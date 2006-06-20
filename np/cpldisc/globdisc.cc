// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* globdisc.cc - global discretization numproc and the corresponding
 * transfer operators numproc.
 * History:
 * Aug. 29, 2002 - created
 * Nov. 13, 2002 - transfer operators
 * May 16, 2003 - number of skip flags can be different from that of equations
 * Jun. 2, 2003 - sparse matrix blocks
 * Jun. 20, 2003 - regular output
 * D. Logashenko
 * Based on ideas of S. Paxion and N. Simus.
 **
 * USAGE:
 * npcreate <name> $c globaldisc;
 * npinit <name>
 *  [$undeclare] # to release the data structure (for ex. from the prev. init.)
 *  {$disc<id> <local discretization numproc>} # the local discretizations
 *  [{$scale<id> <factor>}] # the scaling factor for the loc. disc. <id>
 *  [{$link|copy <export parameter> <import parameter>}] # the correspondence of parameters
 *  [{$set <import parameter> {<value>}}] # explicit specification of values
 *  [{$vd <import parameter> <vec. desc>}] # specification using vec. descriptors
 *  [$sparse] # to assemble sparse matrix blocks (cf. macro __AUTO_SPARSE_DETECT__ below)
 *  [$display no|red|full] # information about the progress of the assembling
 *  [$print_defect] # turnes on the dbg printing of the global defect (finest grid lev. only)
 *  [$print_matrix] # turnes on the dbg printing of the global matrix (finest grid lev. only)
 *  <NP_T_ASSEMBLE class specific parameters>;
 * ...
 * npcreate <transfer name> $c coupled_transfer;
 * npinit <transfer name>
 *  $global <name>
 *  <standard NP_TRANSFER initialization parameters (for ex. - none)>;
 *
 ***REMARK: 1. There are two possible specifications for <export parameter> and
 * <import parameter>: <discretization id>:<export parameter name> or
 * merely <export parameter name>, and the same for the import parameters.
 * <discretization id>'s should be given for both the export and import
 * parameters or omitted for both them. If the discretization id's are not
 * specified, the first found parameter with the given name is used.
 ***2. There are two types of correspondencies between the export and import
 * parameters. The option $link means that the values of the export parameter
 * are copied to the given import parameters and the derivatives are also
 * taken into account when the Jacobian is being computed. The option $copy
 * merely copies the values, the derivatives are then neglected.
 ***3. There are two possibilities to set constant values to import parameters.
 * The first one is to assinge the constants explicitly in the script at the
 * initialization. This is done using the option $set: merely specify the
 * constants. The another way is to assign a vector descriptor. This is allowed
 * only for parameters of the type POSITION_IP, POSITION_BIP, POSITION_NODE and
 * POSITION_ELEMENT. For POSITION_NODE and POSITION_ELEMENT the values are
 * merely read from the vectors. For POSITION_IP and POSITION_BIP the nodal
 * values are read from the v. d. and then interpolated to the integration
 * points.
 */

/* System headers: */
#include <stdio.h>
#include <math.h>
#include <string.h>

/* The UG namespaces: */
#include "namespace.h"
USING_UG_NAMESPACES

/* UG headers: */
#include "compiler.h"
extern "C"
{
  #include "debug.h"
  #include "ugdevices.h"
  #include "general.h"
  #include "gm.h"
  #include "misc.h"
  #include "pcr.h"
  #include "np.h"
  #include "shapes.h"
  #include "fvgeom.h"
  #include "assemble.h"
  #include "transfer.h"
  #include "sm.h" /* for sparse block structures */
}

/* Local discretization header: */
#include "partdisc.hh"

/* Own header: */
extern "C"
{
  #include "globdisc.h"
}

/* Define __AUTO_SPARSE_DETECT__ to switch the automatic detection
 * of sparse blocks using the macrofunction MD_IS_SPARSE. Undefine
 * __AUTO_SPARSE_DETECT__ to use the option $sparse instead:
 */
#undef __AUTO_SPARSE_DETECT__

/** For debugging only: **/
#undef __GLOB_NUM_DIF__
/**/

/*** Connection handling: ***/

struct import_param_list
{
  import_param_list * tl; /* to the next */
  import_param * import; /* the parameter */
  char deriv; /* nonzero iff the derivatives should be computed */
};

struct param_connection
{
  param_connection * tl; /* to the next connection */
  export_param * exp_param; /* the exported param */
  import_param_list * import_list; /* the connected import parameters */
  char deriv; /* nonzero if the derivatives of exp_param should be computed */

  /* Constructor and destructor: */

  param_connection
  (
    param_connection * the_tl, /* the link */
    export_param * the_exp_param /* the export parameter */
  )
  {
    tl = the_tl; exp_param = the_exp_param;
    import_list = NULL;
    deriv = 0;
  };

  ~param_connection ()
  {
    import_param_list * imp_elem;
    while (import_list != NULL)
    {
      imp_elem = import_list;
      import_list = import_list->tl;
      delete imp_elem;
    };
  };

  /* Functions: */

  /* add_import - adds an import parameter to the list. Returns 0 if OK,
   * nonzero if memory allocation failed:
   */
  INT add_import
  (
    import_param * imp_param, /* the import parameter to attach */
    char deriv_imp /* nonzero iff derivatives for this import param. */
  )
  {
    import_param_list * imp_elem;
    if ((imp_elem = new import_param_list) == NULL)
      return __LINE__;
    imp_elem->import = imp_param;
    imp_elem->deriv = deriv_imp;
    deriv |= deriv_imp;
    imp_elem->tl = import_list;
    import_list = imp_elem;
    imp_param->used = 1;
    return 0;
  };

  void copy_values /* copies the values from exp_param to import params */
  (
    FVElementGeometry * fvg
  );

  void display () /* print the information about this connection */
  {
    import_param_list * imp_elem;
    UserWriteF (" - '%d:%s'", exp_param->disc->id, exp_param->name);
    if (! deriv) UserWriteF ("*");
    UserWriteF (" to");
    for (imp_elem = import_list; imp_elem != NULL; imp_elem = imp_elem->tl)
    {
      UserWriteF (" '%d:%s'", imp_elem->import->disc->id,
                  imp_elem->import->name);
      if (! imp_elem->deriv) UserWriteF ("*");
    };
    UserWriteF ("\n");
  };
};

/* param_connection::copy_values - transfers the values through the
 * connection:
 */
void param_connection::copy_values (FVElementGeometry * fvg)
{
  INT i, pos;
  import_param_list * imp_elem;
  INT n_pos = n_position_points (fvg, exp_param->p_type);

  for (imp_elem = import_list; imp_elem != NULL; imp_elem = imp_elem->tl)
    for (pos = 0; pos < n_pos; pos++)
      for (i = 0; i < exp_param->n_comp; i++)
        (* (imp_elem->import))(pos, i) = (* exp_param)(pos, i);
}

/*** Predefined values of import parameters: ***/

/* Generic class for predefined values: */
struct param_predefined
{
  param_predefined * tl; /* the next in the list */
  import_param * param; /* the parameter to assigne the values to */

  /* Constructor: */
  param_predefined (param_predefined * the_tl, import_param * the_param)
  {
    tl = the_tl; param = the_param;
  };

  /* Destructor (virtual): */
  virtual ~param_predefined () {};

  /* set - copies the values kept in the array to the parameter: */
  virtual void set (FVElementGeometry * fvg) = 0;

  /* display - prints the type of the setting and/or the values */
  virtual void display () = 0;
};

/* Generic class for constant values: */
struct param_value : public param_predefined
{
  DOUBLE * value; /* array of the values */

  /* Constructor: */
  param_value
  (
    param_predefined * the_tl, /* to the next element */
    import_param * the_param, /* the import parameter to set */
    DOUBLE * the_val /* the array of the values */
  )
    : param_predefined (the_tl, the_param)
  {
    value = the_val;
  };

  /* set - copies the values kept in the array to the parameter: */
  void set (FVElementGeometry * fvg);

  /* display - prints the values: */
  void display ();
};

/* param_value::set - copies the values kept in the array to the parameter: */
void param_value::set (FVElementGeometry * fvg)
{
  INT i, pos;
  INT n_pos = n_position_points (fvg, param->p_type);

  for (pos = 0; pos < n_pos; pos++)
    for (i = 0; i < param->n_comp; i++)
      (* param)(pos, i) = value [i];
}

/* param_value::display - prints the values from the array 'value': */
void param_value::display ()
{
  INT i;
  UserWriteF (" %s =", param->name);
  for (i = 0; i < param->n_comp; i++)
    UserWriteF (" %g", value [i]);
  UserWriteF (";\n");
};

/* The class for max. DIM values: */
struct param_DIM_value : public param_value
{
  DOUBLE array [DIM];

  /* Constructor: */
  param_DIM_value (param_predefined * the_tl, import_param * the_param)
    : param_value (the_tl, the_param, array) {};
};

/* The class for arbitrary number of values: */
struct param_long_value : public param_value
{
  /* Constructor: */
  param_long_value (param_predefined * the_tl, import_param * the_param)
    : param_value (the_tl, the_param, NULL)
  {
    value = new DOUBLE [param->n_comp];
  };

  /* Destructor: */
  ~param_long_value ()
  {
    if (value != NULL) free (value);
  };
};

/* Class for setting an import parameter from a vector descriptor: */
struct param_vec_desc : public param_predefined
{
  VECDATA_DESC * vd; /* the vector descriptor to take the values from */

  /* Constructor: */
  param_vec_desc
  (
    param_predefined * the_tl, /* to the next element */
    import_param * the_param, /* the parameter to set */
    VECDATA_DESC * the_vd /* the vector descriptor */
  )
    : param_predefined (the_tl, the_param)
  {
    vd = the_vd;
  };

  /* set - copies the values from the v. d. to the parameter: */
  void set (FVElementGeometry * fvg);

  /* display - prints the name of the v. d.: */
  void display ();
};

/* param_vec_desc::set - copies the values from the v. d. to the parameter: */
void param_vec_desc::set (FVElementGeometry * fvg)
{
  INT n;
  INT i, pos;
  SHORT * comp;
  VECTOR * vec;
  DOUBLE nd_val [MAXNC];

  switch (param->p_type)
  {
  case POSITION_IP :
    /* Interpolate from nodal values: */
    comp = VD_ncmp_cmpptr_of_otype (vd, NODEVEC, &n);
    if (n != param->n_comp) return;
    for (i = 0; i < n; i++)
    {
      for (pos = 0; pos < FVG_NSCV (fvg); pos++)
        nd_val[pos] = VVALUE (NVECTOR (CORNER (FVG_ELEM (fvg), pos)), comp [i]);
      for (pos = 0; pos < FVG_NSCVF (fvg); pos++)
        if (InterpolateFEFunction (DIM, FVG_TAG (fvg),
                                   SCVF_LIP (FVG_SCVF (fvg, i)), nd_val, & ((* param)(pos, i))))
          return;
    };
    break;

  case POSITION_BIP :
    /* Interpolate from nodal values: */
    comp = VD_ncmp_cmpptr_of_otype (vd, NODEVEC, &n);
    if (n != param->n_comp) return;
    for (i = 0; i < n; i++)
    {
      for (pos = 0; pos < FVG_NSCV (fvg); pos++)
        nd_val[pos] = VVALUE (NVECTOR (CORNER (FVG_ELEM (fvg), pos)), comp [i]);
      for (pos = 0; pos < FVG_NSCVBF (fvg); pos++)
        if (InterpolateFEFunction (DIM, FVG_TAG (fvg),
                                   SCVBF_LIP (FVG_SCVBF (fvg, i)), nd_val, & ((* param)(pos, i))))
          return;
    };
    break;

  case POSITION_NODE :
    comp = VD_ncmp_cmpptr_of_otype (vd, NODEVEC, &n);
    if (n != param->n_comp) return;
    for (pos = 0; pos < FVG_NSCV (fvg); pos++)
    {
      vec = NVECTOR (CORNER (FVG_ELEM (fvg), pos));
      for (i = 0; i < n; i++)
        (* param)(pos, i) = VVALUE (vec, comp [i]);
    };
    break;

  case POSITION_ELEMENT :
    comp = VD_ncmp_cmpptr_of_otype (vd, ELEMVEC, &n);
    if (n != param->n_comp) return;
    vec = EVECTOR (FVG_ELEM (fvg));
    for (i = 0; i < n; i++)
      (* param)(pos, i) = VVALUE (vec, comp [i]);
    break;

    /* no other types implemented here */
  };
}

/* param_vec_desc::display - prints the name of the v. d.: */
void param_vec_desc::display ()
{
  if (param->p_type == POSITION_IP || param->p_type == POSITION_BIP)
    UserWriteF (" %s is interpolated from v. d. '%s'\n",
                param->name, ENVITEM_NAME (vd));
  else
    UserWriteF (" %s is read from v. d. '%s'\n",
                param->name, ENVITEM_NAME (vd));
}

/*** Local discretization: ***/

struct t_disc_list
{
  t_disc_list * tl; /* the next element */
  np_part_discretization * disc; /* the discretization */
  DOUBLE scale; /* the scaling factor for the discretization */
  char scale_flag; /* zero iff scale == 1 (may be always set to nonzero) */
};

/***** The numproc class: *****/

struct np_coupled_global_disc : public np_t_assemble
{
  t_disc_list * disc_list; /* the local discretizations to couple */
  param_connection * connection_list; /* the declared connections */
  param_predefined * value_list; /* the predefined parameter values */
  INT n_equations; /* number of all equations in all the discretizations */

# ifndef __AUTO_SPARSE_DETECT__
  INT assemble_sparse; /* nonzero to assemble sparse matrix blocks */
# endif

  /* Output control flags: */
  INT display_progress; /* for printing of the processed grid levels */
  INT print_defect; /* debugging output of the assembled defect */
  INT print_matrix; /* debugging output of the assembled preconditioner */

  /* Functions: */

  INT add_discretization
  (
    np_part_discretization * disc, /* the discretization to add */
    INT id, /* the id for this discretization */
    DOUBLE scale /* the scaling factor */
  );
  void sort_discretizations ();
  void assign_VD_ind ();
  void clear_used_flags ();
  INT check_skip_flag_ind ();
  INT enter_discretizations
  (
    INT argc, /* number of strings */
    char * * argv /* the strings */
  );
  void print_discretizations ();
  void undeclare_discretizations ();

  INT create_connection
  (
    import_param * imp_param, /* the import parameter */
    export_param * exp_param, /* the export parameter */
    char deriv /* nonzero iff derivatives for this link */
  );
  INT read_connections
  (
    INT argc, /* number of arguments */
    char * * argv /* the arguments */
  );
  INT sort_connections ();
  void print_connections ();
  void undeclare_connections ();

  INT parse_param_value
  (
    char * command /* the text specification of the values */
  );
  INT parse_param_vec_desc
  (
    char * command /* the text specification of the values */
  );
  INT read_param_predefined
  (
    INT argc, /* number of arguments */
    char * * argv /* the arguments */
  );
  void print_param_predefined ();
  void undeclare_param_predefined ();
  INT verify_sparsity
  (
    MATDATA_DESC * md /* the vector descriptor to check */
  );
};

#define NP_COUPLED_GLOBAL_DISC_CLASS_NAME T_ASSEMBLE_CLASS_NAME ".globaldisc"

/*** Menagment of the local discretizations: ***/

/* add_discretization - adds a discretization to the list. Returns
 * 0 if OK, nonzero on an error:
 */
INT np_coupled_global_disc::add_discretization
(
  np_part_discretization * disc, /* the discretization to add */
  INT id, /* the id for this discretization */
  DOUBLE scale /* the scaling factor */
)
{
  t_disc_list * elem;

  for (elem = disc_list; elem != NULL; elem = elem->tl)
    if (elem->disc->id == id)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "Redeclaration of discretization id %d", id);
      return __LINE__;
    };

  if ((elem = new t_disc_list) == NULL)
  {
    PrintErrorMessage ('E', "globdisc",
                       "Out of memory when declaring new discretization");
    return __LINE__;
  };

  elem->tl = disc_list;
  elem->disc = disc;
  elem->scale = scale;
  elem->scale_flag = (scale != 1);
  disc->id = id;
  disc_list = elem;

  return 0;
}

/* sort_discretizations - sorts the discretizations in the list
 * according to their id's:
 */
void np_coupled_global_disc::sort_discretizations ()
{
  t_disc_list * list, * cur_elem, * place;

  if (disc_list == NULL) return; /* nothing to sort */

  list = disc_list; disc_list = disc_list->tl; list->tl = NULL;
  while (disc_list != NULL)
  {
    cur_elem = disc_list; disc_list = disc_list->tl;
    if (cur_elem->disc->id < list->disc->id)
    {
      cur_elem->tl = list; list = cur_elem;
    }
    else
    {
      for (place = list; place->tl != NULL; place = place->tl)
        if (cur_elem->disc->id < place->tl->disc->id)
          break;
      cur_elem->tl = place->tl; place->tl = cur_elem;
    };
  };

  disc_list = list;
}

/* assign_VD_ind - assignes the indices of components in the
 * solution/rhs vector descriptors to the discretizations.
 * Counts the number of the components:
 */
void np_coupled_global_disc::assign_VD_ind ()
{
  t_disc_list * elem;
  INT cur_cmp;

  cur_cmp = 0;
  for (elem = disc_list; elem != NULL; elem = elem->tl)
  {
    elem->disc->first_vec_cmp = cur_cmp;
    cur_cmp += elem->disc->n_equations;
  };
  n_equations = cur_cmp;
}

/* clear_used_flags - sets the 'used' flags for all import parameters
 * of all the discretizations to zero:
 */
void np_coupled_global_disc::clear_used_flags ()
{
  t_disc_list * elem;
  import_param * param;

  for (elem = disc_list; elem != NULL; elem = elem->tl)
    for (param = elem->disc->import_list; param != NULL; param = param->tl)
      param->used = 0;
}

/* check_skip_flag_ind - verifies that all the discretizations except
 * of the last one have exactly as many skip flags as equations. The
 * function returns zero if this condition is satisfied and the number
 * of the requested skip flags is less or equal sizeof (INT) * 8, nonzero
 * otherwise:
 */
INT np_coupled_global_disc::check_skip_flag_ind ()
{
  t_disc_list * elem;
  INT skip_flag_ind;

  if (disc_list == NULL) return 0;

  /* Assign the indices of the skip flags: */
  skip_flag_ind = 0;
  for (elem = disc_list; elem->tl != NULL; elem = elem->tl)
    if (elem->disc->n_skip_flags != elem->disc->n_equations)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "Disc. %d (%s) is not the last but has different numbers of"
                          " equations and skip flags", elem->disc->id, ENVITEM_NAME (elem->disc));
      return __LINE__;
    }
    else
      skip_flag_ind += elem->disc->n_skip_flags;

  if (skip_flag_ind > sizeof (INT) * 8)
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Number of skip flags required (%d) is too large", skip_flag_ind);
    return __LINE__;
  };

  return 0;
}

/* enter_discretizations - reads the discretizations, sorts them
 * and initialize their indices. Returns 0 if OK, nonzero on an error:
 */
INT np_coupled_global_disc::enter_discretizations
(
  INT argc, /* number of strings */
  char * * argv /* the strings */
)
{
  INT i;
  INT id;
  np_part_discretization * disc;
  char disc_name [64];
  DOUBLE scale;

  for (i = 0; i < argc; i++)
    if (sscanf (argv [i], "disc%d %s", &id, disc_name) == 2)
    {
      if ((disc = (np_part_discretization *) GetNumProcByName (NP_MG (this),
                                                               disc_name, NP_PART_DISCRETIZATION_CLASS_NAME)) == NULL)
      {
        PrintErrorMessageF ('E', "globdisc", "Cannot get discretization '%s'",
                            disc_name);
        return __LINE__;
      };
      sprintf (disc_name, "scale%d", id);
      if (ReadArgvDOUBLE (disc_name, &scale, argc, argv))
        scale = 1;
      if (add_discretization (disc, id, scale))
        return __LINE__;
    };

  sort_discretizations ();
  assign_VD_ind ();
  clear_used_flags ();
  return check_skip_flag_ind ();
}

/* print_discretizations - prints the information about the discretizations: */
void np_coupled_global_disc::print_discretizations ()
{
  t_disc_list * elem;

  if (disc_list == NULL)
  {
    UserWriteF ("No discretizations declared.\n");
    return;
  };

  UserWriteF ("Registered discretizations:\n");
  for (elem = disc_list; elem != NULL; elem = elem->tl)
    UserWriteF (" %d: '%s', %d eq., scaled %lg, bc_id = %d\n",
                elem->disc->id, ENVITEM_NAME (elem->disc), elem->disc->n_equations,
                (double) elem->scale, elem->disc->bc_id);
}

/* undeclare_discretizations - deletes the declaration structures
 * of the discretizations:
 */
void np_coupled_global_disc::undeclare_discretizations ()
{
  t_disc_list * disc_elem;

  while (disc_list != NULL)
  {
    disc_elem = disc_list;
    disc_list = disc_list->tl;
    delete disc_elem;
  };
}

/*** Menagement of the connections: ***/

/* create_connection - creates a connection between an import and an
 * export parameter. The function returns zero if OK, nonzero on an
 * error:
 */
INT np_coupled_global_disc::create_connection
(
  import_param * imp_param, /* the import parameter */
  export_param * exp_param, /* the export parameter */
  char deriv /* nonzero iff derivatives for this link */
)
{
  param_connection * conn;

  /* Check compartibility: */
  if (imp_param->p_type != exp_param->p_type
      || imp_param->n_comp != exp_param->n_comp)
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Export param. '%s' (%s) is incompartible with import param. '%s' (%s)",
                        exp_param->name, ENVITEM_NAME (exp_param->disc),
                        imp_param->name, ENVITEM_NAME (imp_param->disc));
    return __LINE__;
  };

  /* Look for the connection structure; create it if not found: */
  for (conn = connection_list; conn != NULL; conn = conn->tl)
    if (conn->exp_param->disc->id == exp_param->disc->id
        && strcmp (conn->exp_param->name, exp_param->name) == 0)
      break;
  if (conn == NULL)
  {
    if ((conn = new param_connection (connection_list, exp_param)) == NULL)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "Cannot create connection for export param. '%s'", exp_param->name);
      return __LINE__;
    };
    connection_list = conn;
  };
  /* Add the import parameter to the connection: */
  if (conn->add_import (imp_param, deriv))
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Cannot link export param. '%s' to import param '%s'",
                        exp_param->name, imp_param->name);
    return __LINE__;
  };

  return 0;
}

/* find_import_param - looks for an import parameter in a dicsretization.
 * Returns the pointer to this parameter or NULL if not found:
 */
static import_param * find_import_param
(
  np_part_discretization * disc, /* the discretization */
  char * name /* the name of the parameter */
)
{
  import_param * param;

  for (param = disc->import_list; param != NULL; param = param->tl)
    if (strcmp (param->name, name) == 0)
      return param;

  return NULL;
}

/* find_export_param - looks for an export parameter in a dicsretization.
 * Returns the pointer to this parameter or NULL if not found:
 */
static export_param * find_export_param
(
  np_part_discretization * disc, /* the discretization */
  char * name /* the name of the parameter */
)
{
  export_param * param;

  for (param = disc->export_list; param != NULL; param = param->tl)
    if (strcmp (param->name, name) == 0)
      return param;

  return NULL;
}

/* read_connections - reads declared connections from the
 * script command line. The function returns 0 if OK, nonzero on an error:
 */
INT np_coupled_global_disc::read_connections
(
  INT argc, /* number of arguments */
  char * * argv /* the arguments */
)
{
  INT i;
  char imp_name [64], exp_name [64];
  INT exp_id, imp_id;
  import_param * imp_param;
  export_param * exp_param;
  t_disc_list * disc_elem;
  char deriv;

  if (disc_list == NULL)
  {
    PrintErrorMessage ('E', "globdisc", "No discretization declared!");
    return __LINE__;
  };

  for (i = 0; i < argc; i++)
    if ((deriv = 1, sscanf (argv[i], "link %d:%s %d:%s", &exp_id, exp_name,
                            &imp_id, imp_name) == 4)
        || (deriv = 0, sscanf (argv[i], "copy %d:%s %d:%s", &exp_id, exp_name,
                               &imp_id, imp_name) == 4))
    {
      /* Find the discretisation for the import parameter: */
      disc_elem = disc_list;
      while (disc_elem->disc->id != imp_id)
        if ((disc_elem = disc_elem->tl) == NULL)
        {
          PrintErrorMessageF ('E', "globdisc",
                              "No disc. id %d for the decl. '%s'", imp_id, argv[i]);
          return __LINE__;
        };
      /* Find the import parameter: */
      if ((imp_param = find_import_param (disc_elem->disc, imp_name)) == NULL)
      {
        PrintErrorMessageF ('E', "globdisc",
                            "Cannot find import param. '%s' for disc. id %d", imp_name, imp_id);
        return __LINE__;
      };
      /* Find the discretisation for the export parameter: */
      disc_elem = disc_list;
      while (disc_elem->disc->id != exp_id)
        if ((disc_elem = disc_elem->tl) == NULL)
        {
          PrintErrorMessageF ('E', "globdisc",
                              "No disc. id %d for the decl. '%s'", exp_id, argv[i]);
          return __LINE__;
        };
      /* Find the export parameter: */
      if ((exp_param = find_export_param (disc_elem->disc, exp_name)) == NULL)
      {
        PrintErrorMessageF ('E', "globdisc",
                            "Cannot find export param. '%s' for disc. id %d", imp_name, imp_id);
        return __LINE__;
      };
      /* Create the connection: */
      if (create_connection (imp_param, exp_param, deriv))
        return __LINE__;
    }
    else if ((deriv = 1, sscanf (argv[i], "link %s %s", exp_name, imp_name)
              == 2)
             || (deriv = 0, sscanf (argv[i], "copy %s %s", exp_name, imp_name) == 2))
    {
      /* Look for the import parameter: */
      disc_elem = disc_list;
      while ((imp_param = find_import_param (disc_elem->disc, imp_name)) == NULL)
        if ((disc_elem = disc_elem->tl) == NULL)
        {
          PrintErrorMessageF ('E', "globdisc",
                              "Cannot find import parameter '%s'", imp_name);
          return __LINE__;
        };
      /* Look for the export parameter: */
      disc_elem = disc_list;
      while ((exp_param = find_export_param (disc_elem->disc, exp_name)) == NULL)
        if ((disc_elem = disc_elem->tl) == NULL)
        {
          PrintErrorMessageF ('E', "globdisc",
                              "Cannot find export parameter '%s'", exp_name);
          return __LINE__;
        };
      /* Create the connection: */
      if (create_connection (imp_param, exp_param, deriv))
        return __LINE__;
    };

  return 0;
}

/* sort_connections - orders the connections so that the prerequisites
 * of a connection are located earlier in the list than the connection itself.
 * The function returns nonzero if the connections cannot be ordered in
 * this way (i.e. there are cyclic dependencies):
 */
INT np_coupled_global_disc::sort_connections ()
{
  INT i, ret_val;
  param_connection * list, * rest, * conn, * t;
  import_param * prereq;
  import_param_list * i_param;
  char updated, move_conn;

  list = NULL;
  ret_val = 0;
  while (connection_list != NULL)
  {
    rest = connection_list;
    connection_list = NULL;
    updated = 0;
    while (rest != NULL)
    {
      conn = rest; rest = rest->tl;
      move_conn = 1;
      if (conn->exp_param->prerequisite != NULL)
        for (i = 0; (prereq = conn->exp_param->prerequisite [i]) != NULL; i++)
          if (prereq->used)
          {
            for (t = list; t != NULL; t = t->tl)
              for (i_param = t->import_list; i_param != NULL; i_param = i_param->tl)
                if (i_param->import == prereq)
                  goto prereq_found;
            /* We have not found the prerequisite in the new list.
             * This connection should wait:
             */
            move_conn = 0;
            break;
prereq_found:;
          };
      if (move_conn)
      {
        conn->tl = list; list = conn;
        updated = 1;
      }
      else
      {
        conn->tl = connection_list; connection_list = conn;
      };
    };
    if (! updated)
    {  /* We have not moved any connections! Then there are cycles. */
      PrintErrorMessageF ('E', "globdisc",
                          "Cyclic dependence of the parameters! Check export param '%s' (%s)",
                          connection_list->exp_param->name,
                          ENVITEM_NAME (connection_list->exp_param->disc));
      ret_val = __LINE__;
      break;
    };
  };

  /* Write the connections from 'list' to 'connection_list' in the
   * inverse order:
   */
  while (list != NULL)
  {
    conn = list; list = list->tl;
    conn->tl = connection_list; connection_list = conn;
  };

  return ret_val;
}

/* print_connections - prints the declared connections: */
void np_coupled_global_disc::print_connections ()
{
  param_connection * conn;

  if (connection_list == NULL)
  {
    UserWriteF ("No connections declared.\n");
    return;
  };

  UserWriteF ("Connections:\n");
  for (conn = connection_list; conn != NULL; conn = conn->tl)
    conn->display ();
}

/* undeclare_connections - deletes the declaration structures of the
 * connections:
 */
void np_coupled_global_disc::undeclare_connections ()
{
  param_connection * conn;

  while (connection_list != NULL)
  {
    conn = connection_list;
    connection_list = connection_list->tl;
    delete conn;
  };
}

/*** Menagement of the predefined parameter values: ***/

/* The name of the specification of the values in scripts: */
static char set_command [] = "set";
static char vd_command [] = "vd";

/* parse_param_value - parses the specification of constant values for a given
 * import parameter and creates a structure 'param_value'. Returns zero
 * if OK, nonzero on errors. Prints the corresponding warnings:
 */
INT np_coupled_global_disc::parse_param_value
(
  char * command /* the text specification of the values */
)
{
  t_disc_list * disc_elem;
  INT i, disc_id;
  import_param * param;
  char param_name [64];
  param_value * p_val;
  double t;
  int n;
  char * str = command + sizeof (set_command);

  /*** Get the import parameter: ***/
  if (sscanf (str, " %d:%s %n", &disc_id, param_name, &n) == 2)
  {
    /* Find the discretisation for the import parameter: */
    disc_elem = disc_list;
    while (disc_elem->disc->id != disc_id)
      if ((disc_elem = disc_elem->tl) == NULL)
      {
        PrintErrorMessageF ('E', "globdisc",
                            "No disc. id %d for the setting '%s'", disc_id, command);
        return __LINE__;
      };
    /* Find the import parameter: */
    if ((param = find_import_param (disc_elem->disc, param_name)) == NULL)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "Cannot find import param. '%s' for disc. id %d in '%s'",
                          param_name, disc_id, command);
      return __LINE__;
    };
  }
  else if (sscanf (str, " %s %n", param_name, &n) == 1)
  {
    /* Look for the import parameter: */
    disc_elem = disc_list;
    while ((param = find_import_param (disc_elem->disc, param_name)) == NULL)
      if ((disc_elem = disc_elem->tl) == NULL)
      {
        PrintErrorMessageF ('E', "globdisc",
                            "Cannot find import parameter '%s' for '%s'", param_name, command);
        return __LINE__;
      };
  }
  else
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Cannot read the import parameter name in '%s'", command);
    return __LINE__;
  };

  /*** Create the 'param_value' structure: ***/
  if (param->n_comp <= DIM)
    p_val = new param_DIM_value (value_list, param);
  else
    p_val = new param_long_value (value_list, param);

  if (p_val == NULL)
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Cannot allocate memory when processing '%s'", command);
    return __LINE__;
  };

  /*** Read the values: ***/
  str += n;
  for (i = 0; i < param->n_comp; i++)
  {
    if (sscanf (str, " %lg%n", &t, &n) != 1)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "Cannot read the value #%d (of %d) in '%s'", i + 1, param->n_comp,
                          command);
      return __LINE__;
    };
    p_val->value [i] = t;
    str += n;
  };

  /*** Insert the structure into the list: ***/
  value_list = p_val;
  param->used = 2;

  return 0;
}

/* parse_param_vec_desc - parses the command that assignes a vector descriptor
 * to a parameter. The type of the vector descriptor is checked. The
 * function returns 0 if OK, nonzero on an error:
 */
INT np_coupled_global_disc::parse_param_vec_desc
(
  char * command /* the text specification of the v. d. */
)
{
  t_disc_list * disc_elem;
  import_param * param;
  VECDATA_DESC * vd;
  param_vec_desc * param_vd;
  INT disc_id;
  char param_name [64], vd_name [64];
  INT n;
  char * str = command + sizeof (vd_command);

  /* Read the names of the parameter and the vector descriptor: */
  if (sscanf (str, " %d:%s %s", &disc_id, param_name, vd_name) == 3)
  {
    /* Find the discretisation for the import parameter: */
    disc_elem = disc_list;
    while (disc_elem->disc->id != disc_id)
      if ((disc_elem = disc_elem->tl) == NULL)
      {
        PrintErrorMessageF ('E', "globdisc",
                            "No disc. id %d for the command '%s'", disc_id, command);
        return __LINE__;
      };
    /* Find the import parameter: */
    if ((param = find_import_param (disc_elem->disc, param_name)) == NULL)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "Cannot find import param. '%s' for disc. id %d in '%s'",
                          param_name, disc_id, command);
      return __LINE__;
    };
  }
  else if (sscanf (str, " %s %s", param_name, vd_name) == 2)
  {
    /* Look for the import parameter: */
    disc_elem = disc_list;
    while ((param = find_import_param (disc_elem->disc, param_name)) == NULL)
      if ((disc_elem = disc_elem->tl) == NULL)
      {
        PrintErrorMessageF ('E', "globdisc",
                            "Cannot find import parameter '%s' for '%s'", param_name, command);
        return __LINE__;
      };
  }
  else
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Cannot read the parameter or the vector descriptor in '%s'", command);
    return __LINE__;
  };

  /* Get the vector descriptor and check its type: */
  if ((vd = GetVecDataDescByName (NP_MG (this), vd_name)) == NULL)
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Cannot get the v. d. '%s' in '%s'", vd_name, command);
    return __LINE__;
  };

  switch (param->p_type)
  {
  case POSITION_NODE :
  case POSITION_IP :
  case POSITION_BIP :
    VD_ncmp_cmpptr_of_otype (vd, NODEVEC, &n);
    if (n != param->n_comp)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "V. d. '%s' should have %d cmps. at nodes for param. '%s'",
                          vd_name, param->n_comp, param->name);
      return __LINE__;
    };
    break;

  case POSITION_ELEMENT :
    VD_ncmp_cmpptr_of_otype (vd, ELEMVEC, &n);
    if (n != param->n_comp)
    {
      PrintErrorMessageF ('E', "globdisc",
                          "V. d. '%s' should have %d cmps. at elements for param. '%s'",
                          vd_name, param->n_comp, param->name);
      return __LINE__;
    };
    break;

  default :
    PrintErrorMessageF ('E', "globdisc",
                        "Param. '%s' cannot be initialized with a v. d. because of its type",
                        param->name);
    return __LINE__;
  };

  /* Create the param_vec_desc structure: */
  if ((param_vd = new param_vec_desc (value_list, param, vd)) == NULL)
  {
    PrintErrorMessageF ('E', "globdisc",
                        "Cannot allocate memory when processing '%s'", command);
    return __LINE__;
  };
  value_list = param_vd;
  param->used = 3;

  return 0;
}

/* read_param_predefined - reads the specifications of the parameter values
 * from the script. Returns zero if OK, nonzero on an error:
 */
INT np_coupled_global_disc::read_param_predefined
(
  INT argc, char * * argv /* the script option strings */
)
{
  INT i, r;

  for (i = 0; i < argc; i++)
    if (strncmp (argv[i], set_command, sizeof (set_command) - 1) == 0)
    {
      if ((r = parse_param_value (argv[i])) != 0)
        return r;
    }
    else if (strncmp (argv[i], vd_command, sizeof (vd_command) - 1) == 0)
      if ((r = parse_param_vec_desc (argv[i])) != 0)
        return r;
  return 0;
}

/* print_param_predefined - prints the predefined values of parameters: */
void np_coupled_global_disc::print_param_predefined ()
{
  param_predefined * t;

  if (value_list == NULL) return;

  UserWriteF ("Parameter values specified by the user:\n");
  for (t = value_list; t != NULL; t = t->tl)
    t->display ();
}

/* undeclare_param_predefined - deletes the 'param_predefined' structures from
 * the list:
 */
void np_coupled_global_disc::undeclare_param_predefined ()
{
  param_predefined * t;

  while (value_list != NULL)
  {
    t = value_list;
    value_list = t->tl;
    delete t;
  };
}

/* verify_sparsity - checks whether the structure of the sparse
 * matrix suits the patterns and the order of the part discretizations
 * registered in this numproc. Returns 0 if OK, nonzero on an error:
 */
INT np_coupled_global_disc::verify_sparsity
(
  MATDATA_DESC * md /* the vector descriptor to check */
)
{
  INT i, j;
  SHORT * comp;
  SHORT * d_comp;
  np_part_discretization * disc;
  t_disc_list * disc_elem;
  INT d_p_d, p_d;

  /* Check whether the used patterns are available: */
  if (n_equations != MD_ROWS_IN_MTYPE (md, MTP (0, 0))
      || n_equations != MD_COLS_IN_MTYPE (md, MTP (0, 0)))
  {
    PrintErrorMessage ('E', "AssemblePrecond",
                       "Wrong size for the nodal sparse blocks");
    return __LINE__;
  };

  comp = MD_MCMPPTR_OF_MTYPE (md, MTP (0, 0));
  d_comp = MD_MCMPPTR_OF_MTYPE (md, DMTP (0));

  /* Check that there are no other patterns: */
  for (i = 0; i < NMATTYPES; i++)
    if (i != DMTP (0) && i != MTP (0, 0) && MD_SM (md, i) != NULL)
    {
      PrintErrorMessage ('E', "AssemblePrecond",
                         "There are sparse matrix types different from the nodal ones");
      return __LINE__;
    };

  /* Check the Jacobian patterns: */
  for (disc_elem = disc_list; disc_elem != NULL; disc_elem = disc_elem->tl)
  {
    disc = disc_elem->disc;
    for (i = 0; i < disc->n_equations; i++)
      for (j = 0; j < disc->n_equations; j++)
      {
        if (disc->jac_sparse_pattern == NULL)
          p_d = d_p_d = 1;
        else
        {
          p_d = disc->jac_sparse_pattern (disc, 0, i, j);
          d_p_d = disc->jac_sparse_pattern (disc, 1, i, j);
        };
        if (p_d && comp [i * n_equations + j] < 0)
        {
          PrintErrorMessageF ('E', "AssemblePrecond",
                              "The offd. spars. pttrn of m. d. '%s' is incompart. in disc. '%s' (id %d)",
                              ENVITEM_NAME (md), ENVITEM_NAME (disc), disc->id);
          return __LINE__;
        };
        if (d_p_d && d_comp [i * n_equations + j] < 0)
        {
          PrintErrorMessageF ('E', "AssemblePrecond",
                              "The diag. spars. pttrn of m. d. '%s' is incompart. in disc. '%s' (id %d)",
                              ENVITEM_NAME (md), ENVITEM_NAME (disc), disc->id);
          return __LINE__;
        };
      };
  };

  return 0;
}

/*** Standard numproc functions: ***/

static INT Init (NP_BASE * base, INT argc, char * * argv)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) base;
  INT res;

  /* Delete the old initializations: */
  np->undeclare_discretizations ();
  np->undeclare_connections ();
  np->undeclare_param_predefined ();

  if (ReadArgvOption ("undeclare", argc, argv))
  {
    /* If said 'undeclare' than that's all: */
    return NP_ACTIVE;
  };

  /* Initialize the base class: */
  res = NPTAssembleInit (base, argc, argv);
  if (res != NP_ACTIVE && res != NP_EXECUTABLE)
    return res;

  /* Read the discretization specifications: */
  if (np->enter_discretizations (argc, argv))
    return NP_NOT_ACTIVE;

  /* Read the connections: */
  if (np->read_connections (argc, argv))
    return NP_NOT_ACTIVE;

  /* Order connections: */
  if (np->sort_connections ())
    return NP_NOT_ACTIVE;

  /* Read the predefined values: */
  if (np->read_param_predefined (argc, argv))
    return NP_NOT_ACTIVE;

# ifndef __AUTO_SPARSE_DETECT__
  np->assemble_sparse = ReadArgvOption ("sparse", argc, argv);
# endif

  /* Read the output control options: */
  np->display_progress = ReadArgvDisplay (argc, argv);
  np->print_defect = ReadArgvOption ("print_defect", argc, argv);
  np->print_matrix = ReadArgvOption ("print_matrix", argc, argv);

  return res;
}

static INT Display (NP_BASE * base)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) base;
  INT res;

  UserWriteF ("Global discretization assembling numproc.\n");

  if ((res = NPTAssembleDisplay (base)) != 0)
    return res;

  np->print_discretizations ();
  np->print_connections ();
  np->print_param_predefined ();

  return 0;
}

/*** NP_T_ASSEMBLE's functions: ***/

/* PreProcess - calls the preprocess functions of the discretizations.
 * Returns zero if all of them returned zero, nonzero otherwise:
 */
static INT PreProcess (NP_T_ASSEMBLE * ass, INT fl, INT tl, DOUBLE t_new,
                       DOUBLE t, DOUBLE t_old, VECDATA_DESC * u_new, VECDATA_DESC * u,
                       VECDATA_DESC * u_old, INT * res)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) ass;
  t_disc_list * disc_elem;

  if (np->display_progress == PCR_FULL_DISPLAY)
    UserWriteF ("PreProcess assembling: lev. %d-%d, t_new=%g, t=%g, t_old=%g\n",
                fl, tl, t_new, t, t_old);

  /* Call the preprocess's for all the discretizations: */
  for (disc_elem = np->disc_list; disc_elem != NULL;
       disc_elem = disc_elem->tl)
    if (disc_elem->disc->preprocess != NULL)
      if (disc_elem->disc->preprocess (disc_elem->disc, fl, tl,
                                       t_new, t, t_old, u_new, u, u_old))
      {
        * res = __LINE__; return __LINE__;
      };

  * res = 0; return 0;
}

/* PostProcess - calls the postprocess functions of the discretizations.
 * Returns zero if all of them returned zero, nonzero otherwise:
 */
static INT PostProcess (NP_T_ASSEMBLE * ass, INT fl, INT tl, DOUBLE t_new,
                        DOUBLE t, DOUBLE t_old, VECDATA_DESC * u_new, VECDATA_DESC * u,
                        VECDATA_DESC * u_old, INT * res)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) ass;
  t_disc_list * disc_elem;

  if (np->display_progress == PCR_FULL_DISPLAY)
    UserWriteF ("PostProcess assembling: lev. %d-%d, t_new=%g, t=%g, t_old=%g\n",
                fl, tl, t_new, t, t_old);

  /* Call the preprocess's for all the discretizations: */
  for (disc_elem = np->disc_list; disc_elem != NULL;
       disc_elem = disc_elem->tl)
    if (disc_elem->disc->postprocess != NULL)
      if (disc_elem->disc->postprocess (disc_elem->disc, fl, tl,
                                        t_new, t, t_old, u_new, u, u_old))
      {
        * res = __LINE__; return __LINE__;
      };

  * res = 0; return 0;
}

/* AssembleInitial - should assemble the initial condition. Calls the
 * corresponding functions of the part discretizations (if supplied).
 * Returns 0 if OK, nonzero on an error:
 */
static INT AssembleInitial (NP_T_ASSEMBLE * ass, INT fl, INT tl, DOUBLE time,
                            VECDATA_DESC * u, INT * res)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) ass;
  np_part_discretization * disc;
  t_disc_list * disc_elem;

  for (disc_elem = np->disc_list; disc_elem != NULL; disc_elem = disc_elem->tl)
    if ((disc = disc_elem->disc)->assemble_initial != NULL)
      if (disc->assemble_initial (disc, fl, tl, time, u))
      {
        PrintErrorMessageF ('E', "AssembleInitial",
                            "'assemble_initial' of '%s' (id %d) returned error",
                            ENVITEM_NAME (disc), disc->id);
        * res = __LINE__; return __LINE__;
      };

  * res = 0; return 0;
}

/* AssembleSolution - assembles the boundary values, sets the skip bits,
 * ... Returns 0 if OK, nonzero on an error:
 */
static INT AssembleSolution (NP_T_ASSEMBLE * ass, INT fl, INT tl, DOUBLE time,
                             VECDATA_DESC * u, INT * res)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) ass;
  SHORT * comp;
  INT n_comp;
  np_part_discretization * disc;
  INT level;
  GRID * grid;
  NODE * node;
  VECTOR * vec;
  t_disc_list * disc_elem;

  /* Get the components of the vector descriptor: */
  comp = VD_ncmp_cmpptr_of_otype (u, NODEVEC, &n_comp);
  if (n_comp != np->n_equations)
  {
    PrintErrorMessageF ('E', "AssembleSolution",
                        "Wrong vector descriptor or format: only %d comp. specified, but %d needed",
                        n_comp, np->n_equations);
    * res = __LINE__; return __LINE__;
  };

  if (np->display_progress == PCR_FULL_DISPLAY)
    UserWriteF ("AssembleSolution: ");

  /* Loop the nodes: */
  for (level = fl; level <= tl; level++)
  {
    if ((grid = GRID_ON_LEVEL (NP_MG (np), level)) == NULL)
    {
      * res = __LINE__; return __LINE__;
    };
    for (node = FIRSTNODE (grid); node != NULL; node = SUCCN (node))
    {
      vec = NVECTOR (node);
      vec->skip = 0;
      /* Loop the discretizations: */
      for (disc_elem = np->disc_list; disc_elem != NULL; disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        /* Assemble the local solution in the discretization: */
        if (disc->assemble_solution (disc, node, time))
        {
          PrintErrorMessageF ('E', "AssembleSolution",
                              "assemble_solution of '%s' (id %d) returned an error",
                              ENVITEM_NAME (disc), disc->id);
          * res = __LINE__; return __LINE__;
        };
        /* Set the values: */
        for (n_comp = 0; n_comp < disc->n_equations; n_comp++)
          disc->boundary_value (disc, n_comp,
                                & VVALUE (vec, comp [disc->first_vec_cmp + n_comp]));
        /* Set the skip flags: */
        for (n_comp = 0; n_comp < disc->n_skip_flags; n_comp++)
          if (disc->skip_flag (disc, n_comp))
            vec->skip |= 1 << (disc->first_vec_cmp + n_comp);
        if (disc->n_skip_flags < disc->n_equations)
        /* this may happen for the last disc. only! */
        {
          n_comp = disc->first_vec_cmp + disc->n_skip_flags - 1;
          if (vec->skip & (1 << n_comp)) /* if the last flag set */
            vec->skip |= (~ (INT) 1) << n_comp; /* set the other ones */
        };
      };
    };
    /* The (reduced) output: */
    if (np->display_progress >= PCR_RED_DISPLAY)
      UserWriteF ("[S:%d]", level);
  };

  if (np->display_progress >= PCR_RED_DISPLAY)
    UserWriteF ("\n");

  * res = 0; return 0;
}

/* AssembleDefect - assembles the global defect all the discretizations.
 * The defect is added to that kept in the vector descriptor.
 * Returns 0 if OK, nonzero on an error:
 */
static INT AssembleDefect (NP_T_ASSEMBLE * ass, INT fl, INT tl,
                           DOUBLE time, DOUBLE s_m, DOUBLE s_a, VECDATA_DESC * u, VECDATA_DESC * d,
                           MATDATA_DESC * J, INT * res)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) ass;
  INT level;
  GRID * grid;
  ELEMENT * elem;
  INT i, co;
  INT n_comp;
  SHORT * u_comp, * d_comp;
  t_disc_list * disc_elem;
  np_part_discretization * disc;
  FVElementGeometry fvg;
  VECTOR * vec;
  param_connection * conn;
  param_predefined * p_val;

  /* Get access to the vector descriptors: */
  u_comp = VD_ncmp_cmpptr_of_otype (u, NODEVEC, &n_comp);
  if (n_comp != np->n_equations)
  {
    PrintErrorMessageF ('E', "AssembleDefect",
                        "Wrong solution vec. desc. or format: only %d comp. specified, but %d needed",
                        n_comp, np->n_equations);
    * res = __LINE__; return __LINE__;
  };
  d_comp = VD_ncmp_cmpptr_of_otype (d, NODEVEC, &n_comp);
  if (n_comp != np->n_equations)
  {
    PrintErrorMessageF ('E', "AssembleDefect",
                        "Wrong defect vec. desc. or format: only %d comp. specified, but %d needed",
                        n_comp, np->n_equations);
    * res = __LINE__; return __LINE__;
  };

  if (np->display_progress == PCR_FULL_DISPLAY)
    UserWriteF ("AssembleDefect: s_a = %g, s_m = %g ", s_a, s_m);

  /* Assemble the defect: */
  for (level = fl; level <= tl; level++)
  {
    if (np->display_progress >= PCR_RED_DISPLAY) UserWriteF ("[D:%d", level);
    if ((grid = GRID_ON_LEVEL (NP_MG (np), level)) == NULL)
    {
      * res = __LINE__; return __LINE__;
    };
    /* Loop the elements: */
    for (elem = FIRSTELEMENT (grid); elem != NULL; elem = SUCCE (elem))
    {
      /* Initialize the FV geometry: */
      if (EvaluateFVGeometry (elem, &fvg))
      {
        * res = __LINE__; return __LINE__;
      };
      /* Initialize the shapes and their derivatives: */
      if (EvaluateShapesAndDerivatives (&fvg, FILL_ALL))
      {
        * res = __LINE__; return __LINE__;
      };
      /* Write the local solution: */
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        for (co = 0; co < FVG_NSCV (&fvg); co++)
        {
          vec = NVECTOR (CORNER (FVG_ELEM (&fvg), co));
          for (i = 0; i < disc->n_equations; i++)
            disc->unknown (disc, co, i)
              = VVALUE (vec, u_comp [disc->first_vec_cmp + i]);
        };
      };
      /* Set the explicitly specified values of the import parameters: */
      for (p_val = np->value_list; p_val != NULL; p_val = p_val->tl)
        p_val->set (&fvg);
      /* Compute all the export parameters in the connections and
       * transfer the values through the connections:
       */
      for (conn = np->connection_list; conn != NULL; conn = conn->tl)
      {
        if (conn->exp_param->compute (&fvg, time))
        {
          PrintErrorMessageF ('E', "AssembleDefect",
                              "'compute' of export param. '%s' (%s) failed",
                              conn->exp_param->name, ENVITEM_NAME (conn->exp_param->disc));
          * res = __LINE__; return __LINE__;
        };
        conn->copy_values (&fvg);
      };
      /* Compute the local defects and add them to the global one: */
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        /* Compute the local defect: */
        if (disc->assemble_defect (disc, &fvg, time, s_m, s_a))
        {
          PrintErrorMessageF ('E', "AssembleDefect",
                              "assemble_defect of '%s' retuned error", ENVITEM_NAME (disc));
          * res = __LINE__; return __LINE__;
        };
        /* Add it to the global defect: */
        for (co = 0; co < FVG_NSCV (&fvg); co++)
        {
          vec = NVECTOR (CORNER (FVG_ELEM (&fvg), co));
          for (i = 0; i < disc->n_equations; i++)
            VVALUE (vec, d_comp [disc->first_vec_cmp + i])
              += disc->defect (disc, co, i) * disc_elem->scale;
        };
      };
    };
    /* Set the defect to zero at Dirichlet values: */
    for (vec = FIRSTVECTOR (grid); vec != NULL; vec = SUCCVC (vec))
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        for (i = 0; i < disc->n_equations; i++)
          if ((disc->Dirichlet_bc == NULL) ?
              (vec->skip & (1 << (disc->first_vec_cmp + i)))
              : disc->Dirichlet_bc (disc, vec->skip, i))
            VVALUE (vec, d_comp [disc->first_vec_cmp + i]) = 0;
      };
    if (np->display_progress >= PCR_RED_DISPLAY) UserWriteF ("]");
  };

  if (np->display_progress >= PCR_RED_DISPLAY)
    UserWriteF ("\n");

  if (np->print_defect)
  {
    NODE * node;
    UserWriteF ("***** Global defect: level %d, s_a = %g, s_m = %g: *****\n",
                tl, s_a, s_m);
    grid = GRID_ON_LEVEL (NP_MG (np), tl);
    for (node = FIRSTNODE (grid); node != NULL; node = SUCCN (node))
    {
      if (NVECTOR (node)->skip)
        UserWriteF ("%d (0%o): defect", ID (node), NVECTOR (node)->skip);
      else
        UserWriteF ("%d: defect", ID (node));
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        UserWriteF (" (%d:", disc->id);
        for (i = 0; i < disc->n_equations; i++)
          UserWriteF (" %g",
                      VVALUE (NVECTOR (node), d_comp [disc->first_vec_cmp + i]));
        UserWriteF (")");
      };
      UserWriteF (" & solution");
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        UserWriteF (" (%d:", disc->id);
        for (i = 0; i < disc->n_equations; i++)
          UserWriteF (" %g",
                      VVALUE (NVECTOR (node), u_comp [disc->first_vec_cmp + i]));
        UserWriteF (")");
      };
      UserWriteF (";\n");
    };
    UserWriteF ("***** End of the global defect *****\n");
  };

  * res = 0; return 0;
}

/** For debugging only: assembling the defect on one element for
** the numerical differentiation.
**/
#ifdef __GLOB_NUM_DIF__
#define ESOL(co,i) (sol [(co) * MAXNC + (i)])
#define EDEF(co,i) (def [(co) * MAXNC + (i)])
static INT AssembleElementDefect
(
  np_coupled_global_disc * np,
  ELEMENT * elem,
  DOUBLE time,
  DOUBLE s_m,
  DOUBLE s_a,
  DOUBLE * sol,
  DOUBLE * def
)
{
  np_part_discretization * disc;
  t_disc_list * disc_elem;
  param_connection * conn;
  param_predefined * p_val;
  FVElementGeometry fvg;
  INT co, i;

  if (EvaluateFVGeometry (elem, &fvg))
    return __LINE__;
  if (EvaluateShapesAndDerivatives (&fvg, FILL_ALL))
    return __LINE__;

  for (disc_elem = np->disc_list; disc_elem != NULL; disc_elem = disc_elem->tl)
  {
    disc = disc_elem->disc;
    for (co = 0; co < FVG_NSCV (&fvg); co++)
      for (i = 0; i < disc->n_equations; i++)
        disc->unknown (disc, co, i) = ESOL (co, disc->first_vec_cmp + i);
  };

  for (p_val = np->value_list; p_val != NULL; p_val = p_val->tl)
    p_val->set (&fvg);

  for (conn = np->connection_list; conn != NULL; conn = conn->tl)
  {
    if (conn->exp_param->compute (&fvg, time))
      return __LINE__;
    conn->copy_values (&fvg);
  };

  for (disc_elem = np->disc_list; disc_elem != NULL; disc_elem = disc_elem->tl)
  {
    disc = disc_elem->disc;
    if (disc->assemble_defect (disc, &fvg, time, s_m, s_a))
      return __LINE__;
    for (co = 0; co < FVG_NSCV (&fvg); co++)
      for (i = 0; i < disc->n_equations; i++)
        EDEF (co, disc->first_vec_cmp + i) = disc->defect (disc, co, i);
  };

  return 0;
}
#endif
/**/

/* AssemblePrecond - assembles the preconditioner matrix for the jacobian
 * (if all the derivatives are correct than - the jacobian). The function
 * return 0 if OK, nonzero on an error:
 */
static INT AssemblePrecond (NP_T_ASSEMBLE * ass, INT fl, INT tl, DOUBLE time,
                            DOUBLE s_a, VECDATA_DESC * u, VECDATA_DESC * d, VECDATA_DESC * v,
                            MATDATA_DESC * J, INT * res) /* Remark: s_m = 1 */
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) ass;
  INT level;
  GRID * grid;
  ELEMENT * elem;
  INT i, j, co_1, co_2;
  INT n_i_comp, n_j_comp;
  INT k, l;
  SHORT * u_comp, * J_comp, ind;
  SHORT * J_diag_comp, * J_entry_cmp; /* really used only for sparse blocks */
  t_disc_list * disc_elem;
  np_part_discretization * disc;
  import_param * imp_param;
  export_param * exp_param;
  FVElementGeometry fvg;
  VECTOR * vec;
  MATRIX * mat;
  param_connection * conn;
  param_predefined * p_val;
  import_param_list * import_elem;
  DOUBLE s;
  SPARSE_MATRIX * sm, * smd; /* used only for sparse blocks */

# ifdef __AUTO_SPARSE_DETECT__
# define ASSEMBLE_SPARSE MD_IS_SPARSE (J)
# else
# define ASSEMBLE_SPARSE (np->assemble_sparse)
# endif

  /** For debugging only: numerical differentiation: **/
# ifdef __GLOB_NUM_DIF__
# define EORIGSOL(co,i) (orig_sol [(co) * MAXNC + (i)])
# define EORIGDEF(co,i) (orig_def [(co) * MAXNC + (i)])
# define GLOB_NUM_DIF_H 1e-7
  DOUBLE orig_sol [MAXNC * 10];
  DOUBLE orig_def [MAXNC * 10];
  DOUBLE sol [MAXNC * 10];
  DOUBLE def [MAXNC * 10];

  if (np->n_equations > 10)
    return __LINE__;
# endif
  /**/

  /* Get access to the vector and matrix descriptors: */
  u_comp = VD_ncmp_cmpptr_of_otype (u, NODEVEC, &n_i_comp);
  if (n_i_comp != np->n_equations)
  {
    PrintErrorMessageF ('E', "AssemblePrecond",
                        "Wrong solution vec. desc. or format: %d comp. specified, but %d needed",
                        n_i_comp, np->n_equations);
    * res = __LINE__; return __LINE__;
  };

  /* Choose the correct matrix module and verify the matrix: */
  if (ASSEMBLE_SPARSE)
  {
    J_comp = MD_MCMPPTR_OF_MTYPE (J, MTP (0, 0));
    sm = MD_SM (J, MTP (0, 0));
    J_diag_comp = MD_MCMPPTR_OF_MTYPE (J, DMTP (0));
    smd = MD_SM (J, DMTP (0));
    np->verify_sparsity (J);
    if (np->display_progress == PCR_FULL_DISPLAY)
      UserWriteF ("AssemblePrecond (sparse fmt): s_a = %g ", s_a);
  }
  else
  {
    J_diag_comp = J_comp = MD_nr_nc_mcmpptr_of_ro_co
                             (J, NODEVEC, NODEVEC, &k, &n_j_comp);
    smd = sm = NULL;
    if (k != n_i_comp || n_j_comp != n_i_comp)
    {
      PrintErrorMessageF ('E', "AssemblePrecond",
                          "Wrong matrix type: point-block size %dx%d instead of %dx$d",
                          k, n_j_comp, n_i_comp, n_i_comp);
      * res = __LINE__; return __LINE__;
    };
    if (np->display_progress == PCR_FULL_DISPLAY)
      UserWriteF ("AssemblePrecond: s_a = %g ", s_a);
  };

  /* Set the matrix entries to zero: */
  if (dmatclear (NP_MG (np), fl, tl, ALL_VECTORS, J) != NUM_OK)
  {
    * res = __LINE__; return __LINE__;
  };

  /* Loop the elements: */
  for (level = fl; level <= tl; level++)
  {
    if (np->display_progress >= PCR_RED_DISPLAY) UserWriteF ("[M:%d", level);
    if ((grid = GRID_ON_LEVEL (NP_MG (np), level)) == NULL)
    {
      * res = __LINE__; return __LINE__;
    };
    for (elem = FIRSTELEMENT (grid); elem != NULL; elem = SUCCE (elem))
    {
      /** For debugging only: numerical differentiation: **/
#     ifdef __GLOB_NUM_DIF__
      for (co_1 = 0; co_1 < CORNERS_OF_ELEM (elem); co_1++)
        for (i = 0; i < np->n_equations; i++)
          EORIGSOL (co_1, i) = ESOL (co_1, i)
                                 = VVALUE (NVECTOR (CORNER (elem, co_1)), u_comp [i]);
      if (AssembleElementDefect (np, elem, time, 1, s_a, orig_sol, orig_def))
        return __LINE__;
      for (co_2 = 0; co_2 < CORNERS_OF_ELEM (elem); co_2++)
        for (j = 0; j < np->n_equations; j++)
        {
          ESOL (co_2, j) += GLOB_NUM_DIF_H;
          if (AssembleElementDefect (np, elem, time, 1, s_a, sol, def))
            return __LINE__;
          for (co_1 = 0; co_1 < CORNERS_OF_ELEM (elem); co_1++)
          {
            mat = GetMatrix (NVECTOR (CORNER (elem, co_1)),
                             NVECTOR (CORNER (elem, co_2)));
            J_entry_cmp = (MDIAG (mat)) ? J_diag_comp : J_comp;
            for (i = 0; i < np->n_equations; i++)
              if ((ind = J_entry_cmp [i * np->n_equations + j]) >= 0)
                MVALUE (mat, ind) +=
                  (EDEF (co_1, i) - EORIGDEF (co_1, i)) / GLOB_NUM_DIF_H;
          };
          for (co_1 = 0; co_1 < CORNERS_OF_ELEM (elem); co_1++)
            for (i = 0; i < np->n_equations; i++)
              ESOL (co_1, i) = EORIGSOL (co_1, i);
        };
#     else
      /**/
      /* Initialize the FV geometry: */
      if (EvaluateFVGeometry (elem, &fvg))
      {
        * res = __LINE__; return __LINE__;
      };
      /* Initialize the shapes and their derivatives: */
      if (EvaluateShapesAndDerivatives (&fvg, FILL_ALL))
      {
        * res = __LINE__; return __LINE__;
      };
      /* Write the local solution: */
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        for (co_1 = 0; co_1 < FVG_NSCV (&fvg); co_1++)
        {
          vec = NVECTOR (CORNER (FVG_ELEM (&fvg), co_1));
          for (i = 0; i < disc->n_equations; i++)
            disc->unknown (disc, co_1, i)
              = VVALUE (vec, u_comp [disc->first_vec_cmp + i]);
        };
      };
      /* Set the explicitly specified values of the import parameters: */
      for (p_val = np->value_list; p_val != NULL; p_val = p_val->tl)
        p_val->set (&fvg);
      /* Compute all the export parameters in the connections and
       * transfer the values through the connections:
       */
      for (conn = np->connection_list; conn != NULL; conn = conn->tl)
      {
        if (conn->exp_param->compute (&fvg, time))
        {
          PrintErrorMessageF ('E', "AssemblePrecond",
                              "'compute' of export param. '%s' (%s) failed",
                              conn->exp_param->name, ENVITEM_NAME (conn->exp_param->disc));
          * res = __LINE__; return __LINE__;
        };
        conn->copy_values (&fvg);
      };
      /* Assemble the jacobians of the discretizations: */
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        if (disc->assemble_jacobian (disc, &fvg, time, s_a))
        {
          PrintErrorMessageF ('E', "AssemblePrecond",
                              "assemble_jacobian of '%s' returned error", ENVITEM_NAME (disc));
          * res = __LINE__; return __LINE__;
        };
        for (co_1 = 0; co_1 < FVG_NSCV (&fvg); co_1++)
          for (co_2 = 0; co_2 < FVG_NSCV (&fvg); co_2++)
          {
            mat = GetMatrix (NVECTOR (CORNER (FVG_ELEM (&fvg), co_1)),
                             NVECTOR (CORNER (FVG_ELEM (&fvg), co_2)));
            J_entry_cmp = (MDIAG (mat)) ? J_diag_comp : J_comp;
            n_i_comp = disc->first_vec_cmp;
            for (i = 0; i < disc->n_equations; i++)
              for (j = 0; j < disc->n_equations; j++)
                if ((ind = J_entry_cmp
                           [(n_i_comp + i) * np->n_equations + (n_i_comp + j)]) >= 0)
                  MVALUE (mat, ind) += disc->jacobian (disc, co_1, i, co_2, j);
          };
      };
      /* Assemble the other derivatives: */
      for (conn = np->connection_list; conn != NULL; conn = conn->tl)
        if (conn->deriv)
        {
          /* Differentiate the export param. (w.r.t. unknowns of its disc.): */
          exp_param = conn->exp_param;
          if (exp_param->compute_D_param_wrt_unk (&fvg, time))
          {
            PrintErrorMessageF ('E', "AssemblePrecond",
                                "compute_D_param_wrt_unk of export param. '%s' (%s) returned error",
                                exp_param->name, ENVITEM_NAME (exp_param->disc));
            * res = __LINE__; return __LINE__;
          };
          n_i_comp = n_position_points (&fvg, exp_param->p_type);
          /* Follow the connection, compute the derivatives of the defect: */
          for (import_elem = conn->import_list; import_elem != NULL;
               import_elem = import_elem->tl)
            if (import_elem->deriv)
            {
              imp_param = import_elem->import;
              disc = imp_param->disc;
              n_j_comp = disc->first_vec_cmp;
              if (imp_param->compute_D_defect (&fvg, time, s_a))
              {
                PrintErrorMessageF ('E', "AssemblePrecond",
                                    "compute_D_defect of import param. '%s' (%s) retuned error",
                                    imp_param->name, ENVITEM_NAME (disc));
                * res = __LINE__; return __LINE__;
              };
              for (co_1 = 0; co_1 < FVG_NSCV (&fvg); co_1++)
                for (co_2 = 0; co_2 < FVG_NSCV (&fvg); co_2++)
                {
                  mat = GetMatrix (NVECTOR (CORNER (FVG_ELEM (&fvg), co_1)),
                                   NVECTOR (CORNER (FVG_ELEM (&fvg), co_2)));
                  J_entry_cmp = (MDIAG (mat)) ? J_diag_comp : J_comp;
                  for (i = 0; i < disc->n_equations; i++)
                    for (j = 0; j < exp_param->disc->n_equations; j++)
                      if ((ind = J_entry_cmp [(n_j_comp + i) * np->n_equations
                                              + (exp_param->disc->first_vec_cmp + j)]) >= 0)
                      {
                        s = 0;
                        for (k = 0; k < n_i_comp; k++)
                          for (l = 0; l < exp_param->n_comp; l++)
                            s += imp_param->D_defect (co_1, i, k, l)
                                 * exp_param->D_param_wrt_unk (k, l, co_2, j);
                        MVALUE (mat, ind) += s;
                      };
                };
            };
        };
#     endif /** __GLOB_NUM_DIF__ - a macro for debugging purposes **/
    };
    /* Assemble the Dirichlet boundary on this grid level;
     * otherwise scale the matrix:
     */
    for (vec = FIRSTVECTOR (grid); vec != NULL; vec = SUCCVC (vec))
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        for (i = 0; i < disc->n_equations; i++)
        {
          n_i_comp = disc->first_vec_cmp + i;
          if ((disc->Dirichlet_bc == NULL) ?
              (vec->skip & (1 << n_i_comp))
              : disc->Dirichlet_bc (disc, vec->skip, i))
          /* Set the identity matrix for this component: */
          {
            if (ASSEMBLE_SPARSE)
            {
              /* The offdiagonal blocks: */
              for (mat = MNEXT (VSTART (vec)); mat != NULL; mat = MNEXT (mat))
                for (j = sm->row_start [n_i_comp];
                     j < sm->row_start [n_i_comp + 1]; j++)
                  MVALUE (mat, sm->offset [j]) = 0;

              /* The diagonal block: */
              mat = VSTART (vec);
              for (j = smd->row_start [n_i_comp];
                   j < smd->row_start [n_i_comp + 1]; j++)
                MVALUE (mat, smd->offset [j]) = (smd->col_ind [j] == n_i_comp) ?
                                                1 : 0;
            }
            else
            {
              for (mat = VSTART (vec); mat != NULL; mat = MNEXT (mat))
                for (j = 0; j < np->n_equations; j++)
                  MVALUE (mat, J_comp [np->n_equations * n_i_comp + j]) = 0;
              MVALUE (VSTART (vec), J_comp [np->n_equations * n_i_comp + n_i_comp]) = 1;
            };
          }
          else if (disc_elem->scale_flag)
          {
            s = disc_elem->scale;
            if (ASSEMBLE_SPARSE)
            {
              /* The offdiagonal blocks: */
              for (mat = MNEXT (VSTART (vec)); mat != NULL; mat = MNEXT (mat))
                for (j = sm->row_start [n_i_comp];
                     j < sm->row_start [n_i_comp + 1]; j++)
                  MVALUE (mat, sm->offset [j]) *= s;

              /* The diagonal block: */
              mat = VSTART (vec);
              for (j = smd->row_start [n_i_comp];
                   j < smd->row_start [n_i_comp + 1]; j++)
                MVALUE (mat, smd->offset [j]) *= s;
            }
            else
              for (mat = VSTART (vec); mat != NULL; mat = MNEXT (mat))
                for (j = 0; j < np->n_equations; j++)
                  MVALUE (mat, J_comp [np->n_equations * n_i_comp + j]) *= s;
          };
        };
      };

    /* The (reduced) output: */
    if (np->display_progress >= PCR_RED_DISPLAY) UserWriteF ("]");
  };

  if (np->display_progress >= PCR_RED_DISPLAY)
    UserWriteF ("\n");

  /*** The assembling is done, below printing only: ***/
  if (np->print_matrix)
  {
    NODE * node;
    UserWriteF ("***** Global matrix: level %d, s_a = %g *****\n", tl, s_a);
    grid = GRID_ON_LEVEL (NP_MG (np), tl);
    for (node = FIRSTNODE (grid); node != NULL; node = SUCCN (node))
    {
      UserWriteF ("%d (vec %d):", ID (node), VINDEX (NVECTOR (node)));
      for (disc_elem = np->disc_list; disc_elem != NULL;
           disc_elem = disc_elem->tl)
      {
        disc = disc_elem->disc;
        UserWriteF (" (%d:", disc->id);
        for (mat = VSTART (NVECTOR (node)); mat != NULL; mat = MNEXT (mat))
        {
          J_entry_cmp = (MDIAG (mat)) ? J_diag_comp : J_comp;
          if (disc->n_equations == 1)
          {
            if (fabs (s = MVALUE (mat,
                                  J_entry_cmp [disc->first_vec_cmp * np->n_equations
                                               + disc->first_vec_cmp])) < 1e-12)
              continue;
            UserWriteF (" %g", s);
          }
          else
          {
            UserWriteF (" {");
            for (i = 0; i < disc->n_equations; i++)
            {
              n_i_comp = disc->first_vec_cmp + i;
              UserWriteF ("[");
              for (j = 0; j < np->n_equations; j++)
                if ((ind = J_entry_cmp [n_i_comp * np->n_equations + j]) >= 0)
                  UserWriteF (" %g", MVALUE (mat, ind));
                else
                  UserWriteF (" @");
              UserWriteF ("]");
            };
            UserWriteF ("}");
          };
          UserWriteF (" to vec %d", VINDEX (MDEST (mat)));
        };
        UserWriteF (");");
      };
      UserWriteF ("\n");
    };
    UserWriteF ("***** End of the global matrix *****\n");
  };

  * res = 0; return 0;
# undef ASSEMBLE_SPARSE
}

/*** Numproc constructor: ***/
static INT Construct (NP_BASE * base)
{
  np_coupled_global_disc * np = (np_coupled_global_disc *) base;
  NP_T_ASSEMBLE * ass = (NP_T_ASSEMBLE *) base;

  /* Construct NP_BASE: */
  base->Init = Init;
  base->Display = Display;
  base->Execute = NPTAssembleExecute;

  /* Construct NP_T_ASSEMBLE: */
  ass->TAssemblePreProcess = PreProcess;
  ass->TAssembleInitial = AssembleInitial;
  ass->TAssembleSolution = AssembleSolution;
  ass->TAssembleDefect = AssembleDefect;
  ass->TAssembleMatrix = AssemblePrecond;
  ass->TNAssembleMatrix = NULL;
  ass->TAssemblePostProcess = PostProcess;
  ass->TAssembleFinal = NULL;

  /* Construct np_coupled_global_disc: */
  np->disc_list = NULL;
  np->connection_list = NULL;
  np->value_list = NULL;
  np->n_equations = 0;

  return 0;
}

/***** The transfer numproc class: *****/

struct np_coupled_transfer : public np_transfer
{
  np_coupled_global_disc * global_disc; /* attached discretization */
};

/* The transfer preprocess: */
static INT transfer_PreProcess
(
  struct np_transfer * np, /* 'this' pointer */
  INT * fl, /* from level */
  INT tl, /* to level */
  VECDATA_DESC * s, /* solution vector */
  VECDATA_DESC * d, /* defect vector */
  MATDATA_DESC * A, /* system matrix */
  INT * res /* result */
)
{
  np_coupled_transfer * trans = (np_coupled_transfer *) np;
  t_disc_list * disc_elem;
  np_part_discretization * disc;

  if (trans->global_disc == NULL)
  {
    * res = __LINE__; return __LINE__;
  };

  #ifdef ModelP
  if (a_vector_vecskip (NP_MG (np), *fl, tl, s) != NUM_OK)
  {
    * res = __LINE__; return __LINE__;
  };
  #endif

  for (disc_elem = trans->global_disc->disc_list; disc_elem != NULL;
       disc_elem = disc_elem->tl)
    if ((disc = disc_elem->disc)->preprocess_transgrid != NULL)
      if (disc->preprocess_transgrid (disc, *fl, tl, s, d, A))
      {
        * res = __LINE__; return __LINE__;
      };

  * res = 0; return 0;
}

/* The interpolation of the correction: */
static INT transfer_InterpolateCorrection
(
  struct np_transfer * np, /* 'this' pointer */
  INT level, /* grid level */
  VECDATA_DESC * dst, /* destination vector */
  VECDATA_DESC * src, /* source vector */
  MATDATA_DESC * A, /* system matrix */
  VEC_SCALAR damp, /* damping factors */
  INT * res
)
{
  np_coupled_transfer * trans = (np_coupled_transfer *) np;
  t_disc_list * disc_elem;
  np_part_discretization * disc;

  if (trans->global_disc == NULL)
  {
    * res = __LINE__; return __LINE__;
  };

  #ifdef ModelP
  if (l_ghostvector_consistent (GRID_ON_LEVEL (NP_MG (np), level - 1), src)
      != NUM_OK)
  {
    * res = __LINE__; return __LINE__;
  };
  #endif

  for (disc_elem = trans->global_disc->disc_list; disc_elem != NULL;
       disc_elem = disc_elem->tl)
    if ((disc = disc_elem->disc)->interpolate_correction != NULL)
      if (disc->interpolate_correction (disc, level, dst, src, A, damp))
      {
        * res = __LINE__; return __LINE__;
      };

  * res = 0; return 0;
}

/* The restriction of the defect: */
static INT transfer_RestrictDefect
(
  struct np_transfer * np, /* 'this' pointer */
  INT level, /* grid level */
  VECDATA_DESC * dst, /* destination vector */
  VECDATA_DESC * src, /* source vector */
  MATDATA_DESC * A, /* system matrix */
  VEC_SCALAR damp, /* damping factors */
  INT * res
)
{
  np_coupled_transfer * trans = (np_coupled_transfer *) np;
  t_disc_list * disc_elem;
  np_part_discretization * disc;

  if (trans->global_disc == NULL)
  {
    * res = __LINE__; return __LINE__;
  };

  for (disc_elem = trans->global_disc->disc_list; disc_elem != NULL;
       disc_elem = disc_elem->tl)
    if ((disc = disc_elem->disc)->restrict_defect != NULL)
      if (disc->restrict_defect (disc, level, dst, src, A, damp))
      {
        * res = __LINE__; return __LINE__;
      };

  #ifdef ModelP
  if (l_ghostvector_collect (GRID_ON_LEVEL (NP_MG (np), level - 1), dst)
      != NUM_OK)
  {
    * res = __LINE__; return __LINE__;
  };
  #endif

  * res = 0; return 0;
}

/* Interpolation of new vectors: the standard function */
static INT transfer_InterpolateNewVectors
(
  struct np_transfer * np, /* 'this' pointer */
  INT fl, /* from level */
  INT tl, /* to level */
  VECDATA_DESC * vec, /* the vector to interpolate */
  INT * res /* result */
)
{
  INT i;

  for (i = fl + 1; i <= tl; i++)
  {
    #ifdef ModelP
    if (l_ghostvector_consistent (GRID_ON_LEVEL (NP_MG (np), i - 1), vec)
        != NUM_OK)
    {
      * res = __LINE__; return __LINE__;
    };
    #endif
    if ((* res = StandardInterpolateNewVectors (GRID_ON_LEVEL (NP_MG (np), i),
                                                vec)) != 0)
      return __LINE__;
  };

  * res = 0; return 0;
}

/* Project the solution: the standard function */
static INT transfer_ProjectSolution
(
  struct np_transfer * np, /* 'this' pointer */
  INT fl, /* from level */
  INT tl, /* to level */
  VECDATA_DESC * vec, /* the vector to project */
  INT * res /* result */
)
{
  INT i;

  for (i = tl - 1; i >= fl; i--)
  {
    if ((* res = StandardProject (GRID_ON_LEVEL (NP_MG (np), i),
                                  vec, vec)) != 0)
      return __LINE__;
    #ifdef ModelP
    if (l_ghostvector_project (GRID_ON_LEVEL (NP_MG (np), i), vec) != NUM_OK)
    {
      * res = __LINE__; return __LINE__;
    };
    #endif
  };

  * res = 0; return 0;
}

/* The transfer postprocess: */
static INT transfer_PostProcess
(
  struct np_transfer * np, /* 'this' pointer */
  INT * fl, /* from level */
  INT tl, /* to level */
  VECDATA_DESC * s, /* solution vector */
  VECDATA_DESC * d, /* defect vector */
  MATDATA_DESC * A, /* system matrix */
  INT * res /* result */
)
{
  np_coupled_transfer * trans = (np_coupled_transfer *) np;
  t_disc_list * disc_elem;
  np_part_discretization * disc;

  if (trans->global_disc == NULL)
  {
    * res = __LINE__; return __LINE__;
  };

  for (disc_elem = trans->global_disc->disc_list; disc_elem != NULL;
       disc_elem = disc_elem->tl)
    if ((disc = disc_elem->disc)->postprocess_transgrid != NULL)
      if (disc->postprocess_transgrid (disc, *fl, tl, s, d, A))
      {
        * res = __LINE__; return __LINE__;
      };

  return 0;
}

/*** Transfer numproc functions: ***/

static INT transfer_Init (NP_BASE * base, INT argc, char * * argv)
{
  np_coupled_transfer * trans = (np_coupled_transfer *) base;

  if ((trans->global_disc = (np_coupled_global_disc *) ReadArgvNumProc
                              (NP_MG (base), "global", NP_COUPLED_GLOBAL_DISC_CLASS_NAME, argc, argv))
      == NULL)
  {
    PrintErrorMessageF ('E', "coupled_transfer",
                        "Specify the global discretization numproc ($global) of the class '%s'",
                        NP_COUPLED_GLOBAL_DISC_CLASS_NAME);
    return NP_NOT_ACTIVE;
  };

  return NPTransferInit ((NP_TRANSFER *) base, argc, argv);
}

static INT transfer_Display (NP_BASE * base)
{
  np_coupled_transfer * trans = (np_coupled_transfer *) base;

  UserWriteF ("Coupled assembling transfer numproc:\n");
  UserWriteF (" global assembling numproc: %s\n",
              ENVITEM_NAME (trans->global_disc));

  return NPTransferDisplay ((NP_TRANSFER *) base);
}

/*** The constuctor: ***/

static INT transfer_Construct (NP_BASE * base)
{
  NP_TRANSFER * trans = (NP_TRANSFER *) base;

  base->Init = transfer_Init;
  base->Display = transfer_Display;
  base->Execute = NPTransferExecute;

  trans->PreProcess = transfer_PreProcess;
  trans->PreProcessSolution = NULL;
  trans->PreProcessProject = NULL;
  trans->InterpolateCorrection = transfer_InterpolateCorrection;
  trans->RestrictDefect = transfer_RestrictDefect;
  trans->InterpolateNewVectors = transfer_InterpolateNewVectors;
  trans->ProjectSolution = transfer_ProjectSolution;
  trans->AdaptCorrection = NULL;
  trans->PostProcess = transfer_PostProcess;
  trans->PostProcessProject = NULL;
  trans->PostProcessSolution = NULL;

  return 0;
}

/***** Class installer: *****/
INT NS_DIM_PREFIX Install_GlobalDisc ()
{
  if (CreateClass (NP_COUPLED_GLOBAL_DISC_CLASS_NAME,
                   sizeof (np_coupled_global_disc), Construct))
    return __LINE__;

  if (CreateClass (TRANSFER_CLASS_NAME ".coupled_transfer",
                   sizeof (np_coupled_transfer), transfer_Construct))
    return __LINE__;

  return 0;
}

/* End of File */
