// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* partdisc.hh - header for discretizations of model equations
 * in the coupled models.
 * History:
 * Aug. 29, 2002 - created
 * D. Logashenko
 * Based on the ideas of S. Paxion and N. Simus.
 **
 * Notation:
 * The discretization should discretize the model equations
 * only on one grid element. The discretization possesses
 * unknowns assigned to it, requires import parameters and
 * provides export parameters. It is supposed here, that
 * all unknowns are nodal (refere to the corners of the element).
 * Unlike the unknowns the export and input parameters have
 * a position type. Position type says to which positions
 * (in the FV geometry) the values of the parameters belong.
 * This counts the values and serves the compartibility of the
 * parameters. (There are, for ex., nodal parameters, ip-parameters,
 * ...)
 ***REMARK: Export parameter should not depend explicitly import parameters.
 * They are supposed to depend (explicitly) on the unknowns only. Else
 * the derivatives of the export parameters w.r.t. the import parameters
 * are neglected.
 *      Nevertheless, if it is OK to neglect this derivatives, you
 * may implement such dependencies. But if an export parameter explicitly
 * depends on an import parameter than that import parameter should be
 * initialized before the export parameter is computed. Such import parameters
 * are called prerequisites of the export parameter. To make the global
 * discretization to follow this computation order, you must attach to
 * the export parameter an NULL-terminated array of pointers to its
 * prerequisites. (Cf. 'export_param::set_prerequisites'.) Note that no
 * cyclic dependencies are allowed.
 **
 * The numprocs based on the class declared here provide functions that
 * assemble the defect as well as the stiffness and the mass matrices on
 * one element only. The following approach is assumed here: the functions
 * 'assemble_...' compute these data and save somewere in the numproc
 * structure (i. e. these data are not directly written to the grid
 * functions and matrices). Then the global assembling numproc accesses
 * these data using the further functions supplied in this numproc and
 * collects them in the global defect and the global matrix.
 **
 * Exception: The function 'assemble_initial' should do everything itself,
 * including placing the value in the VECTORs. This allows you to decide
 * whether to set the initial conditions at all. Set 'assemble_initial' to
 * NULL to do nothing.
 **
 * When computing the jacobian of the discretization, values of import
 * parameters should be considered as constants not depending on unknowns.
 * The derivatives coming from this dependence are assembled by the
 * global assembling routine automatically.
 **
 * Several discretizations may share the same boundary condition id. Then
 * They get the same boundary data. But they have two different sets of
 * skip flags.
 *      The number of skip flags may differ from that of equations. The skip
 * flags are treated internally in the discretizations, so this does not
 * hinder the assembling routine (cf. the member function 'skip_flag'). Set
 * 'skip_flag' to NULL to treat the skip flags in the standard way.
 *      But this may happen for one discretization only. This discretization
 * should have the greatest id. If the specified number of the skip flags
 * is less than the number of the equations (variables), the set of the
 * skip flags is extended to the full one by replication of the last
 * (for this disc.) skip bit. It is allowed for the number of equations
 * to be greater than the maximum number of skip bits (sizeof (INT)). In this
 * case, the specified number of skip bits should suit this constrait, and
 * the last specified bit defines the skip bits for the other equations.
 *      The number of skip bits can be greater than the number of equations,
 * too. Perhaps, it can be useful at assembling phase.
 * Remark: The standard structure of the skip bits (skip bit #i set <=>
 * the component #i contains a Dirichlet value) is important for the
 * parallelization (for ex. for the computation of the diagonal matrix entries
 * at border nodes, etc). The additional skip bits are never checked.
 **
 * Standards for the boundary conditions:
 * The array 'in' of the boundary should consist of at least
 * MIN_BC_IN_SIZE entries. The first MIN_BC_IN_SIZE are used for
 * the system purposes:
 * in [0..DOM_N_IN_PARAMS-1]: cf. domain.h
 * in [DOM_N_IN_PARAMS] = the bc id of the discretization (converted to DOUBLE)
 * in [DOM_N_IN_PARAMS+1] = time
 * Other entries can be used for the discretization specific information.
 * The parameter 'in' can be set to NULL or omitted (by the
 * 'np_part_discretization::' bc interface functions). Then an array
 * of exactly MIN_BC_IN_SIZE is used.
 *      The 'type' and 'bc_value' arrays have a discretization specific
 * structure. The class does not make any assumptions on it. Typically,
 * the following structure is advised: 'type' is an array of n_equations
 * INTs that keep the bc type for every component, 'bc_value' is
 * an array of n_equations DOUBLEs that are the values of the values of
 * the bc.
 **
 * REMARKS for sparse matrices:
 * 1. To assemble sparse matrix blocks you should supply the sparsity pattern
 * for the jacobian. If you set it to NULL, a full matrix is supposed for
 * this discretization. This pattern is then checked to suit the one of the
 * sparse matrix blocks: the latter ones should be subsets of the former ones.
 * The sparsity patterns for the parameters are never checked. Nonzero values
 * outside of the patterns are merely neglected.
 * 2. Even OUTSIDE of the specified pattern you should CORRECTLY SET the
 * values of the Jacobian: e. g. set them to zero! These values are used
 * for example when full matrices are assembled. Any way, incorrect
 * initialization can lead to numerical errors.
 */

#ifndef __UG_CPL_PART_DISC__
#define __UG_CPL_PART_DISC__

/* UG headers: */
#include "namespace.h"
#include "compiler.h"

START_UGDIM_NAMESPACE

/* The number of the system entries of the input bc array: */
#define MIN_BC_IN_SIZE (DOM_N_IN_PARAMS + 2)

/* The position types: */
enum position_type
{
  POSITION_IP = 0, /* integration point */
  POSITION_BIP = 1, /* boundary integration point */
  POSITION_NODE = 2, /* corner of the element */
  POSITION_ELEMENT = 3, /* the element itself */
  POSITION_SIDE = 4, /* a side of the element */
  N_POSITION_TYPES
};

/* n_position_points - returns the number of position points for
 * the position types:
 */
INT n_position_points
(
  FVElementGeometry * fvg, /* specifies the element */
  position_type type /* the position type to check */
);

/* max_n_position_points - array that keeps the maximum number of
 * position points in a given position type:
 */
extern INT max_n_position_points [N_POSITION_TYPES];

/* Maximum number of position points in all position types: */
#define MAX_POS_PNTS MAX (MAXF, MAX (MAXBF, MAX (MAXNC, MAXS)))

struct np_part_discretization;

/*** The import parameter class: ***/
struct import_param
{
  /* Data provided by the user: */

  import_param * tl; /* next in the list of import params of this disc. */
  np_part_discretization * disc; /* own discretization */

  char * name; /* name of this parameter */
  position_type p_type; /* its position type */
  INT n_comp; /* number of components at every position point */

  /* Data set by the assembling procedure: */
  char used; /* nonzero iff the param. is initialized during assembling */

  /* Constructor: */

  import_param
  (
    import_param * the_tl, /* pointer to the next one in the list, or NULL */
    np_part_discretization * the_disc, /* own discretization */
    char * the_name, /* name of this parameter */
    position_type the_p_type, /* its position type */
    INT the_n_comp = 1 /* number of components at every position point */
  )
  {
    tl = the_tl; disc = the_disc; name = the_name;
    p_type = the_p_type; n_comp = the_n_comp;
  };

  /* Functions for computing the values: */

  virtual INT compute_D_defect /* computes the derivative of the defect */
  (
    FVElementGeometry * fvg, /* the FV geometry data */
    DOUBLE time, /* the time argument */
    DOUBLE s_a /* the factor for the stiffness matrix (s_m = 1) */
  );

  /* Functions to read/write values: */

  virtual DOUBLE & operator () /* to read/write the values of the parameter */
  (
    INT pos, /* index of the position point */
    INT comp /* index of the component */
  ) = 0;

  virtual DOUBLE D_defect /* to read the derivative of the defect */
  (
    INT co, /* corner of the element for the defect */
    INT d_comp, /* component of the defect */
    INT pos, /* position point of the parameter */
    INT p_comp /* component of the parameter */
  );

  /* Other functions: */

  virtual void Display (); /* to display the parameter */

  INT num_D_defect /* num. differentiation of defect */
  (
    FVElementGeometry * fvg,
    DOUBLE time,
    DOUBLE s_a, /* these three are like for compute_D_defect */
    DOUBLE * deriv, /* the array to save the derivative in */
    DOUBLE h = 1e-7 /* the step for the numerical integration */
  );
};

/*** The export parameter class: ***/
struct export_param
{
  export_param * tl; /* next in the list of export params of this disc. */
  np_part_discretization * disc; /* own discretization */

  char * name; /* name of this parameter */
  position_type p_type; /* its position type */
  INT n_comp; /* number of components at every position point */

  import_param * * prerequisite; /* NULL-term. array of prerequisites (or NULL) */

  /* Constructor: */

  export_param
  (
    export_param * the_tl, /* pointer to the next one in the list, or NULL */
    np_part_discretization * the_disc, /* own discretization */
    char * the_name, /* name of this parameter */
    position_type the_p_type, /* its position type */
    INT the_n_comp = 1 /* number of components at every position point */
  )
  {
    tl = the_tl; disc = the_disc; name = the_name;
    p_type = the_p_type; n_comp = the_n_comp;
    prerequisite = NULL;
  };

  /* Functions to compute values: */

  virtual INT compute /* compute the values */
  (
    FVElementGeometry * fvg, /* the FV geometry data */
    DOUBLE time /* the time argument */
  ) = 0;

  virtual INT compute_D_param_wrt_unk /* derive the param. wrt the unknowns */
  (
    FVElementGeometry * fvg, /* the FV geometry data */
    DOUBLE time /* the time argument */
  );

  /* Functions to read/write values: */

  virtual DOUBLE operator () /* to read the values of the param. */
  (
    INT pos, /* index of the position point */
    INT comp /* index of the component */
  ) = 0;

  virtual DOUBLE D_param_wrt_unk /* derivative of the param. wrt unknowns */
  (
    INT pos, /* index of the position point of the param */
    INT p_comp, /* index of the component of the param */
    INT co, /* corner to which the unknown is assigned */
    INT unk_comp /* index of the unknown */
  );

  /* Other functions: */

  void set_prerequisites /* sets the array of prerequisites */
  (
    import_param * * the_prerequisites
  )
  {
    prerequisite = the_prerequisites;
  };

  virtual void Display (); /* to display the parameter */

  INT num_D_param /* numerical differentiation of the parameter */
  (
    FVElementGeometry * fvg,
    DOUBLE time, /* these two are identical to those of compute_D_param_wrt_unk */
    DOUBLE * deriv, /* to save the derivative */
    DOUBLE h = 1e-7 /* the numerical differentiation step */
  );
};

/*** The class of the partial discretization: ***/
struct np_part_discretization : public np_base
{
  /* Fields set by a global discretization numproc (do not rewrite it): */
  INT id; /* id of the discretization */
  INT first_vec_cmp; /* the index of the first component in the VD */
  /* Remark: first_vec_cmp is the index of the first skip flag, too. */

  /* User-defined data: */
  INT n_equations; /* number of equation and unknowns */
  INT bc_id; /* the identifier for the boundary conditions */
  INT n_skip_flags; /* number of skip flags used by this discretization */

  /* Lists of the parameters: */
  import_param * import_list; /* the list of import parameters */
  export_param * export_list; /* the list of export parameters */

  /* Functions for assembling the values: */

  INT (* assemble_initial) /* assembles initial values */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT fl, /* from level */
    INT tl, /* to level */
    DOUBLE time, /* the time argument */
    VECDATA_DESC * u /* where to assemble the initial values */
  );

  INT (* preprocess) /* called during running of the TAssemblePreProcess */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT fl, /* from level */
    INT tl, /* to level */
    DOUBLE new_time, /* the new time level */
    DOUBLE cur_time, /* the current time level */
    DOUBLE old_time, /* the old time level */
    VECDATA_DESC * new_sol, /* vd for the new solution (still not init.) */
    VECDATA_DESC * cur_sol, /* vd of the current solution */
    VECDATA_DESC * old_sol /* vd of the old solution */
  );

  INT (* assemble_solution) /* assembles the boundary and the skip flags */
  (
    np_part_discretization * np, /* 'this' pointer */
    NODE * node, /* the grid node */
    DOUBLE time /* the time argument */
  );

  INT (* assemble_defect) /* assembles the local defect */
  (
    np_part_discretization * np, /* 'this' pointer */
    FVElementGeometry * fvg, /* FV geometry data */
    DOUBLE time, /* the time argument */
    DOUBLE s_m, /* the mass matrix factor */
    DOUBLE s_a /* the stiffness matrix factor */
  );

  INT (* assemble_jacobian) /* assembles the local jacobian */
  (
    np_part_discretization * np, /* 'this' pointer */
    FVElementGeometry * fvg, /* FV geometry data */
    DOUBLE time, /* the time argument */
    DOUBLE s_a /* the stiffness matrix factor (s_m = 1) */
  );

  INT (* postprocess) /* called during running of TAssemblePostProcess */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT fl, /* from level */
    INT tl, /* to level */
    DOUBLE new_time, /* the new time level */
    DOUBLE cur_time, /* the current time level */
    DOUBLE old_time, /* the old time level */
    VECDATA_DESC * new_sol, /* vd for the new solution (still not init.) */
    VECDATA_DESC * cur_sol, /* vd of the current solution */
    VECDATA_DESC * old_sol /* vd of the old solution */
  );

  /* Functions for reading/writing values: */

  void (* boundary_value) /* for init. of boundary values at a given node */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT comp, /* index of the unknown */
    DOUBLE * val /* to write the value to */
  );
  /* As the argument 'val', this function gets the pointer to the position
   * of the value in the corresponding VECTOR. As the function is called,
   * it can decide itself whether or not to set this value. There are no
   * other initializations of the nodal values.
   */

  INT (* skip_flag) /* for reading skip flags for a given node */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT i /* the (discretization's) index of the unknown */
  );
  /* This function should return nonzero to set the skip flag to 1
   * and zero to set it to zero. Other properties of the return value
   * are neglected.
   */

  DOUBLE & (* unknown) /* for writing/reading the local values of the unk. */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT co, /* the corner of the element */
    INT unk /* index of the unknown */
  );

  DOUBLE (* defect) /* returns value of the local defect */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT co, /* the corner of the element */
    INT comp /* index of the component */
  );

  DOUBLE (* jacobian) /* returns the entries of the local jacobian */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT co_1, /* index of the first corner */
    INT unk_1, /* index of the first unknown */
    INT co_2, /* index of the second corner */
    INT unk_2 /* index of the second unknown */
  );

  INT (* Dirichlet_bc) /* returns nonzero for Dirichle BC components */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT skip_flags, /* all skip flags from a VECTOR */
    INT comp /* index of the component */
  );
  /* Remark: 'Dirichle_bc' should interpret the skip flags into the
   * Dirichle BC positions. This function should use the
   * np_part_discretization's fields 'first_skip_flag' and 'n_skip_flags'
   * and return nonzero only for the positions that should contain the
   * Dirichlet values in the grid function. You may set 'dirichlet_bc' to
   * NULL, then the skip flags are treated in the standard way: the i-th
   * flag set iff the i-th component is a Dirichlet BC.
   */

  /*** Sparsity patterns for matrices (if sparse matrices are used): ***/

  /* This function should return nonzero for the entries that
   * can be nonzero and zero for the entries that are always zero.
   * If you specify NULL then nonzero value is supposed for all positions:
   */

  INT (* jac_sparse_pattern) /* the pattern of the blocks of the jacobian */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT diag, /* nonzero for the diagonal blocks, 0 for offdiagonal ones */
    INT i, /* the index of the first unknown */
    INT j /* the index of the second unknown */
  );

  /*** The grid transfer function: projection and restriction ***/

  INT (* preprocess_transgrid) /* the preprocess */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT fl, /* from level */
    INT tl, /* to level */
    VECDATA_DESC * u, /* solution vector */
    VECDATA_DESC * d, /* defect vector */
    MATDATA_DESC * A /* system matrix */
  );

  INT (* interpolate_correction) /* interpolates the correction */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT level, /* grid level */
    VECDATA_DESC * dst, /* destination vector */
    VECDATA_DESC * src, /* source vector */
    MATDATA_DESC * A, /* system matrix */
    VEC_SCALAR damp /* damping factors */
  );

  INT (* restrict_defect) /* restricts the defect */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT level, /* grid level */
    VECDATA_DESC * dst, /* destination vector */
    VECDATA_DESC * src, /* source vector */
    MATDATA_DESC * A, /* system matrix */
    VEC_SCALAR damp /* damping factors */
  );

  INT (* postprocess_transgrid) /* the postprocess */
  (
    np_part_discretization * np, /* 'this' pointer */
    INT fl, /* from level */
    INT tl, /* to level */
    VECDATA_DESC * u, /* solution vector */
    VECDATA_DESC * d, /* defect vector */
    MATDATA_DESC * A /* system matrix */
  );

  /*** Functions to access boundary conditions: ***/

  INT bndp_cond /* analogue of BNDP_BndCond */
  (
    NODE * node, /* node where the bc should be computed */
    INT * n, /* to read the number of segments for this node */
    INT i, /* the index of the segment (0..n-1) */
    DOUBLE time, /* the time */
    INT * type, /* the type of the bc for every component
                 * (array of n_equations values)
                 */
    DOUBLE * bc_value, /* the values returned by the bc function
                        * (array of n_equations values)
                        */
    DOUBLE * in = NULL /* values to transfer to the bc function
                        * (cf. standards above) or NULL if none
                        */
  );

  INT bnds_cond /* analogue of BNDS_BndCond */
  (
    const ELEMENT * elem, /* the element for which the bc should be computed */
    INT side, /* index of the side of an element */
    DOUBLE * local, /* the local coordinates on the side (DIM - 1 values) */
    DOUBLE time, /* the time */
    INT * type, /* the type of the bc for every component
                 * (array of n_equations values)
                 */
    DOUBLE * bc_value, /* the values returned by the bc function
                        * (array of n_equations values)
                        */
    DOUBLE * in = NULL /* values to transfer to the bc function
                        * (cf. standards above) or NULL if none
                        */
  );

  /* Numerical differenciation: */

  INT num_jacobian /* the numerical Jacobian computed using the defect */
  (
    FVElementGeometry * fvg, /* FV geometry data */
    DOUBLE time, /* the time argument */
    DOUBLE s_a, /* the stiffness matrix factor (s_m = 1) */
    DOUBLE * jac, /* where to save the jacobian (cf. partdisc.cc) */
    DOUBLE h = 1e-7 /* the step of the numerical differenciation */
  );

  /* Auxiliary functions: */

  void display_import_export (); /* prints import and export parameters */

  /* The special allocation operator: */
  void * operator new
  (
    size_t s, /* size of the object: specified automatically */
    void * buffer /* the buffer to keep the object in */
  )
  {
    return buffer;
  };

};

/* Macro for adjusting the C++ object initialization with the UG engine.
 * This macro declares two functions: the UG constructor of the class
 * and the UG numproc installer. The first one is static. Declare the
 * installer in a header as
 * extern "C" INT installer_name (void);
 * The arguments of the macro are:
 * class_name - the name of your class derived from 'np_part_discretization',
 * np_name - the name of the numproc (in "..."),
 * installer_name - the installer name.
 */
#define UG_INSTALL_CLASS(class_name, np_name, installer_name) \
  static INT class_name ## _UG_constructor (NP_BASE * base) \
  { \
    new (base) class_name; \
    return 0; \
  } \
  INT installer_name () \
  { \
    if (CreateClass (np_name, sizeof (class_name), class_name ## _UG_constructor)) \
      return __LINE__; \
    return 0; \
  }

/* The standard functions for the transfer operators: */

INT cpl_standard_interpolate_correction /* interpolates the correction */
(
  np_part_discretization * np, /* 'this' pointer */
  INT level, /* grid level */
  VECDATA_DESC * dst, /* destination vector */
  VECDATA_DESC * src, /* source vector */
  MATDATA_DESC * A, /* system matrix */
  VEC_SCALAR damp /* damping factors */
);

INT cpl_standard_restrict_defect /* restricts the defect */
(
  np_part_discretization * np, /* 'this' pointer */
  INT level, /* grid level */
  VECDATA_DESC * dst, /* destination vector */
  VECDATA_DESC * src, /* source vector */
  MATDATA_DESC * A, /* system matrix */
  VEC_SCALAR damp /* damping factors */
);

/* The numproc class name: */
#define NP_PART_DISCRETIZATION_CLASS_NAME "partdisc"

END_UGDIM_NAMESPACE

#endif

/* End of File */
