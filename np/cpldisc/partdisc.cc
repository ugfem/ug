// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* partdisc.cc - common data and functions for the local
 * discretizations.
 * History:
 * Sep. 9, 2002 - created
 * D. Logashenko
 * Based on ideas of S. Paxion and N. Simus.
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
  #include "shapes.h"
  #include "misc.h"
  #include "np.h"
  #include "fvgeom.h"
}

/* Own header: */
#include "partdisc.hh"

/* n_position_points - returns the number of position points for
 * the position types:
 */
INT NS_DIM_PREFIX n_position_points
(
  FVElementGeometry * fvg, /* specifies the element */
  position_type type /* the position type to check */
)
{
  switch (type)
  {
  case POSITION_IP :
    return FVG_NSCVF (fvg);
  case POSITION_BIP :
    return FVG_NSCVBF (fvg);
  case POSITION_NODE :
    return FVG_NSCV (fvg);
  case POSITION_ELEMENT :
    return 1;
  case POSITION_SIDE :
    return SIDES_OF_ELEM (FVG_ELEM (fvg));
  };
  return 0;
}

/* max_n_position_points - array that keeps the maximum number of
 * position points in a given position type:
 */
INT NS_DIM_PREFIX max_n_position_points [N_POSITION_TYPES] =
{
  MAXF, MAXBF, MAXNC, 1, MAXS
};

/*** Dummy virtual functions for the derivatives: zero values ***/

INT import_param::compute_D_defect /* dummy substitution */
(
  FVElementGeometry * fvg, /* the FV geometry data */
  DOUBLE time, /* the time argument */
  DOUBLE s_a /* the factor for the stiffness matrix (s_m = 1) */
)
{
  return 0;
}

DOUBLE import_param::D_defect /* dummy substitution */
(
  INT co, /* corner of the element for the defect */
  INT d_comp, /* component of the defect */
  INT pos, /* position point of the parameter */
  INT p_comp /* component of the parameter */
)
{
  return 0;
}

INT export_param::compute_D_param_wrt_unk /* dummy substitution */
(
  FVElementGeometry * fvg, /* the FV geometry data */
  DOUBLE time /* the time argument */
)
{
  return 0;
}

DOUBLE export_param::D_param_wrt_unk /* dummy substitution */
(
  INT pos, /* index of the position point of the param */
  INT p_comp, /* index of the component of the param */
  INT co, /* corner to which the unknown is assigned */
  INT unk_comp /* index of the unknown */
)
{
  return 0;
}

/*** The standard transgrid functions: ***/

INT NS_DIM_PREFIX cpl_standard_interpolate_correction /* interpolates the correction */
(
  np_part_discretization * np, /* 'this' pointer */
  INT level, /* grid level */
  VECDATA_DESC * dst, /* destination vector */
  VECDATA_DESC * src, /* source vector */
  MATDATA_DESC * A, /* system matrix */
  VEC_SCALAR damp /* damping factors */
)
{
  INT i, k, n;
  SHORT * from_comp, * to_comp;
  GRID * fine_grid;
  NODE * node;
  VECTOR * vec, * c_vec;
  VERTEX * vert;
  ELEMENT * elem;
  DOUBLE c [MAX_CORNERS_OF_ELEM];

  /* Get the components of the vector descriptors: */
  from_comp = VD_ncmp_cmpptr_of_otype (src, NODEVEC, &n);
  if (n < np->n_equations)
    return __LINE__;
  to_comp = VD_ncmp_cmpptr_of_otype (dst, NODEVEC, &n);
  if (n < np->n_equations)
    return __LINE__;

  /* Get the fine grid: */
  if (level < 1) return __LINE__;
  fine_grid = GRID_ON_LEVEL (NP_MG (np), level);

  /* Prolong the correction: */
  for (node = FIRSTNODE (fine_grid); node != NULL; node = SUCCN (node))
  {
    vec = NVECTOR (node);

    /* Set the components to zero: */
    for (k = 0; k < np->n_equations; k++)
      VVALUE (vec, to_comp [np->first_vec_cmp + k]) = 0;

    if (CORNERTYPE (node))
    {
      /* Copy the components of the vector: */
      c_vec = NVECTOR ((NODE *) NFATHER (node));
      for (k = 0; k < np->n_equations; k++)
        if ((np->Dirichlet_bc == NULL) ?
            (vec->skip & (1 << (np->first_vec_cmp + k))) == 0
            : np->Dirichlet_bc (np, vec->skip, k) == 0)
          VVALUE (vec, to_comp [np->first_vec_cmp + k])
            = damp [np->first_vec_cmp + k]
              * VVALUE (c_vec, from_comp [np->first_vec_cmp + k]);
    }
    else
    {
      /* Sum the 'hat'-functions: */
      vert = MYVERTEX (node);
      elem = VFATHER (vert);
      GNs ((n = CORNERS_OF_ELEM (elem)), LCVECT (vert), c);
      for (i = 0; i < n; i++)
      {
        c_vec = NVECTOR (CORNER (elem, i));
        for (k = 0; k < np->n_equations; k++)
          if ((np->Dirichlet_bc == NULL) ?
              (vec->skip & (1 << (np->first_vec_cmp + k))) == 0
              : np->Dirichlet_bc (np, vec->skip, k) == 0)
            VVALUE (vec, to_comp [np->first_vec_cmp + k])
              += damp [np->first_vec_cmp + k] * c [i]
                 * VVALUE (c_vec, from_comp [np->first_vec_cmp + k]);
      };
    };
  };

  return 0;
}

INT NS_DIM_PREFIX cpl_standard_restrict_defect /* restricts the defect */
(
  np_part_discretization * np, /* 'this' pointer */
  INT level, /* grid level */
  VECDATA_DESC * dst, /* destination vector */
  VECDATA_DESC * src, /* source vector */
  MATDATA_DESC * A, /* system matrix */
  VEC_SCALAR damp /* damping factors */
)
{
  INT i, k, n;
  SHORT * from_comp, * to_comp;
  GRID * fine_grid;
  NODE * node;
  VECTOR * vec, * c_vec;
  VERTEX * vert;
  ELEMENT * elem;
  DOUBLE c [MAX_CORNERS_OF_ELEM];

  /* Get the components of the vector descriptors: */
  from_comp = VD_ncmp_cmpptr_of_otype (src, NODEVEC, &n);
  if (n < np->n_equations)
    return __LINE__;
  to_comp = VD_ncmp_cmpptr_of_otype (dst, NODEVEC, &n);
  if (n < np->n_equations)
    return __LINE__;

  /* Get the fine grid: */
  if (level < 1) return __LINE__;
  fine_grid = GRID_ON_LEVEL (NP_MG (np), level);

  /* Set the components to zero: */
  for (vec = PFIRSTVECTOR (GRID_ON_LEVEL (NP_MG (np), level - 1)); vec != NULL;
       vec = SUCCVC (vec))
    if (VNCLASS (vec) >= NEWDEF_CLASS)
      for (k = 0; k < np->n_equations; k++)
        VVALUE (vec, to_comp [np->first_vec_cmp + k]) = 0;

  /* Prolong the correction: */
  for (node = FIRSTNODE (fine_grid); node != NULL; node = SUCCN (node))
    if (VCLASS (vec = NVECTOR (node)) >= NEWDEF_CLASS)
    {
      /* Sum the 'hat'-functions: */
      if (CORNERTYPE (node))
      {
        c_vec = NVECTOR ((NODE *) NFATHER (node));
        for (k = 0; k < np->n_equations; k++)
          if ((np->Dirichlet_bc == NULL) ?
              (c_vec->skip & (1 << (np->first_vec_cmp + k))) == 0
              : np->Dirichlet_bc (np, c_vec->skip, k) == 0)
            VVALUE (c_vec, to_comp [np->first_vec_cmp + k])
              += damp [np->first_vec_cmp + k]
                 * VVALUE (vec, from_comp [np->first_vec_cmp + k]);
      }
      else
      {
        vert = MYVERTEX (node);
        elem = VFATHER (vert);
        GNs ((n = CORNERS_OF_ELEM (elem)), LCVECT (vert), c);
        for (i = 0; i < n; i++)
        {
          c_vec = NVECTOR (CORNER (elem, i));
          for (k = 0; k < np->n_equations; k++)
            if ((np->Dirichlet_bc == NULL) ?
                (c_vec->skip & (1 << (np->first_vec_cmp + k))) == 0
                : np->Dirichlet_bc (np, c_vec->skip, k) == 0)
              VVALUE (c_vec, to_comp [np->first_vec_cmp + k])
                += damp [np->first_vec_cmp + k] * c [i]
                   * VVALUE (vec, from_comp [np->first_vec_cmp + k]);
        };
      };
    };

  return 0;
}

/*** Implementation of the boundary condition interface: ***/

/* bndp_cond - accesses the boundary conditions at a point.
 * The function returns 0 if OK, nonzero on an error:
 */
INT np_part_discretization::bndp_cond /* analogue of BNDP_BndCond */
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
  DOUBLE * in /* values to transfer to the bc function
               * (cf. standards above) or NULL if none
               */
)
{
  DOUBLE loc_in [MIN_BC_IN_SIZE];

  if (in == NULL) in = loc_in;
  in [DOM_N_IN_PARAMS] = bc_id;
  in [DOM_N_IN_PARAMS + 1] = time;

  return BNDP_BndCond (V_BNDP (MYVERTEX (node)), n, i, in, bc_value, type);
}

/* bnds_cond - accesses the boundary condition on an element side.
 * The function returns 0 if OK, nonzero on an error:
 */
INT np_part_discretization::bnds_cond /* analogue of BNDS_BndCond */
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
  DOUBLE * in /* values to transfer to the bc function
               * (cf. standards above) or NULL if none
               */
)
{
  DOUBLE loc_in [MIN_BC_IN_SIZE];

  if (in == NULL) in = loc_in;
  in [DOM_N_IN_PARAMS] = bc_id;
  in [DOM_N_IN_PARAMS + 1] = time;

  return BNDS_BndCond (ELEM_BNDS (elem, side), local, in, bc_value, type);
}

/** Implementations for printing information about import/export parameters: **/

static const char * position_type_name [N_POSITION_TYPES] =
{
  "integration point", "boundary integration point",
  "node", "element", "side"
};

void import_param::Display ()
{
  UserWriteF (" import parameter '%s': %d cmp. at every %s\n",
              name, n_comp, position_type_name [p_type]);
}

void export_param::Display ()
{
  UserWriteF (" export parameter '%s': %d cmp. at every %s\n",
              name, n_comp, position_type_name [p_type]);
}

static char none_str [] = " -- none --\n";

void np_part_discretization::display_import_export ()
{
  import_param * imp_param;
  export_param * exp_param;

  UserWrite ("Import parameters:\n");
  if (import_list == NULL)
    UserWrite (none_str);
  else
    for (imp_param = import_list; imp_param != NULL; imp_param = imp_param->tl)
      imp_param->Display ();

  UserWrite ("Export parameters:\n");
  if (export_list == NULL)
    UserWrite (none_str);
  else
    for (exp_param = export_list; exp_param != NULL; exp_param = exp_param->tl)
      exp_param->Display ();
}

/*** Implementation of the numerical differentiation: ***/

/* num_jacobian - computes the numerical Jacobian of the defect w.r.t. the
 * unknowns. The function is assumed to be called from the 'assemble_jacobian'
 * function and needs further processing of the local jacobian entries.
 * THE FUNCTION WORKS WITH NO MORE THAN 10 EQUATIONS IN THE DISCRETIZATION!
 **
 * The function get the array for the jacobian as an argument (the array
 * 'jac'). The array should have the size MAXNC * <number of equations>
 * * MAXNC * <number of equations>. For the structure of the array s. the
 * mactor JAC in the function. (For ex., if you have only one unknown in
 * the discretization, then 'jac' has type 'DOUBLE [MAXNC] [MAXNC]', and
 * jac [i] [j] denotes the entry (i, j) of the local jacobian.)
 **
 * The function returns 0 if OK, nonzero on an error:
 */
INT np_part_discretization::num_jacobian
(
  FVElementGeometry * fvg, /* FV geometry data */
  DOUBLE time, /* the time argument */
  DOUBLE s_a, /* the stiffness matrix factor (s_m = 1) */
  DOUBLE * jac, /* where to save the jacobian */
  DOUBLE h /* the step of the numerical differenciation */
)
{
  DOUBLE sol [MAXNC * 10];
  DOUBLE def [MAXNC * 10];
  INT n_co = FVG_NSCV (fvg);
  INT i_eq, j_eq, i_co, j_co;

  if (n_equations > 10)
    return __LINE__;

# define SOL(co,eq) (sol [(co) * 10 + (eq)])
# define DEF(co,eq) (def [(co) * 10 + (eq)])
# define JAC(f_co,i,t_co,j) (jac [(((f_co * n_equations) + i) * MAXNC + t_co) \
                                  * n_equations + j])

  /* Save the original solution: */
  for (i_co = 0; i_co < n_co; i_co++)
    for (i_eq = 0; i_eq < n_equations; i_eq++)
      SOL (i_co, i_eq) = unknown (this, i_co, i_eq);

  /* Compute the defect at the original point and save it: */
  if (assemble_defect (this, fvg, time, 1, s_a))
    return __LINE__;
  for (i_co = 0; i_co < n_co; i_co++)
    for (i_eq = 0; i_eq < n_equations; i_eq++)
      DEF (i_co, i_eq) = defect (this, i_co, i_eq);

  /* Compute the derivatives: */
  for (j_co = 0; j_co < n_co; j_co++)
    for (j_eq = 0; j_eq < n_equations; j_eq++)
    {
      /* Increment the solution: */
      unknown (this, j_co, j_eq) += h;
      /* Compute the defect once again: */
      if (assemble_defect (this, fvg, time, 1, s_a))
        return __LINE__;
      /* Compute the derivative: */
      for (i_co = 0; i_co < n_co; i_co++)
        for (i_eq = 0; i_eq < n_equations; i_eq++)
          JAC (i_co, i_eq, j_co, j_eq)
            = (defect (this, i_co, i_eq) - DEF (i_co, i_eq)) / h;
      /* Reset the solution: */
      for (i_co = 0; i_co < n_co; i_co++)
        for (i_eq = 0; i_eq < n_equations; i_eq++)
          unknown (this, i_co, i_eq) = SOL (i_co, i_eq);
    };

# undef JAC
# undef DEF
# undef SOL
  return 0;
}

/* num_D_defect - computes the derivative of the dependence
 * of the defect on an import parameter numerically. This is not a virtual
 * function, you should call it from your function. This allows to allocate
 * (or declare) the arrays of the right size. The function returns 0 if
 * OK, and nonzero on an error.
 * THE FUNCTION WORKS WITH NO MORE THAN 5 COMPONENTS of the parameter
 * and NO MORE THAN 10 EQUATIONS IN THE DISCRETIZATION!
 **
 * This function gets the array for the derivatives as an argument (the
 * argument 'deriv'). This should be an array of size MAXNC * <number of
 * equations> * <number of position points> * <number of components of the
 * parameter>. For the structure of the array s. the macro DER in the function.
 */
INT import_param::num_D_defect
(
  FVElementGeometry * fvg,
  DOUBLE time,
  DOUBLE s_a, /* these three are like for compute_D_defect */
  DOUBLE * deriv, /* the array to save the derivative in */
  DOUBLE h /* the step for the numerical integration */
)
{
  DOUBLE prm [MAX_POS_PNTS * 5];
  DOUBLE def [MAXNC * 10];
  INT n_eq = disc->n_equations;
  INT n_co = FVG_NSCV (fvg);
  INT n_points = n_position_points (fvg, p_type);
  INT max_points = max_n_position_points [p_type];
  INT i_comp, i_point, i, j;

  if (n_comp > 5 || n_eq > 10)
    return __LINE__;

# define PRM(point,comp) (prm [(point) * 5 + (comp)])
# define DEF(co,eq) (def [(co) * 10 + (eq)])
# define DER(co,eq,pos,comp) (deriv \
                              [((((co) * n_eq) + (eq)) * max_points + (pos)) * n_comp + (comp)])

  /* Save the original values of the parameter: */
  for (i = 0; i < n_points; i++)
    for (j = 0; j < n_comp; j++)
      PRM (i, j) = (* this)(i, j);

  /* Compute the defect at the original point and save it: */
  if (disc->assemble_defect (disc, fvg, time, 1, s_a))
    return __LINE__;
  for (i = 0; i < n_co; i++)
    for (j = 0; j < n_eq; j++)
      DEF (i, j) = disc->defect (disc, i, j);

  /* Compute the derivatives: */
  for (i_point = 0; i_point < n_points; i_point++)
    for (i_comp = 0; i_comp < n_comp; i_comp++)
    {
      /* Increment the value of the parameter: */
      (* this)(i_point, i_comp) += h;
      /* Compute the defect once again: */
      if (disc->assemble_defect (disc, fvg, time, 1, s_a))
        return __LINE__;
      /* Compute the derivative: */
      for (i = 0; i < n_co; i++)
        for (j = 0; j < n_eq; j++)
          DER (i, j, i_point, i_comp) = (disc->defect (disc, i, j) - DEF (i, j))
                                        / h;
      /* Reset the values of the parameter: */
      for (i = 0; i < n_points; i++)
        for (j = 0; j < n_comp; j++)
          (* this)(i, j) = PRM (i, j);
    };

# undef PRM
# undef DEF
# undef DER
  return 0;
}

/* num_D_param - numerical differentiation of an export parameter w.r.t.
 * the unknowns. This is not to virtual function, too. You should call
 * it from your function. This function returns 0 if OK, nonzero on an error.
 * THE FUNCTION WORKS WITH NO MORE THAN 5 COMPONENTS of the parameter
 * and NO MORE THAN 10 EQUATIONS IN THE DISCRETIZATION!
 **
 * This function gets the array for the derivatives as an argument (the
 * argument 'deriv'). This should be an array of size MAXNC * <number of
 * equations> * <number of position points> * <number of components of the
 * parameter>. For the structure of the array s. the macro DER in the function.
 */
INT export_param::num_D_param
(
  FVElementGeometry * fvg,
  DOUBLE time, /* these two are identical to those of compute_D_param_wrt_unk */
  DOUBLE * deriv, /* to save the derivative */
  DOUBLE h /* the numerical differentiation step */
)
{
  DOUBLE prm [MAX_POS_PNTS * 5];
  DOUBLE sol [MAXNC * 10];
  INT n_eq = disc->n_equations;
  INT n_co = FVG_NSCV (fvg);
  INT n_points = n_position_points (fvg, p_type);
  INT max_points = max_n_position_points [p_type];
  INT i_co, i_eq, i, j;

  if (n_comp > 5 || n_eq > 10)
    return __LINE__;

# define PRM(point,comp) (prm [(point) * 5 + (comp)])
# define SOL(co,eq) (sol [(co) * 10 + (eq)])
# define DER(pos,comp,co,eq) (deriv \
                              [((((co) * n_eq) + (eq)) * max_points + (pos)) * n_comp + (comp)])

  /* Save the original solution: */
  for (i = 0; i < n_co; i++)
    for (j = 0; j < n_eq; j++)
      SOL (i, j) = disc->unknown (disc, i, j);

  /* Compute the parameter and save it: */
  if (compute (fvg, time))
    return __LINE__;
  for (i = 0; i < n_points; i++)
    for (j = 0; j < n_comp; j++)
      PRM (i, j) = (* this)(i, j);

  /* Numerical differentiation: */
  for (i_co = 0; i_co < n_co; i_co++)
    for (i_eq = 0; i_eq < n_eq; i_eq++)
    {
      /* Increment the solution: */
      disc->unknown (disc, i_co, i_eq) += h;
      /* Compute the parameter again: */
      if (compute (fvg, time))
        return __LINE__;
      /* Compute the finite difference: */
      for (i = 0; i < n_points; i++)
        for (j = 0; j < n_comp; j++)
          DER (i, j, i_co, i_eq) = ((* this)(i, j) - PRM (i, j)) / h;
      /* Reset the solution: */
      for (i = 0; i < n_co; i++)
        for (j = 0; j < n_eq; j++)
          disc->unknown (disc, i, j) = SOL (i, j);
    };

# undef PRM
# undef SOL
# undef DER
  return 0;
}

/* End of File */
