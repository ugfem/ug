// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ls.h                                                                                                          */
/*																			*/
/* Purpose:   definition of the linear solver num proc type                             */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   November 29, 1996                                                                         */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __LS__
#define __LS__

#include "np.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define LINEAR_SOLVER_CLASS_NAME "linear_solver"

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/* a data type for returning the status of the computation                  */
typedef struct {
  INT error_code;                           /* error code                       */
  INT converged;                            /* error code                       */
  VEC_SCALAR first_defect;                  /* first defect                     */
  VEC_SCALAR last_defect;                   /* last defect                      */
  INT number_of_linear_iterations;          /* number of iterations             */
} LRESULT;

struct np_linear_solver {
  NP_BASE base;                              /* inherits base class             */

  /* data (optional, necessary for calling the generic execute routine)   */
  VECDATA_DESC *x;                       /* solution                        */
  VECDATA_DESC *b;                       /* defect                          */
  MATDATA_DESC *A;                       /* matrix                          */
  VEC_SCALAR reduction;                      /* reduction factor                */
  VEC_SCALAR abslimit;                       /* absolute limit for the defect   */

  /* functions */
  INT (*PreProcess)
    (struct np_linear_solver *,              /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *,                                       /* baselevel used by the solver    */
    INT *);                                      /* result                          */
  INT (*Defect)                          /* b := b - Ax                     */
    (struct np_linear_solver *,              /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*Residuum)                        /* computes norm of the defect     */
    (struct np_linear_solver *,              /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    LRESULT *);                                  /* result structure                */
  INT (*Solver)                          /* b := b - Ax                     */
    (struct np_linear_solver *,              /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    VEC_SCALAR,                                  /* reduction factor                */
    VEC_SCALAR,                                  /* absolute limit for the defect   */
    LRESULT *);                                  /* result structure                */
  INT (*PostProcess)
    (struct np_linear_solver *,              /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_linear_solver NP_LINEAR_SOLVER;

typedef INT (*PreProcessLinearSolverProcPtr)                                 \
  (NP_LINEAR_SOLVER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *,   \
  INT *, INT *);
typedef INT (*LinearDefectProcPtr)                                           \
  (NP_LINEAR_SOLVER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *,   \
  INT *);
typedef INT (*LinearResiduumProcPtr)                                         \
  (NP_LINEAR_SOLVER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *,   \
  LRESULT *);
typedef INT (*LinearSolverProcPtr)                                           \
  (NP_LINEAR_SOLVER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *,   \
  VEC_SCALAR, VEC_SCALAR, LRESULT *);
typedef INT (*PostProcessLinearSolverProcPtr)                                \
  (NP_LINEAR_SOLVER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *,   \
  INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for LinearSolver num procs */
INT NPLinearSolverInit (NP_LINEAR_SOLVER *theNP, INT argc , char **argv);

/* generic display function for LinearSolver num procs */
INT NPLinearSolverDisplay (NP_LINEAR_SOLVER *theNP);

/* generic execute function for LinearSolver num procs */
INT NPLinearSolverExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard LinearSolver num proc type */
INT InitLinearSolver (void);

#endif
