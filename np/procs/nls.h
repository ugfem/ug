// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  nls.h                                                                                                         */
/*																			*/
/* Purpose:   definition of the assemble num proc type                                  */
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

#ifndef __NLS__
#define __NLS__

#include "np.h"
#include "ls.h"
#include "assemble.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define NL_SOLVER_CLASS_NAME "nl_solver"

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
  INT number_of_nonlinear_iterations;       /* number of iterations             */
  INT number_of_line_searches;              /* number of line search steps      */
  INT rho_first;                            /* first rho                        */
  INT total_linear_iterations;          /* total number                     */
  INT max_linear_iterations;            /* max number of linear iterations  */
  DOUBLE exec_time;                     /* for this nonlinear solve ...     */
} NLRESULT;

struct np_nl_solver {
  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  VECDATA_DESC *x;                       /* solution                        */
  NP_NL_ASSEMBLE *Assemble;              /* the assemble numproc			*/
  VEC_SCALAR reduction;                      /* reduction factor                */
  VEC_SCALAR abslimit;                       /* absolute limit for the defect   */

  /* functions */
  INT (*PreProcess)
    (struct np_nl_solver *,                  /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
  INT (*Solver)                          /* b := b - Ax                     */
    (struct np_nl_solver *,                  /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    NP_NL_ASSEMBLE *,                                /* the assemble numproc			*/
    VEC_SCALAR,                                  /* absolute limit for the defect   */
    VEC_SCALAR,                                  /* reduction factor                */
    NLRESULT *);                                 /* result structure                */
  INT (*PostProcess)
    (struct np_nl_solver *,                  /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
};
typedef struct np_nl_solver NP_NL_SOLVER;

typedef INT (*PreProcessNLSolverProcPtr)                                    \
  (NP_NL_SOLVER *, INT, INT *);
typedef INT (*Solver)                                                       \
  (NP_NL_SOLVER *, INT, VECDATA_DESC *, NP_NL_ASSEMBLE *, VEC_SCALAR *,      \
  VEC_SCALAR *, NLRESULT *);
typedef INT (*PostProcessNLSolverProcPtr)                                   \
  (NP_NL_SOLVER *, INT, VECDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for LinearSolver num procs */
INT NPNLSolverInit (NP_NL_SOLVER *theNP, INT argc , char **argv);

/* generic display function for LinearSolver num procs */
INT NPNLSolverDisplay (NP_NL_SOLVER *theNP);

/* generic execute function for LinearSolver num procs */
INT NPNLSolverExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard LinearSolver num proc type */
INT InitNonlinearSolver (void);

#endif
