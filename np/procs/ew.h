// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ew.h                                                                                                          */
/*																			*/
/* Purpose:   definition of the eigenvalue solver num proc                              */
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

#ifndef __EW__
#define __EW__

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

#define MAX_NUMBER_EW  20
#define EW_SOLVER_CLASS_NAME "ew"

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/* a data type for returning the status of the computation                  */
typedef struct {
  INT error_code;                           /* error code                       */
  LRESULT lresult[MAX_NUMBER_EW];       /* result of linear solver          */
  VEC_SCALAR defect[MAX_NUMBER_EW];     /* defect                           */
  INT converged[MAX_NUMBER_EW];             /* convergence flag                 */
  INT number_of_iterations[MAX_NUMBER_EW];        /* number of iterations       */
} EWRESULT;

struct np_ew_solver {

  NP_BASE base;                              /* inherits base class             */

  /* data (optional, necessary for calling the generic execute routine)   */
  INT nev;                               /* number of eigenvectors          */
  VECDATA_DESC *ev[MAX_NUMBER_EW];       /* eigenvectors                    */
  DOUBLE ew[MAX_NUMBER_EW];              /* eigenvalues                     */
  NP_NL_ASSEMBLE *Assemble;              /* assembling stiffness matrix     */
                                         /* and right hand side             */
  VEC_SCALAR reduction;                      /* reduction factor                */
  VEC_SCALAR abslimit;                       /* absolute limit for the defect   */

  /* functions */
  INT (*PreProcess)
    (struct np_ew_solver *,                  /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT,                                         /* number of eigenvectors          */
    VECDATA_DESC **,                             /* eigenvectors                    */
    NP_NL_ASSEMBLE *,                            /* matrix and right hand side      */
    INT *);                                      /* result                          */
  INT (*Rayleigh)
    (struct np_ew_solver *,                  /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* eigenvector                     */
    NP_NL_ASSEMBLE *,                            /* matrix and right hand side      */
    DOUBLE *,                                    /* numerator and denominator       */
    DOUBLE *,                                    /* quotient                        */
    INT *);                                      /* result                          */
  INT (*Solver)
    (struct np_ew_solver *,                  /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT,                                         /* number of eigenvectors          */
    VECDATA_DESC **,                             /* eigenvectors                    */
    DOUBLE *,                                    /* eigenvalues                     */
    NP_NL_ASSEMBLE *,                            /* matrix and right hand side      */
    VEC_SCALAR,                                  /* reduction factor                */
    VEC_SCALAR,                                  /* absolut limit for the defect    */
    EWRESULT *);                                 /* result structure                */
  INT (*PostProcess)
    (struct np_ew_solver *,                  /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT,                                         /* number of eigenvectors          */
    VECDATA_DESC **,                             /* eigenvectors                    */
    NP_NL_ASSEMBLE *,                            /* matrix and right hand side      */
    INT *);                                      /* result                          */
};
typedef struct np_ew_solver NP_EW_SOLVER;

typedef INT (*PreProcessEWSolverProcPtr)                                     \
  (NP_EW_SOLVER *, INT, INT, VECDATA_DESC **, NP_NL_ASSEMBLE *, INT *);
typedef INT (*PreProcessRayleighProcPtr)                                     \
  (NP_EW_SOLVER *, INT, VECDATA_DESC *,                                       \
  NP_NL_ASSEMBLE *, DOUBLE *, DOUBLE *, INT *);
typedef INT (*EWSolverProcPtr)                                               \
  (NP_EW_SOLVER *, INT, INT, VECDATA_DESC **, DOUBLE *,                       \
  NP_NL_ASSEMBLE *, VEC_SCALAR, VEC_SCALAR, EWRESULT *);
typedef INT (*PostProcessEWSolverProcPtr)                                    \
  (NP_EW_SOLVER *, INT, INT, VECDATA_DESC **, NP_NL_ASSEMBLE *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for LinearSolver num procs */
INT NPEWSolverInit (NP_EW_SOLVER *theNP, INT argc , char **argv);

/* generic display function for LinearSolver num procs */
INT NPEWSolverDisplay (NP_EW_SOLVER *theNP);

/* generic execute function for LinearSolver num procs */
INT NPEWSolverExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard LinearSolver num proc type */
INT InitEW (void);

#endif
