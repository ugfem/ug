// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ts.h                                                              */
/*                                                                          */
/* Purpose:   time-stepping scheme for (non-)linear time-dependent problems */
/*                                                                          */
/* Author:    Peter Bastian                                                                                             */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   January  09, 1997  begin                                      */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __TS__
#define __TS__

#include "np.h"
#include "nls.h"
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

#define T_SOLVER_CLASS_NAME "ts"

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct np_t_solver {

  NP_NL_ASSEMBLE nlass;                      /* derived from nonlinear assemble */

  /* things to be initialized by generic init */
  VECDATA_DESC *y;                       /* solution vector                             */
  NP_T_ASSEMBLE *tass;                   /* time assemble numproc			*/
  NP_NL_SOLVER *nlsolve;                 /* nonlinear solver numproc		*/
  VEC_SCALAR reduction;                      /* reduction factor per time step  */
  VEC_SCALAR abslimit;                       /* absolute limit for the defect   */

  /* functions */
  INT (*TimePreProcess)                                  /* called before first time step   */
    (struct np_t_solver *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
  INT (*TimeInit)                        /* initialize, set initial values  */
    (struct np_t_solver *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                                      */
  INT (*TimeStep)                        /* b := b - Ax                     */
    (struct np_t_solver *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                                      */
  INT (*TimePostProcess)                                 /* to be called after last timestep*/
    (struct np_t_solver *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_t_solver NP_T_SOLVER;

typedef INT (*TimePreProcessProcPtr)                                         \
  (NP_T_SOLVER *, INT, INT *);
typedef INT (*TimeInitProcPtr)                                               \
  (NP_T_SOLVER *, INT, INT *);
typedef INT (*TimeStepProcPtr)                                               \
  (NP_T_SOLVER *, INT, INT *);
typedef INT (*TimePostProcessProcPtr)                                        \
  (NP_T_SOLVER *, INT, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* default init function for tsolver num procs */
INT NPTSolverInit (NP_T_SOLVER *np, INT argc, char **argv);

/* default display function for tsolver num procs */
INT NPTSolverDisplay (NP_T_SOLVER *np);

/* default execute function for tsolver num procs */
INT NPTSolverExecute (NP_BASE *np, INT argc, char **argv);

/* init */
INT InitTSolver (void);

#endif
