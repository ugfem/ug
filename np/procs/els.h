// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  els.h                                                                                                         */
/*																			*/
/* Purpose:   definition of the extended linear solver num proc type            */
/*																			*/
/* Author:    Klaus Johannsen                                               */
/*            IWR/Technische Simulation                                     */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: klaus.johannsen@iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   19.07.02 begin                                                */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __ELS__
#define __ELS__

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

#define ELINEAR_SOLVER_CLASS_NAME "ext_linear_solver"

/* access macros */
#define NPELS_A(p)                              (((NP_ELINEAR_SOLVER*)(p))->A)
#define NPELS_b(p)                              (((NP_ELINEAR_SOLVER*)(p))->b)
#define NPELS_x(p)                              (((NP_ELINEAR_SOLVER*)(p))->x)
#define NPELS_red(p)                    (((NP_ELINEAR_SOLVER*)(p))->reduction)
#define NPELS_abs(p)                    (((NP_ELINEAR_SOLVER*)(p))->abslimit)
#define NPELS_setbasetime(p)    (((NP_ELINEAR_SOLVER*)(p))->setbasetime)
#define NPELS_printbasetime(p)  (((NP_ELINEAR_SOLVER*)(p))->printbasetime)

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/* a data type for returning the status of the computation                  */
typedef struct {
  INT error_code;                           /* error code                       */
  INT converged;                            /* error code                       */
  EVEC_SCALAR first_defect;                 /* first defect                     */
  EVEC_SCALAR last_defect;                  /* last defect                      */
  INT number_of_linear_iterations;          /* number of iterations             */
} ELRESULT;

struct np_elinear_solver {
  NP_BASE base;                              /* inherits base class             */

  /* data (optional, necessary for calling the generic execute routine)   */
  EVECDATA_DESC *x;                      /* solution                        */
  EVECDATA_DESC *b;                      /* defect                          */
  EMATDATA_DESC *A;                      /* matrix                          */
  EVEC_SCALAR reduction;                     /* reduction factor                */
  EVEC_SCALAR abslimit;                      /* absolute limit for the defect   */
  INT setbasetime;                                               /* collect time portions for base level solver */
  INT printbasetime;                                             /* print collected time for base level solver */

  /* functions */
  INT (*PreProcess)
    (struct np_elinear_solver *,             /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector                 */
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    INT *,                                       /* baselevel used by the solver    */
    INT *);                                      /* result                          */
  INT (*Defect)                          /* b := b - Ax                     */
    (struct np_elinear_solver *,             /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector                 */
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    INT *);                                      /* result                          */
  INT (*Residuum)                        /* computes norm of the defect     */
    (struct np_elinear_solver *,             /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    EVECDATA_DESC *,                             /* solution vector                 */
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    ELRESULT *);                                 /* result structure                */
  INT (*Solver)                          /* b := b - Ax                     */
    (struct np_elinear_solver *,             /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector                 */
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    EVEC_SCALAR,                                 /* reduction factor                */
    EVEC_SCALAR,                                 /* absolute limit for the defect   */
    ELRESULT *);                                 /* result structure                */
  INT (*PostProcess)
    (struct np_elinear_solver *,             /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector                 */
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_elinear_solver NP_ELINEAR_SOLVER;

typedef INT (*PreProcessELinearSolverProcPtr)(NP_ELINEAR_SOLVER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, INT *, INT *);
typedef INT (*ELinearDefectProcPtr)(NP_ELINEAR_SOLVER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, INT *);
typedef INT (*ELinearResiduumProcPtr)(NP_ELINEAR_SOLVER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, ELRESULT *);
typedef INT (*ELinearSolverProcPtr)(NP_ELINEAR_SOLVER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, EVEC_SCALAR, EVEC_SCALAR, ELRESULT *);
typedef INT (*PostProcessELinearSolverProcPtr)(NP_ELINEAR_SOLVER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for LinearSolver num procs */
INT NPELinearSolverInit (NP_ELINEAR_SOLVER *theNP, INT argc , char **argv);

/* generic display function for LinearSolver num procs */
INT NPeLinearSolverDisplay (NP_ELINEAR_SOLVER *theNP);

/* generic execute function for LinearSolver num procs */
INT NPELinearSolverExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard LinearSolver num proc type */
INT InitELinearSolver (void);

#endif
