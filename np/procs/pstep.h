// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pstep.h                                                           */
/*                                                                          */
/* Purpose:   parameter-stepping scheme for nonlin. time-dependent problems */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            IWR/TS                                                        */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheiner Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   August 15, 2000                                               */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __PSTEP__
#define __PSTEP__

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

#define P_STEP_CLASS_NAME "pstep"

#define INIT_PSTEPRESULT(p) p.converged=0;p.number_of_nonlinear_iterations=0;\
  p.number_of_linear_iterations=0;                 \
  p.max_linear_iterations=0;                       \
  p.exec_time=0.0;

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct pstep_result {

  INT converged;                            /* error code                       */
  INT jumped;                                                           /* 1 if nothing had to be done      */
  INT number_of_nonlinear_iterations;       /* number of nonlin. iterations     */
  INT number_of_linear_iterations;          /* number of lin. iterations        */
  INT max_linear_iterations;                /* max lin/nonlin iterations        */
};
typedef struct pstep_result PSTEP_RESULT;

struct np_p_step {

  NP_BASE base;                              /* derived from nonlinear assemble */

  /* things to be initialized by generic init */
  EVECDATA_DESC *sol_p0;                 /* solution vector	0				*/
  EVECDATA_DESC *sol_p1;                 /* solution vector	1				*/

  /* functions */
  INT (*PreProcess)                                              /* called before first time step   */
    (struct np_p_step *,                     /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *sol_p0,               /* solution before step            */
    INT *);                                      /* result                          */
  INT (*Step)                            /* calculate solution at time t1   */
    (struct np_p_step *,                     /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *sol_p0,               /* solution before step            */
    EVECDATA_DESC *sol_p1,               /* solution after step             */
    PSTEP_RESULT *);                             /* result                                      */
  INT (*PostProcess)                                 /* to be called after last timestep*/
    (struct np_p_step *,                     /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_p_step NP_P_STEP;

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

INT InitPStep (void);

#endif
