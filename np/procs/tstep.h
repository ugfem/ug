// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      tstep.h                                                           */
/*                                                                          */
/* Purpose:   time-stepping scheme for (non-)linear time-dependent problems */
/*                                                                          */
/* Author:    Klaus Johannsen                                                                                           */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   Jan 31, 1998  begin                                           */
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

#ifndef __TSTEP__
#define __TSTEP__

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

#define T_STEP_CLASS_NAME "tstep"

#define INIT_TSTEPRESULT(t) t.converged=0;t.number_of_nonlinear_iterations=0;\
  t.number_of_linear_iterations=0;                 \
  t.max_linear_iterations=0;                       \
  t.exec_time=0.0;

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct tstep_result {

  INT converged;                            /* error code                       */
  INT jumped;                                                           /* 1 if nothing had to be done      */
  INT number_of_nonlinear_iterations;       /* number of nonlin. iterations     */
  INT number_of_linear_iterations;          /* number of lin. iterations        */
  INT max_linear_iterations;                /* max lin/nonlin iterations        */
  DOUBLE exec_time;                     /* for this nonlinear solve ...     */
};
typedef struct tstep_result TSTEP_RESULT;

struct np_t_step {

  NP_BASE base;                              /* derived from nonlinear assemble */

  /* things to be initialized by generic init */
  DOUBLE t0;                                 /* time 0                          */
  VECDATA_DESC *sol_t0;                  /* solution vector at time t0		*/
  DOUBLE t1;                                 /* time 1                          */
  VECDATA_DESC *sol_t1;                  /* solution vector at time t1		*/

  /* functions */
  INT (*TimePreProcess)                                  /* called before first time step   */
    (struct np_t_step *,                     /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
  INT (*TimeInit)                        /* set initial values              */
    (struct np_t_step *,                     /* pointer to (derived) object     */
    INT,                                         /* level                           */
    DOUBLE,                                      /* t0                              */
    VECDATA_DESC *sol_t0,                /* initial values                  */
    INT *);                                      /* result                                      */
  INT (*TimeStep)                        /* calculate solution at time t1   */
    (struct np_t_step *,                     /* pointer to (derived) object     */
    INT,                                         /* level                           */
    DOUBLE,                                      /* t0                              */
    VECDATA_DESC *sol_t0,                /* solution at t0                  */
    DOUBLE,                                      /* t1                              */
    VECDATA_DESC *sol_t1,                /* solution at t1                  */
    TSTEP_RESULT *);                             /* result                                      */
  INT (*TimePostProcess)                                 /* to be called after last timestep*/
    (struct np_t_step *,                     /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_t_step NP_T_STEP;

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

INT InitTStep (void);

#endif
