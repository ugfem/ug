// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  bdf.h                                                                                                         */
/*																			*/
/* Purpose:   header file for BDF time solver								*/
/*																			*/
/* Author:	  Peter Bastian                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   JAN 10 1997                                                                                       */
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

#ifndef __BDF__
#define __BDF__

#include "db.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

typedef struct
{
  NP_T_SOLVER tsolver;                   /* derived from class NP_T_SOLVER  */

  /* local variables */
  INT step;                              /* number of time step             */
  DOUBLE dt;                             /* size of time step               */
  DOUBLE t_p1;                           /* time t_k+1                      */
  DOUBLE t_0;                            /* time t_k                        */
  DOUBLE t_m1;                           /* time t_k-1                      */
  NP_ORDERED_LIST *TimeControl;          /* list for time steps             */
  DOUBLE disabled_timestep;              /* storage for timestep disabled   */
                                         /* TimeControl                     */

  /* parameters (to be set with init function */
  INT baselevel;                         /* for nested iteration            */
  INT order;                             /* 1,2 are allowed                 */
  INT predictorder;                      /* 0,1 are allowed                 */
  INT nested;                            /* use nested iteration            */
  INT nlinterpolate;                     /* nonlinear interpolation         */
  INT presteps;                                               /* number of steps for start grid */
  INT optnlsteps;                        /* optimal number of nonlin. steps */
  INT rep;                               /* for repeat solver after grid chg*/
  INT Break;                             /* break after error estimator     */
  INT Continue;                          /* continue after error estimator  */
  INT refarg;                            /* refine copy all                 */
  INT noabort;                           /* take last iterate if not converged */
  DOUBLE tstart;                         /* start time                      */
  DOUBLE dtstart;                        /* time step to begin with         */
  DOUBLE dtmin;                          /* smallest time step allowed      */
  DOUBLE dtmax;                          /* largest time step allowed       */
  DOUBLE dtscale;                        /* scaling factor applied after ste*/
  DOUBLE rhogood;                        /* threshold for step doubling     */
  NP_TRANSFER *trans;                    /* uses transgrid for nested iter  */
  NP_ERROR *error;                       /* error indicator                 */
  NP_ERROR *initerror;                       /* error indicator for error in initial conditions                */
  INT ctn;                               /* change to nested iteration      */
  INT hist;
  INT list_i;
  INT list_n;
  DOUBLE list_dt[50];
  DOUBLE list_work[50];
  INT displayMode;
  char scaleName[NAMELEN];
  DOUBLE scale;

  /* statistics */
  INT number_of_nonlinear_iterations;   /* number of iterations             */
  INT total_linear_iterations;          /* total number                     */
  INT max_linear_iterations;            /* max number of linear iterations  */
  DOUBLE exec_time;                     /* for nonlinear solver ...         */

  /* and XDATA_DESCs */
  VECDATA_DESC *y_p1;                    /* solution y_k+1                  */
  VECDATA_DESC *y_0;                     /* solution y_k                    */
  VECDATA_DESC *y_m1;                    /* solution y_k-1                  */
  VECDATA_DESC *b;                       /* saved nonlinear solution        */

} NP_BDF;                                /*final class implementing BDF(1,2)*/

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* create standard LinearSolver num proc type */
INT BDFPreProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, INT *res);
INT BDFPostProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, VECDATA_DESC *d, MATDATA_DESC *J, INT *res);
INT BDFAssembleSolution (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, INT *res);
INT BDFAssembleDefect (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, VECDATA_DESC *d, MATDATA_DESC *J, INT *res);
INT BDFAssembleMatrix (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *J, INT *res);
INT BDFNAssembleMatrix (NP_NL_ASSEMBLE *ass, INT fl, INT tl, NODE *n, VECDATA_DESC *u, VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *J, INT *res);
INT BDFTimePreProcess (NP_T_SOLVER *ts, INT level, INT *res);
INT BDFTimeInit (NP_T_SOLVER *ts, INT level, INT *res);
INT BDFTimePostProcess (NP_T_SOLVER *ts, INT level, INT *res);
INT InitBDFSolver (void);
INT BDFInit (NP_BASE *base, INT argc, char **argv);
INT BDFDisplay (NP_BASE *theNumProc);

#endif
