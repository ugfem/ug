// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_solve.h													*/
/*																			*/
/* Purpose:   algebraic multigrid solver									*/
/*																			*/
/* Author:	  Peter Bastian                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   05 FEB 1996 Begin												*/
/*			  02 OKT 1997 redesign											*/
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

#ifndef __AMG_SOLVE__
#define __AMG_SOLVE__

#include "amg_coarsen.h"

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
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

#define AMG_DJAC                        1               /* smoother types						*/
#define AMG_SOR                         2
#define AMG_SSOR                        3
#define AMG_ILU                         4
#define AMG_MGC                         5
#define AMG_EX                          6

#define AMG_LS                          5               /* solver types							*/
#define AMG_CG                          6
#define AMG_BCGS                        7

typedef struct {                                        /* parameters for solver				*/
  int verbose;                                          /* be verbose							*/

  /* fine grid solver */
  int solver;                                                   /* type of solver to be used			*/
  int preconditioner;                                   /* type of preconditioner				*/
  int maxit;                                                    /* max number of iterations				*/
  int ex_maxit;                     /* 1 to execute exactly maxit iterations*/
  double red_factor;                                    /* reqired reduction in residual norm	*/
  double dnorm_min;                                     /* convergence limit					*/

  /* coarse grid solver */
  int coarse_smoother;                          /* type of smoother on coarse grid		*/
  int coarse_maxit;                                     /* iteration number for coarse grid sol */
  double coarse_red_factor;                     /* same for coarsest level				*/

  /* multigrid cycle */
  int n1,n2;                                                    /* pre and post smoothing				*/
  int gamma;                                                    /* cycle form							*/
  double omega_p[AMG_MAX_COMP];         /* damping factor for coarse grid corr	*/

  /* smoother */
  int smoother;                                         /* type of smoother to be used			*/
  double omega[AMG_MAX_COMP];                   /* damping factor per component			*/

} AMG_SolverContext;


/****************************************************************************/
/*																			*/
/* functions																*/
/*																			*/
/****************************************************************************/

int AMG_Build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A);
int AMG_Solve (AMG_VECTOR *x, AMG_VECTOR *b);

#endif
