// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      bdf.c                                                             */
/*                                                                          */
/* Purpose:   implement BDF(1) and BDF(2) as a tsolver                                          */
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
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "devices.h"
#include "general.h"
#include "ugstruct.h"
#include "gm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"

#include "nls.h"
#include "ls.h"
#include "nls.h"
#include "assemble.h"
#include "transfer.h"
#include "ts.h"
#include "bdf.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/*  Class Definition                                                                    */
/*                                                                          */
/****************************************************************************/

typedef struct
{
  NP_T_SOLVER tsolver;                                   /* derived from class NP_T_SOLVER	*/

  /* local variables */
  INT step;                                                              /* number of time step				*/
  DOUBLE dt;                                                     /* size of time step				*/
  DOUBLE t_p1;                                                   /* time t_k+1                                  */
  DOUBLE t_0;                                                        /* time t_k                                        */
  DOUBLE t_m1;                                                   /* time t_k-1                                  */

  /* parameters (to be set with init function */
  INT baselevel;                                                 /* for nested iteration		    */
  INT order;                                                             /* 1,2 are allowed					*/
  INT predictorder;                                              /* 0,1 are allowed					*/
  INT nested;                                                            /* use nested iteration                        */
  DOUBLE dtstart;                                                /* time step to begin with			*/
  DOUBLE dtmin;                                                  /* smallest time step allowed		*/
  DOUBLE dtmax;                                                  /* largest time step allowed		*/
  DOUBLE dtscale;                                                /* scaling factor applied after ste*/
  DOUBLE rhogood;                                                /* threshold for step doubling		*/
  NP_TRANSFER *trans;                                            /* uses transgrid for nested iter  */

  /* statistics */
  INT number_of_nonlinear_iterations;       /* number of iterations             */
  INT total_linear_iterations;          /* total number                     */
  INT max_linear_iterations;            /* max number of linear iterations  */
  DOUBLE exec_time;                     /* for nonlinear solver ...             */

  /* and XDATA_DESCs */
  VECDATA_DESC *y_p1;                    /* solution y_k+1                                      */
  VECDATA_DESC *y_0;                     /* solution y_k                                        */
  VECDATA_DESC *y_m1;                    /* solution y_k-1                                      */
  VECDATA_DESC *b;                                               /* saved nonlinear solution		*/

} NP_BDF;                                                                /*final class implementing BDF(1,2)*/


/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* Nonlinear Assemble Interface provided to nonlinear solver				*/
/* REMEMBER: NP_T_SOLVER is derived from NP_NL_ASSEMBLE.                                */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/

static INT BDFPreProcess
  (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, INT *res)
{
  return(0);
}

static INT BDFAssembleSolution
  (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, INT *res)
{
  NP_BDF *bdf;
  NP_T_ASSEMBLE *tass;

  /* get numprocs ... */
  bdf = (NP_BDF *) ass;                 /* this is the trick, tsolver is derived	*/
  /* from nonlinear assemble !				*/
  tass = bdf->tsolver.tass;             /* since we need to access it quite often       */

  /* now call time assemble with correct time value */
  return((*tass->TAssembleSolution)(tass,fl,tl,bdf->t_p1,u,res));
}

static INT BDFAssembleDefect
  (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u,
  VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  NP_BDF *bdf;
  NP_T_ASSEMBLE *tass;
  DOUBLE s_a,s_m;
  DOUBLE dt_p1,dt_0,g_p1,g_0,g_m1;

  /* get numprocs ... */
  bdf = (NP_BDF *) ass;                 /* this is the trick, tsolver is derived	*/
  tass = bdf->tsolver.tass;             /* since we need to access it quite often       */

  /* compute coefficients */
  dt_p1 = bdf->t_p1-bdf->t_0;
  dt_0  = bdf->t_0-bdf->t_m1;
  g_p1  = (dt_0+2*dt_p1)/(dt_0+dt_p1);
  g_0   = -(dt_0+dt_p1)/dt_0;
  g_m1  = dt_p1*dt_p1/(dt_0*dt_0+dt_0*dt_p1);

  /* compute scaling factors depending on selected order */
  switch (bdf->order)
  {
  case 1 : s_m = 1.0; s_a = -dt_p1; break;
  case 2 : s_m = 1.0; s_a = -dt_p1/g_p1; break;
  default :
    UserWrite("BDFAssembleDefect: invalid order\n");
    return(1);
  }

  /* copy precomputed part of defect */
  a_dcopy(ass->base.mg,fl,tl,d,EVERY_CLASS,bdf->b);

  /* call function from time assemble interface */
  return( (*tass->TAssembleDefect)(tass,fl,tl,bdf->t_p1,s_m,s_a,u,d,J,res) );
}

static INT BDFAssembleMatrix
  (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u,
  VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *J, INT *res)
{
  NP_BDF *bdf;
  NP_T_ASSEMBLE *tass;
  DOUBLE s_a;
  DOUBLE dt_p1,dt_0,g_p1,g_0,g_m1;

  /* get numprocs ... */
  bdf = (NP_BDF *) ass;                 /* this is the trick, tsolver is derived	*/
  tass = bdf->tsolver.tass;             /* since we need to access it quite often       */

  /* compute coefficients */
  dt_p1 = bdf->t_p1-bdf->t_0;
  dt_0  = bdf->t_0-bdf->t_m1;
  g_p1  = (dt_0+2*dt_p1)/(dt_0+dt_p1);
  g_0   = -(dt_0+dt_p1)/dt_0;
  g_m1  = dt_p1*dt_p1/(dt_0*dt_0+dt_0*dt_p1);

  /* compute scaling factors depending on selected order */
  switch (bdf->order)
  {
  case 1 : s_a = -dt_p1; break;
  case 2 : s_a = -dt_p1/g_p1; break;
  default :
    UserWrite("BDFAssembleMatrix: invalid order\n");
    return(1);
  }

  /* call function from time assemble interface */
  return( (*tass->TAssembleMatrix)(tass,fl,tl,bdf->t_p1,s_a,u,d,v,J,res) );
}

static INT BDFPostProcess
  (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x,
  VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  return(0);
}

/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* Implementation of NP_T_SOLVER functions                                                              */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/

static INT TimePreProcess (NP_T_SOLVER *ts, INT level, INT *res)
{
  NP_BDF *bdf;

  /* get numprocs ... */
  bdf = (NP_BDF *) ts;

  /* allocate XDATA_DESCs here */
  if (ts->y == NULL)
  {
    UserWrite("solution y is not defined\n");
    return(__LINE__);
  }
  if (AllocVDFromVD(ts->nlass.base.mg,0,level,ts->y,&(bdf->y_p1))) return(__LINE__);
  if (AllocVDFromVD(ts->nlass.base.mg,0,level,ts->y,&(bdf->y_m1))) return(__LINE__);
  if (AllocVDFromVD(ts->nlass.base.mg,0,level,ts->y,&(bdf->b))) return(__LINE__);

  return(0);
}

static INT TimeInit (NP_T_SOLVER *ts, INT level, INT *res)
{
  NP_BDF *bdf;
  NP_T_ASSEMBLE *tass;
  char buffer[128];

  /* get numprocs ... */
  bdf = (NP_BDF *) ts;
  tass = bdf->tsolver.tass;

  /* initialize bdf local variables */
  bdf->dt = bdf->dtstart;
  bdf->step = 0;
  bdf->t_0 = 0.0;       /* we always start with t = 0 ! */
  bdf->t_m1 = - bdf->dt;

  /* set initial values and boundary conditions in y_0 */
  *res = 1;
  if ( (*tass->TAssembleInitial)(tass,0,level,bdf->y_0,res) )
    return(1);
  if ( (*tass->TAssembleSolution)(tass,0,level,bdf->t_0,bdf->y_0,res) )
    return(1);

  /* write time to shell */
  sprintf(buffer,"%12.4lE",bdf->t_0);
  SetStringVar("TIME",buffer);
  sprintf(buffer,"%12.4lE",bdf->dt);
  SetStringVar("TIMESTEP",buffer);

  /* statistics init */
  bdf->number_of_nonlinear_iterations = 0;
  bdf->total_linear_iterations = 0;
  bdf->max_linear_iterations = 0;
  bdf->exec_time = 0.0;

  /* return ok */
  *res = 0;
  return(*res);
}

static INT TimeStep (NP_T_SOLVER *ts, INT level, INT *res)
{
  NP_BDF *bdf;
  NP_T_ASSEMBLE *tass;
  NP_NL_SOLVER *nlsolve;
  DOUBLE dt_p1,dt_0,g_p1,g_0,g_m1;
  DOUBLE Factor[MAX_VEC_COMP];
  INT n_unk;
  INT i,k;
  INT low,bl;
  INT verygood,bad;
  NLRESULT nlresult;
  MULTIGRID *mg;
  char buffer[128];

  /* get numprocs ... */
  bdf = (NP_BDF *) ts;                  /* this is the trick, tsolver is derived	*/
  tass = bdf->tsolver.tass;             /* since we need to access it quite often       */
  nlsolve = ts->nlsolve;
  mg = nlsolve->base.mg;

  /* get number of components */
  n_unk = VD_NCOMP(bdf->y_p1);

  /* initialize strategy flags */
  verygood = bad = 0;

  /* compute solution at new time step ON CURRENT GRID */
  while (1)
  {
    /* advance time level with current time step */
    bdf->t_p1 = bdf->t_0+bdf->dt;

    /* compute coefficients */
    dt_p1 = bdf->t_p1-bdf->t_0;
    dt_0  = bdf->t_0-bdf->t_m1;
    g_p1  = (dt_0+2*dt_p1)/(dt_0+dt_p1);
    g_0   = -(dt_0+dt_p1)/dt_0;
    g_m1  = dt_p1*dt_p1/(dt_0*dt_0+dt_0*dt_p1);

    /* determine level where to predict to new time step */
    if (bdf->nested) low = MIN(bdf->baselevel,level);
    else low = level;

    /* predict to new time step on level low */
    a_dcopy(mg,0,low,bdf->y_p1,EVERY_CLASS,bdf->y_0);
    if (bdf->predictorder==1 && bdf->step>0)
    {
      for (i=0; i<n_unk; i++) Factor[i] = (bdf->t_p1-bdf->t_m1)/dt_0;
      a_dscale(mg,0,low,bdf->y_p1,EVERY_CLASS,Factor);
      for (i=0; i<n_unk; i++) Factor[i] = 1.0-(bdf->t_p1-bdf->t_m1)/dt_0;
      a_daxpy(mg,0,low,bdf->y_p1,EVERY_CLASS,Factor,bdf->y_m1);
    }

    if ( (*tass->TAssemblePreProcess)(tass,0,level,
                                      bdf->t_p1,bdf->t_0,bdf->t_m1,
                                      bdf->y_p1,bdf->y_0,bdf->y_m1,res) )
      return(__LINE__);

    /* set Dirichlet conditions in predicted solution */
    if ( (*tass->TAssembleSolution)(tass,0,low,bdf->t_p1,bdf->y_p1,res) )
      return(__LINE__);

    /* do nested iteration on new time step */
    for (k=low; k<=level; k++)
    {
      if (bdf->nested) UserWriteF("Nested Iteration on level %d\n",k);

      /* prepare constant part of defect */
      a_dset(mg,0,k,bdf->b,EVERY_CLASS,0.0);
      if (bdf->order==1)
      {
        if ( (*tass->TAssembleDefect)(tass,0,k,bdf->t_0,-1.0,0.0,bdf->y_0,bdf->b,NULL,res) )
          return(__LINE__);
      }
      else
      {
        if ( (*tass->TAssembleDefect)(tass,0,k,bdf->t_0,g_0/g_p1,0.0,bdf->y_0,bdf->b,NULL,res) )
          return(__LINE__);
        if ( (*tass->TAssembleDefect)(tass,0,k,bdf->t_m1,g_m1/g_p1,0.0,bdf->y_m1,bdf->b,NULL,res) )
          return(__LINE__);
      }

      /* solve nonlinear problem on level k */
      if (nlsolve->PreProcess!=NULL)
        if ( (*nlsolve->PreProcess)(nlsolve,k,res) )
          return(__LINE__);
      if ( (*nlsolve->Solver)(nlsolve,k,bdf->y_p1,&bdf->tsolver.nlass,nlsolve->abslimit,nlsolve->reduction,&nlresult) )
        return(__LINE__);

      /* update statisitics */
      bdf->number_of_nonlinear_iterations += nlresult.number_of_nonlinear_iterations;
      bdf->total_linear_iterations += nlresult.total_linear_iterations;
      bdf->max_linear_iterations = MAX(bdf->max_linear_iterations,nlresult.max_linear_iterations);
      bdf->exec_time += nlresult.exec_time;

      if (nlsolve->PostProcess!=NULL)
        if ( (*nlsolve->PostProcess)(nlsolve,k,bdf->y_p1,res) )
          return(__LINE__);
      if (!nlresult.converged)
      {
        bdf->dt *= 0.5;                         /* reduce time step     */
        if (bdf->dt<bdf->dtmin) {
          UserWrite("time step too small -- aborting\n");
          return(__LINE__);
        }
        UserWrite("reduced time step\n");
        bad=1;
        break;                               /* and try all over again  */
      }
      else {
        if (nlresult.rho_first<=bdf->rhogood)
          verygood=1;
        else
          verygood=0;
      }

      /* interpolate up */
      if (k<level)
      {
        if (bdf->trans->PreProcessProject!=NULL)
          if ((*bdf->trans->PreProcessProject)(bdf->trans,k+1,&bl,res))
            return(__LINE__);
        for (i=0; i<n_unk; i++) Factor[i] = 1.0;
        if ((*bdf->trans->InterpolateCorrection)(bdf->trans,k+1,bdf->y_p1,bdf->y_p1,NULL,Factor,res))
          return(__LINE__);
        if (bdf->trans->PostProcessProject!=NULL)
          if ((*bdf->trans->PostProcessProject)(bdf->trans,k+1,res))
            return(__LINE__);

        /* set Dirichlet conditions in predicted solution */
        if ( (*tass->TAssembleSolution)(tass,k+1,k+1,bdf->t_p1,bdf->y_p1,res) )
          return(__LINE__);
      }
    }

    /* check convergence */
    if (!nlresult.converged) continue;                  /* start again with smaller time step   */
    else break;                                                                 /* exit while loop						*/
  }

  /* LATER: In the adaptive case now iterate: adapt grid & time step, solve again */

  if ( (*tass->TAssemblePostProcess)(tass,0,level,
                                     bdf->t_p1,bdf->t_0,bdf->t_m1,
                                     bdf->y_p1,bdf->y_0,bdf->y_m1,res) )
    return(__LINE__);

  /* accept new time step */
  a_dcopy(mg,0,level,bdf->y_m1,EVERY_CLASS,bdf->y_0 ); bdf->t_m1 = bdf->t_0;
  a_dcopy(mg,0,level,bdf->y_0 ,EVERY_CLASS,bdf->y_p1); bdf->t_0  = bdf->t_p1;
  bdf->step++;

  /* write time to shell */
  sprintf(buffer,"%12.4lE",bdf->t_0);
  SetStringVar("TIME",buffer);
  sprintf(buffer,"%12.4lE",bdf->dt);
  SetStringVar("TIMESTEP",buffer);

  UserWriteF("TIMESTEP %4d: TIME=%10.4lg DT=%10.4lg EXECT=%10.4lg NLIT=%5d LIT=%5d MAXLIT=%3d\n",
             bdf->step,bdf->t_0,bdf->dt,bdf->exec_time,bdf->number_of_nonlinear_iterations,
             bdf->total_linear_iterations,bdf->max_linear_iterations);

  /* chose new dt for next time step */
  if (verygood && (!bad) && bdf->dt*2<=bdf->dtmax)
  {
    bdf->dt *= 2.0;
    UserWrite("doubling time step\n");
    *res=0; return(*res);
  }
  if ((!bad) && bdf->dt*bdf->dtscale<=bdf->dtmax && bdf->dt*bdf->dtscale>=bdf->dtmin)
  {
    bdf->dt *= bdf->dtscale;
    *res=0; return(*res);
  }


  return(0);
}

static INT TimePostProcess (NP_T_SOLVER *ts, INT level, INT *res)
{
  NP_BDF *bdf;

  /* get numprocs ... */
  bdf = (NP_BDF *) ts;

  /* free XDATA_DESCs here */
  FreeVD(ts->nlass.base.mg,0,level,bdf->y_0);
  FreeVD(ts->nlass.base.mg,0,level,bdf->y_m1);
  FreeVD(ts->nlass.base.mg,0,level,bdf->b);

  return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  Init															*/
/*																			*/
/* Purpose:   init solve cycle												*/
/*																			*/
/* Input:         NumProcType Init Function										*/
/*																			*/
/* Output: INT 0: ok														*/
/*			   1: error														*/
/*																			*/
/****************************************************************************/

static INT BDFInit (NP_BASE *base, INT argc, char **argv)
{
  NP_BDF *bdf;
  INT r;

  /* get numprocs ... */
  bdf = (NP_BDF *) base;

  /* call tsolver init */
  r = NPTSolverInit(&bdf->tsolver,argc,argv);

  /* read data descs */
  bdf->y_0 = bdf->tsolver.y;       /* allocated already in tsolver */
  bdf->y_p1 = ReadArgvVecDesc(base->mg,"yp1",argc,argv);
  bdf->y_m1 = ReadArgvVecDesc(base->mg,"ym1",argc,argv);
  bdf->b = ReadArgvVecDesc(base->mg,"b",argc,argv);

  /* read other numprocs */
  bdf->trans = (NP_TRANSFER *) ReadArgvNumProc(base->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  if (bdf->trans == NULL) return(NP_NOT_ACTIVE);

  /* set configuration parameters */
  if (ReadArgvINT("baselevel",&(bdf->baselevel),argc,argv))
  {
    UserWrite("default: baselevel=0\n");
    bdf->baselevel=0;
  }
  if ((bdf->baselevel<0)||(bdf->baselevel>32)) return(NP_NOT_ACTIVE);

  if (ReadArgvINT("order",&(bdf->order),argc,argv))
  {
    UserWrite("default: order=1\n");
    bdf->order=1;
  }
  if ((bdf->order<1)||(bdf->order>2)) return(NP_NOT_ACTIVE);

  if (ReadArgvINT("predictorder",&(bdf->predictorder),argc,argv))
  {
    UserWrite("default: predictorder=0\n");
    bdf->predictorder=0;
  }
  if ((bdf->predictorder<0)||(bdf->predictorder>1)) return(NP_NOT_ACTIVE);

  if (ReadArgvINT("nested",&(bdf->nested),argc,argv))
  {
    UserWrite("default: nested=0\n");
    bdf->nested=0;
  }
  if ((bdf->nested<0)||(bdf->nested>1)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("dtstart",&(bdf->dtstart),argc,argv))
  {
    UserWrite("dtstart must be specified\n");
    return(NP_NOT_ACTIVE);
  }
  if ((bdf->dtstart<0.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("dtmin",&(bdf->dtmin),argc,argv))
  {
    UserWrite("dtmin must be specified\n");
    return(NP_NOT_ACTIVE);
  }
  if ((bdf->dtmin<0.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("dtmax",&(bdf->dtmax),argc,argv))
  {
    UserWrite("dtmax must be specified\n");
    return(NP_NOT_ACTIVE);
  }
  if ((bdf->dtmax<0.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("dtscale",&(bdf->dtscale),argc,argv))
  {
    UserWrite("dtscale must be specified\n");
    return(NP_NOT_ACTIVE);
  }
  if ((bdf->dtscale<0.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("rhogood",&(bdf->rhogood),argc,argv))
  {
    UserWrite("default: rhogood=0.01\n");
    bdf->rhogood = 0.01;
  }
  if ((bdf->rhogood<0.0) || (bdf->rhogood>1)) return(NP_NOT_ACTIVE);

  return (r);
}

/****************************************************************************/
/*																			*/
/* Function:  Display														*/
/*																			*/
/* Purpose:   display linear multigrid cycle								*/
/*																			*/
/* Input:     none															*/
/*																			*/
/* Output:    INT 0: ok														*/
/*			   1: error														*/
/*																			*/
/****************************************************************************/

static INT BDFDisplay (NP_BASE *theNumProc)
{
  NP_BDF *bdf;

  /* get numprocs ... */
  bdf = (NP_BDF *) theNumProc;

  /* general ts display */
  NPTSolverDisplay(&bdf->tsolver);

  UserWrite("BDF data:\n");
  if (bdf->y_p1 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"y_p1",ENVITEM_NAME(bdf->y_p1));
  if (bdf->y_0  != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"y_0 ",ENVITEM_NAME(bdf->y_0 ));
  if (bdf->y_m1 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"y_m1",ENVITEM_NAME(bdf->y_m1));
  if (bdf->b    != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"b   ",ENVITEM_NAME(bdf->b   ));

  return (0);
}


/****************************************************************************/
/*																			*/
/* Function:  Execute														*/
/*																			*/
/* Purpose:   init solve cycle												*/
/*																			*/
/* Input:         NumProcType Init Function										*/
/*																			*/
/* Output: INT 0: ok														*/
/*			   1: error														*/
/*																			*/
/****************************************************************************/

static INT BDFExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_BDF *bdf;
  NP_T_SOLVER *np;
  INT result,level;

  /* get numprocs ... */
  np = (NP_T_SOLVER *) theNP;
  bdf = (NP_BDF *) theNP;

  /* get level */
  level = CURRENTLEVEL(theNP->mg);

  /* check functions */
  if (np->y == NULL) {
    PrintErrorMessage('E',"BDFExecute","no vector y");
    return (1);
  }
  if (np->tass == NULL) {
    PrintErrorMessage('E',"BDFExecute","no assemble num proc");
    return (1);
  }
  if (np->nlsolve == NULL) {
    PrintErrorMessage('E',"BDFExecute","no solver num proc");
    return (1);
  }

  /* call preprocess function, allocates vectors ! */
  if (ReadArgvOption("pre",argc,argv)) {
    if (np->TimePreProcess != NULL)
      if ((*np->TimePreProcess)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimePreProcess failed, error code %d\n",
                   result);
        return (1);
      }
  }

  /* set initial values */
  if (ReadArgvOption("init",argc,argv)) {
    if (np->TimeInit != NULL)
      if ((*np->TimeInit)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeInit failed, error code %d\n",
                   result);
        return (1);
      }
  }

  /* execute bdf1, nonnested */
  if (ReadArgvOption("bdf1",argc,argv)) {
    bdf->order = 1;
    bdf->nested = 0;
    if (np->TimeStep != NULL)
      if ((*np->TimeStep)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeStep failed, error code %d\n",
                   result);
        return (1);
      }
  }

  /* execute bdf2, nonnested */
  if (ReadArgvOption("bdf2",argc,argv)) {
    bdf->order = 2;
    bdf->nested = 0;
    if (np->TimeStep != NULL)
      if ((*np->TimeStep)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeStep failed, error code %d\n",
                   result);
        return (1);
      }
  }

  /* execute bdf1, nested */
  if (ReadArgvOption("bdf1n",argc,argv)) {
    bdf->order = 1;
    bdf->nested = 1;
    if (np->TimeStep != NULL)
      if ((*np->TimeStep)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeStep failed, error code %d\n",
                   result);
        return (1);
      }
  }

  /* execute bdf2, nested */
  if (ReadArgvOption("bdf2n",argc,argv)) {
    bdf->order = 2;
    bdf->nested = 1;
    if (np->TimeStep != NULL)
      if ((*np->TimeStep)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeStep failed, error code %d\n",
                   result);
        return (1);
      }
  }

  /* postprocess, free vectors */
  if (ReadArgvOption("post",argc,argv)) {
    if (np->TimePostProcess != NULL)
      if ((*np->TimePostProcess)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimePostProcess failed, error code %d\n",
                   result);
        return (1);
      }
  }

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  BDFConstruct													*/
/*																			*/
/* Purpose:	  create an object of class bdf									*/
/*																			*/
/* Input:	  theNP	num proc data structure									*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*			  else : error													*/
/*																			*/
/****************************************************************************/


static INT BDFConstruct (NP_BASE *theNP)
{
  NP_T_SOLVER *ts;
  NP_NL_ASSEMBLE *na;

  /* get numprocs ... */
  ts = (NP_T_SOLVER *) theNP;
  na = &ts->nlass;

  /* set base functions */
  theNP->Init    = BDFInit;
  theNP->Display = BDFDisplay;
  theNP->Execute = BDFExecute;

  /* nonlinear assemble functions */
  na->PreProcess                      = BDFPreProcess;
  na->PostProcess             = BDFPostProcess;
  na->NLAssembleSolution      = BDFAssembleSolution;
  na->NLAssembleDefect        = BDFAssembleDefect;
  na->NLAssembleMatrix        = BDFAssembleMatrix;

  /* and the time solver */
  ts->TimePreProcess              = TimePreProcess;
  ts->TimeInit                    = TimeInit;
  ts->TimeStep                    = TimeStep;
  ts->TimePostProcess             = TimePostProcess;

  return(0);
}


/****************************************************************************/
/*
   InitNewtonSolver  - Init this file

   SYNOPSIS:
   INT InitNewtonSolver (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

INT InitBDFSolver (void)
{
  if (CreateClass (T_SOLVER_CLASS_NAME ".bdf",
                   sizeof(NP_BDF), BDFConstruct))
    return (__LINE__);

  return (0);
}
