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
#include "debug.h"
#include "ugstruct.h"
#include "misc.h"
#include "gm.h"
#include "evm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"

#include "nls.h"
#include "ls.h"
#include "nls.h"
#include "assemble.h"
#include "transfer.h"
#include "error.h"
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

REP_ERR_FILE;

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
  INT nlinterpolate;                                             /* nonlinear interpolation			*/
  INT optnlsteps;                                            /* optimal number of nonlin. steps */
  INT rep;                               /* for repeat solver after grid changed */
  INT Break;                                                     /* break after error estimator         */
  INT Continue;                                              /* continue after error estimator  */
  INT copyall;                                               /* refine copy all                 */
  DOUBLE tstart;                                                 /* start time                          */
  DOUBLE dtstart;                                                /* time step to begin with			*/
  DOUBLE dtmin;                                                  /* smallest time step allowed		*/
  DOUBLE dtmax;                                                  /* largest time step allowed		*/
  DOUBLE dtscale;                                                /* scaling factor applied after ste*/
  DOUBLE rhogood;                                                /* threshold for step doubling		*/
  NP_TRANSFER *trans;                                            /* uses transgrid for nested iter  */
  NP_ERROR *error;                       /* error indicator                 */
  INT err_toplevel;                                              /* toplevel for error estimation	*/
  INT err_baselevel;                                             /* baselevel for error estimation	*/
  INT ctn;                                                               /* change to nested iteration		*/
  INT hist;
  INT list_i;
  INT list_n;
  DOUBLE list_dt[50];
  DOUBLE list_work[50];

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
  dcopy(NP_MG(ass),fl,tl,ALL_VECTORS,d,bdf->b);

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
  /* initialize bdf local variables */
  bdf->dt = bdf->dtstart;
  bdf->step = 0;
  bdf->t_0 = bdf->tstart;
  bdf->t_m1 = - bdf->dt;

  /* set initial values and boundary conditions in y_0 */
  *res = 1;
  if (tass->TAssembleInitial != NULL)
    if ( (*tass->TAssembleInitial)(tass,0,level,bdf->t_0,bdf->y_0,res) )
      return(1);
  if ( (*tass->TAssembleSolution)(tass,0,level,bdf->t_0,bdf->y_0,res) )
    return(1);

  /* write time to shell */
  sprintf(buffer,"%12.4lE",bdf->t_0);
  SetStringVar("TIME",buffer);
  SetStringVar(":BDF:TIME",buffer);
  sprintf(buffer,"%12.4lE",bdf->dt);
  SetStringVar("TIMESTEP",buffer);
  SetStringVar(":BDF:DT",buffer);
  SetStringVar(":BDF:SDT",buffer);

  /* statistics init */
  bdf->number_of_nonlinear_iterations = 0;
  bdf->total_linear_iterations = 0;
  bdf->max_linear_iterations = 0;
  bdf->exec_time = 0.0;

  /* init time-optimization-list */
  bdf->list_i = 0;
  bdf->list_n = 0;

  /* return ok */
  *res = 0;
  return(*res);
}

static INT TimeStep (NP_T_SOLVER *ts, INT level, INT *res)
{
  NP_BDF *bdf;
  NP_T_ASSEMBLE *tass;
  NP_NL_SOLVER *nlsolve;
  DOUBLE dt_p1,dt_0,g_p1,g_0,g_m1,qfm_dt,dtfactor;
  DOUBLE Factor[MAX_VEC_COMP];
  INT n_unk;
  INT i,k,mg_changed,changed,ret;
  INT low,llow,nlinterpolate,last_number_of_nonlinear_iterations;
  INT verygood,bad;
  NLRESULT nlresult;
  ERESULT eresult;
  MULTIGRID *mg;
  char buffer[128];
  static INT qfm;

  /* get numprocs ... */
  bdf = (NP_BDF *) ts;                  /* this is the trick, tsolver is derived	*/
  tass = bdf->tsolver.tass;             /* since we need to access it quite often       */
  nlsolve = ts->nlsolve;
  mg = nlsolve->base.mg;

  /* get number of components */
  n_unk = VD_NCOMP(bdf->y_p1);

  /* initialize strategy flags */
  verygood = bad = 0;
  eresult.step = 0.;

  /* compute solution at new time step */
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
    last_number_of_nonlinear_iterations =
      bdf->number_of_nonlinear_iterations;

    if (bdf->Continue) {
      nlinterpolate = bdf->nlinterpolate - 1;
      nlresult.converged = 1;
      goto Continue;
    }
    /* determine level where to predict to new time step */
    if (bdf->nested) low = MIN(bdf->baselevel,level);
    else low = level;

    /* grid adaption ? */
    if (bdf->error != NULL) nlinterpolate = bdf->nlinterpolate;
    else nlinterpolate = 0;

    /* predict to new time step on level low */
    dcopy(mg,0,low,ALL_VECTORS,bdf->y_p1,bdf->y_0);
    if (bdf->predictorder==1 && bdf->t_0>0.0)
    {
      dscal(mg,0,low,ALL_VECTORS,bdf->y_p1,(bdf->t_p1-bdf->t_m1)/dt_0);
      daxpy(mg,0,low,ALL_VECTORS,bdf->y_p1,1.0-(bdf->t_p1-bdf->t_m1)/dt_0,bdf->y_m1);
    }

    if ( (*tass->TAssemblePreProcess)(tass,0,level,
                                      bdf->t_p1,bdf->t_0,bdf->t_m1,
                                      bdf->y_p1,bdf->y_0,bdf->y_m1,res) )
      return(__LINE__);

    /* set Dirichlet conditions in predicted solution */
    if ( (*tass->TAssembleSolution)(tass,0,low,bdf->t_p1,bdf->y_p1,res) )
      return(__LINE__);

    /* do nested iteration on new time step */
    if (bdf->trans->PreProcessProject!=NULL)
      if ((*bdf->trans->PreProcessProject)(bdf->trans,0,level,res))
        return(__LINE__);
    for (llow=low,changed=0;;)
    {
      for (k=llow; k<=level; k++)
      {
        if (bdf->nested || changed) UserWriteF("Nested Iteration on level %d (%d)\n",k,level);

        /* prepare constant part of defect */
        dset(mg,0,k,ALL_VECTORS,bdf->b,0.0);
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
          if ( (*nlsolve->PreProcess)(nlsolve,k,bdf->y_p1,res) )
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
          if (!bdf->ctn || llow<level || bdf->baselevel>=llow)
          {
            /* halfen time step */
            bdf->dt *= 0.5;
            if (bdf->dt<bdf->dtmin)
            {
              UserWrite("time step too small -- aborting\n");
              return(__LINE__);
            }
            UserWrite("halfen time step\n");
            bad=1;

            /* restart optimization */
            bdf->list_i = 0;
            bdf->list_n = 0;
          }
          break;                                     /* and try all over again  */
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
          for (i=0; i<n_unk; i++) Factor[i] = 1.0;
          if ((*bdf->trans->InterpolateCorrection)(bdf->trans,k+1,bdf->y_p1,bdf->y_p1,NULL,Factor,res))
            return(__LINE__);
          /* set Dirichlet conditions in predicted solution */
          if ( (*tass->TAssembleSolution)(tass,k+1,k+1,bdf->t_p1,bdf->y_p1,res) )
            return(__LINE__);
        }
        else if (nlinterpolate > 0) {
          if (bdf->error->PreProcess != NULL)
            if ((*bdf->error->PreProcess)(bdf->error,level,res))
              NP_RETURN(1,res[0]);
          if (bdf->error->TimeError == NULL)
            NP_RETURN(1,res[0]);
#ifdef ModelP
          a_outervector_consistent(mg,0,level,bdf->y_p1);
          a_outervector_consistent(mg,0,level,bdf->y_0);
#endif
          if ((*bdf->error->TimeError)
                (bdf->error,level,bdf->t_p1,&dt_0,bdf->y_p1,bdf->y_0,
                ts,&eresult))
            NP_RETURN(1,res[0]);
          if (nlsolve->PostProcess!=NULL)
            if ( (*nlsolve->PostProcess)(nlsolve,k,bdf->y_p1,res) )
              NP_RETURN(1,res[0]);
          if (bdf->error->PostProcess != NULL)
            if ((*bdf->error->PostProcess)(bdf->error,level,res))
              NP_RETURN(1,res[0]);
          if (bdf->Break) return(0);
Continue:
          if (eresult.refine + eresult.coarse > 0)
          {
            if (bdf->err_toplevel>=0)
              for (i=bdf->err_toplevel; i<=TOPLEVEL(mg); i++)
                if (ClearMarksOnLevel(GRID_ON_LEVEL(mg,i),1)!=GM_OK)
                  return(__LINE__);
            if (bdf->err_baselevel>=0)
              for (i=0; i<=bdf->err_baselevel; i++)
                if (ClearMarksOnLevel(GRID_ON_LEVEL(mg,i),-1)!=GM_OK)
                  return(__LINE__);

            if (bdf->copyall) {
              if (RefineMultiGrid(mg,GM_COPY_ALL,
                                  GM_REFINE_PARALLEL,
                                  GM_REFINE_NOHEAPTEST) != GM_OK)
                NP_RETURN(1,res[0]);
            }
            else
            if (RefineMultiGrid(mg,GM_REFINE_TRULY_LOCAL,
                                GM_REFINE_PARALLEL,
                                GM_REFINE_NOHEAPTEST) != GM_OK)
              NP_RETURN(1,res[0]);
          }
          if (level != TOPLEVEL(mg)) {
            if (level < TOPLEVEL(mg)) {
              if (InterpolateVDAllocation(mg,bdf->y_m1))
                NP_RETURN(1,res[0]);
              if (InterpolateVDAllocation(mg,bdf->y_0))
                NP_RETURN(1,res[0]);
              if (InterpolateVDAllocation(mg,bdf->y_p1))
                NP_RETURN(1,res[0]);
              if (InterpolateVDAllocation(mg,bdf->b))
                NP_RETURN(1,res[0]);
            }
            level = TOPLEVEL(mg);
            mg_changed = 1;
          }
          else {
            mg_changed = 0;
            for (i=0; i<=level; i++)
              if (GSTATUS(GRID_ON_LEVEL(mg,i),GSTATUS_BDF))
              {
                RESETGSTATUS(GRID_ON_LEVEL(mg,i),GSTATUS_BDF);
                mg_changed = 1;
              }
          }
          if (mg_changed)
          {
            k = level - 1;
            if (bdf->trans->PreProcessSolution != NULL)
              if ((*bdf->trans->PreProcessSolution)
                    (bdf->trans,0,level,bdf->y_p1,res))
                NP_RETURN(1,res[0]);
            if ((*bdf->trans->InterpolateNewVectors)
                  (bdf->trans,0,level,bdf->y_m1,res))
              NP_RETURN(1,res[0]);
            if ((*bdf->trans->InterpolateNewVectors)
                  (bdf->trans,0,level,bdf->y_0,res))
              NP_RETURN(1,res[0]);
            if ((*bdf->trans->InterpolateNewVectors)
                  (bdf->trans,0,level,bdf->y_p1,res))
              NP_RETURN(1,res[0]);
            if (bdf->trans->PostProcessSolution != NULL)
              if ((*bdf->trans->PostProcessSolution)
                    (bdf->trans,0,level,bdf->y_p1,res))
                NP_RETURN(1,res[0]);
            nlinterpolate--;
            if(bdf->rep ==0) {
              k = level;
              nlinterpolate = 0;
              if (nlsolve->PostProcess!=NULL)
                if ((*nlsolve->PostProcess)
                      (nlsolve,k,bdf->y_p1,res))
                  NP_RETURN(1,res[0]);
            }
          }
          else {
            k = level;
            nlinterpolate = 0;
            if (nlsolve->PostProcess!=NULL)
              if ( (*nlsolve->PostProcess)(nlsolve,k,bdf->y_p1,res) )
                NP_RETURN(1,res[0]);
          }
        }
      }

      /* switch to nested? */
      if (!nlresult.converged && bdf->ctn && llow==level && bdf->baselevel<llow)
      {
        llow = bdf->baselevel;
        changed=1;
        UserWrite("Change to Nested Iteration\n");
      }
      else
        break;
    }
    if (bdf->trans->PostProcessProject!=NULL)
      if ((*bdf->trans->PostProcessProject)(bdf->trans,0,level,res))
        return(__LINE__);

    /* check convergence */
    if (!nlresult.converged) continue;                  /* start again with smaller time step   */
    else break;                                                                 /* exit while loop						*/
  }
  if ( (*tass->TAssemblePostProcess)(tass,0,level,
                                     bdf->t_p1,bdf->t_0,bdf->t_m1,
                                     bdf->y_p1,bdf->y_0,bdf->y_m1,res) )
    return(__LINE__);

  /* accept new time step */
  dcopy(mg,0,level,ALL_VECTORS,bdf->y_m1,bdf->y_0 );
  bdf->t_m1 = bdf->t_0;
  dcopy(mg,0,level,ALL_VECTORS,bdf->y_0 ,bdf->y_p1);
  bdf->t_0  = bdf->t_p1;
  bdf->step++;

  /* write time to shell */
  sprintf(buffer,"%12.4lE",bdf->t_0);
  SetStringVar("TIME",buffer);
  SetStringVar(":BDF:TIME",buffer);
  sprintf(buffer,"%12.4lE",bdf->dt);
  SetStringVar("TIMESTEP",buffer);
  SetStringVar(":BDF:DT",buffer);
  sprintf(buffer,"%12.4lE",bdf->t_0/bdf->exec_time);
  SetStringVar(":BDF:AEFF",buffer);
  sprintf(buffer,"%12.4lE",bdf->dt/nlresult.exec_time);
  SetStringVar(":BDF:EFF",buffer);

  UserWriteF("TIMESTEP %4d: TIME=%10.4lg DT=%10.4lg EXECT=%10.4lg NLIT=%5d LIT=%5d MAXLIT=%3d QFM=%d\n",
             bdf->step,bdf->t_0,bdf->dt,bdf->exec_time,bdf->number_of_nonlinear_iterations,
             bdf->total_linear_iterations,bdf->max_linear_iterations,qfm);

  /* chose new dt for next time step */
  if ((bdf->optnlsteps) && (nlresult.converged))
  {
    qfm = 0;
    if (bdf->hist>2)
    {
      /* fill in result */
      bdf->list_dt[bdf->list_i]    = bdf->dt;
      bdf->list_work[bdf->list_i] = nlresult.exec_time/bdf->dt;
      bdf->list_i = (bdf->list_i+1)%bdf->hist;
      bdf->list_n++;
      bdf->list_n = MIN(bdf->list_n,bdf->hist);

      /* get quadratic fitted dt for minimal work */
      ret = QuadraticFittedMin (bdf->list_dt,bdf->list_work,bdf->list_n,&qfm_dt);
      if (ret==0)
      {
        qfm_dt = MAX(0.5*bdf->dt,qfm_dt);
        qfm_dt = MIN(2.0*bdf->dt,qfm_dt);
        qfm_dt = MAX(bdf->dtmin,qfm_dt);
        qfm_dt = MIN(bdf->dtmax,qfm_dt);
        bdf->dt = qfm_dt;
        qfm = 1;
      }
    }

    if (!qfm)
    {
      k = bdf->number_of_nonlinear_iterations - last_number_of_nonlinear_iterations;
      if (k <= 0) bdf->dt *= 2.0;
      else
      {
        if (bdf->optnlsteps==1) dtfactor = QUOT(20+bdf->number_of_nonlinear_iterations,5+bdf->step);
        else dtfactor = bdf->optnlsteps;
        bdf->dt *= SQRT(QUOT(dtfactor,k));
      }
      if (bdf->dt < bdf->dtmin) bdf->dt = bdf->dtmin;
      else if (bdf->dt > bdf->dtmax) bdf->dt = bdf->dtmax;
      PRINTDEBUG(np,1,("new time step %f k %d n %d l %d\n",bdf->dt,k,bdf->number_of_nonlinear_iterations,last_number_of_nonlinear_iterations));
      *res = 0;
    }
  }
  else if (eresult.step ==0.)
  {
    if (verygood && (!bad) && bdf->dt*2<=bdf->dtmax)
    {
      bdf->dt *= 2.0;
      UserWrite("doubling time step\n");
      *res=0;
    }
    else if ((!bad) && bdf->dt*bdf->dtscale<=bdf->dtmax && bdf->dt*bdf->dtscale>=bdf->dtmin)
    {
      bdf->dt *= bdf->dtscale;
      *res=0;
    }
  }
  else if (eresult.step > bdf->dt && (!bad))
  {
    if (eresult.step <=bdf->dtmax)
    {
      bdf->dt = eresult.step;
      UserWrite("time step modified\n");
    }
    if (eresult.step >=bdf->dtmax)
    {
      bdf->dt = bdf->dtmax;
      UserWrite("time step modified, dt= dtmax\n");
    }
  }
  else if (eresult.step < bdf->dt)
  {
    if ( eresult.step < bdf->dtmin)
    {
      bdf->dt = bdf->dtmin;
      UserWrite("time step modified, \n WARNING dt = dtmin\n");
    }
    if ( eresult.step >= bdf->dtmin)
    {
      bdf->dt = eresult.step;
      UserWrite("time step modified\n");
    }
  }

  /* save suggested timestep */
  sprintf(buffer,"%12.4lE",bdf->dt);
  SetStringVar(":BDF:SDT",buffer);

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
  VECDATA_DESC *tmp;
  INT r;

  /* get numprocs ... */
  bdf = (NP_BDF *) base;

  /* call tsolver init */
  r = NPTSolverInit(&bdf->tsolver,argc,argv);

  /* read data descs */
  bdf->y_0 = bdf->tsolver.y;       /* allocated already in tsolver */
  tmp = ReadArgvVecDesc(base->mg,"yp1",argc,argv);
  if (tmp!=NULL) bdf->y_p1 = tmp;
  tmp = ReadArgvVecDesc(base->mg,"ym1",argc,argv);
  if (tmp!=NULL) bdf->y_m1 = tmp;
  tmp = ReadArgvVecDesc(base->mg,"b",argc,argv);
  if (tmp!=NULL) bdf->b = tmp;

  /* read other numprocs */
  bdf->trans = (NP_TRANSFER *) ReadArgvNumProc(base->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  if (bdf->trans == NULL) return(NP_NOT_ACTIVE);
  bdf->error = (NP_ERROR *) ReadArgvNumProc(base->mg,"E",ERROR_CLASS_NAME,argc,argv);
  if (bdf->error == NULL)
  {
    UserWrite("no indicator active\n");
  }

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
  if (ReadArgvINT("ctn",&(bdf->ctn),argc,argv))
  {
    UserWrite("default: change to nested: OFF\n");
    bdf->ctn=0;
  }
  if ((bdf->nested<0)||(bdf->nested>1)) return(NP_NOT_ACTIVE);
  if (ReadArgvINT("optnlsteps",&(bdf->optnlsteps),argc,argv)) bdf->optnlsteps = 0;
  if (bdf->optnlsteps < 0) return(NP_NOT_ACTIVE);
  if (bdf->optnlsteps==1)
  {
    bdf->hist = 4;
  }
  else
  {
    if (ReadArgvINT("hist",&(bdf->hist),argc,argv)) bdf->hist = 0;
  }
  if (bdf->hist < 0 || bdf->hist > 50) return(NP_NOT_ACTIVE);
  if (ReadArgvINT("rep",&(bdf->rep),argc,argv))
    bdf->rep=1;
  if (ReadArgvINT("nlinterpolate",&(bdf->nlinterpolate),argc,argv))
    bdf->nlinterpolate=0;
  if (bdf->nlinterpolate<0) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("tstart",&(bdf->tstart),argc,argv))
    bdf->tstart = 0.0;
  if (ReadArgvDOUBLE("dtstart",&(bdf->dtstart),argc,argv))
  {
    UserWrite("dtstart must be specified\n");
    return(NP_NOT_ACTIVE);
  }
  if ((bdf->dtstart<0.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("dtmin",&(bdf->dtmin),argc,argv))
  {
    bdf->dtmin = bdf->dtstart;
    return(NP_NOT_ACTIVE);
  }
  if ((bdf->dtmin<0.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("dtmax",&(bdf->dtmax),argc,argv))
  {
    bdf->dtmax = bdf->dtstart;
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

  if (bdf->error==NULL || ReadArgvINT("etl",&(bdf->err_toplevel),argc,argv))
  {
    bdf->err_toplevel = -1;
  }

  if (bdf->error==NULL || ReadArgvINT("ebl",&(bdf->err_baselevel),argc,argv))
  {
    bdf->err_baselevel = -1;
  }
  if (bdf->err_toplevel>=0 && bdf->err_baselevel>=0 && bdf->err_baselevel>bdf->err_toplevel) return(NP_NOT_ACTIVE);

  bdf->copyall = ReadArgvOption("copyall",argc,argv);

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

  UserWrite("\nBDF data:\n");
  if (bdf->trans != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(bdf->trans));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (bdf->error != NULL)
  {
    UserWriteF(DISPLAY_NP_FORMAT_SS,"E",ENVITEM_NAME(bdf->error));
    UserWriteF(DISPLAY_NP_FORMAT_SI,"err_toplevel",
               (int)bdf->err_toplevel);
    UserWriteF(DISPLAY_NP_FORMAT_SI,"err_baselevel",
               (int)bdf->err_baselevel);
    UserWriteF(DISPLAY_NP_FORMAT_SI,"copyall",
               (int)bdf->copyall);
  }
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"E","---");
  UserWriteF(DISPLAY_NP_FORMAT_SF,"t_m1",(float)bdf->t_m1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"t_0",(float)bdf->t_0);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"t_p1",(float)bdf->t_p1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"dt",(float)bdf->dt);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"dtmin",(float)bdf->dtmin);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"dtmax",(float)bdf->dtmax);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"nested",(int)bdf->nested);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"ctn",(int)bdf->ctn);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"nlinterpolate",(int)bdf->nlinterpolate);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"optnlsteps",(int)bdf->optnlsteps);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"hist",(int)bdf->hist);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"dtscale",(float)bdf->dtscale);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"rhogood",(float)bdf->rhogood);
  if (bdf->y_p1 != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"y_p1",ENVITEM_NAME(bdf->y_p1));
  if (bdf->y_0  != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"y_0 ",ENVITEM_NAME(bdf->y_0 ));
  if (bdf->y_m1 != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"y_m1",ENVITEM_NAME(bdf->y_m1));
  if (bdf->b    != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"b   ",ENVITEM_NAME(bdf->b   ));

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
  DOUBLE initialtime,dtime;
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
  bdf->Break = ReadArgvOption("Break",argc,argv);
  bdf->Continue = ReadArgvOption("Continue",argc,argv);

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
  if (ReadArgvOption("init",argc,argv))
    if (np->TimeInit != NULL) {
      if (ReadArgvDOUBLE("t",&initialtime,argc,argv) == 0)
        bdf->tstart = initialtime;
      if (ReadArgvDOUBLE("dt",&dtime,argc,argv) == 0)
        bdf->dtstart = dtime;
      if ((*np->TimeInit)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeInit failed, error code %d\n",result);
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
  NP_BDF *bdf;

  /* get numprocs ... */
  ts = (NP_T_SOLVER *) theNP;
  na = &ts->nlass;
  bdf = (NP_BDF *) theNP;

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

  bdf->y_p1 = NULL;
  bdf->y_m1 = NULL;
  bdf->b = NULL;

  return(0);
}

/****************************************************************************/
/*D
   bdf - BDF time stepping scheme

   DESCRIPTION:
   This numproc executes a BDF time stepping scheme of first and second order.
   In every time step it calls a nonlinear solver, where the linearization
   in the nonlinear solver is constructed from a time dependent assembling
   routine of type 'NP_T_ASSEMBLE'.

   .vb
   npinit <name> $y <sol> $A <tnlass> $S <nlsolver>
          $T <transfer> [$baselevel <bl>] $order {1|2}
                  $predictorder {0|1} $nested {0|1}
          $tstart <time> $dtstart <step> $dtmin <min step> $dtmax <max step>
          $dtscale <fac> $rhogood <rho> [$J <mat>];
   .ve

   .  $y~<sol>            - solution vector
   .  $A~<tnlass>         - assemble numproc of type 'NP_T_ASSEMBLE'
   .  $S~<nlsolver>       - nonlinear solver numproc of type 'NP_NL_SOLVER'
   .  $T~<transfer>       - transfer numproc
   .  $baselevel~<bl>     - baselevel for nested iteration
   .  $order~{1|2}        - order of the time stepping scheme
   .  $predictorder~{0|1} - order of the predictor for the next start vector
   .  $nested~{0|1}       - flag nested nonlinear solver in every time step
   .  $tstart~<time>      - start time
   .  $dtstart~<step>     - time step to begin with
   .  $dtmin~<min step>   - smallest time step allowed
   .  $dtmax~<max step>   - largest time step allowed
   .  $dtscale~<fac>      - scaling factor applied after step
   .  $rhogood~<rho>      - threshold for step doubling

   'npexecute <name> [$pre] [$init] [$bdf1|$bdf1n|$bdf2|$bdf2n] [$post];'

   .  $pre - preprocess
   .  $init - initialization
   .  $bdf1 - backward Euler of first order
   .  $bdf1n - nested backward Euler of first order
   .  $bdf1 - backward Euler of second order
   .  $bdf1n - nested backward Euler of second order
   .  $post - postprocess

   EXAMPLE:
   .vb
 # grid transfer numproc
   npcreate transfer $c transfer;
   npinit transfer;

 # assemble numproc
   npcreate tnlass $c boxtnl;
   npinit tnlass;

 # linear solver and iteration numprocs
   npcreate smooth $c gs;
   npinit smooth;

   npcreate basesolver $c ls;
   npinit basesolver $red 0.001 $m 10 $I smooth;

   npcreate lmgc $c lmgc;
   npinit lmgc $S smooth smooth basesolver $T transfer $n1 2 $n2 2;

   npcreate mgs $c ls;
   npinit mgs $m 40 $I lmgc $display full;

 # nonlinear solver numproc to be used by time solver
   npcreate newton $c newton;
   npinit newton $red 1.0E-8 $T transfer $S mgs
              $rhoreass 0.8 $lsteps 6 $maxit 50 $line 1 $linrate 0
              $lambda 1.0 $linminred 1.0E-2 $display full;

 # the time solver
   npcreate ts $c bdf;
   npinit ts $y sol $A tnlass $S newton
          $T transfer $baselevel 0 $order 1 $predictorder 0 $nested 0
          $dtstart 0.00125 $dtmin 0.001 $dtmax 0.00125 $dtscale 1.0
          $rhogood 0.001;

 # first time step
   TIMESTEPS = 100;
   step = 1;
   npexecute ts $pre $init;
   npexecute ts $bdf1;

 # time step loop
   repeat {
    if (step==TIMESTEPS) break;
    step = step+1;
    npexecute ts $bdf2;
   }
   npexecute ts $post;
   .ve
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   InitBDFSolver  - Init this file

   SYNOPSIS:
   INT InitBDFSolver (void);

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
  INT error=0;

  if (MakeStruct(":BDF")!=0)
  {
    SetHiWrd(error,__LINE__);
    return (1);
  }

  if (CreateClass (T_SOLVER_CLASS_NAME ".bdf",sizeof(NP_BDF), BDFConstruct)) return (__LINE__);

  return (0);
}
