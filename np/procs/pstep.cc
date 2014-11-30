// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pstep.c                                                       */
/*                                                                          */
/* Purpose:   parameter-step schemes                                        */
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

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <config.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include "ugdevices.h"
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

#include "pstep.h"
#include "reinit.h"

USING_UG_NAMESPACES

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
/*  Class Definition                                                        */
/*                                                                          */
/****************************************************************************/

typedef struct
{
  NP_P_STEP pstep;                       /* derived from class NP_P_STEP    */

  /* nonlinear assemmlies used */
  NP_NL_ASSEMBLE nlass;
  NP_ENL_ASSEMBLE enlass;

  /* parameters (to be set with init function) */
  INT n;                                 /* number of parameters            */
  char name[EXTENSION_MAX][NAMESIZE];    /* names of parameters             */
  INT baselevel;                         /* for nested iteration            */
  INT nested;                            /* use nested iteration            */
  INT displayMode;

  /* references to numprocs  */
  NP_TRANSFER *trans;                    /* uses transgrid for nested iter  */
  NP_T_ASSEMBLE *tass;                   /* time assemble numproc			*/
  NP_NL_SOLVER *nlsolve;                 /* nonlinear solver numproc		*/
  NP_ENL_SOLVER *enlsolve;               /* extended nonlinear solver np	*/
  NP_REINIT *reinit;                         /* reinit numproc                          */

  /* internals */
  INT n_step;                                /* number of step                  */
  DOUBLE last_nls_dp[EXTENSION_MAX];         /* initial parameter difference    */
  DOUBLE last_nls_nt;                        /* inital norm of sol_t            */
  DOUBLE last_enls_scale;                    /* last sucessfull scale           */
  INT last_solve_mode;                       /* 0: init                         */
  /* >0: nls in sequence             */
  /* <0: enls in sequence            */
  EVECDATA_DESC *sol_t;                      /* tangential                      */
  DOUBLE scale;
  DOUBLE raster;                         /* if >0, std newton is forced to  */
                                         /* raster                          */
  DOUBLE raster_size[EXTENSION_MAX];     /* raster sizes                    */

} NP_SPS;                                /* simple pstep                    */

/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Nonlinear Assemble Interface provided to nonlinear solver				*/
/*                                                                          */
/****************************************************************************/

#define TIME_INFTY                      1e50
static NP_P_STEP *current_pstep;

static INT SPS_NLPreProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, INT *res)
{
  return(0);
}

static INT SPS_NLAssembleSolution (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, INT *res)
{
  INT i;
  NP_SPS *sps;
  NP_T_ASSEMBLE *tass;
  NP_REINIT *reinit;
  REINIT_RESULT reinit_res;

  /* get numprocs */
  sps = (NP_SPS*)current_pstep;
  tass = sps->tass;
  reinit=sps->reinit;

  /* reinit */
  for (i=0; i<sps->pstep.sol_p0->n; i++)
    if ((*reinit->ReinitProblem)(reinit,sps->name[i],EVDD_E(sps->pstep.sol_p1,tl,i),&reinit_res)) NP_RETURN(1,*res);

  return((*tass->TAssembleSolution)(tass,fl,tl,TIME_INFTY,u,res));
}

static INT SPS_NLAssembleDefect (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  INT i;
  NP_SPS *sps;
  NP_T_ASSEMBLE *tass;
  NP_REINIT *reinit;
  REINIT_RESULT reinit_res;

  /* get tass numproc */
  sps = (NP_SPS*)current_pstep;
  tass = sps->tass;
  reinit=sps->reinit;

  /* reinit */
  for (i=0; i<sps->pstep.sol_p0->n; i++)
    if ((*reinit->ReinitProblem)(reinit,sps->name[i],EVDD_E(sps->pstep.sol_p1,tl,i),&reinit_res)) NP_RETURN(1,*res);

  /* assemble defect */
  dset(NP_MG(sps),fl,tl,ALL_VECTORS,d,0.0);
  if ((*tass->TAssembleDefect)(tass,fl,tl,0.0,-1.0,0.0,current_pstep->sol_p1->vd,d,NULL,res)) NP_RETURN(1,*res);
  return ((*tass->TAssembleDefect)(tass,fl,tl,TIME_INFTY,1.0,-TIME_INFTY,u,d,J,res));
}

static INT SPS_NLAssembleMatrix (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *J, INT *res)
{
  INT i;
  NP_SPS *sps;
  NP_T_ASSEMBLE *tass;
  NP_REINIT *reinit;
  REINIT_RESULT reinit_res;

  /* get tass numproc */
  sps = (NP_SPS*)current_pstep;
  tass = sps->tass;
  reinit=sps->reinit;

  /* reinit */
  for (i=0; i<sps->pstep.sol_p0->n; i++)
    if ((*reinit->ReinitProblem)(reinit,sps->name[i],EVDD_E(sps->pstep.sol_p1,tl,i),&reinit_res)) NP_RETURN(1,*res);

  return ((*tass->TAssembleMatrix)(tass,fl,tl,TIME_INFTY,-TIME_INFTY,u,d,v,J,res));
}

static INT SPS_NLPostProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  return(0);
}

/****************************************************************************/
/*                                                                          */
/* Extended nonlinear Assemble Interface provided to extended nonlin solver */
/*                                                                          */
/****************************************************************************/

static INT SPS_ENLPreProcess (NP_ENL_ASSEMBLE *ass, INT fl, INT tl, EVECDATA_DESC *x, INT *res)
{
  return(0);
}

static INT SPS_ENLAssembleSolution (NP_ENL_ASSEMBLE *ass, INT fl, INT tl, EVECDATA_DESC *u, INT *res)
{
  NP_T_ASSEMBLE *tass;

  tass = ((NP_SPS*)current_pstep)->tass;
  return((*tass->TAssembleSolution)(tass,fl,tl,TIME_INFTY,u->vd,res));
}

static INT SPS_ENLAssembleDefect (NP_ENL_ASSEMBLE *ass, INT fl, INT tl, EVECDATA_DESC *u, EVECDATA_DESC *d, EMATDATA_DESC *J, INT *res)
{
  INT i;
  NP_SPS *sps;
  NP_T_ASSEMBLE *tass;
  DOUBLE r;
  VECDATA_DESC *x,*p,*t;
  NP_REINIT *reinit;
  REINIT_RESULT reinit_res;

  /* get tass numproc */
  sps = (NP_SPS*)current_pstep;
  tass = sps->tass;
  reinit = sps->reinit;
  assert(u->n==1);
  x=u->vd;
  p=current_pstep->sol_p0->vd;
  t=sps->sol_t->vd;

  /* reinit problem */
  for (i=0; i<u->n; i++)
    if ((*reinit->ReinitProblem)(reinit,sps->name[i],EVDD_E(u,tl,i),&reinit_res)) return(1);

  /* assemble defect */
  dcopy(NP_MG(sps),fl,tl,ALL_VECTORS,d->vd,x);
  daxpy(NP_MG(sps),fl,tl,ALL_VECTORS,d->vd,-1.0,p);
  ddot(NP_MG(sps),fl,tl,ON_SURFACE,d->vd,t,&r);
  EVDD_E(d,tl,0)=r-sps->scale*sps->last_nls_nt;

  dset(NP_MG(sps),fl,tl,ALL_VECTORS,d->vd,0.0);
  if ((*tass->TAssembleDefect)(tass,fl,tl,0.0,-1.0,0.0,current_pstep->sol_p1->vd,d->vd,NULL,res)) { *res=1; REP_ERR_RETURN(1); }   /* NP_RETURN(1,*res);*/
  if ((*tass->TAssembleDefect)(tass,fl,tl,TIME_INFTY,1.0,-TIME_INFTY,u->vd,d->vd,J->mm,res))              { *res=1; REP_ERR_RETURN(1); }       /* NP_RETURN(1,*res);*/

  return(0);
}

#define PARAMETER_EPS                   1e-08

#ifdef ModelP
static INT l_vector_makeinconsistent_pstep (GRID *g, const VECDATA_DESC *x)
{
  VECTOR *v;
  INT vc,i,type,mask,n;
  const SHORT *Comp;

  if (VD_IS_SCALAR(x)) {
    mask = VD_SCALTYPEMASK(x);
    vc = VD_SCALCMP(x);
    for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
      if ((mask & VDATATYPE(v)) &&
          (DDD_InfoPriority(PARHDR(v)) != PrioMaster))
        VVALUE(v,vc) = 0.0;
  }
  else
    for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v)) {
      type = VTYPE(v);
      n = VD_NCMPS_IN_TYPE(x,type);
      if (n == 0) continue;
      if (DDD_InfoPriority(PARHDR(v)) == PrioMaster) continue;
      Comp = VD_CMPPTR_OF_TYPE(x,type);
      for (i=0; i<n; i++)
        VVALUE(v,Comp[i]) = 0.0;
    }

  return(NUM_OK);
}
#endif

static INT SPS_ENLAssembleMatrix (NP_ENL_ASSEMBLE *ass, INT fl, INT tl, EVECDATA_DESC *u, EVECDATA_DESC *d, EVECDATA_DESC *v, EMATDATA_DESC *J, INT *res)
{
  INT i,j,level;
  NP_SPS *sps;
  NP_T_ASSEMBLE *tass;
  NP_REINIT *reinit;
  REINIT_RESULT reinit_res;

  /* get tass numproc */
  sps = (NP_SPS*)current_pstep;
  reinit = sps->reinit;
  tass = sps->tass;
  assert(sps->n==u->n);

  /* assemble matrix */
  for (level=fl; level<=tl; level++) for (i=0; i<u->n; i++) for (j=0; j<u->n; j++) EMDD_EE(J,level,i*u->n+j)=0.0;
  dset(NP_MG(sps),fl,tl,ALL_VECTORS,J->em[0],0.0);
  for (i=0; i<u->n; i++)
    if ((*reinit->ReinitProblem)(reinit,sps->name[i],EVDD_E(u,tl,i),&reinit_res)) return(1);
  if ((*tass->TAssembleDefect)(tass,fl,tl,0.0,-1.0,0.0,u->vd,J->em[0],NULL,res)) NP_RETURN(1,*res);
  if ((*tass->TAssembleDefect)(tass,fl,tl,TIME_INFTY,1.0,-TIME_INFTY,u->vd,J->em[0],NULL,res)) NP_RETURN(1,*res);
  for (i=0; i<u->n; i++)
  {
    dset(NP_MG(sps),fl,tl,ALL_VECTORS,J->me[i],0.0);
    if ((*reinit->ReinitProblem)(reinit,sps->name[i],(1.0+PARAMETER_EPS)*EVDD_E(u,tl,i),&reinit_res)) return(1);
    if ((*tass->TAssembleDefect)(tass,fl,tl,0.0,-1.0,0.0,u->vd,J->me[i],NULL,res)) return(1);
    if ((*tass->TAssembleDefect)(tass,fl,tl,TIME_INFTY,1.0,-TIME_INFTY,u->vd,J->me[i],NULL,res)) NP_RETURN(1,*res);
    if ((*reinit->ReinitProblem)(reinit,sps->name[i],EVDD_E(u,tl,i),&reinit_res)) return(1);
    if (daxpy(NP_MG(sps),fl,tl,ALL_VECTORS,J->me[i],-1.0,J->em[0])) return(1);
    if (dscal(NP_MG(sps),fl,tl,ALL_VECTORS,J->me[i],1.0/PARAMETER_EPS/EVDD_E(u,tl,i))) return(1);
  }
  for (i=0; i<u->n; i++)
  {
    if (dcopy(NP_MG(sps),fl,tl,ALL_VECTORS,J->em[i],sps->sol_t->vd)) return(1);
#ifdef ModelP
    for (j=fl; j<tl; j++)
      if (l_vector_makeinconsistent_pstep(GRID_ON_LEVEL(NP_MG(sps),j),J->me[i])) return(1);
#endif
  }
  return ((*tass->TAssembleMatrix)(tass,fl,tl,TIME_INFTY,-TIME_INFTY,u->vd,d->vd,v->vd,J->mm,res));
}

static INT SPS_ENLPostProcess (NP_ENL_ASSEMBLE *ass, INT fl, INT tl, EVECDATA_DESC *x, EVECDATA_DESC *d, EMATDATA_DESC *J, INT *res)
{
  return(0);
}

/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* Implementation of NP_P_STEP functions                                                                */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/

static INT SPS_PreProcess (NP_P_STEP *pstep, INT level, EVECDATA_DESC *sol_p0, INT *res)
{
  MULTIGRID *mg;
  NP_SPS *sps;
  NP_REINIT *reinit;
  INT i;

  sps = (NP_SPS *)pstep;
  reinit = sps->reinit;
  mg = pstep->base.mg;

  /* internals */
  sps->n_step=0;
  sps->last_enls_scale=2.0;
  sps->last_solve_mode=0;
  for (i=0; i<sol_p0->n; i++)
    if ((*reinit->GetProblemParameter)(reinit,sps->name[i],&(EVDD_E(sol_p0,level,i)))) return(1);
  if (AllocEVDFromEVD(mg,0,level,sol_p0,&(sps->sol_t))) return(1);

  return(0);
}

static INT SPS_Step (NP_P_STEP *pstep, INT level, EVECDATA_DESC *sol_p0, EVECDATA_DESC *sol_p1, PSTEP_RESULT *pstep_result)
{
  NP_SPS *sps;
  NP_T_ASSEMBLE *tass;
  NP_NL_SOLVER *nlsolve;
  NP_ENL_SOLVER *enlsolve;
  NP_TRANSFER *trans;
  MULTIGRID *mg;
  INT i,ncomp,flevel,k,res;
  NLRESULT nlresult;
  ENLRESULT enlresult;
  DOUBLE Factor[MAX_VEC_COMP];
  INT do_nls,n0,n1;
  char text[128];
  DOUBLE s,nt,diff;

  /* check */
  if (sol_p0->n!=sol_p1->n) return(1);
  if (sol_p0->n!=1) return(1);

  /* get numprocs */
  sps = (NP_SPS *)pstep;
  tass = sps->tass;
  trans = sps->trans;
  nlsolve = sps->nlsolve;
  enlsolve = sps->enlsolve;
  mg = nlsolve->base.mg;

  /* set EVECDATA_DESC for debugging */
  sps->pstep.sol_p0 = sol_p0;
  sps->pstep.sol_p1 = sol_p1;

  /* get number of components */
  ncomp = VD_NCOMP(sol_p0->vd);

  /* reset pstep_result */
  pstep_result->converged = 1;
  pstep_result->number_of_nonlinear_iterations = 0;
  pstep_result->number_of_linear_iterations = 0;
  pstep_result->max_linear_iterations = 0;

  /* determine level where to predict to new time step */
  if (sps->nested) flevel = MIN(sps->baselevel,level);
  else flevel = level;

  /* prepare transfer */
  if (trans->PreProcessProject!=NULL)
    if ((*trans->PreProcessProject)(trans,0,level,&res)) NP_RETURN(1,res);

  /* decide solution scheme */
  if (sps->n_step==0)
    do_nls=1;
  else if (sps->last_solve_mode>0)
    do_nls=1;
  else if (sps->last_solve_mode<0 && sps->last_solve_mode>-5)
    do_nls=0;
  else
    do_nls=1;

  /* assemble */
  current_pstep = pstep;
  if ((*tass->TAssemblePreProcess)(tass,0,level,TIME_INFTY,0.0,0.0,sol_p1->vd,sol_p0->vd,NULL,&res)) NP_RETURN(1,res);
  if ((*tass->TAssembleSolution)(tass,0,flevel,TIME_INFTY,sol_p1->vd,&res)) NP_RETURN(1,res);

  /* do (nested) iteration on new parameter step */
  for (k=flevel; k<=level; k++)
  {
    if (sps->nested) UserWriteF("Nested Iteration on level %d (%d)\n",k,level);

    if (do_nls)
    {
      /* predict to new step on level low */
      if (k==flevel)
      {
        dcopy(mg,0,flevel,ALL_VECTORS,sol_p1->vd,sol_p0->vd);
        if (sps->n_step>0)
        {
          for (i=0; i<sps->pstep.sol_p0->n; i++)
          {
            if (EVDD_E(sps->sol_t,level,i)*sps->last_nls_dp[i]>=0.0)
              EVDD_E(sol_p1,level,i)=EVDD_E(sps->pstep.sol_p0,level,i)+sps->last_nls_dp[i];
            else
              EVDD_E(sol_p1,level,i)=EVDD_E(sps->pstep.sol_p0,level,i)-sps->last_nls_dp[i];
            if (sps->raster>0.0 && EVDD_E(sol_p1,level,i)>EVDD_E(sps->pstep.sol_p0,level,i))
            {
              n0=(INT)floor(EVDD_E(sps->pstep.sol_p0,level,i)/sps->raster_size[i]);
              n1=(INT)floor(EVDD_E(sol_p1,level,i)/sps->raster_size[i]+0.5);
              if (n0>=n1) n1=n0+1;
              EVDD_E(sol_p1,level,i)=n1*sps->raster_size[i];
            }
            else if (sps->raster>0.0 && EVDD_E(sol_p1,level,i)<EVDD_E(sps->pstep.sol_p0,level,i))
            {
              n0=(INT)ceil(EVDD_E(sps->pstep.sol_p0,level,i)/sps->raster_size[i]);
              n1=(INT)ceil(EVDD_E(sol_p1,level,i)/sps->raster_size[i]+0.5);
              if (n0<=n1) n1=n0-1;
              EVDD_E(sol_p1,level,i)=n1*sps->raster_size[i];
            }
          }
        }
      }

      /* prepare nonlinear solver on level k */
      if (nlsolve->PreProcess!=NULL)
        if ((*nlsolve->PreProcess)(nlsolve,k,sol_p1->vd,&res)) NP_RETURN(1,res);

      /* solve nonlinear on level k */
      if ((*nlsolve->Solver)
            (nlsolve,k,sol_p1->vd,&sps->nlass,nlsolve->abslimit,nlsolve->reduction,&nlresult))
        return(__LINE__);

      /* postprocess nonlinear solver on level k */
      if (nlsolve->PostProcess!=NULL)
        if ((*nlsolve->PostProcess)(nlsolve,k,sol_p1->vd,&res)) NP_RETURN(1,res);

      /* update result */
      pstep_result->number_of_nonlinear_iterations += nlresult.number_of_nonlinear_iterations;
      if (nlresult.number_of_nonlinear_iterations==0) pstep_result->jumped=1;
      else pstep_result->jumped=0;
      pstep_result->number_of_linear_iterations += nlresult.total_linear_iterations;
      pstep_result->max_linear_iterations = MAX(pstep_result->max_linear_iterations,nlresult.max_linear_iterations);
      if (!nlresult.converged) { do_nls=0; sps->last_solve_mode=0; }
    }
    if (!do_nls)
    {
      for (sps->scale=MIN(2.0*sps->last_enls_scale,1.0); sps->scale>1e-6; sps->scale*=0.5)
      {
        /* predict to new step on level low */
        if (k==flevel)
        {
          s=0.9*sps->scale;
          if (dcopy(mg,0,flevel,ALL_VECTORS,sol_p1->vd,sol_p0->vd)) return(1);
          if (daxpy(mg,0,flevel,ALL_VECTORS,sol_p1->vd,sps->last_nls_nt*s,sps->sol_t->vd)) return(1);
          for (i=0; i<sps->pstep.sol_p0->n; i++)
            EVDD_E(sol_p1,flevel,i)=EVDD_E(sps->pstep.sol_p0,flevel,i)+s*EVDD_E(sps->sol_t,flevel,i);
        }

        /* prepare nonlinear solver on level k */
        if (enlsolve->PreProcess!=NULL)
          if ((*enlsolve->PreProcess)(enlsolve,k,sol_p1,&res)) NP_RETURN(1,res);

        /* solve nonlinear on level k */
        if ((*enlsolve->Solver)
              (enlsolve,k,sol_p1,&sps->enlass,enlsolve->abslimit,enlsolve->reduction,&enlresult))
          return(__LINE__);

        /* postprocess nonlinear solver on level k */
        if (enlsolve->PostProcess!=NULL)
          if ((*enlsolve->PostProcess)(enlsolve,k,sol_p1,&res)) NP_RETURN(1,res);

        /* update result */
        pstep_result->number_of_nonlinear_iterations += enlresult.number_of_nonlinear_iterations;
        if (enlresult.number_of_nonlinear_iterations==0) pstep_result->jumped=1;
        else pstep_result->jumped=0;
        pstep_result->number_of_linear_iterations += enlresult.total_linear_iterations;
        pstep_result->max_linear_iterations = MAX(pstep_result->max_linear_iterations,enlresult.max_linear_iterations);
        if (enlresult.converged) break;
      }
      if (!enlresult.converged)
      {
        pstep_result->converged = 0;
        break;
      }
    }

    /* are we ready ? */
    if (k==level) break;

    /* interpolate up */
    for (i=0; i<ncomp; i++) Factor[i] = 1.0;
    if ((*trans->InterpolateCorrection)(trans,k+1,sol_p1->vd,sol_p1->vd,NULL,Factor,&res)) NP_RETURN(1,res);

    /* set Dirichlet conditions in predicted solution */
    if ((*tass->TAssembleSolution)(tass,k+1,k+1,TIME_INFTY,sol_p1->vd,&res)) NP_RETURN(1,res);
  }

  /* postprocess assemble */
  if ((*tass->TAssemblePostProcess)(tass,0,level,TIME_INFTY,0.0,0.0,sol_p1->vd,sol_p0->vd,NULL,&res)) NP_RETURN(1,res);

  /* postprocess transfer */
  if (trans->PostProcessProject!=NULL)
    if ((*trans->PostProcessProject)(trans,0,level,&res)) NP_RETURN(1,res);

  if (sps->displayMode == PCR_RED_DISPLAY || sps->displayMode == PCR_FULL_DISPLAY )
  {
    UserWriteF("\n");
    CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(sps),'%',"\n"); UserWrite(text);
    for (i=0; i<sps->pstep.sol_p0->n; i++)
    {
      UserWriteF("%17s:   %13e\n",sps->name[0],EVDD_E(sol_p1,level,0));
      UserWriteF("             step:   %13e\n",EVDD_E(sol_p1,level,0)-EVDD_E(sol_p0,level,0));
      if (do_nls==0)
        UserWriteF("            scale:   %13e\n",sps->scale);
    }
    UserWriteF("\n");
  }
  if (sps->displayMode == PCR_FULL_DISPLAY)
  {
    UserWriteF("               nl: %4d it.\n",(int)pstep_result->number_of_nonlinear_iterations);
    UserWriteF("              lin: %4d it.\n",(int)pstep_result->number_of_linear_iterations);
    UserWriteF("       max lin/nl: %4d it.\n",(int)pstep_result->max_linear_iterations);
    UserWriteF("\n");
  }

  /* update internals */
  if (pstep_result->converged)
  {
    dcopy(mg,0,level,ALL_VECTORS,sps->sol_t->vd,sps->pstep.sol_p1->vd);
    daxpy(mg,0,level,ALL_VECTORS,sps->sol_t->vd,-1.0,sps->pstep.sol_p0->vd);
    dnrm2(mg,0,level,ON_SURFACE,sps->sol_t->vd,&nt);
    assert(nt!=0.0);
    dscal(mg,0,level,ALL_VECTORS,sps->sol_t->vd,1.0/nt);
    for (i=0; i<sps->pstep.sol_p0->n; i++)
    {
      diff=EVDD_E(sps->pstep.sol_p1,level,i)-EVDD_E(sps->pstep.sol_p0,level,i);
      if (diff*EVDD_E(sps->sol_t,level,i)<=0.0) sps->last_solve_mode=0;
      EVDD_E(sps->sol_t,level,i)=diff;
    }
    if (do_nls==1)
    {
      for (i=0; i<sps->pstep.sol_p0->n; i++)
        sps->last_nls_dp[i]=EVDD_E(sps->sol_t,level,i);
      sps->last_nls_nt=nt;
      if (sps->last_solve_mode>0) sps->last_solve_mode++;
      else sps->last_solve_mode=1;
    }
    else
    {
      sps->last_enls_scale=sps->scale;
      if (sps->scale<1.0) sps->last_solve_mode=0;
      if (sps->last_solve_mode<0) sps->last_solve_mode--;
      else sps->last_solve_mode=-1;
    }
    sps->n_step++;
  }

  return(0);
}

static INT SPS_PostProcess (NP_P_STEP *pstep, INT level, INT *res)
{
  MULTIGRID *mg;
  NP_SPS *sps;

  sps = (NP_SPS *)pstep;
  mg = pstep->base.mg;

  /* internals */
  if (FreeEVD(mg,0,level,sps->sol_t)) REP_ERR_RETURN(1);

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

static INT SPS_Init (NP_BASE *base, INT argc, char **argv)
{
  NP_SPS *sps;
  VECDATA_DESC *tmp;
  INT i,ret;
  MULTIGRID *mg;

  /* get numprocs ... */
  sps = (NP_SPS *)base;
  mg = base->mg;

  /* read other numprocs */
  sps->tass = (NP_T_ASSEMBLE *) ReadArgvNumProc(mg,"A",T_ASSEMBLE_CLASS_NAME,argc,argv);
  if (sps->tass==NULL) return (NP_NOT_ACTIVE);
  sps->nlsolve = (NP_NL_SOLVER *) ReadArgvNumProc(mg,"S",NL_SOLVER_CLASS_NAME,argc,argv);
  sps->enlsolve = (NP_ENL_SOLVER *) ReadArgvNumProc(mg,"E",ENL_SOLVER_CLASS_NAME,argc,argv);
  if (sps->nlsolve==NULL && sps->enlsolve==NULL) return (NP_NOT_ACTIVE);
  sps->trans = (NP_TRANSFER *) ReadArgvNumProc(mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  if (sps->trans==NULL) return(NP_NOT_ACTIVE);
  sps->reinit = (NP_REINIT*) ReadArgvNumProc(mg,"R",REINIT_CLASS_NAME,argc,argv);
  if (sps->reinit==NULL) return(NP_NOT_ACTIVE);

  /* set configuration parameters */
  if (ReadArgvINT("baselevel",&(sps->baselevel),argc,argv)) sps->baselevel=0;
  if ((sps->baselevel<0)||(sps->baselevel>32)) return(NP_NOT_ACTIVE);
  if (ReadArgvINT("nested",&(sps->nested),argc,argv)) sps->nested=0;
  if ((sps->nested<0)||(sps->nested>1)) return(NP_NOT_ACTIVE);
  if (ReadArgvDOUBLE("r",&(sps->raster),argc,argv)) sps->raster=0.0;
  if (sps->raster>0.0)
    for (i=0; i<EXTENSION_MAX; i++)
      sps->raster_size[i]=sps->raster;
  sps->displayMode=ReadArgvDisplay(argc,argv);
  ret=NP_EXECUTABLE;

  /* read configuration for executability */
  tmp = ReadArgvVecDesc(base->mg,"sol",argc,argv); if (tmp==NULL) ret=NP_ACTIVE;
  if (AllocEVDForVD(mg,tmp,1,&(sps->pstep.sol_p0))) ret=NP_ACTIVE;
  if (sps->pstep.sol_p0==NULL) ret=NP_ACTIVE;
  if (ReadArgvChar("n0",sps->name[0],argc,argv)) ret=NP_ACTIVE;
  if (ReadArgvDOUBLE("p0",&(EVDD_E(sps->pstep.sol_p0,TOPLEVEL(mg),0)),argc,argv)) ret=NP_ACTIVE;
  sps->n=1;
  sps->scale=1.0;

  return (ret);
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

static INT SPS_Display (NP_BASE *theNumProc)
{
  INT i;
  NP_SPS *sps;

  /* get numprocs */
  sps = (NP_SPS *)theNumProc;

  UserWrite("\nSGS configuration:\n");
  if (sps->tass!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(sps->tass));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"A","---");
  if (sps->nlsolve!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(sps->nlsolve));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"S","---");
  if (sps->enlsolve!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"E",ENVITEM_NAME(sps->enlsolve));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"E","---");
  if (sps->trans!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(sps->trans));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (sps->reinit!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(sps->reinit));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"nested",(int)sps->nested);
  if (sps->displayMode==PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (sps->displayMode==PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (sps->displayMode==PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  if (sps->pstep.sol_p0)
  {
    UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_p0",ENVITEM_NAME(sps->pstep.sol_p0));
  }
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_p0","---");

  if (sps->pstep.sol_p1)
  {
    UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_p1",ENVITEM_NAME(sps->pstep.sol_p1));
  }
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_p1","---");

  for (i=0; i<sps->n; i++)
    UserWriteF(DISPLAY_NP_FORMAT_SF,sps->name[i],(float)EVDD_E(sps->pstep.sol_p0,TOPLEVEL(NP_MG(theNumProc)),i));

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  Execute                                                       */
/*                                                                          */
/* Purpose:   init solve cycle                                              */
/*                                                                          */
/* Input:     NumProcType Init Function                                     */
/*                                                                          */
/* Output: INT 0: ok                                                        */
/*             1: error                                                     */
/*                                                                          */
/****************************************************************************/

static INT PSTEP_Execute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_P_STEP *pstep;
  NP_SPS *sps;
  INT level,result,i;
  PSTEP_RESULT presult;
  MULTIGRID *mg;
  DOUBLE dp;

  /* get numprocs */
  pstep = (NP_P_STEP *)theNP;
  sps = (NP_SPS *)theNP;

  /* get level */
  mg = theNP->mg;
  level = CURRENTLEVEL(mg);

  /* call preprocess function */
  if (ReadArgvOption("pre",argc,argv))
    if (pstep->PreProcess!=NULL)
      if ((*pstep->PreProcess)(pstep,level,pstep->sol_p0,&result))
      {
        UserWriteF("PSTEP_Execute: PreProcess failed, error code %d\n",result);
        return (1);
      }

  /* execute */
  if (ReadArgvOption("step",argc,argv))
    if (pstep->Step!=NULL)
    {
      /* read step configuration */
      if (ReadArgvDOUBLE("dp",&dp,argc,argv)) dp=0.0;

      /* allocate tmp vector */
      if (AllocEVDFromEVD(mg,0,level,pstep->sol_p0,&(pstep->sol_p1))) return (1);
      EVDD_E(pstep->sol_p1,level,0)=EVDD_E(pstep->sol_p0,level,0)+dp;

      /* calculate solution */
      if ((*pstep->Step)(pstep,level,pstep->sol_p0,pstep->sol_p1,&presult))
      {
        UserWriteF("PSTEP_Execute: Step failed, error code\n");
        return (1);
      }
      if (!presult.converged)
      {
        UserWriteF("PSTEP_Execute: Step failed, cannot calculate solution\n");
        return (0);
      }
      else
      {
        dcopy(mg,0,level,ALL_VECTORS,pstep->sol_p0->vd,pstep->sol_p1->vd);
        for (i=0; i<pstep->sol_p0->n; i++)
          EVDD_E(pstep->sol_p0,level,i)=EVDD_E(pstep->sol_p1,level,i);
      }

      /* free tmp-vector */
      if (FreeEVD(mg,0,level,pstep->sol_p1)) return (1);
    }

  /* postprocess */
  if (ReadArgvOption("post",argc,argv))
    if (pstep->PostProcess!=NULL)
      if ((*pstep->PostProcess)(pstep,level,&result))
      {
        UserWriteF("PSTEP_Execute: PostProcess failed, error code %d\n", result);
        return (1);
      }

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  SPSConstruct													*/
/*																			*/
/* Purpose:	  create an object of class pstep								*/
/*																			*/
/* Input:	  theNP	num proc data structure									*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*			  else : error													*/
/*																			*/
/****************************************************************************/

static INT SPSConstruct (NP_BASE *theNP)
{
  NP_P_STEP *pstep;
  NP_NL_ASSEMBLE *na;
  NP_ENL_ASSEMBLE *ena;
  NP_SPS *sps;

  /* get numprocs */
  pstep = (NP_P_STEP *)theNP;
  sps = (NP_SPS *)theNP;
  na = &sps->nlass;
  ena = &sps->enlass;

  /* set base functions */
  theNP->Init    = SPS_Init;
  theNP->Display = SPS_Display;
  theNP->Execute = PSTEP_Execute;

  /* nonlinear assemble functions */
  na->PreProcess                      = SPS_NLPreProcess;
  na->PostProcess             = SPS_NLPostProcess;
  na->NLAssembleSolution      = SPS_NLAssembleSolution;
  na->NLAssembleDefect        = SPS_NLAssembleDefect;
  na->NLAssembleMatrix        = SPS_NLAssembleMatrix;

  /* nonlinear assemble functions */
  ena->PreProcess                     = SPS_ENLPreProcess;
  ena->PostProcess                    = SPS_ENLPostProcess;
  ena->ENLAssembleSolution    = SPS_ENLAssembleSolution;
  ena->ENLAssembleDefect              = SPS_ENLAssembleDefect;
  ena->ENLAssembleMatrix              = SPS_ENLAssembleMatrix;

  /* and the parameter solver */
  pstep->PreProcess       = SPS_PreProcess;
  pstep->Step                     = SPS_Step;
  pstep->PostProcess      = SPS_PostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitPStep - Init this file

   SYNOPSIS:
   INT InitTStep (void);

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

INT NS_DIM_PREFIX InitPStep (void)
{
  if (CreateClass (P_STEP_CLASS_NAME ".sps",sizeof(NP_SPS), SPSConstruct)) return (__LINE__);

  return (0);
}
