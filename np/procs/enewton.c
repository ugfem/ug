// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      enewton.c                                                     */
/*                                                                          */
/* Purpose:   newton method with nonlinear extension                        */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            IWR/TS                                                        */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheiner Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   August 14, 2000                                               */
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

#include "ugdevices.h"
#include "general.h"
#include "debug.h"
#include "gm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"
#include "block.h"

#include "nls.h"
#include "ls.h"
#include "els.h"
#include "assemble.h"
#include "transfer.h"
#include "enewton.h"

#ifdef __cplusplus
#ifdef __TWODIM__
using namespace UG2d;
#else
using namespace UG3d;
#endif
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define CSTART()    clock_start=CURRENT_TIME_LONG;
#define CSTOP(t,c)  t+=(CURRENT_TIME_LONG-clock_start);c++

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* variables for timing measurement		*/
static INT defect_c, newton_c, linear_c;
static DOUBLE defect_t, newton_t, linear_t;
static DOUBLE clock_start;

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

  NP_ENL_SOLVER nlsolver;               /* derived from abstract class NP_ENL_SOLVER		*/

  /* parameters to be set via npinit */
  NP_ELINEAR_SOLVER *esolve;            /* uses linear extended solver					*/
  NP_LINEAR_SOLVER *solve;              /* uses linear solver							*/
  NP_TRANSFER *trans;                           /* uses transgrid								*/
  INT displayMode;                              /* for PCR										*/
  INT maxit;                                            /* maximum number of newton iterations			*/
  INT linearRate;                               /* 1 if nonquadratic nonlin rate assumed                */
  DOUBLE lambda;                                /* nonlinear damp factor in $step and $nmg_step */
  EVEC_SCALAR linMinRed;            /* minimum reduction for linear solver		*/
  EVEC_SCALAR scale;                /* scaling of components						*/
  EVEC_SCALAR divFactor;            /* divergence factor for nonlin iteration	*/
  INT noLastDef;                                        /* no defect check after last iteration			*/
  INT force_iteration;                  /* if 1 at least 1 iteration is carried out     */

  /* and XDATA_DESCs */
  EMATDATA_DESC *J;                             /* the Matrix to be solved						*/
  EVECDATA_DESC *d;                             /* nonlinear defect								*/
  EVECDATA_DESC *dold;                  /* old nonlinear defect			                        */
  EVECDATA_DESC *dsave;                 /* last nonlinear defect			                        */
  EVECDATA_DESC *v;                             /* correction computed by newton step                   */
  VECDATA_DESC *s;                  /* help vector									*/

} NP_ENEWTON;   /* this is a final class */


/****************************************************************************/
/*D
   enewton - nonlinear solver numproc

   DESCRIPTION:
   This numproc executes a newton iteration based on the linearizazion
   provided by an assemble numproc. Depending of the linearized problem,
   the method can be applied as an inexact newton solver resp. a fixpoint
   iteration as well.
   This numproc can be configuered for several damping and line search
   strategies.
   D*/
/****************************************************************************/

static INT ENonLinearDefect (MULTIGRID *mg, INT level, INT init, EVECDATA_DESC *x, NP_ENEWTON *newton, NP_ENL_ASSEMBLE *ass, EVEC_SCALAR defect, INT *error)
{
  LRESULT lr;
  ELRESULT elr;
  INT i,n;

  n=VD_NCOMP(x->vd);

  /* project solution to all grid levels */
  if (newton->trans->PreProcessProject!=NULL)
    if ((*newton->trans->PreProcessProject)(newton->trans,0,level,error))                                               { *error = __LINE__; REP_ERR_RETURN(*error); }
  if ((*newton->trans->ProjectSolution)(newton->trans,0,level,x->vd,error))                                               { *error = __LINE__; REP_ERR_RETURN(*error); }
  if (newton->trans->PostProcessProject!=NULL)
    if ((*newton->trans->PostProcessProject)(newton->trans,0,level,error))                                     { *error = __LINE__; REP_ERR_RETURN(*error); }

  if (init)
  {
    /* preprocess assemble once before all calls */
    if (ass->PreProcess!=NULL)
      if ((*ass->PreProcess)(ass,0,level,x,error))                                                                                    { *error = __LINE__; REP_ERR_RETURN(*error); }

    /* set dirichlet conditions on all grid levels */
    if ((*ass->ENLAssembleSolution)(ass,0,level,x,error))                                                                           { *error = __LINE__; REP_ERR_RETURN(*error); }
  }

  /* compute new nonlinear defect */
  CSTART();
  for (i=0; i<newton->d->n; i++) EVDD_E(newton->d,level,i)=0.0;dset(mg,0,level,ALL_VECTORS,newton->d->vd,0.0);
  *error = 0;
  if ((*ass->ENLAssembleDefect)(ass,0,level,x,newton->d,newton->J,error))                                                 { *error = __LINE__; REP_ERR_RETURN(*error); }
  if (*error) return(0);
  CSTOP(defect_t,defect_c);

  if (UG_math_error)
  {
    UserWrite("math error in NLAssembleDefect\n");
    UG_math_error = 0;
    *error = __LINE__; REP_ERR_RETURN(*error);
  }

  /* compute norm of defect */
  if (newton->solve!=NULL)
  {
    if ((*newton->solve->Residuum)(newton->solve,0,level,newton->v->vd,newton->d->vd,newton->J->mm,&lr))            { *error = __LINE__; REP_ERR_RETURN(*error); }
    for (i=0; i<n; i++) defect[i] = lr.last_defect[i];
    for (i=0; i<x->n; i++) defect[i+n] = ABS(EVDD_E(newton->d,level,i));
  }
  else
  {
    if ((*newton->esolve->Residuum)(newton->esolve,0,level,newton->v,newton->d,newton->J,&elr))             { *error = __LINE__; REP_ERR_RETURN(*error); }
    for (i=0; i<n+x->n; i++) defect[i] = elr.last_defect[i];
  }

  return (0);
}

static INT ENewtonPreProcess  (NP_ENL_SOLVER *solve, INT level, EVECDATA_DESC *x, INT *result)
{
  NP_ENEWTON *newton;
  INT i;

  newton = (NP_ENEWTON *) solve;
  if (AllocEMDFromEVD(solve->base.mg,0,level,x,x,&newton->J)) NP_RETURN(1,result[0]);

  /* check function pointers in numprocs */
  if (newton->trans->base.status < NP_ACTIVE)
  {
    UserWrite("Newton: newton->trans not active\n");
    NP_RETURN(1,result[0]);
  }
  if (newton->trans->ProjectSolution==NULL)
  {
    UserWrite("Newton: newton->trans->ProjectSolution not defined\n");
    NP_RETURN(1,result[0]);
  }
  if (newton->solve!=NULL)
  {
    if (newton->solve->base.status < NP_ACTIVE)
    {
      UserWrite("Newton: newton->solve not active\n");
      NP_RETURN(1,result[0]);
    }
    if (newton->solve->Solver==NULL)
    {
      UserWrite("Newton: newton->solve->Solver not defined\n");
      NP_RETURN(1,result[0]);
    }
    if (newton->solve->Residuum==NULL)
    {
      UserWrite("Newton: newton->solve->Residuum not defined\n");
      NP_RETURN(1,result[0]);
    }
  }
  else
  {
    if (newton->esolve->base.status < NP_ACTIVE)
    {
      UserWrite("Newton: newton->esolve not active\n");
      NP_RETURN(1,result[0]);
    }
    if (newton->esolve->Solver==NULL)
    {
      UserWrite("Newton: newton->esolve->Solver not defined\n");
      NP_RETURN(1,result[0]);
    }
    if (newton->esolve->Residuum==NULL)
    {
      UserWrite("Newton: newton->esolve->Residuum not defined\n");
      NP_RETURN(1,result[0]);
    }
  }
  return(0);
}


static INT ENewtonPostProcess (NP_ENL_SOLVER *solve, INT level, EVECDATA_DESC *x, INT *result)
{
  NP_ENEWTON *newton;
  INT i;

  newton = (NP_ENEWTON *) solve;
  if (FreeEMD(solve->base.mg,0,level,newton->J)) REP_ERR_RETURN(1);

  return(0);
}

static INT ENewtonSolver (NP_ENL_SOLVER *nls, INT level, EVECDATA_DESC *x, NP_ENL_ASSEMBLE *ass, EVEC_SCALAR abslimit, EVEC_SCALAR reduction, ENLRESULT *res)
{
  NP_ENEWTON *newton;                                           /* object pointer						*/
  MULTIGRID *mg;                                                /* multigrid from base class			*/
  INT r;                                                                /* iteration counter			                */
  INT i,j,kk;                                                           /* some loop counters					*/
  char text[DISPLAY_WIDTH+4];                           /* display text in PCR					*/
  INT PrintID;                                                  /* print id for PCR						*/
  EVEC_SCALAR defect, defect2reach;             /* component--wise norm					*/
  EVEC_SCALAR defectmax;                                /* max defect without codivergence		*/
  INT n_unk;                                                    /* number of components in solution		*/
  EVEC_SCALAR linred;                                           /* parameters for linear solver			*/
  EVEC_SCALAR red_factor;                               /* convergence factor for linear iter	*/
  DOUBLE la;                                                            /* damping factor in line search		*/
  INT bl;                                                               /* baselevel returned by preprocess		*/
  INT error;                                                            /* for return value						*/
  LRESULT lr;                                                           /* result of linear solver				*/
  ELRESULT elr;                                                 /* result of ext. linear solver			*/
  DOUBLE s_tmp;                         /* tmp norm                             */
  INT use_second;
  DOUBLE eh[EXTENSION_MAX];                 /* extension help vector                */
  DOUBLE sc[EXTENSION_MAX*EXTENSION_MAX];       /* schur complement                 */
  DOUBLE si[EXTENSION_MAX*EXTENSION_MAX];       /* schur complement inverse         */

  /* get status */
  newton = (NP_ENEWTON *) nls;
  mg = nls->base.mg;
  UG_math_error = 0;

  /* fill result variable with error condition */
  res->error_code = 0;
  res->converged = 0;
  res->rho_first = 0.0;
  res->number_of_nonlinear_iterations = 0;
  res->number_of_line_searches = 0;
  res->total_linear_iterations = 0;
  res->max_linear_iterations = 0;
  res->exec_time = 0.0;

  /* initialize timers and counters */
  defect_c = newton_c = linear_c = 0;
  defect_t = newton_t = linear_t = 0.0;

  /* check function pointers in numprocs */
  if (ass->ENLAssembleSolution==NULL)
  {
    UserWrite("ENewton: ass->ENLAssembleSolution not defined\n");
    res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);
  }
  if (ass->ENLAssembleDefect==NULL)
  {
    UserWrite("ENewton: ass->ENLAssembleDefect not defined\n");
    res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);
  }
  if (ass->ENLAssembleMatrix==NULL)
  {
    UserWrite("ENewton: ass->ENLAssembleMatrix not defined\n");
    res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);
  }

  /* dynamic XDATA_DESC allocation */
  if (ass->A == NULL) ass->A = newton->J;
  if (AllocEVDFromEVD(mg,0,level,x,  &(newton->v)))                                               {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (AllocEVDFromEVD(mg,0,level,x,  &(newton->d)))                                               {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* get number of components */
  n_unk = VD_NCOMP(x->vd)+x->n;

  /* init ass once and compute nonlinear defect */
  if (ENonLinearDefect(mg,level,TRUE,x,newton,ass,defect,&error)!=0)              {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (error) goto exit;

  /* display norm of nonlinear defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(newton),'#',"\n");
  if (PrepareEPCR(newton->d,newton->displayMode,text,&PrintID))               {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (esc_mul(defect2reach,defect,reduction,newton->d))                                   {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (esc_mul(defectmax,defect,newton->divFactor,newton->d))                              {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  use_second=0; for (i=0; i<n_unk; i++) if (defectmax[i]==0.0) use_second=1;
  if (DoPCR(PrintID,defect,PCR_CRATE))                                                                    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  for (i=0; i<n_unk; i++) res->first_defect[i] = defect[i];

  /* check if iteration is necessary */
  if (esc_cmp(defect,abslimit,newton->d) && !newton->force_iteration)
  {
    res->converged = 1;
    for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
    res->error_code = 0;
    goto exit;
  }

  /* initialize reduction factor for linear solver */
  for (i=0; i<n_unk; i++) red_factor[i] = newton->linMinRed[i];

  /* do newton iterations */
  for (r=1; r<=newton->maxit; r++)
  {
    if (UG_math_error)
    {
      UserWrite("math error before newton loop !\n");
      UG_math_error = 0; res->error_code = __LINE__;
      break;
    }
    if (res->converged) break;

    /* compute jacobian */
    CSTART();
    for (i=0; i<newton->v->n; i++) EVDD_E(newton->v,level,i)=0.0;dset(mg,0,level,ALL_VECTORS,newton->v->vd,0.0);
    if ((*ass->ENLAssembleMatrix)(ass,0,level,x,newton->d,newton->v,newton->J,&error))
    {
      res->error_code = __LINE__;
      REP_ERR_RETURN(res->error_code);
    }
    CSTOP(newton_t,newton_c);
    if (UG_math_error)
    {
      UserWrite("math error in ENLAssembleMatrix !\n");
      UG_math_error = 0; res->error_code = __LINE__;
      break;
    }

    /* prepare and solve linear system */
    if (newton->solve!=NULL)
    {
      /* solve extended system using schur-complement */
      CSTART();
      for (i=0; i<n_unk; i++) linred[i] = red_factor[i];bl = 0;
      if (newton->solve->PreProcess!=NULL)
        if ((*newton->solve->PreProcess)(newton->solve,level,newton->v->vd,newton->d->vd,newton->J->mm,&bl,&error))
        {
          UserWriteF("ENewtonSolver: solve->PreProcess failed, error code %d\n",error); res->error_code = __LINE__;
                        #ifndef ModelP
          REP_ERR_RETURN (res->error_code);
                                        #endif
          REP_ERR_INC;
                                        #ifndef Debug
          return (res->error_code);
                                        #endif
        }

      /* build Schur complement */
      if (AllocVDFromVD(mg,0,level,x->vd,&(newton->s))) {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
      for (i=0; i<newton->d->n; i++)
      {
        dcopy(mg,0,level,ALL_VECTORS,newton->s,newton->J->me[i]);
        dset(mg,0,level,ALL_VECTORS,newton->v->vd,0.0);
        if ((*newton->solve->Residuum)(newton->solve,0,level,newton->v->vd,newton->s,newton->J->mm,&lr)) { res->error_code = __LINE__; goto exit; }
        if ((*newton->solve->Solver)(newton->solve,level,newton->v->vd,newton->s,newton->J->mm,newton->solve->abslimit,linred,&lr)) { res->error_code = __LINE__; goto exit; }
        res->total_linear_iterations += lr.number_of_linear_iterations;
        res->max_linear_iterations = MAX(res->max_linear_iterations,lr.number_of_linear_iterations);
        for (j=0; j<newton->d->n; j++)
        {
          if (ddot(mg,0,level,ON_SURFACE,newton->J->em[j],newton->v->vd,&(sc[i*newton->d->n+j])))
          {
            res->error_code = __LINE__;
            goto exit;
          }
          sc[i*newton->d->n+j]=EMDD_EE(newton->J,level,i*newton->d->n+j)-sc[i*newton->d->n+j];
        }
      }

      /* solve */
      dcopy(mg,0,level,ALL_VECTORS,newton->s,newton->d->vd);
      dset(mg,0,level,ALL_VECTORS,newton->v->vd,0.0);
      if ((*newton->solve->Residuum)(newton->solve,0,level,newton->v->vd,newton->s,newton->J->mm,&lr)) { res->error_code = __LINE__; goto exit; }
      if ((*newton->solve->Solver)(newton->solve,level,newton->v->vd,newton->s,newton->J->mm,newton->solve->abslimit,linred,&lr))
      {
        res->error_code = __LINE__;
        goto exit;
      }
      res->total_linear_iterations += lr.number_of_linear_iterations;
      res->max_linear_iterations = MAX(res->max_linear_iterations,lr.number_of_linear_iterations);
      for (i=0; i<newton->d->n; i++)
      {
        if (ddot(mg,0,level,ON_SURFACE,newton->J->em[i],newton->v->vd,&(eh[i])))
        {
          res->error_code = __LINE__;
          goto exit;
        }
        eh[i]=EVDD_E(newton->d,level,i)-eh[i];
      }
      if (SolveFullMatrix(newton->d->n,EVDD_E_PTR(newton->v,level),sc,eh)) { res->error_code = __LINE__; goto exit; }
      dcopy(mg,0,level,ALL_VECTORS,newton->s,newton->d->vd);
      for (i=0; i<newton->d->n; i++)
        if (daxpy(mg,0,level,ALL_VECTORS,newton->s,-1.0*EVDD_E(newton->v,level,i),newton->J->me[i])) { res->error_code = __LINE__; goto exit; }
      dset(mg,0,level,ALL_VECTORS,newton->v->vd,0.0);
      if ((*newton->solve->Residuum)(newton->solve,0,level,newton->v->vd,newton->s,newton->J->mm,&lr)) { res->error_code = __LINE__; goto exit; }
      if ((*newton->solve->Solver)(newton->solve,level,newton->v->vd,newton->s,newton->J->mm,newton->solve->abslimit,linred,&lr))
      {
        res->error_code = 0;
        goto exit;
      }
      res->total_linear_iterations += lr.number_of_linear_iterations;
      res->max_linear_iterations = MAX(res->max_linear_iterations,lr.number_of_linear_iterations);
      if (newton->solve->PostProcess!=NULL)
        if ((*newton->solve->PostProcess)(newton->solve,level,newton->v->vd,newton->d->vd,newton->J->mm,&error))
        {
          res->error_code = __LINE__;
          REP_ERR_RETURN(res->error_code);
        }
      CSTOP(linear_t,linear_c);
      if (UG_math_error) { UG_math_error = 0; res->error_code = __LINE__; }
      if (FreeVD(mg,0,level,newton->s)) REP_ERR_RETURN(1);

      /* if linear solver did not converge, return here */
      if (!lr.converged)
      {
        UserWrite("\nLinear solver did not converge in Newton method\n");
        if (newton->linearRate < 2)
        {
          res->error_code = 0;                           /* no error but exit */
          goto exit;
        }
      }
    }
    else if (newton->esolve!=NULL)
    {
      /* solve extended linear system with extended facilities */

      /* preprocess */
      CSTART();
      for (i=0; i<n_unk; i++) linred[i] = red_factor[i];bl = 0;
      if (newton->esolve->PreProcess!=NULL)
        if ((*newton->esolve->PreProcess)(newton->esolve,level,newton->v,newton->d,newton->J,&bl,&error))
        {
          UserWriteF("ENewtonSolver: esolve->PreProcess failed, error code %d\n",error); res->error_code = __LINE__;
                        #ifndef ModelP
          REP_ERR_RETURN (res->error_code);
                                        #endif
          REP_ERR_INC;
                                        #ifndef Debug
          return (res->error_code);
                                        #endif
        }

      /* solve */
      deset(mg,0,level,ALL_VECTORS,newton->v,0.0);
      if ((*newton->esolve->Residuum)(newton->esolve,0,level,newton->v,newton->d,newton->J,&elr)) { res->error_code = __LINE__; goto exit; }
      if ((*newton->esolve->Solver)(newton->esolve,level,newton->v,newton->d,newton->J,newton->esolve->abslimit,linred,&elr)) { res->error_code = __LINE__; goto exit; }
      res->total_linear_iterations += elr.number_of_linear_iterations;
      res->max_linear_iterations = MAX(res->max_linear_iterations,elr.number_of_linear_iterations);

      /* postprocess */
      if (newton->esolve->PostProcess!=NULL)
        if ((*newton->esolve->PostProcess)(newton->esolve,level,newton->v,newton->d,newton->J,&error))
        {
          res->error_code = __LINE__;
          REP_ERR_RETURN(res->error_code);
        }
      CSTOP(linear_t,linear_c);
      if (UG_math_error) { UG_math_error = 0; res->error_code = __LINE__; }

      /* if linear solver did not converge, return here */
      if (!elr.converged)
      {
        UserWrite("\nLinear solver did not converge in Newton method\n");
        if (newton->linearRate < 2)
        {
          res->error_code = 0;                           /* no error but exit */
          goto exit;
        }
      }
    }
    else
      assert(0);

    /* update nonlinear solution */
    daxpy(mg,0,level,ALL_VECTORS,x->vd,-1.0,newton->v->vd);
    for (i=0; i<x->n; i++)
      EVDD_E(x,level,i)-=EVDD_E(newton->v,level,i);

    /* print norm of defect */
    if (ENonLinearDefect(mg,level,FALSE,x,newton,ass,defect,&error))
      if (error) goto exit;
    if (DoPCR(PrintID,defect,PCR_CRATE)) {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

    /* check convergence of nonlinear iteration */
    if (esc_cmp(defect,abslimit,newton->d)) {res->converged=1; break;}
    if (esc_cmp(defect,defect2reach,newton->d)) {res->converged=1; break;}
    if (use_second)
    {
      use_second=0;
      if (esc_mul(defectmax,defect,newton->divFactor,newton->d))              {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
    }
    else
    {
      if (!esc_cmp(defect,defectmax,newton->d)) break;
    }

    /* compute new reduction factor, assuming quadratic convergence */
    for (i=0; i<n_unk; i++)
    {
      red_factor[i] = MIN(reduction[i]*reduction[i],newton->linMinRed[i]);
      if (newton->linearRate == 1)
        red_factor[i] = MIN(reduction[i],newton->linMinRed[i]);
      if (newton->linearRate == 2)
        red_factor[i] = newton->linMinRed[i];
    }
  }

  /* print norm of defect */
  if (DoPCR(PrintID,defect,PCR_AVERAGE))  {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* report results and mean execution times */
  res->error_code = 0;
  res->number_of_nonlinear_iterations = newton_c;
  res->number_of_line_searches = defect_c;
  for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
  if (!res->converged) UserWriteF("ENL SOLVER: desired convergence not reached\n");
  UserWriteF("AVG EXEC TIMES: DEF[%2d]=%10.4g JAC[%2d]=%10.4g LIN[%2d]=%10.4g\n",
             defect_c,defect_t/defect_c,newton_c,newton_t/newton_c,linear_c,linear_t/linear_c);
  res->exec_time = defect_t+newton_t+linear_t;

  /* postprocess assemble once at the end */
  if (ass->PostProcess!=NULL)
    if ((*ass->PostProcess)(ass,0,level,x,newton->d,newton->J,&error))
    {
      res->error_code = __LINE__;
      REP_ERR_RETURN(res->error_code);
    }

exit:
  if (PostPCR(PrintID,NULL))                              {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* deallocate local XDATA_DESCs */
  if (FreeEVD(mg,0,level,newton->d)) REP_ERR_RETURN(1);
  if (FreeEVD(mg,0,level,newton->v)) REP_ERR_RETURN(1);
  if (res->error_code==0)
    return (0);
  else
    REP_ERR_RETURN (res->error_code);
}

static INT ENewtonInit (NP_BASE *base, INT argc, char **argv)
{
  NP_ENEWTON *newton;
  INT i;

  newton = (NP_ENEWTON *) base;

  /* read  data descs */
  newton->J = ReadArgvEMatDesc(base->mg,"J",argc,argv);
  newton->d = ReadArgvEVecDesc(base->mg,"d",argc,argv);
  newton->v = ReadArgvEVecDesc(base->mg,"v",argc,argv);
  newton->dold = ReadArgvEVecDesc(base->mg,"dold",argc,argv);
  newton->dsave = ReadArgvEVecDesc(base->mg,"dsave",argc,argv);

  /* read other numprocs */
  newton->trans = (NP_TRANSFER *) ReadArgvNumProc(base->mg,"T",TRANSFER_CLASS_NAME,
                                                  argc,argv);
  if (newton->trans == NULL)
  {
    PrintErrorMessage('E',"ENewtonInit","cannot read transfer num proc");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  newton->solve = (NP_LINEAR_SOLVER *) ReadArgvNumProc(base->mg,"S",LINEAR_SOLVER_CLASS_NAME,argc,argv);
  if (newton->solve == NULL)
  {
    newton->esolve = (NP_ELINEAR_SOLVER *) ReadArgvNumProc(base->mg,"S",ELINEAR_SOLVER_CLASS_NAME,argc,argv);
    if (newton->esolve == NULL)
    {
      PrintErrorMessage('E',"ENewtonInit","cannot read neither solve nor esolve num proc");
      REP_ERR_RETURN(NP_NOT_ACTIVE);
    }
  }
  if (ReadArgvINT("fi",&(newton->force_iteration),argc,argv)) newton->force_iteration = 0;
  if (ReadArgvINT("maxit",&(newton->maxit),argc,argv)) newton->maxit = 50;
  if ((newton->maxit<0)||(newton->maxit>1000))
  {
    PrintErrorMessage('E',"ENewtonInit","maxit <= 1000");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("linrate",&(newton->linearRate),argc,argv)) newton->linearRate = 0;
  if ((newton->linearRate<0)||(newton->linearRate>2))
  {
    PrintErrorMessage('E',"ENewtonInit","linrate = 0,1 or 2");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (esc_read(newton->linMinRed,NP_FMT(newton),newton->d,"linminred",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      newton->linMinRed[i] = 1e-4;

  for (i=0; i<MAX_VEC_COMP; i++)
    if ((newton->linMinRed[i]<0.0)||(newton->linMinRed[i]>=1.0))
    {
      PrintErrorMessage('E',"ENewtonInit","linminred must be in (0,1)");
      REP_ERR_RETURN(NP_NOT_ACTIVE);
    }

  if (esc_read(newton->scale,NP_FMT(newton),newton->d,"scale",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      newton->scale[i] = 1.0;

  if (esc_read(newton->divFactor,NP_FMT(newton),newton->d,"divfac",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      newton->divFactor[i] = 1e5;

  for (i=0; i<MAX_VEC_COMP; i++)
    if ((newton->divFactor[i]<=1.0))
    {
      PrintErrorMessage('E',"ENewtonInit","divfac must be in )1,inf(");
      REP_ERR_RETURN(NP_NOT_ACTIVE);
    }

  /* set noLastDef option */
  newton->noLastDef = ReadArgvOption("noLastDef",argc,argv);

  /* set display option */
  newton->displayMode = ReadArgvDisplay(argc,argv);

  /* call general nls init */
  return (NPENLSolverInit(&(newton->nlsolver),argc,argv));
}

static INT ENewtonDisplay (NP_BASE *theNumProc)
{
  NP_ENEWTON *newton;

  newton  = (NP_ENEWTON*) theNumProc;

  /* general nls display */
  NPENLSolverDisplay(&(newton->nlsolver));

  if (newton->J != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"J",ENVITEM_NAME(newton->J));
  if (newton->v != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"v",ENVITEM_NAME(newton->v));
  if (newton->d != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"d",ENVITEM_NAME(newton->d));
  if (newton->dold != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"dold",ENVITEM_NAME(newton->dold));
  if (newton->dsave != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"dsave",ENVITEM_NAME(newton->dsave));

  if (newton->solve != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(newton->solve));
  else if (newton->esolve != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(newton->esolve));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S","---");
  if (newton->trans != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(newton->trans));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (newton->displayMode == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (newton->displayMode == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (newton->displayMode == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  UserWriteF(DISPLAY_NP_FORMAT_SI,"maxit",(int)newton->maxit);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"linrate",(int)newton->linearRate);
  if (newton->d!=NULL) if (esc_disp(newton->linMinRed,newton->d,"linMinRed")) REP_ERR_RETURN (1);
  if (newton->d!=NULL) if (esc_disp(newton->divFactor,newton->d,"divfac")) REP_ERR_RETURN (1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"lambda",(float)newton->lambda);

  return (0);
}

static INT ENewtonConstruct (NP_BASE *theNP)
{
  NP_ENL_SOLVER *np;

  /* set general functions */
  theNP->Init = ENewtonInit;
  theNP->Display = ENewtonDisplay;
  theNP->Execute = NULL;

  np = (NP_ENL_SOLVER *) theNP;
  np->PreProcess = ENewtonPreProcess;
  np->Solver = ENewtonSolver;
  np->PostProcess = ENewtonPostProcess;

  return(0);
}


/****************************************************************************/
/*
   InitENewtonSolver  - Init this file

   SYNOPSIS:
   INT InitENewtonSolver (void);

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

INT InitENewtonSolver (void)
{
  if (CreateClass (ENL_SOLVER_CLASS_NAME ".enewton",sizeof(NP_ENEWTON), ENewtonConstruct))
    REP_ERR_RETURN (__LINE__);

  return (0);
}
