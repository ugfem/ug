// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      newton.c                                                      */
/*                                                                          */
/* Purpose:   newton method (type nls)                                      */
/*                                                                          */
/* Author:    Peter Bastian                                                                                             */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   November 29, 1996                                             */
/*            January  07, 1997 it works !                                  */
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

#include "nls.h"
#include "ls.h"
#include "assemble.h"
#include "transfer.h"
#include "newton.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAX_LINE_SEARCH                 20
#define LINE_SEARCH_REDUCTION   0.5

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
  NP_NL_SOLVER nlsolver;                /* derived from abstract class NP_NL_SOLVER		*/

  /* parameters to be set via npinit */
  NP_LINEAR_SOLVER *solve;              /* uses linear solver							*/
  NP_TRANSFER *trans;                           /* uses transgrid								*/
  INT displayMode;                              /* for PCR										*/
  INT maxit;                                            /* maximum number of newton iterations			*/
  INT linearRate;                               /* 1 if nonquadratic nonlin rate assumed                */
  INT lineSearch;                               /* do line search                                                               */
  INT maxLineSearch;                            /* maximum number of line search steps			*/
  DOUBLE rhoReass;                              /* reassemble if nonlin conv worth than this    */
  DOUBLE lambda;                                /* nonlinear damp factor in $step and $nmg_step */
  DOUBLE linMinRed[MAX_VEC_COMP];      /* minimum reduction for linear solver		*/
  DOUBLE scale[MAX_VEC_COMP];       /* scaling of components						*/
  DOUBLE divFactor[MAX_VEC_COMP];         /* divergence factor for nonlin iteration	*/
  INT noLastDef;                                        /* no defect check after last iteration			*/
  INT force_iteration;                  /* if 1 at least 1 iteration is carried out     */

  /* and XDATA_DESCs */
  MATDATA_DESC *J;                              /* the Matrix to be solved						*/
  VECDATA_DESC *d;                              /* nonlinear defect								*/
  VECDATA_DESC *dold;                           /* old nonlinear defect			                        */
  VECDATA_DESC *dsave;                  /* last nonlinear defect			                        */
  VECDATA_DESC *v;                              /* correction computed by newton step                   */
  VECDATA_DESC *s;                              /* saved nonlinear solution						*/

} NP_NEWTON;    /* this is a final class */


/****************************************************************************/
/*D
   newton - nonlinear solver numproc

   DESCRIPTION:
   This numproc executes a newton iteration based on the linearizazion
   provided by an assemble numproc. Depending of the linearized problem,
   the method can be applied as an inexact newton solver resp. a fixpoint
   iteration as well.
   This numproc can be configuered for several damping and line search
   strategies.

   .vb
   npinit <name> $x <sol>
              $A <assemble> $T <transfer> $S <solver>
              $abslimit <sc double list> $red <sc double list>
              $rhoreass <double> $lsteps <l> $maxit <m> $line {0|1}
                          $linrate {0|1|2} [$lambda {double}]
                          $linminred <sc double list> [$display {no|red|full}] [$J <mat>];
   .ve

   .  $x~<sol>           - solution vector
   .  $A~<assemble>      - assemble numproc of type 'NP_NL_ASSEMBLE'
   .  $S~<linear solver> - linear solver numproc of type 'NP_LINEAR_SOLVER'
   .  $T~<transfer>      - transfer numproc
   .  $abslimit~<sc~double~list>  - absolute limit for the defect (default 1.0E-10)
   .  $red~<sc~double~list> - reduction factor
   .  $rhoreass~<double>    - reassemble if nonlinear convergence rate worse than this
   .  $lsteps~<l>           - maximum number of line search steps
   .  $maxit~<m>            - maximum number of nonlinear iterations
   .  $line~{0|1}           - do line search
   .  $linrate~{0|1|2}      - 0: quadratic nonlin rate assumed, 1: nonquadratic nonlin rate assumed, 2: accept the corretion always
   .  $lambda~{double}      - nonlinear damp factor
   .  $linminred~<sc~double~list> - minimum reduction for linear solver
   .  $display~{no|red|full}] - display mode
   .  $J~<mat>                - Jacobi matrix (optional)

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n   nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata
   .n   if only a single value is specified, this will be used for all components

   'npexecute <name> [$i] [$s] [$p];

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
 # grid transfer numproc
   npcreate transfer $c transfer;
   npinit transfer;

 # assemble numproc
   npcreate nlass $c boxnl;
   npinit nlass;

 # linear solver and iteration numprocs
   npcreate smooth $c ilu;
   npinit smooth $damp 0.92;

   npcreate basesolver $c ls;
   npinit basesolver $red 0.0001 $m 50 $I smooth;

   npcreate lmgc $c lmgc;
   npinit lmgc $S smooth smooth basesolver $T transfer $n1 2 $n2 2 $g 1 $b 0;

   npcreate mgs $c ls;
   npinit mgs $m 40 $I lmgc $display full;

 # nonlinear solver numproc
   npcreate newton $c newton;
   npinit newton $x sol
              $A nlass $T transfer $S mgs
              $abslimit 1.0E-10 $red 1.0E-5
              $rhoreass 0.8 $lsteps 6 $maxit 50 $line 1 $linrate 0
              $lambda 1.0 $linminred 1.0E-2 $display full;

   npexecute newton $i $s $p;
   .ve
   D*/
/****************************************************************************/

static INT NonLinearDefect (MULTIGRID *mg, INT level, INT init, VECDATA_DESC *x,
                            NP_NEWTON *newton, NP_NL_ASSEMBLE *ass, VEC_SCALAR defect, INT *error)
{
  LRESULT lr;                           /* result of linear solver				*/
  INT i,n_unk;

  n_unk = VD_NCOMP(x);

  /* project solution to all grid levels */
  if (newton->trans->PreProcessProject!=NULL)
    if ((*newton->trans->PreProcessProject)
          (newton->trans,0,level,error)) {
      *error = __LINE__;
      REP_ERR_RETURN(*error);
    }
  if ((*newton->trans->ProjectSolution)(newton->trans,0,level,x,error)) {
    *error = __LINE__;
    REP_ERR_RETURN(*error);
  }
  if (newton->trans->PostProcessProject!=NULL)
    if ((*newton->trans->PostProcessProject)
          (newton->trans,0,level,error)) {
      *error = __LINE__;
      REP_ERR_RETURN(*error);
    }

  if (init)
  {
    /* preprocess assemble once before all calls */
    if (ass->PreProcess!=NULL)
      if ((*ass->PreProcess)(ass,0,level,x,error)) {
        *error = __LINE__;
        REP_ERR_RETURN(*error);
      }

    /* set dirichlet conditions on all grid levels */
    if ((*ass->NLAssembleSolution)(ass,0,level,x,error)) {
      *error = __LINE__;
      REP_ERR_RETURN(*error);
    }
  }

  /* compute new nonlinear defect */
  CSTART();
  dset(mg,0,level,ALL_VECTORS,newton->d,0.0);
  *error = 0;
  if ((*ass->NLAssembleDefect)(ass,0,level,x,newton->d,newton->J,error)) {
    *error = __LINE__;
    REP_ERR_RETURN(*error);
  }
  if (*error)
    return(0);
  CSTOP(defect_t,defect_c);

  if (newton->lineSearch == 3)
    dcopy(mg,0,level,ALL_VECTORS,newton->dsave,newton->d);
  if (UG_math_error) {
    UserWrite("math error in NLAssembleDefect\n");
    UG_math_error = 0;
    *error = __LINE__;
    REP_ERR_RETURN(*error);
  }

  IFDEBUG(np,3)
  UserWrite("---- After computation of nonlinear defect\n");
  ListVectorRange(mg,0,level,0,0,1000,FALSE,TRUE,~(INT)1,LV_MOD_DEFAULT);
  ENDDEBUG

  /* compute norm of defect */
  if ((*newton->solve->Residuum)(newton->solve,0,level,newton->v,newton->d,newton->J,&lr)) {
    *error = __LINE__;
    REP_ERR_RETURN(*error);
  }
  for (i=0; i<n_unk; i++) defect[i] = lr.last_defect[i];

  return (0);
}

static INT NewtonPreProcess  (NP_NL_SOLVER *solve, INT level, VECDATA_DESC *x, INT *result)
{
  NP_NEWTON *newton;

  newton = (NP_NEWTON *) solve;
  if (AllocMDFromVD(solve->base.mg,0,level,x,x,&newton->J))
    NP_RETURN(1,result[0]);

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
  return(0);
}


static INT NewtonPostProcess (NP_NL_SOLVER *solve, INT level, VECDATA_DESC *x, INT *result)
{
  NP_NEWTON *newton;

  newton = (NP_NEWTON *) solve;
  if (FreeMD(solve->base.mg,0,level,newton->J)) REP_ERR_RETURN(1);

  return(0);
}



static INT NewtonSolver      (NP_NL_SOLVER *nls, INT level, VECDATA_DESC *x,
                              NP_NL_ASSEMBLE *ass, VEC_SCALAR abslimit, VEC_SCALAR reduction,
                              NLRESULT *res)
{
  NP_NEWTON *newton;                                            /* object pointer						*/
  MULTIGRID *mg;                                                /* multigrid from base class			*/
  INT r;                                                                /* iteration counter			                */
  INT i,kk;                                                             /* some loop counters					*/
  char text[DISPLAY_WIDTH+4];                           /* display text in PCR					*/
  INT PrintID;                                                  /* print id for PCR						*/
  VEC_SCALAR defect, defect2reach;              /* component--wise norm					*/
  VEC_SCALAR defectmax;                                 /* max defect without codivergence		*/
  INT n_unk;                                                            /* number of components in solution		*/
  DOUBLE s,sold,sprime,s2reach,sred;            /* combined defect norm					*/
  INT reassemble=1;                                             /* adaptive computation of jacobian		*/
  VEC_SCALAR linred;                                            /* parameters for linear solver			*/
  DOUBLE red_factor[MAX_VEC_COMP];              /* convergence factor for linear iter	*/
  DOUBLE la;                                                            /* damping factor in line search		*/
  DOUBLE rho[MAX_LINE_SEARCH+1];                /* reduction factors of linesearch		*/
  DOUBLE rhomin;                                                /* best reduction if !accept			*/
  INT best_ls;                                                  /* best ls if !accept					*/
  INT accept;                                                           /* line search accepted					*/
  INT bl;                                                               /* baselevel returned by preprocess		*/
  INT error;                                                            /* for return value						*/
  LRESULT lr;                                                           /* result of linear solver				*/
  DOUBLE lambda_old;                    /* last accepted lamda                  */
  DOUBLE lambda_min;                    /* minimal lamda                        */
  DOUBLE mu;                            /* guess for lambda                     */
  DOUBLE s_tmp;                         /* tmp norm                             */
  INT use_second;

  /* get status */
  newton = (NP_NEWTON *) nls;      /* cast from abstract base class to final class*/
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
  if (newton->lineSearch == 3) {
    lambda_min = lambda_old = 1.0;
    for (kk=1; kk<=newton->maxLineSearch; kk++)
      lambda_min = LINE_SEARCH_REDUCTION * lambda_min;
  }
  /* check function pointers in numprocs */
  if (ass->NLAssembleSolution==NULL)
  {
    UserWrite("Newton: ass->NLAssembleSolution not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  if (ass->NLAssembleDefect==NULL)
  {
    UserWrite("Newton: ass->NLAssembleDefect not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  if (ass->NLAssembleMatrix==NULL)
  {
    UserWrite("Newton: ass->NLAssembleMatrix not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  /* dynamic XDATA_DESC allocation */
  if (ass->A == NULL)
    ass->A = newton->J;
  if (AllocVDFromVD(mg,0,level,x,  &(newton->v)))
  {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(newton->d)))
  {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (newton->lineSearch == 3) {
    if (AllocVDFromVD(mg,0,level,x,  &(newton->dold)))
    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
    if (AllocVDFromVD(mg,0,level,x,  &(newton->dsave)))
    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  }
  /* get number of components */
  n_unk = VD_NCOMP(x);

  /* init ass once and compute nonlinear defect */
  if (NonLinearDefect(mg,level,TRUE,x,newton,ass,defect,&error)!=0)
  {
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  if (error)
    goto exit;
  /* display norm of nonlinear defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(newton),'#',"\n");
  if (PreparePCR(newton->d,newton->displayMode,text,&PrintID))    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (sc_mul(defect2reach,defect,reduction,newton->d))                    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (sc_mul(defectmax,defect,newton->divFactor,newton->d))               {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  use_second=0; for (i=0; i<n_unk; i++) if (defectmax[i]==0.0) use_second=1;
  if (DoPCR(PrintID,defect,PCR_CRATE))                                                    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  for (i=0; i<n_unk; i++) res->first_defect[i] = defect[i];

  /* compute single norm */
  s = 0.0;
  for (i=0; i<n_unk; i++)
    s += newton->scale[i]*newton->scale[i]*defect[i]*defect[i];
  s = sqrt(s);
  sprime = s;
  sold = s * sqrt(2.0);
  sred = 1.0E10; for (i=0; i<n_unk; i++) sred = MIN(sred,reduction[i]);
  s2reach = s*sred;
  if (newton->lineSearch)
    if (newton->displayMode == PCR_FULL_DISPLAY)
      UserWriteF(" ++ s=%12.4E Initial nonlinear residual\n",s);

  /* check if iteration is necessary */
  if (sc_cmp(defect,abslimit,newton->d) && !newton->force_iteration) {
    res->converged = 1;
    for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
    res->error_code = 0;
    goto exit;
  }

  /* initialize reduction factor for linear solver */
  for (i=0; i<n_unk; i++) red_factor[i] = newton->linMinRed[i];
  reassemble = 1;

  /* do newton iterations */
  for (r=1; r<=newton->maxit; r++)
  {
    if (UG_math_error) {
      UserWrite("math error before newton loop !\n");
      UG_math_error = 0;
      res->error_code = __LINE__;
      break;
    }

    if (res->converged) break;             /* solution already found */

    /* compute jacobian */
    CSTART();
    dset(mg,0,level,ALL_VECTORS,newton->v,0.0);
    if (reassemble)
    {
      if ((*ass->NLAssembleMatrix)(ass,0,level,x,newton->d,newton->v,newton->J,&error)) {
        res->error_code = __LINE__;
        REP_ERR_RETURN(res->error_code);
      }
      reassemble = 0;
    }
    CSTOP(newton_t,newton_c);
    if (UG_math_error) {
      UserWrite("math error in NLAssembleMatrix !\n");
      UG_math_error = 0;
      res->error_code = __LINE__;
      break;
    }

    /* solve linear system */
    CSTART();
    for (i=0; i<n_unk; i++) linred[i] = red_factor[i];bl = 0;
    if (newton->solve->PreProcess!=NULL)
      if ((*newton->solve->PreProcess)(newton->solve,level,newton->v,newton->d,newton->J,&bl,&error)) {
        UserWriteF("NewtonSolver: solve->PreProcess failed, error code %d\n",error);
        res->error_code = __LINE__;
                #ifndef ModelP
        REP_ERR_RETURN (res->error_code);
                #endif
        REP_ERR_INC;
                                #ifndef Debug
        return (res->error_code);
                #endif
      }
    if ((*newton->solve->Residuum)(newton->solve,0,level,newton->v,newton->d,newton->J,&lr))
    {
      res->error_code = __LINE__;
      goto exit;
    }
    if ((*newton->solve->Solver)(newton->solve,level,newton->v,newton->d,newton->J,newton->solve->abslimit,linred,&lr))
    {
      res->error_code = __LINE__;
      goto exit;
    }
    if (newton->solve->PostProcess!=NULL)
      if ((*newton->solve->PostProcess)(newton->solve,level,newton->v,newton->d,newton->J,&error)) {
        res->error_code = __LINE__;
        REP_ERR_RETURN(res->error_code);
      }
    CSTOP(linear_t,linear_c);
    if (UG_math_error) {
      UG_math_error = 0;
      res->error_code = __LINE__;
    }

    /* if linear solver did not converge, return here */
    if (!lr.converged) {
      UserWrite("\nLinear solver did not converge in Newton method\n");
      if (newton->linearRate < 2) {
        res->error_code = 0;                     /* no error but exit */
        goto exit;                         /* or goto exit2 ??? */
      }
    }

    IFDEBUG(np,3)
    UserWrite("---- After linear solver\n");
    ListVectorRange(mg,0,level,0,0,1000,FALSE,TRUE,~(INT)1,LV_MOD_DEFAULT);
    ENDDEBUG

    /* linear solver statistics */
    res->total_linear_iterations += lr.number_of_linear_iterations;
    res->max_linear_iterations = MAX(res->max_linear_iterations,lr.number_of_linear_iterations);

    if (newton->lineSearch)
      if (newton->displayMode == PCR_FULL_DISPLAY)
        UserWriteF(" ++ newton step %3d\n",r);

    /* save current solution for line search */
    if (AllocVDFromVD(mg,0,level,x,  &(newton->s)))
    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
    dcopy(mg,0,level,ALL_VECTORS,newton->s,x);

    /* do a line search */
    la = newton->lambda; accept=0;
    for (kk=1; kk<=newton->maxLineSearch; kk++) {

      if (newton->lineSearch == 3) {
        if (lambda_old == 1.0) {
          if (kk == 1)
            mu = sold*sold*0.5 / (s*s);
          else
            mu = sold * 0.5 / sprime;
          UserWriteF("lambda_old %f s %f sold %f sprime %f mu %f\n",
                     lambda_old,s,sold,sprime,mu);
          la = MIN(1.0,mu);
        }
        else {


          if (dnrm2(mg,0,level,ON_SURFACE,newton->dold,&s_tmp)!=NUM_OK)
            NP_RETURN(1,res->error_code);
          UserWriteF("kk %d s_old %f\n",kk,s_tmp);

          if (dnrm2(mg,0,level,ON_SURFACE,newton->dsave,&s_tmp)!=NUM_OK)
            NP_RETURN(1,res->error_code);
          UserWriteF("kk %d s_save %f\n",kk,s_tmp);

          if (dscal(mg,0,level,ALL_VECTORS,newton->dold,lambda_old - 1.0))
            REP_ERR_RETURN (1);
          if (dadd(mg,0,level,ALL_VECTORS,newton->dold,newton->dsave) != NUM_OK)
            REP_ERR_RETURN (1);
          if (dnrm2(mg,0,level,ON_SURFACE,newton->dold,&s_tmp)!=NUM_OK)
            NP_RETURN(1,res->error_code);
          s_tmp *= s_tmp;
          if (dsub(mg,0,level,ALL_VECTORS,newton->dold,newton->dsave) != NUM_OK)
            REP_ERR_RETURN (1);
          if (dscal(mg,0,level,ALL_VECTORS,newton->dold,1.0 / (lambda_old - 1.0)))
            REP_ERR_RETURN (1);
          if (s_tmp <= 0.0)
            la = lambda_old * LINE_SEARCH_REDUCTION;
          else {
            s_tmp = sqrt(s_tmp);
            UserWriteF("kk %d s_tmp %f %f\n",kk,s_tmp,
                       sold*(lambda_old-1)+s);
            if (kk == 1)
              mu=lambda_old*lambda_old*sold*sold*0.5/(s*s_tmp);
            else
              mu = lambda_old*lambda_old*s*0.5 / s_tmp;
            UserWriteF("lambda_old %f s %f sold %f sprime %f mu %f\n",
                       lambda_old,s,sold,sprime,mu);
            la = MIN(1.0,mu);
          }

          mu = sold * 0.5 *lambda_old / sprime;

        }
        if (kk == 1)
          if (newton->lineSearch == 3)
            dcopy(mg,0,level,ALL_VECTORS,newton->dold,newton->dsave);
        la = MAX(la,lambda_min);
        lambda_old = la;
      }

      /* update solution */
      dcopy(mg,0,level,ALL_VECTORS,x,newton->s);
      daxpy(mg,0,level,ALL_VECTORS,x,-la,newton->v);

      if (newton->maxit==r && newton->noLastDef) {
        res->error_code = 0;
        if (FreeVD(mg,0,level,newton->s)) REP_ERR_RETURN(1);
        goto exit2;
      }
      if (NonLinearDefect(mg,level,FALSE,x,newton,ass,defect,&error)!=0)
      {
        res->error_code = __LINE__;
        REP_ERR_RETURN(res->error_code);
      }
      if (error)
        goto exit;

      /* compute single norm */
      sold = sprime;
      sprime = 0.0;
      for (i=0; i<n_unk; i++)
        sprime += newton->scale[i]*newton->scale[i]*defect[i]*defect[i];
      sprime = sqrt(sprime);

      rho[kk] = sprime/s;

      /* print results */
      if (newton->lineSearch)
        if (newton->displayMode == PCR_FULL_DISPLAY)
          UserWriteF(" ++ ls=%2d, s=%12.4E, rho=%8.4g, lambda= %8.4g\n",
                     kk,sprime,rho[kk],fabs(la));

      if (sprime/s<=1-0.25*fabs(la) || !newton->lineSearch) {
        lambda_old = la;
        accept=1;
        break;
      }

      /* else reduce lambda */
      la = LINE_SEARCH_REDUCTION*la;
    }

    /* if not accepted */
    if (!accept)
      switch (newton->lineSearch)
      {
      case 1 :
        /* break iteration */
        UserWrite("line search not accepted, Newton diverged\n");
        res->error_code = 0;
        if (FreeVD(mg,0,level,newton->s)) REP_ERR_RETURN(1);
        goto exit;

      case 2 :
        /* accept best result */
        best_ls = 1;
        rhomin = rho[best_ls];
        for (kk=2; kk<=newton->maxLineSearch; kk++)
          if (rhomin>rho[kk]) {
            rhomin = rho[kk];
            best_ls = kk;
          }

        UserWriteF(" ++ accepting linesearch %d\n",best_ls);

        /* set lambda factor */
        la = newton->lambda * pow(LINE_SEARCH_REDUCTION,best_ls-1);
        lambda_old = la;

        /* update solution */
        dcopy(mg,0,level,ALL_VECTORS,x,newton->s);
        daxpy(mg,0,level,ALL_VECTORS,x,-la,newton->v);

        if (NonLinearDefect(mg,level,FALSE,x,newton,ass,
                            defect,&error)!=0)
        {
          res->error_code = __LINE__;
          REP_ERR_RETURN(res->error_code);
        }
        if (error)
          goto exit;
        break;

        /*default: accept */
      }
    if (FreeVD(mg,0,level,newton->s)) REP_ERR_RETURN(1);

    /* print norm of defect */
    if (DoPCR(PrintID,defect,PCR_CRATE)) {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

    /* save convergence of first step, may be used to increase stepsize */
    if (r==1) res->rho_first = sprime/s;

    /* reassemble if nonlinear convergence bad */
    if (sprime/s >= newton->rhoReass) reassemble = 1;
    if (!newton->linearRate) reassemble = 1;

    /* check convergence of nonlinear iteration */
    if (sprime<s2reach) {res->converged=1; break;}
    if (sc_cmp(defect,abslimit,newton->d)) {res->converged=1; break;}
    if (sc_cmp(defect,defect2reach,newton->d)) {res->converged=1; break;}
    if (use_second)
    {
      use_second=0;
      if (sc_mul(defectmax,defect,newton->divFactor,newton->d))               {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
    }
    else
    {
      if (!sc_cmp(defect,defectmax,newton->d)) break;
    }

    /* compute new reduction factor, assuming quadratic convergence */
    for (i=0; i<n_unk; i++)
    {
      red_factor[i] = MIN((sprime/s)*(sprime/s),newton->linMinRed[i]);
      if (newton->linearRate == 1)
        red_factor[i] = MIN(sprime/s,newton->linMinRed[i]);
      if (newton->linearRate == 2)
        red_factor[i] = newton->linMinRed[i];
    }

    /* accept new iterate */
    sold = s;
    s = sprime;
  }

  /* print norm of defect */
  if (DoPCR(PrintID,defect,PCR_AVERAGE))  {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* report results and mean execution times */
  res->error_code = 0;
  res->number_of_nonlinear_iterations = newton_c;
  res->number_of_line_searches = defect_c;
  for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
  if (!res->converged) UserWriteF("NL SOLVER: desired convergence not reached\n");
exit2:
  UserWriteF("AVG EXEC TIMES: DEF[%2d]=%10.4g JAC[%2d]=%10.4g LIN[%2d]=%10.4g\n",
             defect_c,defect_t/defect_c,newton_c,newton_t/newton_c,linear_c,linear_t/linear_c);
  res->exec_time = defect_t+newton_t+linear_t;

  /* postprocess assemble once at the end */
  if (ass->PostProcess!=NULL)
    if ((*ass->PostProcess)(ass,0,level,x,newton->d,newton->J,&error)) {
      res->error_code = __LINE__;
      REP_ERR_RETURN(res->error_code);
    }

exit:
  if (PostPCR(PrintID,NULL))                              {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* deallocate local XDATA_DESCs */
  if (FreeVD(mg,0,level,newton->d)) REP_ERR_RETURN(1);
  if (FreeVD(mg,0,level,newton->v)) REP_ERR_RETURN(1);
  if (newton->lineSearch == 3) {
    if (FreeVD(mg,0,level,newton->dold)) REP_ERR_RETURN(1);
    if (FreeVD(mg,0,level,newton->dsave)) REP_ERR_RETURN(1);
  }
  if (res->error_code==0)
    return (0);
  else
    REP_ERR_RETURN (res->error_code);
}

static INT NewtonInit (NP_BASE *base, INT argc, char **argv)
{
  NP_NEWTON *newton;                                            /* object pointer						*/
  INT i;

  newton = (NP_NEWTON *) base;

  /* read  data descs */
  newton->J = ReadArgvMatDesc(base->mg,"J",argc,argv);
  newton->d = ReadArgvVecDesc(base->mg,"d",argc,argv);
  newton->v = ReadArgvVecDesc(base->mg,"v",argc,argv);
  newton->s = ReadArgvVecDesc(base->mg,"s",argc,argv);
  newton->dold = ReadArgvVecDesc(base->mg,"dold",argc,argv);
  newton->dsave = ReadArgvVecDesc(base->mg,"dsave",argc,argv);

  /* read other numprocs */
  newton->trans = (NP_TRANSFER *) ReadArgvNumProc(base->mg,"T",TRANSFER_CLASS_NAME,
                                                  argc,argv);
  if (newton->trans == NULL) {
    PrintErrorMessage('E',"NewtonInit","cannot read transfer num proc");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  newton->solve = (NP_LINEAR_SOLVER *) ReadArgvNumProc(base->mg,"S",
                                                       LINEAR_SOLVER_CLASS_NAME,argc,argv);
  if (newton->solve == NULL) {
    PrintErrorMessage('E',"NewtonInit","cannot read solve num proc");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  /* set configuration parameters */
  if (ReadArgvDOUBLE("rhoreass",&(newton->rhoReass),argc,argv))
    newton->rhoReass = 0.8;
  if ((newton->rhoReass<0.0)||(newton->rhoReass>1.0)) {
    PrintErrorMessage('E',"NewtonInit","rhoreass must be in (0,1)");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("lsteps",&(newton->maxLineSearch),argc,argv))
    newton->maxLineSearch=6;
  if ((newton->maxLineSearch<0)||(newton->maxLineSearch>=MAX_LINE_SEARCH)) {
    PrintErrorMessageF('E',"NewtonInit","maxLineSearch < %d",(int)MAX_LINE_SEARCH);
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("line",&(newton->lineSearch),argc,argv)) {
    newton->lineSearch = 0;
    newton->maxLineSearch=1;
  }
  if ((newton->lineSearch<0)||(newton->lineSearch>3)) {
    PrintErrorMessage('E',"NewtonInit","line = 0,1,2 or 3");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("fi",&(newton->force_iteration),argc,argv)) {
    newton->force_iteration = 0;
  }
  if (ReadArgvINT("maxit",&(newton->maxit),argc,argv))
    newton->maxit = 50;
  if ((newton->maxit<0)||(newton->maxit>1000)) {
    PrintErrorMessage('E',"NewtonInit","maxit <= 1000");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("linrate",&(newton->linearRate),argc,argv))
    newton->linearRate = 0;
  if ((newton->linearRate<0)||(newton->linearRate>2)) {
    PrintErrorMessage('E',"NewtonInit","linrate = 0,1 or 2");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvDOUBLE("lambda",&(newton->lambda),argc,argv))
    newton->lambda = 1.0;
  if ((newton->lambda<-2.0)||(newton->lambda>2.0)) {
    PrintErrorMessage('E',"NewtonInit","lambda must be in (-2,2)");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (sc_read(newton->linMinRed,NP_FMT(newton),newton->s,"linminred",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      newton->linMinRed[i] = 0.001;

  for (i=0; i<MAX_VEC_COMP; i++)
    if ((newton->linMinRed[i]<0.0)||(newton->linMinRed[i]>=1.0))
    {
      PrintErrorMessage('E',"NewtonInit","linminred must be in (0,1)");
      REP_ERR_RETURN(NP_NOT_ACTIVE);
    }

  if (sc_read(newton->scale,NP_FMT(newton),newton->s,"scale",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      newton->scale[i] = 1.0;

  if (sc_read(newton->divFactor,NP_FMT(newton),newton->s,"divfac",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      newton->divFactor[i] = 1e5;

  for (i=0; i<MAX_VEC_COMP; i++)
    if ((newton->divFactor[i]<=1.0))
    {
      PrintErrorMessage('E',"NewtonInit","divfac must be in )1,inf(");
      REP_ERR_RETURN(NP_NOT_ACTIVE);
    }
  /* set noLastDef option */
  newton->noLastDef = ReadArgvOption("noLastDef",argc,argv);

  /* set display option */
  newton->displayMode = ReadArgvDisplay(argc,argv);

  /* call general nls init */
  return (NPNLSolverInit(&(newton->nlsolver),argc,argv));
}

static INT NewtonDisplay (NP_BASE *theNumProc)
{
  NP_NEWTON *newton;

  newton  = (NP_NEWTON*) theNumProc;

  /* general nls display */
  NPNLSolverDisplay(&(newton->nlsolver));

  if (newton->J != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"J",ENVITEM_NAME(newton->J));
  if (newton->v != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"v",ENVITEM_NAME(newton->v));
  if (newton->d != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"d",ENVITEM_NAME(newton->d));
  if (newton->dold != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"dold",ENVITEM_NAME(newton->dold));
  if (newton->dsave != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"dsave",ENVITEM_NAME(newton->dsave));
  if (newton->s != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(newton->s));

  if (newton->solve != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(newton->solve));
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
  UserWriteF(DISPLAY_NP_FORMAT_SI,"line",(int)newton->lineSearch);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"lsteps",(int)newton->maxLineSearch);
  if (sc_disp(newton->linMinRed,newton->s,"linMinRed")) REP_ERR_RETURN (1);
  if (sc_disp(newton->divFactor,newton->s,"divfac")) REP_ERR_RETURN (1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"lambda",(float)newton->lambda);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"rhoreass",(float)newton->rhoReass);

  return (0);
}

static INT NewtonConstruct (NP_BASE *theNP)
{
  NP_NL_SOLVER *np;

  /* set general functions */
  theNP->Init = NewtonInit;
  theNP->Display = NewtonDisplay;
  theNP->Execute = NPNLSolverExecute;

  np = (NP_NL_SOLVER *) theNP;
  np->PreProcess = NewtonPreProcess;
  np->Solver = NewtonSolver;
  np->PostProcess = NewtonPostProcess;

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

INT InitNewtonSolver (void)
{
  if (CreateClass (NL_SOLVER_CLASS_NAME ".newton",
                   sizeof(NP_NEWTON), NewtonConstruct))
    REP_ERR_RETURN (__LINE__);

  return (0);
}
