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

#include "devices.h"
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

#define MAX_LINE_SEARCH                 10
#define LINE_SEAARCH_REDUCTION  0.5

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

#undef __DEBUG_NEWTON__

/* variables for timing measurement		*/
static int defect_c, newton_c, linear_c;
static double defect_t, newton_t, linear_t;
static clock_t clock_start;

#define CSTART()        clock_start=clock()
#define CSTOP(t,c)  t+=((double)(clock()-clock_start))/((double)CLOCKS_PER_SEC);c++

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
  DOUBLE linMinRed;                             /* minimum reduction for linear solver			*/

  /* and XDATA_DESCs */
  MATDATA_DESC *J;                              /* the Matrix to be solved						*/
  VECDATA_DESC *d;                              /* nonlinear defect								*/
  VECDATA_DESC *v;                              /* correction computed by newton step                   */
  VECDATA_DESC *s;                              /* saved nonlinear solution						*/

} NP_NEWTON;    /* this is a final class */

static INT NonLinearDefect (MULTIGRID *mg, INT level, INT init, VECDATA_DESC *x, NP_NEWTON *newton, NP_NL_ASSEMBLE *ass, VEC_SCALAR defect)
{
  LRESULT lr;                           /* result of linear solver				*/
  INT i,bl,error,n_unk;

  n_unk = VD_NCOMP(x);

  /* project solution to all grid levels */
  for (i=level; i>0; i--)
  {
    if (newton->trans->PreProcessProject!=NULL)
      if ((*newton->trans->PreProcessProject)(newton->trans,i,&bl,&error)) {
        error = __LINE__;
        REP_ERR_RETURN(error);
      }
    if ((*newton->trans->ProjectSolution)(newton->trans,i,x,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }
    if (newton->trans->PostProcessProject!=NULL)
      if ((*newton->trans->PostProcessProject)(newton->trans,i,&error)) {
        error = __LINE__;
        REP_ERR_RETURN(error);
      }
  }

  if (init)
  {
    /* preprocess assemble once before all calls */
    if (ass->PreProcess!=NULL)
      if ((*ass->PreProcess)(ass,0,level,x,&error)) {
        error = __LINE__;
        REP_ERR_RETURN(error);
      }

    /* set dirichlet conditions on all grid levels */
    if ((*ass->NLAssembleSolution)(ass,0,level,x,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }
  }

  /* compute new nonlinear defect */
  CSTART();
  for (i=0; i<=level; i++) l_dset(GRID_ON_LEVEL(mg,i),newton->d,EVERY_CLASS,0.0);
  if ((*ass->NLAssembleDefect)(ass,0,level,x,newton->d,newton->J,&error)) {
    error = __LINE__;
    REP_ERR_RETURN(error);
  }
  CSTOP(defect_t,defect_c);
  if (UG_math_error) {
    UserWrite("math error in NLAssembleDefect\n");
    UG_math_error = 0;
    error = __LINE__;
    REP_ERR_RETURN(error);
  }

        #ifdef __DEBUG_NEWTON__
  UserWrite("---- After computation of nonlinear defect\n");
  ListVectorRange(mg,0,level,0,1000,FALSE,TRUE);
        #endif

  /* compute norm of defect */
  if ((*newton->solve->Residuum)(newton->solve,0,level,newton->v,newton->d,newton->J,&lr)) {
    error = __LINE__;
    REP_ERR_RETURN(error);
  }
  for (i=0; i<n_unk; i++) defect[i] = lr.last_defect[i];

  return (0);
}

/****************************************************************************/
/*D
   newton - numproc for ...

   DESCRIPTION:
   This numproc executes ...

   .vb
   npinit [$x <sol>] [$A <assemble numproc>] [$red <sc double list>]
       [$abslimit <sc double list>] ....
   .ve

   .  $x~<sol> - the solution vector
   .  $abslimit~<sc~double~list> - absolute limit for the defect
   .  $reduction~<sc~double~list> - reduction factor

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si
   <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate nl $t newton;
   npinit $x sol $A box $red 0.0001 ..............

   npexecute solver $i $d $r $s $p;
   .ve
   D*/
/****************************************************************************/

static INT NewtonPreProcess  (NP_NL_SOLVER *solve, INT level, INT *result)
{
  return(0);
}


static INT NewtonPostProcess (NP_NL_SOLVER *solve, INT level, VECDATA_DESC *x, INT *result)
{
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
  INT n_unk;                                                            /* number of components in solution		*/
  DOUBLE s, sprime, s2reach, sred;              /* combined defect norm					*/
  INT reassemble=1;                                             /* adaptive computation of jacobian		*/
  VEC_SCALAR linred;                                            /* parameters for linear solver			*/
  DOUBLE red_factor;                                            /* convergence factor for linear iter	*/
  DOUBLE la;                                                            /* damping factor in line search		*/
  DOUBLE rho[MAX_LINE_SEARCH+1];                /* reduction factors of linesearch		*/
  DOUBLE rhomin;                                                /* best reduction if !accept			*/
  INT best_ls;                                                  /* best ls if !accept					*/
  INT accept;                                                           /* line search accepted					*/
  DOUBLE Factor[MAX_VEC_COMP];                  /* for damping factor					*/
  INT bl;                                                               /* baselevel returned by preprocess		*/
  INT error;                                                            /* for return value						*/
  LRESULT lr;                                                           /* result of linear solver				*/

  /* get status */
  newton = (NP_NEWTON *) nls;      /* cast from abstract base class to final class*/
  mg = nls->base.mg;

  /* fill result variable with error condition */
  res->error_code = 0;
  res->converged = 0;
  res->rho_first = 0.0;
  res->number_of_nonlinear_iterations = 0;
  res->number_of_line_searches = 0;

  /* initialize timers and counters */
  defect_c = newton_c = linear_c = 0;
  defect_t = newton_t = linear_t = 0.0;

  /* check function pointers in numprocs */
  if (newton->trans->ProjectSolution==NULL)
  {
    UserWrite("Newton: newton->trans->ProjectSolution not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN (res->error_code);
  }
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
  if (newton->solve->Solver==NULL)
  {
    UserWrite("Newton: newton->solve->Solver not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  if (newton->solve->Residuum==NULL)
  {
    UserWrite("Newton: newton->solve->Residuum not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }

  /* dynamic XDATA_DESC allocation */
  if (AllocMDFromVD(mg,0,level,x,x,&(newton->J))) {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(newton->v))) {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(newton->d))) {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(newton->s))) {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* get number of components */
  n_unk = VD_NCOMP(x);

  /* init ass once and compute nonlinear defect */
  if (NonLinearDefect(mg,level,TRUE,x,newton,ass,defect)!=0)
  {
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }

  /* display norm of nonlinear defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(newton),'#',"\n");
  if (PreparePCR(newton->d,newton->displayMode,text,&PrintID))    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (sc_mul(defect2reach,defect,reduction,newton->d))                    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (DoPCR(PrintID,defect,PCR_CRATE))                                                    {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  for (i=0; i<n_unk; i++) res->first_defect[i] = defect[i];

  /* compute single norm */
  s = 0.0; for (i=0; i<n_unk; i++) s += defect[i]*defect[i];s = sqrt(s);
  sred = 1.0E10; for (i=0; i<n_unk; i++) sred = MIN(sred,reduction[i]);
  s2reach = s*sred;
  if (newton->lineSearch)
    UserWriteF(" ++ s=%12.4lE Initial nonlinear residual\n",s);

  /* check if iteration is necessary */
  if (sc_cmp(defect,abslimit,newton->d)) {
    res->converged = 1;
    for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
    res->error_code = 0;
    goto exit;
  }

  /* initialize reduction factor for linear solver */
  red_factor = newton->linMinRed;
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
    for (i=0; i<=level; i++)
      l_dset(GRID_ON_LEVEL(mg,i),newton->v,EVERY_CLASS,0.0);
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
    for (i=0; i<n_unk; i++) linred[i] = red_factor;
    if (newton->solve->PreProcess!=NULL)
      if ((*newton->solve->PreProcess)(newton->solve,level,newton->v,newton->d,newton->J,&bl,&error)) {
        UserWriteF("NewtonSolver: solve->PreProcess failed, error code %d\n",error);
        res->error_code = __LINE__;
        REP_ERR_RETURN(res->error_code);
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
    if ((newton->linearRate < 2) && (!lr.converged)) {
      UserWrite("\nLinear solver did not converge in Newton method\n");
      res->error_code = 0;                   /* no error but exit */
      goto exit;
    }

                #ifdef __DEBUG_NEWTON__
    UserWrite("---- After linear solver\n");
    ListVectorRange(mg,0,level,0,1000,FALSE,TRUE);
                #endif

    if (newton->lineSearch)
      UserWriteF(" ++ newton step %3d\n",r);

    /* save current solution for line search */
    a_dcopy(mg,0,level,newton->s,EVERY_CLASS,x);

    /* do a line search */
    la = newton->lambda; accept=0;
    for (kk=1; kk<=newton->maxLineSearch; kk++) {

      /* set lambda factor */
      for (i=0; i<n_unk; i++) Factor[i] = -la;

      /* update solution */
      a_dcopy(mg,0,level,x,EVERY_CLASS,newton->s);
      a_daxpy(mg,0,level,x,EVERY_CLASS,Factor,newton->v);

      if (NonLinearDefect(mg,level,FALSE,x,newton,ass,defect)!=0)
      {
        res->error_code = __LINE__;
        REP_ERR_RETURN(res->error_code);
      }

      /* compute single norm */
      sprime = 0.0; for (i=0; i<n_unk; i++) sprime += defect[i]*defect[i];
      sprime = sqrt(sprime);

      rho[kk] = sprime/s;

      /* print results */
      if (newton->lineSearch)
        UserWriteF(" ++ ls=%2d, s=%12.4lE, rho=%8.4lg, lambda= %8.4lg\n",kk,sprime,rho[kk],fabs(la));

      if (sprime/s<=1-0.25*fabs(la) || !newton->lineSearch) {
        accept=1;
        break;
      }

      /* else reduce lambda */
      la = LINE_SEAARCH_REDUCTION*la;
    }

    /* if not accepted */
    if (!accept)
      switch (newton->lineSearch)
      {
      case 1 :
        /* break iteration */
        UserWrite("line search not accepted, Newton diverged\n");
        res->error_code = 0;
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
        la = newton->lambda * pow(LINE_SEAARCH_REDUCTION,best_ls-1);
        for (i=0; i<n_unk; i++) Factor[i] = -la;

        /* update solution */
        a_dcopy(mg,0,level,x,EVERY_CLASS,newton->s);
        a_daxpy(mg,0,level,x,EVERY_CLASS,Factor,newton->v);

        if (NonLinearDefect(mg,level,FALSE,x,newton,ass,defect)!=0)
        {
          res->error_code = __LINE__;
          REP_ERR_RETURN(res->error_code);
        }
        break;

        /*default: accept */
      }


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

    /* compute new reduction factor, assuming quadratic convergence */
    red_factor = MIN((sprime/s)*(sprime/s),newton->linMinRed);
    if (newton->linearRate == 1)
      red_factor = MIN(sprime/s,newton->linMinRed);
    if (newton->linearRate == 2)
      red_factor = newton->linMinRed;

    /* accept new iterate */
    s = sprime;
  }

  /* print norm of defect */
  if (DoPCR(PrintID,defect,PCR_AVERAGE))  {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* if converged, then report results and mean execution times */
  if (res->converged) {
    res->error_code = 0;
    res->number_of_nonlinear_iterations = newton_c;
    res->number_of_line_searches = defect_c;
    for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
    UserWriteF("AVG EXEC TIMES: DEF[%2d]=%10.4lg JAC[%2d]=%10.4lg LIN[%2d]=%10.4lg\n",
               defect_c,defect_t/defect_c,newton_c,newton_t/newton_c,linear_c,linear_t/linear_c);
  }

  /* postprocess assemble once at the end */
  if (ass->PostProcess!=NULL)
    if ((*ass->PostProcess)(ass,0,level,x,newton->d,newton->J,&error)) {
      res->error_code = __LINE__;
      REP_ERR_RETURN(res->error_code);
    }

exit:
  if (PostPCR(PrintID,NULL))                              {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* deallocate local XDATA_DESCs */
  FreeVD(mg,0,level,newton->d);
  FreeVD(mg,0,level,newton->v);
  FreeVD(mg,0,level,newton->s);
  FreeMD(mg,0,level,newton->J);

  if (res->error_code==0)
    return (0);
  else
    REP_ERR_RETURN (res->error_code);
}

static INT NewtonInit (NP_BASE *base, INT argc, char **argv)
{
  NP_NEWTON *newton;                                            /* object pointer						*/

  newton = (NP_NEWTON *) base;

  /* read  data descs */
  newton->J = ReadArgvMatDesc(base->mg,"J",argc,argv);
  newton->d = ReadArgvVecDesc(base->mg,"d",argc,argv);
  newton->v = ReadArgvVecDesc(base->mg,"v",argc,argv);
  newton->s = ReadArgvVecDesc(base->mg,"s",argc,argv);

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
  if ((newton->lineSearch<0)||(newton->lineSearch>2)) {
    PrintErrorMessage('E',"NewtonInit","line = 0 or 2");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
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
  if (ReadArgvDOUBLE("linminred",&(newton->linMinRed),argc,argv))
    newton->linMinRed = 0.001;
  if ((newton->linMinRed<0.0)||(newton->linMinRed>=1.0)) {
    PrintErrorMessage('E',"NewtonInit","linminred must be in (0,1)");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
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
  UserWriteF(DISPLAY_NP_FORMAT_SF,"linminred",(float)newton->linMinRed);
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
