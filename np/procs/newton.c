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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

#undef __DEBUG_NEWTON__

static clock_t clock_start;                     /* used for internal timing					*/

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
  int defect_c, newton_c, linear_c;             /* variables for timing measurement		*/
  double defect_t, newton_t, linear_t;      /* variables for timing measurement		*/
  char text[DISPLAY_WIDTH+4];                           /* display text in PCR					*/
  INT PrintID;                                                  /* print id for PCR						*/
  VEC_SCALAR defect, defect2reach;              /* component--wise norm					*/
  INT n_unk;                                                            /* number of components in solution		*/
  DOUBLE s, sprime, s2reach, sred;              /* combined defect norm					*/
  char buffer[128];                                             /* for messages							*/
  INT reassemble=1;                                             /* adaptive computation of jacobian		*/
  VEC_SCALAR linred;                                            /* parameters for linear solver			*/
  DOUBLE red_factor;                                            /* convergence factor for linear iter	*/
  DOUBLE la;                                                            /* damping factor in line search		*/
  INT accept;                                                           /* line search accepted					*/
  DOUBLE Factor[MAX_VEC_COMP];                  /* for damping factor					*/
  INT bl;                                                               /* baselevel returned by preprocess		*/
  INT error;                                                            /* for return value						*/
  LRESULT lr;                                                           /* result of linear solver				*/

  /* get status */
  newton = (NP_NEWTON *) nls;      /* cast from abstract base class to final class*/
  mg = nls->base.mg;

  /* fill result variable with error condition */
  res->error_code = __LINE__;
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
    return(res->error_code);
  }
  if (ass->NLAssembleSolution==NULL)
  {
    UserWrite("Newton: ass->NLAssembleSolution not defined\n");
    res->error_code = __LINE__;
    return(res->error_code);
  }
  if (ass->NLAssembleDefect==NULL)
  {
    UserWrite("Newton: ass->NLAssembleDefect not defined\n");
    res->error_code = __LINE__;
    return(res->error_code);
  }
  if (ass->NLAssembleMatrix==NULL)
  {
    UserWrite("Newton: ass->NLAssembleMatrix not defined\n");
    res->error_code = __LINE__;
    return(res->error_code);
  }
  if (newton->solve->Solver==NULL)
  {
    UserWrite("Newton: newton->solve->Solver not defined\n");
    res->error_code = __LINE__;
    return(res->error_code);
  }
  if (newton->solve->Residuum==NULL)
  {
    UserWrite("Newton: newton->solve->Residuum not defined\n");
    res->error_code = __LINE__;
    return(res->error_code);
  }

  /* dynamic XDATA_DESC allocation */
  if (AllocMDFromVD(mg,0,level,x,x,&(newton->J))) {res->error_code = __LINE__; return(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(newton->v))) {res->error_code = __LINE__; return(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(newton->d))) {res->error_code = __LINE__; return(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(newton->s))) {res->error_code = __LINE__; return(res->error_code);}

  /* project initial guess to all grid levels */
  for (i=level; i>0; i--)
  {
    if (newton->trans->PreProcessProject!=NULL)
      if ((*newton->trans->PreProcessProject)(newton->trans,i,&bl,&error)) {
        res->error_code = __LINE__;
        return(res->error_code);
      }
    if ((*newton->trans->ProjectSolution)(newton->trans,i,x,&error)) {
      res->error_code = __LINE__;
      return(res->error_code);
    }
    if (newton->trans->PostProcessProject!=NULL)
      if ((*newton->trans->PostProcessProject)(newton->trans,i,&error)) {
        res->error_code = __LINE__;
        return(res->error_code);
      }
  }

  /* preprocess assemble once before all calls */
  if (ass->PreProcess!=NULL)
    if ((*ass->PreProcess)(ass,0,level,x,&error)) {
      res->error_code = __LINE__;
      return(res->error_code);
    }

  /* set dirichlet conditions on all grid levels */
  if ((*ass->NLAssembleSolution)(ass,0,level,x,&error)) {
    res->error_code = __LINE__;
    return(res->error_code);
  }

  /* get number of components */
  n_unk = VD_NCOMP(x);

  /* assemble defect from current solution on all levels */
  CSTART();
  for (i=0; i<=level; i++) l_dset(GRID_ON_LEVEL(mg,i),newton->d,EVERY_CLASS,0.0);
  if ((*ass->NLAssembleDefect)(ass,0,level,x,newton->d,newton->J,&error)) {
    res->error_code = __LINE__;
    return(res->error_code);
  }
  CSTOP(defect_t,defect_c);

        #ifdef __DEBUG_NEWTON__
  UserWrite("===== After computation of defect\n");
  /* ListVectorRange(mg,0,level,0,1000,TRUE,TRUE); */
        #endif

  /* compute norm of nonlinear defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(newton),'*',"\n");
  if (PreparePCR(newton->d,newton->displayMode,text,&PrintID))    {res->error_code = __LINE__; return(res->error_code);}
  if (s_eunorm(mg,0,level,newton->d,defect))                                              {res->error_code = __LINE__; return(res->error_code);}
  if (sc_mul(defect2reach,defect,reduction,newton->d))                    {res->error_code = __LINE__; return(res->error_code);}
  if (DoPCR(PrintID,defect,PCR_CRATE))                                                    {res->error_code = __LINE__; return(res->error_code);}
  for (i=0; i<n_unk; i++) res->first_defect[i] = defect[i];

  /* compute single norm */
  s = 0.0; for (i=0; i<n_unk; i++) s += defect[i]*defect[i];s = sqrt(s);
  sred = 1.0E10; for (i=0; i<n_unk; i++) sred = MIN(sred,nls->reduction[i]);
  s2reach = s*sred;
  if (newton->lineSearch)
  {
    sprintf(buffer," ++ s=%12.4lE Initial nonlinear residual\n",s);
    UserWrite(buffer);
  }

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
      for (i=0; i<=level; i++)
        l_dmatset(GRID_ON_LEVEL(mg,i),newton->J,0.0);
      if ((*ass->NLAssembleMatrix)(ass,0,level,x,newton->d,newton->v,newton->J,&error)) {
        res->error_code = __LINE__;
        return(res->error_code);
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
        res->error_code = __LINE__;
        return(res->error_code);
      }
    if ((*newton->solve->Residuum)(newton->solve,level,newton->v,newton->d,newton->J,&lr))
    {
      res->error_code = __LINE__;
      goto exit;
    }
    if ((*newton->solve->Solver)(newton->solve,level,newton->v,newton->d,newton->J,abslimit,linred,&lr))
    {
      res->error_code = __LINE__;
      goto exit;
    }
    if (newton->solve->PostProcess!=NULL)
      if ((*newton->solve->PostProcess)(newton->solve,level,newton->v,newton->d,newton->J,&error)) {
        res->error_code = __LINE__;
        return(res->error_code);
      }
    CSTOP(linear_t,linear_c);
    if (UG_math_error) {
      UG_math_error = 0;
      res->error_code = __LINE__;
    }

    /* if linear solver did not converge, return here */
    if (!lr.converged) {
      UserWrite("\nLinear solver did not converge in Newton method\n");
      res->error_code = __LINE__;
      goto exit;
    }

                #ifdef __DEBUG_NEWTON__
    UserWrite("===== After linear solver\n");
    ListVectorRange(mg,0,level,0,1000,FALSE,TRUE);
                #endif

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

      /* project solution to all grid levels */
      for (i=level; i>0; i--)
      {
        if (newton->trans->PreProcessProject!=NULL)
          if ((*newton->trans->PreProcessProject)(newton->trans,i,&bl,&error)) {
            res->error_code = __LINE__;
            return(res->error_code);
          }
        if ((*newton->trans->ProjectSolution)(newton->trans,i,x,&error)) {
          res->error_code = __LINE__;
          return(res->error_code);
        }
        if (newton->trans->PostProcessProject!=NULL)
          if ((*newton->trans->PostProcessProject)(newton->trans,i,&error)) {
            res->error_code = __LINE__;
            return(res->error_code);
          }
      }

      /* compute new nonlinear defect with given lambda */
      CSTART();
      for (i=0; i<=level; i++) l_dset(GRID_ON_LEVEL(mg,i),newton->d,EVERY_CLASS,0.0);
      if ((*ass->NLAssembleDefect)(ass,0,level,x,newton->d,newton->J,&error)) {
        res->error_code = __LINE__;
        return(res->error_code);
      }
      CSTOP(defect_t,defect_c);
      if (UG_math_error) {
        UserWrite("math error in NLAssembleDefect\n");
        UG_math_error = 0;
        res->error_code = __LINE__;
        break;
      }

                        #ifdef __DEBUG_NEWTON__
      UserWrite("===== After computation of nonlinear defect\n");
      ListVectorRange(mg,0,level,0,1000,FALSE,TRUE);
                        #endif

      /* compute norm of defect */
      if (s_eunorm(mg,0,level,newton->d,defect))      {res->error_code = __LINE__; return(res->error_code);}

      /* compute single norm */
      sprime = 0.0; for (i=0; i<n_unk; i++) sprime += defect[i]*defect[i];
      sprime = sqrt(sprime);

      /* print results */
      if (newton->lineSearch)
      {
        sprintf(buffer," ++ ls=%2d, s=%12.4lE, rho=%12.4lg, newton step %3d\n",kk,sprime,sprime/s,r);
        UserWrite(buffer);
      }

      if (sprime/s<=1-0.25*la || !newton->lineSearch) {
        accept=1;
        break;
      }

      /* else reduce lambda */
      la = 0.5*la;
    }

    /* reduce time step if not accepted */
    if (!accept) {
      UserWrite("line search not accepted, Newton diverged\n");
      res->error_code = __LINE__;
      break;
    }

    /* print norm of defect */
    if (DoPCR(PrintID,defect,PCR_CRATE)) {res->error_code = __LINE__; return(res->error_code);}

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
    if (newton->linearRate) red_factor = MIN(sprime/s,newton->linMinRed);

    /* accept new iterate */
    s = sprime;
  }

  /* print norm of defect */
  if (DoPCR(PrintID,defect,PCR_AVERAGE))  {res->error_code = __LINE__; return(res->error_code);}
  if (PostPCR(PrintID,NULL))                              {res->error_code = __LINE__; return(res->error_code);}

  /* if converged, then accept new time step */
  if (res->converged) {
    res->error_code = 0;
    res->number_of_nonlinear_iterations = newton_c;
    res->number_of_line_searches = defect_c;
    for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
    sprintf(buffer,"AVG EXEC TIMES: DEF[%2d]=%10.4lg JAC[%2d]=%10.4lg LIN[%2d]=%10.4lg\n",
            defect_c,defect_t/defect_c,newton_c,newton_t/newton_c,linear_c,linear_t/linear_c);
    UserWrite(buffer);
  }

  /* postprocess assemble once at the end */
  if (ass->PostProcess!=NULL)
    if ((*ass->PostProcess)(ass,0,level,x,newton->d,newton->J,&error)) {
      res->error_code = __LINE__;
      return(res->error_code);
    }

exit:
  /* deallocate local XDATA_DESCs */
  FreeVD(mg,0,level,newton->d);
  FreeVD(mg,0,level,newton->v);
  FreeVD(mg,0,level,newton->s);
  FreeMD(mg,0,level,newton->J);

  return(res->error_code);
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
  if (newton->trans == NULL) return(NP_NOT_ACTIVE);
  newton->solve = (NP_LINEAR_SOLVER *) ReadArgvNumProc(base->mg,"S",
                                                       LINEAR_SOLVER_CLASS_NAME,argc,argv);
  if (newton->solve == NULL) return(NP_NOT_ACTIVE);

  /* set configuration parameters */
  if (ReadArgvDOUBLE("rhoreass",&(newton->rhoReass),argc,argv))
    newton->rhoReass = 0.8;
  if ((newton->rhoReass<0.0)||(newton->rhoReass>1.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvINT("lsteps",&(newton->maxLineSearch),argc,argv))
    newton->maxLineSearch=6;
  if ((newton->maxLineSearch<0)||(newton->maxLineSearch>20)) return(NP_NOT_ACTIVE);

  if (ReadArgvINT("maxit",&(newton->maxit),argc,argv))
    newton->maxit = 50;
  if ((newton->maxit<0)||(newton->maxit>1000)) return(NP_NOT_ACTIVE);

  if (ReadArgvINT("line",&(newton->lineSearch),argc,argv))
    newton->lineSearch = 1;
  if ((newton->lineSearch<0)||(newton->lineSearch>1)) return(NP_NOT_ACTIVE);

  if (ReadArgvINT("linrate",&(newton->linearRate),argc,argv))
    newton->linearRate = 0;
  if ((newton->linearRate<0)||(newton->linearRate>1)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("lambda",&(newton->lambda),argc,argv))
    newton->lambda = 1.0;
  if ((newton->lambda<0.0)||(newton->lambda>2.0)) return(NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("linminred",&(newton->linMinRed),argc,argv))
    newton->linMinRed = 0.001;
  if ((newton->linMinRed<0.0)||(newton->linMinRed>=1.0)) return(NP_NOT_ACTIVE);


  /* set display option */
  newton->displayMode = ReadArgvDisplay(argc,argv);

  /* call general nls init */
  return (NPNLSolverInit(&(newton->nlsolver),argc,argv));
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

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  Execute														*/
/*																			*/
/* Purpose:	  execute linear solver											*/
/*																			*/
/* Input:	  theNumProc	num proc data structure							*/
/*			  theMG			data newton										*/
/*			  argc,argv		command options									*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*			  else : error													*/
/*																			*/
/****************************************************************************/


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
    return (__LINE__);

  return (0);
}
