// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ls.c	                                                                                                        */
/*																			*/
/* Purpose:   linear solver num procs                                           */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   November 29, 1996                                                                         */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "devices.h"
#include "general.h"
#include "debug.h"
#include "ugstruct.h"
#include "gm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"

#include "iter.h"
#include "ls.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define ABS_LIMIT 1e-10

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct np_ls
{
  NP_LINEAR_SOLVER ls;

  NP_ITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;

  VECDATA_DESC *c;

  INT (*Prepare)
    (struct np_ls *,                         /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*Update)
    (struct np_ls *,                         /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* correction vector               */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*Close)
    (struct np_ls *,                         /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_ls NP_LS;

typedef struct
{
  NP_LS ls;

  DOUBLE rho;
  VECDATA_DESC *p;
  VECDATA_DESC *t;

} NP_CG;

typedef struct
{
  NP_LINEAR_SOLVER ls;

  NP_ITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;
  INT restart;

  DOUBLE sp;
  VEC_SCALAR weight;
  VECDATA_DESC *p;
  VECDATA_DESC *pp;
  VECDATA_DESC *t;
  VECDATA_DESC *h1;
  VECDATA_DESC *h2;
  VECDATA_DESC *h3;

} NP_CR;

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static VEC_SCALAR Factor_One;

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   NP_LINEAR_SOLVER - type definition for linear solvers

   DESCRIPTION:
   This numproc type is used for the description of linear solvers.
   It can be called by the given interface from a nonlinearsolver.
   Initializing the data is optional; it can be done with

   'INT NPLinearSolverInit (NP_LINEAR_SOLVER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPLinearSolverDisplay (NP_LINEAR_SOLVER *theNP);'
   'INT NPLinearSolverExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb


   ..... fill in data structure here when the realizition is finished


   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPLinearSolverInit (NP_LINEAR_SOLVER *np, INT argc , char **argv)
{
  INT i;

  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);
  if (sc_read(np->abslimit,np->x,"abslimit",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->abslimit[i] = ABS_LIMIT;
  if (sc_read(np->reduction,NULL,"red",argc,argv))
    return(NP_ACTIVE);

  if ((np->x == NULL) || (np->b == NULL) || (np->A == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPLinearSolverDisplay (NP_LINEAR_SOLVER *np)
{
  if ((np->x != NULL) || (np->b != NULL) || (np->A != NULL)) {
    UserWrite("symbolic user data:\n");
    if (np->A != NULL)
      UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
    if (np->x != NULL)
      UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
    if (np->b != NULL)
      UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->b));
    UserWrite("\n");
  }
  UserWrite("configuration parameters:\n");
  if (np->x != NULL)
    if (sc_disp(np->reduction,np->x,"red"))
      REP_ERR_RETURN (1);
  if (sc_disp(np->abslimit,np->x,"abslimit"))
    REP_ERR_RETURN (1);

  return(0);
}

INT NPLinearSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LINEAR_SOLVER *np;
  LRESULT lresult;
  INT result,level,bl;

  np = (NP_LINEAR_SOLVER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no vector x");
    REP_ERR_RETURN (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no vector b");
    REP_ERR_RETURN (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no matrix A");
    REP_ERR_RETURN (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no PreProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PreProcess)(np,level,np->x,np->b,np->A,&bl,&result)) {
      UserWriteF("NPLinearSolverExecute: PreProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("d",argc,argv)) {
    if (np->Defect == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Defect");
      REP_ERR_RETURN (1);
    }
    if ((*np->Defect)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPLinearSolverExecute: Defect failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("r",argc,argv)) {
    if (np->Residuum == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Residuum");
      REP_ERR_RETURN (1);
    }
    if ((*np->Residuum)(np,bl,level,np->x,np->b,np->A,&lresult)) {
      UserWriteF("NPLinearSolverExecute: Residuum failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->Solver == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Solver");
      REP_ERR_RETURN (1);
    }
    if ((*np->Solver)(np,level,np->x,np->b,np->A,
                      np->abslimit,np->reduction,&lresult)) {
      UserWriteF("NPLinearSolverExecute: Solver failed, error code %d\n",
                 lresult.error_code);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no PostProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PostProcess)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPLinearSolverExecute: PostProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }
  return(0);
}

/* tools for linear solvers */

static INT LinearSolverInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    REP_ERR_RETURN(NP_NOT_ACTIVE);

  np->display = ReadArgvDisplay(argc,argv);
  np->Iter = (NP_ITER *)
             ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);
  if (np->Iter == NULL)
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->baselevel = 0;
  np->c = ReadArgvVecDesc(theNP->mg,"c",argc,argv);

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT LinearSolverDisplay (NP_BASE *theNP)
{
  NP_LS *np;

  np = (NP_LS *) theNP;
  NPLinearSolverDisplay(&np->ls);

  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->Iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");
  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->c));

  return (0);
}

static INT LinearSolverPreProcess (NP_LINEAR_SOLVER *theNP, INT level,
                                   VECDATA_DESC *x, VECDATA_DESC *b,
                                   MATDATA_DESC *A,
                                   INT *baselevel, INT *result)
{
  NP_LS *np;

  /* store passed XXXDATA_DESCs */
  NPLS_A(theNP) = A;
  NPLS_x(theNP) = x;
  NPLS_b(theNP) = b;

  np = (NP_LS *) theNP;
  if (np->Iter->PreProcess != NULL)
    if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result))
      REP_ERR_RETURN(1);
  np->baselevel = MIN(*baselevel,level);

  return(0);
}

static INT LinearDefect (NP_LINEAR_SOLVER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *result)
{
  NP_LS *np;

  np = (NP_LS *) theNP;
  if (s_dmatmul_minus(theNP->base.mg,np->baselevel,level,b,A,x,EVERY_CLASS)
      != NUM_OK) NP_RETURN(1,result[0]);
  return (0);
}

static INT LinearResiduum (NP_LINEAR_SOLVER *theNP, INT bl, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           LRESULT *lresult)
{
  NP_LS *np;

  np = (NP_LS *) theNP;
        #ifdef ModelP
  if (a_vector_collect(theNP->base.mg,bl,level,b)) NP_RETURN(1,lresult->error_code);
        #endif
  if (s_eunorm(theNP->base.mg,bl,level,b,lresult->last_defect)) NP_RETURN(1,lresult->error_code);

  return(0);
}

static INT LinearSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_LS *np;
  VEC_SCALAR defect2reach;
  INT i,bl,PrintID;
  char text[DISPLAY_WIDTH+4];

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  np = (NP_LS *) theNP;
  bl = np->baselevel;
  if (np->Iter->Iter == NULL) NP_RETURN(1,lresult->error_code);
  if (np->Update == NULL) NP_RETURN(1,lresult->error_code);
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->c)) NP_RETURN(1,lresult->error_code);
  if (np->Prepare != NULL)
    if ((*np->Prepare)(np,level,x,&lresult->error_code)) REP_ERR_RETURN (1);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY)
    if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY)
    if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;
    if (l_dset(GRID_ON_LEVEL(theNP->base.mg,level),np->c,EVERY_CLASS,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    if ((*np->Iter->Iter)(np->Iter,level,np->c,b,A,&lresult->error_code)) REP_ERR_RETURN (1);
    if ((*np->Update)(np,level,x,np->c,b,A,&lresult->error_code)) REP_ERR_RETURN (1);
    if (LinearResiduum(theNP,bl,level,x,b,A,lresult)) REP_ERR_RETURN(1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  FreeVD(theNP->base.mg,bl,level,np->c);
  if (np->Close != NULL)
    if ((*np->Close)(np,level,&lresult->error_code))
      REP_ERR_RETURN (1);
  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
  }

  return (0);
}

static INT LinearSolverPostProcess (NP_LINEAR_SOLVER *theNP,
                                    INT level,
                                    VECDATA_DESC *x, VECDATA_DESC *b,
                                    MATDATA_DESC *A,
                                    INT *result)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

  if (np->Iter->PostProcess == NULL)
    return(0);
  return((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result));
}

/****************************************************************************/
/*D
   ls - numproc for linear solvers

   DESCRIPTION:
   This numproc executes a linear solver: it performs an iteration
   up to convergence.

   .vb
   npinit [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $m~<maxit> - maximal number of iterations
   .  $I~<iteration> - iteration numproc
   .  $d - display modus

   'npexecute <name> [$i] [$d] [$r] [$s] [$p]'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate pre $c ilu;           npinit pre;
   npcreate post $c ilu;          npinit post;
   npcreate base $c ilu;          npinit base $n 3;
   npcreate basesolver $c ls;     npinit basesolver $red 0.001 $I base;
   npcreate transfer $c transfer; npinit transfer;
   npcreate lmgc $c lmgc;         npinit lmgc $S pre post basesolver $T transfer;
   npcreate mgs $c ls;            npinit mgs $A MAT $x sol $b rhs
                                          $red 0.00001 $I lmgc $d full;
   .ve
   D*/
/****************************************************************************/

static INT LSUpdate (NP_LS *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *c,
                     VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  if (a_daxpy (theNP->ls.base.mg,theNP->baselevel,level,
               x,EVERY_CLASS,Factor_One,c) != NUM_OK) NP_RETURN(1,result[0]);

  return(0);
}

static INT LSConstruct (NP_BASE *theNP)
{
  NP_LS *np;

  theNP->Init = LinearSolverInit;
  theNP->Display = LinearSolverDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_LS *) theNP;
  np->ls.PreProcess = LinearSolverPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = LinearSolver;
  np->ls.PostProcess = LinearSolverPostProcess;

  np->Prepare = NULL;
  np->Update = LSUpdate;
  np->Close = NULL;

  return(0);
}

/****************************************************************************/
/*D
   cg - numproc for the conjugate gradient method

   DESCRIPTION:
   This numproc executes a conjugate gradient step. It is preconditioned
   by an iteration numproc, e. g. a multi grid cycle or a smoother.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $I <iteration> [$d {full|red|no}]
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $P~<iteration> - preconditioner
   .  $d - display modus

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate pre $c ilu;           npinit pre;
   npcreate post $c ilu;          npinit post;
   npcreate base $c ilu;          npinit base $n 3;
   npcreate basesolver $c cg;     npinit basesolver $red 0.001 $I base;
   npcreate transfer $c transfer; npinit transfer;
   npcreate lmgc $c lmgc;         npinit lmgc $S pre post basesolver $T transfer;
   npcreate mgs $c cg;            npinit mgs $A MAT $x sol $b rhs
                                          $red 0.00001 $I lmgc $d full;
   .ve
   D*/
/****************************************************************************/

static INT CGInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  np->p = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);

  return (LinearSolverInit(theNP,argc,argv));
}

static INT CGDisplay (NP_BASE *theNP)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  LinearSolverDisplay(theNP);
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->p != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));

  return (0);
}

static INT CGPrepare (NP_LS *theNP, INT level, VECDATA_DESC *x, INT *result)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  if (AllocVDFromVD(theNP->ls.base.mg,theNP->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (a_dset(theNP->ls.base.mg,theNP->baselevel,level,np->p,EVERY_CLASS,0.0)
      != NUM_OK) NP_RETURN(1,result[0]);
  np->rho = 1.0;

  return(0);
}

static INT CGUpdate (NP_LS *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *c,
                     VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_CG *np;
  MULTIGRID *theMG;
  VEC_SCALAR scal;
  DOUBLE lambda;
  INT ncomp,j;

  np = (NP_CG *) theNP;
  theMG = theNP->ls.base.mg;
  ncomp = VD_NCOMP(x);
  if (AllocVDFromVD(theMG,theNP->baselevel,level,x,&np->t)) NP_RETURN(1,result[0]);
  if (a_dset(theMG,theNP->baselevel,level,np->t,EVERY_CLASS,0.0) != NUM_OK) NP_RETURN(1,result[0]);
  for (j=theNP->baselevel; j<=level; j++)
    if (l_dmatmul(GRID_ON_LEVEL(theMG,j),np->t,EVERY_CLASS,A,c,EVERY_CLASS)
        !=NUM_OK) NP_RETURN(1,result[0]);
  if (a_daxpy(theMG,theNP->baselevel,level,b,EVERY_CLASS,Factor_One,np->t)) NP_RETURN(1,result[0]);
  if (s_ddot(theMG,theNP->baselevel,level,c,b,scal) !=NUM_OK) NP_RETURN(1,result[0]);
  lambda = 0.0;
  for (j=0; j<ncomp; j++) lambda += scal[j];
  for (j=0; j<ncomp; j++) scal[j] = lambda / np->rho;
  np->rho = lambda;
  if (a_dscale(theMG,theNP->baselevel,level,np->p,EVERY_CLASS,scal)
      != NUM_OK) NP_RETURN(1,result[0]);
  if (a_daxpy (theMG,theNP->baselevel,level,np->p,EVERY_CLASS,Factor_One,c)
      != NUM_OK) NP_RETURN(1,result[0]);
  if (a_dset(theMG,theNP->baselevel,level,np->t,EVERY_CLASS,0.0) != NUM_OK) NP_RETURN(1,result[0]);
  for (j=theNP->baselevel; j<=level; j++)
    if (l_dmatmul (GRID_ON_LEVEL(theMG,j),np->t,EVERY_CLASS,
                   A,np->p,EVERY_CLASS) != NUM_OK) NP_RETURN(1,result[0]);
  if (s_ddot (theMG,theNP->baselevel,level,np->t,np->p,scal) != NUM_OK) NP_RETURN(1,result[0]);
  lambda = 0.0;
  for (j=0; j<ncomp; j++) lambda += scal[j];
  for (j=0; j<ncomp; j++) scal[j] = np->rho / lambda;
  if (a_daxpy(theMG,theNP->baselevel,level,x,EVERY_CLASS,scal,np->p)
      != NUM_OK) NP_RETURN(1,result[0]);
  for (j=0; j<ncomp; j++) scal[j] = - np->rho / lambda;
  if (a_daxpy (theMG,theNP->baselevel,level,b,EVERY_CLASS,scal,np->t)
      != NUM_OK) NP_RETURN(1,result[0]);
  FreeVD(theNP->ls.base.mg,theNP->baselevel,level,np->t);
  if (theNP->display == PCR_FULL_DISPLAY)
    UserWriteF("      rho %-.4g \n",np->rho);

  return(0);
}

static INT CGClose (NP_LS *theNP, INT level, INT *result)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  FreeVD(theNP->ls.base.mg,theNP->baselevel,level,np->p);

  return(0);
}

static INT CGConstruct (NP_BASE *theNP)
{
  NP_LS *np;

  theNP->Init = CGInit;
  theNP->Display = CGDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_LS *) theNP;
  np->ls.PreProcess = LinearSolverPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = LinearSolver;
  np->ls.PostProcess = LinearSolverPostProcess;

  np->Prepare = CGPrepare;
  np->Update = CGUpdate;
  np->Close = CGClose;

  return(0);
}

/****************************************************************************/
/*D
   cr - numproc for the conjugate residuum method

   DESCRIPTION:
   This numproc executes a conjugate residuum step. It is preconditioned
   by an iteration numproc, e. g. a multi grid cycle or a smoother.

   .vb
   npinit [$x <sol>] [$b <rhs>] [$A <mat>]
       [$p <condir>] [$h1 <help1>] [$h2 <help2>]
       $I <iteration> [$d {full|red|no}]
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $P~<iteration> - preconditioner
   .  $d - display modus

   'npexecute <name> [$i] [$d] [$r] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate pre $c ilu;           npinit pre;
   npcreate post $c ilu;          npinit post;
   npcreate base $c ilu;          npinit base $n 3;
   npcreate basesolver $c cg;     npinit basesolver $red 0.001 $I base;
   npcreate transfer $c transfer; npinit transfer;
   npcreate lmgc $c lmgc;         npinit lmgc $S pre post basesolver $T transfer;
   npcreate mgs $c cr;            npinit mgs $A MAT $x sol $b rhs
                                          $red 0.00001 $I lmgc $d full;
   .ve
   D*/
/****************************************************************************/

static INT CRInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_CR *np;
  INT i;

  np = (NP_CR *) theNP;
  if (sc_read (np->weight,NULL,"weight",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] = 1.0;
  np->p  = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->pp = ReadArgvVecDesc(theNP->mg,"pp",argc,argv);
  np->h1 = ReadArgvVecDesc(theNP->mg,"h1",argc,argv);
  np->h2 = ReadArgvVecDesc(theNP->mg,"h2",argc,argv);
  np->h3 = ReadArgvVecDesc(theNP->mg,"h3",argc,argv);
  np->t  = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("r",&(np->restart),argc,argv))
    np->restart = 0;
  if (np->restart<0) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->Iter = (NP_ITER *) ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);
  if (np->Iter == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->baselevel = 0;

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT CRDisplay (NP_BASE *theNP)
{
  NP_CR *np;

  np = (NP_CR *) theNP;
  NPLinearSolverDisplay(&np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"r",(int)np->restart);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->Iter != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");
  if (np->display == PCR_NO_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->p != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  if (np->pp != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"pp",ENVITEM_NAME(np->pp));
  if (np->h1 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h1",ENVITEM_NAME(np->h1));
  if (np->h2 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h2",ENVITEM_NAME(np->h2));
  if (np->h3 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h3",ENVITEM_NAME(np->h3));
  if (np->t != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h3",ENVITEM_NAME(np->t));
  if (np->p != NULL) if (sc_disp(np->weight,np->p,"weight")) REP_ERR_RETURN (1);

  return (0);
}

static INT CRPreProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_CR *np;
  INT i;

  np = (NP_CR *) theNP;
  if (np->Iter->PreProcess != NULL)
    if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);
  np->baselevel = MIN(*baselevel,level);

  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->pp)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->t)) NP_RETURN(1,result[0]);

  return(0);
}

static INT CRPostProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_CR *np;

  np = (NP_CR *) theNP;
  FreeVD(np->ls.base.mg,np->baselevel,level,np->p);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->pp);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->t);

  if (np->Iter->PostProcess == NULL) return(0);
  return((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result));

  return(0);
}

void PrintEunorm (MULTIGRID *theMG, VECDATA_DESC *v, char *name)
{
  DOUBLE eu;

  s_eunorm(theMG,0,0,v,&eu);
  UserWriteF("EUNORM(%s): %f\n",name,(float)eu);

  return;
}

static INT CRSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_CR *np;
  VEC_SCALAR defect2reach,scal,alpha0, alpha1;
  INT i,j,bl,PrintID;
  char text[DISPLAY_WIDTH+4];
  DOUBLE s,t;

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  /* prepare */
  np = (NP_CR *) theNP;
  bl = np->baselevel;
  if (np->Iter->Iter == NULL) NP_RETURN(1,lresult->error_code);
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->h1)) NP_RETURN(1,lresult->error_code);
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->h2)) NP_RETURN(1,lresult->error_code);
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->h3)) NP_RETURN(1,lresult->error_code);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY) if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY) if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;

  if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  if (s_dset(theNP->base.mg,np->baselevel,level,np->p,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  if ((*np->Iter->Iter)(np->Iter,level,np->p,np->h1,A,&lresult->error_code)) REP_ERR_RETURN (1);
  if (s_dset(theNP->base.mg,np->baselevel,level,np->pp,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  if (s_dset(theNP->base.mg,np->baselevel,level,np->t,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  np->sp = 1.0;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;

    /* update x, b */
    if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->h2,A,np->p,EVERY_CLASS)) REP_ERR_RETURN (1);
    if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,np->h2)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    if (s_dset(theNP->base.mg,np->baselevel,level,np->h3,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    if ((*np->Iter->Iter)(np->Iter,level,np->h3,np->h2,A,&lresult->error_code)) REP_ERR_RETURN (1);
    if (s_ddot (theNP->base.mg,np->baselevel,level,np->h1,np->h3,scal)!=NUM_OK) REP_ERR_RETURN (1);
    s=0.0; for (j=0; j<VD_NCOMP(x); j++) s += np->weight[j]*scal[j];
    if (s_ddot (theNP->base.mg,np->baselevel,level,b,np->h3,scal)!=NUM_OK) REP_ERR_RETURN (1);
    t=0.0; for (j=0; j<VD_NCOMP(x); j++) t += np->weight[j]*scal[j];
    for (j=0; j<VD_NCOMP(x); j++) scal[j] = t/s;
    if (s_daxpy (theNP->base.mg,np->baselevel,level,x,scal,np->p)!= NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) scal[j] = -scal[j];
    if (s_daxpy (theNP->base.mg,np->baselevel,level,b,scal,np->h1)!= NUM_OK) REP_ERR_RETURN (1);

    /* update p */
    if (np->restart>0 && i%np->restart==0)
    {
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dset(theNP->base.mg,np->baselevel,level,np->p,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if ((*np->Iter->Iter)(np->Iter,level,np->p,np->h1,A,&lresult->error_code)) REP_ERR_RETURN (1);
      if (s_dset(theNP->base.mg,np->baselevel,level,np->pp,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dset(theNP->base.mg,np->baselevel,level,np->t,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      np->sp = 1.0;
    }
    else
    {
      if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->h2,A,np->h3,EVERY_CLASS)) REP_ERR_RETURN (1);
      if (s_ddot (theNP->base.mg,np->baselevel,level,np->h2,np->h3,scal)!=NUM_OK) REP_ERR_RETURN (1);
      t=0.0; for (j=0; j<VD_NCOMP(x); j++) t += np->weight[j]*scal[j];
      for (j=0; j<VD_NCOMP(x); j++) alpha0[j] = -t/s;
      if (s_ddot (theNP->base.mg,np->baselevel,level,np->h2,np->t,scal)!=NUM_OK) REP_ERR_RETURN (1);
      t=0.0; for (j=0; j<VD_NCOMP(x); j++) t += np->weight[j]*scal[j];
      for (j=0; j<VD_NCOMP(x); j++) alpha1[j] = -t/np->sp;
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->t,np->h3)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,np->pp)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->pp,np->p)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->p,np->h3)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,alpha0,np->pp)!= NUM_OK) REP_ERR_RETURN (1);
      if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,alpha1,np->h1)!= NUM_OK) REP_ERR_RETURN (1);
      np->sp = s;
    }

    /* redisuum */
    if (LinearResiduum(theNP,bl,level,x,b,A,lresult)) REP_ERR_RETURN (1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  FreeVD(theNP->base.mg,bl,level,np->h1);
  FreeVD(theNP->base.mg,bl,level,np->h2);
  FreeVD(theNP->base.mg,bl,level,np->h3);
  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
  }

  return (0);
}

static INT CRConstruct (NP_BASE *theNP)
{
  NP_CR *np;

  theNP->Init = CRInit;
  theNP->Display = CRDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_CR *) theNP;
  np->ls.PreProcess = CRPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = CRSolver;
  np->ls.PostProcess = CRPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitLinearSolver	- Init this file

   SYNOPSIS:
   INT InitLinearSolver ();

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

INT InitLinearSolver ()
{
  INT i;

  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".ls",sizeof(NP_LS),LSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".cg",sizeof(NP_CG),CGConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".cr",sizeof(NP_CR),CRConstruct))
    REP_ERR_RETURN (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;
  if (MakeStruct(":ls")) REP_ERR_RETURN(__LINE__);
  if (MakeStruct(":ls:avg")) REP_ERR_RETURN(__LINE__);

  return (0);
}
