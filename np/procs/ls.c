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

#include "general.h"
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

typedef struct
{
  NP_LINEAR_SOLVER ls;

  NP_ITER *Iter;
  INT maxiter;
  INT baselevel;
  INT display;
  VECDATA_DESC *c;

} NP_LS;

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
      np->abslimit[i] = 1.0;
  if (sc_read(np->reduction,NULL,"red",argc,argv))
    return(NP_ACTIVE);

  if ((np->x == NULL) || (np->b == NULL) || (np->A == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPLinearSolverDisplay (NP_LINEAR_SOLVER *np)
{
  UserWrite("symbolic user data:\n");
  if (np->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  if (np->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
  if (np->b != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->x));
  UserWrite("\n");

  UserWrite("configuration parameters:\n");
  if (np->x != NULL)
    if (sc_disp(np->reduction,np->x,"red"))
      return (1);
  if (sc_disp(np->abslimit,np->x,"abslimit"))
    return (1);

  return(0);
}

INT NPLinearSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LINEAR_SOLVER *np;
  VECDATA_DESC *x,*b;
  MATDATA_DESC *A;
  LRESULT lresult;
  INT result,level,bl;

  np = (NP_LINEAR_SOLVER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no vector x");
    return (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no vector b");
    return (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no matrix A");
    return (1);
  }

  if (ReadOption("i",argc,argv)) {
    if (*np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,np->x,np->b,np->A,&bl,&result)) {
      UserWriteF("NPLinearSolverExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("d",argc,argv)) {
    if (*np->Defect == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Defect");
      return (1);
    }
    if ((*np->Defect)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPLinearSolverExecute: Defect failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("r",argc,argv)) {
    if (*np->Residuum == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Residuum");
      return (1);
    }
    if ((*np->Residuum)(np,level,np->x,np->b,np->A,&lresult)) {
      UserWriteF("NPLinearSolverExecute: Residuum failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("s",argc,argv)) {
    if (*np->Solver == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Solver");
      return (1);
    }
    if ((*np->Solver)(np,level,np->x,np->b,np->A,
                      np->reduction,np->abslimit,&lresult)) {
      UserWriteF("NPLinearSolverExecute: Solver failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("p",argc,argv)) {
    if (*np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPLinearSolverExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }
  return(0);
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
   npcreate sm $t ilu;
    npinit $n 1 $damp 1.0 $beta 0.0;
   npcreate base $t ls;
    npinit $m 10 $i ilu $d no $red 1e-4;
   npcreate mgc $t lmgc;
    npinit $g 1 $s sm sm base;
   npcreate solver $t ls;
    npinit $x x $b b $A MAT $m 20 $I mgc $d full $red 1e-4;

   npexecute solver $i $d $r $s;
   npexecute solver $s;
   npexecute solver $p;
   .ve
   D*/
/****************************************************************************/

static INT LinearSolverInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    return(NP_NOT_ACTIVE);

  np->display = ReadArgvDisplay(argc,argv);
  np->Iter = (NP_ITER *)ReadArgvNumProc(theNP->mg,"I",argc,argv);
  if (np->Iter == NULL)
    return(NP_NOT_ACTIVE);
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
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->c));

  return (0);
}

static INT LinearSolverPreProcess (NP_LINEAR_SOLVER *theNP, INT level,
                                   VECDATA_DESC *x, VECDATA_DESC *b,
                                   MATDATA_DESC *A,
                                   INT *baselevel, INT *result)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

  if (np->Iter->PreProcess != NULL)
    if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result))
      return(1);
  np->baselevel = *baselevel;

  return(0);
}

static INT LinearDefect (NP_LINEAR_SOLVER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *result)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

  if (s_dmatmul_minus(theNP->base.mg,np->baselevel,level,b,A,x,EVERY_CLASS)
      != NUM_OK) {
    result[0] = __LINE__;
    return(1);
  }
}

static INT LinearResiduum (NP_LINEAR_SOLVER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           LRESULT *lresult)
{
  NP_LS *np;
  INT bl;

  np = (NP_LS *) theNP;
  bl = np->baselevel;

        #ifdef ModelP
  if (a_vector_collect(theNP->base.mg,bl,level,b)) {
    lresult->error_code = __LINE__;
    return(1);
  }
        #endif

  if (s_eunorm(theNP->base.mg,bl,level,b,lresult->last_defect)) {
    lresult->error_code = __LINE__;
    return(1);
  }

  return(0);
}

static INT LinearSolver (NP_LINEAR_SOLVER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         VEC_SCALAR abslimit, VEC_SCALAR reduction,
                         LRESULT *lresult)
{
  NP_LS *np;
  VEC_SCALAR defect2reach;
  INT i,bl,PrintID;
  char text[DISPLAY_WIDTH+4];

  np = (NP_LS *) theNP;
  bl = np->baselevel;
  if (np->Iter->Iter == NULL) {
    lresult->error_code = __LINE__;
    return(1);
  }
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->c)) {
    lresult->error_code = __LINE__;
    return(1);
  }
  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH," Linear Solver ",'*',"\n");
  if (PreparePCR(x,np->display,text,&PrintID)) {
    lresult->error_code = __LINE__;
    return(1);
  }
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) {
    lresult->error_code = __LINE__;
    return(1);
  }
  if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) {
    lresult->error_code = __LINE__;
    return(1);
  }
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;
    if (l_dset(GRID_ON_LEVEL(theNP->base.mg,level),np->c,EVERY_CLASS,0.0)
        != NUM_OK) {
      lresult->error_code = __LINE__;
      return(1);
    }
    if ((*np->Iter->Iter)(np->Iter,level,np->c,b,A,&lresult->error_code))
      return (1);
    if (a_daxpy (theNP->base.mg,bl,level,x,EVERY_CLASS,Factor_One,np->c)
        != NUM_OK) {
      lresult->error_code = __LINE__;
      return(1);
    }
    if (LinearResiduum(theNP,level,x,b,A,lresult))
      return(1);
    if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) {
      lresult->error_code = __LINE__;
      return (1);
    }
    if (sc_cmp(lresult->last_defect,abslimit,b) ||
        sc_cmp(lresult->last_defect,defect2reach,b)) {
      lresult->converged = 1;
      break;
    }
  }
  FreeVD(theNP->base.mg,np->c);
  if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) {
    lresult->error_code = __LINE__;
    return (1);
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

static INT LSConstruct (NP_BASE *theNP)
{
  NP_LINEAR_SOLVER *np;

  theNP->Init = LinearSolverInit;
  theNP->Display = LinearSolverDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_LINEAR_SOLVER *) theNP;
  np->PreProcess = LinearSolverPreProcess;
  np->Defect = LinearDefect;
  np->Residuum = LinearResiduum;
  np->Solver = LinearSolver;
  np->PostProcess = LinearSolverPostProcess;

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

  if (CreateClass ("ls", sizeof(NP_LS), LSConstruct))
    return (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  return (0);
}
