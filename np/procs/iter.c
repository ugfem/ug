// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  iter.c			                                                                                */
/*																			*/
/* Purpose:   iteration num procs                                               */
/*																			*/
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
#include "np.h"

#include "transfer.h"
#include "ls.h"
#include "iter.h"

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
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  NP_ITER iter;

  VEC_SCALAR damp;
  VEC_SCALAR beta;
  INT n;

  MATDATA_DESC *L;
  VECDATA_DESC *t;

} NP_ILU;

typedef struct
{
  NP_ITER iter;

  INT gamma;
  INT nu1;
  INT nu2;
  INT baselevel;

  NP_TRANSFER *Transfer;
  NP_ITER *PreSmooth;
  NP_ITER *PostSmooth;
  NP_LINEAR_SOLVER *BaseSolver;

  VECDATA_DESC *t;

} NP_LMGC;

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
   NP_ITER - type definition for iterations

   DESCRIPTION:
   This numproc type is used for the description of linear iterations.
   It can be called by the given interface from a linear solver.
   Initializing the data is optional; it can be done with

   'INT NPIterInit (NP_ITER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPIterDisplay (NP_ITER *theNP);'
   'INT NPIterExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb


   ..... fill in data structure here when the realizition is finished


   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPIterInit (NP_ITER *np, INT argc , char **argv)
{
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->c = ReadArgvVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);

  if ((np->A == NULL) || (np->b == NULL) || (np->c == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPIterDisplay (NP_ITER *np)
{
  if ((np->A == NULL) && (np->b == NULL) && (np->c == NULL))
    return(0);
  UserWrite("symbolic user data:\n");
  if (np->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  if (np->b != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->b));
  if (np->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->c));
  UserWrite("\n");

  return(0);
}

INT NPIterExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ITER *np;
  INT result,bl,level;

  np = (NP_ITER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->c == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no vector c");
    return (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no vector b");
    return (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no matrix A");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,np->c,np->b,np->A,&bl,&result)) {
      UserWriteF("NPIterExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->Iter == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no Iter");
      return (1);
    }
    if ((*np->Iter)(np,level,np->c,np->b,np->A,&result)) {
      UserWriteF("NPIterExecute: Iter failed, error code %d\n", result);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->c,np->b,np->A,&result)) {
      UserWriteF("NPIterExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  return(0);
}

/****************************************************************************/
/*D
   ilu - numproc for point block beta-modified ilu smoother

   DESCRIPTION:
   This numproc executes a point block ilu smoother, using the blas routines
   'l_ilubthdecomp' and 'l_luiter'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list> $beta <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate sm $t ilu;
    npinit $n 1 $damp 1.0 $beta 0.0;
   npcreate base $t ls;
    npinit $m 10 $i ilu $d no $red 1e-4;
   npcreate mgc $t lmgc;
    npinit $g 1 $S sm sm base;
   npcreate solver $t ls;
    npinit $x x $b b $A MAT $m 20 $i mgc $d full $red 1e-4;

   npexecute solver $i $d $r $s $p;
   .ve
   D*/
/****************************************************************************/

static INT ILU_Init (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ILU *np;
  INT i;

  np = (NP_ILU *) theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i] = 1.0;
  sc_read(np->damp,np->iter.b,"damp",argc,argv);
  for (i=0; i<MAX_VEC_COMP; i++) np->beta[i] = 0.0;
  sc_read(np->damp,np->iter.b,"beta",argc,argv);
  if (ReadArgvINT("n",&(np->n),argc,argv))
    np->n = 1;
  np->L = ReadArgvMatDesc(theNP->mg,"L",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);

  return (NPIterInit(&np->iter,argc,argv));
}

static INT ILU_Display (NP_BASE *theNP)
{
  NP_ILU *np;

  np = (NP_ILU *) theNP;

  NPIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) return (1);
  if (sc_disp(np->damp,np->iter.b,"beta")) return (1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",(int)np->n);
  if (np->L != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(np->L));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));

  return (0);
}

static INT ILU_PreProcess  (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *baselevel, INT *result)
{
  NP_ILU *np;
  GRID *theGrid;

  np = (NP_ILU *) theNP;

  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->L)) {
    result[0] = __LINE__;
    return (1);
  }
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_setindex(theGrid)) {
    result[0] = __LINE__;
    return (1);
  }
  if (l_dmatcopy(theGrid,np->L,A) != NUM_OK) {
    result[0] = __LINE__;
    return (1);
  }
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->L,TRUE) != NUM_OK) {
    result[0] = __LINE__;
    return (1);
  }
        #endif
  if (l_ilubthdecomp(theGrid,np->L,np->beta,NULL,NULL,NULL) !=NUM_OK) {
    result[0] = __LINE__;
    return (1);
  }
  *baselevel = level;

  return (0);
}

static INT ILU_Iter (NP_ITER *theNP, INT level,
                     VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                     INT *result)
{
  NP_ILU *np;
  GRID *theGrid;
  INT i;

  np = (NP_ILU *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);

  if (AllocVDFromVD(theNP->base.mg,level,level,x,&np->t)) {
    result[0] = __LINE__;
    return(1);
  }

  for (i=0; i<np->n; i++) {
        #ifdef ModelP
    if (l_vector_collect(theGrid,b) != NUM_OK) {
      result[0] = __LINE__;
      return (1);
    }
            #endif
    if (l_luiter(theGrid,np->t,np->L,b) != NUM_OK) {
      result[0] = __LINE__;
      return (1);
    }
        #ifdef ModelP
    if (l_vector_consistent(theGrid,np->t) != NUM_OK) {
      result[0] = __LINE__;
      return (1);
    }
        #endif
    if (l_luiter(theGrid,np->t,np->L,b) != NUM_OK) {
      result[0] = __LINE__;
      return (1);
    }
    if (l_dscale(theGrid,np->t,ACTIVE_CLASS,np->damp) != NUM_OK) {
      result[0] = __LINE__;
      return (1);
    }
    if (l_daxpy(theGrid,x,ACTIVE_CLASS,Factor_One,np->t) != NUM_OK) {
      result[0] = __LINE__;
      return (1);
    }
    if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,np->t,ACTIVE_CLASS)
        != NUM_OK) {
      result[0] = __LINE__;
      return (1);
    }
  }
  FreeVD(theNP->base.mg,level,level,np->t);

  return (0);
}

static INT ILU_PostProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *result)
{
  NP_ILU *np;

  np = (NP_ILU *) theNP;

  FreeMD(theNP->base.mg,level,level,np->L);

  return(0);
}

static INT ILUConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = ILU_Init;
  theNP->Display = ILU_Display;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = ILU_PreProcess;
  np->Iter = ILU_Iter;
  np->PostProcess = ILU_PostProcess;

  return(0);
}

/****************************************************************************/
/*D
   lmgc - numproc for linear multigrid cycle

   DESCRIPTION:
   This numproc executes

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $S <pre post base> $T <transfer>
       $b <baselevel> $g <gamma> [$n1 <it>] [$n2 <it>]
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $T~<transfer> - transfer numproc
   .  $S~<pre~post~base> - numprocs for pre- and postsmoother, base solver
   .  $b~<baselevel> - baselevel where the base solver is called
   .  $g~<gamma> - number of iterations of Lmgc per level
   .  $n1~<it> - number of iterations of the presmoother
   .  $n2~<it> - number of iteration of the postsmoother

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate sm $t ilu;
    npinit $n 1 $damp 1.0 $beta 0.0;
   npcreate base $t ls;
    npinit $m 10 $i ilu $d no $red 1e-4;
   npcreate mgc $t lmgc;
    npinit $g 1 $S sm sm base;
   npcreate solver $t ls;
    npinit $x x $b b $A MAT $m 20 $i mgc $d full $red 1e-4;

   npexecute solver $i $d $r $s $p;
   .ve
   D*/
/****************************************************************************/

static INT LmgcInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LMGC *np;
  INT i;
  char post[VALUELEN],pre[VALUELEN],base[VALUELEN];

  np = (NP_LMGC *) theNP;

  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->Transfer = (NP_TRANSFER *)ReadArgvNumProc(theNP->mg,"T",argc,argv);
  for (i=1; i<argc; i++)
    if (argv[i][0]=='s') {
      if (sscanf(argv[i],"s %s %s %s",pre,post,base)!=3)
        continue;
      np->PreSmooth = (NP_ITER *)GetNumProcFromName(pre);
      np->PostSmooth = (NP_ITER *)GetNumProcFromName(post);
      np->BaseSolver = (NP_LINEAR_SOLVER *)GetNumProcFromName(base);
      break;
    }

  if (ReadArgvINT("g",&(np->gamma),argc,argv))
    np->gamma = 1;
  if (ReadArgvINT("n1",&(np->nu1),argc,argv))
    np->nu1 = 1;
  if (ReadArgvINT("n2",&(np->nu2),argc,argv))
    np->nu2 = 1;
  if (ReadArgvINT("b",&(np->baselevel),argc,argv))
    np->baselevel = 0;

  if (np->Transfer == NULL) return(NP_NOT_ACTIVE);
  if (np->PreSmooth == NULL) return(NP_NOT_ACTIVE);
  if (np->PostSmooth == NULL) return(NP_NOT_ACTIVE);
  if (np->BaseSolver == NULL) return(NP_NOT_ACTIVE);

  return (NPIterInit(&np->iter,argc,argv));
}

static INT LmgcDisplay (NP_BASE *theNP)
{
  NP_LMGC *np;

  np = (NP_LMGC *) theNP;

  NPIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->gamma);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)np->nu1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n2",(int)np->nu2);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"b",(int)np->baselevel);

  if (np->Transfer != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(np->Transfer));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (np->PreSmooth != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"pre",ENVITEM_NAME(np->PreSmooth));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"pre","---");
  if (np->PostSmooth != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"post",ENVITEM_NAME(np->PostSmooth));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"post","---");
  if (np->BaseSolver != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"base",ENVITEM_NAME(np->BaseSolver));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"base","---");

  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));

  return (0);
}

static INT LmgcPreProcess  (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *baselevel, INT *result)
{
  NP_LMGC *np;
  INT i;

  np = (NP_LMGC *) theNP;

  if (np->Transfer->PreProcess != NULL)
    for (i = np->baselevel; i <= level; i++)
      if ((*np->Transfer->PreProcess)
            (np->Transfer,i,x,b,A,baselevel,result))
        return(1);

  if (np->PreSmooth->PreProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PreProcess)
            (np->PreSmooth,i,x,b,A,baselevel,result))
        return(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PreProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PreSmooth->PreProcess)
              (np->PostSmooth,i,x,b,A,baselevel,result))
          return(1);

  if (np->BaseSolver->PreProcess != NULL)
    if ((*np->BaseSolver->PreProcess)
          (np->BaseSolver,np->baselevel,x,b,A,baselevel,result))
      return(1);
  *baselevel = np->baselevel;

  return (0);
}

static INT Lmgc (NP_ITER *theNP, INT level,
                 VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A,
                 INT *result)
{
  NP_LMGC *np;
  MULTIGRID *theMG;
  LRESULT lresult;
  INT i;

  np = (NP_LMGC *) theNP;

  if (level <= np->baselevel) {
    if ((*np->BaseSolver->Residuum)
          (np->BaseSolver,level,c,b,A,&lresult))
      return(1);
    if ((*np->BaseSolver->Solver)(np->BaseSolver,level,c,b,A,
                                  np->BaseSolver->abslimit,
                                  np->BaseSolver->reduction,&lresult)) {
      result[0] = lresult.error_code;
      return(1);
    }
    if (!lresult.converged)
      PrintErrorMessage('W',"Lmgc","no convergence of BaseSolver");
    return(0);
  }

  theMG = theNP->base.mg;
  for (i=0; i<np->nu1; i++)
    if ((*np->PreSmooth->Iter)(np->PreSmooth,level,c,b,A,result))
      return(1);
  if ((*np->Transfer->RestrictDefect)
        (np->Transfer,level,b,b,A,Factor_One,result))
    return(1);
  if (l_dset(GRID_ON_LEVEL(theMG,level-1),c,EVERY_CLASS,0.0) != NUM_OK) {
    result[0] = __LINE__;
    return(1);
  }
  for (i=0; i<np->gamma; i++)
    if (Lmgc(theNP,level-1,c,b,A,result))
      return(1);
  if (AllocVDFromVD(theMG,level,level,c,&np->t)) {
    result[0] = __LINE__;
    return(1);
  }
  if ((*np->Transfer->InterpolateCorrection)
        (np->Transfer,level,np->t,c,A,Factor_One,result))
    return(1);
  if (l_daxpy  (GRID_ON_LEVEL(theMG,level),c,EVERY_CLASS,Factor_One,np->t)
      != NUM_OK) {
    result[0] = __LINE__;
    return(1);
  }
  if (l_dmatmul_minus(GRID_ON_LEVEL(theMG,level),b,2,A,np->t,EVERY_CLASS)
      != NUM_OK) {
    result[0] = __LINE__;
    return(1);
  }
  FreeVD(theMG,level,level,np->t);
  for (i=0; i<np->nu1; i++)
    if ((*np->PostSmooth->Iter)(np->PostSmooth,level,c,b,A,result))
      return(1);
  if (np->Transfer->AdaptCorrection != NULL)
    if ((*np->Transfer->AdaptCorrection)(np->Transfer,level,c,b,A,result))
      return(1);

  return (0);
}

static INT LmgcPostProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *result)
{
  NP_LMGC *np;
  INT i;

  np = (NP_LMGC *) theNP;

  if (np->Transfer->PostProcess != NULL)
    for (i = np->baselevel; i <= level; i++)
      if ((*np->Transfer->PostProcess)(np->Transfer,i,x,b,A,result))
        return(1);

  if (np->PreSmooth->PostProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PostProcess)
            (np->PreSmooth,i,x,b,A,result))
        return(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PostProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PreSmooth->PostProcess)
              (np->PostSmooth,i,x,b,A,result))
          return(1);

  if (np->BaseSolver->PostProcess != NULL)
    if ((*np->BaseSolver->PostProcess)
          (np->BaseSolver,np->baselevel,x,b,A,result))
      return(1);

  return (0);
}

static INT LmgcConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = LmgcInit;
  theNP->Display = LmgcDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = LmgcPreProcess;
  np->Iter = Lmgc;
  np->PostProcess = LmgcPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitIter	- Init this file

   SYNOPSIS:
   INT InitIter ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

INT InitIter ()
{
  INT i;

  if (CreateClass ("ilu", sizeof(NP_ILU), ILUConstruct))
    return (__LINE__);
  if (CreateClass ("lmgc", sizeof(NP_LMGC), LmgcConstruct))
    return (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  return (0);
}
