// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  eiter.c			                                                                                */
/*																			*/
/* Purpose:   extended iter num proc type                                   */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            IWR/Technische Simulation                                     */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: klaus.johannsen@iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   22.07.02 begin                                                */
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
#include <math.h>

#include "general.h"
#include "debug.h"
#include "dlmgr.h"
#include "gm.h"
#include "ugm.h"
#include "algebra.h"
#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "ugdevices.h"
#include "udm.h"
#include "pcr.h"
#include "debug.h"
#include "fifo.h"
#include "evm.h"
#include "misc.h"
#include "ugstruct.h"

#include "transfer.h"
#include "ls.h"
#include "els.h"
#include "iter.h"
#include "eiter.h"
#include "project.h"
#include "disctools.h"
#include "block.h"

#include "ff_gen.h"
#include "ff.h"
#include "ugeblas.h"
#include "order.h"

#ifdef __cplusplus
#ifdef __TWODIM__
using namespace UG2d;
#else
using namespace UG3d;
#endif
#endif

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

  NP_EITER eiter;

  VECDATA_DESC *Af[EXTENSION_MAX];
  DOUBLE sc[EXTENSION_MAX*EXTENSION_MAX];
  NP_ITER *iter;

} NP_SCITER;

typedef struct
{
  NP_EITER iter;

  INT gamma;
  INT nu1;
  INT nu2;
  INT baselevel;

  NP_TRANSFER *Transfer;
  NP_EITER *PreSmooth;
  NP_EITER *PostSmooth;
  NP_ELINEAR_SOLVER *BaseSolver;

  EVECDATA_DESC *t;

  EVEC_SCALAR damp;

} NP_ELMGC;

typedef struct
{
  NP_EITER eiter;

  INT n;                            /* # vectors                            */
  INT MarkKey[MAXLEVEL];            /* key for Mark/Release                 */
  INT count;                        /* counter for MarkKey                  */
  DOUBLE *mat[MAXLEVEL];            /* matrix                               */
  DOUBLE *scale[MAXLEVEL];          /* diagonal scaling of matrix           */
  INT mem;                          /* memory used temporary (bytes)        */

  DOUBLE *sol;                      /* sol vector                           */
  DOUBLE *rhs;                      /* rhs vector                           */

} NP_EEX;

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

INT NS_PREFIX NPEIterInit (NP_EITER *np, INT argc , char **argv)
{
  np->A = ReadArgvEMatDesc(np->base.mg,"A",argc,argv);
  np->c = ReadArgvEVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvEVecDesc(np->base.mg,"r",argc,argv);

  if ((np->A == NULL) || (np->b == NULL) || (np->c == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NS_PREFIX NPEIterDisplay (NP_EITER *np)
{
  if ((np->A == NULL) && (np->b == NULL) && (np->c == NULL))
    return(0);
  UserWrite("symbolic user data:\n");
  if (np->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  if (np->b != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->b));
  if (np->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->c));
  UserWrite("\n");

  return(0);
}


INT NS_PREFIX NPEIterExecute (NP_BASE *theNP, INT argc , char **argv)
{
  REP_ERR_RETURN (1);
}

/****************************************************************************/
/*D
   sciter - numproc for schur complement iterations

   DESCRIPTION:
   This numproc realizes the schur complement elimination

   .vb
   npinit <name>
   .ve

   .  $d~<dummy> - correction vector

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT SCITER_Init (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SCITER *np = (NP_SCITER *) theNP;

  np->iter=(NP_ITER *)ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);
  if (np->iter==NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);

  return (NPEIterInit(&np->eiter,argc,argv));
}

static INT SCITER_Display (NP_BASE *theNP)
{
  NP_SCITER *np = (NP_SCITER *) theNP;

  NPEIterDisplay(&np->eiter);

  if (np->iter != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"iter",ENVITEM_NAME(np->iter));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"iter","---");

  return (0);
}

static INT SCITER_PreProcess  (NP_EITER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *baselevel, INT *result)
{
  INT i,j;
  DOUBLE a;
  NP_SCITER *np = (NP_SCITER *) theNP;
  VECDATA_DESC *tmp=NULL;

  /* prepare internal iteration */
  if (np->iter->PreProcess!=NULL)
    if ((*np->iter->PreProcess)(np->iter,level,x->vd,b->vd,A->mm,baselevel,result)) REP_ERR_RETURN(1);

  /* allocate Af-vectors */
  if (AllocVDFromVD(NP_MG(theNP),level,level,x->vd,&tmp)) NP_RETURN(1,result[0]);
  for (i=0; i<x->n; i++)
    if (AllocVDFromVD(NP_MG(theNP),level,level,x->vd,&(np->Af[i]))) NP_RETURN(1,result[0]);

  /* calculate Af and schur complement */
  for (i=0; i<x->n; i++)
  {
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,tmp,A->me[i])!= NUM_OK) REP_ERR_RETURN (1);
    if ((*np->iter->Iter)(np->iter,level,np->Af[i],tmp,A->mm,result)) NP_RETURN(1,result[0]);
    for (j=0; j<x->n; j++)
    {
      if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,A->em[j],np->Af[i],&a)!=NUM_OK) REP_ERR_RETURN (1);
      np->sc[j*x->n+i]=EMDD_EE(A,level,j*x->n+i)-a;
    }
  }

  /* free one */
  if (FreeVD(NP_MG(theNP),level,level,tmp)) NP_RETURN(1,result[0]);

  return (0);
}

static INT SolveFullMatrix_copy (INT n, DOUBLE *x, DOUBLE *mat, DOUBLE *b)
{
  DOUBLE copy[EXTENSION_MAX*EXTENSION_MAX];
  INT i,j;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      copy[n*i+j]=mat[n*i+j];
  return (SolveFullMatrix(n,x,copy,b));
}

static INT SCITER_Iter (NP_EITER *theNP, INT level, EVECDATA_DESC *c, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  INT i;
  DOUBLE s[EXTENSION_MAX];
  NP_SCITER *np = (NP_SCITER *) theNP;
  VECDATA_DESC *tmp=NULL;
  DOUBLE norm;

  /* inner iteration */
  if (AllocVDFromVD(NP_MG(theNP),level,level,c->vd,&tmp)) NP_RETURN(1,result[0]);
  if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,tmp,b->vd)!= NUM_OK) REP_ERR_RETURN (1);
  if ((*np->iter->Iter)(np->iter,level,c->vd,tmp,A->mm,result)) NP_RETURN(1,result[0]);
  if (FreeVD(NP_MG(theNP),level,level,tmp)) NP_RETURN(1,result[0]);

  /* outter iteration */
  for (i=0; i<c->n; i++)
  {
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,A->em[i],c->vd,s+i)!=NUM_OK) REP_ERR_RETURN (1);
    s[i]=EVDD_E(b,level,i)-s[i];
  }
  if (SolveFullMatrix_copy(c->n,EVDD_E_PTR(c,level),np->sc,s)) NP_RETURN(1,result[0]);
  for (i=0; i<c->n; i++)
    if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,c->vd,-EVDD_E(c,level,i),np->Af[i])) REP_ERR_RETURN (1);

  /* update defect */
  if (dematmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,c)!= NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SCITER_PostProcess (NP_EITER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  INT i;
  NP_SCITER *np = (NP_SCITER *) theNP;

  /* free Af-vectors */
  for (i=0; i<x->n; i++)
    if (FreeVD(NP_MG(theNP),level,level,np->Af[i])) NP_RETURN(1,result[0]);

  /* postprocess internal iteration */
  if (np->iter->PostProcess!=NULL)
    if ((*np->iter->PostProcess)(np->iter,level,x->vd,b->vd,A->mm,result)) REP_ERR_RETURN(1);

  return (0);
}

static INT SCITER_Construct (NP_BASE *theNP)
{
  NP_EITER *np = (NP_EITER *) theNP;

  theNP->Init             = SCITER_Init;
  theNP->Display          = SCITER_Display;
  theNP->Execute          = NPIterExecute;
  np->PreProcess          = SCITER_PreProcess;
  np->Iter                        = SCITER_Iter;
  np->PostProcess         = SCITER_PostProcess;

  return(0);
}

/****************************************************************************/
/*D
   elmgc - numproc for extended linear multigrid cycle

   DESCRIPTION:
   This numproc executes an extended linear multigrid cycle.

   .vb
   npinit <name> [$c <cor>] [$r <rhs>] [$A <mat>]
       $S <pre post base> $T <transfer>
       [$b <baselevel>] [$g <gamma>] [$n1 <it>] [$n2 <it>]
   .ve

   .  $c~<cor> - correction vector
   .  $r~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $T~<transfer> - transfer numproc
   .  $S~<pre~post~base> - numprocs for pre- and postsmoother, base solver
   .  $b~<baselevel> - baselevel where the base solver is called
   .  $g~<gamma> - number of iterations of Lmgc per level (default gamma = 1)
   .  $n1~<it> - number of iterations of the presmoother (default n1 = 1)
   .  $n2~<it> - number of iteration of the postsmoother (default n2 = 1)

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT ELmgcInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ELMGC *np;
  INT i,ret;
  char post[VALUELEN],pre[VALUELEN],base[VALUELEN];

  np = (NP_ELMGC *) theNP;

  np->t = ReadArgvEVecDesc(theNP->mg,"t",argc,argv);
  np->Transfer = (NP_TRANSFER *)ReadArgvNumProc(theNP->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  for (i=1; i<argc; i++)
    if (argv[i][0]=='S')
    {
      if (sscanf(argv[i],"S %s %s %s",pre,post,base)!=3) continue;
      np->PreSmooth = (NP_EITER *)GetNumProcByName(theNP->mg,pre,EITER_CLASS_NAME);
      np->PostSmooth = (NP_EITER *)GetNumProcByName(theNP->mg,post,EITER_CLASS_NAME);
      np->BaseSolver = (NP_ELINEAR_SOLVER *)GetNumProcByName(theNP->mg,base,ELINEAR_SOLVER_CLASS_NAME);
      break;
    }

  if (ReadArgvINT("g",&(np->gamma),argc,argv)) np->gamma = 1;
  if (ReadArgvINT("n1",&(np->nu1),argc,argv)) np->nu1 = 1;
  if (ReadArgvINT("n2",&(np->nu2),argc,argv)) np->nu2 = 1;
  if (ReadArgvINT("b",&(np->baselevel),argc,argv)) np->baselevel = 0;
  if (np->baselevel<0)
  {
    for (i=FULLREFINELEVEL(NP_MG(theNP)); i>0; i--)
      if (NVEC(GRID_ON_LEVEL(NP_MG(theNP),i))<=-np->baselevel)
        break;
    np->baselevel=i;
  }

  if (np->Transfer == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PreSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PostSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->BaseSolver == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);

  ret = NPEIterInit(&np->iter,argc,argv);

  if (esc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->damp[i] = 1.0;

  return (ret);
}

static INT ELmgcDisplay (NP_BASE *theNP)
{
  NP_ELMGC *np;

  np = (NP_ELMGC *) theNP;

  NPEIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->gamma);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)np->nu1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n2",(int)np->nu2);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);

  if (np->Transfer != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(np->Transfer));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (np->PreSmooth != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"pre",ENVITEM_NAME(np->PreSmooth));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"pre","---");
  if (np->PostSmooth != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"post",ENVITEM_NAME(np->PostSmooth));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"post","---");
  if (np->BaseSolver != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"base",ENVITEM_NAME(np->BaseSolver));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"base","---");

  if (np->t != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));

  if (np->iter.b!=NULL)           { if (esc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);}
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"damp","--- (cannot display)");

  return (0);
}

static INT ELmgcPreProcess (NP_EITER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_ELMGC *np = (NP_ELMGC *) theNP;
  INT i;

  if (np->Transfer->PreProcess != NULL)
    if ((*np->Transfer->PreProcess)(np->Transfer,&(np->baselevel),level,x->vd,b->vd,A->mm,result)) REP_ERR_RETURN(1);

  if (np->PreSmooth->PreProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PreProcess)(np->PreSmooth,i,x,b,A,baselevel,result)) REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PreProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PostSmooth->PreProcess)(np->PostSmooth,i,x,b,A,baselevel,result)) REP_ERR_RETURN(1);

  *baselevel = MIN(np->baselevel,level);
  if (np->gamma>0 && np->BaseSolver->PreProcess != NULL)
    if ((*np->BaseSolver->PreProcess)(np->BaseSolver,*baselevel,x,b,A,baselevel,result)) REP_ERR_RETURN(1);

  return (0);
}

static INT ELmgc (NP_EITER *theNP, INT level, EVECDATA_DESC *c, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  NP_ELMGC *np=(NP_ELMGC *) theNP;
  MULTIGRID *theMG=NP_MG(theNP);
  GRID *theGrid;
  ELRESULT lresult;
  INT i;
  DOUBLE eunorm;

  /* admin */
  NPEIT_A(theNP) = A;
  NPEIT_c(theNP) = c;
  NPEIT_b(theNP) = b;
  theGrid = GRID_ON_LEVEL(theMG,level);

  if (level <= np->baselevel)
  {
    /* solve on baselevel */
    if ((*np->BaseSolver->Residuum)(np->BaseSolver,MIN(level,np->baselevel),level,c,b,A,&lresult)) REP_ERR_RETURN(1);
    if ((*np->BaseSolver->Solver)(np->BaseSolver,level,c,b,A, np->BaseSolver->abslimit, np->BaseSolver->reduction,&lresult)) NP_RETURN(1,result[0]);
    return(0);
  }

  /* presmooth */
  if (AllocEVDFromEVD(theMG,level,level,c,&np->t)) NP_RETURN(1,result[0]);
  for (i=0; i<np->nu1; i++)
  {
    if ((*np->PreSmooth->Iter)(np->PreSmooth,level,np->t,b,A,result)) REP_ERR_RETURN(1);
    if (deadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }

  /* coarse grid correction */
  if ((*np->Transfer->RestrictDefect)(np->Transfer,level,b->vd,b->vd,A->mm,Factor_One,result)) REP_ERR_RETURN(1);
  EVDD_E(b,level-1,0)=EVDD_E(b,level,0);
  if (deset(theMG,level-1,level-1,ALL_VECTORS,c,0.0) != NUM_OK) NP_RETURN(1,result[0]);
  for (i=0; i<np->gamma; i++)
    if (ELmgc(theNP,level-1,c,b,A,result)) REP_ERR_RETURN(1);

  if ((*np->Transfer->InterpolateCorrection)(np->Transfer,level,np->t->vd,c->vd,A->mm,np->damp,result)) REP_ERR_RETURN(1);
  EVDD_E(b,level,0)=EVDD_E(b,level-1,0);
  if (deadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  if (dematmul_minus(theMG,level,level,ALL_VECTORS,b,A,np->t) != NUM_OK) NP_RETURN(1,result[0]);

  /* postsmooth */
  for (i=0; i<np->nu2; i++)
  {
    if ((*np->PostSmooth->Iter)(np->PostSmooth,level,np->t,b,A,result)) REP_ERR_RETURN(1);
    if (deadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }

  /* admin */
  if (FreeEVD(NP_MG(theNP),level,level,np->t)) REP_ERR_RETURN(1);

  return (0);
}

static INT ELmgcPostProcess (NP_EITER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  NP_ELMGC *np;
  INT i;

  np = (NP_ELMGC *) theNP;

  if (np->gamma>0 && np->BaseSolver->PostProcess != NULL)
    if ((*np->BaseSolver->PostProcess)(np->BaseSolver,np->baselevel,x,b,A,result)) REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PostProcess != NULL)
      for (i = level; i >= np->baselevel+1; i--)
        if ((*np->PostSmooth->PostProcess)(np->PostSmooth,i,x,b,A,result)) REP_ERR_RETURN(1);

  if (np->PreSmooth->PostProcess != NULL)
    for (i = level; i >= np->baselevel+1; i--)
      if ((*np->PreSmooth->PostProcess)(np->PreSmooth,i,x,b,A,result)) REP_ERR_RETURN(1);

  if (np->Transfer->PostProcess != NULL)
    if ((*np->Transfer->PostProcess)(np->Transfer,&(np->baselevel),level,x->vd,b->vd,A->mm,result)) REP_ERR_RETURN(1);

  return (0);
}

static INT ELmgcConstruct (NP_BASE *theNP)
{
  NP_EITER *np;

  theNP->Init = ELmgcInit;
  theNP->Display = ELmgcDisplay;
  theNP->Execute = NPEIterExecute;

  np = (NP_EITER *) theNP;
  np->PreProcess = ELmgcPreProcess;
  np->Iter = ELmgc;
  np->PostProcess = ELmgcPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   eex - numproc for extended exact solver

   DESCRIPTION:
   This numproc solves an extended linear system

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>];
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix

   'npexecute <name> [$i] [$s] [$p];'

   .  $p - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT EEXCopyMatrix (GRID *theGrid, EVECDATA_DESC *x, EMATDATA_DESC *A, INT n, DOUBLE *mat)
{
  INT ment,index,rindex,rtype,rcomp,cindex,ctype,ccomp,i,j;
  VECTOR *theV,*theW;
  MATRIX *theM;
  SHORT *comp;

    #ifdef ModelP
  if (FIRSTVECTOR(theGrid) == NULL)
    return(0);
    #endif

  for (i=0; i<n*n; i++) mat[i]=0.0;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    rindex = VINDEX(theV);
    rtype = VTYPE(theV);
    rcomp = VD_NCMPS_IN_TYPE(x->vd,rtype);
    for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
    {
      theW = MDEST(theM);
      cindex = VINDEX(theW);
      ctype = VTYPE(theW);
      ccomp = VD_NCMPS_IN_TYPE(x->vd,ctype);
      comp = MD_MCMPPTR_OF_RT_CT(A->mm,rtype,ctype);
      for (i=0; i<rcomp; i++)
        for (j=0; j<ccomp; j++)
          mat[(rindex+i)*n+(cindex+j)]=MVALUE(theM,comp[i*ccomp+j]);
    }
    for (i=0; i<A->n; i++)
    {
      for (j=0; j<VD_NCMPS_IN_TYPE(A->me[i],rtype); j++)
        mat[(rindex+j)*n+(n-A->n+i)]=VVALUE(theV,VD_CMP_OF_TYPE(A->me[i],rtype,j));
      for (j=0; j<VD_NCMPS_IN_TYPE(A->em[i],rtype); j++)
        mat[(n-A->n+i)*n+(rindex+j)]=VVALUE(theV,VD_CMP_OF_TYPE(A->em[i],rtype,j));
    }
  }
  for (i=0; i<A->n; i++)
    for (j=0; j<A->n; j++)
      mat[(n-A->n+i)*n+(n-A->n+j)]=EMDD_EE(A,GLEVEL(theGrid),A->n*i+j);

  return (0);
}


static INT EEXInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EEX *np;

  np = (NP_EEX *) theNP;
  np->n = -1;
  np->count = -1;

  return (NPEIterInit(&np->eiter,argc,argv));
}

static INT EEXDisplay (NP_BASE *theNP)
{
  NP_EEX *np;
  char name [32];

  np = (NP_EEX *) theNP;
  NPEIterDisplay(&np->eiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",(int)np->n);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"count",(int)np->count);

  return (0);
}

void WriteScilabVector (char *name, INT n, DOUBLE *vec)
{
  INT i,j;
  FILE *f=fopen(name,"w");

  for (i=0; i<n; i++)
    fprintf(f,"%e\n",vec[i]);

  fclose(f);
}

void WriteScilabMatrix (char *name, INT n, DOUBLE *mat)
{
  INT i,j;
  FILE *f=fopen(name,"w");

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      fprintf(f,"%e ",mat[i*n+j]);
    fprintf(f,"\n");
  }

  fclose(f);
}

static INT EEX_ScaleMatrix (INT n, DOUBLE *mat, DOUBLE *scale)
{
  INT i,j;
  DOUBLE s;

  for (i=0; i<n; i++)
  {
    s=0.0; for (j=0; j<n; j++) s+=mat[i*n+j]*mat[i*n+j];s=sqrt(s);
    if (s==0.0) return 1;
    s=1.0/s; scale[i]=s;
    for (j=0; j<n; j++)
      mat[n*i+j]*=s;
  }

  return 0;
}

static INT EEX_ScaleVector (INT n, DOUBLE *scale, DOUBLE *vec)
{
  INT i;

  for (i=0; i<n; i++)
    vec[i]*=scale[i];

  return 0;
}

static INT EEXPreProcess  (NP_EITER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_EEX *np=(NP_EEX *) theNP;
  GRID *theGrid = NP_GRID(theNP,level);
  HEAP *theHeap=MGHEAP(NP_MG(theNP));
  VECTOR *theV;
  MATRIX *theM;
  INT n;

  /* admin */
  for (theV=FIRSTVECTOR(theGrid),n=0; theV!=NULL; theV=SUCCVC(theV)) { VINDEX(theV)=n; n+=VD_NCMPS_IN_TYPE(x->vd,VTYPE(theV)); }
  n+=x->n; np->n=n; *baselevel=level;

  /* get memory */
  np->count++;
  if (MarkTmpMem(theHeap,&(np->MarkKey[np->count]))) REP_ERR_RETURN(1);
  if (np->count==0)
  {
    np->sol=(DOUBLE*)GetTmpMem(theHeap,np->n*sizeof(DOUBLE),np->MarkKey[np->count]);
    if (np->sol==NULL) REP_ERR_RETURN(1);
    np->rhs=(DOUBLE*)GetTmpMem(theHeap,np->n*sizeof(DOUBLE),np->MarkKey[np->count]);
    if (np->rhs==NULL) REP_ERR_RETURN(1);
    np->mat[np->count]=(DOUBLE*)GetTmpMem(theHeap,(np->n*np->n+np->n)*sizeof(DOUBLE),np->MarkKey[np->count]);
    if (np->mat==NULL) REP_ERR_RETURN(1);
    np->scale[np->count]=(DOUBLE*)GetTmpMem(theHeap,(np->n)*sizeof(DOUBLE),np->MarkKey[np->count]);
    if (np->scale==NULL) REP_ERR_RETURN(1);
  }

  /* copy matrix */
  if (EEXCopyMatrix(theGrid,x,A,n,np->mat[np->count])) REP_ERR_RETURN(1);
  if (EEX_ScaleMatrix(n,np->mat[np->count],np->scale[np->count])) REP_ERR_RETURN(1);

  /* decompose matrix */
  if (Yams(n,NULL,np->mat[np->count],NULL)) REP_ERR_RETURN(1);

  return 0;
}

static INT EEXSmoother (NP_EITER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  NP_EEX *np=(NP_EEX *)theNP;
  INT i,j,n,vent,type;
  GRID *theGrid;
  DOUBLE *sol,*rhs;
  VECTOR *theV;
  SHORT *comp;

  /* admin */
  theGrid = NP_GRID(theNP,level);
  NPEIT_A(theNP) = A;
  NPEIT_c(theNP) = x;
  NPEIT_b(theNP) = b;

  /* init */
  n=np->n;
  if (n==0) return(0);
  sol=np->sol; rhs=np->rhs;

  /* copy b to rhs */
  if (MD_IS_SCALAR(A->mm))
  {
    vent=VD_SCALCMP(b->vd);
    j=0;
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV),i++)
      if (VD_NCMPS_IN_TYPE(b->vd,VTYPE(theV))>0)
        rhs[j++]=VVALUE(theV,vent);
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV))
    {
      type = VTYPE(theV);
      comp = VD_CMPPTR_OF_TYPE(b->vd,type);
      for (j=0; j<VD_NCMPS_IN_TYPE(b->vd,type); j++)
        rhs[i++]=VVALUE(theV,comp[j]);
    }
  }
  for (i=0; i<A->n; i++)
    rhs[n-A->n+i]=EVDD_E(b,level,i);

  /* solve */
  if (EEX_ScaleVector(n,np->scale[np->count],np->rhs)) REP_ERR_RETURN(1);
  if (Yams(n,sol,np->mat[np->count],rhs)) REP_ERR_RETURN(1);

  /* copy sol to x */
  if (MD_IS_SCALAR(A->mm))
  {
    vent=VD_SCALCMP(x->vd);
    j=0;
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV),i++)
      if (VD_NCMPS_IN_TYPE(x->vd,VTYPE(theV))>0)
        VVALUE(theV,vent)=rhs[j++];
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV))
    {
      type = VTYPE(theV);
      comp = VD_CMPPTR_OF_TYPE(x->vd,type);
      for (j=0; j<VD_NCMPS_IN_TYPE(x->vd,type); j++)
        VVALUE(theV,comp[j])=sol[i++];
    }
  }
  for (i=0; i<A->n; i++)
    EVDD_E(x,level,i)=sol[n-A->n+i];

  /* update defect */
    #ifdef ModelP
  if (l_vector_consistent(theGrid,x->vd) != NUM_OK) NP_RETURN(1,result[0]);
    #endif

  if (dematmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK) NP_RETURN(1,result[0]);

  return 0;
}

static INT EEXPostProcess (NP_EITER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  NP_EEX *np=(NP_EEX*)theNP;
  HEAP *theHeap=MGHEAP(NP_MG(theNP));

    #ifdef ModelP
  if (FIRSTVECTOR(NP_GRID(theNP,level)) == NULL)
    return(0);
    #endif

  ReleaseTmpMem(theHeap,np->MarkKey[np->count]);
  np->mat[np->count]=NULL;
  if (np->count == 0) { np->sol=np->rhs=NULL; }
  np->count--;

  return 0;
}

static INT EEXConstruct (NP_BASE *theNP)
{
  NP_EITER *np;

  theNP->Init     = EEXInit;
  theNP->Display  = EEXDisplay;
  theNP->Execute  = NPEIterExecute;

  np = (NP_EITER *) theNP;
  np->PreProcess  = EEXPreProcess;
  np->Iter        = EEXSmoother;
  np->PostProcess = EEXPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitEIter	- Init this file

   SYNOPSIS:
   INT InitEIter ();

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

INT InitEIter ()
{
  INT i;

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  if (CreateClass(EITER_CLASS_NAME ".sciter",sizeof(NP_SCITER),SCITER_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(EITER_CLASS_NAME ".elmgc",sizeof(NP_ELMGC),ELmgcConstruct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(EITER_CLASS_NAME ".eex",sizeof(NP_EEX),EEXConstruct)) REP_ERR_RETURN (__LINE__);

  return (0);
}
