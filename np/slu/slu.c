// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  slu.c	                                                                                                        */
/*																			*/
/* Purpose:   interface to SuperLU from:									*/
/*                James W. Demmel, UCB										*/
/*                John R. Gilbert, Xerox Paolo Alto Research Center			*/
/*                Xiaoye S. Li, NERSC										*/
/*																			*/
/* Author:	  Klaus Johannsen	                                                                                */
/*			  IWR/Technische Simulation										*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69117 Heidelberg			                                                                */
/*																			*/
/* History:   October 15, 2001                                                                          */
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
#include <assert.h>

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
#include "iter.h"
#include "slu.h"

#include "namespace.h"

USING_UG_NAMESPACES

extern void    SLU_GetStat (float *ftime, float *fflops, float *stime, float *sflops);
extern int     lsame_ (char *, char *);
extern int     sp_ienv (int);

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

typedef enum {
  SLU_NC,
  SLU_NR,
  SLU_SC,
  SLU_SR,
  SLU_NCP,
  SLU_DN
} Stype_t;

typedef enum {
  SLU__S,
  SLU__D,
  SLU__C,
  SLU__Z
} Dtype_t;

typedef enum {
  SLU_GE,
  SLU_TRLU,
  SLU_TRUU,
  SLU_TRL,
  SLU_TRU,
  SLU_SYL,
  SLU_SYU,
  SLU_HEL,
  SLU_HEU
} Mtype_t;

typedef struct {
  Stype_t Stype;
  Dtype_t Dtype;
  Mtype_t Mtype;
  int nrow;
  int ncol;
  void *Store;
} SuperMatrix;

typedef struct {
  int panel_size;
  int relax;
  double diag_pivot_thresh;
  double drop_tol;
} factor_param_t;

typedef struct {
  float for_lu;
  float total_needed;
  int expansions;
} mem_usage_t;

typedef struct {
  int nnz;
  void *nzval;
  int  *colind;
  int  *rowptr;
} NRformat;

typedef struct {
  int lda;
  void *nzval;
} DNformat;

typedef struct {
  int nnz;
  void *nzval;
  int  *rowind;
  int  *colptr;
} NCformat;

typedef struct {
  int nnz;
  int nsuper;
  void *nzval;
  int  *nzval_colptr;
  int  *rowind;
  int *rowind_colptr;
  int *col_to_sup;
  int *sup_to_col;
} SCformat;

typedef struct {
  NP_ITER iter;

  INT optimize;                         /* 0: natural ordering							*/
  /* 1: minimum degree on structure of A'*A		*/
  /* 2: minimum degree on structure of A'+A		*/
  /* 3: approx. min. degree for unsym. matrices	*/
  INT n;                                        /* size of problem								*/
  INT Annz;                                     /* nonzeros of A								*/
  double rpg;                                   /* growth rate of pivot element					*/
  double rcond;                         /* reciprocal condition number					*/
  INT display;

  SuperMatrix A,L,U,B,X;    /* SuperLU specific stuff                                           */
  mem_usage_t mem_usage;
  factor_param_t iparam;
  int *perm_r,*perm_c,*etree;
  double *rhsb,*rhsx,*R,*C,ferr,berr;
} NP_SLU;

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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

static INT SLU_MarkKey_nb;
static INT SLU_MarkKey[2];
static HEAP *SLU_Heap;

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

int *intMalloc (int);
float *floatMalloc (int);
double *doubleMalloc (int);
void dCreate_Dense_Matrix (SuperMatrix *, int, int, double *, int, Stype_t, Dtype_t, Mtype_t);
void get_perm_c (int, SuperMatrix *, int *);
void dgssvx (int decompose, char *, char *, char *, SuperMatrix *, factor_param_t *, int *, int *, int *, char *, double *, double *,
             SuperMatrix *, SuperMatrix *, void *, int, SuperMatrix *, SuperMatrix *, double *, double *, double *, double *, mem_usage_t *, int *);
void superlu_free (void*);
void *superlu_malloc (int);

/****************************************************************************/
/*																			*/
/* interface																*/
/*																			*/
/****************************************************************************/

static int ReadMatrix (SuperMatrix *A, GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *J)
{
  INT i,j,n,rindex,rtype,rcomp,nnz,ctype,ccomp,cindex;
  SHORT *comp;
  int *fe,*dest;
  double *a;
  VECTOR *theV,*theW;
  MATRIX *theM;
  NRformat *Astore;

  for (theV=FIRSTVECTOR(theGrid),rindex=0; theV!=NULL; theV=SUCCVC(theV))
  {
    rtype=VTYPE(theV);
    rcomp=VD_NCMPS_IN_TYPE(x,rtype);
    VINDEX(theV)=rindex;
    rindex+=rcomp;
  }
  n=rindex;
  fe=intMalloc(n+1); if (fe==NULL) REP_ERR_RETURN(1);
  for (theV=FIRSTVECTOR(theGrid),nnz=0; theV!=NULL; theV=SUCCVC(theV))
  {
    rindex=VINDEX(theV);
    rtype=VTYPE(theV);
    rcomp=VD_NCMPS_IN_TYPE(x,rtype);
    for (i=0; i<rcomp; i++)
    {
      fe[rindex+i]=nnz;
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      {
        theW=MDEST(theM);
        ctype=VTYPE(theW);
        ccomp=VD_NCMPS_IN_TYPE(x,ctype);
        nnz+=ccomp;
      }
    }
  }
  fe[n]=nnz;
  dest=intMalloc(nnz); if (dest==NULL) REP_ERR_RETURN(1);
  a=doubleMalloc(nnz); if (a==NULL) REP_ERR_RETURN(1);
  for (theV=FIRSTVECTOR(theGrid),nnz=0; theV!=NULL; theV=SUCCVC(theV))
  {
    rindex=VINDEX(theV);
    rtype=VTYPE(theV);
    rcomp=VD_NCMPS_IN_TYPE(x,rtype);
    for (i=0; i<rcomp; i++)
    {
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      {
        theW=MDEST(theM);
        cindex=VINDEX(theW);
        ctype=VTYPE(theW);
        ccomp=VD_NCMPS_IN_TYPE(x,ctype);
        comp=MD_MCMPPTR_OF_RT_CT(J,rtype,ctype);
        for (j=0; j<ccomp; j++)
        {
          dest[nnz]=cindex+j;
          a[nnz]=MVALUE(theM,comp[i*ccomp+j]);
          nnz++;
        }
      }
    }
  }

  /* creating supermatrix */
  A->Stype=SLU_NR;    A->Dtype=SLU__D;    A->Mtype=SLU_GE;
  A->nrow=n;      A->ncol=n;
  A->Store=(void*)superlu_malloc(sizeof(NRformat));
  if (A->Store==NULL) REP_ERR_RETURN(1);
  Astore = (NRformat*)A->Store;
  Astore->nnz=fe[n];
  Astore->nzval=a;
  Astore->colind=dest;
  Astore->rowptr=fe;

  /* return */
  return(0);
}

int ReadVector (INT n, SuperMatrix *B, GRID *theGrid, VECDATA_DESC *x)
{
  int i,index,type;
  double *v;
  VECTOR *theV;
  DNformat *Bstore;
  SHORT *comp;

  Bstore = (DNformat*)B->Store;
  v=(double*)Bstore->nzval;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    type=VTYPE(theV);
    comp=VD_CMPPTR_OF_TYPE(x,type);
    index=VINDEX(theV);
    for (i=0; i<VD_NCMPS_IN_TYPE(x,type); i++)
      v[index+i]=VVALUE(theV,comp[i]);
  }

  return(0);
}

int WriteVector (INT n, SuperMatrix *B, GRID *theGrid, VECDATA_DESC *x)
{
  int i,index,type;
  double *v;
  VECTOR *theV;
  DNformat *Bstore;
  SHORT *comp;

  Bstore = (DNformat*)B->Store;
  v = (double*)Bstore->nzval;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    type=VTYPE(theV);
    comp=VD_CMPPTR_OF_TYPE(x,type);
    index=VINDEX(theV);
    for (i=0; i<VD_NCMPS_IN_TYPE(x,type); i++)
      VVALUE(theV,comp[i])=v[index+i];
  }

  return(0);
}

void *SLU_Malloc (int size)
{
  return((void*)GetTmpMem(SLU_Heap,size,SLU_MarkKey[SLU_MarkKey_nb]));
}

static INT SLU_Output (NP_BASE* theNP)
{
  NP_SLU *np;
  SCformat *Lstore;
  NCformat *Ustore;
  float ftime,fflops,stime,sflops;
  char buffer[128];

  np = (NP_SLU *) theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",np->n);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Rpg",np->rpg);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Rcn",np->rcond);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Ferror",np->ferr);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Berror",np->berr);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"A entries",np->Annz);
  Lstore = (SCformat *) np->L.Store;
  if (Lstore!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SI,"L entries",Lstore->nnz);
  Ustore = (NCformat *) np->U.Store;
  if (Ustore!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SI,"U entries",Ustore->nnz);
  if (Lstore!=NULL && Ustore!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SI,"L+U entries",Lstore->nnz + Ustore->nnz - np->n);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Mem L/U (MB)",np->mem_usage.for_lu/1e6);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Mem tot (MB)",np->mem_usage.total_needed/1e6);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"Expansions",np->mem_usage.expansions);
  SLU_GetStat(&ftime,&fflops,&stime,&sflops);
  UserWriteF("\n");
  sprintf(buffer,"%.2f s",ftime);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"Factor time",buffer);
  sprintf(buffer,"%.2f MFlops",fflops*1e-6);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"Factor flops",buffer);
  sprintf(buffer,"%.2f s",stime);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"Solve time",buffer);
  sprintf(buffer,"%.2f MFlops",sflops*1e-6);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"Solve flops",buffer);

  return (0);
}

static INT SLUPreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_SLU *np = (NP_SLU*)theNP;
  HEAP *theHeap = MGHEAP(NP_MG(theNP));
  GRID *theGrid = NP_GRID(theNP,level);
  VECTOR *theV;
  MATRIX *theM;
  NCformat *Ustore;
  SCformat *Lstore;
  char text[DISPLAY_WIDTH+4];
  char fact[1],equed[1],trans[1],refact[1];
  int firstfact,info,n,m,i;
  double *a;
  void *work;

        #ifdef ModelP
  if (me != master) return(0);
        #endif

  /* init heap management */
  SLU_MarkKey_nb=0;
  SLU_Heap=theHeap;
  MarkTmpMem(SLU_Heap,&(SLU_MarkKey[SLU_MarkKey_nb]));

  /* control parameters */
  *fact='E'; *equed='B'; *trans='N'; *refact='N';
  firstfact=lsame_(fact,"F") || lsame_(refact,"Y");
  np->iparam.panel_size=sp_ienv(1); np->iparam.relax=sp_ienv(2);
  np->iparam.diag_pivot_thresh=1.0; np->iparam.drop_tol=-1;

  /* cp matrix: ug --> SuperLU */
  if (ReadMatrix(&np->A,theGrid,b,A)) REP_ERR_RETURN(1);
  np->n=n=np->A.nrow; m=np->A.ncol;

  /* cp vector: ug --> SuperLU */
  if (!(np->rhsb=doubleMalloc(m))) REP_ERR_RETURN(1);
  dCreate_Dense_Matrix(&np->B,m,1,np->rhsb,m,SLU_DN,SLU__D,SLU_GE);
  if (ReadVector(n,&np->B,theGrid,b)) REP_ERR_RETURN(1);
  if (!(np->rhsx=doubleMalloc(m))) REP_ERR_RETURN(1);
  dCreate_Dense_Matrix(&np->X,m,1,np->rhsx,m,SLU_DN,SLU__D,SLU_GE);
  np->Annz=((NCformat*)np->A.Store)->nnz;

  /* solve */
  if (!(np->etree =intMalloc(n))) REP_ERR_RETURN(1);
  if (!(np->perm_r=intMalloc(m))) REP_ERR_RETURN(1);
  if (!(np->perm_c=intMalloc(n))) REP_ERR_RETURN(1);
  get_perm_c(np->optimize,&np->A,np->perm_c);
  if (!(np->R=(double *)superlu_malloc(np->A.nrow*sizeof(double)))) REP_ERR_RETURN(1);
  if (!(np->C=(double *)superlu_malloc(np->A.ncol*sizeof(double)))) REP_ERR_RETURN(1);
  dgssvx(1,fact,trans,refact,&np->A,&np->iparam,np->perm_c,np->perm_r,np->etree,equed,np->R,np->C,&np->L,&np->U,work,0,
         &np->B,&np->X,&np->rpg,&np->rcond,&np->ferr,&np->berr,&np->mem_usage,&info);
  if (info!=0) REP_ERR_RETURN(1);

  if (np->display>PCR_NO_DISPLAY)
  {
    CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
    UserWriteF(text);
    if (SLU_Output((NP_BASE*)theNP)) REP_ERR_RETURN(1);
  }

  return (0);
}

static INT SLUIter (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SLU *np = (NP_SLU*)theNP;
  GRID *theGrid = NP_GRID(theNP,level);
  int info,bl,firstfact;
  char fact[1],equed[1],trans[1],refact[1];
  void *work;
  VEC_SCALAR defect;

        #ifdef ModelP
  if (me != master) return(0);
        #endif

  /* control parameters */
  *fact='F'; *equed='B'; *trans='N'; *refact='N';

  /* solve */
  if (ReadVector(np->n,&np->B,theGrid,b)) REP_ERR_RETURN(1);
  SLU_MarkKey_nb=1;
  MarkTmpMem(SLU_Heap,&(SLU_MarkKey[SLU_MarkKey_nb]));
  dgssvx(0,fact,trans,refact,&np->A,&np->iparam,np->perm_c,np->perm_r,np->etree,equed,np->R,np->C,&np->L,&np->U,work,0,
         &np->B,&np->X,&np->rpg,&np->rcond,&np->ferr,&np->berr,&np->mem_usage,&info);
  ReleaseTmpMem(SLU_Heap,SLU_MarkKey[SLU_MarkKey_nb]);
  SLU_MarkKey_nb=0;
  if (info!=0) REP_ERR_RETURN(1);
  if (WriteVector(np->n,&np->X,theGrid,x)) REP_ERR_RETURN(1);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK) REP_ERR_RETURN(1);

  return(0);
}

static INT SLUPostProcess (NP_ITER *theNP, INT level,VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SLU *np;

        #ifdef ModelP
  if (me != master) return(0);
        #endif

  np = (NP_SLU *) theNP;
  ReleaseTmpMem(SLU_Heap,SLU_MarkKey[SLU_MarkKey_nb]);

  return(0);
}

static INT SLUInit (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_SLU *np;

  np = (NP_SLU *) theNP;
  if (ReadArgvINT ("o",&np->optimize,argc,argv)) np->optimize=3;
  if (np->optimize<0 || np->optimize>3) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display=ReadArgvDisplay(argc,argv);

  return (NP_ACTIVE);
}

static INT SLUDisplay (NP_BASE *theNP)
{
  NP_SLU *np;

  np = (NP_SLU *) theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SI,"o",np->optimize);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"display",np->display);
  UserWriteF("\n");

  return (SLU_Output(theNP));
}

INT SLUExecute (NP_BASE *theNP, INT argc , char **argv)
{
  return(1);
}

static INT SLUConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init     = SLUInit;
  theNP->Display  = SLUDisplay;
  theNP->Execute  = SLUExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess  = SLUPreProcess;
  np->Iter        = SLUIter;
  np->PostProcess = SLUPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitSLU	- Init this file

   SYNOPSIS:
   INT InitSLU ();

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

INT InitSLU ()
{
  if (CreateClass(ITER_CLASS_NAME ".slu",sizeof(NP_SLU),SLUConstruct)) REP_ERR_RETURN (__LINE__);

  return(0);
}
