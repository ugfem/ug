// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  block.c			                                                                                */
/*																			*/
/* Purpose:   num procs to create blocks                                    */
/*																			*/
/*																			*/
/* Author: Klaus Johannsen                                                  */
/*         Sit                                                              */
/*         Universitaet Heidelberg                                          */
/*         INF 368                                                          */
/*         69120 Heidelberg                                                 */
/*         email: ug@ica3.uni-stuttgart.de                                  */
/*                                                                          */
/* History:   Sep 27, 2004 begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
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
#include "iter.h"
#include "project.h"
#include "disctools.h"
#include "block.h"

#include "ff_gen.h"
#include "ff.h"
#include "ugblas.h"
#include "blocking.h"

USING_UG_NAMESPACES

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
  NP_BLOCKING blocking;
} NP_ELEM_BLOCK;

typedef struct
{
  NP_BLOCKING blocking;

  INT depth;
} NP_SAB;

typedef struct
{
  NP_BLOCKING blocking;

  INT n;                                                                        /* desired block size           */
  INT r[MAXLEVEL];                                                      /* realized block size          */
} NP_DD;

typedef struct
{
  NP_BLOCKING blocking;
} NP_UB;

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

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   elemblock - numproc for element blocking

   DESCRIPTION:
   This numproc returns a element-wise blocking

   .vb
   npinit <name>;
   .ve

   'npexecute not possible'

   D*/
/****************************************************************************/

static INT ELEM_BLOCK_Init (NP_BASE *theNP, INT argc , char **argv)
{
  return (NP_ACTIVE);
}

static INT ELEM_BLOCK_Display (NP_BASE *theNP)
{
  return (0);
}

static INT ELEM_BLOCK_Blocking (NP_BLOCKING *theNP, GetMemProcPtr GetMem, INT level, MATDATA_DESC *A, BLOCKING_STRUCTUR *bs, INT *result)
{
  GRID *theGrid=NP_GRID(theNP,level);
  INT i,n,nt;
  ELEMENT *theElement;
  VECTOR **vl;

  bs->n=NT(theGrid);
  bs->nb=(INT *)(*GetMem)(sizeof(INT)*bs->n);
  for (theElement=FIRSTELEMENT(theGrid),n=nt=0; theElement!=NULL; theElement=SUCCE(theElement),n++)
  {
    bs->nb[n]=CORNERS_OF_ELEM(theElement);
    nt+=bs->nb[n];
  }
  bs->vb=(VECTOR***)(*GetMem)(sizeof(VECTOR**)*bs->n);
  vl=(VECTOR**)(*GetMem)(sizeof(VECTOR*)*nt);
  for (theElement=FIRSTELEMENT(theGrid),n=nt=0; theElement!=NULL; theElement=SUCCE(theElement),n++)
  {
    bs->vb[n]=vl+nt;
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      vl[nt++]=NVECTOR(CORNER(theElement,i));
  }

  return (0);
}

static INT ELEM_BLOCK_Construct (NP_BASE *theNP)
{
  NP_BLOCKING *np;

  theNP->Init=ELEM_BLOCK_Init;
  theNP->Display=ELEM_BLOCK_Display;
  theNP->Execute=NULL;

  np=(NP_BLOCKING *)theNP;
  np->PreProcess=NULL;
  np->Blocking=ELEM_BLOCK_Blocking;
  np->PostProcess=NULL;

  return(0);
}

/****************************************************************************/
/*D
   sab - numproc for simple algebraic blocking

   DESCRIPTION:
   This numproc returns a element-wise blocking

   .vb
   npinit <name>;
   .ve

   'npexecute not possible'

   D*/
/****************************************************************************/

static INT SAB_Init (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SAB *np=(NP_SAB *)theNP;

  if (ReadArgvINT ("depth",&np->depth,argc,argv)) np->depth=1;
  if (np->depth<0) return(NP_NOT_ACTIVE);

  return (NP_ACTIVE);
}

static INT SAB_Display (NP_BASE *theNP)
{
  NP_SAB *np=(NP_SAB *)theNP;

  UserWriteF(DISPLAY_NP_FORMAT_SS,"depth",BOOL_2_YN(np->depth));

  return (0);
}

static INT SAB_Reset (VECTOR *theV, INT depth)
{
  MATRIX *theM;

  SETVCUSED(theV,0);
  if (depth>0)
    for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      SAB_Reset(MDEST(theM),depth-1);
  return(0);
}

static INT SAB_GetN (VECTOR *theV, INT depth)
{
  INT n;
  MATRIX *theM;

  n=0;
  if (!VCUSED(theV))
  {
    n=1;
    SETVCUSED(theV,1);
  }
  if (depth>0)
    for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      if (!VCUSED(MDEST(theM)))
        n+=SAB_GetN(MDEST(theM),depth-1);

  return(n);
}

static INT SAB_SetBV (VECTOR *theV, INT depth, VECTOR **vl)
{
  INT i,n;
  MATRIX *theM;

  n=0;
  if (!VCUSED(theV))
  {
    n=1;
    SETVCUSED(theV,1);
    vl[0]=theV;
    vl++;
  }
  if (depth>0)
    for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      if (!VCUSED(MDEST(theM)))
      {
        i=SAB_SetBV(MDEST(theM),depth-1,vl);
        vl+=i;
        n+=i;
      }

  return(n);
}

static INT SAB_Blocking (NP_BLOCKING *theNP, GetMemProcPtr GetMem, INT level, MATDATA_DESC *A, BLOCKING_STRUCTUR *bs, INT *result)
{
  NP_SAB *np=(NP_SAB *)theNP;
  GRID *theGrid=NP_GRID(theNP,level);
  VECTOR **vl,*theV;
  INT n,i,nt;

  bs->n=NVEC(theGrid);
  bs->nb=(INT *)(*GetMem)(sizeof(INT)*bs->n);
  bs->vb=(VECTOR***)(*GetMem)(sizeof(VECTOR**)*bs->n);
  for (theV=FIRSTVECTOR(theGrid),n=nt=0; theV!=NULL; theV=SUCCVC(theV),n++)
  {
    SAB_Reset(theV,np->depth);
    i=SAB_GetN(theV,np->depth);
    nt+=i;
    bs->nb[n]=i;
  }
  vl=(VECTOR**)(*GetMem)(sizeof(VECTOR*)*nt);
  for (theV=FIRSTVECTOR(theGrid),n=nt=0; theV!=NULL; theV=SUCCVC(theV),n++)
  {
    bs->vb[n]=vl+nt;
    SAB_Reset(theV,np->depth);
    i=SAB_SetBV(theV,np->depth,vl+nt);
    nt+=i;
  }

  return (0);
}

static INT SAB_Construct (NP_BASE *theNP)
{
  NP_BLOCKING *np;

  theNP->Init=SAB_Init;
  theNP->Display=SAB_Display;
  theNP->Execute=NULL;

  np=(NP_BLOCKING *)theNP;
  np->PreProcess=NULL;
  np->Blocking=SAB_Blocking;
  np->PostProcess=NULL;

  return(0);
}

/****************************************************************************/
/*D
   dd - numproc for dd blocking on algebraic level

   DESCRIPTION:
   This numproc returns a element-wise blocking

   .vb
   npinit <name>;
   .ve

   'npexecute not possible'

   D*/
/****************************************************************************/

static INT DD_Init (NP_BASE *theNP, INT argc , char **argv)
{
  NP_DD *np=(NP_DD *)theNP;
  INT i;

  if (ReadArgvINT ("n",&np->n,argc,argv)) np->n=1;
  if (np->n<0) return(NP_NOT_ACTIVE);
  for (i=0; i<MAXLEVEL; i++) np->r[i]=0;

  return (NP_ACTIVE);
}

static INT DD_Display (NP_BASE *theNP)
{
  NP_DD *np=(NP_DD *)theNP;
  INT i;
  char buffer[32];

  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",(int)np->n);
  for (i=0; i<MAXLEVEL; i++)
    if (np->r[i]>0)
    {
      sprintf(buffer,"r[%d]",i);
      UserWriteF(DISPLAY_NP_FORMAT_SI,buffer,(int)np->r[i]);
    }

  return (0);
}

static INT DD_Blocking (NP_BLOCKING *theNP, GetMemProcPtr GetMem, INT level, MATDATA_DESC *A, BLOCKING_STRUCTUR *bs, INT *result)
{
  NP_DD *np=(NP_DD *)theNP;
  GRID *theGrid=NP_GRID(theNP,level);
  void *buffer;
  VECTOR **vlist;
  VECTOR *theV;
  MATRIX *theM;
  FIFO fifo;
  INT i,n,b,v_idx;

  /* perform reverse CutHill McGee ordering of vectors */
  n=NVEC(theGrid);
  buffer=(void *)(*GetMem)(sizeof(VECTOR*)*n); assert(buffer!=NULL);
  vlist = (VECTOR**)(*GetMem)(sizeof(VECTOR*)*n); assert(vlist!=NULL);
  fifo_init(&fifo,buffer,sizeof(VECTOR*)*n);
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV)) SETVCUSED(theV,0);
  fifo_in(&fifo,FIRSTVECTOR(theGrid)); SETVCUSED(FIRSTVECTOR(theGrid),1);
  while(!fifo_empty(&fifo))
  {
    theV=(VECTOR *)fifo_out(&fifo);
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      if (!VCUSED(MDEST(theM)))
      {
        fifo_in(&fifo,(void *)MDEST(theM));
        SETVCUSED(MDEST(theM),1);
      }
  }
  fifo_in(&fifo,(void *)theV);
  SETVCUSED(theV,0); i=0;
  while(!fifo_empty(&fifo))
  {
    theV=(VECTOR *)fifo_out(&fifo);
    vlist[i++]=theV;
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      if (VCUSED(MDEST(theM)))
      {
        fifo_in(&fifo,(void *)MDEST(theM));
        SETVCUSED(MDEST(theM),0);
      }
  }
  assert(i==n);
  for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV)) vlist[i++]=theV;
  for (i=0; i<n; i++) GRID_UNLINK_VECTOR(theGrid,vlist[i]);
  for (i=0; i<n; i++) GRID_LINK_VECTOR(theGrid,vlist[i],PRIO(vlist[i]));

  /* determine the blocking */
  b=(INT)ceil(((DOUBLE)n/(DOUBLE)np->n));
  np->r[level]=(INT)floor(((DOUBLE)n/(DOUBLE)b)+0.5);
  bs->n=(INT)ceil((DOUBLE)n/(DOUBLE)np->r[level]);
  bs->nb=(INT *)(*GetMem)(sizeof(INT)*bs->n);
  bs->vb=(VECTOR***)(*GetMem)(sizeof(VECTOR**)*bs->n);
  for(i=v_idx=0; i<bs->n; i++)
  {
    if (i<bs->n-1) bs->nb[i]=np->r[level];
    else bs->nb[i]=n-v_idx;
    assert(v_idx<n);
    bs->vb[i]=vlist+v_idx;
    v_idx+=np->r[level];
  }

  return (0);
}

static INT DD_Construct (NP_BASE *theNP)
{
  NP_BLOCKING *np;

  theNP->Init=DD_Init;
  theNP->Display=DD_Display;
  theNP->Execute=NULL;

  np=(NP_BLOCKING *)theNP;
  np->PreProcess=NULL;
  np->Blocking=DD_Blocking;
  np->PostProcess=NULL;

  return(0);
}

/****************************************************************************/
/*D
   ub - numproc for ub blocking on algebraic level

   DESCRIPTION:
   This numproc returns the (u)ltimate (b)locking

   .vb
   npinit <name>;
   .ve

   'npexecute not possible'

   D*/
/****************************************************************************/

static INT UB_Init (NP_BASE *theNP, INT argc , char **argv)
{
  return (NP_ACTIVE);
}

static INT UB_Display (NP_BASE *theNP)
{
  return (0);
}

static INT UB_Block (MATRIX *m)
{
  INT n;
  DOUBLE l,ll;
  DOUBLE_VECTOR p1,p2,q;
  VECTOR *v1,*v2;
  MATRIX *theM;

  v1=MDEST(m); v2=MDEST(MADJ(m));
  VectorPosition(v1,p1); VectorPosition(v2,p2); V_DIM_EUKLIDNORM_OF_DIFF(p1,p2,l);
  n=0;
  for (theM=MNEXT(VSTART(v1)); theM!=NULL; theM=MNEXT(theM))
  {
    VectorPosition(MDEST(theM),q);
    V_DIM_EUKLIDNORM_OF_DIFF(p1,q,ll);
    if (ll>3.0*l) n++;
  }
  for (theM=MNEXT(VSTART(v2)); theM!=NULL; theM=MNEXT(theM))
  {
    VectorPosition(MDEST(theM),q);
    V_DIM_EUKLIDNORM_OF_DIFF(p2,q,ll);
    if (ll>3.0*l) n++;
  }
  return(n);
}

static INT UB_WeiredElem (ELEMENT *e)
{
  DOUBLE min,max;

  min=3.14159265; max=0.0;
  if (MinMaxAngle(e,&min,&max)) assert(0);
  if (max>=0.666*3.14159265) return(1);
  return(0);
}

static INT UB_Blocking (NP_BLOCKING *theNP, GetMemProcPtr GetMem, INT level, MATDATA_DESC *A, BLOCKING_STRUCTUR *bs, INT *result)
{
  GRID *theGrid=NP_GRID(theNP,level);
  void *buffer;
  INT i,n,ne;
  VECTOR **vlist, *theV;
  MATRIX *theM;
  FIFO fifo;
  ELEMENT *e;

  /* prepare */
  n=NVEC(theGrid); ne=NT(theGrid);
  vlist = (VECTOR**)(*GetMem)(sizeof(VECTOR*)*n); assert(vlist!=NULL);
  bs->nb=(INT *)(*GetMem)(sizeof(INT)*(n+ne));
  bs->vb=(VECTOR ***)(*GetMem)(sizeof(VECTOR**)*(n+ne));
  buffer=(void *)(*GetMem)(sizeof(VECTOR*)*n); assert(buffer!=NULL);
  fifo_init(&fifo,buffer,sizeof(VECTOR*)*n);
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV)) SETVCUSED(theV,0);

  /* get strongly coupled line, defined by UG_Block(theM)==0  */
  for (bs->n=0; FIRSTVECTOR(theGrid)!=NULL; bs->n++)
  {
    fifo_in(&fifo,FIRSTVECTOR(theGrid)); SETVCUSED(FIRSTVECTOR(theGrid),1);
    n=0;
    while(!fifo_empty(&fifo))
    {
      theV=(VECTOR *)fifo_out(&fifo);
      vlist[n++]=theV;
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
        if (UB_Block(theM) && !VCUSED(MDEST(theM)))
        {
          fifo_in(&fifo,(void *)MDEST(theM));
          SETVCUSED(MDEST(theM),1);
        }
    }
    bs->nb[bs->n]=n;
    bs->vb[bs->n]=(VECTOR **)(*GetMem)(sizeof(VECTOR*)*n);
    for (i=0; i<n; i++)
    {
      bs->vb[bs->n][i]=vlist[i];
      GRID_UNLINK_VECTOR(theGrid,vlist[i]);
    }
  }
  for (n=0; n<bs->n; n++)
    for (i=0; i<bs->nb[n]; i++)
      GRID_LINK_VECTOR(theGrid,bs->vb[n][i],PRIO(bs->vb[n][i]));

  /* add vectors of weired elements */
  for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
    if (UB_WeiredElem(e))
    {
      GetVectorsOfNodes(e,&n,vlist);
      bs->nb[bs->n]=n;
      bs->vb[bs->n]=(VECTOR **)(*GetMem)(sizeof(VECTOR*)*n);
      for (i=0; i<n; i++)
        bs->vb[bs->n][i]=vlist[i];
      bs->n++;
    }

  return (0);
}

static INT UB_Construct (NP_BASE *theNP)
{
  NP_BLOCKING *np;

  theNP->Init=UB_Init;
  theNP->Display=UB_Display;
  theNP->Execute=NULL;

  np=(NP_BLOCKING *)theNP;
  np->PreProcess=NULL;
  np->Blocking=UB_Blocking;
  np->PostProcess=NULL;

  return(0);
}

/****************************************************************************/
/*
   InitBlock	- Init this file

   SYNOPSIS:
   INT InitBlock ();

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

INT NS_DIM_PREFIX InitBlocking ()
{
  if (CreateClass(BLOCKING_CLASS_NAME ".elemblock",sizeof(NP_ELEM_BLOCK),ELEM_BLOCK_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(BLOCKING_CLASS_NAME ".sab",sizeof(NP_SAB),SAB_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(BLOCKING_CLASS_NAME ".dd",sizeof(NP_DD),DD_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(BLOCKING_CLASS_NAME ".ub",sizeof(NP_UB),UB_Construct)) REP_ERR_RETURN (__LINE__);

  return (0);
}
