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
  NP_ELEM_BLOCK *np=(NP_ELEM_BLOCK *)theNP;

  return (NP_ACTIVE);
}

static INT ELEM_BLOCK_Display (NP_BASE *theNP)
{
  return (0);
}

static INT ELEM_BLOCK_Blocking (NP_BLOCKING *theNP, GetMemProcPtr GetMem, INT level, MATDATA_DESC *A, BLOCKING_STRUCTUR *bs, INT *result)
{
  NP_ELEM_BLOCK *np=(NP_ELEM_BLOCK *)theNP;
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
  VECTOR **vl,*theV,*theW;
  MATRIX *theM;
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
  if (CreateClass(BLOCKING_CLASS_NAME ".sab",sizeof(NP_ELEM_BLOCK),SAB_Construct)) REP_ERR_RETURN (__LINE__);

  return (0);
}
