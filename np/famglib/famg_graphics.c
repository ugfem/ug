// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  famg_graphics.c												*/
/*																			*/
/* Purpose:   graphic functions                                                                         */
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   August 97 begin, Stuttgart									*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#ifdef UG_DRAW

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>


#include "wop.h"
#include "wpm.h"
#include "misc.h"
#include "evm.h"
#include "cw.h"
#include "graph.h"
#include "gm.h"

#include "famginterface.h"

/* RCS_ID
   $Header$
 */

static long BlackColor; /* Black  */
static long RedColor; /* Red  */
static int GlobalCounter;
static int GlobalN;
static int GlobalCGCounter;
static int GlobalCGN;
static int GlobalNF;
static int GlobalLevel;
static int GlobalMaxLevel;
static void **GlobalNodePtr;
static void **GlobalCGNodePtr;


struct FAMGPlotObject
{
  struct PlotObjHead theHead;      /* the head */
  int level;
};


static INT SetFAMGGraph (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  struct FAMGPlotObject *theObj;
  int l;

  theObj = (struct FAMGPlotObject *) &(thePlotObj->theExternObject);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

  if (ReadArgvINT("l",&l,argc,argv)) l = 0;
  theObj->level = l;

  return (ACTIVE);
}

static INT  DisplayFAMGGraph (PLOTOBJ *thePlotObj)
{
  struct FAMGPlotObject *theObj;

  theObj = (struct FAMGPlotObject *) &(thePlotObj->theExternObject);

  UserWriteF("level: %d \n",theObj->level);
  return (0);
}

static INT PreProcessFAMGGraph (PICTURE *thePicture, WORK *theWork)
{
  OUTPUTDEVICE *theOD;
  struct FAMGPlotObject *theObj;
  int level,maxlevel;

  theObj = (struct FAMGPlotObject *) &(PIC_PO(thePicture)->theExternObject);
  theOD  = PIC_OUTPUTDEV(thePicture);

  BlackColor = theOD->black;
  RedColor = theOD->red;
  GlobalCounter = 0;
  GlobalCGCounter = 0;
  GlobalMaxLevel = maxlevel = FAMG_GetMaxLevel();
  level = theObj->level;
  if(level > maxlevel)
  {
    level = maxlevel;
    UserWriteF("plotting level %d \n",level);
  }

  GlobalNodePtr = FAMG_GetExtraPtr(level);
  GlobalN = FAMG_GetN(level);
  GlobalNF = FAMG_GetNF(level);
  GlobalLevel = level;


  if(level < maxlevel)
  {
    GlobalCGNodePtr = FAMG_GetExtraPtr(level+1);
    GlobalCGN = FAMG_GetN(level+1);
  }
  return(0);
}




static INT EvalFAMGGraph (DRAWINGOBJ *theDO, INT *end)
{
  VERTEX *vertex, *nbvertex;
  FAMG_MatrixPtr mat, *matptr;
  FAMG_TransferEntry *trans;
  DOUBLE_VECTOR mypos,nbpos;
  int j;

  UgSetLineWidth(1);

  if(GlobalCounter == GlobalN) {*end = 1; return (0);}

  vertex = (VERTEX *) GlobalNodePtr[GlobalCounter];
  V_DIM_COPY(CVECT(vertex),mypos);
  matptr = FAMG_GetMatrixPtr(GlobalLevel,GlobalCounter);
  mat = *matptr;

  if((mat.index)->type)
  {
    DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO);
    DO_2c(theDO) = 1; DO_inc(theDO);
    DO_2l(theDO) = BlackColor; DO_inc(theDO);
    DO_2s(theDO) = EMPTY_CIRCLE_MARKER; DO_inc(theDO);
    DO_2s(theDO) = 8; DO_inc(theDO);
    V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
  }
  else
  {
    DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO);
    DO_2c(theDO) = 1; DO_inc(theDO);
    DO_2l(theDO) = BlackColor; DO_inc(theDO);
    DO_2s(theDO) = FILLED_CIRCLE_MARKER; DO_inc(theDO);
    DO_2s(theDO) = 8; DO_inc(theDO);
    V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
  }


  /* print id */
  DO_2c(theDO) = DO_TEXT; DO_inc(theDO);
  DO_2l(theDO) = BlackColor; DO_inc(theDO);
  DO_2c(theDO) = TEXT_REGULAR; DO_inc(theDO);
  DO_2c(theDO) = TEXT_NOT_CENTERED; DO_inc(theDO);
  DO_2s(theDO) = 6; DO_inc(theDO);
  V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
  sprintf(DO_2cp(theDO),"%d",GlobalCounter);
  DO_inc_str(theDO);

  /* print edges */
  while(mat.nc > 0)
  {
    (mat.index)++; (mat.nc)--;
    j = (mat.index)->id;
    nbvertex = (VERTEX *) GlobalNodePtr[j];
    V_DIM_COPY(CVECT(nbvertex),nbpos);
    DO_2c(theDO) = DO_LINE; DO_inc(theDO);
    DO_2l(theDO) =BlackColor; DO_inc(theDO);
    V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    V2_COPY(nbpos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
  }

  trans = FAMG_GetTransferEntry(GlobalLevel,GlobalCounter);
  if (trans != NULL)
  {
    for(trans = trans->next; trans != NULL; trans = trans->next)
    {
      j = trans->id.f0;
      nbvertex = (VERTEX *) GlobalNodePtr[j];
      V_DIM_COPY(CVECT(nbvertex),nbpos);
      DO_2c(theDO) = DO_LINE; DO_inc(theDO);
      DO_2l(theDO) =RedColor; DO_inc(theDO);
      V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
      V2_COPY(nbpos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    }
  }

  GlobalCounter++;
  DO_2c(theDO) = DO_NO_INST;

  if(GlobalCounter == GlobalN) *end = 1;
  else *end = 0;

  return(0);
}


static INT PostProcessFAMGGraph (PICTURE *thePicture, WORK *theWork)
{
  return(0);
}


INT InitFAMGGraph (void)
{
  PLOTOBJHANDLING *thePOH;
  WORKPROCS *theWP;
  EXTERNWORK *theEXW;
  PLOTOBJTYPE *thePOT;

  /* create WorkHandling for 'FAMGGraph' */
  if ((thePOH=CreatePlotObjHandling ("FAMGGraph"))    == NULL)
    return (1);

  POH_DYNAMIC_INFO(thePOH) = NULL;
  POH_CLICKACTION(thePOH)  = NULL;

  /* draw work */
  POH_NBCYCLES(thePOH,DRAW_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
  WP_WORKMODE(theWP) = EXTERN;
  theEXW = WP_EXTERNWISE(theWP);
  theEXW->EXT_PreProcessProc = PreProcessFAMGGraph;
  theEXW->EXT_EvaluateProc   = EvalFAMGGraph;
  theEXW->EXT_ExecuteProc        = Draw2D;
  theEXW->EXT_PostProcessProc = NULL;   /* PostProcessFAMGGraph; */

  if ((thePOT=GetPlotObjType("FAMGGraph"))    == NULL)
    return (1);
  thePOT->Dimension = TYPE_2D;
  thePOT->SetPlotObjProc = SetFAMGGraph;
  thePOT->DispPlotObjProc = DisplayFAMGGraph;

  return(0);
}

#endif /* UG_DRAW */
