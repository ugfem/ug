// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  wop.c                                                                                                                 */
/*																			*/
/* Purpose:   work functions on pictures									*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: johannse@iwr1.iwr.uni-heidelberg.de                                 */
/*																			*/
/* History:   21.06.93 begin, ug version ug21Xmas3d                                             */
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

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compiler.h"
#include "wop.h"
#include "wpm.h"
#include "plotproc.h"
#include "misc.h"
#include "evm.h"
#include "simplex.h"
#include "cw.h"
#include "graph.h"
#include "gm.h"

#ifdef __TWODIM__
#include "shapes2d.h"
#endif
#ifdef __THREEDIM__
#include "shapes3d.h"
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

/* definition of the arrow */
#define ARR_ALPHA                                               0.7
#define ARR_SIN                                                 0.5
#define ARR_COS                                            -0.866

/* ctrl entries of element ctrl */
INT ce_VSIDES;
#define VSIDES_LEN                                              4
#define VSIDES(p)                                               CW_READ(p,ce_VSIDES)
#define SETVSIDES(p,n)                                  CW_WRITE(p,ce_VSIDES,n)
#define VIEWABLE(p,i)                                   (VSIDES(p) & (1<<i))

INT ce_NORDER;
#define NORDER_LEN                                              5
#define NORDER(p)                                               CW_READ(p,ce_NORDER)
#define SETNORDER(p,n)                                  CW_WRITE(p,ce_NORDER,n)

INT ce_COUNT;
#define COUNT_LEN                                               4
#define COUNT(p)                                                CW_READ(p,ce_COUNT)
#define SETCOUNT(p,n)                                   CW_WRITE(p,ce_COUNT,n)

INT ce_CUTMODE;
#define CUTMODE_LEN                                     2
#define CUTMODE(p)                                              CW_READ(p,ce_CUTMODE)
#define SETCUTMODE(p,n)                                 CW_WRITE(p,ce_CUTMODE,n)

/* values for CUTMODE */
#define CM_BEHIND                                               0
#define CM_INTERSECT                                    1
#define CM_INFRONT                                              2


/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef void (*ProjectionProcPtr)(COORD *, COORD_POINT *);

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/************************************************************************/
/************ ordinary static variables   *******************************/
/************************************************************************/

#ifdef __DO_HEAP_USED__
static INT Heap_Used_Min;
static INT Heap_Used_Max;
#endif


/* internal function used to create value of NodeOrder from code determined in function 'CalcNodeOrder' */
/* if ">=" means nearer, than OrderIndex[111111] = 0 is the index for x1>=x2>=x3>=x4					*/
static INT OrderIndex[64] = {23,17,-1,15,  21,-1,11,9,   -1,-1,-1,14,  -1,-1,-1,8,
                             -1,-1,-1,-1,  20,-1,10,-1,  -1,-1,-1,-1,  -1,-1,7,6,
                             22,16,-1,-1,  -1,-1,-1,-1,  -1,13,-1,12,  -1,-1,-1,-1,
                             19,-1,-1,-1,  18,-1,-1,-1,   5,3,-1,2,         4,-1,1,0        };

static INT NoOfViewableSides[16] =      {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4}; /* get number of viewable sides	*/

/* unit vectors */
static COORD ex[3] = {1.0, 0.0, 0.0};
static COORD ey[3] = {0.0, 1.0, 0.0};
static COORD ez[3] = {0.0, 0.0, 1.0};

/*----------- variables describing taransformations --------------------*/
static COORD ObsTrafo[16], InvObsTrafo[16];
static COORD CutTrafo[16], InvCutTrafo[16];
static INT CUT_CutExisting;
static COORD CUT_CutNormal[3];
static INT CUT_CutAtFront;
static ProjectionProcPtr OBS_ProjectProc;
static INT OBS_Perspective;
static COORD OBS_PerspCorr[2];
static COORD OBS_ViewDirection[3];
static COORD OBS_ViewPlaneDist;



/****************************************************************************/
/************ variables used for communication of functions *****************/
/****************************************************************************/

/*---------- variables use by OrderElements etc ----------------------------*/
static VIEWEDOBJ                        *OE_ViewedObj;

/*---------- variables use by GetFirst/NextElement... ----------------------*/
static MULTIGRID                                        *GE_MG;
static INT GE_fromLevel,GE_toLevel;


/*---------- input variables of 'EW_ElementEval2D/3D' ----------------------*/
/* defines for 2D/3D */
#define PLOT_COPY                       YELLOW                  /* values for 'Elem2Plot'		*/
#define PLOT_IRR                        GREEN
#define PLOT_REG                        REGULAR_CLASS
#define PLOT_ALL                        (REGULAR_CLASS+1) /* this means all levels		*/

#define COLOR_COPY                      YELLOW                  /* values for 'NoColor' and     */
#define COLOR_IRR                       GREEN                   /* 'Color'						*/
#define COLOR_REG                       REGULAR_CLASS
#define COLOR_LOWER_LEVEL       (REGULAR_CLASS+1)
#define COLOR_EDGE                      (REGULAR_CLASS+2)

/* defines 2D */
#define EE2D_TEXTSIZE                   8

#define COLOR_BND                       (REGULAR_CLASS+3)
#define COLOR_ELEMID            (REGULAR_CLASS+4)

/* defines 3D */
#define EE3D_TEXTSIZE                   8

#define COLOR_CUT_EDGE          (REGULAR_CLASS+3)

/* 2D */
static INT EE2D_Elem2Plot[10];          /* 1 if element has to be plotted			*/
static long EE2D_NoColor[10];   /* 1 if no color (background color) used	*/
static long EE2D_Color[10];             /* colors used								*/
static INT EE2D_MaxLevel;               /* level considered to be the top level         */
static INT EE2D_ElemID;                 /* 1 if element ID has to be plotted		*/


/* 3D */
static INT EE3D_Elem2Plot[10];          /* 1 if element has to be plotted			*/
static long EE3D_NoColor[10];   /* 1 if no color (background color) used	*/
static long EE3D_Color[10];             /* colors used								*/
static INT EE3D_MaxLevel;               /* level considered to be the top level         */


/*---------- working variables of 'NW_NodesEval2D' -------------------------*/
static long NE_IDColor;                 /* color of node ID's                       */
static long NE_BndMarkerColor;  /* color of bnd marks						*/
static long NE_InnerMarkerColor; /* color of inner marks                                        */
static INT NE_EvalNodeID;               /* 1 if to evaluate                                             */
static INT NE_EvalInnerNode;    /* 1 if to evaluate                                             */
static INT NE_EvalBndNode;              /* 1 if to evaluate                                             */
static short NE_InnerMarker;    /* marker for inner nodes					*/
static short NE_BndMarker;              /* marker for bnd nodes                                         */
static short NE_InnerMarkerSize; /* markersize for inner nodes				*/
static short NE_BndMarkerSize;  /* markersize for bnd nodes                             */


/*---------- working variables of 'EW_FindElement3D' -----------------------*/
static COORD_POINT FE2D_MousePos;


/*---------- working variables of 'EW_FindNode2D' --------------------------*/
#define FN2D_INVSIZE                    3
#define FN2D_ACC                                3

static COORD_POINT FN2D_MousePos;
static INT FN2D_found;

/*---------- working variables of 'EW_FindElement2D' -----------------------*/
static INT FE2D_found;

/*---------- working variables of 'EW_FindElement3D' -----------------------*/
static COORD_POINT FE3D_MousePos;


/*---------- working variables of 'EW_EScalar2D' ---------------------------*/
static PreprocessingProcPtr EScalar2D_PreProcess;
static ElementPlotProcPtr EScalar2D_EvalFct;
static COORD EScalar2D_V2C_factor;
static COORD EScalar2D_V2C_offset;
static INT EScalar2D_mode;
static INT EScalar2D_depth;
static INT EScalar2D_numOfContours;
static DOUBLE EScalar2D_ContValues[PO_MAXCONTOURS];
static DOUBLE EScalar2D_minValue;
static DOUBLE EScalar2D_maxValue;


/*---------- working variables of 'EW_EScalar3D' ---------------------------*/
static PreprocessingProcPtr EScalar3D_PreProcess;
static ElementPlotProcPtr EScalar3D_EvalFct;
static COORD EScalar3D_V2C_factor;
static COORD EScalar3D_V2C_offset;
static INT EScalar3D_mode;
static INT EScalar3D_depth;
static INT EScalar3D_numOfContours;
static DOUBLE EScalar3D_ContValues[PO_MAXCONTOURS];
static DOUBLE EScalar3D_minValue;
static DOUBLE EScalar3D_maxValue;


/*---------- working variables of 'EW_EVector3D' ---------------------------*/
#define RASTERPOINTS_MAX                200

static ElementVectorProcPtr EVector_EvalFct;
static COORD EVector_rastersize;
static INT EVector_cutvector;
static COORD EVector_V2L_factor;
static COORD EVector_CutLenFactor;
static long EVector_ColorCut;

/* 2D */
static long EVector2D_ColorNormal;

/* 3D */
static INT EVector3D_projectvector;
static COORD EVector3D_V2C_factor;
static COORD EVector3D_V2C_offset;


/*---------- working variables of 'GetNode...' routines --------------------*/
static MULTIGRID        *GNode_MG;
static INT GNode_fromLevel;
static INT GNode_toLevel;

/*---------- working variables of 'GetElement...' routines -----------------*/
static MULTIGRID        *GElem_MG;
static INT GElem_fromLevel;
static INT GElem_toLevel;

/*---------- working variables of 'FindRange' routines ---------------------*/
static INT GEN_FR_put;
static DOUBLE GEN_FR_min;
static DOUBLE GEN_FR_max;


/*---------- working variables of 'WorkOnPicture' routines -----------------*/
static OUTPUTDEVICE                             *WOP_OutputDevice;
static PICTURE                                          *WOP_Picture;
static PLOTOBJ                                          *WOP_PlotObj;
static VIEWEDOBJ                                        *WOP_ViewedObj;
static WORK                                             *WOP_Work;
static PLOTOBJHANDLING                          *WOP_PlotObjHandling;
static MULTIGRID                                        *WOP_MG;
static WORKPROCS                                        *WOP_WorkProcs;
static INT WOP_ViewDim;
static INT WOP_WorkMode;
static ELEMENT                                          *WOP_Element;
static NODE                                             *WOP_Node;
static VECTOR                                           *WOP_Vector;
static DRAWINGOBJ WOP_DrawingObject[DO_SIZE];

static GEN_ExecuteProcPtr WOP_GEN_ExecuteProc;
static GEN_PostProcessProcPtr WOP_GEN_PostProcessProc;
static GEN_PreProcessProcPtr WOP_GEN_PreProcessProc;
static EW_GetFirstElementProcPtr WOP_EW_GetFirstElementProc;
static EW_GetNextElementProcPtr WOP_EW_GetNextElementProc;
static EW_EvaluateProcPtr WOP_EW_EvaluateProc;
static NW_GetFirstNodeProcPtr WOP_NW_GetFirstNodeProc;
static NW_GetNextNodeProcPtr WOP_NW_GetNextNodeProc;
static NW_EvaluateProcPtr WOP_NW_EvaluateProc;
static VW_GetFirstVectorProcPtr WOP_VW_GetFirstVectorProc;
static VW_GetNextVectorProcPtr WOP_VW_GetNextVectorProc;
static VW_EvaluateProcPtr WOP_VW_EvaluateProc;
static EXT_EvaluateProcPtr WOP_EXT_EvaluateProc;

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   CreatePlotObjHandling - Allocate a new PLOTOBJHANDLING

   SYNOPSIS:
   static PLOTOBJHANDLING	*CreatePlotObjHandling (char *PlotObjTypeName)

   PARAMETERS:
   .  PlotObjTypeName - the name

   DESCRIPTION:
   This function allocates a new 'PLOTOBJHANDLING'.

   RETURN VALUE:
   pointer to 'PLOTOBJHANDLING'

   NULL
   D*/
/****************************************************************************/

static PLOTOBJHANDLING  *CreatePlotObjHandling (char *PlotObjTypeName)
{
  /* allocate PLOTOBJHANDLING envItem */
  return ((PLOTOBJHANDLING*) CreatePlotObjType (PlotObjTypeName,sizeof(PLOTOBJHANDLING)));
}

/****************************************************************************/
/*D
   GetPlotObjHandling - Get PLOTOBJHANDLING

   SYNOPSIS:
   static PLOTOBJHANDLING *GetPlotObjHandling (char *PlotObjHandlingName);

   PARAMETERS:
   .  PlotObjHandlingName - the name

   DESCRIPTION:
   This function gets 'PLOTOBJHANDLING' by name.

   RETURN VALUE:
   PLOTOBJHANDLING *
   .n     pointer to 'PLOTOBJHANDLING'
   .n     NULL if error occured.
   D*/
/****************************************************************************/

static PLOTOBJHANDLING  *GetPlotObjHandling (char *PlotObjHandlingName)
{
  return ((PLOTOBJHANDLING*)GetPlotObjType(PlotObjHandlingName));
}

/****************************************************************************/
/*
   PerspectiveProjection -

   SYNOPSIS:
   static void PerspectiveProjection (COORD *in, COORD_POINT *ScreenPoint);

   PARAMETERS:
   .  in - input vect (3d) to transform
   .  ScreenPoint -

   DESCRIPTION:
   This function project a screen vector perspectively onto the (2d) screen.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void PerspectiveProjection (COORD *in, COORD_POINT *ScreenPoint)
{
  COORD k;

  k = OBS_ViewPlaneDist/(OBS_ViewPlaneDist-in[2]);
  (*ScreenPoint).x = k*in[0] + (1-k)*OBS_PerspCorr[0];
  (*ScreenPoint).y = k*in[1] + (1-k)*OBS_PerspCorr[1];
}

/****************************************************************************/
/*
   NormalProjection - Project a screen vector parallel onto the (2d) screen

   SYNOPSIS:
   static void NormalProjection (COORD *in, COORD_POINT *ScreenPoint);

   PARAMETERS:
   .  in - input vect (3d) to transform
   .  ScreenPoint -

   DESCRIPTION:
   This function projects a screen vector parallel onto the (2d) screen.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void NormalProjection (COORD *in, COORD_POINT *ScreenPoint)
{
  (*ScreenPoint).x = in[0];
  (*ScreenPoint).y = in[1];
}

/****************************************************************************/
/*
   BuildObsTrafo - Build observer transformation

   SYNOPSIS:
   static INT BuildObsTrafo (PICTURE *thePicture);

   PARAMETERS:
   .  thePicture -

   DESCRIPTION:
   This function builds observer transformation.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT BuildObsTrafo (PICTURE *thePicture)
{
  VIEWEDOBJ *theViewedObj;
  PLOTOBJ *thePlotObj;
  COORD VRS_2_PHS[16], PHS_2_VRS[16], VRS_2_SCS[16];
  COORD ZD[3];
  COORD *MP, *XD, *YD;
  INT *LL, *UR;

  theViewedObj = PIC_VO(thePicture);
  if (VO_STATUS(theViewedObj) != ACTIVE) return (1);
  thePlotObj = VO_PO(theViewedObj);
  LL = PIC_GLL(thePicture);
  UR = PIC_GUR(thePicture);
  MP = VO_PMP(theViewedObj);
  XD = VO_PXD(theViewedObj);
  YD = VO_PYD(theViewedObj);

  switch (PO_DIM(thePlotObj))
  {
  case TYPE_2D :

    /* set trafo from phys. space to view reference system */
    VRS_2_PHS[0] = 2.0*XD[0];               VRS_2_PHS[3] = 2.0*YD[0];               VRS_2_PHS[6] = MP[0]-XD[0]-YD[0];
    VRS_2_PHS[1] = 2.0*XD[1];               VRS_2_PHS[4] = 2.0*YD[1];               VRS_2_PHS[7] = MP[1]-XD[1]-YD[1];
    VRS_2_PHS[2] = 0.0;                     VRS_2_PHS[5] = 0.0;                     VRS_2_PHS[8] = 1.0;
    if (M3_Invert(PHS_2_VRS,VRS_2_PHS)) return (1);

    /* set trafo from view reference system to screen system */
    VRS_2_SCS[0] = UR[0]-LL[0];     VRS_2_SCS[3] = 0.0;                     VRS_2_SCS[6] = LL[0];
    VRS_2_SCS[1] = 0.0;                     VRS_2_SCS[4] = UR[1]-LL[1];     VRS_2_SCS[7] = LL[1];
    VRS_2_SCS[2] = 0.0;                     VRS_2_SCS[5] = 0.0;                     VRS_2_SCS[8] = 1.0;

    /* set ObsTrafo and its inverse */
    M3_TIMES_M3(VRS_2_SCS,PHS_2_VRS,ObsTrafo)
    if (M3_Invert(InvObsTrafo,ObsTrafo)) return (1);
    OBS_ProjectProc = NormalProjection;
    break;

  case TYPE_3D :

    /* set trafo from phys. space to view reference system */
    V3_VECTOR_PRODUCT(XD,YD,ZD)
    if (V3_Normalize(ZD)) return (1);
    VRS_2_PHS[0] = XD[0];           VRS_2_PHS[4] = YD[0];           VRS_2_PHS[8] = ZD[0];           VRS_2_PHS[12]= MP[0];
    VRS_2_PHS[1] = XD[1];           VRS_2_PHS[5] = YD[1];           VRS_2_PHS[9] = ZD[1];           VRS_2_PHS[13]= MP[1];
    VRS_2_PHS[2] = XD[2];           VRS_2_PHS[6] = YD[2];           VRS_2_PHS[10]= ZD[2];           VRS_2_PHS[14]= MP[2];
    VRS_2_PHS[3] = 0.0;             VRS_2_PHS[7] = 0.0;             VRS_2_PHS[11]= 0.0;             VRS_2_PHS[15]= 1.0;
    if (M4_Invert(PHS_2_VRS,VRS_2_PHS)) return (1);

    /* set trafo from view reference system to screen system */
    VRS_2_SCS[0] = 0.5*(UR[0]-LL[0]);       VRS_2_SCS[4] = 0.0;                             VRS_2_SCS[8] = 0.0;             VRS_2_SCS[12]= 0.5*(UR[0]+LL[0]);
    VRS_2_SCS[1] = 0.0;                             VRS_2_SCS[5] = 0.5*(UR[1]-LL[1]);       VRS_2_SCS[9] = 0.0;             VRS_2_SCS[13]= 0.5*(UR[1]+LL[1]);
    VRS_2_SCS[2] = 0.0;                             VRS_2_SCS[6] = 0.0;                             VRS_2_SCS[10]= 1.0;             VRS_2_SCS[14]= 0.0;
    VRS_2_SCS[3] = 0.0;                             VRS_2_SCS[7] = 0.0;                             VRS_2_SCS[11]= 0.0;             VRS_2_SCS[15]= 1.0;

    /* set ObsTrafo and its inverse */
    M4_TIMES_M4(VRS_2_SCS,PHS_2_VRS,ObsTrafo)
    if (M4_Invert(InvObsTrafo,ObsTrafo)) return (1);
    if (VO_PERSPECTIVE(theViewedObj) == YES)
    {
      OBS_Perspective = YES;
      OBS_ProjectProc = PerspectiveProjection;
    }
    else
    {
      OBS_Perspective = NO;
      OBS_ProjectProc = NormalProjection;
    }
    OBS_PerspCorr[0] = VRS_2_SCS[12]; OBS_PerspCorr[1] = VRS_2_SCS[13];
    V3_SUBTRACT(VO_VP(theViewedObj),VO_VT(theViewedObj),OBS_ViewDirection)
    V3_EUKLIDNORM(OBS_ViewDirection,OBS_ViewPlaneDist)
    if (OBS_ViewPlaneDist<SMALL_C) return (1);
    break;

  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*
   BuildCutTrafo - Build cut transformation

   SYNOPSIS:
   static INT BuildCutTrafo (CUT *theCut, COORD *theViewDir);

   PARAMETERS:
   .  theCut -
   .  theViewDir -

   DESCRIPTION:
   This function builds cut transformation.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT BuildCutTrafo (CUT *theCut, COORD *theViewDir)
{
  COORD XD[3], YD[3], ZD[3];
  COORD *PN, *PP;
  COORD scpr;

  CUT_CutExisting = 0;
  if (theCut == NULL) return (1);
  if (CUT_STATUS(theCut)!=ACTIVE) return (0);
  CUT_CutExisting = 1;
  PN = CUT_PN(theCut);
  PP = CUT_PP(theCut);

  /* set trafo from phys. space to view reference system */
  V3_COPY(PN,ZD)
  if (V3_Normalize(ZD)) return (1);
  V3_Orthogonalize(ex,ZD,XD);
  if (V3_Normalize(XD))
  {
    V3_Orthogonalize(ey,ZD,XD);
    if (V3_Normalize(XD)) return (1);
  }
  V3_VECTOR_PRODUCT(ZD,XD,YD)
  if (V3_Normalize(YD)) return (1);
  InvCutTrafo[0] = XD[0];         InvCutTrafo[4] = YD[0];         InvCutTrafo[8] = ZD[0];         InvCutTrafo[12]= PP[0];
  InvCutTrafo[1] = XD[1];         InvCutTrafo[5] = YD[1];         InvCutTrafo[9] = ZD[1];         InvCutTrafo[13]= PP[1];
  InvCutTrafo[2] = XD[2];         InvCutTrafo[6] = YD[2];         InvCutTrafo[10]= ZD[2];         InvCutTrafo[14]= PP[2];
  InvCutTrafo[3] = 0.0;           InvCutTrafo[7] = 0.0;           InvCutTrafo[11]= 0.0;           InvCutTrafo[15]= 1.0;
  if (M4_Invert(CutTrafo,InvCutTrafo)) return (1);

  /* is cut plane at front? */
  V3_SCALAR_PRODUCT(theViewDir,PN,scpr)
  CUT_CutAtFront = 0;
  if (scpr>0.0) CUT_CutAtFront = 1;
  V3_COPY(ZD,CUT_CutNormal)

  return (0);
}

/****************************************************************************/
/*
   EW_GetNextElement_hor_fw_up - Return next Element in horizontal, forward, upward direction

   SYNOPSIS:
   static ELEMENT *EW_GetNextElement_hor_fw_up (ELEMENT *theElement);

   PARAMETERS:
   .  theElement - Element for which the successor is to be determined

   DESCRIPTION:
   This function returns next Element in horizontal, forward, upward direction.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to
   .n     NULL if the argument is the last Element.
 */
/****************************************************************************/

static ELEMENT *EW_GetNextElement_hor_fw_up (ELEMENT *theElement)
{
  INT currLevel;

  currLevel = LEVEL(theElement);

  do
  {
    if (SUCCE(theElement)==NULL)
    {
      if (++currLevel>GElem_toLevel)
        return (NULL);
      else
        theElement = FIRSTELEMENT(GRID_ON_LEVEL(GElem_MG,currLevel));
    }
    else
      theElement = SUCCE(theElement);
  }
  while (!USED(theElement));

  return (theElement);
}

/****************************************************************************/
/*
   EW_GetNextElement_hor_fw_down - Return next Element in horizontal, forward, downward direction

   SYNOPSIS:
   static ELEMENT *EW_GetNextElement_hor_fw_down (ELEMENT *theElement);

   PARAMETERS:
   .  theElement - Element for which the successor is to be determined

   DESCRIPTION:
   This function return next Element in horizontal, forward, downward direction.

   RETURN VALUE:
   ELEMENT *
   .n      pointer to
   .n      NULL if the argument is the last Element.
 */
/****************************************************************************/

static ELEMENT *EW_GetNextElement_hor_fw_down (ELEMENT *theElement)
{
  INT currLevel;

  currLevel = LEVEL(theElement);

  do
  {
    if (SUCCE(theElement)==NULL)
    {
      if (--currLevel<GElem_fromLevel)
        return (NULL);
      else
        theElement = FIRSTELEMENT(GRID_ON_LEVEL(GElem_MG,currLevel));
    }
    else
      theElement = SUCCE(theElement);
  }
  while (!USED(theElement));

  return (theElement);
}

/****************************************************************************/
/*
   EW_GetNextElement_hor_bw_up - Return next Element in horizontal, forward, upward direction

   SYNOPSIS:
   static ELEMENT *EW_GetNextElement_hor_bw_up (ELEMENT *theElement);

   PARAMETERS:
   .  theElement - Element for which the successor is to be determined

   DESCRIPTION:
   This function return next Element in horizontal, forward, upward direction.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to
   .n     NULL if the argument is the last Element.
 */
/****************************************************************************/

static ELEMENT *EW_GetNextElement_hor_bw_up (ELEMENT *theElement)
{
  INT currLevel;

  currLevel = LEVEL(theElement);

  do
  {
    if (PREDE(theElement)==NULL)
    {
      if (++currLevel>GElem_toLevel)
        return (NULL);
      else
        theElement = LASTELEMENT(GRID_ON_LEVEL(GElem_MG,currLevel));
    }
    else
      theElement = PREDE(theElement);
  }
  while (!USED(theElement));

  return (theElement);
}

/****************************************************************************/
/*
   EW_GetNextElement_hor_bw_down - return next Element in horizontal, forward, downward direction

   SYNOPSIS:
   static ELEMENT *EW_GetNextElement_hor_bw_down (ELEMENT *theElement)

   PARAMETERS:
   .  theElement - Element for which the successor is to be determined

   DESCRIPTION:
   This function return next Element in horizontal, forward, downward direction.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to ELEMENT *
   .n     NULL if the argument is the last Element.
 */
/****************************************************************************/

static ELEMENT *EW_GetNextElement_hor_bw_down (ELEMENT *theElement)
{
  INT currLevel;

  currLevel = LEVEL(theElement);

  do
  {
    if (PREDE(theElement)==NULL)
    {
      if (--currLevel<GElem_fromLevel)
        return (NULL);
      else
        theElement = LASTELEMENT(GRID_ON_LEVEL(GElem_MG,currLevel));
    }
    else
      theElement = PREDE(theElement);
  }
  while (!USED(theElement));

  return (theElement);
}
/****************************************************************************/
/*
   EW_GetNextElement_vert_fw_up - Return next Element in vertical, forward, upward direction

   SYNOPSIS:
   static ELEMENT *EW_GetNextElement_vert_fw_up (ELEMENT *theElement);

   PARAMETERS:
   .  theElement - Element for which the successor is to be determined

   DESCRIPTION:
   This function return next Element in vertical, forward, upward direction.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to ELEMENT *
   .n     NULL if the argument is the last Element.
 */
/****************************************************************************/

static ELEMENT *EW_GetNextElement_vert_fw_up (ELEMENT *theElement)
{
  do
  {
    /* something above? */
    if ((LEVEL(theElement)<GElem_toLevel) && (SON(theElement,0)!=NULL))
    {
      /* first go up */
      theElement = SON(theElement,0);

      /* go back while father not left */
      while ((PREDE(theElement)!=NULL) && (EFATHER(PREDE(theElement))==EFATHER(theElement)))
        theElement = PREDE(theElement);
    }
    else
    {
      /* descent while all sons of the father are processed */
      while ((LEVEL(theElement)>GElem_fromLevel)
             && ((SUCCE(theElement)==NULL) || (EFATHER(SUCCE(theElement))!=EFATHER(theElement))))
        theElement = EFATHER(theElement);

      theElement = SUCCE(theElement);
      if (theElement==NULL)
        return (NULL);
    }
  }
  while (!USED(theElement));

  return (theElement);
}
/****************************************************************************/
/*
   EW_GetNextElement_vert_bw_up - Return next Element in vertical, backward, upward direction

   SYNOPSIS:
   static ELEMENT *EW_GetNextElement_vert_bw_up (ELEMENT *theElement);

   PARAMETERS:
   .  theElement - Element for which the successor is to be determined

   DESCRIPTION:
   This function returns next Element in vertical, backward, upward direction.

   RETURN VALUE:
   ELEMENT *
   .n      pointer to ELEMENT *
   .n      NULL if the argument is the last Element.
 */
/****************************************************************************/

static ELEMENT *EW_GetNextElement_vert_bw_up (ELEMENT *theElement)
{
  do
  {
    if ((PREDE(theElement)!=NULL) && (EFATHER(PREDE(theElement))==EFATHER(theElement)))
    {
      /* predecessor has the same father */
      theElement = PREDE(theElement);

      /* go up as long as possible */
      while ((LEVEL(theElement)<GElem_toLevel) && (SON(theElement,0)!=NULL))
      {
        theElement = SON(theElement,0);

        /* go forward while father not left */
        while ((SUCCE(theElement)!=NULL) && (EFATHER(SUCCE(theElement))==EFATHER(theElement)))
          theElement = SUCCE(theElement);
      }
    }
    else
    {
      /* descent to father */
      if (LEVEL(theElement)>GElem_fromLevel)
        theElement = EFATHER(theElement);
      else
        return (NULL);
    }
  }
  while (!USED(theElement));

  return (theElement);
}

/****************************************************************************/
/*																			*/
/* Function:  EW_GetFirstElement_hor_fw_up									*/
/*																			*/
/* Purpose:   return first Element in horizontal, forward, upward direction */
/*																			*/
/* Input:	  the corresponding multigrid, from level, to level                     */
/*																			*/
/* Output:	  first Element, NULL if there is none to return				*/
/*																			*/
/****************************************************************************/

static ELEMENT *EW_GetFirstElement_hor_fw_up (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GElem_MG                = theMG;
  GElem_fromLevel = fromLevel;
  GElem_toLevel   = toLevel;

  theElement = FIRSTELEMENT(GRID_ON_LEVEL(GElem_MG,fromLevel));

  if (theElement==NULL)
    return (NULL);
  else if (USED(theElement))
    return (theElement);
  else
    return (EW_GetNextElement_hor_fw_up(theElement));
}

/****************************************************************************/
/*
   EW_GetFirstElement_hor_fw_down - return first Element in horizontal, forward, downward direction

   SYNOPSIS:
   static ELEMENT *EW_GetFirstElement_hor_fw_down (MULTIGRID *theMG, INT fromLevel,
   INT toLevel);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function returns first Element in horizontal, forward, downward direction.

   RETURN VALUE:
   ELEMENT *
   .n      pointer to first Element
   .n      NULL if there is none to return.
 */
/****************************************************************************/

static ELEMENT *EW_GetFirstElement_hor_fw_down (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GElem_MG                = theMG;
  GElem_fromLevel = fromLevel;
  GElem_toLevel   = toLevel;

  theElement = FIRSTELEMENT(GRID_ON_LEVEL(GElem_MG,toLevel));

  if (theElement==NULL)
    return (NULL);
  else if (USED(theElement))
    return (theElement);
  else
    return (EW_GetNextElement_hor_fw_down(theElement));
}
/****************************************************************************/
/*
   EW_GetFirstElement_hor_bw_up	- Return first Element in horizontal, backward, upward direction

   SYNOPSIS:
   static ELEMENT *EW_GetFirstElement_hor_bw_up (MULTIGRID *theMG,
   INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG - the pointer to multigrid
   .  from Level -
   .  toLevel -

   DESCRIPTION:
   This function returns first Element in horizontal, backward, upward direction.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to first Element
   .n     NULL if there is none to return.
 */
/****************************************************************************/

static ELEMENT *EW_GetFirstElement_hor_bw_up (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GElem_MG                = theMG;
  GElem_fromLevel = fromLevel;
  GElem_toLevel   = toLevel;

  theElement = LASTELEMENT(GRID_ON_LEVEL(GElem_MG,fromLevel));

  if (theElement==NULL)
    return (NULL);
  else if (USED(theElement))
    return (theElement);
  else
    return (EW_GetNextElement_hor_bw_up(theElement));
}

/****************************************************************************/
/*
   EW_GetFirstElement_hor_bw_down - Return first Element in horizontal, backward, downward direction

   SYNOPSIS:
   static ELEMENT *EW_GetFirstElement_hor_bw_down (MULTIGRID *theMG,
   INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function returns first Element in horizontal, backward, downward direction.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to first Element
   .n     NULL if there is none to return.
 */
/****************************************************************************/

static ELEMENT *EW_GetFirstElement_hor_bw_down (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GElem_MG                = theMG;
  GElem_fromLevel = fromLevel;
  GElem_toLevel   = toLevel;

  theElement = LASTELEMENT(GRID_ON_LEVEL(GElem_MG,toLevel));

  if (theElement==NULL)
    return (NULL);
  else if (USED(theElement))
    return (theElement);
  else
    return (EW_GetNextElement_hor_bw_down(theElement));
}

/****************************************************************************/
/*
   EW_GetFirstElement_vert_fw_up - Return first Element in vertical, forward, upward direction

   SYNOPSIS:
   static ELEMENT *EW_GetFirstElement_vert_fw_up (MULTIGRID *theMG, INT fromLevel,
   INT toLevel)

   PARAMETERS:
   .  theMG -
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function  returns first Element in vertical, forward, upward direction.

   RETURN VALUE:
   ELEMENT *
   .n      pointer to first Element
   .n      NULL if there is none to return.
 */
/****************************************************************************/

static ELEMENT *EW_GetFirstElement_vert_fw_up (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GElem_MG                = theMG;
  GElem_fromLevel = fromLevel;
  GElem_toLevel   = toLevel;

  theElement = FIRSTELEMENT(GRID_ON_LEVEL(GElem_MG,fromLevel));

  if (theElement==NULL)
    return (NULL);
  else if (USED(theElement))
    return (theElement);
  else
    return (EW_GetNextElement_vert_fw_up(theElement));
}

/****************************************************************************/
/*
   EW_GetFirstElement_vert_bw_up -  Return first Element in vertical, backward, upward direction

   SYNOPSIS:
   static ELEMENT *EW_GetFirstElement_vert_bw_up (MULTIGRID *theMG,
   INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG -
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function returns first Element in vertical, backward, upward direction.

   RETURN VALUE:
   ELEMENT *
   .n      pointer to first Element
   .n      NULL if there is none to return.
 */
/****************************************************************************/

static ELEMENT *EW_GetFirstElement_vert_bw_up (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GElem_MG                = theMG;
  GElem_fromLevel = fromLevel;
  GElem_toLevel   = toLevel;

  theElement = LASTELEMENT(GRID_ON_LEVEL(GElem_MG,fromLevel));

  if (theElement==NULL)
    return (NULL);

  /* go up as long as possible */
  while ((LEVEL(theElement)<GElem_toLevel) && (SON(theElement,0)!=NULL))
  {
    theElement = SON(theElement,0);

    /* go foreward while father not left */
    while ((SUCCE(theElement)!=NULL) && (EFATHER(SUCCE(theElement))==EFATHER(theElement)))
      theElement = SUCCE(theElement);
  }

  if (USED(theElement))
    return (theElement);
  else
    return (EW_GetNextElement_vert_bw_up(theElement));
}

/****************************************************************************/
/*
   MarkElements_MGS - Mark elements on surface of multigrid

   SYNOPSIS:
   static INT MarkElements_MGS (MULTIGRID *theMG, INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function marks elements on surface of multigrid.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MarkElements_MGS (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  for (i=fromLevel; i<toLevel; i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      if (NSONS(theElement)==0)
        SETUSED(theElement,1);
      else
        SETUSED(theElement,0);

  for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
    SETUSED(theElement,1);

  return (0);
}

/****************************************************************************/
/*
   MarkElements_MGS - Mark elements on surface of multigrid

   SYNOPSIS:
   static INT MarkElements_ID (MULTIGRID *theMG, INT id);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  id -

   DESCRIPTION:
   This function marks elements on surface of multigrid.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT MarkElements_ID (MULTIGRID *theMG, INT id)
{
  ELEMENT *theElement;
  INT i;

  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      if (ID(theElement)==id)
        SETUSED(theElement,1);
      else
        SETUSED(theElement,0);

  return (0);
}

/****************************************************************************/
/*
   MarkElements2D - Mark elements on surface of multigrid specified by EE2D_Elem2P

   SYNOPSIS:
   static INT MarkElements2D (MULTIGRID *theMG, INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG -
   .  fromlevel -
   .  toLevel -

   DESCRIPTION:
   This function marks elements on surface of multigrid specified by EE2D_Elem2P.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MarkElements2D (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  if (EE2D_Elem2Plot[PLOT_ALL])
  {
    for (i=fromLevel; i<toLevel; i++)
      for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
        if (NSONS(theElement)==0 && EE2D_Elem2Plot[ECLASS(theElement)])
          SETUSED(theElement,1);
        else
          SETUSED(theElement,0);
  }
  else
  {
    for (i=fromLevel; i<toLevel; i++)
      for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
        SETUSED(theElement,0);
  }

  for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
    if (EE2D_Elem2Plot[ECLASS(theElement)])
      SETUSED(theElement,1);
    else
      SETUSED(theElement,0);

  for (i=toLevel+1; i<=TOPLEVEL(theMG); i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      SETUSED(theElement,0);

  return (0);
}

/****************************************************************************/
/*
   MarkElements3D - Mark elements on surface of multigrid specified by EE2D_Elem2P

   SYNOPSIS:
   static INT MarkElements3D (MULTIGRID *theMG, INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG -
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function marks elements on surface of multigrid specified by EE2D_Elem2P.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
 */
/****************************************************************************/

static INT MarkElements3D (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  if (EE3D_Elem2Plot[PLOT_ALL])
  {
    for (i=fromLevel; i<toLevel; i++)
      for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
        if (NSONS(theElement)==0 && EE3D_Elem2Plot[ECLASS(theElement)])
          SETUSED(theElement,1);
        else
          SETUSED(theElement,0);
  }
  else
  {
    for (i=fromLevel; i<toLevel; i++)
      for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
        SETUSED(theElement,0);
  }

  for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
    if (EE3D_Elem2Plot[ECLASS(theElement)])
      SETUSED(theElement,1);
    else
      SETUSED(theElement,0);

  for (i=toLevel+1; i<=TOPLEVEL(theMG); i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      SETUSED(theElement,0);

  return (0);
}

/****************************************************************************/
/*
   MarkElements_Level -  Mark elements on level of multigrid

   SYNOPSIS:
   static INT MarkElements_Level (MULTIGRID *theMG, INT level);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  level -

   DESCRIPTION:
   This function marks elements on level of multigrid.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 when error occured.
 */
/****************************************************************************/

static INT MarkElements_Level (MULTIGRID *theMG, INT level)
{
  ELEMENT *theElement;
  INT i;

  level = MAX(level,0);
  level = MIN(level,CURRENTLEVEL(theMG));

  for (i=0; i<level; i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      SETUSED(theElement,0);

  for (theElement=FIRSTELEMENT(theMG->grids[level]); theElement!=NULL; theElement=SUCCE(theElement))
    SETUSED(theElement,1);

  return (0);
}

/****************************************************************************/
/*
   MarkElements_MGS_Bnd - Mark bnd elements on surface of multigrid

   SYNOPSIS:
   static INT MarkElements_MGS_Bnd (MULTIGRID *theMG, INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function marks bnd elements on surface of multigrid.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MarkElements_MGS_Bnd (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  for (i=fromLevel; i<toLevel; i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      if (NSONS(theElement)==0 && OBJT(theElement)==BEOBJ)
        SETUSED(theElement,1);
      else
        SETUSED(theElement,0);

  for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
    if (OBJT(theElement)==BEOBJ)
      SETUSED(theElement,1);
    else
      SETUSED(theElement,0);

  return (0);
}

/****************************************************************************/
/*
   MarkElements_MGS_Bnd_Cut - Mark elements on surface of multigrid lying on bnd or cut

   SYNOPSIS:
   static INT MarkElements_MGS_Bnd_Cut (MULTIGRID *theMG, INT fromLevel,
   INT toLevel);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function marks elements on surface of multigrid lying on bnd or cut.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MarkElements_MGS_Bnd_Cut (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  for (i=fromLevel; i<toLevel; i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      if (NSONS(theElement)==0 && (OBJT(theElement)==BEOBJ || CUTMODE(theElement)==CM_INTERSECT))
        SETUSED(theElement,1);
      else
        SETUSED(theElement,0);

  for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
    if (OBJT(theElement)==BEOBJ || CUTMODE(theElement)==CM_INTERSECT)
      SETUSED(theElement,1);
    else
      SETUSED(theElement,0);

  return (0);
}

/****************************************************************************/
/*
   MarkElements_MGS_Bnd_and_Cut	- Mark elements on surface of multigrid lying on bnd or cut

   SYNOPSIS:
   static INT MarkElements_MGS_Bnd_and_Cut (MULTIGRID *theMG, INT fromLevel,
   INT toLevel);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function marks elements on surface of multigrid lying on bnd or cut.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MarkElements_MGS_Bnd_and_Cut (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  for (i=fromLevel; i<toLevel; i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      if (NSONS(theElement)==0 && OBJT(theElement)==BEOBJ && CUTMODE(theElement)==CM_INTERSECT)
        SETUSED(theElement,1);
      else
        SETUSED(theElement,0);

  for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
    if (OBJT(theElement)==BEOBJ && CUTMODE(theElement)==CM_INTERSECT)
      SETUSED(theElement,1);
    else
      SETUSED(theElement,0);

  return (0);
}

/****************************************************************************/
/*
   MarkElements_MGS_Cut - Mark elements on surface of multigrid being intersected

   SYNOPSIS:
   static INT MarkElements_MGS_Cut (MULTIGRID *theMG, INT fromLevel, INT toLevel);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  fromLevel -
   .  toLevel -

   DESCRIPTION:
   This function marks elements on surface of multigrid being intersected.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MarkElements_MGS_Cut (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  ELEMENT *theElement;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  for (i=fromLevel; i<toLevel; i++)
    for (theElement=FIRSTELEMENT(theMG->grids[i]); theElement!=NULL; theElement=SUCCE(theElement))
      if (NSONS(theElement)==0 && CUTMODE(theElement)==CM_INTERSECT)
        SETUSED(theElement,1);
      else
        SETUSED(theElement,0);

  for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
    if (CUTMODE(theElement)==CM_INTERSECT)
      SETUSED(theElement,1);
    else
      SETUSED(theElement,0);

  return (0);
}

/****************************************************************************/
/*
   EW_GetFirstElement_vert_fw_up_Proc  - Get the GetFirstElementProc

   SYNOPSIS:
   static EW_GetFirstElementProcPtr EW_GetFirstElement_vert_fw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetFirstElementProc.

   RETURN VALUE:
   EW_GetFirstElementProcPtr
   .n    pointer
   .n    NULL
 */
/****************************************************************************/

static EW_GetFirstElementProcPtr EW_GetFirstElement_vert_fw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetFirstElement_vert_fw_up);
}
/****************************************************************************/
/*
   EW_GetFirstElement_vert_bw_up_Proc - Get the GetFirstElementProc

   SYNOPSIS:
   static EW_GetFirstElementProcPtr EW_GetFirstElement_vert_bw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetFirstElementProc.

   RETURN VALUE:
   EW_GetFirstElementProcPtr
   .n      pointer to
   .n      Null if error occured.
 */
/****************************************************************************/
static EW_GetFirstElementProcPtr EW_GetFirstElement_vert_bw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetFirstElement_vert_bw_up);
}
/****************************************************************************/
/*
   EW_GetFirstElement_hor_fw_up_Proc - Get the GetFirstElementProc

   SYNOPSIS:
   static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_fw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetFirstElementProc.

   RETURN VALUE:
   EW_GetFirstElementProcPtr
   .n      pointer to
   .n      Null if error occured.
 */
/****************************************************************************/
static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_fw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetFirstElement_hor_fw_up);
}
/****************************************************************************/
/*
   EW_GetFirstElement_hor_fw_down_Proc - Get the GetFirstElementProc

   SYNOPSIS:
   static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_fw_down_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetFirstElementProc.

   RETURN VALUE:
   EW_GetFirstElementProcPtr
   .n      pointer to
   .n      Null if error occured
 */
/****************************************************************************/
static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_fw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetFirstElement_hor_fw_down);
}
/****************************************************************************/
/*
   EW_GetFirstElement_hor_bw_up_Proc  - Get the GetFirstElementProc

   SYNOPSIS:
   static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_bw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetFirstElementProc.

   RETURN VALUE:
   EW_GetFirstElementProcPtr
   .n      pointer to
   .n      Null if error occured.
 */
/****************************************************************************/
static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_bw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetFirstElement_hor_bw_up);
}
/****************************************************************************/
/*
   EW_GetFirstElement_hor_bw_down_Proc  - Get the GetFirstElementProc

   SYNOPSIS:
   static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_bw_down_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetFirstElementProc.

   RETURN VALUE:
   EW_GetFirstElementProcPtr
   .n       pointer to
   .n       Null if error occured.
 */
/****************************************************************************/
static EW_GetFirstElementProcPtr EW_GetFirstElement_hor_bw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetFirstElement_hor_bw_down);
}

/****************************************************************************/
/*
   EW_GetNextElement_vert_fw_up_Proc - Get the GetNextElementProc

   SYNOPSIS:
   static EW_GetNextElementProcPtr EW_GetNextElement_vert_fw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetNextElementProc.

   RETURN VALUE:
   EW_GetNextElementProcPtr
   .n      pointer to
   .n      NULL, when error occured.
 */
/****************************************************************************/

static EW_GetNextElementProcPtr EW_GetNextElement_vert_fw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetNextElement_vert_fw_up);
}
/****************************************************************************/
/*
   EW_GetNextElement_vert_bw_up_Proc - Get the GetNextElementProc

   SYNOPSIS:
   static EW_GetNextElementProcPtr EW_GetNextElement_vert_bw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetNextElementProc.

   RETURN VALUE:
   EW_GetNextElementProcPtr
   .n      pointer to
   .n      NULL, when error occured.
 */
/****************************************************************************/
static EW_GetNextElementProcPtr EW_GetNextElement_vert_bw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetNextElement_vert_bw_up);
}
/****************************************************************************/
/*
   EW_GetNextElement_hor_fw_up_Proc - Get the GetNextElementProc

   SYNOPSIS:
   static EW_GetNextElementProcPtr EW_GetNextElement_hor_fw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetNextElementProc.

   RETURN VALUE:
   EW_GetNextElementProcPtr
   .n      pointer to
   .n      NULL if error occured.
 */
/****************************************************************************/
static EW_GetNextElementProcPtr EW_GetNextElement_hor_fw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetNextElement_hor_fw_up);
}
/****************************************************************************/
/*
   EW_GetNextElement_hor_fw_down_Proc - Get the GetNextElementProc

   SYNOPSIS:
   static EW_GetNextElementProcPtr EW_GetNextElement_hor_fw_down_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetNextElementProc.

   RETURN VALUE:
   EW_GetNextElementProcPtr
   .n     pointer to
   .n     NULL if error occured.
 */
/****************************************************************************/
static EW_GetNextElementProcPtr EW_GetNextElement_hor_fw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetNextElement_hor_fw_down);
}

/****************************************************************************/
/*
   EW_GetNextElement_hor_bw_up_Proc - Get the GetNextElementProc

   SYNOPSIS:
   static EW_GetNextElementProcPtr EW_GetNextElement_hor_bw_up_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetNextElementProc.

   RETURN VALUE:
   EW_GetNextElementProcPtr
   .n     pointer to
   .n     NULL if error occured.
 */
/****************************************************************************/
static EW_GetNextElementProcPtr EW_GetNextElement_hor_bw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetNextElement_hor_bw_up);
}

/****************************************************************************/
/*
   EW_GetNextElement_hor_bw_down_Proc - Get the GetNextElementProc

   SYNOPSIS:
   static EW_GetNextElementProcPtr EW_GetNextElement_hor_bw_down_Proc
   (VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theViewedObj -

   DESCRIPTION:
   This function gets the GetNextElementProc.

   RETURN VALUE:
   EW_GetNextElementProcPtr
   .n      pointer to
   .n      NULL if error occured.
 */
/****************************************************************************/

static EW_GetNextElementProcPtr EW_GetNextElement_hor_bw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (EW_GetNextElement_hor_bw_down);
}

/****************************************************************************/
/*																			*/
/* Function:  MarkNodes_MGS                                                                                             */
/*																			*/
/* Purpose:   mark nodes on surface of multigrid							*/
/*																			*/
/* Input:	  MULTIGRID *theMG												*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT MarkNodes_MGS (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  NODE *theNode;
  INT i;

  fromLevel = MAX(fromLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  for (i=fromLevel; i<toLevel; i++)
    for (theNode=FIRSTNODE(theMG->grids[i]); theNode!=NULL; theNode=SUCCN(theNode))
      if (SONNODE(theNode)==NULL)
        SETUSED(theNode,1);
      else
        SETUSED(theNode,0);

  for (theNode=FIRSTNODE(theMG->grids[toLevel]); theNode!=NULL; theNode=SUCCN(theNode))
    SETUSED(theNode,1);

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  MarkNodes_OfMarkedElem										*/
/*																			*/
/* Purpose:   mark nodes on level of multigrid								*/
/*																			*/
/* Input:	  MULTIGRID *theMG												*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

#define MARKMODE_ALL            0
#define MARKMODE_COPY           1
#define MARKMODE_IRREG          2
#define MARKMODE_REG            3

static INT MarkNodes_OfMarkedElem (MULTIGRID *theMG, INT fromLevel, INT toLevel, INT mode)
{
  NODE *theNode;
  ELEMENT *theElement;
  INT i, j, limit;

  fromLevel = MAX(fromLevel,0);
  fromLevel = MIN(fromLevel,CURRENTLEVEL(theMG));
  toLevel = MAX(toLevel,0);
  toLevel = MIN(toLevel,CURRENTLEVEL(theMG));

  /* set marks from elements */
  limit = 0;
  switch (mode)
  {
  case MARKMODE_ALL :
    for (i=fromLevel; i<=toLevel; i++)
      for (theNode=FIRSTNODE(theMG->grids[i]); theNode!=NULL; theNode=SUCCN(theNode))
        SETUSED(theNode,1);
    break;
  case MARKMODE_REG :
    limit = MAX(limit,REGULAR_CLASS);
  case MARKMODE_IRREG :
    limit = MAX(limit,GREEN);
  case MARKMODE_COPY :
    limit = MAX(limit,YELLOW);
    for (i=fromLevel; i<=toLevel; i++)
      for (theNode=FIRSTNODE(theMG->grids[i]); theNode!=NULL; theNode=SUCCN(theNode))
        SETUSED(theNode,0);
    for (theElement=FIRSTELEMENT(theMG->grids[toLevel]); theElement!=NULL; theElement=SUCCE(theElement))
      if (ECLASS(theElement)>=limit)
        for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
          SETUSED(CORNER(theElement,j),1);
    break;
  default :
    return (1);
  }

  /* only marks on surface count */
  for (i=fromLevel; i<toLevel; i++)
    for (theNode=FIRSTNODE(theMG->grids[i]); theNode!=NULL; theNode=SUCCN(theNode))
      if (SONNODE(theNode)!=NULL)
        SETUSED(theNode,0);

  /* skip upper levels */
  for (i=toLevel+1; i<=TOPLEVEL(theMG); i++)
    for (theNode=FIRSTNODE(theMG->grids[i]); theNode!=NULL; theNode=SUCCN(theNode))
      SETUSED(theNode,0);

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetNextNode_hor_fw_up										*/
/*																			*/
/* Purpose:   return next Node in horizontal, forward, upward direction         */
/*																			*/
/* Input:	  Node for which the successor is to be determined				*/
/*																			*/
/* Output:	  next Node, NULL if the argument is the last Node				*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetNextNode_hor_fw_up (NODE *theNode)
{
  INT currLevel;

  currLevel = LEVEL(theNode);

  do
  {
    if (SUCCN(theNode)==NULL)
    {
      if (++currLevel>GNode_toLevel)
        return (NULL);
      else
        theNode = FIRSTNODE(GRID_ON_LEVEL(GNode_MG,currLevel));
    }
    else
      theNode = SUCCN(theNode);
  }
  while (!USED(theNode));

  return (theNode);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetNextNode_hor_fw_down									*/
/*																			*/
/* Purpose:   return next Node in horizontal, forward, downward direction	*/
/*																			*/
/* Input:	  Node for which the successor is to be determined				*/
/*																			*/
/* Output:	  next Node, NULL if the argument is the last Node				*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetNextNode_hor_fw_down (NODE *theNode)
{
  INT currLevel;

  currLevel = LEVEL(theNode);

  do
  {
    if (SUCCN(theNode)==NULL)
    {
      if (--currLevel<GNode_fromLevel)
        return (NULL);
      else
        theNode = FIRSTNODE(GRID_ON_LEVEL(GNode_MG,currLevel));
    }
    else
      theNode = SUCCN(theNode);
  }
  while (!USED(theNode));

  return (theNode);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetNextNode_hor_bw_up										*/
/*																			*/
/* Purpose:   return next Node in horizontal, forward, upward direction         */
/*																			*/
/* Input:	  Node for which the successor is to be determined				*/
/*																			*/
/* Output:	  next Node, NULL if the argument is the last Node				*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetNextNode_hor_bw_up (NODE *theNode)
{
  INT currLevel;

  currLevel = LEVEL(theNode);

  do
  {
    if (PREDN(theNode)==NULL)
    {
      if (++currLevel>GNode_toLevel)
        return (NULL);
      else
        theNode = LASTNODE(GRID_ON_LEVEL(GNode_MG,currLevel));
    }
    else
      theNode = PREDN(theNode);
  }
  while (!USED(theNode));

  return (theNode);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetNextNode_hor_bw_down									*/
/*																			*/
/* Purpose:   return next Node in horizontal, forward, downward direction	*/
/*																			*/
/* Input:	  Node for which the successor is to be determined				*/
/*																			*/
/* Output:	  next Node, NULL if the argument is the last Node				*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetNextNode_hor_bw_down (NODE *theNode)
{
  INT currLevel;

  currLevel = LEVEL(theNode);

  do
  {
    if (PREDN(theNode)==NULL)
    {
      if (--currLevel<GNode_fromLevel)
        return (NULL);
      else
        theNode = LASTNODE(GRID_ON_LEVEL(GNode_MG,currLevel));
    }
    else
      theNode = PREDN(theNode);
  }
  while (!USED(theNode));

  return (theNode);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetNextNode_leave_fw										*/
/*																			*/
/* Purpose:   return next leave-Node in forward direction					*/
/*																			*/
/* Input:	  Node for which the successor is to be determined				*/
/*																			*/
/* Output:	  next Node, NULL if the argument is the last Node				*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetNextNode_leave_fw (NODE *theNode)
{
  INT currLevel;
  VERTEX *theVertex;

  do
  {
    if ((theVertex=SUCCV(MYVERTEX(theNode)))==NULL)
    {
      currLevel = LEVEL(MYVERTEX(theNode));
      if (++currLevel>GNode_toLevel)
        return (NULL);
      else
        theNode = TOPNODE(FIRSTVERTEX(GRID_ON_LEVEL(GNode_MG,currLevel)));
    }
    else
      theNode=TOPNODE(theVertex);
    assert(theNode!=NULL);
    while (LEVEL(theNode)>GNode_toLevel && theNode!=NULL) theNode=NFATHER(theNode);
    if (theNode==NULL) return (NULL);
  }
  while (!USED(theNode));

  return (theNode);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetNextNode_leave_bw										*/
/*																			*/
/* Purpose:   return next leave-Node in backward direction					*/
/*																			*/
/* Input:	  Node for which the successor is to be determined				*/
/*																			*/
/* Output:	  next Node, NULL if the argument is the last Node				*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetNextNode_leave_bw (NODE *theNode)
{
  INT currLevel;
  VERTEX *theVertex;

  do
  {
    if ((theVertex=PREDV(MYVERTEX(theNode)))==NULL)
    {
      currLevel = LEVEL(MYVERTEX(theNode));
      if (--currLevel<GNode_fromLevel)
        return (NULL);
      else
        theNode = TOPNODE(LASTVERTEX(GRID_ON_LEVEL(GNode_MG,currLevel)));
    }
    else
      theNode = TOPNODE(theVertex);
    while (LEVEL(theNode)>GNode_toLevel && theNode!=NULL) theNode=NFATHER(theNode);
    if (theNode==NULL) return (NULL);
  }
  while (!USED(theNode));

  return (theNode);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetFirstNode_hor_fw_up                                                                         */
/*																			*/
/* Purpose:   return first Node in horizontal, forward, upward direction	*/
/*																			*/
/* Input:	  the corresponding multigrid, from level, to level                     */
/*																			*/
/* Output:	  first Node, NULL if there is none to return					*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetFirstNode_hor_fw_up (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  NODE *theNode;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GNode_MG                = theMG;
  GNode_fromLevel = fromLevel;
  GNode_toLevel   = toLevel;

  theNode = FIRSTNODE(GRID_ON_LEVEL(GNode_MG,fromLevel));

  if (theNode==NULL)
    return (NULL);
  else if (USED(theNode))
    return (theNode);
  else
    return (NW_GetNextNode_hor_fw_up(theNode));
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetFirstNode_hor_fw_down									*/
/*																			*/
/* Purpose:   return first Node in horizontal, forward, downward directio	*/
/*																			*/
/* Input:	  the corresponding multigrid, from level, to level                     */
/*																			*/
/* Output:	  first Node, NULL if there is none to return					*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetFirstNode_hor_fw_down (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  NODE *theNode;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GNode_MG                = theMG;
  GNode_fromLevel = fromLevel;
  GNode_toLevel   = toLevel;

  theNode = FIRSTNODE(GRID_ON_LEVEL(GNode_MG,toLevel));

  if (theNode==NULL)
    return (NULL);
  else if (USED(theNode))
    return (theNode);
  else
    return (NW_GetNextNode_hor_fw_down(theNode));
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetFirstNode_hor_bw_up                                                                         */
/*																			*/
/* Purpose:   return first Node in horizontal, backward, upward direction	*/
/*																			*/
/* Input:	  the corresponding multigrid, from level, to level                     */
/*																			*/
/* Output:	  first Node, NULL if there is none to return					*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetFirstNode_hor_bw_up (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  NODE *theNode;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GNode_MG                = theMG;
  GNode_fromLevel = fromLevel;
  GNode_toLevel   = toLevel;

  theNode = LASTNODE(GRID_ON_LEVEL(GNode_MG,fromLevel));

  if (theNode==NULL)
    return (NULL);
  else if (USED(theNode))
    return (theNode);
  else
    return (NW_GetNextNode_hor_bw_up(theNode));
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetFirstNode_hor_bw_down									*/
/*																			*/
/* Purpose:   return first Node in horizontal, backward, downward directi	*/
/*																			*/
/* Input:	  the corresponding multigrid, from level, to level                     */
/*																			*/
/* Output:	  first Node, NULL if there is none to return					*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetFirstNode_hor_bw_down (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  NODE *theNode;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GNode_MG                = theMG;
  GNode_fromLevel = fromLevel;
  GNode_toLevel   = toLevel;

  theNode = LASTNODE(GRID_ON_LEVEL(GNode_MG,toLevel));

  if (theNode==NULL)
    return (NULL);
  else if (USED(theNode))
    return (theNode);
  else
    return (NW_GetNextNode_hor_bw_down(theNode));
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetFirstNode_leave_fw										*/
/*																			*/
/* Purpose:   return first leave-Node in forward direction					*/
/*																			*/
/* Input:	  the corresponding multigrid, from level, to level                     */
/*																			*/
/* Output:	  first Node, NULL if there is none to return					*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetFirstNode_leave_fw (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  NODE *theNode;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GNode_MG                = theMG;
  GNode_fromLevel = fromLevel;
  GNode_toLevel   = toLevel;

  theNode = TOPNODE(FIRSTVERTEX(GRID_ON_LEVEL(theMG,fromLevel)));
  while (LEVEL(theNode)>GNode_toLevel && theNode!=NULL) theNode=NFATHER(theNode);
  if (theNode==NULL) return (NULL);

  if (theNode==NULL)
    return (NULL);
  else if (USED(theNode))
    return (theNode);
  else
    return (NW_GetNextNode_leave_fw(theNode));
}

/****************************************************************************/
/*																			*/
/* Function:  NW_GetFirstNode_leave_bw										*/
/*																			*/
/* Purpose:   return first leave-Node in backward direction                             */
/*																			*/
/* Input:	  the corresponding multigrid, from level, to level                     */
/*																			*/
/* Output:	  first Node, NULL if there is none to return					*/
/*																			*/
/****************************************************************************/

static NODE *NW_GetFirstNode_leave_bw (MULTIGRID *theMG, INT fromLevel, INT toLevel)
{
  NODE *theNode;

  if (theMG==NULL)
    return (NULL);

  if (fromLevel<0)
    return (NULL);

  if (toLevel>TOPLEVEL(theMG))
    return (NULL);

  if (toLevel<fromLevel)
    return (NULL);

  GNode_MG                = theMG;
  GNode_fromLevel = fromLevel;
  GNode_toLevel   = toLevel;

  theNode = LASTNODE(GRID_ON_LEVEL(GNode_MG,fromLevel));
  while (LEVEL(theNode)>GNode_toLevel && theNode!=NULL) theNode=NFATHER(theNode);
  if (theNode==NULL) return (NULL);

  if (theNode==NULL)
    return (NULL);

  theNode = TOPNODE(LASTVERTEX(GRID_ON_LEVEL(theMG,toLevel)));

  if (USED(theNode))
    return (theNode);
  else
    return (NW_GetNextNode_leave_bw(theNode));
}

/****************************************************************************/
/*																			*/
/* Function:  GetFirstNode-Proc                                                                                         */
/*																			*/
/* Purpose:   get the GetFirstNodeProc										*/
/*																			*/
/* Input:	  NODE *theNode, char *theDrawingObject                                                 */
/*																			*/
/* Output:	  NULL: error													*/
/*																			*/
/****************************************************************************/

static NW_GetFirstNodeProcPtr NW_GetFirstNode_leave_fw_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetFirstNode_leave_fw);
}

static NW_GetFirstNodeProcPtr NW_GetFirstNode_leave_bw_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetFirstNode_leave_bw);
}

static NW_GetFirstNodeProcPtr NW_GetFirstNode_hor_fw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetFirstNode_hor_fw_up);
}

static NW_GetFirstNodeProcPtr NW_GetFirstNode_hor_fw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetFirstNode_hor_fw_down);
}

static NW_GetFirstNodeProcPtr NW_GetFirstNode_hor_bw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetFirstNode_hor_bw_up);
}

static NW_GetFirstNodeProcPtr NW_GetFirstNode_hor_bw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetFirstNode_hor_bw_down);
}

/****************************************************************************/
/*																			*/
/* Function:  GetNextNode-Proc												*/
/*																			*/
/* Purpose:   get the GetNextNodeProc										*/
/*																			*/
/* Input:	  NODE *theNode, char *theDrawingObject                                                 */
/*																			*/
/* Output:	  NULL: error													*/
/*																			*/
/****************************************************************************/

static NW_GetNextNodeProcPtr NW_GetNextNode_leave_fw_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetNextNode_leave_fw);
}

static NW_GetNextNodeProcPtr NW_GetNextNode_leave_bw_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetNextNode_leave_bw);
}

static NW_GetNextNodeProcPtr NW_GetNextNode_hor_fw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetNextNode_hor_fw_up);
}

static NW_GetNextNodeProcPtr NW_GetNextNode_hor_fw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetNextNode_hor_fw_down);
}

static NW_GetNextNodeProcPtr NW_GetNextNode_hor_bw_up_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetNextNode_hor_bw_up);
}

static NW_GetNextNodeProcPtr NW_GetNextNode_hor_bw_down_Proc (VIEWEDOBJ *theViewedObj)
{
  return (NW_GetNextNode_hor_bw_down);
}

/****************************************************************************/
/*
   EW_DoNothing0D - Do nothing

   SYNOPSIS:
   static INT EW_DoNothing0D (DRAWINGOBJ *q)

   PARAMETERS:
   .  q - the drawing object

   DESCRIPTION:
   This function does nothing (just for programming new ones).

   RETURN VALUE:
   INT

   0 when ok

   1 when erro occured
 */
/****************************************************************************/

static INT EW_DoNothing0D (DRAWINGOBJ *q)
{
  INT end;

  end = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc_RANGE(q);
      break;
    case DO_LINE :
      DO_inc_LINE(q,0);
      break;
    case DO_ARROW :
      DO_inc_ARROW(q,0);
      break;
    case DO_INVERSE_LINE :
      DO_inc_INVERSE_LINE(q,0);
      break;
    case DO_POLYLINE :
      DO_inc_POLYLINE(q,0);
      break;
    case DO_TEXT :
      DO_inc_TEXT(q,0);
      break;
    case DO_POLYMARK :
      DO_inc_POLYMARK(q,0);
      break;
    case DO_POLYGON :
      DO_inc_POLYGON(q,0);
      break;
    case DO_ERASE_SURRPOLYGON :
      DO_inc_ERASE_SURRPOLYGON(q,0);
      break;
    case DO_INVERSE_POLYGON :
      DO_inc_INVERSE_POLYGON(q,0);
      break;
    case DO_ERASE_POLYGON :
      DO_inc_ERASE_POLYGON(q,0);
      break;
    case DO_SURRPOLYGON :
      DO_inc_SURRPOLYGON(q,0);
      break;
    default :
      return (1);
    }

  }

  return (0);
}

/****************************************************************************/
/*
   Draw2D - Draw content of a 2D drawing object

   SYNOPSIS:
   static INT Draw2D (DRAWINGOBJ *q)

   PARAMETERS:
   .  q - the drawing object

   DESCRIPTION:
   This function draws content of a 2D drawing object.

   RETURN VALUE:
   INT

   0 when ok

   1 when error occured
 */
/****************************************************************************/

static INT Draw2D (DRAWINGOBJ *q)
{
  INT j, n, centered, end, mode;
  COORD help[2];
  COORD_POINT a, b, point[MAX_POINTS_OF_POLY];
  long color;

  end = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc_RANGE(q);
      break;
    case DO_LINE :
      DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,&a);
      UgMove(a);
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,&a);
      UgDraw(a);
      break;
    case DO_ARROW :
      DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,point);
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,point+3);
      point[1].x = ARR_ALPHA*point[3].x + (1.0-ARR_ALPHA)*point[0].x;
      point[1].y = ARR_ALPHA*point[3].y + (1.0-ARR_ALPHA)*point[0].y;
      a.x = point[3].x-point[1].x; a.y = point[3].y-point[1].y;
      point[2].x = point[1].x + ARR_COS*a.x - ARR_SIN*a.y;
      point[2].y = point[1].y + ARR_SIN*a.x + ARR_COS*a.y;
      point[4].x = point[1].x + ARR_COS*a.x + ARR_SIN*a.y;
      point[4].y = point[1].y - ARR_SIN*a.x + ARR_COS*a.y;
      point[5].x = point[1].x; point[5].y = point[1].y;
      UgPolyLine(point,6);
      break;
    case DO_INVERSE_LINE :
      DO_inc(q)
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,&a);
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,&b);
      UgInverseLine(a,b);
      break;
    case DO_POLYLINE :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q)
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,&a);
      UgMove(a);
      for (j=1; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,&a);
        UgDraw(a);
      }
      break;
    case DO_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgPolygon(point,j);
      break;
    case DO_ERASE_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgErasePolygon(point,j);
      break;
    case DO_INVERSE_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgInversePolygon(point,j);
      break;
    case DO_ERASE_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgErasePolygon(point,j);
      point[j].x=point[0].x; point[j].y=point[0].y;
      UgPolyLine(point,j+1);
      break;
    case DO_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      color = DO_2l(q); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgPolygon(point,j);
      UgSetColor(color);
      point[j].x=point[0].x; point[j].y=point[0].y;
      UgPolyLine(point,j+1);
      break;
    case DO_TEXT :
      DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      mode = TEXT_REGULAR; if (DO_2c(q)==TEXT_INVERSE) mode = TEXT_INVERSE;DO_inc(q)
      centered = 0; if (DO_2c(q)==TEXT_CENTERED) centered = 1;DO_inc(q)
      UgSetTextSize(DO_2s(q)); DO_inc(q);
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,&a);
      if (centered)
        UgCenteredText(a,(char*)q,mode);
      else
      {
        UgMove(a);
        UgText((char*)q,mode);
      }
      DO_inc_str(q);
      break;
    case DO_POLYMARK :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      UgSetMarker(DO_2s(q)); DO_inc(q);
      UgSetMarkerSize(DO_2s(q)); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgPolymark(point,n);
      break;
    default :
      return (1);
    }
  }

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_SelectNode2D												*/
/*																			*/
/* Purpose:   select a node                                                                                             */
/*																			*/
/* Input:	  char *q: the drawing object									*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT NW_SelectNode2D (DRAWINGOBJ *q)
{
  INT end, found;
  COORD help[2];
  COORD_POINT a, point[4];

  found = 0;
  end = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc_RANGE(q);
      break;
    case DO_LINE :
      DO_inc_LINE(q,2);
      break;
    case DO_ARROW :
      DO_inc_ARROW(q,2);
      break;
    case DO_INVERSE_LINE :
      DO_inc_INVERSE_LINE(q,2);
      break;
    case DO_POLYLINE :
      DO_inc_POLYLINE(q,2);
      break;
    case DO_TEXT :
      DO_inc_TEXT(q,2);
      break;
    case DO_POLYMARK :
      DO_inc(q)
      if (DO_2c(q)!=1) return (1);DO_inc_n(q,4);
      V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
      (*OBS_ProjectProc)(help,&a);
      if (ABS(FN2D_MousePos.x - a.x) < FN2D_ACC)
        if (ABS(FN2D_MousePos.y - a.y) < FN2D_ACC)
          found = 1;
      break;
    case DO_POLYGON :
      DO_inc_POLYGON(q,2);
      break;
    case DO_ERASE_SURRPOLYGON :
      DO_inc_ERASE_SURRPOLYGON(q,2);
      break;
    case DO_INVERSE_POLYGON :
      DO_inc_INVERSE_POLYGON(q,2);
      break;
    case DO_ERASE_POLYGON :
      DO_inc_ERASE_POLYGON(q,2);
      break;
    case DO_SURRPOLYGON :
      DO_inc_SURRPOLYGON(q,2);
      break;
    default :
      return (1);
    }

  }

  /* if found, put in selection list and invert */
  if (found)
  {
    /* put in/delete from selection list */
    if (SELECTIONMODE(WOP_MG)!=nodeSelection)
      ClearSelection(WOP_MG);
    if (AddNodeToSelection(WOP_MG,WOP_Node) == GM_ERROR)
      if (RemoveNodeFromSelection(WOP_MG,WOP_Node) == GM_ERROR)
        return (1);

    /* invert surrounding of node */
    point[0].x = point[3].x = a.x-FN2D_INVSIZE;
    point[0].y = point[1].y = a.y-FN2D_INVSIZE;
    point[2].x = point[1].x = a.x+FN2D_INVSIZE;
    point[2].y = point[3].y = a.y+FN2D_INVSIZE;
    UgInversePolygon(point,4);

    /* we have found a node */
    FN2D_found = 1;
  }

  return (0);
}

/****************************************************************************/
/*
   Draw3D - Draw content of a 2D drawing object

   SYNOPSIS:
   static INT Draw3D (DRAWINGOBJ *q);

   PARAMETERS:
   .  q - the drawing object

   DESCRIPTION:
   This function draws content of a 2D drawing object.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT Draw3D (DRAWINGOBJ *q)
{
  INT j, n, centered, end, mode;
  COORD help[3];
  COORD_POINT a, b, point[MAX_POINTS_OF_POLY];
  long color;

  end = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc_RANGE(q);
      break;
    case DO_LINE :
      DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,&a);
      UgMove(a);
      V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,&a);
      UgDraw(a);
      break;
    case DO_ARROW :
      DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,point);
      V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,point+3);
      point[1].x = ARR_ALPHA*point[3].x + (1.0-ARR_ALPHA)*point[0].x;
      point[1].y = ARR_ALPHA*point[3].y + (1.0-ARR_ALPHA)*point[0].y;
      a.x = point[3].x-point[1].x; a.y = point[3].y-point[1].y;
      point[2].x = point[1].x + ARR_COS*a.x - ARR_SIN*a.y;
      point[2].y = point[1].y + ARR_SIN*a.x + ARR_COS*a.y;
      point[4].x = point[1].x + ARR_COS*a.x + ARR_SIN*a.y;
      point[4].y = point[1].y - ARR_SIN*a.x + ARR_COS*a.y;
      point[5].x = point[1].x; point[5].y = point[1].y;
      UgPolyLine(point,6);
      break;
    case DO_INVERSE_LINE :
      DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,&a);
      V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,&b);
      UgInverseLine(a,b);
      break;
    case DO_POLYLINE :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,&a);
      UgMove(a);
      for (j=1; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,&a);
        UgDraw(a);
      }
      break;
    case DO_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgPolygon(point,j);
      break;
    case DO_ERASE_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgErasePolygon(point,j);
      break;
    case DO_INVERSE_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgInversePolygon(point,j);
      break;
    case DO_ERASE_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgErasePolygon(point,j);
      point[j].x=point[0].x; point[j].y=point[0].y;
      UgPolyLine(point,j+1);
      break;
    case DO_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      color = DO_2l(q); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgPolygon(point,j);
      UgSetColor(color);
      point[j].x=point[0].x; point[j].y=point[0].y;
      UgPolyLine(point,j+1);
      break;
    case DO_TEXT :
      DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      mode = TEXT_REGULAR; if (DO_2c(q)==TEXT_INVERSE) mode = TEXT_INVERSE;DO_inc(q)
      centered = 0; if (DO_2c(q) == TEXT_CENTERED) centered = 1;DO_inc(q)
      UgSetTextSize(DO_2s(q)); DO_inc(q);
      if ((DO_2c(q)) == TEXT_CENTERED)
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help);DO_inc_n(q,3);
      (*OBS_ProjectProc)(help,&a);
      if (centered)
        UgCenteredText(a,(char*)q,mode);
      else
      {
        UgMove(a);
        UgText((char*)q,mode);
      }
      DO_inc_str(q);
      break;
    case DO_POLYMARK :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      UgSetColor(DO_2l(q)); DO_inc(q);
      UgSetMarker(DO_2s(q)); DO_inc(q);
      UgSetMarkerSize(DO_2s(q)); DO_inc(q);
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      UgPolymark(point,n);
      break;
    default :
      return (1);
    }
  }

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  EW_FindElement2D												*/
/*																			*/
/* Purpose:   find element in 2D drawing object                                                         */
/*																			*/
/* Input:	  DRAWINGOBJ *q: the drawing object							*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT EW_SelectElement2D (DRAWINGOBJ *q)
{
  INT j, n, end, found;
  COORD help[3];
  COORD_POINT point[MAX_POINTS_OF_POLY];

  end = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc_RANGE(q);
      break;
    case DO_LINE :
      DO_inc_LINE(q,2);
      break;
    case DO_ARROW :
      DO_inc_ARROW(q,2);
      break;
    case DO_INVERSE_LINE :
      DO_inc_INVERSE_LINE(q,2);
      break;
    case DO_POLYLINE :
      DO_inc_POLYLINE(q,2);
      break;
    case DO_TEXT :
      DO_inc_TEXT(q,2);
      break;
    case DO_POLYMARK :
      DO_inc_POLYMARK(q,2);
      break;
    case DO_POLYGON :
    case DO_ERASE_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc_n(q,2)
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      found = PointInPolygon(point,j,FE2D_MousePos);
      break;
    case DO_INVERSE_POLYGON :
    case DO_ERASE_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      found = PointInPolygon(point,j,FE2D_MousePos);
      break;
    case DO_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc_n(q,3)
      for (j=0; j<n; j++)
      {
        V2_TRAFOM3_V2(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,2);
        (*OBS_ProjectProc)(help,point+j);
      }
      found = PointInPolygon(point,j,FE2D_MousePos);
      break;
    default :
      return (1);
    }

  }

  /* if found ... */
  if (found)
  {
    /* put in/delete from selection list */
    if (SELECTIONMODE(WOP_MG)!=elementSelection)
      ClearSelection(WOP_MG);
    if (AddElementToSelection(WOP_MG,WOP_Element) == GM_ERROR)
      if (RemoveElementFromSelection(WOP_MG,WOP_Element) == GM_ERROR)
        return (1);

    /* invert element */
    UgInversePolygon(point,n);

    /* something (i.e. an element) was found */
    FE2D_found = 1;
  }

  return (0);
}


/****************************************************************************/
/*
   FindRange2D - Find range for 2D

   SYNOPSIS:
   static INT FindRange2D (DRAWINGOBJ *q);

   PARAMETERS:
   q - the drawing object

   DESCRIPTION:
   This function finds range for 2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT FindRange2D (DRAWINGOBJ *q)
{
  INT end;

  end = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc(q)
      GEN_FR_min = MIN(GEN_FR_min,DO_2C(q)); DO_inc(q);
      GEN_FR_max = MAX(GEN_FR_max,DO_2C(q)); DO_inc(q);
      end = 1;
      break;
    case DO_LINE :
      DO_inc_LINE(q,2);
      break;
    case DO_ARROW :
      DO_inc_ARROW(q,2);
      break;
    case DO_INVERSE_LINE :
      DO_inc_INVERSE_LINE(q,2);
      break;
    case DO_POLYLINE :
      DO_inc_POLYLINE(q,2);
      break;
    case DO_TEXT :
      DO_inc_TEXT(q,2);
      break;
    case DO_POLYMARK :
      DO_inc_POLYMARK(q,2);
      break;
    case DO_POLYGON :
      DO_inc_POLYGON(q,2);
      break;
    case DO_ERASE_SURRPOLYGON :
      DO_inc_ERASE_SURRPOLYGON(q,2);
      break;
    case DO_INVERSE_POLYGON :
      DO_inc_INVERSE_POLYGON(q,2);
      break;
    case DO_ERASE_POLYGON :
      DO_inc_ERASE_POLYGON(q,2);
      break;
    case DO_SURRPOLYGON :
      DO_inc_SURRPOLYGON(q,2);
      break;
    default :
      return (1);
    }

  }

  return (0);
}

/****************************************************************************/
/*
   FindRange3D -  Find range for 3D

   SYNOPSIS:
   static INT FindRange3D (DRAWINGOBJ *q);

   PARAMETERS:
   .  q - the drawing object

   DESCRIPTION:
   This function finds range for 3D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT FindRange3D (DRAWINGOBJ *q)
{
  INT end;

  end = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc(q)
      GEN_FR_min = MIN(GEN_FR_min,DO_2C(q)); DO_inc(q);
      GEN_FR_max = MAX(GEN_FR_max,DO_2C(q)); DO_inc(q);
      break;
    case DO_LINE :
      DO_inc_LINE(q,3);
      break;
    case DO_ARROW :
      DO_inc_ARROW(q,3);
      break;
    case DO_INVERSE_LINE :
      DO_inc_INVERSE_LINE(q,3);
      break;
    case DO_POLYLINE :
      DO_inc_POLYLINE(q,3);
      break;
    case DO_TEXT :
      DO_inc_TEXT(q,3);
      break;
    case DO_POLYMARK :
      DO_inc_POLYMARK(q,3);
      break;
    case DO_POLYGON :
      DO_inc_POLYGON(q,3);
      break;
    case DO_ERASE_SURRPOLYGON :
      DO_inc_ERASE_SURRPOLYGON(q,3);
      break;
    case DO_INVERSE_POLYGON :
      DO_inc_INVERSE_POLYGON(q,3);
      break;
    case DO_ERASE_POLYGON :
      DO_inc_ERASE_POLYGON(q,3);
      break;
    case DO_SURRPOLYGON :
      DO_inc_SURRPOLYGON(q,3);
      break;
    default :
      return (1);
    }

  }

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  GEN_PostProcess_Scalar_FR                                                                         */
/*																			*/
/* Purpose:   postprocess for findrange of scalar plot						*/
/*																			*/
/* Input:	  PICTURE *thePicture, WORK *theWork							*/
/*																			*/
/* Return:	  INT 0: ok                                                                                                     */
/*			  INT 1: an error occurred										*/
/*																			*/
/****************************************************************************/

static INT GEN_PostProcess_Scalar_FR (PICTURE *thePicture, WORK *theWork)
{
  struct FindRange_Work *FR_Work;
  DOUBLE m,l;
  INT i;

  FR_Work = W_FINDRANGE_WORK(theWork);

  if (GEN_FR_min>GEN_FR_max)
  {
    UserWrite("findrange failed\n");
    return (0);
  }

  /* postprocess findrange */
  if (FR_Work->symmetric==YES)
  {
    GEN_FR_max = MAX(ABS(GEN_FR_min),ABS(GEN_FR_max));
    GEN_FR_min = -GEN_FR_max;
  }
  if (FR_Work->zoom!=1.0)
  {
    m = 0.5*(GEN_FR_max + GEN_FR_min);
    l = 0.5*(GEN_FR_max - GEN_FR_min);
    GEN_FR_min = m - FR_Work->zoom*l;
    GEN_FR_max = m + FR_Work->zoom*l;
  }
  FR_Work->min = GEN_FR_min;
  FR_Work->max = GEN_FR_max;

  /* store if */
  if (GEN_FR_put == YES)
  {
    PIC_PO(thePicture)->theEspo.min = GEN_FR_min;
    PIC_PO(thePicture)->theEspo.max = GEN_FR_max;
    if (PIC_PO(thePicture)->theEspo.mode==PO_CONTOURS_EQ)
    {
      for (i=0; i<PIC_PO(thePicture)->theEspo.numOfContours; i++)
        PIC_PO(thePicture)->theEspo.contValues[i] = GEN_FR_min + (DOUBLE)i * (GEN_FR_max - GEN_FR_min) / (DOUBLE)(PIC_PO(thePicture)->theEspo.numOfContours-1);
    }
  }

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  EW_PostProcess_EVector_FR                                                                         */
/*																			*/
/* Purpose:   postprocess for findrange of vector plot						*/
/*																			*/
/* Input:	  PICTURE *thePicture, WORK *theWork							*/
/*																			*/
/* Return:	  INT 0: ok                                                                                                     */
/*			  INT 1: an error occurred										*/
/*																			*/
/****************************************************************************/

static INT GEN_PostProcess_Vector_FR (PICTURE *thePicture, WORK *theWork)
{
  struct FindRange_Work *FR_Work;

  FR_Work = W_FINDRANGE_WORK(theWork);

  /* postprocess findrange */
  FR_Work->min = 0.0;
  if (FR_Work->zoom>0.0)
    FR_Work->max = GEN_FR_max*FR_Work->zoom;
  else
    FR_Work->max = GEN_FR_max;

  /* store if */
  if (GEN_FR_put == YES)
    PIC_PO(thePicture)->theEvpo.max = FR_Work->max;

  return (0);
}

/**********************************************************************************************************/
/************************************ Part only for 2D Version ********************************************/
/**********************************************************************************************************/

#ifdef __TWODIM__

/****************************************************************************/
/*
   OrderElements - Order elements w.r.t. theViewedObject

   SYNOPSIS:
   static INT OrderElements_2D (MULTIGRID *theMG, VIEWEDOBJ *theViewedObj)

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  heViewedObj -

   DESCRIPTION:
   This function order elements w.r.t. theViewedObject.

   RETURN VALUE:
   INT

   0 when ok

   1 when error occured
 */
/****************************************************************************/

static INT OrderElements_2D (MULTIGRID *theMG, VIEWEDOBJ *theViewedObj)
{
  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_PlotBndOfElem2D - Initialize input variables of EW_BndOfElemEval2D

   SYNOPSIS:
   static INT EW_PreProcess_PlotBndOfElem2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes input variables of EW_BndOfElemEval2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_PlotBndOfElem2D (PICTURE *thePicture, WORK *theWork)
{
  struct GridPlotObj2D *theGpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theGpo = &(PIC_PO(thePicture)->theGpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  /* see if boundary has to be plotted */
  if (theGpo->PlotBoundary == NO)
    return (1);

  EE2D_NoColor[COLOR_BND]                 = 0;
  EE2D_Color[COLOR_BND]                   = theOD->blue;
  EE2D_Elem2Plot[PLOT_ALL]                = 1;
  EE2D_Elem2Plot[PLOT_COPY]               = 1;
  EE2D_Elem2Plot[PLOT_IRR]                = 1;
  EE2D_Elem2Plot[PLOT_REG]                = 1;

  /* mark surface elements */
  if (MarkElements_MGS_Bnd(theMG,0,CURRENTLEVEL(theMG))) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_PlotBlackBnd2D	- Initialize input variables of EW_BndOfElemEval2D

   SYNOPSIS:
   static INT EW_PreProcess_PlotBlackBnd2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes input variables of EW_BndOfElemEval2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_PlotBlackBnd2D (PICTURE *thePicture, WORK *theWork)
{
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  EE2D_NoColor[COLOR_BND]                 = 0;
  EE2D_Color[COLOR_BND]                   = theOD->black;
  EE2D_Elem2Plot[PLOT_ALL]                = 1;
  EE2D_Elem2Plot[PLOT_COPY]               = 1;
  EE2D_Elem2Plot[PLOT_IRR]                = 1;
  EE2D_Elem2Plot[PLOT_REG]                = 1;

  /* mark surface elements */
  if (MarkElements_MGS_Bnd(theMG,0,CURRENTLEVEL(theMG))) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_PlotElements2D - Initialize input variables of EW_ElementEval2D for GridPlot2D

   SYNOPSIS:
   static INT EW_PreProcess_PlotElements2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes input variables of EW_ElementEval2D for GridPlot2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_PlotElements2D (PICTURE *thePicture, WORK *theWork)
{
  struct GridPlotObj2D *theGpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theGpo = &(PIC_PO(thePicture)->theGpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  if (theGpo->WhichElem == PO_NO && theGpo->PlotElemID == NO)
    return (1);

  EE2D_NoColor[COLOR_EDGE]                        = 0;
  EE2D_NoColor[COLOR_LOWER_LEVEL]         = 1;
  EE2D_NoColor[COLOR_BND]                         = 1;
  if (theGpo->ElemColored == YES)
  {
    EE2D_NoColor[COLOR_COPY]                = 0;
    EE2D_NoColor[COLOR_IRR]                 = 0;
    EE2D_NoColor[COLOR_REG]                 = 0;
  }
  else
  {
    EE2D_NoColor[COLOR_COPY]                = 1;
    EE2D_NoColor[COLOR_IRR]                 = 1;
    EE2D_NoColor[COLOR_REG]                 = 1;
  }
  if (theGpo->PlotBoundary == YES)
    EE2D_NoColor[COLOR_BND]                 = 0;


  EE2D_Color[COLOR_COPY]                  = theOD->yellow;
  EE2D_Color[COLOR_IRR]                   = theOD->green;
  EE2D_Color[COLOR_REG]                   = theOD->red;
  EE2D_Color[COLOR_LOWER_LEVEL]   = theOD->white;
  EE2D_Color[COLOR_EDGE]                  = theOD->black;
  EE2D_Color[COLOR_BND]                   = theOD->blue;
  EE2D_Color[COLOR_ELEMID]                = theOD->orange;

  EE2D_Elem2Plot[PLOT_ALL]                = 0;
  EE2D_Elem2Plot[PLOT_COPY]               = 0;
  EE2D_Elem2Plot[PLOT_IRR]                = 0;
  EE2D_Elem2Plot[PLOT_REG]                = 0;


  switch (theGpo->WhichElem)
  {
  case PO_NO :
    break;
  case PO_ALL :
    EE2D_Elem2Plot[PLOT_ALL] = 1;
  case PO_COPY :
    EE2D_Elem2Plot[PLOT_COPY] = 1;
  case PO_IRR :
    EE2D_Elem2Plot[PLOT_IRR] = 1;
  case PO_REG :
    EE2D_Elem2Plot[PLOT_REG] = 1;
  }

  EE2D_ElemID                                     = 0;
  if (theGpo->PlotElemID == YES)
    EE2D_ElemID                             = 1;

  /* mark surface elements */
  EE2D_MaxLevel = CURRENTLEVEL(theMG);
  if (MarkElements2D(theMG,0,EE2D_MaxLevel)) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_SelectElement2D - Initialize input variables finding element 2D

   SYNOPSIS:
   static INT EW_PreProcess_SelectElement2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes input variables finding element 2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_SelectElement2D (PICTURE *thePicture, WORK *theWork)
{
  struct GridPlotObj2D *theGpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theGpo = &(PIC_PO(thePicture)->theGpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  EE2D_NoColor[COLOR_EDGE]                        = 1;
  EE2D_NoColor[COLOR_LOWER_LEVEL]         = 1;
  EE2D_NoColor[COLOR_BND]                         = 1;
  EE2D_NoColor[COLOR_COPY]                        = 1;
  EE2D_NoColor[COLOR_IRR]                         = 1;
  EE2D_NoColor[COLOR_REG]                         = 1;
  EE2D_NoColor[COLOR_BND]                         = 1;

  EE2D_Elem2Plot[PLOT_ALL]                = 0;
  EE2D_Elem2Plot[PLOT_COPY]               = 0;
  EE2D_Elem2Plot[PLOT_IRR]                = 0;
  EE2D_Elem2Plot[PLOT_REG]                = 0;
  switch (theGpo->WhichElem)
  {
  case PO_NO :
    break;
  case PO_ALL :
    EE2D_Elem2Plot[PLOT_ALL] = 1;
  case PO_COPY :
    EE2D_Elem2Plot[PLOT_COPY] = 1;
  case PO_IRR :
    EE2D_Elem2Plot[PLOT_IRR] = 1;
  case PO_REG :
    EE2D_Elem2Plot[PLOT_REG] = 1;
  }

  FE2D_found                                                      = 0;

  /* store mouse position */
  FE2D_MousePos.x = W_SELECTELEMENT_WORK(theWork)->PixelX;
  FE2D_MousePos.y = W_SELECTELEMENT_WORK(theWork)->PixelY;

  if (MarkElements2D(theMG,0,EE2D_MaxLevel)) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_BndOfElemEval2D -  Evaluate bnd of 2D element

   SYNOPSIS:
   static INT EW_BndOfElemEval2D (ELEMENT *theElement, DRAWINGOBJ *theDO);

   PARAMETERS:
   .  theElement -
   .  theDO -

   DESCRIPTION:
   This function evaluates bnd of 2D element (triangle/quadrilateral).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_BndOfElemEval2D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, n;
  COORD *x[MAX_CORNERS_OF_ELEM];

  /* get coordinates of corners of the element */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* store bnd sides on drawing obj */
  n = CORNERS_OF_ELEM(theElement);
  for (i=0; i<n; i++)
  {
    if (SIDE(theElement,i) == NULL) continue;
    DO_2c(theDO) = DO_LINE; DO_inc(theDO)
    DO_2l(theDO) = EE2D_Color[COLOR_BND]; DO_inc(theDO);
    V2_COPY(x[i],DO_2Cp(theDO)); DO_inc_n(theDO,2);
    V2_COPY(x[(i+1)%n],DO_2Cp(theDO)); DO_inc_n(theDO,2);
  }

  DO_2c(theDO) = DO_NO_INST;

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_PreProcess_PlotNodes2D                                                                         */
/*																			*/
/* Purpose:   initialize input variables of NW_PlotNodes2D					*/
/*																			*/
/* Input:	  PICTURE *thePicture, WORK *theWork							*/
/*																			*/
/* Return:	  INT 0: ok                                                                                                     */
/*			  INT 1: an error occurred										*/
/*																			*/
/****************************************************************************/

static INT NW_PreProcess_PlotNodes2D (PICTURE *thePicture, WORK *theWork)
{
  struct GridPlotObj2D *theGpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;
  INT mode;

  theGpo = &(PIC_PO(thePicture)->theGpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  NE_IDColor                                      = theOD->black;
  NE_BndMarkerColor                       = theOD->red;
  NE_InnerMarkerColor             = theOD->red;
  NE_InnerMarker                          = FILLED_CIRCLE_MARKER;
  NE_BndMarker                            = FILLED_SQUARE_MARKER;
  NE_InnerMarkerSize                      = 4;
  NE_BndMarkerSize                        = 4;

  NE_EvalNodeID                           = 0;
  NE_EvalInnerNode                        = 0;
  NE_EvalBndNode                          = 0;
  if (theGpo->PlotNodeID == YES)
    NE_EvalNodeID                   = 1;
  if (theGpo->PlotNodes == YES)
  {
    NE_EvalInnerNode                = 1;
    NE_EvalBndNode                  = 1;
  }

  /* mark nodes */
  switch (theGpo->WhichElem)
  {
  case PO_ALL :
    mode = MARKMODE_ALL;
    break;
  case PO_COPY :
    mode = MARKMODE_COPY;
    break;
  case PO_IRR :
    mode = MARKMODE_IRREG;
    break;
  case PO_REG :
    mode = MARKMODE_REG;
    break;
  default :
    return (1);
  }
  if (MarkNodes_OfMarkedElem(theMG,0,CURRENTLEVEL(theMG),mode)) return (1);

  return (0);
}

/****************************************************************************/
/*
   NW_PreProcess_SelectNode2D - Initialize input variables of EW_BndOfElemEval2D

   SYNOPSIS:
   static INT NW_PreProcess_SelectNode2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes input variables of EW_BndOfElemEval2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT NW_PreProcess_SelectNode2D (PICTURE *thePicture, WORK *theWork)
{
  struct GridPlotObj2D *theGpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;
  INT mode;

  theGpo = &(PIC_PO(thePicture)->theGpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  NE_InnerMarkerColor             = theOD->red;
  NE_BndMarkerColor                       = theOD->red;
  NE_InnerMarker                          = FILLED_CIRCLE_MARKER;
  NE_BndMarker                            = FILLED_SQUARE_MARKER;
  NE_InnerMarkerSize                      = 4;
  NE_BndMarkerSize                        = 4;

  NE_EvalNodeID                           = 0;
  NE_EvalInnerNode                        = 1;
  NE_EvalBndNode                          = 1;

  FN2D_found                                      = 0;

  /* store mouse position */
  FN2D_MousePos.x = W_SELECTNODE_WORK(theWork)->PixelX;
  FN2D_MousePos.y = W_SELECTNODE_WORK(theWork)->PixelY;

  /* mark nodes */
  switch (theGpo->WhichElem)
  {
  case PO_ALL :
    mode = MARKMODE_ALL;
    break;
  case PO_COPY :
    mode = MARKMODE_COPY;
    break;
  case PO_IRR :
    mode = MARKMODE_IRREG;
    break;
  case PO_REG :
    mode = MARKMODE_REG;
    break;
  default :
    return (1);
  }
  if (MarkNodes_OfMarkedElem(theMG,0,CURRENTLEVEL(theMG),mode)) return (1);

  return (0);
}

/****************************************************************************/
/*
   InvertElementSelection2D - invert element selection

   SYNOPSIS:
   static INT InvertElementSelection2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function inverts element selection.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InvertElementSelection2D (PICTURE *thePicture, WORK *theWork)
{
  ELEMENT *theElement;
  COORD_VECTOR help;
  COORD_POINT points[4];
  INT i, j;

  /* evaluate and execute */
  if (SELECTIONMODE(WOP_MG)==elementSelection)
    for (i=0; i<SELECTIONSIZE(WOP_MG); i++)
    {
      theElement = (ELEMENT *)SELECTIONOBJECT(WOP_MG,i);
      if (!USED(theElement)) continue;
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      {
        V2_TRAFOM3_V2(CVECT(MYVERTEX(CORNER(theElement,j))),ObsTrafo,help);
        (*OBS_ProjectProc)(help,points+j);
      }
      UgInversePolygon(points,j);
    }

  return (0);
}

/****************************************************************************/
/*
   /* Function: InvertNodeSelection2D                                                                                   */
/*																			*/
/* Purpose:   invert node selection                                                                             */
/*																			*/
/* Input:	  PICTURE *thePicture, WORK *theWork							*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT InvertNodeSelection2D (PICTURE *thePicture, WORK *theWork)
{
  NODE *theNode;
  COORD_VECTOR help;
  COORD_POINT a, point[4];
  INT i;

  /* evaluate and execute */
  if (SELECTIONMODE(WOP_MG)==nodeSelection)
    for (i=0; i<SELECTIONSIZE(WOP_MG); i++)
    {
      theNode = (NODE *)SELECTIONOBJECT(WOP_MG,i);
      if (!USED(theNode)) continue;
      V2_TRAFOM3_V2(CVECT(MYVERTEX(theNode)),ObsTrafo,help);
      (*OBS_ProjectProc)(help,&a);

      /* invert surrounding of node */
      point[0].x = point[3].x = a.x-FN2D_INVSIZE;
      point[0].y = point[1].y = a.y-FN2D_INVSIZE;
      point[2].x = point[1].x = a.x+FN2D_INVSIZE;
      point[2].y = point[3].y = a.y+FN2D_INVSIZE;
      UgInversePolygon(point,4);
    }

  return (0);
}

/****************************************************************************/
/*
   EW_ElementEval2D -

   SYNOPSIS:
   static INT EW_ElementEval2D (ELEMENT *theElement, DRAWINGOBJ *theDO);

   PARAMETERS:
   .  theElement -
   .  theDO -

   DESCRIPTION:
   This function evaluates geometry of 2D element (triangle/quadrilateral).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_ElementEval2D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, j;
  COORD *x[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoint;

  /* get coordinates of corners of the element */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* store viewable sides on drawing obj */
  if (LEVEL(theElement)<EE2D_MaxLevel)
  {
    if (EE2D_NoColor[COLOR_LOWER_LEVEL])
    {
      DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
      DO_2c(theDO) = CORNERS_OF_ELEM(theElement); DO_inc(theDO)
    }
    else
    {
      DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
      DO_2c(theDO) = CORNERS_OF_ELEM(theElement); DO_inc(theDO)
      DO_2l(theDO) = EE2D_Color[COLOR_LOWER_LEVEL]; DO_inc(theDO);
    }
    DO_2l(theDO) = EE2D_Color[COLOR_EDGE]; DO_inc(theDO);
    for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
    {
      V2_COPY(x[j],DO_2Cp(theDO));
      DO_inc_n(theDO,2);
    }
  }
  else
  {
    if (EE2D_NoColor[ECLASS(theElement)])
    {
      DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
      DO_2c(theDO) = CORNERS_OF_ELEM(theElement); DO_inc(theDO)
    }
    else
    {
      DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
      DO_2c(theDO) = CORNERS_OF_ELEM(theElement); DO_inc(theDO)
      DO_2l(theDO) = EE2D_Color[ECLASS(theElement)]; DO_inc(theDO);
    }
    DO_2l(theDO) = EE2D_Color[COLOR_EDGE]; DO_inc(theDO);
    for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
    {
      V2_COPY(x[j],DO_2Cp(theDO));
      DO_inc_n(theDO,2);
    }
  }

  /* plot element ID */
  if (EE2D_ElemID)
  {
    V2_CLEAR(MidPoint)
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      V2_ADD(MidPoint,x[i],MidPoint)
      V2_SCALE(1.0/(COORD)i,MidPoint)

      DO_2c(theDO) = DO_TEXT;DO_inc(theDO)
    DO_2l(theDO) = EE2D_Color[COLOR_ELEMID]; DO_inc(theDO);
    DO_2c(theDO) = TEXT_REGULAR; DO_inc(theDO)
    DO_2c(theDO) = TEXT_CENTERED; DO_inc(theDO)
    DO_2s(theDO) = EE2D_TEXTSIZE; DO_inc(theDO);
    V2_COPY(MidPoint,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    sprintf(DO_2cp(theDO),"%d",(int)ID(theElement)); DO_inc_str(theDO);
  }

  DO_2c(theDO) = DO_NO_INST;

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  NW_NodesEval2D												*/
/*																			*/
/* Purpose:   evaluate node                                                                                             */
/*																			*/
/* Input:	  NODE *theNode, char *theDrawingObject                                                 */
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT NW_NodesEval2D (NODE *theNode, DRAWINGOBJ *theDO)
{
  /* plot marks of boundary nodes */
  if (NE_EvalInnerNode)
  {
    DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO)
    DO_2c(theDO) = 1; DO_inc(theDO)
    DO_2l(theDO) = NE_BndMarkerColor; DO_inc(theDO);
    DO_2s(theDO) = NE_BndMarker; DO_inc(theDO);
    DO_2s(theDO) = NE_BndMarkerSize; DO_inc(theDO);
    V2_COPY(CVECT(MYVERTEX(theNode)),DO_2Cp(theDO)); DO_inc_n(theDO,2);
  }

  /* plot marks of inner nodes */
  if (NE_EvalInnerNode)
  {
    DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO)
    DO_2c(theDO) = 1; DO_inc(theDO)
    DO_2l(theDO) = NE_InnerMarkerColor; DO_inc(theDO);
    DO_2s(theDO) = NE_InnerMarker; DO_inc(theDO);
    DO_2s(theDO) = NE_InnerMarkerSize; DO_inc(theDO);
    V2_COPY(CVECT(MYVERTEX(theNode)),DO_2Cp(theDO)); DO_inc_n(theDO,2);
  }

  /* plot node ID */
  if (NE_EvalNodeID)
  {
    DO_2c(theDO) = DO_TEXT; DO_inc(theDO)
    DO_2l(theDO) = NE_IDColor; DO_inc(theDO)
    DO_2c(theDO) = TEXT_REGULAR; DO_inc(theDO)
    DO_2c(theDO) = TEXT_NOT_CENTERED; DO_inc(theDO)
    DO_2s(theDO) = EE2D_TEXTSIZE; DO_inc(theDO);
    V2_COPY(CVECT(MYVERTEX(theNode)),DO_2Cp(theDO)); DO_inc_n(theDO,2);
    sprintf(DO_2cp(theDO),"%d",(int)ID(theNode)); DO_inc_str(theDO);
  }

  DO_2c(theDO) = DO_NO_INST;

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EScalar3D - Initialize for C(olor)C(ontour) plot

   SYNOPSIS:
   static INT EW_PreProcess_EScalar2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for C(olor)C(ontour) plot.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_EScalar2D (PICTURE *thePicture, WORK *theWork)
{
  struct ElemScalarPlotObj2D *theEspo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;
  INT i;

  theEspo = &(PIC_PO(thePicture)->theEspo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  /* set value->color fct, eval fct */
  if (theEspo->max - theEspo->min < SMALL_D)
    if (W_ID(theWork) != FINDRANGE_WORK)
    {
      UserWrite("maxValue has to be larger than minValue\n");
      return (1);
    }
  EScalar2D_EvalFct        = theEspo->EvalFct->PlotProc;
  if ((theEspo->max - theEspo->min)==0)
    EScalar2D_V2C_factor = 0;
  else
    EScalar2D_V2C_factor = (theOD->spectrumEnd - theOD->spectrumStart)/(theEspo->max - theEspo->min);
  EScalar2D_V2C_offset = theOD->spectrumStart - EScalar2D_V2C_factor*theEspo->min;
  EScalar2D_mode           = theEspo->mode;
  if (EScalar2D_mode == PO_CONTOURS_EQ)
  {
    EScalar2D_numOfContours = theEspo->numOfContours;
    for (i=0; i<theEspo->numOfContours; i++)
      EScalar2D_ContValues[i] = theEspo->contValues[i];
  }
  EScalar2D_depth                 = theEspo->depth;

  /* mark suface elements on boundary */
  if (MarkElements_MGS(theMG,0,CURRENTLEVEL(theMG))) return (1);

  /* prepare evaluation routine */
  if (theEspo->EvalFct->PreprocessProc!=NULL)
    if ((*theEspo->EvalFct->PreprocessProc)(theMG))
      return (1);;

  return (0);
}

/****************************************************************************/
/*
   PlotColorTriangle2D -  Plot on triangle color(2D coord) with depth

   SYNOPSIS:
   static INT PlotColorTriangle2D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *TP0, COORD *TP1, COORD *TP2,
   INT depth, DRAWINGOBJ **PtrDO);

   PARAMETERS:
   .  theElement -
   .  CornersOfElem - its corners
   .  TP0 -
   .  TP1 -
   .  TP2 -
   .  depth -
   .  PtrDO -

   DESCRIPTION:
   This function plots on triangle color(2D coord) with depth.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PlotColorTriangle2D (ELEMENT *theElement, const COORD **CornersOfElem, const COORD *TP0, const COORD *TP1, const COORD *TP2, INT depth, DRAWINGOBJ **PtrDO)
{
  COORD_VECTOR EvalPoint, LocalCoord, MP0, MP1, MP2;
  INT i;
  long Color;
  DOUBLE value;

  if (depth<=0)
  {
    /* get values */
    for (i=0; i<DIM; i++)
      EvalPoint[i] = (TP0[i]+TP1[i]+TP2[i])/3.0;
    if (GlobalToLocal2d(3,CornersOfElem,EvalPoint,LocalCoord)) return (1);
    value = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    Color = (long)(EScalar2D_V2C_factor*value+EScalar2D_V2C_offset);
    Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
    Color = MAX(Color,WOP_OutputDevice->spectrumStart);

    /* draw */
    DO_2c(*PtrDO) = DO_POLYGON; DO_inc(*PtrDO)
    DO_2c(*PtrDO) = 3; DO_inc(*PtrDO)
    DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
    V2_COPY(TP0,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
    V2_COPY(TP1,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
    V2_COPY(TP2,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);

    /* store range */
    EScalar2D_minValue = MIN(EScalar2D_minValue,value);
    EScalar2D_maxValue = MAX(EScalar2D_maxValue,value);
  }
  else
  {
    /* find corners of subdivided triangles */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(TP0[i]+TP1[i]);
      MP1[i] = 0.5*(TP1[i]+TP2[i]);
      MP2[i] = 0.5*(TP2[i]+TP0[i]);
    }
    if (PlotColorTriangle2D(theElement,CornersOfElem,TP0,MP0,MP2,depth-1,PtrDO)) return (1);
    if (PlotColorTriangle2D(theElement,CornersOfElem,MP0,TP1,MP1,depth-1,PtrDO)) return (1);
    if (PlotColorTriangle2D(theElement,CornersOfElem,TP2,MP2,MP1,depth-1,PtrDO)) return (1);
    if (PlotColorTriangle2D(theElement,CornersOfElem,MP0,MP1,MP2,depth-1,PtrDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   PlotColorQuadrilateral2D - Plot on quadrilateral color(2D coord) with depth

   SYNOPSIS:
   static INT PlotColorQuadrilateral2D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *QP0, COORD *QP1, COORD *QP2,
   COORD *QP3, INT depth, DRAWINGOBJ **PtrDO);

   PARAMETERS:
   .  theElement -
   .  CornersOfElem -
   .  QP0 -
   .  QP1 -
   .  QP2 -
   .  QP3 -
   .  depth -
   .  PtrDO - the drawing object to draw on

   DESCRIPTION:
   This function plots on quadrilateral color(2D coord) with depth.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PlotColorQuadrilateral2D (ELEMENT *theElement, const COORD **CornersOfElem, const COORD *QP0, const COORD *QP1, const COORD *QP2, const COORD *QP3, INT depth, DRAWINGOBJ **PtrDO)
{
  COORD_VECTOR EVP, LocalCoord, MP0, MP1, MP2, MP3;
  INT i;
  long Color;
  DOUBLE value;

  for (i=0; i<DIM; i++)
    EVP[i] = (QP0[i]+QP1[i]+QP2[i]+QP3[i])*0.25;
  if (depth<=0)
  {
    if (ID(theElement)==98)
      i=0;
    /* get values */
    if (GlobalToLocal2d(4,CornersOfElem,EVP,LocalCoord)) return (1);
    value = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    Color = (long)(EScalar2D_V2C_factor*value+EScalar2D_V2C_offset);
    Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
    Color = MAX(Color,WOP_OutputDevice->spectrumStart);

    /* draw */
    DO_2c(*PtrDO) = DO_POLYGON; DO_inc(*PtrDO)
    DO_2c(*PtrDO) = 4; DO_inc(*PtrDO)
    DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
    V2_COPY(QP0,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
    V2_COPY(QP1,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
    V2_COPY(QP2,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
    V2_COPY(QP3,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);

    /* store range */
    EScalar2D_minValue = MIN(EScalar2D_minValue,value);
    EScalar2D_maxValue = MAX(EScalar2D_maxValue,value);
  }
  else
  {
    /* find corners of subdivided quadrilaterals */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(QP0[i]+QP1[i]);
      MP1[i] = 0.5*(QP1[i]+QP2[i]);
      MP2[i] = 0.5*(QP2[i]+QP3[i]);
      MP3[i] = 0.5*(QP3[i]+QP0[i]);
    }
    if (PlotColorQuadrilateral2D(theElement,CornersOfElem,QP0,MP0,EVP,MP3,depth-1,PtrDO)) return (1);
    if (PlotColorQuadrilateral2D(theElement,CornersOfElem,MP0,QP1,MP1,EVP,depth-1,PtrDO)) return (1);
    if (PlotColorQuadrilateral2D(theElement,CornersOfElem,EVP,MP1,QP2,MP2,depth-1,PtrDO)) return (1);
    if (PlotColorQuadrilateral2D(theElement,CornersOfElem,EVP,MP2,QP3,MP3,depth-1,PtrDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   PointOnLine2D - Cals point between two points with contourValue

   SYNOPSIS:
   static INT PointOnLine2D (DOUBLE contourValue, DOUBLE value0,
   DOUBLE value1, COORD_VECTOR vec0, COORD_VECTOR vec1,
   COORD_VECTOR p);

   PARAMETERS:
   .  contourValue -
   .  value0 -
   .  value1 -
   .  vec0 -
   .  vec1 -
   .  VECTOR p -

   DESCRIPTION:
   This function cals point between two points with contourValue.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PointOnLine2D (DOUBLE contourValue, DOUBLE value0, DOUBLE value1, const COORD_VECTOR vec0, const COORD_VECTOR vec1, COORD_VECTOR p)
{
  DOUBLE alpha;

  if (value0==value1)
  {
    if (value0==contourValue)
    {
      V2_LINCOMB(0.5,vec0,0.5,vec1,p);
      return (1);
    }
    return (0);
  }
  else
  {
    alpha = (contourValue-value0)/(value1-value0);
    if (alpha<0.0 || alpha>1.0) return (0);
    V2_LINCOMB(1.0-alpha,vec0,alpha,vec1,p);
    return (1);
  }
}

/****************************************************************************/
/*
   PlotContourTriangle2D - plot on triangle contourlines (2D coord) with depth

   SYNOPSIS:
   static INT PlotContourTriangle2D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *TP0, COORD *TP1, COORD *TP2,
   INT depth, DRAWINGOBJ **PtrDO);

   PARAMETER:
   .  theElement -
   .  CornersOfElem -
   .  TP0 -
   .  TP1 -
   .  TP2 -
   .  depth -
   .  PtrDO - the drawing object to draw on

   DESCRIPTION:
   This function plots on triangle contourlines (2D coord) with depth.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PlotContourTriangle2D (ELEMENT *theElement, const COORD **CornersOfElem, const COORD *TP0, const COORD *TP1, const COORD *TP2, INT depth, DRAWINGOBJ **PtrDO)
{
  COORD_VECTOR LocalCoord, MP0, MP1, MP2, PointMid, Point[3];
  INT i, j, n, min, max;
  long Color;
  DOUBLE v0, v1, v2, vmin, vmax;

  if (depth<=0)
  {
    /* get values at the corners */
    if (GlobalToLocal2d(3,CornersOfElem,TP0,LocalCoord)) return (1);
    v0      = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal2d(3,CornersOfElem,TP1,LocalCoord)) return (1);
    v1      = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal2d(3,CornersOfElem,TP2,LocalCoord)) return (1);
    v2      = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    vmin = MIN(v0,v1); vmin = MIN(vmin,v2);
    vmax = MAX(v0,v1); vmax = MAX(vmax,v2);

    /* store range */
    EScalar2D_minValue = MIN(EScalar2D_minValue,vmin);
    EScalar2D_maxValue = MAX(EScalar2D_maxValue,vmax);

    /* find contours to be plotted */
    for (min=0; min<EScalar2D_numOfContours; min++)
      if (EScalar2D_ContValues[min]>=vmin)
        break;
    for (max=EScalar2D_numOfContours-1; max>=0; max--)
      if (EScalar2D_ContValues[max]<=vmax)
        break;

    /* draw contours */
    for (i=min; i<=max; i++)
    {
      /* set color */
      Color = (long)(EScalar2D_V2C_factor*EScalar2D_ContValues[i]+EScalar2D_V2C_offset);
      Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
      Color = MAX(Color,WOP_OutputDevice->spectrumStart);

      /* calculate points on each side of triangle having the right value */
      n=0;
      if (PointOnLine2D(EScalar2D_ContValues[i],v0,v1,TP0,TP1,Point[n])) n++;
      if (PointOnLine2D(EScalar2D_ContValues[i],v1,v2,TP1,TP2,Point[n])) n++;
      if (PointOnLine2D(EScalar2D_ContValues[i],v2,v0,TP2,TP0,Point[n])) n++;

      /* draw */
      switch (n)
      {
      case 2 :
        DO_2c(*PtrDO) = DO_LINE; DO_inc(*PtrDO)
        DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
        V2_COPY(Point[0],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[1],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        break;
      case 3 :
        DO_2c(*PtrDO) = DO_POLYLINE; DO_inc(*PtrDO)
        DO_2c(*PtrDO) = 5; DO_inc(*PtrDO)
        DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
        for (j=0; j<DIM; j++)
          PointMid[j] = (Point[0][j]+Point[1][j]+Point[2][j])/3.0;
        V2_COPY(Point[0],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(PointMid,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[1],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(PointMid,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[2],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        break;
      }
    }
  }
  else
  {
    /* find corners of subdivided triangles */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(TP0[i]+TP1[i]);
      MP1[i] = 0.5*(TP1[i]+TP2[i]);
      MP2[i] = 0.5*(TP2[i]+TP0[i]);
    }
    if (PlotContourTriangle2D(theElement,CornersOfElem,TP0,MP0,MP2,depth-1,PtrDO)) return (1);
    if (PlotContourTriangle2D(theElement,CornersOfElem,MP0,TP1,MP1,depth-1,PtrDO)) return (1);
    if (PlotContourTriangle2D(theElement,CornersOfElem,TP2,MP2,MP1,depth-1,PtrDO)) return (1);
    if (PlotContourTriangle2D(theElement,CornersOfElem,MP0,MP1,MP2,depth-1,PtrDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   PlotContourQuadrilateral2D - Plot on quadrilateral contourlines (2D coord) with depth

   SYNOPSIS:
   static INT PlotContourQuadrilateral2D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *QP0, COORD *QP1, COORD *QP2,
   COORD *QP3, INT depth, DRAWINGOBJ **PtrDO);

   PARAMETERS:
   .  theElement -
   .  CornersOfElem -
   .  QP0 -
   .  QP1 -
   .  QP2 -
   .  QP3 -
   .  depth -
   .  PtrDO -  the drawing object to draw on

   DESCRIPTION:
   This function plots on quadrilateral contourlines (2D coord) with depth.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PlotContourQuadrilateral2D (ELEMENT *theElement, const COORD **CornersOfElem, const COORD *QP0, const COORD *QP1, const COORD *QP2, const COORD *QP3, INT depth, DRAWINGOBJ **PtrDO)
{
  COORD_VECTOR EVP, LocalCoord, MP0, MP1, MP2, MP3, PointMid, Point[4];
  INT i, j, n, min, max;
  long Color;
  DOUBLE v0, v1, v2, v3, vmin, vmax;

  for (i=0; i<DIM; i++)
    EVP[i] = (QP0[i]+QP1[i]+QP2[i]+QP3[i])*0.25;
  if (depth<=0)
  {
    /* get values at the corners */
    if (GlobalToLocal2d(4,CornersOfElem,QP0,LocalCoord)) return (1);
    v0      = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal2d(4,CornersOfElem,QP1,LocalCoord)) return (1);
    v1      = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal2d(4,CornersOfElem,QP2,LocalCoord)) return (1);
    v2      = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal2d(4,CornersOfElem,QP3,LocalCoord)) return (1);
    v3      = (*EScalar2D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    vmin = MIN(v0,v1); vmin = MIN(vmin,v2); vmin = MIN(vmin,v3);
    vmax = MAX(v0,v1); vmax = MAX(vmax,v2); vmax = MAX(vmax,v3);

    /* store range */
    EScalar2D_minValue = MIN(EScalar2D_minValue,vmin);
    EScalar2D_maxValue = MAX(EScalar2D_maxValue,vmax);

    /* find contours to be plotted */
    for (min=0; min<EScalar2D_numOfContours; min++)
      if (EScalar2D_ContValues[min]>=vmin)
        break;
    for (max=EScalar2D_numOfContours-1; max>=0; max--)
      if (EScalar2D_ContValues[max]<=vmax)
        break;

    /* draw contours */
    for (i=min; i<=max; i++)
    {
      /* set color */
      Color = (long)(EScalar2D_V2C_factor*EScalar2D_ContValues[i]+EScalar2D_V2C_offset);
      Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
      Color = MAX(Color,WOP_OutputDevice->spectrumStart);

      /* calculate points on each side of triangle having the right value */
      n=0;
      if (PointOnLine2D(EScalar2D_ContValues[i],v0,v1,QP0,QP1,Point[n])) n++;
      if (PointOnLine2D(EScalar2D_ContValues[i],v1,v2,QP1,QP2,Point[n])) n++;
      if (PointOnLine2D(EScalar2D_ContValues[i],v2,v3,QP2,QP3,Point[n])) n++;
      if (PointOnLine2D(EScalar2D_ContValues[i],v3,v0,QP3,QP0,Point[n])) n++;

      /* draw */
      switch (n)
      {
      case 1 :
        DO_2c(*PtrDO) = DO_LINE; DO_inc(*PtrDO)
        DO_2l(*PtrDO) = WOP_OutputDevice->black; DO_inc(*PtrDO);
        V2_COPY(Point[0],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[0],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        break;
      case 2 :
        DO_2c(*PtrDO) = DO_LINE; DO_inc(*PtrDO)
        DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
        V2_COPY(Point[0],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[1],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        break;
      case 3 :
        DO_2c(*PtrDO) = DO_POLYLINE; DO_inc(*PtrDO)
        DO_2c(*PtrDO) = 5; DO_inc(*PtrDO)
        DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
        for (j=0; j<DIM; j++)
          PointMid[j] = (Point[0][j]+Point[1][j]+Point[2][j])/3.0;
        V2_COPY(Point[0],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(PointMid,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[1],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(PointMid,DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[2],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        break;
      case 4 :
        DO_2c(*PtrDO) = DO_LINE; DO_inc(*PtrDO)
        DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
        V2_COPY(Point[0],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[2],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);

        DO_2c(*PtrDO) = DO_LINE; DO_inc(*PtrDO)
        DO_2l(*PtrDO) = Color; DO_inc(*PtrDO);
        V2_COPY(Point[1],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
        V2_COPY(Point[3],DO_2Cp(*PtrDO)); DO_inc_n(*PtrDO,2);
      }
    }
  }
  else
  {
    /* find corners of subdivided quadrilaterals */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(QP0[i]+QP1[i]);
      MP1[i] = 0.5*(QP1[i]+QP2[i]);
      MP2[i] = 0.5*(QP2[i]+QP3[i]);
      MP3[i] = 0.5*(QP3[i]+QP0[i]);
    }
    if (PlotContourQuadrilateral2D(theElement,CornersOfElem,QP0,MP0,EVP,MP3,depth-1,PtrDO)) return (1);
    if (PlotContourQuadrilateral2D(theElement,CornersOfElem,MP0,QP1,MP1,EVP,depth-1,PtrDO)) return (1);
    if (PlotContourQuadrilateral2D(theElement,CornersOfElem,EVP,MP1,QP2,MP2,depth-1,PtrDO)) return (1);
    if (PlotContourQuadrilateral2D(theElement,CornersOfElem,EVP,MP2,QP3,MP3,depth-1,PtrDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   EW_EScalar2D	- C(olor)C(ontour) plot

   SYNOPSIS:
   static INT EW_EScalar2D (ELEMENT *theElement, DRAWINGOBJ *theDO);

   PARAMETERS:
   .  theElement -
   .  theDO -

   DESCRIPTION:
   This function plots C(olor)C(ontour).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_EScalar2D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, n;
  const COORD *x[MAX_CORNERS_OF_ELEM];
  DRAWINGOBJ *p, *range;

  n = CORNERS_OF_ELEM(theElement);

  /* get coordinates of corners of the element */
  for (i=0; i<n; i++)
    x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* draw polygon with depth */
  p = theDO;
  EScalar2D_minValue = MAX_D; EScalar2D_maxValue = -MAX_D;
  DO_2c(theDO) = DO_RANGE; DO_inc(theDO); range = theDO; DO_inc_n(theDO,2);
  switch (EScalar2D_mode)
  {
  case PO_COLOR :
    if (n==TRIANGLE)
    {
      if (PlotColorTriangle2D(theElement,x,x[0],x[1],x[2],EScalar2D_depth,&theDO)) return (1);
    }
    else
    {
      if (PlotColorQuadrilateral2D(theElement,x,x[0],x[1],x[2],x[3],EScalar2D_depth,&theDO)) return (1);
    }
    break;
  case PO_CONTOURS_EQ :
    if (n==TRIANGLE)
    {
      if (PlotContourTriangle2D(theElement,x,x[0],x[1],x[2],EScalar2D_depth,&theDO)) return (1);
    }
    else
    {
      if (PlotContourQuadrilateral2D(theElement,x,x[0],x[1],x[2],x[3],EScalar2D_depth,&theDO)) return (1);
    }
    break;
  default :
    return (1);
  }

  DO_2c(theDO) = DO_NO_INST;
  DO_2C(range) = EScalar2D_minValue; DO_inc(range);
  DO_2C(range) = EScalar2D_maxValue;

#ifdef __DO_HEAP_USED__
  n = (INT)theDO - (INT)p;
  Heap_Used_Min = MIN(Heap_Used_Min,n);
  Heap_Used_Max = MAX(Heap_Used_Max,n);
#endif

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EScalar2D_FR - Initialize for findrange of scalar plot

   SYNOPSIS:
   static INT EW_PreProcess_EScalar2D_FR (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for findrange of scalar plot.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_EScalar2D_FR (PICTURE *thePicture, WORK *theWork)
{
  if (EW_PreProcess_EScalar2D (thePicture,theWork))
    return (1);

  /* reset min and max values */
  GEN_FR_put = W_FINDRANGE_WORK(theWork)->put;
  GEN_FR_min = MAX_D;
  GEN_FR_max = -MAX_D;

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EVector2D - Initialize for vector plot 2D

   SYNOPSIS:
   static INT EW_PreProcess_EVector2D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for vector plot 2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_EVector2D (PICTURE *thePicture, WORK *theWork)
{
  struct ElemVectorPlotObj2D *theEvpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theEvpo = &(PIC_PO(thePicture)->theEvpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  /* set value->length fct, eval fct */
  if (theEvpo->max < SMALL_D)
    if (W_ID(theWork) != FINDRANGE_WORK)
    {
      UserWrite("maxValue has to be larger than zero\n");
      return (1);
    }
  EVector_rastersize        = theEvpo->RasterSize;
  EVector_cutvector         = theEvpo->CutVectors;
  EVector_EvalFct           = theEvpo->EvalFct->PlotProc;
  EVector_V2L_factor        = EVector_rastersize/theEvpo->max;                                  /* scale length of vectors			*/
  EVector_CutLenFactor  = theEvpo->CutLenFactor;
  EVector_ColorCut          = theOD->red;
  EVector2D_ColorNormal = theOD->black;

  /* mark suface elements on boundary */
  if (MarkElements_MGS(theMG,0,CURRENTLEVEL(theMG))) return (1);

  /* prepare evaluation routine */
  if (theEvpo->EvalFct->PreprocessProc!=NULL)
    if ((*theEvpo->EvalFct->PreprocessProc)(theMG))
      return (1);

  return (0);
}

/****************************************************************************/
/*
   FindRasterPoints2D - Find rasterpoints in 2D

   SYNOPSIS:
   static INT FindRasterPoints2D (COORD RasterSize, COORD **Polygon,
   INT Number, COORD_VECTOR *RasterPoints, INT *RPNumber);

   PARAMETERS:
   .  RssterSize -
   .  Polygon -
   .  Number -
   .  RasterPoints -
   .  RPNumber -

   DESCRIPTION:
   This function finds rasterpoints in 2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT FindRasterPoints2D (COORD RasterSize, const COORD **Polygon, INT Number, COORD_VECTOR *RasterPoints, INT *RPNumber)
{
  INT i, j, k, i0, i1, j0, j1, c0, c1;
  COORD xmin, xmax, ymin, ymax;
  COORD diff[MAX_POINTS_OF_POLY][2], test[2];

  *RPNumber = 0;
  if (Number<2) return (0);

  xmin = ymin = MAX_C;
  xmax = ymax = -MAX_C;
  for (i=0; i<Number; i++)
  {
    xmin = MIN(xmin,Polygon[i][0]);
    xmax = MAX(xmax,Polygon[i][0]);
    ymin = MIN(ymin,Polygon[i][1]);
    ymax = MAX(ymax,Polygon[i][1]);
    diff[i][0] = Polygon[(i+1)%Number][0] - Polygon[i][0];
    diff[i][1] = Polygon[(i+1)%Number][1] - Polygon[i][1];
  }
  i0 = (INT)ceil(xmin/RasterSize);
  i1 = (INT)floor(xmax/RasterSize);
  j0 = (INT)ceil(ymin/RasterSize);
  j1 = (INT)floor(ymax/RasterSize);

  for (i=i0; i<=i1; i++)
    for (j=j0; j<=j1; j++)
    {
      c0 = c1 = 0;
      for (k=0; k<Number; k++)
      {
        test[0] = RasterSize*(COORD)(i) - Polygon[k][0];
        test[1] = RasterSize*(COORD)(j) - Polygon[k][1];
        if (diff[k][0]*test[1]>=diff[k][1]*test[0]) c0++;
        if (diff[k][0]*test[1]<=diff[k][1]*test[0]) c1++;
      }
      if (c0==Number || c1==Number)
      {
        RasterPoints[*RPNumber][0] = RasterSize*(COORD)(i);
        RasterPoints[*RPNumber][1] = RasterSize*(COORD)(j);
        (*RPNumber)++;
      }
      if (*RPNumber>=RASTERPOINTS_MAX)
        return (0);
    }

  return (0);
}

/****************************************************************************/
/*
   EW_EVector2D	- Evaluate elements for vector drawing

   SYNOPSIS:
   static INT EW_EVector2D (ELEMENT *theElement, DRAWINGOBJ *theDO);

   PARAMETERS:
   .  theElement -
   .  theDO -

   DESCRIPTION:
   This function evaluates elements for vector drawing.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_EVector2D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, nr;
  COORD_VECTOR LocalCoord, Poly[MAX_POINTS_OF_POLY], RasterPoint[RASTERPOINTS_MAX];
  const COORD *x[MAX_CORNERS_OF_ELEM];
  COORD norm;
  long Color;
  DOUBLE min, max;
  DOUBLE_VECTOR Arrow;

  /* get coordinates of corners of the element and their z coordinates in cut system */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));
    V2_COPY(x[i],Poly[i]);
  }

  /* get arrows with rastersize */
  if (FindRasterPoints2D(EVector_rastersize,x,i,RasterPoint,&nr)) return (1);

  /* handle arrows */
  min = MAX_D; max = -MAX_D;
  for (i=0; i<nr; i++)
  {
    if (GlobalToLocal2d(CORNERS_OF_ELEM(theElement),x,RasterPoint[i],LocalCoord)) return (1);
    (*EVector_EvalFct)(theElement,x,LocalCoord,Arrow);
    V2_SCALE(EVector_V2L_factor,Arrow)

    /* find color and size of arrow, define its endpoint on the cutplane */
    V2_EUKLIDNORM(Arrow,norm)
    max = MAX(max,norm/EVector_V2L_factor); min = MIN(min,norm/EVector_V2L_factor);

    if ((norm>EVector_rastersize*EVector_CutLenFactor) && EVector_cutvector)
    {
      Color = EVector_ColorCut;
      V2_SCALE(EVector_rastersize*EVector_CutLenFactor/norm,Arrow)
    }
    else
      Color = EVector2D_ColorNormal;
    V2_ADD(RasterPoint[i],Arrow,Arrow);

    /* draw arrow */
    DO_2c(theDO) = DO_ARROW; DO_inc(theDO)
    DO_2l(theDO) = Color; DO_inc(theDO);
    V2_COPY(RasterPoint[i],DO_2Cp(theDO)); DO_inc_n(theDO,2);
    V2_COPY(Arrow,DO_2Cp(theDO)); DO_inc_n(theDO,2);
  }

  if (nr>0)
  {
    /* store range */
    DO_2c(theDO) = DO_RANGE; DO_inc(theDO)
    DO_2C(theDO) = min; DO_inc(theDO);
    DO_2C(theDO) = max; DO_inc(theDO);
  }

  DO_2c(theDO) = DO_NO_INST;

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EVector2D_FR - Initialize for findrange of vector plot in 2D

   SYNOPSIS:
   static INT EW_PreProcess_EVector2D_FR (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for findrange of vector plot in 2D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_EVector2D_FR (PICTURE *thePicture, WORK *theWork)
{
  if (EW_PreProcess_EVector2D (thePicture,theWork))
    return (1);

  /* reset min and max values */
  GEN_FR_put = W_FINDRANGE_WORK(theWork)->put;
  GEN_FR_min = MAX_D;
  GEN_FR_max = -MAX_D;

  return (0);
}

#endif

/**********************************************************************************************************/
/************************************ Part only for 3D Version ********************************************/
/**********************************************************************************************************/

#ifdef __THREEDIM__

/****************************************************************************/
/*
   GetPolyElemSideISHalfSpace -  Get polygon

   SYNOPSIS:
   static INT GetPolyElemSideISHalfSpace (COORD **Corners,
   COORD *CutZCoord, INT NodeOrder, INT side, COORD_VECTOR *Poly,
   INT *Number);

   PARAMETERS:
   .  Corners -
   .  CutZCoord -
   .  NodeOrder -
   .  side -
   .  Poly -
   .  Number -

   DESCRIPTION:
   This function gets polygon being the intersection of an elementside of an
   tetrahedron an the half space behind cut plane.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT GetPolyElemSideISHalfSpace (COORD **Corners, COORD *CutZCoord, INT NodeOrder, INT side, COORD_VECTOR *Poly, INT *Number)
{
  INT i, j, count1, count2, count3;
  COORD *x[CORNERS_OF_TETRASIDE], Z[CORNERS_OF_TETRASIDE];

  /* get corners of the element side */
  count1 = 0;
  count2 = 0;
  for (i=0; i<CORNERS_OF_TETRASIDE; i++)
  {
    j = SideCornerIndex[side][NodeOrder][i];
    x[i] = Corners[j];
    Z[i] = CutZCoord[j];
    if (Z[i] > 0.0)
      count1++;
    if (Z[i] < 0.0)
      count2++;
  }
  count3 = CORNERS_OF_TETRASIDE-count1-count2;

  /* build the polygon */
  *Number = 0;
  if (count1+count3 >= CORNERS_OF_TETRASIDE)
    return (0);
  for (i=0; i<count2+count3; i++)
    for (j=0; j<DIM; j++)
      Poly[i][j] = x[CORNERS_OF_TETRASIDE-i-1][j];
  *Number = count2+count3;

  if (*Number == CORNERS_OF_TETRASIDE)
    return (0);
  switch (count1)
  {
  case (1) :
    if (*Number != 2) return (1);
    V3_LINCOMB(Z[0]/(Z[0]-Z[1]), x[1], -Z[1]/(Z[0]-Z[1]), x[0], Poly[2])
    V3_LINCOMB(Z[0]/(Z[0]-Z[2]), x[2], -Z[2]/(Z[0]-Z[2]), x[0], Poly[3])
    *Number = 4;
    return (0);
  case (2) :
    if (*Number != 1) return (1);
    V3_LINCOMB(Z[0]/(Z[0]-Z[2]), x[2], -Z[2]/(Z[0]-Z[2]), x[0], Poly[1])
    V3_LINCOMB(Z[1]/(Z[1]-Z[2]), x[2], -Z[2]/(Z[1]-Z[2]), x[1], Poly[2])
    *Number = 3;
    return (0);
  }

  return (1);
}

/****************************************************************************/
/*
   GetLineElemSideISCutPlane - Get line

   SYNOPSIS:
   static INT GetLineElemSideISCutPlane (COORD **Corners,
   COORD *CutZCoord, INT NodeOrder, INT side, COORD_VECTOR *Line,
   INT *Number);

   PARAMETERS:
   .  Corners -
   .  CutZCoord -
   .  NodeOrder -
   .  side -
   .  Line -
   .  Number -

   DESCRIPTION:
   This function gets line being the intersection of an elementside of an
   tetrahedron an the cut plane.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT GetLineElemSideISCutPlane (COORD **Corners, COORD *CutZCoord, INT NodeOrder, INT side, COORD_VECTOR *Line, INT *Number)
{
  INT i, j, c;
  COORD *x[CORNERS_OF_TETRASIDE], Z[CORNERS_OF_TETRASIDE];

  /* get corners of the element side */
  c = 0;
  for (i=0; i<CORNERS_OF_TETRASIDE; i++)
  {
    j = SideCornerIndex[side][NodeOrder][i];
    x[i] = Corners[j];
    Z[i] = CutZCoord[j];
    if (Z[i] > 0.0) c++;
  }

  /* build the polygon */
  switch (c)
  {
  case (1) :
    V3_LINCOMB(Z[0]/(Z[0]-Z[1]), x[1], -Z[1]/(Z[0]-Z[1]), x[0], Line[0])
    V3_LINCOMB(Z[0]/(Z[0]-Z[2]), x[2], -Z[2]/(Z[0]-Z[2]), x[0], Line[1])
    *Number = 2;
    return (0);
  case (2) :
    V3_LINCOMB(Z[0]/(Z[0]-Z[2]), x[2], -Z[2]/(Z[0]-Z[2]), x[0], Line[0])
    V3_LINCOMB(Z[1]/(Z[1]-Z[2]), x[2], -Z[2]/(Z[1]-Z[2]), x[1], Line[1])
    *Number = 2;
    return (0);
  default :
    *Number = 0;
    return (0);
  }
}

/****************************************************************************/
/*
   GetPolyElemISCutPlane - Get polygon

   SYNOPSIS:
   static INT GetPolyElemISCutPlane (COORD **CornerDC,
   COORD *CutZCoord, INT NodeOrder, COORD_VECTOR *Poly, INT *Number);

   PARAMETERS:
   .  CornerDC -
   .  CutZCoord -
   .  NodeOrder -
   .  Poly -
   .  Number -

   DESCRIPTION:
   This function gets polygon being the intersection of an
   tetrahedron and the cut plane.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT GetPolyElemISCutPlane (COORD **CornerDC, COORD *CutZCoord, INT NodeOrder, COORD_VECTOR *Poly, INT *Number)
{
  INT i, j, count1, count2;
  COORD *x[MAX_CORNERS_OF_ELEM], Z[MAX_CORNERS_OF_ELEM];

  /* get corners of the element side */
  count1 = 0;
  count2 = 0;
  for (i=0; i<4; i++)           /* 4=Corners of tetrahedron */
  {
    j = CornerIndex[NodeOrder][i];
    x[i] = CornerDC[j];
    Z[i] = CutZCoord[j];
    if (Z[i] > SMALL_C)
      count1++;
    if (Z[i] < -SMALL_C)
      count2++;
  }

  /* build the polygon */
  *Number = 0;
  switch (count1)
  {
  case (0) :
    switch (count2)
    {
    case (0) :
      return (1);
    case (1) :
      V3_COPY(x[0], Poly[0])
      V3_COPY(x[1], Poly[1])
      V3_COPY(x[2], Poly[2])
      *Number = 3;
      return (0);
    case (2) :
    case (3) :
    case (4) :
      return (0);
    }
  case (1) :
    switch (count2)
    {
    case (0) :
      V3_COPY(x[1], Poly[0])
      V3_COPY(x[2], Poly[1])
      V3_COPY(x[3], Poly[2])
      *Number = 3;
      return (0);
    case (1) :
      V3_COPY(x[1], Poly[0])
      V3_COPY(x[2], Poly[1])
      V3_LINCOMB(Z[0]/(Z[0]-Z[3]), x[3], -Z[3]/(Z[0]-Z[3]), x[0], Poly[2])
      *Number = 3;
      return (0);
    case (2) :
      V3_COPY(x[1], Poly[0])
      V3_LINCOMB(Z[0]/(Z[0]-Z[2]), x[2], -Z[2]/(Z[0]-Z[2]), x[0], Poly[1]);
      V3_LINCOMB(Z[0]/(Z[0]-Z[3]), x[3], -Z[3]/(Z[0]-Z[3]), x[0], Poly[2]);
      *Number = 3;
      return (0);
    case (3) :
      V3_LINCOMB(Z[0]/(Z[0]-Z[1]), x[1], -Z[1]/(Z[0]-Z[1]), x[0], Poly[0])
      V3_LINCOMB(Z[0]/(Z[0]-Z[2]), x[2], -Z[2]/(Z[0]-Z[2]), x[0], Poly[1])
      V3_LINCOMB(Z[0]/(Z[0]-Z[3]), x[3], -Z[3]/(Z[0]-Z[3]), x[0], Poly[2])
      *Number = 3;
      return (0);
    default :
      return (1);
      break;
    }
  case (2) :
    switch (count2)
    {
    case (0) :
      return (0);
    case (1) :
      V3_COPY(x[2], Poly[0])
      V3_LINCOMB(Z[0]/(Z[0]-Z[3]), x[3], -Z[3]/(Z[0]-Z[3]), x[0], Poly[1])
      V3_LINCOMB(Z[1]/(Z[1]-Z[3]), x[3], -Z[3]/(Z[1]-Z[3]), x[1], Poly[2])
      *Number = 3;
      return (0);
    case (2) :
      V3_LINCOMB(Z[0]/(Z[0]-Z[2]), x[2], -Z[2]/(Z[0]-Z[2]), x[0], Poly[0])
      V3_LINCOMB(Z[0]/(Z[0]-Z[3]), x[3], -Z[3]/(Z[0]-Z[3]), x[0], Poly[1])
      V3_LINCOMB(Z[1]/(Z[1]-Z[3]), x[3], -Z[3]/(Z[1]-Z[3]), x[1], Poly[2])
      V3_LINCOMB(Z[1]/(Z[1]-Z[2]), x[2], -Z[2]/(Z[1]-Z[2]), x[1], Poly[3])
      *Number = 4;
      return (0);
    default :
      return (1);
      break;
    }
  case (3) :
    switch (count2)
    {
    case (0) :
      return (0);
    case (1) :
      V3_LINCOMB(Z[0]/(Z[0]-Z[3]), x[3], -Z[3]/(Z[0]-Z[3]), x[0], Poly[0])
      V3_LINCOMB(Z[1]/(Z[1]-Z[3]), x[3], -Z[3]/(Z[1]-Z[3]), x[1], Poly[1])
      V3_LINCOMB(Z[2]/(Z[2]-Z[3]), x[3], -Z[3]/(Z[2]-Z[3]), x[2], Poly[2])
      *Number = 3;
      return (0);
    default :
      return (1);
    }
  case (4) :
    return (0);
  }

  return (1);
}

/****************************************************************************/
/*
   EW_ElementEval3D - Evaluate geometry of TetraHedron

   SYNOPSIS:
   static INT EW_ElementEval3D (ELEMENT *theElement, DRAWINGOBJ *theDO);

   PARAMETERS:
   .  theElement -
   .  theDO -

   DESCRIPTION:
   This function evaluates geometry of TetraHedron.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_ElementEval3D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, j, NodeOrder, n;
  COORD *x[MAX_CORNERS_OF_ELEM], z[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR Polygon[MAX_POINTS_OF_POLY];
  INT Viewable[MAX_SIDES_OF_ELEM];

  DO_2c(theDO) = DO_NO_INST;

  if (ID(theElement)==3231)
    i = 0;

  if (!CUT_CutExisting || CUTMODE(theElement)==CM_BEHIND)
  {
    /* plot full element */

    /* determine viewable sides */
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      Viewable[i] = VIEWABLE(theElement,i);
    if (EE3D_Elem2Plot[PLOT_ALL])
    {
      /* plot only parts lying on the boundary */
      if (OBJT(theElement)==BEOBJ)
      {
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
          if (SIDE(theElement,i)==NULL)
            Viewable[i] = 0;
      }
      else
        return (0);
    }
    else
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
        if (NBELEM(theElement,i) != NULL)
          if (EE3D_Elem2Plot[ECLASS(NBELEM(theElement,i))])
            Viewable[i] = 0;

    /* get coordinates of corners of the element */
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

    /* store viewable sides on drawing obj */
    if (LEVEL(theElement)<EE3D_MaxLevel)
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        if (!Viewable[i]) continue;
        if (EE3D_NoColor[COLOR_LOWER_LEVEL])
        {
          DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = 3; DO_inc(theDO)
        }
        else
        {
          DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = 3; DO_inc(theDO)
          DO_2l(theDO) = EE3D_Color[COLOR_LOWER_LEVEL]; DO_inc(theDO)
        }
        DO_2l(theDO) = EE3D_Color[COLOR_EDGE]; DO_inc(theDO);
        for (j=0; j<CORNERS_OF_TETRASIDE; j++)
        {
          V3_COPY(x[CornerOfSide[i][j]],DO_2Cp(theDO));
          DO_inc_n(theDO,3);
        }
      }
    else
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        if (!Viewable[i]) continue;

        if (EE3D_NoColor[ECLASS(theElement)])
        {
          DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = 3; DO_inc(theDO)
        }
        else
        {
          DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = 3; DO_inc(theDO)
          DO_2l(theDO) = EE3D_Color[ECLASS(theElement)]; DO_inc(theDO);
        }
        DO_2l(theDO) = EE3D_Color[COLOR_EDGE]; DO_inc(theDO);
        for (j=0; j<CORNERS_OF_TETRASIDE; j++)
        {
          V3_COPY(x[CornerOfSide[i][j]],DO_2Cp(theDO));
          DO_inc_n(theDO,3);
        }

        /* inverse if selected */
        if (IsElementSelected(GElem_MG,theElement))
        {
          DO_2c(theDO) = DO_INVERSE_POLYGON; DO_inc(theDO)
          DO_2c(theDO) = 3; DO_inc(theDO)
          for (j=0; j<CORNERS_OF_TETRASIDE; j++)
          {
            V3_COPY(x[CornerOfSide[i][j]],DO_2Cp(theDO));
            DO_inc_n(theDO,3);
          }
        }
      }
  }
  else if (CUTMODE(theElement)==CM_INTERSECT)
  {
    /* plot cutted element */

    /* determine viewable sides */
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      Viewable[i] = VIEWABLE(theElement,i);
    if (EE3D_Elem2Plot[PLOT_ALL])
    {
      /* only sides lying on the boundary are visible */
      if (OBJT(theElement)==BEOBJ)
      {
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
          if (SIDE(theElement,i)==NULL)
            Viewable[i] = 0;
      }
      else
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
          Viewable[i] = 0;
    }
    else
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
        if (NBELEM(theElement,i) != NULL)
          if (EE3D_Elem2Plot[ECLASS(NBELEM(theElement,i))])
            Viewable[i] = 0;

    /* get coordinates of corners of the element and their z coordinates in cut system */
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));
      V3_TRAFO4_SC(x[i],CutTrafo,z[i])
    }

    /* get node order */
    NodeOrder = NORDER(theElement);

    /* plot that parts of the viewable sides of the element lying behind cut plane */
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      if (!Viewable[i]) continue;

      /* determine polygon arising from intersection of triangle with half space behind cut plane */
      if (GetPolyElemSideISHalfSpace (x,z,NodeOrder,i,Polygon,&n))
        return (1);
      if (n<=2) continue;

      /* store on drawing object */
      if (LEVEL(theElement)<EE3D_MaxLevel)
      {
        if (EE3D_NoColor[COLOR_LOWER_LEVEL])
        {
          DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = n; DO_inc(theDO)
        }
        else
        {
          DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = n; DO_inc(theDO)
          DO_2l(theDO) = EE3D_Color[COLOR_LOWER_LEVEL]; DO_inc(theDO);
        }
        DO_2l(theDO) = EE3D_Color[COLOR_EDGE]; DO_inc(theDO);
        for (j=0; j<n; j++)
        {
          V3_COPY(Polygon[j],DO_2Cp(theDO));
          DO_inc_n(theDO,3);
        }
      }
      else
      {
        if (EE3D_NoColor[ECLASS(theElement)])
        {
          DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = n; DO_inc(theDO)
        }
        else
        {
          DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
          DO_2c(theDO) = n; DO_inc(theDO)
          DO_2l(theDO) = EE3D_Color[ECLASS(theElement)]; DO_inc(theDO);
        }
        DO_2l(theDO) = EE3D_Color[COLOR_EDGE]; DO_inc(theDO);
        for (j=0; j<n; j++)
        {
          V3_COPY(Polygon[j],DO_2Cp(theDO));
          DO_inc_n(theDO,3);
        }

        /* inverse if selected */
        if (IsElementSelected(GElem_MG,theElement))
        {
          DO_2c(theDO) = DO_INVERSE_POLYGON; DO_inc(theDO)
          DO_2c(theDO) = n; DO_inc(theDO)
          for (j=0; j<n; j++)
          {
            V3_COPY(Polygon[j],DO_2Cp(theDO));
            DO_inc_n(theDO,3);
          }
        }
      }
    }

    /* plot intersection of element with cut plane if */
    if (CUT_CutAtFront)
    {
      if (GetPolyElemISCutPlane(x,z,NodeOrder,Polygon,&n))
        return (1);

      /* store on drawing object */
      if (n>2)
        if (LEVEL(theElement)<EE3D_MaxLevel)
        {
          switch (EE3D_NoColor[COLOR_LOWER_LEVEL] | (EE3D_NoColor[COLOR_CUT_EDGE]<<1))
          {
          case 0 :
            DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[COLOR_LOWER_LEVEL]; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[COLOR_CUT_EDGE]; DO_inc(theDO)
            break;
          case 1 :
            DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[COLOR_CUT_EDGE]; DO_inc(theDO)
            break;
          case 2 :
            DO_2c(theDO) = DO_POLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[COLOR_LOWER_LEVEL]; DO_inc(theDO)
            break;
          case 3 :
            DO_2c(theDO) = DO_ERASE_POLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            break;
          }
          for (j=0; j<n; j++)
          {
            V3_COPY(Polygon[j],DO_2Cp(theDO));
            DO_inc_n(theDO,3);
          }
        }
        else
        {
          switch (EE3D_NoColor[ECLASS(theElement)] | (EE3D_NoColor[COLOR_CUT_EDGE]<<1))
          {
          case 0 :
            DO_2c(theDO) = DO_SURRPOLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[ECLASS(theElement)]; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[COLOR_CUT_EDGE]; DO_inc(theDO)
            break;
          case 1 :
            DO_2c(theDO) = DO_ERASE_SURRPOLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[COLOR_CUT_EDGE]; DO_inc(theDO)
            break;
          case 2 :
            DO_2c(theDO) = DO_POLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            DO_2l(theDO) = EE3D_Color[COLOR_LOWER_LEVEL]; DO_inc(theDO)
            break;
          case 3 :
            DO_2c(theDO) = DO_ERASE_POLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            break;
          }
          for (j=0; j<n; j++)
          {
            V3_COPY(Polygon[j],DO_2Cp(theDO));
            DO_inc_n(theDO,3);
          }


          /* inverse if selected */
          if (IsElementSelected(GElem_MG,theElement))
          {
            DO_2c(theDO) = DO_INVERSE_POLYGON; DO_inc(theDO)
            DO_2c(theDO) = n; DO_inc(theDO)
            for (j=0; j<n; j++)
            {
              V3_COPY(Polygon[j],DO_2Cp(theDO));
              DO_inc_n(theDO,3);
            }
          }
        }
    }
  }
  DO_2c(theDO) = DO_NO_INST;

  return (0);
}

/****************************************************************************/
/*
   EW_ECutBnd3D	-  Evaluate geometry of TetraHedron

   SYNOPSIS:
   static INT EW_ECutBnd3D (ELEMENT *theElement, DRAWINGOBJ *theDO);

   PARAMETERS:
   .  theElement -
   .  theDO -

   DESCRIPTION:
   This function evaluates geometry of TetraHedron.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_ECutBnd3D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, NodeOrder, n;
  COORD *x[MAX_CORNERS_OF_ELEM], z[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR Line[2];

  DO_2c(theDO) = DO_NO_INST;

  if (CUTMODE(theElement)==CM_INTERSECT && OBJT(theElement)==BEOBJ && CUT_CutAtFront)
  {
    /* get coordinates of corners of the element and their z coordinates in cut system */
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));
      V3_TRAFO4_SC(x[i],CutTrafo,z[i])
    }

    /* get node order */
    NodeOrder = NORDER(theElement);

    /* plot boundary side intersection with cut-plane */
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      if (SIDE(theElement,i)==NULL) continue;

      /* determine line arising from intersection of side with cut plane */
      if (GetLineElemSideISCutPlane (x,z,NodeOrder,i,Line,&n))
        return (1);
      if (n<2) continue;

      /* store line on drawing object */
      DO_2c(theDO) = DO_LINE; DO_inc(theDO)
      DO_2l(theDO) = EE3D_Color[COLOR_CUT_EDGE]; DO_inc(theDO);
      V3_COPY(Line[0],DO_2Cp(theDO)); DO_inc_n(theDO,3);
      V3_COPY(Line[1],DO_2Cp(theDO)); DO_inc_n(theDO,3);
    }
  }
  DO_2c(theDO) = DO_NO_INST;

  return (0);
}

/****************************************************************************/
/*
   EW_FindElement3D - Find element in 3D drawing object

   SYNOPSIS:
   static INT EW_SelectElement3D (DRAWINGOBJ *q);

   PARAMETERS:
   .  q - the drawing object

   DESCRIPTION:
   This function finds element in 3D drawing object.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_SelectElement3D (DRAWINGOBJ *q)
{
  INT j, n, end, found;
  COORD help[3];
  COORD_POINT point[MAX_POINTS_OF_POLY];

  end = found = 0;
  while (!end)
  {
    switch (DO_2c(q))
    {
    case DO_NO_INST :
      end = 1;
      break;
    case DO_RANGE :
      DO_inc_RANGE(q);
      break;
    case DO_LINE :
      DO_inc_LINE(q,3);
      break;
    case DO_ARROW :
      DO_inc_ARROW(q,3);
      break;
    case DO_INVERSE_LINE :
      DO_inc_INVERSE_LINE(q,3);
      break;
    case DO_POLYLINE :
      DO_inc_POLYLINE(q,3);
      break;
    case DO_TEXT :
      DO_inc_TEXT(q,3);
      break;
    case DO_POLYMARK :
      DO_inc_POLYMARK(q,3);
      break;
    case DO_POLYGON :
    case DO_ERASE_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc_n(q,2)
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      found |= PointInPolygon(point,j,FE3D_MousePos);
      break;
    case DO_INVERSE_POLYGON :
    case DO_ERASE_POLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc(q)
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      found |= PointInPolygon(point,j,FE3D_MousePos);
      break;
    case DO_SURRPOLYGON :
      DO_inc(q)
      n = DO_2c(q); DO_inc_n(q,3)
      for (j=0; j<n; j++)
      {
        V3_TRAFOM4_V3(DO_2Cp(q),ObsTrafo,help); DO_inc_n(q,3);
        (*OBS_ProjectProc)(help,point+j);
      }
      found |= PointInPolygon(point,j,FE3D_MousePos);
      break;
    default :
      return (1);
    }

  }

  /* if found ... */
  if (found)
  {
    /* put in/delete from selection list */
    if (SELECTIONMODE(WOP_MG)!=elementSelection)
      ClearSelection(WOP_MG);
    if (AddElementToSelection(WOP_MG,WOP_Element) == GM_ERROR)
      if (RemoveElementFromSelection(WOP_MG,WOP_Element) == GM_ERROR)
        return (1);

    /* plot part lying in front */
    if (EW_ElementEval3D(WOP_Element,WOP_DrawingObject)) return (1);
    if (Draw3D(WOP_DrawingObject)) return (1);
    WOP_EW_GetNextElementProc       = EW_GetNextElement_vert_fw_up;
    WOP_GEN_ExecuteProc             = Draw3D;
  }

  return (0);
}

/****************************************************************************/
/*
   CalcViewableSides - Determination of viewable tetrahedra sides

   SYNOPSIS:
   static void CalcViewableSides (ELEMENT *theElement);

   PARAMETERS:
   .  theElement -

   DESCRIPTION:
   This function determines viewable tetrahedra sides.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void CalcViewableSides (ELEMENT *theElement)
{
  COORD_VECTOR Vector, Vector01, Vector02, Vector03, ViewDirection;
  INT Viewablility;
  INT i;
  COORD *x[MAX_CORNERS_OF_ELEM];
  COORD ScalarPrd;

  /* load corners of the element */
  for( i=0; i<CORNERS_OF_ELEM(theElement); i++)
    x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* set view direction */
  if (OBS_Perspective == YES)
  {
    Viewablility = 0;
    for( i=0; i<SIDES_OF_ELEM(theElement); i++ )
    {
      V3_SUBTRACT(VO_VP(OE_ViewedObj),x[CornerOfSide[i][0]],ViewDirection)
      V3_SUBTRACT(x[CornerOfSide[i][1]],x[CornerOfSide[i][0]],Vector01)
      V3_SUBTRACT(x[CornerOfSide[i][2]],x[CornerOfSide[i][0]],Vector02)
      V3_SUBTRACT(x[ OppositeCorner[i]],x[CornerOfSide[i][0]],Vector03)
      V3_VECTOR_PRODUCT(Vector01,Vector02,Vector)
      V3_SCALAR_PRODUCT(Vector,Vector03,ScalarPrd)
      if (ScalarPrd>0)
        V3_SCALE(-1.0, Vector)
        if (Vector[0]*ViewDirection[0]+Vector[1]*ViewDirection[1]+Vector[2]*ViewDirection[2]>0)
          Viewablility |= (1<<i);
    }
  }
  else
  {
    V3_SUBTRACT(VO_VP(OE_ViewedObj),VO_VT(OE_ViewedObj),ViewDirection);
    Viewablility = 0;
    for( i=0; i<SIDES_OF_ELEM(theElement); i++ )
    {
      V3_SUBTRACT(x[CornerOfSide[i][1]], x[CornerOfSide[i][0]], Vector01)
      V3_SUBTRACT(x[CornerOfSide[i][2]], x[CornerOfSide[i][0]], Vector02)
      V3_SUBTRACT(x[ OppositeCorner[i]], x[CornerOfSide[i][0]], Vector03)
      V3_VECTOR_PRODUCT(Vector01, Vector02, Vector)
      V3_SCALAR_PRODUCT(Vector,Vector03,ScalarPrd)
      if (ScalarPrd>0)
        V3_SCALE(-1.0, Vector)
        if (Vector[0]*ViewDirection[0]+Vector[1]*ViewDirection[1]+Vector[2]*ViewDirection[2]>0)
          Viewablility |= (1<<i);
    }
  }
  SETVSIDES(theElement,Viewablility);
}

/****************************************************************************/
/*
   CalcViewableSidesOnGrid - Determination of viewable tetrahedra sides

   SYNOPSIS:
   static void CalcViewableSidesOnGrid (GRID *theGrid);

   PARAMETERS:
   theGrid - pointer to grid

   DESCRIPTION:
   This function determines viewable tetrahedra sides.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void CalcViewableSidesOnGrid (GRID *theGrid)
{
  ELEMENT *theElement, *theNeighbor;
  INT i, j;

  /* calc viewable sides for each element */
  for (theElement=theGrid->elements; theElement!= NULL; theElement=SUCCE(theElement))
    CalcViewableSides(theElement);

  /* make the viewable sides consistent */
  for (theElement=theGrid->elements; theElement!= NULL; theElement=SUCCE(theElement))
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      if ((theNeighbor=NBELEM(theElement,i))==NULL) continue;
      if (ID(theElement) < ID(theNeighbor))
      {
        for(j=0; j<SIDES_OF_ELEM(theElement); j++)
          if (NBELEM(theNeighbor,j)==theElement)
            break;
        if (VIEWABLE(theElement,i) == VIEWABLE(theNeighbor,j))
          if (VIEWABLE(theElement,i))
            SETVSIDES(theElement,VSIDES(theElement)&(~(1<<i)));
          else
            SETVSIDES(theElement,VSIDES(theElement)|(1<<i));
      }
    }

}

/****************************************************************************/
/*
   OrderSons - Order elements with respect to view orientation on all levels

   SYNOPSIS:
   static INT OrderSons (ELEMENT **table,ELEMENT *theElement);

   PARAMETERS:
   .  table -
   .  theElement -

   DESCRIPTION:
   This function orders elements with respect to view orientation on all levels.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT OrderSons (ELEMENT **table,ELEMENT *theElement)
{
  INT i, j, Count, nsons;
  INT LastShellBegin, NewShellBegin, ActualPosition;
  ELEMENT *NbElement, *SonElement, *SonList[MAX_SONS];

  /* get son list (not stored in element) */
  if (GetSons(theElement,SonList)!=0) return(1);

  /* init list and numbers */
  LastShellBegin = 0;
  ActualPosition = 0;
  nsons = NSONS(theElement);
  for  (i=0; i<nsons; i++)
  {
    SonElement = SonList[i];

    /* count how many neighbor-sons are overlapped by SonElement */
    Count = 0;
    for (j=0; j<SIDES_OF_ELEM(SonElement); j++)
      if (NBELEM(SonElement,j)!=NULL && EFATHER(NBELEM(SonElement,j))==theElement && (!VIEWABLE(SonElement,j))) Count++;
    if (Count)
    {
      SETCOUNT(SonElement,Count);
      SETUSED(SonElement,0);
    }
    else
    {
      table[ActualPosition++] = SonElement;
      SETUSED(SonElement,1);
    }
  }
  NewShellBegin = ActualPosition;

  /* create list */
  while (ActualPosition<nsons)
  {
    /* create a new shell */
    for (i=LastShellBegin; i<NewShellBegin; i++)
    {
      for (j=0; j<SIDES_OF_ELEM(table[i]); j++)
      {
        if ((NbElement=NBELEM(table[i],j))==NULL) continue;
        if (EFATHER(NbElement)!=theElement || USED(NbElement)) continue;
        if ((Count = COUNT(NbElement)-1)==0)
        {
          table[ActualPosition++] = NbElement;
          SETUSED(NbElement,1);
        }
        else
          SETCOUNT(NbElement,Count);
      }
    }

    /* set shell pointers */
    LastShellBegin = NewShellBegin;
    NewShellBegin = ActualPosition;
  }
  return (0);
}

/****************************************************************************/
/*
   CalcCrossingPoint - Calculate crossing point of two lines

   SYNOPSIS:
   static INT CalcCrossingPoint (COORD_POINT P1, COORD_POINT P2,
   COORD_POINT P3, COORD_POINT P4, COORD *alpha, COORD *beta);

   PARAMETERS:
   .  P1 -
   .  P2 -
   .  P3 -
   .  P4 -
   .  alpha -
   .  beta -

   DESCRIPTION:
   This function calculates crossing point of two lines.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT CalcCrossingPoint (COORD_POINT P1, COORD_POINT P2, COORD_POINT P3, COORD_POINT P4, COORD *alpha, COORD *beta)
{
  INT flags1;
  COORD determinante, c1, c2;

  /* check if one endpoint of line0 coincide with one endpoint of line1 */
  if (ABS(P1.x - P3.x)<SMALL_C && ABS(P1.y - P3.y)<SMALL_C) return(0);
  if (ABS(P1.x - P4.x)<SMALL_C && ABS(P1.y - P4.y)<SMALL_C) return(0);
  if (ABS(P2.x - P3.x)<SMALL_C && ABS(P2.y - P3.y)<SMALL_C) return(0);
  if (ABS(P2.x - P4.x)<SMALL_C && ABS(P2.y - P4.y)<SMALL_C) return(0);

  flags1 = 0;
  if (ABS(P1.x - P2.x)<SMALL_C) flags1 |= 1;
  if (ABS(P1.y - P2.y)<SMALL_C) flags1 |= 2;
  if (ABS(P3.x - P4.x)<SMALL_C) flags1 |= 4;
  if (ABS(P3.y - P4.y)<SMALL_C) flags1 |= 8;

  switch (flags1)
  {
  case (0) :
    /* the natural case */
    determinante = (P2.y-P1.y)*(P4.x-P3.x) - (P2.x-P1.x)*(P4.y-P3.y);
    if (ABS(determinante)<SMALL_C)
    {
      /* the lines are parallel */
      /* check if P1 (or P2) is on line1 */
      c1 = (P1.y-P3.y)/(P4.y-P3.y);
      c2 = (P2.y-P3.y)/(P4.y-P3.y);
      if (ABS((1.0-c1)*P3.x + c1*P4.x - P1.x)<SMALL_C)
      {
        if (((c1<=0.0) && (c2<=0.0)) || ((1.0<=c1) && (1.0<=c2)))
          return (0);
        if (0.0<c1 && c1<1.0)
        {
          *beta = c1;
          *alpha = 0.0;
          return (1);
        }
        if (0.0<c2 && c2<1.0)
        {
          *beta = c2;
          *alpha = 1.0;
          return (1);
        }

        /* line1 is contained in line0: calculate position of P3 on line0 */
        c1 = (P3.y-P1.y)/(P2.y-P1.y);
        if (c1<=0.0 || 1.0<=c1)
          /* this is impossible! */
          return (0);
        *alpha = c1;
        *beta = 0.0;
      }
      return (0);
    }
    else
    {
      /* the lines are not parallel */
      *alpha = ((P4.x-P3.x)*(P3.y-P1.y) - (P4.y-P3.y)*(P3.x-P1.x))/determinante;
      *beta  = ((P2.x-P1.x)*(P3.y-P1.y) - (P2.y-P1.y)*(P3.x-P1.x))/determinante;
      if (0.0<*alpha && *alpha<1.0 && 0.0<*beta && *beta<1.0)
        return (1);
      return (0);
    }
  case (1) :
    /* line0 is vertical */
    *beta = (P1.x-P3.x)/(P4.x-P3.x);
    *alpha = (((1.0-(*beta))*P3.y + (*beta)*P4.y)-P1.y)/(P2.y-P1.y);
    if (0.0<*alpha && *alpha<1.0 && 0.0<*beta && *beta<1.0)
      return (1);
    return (0);
  case (2) :
    /* line0 is horizontal */
    *beta = (P1.y-P3.y)/(P4.y-P3.y);
    *alpha = (((1.0-(*beta))*P3.x + (*beta)*P4.x)-P1.x)/(P2.x-P1.x);
    if (0.0<*alpha && *alpha<1.0 && 0.0<*beta && *beta<1.0)
      return (1);
    return (0);
  case (3) :
    /* line0 is degenerated */
    *beta = (P1.y-P3.y)/(P4.y-P3.y);
    if (ABS((1.0-(*beta))*P3.x + (*beta)*P4.x - P1.x)<SMALL_C)
      if (0.0<*beta && *beta<1.0)
      {
        *alpha = 0.5;
        return (1);
      }
    return (0);
  case (4) :
    /* line1 is vertical */
    *alpha = (P3.x-P1.x)/(P2.x-P1.x);
    *beta  = (((1.0-(*alpha))*P1.y + (*alpha)*P2.y)-P3.y)/(P4.y-P3.y);
    if (0.0<*alpha && *alpha<1.0 && 0.0<*beta && *beta<1.0)
      return (1);
    return (0);
  case (5) :
    /* both lines are vertical */
    if (ABS(P1.x - P3.x)<SMALL_C)
    {
      c1 = (P1.y-P3.y)/(P4.y-P3.y);
      if (0.0<c1 && c1<1.0)
      {
        *alpha = 0.0;
        *beta  = c1;
        return (1);
      }
      c2 = (P2.y-P3.y)/(P4.y-P3.y);
      if (0.0<c2 && c2<1.0)
      {
        *alpha = 1.0;
        *beta  = c2;
        return (1);
      }
      c1 = (P3.y-P1.y)/(P2.y-P1.y);
      if (0.0<c1 && c1<1.0)
      {
        *alpha = c1;
        *beta  = 0.0;
        return (1);
      }
      c2 = (P4.y-P1.y)/(P2.y-P1.y);
      if (0.0<c2 && c2<1.0)
      {
        *alpha = c2;
        *beta  = 1.0;
        return (1);
      }
      /* impossible case */
      return (0);
    }
    return (0);
  case (6) :
    /* line0 is horizontal, line1 is vertical */
    *alpha = (P3.x-P1.x)/(P2.x-P1.x);
    *beta  = (P1.y-P3.y)/(P4.y-P3.y);
    if (0.0<*alpha && *alpha<1.0 && 0.0<*beta && *beta<1.0)
      return (1);
    return (0);
  case (7) :
    /* line0 is degenerated, line1 is vertical */
    if (ABS(P1.x - P3.x)<SMALL_C)
    {
      *alpha = 0.5;
      *beta  = (P1.y-P3.y)/(P4.y-P3.y);
      if ( 0.0<*beta && *beta<1.0)
        return (1);
    }
    return (0);
  case (8) :
    /* line1 is horizontal */
    *alpha = (P3.y-P1.y)/(P2.y-P1.y);
    *beta  = (((1.0-(*alpha))*P1.x + (*alpha)*P2.x)-P3.x)/(P4.x-P3.x);
    if (0.0<*alpha && *alpha<1.0 && 0.0<*beta && *beta<1.0)
      return (1);
    return (0);
  case (9) :
    /* line0 is vertical, line1 is horizontal */
    *alpha = (P3.y-P1.y)/(P2.y-P1.y);
    *beta  = (P1.x-P3.x)/(P4.x-P3.x);
    if (0.0<*alpha && *alpha<1.0 && 0.0<*beta && *beta<1.0)
      return (1);
    return (0);
  case (10) :
    /* both lines are horizontal */
    if (ABS(P1.y - P3.y)<SMALL_C)
    {
      c1 = (P1.x-P3.x)/(P4.x-P3.x);
      if (0.0<c1 && c1<1.0)
      {
        *alpha = 0.0;
        *beta  = c1;
        return (1);
      }
      c2 = (P2.x-P3.x)/(P4.x-P3.x);
      if (0.0<c2 && c2<1.0)
      {
        *alpha = 1.0;
        *beta  = c2;
        return (1);
      }
      c1 = (P3.x-P1.x)/(P2.x-P1.x);
      if (0.0<c1 && c1<1.0)
      {
        *alpha = c1;
        *beta  = 0.0;
        return (1);
      }
      c2 = (P4.x-P1.x)/(P2.x-P1.x);
      if (0.0<c2 && c2<1.0)
      {
        *alpha = c2;
        *beta  = 1.0;
        return (1);
      }
      /* impossible case */
      return (0);
    }
    return (0);
  case (11) :
    /* line0 is degenerated, line1 is horizontal */
    if (ABS(P1.y - P3.y)<SMALL_C)
    {
      *alpha = 0.5;
      *beta  = (P1.x-P3.x)/(P4.x-P3.x);
      if ( 0.0<*beta && *beta<1.0)
        return (1);
    }
    return (0);
  case (12) :
    /* line1 is degenerated */
    *alpha = (P3.y-P1.y)/(P2.y-P1.y);
    if (ABS((1.0-(*alpha))*P1.x + (*alpha)*P2.x - P3.x)<SMALL_C)
      if (0.0<*alpha && *alpha<1.0)
      {
        *beta = 0.5;
        return (1);
      }
    return (0);
  case (13) :
    /* line0 is vertical, lin1 is degenerated */
    if (ABS(P1.x - P3.x)<SMALL_C)
    {
      *alpha = (P3.y-P1.y)/(P2.y-P1.y);
      *beta  = 0.5;
      if ( 0.0<*alpha && *alpha<1.0)
        return (1);
    }
    return (0);
  case (14) :
    /* line0 is horizontal, lin1 is degenerated */
    if (ABS(P1.y - P3.y)<SMALL_C)
    {
      *alpha = (P3.x-P1.x)/(P2.x-P1.x);
      *beta  = 0.5;
      if ( 0.0<*alpha && *alpha<1.0)
        return (1);
    }
    return (0);
  case (15) :
    /* both lines are degenerated */
    if (ABS(P1.x - P3.x)<SMALL_C && ABS(P1.y - P3.y)<SMALL_C)
    {
      *alpha = 0.5;
      *beta  = 0.5;
      return (1);
    }
    return (0);
  }
  return (-1);
}

/****************************************************************************/
/*
   CompareTriangles - Test, whether an elements side hides another elements side

   SYNOPSIS:
   static INT CompareTriangles (COORD_VECTOR Triangle[2][3],
   COORD_POINT ScreenPoint[2][3]);

   PARAMETERS:
   .  Triangle[2][3] -
   .  ScreenPoint[2][3] -

   DESCRIPTION:
   This function tests, whether an elements side hides another elements side.

   RETURN VALUE:
   INT
   .n     0 if triangles do not hide each other
   .n     1 if triangle 0 hides triangle 1
   .n     -1 if triangle 1 hides triangle 0.
 */
/****************************************************************************/

static INT CompareTriangles (COORD_VECTOR Triangle[2][3], COORD_POINT ScreenPoint[2][3])
{
  COORD alpha, beta, z[2];
  COORD lambda[3], rhs[3], M[9], Inverse[9];
  int i, i1, j, j1;

  /* decide if sides of triangle are crossing */
  for (i=0; i<CORNERS_OF_TETRASIDE; i++)
  {
    i1 = (i+1)%CORNERS_OF_TETRASIDE;
    for( j=0; j<CORNERS_OF_TETRASIDE; j++ )
    {
      j1 = (j+1)%CORNERS_OF_TETRASIDE;
      if (!CalcCrossingPoint(ScreenPoint[0][i],ScreenPoint[0][i1],ScreenPoint[1][j], ScreenPoint[1][j1],&alpha,&beta)) continue;

      /* crossing found, now decide which hides the other */
      z[0] = (1.0-alpha)*Triangle[0][i][2] + alpha*Triangle[0][i1][2];
      z[1] = (1.0-beta)*Triangle[1][j][2] + beta*Triangle[1][j1][2];
      if (ABS(z[0]-z[1])<SMALL_C) continue;
      if (z[0]>z[1])
        return (1);
      else
        return (-1);
    }
  }

  /* now check, if one side contains the other. In that case, the */
  /* covering side also covers the other one's center of mass.    */
  for (i=0; i<2; i++)
  {
    /* now invert a 3x3 system */
    M[0]=(COORD)ScreenPoint[i][0].x, M[3]=(COORD)ScreenPoint[i][1].x, M[6] =(COORD)ScreenPoint[i][2].x;
    M[1]=(COORD)ScreenPoint[i][0].y, M[4]=(COORD)ScreenPoint[i][1].y, M[7] =(COORD)ScreenPoint[i][2].y;
    M[2]=1.0,                                                M[5]=1.0,                                                M[8] =1.0;

    if (M3_Invert(Inverse, M)) continue;
    j=1-i;
    rhs[0] =(COORD) ((ScreenPoint[j][0].x + ScreenPoint[j][1].x + ScreenPoint[j][2].x)/3.0);
    rhs[1] =(COORD) ((ScreenPoint[j][0].y + ScreenPoint[j][1].y + ScreenPoint[j][2].y)/3.0);
    rhs[2] = 1.0;
    M3_TIMES_V3(Inverse,rhs,lambda)

    /* decide if MidPoint of Triangle[j] lies in the interior of Triangle[i] */
    if (lambda[0]>=0.0 && lambda[1]>=0.0 && lambda[2]>=0.0)
    {
      z[i] = lambda[0]*Triangle[i][0][2] + lambda[1]*Triangle[i][1][2] + lambda[2]*Triangle[i][2][2];
      z[j] = (Triangle[j][0][2] + Triangle[j][1][2] + Triangle[j][2][2])/3.0;
      if (z[i]>z[j])
        return ((i>j) ? (-1) : (1));
      else
        return ((i>j) ? (1) : (-1));
    }
  }
  return (0);
}

/****************************************************************************/
/*
   CompareElements - Test, whether a element hides another

   SYNOPSIS:
   static int CompareElements (const void *ElementHandle0, const void *ElementHandle1);

   PARAMETERS:
   .  ElementHandle0 -
   .  ElementHandle1 -

   DESCRIPTION:
   This function tests, whether a element hides another.

   RETURN VALUE:
   INT
   .n     1 whenelement0 hides element1
   .n      -1 when element1 hides element0
   .n     0 when elements do not hide each other.
 */
/****************************************************************************/

static int CompareElements (const void *ElementHandle0, const void *ElementHandle1)
{
  ELEMENT *theElement0, *theElement1;
  INT i, j, k, i1, k1, b0, b1, found, view0, view1, num0, num1;
  COORD radius[2], norm, alpha, beta[2];
  COORD *Corners[2][4];
  COORD_VECTOR VectorMid[2], ViewDir[2], Vector0, Vector1, Triangle[2][3];
  COORD_POINT ScreenPoints[2][3];

  theElement0 = *((ELEMENT **)ElementHandle0);
  theElement1 = *((ELEMENT **)ElementHandle1);

  /* test, if elements have a common side */
  for (i=0; i<SIDES_OF_ELEM(theElement0); i++)
    if( NBELEM(theElement0,i) == theElement1 )
    {
      if (VIEWABLE(theElement0,i))
        return(-1);
      for (j=0; j<SIDES_OF_ELEM(theElement1); j++)
        if( NBELEM(theElement1,j) == theElement0 )
          if (VIEWABLE(theElement1,j))
            return(1);
      return(0);
    }

  /* do some initializing */
  V3_CLEAR(VectorMid[0])
  V3_CLEAR(VectorMid[1])
  for (i=0; i<4; i++)
  {
    Corners[0][i] = CVECT(MYVERTEX(CORNER(theElement0,i)));
    Corners[1][i] = CVECT(MYVERTEX(CORNER(theElement1,i)));
    V3_ADD(Corners[0][i],VectorMid[0],VectorMid[0])
    V3_ADD(Corners[1][i],VectorMid[1],VectorMid[1])
  }
  V3_SCALE(0.25,VectorMid[0])
  V3_SCALE(0.25,VectorMid[1])

  radius[0] = radius[1] = 0.0;
  for (i=0; i<EDGES_OF_ELEM(theElement0); i++)
    for (j=0; j<2; j++)
    {
      V3_SUBTRACT(Corners[j][CornerOfEdge[i][0]],Corners[j][CornerOfEdge[i][1]],Vector0)
      V3_EUKLIDNORM(Vector0,norm)
      radius[j] = MAX(radius[j],norm);
    }

  /* check if tetrahedrons are contained in spheres which do not overlap */
  switch(OBS_Perspective)
  {
  case (YES) :
    for (i=0; i<2; i++)
    {
      V3_SUBTRACT(VectorMid[i],VO_VP(OE_ViewedObj),ViewDir[i])
      V3_EUKLIDNORM(ViewDir[i],norm);
      beta[i] = SQRT(1.0-radius[i]*radius[i]/norm/norm);
      V3_Normalize(ViewDir[i]);
    }
    V3_SCALAR_PRODUCT(ViewDir[0],ViewDir[1],alpha)
    if (beta[0]*beta[1]-SQRT(1.0-beta[0]*beta[0])*SQRT(1.0-beta[1]*beta[1]) > alpha)
      return(0);
    break;
  case (NO) :
    V3_SUBTRACT(VectorMid[0],VectorMid[1],Vector0)
    V3_Orthogonalize(Vector0,OBS_ViewDirection,Vector1);
    V3_EUKLIDNORM(Vector1,norm)
    if (radius[0]+radius[1]<norm)
      return (0);
    break;
  default :
    return (1);
  }

  /* determine the viewable sides and its numbers */
  view0  = VSIDES(theElement0);
  view1  = VSIDES(theElement1);
  num0   = NoOfViewableSides[view0];
  num1   = NoOfViewableSides[view1];

  /* use visible or unvisible sides, depending on which are less */
  b0 = (num0>2);
  b1 = (num1>2);
  if (b0) num0 = 4 - num0;
  if (b1) num1 = 4 - num1;

  /* test the tetrahedrons by testing triangles */
  i1=0;
  for (i=0; i<num0; i++)
  {
    /* determine triangle of theElement0 */
    while( ((view0>>i1)&1) == b0) i1++;
    for (j=0; j<CORNERS_OF_TETRASIDE; j++)
    {
      V3_TRAFOM4_V3(Corners[0][CornerOfSide[i1][j]],ObsTrafo,Triangle[0][j])
      (*OBS_ProjectProc)(Triangle[0][j],&(ScreenPoints[0][j]));
    }

    /* determine triangle of theElement1 and compare triangles */
    k1=0;
    for (k=0; k<num1; k++)
    {
      while ( ((view1>>k1)&1) == b1 ) k1++;
      for (j=0; j<CORNERS_OF_TETRASIDE; j++)
      {
        V3_TRAFOM4_V3(Corners[1][CornerOfSide[k1][j]],ObsTrafo,Triangle[1][j])
        (*OBS_ProjectProc)(Triangle[1][j],&(ScreenPoints[1][j]));
      }
      found = CompareTriangles(Triangle,ScreenPoints);
      if (found)
        return (found);
      k1++;
    }
    i1++;
  }

  return(0);
}

/****************************************************************************/
/*
   OrderElements - order elements w.r.t. theViewedObject

   SYNOPSIS:
   static INT OrderElements_3D (MULTIGRID *theMG, VIEWEDOBJ *theViewedObj);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theViewObject -

   DESCRIPTION:
   This function orders elements w.r.t. theViewedObject.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT OrderElements_3D (MULTIGRID *theMG, VIEWEDOBJ *theViewedObj)
{
  HEAP *theHeap;
  ELEMENT **table, *theElement;
  GRID *theGrid;
  INT i;

  /* check if multigrid is allready ordered */
  /*.....*/

  /* inits */
  OE_ViewedObj = theViewedObj;

  /* calculate the viewable sides of all elements on all levels */
  for (i=0; i<=theMG->topLevel; i++)
    CalcViewableSidesOnGrid(theMG->grids[i]);

  /* allocate memory for the element list on level 0 (but at least for maximal number of sons of an element) */
  theGrid = theMG->grids[0];
  if (theGrid->nElem<2)
  {
    UserWrite("elements need not to be ordered\n");
    return(0);
  }
  theHeap = theMG->theHeap;
  Mark(theHeap,FROM_TOP);
  if ( (table=(ELEMENT **)GetMem(theHeap,MAX(theGrid->nElem,MAX_SONS)*sizeof(ELEMENT *),FROM_TOP)) == NULL )
  {
    Release(theHeap,FROM_TOP);
    UserWrite("ERROR: could not allocate memory from the MGHeap\n");
    return (1);
  }

  /* order elements on level zero */
  i=0;
  for (theElement=theGrid->elements; theElement!= NULL; theElement=SUCCE(theElement))
    table[i++] = theElement;
  if (i!=NT(theGrid)) return (1);
  SelectionSort((void *)table,theGrid->nElem,sizeof(*table),CompareElements);
  if (PutAtStartOfList(theGrid,i,table)!=GM_OK) return (1);

  /* now order level 1 to toplevel hirarchically */
  for (i=0; i<theMG->topLevel; i++)
  {
    theGrid = theMG->grids[i];
    for (theElement=LASTELEMENT(theGrid); theElement!= NULL; theElement=PREDE(theElement))
    {
      if (NSONS(theElement)<=0) continue;
      OrderSons(table,theElement);
      if (PutAtStartOfList(UPGRID(theGrid),NSONS(theElement),table)!=GM_OK) return (1);
    }
  }

  /* release heap */
  Release(theHeap,FROM_TOP);

  /* store the view to which elements are ordered */
  /*...*/

  return (0);
}

/****************************************************************************/
/*
   OrderNodes - order nodes in elements

   SYNOPSIS:
   static INT OrderNodes (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - pointer to multigrid

   DESCRIPTION:
   This function orders nodes in elements with respect to cutting orientation on all levels.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT OrderNodes (MULTIGRID *theMG)
{
  ELEMENT *theElement;
  COORD z[MAX_CORNERS_OF_ELEM];
  INT i, j;

  /* check if nodes are allready ordered w.r.t current cut */
  /*...*/

  /* calculate the node order of all elements on all levels */
  for (j=0; j<=theMG->topLevel; j++)
    for (theElement=FIRSTELEMENT(theMG->grids[j]); theElement!= NULL; theElement=SUCCE(theElement))
    {
      /* calculate Z-coordinates in the cutting VRS of corners */
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
        V3_TRAFO4_SC(CVECT(MYVERTEX(CORNER(theElement,i))),CutTrafo,z[i])

        i = ((z[0]>=z[1])<<5) + ((z[0]>=z[2])<<4) + ((z[0]>=z[3])<<3)
            + ((z[1]>=z[2])<<2) + ((z[1]>=z[3])<<1) +  (z[2]>=z[3]);

      if (OrderIndex[i] == -1)
        return(1);

      /* set flags in element */
      SETNORDER(theElement,OrderIndex[i]);
      if (z[CornerIndex[OrderIndex[i]][0]] < -SMALL_C)
        SETCUTMODE(theElement,CM_BEHIND);
      else if (z[CornerIndex[OrderIndex[i]][3]] > SMALL_C)
        SETCUTMODE(theElement,CM_INFRONT);
      else
        SETCUTMODE(theElement,CM_INTERSECT);
    }

  /* store the theCut to which nodes are ordered */
  /*...*/

  return(0);
}

/****************************************************************************/
/*
   EW_PreProcess_PlotGrid3D - Initialize input variables of EW_ElementEval3D for GridPlot3D

   SYNOPSIS:
   static INT EW_PreProcess_PlotGrid3D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initialize input variables of EW_ElementEval3D for GridPlot3D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_PlotGrid3D (PICTURE *thePicture, WORK *theWork)
{
  struct GridPlotObj3D *theGpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theGpo = &(PIC_PO(thePicture)->theGpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  EE3D_NoColor[COLOR_LOWER_LEVEL]         = 1;
  EE3D_NoColor[COLOR_CUT_EDGE]            = 0;
  if (theGpo->ElemColored == YES)
  {
    EE3D_NoColor[COLOR_COPY]                = 0;
    EE3D_NoColor[COLOR_IRR]                 = 0;
    EE3D_NoColor[COLOR_REG]                 = 0;
    EE3D_NoColor[COLOR_REG]                 = 0;
  }
  else
  {
    EE3D_NoColor[COLOR_COPY]                = 1;
    EE3D_NoColor[COLOR_IRR]                 = 1;
    EE3D_NoColor[COLOR_REG]                 = 1;
  }

  EE3D_Color[COLOR_COPY]                  = theOD->yellow;
  EE3D_Color[COLOR_IRR]                   = theOD->green;
  EE3D_Color[COLOR_REG]                   = theOD->red;
  EE3D_Color[COLOR_LOWER_LEVEL]   = theOD->white;
  EE3D_Color[COLOR_EDGE]                  = theOD->black;
  EE3D_Color[COLOR_CUT_EDGE]              = theOD->orange;

  EE3D_Elem2Plot[PLOT_ALL]                = 0;
  EE3D_Elem2Plot[PLOT_COPY]               = 0;
  EE3D_Elem2Plot[PLOT_IRR]                = 0;
  EE3D_Elem2Plot[PLOT_REG]                = 0;


  switch (theGpo->WhichElem)
  {
  case PO_NO :
    break;
  case PO_ALL :
    EE3D_Elem2Plot[PLOT_ALL] = 1;
  case PO_COPY :
    EE3D_Elem2Plot[PLOT_COPY] = 1;
  case PO_IRR :
    EE3D_Elem2Plot[PLOT_IRR] = 1;
  case PO_REG :
    EE3D_Elem2Plot[PLOT_REG] = 1;
  }

  /* build cut trafo */
  if (BuildCutTrafo(&(theGpo->theCut),OBS_ViewDirection)) return (1);

  /* order nodes if */
  if (theGpo->theCut.status==ACTIVE)
    if (OrderNodes(theMG)) return (1);

  /* mark surface elements */
  EE3D_MaxLevel = CURRENTLEVEL(theMG);
  if (MarkElements3D(theMG,0,CURRENTLEVEL(theMG))) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_SelectElement3D - Initialize input variables of EW_ElementEval3D for GridPlot3D

   SYNOPSIS:
   static INT EW_PreProcess_SelectElement3D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes input variables of EW_ElementEval3D for GridPlot3D.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if an error occured.
 */
/****************************************************************************/

static INT EW_PreProcess_SelectElement3D (PICTURE *thePicture, WORK *theWork)
{
  FE3D_MousePos.x = W_SELECTELEMENT_WORK(theWork)->PixelX;
  FE3D_MousePos.y = W_SELECTELEMENT_WORK(theWork)->PixelY;

  if (EW_PreProcess_PlotGrid3D (thePicture,theWork)) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EScalar3D_BackGrid - Initialize for plot of backgrid for
   C(olor)C(ontour) plot

   SYNOPSIS:
   static INT EW_PreProcess_EScalar3D_BackGrid (PICTURE *thePicture,
   WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for plot of backgrid for C(olor)C(ontour) plot.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if an error occured
 */
/****************************************************************************/

static INT EW_PreProcess_EScalar3D_BackGrid (PICTURE *thePicture, WORK *theWork)
{
  struct ElemScalarPlotObj3D *theEspo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theEspo = &(PIC_PO(thePicture)->theEspo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  EE3D_NoColor[COLOR_CUT_EDGE]            = 1;
  EE3D_NoColor[COLOR_LOWER_LEVEL]         = 1;
  EE3D_NoColor[COLOR_COPY]                        = 1;
  EE3D_NoColor[COLOR_IRR]                         = 1;
  EE3D_NoColor[COLOR_REG]                         = 1;

  EE3D_Color[COLOR_COPY]                  = theOD->yellow;
  EE3D_Color[COLOR_IRR]                   = theOD->green;
  EE3D_Color[COLOR_REG]                   = theOD->red;
  EE3D_Color[COLOR_LOWER_LEVEL]   = theOD->white;
  EE3D_Color[COLOR_EDGE]                  = theOD->black;
  EE3D_Color[COLOR_CUT_EDGE]              = theOD->orange;

  EE3D_Elem2Plot[PLOT_ALL]                = 1;
  EE3D_Elem2Plot[PLOT_COPY]               = 1;
  EE3D_Elem2Plot[PLOT_IRR]                = 1;
  EE3D_Elem2Plot[PLOT_REG]                = 1;

  /* build cut trafo */
  if (BuildCutTrafo(&(theEspo->theCut),OBS_ViewDirection)) return (1);

  /* order nodes if */
  if (theEspo->theCut.status==ACTIVE)
    if (OrderNodes(theMG)) return (1);

  /* mark suface elements on boundary and cut if */
  if (theEspo->mode == PO_COLOR)
  {
    if (MarkElements_MGS_Bnd(theMG,0,CURRENTLEVEL(theMG))) return (1);
  }
  else
  if (MarkElements_MGS_Bnd_Cut(theMG,0,CURRENTLEVEL(theMG))) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_CutBnd3D - Initialize for plot of backgrid for C(olor)C(ontour) plot

   SYNOPSIS:
   static INT EW_PreProcess_CutBnd3D (PICTURE *thePicture,
   WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for plot of backgrid for C(olor)C(ontour) plot.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if an error occured
 */
/****************************************************************************/

static INT EW_PreProcess_CutBnd3D (PICTURE *thePicture, WORK *theWork)
{
  struct ElemScalarPlotObj3D *theEspo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theEspo = &(PIC_PO(thePicture)->theEspo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  if (theEspo->theCut.status!=ACTIVE)
    return (1);

  EE3D_Color[COLOR_CUT_EDGE]              = theOD->black;

  /* build cut trafo */
  if (BuildCutTrafo(&(theEspo->theCut),OBS_ViewDirection)) return (1);

  /* order nodes */
  if (OrderNodes(theMG)) return (1);

  /* mark suface elements on boundary which are cut */
  if (MarkElements_MGS_Bnd_and_Cut(theMG,0,CURRENTLEVEL(theMG))) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EVector3D_BackGrid - Initialize for plot of backgrid for
   C(olor)C(ontour) plot

   SYNOPSIS:
   static INT EW_PreProcess_EVector3D_BackGrid (PICTURE *thePicture,
   WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for plot of backgrid for C(olor)C(ontour) plot.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if an error occured
 */
/****************************************************************************/

static INT EW_PreProcess_EVector3D_BackGrid (PICTURE *thePicture, WORK *theWork)
{
  struct ElemVectorPlotObj3D *theEvpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theEvpo = &(PIC_PO(thePicture)->theEvpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  EE3D_NoColor[COLOR_CUT_EDGE]            = 1;
  EE3D_NoColor[COLOR_LOWER_LEVEL]         = 1;
  EE3D_NoColor[COLOR_COPY]                        = 1;
  EE3D_NoColor[COLOR_IRR]                         = 1;
  EE3D_NoColor[COLOR_REG]                         = 1;

  EE3D_Color[COLOR_COPY]                  = theOD->yellow;
  EE3D_Color[COLOR_IRR]                   = theOD->green;
  EE3D_Color[COLOR_REG]                   = theOD->red;
  EE3D_Color[COLOR_LOWER_LEVEL]   = theOD->white;
  EE3D_Color[COLOR_EDGE]                  = theOD->black;
  EE3D_Color[COLOR_CUT_EDGE]              = theOD->orange;

  EE3D_Elem2Plot[PLOT_ALL]                = 1;
  EE3D_Elem2Plot[PLOT_COPY]               = 1;
  EE3D_Elem2Plot[PLOT_IRR]                = 1;
  EE3D_Elem2Plot[PLOT_REG]                = 1;

  /* build cut trafo */
  if (BuildCutTrafo(&(theEvpo->theCut),OBS_ViewDirection)) return (1);

  /* order nodes if */
  if (theEvpo->theCut.status==ACTIVE)
    if (OrderNodes(theMG)) return (1);

  /* mark suface elements on boundary and cut */
  if (MarkElements_MGS_Bnd_Cut(theMG,0,CURRENTLEVEL(theMG))) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EScalar3D - Initialize for C(olor)C(ontour) plot

   SYNOPSIS:
   static INT EW_PreProcess_EScalar3D (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for C(olor)C(ontour) plot.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT EW_PreProcess_EScalar3D (PICTURE *thePicture, WORK *theWork)
{
  struct ElemScalarPlotObj3D *theEspo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;
  INT i;

  theEspo = &(PIC_PO(thePicture)->theEspo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  /* set value->color fct, eval fct */
  if (theEspo->max - theEspo->min < SMALL_D)
    if (W_ID(theWork) != FINDRANGE_WORK)
    {
      UserWrite("maxValue has to be larger than minValue\n");
      return (1);
    }

  /* do not plot if cut plane is on the back */
  if (!CUT_CutAtFront) return (1);

  EScalar3D_EvalFct        = theEspo->EvalFct->PlotProc;
  EScalar3D_V2C_factor = (theOD->spectrumEnd - theOD->spectrumStart)/(theEspo->max - theEspo->min);
  EScalar3D_V2C_offset = theOD->spectrumStart - EScalar3D_V2C_factor*theEspo->min;
  EScalar3D_mode           = theEspo->mode;
  if (EScalar3D_mode == PO_CONTOURS_EQ)
  {
    EScalar3D_numOfContours = theEspo->numOfContours;
    for (i=0; i<theEspo->numOfContours; i++)
      EScalar3D_ContValues[i] = theEspo->contValues[i];
  }
  EScalar3D_depth                 = theEspo->depth;

  /* mark suface elements on boundary */
  if (MarkElements_MGS_Cut(theMG,0,CURRENTLEVEL(theMG))) return (1);

  /* prepare evaluation routine */
  if (theEspo->EvalFct->PreprocessProc != NULL)
    if ((*theEspo->EvalFct->PreprocessProc)(theMG)) return (1);

  return (0);
}

/****************************************************************************/
/*
   EW_PreProcess_EScalar3D_FR - Initialize for findrange of scalar plot

   SYNOPSIS:
   static INT EW_PreProcess_EScalar3D_FR (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture -
   .  theWork -

   DESCRIPTION:
   This function initializes for findrange of scalar plot.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT EW_PreProcess_EScalar3D_FR (PICTURE *thePicture, WORK *theWork)
{
  struct ElemScalarPlotObj3D *theEspo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theEspo = &(PIC_PO(thePicture)->theEspo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  if (theEspo->theCut.status!=ACTIVE) return (1);

  EScalar3D_EvalFct        = theEspo->EvalFct->PlotProc;
  EScalar3D_V2C_factor = (theOD->spectrumEnd - theOD->spectrumStart);
  EScalar3D_V2C_offset = theOD->spectrumStart;
  EScalar3D_mode           = PO_COLOR;
  EScalar3D_depth                 = theEspo->depth;

  /* build cut trafo */
  if (BuildCutTrafo(&(theEspo->theCut),OBS_ViewDirection)) return (1);

  /* order nodes */
  if (OrderNodes(theMG)) return (1);

  /* mark suface elements on boundary */
  if (MarkElements_MGS_Cut(theMG,0,CURRENTLEVEL(theMG))) return (1);

  /* prepare evaluation routine */
  if (theEspo->EvalFct->PreprocessProc != NULL)
    if ((*theEspo->EvalFct->PreprocessProc)(theMG)) return (1);

  /* reset min and max values */
  GEN_FR_put = W_FINDRANGE_WORK(theWork)->put;
  GEN_FR_min = MAX_D;
  GEN_FR_max = -MAX_D;

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  EW_PreProcess_EVector3D										*/
/*																			*/
/* Purpose:   initialize for vector plot 3D                                                             */
/*																			*/
/* Input:	  PICTURE *thePicture, WORK *theWork							*/
/*																			*/
/* Return:	  INT 0: ok                                                                                                     */
/*			  INT 1: an error occurred										*/
/*																			*/
/****************************************************************************/

static INT EW_PreProcess_EVector3D (PICTURE *thePicture, WORK *theWork)
{
  struct ElemVectorPlotObj3D *theEvpo;
  OUTPUTDEVICE *theOD;
  MULTIGRID *theMG;

  theEvpo = &(PIC_PO(thePicture)->theEvpo);
  theOD  = PIC_OUTPUTDEV(thePicture);
  theMG  = PO_MG(PIC_PO(thePicture));

  /* set value->length fct, eval fct */
  if (theEvpo->max < SMALL_D)
    if (W_ID(theWork) != FINDRANGE_WORK)
    {
      UserWrite("maxValue has to be larger than zero\n");
      return (1);
    }

  /* do not plot if cut plane is on the back */
  if (!CUT_CutAtFront) return (1);

  EVector_rastersize              = theEvpo->RasterSize;
  EVector_cutvector               = theEvpo->CutVector;
  EVector_CutLenFactor    = theEvpo->CutLenFactor;
  EVector3D_projectvector = theEvpo->ProjectVector;
  EVector_EvalFct                 = theEvpo->EvalFct->PlotProc;
  EVector_V2L_factor              = EVector_rastersize/theEvpo->max;                                    /* scale length of vectors			*/
  EVector3D_V2C_factor    = 0.5*(theOD->spectrumEnd - theOD->spectrumStart);       /* transformation from (-1,1) to     */
  EVector3D_V2C_offset    = theOD->spectrumStart + EVector3D_V2C_factor;                /* color spectrum					*/
  EVector_ColorCut                = theOD->black;

  /* mark suface elements on boundary */
  if (MarkElements_MGS_Cut(theMG,0,CURRENTLEVEL(theMG))) return (1);

  /* prepare evaluation routine */
  if (theEvpo->EvalFct->PreprocessProc!=NULL)
    if ((*theEvpo->EvalFct->PreprocessProc)(theMG)) return (1);;

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  EW_PreProcess_EVector3D_FR									*/
/*																			*/
/* Purpose:   initialize for findrange of vector plot						*/
/*																			*/
/* Input:	  PICTURE *thePicture, WORK *theWork							*/
/*																			*/
/* Return:	  INT 0: ok                                                                                                     */
/*			  INT 1: an error occurred										*/
/*																			*/
/****************************************************************************/

static INT EW_PreProcess_EVector3D_FR (PICTURE *thePicture, WORK *theWork)
{
  if (EW_PreProcess_EVector3D(thePicture,theWork))
    return (1);

  /* reset min and max values */
  GEN_FR_put = W_FINDRANGE_WORK(theWork)->put;
  GEN_FR_min = MAX_D;
  GEN_FR_max = -MAX_D;

  return (0);
}

/****************************************************************************/
/*
   PlotColorTriangle3D - Plot on triangle color(3D coord) with depth

   SYNOPSIS:
   static INT PlotColorTriangle3D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *TP0, COORD *TP1, COORD *TP2,
   INT depth, DRAWINGOBJ **theDO);

   PARAMETERS:
   .  theElement -
   .  CornersOfElem -
   .  TP0 -
   .  TP1 -
   .  TP2 -
   .  depth -
   .  theDO - the drawing object to draw on

   DESCRIPTION:
   This function plot on triangle color(3D coord) with depth.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT PlotColorTriangle3D (ELEMENT *theElement, COORD **CornersOfElem, COORD *TP0, COORD *TP1, COORD *TP2, INT depth, DRAWINGOBJ **theDO)
{
  COORD_VECTOR EvalPoint, LocalCoord, MP0, MP1, MP2;
  INT i;
  long Color;
  DOUBLE value;

  if (depth<=0)
  {
    /* get values */
    for (i=0; i<DIM; i++)
      EvalPoint[i] = (TP0[i]+TP1[i]+TP2[i])/3.0;
    if (GlobalToLocal3d(CornersOfElem,EvalPoint,LocalCoord)) return (1);
    value = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    Color = (long)(EScalar3D_V2C_factor*value+EScalar3D_V2C_offset);
    Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
    Color = MAX(Color,WOP_OutputDevice->spectrumStart);

    /* draw */
    DO_2c(*theDO) = DO_POLYGON; DO_inc(*theDO)
    DO_2c(*theDO) = 3; DO_inc(*theDO)
    DO_2l(*theDO) = Color; DO_inc(*theDO)
    V3_COPY(TP0,DO_2Cp(*theDO)); DO_inc_n(*theDO,3)
    V3_COPY(TP1,DO_2Cp(*theDO)); DO_inc_n(*theDO,3)
    V3_COPY(TP2,DO_2Cp(*theDO)); DO_inc_n(*theDO,3)

    /* store range */
    EScalar3D_minValue = MIN(EScalar3D_minValue,value);
    EScalar3D_maxValue = MAX(EScalar3D_maxValue,value);
  }
  else
  {
    /* find corners of subdivided triangles */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(TP0[i]+TP1[i]);
      MP1[i] = 0.5*(TP1[i]+TP2[i]);
      MP2[i] = 0.5*(TP2[i]+TP0[i]);
    }
    if (PlotColorTriangle3D(theElement,CornersOfElem,TP0,MP0,MP2,depth-1,theDO)) return (1);
    if (PlotColorTriangle3D(theElement,CornersOfElem,MP0,TP1,MP1,depth-1,theDO)) return (1);
    if (PlotColorTriangle3D(theElement,CornersOfElem,TP2,MP2,MP1,depth-1,theDO)) return (1);
    if (PlotColorTriangle3D(theElement,CornersOfElem,MP0,MP1,MP2,depth-1,theDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   PlotColorQuadrilateral3D - Plot on quadrilateral color(3D coord) with depth

   SYNOPSIS:
   static INT PlotColorQuadrilateral3D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *QP0, COORD *QP1, COORD *QP2,
   COORD *QP3, INT depth, DRAWINGOBJ **theDO);

   PARAMETERS:
   .  theElement -
   .  CornersOfElem -
   .  QP0 -
   .  QP1 -
   .  QP2 -
   .  QP3 -
   .  depth -
   .  theDO - the drawing object to draw on

   DESCRIPTION:
   This function plots on quadrilateral color(3D coord) with depth.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT PlotColorQuadrilateral3D (ELEMENT *theElement, COORD **CornersOfElem, COORD *QP0, COORD *QP1, COORD *QP2, COORD *QP3, INT depth, DRAWINGOBJ **theDO)
{
  COORD_VECTOR EVP, LocalCoord, MP0, MP1, MP2, MP3;
  INT i;
  long Color;
  DOUBLE value;

  for (i=0; i<DIM; i++)
    EVP[i] = (QP0[i]+QP1[i]+QP2[i]+QP3[i])*0.25;
  if (depth<=0)
  {
    /* get values */
    if (GlobalToLocal3d(CornersOfElem,EVP,LocalCoord)) return (1);
    value = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    Color = (long)(EScalar3D_V2C_factor*value+EScalar3D_V2C_offset);
    Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
    Color = MAX(Color,WOP_OutputDevice->spectrumStart);

    /* draw */
    DO_2c(*theDO) = DO_POLYGON; DO_inc(*theDO)
    DO_2c(*theDO) = 4; DO_inc(*theDO)
    DO_2l(*theDO) = Color; DO_inc(*theDO)
    V3_COPY(QP0,DO_2Cp(*theDO)); DO_inc_n(*theDO,3)
    V3_COPY(QP1,DO_2Cp(*theDO)); DO_inc_n(*theDO,3)
    V3_COPY(QP2,DO_2Cp(*theDO)); DO_inc_n(*theDO,3)
    V3_COPY(QP3,DO_2Cp(*theDO)); DO_inc_n(*theDO,3)

    /* store range */
    EScalar3D_minValue = MIN(EScalar3D_minValue,value);
    EScalar3D_maxValue = MAX(EScalar3D_maxValue,value);
  }
  else
  {
    /* find corners of subdivided quadrilaterals */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(QP0[i]+QP1[i]);
      MP1[i] = 0.5*(QP1[i]+QP2[i]);
      MP2[i] = 0.5*(QP2[i]+QP3[i]);
      MP3[i] = 0.5*(QP3[i]+QP0[i]);
    }
    if (PlotColorQuadrilateral3D(theElement,CornersOfElem,QP0,MP0,EVP,MP3,depth-1,theDO)) return (1);
    if (PlotColorQuadrilateral3D(theElement,CornersOfElem,MP0,QP1,MP1,EVP,depth-1,theDO)) return (1);
    if (PlotColorQuadrilateral3D(theElement,CornersOfElem,EVP,MP1,QP2,MP2,depth-1,theDO)) return (1);
    if (PlotColorQuadrilateral3D(theElement,CornersOfElem,EVP,MP2,QP3,MP3,depth-1,theDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   PointOnLine3D - Cals point between two points with contourValue

   SYNOPSIS:
   static INT PointOnLine3D (DOUBLE contourValue, DOUBLE value0,
   DOUBLE value1, COORD_VECTOR vec0, COORD_VECTOR vec1, COORD_VECTOR p);

   PARAMETERS:
   .  contourValue -
   .  value0 -
   .  value1 -
   .  vec0 -
   .  vec1 -
   .  p -

   DESCRIPTION:
   This function cals point between two points with contourValue.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT PointOnLine3D (DOUBLE contourValue, DOUBLE value0, DOUBLE value1, COORD_VECTOR vec0, COORD_VECTOR vec1, COORD_VECTOR p)
{
  DOUBLE alpha;

  if (ABS(value0-value1)<SMALL_D)
  {
    if (value0==contourValue)
    {
      V3_LINCOMB(0.5,vec0,0.5,vec1,p);
      return (1);
    }
    return (0);
  }
  else
  {
    alpha = (contourValue-value0)/(value1-value0);
    if (alpha<0.0 || alpha>=1.0) return (0);
    V3_LINCOMB(1.0-alpha,vec0,alpha,vec1,p);
    return (1);
  }
}

/****************************************************************************/
/*
   PlotContourTriangle3D - Plot on triangle contourlines (3D coord) with depth

   SYNOPSIS:
   static INT PlotContourTriangle3D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *TP0, COORD *TP1, COORD *TP2,
   INT depth, DRAWINGOBJ **theDO);

   PARAMETERS:
   .  theElement -
   .  CornerOfElem -
   .  TP) -
   .  TP1 -
   .  TP2 -
   .  depth -
   .  theDO - the drawing object to draw on

   DESCRIPTION:
   This function plots on triangle contourlines (3D coord) with depth.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT PlotContourTriangle3D (ELEMENT *theElement, COORD **CornersOfElem, COORD *TP0, COORD *TP1, COORD *TP2, INT depth, DRAWINGOBJ **theDO)
{
  COORD_VECTOR LocalCoord, MP0, MP1, MP2, PointMid, Point[3];
  INT i, j, n, min, max;
  long Color;
  DOUBLE v0, v1, v2, vmin, vmax;

  if (depth<=0)
  {
    /* get values at the corners */
    if (GlobalToLocal3d(CornersOfElem,TP0,LocalCoord)) return (1);
    v0      = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal3d(CornersOfElem,TP1,LocalCoord)) return (1);
    v1      = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal3d(CornersOfElem,TP2,LocalCoord)) return (1);
    v2      = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    vmin = MIN(v0,v1); vmin = MIN(vmin,v2);
    vmax = MAX(v0,v1); vmax = MAX(vmax,v2);

    /* store range */
    EScalar3D_minValue = MIN(EScalar3D_minValue,vmin);
    EScalar3D_maxValue = MAX(EScalar3D_maxValue,vmax);

    /* find contours to be plotted */
    for (min=0; min<EScalar3D_numOfContours; min++)
      if (EScalar3D_ContValues[min]>=vmin)
        break;
    for (max=EScalar3D_numOfContours-1; max>=0; max--)
      if (EScalar3D_ContValues[max]<=vmax)
        break;

    /* draw contours */
    for (i=min; i<=max; i++)
    {
      /* set color */
      Color = (long)(EScalar3D_V2C_factor*EScalar3D_ContValues[i]+EScalar3D_V2C_offset);
      Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
      Color = MAX(Color,WOP_OutputDevice->spectrumStart);

      /* calculate points on each side of triangle having the right value */
      n=0;
      if (PointOnLine3D(EScalar3D_ContValues[i],v0,v1,TP0,TP1,Point[n])) n++;
      if (PointOnLine3D(EScalar3D_ContValues[i],v1,v2,TP1,TP2,Point[n])) n++;
      if (PointOnLine3D(EScalar3D_ContValues[i],v2,v0,TP2,TP0,Point[n])) n++;

      /* draw */
      switch (n)
      {
      case 2 :
        DO_2c(*theDO) = DO_LINE; DO_inc(*theDO)
        DO_2l(*theDO) = Color; DO_inc(*theDO)
        V3_COPY(Point[0],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[1],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        break;
      case 3 :
        DO_2c(*theDO) = DO_POLYLINE; DO_inc(*theDO)
        DO_2c(*theDO) = 5; DO_inc(*theDO)
        DO_2l(*theDO) = Color; DO_inc(*theDO);
        for (j=0; j<3; j++)
          PointMid[j] = (Point[0][j]+Point[1][j]+Point[2][j])/3.0;
        V3_COPY(Point[0],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(PointMid,DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[1],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(PointMid,DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[2],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        break;
      }
    }
  }
  else
  {
    /* find corners of subdivided triangles */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(TP0[i]+TP1[i]);
      MP1[i] = 0.5*(TP1[i]+TP2[i]);
      MP2[i] = 0.5*(TP2[i]+TP0[i]);
    }
    if (PlotContourTriangle3D(theElement,CornersOfElem,TP0,MP0,MP2,depth-1,theDO)) return (1);
    if (PlotContourTriangle3D(theElement,CornersOfElem,MP0,TP1,MP1,depth-1,theDO)) return (1);
    if (PlotContourTriangle3D(theElement,CornersOfElem,TP2,MP2,MP1,depth-1,theDO)) return (1);
    if (PlotContourTriangle3D(theElement,CornersOfElem,MP0,MP1,MP2,depth-1,theDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   PlotContourQuadrilateral3D - Plot on quadrilateral contourlines (3D coord) with depth

   SYNOPSIS:
   static INT PlotContourQuadrilateral3D (ELEMENT *theElement,
   COORD **CornersOfElem, COORD *QP0, COORD *QP1, COORD *QP2,
   COORD *QP3, INT depth, DRAWINGOBJ **theDO);

   PARAMETERS:
   .  theElement -
   .  CornersOfElem -
   .  QP1 -
   .  QP2 -
   .  QP3 -
   .  depth -
   .  theDO - the drawing object to draw on

   DESCRIPTION:
   This function plots on quadrilateral contourlines (3D coord) with depth.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT PlotContourQuadrilateral3D (ELEMENT *theElement, COORD **CornersOfElem, COORD *QP0, COORD *QP1, COORD *QP2, COORD *QP3, INT depth, DRAWINGOBJ **theDO)
{
  COORD_VECTOR EVP, LocalCoord, MP0, MP1, MP2, MP3, PointMid, Point[4];
  INT i, j, n, min, max;
  long Color;
  DOUBLE v0, v1, v2, v3, vmin, vmax;

  for (i=0; i<DIM; i++)
    EVP[i] = (QP0[i]+QP1[i]+QP2[i]+QP3[i])*0.25;
  if (depth<=0)
  {
    /* get values at the corners */
    if (GlobalToLocal3d(CornersOfElem,QP0,LocalCoord)) return (1);
    v0      = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal3d(CornersOfElem,QP1,LocalCoord)) return (1);
    v1      = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal3d(CornersOfElem,QP2,LocalCoord)) return (1);
    v2      = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    if (GlobalToLocal3d(CornersOfElem,QP3,LocalCoord)) return (1);
    v3      = (*EScalar3D_EvalFct)(theElement,CornersOfElem,LocalCoord);
    vmin = MIN(v0,v1); vmin = MIN(vmin,v2); vmin = MIN(vmin,v3);
    vmax = MAX(v0,v1); vmax = MAX(vmax,v2); vmax = MAX(vmin,v3);

    /* store range */
    EScalar3D_minValue = MIN(EScalar3D_minValue,vmin);
    EScalar3D_maxValue = MAX(EScalar3D_maxValue,vmax);

    /* find contours to be plotted */
    for (min=0; min<EScalar3D_numOfContours; min++)
      if (EScalar3D_ContValues[min]>=vmin)
        break;
    for (max=EScalar3D_numOfContours-1; max>=0; max--)
      if (EScalar3D_ContValues[max]<=vmax)
        break;

    /* draw contours */
    for (i=min; i<=max; i++)
    {
      /* set color */
      Color = (long)(EScalar3D_V2C_factor*EScalar3D_ContValues[i]+EScalar3D_V2C_offset);
      Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
      Color = MAX(Color,WOP_OutputDevice->spectrumStart);

      /* calculate points on each side of triangle having the right value */
      n=0;
      if (PointOnLine3D(EScalar3D_ContValues[i],v0,v1,QP0,QP1,Point[n])) n++;
      if (PointOnLine3D(EScalar3D_ContValues[i],v1,v2,QP1,QP2,Point[n])) n++;
      if (PointOnLine3D(EScalar3D_ContValues[i],v2,v3,QP2,QP3,Point[n])) n++;
      if (PointOnLine3D(EScalar3D_ContValues[i],v3,v0,QP3,QP0,Point[n])) n++;

      /* draw */
      switch (n)
      {
      case 2 :
        DO_2c(*theDO) = DO_LINE; DO_inc(*theDO)
        DO_2l(*theDO) = Color; DO_inc(*theDO);
        V3_COPY(Point[0],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[1],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        break;
      case 3 :
        DO_2c(*theDO) = DO_POLYLINE; DO_inc(*theDO)
        DO_2c(*theDO) = 5; DO_inc(*theDO)
        DO_2l(*theDO) = Color; DO_inc(*theDO);
        for (j=0; j<3; j++)
          PointMid[j] = (Point[0][j]+Point[1][j]+Point[2][j])/3.0;
        V3_COPY(Point[0],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(PointMid,DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[1],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(PointMid,DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[2],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        break;
      case 4 :
        DO_2c(*theDO) = DO_LINE; DO_inc(*theDO)
        DO_2l(*theDO) = Color; DO_inc(*theDO);
        V3_COPY(Point[0],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[2],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);

        DO_2c(*theDO) = DO_LINE; DO_inc(*theDO)
        DO_2l(*theDO) = Color; DO_inc(*theDO);
        V3_COPY(Point[1],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
        V3_COPY(Point[3],DO_2Cp(*theDO)); DO_inc_n(*theDO,3);
      }
    }
  }
  else
  {
    /* find corners of subdivided quadrilaterals */
    for (i=0; i<DIM; i++)
    {
      MP0[i] = 0.5*(QP0[i]+QP1[i]);
      MP1[i] = 0.5*(QP1[i]+QP2[i]);
      MP2[i] = 0.5*(QP2[i]+QP3[i]);
      MP3[i] = 0.5*(QP3[i]+QP0[i]);
    }
    if (PlotContourQuadrilateral3D(theElement,CornersOfElem,QP0,MP0,EVP,MP3,depth-1,theDO)) return (1);
    if (PlotContourQuadrilateral3D(theElement,CornersOfElem,MP0,QP1,MP1,EVP,depth-1,theDO)) return (1);
    if (PlotContourQuadrilateral3D(theElement,CornersOfElem,EVP,MP1,QP2,MP2,depth-1,theDO)) return (1);
    if (PlotContourQuadrilateral3D(theElement,CornersOfElem,EVP,MP2,QP3,MP3,depth-1,theDO)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   EW_EScalar3D	- Initialize for C(olor)C(ontour) plot

   SYNOPSIS:
   static INT EW_EScalar3D (ELEMENT *theElement, DRAWINGOBJ *theDO);

   PARAMETERS:
   .  theElement -
   .  theDO -

   DESCRIPTION:
   This function initializes for C(olor)C(ontour) plot.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/

static INT EW_EScalar3D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, n, NodeOrder;
  COORD_VECTOR Poly[MAX_POINTS_OF_POLY];
  COORD *x[MAX_CORNERS_OF_ELEM], z[MAX_CORNERS_OF_ELEM];
  DRAWINGOBJ *range;

  DO_2c(theDO) = DO_NO_INST;

  /* get node order */
  NodeOrder = NORDER(theElement);

  /* get coordinates of corners of the element and their z coordinates in cut system */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));
    V3_TRAFO4_SC(x[i],CutTrafo,z[i])
  }

  /* determine polygon being intersection of element wth cut plane */
  if (GetPolyElemISCutPlane(x,z,NodeOrder,Poly,&n))
    return (1);
  if (n<=2) return (0);

  /* draw polygon with depth */
  EScalar3D_minValue = MAX_D; EScalar3D_maxValue = -MAX_D;
  DO_2c(theDO) = DO_RANGE; DO_inc(theDO); range = theDO; DO_inc_n(theDO,2);
  switch (EScalar3D_mode)
  {
  case PO_COLOR :
    if (n==TRIANGLE)
    {
      if (PlotColorTriangle3D(theElement,x,Poly[0],Poly[1],Poly[2],EScalar3D_depth,&theDO)) return (1);
    }
    else
    {
      if (PlotColorQuadrilateral3D(theElement,x,Poly[0],Poly[1],Poly[2],Poly[3],EScalar3D_depth,&theDO)) return (1);
    }
    break;
  case PO_CONTOURS_EQ :
    if (n==TRIANGLE)
    {
      if (PlotContourTriangle3D(theElement,x,Poly[0],Poly[1],Poly[2],EScalar3D_depth,&theDO)) return (1);
    }
    else
    {
      if (PlotContourQuadrilateral3D(theElement,x,Poly[0],Poly[1],Poly[2],Poly[3],EScalar3D_depth,&theDO)) return (1);
    }
    break;
  default :
    return (1);
  }

  DO_2c(theDO) = DO_NO_INST;
  DO_2C(range) = EScalar3D_minValue; DO_inc(range);
  DO_2C(range) = EScalar3D_maxValue;


  return (0);
}

/****************************************************************************/
/*
   FindRasterPoints3D - Find rasterpoints in 3D

   SYNOPSIS:
   static INT FindRasterPoints3D (COORD RasterSize, COORD_VECTOR *Polygon,
   INT Number, COORD_VECTOR *RasterPoints, INT *RPNumber);

   PARAMETERS:
   .  Rastersize -
   .  Polygon -
   .  Number -
   .  RasterPoints -
   .  RPNumber -

   DESCRIPTION:
   This function finds rasterpoints in 3D.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
 */
/****************************************************************************/


static INT FindRasterPoints3D (COORD RasterSize, COORD_VECTOR *Polygon, INT Number, COORD_VECTOR *RasterPoints, INT *RPNumber)
{
  INT i, j, k, i0, i1, j0, j1, c0, c1;
  COORD xmin, xmax, ymin, ymax;
  COORD diff[MAX_POINTS_OF_POLY][2], test[2];

  *RPNumber = 0;
  if (Number<2) return (0);

  xmin = ymin = MAX_C;
  xmax = ymax = -MAX_C;
  for (i=0; i<Number; i++)
  {
    xmin = MIN(xmin,Polygon[i][0]);
    xmax = MAX(xmax,Polygon[i][0]);
    ymin = MIN(ymin,Polygon[i][1]);
    ymax = MAX(ymax,Polygon[i][1]);
    diff[i][0] = Polygon[(i+1)%Number][0] - Polygon[i][0];
    diff[i][1] = Polygon[(i+1)%Number][1] - Polygon[i][1];
  }
  i0 = (INT)ceil(xmin/RasterSize);
  i1 = (INT)floor(xmax/RasterSize);
  j0 = (INT)ceil(ymin/RasterSize);
  j1 = (INT)floor(ymax/RasterSize);

  for (i=i0; i<=i1; i++)
    for (j=j0; j<=j1; j++)
    {
      c0 = c1 = 0;
      for (k=0; k<Number; k++)
      {
        test[0] = RasterSize*(COORD)(i) - Polygon[k][0];
        test[1] = RasterSize*(COORD)(j) - Polygon[k][1];
        if (diff[k][0]*test[1]>=diff[k][1]*test[0]) c0++;
        if (diff[k][0]*test[1]<=diff[k][1]*test[0]) c1++;
      }
      if (c0==Number || c1==Number)
      {
        RasterPoints[*RPNumber][0] = RasterSize*(COORD)(i);
        RasterPoints[*RPNumber][1] = RasterSize*(COORD)(j);
        RasterPoints[(*RPNumber)++][2] = 0.0;
      }
      if (*RPNumber==RASTERPOINTS_MAX)
        return (0);
    }

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  EW_EVector3D													*/
/*																			*/
/* Purpose:   evaluate elements for vector drawing							*/
/*																			*/
/* Input:	  PICTURE *thePicture, WORK *theWork							*/
/*																			*/
/* Return:	  INT 0: ok                                                                                                     */
/*			  INT 1: an error occurred										*/
/*																			*/
/****************************************************************************/

static INT EW_EVector3D (ELEMENT *theElement, DRAWINGOBJ *theDO)
{
  INT i, n, NodeOrder, nr;
  COORD_VECTOR LocalCoord, Poly[MAX_POINTS_OF_POLY], Poly2[MAX_POINTS_OF_POLY], RasterPoint[RASTERPOINTS_MAX];
  COORD *x[MAX_CORNERS_OF_ELEM], z[MAX_CORNERS_OF_ELEM];
  COORD scprd, norm, value;
  long Color;
  DOUBLE min, max;
  DOUBLE_VECTOR Arrow;

  /* get node order */
  NodeOrder = NORDER(theElement);

  /* get coordinates of corners of the element and their z coordinates in cut system */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));
    V3_TRAFO4_SC(x[i],CutTrafo,z[i])
  }

  /* determine polygon being intersection of element wth cut plane */
  if (GetPolyElemISCutPlane(x,z,NodeOrder,Poly,&n))
    return (1);
  if (n<=2)
  {
    DO_2c(theDO) = DO_NO_INST;
    return (0);
  }

  /* get arrows with rastersize (transform to cutsystem and back) */
  for (i=0; i<n; i++)
    V3_TRAFOM4_V3(Poly[i],CutTrafo,Poly2[i])
    if (FindRasterPoints3D(EVector_rastersize,Poly2,n,RasterPoint,&nr)) return (1);
  for (i=0; i<nr; i++)
  {
    V3_TRAFOM4_V3(RasterPoint[i],InvCutTrafo,Arrow)
    V3_COPY(Arrow,RasterPoint[i])
  }

  /* handle arrows */
  min = MAX_D; max = -MAX_D;
  for (i=0; i<nr; i++)
  {
    if (GlobalToLocal3d(x,RasterPoint[i],LocalCoord)) return (1);
    (*EVector_EvalFct)(theElement,x,LocalCoord,Arrow);
    V3_SCALE(EVector_V2L_factor,Arrow)

    /* find color and size of arrow, define its endpoint on the cutplane */
    V3_SCALAR_PRODUCT(Arrow,CUT_CutNormal,scprd)
    V3_EUKLIDNORM(Arrow,norm)

    if (norm!=0.0) value = scprd/norm;
    else value = 0.0;
    if (EVector3D_projectvector==YES && norm>SMALL_C)
    {
      V3_LINCOMB(1.0,Arrow,-scprd,CUT_CutNormal,Arrow)
      V3_EUKLIDNORM(Arrow,norm)
    }
    max = MAX(max,norm); min = MIN(min,norm);
    if ((norm>EVector_rastersize*EVector_CutLenFactor) && EVector_cutvector)
    {
      Color = EVector_ColorCut;
      V3_SCALE(EVector_rastersize*EVector_CutLenFactor/norm,Arrow)
    }
    else
    {
      Color = EVector3D_V2C_factor*value + EVector3D_V2C_offset;
      Color = MIN(Color,WOP_OutputDevice->spectrumEnd);
      Color = MAX(Color,WOP_OutputDevice->spectrumStart);
    }
    V3_ADD(RasterPoint[i],Arrow,Arrow);

    /* draw arrow */
    DO_2c(theDO) = DO_ARROW; DO_inc(theDO)
    DO_2l(theDO) = Color; DO_inc(theDO);
    V3_COPY(RasterPoint[i],DO_2Cp(theDO)); DO_inc_n(theDO,3);
    V3_COPY(Arrow,DO_2Cp(theDO)); DO_inc_n(theDO,3);
  }

  /* store range */
  DO_2c(theDO) = DO_RANGE; DO_inc(theDO)
  DO_2C(theDO) = min/EVector_V2L_factor; DO_inc(theDO);
  DO_2C(theDO) = max/EVector_V2L_factor; DO_inc(theDO);

  DO_2c(theDO) = DO_NO_INST;

  return (0);
}

#endif

/****************************************************************************/
/*D
   DrawPictureFrame - Draw a frame around picture

   SYNOPSIS:
   INT DrawPictureFrame (PICTURE *thePicture, INT mode);

   PARAMETERS:
   .  thePicture - draw a frame of the picture
   .  mode - mode of the frame

   DESCRIPTION:
   This function draw a frame around 'thePicture'. The mode determines the color of the
   frame. 'mode'==WOP_ACTIVE results in a orange frame, 'mode'==WOP_NOT_ACTIVE in a black
   and 'mode'==WOP_WORKING in a red one.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
   D*/
/****************************************************************************/

INT DrawPictureFrame (PICTURE *thePicture, INT mode)
{
  OUTPUTDEVICE *theOD;
  long color;
  COORD_POINT p[5];

  /* prepare graph for plot */
  if (PrepareGraph (thePicture)) return (1);

  /* set color */
  theOD  = PIC_OUTPUTDEV(thePicture);
  switch (mode)
  {
  case WOP_ACTIVE :
    color = theOD->orange;
    break;
  case WOP_NOT_ACTIVE :
    color = theOD->black;
    break;
  case WOP_WORKING :
    color = theOD->red;
    break;
  }

  /* plot invalidMode */
  UgSetLineWidth(1);
  UgSetColor(color);
  p[0].x = PIC_GLL(thePicture)[0]; p[0].y = PIC_GLL(thePicture)[1];
  p[1].x = PIC_GUR(thePicture)[0]; p[1].y = PIC_GLL(thePicture)[1];
  p[2].x = PIC_GUR(thePicture)[0]; p[2].y = PIC_GUR(thePicture)[1];
  p[3].x = PIC_GLL(thePicture)[0]; p[3].y = PIC_GUR(thePicture)[1];
  p[4].x = PIC_GLL(thePicture)[0]; p[4].y = PIC_GLL(thePicture)[1];
  UgPolyLine(p,5);

  return (0);
}

/****************************************************************************/
/*D
   ErasePicture	- Erase area of the Picture

   SYNOPSIS:
   INT ErasePicture (PICTURE *thePicture);

   PARAMETERS:
   .  thePicture - to be erased

   DESCRIPTION:
   This function erases area of the Picture.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
   D*/
/****************************************************************************/

INT ErasePicture (PICTURE *thePicture)
{
  COORD_POINT p[4];

  if (PrepareGraph (thePicture)) return (1);
  p[0].x = PIC_GLL(thePicture)[0]; p[0].y = PIC_GLL(thePicture)[1];
  p[1].x = PIC_GUR(thePicture)[0]; p[1].y = PIC_GLL(thePicture)[1];
  p[2].x = PIC_GUR(thePicture)[0]; p[2].y = PIC_GUR(thePicture)[1];
  p[3].x = PIC_GLL(thePicture)[0]; p[3].y = PIC_GUR(thePicture)[1];
  UgErasePolygon(p,4);

  return (0);
}

/****************************************************************************/
/*D
   WorkOnPicture - Work on picture

   SYNOPSIS:
   INT WorkOnPicture (PICTURE *thePicture, WORK *theWork);

   PARAMETERS:
   .  thePicture - the picture to work on
   .  theWork - the work to be performed

   DESCRIPTION:
   This function executes the specified 'WORK' on the specified 'PICTURE'. The 'PICTURE'
   has to have the 'status' 'ACTIVE' (completely initialized). An attemp to perform a
   'WORK' which is not executable results in an output on ug shell "action not executable"
   without errormessage.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
   D*/
/****************************************************************************/

INT WorkOnPicture (PICTURE *thePicture, WORK *theWork)
{
  INT i, end;

#ifdef __DO_HEAP_USED__
  char buffer[128];

  /* inits */
  Heap_Used_Min = MAX_I;
  Heap_Used_Max = -MAX_I;
#endif

  if (thePicture==NULL || theWork==NULL) return (1);
  WOP_Picture = thePicture;
  WOP_ViewedObj = PIC_VO(WOP_Picture);
  if (VO_STATUS(WOP_ViewedObj) != ACTIVE)
  {
    UserWrite("PlotObject and View have to be initialized\n");
    return (0);
  }
  WOP_Work                                        = theWork;
  WOP_OutputDevice                        = UGW_OUTPUTDEV(PIC_UGW(WOP_Picture));
  WOP_PlotObjHandling             = (PLOTOBJHANDLING*)PO_POT(PIC_PO(WOP_Picture));
  WOP_MG                                          = PO_MG(PIC_PO(WOP_Picture));
  if (WOP_MG == NULL) return (1);
  WOP_ViewDim                             = PO_DIM(PIC_PO(WOP_Picture));
  if (WOP_ViewDim == NOT_DEFINED) return (1);

  /* if FINDWORK: is plot valid? */
  if (W_ISSELECTWORK(WOP_Work))
    if (PIC_VALID(WOP_Picture) == NO)
    {
      UserWrite("cannot execute find-work: picture is not valid\n");
      return (0);
    }

  /* build transformation */
  if (BuildObsTrafo(WOP_Picture))
  {
    UserWrite("cannot build transformation\n");
    return (1);
  }

  /* activate low level grahic */
  if (PrepareGraph(WOP_Picture))
  {
    UserWrite("cannot activate low level graphic\n");
    return (1);
  }

  if (POH_NBCYCLES(WOP_PlotObjHandling,W_ID(WOP_Work)) <= 0)
  {
    UserWrite("action not executable\n");
    return (0);
  }

  /* clear if if DRAW_WORK */
  if (W_ID(theWork) == DRAW_WORK)
  {
    if (PO_CBD(PIC_PO(WOP_Picture)) == YES)
      if (ErasePicture(WOP_Picture))
        return (1);
    if (DrawPictureFrame(WOP_Picture,WOP_WORKING)) return (1);
  }

  for (i=0; i<POH_NBCYCLES(WOP_PlotObjHandling,W_ID(WOP_Work)); i++)
  {
    WOP_WorkProcs = POH_WORKPROGS(WOP_PlotObjHandling,W_ID(WOP_Work),i);
    WOP_WorkMode = WP_WORKMODE(WOP_WorkProcs);
    switch (WOP_WorkMode)
    {
    case ELEMENTWISE :

      /* order elements if */
      if (WOP_ViewDim == TYPE_3D)
      {
                                        #ifdef __TWODIM__
        if (OrderElements_2D(WOP_MG,WOP_ViewedObj))
                                        #endif
                                        #ifdef __THREEDIM__
        if (OrderElements_3D(WOP_MG,WOP_ViewedObj))
                                        #endif
        {
          UserWrite("ording of elements failed\n");
          return (1);
        }
      }

      /* set execution functions */
      WOP_GEN_PreProcessProc          =       WP_ELEMWISE(WOP_WorkProcs)->EW_PreProcessProc;
      WOP_EW_GetFirstElementProc      = (*WP_ELEMWISE(WOP_WorkProcs)->EW_GetFirstElementProcProc)(WOP_ViewedObj);
      WOP_EW_GetNextElementProc       = (*WP_ELEMWISE(WOP_WorkProcs)->EW_GetNextElementProcProc)(WOP_ViewedObj);
      WOP_EW_EvaluateProc             =       WP_ELEMWISE(WOP_WorkProcs)->EW_EvaluateProc;
      WOP_GEN_ExecuteProc             =       WP_ELEMWISE(WOP_WorkProcs)->EW_ExecuteProc;
      WOP_GEN_PostProcessProc         =       WP_ELEMWISE(WOP_WorkProcs)->EW_PostProcessProc;
      if (WOP_EW_EvaluateProc==NULL || WOP_GEN_ExecuteProc==NULL) return (1);

      /* work */
      if (WOP_GEN_PreProcessProc!=NULL)
        if ((*WOP_GEN_PreProcessProc)(WOP_Picture,WOP_Work))
          break;
      for (WOP_Element=(*WOP_EW_GetFirstElementProc)(WOP_MG,0,WOP_MG->currentLevel); WOP_Element!=NULL; WOP_Element=(*WOP_EW_GetNextElementProc)(WOP_Element))
      {
        if ((*WOP_EW_EvaluateProc)(WOP_Element,WOP_DrawingObject)) return (1);
        if ((*WOP_GEN_ExecuteProc)(WOP_DrawingObject)) return (1);
      }
      if (WOP_GEN_PostProcessProc!=NULL)
        if ((*WOP_GEN_PostProcessProc)(WOP_Picture,WOP_Work)) return (1);
      break;

    case NODEWISE :

      /* order elements if */
      if (WOP_ViewDim == TYPE_3D)
      {
        /* still missing */
      }

      /* set execution functions */
      WOP_GEN_PreProcessProc          =       WP_NODEWISE(WOP_WorkProcs)->NW_PreProcessProc;
      WOP_NW_GetFirstNodeProc         = (*WP_NODEWISE(WOP_WorkProcs)->NW_GetFirstNodeProcProc)(WOP_ViewedObj);
      WOP_NW_GetNextNodeProc          = (*WP_NODEWISE(WOP_WorkProcs)->NW_GetNextNodeProcProc)(WOP_ViewedObj);
      WOP_NW_EvaluateProc             =       WP_NODEWISE(WOP_WorkProcs)->NW_EvaluateProc;
      WOP_GEN_ExecuteProc             =       WP_NODEWISE(WOP_WorkProcs)->NW_ExecuteProc;
      WOP_GEN_PostProcessProc         =       WP_NODEWISE(WOP_WorkProcs)->NW_PostProcessProc;
      if (WOP_NW_EvaluateProc==NULL || WOP_GEN_ExecuteProc==NULL) return (1);

      /* work */
      if (WOP_GEN_PreProcessProc!=NULL)
        if ((*WOP_GEN_PreProcessProc)(WOP_Picture,WOP_Work))
          break;
      for (WOP_Node=(*WOP_NW_GetFirstNodeProc)(WOP_MG,0,WOP_MG->currentLevel); WOP_Node!=NULL; WOP_Node=(*WOP_NW_GetNextNodeProc)(WOP_Node))
      {
        if ((*WOP_NW_EvaluateProc)(WOP_Node,WOP_DrawingObject)) return (1);
        if ((*WOP_GEN_ExecuteProc)(WOP_DrawingObject)) return (1);
      }
      if (WOP_GEN_PostProcessProc!=NULL)
        if ((*WOP_GEN_PostProcessProc)(WOP_Picture,WOP_Work)) return (1);
      break;

    case VECTORWISE :

      /* set execution functions */
      WOP_GEN_PreProcessProc          =       WP_VECTORWISE(WOP_WorkProcs)->VW_PreProcessProc;
      WOP_VW_GetFirstVectorProc       = (*WP_VECTORWISE(WOP_WorkProcs)->VW_GetFirstVectorProcProc)(WOP_ViewedObj);
      WOP_VW_GetNextVectorProc        = (*WP_VECTORWISE(WOP_WorkProcs)->VW_GetNextVectorProcProc)(WOP_ViewedObj);
      WOP_VW_EvaluateProc             =       WP_VECTORWISE(WOP_WorkProcs)->VW_EvaluateProc;
      WOP_GEN_ExecuteProc             =       WP_VECTORWISE(WOP_WorkProcs)->VW_ExecuteProc;
      WOP_GEN_PostProcessProc         =       WP_VECTORWISE(WOP_WorkProcs)->VW_PostProcessProc;
      if (WOP_VW_EvaluateProc==NULL || WOP_GEN_ExecuteProc==NULL) return (1);

      /* work */
      if (WOP_GEN_PreProcessProc!=NULL)
        if ((*WOP_GEN_PreProcessProc)(WOP_Picture,WOP_Work))
          break;
      for (WOP_Vector=(*WOP_VW_GetFirstVectorProc)(WOP_MG,0,WOP_MG->currentLevel); WOP_Vector!=NULL; WOP_Vector=(*WOP_VW_GetNextVectorProc)(WOP_Vector))
      {
        if ((*WOP_VW_EvaluateProc)(WOP_Vector,WOP_DrawingObject)) return (1);
        if ((*WOP_GEN_ExecuteProc)(WOP_DrawingObject)) return (1);
      }
      if (WOP_GEN_PostProcessProc!=NULL)
        if ((*WOP_GEN_PostProcessProc)(thePicture,WOP_Work)) return (1);
      break;

    case EXTERN :

      /* set execution functions */
      WOP_GEN_PreProcessProc                  = WP_EXTERNWISE(WOP_WorkProcs)->EXT_PreProcessProc;
      WOP_EXT_EvaluateProc                    = WP_EXTERNWISE(WOP_WorkProcs)->EXT_EvaluateProc;
      WOP_GEN_ExecuteProc                     = WP_EXTERNWISE(WOP_WorkProcs)->EXT_ExecuteProc;
      WOP_GEN_PostProcessProc                 = WP_EXTERNWISE(WOP_WorkProcs)->EXT_PostProcessProc;
      if (WOP_EXT_EvaluateProc==NULL || WOP_GEN_ExecuteProc==NULL)
      {
        UserWrite("evaluation or execution procedure is missing\n");
        return (1);
      }

      /* work */
      if (WOP_GEN_PreProcessProc!=NULL)
        if ((*WOP_GEN_PreProcessProc)(thePicture,WOP_Work))
          break;
      end = 0;
      while (!end)
      {
        if ((*WOP_EXT_EvaluateProc)(WOP_DrawingObject,&end)) return (1);
        if ((*WOP_GEN_ExecuteProc)(WOP_DrawingObject)) return (1);
      }
      if (WOP_GEN_PostProcessProc!=NULL)
        if ((*WOP_GEN_PostProcessProc)(thePicture,WOP_Work)) return (1);
      break;

    default :
      return (1);
    }
  }

  /* may be picture is valid now */
  if (W_ID(theWork) == DRAW_WORK)
    PIC_VALID(WOP_Picture) = YES;

  /* flush cash */
  UgFlushCash();

  /* print heap used */
#ifdef __DO_HEAP_USED__
  sprintf(buffer,"Heap_min = %d\nHeap_max = %d\n",(int)Heap_Used_Min,(int)Heap_Used_Max);
  UserWrite(buffer);
#endif

  return (0);
}

/****************************************************************************/
/*D
   DrawUgPicture - Draw the picture

   SYNOPSIS:
   INT DrawUgPicture (PICTURE *thePicture);

   PARAMETERS:
   .  thePicture - the picture to draw

   DESCRIPTION:
   This function draws the picture. It initializes a structure 'WORK' to be a
   'DRAW_WORK', and calls the function 'WOorkOnPicture'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
   D*/
/****************************************************************************/

INT DrawUgPicture (PICTURE *thePicture)
{
  WORK theWork;

  theWork.WorkID = DRAW_WORK;
  if (WorkOnPicture(thePicture,&theWork)) return (1);

  return (0);
}

/****************************************************************************/
/*D
   InitWOP - Initialization

   SYNOPSIS:
   INT InitWOP (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function enrols all 'PLOTOBJHANDLING's. First all ptr. to functions which
   perform the work are set and then the substructures 'PLOTOBJTYPE' of the
   'PLOTOBJHANDLING' are initialized by the call of 'InitPlotObjTypes'.
   All functions used for working on the 'PICTURE's are defined in this file.
   See 'PLOTOBJHANDLING', 'PLOTOBJTYPE'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if en arror occured.
   D*/
/****************************************************************************/

INT InitWOP (void)
{
  PLOTOBJHANDLING *thePOH;
  WORKPROCS *theWP;
  ELEMWISEWORK *theEWW;

        #ifdef __TWODIM__
  NODEWISEWORK *theNWW;

  /* create WorkHandling for 'Grid' */
  if ((thePOH=CreatePlotObjHandling ("Grid"))     == NULL) return (__LINE__);

  /* draw work */
  POH_NBCYCLES(thePOH,DRAW_WORK) = 3;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_PlotElements2D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ElementEval2D;
  theEWW->EW_ExecuteProc                                  = Draw2D;
  theEWW->EW_PostProcessProc                              = InvertElementSelection2D;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,1);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_PlotBndOfElem2D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_BndOfElemEval2D;
  theEWW->EW_ExecuteProc                                  = Draw2D;
  theEWW->EW_PostProcessProc                              = NULL;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,2);
  WP_WORKMODE(theWP) = NODEWISE;
  theNWW = WP_NODEWISE(theWP);
  theNWW->NW_PreProcessProc                               = NW_PreProcess_PlotNodes2D;
  theNWW->NW_GetFirstNodeProcProc                 = NW_GetFirstNode_hor_fw_up_Proc;
  theNWW->NW_GetNextNodeProcProc                  = NW_GetNextNode_hor_fw_up_Proc;
  theNWW->NW_EvaluateProc                                 = NW_NodesEval2D;
  theNWW->NW_ExecuteProc                                  = Draw2D;
  theNWW->NW_PostProcessProc                              = InvertNodeSelection2D;

  /* findrange work */
  POH_NBCYCLES(thePOH,FINDRANGE_WORK) = 0;

  /* selectnode work */
  POH_NBCYCLES(thePOH,SELECTNODE_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,SELECTNODE_WORK,0);
  WP_WORKMODE(theWP) = NODEWISE;
  theNWW = WP_NODEWISE(theWP);
  theNWW->NW_PreProcessProc                               = NW_PreProcess_SelectNode2D;
  theNWW->NW_GetFirstNodeProcProc                 = NW_GetFirstNode_hor_fw_up_Proc;
  theNWW->NW_GetNextNodeProcProc                  = NW_GetNextNode_hor_fw_up_Proc;
  theNWW->NW_EvaluateProc                                 = NW_NodesEval2D;
  theNWW->NW_ExecuteProc                                  = NW_SelectNode2D;
  theNWW->NW_PostProcessProc                              = NULL;

  /* selectelement work */
  POH_NBCYCLES(thePOH,SELECTELEMENT_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,SELECTELEMENT_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_SelectElement2D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ElementEval2D;
  theEWW->EW_ExecuteProc                                  = EW_SelectElement2D;
  theEWW->EW_PostProcessProc                              = NULL;

  /* markelement work */
  POH_NBCYCLES(thePOH,MARKELEMENT_WORK) = 0;

  /* insertnode work */
  POH_NBCYCLES(thePOH,INSERTNODE_WORK) = 0;

  /* movenode work */
  POH_NBCYCLES(thePOH,MOVENODE_WORK) = 0;

  /* insertbndnode work */
  POH_NBCYCLES(thePOH,INSERTBNDNODE_WORK) = 0;


  /* create WorkHandling for 'EScalar' */
  if ((thePOH=CreatePlotObjHandling ("EScalar"))     == NULL) return (__LINE__);

  /* draw work */
  POH_NBCYCLES(thePOH,DRAW_WORK) = 2;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EScalar2D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EScalar2D;
  theEWW->EW_ExecuteProc                                  = Draw2D;
  theEWW->EW_PostProcessProc                              = NULL;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,1);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_PlotBlackBnd2D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_BndOfElemEval2D;
  theEWW->EW_ExecuteProc                                  = Draw2D;
  theEWW->EW_PostProcessProc                              = NULL;

  /* findrange work */
  POH_NBCYCLES(thePOH,FINDRANGE_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,FINDRANGE_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EScalar2D_FR;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EScalar2D;
  theEWW->EW_ExecuteProc                                  = FindRange2D;
  theEWW->EW_PostProcessProc                              = GEN_PostProcess_Scalar_FR;


  /* selectnode work */
  POH_NBCYCLES(thePOH,SELECTNODE_WORK) = 0;

  /* selectelement work */
  POH_NBCYCLES(thePOH,SELECTELEMENT_WORK) = 0;

  /* markelement work */
  POH_NBCYCLES(thePOH,MARKELEMENT_WORK) = 0;

  /* insertnode work */
  POH_NBCYCLES(thePOH,INSERTNODE_WORK) = 0;

  /* movenode work */
  POH_NBCYCLES(thePOH,MOVENODE_WORK) = 0;

  /* insertbndnode work */
  POH_NBCYCLES(thePOH,INSERTBNDNODE_WORK) = 0;



  /* create WorkHandling for 'EVector' */
  if ((thePOH=CreatePlotObjHandling ("EVector"))  == NULL) return (__LINE__);

  /* draw work */
  POH_NBCYCLES(thePOH,DRAW_WORK) = 2;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_PlotBlackBnd2D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_BndOfElemEval2D;
  theEWW->EW_ExecuteProc                                  = Draw2D;
  theEWW->EW_PostProcessProc                              = NULL;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,1);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EVector2D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EVector2D;
  theEWW->EW_ExecuteProc                                  = Draw2D;
  theEWW->EW_PostProcessProc                              = NULL;

  /* findrange work */
  POH_NBCYCLES(thePOH,FINDRANGE_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,FINDRANGE_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EVector2D_FR;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EVector2D;
  theEWW->EW_ExecuteProc                                  = FindRange2D;
  theEWW->EW_PostProcessProc                              = GEN_PostProcess_Vector_FR;

  /* selectnode work */
  POH_NBCYCLES(thePOH,SELECTNODE_WORK) = 0;

  /* selectelement work */
  POH_NBCYCLES(thePOH,SELECTELEMENT_WORK) = 0;

  /* markelement work */
  POH_NBCYCLES(thePOH,MARKELEMENT_WORK) = 0;

  /* insertnode work */
  POH_NBCYCLES(thePOH,INSERTNODE_WORK) = 0;

  /* movenode work */
  POH_NBCYCLES(thePOH,MOVENODE_WORK) = 0;

  /* insertbndnode work */
  POH_NBCYCLES(thePOH,INSERTBNDNODE_WORK) = 0;

        #endif



        #ifdef __THREEDIM__

  /* allocate VSIDES, COUNT, NORDER and CUTMODE in element ctrl or element flag */
  if (AllocateControlEntry(ELEMENT_CW,VSIDES_LEN,&ce_VSIDES) != GM_OK)
    if (AllocateControlEntry(FLAG_CW,VSIDES_LEN,&ce_VSIDES) != GM_OK)
      return (__LINE__);
  if (AllocateControlEntry(ELEMENT_CW,NORDER_LEN,&ce_NORDER) != GM_OK)
    if (AllocateControlEntry(FLAG_CW,NORDER_LEN,&ce_NORDER) != GM_OK)
      return (__LINE__);
  if (AllocateControlEntry(ELEMENT_CW,COUNT_LEN,&ce_COUNT) != GM_OK)
    if (AllocateControlEntry(FLAG_CW,COUNT_LEN,&ce_COUNT) != GM_OK)
      return (__LINE__);
  if (AllocateControlEntry(ELEMENT_CW,CUTMODE_LEN,&ce_CUTMODE) != GM_OK)
    if (AllocateControlEntry(FLAG_CW,CUTMODE_LEN,&ce_CUTMODE) != GM_OK)
      return (__LINE__);

  /* create WorkHandling for 'Grid' */
  if ((thePOH=CreatePlotObjHandling ("Grid"))     == NULL) return (__LINE__);

  /* draw work */
  POH_NBCYCLES(thePOH,DRAW_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_PlotGrid3D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ElementEval3D;
  theEWW->EW_ExecuteProc                                  = Draw3D;
  theEWW->EW_PostProcessProc                              = NULL;

  /* findrange work */
  POH_NBCYCLES(thePOH,FINDRANGE_WORK) = 0;

  /* selectnode work */
  POH_NBCYCLES(thePOH,SELECTNODE_WORK) = 0;

  /* selectelement work */
  POH_NBCYCLES(thePOH,SELECTELEMENT_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,SELECTELEMENT_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_SelectElement3D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_bw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_bw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ElementEval3D;
  theEWW->EW_ExecuteProc                                  = EW_SelectElement3D;
  theEWW->EW_PostProcessProc                              = NULL;

  /* markelement work */
  POH_NBCYCLES(thePOH,MARKELEMENT_WORK) = 0;

  /* insertnode work */
  POH_NBCYCLES(thePOH,INSERTNODE_WORK) = 0;

  /* movenode work */
  POH_NBCYCLES(thePOH,MOVENODE_WORK) = 0;

  /* insertbndnode work */
  POH_NBCYCLES(thePOH,INSERTBNDNODE_WORK) = 0;




  /* create WorkHandling for 'EScalar' */
  if ((thePOH=CreatePlotObjHandling ("EScalar"))  == NULL) return (__LINE__);

  /* draw work */
  POH_NBCYCLES(thePOH,DRAW_WORK) = 3;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EScalar3D_BackGrid;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ElementEval3D;
  theEWW->EW_ExecuteProc                                  = Draw3D;
  theEWW->EW_PostProcessProc                              = NULL;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,1);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EScalar3D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EScalar3D;
  theEWW->EW_ExecuteProc                                  = Draw3D;
  theEWW->EW_PostProcessProc                              = NULL;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,2);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_CutBnd3D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ECutBnd3D;
  theEWW->EW_ExecuteProc                                  = Draw3D;
  theEWW->EW_PostProcessProc                              = NULL;

  /* findrange work */
  POH_NBCYCLES(thePOH,FINDRANGE_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,FINDRANGE_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EScalar3D_FR;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EScalar3D;
  theEWW->EW_ExecuteProc                                  = FindRange3D;
  theEWW->EW_PostProcessProc                              = GEN_PostProcess_Scalar_FR;

  /* selectnode work */
  POH_NBCYCLES(thePOH,SELECTNODE_WORK) = 0;

  /* selectelement work */
  POH_NBCYCLES(thePOH,SELECTELEMENT_WORK) = 0;

  /* markelement work */
  POH_NBCYCLES(thePOH,MARKELEMENT_WORK) = 0;

  /* insertnode work */
  POH_NBCYCLES(thePOH,INSERTNODE_WORK) = 0;

  /* movenode work */
  POH_NBCYCLES(thePOH,MOVENODE_WORK) = 0;

  /* insertbndnode work */
  POH_NBCYCLES(thePOH,INSERTBNDNODE_WORK) = 0;




  /* create WorkHandling for 'EVector' */
  if ((thePOH=CreatePlotObjHandling ("EVector"))  == NULL) return (__LINE__);

  /* draw work */
  POH_NBCYCLES(thePOH,DRAW_WORK) = 3;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EVector3D_BackGrid;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ElementEval3D;
  theEWW->EW_ExecuteProc                                  = Draw3D;
  theEWW->EW_PostProcessProc                              = NULL;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,1);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_CutBnd3D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_ECutBnd3D;
  theEWW->EW_ExecuteProc                                  = Draw3D;
  theEWW->EW_PostProcessProc                              = NULL;

  theWP = POH_WORKPROGS(thePOH,DRAW_WORK,2);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EVector3D;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EVector3D;
  theEWW->EW_ExecuteProc                                  = Draw3D;
  theEWW->EW_PostProcessProc                              = NULL;

  /* findrange work */
  POH_NBCYCLES(thePOH,FINDRANGE_WORK) = 1;

  theWP = POH_WORKPROGS(thePOH,FINDRANGE_WORK,0);
  WP_WORKMODE(theWP) = ELEMENTWISE;
  theEWW = WP_ELEMWISE(theWP);
  theEWW->EW_PreProcessProc                               = EW_PreProcess_EVector3D_FR;
  theEWW->EW_GetFirstElementProcProc              = EW_GetFirstElement_vert_fw_up_Proc;
  theEWW->EW_GetNextElementProcProc               = EW_GetNextElement_vert_fw_up_Proc;
  theEWW->EW_EvaluateProc                                 = EW_EVector3D;
  theEWW->EW_ExecuteProc                                  = FindRange3D;
  theEWW->EW_PostProcessProc                              = GEN_PostProcess_Vector_FR;

  /* selectnode work */
  POH_NBCYCLES(thePOH,SELECTNODE_WORK) = 0;

  /* selectelement work */
  POH_NBCYCLES(thePOH,SELECTELEMENT_WORK) = 0;

  /* markelement work */
  POH_NBCYCLES(thePOH,MARKELEMENT_WORK) = 0;

  /* insertnode work */
  POH_NBCYCLES(thePOH,INSERTNODE_WORK) = 0;

  /* movenode work */
  POH_NBCYCLES(thePOH,MOVENODE_WORK) = 0;

  /* insertbndnode work */
  POH_NBCYCLES(thePOH,INSERTBNDNODE_WORK) = 0;

        #endif

  /* set PlotObjTypes of PlotObjHandlings */
  if (InitPlotObjTypes()) return (__LINE__);

  return(0);
}
