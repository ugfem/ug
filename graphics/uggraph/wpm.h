// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  wpm.h                                                                                                                 */
/*																			*/
/* Purpose:   defines data structure for wpm.c								*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Institut fuer Computeranwendungen                                                     */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de                                                */
/*																			*/
/* History:   8.12.94 begin, ug3-version									*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __WPM__
#define __WPM__

#include "compiler.h"
#include "devices.h"
#include "ugenv.h"
#include "gm.h"
#include "num.h"

/****************************************************************************/
/*																			*/
/* defines in the arbitrary order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Defines and Macros, misc                                                                                             */
/*																			*/
/****************************************************************************/

/* formats for display routines */
#define DISPLAY_PO_FORMAT_SS            "%-15.12s = %-25.22s\n"
#define DISPLAY_PO_FORMAT_SF            "%-15.12s = %-7.4g\n"
#define DISPLAY_PO_FORMAT_SFF           "%-15.12s = %-7.4g  %-7.4g\n"
#define DISPLAY_PO_FORMAT_SFFF          "%-15.12s = %-7.4g  %-7.4g  %-7.4g\n"
#define DISPLAY_PO_FORMAT_SI            "%-15.12s = %-2d\n"
#define DISPLAY_PO_FORMAT_SII           "%-15.12s = %-2d  %-2d\n"
#define DISPLAY_PO_FORMAT_SIII          "%-15.12s = %-2d  %-2d  %-2d\n"

/* dimension of PLOTOBJ, VIEWEDOBJ */
#define NOT_DEFINED                     0
#define TYPE_2D                                 1
#define TYPE_3D                                 2

/* status field of VIEW, CUT, PLOTOBJTYPE */
#define NOT_INIT                                0
#define NOT_ACTIVE                              1
#define ACTIVE                                  2

/* info box status in UGWINDOW */
#define BOX_INVALID             -2      /* indicates that info-box has to be redrawn*/
#define NO_INFO_AVAILABLE       -3      /* no information avaiable to be printed	*/
#define MOUSE_IN_CURR_PIC       -4      /* mouse in current picture of active graphw*/
#define MOUSE_OUT_CURR_PIC      -5      /* mouse in current picture of active graphw*/
#define STATIC_TEXT                     -6      /* static text instead of dynamic info		*/

/****************************************************************************/
/*																			*/
/* Defines and Macros for CUT												*/
/*																			*/
/****************************************************************************/

#define CUT_STATUS(p)                   ((p)->status)
#define CUT_PN(p)                               ((p)->PlaneNormal)
#define CUT_PP(p)                               ((p)->PlanePoint)

/****************************************************************************/
/*																			*/
/* Defines and Macros for PLOTOBJ											*/
/*																			*/
/****************************************************************************/

#define EXT_DESC_SIZE                   200

#define PO_NO                                   0
#define PO_COPY                                 1
#define PO_IRR                                  2
#define PO_REG                                  3
#define PO_ALL                                  4

#define PO_COLOR                                0
#define PO_CONTOURS_EQ                  1

#define PO_MAXCONTOURS                  50

#define PO_STATUS(p)                    ((p)->theHead.status)
#define PO_CBD(p)                               ((p)->theHead.clearBeforeDraw)
#define PO_POT(p)                               ((p)->theHead.thePlotObjType)
#define PO_MIDPOINT(p)                  ((p)->theHead.theMidPoint)
#define PO_RADIUS(p)                    ((p)->theHead.theRadius)
#define PO_MG(p)                                ((p)->theHead.theMG)
#define PO_NAME(p)                              ((p)->theHead.name)
#define PO_USESCUT(p)                   ((p)->theHead.UsesCut)
#define PO_DIM(p)                               (((p)->theHead.thePlotObjType==NULL) ? (NOT_DEFINED) : ((p)->theHead.thePlotObjType->Dimension))

/****************************************************************************/
/*																			*/
/* Defines and Macros for VIEWEDOBJ                                                                             */
/*																			*/
/****************************************************************************/

#define VO_PO(p)                                (&((p)->thePlotObj))
#define VO_DIM(p)                               (PO_DIM(&((p)->thePlotObj)))
#define VO_STATUS(p)                    ((p)->status)
#define VO_CUT(p)                               (&((p)->theCut))
#define VO_PERSPECTIVE(p)               ((p)->perspective)
#define VO_VP(p)                                ((p)->ViewPoint)
#define VO_VT(p)                                ((p)->ViewTarget)
#define VO_PMP(p)                               ((p)->PlaneMidpoint)
#define VO_PXD(p)                               ((p)->PlaneXDir)
#define VO_PYD(p)                               ((p)->PlaneYDir)

#define VO_TRAFO(p)                             ((p)->ObsTrafo)
#define VO_INVTRAFO(p)                  ((p)->InvObsTrafo)

/****************************************************************************/
/*																			*/
/* Defines and Macros for PICTURE											*/
/*																			*/
/****************************************************************************/

#define PIC_UGW(p)                                      ((p)->theUgWindow)
#define PIC_OUTPUTDEV(p)                        ((p)->theUgWindow->theOutputDevice)
#define PIC_VO(p)                                       (&((p)->theViewedObj))
#define PIC_PO(p)                                       (&((p)->theViewedObj.thePlotObj))
#define PIC_MG(p)                                       ((p)->theViewedObj.thePlotObj.theHead.theMG)
#define PIC_VALID(p)                            ((p)->Valid)
#define PIC_GLL(p)                                      ((p)->Global_LL)
#define PIC_GUR(p)                                      ((p)->Global_UR)
#define PIC_NAME(p)                             ((p)->v.name)
#define PIC_POT(p)                                      ((p)->theViewedObj.thePlotObj.theHead.thePlotObjType)
#define PIC_SIGN_X(p)                           ((p)->sign_x)
#define PIC_SIGN_Y(p)                           ((p)->sign_y)

/****************************************************************************/
/*																			*/
/* Defines and Macros for UGWINDOW											*/
/*																			*/
/****************************************************************************/

#define UGW_OUTPUTDEV(p)                        ((p)->theOutputDevice)
#define UGW_IFWINDOW(p)                         ((p)->theIFWindow)
#define UGW_GLL(p)                                      ((p)->Global_LL)
#define UGW_GUR(p)                                      ((p)->Global_UR)
#define UGW_LLL(p)                                      ((p)->Local_LL)
#define UGW_LUR(p)                                      ((p)->Local_UR)
#define UGW_NPIC(p)                             ((p)->NbPictures)
#define UGW_NAME(p)                             ((p)->d.name)
#define UGW_VALID(p)                            ((p)->Valid)
#define UGW_CURRTOOL(p)                         ((p)->currTool)
#define UGW_CURRFUNC(p)                         ((p)->currFunc)
#define UGW_INFOTEXT(p)                         ((p)->info)
#define UGW_BOXSTATE(p)                         ((p)->InfoBoxState)

/****************************************************************************/
/*																			*/
/* Defines and Macros for PLOTOBJTYPE										*/
/*																			*/
/****************************************************************************/

#define POT_DIM(p)                                      ((p)->Dimension)

/****************************************************************************/
/*																			*/
/* structures and ...														*/
/*																			*/
/****************************************************************************/

/*----------- definition of structs ----------------------------------------*/

struct Cut {

  INT status;                                                           /* see above									*/
  DOUBLE PlaneNormal[3];                                        /* normal vector of the plane					*/
  DOUBLE PlanePoint[3];                                         /* point on the plane							*/
};

/*---------------------------- PlotObj head ---------------------------------*/

struct PlotObjHead {                                            /* head of all PlotObjs                                                 */

  INT status;                                                           /* see above									*/
  struct PlotObjType *thePlotObjType;           /* type of PlotObj								*/
  MULTIGRID *theMG;                                                     /* the data base								*/
  DOUBLE theMidPoint[3];                                        /* MidPoint of covering 2/3D sphere                     */
  DOUBLE theRadius;                                                     /* radius of covering 2/3D sphere				*/
  INT clearBeforeDraw;                                          /* YES or NO									*/
  char name[NAMESIZE];
  INT UsesCut;                                                          /* YES or NO									*/
};

/*----------- application dimension independent PlotObj ---------------------*/

struct MatrixPlotObj {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of matrix */
  MVALUES *EvalFct;                                                     /* evaluation proceedure						*/
  INT log;                                                                      /* use log of absolute value					*/
  INT rel;                                                                      /* values relative to diagonal entry			*/
  INT BV;                                                                       /* plot blockvector blocks						*/
  DOUBLE thresh;                                                        /* plot entries only if |.|>thresh				*/
  INT conn;                                                                     /* plot connections								*/
  INT extra;                                                                    /* plot extra connections						*/
  DOUBLE min, max;                                                      /* range										*/
  MATDATA_DESC *Matrix;                                         /* matrix                                                               */
  DOUBLE dash;                                                          /* length of the line segments in dashed lines */
  DOUBLE space;                                                         /* gap between line segments in dashed lines*/
};

struct ExternPlotObject {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for extern object */
  char ExternObjDesc[EXT_DESC_SIZE];                    /* description of extern plotobj				*/
};

/*----------- application dimension 2 PlotObjs ------------------------------*/

struct ElemScalarPlotObj2D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of elem scalar field */
  EVALUES *EvalFct;                                                     /* evaluation proceedure						*/
  DOUBLE min, max;                                                      /* range										*/
  INT mode;                                                                     /* COLOR or CONTOURS							*/
  INT PlotGrid;                                                         /* plot grid together with scalar field			*/
  INT depth;                                                                    /* depth of recoursive subdevision of elements	*/
  INT numOfContours;                                                    /* nb of contourlines if used					*/
  DOUBLE contValues[PO_MAXCONTOURS];                    /* contour values if used						*/
};

struct ElemVectorPlotObj2D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of elem vector field */
  EVECTOR *EvalFct;                                                     /* evaluation proceedure						*/
  INT PlotGrid;                                                         /* plot grid together with scalar field			*/
  DOUBLE max;                                                           /* range										*/
  DOUBLE RasterSize;                                                    /* size of raster used for arrows				*/
  INT CutVectors;                                                        /* YES or NO									*/
  DOUBLE CutLenFactor;                                                  /* vector will be cut if longer then                    */
  /*	'CutLenFactor*RasterSize'					*/
};

struct GridPlotObj2D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of grid */
  DOUBLE ShrinkFactor;                                          /* YES or NO									*/
        #ifdef ModelP
  DOUBLE PartShrinkFactor;                                      /* YES or NO									*/
        #endif
  INT PlotElemID;                                                       /* YES or NO									*/
  INT PlotRefMarks;                                                     /* YES or NO									*/
  INT PlotIndMarks;                                                     /* YES or NO									*/
  INT PlotNodeID;                                                       /* YES or NO									*/
  INT PlotNodes;                                                        /* YES or NO									*/
  INT PlotBoundary;                                                     /* YES or NO									*/
  INT PlotSubdomain;                                                    /* YES or NO									*/
  INT WhichElem;                                                        /* see above									*/
  INT ElemColored;                                                      /* YES or NO									*/
  VECDATA_DESC *FreeBnd;                                        /* global coords of new free boundary			*/
};

struct VecMatPlotObj2D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of vector-matrix graph */
  INT Marker;                                                                   /* YES or NO									*/
  INT Type[MAXVECTORS];                                         /* YES or NO									*/
  INT Connections;                                                      /* YES or NO									*/
  INT Extra;                                                                    /* YES or NO									*/
  INT Idx;                                                                      /* YES or NO									*/
  INT Part;                                                                     /* YES or NO									*/
  INT Order;                                                                    /* YES or NO									*/
  INT Dependency;                                                       /* YES or NO									*/
  INT ConnectVectors;                                                   /* YES or NO									*/
  INT Boundary;                                                         /* YES or NO									*/
  VECDATA_DESC *vd;                                                     /* NULL or vector                                                       */
  MATDATA_DESC *md;                                                     /* NULL or matrix                                                       */
};

struct LinePlotObj2D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of line field */
  EVALUES *EvalFct;                                                     /* evaluation proceedure						*/
  DOUBLE min, max;                                                      /* range										*/
  INT yLog;                                                                     /* draw y-axis in logarithmic scale				*/
  DOUBLE_VECTOR left, right;                                    /* line in 2D physical space					*/
  INT depth;                                                                    /* depth of recoursive subdevision of elements	*/
  DOUBLE color;                                                         /* value between 0 and 1 specifiing the color   */
  DOUBLE aspectratio;                                                   /* ratio of the picture							*/
  INT nHit;                                                                     /* # elements hit by the line					*/
  DOUBLE xmin;                                                          /* min intersection between grid and line	    */
  DOUBLE xmax;                                                          /* max intersection between grid and line	    */
};

/*----------- application dimension 3 PlotObjs -----------------------------*/

struct DomainPlotObj3D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 3D-View of domain */
  INT NbOfSteps;                                                        /* number of lines in each direction of patch	*/
};

struct ElemScalarPlotObj3D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 3D-View of elem scalar field */
  EVALUES *EvalFct;                                                     /* evaluation proceedure						*/
  DOUBLE min, max;                                                      /* range										*/
  INT mode;                                                                     /* COLOR or CONTOURS							*/
  INT depth;                                                                    /* depth of recoursive subdevision of elements	*/
  INT numOfContours;                                                    /* nb of contourlines if used					*/
  DOUBLE contValues[PO_MAXCONTOURS];                    /* contour values if used						*/
};

struct ElemVectorPlotObj3D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 3D-View of elem vector field */
  EVECTOR *EvalFct;                                                     /* evaluation proceedure						*/
  DOUBLE max;                                                           /* range										*/
  DOUBLE RasterSize;                                                    /* size of raster used for arrows				*/
  INT CutVector;                                                        /* YES or NO									*/
  INT ProjectVector;                                                    /* YES or NO									*/
  DOUBLE CutLenFactor;                                          /* vector will be cut if longer then                    */
  /*	'CutLenFactor*RasterSize'					*/
};

struct VecMatPlotObj3D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 3D-View of VECTOR-MATRIX-data */
  INT Type[MAXVECTORS];                                         /* which types if Vectors set					*/
  INT Idx;                                                                      /* YES or NO									*/
  VECDATA_DESC *vd;                                                     /* NULL or vector								*/
  MATDATA_DESC *md;                                                     /* NULL or matrix                                                       */
};

struct GridPlotObj3D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 3D-View of grid */
  DOUBLE ShrinkFactor;                                          /* YES or NO									*/
        #ifdef ModelP
  DOUBLE PartShrinkFactor;                                      /* YES or NO									*/
        #endif
  INT NodeMarkers;                                                      /* plot node markers (will only for shrink<1)	*/
  INT NodeIndex;                                                        /* plot node indices (only together with marker)*/
  INT Vectors;                                                          /* plot vectors (mutually exclusive with nodes)	*/
  INT VecIndex;                                                         /* plot vector indices							*/
  INT OType[MAXVOBJECTS];                                       /* which object types if Vectors set			*/
  INT ElemColored;                                                      /* YES or NO									*/
  INT WhichElem;                                                        /* see above									*/
  INT PlotSelection;                                                            /* see above									*/
};

union PlotObj {

  struct PlotObjHead theHead;

#ifdef __TWODIM__
  struct MatrixPlotObj theMpo;
  struct ExternPlotObject theExternObject;
  struct ElemScalarPlotObj2D theEspo;
  struct ElemVectorPlotObj2D theEvpo;
  struct GridPlotObj2D theGpo;
  struct VecMatPlotObj2D theVmo;
  struct LinePlotObj2D theLpo;
#endif

#ifdef __THREEDIM__
  struct MatrixPlotObj theMpo;
  struct ExternPlotObject theExternObject;
  struct DomainPlotObj3D theDpo;
  struct ElemScalarPlotObj3D theEspo;
  struct ElemVectorPlotObj3D theEvpo;
  struct VecMatPlotObj3D theVmo;
  struct GridPlotObj3D theGpo;
#endif

};

struct ViewedObj {

  /* the object */
  union PlotObj thePlotObj;                                     /* plot object									*/

  /* the view */
  INT status;                                                           /* see above									*/
  INT perspective;                                                      /* YES or NO									*/
  DOUBLE ViewPoint[3];                                          /* Observer Stand								*/
  DOUBLE ViewTarget[3];                                         /* View target point							*/
  DOUBLE PlaneMidpoint[3];                                      /* description of projection plane (the infinite*/
  DOUBLE PlaneXDir[3], PlaneYDir[3];            /* extension touches the ViewTarget)			*/

  DOUBLE ObsTrafo[16];
  DOUBLE InvObsTrafo[16];

  struct Cut theCut;
};

struct PICture {

  ENVVAR v;                                                                             /* envitem of the picture										*/

  /* specification of the picture */
  struct UgWindow *theUgWindow;                                 /* window of that canvas										*/
  INT Global_LL[2], Global_UR[2];                               /* size of picture w.r.t. parent i.e.the ugwindow local pixel sp*/
  INT sign_x, sign_y;                                                           /* 1 if right/up system, -1 if left/down system					*/
  INT Valid;                                                                            /* YES or NO													*/
  struct ViewedObj theViewedObj;                                /* the object as you see it                                                                     */
};

struct UgWindow {

  ENVDIR d;                                                                             /* envitem of the UgWindow										*/

  struct outputdevice *theOutputDevice;                 /* corresponding Output Device									*/
  WINDOWID theIFWindow;                                                 /* identification of interface window							*/
  INT Valid;                                                                            /* YES or NO													*/
  INT NbPictures;                                                               /* number of pictures for that ugwindow                                                 */
  INT Global_LL[2], Global_UR[2];                               /* size of Ugwindow w.r.t. parent i.e.the ug-screen pixelspace	*/
  INT Local_LL[2], Local_UR[2];                                 /* real pixelrange of UgWindow, given by LowerLeft, UpperRight	*/

  /* info and tool box */
  INT currTool;
  INT currFunc;
  INT InfoBoxState;
  char info[INFO_SIZE];
};

/*----------- typedef for functions ----------------------------------------*/

typedef INT (*DispPlotObjProcPtr)(union PlotObj *thePlotObj);
typedef INT (*SetPlotObjProcPtr)(union PlotObj *thePlotObj, INT argc, char **argv);

struct PlotObjType {

  ENVVAR v;                                                                             /* envitem of the UgWindow					*/

  INT Dimension;                                                                /* see above								*/
  SetPlotObjProcPtr SetPlotObjProc;                             /* proc for initializing the PlotObj		*/
  DispPlotObjProcPtr DispPlotObjProc;                   /* proc for displaying the PlotObj			*/
};

/****************************************************************************/
/*																			*/
/*					typedef for structs                                                                     */
/*																			*/
/****************************************************************************/

typedef struct UgWindow UGWINDOW;
typedef struct PICture PICTURE;
typedef union  PlotObj PLOTOBJ;
typedef struct ViewedObj VIEWEDOBJ;
typedef struct PlotObjType PLOTOBJTYPE;
typedef struct Cut CUT;
typedef struct PlotObjHead PO_HEAD;

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* create/dispose/first/next ... */
UGWINDOW           *CreateUgWindow                                      (OUTPUTDEVICE *theOutputDevice, const char *UgWindowName, INT x, INT y, INT width, INT height);
INT                     DisposeUgWindow                                 (UGWINDOW *theUgWindow);
UGWINDOW           *GetFirstUgWindow                            (void);
UGWINDOW           *GetNextUgWindow                             (const UGWINDOW *theUgWindow);
PICTURE            *CreatePicture                                       (const char *PictureName, UGWINDOW *theUgWindow, const INT *Global_LL, const INT *Global_UR);
INT                     DisposePicture                                  (PICTURE *thePicture);
PICTURE            *GetFirstPicture                             (const UGWINDOW *theUgWindow);
PICTURE            *GetNextPicture                                      (const PICTURE *thePicture);

/* validate/invalidate ... */
INT                     InvalidatePicture                               (PICTURE *thePicture);
INT                     InvalidatePicturesOfMG                  (MULTIGRID *theMG);
INT                     InvalidatePicturesOfUgWindow    (UGWINDOW *theUgW);
INT                     UpdateUgWindow                                  (UGWINDOW *theUgWindow, const PICTURE *EvalPicture);
INT                     InvalidateUgWindow                              (UGWINDOW *theUgW);
INT                             InvalidateUgWindowsOfMG                 (MULTIGRID *theMG);
void                    ResetToolBoxState                               (UGWINDOW *ugw);

/* copy ViewedObject */
/*INT                   CopyViewedObjToPicture			(PICTURE *thePicture, VIEWEDOBJ *theViewedObj);*/

/* initializing/changing view */
INT                     SetView                                                 (PICTURE *thePicture, const DOUBLE *viewPoint, const DOUBLE *targetPoint, const DOUBLE *xAxis, const INT *perspective,
                                                                                 INT RemoveCut, const DOUBLE *cutPoint, const DOUBLE *cutNormal);
INT                             CopyView                                                (const PICTURE *mypic, INT all, INT cut);
INT                             PrintViewSettings                               (const PICTURE *thePicture);
INT                     DisplayViewOfViewedObject               (const PICTURE *thePicture);
INT                     Walk                                                    (PICTURE *thePicture, const DOUBLE *vrsDelta);
INT                     RunAroundTargetPoint                    (PICTURE *thePicture, DOUBLE vrsDirectionAngle, DOUBLE vrsAngle);
INT                     Zoom                                                    (PICTURE *thePicture, DOUBLE factor);
INT                     DragProjectionPlane                     (PICTURE *thePicture, DOUBLE vrsDeltaX, DOUBLE vrsDeltaY);
INT                     RotateProjectionPlane                   (PICTURE *thePicture, DOUBLE vrsAngle);

/* operations on PlotObjs */
INT                     SpecifyPlotObjOfViewedObject    (PICTURE *thePicture, MULTIGRID *theMG, const char *PlotObjTypeName, INT argc, char **argv);
INT                     DisplayObject                                   (PLOTOBJ *thePlotObj);

/* miscellaneous */
UGWINDOW           *WinID2UgWindow                                      (WINDOWID theIFWindowID);
UGWINDOW           *GetUgWindow                                         (const char *name);
void                    ListUgWindow                                    (const UGWINDOW *theUgWindow, INT current);
PICTURE            *GetUgPicture                                        (const UGWINDOW *theUgWindow, const char *name);
void                    ListPicture                                             (const PICTURE *thePicture, INT current);
INT                             MovePictureToNewWindow                  (PICTURE *pic);
PICTURE            *Mouse2Picture                                       (const UGWINDOW *theUgWindow, INT *MousePosition);
PLOTOBJTYPE    *GetPlotObjType                                  (const char *PlotObjTypeName);
PLOTOBJTYPE    *CreatePlotObjType                               (const char *PlotObjTypeName, INT size);
PLOTOBJTYPE    *GetFirstPlotObjType                             (void);
PLOTOBJTYPE    *GetNextPlotObjType                              (const PLOTOBJTYPE *thePlotObj);
INT                     ResizeViewPlane                                 (VIEWEDOBJ *theVO, const INT *Pix_LL_old, const INT *Pix_UR_old, const INT *Pix_LL_new, const INT *Pix_UR_new);
void                    ListWindowPictureHeader                 (void);
INT                     InitPlotObjTypes                                (void);
INT                     InitWPM                                                 (void);


#endif
