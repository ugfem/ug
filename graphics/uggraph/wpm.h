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
#define PO_DIM(p)                               (((p)->theHead.thePlotObjType==NULL) ? (NOT_DEFINED) : ((p)->theHead.thePlotObjType->Dimension))

/****************************************************************************/
/*																			*/
/* Defines and Macros for VIEWEDOBJ                                                                             */
/*																			*/
/****************************************************************************/

#define VO_PO(p)                                (&((p)->thePlotObj))
#define VO_DIM(p)                               (PO_DIM(&((p)->thePlotObj)))
#define VO_STATUS(p)                    ((p)->status)
#define VO_PERSPECTIVE(p)               ((p)->perspective)
#define VO_VP(p)                                ((p)->ViewPoint)
#define VO_VT(p)                                ((p)->ViewTarget)
#define VO_PMP(p)                               ((p)->PlaneMidpoint)
#define VO_PXD(p)                               ((p)->PlaneXDir)
#define VO_PYD(p)                               ((p)->PlaneYDir)

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
#define UGW_CURRTOOL(p)                         ((p)->currTool)
#define UGW_NAME(p)                             ((p)->d.name)
#define UGW_VALID(p)                            ((p)->Valid)

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
  COORD PlaneNormal[3];                                         /* normal vector of the plane					*/
  COORD PlanePoint[3];                                          /* point on the plane							*/
};

/*---------------------------- PlotObj head ---------------------------------*/

struct PlotObjHead {                                            /* head of all PlotObjs                                                 */

  INT status;                                                           /* see above									*/
  struct PlotObjType *thePlotObjType;           /* type of PlotObj								*/
  MULTIGRID *theMG;                                                     /* the data base								*/
  COORD theMidPoint[3];                                         /* MidPoint of covering 2/3D sphere                     */
  COORD theRadius;                                                      /* radius of covering 2/3D sphere				*/
  INT clearBeforeDraw;                                          /* YES or NO									*/
  char name[NAMESIZE];
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
  SYMBOL *Matrix;                                                       /* matrix symbol iff							*/
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
  COORD RasterSize;                                                     /* size of raster used for arrows				*/
  INT CutVectors;                                                        /* YES or NO									*/
  COORD CutLenFactor;                                                   /* vector will be cut if longer then                    */
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
  INT WhichElem;                                                        /* see above									*/
  INT ElemColored;                                                      /* YES or NO									*/
};

struct VecMatPlotObj2D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of vector-matrix graph */
  INT Marker;                                                                   /* YES or NO									*/
  INT Type[MAXVECTORS];                                         /* YES or NO									*/
  INT Connections;                                                      /* YES or NO									*/
  INT Extra;                                                                    /* YES or NO									*/
  INT Idx;                                                                      /* YES or NO									*/
  INT Order;                                                                    /* YES or NO									*/
  INT Dependency;                                                       /* YES or NO									*/
  INT ConnectVectors;                                                   /* YES or NO									*/
  INT Boundary;                                                         /* YES or NO									*/
  SYMBOL *vs;                                                                   /* NULL or vector symbol						*/
  SYMBOL *ms;                                                                   /* NULL or matrix symbol						*/
};

struct LinePlotObj2D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 2D-View of line field */
  EVALUES *EvalFct;                                                     /* evaluation proceedure						*/
  DOUBLE min, max;                                                      /* range										*/
  COORD_VECTOR left, right;                                     /* line in 2D physical space					*/
  INT depth;                                                                    /* depth of recoursive subdevision of elements	*/
  DOUBLE color;                                                         /* value between 0 and 1 specifiing the color   */
  COORD aspectratio;                                                    /* ratio of the picture							*/
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
  struct Cut theCut;                                                    /* description of the cut						*/
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
  struct Cut theCut;                                                    /* description of the cut						*/
  EVECTOR *EvalFct;                                                     /* evaluation proceedure						*/
  DOUBLE max;                                                           /* range										*/
  COORD RasterSize;                                                     /* size of raster used for arrows				*/
  INT CutVector;                                                        /* YES or NO									*/
  INT ProjectVector;                                                    /* YES or NO									*/
  COORD CutLenFactor;                                                   /* vector will be cut if longer then                    */
  /*	'CutLenFactor*RasterSize'					*/
};

struct GridPlotObj3D {

  struct PlotObjHead theHead;                           /* the head                                                                     */

  /* data for 3D-View of grid */
  DOUBLE ShrinkFactor;                                          /* YES or NO									*/
        #ifdef ModelP
  DOUBLE PartShrinkFactor;                                      /* YES or NO									*/
        #endif
  struct Cut theCut;                                                    /* description of the cut						*/
  INT ElemColored;                                                      /* YES or NO									*/
  INT WhichElem;                                                        /* see above									*/
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
  struct GridPlotObj3D theGpo;
#endif

};

struct ViewedObj {

  /* the object */
  union PlotObj thePlotObj;                                     /* plot object									*/

  /* the view */
  INT status;                                                           /* see above									*/
  INT perspective;                                                      /* YES or NO									*/
  COORD ViewPoint[3];                                           /* Observer Stand								*/
  COORD ViewTarget[3];                                          /* View target point							*/
  COORD PlaneMidpoint[3];                                       /* description of projection plane (the infinite*/
  COORD PlaneXDir[3], PlaneYDir[3];             /* extension touches the ViewTarget)			*/
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
  INT currTool;                                                                 /* current tool chosen											*/
  INT Global_LL[2], Global_UR[2];                               /* size of Ugwindow w.r.t. parent i.e.the ug-screen pixelspace	*/
  INT Local_LL[2], Local_UR[2];                                 /* real pixelrange of UgWindow, given by LowerLeft, UpperRight	*/
};

/*----------- typedef for functions ----------------------------------------*/

typedef INT (*DispPlotObjProcPtr)(union PlotObj *thePlotObj);
typedef INT (*SetPlotObjProcPtr)(union PlotObj *thePlotObj, INT argc, char **argv);

struct PlotObjType {

  ENVVAR v;                                                                             /* envitem of the UgWindow										*/

  INT Dimension;                                                                /* see above													*/
  SetPlotObjProcPtr SetPlotObjProc;                             /* proc for initializing the PlotObj							*/
  DispPlotObjProcPtr DispPlotObjProc;                   /* proc for displaying the PlotObj								*/
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

/* copy ViewedObject */
/*INT                   CopyViewedObjToPicture			(PICTURE *thePicture, VIEWEDOBJ *theViewedObj);*/

/* initializing/changing view */
INT                     SetView                                                 (PICTURE *thePicture, const COORD *viewPoint, const COORD *targetPoint, const COORD *xAxis, const INT *perspective);
INT                             PrintViewSettings                               (const PICTURE *thePicture);
INT                     DisplayViewOfViewedObject               (const PICTURE *thePicture);
INT                     Walk                                                    (PICTURE *thePicture, const COORD *vrsDelta);
INT                     RunAroundTargetPoint                    (PICTURE *thePicture, COORD vrsDirectionAngle, COORD vrsAngle);
INT                     Zoom                                                    (PICTURE *thePicture, COORD factor);
INT                     DragProjectionPlane                     (PICTURE *thePicture, COORD vrsDeltaX, COORD vrsDeltaY);
INT                     RotateProjectionPlane                   (PICTURE *thePicture, COORD vrsAngle);

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
