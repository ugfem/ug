// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  shapes.h														*/
/*																			*/
/* Purpose:   header file for shape functions								*/
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   28.11.95 begin, ug version 3.1								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __SHAPES__
#define __SHAPES__

#ifndef __GM__
#include "gm.h"
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

#ifdef __TWODIM__
#define GlobalToLocal(n,c,e,l)           GlobalToLocal2d (n,c,e,l)
#endif
#ifdef __THREEDIM__
#define GlobalToLocal(n,c,e,l)           GlobalToLocal3d (c,e,l)
#endif

#ifdef __TWODIM__
#define DNDS(n,i,s,t,r)                 if (n==3)\
  {\
    switch (i)\
    {\
    case 0 : r=-1.0; break;\
    case 1 : r=1.0; break;\
    case 2 : r=0.0; break;\
    }\
  }\
  else if (n==4)\
  {\
    switch (i)\
    {\
    case 0 : r=-0.25*(1-t); break;\
    case 1 : r=0.25*(1-t); break;\
    case 2 : r=0.25*(1+t); break;\
    case 3 : r=-0.25*(1+t); break;\
    }\
  }

#define DNDT(n,i,s,t,r)             if (n==3)\
  {\
    switch (i)\
    {\
    case 0 : r=-1; break;\
    case 1 : r=0; break;\
    case 2 : r=1; break;\
    }\
  }\
  else if (n==4)\
  {\
    switch (i)\
    {\
    case 0 : r=-0.25*(1-s); break;\
    case 1 : r=-0.25*(1+s); break;\
    case 2 : r=0.25*(1+s); break;\
    case 3 : r=0.25*(1-s); break;\
    }\
  }
#endif

#ifdef __TWODIM__
#define CORNER_COORDINATES(e,n,x) {(n)=CORNERS_OF_ELEM((e));              \
                                   (x)[0]=CVECT(MYVERTEX(CORNER((e),0))); \
                                   (x)[1]=CVECT(MYVERTEX(CORNER((e),1))); \
                                   (x)[2]=CVECT(MYVERTEX(CORNER((e),2))); \
                                   if(CORNERS_OF_ELEM((e))==4)                \
                                     (x)[3]=CVECT(MYVERTEX(CORNER((e),3)));}
#define LOCAL_2_GLOBAL(x,local,global)                    \
  {(global)[0] = (1.0-(local)[0]-(local)[1])*(x)[0][0]     \
                 +(local)[0]*(x)[1][0] + (local)[1]*(x)[2][0];          \
   (global)[1] = (1.0-(local)[0]-(local)[1])*(x)[0][1]     \
                 +(local)[0]*(x)[1][1] + (local)[1]*(x)[2][1]; }
#endif

#ifdef __THREEDIM__
#define CORNER_COORDINATES(e,n,x) {(n)=CORNERS_OF_ELEM((e));               \
                                   (x)[0]=CVECT(MYVERTEX(CORNER((e),0)));  \
                                   (x)[1]=CVECT(MYVERTEX(CORNER((e),1)));  \
                                   (x)[2]=CVECT(MYVERTEX(CORNER((e),2)));  \
                                   (x)[3]=CVECT(MYVERTEX(CORNER((e),3)));}
/*
                                     if(CORNERS_OF_ELEM((e))==5)             \
                                      (x)[4]=CVECT(MYVERTEX(CORNER((e),4)));}\
   if(CORNERS_OF_ELEM((e))==8){             \
                                      (x)[4]=CVECT(MYVERTEX(CORNER((e),4))); \
                                      (x)[5]=CVECT(MYVERTEX(CORNER((e),5))); \
                                      (x)[6]=CVECT(MYVERTEX(CORNER((e),6))); \
                                      (x)[7]=CVECT(MYVERTEX(CORNER((e),7)));}}*/
#define LOCAL_2_GLOBAL(x,local,global)                                    \
  {(global)[0] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][0]          \
                 +(local)[0]*(x)[1][0] + (local)[1]*(x)[2][0] + (local)[2]*(x)[3][0];   \
   (global)[1] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][1]          \
                 +(local)[0]*(x)[1][1] + (local)[1]*(x)[2][1] + (local)[2]*(x)[3][1];   \
   (global)[1] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][2]          \
                 +(local)[0]*(x)[1][2] + (local)[1]*(x)[2][2] + (local)[2]*(x)[3][2]; }
#define LOCAL_TO_GLOBAL(A,B,C)          {C[0] = (1.0-B[0]-B[1]-B[2])*A[0][0]+B[0]*A[1][0]+B[1]*A[2][0]+B[2]*A[3][0];\
                                         C[1] = (1.0-B[0]-B[1]-B[2])*A[0][1]+B[0]*A[1][1]+B[1]*A[2][1]+B[2]*A[3][1];\
                                         C[2] = (1.0-B[0]-B[1]-B[2])*A[0][2]+B[0]*A[1][2]+B[1]*A[2][2]+B[2]*A[3][2];}
#endif

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

#ifdef __THREEDIM__
extern COORD HexRefCoord[MAX_CORNERS_OF_ELEM][DIM];
extern COORD TransfCoeff[MAX_CORNERS_OF_ELEM][DIM];
extern COORD CenterOfIntergrSurf[MAX_EDGES_OF_ELEM][DIM];
#endif

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

DOUBLE          GN              (INT n, INT i, COORD *local);
COORD      *LMP         (INT n);


#ifdef __TWODIM__
DOUBLE  N                               (INT n, INT i, DOUBLE s, DOUBLE t);
DOUBLE  dNds                    (INT n, INT i, DOUBLE s, DOUBLE t);
DOUBLE  dNdt                    (INT n, INT i, DOUBLE s, DOUBLE t);

INT     Derivatives     (INT n, const DOUBLE *px, const DOUBLE *py, DOUBLE ips, DOUBLE ipt, DOUBLE *dNdx, DOUBLE *dNdy, DOUBLE *detJ);
INT             Gradients               (INT n, const COORD **theCorners, DOUBLE ips, DOUBLE ipt, DOUBLE_VECTOR Gradient[MAX_CORNERS_OF_ELEM], DOUBLE *DetJ);

INT             L2GDerivative2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD *Derivative);
INT             LocalToGlobal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR GlobalCoord);
INT     GlobalToLocal2d (INT n,const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);
INT             specialGlobalToLocal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);
#endif

#ifdef __THREEDIM__
INT TransformGlobalToLocal3D(ELEMENT *theElement, COORD_VECTOR Global, COORD_VECTOR Local);
DOUBLE  N                   (const INT i, const COORD *LocalCoord);
INT     GlobalToLocal3d     (const COORD **Corners, const COORD *EvalPoint, COORD *LocalCoord);
INT     TetraDerivative     (ELEMENT *theElement, const COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM]);
INT     TetraVolume         (const COORD **theCorners, COORD *volume);
INT     FV_TetInfo          (const COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM]);
INT     Side_TetInfo            (COORD **theCorners, INT side, COORD_VECTOR Area, COORD_VECTOR GIP[3]);
INT     TetraSideNormals    (ELEMENT *theElement, COORD **theCorners, COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM]);
INT     TetMaxSideAngle     (ELEMENT *theElement, const COORD **theCorners, COORD *MaxAngle);
INT     TetAngleAndLength   (ELEMENT *theElement, const COORD **theCorners, COORD *Angle, COORD *Length);
INT     FV_AliTetInfo       (const COORD **CornerPoints, COORD_VECTOR Area[6], DOUBLE_VECTOR conv, COORD_VECTOR GIP[6], COORD_VECTOR LIP[6]);
INT     FV_TetInfo_for_conv (ELEMENT *theElement, const COORD **CornerPoints, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM], COORD_VECTOR conv);
INT     GFUIP               (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     GCUIP               (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     COPYIP              (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);

INT TransformCoefficients (ELEMENT *theElement);
INT JacobiMatrix (COORD_VECTOR alpha, COORD *Matrix);
DOUBLE JacobiMatrixDet (COORD *Matrix);
INT GradientsOfTransormation (COORD_VECTOR alpha, COORD_VECTOR grad, int i);
INT LocalToGlobalHEX (COORD_VECTOR Local, COORD_VECTOR Global);
INT GlobalToLocalHEX (COORD_VECTOR Global, COORD_VECTOR Local);
INT HexaDerivatives (COORD *alpha, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM]);
DOUBLE NiHex (int i, COORD *alpha);
DOUBLE NiHexDer (int i, int j, COORD *alpha);
INT FV_HexInfo (ELEMENT *theElement, COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM]);
INT HexaVolume (COORD **theCorners, COORD *volume);
#endif



#endif
