// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  shapes.h														*/
/*																			*/
/* Purpose:   header file for shape functions								*/
/*																			*/
/* Author:	  Klaus Johannsen 												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   28.11.95 begin, ug version 3.1								*/
/*																			*/
/* Remarks: 																*/
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
#define GlobalToLocal(n,c,e,l)		 GlobalToLocal2d (n,c,e,l)
#endif
#ifdef __THREEDIM__
#define GlobalToLocal(n,c,e,l)		 GlobalToLocal3d (c,e,l)
#endif

#ifdef __TWODIM__
#define DNDS(n,i,s,t,r) 		if (n==3)\
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

#define DNDT(n,i,s,t,r)		    if (n==3)\
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
#define CORNER_COORDINATES_TRIANGLE(e,n,x)                                \
                                   {(n)=3;                                \
								   (x)[0]=CVECT(MYVERTEX(CORNER((e),0))); \
								   (x)[1]=CVECT(MYVERTEX(CORNER((e),1))); \
								   (x)[2]=CVECT(MYVERTEX(CORNER((e),2)));}

#define CORNER_COORDINATES_QUADRILATERAL(e,n,x)                           \
                                   {(n)=4;                                \
								   (x)[0]=CVECT(MYVERTEX(CORNER((e),0))); \
								   (x)[1]=CVECT(MYVERTEX(CORNER((e),1))); \
								   (x)[2]=CVECT(MYVERTEX(CORNER((e),2))); \
								   (x)[3]=CVECT(MYVERTEX(CORNER((e),3)));}

#define CORNER_COORDINATES(e,n,x)                               \
 {if (TAG((e))==TRIANGLE)                                       \
                 CORNER_COORDINATES_TRIANGLE((e),(n),(x))       \
  else if (TAG((e))==QUADRILATERAL)                             \
                 CORNER_COORDINATES_QUADRILATERAL((e),(n),(x))}

#define LOCAL_TO_GLOBAL_TRIANGLE(x,local,global)          \
 {(global)[0] = (1.0-(local)[0]-(local)[1])*(x)[0][0]     \
   +(local)[0]*(x)[1][0] + (local)[1]*(x)[2][0];          \
  (global)[1] = (1.0-(local)[0]-(local)[1])*(x)[0][1]     \
   +(local)[0]*(x)[1][1] + (local)[1]*(x)[2][1]; }  

#define LOCAL_TO_GLOBAL_QUADRILATERAL(x,local,global)                   \
 {(global)[0] = 0.25 * ((1.0-(local)[0])*(1.0-(local)[1])*(x)[1][0]     \
                      + (1.0+(local)[0])*(1.0-(local)[1])*(x)[2][0]     \
                      + (1.0+(local)[0])*(1.0+(local)[1])*(x)[3][0]     \
                      + (1.0-(local)[0])*(1.0+(local)[1])*(x)[4][0]);   \
  (global)[1] = 0.25 * ((1.0-(local)[0])*(1.0-(local)[1])*(x)[1][1]     \
                      + (1.0+(local)[0])*(1.0-(local)[1])*(x)[2][1]     \
                      + (1.0+(local)[0])*(1.0+(local)[1])*(x)[3][1]     \
                      + (1.0-(local)[0])*(1.0+(local)[1])*(x)[4][1]);}

#define LOCAL_TO_GLOBAL(n,x,local,global)                                 \
 {if ((n) == 3)      LOCAL_TO_GLOBAL_TRIANGLE((x),(local),(global))       \
  else if ((n) == 4) LOCAL_TO_GLOBAL_QUADRILATERAL((x),(local),(global))}

#define AREA_OF_TRIANGLE(x,area)     {DOUBLE detJ; DOUBLE_VECTOR a,b;      \
									  V2_SUBTRACT((x)[1],(x)[0],a);        \
									  V2_SUBTRACT((x)[2],(x)[0],b);        \
									  V2_VECTOR_PRODUCT(a,b,detJ);         \
									  (area) = ABS(detJ) * 0.5;}

#define AREA_OF_QUADRILATERAL(x,area)   {DOUBLE detJ,ar; DOUBLE_VECTOR a,b;  \
										 V2_SUBTRACT((x)[1],(x)[0],a);       \
										 V2_SUBTRACT((x)[2],(x)[0],b);       \
										 V2_VECTOR_PRODUCT(a,b,detJ);        \
									     ar = ABS(detJ) * 0.5;               \
										 V2_SUBTRACT((x)[3],(x)[0],a);       \
										 V2_VECTOR_PRODUCT(a,b,detJ);        \
									     (area) = ABS(detJ) * 0.5 + ar;}

#define AREA_OF_ELEMENT(n,x,area)                       \
 {if ((n) == 3)      AREA_OF_TRIANGLE((x),(area))       \
  else if ((n) == 4) AREA_OF_QUADRILATERAL((x),(area)) }

#endif /* __TWODIM__ */

#ifdef __THREEDIM__
#define CORNER_COORDINATES_TETRAHEDRON(e,n,x)                              \
                                  {(n) = 4;                                \
								   (x)[0]=CVECT(MYVERTEX(CORNER((e),0)));  \
								   (x)[1]=CVECT(MYVERTEX(CORNER((e),1)));  \
								   (x)[2]=CVECT(MYVERTEX(CORNER((e),2)));  \
				                   (x)[3]=CVECT(MYVERTEX(CORNER((e),3)));}

#define CORNER_COORDINATES_PYRAMID(e,n,x)                                  \
                                  {(n) = 5;                                \
								   (x)[0]=CVECT(MYVERTEX(CORNER((e),0)));  \
								   (x)[1]=CVECT(MYVERTEX(CORNER((e),1)));  \
								   (x)[2]=CVECT(MYVERTEX(CORNER((e),2)));  \
                                   (x)[3]=CVECT(MYVERTEX(CORNER((e),3)));  \
                                   (x)[4]=CVECT(MYVERTEX(CORNER((e),4)));}

#define CORNER_COORDINATES_HEXAHEDRON(e,n,x)                               \
                                  {(n) = 8;                                \
								   (x)[0]=CVECT(MYVERTEX(CORNER((e),0)));  \
								   (x)[1]=CVECT(MYVERTEX(CORNER((e),1)));  \
								   (x)[2]=CVECT(MYVERTEX(CORNER((e),2)));  \
				                   (x)[3]=CVECT(MYVERTEX(CORNER((e),3)));  \
								   (x)[4]=CVECT(MYVERTEX(CORNER((e),4)));  \
								   (x)[5]=CVECT(MYVERTEX(CORNER((e),5)));  \
								   (x)[6]=CVECT(MYVERTEX(CORNER((e),6)));  \
								   (x)[7]=CVECT(MYVERTEX(CORNER((e),7)));}

#define CORNER_COORDINATES(e,n,x)                                            \
  {if (TAG((e))==TETRAHEDRON)     CORNER_COORDINATES_TETRAHEDRON((e),(n),(x))\
   else if (TAG((e))==PYRAMID)    CORNER_COORDINATES_PYRAMID((e),(n),(x))    \
   else if (TAG((e))==HEXAHEDRON) CORNER_COORDINATES_HEXAHEDRON((e),(n),(x))}
 
#define LOCAL_TO_GLOBAL_TETRAHEDRON(x,local,global)                       \
 {(global)[0] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][0]          \
   +(local)[0]*(x)[1][0] + (local)[1]*(x)[2][0] + (local)[2]*(x)[3][0];   \
  (global)[1] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][1]          \
   +(local)[0]*(x)[1][1] + (local)[1]*(x)[2][1] + (local)[2]*(x)[3][1];   \
  (global)[2] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][2]          \
   +(local)[0]*(x)[1][2] + (local)[1]*(x)[2][2] + (local)[2]*(x)[3][2]; }  

#define LOCAL_TO_GLOBAL_PYRAMID(x,local,global)                                 \
 {(global)[0] = 0.25 * ((1.0-(local)[0])*(1.0-(local)[1])*(1.0-(local)[2])*(x)[0][0]     \
                  + (local)[0]*(1.0-(local)[1])*(1.0-(local)[2])*(x)[1][0]          \
                  + (local)[0]*(local)[1]*(1.0-(local)[2])*(x)[2][0]                \
                  + (1.0-(local)[0])*(local)[1]*(1.0-(local)[2])*(x)[3][0]          \
                  + (1.0-(local)[0]-(local)[1])*(local)[2]*(x)[4][0]);              \
  (global)[1] = 0.25 * ((1.0-(local)[0])*(1.0-(local)[1])*(1.0-(local)[2])*(x)[0][1]     \
                  + (local)[0]*(1.0-(local)[1])*(1.0-(local)[2])*(x)[1][1]          \
                  + (local)[0]*(local)[1]*(1.0-(local)[2])*(x)[2][1]                \
                  + (1.0-(local)[0])*(local)[1]*(1.0-(local)[2])*(x)[3][1]          \
                  + (1.0-(local)[0]-(local)[1])*(local)[2]*(x)[4][1]);              \
  (global)[2] = 0.25 * ((1.0-(local)[0])*(1.0-(local)[1])*(1.0-(local)[2])*(x)[0][2]     \
                  + (local)[0]*(1.0-(local)[1])*(1.0-(local)[2])*(x)[1][2]          \
                  + (local)[0]*(local)[1]*(1.0-(local)[2])*(x)[2][2]                \
                  + (1.0-(local)[0])*(local)[1]*(1.0-(local)[2])*(x)[3][2]          \
                  + (1.0-(local)[0]-(local)[1])*(local)[2]*(x)[4][2]);}
#define LOCAL_TO_GLOBAL_HEXAHEDRON(x,local,global)                               \
 {(global)[0] = 0.125 * ((1.0-(local)[0])*(1.0-(local)[1])*(1.0-(local)[2])*(x)[0][0]     \
                  + (local)[0]*(1.0-(local)[1])*(1.0-(local)[2])*(x)[1][0]           \
                  + (local)[0]*(local)[1]*(1.0-(local)[2])*(x)[2][0]                 \
                  + (1.0-(local)[0])*(local)[1]*(1.0-(local)[2])*(x)[3][0]           \
                  + (1.0-(local)[0])*(1.0-(local)[1])*(local)[2]*(x)[4][0]           \
                  + (local)[0]*(1.0-(local)[1])*(local)[2]*(x)[5][0]                 \
                  + (local)[0]*(local)[1]*(local)[2]*(x)[6][0]                       \
                  + (1.0-(local)[0])*(local)[1]*(local)[2]*(x)[7][0]);               \
  (global)[1] = 0.125 * ((1.0-(local)[0])*(1.0-(local)[1])*(1.0-(local)[2])*(x)[0][1]     \
                  + (local)[0]*(1.0-(local)[1])*(1.0-(local)[2])*(x)[1][1]           \
                  + (local)[0]*(local)[1]*(1.0-(local)[2])*(x)[2][1]                 \
                  + (1.0-(local)[0])*(local)[1]*(1.0-(local)[2])*(x)[3][1]           \
                  + (1.0-(local)[0])*(1.0-(local)[1])*(local)[2]*(x)[4][1]           \
                  + (local)[0]*(1.0-(local)[1])*(local)[2]*(x)[5][1]                 \
                  + (local)[0]*(local)[1]*(local)[2]*(x)[6][1]                       \
                  + (1.0-(local)[0])*(local)[1]*(local)[2]*(x)[7][1]);               \
  (global)[2] = 0.125 * ((1.0-(local)[0])*(1.0-(local)[1])*(1.0-(local)[2])*(x)[0][2]     \
                  + (local)[0]*(1.0-(local)[1])*(1.0-(local)[2])*(x)[1][2]           \
                  + (local)[0]*(local)[1]*(1.0-(local)[2])*(x)[2][2]                 \
                  + (1.0-(local)[0])*(local)[1]*(1.0-(local)[2])*(x)[3][2]           \
                  + (1.0-(local)[0])*(1.0-(local)[1])*(local)[2]*(x)[4][2]           \
                  + (local)[0]*(1.0-(local)[1])*(local)[2]*(x)[5][2]                 \
                  + (local)[0]*(local)[1]*(local)[2]*(x)[6][2]                       \
                  + (1.0-(local)[0])*(local)[1]*(local)[2]*(x)[7][2]);}

#define LOCAL_TO_GLOBAL(n,x,local,global)                                 \
 {if ((n) == 4)      LOCAL_TO_GLOBAL_TETRAHEDRON((x),(local),(global))    \
  else if ((n) == 5) LOCAL_TO_GLOBAL_PYRAMID((x),(local),(global))        \
  else if ((n) == 8) LOCAL_TO_GLOBAL_HEXAHEDRON((x),(local),(global))}

#define AREA_OF_TETRAHEDRON(x,area)  {DOUBLE detJ; DOUBLE_VECTOR a,b,c;    \
									  V3_SUBTRACT((x)[1],(x)[0],a);        \
									  V3_SUBTRACT((x)[2],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[3],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  (area) = ABS(detJ)/6.0;}

#define AREA_OF_PYRAMID(x,area)      {DOUBLE detJ,ar; DOUBLE_VECTOR a,b,c; \
									  V3_SUBTRACT((x)[1],(x)[0],a);        \
									  V3_SUBTRACT((x)[2],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[4],(x)[0],b);        \
									  V3_SCALAR_PRODUCT(b,c,detJ);         \
									  ar = ABS(detJ)/6.0;                  \
									  V3_SUBTRACT((x)[3],(x)[0],a);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SCALAR_PRODUCT(b,c,detJ);         \
									  (area) = ABS(detJ)/6.0 + ar;}

#define AREA_OF_HEXAHEDRON(x,area)   {DOUBLE detJ,ar; DOUBLE_VECTOR a,b,c; \
									  V3_SUBTRACT((x)[1],(x)[0],a);        \
									  V3_SUBTRACT((x)[2],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[5],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  ar = ABS(detJ)/6.0;                  \
									  V3_SUBTRACT((x)[2],(x)[0],a);        \
									  V3_SUBTRACT((x)[5],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[6],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  ar += ABS(detJ)/6.0;                 \
									  V3_SUBTRACT((x)[4],(x)[0],a);        \
									  V3_SUBTRACT((x)[5],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[6],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  ar += ABS(detJ)/6.0;                 \
									  V3_SUBTRACT((x)[2],(x)[0],a);        \
									  V3_SUBTRACT((x)[3],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[6],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  ar += ABS(detJ)/6.0;                 \
									  V3_SUBTRACT((x)[3],(x)[0],a);        \
									  V3_SUBTRACT((x)[4],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[6],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  ar+= ABS(detJ)/6.0;                  \
									  V3_SUBTRACT((x)[3],(x)[7],a);        \
									  V3_SUBTRACT((x)[4],(x)[7],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[6],(x)[7],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  (area) = ABS(detJ)/6.0 + ar;}

#define AREA_OF_ELEMENT(n,x,area)                        \
 {if ((n) == 4)      {AREA_OF_TETRAHEDRON((x),(area));}  \
  else if ((n) == 5) {AREA_OF_PYRAMID((x),(area));}      \
  else if ((n) == 8) {AREA_OF_HEXAHEDRON((x),(area));}}

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
extern	COORD	HexRefCoord[MAX_CORNERS_OF_ELEM][DIM];
extern  COORD   TransfCoeff[MAX_CORNERS_OF_ELEM][DIM];
extern  COORD   CenterOfIntergrSurf[MAX_EDGES_OF_ELEM][DIM];
#endif

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

DOUBLE 		GN 		(INT n, INT i, COORD *local);
COORD 	   *LMP 	(INT n);


#ifdef __TWODIM__
DOUBLE	N				(INT n, INT i, DOUBLE s, DOUBLE t);
DOUBLE	dNds			(INT n, INT i, DOUBLE s, DOUBLE t);
DOUBLE	dNdt			(INT n, INT i, DOUBLE s, DOUBLE t);

INT 	Derivatives 	(INT n, const DOUBLE *px, const DOUBLE *py, DOUBLE ips, DOUBLE ipt, DOUBLE *dNdx, DOUBLE *dNdy, DOUBLE *detJ);
INT		Gradients		(INT n, const COORD **theCorners, DOUBLE ips, DOUBLE ipt, DOUBLE_VECTOR Gradient[MAX_CORNERS_OF_ELEM], DOUBLE *DetJ);

INT		L2GDerivative2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD *Derivative);
INT		LocalToGlobal2d	(INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR GlobalCoord);
INT 	GlobalToLocal2d (INT n,const  COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);
INT		specialGlobalToLocal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);
#endif

#ifdef __THREEDIM__
INT TransformGlobalToLocal3D(ELEMENT *theElement, COORD_VECTOR Global, COORD_VECTOR Local);
DOUBLE  N                   (const INT i, const COORD *LocalCoord);
INT     GlobalToLocal3d     (const COORD **Corners, const COORD *EvalPoint, COORD *LocalCoord);
INT     TetraDerivative     (ELEMENT *theElement, const COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM]);
INT     TetraVolume         (const COORD **theCorners, COORD *volume);
INT     FV_TetInfo          (const COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM]);
INT 	Side_TetInfo		(COORD **theCorners, INT side, COORD_VECTOR Area, COORD_VECTOR GIP[3]);
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
