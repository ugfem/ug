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
										case 0 : r=t-1; break;\
										case 1 : r=1-t; break;\
										case 2 : r=t;   break;\
										case 3 : r=-t;  break;\
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
										case 0 : r=s-1; break;\
										case 1 : r=-s;  break;\
										case 2 : r=s;   break;\
										case 3 : r=1-s; break;\
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

#define LOCAL_TO_GLOBAL_QUADRILATERAL(x,local,global)         \
{(global)[0] = (1.0-(local)[0])*(1.0-(local)[1])*(x)[0][0]    \
   + (local)[0]*(1.0-(local)[1])*(x)[1][0]                    \
   + (local)[0]*(local)[1]*(x)[2][0]                          \
   + (1.0-(local)[0])*(local)[1]*(x)[3][0];                   \
 (global)[1] = (1.0-(local)[0])*(1.0-(local)[1])*(x)[0][1]    \
   + (local)[0]*(1.0-(local)[1])*(x)[1][1]                    \
   + (local)[0]*(local)[1]*(x)[2][1]                          \
   + (1.0-(local)[0])*(local)[1]*(x)[3][1];}

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

#define TRANSFORMATION_OF_TRIANGLE(x,M)     \
    { V2_SUBTRACT((x)[1],(x)[0],(M)[0]);    \
	  V2_SUBTRACT((x)[2],(x)[0],(M)[1]); }

#define TRANSFORMATION_OF_QUADRILATERAL(x,local,M)                        \
	{ DOUBLE a;                                                           \
	  a = 1.0 - (local)[1];                                               \
	  (M)[0][0] = a*((x)[1][0]-x[0][0])+(local)[1]*((x)[2][0]-(x)[3][0]); \
	  (M)[0][1] = a*((x)[1][1]-x[0][1])+(local)[1]*((x)[2][1]-(x)[3][1]); \
	  a = 1.0 - (local)[0];                                               \
	  (M)[1][0] = a*((x)[3][0]-x[0][0])+(local)[0]*((x)[2][0]-(x)[1][0]); \
	  (M)[1][1] = a*((x)[3][1]-x[0][1])+(local)[0]*((x)[2][1]-(x)[1][1]); }

#define TRANSFORMATION(n,x,local,M)                                      \
 {if ((n) == 3)      {TRANSFORMATION_OF_TRIANGLE((x),(M));}              \
  else if ((n) == 4) {TRANSFORMATION_OF_QUADRILATERAL((x),(local),(M));} }

#define SIDE_NORMAL(n,i,x,normal)                 \
   { DOUBLE s; DOUBLE_VECTOR y;                   \
   	 V2_SUBTRACT(x[(i+1)%n],x[i],y);              \
 	 V2_EUKLIDNORM(y,s);                          \
	 V2_SCALE(1.0/s,y);                           \
	 V2_SUBTRACT(x[(i+1)%n],x[(i+2)%n],normal);   \
	 V2_SCALAR_PRODUCT(normal,y,s);               \
	 V2_LINCOMB(1.0,normal,-s,y,normal);          \
	 V2_EUKLIDNORM(normal,s);                     \
	 V2_SCALE(1.0/s,normal);}         

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

#define CORNER_COORDINATES_PRISM(e,n,x)                                    \
                                  {(n) = 6;                                \
								   (x)[0]=CVECT(MYVERTEX(CORNER((e),0)));  \
								   (x)[1]=CVECT(MYVERTEX(CORNER((e),1)));  \
								   (x)[2]=CVECT(MYVERTEX(CORNER((e),2)));  \
                                   (x)[3]=CVECT(MYVERTEX(CORNER((e),3)));  \
                                   (x)[4]=CVECT(MYVERTEX(CORNER((e),3)));  \
                                   (x)[5]=CVECT(MYVERTEX(CORNER((e),4)));}

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
   else if (TAG((e))==PRISM)      CORNER_COORDINATES_PRISM((e),(n),(x))      \
   else if (TAG((e))==HEXAHEDRON) CORNER_COORDINATES_HEXAHEDRON((e),(n),(x))}
 
#define LOCAL_TO_GLOBAL_TETRAHEDRON(x,local,global)                       \
 {(global)[0] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][0]          \
   +(local)[0]*(x)[1][0] + (local)[1]*(x)[2][0] + (local)[2]*(x)[3][0];   \
  (global)[1] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][1]          \
   +(local)[0]*(x)[1][1] + (local)[1]*(x)[2][1] + (local)[2]*(x)[3][1];   \
  (global)[2] = (1.0-(local)[0]-(local)[1]-(local)[2])*(x)[0][2]          \
   +(local)[0]*(x)[1][2] + (local)[1]*(x)[2][2] + (local)[2]*(x)[3][2]; }  

#define LOCAL_TO_GLOBAL_PYRAMID(x,local,global)                             \
 {DOUBLE a,b,a0,a1,a2,a3;                                                   \
  a = 1.0 - (local)[0];                                                     \
  b = 1.0 - (local)[1];                                                     \
  if ((local)[0] > (local)[1]) {                                            \
  a0 = a * b - (local)[2] * b;                                              \
  a1 = (local)[0] * b - (local)[2]*(local)[1];                              \
  a2 = (local)[0] * (local)[1] + (local)[2]*(local)[1];                     \
  a3 = a * (local)[1] - (local)[2] * (local)[1];                            \
  (global)[0] =                                                               \
	a0*(x)[0][0]+a1*(x)[1][0]+a2*(x)[2][0]+a3*(x)[3][0]+(local)[2]*(x)[4][0]; \
  (global)[1] =                                                               \
	a0*(x)[0][1]+a1*(x)[1][1]+a2*(x)[2][1]+a3*(x)[3][1]+(local)[2]*(x)[4][1]; \
  (global)[2] =                                                               \
	a0*(x)[0][2]+a1*(x)[1][2]+a2*(x)[2][2]+a3*(x)[3][2]+(local)[2]*(x)[4][2];}\
  else {                                                                      \
  a0 = a * b - (local)[2] * a;                                                \
  a1 = (local)[0] * b - (local)[2]*(local)[0];                                \
  a2 = (local)[0] * (local)[1] + (local)[2]*(local)[0];                       \
  a3 = a * (local)[1] - (local)[2] * (local)[0];                              \
  (global)[0] =                                                               \
	a0*(x)[0][0]+a1*(x)[1][0]+a2*(x)[2][0]+a3*(x)[3][0]+(local)[2]*(x)[4][0]; \
  (global)[1] =                                                               \
	a0*(x)[0][1]+a1*(x)[1][1]+a2*(x)[2][1]+a3*(x)[3][1]+(local)[2]*(x)[4][1]; \
  (global)[2] =                                                               \
	a0*(x)[0][2]+a1*(x)[1][2]+a2*(x)[2][2]+a3*(x)[3][2]+(local)[2]*(x)[4][2];}}

#define LOCAL_TO_GLOBAL_PRISM(x,local,global)                               \
 {DOUBLE a,b,a0,a1,a2,a3,a4,a5,a6,a7;                                       \
  a = 1.0 - (local)[0] - (local)[1];                                        \
  b = 1.0 - (local)[2];                                                     \
  a0 = a * b;                                                               \
  a1 = (local)[0] * b;                                                      \
  a2 = (local)[1] * b;                                                      \
  a3 = a * (local)[2];                                                      \
  a4 = (local)[0] * (local)[2];                                             \
  a5 = (local)[1] * (local)[2];                                             \
  (global)[0] =                                                             \
	a0*(x)[0][0]+a1*(x)[1][0]+a2*(x)[2][0]+a3*(x)[3][0]+                    \
    a4*(x)[4][0]+a5*(x)[5][0];                                              \
  (global)[1] =                                                             \
	a0*(x)[0][1]+a1*(x)[1][1]+a2*(x)[2][1]+a3*(x)[3][1]+                    \
	a4*(x)[4][1]+a5*(x)[5][1];                                              \
  (global)[2] =                                                             \
	a0*(x)[0][2]+a1*(x)[1][2]+a2*(x)[2][2]+a3*(x)[3][2]+                    \
	a4*(x)[4][2]+a5*(x)[5][2]; }

#define LOCAL_TO_GLOBAL_HEXAHEDRON(x,local,global)                          \
 {DOUBLE a,b,c,a0,a1,a2,a3,a4,a5,a6,a7;                                     \
  a = 1.0 - (local)[0];                                                     \
  b = 1.0 - (local)[1];                                                     \
  c = 1.0 - (local)[2];                                                     \
  a0 = a * b * c;                                                           \
  a1 = (local)[0] * b * c;                                                  \
  a2 = (local)[0] * (local)[1] * c;                                         \
  a3 = a * (local)[1] * c;                                                  \
  a4 = a * b * (local)[2];                                                  \
  a5 = (local)[0] * b * (local)[2];                                         \
  a6 = (local)[0] * (local)[1] * (local)[2];                                \
  a7 = a * (local)[1] * (local)[2];                                         \
  (global)[0] =                                                             \
	a0*(x)[0][0]+a1*(x)[1][0]+a2*(x)[2][0]+a3*(x)[3][0]+                    \
	a4*(x)[4][0]+a5*(x)[5][0]+a6*(x)[6][0]+a7*(x)[7][0];                    \
  (global)[1] =                                                             \
	a0*(x)[0][1]+a1*(x)[1][1]+a2*(x)[2][1]+a3*(x)[3][1]+                    \
	a4*(x)[4][1]+a5*(x)[5][1]+a6*(x)[6][1]+a7*(x)[7][1];                    \
  (global)[2] =                                                             \
	a0*(x)[0][2]+a1*(x)[1][2]+a2*(x)[2][2]+a3*(x)[3][2]+                    \
	a4*(x)[4][2]+a5*(x)[5][2]+a6*(x)[6][2]+a7*(x)[7][2]; }

#define LOCAL_TO_GLOBAL(n,x,local,global)                                 \
 {if ((n) == 4)      LOCAL_TO_GLOBAL_TETRAHEDRON((x),(local),(global))    \
  else if ((n) == 5) LOCAL_TO_GLOBAL_PYRAMID((x),(local),(global))        \
  else if ((n) == 6) LOCAL_TO_GLOBAL_PRISM((x),(local),(global))          \
  else if ((n) == 8) LOCAL_TO_GLOBAL_HEXAHEDRON((x),(local),(global))}

#define AREA_OF_TETRAHEDRON(x,area)  {DOUBLE detJ; DOUBLE_VECTOR a,b,c;    \
									  V3_SUBTRACT((x)[1],(x)[0],a);        \
									  V3_SUBTRACT((x)[2],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[3],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  (area) = ABS(detJ)/6.0;}

#define AREA_OF_PYRAMID(x,area)      {DOUBLE detJ,ar; DOUBLE_VECTOR a,b,c,d;\
									  V3_SUBTRACT((x)[1],(x)[0],a);         \
									  V3_SUBTRACT((x)[2],(x)[0],b);         \
									  V3_VECTOR_PRODUCT(a,b,c);             \
									  V3_SUBTRACT((x)[4],(x)[0],d);         \
									  V3_SCALAR_PRODUCT(c,d,detJ);          \
									  ar = ABS(detJ)/6.0;                   \
									  V3_SUBTRACT((x)[3],(x)[0],a);         \
									  V3_VECTOR_PRODUCT(a,b,c);             \
									  V3_SCALAR_PRODUCT(c,d,detJ);          \
									  (area) = ABS(detJ)/6.0 + ar;}

#define AREA_OF_PRISM(x,area)   {DOUBLE detJ,ar; DOUBLE_VECTOR a,b,c;      \
									  V3_SUBTRACT((x)[1],(x)[0],a);        \
									  V3_SUBTRACT((x)[2],(x)[0],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[3],(x)[0],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  ar = ABS(detJ)/6.0;                  \
									  V3_SUBTRACT((x)[2],(x)[1],a);        \
									  V3_SUBTRACT((x)[3],(x)[1],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[4],(x)[1],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  ar += ABS(detJ)/6.0;                 \
									  V3_SUBTRACT((x)[2],(x)[5],a);        \
									  V3_SUBTRACT((x)[3],(x)[5],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[4],(x)[5],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
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
									  ar += ABS(detJ)/6.0;                 \
									  V3_SUBTRACT((x)[3],(x)[7],a);        \
									  V3_SUBTRACT((x)[4],(x)[7],b);        \
									  V3_VECTOR_PRODUCT(a,b,c);            \
									  V3_SUBTRACT((x)[6],(x)[7],a);        \
									  V3_SCALAR_PRODUCT(a,c,detJ);         \
									  (area) = ABS(detJ)/6.0 + ar;}

#define AREA_OF_ELEMENT(n,x,area)                        \
 {if ((n) == 4)      {AREA_OF_TETRAHEDRON((x),(area));}  \
  else if ((n) == 5) {AREA_OF_PYRAMID((x),(area));}      \
  else if ((n) == 6) {AREA_OF_PRISM((x),(area));}        \
  else if ((n) == 8) {AREA_OF_HEXAHEDRON((x),(area));}}

#define TRANSFORMATION_OF_TETRAHEDRON(x,M)   \
	{ V3_SUBTRACT((x)[1],(x)[0],(M)[0]);     \
	  V3_SUBTRACT((x)[2],(x)[0],(M)[1]);     \
	  V3_SUBTRACT((x)[3],(x)[0],(M)[2]);}

#define TRANSFORMATION_OF_PYRAMID(x,local,M)                              \
	{ DOUBLE a,b,c;                                                       \
	  a = (x)[0][0]-(x)[1][0]+(x)[2][0]-(x)[3][0];                        \
	  b = (x)[0][1]-(x)[1][1]+(x)[2][1]-(x)[3][1];                        \
	  c = (x)[0][2]-(x)[1][2]+(x)[2][2]-(x)[3][2];                        \
	  if ((local)[0] > (local)[1]) {                                      \
	  (M)[0][0] = (x)[1][0]-(x)[0][0]+(local)[1]*a;                       \
	  (M)[0][1] = (x)[1][1]-(x)[0][1]+(local)[1]*b;                       \
	  (M)[0][2] = (x)[1][2]-(x)[0][2]+(local)[1]*c;                       \
	  (M)[1][0] = (x)[3][0]-(x)[0][0]+((local)[0]+(local)[2])*a;          \
	  (M)[1][1] = (x)[3][1]-(x)[0][1]+((local)[0]+(local)[2])*b;          \
	  (M)[1][2] = (x)[3][2]-(x)[0][2]+((local)[0]+(local)[2])*c;          \
	  (M)[2][0] = (x)[4][0]-(x)[0][0]+(local)[1]*a;                       \
	  (M)[2][1] = (x)[4][1]-(x)[0][1]+(local)[1]*b;                       \
	  (M)[2][2] = (x)[4][2]-(x)[0][2]+(local)[1]*c;}                      \
	  else {                                                              \
	  (M)[0][0] = (x)[1][0]-(x)[0][0]+((local)[1]+(local)[2])*a;          \
	  (M)[0][1] = (x)[1][1]-(x)[0][1]+((local)[1]+(local)[2])*b;          \
	  (M)[0][2] = (x)[1][2]-(x)[0][2]+((local)[1]+(local)[2])*c;          \
	  (M)[1][0] = (x)[3][0]-(x)[0][0]+(local)[0]*a;                       \
	  (M)[1][1] = (x)[3][1]-(x)[0][1]+(local)[0]*b;                       \
	  (M)[1][2] = (x)[3][2]-(x)[0][2]+(local)[0]*c;                       \
	  (M)[2][0] = (x)[4][0]-(x)[0][0]+(local)[0]*a;                       \
	  (M)[2][1] = (x)[4][1]-(x)[0][1]+(local)[0]*b;                       \
	  (M)[2][2] = (x)[4][2]-(x)[0][2]+(local)[0]*c;}  }

#define TRANSFORMATION_OF_PRISM(x,local,M)                                \
	{ DOUBLE a0,a1,a2,b0,b1,b2;                                           \
	  a0 = (x)[0][0]-(x)[1][0]+(x)[3][0]-(x)[4][0];                       \
	  a1 = (x)[0][1]-(x)[1][1]+(x)[3][1]-(x)[4][1];                       \
	  a2 = (x)[0][2]-(x)[1][2]+(x)[3][2]-(x)[4][2];                       \
	  b0 = (x)[0][0]-(x)[2][0]+(x)[3][0]-(x)[5][0];                       \
	  b1 = (x)[0][1]-(x)[2][1]+(x)[3][1]-(x)[5][1];                       \
	  b2 = (x)[0][2]-(x)[2][2]+(x)[3][2]-(x)[5][2];                       \
	  (M)[0][0] = (x)[1][0]-(x)[0][0]+(local)[2]*a0;                      \
	  (M)[0][1] = (x)[1][1]-(x)[0][1]+(local)[2]*a1;                      \
	  (M)[0][2] = (x)[1][2]-(x)[0][2]+(local)[2]*a2;                      \
	  (M)[1][0] = (x)[3][0]-(x)[0][0]+(local)[2]*b0;                      \
	  (M)[1][1] = (x)[3][1]-(x)[0][1]+(local)[2]*b1;                      \
	  (M)[1][2] = (x)[3][2]-(x)[0][2]+(local)[2]*b2;                      \
	  (M)[2][0] = (x)[4][0]-(x)[0][0]+(local)[0]*a0+(local)[1]*b0;        \
	  (M)[2][1] = (x)[4][1]-(x)[0][1]+(local)[0]*a1+(local)[1]*b1;        \
	  (M)[2][2] = (x)[4][2]-(x)[0][2]+(local)[0]*a2+(local)[1]*b2;}  

#define TRANSFORMATION_OF_HEXAHEDRON(x,local,M)                           \
	{ DOUBLE a,b,c,a0,a1,a2,a3;                                           \
	  a = 1.0 - (local)[0];                                               \
	  b = 1.0 - (local)[1];                                               \
	  c = 1.0 - (local)[2];                                               \
      a0 = b * c;                                                         \
      a1 = (local)[1] * c;                                                \
      a2 = (local)[1] * (local)[2];                                       \
      a3 = b * (local)[2];                                                \
	  (M)[0][0] = a0*((x)[1][0]-x[0][0])+a1*((x)[2][0]-(x)[3][0])         \
	            + a2*((x)[6][0]-x[7][0])+a3*((x)[5][0]-(x)[4][0]);        \
	  (M)[0][1] = a0*((x)[1][1]-x[0][1])+a1*((x)[2][1]-(x)[3][1])         \
	            + a2*((x)[6][1]-x[7][1])+a3*((x)[5][1]-(x)[4][1]);        \
	  (M)[0][2] = a0*((x)[1][2]-x[0][2])+a1*((x)[2][2]-(x)[3][2])         \
	            + a2*((x)[6][2]-x[7][2])+a3*((x)[5][2]-(x)[4][2]);        \
      a0 = a * c;                                                         \
      a1 = (local)[0] * c;                                                \
      a2 = (local)[0] * (local)[2];                                       \
      a3 = a * (local)[2];                                                \
	  (M)[1][0] = a0*((x)[3][0]-x[0][0])+a1*((x)[2][0]-(x)[1][0])         \
	            + a2*((x)[6][0]-x[5][0])+a3*((x)[7][0]-(x)[4][0]);        \
	  (M)[1][1] = a0*((x)[3][1]-x[0][1])+a1*((x)[2][1]-(x)[1][1])         \
	            + a2*((x)[6][1]-x[5][1])+a3*((x)[7][1]-(x)[4][1]);        \
	  (M)[1][2] = a0*((x)[3][2]-x[0][2])+a1*((x)[2][2]-(x)[1][2])         \
	            + a2*((x)[6][2]-x[5][2])+a3*((x)[7][2]-(x)[4][2]);        \
      a0 = a * b;                                                         \
      a1 = (local)[0] * b;                                                \
      a2 = (local)[0] * (local)[1];                                       \
      a3 = a * (local)[1];                                                \
	  (M)[2][0] = a0*((x)[4][0]-x[0][0])+a1*((x)[5][0]-(x)[1][0])         \
	            + a2*((x)[6][0]-x[2][0])+a3*((x)[7][0]-(x)[3][0]);        \
	  (M)[2][1] = a0*((x)[4][1]-x[0][1])+a1*((x)[5][1]-(x)[1][1])         \
	            + a2*((x)[6][1]-x[2][1])+a3*((x)[7][1]-(x)[3][1]);        \
	  (M)[2][2] = a0*((x)[4][2]-x[0][2])+a1*((x)[5][2]-(x)[1][2])         \
	            + a2*((x)[6][2]-x[2][2])+a3*((x)[7][2]-(x)[3][2]); }

#define TRANSFORMATION(n,x,local,M)                                     \
 {if ((n) == 4)      {TRANSFORMATION_OF_TETRAHEDRON((x),(M));}          \
  else if ((n) == 5) {TRANSFORMATION_OF_PYRAMID((x),(local),(M));}      \
  else if ((n) == 6) {TRANSFORMATION_OF_PRISM((x),(local),(M));}        \
  else if ((n) == 8) {TRANSFORMATION_OF_HEXAHEDRON((x),(local),(M));}}

#define SIDE_NORMAL(n,i,x,normal)                                \
  { DOUBLE s; DOUBLE_VECTOR a,b; ELEMENT e;                      \
	V3_SUBTRACT(x[CORNER_OF_SIDE_REF((n),(i),1)],                \
                x[CORNER_OF_SIDE_REF((n),(i),0)],a);             \
	V3_SUBTRACT(x[CORNER_OF_SIDE_REF((n),(i),2)],                \
				x[CORNER_OF_SIDE_REF((n),(i),0)],b);             \
	V3_VECTOR_PRODUCT(a,b,(normal));                             \
	V3_EUKLIDNORM((normal),s);                                   \
	V3_SCALE(1.0/s,(normal));}

#endif

#ifdef __TWODIM__
#define GlobalToLocal(n,c,e,l)		 GlobalToLocal2d (n,c,e,l)
#define LocalToGlobal(n,c,e,l)		 LocalToGlobal2d (n,c,e,l)
#endif
#ifdef __THREEDIM__
#define GlobalToLocal(n,c,e,l)		 GlobalToLocal3d (n,c,e,l)
#define LocalToGlobal(n,c,e,l)		 LocalToGlobal3d (n,c,e,l)
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

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

DOUBLE 		GN 		(INT n, INT i, COORD *ip_local);
INT 		GNs 	(INT n, COORD *ip_local, DOUBLE *result);
INT 		D_GN 	(INT n, INT i, COORD *ip_local, DOUBLE *derivative);
COORD 	   *LMP 	(INT n);

INT LocalCornerCoordinates (INT dim, INT tag, INT corner, DOUBLE *result);
/*****************************************************************************
 * This function delivers the coordinates of the corners in the refernence
 * element. 
 * Parameters are:
 * dim          space dimension, 2D can be called in 3D
 * tag 			identifies element type, i.e. 3=TRIANGLE, 4=QUADRILATERAL here
 * corner		number of the corner
 * result		array to place result
 */

INT InterpolateFEFunction (INT dim, INT tag, DOUBLE ip_local[DIM], 
		DOUBLE nodal_values[MAX_CORNERS_OF_ELEM], DOUBLE *result);
/*****************************************************************************
 * This function interpolates  a finite element function given by
 * nodal values at some point in the interior of the element. The coordinates
 * of this point must be given in local coordinates.
 * Parameters are:
 * dim          space dimension, 2D can be called in 3D
 * tag 			identifies element type, i.e. 3=TRIANGLE, 4=QUADRILATERAL here
 * ip_local		local coordinates of interpolation point
 * nodal_values array of nodal values of function to interpolate
 * result		pointer where to place result (one DOUBLE value)
 */


INT LinearTrafo (INT dim, INT tag);
/*****************************************************************************
 * This function returns true if transformation from reference element
 * to an arbitrary element is linear. This can be used to avoid
 * recomputation of the Jacobian of the transformation.
 * Parameters are:
 * dim          space dimension, 2D can be called in 3D
 * tag 			identifies element type, i.e. 3=TRIANGLE, 4=QUADRILATERAL here
 */


INT JacobianInverse (INT dim, INT tag, DOUBLE co_global[MAX_CORNERS_OF_ELEM][DIM],
		DOUBLE ip_local[DIM], DOUBLE Jinv[DIM][DIM], DOUBLE *detJ);
/*****************************************************************************
 * Compute inverse of the jacobian of transformation of reference element
 * to some element given by global coordinates of corners. The determinant
 * of the jacobian (no its inverse!) is also provided as a result.
 * Parameters are:
 * dim          space dimension, 2D can be called in 3D
 * tag 			identifies element type, i.e. 3=TRIANGLE, 4=QUADRILATERAL here
 * co_global	global coordinates of the corners of the element
 * ip_local		local coordinates of interpolation point
 * Jinv			place to store the inverse of jacobian
 * detJ			place to store determinant
 */


INT GradientFEFunction (INT dim, INT tag, DOUBLE ip_local[DIM], DOUBLE Jinv[DIM][DIM], 
		DOUBLE nodal_values[MAX_CORNERS_OF_ELEM], DOUBLE result[DIM]);
/*****************************************************************************
 * Compute gradient in global coordinates of some finite element function
 * given by nodal values at some point within the element.
 * Parameters are:
 * dim          space dimension, 2D can be called in 3D
 * tag 			identifies element type, i.e. 3=TRIANGLE, 4=QUADRILATERAL here
 * ip_local		local coordinates of interpolation point
 * Jinv			inverse of jacobian of transformation computed by fct above
 * nodal_values array of nodal values of function to interpolate
 * result		pointer where to place result (a vector)
 */




#ifdef __TWODIM__
DOUBLE	N				(INT n, INT i, DOUBLE s, DOUBLE t);
DOUBLE	dNds			(INT n, INT i, DOUBLE s, DOUBLE t);
DOUBLE	dNdt			(INT n, INT i, DOUBLE s, DOUBLE t);

INT 	Derivatives 	(INT n, const DOUBLE *px, const DOUBLE *py, DOUBLE ips, DOUBLE ipt, DOUBLE *dNdx, DOUBLE *dNdy, DOUBLE *detJ);
INT		Gradients		(INT n, const COORD **theCorners, DOUBLE ips, DOUBLE ipt, DOUBLE_VECTOR Gradient[MAX_CORNERS_OF_ELEM], DOUBLE *DetJ);

INT		L2GDerivative2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD *Derivative);
INT		LocalToGlobal2d	(INT n, const COORD **Corners, const COORD *EvalPoint, COORD *GlobalCoord);
INT 	GlobalToLocal2d (INT n, const COORD **Corners, const COORD *EvalPoint, COORD *LocalCoord);
#endif

#ifdef __THREEDIM__						 
INT		LocalToGlobal3d	(INT n, const COORD **Corners, const COORD *EvalPoint, COORD *GlobalCoord);
INT     GlobalToLocal3d     (INT n, const COORD **Corners, const COORD *EvalPoint, COORD *LocalCoord);
INT GetSkewedUIP (const COORD_VECTOR *theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], const DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
DOUBLE  N                   (const INT i, const COORD *LocalCoord);
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
#endif

#endif
