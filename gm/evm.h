// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  evm.h 														*/
/*																			*/
/* Purpose:   elementary vector manipulations, header for evm.c 			*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Institut fuer Computeranwendungen 							*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de 					    	*/
/*																			*/
/* History:   8.12.94 begin, ug3-version									*/
/*																			*/
/* Remarks: 																*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __EVM__
#define __EVM__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __GM__
#include "gm.h"
#endif

#include "debug.h"

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
/* macros																	*/
/*																			*/
/****************************************************************************/

/* space dimension indices */
#define _X_		0
#define _Y_		1
#define _Z_		2

/* misc macros */
#define SQRT(a)							sqrt((double)(a))
#define POW(a,b)						pow((double)(a),(double)(b))
#define ISNaN(x)						(!((x)-(x)==0) || ((x)-(x)!=0))

/* macros for coord points */
#define COPY_SC_TO_SH(p1,p2)			(p2).x=(short)((p1).x);(p2).y=(short)((p1).y)
#define CP_SUBTRACT(A,B,C)			   {(C).x = (A).x - (B).x;\
										(C).y = (A).y - (B).y;}
#define CP_LIMCOMB(a,A,b,B,C)		   {(C).x = (COORD)(a)*(A).x + (COORD)(b)*(B).x;\
										(C).y = (COORD)(a)*(A).y + (COORD)(b)*(B).y;}
#define CP_SCALARPRODUCT(A,B,c) 		(c) = (A).x*(B).x + (A).y*(B).y;
#define CP_EUKLIDNORM(A,b)				(b) = (COORD)sqrt((double)((A).x*(A).x+(A).y*(A).y));

/* macros for 2D vector operations */
#define V2_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
										(C)[1] = (a)*(A)[1] + (b)*(B)[1];}
#define V2_COPY(A,C)				   {(C)[0] = (A)[0];\
										(C)[1] = (A)[1];}
#define V2_SUBTRACT(A,B,C)			   {(C)[0] = (A)[0] - (B)[0];\
										(C)[1] = (A)[1] - (B)[1];}
#define V2_ADD(A,B,C)				   {(C)[0] = (A)[0] + (B)[0];\
										(C)[1] = (A)[1] + (B)[1];}
#define V2_ADD1(A,C)				   {(C)[0] += (A)[0];\
										(C)[1] += (A)[1];}
#define V2_AVG2(A,B,C)				   {(C)[0] = 0.5*((A)[0] + (B)[0]);\
	                                    (C)[1] = 0.5*((A)[1] + (B)[1]);}
#define V2_AVG4(A,B,C,D,E)			   {(E)[0] = 0.25*((A)[0]+(B)[0]+(C)[0]+(D)[0]);\
                                     	(E)[1] = 0.25*((A)[1]+(B)[1]+(C)[1]+(D)[1]);}
#define V2_SCALE(c,C)				   {(C)[0] = (c)*(C)[0];\
										(C)[1] = (c)*(C)[1];}
#define V2_NORMAL(A,B)				   {(A)[0] =  (B)[1];\
										(A)[1] = -(B)[0];}
#define V2_SCALEADD1(c,A,C)			   {(C)[0] += (c)*(A)[0];\
										(C)[1] += (c)*(A)[1];}
#define V2_SCALESET(c,A,C)			   {(C)[0] = (c)*(A)[0];\
										(C)[1] = (c)*(A)[1];}
#define V2_VECTOR_PRODUCT(A,B,c)		(c) = (A)[0]*(B)[1] - (A)[1]*(B)[0];
#define V2_COMPARE(A,B,c)				(c) = (((ABS((A)[0]-(B)[0])<SMALL_C)&&(ABS((A)[1]-(B)[1])<SMALL_C))?(0):(1));
#define V2_ISEQUAL(A,B) 				((ABS((A)[0]-(B)[0])<SMALL_C)&&(ABS((A)[1]-(B)[1])<SMALL_C))
#define V2_EUKLIDNORM(A,b)				(b) = sqrt((double)((A)[0]*(A)[0]+(A)[1]*(A)[1]));
#define V2_EUKLIDNORM_OF_DIFF(A,B,b)	(b) = sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])));
#define V2_CLEAR(A) 				   {(A)[0] = 0.0; (A)[1]= 0.0;}
#define V2_SCALAR_PRODUCT(A,B,c)		(c) = (A)[0]*(B)[0]+(A)[1]*(B)[1];
#define V2_ISZERO(A)					(ABS((A)[0])<SMALL_C && ABS((A)[1])<SMALL_C)
#define V2_SUP(v,s)                    {s = MAX(ABS(v[0]),ABS(v[1]));}


/* macros for 2D matrix-vector operations */
#define M2_TIMES_V2(M,A,B)			   {(B)[0] = (M)[0]*(A)[0] + (M)[2]*(A)[1];\
										(B)[1] = (M)[1]*(A)[0] + (M)[3]*(A)[1];}
#define MM2_TIMES_V2(M,A,B)			   {(B)[0] = (M)[0][0]*(A)[0] + (M)[0][1]*(A)[1];\
										(B)[1] = (M)[1][0]*(A)[0] + (M)[1][1]*(A)[1];}

#define MT2_TIMES_V2(M,A,B)			   {(B)[0] = (M)[0][0]*(A)[0] + (M)[1][0]*(A)[1];\
										(B)[1] = (M)[0][1]*(A)[0] + (M)[1][1]*(A)[1];}

#define MD2_TIMES_V2(M,A,B)			   {(B)[0] = (M)[0]*(A)[0];\
										(B)[1] = (M)[1]*(A)[1];}

/* macros for matrix operations */
#define M2_SCALE(c,M)				   {(M)[0] = (c)*(M)[0];\
										(M)[1] = (c)*(M)[1];\
										(M)[2] = (c)*(M)[2];\
										(M)[3] = (c)*(M)[3];}
#define M2_ADD(A,B,C)				   {(C)[0] = (A)[0]+(B)[0];\
										(C)[1] = (A)[1]+(B)[1];\
										(C)[2] = (A)[2]+(B)[2];\
										(C)[3] = (A)[3]+(B)[3];}
#define M2_COPY(A,C) {(C)[0] = (A)[0]; (C)[1] = (A)[1];\
	(C)[2] = (A)[2]; (C)[3] = (A)[3];}
#define MM2_COPY(A,C) {(C)[0][0] = (A)[0][0]; (C)[1][0] = (A)[1][0];\
	(C)[0][1] = (A)[0][1]; (C)[1][1] = (A)[1][1];}
#define M2_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0]+(b)*(B)[0];\
										(C)[1] = (a)*(A)[1]+(b)*(B)[1]);\
										(C)[2] = (a)*(A)[2]+(b)*(B)[2]);\
										(C)[3] = (a)*(A)[3]+(b)*(B)[3]);}
#define M2_TIMES_M2(A,B,C)			   {(C)[0] = (A)[0]*(B)[0]+(A)[2]*(B)[1];\
										(C)[1] = (A)[1]*(B)[0]+(A)[3]*(B)[1];\
										(C)[2] = (A)[0]*(B)[2]+(A)[2]*(B)[3];\
										(C)[3] = (A)[1]*(B)[2]+(A)[3]*(B)[3];}
#define M2_INVERT(M,IM,det)                       \
{ DOUBLE invdet;                                  \
  det = (M)[0][0]*(M)[1][1]-(M)[1][0]*(M)[0][1];  \
  if (ABS((det))<SMALL_D*SMALL_D)                 \
	RETURN(1);                                    \
  invdet = 1.0 / (det);                           \
  (IM)[0][0] =  (M)[1][1]*invdet;                 \
  (IM)[1][0] = -(M)[0][1]*invdet;                 \
  (IM)[0][1] = -(M)[1][0]*invdet;                 \
  (IM)[1][1] =  (M)[0][0]*invdet;}

/* macros for vector operations */
#define V3_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
										(C)[1] = (a)*(A)[1] + (b)*(B)[1];\
										(C)[2] = (a)*(A)[2] + (b)*(B)[2];}
#define V3_COPY(A,C)				   {(C)[0] = (A)[0];\
										(C)[1] = (A)[1];\
										(C)[2] = (A)[2];}
#define V3_SUBTRACT(A,B,C)			   {(C)[0] = (A)[0] - (B)[0];\
										(C)[1] = (A)[1] - (B)[1];\
										(C)[2] = (A)[2] - (B)[2];}
#define V3_AVG2(A,B,C)				   {(C)[0] = 0.5*((A)[0] + (B)[0]);\
										(C)[1] = 0.5*((A)[1] + (B)[1]);\
										(C)[2] = 0.5*((A)[2] + (B)[2]);}
#define V3_AVG4(A,B,C,D,E)			   {(E)[0] = 0.25*((A)[0]+(B)[0]+(C)[0]+(D)[0]);\
										(E)[1] = 0.25*((A)[1]+(B)[1]+(C)[1]+(D)[1]);\
										(E)[2] = 0.25*((A)[2]+(B)[2]+(C)[2]+(D)[2]);}
#define V3_ADD(A,B,C)				   {(C)[0] = (A)[0] + (B)[0];\
										(C)[1] = (A)[1] + (B)[1];\
										(C)[2] = (A)[2] + (B)[2];}
#define V3_ADD1(A,C)				   {(C)[0] += (A)[0];\
										(C)[1] += (A)[1];\
										(C)[2] += (A)[2];}
#define V3_SCALE(c,C)				   {(C)[0] = (c)*(C)[0];\
										(C)[1] = (c)*(C)[1];\
										(C)[2] = (c)*(C)[2];}
#define V3_SCALEADD1(c,A,C)			   {(C)[0] += (c)*(A)[0];\
										(C)[1] += (c)*(A)[1];\
										(C)[2] += (c)*(A)[2];}
#define V3_SCALESET(c,A,C)			   {(C)[0] = (c)*(A)[0];\
										(C)[1] = (c)*(A)[1];\
										(C)[2] = (c)*(A)[2];}
#define V3_VECTOR_PRODUCT(A,B,C)	   {(C)[0] = (A)[1]*(B)[2] - (A)[2]*(B)[1];\
										(C)[1] = (A)[2]*(B)[0] - (A)[0]*(B)[2];\
										(C)[2] = (A)[0]*(B)[1] - (A)[1]*(B)[0];}
#define V3_EUKLIDNORM(A,b)				(b) = (sqrt((double)((A)[0]*(A)[0]+(A)[1]*(A)[1]+(A)[2]*(A)[2])));
#define V3_COMPARE(A,B,c)				(c) = (((ABS((A)[0]-(B)[0])<SMALL_C)&&(ABS((A)[1]-(B)[1])<SMALL_C)&&(ABS((A)[2]-(B)[2])<SMALL_C))?(0):(1));
#define V3_ISEQUAL(A,B) 				((ABS((A)[0]-(B)[0])<SMALL_C)&&(ABS((A)[1]-(B)[1])<SMALL_C)&&(ABS((A)[2]-(B)[2])<SMALL_C))
#define V3_EUKLIDNORM_OF_DIFF(A,B,b)	(b) = (sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])+((A)[2]-(B)[2])*((A)[2]-(B)[2]))));
#define V3_CLEAR(A) 				   {(A)[0] = 0.0; (A)[1]= 0.0; (A)[2] = 0.0;}
#define V3_SCALAR_PRODUCT(A,B,c)		(c) = ((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]);
#define V3_SCAL_PROD(A,B)		        ((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2])
#define V3_ISZERO(A)					((A)[0]==0.0 && (A)[1]==0.0 && (A)[2]==0.0)
#define V3_SUP(v,s)                    {s = MAX(ABS(v[0]),MAX(ABS(v[1]),ABS(v[2])));}

/* macros for matrix-vector operations */
#define M3_TIMES_V3(M,A,B)			   {(B)[0] = (M)[0]*(A)[0] + (M)[3]*(A)[1] + (M)[6]*(A)[2];\
										(B)[1] = (M)[1]*(A)[0] + (M)[4]*(A)[1] + (M)[7]*(A)[2];\
										(B)[2] = (M)[2]*(A)[0] + (M)[5]*(A)[1] + (M)[8]*(A)[2];}

#define MM3_TIMES_V3(M,A,B)			   {(B)[0] = (M)[0][0]*(A)[0] + (M)[0][1]*(A)[1] + (M)[0][2]*(A)[2];\
										(B)[1] = (M)[1][0]*(A)[0] + (M)[1][1]*(A)[1] + (M)[1][2]*(A)[2];\
										(B)[2] = (M)[2][0]*(A)[0] + (M)[2][1]*(A)[1] + (M)[2][2]*(A)[2];}

#define MT3_TIMES_V3(M,A,B)			   {(B)[0] = (M)[0][0]*(A)[0] + (M)[1][0]*(A)[1] + (M)[2][0]*(A)[2];\
										(B)[1] = (M)[0][1]*(A)[0] + (M)[1][1]*(A)[1] + (M)[2][1]*(A)[2];\
										(B)[2] = (M)[0][2]*(A)[0] + (M)[1][2]*(A)[1] + (M)[2][2]*(A)[2];}

#define MD3_TIMES_V3(M,A,B)			   {(B)[0] = (M)[0]*(A)[0];\
										(B)[1] = (M)[1]*(A)[1];\
										(B)[2] = (M)[2]*(A)[2];}

/* macros for matrix operations */
#define M3_SCALE(c,M)				   {(M)[0] = (c)*(M)[0];\
										(M)[1] = (c)*(M)[1];\
										(M)[2] = (c)*(M)[2];\
										(M)[3] = (c)*(M)[3];\
										(M)[4] = (c)*(M)[4];\
										(M)[5] = (c)*(M)[5];\
										(M)[6] = (c)*(M)[6];\
										(M)[7] = (c)*(M)[7];\
										(M)[8] = (c)*(M)[8];}
#define M3_ADD(A,B,C) 		   {(C)[0] = (A)[0]+(B)[0];\
										(C)[1] = (A)[1]+(B)[1];\
										(C)[2] = (A)[2]+(B)[2];\
										(C)[3] = (A)[3]+(B)[3];\
										(C)[4] = (A)[4]+(B)[4];\
										(C)[5] = (A)[5]+(B)[5];\
										(C)[6] = (A)[6]+(B)[6];\
										(C)[7] = (A)[7]+(B)[7];\
										(C)[8] = (A)[8]+(B)[8];}
#define M3_COPY(A,C) {(C)[0] = (A)[0]; (C)[1] = (A)[1]; (C)[2] = (A)[2];\
(C)[3] = (A)[3]; (C)[4] = (A)[4]; (C)[5] = (A)[5];\
(C)[6] = (A)[6]; (C)[7] = (A)[7]; (C)[8] = (A)[8];}
#define MM3_COPY(A,C) {(C)[0][0] = (A)[0][0]; (C)[1][0] = (A)[1][0]; (C)[2][0] = (A)[2][0];\
(C)[0][1] = (A)[0][1]; (C)[1][1] = (A)[1][1]; (C)[2][1] = (A)[2][1];\
(C)[0][2] = (A)[0][2]; (C)[1][2] = (A)[1][2]; (C)[2][2] = (A)[2][2];}
#define M3_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0]+(b)*(B)[0];\
										(C)[1] = (a)*(A)[1]+(b)*(B)[1];\
										(C)[2] = (a)*(A)[2]+(b)*(B)[2];\
										(C)[3] = (a)*(A)[3]+(b)*(B)[3];\
										(C)[4] = (a)*(A)[4]+(b)*(B)[4];\
										(C)[5] = (a)*(A)[5]+(b)*(B)[5];\
										(C)[6] = (a)*(A)[6]+(b)*(B)[6];\
										(C)[7] = (a)*(A)[7]+(b)*(B)[7];\
										(C)[8] = (a)*(A)[8]+(b)*(B)[8];}
#define M3_TIMES_M3(A,B,C)			   {(C)[0] = (A)[0]*(B)[0]+(A)[3]*(B)[1]+(A)[6]*(B)[2];\
										(C)[1] = (A)[1]*(B)[0]+(A)[4]*(B)[1]+(A)[7]*(B)[2];\
										(C)[2] = (A)[2]*(B)[0]+(A)[5]*(B)[1]+(A)[8]*(B)[2];\
										(C)[3] = (A)[0]*(B)[3]+(A)[3]*(B)[4]+(A)[6]*(B)[5];\
										(C)[4] = (A)[1]*(B)[3]+(A)[4]*(B)[4]+(A)[7]*(B)[5];\
										(C)[5] = (A)[2]*(B)[3]+(A)[5]*(B)[4]+(A)[8]*(B)[5];\
										(C)[6] = (A)[0]*(B)[6]+(A)[3]*(B)[7]+(A)[6]*(B)[8];\
										(C)[7] = (A)[1]*(B)[6]+(A)[4]*(B)[7]+(A)[7]*(B)[8];\
										(C)[8] = (A)[2]*(B)[6]+(A)[5]*(B)[7]+(A)[8]*(B)[8];}

#define M3_INVERT(M,IM,det)                       \
{ DOUBLE invdet;                                  \
  (det) = (M)[0][0]*(M)[1][1]*(M)[2][2]           \
	+ (M)[0][1]*(M)[1][2]*(M)[2][0]               \
	  + (M)[0][2]*(M)[1][0]*(M)[2][1]             \
		- (M)[0][2]*(M)[1][1]*(M)[2][0]           \
		  - (M)[0][0]*(M)[1][2]*(M)[2][1]         \
			- (M)[0][1]*(M)[1][0]*(M)[2][2];      \
  if (ABS((det))<SMALL_D*SMALL_D)                 \
	RETURN(1);                                    \
  invdet = 1.0 / (det);                           \
  (IM)[0][0] = ( (M)[1][1]*(M)[2][2] - (M)[1][2]*(M)[2][1]) * invdet;  \
  (IM)[0][1] = (-(M)[1][0]*(M)[2][2] + (M)[1][2]*(M)[2][0]) * invdet;  \
  (IM)[0][2] = ( (M)[1][0]*(M)[2][1] - (M)[1][1]*(M)[2][0]) * invdet;  \
  (IM)[1][0] = (-(M)[0][1]*(M)[2][2] + (M)[0][2]*(M)[2][1]) * invdet;  \
  (IM)[1][1] = ( (M)[0][0]*(M)[2][2] - (M)[0][2]*(M)[2][0]) * invdet;  \
  (IM)[1][2] = (-(M)[0][0]*(M)[2][1] + (M)[0][1]*(M)[2][0]) * invdet;  \
  (IM)[2][0] = ( (M)[0][1]*(M)[1][2] - (M)[0][2]*(M)[1][1]) * invdet;  \
  (IM)[2][1] = (-(M)[0][0]*(M)[1][2] + (M)[0][2]*(M)[1][0]) * invdet;  \
  (IM)[2][2] = ( (M)[0][0]*(M)[1][1] - (M)[0][1]*(M)[1][0]) * invdet;} 

#ifdef __MPW32__
#define M4_TIMES_M4(A,B,C)				{int i,j,k; for (i=0;i<4;i++) for (j=0;j<4;j++) for (k=0;k<4;k++) \
																(C)[i+4*j] = (A)[i+4*k] * (B)[k+4*j];}
#else
#define M4_TIMES_M4(A,B,C)			   {(C)[0]	= (A)[ 0]*(B)[ 0]+(A)[ 4]*(B)[ 1]+(A)[ 8]*(B)[ 2]+(A)[12]*(B)[ 3];\
										(C)[1]	= (A)[ 1]*(B)[ 0]+(A)[ 5]*(B)[ 1]+(A)[ 9]*(B)[ 2]+(A)[13]*(B)[ 3];\
										(C)[2]	= (A)[ 2]*(B)[ 0]+(A)[ 6]*(B)[ 1]+(A)[10]*(B)[ 2]+(A)[14]*(B)[ 3];\
										(C)[3]	= (A)[ 3]*(B)[ 0]+(A)[ 7]*(B)[ 1]+(A)[11]*(B)[ 2]+(A)[15]*(B)[ 3];\
										(C)[4]	= (A)[ 0]*(B)[ 4]+(A)[ 4]*(B)[ 5]+(A)[ 8]*(B)[ 6]+(A)[12]*(B)[ 7];\
										(C)[5]	= (A)[ 1]*(B)[ 4]+(A)[ 5]*(B)[ 5]+(A)[ 9]*(B)[ 6]+(A)[13]*(B)[ 7];\
										(C)[6]	= (A)[ 2]*(B)[ 4]+(A)[ 6]*(B)[ 5]+(A)[10]*(B)[ 6]+(A)[14]*(B)[ 7];\
										(C)[7]	= (A)[ 3]*(B)[ 4]+(A)[ 7]*(B)[ 5]+(A)[11]*(B)[ 6]+(A)[15]*(B)[ 7];\
										(C)[8]	= (A)[ 0]*(B)[ 8]+(A)[ 4]*(B)[ 9]+(A)[ 8]*(B)[10]+(A)[12]*(B)[11];\
										(C)[9]	= (A)[ 1]*(B)[ 8]+(A)[ 5]*(B)[ 9]+(A)[ 9]*(B)[10]+(A)[13]*(B)[11];\
										(C)[10] = (A)[ 2]*(B)[ 8]+(A)[ 6]*(B)[ 9]+(A)[10]*(B)[10]+(A)[14]*(B)[11];\
										(C)[11] = (A)[ 3]*(B)[ 8]+(A)[ 7]*(B)[ 9]+(A)[11]*(B)[10]+(A)[15]*(B)[11];\
										(C)[12] = (A)[ 0]*(B)[12]+(A)[ 4]*(B)[13]+(A)[ 8]*(B)[14]+(A)[12]*(B)[15];\
										(C)[13] = (A)[ 1]*(B)[12]+(A)[ 5]*(B)[13]+(A)[ 9]*(B)[14]+(A)[13]*(B)[15];\
										(C)[14] = (A)[ 2]*(B)[12]+(A)[ 6]*(B)[13]+(A)[10]*(B)[14]+(A)[14]*(B)[15];\
										(C)[15] = (A)[ 3]*(B)[12]+(A)[ 7]*(B)[13]+(A)[11]*(B)[14]+(A)[15]*(B)[15];}
#endif

/****************************************************************************/
/*																			*/
/* typedef of DIM-routines													*/
/*																			*/
/****************************************************************************/

#ifdef __TWODIM__

#define V_DIM_LINCOMB(a,A,b,B,C)		V2_LINCOMB(a,A,b,B,C)
#define V_DIM_COPY(A,C)				    V2_COPY(A,C)
#define V_DIM_SUBTRACT(A,B,C)			V2_SUBTRACT(A,B,C)
#define V_DIM_ADD(A,B,C)			    V2_ADD(A,B,C)
#define V_DIM_ADD1(A,C)		            V2_ADD1(A,C)
#define V_DIM_AVG2(A,B,C)			    V2_AVG2(A,B,C)
#define V_DIM_AVG4(A,B,C,D,E)			V2_AVG4(A,B,C,D,E)
#define V_DIM_SCALE(c,C)			    V2_SCALE(c,C)
#define V_DIM_SCALEADD1(c,A,C)			V2_SCALEADD1(c,A,C)
#define V_DIM_SCALESET(c,A,C)			V2_SCALESET(c,A,C)
#define V_DIM_VECTOR_PRODUCT(A,B,c)		V2_VECTOR_PRODUCT(A,B,c)
#define V_DIM_COMPARE(A,B,c)			V2_COMPARE(A,B,c)
#define V_DIM_ISEQUAL(A,B) 			    V2_ISEQUAL(A,B)
#define V_DIM_EUKLIDNORM(A,b)			V2_EUKLIDNORM(A,b)
#define V_DIM_EUKLIDNORM_OF_DIFF(A,B,b)	V2_EUKLIDNORM_OF_DIFF(A,B,b)
#define V_DIM_CLEAR(A) 				    V2_CLEAR(A)
#define V_DIM_SCALAR_PRODUCT(A,B,c)		V2_SCALAR_PRODUCT(A,B,c)
#define V_DIM_ISZERO(A)				    V2_ISZERO(A)
#define M_TIMES_V_DIM(M,A,B)	        M2_TIMES_V2(M,A,B)
#define MM_TIMES_V_DIM(M,A,B)	      	MM2_TIMES_V2(M,A,B)
#define MT_TIMES_V_DIM(M,A,B)	      	MT2_TIMES_V2(M,A,B)
#define MD_TIMES_V_DIM(M,A,B)	      	MD2_TIMES_V2(M,A,B)
#define M_DIM_ADD(A,B,C)			    M2_ADD(A,B,C)
#define M_DIM_COPY(A,C)				    M2_COPY(A,C)
#define MM_DIM_COPY(A,C)				MM2_COPY(A,C)
#define M_DIM_SCALE(c,M)			    M2_SCALE(c,M)
#define V_DIM_Normalize(a)			    V2_Normalize(a)
#define M_DIM_INVERT(M,IM,det)			M2_INVERT(M,IM,det)
#define V_DIM_SUP(v,s)			        V2_SUP(v,s)

#endif

#ifdef __THREEDIM__

#define V_DIM_LINCOMB(a,A,b,B,C)		V3_LINCOMB(a,A,b,B,C)
#define V_DIM_COPY(A,C)				    V3_COPY(A,C)
#define V_DIM_SUBTRACT(A,B,C)			V3_SUBTRACT(A,B,C)
#define V_DIM_ADD(A,B,C)			    V3_ADD(A,B,C)
#define V_DIM_AVG2(A,B,C)			    V3_AVG2(A,B,C)
#define V_DIM_AVG4(A,B,C,D,E)			V3_AVG4(A,B,C,D,E)
#define V_DIM_ADD1(A,C)			        V3_ADD1(A,C)
#define V_DIM_SCALE(c,C)			    V3_SCALE(c,C)
#define V_DIM_SCALEADD1(c,A,C)			V3_SCALEADD1(c,A,C)
#define V_DIM_SCALESET(c,A,C)			V3_SCALESET(c,A,C)
#define V_DIM_VECTOR_PRODUCT(A,B,c)		V3_VECTOR_PRODUCT(A,B,c)
#define V_DIM_COMPARE(A,B,c)			V3_COMPARE(A,B,c)
#define V_DIM_ISEQUAL(A,B) 			    V3_ISEQUAL(A,B)
#define V_DIM_EUKLIDNORM(A,b)			V3_EUKLIDNORM(A,b)
#define V_DIM_EUKLIDNORM_OF_DIFF(A,B,b)	V3_EUKLIDNORM_OF_DIFF(A,B,b)
#define V_DIM_CLEAR(A) 				    V3_CLEAR(A)
#define V_DIM_SCALAR_PRODUCT(A,B,c)		V3_SCALAR_PRODUCT(A,B,c)
#define V_DIM_ISZERO(A)				    V3_ISZERO(A)
#define M_TIMES_V_DIM(M,A,B)			M3_TIMES_V3(M,A,B)
#define MM_TIMES_V_DIM(M,A,B)			MM3_TIMES_V3(M,A,B)
#define MT_TIMES_V_DIM(M,A,B)			MT3_TIMES_V3(M,A,B)
#define MD_TIMES_V_DIM(M,A,B)			MD3_TIMES_V3(M,A,B)
#define M_DIM_ADD(A,B,C) 		   	    M3_ADD(A,B,C)
#define M_DIM_COPY(A,C)				    M3_COPY(A,C)
#define MM_DIM_COPY(A,C)				MM3_COPY(A,C)
#define M_DIM_SCALE(c,M)			    M3_SCALE(c,M)
#define V_DIM_Normalize(a)			    V3_Normalize(a)
#define M_DIM_INVERT(M,IM,det)			M3_INVERT(M,IM,det)
#define V_DIM_SUP(v,s)			        V3_SUP(v,s)

#endif

/****************************************************************************/
/*																			*/
/* typedef of 2d points for screen coordinates								*/
/*																			*/
/****************************************************************************/

struct coord_point
{
	COORD x;
	COORD y;
};

typedef struct coord_point COORD_POINT;

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* general routines */
INT 		ClipRectangleAgainstRectangle		(const COORD *r1min, const COORD *r1max, COORD *r2min, COORD *r2max);
INT 		CheckRectagleIntersection			(const COORD *r1min, const COORD *r1max, const COORD *r2min, const COORD *r2max);
INT 		CheckRectangle						(const COORD *rmin, const COORD *rmax, const COORD minsize);
INT 		PointInTriangle 					(const COORD_POINT *Points, const COORD_POINT Point);
INT 		PointInPolygon						(const COORD_POINT *Points, INT n, COORD_POINT Point);
INT 		PointInPolygonC 					(const COORD_VECTOR_2D *Points, INT n, const COORD_VECTOR_2D Point);
INT 		PolyArea 							(INT n, COORD_VECTOR_2D *Polygon, COORD *Area);


/* 2D routines */
INT 		M2_Invert							(COORD *Inverse, const COORD *Matrix);
DOUBLE		vp									(const DOUBLE x1, const DOUBLE y1, const DOUBLE x2, const DOUBLE y2);
INT 		V2_Normalize						(COORD *a);
INT 		V2_Rotate							(COORD *vector, COORD alpha);
DOUBLE		tarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2);
DOUBLE		qarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,DOUBLE x3,DOUBLE y3);
DOUBLE		c_tarea								(const COORD *x0, const COORD *x1, const COORD *x2);
DOUBLE		c_qarea								(const COORD *x0, const COORD *x1, const COORD *x2, const COORD *x3);
DOUBLE		ctarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2);
DOUBLE		cqarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,DOUBLE x3,DOUBLE y3);


/* 3D routines */
INT 		M3_Invert							(COORD *Inverse, const COORD *Matrix);
INT 		V3_Normalize						(COORD *a);
INT 		V3_NormVectorProduct				(const COORD *a, const COORD *b, COORD *result);
INT 		V3_Rotate							(COORD *vector, const COORD *axis, COORD alpha);
INT 		V3_Angle							(const COORD *a, const COORD *b, COORD *result);
INT 		V3_Orthogonalize					(const COORD *a, const COORD *b, COORD *r);
INT 		V3_Project 							(const COORD *a, const COORD *b, COORD *r);


/* 4D routines */
INT 		M4_Invert							(COORD *Inverse, const COORD *Matrix);

#endif
