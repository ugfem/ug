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


/* RCS_ID
$Header$
*/

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
#define CP_LIMCOMB(a,A,b,B,C)		   {(C).x = (DOUBLE)(a)*(A).x + (DOUBLE)(b)*(B).x;\
										(C).y = (DOUBLE)(a)*(A).y + (DOUBLE)(b)*(B).y;}
#define CP_SCALARPRODUCT(A,B,c) 		(c) = (A).x*(B).x + (A).y*(B).y;
#define CP_EUKLIDNORM(A,b)				(b) = (DOUBLE)sqrt((double)((A).x*(A).x+(A).y*(A).y));

/* macros for 1D vector operations */
#define V1_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0] + (b)*(B)[0];}
#define V1_SET(a,A)					   {(A)[0] = a;}
#define V1_COPY(A,C)				   {(C)[0] = (A)[0];}
#define V1_SUBTRACT(A,B,C)			   {(C)[0] = (A)[0] - (B)[0];}
#define V1_ADD(A,B,C)				   {(C)[0] = (A)[0] + (B)[0];}
#define V1_ADD1(A,C)				   {(C)[0] += (A)[0];}
#define V1_AVG2(p1,p2,p)			   {(p)[0] = 0.5*((p1)[0]+(p2)[0]);}
#define V1_AVG3(p1,p2,p3,p)			   {(p)[0] = ((p1)[0]+(p2)[0]+(p3)[0])/3.0;
#define V1_AVG4(p1,p2,p3,p4,p)		   {(p)[0] = 0.25*((p1)[0]+(p2)[0]+(p3)[0]+(p4)[0]);}
#define V1_SCALE(c,C)				   {(C)[0] = (c)*(C)[0];}
#define V1_SCALEADD1(c,A,C)			   {(C)[0] += (c)*(A)[0];}
#define V1_SCALESET(c,A,C)			   {(C)[0] = (c)*(A)[0];}
#define V1_VECTOR_PRODUCT(A,B,c)		(c) = (A)[0]*(B)[1];
#define V1_COMPARE(A,B,c)				(c) = (ABS((A)[0]-(B)[0])<SMALL_C);
#define V1_ISEQUAL(A,B) 				(ABS((A)[0]-(B)[0])<SMALL_C)
#define V1_EUKLIDNORM(A,b)				(b) = (A)[0]);
#define V1_EUKLIDNORM_OF_DIFF(A,B,b)	(b) = sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])));
#define V1_CLEAR(A) 				   {(A)[0] = 0.0;}
#define V1_SCALAR_PRODUCT(A,B,c)		(c) = (A)[0]*(B)[0];
#define V1_SCAL_PROD(A,B)				((A)[0]*(B)[0])
#define V1_ISZERO(A)					(ABS((A)[0])<SMALL_C)
#define V1_SUP(v,s)                    {s = v[0]);}
#define V1_NORMAL(s,n)				   {(n)[0] = (s)[0];}

#define M1_INVERT(M,IM,det)				{det = (IM)[0];							\
										 if (ABS((det))<SMALL_D*SMALL_D)		\
										 	return(1);							\
										 (M)[0] = 1./det;}
#define MT1_TIMES_V1(M,A,B)			   {(B)[0] = (M)[0][0]*(A)[0];}

/* macros for 2D vector operations */
#define V2_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
										(C)[1] = (a)*(A)[1] + (b)*(B)[1];}
#define V2_SET(a,A)					   {(A)[0] = a;\
										(A)[1] = a;}
#define V2_COPY(A,C)				   {(C)[0] = (A)[0];\
										(C)[1] = (A)[1];}
#define V2_SUBTRACT(A,B,C)			   {(C)[0] = (A)[0] - (B)[0];\
										(C)[1] = (A)[1] - (B)[1];}
#define V2_ADD(A,B,C)				   {(C)[0] = (A)[0] + (B)[0];\
										(C)[1] = (A)[1] + (B)[1];}
#define V2_ADD1(A,C)				   {(C)[0] += (A)[0];\
										(C)[1] += (A)[1];}
#define V2_AVG2(p1,p2,p)			   {(p)[0] = 0.5*((p1)[0]+(p2)[0]); \
										(p)[1] = 0.5*((p1)[1]+(p2)[1]);}
#define V2_AVG3(p1,p2,p3,p)			   {(p)[0] = ((p1)[0]+(p2)[0]+(p3)[0])/3.0; \
										(p)[1] = ((p1)[1]+(p2)[1]+(p3)[1])/3.0;}
#define V2_AVG4(p1,p2,p3,p4,p)		   {(p)[0] = 0.25*((p1)[0]+(p2)[0]+(p3)[0]+(p4)[0]); \
										(p)[1] = 0.25*((p1)[1]+(p2)[1]+(p3)[1]+(p4)[1]);}
#define V2_SCALE(c,C)				   {(C)[0] = (c)*(C)[0];\
										(C)[1] = (c)*(C)[1];}
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
#define V2_SCAL_PROD(A,B)				((A)[0]*(B)[0]+(A)[1]*(B)[1])
#define V2_ISZERO(A)					(ABS((A)[0])<SMALL_C && ABS((A)[1])<SMALL_C)
#define V2_SUP(v,s)                    {s = MAX(ABS(v[0]),ABS(v[1]));}
#define V2_NORMAL(s,n)				   {(n)[0] = (s)[1]; (n)[1] = -(s)[0];}


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
#define M2_COPY(A,C)				   {(C)[0] = (A)[0];\
										(C)[1] = (A)[1];\
										(C)[2] = (A)[2];\
										(C)[3] = (A)[3];}
#define MM2_COPY(A,C)				   {(C)[0][0] = (A)[0][0];\
										(C)[0][1] = (A)[0][1];\
										(C)[1][0] = (A)[1][0];\
										(C)[1][1] = (A)[1][1];}
#define M2_DET(M)						((M)[0]*(M)[3] - (M)[1]*(M)[2])
#define MM2_DET(M)						((M)[0][0]*(M)[1][1] - (M)[1][0]*(M)[0][1])
#define M2_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0]+(b)*(B)[0];\
										(C)[1] = (a)*(A)[1]+(b)*(B)[1]);\
										(C)[2] = (a)*(A)[2]+(b)*(B)[2]);\
										(C)[3] = (a)*(A)[3]+(b)*(B)[3]);}
#define M2_TIMES_M2(A,B,C)			   {(C)[0] = (A)[0]*(B)[0]+(A)[2]*(B)[1];\
										(C)[1] = (A)[1]*(B)[0]+(A)[3]*(B)[1];\
										(C)[2] = (A)[0]*(B)[2]+(A)[2]*(B)[3];\
										(C)[3] = (A)[1]*(B)[2]+(A)[3]*(B)[3];}
#define M2_INVERT_STD(M,IM,det)                       \
{ DOUBLE invdet;                                  \
  det = (M)[0]*(M)[3]-(M)[1]*(M)[2];  \
	if (ABS((det))<SMALL_D*SMALL_D) det= 0.;  \
	else {                                      \
	invdet = 1.0 / (det);                       \
	(IM)[0] =  (M)[3]*invdet;             \
	(IM)[1] = -(M)[1]*invdet;             \
	(IM)[2] = -(M)[2]*invdet;             \
	(IM)[3] =  (M)[0]*invdet;}}

#define M2_INVERT(M,IM,det)                   \
{ DOUBLE invdet;                                  \
  det = (M)[0][0]*(M)[1][1]-(M)[1][0]*(M)[0][1];  \
	if (ABS((det))<SMALL_D*SMALL_D) det= 0.;  \
	else {                                      \
	invdet = 1.0 / (det);                       \
	(IM)[0][0] =  (M)[1][1]*invdet;             \
	(IM)[1][0] = -(M)[1][0]*invdet;             \
	(IM)[0][1] = -(M)[0][1]*invdet;             \
	(IM)[1][1] =  (M)[0][0]*invdet;}}

#define M2_MAXNORM(M,n)					(n)=MAX(ABS((M)[0])+ABS((M)[2]),ABS((M)[1])+ABS((M)[3]))

/* macros for vector operations */
#define V3_LINCOMB(a,A,b,B,C)		   {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
										(C)[1] = (a)*(A)[1] + (b)*(B)[1];\
										(C)[2] = (a)*(A)[2] + (b)*(B)[2];}
#define V3_SET(a,A)					   {(A)[0] = a;\
										(A)[1] = a;\
										(A)[2] = a;}
#define V3_COPY(A,C)				   {(C)[0] = (A)[0];\
										(C)[1] = (A)[1];\
										(C)[2] = (A)[2];}
#define V3_SUBTRACT(A,B,C)			   {(C)[0] = (A)[0] - (B)[0];\
										(C)[1] = (A)[1] - (B)[1];\
										(C)[2] = (A)[2] - (B)[2];}
#define V3_ADD(A,B,C)				   {(C)[0] = (A)[0] + (B)[0];\
										(C)[1] = (A)[1] + (B)[1];\
										(C)[2] = (A)[2] + (B)[2];}
#define V3_ADD1(A,C)				   {(C)[0] += (A)[0];\
										(C)[1] += (A)[1];\
										(C)[2] += (A)[2];}
#define V3_AVG2(p1,p2,p)			   {(p)[0] = 0.5*((p1)[0]+(p2)[0]); \
										(p)[1] = 0.5*((p1)[1]+(p2)[1]); \
										(p)[2] = 0.5*((p1)[2]+(p2)[2]);}
#define V3_AVG3(p1,p2,p3,p)			   {(p)[0] = ((p1)[0]+(p2)[0]+(p3)[0])/3.0; \
										(p)[1] = ((p1)[1]+(p2)[1]+(p3)[1])/3.0; \
										(p)[2] = ((p1)[2]+(p2)[2]+(p3)[2])/3.0;}
#define V3_AVG4(p1,p2,p3,p4,p)		   {(p)[0] = 0.25*((p1)[0]+(p2)[0]+(p3)[0]+(p4)[0]); \
										(p)[1] = 0.25*((p1)[1]+(p2)[1]+(p3)[1]+(p4)[1]); \
										(p)[2] = 0.25*((p1)[2]+(p2)[2]+(p3)[2]+(p4)[2]);}
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
#define M3_COPY(A,C)		 		   {(C)[0] = (A)[0];\
										(C)[1] = (A)[1];\
										(C)[2] = (A)[2];\
										(C)[3] = (A)[3];\
										(C)[4] = (A)[4];\
										(C)[5] = (A)[5];\
										(C)[6] = (A)[6];\
										(C)[7] = (A)[7];\
										(C)[8] = (A)[8];}
#define MM3_COPY(A,C)		 		   {(C)[0][0] = (A)[0][0]; \
										(C)[0][1] = (A)[0][1]; \
										(C)[0][2] = (A)[0][2]; \
										(C)[1][0] = (A)[1][0]; \
										(C)[1][1] = (A)[1][1]; \
										(C)[1][2] = (A)[1][2]; \
										(C)[2][0] = (A)[2][0]; \
										(C)[2][1] = (A)[2][1]; \
										(C)[2][2] = (A)[2][2];}
#define M3_DET(M)						((M)[0]*(M)[4]*(M)[8] + (M)[1]*(M)[5]*(M)[6] + (M)[2]*(M)[3]*(M)[7] \
		    							-(M)[2]*(M)[4]*(M)[6] - (M)[0]*(M)[5]*(M)[7] - (M)[1]*(M)[3]*(M)[8])
#define MM3_DET(M)						((M)[0][0]*(M)[1][1]*(M)[2][2] + (M)[0][1]*(M)[1][2]*(M)[2][0] + (M)[0][2]*(M)[1][0]*(M)[2][1] \
										-(M)[0][2]*(M)[1][1]*(M)[2][0] - (M)[0][0]*(M)[1][2]*(M)[2][1] - (M)[0][1]*(M)[1][0]*(M)[2][2])
#define M3_ADD(A,B,C) 		   			{(C)[0] = (A)[0]+(B)[0];\
										(C)[1] = (A)[1]+(B)[1];\
										(C)[2] = (A)[2]+(B)[2];\
										(C)[3] = (A)[3]+(B)[3];\
										(C)[4] = (A)[4]+(B)[4];\
										(C)[5] = (A)[5]+(B)[5];\
										(C)[6] = (A)[6]+(B)[6];\
										(C)[7] = (A)[7]+(B)[7];\
										(C)[8] = (A)[8]+(B)[8];}
#define M3_ADDMATRIX(A,B,C) 		    M3_ADD(A,B,C)
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
		return(1);                                    \
	invdet = 1.0 / (det);                           \
	(IM)[0][0] = ( (M)[1][1]*(M)[2][2] - (M)[1][2]*(M)[2][1]) * invdet;  \
	(IM)[0][1] = (-(M)[0][1]*(M)[2][2] + (M)[0][2]*(M)[2][1]) * invdet;  \
	(IM)[0][2] = ( (M)[0][1]*(M)[1][2] - (M)[0][2]*(M)[1][1]) * invdet;  \
	(IM)[1][0] = (-(M)[1][0]*(M)[2][2] + (M)[1][2]*(M)[2][0]) * invdet;  \
	(IM)[1][1] = ( (M)[0][0]*(M)[2][2] - (M)[0][2]*(M)[2][0]) * invdet;  \
	(IM)[1][2] = (-(M)[0][0]*(M)[1][2] + (M)[0][2]*(M)[1][0]) * invdet;  \
	(IM)[2][0] = ( (M)[1][0]*(M)[2][1] - (M)[1][1]*(M)[2][0]) * invdet;  \
	(IM)[2][1] = (-(M)[0][0]*(M)[2][1] + (M)[0][1]*(M)[2][0]) * invdet;  \
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
#define M4_COPY(A,B)					{(B)[0] = (A)[0]; (B)[1] = (A)[1]; (B)[2] = (A)[2]; (B)[3] = (A)[3]; \
										 (B)[4] = (A)[4]; (B)[5] = (A)[5]; (B)[6] = (A)[6]; (B)[7] = (A)[7]; \
										 (B)[8] = (A)[8]; (B)[9] = (A)[9]; (B)[10] = (A)[10]; (B)[11] = (A)[11]; \
										 (B)[12] = (A)[12]; (B)[13] = (A)[13]; (B)[14] = (A)[14]; (B)[15] = (A)[15];}

#define M4_CLEAR(B)						{(B)[0] = 0.0; (B)[1] = 0.0; (B)[2] = 0.0; (B)[3] = 0.0; \
										 (B)[4] = 0.0; (B)[5] = 0.0; (B)[6] = 0.0; (B)[7] = 0.0; \
										 (B)[8] = 0.0; (B)[9] = 0.0; (B)[10] = 0.0; (B)[11] = 0.0; \
										 (B)[12] = 0.0; (B)[13] = 0.0; (B)[14] = 0.0; (B)[15] = 0.0;}

/* macros for exact solver (EX) */
#define EX_MAT(m,b,i,j)			((m)[2*(b)*(i) + (j)])

/****************************************************************************/
/*																			*/
/* typedef of DIM-routines													*/
/*																			*/
/****************************************************************************/

#ifdef __TWODIM__

#define V_BDIM_LINCOMB(a,A,b,B,C)		V1_LINCOMB(a,A,b,B,C)
#define V_BDIM_SET(a,A)				    V1_SET(a,A)
#define V_BDIM_COPY(A,C)				V1_COPY(A,C)
#define V_BDIM_SUBTRACT(A,B,C)			V1_SUBTRACT(A,B,C)
#define V_BDIM_ADD(A,B,C)			    V1_ADD(A,B,C)
#define V_BDIM_ADD1(A,C)		        V1_ADD1(A,C)
#define V_BDIM_AVG2(A,B,X)		        V1_AVG2(A,B,X)
#define V_BDIM_AVG3(A,B,C,X)		    V1_AVG3(A,B,C,X)
#define V_BDIM_AVG4(A,B,C,D,X)	        V1_AVG4(A,B,C,D,X)
#define V_BDIM_SCALE(c,C)			    V1_SCALE(c,C)
#define V_BDIM_SCALEADD1(c,A,C)			V1_SCALEADD1(c,A,C)
#define V_BDIM_SCALESET(c,A,C)			V1_SCALESET(c,A,C)
#define V_BDIM_VECTOR_PRODUCT(A,B,c)	V1_VECTOR_PRODUCT(A,B,c)
#define V_BDIM_COMPARE(A,B,c)			V1_COMPARE(A,B,c)
#define V_BDIM_ISEQUAL(A,B) 		    V1_ISEQUAL(A,B)
#define V_BDIM_EUKLIDNORM(A,b)			V1_EUKLIDNORM(A,b)
#define V_BDIM_EUKLIDNORM_OF_DIFF(A,B,b)V1_EUKLIDNORM_OF_DIFF(A,B,b)
#define V_BDIM_CLEAR(A) 			    V1_CLEAR(A)
#define V_BDIM_SCALAR_PRODUCT(A,B,c)	V1_SCALAR_PRODUCT(A,B,c)
#define V_BDIM_SCAL_PROD(A,B)			V1_SCAL_PROD(A,B)
#define V_BDIM_ISZERO(A)				V1_ISZERO(A)
#define V_BDIM_SUP(v,s)			        V1_SUP(v,s)
#define V_BDIM_Normalize(a)			    V1_Normalize(a)

#define M_BDIM_INVERT(M,IM,det)			M1_INVERT(M,IM,det)
#define MT_TIMES_V_BDIM(M,A,B)			MT1_TIMES_V1(M,A,B)

#define V_DIM_LINCOMB(a,A,b,B,C)		V2_LINCOMB(a,A,b,B,C)
#define V_DIM_SET(a,A)				    V2_SET(a,A)
#define V_DIM_COPY(A,C)				    V2_COPY(A,C)
#define V_DIM_SUBTRACT(A,B,C)			V2_SUBTRACT(A,B,C)
#define V_DIM_ADD(A,B,C)			    V2_ADD(A,B,C)
#define V_DIM_ADD1(A,C)		            V2_ADD1(A,C)
#define V_DIM_AVG2(A,B,X)		        V2_AVG2(A,B,X)
#define V_DIM_AVG3(A,B,C,X)		        V2_AVG3(A,B,C,X)
#define V_DIM_AVG4(A,B,C,D,X)	        V2_AVG4(A,B,C,D,X)
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
#define V_DIM_SCAL_PROD(A,B)			V2_SCAL_PROD(A,B)
#define V_DIM_SP(A,B)	        		V2_SCAL_PROD(A,B)
#define V_DIM_ISZERO(A)				    V2_ISZERO(A)
#define V_DIM_SUP(v,s)			        V2_SUP(v,s)
#define V_DIM_Normalize(a)			    V2_Normalize(a)
#define V_DIM_NORMAL(a)				    V2_NORMAL(a)

#define M_TIMES_V_DIM(M,A,B)	        M2_TIMES_V2(M,A,B)
#define MM_TIMES_V_DIM(M,A,B)	      	MM2_TIMES_V2(M,A,B)
#define MT_TIMES_V_DIM(M,A,B)	      	MT2_TIMES_V2(M,A,B)
#define MD_TIMES_V_DIM(M,A,B)	      	MD2_TIMES_V2(M,A,B)
#define M_DIM_ADD(A,B,C)			    M2_ADD(A,B,C)
#define M_DIM_COPY(A,C)				    M2_COPY(A,C)
#define MM_DIM_COPY(A,C)				MM2_COPY(A,C)
#define M_DIM_DET(M)					M2_DET(M)
#define MM_DIM_DET(M)					MM2_DET(M)
#define M_DIM_SCALE(c,M)			    M2_SCALE(c,M)
#define M_DIM_INVERT(M,IM,det)			M2_INVERT(M,IM,det)

#endif

#ifdef __THREEDIM__

#define V_BDIM_LINCOMB(a,A,b,B,C)		V2_LINCOMB(a,A,b,B,C)
#define V_BDIM_SET(a,A)				    V2_SET(a,A)
#define V_BDIM_COPY(A,C)				V2_COPY(A,C)
#define V_BDIM_SUBTRACT(A,B,C)			V2_SUBTRACT(A,B,C)
#define V_BDIM_ADD(A,B,C)			    V2_ADD(A,B,C)
#define V_BDIM_ADD1(A,C)		        V2_ADD1(A,C)
#define V_BDIM_AVG2(A,B,X)		        V2_AVG2(A,B,X)
#define V_BDIM_AVG3(A,B,C,X)		    V2_AVG3(A,B,C,X)
#define V_BDIM_AVG4(A,B,C,D,X)	        V2_AVG4(A,B,C,D,X)
#define V_BDIM_SCALE(c,C)			    V2_SCALE(c,C)
#define V_BDIM_SCALEADD1(c,A,C)			V2_SCALEADD1(c,A,C)
#define V_BDIM_SCALESET(c,A,C)			V2_SCALESET(c,A,C)
#define V_BDIM_VECTOR_PRODUCT(A,B,c)	V2_VECTOR_PRODUCT(A,B,c)
#define V_BDIM_COMPARE(A,B,c)			V2_COMPARE(A,B,c)
#define V_BDIM_ISEQUAL(A,B) 		    V2_ISEQUAL(A,B)
#define V_BDIM_EUKLIDNORM(A,b)			V2_EUKLIDNORM(A,b)
#define V_BDIM_EUKLIDNORM_OF_DIFF(A,B,b)V2_EUKLIDNORM_OF_DIFF(A,B,b)
#define V_BDIM_CLEAR(A) 			    V2_CLEAR(A)
#define V_BDIM_SCALAR_PRODUCT(A,B,c)	V2_SCALAR_PRODUCT(A,B,c)
#define V_BDIM_SCAL_PROD(A,B)			V2_SCAL_PROD(A,B)
#define V_BDIM_ISZERO(A)				V2_ISZERO(A)
#define V_BDIM_SUP(v,s)			        V2_SUP(v,s)
#define V_BDIM_Normalize(a)			    V2_Normalize(a)

#define M_BDIM_INVERT(M,IM,det)			M2_INVERT(M,IM,det)
#define MT_TIMES_V_BDIM(M,A,B)			MT2_TIMES_V2(M,A,B)

#define V_DIM_LINCOMB(a,A,b,B,C)		V3_LINCOMB(a,A,b,B,C)
#define V_DIM_SET(a,A)				    V3_SET(a,A)
#define V_DIM_COPY(A,C)				    V3_COPY(A,C)
#define V_DIM_SUBTRACT(A,B,C)			V3_SUBTRACT(A,B,C)
#define V_DIM_ADD(A,B,C)			    V3_ADD(A,B,C)
#define V_DIM_ADD1(A,C)			        V3_ADD1(A,C)
#define V_DIM_AVG2(A,B,X)		        V3_AVG2(A,B,X)
#define V_DIM_AVG3(A,B,C,X)		        V3_AVG3(A,B,C,X)
#define V_DIM_AVG4(A,B,C,D,X)	        V3_AVG4(A,B,C,D,X)
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
#define V_DIM_SCAL_PROD(A,B)			V3_SCAL_PROD(A,B)
#define V_DIM_SP(A,B)	         		V3_SCAL_PROD(A,B)
#define V_DIM_ISZERO(A)				    V3_ISZERO(A)
#define M_TIMES_V_DIM(M,A,B)			M3_TIMES_V3(M,A,B)
#define MM_TIMES_V_DIM(M,A,B)			MM3_TIMES_V3(M,A,B)
#define MT_TIMES_V_DIM(M,A,B)			MT3_TIMES_V3(M,A,B)
#define MD_TIMES_V_DIM(M,A,B)			MD3_TIMES_V3(M,A,B)
#define M_DIM_ADD(A,B,C) 		   	    M3_ADD(A,B,C)
#define M_DIM_COPY(A,C)				    M3_COPY(A,C)
#define MM_DIM_COPY(A,C)				MM3_COPY(A,C)
#define M_DIM_DET(M)					M3_DET(M)
#define MM_DIM_DET(M)					MM3_DET(M)
#define M_DIM_SCALE(c,M)			    M3_SCALE(c,M)
#define V_DIM_Normalize(a)			    V3_Normalize(a)
#define V_DIM_NORMAL(a)				    V3_NORMAL(a)
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
	DOUBLE x;
	DOUBLE y;
};

typedef struct coord_point COORD_POINT;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

extern const DOUBLE unit_vec[DIM][DIM];

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* general routines */
INT 		ClipRectangleAgainstRectangle		(const DOUBLE *r1min, const DOUBLE *r1max, DOUBLE *r2min, DOUBLE *r2max);
INT 		CheckRectagleIntersection			(const DOUBLE *r1min, const DOUBLE *r1max, const DOUBLE *r2min, const DOUBLE *r2max);
INT 		CheckRectangle						(const DOUBLE *rmin, const DOUBLE *rmax, const DOUBLE minsize);
INT 		PointInTriangle 					(const COORD_POINT *Points, const COORD_POINT Point);
INT 		PointInPolygon						(const COORD_POINT *Points, INT n, COORD_POINT Point);
INT 		PointInPolygonC 					(const DOUBLE_VECTOR_2D *Points, INT n, const DOUBLE_VECTOR_2D Point);
INT 		PolyArea 							(INT n, DOUBLE_VECTOR_2D *Polygon, DOUBLE *Area);
INT 		QuadraticFittedMin 					(DOUBLE *x, DOUBLE *y, INT n, DOUBLE *minx);
INT 		EXDecomposeMatrixFLOAT 				(FLOAT *Mat, INT bw, INT n);
INT 		EXDecomposeMatrixDOUBLE 			(DOUBLE *Mat, INT bw, INT n);
INT 		EXApplyLUFLOAT 						(FLOAT *Mat, INT bw, INT n, DOUBLE *Vec);
INT 		EXApplyLUDOUBLE 					(DOUBLE *Mat, INT bw, INT n, DOUBLE *Vec);
INT 		LineISTriangle3D 					(const DOUBLE *c1, const DOUBLE *c2, const DOUBLE *c3, const DOUBLE *p1, const DOUBLE *p2, DOUBLE *lambda);


/* 2D routines */
INT 		M2_Invert							(DOUBLE *Inverse, const DOUBLE *Matrix);
DOUBLE		vp									(const DOUBLE x1, const DOUBLE y1, const DOUBLE x2, const DOUBLE y2);
INT 		V2_Normalize						(DOUBLE *a);
INT 		V2_Rotate							(DOUBLE *vector, DOUBLE alpha);
INT			V2_IntersectLineSegments			(const DOUBLE_VECTOR a0, const DOUBLE_VECTOR a1, const DOUBLE_VECTOR b0, const DOUBLE_VECTOR b1, DOUBLE *lambda);
DOUBLE		tarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2);
DOUBLE		qarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,DOUBLE x3,DOUBLE y3);
DOUBLE		c_tarea								(const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2);
DOUBLE		c_qarea								(const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2, const DOUBLE *x3);
DOUBLE		ctarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2);
DOUBLE		cqarea								(DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,DOUBLE x3,DOUBLE y3);
DOUBLE		V_te								(const DOUBLE *x0, const DOUBLE *x1,
												 const DOUBLE *x2, const DOUBLE *x3);
DOUBLE		V_py								(const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2,
												 const DOUBLE *x3, const DOUBLE *x4);
DOUBLE		V_pr								(const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2,
												 const DOUBLE *x3, const DOUBLE *x4, const DOUBLE *x5);
DOUBLE		V_he								(const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2, const DOUBLE *x3,
												 const DOUBLE *x4, const DOUBLE *x5, const DOUBLE *x6, const DOUBLE *x7);


/* 3D routines */
INT 		M3_Invert							(DOUBLE *Inverse, const DOUBLE *Matrix);
INT 		V3_Normalize						(DOUBLE *a);
INT 		V3_NormVectorProduct				(const DOUBLE *a, const DOUBLE *b, DOUBLE *result);
INT 		V3_Rotate							(DOUBLE *vector, const DOUBLE *axis, DOUBLE alpha);
INT 		V3_Angle							(const DOUBLE *a, const DOUBLE *b, DOUBLE *result);
INT 		V3_Orthogonalize					(const DOUBLE *a, const DOUBLE *b, DOUBLE *r);
INT 		V3_Project 							(const DOUBLE *a, const DOUBLE *b, DOUBLE *r);


/* 4D routines */
INT 		M4_Invert							(DOUBLE *Inverse, const DOUBLE *Matrix);

/* volume calculations*/
DOUBLE		GeneralElementVolume				(INT tag, DOUBLE *x_co[]);
DOUBLE		ElementVolume						(const ELEMENT *elem);

#endif
