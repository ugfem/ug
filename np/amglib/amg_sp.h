// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_sp.h														*/
/*																			*/
/* Purpose:   interface to sparse matrix/vector data structure for amg          */
/*																			*/
/* Author:	  Peter Bastian                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   01 FEB 1996 Begin												*/
/*			  30 SEP 1997 redesign											*/
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

#ifndef __AMG_SP__
#define __AMG_SP__

#include "amg_header.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define AMG_LINES_PER_PAGE              60      /* repeat legend every this line		*/
#define AMG_COLS_PER_LINE               3       /* matrix entries per line                              */

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/******** a block vector, i.e. n*blocksize entries    ***********************/
/****************************************************************************/

struct amg_vector {                                     /* data + ptr to nodes structure		*/
  char name[AMG_NAME_SIZE];                     /* name of this vector					*/
  int n;                                                        /* dimension of vector in blocks		*/
  int b;                                                        /* dimension of block					*/
  double *x;                                                    /* n*b doubles							*/
} ;

typedef
struct amg_vector AMG_VECTOR;           /* a whole vector						*/

/****************************************************************************/
/******** functions for vectors                             *****************/
/****************************************************************************/

#define AMG_VECTOR_NAME(p)                      ((p)->name)
#define AMG_VECTOR_N(p)                         ((p)->n)
#define AMG_VECTOR_B(p)                         ((p)->b)
#define AMG_VECTOR_X(p)                         ((p)->x)
#define AMG_VECTOR_ENTRY(p,i,ii)        ((p)->x[(i)*(p)->b+(ii)])

AMG_VECTOR *AMG_NewVector (int n, int b, char *name);

/****************************************************************************/
/******** pattern of a general m x n sparse block matrix    *****************/
/****************************************************************************/

struct amg_matrix {                                     /* combines all links					*/
  char name[AMG_NAME_SIZE];                     /* name of this matrix					*/
  int n;                                                        /* dimension of matrix (always square)	*/
  int b;                                                        /* dimension of blocks: b x b			*/
  int bb;                                                       /* block size in doubles: ie bs=b*b		*/
  int system_as_scalar;                         /* system treated as scalars if > 1		*/
  int bandwidth;                                        /* bandwidth of the matrix				*/
  int nonzeros;                                         /* number of nonzero blocks allocated	*/
  int connections;                                      /* nonzeros actually used				*/
  int *ra;                                                      /* ra[i]: index of first entry of row i	*/
  int *ja;                                                      /* aj[k]: col index of entry k			*/
  double *a;                                                    /* the matrix							*/
} ;

typedef
struct amg_matrix AMG_MATRIX;                   /* a matrix                                                     */

/****************************************************************************/
/******** functions for vectors                             *****************/
/****************************************************************************/

#define AMG_MATRIX_NAME(p)                      ((p)->name)
#define AMG_MATRIX_N(p)                         ((p)->n)
#define AMG_MATRIX_B(p)                         ((p)->b)
#define AMG_MATRIX_BB(p)                        ((p)->bb)
#define AMG_MATRIX_SAS(p)                       ((p)->system_as_scalar)
#define AMG_MATRIX_NONZEROS(p)          ((p)->nonzeros)
#define AMG_MATRIX_CONNECTIONS(p)       ((p)->connections)
#define AMG_MATRIX_BW(p)                        ((p)->bandwidth)
#define AMG_MATRIX_RA(p)                        ((p)->ra)
#define AMG_MATRIX_JA(p)                        ((p)->ja)
#define AMG_MATRIX_A(p)                         ((p)->a)

/* Construction */
AMG_MATRIX *AMG_NewMatrix (int n, int b, int nonzeros, int system_as_scalar, char *name);
AMG_MATRIX *AMG_CopyMatrix (AMG_MATRIX *A, char *name);
int             AMG_SetRowLength (AMG_MATRIX *A, int i, int l);
int             AMG_FindEntry (AMG_MATRIX *A, int i, int j);
int             AMG_InsertEntry (AMG_MATRIX *A, int i, int j);
int             AMG_InsertValues (AMG_MATRIX *A, int i, int j, double *aij);
int             AMG_AddValues (AMG_MATRIX *A, int i, int j, double *aij);

/* input / output */
int       AMG_PrintVector   (int k, AMG_VECTOR **vlist, char *text);
int       AMG_PrintMatrix   (AMG_MATRIX *A, char *text);

#endif
