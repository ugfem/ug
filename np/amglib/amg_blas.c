// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_blas.c													*/
/*																			*/
/* Purpose:   BLAS (Basic Linear Algebra Subroutines) for amg package		*/
/*			  Contains Level 1 and 2										*/
/*																			*/
/* Author:	  Peter Bastian					                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   31 JAN 1996 Begin												*/
/*			  01 OKT 1997 redesign											*/
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
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "amg_header.h"
#include "amg_low.h"
#include "amg_sp.h"
#include "amg_blas.h"

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
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   blas - BLAS routines based on the amg data structure

   SYNOPSIS:
   int    AMG_dset      (AMG_VECTOR *x, double a                        );
   int    AMG_randomize (AMG_VECTOR *x						                );
   int    AMG_dcopy     (AMG_VECTOR *x, AMG_VECTOR *y               );
   int    AMG_dscale    (AMG_VECTOR *x, double a                        );
   int    AMG_daxpy     (AMG_VECTOR *x, double a,      AMG_VECTOR *y);
   double AMG_ddot      (AMG_VECTOR *x, AMG_VECTOR *y               );
   int    AMG_dmatset   (AMG_MATRIX *A, double a                        );
   int    AMG_dmatcopy  (AMG_MATRIX *A, AMG_MATRIX *B               );
   int    AMG_dmatmul   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_VECTOR *y);
   int    AMG_dmatminus (AMG_VECTOR *x, AMG_MATRIX *A, AMG_VECTOR *y);


   PARAMETERS:
   .  x,y - vectors
   .  A,B - sparse matrices
   .  a - scalar

   DESCRIPTION:
   .n AMG_dset sets all values of vector x to a.
   .n AMG_randomize sets random values in vector x.
   .n AMG_dcopy copies vector x to vector y.
   .n AMG_dscale multiplies each entry of x with a.
   .n AMG_daxpy computes x = x + ay
   .n AMG_ddot computes euclidean scalar product of two vectors
   .n AMG_dmatset sets al entries of matrix A to a.
   .n AMG_dmatcopy copies matrix B to A (same structure assumed)
   .n AMG_dmatmul computes x = Ay
   .n AMG_dmatminus computes x = x - Ay

   All routines check the compatibility of their arguments and return
   errors if necessary. Where appropriate the routines have been
   implemented in a scalar and a block version to enhance efficiency
   in the scalar case.

   RETURN VALUE:
   .n AMG_OK
   .n AMG_FATAL

   D*/
/****************************************************************************/

int AMG_dset (AMG_VECTOR *x, double a)
{
  register int i,n;
  register double *values;

  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  values = AMG_VECTOR_X(x);
  for (i=0; i<n; i++) *values++ = a;

  return(AMG_OK);
}

int AMG_randomize (AMG_VECTOR *x)
{
  register int i,n;
  register double *values;

  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  values = AMG_VECTOR_X(x);
  for (i=0; i<n; i++) *values++ = rand();

  return(AMG_OK);
}

int AMG_dcopy (AMG_VECTOR *x, AMG_VECTOR *y)
{
  register int i,n;
  register double *values_x, *values_y;

  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(y)) return(AMG_FATAL);
  if (AMG_VECTOR_B(x)!=AMG_VECTOR_B(y)) return(AMG_FATAL);

  values_x = AMG_VECTOR_X(x);
  values_y = AMG_VECTOR_X(y);
  for (i=0; i<n; i++) *values_x++ = *values_y++;

  return(AMG_OK);
}

int AMG_dscale (AMG_VECTOR *x, double a)
{
  register int i,n;
  register double *values;

  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  values = AMG_VECTOR_X(x);
  for (i=0; i<n; i++) *values++ *= a;

  return(AMG_OK);
}

int AMG_daxpy (AMG_VECTOR *x, double a, AMG_VECTOR *y)
{
  register int i,n;
  register double *values_x, *values_y;

  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(y)) return(AMG_FATAL);
  if (AMG_VECTOR_B(x)!=AMG_VECTOR_B(y)) return(AMG_FATAL);

  values_x = AMG_VECTOR_X(x);
  values_y = AMG_VECTOR_X(y);
  for (i=0; i<n; i++) *values_x++ += a * (*values_y++);

  return(AMG_OK);
}

double AMG_ddot (AMG_VECTOR *x, AMG_VECTOR *y)
{
  register int i,n;
  register double *values_x, *values_y;
  register double s=0.0;

  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(y)) return(AMG_FATAL);
  if (AMG_VECTOR_B(x)!=AMG_VECTOR_B(y)) return(AMG_FATAL);

  values_x = AMG_VECTOR_X(x);
  values_y = AMG_VECTOR_X(y);
  for (i=0; i<n; i++) s += (*values_x++) * (*values_y++);

  return(s);
}



int AMG_dmatset (AMG_MATRIX *A, double a)
{
  int size,i;
  double *values;

  size = AMG_MATRIX_N(A)*AMG_MATRIX_BB(A);
  values = AMG_MATRIX_A(A);
  for (i=0; i<size; i++) *values++ = a;

  return(AMG_OK);
}


int AMG_dmatcopy (AMG_MATRIX *A, AMG_MATRIX *B)
{
  int size_a,size_b,i;
  double *values_a, *values_b;

  size_a = AMG_MATRIX_N(A)*AMG_MATRIX_BB(A);
  size_b = AMG_MATRIX_N(B)*AMG_MATRIX_BB(B);
  if (size_a!=size_b) return(AMG_FATAL);

  values_a = AMG_MATRIX_A(A);
  values_b = AMG_MATRIX_A(B);
  for (i=0; i<size_a; i++) *values_a++ = *values_b++;

  return(AMG_OK);
}

#define ZERO2(x)    x[0] = x[1] = 0.0;
#define ZERO3(x)    x[0] = x[1] = x[2] = 0.0;
#define ZERO4(x)    x[0] = x[1] = x[2] = x[3] = 0.0;


#define XAY2(x,a,y)     {x[0] += a[0]*y[0] + a[1]*y[1]; \
                         x[1] += a[2]*y[0] + a[3]*y[1];}

#define XAY3(x,a,y)     {x[0] += a[0]*y[0] + a[1]*y[1] + a[2]*y[2]; \
                         x[1] += a[3]*y[0] + a[4]*y[1] + a[5]*y[2]; \
                         x[2] += a[6]*y[0] + a[7]*y[1] + a[8]*y[2];}

#define XAY4(x,a,y)     {x[0] += a[0]*y[0] + a[1]*y[1] + a[2]*y[2] + a[3]*y[3]; \
                         x[1] += a[4]*y[0] + a[5]*y[1] + a[6]*y[2] + a[7]*y[3]; \
                         x[2] += a[8]*y[0] + a[9]*y[1] + a[10]*y[2] + a[11]*y[3]; \
                         x[3] += a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];} \

int AMG_dmatmul (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_VECTOR *y_)
{
  register int n,i,k,start,end;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;

  /* plausi */
  if (AMG_VECTOR_N(x_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_N(y_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(x_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(y_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);

  /* prepare data */
  n = AMG_VECTOR_N(x_);
  b = AMG_VECTOR_B(x_);
  bb = AMG_MATRIX_BB(A);
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* loop */
  switch (b)
  {
  case 1 :
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      s = a[start]*y[i];
      for (k=start+1; k<end; k++) s += a[k]*y[ja[k]];
      x[i] = s;
    }
    break;

  case 2 :
    xx = x; aa = a;
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      ZERO2(xx);
      yy = y+(i*b); XAY2(xx,aa,yy); aa+=bb;
      for (k=start+1; k<end; k++)
      {
        yy = y+(ja[k]*b); XAY2(xx,aa,yy); aa+=bb;
      }
      xx+=b;
    }
    break;

  case 3 :
    xx = x; aa = a;
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      ZERO3(xx);
      yy = y+(i*b); XAY3(xx,aa,yy); aa+=bb;
      for (k=start+1; k<end; k++)
      {
        yy = y+(ja[k]*b); XAY3(xx,aa,yy); aa+=bb;
      }
      xx+=b;
    }
    break;

  case 4 :
    xx = x; aa = a;
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      ZERO4(xx);
      yy = y+(i*b); XAY4(xx,aa,yy); aa+=bb;
      for (k=start+1; k<end; k++)
      {
        yy = y+(ja[k]*b); XAY4(xx,aa,yy); aa+=bb;
      }
      xx+=b;
    }
    break;

  default :
    AMG_Print("dmatmul: blocksize>4 not implemented yet\n");
    break;
  }

  return(AMG_OK);
}

#define XAYM2(x,a,y)   {x[0] -= a[0]*y[0] + a[1]*y[1]; \
                        x[1] -= a[2]*y[0] + a[3]*y[1];}

#define XAYM3(x,a,y)   {x[0] -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2]; \
                        x[1] -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2]; \
                        x[2] -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];}

#define XAYM4(x,a,y)   {x[0] -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2] + a[3]*y[3]; \
                        x[1] -= a[4]*y[0] + a[5]*y[1] + a[6]*y[2] + a[7]*y[3]; \
                        x[2] -= a[8]*y[0] + a[9]*y[1] + a[10]*y[2] + a[11]*y[3]; \
                        x[3] -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];} \


int AMG_dmatminus (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_VECTOR *y_)
{
  register int n,i,k,start,end;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;

  /* plausi */
  if (AMG_VECTOR_N(x_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_N(y_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(x_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(y_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);

  /* prepare data */
  n = AMG_VECTOR_N(x_);
  b = AMG_VECTOR_B(x_);
  bb = AMG_MATRIX_BB(A);
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* loop */
  switch (b)
  {
  case 1 :
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      s = a[start]*y[i];
      for (k=start+1; k<end; k++) s += a[k]*y[ja[k]];
      x[i] -= s;
    }
    break;

  case 2 :
    xx = x; aa = a;
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      ZERO2(xx);
      yy = y+(i*b); XAYM2(xx,aa,yy); aa+=bb;
      for (k=start+1; k<end; k++)
      {
        yy = y+(ja[k]*b); XAYM2(xx,aa,yy); aa+=bb;
      }
      xx+=b;
    }
    break;

  case 3 :
    xx = x; aa = a;
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      ZERO3(xx);
      yy = y+(i*b); XAYM3(xx,aa,yy); aa+=bb;
      for (k=start+1; k<end; k++)
      {
        yy = y+(ja[k]*b); XAYM3(xx,aa,yy); aa+=bb;
      }
      xx+=b;
    }
    break;

  case 4 :
    xx = x; aa = a;
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      ZERO4(xx);
      yy = y+(i*b); XAYM4(xx,aa,yy); aa+=bb;
      for (k=start+1; k<end; k++)
      {
        yy = y+(ja[k]*b); XAYM4(xx,aa,yy); aa+=bb;
      }
      xx+=b;
    }
    break;

  default :
    AMG_Print("dmatmul: blocksize>4 not implemented yet\n");
    break;
  }

  return(AMG_OK);
}
