// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_iter.c													*/
/*																			*/
/* Purpose:   Simple iterative schemes for amg package						*/
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
/* History:   04 FEB 1996 Begin												*/
/*			  01 OKT 1997 redesign											*/
/*			  21 OKT 1997 EX code due to Klaus Johannsen					*/
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
#include <math.h>

#include "amg_header.h"
#include "amg_low.h"
#include "amg_sp.h"
#include "amg_blas.h"
#include "amg_iter.h"

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
   AMG_jac - jacobi kernel

   SYNOPSIS:
   int AMG_jac (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d)

   PARAMETERS:
   .  A - matrix
   .  d - defect
   .  v - correction

   DESCRIPTION:
   Solves the system Dv=d where D=diag(A).

   RETURN VALUE:
   .n AMG_OK
   .n AMG_FATAL

   D*/
/****************************************************************************/

int AMG_jac (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i;
  register double *v, *d, *a;
  register int *ra;
  register int b,bb;
  double om;

  /* plausi */
  if (AMG_VECTOR_N(v_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_N(d_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(v_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(d_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);

  /* prepare data */
  n  = AMG_MATRIX_N(A);
  b  = AMG_MATRIX_B(A);
  bb = AMG_MATRIX_BB(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);

  /* loop */
  switch (b)
  {
  case 1 :
    om = omega[0];
    for (i=0; i<n; i++) v[i] = om*d[i]/a[ra[i]];
    break;

  default :
    AMG_Print("jac: blocksize>1 not implemented yet\n");
    break;
  }


  return(AMG_FATAL);
}



/****************************************************************************/
/*D
   AMG_sorf - SOR kernels forward and backward

   SYNOPSIS:
   int AMG_sorf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
   int AMG_sorb (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)

   PARAMETERS:
   .  A - matrix
   .  d - defect
   .  v - correction
   .  omega - damping factor

   DESCRIPTION:
   Solve the system Lv=d or Uv=d with SOR damping where L is the
   lower triangle of A and U the upper triangle of A. The damping factor
   is an array containing one value per component.

   RETURN VALUE:
   .n AMG_OK
   .n AMG_FATAL

   D*/
/****************************************************************************/

int AMG_sorf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,end,start;
  register double *v, *d, *a;
  register int *ra, *ja;
  register int b,bb;
  register double s,om;

  /* plausi */
  if (AMG_VECTOR_N(v_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_N(d_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(v_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(d_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);

  /* prepare data */
  n  = AMG_MATRIX_N(A);
  b  = AMG_MATRIX_B(A);
  bb = AMG_MATRIX_BB(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* loop */
  switch (b)
  {
  case 1 :
    om=omega[0];
    for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      s = 0.0;
      for (k=start+1; k<end; k++)
        if (ja[k]<i) s += a[k]*d[ja[k]];
      v[i] = om*(d[i]-s)/a[ra[i]];
    }
    break;

  default :
    AMG_Print("sor: blocksize>1 not implemented yet\n");
    break;
  }


  return(AMG_FATAL);
}

int AMG_sorb (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end;
  register double *v, *d, *a;
  register int *ra, *ja;
  register int b,bb;
  register double s,om;

  /* plausi */
  if (AMG_VECTOR_N(v_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_N(d_)!=AMG_MATRIX_N(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(v_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
  if (AMG_VECTOR_B(d_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);

  /* prepare data */
  n  = AMG_MATRIX_N(A);
  b  = AMG_MATRIX_B(A);
  bb = AMG_MATRIX_BB(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* loop */
  switch (b)
  {
  case 1 :
    om=omega[0];
    for (i=n-1; i>=0; i--)
    {
      start = ra[i]; end = start+ja[start];
      s = 0.0;
      for (k=start+1; k<end; k++)
        if (ja[k]>i) s += a[k]*d[ja[k]];
      v[i] = om*(d[i]-s)/a[ra[i]];
    }
    break;

  default :
    AMG_Print("sor: blocksize>1 not implemented yet\n");
    break;
  }


  return(AMG_FATAL);
}

/****************************************************************************/
/*D
   AMG_EXDecomposeMatrixdouble - LU decompose a band matrix (double numbers)

   SYNOPSIS:
   int AMG_EXDecomposeMatrixdouble (double *Mat, int bw, int n);

   PARAMETERS:
   .  Mat - pointer to double array containing the bandmatrix
   .  bw - bandwidth
   .  n - number of rows (==columns) of the matrix

   DESCRIPTION:
   This function calculates the Gauss decomposition of the given band matrix;
   the L and U factors are stored instead of the band matrix.

   RETURN VALUE:
   int  0: o.k.
        1: main diagonal element to small (0.0); can not divide

   SEE ALSO:
   EXDecomposeMatrixFLOAT, EXCopyMatrixdouble, EXApplyLUdouble
   D*/
/****************************************************************************/

int AMG_EXDecomposeMatrix (double *Mat, int bw, int n)
{
  int i,j,k;
  double f,d;

  for (i=0; i<n-1; i++)
  {
    d = AMG_EX_MAT(Mat,bw,i,i);
    if (AMG_ABS(d)<=1.0E-80) return (1);
    for (j=i+1; j<=AMG_MIN(i+bw,n-1); j++)
    {
      f = AMG_EX_MAT(Mat,bw,j,i)/d;
      AMG_EX_MAT(Mat,bw,j,i) = f;
      for (k=i+1; k<=AMG_MIN(i+bw,n-1); k++)
        AMG_EX_MAT(Mat,bw,j,k) -= f*AMG_EX_MAT(Mat,bw,i,k);
    }
  }
  return (0);
}

/****************************************************************************/
/*D
   AMG_EXApplyLUdouble - applies a LU decomposed band matrix (double numbers)

   SYNOPSIS:
   int AMG_EXApplyLUdouble (double *Mat, int bw, int n, double *Vec);

   PARAMETERS:
   .  Mat - pointer to double array containing the bandmatrix
   .  bw - bandwidth
   .  n - number of rows (==columns) of the matrix
   .  Vec - pointer to double array containing the vector

   DESCRIPTION:
   This function solves for the LU decomposed band matrix 'Mat' the equation
   L*U x = f.
   f is provided in 'Vec' and the result x is returned again in 'Vec'.

   RETURN VALUE:
   int  0: o.k.

   SEE ALSO:
   EXApplyLUFLOAT, EXCopyMatrixdouble, EXDecomposeMatrixdouble
   D*/
/****************************************************************************/

int AMG_EXApplyLU (double *Mat, int bw, int n, double *Vec)
{
  int i,j;

  /* invert lower */
  for (i=1; i<n; i++)
    for (j=AMG_MAX(i-bw,0); j<i; j++)
      Vec[i] -= AMG_EX_MAT(Mat,bw,i,j)*Vec[j];

  /* invert upper */
  for (i=n-1; i>=0; i--)
  {
    for (j=i+1; j<=AMG_MIN(i+bw,n-1); j++)
      Vec[i] -= AMG_EX_MAT(Mat,bw,i,j)*Vec[j];
    Vec[i] /= AMG_EX_MAT(Mat,bw,i,i);
  }
  return (0);
}
