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
  register int n,i,k,start,end;
  register double *v, *d, *a, *xx, *aa, *yy;
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
  register int n,i,k,start,end;
  register double *v, *d, *a, *xx, *aa, *yy;
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
  register double *v, *d, *a, *xx, *aa, *yy;
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
