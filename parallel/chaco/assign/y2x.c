// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/defs.h"


void x2y(xvecs, ndims, nmyvtxs, wsqrt)
double **xvecs;                 /* pointer to list of x-vectors */
int ndims;                      /* number of divisions to make (# xvecs) */
int nmyvtxs;                    /* number of vertices I own (lenght xvecs) */
double *wsqrt;                  /* sqrt of vertex weights */
/* Convert from x to y by dividing by wsqrt. */
{
  double *wptr;                 /* loops through wsqrt */
  double *xptr;                 /* loops through elements of a xvec */
  int i, j;                     /* loop counters */

  if (wsqrt == NULL) return;

  for (i=1; i<=ndims; i++) {
    xptr = xvecs[i];
    wptr = wsqrt+1;
    for (j=nmyvtxs; j; j--) {
      *xptr++ /= *wptr++;
    }
  }
}


void y2x(yvecs, ndims, nmyvtxs, wsqrt)
double **yvecs;                 /* pointer to list of y-vectors */
int ndims;                      /* number of divisions to make (# yvecs) */
int nmyvtxs;                    /* number of vertices I own (lenght yvecs) */
double *wsqrt;                  /* sqrt of vertex weights */
/* Convert from y to x by multiplying by wsqrt. */
{
  double *wptr;                 /* loops through wsqrt */
  double *yptr;                 /* loops through elements of a yvec */
  int i, j;                     /* loop counters */

  if (wsqrt == NULL) return;

  for (i=1; i<=ndims; i++) {
    yptr = yvecs[i]+1;
    wptr = wsqrt+1;
    for (j=nmyvtxs; j; j--) {
      /*
               printf("y2x(): i= %d j=%d *yptr=%10.6f *wptr=%10.6f\n",i,j,*yptr,*wptr);
       */
      *yptr++ *= *wptr++;
    }
  }
}
