// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"

/* Finds some eigenvalues of T using Sturm sequence bisection.
   Based on Wilkinson's algorithm, AEP, p.302. Beta must be offset
   as it is passed in so that indices are consistent with AEP. */
void bisect(alpha,beta,j,Anorm,workj,ritz,nevals_left, nevals_right,tol)
double *alpha;                  /* vector of Lanczos scalars */
double *beta;                   /* vector of Lanczos scalars */
int j;                          /* number of Lanczos iterations taken */
double Anorm;                   /* Gershgorin estimate */
double *workj;                  /* work vector for Sturm sequence */
int nevals_left;                /* number of evals on right to find */
int nevals_right;               /* number of evals on left to find */
double *ritz;                   /* array holding evals */
double tol;                     /* tolerance on bracket width */
{
  extern int DEBUG_EVECS;       /* debug flag for eigen computation */
  extern double DOUBLE_MAX;             /* largest double value */
  int index;                    /* index of sturm polynomial */
  int i;                        /* loop index */
  double *pntr;                 /* pntr to double array */
  double x1,x2;                 /* the bracketing interval */
  int x1cnt, x2cnt;             /* Sturm counts at x1 and x2 */
  double x;                     /* the inserted point */
  int xcnt;                     /* the Sturm count at x */
  int steps;                    /* number of bisection steps for a Ritzval */
  int tot_steps;                /* number of bisection steps for all Ritzvals */
  int numbracketed;             /* number of evals between x1 and x2 */
  int sturmcnt();

  /* Initialize portion of ritz we will use (use max double so scanmin
     will work properly when called later on) */
  pntr = &ritz[1];
  for (i=j; i; i--) {
    *pntr++ = DOUBLE_MAX;
  }

  tot_steps = 0;

  /* find evals on left in decreasing index order */
  x2 = Anorm;
  x2cnt = j;
  numbracketed = j;
  for (index = nevals_left; index >= 1; index--) {
    x1 = 0;
    x1cnt = 0;
    steps = 0;
    while ((x2 - x1) > tol || numbracketed > 1) {
      x = 0.5*(x1 + x2);
      xcnt = sturmcnt(alpha,beta,j,x,workj);
      if (xcnt >= index) {
        x2 = x;
        x2cnt = xcnt;
      }
      else {
        x1 = x;
        x1cnt = xcnt;
      }
      numbracketed = x2cnt - x1cnt;
      steps++;
    }
    ritz[index] = 0.5*(x1 + x2);
    if (DEBUG_EVECS > 2) {
      {char buf[150]; sprintf(buf,"index %d, bisection steps %d, root %20.16f\n", index,steps,ritz[index]);UserWrite(buf);}
    }
    tot_steps += steps;
  }

  /* find evals on right in increasing index order */
  x1 = 0;
  x1cnt = 0;
  for (index = j - nevals_right + 1; index <= j; index++) {
    x2 = Anorm;
    x2cnt = j;
    steps = 0;
    while ((x2 - x1) > tol || numbracketed > 1) {
      x = 0.5*(x1 + x2);
      xcnt = sturmcnt(alpha,beta,j,x,workj);
      if (xcnt >= index) {
        x2 = x;
        x2cnt = xcnt;
      }
      else {
        x1 = x;
        x1cnt = xcnt;
      }
      numbracketed = x2cnt - x1cnt;
      steps++;
    }
    ritz[index] = 0.5*(x1 + x2);
    if (DEBUG_EVECS > 2) {
      {char buf[150]; sprintf(buf,"index %d, bisection steps %d, root %20.16f\n", index,steps,ritz[index]);UserWrite(buf);}
    }
    tot_steps += steps;
  }
  if (DEBUG_EVECS > 2) {
    {char buf[150]; sprintf(buf,"Total number of bisection steps %d\n",tot_steps);UserWrite(buf);}
  }
}
