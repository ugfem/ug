// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>
#include <stdio.h>
#include "../main/defs.h"

/* Finds needed eigenvalues of tridiagonal T using either the QL
   algorithm or Sturm sequence bisection, whichever is predicted
   to be faster based on a simple complexity model. */
void get_ritzvals(alpha,beta,j,Anorm,workj,ritz,d,left_goodlim,right_goodlim,eigtol)
double *alpha;                  /* vector of Lanczos scalars */
double *beta;                   /* vector of Lanczos scalars */
int j;                          /* number of Lanczos iterations taken */
double Anorm;                   /* Gershgorin estimate */
double *workj;                  /* work vector for Sturm sequence */
double *ritz;                   /* array holding evals */
int d;                          /* problem dimension = num. eigenpairs needed */
int left_goodlim;               /* number of ritz pairs checked on left end */
int right_goodlim;              /* number of ritz pairs checked on right end */
double eigtol;                  /* tolerance on eigenpair */
{
  extern double BISECTION_SAFETY;       /* bisection tolerance function divided by this */
  extern int DEBUG_EVECS;               /* debug flag for eigen computation */
  int nvals_left;                       /* numb. evals to find on left end of spectrum */
  int nvals_right;                      /* numb. evals to find on right end of spectrum */
  double bisection_tol;                 /* width of interval bisection should converge to */
  double predicted_steps;               /* predicts number of required bisection steps */
  void ql(), shell_sort(), bisect();

  /* Determine number of ritzvals to find on left and right ends */
  nvals_left = max(d,left_goodlim);
  nvals_right = min(j-nvals_left,right_goodlim);

  /* Predict work for bisection and ql assuming bisection takes roughly 5j flops per
     step, ql takes roughly 30j^j flops per call. (Ignored sort following ql.) */
  bisection_tol = eigtol * eigtol / BISECTION_SAFETY;
  predicted_steps = (nvals_left + nvals_right)*
                    (log10(Anorm/bisection_tol)/log10(2.0));

  if (5*predicted_steps < 30*j*j) {
    bisect(alpha,beta-1,j,Anorm,workj,ritz,nvals_left,nvals_right,bisection_tol);
    if (DEBUG_EVECS > 2) {{char buf[150]; sprintf(buf,"  bisection\n");UserWrite(buf);}}
  }
  else {
    if (DEBUG_EVECS > 2) {{char buf[150]; sprintf(buf,"  ql\n");UserWrite(buf);}}
    ql(ritz,workj-1,j);
    shell_sort(j,ritz);
  }
}
