// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Finds eigenvector s of T and returns residual norm. */
double Tevec(alpha,beta,j,ritz,s)
double *alpha;                          /* vector of Lanczos scalars */
double *beta;                           /* vector of Lanczos scalars */
int j;                                  /* number of Lanczos iterations taken */
double ritz;                            /* approximate eigenvalue  of T */
double *s;                              /* approximate eigenvector of T */
{
  int i;                                /* index */
  double residual;                      /* how well recurrence gives eigenvector */

  double normalize();                   /* normalizes vector, returns norm */
  double fabs();                        /* intrinsic absolute value */

  s[j] = 1.0;

  if (j == 1)  {
    residual = (alpha[1]-ritz);
  }

  if (j >= 2)  {
    s[j-1] = - (alpha[j]-ritz)/beta[j];
    for (i=j; i>=3; i--)  {
      s[i-2] = -((alpha[i-1]-ritz)*s[i-1] + beta[i]*s[i])/beta[i-1];
    }
    residual = (alpha[1]-ritz)*s[1] + beta[2]*s[2];
  }

  residual = fabs(residual)/normalize(s,1,j);
  return(residual);
}
