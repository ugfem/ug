// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Evaluates principal minor polynomial */
/* Based on Wilkinson's algorithm, AEP, p.302 */
int sturmcnt(alpha,beta,j,mu,p)
double *alpha;                  /* vector of Lanczos scalars */
double *beta;                   /* vector of Lanczos scalars */
int j;                          /* index of Lanczos step we're on */
double mu;                      /* argument to the sequence generating polynomial */
double *p;                      /* work vector for sturm sequence */
{
  int i;                        /* loop index */
  int cnt;                      /* number of sign changes in the sequence */
  int sign;                     /* algebraic sign of current sequence value */
  int last_sign;                /* algebraic sign of previous sequence value */
  double *pntr_p;               /* for stepping through sequence array */
  double *pntr_p1;              /* for stepping through sequence array one behind */
  double *pntr_p2;              /* for stepping through sequence array two behind */
  double *pntr_alpha;           /* for stepping through alpha array */
  double *pntr_beta;            /* for stepping through beta array */

  if (j == 1) {
    /* have to special case this one */
    if (alpha[1] > mu) {cnt = 1;}
    else {cnt = 0;}
  }
  else {
    /* compute the Sturm sequence */
    p[0] = 1;
    p[1] = alpha[1] - mu;
    pntr_p = &p[2];
    pntr_p1 = &p[1];
    pntr_p2 = &p[0];
    pntr_alpha = &alpha[2];
    pntr_beta = &beta[2];
    for (i=2; i<=j; i++) {
      *pntr_p++ = (*pntr_alpha - mu)*(*pntr_p1) -
                  (*pntr_beta)*(*pntr_beta)*(*pntr_p2);
      pntr_alpha++;
      pntr_beta++;
      pntr_p1++;
      pntr_p2++;
    }
    /* count sign changes */
    cnt = 0;
    last_sign = 1;
    pntr_p = &p[1];
    for (i=j; i; i--) {
      if (*pntr_p > 0) {sign = 1;}
      else if (*pntr_p < 0) {sign = -1;}
      else {sign = -last_sign;}
      if (sign == last_sign) {cnt++;}
      last_sign = sign;
      pntr_p++;
    }
  }

  /* cnt is number of evals strictly > than mu. Instead
     send back index of the eigenvalue just left of mu. */
  return(j-cnt);
}
