// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Eigensolution of real symmetric tridiagonal matrix. */
/* After Numerical Recipies p. 380, removed eigenvector calc. */
/* CAREFUL: e is used as workspace, d returns eigenvals. */

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

void ql(d,e,n)
double d[],e[];
int n;

{
  int m,l,iter,i;
  double s,r,p,g,f,dd,c,b;
  void nrerror();
  double fabs();
  double sqrt();

  /* off diagonal data starts at second vector location, so renumber */
  for (i = 2; i <= n; i++) e[i-1] = e[i];
  e[n] = 0.0;

  for (l = 1; l <= n; l++)  {
    iter = 0;
    do {
      for (m = l; m <= n-1; m++)  {
        dd = fabs(d[m]) + fabs(d[m+1]);
        if (fabs(e[m]) + dd == dd) break;
      }
      if (m != l)  {
        if (iter++ == 50) nrerror("ql: Too many iterations (max 50).");
        g = (d[l+1] - d[l])/(2.0 * e[l]);
        r = sqrt((g*g) + 1.0);
        g = d[m] - d[l] + e[l]/(g + SIGN(r,g));
        s = c = 1.0;
        p = 0.0;
        for (i = m-1; i >= l; i--)  {
          f = s * e[i];
          b = c * e[i];
          if (fabs(f) >= fabs(g))  {
            c = g/f;
            r = sqrt((c*c) + 1.0);
            e[i+1] = f * r;
            c *= (s = 1.0/r);
          }
          else  {
            s = f/g;
            r = sqrt((s*s) + 1.0);
            e[i+1] = g * r;
            s *= (c = 1.0/r);
          }
          g = d[i+1] - p;
          r = (d[i] - g) * s + 2.0*c*b;
          p = s * r;
          d[i+1] = g + p;
          g = c*r - b;
        }
        d[l] = d[l] - p;
        e[l] = g;
        e[m] = 0.0;
      }
    }
    while (m != l);

  }
}
