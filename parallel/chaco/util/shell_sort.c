// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>

#define ALN2I 1.442695022
#define TINY 1.0e-5

/* Sort a double array using Shell's method. Modified algorithm
   from p. 245, Numerical Recipies (changed float to double). */
void shell_sort(n,arr)
double arr[];
int n;

{
  int nn,m,j,i,lognb2;
  double t;

  lognb2 = (log((double) n)*ALN2I+TINY);
  m = n;
  for (nn=1; nn<=lognb2; nn++) {
    m >>= 1;
    for (j=m+1; j<=n; j++) {
      i = j-m;
      t = arr[j];
      while (i >= 1 && arr[i] > t) {
        arr[i+m] = arr[i];
        i -= m;
      }
      arr[i+m] = t;
    }
  }
}
