// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Orthogonalize to all one's */
void orthog1(x,beg,end)
double *x;
int beg, end;
{
  int i;
  double *pntr;
  double sum;
  int len;

  len = end - beg + 1;
  sum = 0.0;
  pntr = x + beg;
  for (i = len; i; i--)  {
    sum += *pntr++;
  }
  sum /= len;
  pntr = x + beg;
  for (i = len; i; i--)  {
    *pntr++ -= sum;

  }
}
