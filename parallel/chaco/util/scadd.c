// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Scaled add - fills vec1 with vec1 + alpha*vec2 over range*/
void scadd(vec1,beg,end,fac,vec2)
double *vec1;
int beg, end;
double fac;
double *vec2;
{
  int i;

  vec1 = vec1 + beg;
  vec2 = vec2 + beg;
  for (i = end - beg + 1; i; i--)  {
    (*vec1++) += fac * (*vec2++);
  }
}
