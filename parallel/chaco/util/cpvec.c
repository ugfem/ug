// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Copy a range of a double vector*/
void cpvec(copy,beg,end,vec)
double *copy;
int beg, end;
double *vec;
{
  int i;

  copy = copy + beg;
  vec = vec + beg;
  for (i = end - beg + 1; i; i--)  {
    (*copy++) = *vec++;
  }
}
