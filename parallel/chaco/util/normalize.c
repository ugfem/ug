// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Normalizes an n-vector over range. */
double normalize(vec,beg,end)
double *vec;
int beg, end;
{
  int i;
  double scale;
  double norm();

  scale = norm(vec,beg,end);
  vec = vec + beg;
  for (i = end - beg + 1; i; i--)  {
    *vec = *vec/scale;
    vec++;
  }
  return(scale);
}
