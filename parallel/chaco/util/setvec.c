// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Set a double precision vector to constant over range. */
void setvec(vec,beg,end,setval)
double vec[];
int beg, end;
double setval;
{
  int i;

  vec = vec + beg;
  for(i = end - beg + 1; i; i--)   {
    (*vec++) = setval;
  }
}
