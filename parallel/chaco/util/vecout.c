// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        <math.h>

/* Print vertically range of double vector. */
void vecout(vec,beg,end,tag)
double *vec;
int beg,end;
char *tag;
{
  int i;

  {char buf[150]; sprintf(buf,"%s:\n", tag);UserWrite(buf);}
  for(i = beg; i <= end; i++) {
    if (fabs(vec[i]) >= 1.0e-16)
    {char buf[150]; sprintf(buf,"%2d.   %24.16f\n", i, vec[i]);UserWrite(buf);}
    else
    {char buf[150]; sprintf(buf,"%2d.         %g \n", i, vec[i]);UserWrite(buf);}

  }
}
