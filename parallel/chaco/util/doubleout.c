// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include  <stdio.h>

/* Print a double precesion number with filtering format. */
void doubleout(number,mode)
double number;
int mode;
{
  double fabs();                /* intrinsic absolute value function */

  if (mode == 1) {
    if (fabs(number) < 100) {
      {char buf[150]; sprintf(buf,"  %19.16f",number);UserWrite(buf);}
    }
    else {
      {char buf[150]; sprintf(buf,"  %19g",number);UserWrite(buf);}
    }
  }
}
