// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

void orthogvec(vec1, beg, end, vec2)
double *vec1;           /* vector to be orthogonalized */
int beg, end;           /* start and stop range for vector */
double *vec2;           /* vector to be orthogonalized against */
{
  double alpha;
  double dot();
  void scadd();

  alpha = -dot(vec1, beg, end, vec2)/dot(vec2, beg, end, vec2);
  scadd(vec1, beg, end, alpha, vec2);
}
