// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"


void make_maps(assignment, nvtxs, set, glob2loc, loc2glob, ready)
short *assignment;              /* set assignments for graph */
int nvtxs;                      /* length of assignment */
int set;                        /* set value denoting subgraph */
int *glob2loc;                  /* graph -> subgraph numbering map */
int *loc2glob;                  /* subgraph -> graph numbering map */
short *ready;
{
  int i, j;                     /* loop counter */

  j = 1;
  for (i=1; i<=nvtxs; i++) {
    if ((assignment[i] == set) && (ready[i]==FALSE)) {
      if (glob2loc != NULL) glob2loc[i] = j;
      loc2glob[j] = i;
      j++;
    }
  }
}
