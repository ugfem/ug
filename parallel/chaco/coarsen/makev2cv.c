// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


void makev2cv(mflag, nvtxs, v2cv)
/* Construct mapping from original graph vtxs to coarsened graph vtxs. */
int *mflag;                     /* flag indicating vtx selected or not */
int nvtxs;                      /* number of vtxs in original graph */
int *v2cv;                      /* mapping from vtxs to coarsened vtxs */
{
  int i, j;                     /* loop counters */

  j = 1;
  for (i=1; i<=nvtxs; i++) {
    if (mflag[i] == 0 || mflag[i] > i) v2cv[i] = j++;
    else v2cv[i] = v2cv[mflag[i]];
  }
}


void makemflag(merged, nmerged, mflag, nvtxs)
struct ipairs *merged;          /* vertex pairs that got merged */
int nmerged;                    /* number of merged pairs */
int *mflag;                     /* vector denoting matching vertices */
int nvtxs;                      /* total number of vertices */
{
  int *iptr;                    /* loops through mflag */
  int i;                        /* loop counter */

  for (iptr=mflag, i=nvtxs; i; i--) *(++iptr) = 0;

  for (i=0; i<nmerged; i++) {
    mflag[merged[i].val1] = merged[i].val2;
    mflag[merged[i].val2] = merged[i].val1;
  }
}
