// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>
#include <stdio.h>
#include "../main/structs.h"
#include "../main/defs.h"


void make_subgeom(igeom, coords, subcoords, subnvtxs, loc2glob)
int igeom;                      /* 1, 2 or 3 dimensional geometry? */
float **coords;                 /* x, y and z coordinates of vertices */
float **subcoords;              /* x, y ans z coordinates in subgraph */
int subnvtxs;                   /* number of vertices in subgraph */
int *loc2glob;                  /* maps from subgraph to graph numbering */
{
  int i;                        /* loop counter */

  for (i=1; i<=subnvtxs; i++) {
    subcoords[0][i] = coords[0][loc2glob[i]];
    if (igeom > 1) subcoords[1][i] = coords[1][loc2glob[i]];
    if (igeom > 2) subcoords[2][i] = coords[2][loc2glob[i]];
  }
}
