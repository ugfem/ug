// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>
#include <stdio.h>
#include "../main/structs.h"
#include "../main/defs.h"


void coordsout(coords, igeom, nvtxs, file)
float **coords;                 /* x, y and z coordinates of vertices */
int igeom;                      /* 1, 2 or 3 dimensional geometry? */
int nvtxs;                      /* number of vtxs in graph */
FILE *file;                     /* output file if not null */
{
  int i, j;                     /* loop counters */

  if (file == NULL) {
    for (i=1; i<=nvtxs; i++) {
      for (j=0; j<igeom; j++) {char buf[150]; sprintf(buf,"%.8g ", coords[j][i]);UserWrite(buf);}
      {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
    }
  }
  else {
    for (i=1; i<=nvtxs; i++) {
      for (j=0; j<igeom; j++) fprintf(file, "%.8g ", coords[j][i]);
      fprintf(file, "\n");
    }
  }
}
