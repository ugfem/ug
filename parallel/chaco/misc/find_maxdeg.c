// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

/* Find the maximum weighted degree of a vertex. */
double find_maxdeg(graph, nvtxs)
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vertices */
{
  double maxdeg;                /* maximum degree of a vertex */
  int using_ewgts;              /* are edge weights being used? */
  int i;                        /* loop counter */

  using_ewgts = (graph[1]->ewgts != NULL);

  /* Find the maximum weighted degree of a vertex. */
  maxdeg = 0;
  if (using_ewgts) {
    for (i=1; i<=nvtxs; i++) {
      if (-graph[i]->ewgts[0] > maxdeg) maxdeg = -graph[i]->ewgts[0];
    }
  }
  else {
    for (i=1; i<=nvtxs; i++) {
      if (graph[i]->nedges > maxdeg) maxdeg = graph[i]->nedges - 1;
    }
  }
  return(maxdeg);
}
