// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


void makecewgts(graph, cgraph, cnvtxs, merged, reduction, nmerged,
                v2cv, seenflag, using_ewgts)
/* Compute the edge weights of coarsened graph. */
struct vtx_data **graph;        /* array of vtx data for graph */
struct vtx_data **cgraph;       /* coarsened version of graph */
int cnvtxs;                     /* number of vertices in coarsened graph */
struct ipairs *merged;          /* vtx pairs that get merged */
struct fpairs *reduction;       /* multipliers reducing edge weights */
int nmerged;                    /* number of vertices merged */
int *v2cv;                      /* mapping from vtxs to coarsened vtxs */
int using_ewgts;                /* are edge weights being used? */
short *seenflag;                /* flags for cnvtxs, initialized to zero */
{
  extern Heap   *heap;      /* pointer to heap of multigrid */
  extern double *MEM_OK;    /* variable for memory overflow exeception */
  extern double adj_cewgts_time;        /* time spent in adjcewgts */
  double time;                  /* timing value */
  float *eptr;                  /* loops through ewgts */
  double sum;                   /* sum of edge weights */
  int temp;                     /* temporary value */
  int neighbor;                 /* neighboring vertex in coarse graph */
  int i, j, k;                  /* loop counter */
  double drandom();
  double seconds();
  void adjcewgts();

  /* Randomly permute the merged edges to avoid ordering effects. */
  for (i=0; i<nmerged; i++) {
    j = nmerged*drandom();
    temp = merged[i].val1;
    merged[i].val1 = merged[j].val1;
    merged[j].val1 = temp;
    temp = merged[i].val2;
    merged[i].val2 = merged[j].val2;
    merged[j].val2 = temp;
  }
  /*pmerges(merged, nmerged);*/

  /* Initialize the edge weights of the coarsened graph */
  /* For each non-matched edge in original edge, add its weight to the */
  /* corresonding edge in coarsened graph. */
  /* Only work with one instance of edges, and copy to transpose later. */
  /* Note: this approach only works if later transformations to edge */
  /* weights work identically on any edges that overlap in coarse graph. */

  /* Now loop through the merged edges, adjusting neighboring edge weights. */

  time = seconds();
  adjcewgts(graph, cgraph, cnvtxs, merged, reduction, nmerged, v2cv, seenflag,
            using_ewgts);
  if (!MEM_OK) return;
  adj_cewgts_time += seconds() - time;

  /* Since only one instance of each edge has been handled, the weight */
  /* of its transpose must be set as well. */
  for (i=1; i<=cnvtxs; i++) {
    for (j=1; j<cgraph[i]->nedges; j++) {
      neighbor = cgraph[i]->edges[j];
      if (i < neighbor) {
        /* Find weight of edge neighbor->i. */
        for (k=1;
             k<cgraph[neighbor]->nedges && cgraph[neighbor]->edges[k]!=i;
             k++) ;
        cgraph[i]->ewgts[j] = cgraph[neighbor]->ewgts[k];
      }
    }
  }


  /* Finally, set the matrix diagonals correctly. */
  for (i=1; i<=cnvtxs; i++) {
    sum = 0;
    eptr = cgraph[i]->ewgts;
    for (j=cgraph[i]->nedges-1; j; j--) {
      sum += *(++eptr);
    }
    cgraph[i]->ewgts[0] = -sum;
  }
}
