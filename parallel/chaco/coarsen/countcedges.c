// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


void countcedges(graph, nvtxs, start, seenflag, mflag, v2cv, pcnedges)
/* Count edges in coarsened graph. */
struct vtx_data **graph;        /* array of vtx data for graph */
int nvtxs;                      /* number of vertices in graph */
int *start;                     /* start of edgevals list for each vertex */
short *seenflag;                /* flags indicating vtxs already counted */
int *mflag;                     /* flag indicating vtx matched or not */
int *v2cv;                      /* mapping from fine to coarse vertices */
int *pcnedges;                  /* number of edges in coarsened graph */
{
  int *jptr;                    /* loops through edge list */
  int cnedges;                  /* twice number of edges in coarsened graph */
  int neighbor;                 /* neighboring vertex */
  int cneighbor;                /* neighboring vertex in coarse graph */
  int nneighbors;               /* number of neighboring vertices */
  int newi;                     /* loops over vtxs in coarsened graph */
  int i, j;                     /* loop counters */

  cnedges = 0;
  newi = 1;
  start[1] = 0;
  for (i=1; i<=nvtxs; i++) {
    if (mflag[i] == 0 || mflag[i] > i) {
      nneighbors = 0;
      jptr = graph[i]->edges;
      for (j=graph[i]->nedges-1; j; j--) {
        /* Has edge already been added? */
        neighbor = *(++jptr);
        if (neighbor != mflag[i]) {
          cneighbor = v2cv[neighbor];

          if (!seenflag[cneighbor]) {
            nneighbors++;
            seenflag[cneighbor] = 1;
          }
        }
      }
    }
    if (mflag[i] > i) {
      jptr = graph[mflag[i]]->edges;
      for (j=graph[mflag[i]]->nedges-1; j; j--) {
        neighbor = *(++jptr);
        if (neighbor != i) {
          cneighbor = v2cv[neighbor];

          if (!seenflag[cneighbor]) {
            nneighbors++;
            seenflag[cneighbor] = 1;
          }
        }
      }
    }
    if (mflag[i] == 0 || mflag[i] > i) {
      start[newi+1] = start[newi] + nneighbors;
      newi++;
      cnedges += nneighbors;

      /* Now clear the seenflags */
      jptr = graph[i]->edges;
      for (j=graph[i]->nedges-1; j; j--) seenflag[v2cv[*(++jptr)]] = 0;
      if (mflag[i] > i) {
        jptr = graph[mflag[i]]->edges;
        for (j=graph[mflag[i]]->nedges-1; j; j--) seenflag[v2cv[*(++jptr)]] = 0;
      }
    }
  }
  *pcnedges = cnedges/2;
}
