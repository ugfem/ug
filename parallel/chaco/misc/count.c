// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/params.h"
#include        "../main/structs.h"
#include        "../main/defs.h"


void count(graph, nvtxs, sets, nsets, hops, dump)
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vtxs in graph */
short *sets;                    /* processor each vertex is assigned to */
int nsets;                      /* number of sets partitioned into */
short (*hops)[MAXSETS];         /* hops metric between sets */
int dump;                       /* flag for extended output */
{
  int nguys[MAXSETS_TOT];       /* number of vtxs in each set */
  int using_ewgts;              /* are edge weights being used? */
  int ncross;                   /* number of outgoing edges */
  int nhops;                    /* number of hops */
  int neighbor;                 /* neighbor of a vertex */
  int nmax, nmin;               /* largest and smallest set sizes */
  int i, j;                     /* loop counters */

  for (i=0; i<MAXSETS_TOT; i++) nguys[i] = 0;

  using_ewgts = (graph[1]->ewgts != NULL);

  ncross = nhops = 0;
  for (i=1; i<=nvtxs; i++) {
    nguys[sets[i]] += graph[i]->vwgt;

    for (j=1; j<graph[i]->nedges; j++) {
      neighbor = graph[i]->edges[j];
      if (sets[neighbor] != sets[i]) {
        if (using_ewgts) {
          ncross += graph[i]->ewgts[j];
          nhops += graph[i]->ewgts[j]*hops[sets[i]][sets[neighbor]];
        }
        else {
          ncross++;
          nhops += hops[sets[i]][sets[neighbor]];
        }
      }
    }
  }

  ncross /= 2;
  nhops /= 2;

  nmax = 0;
  nmin = nvtxs;
  for (i=0; i<nsets; i++) {
    if (nguys[i] > nmax) nmax = nguys[i];
    if (nguys[i] < nmin) nmin = nguys[i];
  }
  {char buf[150]; sprintf(buf,"In subgraph: Cuts=%d, Hops=%d; Max=%d, Min=%d (nvtxs=%d).\n",
                          ncross, nhops, nmax, nmin, nvtxs);UserWrite(buf);}

  if (dump) {
    for (i=0; i<nsets; i++) {char buf[150]; sprintf(buf," Size of %d = %d\n", i, nguys[i]);UserWrite(buf);}

    for (i=0; i<nvtxs; i++) {char buf[150]; sprintf(buf,"%d\n", sets[i]);UserWrite(buf);}
    {char buf[150]; sprintf(buf,"\n\n");UserWrite(buf);}
  }
}
