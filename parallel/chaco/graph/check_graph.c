// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <math.h>
#include "../main/defs.h"
#include "../main/structs.h"
#include "general.h"

/* Check graph for errors */

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

static int is_an_edge();

int check_graph(graph, nvtxs, nedges)
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vertices */
int nedges;                     /* number of edges */
{
  float eweight;                /* edge weight */
  float wgt_sum;                /* sum of edge weights */
  int flag;                     /* flag for error free graph */
  int no_edge_count;            /* warning flag for isolated vertex */
  int using_ewgts;              /* are edge weights being used? */
  int narcs;                    /* number of neighbors of a vertex */
  int neighbor;                 /* neighbor of a vertex */
  int i, j;                     /* loop counters */

  flag = TRUE;
  no_edge_count = 0;
  using_ewgts = (graph[1]->ewgts != NULL);
  narcs = 0;
  for (i=1; i<=nvtxs; i++) {
    narcs += graph[i]->nedges - 1;

    if (graph[i]->edges[0] != i) {
      {char buf[150]; sprintf(buf," Self edge wrong for vtx %d\n", i);UserWrite(buf);}
      flag = FALSE;
    }

    if (graph[i]->nedges == 1) {
      if (1) {char buf[150]; sprintf(buf,"WARNING: Vertex %d has no neighbors\n", i);UserWrite(buf);}
      ++no_edge_count;
    }

    if (using_ewgts) wgt_sum = graph[i]->ewgts[0];

    for (j=1; j<graph[i]->nedges; j++) {
      neighbor = graph[i]->edges[j];
      if (using_ewgts) wgt_sum += graph[i]->ewgts[j];

      if (neighbor < 1 || neighbor > nvtxs) {
        {char buf[150]; sprintf(buf," Edge (%d,%d) included, but nvtxs = %d\n",
                                i, neighbor, nvtxs);UserWrite(buf);}
        flag = FALSE;
      }

      if (!is_an_edge(graph[neighbor], i, &eweight)) {
        {char buf[150]; sprintf(buf," Edge (%d,%d) included but not (%d,%d)\n",
                                i, neighbor, neighbor, i);UserWrite(buf);}
        flag = FALSE;
      }
      else if (using_ewgts && eweight != graph[i]->ewgts[j]) {
        {char buf[150]; sprintf(buf," Weight of (%d,%d)=%g, but weight of (%d,%d)=%g\n",
                                i, neighbor, graph[i]->ewgts[j], neighbor, i, eweight);UserWrite(buf);}
        flag = FALSE;
      }
    }

    if (using_ewgts && fabs((double) wgt_sum) > 1.0e-6) {
      {char buf[150]; sprintf(buf," Sum of edge weights for vertex %d = %g\n", i, wgt_sum);UserWrite(buf);}
      flag = FALSE;
    }
  }

  if (no_edge_count > 1) {char buf[150]; sprintf(buf,"WARNING: %d vertices have no neighbors\n", no_edge_count);UserWrite(buf);}

  if (narcs != 2*nedges) {
    {char buf[150]; sprintf(buf," twice nedges = %d, but I count %d\n", 2*nedges, narcs);UserWrite(buf);}
    flag = FALSE;
  }
  return(flag);
}


int check(graph, nvtxs, nedges)
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vertices */
int nedges;                     /* number of edges */
{
  float eweight;                /* edge weight */
  float wgt_sum;                /* sum of edge weights */
  int flag;                     /* flag for error free graph */
  int no_edge_count;            /* warning flag for isolated vertex */
  int bad_edge_count;           /* warning flag for isolated vertex */
  int using_ewgts;              /* are edge weights being used? */
  int narcs;                    /* number of neighbors of a vertex */
  int neighbor;                 /* neighbor of a vertex */
  int i, j;                     /* loop counters */
  int max,min,vmax,vmin;

  flag = TRUE;
  no_edge_count = bad_edge_count = 0;
  if (graph[1]->ewgts != NULL)
  {
    {char buf[150]; sprintf(buf,"GAU: now graph has edges!! ptraddr=%d\n", graph[1]->ewgts);UserWrite(buf);}
    using_ewgts = (graph[1]->ewgts != NULL);
  }
  else using_ewgts = 0;
  narcs = 0;
  max=vmax=vmin=0;
  min=vmin=64;

  for (i=1; i<=nvtxs; i++) {
    narcs += graph[i]->nedges - 1;

    if (graph[i]->edges[0] != i) {
      {char buf[150]; sprintf(buf," Self edge wrong for vtx %d\n", i);UserWrite(buf);}
      flag = FALSE;
    }

    if (graph[i]->nedges < 1) {
      if (1) {char buf[150]; sprintf(buf,"FATAL: Vertex %d has %d neighbors\n", i, graph[i]->nedges);UserWrite(buf);}
      ++bad_edge_count;
    }

    if (graph[i]->nedges == 1) {
      if (1) {char buf[150]; sprintf(buf,"WARNING: Vertex %d has no neighbors\n", i);UserWrite(buf);}
      ++no_edge_count;
    }

    if (graph[i]->vwgt<1 || graph[i]->vwgt>512)
    {char buf[150]; sprintf(buf,"WARNING: vertex %d has weight %d\n",i,graph[i]->vwgt);UserWrite(buf);}

    if (graph[i]->vwgt<min)
    {
      vmin = i;
      min = graph[i]->vwgt;
    }

    if (graph[i]->vwgt>max)
    {
      vmax = i;
      max = graph[i]->vwgt;
    }

    if (using_ewgts) wgt_sum = graph[i]->ewgts[0];

    for (j=1; j<graph[i]->nedges; j++) {
      neighbor = graph[i]->edges[j];
      if (using_ewgts) wgt_sum += graph[i]->ewgts[j];

      if (neighbor < 1 || neighbor > nvtxs) {
        {char buf[150]; sprintf(buf,"FATAL: Edge (%d,%d) included, but nvtxs = %d\n",
                                i, neighbor, nvtxs);UserWrite(buf);}
        flag = FALSE;
      }

      if (!is_an_edge(graph[neighbor], i, &eweight)) {
        {char buf[150]; sprintf(buf,"FATAL: Edge (%d,%d) included but not (%d,%d)\n",
                                i, neighbor, neighbor, i);UserWrite(buf);}
        flag = FALSE;
      }
      else if (using_ewgts && eweight != graph[i]->ewgts[j]) {
        {char buf[150]; sprintf(buf," Weight of (%d,%d)=%g, but weight of (%d,%d)=%g\n",
                                i, neighbor, graph[i]->ewgts[j], neighbor, i, eweight);UserWrite(buf);}
        flag = FALSE;
      }
    }

    if (using_ewgts && fabs((double) wgt_sum) > 1.0e-6) {
      {char buf[150]; sprintf(buf," Sum of edge weights for vertex %d = %g\n", i, wgt_sum);UserWrite(buf);}
      flag = FALSE;
    }
  }

  if (bad_edge_count > 0) {char buf[150]; sprintf(buf,"FATAL: %d vertices have no neighbors\n", bad_edge_count);UserWrite(buf);}
  if (no_edge_count > 0) {char buf[150]; sprintf(buf,"WARNING: %d vertices have no neighbors\n", no_edge_count);UserWrite(buf);}


  if (narcs != 2*nedges) {
    if (nedges>=0)
    {
      {char buf[150]; sprintf(buf,"FATAL: twice nedges = %d, but I count %d\n", 2*nedges, narcs);UserWrite(buf);}
      flag = FALSE;
    }
  }
  {char buf[150]; sprintf(buf,"checkgraph: vmin=%d min=%d vmax=%d max=%d\n",vmin,min,vmax,max);UserWrite(buf);}
  return(flag);
}

static int is_an_edge(vertex, v2, weight2)
struct vtx_data *vertex;        /* data for a vertex */
int v2;                         /* neighbor to look for */
float *weight2;                 /* weight of edge if found */
{
  int i;                        /* loop counter */

  for (i=1; i<vertex->nedges; i++) {
    if (vertex->edges[i] == v2) {
      if (vertex->ewgts != NULL) *weight2 = vertex->ewgts[i];
      else *weight2 = 1;
      return(TRUE);
    }
  }

  return(FALSE);
}
