// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/structs.h"
#include        "../main/defs.h"

/* Free a graph data structure. */

void free_graph(graph)
struct vtx_data **graph;
{
  extern Heap   *heap;      /* pointer to heap of multigrid */


  if (graph[1]->ewgts != NULL) sfree((char *) graph[1]->ewgts);
  sfree((char *) graph[1]->edges);
  sfree((char *) graph[1]);
  sfree((char *) graph);
}
