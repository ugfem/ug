// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include "../main/defs.h"
#include "../main/structs.h"

/* Print out subgraph in program-readable form. */
void graphout(graph, nvtxs, tag, file_name)
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vtxs in graph */
char *tag;                      /* message to include */
char *file_name;                /* output file name if not null */
{
  FILE *file;                   /* output file */
  int using_vwgts;              /* Are vertices weighted? */
  int using_ewgts;              /* Are edges weighted? */
  int nedges;                   /* number of edges in graph */
  int option;                   /* output option */
  int i, j;                     /* loop counter */


  /* Determine all the appropriate parameters. */
  using_ewgts = (graph[1]->ewgts != NULL);
  using_vwgts = FALSE;
  nedges = 0;
  for (i=1; i<=nvtxs; i++) {
    if (graph[i]->vwgt != 1) using_vwgts = TRUE;
    nedges += graph[i]->nedges - 1;
  }

  option = 0;
  if (using_ewgts) option += 1;
  if (using_vwgts) option += 10;

  if (file_name == NULL) {
    if (tag != NULL) {char buf[150]; sprintf(buf,"%% graphout: %s\n", tag);UserWrite(buf);}
    {char buf[150]; sprintf(buf," %d %d", nvtxs, nedges/2, option);UserWrite(buf);}
    if (option != 0) {char buf[150]; sprintf(buf,"  %d", option);UserWrite(buf);}
    {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
    for (i=1; i<=nvtxs; i++) {
      if (using_vwgts) {char buf[150]; sprintf(buf,"%d ", graph[i]->vwgt);UserWrite(buf);}
      for (j=1; j<graph[i]->nedges; j++) {
        {char buf[150]; sprintf(buf," %d", graph[i]->edges[j]);UserWrite(buf);}
        if (using_ewgts) {char buf[150]; sprintf(buf," %g ", graph[i]->ewgts[j]);UserWrite(buf);}
      }
      {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
    }
  }

  else {
    file = fopen(file_name, "w");
    if (tag != NULL) fprintf(file, "%% graphout: %s\n", tag);
    fprintf(file, " %d %d", nvtxs, nedges/2, option);
    if (option != 0) fprintf(file, "  %d", option);
    fprintf(file, "\n");
    for (i=1; i<=nvtxs; i++) {
      if (using_vwgts) fprintf(file, "%d ", graph[i]->vwgt);
      for (j=1; j<graph[i]->nedges; j++) {
        fprintf(file, " %d", graph[i]->edges[j]);
        if (using_ewgts) fprintf(file, " %g ", graph[i]->ewgts[j]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }
  return;
}
