// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


/* Change from a FORTRAN graph style to our graph data structure. */

void reformat(start, adjacency, nvtxs, pnedges, vwgts, ewgts, pgraph)
int *start;			/* start of edge list for each vertex */
int *adjacency;			/* edge list data */
int nvtxs;			/* number of vertices in graph */
int *pnedges;			/* ptr to number of edges in graph */
int *ewgts;			/* weights for all edges */
int *vwgts;			/* weights for all vertices */
struct vtx_data ***pgraph;	/* ptr to array of vtx data for graph */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   struct vtx_data **graph;	/* array of vtx data for graph */
   struct vtx_data *links;	/* ptr to space for data for all vtxs */
   int *edges;			/* space for all adjacency lists */
   float *eweights;		/* space for all edge weights */
   int *eptr;			/* steps through adjacency list */
   int *wptr;			/* steps through edge weights list */
   int size;			/* length of all edge lists */
   float sum;			/* sum of edge weights for a vtx */
   int using_ewgts;		/* are edge weights being used? */
   int using_vwgts;		/* are vertex weights being used? */
   int i, j;			/* loop counters */

   using_ewgts = (ewgts != NULL);
   using_vwgts = (vwgts != NULL);

   graph = (struct vtx_data **) (MEM_OK = smalloc((unsigned)(nvtxs+1)*sizeof(struct vtx_data *));
   if (!MEM_OK) return;
   *pgraph = graph;

   /* Set up all the basic data structure for the vertices. */
   /* Replace many small mallocs by a few large ones. */
   links = (struct vtx_data *) (MEM_OK = smalloc((unsigned)(nvtxs)*sizeof(struct vtx_data));
   if (!MEM_OK) return;
   for (i=1; i<=nvtxs; i++) {
      graph[i] = links++;
   }

   /* Now fill in all the data fields. */
   *pnedges = start[nvtxs+1]/2;
   size = start[nvtxs+1] + nvtxs;
   edges = (int *) (MEM_OK = smalloc((unsigned)size*sizeof(int));
   if (!MEM_OK) return;
   if (using_ewgts) eweights = (float *) (MEM_OK = smalloc((unsigned)size*sizeof(float));
   if (!MEM_OK) return;

   eptr = adjacency;
   wptr = ewgts;

   for (i=1; i<=nvtxs; i++) {
      size = start[i+1] - start[i];
      if (using_vwgts) graph[i]->vwgt = *(++vwgts);
      else graph[i]->vwgt = 1;
      graph[i]->nedges = size+1;
      graph[i]->edges = edges;
      *edges++ = i;
      for (j=size; j; j--) *edges++ = *eptr++;
      if (using_ewgts) {
         graph[i]->ewgts = eweights;
         eweights++;
         sum = 0;
         for (j=size; j; j--) {
	    sum += *wptr;
	    *eweights++ = *wptr++;
         }
         graph[i]->ewgts[0] = -sum;
      }
      else graph[i]->ewgts = NULL;
   }
}
