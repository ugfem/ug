// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


/* Find a maximal matching in a graph using simple greedy algorithm. */
/* Randomly permute vertices, and then have each select an unmatched */
/* neighbor. */

int maxmatch(graph, nvtxs, mflag, using_ewgts)
/* Find a maximal matching in the graph, using variant of Luby's algorithm. */
/* Assign a random value to each edge, and include an edge in */
/* the matching if its value larger than all neighboring edges. */
struct vtx_data **graph;	/* array of vtx data for graph */
int nvtxs;			/* number of vertices in graph */
int *mflag;			/* flag indicating vtx selected or not */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   int *order;			/* random ordering of vertices */
   int *iptr, *jptr;		/* loops through integer arrays */
   double prob_sum;		/* sum of probabilities to select from */
   double val;			/* random value for selecting neighbor */
   float ewgt;			/* edge weight */
   int save;			/* neighbor vertex if only one active */
   int vtx;			/* vertex to process next */
   int neighbor;		/* neighbor of a vertex */
   int nmerged;			/* number of edges in matching */
   int i, j;			/* loop counters */
   
   double drandom();
   void randomize();

   /* First, randomly permute the vertices. */
   iptr = order = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
   if (!MEM_OK) return;
   jptr = mflag;
   for (i=1; i<=nvtxs; i++) {
      *(++iptr) = i;
      *(++jptr) = 0;
   }
   randomize(order, nvtxs);

   ewgt = 1;
   nmerged = 0;
   for (i=1; i<=nvtxs; i++) {
      vtx = order[i];
      if (mflag[vtx] == 0) {	/* Not already matched. */
         /* Add up sum of edge weights of neighbors. */
         prob_sum = 0;
	 save = 0;
         for (j=1; j<graph[vtx]->nedges; j++) {
	    neighbor = graph[vtx]->edges[j];
	    if (mflag[neighbor] == 0) {
	       if (using_ewgts) ewgt = graph[vtx]->ewgts[j];
	       /* Set flag for single possible neighbor. */
	       if (prob_sum == 0) save = neighbor;
	       else save = 0;
	       prob_sum += ewgt/graph[neighbor]->vwgt;
	    }
	 }

	 if (prob_sum != 0) {	/* Does vertex have contractible edges? */
	    nmerged++;
	    if (save != 0) {	/* Only one neighbor, special case. */
	       mflag[vtx] = save;
	       mflag[save] = vtx;
	    }
	    else {	/* Pick randomly neighbor, skewed by edge weights. */
	       val = drandom()*prob_sum*.999999;
	       prob_sum = 0;
               for (j=1; !mflag[vtx]; j++) {
	          neighbor = graph[vtx]->edges[j];
	          if (mflag[neighbor] == 0) {
	             if (using_ewgts) ewgt = graph[vtx]->ewgts[j];
	             prob_sum += ewgt/graph[neighbor]->vwgt;
		     if (prob_sum >= val) {
			mflag[vtx] = neighbor;
			mflag[neighbor] = vtx;
		     }
		  }
	       }
	    }
	 }
      }
   }

   sfree((char *) order);
   return(nmerged);
}
