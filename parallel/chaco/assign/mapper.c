// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"../main/params.h"
#include	"../main/defs.h"
#include	"../main/structs.h"


void mapper(graph, xvecs, nvtxs, active, sets, ndims, mediantype, goal,
       vwgt_max, part_type)
struct vtx_data **graph;	/* data structure with vertex weights */
double **xvecs;			/* continuous indicator vectors */
int nvtxs;			/* number of vtxs in graph */
int *active;			/* space for nmyvals ints */
short *sets;			/* returned processor assignment for my vtxs */
int ndims;			/* number of dimensions being divided into */
int mediantype;			/* type of eigenvector partitioning to use */
double *goal;			/* desired set sizes */
int vwgt_max;			/* largest vertex weight */
int part_type;
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;  /* variable for memory overflow exception */
   double temp_goal[2];		/* combined set goals if using option 1. */
   double wbelow;		/* weight of vertices with negative values */
   short *temp_sets;		/* sets vertices get assigned to */
   int vweight;			/* weight of a vertex */
   int using_vwgts;		/* are vertex weights being used? */
   int bits;			/* bits for assigning set numbers */
   int i, j;			/* loop counters */
   void median(), rec_median_k(), median_assign();
   void map2d(), map3d();

   /* NOTE: THIS EXPECTS XVECS, NOT YVECS! */

   using_vwgts = (vwgt_max != 1);

   if (ndims == 1 && mediantype == 3) mediantype = 1;

   if (mediantype == 3) {	/* Divide using bipartite matching. */
      if (ndims == 2) {map2d(graph, xvecs, nvtxs, sets, goal, vwgt_max);
                       if (!MEM_OK) return;
                      }
      else if (ndims == 3) {map3d(graph, xvecs, nvtxs, sets, goal, vwgt_max);
                            if (!MEM_OK) return;
                           } 
   }

   else if (mediantype == 2) {	/* Divide recursively using medians. */
      rec_median_k(graph, xvecs, nvtxs, active, ndims, goal, using_vwgts, sets, 
	               part_type);
      if (!MEM_OK) return;
   }


   else if (mediantype == 1) {	/* Cut with independent medians => unbalanced. */
      bits = 1;
      temp_sets = (short *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(short));
      if (!MEM_OK) return;
      for (j=1; j<=nvtxs; j++) sets[j] = 0;

      for (i=1; i<=ndims; i++) {
	 temp_goal[0] = temp_goal[1] = 0;
	 for (j=0; j<(1<<ndims); j++) {
	    if (bits & j) temp_goal[0] += goal[j];
	    else temp_goal[1] += goal[j];
	 }
	 bits <<= 1;

	 median(graph, xvecs[i], nvtxs, active, temp_goal, using_vwgts, temp_sets);
	 for (j=1; j<=nvtxs; j++) sets[j] = (sets[j]<<1) + temp_sets[j];
      }
      sfree((char *) temp_sets);
   }

   else if (mediantype == 0) {	/* Divide at zero instead of median. */
      bits = 1;
      temp_sets = (short *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(short));
      if (!MEM_OK) return;
      for (j=1; j<=nvtxs; j++) sets[j] = 0;

      for (i=1; i<=ndims; i++) {
	 temp_goal[0] = temp_goal[1] = 0;
	 for (j=0; j<(1<<ndims); j++) {
	    if (bits & j) temp_goal[0] += goal[j];
	    else temp_goal[1] += goal[j];
	 }
	 bits <<= 1;

	 wbelow = 0;
	 vweight = 1;
	 for (j=1; j<=nvtxs; j++) {
	    if (using_vwgts) vweight = graph[i]->vwgt;
	    if (xvecs[i][j] < 0) wbelow += vweight;
	 }

	 median_assign(graph, xvecs[i], nvtxs, temp_goal, using_vwgts, temp_sets,
		       wbelow, (double) 0.0);

	 for (j=1; j<=nvtxs; j++) sets[j] = (sets[j]<<1) + temp_sets[j];
      }
   }
}
