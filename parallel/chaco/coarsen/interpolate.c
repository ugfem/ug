// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


/* Interpolate the eigenvector from the coarsened to the original graph.
   This may require regenerating mflag and v2cv arrays from merged array.

   I also need to undo the merged edges in the reverse order from that in
   which they were collapsed.
*/


void interpolate(vecs, cvecs, ndims, graph, nvtxs, cgraph, cnvtxs, merged,
	    reduction, nmerged, using_ewgts)
double **vecs;			/* approximate eigenvectors for graph */
double **cvecs;			/* exact eigenvectors for coarse graph */
int ndims;			/* number of vectors to interpolate */
struct vtx_data **graph;	/* array of vtx data for graph */
int nvtxs;			/* number of vertices in graph */
struct vtx_data **cgraph;	/* coarsened version of graph */
int cnvtxs;			/* number of vertices in coarsened graph */
struct ipairs *merged;		/* vtx pairs that got merged */
struct fpairs *reduction;	/* multipliers for adjusting edge weights */
int nmerged;			/* number of edges contracted */
int using_ewgts;		/* are edge weights being used in fine graph? */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   short *space;		/* space needed for flags in unmerge */
   int *mflag;			/* flag indicating vtx matched or not */
   int *v2cv;			/* mapping from fine to coarse vertices */
   int i, j;			/* loop counters */
   
   void makemflag(), makev2cv(), unmerge();

   /* Allocate and initialize space. */
   mflag = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
   if (!MEM_OK) return;
   v2cv = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
   if (!MEM_OK) return;

   /* Make mflag. */
   makemflag(merged, nmerged, mflag, nvtxs);

   /* Make v2cv. */
   makev2cv(mflag, nvtxs, v2cv);

   /* Now uncompress the coarse eigenvector by replicating matched values. */
   for (i=1; i<=ndims; i++) {
      for (j=1; j<=nvtxs; j++) vecs[i][j] = cvecs[i][v2cv[j]];
   }

   /* Really only need space for maximum degree in cgraph + 1. */
   space = (short *) (MEM_OK = smalloc((unsigned) cnvtxs*sizeof(short));
   if (!MEM_OK) return;

   /* Now unmerge all the edges, generating new heights for merged vertices. */
   unmerge(vecs, ndims, graph, cgraph, merged, reduction, nmerged, v2cv, space,
	   using_ewgts);

   sfree((char *) space);
   sfree((char *) v2cv);
   sfree((char *) mflag);
}
