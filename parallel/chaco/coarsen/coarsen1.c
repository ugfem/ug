// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

/*
  Generate matching.
  Construct mapping from vtxs to coarsened vtxs (v2cv).
  Count edges in coarsened graph.
  Allocate space for new graph.
  Fill in adjacency fields and vertex weights.
  Compute edge weights.
    For all original edges {
      sum weights into corresponding coarsened edges.
      (Note: contracted edges don't get included.)
      adjust weights of coarsened edges to account for missing contracted edges.
    }

    The adjustment of coarsened edge weights proceeds as follows:
    (i) Randomly permute list of matching edges
        (Since weight adjustment is order dependent.)
    (iii) For all contracted edges:
	 Use mapping from old to new edges to identify nonparallel edges.
	 Adjust weights as if these were the only edges.

    This final numerical adjustment involves the following process:
    (a) Compute z values that minimize the sum of kz^2.
    (b) Assume merged vertices end up at center of mass.
    (c) Choose new k's to preserve forces on non-merged vertices.

    (The implicit assumption here is that the z values for the vtxs adjacent
     to a merged vertex are equal.  I only apply this to edges that won't
     collapse onto each other in the coarsened graph.)

*/

void coarsen1(graph, nvtxs, nedges, pcgraph, pcnvtxs, pcnedges, pmerged, preduction,
	 pmflag, pv2cv, igeom, coords, ccoords, using_ewgts)
struct vtx_data **graph;	/* array of vtx data for graph */
int nvtxs;			/* number of vertices in graph */
int nedges;			/* number of edges in graph */
struct vtx_data ***pcgraph;	/* coarsened version of graph */
int *pcnvtxs;			/* number of vtxs in coarsened graph */
int *pcnedges;			/* number of edges in coarsened graph */
struct ipairs **pmerged;	/* vtx pairs that get merged */
struct fpairs **preduction;	/* multipliers for reducing edge weights */
int **pmflag;			/* pointer to mflag */
int **pv2cv;			/* pointer to v2cv */
int igeom;			/* dimension for geometric information */
float **coords;			/* coordinates for vertices */
float **ccoords;		/* coordinates for coarsened vertices */
int using_ewgts;		/* are edge weights being used? */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern double coarsen_time;
   extern double match_time;
   struct ipairs *merged;	/* vtx pairs that get merged */
   double time;			/* time routine is entered */
   int *mflag;			/* flag indicating vtx matched or not */
   int nmerged;			/* number of edges contracted */
   int i, j;			/* loop counters */
   double seconds();
   int maxmatch();
   void makecgraph();

   time = seconds();

   /* Allocate and initialize space. */
   *pmflag = mflag = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
   if (!MEM_OK) return;

   /* Find a maximal matching in the graph. */
   /* nmerged = maxmatch2(graph, nvtxs, nedges, mflag, using_ewgts);*/
   nedges = nedges; /* to placate alint */
   nmerged = maxmatch(graph, nvtxs, mflag, using_ewgts);
   if (!MEM_OK) return;
   /*checkmflag(mflag, nvtxs);*/
   match_time += seconds() - time;

   /* Use mflag array to generate list of merged edges. */
   *pmerged = merged = (struct ipairs *) (MEM_OK = smalloc((unsigned) nmerged*sizeof(struct ipairs));
   if (!MEM_OK) return;
   *preduction = (struct fpairs *) (MEM_OK = smalloc((unsigned) nmerged*sizeof(struct fpairs));
   if (!MEM_OK) return;
   j = 0;
   for (i=1; i<=nvtxs; i++) {
      if (mflag[i] > i) {
	 merged[j].val1 = i;
	 merged[j].val2 = mflag[i];
	 j++;
      }
   }

   /* Now construct coarser graph by contracting along matching edges. */
   /* Pairs of values in mflag array indicate matched vertices. */
   /* A zero value indicates that vertex is unmatched. */

   *pv2cv = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
   if (!MEM_OK) return;

   makecgraph(graph, nvtxs, pcgraph, pcnvtxs, pcnedges, mflag, merged,
	      *preduction, *pv2cv, nmerged, using_ewgts, igeom, coords, ccoords);
   if (!MEM_OK) return;
   coarsen_time += seconds() - time;
}
