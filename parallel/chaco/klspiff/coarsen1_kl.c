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
    }

*/



void coarsen1_kl(graph, nvtxs, nedges, pcgraph, pcnvtxs, pcnedges,
	    pv2cv, igeom, coords, ccoords, using_ewgts)
struct vtx_data **graph;	/* array of vtx data for graph */
int nvtxs;			/* number of vertices in graph */
int nedges;			/* number of edges in graph */
struct vtx_data ***pcgraph;	/* coarsened version of graph */
int *pcnvtxs;			/* number of vtxs in coarsened graph */
int *pcnedges;			/* number of edges in coarsened graph */
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
   double time;			/* time routine is entered */
   int *mflag;			/* flag indicating vtx matched or not */
   int nmerged;			/* number of edges contracted */
   double seconds();
   int maxmatch();
   void makecgraph_kl();

   time = seconds();

   /* to placate alint */
   nedges = nedges;

   /* Allocate and initialize space. */
   mflag = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
   if (!MEM_OK) return;

   /* Find a maximal matching in the graph. */
   nmerged = maxmatch(graph, nvtxs, mflag, using_ewgts);
   if (!MEM_OK) return;
   match_time += seconds() - time;

   /* Now construct coarser graph by contracting along matching edges. */
   /* Pairs of values in mflag array indicate matched vertices. */
   /* A zero value indicates that vertex is unmatched. */

   *pv2cv = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
   if (!MEM_OK) return;

   makecgraph_kl(graph, nvtxs, pcgraph, pcnvtxs, pcnedges, mflag, 
	         *pv2cv, nmerged, using_ewgts, igeom, coords, ccoords);
   if (!MEM_OK) return;

   sfree((char *) mflag);
   coarsen_time += seconds() - time;
}
