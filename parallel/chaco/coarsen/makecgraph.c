// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


void makecgraph(graph, nvtxs, pcgraph, pcnvtxs, pcnedges, mflag, merged, reduction,
	   v2cv, nmerged, using_ewgts, igeom, coords, ccoords)
struct vtx_data **graph;	/* array of vtx data for graph */
int nvtxs;			/* number of vertices in graph */
struct vtx_data ***pcgraph;	/* coarsened version of graph */
int *pcnvtxs;			/* number of vtxs in coarsened graph */
int *pcnedges;			/* number of edges in coarsened graph */
int *mflag;			/* flag indicating vtx matched or not */
struct ipairs *merged;		/* vtx pairs that get merged */
struct fpairs *reduction;	/* multipliers reducing edge weights */
int *v2cv;			/* mapping from vtxs to coarsened vtxs */
int nmerged;			/* number of merged vertices */
int using_ewgts;		/* are edge weights being used? */
int igeom;			/* dimensions of geometric data */
float **coords;			/* coordinates for vertices */
float **ccoords;		/* coordinates for coarsened vertices */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern double make_cgraph_time;
   extern double make_cewgts_time;
   extern int DEBUG_COARSEN;	/* debug flag for coarsening output */
   struct vtx_data **cgraph;	/* coarsened version of graph */
   struct vtx_data *links;	/* space for all the vertex data */
   short *seenflag;		/* flags vtxs already in edge list */
   short *sptr;			/* loops through seenflag */
   double time, time2;		/* timing parameters */
   float *eweights;		/* space for edge weights in coarsened graph */
   float *eptr;			/* loops through edge weights */
   float ewgt;			/* edge weight */
   int *start;			/* start of edgevals list for each vertex */
   int cnvtxs;			/* number of vtxs in coarsened graph */
   int cnedges;			/* twice number of edges in coarsened graph */
   int neighbor;		/* neighboring vertex */
   int cneighbor;		/* numbering in coarse graph of neighbor vertex */
   int nseen;			/* number of neighbors of a coarse vertex */
   int size;			/* space needed for coarsened graph */
   int *edges;			/* space for edges in coarsened graph */
   int *iptr;			/* loops through edge lists */
   int cvtx;			/* vertex number in coarsened graph */
   double m1, m2;		/* vertex weights of vertices being merged */
   int v1, v2;			/* vertices being merged */
   int i, j;			/* loop counters */
   
   double seconds();
   void makev2cv(), countcedges(), makecewgts();

   /* Compute the number of vertices and edges in the coarsened graph, */
   /* and construct start pointers into coarsened edge array. */
   time = seconds();

   *pcnvtxs = cnvtxs = nvtxs - nmerged;

   /* Construct mapping from original graph vtxs to coarsened graph vtxs. */
   makev2cv(mflag, nvtxs, v2cv);

   start = (int *) (MEM_OK = smalloc((unsigned) (cnvtxs+2)*sizeof(int));
   if (!MEM_OK) return;
   seenflag = (short *) (MEM_OK = smalloc((unsigned) (cnvtxs+1)*sizeof(short));
   if (!MEM_OK) return;
   sptr = seenflag;
   for (i=cnvtxs; i; i--) *(++sptr) = 0;
   countcedges(graph, nvtxs, start, seenflag, mflag, v2cv, pcnedges);
   cnedges = *pcnedges;

   if (DEBUG_COARSEN > 0) {
     {char buf[150]; sprintf(buf," Coarse graph has %d vertices and %d edges\n", cnvtxs, cnedges);UserWrite(buf);}
   }

   /* Now allocate space for the new graph. */
   *pcgraph = cgraph = (struct vtx_data **) (MEM_OK = smalloc((unsigned) (cnvtxs+1)*sizeof(struct vtx_data *));
   if (!MEM_OK) return;
   links = (struct vtx_data *) (MEM_OK = smalloc((unsigned) cnvtxs*sizeof(struct vtx_data));
   if (!MEM_OK) return;
   for (i=1; i<=cnvtxs; i++) cgraph[i] = links++;

   size = 2*cnedges + cnvtxs;
   edges = (int *) (MEM_OK = smalloc((unsigned) size*sizeof(int));
   if (!MEM_OK) return;
   eweights = (float *) (MEM_OK = smalloc((unsigned) size*sizeof(float));
   if (!MEM_OK) return;

   /* Fill in simple data fields for coarsened graph. */
   /* Edges and weights are put in later. */
   for (i=1; i<=cnvtxs; i++) {
      size = start[i+1] - start[i] + 1;
      cgraph[i]->nedges = size;
      cgraph[i]->edges = edges;
      cgraph[i]->ewgts = eweights;
      cgraph[i]->edges[0] = i;
      edges += size;
      eweights += size;
   }

   /* Use the renumbering to fill in the edge fields for the new graph. */
   /* Use the start array as an index into the edges fields. */
   /* Initialize the edge weights as the sum of those of their constituent */
   /* edges.  However, only set weights for half the edges and fix their */
   /* counterparts after adjustment for the spring analogy. */
   ewgt = 1;
   for (i=1; i<=nvtxs; i++) {
      /* Deal with unmatched vertex, or first occurence of matched edge. */
      if (mflag[i] > i || mflag[i] == 0) {
         cvtx = v2cv[i];
	 cgraph[cvtx]->vwgt = graph[i]->vwgt;
         nseen = 1;
	 iptr = graph[i]->edges;
	 if (using_ewgts) eptr = graph[i]->ewgts;
         for (j=graph[i]->nedges-1; j; j--) {
	    neighbor = *(++iptr);
	    cneighbor = v2cv[neighbor];
	    if (using_ewgts) ewgt = *(++eptr);
	    if (cvtx > cneighbor) {
	       /* Only set edge weights for edges from big to small vtxs. */
	       if (!seenflag[cneighbor]) {	/* New neighbor. */
	          cgraph[cvtx]->edges[nseen] = cneighbor;
		  cgraph[cvtx]->ewgts[nseen] = ewgt;
	          seenflag[cneighbor] = nseen++;
	       }
	       else {
		  cgraph[cvtx]->ewgts[seenflag[cneighbor]] += ewgt;
	       }
	    }
	    else if (cvtx < cneighbor) {
	       if (!seenflag[cneighbor]) {	/* New neighbor. */
	          cgraph[cvtx]->edges[nseen] = cneighbor;
	          seenflag[cneighbor] = nseen++;
	       }
	    }
	 }
      }

      /* Now handle matching vertex. */
      if (mflag[i] > i) {
	 v1 = mflag[i];
	 cgraph[cvtx]->vwgt += graph[v1]->vwgt;
	 iptr = graph[v1]->edges;
	 if (using_ewgts) eptr = graph[v1]->ewgts;
         for (j=graph[v1]->nedges-1; j; j--) {
	    neighbor = *(++iptr);
	    cneighbor = v2cv[neighbor];
	    if (using_ewgts) ewgt = *(++eptr);
	    if (cvtx > cneighbor) {
	       if (!seenflag[cneighbor]) {	/* New neighbor. */
	          cgraph[cvtx]->edges[nseen] = cneighbor;
		  cgraph[cvtx]->ewgts[nseen] = ewgt;
	          seenflag[cneighbor] = nseen++;
	       }
	       else {
		  cgraph[cvtx]->ewgts[seenflag[cneighbor]] += ewgt;
	       }
	    }
	    else if (cvtx < cneighbor) {
	       if (seenflag[cneighbor] == 0) {	/* New neighbor. */
	          cgraph[cvtx]->edges[nseen] = cneighbor;
	          seenflag[cneighbor] = nseen++;
	       }
	    }
	 }
      }
      if (mflag[i] > i || mflag[i] == 0) {
	 start[cvtx] = nseen;
         /* Now clear the seenflag values */
	 iptr = cgraph[cvtx]->edges;
         for (j=cgraph[cvtx]->nedges-1; j; j--) {
	    seenflag[*(++iptr)] = 0;
	 }
      }
   }

   sfree((char *) start);

   /* Finally, generate the edge weights for the coarsened graph */
   time2 = seconds();
   makecewgts(graph, cgraph, cnvtxs, merged, reduction, nmerged,
	      v2cv, seenflag, using_ewgts);
   if (!MEM_OK) return;
   make_cewgts_time += seconds() - time2;

   /* If desired, make new vtx coordinates = center-of-mass of their parents. */
   if (coords != NULL && ccoords != NULL && igeom > 0) {
      for (i=0; i<igeom; i++) {
	 ccoords[i] = (float *) (MEM_OK = smalloc((unsigned) (cnvtxs+1)*sizeof(float));
     if (!MEM_OK) return;
      }
      for (j=1; j<=nvtxs; j++) {
	 if (mflag[j] == 0) {	/* If unmatched, leave it alone. */
            for (i=0; i<igeom; i++) {
	       ccoords[i][v2cv[j]] = coords[i][j];
            }
	 }
	 else if (j < mflag[j]) {	/* If matched, use center of mass. */
	    v1 = j;
	    v2 = mflag[j];
	    m1 = graph[v1]->vwgt;
	    m2 = graph[v2]->vwgt;
            for (i=0; i<igeom; i++) {
	       ccoords[i][v2cv[j]] = (m1*coords[i][v1]+m2*coords[i][v2])/(m1+m2);
            }
	 }
      }
   }

   sfree((char *) seenflag);

   make_cgraph_time += seconds() - time;
}
