// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


void makecgraph_kl(graph, nvtxs, pcgraph, pcnvtxs, pcnedges, mflag,
	      v2cv, nmerged, using_ewgts, igeom, coords, ccoords)
struct vtx_data **graph;	/* array of vtx data for graph */
int nvtxs;			/* number of vertices in graph */
struct vtx_data ***pcgraph;	/* coarsened version of graph */
int *pcnvtxs;			/* number of vtxs in coarsened graph */
int *pcnedges;			/* number of edges in coarsened graph */
int *mflag;			/* flag indicating vtx matched or not */
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
   extern int DEBUG_COARSEN;	/* debug flag for coarsening output */
   struct vtx_data **cgraph;	/* coarsened version of graph */
   struct vtx_data *links;	/* space for all the vertex data */
   struct vtx_data **gptr;	/* loops through cgraph */
   struct vtx_data *cgptr;	/* loops through cgraph */
   int *start;			/* start of edgevals list for each vertex */
   int *iptr;			/* loops through integer arrays */
   short *seenflag;		/* flags for vtxs already put in edge list */
   short *sptr;			/* loops through seenflags */
   float *eweights;		/* space for edge weights in coarsened graph */
   float *fptr;			/* loops through eweights */
   float ewgt;			/* edge weight */
   double ewgt_sum;		/* sum of edge weights */
   double time;			/* timing parameters */
   int nseen;			/* number of edges of coarse graph seen so far */
   int cnvtxs;			/* number of vtxs in coarsened graph */
   int cnedges;			/* twice number of edges in coarsened graph */
   int neighbor;		/* neighboring vertex */
   int size;			/* space needed for coarsened graph */
   int *edges;			/* space for edges in coarsened graph */
   int newvtx;			/* vertex number in coarsened graph */
   int newneighbor;		/* neighboring vertex number in coarsened graph */
   double m1, m2;		/* vertex weights of vertices being merged */
   int v1, v2;			/* vertices being merged */
   int i, j;			/* loop counters */
   double seconds();
   
   void makev2cv(), countcedges();

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

   size = 2*cnedges + cnvtxs;
   edges = (int *) (MEM_OK = smalloc((unsigned) size*sizeof(int));
   if (!MEM_OK) return;
   eweights = (float *) (MEM_OK = smalloc((unsigned) size*sizeof(float));
   if (!MEM_OK) return;

   /* Fill in simple data fields for coarsened graph. */
   /* Edges and weights are put in later. */
   gptr = cgraph;
   for (i=1; i<=cnvtxs; i++) {
      links->vwgt = 0;
      size = start[i+1] - start[i] + 1;
      links->nedges = size;
      links->edges = edges;
      links->ewgts = eweights;
      links->edges[0] = i;
      edges += size;
      eweights += size;
      *(++gptr) = links++;
   }
   sfree((char *) start);

   /* Now form new vertex weights by adding those from contracted edges. */
   gptr = graph;
   for (i=1; i<=nvtxs; i++) {
      cgraph[v2cv[i]]->vwgt += (*(++gptr))->vwgt;
   }

   /* Use the renumbering to fill in the edge lists for the new graph. */
   ewgt = 1;
   for (i=1; i<=nvtxs; i++) {
      nseen = 1;
      newvtx = v2cv[i];
      cgptr = cgraph[newvtx];
      ewgt_sum = 0;

      /* Unmatched edge, or first appearance of matched edge. */
      if (mflag[i] > i || mflag[i] == 0) {
         iptr = graph[i]->edges;
         if (using_ewgts) fptr = graph[i]->ewgts;
         for (j=graph[i]->nedges-1; j; j--) {
	    neighbor = *(++iptr);
	    if (v2cv[neighbor] != newvtx) {
	       newneighbor = v2cv[neighbor];
	       if (using_ewgts) ewgt = *(++fptr);
	       ewgt_sum += ewgt;

	       if (seenflag[newneighbor] == 0) {	/* New neighbor. */
	          cgptr->edges[nseen] = newneighbor;
	          cgptr->ewgts[nseen] = ewgt;
		  seenflag[newneighbor] = nseen++;
	       }
	       else {				/* Already seen neighbor. */
	          cgptr->ewgts[seenflag[newneighbor]] += ewgt;
	       }

	    }
	    else if (using_ewgts) ++fptr;
         }
      }

      /* Now handle the matched vertex. */
      if (mflag[i] > i) {
         iptr = graph[mflag[i]]->edges;
         if (using_ewgts) fptr = graph[mflag[i]]->ewgts;
         for (j=graph[mflag[i]]->nedges-1; j; j--) {
	    neighbor = *(++iptr);
	    if (v2cv[neighbor] != newvtx) {
	       newneighbor = v2cv[neighbor];
	       if (using_ewgts) ewgt = *(++fptr);
	       ewgt_sum += ewgt;

	       if (seenflag[newneighbor] == 0) {	/* New neighbor. */
	          cgptr->edges[nseen] = newneighbor;
	          cgptr->ewgts[nseen] = ewgt;
		  seenflag[newneighbor] = nseen++;
	       }
	       else {				/* Already seen neighbor. */
	          cgptr->ewgts[seenflag[newneighbor]] += ewgt;
	       }
	    }
	    else if (using_ewgts) ++fptr;
         }
      }
      if (mflag[i] > i || mflag[i] == 0) {
	 cgptr->ewgts[0] = -ewgt_sum;
	 /* Now clear the seenflag values. */
         iptr = cgraph[newvtx]->edges;
	 for (j=cgraph[newvtx]->nedges-1; j; j--) {
	    seenflag[*(++iptr)] = 0;
	 }
      }
   }


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
