// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

static void reduce();

void adjcewgts(graph, cgraph, cnvtxs, merged, reduction, nmerged, v2cv, seenflag,
	       using_ewgts)
/* Adjust the weights of edges adjacent to merged edges. */
struct vtx_data **graph;	/* array of vtx data for graph */
struct vtx_data **cgraph;	/* coarsened version of graph */
int cnvtxs;			/* number of vertices in coarse graph */
struct ipairs *merged;		/* vtx pairs that get merged */
struct fpairs *reduction;	/* multipliers reducing each edge weight */
int nmerged;			/* number of vertices that were merged */
int *v2cv;			/* mapping from original to coarse vtxs */
short *seenflag;		/* space for marking cvtxs (initially zero) */
int using_ewgts;		/* are edge weights being used in graph? */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   int *flag;			/* space for cneighbors+1 values */
   int *fptr;			/* loops through flag */
   int *vptr, *vptr2;		/* loops through edge list */
   double r12, r34;		/* fractional reduction in edge weights */
   double k12, k23, k34;	/* (summed) edge wgts for nonparallel edges */
   double ewgt1, ewgt2;		/* edge weights */
   double ratio;		/* sum of inverse edge weights */
   int v2, v3;			/* vertex numbers */
   int cvmin, cvmax;		/* min/max of two vtxs defining an edge */
   int m2, m3;			/* vertex weights of contracting vertices */
   int cv, cv0;			/* contracted vertex number */
   int neighbor;		/* neighboring vertex */
   int cmaxdeg;			/* maximum number of neighbors in coarse graph */
   int i, j, k, l;		/* loop counters */
   

   /* Find maximum degree in coarse graph. */
   cmaxdeg = 1;
   for (i=1; i<=cnvtxs; i++) {
      if (cgraph[i]->nedges > cmaxdeg) cmaxdeg = cgraph[i]->nedges;
   }
   flag = (int *) (MEM_OK = smalloc((unsigned) cmaxdeg*sizeof(int));
   if (!MEM_OK) return;

   /* Loop through all merged edges in randomly permuted order. */
   ewgt1 = ewgt2 = 1;
   ratio = .5;
   k23 = 1;
   for (i=0; i<nmerged; i++) {
      v2 = min(merged[i].val1, merged[i].val2);
      v3 = max(merged[i].val1, merged[i].val2);
      cv = v2cv[v2];

      /* Find spring constant and masses. */
      m2 = graph[v2]->vwgt;
      m3 = graph[v3]->vwgt;
      if (using_ewgts) k23 = 0;

      /* Treat three classes of edges differently.  Edges forming a square */
      /* with another matching edge are ignored.  Those forming a triangle  */
      /* involving current matching edge are used to modify k23.  And the */
      /* rest are used to compute forces and spring extensions. */

      /* Consider all adjacent edges in graph and determine where they get */
      /* mapped to in coarse graph.  By looking at neighbors of both merged */
      /* vertices, I can distinguish three types of edges. */
      fptr = flag;
      for (k=cgraph[cv]->nedges-1; k; k--) *(++fptr) = 0;

      vptr2 = graph[v2]->edges;
      k = 1;
      for (j=graph[v2]->nedges-1; j; j--) {
	 neighbor = *(++vptr2);
	 if (neighbor != v3) {		/* If == v3, edge disappears. */
	    cv0 = v2cv[neighbor];
	    if (!seenflag[cv0]) {
	       flag[k] = 1;
	       seenflag[cv0] = k++;
	    }
	 }
      }

      for (j=1; j<graph[v3]->nedges; j++) {
	 neighbor = graph[v3]->edges[j];
	 if (neighbor == v2) {
	    if (using_ewgts) k23 += graph[v3]->ewgts[j];
	 }
	 else {
	    cv0 = v2cv[neighbor];

	    if (!seenflag[cv0]) {
	       flag[k] = 2;
	       seenflag[cv0] = k++;
	    }
	    else if (flag[seenflag[cv0]] == 1) {
	       /* Flag edge as not from either end. */
	       flag[seenflag[cv0]] = 3;

	       /* Collapsing edge; is it triangle or matching square? */
	       vptr = graph[v2]->edges;
	       for (l=1; l<graph[v2]->nedges && *(++vptr)!=neighbor; l++) ;
	       if (l < graph[v2]->nedges) {	/* Triangle!  Modify k23. */
		  if (using_ewgts) {
		     ewgt1 = graph[v3]->ewgts[j];
		     ewgt2 = graph[v2]->ewgts[l];
		     ratio = ewgt1*ewgt2/(ewgt1+ewgt2);
		  }
		  k23 += ratio;
	       }
	    }
	 }
      }

      /* Now add the spring constants in sets 1 and 2. */
      k12 = k34 = 0;
      for (j=1; j<cgraph[cv]->nedges; j++) {
	 if (flag[j] == 3) flag[j] = 0;
	 if (flag[j]) {
	    cvmin = min(cv, cgraph[cv]->edges[j]);
	    cvmax = max(cv, cgraph[cv]->edges[j]);
	    if (cvmax == cv) k = j;
	    else {
	       vptr = cgraph[cvmax]->edges;
	       for (k=1; *(++vptr)!=cvmin; k++) ;
	    }
	    if (flag[j] == 1) {
	       flag[j] = k;
	       k12 += cgraph[cvmax]->ewgts[k];
	    }
	    else if (flag[j] == 2) {
	       flag[j] = -k;
	       k34 += cgraph[cvmax]->ewgts[k];
	    }
	 }
      }

      /* Now use spring constants and weights to compute proportional */
      /* reduction in spring values. */
      reduce(k12, k23, k34, m2, m3, &r12, &r34);
      reduction[i].val1 = r12;
      reduction[i].val2 = r34;

/*
{char buf[150]; sprintf(buf,"Matched edge (%d,%d) has constants %g, %g, %g => reductions %g, %g\n", v2, v3, k12, k23, k34, r12, r34);UserWrite(buf);}
*/

      /* Finally, use reduction to adjust spring constant values. */
      for (j=1; j<cgraph[cv]->nedges; j++) {
	 if (flag[j]) {
	    cvmax = max(cv, cgraph[cv]->edges[j]);
	    if (flag[j] > 0) cgraph[cvmax]->ewgts[flag[j]] *= r12;
	    else cgraph[cvmax]->ewgts[-flag[j]] *= r34;
	 }
	 /* Now clear the seenflag. */
	 seenflag[cgraph[cv]->edges[j]] = 0;
      }
   }
   sfree((char *) flag);
}


static void reduce(k12, k23, k34, m2, m3, r12, r34)
double k12, k23, k34;		/* summed edge weights at left and right */
int m2, m3;			/* vertex weights in middle */
double *r12, *r34;		/* returned multipliers for left and right */
{
   double z1, z2, z3;		/* z-displacement of each segment */
   double y1, y3;		/* displacement after merging center vtxs */

   /* Determine elevations that minimize  \sum k z^2. */
   /* Handle special cases first. */
   if (k12 == 0 && k34 == 0) {
      /*z1 = z3 = 0; z2 = 1;*/
      *r12 = 0;
      *r34 = 0;
   }
   else if (k12 == 0) {
      /* z1 = 0;*/
      *r12 = 0;
      z2 = k34;
      z3 = k23;
      y3 = z3 + z2*m2/(m2+m3);
      *r34 = z3/y3;
   }
   else if (k34 == 0) {
      /*z3 = 0;*/
      *r34 = 0;
      z2 = k12;
      z1 = k23;
      y1 = z1 + z2*m3/(m2+m3);
      *r12 = z1/y1;
   }
   else {	/* General case */
      z1 = k23*k34;
      z2 = k12*k34;
      z3 = k12*k23;
      /* Find center-of-mass location for merged center vertices. */
      y1 = z1 + z2*m3/(m2+m3);
      y3 = z3 + z2*m2/(m2+m3);

      /* Determine reduction in spring constant to preserve end forces. */
      *r12 = z1/y1;
      *r34 = z3/y3;
   }
}
