// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


void countup(graph, nvtxs, assignment, ndims, nsets_glob, print_lev, outfile)
struct vtx_data **graph;	/* graph data structure */
int nvtxs;			/* number of vtxs in graph */
short *assignment;		/* set number of each vtx (length nvtxs+1) */
int ndims;			/* number of cuts at each level */
int nsets_glob;			/* total number of sets for division of graph */
int print_lev;			/* level of output */
FILE *outfile;			/* output file if not NULL */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   float *hopsize;		/* number of hops for each set */
   float *cutsize;		/* number of cuts for each set */
   int *setsize;		/* number of vtxs in each set */
   int *setseen;		/* flags for sets adjacent to a particular set */
   int *inorder;		/* list of vtxs in each set */
   int *startptr;		/* indices into inorder array */
   float ncuts;			/* total number of edges connecting sets */
   float nhops;			/* total cuts weighted by hypercube hops */
   float ewgt;			/* edge weight */
   int nsets;			/* number of sets after a level */
   int vtx;			/* vertex in graph */
   int set, set2, set3;		/* sets neighboring vtxs are assigned to */
   int onbdy;			/* counts number of neighboring set for a vtx */
   int bdyvtxs;			/* sum of onbdy values for a set */
   int neighbor_sets;		/* number of neighboring sets for a set */
   int total_bdyvtxs;		/* sum of all onbdy values in whole graph  */
   int total_neighbors;		/* number of neighboring sets in graph */
   int neighbor;		/* neighbor of a vertex */
   int bits;			/* bit pattern for counting hops */
   int start_dims;		/* starting dimension for output loop */
   int using_ewgts;		/* are edge weights being used? */
   int level;			/* recursion level of partition */
   int print2file;		/* should I print to a file? */
   int i, j, k, l;		/* loop counters */
   

   print2file = (outfile != NULL);
   using_ewgts = (graph[1]->ewgts != NULL);
   ewgt = 1;

   cutsize = (float *) (MEM_OK = smalloc((unsigned) nsets_glob*sizeof(float));
   if (!MEM_OK) return;
   hopsize = (float *) (MEM_OK = smalloc((unsigned) nsets_glob*sizeof(float));
   if (!MEM_OK) return;
   setsize = (int *) (MEM_OK = smalloc((unsigned) nsets_glob*sizeof(int));
   if (!MEM_OK) return;
   setseen = (int *) (MEM_OK = smalloc((unsigned) nsets_glob*sizeof(int));
   if (!MEM_OK) return;
   startptr = (int *) (MEM_OK = smalloc((unsigned) (nsets_glob+1)*sizeof(int));
   if (!MEM_OK) return;
   inorder = (int *) (MEM_OK = smalloc((unsigned) nvtxs*sizeof(int));
   if (!MEM_OK) return;
   for (j=0; j<nsets_glob; j++) setsize[j] = 0;
   for (i=1; i<=nvtxs; i++) ++setsize[assignment[i]];

   /* Modify setsize to become index into vertex list. */
   for (j=1; j<nsets_glob; j++) setsize[j] += setsize[j-1];
   for (j=nsets_glob-1; j>0; j--) startptr[j] = setsize[j] = setsize[j-1];
   startptr[0] = setsize[0] = 0;
   startptr[nsets_glob] = nvtxs;
   for (i=1; i<=nvtxs; i++) {
      set = assignment[i];
      inorder[setsize[set]] = i;
      setsize[set]++;
   }

     for (j=0; j<nsets_glob; j++)
     {  
	    cutsize[j] = 0;
	    hopsize[j] = 0;
	    setsize[j] = 0;
	 }

	 for (i=1; i<=nvtxs; i++) {
	    set = assignment[i];
	    setsize[set] += graph[i]->vwgt;
	    for (j=1; j<graph[i]->nedges; j++) {
	       neighbor = graph[i]->edges[j];
	       set2 = assignment[neighbor];
	       if (set != set2) {
		  if (using_ewgts) ewgt = graph[i]->ewgts[j];
		  cutsize[set] += ewgt;
		  bits = set^set2;
		  for (l=bits; l; l>>=1) if (l&1) hopsize[set] += ewgt;
	       }
	    }
	 }
	 ncuts = nhops = 0;
	 total_bdyvtxs = total_neighbors = 0;

	 for (set=0; set<nsets_glob; set++) {
	    /* Compute number of set neighbors, and number of vtxs on boundary. */
	    bdyvtxs = 0;
	    for (i=0; i<nsets_glob; i++) setseen[i] = 0;
	       set2 = set;
	       onbdy = 0;
	       for (i=startptr[set2]; i<startptr[set2+1]; i++) {
		  vtx = inorder[i];
	          for (j=1; j<graph[vtx]->nedges; j++) {
	             neighbor = graph[vtx]->edges[j];
	             set3 = assignment[neighbor];
		     if (set3 != set) {		/* Is vtx on boundary? */
			/* Has this neighboring set been seen already? */
			if (setseen[set3] >= 0) {
			   ++onbdy;
			   setseen[set3] = -setseen[set3] - 1;
			}
		     }
		  }
		  /* Now reset all the setseen values to be positive. */
		  if (onbdy != 0) for (j=1; j<graph[vtx]->nedges; j++) {
	             neighbor = graph[vtx]->edges[j];
	             set3 = assignment[neighbor];
		     if (setseen[set3] < 0) setseen[set3] = -setseen[set3];
		  }
       }
	       bdyvtxs += onbdy;
	    /* Now count up the number of neighboring sets. */
	    neighbor_sets = 0;
	    for (i=0; i<nsets_glob; i++) {
	       if (setseen[i] != 0) ++neighbor_sets;
	    }

	    ncuts += cutsize[set];
	    nhops += hopsize[set];
	    total_bdyvtxs += bdyvtxs;
	    total_neighbors += neighbor_sets;
	 }
         ncuts /= 2;
         nhops /= 2;
        {char buf[150]; sprintf(buf,"Total distinct  cuts: %g   hops: %g   bndy_vtxs: %d   adj_sets: %d\n", ncuts, nhops, total_bdyvtxs, total_neighbors);UserWrite(buf);}
         if (print2file) fprintf(outfile, "Total distinct  cuts: %g   hops: %g   bndy_vtxs: %d   adj_sets: %d\n", 
            ncuts, nhops, total_bdyvtxs, total_neighbors);

   sfree((char *) cutsize);
   sfree((char *) hopsize);
   sfree((char *) setsize);
   sfree((char *) setseen);
   sfree((char *) startptr);
   sfree((char *) inorder);
}
