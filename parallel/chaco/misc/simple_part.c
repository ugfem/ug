// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"../main/structs.h"
#include	"../main/defs.h"


/* Partition vertices into sets in one of several simplistic ways. */
void simple_part(graph, nvtxs, sets, nsets, simple_type, goal)
struct vtx_data **graph;	/* data structure for graph */
int nvtxs;			/* total number of vtxs in graph */
short *sets;			/* sets vertices get assigned to */
int nsets;			/* number of sets at each division */
int simple_type;		/* type of decomposition */
double *goal;			/* desired set sizes */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int DEBUG_TRACE;	/* trace the execution of the code */
   double cutoff;		/* ending weight for a partition */
   int using_vwgts;		/* are vertex weights active? */
   int *order;			/* random ordering of vertices */
   int weight;			/* sum of vertex weights in a partition */
   int i, j;			/* loop counters */
   double *assigned;       /* pointer to array for number of nvtxs */
   
   void exit(), randomize();

   assigned = (double *) (MEM_OK = smalloc((unsigned) nsets*sizeof(double));
   if (!MEM_OK) return;
   for (i=0; i<nsets; i++) assigned[i] = 0;
   using_vwgts = (graph != NULL);

   /* Scattered initial decomposition. */
   if (simple_type == 1) {
      if (DEBUG_TRACE > 0) {
	{char buf[150]; sprintf(buf,"Generating scattered partition, nvtxs = %d\n", nvtxs);UserWrite(buf);}
      }
      j = 0;
      for (i=1; i<=nvtxs; ) {
         if (goal[j] > assigned[j])
         {
            sets[i] = j;
            assigned[j]++;
            i++;
         }
         j++;
         if (j == nsets) j = 0;
      }
   }

   /* Random initial decomposition. */
   if (simple_type == 2) {
      if (DEBUG_TRACE > 0) {
	{char buf[150]; sprintf(buf,"Generating random partition, nvtxs = %d\n", nvtxs);UserWrite(buf);}
      }
      /* Construct random order in which to step through graph. */
      order = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
      if (!MEM_OK) return;
      for (i=1; i<=nvtxs; i++) order[i] = i;
      randomize(order, nvtxs);

      j = 0;
      for (i=1; i<=nvtxs; ) {
         if (goal[j] > assigned[j])
         {
            sets[order[i]] = j;
            if (using_vwgts) assigned[j] += graph[order[i]]->vwgt;
	        else assigned[j]++;
            i++;
         }
         else j++;
         if (j == nsets) j = 0;
      }
      sfree((char *) order);
   }

   /* Linearly ordered initial decomposition. */
   if (simple_type == 3) {
      if (DEBUG_TRACE > 0) {
	{char buf[150]; sprintf(buf,"Generating linear partition, nvtxs = %d\n", nvtxs);UserWrite(buf);}
      }
      j = 0;
      for (i=1; i<=nvtxs; ) {
         if (goal[j] > assigned[j])
         {
            sets[i] = j;
            if (using_vwgts) assigned[j] += graph[i]->vwgt;
	        else assigned[j]++;
            i++;
         }
         else j++;
         if (j == nsets) j = 0;
      }
   }

   /* Read initialization from file or array. */
   if (simple_type == 4) {
     {char buf[150]; sprintf(buf," Sorry, (simple_type = 4) option not yet functional.\n");UserWrite(buf);}
      exit(1);
   }
}
