// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


void inertial1d(graph, nvtxs, ndims, x, sets, goal, using_vwgts, part_type)
struct vtx_data **graph;	/* graph data structure */
int nvtxs;			/* number of vtxs in graph */
int ndims;			/* number of cuts to make at once */
float *x;			/* x coordinates of vertices */
short *sets;			/* set each vertex gets assigned to */
double *goal;			/* desired set sizes */
int using_vwgts;		/* are vertex weights being used? */
int part_type;
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern double median_time;	/* time to find medians */
   double *value;		/* values passed to median routine */
   double time;			/* timing variables */
   int *space;			/* space required by median routine */
   int i;			/* loop counter */
   void rec_median_1();
   double seconds();
   

   value = (double *) (MEM_OK = smalloc((unsigned)(nvtxs+1)*sizeof(double));
   if (!MEM_OK) return;

   /* Copy values into double precision array. */
   for (i=1; i<=nvtxs; i++) value[i] = x[i];

   /* Now find the median value and partition based upon it. */
   space = (int *) (MEM_OK = smalloc((unsigned)nvtxs*sizeof(int));
   if (!MEM_OK) return;

   time = seconds();
   rec_median_1(graph, value, nvtxs, space, ndims, goal, using_vwgts, sets, 
	            TRUE, part_type);
   if (!MEM_OK) return;
   median_time += seconds() - time;

   sfree((char *) space);
   sfree((char *) value);
}
