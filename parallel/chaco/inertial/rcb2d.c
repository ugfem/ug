// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>
#include <stdio.h>
#include "../main/structs.h"
#include "../main/defs.h"

void rcb2d(graph, nvtxs, ndims, x, y, sets, goal, using_vwgts, part_type)
struct vtx_data **graph;	/* graph data structure for weights */
int nvtxs;			/* number of vtxs in graph */
int ndims;			/* number of cuts to make at this step */
float *x, *y;			/* x and y coordinates of vertices */
short *sets;			/* set each vertex gets assigned to */
double *goal;			/* desired set sizes */
int using_vwgts;		/* are vertex weights being used? */
int part_type;
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int DEBUG_INERTIAL;	/* debug flag for inertial method */
   extern double inertial_axis_time;	/* time spent finding inertial axis */
   extern double median_time;	/* time spent computing medians */
   double tensor[2][2];		/* inertial tensor */
   double evec[2];		/* eigenvector of tensor */
   double *value;		/* values along selected direction to sort */
   double xcm, ycm;		/* center of mass in each direction */
   double xx, yy, xy;		/* elements of inertial tensor */
   double xdif, ydif;		/* deviation from center of mass */
   double eval, res;		/* eigenvalue and error in eval calculation */
   double vwgt_sum;		/* sum of all the vertex weights */
   double time;			/* timing parameters */
   int *space;			/* space required by median routine */
   int i;			/* loop counter */
   double seconds();
   int rcbdir = 0;
   
   void evals2(), eigenvec2(), rec_median_1();

   /* Compute center of mass and total mass. */

   time = seconds();

   /* Allocate space for value array. */

   value = (double *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(double));
   if (!MEM_OK) return;

   /* calculate the bisection coordinate direction */
{
	int max_x, min_x;
	int max_y, min_y;

	max_x = min_x = x[1];
	max_y = min_y = y[1];

	for (i=1; i<=nvtxs; i++) {
		if (max_x<x[i]) max_x = x[i];
		if (min_x>x[i]) min_x = x[i];

		if (max_y<y[i]) max_y = y[i];
		if (min_y>y[i]) min_y = y[i];
	}
	if ((max_x-min_x) < (max_y-min_y)) rcbdir = 1;
}

	/* Calculate value to sort/split on for each cell. */
	/* This is inner product with eigenvector. */
	for (i=1; i<=nvtxs; i++) {
		if (rcbdir == 0)
			value[i] = x[i];
		else
			value[i] = y[i];
	}


   /* Now find the median value and partition based upon it. */
   space = (int *) (MEM_OK = smalloc((unsigned) nvtxs*sizeof(int));
   if (!MEM_OK) return;
   time = seconds();
   rec_median_1(graph, value, nvtxs, space, ndims, goal, using_vwgts, sets, 
	            FALSE, part_type);
   if (!MEM_OK) return;
   median_time += seconds() - time;

   sfree((char *) space);
   sfree((char *) value);
}
