// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"../main/structs.h"
#include	"../main/params.h"
#include	"../main/defs.h"


/* Recursively apply median to a SINGLE vector of values */

void rec_median_1(graph, vals, nvtxs, active, ndims, goal, using_vwgts, 
                  assign, top, part_type)
struct vtx_data **graph;	/* data structure with vertex weights */
double *vals;			/* values of which to find median */
int nvtxs;			/* number of values I own */
int *active;			/* space for list of nmyvals ints */
int ndims;			/* number of dimensions to divide */
double *goal;			/* desired sizes for sets */
int using_vwgts;		/* are vertex weights being used? */
short *assign;			/* set each vertex gets assigned to */
int top;			/* is this the top call in the recursion? */
int part_type;
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int ARCH_GOAL;
   struct vtx_data **sub_graph;	/* subgraph data structure with vertex weights */
   double *sub_vals;		/* subgraph entries in vals vector */
   double merged_goal[MAXSETS/2];	/* combined goal values */
   double sub_vwgt_sum;		/* sum of vertex weights in subgraph */
   int *loc2glob;		/* mapping from subgraph to graph numbering */
   int inv_gray[MAXSETS];	/* inverse gray code of set number */
   short *sub_assign;		/* assignment returned from subgraph */
   int sub_nvtxs;		/* number of vertices in subgraph */
   int setsize[2];		/* number of vertices in two subsets */
   int maxsize;			/* number of vertices in larger subset */
   int nsets;			/* number of sets we are dividing into */
   int i, j;			/* loop counters */
   int _gray_(), bit_reverse();
   void rec_median_1(), median(), make_subgoal(), make_subvector();
   void make_maps(), merge_assignments();

   /* Note: This is HYPERCUBE/MESH dependent.  We'll want to combine the */
   /* sizes of different sets on the different architectures. */
   nsets = 1<<ndims;
	if (ARCH_GOAL==0)
	{
		merge_goals(merged_goal, goal, ndims, part_type);
	}
	else 
	{
		for (i=0; i<2; i++) 
		{
		merged_goal[i] = 0;
		for (j=i; j<nsets; j+=2) merged_goal[i] += goal[j];
		}
	}

   median(graph, vals, nvtxs, active, merged_goal, using_vwgts, assign);

   if (ndims > 1) {
      /* Find size of largest subproblem. */
      setsize[0] = setsize[1] = 0;
      for (i=1; i<=nvtxs; i++) ++setsize[assign[i]];
      maxsize = max(setsize[0], setsize[1]);

      sub_assign = (short *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(short));
      if (!MEM_OK) return;
      sub_vals = (double *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(double));
      if (!MEM_OK) return;
      loc2glob = (int *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(int));
      if (!MEM_OK) return;
      if (!using_vwgts) sub_graph = NULL;
      else 
      {
           sub_graph = (struct vtx_data **) (MEM_OK =
			        smalloc((unsigned) (maxsize+1)*sizeof(struct vtx_data *));
           if (!MEM_OK) return;
      }
      

      for (i=0; i<2; i++) {
         /* Construct subproblem. */
	 sub_nvtxs = setsize[i];
	 make_maps(assign, nvtxs, i, (int *) NULL, loc2glob);

	 if (!using_vwgts) {
	    sub_vwgt_sum = sub_nvtxs;
	 }
	 else {
	    sub_vwgt_sum = 0;
	    for (j=1; j<=sub_nvtxs; j++) {
	       sub_graph[j] = graph[loc2glob[j]];
	       sub_vwgt_sum += sub_graph[j]->vwgt;
	    }
	 }

	 make_subvector(vals, sub_vals, sub_nvtxs, loc2glob);

	make_subgoal_med(goal, merged_goal, 1, 1<<ndims, i, 
		                   sub_vwgt_sum, part_type);

         rec_median_1(sub_graph, sub_vals, sub_nvtxs, active, ndims-1, merged_goal,
		     using_vwgts, sub_assign, FALSE,part_type);
         if (!MEM_OK) return;

	 /* Merge subassignment with current assignment. */
	 /* Note: Following is HYPERCUBE/MESH dependent. */
	 merge_assignments_med(assign, sub_assign, 1, sub_nvtxs, loc2glob, part_type);
      }
      /* Now I have the set stiped, with set assignments in bit-reversed order. */
      /* For hypercubes, it's better to Grey code. */
      /* Hence the following is HYPERCUBE/MESH dependent. */
      /* For a mesh, I probably want bit-reverse, but not gray coded. */
      if (top) {
         for (i=0; i<nsets; i++) inv_gray[i] = _gray_(bit_reverse(i, ndims));
         for (i=1; i<=nvtxs; i++) assign[i] = inv_gray[assign[i]];
      }

      if (sub_graph != NULL) sfree((char *) sub_graph);
      sfree((char *) loc2glob);
      sfree((char *) sub_vals);
      sfree((char *) sub_assign);
   }
}



/* Recursively apply median to a SEVERAL vectors of values */
/* Divide with first, and use result to divide with second, etc. */

void rec_median_k(graph, vals, nvtxs, active, ndims, goal, using_vwgts, 
                  assign, part_type)
struct vtx_data **graph;	/* data structure with vertex weights */
double **vals;			/* values of which to find median */
int nvtxs;			/* number of values I own */
int *active;			/* space for list of nmyvals ints */
int ndims;			/* number of dimensions to divide */
double *goal;			/* desired sizes for sets */
int using_vwgts;		/* are vertex weights being used? */
short *assign;			/* set each vertex gets assigned to */
int part_type;
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */ 
   extern int  ARCH_GOAL;
   struct vtx_data **sub_graph;	/* subgraph data structure with vertex weights */
   double *sub_vals[MAXDIMS];	/* subgraph entries in vals vectors */
   double merged_goal[MAXSETS/2];	/* combined goal values */
   double sub_vwgt_sum;		/* sum of vertex weights in subgraph */
   int *loc2glob;		/* mapping from subgraph to graph numbering */
   short *sub_assign;		/* assignment returned from subgraph */
   int sub_nvtxs;		/* number of vertices in subgraph */
   int setsize[2];		/* number of vertices in two subsets */
   int maxsize;			/* number of vertices in larger subset */
   int nsets;			/* number of sets we are dividing into */
   int i, j;			/* loop counters */
   short *ready;
   void rec_median_k(), median(), make_subgoal(), make_subvector();
   void make_maps(), merge_assignments();

   /* Note: This is HYPERCUBE/MESH dependent.  We'll want to combine the */
   /* sizes of different sets on the different architectures. */
   nsets = 1<<ndims;
	if (ARCH_GOAL==0)
	{
		merge_goals(merged_goal, goal, ndims, part_type);
	}
	else 
	{
		for (i=0; i<2; i++) 
		{
		merged_goal[i] = 0;
		for (j=i; j<nsets; j+=2) merged_goal[i] += goal[j];
		}
	}

   median(graph, vals[1], nvtxs, active, merged_goal, using_vwgts, assign);

   if (ndims > 1) {
      /* Find size of largest subproblem. */
      setsize[0] = setsize[1] = 0;
      for (i=1; i<=nvtxs; i++) ++setsize[assign[i]];
      maxsize = max(setsize[0], setsize[1]);

      ready = (short *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(short));
      for (i=0; i<=nvtxs; i++) ready[i] = FALSE;
      sub_assign = (short *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(short));
      if (!MEM_OK) return;
      for (i=1; i<ndims; i++) {
         sub_vals[i] = (double *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(double));
         if (!MEM_OK) return;
      }
      loc2glob = (int *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(int));
      if (!MEM_OK) return;
      if (!using_vwgts) sub_graph = NULL;
      else 
      {
          sub_graph = (struct vtx_data **) (MEM_OK =
			         smalloc((unsigned) (maxsize+1)*sizeof(struct vtx_data *));
          if (!MEM_OK) return;
      }

      for (i=0; i<2; i++) {
         /* Construct subproblem. */
	 sub_nvtxs = setsize[i];
	 make_maps(assign, nvtxs, i, (int *) NULL, loc2glob, ready);

	 if (!using_vwgts) {
	    sub_vwgt_sum = sub_nvtxs;
	 }
	 else {
	    sub_vwgt_sum = 0;
	    for (j=1; j<=sub_nvtxs; j++) {
	       sub_graph[j] = graph[loc2glob[j]];
	       sub_vwgt_sum += sub_graph[j]->vwgt;
	    }
	 }

	 for (j=2; j<=ndims; j++) {
	    make_subvector(vals[j], sub_vals[j-1], sub_nvtxs, loc2glob);
	 }

	make_subgoal_med(goal, merged_goal, 1, 1<<ndims, i, 
		                   sub_vwgt_sum, part_type);

         rec_median_k(sub_graph, sub_vals, sub_nvtxs, active, ndims-1, 
		              merged_goal, using_vwgts, sub_assign, part_type>>1);
         if (!MEM_OK) return;

	 /* Merge subassignment with current assignment. */
	 /* Note: Following is HYPERCUBE/MESH dependent. */
		merge_assignments_med(assign, sub_assign, ndims, 1, sub_nvtxs, 
		                  loc2glob, part_type, ready);
      }
      if (sub_graph != NULL) sfree((char *) sub_graph);
      sfree((char *) loc2glob);
      for (i=1; i<ndims; i++) sfree((char *) sub_vals[i]);
      sfree((char *) sub_assign);
   }
}
