// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"
#include "../misc/array.h"

int nsets_glob;         /* total number of sets to partition into */

void recurse(graph, nvtxs, nedges, ndims, ndims_tot, hop_mtx, vwsqrt, spec,
	inert, rcb, KL, mediantype, mkconnected, solver_flag,
	coarse_flag, vmax, eigtol, igeom, coords, assignment, goal,
        scatt, randm, lin, assigned, xlen, ylen)
struct vtx_data **graph;	/* data structure for graph */
int nvtxs;			/* number of vertices in full graph */
int nedges;			/* number of edges in full graph */
int ndims;			/* number of eigenvectors (2^d sets) */
int ndims_tot;			/* number of cuts made in total */
short (*hop_mtx)[MAXSETS];	/* between-set hop cost for KL */
double *vwsqrt;			/* sqrt of vertex weights (length nvtxs+1) */
int spec;			/* run spectral global decomposition? */
int inert;			/* run inertial global decomposition? */
int rcb;			/* run coordinate global decomposition? */
int KL;				/* run Kernighan-Lin local optimization? */
int mediantype;			/* method for partitioning eigenvector */
int mkconnected;		/* check for connectivity & add phantom edges? */
int solver_flag;		/* which eigensolver should I use? */
int coarse_flag;		/* should I use multilevel techniques? */
int vmax;			/* if so, how many vertices to coarsen down to? */
double eigtol;			/* tolerance on eigenvectors */
int igeom;			/* geometric dimension for inertial method */
float **coords;			/* coordinates for inertial method */
short *assignment;		/* set number of each vtx (length n) */
double *goal;			/* desired set sizes */
int scatt;			/* scattered global decomposition */
int randm;			/* random global decomposition */
int lin;			/* linear global decomposition */
int assigned;       		/* bit assigned so far */
int xlen;			/* x length of sub processor array */
int ylen;			/* y length of sub processor array */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int ARCH_GOAL;	/* trace the execution of the code */
   extern int DEBUG_TRACE;	/* trace the execution of the code */
   extern int DEBUG_RECURSE;/* trace the path of recursion */
   struct vtx_data **subgraph;	/* data structure for subgraph */
   int *loc2glob;		/* mapping from subgraph to graph numbering */
   int *glob2loc;		/* mapping from graph to subgraph numbering */
   short *degree;		/* degrees of graph vertices from a subgraph */
   int setsize[MAXSETS];	/* sizes of sets created by division */
   double merged_goal[MAXSETS];	/* sizes of sets at this partition level */
   double *subgoal;		/* sizes of sets needed to recurse */
   double sub_vwgt_sum;		/* sum of subgraph vertex weights */
   int maxsize;			/* size of largest subgraph */
   int nsets;			/* number of sets created by a division */
   /* int nsets_tot; */		/* total sets to divide subgraph into */
   int set;			/* set that a vertex is assigned to */
   int sub_ndims;		/* number of dimensions left for subgraphs */
   int subnvtxs;		/* number of vertices in subgraph */
   int subnedges;		/* number of edgess in subgraph */
   double *subvwsqrt;		/* vwsqrt array for subgraph */
   short *subassign;		/* set assignments for subgraph */
   float **subcoords;		/* coordinates for subgraph */
   int using_vwgts;		/* are vertex weights being used? */
   int i, j, k, l;			/* loop counter */
   int next_set;        /* next set number generated */
   int maxgoal;        /* goal with highest number of processor numbers */
   int set_nr;    
   int subassigned;
   int part_type;
   int recurse_flag;
   int array2proc[MAXSETS];
   int xlens[MAXSETS];
   int ylens[MAXSETS];
   int set_count[MAXSETS];
   short  *ready;
   short  suba;
   void divide(), make_maps(), make_subgraph(), make_subgoal_rec();
   void make_subvector(), make_subgeom(), remake_graph();
   void merge_assignments();

   if (DEBUG_TRACE > 0) { 
     {char buf[150]; sprintf(buf,"Entering recurse with nvtxs=%d, nedges=%d, nsets_tot=%d, assigned=%d\n", nvtxs, nedges, ndims_tot, assigned);UserWrite(buf);}
   }

	recurse_flag=FALSE;

   /* Compute some simple parameters. */
   using_vwgts = (vwsqrt != NULL);

   if (ndims == 0) return;

	nsets = 1<<ndims;

	if (ARCH_GOAL==0)
	{
		if (xlen*ylen <= 1) return;
		/* initialize partition information */
		for (i=0; i<MAXSETS; i++)
			merged_goal[i]=array2proc[i]=xlens[i]=ylens[i]=set_count[i]=0;
	
		/* determine part_type and information about subpartitions */
		part_type = partition_type(xlen,ylen,ndims,array2proc,xlens,ylens);
	}
       
   for (i=1; i<=nvtxs; i++) assignment[i] = 0;

   if (nvtxs <= 1) return;

   /* Construct desired set sizes for this division step. */
   /* Note: This is HYPERCUBE/MESH dependent.  You want to combine */
   /* different sets on the different architectures. */
	maxgoal = 0;
	if (ARCH_GOAL == 0)
	{
		for (j=0; j<xlen; j++)
			for (l=0; l<ylen; l++)
			{
				set_nr = map(j,l,ndims,part_type,xlens,ylens);
				merged_goal[set_nr] += goal[l*xlen+j];
				set_count[set_nr]++;
			}
		for (k=0; k<nsets; k++)
			maxgoal = max(maxgoal,set_count[k]);
	}	
	else 
	{
		for (i=0; i<nsets; i++) 
		{
			sub_vwgt_sum = k = 0;
			for (j=(i<<ndims_tot)|assigned; j<nsets_glob; j+=(nsets<<ndims_tot)) 
			{
				sub_vwgt_sum += goal[j>>ndims_tot];
				k++;
			}
			/* remember goal with highest number of processor numbers */ 
			if (maxgoal<k) maxgoal = k;
			merged_goal[i] = sub_vwgt_sum;
		}
	}

   /* Perform a single division step. */
   divide(graph, nvtxs, nedges, ndims, vwsqrt, spec, inert, rcb, KL, mediantype,
	  mkconnected, solver_flag, coarse_flag, vmax, eigtol, hop_mtx, igeom,
	  coords, assignment, merged_goal, scatt, randm, lin, part_type);
   if (!MEM_OK) return;

	if (ARCH_GOAL == 0)
	{
		/* determine whether to recurse */  
		if (xlen*ylen>=2) recurse_flag = TRUE;
		else recurse_flag = FALSE;
	}
	else 
	{
		/* compute minimal setnumber that would be generated by next */
		/* minimal (=bisection) subdivision step                     */
		next_set = 1<<(ndims+ndims_tot)|assigned;

		if (!((nsets_glob-1) >= next_set))
			recurse_flag = FALSE;
		else
			recurse_flag = TRUE;
	}


   /* no more recursion step,then print out setnumbers */
   if (recurse_flag==FALSE)
   {
#    ifdef __DEBUG__
     if (DEBUG_RECURSE > 0)
     {
	   for (set=0; set<nsets; set++)
          {char buf[150]; sprintf(buf,"Recursion terminating for set=%d\n", assigned);UserWrite(buf);}
     }
#    endif 
   }

   if (recurse_flag==TRUE) {		/* Need to recurse */
      /* Find size of largest subgraph. */
      for (i=0; i<nsets; i++) setsize[i] = 0;
      for (i=1; i<=nvtxs; i++)  ++setsize[assignment[i]];
      maxsize = 0;
      for (i=0; i<nsets; i++) if (setsize[i] > maxsize) maxsize = setsize[i];

      /* Allocate space for recursing based on size of largest subgraph. */
      subgoal = (double *) (MEM_OK = smalloc((unsigned) maxgoal*sizeof(double));
      if (!MEM_OK) return;
      glob2loc = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
      if (!MEM_OK) return;
      loc2glob = (int *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(int));
      if (!MEM_OK) return;
      if (graph != NULL) {
	 subgraph =(struct vtx_data **) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(struct vtx_data *));
         if (!MEM_OK) return;
         degree = (short *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(short));
         if (!MEM_OK) return;
      }
      else {
	 subgraph = NULL;
      }
      subassign = (short *) (MEM_OK = smalloc((unsigned) (maxsize+1)*sizeof(short));
      if (!MEM_OK) return;
      ready = (short *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(short));
      if (!MEM_OK) return;
      for (i=0; i<=nvtxs; i++)  ready[i] = FALSE;
      if (!using_vwgts) subvwsqrt = NULL;
      else 
      {
            subvwsqrt = (double *) (MEM_OK = smalloc((unsigned)(maxsize+1)*sizeof(double));
            if (!MEM_OK) return;
      }

      subcoords = NULL;
      if (inert || rcb) {
	 subcoords = (float **) (MEM_OK = smalloc((unsigned)3*sizeof(float *));
     if (!MEM_OK) return;
	 subcoords[0] = subcoords[1] = subcoords[2] = NULL;
	 subcoords[0] = (float *) (MEM_OK = smalloc((unsigned)(maxsize+1)*sizeof(float));
     if (!MEM_OK) return;
	 if (igeom > 1) subcoords[1] = (float *) (MEM_OK = smalloc((unsigned)(maxsize+1)*sizeof(float));
     if (!MEM_OK) return;
	 if (igeom > 2) subcoords[2] = (float *) (MEM_OK = smalloc((unsigned)(maxsize+1)*sizeof(float));
     if (!MEM_OK) return;
      }

      /* Partition each of the subgraphs */
      for (set=0; set<nsets; set++) {

	/* determine partitioning dimension */
	if (ARCH_GOAL==0)
	{
		sub_ndims = compute_sub_ndims(xlens[set],ylens[set],ndims);
	}
	else
	{
		sub_ndims = ndims;
		/* loop until the maximal setnumber generated by the this partitioning */
		/* dimension is a processor number which is needed                     */
		next_set = ((((1<<sub_ndims)-1)<<ndims|set)<<ndims_tot|assigned);
		while (nsets_glob<=next_set && sub_ndims>0)
		{
			sub_ndims-- ;
			/* max. setnumber = max. generated bitcombination + bits set by this */
			/*					 division step + bits set before is division step */     
			next_set = ((((1<<sub_ndims)-1)<<ndims|set)<<ndims_tot|assigned);
		}
	}

	 /* generated setnumbers would exceed number of processors needed */
	 if (sub_ndims==0) 
     {
#      ifdef __DEBUG__
       if (DEBUG_RECURSE > 0)
       {
        {char buf[150]; sprintf(buf,"Recursion terminating for set=%d\n", set<<ndims_tot|assigned);UserWrite(buf);}
       }
#      endif
		for (i=1; i<=nvtxs; i++)
		{
			suba = 0;
			if (assignment[i]==set && ready[i]==FALSE) 
			{
				assignment[i] = suba2a(NULL, NULL, i, &suba, ndims, 
				                       xlen, ylen, set, xlens, ylens, part_type);
				ready[i] = TRUE;
			}
		}
       continue;
     }

#    ifdef __DEBUG__
     if (DEBUG_RECURSE > 0)
     {
       {char buf[150]; sprintf(buf,"Recursion continues: next_set=%d sub_ndims=%d\n",next_set,sub_ndims);UserWrite(buf);} 
     }
#    endif 

     /* determine subgraph size */
	 subnvtxs = setsize[set];

	 /* Construct mappings between local and global vertex numberings. */
	 make_maps(assignment, nvtxs, set, glob2loc, loc2glob, ready);

	 /* Form the subgraph in our graph format. */
	 if (graph != NULL) make_subgraph(graph, subgraph, subnvtxs, &subnedges,
				  assignment, set, glob2loc, loc2glob, degree, ready);
         else subnedges = 0;		/* Otherwise some output is garbage */

	 if (!using_vwgts) sub_vwgt_sum = subnvtxs;
	 else {
	    sub_vwgt_sum = 0;
	    for (i=1; i<=subnvtxs; i++) sub_vwgt_sum += subgraph[i]->vwgt;
	 }

	 make_subgoal_rec(goal, subgoal, ndims, ndims_tot, set, sub_vwgt_sum, 
	                  assigned, xlen, ylen, xlens, ylens, part_type);

	 /* Condense the relevant vertex weight array. */
         if (using_vwgts) make_subvector(vwsqrt, subvwsqrt, subnvtxs, loc2glob);

	 if (inert || rcb) make_subgeom(igeom, coords, subcoords, subnvtxs, loc2glob);

	if (ARCH_GOAL==0)
	{
		subassigned = 0;
	}
	else
	{
		subassigned = (set<<ndims_tot)|assigned;
	}
         recurse(subgraph, subnvtxs, subnedges, sub_ndims, ndims_tot+ndims, hop_mtx,
		 subvwsqrt, spec, inert, rcb, KL, mediantype, mkconnected, solver_flag, coarse_flag, vmax, eigtol, igeom, subcoords,
		 subassign, subgoal, scatt, randm, lin, subassigned,
		 xlens[set],ylens[set]);

         if (!MEM_OK) return;

	 /* Undo the subgraph construction. */

	 if (graph != NULL) remake_graph(subgraph, subnvtxs, loc2glob, degree);

	 /* Merge the subgraph partitioning with the graph partitioning. */
         merge_assignments_rec(assignment, subassign, ndims, subnvtxs, loc2glob, 
		                       set, xlen, ylen, xlens, ylens, part_type, ready);
      }
      /* Free everything allocated for subgraphs. */
      sfree((char *) ready);
      sfree((char *) subassign);
      sfree((char *) loc2glob);
      sfree((char *) glob2loc);
      if (graph != NULL) sfree((char *) degree);
      if (subgraph != NULL) sfree((char *) subgraph);
      if (subvwsqrt != NULL) sfree((char *) subvwsqrt);
      if (subcoords != NULL) {
         if (subcoords[0] != NULL) sfree((char *) subcoords[0]);
         if (subcoords[1] != NULL) sfree((char *) subcoords[1]);
         if (subcoords[2] != NULL) sfree((char *) subcoords[2]);
	 sfree((char *) subcoords);
      }
      sfree((char *) subgoal);
   }
}
