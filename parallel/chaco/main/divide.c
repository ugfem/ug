// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"


void divide(graph, nvtxs, nedges, ndims, vwsqrt, spec, inert, KL, 
       mediantype, mkconnected, solver_flag, coarse_flag, vmax, eigtol,
       hop_mtx, igeom, coords, assignment, goal, scatt, randm, lin, part_type)
struct vtx_data **graph;	/* graph data structure */
int nvtxs;			/* number of vertices in graph */
int nedges;			/* number of edges in graph */
int ndims;			/* number of eigenvectors (2^d sets) */
double *vwsqrt;			/* sqrt of vertex weights (length nvtxs+1) */
int spec;			/* run spectral global decomposition? */
int inert;			/* run inertial global decomposition? */
int KL;				/* run Kernighan-Lin local optimization? */
int mediantype;			/* method for partitioning eigenvector */
int mkconnected;		/* check connectivity & add phantom edges? */
int solver_flag;		/* which eigensolver should I use? */
int coarse_flag;		/* should I use multilevel techniques? */
int vmax;			/* if so, # vertices to coarsen down to */
double eigtol;			/* tolerance on eigenvectors */
short (*hop_mtx)[MAXSETS];	/* between-set hop costs for KL */
int igeom;			/* geometry dimension for inertial method */
float **coords;			/* coordinates for inertial method */
short *assignment;		/* set number of each vtx (length n) */
double *goal;			/* desired set sizes */
int scatt;			/* indicates scattered global decomposition */
int randm;			/* indicates random global decomposition */
int lin;			/* indicates lin global decomposition */
int part_type;
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int DEBUG_CONNECTED;	/* debug flag for connected components */
   extern int DEBUG_EVECS;	/* debug flag for eigenvector calculation */
   extern int COARSE_NLEVEL_KL;	/* how often to invoke KL while uncoarsening */
   struct connect_data *cdata;	/* data for enforcing connectivity */
   double *yvecs[MAXDIMS+1];	/* space for pointing to eigenvectors */
   double evals[MAXDIMS+1];	/* corresponding eigenvalues */
   double maxdeg;		/* maximum weighted degree of a vertex */
   int *active;			/* keeping track of vtxs in BFS (length nvtxs)*/
   short *mark;			/* for marking vtxs in BFS (length nvtxs+1) */
   int using_vwgts;		/* are vertex weights being used? */
   int using_ewgts;		/* are edge weights being used? */
   int nsets;			/* number of sets being divided into */
   int vwgt_max;		/* largest vertex weight */
   int i;			/* loop counters */
   int simple_type;		/* which type of simple partitioning to use */
   
   void make_connected(), print_connected(), coarsen_kl(), eigensolve();
   void assign(), make_unconnected(), inertial(), klspiff(), simple_part();
   double find_maxdeg();

   using_vwgts = (vwsqrt != NULL);
   if (graph == NULL) using_ewgts = FALSE;
   else using_ewgts = (graph[1]->ewgts != NULL);
   nsets = (1<<ndims);

   if (using_vwgts) {
      vwgt_max = 0;
      for (i=1; i<=nvtxs; i++) {
	 if (graph[i]->vwgt > vwgt_max) vwgt_max = graph[i]->vwgt;
      }
   }
   else vwgt_max = 1;

   /* Perform one of the global partitionings on this sub-graph*/
   if(spec) {
      if (DEBUG_EVECS > 0){char buf[150]; sprintf(buf,"Eigen pair tolerance %g\n",eigtol);UserWrite(buf);}

      active = (int *) (MEM_OK = smalloc((unsigned) nvtxs*sizeof(int));
      if (!MEM_OK) return;
      if (coarse_flag != 2) for (i=1; i<=ndims; i++)  {
	 yvecs[i] = (double *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(double));
      if (!MEM_OK) return;
      }

# ifdef __DEBUG__
/* debug_graph */
{
extern int DEBUG_GRAPH;
	if (DEBUG_GRAPH>0)
	{
		int flag;
        char buf[100];

		flag = check(graph,nvtxs,nedges);
		if (flag == FALSE)
		{
			sprintf(buf,"FATAL: check_graph returned FALSE in divide() before mk_cnt\n");
			UserWrite(buf);
			sprintf(buf,"nvtxs=%d,nedges=%d,ndims=%d\n",nvtxs,nedges,ndims);
			UserWrite(buf);
		}
        else
        {
            sprintf(buf,"OK: check_graph = TRUE in divide() before mk_cnt,nvtxs=%d,nedges=%d\n",nvtxs,nedges);
            UserWrite(buf);
        }

	}
}
# endif


      if (mkconnected) {
	 if (coarse_flag != 2) {
	    /* If doing multi-level KL, this happens on coarse graph. */
	    mark = (short *) &(yvecs[1][0]);
      	    make_connected(graph, nvtxs, &nedges, mark, active, &cdata);
            if (!MEM_OK) return;
	    if (DEBUG_CONNECTED > 0) {
   	      {char buf[150]; sprintf(buf,"Enforcing connectivity\n");UserWrite(buf);}
	       print_connected(cdata);
	    }
	 }
      }
      else {
	 if (DEBUG_CONNECTED > 0) {
   	   {char buf[150]; sprintf(buf,"Not enforcing connectivity\n");UserWrite(buf);}
	 }
      }

# ifdef __DEBUG__
/* debug_graph */
{
extern int DEBUG_GRAPH;
	if (DEBUG_GRAPH>0)
	{
		int flag;
        char buf[100];

		flag = check(graph,nvtxs,nedges);
		if (flag == FALSE)
		{
			sprintf(buf,"FATAL: check_graph returned FALSE in divide() after mk_cnt\n");
			UserWrite(buf);
			sprintf(buf,"nvtxs=%d,nedges=%d,ndims=%d\n",nvtxs,nedges,ndims);
			UserWrite(buf);
		}
        else
        {
			sprintf(buf,"OK: check_graph = TRUE in divide() after mk_cnt,nvtxs=%d,nedges=%d\n",nvtxs,nedges);
			UserWrite(buf);
        }

	}
}
# endif 


      if (coarse_flag == 2) {	/* Divide coarse graph, smooth with KL. */
         coarsen_kl(graph, nvtxs, nedges, ndims, COARSE_NLEVEL_KL, 0, vmax,
		    eigtol, assignment, nsets, hop_mtx, mediantype,
		    mkconnected, using_ewgts, using_vwgts, solver_flag,
		    goal, vwgt_max, FALSE, part_type);
         if (!MEM_OK) return;
      }
      else {
         maxdeg = find_maxdeg(graph, nvtxs);
         eigensolve(graph, nvtxs, nedges, maxdeg, yvecs, evals, ndims, eigtol,
		   vwsqrt, solver_flag, coarse_flag, vmax, igeom, coords,
		   using_ewgts, using_vwgts);
         if (!MEM_OK) return;

         assign(graph, yvecs, nvtxs, ndims, vwsqrt, assignment, active,
	        mediantype, goal, vwgt_max, part_type);
         if (!MEM_OK) return;
         for (i=1; i<=ndims; i++) sfree((char *) yvecs[i]);
      }

# ifdef __DEBUG__
/* debug_graph */
{
extern int DEBUG_GRAPH;

    /* debug_graph */
    if (DEBUG_GRAPH>0)
    {		
		int flag;
        char buf[100];

        flag = check(graph,nvtxs,nedges);
        if (flag == FALSE)
        {
            sprintf(buf,"FATAL: check_graph returned FALSE in divide before un_cnt\n");
            UserWrite(buf);
            sprintf(buf,"nvtxs=%d,nedges=%d,ndims=%d\n",nvtxs,nedges,ndims);
            UserWrite(buf);
        }
        else
        {
            sprintf(buf,"OK: check_graph = TRUE in divide() before un_cnt,nvtxs=%d,nedges=%d\n",nvtxs,nedges);
            UserWrite(buf);
        }

    }
}
# endif

      if (mkconnected) {
	 if (coarse_flag != 2) make_unconnected(graph, &nedges, &cdata);
      }

# ifdef __DEBUG__
/* debug_graph */
{
extern int DEBUG_GRAPH;

    /* debug_graph */
    if (DEBUG_GRAPH>0)
    {
		int flag;
        char buf[100];

        flag = check(graph,nvtxs,nedges);
        if (flag == FALSE)
        {
            sprintf(buf,"FATAL: check_graph returned FALSE in divide after un_cnt\n");
            UserWrite(buf);
            sprintf(buf,"nvtxs=%d,nedges=%d,ndims=%d\n",nvtxs,nedges,ndims);
            UserWrite(buf);
        }
        else
        {
            sprintf(buf,"OK: check_graph = TRUE in divide() after un_cnt,nvtxs=%d,nedges=%d\n",nvtxs,nedges);
            UserWrite(buf);
        }
    }
}
# endif

      sfree ((char *) active);
   }
   else if (inert) {
      inertial(graph, nvtxs, ndims, igeom, coords, assignment, goal, 
	           using_vwgts, part_type);
      if (!MEM_OK) return;
   }
   else if (scatt) {
      simple_type = 1;
      simple_part(graph, nvtxs, assignment, nsets, simple_type, goal);
      if (!MEM_OK) return;
   }
   else if (randm) {
      simple_type = 2;
      simple_part(graph, nvtxs, assignment, nsets, simple_type, goal);
      if (!MEM_OK) return;
   }
   else if (lin) {
      simple_type = 3;
      simple_part(graph, nvtxs, assignment, nsets, simple_type, goal);
      if (!MEM_OK) return;
   }

   /* Perform a local refinement, if specified, on the global partitioning. */
   if (KL && coarse_flag != 2)  {
      /* Find the maximum weighted degree of a vertex. */
      if (!spec) {
         maxdeg = find_maxdeg(graph, nvtxs);
      }
      /* If coarse_flag = 2 already did KL as part of multilevel KL/spec. */
      klspiff(graph, nvtxs, assignment, nsets, hop_mtx, goal, vwgt_max, maxdeg);
      if (!MEM_OK) return;
   }
}
