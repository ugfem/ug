// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include "../main/defs.h"
#include "../main/params.h"
#include "../main/structs.h"
#include "ppif.h"

double *SQRTS;			/* precomputed square roots for efficiency */

double DOUBLE_EPSILON;		/* machine precision */
double DOUBLE_MAX;		/* largest double precision value */


int balance(graph, nvtxs, nedges, ndims, nsets_tot, using_vwgts, spec, inert,
       KL, mediantype, mkconnected, solver_flag, coarse_flag, vmax,
       eigtol, seed, igeom, coords, assignment, goal, scatt, randm, lin,
       graphname, geomname, assignname, outfile, dimx, dimy)
struct vtx_data **graph;	/* data structure for graph */
int nvtxs;			/* number of vertices in full graph */
int nedges;			/* number of edges in graph */
int ndims;			/* number of eigenvectors (2^d sets) */
int nsets_tot;			/* total number of sets to divide into */
int using_vwgts;		/* are vertex weights being used? */
int spec;			/* run spectral global decomposition? */
int inert;			/* run inertial method global decomposition? */
int KL;				/* run Kernighan-Lin local optimization? */
int mediantype;			/* method for partitioning eigenvector */
int mkconnected;		/* check for connectivity & add phantom edges? */
int solver_flag;		/* which eigensolver should I use? */
int coarse_flag;		/* should I apply multilevel techniques? */
int vmax;			/* if so, how many vtxs to coarsen down to */
double eigtol;			/* tolerance on eigenvectors */
long seed;			/* for random graph mutations */
int igeom;			/* geometry dimension if using inertial method */
float **coords;			/* coordinates of vertices if used */
short *assignment;		/* set number of each vtx (length n) */
double *goal;			/* desired sizes for each set */
int scatt;			/* run scattered global partitioner? */
int randm;			/* run random global partitioner? */
int lin;			/* run linear global partitioner? */
char *graphname, *geomname;	/* names of input files */
char *assignname;		/* name of assignment output file */
FILE *outfile;			/* in which to print output metrics */
int dimx;				/* x dimension of processor array */
int dimy;				/* y dimension of processor array */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int CHECK_INPUT;	/* check the input for consistency? */
   extern int ECHO;		/* print out input and user parameters? */
   extern int OUTPUT_METRICS;	/* whether to tally and display cost metrics */
   extern int OUTPUT_ASSIGN;	/* print assignment to a file? */
   extern int NSQRTS;		/* number of square roots to precompute */
   extern int KL_METRIC;	/* KL interset cost: 1=>cuts, 2=>hops */
   extern int RAND_MAXIMUM;	/* largest value returnable by rand() */
   extern double check_input_time;	/* time spent checking input */
   extern double partition_time;/* time spent partitioning graph */
   extern double count_time;	/* time spent evaluating the answer */
   extern double print_assign_time;	/* time spent writing output file */
   extern int nsets_glob;       /* total number of sets to generate */
   struct vtx_data **graph2;	/* data structure for graph */
   short hop_mtx[MAXSETS][MAXSETS];	/* between-set hop cost for KL */
   double *vwsqrt;		/* sqrt of vertex weights (length nvtxs+1) */
   double time;
   double temp;			/* dummy parameter */
   int nsets;			/* number of sets created by each divide */
   int bits;			/* used in computing hops */
   int flag;			/* return code from checking input */
   int i, j, k;			/* loop counters */
   double seconds();
   void srand48();
   double sqrt();
   int check_input();
   void setrandom(), makevwsqrt(), recurse(), countup();
   void reflect_input(), machine_params(), assign_out();

   /* Check the input for inconsistencies. */
   if (CHECK_INPUT) {
      time = seconds();
      flag = check_input(graph, nvtxs, nedges, spec, inert, KL,
	 mediantype, solver_flag, coarse_flag, vmax, 1, 2,
	 igeom, coords, scatt, randm, lin, graphname);
      check_input_time = seconds() - time;
      if (!flag) {
        {char buf[150]; sprintf(buf," ERROR: Inconsistent input.\n");UserWrite(buf);}
         return(1);
      }
   }

   if (ECHO != 0) {
      reflect_input(nvtxs, spec, inert, KL, solver_flag, coarse_flag, vmax,
	 ndims, nsets_tot, igeom, scatt, randm, lin, graphname, geomname,
	 assignname, outfile);
   }

   time = seconds();

   /* Perform some one-time initializations. */
   srand48(seed);
   setrandom(seed);
   machine_params(&DOUBLE_EPSILON, &temp, &DOUBLE_MAX, &RAND_MAXIMUM);
   nsets = (1<<ndims);
   SQRTS = (double *) (MEM_OK = smalloc((unsigned) (NSQRTS+1)*sizeof(double));
   if (!MEM_OK) return;
   for (i=1; i<=NSQRTS; i++) SQRTS[i] = sqrt((double) i);

   if (using_vwgts) {
      vwsqrt = (double *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(double));
      if (!MEM_OK) return;
      makevwsqrt(vwsqrt, graph, nvtxs);
   }
   else vwsqrt = NULL;

   /* Initialize cost function for KL-spiff */
   if (KL) {
      for (i=0; i<nsets; i++) {
         for (j=0; j<nsets; j++) {
      	    if (KL_METRIC == 2) {		/* Count hypercube hops */
               hop_mtx[i][j] = 0;
      	       bits = i^j;
	       for (k=bits; k; k>>=1) if (k&1) ++hop_mtx[i][j];
	    }
      	    if (KL_METRIC == 3) {		/* Count array hops */
		if (i==j)
               		hop_mtx[i][j] = 0;
		else 
               		hop_mtx[i][j] = absval(i%DimX - j%DimX) + absval(i/DimX - j/DimX);
		}
	    else if (KL_METRIC == 1) {	/* Count cut edges */
	       hop_mtx[i][j] = (i==j) ? 0 : 1;
	    }
	 }
      }
   }

   graph2 = graph;
   if (inert && !KL && !using_vwgts) graph2 = NULL;
   recurse(graph2, nvtxs, nedges, ndims, 0, hop_mtx, vwsqrt, spec, inert,
         KL, mediantype, mkconnected, solver_flag, coarse_flag,
         vmax, eigtol, igeom, coords, assignment, goal, scatt, randm, lin, 0, 
         dimx, dimy);
   if (!MEM_OK) return;

   partition_time += seconds() - time;

   /* Count the cuts and hops. */
   if (OUTPUT_METRICS > 0) {
      time = seconds();
      if (graph != NULL) 
      { 
          countup(graph, nvtxs, assignment, ndims, nsets_glob,
				 OUTPUT_METRICS, outfile);
          if (!MEM_OK) return;
      }

      count_time += seconds() - time;
   }

   if (OUTPUT_ASSIGN > 0) {
      time = seconds();
      assign_out(nvtxs, assignment, assignname);
      print_assign_time += seconds() - time;
   }

   if (vwsqrt != NULL ) sfree((char *) vwsqrt);
   if (SQRTS != NULL) sfree((char *) SQRTS);

   return(0);
}
