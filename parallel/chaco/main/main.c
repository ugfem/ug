// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/params.h"
#include "../main/structs.h"

extern Heap   *heap;     /* pointer to heap of multigrid */
double *MEM_OK;   /* variable for memory overflow exeception */

int main()
{
   extern double EIGEN_TOLERANCE;	/* tolerance for eigen calculations */
   extern int MAKE_CONNECTED;	/* connect graph if using spectral method? */
   extern int MAPPING_TYPE;	/* how to map from evecs to partition */
   extern int OUTPUT_ASSIGN;	/* whether to write assignments to file */
   extern int DEBUG_MEMORY;	/* debug memory allocation and freeing? */
   extern double input_time;	/* times data file input */
   extern long RANDOM_SEED;	/* seed for random number generators */
   FILE *fin;			/* input file */
   FILE *fingeom;		/* geometry input file (for inertial method) */
   FILE *outfile;		/* file for outputing run results */
   float *x, *y, *z;		/* coordinates for inertial method */
   int *start;			/* start of edge list for each vertex */
   int *adjacency;		/* edge list data */
   int *ewgts;			/* weights for all edges */
   int *vwgts;			/* weights for all vertices */
   short *assignment;		/* set number of each vtx (length nvtxs+1) */
   int *free_ptr;		/* adjacency ptr to free copy of graph */
   double eigtol;		/* tolerance in eigenvector calculation */
   int nvtxs;			/* number of vertices in graph */
   int ndims;			/* dimension of recursive partitioning */
   int ndims_tot;		/* total number dimensions to partition into */
   int spec;			/* run spectral global decomposition? */
   int inert;			/* run inertial global decomposition? */
   int KL;			/* run Kernighan-Lin local optimization? */
   long seed;			/* for random graph mutations */
   int mediantype;		/* method for partitioning eigenvector */
   int mkconnected;		/* check connectivity & add phantom edges? */
   int solver_flag;		/* which eigensolver should I use? */
   int coarse_flag;		/* should I use multilevel techniques? */
   int vmax;			/* if so, how many vertices to coarsen down to? */
   int igeom;			/* geometry dimension if inertial method */
   char graphname[NAME_LENGTH];	/* name of graph input file */
   char geomname[NAME_LENGTH];	/* name of geometry input file */
   char assignname[NAME_LENGTH];/* name of output file */
   char *geomptr, *assignptr;	/* names or null pointers */
   int scatt;			/* run scatt global partitioner? */
   int randm;			/* run randm global partitioner? */
   int lin;			/* run linear global partitioner? */
   double time;			/* timing marker */

   double seconds();		/* returns elapsed time in seconds */
   int interface();
   void input_queries(), input_graph();
   void input_geom(), smalloc_stats();

   time = seconds();

   /*malloc_debug(2);*/

   input_queries(&fin, &fingeom, graphname, geomname, &ndims, &ndims_tot,
		 assignname, &spec, &inert, &KL, &solver_flag, &coarse_flag,
		 &vmax, &scatt, &randm, &lin, &outfile);
   if (OUTPUT_ASSIGN > 0) assignptr = assignname;
   else assignptr = NULL;

   input_graph(fin, graphname, &start, &adjacency, &nvtxs, &vwgts, &ewgts);
   if (!MEM_OK) return;

   if (inert) {
      input_geom(fingeom, geomname, nvtxs, &igeom, &x, &y, &z);
      if (!MEM_OK) return;
      geomptr = geomname;
   }
   else {
      geomptr = NULL;
      igeom = 0;
      fingeom = NULL;
      x = y = z = NULL;
   }

   input_time = seconds() - time;

   eigtol = EIGEN_TOLERANCE;
   mediantype = MAPPING_TYPE;
   seed = RANDOM_SEED;
   assignment = (short *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(short));
   if (!MEM_OK) return;

   /* Initialize free_ptr to delete one copy of graph. */
   /* To not delete anything, set free_ptr = NULL. */
   free_ptr = adjacency;

   /* Initialize the connectivity flag. */
   if (spec == TRUE) mkconnected = MAKE_CONNECTED;
   else mkconnected = FALSE;

   interface(start, adjacency, nvtxs, vwgts, ewgts, ndims, ndims_tot,
        spec, inert, KL, mediantype, mkconnected, solver_flag,
        coarse_flag, vmax, eigtol, seed, igeom, x, y, z, assignment, free_ptr, 
        scatt, randm, lin, graphname, geomptr, assignptr, outfile);
   if (!MEM_OK) return;

   if (inert) {
      if (z != NULL) sfree((char *) z);
      if (y != NULL) sfree((char *) y);
      if (x != NULL) sfree((char *) x);
   }
   sfree((char *) assignment);

   if (DEBUG_MEMORY > 0) {
     {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
      smalloc_stats();
  }

  {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
   return(0);
}
