// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*                                                                          */
/* File:      interface.c                                                   */
/*                                                                          */
/* Purpose:   interface for clustered partitioning using Chaco's            */
/*            partitioning methods in lb4.c of ugp.                         */
/*                                                                          */
/*            Chaco's partioning schemes are one of 5 spectral methods,     */
/*            a high-end multilevel strategy, inertial partitioning or      */
/*            one of 3 simple method's. The partioning may be optimized     */
/*            by a local KL method. Altough one can choose betwween         */
/*            bisection, quadrsection and octasection.                      */
/*            Currently n x m  processor array configurations are possible. */
/*                                                                          */
/*            This file contains multiple references to files, which are    */
/*            included in Chaco's original source code. Probably they are   */
/*            of interest when changing the code, because Chaco consists    */
/*            of a more than a hundred source files.                        */
/*                                                                          */
/* Author:    Stefan Lang                                                   */
/*            Institut fuer Mathematische Maschinen und                     */
/*            Datenverarbeitung III                                         */
/*            Universitaet Erlangen-Nuernberg                               */
/*            Martensstrasse 3                                              */
/*            91058 Erlangen                                                */
/*                                                                          */
/* History:   3 Jan 94 begin                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"
#include "../main/params.h"
#include "../lb4.h"
#include "gm.h"

/* defines to extract flags from  mode */

#define CHECK_INPUT_MASK      01
#define OUTPUT_METRICS_MASK   06
 
#define OUTPUT_ASSIGN_MASK   010
#define DEBUG_MEMORY_MASK    060
#define OUTPUT_TIME_MASK    0300

#define CHECK_INPUT_SHIFT      0
#define OUTPUT_METRICS_SHIFT   1 
#define OUTPUT_ASSIGN_SHIFT    3
#define DEBUG_MEMORY_SHIFT     4
#define OUTPUT_TIME_SHIFT      6

#define GET_CHECK_INPUT(x)      (((x)&CHECK_INPUT_MASK)>>CHECK_INPUT_SHIFT)
#define GET_OUTPUT_METRICS(x)   (((x)&OUTPUT_METRICS_MASK)>>OUTPUT_METRICS_SHIFT)
#define GET_OUTPUT_ASSIGN(x)   (((x)&OUTPUT_ASSIGN_MASK)>>OUTPUT_ASSIGN_SHIFT)
#define GET_DEBUG_MEMORY(x)     (((x)&DEBUG_MEMORY_MASK)>>DEBUG_MEMORY_SHIFT)
#define GET_OUTPUT_TIME(x)      (((x)&OUTPUT_TIME_MASK)>>OUTPUT_TIME_SHIFT)

/* global variables defined in other modules */
extern int nsets_glob;           /* total number of processors to partition into */

/* global variables for this source file. Static! */
static double start_time;		/* time routine is entered */

/* reformats the graph in adding edges to neighbors, which have not the */
/* same cluster depth and are therefore not partitioned here, at the    */
/* end of the edge list. Thus the old graph can be regained.            */

void reformat_graph (CLUSTER **clusters, struct vtx_data **graph, int nvtxs,
                     int *nedges, int *glob2loc, int *loc2glob, int using_ewgts)
{
	struct vtx_data *gptr;
	float *fptr;
	int *iptr,*pvtx;
	double ewgtsum;
	int subnedges;
	int neighbor;
	int newnedges;
	int tempvtx;
	int tempwgt;
    int i,j;    
        
	subnedges = 0;

	/* begin part of Chaco's make_maps() */
	for (i=1; i<=nvtxs; i++)
	{
		pvtx = graph[i]->edges;
		if (glob2loc != NULL) glob2loc[*pvtx] = i;
		loc2glob[i] = *pvtx;
		*pvtx = i;
	}
	/* end part of Chaco's make_maps() */

	for (i=1; i<=nvtxs; i++)
	{
		gptr = graph[i];
		ewgtsum = 0;
		newnedges = gptr->nedges;
		/* Move all deactivated edges to the end of the list. */
		iptr = gptr->edges + 1;
		if (using_ewgts) fptr = gptr->ewgts + 1;
		for (j=1; j<newnedges; ) 
		{
			neighbor = *iptr;
			if (glob2loc[neighbor] != 0) 
			{   /* Keep vertex in edge list. */
				*iptr= glob2loc[neighbor];
				if (using_ewgts) ewgtsum += *fptr++;
				j++;
				iptr++;
			}
			else 
			{   /* Move vertex to back of edge list. */
				--newnedges;
				tempvtx = gptr->edges[newnedges];
				gptr->edges[newnedges] = neighbor;
				*iptr = tempvtx;
				if (using_ewgts) 
				{
					tempwgt = gptr->ewgts[newnedges];
					gptr->ewgts[newnedges] = *fptr;
					*fptr = tempwgt;
				}
			}
		}

	gptr->nedges = newnedges;
	subnedges += newnedges;
	if (using_ewgts) gptr->ewgts[0] = -ewgtsum;
	}
	*nedges = (subnedges-nvtxs)/2;

}

/****************************************************************************/
/*                                                                          */
/* Function:  interface                                                     */
/*                                                                          */
/* Purpose:   interface between ugp and the Chaco partitioning tool         */
/*            on level of clustering.                                       */
/*            interface expects the clusterformat of lb4.h. It reformats    */
/*            the clusters into an Chaco internal graph format (struct      */
/*            vtx_data). Furthermore the Chaco specific parameter are       */
/*            set for controlling e.g. the strategy used for partitioning   */
/*            and a lot of other things like input_checks, dimension of     */
/*            partitioning, which eigensolver to use.                       */
/*            Finially interface does the whole timing of the partitioning  */
/*            process.                                                      */
/*                                                                          */
/* Input:     CLUSTER **clusters    the clusters to partition               */
/*            int     nvtxs         number of clusters to partition         */
/*                                  (length of *clusters and assign)        */
/*            short   *assign       to which processor each cluster gets    */
/*                                  assigned to                             */
/*            double  *goal         how many elements should be assign to   */
/*                                  the processors                          */
/*            int     nsets_tot     number of processors (length of goal)   */
/*            int     glob_method   which partitioning strategy to apply    */
/*            int     eigen         which eigensolver to use                */
/*            int     loc           refine partitioning with local KL?      */
/*            int     dims          dimension of partitioning (number of    */
/*                                  sets created in each part. step 2-8)    */
/*            int     weights       use edge and/or vertex weight?          */
/*            int     coarse        number of vertices of coarse graph      */
/*            int     mode          contains flags for use of several Chaco */
/*                                  ingredients                             */
/*            int     *glob2loc     map of clusters from global to local    */
/*                                  graph numbering scheme                  */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if an error occured                                        */
/*                                                                          */
/****************************************************************************/

int interface (CLUSTER **clusters, int nvtxs, short *assign, double *goal,
				int nsets_tot, int glob_method, int eigen, int loc, int dims, 
				int weights, int coarse, int mode, int *glob2loc, 
                int dimx, int dimy)
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */

	/* variables of Chaco's main() */
	extern double EIGEN_TOLERANCE;/* tolerance for eigen calculations */
	extern int MAKE_CONNECTED;   /* connect graph if using spectral method? */
	extern int MAPPING_TYPE;     /* how to map from evecs to partition */
	extern int OUTPUT_ASSIGN;    /* whether to write assignments to file */
	extern int DEBUG_MEMORY;     /* debug memory allocation and freeing? */
	extern double input_time;    /* times data file input */
   	extern long RANDOM_SEED;     /* seed for random number generators */
	extern int CHECK_INPUT;      /* check graph for consistency? */
	extern int OUTPUT_METRICS;   /* analyse quality of partitioning? */
   	FILE *fin;                   /* input file */
   	FILE *fingeom;               /* geometry input file (for inertial method) */
   	FILE *outfile;               /* file for outputing run results */
   	float *x, *y, *z;            /* coordinates for inertial method */
   	int *start;                  /* start of edge list for each vertex */
   	int *vwgts;                  /* weights for all vertices */
   	short *assignment;           /* set number of each vtx (length nvtxs+1) */
   	int *free_ptr;               /* adjacency ptr to free copy of graph */
   	double eigtol;               /* tolerance in eigenvector calculation */
	int ndims;                   /* dimension of recursive partitioning */
   	int spec;                    /* run spectral global decomposition? */
   	int inert;                   /* run inertial global decomposition? */
   	int rcb;                     /* run coordinate global decomposition? */
   	int KL;                      /* run Kernighan-Lin local optimization? */
   	long seed;                   /* for random graph mutations */
   	int mediantype;              /* method for partitioning eigenvector */
   	int mkconnected;             /* check connectivity & add phantom edges? */
   	int solver_flag;             /* which eigensolver should I use? */
   	int coarse_flag;             /* should I use multilevel techniques? */
   	int vmax;                    /* if so, how many vertices to coarsen down to? */
   	int igeom;                   /* geometry dimension if inertial method */
   	char graphname[NAME_LENGTH]; /* name of graph input file */
   	char geomname[NAME_LENGTH];  /* name of geometry input file */
   	char assignname[NAME_LENGTH];/* name of output file */
	char outname[NAME_LENGTH];   /* name of file for outputting run results */ 
   	char *geomptr, *assignptr;   /* names or null pointers */
   	int scatt;                   /* run scatt global partitioner? */
   	int randm;                   /* run randm global partitioner? */
   	int lin;                     /* run linear global partitioner? */
   	double time;                 /* timing marker */
	
	/* variables of Chaco's interface() */
    extern int ECHO;                /* print timing output to file? */
    extern int OUTPUT_TIME; /* following for output options */
   	double total_time;		/* total time spent in partitioning code */
   	struct vtx_data **graph;	/* graph data structure */
   	float **coords;		/* coordinates for vertices if used */
   	int *vptr;			/* loops through vwgts */
   	double vwgt_sum;		/* sum of all the vertex weights */
   	int flag;			/* return code from balance */
   	int nedges;			/* number of edges in graph */
   	int using_vwgts;		/* are vertex weights being used? */
   	int using_ewgts;     /* are edge weights being used? */
   	extern double start_time;		/* time routine is entered */
   	int i;			/* loop counter */
   	double time_tol;		/* effective resolution of timer */
   	double seconds();
   	int balance();
   	void reformat(), free_graph();
	void do_timing();
	
	/* variables of Chaco's input_graph() */
	int vtxwts;                  /* are vertex weights in input file? */
	int edgewts;                 /* are edge weights in input file? */
	int eweight;                 /* edgeweight being read */
	int narcs;                   /* number of edges in graph */
	int size;					 /* length of all edge lists */

	/* variables of Chaco's reformat() */
	struct vtx_data *links;      /* ptr to space for data for all vtx */
	int *edges;                  /* space for all adjacency lists */
	float *eweights;             /* space for all edge weights */
	float sum;                   /* sum of edge weights for a vtx */
	int j;                       /* loop counter */

	/* additionally declared variables */
	CLUSTER *cptr;               /* loops through clusters */
	int *pvtx;                   /* loops through vtxs */
	int *loc2glob;               /* map from local to global numbering */
	int total_edges;
	int count;

	start_time = seconds();

	/* begin part of Chaco's input_queries() */

	/* Is an outputfile wanted */
	fin = NULL;
	fingeom = NULL;
	if (ECHO)
	{
		strcpy(outname,"output.lb4");
		outfile = fopen(outname, "w");
	}
	else outfile = NULL;
	
	/* Initialize then set the flags */
	scatt = randm = lin = inert = spec = KL = FALSE;
	coarse_flag = solver_flag = 0;

	/* which partitioning strategy to use */
	if (glob_method == 6) 
	{ 
		scatt = TRUE;
	}
	else if (glob_method == 5) 
	{ 
		randm = TRUE;
	}
	else if (glob_method == 4) 
	{ 
		lin = TRUE;
	}
	else if (glob_method == 3) 
	{ 
		inert = TRUE;
		fingeom = NULL;
	}
	else if (glob_method == 0)
	{
		rcb = TRUE;
		fingeom = NULL;
	}
	else if (glob_method == 2) 
	{ 
		spec = TRUE;
		fingeom = NULL;
		if (eigen >= 5)
		{
			solver_flag =	eigen - 4;
			coarse_flag = 1;
			vmax = coarse;
		}
		else solver_flag = eigen;
	}
	else if (glob_method == 1) 
	{
		spec = TRUE;
		KL = TRUE;
		fingeom = NULL;
		solver_flag = eigen;
		coarse_flag = 2;
		vmax = coarse;
	}

	/* Get local method, if any */
	if (loc || glob_method == 1)
   	{
			KL = TRUE;
   	}

   	/* Get number of dimensions in which to partition at each level. For inertial
   	   must do bisection hence ndims = 1. Also restricted by nsets_tot. */
   	nsets_glob = nsets_tot;

	if (nsets_tot < 4) {ndims = 1;}
	else if (nsets_tot < 8) 
	{
			if (dims <= 2) 
			{
				ndims = dims; 
			}
			else ndims = 2;
	}
	else if (nsets_tot >= 8) 
	{
		ndims = dims;
	}

	/* end part of Chaco's input_queries() */

	/* analyse the mode parameter */
	if (mode)
	{
		CHECK_INPUT    = GET_CHECK_INPUT(mode);
		OUTPUT_METRICS = GET_OUTPUT_METRICS(mode);
		OUTPUT_ASSIGN  = GET_OUTPUT_ASSIGN(mode);
		DEBUG_MEMORY   = GET_DEBUG_MEMORY(mode);
		OUTPUT_TIME    = GET_OUTPUT_TIME(mode);
	}		

	/* begin part of Chaco's main() */

	if (OUTPUT_ASSIGN > 0) 
	{
		strcpy(assignname,"assign.lb4"); 
		assignptr = assignname; 
	} 
	else 
		assignptr = NULL;

	if (inert || rcb)
	{

		/* begin part of Chaco's input_geom */
		#ifdef __TWODIM__
		igeom = 2;
		#endif
		#ifdef __THREEDIM__
		igeom = 3;
		#endif
		x = (float *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(float));
		if (!MEM_OK) return;
		y = (float *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(float));
		if (!MEM_OK) return;
		#ifdef __THREEDIM__
		z = (float *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(float));
		if (!MEM_OK) return;
		#endif

		for (i=1; i<=nvtxs; i++)
		{ 
			cptr = clusters[i-1];
			x[i] = (float) cptr->sx;
			y[i] = (float) cptr->sy;
			#ifdef __THREEDIM__
			z[i] = (float) cptr->sz;
			#endif
		} 

		/* end part of Chaco's input_geom */

	}
	else
	{
		geomptr = NULL;
		igeom   = 0;
		fingeom = NULL;
		x = y = z = NULL;
	}

	eigtol = EIGEN_TOLERANCE;
	mediantype = MAPPING_TYPE;
	seed = RANDOM_SEED;
	assignment = assign;

	/* Initialize the connectivity flag. */
	if (spec == TRUE) mkconnected = MAKE_CONNECTED;
	else mkconnected = FALSE;

	/* end part of Chaco's main */


	/* begin part of Chaco's interface() */

	/* begin part of Chaco's reformat() */

	if (weights == 0)
	{
		using_vwgts = 0;
		using_ewgts = 0;
	}
	else if (weights == 1)
	{
		using_vwgts = 1;
		using_ewgts = 0;
	}
	else if (weights == 2)
	{
		using_vwgts = 0;
		using_ewgts = 1;
    }
    else if (weights == 3)
    {
		using_vwgts = 1;
		using_ewgts = 1;
	}        

	graph = (struct vtx_data **) (MEM_OK = smalloc((unsigned)(nvtxs+1)
	                                     *sizeof(struct vtx_data *));
	if (!MEM_OK) return;

	/* Set up all the basic data structure for the vertices. */
	/* Replace many small mallocs by a few large ones. */
	links = (struct vtx_data *) (MEM_OK = smalloc((unsigned)(nvtxs)*sizeof(struct vtx_data));
	if (!MEM_OK) return;
	for (i=1; i<=nvtxs; i++) 
	{
		graph[i] = links++;
	}

	/* count number of edges self-edges included */
	size = 0;
	for (j=0; j<nvtxs; j++)
	{
		cptr = clusters[j];
		i = 0;
		while (cptr->edges[i++]!=0)
		{
			size++;
		}
	}
	total_edges = size;

	edges = (int *) (MEM_OK = smalloc((unsigned)size*sizeof(int));
	if (!MEM_OK) return;
	if (using_ewgts) 
	{
		eweights = (float *) (MEM_OK = smalloc((unsigned)size*sizeof(float));
		if (!MEM_OK) return;
	}
	else eweights = NULL;

	/* Now fill in all the data fields. */
	for (i=1; i<=nvtxs; i++) {
		cptr = clusters[i-1];
		for (j=0,size=0; cptr->edges[j]!=0 && j<MAX_SIDES_OF_ELEM+1; j++) 
				size++;
		if (using_vwgts) 
			graph[i]->vwgt = cptr->level_size[cptr->depth];
		else 
			graph[i]->vwgt = 1;
		graph[i]->nedges = size;
		graph[i]->edges = edges;
		*edges = cptr->edges[0];
		edges++;
		count++;
		for (j=size-1; j>0; j--) 
		{
			*edges = cptr->edges[j];
			edges++;
			count++;
		}
		if (using_ewgts) {
			graph[i]->ewgts = eweights;
			eweights++;
			sum = 0;
			for (j=size-1; j; j--) 
			{
				sum += cptr->edges[j];
				*eweights++ = cptr->edges[j];
			}
			graph[i]->ewgts[0] = -sum;
		}
		else graph[i]->ewgts = NULL;
	}

	/* end part of Chaco's reformat() */

	if (x == NULL) coords = NULL;
	else {
		coords = (float **) (MEM_OK = smalloc((unsigned) 3*sizeof(float *));
		if (!MEM_OK) return;
		coords[0] = x;
		coords[1] = y;
		coords[2] = z;
	}

	/* allocate map from local back to global numbering */
	loc2glob = (int *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(int));
	if (!MEM_OK) return;

	reformat_graph(clusters,graph,nvtxs,&nedges,glob2loc,loc2glob,using_ewgts);

	flag = balance(graph, nvtxs, nedges, ndims, nsets_tot, using_vwgts, spec, 
				inert, rcb, KL, mediantype, mkconnected, solver_flag, coarse_flag,
				vmax, eigtol, seed, igeom, coords, assignment, goal, scatt, 
				randm, lin, graphname, geomname, assignname, outfile,dimx,dimy); 
	if (!MEM_OK) return;

	/* reset glob2loc map */ 
	for (i=1; i<=nvtxs; i++)
	{
		j = loc2glob[i];
		glob2loc[j] = 0;
	}
	sfree((char *) loc2glob);

	if (flag){char buf[150]; sprintf(buf," balance returned with flag %d\n", flag);UserWrite(buf);}

	/* free allocated space */
	if (x != NULL) sfree((char*) x);
	if (y != NULL) sfree((char*) y);
	/* free eweights, edges, links and graph */
	free_graph(graph);
	if (coords != NULL) sfree((char *) coords);

	/* now do timing */
	do_timing(outfile);

	if (DEBUG_MEMORY > 0) {
{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
	smalloc_stats();
	}

{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
	return(0);
}


/* timing for the different parts during the partitioning */

void do_timing(FILE *outfile)
{
   extern int ECHO;                /* print timing output to file? */
   extern int OUTPUT_TIME; /* following for output options */

   extern double input_time;
   extern double partition_time;
   extern double reformat_time;
   extern double check_input_time;
   extern double count_time;
   extern double print_assign_time;
   
   extern double inertial_time;
   extern double rcb_time;
   extern double inertial_axis_time;
   extern double median_time;

   extern double coarsen_time;     /* following for timing coarsening */
								   /* algorithm */
   extern double match_time;
   extern double make_cgraph_time;
   extern double make_cewgts_time;
   extern double adj_cewgts_time;

   extern double lanczos_time;     /* following for timing Lanczos */
								   /* algorithms */
   extern double splarax_time;
   extern double orthog_time;
   extern double ql_time;
   extern double tevec_time;
   extern double ritz_time;
   extern double evec_time;
   extern double check_time;
   extern double blas_time;
   extern double init_time;
   extern double scan_time;
   extern double debug_time;
   extern double probe_time;
   extern double pause_time;
   double other_time;
   
   extern double rqi_symmlq_time;/* following for timing ML RQI/Symmlq */
								 /* algorithm */
   extern double refine_time;
   
   extern double kl_total_time;    /* following for timing KL algorithm */
   extern double kl_init_time;
   extern double nway_kl_time;
   extern double start_time;		/* time routine is entered */

   double time_tol;		/* effective resolution of timer */
   double total_time;           /* total time spent in partitioning code*/
   int print2file;		/* should I print to a file? */
   double seconds();

   time_tol = 0.005;

   print2file = (ECHO < 0);
   if (OUTPUT_TIME > 0) {
      total_time = seconds() - start_time + input_time;
      if (total_time != 0) {
        {char buf[150]; sprintf(buf,"\nTotal time: %g sec.\n", total_time);UserWrite(buf);}
         if (print2file) fprintf(outfile, "\nTotal time: %g sec.\n", total_time);
         if (input_time != 0) {
	   {char buf[150]; sprintf(buf,"  file input %g\n", input_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  file input %g\n", input_time);
	 }
         if (reformat_time != 0){
	   {char buf[150]; sprintf(buf,"  reformatting %g\n", reformat_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  reformatting %g\n", reformat_time);
	 }
         if (check_input_time != 0){
	   {char buf[150]; sprintf(buf,"  checking input %g\n", check_input_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  checking input %g\n", check_input_time);
	 }
         if (partition_time != 0){
	   {char buf[150]; sprintf(buf,"  partitioning %g\n", partition_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  partitioning %g\n", partition_time);
	 }
         if (count_time != 0){
	   {char buf[150]; sprintf(buf,"  evaluation %g\n", count_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  evaluation %g\n", count_time);
	 }
         if (print_assign_time != 0){
	   {char buf[150]; sprintf(buf,"  printing assignment file %g\n", print_assign_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  printing assignment file %g\n", print_assign_time);
	 }
         other_time = total_time - input_time - reformat_time 
            - check_input_time - partition_time - count_time - print_assign_time;
         if (other_time > time_tol){
	   {char buf[150]; sprintf(buf,"  other %g\n", other_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  other %g\n", other_time);
	 }
      }
   }

   if (OUTPUT_TIME > 1) {
      if (inertial_time != 0) {
         if (inertial_time != 0) {
	   {char buf[150]; sprintf(buf,"\nInertial time: %g sec.\n", inertial_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "\nInertial time: %g sec.\n", inertial_time);
	 }
         if (inertial_axis_time != 0) {
	   {char buf[150]; sprintf(buf,"  inerial axis %g\n", inertial_axis_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  inerial axis %g\n", inertial_axis_time);
	 }
         if (median_time != 0) {
	   {char buf[150]; sprintf(buf,"  median finding %g\n", median_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  median finding %g\n", median_time);
	 }
         other_time = inertial_time - inertial_axis_time - median_time;
         if (other_time > time_tol) {
	   {char buf[150]; sprintf(buf,"  other %g\n", other_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  other %g\n", other_time);
	 }
      }
      if (rcb_time != 0) {
         if (rcb_time != 0) {
	   {char buf[150]; sprintf(buf,"\nRcb time: %g sec.\n", rcb_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "\nRcb time: %g sec.\n", rcb_time);
	 }
	  }

      if (kl_total_time != 0) {
         if (kl_total_time != 0) {
	   {char buf[150]; sprintf(buf,"\nKL time: %g sec.\n", kl_total_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "\nKL time: %g sec.\n", kl_total_time);
	 }
         if (kl_init_time != 0) {
	   {char buf[150]; sprintf(buf,"  initialization %g\n", kl_init_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  initialization %g\n", kl_init_time);
	 }
         if (nway_kl_time != 0) {
	   {char buf[150]; sprintf(buf,"  nway refinement %g\n", nway_kl_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  nway refinement %g\n", nway_kl_time);
	 }
         other_time = kl_total_time - kl_init_time - nway_kl_time;
         if (other_time > time_tol) {
	   {char buf[150]; sprintf(buf,"  other %g\n", other_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  other %g\n", other_time);
	 }
      }

      if (coarsen_time != 0) {
        {char buf[150]; sprintf(buf,"\nCoarsening time: %g sec.\n", coarsen_time);UserWrite(buf);}
         if (print2file) fprintf(outfile, "\nCoarsening time: %g sec.\n", coarsen_time);
         if (match_time != 0) {
	   {char buf[150]; sprintf(buf,"  maxmatch %g\n", match_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  maxmatch %g\n", match_time);
	 }
         if (make_cgraph_time != 0) {
	   {char buf[150]; sprintf(buf,"  makecgraph %g\n", make_cgraph_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  makecgraph %g\n", make_cgraph_time);
	 }
         if (make_cewgts_time != 0) {
	   {char buf[150]; sprintf(buf,"    makecewgts %g\n", make_cewgts_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "    makecewgts %g\n", make_cewgts_time);
	 }
         if (adj_cewgts_time != 0) {
	   {char buf[150]; sprintf(buf,"      adjcewgts %g\n", adj_cewgts_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "      adjcewgts %g\n", adj_cewgts_time);
	 }
      }

      if (lanczos_time != 0) {
        {char buf[150]; sprintf(buf,"\nLanczos time: %g sec.\n", lanczos_time);UserWrite(buf);}
         if (print2file) fprintf(outfile, "\nLanczos time: %g sec.\n", lanczos_time);
         if (splarax_time != 0) {
	   {char buf[150]; sprintf(buf,"  matvec/solve %g\n", splarax_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  matvec/solve %g\n", splarax_time);
	 }
         if (blas_time != 0) {
	   {char buf[150]; sprintf(buf,"  vector ops %g\n", blas_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  vector ops %g\n", blas_time);
	 }
         if (evec_time != 0) {
	   {char buf[150]; sprintf(buf,"  assemble evec %g\n", evec_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  assemble evec %g\n", evec_time);
	 }
         if (init_time != 0) {
	   {char buf[150]; sprintf(buf,"  malloc/init/free %g\n", init_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  malloc/init/free %g\n", init_time);
	 }
         if (check_time != 0) {
	   {char buf[150]; sprintf(buf,"  check results %g\n", check_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  check results %g\n", check_time);
	 }
         if (orthog_time != 0) {
	   {char buf[150]; sprintf(buf,"  maintain orthog %g\n", orthog_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  maintain orthog %g\n", orthog_time);
	 }
         if (scan_time != 0) {
	   {char buf[150]; sprintf(buf,"  scan %g\n", scan_time);UserWrite(buf);} 
	    if (print2file) fprintf(outfile, "  scan %g\n", scan_time); 
	 }
         if (debug_time != 0) {
	   {char buf[150]; sprintf(buf,"  debug/warning/check %g\n", debug_time);UserWrite(buf);} 
	    if (print2file) fprintf(outfile, "  debug/warning/check %g\n", debug_time); 
	 }
         if (ql_time != 0) {
	   {char buf[150]; sprintf(buf,"  ql/bisection %g\n", ql_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  ql/bisection %g\n", ql_time);
	 }
         if (tevec_time != 0) {
	   {char buf[150]; sprintf(buf,"  tevec %g\n", tevec_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  tevec %g\n", tevec_time);
	 }
         if (ritz_time != 0) {
	   {char buf[150]; sprintf(buf,"  compute ritz %g\n", ritz_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  compute ritz %g\n", ritz_time);
	 }
         if (pause_time != 0) {
	   {char buf[150]; sprintf(buf,"  pause %g\n", pause_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  pause %g\n", pause_time);
	 }
         if (probe_time != 0) {
	   {char buf[150]; sprintf(buf,"\n  probe %g\n", probe_time);UserWrite(buf);} 
	    if (print2file) fprintf(outfile, "\n  probe %g\n", probe_time); 
	 }
         other_time = lanczos_time - splarax_time - orthog_time 
   	   - ql_time - tevec_time - ritz_time - evec_time - check_time 
   	   - blas_time - init_time - scan_time - debug_time - pause_time;
         if (other_time > time_tol && other_time != lanczos_time) {
	   {char buf[150]; sprintf(buf,"  other %g\n", other_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  other %g\n", other_time);
	 }
      } 

      if (rqi_symmlq_time != 0) {
        {char buf[150]; sprintf(buf,"\nRQI/Symmlq time: %g sec.\n", rqi_symmlq_time);UserWrite(buf);}
         if (print2file) fprintf(outfile, "\nRQI/Symmlq time: %g sec.\n", rqi_symmlq_time);
         if (coarsen_time != 0) {
	   {char buf[150]; sprintf(buf,"  coarsening %g\n", coarsen_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  coarsening %g\n", coarsen_time);
	 }
         if (refine_time != 0) {
	   {char buf[150]; sprintf(buf,"  refinement %g\n", refine_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  refinement %g\n", refine_time);
	 }
         if (lanczos_time != 0) {
	   {char buf[150]; sprintf(buf,"  lanczos %g\n", lanczos_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  lanczos %g\n", lanczos_time);
	 }
         other_time = rqi_symmlq_time - coarsen_time - refine_time - lanczos_time; 
         if (other_time > time_tol) {
	   {char buf[150]; sprintf(buf,"  other %g\n", other_time);UserWrite(buf);}
	    if (print2file) fprintf(outfile, "  other %g\n", other_time);
	 }
      }
   }
}

