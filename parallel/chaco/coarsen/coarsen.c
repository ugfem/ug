// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"

void coarsen(graph, nvtxs, nedges, yvecs, ndims, igeom, coords, nstep, step,
	vmax, eigtol, using_ewgts, using_vwgts, solver_flag, give_up)
/* Coarsen until nvtxs < vmax, compute and uncoarsen. */
struct vtx_data **graph;	/* array of vtx data for graph */
int nvtxs;			/* number of vertices in graph */
int nedges;			/* number of edges in graph */
double **yvecs;			/* eigenvectors returned */
int ndims;			/* number of eigenvectors to calculate */
int igeom;			/* dimension for geometric information */
float **coords;			/* coordinates for vertices */
int nstep;			/* number of coarsenings between RQI steps */
int step;			/* current step number */
int vmax;			/* largest subgraph to stop coarsening */
double eigtol;			/* tolerence in eigen calculation */
int using_ewgts;		/* are edge weights being used? */
int using_vwgts;		/* are vertices weights being used? */
int solver_flag;		/* which eigensolver to use */
int give_up;			/* has coarsening bogged down? */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int DEBUG_COARSEN;	/* debug flag for coarsening */
   extern int PERTURB;		/* was matrix perturbed in Lanczos? */
   extern float COARSEN_RATIO_MIN;	/* min vtx reduction for coarsening */
   extern double refine_time;	/* time for RQI/Symmlq iterative refinement */
   struct vtx_data **cgraph;	/* array of vtx data for coarsened graph */
   struct ipairs *merged;	/* vtx pairs that get merged */
   struct fpairs *reduction;	/* multipliers for reducing edge weights */
   struct orthlink *orthlist;	/* list of lower evecs to suppress */
   struct orthlink *newlink;	/* lower evec to suppress */
   double *cyvecs[MAXDIMS+1];	/* eigenvectors for subgraph */
   double evals[MAXDIMS+1];	/* eigenvalues returned */
   double *r1, *r2, *work;	/* space needed by symmlq/RQI */
   double *v, *w, *x, *y;	/* space needed by symmlq/RQI */
   double evalest;		/* eigenvalue estimate returned by RQI */
   double maxdeg;		/* maximum weighted degree of a vertex */
   float **ccoords;		/* coordinates for coarsened graph */
   double *vwsqrt = NULL;	/* square root of vertex weights */
   double norm, alpha;		/* values used for orthogonalization */
   double initshift;		/* initial shift for RQI */
   int *mflag, *v2cv;		/* unused parameters */
   int oldperturb;		/* saves PERTURB value */
   int cnvtxs;			/* number of vertices in coarsened graph */
   int cnedges;			/* number of edges in coarsened graph */
   int nmerged;			/* number of vertices that get merged */
   int nextstep;		/* next step in RQI test */
   int i, j;			/* loop counters */
   double time;			/* time marker */
   
   double dot(), normalize(), find_maxdeg(), seconds();
   struct orthlink *makeorthlnk();
   void makevwsqrt(), eigensolve(), coarsen1(), interpolate(), orthogvec();
   void orthog1(), RQI(), scadd(), free_graph();

   if (DEBUG_COARSEN > 0) {
     {char buf[150]; sprintf(buf," Entering coarsen, step=%d, nvtxs=%d, nedges=%d, vmax=%d\n", step, nvtxs, nedges, vmax);UserWrite(buf);}
   }

   /* Is problem small enough to solve? */
   if (nvtxs <= vmax || give_up) {
      if (using_vwgts) {
         vwsqrt = (double *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(double));
         if (!MEM_OK) return;
	 makevwsqrt(vwsqrt, graph, nvtxs);
      }
      else vwsqrt = NULL;
      maxdeg = find_maxdeg(graph, nvtxs);

      eigensolve(graph, nvtxs, nedges, maxdeg, yvecs, evals, ndims, eigtol,
		 vwsqrt, solver_flag, FALSE, 0, igeom, coords, using_ewgts,
		using_vwgts);
      if (!MEM_OK) return;

      if (vwsqrt != NULL) sfree((char *) vwsqrt);
      return;
   }

   /* Otherwise I have to coarsen. */
   if (coords == NULL || igeom <= 0) ccoords = NULL;
   else {ccoords = (float **) (MEM_OK = smalloc((unsigned) (igeom+1)*sizeof(float *));
         if (!MEM_OK) return;
        }
   coarsen1(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &merged, 
	   &reduction, &mflag, &v2cv, igeom, coords, ccoords, using_ewgts);
   if (!MEM_OK) return;
   sfree((char *) mflag);
   sfree((char *) v2cv);

   /* If coarsening isn't working very well, give up and partition. */
   give_up = FALSE;
   if (nvtxs*COARSEN_RATIO_MIN < cnvtxs) {
     {char buf[150]; sprintf(buf,"WARNING: Coarsening not making enough progress, nvtxs = %d, cnvtxs = %d.\n", nvtxs, cnvtxs);UserWrite(buf);}
     {char buf[150]; sprintf(buf,"         Recursive coarsening being stopped prematurely.\n");UserWrite(buf);}
      give_up = TRUE;
   }

   /* Create space for subgraph yvecs. */
   for (i=1; i<=ndims; i++) {
      cyvecs[i] = (double *) (MEM_OK = smalloc((unsigned) (cnvtxs+1)*sizeof(double));
      if (!MEM_OK) return;
   }

   /* Now recurse on coarse subgraph. */
   nextstep = step - 1;
   if (step <= 0) nextstep = nstep - 1; 
   coarsen(cgraph, cnvtxs, cnedges, cyvecs, ndims, igeom, ccoords, nstep,
	   nextstep, vmax, eigtol, TRUE, TRUE, solver_flag, give_up);
   if (!MEM_OK) return;

   nmerged = nvtxs - cnvtxs;
   interpolate(yvecs, cyvecs, ndims, graph, nvtxs, cgraph, cnvtxs,
	       merged, reduction, nmerged, using_ewgts);
   if (!MEM_OK) return;

   /* I need to do Rayleigh Quotient Iteration each nstep stages. */
   time = seconds();
   if (step == 0 || step == nstep) {
      oldperturb = PERTURB;
      PERTURB = FALSE;
      /* Should I do some orthogonalization here against vwsqrt? */
      if (using_vwgts) {
         vwsqrt = (double *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(double));
         if (!MEM_OK) return;
	 makevwsqrt(vwsqrt, graph, nvtxs);

         for (i=1; i<=ndims; i++) orthogvec(yvecs[i], 1, nvtxs, vwsqrt);
      }
      else for (i=1; i<=ndims; i++) orthog1(yvecs[i], 1, nvtxs);

      /* Allocate space that will be needed in RQI. */
      r1 = (double *) (MEM_OK = smalloc((unsigned) 7*(nvtxs+1)*sizeof(double));
      if (!MEM_OK) return;
      r2 = &r1[nvtxs+1];
      v = &r1[2*(nvtxs+1)];
      w = &r1[3*(nvtxs+1)];
      x = &r1[4*(nvtxs+1)];
      y = &r1[5*(nvtxs+1)];
      work = &r1[6*(nvtxs+1)];

      initshift = 0;
      orthlist = NULL;
      for (i=1; i<ndims; i++) {
	 normalize(yvecs[i], 1, nvtxs);
         RQI(graph, yvecs[i], nvtxs, r1, r2, v, w, x, y, work, 
	     eigtol, initshift, &evalest, vwsqrt, orthlist);
         if (!MEM_OK) return;

         /* Now orthogonalize higher yvecs against this one. */
         norm = dot(yvecs[i], 1, nvtxs, yvecs[i]);
         for (j=i+1; j<=ndims; j++) {
	    alpha = -dot(yvecs[j], 1, nvtxs, yvecs[i])/norm;
	    scadd(yvecs[j], 1, nvtxs, alpha, yvecs[i]);
	 }

	 /* Now prepare for next pass through loop. */
	 initshift = evalest;
	 newlink = makeorthlnk();
     if (!MEM_OK) return;
	 newlink->vec = yvecs[i];
	 newlink->pntr = orthlist;
	 orthlist = newlink;

      }
      normalize(yvecs[ndims], 1, nvtxs);

      RQI(graph, yvecs[ndims], nvtxs, r1, r2, v, w, x, y, work,
	  eigtol, initshift, &evalest, vwsqrt, orthlist);
      if (!MEM_OK) return;
      refine_time += seconds() - time;

      /* Free the space allocated for RQI. */
      while (orthlist != NULL) {
	 newlink = orthlist->pntr;
	 sfree((char *) orthlist);
	 orthlist = newlink;
      }
      sfree((char *) r1);
      if (vwsqrt != NULL)  sfree((char *) vwsqrt);
      PERTURB = oldperturb;
   }
   if (DEBUG_COARSEN > 0) {
     {char buf[150]; sprintf(buf," Leaving coarsen, step=%d\n", step);UserWrite(buf);}
   }

   /* Free the space that was allocated. */
   if (ccoords != NULL) {
      for (i=0; i<igeom; i++) sfree((char *) ccoords[i]);
      sfree((char *) ccoords);
   }
   for (i=ndims; i>0; i--) sfree((char *) cyvecs[i]);
   sfree((char *) reduction);
   sfree((char *) merged);
   free_graph(cgraph);
}
