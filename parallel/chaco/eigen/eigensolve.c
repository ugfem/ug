// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"

void eigensolve(graph, nvtxs, nedges, maxdeg, yvecs, evals, ndims, eigtol,
                vwsqrt, solver_flag, coarse_flag, vmax, igeom, coords, using_ewgts,
                using_vwgts)
/* Invoke the eigenvector calculation */
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vertices in graph */
int nedges;                     /* number of edges in graph */
double maxdeg;                  /* largest (weighted) degree of a vertex */
double **yvecs;                 /* space for pointing to eigenvectors */
double *evals;                  /* eigenvalues associated with eigenvectors */
int ndims;                      /* number of eigenvectors (2^d sets) */
double eigtol;                  /* tolerance on eigenvectors */
double *vwsqrt;                 /* sqrt of vertex weights (length nvtxs+1) */
int solver_flag;                /* flag indicating which solver to use */
int coarse_flag;                /* use multi-level techniques? */
int vmax;                       /* if so, how many vtxs to coarsen down to? */
int igeom;                      /* geometric dimensionality if given coords */
float **coords;                 /* coordinates of vertices */
int using_ewgts;                /* are edge weights being used? */
int using_vwgts;                /* are vertex weights being used? */
{
  extern Heap   *heap;      /* pointer to heap of multigrid */
  extern double *MEM_OK;    /* variable for memory overflow exeception */
  extern int DEBUG_TRACE;       /* trace the execution of the code */
  extern int DEBUG_EVECS;       /* debug flag for eigenvector generation */
  extern int DEBUG_PERTURB;     /* debug flag for matrix perturbation */
  extern int LANCZOS_SO_TIME;   /* perform detailed timing on Lanczos_SO? */
  extern int PERTURB;           /* randomly perturb to break symmetry? */
  extern int NPERTURB;          /* number of edges to perturb */
  extern double PERTURB_MAX;    /* maximum size of perturbation */
  extern int COARSE_NLEVEL_RQI; /* how often to invoke RQI while uncoarsening */
  extern double lanczos_time;   /* time spent in Lanczos algorithm */
  extern double rqi_symmlq_time; /* time spent in Multilevel RQI/Symmlq method */
  double bound[MAXDIMS+1];      /* ritz approximation bounds to eigenpairs */
  double time;                  /* time marker */
  int step;                     /* current step in RQI counting */
  int nstep;                    /* number of uncoarsening levels between RQIs */
  int version;                  /* which version of sel. orth. to use */
  double seconds();
  void coarsen(), lanczos_FO(), lanczos_SO(), lanczos_SO_no_time(), vecout();
  void perturb_init(), perturb_clear(), lanczos_IO();

  if (DEBUG_TRACE > 0) {
    {char buf[150]; sprintf(buf,"Entering eigensolve, nvtxs = %d, nedges = %d\n", nvtxs, nedges);UserWrite(buf);}
  }

  if (coarse_flag == 1) {       /* Solve using multi-level scheme SYMMLQ/RQI. */
    time = seconds();
    nstep = COARSE_NLEVEL_RQI;
    step = 0;
    coarsen(graph, nvtxs, nedges, yvecs, ndims, igeom, coords, nstep, step,
            vmax, eigtol, using_ewgts, using_vwgts, solver_flag, FALSE);
    if (!MEM_OK) return;
    rqi_symmlq_time += seconds() - time;
  }
  else {                        /* Compute eigenvectors directly. */
    if (PERTURB) {
      if (NPERTURB > 0 && PERTURB_MAX > 0.0) {
        perturb_init(nvtxs);
        if (DEBUG_PERTURB > 0) {
          {char buf[150]; sprintf(buf,"Matrix being perturbed with scale %e\n", PERTURB_MAX);UserWrite(buf);}
        }
      }
      else if (DEBUG_PERTURB > 0) {
        {char buf[150]; sprintf(buf,"Matrix not being perturbed\n");UserWrite(buf);}
      }
    }

    if (solver_flag == 1) {
      time = seconds();
      lanczos_FO(graph,nvtxs,ndims,yvecs,evals,bound,eigtol,vwsqrt);
      if (!MEM_OK) return;
      lanczos_time += seconds() - time;
    }
    if (solver_flag == 2) {
      time = seconds();
      lanczos_IO(graph,nvtxs,ndims,yvecs,evals,bound,eigtol,vwsqrt);
      if (!MEM_OK) return;
      lanczos_time += seconds() - time;
    }
    else if (solver_flag == 3) {
      version = 1;
      time = seconds();
      if (LANCZOS_SO_TIME) {
        lanczos_SO(graph,nvtxs,ndims,yvecs,evals,bound,eigtol,vwsqrt,maxdeg,version);
        if (!MEM_OK) return;
      }
      else {
        lanczos_SO_no_time(graph,nvtxs,ndims,yvecs,evals,bound,eigtol,vwsqrt,maxdeg,version);
        if (!MEM_OK) return;
      }
      lanczos_time += seconds() - time;
    }
    else if (solver_flag == 4) {
      version = 2;
      time = seconds();
      if (LANCZOS_SO_TIME) {
        lanczos_SO(graph,nvtxs,ndims,yvecs,evals,bound,eigtol,vwsqrt,maxdeg,version);
        if (!MEM_OK) return;
      }
      else {
        lanczos_SO_no_time(graph,nvtxs,ndims,yvecs,evals,bound,eigtol,vwsqrt,maxdeg,version);
        if (!MEM_OK) return;
      }
      lanczos_time += seconds() - time;
    }

    if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0) {
      perturb_clear();
    }
  }
  if (DEBUG_EVECS > 4) {
    for (nstep=1; nstep<=ndims; nstep++) {
      vecout(yvecs[nstep], 1, nvtxs, "Eigenvector");
    }
  }
}
