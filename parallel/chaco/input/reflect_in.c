// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"

/* Print out the input options. */
void reflect_input(nvtxs, spec, inert, KL, solver_flag, coarse_flag, vmax,
                   ndims, ndims_tot, igeom, scatt, randm, lin, graphname, geomname,
                   assignname, outfile)
int nvtxs;                      /* number of vertices */
int spec;                       /* invoke spectral method? */
int inert;                      /* invoke inertial method? */
int KL;                         /* invoke Kernighan-Lin? */
int solver_flag;                /* type of eigensolver */
int coarse_flag;                /* flag for recursive coarsening  */
int vmax;                       /* smallest acceptable coarsened nvtxs */
int ndims;                      /* partitioning level */
int ndims_tot;                  /* total number of cuts to make */
int igeom;                      /* geometric dimension for inertial method */
int scatt;                      /* use scattered global partitioner? */
int randm;                      /* use randm global partitioner? */
int lin;                        /* use linear global partitioner? */
char *graphname, *geomname;     /* input file names */
char *assignname;               /* assignment output file name */
FILE *outfile;                  /* file to write output to */
{
  extern double EIGEN_TOLERANCE;        /* eigen-tolerance convergence criteria */
  extern double BISECTION_SAFETY;       /* divides Lanczos bisection tolerance  */
  extern int MAKE_CONNECTED;    /* connect graph if using spectral method? */
  extern int LANCZOS_SO_INTERVAL;       /* interval between orthog checks in SO */
  extern int COARSE_NLEVEL_RQI; /* # levels between RQI calls in uncoarsening */
  extern int KL_RANDOM;         /* use randomness in Kernighan-Lin? */
  extern int KL_METRIC;         /* KL interset cost: 1=>cuts, 2=>hops */
  extern int KL_NTRIES_BAD;     /* # unhelpful passes before quitting KL */
  extern int KL_BAD_MOVES;      /* number of unhelpful moves in a row allowed */
  extern int COARSE_NLEVEL_KL;  /* # levels between KL calls in uncoarsening */
  extern int PERTURB ;          /* randomly perturb matrix in spectral method? */
  extern int NPERTURB;          /* if so, how many edges to modify? */
  extern double PERTURB_MAX;    /* largest value for perturbation */
  extern int ECHO;              /* copy input parameters back to screen? */
  int print2file;               /* should I print to a file as well? */
  char *true_or_false();

  print2file = (ECHO < 0);

  {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
  if (graphname != NULL) {
    {char buf[150]; sprintf(buf,"Graph file: %s; ", graphname);UserWrite(buf);}
    if (print2file) fprintf(outfile, "Graph file: %s; ", graphname);
  }
  {char buf[150]; sprintf(buf,"Number of vertices = %d\n", nvtxs);UserWrite(buf);}
  if (print2file) fprintf(outfile, "Number of vertices = %d\n", nvtxs);
  if (assignname != NULL) {
    {char buf[150]; sprintf(buf,"Assignment output file: %s\n", assignname);UserWrite(buf);}
    if (print2file) fprintf(outfile, "Assignment output file: %s\n", assignname);
  }

  /* Print global partitioning strategy. */
  {char buf[150]; sprintf(buf,"Global method: ");UserWrite(buf);}
  if (print2file) fprintf(outfile, "Global method: ");
  if (coarse_flag == 2) {
    {char buf[150]; sprintf(buf,"Multilevel\n");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Multilevel\n");
  }
  else if (spec) {
    {char buf[150]; sprintf(buf,"Spectral\n");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Spectral\n");
  }
  else if (inert) {
    {char buf[150]; sprintf(buf,"Inertial\n");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Inertial\n");
  }
  else if (lin) {
    {char buf[150]; sprintf(buf,"Linear\n");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Linear\n");
  }
  else if (randm) {
    {char buf[150]; sprintf(buf,"Random\n");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Random\n");
  }
  else if (scatt) {
    {char buf[150]; sprintf(buf,"Scattered\n");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Scattered\n");
  }

  if (coarse_flag == 2) {
    {char buf[150]; sprintf(buf,"Number of vertices to coarsen down to: %d\n", vmax);UserWrite(buf);}
    if (print2file) fprintf(outfile, "Number of vertices to coarsen down to: %d\n", vmax);
    {char buf[150]; sprintf(buf,"Coarse graph Lanczos: ");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Coarse graph Lanczos: ");
    if (solver_flag == 1) {
      {char buf[150]; sprintf(buf,"Full orthogonalization\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, " Full orthogonalization\n");
    }
    else if (solver_flag == 2) {
      {char buf[150]; sprintf(buf,"Full orthogonalization, inverse operator\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Full orthogonalization, inverse operator\n");
    }
    else if (solver_flag == 3) {
      {char buf[150]; sprintf(buf,"Selective orthogonalization, both ends\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Selective orthogonalization, both ends\n");
    }
    else if (solver_flag == 4) {
      {char buf[150]; sprintf(buf,"Selective orthogonalization, left end only\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Selective orthogonalization, left end only\n");
    }
  }

  else if (spec) {
    {char buf[150]; sprintf(buf,"Eigensolver: ");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Eigensolver: ");
    if (coarse_flag == 1) {
      {char buf[150]; sprintf(buf,"Multilevel RQI/Symmlq\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Multilevel RQI/Symmlq\n");
      {char buf[150]; sprintf(buf,"On coarse graph: ");UserWrite(buf);}
      if (print2file) fprintf(outfile, "On coarse graph: ");
    }

    if (solver_flag == 1) {
      {char buf[150]; sprintf(buf,"Lanczos with full orthogonalization\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Lanczos with full orthogonalization\n");
    }
    else if (solver_flag == 2) {
      {char buf[150]; sprintf(buf,"Lanczos with full orthogonalization and inverted operator\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Lanczos with full orthogonalization and inverted operator\n");
    }
    else if (solver_flag == 3) {
      {char buf[150]; sprintf(buf,"Lanczos with selective orthogonalization at both ends\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Lanczos with selective orthogonalization at both ends\n");
    }
    else if (solver_flag == 4) {
      {char buf[150]; sprintf(buf,"Lanczos with selective orthogonalization at left end only\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Lanczos with selective orthogonalization at left end only\n");
    }
    else if (solver_flag == 5) {
      {char buf[150]; sprintf(buf,"Multilevel RQI/Symmlq\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Multilevel RQI/Symmlq\n");
    }

    if (coarse_flag == 1) {
      {char buf[150]; sprintf(buf,"Number of vertices to coarsen down to: %d\n", vmax);UserWrite(buf);}
      if (print2file) fprintf(outfile, "Number of vertices to coarsen down to: %d\n", vmax);
    }
  }



  else if (inert) {
    if (geomname != NULL) {
      {char buf[150]; sprintf(buf,"Geometry input file: %s. Dimensionality = %d\n", geomname, igeom);UserWrite(buf);}
      if (print2file) fprintf(outfile, "Geometry input file: %s. Dimensionality = %d\n", geomname, igeom);
    }
  }

  if (coarse_flag != 2) {
    if (KL) {
      {char buf[150]; sprintf(buf,"Local method: Kernighan-Lin\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Local method: Kernighan-Lin\n");
    }
    else {
      {char buf[150]; sprintf(buf,"No local method\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "No local method\n");
    }
  }




  {char buf[150]; sprintf(buf,"Total number of dimensions = %d\n", ndims_tot);UserWrite(buf);}
  if (print2file) fprintf(outfile, "Total number of dimensions = %d\n", ndims_tot);
  if (ndims_tot > 1) {
    if (ndims == 1) {
      {char buf[150]; sprintf(buf,"Using bisection\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Using bisection\n");
    }
    else if (ndims == 2) {
      {char buf[150]; sprintf(buf,"Using quadrisection\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Using quadrisection\n");
    }
    else if (ndims == 3) {
      {char buf[150]; sprintf(buf,"Using octasection\n");UserWrite(buf);}
      if (print2file) fprintf(outfile, "Using octasection\n");
    }
  }

  if (abs(ECHO) > 1) {
    {char buf[150]; sprintf(buf,"Parameters:\n");UserWrite(buf);}
    if (print2file) fprintf(outfile, "Parameters:\n");
    if (spec) {
      {char buf[150]; sprintf(buf,"  MAKE_CONNECTED = %s\n", true_or_false(MAKE_CONNECTED));UserWrite(buf);}
      if (print2file) fprintf(outfile, "  MAKE_CONNECTED = %s\n", true_or_false(MAKE_CONNECTED));
      {char buf[150]; sprintf(buf,"  EIGEN_TOLERANCE = %g\n", EIGEN_TOLERANCE);UserWrite(buf);}
      if (print2file) fprintf(outfile, "  EIGEN_TOLERANCE = %g\n", EIGEN_TOLERANCE);
      {char buf[150]; sprintf(buf,"  LANCZOS_SO_INTERVAL = %d\n", LANCZOS_SO_INTERVAL);UserWrite(buf);}
      if (print2file) fprintf(outfile, "  LANCZOS_SO_INTERVAL = %d\n", LANCZOS_SO_INTERVAL);
      {char buf[150]; sprintf(buf,"  BISECTION_SAFETY = %g\n", BISECTION_SAFETY);UserWrite(buf);}
      if (print2file) fprintf(outfile, "  BISECTION_SAFETY = %g\n", BISECTION_SAFETY);

      if (PERTURB) {
        {char buf[150]; sprintf(buf,"  PERTURB = TRUE, NPERTURB = %d, PERTURB_MAX = %g\n", NPERTURB, PERTURB_MAX);UserWrite(buf);}
        if (print2file) fprintf(outfile, "  PERTURB = TRUE, NPERTURB = %d, PERTURB_MAX = %g\n",
                                NPERTURB, PERTURB_MAX);
      }
      else {
        {char buf[150]; sprintf(buf,"  PERTURB = FALSE\n");UserWrite(buf);}
        if (print2file) fprintf(outfile, "  PERTURB = FALSE\n");
      }
    }

    if (coarse_flag == 1) {
      {char buf[150]; sprintf(buf,"  COARSE_NLEVEL_RQI = %d\n", COARSE_NLEVEL_RQI);UserWrite(buf);}
      if (print2file) fprintf(outfile, "  COARSE_NLEVEL_RQI = %d\n", COARSE_NLEVEL_RQI);
    }

    if (KL) {
      {char buf[150]; sprintf(buf,"  KL_RANDOM = %s\n", true_or_false(KL_RANDOM));UserWrite(buf);}
      if (print2file) fprintf(outfile, "  KL_RANDOM = %s\n", true_or_false(KL_RANDOM));
      if (KL_METRIC == 1) {
        {char buf[150]; sprintf(buf,"  KL_METRIC = Cuts\n");UserWrite(buf);}
        if (print2file) fprintf(outfile, "  KL_METRIC = Cuts\n");
      }
      else if (KL_METRIC == 2) {
        {char buf[150]; sprintf(buf,"  KL_METRIC = Hops\n");UserWrite(buf);}
        if (print2file) fprintf(outfile, "  KL_METRIC = Hops\n");
      }
      {char buf[150]; sprintf(buf,"  KL_NTRIES_BAD = %d\n", KL_NTRIES_BAD);UserWrite(buf);}
      if (print2file) fprintf(outfile, "  KL_NTRIES_BAD = %d\n", KL_NTRIES_BAD);
      {char buf[150]; sprintf(buf,"  KL_BAD_MOVES = %d\n", KL_BAD_MOVES);UserWrite(buf);}
      if (print2file) fprintf(outfile, "  KL_BAD_MOVES = %d\n", KL_BAD_MOVES);
    }

    if (coarse_flag == 2) {
      {char buf[150]; sprintf(buf,"  COARSE_NLEVEL_KL = %d\n", COARSE_NLEVEL_KL);UserWrite(buf);}
      if (print2file) fprintf(outfile, "  COARSE_NLEVEL_KL = %d\n", COARSE_NLEVEL_KL);
    }
  }
  {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
  if (print2file) fprintf(outfile, "\n");
}
