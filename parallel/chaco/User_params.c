// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "main/defs.h"

int CHECK_INPUT = TRUE;         /* check the input for consistency? */

int ECHO = 0;                   /* print input/param options? to file? (-2..2) */

long RANDOM_SEED = 7654321L;    /* seed for random number generators */

double EIGEN_TOLERANCE = 1e-3;  /* convergence criteria for eigen-tolerance */
double BISECTION_SAFETY = 10;   /* divides Lanczos bisection tolerance function */

int MAKE_CONNECTED = TRUE;      /* connect graph if using spectral method? */

int NSQRTS = 1000;              /* # square roots to precompute if coarsening */

int LANCZOS_SO_INTERVAL = 10;   /* interval between orthog checks in SO */
int LANCZOS_SO_TIME = 0;        /* perform detailed timing of Lanczos_SO? (0..1) */

int MAPPING_TYPE = 2;           /* how to map from eigenvectors to partition */
                                /* 0 => cut at origin, 1 => cut at medians */
                                /* 2 => recursive median, 3 => bipartite match */

int ARCH_GOAL = 0;              /* which sets to choose for subgoals?    */
/* 0 =>array  1 =>hypercube architecture */

float COARSEN_RATIO_MIN = .7;   /* min vtx reduction at each coarsen stage */
int COARSE_NLEVEL_RQI = 2;      /* # levels between RQI calls in uncoarsening */

int KL_RANDOM = TRUE;           /* use randomness in Kernighan-Lin? */
int KL_ONLY_BNDY = FALSE;       /* only consider moving vtxs on boundary? */
int KL_UNDO_LIST = TRUE;        /* only re-sort vertices changed in prior pass? */
int KL_METRIC = 1;              /* KL interset cost: 1=>cuts, 2=>hypercube hops */
/* 3=>array hops                                */
int KL_NTRIES_BAD = 2;          /* # unhelpful passes before quitting KL */
int KL_BAD_MOVES = 40;          /* number of unhelpful moves in a row allowed */
int COARSE_NLEVEL_KL = 3;       /* # levels between KL calls in uncoarsening */

int PERTURB = TRUE;             /* randomly perturb matrix in spectral method? */
int NPERTURB = 2;               /* if so, how many edges to modify? */
double PERTURB_MAX = 1.0e-4;    /* largest value for perturbation */

int OPT3D_NTRIES = 5;           /* # local opts to look for global min in opt3d */

int DEBUG_CONNECTED = 0;        /* debug flag for connected components (0..1) */
int DEBUG_PERTURB = 0;          /* debug flag for matrix perturbation (0..1) */
int DEBUG_ASSIGN = 0;           /* debug flag for assignment to sets (0..1) */
int DEBUG_INERTIAL = 0;         /* debug flag for inertial method (0..1) */
int DEBUG_OPTIMIZE = 0;         /* debug flag for optimization/rotation (0..2) */
int DEBUG_BPMATCH = 0;          /* debug flag for bipartite matching code (0..2) */
int DEBUG_COARSEN = 0;          /* debug flag for coarsening/uncoarsening (0..1) */
int DEBUG_EVECS = 0;            /* debug flag for eigenvector generation (0..5) */
int DEBUG_KL = 0;               /* debug flag for Kernighan-Lin (0..3) */
int DEBUG_MEMORY = 0;           /* debug flag for smalloc/sfree (0..2) */
int DEBUG_INPUT = 0;            /* debug flag for having read input files (0..1) */
int DEBUG_TRACE = 0;            /* trace main execution path (0..1) */
int DEBUG_RECURSE = 0;          /* trace path of recursion (0..1) */
int DEBUG_GRAPH = 0;            /* trace consistency of graph (0..1) */
int DEBUG_ARRAY = 0;        /* debug flag for array architecture */

int WARNING_EVECS = 2;          /* warnings in eigenvector generation (0..3)*/
double WARNING_ORTHTOL = 2;     /* warning if Ares and bjitol have this ratio */
double WARNING_MISTOL = 100;    /* warning if Ares and bjitol have this ratio */
double WARNING_SRESTOL = 5;     /* warning if evec of T is inaccurate. */

int OUTPUT_TIME = 2;            /* At what level to display timings (0..2) */
int OUTPUT_METRICS = 0;         /* controls displaying of results (0..3) */
int OUTPUT_ASSIGN = 0;          /* whether to write assignments to file (0..1) */

double *MEM_OK = (double *)0x1L; /* variable for memory overflow exeception */
Heap   *heap = NULL;     /* pointer to heap of multigrid */
