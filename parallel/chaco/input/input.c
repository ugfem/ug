// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        <string.h>
#include        "../main/defs.h"
#include        "../main/params.h"


static void get_eigensolver();

void input_queries(fin, fingeom, inname, geomname, ndims, ndims_tot, assignname,
                   spec, inert, KL, solver_flag, coarse_flag, vmax, scatt, randm,
                   lin, outfile)

FILE **fin;                     /* input file */
FILE **fingeom;                 /* geometry input file (for inertial method)  */
char *inname;                   /* name of graph input file */
char *geomname;                 /* name of geometry input file */
int *ndims;                     /* number of divisions at each stage */
int *ndims_tot;                 /* target number of hypercube dimensions */
char *assignname;               /* name of assignment output file */
int *spec;                      /* invoke spectral method? */
int *inert;                     /* invoke inertial method? */
int *KL;                        /* invoke Kernighan-Lin? */
int *solver_flag;               /* which eigen-solver should I use? */
int *coarse_flag;               /* should I use multilevel solver? */
int *vmax;                      /* if so, how far should I coarsen? */
int *scatt;                     /* invoke scattered partitioning? */
int *randm;                     /* invoke random partitioning? */
int *lin;                       /* invoke linear partitioning? */
FILE **outfile;                 /* file for outputing run results */
{
  extern int OUTPUT_ASSIGN;     /* whether to write assignments to file */
  extern int ECHO;              /* copy input to screen? results to file? */
  char outname[NAME_LENGTH];    /* name of file for outputting run results */
  int glob_method;              /* what global partitioning strategy to use? */
  int loc_method;               /* what local refinement strategy to use? */
  int input_int();

  /* Name and open input file. */
  *fin = NULL;
  while (*fin == NULL) {
    {char buf[150]; sprintf(buf,"Graph input file: ");UserWrite(buf);}
    scanf("%s", inname);

    *fin = fopen(inname, "r");
    if (*fin == NULL) {
      {char buf[150]; sprintf(buf,"Graph file `%s' not found.\n", inname);UserWrite(buf);}
    }
  }

  /* Name output file. */
  if (OUTPUT_ASSIGN) {
    {char buf[150]; sprintf(buf,"Assignment output file: ");UserWrite(buf);}
    scanf("%s", assignname);
  }

  /* Name output file. */
  if (ECHO < 0) {
    {char buf[150]; sprintf(buf,"File name for saving run results: ");UserWrite(buf);}
    scanf("%s", outname);
    *outfile = fopen(outname, "w");
  }
  else *outfile = NULL;

  /* Get global method, if any. */
  glob_method = 0;
  while (glob_method < 1 || glob_method > 6) {
    {char buf[150]; sprintf(buf,"Global partitioning method:\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (1) Multilevel\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (2) Spectral\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (3) Inertial\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (4) Linear\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (5) Random\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (6) Scattered\n");UserWrite(buf);}
    glob_method = input_int();
  }

  /* Initialize then set the flags */
  *scatt = *randm = *lin = *inert = *spec = *KL = FALSE;
  *coarse_flag = *solver_flag = 0;

  if (glob_method == 6) {
    *scatt = TRUE;
  }
  else if (glob_method == 5) {
    *randm = TRUE;
  }
  else if (glob_method == 4) {
    *lin = TRUE;
  }
  else if (glob_method == 3) {
    *inert = TRUE;
    *fingeom = NULL;
    while (*fingeom == NULL) {
      {char buf[150]; sprintf(buf,"Geometry input file name: ");UserWrite(buf);}
      scanf("%s", geomname);

      *fingeom = fopen(geomname, "r");
      if (*fingeom == NULL) {
        {char buf[150]; sprintf(buf,"Geometry input file `%s' not found.\n", geomname);UserWrite(buf);}
      }
    }
  }
  else if (glob_method == 2) {
    *spec = TRUE;
    *fingeom = NULL;
    get_eigensolver(solver_flag, coarse_flag, vmax);
  }
  else if (glob_method == 1) {
    *spec = TRUE;
    *KL = TRUE;
    *fingeom = NULL;
    *coarse_flag = 2;
    *vmax = 0;
    while (*vmax <= 1) {
      {char buf[150]; sprintf(buf,"Number of vertices to coarsen down to: ");UserWrite(buf);}
      *vmax = input_int();
    }
    get_eigensolver(solver_flag, coarse_flag, vmax);
  }

  /* Get local method, if any */
  loc_method = 0;
  if (glob_method == 1) loc_method = 1;
  while (loc_method < 1 || loc_method > 2) {
    {char buf[150]; sprintf(buf,"Local refinement method:\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (1) Kernighan-Lin\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (2) None\n");UserWrite(buf);}
    loc_method = input_int();
  }

  if (loc_method == 1) {
    *KL = TRUE;
  }

  /* Get total number of dimensions in which to partition. */
  *ndims_tot = 0;
  while (*ndims_tot < 1 || *ndims_tot > MAXDIMS_TOT) {
    {char buf[150]; sprintf(buf,"Total number of target hypercube dimensions: ");UserWrite(buf);}
    *ndims_tot = input_int();
    if (*ndims_tot < 1 || *ndims_tot > MAXDIMS_TOT) {
      {char buf[150]; sprintf(buf," Number of divisions must be between 1 and %d\n", MAXDIMS_TOT);UserWrite(buf);}
    }
  }

  /* Get number of dimensions in which to partition at each level. For inertial
     must do bisection hence ndims = 1. Also restricted by ndims_tot. */
  if (*ndims_tot == 1) {*ndims = 1;}
  else if (*ndims_tot == 2) {
    {char buf[150]; sprintf(buf,"Partitioning dimension: \n");UserWrite(buf);}
    while (*ndims < 1 || *ndims > 2) {
      {char buf[150]; sprintf(buf,"   (1) Bisection\n");UserWrite(buf);}
      {char buf[150]; sprintf(buf,"   (2) Quadrisection\n");UserWrite(buf);}
      *ndims = input_int();
    }
  }
  else if (*ndims_tot >= 3) {
    {char buf[150]; sprintf(buf,"Partitioning dimension: \n");UserWrite(buf);}
    while (*ndims < 1 || *ndims > 3) {
      {char buf[150]; sprintf(buf,"   (1) Bisection\n");UserWrite(buf);}
      {char buf[150]; sprintf(buf,"   (2) Quadrisection\n");UserWrite(buf);}
      {char buf[150]; sprintf(buf,"   (3) Octasection\n");UserWrite(buf);}
      *ndims = input_int();
    }
  }
}


/* Routine to accept appropriate flags/parameters for eigensolver */
static void get_eigensolver(solver_flag, coarse_flag, vmax)
int *solver_flag;               /* which eigen-solver should I use? */
int *coarse_flag;               /* should I use multilevel solver? */
int *vmax;                      /* if so, how far should I coarsen? */
{
  int max_flag;                 /* largest allowed value for solver_flag */
  void get_eigensolver();

  if (*coarse_flag == 0) max_flag = 5;
  else max_flag = 4;

  *solver_flag = 0;
  while (*solver_flag < 1 || *solver_flag > max_flag) {
    if (max_flag == 4) {char buf[150]; sprintf(buf,"Coarse graph eigensolver:\n");UserWrite(buf);}
    else{char buf[150]; sprintf(buf,"Eigensolver:\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (1) Lanczos with full orthogonalization\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (2) Lanczos with full orthogonalization and inverted operator\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (3) Lanczos with selective orthogonalization at both ends\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"  (4) Lanczos with selective orthogonalization at left end only\n");UserWrite(buf);}
    if (max_flag > 4) {char buf[150]; sprintf(buf,"  (5) Multilevel RQI/Symmlq\n");UserWrite(buf);}
    *solver_flag = input_int();
  }

  if (*solver_flag == 5) {
    *coarse_flag = 1;
    *vmax = 0;
    while (*vmax <= 0) {
      {char buf[150]; sprintf(buf,"Number of vertices to coarsen down to: ");UserWrite(buf);}
      *vmax = input_int();
    }
    get_eigensolver(solver_flag, coarse_flag, vmax);
  }
}
