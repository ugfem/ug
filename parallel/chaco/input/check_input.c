// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <math.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"

/* Check input parameters and graph for consistency. */
int check_input(graph, nvtxs, nedges, spec, inert, KL, mediantype,
                solver_flag, coarse_flag, vmax, ndims, ndims_tot, igeom, coords,
                scatt, randm, lin, graphname)
struct vtx_data **graph;        /* linked lists of vertex data */
int nvtxs;                      /* number of vertices */
int nedges;                     /* number of edges */
int spec;                       /* invoke spectral method? */
int inert;                      /* invoke inertial method? */
int KL;                         /* invoke Kernighan-Lin? */
int mediantype;                 /* median-finding option for spec */
int solver_flag;                /* type of eigensolver */
int coarse_flag;                /* flag for recursive coarsening  */
int vmax;                       /* smallest acceptable coarsened nvtxs */
int ndims;                      /* partitioning level */
int ndims_tot;                  /* total number of cuts to make */
int igeom;                      /* geometric dimension for inertial method */
float **coords;                 /* coordinates for inertial method */
int scatt;                      /* use scattered global partitioner? */
int randm;                      /* use randm global partitioner? */
int lin;                        /* use linear global partitioner? */
char *graphname;                /* graph input file name */
{
  int nmethods;                 /* number of methods selected */
  int flag;                     /* does input data pass all the tests? */
  int flag_graph;               /* does graph check out OK? */
  int check_graph();

  flag = TRUE;

  if (ndims_tot < 1 || ndims_tot > MAXDIMS_TOT) {
    {char buf[150]; sprintf(buf,"Number of hypercube dimensions = %d, should be in [1,%d].\n", ndims_tot, MAXDIMS_TOT);UserWrite(buf);}
    flag = FALSE;
  }

  if (ndims < 1 || ndims > MAXDIMS) {
    {char buf[150]; sprintf(buf,"Partitioning at each step = %d, should be in [1,%d].\n", ndims, MAXDIMS);UserWrite(buf);}
    flag = FALSE;
  }

  if (ndims > ndims_tot) {
    {char buf[150]; sprintf(buf,"Partitioing step %d is greater than hypercube size %d.\n", ndims, ndims_tot);UserWrite(buf);}
    flag = FALSE;
  }

  nmethods = 0;
  if (spec) ++nmethods;
  if (inert) ++nmethods;
  if (scatt) ++nmethods;
  if (randm) ++nmethods;
  if (lin) ++nmethods;

  if (nmethods == 0) {
    {char buf[150]; sprintf(buf,"No global partitioning option selected.\n");UserWrite(buf);}
    flag = FALSE;
  }

  else if (nmethods > 1) {
    {char buf[150]; sprintf(buf,"Multiple global partition options selected.\n");UserWrite(buf);}
    flag = FALSE;
  }


  if (coarse_flag == 1 && !spec) {
    {char buf[150]; sprintf(buf,"Spectral partitioning not active, but coarse_flag = %d\n", coarse_flag);UserWrite(buf);}
    flag = FALSE;
  }

  if (coarse_flag == 2 && (!spec || !KL)) {
    {char buf[150]; sprintf(buf,"When using multilevel KL/spec method, both spectral and KL must be active\n");UserWrite(buf);}
    flag = FALSE;
  }

  if (coarse_flag == 2 && vmax <= ndims) {
    {char buf[150]; sprintf(buf,"When using multilevel KL/spec method, both spectral and KL must be active\n");UserWrite(buf);}
    flag = FALSE;
  }

  if (spec) {
    if (mediantype < 0 || mediantype > 3) {
      {char buf[150]; sprintf(buf,"Bad median finding option (%d).\n", mediantype);UserWrite(buf);}
      flag = FALSE;
    }
    if (solver_flag < 1 || solver_flag > 8) {
      {char buf[150]; sprintf(buf,"Invalid eigensolver selected (%d)\n", solver_flag);UserWrite(buf);}
      flag = FALSE;
    }
    if (coarse_flag == 1 || coarse_flag == 2) {
      if (vmax < 2) {
        {char buf[150]; sprintf(buf,"Bad number of vertices to coarsen down to (%d)\n", vmax);UserWrite(buf);}
        flag = FALSE;
      }
      else if (vmax <= ndims) {
        {char buf[150]; sprintf(buf,"Number of vertices of coarse graph (%d) must be greater than\n", vmax);UserWrite(buf);}
        {char buf[150]; sprintf(buf,"  the number of eigenvectors (%d)\n", ndims);UserWrite(buf);}
        flag = FALSE;
      }
    }
  }

  if (inert) {
    if (igeom < 1 || igeom > 3) {
      {char buf[150]; sprintf(buf,"Geometry must be 1-, 2- or 3-dimensional for inertial method\n");UserWrite(buf);}
      flag = FALSE;
    }
    if (igeom > 0 && coords == NULL) {
      {char buf[150]; sprintf(buf,"No coordinates given for inertial method\n");UserWrite(buf);}
      flag = FALSE;
    }
    else if (igeom > 0 && coords[0] == NULL) {
      {char buf[150]; sprintf(buf,"No X-coordinates given for inertial method\n");UserWrite(buf);}
      flag = FALSE;
    }
    else if (igeom > 1 && coords[1] == NULL) {
      {char buf[150]; sprintf(buf,"No Y-coordinates given for inertial method\n");UserWrite(buf);}
      flag = FALSE;
    }
    else if (igeom > 2 && coords[2] == NULL) {
      {char buf[150]; sprintf(buf,"No Z-coordinates given for inertial method\n");UserWrite(buf);}
      flag = FALSE;
    }
  }


  /* Now check for consistency in the graph.  Every edge should appear twice. */
  if (graph != NULL)  {
    flag_graph = check_graph(graph, nvtxs, nedges);
    if (!flag_graph) {
      if (graphname != NULL) {char buf[150]; sprintf(buf,"Errors in graph input file %s\n", graphname);UserWrite(buf);}
      else{char buf[150]; sprintf(buf,"Errors in graph\n");UserWrite(buf);}

    }
  }

  else {
    /* Only allowed if simple or inertial w/o KL and no weights. */
    flag_graph = TRUE;
    if ((!inert & !scatt & !randm & !lin) || KL) {
      {char buf[150]; sprintf(buf,"No graph input.  Only allowed for inertial method without KL.\n");UserWrite(buf);}
      flag_graph = FALSE;
    }
  }

  flag = flag && flag_graph;

  return(flag);
}
