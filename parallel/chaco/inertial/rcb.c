// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>
#include <stdio.h>
#include "../main/structs.h"
#include "../main/defs.h"


void rcbi(graph, nvtxs, ndims, igeom, coords, sets, goal, using_vwgts,
          part_type)
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vtxs in graph */
int ndims;                      /* number of cuts to make at once */
int igeom;                      /* 2 or 3 dimensional geometry? */
float **coords;                 /* x, y and z coordinates of vertices */
short *sets;                    /* set each vertex gets assigned to */
double *goal;                   /* desired set sizes */
int using_vwgts;                /* are vertex weights being used? */
int part_type;
{
  extern DEBUG_TRACE;           /* trace the execution of the code */
  extern double rcb_time;       /* time spend in inertial calculations */
  double time;                  /* timing parameter */
  double seconds();
  void rcb2d(), rcb3d();

  time = seconds();

  if (DEBUG_TRACE > 0) {
    {char buf[150]; sprintf(buf,"Entering rcb, nvtxs = %d\n", nvtxs);UserWrite(buf);}
  }

  else if (igeom == 2) rcb2d(graph, nvtxs, ndims, coords[0], coords[1],
                             sets, goal, using_vwgts, part_type);

  else if (igeom == 3) rcb3d(graph, nvtxs, ndims, coords[0], coords[1],
                             coords[2], sets, goal, using_vwgts, part_type);
  rcb_time += seconds() - time;
}
