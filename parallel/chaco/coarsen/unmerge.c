// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"

static void displace();

void unmerge(vecs, ndims, graph, cgraph, merged, reduction, nmerged, v2cv, space,
             using_ewgts)
double **vecs;                  /* eigenvectors being modified */
int ndims;                      /* number of eigenvectors */
struct vtx_data **graph;        /* array of vtx data for graph */
struct vtx_data **cgraph;       /* coarsened version of graph */
struct ipairs *merged;          /* vtx pairs that get merged */
struct fpairs *reduction;       /* multipliers used to reduce edge weights */
int nmerged;                    /* number of vertices that were merged */
int *v2cv;                      /* mapping from original to coarse vtxs */
short *space;                   /* space for max_cneighbors+1 shorts */
int using_ewgts;                /* are edge weights being used in graph? */
{
  double k1, k2, k12;           /* (summed) edge weights for nonparallel edges */
  double k1z[MAXDIMS+1];        /* sums of k*vecs for one vertex in edge */
  double k2z[MAXDIMS+1];        /* sums of k*vecs for other vertex in edge */
  double z0[MAXDIMS+1];         /* vec value of vertex being unmerged */
  double wt1, wt2, wtsum;       /* edge weights and their sum */
  double ewgt;                  /* edge weight */
  int v1, v2;                   /* vertex numbers */
  int m1, m2;                   /* vertex weights of contracting vertices */
  int cv, cv0;                  /* contracted vertex number */
  int neighbor;                 /* neighboring vertex */
  int square;                   /* square w/ two matching edges been found? */
  int i, j, k, l, idim;         /* loop counter */

  /* Unmerge the edges in reverse order. */
  k12 = 1;
  for (i=nmerged-1; i>=0; i--) {
    /* Read all the easy information about the edge. */
    v1 = min(merged[i].val1, merged[i].val2);
    v2 = max(merged[i].val1, merged[i].val2);
    cv = v2cv[v1];
    m1 = graph[v1]->vwgt;
    m2 = graph[v2]->vwgt;

    if (using_ewgts) {
      for (j=1; j<graph[v2]->nedges && graph[v2]->edges[j] != v1; j++) ;
      k12 = graph[v2]->ewgts[j];
    }

    for (idim=1; idim<=ndims; idim++) z0[idim] = vecs[idim][v1];

    /* Figure out what the spring constants were before merging these */
    /* vertices.  Then use this information to compute z values for */
    /* the unmerging vertices. */

    /* Treat three classes of edges differently.  Edges forming a square */
    /* with another matching edge and those forming a triangle are unchanged */
    /* in coarsening process, so they can be factored back out. */
    /* The rest are modified by the multipliers.  */

    /* Consider all adjacent edges in graph and determine where they get */
    /* mapped to in coarse graph.  By looking at neighbors of both merged */
    /* vertices, I can distinguish these three types of edges. */
    for (k=1; k<cgraph[cv]->nedges; k++) space[k] = 0;
    for (j=1; j<graph[v1]->nedges; j++) {
      neighbor = graph[v1]->edges[j];
      if (neighbor != v2) {             /* If == v2, edge disappears. */
        cv0 = v2cv[neighbor];
        for (k=1; k<cgraph[cv]->nedges && cgraph[cv]->edges[k]!=cv0; k++) ;
        space[k] = 1;
      }
    }

    for (j=1; j<graph[v2]->nedges; j++) {
      neighbor = graph[v2]->edges[j];
      if (neighbor != v1) {
        cv0 = v2cv[neighbor];
        for (k=1; k<cgraph[cv]->nedges && cgraph[cv]->edges[k]!=cv0; k++) ;
        if (space[k] == 0) space[k] = 2;
        else if (space[k] == 1) space[k] = 3;
      }
    }

    /* space=1 or 2 => adjacent only to v1 or v2 in coarsened graph. */
    /* space=3 => adjacent to both (so triangle or square). */

    /* Now compute unmerged spring constants of all the neighboring edges, */
    /* and sum into different sets. */

    /* If space = 1 or 2, then edge weight goes to either cwgt/reduction or
       it goes to weight in original graph.
       If vwgt < 0, then use weight in original graph.
       If space = 3, then either triangle or square.
       Triangle separates into two edges with proportions = original proportions.
       Square separates into original edges.
     */
    k1 = k2 = 0;
    for (idim=1; idim<=ndims; idim++) k1z[idim] = k2z[idim] = 0;
    cgraph[cv]->vwgt = -1;
    for (j=1; j<graph[v1]->nedges; j++) {
      neighbor = graph[v1]->edges[j];
      if (neighbor != v2) {             /* If == v2, edge disappears. */
        cv0 = v2cv[neighbor];
        if (cgraph[cv0]->vwgt < 0) {
          /* Neighboring vtx already modified, revert to original edge wt. */
          if (using_ewgts) ewgt = graph[v1]->ewgts[j];
          else ewgt = 1;
        }
        else {
          /* Find coarse edge weight for possible adjustment */
          for (k=1; k<cgraph[cv]->nedges && cgraph[cv]->edges[k]!=cv0; k++) ;
          if (space[k] == 1) {          /* Just adjacent to v1. */
            ewgt = cgraph[cv]->ewgts[k]/reduction[i].val1;
          }
          else if (space[k] == 3) {             /* Square or triangle edge. */
            /* If square, revert to original edge weight. */
            /* Else if triangle, keep original ratio, but unscale. */
            square = FALSE;
            for (l=1; l<graph[v2]->nedges; l++) {       /* Look for square. */
              if (v2cv[graph[v2]->edges[l]] == cv0 &&
                  graph[v2]->edges[l] != neighbor) break;
            }
            if (l < graph[v2]->nedges) square = TRUE;
            else {              /* Look for other square. */
              for (l=1; l<graph[v1]->nedges; l++) {
                if (v2cv[graph[v1]->edges[l]] == cv0 &&
                    graph[v1]->edges[l] != neighbor) break;
              }
              if (l < graph[v1]->nedges) square = TRUE;
            }
            if (square) {               /* Square! */
              if (using_ewgts) ewgt = graph[v1]->ewgts[j];
              else ewgt = 1;
            }
            else {              /* Triangle. */
              for (l=1; l<graph[v2]->nedges &&
                   graph[v2]->edges[l]!=neighbor; l++) ;
              if (using_ewgts) {
                wt1 = graph[v1]->ewgts[j];
                wt2 =  graph[v2]->ewgts[l];
              }
              else {
                wt1 = wt2 = 1;
              }
              wtsum = wt1 + wt2;
              /* Might as well do other edge while I'm here. */
              /* Then alter space[k] to indicate that it's done. */
              ewgt = wt2*cgraph[cv]->ewgts[k]/wtsum;
              /*
                 if (using_ewgts) ewgt = graph[v2]->ewgts[l];
                 else ewgt = 1;
               */
              k2 += ewgt;
              for (idim=1; idim<=ndims; idim++) {
                k2z[idim] += ewgt*vecs[idim][neighbor];
              }
              space[k] = 0;
              ewgt = wt1*cgraph[cv]->ewgts[k]/wtsum;
              /*
                 if (using_ewgts) ewgt = graph[v1]->ewgts[j];
                 else ewgt = 1;
               */
            }
          }
        }
        k1 += ewgt;
        for (idim=1; idim<=ndims; idim++) {
          k1z[idim] += ewgt*vecs[idim][neighbor];
        }
      }
    }
    /* Now do the same for edges incident to v2. */
    for (j=1; j<graph[v2]->nedges; j++) {
      neighbor = graph[v2]->edges[j];
      if (neighbor != v1) {             /* If == v1, edge disappears. */
        cv0 = v2cv[neighbor];
        if (cgraph[cv0]->vwgt < 0) {
          /* Neighboring vtx already modified, revert to original edge wt. */
          if (using_ewgts) ewgt = graph[v2]->ewgts[j];
          else ewgt = 1;
        }
        else {
          /* Find coarse edge weight for possible adjustment */
          for (k=1; k<cgraph[cv]->nedges && cgraph[cv]->edges[k]!=cv0; k++) ;
          if (space[k] == 2) {          /* Just adjacent to v2. */
            ewgt = cgraph[cv]->ewgts[k]/reduction[i].val2;
          }
          else if (space[k] == 3) {
            /* Must be square, so revert to original edge weight. */
            if (using_ewgts) ewgt = graph[v2]->ewgts[j];
            else ewgt = 1;
          }
          else if (space[k] == 0) ewgt = 0;             /* Already done. */
        }
        k2 += ewgt;
        for (idim=1; idim<=ndims; idim++) {
          k2z[idim] += ewgt*vecs[idim][neighbor];
        }
      }
    }

    /* Now use spring constants and weights to compute proportional */
    /* new z displacements. */
    for (idim=1; idim<=ndims; idim++) {
      displace(m1, m2, k1, k2, k1z[idim], k2z[idim], z0[idim], k12,
               &vecs[idim][v1], &vecs[idim][v2]);
    }
  }
}



static void displace(m1, m2, k1, k2, k1z, k2z, z0, k12, z1, z2)
/* Compute the displacement of the unmerged vertices. */
/* Assume center-of-mass stays unchanged, and sum  k*(dz)^2 is minimized. */
/* Note that the spring constants should be those in the unmerged graph. */
int m1, m2;             /* masses of vtxs being unjoined */
double k1, k2;          /* sum of spring constants from v1 and v2 */
double k1z, k2z;        /* sum of k*z values from v1 and v2 */
double z0;              /* displacement of merged vertex */
double k12;             /* spring constant of merged edge */
double *z1, *z2;        /* computed displacements for v1 and v2 */
{
  double r1, r2;        /* ratios of masses */
  double denom;         /* demoninator */

  r1 =  m1/((double) m1+m2);
  r2 = 1.0 - r1;

  denom = 1/(r2*r2*k1 + r1*r1*k2 + k12);
  *z1 = denom*(r2*r2*k1z - r1*r2*k2z + r1*z0*k2 + z0*k12);
  *z2 = denom*(r1*r1*k2z - r1*r2*k1z + r2*z0*k1 + z0*k12);
}
