// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/params.h"
#include        "../main/structs.h"
#include        "../main/defs.h"

/* Idea:
   'buckets[i][j]' is a set of buckets to sort moves from i to j.
   listspace[i] is space for lists in buckets[i][j].
   Loop through all nonequal pairs [i][j], taking the first element
   in each list.  Compare them all to find the largest allowed move.
   Make that move, and save it in movelist.
 */


void bucketsorts(graph, nvtxs, buckets, listspace, dvals, sets, maxdval,
                 nsets, parity, hops, bspace, list_length, npass)
struct vtx_data **graph;        /* graph data structure */
int nvtxs;                      /* number of vertices */
struct bilist ****buckets;      /* array of lists for bucket sort */
struct bilist **listspace;      /* list data structure for each vertex */
int **dvals;                    /* d-values for each vertex for removing */
short *sets;                    /* processor each vertex is assigned to */
int maxdval;                    /* maximum possible dvalue for a vertex */
int nsets;                      /* number of sets being divided into */
int parity;                     /* work in forward or backward direction? */
short (*hops)[MAXSETS];                 /* hop cost between sets */
int *bspace;                    /* indices for randomly ordering vtxs */
int list_length;                /* number of values in bspace to work with */
int npass;                      /* which pass through KL is this? */
{
  extern int KL_RANDOM;         /* use randomness in KL? */
  extern int KL_ONLY_BNDY;      /* only consider edges on set boundaries? */
  extern int KL_UNDO_LIST;      /* only sort vertices who have moved. */
  struct bilist **bptr;         /* loops through set of buckets */
  struct bilist *lptr;          /* pointer to an element in listspace */
  float *ewptr;                 /* loops through edge weights */
  int *bsptr;                   /* loops through bspace */
  int *edges;                   /* edge list for a vertex */
  int *sptr;                    /* loops through dvals */
  int myset;                    /* set that current vertex belongs to */
  int newset;                   /* set current vertex could move to */
  int set;                      /* set that neighboring vertex belongs to */
  int weight;                   /* edge weight for a particular edge */
  int using_ewgts;              /* are edge weights being used? */
  int vtx;                      /* vertex in graph */
  short myhop;                  /* hops associated with current vertex */
  int internal;                 /* vertex only has neighbors in same set? */
  int i, j, l;                  /* loop counters */
  void randomize(), add2bilist();

  /* For each vertex, compute d-values for each possible transition. */
  /* Then store them in each appropriate bucket. */

  if (!KL_UNDO_LIST || (nvtxs == list_length)) {
    /* Empty all the buckets. */
    bptr = buckets[0][1];
    for (i=nsets*(nsets-1)*(2*maxdval+1); i; i--) *bptr++ = NULL;
  }

  /* Randomize the order of the vertices */

  if (npass == 1 || (KL_UNDO_LIST && nvtxs == list_length) ||
      (!KL_UNDO_LIST && !KL_RANDOM)) {
    /* Don't need to reoder if about to randomize. */
    bsptr = bspace;
    if (parity) for (i=nvtxs; i; i--) *bsptr++ = nvtxs-i+1;
    else for (i=nvtxs; i; i--) *bsptr++ = i;
  }
  if (KL_RANDOM) randomize(bspace-1, list_length);

  if (!KL_ONLY_BNDY && !KL_UNDO_LIST) {
    sptr = &dvals[1][0];
    for (j=(nsets-1)*nvtxs; j; j--) *sptr++ = 0;
  }

  /* Now compute d-vals by seeing which sets neighbors belong to. */
  using_ewgts = (graph[1]->ewgts != NULL);
  weight = 1;
  bsptr = bspace;
  for (i=0; i<list_length; i++) {
    vtx = *bsptr++;
    myset = sets[vtx];

    if (KL_ONLY_BNDY) {         /* Check to see if on set boundary. */
      edges = graph[vtx]->edges;
      for (j=graph[vtx]->nedges-1; j; j--) {
        if (sets[*(++edges)] != myset) break;
      }
      internal = !j;
    }

    if (!KL_ONLY_BNDY || !internal) {
      if (KL_ONLY_BNDY || KL_UNDO_LIST) {
        for (j=0; j<nsets-1; j++) dvals[vtx][j] = 0;
      };

      /* First count the neighbors in each set. */
      edges = graph[vtx]->edges;
      if (using_ewgts) ewptr = graph[vtx]->ewgts;
      for (j=graph[vtx]->nedges-1; j; j--) {
        set = sets[*(++edges)];
        if (set < 0) set = -set - 1;
        if (using_ewgts) weight = .5 + *(++ewptr);
        myhop = hops[myset][set];

        l = 0;
        for (newset=0; newset<nsets; newset++) {
          if (newset != myset) {
            dvals[vtx][l] += weight*(myhop - hops[newset][set]);
            l++;
          }
        }
      }

      /* Now add to appropriate buckets. */
      l = 0;
      for (newset=0; newset<nsets; newset++) {
        if (newset != myset) {
          lptr = listspace[l];
          add2bilist(&lptr[vtx], &buckets[myset][newset][dvals[vtx][l]+maxdval]);
          ++l;
        }
      }
    }
    else {      /* Flag vertex as if already moved. */
      sets[vtx] = -sets[vtx] - 1;
    }
  }
}
