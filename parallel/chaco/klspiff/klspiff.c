// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/structs.h"
#include        "../main/params.h"
#include        "../main/defs.h"

/* Think hard about space.
   Put active list into each routine.

   It may be possible to overlap dvals with active in kl, requiring
   a total of nvtxs space.
 */
static void free_kl();


void klspiff(graph, nvtxs, sets, nsets, hops, goal, vwgt_max, maxdeg)
struct vtx_data **graph;        /* list of graph info for each vertex */
int nvtxs;                      /* number of vertices in graph */
short *sets;                    /* local partitioning of vtxs */
int nsets;                      /* number of sets at each level */
short (*hops)[MAXSETS];         /* hop cost between sets */
double *goal;                   /* desired set sizes */
int vwgt_max;                   /* largest vertex weight */
double maxdeg;                  /* largest weighted vertex degree */
{
  extern Heap   *heap;      /* pointer to heap of multigrid */
  extern double *MEM_OK;    /* variable for memory overflow exeception */
  extern int DEBUG_TRACE;       /* debug flag for Kernighan-Lin */
  extern int DEBUG_KL;          /* debug flag for Kernighan-Lin */
  extern double kl_total_time;
  extern double kl_init_time;
  extern double nway_kl_time;
  struct bilist ****buckets;    /* space for bucket sorts */
  struct bilist **listspace;    /* space for all bidirectional elements */
  int **dvals;                  /* change in penalty for each possible move */
  int **tops;                   /* starting dval for each type of move */
  double time, time1;           /* timing variables */
  int maxhop;                   /* maximum hops between sets */
  int maxdval;                  /* largest transition cost for a vertex */
  int i, j;                     /* loop counters */
  double seconds();
  void initialize(), count(), nway_kl();

  time = seconds();

  if (DEBUG_TRACE > 0) {
    {char buf[150]; sprintf(buf,"Entering klspiff, nvtxs = %d\n", nvtxs);UserWrite(buf);}
  }

# ifdef __DEBUG__
  /* debug_graph*/
  {
    extern int DEBUG_GRAPH;

    /* debug_graph */
    if (DEBUG_GRAPH>0)
    {
      int flag;
      char buf[100];
      flag = check(graph,nvtxs,-1);
      if (flag == FALSE)
      {
        sprintf(buf,"FATAL: check_graph returned FALSE in klspiff before kl()\n");
        UserWrite(buf);
        sprintf(buf,"nvtxs=%d,nsets=%d\n",nvtxs,nsets);
        UserWrite(buf);
      }
      else
      {
        sprintf(buf,"OK: check_graph = TRUE in klspiff() before kl,nvtxs=%d,nsets=%d\n",nvtxs,nsets);
        UserWrite(buf);
      }
    }
  }
# endif


  /* Find the largest hop value. */
  maxhop = 0;
  for (i=0; i<nsets; i++) {
    for (j=0; j<nsets; j++) {
      if (hops[i][j] > maxhop) maxhop = hops[i][j];
    }
  }

  maxdval = maxhop*maxdeg;

  /* Allocate and initialize a bunch of space for KL. */
  time1 = seconds();
  initialize(&buckets, &listspace, &dvals, &tops, nvtxs, nsets, maxdval);
  if (!MEM_OK) return;
  kl_init_time += seconds() - time1;

  if (DEBUG_KL > 1) {
    {char buf[150]; sprintf(buf," Before KL: ");UserWrite(buf);}
    count(graph, nvtxs, sets, nsets, hops, FALSE);
  }

  time1 = seconds();
  nway_kl(graph, nvtxs, buckets, listspace, tops, dvals, sets,
          maxdval, nsets, goal, hops, vwgt_max);
  if (!MEM_OK) return;
  nway_kl_time += seconds() - time1;

  if (DEBUG_KL > 1) {
    {char buf[150]; sprintf(buf," After KL:");UserWrite(buf);}
    count(graph, nvtxs, sets, nsets, hops, FALSE);
  }

  free_kl(buckets, listspace, dvals, tops);

# ifdef __DEBUG__
  /* debug_graph */
  {
    extern int DEBUG_GRAPH;

    /* debug_graph */
    if (DEBUG_GRAPH>0)
    {
      int flag;
      char buf[100];

      flag = check(graph,nvtxs,-1);
      if (flag == FALSE)
      {
        sprintf(buf,"FATAL: check_graph returned with FALSE after kl\n");
        UserWrite(buf);
        sprintf(buf,"nvtxs=%d,nsets=%d\n",nvtxs,nsets);
        UserWrite(buf);
      }
      else
      {
        sprintf(buf,"OK: check_graph = TRUE in klspiff() after kl,nvtxs=%d,nsets=%d\n",nvtxs,nsets);
        UserWrite(buf);
      }

    }
  }
# endif


  kl_total_time += seconds() - time;
}


static void free_kl(buckets, listspace, dvals, tops)
/* Free everything malloc'd for KL. */
struct bilist ****buckets;      /* space for bucket sorts */
struct bilist **listspace;      /* space for all bidirectional elements */
int **dvals;                    /* change in penalty for each possible move */
int **tops;                     /* starting dval for each type of move */
{
  extern Heap   *heap;      /* pointer to heap of multigrid */



  sfree((char *) dvals);
  sfree((char *) tops);

  sfree((char *) listspace[0]);
  sfree((char *) buckets[0][1]);
  sfree((char *) listspace);
  sfree((char *) buckets);
}
