// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/structs.h"
#include        "../main/defs.h"


/* Find the median of set of values. */
/* Can also find nested medians of several sets of values */
/* Routine works by repeatedly guessing a value, and discarding those */
/* values which are on the wrong side of the guess. */

void median(graph, vals, nvtxs, active, goal, using_vwgts, sets)
struct vtx_data **graph;        /* data structure with vertex weights */
double *vals;                   /* values of which to find median */
int nvtxs;                      /* number of values I own */
int *active;                    /* space for list of nvtxs ints */
double *goal;                   /* desired sizes for sets */
int using_vwgts;                /* are vertex weights being used? */
short *sets;                    /* set each vertex gets assigned to */
{
  extern int DEBUG_ASSIGN;  /* debug flag for assignment */
  double *vptr;                 /* loops through vals array */
  double val;                   /* value in vals array */
  double maxval;                /* largest active value */
  double max;               /* largest value */
  double minval;                /* smallest active value */
  double min;               /* smallest value */
  double guess;                 /* approximate median value */
  double nearup;                /* lowest guy above guess */
  double neardown;              /* highest guy below guess */
  double whigh;                 /* total weight of values above maxval */
  double wlow;                  /* total weight of values below minval */
  double wabove;                /* total weight of active values above guess */
  double wbelow;                /* total weight of active values below guess */
  double lweight;               /* desired weight of lower values in set */
  double uweight;               /* desired weight of upper values in set */
  int flags;                    /* logical check for each set */
  double frac;                  /* fraction of values I want less than guess */
  int *aptr;                    /* loops through active array */
  int *aptr2;                   /* helps update active array */
  int myactive;                 /* number of active values I own */
  double wfree;                 /* weight of vtxs not yet divided */
  int removed;                  /* number of my values eliminated */
  int npass = 0;                /* counts passes required to find median */
  int done;                     /* check for termination criteria */
  int vtx;                      /* vertex being considered */
  int i;                        /* loop counters */
  void median_assign();

  /* Initialize. */
  if (DEBUG_ASSIGN > 0)
  {char buf[150]; sprintf(buf,"   Entering median: nvtxs=%d goal[0]=%f goal[1]=%f \n",nvtxs,goal[0],goal[1]);UserWrite(buf);}

  /* Determine the desired weight sums for the two different sets. */
  lweight = goal[0];
  uweight = goal[1];

  myactive = nvtxs;
  whigh = wlow = 0;

  /* Find largest and smallest values in vector, and construct active list. */
  vptr = vals;
  aptr = active;
  minval = maxval = *(++vptr);
  *aptr++ = 1;
  for (i=2; i<=nvtxs; i++) {
    *aptr++ = i;
    val = *(++vptr);
    if (val > maxval) maxval = val;
    if (val < minval) minval = val;
  }

  /* remember extrem values */
  max = maxval;
  min = minval;

  done = FALSE;
  /* Extreme partitioning necessary? */
  if (goal[0] == 0)
  {
    guess = min-0.1;
    wlow = 0;
    done = TRUE;
  }
  if (goal[1] == 0)
  {
    guess = max+0.1;
    for (vtx=1; vtx<=nvtxs; vtx++)
    {
      if (using_vwgts) wlow += graph[vtx]->vwgt;
      else wlow++;
    }
    done = TRUE;
  }
  /* Loop until all sets are partitioned correctly. */
  while (!done) {
    npass++;

    /* Select a potential dividing value. */
    /* Currently, this assumes a linear distribution. */
    wfree = lweight + uweight - (wlow + whigh);
    frac = (lweight - wlow)/wfree;
    if (frac < .5) frac += .1;
    if (frac > .5) frac -= .1;
    guess = minval + frac*(maxval - minval);

    /* Now count the guys above and below this guess. */
    /* Also find nearest values on either side of guess. */
    wabove = wbelow = 0;
    nearup = maxval;
    neardown = minval;

    aptr = active;
    for (i=0; i<myactive; i++) {
      vtx = *aptr++;
      val = vals[vtx];
      if (val > guess) {
        if (using_vwgts) wabove += graph[vtx]->vwgt;
        else wabove++;
        if (val < nearup) nearup = val;
      }
      else if (val < guess) {
        if (using_vwgts) wbelow += graph[vtx]->vwgt;
        else wbelow++;
        if (val > neardown) neardown = val;
      }
      else {
        nearup = guess;
        neardown = guess;
      }
    }

    /* Select a half to discard. */
    done = TRUE;
    if (wlow + wbelow > lweight) {              /* Discard upper set. */
      flags = 1;
      whigh += wabove;
      maxval = neardown;
      done = FALSE;
    }
    else if (whigh+wabove > uweight) {          /* Discard lower set. */
      flags = -1;
      wlow += wbelow;
      minval = nearup;
      done = FALSE;
    }
    else {                                      /* Perfect partition! */
      flags = 0;
      wlow += wbelow;
      whigh += wabove;
      minval = nearup;
      maxval = neardown;
    }


    /* Remove discarded vertices from active list. */
    if (!done) {
      removed = 0;
      aptr = aptr2 = active;
      if (flags == 1) {
        for (i=0; i<myactive; i++) {
          if (vals[*aptr] <= neardown) *aptr2++ = *aptr;
          else ++removed;
          aptr++;
        }
      }
      else if (flags == -1) {
        for (i=0; i<myactive; i++) {
          if (vals[*aptr] >= nearup) *aptr2++ = *aptr;
          else ++removed;
          aptr++;
        }
      }
      myactive -= removed;
    }
    else myactive = 0;

    /* Check for alternate termination criteria. */
    if (!done && (maxval == minval))
    {
      done = TRUE;
      if (DEBUG_ASSIGN > 0)
      {
        {char buf[150]; sprintf(buf,"   alternate termination: maxval=%f minval=%f\n",maxval,minval);UserWrite(buf);}
      }
    }
  }
  if (DEBUG_ASSIGN > 0)
  {char buf[150]; sprintf(buf,"   entering medianassign: npass=%d guess=%f \n", npass, guess);UserWrite(buf);}

  median_assign(graph, vals, nvtxs, goal, using_vwgts, sets, wlow, guess);
}


void median_assign(graph, vals, nvtxs, goal, using_vwgts, sets,
                   wlow, guess)
struct vtx_data **graph;        /* data structure with vertex weights */
double *vals;                   /* values of which to find median */
int nvtxs;                      /* number of values I own */
double *goal;                   /* desired sizes for sets */
int using_vwgts;                /* are vertex weights being used? */
short *sets;                    /* assigned set for each vertex */
double wlow;                    /* sum of weights below guess */
double guess;                   /* median value */
{
  int i;                        /* loop counter */

  for (i=1; i<=nvtxs; i++) {
    if (vals[i] < guess) sets[i] = 0;
    else if (vals[i] > guess) sets[i] = 1;
    else {
      if (wlow < goal[0]) {
        sets[i] = 0;
        if (using_vwgts) wlow += graph[i]->vwgt;
        else wlow++;
      }
      else sets[i] = 1;
    }
  }
}
