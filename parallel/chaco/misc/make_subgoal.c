// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include "../main/defs.h"

extern int nsets_glob;   /* total number of sets to partition into */

/* Use same ratios as original goals, but adjust based on set sizes. */
void make_subgoal(goal, subgoal, ndims, nsets, set, sub_vwgt_sum)
double *goal;                   /* goals for sets */
double *subgoal;                /* goals for subset of sets */
int ndims;                      /* number of bits to set at once */
int nsets;                      /* length of goal array */
int set;                        /* which set am I in? */
double sub_vwgt_sum;            /* sum of subgraph vertex weights */
{
  double tweight;               /* total weight among all subgoals */
  double ratio;                 /* scaling factor */
  int step;                     /* steps size in walking through array */
  int i, j;                     /* loop counter */
  double threshold;     /* threshold for minimal subgoal sizes */

  step = 1<<ndims;
  threshold = 1e-4;

  tweight = 0;
  j = 0;
  for (i=set; i<nsets; i+=step, j++) {
    subgoal[j] = goal[i];
    tweight += goal[i];
  }

  /* test for tweight with value 0 */
  if (tweight<threshold)
  {
    if (sub_vwgt_sum<threshold)
    {
      ratio = 1;
    }
    else
    {
      {char buf[150]; sprintf(buf,"make_subgoal failed: tweight=%f, sub_vwgt_sum=%f\n",tweight, sub_vwgt_sum);UserWrite(buf);}
      {char buf[150]; sprintf(buf,"Sorry, don't know how to partition subgraph.\n");UserWrite(buf);}
      exit(1);
    }
  }
  else ratio = sub_vwgt_sum/tweight;

  for (j=0; j<nsets/step; j++) {
    subgoal[j] *= ratio;
  }
}


void make_subgoal_med(goal, subgoal, ndims, nsets, set, sub_vwgt_sum,
                      part_type)
double *goal;                   /* goals for sets */
double *subgoal;                /* goals for subset of sets */
int ndims;                      /* number of bits to set at once */
int nsets;                      /* length of goal array */
int set;                        /* which set am I in? */
double sub_vwgt_sum;            /* sum of subgraph vertex weights */
int part_type;
{
  extern int ARCH_GOAL;
  double tweight;               /* total weight among all subgoals */
  double ratio;                 /* scaling factor */
  int step;                     /* steps size in walking through array */
  int i, j;                     /* loop counter */
  double threshold;     /* threshold for minimal subgoal sizes */

  step = 1<<ndims;
  threshold = 1e-4;

  tweight = 0;
  j = 0;

  if (ARCH_GOAL==0)
  {
    if (part_type==QUAD_1)
    {
      if (set==0)
      {
        subgoal[0] = goal[0];
        subgoal[1] = goal[2];
        tweight    = goal[0]+goal[2];
      }
      else                    /* set==1 */
      {
        subgoal[0] = goal[1];
        subgoal[1] = goal[3];
        tweight    = goal[1]+goal[3];
      }
    }
    else if (part_type==QUAD_2 || part_type==QUAD_3 || part_type==QUAD_4)
    {
      if (set==0)
      {
        subgoal[0] = goal[0];
        subgoal[1] = goal[1];
        tweight    = goal[0]+goal[1];
      }
      else                    /* set==1 */
      {
        subgoal[0] = goal[2];
        subgoal[1] = goal[3];
        tweight    = goal[2]+goal[3];
      }
    }
  }
  else
  {
    for (i=set; i<nsets; i+=step, j++)
    {
      subgoal[j] = goal[i];
      tweight += goal[i];
    }
  }

  /* test for tweight with value 0 */
  if (tweight<threshold)
  {
    if (sub_vwgt_sum<threshold)
    {
      ratio = 1;
    }
    else
    {
      {char buf[150]; sprintf(buf,"make_subgoal failed: tweight=%f, sub_vwgt_sum=%f\n",tweight, sub_vwgt_sum);UserWrite(buf);}
      {char buf[150]; sprintf(buf,"Sorry, don't know how to partition subgraph.\n");UserWrite(buf);}
      exit(1);
    }
  }
  else ratio = sub_vwgt_sum/tweight;

  for (j=0; j<nsets/step; j++) {
    subgoal[j] *= ratio;
  }
}


/* this version of make_subgoal is for use in recurse, when partitioning */
/* in n xm sets.                                                         */
/* Use same ratios as original goals, but adjust based on set sizes.     */
void make_subgoal_rec(goal, subgoal, ndims, ndims_tot, set, sub_vwgt_sum, assigned,
                      xlen, ylen, xlens, ylens, part_type)
double *goal;                   /* goals for sets */
double *subgoal;                /* goals for subset of sets */
int ndims;                      /* number of bits to set at once */
int ndims_tot;                  /* bits set up to this time */
int set;                        /* which set am I in? */
double sub_vwgt_sum;            /* sum of subgraph vertex weights */
int assigned;       /* bits presently assigned in this path */
int xlen;
int ylen;
int *xlens;
int *ylens;
int part_type;
{
  extern int DEBUG_ARRAY;
  extern int ARCH_GOAL;
  double tweight;               /* total weight among all subgoals */
  double ratio;                 /* scaling factor */
  int step;                     /* steps size in walking through array */
  int i, j, l;                  /* loop counter */
  double threshold;     /* threshold for minimal subgoal sizes */
  char buffer[100];
  int map2set;

  step = 1<<ndims;
  threshold = 1e-4;

  tweight = 0;
  j = 0;
  if (ARCH_GOAL==0)
  {
    for (l=0; l<ylen; l++)
      for (i=0; i<xlen; i++)
      {
        map2set = map(i,l,ndims,part_type,xlens,ylens);
        if (map2set==set)
        {
          subgoal[j] = goal[l*xlen+i];
          tweight += goal[l*xlen+i];
          j++;
          if (DEBUG_ARRAY>0)
          {
            sprintf(buffer,"subgoal: set=%d j=%d proc=%d subgoal=%f\n",
                    set, j-1, l*xlen+i, subgoal[j-1]);
            UserWrite(buffer);
          }
        }
      }
  }
  else
  {
    for (i=set<<ndims_tot|assigned; i<nsets_glob; i+=step<<ndims_tot, j++)
    {
      subgoal[j] = goal[i>>ndims_tot];
      tweight += goal[i>>ndims_tot];
    }
  }

  /* test for tweight with value 0 */
  if (tweight<threshold)
  {
    if (sub_vwgt_sum<threshold)
    {
      ratio = 1;
    }
    else
    {
      {char buf[150]; sprintf(buf,"make_subgoal_rec failed: tweight=%f, sub_vwgt_sum=%f\n",tweight, sub_vwgt_sum);UserWrite(buf);}
      {char buf[150]; sprintf(buf,"Sorry, don't know how to partition subgraph.\n");UserWrite(buf);}
      exit(1);
    }
  }
  else ratio = sub_vwgt_sum/tweight;
  for (i=0; i<j; i++) {
    subgoal[i] *= ratio;

#     ifdef __DEBUG__
    if (DEBUG_ARRAY>0)
    {
      sprintf(buffer,"new subgoal: i=%d subgoal=%f sub_vwgt_sum=%f, tweight=%f\n",
              i, subgoal[i], sub_vwgt_sum, tweight);
      UserWrite(buffer);
    }
#     endif
  }
}
