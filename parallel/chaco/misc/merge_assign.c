// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include "../main/defs.h"

/* Combine the old assignment value with the new partition. */
void merge_assignments(assignment, subassign, ndims, subnvtxs, loc2glob, part_type)
short *assignment;              /* assignment list for graph */
short *subassign;               /* subgraph assignment list */
int ndims;                      /* number of bits to shift over */
int subnvtxs;                   /* number of vtxs in subgraph */
int *loc2glob;                  /* subgraph -> graph numbering map */
{
  extern int ARCH_GOAL;
  int i;                        /* loop counters */
  char buffer[100];

  for (i=1; i<=subnvtxs; i++)
  {
    subassign++;
    assignment[loc2glob[i]] |=  (*subassign << ndims);
  }
}

/* Combine the old assignment value with the new partition. */
void merge_assignments_med(assignment, subassign, ndims, shift, subnvtxs, loc2glob, part_type)
short *assignment;              /* assignment list for graph */
short *subassign;               /* subgraph assignment list */
int ndims;                      /* number of bits to shift over */
int shift;
int subnvtxs;                   /* number of vtxs in subgraph */
int *loc2glob;                  /* subgraph -> graph numbering map */
{
  extern int ARCH_GOAL;
  int i;                        /* loop counters */
  char buffer[100];

  for (i=1; i<=subnvtxs; i++)
  {
    subassign++;
    if (ARCH_GOAL==0)
    {
      if (ndims==2)
      {
        if (part_type==QUAD_1)
          assignment[loc2glob[i]] = ((*subassign)<<shift)
                                    |assignment[loc2glob[i]];
        else if (part_type==QUAD_2 || part_type==QUAD_3 || part_type==QUAD_4)
          assignment[loc2glob[i]] = (assignment[loc2glob[i]]<<shift)
                                    |(*subassign);
      }
      else
      {
        sprintf(buffer,"FATAL: merge_assignments() failed, ndims=%d\n", ndims);
        UserWrite(buffer);
      }
    }
    else
      assignment[loc2glob[i]] |=  (*subassign << shift);
  }
}

void merge_assignments_rec(assignment, subassign, ndims, subnvtxs, loc2glob,
                           set, xlen, ylen, xlens, ylens, part_type, ready)
short *assignment;              /* assignment list for graph */
short *subassign;               /* subgraph assignment list */
int ndims;                      /* number of bits to shift over */
int subnvtxs;                   /* number of vtxs in subgraph */
int *loc2glob;                  /* subgraph -> graph numbering map */
int set;
int xlen;
int ylen;
int *xlens;
int *ylens;
int part_type;
short *ready;
{
  extern int DEBUG_ARRAY;
  extern int ARCH_GOAL;
  int i;                        /* loop counters */
  char buffer[100];

  for (i=1; i<=subnvtxs; i++)
  {
    subassign++;
    if (ARCH_GOAL==0)
    {
      if (ndims==1 || ndims==2)
      {

#               ifdef __DEBUG__
        if (DEBUG_ARRAY>0)
        {
          if (DEBUG_ARRAY>1)
          {
            if (assignment[loc2glob[i]]==set)
            {
              sprintf(buffer,"OK: set=%d assignment=%d loc2glob=%d i=%d subassign=%d in merge_assignments_rec()\n", set,
                      assignment[loc2glob[i]], loc2glob[i], i,
                      *subassign);
              UserWrite(buffer);
            }
          }
          if (assignment[loc2glob[i]]!=set)
          {
            sprintf(buffer,"FATAL: set=%d assignment=%d loc2glob=%d i=%d\n",
                    set, assignment[loc2glob[i]],
                    loc2glob[i], i);
            UserWrite(buffer);
          }
        }
#               endif
        assignment[loc2glob[i]] =  suba2a(assignment,loc2glob,
                                          i, subassign, ndims,
                                          xlen,ylen,set,xlens,
                                          ylens,part_type);
        ready[loc2glob[i]] = TRUE;
      }
      else
      {
        sprintf(buffer,"FATAL: merge_assignments_rec() failed, ndims=%d\n", ndims);
      }
    }
    else
      assignment[loc2glob[i]] |=  (*subassign << ndims);
  }
}
