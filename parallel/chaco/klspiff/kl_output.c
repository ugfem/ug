// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/structs.h"
#include        "../main/defs.h"

static void p1bucket();

void pbuckets(buckets, listspace, maxdeg, nsets)
struct bilist ****buckets;      /* pointers to bucket lists */
struct bilist **listspace;      /* elements within buckets */
int maxdeg;                     /* maximum degree of a vertex */
int nsets;                      /* number of sets being divided into */
{
  struct bilist *lptr;          /* points to correct listspace */
  int i, j;                     /* loop counter */

  {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
  for (i=0; i<nsets; i++) {
    for (j=0; j<nsets; j++) {
      if (i!=j) {
        {char buf[150]; sprintf(buf,"For transition %d -> %d\n", i, j);UserWrite(buf);}
        if (j > i) lptr = listspace[j-1];
        else lptr = listspace[j];
        p1bucket(buckets[i][j], lptr, maxdeg);
        {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
      }
    }
  }
  {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
}


static void p1bucket(bucket, lptr, maxdeg)
struct bilist **bucket;         /* buckets holding bucket list */
struct bilist *lptr;            /* elements within bucket */
int maxdeg;                     /* maximum degree of a vertex */
{
  struct bilist *bptr;          /* loops through list at a bucket */
  int val;                      /* element in a bucket */
  int i;                        /* loop counter */

  for (i=2*maxdeg; i>=0; i--) {
    if (bucket[i] != NULL) {
      {char buf[150]; sprintf(buf,"  Bucket %d:", i-maxdeg);UserWrite(buf);}
      for (bptr=bucket[i]; bptr!=NULL; bptr=bptr->next) {
        val = ((int) bptr - (int) lptr)/sizeof(struct bilist);
        {char buf[150]; sprintf(buf," %d", val);UserWrite(buf);}
      }
      {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
    }
  }
}
