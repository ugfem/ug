// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

void solistout(solist,n,ngood,j)
struct orthlink **solist;       /* vector of pntrs to orthlnks */
int n;                          /* length of vecs to orth. against */
int ngood;                      /* number of good vecs on list*/
int j;                          /* current number of Lanczos steps */
{
  int i;                                /* index */
  extern int DEBUG_EVECS;               /* debugging output level for eigen computations */

  /* to placate alint */
  n = n;

  for (i=1; i<=ngood; i++) {
    if ( (solist[i])->index <= (int)(j/2) ) {
      {char buf[150]; sprintf(buf,".");UserWrite(buf);}
    }
    else {
      {char buf[150]; sprintf(buf,"+");UserWrite(buf);}
    }
    /*  Really detailed output:
       {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
       {char buf[150]; sprintf(buf,"depth %d\n",(solist[i])->depth);UserWrite(buf);}
       {char buf[150]; sprintf(buf,"index %d\n",(solist[i])->index);UserWrite(buf);}
       {char buf[150]; sprintf(buf,"ritzval %g\n",(solist[i])->ritzval);UserWrite(buf);}
       {char buf[150]; sprintf(buf,"betaji %g\n",(solist[i])->betaji);UserWrite(buf);}
       {char buf[150]; sprintf(buf,"tau %g\n",(solist[i])->tau);UserWrite(buf);}
       {char buf[150]; sprintf(buf,"prevtau %g\n",(solist[i])->prevtau);UserWrite(buf);}
       vecout((solist[i])->vec,1,n,"vec");
     */
  }
  {char buf[150]; sprintf(buf,"%d\n",ngood);UserWrite(buf);}

  if (DEBUG_EVECS > 2) {
    {char buf[150]; sprintf(buf,"  actual indicies: ");UserWrite(buf);}
    for (i=1; i<=ngood; i++) {
      {char buf[150]; sprintf(buf," %2d",solist[i]->index);UserWrite(buf);}
    }
    {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
  }
}
