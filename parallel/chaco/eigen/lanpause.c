// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

/* Determine whether to pause in Lanczos */
int lanpause(j,lastpause,interval,q,n,pausemode,version)
int j;                  /* current step */
int lastpause;          /* when last paused */
int interval;           /* interval between pauses */
double **q;             /* the Lanczos vectors */
int n;                  /* length of Lanczos vectors */
int *pausemode;         /* which pausing criterion to use */
int version;            /* which version of sel. orth. we are using */
{
  extern int DEBUG_EVECS;               /* debugging level for eigen computation */
  double paige_dot;                     /* q[j]^T q[1] */
  double paigetol;                      /* pause if paigedot > paigetol */
  double dot();                         /* standard dot product */
  double fabs();                        /* intrinsic abs. value */
  void checkorth();


  /* Check orthogonality of last Lanczos vector against previous ones */
  if (DEBUG_EVECS > 3) {
    checkorth(q,n,j);
  }

  /* periodic reorthogonalization */
  if (version == 1 || version == 2) {
    if ((j - lastpause) == interval) {return(TRUE);}
    else                             {return(FALSE);}
  }

  /* Run until orthogonality with first Lanczos vector deteriorates, then switch
     switch to periodic reorthog. */
  if (version == 3) {
    paigetol = 1.0e-3;
    if (*pausemode == 1) {
      paige_dot = fabs(dot(q[1],1,n,q[j]));
      if (paige_dot > paigetol && j > 1) {
        if (DEBUG_EVECS > 1) {
          {char buf[150]; sprintf(buf,"  Pausing on step %3d with Paige prod. = %g\n",j,paige_dot);UserWrite(buf);}
        }
        *pausemode = 2;
        return(TRUE);
      }
      else {return(FALSE);}
    }
    if (*pausemode == 2) {
      if ((j - lastpause) == interval) {return(TRUE);}
      else                             {return(FALSE);}
    }
  }

  /* shouldn't ever get this far, but alint really wants a return value */
  return(FALSE);
}
