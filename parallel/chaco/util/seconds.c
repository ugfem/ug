// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include   <sys/time.h>
#include   "ppif.h"

/* Faster and more robust timer. */
/*
   double seconds()
   {
   double time;
   struct rusage rusage;
   int getrusage();

   getrusage(RUSAGE_SELF, &rusage);
   time = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
          1.0e-6*(rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));
   return(time);
   }
 */

/* The old timer, which wrapped around after ~36 minutes. */
double seconds()
{
  return((double) CurrentTime());
}
