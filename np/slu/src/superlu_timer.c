// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*******************************************************/
/* This file was modified according to be integrated   */
/* in the ug-software-package. The changes are due to  */
/*                                                     */
/* Klaus Johannsen,                                    */
/* University of Heidelberg,                           */
/* Im Neuenheimer Feld 368,                            */
/* Germany.                                            */
/*******************************************************/

/*
 * Purpose
 * =======
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */


#ifdef SUN
/*
 *      It uses the system call gethrtime(3C), which is accurate to
 *	nanoseconds.
 */
#include <sys/time.h>

double SuperLU_timer_() {
  return ( (double)gethrtime() / 1e9 );
}

#else

#ifndef __MWCW__
#include <time.h>
#include <sys/times.h>
#endif

double SuperLU_timer_()
{
#ifndef __MWCW__
  struct tms use;
  double tmp;
  times(&use);
  tmp = use.tms_utime;
  tmp += use.tms_stime;
  return (double)tmp / (double)CLK_TCK;
#else
  /* do nothing */
  return 0.0;
#endif
}

#endif
