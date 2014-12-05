// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*

   include relevant headers and define two macros to access the time

 */

#ifndef UG_UGTIME_H
#define UG_UGTIME_H

#ifndef UGLIB
#error Internal UG-lib header, must not be used in applications!
#endif

/* stolen from autoconf-docs */
#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

/* current time as DOUBLE value
   CURRENT_TIME should be the most accurate time (usually in micro seconds)
   CURRENT_TIME_LONG should be a time which measures some days without overflow
 */

/* !!! test for the time-functions */
#define CURRENT_TIME   (((DOUBLE)clock())/((DOUBLE)CLOCKS_PER_SEC))
#define CURRENT_TIME_LONG       CURRENT_TIME

/* !! test for difftime */
#define ARCH_DIFF_TIMER(x,y) (difftime((time_t)(x),(time_t)(y)))

#endif
