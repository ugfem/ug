// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      stat.h                                                        */
/*                                                                          */
/* Purpose:   DDD statistical evaluation, header file                       */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/*                                                                          */
/* History:   95/01/16 kb  begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD_STAT_H__
#define __DDD_STAT_H__



/****************************************************************************/
/*                                                                          */
/* compile time constants defining array sizes                              */
/*                                                                          */
/****************************************************************************/

#define STAT_NDURS        64
#define STAT_NITEMS       32


/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/* STATISTICS: global record for statistical data                           */
/****************************************************************************/

typedef struct
{
  double starttime[4];                /* temporary starttime timer 1-4      */
  double durations[STAT_NDURS];       /* miscellaneous time measurements    */
  int items[STAT_NITEMS];             /* miscellaneous count measurements   */
} STATISTICS;




/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/* from stat.c
      statistical data (compile with option -DStatistics) */

#ifdef Statistics
extern STATISTICS statistics;
#endif


/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/



#ifdef Statistics

/*** macros for "DDD with statistic evaluation" ***/

#       define STAT_RESETX(n)       statistics.starttime[(n)] = CurrentTime()
#       define STAT_RESET           STAT_RESETX(0)
#       define STAT_RESET1          STAT_RESETX(0)
#       define STAT_RESET2          STAT_RESETX(1)
#       define STAT_RESET3          STAT_RESETX(2)
#       define STAT_RESET4          STAT_RESETX(3)

#       define STAT_TIMERX(n,d)     statistics.durations[(d)] =  \
  CurrentTime()-statistics.starttime[(n)]
#       define STAT_TIMER(d)        STAT_TIMERX(0,d)
#       define STAT_TIMER1(d)       STAT_TIMERX(0,d)
#       define STAT_TIMER2(d)       STAT_TIMERX(1,d)
#       define STAT_TIMER3(d)       STAT_TIMERX(2,d)
#       define STAT_TIMER4(d)       STAT_TIMERX(3,d)

#       define STAT_INCTIMERX(n,d)  statistics.durations[(d)] +=  \
  CurrentTime()-statistics.starttime[(n)]
#       define STAT_INCTIMER(d)     STAT_INCTIMERX(0,d)
#       define STAT_INCTIMER1(d)    STAT_INCTIMERX(0,d)
#       define STAT_INCTIMER2(d)    STAT_INCTIMERX(1,d)
#       define STAT_INCTIMER3(d)    STAT_INCTIMERX(2,d)
#       define STAT_INCTIMER4(d)    STAT_INCTIMERX(3,d)

#       define STAT_SETTIMER(d,v)   statistics.durations[(d)] = (v)
#       define STAT_GETTIMER(d)     statistics.durations[(d)]


#       define STAT_SETCOUNT(d,v)   statistics.items[(d)] = (v)
#       define STAT_GETCOUNT(d)     statistics.items[(d)]


#       define STAT_ZEROTIMER           { int i; for (i=0; i<STAT_NDURS; i++)  \
                                            STAT_SETTIMER(i,0.0);}
#       define STAT_ZEROCOUNT           { int i; for (i=0; i<STAT_NITEMS; i++)  \
                                            STAT_SETCOUNT(i,0);}
#       define STAT_ZEROALL         STAT_ZEROTIMER; STAT_ZEROCOUNT


#else

/*** macros for "DDD without statistic evaluation" ***/

#       define STAT_RESET
#       define STAT_RESET0
#       define STAT_RESET1
#       define STAT_RESET2
#       define STAT_RESET3
#       define STAT_RESET4

#       define STAT_TIMER(d)
#       define STAT_TIMER1(d)
#       define STAT_TIMER2(d)
#       define STAT_TIMER3(d)
#       define STAT_TIMER4(d)

#       define STAT_INCTIMER(d)
#       define STAT_INCTIMER1(d)
#       define STAT_INCTIMER2(d)
#       define STAT_INCTIMER3(d)
#       define STAT_INCTIMER4(d)

#       define STAT_INCTIMER(d)
#       define STAT_SETTIMER(d,v)
#       define STAT_GETTIMER(d)    0.0
#       define STAT_SETCOUNT(d,v)
#       define STAT_GETCOUNT(d)    0

#       define STAT_ZEROTIMER
#       define STAT_ZEROCOUNT
#       define STAT_ZEROALL

#endif


#endif
