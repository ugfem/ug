// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugtimer.h														*/
/*																			*/
/* Purpose:   implements a simple timing facility							*/
/*																			*/
/* Author:	  Stefan Lang                                                   */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   980805   begin                                                */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UGTIMER__
#define __UGTIMER__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAX_TIMER       30

#define NEW_TIMER(n)    { new_timer(&(n)); }
#define DEL_TIMER(n)    { ug_timer[(n)].used = 0; }
#define START_TIMER(n)  { ug_timer[(n)].start = CURRENT_TIME; }
#define STOP_TIMER(n)   { ug_timer[(n)].stop = CURRENT_TIME; }
#define DIFF_TIMER(n)   (ug_timer[(n)].stop-ug_timer[(n)].start)
#define SUM_TIMER(n)    { STOP_TIMER(n) ug_timer[(n)].sum += DIFF_TIMER(n); }
#define EVAL_TIMER(n)   (ug_timer[(n)].sum)


/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef struct {
  char used;
  double start;
  double stop;
  double sum;
} UG_TIMER;


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

extern UG_TIMER ug_timer[MAX_TIMER];


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/


#endif
