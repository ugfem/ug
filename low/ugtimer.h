// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  ugtimer.h														*/
/*																			*/
/* Purpose:   implements a simple timing facility							*/
/*																			*/
/* Author:	  Stefan Lang                                                   */ 
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   980805   begin                                                */ 
/*																			*/
/* Remarks: 																*/
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

#define NEW_TIMER(n)	{ new_timer(&(n)); }
#define DEL_TIMER(n)	{ ug_timer[(n)].used = 0; }
#define START_TIMER(n)  { ug_timer[(n)].start = CURRENT_TIME; }
#define STOP_TIMER(n)	{ ug_timer[(n)].stop = CURRENT_TIME; }
#define DIFF_TIMER(n)	(ug_timer[(n)].stop-ug_timer[(n)].start}
#define SUM_TIMER(n)	{ STOP_TIMER(n) ug_timer[(n)].sum += DIFF_TIMER(n); }
#define EVAL_TIMER(n)	(ug_timer[(n)].sum)


/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/


#endif
