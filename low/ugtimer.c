// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugtimer.c														*/
/*																			*/
/* Purpose:   implements a simple timimg facility			                */
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdio.h>
#include <assert.h>

#include "general.h"

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


/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {
  char used;
  double start_timer;
  double stop_timer;
  double sum_timer;
} UG_TIMER;

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

UG_TIMER ug_timer[MAX_TIMER];


/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   new_timer - allocate a new timer

   SYNOPSIS:
   void new_timer (int *n)

   PARAMETERS:

   DESCRIPTION:
   allocate a new timer

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

void new_timer (int *n)
{
  int i;

  *n = -1;

  for (i=0; i<MAX_TIMER; i++)
    if (ug_timer[i].used == 0)
    {
      *n = i;
      ug_timer[i].used = 1;
      ug_timer[i].start_timer = 0.0;
      ug_timer[i].stop_timer = 0.0;
      ug_timer[i].sum_timer = 0.0;
      break;
    }
  if (*n == -1)
  {
    printf("NEW_TIMER(): couldn't allocate new timer!\n");
    fflush(stdout);
    assert(0);
  }
  return;
}
