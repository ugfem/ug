// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  sif.c	                                                                                                        */
/*																			*/
/* Purpose:   standard interface for ug                                                                         */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*																			*/
/* History:   22 Sep 95 begin, ug31                                                                             */
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

/* standard C includes */
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* interface includes */
#include "compiler.h"
#include "devices.h"
#include "initdev.h"
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

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* export global variables													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* export global variables per function call								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*
   GetScreenSize - Return sreen size

   SYNOPSIS:
   INT GetScreenSize (INT size[2]);

   PARAMETERS:
   .  size[2] - pointer to the size of screen

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

INT GetScreenSize (INT size[2])
{
  size[0] = 0;
  size[1] = 0;

  return (0);
}

/****************************************************************************/
/*
   GUI_GetNextEvent - Process an event from the system and pass it to ug

   SYNOPSIS:
   INT GetNextUGEvent (EVENT *theEvent, INT EventMask);

   PARAMETERS:
   .  theEvent - pointer to ug event
   .  EventMask - determines which events are reported

   DESCRIPTION:
   This function processes an event from the system and passes it to ug if necessary.
   Here we only return the string event.

   RETURN VALUE:
   INT
   .n     0 if no event occurred (ug or system)
   .n     1 if an event occurred (ug or system).
 */
/****************************************************************************/

INT GetNextUGEvent (EVENT *theEvent, INT EventMask)
{
  char *s;
  int cmdKey, onlyCmdKey;

  /* no event as default */
  theEvent->Type = NO_EVENT;
  theEvent->NoEvent.InterfaceEvent = 0;

  /* we have no command keys */
  if (EventMask==TERM_CMDKEY) return(0);

  /* read in a string from the user and store it in event structure */
  gets(theEvent->TermString.String);
  theEvent->Type = TERM_STRING;

  /* ready */
  return(0);
}


/****************************************************************************/
/*
   InitScreen - Init rest of GUI and return pointer to screen outputdevice

   SYNOPSIS:
   OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error);

   PARAMETERS:
   .  argcp - pointer to argument counter
   .  argv  - argument vector
   .  error - errorcode

   DESCRIPTION:
   Simply return NULL, since we do not have graphics capability.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      POINTER if all is o.k.
   .n      NULL if an error occurred.

 */
/****************************************************************************/

OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error)
{
  *error=0;
  return(NULL);
}


void ExitScreen (void)
{}


/****************************************************************************/
/*
   WriteString - write a string to a terminal window

   SYNOPSIS:
   void WriteString (const char *s);

   PARAMETERS:
   .  s - string to write

   DESCRIPTION:
   This function writes a string to a terminal window. Simply map it
   to printf().

   RETURN VALUE:
   void
 */
/****************************************************************************/

void WriteString (const char *s)
{
  printf("%s",s);
  return;
}

/****************************************************************************/
/*
   MousePosition - Get current mouse position


   SYNOPSIS:
   void MousePosition (INT *ScreenPoint);

   PARAMETERS:
   .  ScreenPoint - return result in this vector

   DESCRIPTION:
   This function gets current mouse position.


   RETURN VALUE:
   void
 */
/****************************************************************************/

void MousePosition (INT *ScreenPoint)
{
  return;
}


/****************************************************************************/
/*
   MouseStillDown - Determine if mouse button is still pressed


   SYNOPSIS:
   INT MouseStillDown (void);

   PARAMETERS:
   no parameters

   DESCRIPTION:
   This function returns true (1) if the mouse button is still pressed.
   The function should only be called after a button pressed event has been
   reported.


   RETURN VALUE:
   INT
   .n 0 if mouse button has been released
   .n 1 if mouse button is still pressed
 */
/****************************************************************************/

INT MouseStillDown (void)
{
  return (0);
}

void DrawInfoBox (WINDOWID win, const char *info)
{
  return;
}

INT WhichTool (WINDOWID win, const INT mouse[2], INT *tool)
{
  return (0);
}
