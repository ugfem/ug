/****************************************************************************/
/*																			*/
/* File:	  NeXTMain.c 													*/
/*																			*/
/* Purpose:   NEXTSTEP graphical user interface for ug 3.0 					*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						    	*/
/*																			*/
/*	History:  September 12, 1997 begin										*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files 									*/
/*																			*/
/****************************************************************************/

/* NeXT specific includes */
#import <appkit/Application.h>

/* standard C includes */
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* interface includes */
#include "heaps.h"
#include "defaults.h"
#include "compiler.h"
#include "devices.h"
#include "initdev.h"
#include "general.h"

#include "UGshell.h"
#include "NeXTGraph.h"
#include "NeXTMain.h"

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


static HEAP *guiHeap=NULL;
static long guiHeapSize=32000;

/****************************************************************************/
/*																			*/
/* export global variables													*/
/*																			*/
/****************************************************************************/


id		theUGshell;			/* the shell object */


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
	NXSize screenSize;
	
	[NXApp getScreenSize:&screenSize];
	printf ("Screen: Width = %f, height = %f\n", screenSize.width, screenSize.height);

	size[0] = screenSize.width;
	size[1] = screenSize.height;
	
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
	char buffer[256];
	NXSize screenSize;
	OUTPUTDEVICE *outdev;
		
	[Application new];

	theUGshell = [UGshell alloc];
	[NXApp setDelegate:theUGshell];
	
	[NXApp getScreenSize:&screenSize];
	printf ("Width = %f, height = %f\n", screenSize.width, screenSize.height);

	if ( [theUGshell setUp] == 0 )	{
		printf("UGshell setUp method failed\n");
		*error = 1;
		[NXApp free];
		}
	
	/* create output device */
	if ((outdev=InitNeXTOutputDevice())==NULL) return(NULL);

	/* get gui heapsize */
	if (GetDefaultValue(DEFAULTSFILENAME,"guimemory",buffer)==0)
		sscanf(buffer," %d ",&guiHeapSize);

	/* allocate gui heap */
	if ((guiHeap=NewHeap(GENERAL_HEAP,guiHeapSize,malloc(guiHeapSize)))==NULL) return(NULL);
	
	return (outdev);

	*error=0;
	return(NULL);
}


void ExitScreen (void)
{
	[NXApp free];
}


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
	//printf("%s",s);
	[theUGshell appendToText:s];
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
