/****************************************************************************/
/*																			*/
/* File:	  NeXTSurface.c													*/
/*																			*/
/* Purpose: 																*/
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

#import <app/appkit.h>

/* standard C includes */
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* dev interface */
#include "devices.h"

/* interface includes */
#include "compiler.h"
#include "general.h"


/* mif includes */
#include "MacMain.m"
#include "MacSurface.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define AboutMargin 	0	/* margin in pixels seperating pict and window	*/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* menu handles */
static MenuHandle myMenus[menuCount];	/* menuCount defined in MacGui.m	*/

static short currentCurs = 128; 		/* resource id of current cursor	*/

static Rect dragRect;					/* rect where windows can be moved	*/


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

Rect *DragRect ()
{
	return(&dragRect);
}

void SetMyCursor (short id)
{
	CursHandle cursor;
	
	cursor = GetCursor(id);
	if (cursor==NULL) return;
	SetCursor(*cursor);
	currentCurs = id;
}

void MousePosition (INT *ScreenPoint)
{
	Point pt;
	
	GetMouse(&pt);
	ScreenPoint[0] = (INT) pt.h;
	ScreenPoint[1] = (INT) pt.v;
	return;
}

INT MouseStillDown(void)
{
	/*if (StillDown())
		return (1);
	else
	{
		char buffer[8];
		
		UserIn(buffer);
		return (0);
	}*/
	return (StillDown());
}


/****************************************************************************/
/*																			*/
/* Function:  ScheduleCommand												*/
/*																			*/
/* Purpose:   pass control to the command selected from a menu				*/
/*																			*/
/* Input:	  short theMenu: menu number									*/
/*			  short theItem: item number									*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

void ScheduleCommand (short theMenu,short theItem)
{	
	switch (theMenu)
	{
		case appleID:
			switch (theItem)
			{
				case aboutCommand:	About(); break;
			}
			break;
	
		default:
			break;
	}
	return;
}

/****************************************************************************/
/*																			*/
/* Function:  DoCommand 													*/
/*																			*/
/* Purpose:   process command selection 									*/
/*																			*/
/* Input:	  long: menu and item information								*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

void DoCommand (long mResult)
{
	short temp,theItem,theMenu;
	Str255 name;
	
	theItem = LoWrd(mResult);
	theMenu = HiWrd(mResult);
	
	switch (theMenu)
	{
	
		case appleID :
			if (theItem==1)
				About();
			else
			{
				GetItem(myMenus[appleM],theItem,name);
				temp = OpenDeskAcc(name);
			}
			break;
		default:
			break;
	}
	HiliteMenu(0);
}

/****************************************************************************/
/*																			*/
/* Function:  SetUpMenus													*/
/*																			*/
/* Purpose:   read menus from resource file 								*/
/*																			*/
/* Input:	  none															*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

void SetUpMenus ()
{
	int i;
	
	for (i=0; i<menuCount; i++)
	{
		myMenus[i] = GetMenu(i+appleID);
		assert(myMenus[i]!=NULL);
		AddResMenu(myMenus[i],'DRVR');
	}
	for (i=0; i<menuCount; i++)
		InsertMenu(myMenus[i],0);				/* Add menus to menu bar */
	
	DrawMenuBar();
}

/****************************************************************************/
/*																		
   GetScreenSize - Return the available size of the screen												

   SYNOPSIS:
   INT GetScreenSize (INT size[2]);

   PARAMETERS:
.  size[2] - pointer to size vector

   DESCRIPTION:
   This function returns the available size of the screen.	
   
   RETURN VALUE:
   INT
.n     1 if screen available										
.n     0 if not.													
*/																			

/****************************************************************************/

INT GetScreenSize (INT size[2])
{
	size[0] = SCREEN_WIDTH;
	size[1] = SCREEN_HEIGHT-MENU_BAR;
	
	return (1);
}

/****************************************************************************/
/*																			
   InitMacGuiCommands - Load mac specific commands into environment											

   SYNOPSIS:
   int InitMacSurface (void)

   PARAMETERS:
.  void
 
   DESCRIPTION:
   This function loads mac specific commands into environment,
   makes Macintosh menus.

   RETURN VALUE:
   int
   0, if all is o.k.
												
   1, if an error occurred										
*/																			
/****************************************************************************/


int InitMacSurface (void)
{	
	/* reactangle where windows can be moved */
	SetRect(&dragRect,	-1,
						MENU_BAR-SCROLL_BAR-1,
						SCREEN_WIDTH  + SCROLL_BAR,
						SCREEN_HEIGHT + SCROLL_BAR);
	
	return(0);
}

