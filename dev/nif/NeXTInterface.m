/****************************************************************************/
/*																			*/
/* File:	  NeXTInterface.c												*/
/*																			*/
/* Purpose:   NEXTSTEP graphical user interface for ug 3.0 					*/
/*            All functions UG needs for its graphical output are			*/
/*            defined here.	See below for comments.							*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						    	*/
/*																			*/
/*	History:  January 16, 1998 begin										*/
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
#include "compiler.h"
#include "devices.h"
#include "initdev.h"
#include "debug.h"
#include "general.h"
#include "misc.h"
#include "heaps.h"
#include "defaults.h"
#include "cmdint.h"

#include "UGShellWindow.h"
#include "UGGraphWindow.h"
#include "NeXTInterface.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/
#define VERSION 	"This is "UG_VERSION" from $Date$\n"

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/



static int doneFlag=1;
static int nbuf=0;
static char buf[INPUTBUFFERLEN];

static HEAP *guiHeap=NULL;
static long guiHeapSize=32000;
SHORT_POINT moveto_point;
static GRAPH_WINDOW *currgw;
static GRAPH_WINDOW *windowList=NULL;
static OUTPUTDEVICE *NeXTOutputDevice;

struct ugcolortable	{
	int red;
	int green;
	int blue;
};
typedef struct ugcolortable UGColorTable;
static UGColorTable ColorTable[256];



UGShellWindow		*theUGshell;			/* the shell object */

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

/****************************************************************************
	The following functions are connected with UGs graphics system in 
	InitScreen.
*****************************************************************************/

static void NeXTMove (SHORT_POINT point)
{
	moveto_point.x = point.x;
	moveto_point.y = point.y;
}

static void NeXTDraw (SHORT_POINT point)
{
	[currgw->theView lockFocus];
	PSmoveto(moveto_point.x, moveto_point.y);
	PSlineto(point.x, point.y);
	PSstroke();
	[currgw->theView unlockFocus];
}

static void NeXTPolyline (SHORT_POINT *points, INT n)
{
	int i;
	
	if (n<2) return;

	[currgw->theView lockFocus];
	PSmoveto(points[0].x, points[0].y);
	for (i=1; i<n; i++) 
		PSlineto(points[i].x, points[i].y);
	PSstroke();
	[currgw->theView unlockFocus];
}

static void NeXTInversePolyline (SHORT_POINT *points, INT n)
{
	return;
}

static void NeXTPolygon (SHORT_POINT *points, INT n)
{
	int i;
	
	if (n<3) return;
	
	[currgw->theView lockFocus];
	PSmoveto(points[0].x, points[0].y);
	for (i=1; i<n; i++) PSlineto (points[i].x, points[i].y);
	PSclosepath();
	PSfill();
	PSstroke();
	[currgw->theView unlockFocus];
	return;
}

static void NeXTInversePolygon (SHORT_POINT *points, INT n)
{
	return;
}

static void NeXTErasePolygon (SHORT_POINT *points, INT n)
{
	return;
}

static void Marker (short n, short s, SHORT_POINT point)
{
	return;
}

		
static void NeXTPolymark (short n, SHORT_POINT *points)
{
	return;
}

static void InvMarker (short n, short s, SHORT_POINT point)
{
	return;
}

		
static void NeXTInvPolymark (short n, SHORT_POINT *points)
{
	return;
}

static void NeXTDrawText (const char *s, INT mode)
{
	return;
}

static void NeXTCenteredText (SHORT_POINT point, const char *s, INT mode)
{
	return;
}

static void NeXTClearViewPort (void)
{
	return;
}

static void NeXTSetLineWidth (short w)
{
	PSsetlinewidth((float)w);
	return;
}

static void NeXTSetTextSize (short s)
{
	id theFont;
    theFont = [Font newFont:"Ohlfs" size:(float)s];
	[[theUGshell shellText] setFont:theFont];
	return;
}

static void NeXTSetMarkerSize (short s)
{
	return;
}

static void NeXTSetMarker (short s)
{
	return;
}

static void NeXTSetColor (long index)
{	
	NXSetColor(NXConvertRGBToColor((float)ColorTable[index].red/(float)0xFFFF,
									(float)ColorTable[index].green/(float)0xFFFF,
									(float)ColorTable[index].blue/(float)0xFFFF));
	return;
}

static void NeXTSetPaletteEntry (long index, short r, short g, short b)
{
	return;
}

static void NeXTSetNewPalette (long start, long count, short *r, short *g, short *b)
{
	return;
}

static void NeXTGetPaletteEntry (long index, short *r, short *g, short *b)
{
	return;
}

static void NeXTFlush (void)
{
	[currgw->theView lockFocus];
	PSstroke();
	NXPing();
	[currgw->theWindow flushWindow];
	[currgw->theView unlockFocus];
	return;
}




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
	NXEvent *event;
	char *s;
	int cmdKey, onlyCmdKey;
	static int count = 0;
	int		NeXTEventMask;
	
	/* no event as default */
	theEvent->Type = NO_EVENT;
	theEvent->NoEvent.InterfaceEvent = 0;
	theEvent->NoEvent.GraphWinActive = (WINDOWID) 0;
return 0;
	/* event loop */
	switch (EventMask)
	{
		case EVERY_EVENT:
			NeXTEventMask = NX_ALLEVENTS;
			//printf("EVERY_EVENT passed to GetNextUGEvent()\n");
			break;
		case TERM_STRING:
			NeXTEventMask = NX_KEYDOWNMASK;
			//printf("TERM_STRING passed to GetNextUGEvent()\n");
			break;
		case TERM_CMDKEY:
			NeXTEventMask = NX_KEYDOWNMASK | NX_FLAGSCHANGEDMASK;
			//printf("TERM_CMDKEY passed to GetNextUGEvent()\n");
			break;
		default:
			printf("Error in GetNextUGEvent:\tswitch (EventMask=%d)\n",EventMask);
			return (EVENT_ERROR);
	}

	//printf("EventMask: %d,\tNeXTEventMask: %d,\t",EventMask,NeXTEventMask);

	// NX_MODALRESPTHRESHOLD ?? NX_BASETHRESHOLD
	event = [NXApp getNextEvent:NeXTEventMask waitFor:0.0 threshold:NX_BASETHRESHOLD];
	if (event == NULL) {
		//printf("event==NULL\n");
		return 0;
		}

	//printf("event type: %d\n",event->type);

	// Send the event to where it belongs
	[theUGshell sendEvent:event];
			
	switch (event->type)
	{
		case NX_LMOUSEDOWN:
			//printf("NX_LMOUSEDOWN (%d)\n", NX_LMOUSEDOWN);
			break;

		case NX_LMOUSEUP:
			//printf("NX_LMOUSEUP (%d)\n", NX_LMOUSEUP);
			break;

		case NX_KEYDOWN:
		case NX_FLAGSCHANGED:
			//printf("NX_KEYDOWN (%d) NX_FLAGSCHANGED (%d)\n",
			//		NX_KEYDOWN, NX_FLAGSCHANGED);
			/*s = ShellHandleKeybordEvent(event);
			if (s==NULL)
			{
				theEvent->NoEvent.InterfaceEvent = 1;
				break;
			}
			if (event->type==NX_FLAGSCHANGED)
			{
				//printf("command key was pressed\n");
				theEvent->Type = TERM_CMDKEY;
				theEvent->TermCmdKey.CmdKey = s[0];
				break;
			}*/
			theEvent->Type = TERM_STRING;
			strcpy(theEvent->TermString.String,s);
			break;
		
		case NX_MOUSEMOVED:
			//printf("NX_MOUSEMOVED (%d)\n");
			break;

		default:
			//printf("default %d\n",event->type);
			break;
	}

	return 0;
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
	printf("MouseStillDown ausgefuehrt\n");
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



static WINDOWID NeXT_OpenOutput (
	const char *title,						/* title of the window	 		*/
	INT x, INT y, INT width, INT height,	/* plot rgn in standard coord.	*/
	INT *Global_LL, INT *Global_UR, 		/* global machine coordinates	*/
	INT *Local_LL, INT *Local_UR,			/* local machine coordinates	*/
	INT *error) 							/* error code					*/
{
	id theWindow;
	NXRect graphicsRect;

	GRAPH_WINDOW	*gw;
	
	*error = 0;
	/* create GRAPH_WINDOW structure and put in list */
	gw = (GRAPH_WINDOW*) GetMem(guiHeap,sizeof(GRAPH_WINDOW),GENERAL_HEAP);
	if (gw==NULL)	{*error=1; return(0);}
	gw->next   = windowList;
	windowList = gw;
		
	NXSetRect(&graphicsRect, x, y, width, height);
	gw->theWindow = [ [UGGraphWindow alloc]
				initContent: &graphicsRect
				style:		NX_TITLEDSTYLE
				backing:	NX_BUFFERED
				buttonMask:	NX_MINIATURIZEBUTTONMASK
				defer:		NO];
	[gw->theWindow setTitle:title];
	[gw->theWindow display];
	
	NXSetRect(&graphicsRect, 0, 0, width, height);
	gw->theView = [[View alloc] initFrame:&graphicsRect];
	[gw->theWindow setContentView:gw->theView];
	//[gw->theView setBackgroundGray:NX_WHITE];
	
	/* fill global and local lower left and upper right in the devices coordinate system */
	gw->Global_LL[0] = Global_LL[0] = x;
	gw->Global_LL[1] = Global_LL[1] = y+height;
	gw->Global_UR[0] = Global_UR[0] = x+width;
	gw->Global_UR[1] = Global_UR[1] = y;
	
	gw->Local_LL[0] = Local_LL[0] = 0;
	gw->Local_LL[1] = Local_LL[1] = height;
	gw->Local_UR[0] = Local_UR[0] = width;
	gw->Local_UR[1] = Local_UR[1] = 0;
	
	return((WINDOWID)gw);
}


static INT NeXT_CloseOutput (WINDOWID win)
{
	[((GRAPH_WINDOW *)win)->theView unlockFocus];
	[((GRAPH_WINDOW *)win)->theWindow close];
}


INT NeXT_ActivateOutput (WINDOWID win)
{
	currgw = (GRAPH_WINDOW *)(win);

	[currgw->theWindow makeKeyAndOrderFront:nil];
	//[currgw->theView lockFocus];
	return 0;
}

INT NeXT_UpdateOutput (WINDOWID win, INT tool)
{
	return 0;
}



void MousePosition (INT *ScreenPoint)
{
	NXPoint mousePoint;
	
	[theUGshell getMouseLocation:&mousePoint];
	ScreenPoint[0] = (INT) mousePoint.x;
	ScreenPoint[1] = (INT) mousePoint.y;
	return;
}

char *ShellHandleKeybordEvent (NXEvent *theEvent)
{
	char string[2];
	string[1]='\0';
	
	string[0] = theEvent->data.key.charCode;
	[theUGshell appendToText:string];

	//printf("ShellHandleKeybordEvent char ->%s<- (%d,%d)\n",string,string[0],'\n');
	
	if ((int)string[0] == 13)	{
		//printf("return getippt\n");
		buf[nbuf] = '\0';
		nbuf = 0;
		return buf;
	}
	buf[nbuf++] = string[0];
	
	return NULL;
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
	NXRect	aRect;
	NXSize	screenSize;
	OUTPUTDEVICE *outdev;
	int i,j;
	unsigned short res,delta,max,r,g,b;
	int scrollback;
	int charsperline;
	int fontNum;
	int fontSize;
	int TermWinH=0,TermWinV=0,TermWinDH=400,TermWinDV=300;
	
	/*
	 * Get some default values
	 */
	if (GetDefaultValue(DEFAULTSFILENAME,"scrollback",buffer)==0)
		sscanf(buffer," %d ",&scrollback);
	if (GetDefaultValue(DEFAULTSFILENAME,"charsperline",buffer)==0)
		sscanf(buffer," %d ",&charsperline);
	if (GetDefaultValue(DEFAULTSFILENAME,"size",buffer)==0)
		sscanf(buffer," %d ",&fontSize);
	if (GetDefaultValue(DEFAULTSFILENAME,"font",buffer)==0)
	{
		if (strcmp(buffer,"Monaco")==0) 	fontNum = 0;
	}
	
	/* read default size of shell window */
	if (GetDefaultValue(DEFAULTSFILENAME,"TermWinH",buffer)==0)
		sscanf(buffer," %d ",&TermWinH);
	if (GetDefaultValue(DEFAULTSFILENAME,"TermWinV",buffer)==0)
		sscanf(buffer," %d ",&TermWinV);
	if (GetDefaultValue(DEFAULTSFILENAME,"TermWinDH",buffer)==0)
		sscanf(buffer," %d ",&TermWinDH);
	if (GetDefaultValue(DEFAULTSFILENAME,"TermWinDV",buffer)==0)
		sscanf(buffer," %d ",&TermWinDV);
	
		
	[Application new];

	/*
	 * Open up a Window
	 */
	NXSetRect(&aRect, 100.0, 350.0, 600.0, 300.0);
	theUGshell = [[UGShellWindow alloc] initContent:&aRect
					style:NX_RESIZEBARSTYLE
					backing:NX_BUFFERED
					buttonMask:NX_MINIATURIZEBUTTONMASK
					defer:NO];
	[theUGshell setTitle:"UG shell"];
	[theUGshell display];
	[NXApp setDelegate:theUGshell];
	
	if ( [theUGshell setUp] == 0 )	{
		printf("UGshell setUp method failed\n");
		*error = 1;
		[NXApp free];
		return NULL;
		}
	
	/* create output device */
	if ((NeXTOutputDevice=CreateOutputDevice("screen"))==NULL) return(NULL);

	/* init output device 'NeXT' */
	NeXTOutputDevice->OpenOutput	= NeXT_OpenOutput;
	NeXTOutputDevice->CloseOutput	= NeXT_CloseOutput;
	NeXTOutputDevice->ActivateOutput	= NeXT_ActivateOutput;
	NeXTOutputDevice->UpdateOutput		= NeXT_UpdateOutput;
	NeXTOutputDevice->v.locked		= 1;
	

	/*
	 * Now do a lot of setup stuff 
	 * (connect drawing functions, do color initialization, ...)
	 */
	
	res = 31;
	delta = 2048;
	max = 63488;

	/* fixed colors */
	i = 0;
	ColorTable[i].red = 0xFFFF;	ColorTable[i].green = 0xFFFF; ColorTable[i++].blue = 0xFFFF;
	ColorTable[i].red = 0xD000; ColorTable[i].green = 0xD000; ColorTable[i++].blue = 0xD000;
	ColorTable[i].red = 0xFFFF; ColorTable[i].green = 0xFFFF; ColorTable[i++].blue = 0x0   ;
	ColorTable[i].red = 0xFFFF; ColorTable[i].green = 0x0	; ColorTable[i++].blue = 0xFFFF;
	ColorTable[i].red = 0xFFFF; ColorTable[i].green = 0x0	; ColorTable[i++].blue = 0x0   ;
	ColorTable[i].red = 0x0;	ColorTable[i].green = 0xFFFF; ColorTable[i++].blue = 0xFFFF;
	ColorTable[i].red = 0x0;	ColorTable[i].green = 0xFFFF; ColorTable[i++].blue = 0x0   ;
	ColorTable[i].red = 0x0;	ColorTable[i].green = 0x0	; ColorTable[i++].blue = 0xFFFF;
	ColorTable[i].red = 0x0;	ColorTable[i].green = 0x0	; ColorTable[i++].blue = 0x0   ;
	ColorTable[i].red = 65520 ; ColorTable[i].green = 32240 ; ColorTable[i++].blue = 0x0   ;
	ColorTable[i].red = 65520 ; ColorTable[i].green = 60000 ; ColorTable[i++].blue = 0x0   ;

	/* color spectrum */
	/* TODO: This can be done much nicer on a NeXT */
	r = g = 0; b = max;
	ColorTable[i].red = r; ColorTable[i].green = g; ColorTable[i].blue = b; i++;

	/* blue to cyan */
	for (j=0; j<res; j++)
	{
		g += delta;
		ColorTable[i].red = r; ColorTable[i].green = g; ColorTable[i].blue = b; i++;
	}
	/* cyan to green */
	for (j=0; j<res; j++)
	{
		b -= delta;
		ColorTable[i].red = r; ColorTable[i].green = g; ColorTable[i].blue = b; i++;
	}
	/* green to yellow */
	for (j=0; j<res; j++)
	{
		r += delta;
		ColorTable[i].red = r; ColorTable[i].green = g; ColorTable[i].blue = b; i++;
	}
	/* yellow to red */
	for (j=0; j<res; j++)
	{
		g -= delta;
		ColorTable[i].red = r; ColorTable[i].green = g; ColorTable[i].blue = b; i++;
	}

	NeXTOutputDevice->black = 8;
	NeXTOutputDevice->gray = 1;
	NeXTOutputDevice->white = 0;
	NeXTOutputDevice->red = 4;
	NeXTOutputDevice->green = 6;
	NeXTOutputDevice->blue = 7;
	NeXTOutputDevice->cyan = 5;
	NeXTOutputDevice->yellow = 2;
	NeXTOutputDevice->darkyellow = 10;
	NeXTOutputDevice->magenta = 3;
	NeXTOutputDevice->orange = 9;
	NeXTOutputDevice->hasPalette = 1;
	NeXTOutputDevice->range = i;
	NeXTOutputDevice->spectrumStart = 11;
	NeXTOutputDevice->spectrumEnd = i-1;

	if (i>0) printf("Using color map with %d entries\n",i);

	NeXTOutputDevice->signx		 		= 1;
	NeXTOutputDevice->signy		 		= -1;
	
	/* init pointers to basic drawing functions */
	NeXTOutputDevice->Move				= NeXTMove;
	NeXTOutputDevice->Draw				= NeXTDraw;
	NeXTOutputDevice->Polyline			= NeXTPolyline;
	NeXTOutputDevice->InversePolyline	= NeXTInversePolyline;
	NeXTOutputDevice->Polygon			= NeXTPolygon;
	NeXTOutputDevice->InversePolygon 	= NeXTInversePolygon;
	NeXTOutputDevice->ErasePolygon		= NeXTErasePolygon;
	NeXTOutputDevice->Polymark			= NeXTPolymark;
	NeXTOutputDevice->InvPolymark		= NeXTInvPolymark;
	NeXTOutputDevice->DrawText			= NeXTDrawText;
	NeXTOutputDevice->CenteredText		= NeXTCenteredText;
	NeXTOutputDevice->ClearViewPort		= NeXTClearViewPort; 					
	
	/* init pointers to set functions */
	NeXTOutputDevice->SetLineWidth		= NeXTSetLineWidth;
	NeXTOutputDevice->SetTextSize		= NeXTSetTextSize;
	NeXTOutputDevice->SetMarker			= NeXTSetMarker;
	NeXTOutputDevice->SetMarkerSize		= NeXTSetMarkerSize;
	NeXTOutputDevice->SetColor			= NeXTSetColor;
	NeXTOutputDevice->SetPaletteEntry	= NeXTSetPaletteEntry;
	NeXTOutputDevice->SetNewPalette		= NeXTSetNewPalette;
	
	/* init pointers to miscellaneous functions */
	NeXTOutputDevice->GetPaletteEntry	= NeXTGetPaletteEntry;
	NeXTOutputDevice->Flush				= NeXTFlush;

	printf("output device 'screen' for NeXT window manager created\n");

	/* get gui heapsize */
	if (GetDefaultValue(DEFAULTSFILENAME,"guimemory",buffer)==0)
		sscanf(buffer," %d ",&guiHeapSize);

	/* allocate gui heap */
	if ((guiHeap=NewHeap(GENERAL_HEAP,guiHeapSize,malloc(guiHeapSize)))==NULL)
		return(NULL);
		
	return (NeXTOutputDevice);

	*error=0;
	return(NULL);
}


void ExitScreen (void)
{
	[theUGshell free];
	[NXApp free];
}

/****************************************************************************/
/*D
   NeXTCommandLoop - NEXTSTEP version of CommandLoop 

   SYNOPSIS:
   void CommandLoop (int argc, char **argv);

   PARAMETERS:
.  argc - argument counter
.  argv - argument vector

   DESCRIPTION:
   This command is used 

   RETURN VALUE:
   void
D*/
/****************************************************************************/

static void PrintVersionString (void)
{
	char ver[128];
	int  i,j,k;

	
	strcpy(ver,VERSION);
	for (i=0; i<100; i++)
  	{
		if (ver[i] == '\0') break;
		if (ver[i] == '$') break;
  	}
	k = 0;
	for (j=i+6; j<100; j++)
  	{
		if (ver[j] == '$') 
	  		k = 1;
		else
	  		ver[j-6-k] = ver[j];
		if (ver[j] == '\0') break;
  	}

   	/* print version */ 
	UserWrite(ver);
}

void NeXTCommandLoop (int argc, char **argv)
{
	INT error;
	int i,kerr;
	char c,errLine[256],spcLine[256],buffer[256];
	char *inpLine;
	char *strStart;
	int batch = FALSE;
	
	for (i=1; i<argc; i++)
		if (argv[i][0]!='-')
			batch = TRUE;

	/* alloc input line buffer */
	if ((inpLine=(char *)malloc(cmdintbufsize))==NULL)
	{
		PrintErrorMessage('F',"NeXTCommandLoop()","could not allocate inpLine buffer");
		return;
	}
	inpLine[0] = (char) 0;

	#ifdef ModelP
	if (me==master)
	{
		/* FOR MASTER PROCESSOR */
	#endif

		PrintVersionString();

		/* reset doneFlag */
		doneFlag=FALSE;


		/* if (argc==-1): second start of CommandLoop */
		if (argc != -1)
		{
			/* execute init script */
	    	if (GetDefaultValue(DEFAULTSFILENAME,"initscript",buffer)==0) {
 		    	strcpy(inpLine,"execute ");
				strcat(inpLine,buffer);
				error = InterpretCommand(inpLine);
				if (error==QUITCODE)
			    	doneFlag=TRUE;
			}
		}

		if (batch)
		{
        	i = 1;
        	while (i<argc && !doneFlag)
        	{
            	/* execute batch file */
            	if (argv[i][0]!='-')
            	{
                	sprintf(inpLine,"execute %s\n",argv[i]);

                	InterpretCommand(inpLine); /* execute command line argument */
					if (i + 1 < argc)
				    	if (strcmp(argv[i+1],"-noquit") == 0) {
					    	CommandLoop(-1,NULL);
							return;
						}
                	InterpretCommand("quit\n");/* end program */
            		i++;
                	continue;
            	}
            	/* set command from command line */
            	if ((argv[i][0]=='-')&&(argv[i][1]=='S'))
            	{
                	if (i+1<argc)
                	{
                    	sprintf(inpLine,"%s\n",(argv[i+1]));
                    	InterpretCommand(inpLine);
                    	i++;
                	}
                	else
                		UserWrite("Error in command line option -S\n");
            		i++;
                	continue;
            	}
            	/* logon command from command line */
            	if ((argv[i][0]=='-')&&(argv[i][1]=='L'))
            	{
                	if (i+1<argc)
                	{
                    	sprintf(inpLine,"logon %s\n",(argv[i+1]));
                    	InterpretCommand(inpLine);
                    	i++;
                	}
                	else
                		UserWrite("Error in command line option -L\n");
            		i++;
                	continue;
            	}
            	i++;
			}
		}
		
		// Print the prompt
		[theUGshell appendToText:PROMPT];
		[NXApp run];

	#ifdef ModelP
	}
	else
	{
		/* FOR PROCESSORS WITH ME!=MASTER */

		/* reset doneFlag */
		doneFlag = FALSE;


		while (!doneFlag)
		{
			if (doneFlag) break;
			error=ParExecCommand(inpLine);
			if (error==QUITCODE)
				doneFlag=TRUE;
		}
	}
	#endif


	/* call ExitUg() at the end of CommandLoop in order to avoid that
	   the application programmer will forget to call it at the end of
	   the application. */
	ExitUg();
}

