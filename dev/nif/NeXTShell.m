/****************************************************************************/
/*																			*/
/* File:	  NeXTShell.c													*/
/*																			*/
/* Purpose:   handle the shell window (text i/o)							*/
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

/* standard C library */
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* low module */
#include "misc.h"
#include "defaults.h"
#include "general.h"

/* interface includes */
#include "devices.h"

/* mif includes */
#include "MacMain.m"
#include "MacMain.h"
#include "MacSurface.h"
#include "MacShell.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define SHELLTITLE	"<<<<< "UG_VERSION" shell >>>>>"

#define TERMHMIN	200 		/* minimum horizontal size of term window	*/
#define TERMVMIN	100 		/* minimum vertical size of term window 	*/

#define HiWrd(aLong) (((aLong) >> 16) & 0xFFFF)
#define LoWrd(aLong) ((aLong) & 0xFFFF)

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static ShellWindow *shell;					/* my shell window				*/
static char inpline[INPUTBUFFERLEN];		/* input line buffer			*/
static int nchars=0;						/* number of charcters in input */

static int scrollback=100;		/* default number of lines in scrollback bu */
static int charsperline=256;	/* default line width for shell window		*/
static int fontNum=monaco;		/* default font in shell window (Monaco)	*/
static int fontSize=9;			/* default font size						*/

/**** the universal headers contain the following defines: ******************/

#ifndef USES_UNIVERSAL_HEADERS
#ifndef USESROUTINEDESCRIPTORS
	#if USESROUTINEDESCRIPTORS
	typedef UniversalProcPtr ControlActionUPP;
	
	#define CallControlActionProc(userRoutine, theControl, partCode)		\
			CallUniversalProc((UniversalProcPtr)(userRoutine), uppControlActionProcInfo, (theControl), (partCode))
	#define NewControlActionProc(userRoutine)		\
			(ControlActionUPP) NewRoutineDescriptor((ProcPtr)(userRoutine), uppControlActionProcInfo, GetCurrentISA())
	#else 
	/*typedef ControlActionProcPtr ControlActionUPP; this doesn't run on Peter's Mac */
	typedef ProcPtr ControlActionUPP; /* this one works */
	
	#define CallControlActionProc(userRoutine, theControl, partCode)		\
			(*(userRoutine))((theControl), (partCode))
	#define NewControlActionProc(userRoutine)		\
			(ControlActionUPP)(userRoutine)
	#endif
#endif
#endif


static ControlActionUPP MyShellScrollActionPtr;

/****************************************************************************/
/*																			*/
/* Functions: GrowShellWindow												*/
/*																			*/
/* Purpose:   grow shell window 											*/
/*																			*/
/* Input:	  void *Data: the event 										*/
/*																			*/
/* Output:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT GrowShellWindow (EventRecord *theEvent)
{
	long growResult;
	Rect r,sizeRect,beforeRect,afterRect,viewRect;
	WindowPtr theWindow;

	/* get event */
	theWindow = MAC_WIN(shell);
	SetPort(theWindow);

	/* grow window */
	beforeRect = theWindow->portRect;
	SetRect(&sizeRect,TERMHMIN,TERMVMIN,SCREEN_WIDTH,SCREEN_HEIGHT);
	growResult = GrowWindow(theWindow,theEvent->where,&sizeRect);
	if (growResult==0)
		return(NO_POS_CHANGE);

	/* actually change windows size */
	SetPort(theWindow);
	SizeWindow(theWindow,LoWrd(growResult),HiWrd(growResult),true);
	afterRect = theWindow->portRect;

	/* make the new regions invalid */
	if (afterRect.right>beforeRect.right)
	{
		SetRect(&r,beforeRect.right-16,afterRect.top,afterRect.right,afterRect.bottom);
		InvalRect(&r);
	}
	if (afterRect.bottom>beforeRect.bottom)
	{
		SetRect(&r,afterRect.left,beforeRect.bottom-16,afterRect.right,afterRect.bottom);
		InvalRect(&r);
	}

	/* redraw resize */
	DrawGrowIcon(theWindow);
	r = theWindow->portRect;
	SetRect(&r,r.right-15,r.bottom-15,r.right,r.bottom);
	ValidRect(&r);
		
	/* move & resize scrollbars */
	r = theWindow->portRect;
	MoveControl(shell->vScrollBar,r.right-15,0);
	SizeControl(shell->vScrollBar,16,r.bottom-14);
	SetRect(&r,r.right-15,0,r.right,r.bottom-14);
	ValidRect(&r);
	r = theWindow->portRect;
	MoveControl(shell->hScrollBar,0,r.bottom-15);
	SizeControl(shell->hScrollBar,r.right-14,16);
	SetRect(&r,0,r.bottom-15,r.right-14,r.bottom);
	ValidRect(&r);
	
	/* update text view */
	r = theWindow->portRect;
	SetRect(&viewRect,0,0,r.right-16,r.bottom-16);
	shell->lines = (viewRect.bottom-viewRect.top-2)/shell->lineHeight;
	shell->cols = (viewRect.right-viewRect.left-4)/shell->charWidth;
	(**(shell->textH)).viewRect = viewRect;
	
	return (POS_CHANGE);
}

/****************************************************************************/
/*																			*/
/* Functions: DragShellWin													*/
/*																			*/
/* Purpose:   drag shell window 											*/
/*																			*/
/* Input:	  void *Data: the event 										*/
/*																			*/
/* return:	  POS_CHANGE													*/
/*			  NO_POS_CHANGE													*/
/*																			*/
/****************************************************************************/

INT DragShellWin (EventRecord *theEvent)
{
	WindowPtr theWindow;
	
	/* set pointers */
	theWindow = MAC_WIN(shell);
	
	/* drag window */
	DragWindow (theWindow,theEvent->where,DragRect());
	
	return(POS_CHANGE);
}

/****************************************************************************/
/*																			*/
/* Functions: MyShellScrollAction											*/
/*																			*/
/* Purpose:   do scroll action												*/
/*																			*/
/* Input:	  ControlHandle, partcode										*/
/*																			*/
/* Output:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static pascal void MyShellScrollAction (ControlHandle c, short p)
{
	short new;

	/* outside original partcode ? */
	if (p==0) return;
	
	if (c==shell->vScrollBar)
	{
		switch (p)
		{
			case inUpButton :
				new = MAX(GetCtlMin(c),shell->line-1);
				break;
			case inDownButton :
				new = MIN(GetCtlMax(c),shell->line+1);
				break;
			case inPageUp :
				new = MAX(GetCtlMin(c),shell->line-(shell->lines+1));
				break;
			case inPageDown :
				new = MIN(GetCtlMax(c),shell->line+(shell->lines+1));
				break;
		}
		if (shell->line!=new)
		{
			TEScroll(0,shell->lineHeight*(shell->line-new),shell->textH);
			shell->line = new;
		}
		SetCtlValue(c,new);
		DrawControls(MAC_WIN(shell));
	}
	if (c==shell->hScrollBar)
	{
		switch (p)
		{
			case inUpButton :
				new = MAX(GetCtlMin(c),shell->col-1);
				break;
			case inDownButton :
				new = MIN(GetCtlMax(c),shell->col+1);
				break;
			case inPageUp :
				new = MAX(GetCtlMin(c),shell->col-(shell->cols+1));
				break;
			case inPageDown :
				new = MIN(GetCtlMax(c),shell->col+(shell->cols+1));
				break;
		}
		if (shell->col!=new)
		{
			TEScroll(shell->charWidth*(shell->col-new),0,shell->textH);
			shell->col = new;
		}
		SetCtlValue(c,new);
		DrawControls(MAC_WIN(shell));
	}
}

/****************************************************************************/
/*																			*/
/* Functions: ShellWinContentClick											*/
/*																			*/
/* Purpose:   handle click into shell window								*/
/*																			*/
/* Input:	  void *Data: the event 										*/
/*																			*/
/* Output:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

void ShellWinContentClick (WindowPtr theWindow,EventRecord *theEvent)
{
	ControlHandle whichControl;
	short part1,part2;
	short new;

	SetPort(theWindow);
	GlobalToLocal(&(theEvent->where));
		
	/* check if a control has been hit */
	part1 = FindControl(theEvent->where,theWindow,&whichControl);
	switch (part1)
	{
		case 0:
			/* cut and paste mechanism here */
			return;
			
		case inThumb:
			if (TrackControl(whichControl,theEvent->where,(ControlActionUPP)-1)!=0)
			{
				new = GetCtlValue(whichControl);
				if ((whichControl==shell->vScrollBar)&&(shell->line!=new))
				{
					TEScroll(0,shell->lineHeight*(shell->line-new),shell->textH);
					shell->line = new;
				}
				if ((whichControl==shell->hScrollBar)&&(shell->col!=new))
				{
					TEScroll(shell->charWidth*(shell->col-new),0,shell->textH);
					shell->col = new;
				}
			}
			break;
		
		case inUpButton :
		case inDownButton :
		case inPageUp :
		case inPageDown :
			part2 = TrackControl(whichControl,theEvent->where,MyShellScrollActionPtr);
			break;
	}
	return;
}


/****************************************************************************/
/*																			*/
/* Functions: ActivateShellWin												*/
/*																			*/
/* Purpose:   process activate event for shell window						*/
/*																			*/
/* Input:	  WIN *w: pointer to application window data structure			*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

INT ActivateShellWin ()
{
	WindowPtr theWindow;
	
	theWindow = MAC_WIN(shell);
	SetPort(theWindow);
	DrawGrowIcon(theWindow);
	TEActivate(shell->textH);
	ShowControl(shell->vScrollBar);
	ShowControl(shell->hScrollBar);
	SetMyCursor(textCurs);
	
	return (0);
}

/****************************************************************************/
/*																			*/
/* Functions: DeactivateShellWin											*/
/*																			*/
/* Purpose:   process deactivate event for shell window 					*/
/*																			*/
/* Input:	  WIN *w: pointer to application window data structure			*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

INT DeactivateShellWin ()
{
	SetPort(MAC_WIN(shell));
	TEDeactivate(shell->textH);
	HideControl(shell->vScrollBar);
	HideControl(shell->hScrollBar);
	DrawGrowIcon(MAC_WIN(shell));
	
	return (0);
}

/****************************************************************************/
/*																			*/
/* Functions: UpdateShellWin												*/
/*																			*/
/* Purpose:   do update event for shell window								*/
/*																			*/
/* Input:	  WIN *w: pointer to application window data structure			*/
/*																			*/
/* Output:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT UpdateShellWin (void)
{
	WindowPtr theWindow;
	
	theWindow = MAC_WIN(shell);
	
	SetPort(theWindow);
	BeginUpdate(theWindow);
	EraseRgn(((GrafPtr)(theWindow))->visRgn);	  /* only what is necessary */
	DrawGrowIcon(theWindow);
	DrawControls(theWindow);
	TEUpdate(&(theWindow->portRect),shell->textH);
	EndUpdate(theWindow);
	
	return (0);
}



/****************************************************************************/
/*																			*/
/* Functions: IdleShellWindow												*/
/*																			*/
/* Purpose:   periodic task for shell window								*/
/*																			*/
/* Input:	  WIN *w: pointer to application window data structure			*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

void IdleShellWindow ()
{
	TEIdle(shell->textH);
}

/****************************************************************************/
/*																			*/
/* Functions: InsertInShellWin												*/
/*																			*/
/* Purpose:   insert a single character into shell window					*/
/*																			*/
/* Input:	  TERM_WIN *shell: window to insert 							*/
/*			  char c: character to insert									*/
/*																			*/
/* Output:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT InsertInShellWin (char c)
{
	int nLinesBefore,nLinesAfter,linesToDelete,linesToScroll;
	
	nLinesBefore = (**(shell->textH)).nLines;

	/* insert text into text edit record */
	TESetSelect(32760,32760,shell->textH);
	TEKey(c,shell->textH);
	
	/* delete lines past scroll back buffer */
	nLinesAfter = (**(shell->textH)).nLines;
	if (nLinesAfter>scrollback)
	{
		linesToDelete = nLinesAfter-scrollback;
		TESetSelect((long)(**(shell->textH)).lineStarts[0],
					(long)(**(shell->textH)).lineStarts[linesToDelete],shell->textH);
		TEDelete(shell->textH);
		TESetSelect(32760,32760,shell->textH);
		TECalText(shell->textH);
		nLinesAfter -= linesToDelete;
	}
	
	/* adjust view */
	if (nLinesAfter>shell->line+shell->lines)
	{
		linesToScroll = nLinesAfter-(shell->line+shell->lines);
		TEScroll(0,-shell->lineHeight*linesToScroll,shell->textH);
		TECalText(shell->textH);
		shell->line += linesToScroll;
	}
	
	/* adjust scroll bars */
	if (nLinesAfter!=nLinesBefore)
	{
		SetCtlMax(shell->vScrollBar,(short)nLinesAfter);
		SetCtlValue(shell->vScrollBar,(short)shell->line);
	}
	
	return (0);
}
	
/****************************************************************************/
/*																			*/
/* Functions: MacWriteString												*/
/*																			*/
/* Purpose:   write a string to a shell window								*/
/*																			*/
/* Input:	  TERM_WIN *shell: window to write to							*/
/*			  char *s: string to write										*/
/*																			*/
/* Output:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

void MacWriteString (char *s)
{
	int nLinesBefore,nLinesAfter,linesToDelete,linesToScroll;
	GrafPtr savePort;

#	if GENERATING68K 

	/*
		NB: this is not necessary for CW C-compiler version 1.2 if 'MPW newlines' is checked
			in the 'Language' preferences
	*/
	
	char *cp;
	
	/* replace \n by \r for Metrowerks CodeWarrior C compiler */
	for (cp=s; *cp!='\0'; cp++)
		if (*cp=='\n')
			*cp = '\r';
#	endif
	
	nLinesBefore = (**(shell->textH)).nLines;

	/* save previous window pointer */
	GetPort(&savePort);

	/* insert text into text edit record */
	SetPort(MAC_WIN(shell));
	TESetSelect(32760,32760,shell->textH);
	TEInsert((Ptr) s,(long) strlen(s),shell->textH);
	
	/* delete lines past scroll back buffer */
	nLinesAfter = (**(shell->textH)).nLines;
	if (nLinesAfter>scrollback)
	{
		linesToDelete = nLinesAfter-scrollback;
		TESetSelect((long)(**(shell->textH)).lineStarts[0],
					(long)(**(shell->textH)).lineStarts[linesToDelete],shell->textH);
		TEDelete(shell->textH);
		TESetSelect(32760,32760,shell->textH);
		TECalText(shell->textH);
		nLinesAfter -= linesToDelete;
	}
	
	/* adjust view */
	if (nLinesAfter>shell->line+shell->lines)
	{
		linesToScroll = nLinesAfter-(shell->line+shell->lines);
		TEScroll(0,-shell->lineHeight*linesToScroll,shell->textH);
		TECalText(shell->textH);
		shell->line += linesToScroll;
	}
	
	/* adjust scroll bars */
	if (nLinesAfter!=nLinesBefore)
	{
		SetCtlMax(shell->vScrollBar,(short)nLinesAfter);
		SetCtlValue(shell->vScrollBar,(short)shell->line);
	}
	
	/* reset saved window pointer */
	SetPort(savePort);
	return;
}

/****************************************************************************/
/*																			*/
/* Function:   ShellHandleKeybordEvent										*/
/*																			*/
/* Purpose:    handle keypress event for shell window						*/
/*																			*/
/* Input:																	*/
/*																			*/
/* Output:	   pointer to buffer if \n										*/
/*			   else NULL													*/
/*																			*/
/****************************************************************************/

char *ShellHandleKeybordEvent (char key)
{
	switch (key)
	{
		case '\n': /* carriage return */
		case '\r': /* carriage return */
			InsertInShellWin(key);
			inpline[nchars] = 0; /* zero terminated string */
			nchars = 0;
			return (inpline);
		
		case '\b': /* backspace */
			if (nchars>0)
			{
				InsertInShellWin(key);
				nchars--;
			}
			return (NULL);
		
		case '\f':	/* page down */
			MyShellScrollAction(shell->vScrollBar,inPageDown);
			return (NULL);
		
		case '\v':	/* page up */
			MyShellScrollAction(shell->vScrollBar,inPageUp);
			return (NULL);
		
		default:
			if ((key>=' ')&&(key<='~')||(key=='\t'))
			{
				InsertInShellWin(key);
				inpline[nchars++] = key;
			}
			return (NULL);;
	}
}

/****************************************************************************/
/*																			*/
/* Function:   ShellInitAndOpen 											*/
/*																			*/
/* Purpose:    Open a shell window and initialize a ShellWindow data str.	*/
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*																			*/
/* Output:	   0: OK														*/
/*			   1: error, could not complete 								*/
/*																			*/
/****************************************************************************/

int ShellInitAndOpen (ShellWindow *sh)
{
	Rect r,destRect,viewRect;
	GrafPtr myPort;
	FontInfo info;
	char buffer[256];
	int TermWinH=0,TermWinV=0,TermWinDH=400,TermWinDV=300;
	short maxX,maxY;
	char ShellTitle[256];
	
	shell = sh;
	
	strcpy(ShellTitle,SHELLTITLE);
	
	/* init the Macintosh toolbox */
	InitGraf(&qd.thePort);
	InitFonts();
	FlushEvents(everyEvent,0);
	InitWindows();
	TEInit();
	InitDialogs(nil);
	InitCursor();
	InitMenus();

	/* show about box */
	About();
	
	/* draw menu bar */
	SetUpMenus();
	
	/* create routine descriptor for running in mixed mode on Power PC
	   (CAUTION: usage of the new universal headers required!, s.a.) */
	MyShellScrollActionPtr = NewControlActionProc(MyShellScrollAction);
	
	/* read default values for text edit */
	if (GetDefaultValue(DEFAULTSFILENAME,"scrollback",buffer)==0)
		sscanf(buffer," %d ",&scrollback);
	if (GetDefaultValue(DEFAULTSFILENAME,"charsperline",buffer)==0)
		sscanf(buffer," %d ",&charsperline);
	if (GetDefaultValue(DEFAULTSFILENAME,"size",buffer)==0)
		sscanf(buffer," %d ",&fontSize);
	if (GetDefaultValue(DEFAULTSFILENAME,"font",buffer)==0)
	{
		if (strcmp(buffer,"Monaco")==0) 	fontNum = monaco;
		if (strcmp(buffer,"Courier")==0)	fontNum = courier;
		if (strcmp(buffer,"Times")==0)		fontNum = times;
		if (strcmp(buffer,"New York")==0)	fontNum = newYork;
		if (strcmp(buffer,"Helvetica")==0)	fontNum = helvetica;
		if (strcmp(buffer,"Geneva")==0) 	fontNum = geneva;
		if (strcmp(buffer,"Athens")==0) 	fontNum = athens;
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

	maxX = SCREEN_WIDTH-10;
	maxY = SCREEN_HEIGHT-40;
	
	TermWinH  = MIN(maxX-200,MAX(0,TermWinH));
	TermWinV  = MIN(maxY-150,MAX(0,TermW