// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  xmain.c														*/
/*																			*/
/* Purpose:   main file for X11 interface									*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*																			*/
/* History:   15.02.94 begin, ug3.0                                                                             */
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

/* X11 includes */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/Xresource.h>
#include <X11/keysym.h>

#define USE_XAW     /* enables use of Athena Text Widget in libXaw, libXt */
/* undefine USE_XAW, if you want pure X (only libX11) */

#ifdef USE_XAW
/* Xt & Xaw includes */
#include <X11/Intrinsic.h>
#include <X11/Vendor.h>
#include <X11/StringDefs.h>
#include <X11/Xaw/XawInit.h>
#include <X11/Xaw/AsciiText.h>
#include <X11/Xaw/Text.h>
#endif /* USE_XAW */

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

/* Xif includes */
#include "xshell.h"
#include "xgraph.h"


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

static ShellWindow shell;                                       /* our only shell window		*/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*																			*/
/* export global variables													*/
/*																			*/
/****************************************************************************/

/* global data needed everywhere */
#ifdef USE_XAW
XtAppContext context;   /* application context */
Widget applShell,
       ugshell;
XawTextPosition CursorPos,
                CutBeginPos = 0;
#endif

Display *display;                                                       /* the display					*/
int screen_num;                                                         /* screen on display			*/
char *prog_name;                                                        /* my own name					*/
Screen *screen_ptr;                                             /* dont know for what			*/
unsigned int display_width;                             /* size of screen if needed     */
unsigned int display_height;
int if_argc;                                                            /* command line args			*/
char **if_argv;

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
  size[0] = display_width;
  size[1] = display_height;

  return (0);
}

/****************************************************************************/
/*										*/
/* Function:  GUI_GetNextEvent							*/
/*										*/
/* Purpose:   process an event from the system and pass it to ug if nec.	*/
/*										*/
/* Input:	  EVENT *theEvent	 pointer to ug event to be filled	*/
/*										*/
/* return: INT 0: no event occurred (ug or system)				*/
/*			   1: an event occurred (ug or system)			*/
/*										*/
/****************************************************************************/

static Bool callback (Display *d, XEvent *report, char *arg)
{
  return(True);
}

/****************************************************************************/
/*
   GUI_GetNextEvent - Process an event from the system and pass it to ug

   SYNOPSIS:
   INT GetNextUGEvent (EVENT *theEvent, INT EventMask);

   PARAMETERS:
   .  theEvent - pointer to ug event

   .  EventMask -

   DESCRIPTION:
   This function processes an event from the system and passes it to ug if necessary.

   RETURN VALUE:
   INT
   .n     0 if no event occurred (ug or system)
   .n     1 if an event occurred (ug or system).
 */
/****************************************************************************/

INT GetNextUGEvent (EVENT *theEvent, INT Eventmask)
{
  XEvent report;
  XWindowAttributes xwa;
  char *s;
  GraphWindow *gw;
  int where_x,where_y,tool;
  int x,y,w,h;
  int cmdKey, onlyCmdKey;
  int flag;

  /* no event as default */
  theEvent->Type = NO_EVENT;
  theEvent->NoEvent.InterfaceEvent = 0;

  /* event loop */
  switch (Eventmask)
  {
  case EVERY_EVENT :
    if (!XCheckIfEvent(display,&report,callback,s)) return(0);
    PRINTDEBUG(dev,1,("XCheckIfEvent(): matching event found\n"));
    onlyCmdKey = 0;
    break;
  case TERM_STRING :
    if (!XCheckMaskEvent(display,KeyPressMask,&report)) return(0);
    onlyCmdKey = 0;
    break;
  case TERM_CMDKEY :
    if (!XCheckMaskEvent(display,KeyPressMask,&report)) return(0);
    onlyCmdKey = 1;
    break;
  default :
    return (0);
  }

  switch (report.type)
  {
  case Expose :
    if (report.xexpose.window==shell.win)
    {
                                #ifndef USE_XAW
      ShellHandleExposeEvent(&shell,&report);
                                #endif
      theEvent->NoEvent.InterfaceEvent = 1;
      break;
    }
    if (report.xexpose.count!=0) return(0);
    gw = WhichGW(report.xexpose.window);
    if (gw==NULL) break;
    theEvent->Type = DOC_UPDATE;
    theEvent->DocUpdate.win = (WINDOWID) gw;
    break;

  case EnterNotify :
    gw = WhichGW(report.xcrossing.window);
    if (gw==NULL) break;
    theEvent->Type = DOC_ACTIVATE;
    theEvent->DocActivate.win = (WINDOWID) gw;
    SetCurrentGW(gw);
    /*
       IFDEBUG(dev,1)
       printf("reporting DOC_ACTIVATE for view %s\n",MY_VIEW(gw)->name);
       ENDDEBUG
     */
    break;

  case ConfigureNotify :
    if (report.xconfigure.window==shell.win)
    {
                                #ifdef USE_XAW
      flag=XtDispatchEvent(&report);
      IFDEBUG(dev,1)
      if (flag==FALSE)
        PRINTDEBUG(dev,1,("XtDispatchEvent(): NO handler for this event found\n"))
        else
          PRINTDEBUG(dev,1,("XtDispatchEvent(): handler for this event found\n"))
          ENDDEBUG
                                #else /* USE_XAW */
      ShellHandleResizeEvent(&shell,&report);
      theEvent->NoEvent.InterfaceEvent = 1;
                                #endif /* USE_XAW */
          break;
    }
    gw = WhichGW(report.xconfigure.window);
    if (gw==NULL) break;
    /* check if size of window changed */
    if (  (report.xconfigure.width!=gw->window_width)
          ||(report.xconfigure.height!=gw->window_height))
    {
      /* new size, origin may have changed also ! */
      /* CAUTION: windows new position could not be determined ! */
      x = report.xconfigure.x;
      y = report.xconfigure.y;
      w = report.xconfigure.width;
      h = report.xconfigure.height;
      theEvent->Type = DOC_GROW;
      theEvent->DocGrow.win = (WINDOWID) gw;
      /*theEvent->DocGrow.Global_LL[0] = x;
         theEvent->DocGrow.Global_LL[1] = display_height-y;
         theEvent->DocGrow.Global_UR[0] = x+w;
         theEvent->DocGrow.Global_UR[1] = display_height-(y+h);*/
      theEvent->DocGrow.Global_LL[0] = x;
      theEvent->DocGrow.Global_LL[1] = y+h;
      theEvent->DocGrow.Global_UR[0] = x+w;
      theEvent->DocGrow.Global_UR[1] = y;
      theEvent->DocGrow.Local_LL[0] = 0;
      theEvent->DocGrow.Local_LL[1] = h-CONTROLSIZE-2;
      theEvent->DocGrow.Local_UR[0] = w-1;
      theEvent->DocGrow.Local_UR[1] = 0;
      gw->window_x = x;
      gw->window_y = y;
      gw->window_width = w;
      gw->window_height = h;
      /*
         IFDEBUG(dev,1)
         printf("reporting DOC_GROW dxmin=%g dymin=%g dxmax=%g dymax=%g\n",
              theEvent->DocGrow.LowerLeft[0],theEvent->DocGrow.LowerLeft[1],
              theEvent->DocGrow.UpperRight[0],theEvent->DocGrow.UpperRight[1]);
         ENDDEBUG
       */
      break;
    }
    if (  (report.xconfigure.x!=gw->window_x)
          ||(report.xconfigure.y!=gw->window_y))
    {
      /* it's a drag event */
      x = report.xconfigure.x;
      y = report.xconfigure.y;
      w = gw->window_width;
      h = gw->window_height;
      theEvent->Type = DOC_DRAG;
      theEvent->DocDrag.win = (WINDOWID) gw;
      theEvent->DocDrag.Global_LL[0] = x;
      theEvent->DocDrag.Global_LL[1] = display_height-y;
      theEvent->DocDrag.Global_UR[0] = x+w;
      theEvent->DocDrag.Global_UR[1] = display_height-(y+h);
      gw->window_x = report.xconfigure.x;
      gw->window_y = report.xconfigure.y;
      /*
         IFDEBUG(dev,1)
         printf("reporting DOC_DRAG dxmin=%g dymin=%g dxmax=%g dymax=%g\n",
              theEvent->DocDrag.LowerLeft[0],theEvent->DocDrag.LowerLeft[1],
              theEvent->DocDrag.UpperRight[0],theEvent->DocDrag.UpperRight[1]);
         ENDDEBUG
       */
      break;
    }
    break;

  case ButtonPress :
                        #ifdef USE_XAW
    if (report.xbutton.window == shell.win)
    {
      if (report.xbutton.button == Button1 ||
          report.xbutton.button == Button3 )
      {
        CursorPos = XawTextGetInsertionPoint(shell.wid);
        XawTextDisplayCaret(shell.wid,FALSE);
        if (report.xbutton.button == Button3)
          XawTextSetInsertionPoint(shell.wid,CutBeginPos);
      }
      break;
    }
                        #endif /* USE_XAW */
    gw = WhichGW(report.xbutton.window);
    if (gw==NULL) break;
    where_x = report.xbutton.x;
    where_y = report.xbutton.y;
    if (WhichTool(gw,where_x,where_y,&tool))
    {
      theEvent->Type = DOC_CHANGETOOL;
      theEvent->DocChangeTool.win = (WINDOWID) gw;
      theEvent->DocChangeTool.Tool = tool;
      IFDEBUG(dev,1)
      printf("reporting DOC_CHANGETOOL tool=%d\n",tool);
      ENDDEBUG
      break;
    }
    if (DrawRegion(gw,where_x,where_y))
    {
      theEvent->Type = DOC_CONTENTCLICK;
      theEvent->DocContentClick.win = (WINDOWID) gw;
      theEvent->DocContentClick.MousePosition[0] = where_x;
      theEvent->DocContentClick.MousePosition[1] = where_y;
    }

    break;

                #ifdef USE_XAW
  case ButtonRelease :
    if (report.xbutton.window == shell.win)
    {
      if (report.xbutton.button == Button2)
      {
        char *cutbuffer;
        int cnt;

        cutbuffer = XFetchBytes(display,&cnt);
        if (cutbuffer == NULL)
        {
          IFDEBUG(dev,1)
          printf("cut buffer empty\n");
          ENDDEBUG
          break;
        }
        IFDEBUG(dev,1)
        printf("ButtonRelease with cnt=%d cutbuffer=%s\n",cnt,cutbuffer);
        ENDDEBUG

        AppendOrInsertCutbuffer(&shell,cutbuffer,cnt);
        XFree(cutbuffer);
      }
    }
    break;
                #endif /* USE_XAW */

  case KeyPress :
    s = ShellHandleKeybordEvent(&shell,&report,&cmdKey,onlyCmdKey);
    if (s==NULL)
    {
      theEvent->NoEvent.InterfaceEvent = 1;
      break;
    }
    if (cmdKey)
    {
      theEvent->Type = TERM_CMDKEY;
      theEvent->TermCmdKey.CmdKey = s[0];
      break;
    }
    theEvent->Type = TERM_STRING;
    strcpy(theEvent->TermString.String,s);
    break;

  default :
    break;
  }

        #ifdef USE_XAW
  /* Send all events to shell widget */
  flag=XtDispatchEvent(&report);
  IFDEBUG(dev,1)
  if (flag==FALSE)
    PRINTDEBUG(dev,1,("XtDispatchEvent(): NO handler for this event found\n"))
    else
      PRINTDEBUG(dev,1,("XtDispatchEvent(): Handler for this event found\n"))
      ENDDEBUG

      switch (report.type)
      {
      case ButtonRelease :
        if (report.xbutton.window == shell.win)
        {
          if (report.xbutton.button == Button1)
            CutBeginPos = XawTextGetInsertionPoint(shell.wid);

          if (report.xbutton.button == Button1 ||
              report.xbutton.button == Button3 )
          {
            XawTextSetInsertionPoint(shell.wid,CursorPos);
            XawTextDisplayCaret(shell.wid,TRUE);
          }
        }
        break;
      }
        #endif

  return(0);
}

/****************************************************************************/
/*
   InitScreen - Init rest of GUI and return pointer to screen outputdevice

   SYNOPSIS:
   OUTPUTDEVICE *InitScreen (int argc, char **argv, INT *error);

   PARAMETERS:
   .  argc - argument counter
   .  argv - argument vector
   .  error - errorcode

   DESCRIPTION:
   This function inits rest of GUI and return ptr to screen outputdevice.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      POINTER if all is o.k.
   .n      NULL if an error occurred.

 */
/****************************************************************************/

OUTPUTDEVICE *InitScreen (int argc, char **argv, INT *error)
{
  OUTPUTDEVICE *d;
  char buf[128];

  /* copy parameters to globals */
  if_argc = argc;
  if_argv = argv;

  /* connect procedure */
  prog_name = argv[0];
  strcpy(buf,SHELLWINNAME);
  argv[0] = buf;

        #ifdef USE_XAW
  applShell = XtAppInitialize (&context, "Xug3",
                               (XrmOptionDescRec*)NULL, 0,
                               &argc, argv,
                               (String*)NULL,
                               (Arg*)NULL, 0);

  display = XtDisplay (applShell);
        #else /* USE_XAW */
  display=XOpenDisplay(NULL);
        #endif /* USE_XAW */
  if (display == NULL)
  {
    fprintf(stderr,"%s: cannot connect to X server %s\n",prog_name,XDisplayName(NULL));
    exit(-1);
  }
  screen_num = DefaultScreen(display);
  screen_ptr = DefaultScreenOfDisplay(display);
  display_width = DisplayWidth(display,screen_num);
  display_height = DisplayHeight(display,screen_num);

  /* make output device first to have color table */
  d = InitXOutputDevice();
  if (d==NULL) {*error=1; return(NULL);}

  /* create shell window */
  ShellOpen(&shell);

  /* init cursors and tools */
  InitControls(shell.win);

  *error = 0;
  return(d);
}

/****************************************************************************/
/*
   WriteString - write a string to a terminal window

   SYNOPSIS:
   void WriteString (const char *s);

   PARAMETERS:
   .  s -

   DESCRIPTION:
   This function writes a string to a terminal window.

   RETURN VALUE:
   void
 */
/****************************************************************************/

void WriteString (const char *s)
{
        #ifdef STDIF
  fputs(s,stdout);
        #ifdef Debug
  fflush(stdout);
        #endif
  /* printf("%s",s); */
        #else
  ShellInsertString(&shell,s);
  return;
        #endif
}
