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
#include <X11/Shell.h>
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
#include "xmain.h"

#ifdef __GUI__
#include "gui.h"
#endif


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define DBG_LEVEL               4

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static ShellWindow shell;                                       /* our only shell window		*/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* export global variables													*/
/*																			*/
/****************************************************************************/

/* global data needed everywhere */
#ifdef USE_XAW
XtAppContext context;   /* application context */
Widget toplevel,
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
#ifndef __GUI__
int user_interface = XUI;                                       /* user interface to open       */
#else
int user_interface = GUI;                                       /* user interface to open       */
#endif
int cui = 0;                                /* reset toggle for cui         */

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

INT GetNextUGEvent_CUI (EVENT *theEvent, INT EventMask)
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

INT GetNextUGEvent_XUI (EVENT *theEvent, INT Eventmask)
{
  XEvent report;
  XWindowAttributes xwa;
  Window root,child;
  INT root_x,root_y;
  char *s;
  GraphWindow *gw;
  int where_x,where_y,tool;
  INT pt[2];
  int x,y,w,h;
  int cmdKey, onlyCmdKey;
  int flag;
  unsigned int keys_buttons;
  static int count = 0;

  /* no event as default */
  theEvent->Type = NO_EVENT;
  theEvent->NoEvent.InterfaceEvent = 0;
  theEvent->NoEvent.GraphWinActive = 0;
  /*if (!XQueryPointer(display,report.xmotion.window,&root,&child,
                                     &root_x,&root_y,&where_x,&where_y,&keys_buttons))*/

  /* event loop */
  switch (Eventmask)
  {
  case EVERY_EVENT :
    if (!XCheckIfEvent(display,&report,callback,NULL)) return(0);
    PRINTDEBUG(dev,DBG_LEVEL,("XCheckIfEvent(): matching event found count=%d\n",count++));
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
       IFDEBUG(dev,DBG_LEVEL)
       printf("reporting DOC_ACTIVATE for view %s\n",MY_VIEW(gw)->name);
       ENDDEBUG
     */
    break;

  case ConfigureNotify :
    if (report.xconfigure.window==shell.win)
    {
                                #ifdef USE_XAW
      flag=XtDispatchEvent(&report);
      IFDEBUG(dev,DBG_LEVEL)
      if (flag==FALSE) {
        PRINTDEBUG(dev,DBG_LEVEL,("XtDispatchEvent(): NO handler for this event found\n"));
      }
      else {
        PRINTDEBUG(dev,DBG_LEVEL,("XtDispatchEvent(): handler for this event found\n"));
      }
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
         IFDEBUG(dev,DBG_LEVEL)
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
         IFDEBUG(dev,DBG_LEVEL)
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
    pt[0] = where_x = report.xbutton.x;
    pt[1] = where_y = report.xbutton.y;
    if (WhichTool((WINDOWID)gw,pt,&tool))
    {
      theEvent->Type = DOC_CHANGETOOL;
      theEvent->DocChangeTool.win = (WINDOWID) gw;
      theEvent->DocChangeTool.Tool = tool;
      theEvent->DocChangeTool.MousePosition[0] = where_x;
      theEvent->DocChangeTool.MousePosition[1] = where_y;
      IFDEBUG(dev,DBG_LEVEL)
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
          IFDEBUG(dev,DBG_LEVEL)
          printf("cut buffer empty\n");
          ENDDEBUG
          break;
        }
        IFDEBUG(dev,DBG_LEVEL)
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

  case MotionNotify :
    while (XCheckMaskEvent(display,PointerMotionMask,&report)) ;
    if (!XQueryPointer(display,report.xmotion.window,&root,&child,
                       &root_x,&root_y,&where_x,&where_y,&keys_buttons))
      break;

    gw = WhichGW(report.xmotion.window);
    if (gw!=NULL)
    {
      theEvent->NoEvent.GraphWinActive = (WINDOWID) gw;
      theEvent->NoEvent.Mouse[0] = where_x;
      theEvent->NoEvent.Mouse[1] = where_y;
    }
    break;

  default :
    break;
  }

        #ifdef USE_XAW
  /* Send all events to shell widget */
  flag=XtDispatchEvent(&report);
  IFDEBUG(dev,DBG_LEVEL)
  if (flag==FALSE) {
    PRINTDEBUG(dev,DBG_LEVEL,("XtDispatchEvent(): NO handler for this event found\n"));
  }
  else {
    PRINTDEBUG(dev,DBG_LEVEL,("XtDispatchEvent(): Handler for this event found\n"));
  }
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

INT GetNextUGEvent (EVENT *theEvent, INT EventMask)
{
  if (CUI_ON) GetNextUGEvent_CUI (theEvent,EventMask);
  else if (!NUI_ON) GetNextUGEvent_XUI (theEvent,EventMask);
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
   This function inits rest of GUI and return ptr to screen outputdevice.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      POINTER if all is o.k.
   .n      NULL if an error occurred.

 */
/****************************************************************************/

OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error)
{
  OUTPUTDEVICE *d;
  char buf[128];
  char buffer[128];
  int i,j;
        #ifdef USE_XAW
  int n;
  Arg args[20];
        #endif

  /* copy parameters to globals */
  if_argc = *argcp;
  if_argv = argv;

  /* connect procedure */
  prog_name = argv[0];
  strcpy(buf,SHELLWINNAME);
  argv[0] = buf;

  /* now set user interface */
  for (i=1; i<*argcp; i++)
    if (strcmp(argv[i],"-ui")==0)
    {
      int ok = 0;
      if (sscanf(argv[i+1],"%s",buffer)!=1)
      {
        fprintf(stderr,"%s: invalid use of option -ui [selected_ui]\n",prog_name);
        fprintf(stderr,"%s: choose for option -ui [" XUI_STRING "|" CUI_STRING "|"
                NUI_STRING "|" GUI_STRING "|" XGUI_STRING "|"
                CGUI_STRING "|" CNUI_STRING "]\n",prog_name);
        exit(-1);
      }

      if (strcmp(buffer,XUI_STRING)==0)       { user_interface = XUI; ok = 1; }
      if (strcmp(buffer,CUI_STRING)==0)       { user_interface = CUI; ok = 1; }
      if (strcmp(buffer,GUI_STRING)==0)       { user_interface = GUI; ok = 1; }
      if (strcmp(buffer,NUI_STRING)==0)       { user_interface = NUI; ok = 1; }
      if (strcmp(buffer,XGUI_STRING)==0)      { user_interface = XGUI; ok = 1; }
      if (strcmp(buffer,CGUI_STRING)==0)      { user_interface = CGUI; ok = 1; }
      if (strcmp(buffer,CNUI_STRING)==0)      { user_interface = CNUI; ok = 1; }

      if (!ok)
      {
        fprintf(stderr,"%s: invalid use of option -ui [selected_ui]\n",prog_name);
        fprintf(stderr,"%s: choose for option -ui [" XUI_STRING "|" CUI_STRING "|"
                NUI_STRING "|" GUI_STRING "|" XGUI_STRING "|"
                CGUI_STRING "|" CNUI_STRING "]\n",prog_name);
        exit(-1);
      }

      /* erase arguments from arglist */
      for (j=i+2; j<*argcp; j++) argv[j-2] = argv[j];
      i -= 1;
      *argcp -= 2;
    }

  if (NUI_ON) { *error = 0; return(NULL); }

        #ifdef USE_XAW
  /* set input focus to true, due to problems with some
     window manager implemtations, DEC, LINUX, SUN      */
  n=0;
  XtSetArg(args[n], XtNinput, TRUE); n++;

  toplevel = XtAppInitialize (&context, "Xug3",
                              (XrmOptionDescRec*)NULL, 0,
                              argcp, argv,
                              (String*)NULL,
                              args, n);

  display = XtDisplay(toplevel);
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

    #ifdef __GUI__
  /* init gui */
  if(GUI_Init(display,context))
  {
    printf("InitScreen(): GUI_Init() returned error");
    return(NULL);
  }
    #endif

  *error = 0;
  return(d);
}



void ExitScreen (void)
{}



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
  if (CUI_ON)
  {
    fputs(s,stdout);
                #ifdef Debug
    fflush(stdout);
                #endif
    /* printf("%s",s); */
  }
  else if (!NUI_ON)
  {
    ShellInsertString(&shell,s);
    return;
  }
}
