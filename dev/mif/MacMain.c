// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MacMain.c                                                                                                     */
/*																			*/
/* Purpose:   Macintosh graphical user interface for ug 3.0                             */
/*																			*/
/* Author:	  Peter Bastian/Henrik Rentz-Reichert							*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							*/
/*					 ug@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*			  13.02.95 begin, ug version 3.0								*/
/*																			*/
/* Remarks:   former MacGui.c/.h											*/
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment mif
#endif

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* standard C includes */
/*  #include <string.h> */
/*  #include <math.h> */
/*  #include <stddef.h> */
/*  #include <stdlib.h> */
/*  #include <stdio.h> */
/*  #include <assert.h> */

/* mac toolbox includes */
#include <ToolUtils.h>
/*  #include <Types.h> */
/*  #include <Memory.h> */
/*  #include <Quickdraw.h> */
/*  #include <Fonts.h> */
/*  #include <Events.h> */
/*  #include <Menus.h> */
/*  #include <Windows.h> */
/*  #include <Palettes.h> */
/*  #include <TextEdit.h> */
/*  #include <Dialogs.h> */
/*  #include <OSUtils.h> */
/*  #include <SegLoad.h> */
/*  #include <Packages.h> */
/*  #include <Files.h> */

/* interface includes */
#include "devices.h"
#include "debug.h"
#include "initdev.h"

/* mif includes */
#include "MacMain.m"
#include "MacSurface.h"
#include "MacShell.h"
#include "MacGraph.h"

/*  #include "MacMain.h" */

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define OUTBUFFSIZE 1024        /* for buffering of the shell output			*/

#define DATA_MAX        2               /* max size of data returned by GUI_GetNextEvent*/

/* quick fix: printf (called by debug.c fcts) will stimulate the SIOUX-interface
   to open another terminal window wich belongs to the application but is
   unknown to ug */
#define my_assert(ass)          if ((ass)==FALSE) return (0)

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static int outbuffpos = 0;
static char outbuff[OUTBUFFSIZE] = "";

static ShellWindow shell;                                       /* our only shell window		*/

static EventRecord theEvent;    /* last event handled						*/
static INT MouseLocation[2];    /* last mouse location						*/
static char Char;                               /* last character from keyboard                         */
static INT ChosenTool;                  /* last chosen tool                                             */

static OUTPUTDEVICE *MacOutputDevice;   /* ptr to MacOutputDevice			*/

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


static void MacFlushOutbuff ()
{
  if (outbuffpos==0)
    return;

  MacWriteString(outbuff);

  /* reset outbuff */
  outbuff[0] = '\0';
  outbuffpos = 0;

  return;
}

void WriteString (const char *s)
{
  int len;

  assert(outbuff[outbuffpos]=='\0');

  if ((outbuffpos+(len=strlen(s)))<OUTBUFFSIZE-1)
  {
    strcpy(outbuff+outbuffpos,s);
    outbuffpos += len;
  }
  else
  {
    len = OUTBUFFSIZE-outbuffpos-1;
    strncpy(outbuff+outbuffpos,s,len);
    outbuffpos += len;
    MacFlushOutbuff();
    WriteString(s+len);
  }

  return;
}


/****************************************************************************/
/*
   GUI_GetNextEvent - Do system tasks, system events and pass next ug event

   SYNOPSIS:
   INT GetNextUGEvent (EVENT *reportEvent, INT EventMask);

   PARAMETERS:
   .  reportEvent -
   .  EventMask -

   DESCRIPTION:
   This function does system tasks, system events and passe next ug event.

   RETURN VALUE:
   INT
   .n    0 if no event occurred (ug or system)
   .n    1 if an event occurred (ug or system)
 */
/****************************************************************************/

INT GetNextUGEvent (EVENT *reportEvent, INT EventMask)
{
  WindowPtr whichWindow;
  GRAPH_WINDOW *gw;
  Point pt;
  INT WhereIn,rv;
  char *s;
  short MacEventMask;

  /* flush outbuff */
  if (outbuffpos!=0)
    MacFlushOutbuff();

  /* translate EventMask into machine specific one */
  switch (EventMask)
  {
  case EVERY_EVENT :
    MacEventMask = everyEvent;
    break;

  case TERM_STRING :
    MacEventMask = keyDownMask | autoKeyMask;
    break;


  case TERM_CMDKEY :
    MacEventMask = keyDownMask;
    break;

  default :
    return (EVENT_ERROR);
  }

  /* no event as default */
  reportEvent->Type = NO_EVENT;
  reportEvent->NoEvent.InterfaceEvent = 0;
  reportEvent->NoEvent.GraphWinActive = (WINDOWID) 0;

  /* do system tasks */
  SystemTask();

  /* periodic tasks for windows */
  if ((whichWindow=FrontWindow())!=NULL)
  {
    if (whichWindow==(WindowPtr)&(shell.theRecord))
      IdleShellWindow();
    gw = WhichGW(whichWindow);
    if (gw!=NULL)
    {
      GetMouse(&pt);
      reportEvent->NoEvent.GraphWinActive = (WINDOWID) gw;
      reportEvent->NoEvent.Mouse[0] = pt.h;
      reportEvent->NoEvent.Mouse[1] = pt.v;
    }
  }

  /* get event */
  if (GetNextEvent(MacEventMask,&theEvent)==TRUE)
  {
    WhereIn = FindWindow(theEvent.where,&whichWindow);

    /* handle events */
    switch (theEvent.what)
    {
    case osEvt :
      if (theEvent.message>>24==suspendResumeMessage)
      {
        if (theEvent.message & resumeFlag)
        {
          whichWindow = FrontWindow();
          if (whichWindow==(WindowPtr)&(shell.theRecord))
          {
            /* shell window */
            ActivateShellWin();
            break;
          }
          else
          {
            /* graph window */
            gw = WhichGW(whichWindow);
            my_assert(gw!=NULL);
            Mac_ActivateOutput((WINDOWID)gw);
            break;
          }
        }
        else
        {
          whichWindow = FrontWindow();
          if (whichWindow==(WindowPtr)&(shell.theRecord))
          {
            /* shell window */
            DeactivateShellWin();
            break;
          }
        }
      }
      break;

    case mouseDown :
      switch (WhereIn)
      {
      /* do system events */
      case inSysWindow :
        SystemClick(&theEvent,whichWindow);
        break;
      case inMenuBar :
        DoCommand(MenuSelect(theEvent.where));
        break;

      /* do ug events */
      case inGrow :
        if (whichWindow==(WindowPtr)&(shell.theRecord))
        {
          /* shell window */
          if ((rv=GrowShellWindow (&theEvent))==NO_POS_CHANGE)
          {
            SelectWindow(whichWindow);
            ActivateShellWin();
          }
        }
        else
        {
          /* graph window */
          gw = WhichGW(whichWindow);
          my_assert(gw!=NULL);
          if (GrowGraphWindow(gw,&theEvent,&(reportEvent->DocGrow))==POS_CHANGE)
          {
            reportEvent->Type                       = DOC_GROW;
            reportEvent->DocGrow.win        = (WINDOWID) gw;
          }
          else if (whichWindow!=FrontWindow())
          {
            SelectWindow(whichWindow);

            gw = WhichGW(whichWindow);
            my_assert(gw!=NULL);
            Mac_ActivateOutput((WINDOWID)gw);
            reportEvent->Type                               = DOC_ACTIVATE;
            reportEvent->DocActivate.win    = (WINDOWID) gw;
          }
        }
        break;
      case inDrag :
        if (whichWindow==(WindowPtr)&(shell.theRecord))
        {
          /* shell window */
          DragShellWin(&theEvent);
        }
        else
        {
          /* graph window */
          gw = WhichGW(whichWindow);
          my_assert(gw!=NULL);
          if (DragGraphWindow(gw,&theEvent,&(reportEvent->DocDrag))==POS_CHANGE)
          {
            reportEvent->Type                       = DOC_DRAG;
            reportEvent->DocDrag.win        = (WINDOWID) gw;
          }
          else if (whichWindow!=FrontWindow())
          {
            SelectWindow(whichWindow);

            gw = WhichGW(whichWindow);
            my_assert(gw!=NULL);
            Mac_ActivateOutput((WINDOWID)gw);
            reportEvent->Type                               = DOC_ACTIVATE;
            reportEvent->DocActivate.win    = (WINDOWID) gw;
          }
        }
        break;
      case inGoAway :
        if (TrackGoAway(whichWindow,theEvent.where))
          if (whichWindow==(WindowPtr)&(shell.theRecord))
          {
            /* shell window */
            reportEvent->Type                               = TERM_GOAWAY;
          }
          else
          {
            /* graph window */
            gw = WhichGW(whichWindow);
            my_assert(gw!=NULL);
            reportEvent->Type                               = DOC_GOAWAY;
            reportEvent->DocGoAway.win              = (WINDOWID) gw;
          }
        break;
      case inContent :
        if (whichWindow!=FrontWindow())
        {
          SelectWindow(whichWindow);

          if (whichWindow==(WindowPtr)&(shell.theRecord))
          {
            /* shell window */
            ActivateShellWin();
            break;
          }
          else
          {
            /* graph window */
            gw = WhichGW(whichWindow);
            my_assert(gw!=NULL);
            Mac_ActivateOutput((WINDOWID)gw);
            reportEvent->Type                               = DOC_ACTIVATE;
            reportEvent->DocActivate.win    = (WINDOWID) gw;
            break;
          }
        }
        else
        {
          if (whichWindow==(WindowPtr)&(shell.theRecord))
          {
            ShellWinContentClick(whichWindow,&theEvent);
          }
          else
          {
            /* graph window */
            gw = WhichGW(whichWindow);
            my_assert(gw!=NULL);
            SetPort(whichWindow);
            GlobalToLocal(&(theEvent.where));
            MouseLocation[0] = (INT)theEvent.where.h;
            MouseLocation[1] = (INT)theEvent.where.v;
            if (WhichTool((WINDOWID) gw,MouseLocation,&ChosenTool))
            {
              reportEvent->Type                               = DOC_CHANGETOOL;
              reportEvent->DocChangeTool.win  = (WINDOWID) gw;
              reportEvent->DocChangeTool.Tool = ChosenTool;
              reportEvent->DocChangeTool.MousePosition[0]   = MouseLocation[0];
              reportEvent->DocChangeTool.MousePosition[1]   = MouseLocation[1];
            }
            else
            {
              reportEvent->Type = DOC_CONTENTCLICK;
              reportEvent->DocContentClick.win                                = (WINDOWID) gw;
              reportEvent->DocContentClick.MousePosition[0]   = (DOUBLE) MouseLocation[0];
              reportEvent->DocContentClick.MousePosition[1]   = (DOUBLE) MouseLocation[1];
            }
          }
        }
        break;
      }
      break;

    case activateEvt :
      whichWindow = (WindowPtr)theEvent.message;
      if (whichWindow==(WindowPtr)&(shell.theRecord))
      {
        /* shell window */
        if (BitAnd(theEvent.modifiers,activeFlag)!=0)
          ActivateShellWin();
        else
          DeactivateShellWin();
      }
      else
      {
        /* graph window */
        gw = WhichGW(whichWindow);
        my_assert(gw!=NULL);

        if (BitAnd(theEvent.modifiers,activeFlag)!=0)
        {
          Mac_ActivateOutput((WINDOWID)gw);
          reportEvent->Type                               = DOC_ACTIVATE;
          reportEvent->DocActivate.win    = (WINDOWID) gw;
        }
      }
      break;

    case updateEvt :
      whichWindow = (WindowPtr)theEvent.message;
      if (whichWindow==(WindowPtr)&(shell.theRecord))
      {
        /* shell window */
        if (whichWindow==FrontWindow())
          SetMyCursor(textCurs);
        UpdateShellWin();
      }
      else
      {
        /* graph window */
        gw = WhichGW(whichWindow);
        my_assert(gw!=NULL);

        if (whichWindow==FrontWindow())
          SetMyCursor(gw->currTool);
        Mac_UpdateOutput((WINDOWID)gw,arrowTool);
        reportEvent->Type                       = DOC_UPDATE;
        reportEvent->DocUpdate.win      = (WINDOWID) gw;
      }
      break;

    case keyDown :
    case autoKey :
    {
      char c=(char)BitAnd(theEvent.message,charCodeMask);
      if (BitAnd(theEvent.modifiers,cmdKey))
      {
        if (c=='\34')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_PAGE_LEFT,0);
          break;
        }
        if (c=='\35')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_PAGE_RIGHT,0);
          break;
        }
        if (c=='\36')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_PAGE_UP,0);
          break;
        }
        if (c=='\37')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_PAGE_DOWN,0);
          break;
        }
        reportEvent->Type = TERM_CMDKEY;
        reportEvent->TermCmdKey.CmdKey = c;

        /* cmd + shift + 1 --> cmd + ! */
        if (reportEvent->TermCmdKey.CmdKey=='1')
          if (BitAnd(theEvent.modifiers,shiftKey))
            reportEvent->TermCmdKey.CmdKey = '!';
        break;
      }
      else if (BitAnd(theEvent.modifiers,optionKey))
      {
        if (c=='\34')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_COL_LEFT,0);
          break;
        }
        if (c=='\35')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_COL_RIGHT,0);
          break;
        }
        if (c=='\36')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_LINE_UP,0);
          break;
        }
        if (c=='\37')
        {
          reportEvent->NoEvent.InterfaceEvent = 1;
          ShellHandleKeybordEvent(SK_LINE_DOWN,0);
          break;
        }
      }
      if (c=='\f')
      {
        reportEvent->NoEvent.InterfaceEvent = 1;
        ShellHandleKeybordEvent(SK_PAGE_DOWN,0);
        break;
      }
      else if (c=='\v')
      {
        reportEvent->NoEvent.InterfaceEvent = 1;
        ShellHandleKeybordEvent(SK_PAGE_UP,0);
        break;
      }
      else if ((s=ShellHandleKeybordEvent(SK_NONE,c))==NULL)
      {
        reportEvent->NoEvent.InterfaceEvent = 1;
        break;
      }
      else
      {
        reportEvent->Type = TERM_STRING;
        strcpy(reportEvent->TermString.String,s);
      }
      break;
    }
    }
  }
  return (0);
}

/****************************************************************************/
/*
   InitScreen - Init rest of GUI and return ptr to screen outputdevice

   SYNOPSIS:
   OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error);

   PARAMETERS:
   .  argcp - pointer to argument counter
   .  argv  - argument vector
   .  error -

   DESCRIPTION:
   This function inits rest of GUI and returns ptr to screen outputdevice.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      POINTER to
   .n      0 if an error occurred.
 */

/****************************************************************************/

OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error)
{
  /* init shell etc */
  if (ShellInitAndOpen(&shell)!=0)
  {
    *error = 2;
    return (NULL);
  }
  if ((MacOutputDevice=InitMacOutputDevice())==NULL)
  {
    *error = 3;
    return (NULL);
  }
  if (InitMacSurface())
  {
    *error = 1;
    return (NULL);
  }

  *error = 0;
  return (MacOutputDevice);
}



void ExitScreen (void)
{}
