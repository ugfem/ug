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
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*					 henrik@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7007									*/
/*			  fax  : 0049-(0)711-685-7000									*/
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
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifndef __MWCW__        /* don't need that: included MacHeadersPPC */

/* mac toolbox includes */
#include <desk.h>
#include <ToolUtils.h>
#include <Types.h>
#include <Memory.h>
#include <Quickdraw.h>
#include <Fonts.h>
#include <Events.h>
#include <Menus.h>
#include <Windows.h>
#include <Palettes.h>
#include <TextEdit.h>
#include <Dialogs.h>
#include <OSUtils.h>
#include <SegLoad.h>
#include <osevents.h>
#include <Packages.h>
#include <Files.h>

#endif

/* interface includes */
#include "compiler.h"
#include "misc.h"
#include "devices.h"
#include "initdev.h"
#include "heaps.h"

/* mif includes */
#include "MacMain.m"
#include "MacSurface.h"
#include "MacShell.h"
#include "MacGraph.h"

#include "MacMain.h"

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

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static int outbuffpos;
static char outbuff[OUTBUFFSIZE];

static ShellWindow shell;                                       /* our only shell window		*/

static EventRecord theEvent;    /* last event handled						*/
static INT MouseLocation[2];    /* last mouse location						*/
static char Char;                               /* last character from keyboard                         */
static INT ChosenTool;                  /* last chosen tool                                             */

OUTPUTDEVICE *MacOutputDevice;  /* ptr to MacOutputDevice					*/

/* data for CVS */
static char rcsid[] = "$Header$";

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

  if ((outbuffpos+(len=strlen(s)))<OUTBUFFSIZE-1)
  {
    outbuffpos += len;
    strcat(outbuff,s);
  }
  else
  {
    len = OUTBUFFSIZE-outbuffpos-1;
    outbuffpos += len;
    strncat(outbuff,s,len);
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

  /* do system tasks */
  SystemTask();

  /* periodic tasks for windows */
  if ((whichWindow=FrontWindow())!=NULL)
  {
    if (whichWindow==(WindowPtr)&(shell.theRecord))
      IdleShellWindow();
  }

  /* get event */
  if (GetNextEvent(MacEventMask,&theEvent)==TRUE)
  {
    WhereIn = FindWindow(theEvent.where,&whichWindow);

    /* handle events */
    switch (theEvent.what)
    {
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
          assert(gw!=NULL);
          if (GrowGraphWindow(gw,&theEvent,&(reportEvent->DocGrow))==POS_CHANGE)
          {
            reportEvent->Type                       = DOC_GROW;
            reportEvent->DocGrow.win        = (WINDOWID) gw;
          }
          else if (whichWindow!=FrontWindow())
          {
            SelectWindow(whichWindow);

            gw = WhichGW(whichWindow);
            assert(gw!=NULL);
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
          assert(gw!=NULL);
          if (DragGraphWindow(gw,&theEvent,&(reportEvent->DocDrag))==POS_CHANGE)
          {
            reportEvent->Type                       = DOC_DRAG;
            reportEvent->DocDrag.win        = (WINDOWID) gw;
          }
          else if (whichWindow!=FrontWindow())
          {
            SelectWindow(whichWindow);

            gw = WhichGW(whichWindow);
            assert(gw!=NULL);
            Mac_ActivateOutput((WINDOWID)gw);
            reportEvent->Type                               = DOC_ACTIVATE;
            reportEvent->DocActivate.win    = (WINDOWID) gw;
          }
        }
        break;
      case inGoAway :
        if (whichWindow==(WindowPtr)&(shell.theRecord))
        {
          /* shell window */
          reportEvent->Type                               = TERM_GOAWAY;
        }
        else
        {
          /* graph window */
          gw = WhichGW(whichWindow);
          assert(gw!=NULL);
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
            assert(gw!=NULL);
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
            assert(gw!=NULL);
            SetPort(whichWindow);
            GlobalToLocal(&(theEvent.where));
            MouseLocation[0] = (INT)theEvent.where.h;
            MouseLocation[1] = (INT)theEvent.where.v;
            if (GetTool(whichWindow,MouseLocation,&ChosenTool))
            {
              reportEvent->Type                               = DOC_CHANGETOOL;
              reportEvent->DocChangeTool.win  = (WINDOWID) gw;
              reportEvent->DocChangeTool.Tool = ChosenTool;
            }
            else
            {
              reportEvent->Type = DOC_CONTENTCLICK;
              reportEvent->DocContentClick.win                                = (WINDOWID) gw;
              reportEvent->DocContentClick.MousePosition[0]   = (COORD) MouseLocation[0];
              reportEvent->DocContentClick.MousePosition[1]   = (COORD) MouseLocation[1];
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
        assert(gw!=NULL);

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
        UpdateShellWin();
      }
      else
      {
        /* graph window */
        gw = WhichGW(whichWindow);
        assert(gw!=NULL);

        Mac_UpdateOutput((WINDOWID)gw,NULL,arrowTool);
        reportEvent->Type                       = DOC_UPDATE;
        reportEvent->DocUpdate.win      = (WINDOWID) gw;
      }
      break;

    case keyDown :
    case autoKey :
      if (BitAnd(theEvent.modifiers,cmdKey))
      {
        reportEvent->Type = TERM_CMDKEY;
        reportEvent->TermCmdKey.CmdKey = (char)BitAnd(theEvent.message,charCodeMask);

        /* special: cmd + shift doesn't work */
        if (reportEvent->TermCmdKey.CmdKey=='Á')
          reportEvent->TermCmdKey.CmdKey = '!';
        break;
      }
      else if ((s=ShellHandleKeybordEvent((char)BitAnd(theEvent.message,charCodeMask)))==NULL)
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
  return (0);
}

/****************************************************************************/
/*
   InitScreen - Init rest of GUI and return ptr to screen outputdevice

   SYNOPSIS:
   OUTPUTDEVICE *InitScreen (int argc, char **argv, INT *error);

   PARAMETERS:
   .  argc - argument counter
   .  argv - argument vector
   .  error -

   DESCRIPTION:
   This function inits rest of GUI and returns ptr to screen outputdevice.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      POINTER to
   .n      0 if an error occurred.
 */

/****************************************************************************/

OUTPUTDEVICE *InitScreen (int argc, char **argv, INT *error)
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
