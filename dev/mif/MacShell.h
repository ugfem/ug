// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MacTerminalWindows.h											*/
/*																			*/
/* Purpose:   header file for scrollable text windows						*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					*/
/*																			*/
/* History:   27.02.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __MACSHELL__
#define __MACSHELL__

#include <Windows.h>
#include <TextEdit.h>
#include <Controls.h>

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/********************************************************************************/
/*																				*/
/* data structures exported by the corresponding source file					*/
/*																				*/
/********************************************************************************/

typedef struct {

  WindowRecord theRecord;                       /* used by mac's window manager             */

  /* now terminal window specific stuff */
  TEHandle textH;                                       /* text edit handle                                             */
  int lineHeight;                                       /* height of one line in pixels                         */
  int charWidth;                                        /* max width of a character in pixels		*/
  int line;                                                     /* line number of first line in window		*/
  int col;                                                      /* col number of first col in window		*/
  int lines;                                                    /* no of lines the window can display		*/
  int cols;                                                     /* no of cols the window can display		*/
  ControlHandle vScrollBar;                     /* vertical scroll bar						*/
  ControlHandle hScrollBar;                     /* horizontal scroll bar					*/

} ShellWindow;

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT              GrowShellWindow                        (EventRecord *theEvent);
INT              DragShellWin                           (EventRecord *theEvent);
void     ShellWinContentClick           (WindowPtr theWindow,EventRecord *theEvent);
INT      ActivateShellWin                       ();
INT      DeactivateShellWin             ();
INT      UpdateShellWin                         (void);
void     IdleShellWindow                        ();
char    *ShellHandleKeybordEvent        (char key);
void     MacWriteString                         (char *s);

int      ShellInitAndOpen                       (ShellWindow *sh);

#endif
