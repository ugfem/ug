// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MacSurface.h													*/
/*																			*/
/* Purpose:   header file for MacSurface.c									*/
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
/* Remarks:   former MacGuSurfacei.c/.h                                                                         */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __MACSURFACE__
#define __MACSURFACE__

#include <Quickdraw.h>

#ifndef __COMPILER__
#include "compiler.h"
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

/* screen's size from the Quidraw globals */
#define SCREEN_WIDTH            (((qd.screenBits).bounds).right)
#define SCREEN_HEIGHT           (((qd.screenBits).bounds).bottom)

#define MENU_BAR                        39      /* the menu-bar has a hight of 39 pixels	*/
#define MARGIN_TO_SCREEN        2       /* leave a margin between window and screen */
#define TITLE_BAR                       19      /* height of the tiltel bar of a window         */
#define SCROLL_BAR                      15      /* width of h/v scrollbars					*/

#define NO_POS_CHANGE           0       /* drag and grow return code				*/
#define POS_CHANGE                      1       /* drag and grow return code				*/

#define HiWrd(aLong)            (((aLong) >> 16) & 0xFFFF)
#define LoWrd(aLong)            ((aLong) & 0xFFFF)

#define MAC_WIN(w)                      ((WindowPtr) & ((w)->theRecord))

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

Rect *DragRect                  ();
void SetMyCursor                (short id);
void About                              ();
void ScheduleCommand    (short theMenu,short theItem);
void DoCommand                  (long mResult);
void SetUpMenus                 ();

int  InitMacSurface     (void);

#endif
