// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  UGShellWindow.h												*/
/*																			*/
/* Purpose:   Interface to UGShellWindow object.							*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/*	History:  January 21, 1998 begin										*/
/*																			*/
/****************************************************************************/

#ifndef _UGSHELLWINDOW_
#define _UGSHELLWINDOW_

#import <appkit/appkit.h>

@interface UGShellWindow : Window
{
  ScrollView      *theScrollView;
  Text            *theText;
  Panel           *infoPanel;
  id theMenu;
  id shellFont;
}

- setUp;
- free;
- shellText;
- appendToText:(const char *)val;
@end

#endif
