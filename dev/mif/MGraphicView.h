// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MGraphicView.h												*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Interdisziplin"ares Zentrum f"ur Wissenschaftliches Rechnen	*/
/*			  Universit"at Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69210 Heidelberg												*/
/*			  email: Volker.Reichenberger@IWR.Uni-Heidelberg.DE		        */
/*																			*/
/*	History:  June 4, 1999 begin (based on OPENSTEP code)					*/
/*																			*/
/****************************************************************************/

#ifndef _MGRAPHICVIEW_
#define _MGRAPHICVIEW_

#import <AppKit/AppKit.h>

#import "MShell.h"

@interface MGraphicView : NSView
{
  MShellTextView *theShellTextView;
}

- (void) drawRect:(NSRect)aRect;
- (void) keyDown: (NSEvent *)theEvent;

@end

#endif
