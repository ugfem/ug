// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MShellWindow.h												*/
/*																			*/
/* Purpose:   Interface to MShellWindow class.								*/
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

#ifndef _MShellWindow_
#define _MShellWindow_

#import <AppKit/AppKit.h>

@interface MShellWindow : NSWindow
{
  NSCharacterSet  *finalCharSet;
  NSMutableString     *commandString;
}

- (void)keyDown:(NSEvent *)theEvent;

- (id)initWithContentRect:(NSRect) contentRect
 styleMask:(unsigned int)styleMask
 backing:(NSBackingStoreType) backingType
 defer:(BOOL)flag;
@end

#endif
