// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MShell.h														*/
/*																			*/
/* Purpose:   Interface to MShell class.									*/
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

#ifndef _MSHELL_
#define _MSHELL_

#import <AppKit.h>

#import "MShellView.h"
#import "MShellWindow.h"
#import "MShellTextView.h"

@interface MShell : NSObject
{
  MShellWindow    *theWindow;
  MShellView          *theShellView;
  NSScrollView    *theScrollView;
  MShellTextView  *theTextView;
  NSMenu                  *theMenu;
  id shellFont;
  char                    *inpLine;
}

+ (MShell *)    instantiate;
- (BOOL)        setUp;
- (void)        applicationDidFinishLaunching:(NSNotification *)aNotification;
- (void)        appendToText:(const char *)val;
- (void)        appendPrompt;
- (BOOL)        executeFile:(NSString *)fileName;
- (NSPoint)     getMouseLocation;
- (void)        setScrollback:(int)sb;
- (void)        setCharactersPerLine:(int)cpl;
- (void)        setFontSize:(float)fsize;
- (BOOL)        interpretCommand:(NSString*)command;

- (MShellWindow*) window;
/*- (BOOL)textView:(NSTextView *)aTextView shouldChangeTextInRange:(NSRange)affectedCharRange replacementString:(NSString *)replacementString;
   - (NSRange)textView:(NSTextView *)aTextView willChangeSelectionFromCharacterRange:(NSRange)oldSelectedCharRange toCharacterRange:(NSRange)newSelectedCharRange;*/
@end

#endif
