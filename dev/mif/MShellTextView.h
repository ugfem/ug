// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MShellTextView.h												*/
/*																			*/
/* Purpose:   Interface for the MShellTextView class.						*/
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


#import <AppKit/AppKit.h>


#define MAX_BUFFERED_COMMANDS 128


@interface MShellTextView : NSTextView
{
  int startPosition;                                            // This is the start of the current command.
  NSTextStorage *theTextStore;

  NSMutableArray      *theCommandList;          // store the history of entered commands
  unsigned int commandCount;                    // count the number of entered commands
  unsigned int currCommandNumber;               // for scrolling through command list
  NSMutableString     *currentLine;             // the current input text

  NSDictionary *textAttributes;                 // Normal text attributes
  NSDictionary *errorAttributes;                // Error highlighting attributes
  NSDictionary *promptAttributes;               // Prompt text attributes
  NSDictionary *cmdlineAttributes;      // Command line text attributes
}

- (id) initWithFrame:(NSRect)frameRect;

- (void) appendPrompt;

- (void) appendText:(NSString *)str;

- (void)keyDown:(NSEvent *)theEvent;

- (void) replaceCurrentCommandWith:(NSString *)command;

- (void)highlightError;

@end
