/****************************************************************************/
/*																			*/
/* File:	  MGraphicView.m												*/
/*																			*/
/* Purpose:   Implementation of the MGraphicView class.						*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Interdisziplin"ares Zentrum f"ur Wissenschaftliches Rechnen	*/
/*			  Universit"at Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69210 Heidelberg												*/
/*			  email: Volker.Reichenberger@IWR.Uni-Heidelberg.DE		    	*/
/*																			*/
/*	History:  June 4, 1999 begin (based on OPENSTEP code)					*/
/*																			*/
/****************************************************************************/


#import "MGraphicView.h"

#include "cmdint.h"

@implementation MGraphicView

- (id) initWithFrame:(NSRect)frameRect
{
    [super initWithFrame:frameRect];
	theShellTextView = [[MShell instantiate] textview];
    return self;
}

- (void) keyDown:(NSEvent *)theEvent
{
	[theShellTextView keyDown:theEvent];
}

- (void)drawRect:(NSRect)aRect
{
	PSsetgray(NSWhite);
	NSRectFill(aRect);

    //InterpretCommand("plot");
}
@end

