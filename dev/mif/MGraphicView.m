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

	[self addTrackingRect:frameRect owner:self userData:nil assumeInside:NO];
	
    arrowImage = [NSImage imageNamed:@"Arrow"];
    magnifyImage = [NSImage imageNamed:@"Magnify"];
	
    arrowCursor		= [[NSCursor alloc] initWithImage:arrowImage
                                              hotSpot:NSMakePoint(4,0)];
    magnifyCursor	= [[NSCursor alloc] initWithImage:magnifyImage
                                              hotSpot:NSMakePoint(4,4)];

    if ( [magnifyCursor image]==nil )
        printf("No image for magnify cursor\n");
    
	[self addCursorRect:frameRect cursor:magnifyCursor];
    currentCursor = arrowCursor;
    
    [magnifyCursor setOnMouseEntered:YES];
    
    return self;
}

- (void) keyDown:(NSEvent *)theEvent
{
	[theShellTextView keyDown:theEvent];
}

- (void)drawRect:(NSRect)aRect
{
	//PSsetgray(NSWhite);
	//NSRectFill(aRect);
	printf("drawRect called\n");
    //InterpretCommand("plot");
}

- (void) mouseEntered:(NSEvent *)anEvent
{
    printf("MouseEntered\n");
}

- (void) mouseExited:(NSEvent *)theEvent
{
    printf("MouseExited\n");    
}

- (void)mouseDown:(NSEvent *)theEvent
{
   BOOL keepOn = YES;
   BOOL isInside = YES;
   NSPoint mouseLoc;
   
   do {
       mouseLoc = [self convertPoint:[theEvent locationInWindow] fromView:nil];
       isInside = [self mouse:mouseLoc inRect:[self bounds]];
       switch ([theEvent type]) {
         case NSLeftMouseDragged:
               //[self highlight:isInside];
               break;
         case NSLeftMouseUp:
               //if (isInside) [self doSomethingSignificant];
               //[self highlight:NO];
               keepOn = NO;
               break;
         default:
               /* Ignore any other kind of event. */
               break;
      }
      theEvent = [[self window] nextEventMatchingMask: NSLeftMouseUpMask |
            NSLeftMouseDraggedMask];
    }	while (keepOn);
	return;
}

- (BOOL) acceptsFirstResponder	{ return YES; }

- (void) resetCursorRects
{
    printf("resetCursorRects called\n");
    [self addCursorRect:[self frame] cursor:currentCursor];
}

@end

