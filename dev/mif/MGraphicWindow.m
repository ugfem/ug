/****************************************************************************/
/*																			*/
/* File:	  MGraphicWindow.m												*/
/*																			*/
/* Purpose:   Implementation of MGraphicWindow class.						*/
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


#import "MGraphicWindow.h"

#define BORDERSIZE 2
#define BOTTOMBORDERSIZE 1
#define SEPARATORSIZE 3
#define STATUSTEXTHEIGHT 21

@implementation MGraphicWindow

- (id)initWithContentRect:(NSRect)contentRect
                styleMask:(unsigned int)styleMask
                backing:(NSBackingStoreType)backingType
                defer:(BOOL)flag
{
    NSRect viewRect;

    [super initWithContentRect:NSMakeRect(contentRect.origin.x,
                                          contentRect.origin.y,
                                          contentRect.size.width+BORDERSIZE+BORDERSIZE,
                                          contentRect.size.height+BORDERSIZE+BOTTOMBORDERSIZE+SEPARATORSIZE+STATUSTEXTHEIGHT)
                     styleMask:styleMask backing:backingType defer:flag];

    [self setBackgroundColor:[NSColor whiteColor]];

    // Now set up the MGraphicView inside the GraphicWindow
    viewRect = [[self contentView] bounds];
    theView = [[MGraphicView alloc]
            initWithFrame:NSMakeRect(BORDERSIZE,
                                     BORDERSIZE+BOTTOMBORDERSIZE+SEPARATORSIZE+STATUSTEXTHEIGHT,
                                     contentRect.size.width,
                                     contentRect.size.height)];
    [[self contentView] addSubview:theView];

    // Set up the status text line
    statusText = [[NSTextView alloc]
			initWithFrame:NSMakeRect(BORDERSIZE,
                                     BOTTOMBORDERSIZE,
                                     contentRect.size.width,
                                     contentRect.size.height)];
    [statusText setEditable:NO];
    [statusText setSelectable:YES];
    [statusText setBackgroundColor:[NSColor whiteColor]];
    [statusText insertText:@"Status"];
    
	[self display];
    [self makeKeyAndOrderFront:nil];
    [self makeFirstResponder:theView];

    /* Let the graphics window accept mouse moved events (for tracking of
       x and y coordinates).  Default is NO for NSWindow objects.*/
    [self setAcceptsMouseMovedEvents:YES];

    return self;
}

- (void) close
{
    [theView dealloc];
	[self dealloc];
    return;
}

- (void) activateOutput
{
	[self makeKeyAndOrderFront:nil];
}

- (void) updateOutput
{
    [theView drawRect:[theView frame]];
}


- (void) moveToPoint:(SHORT_POINT)point
{
	moveto_x = point.x;
	moveto_y = point.y;
}

- (void) drawLineTo:(SHORT_POINT)point
{
	[theView lockFocus];
	PSmoveto(moveto_x, moveto_y);
	PSlineto(point.x, point.y);
	PSstroke();
	[theView unlockFocus];
}


- (void) drawPolyLine:(SHORT_POINT *)points noOfPoints:(INT)n
{
	int i;
	
	if (n<2) return;

	[theView lockFocus];
	PSmoveto(points[0].x, points[0].y);
	for (i=1; i<n; i++) 
		PSlineto(points[i].x, points[i].y);
	PSstroke();
	[theView unlockFocus];
}


- (void) drawInversePolyLine:(SHORT_POINT *)points noOfPoints:(INT)n
{
	return;
}


- (void) drawPolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
	int i;
	
	if (n<3) return;
	
	[theView lockFocus];
	PSmoveto(points[0].x, points[0].y);
	for (i=1; i<n; i++) PSlineto (points[i].x, points[i].y);
	PSclosepath();
	PSfill();
	PSstroke();
	[theView unlockFocus];
	return;
}


- (void) drawInversePolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
	return;
}


- (void) erasePolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
	return;
}


- (void) marker
{
	return;
}


- (void) polyMark
{
	return;
}


- (void) invMarker
{
	return;
}


- (void) invPolyMark
{
	return;
}


- (void) drawText:(NSString) text
{
    //printf("drawText: <%s>\n", text);
	return;
}

- (void) drawCenteredText:(NSString) text
{
    //printf("drawCenteredText: <%s>\n", text);
	return;
}


- (void) clearView
{
	return;
}

- (void) setLineWidth:(short)width
{
	PSsetlinewidth((float)width);
	return;
}

- (void) setTextSize:(short)width
{
	return;
}

- (void) setMarkerSize:(short)width
{
	return;
}

- (void) setMarker:(short)width
{
	return;
}

- (void) setPaletteEntryWithIndex:(long)index red:(short)r green:(short)g blue:(short)b
{
	return;
}

- (void) getPaletteEntryWithIndex:(long)index red:(short *)r green:(short *)g blue:(short *)b
{
	return;
}

- (void) flush
{
	[theView lockFocus];
	PSstroke();
	[self flushWindow];
	[theView unlockFocus];
	return;
}

@end
