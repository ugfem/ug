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

#import "MShell.h"
#import "MGraphicWindow.h"

#define BORDERSIZE 2
#define BOTTOMBORDERSIZE 2
#define SEPARATORSIZE 3
#define STATUSTEXTHEIGHT 21

static NSMutableDictionary *textAttributes;



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

    [self setBackgroundColor:[NSColor grayColor]];

    // Now set up the MGraphicView inside the GraphicWindow
    viewRect = [[self contentView] bounds];
    theView = [[MGraphicView alloc]
            initWithFrame:NSMakeRect(BORDERSIZE,
                                     BOTTOMBORDERSIZE+SEPARATORSIZE+STATUSTEXTHEIGHT,
                                     contentRect.size.width,
                                     contentRect.size.height)];
    [[self contentView] addSubview:theView];

	// set default colors
    color = [NSColor blackColor];
    backgroundColor = [NSColor whiteColor];
    textColor = [NSColor whiteColor];

    // Set up the status text line
    statusText = [[NSTextView alloc]
			initWithFrame:NSMakeRect(BORDERSIZE+9,
                            BOTTOMBORDERSIZE,
                            contentRect.size.width-18,
                            STATUSTEXTHEIGHT)];

    [[self contentView] addSubview:statusText];
    [statusText alignCenter:statusText]; 
    
    [statusText setSelectable:YES];
    [statusText setDrawsBackground:YES];
    [statusText setBackgroundColor:[self backgroundColor]];
    [statusText setTextColor:textColor];
    [statusText setEditable:YES];
    [statusText insertText:@"Status"];
    [statusText setEditable:NO];    
    [statusText display];
    
	[self display];
    [self orderFront:self];
    [self makeFirstResponder:theView];
    [self setDelegate:[[MShell instantiate] window]];

    /* Let the graphics window accept mouse moved events (for tracking of
       x and y coordinates).  Default is NO for NSWindow objects.*/
    [self setAcceptsMouseMovedEvents:YES];

    textAttributes = [NSDictionary dictionaryWithObjectsAndKeys:
                            [NSColor blackColor],NSForegroundColorAttributeName,
                            [NSColor whiteColor],NSBackgroundColorAttributeName,
                            [NSFont boldSystemFontOfSize:10.0], NSFontAttributeName,
                            nil];
    textsize = 10.0;
    
    return self;
}

- (void) close
{
    [theView dealloc];
    [super close];
    return;
}

- (void) activateOutput
{
    [self orderFront:self];
}

- (void) updateOutput
{
    printf("Update Output\n");
}


- (void) moveToPoint:(SHORT_POINT)point
{
	moveto_x = point.x;
	moveto_y = point.y;
}

- (void) drawLineTo:(SHORT_POINT)point
{
	[theView lockFocus];
    [NSBezierPath strokeLineFromPoint:(NSPoint){moveto_x, moveto_y}
                      toPoint:(NSPoint){point.x, point.y}];
	//PSmoveto(moveto_x, moveto_y);
	//PSlineto(point.x, point.y);
	//PSstroke();
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
	int i;
    printf("drawInversePolyLine\n");
    [theView lockFocus];
    [[NSColor blueColor] set];
    [path moveToPoint:(NSPoint){points[0].x, points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, points[i].y}];
    [path closePath];
    [path fill];
    [path removeAllPoints];	// There is no 'reset' method, though the manual says so...
    //[color set];
    [theView unlockFocus];
	return;
}


- (void) drawPolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
	int i;    
	if (n<3) return;

    [theView lockFocus];
    [path moveToPoint:(NSPoint){points[0].x, points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, points[i].y}];
    [path closePath];
    [path fill];
    [path removeAllPoints];	// There is no 'reset' method, though the manual says so...
    [theView unlockFocus];
	return;
}


- (void) drawShadedPolygon:(SHORT_POINT *)points noOfPoints:(INT)n intensity:(DOUBLE)it
{
    printf("drawShadedPolygon\n");
	return;
}

- (void) drawInversePolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
    int i;
    printf("drawInversePolygon\n");
    [theView lockFocus];
    [[NSColor greenColor] set];
    [path moveToPoint:(NSPoint){points[0].x, points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, points[i].y}];
    [path closePath];
    [path fill];
    [path removeAllPoints];	// There is no 'reset' method, though the manual says so...
    //[color set];
    [theView unlockFocus];
    return;
}


- (void) erasePolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
    int i;
    [theView lockFocus];
    [backgroundColor set];
    /*[path moveToPoint:(NSPoint){points[0].x, points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, points[i].y}];
    [path closePath];
    [path stroke];
    [path removeAllPoints];	// There is no 'reset' method, though the manual says so...
	*/
    PSmoveto(points[0].x, points[0].y);
    for (i=1; i<n; i++) PSlineto (points[i].x, points[i].y);
    PSclosepath();
    PSfill();
	//PSstroke();
    [theView unlockFocus];
    //[color set];
	return;
}


- (void) setMarkerSize:(short)width
{
    printf("setMarkerSize\n");
    return;
}

- (void) setMarker:(short)width
{
    printf("setMarker\n");
    return;
}

- (void) marker
{
    printf("marker\n");
	return;
}


- (void) polyMark:(SHORT_POINT *)points noOfPoints:(INT)n
{
    printf("polyMark\n");
	return;
}


- (void) invMarker:(SHORT_POINT)point
{
    printf("invMarker\n");
	return;
}


- (void) invPolyMark:(SHORT_POINT *)points noOfPoints:(INT)n
{
    printf("invPolyMark\n");
	return;
}


- (void) drawText:(const char*) text mode:(INT)m
{
    NSAttributedString *s = [[NSAttributedString alloc] initWithString:[NSString stringWithCString:text]
                                   attributes:textAttributes];
    [theView lockFocus];
	[s drawAtPoint:(NSPoint){moveto_x,moveto_y}];
    [theView unlockFocus];
    //printf("drawText, mode %d: <%s>\n", m, text);
	return;
}

- (void) drawCenteredText:(const char*) text atPoint:(SHORT_POINT)point mode:(INT)m
{
    NSAttributedString *s = [[NSAttributedString alloc] initWithString:[NSString stringWithCString:text]
                                   attributes:textAttributes];
    [theView lockFocus];
    [s drawAtPoint:(NSPoint){moveto_x+40,moveto_y}];
    [theView unlockFocus];
    //printf("drawCenteredText: <%s>\n", text);
	return;
}


- (void) clearView
{
    printf("clearView\n");
    [theView lockFocus];
    [[NSColor blueColor] set];
    [[NSBezierPath bezierPathWithRect:[theView frame]] fill];
    [theView unlockFocus];
	return;
}

- (void) setColor:(NSColor*)col
{
	/* dummy */
}

- (void) setLineWidth:(short)width
{
	//PSsetlinewidth((float)width);
    [path setLineWidth:(float)width];
	return;
}

- (void) setTextSize:(short)h
{
    textsize = (float)h;
	return;
}

- (void) setPaletteEntryWithIndex:(long)index red:(short)r green:(short)g blue:(short)b
{
    printf("Graphic: setPaletteEntryWithIndex %ld (%d,%d,%d)\n",index,r,g,b);
	return;
}

- (void) setNewPaletteEntryWithIndex:(long)index withCount:(long)count red:(short)r green:(short)g blue:(short)b
{
    printf("Graphic: setNewPaletteEntryWithIndex %ld count %ld (%d,%d,%d)\n",index,count,r,g,b);
    return;
}

- (void) getPaletteEntryWithIndex:(long)index red:(short *)r green:(short *)g blue:(short *)b
{
    printf("Graphic: getPaletteEntryWithIndex %ld = (%d,%d,%d)\n",index,*r,*g,*b);
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
