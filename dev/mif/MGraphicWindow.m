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

    viewHeight = contentRect.size.height;
    
	// set default colors
    currentColor = [NSColor blackColor];
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
    //[statusText setBackgroundColor:[self backgroundColor]];
    [statusText setTextColor:textColor];
    [statusText setEditable:YES];
    [statusText insertText:@""];
    [statusText setEditable:NO];    
    [statusText display];
        
	[self display];
    [self orderFront:self];
	[self setReleasedWhenClosed:YES];
    [self makeFirstResponder:theView];

    /* Let the graphics window accept mouse moved events (for tracking of
       x and y coordinates).  Default is NO for NSWindow objects.*/
    [self setAcceptsMouseMovedEvents:YES];

    textAttributes = [NSDictionary dictionaryWithObjectsAndKeys:
                            [NSColor blackColor],NSForegroundColorAttributeName,
                            [NSColor whiteColor],NSBackgroundColorAttributeName,
                            [NSFont boldSystemFontOfSize:10.0], NSFontAttributeName,
                            nil];
    textsize = 10.0;

    path = [NSBezierPath bezierPath];
    
	updateDate = [NSDate date];
	timeBetweenUpdates = 2;		// wait this many seconds til update

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
    [self makeKeyWindow];
	[self makeKeyAndOrderFront:self];
	//[self orderFont:self];
}

- (void) updateOutput
{
    printf("Update Output\n");
}


- (void) moveToPoint:(SHORT_POINT)point
{
	moveto_x = point.x;
	moveto_y = viewHeight-point.y;
}

- (void) drawLineTo:(SHORT_POINT)point
{
	NSTimeInterval d;
	
	[theView lockFocus];
	[currentColor set];
	[NSBezierPath setDefaultLineWidth:linewidth];
	[NSBezierPath strokeLineFromPoint:NSMakePoint(moveto_x, moveto_y) 
							  toPoint:NSMakePoint(point.x, viewHeight-point.y)];
	
	/* use this to update every 'timeBetweenUpdates' seconds */
	/*
	d = [[NSDate date] timeIntervalSinceDate:updateDate];
	
	if ( fabs(d) > (float)timeBetweenUpdates ) 	{	
		// if the update date is earlier than the current date update graphics
		printf("update %g\n",(float)d);
		[self flushWindow];
		[updateDate release];
		updateDate = [NSDate date];
	}
	*/
	
	if ( ups>500 )	{	// don't do this very often
		[self flushWindow];
		ups = 0;
	}
	ups++;
	
	[theView unlockFocus];
}


- (void) drawPolyLine:(SHORT_POINT *)point noOfPoints:(INT)n
{
	int i;
	
	if (n<2) return;

	[theView lockFocus];
	[currentColor set];
	[NSBezierPath setDefaultLineWidth:0.0];
	for (i=0; i<n-1; i++)
		[NSBezierPath strokeLineFromPoint:NSMakePoint(point[i].x, viewHeight-point[i].y) 
								  toPoint:NSMakePoint(point[i].x, viewHeight-point[i].y)];
	[theView unlockFocus];
}


- (void) drawInversePolyLine:(SHORT_POINT *)points noOfPoints:(INT)n
{
	int i;
    // EXPERIMENTAL! (TODO: Is this still used?)
    printf("drawInversePolyLine\n");
    [theView lockFocus];
    [[NSColor blueColor] set];
    [path moveToPoint:(NSPoint){points[0].x, viewHeight-points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, viewHeight-points[i].y}];
    [path closePath];
    [path stroke];
    [path removeAllPoints];	// There is no 'reset' method, though the manual says so...
    [currentColor set];
    [theView unlockFocus];
	return;
}


- (void) drawPolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
	int i;    
	if (n<3) return;

    [theView lockFocus];
	[currentColor set];
    [path moveToPoint:(NSPoint){points[0].x, viewHeight-points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, viewHeight-points[i].y}];
    [path closePath];
    [path fill];
    [path removeAllPoints];
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
    [path moveToPoint:(NSPoint){points[0].x, viewHeight-points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, viewHeight-points[i].y}];
    [path closePath];
    [path fill];
    [path removeAllPoints];	// There is no 'reset' method, though the manual says so...
    [currentColor set];
    [theView unlockFocus];
    return;
}


- (void) erasePolygon:(SHORT_POINT *)points noOfPoints:(INT)n
{
    int i;
    [theView lockFocus];
    /*PSmoveto(points[0].x, points[0].y);
    for (i=1; i<n; i++) PSlineto (points[i].x, points[i].y);
    PSclosepath();
    PSfill();
    PSstroke();*/
    [backgroundColor set];
    [path moveToPoint:(NSPoint){points[0].x, viewHeight-points[0].y}];
    for (i=1; i<n; i++) [path lineToPoint:(NSPoint){points[i].x, viewHeight-points[i].y}];
    [path closePath];
    [path fill];
    [path removeAllPoints];	// There is no 'reset' method, though the manual says so...
    [theView unlockFocus];
    [currentColor set];
	return;
}


- (void) setMarkerSize:(short)width
{
    printf("setMarkerSize\n");
    marker_size = width;
    return;
}

- (void) setMarker:(short)n
{
    printf("setMarker\n");
    marker_id = n;
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
    [textColor set];
    [currentColor set];
    [theView lockFocus];
    [s drawAtPoint:(NSPoint){moveto_x,moveto_y}];
    [theView unlockFocus];
	return;
}

- (void) drawCenteredText:(const char*) text atPoint:(SHORT_POINT)point mode:(INT)m
{
    float w,h;
    NSAttributedString *s = [[NSAttributedString alloc] initWithString:[NSString stringWithCString:text]
                                   attributes:textAttributes];
    [theView lockFocus];
    w = [s size].width;
    h = [s size].height;
    [currentColor set];
    [s drawInRect:NSMakeRect(point.x-w/2,viewHeight-point.y-h/2,w,h)];
    //[s drawInRect:NSMakeRect(point.x,viewHeight-point.y,w,h)];
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

- (void) setColorRed:(float)r green:(float)g blue:(float)b
{
    [theView lockFocus];
	[currentColor release];
    currentColor = [NSColor colorWithDeviceRed:r green:g blue:b alpha:1.0];
	[currentColor set];
    [theView unlockFocus];
}

- (void) setLineWidth:(short)width
{
	linewidth = (float)width;
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
	[self flushWindow];
	[theView unlockFocus];
	return;
}

@end
