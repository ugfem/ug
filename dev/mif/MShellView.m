/****************************************************************************/
/*																			*/
/* File:	  MShellView.m													*/
/*																			*/
/* Purpose:   Implementation for the MShellView class.						*/
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

/* based on CLIView.m by Philippe Mougin */

/*

The organisation of the view hierarchy:

      ----- MShellView -----------------------------------
      |                                                  |
      | ----- NSScrollView ----------------------------  |
      | |                                             |  |
      | |  ----- MShellTextView --------------------  |  |
      | |  | prompt>                               |  |  |
      | |  |                                       |  |  |
      | |  |                                       |  |  |
      | |  |                                       |  |  |
      | |  |                                       |  |  |
      | |  |                                       |  |  |
      | |  -----------------------------------------  |  |
      | |                                             |  |
      | -----------------------------------------------  |
      |                                                  |
      ----------------------------------------------------

A MShellView has one subview: a NSScrollView.

This NSScrollView has a MShellTextView as document view.

The MShellTextView is the view that displays the prompt, receive the keyboard events from the user, display
the commands entered by the user and the results of those commands etc.  

*/


#import "MShellView.h"
#import "MShellTextView.h"
#import "MShell.h"

@interface MShellView(MShellViewPrivate)
- (MShellTextView *)textView;
@end

extern MShell *theUGshell;


@implementation MShellView

- (float)fontSize
{ return [[[self textView] font] pointSize]; }

- (id) _init  // construction and configuration of the view hierarchy
{
    NSSize contentSize;
    theScrollView =[[[NSScrollView alloc] initWithFrame:[self bounds]] autorelease];
    contentSize = [theScrollView contentSize];

    [theScrollView setBorderType:NSNoBorder];
    [theScrollView setHasVerticalScroller:YES];
    [theScrollView setHasHorizontalScroller:YES];
    [theScrollView setAutoresizingMask:(NSViewWidthSizable | NSViewHeightSizable)]; 

    /* Create the MShellTextView */
    theTextView = [[[MShellTextView alloc] initWithFrame:NSMakeRect(0, 0,
           [theScrollView contentSize].width, [theScrollView contentSize].height)] autorelease];
    // 2000.0 instead of [theScrollView contentSize].width?

    [theTextView setMinSize:(NSSize){50.0, 50.0}];
    [theTextView setMaxSize:(NSSize){1e7, 1e7}];
    [theTextView setVerticallyResizable:YES];
    [theTextView setHorizontallyResizable:YES];
    [theTextView setAutoresizingMask:NSViewWidthSizable];
    [[theTextView textContainer] setWidthTracksTextView:NO];

    [theTextView setEditable:YES];
    [theTextView setSelectable:YES];
    [[NSFont userFixedPitchFontOfSize:10.0] set];
    //[theTextView setFont:[NSFont userFixedPitchFontOfSize:10.0]];
    [self addSubview:theScrollView];
    [theScrollView setDocumentView:theTextView];

    /*
     * Get notification if text input ended:
     * the MShellWindow object will from now on get a notification about
     * all textDidEndEditing events in theTextView.
     
    [[NSNotificationCenter defaultCenter]
        addObserver:self
        selector:@selector(textDidEndEditing:)
        name:NSTextDidEndEditingNotification
        object:theTextView];
     */

    return self;
}

- (id)initWithCoder:(NSCoder *)coder
{
  self = [super initWithCoder:coder];
  [self _init];
  [self setAutoresizesSubviews:YES];
  return self;
}

- (id)initWithFrame:(NSRect)frameRect
{
  self = [super initWithFrame:frameRect];
  if ( self!=nil )
  {
    [self _init];
    [self setAutoresizesSubviews:YES];    
    return self;
  }
  return nil;
}

- (MShellTextView *)textView
{
  return theTextView;
}

- (NSScrollView *)scrollView
{
  return theScrollView;
}

@end
