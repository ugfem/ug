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

@interface MShellView(MShellViewPrivate)
- (MShellTextView *)shellView;
@end

@implementation MShellView

- (float)fontSize
{ return [[[self shellView] font] pointSize]; }

- (id) _init  // construction and configuration of the view hierarchy
{
    NSScrollView *scrollview =[[[NSScrollView alloc] initWithFrame:[self bounds]] autorelease];
    NSSize contentSize = [scrollview contentSize];
    MShellTextView *shellView;

    [scrollview setBorderType:NSNoBorder];
    [scrollview setHasVerticalScroller:YES];
    [scrollview setHasHorizontalScroller:NO];
    [scrollview setAutoresizingMask:NSViewWidthSizable | NSViewHeightSizable]; 

    shellView = [[[MShellTextView alloc] initWithFrame:NSMakeRect(0, 0,[scrollview contentSize].width, [scrollview contentSize].height)] autorelease];
    [shellView setMinSize:(NSSize){0.0, contentSize.height}];
    [shellView setMaxSize:(NSSize){1e7, 1e7}];
    [shellView setVerticallyResizable:YES];
    [shellView setHorizontallyResizable:NO];
    [shellView setAutoresizingMask:NSViewWidthSizable ];
    [[shellView textContainer] setWidthTracksTextView:YES];

    [scrollview setDocumentView:shellView];
    [self addSubview:scrollview];
    [self setFontSize:[[NSUserDefaults standardUserDefaults] integerForKey:@"fontSize"]];
    
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
  if (self = [super initWithFrame:frameRect])
  {
    [self _init];
    [self setAutoresizesSubviews:YES];    
    return self;
  }
  return nil;
}


// Private

- (MShellTextView *)shellView
{
  return [[[self subviews] objectAtIndex:0] documentView];
}

@end
