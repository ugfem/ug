/****************************************************************************/
/*																			*/
/* File:	  MShell.m														*/
/*																			*/
/* Purpose:   Implementation of the MShell class.							*/
/*            See information below.										*/
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

/* UG includes */
#include "compiler.h"
#include "misc.h"
#include "defaults.h"
#include "general.h"
#include "devices.h"
#include "cmdint.h"

static MShell	*singleInstance = nil;

@implementation MShell

/*
	Create a unique instance of MShell.  If it doesn't exist yet, set it up.  If it
 	already exists return the already existing copy.
 */
+ (MShell *) instantiate;
{
    printf(ARCHNAME);
    printf(" instantiate...\n");

    // Ensure that only one instance exists
    if (singleInstance == nil)	{
        singleInstance = [[self alloc] init];
        if ([singleInstance setUp]==NO) return nil;
        }
    
    return singleInstance;
}

- (BOOL) setUp
{
	NSSize	contentSize;

    printf(ARCHNAME);
    printf(" setUp...\n");
    
	/*
	 * Open up a Window
	 */
	theWindow = [[NSWindow alloc] initWithContentRect:NSMakeRect(100.0, 350.0, 600.0, 300.0)
					styleMask:(NSTitledWindowMask | NSResizableWindowMask | NSMiniaturizableWindowMask)
					backing:NSBackingStoreRetained
					defer:YES];
	[theWindow setTitle:@"UG shell"];
	[NSApp setDelegate:theWindow];
	
    
	/* 
	 * Put a scrollView into theWindow
	 */
	theScrollView = [[NSScrollView alloc]
    	initWithFrame:[[theWindow contentView] bounds]];

	contentSize = [theScrollView contentSize];
	[theScrollView setBorderType:NSNoBorder];
	[theScrollView setHasVerticalScroller:YES];
	[theScrollView setHasHorizontalScroller:YES];
	[theScrollView setAutoresizingMask:(NSViewWidthSizable | NSViewHeightSizable)];

    
    /* Create the NSTextView and add it to the window */
	theTextView = [[MShellTextView alloc] initWithFrame:NSMakeRect(0, 0,
    	2000, [theScrollView contentSize].height)];
	[theTextView setMinSize:(NSSize){50.0, 50.0}];
	[theTextView setMaxSize:(NSSize){1e7, 1e7}];
	[theTextView setVerticallyResizable:YES];
	[theTextView setHorizontallyResizable:YES];
	[theTextView setAutoresizingMask:NSViewWidthSizable];
	
	[[theTextView textContainer] setWidthTracksTextView:YES];

	[theScrollView setDocumentView:theTextView];
    [[theWindow contentView] addSubview:theScrollView];
    [theWindow makeKeyAndOrderFront:nil];
    [theWindow makeFirstResponder:theTextView];
    [theWindow setDelegate:theTextView];
	    
    [theTextView setEditable:YES];
    [theTextView setSelectable:YES];
    [theTextView setFieldEditor:YES];
    [[NSFont userFixedPitchFontOfSize:10.0] set];

	/* 
	 * Get notification if text input ended: 
	 * the MShellWindow object will from now on get a notification about
     * all textDidEndEditing events in theTextView.
     */
	[[NSNotificationCenter defaultCenter]
		addObserver:self
     	selector:@selector(textDidEndEditing:)
		name:NSTextDidEndEditingNotification
       	object:theTextView];
    
	/* alloc input line buffer */
	if ((inpLine=(char *)malloc(cmdintbufsize))==NULL)
	{
		PrintErrorMessage('F',"MShell setUp method","could not allocate inpLine buffer");
		return NO;
	}
	inpLine[0] = (char) 0;

    /*[[NSNotificationCenter defaultCenter]
		addObserver:self
     	selector:@selector(textView:shouldChangeTextInRange:replacementString:)
		name:NSTextDidChangeNotification
       	object:theTextView];
    
	[[NSNotificationCenter defaultCenter]
		addObserver:self
     	selector:@selector(textView:willChangeSelectionFromCharacterRange:toCharacterRange:)
		name:NSTextViewDidChangeSelectionNotification
       	object:theTextView];*/

    theCommandList = [NSMutableArray arrayWithCapacity:100];
    
    printf("MShellWindow: setUp done\n");

	/* Don't display it, this is done later after the correct size is set */
	/* [theWindow display]; */
	
	return YES;
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
	printf ("applicationDidFinishLaunching called\n");
    [self appendToText:"applicationDidFinishLaunching called\n"];
}

- (void)appendToText:(const char *)val
{
    [[theTextView textStorage]
        replaceCharactersInRange:NSMakeRange([[theTextView textStorage] length], 0)
        withString:[NSString stringWithCString:val]];
    [theTextView scrollRangeToVisible:NSMakeRange([[theTextView textStorage] length], 0)];
    [theTextView didChangeText];
    return;
}

- (NSPoint) getMouseLocation
{
    return [theWindow mouseLocationOutsideOfEventStream];
}

- (void)textDidEndEditing:(NSNotification *)notification
{
    int error;
    char tab[2], newline[2];
    NSRange	lineRange;
    NSString *lastLine;
    NSCharacterSet *newlineSet;
    unsigned whyEnd = [[[notification userInfo] objectForKey:@"NSTextMovement"] unsignedIntValue];

    tab[0]='\t';		tab[1]='\0';
    newline[0]='\n';	newline[1]='\0';
    newlineSet = [NSCharacterSet characterSetWithCharactersInString:[NSString stringWithCString:newline]];
    
    // Unscroll the previous text.
    [theTextView scrollRangeToVisible:NSMakeRange([[theTextView textStorage] length], 0)];

    // If the text editing was ended because of a tab then append the tab.
    // TODO: Replace this with command completion.
    if (whyEnd == NSTabTextMovement)	{
        [self appendToText:tab];
        [theWindow makeFirstResponder:theTextView];
    }
    // If the text was ended by a return, then get the last line (the command line) and
    // process it.
    else if ( whyEnd== NSReturnTextMovement )	{
        // Find the last newline character
        lineRange = [[[theTextView textStorage] 
        					string] 
        					rangeOfCharacterFromSet:newlineSet options:NSBackwardsSearch];
        // Get the command line (assume prompt-length of 2)
        lastLine = [[[theTextView textStorage] string]
			substringWithRange:NSMakeRange(lineRange.location+2,
                                  [[theTextView textStorage] length]-lineRange.location-2)];
        // now append the newline character
        [self appendToText:newline];
        // printf("--%s--\n",[lastLine lossyCString]);
		
        // Add the command to the history of commands
        [theCommandList insertObject:lastLine atIndex:(commandCount%MAX_BUFFERED_COMMANDS)];
        commandCount++;
		currCommandNumber = commandCount;
        
        // Process the command
        inpLine = strcpy(inpLine,[lastLine lossyCString]);
        if (error = InterpretCommand(inpLine))	{
        	if (error==QUITCODE)	{
                [NSApp stop:nil];
                return;
            	}
        }
        
        // Set everything up to get the next command
        [theWindow makeFirstResponder:theTextView];
        [self appendToText:PROMPT];
    }
    else if (whyEnd == NSUpTextMovement)	{
        // Find the last newline character
        lineRange = [[[theTextView textStorage] string]
			rangeOfCharacterFromSet:newlineSet options:NSBackwardsSearch];
		lineRange.length += strlen(PROMPT);
        // Replace the current text on the command line if selected command exists
        --currCommandNumber;
        if (currCommandNumber>=0)
			[[theTextView textStorage]
				replaceCharactersInRange:lineRange
            	withString:[theCommandList objectAtIndex:currCommandNumber]];
    }
    else if (whyEnd == NSDownTextMovement)	{
        // Find the last newline character
        lineRange = [[[theTextView textStorage] string]
			rangeOfCharacterFromSet:newlineSet options:NSBackwardsSearch];
        // Replace the current text on the command line
        ++currCommandNumber;
        if (currCommandNumber < commandCount)
			[[theTextView textStorage]
				replaceCharactersInRange:lineRange
            	withString:[theCommandList objectAtIndex:currCommandNumber]];
    }
        
    return;
}

- (void)textDidChange:(NSNotification *)notification
{
	return;
}


- (void)	setScrollback:(int)sb
{

}

- (void)	setCharactersPerLine:(int)cpl
{

}

- (void)	setFontSize:(float)fsize
{
	[[NSFont userFixedPitchFontOfSize:fsize] set];
}


/*- (BOOL)textView:(NSTextView *)aTextView shouldChangeTextInRange:(NSRange)affectedCharRange replacementString:(NSString *)replacementString
{
	printf("===range %d %d===\n",
        affectedCharRange.location,affectedCharRange.length);
    return YES;
}


- (NSRange)textView:(NSTextView *)aTextView willChangeSelectionFromCharacterRange:(NSRange)oldSelectedCharRange toCharacterRange:(NSRange)newSelectedCharRange
{
	printf("---old %d %d new %d %d---\n",
        oldSelectedCharRange.location,oldSelectedCharRange.length, 
        newSelectedCharRange.location,newSelectedCharRange.length);
    return oldSelectedCharRange;
}*/

- (BOOL) executeFile:(NSString *)fileName
{
    char command[160];
    char fname[128];
	[fileName getCString:fname maxLength:128];
    
    sprintf(command,"execute %s",fname);
    InterpretCommand(command);
    return YES;
}

- (MShellWindow*) window
{
	return theWindow;
}

@end

