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
    // Ensure that only one instance exists
    if (singleInstance == nil)	{
        singleInstance = [[self alloc] init];
        if ([singleInstance setUp]==NO) return nil;
        }
    
    return singleInstance;
}

- (BOOL) setUp
{
	/*
	 * Open up a Window
	 */
    int TermWinX, TermWinY, TermWinW, TermWinH;
    char buffer[256];
    
    /* read default size of shell window */
    if (GetDefaultValue(DEFAULTSFILENAME,"TermWinH",buffer)==0)
        sscanf(buffer," %d ",&TermWinX);
    if (GetDefaultValue(DEFAULTSFILENAME,"TermWinV",buffer)==0)
        sscanf(buffer," %d ",&TermWinY);
    if (GetDefaultValue(DEFAULTSFILENAME,"TermWinDH",buffer)==0)
        sscanf(buffer," %d ",&TermWinW);
    if (GetDefaultValue(DEFAULTSFILENAME,"TermWinDV",buffer)==0)
        sscanf(buffer," %d ",&TermWinH);
    theWindow = [[MShellWindow alloc] initWithContentRect:NSMakeRect(TermWinX, TermWinY, TermWinW, TermWinH)
					styleMask:(NSTitledWindowMask | NSResizableWindowMask | NSMiniaturizableWindowMask)
                    backing:NSBackingStoreBuffered
					defer:NO];
	[theWindow setTitle:@"UG shell"];
	[NSApp setDelegate:theWindow];

	/*
     *  Put a MShellView into the window
     */
    theShellView = [[MShellView alloc] initWithFrame:[[theWindow contentView] bounds]];
	theTextView = [theShellView textView];
        
    [[theWindow contentView] addSubview:[theShellView scrollView]];
    [theWindow makeKeyAndOrderFront:self];
    [theWindow makeFirstResponder:theTextView];
	//[theWindow setDelegate:theTextView];
	    
    /* allocate input line buffer */
    if ((inpLine=(char *)malloc(cmdintbufsize))==NULL)
    {
        PrintErrorMessage('F',"MShell setUp method","could not allocate inpLine buffer");
        return NO;
    }
    inpLine[0] = (char) 0;

	[theWindow display];
	
    return YES;
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
	//printf ("applicationDidFinishLaunching called\n");
    //[self appendToText:"applicationDidFinishLaunching called\n"];
}

- (void) appendPrompt
{
	[theTextView appendPrompt];
}

- (void) appendToText:(const char *)val
{
    [theTextView appendText:[NSString stringWithCString:val]];
    return;
}

- (NSPoint) getMouseLocation
{
    return [theWindow mouseLocationOutsideOfEventStream];
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
    fontsize = fsize;
    printf ("6 Current font is %s", [[[theTextView font] displayName] lossyCString]);
    printf (" at size %f\n", [[theTextView font] pointSize]);
    //[[NSFont userFixedPitchFontOfSize:24.0] set];
    [theTextView setFont:[NSFont userFixedPitchFontOfSize:fontsize]];
    //[theTextView setTextColor:[NSColor redColor]];
    //[theTextView textDidChange];
    printf ("7 Current font is %s", [[[theTextView font] displayName] lossyCString]);
    printf (" at size %f\n", [[theTextView font] pointSize]);

    //[[theTextView text] beginEditing];
    //[myTextView setFont:[NSFont fontWithName:@"Helvetica-Oblique"
//    size:12.0] range:changeCharRange];
 //   [theTextStore endEditing];
  //  [theTextView didChangeText];

}


/*- (BOOL)theTextView:(NSTextView *)aTextView shouldChangeTextInRange:(NSRange)affectedCharRange replacementString:(NSString *)replacementString
{
	printf("===range %d %d===\n",
        affectedCharRange.location,affectedCharRange.length);
    return YES;
}


- (NSRange)theTextView:(NSTextView *)aTextView willChangeSelectionFromCharacterRange:(NSRange)oldSelectedCharRange toCharacterRange:(NSRange)newSelectedCharRange
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

- (BOOL) interpretCommand:(NSString*)command
{
    int error;

    inpLine = strcpy(inpLine,[command lossyCString]);
    error = InterpretCommand(inpLine);
    if ( error )	{
        if (error==QUITCODE)	{
            [NSApp stop:nil];
            return NO;
        }
        else	{
            [theTextView highlightError];
		}
    }

    return YES;
}

- (MShellWindow*) window
{
	return theWindow;
}

- (MShellTextView*)	textview
{
    return theTextView;
}
@end

