/****************************************************************************/
/*																			*/
/* File:	  MShellWindow.m												*/
/*																			*/
/* Purpose:   Implementation of the MShellWindow class.						*/
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


#import "MShellWindow.h"

/* low module */
#include "misc.h"
#include "defaults.h"
#include "general.h"

/* interface includes */
#include "devices.h"
#include "cmdint.h"


unsigned int ugTextLength;

@implementation MShellWindow

- (id)initWithContentRect:(NSRect)contentRect
                styleMask:(unsigned int)styleMask
                backing:(NSBackingStoreType)backingType
                defer:(BOOL)flag
/*
 The init method is overwritten here to set up the NSCharacterSet called finalCharSet,
 which is used later to determine which keyDown events can simply be passed to the
 TextView - all other keyDown events have to be treated differently.
 */
{
    
	if ( self = [super initWithContentRect:contentRect styleMask:styleMask
                             backing:backingType defer:flag] )
    {
        NSMutableCharacterSet *workingSet = [[NSCharacterSet alphanumericCharacterSet] mutableCopy];
        [workingSet addCharactersInString:@";:,."];
        finalCharSet = [workingSet copy];
        [workingSet release];
        commandString = [NSMutableString stringWithString:@""];
    }
    
	return self;
}

- (void)keyDown:(NSEvent *)theEvent
{
    char newline = '\n';
	NSString *theChars = [theEvent characters];		// get Character from theEvent

	[super keyDown:theEvent];
    printf("keyDown>>%c..%i<<\n",
           (char)[theChars characterAtIndex:0],(int)[theChars characterAtIndex:0]);
	return;
    
    if ([theChars length] > 1)
        printf("MShellWindow, keyDown: theChar length is %d\n",[theChars length]);
    
    if ( [finalCharSet characterIsMember:[theChars characterAtIndex:0]] )	{
        [super keyDown:theEvent];
        }
    else	{
        printf ("not a character event\n");
        }
    
    // append the key to the commandString
    [commandString appendString:theChars];

    if ([theChars characterAtIndex:0] == '\n')	{
		/* Process command */
        printf("\\n gepresst\n");
        }

    return;
}

@end

