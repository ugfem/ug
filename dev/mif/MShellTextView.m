/****************************************************************************/
/*																			*/
/* File:	  MShellTextView.												*/
/*																			*/
/* Purpose:   Implementation of the MShellTextView class.					*/
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
#import "MShellTextView.h"

#include "cmdint.h"


#define RETURN    @"\015"
#define BACKSPACE @"\010"


extern MShell *theUGshell;

@implementation MShellTextView

- (id) initWithFrame:(NSRect)frameRect
{
    unsigned int i;
    
    [super initWithFrame:frameRect];
    if ( self==nil ) return nil;
    
    theCommandList = [NSMutableArray arrayWithCapacity:MAX_BUFFERED_COMMANDS];
    for ( i=0; i<MAX_BUFFERED_COMMANDS; i++ )
        [theCommandList insertObject: @"" atIndex:i];

    currentLine = [NSMutableString string];
    
    return self;
}

- (void) appendPrompt
{
    [[self textStorage]
        replaceCharactersInRange:NSMakeRange([[self textStorage] length], 0)
        withString:[NSString stringWithCString:PROMPT]];

    startPosition = [[self textStorage] length];
    [self scrollRangeToVisible:NSMakeRange(startPosition, 0)];
    [self didChangeText];
    return;
}

- (void) keyDown:(NSEvent *)theEvent
{
	if ([theEvent type] != NSKeyDown)	{
        [super keyDown:theEvent];
        return;
    }
    
    if ([self selectedRange].location < startPosition)
    {
        [self moveToEndOfDocument:self]; 	// ...not implemented yet
        [self setSelectedRange:NSMakeRange([[self string] length],0)];
        [self scrollRangeToVisible:[self selectedRange]];  
    }
	/*
     *	Depending on which key was pressed we do:
     *	1. Return -> execute command
     *
     *
     */
    if ( [[theEvent characters] isEqualToString:RETURN] )
    {
        NSString *command = [[self string] substringFromIndex:startPosition];

        // Add command to command list.
        if ( [command length] > 0 )	{
            if (commandCount>0 )	{
                if ( ![command isEqualToString:[theCommandList objectAtIndex:(commandCount%MAX_BUFFERED_COMMANDS)-1]] )	{
                    [theCommandList replaceObjectAtIndex:commandCount%MAX_BUFFERED_COMMANDS  withObject:command];
                    commandCount++;
                }
            }
            else {
                [theCommandList replaceObjectAtIndex:commandCount%MAX_BUFFERED_COMMANDS  withObject:command];
                commandCount++;
            }
            currCommandNumber=commandCount;
        }
        else {
            [super keyDown:theEvent];
            return;
        }
        [self moveToEndOfDocument:self];
        [self setSelectedRange:NSMakeRange([[self string] length],0)];
        [super keyDown:theEvent];
        [theUGshell interpretCommand:command];
        [self appendPrompt];
        [self scrollRangeToVisible:[self selectedRange]];
        startPosition = [[self string] length];
        [currentLine setString:@""];
    }
	
    else if ([[theEvent characters] characterAtIndex:0] == NSLeftArrowFunctionKey)
    {
      if ([self selectedRange].location > startPosition)
          [super keyDown:theEvent];
    }
    
    else if ( [[theEvent characters] characterAtIndex:0] == NSUpArrowFunctionKey )
    {
        NSString *command = [[self string] substringFromIndex:startPosition];
		
		// if something was typed in save it
		if ( [command length]>0 && commandCount==currCommandNumber && currCommandNumber>0)	{
            if ( ![command isEqual:currentLine] ) {
                [currentLine setString:command];
            }
		}

        //if ( commandCount==currCommandNumber ) currCommandNumber--;
        if ( currCommandNumber>0 )	{
            currCommandNumber--;
            [self replaceCurrentCommandWith:[theCommandList objectAtIndex:currCommandNumber%MAX_BUFFERED_COMMANDS]];
        }
    }

    else if ( [[theEvent characters] characterAtIndex:0] == NSDownArrowFunctionKey )
    {
        currCommandNumber++;
        if ( currCommandNumber<commandCount )	{
            [self replaceCurrentCommandWith:[theCommandList objectAtIndex:(currCommandNumber%MAX_BUFFERED_COMMANDS)]];
        }
		else	{
            [self replaceCurrentCommandWith:currentLine];
            currCommandNumber = commandCount;
        }
    }

    else
        [super keyDown:theEvent];

	return;
}

- (void) replaceCurrentCommandWith:(NSString *)command
{
    [self setSelectedRange:NSMakeRange(startPosition,[[self string] length])];
    [self insertText:command];
    [self moveToEndOfDocument:self];
    [self scrollRangeToVisible:[self selectedRange]];

    return;
}

@end
