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

    commandCount = currCommandNumber = 0;
    theCommandList = [NSMutableArray arrayWithCapacity:MAX_BUFFERED_COMMANDS];
    for ( i=0; i<MAX_BUFFERED_COMMANDS; i++ )
        [theCommandList insertObject: @"" atIndex:i];

    currentLine = [NSMutableString string];
	
    [self setRichText:YES];
	[self setFont:[NSFont userFixedPitchFontOfSize:10.0]];
        
    /*[self setFont:[NSFont userFixedPitchFontOfSize:10]
          range:NSMakeRange(0,[[self string] length])];

    printf ("4 Desired font is %s, ", [[[NSFont userFixedPitchFontOfSize:18] displayName] lossyCString]);
    printf ("Current font is %s ", [[[self font] displayName] lossyCString]);
    printf ("at size %f\n", [[self font] pointSize]);*/

    textAttributes = [NSDictionary dictionaryWithObjectsAndKeys:
        					[NSColor blackColor],NSForegroundColorAttributeName,
                            [NSColor whiteColor],NSBackgroundColorAttributeName,
        					[NSFont userFixedPitchFontOfSize:10.0], NSFontAttributeName,
                            nil];
   	
    errorAttributes = [NSDictionary dictionaryWithObjectsAndKeys:
                            [NSColor whiteColor],NSForegroundColorAttributeName,
                            [NSColor blackColor],NSBackgroundColorAttributeName,
        					[NSFont userFixedPitchFontOfSize:10.0], NSFontAttributeName,
                            nil];
   	
    promptAttributes = [NSDictionary dictionaryWithObjectsAndKeys:
        					[NSColor blackColor],NSForegroundColorAttributeName,
                            [NSColor whiteColor],NSBackgroundColorAttributeName,
                            [NSFont userFixedPitchFontOfSize:10.0], NSFontAttributeName,
                            nil];

    cmdlineAttributes = [NSDictionary dictionaryWithObjectsAndKeys:
                            [NSColor blackColor],NSForegroundColorAttributeName,
                            [NSColor whiteColor],NSBackgroundColorAttributeName,
                            [NSFont userFixedPitchFontOfSize:10.0], NSFontAttributeName,
                            nil];

    theTextStore = [self textStorage];

    if ([self shouldChangeTextInRange:NSMakeRange(0,[[self string] length])
                    replacementString:nil])
    {
        [theTextStore beginEditing];
        [theTextStore addAttributes:textAttributes range:NSMakeRange(0,[[self string] length])];
        [theTextStore fixAttributesInRange:NSMakeRange(0,[[self string] length])];
        [theTextStore endEditing];
        [self didChangeText];
	}
    
    printf ("5 Current font is %s ", [[[self font] displayName] lossyCString]);
    printf ("at size %f\n", [[self font] pointSize]);

    return self;
}

- (void) appendPrompt
{
	int plen = strlen(PROMPT);
    
    [theTextStore
        replaceCharactersInRange:NSMakeRange([theTextStore length], 0)
        withString:[NSString stringWithCString:PROMPT]];

    if ([self shouldChangeTextInRange:NSMakeRange([theTextStore length]-plen, plen)
                    replacementString:nil])
    {
        [theTextStore beginEditing];
        [theTextStore addAttributes:promptAttributes range:NSMakeRange([theTextStore length]-plen, plen-1)];
        [theTextStore addAttributes:cmdlineAttributes range:NSMakeRange([theTextStore length]-1, 1)];
        [theTextStore endEditing];
        [self didChangeText];
    }

    startPosition = [theTextStore length];
    [self scrollRangeToVisible:NSMakeRange(startPosition, 0)];

	return;
}

- (void) appendText:(NSString *)str
{
    static int i;

    [self moveToEndOfDocument:self];
    [self insertText:str];
    [self didChangeText];

    if (++i%32==0 ) [self displayIfNeeded];
    
    return;
}

- (void) keyDown:(NSEvent *)theEvent
{
    unichar key;
    
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

    key = [[theEvent characters] characterAtIndex:0];

    if ( (key==NSNewlineCharacter || key==NSEnterCharacter || key==NSCarriageReturnCharacter)
         && !([theEvent modifierFlags] & NSCommandKeyMask))
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
	
    else if ( key==NSLeftArrowFunctionKey
              && !([theEvent modifierFlags] & NSCommandKeyMask))
   {
      	if ([self selectedRange].location > startPosition)
       		[super keyDown:theEvent];
    }
    
    else if ( key==NSUpArrowFunctionKey
              && !([theEvent modifierFlags] & NSCommandKeyMask))
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

    else if ( key==NSDownArrowFunctionKey 
              && !([theEvent modifierFlags] & NSCommandKeyMask))
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

    else if ( (key==NSBackspaceCharacter || key==NSBackTabCharacter || key==NSDeleteCharacter)
              && !([theEvent modifierFlags] & NSCommandKeyMask))
    {
        if ([self selectedRange].location > startPosition)	{
            [super keyDown:theEvent];
       	}
	}

    else if ( key==NSTabCharacter
              && !([theEvent modifierFlags] & NSCommandKeyMask))	
	{
		// Command completion based on command history
     	int k;
        NSString *command = [[self string] substringFromIndex:startPosition];
		NSMutableString *str = [NSMutableString string];
        for ( k=currCommandNumber-1; k>=0; k-- )	{
            [str setString:[theCommandList objectAtIndex:(k%MAX_BUFFERED_COMMANDS)]];
            if ( [str length]>[command length] &&
                		[command isEqual:[str substringToIndex:[command length]]] )
                [self replaceCurrentCommandWith:str];
        }
    }

    else if ( ([theEvent modifierFlags] & NSCommandKeyMask) )	{
        // Command key
        char keyCommand[256];
        int k;
        NSString *command;
        //DoCmdKey (key, keyCommand);
        printf ("command key %c = %s\n",key, "[not implemented]" );
        [super keyDown:theEvent];
        return;
        /*command = [NSString withCString:keyCommand];
        [self replaceCurrentCommandWith:[NSString withCString:keyCommand]];

        [self moveToEndOfDocument:self];
        [self setSelectedRange:NSMakeRange([[self string] length],0)];

        [super keyDown:theEvent];
        [theUGshell interpretCommand:command];
        [self appendPrompt];
        [self scrollRangeToVisible:[self selectedRange]];
        startPosition = [[self string] length];
        [currentLine setString:@""];*/
    }
    
	else	{
        [super keyDown:theEvent];
    }

	return;
}

- (void) replaceCurrentCommandWith:(NSString *)command
{
    [self setSelectedRange:NSMakeRange(startPosition,[[self string] length]-startPosition)];
    [theTextStore addAttributes:cmdlineAttributes range:NSMakeRange(startPosition,[[self string] length]-startPosition)];
    [self insertText:command];
    [self moveToEndOfDocument:self];
    [self scrollRangeToVisible:[self selectedRange]];

    return;
}

- (void)highlightError
{
	NSRange range=[[theTextStore string] rangeOfString:@"\n" options:NSLiteralSearch
				range:NSMakeRange(startPosition,[[self string] length]-startPosition)];

	if ( range.length>0 )
        if ([self shouldChangeTextInRange:NSMakeRange(startPosition,range.location-startPosition) replacementString:nil])
        {
            [theTextStore beginEditing];
            [theTextStore addAttributes:errorAttributes range:NSMakeRange(startPosition,range.location-startPosition)];
            [theTextStore endEditing];
            [self didChangeText];
        }
}


@end
