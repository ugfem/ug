/****************************************************************************/
/*																			*/
/* File:	  UGShellWindow.m												*/
/*																			*/
/* Purpose:   Implementation of the UGShellWindow Object.					*/
/*            See information below.										*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						    	*/
/*																			*/
/*	History:  January 16, 1998 begin										*/
/*																			*/
/****************************************************************************/


#import "UGShellWindow.h"

/* low module */
#include "misc.h"
#include "defaults.h"
#include "general.h"

/* interface includes */
#include "devices.h"
#include "cmdint.h"


unsigned int ugTextLength;
NXStream *ugTextStream;

@implementation UGShellWindow

- setUp
{
	NXRect	aRect,textRect;
	NXSize	contentSize;
	char	buffer[256];

	/*
	 * Set up a Panel
	 */
	NXSetRect(&aRect, 100.0, 700.0, 300.0, 40.0);
	infoPanel = [[Panel alloc] initContent:&aRect
					style:NX_TITLEDSTYLE
					backing:NX_BUFFERED
					buttonMask:NX_CLOSEBUTTONMASK
					defer:YES];
	[infoPanel setTitle:"About UG..."];
	[infoPanel removeFromEventMask:(NX_KEYDOWNMASK | NX_KEYUPMASK)];

	/*
	 * Set up a Menu
	 */
	theMenu = [[Menu alloc] initTitle:"UG Shell"];
	[[theMenu addItem:"Info..."
			action:@selector(orderFront:)
			keyEquivalent:'\0']
			setTarget:infoPanel];
	[theMenu addItem:"Hide"
			action:@selector(hide:)
			keyEquivalent:'h'];
	[theMenu addItem:"Quit"
			action:@selector(terminate:)
			keyEquivalent:'q'];
	[theMenu sizeToFit];
	//[theMenu display];
	[NXApp setMainMenu:theMenu];

	/* 
	 * Put a scrollView into theWindow
	 */
	//[self getFrame:&aRect];
	//aRect.origin.x = aRect.origin.y = 0.0;
	//printf("Windowframe size %d %d\n",aRect.size.width,aRect.size.height);
	NXSetRect(&aRect, 0.0, 0.0, 600.0, 300.0);
	theScrollView = [[ScrollView alloc] initFrame:&aRect];

	[theScrollView setBorderType:NX_BEZEL];
	[theScrollView setVertScrollerRequired:YES];
	[theScrollView setHorizScrollerRequired:YES];
	[theScrollView setDynamicScrolling:NX_WIDTHSIZABLE | NX_HEIGHTSIZABLE];
	[theScrollView setOpaque:YES];
	[theScrollView setBackgroundGray:NX_WHITE];

	/*
	 * Create and configure a NSTextView to fit into the NSScrollView
	 */
	[theScrollView getContentSize:&contentSize];
	NXSetRect(&aRect, 0.0, 0.0, contentSize.width, contentSize.height);
	theText = [[Text alloc] initFrame: &aRect  
					text:"" 
					alignment:NX_LEFTALIGNED];
	[theText setOpaque:YES];
	[theText setBackgroundGray:NX_WHITE];
	[theText setCharFilter:NXFieldFilter];
	[theText setDelegate:self];
	[[self contentView] addSubview:theText];
    [theText notifyAncestorWhenFrameChanged:YES];
	[theText setEditable:YES];
    shellFont = [Font newFont:"Ohlfs" size:10];
	[theText setFont:shellFont];

    [theText setMinSize:&(contentSize)];
	textRect.size.width  = 3000.0;
    textRect.size.height = 3000.0;
    [theText setMaxSize:&(textRect.size)];

    [theText setHorizResizable:YES];
    [theText setVertResizable:YES];
    [[theText superview] setAutoresizeSubviews:YES];
    [[theText superview] setAutosizing: NX_HEIGHTSIZABLE | 
										  NX_WIDTHSIZABLE];
	
	[theScrollView setDocView:theText];
	[self setContentView:theScrollView];
	[self makeKeyAndOrderFront:nil];
	[self orderFrontRegardless];
	[self makeFirstResponder:theText];

	[self makeKeyWindow];
	[theText selectAll:nil];
	
	ugTextLength = 0;
	ugTextStream = NXOpenMemory(0, 0, NX_READWRITE);

	printf("UGShellWindow: setUp done\n");
	
	return self;
}

- shellText
{
	return theText;
}

- appendToText:(const char *)val
{
	int length = [theText textLength];
	[self setDocEdited:YES];
	
	[theText setSel:length :length];

	[theText replaceSel:val];
	[theText scrollSelToVisible];
	[theText display];
	[theText setEditable: YES];
	return self;
}

- free
{
	NXCloseMemory(ugTextStream, NX_FREEBUFFER);
	return [super free];
}

@end

@implementation UGShellWindow(TextDelegate)

- textDidEnd:sender endChar:(unsigned short)whyEnd
{
	int		len,maxlen;
	char	*lastLine;
	char	*buf;
	NXStream	*str;
	INT error;
	int i,kerr;
	char c,errLine[256],spcLine[256],buffer[256];
	char *inpLine;
	char *strStart;
	int batch = FALSE;

	printf("In textDidEnd, whyEnd=%d, tab=%d\n",whyEnd,'\t');
	if (whyEnd==17)	{
		printf("Appending tab\n");
		[self appendToText:"\t"];
		return self;
		}
	str = NXOpenMemory(0, 0, NX_READWRITE);
	[sender writeText:str];
	NXGetMemoryBuffer(str, &buf, &len, &maxlen);
	lastLine = rindex(buf, '\n');
	if (lastLine)
		lastLine++;								/* skip past '\n' */
	else
		lastLine = buf;								/* get first line */

	if (lastLine=='\n')	{
		[self appendToText:"\n"];
		NXCloseMemory(str, NX_FREEBUFFER);
		return self;
		}
	
	printf("lastLine=%s\n", lastLine);
	if (strlen(lastLine) > 0)	{
		// Skip the prompt
		lastLine[0] = ' ';
		lastLine[1] = ' ';
		[self appendToText:"\n"];
		printf("Calling InterpretCommand(%s)\n", lastLine);
		if ((error=InterpretCommand(lastLine))!=DONE)	{
			if (error==QUITCODE)
				ExitUg();
			else
				UserWrite("Error position: \n");
			}
		}
	// Print the prompt
	[self appendToText:PROMPT];
	NXCloseMemory(str, NX_FREEBUFFER);

printf("\ttextDidEnd finished\n");

	return self;
}
@end
