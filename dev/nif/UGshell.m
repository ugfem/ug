/*
 * UGshell.m
 */

#import "UGshell.h"

@implementation UGshell

-setUp
{
	NXRect	 	aRect,textRect;
	NXSize 		contentSize;
	
	printf("In setUp\n");

	/*
	 * Set up a Panel
	 */
	NXSetRect(&aRect, 100.0, 700.0, 300.0, 40.0);
	infoPanel = [[Panel alloc] initContent:&aRect
					style:NX_TITLEDSTYLE
					backing:NX_BUFFERED
					buttonMask:NX_CLOSEBUTTONMASK
					defer:YES];
	[infoPanel setTitle:"About UG"];
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
	[NXApp setMainMenu:theMenu];

	/*
	 * Set up a Window
	 */
	NXSetRect(&aRect, 100.0, 350.0, 600.0, 300.0);
	shellWindow = [[Window alloc] initContent:&aRect
					style:NX_RESIZEBARSTYLE
					backing:NX_BUFFERED
					buttonMask:NX_MINIATURIZEBUTTONMASK
					defer:NO];
	[shellWindow setTitle:"UG shell"];
	[shellWindow display];
	
	/* 
	 * Put a scrollView into theWindow
	 */
	NXSetRect(&aRect, 0.0, 0.0, 600.0, 300.0);
	shellScrollView = [[ScrollView alloc] initFrame:&aRect];

	[shellScrollView setBorderType:NX_BEZEL];
	[shellScrollView setVertScrollerRequired:YES];
	[shellScrollView setHorizScrollerRequired:YES];
	[shellScrollView setDynamicScrolling:NX_WIDTHSIZABLE | NX_HEIGHTSIZABLE];
	[shellScrollView setOpaque:YES];
	[shellScrollView setBackgroundGray:NX_WHITE];

	/*
	 * Create and configure a NSTextView to fit into the NSScrollView
	 */
	[shellScrollView getContentSize:&contentSize];
	NXSetRect(&aRect, 0.0, 0.0, contentSize.width, contentSize.height);
	shellText = [[Text alloc] initFrame: &aRect  
					text:"" 
					alignment:NX_LEFTALIGNED];
	[shellText setOpaque:YES];
	[shellText setBackgroundGray:NX_WHITE];
	[shellText setCharFilter:NXFieldFilter];
	[shellText setDelegate:self];
	[[shellWindow contentView] addSubview:shellText];
    [shellText notifyAncestorWhenFrameChanged:YES];
	[shellText setEditable:YES];
    shellFont = [Font newFont:"Ohlfs" size:10];
	[shellText setFont:shellFont];

    [shellText setMinSize:&(contentSize)];
	textRect.size.width  = 3000.0;
    textRect.size.height = 3000.0;
    [shellText setMaxSize:&(textRect.size)];

    [shellText setHorizResizable:YES];
    [shellText setVertResizable:YES];
    [[shellText superview] setAutoresizeSubviews:YES];
    [[shellText superview] setAutosizing: NX_HEIGHTSIZABLE | 
										  NX_WIDTHSIZABLE];
	
	[shellScrollView setDocView:shellText];
	[shellWindow setContentView:shellScrollView];
	[shellWindow makeKeyAndOrderFront:nil];
	[shellWindow makeFirstResponder:shellText];

	[shellWindow orderFront:nil];
	[shellWindow makeKeyWindow];
	[shellText selectAll:nil];

	printf("setUp done\n");
	
	return self;
}

- shellWindow
{
	return shellWindow;
}

- shellText
{
	return shellText;
}

- appendToText:(const char *)val
{
	int length = [shellText textLength];
	[shellWindow setDocEdited:YES];
	
	[shellText setSel:length :length];

	[shellText replaceSel:val];

	[shellText scrollSelToVisible];
	[shellText display];

	return self;
}

@end


@implementation UGshell(TextDelegate)

/* Invoked when the text ends.
 * Find the last line and send it down the pipe...
 */

- textDidEnd:sender endChar:(unsigned short)whyEnd
{
	NXStream	*str;
	char		*lastLine;
	char		*buf;
	int			len, maxlen;

	printf("Entering textDidEnd\n");
	
	str = NXOpenMemory(0, 0, NX_READWRITE);
	[sender writeText:str];
	NXGetMemoryBuffer(str, &buf, &len, &maxlen);
	lastLine = strrchr(buf, '\n');
	
	printf("1. buf = %s, lastLine = %s\n",buf,lastLine);
	
	if (lastLine) {
		lastLine++;						/* skip past '\n' */
	}
	else{
		lastLine = buf;					/* get first line */
	}
	
	printf("2. buf = %s, lastLine = %s\n",buf,lastLine);
	
	if (strlen(lastLine) > 0){
		[self appendToText:"\n"];
		[self appendToText:": "];
		[self appendToText:lastLine];
		[self appendToText:"\n"];
		//[sender setEditable:NO]; /* wait for response */
	}
	NXCloseMemory(str, NX_FREEBUFFER);

	return self;
}

@end
