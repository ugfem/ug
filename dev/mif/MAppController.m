/****************************************************************************/
/*																			*/
/* File:	  MAppController.m												*/
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

#include <unistd.h>
#include <sys/param.h>

#import "MAppController.h"
#import "MShell.h"

extern MShell *theUGshell;

@implementation MAppController

- (void)appendText:(id)sender
{
    [theUGshell appendToText:[sender lossyCString]];
}

- (void)closeGraphic:(id)sender
{
	
}

- (void)executeCommand:(id)sender
{
    printf("execute command\n");
    //[theUGshell interpretCommand:sender];
}

- (void)runScript:(id)sender
{
    NSOpenPanel *openPanel = [NSOpenPanel openPanel];
    int runModalResult;
    char *buffer = malloc(MAXPATHLEN*sizeof(char));
    buffer = getwd(buffer);
    if (!_openPanelPath) _openPanelPath = [NSString stringWithCString:buffer];
    
    // Configure the OpenPanel
    [openPanel setTitle:NSLocalizedString(@"Run UG script", @"Select an UG script file to run.")];
    [openPanel setTreatsFilePackagesAsDirectories:NO];
    [openPanel setAllowsMultipleSelection:NO];
    [openPanel setCanChooseDirectories:NO];
    [openPanel setCanChooseFiles:YES];
	[openPanel setDirectory:_openPanelPath];
    //printf("%s\n", [_openPanelPath lossyCString]);
    free(buffer);

    // Run the panel
    runModalResult = [openPanel runModal];
    if (runModalResult == NSOKButton) {
        NSString *path = [openPanel filename];
        NSString *execute = @"execute ";
        NSString *command = [execute stringByAppendingString:path];
        printf ("%s\n", [command lossyCString]);
        [theUGshell interpretCommand:command];
    }
    [_openPanelPath release];
    _openPanelPath = [[openPanel directory] copyWithZone:[self zone]];

    return;
}

- (void)newGraphic:(id)sender
{
}

- (void)saveDocument:(id)sender
{
}

- (void)saveDocumentAs:(id)sender
{
}

- (void)showInfo:(id)sender
{
}

- (void)showInspector:(id)sender
{
}

- (void)showPreferences:(id)sender
{
}

- (void)dealloc {
    [_openPanelPath release];
    [super dealloc];
}

@end
