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

#import "MAppController.h"
#import "MShell.h"

@implementation MAppController

- (void)appendText:(id)sender
{
	[[MShell instantiate] appendToText:[sender lossyCString]];
}

- (void)closeGraphic:(id)sender
{
	
}

- (void)executeCommand:(id)sender
{
    [[MShell instantiate] interpretCommand:sender];
}

- (void)newGraphic:(id)sender
{
}

- (void)saveDocument:(id)sender
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

@end
