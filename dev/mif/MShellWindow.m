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
    self = [super initWithContentRect:contentRect styleMask:styleMask
                                 backing:backingType defer:flag];

    return self;
}

@end

