/****************************************************************************/
/*																			*/
/* File:	  NeXTGraph.m													*/
/*																			*/
/* Purpose:   handling windows for graphical output 						*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						    	*/
/*																			*/
/*	History:  September 12, 1997 begin										*/
/*																			*/
/****************************************************************************/

#import <dpsclient/psops.h>



#import <appkit/appkit.h>


/* NeXT specific includes */
#import <appkit/Application.h>

/* standard C includes */
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* low includes */
#include "compiler.h"
#include "misc.h"
#include "heaps.h"
#include "defaults.h"
#include "general.h"

/* interface includes */
#include "devices.h"
#include "initdev.h"

/* mif includes */
#import "NeXTMain.h"
#import "NeXTGraph.h"
#import "UGshell.h"


extern id		theUGshell;	

static OUTPUTDEVICE *NeXTOutputDevice;

static void NeXTMove (SHORT_POINT point)
{
	PSmoveto(point.x, point.y);
}

static void NeXTDraw (SHORT_POINT point)
{
	PSlineto(point.x, point.y);
}

static void NeXTPolyline (SHORT_POINT *points, INT n)
{
	int i;
	
	if (n<2) return;

	PSmoveto(points[0].x, points[0].y);
	for (i=1; i<n; i++) 
		PSlineto(points[i].x, points[i].y);
}

static void NeXTInversePolyline (SHORT_POINT *points, INT n)
{
	return;
}

static void NeXTPolygon (SHORT_POINT *points, INT n)
{
	return;
}

static void NeXTInversePolygon (SHORT_POINT *points, INT n)
{
	return;
}

static void NeXTErasePolygon (SHORT_POINT *points, INT n)
{
	return;
}

static void Marker (short n, short s, SHORT_POINT point)
{
	return;
}

		
static void NeXTPolymark (short n, SHORT_POINT *points)
{
	return;
}

static void InvMarker (short n, short s, SHORT_POINT point)
{
	return;
}

		
static void NeXTInvPolymark (short n, SHORT_POINT *points)
{
	return;
}

static void NeXTText (const char *s, INT mode)
{
	return;
}

static void NeXTCenteredText (SHORT_POINT point, const char *s, INT mode)
{
	return;
}

static void NeXTClearViewPort (void)
{
	return;
}

static void NeXTSetLineWidth (short w)
{
	return;
}

static void NeXTSetTextSize (short s)
{
	return;
}

static void NeXTSetMarkerSize (short s)
{
	return;
}

static void NeXTSetMarker (short s)
{
	return;
}

static void NeXTSetColor (long index)
{
	return;
}

static void NeXTSetPaletteEntry (long index, short r, short g, short b)
{
	return;
}

static void NeXTSetNewPalette (long start, long count, short *r, short *g, short *b)
{
	return;
}

static void NeXTGetPaletteEntry (long index, short *r, short *g, short *b)
{
	return;
}

static void NeXTFlush (void)
{
	return;
}


static WINDOWID NeXT_OpenOutput (
	const char *title,						/* title of the window	 		*/
	INT x, INT y, INT width, INT height,	/* plot rgn in standard coord.	*/
	INT *Global_LL, INT *Global_UR, 		/* global machine coordinates	*/
	INT *Local_LL, INT *Local_UR,			/* local machine coordinates	*/
	INT *error) 							/* error code					*/
{
	return -1;
}


static INT NeXT_CloseOutput (WINDOWID win)
{
	return -1;
}


INT NeXT_ActivateOutput (WINDOWID win)
{
	return -1;
}

INT NeXT_UpdateOutput (WINDOWID win, INT tool)
{
	return -1;
}


static void InitNeXTPort()
{
	/*	init colors */
	NeXTOutputDevice->black				= 1;
	NeXTOutputDevice->white				= 0;
	NeXTOutputDevice->red				= 254;
	NeXTOutputDevice->green				= 170;
	NeXTOutputDevice->blue				= 3;
	NeXTOutputDevice->cyan				= 66;
	NeXTOutputDevice->orange 			= 223;
	NeXTOutputDevice->yellow 			= 192;
	NeXTOutputDevice->magenta			= 2;
	NeXTOutputDevice->hasPalette 		= 1;
	NeXTOutputDevice->range				= 256;
	NeXTOutputDevice->spectrumStart		= 3;
	NeXTOutputDevice->spectrumEnd		= 254;
	NeXTOutputDevice->signx		 		= 1;
	NeXTOutputDevice->signy		 		= -1;
	
	/* init pointers to basic drawing functions */
	NeXTOutputDevice->Move				= NeXTMove;
	NeXTOutputDevice->Draw				= NeXTDraw;
	NeXTOutputDevice->Polyline			= NeXTPolyline;
	NeXTOutputDevice->InversePolyline	= NeXTInversePolyline;
	NeXTOutputDevice->Polygon			= NeXTPolygon;
	NeXTOutputDevice->InversePolygon 	= NeXTInversePolygon;
	NeXTOutputDevice->ErasePolygon		= NeXTErasePolygon;
	NeXTOutputDevice->Polymark			= NeXTPolymark;
	NeXTOutputDevice->InvPolymark		= NeXTInvPolymark;
	NeXTOutputDevice->UGText				= NeXTText;
	NeXTOutputDevice->CenteredText		= NeXTCenteredText;
	NeXTOutputDevice->ClearViewPort		= NeXTClearViewPort; 					
	
	/* init pointers to set functions */
	NeXTOutputDevice->SetLineWidth		= NeXTSetLineWidth;
	NeXTOutputDevice->SetTextSize		= NeXTSetTextSize;
	NeXTOutputDevice->SetMarker			= NeXTSetMarker;
	NeXTOutputDevice->SetMarkerSize		= NeXTSetMarkerSize;
	NeXTOutputDevice->SetColor			= NeXTSetColor;
	NeXTOutputDevice->SetPaletteEntry	= NeXTSetPaletteEntry;
	NeXTOutputDevice->SetNewPalette		= NeXTSetNewPalette;
	
	/* init pointers to miscellaneous functions */
	NeXTOutputDevice->GetPaletteEntry	= NeXTGetPaletteEntry;
	NeXTOutputDevice->Flush				= NeXTFlush;
}



OUTPUTDEVICE *InitNeXTOutputDevice (void)
{
	/* create output device */
	if ((NeXTOutputDevice=CreateOutputDevice("screen"))==NULL) return(NULL);

	/* init output device 'x11' */
	NeXTOutputDevice->OpenOutput  = NeXT_OpenOutput;
	NeXTOutputDevice->CloseOutput  = NeXT_CloseOutput;
	NeXTOutputDevice->ActivateOutput  = NeXT_ActivateOutput;
	NeXTOutputDevice->UpdateOutput  = NeXT_UpdateOutput;

	NeXTOutputDevice->v.locked = 1;
	InitNeXTPort();

	printf("output device 'screen' for NeXT window manager created\n");

	return(NeXTOutputDevice);
}



void MousePosition (INT *ScreenPoint)
{
	NXPoint mousePoint;
	
	[[theUGshell shellWindow] getMouseLocation:&mousePoint];
	ScreenPoint[0] = (INT) mousePoint.x;
	ScreenPoint[1] = (INT) mousePoint.y;
	return;
}