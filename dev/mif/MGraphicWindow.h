// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MGraphicWindow.h												*/
/*																			*/
/* Purpose:   Interface for MGraphicWindow class.							*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Interdisziplin"ares Zentrum f"ur Wissenschaftliches Rechnen	*/
/*			  Universit"at Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69210 Heidelberg												*/
/*			  email: Volker.Reichenberger@IWR.Uni-Heidelberg.DE		        */
/*																			*/
/*	History:  June 4, 1999 begin (based on OPENSTEP code)					*/
/*																			*/
/****************************************************************************/

#ifndef _MGRAPHICWINDOW_
#define _MGRAPHICWINDOW_

#import <AppKit.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DEVICESH__
#include "devices.h"
#endif

#include "MGraphicView.h"

@interface MGraphicWindow : NSWindow
{
  MGraphicView    *theView;                             // The graphic view within this window.
  NSTextView          *statusText;              // Information about window status goes here.
  /* drawing parameters */
  short marker_size;                                            // size of markers in pixels
  short marker_id;                                              // number of marker
  short textSize;                                               // text size

  INT currTool;                                                 // the current tool
  int moveto_x, moveto_y;                               // value is set in moveToPoint method

}

- (id)initWithContentRect:(NSRect) contentRect
 styleMask:(unsigned int)styleMask
 backing:(NSBackingStoreType) backingType
 defer:(BOOL)flag;
- (void) close;
- (void) activateOutput;
- (void) updateOutput;
- (void) moveToPoint:(SHORT_POINT)point;
- (void) drawLineTo:(SHORT_POINT)point;
- (void) drawPolyLine:(SHORT_POINT *)points noOfPoints:(INT)n;
- (void) drawInversePolyLine:(SHORT_POINT *)points noOfPoints:(INT)n;
- (void) drawPolygon:(SHORT_POINT *)points noOfPoints:(INT)n;
- (void) drawInversePolygon:(SHORT_POINT *)points noOfPoints:(INT)n;
- (void) erasePolygon:(SHORT_POINT *)points noOfPoints:(INT)n;
- (void) marker;
- (void) polyMark;
- (void) invMarker;
- (void) invPolyMark;
- (void) drawText:(NSString) text;
- (void) drawCenteredText:(NSString) text;
- (void) clearView;
- (void) setLineWidth:(short)width;
- (void) setTextSize:(short)width;
- (void) setMarkerSize:(short)width;
- (void) setMarker:(short)width;
- (void) setPaletteEntryWithIndex:(long)index red:(short)r green:(short)g blue:(short)b;
- (void) getPaletteEntryWithIndex:(long)index red:(short *)r green:(short *)g blue:(short *)b;
- (void) flush;
@end

#endif
