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

  NSBezierPath        *path;                            // Bezier path for drawing primitives.
  NSColor                     *currentColor;            // Current color for drawing primitives.
  NSColor                     *backgroundColor;         // Background color of graphic windows.
  NSColor                     *textColor;                       // Color of status text.
  float textsize;                                               // Current text size.

  NSPanel                 *toolPanel;                           // Select tools from this panel

  /* drawing parameters */
  short marker_size;                                            // size of markers in pixels
  short marker_id;                                              // number of marker
  short textSize;                                               // text size

  INT currTool;                                                 // The current tool.
  int moveto_x, moveto_y;                               // Value is set in moveToPoint method.
  int viewHeight;                                               // Height of the view.

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
- (void) drawShadedPolygon:(SHORT_POINT *)points noOfPoints:(INT) n intensity:(DOUBLE)it;
- (void) drawInversePolygon:(SHORT_POINT *)points noOfPoints:(INT)n;
- (void) erasePolygon:(SHORT_POINT *)points noOfPoints:(INT)n;

- (void) setMarkerSize:(short)width;
- (void) setMarker:(short)width;
- (void) polyMark:(SHORT_POINT *)points noOfPoints:(INT)n;
- (void) invMarker:(SHORT_POINT)point;
- (void) invPolyMark:(SHORT_POINT *)points noOfPoints:(INT)n;

- (void) drawText:(const char*)text mode:(INT)m;
- (void) drawCenteredText:(const char*)text atPoint:(SHORT_POINT) point mode:(INT)m;

- (void) clearView;

- (void) setColorRed:(float)r green:(float)g blue:(float)b;
- (void) setLineWidth:(short)width;
- (void) setTextSize:(short)h;

- (void) setPaletteEntryWithIndex:(long)index red:(short)r green:(short)g blue:(short)b;
- (void) setNewPaletteEntryWithIndex:(long)index withCount:(long)count red:(short)r green:(short)g blue:(short)b;
- (void) getPaletteEntryWithIndex:(long)index red:(short *)r green:(short *)g blue:(short *)b;
- (void) flush;

@end

#endif
