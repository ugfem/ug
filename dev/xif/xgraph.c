// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      xgraph.c                                                      */
/*                                                                          */
/* Purpose:   graphical output for ug based on X11                          */
/*                                                                          */
/* Author:    Peter Bastian / Klaus Johannsen                               */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/* email:     ug@ica3.uni-stuttgart.de                                      */
/*                                                                          */
/* History:   17.02.94 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*			  this version uses read only color cells, therefore no                 */
/*			  it does not allow to change the color table during run time	*/
/*																			*/
/*			  Changing text size is not supported yet						*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* X11 includes */
#include <config.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/Xresource.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#include <X11/X.h>

/* standard C includes */
#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <string.h>

/* interface includes */
#include "ugtypes.h"
#include "ugdevices.h"
#include "initdev.h"
#include "general.h"
#include "debug.h"
#include "commands.hh"
#include "pixel.h"

/* Xif includes */
#include "xmain.h"
#include "xgraph.h"

USING_UG_NAMESPACES

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define CSIZE                   256                     /* ug supports at most 256 colors	*/

#define DEFAULTXPOS     50
#define DEFAULTYPOS     100
#define DEFAULTFONT     "9x15"
#define DEFAULTMINX     100
#define DEFAULTMINY     100
#define DEFAULTBORDER   4
#define GRAPHICONNAME   "ug3 view"
#define GRAPHWINNAME    "ug3 view"
#define RESOURCENAME    "ug3"

#define NTOOLS                  7                               /* this is the number of tools		*/

/* bitmaps for icons */
#include "view-icon"
#include "tool0"
#include "tool1"
#include "tool2"
#include "tool3"
#include "tool4"
#include "tool5"
#include "tool6"
#include "shades.h"

/* hack for C++ */
#ifdef __cplusplus
#define CLASS c_class
#else
#define CLASS class
#endif

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* current context */
static GraphWindow *gw;                                 /* active output					*/
static XColor ctab[CSIZE];                              /* one static color table for all wi*/
static int ncolors=0;                                   /* number of entries used			*/
static unsigned long pixels[CSIZE];     /* array to free color cells		*/
static OUTPUTDEVICE *X11OutputDevice;   /* outputdevice that has been initi */
static GraphWindow *windowList=NULL;    /* list of our windows				*/
static Visual *default_visual;
static int default_depth;

static Pixmap tools[NTOOLS];                    /* bitmaps for toolbox				*/
static Cursor cursors[NTOOLS];                  /* id in cursor font				*/

static int MoveMouse = 1;                               /* some local vars for mouse ops	*/
static int StoredMousePos = 0;
static int MouseX,MouseY;

static Colormap our_cmap;             /* may be a private one               */
static int true_color;                /* whether we have true color visual  */
static int current_color;             /* index to current foreground color  */

/* pixmaps for shading patterns */
static Pixmap pattern[NO_PATTERNS];

/* pixel structure for true color visual */
static unsigned long red_mask, green_mask, blue_mask;
static int red_shift, green_shift, blue_shift;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

unsigned long NS_PREFIX UGBlack (void)
{
  return(ctab[X11OutputDevice->black].pixel);
}

unsigned long NS_PREFIX UGWhite (void)
{
  return(ctab[X11OutputDevice->white].pixel);
}

/*==========================================================================*/
/*																			*/
/* Initialization and port functions										*/
/*																			*/
/*==========================================================================*/

static void IFMove (SHORT_POINT point)
{
  gw->x = (int) point.x;
  gw->y = (int) point.y;
}

static void IFDraw (SHORT_POINT point)
{
  int x1, y1, x2, y2;

  /* always draw from left to right to get the same pixels */
  if (point.x <= gw->x) {
    x1 = point.x;  y1 = point.y;
    x2 = gw->x;    y2 = gw->y;
  }
  else {
    x1 = gw->x;    y1 = gw->y;
    x2 = point.x;  y2 = point.y;
  }
  XDrawLine(display,gw->win,gw->gc,x1,y1,x2,y2);
  if (!gw->backing_store)
    XDrawLine(display,gw->pixmap,gw->gc,x1,y1,x2,y2);
  gw->x = (int) point.x;
  gw->y = (int) point.y;
}

static void IFPolyline (SHORT_POINT *points, INT n)
{
  XDrawLines(display,gw->win,gw->gc,(XPoint *)points,n,CoordModeOrigin);
  if (!gw->backing_store)
    XDrawLines(display,gw->pixmap,gw->gc,(XPoint *)points,n,CoordModeOrigin);
}

static void IFInversePolyline (SHORT_POINT *points, INT n)
{
  XGCValues values_return;
  unsigned long valuemask,plane_mask;
  int function;

  valuemask = GCPlaneMask;
  XGetGCValues(display,gw->gc,valuemask,&values_return);
  plane_mask = values_return.plane_mask;
  valuemask = GCFunction;
  XGetGCValues(display,gw->gc,valuemask,&values_return);
  function = values_return.function;
  XSetFunction(display,gw->gc,GXinvert);
  XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
  XDrawLines(display,gw->win,gw->gc,(XPoint *)points,n,CoordModeOrigin);
  if (!gw->backing_store)
    XDrawLines(display,gw->pixmap,gw->gc,(XPoint *)points,n,CoordModeOrigin);
  XSetFunction(display,gw->gc,function);
  XSetPlaneMask(display,gw->gc,plane_mask);
}

static void IFPolygon (SHORT_POINT *points, INT n)
{
  XFillPolygon(display,gw->win,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
  if (!gw->backing_store)
    XFillPolygon(display,gw->pixmap,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
}

static void IFShadedPolygon(SHORT_POINT *points, INT n, DOUBLE intensity)
{
  if (!true_color) {
    int s;

    XSetBackground(display,gw->gc,ctab[X11OutputDevice->black].pixel);
    s = (int)(0.5 + intensity * (DOUBLE)(NO_PATTERNS-1));
    XSetFillStyle(display, gw->gc, FillOpaqueStippled);
    XSetStipple(display, gw->gc, pattern[s]);
    XFillPolygon(display,gw->win,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
    if (!gw->backing_store)
      XFillPolygon(display,gw->pixmap,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
    XSetFillStyle(display, gw->gc, FillSolid);
    XSetBackground(display,gw->gc,ctab[X11OutputDevice->white].pixel);
  }
  else {
    unsigned long pixel;
    DOUBLE red, green, blue;

    pixel = ctab[current_color].pixel;
    red   = (pixel & red_mask  ) >> red_shift;
    green = (pixel & green_mask) >> green_shift;
    blue  = (pixel & blue_mask ) >> blue_shift;
    red   *= intensity;
    green *= intensity;
    blue  *= intensity;
    pixel = ((unsigned long)(red  +0.5) << red_shift  ) +
            ((unsigned long)(green+0.5) << green_shift) +
            ((unsigned long)(blue +0.5) << blue_shift );
    XSetForeground(display, gw->gc, pixel);
    XFillPolygon(display, gw->win, gw->gc, (XPoint *)points, n, Convex, CoordModeOrigin);
    if (!gw->backing_store)
      XFillPolygon(display,gw->pixmap,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
    XSetForeground(display, gw->gc, ctab[current_color].pixel);
  }
}

static void IFInversePolygon (SHORT_POINT *points, INT n)
{
  XGCValues values_return;
  unsigned long valuemask,plane_mask;
  int function;

  valuemask = GCPlaneMask;
  XGetGCValues(display,gw->gc,valuemask,&values_return);
  plane_mask = values_return.plane_mask;
  valuemask = GCFunction;
  XGetGCValues(display,gw->gc,valuemask,&values_return);
  function = values_return.function;
  XSetFunction(display,gw->gc,GXinvert);
  XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
  XFillPolygon(display,gw->win,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
  if (!gw->backing_store)
    XFillPolygon(display,gw->pixmap,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
  XSetFunction(display,gw->gc,function);
  XSetPlaneMask(display,gw->gc,plane_mask);
}

static void IFErasePolygon (SHORT_POINT *points, INT n)
{
  XGCValues gcv;

  XGetGCValues(display, gw->gc, GCForeground|GCBackground, &gcv);
  XSetForeground(display,gw->gc,gcv.background);
  XFillPolygon(display,gw->win,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
  XDrawLines(display,gw->win,gw->gc,(XPoint *)points,n,CoordModeOrigin);
  XDrawLine(display,gw->win,gw->gc,(int)points[n-1].x,(int)points[n-1].y,(int)points[0].x,(int)points[0].y);
  if (!gw->backing_store) {
    XFillPolygon(display,gw->pixmap,gw->gc,(XPoint *)points,n,Convex,CoordModeOrigin);
    XDrawLines(display,gw->pixmap,gw->gc,(XPoint *)points,n,CoordModeOrigin);
    XDrawLine(display,gw->pixmap,gw->gc,(int)points[n-1].x,(int)points[n-1].y,(int)points[0].x,(int)points[0].y);
  }
  XSetForeground(display,gw->gc,gcv.foreground);
}

static void Marker (short n, short s, SHORT_POINT point)
{
  short x,y;
  short top, left, bottom, right;

  x = point.x;
  y = point.y;
  top = y-s/2; bottom = y+s/2;
  left = x-s/2; right = x+s/2;
  n = n%11;

  switch (n)
  {
  case 0 :
    XDrawRectangle(display, gw->win, gw->gc, left, top, s, s);
    if (!gw->backing_store)
      XDrawRectangle(display, gw->pixmap, gw->gc, left, top, s, s);
    break;
  case 1 :
    XDrawRectangle(display, gw->win, gw->gc, left, top, s, s);
    if (!gw->backing_store)
      XDrawRectangle(display, gw->pixmap, gw->gc, left, top, s, s);
    break;
  case 2 :
    XFillRectangle(display, gw->win, gw->gc, left, top, s, s);
    if (!gw->backing_store)
      XFillRectangle(display, gw->pixmap, gw->gc, left, top, s, s);
    break;
  case 3 :
    XDrawArc( display, gw->win, gw->gc, left, top, s, s, 0, 360*64);
    if (!gw->backing_store)
      XDrawArc( display, gw->pixmap, gw->gc, left, top, s, s, 0, 360*64);
    break;
  case 4 :
    XDrawArc( display, gw->win, gw->gc, left, top, s, s, 0, 360*64);
    if (!gw->backing_store)
      XDrawArc( display, gw->pixmap, gw->gc, left, top, s, s, 0, 360*64);
    break;
  case 5 :
    XFillArc( display, gw->win, gw->gc, left, top, s+1, s+1, 0, 360*64);
    if (!gw->backing_store)
      XFillArc( display, gw->pixmap, gw->gc, left, top, s+1, s+1, 0, 360*64);
    break;
  case 6 :
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, gw->win, gw->gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x+s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x+s/2, y, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x, y-s/2, x-s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x, y+s/2);
    }
    break;
  case 7 :
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, gw->win, gw->gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x+s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x+s/2, y, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x, y-s/2, x-s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x, y+s/2);
    }
    break;
  case 8 :
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, gw->win, gw->gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x+s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x+s/2, y, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x, y-s/2, x-s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x, y+s/2);
    }
    break;
  case 9 :
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x+s/2, y);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x+s/2, y);
    }
    break;
  case 10 :
    XDrawLine( display, gw->win, gw->gc, x-s/2, y+s/2, x+s/2, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y-s/2, x+s/2, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y+s/2, x+s/2, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y-s/2, x+s/2, y+s/2);
    }
    break;
  }
}


static void IFPolymark (short n, SHORT_POINT *points)
{
  int i;

  for (i=0; i<n; i++) Marker(gw->marker_id,gw->marker_size,points[i]);
}

static void InvMarker (short n, short s, SHORT_POINT point)
{
  XGCValues values_return;
  unsigned long valuemask,plane_mask;
  int function;
  short x,y;
  short top, left, bottom, right;

  x = point.x;
  y = point.y;
  top = y-s/2; bottom = y+s/2;
  left = x-s/2; right = x+s/2;
  n = n%11;

  /* save values */
  valuemask = GCPlaneMask;
  XGetGCValues(display,gw->gc,valuemask,&values_return);
  plane_mask = values_return.plane_mask;
  valuemask = GCFunction;
  XGetGCValues(display,gw->gc,valuemask,&values_return);
  function = values_return.function;
  XSetFunction(display,gw->gc,GXinvert);

  switch (n)
  {
  case EMPTY_SQUARE_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawRectangle(display, gw->win, gw->gc, left, top, s, s);
    if (!gw->backing_store)
      XDrawRectangle(display, gw->pixmap, gw->gc, left, top, s, s);
    break;
  case GRAY_SQUARE_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawRectangle(display, gw->win, gw->gc, left, top, s, s);
    if (!gw->backing_store)
      XDrawRectangle(display, gw->pixmap, gw->gc, left, top, s, s);
    break;
  case FILLED_SQUARE_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XFillRectangle(display, gw->win, gw->gc, left, top, s, s);
    if (!gw->backing_store)
      XFillRectangle(display, gw->pixmap, gw->gc, left, top, s, s);
    break;
  case EMPTY_CIRCLE_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawArc( display, gw->win, gw->gc, left, top, s, s, 0, 360*64);
    if (!gw->backing_store)
      XDrawArc( display, gw->pixmap, gw->gc, left, top, s, s, 0, 360*64);
    break;
  case GRAY_CIRCLE_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawArc( display, gw->win, gw->gc, left, top, s, s, 0, 360*64);
    if (!gw->backing_store)
      XDrawArc( display, gw->pixmap, gw->gc, left, top, s, s, 0, 360*64);
    break;
  case FILLED_CIRCLE_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XFillArc( display, gw->win, gw->gc, left, top, s+1, s+1, 0, 360*64);
    if (!gw->backing_store)
      XFillArc( display, gw->pixmap, gw->gc, left, top, s+1, s+1, 0, 360*64);
    break;
  case EMPTY_RHOMBUS_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, gw->win, gw->gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x+s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x+s/2, y, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x, y-s/2, x-s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x, y+s/2);
    }
    break;
  case GRAY_RHOMBUS_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, gw->win, gw->gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x+s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x+s/2, y, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x, y-s/2, x-s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x, y+s/2);
    }
    break;
  case FILLED_RHOMBUS_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, gw->win, gw->gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x+s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x+s/2, y, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x, y-s/2, x-s/2, y);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x, y+s/2);
    }
    break;
  case PLUS_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawLine( display, gw->win, gw->gc, x, y+s/2, x, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y, x+s/2, y);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x, y+s/2, x, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y, x+s/2, y);
    }
    break;
  case CROSS_MARKER :
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y+s/2, x+s/2, y-s/2);
    XDrawLine( display, gw->win, gw->gc, x-s/2, y-s/2, x+s/2, y+s/2);
    if (!gw->backing_store) {
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y+s/2, x+s/2, y-s/2);
      XDrawLine( display, gw->pixmap, gw->gc, x-s/2, y-s/2, x+s/2, y+s/2);
    }
    break;
  }
  /* restore values */
  XSetFunction(display,gw->gc,function);
  XSetPlaneMask(display,gw->gc,plane_mask);
}


static void IFInvPolymark (short n, SHORT_POINT *points)
{
  int i;

  for (i=0; i<n; i++) Marker(gw->marker_id,gw->marker_size,points[i]);
}

static void IFDrawText (const char *s, INT mode)
{
  XGCValues values_return;
  unsigned long valuemask,plane_mask;
  int function;

  if (mode==TEXT_REGULAR)
  {
    XDrawString(display,gw->win,gw->gc,gw->x,gw->y,s,strlen(s));
    if (!gw->backing_store)
      XDrawString(display,gw->pixmap,gw->gc,gw->x,gw->y,s,strlen(s));
  }
  else
  {
    valuemask = GCPlaneMask;
    XGetGCValues(display,gw->gc,valuemask,&values_return);
    plane_mask = values_return.plane_mask;
    valuemask = GCFunction;
    XGetGCValues(display,gw->gc,valuemask,&values_return);
    function = values_return.function;
    XSetFunction(display,gw->gc,GXinvert);
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawString(display,gw->win,gw->gc,gw->x,gw->y,s,strlen(s));
    if (!gw->backing_store)
      XDrawString(display,gw->pixmap,gw->gc,gw->x,gw->y,s,strlen(s));
    XSetFunction(display,gw->gc,function);
    XSetPlaneMask(display,gw->gc,plane_mask);
  }
}

static void IFCenteredText (SHORT_POINT point, const char *s, INT mode)
{
  int ts,w;
  XGCValues values_return;
  unsigned long valuemask,plane_mask;
  int function;

  w = XTextWidth(gw->font_info,s,strlen(s));
  ts = gw->font_height;
  if (mode==TEXT_REGULAR)
  {
    XDrawString(display,gw->win,gw->gc,((int)point.x)-w/2,((int)point.y)+ts/2,s,strlen(s));
    if (!gw->backing_store)
      XDrawString(display,gw->pixmap,gw->gc,((int)point.x)-w/2,((int)point.y)+ts/2,s,strlen(s));
  }
  else
  {
    valuemask = GCPlaneMask;
    XGetGCValues(display,gw->gc,valuemask,&values_return);
    plane_mask = values_return.plane_mask;
    valuemask = GCFunction;
    XGetGCValues(display,gw->gc,valuemask,&values_return);
    function = values_return.function;
    XSetFunction(display,gw->gc,GXinvert);
    XSetPlaneMask(display,gw->gc,0xFFFFFFFF);
    XDrawString(display,gw->win,gw->gc,((int)point.x)-w/2,((int)point.y)+ts/2,s,strlen(s));
    if (!gw->backing_store)
      XDrawString(display,gw->pixmap,gw->gc,((int)point.x)-w/2,((int)point.y)+ts/2,s,strlen(s));
    XSetFunction(display,gw->gc,function);
    XSetPlaneMask(display,gw->gc,plane_mask);
  }
}

static void IFClearViewPort (void)
{
  XClearArea(display,gw->win,0,0,gw->window_width,gw->window_height-CONTROLSIZE-1,0);
}

static void IFSetLineWidth (short w)
{
  XGCValues gcv;
  if (w==1)
    gcv.line_width = 0;
  else
    gcv.line_width = (int)w;
  XChangeGC( display, gw->gc, GCLineWidth, &gcv);
}

static void IFSetTextSize (short s)
{
  return;
}

static void IFSetMarkerSize (short s)
{
  gw->marker_size = s;
  return;
}

static void IFSetMarker (short s)
{
  gw->marker_id = s;
  return;
}

static void IFSetColor (long index)
{
  if (index<0 || index>=CSIZE)
    return;
  XSetForeground(display,gw->gc,ctab[(int)index].pixel);
  current_color = index;
}

static void IFSetPaletteEntry (long index, short r, short g, short b)
{
  return;
}

static void IFSetNewPalette (long start, long count, short *r, short *g, short *b)
{
  return;
}

static void IFGetPaletteEntry (long index, short *r, short *g, short *b)
{
  *r = (short) ctab[(int)index].red>>8;
  *g = (short) ctab[(int)index].green>>8;
  *b = (short) ctab[(int)index].blue>>8;
}

static void IFFlush (void)
{
  if (gw->gc==NULL) return;
  XFlushGC(display,gw->gc);
  return;
}

/*--------------------------------------------------------------------------*
*   Convert pixel buffer to XImage data: one has to copy the 16 or 32      *
*   lower bits of the pixel values keeping the clients byte order.         *
*--------------------------------------------------------------------------*/

static void MakeXImageData(void *buffer, void *data, int len, int pad)
{
  int i, c;
  unsigned long pixel, r, g, b, f;
  char *d;
  PIXEL *p;

  p = (PIXEL *)buffer;
  d = (char  *)data;

  for (i = 0; i < len; i++) {
    c = p->cindex;
    f = p->intensity;
    p++;
    r = (ctab[c].pixel & red_mask  ) >> red_shift;
    g = (ctab[c].pixel & green_mask) >> green_shift;
    b = (ctab[c].pixel & blue_mask ) >> blue_shift;
    r = r*f/255;
    g = g*f/255;
    b = b*f/255;
    pixel = (r << red_shift) + (g << green_shift) + (b << blue_shift);
    switch (pad)
    {
    case 16 :
#ifdef __SWAPBYTES__
      d[0] = ((char*)&pixel)[0];
      d[1] = ((char*)&pixel)[1];
#else
      d[0] = ((char*)&pixel)[sizeof(pixel)-2];
      d[1] = ((char*)&pixel)[sizeof(pixel)-1];
#endif
      d += 2;
      break;
    case 32 :
#ifdef __SWAPBYTES__
      d[0] = ((char*)&pixel)[0];
      d[1] = ((char*)&pixel)[1];
      d[2] = ((char*)&pixel)[2];
      d[3] = ((char*)&pixel)[3];
#else
      d[0] = ((char*)&pixel)[sizeof(pixel)-4];
      d[1] = ((char*)&pixel)[sizeof(pixel)-3];
      d[2] = ((char*)&pixel)[sizeof(pixel)-2];
      d[3] = ((char*)&pixel)[sizeof(pixel)-1];
#endif
      d += 4;
    }

  }
}

/*--------------------------------------------------------------------------*
*   write the pixel buffer to the screen                                   *
*--------------------------------------------------------------------------*/

static void IFPlotPixelBuffer(void *buffer, void *data,
                              int x, int y, int w, int h)
{
  HEAP *heap;
  XImage *image;
  int bitmap_pad;
  INT mem_allocated, key;

  /* We support only TrueColor, this is 2004! */
  if (default_visual->CLASS != TrueColor) {
    UserWrite("IFPlotPixelBuffer: Sorry, only TrueColor visuals supported.\n");
    return;
  }

  /* how many bits per pixel do we need? */
  if (default_depth <= 16)
    bitmap_pad = 16;
  else
    bitmap_pad = 32;

  /* allocate mem for XImage data if necessary */
  mem_allocated = NO;
  if (data == NULL) {
    heap = GetCurrentMultigrid()->theHeap;
    MarkTmpMem(heap, &key);
    if ((data = GetTmpMem(heap, w*h*bitmap_pad/8, key)) == NULL) {
      UserWrite("IFPlotPixelBuffer: Sorry, not enough memory.\n");
      ReleaseTmpMem(heap, key);
      return;
    }
    mem_allocated = YES;
  }

  MakeXImageData(buffer, data, w*h, bitmap_pad);

  /* create XImage structure and write image to the screen */
  image = XCreateImage(display, default_visual, default_depth, ZPixmap, 0,
                       (char*)data, w, h, bitmap_pad, 0);
#ifdef __SWAPBYTES__
  image->byte_order = LSBFirst;
#else
  image->byte_order = MSBFirst;
#endif
  XPutImage(display, gw->win, gw->gc, image, 0, 0, x, y, w, h);
  if (!gw->backing_store)
    XPutImage(display, gw->pixmap, gw->gc, image, 0, 0, x, y, w, h);
  image->data = NULL;
  XDestroyImage(image);

  if (mem_allocated) ReleaseTmpMem(heap, key);
}

/*--------------------------------------------------------------------------*/

static int get_component_shift(unsigned long mask)
{
  int shift;

  shift=0;
  while (!(mask & 1)) {
    shift++;
    mask >>= 1;
  }
  return shift;
}

/****************************************************************************/
/*
   InitXPort - implement basic drawing functions by X11

   SYNOPSIS:
   void InitXPort (OUTPUTDEVICE *thePort);

   PARAMETERS:
   .  thePort - PORT structure to initialize

   DESCRIPTION:
   This function  implements basic drawing functions by X11.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void InitXPort (OUTPUTDEVICE *thePort)
{
  Colormap default_cmap;
  int i,j,Private;
  unsigned short res,delta,max,r,g,b;

  if (ncolors!=0) return;

  thePort->signx          = 1;
  thePort->signy          = -1;

  /* init color table */
  ncolors = 0;

  /* get default values */
  default_depth =  DefaultDepth(display,screen_num);
  default_visual = DefaultVisual(display,screen_num);
  default_cmap =   DefaultColormap(display,screen_num);

  /* print display type */
  switch (default_visual->CLASS)
  {
  case PseudoColor :
    printf("visual=%s depth=%d\n","PseudoColor",default_depth);
    break;
  case GrayScale :
    printf("visual=%s depth=%d\n","GrayScale",default_depth);
    break;
  case DirectColor :
    printf("visual=%s depth=%d\n","DirectColor",default_depth);
    break;
  case TrueColor :
    printf("visual=%s depth=%d\n","TrueColor",default_depth);
    break;
  case StaticGray :
    printf("visual=%s depth=%d\n","StaticGray",default_depth);
    break;
  case StaticColor :
    printf("visual=%s depth=%d\n","StaticColor",default_depth);
    break;
  default :
    printf("visual=%s depth=%d\n","unknown",default_depth);
    break;
  }

  /* get info on pixel structure if true color */
  if (true_color = (default_visual->CLASS == TrueColor)) {
    red_mask   = default_visual->red_mask;
    green_mask = default_visual->green_mask;
    blue_mask  = default_visual->blue_mask;

    red_shift  = get_component_shift(red_mask);
    green_shift= get_component_shift(green_mask);
    blue_shift = get_component_shift(blue_mask);
  }

  our_cmap = default_cmap;

  /* check for B&W */
  if (default_depth==1)
  {
    ncolors = 2;
    ctab[0].pixel = WhitePixel(display,screen_num);
    ctab[0].red = 0xFFFF;
    ctab[0].green = 0xFFFF;
    ctab[0].blue = 0xFFFF;
    ctab[1].pixel = BlackPixel(display,screen_num);
    ctab[1].red = 0;
    ctab[1].green = 0;
    ctab[1].blue = 0;
    thePort->black                          = 1;
    thePort->gray               = 1;
    thePort->white                          = 0;
    thePort->red                            = 1;
    thePort->green                          = 1;
    thePort->blue                           = 1;
    thePort->cyan                           = 1;
    thePort->orange                         = 1;
    thePort->yellow                         = 1;
    thePort->darkyellow             = 1;
    thePort->magenta                        = 1;
    thePort->hasPalette             = 0;
    thePort->range                          = 2;
    thePort->spectrumStart          = 0;
    thePort->spectrumEnd            = 1;
    printf("Using B&W color map\n");
  }

  if ((default_depth>=6)&&(ncolors==0)&&((default_visual->CLASS==GrayScale)||(default_visual->CLASS==StaticGray)))
  {
    /* lets try 64 gray scale map */
    ncolors = 64;
    for (i=0; i<64; i++)
    {
      ctab[63-i].red = 1024*i;
      ctab[63-i].green = 1024*i;
      ctab[63-i].blue = 1024*i;
      if (XAllocColor(display,our_cmap,&(ctab[i]))==0)
      {
        ncolors = 0;
        for (j=0; j<i; j++) pixels[j] = ctab[j].pixel;
        XFreeColors(display,our_cmap,pixels,i,0);
        break;
      }
    }
    thePort->black = 63;
    thePort->gray = 48;
    thePort->white = 0;
    thePort->red = 36;
    thePort->green = 27;
    thePort->blue = 54;
    thePort->cyan = 18;
    thePort->yellow = 9;
    thePort->darkyellow  = 10;
    thePort->magenta = 42;
    thePort->orange = 29;
    thePort->hasPalette = 1;
    thePort->range = 64;
    thePort->spectrumStart = 6;
    thePort->spectrumEnd = 63;
    if (ncolors>0) printf("Using 64 grayscale color map\n");
  }

  if ((default_depth>=4)&&(ncolors==0)&&((default_visual->CLASS==GrayScale)||(default_visual->CLASS==StaticGray)))
  {
    /* lets try a 16 gray scale map */
    ncolors = 16;
    for (i=0; i<16; i++)
    {
      ctab[15-i].red = 4096*i;
      ctab[15-i].green = 4096*i;
      ctab[15-i].blue = 4096*i;
      if (XAllocColor(display,our_cmap,&(ctab[i]))==0)
      {
        ncolors = 0;
        for (j=0; j<i; j++) pixels[j] = ctab[j].pixel;
        XFreeColors(display,our_cmap,pixels,i,0);
        break;
      }
    }
    thePort->black = 15;
    thePort->gray = 13;
    thePort->white = 0;
    thePort->red = 11;
    thePort->green = 9;
    thePort->blue = 13;
    thePort->cyan = 5;
    thePort->yellow = 5;
    thePort->darkyellow = 5;
    thePort->magenta = 7;
    thePort->orange = 10;
    thePort->hasPalette = 1;
    thePort->range = 16;
    thePort->spectrumStart = 1;
    thePort->spectrumEnd = 15;
    if (ncolors>0) printf("Using 16 grayscale color map\n");
  }

  if ((default_depth>=8)&&(ncolors==0))
  {
    /* should be a color device */
    res = 31;
    delta = 2048;
    max = 63488;

    /* fixed colors */
    i = 0;
    ctab[i].red = 0xFFFF; ctab[i].green = 0xFFFF; ctab[i].blue = 0xFFFF; i++;
    ctab[i].red = 0xD000; ctab[i].green = 0xD000; ctab[i].blue = 0xD000; i++;
    ctab[i].red = 0xFFFF; ctab[i].green = 0xFFFF; ctab[i].blue = 0x0   ; i++;
    ctab[i].red = 0xFFFF; ctab[i].green = 0x0       ; ctab[i].blue = 0xFFFF; i++;
    ctab[i].red = 0xFFFF; ctab[i].green = 0x0       ; ctab[i].blue = 0x0   ; i++;
    ctab[i].red = 0x0       ; ctab[i].green = 0xFFFF; ctab[i].blue = 0xFFFF; i++;
    ctab[i].red = 0x0       ; ctab[i].green = 0xFFFF; ctab[i].blue = 0x0   ; i++;
    ctab[i].red = 0x0       ; ctab[i].green = 0x0   ; ctab[i].blue = 0xFFFF; i++;
    ctab[i].red = 0x0       ; ctab[i].green = 0x0   ; ctab[i].blue = 0x0   ; i++;
    ctab[i].red = 65520 ; ctab[i].green = 32240 ; ctab[i].blue = 0x0   ; i++;
    ctab[i].red = 65520 ; ctab[i].green = 60000 ; ctab[i].blue = 0x0   ; i++;

    /* color spectrum */
    r = g = 0; b = max;
    ctab[i].red = r; ctab[i].green = g; ctab[i].blue = b; i++;

    /* blau nach cyan */
    for (j=0; j<res; j++)
    {
      g += delta;
      ctab[i].red = r; ctab[i].green = g; ctab[i].blue = b; i++;
    }
    /* cyan nach gruen */
    for (j=0; j<res; j++)
    {
      b -= delta;
      ctab[i].red = r; ctab[i].green = g; ctab[i].blue = b; i++;
    }
    /* gruen nach gelb */
    for (j=0; j<res; j++)
    {
      r += delta;
      ctab[i].red = r; ctab[i].green = g; ctab[i].blue = b; i++;
    }
    /* gelb nach rot */
    for (j=0; j<res; j++)
    {
      g -= delta;
      ctab[i].red = r; ctab[i].green = g; ctab[i].blue = b; i++;
    }

    thePort->black = 8;
    thePort->gray = 1;
    thePort->white = 0;
    thePort->red = 4;
    thePort->green = 6;
    thePort->blue = 7;
    thePort->cyan = 5;
    thePort->yellow = 2;
    thePort->darkyellow = 10;
    thePort->magenta = 3;
    thePort->orange = 9;
    thePort->hasPalette = 1;
    thePort->range = i;
    thePort->spectrumStart = 11;
    thePort->spectrumEnd = i-1;

    ncolors = i;
    Private = NO;

    for (i=0; i<ncolors; i++)
      if (!XAllocColor(display, our_cmap, &ctab[i])) {
        Private = YES;
        our_cmap = XCopyColormapAndFree(display, our_cmap);
        i--;                          /* allocate the color that failed again! */
      }

    if (Private)
      printf("Using private color map with %d entries\n",ncolors);
    else
      printf("Using default color map with %d entries\n",ncolors);
  }

  if (ncolors==0)
  {
    /* settle with b&w */
    ncolors = 2;
    ctab[0].pixel = WhitePixel(display,screen_num);
    ctab[0].red = 0xFFFF;
    ctab[0].green = 0xFFFF;
    ctab[0].blue = 0xFFFF;
    ctab[1].pixel = BlackPixel(display,screen_num);
    ctab[1].red = 0;
    ctab[1].green = 0;
    ctab[1].blue = 0;
    thePort->black                          = 1;
    thePort->gray               = 1;
    thePort->white                          = 0;
    thePort->red                            = 1;
    thePort->green                          = 1;
    thePort->blue                           = 1;
    thePort->cyan                           = 1;
    thePort->orange                         = 1;
    thePort->yellow                         = 1;
    thePort->darkyellow             = 1;
    thePort->magenta                        = 1;
    thePort->hasPalette             = 0;
    thePort->range                          = 2;
    thePort->spectrumStart          = 0;
    thePort->spectrumEnd            = 1;
    printf("Using B&W color map\n");
  }

  thePort->PixelRatio = 1;

  /* init pointers to basic drawing functions */
  thePort->Move                   = IFMove;
  thePort->Draw                   = IFDraw;
  thePort->Polyline               = IFPolyline;
  thePort->InversePolyline= IFInversePolyline;
  thePort->Polygon                = IFPolygon;
  thePort->ShadedPolygon  = IFShadedPolygon;
  thePort->InversePolygon = IFInversePolygon;
  thePort->ErasePolygon   = IFErasePolygon;
  thePort->Polymark               = IFPolymark;
  thePort->InvPolymark    = IFInvPolymark;
  thePort->DrawText               = IFDrawText;
  thePort->CenteredText   = IFCenteredText;
  thePort->ClearViewPort  = IFClearViewPort;

  /* init pointers to set functions */
  thePort->SetLineWidth           = IFSetLineWidth;
  thePort->SetTextSize            = IFSetTextSize;
  thePort->SetMarker                      = IFSetMarker;
  thePort->SetMarkerSize          = IFSetMarkerSize;
  thePort->SetColor                       = IFSetColor;
  thePort->SetPaletteEntry        = IFSetPaletteEntry;
  thePort->SetNewPalette          = IFSetNewPalette;

  /* init pointers to miscellaneous functions */
  thePort->GetPaletteEntry        = IFGetPaletteEntry;
  thePort->Flush                          = IFFlush;
  thePort->PlotPixelBuffer    = IFPlotPixelBuffer;
}

/*==========================================================================*/
/*																			*/
/* The windows controls (statistics, toolbox, ...)							*/
/*																			*/
/*==========================================================================*/

/****************************************************************************/
/*
   InitControls - Load bitmaps from "tooli" include files

   SYNOPSIS:
   int InitControls (Window win);

   PARAMETERS:
   .  win -

   DESCRIPTION:
   This function loads bitmaps from "tooli" include files.

   RETURN VALUE:
   int
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

int NS_PREFIX InitControls (Window win)
{
  int i;

  /* make bitmaps from included files */
  tools[0] = XCreateBitmapFromData(display,win,(const char*)tool0_bits,tool0_width,tool0_height);
  tools[1] = XCreateBitmapFromData(display,win,(const char*)tool1_bits,tool1_width,tool1_height);
  tools[2] = XCreateBitmapFromData(display,win,(const char*)tool2_bits,tool2_width,tool2_height);
  tools[3] = XCreateBitmapFromData(display,win,(const char*)tool3_bits,tool3_width,tool3_height);
  tools[4] = XCreateBitmapFromData(display,win,(const char*)tool4_bits,tool4_width,tool4_height);
  tools[5] = XCreateBitmapFromData(display,win,(const char*)tool5_bits,tool5_width,tool5_height);
  tools[6] = XCreateBitmapFromData(display,win,(const char*)tool6_bits,tool6_width,tool6_height);

  /* set cursor id's */
  cursors[0] = XCreateFontCursor(display,XC_top_left_arrow);
  cursors[1] = XCreateFontCursor(display,XC_X_cursor);
  cursors[2] = XCreateFontCursor(display,XC_crosshair);
  cursors[3] = XCreateFontCursor(display,XC_circle);
  cursors[4] = XCreateFontCursor(display,XC_hand1);
  cursors[5] = XCreateFontCursor(display,XC_heart);
  cursors[6] = XCreateFontCursor(display,XC_box_spiral);

  /* make stipples for shading */
  if (!true_color)
    for (i = 0; i < NO_PATTERNS; i++)
      pattern[i] = XCreateBitmapFromData(display, win, (const char*)pattern_data[i],
                                         PATTERN_SIZE, PATTERN_SIZE);
  return 0;
}

/****************************************************************************/
/*
   WhichTool - check if a tool has been hit

   SYNOPSIS:
   int WhichTool (GraphWindow *gwin, int x, int y, int *tool);

   PARAMETERS:
   .  gwin -
   .  x - pointer to the array of coordinates
   .  y - pointer to the array of coordinates
   .  tool - return tool number

   RETURN VALUE:
   INT
   .n     0 if no tool has been hit
   .n     1 if a valid tool number is returned in tool
 */
/****************************************************************************/

INT NS_PREFIX WhichTool (WINDOWID win, const INT mouse[2], INT *tool)
{
  GraphWindow *gwin;
  int xx,yy,w,h,i;

  gwin = (GraphWindow *)(win);

  xx = gwin->window_width-NTOOLS*CONTROLSIZE;
  yy = gwin->window_height-CONTROLSIZE;
  w = mouse[0]-xx;
  h = mouse[1]-yy;

  if ((h<0)||(h>=CONTROLSIZE)||(w<0)) return(0);
  i = w/CONTROLSIZE;
  if ((i<0)||(i>=NTOOLS)) return(0);
  *tool = i;
  return(1);
}

/****************************************************************************/
/*																			*/
/* Function:   DrawRegion													*/
/*																			*/
/* Purpose:    check if draw region has been hit							*/
/*																			*/
/* Input:	   GraphWindow *gwin the window                                                                 */
/*			   int x			 where										*/
/*			   int y														*/
/*																			*/
/* Output:	   0: draw region has not been hit								*/
/*			   1: draw region has been hit									*/
/*																			*/
/****************************************************************************/

int NS_PREFIX DrawRegion (GraphWindow *gwin, int x, int y)
{
  int xx,yy;

  xx = gwin->window_width;
  yy = gwin->window_height-CONTROLSIZE;

  if ((x<0)||(x>xx)) return(0);
  if ((y<0)||(y>yy)) return(0);
  return(1);
}

void NS_PREFIX DrawInfoBox (WINDOWID win, const char *info)
{
  GraphWindow *gwin;
  XRectangle rect;
  int h,x,y,w;
  int ts;

  /*char buffer[INFO_SIZE];
     strncpy(buffer,info,INFO_SIZE);
     buffer[INFO_LEN] = '\0';*/

  gw = gwin = (GraphWindow *)(win);

  /* set clipping region to info box */
  rect.x = x = 0;
  rect.y = y = gwin->window_height-CONTROLSIZE;
  rect.width = w = gwin->window_width-NTOOLS*CONTROLSIZE;
  rect.height = h = CONTROLSIZE-1;
  XSetClipRectangles(display,gwin->gc,0,0,&rect,1,YSorted);

  /* erase info box */
  XClearArea(display,gwin->win,x,y,w,h,0);

  /* draw info string */
  IFSetColor(X11OutputDevice->black);
  x = (gwin->window_width-NTOOLS*CONTROLSIZE)/2;
  y = gwin->window_height-(CONTROLSIZE/2);
  w = XTextWidth(gwin->font_info,info,strlen(info));
  ts = gwin->font_height;
  XDrawString(display,gwin->win,gwin->gc,x-w/2,y+ts/2,info,strlen(info));
  if (!gwin->backing_store)
    XDrawString(display,gwin->pixmap,gwin->gc,x-w/2,y+ts/2,info,strlen(info));
  /* restore clipping region */
  rect.x = 0;
  rect.y = 0;
  rect.width = gwin->window_width;
  rect.height = gwin->window_height-CONTROLSIZE-1;
  XSetClipRectangles(display,gwin->gc,0,0,&rect,1,YSorted);
}

/*==========================================================================*/
/*																			*/
/* Outputdevice functions													*/
/*																			*/
/*==========================================================================*/

/****************************************************************************/
/*																			*/
/* Function:   WhichGW														*/
/*																			*/
/* Purpose:    return GraphWindow with a given window						*/
/*																			*/
/* Input:	   Window win	  window to find								*/
/*																			*/
/* Output:	   GraphWIndow *  if found										*/
/*			   NULL                   if not found									*/
/*																			*/
/****************************************************************************/

GraphWindow *NS_PREFIX WhichGW (Window win)
{
  GraphWindow *g;

  for (g=windowList; g!=NULL; g=g->next)
    if (g->win==win) return(g);
  return(NULL);
}


/****************************************************************************/
/*																			*/
/* Function:   SetCurrentGW                                                                                             */
/*																			*/
/* Purpose:    set active output window                                                                         */
/*																			*/
/* Input:	   GraphWindow *g												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void NS_PREFIX SetCurrentGW (GraphWindow *g)
{
  gw = g;
}

/****************************************************************************/
/*																			*/
/* Function:   GraphOpen													*/
/*																			*/
/* Purpose:    Open a graph window and initialize a GraphWindow data str.	*/
/*																			*/
/* Input:	   GraphWindow *gw												*/
/*			   char *window_name											*/
/*																			*/
/* Output:	   0: OK														*/
/*			   1: error, could not complete                                                                 */
/*																			*/
/****************************************************************************/

static int GraphOpen (GraphWindow *gw, char *window_name, int x, int y, int width, int height)
{
  unsigned int border_width = DEFAULTBORDER;
  char *icon_name = (char *) GRAPHICONNAME;
  const char *Fontname;
  const char *ug_name = RESOURCENAME;
  unsigned long valuemask = 0;
  static char dash_list[] = {12,24};
  XGCValues values;
  XSizeHints size_hints;
  XWMHints wm_hints;
  XClassHint class_hints;
  XEvent report;

  /* first load font to know size of window */
  Fontname = XGetDefault(display,ug_name,"viewfont");
  if (Fontname==NULL)
    strcpy(gw->font_name,DEFAULTFONT);
  else
    strcpy(gw->font_name,Fontname);
  if ( (gw->font_info=XLoadQueryFont(display,gw->font_name))==NULL)
  {
    /* font could not be loaded, try default font */
    strcpy(gw->font_name,DEFAULTFONT);
    if ( (gw->font_info=XLoadQueryFont(display,gw->font_name))==NULL)
    {
      fprintf(stderr,"%s could not load font %s\n",prog_name,gw->font_name);
      return(1);
    }
  }
  gw->font_ascent = gw->font_info->ascent;
  gw->font_height = gw->font_info->ascent+gw->font_info->descent;
  gw->font_width = gw->font_info->max_bounds.width;

  /* open a window */
  gw->window_width = width;
  gw->window_height = height;
  gw->window_x = x;
  gw->window_y = y;
  gw->win = XCreateSimpleWindow(display,RootWindow(display,screen_num),
                                x,y,width,height,border_width,ctab[X11OutputDevice->black].pixel,
                                ctab[X11OutputDevice->white].pixel);

  /* use backing store for window if possible */
  {
    XSetWindowAttributes attr;
    unsigned long mask;

    if (gw->backing_store = DoesBackingStore(DefaultScreenOfDisplay(display))) {
      attr.backing_store = Always;
      mask = CWBackingStore;
      XChangeWindowAttributes(display, gw->win, mask, &attr);
    }
  }

  /* set the window's color map (may be a private one) */
  XSetWindowColormap(display, gw->win, our_cmap);

  /* create a pixmap with white background	*/
  if (!gw->backing_store) {
    gw->pixmap_depth=DefaultDepth(display,screen_num);
    gw->pixmap_width=gw->window_width;
    gw->pixmap_height=gw->window_height;
    gw->pixmap=XCreatePixmap(display,gw->win,gw->pixmap_width,gw->pixmap_height,gw->pixmap_depth);
    gw->gc=XCreateGC(display,gw->win, valuemask,&values);
    XSetForeground(display,gw->gc,ctab[X11OutputDevice->white].pixel);
    XFillRectangle(display,gw->pixmap,gw->gc,0,0,gw->pixmap_width,gw->pixmap_height);
    XFreeGC(display,gw->gc);
  }

  /* generate icon, needed for properties */
  gw->icon_pixmap = XCreateBitmapFromData(display,gw->win,(const char*)view_icon_bits,view_icon_width,
                                          view_icon_height);

  /* set standard properties */
  size_hints.flags = USPosition | USSize | PMinSize;
  size_hints.min_width = DEFAULTMINX;
  size_hints.min_height = DEFAULTMINY;
  if (XStringListToTextProperty(&window_name,1,&(gw->window_name))==0)
  {
    fprintf(stderr,"%s: structure alloc for window_name failed.\n",prog_name);
    exit(-1);
  }
  if (XStringListToTextProperty(&icon_name,1,&(gw->icon_name))==0)
  {
    fprintf(stderr,"%s: structure alloc for icon_name failed.\n",prog_name);
    exit(-1);
  }
  wm_hints.initial_state = NormalState;
  wm_hints.input = True;
  wm_hints.icon_pixmap = gw->icon_pixmap;
  wm_hints.flags = StateHint | IconPixmapHint | InputHint;
  class_hints.res_name = prog_name;
  class_hints.res_class = (char*) ug_name;
  XSetWMProperties(display,gw->win,&gw->window_name,&gw->icon_name,if_argv,if_argc,
                   &size_hints,&wm_hints,&class_hints);

  /* select event types that will be received */
  XSelectInput(display,gw->win,EnterWindowMask|ExposureMask|KeyPressMask|PointerMotionMask|
               ButtonPressMask|StructureNotifyMask|ButtonReleaseMask|ButtonMotionMask|PointerMotionHintMask);

  /* prepare graphics context */
  gw->gc = XCreateGC(display,gw->win,valuemask,&values);
  XSetFont(display,gw->gc,gw->font_info->fid);
  XSetForeground(display,gw->gc,ctab[X11OutputDevice->black].pixel);
  XSetBackground(display,gw->gc,ctab[X11OutputDevice->white].pixel);
  XSetLineAttributes(display,gw->gc,0,LineSolid,CapRound,JoinRound);
  XSetDashes(display,gw->gc,0,dash_list,2);

  /* now map the window */
  XMapWindow(display,gw->win);

  /* create a region to accumulate update region */
  gw->region = XCreateRegion();

  /* wait here for window to become visible, i.e. wait for			*/
  /* first expose event. I'm not sure that this is the best way, but  */
  /* if opened from a script we will not be in the event loop for       */
  /* a long time.                                                                                                       */
  while (1) {
    XNextEvent(display,&report);
    if (report.xexpose.window==gw->win) break;
  }

  /* und tschuess */
  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  OpenDocumentWindow											*/
/*																			*/
/* Purpose:   open a X window and associate it with theView                             */
/*																			*/
/* Input:	  char *Windowtitle: window title								*/
/*			  UGWINDOW *theUgWindow: view for that window					*/
/*																			*/
/* Output:	  void *: pointer the window struct or							*/
/*					  NULL if an error occured								*/
/*																			*/
/****************************************************************************/

static WINDOWID X11_OpenOutput (const char *title, INT rename, INT x, INT y, INT width, INT height, INT *Global_LL, INT *Global_UR, INT *Local_LL, INT *Local_UR, INT *error)
{
  GraphWindow *gw;

  *error = 0;

  /* create GraphWindow structure */
  gw = (GraphWindow *) malloc(sizeof(GraphWindow));
  if (gw==NULL) return((WINDOWID) NULL);
  gw->next = windowList;
  windowList = gw;

  /* x,y,width,height is in coordinate system with x to the right and y up */
  /* clip requested window against screen and return clipped values		 */
  if (x<0) x=0;
  if (y<0) y=0;
  if ((x)+(width)>display_width) width = display_width-x;
  if ((y)+(height)+CONTROLSIZE>display_height) height = display_height-y-CONTROLSIZE;
  if ((width<DEFAULTMINX)||(height<DEFAULTMINY))
  {
    *error=1;
    UserWriteF("X11_OpenOutput(): ERROR inappropriate window size\n");
    return(0);
  }

  /* open new window */
  /* the following (char *) cast is ugly, but can not really be avoided, since the
          XStringListToTextProperty in GraphOpen wansts a (char **) pointer */
  if (GraphOpen(gw,(char *)title,x,display_height-y-height-CONTROLSIZE,width+1,height+CONTROLSIZE+2)>0)
  {
    *error=1;
    UserWriteF("X11_OpenOutput(): ERROR failed to open window\n");
    return(0);
  }

  /* fill lower left and upper right in the devices coordinate system */
  Global_LL[0] = x; Global_LL[1] = display_height-y;
  Global_UR[0] = Global_LL[0]+width; Global_UR[1] = Global_LL[1]+height;
  Local_LL[0] = 0; Local_LL[1] = height;
  Local_UR[0] = width; Local_UR[1] = 0;

  /* return window ptr */
  return((WINDOWID)gw);
}


/****************************************************************************/
/*																			*/
/* Function:  CloseDocumentWindow											*/
/*																			*/
/* Purpose:   close the X11 window associated with theView					*/
/*																			*/
/* Input:	  VIEW *theView: Port of that View will be inited				*/
/*																			*/
/* Output:	  INT: 0 if all was done well									*/
/*				   1 if an error ocurred									*/
/*																			*/
/****************************************************************************/

static INT X11_CloseOutput (WINDOWID win)
{
  GraphWindow *old,*g;

  /* get window */
  old = (GraphWindow *)(win);

  /* remove from window list */
  if (windowList==old)
    windowList = old->next;
  else
  {
    for (g=windowList; g!=NULL; g=g->next)
      if (g->next==old) break;
    if (g==NULL) return(1);             /* not found */
    g->next = old->next;
  }

  /* close window on screen */
  XDestroyWindow(display,old->win);

  /* free memory */
  free(old);

  /* no error */
  return(0);
}


/****************************************************************************/
/*																			*/
/* Functions: X11_ActivateOutput											*/
/*																			*/
/* Purpose:   activate the window win										*/
/*																			*/
/* Input:	  WINDOWID win													*/
/*																			*/
/* Output:	  0 is OK														*/
/*																			*/
/****************************************************************************/

static INT X11_ActivateOutput (WINDOWID win)
{
  XRectangle rect;

  /* set current output window */
  gw = (GraphWindow *)(win);

  /* set proper clipping region without toolbox */
  rect.x = 0;
  rect.y = 0;
  rect.width = gw->window_width;
  rect.height = gw->window_height-CONTROLSIZE-1;
  XSetClipRectangles(display,gw->gc,0,0,&rect,1,YSorted);

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:   UpdateOutput                                                                                             */
/*																			*/
/* Purpose:    Draws all controls and highlights active tool				*/
/*																			*/
/* Input:	   GraphWindow *gwin											*/
/*																			*/
/* Output:	   0: OK														*/
/*			   1: error, could not complete                                                                 */
/*																			*/
/****************************************************************************/

static INT X11_UpdateOutput (WINDOWID win, INT tool)
{
  int x,y,w,h,i;
  int lw;
  unsigned long fg,bg;
  XGCValues gcv;
  XRectangle rect;
  GraphWindow *gwin;

  IFDEBUG(dev,1)
  printf("Draw Controls\n");
  ENDDEBUG

    gwin = (GraphWindow *) win;

  /* save some values of graphics context */
  XGetGCValues(display, gwin->gc, GCForeground|GCBackground|GCLineWidth, &gcv);
  lw = gcv.line_width;
  fg = gcv.foreground;
  bg = gcv.background;

  /* check window and pixmap size */
  if (!gwin->backing_store) {
    if((gwin->window_width!=gwin->pixmap_width)||(gwin->window_height!=gwin->pixmap_height))
    {                   /* window size has changed, create a new empty pixmap with the new size */
      XFreePixmap(display,gwin->pixmap);
      gwin->pixmap_width=gwin->window_width;
      gwin->pixmap_height=gwin->window_height;
      gwin->pixmap=XCreatePixmap(display,gwin->win,gwin->pixmap_width,gwin->pixmap_height,gwin->pixmap_depth);
      XSetForeground(display,gwin->gc,ctab[X11OutputDevice->white].pixel);
      XFillRectangle(display,gwin->pixmap,gwin->gc,0,0,gwin->pixmap_width,gwin->pixmap_height);
      XSetForeground(display,gwin->gc,ctab[X11OutputDevice->black].pixel);
    }
    else
    {                   /* no changes in size */
      XClearWindow(display,gwin->win);
      XCopyArea(display, gwin->pixmap, gwin->win,gwin->gc, 0, 0,gwin->pixmap_width,gwin->pixmap_height,0,0);
    }
  }

  /* set clipping to whole window */
  rect.x = 0;
  rect.y = 0;
  rect.width = gwin->window_width;
  rect.height = gwin->window_height;
  XSetClipRectangles(display,gwin->gc,0,0,&rect,1,YSorted);

  /* set to black onto white and lw=1, important for copy plane */
  gcv.line_width = 0;
  gcv.foreground = ctab[X11OutputDevice->black].pixel;
  gcv.background = ctab[X11OutputDevice->white].pixel;
  XChangeGC(display, gwin->gc, GCForeground|GCBackground|GCLineWidth, &gcv);

  /* clear controls area including line */
  x = 0;
  y = gwin->window_height-1-CONTROLSIZE;
  w = gwin->window_width;
  h = CONTROLSIZE;
  XClearArea(display,gwin->win,x,y,w,h,0);

  /* draw tools */
  y = gwin->window_height-CONTROLSIZE;
  w = h = CONTROLSIZE;
  for (i=0; i<NTOOLS; i++)
  {
    x = gwin->window_width-(NTOOLS-i)*CONTROLSIZE;
    XCopyPlane(display,tools[i],gwin->win,gwin->gc,0,0,w,h,x,y,1);
    if (!gwin->backing_store)
      XCopyPlane(display,tools[i],gwin->pixmap,gwin->gc,0,0,w,h,x,y,1);
  }

  /* highlite active tool */
  i = tool;
  if ((i<0)||(i>=NTOOLS)) i=0;
  x = gwin->window_width-(NTOOLS-i)*CONTROLSIZE;
  y = gwin->window_height-CONTROLSIZE;
  w = h = CONTROLSIZE;
  gcv.foreground = ctab[X11OutputDevice->white].pixel;
  gcv.background = ctab[X11OutputDevice->black].pixel;
  XChangeGC(display, gwin->gc, GCForeground|GCBackground, &gcv);
  XCopyPlane(display,tools[i],gwin->win,gwin->gc,0,0,w,h,x,y,1);
  gcv.foreground = ctab[X11OutputDevice->black].pixel;
  gcv.background = ctab[X11OutputDevice->white].pixel;
  XChangeGC(display, gwin->gc, GCForeground|GCBackground, &gcv);

  /* set cursor */
  XDefineCursor(display,gwin->win,cursors[i]);

  /* draw separation line */
  y = gwin->window_height-CONTROLSIZE-1;
  w = gwin->window_width;
  XDrawLine(display,gwin->win,gwin->gc,0,y,w-1,y);
  if (!gwin->backing_store)
    XDrawLine(display,gwin->pixmap,gwin->gc,0,y,w-1,y);

  /* restore GC */
  gcv.line_width = lw;
  gcv.foreground = fg;
  gcv.background = bg;
  XChangeGC(display, gwin->gc, GCForeground|GCBackground|GCLineWidth, &gcv);

  /* restore clipping region */
  rect.x = 0;
  rect.y = 0;
  rect.width = gwin->window_width;
  rect.height = gwin->window_height-CONTROLSIZE-1;
  XSetClipRectangles(display,gwin->gc,0,0,&rect,1,YSorted);

  return(0);
}

/****************************************************************************/
/*
   InitOutputDevice - Install output device 'x11'

   SYNOPSIS:
   OUTPUTDEVICE *InitXOutputDevice (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function installs output device 'x11'.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      pointer to
   .n      NULL if an error occured.
 */
/****************************************************************************/

OUTPUTDEVICE *NS_PREFIX InitXOutputDevice (void)
{
  /* create output device */
  if ((X11OutputDevice=CreateOutputDevice("screen"))==NULL) return(NULL);

  /* init output device 'x11' */
  X11OutputDevice->OpenOutput  = X11_OpenOutput;
  X11OutputDevice->CloseOutput  = X11_CloseOutput;
  X11OutputDevice->ActivateOutput  = X11_ActivateOutput;
  X11OutputDevice->UpdateOutput  = X11_UpdateOutput;

  X11OutputDevice->v.locked = 1;
  InitXPort(X11OutputDevice);

  printf("output device 'screen' for x11 window manager created\n");

  return(X11OutputDevice);
}


/****************************************************************************/
/*
   MousePosition - Get current mouse position

   SYNOPSIS:
   void MousePosition (INT *ScreenPoint);

   PARAMETERS:
   .  ScreenPoint - return result in this vector

   DESCRIPTION:
   This function gets current mouse position.


   RETURN VALUE:
   void
 */
/****************************************************************************/

void NS_PREFIX MousePosition (INT *ScreenPoint)
{
  XEvent report;
  int where_x,where_y;
  int exitLoop;
  int root_x,root_y;
  Window root,child;
  unsigned int keys_buttons;

  /* if a mouse pos is stored, return it */
  if (StoredMousePos)
  {
    ScreenPoint[0] = MouseX;
    ScreenPoint[1] = MouseY;
    StoredMousePos = 0;
    return;
  }

  /* get events until the mouse has been moved or button released */
  exitLoop=0;
  while (!exitLoop)
  {

    /* get next event, all events except button and motion are discarded */
    XNextEvent(display,&report);

    /* examine event */
    switch (report.type)
    {
    case ButtonRelease :
      MoveMouse = 0;
      where_x = report.xbutton.x;
      where_y = report.xbutton.y;
      exitLoop=1;
      break;

    case MotionNotify :
      while (XCheckMaskEvent(display,ButtonMotionMask,&report)) ;
      if (!XQueryPointer(display,report.xmotion.window,&root,&child,
                         &root_x,&root_y,&where_x,&where_y,&keys_buttons))
        break;
      exitLoop=1;
      break;

    default :
      break;
    }
  }
  ScreenPoint[0] = where_x;
  ScreenPoint[1] = where_y;

  return;
}

/****************************************************************************/
/*
   MouseStillDown - Determine if mouse button is still pressed


   SYNOPSIS:
   INT MouseStillDown (void);

   PARAMETERS:
   no parameters

   DESCRIPTION:
   This function returns true (1) if the mouse button is still pressed.
   The function should only be called after a button pressed event has been
   reported.


   RETURN VALUE:
   INT
   .n 0 if mouse button has been released
   .n 1 if mouse button is still pressed
 */
/****************************************************************************/

static Bool callback (Display *d, XEvent *report, char *arg)
{
  return(True);
}

INT NS_PREFIX MouseStillDown (void)
{
  XEvent report;
  int where_x,where_y;
  int root_x,root_y;
  Window root,child;
  unsigned int keys_buttons;

  /* if button has been released already */
  if (MoveMouse==0)
  {
    MoveMouse = 1;
    return(0);
  }

  /* get next event, all events except button and motion are discarded */
  if (!XCheckIfEvent(display,&report,callback,NULL)) return(1);

  /* examine event */
  switch (report.type)
  {
  case ButtonRelease :
    MoveMouse = 1;
    return(0);

  case MotionNotify :
    while (XCheckMaskEvent(display,ButtonMotionMask,&report)) ;
    if (!XQueryPointer(display,report.xmotion.window,&root,&child,
                       &root_x,&root_y,&where_x,&where_y,&keys_buttons))
      break;
    StoredMousePos = 1;
    MouseX = where_x;
    MouseY = where_y;
    return(1);

  default :
    break;
  }

  return (1);
}
