// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MacGraph.c													*/
/*																			*/
/* Purpose:   handling windows for graphical output                                             */
/*																			*/
/* Author:	  Peter Bastian/Henrik Rentz-Reichert							*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   26.03.92 begin, ug version 2.0								*/
/*			  13.02.95 begin, ug version 3.0								*/
/*																			*/
/* Remarks:   former MacDocumentWindows.c/.h								*/
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment mif
#endif

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include "MWCW.cmdlinedefs"

#ifndef __MWCW__        /* don't need that: included MacHeadersPPC */

/* Macintosh toolbox includes */
#include <Types.h>
#include <Memory.h>
#include <Quickdraw.h>
#include <Fonts.h>
#include <Events.h>
#include <Menus.h>
#include <Windows.h>
#include <Palettes.h>
#include <TextEdit.h>
#include <Dialogs.h>
#include <OSUtils.h>
#include <ToolUtils.h>
#include <SegLoad.h>
#include <desk.h>
#include <osevents.h>
#include <Packages.h>
#include <Files.h>

#endif

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
#include "MacMain.m"
#include "MacSurface.h"
#include "MacGraph.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#ifdef __MWCW__
#define Str255          ConstStr255Param
#define qdgray          &qd.black
#else
#define qdgray          &qd.black
#endif

#define GET_CHOSEN_TOOL(val,rect,m)     {\
    if (((int)m[1]<rect.bottom-16) || ((int)m[1]>rect.bottom) || ((int)m[0]-(rect.right-119)<0) \
        || (((val=floor(m[0]-(rect.right-119))/15)<0) || (val>=nboftools))) val =       NO_TOOL_CHOSEN;}

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static GRAPH_WINDOW *currgw;                    /* active output					*/
static GRAPH_WINDOW *windowList=NULL;   /* list of our windows				*/
static OUTPUTDEVICE *MacOutputDevice;   /* outputdevice that has been initi */

static long guiHeapSize=32000;  /* size of system heap						*/
static HEAP *guiHeap=NULL;              /* system heap structure					*/
static char buffer[255];

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* Functions used in InitMacPort											*/
/*																			*/
/* Purpose:   basic drawing functions by QuickDraw							*/
/*																			*/
/****************************************************************************/

static void MacMove (SHORT_POINT point)
{
  MoveTo (point.x, point.y);
}

static void MacDraw (SHORT_POINT point)
{
  LineTo (point.x, point.y);
}

static void MacPolyline (SHORT_POINT *points, INT n)
{
  int i;

  if (n<2) return;

  MoveTo (points[0].x, points[0].y);
  for (i=1; i<n; i++)
    LineTo (points[i].x, points[i].y);
}

static void MacInversePolyline (SHORT_POINT *points, INT n)
{
  int i;

  if (n<2) return;

  PenMode(patXor);
  PenPat(qdgray);
  MoveTo (points[0].x, points[0].y);
  for (i=1; i<n; i++)
    LineTo (points[i].x, points[i].y);
  PenNormal();
}

static void MacPolygon (SHORT_POINT *points, INT n)
{
  int i;
  PolyHandle myPoly;

  if (n<3) return;

  myPoly = OpenPoly();
  MoveTo (points[0].x, points[0].y);
  for (i=1; i<n; i++) LineTo (points[i].x, points[i].y);
  LineTo (points[0].x, points[0].y);
  ClosePoly();
  PaintPoly(myPoly);
  FramePoly(myPoly);
  KillPoly(myPoly);
}

static void MacInversePolygon (SHORT_POINT *points, INT n)
{
  int i;
  PolyHandle myPoly;

  if (n<3) return;

  myPoly = OpenPoly();
  MoveTo (points[0].x, points[0].y);
  for (i=1; i<n; i++) LineTo (points[i].x, points[i].y);
  LineTo (points[0].x, points[0].y);
  ClosePoly();
  InvertPoly(myPoly);
  KillPoly(myPoly);
}

static void MacErasePolygon (SHORT_POINT *points, INT n)
{
  int i;
  PolyHandle myPoly;
  RGBColor ForeColor, BackColor;

  if (n<3) return;

  GetForeColor(&ForeColor);
  GetBackColor(&BackColor);
  myPoly = OpenPoly();
  MoveTo (points[0].x, points[0].y);
  for (i=1; i<n; i++) LineTo (points[i].x, points[i].y);
  LineTo (points[0].x, points[0].y);
  ClosePoly();
  ErasePoly(myPoly);
  RGBForeColor(&BackColor);
  FramePoly(myPoly);
  KillPoly(myPoly);
  RGBForeColor(&ForeColor);
}

static void MacSetMarkerSize (short s)
{
  currgw->marker_size = s;
  return;
}

static void MacSetMarker (short s)
{
  currgw->marker_id = s;
  return;
}

static void Marker (SHORT_POINT point)
{
  Rect r;
  PolyHandle myPoly;
  short n,s,S,x,y;

  s = currgw->marker_size/2;
  S = 1.41*s;
  n = currgw->marker_id;

  x = point.x;
  y = point.y;
  r.top = y-s; r.bottom = y+s;
  r.left = x-s; r.right = x+s;

  switch (n)
  {
  case EMPTY_SQUARE_MARKER :
    FrameRect(&r);
    break;
  case GRAY_SQUARE_MARKER :
    PenPat(qdgray);
    PaintRect(&r);
    PenNormal();
    break;
  case FILLED_SQUARE_MARKER :
    PaintRect(&r);
    break;
  case EMPTY_CIRCLE_MARKER :
    FrameOval(&r);
    break;
  case GRAY_CIRCLE_MARKER :
    PenPat(qdgray);
    PaintOval(&r);
    PenNormal();
    break;
  case FILLED_CIRCLE_MARKER :
    PaintOval(&r);
    break;
  case EMPTY_RHOMBUS_MARKER :
    MoveTo(x,y+S);
    LineTo(x+S,y);
    LineTo(x,y-S);
    LineTo(x-S,y);
    LineTo(x,y+S);
    break;
  case GRAY_RHOMBUS_MARKER :
    PenPat(qdgray);
    myPoly = OpenPoly();
    MoveTo(x,y+S);
    LineTo(x+S,y);
    LineTo(x,y-s);
    LineTo(x-S,y);
    LineTo(x,y+S);
    ClosePoly();
    PaintPoly(myPoly);
    KillPoly(myPoly);
    PenNormal();
    break;
  case FILLED_RHOMBUS_MARKER :
    myPoly = OpenPoly();
    MoveTo(x,y+S);
    LineTo(x+S,y);
    LineTo(x,y-S);
    LineTo(x-S,y);
    LineTo(x,y+S);
    ClosePoly();
    PaintPoly(myPoly);
    KillPoly(myPoly);
    PenNormal();
    break;
  case PLUS_MARKER :
    MoveTo(x,y+s);
    LineTo(x,y-s);
    MoveTo(x-s,y);
    LineTo(x+s,y);
    break;
  case CROSS_MARKER :
    MoveTo(x-s,y+s);
    LineTo(x+s,y-s);
    MoveTo(x-s,y-s);
    LineTo(x+s,y+s);
    break;
  default :
    break;
  }
}

static void MacPolymark (short n, SHORT_POINT *points)
{
  int i;

  for (i=0; i<n; i++)
    Marker(points[i]);
}

static void InvMarker (SHORT_POINT point)
{
  Rect r;
  PolyHandle myPoly;
  short n,s,S,x,y;

  s = currgw->marker_size/2;
  S = 1.41*s;
  n = currgw->marker_id;

  x = point.x;
  y = point.y;
  r.top = y-s; r.bottom = y+s;
  r.left = x-s; r.right = x+s;

  switch (n)
  {
  case EMPTY_SQUARE_MARKER :
    PenMode(patXor);
    FrameRect(&r);
    break;
  case GRAY_SQUARE_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    PaintRect(&r);
    PenNormal();
    break;
  case FILLED_SQUARE_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    PaintRect(&r);
    break;
  case EMPTY_CIRCLE_MARKER :
    PenMode(patXor);
    FrameOval(&r);
    break;
  case GRAY_CIRCLE_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    PaintOval(&r);
    PenNormal();
    break;
  case FILLED_CIRCLE_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    PaintOval(&r);
    break;
  case EMPTY_RHOMBUS_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    MoveTo(x,y+S);
    LineTo(x+S,y);
    LineTo(x,y-S);
    LineTo(x-S,y);
    LineTo(x,y+S);
    break;
  case GRAY_RHOMBUS_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    myPoly = OpenPoly();
    MoveTo(x,y+S);
    LineTo(x+S,y);
    LineTo(x,y-s);
    LineTo(x-S,y);
    LineTo(x,y+S);
    ClosePoly();
    InvertPoly(myPoly);
    KillPoly(myPoly);
    PenNormal();
    break;
  case FILLED_RHOMBUS_MARKER :
    myPoly = OpenPoly();
    MoveTo(x,y+S);
    LineTo(x+S,y);
    LineTo(x,y-S);
    LineTo(x-S,y);
    LineTo(x,y+S);
    ClosePoly();
    InvertPoly(myPoly);
    KillPoly(myPoly);
    PenNormal();
    break;
  case PLUS_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    MoveTo(x,y+s);
    LineTo(x,y-s);
    MoveTo(x-s,y);
    LineTo(x+s,y);
    break;
  case CROSS_MARKER :
    PenMode(patXor);
    PenPat(qdgray);
    MoveTo(x-s,y+s);
    LineTo(x+s,y-s);
    MoveTo(x-s,y-s);
    LineTo(x+s,y+s);
    break;
  default :
    break;
  }
  PenNormal();
}

static void MacInvPolymark (short n, SHORT_POINT *points)
{
  int i;

  for (i=0; i<n; i++)
    InvMarker(points[i]);
}

static void MacDrawText (const char *s, INT mode)
{
  switch(mode)
  {
  case TEXT_REGULAR : TextMode(srcOr); break;
  case TEXT_INVERSE : TextMode(srcXor); break;
  default : return;
  }
  strcpy(buffer,s);
  DrawString(((Str255)c2pstr(buffer)));
  TextMode(srcOr);
}

static void MacCenteredText (SHORT_POINT point, const char *s, INT mode)
{
  short ts,w;

  ts = currgw->textSize;

  switch(mode)
  {
  case TEXT_REGULAR : TextMode(srcOr); break;
  case TEXT_INVERSE : TextMode(srcXor); break;
  default : return;
  }
  strcpy(buffer,s);
  c2pstr(buffer);
  w = StringWidth((Str255)buffer);
  MoveTo(point.x-w/2,point.y+ts/2);
  DrawString((Str255)buffer);
  TextMode(srcOr);
}

static void MacClearViewPort (void)
{
  Rect PlotRgnRect;

  PlotRgnRect.left        = currgw->Local_LL[0];
  PlotRgnRect.bottom      = currgw->Local_LL[1];
  PlotRgnRect.right       = currgw->Local_UR[0];
  PlotRgnRect.top         = currgw->Local_UR[1];

  EraseRect(&PlotRgnRect);
}

static void MacSetLineWidth (short w)
{
  PenSize(w,w);
}

static void MacSetTextSize (short s)
{
  TextSize(s);
  currgw->textSize = s;
}

static void MacSetColor (long index)
{
  PmForeColor((short)index);
}

static void MacSetPaletteEntry (long index, short r, short g, short b)
{
  RGBColor newColor;
  PaletteHandle myPalette;

  myPalette = GetPalette(MAC_WIN(currgw));
  newColor.red = r;
  newColor.green = g;
  newColor.blue = b;
  SetEntryColor(myPalette,(short) index,&newColor);
  ActivatePalette(MAC_WIN(currgw));
}

static void MacSetNewPalette (long start, long count, short *r, short *g, short *b)
{
  RGBColor newColor;
  PaletteHandle myPalette;
  long i;

  myPalette = GetPalette(MAC_WIN(currgw));
  for (i=start; i<start+count; i++)
  {
    newColor.red = r[i];
    newColor.green = g[i];
    newColor.blue = b[i];
    SetEntryColor(myPalette,(short) i,&newColor);
  }
  ActivatePalette(MAC_WIN(currgw));
}

static void MacGetPaletteEntry (long index, short *r, short *g, short *b)
{
  RGBColor theColor;
  PaletteHandle myPalette;

  myPalette = GetPalette(MAC_WIN(currgw));
  GetEntryColor(myPalette,(short) index,&theColor);
  *r = theColor.red;
  *g = theColor.green;
  *b = theColor.blue;
}

static void MacFlush (void)
{
  return;
}

/****************************************************************************/
/*
   InitMacPort - Implement basic drawing functions by QuickDraw

   SYNOPSIS:
   static void InitMacPort ();

   PARAMETERS:
   .  none

   DESCRIPTION:
   This function implements basic drawing functions by QuickDraw.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void InitMacPort ()
{

  /*	init colors */
  MacOutputDevice->black                          = 1;
  MacOutputDevice->white                          = 0;
  MacOutputDevice->red                            = 254;
  MacOutputDevice->green                          = 170;
  MacOutputDevice->blue                           = 3;
  MacOutputDevice->cyan                           = 66;
  MacOutputDevice->orange                         = 223;
  MacOutputDevice->yellow                         = 192;
  MacOutputDevice->magenta                        = 2;
  MacOutputDevice->hasPalette             = 1;
  MacOutputDevice->range                          = 256;
  MacOutputDevice->spectrumStart          = 3;
  MacOutputDevice->spectrumEnd            = 254;
  MacOutputDevice->signx                          = 1;
  MacOutputDevice->signy                          = -1;

  /* init pointers to basic drawing functions */
  MacOutputDevice->Move                           = MacMove;
  MacOutputDevice->Draw                           = MacDraw;
  MacOutputDevice->Polyline                       = MacPolyline;
  MacOutputDevice->InversePolyline        = MacInversePolyline;
  MacOutputDevice->Polygon                        = MacPolygon;
  MacOutputDevice->InversePolygon         = MacInversePolygon;
  MacOutputDevice->ErasePolygon           = MacErasePolygon;
  MacOutputDevice->Polymark                       = MacPolymark;
  MacOutputDevice->InvPolymark            = MacInvPolymark;
  MacOutputDevice->DrawText                       = MacDrawText;
  MacOutputDevice->CenteredText           = MacCenteredText;
  MacOutputDevice->ClearViewPort          = MacClearViewPort;

  /* init pointers to set functions */
  MacOutputDevice->SetLineWidth           = MacSetLineWidth;
  MacOutputDevice->SetTextSize            = MacSetTextSize;
  MacOutputDevice->SetMarker                      = MacSetMarker;
  MacOutputDevice->SetMarkerSize          = MacSetMarkerSize;
  MacOutputDevice->SetColor                       = MacSetColor;
  MacOutputDevice->SetPaletteEntry        = MacSetPaletteEntry;
  MacOutputDevice->SetNewPalette          = MacSetNewPalette;

  /* init pointers to miscellaneous functions */
  MacOutputDevice->GetPaletteEntry        = MacGetPaletteEntry;
  MacOutputDevice->Flush                          = MacFlush;
}

/****************************************************************************/
/*																			*/
/* Function:   WhichGW														*/
/*																			*/
/* Purpose:    return GRAPH_WINDOW with a given window						*/
/*																			*/
/* Input:	   Window win	  window to find								*/
/*																			*/
/* Output:	   GraphWIndow *  if found										*/
/*			   NULL                   if not found									*/
/*																			*/
/****************************************************************************/

GRAPH_WINDOW *WhichGW (WindowPtr theWindow)
{
  GRAPH_WINDOW *g;

  for (g=windowList; g!=NULL; g=g->next)
    if (MAC_WIN(g)==theWindow) return(g);
  return(NULL);
}


/****************************************************************************/
/*																			*/
/* Function:   SetCurrentGW                                                                                             */
/*																			*/
/* Purpose:    set active output window                                                                         */
/*																			*/
/* Input:	   GRAPH_WINDOW *g												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void SetCurrentGW (GRAPH_WINDOW *g)
{
  currgw = g;
}

/****************************************************************************/
/*																			*/
/* Function:  DrawToolBox													*/
/*																			*/
/* Purpose:   draw Tool box in theWindow									*/
/*																			*/
/* Input:	  DOC_WIN *theWindow											*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

static void DrawToolBox (GRAPH_WINDOW *gw, INT tool)
{
  Rect r,dstRect;
  GrafPtr myPort;
  PicHandle toolBox;
  WindowPtr theWindow;

  theWindow = MAC_WIN(gw);
  myPort = (GrafPtr) theWindow;

  SetPort(myPort);
  r = myPort->portRect;

  toolBox = GetPicture(TOOLBOX_RSRC_ID);
  if (toolBox!=NULL)
  {
    SetRect(&dstRect,r.right-120,r.bottom-15,r.right-15,r.bottom);
    DrawPicture(toolBox,&dstRect);
  }

  SetRect(&dstRect,r.right-119+tool*15,r.bottom-14,r.right-104+tool*15,r.bottom);
  InvertRect(&dstRect);
}

/****************************************************************************/
/*																			*/
/* Function:  DrawInfoBox													*/
/*																			*/
/* Purpose:   draw info box													*/
/*																			*/
/* Input:	  WINDOWID win, info string										*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

void DrawInfoBox (WINDOWID win, const char *info)
{
  WindowPtr theMacWindow;
  Rect r,box,WinRect;
  GrafPtr myPort;
  Rect myClipRect;
  char buffer[INFO_SIZE];
  short w;
  GRAPH_WINDOW *gw;

  gw = (GRAPH_WINDOW *)win;
  theMacWindow = MAC_WIN(gw);
  WinRect = theMacWindow->portRect;

  strncpy(buffer,info,INFO_SIZE);
  buffer[INFO_LEN] = '\0';

  /* set port */
  myPort = (GrafPtr) theMacWindow;
  SetPort(myPort);

  PmForeColor(1);               /* black */
  r = myPort->portRect;

  /* info box */
  SetRect(&box,r.left+1,r.bottom-14,r.right-121,r.bottom);

  /* set clipping region to info box */
  ClipRect(&box);

  EraseRect(&box);

  c2pstr(buffer);
  w = StringWidth((Str255) buffer);
  MoveTo((box.left+box.right)/2-w/2,box.bottom-3);
  TextSize(9);
  DrawString((Str255) buffer);

  /* adjust clipping rectangle again to plottable area */
  myClipRect.left   = gw->Local_LL[0];
  myClipRect.right  = gw->Local_UR[0];
  myClipRect.bottom = gw->Local_LL[1];
  myClipRect.top    = gw->Local_UR[1];

  ClipRect(&myClipRect);
}

/****************************************************************************/
/*																			*/
/* Functions: ActivateMacWin												*/
/*																			*/
/* Purpose:   process activate event for graph window						*/
/*																			*/
/* Input:	  VIEW *theView: pointer to corresponding view					*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT ActivateMacWin (GRAPH_WINDOW *gw, INT tool)
{
  WindowPtr MacWindow;

  MacWindow = MAC_WIN(gw);

  DrawGrowIcon(MacWindow);
  SetPort(MacWindow);

  switch (tool)
  {
  case arrowTool :         SetMyCursor(arrowCurs);         break;
  case crossTool :         SetMyCursor(crossCurs);         break;
  case choiceTool :        SetMyCursor(choiceCurs);        break;
  case circleTool :        SetMyCursor(circleCurs);        break;
  case handTool :          SetMyCursor(handCurs);          break;
  case heartTool :         SetMyCursor(heartCurs);         break;
  case gnoedelTool :       SetMyCursor(gnoedelCurs);       break;
  }

  return (0);
}

/****************************************************************************/
/*																			*/
/* Functions: GetTool														*/
/*																			*/
/* Purpose:   get tool from mouselocation									*/
/*																			*/
/* Input:	  WindowPtr *theWindow: get tool of THIS MacWindow				*/
/*			  INT *MouseLocation:	                                                                                */
/*																			*/
/* Output:	  INT *ChosenToolPtr: chosen tool								*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/****************************************************************************/

INT WhichTool (WINDOWID win, const INT mp[2], INT *ChosenToolPtr)
{
  INT i;
  Rect WinRect;
  WindowPtr theWindow;

  theWindow = MAC_WIN((GRAPH_WINDOW*)win);

  WinRect = theWindow->portRect;

  GET_CHOSEN_TOOL(i,WinRect,mp);

  *ChosenToolPtr = i;

  if (i==NO_TOOL_CHOSEN)
    return (0);

  return(1);
}

/****************************************************************************/
/*																			*/
/* Functions: GrowGraphWindow												*/
/*																			*/
/* Purpose:   grow a graph window											*/
/*																			*/
/* Input:	  WindowPtr *theWindow: window to update						*/
/*			  EventRecord *theEvent: pointer to keyboard event				*/
/*																			*/
/* return:	  POS_CHANGE													*/
/*			  NO_POS_CHANGE													*/
/*																			*/
/****************************************************************************/

INT GrowGraphWindow (GRAPH_WINDOW *gw, EventRecord *theEvent, DOC_GROW_EVENT *docGrow)
{
  long growResult;
  Rect r,sizeRect,beforeRect,afterRect;
  WindowPtr theWindow;
  Rect myClipRect;

  theWindow = MAC_WIN(gw);
  SetPort(theWindow);

  /* grow window */
  beforeRect = theWindow->portRect;
  SetRect(&sizeRect,GRAPHWIN_MINSIZE,GRAPHWIN_MINSIZE,SCREEN_WIDTH-2*MARGIN_TO_SCREEN,SCREEN_HEIGHT-2*MARGIN_TO_SCREEN-MENU_BAR);
  growResult = GrowWindow(theWindow,theEvent->where,&sizeRect);
  if (growResult!=0)
  {
    /* actually change windows size */
    SizeWindow(theWindow,LoWrd(growResult),HiWrd(growResult),true);
    afterRect = theWindow->portRect;

    /* make the new regions invalid */
    if (afterRect.right>beforeRect.right)
    {
      SetRect(&r,beforeRect.right-15,beforeRect.top,beforeRect.right,beforeRect.bottom-15);
      InvalRect(&r);
    }
    else if (afterRect.right<beforeRect.right)
    {
      SetRect(&r,afterRect.right-15,afterRect.top,afterRect.right,afterRect.bottom-15);
      InvalRect(&r);
    }
    if (afterRect.bottom>beforeRect.bottom)
    {
      SetRect(&r,beforeRect.left,beforeRect.bottom-15,beforeRect.right,beforeRect.bottom);
      InvalRect(&r);
    }
    else if (afterRect.bottom<beforeRect.bottom)
    {
      SetRect(&r,afterRect.left,afterRect.bottom-15,afterRect.right,afterRect.bottom);
      InvalRect(&r);
    }

    /* store and report new size */
    docGrow->Global_LL[0] = gw->Global_LL[0];
    docGrow->Global_LL[1] = gw->Global_LL[1] += afterRect.bottom - SCROLL_BAR - gw->Local_LL[1];
    docGrow->Global_UR[0] = gw->Global_UR[0] += afterRect.right  - SCROLL_BAR - gw->Local_UR[0];
    docGrow->Global_UR[1] = gw->Global_UR[1];

    docGrow->Local_LL[0] = gw->Local_LL[0] = afterRect.left;
    docGrow->Local_LL[1] = gw->Local_LL[1] = afterRect.bottom - SCROLL_BAR;
    docGrow->Local_UR[0] = gw->Local_UR[0] = afterRect.right  - SCROLL_BAR;
    docGrow->Local_UR[1] = gw->Local_UR[1] = afterRect.top;

    ClipRect(&afterRect);
    EraseRect(&afterRect);

    /* adjust clipping rectangle */
    myClipRect.left   = gw->Local_LL[0];
    myClipRect.right  = gw->Local_UR[0];
    myClipRect.bottom = gw->Local_LL[1];
    myClipRect.top    = gw->Local_UR[1];

    ClipRect(&myClipRect);

    return (POS_CHANGE);
  }
  else
    return (NO_POS_CHANGE);
}

/****************************************************************************/
/*																			*/
/* Functions: DragGraphWindow												*/
/*																			*/
/* Purpose:   drag window of a view                                                                             */
/*																			*/
/* Input:	  VIEW *theView: the view										*/
/*			  void *Data: the event                                                                                 */
/*																			*/
/* Output:	  INT *UgWindowPositionChanges: changes of the position of the w*/
/*																			*/
/* return:	  POS_CHANGE													*/
/*			  NO_POS_CHANGE													*/
/*																			*/
/****************************************************************************/

INT DragGraphWindow (GRAPH_WINDOW *gw, EventRecord *theEvent, DOC_DRAG_EVENT *docDrag)
{
  WindowPtr theWindow;
  int TopOld, BottomOld, LeftOld, RightOld, Left, Right, Top, Bottom, DelLeft, DelRight, DelTop, DelBottom;
  Rect theDragRect;
  Point MouseLoc;

  theWindow = MAC_WIN(gw);

  /* store old corners of the window */
  TopOld     = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.top;
  BottomOld  = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.bottom;
  LeftOld    = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.left;
  RightOld   = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.right;

  /* set drag rect */
  GetMouse(&MouseLoc);
  LocalToGlobal(&MouseLoc);
  Left   = DragRect()->left       + MouseLoc.h - LeftOld;
  Right  = DragRect()->right      + MouseLoc.h - RightOld;
  Top    = DragRect()->top        + MouseLoc.v - TopOld;
  Bottom = DragRect()->bottom + MouseLoc.v - BottomOld;
  SetRect(&theDragRect,Left,Top,Right,Bottom);

  /* drag window */
  DragWindow (theWindow,theEvent->where,&theDragRect);

  /* report new size */
  Left   = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.left;
  Right  = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.right;
  Top    = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.top;
  Bottom = (*(((WindowPeek)theWindow)->strucRgn))->rgnBBox.bottom;
  DelLeft         = Left  - LeftOld;
  DelRight        = Right - RightOld;
  DelTop          = Top   - TopOld;
  DelBottom       = Bottom- BottomOld;

  if (DelLeft==0 && DelRight==0 && DelTop==0 && DelBottom)
    return (NO_POS_CHANGE);

  docDrag->Global_LL[0] = gw->Global_LL[0] += DelLeft;
  docDrag->Global_LL[1] = gw->Global_LL[1] += DelBottom;
  docDrag->Global_UR[0] = gw->Global_UR[0] += DelRight;
  docDrag->Global_UR[1] = gw->Global_UR[1] += DelTop;

  return (POS_CHANGE);
}

/****************************************************************************/
/*																			*/
/* Function:  GraphOpen                                                                                                         */
/*																			*/
/* Purpose:   create a new Macintosh window                                                             */
/*																			*/
/* Input:	  window title and size                                                                                 */
/*																			*/
/* Output:	  0: ok                                                                                                                 */
/*																			*/
/****************************************************************************/

static INT GraphOpen (GRAPH_WINDOW *gw, const char *title, short h, short v, short dh, short dv)
{
  WindowPtr MacWin;
  GrafPtr myPort;
  PaletteHandle myPalette;
  char WindowTitle[MAXTITLELENGTH];

  /* copy title because of conversion to pascal string */
  strncpy(WindowTitle,title,MAXTITLELENGTH);
  WindowTitle[MAXTITLELENGTH-1] = '\0';

  MacWin = MAC_WIN(gw);

  /* read in resources */
  if (GetNewCWindow(GRAPH_RSRC_ID,(Ptr) MacWin,(WindowPtr) -1)==NULL)
    return (1);

  myPalette = GetNewPalette(PALETTE_RSRC_ID);
  SetPalette(MacWin,myPalette,FALSE);

  /* move and size window */
  myPort = (GrafPtr) MacWin;
  SetPort(myPort);
  MoveWindow(MacWin,h,v,false);
  SizeWindow(MacWin,dh,dv,false);
  SetWTitle(MacWin,c2pstr(WindowTitle));
  ShowWindow(MacWin);
  SelectWindow(MacWin);
  DrawGrowIcon(MacWin);

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  Mac_OpenOutput												*/
/*																			*/
/* Purpose:   open a Mac window                                                                                         */
/*																			*/
/* Input:	  window title, size and lower left + upper right corners		*/
/*																			*/
/* Output:	  >0 if ok														*/
/*																			*/
/****************************************************************************/

static WINDOWID Mac_OpenOutput (
  const char *title,                                                    /* tiltle of the window                 */
  INT x, INT y, INT width, INT height,          /* plot rgn in standard coord.	*/
  INT *Global_LL, INT *Global_UR,                       /* global machine coordinates	*/
  INT *Local_LL, INT *Local_UR,                         /* local machine coordinates	*/
  INT *error)                                                           /* error code					*/
{
  GRAPH_WINDOW *gw;
  Rect myClipRect,*portRect;
  short h,v,dh,dv;

  /* create GRAPH_WINDOW structure and put in list */
  gw = (GRAPH_WINDOW*) GetMem(guiHeap,sizeof(GRAPH_WINDOW),GENERAL_HEAP);
  if (gw==NULL)
  {*error=1; return(0);}
  gw->next   = windowList;
  windowList = gw;

  gw->currTool = arrowTool;

  /* x,y,width,height is in coordinate system with x to the right and y up*/
  /* clip requested window against screen                                                               */

  /* clip lower left */
  x = MAX(x,0);
  y = MAX(y,0);
  /* clip min size */
  width  = MAX(width ,GRAPHWIN_MINSIZE);
  height = MAX(height,GRAPHWIN_MINSIZE);
  /* clip upper right */
  width  = MIN(width,SCREEN_WIDTH-SCROLL_BAR-x);
  height = MIN(height,SCREEN_HEIGHT-MENU_BAR-SCROLL_BAR-y);

  /* dh,dv: width and height of the port(drawable region): contains scroll bars but not the title bar */
  dh = width  + SCROLL_BAR;
  dv = height + SCROLL_BAR;

  /* h,v: origin of the window = upper left corner of the title bar in global coord */
  h = x;
  v = SCREEN_HEIGHT-dv-y;

  /* open new window */
  if (GraphOpen(gw,title,h,v,dh,dv)!=0)
  {*error=3; return(0);}

  /* set clipping rectangle */
  myClipRect.left   = 0;
  myClipRect.right  = width;
  myClipRect.bottom = height;
  myClipRect.top    = 0;

  /* set clipping rectangle */
  portRect = &(((GrafPtr)MAC_WIN(gw))->portRect);
  myClipRect.left   = portRect->left;
  myClipRect.right  = portRect->right  - SCROLL_BAR;
  myClipRect.bottom = portRect->bottom - SCROLL_BAR;
  myClipRect.top    = portRect->top;

  ClipRect(&myClipRect);

  /* fill global and local lower left and upper right in the devices coordinate system */
  gw->Global_LL[0] = Global_LL[0] = h;
  gw->Global_LL[1] = Global_LL[1] = v+dv-SCROLL_BAR;
  gw->Global_UR[0] = Global_UR[0] = h+dh-SCROLL_BAR;
  gw->Global_UR[1] = Global_UR[1] = v+TITLE_BAR;

  gw->Local_LL[0] = Local_LL[0] = 0;
  gw->Local_LL[1] = Local_LL[1] = height;
  gw->Local_UR[0] = Local_UR[0] = width;
  gw->Local_UR[1] = Local_UR[1] = 0;

  TextFont(monaco);

  /* return window ptr */
  *error = 0;
  return((WINDOWID)gw);
}


/****************************************************************************/
/*																			*/
/* Function:  Mac_CloseOutput												*/
/*																			*/
/* Purpose:   close the Mac window associated with theView					*/
/*																			*/
/* Input:	  VIEW *theView: Port of that View will be inited				*/
/*																			*/
/* Output:	  INT: 0 if all was done well									*/
/*				   1 if an error ocurred									*/
/*																			*/
/****************************************************************************/

static INT Mac_CloseOutput (WINDOWID win)
{
  GRAPH_WINDOW *old,*g;

  /* get window */
  old = (GRAPH_WINDOW *)(win);

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
  CloseWindow(MAC_WIN(old));

  /* free memory */
  DisposeMem(guiHeap,old);

  /* no error */
  return(0);
}


/****************************************************************************/
/*																			*/
/* Functions: Mac_ActivateOutput											*/
/*																			*/
/* Purpose:   activate the window win										*/
/*																			*/
/* Input:	  WINDOWID win													*/
/*																			*/
/* Output:	  0 is OK														*/
/*																			*/
/****************************************************************************/

INT Mac_ActivateOutput (WINDOWID win)
{
  WindowPtr theWindow;

  /* set current output window */
  currgw = (GRAPH_WINDOW *)(win);
  theWindow = MAC_WIN(currgw);

  SetPort(theWindow);

  ActivateMacWin (currgw,currgw->currTool);

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:   Mac_UpdateOutput                                                                                         */
/*																			*/
/* Purpose:    Draws all controls and highlights active tool				*/
/*																			*/
/* Input:	   window id													*/
/*																			*/
/* Output:	   0: OK														*/
/*			   1: error, could not complete                                                                 */
/*																			*/
/****************************************************************************/

INT Mac_UpdateOutput (WINDOWID win, INT tool)
{
  GRAPH_WINDOW *gw;
  WindowPtr theWindow;
  Rect myClipRect,*portRect;

  gw = (GRAPH_WINDOW *) win;
  theWindow = MAC_WIN(gw);

  gw->currTool = tool;

  SetPort(theWindow);

  /* identify clipping rect with port */
  portRect = &(((GrafPtr)theWindow)->portRect);
  ClipRect(portRect);

  BeginUpdate(theWindow);
  /* leave this out: uginterface will manage
     EraseRgn(((GrafPtr)theWindow)->visRgn);*/
  EndUpdate(theWindow);

  /* adjust cursor */
  if (theWindow == FrontWindow())
    ActivateMacWin (gw,tool);

  DrawGrowIcon(theWindow);
  DrawToolBox(gw,tool);

  /* reset clipping rectangle */
  myClipRect.left   = portRect->left;
  myClipRect.right  = portRect->right  - SCROLL_BAR;
  myClipRect.bottom = portRect->bottom - SCROLL_BAR;
  myClipRect.top    = portRect->top;

  ClipRect(&myClipRect);

  return (0);
}

/****************************************************************************/
/*
   InitMacOutputDevice - Install output device 'screen'

   SYNOPSIS:
   OUTPUTDEVICE *InitMacOutputDevice (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function installs output device 'screen'.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n       pointer to
   .n       Null if operation failed.
 */
/****************************************************************************/

OUTPUTDEVICE *InitMacOutputDevice (void)
{
  char buffer[32];

  /* create output device */
  if ((MacOutputDevice=CreateOutputDevice("screen"))==NULL) return(NULL);

  /* init output device 'screen' */
  MacOutputDevice->OpenOutput     = Mac_OpenOutput;
  MacOutputDevice->CloseOutput    = Mac_CloseOutput;
  MacOutputDevice->ActivateOutput = Mac_ActivateOutput;
  MacOutputDevice->UpdateOutput   = Mac_UpdateOutput;

  MacOutputDevice->v.locked               = 1;

  InitMacPort ();

  /* get gui heapsize */
  if (GetDefaultValue(DEFAULTSFILENAME,"guimemory",buffer)==0)
    sscanf(buffer," %d ",&guiHeapSize);

  /* allocate gui heap */
  if ((guiHeap=NewHeap(GENERAL_HEAP,guiHeapSize,malloc(guiHeapSize)))==NULL) return(NULL);

  return (MacOutputDevice);
}
