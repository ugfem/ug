// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugView.c                                                      */
/*                                                                          */
/* Purpose:   ug metafile display program for the Macintosh                     */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr1.iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   03.11.92 begin                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#define __MWCW__

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

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
#include <Packages.h>
#include <Files.h>
#include <errors.h>

#include <FileTransfers.h>
#include <FixMath.h>
#include <Devices.h>

#include <string.h>
#include <strings.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "ugView.m"

#ifndef __MWCW__
        #define MAC 1           /* for HDF */
        #include <df.h>
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define HiWrd(aLong) (((aLong) >> 16) & 0xFFFF)
#define LoWrd(aLong) ((aLong) & 0xFFFF)
#define malloc(n) ((void *) NewPtr((Size) n))
#define free(p) DisposePtr((Ptr) p)

#define opNop                           0
#define opMove                          1
#define opDraw                          2
#define opPolyline                      3
#define opPolygon                       4
#define opPolymark                      5
#define opText                          6
#define opCenteredText          7
#define opSetLineWidth          8
#define opSetMarker                     9
#define opSetMarkerSize         10
#define opSetTextSize           11
#define opSetColor                      12
#define opSetEntry                      13
#define opSetPalette            14
#define opNewLine                       15
#define opNewPolyline           16
#define opNewPolygon            17
#define opNewPolymark           18
#define opNewText                       19
#define opNewCenteredText       20

#define SIZE 50
#define CSIZE 256

#define TRFMX(x) l=x; l=l<<16; x = HiWrd(FixMul(l,sx))          /* use only in Rasterize ! */
#define TRFMY(y) l=y; l=l<<16; y = wy-1-HiWrd(FixMul(l,sy))

/* global states */
#define CLOSED 0
#define OPENED 1

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

typedef struct {
  CWindowRecord theWindowRecord;        /* storage for window information		*/
  WindowPtr theWindow;                          /* only WindowPtr although its a CWin.. */
  ControlHandle vScrollBar;                     /* vertical scroll bar                                  */
  ControlHandle hScrollBar;                     /* horizontal scroll bar                                */
  Rect usableRect;                                      /* window without scroll bars			*/
  short red[256];                                       /* color table							*/
  short green[256];
  short blue[256];
} AWindowRecord ;


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static short screenx, screeny;
static Rect dragRect;
static AWindowRecord myWindow;

/* menu handles */
static MenuHandle myMenus[menuCount];           /* menuCount defined in MacGui.m        */

static int globalState=CLOSED;
static char fileName[64];
static short vRefNum;
static short fx,fy;
static short quitFlag=0;

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* Function:  OneButtonBox                                                                      */
/*                                                                          */
/* Purpose:   display text box with one button								*/
/*                                                                          */
/* Input:     char *text: static text to display (<256 characters)          */
/*            char *button1: text for button 1                                          */
/*                                                                          */
/* Output:    int: 1: button 1 has been hit                                     */
/*                                                                          */
/****************************************************************************/

static int OneButtonBox (char *text, char *button1)
{
  DialogPtr theDialog;
  short itemHit,itemType;
  Handle item;
  Rect box;
  char buffer[256];

  theDialog = GetNewDialog(oneButtonDialog,NULL,(WindowPtr) -1);

  /* draw OK Button */
  GetDialogItem(theDialog,OKButton,&itemType,&item,&box);
  SetPort((GrafPtr) theDialog);
  PenSize(3,3);
  InsetRect(&box,-4,-4);
  FrameRoundRect(&box,16,16);

  /* change text of button 1 */
  strcpy(buffer,button1);
  SetControlTitle((ControlHandle) item,c2pstr(buffer));

  GetDialogItem(theDialog,2,&itemType,&item,&box);
  strcpy(buffer,text);
  SetDialogItemText(item,c2pstr(buffer));
  itemHit=0;
  while (itemHit!=OKButton)
  {
    ModalDialog(NULL,&itemHit);
  }
  DisposeDialog(theDialog);
  return(itemHit);
}

/****************************************************************************/
/*                                                                          */
/* Function:  TwoButtonBox                                                                      */
/*                                                                          */
/* Purpose:   display text box with two buttons								*/
/*                                                                          */
/* Input:     char *text: static text to display (<256 characters)          */
/*            char *button1: text for button 1                                          */
/*            char *button2: text for button 2                                          */
/*                                                                          */
/* Output:    int: 1: button 1 has been hit                                     */
/*            int: 2: button 2 has been hit                                     */
/*                                                                          */
/****************************************************************************/

static int TwoButtonBox (char *text, char *button1, char *button2)
{
  DialogPtr theDialog;
  short itemHit,itemType;
  Handle item;
  Rect box;
  char buffer[256];

  theDialog = GetNewDialog(twoButtonDialog,NULL,(WindowPtr) -1);

  /* draw OK Button */
  GetDialogItem(theDialog,OKButton,&itemType,&item,&box);
  SetPort((GrafPtr) theDialog);
  PenSize(3,3);
  InsetRect(&box,-4,-4);
  FrameRoundRect(&box,16,16);

  /* change text of button 1 */
  strcpy(buffer,button1);
  SetControlTitle((ControlHandle) item,c2pstr(buffer));

  /* change text of button 2 */
  GetDialogItem(theDialog,2,&itemType,&item,&box);
  strcpy(buffer,button2);
  SetControlTitle((ControlHandle) item,c2pstr(buffer));

  GetDialogItem(theDialog,3,&itemType,&item,&box);
  strcpy(buffer,text);
  SetDialogItemText(item,c2pstr(buffer));
  itemHit=0;
  while ((itemHit<=0)||(itemHit>=3))
  {
    ModalDialog(NULL,&itemHit);
  }
  DisposeDialog(theDialog);
  return(itemHit);
}

/****************************************************************************/
/*                                                                          */
/* Function:  ThreeButtonBox                                                                    */
/*                                                                          */
/* Purpose:   display text box with three buttons							*/
/*                                                                          */
/* Input:     char *text: static text to display (<256 characters)          */
/*            char *button1: text for button 1                                          */
/*            char *button2: text for button 2                                          */
/*            char *button3: text for button 3                                          */
/*                                                                          */
/* Output:    int: 1: button 1 has been hit                                     */
/*            int: 2: button 2 has been hit                                     */
/*            int: 3: button 3 has been hit                                     */
/*                                                                          */
/****************************************************************************/

static int ThreeButtonBox (char *text, char *button1, char *button2, char *button3)
{
  DialogPtr theDialog;
  short itemHit,itemType;
  Handle item;
  Rect box;
  char buffer[256];

  theDialog = GetNewDialog(threeButtonDialog,NULL,(WindowPtr) -1);

  /* draw OK Button */
  GetDialogItem(theDialog,OKButton,&itemType,&item,&box);
  SetPort((GrafPtr) theDialog);
  PenSize(3,3);
  InsetRect(&box,-4,-4);
  FrameRoundRect(&box,16,16);

  /* change text of button 1 */
  strcpy(buffer,button1);
  SetControlTitle((ControlHandle) item,c2pstr(buffer));

  /* change text of button 2 */
  GetDialogItem(theDialog,2,&itemType,&item,&box);
  strcpy(buffer,button2);
  SetControlTitle((ControlHandle) item,c2pstr(buffer));

  /* change text of button 3 */
  GetDialogItem(theDialog,3,&itemType,&item,&box);
  strcpy(buffer,button3);
  SetControlTitle((ControlHandle) item,c2pstr(buffer));

  GetDialogItem(theDialog,4,&itemType,&item,&box);
  strcpy(buffer,text);
  SetDialogItemText(item,c2pstr(buffer));
  itemHit=0;
  while ((itemHit<=0)||(itemHit>=4))
  {
    ModalDialog(NULL,&itemHit);
  }
  DisposeDialog(theDialog);
  return(itemHit);
}

/****************************************************************************/
/*                                                                          */
/* Function:  @                                                             */
/*                                                                          */
/* Purpose:                                                                 */
/*                                                                          */
/* Input:                                                                   */
/*                                                                          */
/* Output:                                                                  */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

static int CreateApplicationWindow (AWindowRecord *wr, char *fname, short h, short v, short dh, short dv)
{
  Rect r;
  GrafPtr myPort;
  PaletteHandle myPalette;
  char name[80];

  /* init AWindowRecord */
  wr->theWindow = (WindowPtr) wr;

  /* read in resources */
  if (GetNewCWindow(appWinId,(Ptr)wr,(WindowPtr) -1)==NULL)
  {
    return(1);
  }
  myPalette = GetNewPalette(defaultPaletteId);
  SetPalette(wr->theWindow,myPalette,false);

  /* move and size window */
  myPort = (GrafPtr) wr->theWindow;
  SetPort(myPort);
  MoveWindow(wr->theWindow,h,v,false);
  SizeWindow(wr->theWindow,dh+15,dv+15,false);
  strcpy(name,fname);
  SetWTitle(wr->theWindow,c2pstr(name));
  ShowWindow(wr->theWindow);
  SelectWindow(wr->theWindow);
  DrawGrowIcon(wr->theWindow);
  r = myPort->portRect;

  TextFont(kFontIDMonaco);

  /* get the scroll bars */
  wr->vScrollBar = GetNewControl(vScrollBarId,wr->theWindow);
  wr->hScrollBar = GetNewControl(hScrollBarId,wr->theWindow);

  /* set correct size of the scroll bars */
  MoveControl(wr->vScrollBar,r.right-15,-1);
  SizeControl(wr->vScrollBar,16,r.bottom-13);
  SetControlMinimum(wr->vScrollBar,0);
  SetControlMaximum(wr->vScrollBar,0);
  SetControlValue(wr->vScrollBar,0);
  ShowControl(wr->vScrollBar);
  MoveControl(wr->hScrollBar,-1,r.bottom-15);
  SizeControl(wr->hScrollBar,r.right-13,16);
  SetControlMinimum(wr->hScrollBar,0);
  SetControlMaximum(wr->hScrollBar,0);
  SetControlValue(wr->hScrollBar,0);
  ShowControl(wr->hScrollBar);
  DrawControls(wr->theWindow);

  SetRect(&(wr->usableRect),0,0,dh,dv);

  return(0);
}

static void DisposeApplicationWindow (AWindowRecord *wr)
{
  DisposeControl(wr->vScrollBar);
  DisposeControl(wr->hScrollBar);
  CloseWindow(wr->theWindow);
  return;
}


static void Marker (short n, short s, short x, short y)
{
  Rect r;
  PolyHandle myPoly;

  s = s/2;

  r.top = y-s; r.bottom = y+s;
  r.left = x-s; r.right = x+s;

  n = n%11;

  switch (n)
  {
  case 0 :
    FrameRect(&r);
    break;
  case 1 :
    PenPat(&qd.gray);
    PaintRect(&r);
    PenNormal();
    break;
  case 2 :
    PaintRect(&r);
    break;
  case 3 :
    FrameOval(&r);
    break;
  case 4 :
    PenPat(&qd.gray);
    PaintOval(&r);
    PenNormal();
    break;
  case 5 :
    PaintOval(&r);
    break;
  case 6 :
    MoveTo(x,y+s);
    LineTo(x+s,y);
    LineTo(x,y-s);
    LineTo(x-s,y);
    LineTo(x,y+s);
    break;
  case 7 :
    PenPat(&qd.gray);
    myPoly = OpenPoly();
    MoveTo(x,y+s);
    LineTo(x+s,y);
    LineTo(x,y-s);
    LineTo(x-s,y);
    LineTo(x,y+s);
    ClosePoly();
    PaintPoly(myPoly);
    KillPoly(myPoly);
    PenNormal();
    break;
  case 8 :
    myPoly = OpenPoly();
    MoveTo(x,y+s);
    LineTo(x+s,y);
    LineTo(x,y-s);
    LineTo(x-s,y);
    LineTo(x,y+s);
    ClosePoly();
    PaintPoly(myPoly);
    KillPoly(myPoly);
    PenNormal();
    break;
  case 9 :
    MoveTo(x,y+s);
    LineTo(x,y-s);
    MoveTo(x-s,y);
    LineTo(x+s,y-s);
    break;
  case 10 :
    MoveTo(x-s,y+s);
    LineTo(x+s,y-s);
    MoveTo(x-s,y-s);
    LineTo(x+s,y+s);
    break;
  }
}

static int RasterizeFile (FILE *stream, AWindowRecord *myWindow, short wx, short wy)
{
  char *buffer;                                                 /* input buffer						*/
  long blockSize;                                               /* METABUFFERSIZE					*/
  long blockUsed;                                               /* actual buffer size used			*/
  long itemCounter;                                             /* number of commands in buffer		*/
  char *data;                                                           /* data pointer in buffer			*/
  short fx,fy;                                                  /* file screen size					*/
  Fixed sx,sy;                                                  /* scaling factors					*/
  int i,error,j,size;
  char opCode;
  short x,y,r,g,b,n,lw,ts,m,ms,w;
  short x1,y1,x2,y2;
  short xx[SIZE],yy[SIZE];
  PolyHandle myPoly;
  char s[CSIZE];
  unsigned char c;
  RGBColor newColor;
  PaletteHandle myPalette;
  long l;

  /* get file parameters */
  rewind(stream);
  error = fread(&blockSize,4,1,stream); if (error!=1) return(1);       /* block size */
  error = fread(&fx,2,1,stream);            if (error!=1) return(1);       /* x size */
  error = fread(&fy,2,1,stream);        if (error!=1) return(1);       /* y size */

  /* compute scaling factors */
  sx = FixRatio(wx-1,fx-1);
  sy = FixRatio(wy-1,fy-1);

  /* default values */
  lw = 1;
  ts = 12;
  m = 0;
  ms = 6;

  /* allocate input buffer */
  buffer = malloc(blockSize);
  if (buffer==NULL) return(1);

  SetPort((GrafPtr)(myWindow->theWindow));
  EraseRect(&(myWindow->usableRect));

  /* loop through the blocks */
  while (!feof(stream))
  {
    /* read block parameters */
    error = fread(&blockUsed,4,1,stream);    if (error!=1) {free(buffer); return(1);}
    error = fread(&itemCounter,4,1,stream);  if (error!=1) {free(buffer); return(1);}
    error = fread(buffer,blockUsed,1,stream);if (error!=1) {free(buffer); return(1);}

    /* init pointer to next item */
    data = buffer;

    /* for all items */
    for (i=0; i<itemCounter; i++)
    {
      /* get op code */
      opCode = *(data++);

      switch (opCode)
      {
      case opMove :
        x = *((short *)data);
        data += 2;
        y = *((short *)data);
        data += 2;
        TRFMX(x);
        TRFMY(y);
        MoveTo(x,y);
        break;

      case opDraw :
        x = *((short *)data);
        data += 2;
        y = *((short *)data);
        data += 2;
        TRFMX(x);
        TRFMY(y);
        LineTo(x,y);
        break;

      case opPolyline :
        n = *((short *)data);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        for (j=0; j<n; j++)
        {
          TRFMX(xx[j]);
          TRFMY(yy[j]);
        }
        MoveTo(xx[0],yy[0]);
        for (j=1; j<n; j++) LineTo(xx[j],yy[j]);
        break;

      case opPolygon :
        n = *((short *)data);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        for (j=0; j<n; j++)
        {
          TRFMX(xx[j]);
          TRFMY(yy[j]);
        }
        if (n<3) break;
        myPoly = OpenPoly();
        MoveTo(xx[0],yy[0]);
        for (j=1; j<n; j++) LineTo(xx[j],yy[j]);
        LineTo(xx[0],yy[0]);
        ClosePoly();
        PaintPoly(myPoly);
        FramePoly(myPoly);
        KillPoly(myPoly);
        break;

      case opPolymark :
        n = *((short *)data);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        for (j=0; j<n; j++)
        {
          TRFMX(xx[j]);
          TRFMY(yy[j]);
        }
        for (j=0; j<n; j++) Marker(m,ms,xx[j],yy[j]);
        break;

      case opText :
        n = *((short *)data);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        DrawString((ConstStr255Param)c2pstr(s));
        break;

      case opCenteredText :
        x = *((short *)data);
        data += 2;
        y = *((short *)data);
        data += 2;
        TRFMX(x);
        TRFMY(y);
        n = *((short *)data);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        c2pstr(s);
        w = StringWidth((ConstStr255Param)s);
        MoveTo(x-w/2,y+ts/2);
        DrawString((ConstStr255Param)s);
        break;

      case opSetLineWidth :
        n = *((short *)data);
        data += 2;
        lw = n;
        PenSize(n,n);
        break;

      case opSetTextSize :
        n = *((short *)data);
        data += 2;
        ts = n;
        TextSize(n);
        break;

      case opSetMarker :
        n = *((short *)data);
        data += 2;
        m = n;
        break;

      case opSetMarkerSize :
        n = *((short *)data);
        data += 2;
        ms = n;
        break;

      case opSetColor :
        c = *((unsigned char *)data);
        data++;
        PmForeColor((short)c);
        break;

      case opSetEntry :
        c = *((unsigned char *)data);
        data++;
        r = (short) (*((unsigned char *)data));
        data++;
        g = (short) (*((unsigned char *)data));
        data++;
        b = (short) (*((unsigned char *)data));
        data++;
        myPalette = GetPalette(myWindow->theWindow);
        myWindow->red[c]   = newColor.red   = r<<8;
        myWindow->green[c] = newColor.green = g<<8;
        myWindow->blue[c]  = newColor.blue  = b<<8;
        SetEntryColor(myPalette,(short) c,&newColor);
        ActivatePalette(myWindow->theWindow);
        break;

      case opSetPalette :
        x = (short) (*((unsigned char *)data));
        data++;
        y = (short) (*((unsigned char *)data));
        data++;
        myPalette = GetPalette(myWindow->theWindow);
        for (j=x; j<=y; j++)
        {
          r = (short) (*((unsigned char *)data));
          data++;
          g = (short) (*((unsigned char *)data));
          data++;
          b = (short) (*((unsigned char *)data));
          data++;
          myWindow->red[j]   = newColor.red   = r<<8;
          myWindow->green[j] = newColor.green = g<<8;
          myWindow->blue[j]  = newColor.blue  = b<<8;
          SetEntryColor(myPalette,(short) j,&newColor);
        }
        ActivatePalette(myWindow->theWindow);
        break;

      case opNewLine :
        lw = *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        x1 = *((short *)data);
        data += 2;
        y1 = *((short *)data);
        data += 2;
        x2 = *((short *)data);
        data += 2;
        y2 = *((short *)data);
        data += 2;
        TRFMX(x1);
        TRFMY(y1);
        TRFMX(x2);
        TRFMY(y2);
        PenSize(lw,lw);
        PmForeColor((short)c);
        MoveTo(x1,y1);
        LineTo(x2,y2);
        break;

      case opNewPolyline :
        lw = *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        n = *((short *)data);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        for (j=0; j<n; j++)
        {
          TRFMX(xx[j]);
          TRFMY(yy[j]);
        }
        PenSize(lw,lw);
        PmForeColor((short)c);
        MoveTo(xx[0],yy[0]);
        for (j=1; j<n; j++) LineTo(xx[j],yy[j]);
        break;

      case opNewPolygon :
        c = *((unsigned char *)data);
        data++;
        n = *((short *)data);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        for (j=0; j<n; j++)
        {
          TRFMX(xx[j]);
          TRFMY(yy[j]);
        }
        if (n<3) break;
        PmForeColor((short)c);
        myPoly = OpenPoly();
        MoveTo(xx[0],yy[0]);
        for (j=1; j<n; j++) LineTo(xx[j],yy[j]);
        LineTo(xx[0],yy[0]);
        ClosePoly();
        PaintPoly(myPoly);
        FramePoly(myPoly);
        KillPoly(myPoly);
        break;

      case opNewPolymark :
        m = *((unsigned char *)data);
        data++;
        ms = *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        n = *((short *)data);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        for (j=0; j<n; j++)
        {
          TRFMX(xx[j]);
          TRFMY(yy[j]);
        }
        PmForeColor((short)c);
        for (j=0; j<n; j++) Marker(m,ms,xx[j],yy[j]);
        break;

      case opNewText :
        ts = *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        x = *((short *)data);
        data += 2;
        y = *((short *)data);
        data += 2;
        TRFMX(x);
        TRFMY(y);
        n = *((short *)data);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        MoveTo(x,y);
        TextSize(ts);
        PmForeColor((short)c);
        DrawString((ConstStr255Param)c2pstr(s));
        break;

      case opNewCenteredText :
        ts = *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        x = *((short *)data);
        data += 2;
        y = *((short *)data);
        data += 2;
        TRFMX(x);
        TRFMY(y);
        n = *((short *)data);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        c2pstr(s);
        w = StringWidth((ConstStr255Param)s);
        TextSize(ts);
        PmForeColor((short)c);
        MoveTo(x-w/2,y+ts/2);
        DrawString((ConstStr255Param)c2pstr(s));
        break;

      default :
        break;
      }
    }
  }

  return(0);
}

static int GetFileScreen (FILE *stream, short *fx, short *fy)
{
  long blockSize;                                               /* METABUFFERSIZE					*/
  int error;

  /* get file parameters */
  rewind(stream);
  error = fread(&blockSize,4,1,stream); if (error!=1) return(1);       /* block size */
  error = fread(fx,2,1,stream);             if (error!=1) return(1);       /* x size */
  error = fread(fy,2,1,stream);        if (error!=1) return(1);       /* y size */
  return(0);
}



static void SetUpMenus ()
{
  int i;

  for (i=0; i<menuCount; i++)
    myMenus[i] = GetMenu(i+appleID);
  AppendResMenu(myMenus[0],'DRVR');
  for (i=0; i<menuCount; i++)
    InsertMenu(myMenus[i],0);                                           /* Add menus to menu bar */
  DrawMenuBar();
}

static void SetState (int i)
{
  globalState = i;
  switch (i)
  {
  case CLOSED :
    EnableItem(myMenus[fileM],openCommand);
    DisableItem(myMenus[fileM],closeCommand);
    EnableItem(myMenus[fileM],quitCommand);
    DisableItem(myMenus[displayM],refreshCommand);
    DisableItem(myMenus[displayM],hdfCommand);
    DisableItem(myMenus[displayM],pictCommand);
    break;

  case OPENED :
    DisableItem(myMenus[fileM],openCommand);
    EnableItem(myMenus[fileM],closeCommand);
    EnableItem(myMenus[fileM],quitCommand);
    EnableItem(myMenus[displayM],refreshCommand);
    EnableItem(myMenus[displayM],hdfCommand);
    EnableItem(myMenus[displayM],pictCommand);
    break;
  }
}

static void AboutCommand (void)
{}

static void RefreshCommand (void)
{
  int error;
  FILE *stream;

  if (globalState!=OPENED) return;

  /* set volume */
  error = SetVol(NULL,vRefNum);
  if (((int)error)!=0)
  {
    OneButtonBox("Could not set volume","OK");
    return;
  }

  stream = fopen(fileName,"rb");
  if (stream==NULL)
  {
    OneButtonBox("Could not open file","OK");
    return;
  }

  /* draw something */
  RasterizeFile(stream,&myWindow,fx,fy);
  fclose(stream);

  DrawGrowIcon(myWindow.theWindow);
  DrawControls(myWindow.theWindow);

  return;
}

static void OpenCommand (void)
{
  int error;
  FILE *stream;
  SFReply reply;
  Point where;
  char prompt[32];
  SFTypeList typeList;

  if (globalState==OPENED) return;

  /* get file name */
  strcpy(prompt,"choose metafile:");
  where.h = 20; where.v = 40;
  SFGetFile(where,c2pstr(prompt),NULL,-1,typeList,NULL,&reply);
  if (!reply.good)
  {
    return;
  }

  /* save name and path */
  vRefNum = reply.vRefNum;
  strcpy(fileName,p2cstr(reply.fName));

  /* set volume */
  error = SetVol(NULL,vRefNum);
  if (((int)error)!=0)
  {
    OneButtonBox("Could not set volume","OK");
    return;
  }

  stream = fopen(fileName,"rb");
  if (stream==NULL)
  {
    OneButtonBox("Could not open file","OK");
    return;
  }

  GetFileScreen(stream,&fx,&fy);

  /* open a window */
  error = CreateApplicationWindow(&myWindow,fileName,20,40,fx,fy);
  fclose(stream);
  SetState(OPENED);
  RefreshCommand();

  return;
}

static void CloseCommand (void)
{
  if (globalState==CLOSED) return;

  DisposeApplicationWindow(&myWindow);
  SetState(CLOSED);
  return;
}

static void QuitCommand (void)
{
  CloseCommand();
  quitFlag = 1;
}

#ifndef __MWCW__

/****************************************************************************/
/*                                                                          */
/* Function:  HDFCommand                                                                        */
/*                                                                          */
/* Purpose:   save current view to a hdf image file							*/
/*            currently only (8 bit mode)                                                       */
/*                                                                          */
/* Input:     none                                                                                      */
/*                                                                          */
/* Output:    none                                                          */
/*                                                                          */
/****************************************************************************/

static void SaveToHDF (char *name)
{
  CGrafPtr theCGrafPort;
  WindowPtr theWindow;
  PixMapPtr thePixMap;
  char buffer[80];
  GDHandle gd;
  GDPtr theGDevice;
  unsigned char *theImage,*srcPtr,*dstPtr,*baseAddr;
  Rect screen,window;
  unsigned long size;
  int rowBytes,row,col,i,j;
  PaletteHandle thePalette;
  unsigned char transTable[256],hdfPalette[768];
  RGBColor theColor;

  theWindow = myWindow.theWindow;
  theCGrafPort = (CGrafPtr) theWindow;
  gd = (GDHandle) GetGDevice();
  theGDevice = *gd;

  /* check frontwindow */
  if (theWindow!=FrontWindow())
  {
    OneButtonBox("window to save must be the front window","OK");
    return;
  }

  /* check screen depth and mode */
  thePixMap = *(theGDevice->gdPMap);
  rowBytes = thePixMap->rowBytes&0x1FFF;
  baseAddr =  (unsigned char *) thePixMap->baseAddr;
  if ((thePixMap->pixelType!=0)||(thePixMap->pixelSize!=8))
  {
    OneButtonBox("Sorry, 8 bits/pixel and chunky mode required","OK");
    return;
  }

  /* compute window contents coordinates */
  thePixMap = *(theGDevice->gdPMap);
  SetRect(&screen,thePixMap->bounds.left,thePixMap->bounds.top,thePixMap->bounds.right,thePixMap->bounds.bottom);
  thePixMap = *(theCGrafPort->portPixMap);
  SetRect(&window,-thePixMap->bounds.left,-thePixMap->bounds.top,theCGrafPort->portRect.right-16-thePixMap->bounds.left,theCGrafPort->portRect.bottom-16-thePixMap->bounds.top);
  SectRect(&window,&screen,&window);

  /* allocate image buffer */
  size = ((unsigned long)(window.right-window.left))*((unsigned long)(window.bottom-window.top));
  theImage = (unsigned char *) malloc(size);
  if (theImage==NULL)
  {
    OneButtonBox("Sorry, not enough memory for image buffer","OK");
    return;
  }

  /* create translation table and palette */
  thePalette = GetPalette(theWindow);
  for (i=0; i<256; i++)
  {
    GetEntryColor(thePalette,i,&theColor);
    hdfPalette[i*3] = (theColor.red)>>8;
    hdfPalette[i*3+1] = (theColor.green)>>8;
    hdfPalette[i*3+2] = (theColor.blue)>>8;
    j = Color2Index(&theColor);
    transTable[j] = (unsigned char) i;
  }

  /* copy image */
  dstPtr = theImage;
  for (row=window.top; row<window.bottom; row++)
  {
    srcPtr = baseAddr+row*rowBytes+window.left;
    for (col=window.left; col<window.right; col++)
      *(dstPtr++) = transTable[*(srcPtr++)];
  }

  /* set palette for subsequent image */
  if (DFR8setpalette(hdfPalette))
    OneButtonBox("setpalette failed","OK");

  /* write image to disk */
  if(DFR8putimage(name,theImage,window.right-window.left,window.bottom-window.top,DFTAG_RLE))
  {
    free(theImage);
    OneButtonBox("Image not written correctly","OK");
    return;
  }

  /* free memory */
  free(theImage);

  return;
}

static void HDFCommand (void)
{
  SFReply reply;
  Point where;
  char prompt[128];
  OSErr error;
  char buffer[80];
  char *s;

  if (globalState!=OPENED) return;

  /* get file name */
  strcpy(buffer,"Save HDF image as:");
  strcpy(prompt,fileName);
  if ((s=strstr(prompt,".meta"))!=NULL)
    *s = '\0';
  strcat(prompt,".hdf");
  where.h = 20; where.v = 40;
  SFPutFile(where,c2pstr(buffer),c2pstr(prompt),NULL,&reply);
  if (!reply.good) return;

  /* set volume */
  error = SetVol(NULL,reply.vRefNum);
  if (((int)error)!=0)
  {
    OneButtonBox("Could not set volume","OK");
    return;
  }
  /* change pascal string to C string */
  p2cstr(reply.fName);

  RefreshCommand();
  SaveToHDF(reply.fName);
  return;
}
#endif          /* ifndef __MWCW__ */

/****************************************************************************/
/*                                                                          */
/* Function:  PICTCommand                                                                       */
/*                                                                          */
/* Purpose:   save current view in a PICT image file						*/
/*                                                                          */
/* Input:     none                                                                                      */
/*                                                                          */
/* Output:    none                                                          */
/*                                                                          */
/****************************************************************************/

static int PicCounter ;
static short FileSys_RefNR ;
static PicHandle MyPicHand ;

static pascal void MyPutPicData ( void            *Data_Ptr ,
                                  short ByteCount )
{
  long LongCount ;
  OSErr error;

  LongCount  = ByteCount ;
  PicCounter = PicCounter + ByteCount ;
  error          = FSWrite ( FileSys_RefNR, &LongCount, Data_Ptr ) ;
  if (MyPicHand != NULL)
    (**MyPicHand).picSize = PicCounter ;
}  /* MyPutPicData */

#define cHRes           0x00480000
#define cVRes           0x00480000

static void SaveToPICT (SFReply *reply)
{
  Rect MyPicFrame ;
  CQDProcs MyPicProcs,*SavePtr ;
  OSErr error ;
  long LongZero, LongCount ;
  short Counter ;
  CGrafPtr theCGrafPort;
  WindowPtr theWindow;
  OpenCPicParams myOpenCPicParams;

  /* get my window */
  theWindow = myWindow.theWindow;
  theCGrafPort = (CGrafPtr) theWindow;

  error = FSOpen(reply->fName, reply->vRefNum, &FileSys_RefNR ) ;
  SetStdCProcs ( &MyPicProcs ) ;
  SavePtr = (CQDProcs *) theWindow->grafProcs;
  theWindow->grafProcs   = (QDProcs *) (&MyPicProcs) ;
  MyPicProcs.putPicProc  = NewQDPutPicProc(MyPutPicData) ;

  LongZero   = 0 ;
  LongCount  = 4 ;

  PicCounter = sizeof(Picture) ;
  for (Counter=1 ; Counter <= 128+sizeof(Picture) ; Counter++ )
    error = FSWrite(FileSys_RefNR, &LongCount, &LongZero );
  error = SetFPos(FileSys_RefNR, fsFromStart, 512+sizeof(Picture)) ;

  SetPort((GrafPtr)(theWindow));
  MyPicFrame = ((GrafPtr)(theWindow))->portRect;
  MyPicFrame.right = MyPicFrame.right-16;
  MyPicFrame.bottom = MyPicFrame.bottom-16;
  MyPicHand  = NULL ;
  myOpenCPicParams.srcRect        = MyPicFrame;
  myOpenCPicParams.hRes           = cHRes;
  myOpenCPicParams.vRes           = cVRes;
  myOpenCPicParams.version        = -2;
  myOpenCPicParams.reserved1      = 0;
  myOpenCPicParams.reserved2      = 0;
  MyPicHand  = OpenCPicture ( &myOpenCPicParams );
  ClipRect(&MyPicFrame);

  /* draw picture */
  RefreshCommand();

  ClosePicture() ;

  /* close File */
  error = SetFPos(FileSys_RefNR, fsFromStart, 512) ;
  LongCount = sizeof(Picture) ;
  error = FSWrite(FileSys_RefNR, &LongCount, (void *) (*MyPicHand) );
  error = FSClose(FileSys_RefNR) ;

  KillPicture(MyPicHand);
  theWindow->grafProcs   = (QDProcs *) SavePtr ;

  return;
}


static void PICTCommand (void)
{
  OSErr error ;
  SFReply reply ;
  Point wher ;
  char buffer1[32],buffer2[128];
  char *s;

  if (globalState!=OPENED) return;

  /* get file name */
  wher.v = 40 ;
  wher.h = 40 ;
  strcpy(buffer1,"Name of Pict-File :");
  strcpy(buffer2,fileName);
  if ((s=strstr(buffer2,".meta"))!=NULL)
    *s = '\0';
  strcat(buffer2,".pict");
  SFPutFile ( wher, c2pstr(buffer1), c2pstr(buffer2), NULL, &reply ) ;
  if (!reply.good) return;

  error = Create ( reply.fName, reply.vRefNum, *((OSType *)"DAD2"), *((OSType *)"PICT") ) ;
  if ((error!=noErr) && (error!=dupFNErr))
  {
    OneButtonBox("Could not create file","OK");
    return;
  }

  SaveToPICT(&reply);
  return;
}


/****************************************************************************/
/*                                                                          */
/* Function:  ScheduleCommand                                                           */
/*                                                                          */
/* Purpose:   pass control to the command selected from a menu				*/
/*                                                                          */
/* Input:     short theMenu: menu number                                        */
/*            short theItem: item number                                                */
/*                                                                          */
/* Output:    none                                                          */
/*                                                                          */
/****************************************************************************/

static void ScheduleCommand (short theMenu,short theItem)
{
  switch (theMenu)
  {
  case appleID :
    switch (theItem)
    {
    case aboutCommand :              AboutCommand(); break;
    }
    break;

  case fileID :
    switch (theItem)
    {
    case openCommand :               OpenCommand(); break;
    case closeCommand :              CloseCommand(); break;
    case quitCommand :               QuitCommand(); break;
    }
    break;

  case displayID :
    switch (theItem)
    {
    case refreshCommand :     RefreshCommand(); break;
                                #ifndef __MWCW__
    case hdfCommand :         HDFCommand(); break;
                                #endif
    case pictCommand :        PICTCommand(); break;
    }
    break;

  default :
    break;
  }
  return;
}

/****************************************************************************/
/*                                                                          */
/* Function:  DoCommand                                                         */
/*                                                                          */
/* Purpose:   process command selection										*/
/*                                                                          */
/* Input:     long: menu and item information                                                   */
/*                                                                          */
/* Output:    none                                                                              */
/*                                                                          */
/****************************************************************************/

static void DoCommand (long mResult)
{
  short temp,theItem,theMenu;
  Str255 name;

  theItem = LoWrd(mResult);
  theMenu = HiWrd(mResult);

  switch (theMenu)
  {

  case appleID :
    if (theItem==1)
    {
      ScheduleCommand(theMenu,theItem);
    }
    else
    {
      GetMenuItemText(myMenus[appleM],theItem,name);
      temp = OpenDeskAcc(name);
    }
    break;

  default :
    ScheduleCommand(theMenu,theItem);
    break;
  }
  HiliteMenu(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  DoKey                                                             */
/*                                                                          */
/* Purpose:   process keyboard events                                                                           */
/*                                                                          */
/* Input:     EventRecord *theEvent: pointer to keyboard event              */
/*                                                                          */
/* Output:    none                                                                              */
/*                                                                          */
/****************************************************************************/

static void DoKey (EventRecord *theEvent)
{
  long theChar;

  theChar = BitAnd(theEvent->message,charCodeMask);
  if (BitAnd(theEvent->modifiers,cmdKey)!=0)
    DoCommand(MenuKey(theChar));
}


/****************************************************************************/
/*                                                                          */
/* Function:  ProcessEvent                                                  */
/*                                                                          */
/* Purpose:   get next event from the queue and call the appropriate proc.  */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    none                                                          */
/*                                                                          */
/****************************************************************************/

static void ProcessEvent ()
{
  EventRecord myEvent;
  WindowPtr whichWindow;

  /* do system tasks */
  SystemTask();

  /* get event */
  if (GetNextEvent(everyEvent,&myEvent))
  {
    switch (myEvent.what)
    {
    case mouseDown :
      switch (FindWindow(myEvent.where,&whichWindow))
      {
      case inSysWindow :
        SystemClick(&myEvent,whichWindow);
        break;
      case inMenuBar :
        DoCommand(MenuSelect(myEvent.where));
        break;
      case inDrag :
        DragWindow(whichWindow,myEvent.where,&dragRect);
        break;
      case inContent :
        if (whichWindow!=FrontWindow())
          SelectWindow(whichWindow);
        else
        {
          /*DoContentClick(whichWindow,&myEvent);*/
        }
        break;
      case inGrow :
        /*DoGrowWindow(whichWindow,&myEvent);*/
        break;
      case inGoAway :
        /*DoGoAway(whichWindow,&myEvent);*/
        break;
      }
      break;                           /* mouseDown */

    case keyDown :
    case autoKey :
      DoKey(&myEvent);
      break;

    case activateEvt :
      /*DoActivate(&myEvent);*/
      break;

    case updateEvt :
      whichWindow = (WindowPtr) myEvent.message;
      /*DoUpdate(whichWindow,&myEvent);*/
      break;
    }
  }
} /* ProcessEvent */



int main (void)
{
  /* init the Macintosh toolbox */
  InitGraf(&qd.thePort);
  InitFonts();
  FlushEvents(everyEvent,0);
  InitWindows();
  InitMenus();
  TEInit();
  InitDialogs(nil);
  InitCursor();

  /* read the screen's size */
  screenx = ((qd.screenBits).bounds).right;
  screeny = ((qd.screenBits).bounds).bottom;
  SetRect(&dragRect,4,24,((qd.screenBits).bounds).right-4,((qd.screenBits).bounds).bottom-4);


  SetUpMenus ();
  SetState(CLOSED);
  OpenCommand();

  /* handle events till the end */
  while (!quitFlag)
  {
    ProcessEvent();
  }

  return(0);
}
