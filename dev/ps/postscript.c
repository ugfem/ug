// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  postscript.c													*/
/*																			*/
/* Purpose:   write graphics postscript files                                                           */
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   07.05.96 begin, ug version 3.2								*/
/*																			*/
/* Remarks:   used most of the code of m2ps by P. Bastian					*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "defaults.h"
#include "fileopen.h"
#include "devices.h"
#include "initdev.h"
#include "general.h"


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define TEXTSIZE                        8

#define CREATOR "ug postscript output"


/* postscript defaults */
#define  PS_DEFAULT_LINE_WIDTH     1
#define  PS_DEFAULT_FONT           "Monaco"
#define  PS_FONT_FACTOR            10
#define  LW_FACTOR                 0.03
#define  LW_SCALE                  50.0
#define  PAPER_X                   600
#define  PAPER_Y                   850

#define REL_CHARWIDTH                      0.7

#define COLORS          256
#define GRAY                    0.5
#define GRAY_CC                 (~(unsigned long)0)

#define TRFMX(pt) (((float)(pt.x))*mxx + ((float)(pt.y))*mxy + tx)
#define TRFMY(pt) (((float)(pt.x))*myx + ((float)(pt.y))*myy + ty)

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

typedef struct {

  FILE *psfile;                                         /* output file							*/

  short landscape;                                      /* landscape format						*/

  /* transformation */
  float tx,ty;                                          /* translation							*/
  float mxx,mxy,myx,myy;                        /* matrix								*/

} PSWINDOW ;

static OUTPUTDEVICE *PSOutputDevice=NULL; /* outputdevice that has been ini*/
static PSWINDOW *currPSW=NULL;          /* current postscript window			*/
static FILE *currPSF=NULL;                      /* current postscript file				*/
static short PSmarker;
static short PSmarkersize;

static SHORT_POINT PScp={0,0};
static short PSlw=-1, PSts=-1;
static short PScc=0;
static float tx,ty,mxx,myy,mxy,myx;
static short landscape = 0;

/* static color table for all postscript files */
static float red[COLORS];
static float green[COLORS];
static float blue[COLORS];

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static void PSMoveTo (SHORT_POINT point)
{
  PScp.x = point.x; PScp.y = point.y;
  return;
}

static void PSDrawTo (SHORT_POINT point)
{
  fprintf(currPSF,"%g %g M %g %g S\n",
          TRFMX(PScp),TRFMY(PScp),TRFMX(point),TRFMY(point));
  PScp.x = point.x; PScp.y = point.y;
  return;
}

static void PSCircle (SHORT_POINT point, short r)
{
  SHORT_POINT aux;
  short xr,yr;

  aux.x = 0;
  aux.y = r;
  xr = TRFMX(aux);
  yr = TRFMY(aux);
  r  = sqrt(xr*xr+yr*yr);
  fprintf(currPSF,"N\n");
  fprintf(currPSF,"%g %g M\n",TRFMX(point)+r,TRFMY(point));
  fprintf(currPSF,"%g %g %g %g %g arc\n",TRFMX(point),TRFMY(point),(float)r,(float)0,(float)360);

  fprintf(currPSF,"stroke\n");

  return;
}

static void PSFilledCircle (SHORT_POINT point, short r)
{
  SHORT_POINT aux;
  short xr,yr;

  aux.x = 0;
  aux.y = r;
  xr = TRFMX(aux);
  yr = TRFMY(aux);
  r  = sqrt(xr*xr+yr*yr);
  fprintf(currPSF,"N\n");
  fprintf(currPSF,"%g %g M\n",TRFMX(point)+r,TRFMY(point));
  fprintf(currPSF,"%g %g %g %g %g arc\n",TRFMX(point),TRFMY(point),(float)r,(float)0,(float)360);

  fprintf(currPSF,"C\n");

  return;
}

static void PSPolyline (SHORT_POINT *points, INT nb)
{
  int i;

  fprintf(currPSF,"N\n");
  fprintf(currPSF,"%g %g M\n",TRFMX(points[0]),TRFMY(points[0]));
  for (i=1; i<nb; i++)
    fprintf(currPSF,"%g %g L\n",TRFMX(points[i]),TRFMY(points[i]));
  fprintf(currPSF,"stroke\n");

  return;
}

static void PSPolygon (SHORT_POINT *points, INT nb)
{
  int i;

  fprintf(currPSF,"N\n");
  fprintf(currPSF,"%g %g M\n",TRFMX(points[0]),TRFMY(points[0]));
  for (i=1; i<nb; i++)
    fprintf(currPSF,"%g %g L\n",TRFMX(points[i]),TRFMY(points[i]));
  fprintf(currPSF,"C\n");

  return;
}

static void PSInversePolygon (SHORT_POINT *points, INT nb)
{
  return;
}

static void PSErasePolygon (SHORT_POINT *points, INT nb)
{
  return;
}

void PSPrintColor (float color)
{
  if (color==0.0) fprintf(currPSF,"%d ",0);
  else if (color==1.0) fprintf(currPSF,"%d ",1);
  else fprintf(currPSF,"%.3f ",color);
  return;
}

static void PSGrayForeground (void)
{
  if (PScc==GRAY_CC) return;

  fprintf(currPSF,"%.1f %.1f %.1f R\n",GRAY,GRAY,GRAY);
  PScc = GRAY_CC;
  return;
}

static void PSForeground (long index)
{
  if (PScc==index) return;

  PSPrintColor((float)red[index]);
  PSPrintColor((float)green[index]);
  PSPrintColor((float)blue[index]);
  fprintf(currPSF,"R\n");
  PScc = index;
  return;
}

static void PSMarker (short n, short s, SHORT_POINT point)
{
  short top, left, bottom, right;
  SHORT_POINT points[10],pt;
  unsigned long cc;

  s = s/2;

  top = point.y+s; bottom = point.y-s;
  left = point.x-s; right = point.x+s;

  n = n%NMARKERS;

  switch (n)
  {
  case EMPTY_SQUARE_MARKER :
    points[0].x = left;  points[0].y = bottom;
    points[1].x = right; points[1].y = bottom;
    points[2].x = right; points[2].y = top;
    points[3].x = left;  points[3].y = top;
    points[4].x = left;  points[4].y = bottom;
    PSPolyline(points,5);
    break;
  case GRAY_SQUARE_MARKER :
    points[0].x = left;  points[0].y = bottom;
    points[1].x = right; points[1].y = bottom;
    points[2].x = right; points[2].y = top;
    points[3].x = left;  points[3].y = top;
    cc = PScc;
    PSGrayForeground();
    PSPolygon(points,4);
    PSForeground(cc);
    break;
  case FILLED_SQUARE_MARKER :
    points[0].x = left;  points[0].y = bottom;
    points[1].x = right; points[1].y = bottom;
    points[2].x = right; points[2].y = top;
    points[3].x = left;  points[3].y = top;
    PSPolygon(points,4);
    break;
  case EMPTY_CIRCLE_MARKER :
    PSCircle(point,s);
    break;
  case GRAY_CIRCLE_MARKER :
    PSGrayForeground();
    PSFilledCircle(point,s);
    cc = PScc;
    PSForeground(cc);
    break;
  case FILLED_CIRCLE_MARKER :
    PSFilledCircle(point,s);
    break;
  case EMPTY_RHOMBUS_MARKER :
    points[0].x = point.x;      points[0].y = bottom;
    points[1].x = right;        points[1].y = point.y;
    points[2].x = point.x;      points[2].y = top;
    points[3].x = left;         points[3].y = point.y;
    points[4].x = point.x;      points[4].y = bottom;
    PSPolyline(points,5);
    break;
  case GRAY_RHOMBUS_MARKER :
    points[0].x = point.x;      points[0].y = bottom;
    points[1].x = right;        points[1].y = point.y;
    points[2].x = point.x;      points[2].y = top;
    points[3].x = left;         points[3].y = point.y;
    cc = PScc;
    PSGrayForeground();
    PSPolygon(points,4);
    PSForeground(cc);
    break;
  case FILLED_RHOMBUS_MARKER :
    points[0].x = point.x;      points[0].y = bottom;
    points[1].x = right;        points[1].y = point.y;
    points[2].x = point.x;      points[2].y = top;
    points[3].x = left;         points[3].y = point.y;
    PSPolygon(points,4);
    break;
  case PLUS_MARKER :
    pt.x = point.x; pt.y = bottom;
    PSMoveTo(pt);
    pt.x = point.x; pt.y = top;
    PSDrawTo(pt);
    pt.x = right; pt.y = point.y;
    PSMoveTo(pt);
    pt.x = left; pt.y = point.y;
    PSDrawTo(pt);
    break;
  case CROSS_MARKER :
    pt.x = left; pt.y = bottom;
    PSMoveTo(pt);
    pt.x = right; pt.y = top;
    PSDrawTo(pt);
    pt.x = right; pt.y = bottom;
    PSMoveTo(pt);
    pt.x = left; pt.y = top;
    PSDrawTo(pt);
    break;
  }
}

static void PSPolymark (short n, SHORT_POINT *points)
{
  int i;

  for (i=0; i<n; i++)
    PSMarker(PSmarker,PSmarkersize,points[i]);

  return;
}

static void PSInvPolymark (short n, SHORT_POINT *points)
{}

static void PrintPSString (const char *s)
{
  fputc('(',currPSF);

  while (*s != '\0')
  {
    switch (*s)
    {
    case '\\' :
    case '(' :
    case ')' : fputc('\\',currPSF);
    default  : fputc(*s,currPSF);
    }
    s++;
  }
  fputc(')',currPSF);
}

static void PSText (const char *s, INT mode)
{
  fprintf(currPSF,"%g %g M\n",TRFMX(PScp),TRFMY(PScp));
  if (landscape) fprintf(currPSF,"90 rotate\n");
  PrintPSString(s);
  fprintf(currPSF," show N\n");
  if (landscape) fprintf(currPSF,"-90 rotate\n");
  return;
}

static void PSCenteredText (SHORT_POINT point, const char *s, INT mode)
{
  /* centering?? */
  point.x -= 0.5*REL_CHARWIDTH*PSts*strlen(s);
  PSMoveTo(point);
  PSText(s,mode);
}

static void PSClearViewPort (void)
{
  return;
}

static void PSSetLineWidth (short n)
{
  if (n<1) n=1;
  if (n==PSlw) return;

  fprintf(currPSF,"%.3f W\n",LW_FACTOR+((float)(n-1))*LW_SCALE*LW_FACTOR);
  PSlw = n;
  return;
}

static void PSSetTextSize (short n)
{
  if (n==PSts) return;

  fprintf(currPSF,"/%s findfont %d scalefont setfont\n",PS_DEFAULT_FONT,(int)n);
  PSts = n;
  return;
}

static void PSSetMarker (short n)
{
  PSmarker = n;
  return;
}

static void PSSetMarkerSize (short n)
{
  PSmarkersize = n;
  return;
}

static void PSSetPaletteEntry (long index, short r, short g, short b)
{
  red[index]   = ((float)r)/255.0;
  green[index] = ((float)g)/255.0;
  blue[index]  = ((float)b)/255.0;
  PSPrintColor((float)red[index]);
  PSPrintColor((float)green[index]);
  PSPrintColor((float)blue[index]);
  fprintf(currPSF,"R\n");
  PScc = index;
  return;
}

static void PSSetPalette (long x, long count, short *r, short *g, short *b)
{
  int i,y;

  y = x+count;
  for (i=x; i<y; i++)
  {
    red[i]   = ((float)r[i-x])/255.0;
    green[i] = ((float)g[i-x])/255.0;
    blue[i]  = ((float)b[i-x])/255.0;
  }
  PSPrintColor((float)red[x]);
  PSPrintColor((float)green[x]);
  PSPrintColor((float)blue[x]);
  fprintf(currPSF,"R\n");
  PScc = (unsigned char) x;
  return;
}

static void PSGetPaletteEntry (long index, short *r, short *g, short *b)
{
  return;
}

static void PSFlush (void)
{
  return;
}

/****************************************************************************/
/*
   InitPSPort - init port structure of output device 'ps'

   SYNOPSIS:
   static void InitPSPort (OUTPUTDEVICE *thePort);

   PARAMETERS:
   .  thePort - port structure to initialize

   DESCRIPTION:
   This function inits port structure of output device 'ps'

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void InitPSPort (OUTPUTDEVICE *thePort)
{
  short r,g,b,delta,i,j,max,res;

  /* init pointers to basic drawing functions */
  thePort->Move                   = PSMoveTo;
  thePort->Draw                   = PSDrawTo;
  thePort->Polyline               = PSPolyline;
  thePort->Polygon                = PSPolygon;
  thePort->InversePolygon = PSInversePolygon;
  thePort->ErasePolygon   = PSErasePolygon;
  thePort->Polymark               = PSPolymark;
  thePort->InvPolymark    = PSInvPolymark;
  thePort->Text                   = PSText;
  thePort->CenteredText   = PSCenteredText;
  thePort->ClearViewPort  = PSClearViewPort;

  /* init pointers to set functions */
  thePort->SetLineWidth   = PSSetLineWidth;
  thePort->SetTextSize    = PSSetTextSize;
  thePort->SetMarker              = PSSetMarker;
  thePort->SetMarkerSize  = PSSetMarkerSize;
  thePort->SetColor               = PSForeground;
  thePort->SetPaletteEntry= PSSetPaletteEntry;
  thePort->SetNewPalette  = PSSetPalette;

  /* init pointers to miscellaneous functions */
  thePort->GetPaletteEntry        = PSGetPaletteEntry;
  thePort->Flush                          = PSFlush;

  /* fill port */
  thePort->black = 255;
  thePort->white = 0;
  thePort->red = 254;
  thePort->green = 128;
  thePort->blue = 2;
  thePort->cyan = 65;
  thePort->orange = 220;
  thePort->yellow = 191;
  thePort->darkyellow = 205;
  thePort->magenta = 1;
  thePort->hasPalette = 1;
  thePort->range = 256;
  thePort->spectrumStart = 2;
  thePort->spectrumEnd = 254;
  thePort->signx = 1;
  thePort->signy = 1;

  /* initialize color table */
  res = 63;
  delta = 4;
  max = 252;
  i = 0;

  /* fixed colors */
  red[i] = 255; green[i] = 255; blue[i++] = 255;        /* 0 = white */
  red[i] = 255; green[i] = 0      ; blue[i++] = 255;            /* 1 = magenta */

  /* color spectrum */
  r = g = 0; b = max;
  red[i] = r; green[i] = g; blue[i++] = b;                      /* 2 = blue */

  /* blue to cyan */
  for (j=0; j<res; j++)
  {
    g += delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 65 = cyan */
  /* cyan to green */
  for (j=0; j<res; j++)
  {
    b -= delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 128 = green */
  /* green to yellow */
  for (j=0; j<res; j++)
  {
    r += delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 191 = yellow */
  /* yellow to rot */
  for (j=0; j<res; j++)
  {
    g -= delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 254 = red */
  red[i] = 0; green[i] = 0  ; blue[i++] = 0;                    /* 255 = black */

  for (j=0; j<i; j++)
  {
    red[j]   /= 255.0;
    green[j] /= 255.0;
    blue[j]  /= 255.0;
  }

  return;
}

static INT WritePSHeader (FILE *file, const char *title, INT x, INT y, INT w, INT h)
{
  time_t tp;
  char date[64];

  /* read date */
  if (time(&tp)!=-1)
    strcpy(date,ctime(&tp));
  else
    strcpy(date,"\n");

  /* write header */
  fprintf(file,"%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(file,"%%%%Title: %s\n", title);
  fprintf(file,"%%%%Creator: %s\n",CREATOR);
  fprintf(file,"%%%%CreationDate: %s", date);
  fprintf(file,"%%%%BoundingBox: %d %d %d %d\n",
          (int)x, (int)y, (int)w, (int)h );
  fprintf(file,"%%%%Pages: 1\n");
  fprintf(file,"%%%%DocumentsFonts: %s\n",PS_DEFAULT_FONT);
  fprintf(file,"%%%%Copyright 1994 ug-group - All Rights Reserved Worldwide\n");
  fprintf(file,"%%%%EndComments\n\n");

  /* write initialisation */
  fprintf(file,"1 setlinejoin\n");
  fprintf(file,"1 setlinecap\n");

  fprintf(file,"/%s findfont %d scalefont setfont\n", PS_DEFAULT_FONT, PS_FONT_FACTOR);

  fprintf(file,"\n");

  /* defines */
  fprintf(file,"/M {moveto} def\n");
  fprintf(file,"/S {lineto stroke} def\n");
  fprintf(file,"/L {lineto} def\n");
  fprintf(file,"/C {closepath fill} def\n");
  fprintf(file,"/N {newpath} def\n");
  fprintf(file,"/R {setrgbcolor} def\n");
  fprintf(file,"/W {setlinewidth} def\n");

  fprintf(file,"\n");

  /* write end of prolog */
  fprintf(file,"%%%%Endprolog\n%%\n");
  fprintf(file,"%%%%Page: 1 1\n%%\n\n");

  return (0);
}

static INT ComputeTransformation (PSWINDOW *PSWindow,INT *ox,INT *oy,INT *sx,INT *sy, INT ls)
{
  short j;

  /* compute transformation */
  if (!ls)
  {
    /* regular format */
    tx  = (float) (*ox);
    ty  = (float) (*oy);
    mxx = (float) (*sx);
    myy = (float) (*sy);
    mxy = myx = 0.0;
  }
  else
  {
    /* landscape format */
    tx  = (float) (PAPER_X-(*ox));
    ty  = (float) (*oy);
    mxy =-(float) (*sy);
    myx = (float) (*sx);
    mxx = myy = 0.0;
    *ox = PAPER_X-(*ox)-(*sy);
    j   = *sx;
    *sx = *sy;
    *sy = j;
  }
  PSWindow->landscape = landscape = ls;
  PSWindow->tx  = tx;
  PSWindow->ty  = ty;
  PSWindow->mxx = mxx;
  PSWindow->mxy = mxy;
  PSWindow->myx = myx;
  PSWindow->myy = myy;

  return (0);
}

/****************************************************************************/
/*
   OpenDocumentWindow - Open a postscript file

   SYNOPSIS:
   static WINDOWID OpenPSWindow (char *title, INT x, INT y, INT width,
   INT height, INT *Global_LL, INT *Global_UR, INT *Local_LL, INT *Local_UR,
   INT *error);

   PARAMETERS:
   .  title -
   .  x -
   .  y -
   .  width -
   .  height -
   .  Global_LL -
   .  Global_UR -
   .  Local_LL -
   .  Local_UR -
   .  error -

   DESCRIPTION:
   This function opens a psfile.

   RETURN VALUE:
   WINDOWID
   .n   pointer to the window struct
   .n   NULL if an error occured.
 */
/****************************************************************************/

static WINDOWID OpenPSWindow (const char *title, INT x, INT y, INT width, INT height, INT *Global_LL, INT *Global_UR, INT *Local_LL, INT *Local_UR, INT *error)
{
  PSWINDOW *PSWindow;
  char pspath[BUFFLEN];
  INT sx,sy;

  *error = 0;

  /* create PSWINDOW structure */
  PSWindow = malloc(sizeof(PSWINDOW));
  if (PSWindow==NULL) {*error=1; return(0);}

  /* init postscript window */
  if (GetDefaultValue(DEFAULTSFILENAME,"psfilesdir",pspath)==0)
    PSWindow->psfile = FileOpenUsingSearchPath(title,"w",pspath);
  else
    PSWindow->psfile = fileopen(title,"w");
  if (PSWindow->psfile==NULL)
  {
    free(PSWindow);
    *error = 1;
    return (0);
  }

  /* set currents */
  currPSW = PSWindow;
  currPSF = PSWindow->psfile;

  /* return corners in devices coordinate system (before they are changed by ComputeTransformation) */
  Global_LL[0] = Local_LL[0] = x;  Global_LL[1] = Local_LL[1] = y;
  Global_UR[0] = Local_UR[0] = x + width;
  Global_UR[1] = Local_UR[1] = y + height;

  /* transformation with default origin, scale and format (how could they be made to be changable?) */
  sx = sy = 1;
  ComputeTransformation(PSWindow,&x,&y,&sx,&sy,0);

  /* init psfile */
  WritePSHeader(currPSF,title,x,y,sx*width,sy*height);

  /* write pallette */
  /* PSSetPalette(0,256,red,green,blue); */
  /* default palette */

  PSSetLineWidth(1);
  PSSetTextSize(10);

  /* return window ptr */
  return((WINDOWID)currPSW);
}

/****************************************************************************/
/*
   ClosePSPort - write trailer and close file

   SYNOPSIS:
   static INT ClosePSWindow (WINDOWID win);

   PARAMETERS:
   .  win -

   DESCRIPTION:
   This function writes the trailer and closes the file.

   RETURN VALUE:
   INT
   .n    0 if operation ok
   .n    1 if an error occured.
 */
/****************************************************************************/

static INT ClosePSWindow (WINDOWID win)
{
  currPSW = (PSWINDOW *) win;
  if (currPSW==NULL) return (1);
  currPSF = currPSW->psfile;
  if (currPSF==NULL) return (0);

  /* write trailer */
  fprintf(currPSF,"\nshowpage\n\n");
  fprintf(currPSF,"%%%%Trailer\n");

  fclose(currPSF);
  free(currPSW);

  currPSW = NULL;
  currPSF = NULL;

  return (0);
}


/****************************************************************************/
/*
   SetPSOutput - Activate the window associated with theView

   SYNOPSIS:
   static INT SetPSOutput (WINDOWID win);

   PARAMETERS:
   .  win -

   DESCRIPTION:
   This function activates the window associated with theView.

   RETURN VALUE:
   INT
   .n    0 if operation ok
   .n    1 if an error occured.
 */
/****************************************************************************/

static INT SetPSOutput (WINDOWID win)
{
  /* set current output window */
  currPSW = (PSWINDOW *) win;
  currPSF = currPSW->psfile;

  /* set current transformation */
  tx  = currPSW->tx;
  ty  = currPSW->ty;
  mxx = currPSW->mxx;
  mxy = currPSW->mxy;
  myx = currPSW->myx;
  myy = currPSW->myy;

  return(0);
}

/****************************************************************************/
/*
   UpdateOutput - Draws all controls and highlights active tool

   SYNOPSIS:
   static INT UpdatePSOutput (WINDOWID win, char *s, INT tool)

   PARAMETERS:
   .  win -
   .  s -
   .  tool

   DESCRIPTION:
   This function draws all controls and highlights active tool.

   RETURN VALUE:
   INT

   0, when OK

   1, when error
 */
/****************************************************************************/

static INT UpdatePSOutput (WINDOWID win, INT tool)
{
  return(0);
}

/****************************************************************************/
/*
   InitPSOutputDevice	- Create psfile output device

   SYNOPSIS:
   static OUTPUTDEVICE *InitPSOutputDevice (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This file creates psfile output device.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n    pointer to output device
   .n    NULL if an error occured.

 */
/****************************************************************************/

static OUTPUTDEVICE *InitPSOutputDevice (void)
{
  /* create output device */
  if ( (PSOutputDevice=CreateOutputDevice("ps")) == NULL ) return(NULL);

  /* init output device 'meta' */
  PSOutputDevice->OpenOutput       = OpenPSWindow;
  PSOutputDevice->CloseOutput      = ClosePSWindow;
  PSOutputDevice->ActivateOutput = SetPSOutput;
  PSOutputDevice->UpdateOutput     = UpdatePSOutput;

  PSOutputDevice->v.locked                 = 1;
  PSOutputDevice->PixelRatio       = (DOUBLE) 1.0;
  InitPSPort (PSOutputDevice);

  UserWrite("output device 'ps' created\n");

  return (PSOutputDevice);
}

/****************************************************************************/
/*
   InitPostScript - Initialize psfile

   SYNOPSIS:
   INT InitPS (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes psfile.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if an error occured.
 */
/****************************************************************************/

INT InitPostScript (void)
{
  if ((InitPSOutputDevice()) == NULL) return (1);

  return (0);
}
