// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  metafile.c													*/
/*																			*/
/* Purpose:   write graphics meta files                                                                         */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					        */
/*																			*/
/* History:   01.11.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
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

#define METABUFFERSIZE 16384            /* buffering graphic commmands			*/

/* meta command opcodes */
#define opNop                           0
#define opMove                          1
#define opDraw                          2
#define opPolyline                      3
#define opPolygon                       4
#define opPolymark                      5
#define opText                          6
#define opCenteredText          7
#define opSetLineWidth          8
#define opSetMarker             9
#define opSetMarkerSize         10
#define opSetTextSize           11
#define opSetColor                      12
#define opSetEntry                      13
#define opSetPalette            14
#define opInvPolygon            15
#define opErasePolygon          16

/* The following macro is only valid for long and short data types */

/* Adjust a little endian variable from the native data length to
   the standard legth */

#define ADJLEN(VALP,STD_SIZE) ((void*) (((char *) (VALP)) + (sizeof(*(VALP)) - (STD_SIZE))))

#define NAT2STDS(VAL) ((myshort = swap_short(VAL)),(ADJLEN(&(myshort),STD_SHORT)))

#define NAT2STDL(VAL) ((mylong  = swap_long(VAL)),(ADJLEN(&(mylong ) ,STD_LONG )))

#define MEMCPYS(data,val) {memcpy(data,NAT2STDS(val),STD_SHORT); data += STD_SHORT; }

/* STD_LONG and STD_SHORT are chose n the minimum possible values for   */
/* reasons of efficiency, the code relies on this fact, they may        */
/* NEVER be larger than sizeof(long) resp. sizeof(short)                */
/* This code relies on sizeof(char) being exactly one byte.             */

#define STD_LONG  4
#define STD_SHORT 2

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT littleEndian = 1; /* needed for check LITTLE/BIG-ENDIAN */
static short myshort;
static long mylong;

typedef struct {

  /* File */
  FILE *metafile;                                       /* output file							*/

  /* data management */
  char metabuffer[METABUFFERSIZE];      /* output buffer						*/
  long blockSize;                                       /* METABUFFERSIZE						*/
  long blockUsed;                                       /* actual buffer size used				*/
  long itemCounter;                                     /* number of commands in buffer                 */
  char *data;                                           /* data pointer in buffer				*/
  short xdim, ydim;                                     /* bounding box (screen)				*/
} METAWINDOW ;

static OUTPUTDEVICE *MetaOutputDevice=NULL; /* outputdevice that has been ini*/
static METAWINDOW *currMW=NULL;         /* current meta window					*/
static FILE *currMF=NULL;                       /* current meta file					*/

/* static color table for all metafiles */
static short red[256];
static short green[256];
static short blue[256];
static long saveColor;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static short swap_short (short data)
{
  char *s;
  int i;
  short res;

  if (littleEndian) return(data);

  s = (char *) &res;

  /* the condition i>0 is equal to just "i" */
  for (i = sizeof(short); i; i--)
    s[i-1] = ((char *) &data)[sizeof(short)-i];

  return(res);
}

static long swap_long (long data)
{
  char *s;
  long res;
  int i;

  if (littleEndian) return(data);

  s = (char *) &res;

  for (i = sizeof(long); i; i--)       /* the condition i>0 is equal to just "i" */
    s[i-1] = ((char *) &data)[sizeof(long)-i];

  return(res);
}

static void flush_block (void)
{
  int error;

  if (currMF==NULL) return;

  if (currMW->blockUsed>0)
  {
    error = fwrite(NAT2STDL(currMW->blockUsed),STD_LONG,1,currMF);
    if (error!=1) return;
    error = fwrite(NAT2STDL(currMW->itemCounter),STD_LONG,1,currMF);
    if (error!=1) return;
    error = fwrite(currMW->metabuffer,currMW->blockUsed,1,currMF);
    if (error!=1) return;
  }
  currMW->blockUsed = currMW->itemCounter = 0;
  currMW->data = currMW->metabuffer;
  return;
}

static void MetaMove (SHORT_POINT point)
{
  if (currMW->blockUsed+5>METABUFFERSIZE) flush_block();
  *currMW->data = opMove;
  currMW->data++;

  MEMCPYS(currMW->data,point.x);
  MEMCPYS(currMW->data,point.y);
  currMW->itemCounter++;
  currMW->blockUsed += 5;
  return;
}

static void MetaDraw (SHORT_POINT point)
{
  if (currMW->blockUsed+5>METABUFFERSIZE) flush_block();
  *currMW->data = opDraw;
  currMW->data++;
  MEMCPYS(currMW->data,point.x);
  MEMCPYS(currMW->data,point.y);
  currMW->itemCounter++;
  currMW->blockUsed += 5;
  return;
}

static void MetaPolyline (SHORT_POINT *points, INT nb)
{
  int i,size1;
  short n;

  n = (short)nb;
  if (n<2) return;
  size1 = 3+n*4;
  if (currMW->blockUsed+size1>METABUFFERSIZE) flush_block();

  *currMW->data = opPolyline;
  currMW->data++;
  MEMCPYS(currMW->data,n);
  for (i=0; i<n; i++)
    MEMCPYS(currMW->data,points[i].x);
  for (i=0; i<n; i++)
    MEMCPYS(currMW->data,points[i].y);
  currMW->itemCounter++;
  currMW->blockUsed += size1;
  return;
}

static void MetaPolygon (SHORT_POINT *points, INT nb)
{
  int i,size1;
  short n;

  n = (short)nb;
  if (n<2) return;
  size1 = 3+n*4;
  if (currMW->blockUsed+size1>METABUFFERSIZE) flush_block();

  *currMW->data = opPolygon;
  currMW->data++;
  MEMCPYS(currMW->data,n);

  for (i=0; i<n; i++)
    MEMCPYS(currMW->data,points[i].x);
  for (i=0; i<n; i++)
    MEMCPYS(currMW->data,points[i].y);
  currMW->itemCounter++;
  currMW->blockUsed += size1;
  return;
}

static void MetaInversePolygon (SHORT_POINT *points, INT nb)
{
  return;
}

static void MetaSetColor (long index);

static void MetaErasePolygon (SHORT_POINT *points, INT nb)
{
  long saved;

  saved = saveColor;
  MetaSetColor(0);
  MetaPolygon(points,nb);
  MetaSetColor(saved);
  return;
}

static void MetaPolymark (short n, SHORT_POINT *points)
{
  int i,size1;

  if (n<1) return;
  size1 = 3+n*4;
  if (currMW->blockUsed+size1>METABUFFERSIZE) flush_block();

  *currMW->data = opPolymark;
  currMW->data++;
  MEMCPYS(currMW->data,n);

  for (i=0; i<n; i++)
    MEMCPYS(currMW->data,points[i].x);
  for (i=0; i<n; i++)
    MEMCPYS(currMW->data,points[i].y);
  currMW->itemCounter++;
  currMW->blockUsed += size1;
  return;
}

static void MetaText (const char *s, INT mode)
{
  short n,size;

  n = strlen(s);
  size = 3+n;
  if (currMW->blockUsed+size>METABUFFERSIZE) flush_block();

  *currMW->data = opText;
  currMW->data++;
  MEMCPYS(currMW->data,n);
  memcpy(currMW->data,s,n);
  currMW->data += n;
  currMW->itemCounter++;
  currMW->blockUsed += size;
  return;
}

static void MetaCenteredText (SHORT_POINT point, const char *s, INT mode)
{
  short n,size;

  n = strlen(s);
  size = 7+n;
  if (currMW->blockUsed+size>METABUFFERSIZE) flush_block();

  *currMW->data = opCenteredText;
  currMW->data++;
  MEMCPYS(currMW->data,point.x);
  MEMCPYS(currMW->data,point.y);
  MEMCPYS(currMW->data,n);
  memcpy(currMW->data,s,n);
  currMW->data += n;
  currMW->itemCounter++;
  currMW->blockUsed += size;
  return;
}

static void MetaClearViewPort (void)
{
  return;
}

static void MetaSetLineWidth (short w)
{
  if (currMW->blockUsed+3>METABUFFERSIZE) flush_block();
  *currMW->data = opSetLineWidth;
  currMW->data++;
  MEMCPYS(currMW->data,w);
  currMW->itemCounter++;
  currMW->blockUsed += 3;
  return;
}

static void MetaSetTextSize (short s)
{
  if (currMW->blockUsed+3>METABUFFERSIZE) flush_block();
  *currMW->data = opSetTextSize;
  currMW->data++;
  MEMCPYS(currMW->data,s);
  currMW->itemCounter++;
  currMW->blockUsed += 3;
  return;
}

static void MetaSetMarker (short n)
{
  if (currMW->blockUsed+3>METABUFFERSIZE) flush_block();
  *currMW->data = opSetMarker;
  currMW->data++;
  MEMCPYS(currMW->data,n);
  currMW->itemCounter++;
  currMW->blockUsed += 3;
  return;
}

static void MetaSetMarkerSize (short n)
{
  if (currMW->blockUsed+3>METABUFFERSIZE) flush_block();
  *currMW->data = opSetMarkerSize;
  currMW->data++;
  MEMCPYS(currMW->data,n);
  currMW->itemCounter++;
  currMW->blockUsed += 3;
  return;
}


static void MetaSetColor (long index)
{
  unsigned char i;

  i = index%256;
  saveColor=index;
  if (currMW->blockUsed+2>METABUFFERSIZE) flush_block();
  *currMW->data = opSetColor;
  currMW->data++;
  memcpy(currMW->data,&i,1);
  currMW->data++;
  currMW->itemCounter++;
  currMW->blockUsed += 2;
  return;
}

static void MetaSetPaletteEntry (long index, short r, short g, short b)
{
  return;
}

static void MetaSetPalette (long start, long count, short *r, short *g, short *b)
{
  short size;
  long end;
  unsigned char i,j;
  short k;

  i = start%256;
  end = start+count-1;
  j = end%256;
  if (end<start) return;
  size = 3+count*3;

  if (currMW->blockUsed+size>METABUFFERSIZE) flush_block();
  *currMW->data = opSetPalette;
  currMW->data++;
  *((unsigned char *)currMW->data) = i;
  currMW->data++;
  *((unsigned char *)currMW->data) = j;
  currMW->data++;
  for (k=0; k<count; k++)
  {
    *((unsigned char *)currMW->data++) = (unsigned char) r[k];
    *((unsigned char *)currMW->data++) = (unsigned char) g[k];
    *((unsigned char *)currMW->data++) = (unsigned char) b[k];
  }
  currMW->itemCounter++;
  currMW->blockUsed += size;

  return;
}

static void MetaGetPaletteEntry (long index, short *r, short *g, short *b)
{
  return;
}

static void MetaFlush (void)
{
  return;
}

/****************************************************************************/
/*
   InitMetaPort - init port structure of output device 'meta'

   SYNOPSIS:
   static void InitMetaPort (OUTPUTDEVICE *thePort);

   PARAMETERS:
   .  thePort - port structure to initialize

   DESCRIPTION:
   This function inits port structure of output device 'meta'

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void InitMetaPort (OUTPUTDEVICE *thePort)
{
  short r,g,b,delta,i,j,max,res;

  /* init pointers to basic drawing functions */
  thePort->Move                   = MetaMove;
  thePort->Draw                   = MetaDraw;
  thePort->Polyline               = MetaPolyline;
  thePort->Polygon                = MetaPolygon;
  thePort->InversePolygon = MetaInversePolygon;
  thePort->ErasePolygon   = MetaErasePolygon;
  thePort->Polymark               = MetaPolymark;
  thePort->Text                   = MetaText;
  thePort->CenteredText   = MetaCenteredText;
  thePort->ClearViewPort  = MetaClearViewPort;

  /* init pointers to set functions */
  thePort->SetLineWidth   = MetaSetLineWidth;
  thePort->SetTextSize    = MetaSetTextSize;
  thePort->SetMarker              = MetaSetMarker;
  thePort->SetMarkerSize  = MetaSetMarkerSize;
  thePort->SetColor               = MetaSetColor;
  thePort->SetPaletteEntry= MetaSetPaletteEntry;
  thePort->SetNewPalette  = MetaSetPalette;

  /* init pointers to miscellaneous functions */
  thePort->GetPaletteEntry        = MetaGetPaletteEntry;
  thePort->Flush                          = MetaFlush;

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

  /* blau nach cyan */
  for (j=0; j<res; j++)
  {
    g += delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 65 = cyan */
  /* cyan nach green */
  for (j=0; j<res; j++)
  {
    b -= delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 128 = green */
  /* grüeen nach gelb */
  for (j=0; j<res; j++)
  {
    r += delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 191 = yellow */
  /* gelb nach rot */
  for (j=0; j<res; j++)
  {
    g -= delta;
    red[i] = r; green[i] = g; blue[i++] = b;
  }                                                             /* 254 = red */
  red[i] = 0; green[i] = 0  ; blue[i++] = 0;                    /* 255 = black */

  return;
}



/****************************************************************************/
/*
   OpenDocumentWindow - Open a metafile

   SYNOPSIS:
   static WINDOWID OpenMetaWindow (char *title, INT x, INT y, INT width,
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
   This function opens a metafile.

   RETURN VALUE:
   WINDOWID
   .n   pointer to the window struct
   .n   NULL if an error occured.
 */
/****************************************************************************/

static WINDOWID OpenMetaWindow (const char *title, INT x, INT y, INT width, INT height, INT *Global_LL, INT *Global_UR, INT *Local_LL, INT *Local_UR, INT *error)
{
  METAWINDOW *MetaWindow;
  char metapath[BUFFLEN];

  *error = 0;

  /* create METAWINDOW structure */
  MetaWindow = malloc(sizeof(METAWINDOW));
  if (MetaWindow==NULL) {*error=1; return(0);}

  /* init metawindow */
  MetaWindow->blockSize = METABUFFERSIZE;
  MetaWindow->data = MetaWindow->metabuffer;
  MetaWindow->blockUsed = 0;
  MetaWindow->itemCounter = 0;
  if (GetDefaultValue(DEFAULTSFILENAME,"metafilesdir",metapath)==0)
    MetaWindow->metafile = FileOpenUsingSearchPath(title,"wb",metapath);
  else
    MetaWindow->metafile = fileopen(title,"wb");
  if (MetaWindow->metafile==NULL)
  {
    free(MetaWindow);
    *error = 1;
    return (0);
  }
  MetaWindow->xdim = width;
  MetaWindow->ydim = height;

  /* set currents */
  currMW = MetaWindow;
  currMF = MetaWindow->metafile;

  /* init metafile */
  /* block size */
  fwrite(NAT2STDL(MetaWindow->blockSize),STD_LONG,1,MetaWindow->metafile);
  /* x size     */
  fwrite(NAT2STDS(MetaWindow->xdim),STD_SHORT,1,MetaWindow->metafile);
  /* y size     */
  fwrite(NAT2STDS(MetaWindow->ydim),STD_SHORT,1,MetaWindow->metafile);

  /* write pallette */
  MetaSetPalette(0,256,red,green,blue);                         /* default palette */

  /* return corners in devices coordinate system */
  Global_LL[0] = Local_LL[0] = x;  Global_LL[1] = Local_LL[1] = y;
  Global_UR[0] = Local_UR[0] = x + width;
  Global_UR[1] = Local_UR[1]= y + height;

  /* return window ptr */
  return((WINDOWID)currMW);
}

/****************************************************************************/
/*
   CloseMetaPort - Flush last block and close file

   SYNOPSIS:
   static INT CloseMetaWindow (WINDOWID win);

   PARAMETERS:
   .  win -

   DESCRIPTION:
   This function flushes last block and closes file.

   RETURN VALUE:
   INT
   .n    0 if operation ok
   .n    1 if an error occured.
 */
/****************************************************************************/

static INT CloseMetaWindow (WINDOWID win)
{
  currMW = (METAWINDOW *) win;
  if (currMW==NULL) return (1);
  currMF = currMW->metafile;
  if (currMF==NULL) return (1);

  flush_block();
  fclose(currMF);
  free(currMW);

  currMW = NULL;
  currMF = NULL;

  return (0);
}


/****************************************************************************/
/*
   SetMetaOutput - Activate the window associated with theView

   SYNOPSIS:
   static INT SetMetaOutput (WINDOWID win);

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

static INT SetMetaOutput (WINDOWID win)
{
  /* set current output window */
  currMW = (METAWINDOW *) win;
  currMF = currMW->metafile;

  return(0);
}

/****************************************************************************/
/*
   UpdateOutput - Draws all controls and highlights active tool

   SYNOPSIS:
   static INT UpdateMetaOutput (WINDOWID win, char *s, INT tool)

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

static INT UpdateMetaOutput (WINDOWID win, char *s, INT tool)
{
  return(0);
}

/****************************************************************************/
/*
   InitMetaOutputDevice	- Create metafile output device

   SYNOPSIS:
   static OUTPUTDEVICE *InitMetaOutputDevice (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This file creates metafile output device.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n    pointer to output device
   .n    NULL if an error occured.

 */
/****************************************************************************/

static OUTPUTDEVICE *InitMetaOutputDevice (void)
{
  /* create output device */
  if ( (MetaOutputDevice=CreateOutputDevice("meta")) == NULL ) return(NULL);

  /* init output device 'meta' */
  MetaOutputDevice->OpenOutput     = OpenMetaWindow;
  MetaOutputDevice->CloseOutput    = CloseMetaWindow;
  MetaOutputDevice->ActivateOutput = SetMetaOutput;
  MetaOutputDevice->UpdateOutput   = UpdateMetaOutput;

  MetaOutputDevice->v.locked               = 1;
  MetaOutputDevice->PixelRatio     = (COORD) 1.0;
  InitMetaPort (MetaOutputDevice);

  UserWrite("output device 'meta' created\n");

  return (MetaOutputDevice);
}

/****************************************************************************/
/*
   InitMeta - Initialize metafile

   SYNOPSIS:
   INT InitMeta (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes metafile.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if an error occured.
 */
/****************************************************************************/

INT InitMeta (void)
{
  /* check for little/big endian storage type */
  littleEndian = !( *((char *) &littleEndian));

  if ((InitMetaOutputDevice()) == NULL) return (1);

  return (0);
}
