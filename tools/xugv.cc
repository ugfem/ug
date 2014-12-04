// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*									                                                                                */
/* File:      xugv.c														*/
/*									                                                                                */
/* Purpose:   ug metafile display program for X11                                       */
/*									                                                                                */
/* Author:    Peter Bastian						                                                        */
/*	    Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen                 */
/*	    Universitaet Heidelberg				                                                        */
/*	    Im Neuenheimer Feld 368				                                                        */
/*	    6900 Heidelberg					                                                                */
/*	    internet: bastian@iwr1.iwr.uni-heidelberg.de		                                */
/*									                                                                                */
/*	    Stefan Lang						                                                                        */
/*	    Institut fuer Mathematische Maschinen und		                                */
/*	    Datenverarbeitung III					                                                        */
/*	    Martensstrasse 3					                                                        */
/*	    91058 Erlangen														*/
/*									                                                                                */
/*									                                                                                */
/*									                                                                                */
/* History:   22.10.93 begin												*/
/*									                                                                                */
/* Remarks:   up to this time xugv needs a color screen		                        */
/*									                                                                                */
/****************************************************************************/

/* ugview is a program for displaying ug metafont files on sun workstations */
/* with color screen.							    */


/* general header files */
#include <config.h>
#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <cmath>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <X11/Intrinsic.h>
#include <X11/IntrinsicP.h>
#include <X11/Vendor.h>
#include <X11/StringDefs.h>
#include <X11/keysym.h>
#include <X11/Xaw/XawInit.h>
#include <X11/Xaw/Viewport.h>
#include <X11/Xaw/Simple.h>
#include <X11/Xaw/Label.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Xaw/Command.h>
#include "dev/xif/shades.h" /* !!! */
#include "general.h"

/* defines for opCodes */
#define opNop                           0
#define opMove                  1
#define opDraw                  2
#define opPolyline              3
#define opPolygon               4
#define opPolymark              5
#define opText                  6
#define opCenteredText          7
#define opSetLineWidth          8
#define opSetMarker                     9
#define opSetMarkerSize    10
#define opSetTextSize      11
#define opSetColor             12
#define opSetEntry             13
#define opSetPalette       14
#define opNewLine              15
#define opNewPolyline      16
#define opNewPolygon       17
#define opNewPolymark      18
#define opNewText              19
#define opNewCenteredText  20
#define opShadedPolygon    21

#define SIZE                50
#define CSIZE              256

#define VERBOSE 0

#define MAX_FILES       5

/* needed only for CRAY C90!!					*/
/* because there is in error in fread() for short */
#ifdef __C90__
#define short long
#endif

/* The following macro is only valid for long and short data types */

#define SWAP(X) (littleEndian ? (X) : \
                 ((sizeof(X)!=sizeof(short)) ? swap_long(X) : swap_short(X)))


/* Adjust a little endian var. from the native data length to the standard length */

#define ADJLEN(VALP,STD_SIZE) \
  ((void*) ((char *) (VALP))[sizeof(*(VALP)) - STD_SIZE])

#define INVADJLEN(VALP,STD_SIZE) \
  ((*(VALP))>>((sizeof(*(VALP)) - STD_SIZE)*8))

#define STD2NAT(VAL,STDSIZE) VAL = INVADJLEN(&(VAL),STDSIZE);VAL = SWAP(VAL)

/* STD_LONG and STD_SHORT are chosen the minimum possible values for */
/* reasons of efficiency, the code relies on this fact, they may */
/* NEVER be larger than sizeof(long) resp. sizeof(short) */
/* This code relies on sizeof(char) being exactly one byte */

#define STD_LONG 4
#define STD_SHORT 2

#define XUGV_MAX(a,b) ((a)<(b)) ? (b) : (a)
#define XUGV_MIN(a,b) ((a)<(b)) ? (a) : (b)

/* defines for max. initial windowsize */
#define WIN_HEIGHT                864
#define WIN_WIDTH                1152

/* class name of the application */
#define APPL_KLASSE    "Xugview"

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


XtAppContext kontext;                     /* application kontext */
Display      *display;                  /* display of the application */
int screen;
Colormap cmap, def_colormap;       /* colormap used inside window */
XColor colortable[CSIZE];          /* colortable for indirection  */
short used[CSIZE];                              /* those colors which are used but not allocated */
GC gc;                                                     /* graphics context */
Widget applShell,                               /* Widgets */
       viewport,
       dial,
       xexit,
       picture;
Pixmap pixmap;                 /* pixmap to store picture */
FILE     *stream;
FILE *mstream[MAX_FILES];
char     *newfile;
const char       *option;
short fx, fy;          /* dimensions of drawing object */
short mfx[MAX_FILES], mfy[MAX_FILES];   /* dimensions of drawing object */
short fx_max, fy_max;
int error;
int expose;
Cardinal n;
String file;
short vwidth,
      vheight,
      pwidth,
      pheight;
int dispcells;            /* number of available displaycells */
int incr=1,first,last,film,frame_number;
int _wait=0;                      /* wait for file creation, if file does not exist */
int ignore=0;
int verbose=0;
char frame[50];
int p_opt=0, g_opt=0;
char outext[80];
int through=0;
int stoploop = 0;
int count = 0;
int n_pic, i_pic, f_offset, nbreak;
char *mfile[MAX_FILES];
static int run_max, run_count;
static unsigned int sleep_seconds = 1;
static int auto_nb;

static int littleEndian = 1; /* needed for check LITTLE/BIG-ENDIAN */

/* pixel structure for true color visual */
unsigned long red_mask, green_mask, blue_mask;
int red_shift, green_shift, blue_shift;
int true_color;

/* pixmaps for shading patterns */
static Pixmap pattern[NO_PATTERNS];

/* forward declarations */
static void Marker(short n, short s, short x, short y);
void createGraphics(void);

/* reverse byte order of a short variable independent of sizeof(short) */
static short swap_short (short data)
{
  char *s;
  int i;
  short res;

  s = (char *) &res;

  for (i = sizeof(short); i; i--)       /* the condition i>0 is equal to just "i" */
    s[i-1] = ((char *) &data)[sizeof(short)-i];

  return(res);
}

/* reverse byte order of a long variable independent of sizeof(long) */
static long swap_long (long data)
{
  char *s;
  int i;
  long res;

  s = (char *) &res;

  for (i = sizeof(long); i; i--)       /* the condition i>0 is equal to just "i" */
    s[i-1] = ((char *) &data)[sizeof(long)-i];

  return(res);
}

FILE *fopen_with_wait (const char *name, const char *type)
{
  FILE *stream = NULL;

  do
  {
    stream = fopen(name,type);
    if (stream == NULL) sleep(sleep_seconds);
    else
    {
      if (verbose)
      {
        printf ("xugv: opened file %s\n",name);
        fflush(stdout);
      }
      break;
    }
  }
  while (1);

  return(stream);
}

int GetFileScreen (FILE *stream, short *fx, short *fy)
{
  long blockSize;          /* METABUFFERSIZE */
  int error;
  short newfx;

  if (VERBOSE)
    printf("GetFileScreen():\n");

  /* get file parameters */
  rewind(stream);
  error = fread(&blockSize,STD_LONG,1,stream);
  if (VERBOSE)
    printf("error=%d blockSize=%d hex=%x\n",error,blockSize,blockSize);
  STD2NAT(blockSize,STD_LONG);
  if (error!=1) return(1);                    /* block size */
  error = fread(fx,STD_SHORT,1,stream);
  if (VERBOSE)
    printf("error=%d fx=%d hex=%x\n",error,*fx,*fx);
  if (VERBOSE>2)
  {
    newfx = *fx;
    newfx = swap_short(newfx);
    printf("swap_short: newfx=%d\n",newfx);
    newfx = *fx;
    newfx = SWAP(newfx);
    printf("SWAP: newfx=%d\n",newfx);
    newfx = INVADJLEN(&newfx,STD_SHORT);
    printf("SWAP and INVADJLEN: newfx=%d\n",newfx);
    newfx = *fx;
    newfx = INVADJLEN(&newfx,STD_SHORT);
    printf("INVADJLEN: newfx=%d\n",newfx);
    newfx = SWAP(newfx);
    printf("INVADJLEN and SWAP: newfx=%d\n",newfx);
    newfx = *fx;
    STD2NAT(newfx,STD_SHORT);
    printf("STD2NAT: newfx=%d\n",newfx);
  }

  STD2NAT(*fx,STD_SHORT);
  if (error!=1) return(1);                    /* x size */
  error = fread(fy,STD_SHORT,1,stream);
  if (VERBOSE)
    printf("error=%d fy=%d hex=%x\n",error,*fy,*fy);
  STD2NAT(*fy,STD_SHORT);
  if (error!=1) return(1);                    /* y size */

  if (VERBOSE)
  {
    printf("sizeof blockSize=%d fx=%d fy=%d\n",sizeof(blockSize),sizeof(*fx),sizeof(*fy));
    printf("blockSize=%d fx=%d fy=%d\n",blockSize,*fx,*fy);
    printf("blockSize=%x fx=%x fy=%x\n",blockSize,*fx,*fy);
  }

  return (0);
}

int get_component_shift(unsigned long mask)
{
  int shift;

  shift=0;
  while (!(mask & 1)) {
    shift++;
    mask >>= 1;
  }
  return shift;
}

void createGraphics(void)
{
  int i;
  int ncolors;
  XColor colors[CSIZE];
  XGCValues gcv;
  Visual *visual;
  XPoint edges[4];

  if (VERBOSE)
    printf("createGraphics():");

  /* get default screen */
  screen = DefaultScreen(display);

  /* test if true color visual */
  visual = DefaultVisual(display, screen);
#ifdef __cplusplus
  true_color = (visual->c_class == TrueColor);
#else
  true_color = (visual->class == TrueColor);
#endif

  /* get info on pixel structure */
  if (true_color) {
    red_mask   = visual->red_mask;
    green_mask = visual->green_mask;
    blue_mask  = visual->blue_mask;
    red_shift  = get_component_shift(red_mask);
    green_shift= get_component_shift(green_mask);
    blue_shift = get_component_shift(blue_mask);
  }

  /* get one supported visual type and create colormap */
  ncolors = DisplayCells(display, screen);
  if (ncolors > CSIZE) ncolors = CSIZE;
  def_colormap = DefaultColormap(display, screen);
  for( i=0; i<ncolors; i++) {
    colors[i].pixel = i;
    colors[i].flags = DoRed|DoGreen|DoBlue;
  }

  /* create graphic context */
  gcv.foreground = WhitePixel(display, screen);
  gcv.background = WhitePixel(display, screen);
  gcv.line_width = 1;

  gc = XCreateGC( XtDisplay(picture) , XtWindow(picture),
                  GCForeground|GCBackground|GCLineWidth, &gcv);

  /* create a pixmap to store the picture */
  pixmap = XCreatePixmap(display, XtWindow(picture), pwidth, pheight,
                         DefaultDepthOfScreen(XtScreen(picture)));


  /* clear pixmap's area */
  edges[0].x = 0;
  edges[0].y = 0;
  edges[1].x = 0;
  edges[1].y = pheight;
  edges[2].x = pwidth;
  edges[2].y = pheight;
  edges[3].x = pwidth;
  edges[3].y = 0;

  XFillPolygon(display, pixmap, gc, edges, 4,
               Convex, CoordModeOrigin);

  /* make stipples for shading */
  if (!true_color)
    for (i = 0; i < NO_PATTERNS; i++)
      pattern[i] = XCreateBitmapFromData(display,XtWindow(picture),
                                         pattern_data[i], PATTERN_SIZE, PATTERN_SIZE);
}


/* exposeCB is called at any expose event */
void exposeCB (Widget widget, XtPointer client_data, XEvent *event, Boolean *continue_to_dispatch)
{
  XCopyArea(display, pixmap, XtWindow(picture), gc, 0, 0, pwidth, pheight, 0, 0);
}

/* createcolors creates or changes the colormapping table. For each color  */
/* requested the pixelvalue of the corresponding color in the colormap is  */
/* stored in the colormapping table.					                   */

void createColors(XColor colors[CSIZE], short x, short y)
{
  short i, j, strip, alloc;
  unsigned short mask;
  long d;
  long maxdist, actdist;
  XColor ctab[CSIZE];

  strip = 0;
  mask = ~0x0<<8;
  maxdist = 0x1L<<9;

  actdist = maxdist;

  /* get current colormap values */
  for (i=0; i<CSIZE; i++) ctab[i].pixel = i;
  XQueryColors(display, def_colormap, ctab, dispcells);

  /* loop until very modest color request */
  while (mask)
  {
    /* printf("trying to allocate colors with %d stripped bits:\n", strip); */

    /* all current color requests */
    for (i=x; i<=y; i++)
    {
      /* try to allocate color in colormap */
      alloc = XAllocColor(display, def_colormap, &colors[i]);
      if (alloc)
      {
        if (option[1] == 'v')
          printf("colortable entry %d allocated colormap pixel %d\n", i, colors[i].pixel);
      }
      else
      {
        /* look for color with similar rgb values */
        for (j=0; j<dispcells; j++)
        {
          /* calculate distance */
          d = abs(ctab[j].red - colors[i].red)
              + abs(ctab[j].green - colors[i].green)
              + abs(ctab[j].blue - colors[i].blue);

          /* accept color */
          if (d<actdist)
          {
            colors[i].pixel = j;
            used[i] = 1;
            actdist = d;
          }
          else
          if (option[1] == 'v' && option[2] == '2')
            printf("color request %d for pixel %d failed with distance %ld\n:", i, j, d);

          /* print pixel allocated for current tableentry */
          if ( used[i] && j == dispcells-1 && option[1] == 'v')
            printf("colortable entry %d uses colormap pixel %d with distance %ld\n",
                   i, colors[i].pixel, actdist);
        }

        /* reset distance */
        actdist = maxdist;

        /* one color request fails. free all allocated colors,  */
        /* strip least significant bit of all color requests    */
        /* and try to allocate colors again.			*/
        if (!used[i])
        {
          for (j=x; j<i; j++)
          {
            if (!used[j])
              XFreeColors(display, def_colormap, &(colors[j].pixel), 1, 0L);
            else
              used[j] = 0;
            colors[j].pixel = 0;
          }

          strip++;
          mask = mask<<1;
          maxdist = maxdist<<1;
          actdist = maxdist;

          /* strip least significant bit of all color requests*/
          printf("striping bit positon %d\n", strip);
          for (j=x; j<=y; j++)
          {
            colors[j].red = colors[j].red & mask;
            colors[j].green = colors[j].green & mask;
            colors[j].blue = colors[j].blue & mask;
          }
          break;
        }
      }
    }

    if (i>=y+1)
    {
      /* color request was sucessfull, store color values
                 in colortable */
      if (strip)
        printf("colorrequest was sucessfull:\n %d bits stripped...\n", strip);
      /* printf("storing allocated colors in internal colortable!\n"); */

      for(--i; i>=0; i--)
      {
        colortable[i].pixel = colors[i].pixel;
        if (!true_color) {
          colortable[i].red = ctab[colors[i].pixel].red;
          colortable[i].green = ctab[colors[i].pixel].green;
          colortable[i].blue = ctab[colors[i].pixel].blue;
        }
      }
      return;
    }
  }
  printf("Sorry couldn't allocate requested color(s),\n\
		please close other X applications!\n");
}




/* RasterizeFile reads the file to draw and executes the */
/* drawing commands.					 */
int RasterizeFile(FILE *stream)
{
  char *buffer;                         /* input buffer */
  long blockSize;                     /* METABUFFERSIZE */
  long blockUsed;                     /* actual buffer size used */
  long itemCounter;                         /* number of commands in buffer */
  char *data;                             /* data pointer in buffer */
  short fx, fy;                          /* file screen size */
  int i,error,j,k,size,snb;
  char opCode;
  short xx[SIZE],yy[SIZE];
  short x_cur, y_cur, x, y,r,g,b,n,lw,ts,m,ms,w,x1,y1,x2,y2,shd;
  XPoint xy[SIZE];
  char s[CSIZE];
  unsigned char c,ac;
  const char *fn;
  char string[10];
  char **list;
  XGCValues gcv;
  XColor color;
  XColor colors[CSIZE];
  XFontStruct *font=NULL;

  if (VERBOSE)
    printf("RasterizeFile():");

  /* get file parameters */
  rewind(stream);
  error = fread(&blockSize,STD_LONG,1,stream);
  STD2NAT(blockSize,STD_LONG);
  if (error!=1) return(1);                    /* block size */
  error = fread(&fx,STD_SHORT,1,stream);
  STD2NAT(fx,STD_SHORT);
  if (error!=1) return(1);                    /* x size */
  error = fread(&fy,STD_SHORT,1,stream);
  STD2NAT(fy,STD_SHORT);
  if (error!=1) return(1);                    /* y size */

  /* default values */
  lw = 1;
  ts = 12;
  m = 0;
  ms = 6;
  fn = string;


  /* allocate input buffer */
  buffer = (char*)malloc(blockSize);
  if (buffer==NULL) return(1);

  /* loop through the blocks */
  while (!feof(stream))
  {
    /* read block parameters */
    error = fread(&blockUsed,STD_LONG,1,stream);
    STD2NAT(blockUsed,STD_LONG);
    if (error!=1) {free(buffer); return(1);}
    error = fread(&itemCounter,STD_LONG,1,stream);
    STD2NAT(itemCounter,STD_LONG);
    if (error!=1) {free(buffer); return(1);}
    error = fread(buffer,blockUsed,1,stream);
    if (error!=1) {free(buffer); return(1);}

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
        memcpy(&x_cur,data,STD_SHORT);
        STD2NAT(x_cur,STD_SHORT);
        data += STD_SHORT;
        memcpy(&y_cur,data,STD_SHORT);
        STD2NAT(y_cur,STD_SHORT);
        y_cur = fy - y_cur;
        data += STD_SHORT;
        break;

      case opDraw :
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy - y;
        data += STD_SHORT;

        XDrawLine(display, pixmap, gc,
                  x_cur, y_cur, x, y);
        x_cur = x;
        y_cur = y;
        break;

      case opPolyline :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(STD_SHORT);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          data += STD_SHORT;
        }

        for (j=1; j<n; j++)
          XDrawLine(display, pixmap, gc,
                    xx[j-1],yy[j-1],
                    yy[j],yy[j]);
        x_cur = xx[n-1];
        y_cur = yy[n-1];
        break;

      case opPolygon :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          data += STD_SHORT;
        }

        if (n<3) break;
        for (j=0; j<n; j++) {
          xy[j].x = xx[j];
          xy[j].y = yy[j];
        }

        XFillPolygon(display, pixmap, gc,
                     xy, n, Convex, CoordModeOrigin);
        x_cur = xy[n-1].x;
        y_cur = xy[n-1].y;
        break;

      case opPolymark :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          data += STD_SHORT;
        }

        for (j=0; j<n; j++) Marker(m,ms,xx[j],yy[j]);
        x_cur = xx[n-1];
        y_cur = yy[n-1];

        break;

      case opText :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opCenteredText :
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy - y;
        data += STD_SHORT;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;

        w = XTextWidth( font, s, n);
        x_cur = x - w/2;
        y_cur = y + ts/2;
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opSetLineWidth :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        lw = n;
        gcv.line_width = n;
        XChangeGC( display, gc, GCLineWidth, &gcv);
        break;

      case opSetTextSize :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        ts = n;

        fn = "?x";
        if (n<10) string[4] = 0;else string[5] = 0;
        list = XListFonts(display, "?x", 1, &k);
        if (NULL != list) XLoadFont(display, list[0]);

        break;

      case opSetMarker :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        m = n;
        break;



      case opSetMarkerSize :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        ms = n;
        break;

      case opSetColor :
        ac = *((unsigned char *)data);
        data++;
        gcv.foreground = colortable[ac].pixel;
        XChangeGC( display, gc, GCForeground, &gcv);
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
        color.pixel = 0xffff;
        color.red = r<<8;
        color.green = g<<8;
        color.blue = b<<8;
        color.flags = DoRed|DoGreen|DoBlue;
        if (!ignore) createColors(&color, c, c);
        break;

      case opSetPalette :
        x = (short) (*((unsigned char *)data));
        data++;
        y = (unsigned short) (*((unsigned char *)data));
        data++;

        for (j=x; j<=y; j++)
        {
          r = (unsigned short) (*((unsigned char *)data));
          data++;
          g = (unsigned short) (*((unsigned char *)data));
          data++;
          b = (unsigned short) (*((unsigned char *)data));
          data++;
          colors[j].pixel = 0;
          colors[j].red = r<<8;
          colors[j].green = g<<8;
          colors[j].blue = b<<8;
          colors[j].flags = DoRed|DoGreen|DoBlue;
        }

        if (!ignore) createColors(colors, x, y);
        break;

      case opNewLine :
        lw = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x1,data,STD_SHORT);
        STD2NAT(x1,STD_SHORT);
        data += STD_SHORT;
        memcpy(&y1,data,STD_SHORT);
        STD2NAT(y1,STD_SHORT);
        y1 = fy - y1;
        data += STD_SHORT;
        memcpy(&x2,data,STD_SHORT);
        STD2NAT(x2,STD_SHORT);
        data += STD_SHORT;
        memcpy(&y2,data,STD_SHORT);
        STD2NAT(y2,STD_SHORT);
        y2 = fy - y2;
        data += STD_SHORT;
        gcv.line_width = lw;
        gcv.foreground = colortable[c].pixel;
        XChangeGC( display, gc,
                   (GCForeground|GCLineWidth), &gcv);
        XDrawLine( display, pixmap, gc, x1, y1, x2, y2);
        x_cur = x2;
        y_cur = y2;
        break;

      case opNewPolyline :
        lw = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          data += STD_SHORT;
        }
        gcv.line_width = lw;
        gcv.foreground = colortable[c].pixel;
        XChangeGC( display, gc,
                   (GCForeground|GCLineWidth),  &gcv);
        gcv.foreground = 0;
        XGetGCValues(display, gc, GCForeground, &gcv);

        for (j=1; j<n; j++)
          XDrawLine(display, pixmap, gc,
                    xx[j-1],fy-yy[j-1],xx[j],fy-yy[j]);
        x_cur = xx[n-1];
        y_cur = fy-yy[n-1];
        break;

      case opNewPolygon :
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          data += STD_SHORT;
        }
        if (n<3) break;

        gcv.line_width = lw;
        gcv.foreground = colortable[c].pixel;
        gcv.fill_style = FillSolid;
        XChangeGC( display, gc,
                   GCForeground|GCLineWidth|GCFillStyle, &gcv);

        for (j=0; j<n; j++) {
          xy[j].x = xx[j];
          xy[j].y = fy-yy[j];
        }

        XFillPolygon(display, pixmap, gc,
                     xy, n, Convex, CoordModeOrigin);
        x_cur = xy[n-1].x;
        y_cur = xy[n-1].y;
        break;

      case opNewPolymark :
        m = (short) *((unsigned char *)data);
        data++;
        ms = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          data += STD_SHORT;
        }
        gcv.line_width = lw;
        gcv.foreground =colortable[c].pixel;
        XChangeGC( display, gc,
                   (GCForeground|GCLineWidth), &gcv);
        for (j=0; j<n; j++) Marker(m,ms,xx[j],fy-yy[j]);
        x_cur = xx[n-1];
        y_cur = fy-yy[n-1];

        break;

      case opNewText :
        ts = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy - y;
        data += STD_SHORT;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        gcv.foreground =colortable[c].pixel;
        XChangeGC( display, gc,
                   GCForeground, &gcv);
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opNewCenteredText :
        ts = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy - y;
        data += STD_SHORT;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        gcv.foreground =colortable[c].pixel;
        XChangeGC( display, gc,
                   GCForeground , &gcv);
        w = XTextWidth( font, s, n);
        x_cur = x - w/2;
        y_cur = y + ts/2;
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opShadedPolygon :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        memcpy(&shd,data,STD_SHORT);
        STD2NAT(shd,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          data += STD_SHORT;
        }

        if (n<3) break;
        for (j=0; j<n; j++) {
          xy[j].x = xx[j];
          xy[j].y = yy[j];
        }
        if (!true_color) {
          snb = (int)(0.5+(NO_PATTERNS-1)*shd/1000.0);
          XSetBackground(display ,gc, BlackPixel(display, screen));
          XSetFillStyle(display, gc, FillOpaqueStippled);
          XSetStipple(display, gc, pattern[snb]);
          XFillPolygon(display, pixmap, gc,
                       xy, n, Convex, CoordModeOrigin);
          XSetFillStyle(display, gc, FillSolid);
          XSetBackground(display ,gc, WhitePixel(display, screen));
        }
        else {
          unsigned long pixel;
          float red, green, blue;

          pixel = colortable[ac].pixel;
          red   = (pixel & red_mask  ) >> red_shift;
          green = (pixel & green_mask) >> green_shift;
          blue  = (pixel & blue_mask ) >> blue_shift;
          red   *= shd/1000.0;
          green *= shd/1000.0;
          blue  *= shd/1000.0;
          pixel = ((unsigned long)(red  +0.5) << red_shift  ) +
                  ((unsigned long)(green+0.5) << green_shift) +
                  ((unsigned long)(blue +0.5) << blue_shift );
          XSetForeground(display, gc, pixel);
          XFillPolygon(display, pixmap, gc, xy, n, Convex, CoordModeOrigin);
          XSetForeground(display, gc, colortable[ac].pixel);
        }
        x_cur = xy[n-1].x;
        y_cur = xy[n-1].y;
        break;

      default :
        break;
      }
    }
  }
  fclose(stream);
  return(0);
}


/* RasterizeFile reads the file to draw and executes the */
/* drawing commands.					 */
int RasterizePositionedFile(FILE *stream, long x_offset, long y_offset)
{
  char *buffer;                         /* input buffer */
  long blockSize;                     /* METABUFFERSIZE */
  long blockUsed;                     /* actual buffer size used */
  long itemCounter;                         /* number of commands in buffer */
  char *data;                             /* data pointer in buffer */
  short fx, fy;                          /* file screen size */
  int i,error,j,k,size,snb;
  char opCode;
  short xx[SIZE],yy[SIZE];
  short x_cur, y_cur, x, y,r,g,b,n,lw,ts,m,ms,w,x1,y1,x2,y2,shd;
  XPoint xy[SIZE];
  char s[CSIZE];
  unsigned char c,ac;
  const char *fn;
  char string[10];
  char **list;
  XGCValues gcv;
  XColor color;
  XColor colors[CSIZE];
  XFontStruct *font=NULL;

  if (VERBOSE)
    printf("RasterizeFile():");

  /* get file parameters */
  rewind(stream);
  error = fread(&blockSize,STD_LONG,1,stream);
  STD2NAT(blockSize,STD_LONG);
  if (error!=1) return(1);                    /* block size */
  error = fread(&fx,STD_SHORT,1,stream);
  STD2NAT(fx,STD_SHORT);
  if (error!=1) return(1);                    /* x size */
  error = fread(&fy,STD_SHORT,1,stream);
  STD2NAT(fy,STD_SHORT);
  if (error!=1) return(1);                    /* y size */

  /* default values */
  lw = 1;
  ts = 12;
  m = 0;
  ms = 6;
  fn = string;


  /* allocate input buffer */
  buffer = (char*)malloc(blockSize);
  if (buffer==NULL) return(1);

  /* loop through the blocks */
  while (!feof(stream))
  {
    /* read block parameters */
    error = fread(&blockUsed,STD_LONG,1,stream);
    STD2NAT(blockUsed,STD_LONG);
    if (error!=1) {free(buffer); return(1);}
    error = fread(&itemCounter,STD_LONG,1,stream);
    STD2NAT(itemCounter,STD_LONG);
    if (error!=1) {free(buffer); return(1);}
    error = fread(buffer,blockUsed,1,stream);
    if (error!=1) {free(buffer); return(1);}

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
        memcpy(&x_cur,data,STD_SHORT);
        STD2NAT(x_cur,STD_SHORT);
        x_cur += x_offset;
        data += STD_SHORT;
        memcpy(&y_cur,data,STD_SHORT);
        STD2NAT(y_cur,STD_SHORT);
        y_cur = fy-y_cur;
        y_cur += y_offset;
        data += STD_SHORT;
        break;

      case opDraw :
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        x += x_offset;
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy-y;
        y += y_offset;
        data += STD_SHORT;

        XDrawLine(display, pixmap, gc,
                  x_cur, y_cur, x, y);
        x_cur = x;
        y_cur = y;
        break;

      case opPolyline :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(STD_SHORT);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          xx[k] += x_offset;
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          yy[k] += y_offset;
          data += STD_SHORT;
        }

        for (j=1; j<n; j++)
          XDrawLine(display, pixmap, gc,
                    xx[j-1],yy[j-1],
                    xx[j],yy[j]);
        x_cur = xx[n-1];
        y_cur = yy[n-1];
        break;

      case opPolygon :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          xx[k] += x_offset;
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          yy[k] += y_offset;
          data += STD_SHORT;
        }

        if (n<3) break;
        for (j=0; j<n; j++) {
          xy[j].x = xx[j];
          xy[j].y = yy[j];
        }

        XFillPolygon(display, pixmap, gc,
                     xy, n, Convex, CoordModeOrigin);
        x_cur = xy[n-1].x;
        y_cur = xy[n-1].y;
        break;

      case opPolymark :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          xx[k] += x_offset;
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          yy[k] += y_offset;
          data += STD_SHORT;
        }

        for (j=0; j<n; j++) Marker(m,ms,xx[j],yy[j]);
        x_cur = xx[n-1];
        y_cur = yy[n-1];

        break;

      case opText :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opCenteredText :
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        x += x_offset;
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy-y;
        y += y_offset;
        data += STD_SHORT;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;

        w = XTextWidth( font, s, n);
        x_cur = x - w/2;
        y_cur = y + ts/2;
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opSetLineWidth :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        lw = n;
        gcv.line_width = n;
        XChangeGC( display, gc, GCLineWidth, &gcv);
        break;

      case opSetTextSize :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        ts = n;

        fn = "?x";
        if (n<10) string[4] = 0;else string[5] = 0;
        list = XListFonts(display, "?x", 1, &k);
        if (NULL != list) XLoadFont(display, list[0]);

        break;

      case opSetMarker :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        m = n;
        break;

      case opSetMarkerSize :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        ms = n;
        break;

      case opSetColor :
        ac = *((unsigned char *)data);
        data++;
        gcv.foreground = colortable[ac].pixel;
        XChangeGC( display, gc, GCForeground, &gcv);
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
        color.pixel = 0xffff;
        color.red = r<<8;
        color.green = g<<8;
        color.blue = b<<8;
        color.flags = DoRed|DoGreen|DoBlue;
        if (!ignore) createColors(&color, c, c);
        break;

      case opSetPalette :
        x = (short) (*((unsigned char *)data));
        data++;
        y = (unsigned short) (*((unsigned char *)data));
        data++;

        for (j=x; j<=y; j++)
        {
          r = (unsigned short) (*((unsigned char *)data));
          data++;
          g = (unsigned short) (*((unsigned char *)data));
          data++;
          b = (unsigned short) (*((unsigned char *)data));
          data++;
          colors[j].pixel = 0;
          colors[j].red = r<<8;
          colors[j].green = g<<8;
          colors[j].blue = b<<8;
          colors[j].flags = DoRed|DoGreen|DoBlue;
        }

        if (!ignore) createColors(colors, x, y);
        break;

      case opNewLine :
        lw = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x1,data,STD_SHORT);
        STD2NAT(x1,STD_SHORT);
        x1 += x_offset;
        data += STD_SHORT;
        memcpy(&y1,data,STD_SHORT);
        STD2NAT(y1,STD_SHORT);
        y1 = fy - y1;
        y1 += y_offset;
        data += STD_SHORT;
        memcpy(&x2,data,STD_SHORT);
        STD2NAT(x2,STD_SHORT);
        x2 += x_offset;
        data += STD_SHORT;
        memcpy(&y2,data,STD_SHORT);
        STD2NAT(y2,STD_SHORT);
        y2 = fy - y2;
        y2 += y_offset;
        data += STD_SHORT;
        gcv.line_width = lw;
        gcv.foreground = colortable[c].pixel;
        XChangeGC( display, gc,
                   (GCForeground|GCLineWidth), &gcv);
        XDrawLine( display, pixmap, gc, x1, y1, x2, y2);
        x_cur = x2;
        y_cur = y2;
        break;

      case opNewPolyline :
        lw = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          xx[k] += x_offset;
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] += y_offset;
          data += STD_SHORT;
        }
        gcv.line_width = lw;
        gcv.foreground = colortable[c].pixel;
        XChangeGC( display, gc,
                   (GCForeground|GCLineWidth),  &gcv);
        gcv.foreground = 0;
        XGetGCValues(display, gc, GCForeground, &gcv);

        for (j=1; j<n; j++)
          XDrawLine(display, pixmap, gc,
                    xx[j-1],fy-yy[j-1],xx[j],fy-yy[j]);
        x_cur = xx[n-1];
        y_cur = fy-yy[n-1];
        break;

      case opNewPolygon :
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          xx[k] += x_offset;
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] += y_offset;
          data += STD_SHORT;
        }
        if (n<3) break;

        gcv.line_width = lw;
        gcv.foreground = colortable[c].pixel;
        gcv.fill_style = FillSolid;
        XChangeGC( display, gc,
                   GCForeground|GCLineWidth|GCFillStyle, &gcv);

        for (j=0; j<n; j++) {
          xy[j].x = xx[j];
          xy[j].y = fy-yy[j];
        }

        XFillPolygon(display, pixmap, gc,
                     xy, n, Convex, CoordModeOrigin);
        x_cur = xy[n-1].x;
        y_cur = xy[n-1].y;
        break;

      case opNewPolymark :
        m = (short) *((unsigned char *)data);
        data++;
        ms = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          xx[k] += x_offset;
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] += y_offset;
          data += STD_SHORT;
        }
        gcv.line_width = lw;
        gcv.foreground =colortable[c].pixel;
        XChangeGC( display, gc,
                   (GCForeground|GCLineWidth), &gcv);
        for (j=0; j<n; j++) Marker(m,ms,xx[j],fy-yy[j]);
        x_cur = xx[n-1];
        y_cur = fy-yy[n-1];

        break;

      case opNewText :
        ts = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        x += x_offset;
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy - y;
        y += y_offset;
        data += STD_SHORT;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        gcv.foreground =colortable[c].pixel;
        XChangeGC( display, gc,
                   GCForeground, &gcv);
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opNewCenteredText :
        ts = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x,data,STD_SHORT);
        STD2NAT(x,STD_SHORT);
        x += x_offset;
        data += STD_SHORT;
        memcpy(&y,data,STD_SHORT);
        STD2NAT(y,STD_SHORT);
        y = fy - y;
        y += y_offset;
        data += STD_SHORT;
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        gcv.foreground =colortable[c].pixel;
        XChangeGC( display, gc,
                   GCForeground , &gcv);
        w = XTextWidth( font, s, n);
        x_cur = x - w/2;
        y_cur = y + ts/2;
        XDrawString( display, pixmap, gc,
                     x_cur, y_cur, s, n);
        break;

      case opShadedPolygon :
        memcpy(&n,data,STD_SHORT);
        STD2NAT(n,STD_SHORT);
        data += STD_SHORT;
        memcpy(&shd,data,STD_SHORT);
        STD2NAT(shd,STD_SHORT);
        data += STD_SHORT;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n;
        for (k=0; k<n; k++)
        {
          memcpy(xx+k,data,STD_SHORT);
          STD2NAT(xx[k],STD_SHORT);
          xx[k] += x_offset;
          data += STD_SHORT;
        }
        for (k=0; k<n; k++)
        {
          memcpy(yy+k,data,STD_SHORT);
          STD2NAT(yy[k],STD_SHORT);
          yy[k] = fy-yy[k];
          yy[k] += y_offset;
          data += STD_SHORT;
        }

        if (n<3) break;
        for (j=0; j<n; j++) {
          xy[j].x = xx[j];
          xy[j].y = yy[j];
        }
        if (!true_color) {
          snb = (int)(0.5+(float)(NO_PATTERNS-1)*(float)(shd)/1000.0);
          XSetBackground(display ,gc, BlackPixel(display, screen));
          XSetFillStyle(display, gc, FillOpaqueStippled);
          XSetStipple(display, gc, pattern[snb]);
          XFillPolygon(display, pixmap, gc,
                       xy, n, Convex, CoordModeOrigin);
          XSetFillStyle(display, gc, FillSolid);
          XSetBackground(display ,gc, WhitePixel(display, screen));
        }
        else {
          unsigned long pixel;
          float red, green, blue;

          pixel = colortable[ac].pixel;
          red   = (pixel & red_mask  ) >> red_shift;
          green = (pixel & green_mask) >> green_shift;
          blue  = (pixel & blue_mask ) >> blue_shift;
          red   *= shd/1000.0;
          green *= shd/1000.0;
          blue  *= shd/1000.0;
          pixel = ((unsigned long)(red  +0.5) << red_shift  ) +
                  ((unsigned long)(green+0.5) << green_shift) +
                  ((unsigned long)(blue +0.5) << blue_shift );
          XSetForeground(display, gc, pixel);
          XFillPolygon(display, pixmap, gc, xy, n, Convex, CoordModeOrigin);
          XSetForeground(display, gc, colortable[ac].pixel);
        }
        x_cur = xy[n-1].x;
        y_cur = xy[n-1].y;
        break;

      default :
        break;
      }
    }
  }
  fclose(stream);
  return(0);
}



static void Marker (short n,short s,short x,short y)
{
  short top, left, bottom, right;
  XColor color, color2;
  XGCValues gcv;


  XGetGCValues(display, gc, GCForeground, &gcv);
  color.pixel = gcv.foreground;
  XQueryColor( display, cmap, &color);
  color2.pixel = gcv.foreground;
  color2.red = color.red/2;
  color2.green = color.green/2;
  color2.blue = color.blue/2;

  top = y+s; bottom = y-s;
  left = x-s; right = x+s;

  n = n%11;

  switch (n)
  {
  case 0 :
    XDrawRectangle(display, pixmap, gc, left, top, s, s);
    break;
  case 1 :
    XStoreColor(display, cmap, &color2);
    XDrawRectangle(display, pixmap, gc, left, top, s, s);
    XStoreColor(display, cmap, &color);
    break;
  case 2 :
    XDrawRectangle(display, pixmap, gc, left, top, s, s);
    break;
  case 3 :
    XDrawArc( display, pixmap, gc, x, y, 2*s, 2*s, 0, 360*64);
    break;
  case 4 :
    XStoreColor(display, cmap, &color2);
    XDrawArc( display, pixmap, gc, x, y, 2*s, 2*s, 0, 360*64);
    XStoreColor(display, cmap, &color);
    break;
  case 5 :
    XDrawArc( display, pixmap, gc, x, y, 2*s, 2*s, 0, 360*64);
    break;
  case 6 :
    XDrawLine( display, pixmap, gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, pixmap, gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, pixmap, gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, pixmap, gc, x-s/2, y, x, y+s/2);
    break;
  case 7 :
    XStoreColor(display, cmap, &color2);
    XDrawLine( display, pixmap, gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, pixmap, gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, pixmap, gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, pixmap, gc, x-s/2, y, x, y+s/2);
    XStoreColor(display, cmap, &color);
    break;
  case 8 :
    XDrawLine( display, pixmap, gc, x, y+s/2, x+s/2, y);
    XDrawLine( display, pixmap, gc, x+s/2, y, x, y-s/2);
    XDrawLine( display, pixmap, gc, x, y-s/2, x-s/2, y);
    XDrawLine( display, pixmap, gc, x-s/2, y, x, y+s/2);
    break;
  case 9 :

    XDrawLine( display, pixmap, gc, x, y+s/2, x, y-s/2);
    XDrawLine( display, pixmap, gc, x-s/2, y, x+s/2, y);
    break;
  case 10 :
    XDrawLine( display, pixmap, gc, x-s/2, y+s/2, x+s/2, y-s/2);
    XDrawLine( display, pixmap, gc, x-s/2, y-s/2, x+s/2, y+s/2);
    break;
  }
}


void exit_dial (Widget w, XtPointer data, XEvent *event, Boolean *continue_to_dispatch)
{
  if (event->xbutton.button == Button3)
  {
    XtManageChild(xexit);
    XRaiseWindow(display, XtWindow(xexit));
  }
}

void exit_confirm (Widget w, XtPointer clientdata, XtPointer calldata)
{
  XtUnmanageChild(viewport);
  XtDestroyApplicationContext(kontext);
  exit(0);
}

void exit_cancel (Widget w, XtPointer clientdata, XtPointer calldata)
{
  XtUnmanageChild(xexit);
}

void manage_dial (Widget w, XtPointer data, XEvent *event, Boolean *continue_to_dispatch)
{
  if (event->xbutton.button == Button2)
  {
    XtManageChild(dial);
    XRaiseWindow(display, XtWindow(dial));
  }
}


void dialog_confirm(Widget widget, XtPointer clientdata, XtPointer calldata)
{
  XGCValues gcv;
  XPoint edges[4];
  Arg args[10];
  Cardinal n;
  int j;

  XtUnmanageChild(dial);

  file = XawDialogGetValueString(dial);

  if (NULL == (stream = fopen(file, "r")))
  {
    XtManageChild(dial);
    printf("Cannot open file %s\n", file);

    /* initialize dialog widget */
    n = 0;
    XtSetArg(args[n], XtNlabel, "bad filename!\n type new:\n"); n++;
    XtSetArg(args[n], XtNvalue, ".meta"); n++;
    XtSetValues(dial, args, n);

    return;
  }

  /* initialize dialog widget */
  n = 0;
  XtSetArg(args[n], XtNlabel, "new filename:"); n++;
  XtSetArg(args[n], XtNvalue, ".meta"); n++;
  XtSetValues(dial, args, n);

  XClearWindow(display, XtWindow(picture));

  GetFileScreen(stream, &fx, &fy);

  /* change size of drawing window */
  pwidth = (fx < WIN_WIDTH) ? WIN_WIDTH : fx;
  pheight = (fy < WIN_HEIGHT) ? WIN_HEIGHT : fy;

  /* allow Scrollbars */
  n = 0;
  XtSetArg(args[n], XtNallowHoriz, True); n++;
  XtSetArg(args[n], XtNallowVert, True); n++;
  XtSetArg(args[n], XtNwidth, vwidth); n++;
  XtSetArg(args[n], XtNheight, vheight); n++;

  XtSetValues(viewport, args, n);


  /* change size of drawing window */
  XtResizeWidget(picture, pwidth, pheight, 0);


  /* free allocated colorcells */
  for (j=0; j<dispcells; j++)
  {
    if (!used[j])
      XFreeColors(display, def_colormap, &(colortable[j].pixel), 1, 0L);
    else
      used[j] = 0;
    colortable[j].pixel = 0;
  }

  XFreePixmap(display, pixmap);
  pixmap = XCreatePixmap(display, XtWindow(picture), pwidth, pheight,
                         DefaultDepthOfScreen(XtScreen(picture)));

  /* create graphic context */
  gcv.foreground = WhitePixel(display, screen);
  gcv.background = WhitePixel(display, screen);

  gc = XCreateGC( XtDisplay(picture) , XtWindow(picture),
                  GCForeground|GCBackground, &gcv);

  /* clear pixmap's area */
  edges[0].x = 0;
  edges[0].y = 0;
  edges[1].x = 0;
  edges[1].y = pheight;
  edges[2].x = pwidth;
  edges[2].y = pheight;
  edges[3].x = pwidth;
  edges[3].y = 0;

  XFillPolygon(display, pixmap, gc, edges, 4,
               Convex, CoordModeOrigin);

  RasterizeFile(stream);
  fclose(stream);

  XCopyArea(display, pixmap, XtWindow(picture), gc, 0, 0, pwidth, pheight , 0, 0);
}


void dialog_cancel(Widget widget, XtPointer clientdata, XtPointer calldata)
{
  Arg args[10];

  XtUnmanageChild(dial);

  /* initialize dialog widget */
  n = 0;
  XtSetArg(args[n], XtNlabel, "new filename:"); n++;
  XtSetArg(args[n], XtNvalue, ".meta"); n++;
  XtSetValues(dial, args, n);

}


static Boolean run_film (void)
{
  XGCValues gcv;
  XPoint edges[4];
  char command[200];
  short xshift, yshift;


  /* the first frame has already been displayed by main() */

  if (p_opt) {
    sprintf(command, "xwd -name xugv -silent | xwdtopnm >%s_.%04d 2>/dev/null",
            file, frame_number);
    system(command);
  }
  else if (g_opt) {
    sprintf(command, "xwd -name xugv -silent | xwdtopnm 2>/dev/null | ppmtogif >%s_.%04d 2>/dev/null", file, frame_number);
    system(command);
  }

  if (VERBOSE)
  {
    printf("Frame %s done, press <RETURN> to continue",frame);
    getchar();
  }

  if (count) {printf("[%d] ",frame_number); fflush(stdout);}

  /* next frame, please */
  frame_number+=incr;

  if (frame_number>last)
  {
    if (p_opt || g_opt) exit(0);
    if (stoploop) frame_number-=incr;
    else frame_number=first;
    if (count)
    {
      if (stoploop)
      {
        printf("\nFilm %s done, press <RETURN> to start again",file);
        getchar();
        frame_number = first;
      }
      else
        printf("\n");
    }
    run_count++;
    if (run_count==run_max) exit(0);
  }

  for (i_pic=0; i_pic<n_pic; i_pic++)
  {
    sprintf(frame,"%s.%04d",mfile[i_pic],frame_number);
    if (!_wait)
    {
      mstream[i_pic] = fopen(frame,"r");
      if (mstream[i_pic]==NULL) exit(-1);
    }
    else
    {
      mstream[i_pic] = fopen_with_wait(frame,"r");
      break;
    }
  }

  for (i_pic=0; i_pic<n_pic; i_pic++)
  {
    GetFileScreen (mstream[i_pic], mfx+i_pic, mfy+i_pic);
  }

  /* create graphic context */
  gcv.foreground = WhitePixel(display, screen);
  gcv.background = WhitePixel(display, screen);
  gc = XCreateGC( XtDisplay(picture) , XtWindow(picture), GCForeground|GCBackground, &gcv);

  /* clear pixmap's area */
  edges[0].x = 0;
  edges[0].y = 0;
  edges[1].x = 0;
  edges[1].y = pheight;
  edges[2].x = pwidth;
  edges[2].y = pheight;
  edges[3].x = pwidth;
  edges[3].y = 0;

  XFillPolygon(display, pixmap, gc, edges, 4, Convex, CoordModeOrigin);

  /* draw picture into pixmap */
  xshift = 0; yshift = -fy_max;
  for (i_pic=0; i_pic<n_pic; i_pic++)
  {
    if (i_pic%nbreak == 0)
    {
      xshift = 0;
      yshift += fy_max;
    }
    else
      xshift += fx_max;

    stream = mstream[i_pic];
    RasterizePositionedFile(stream,xshift,yshift);
    fclose(stream);
  }

  XCopyArea(display, pixmap, XtWindow(picture), gc, 0, 0, pwidth, pheight , 0, 0);

  return(False);
}

FILE *auto_fopen (char *name)
{
  int i;
  char name_ext[1024],*home,lock[256];
  FILE *f;

  home=getenv("HOME"); sprintf(lock,"%s/.xugv_tail",home);
  if (auto_nb==-2)
  {
    /* initialize 'auto_nb' */
    for (i=0; i<10000; i++)
    {
      sprintf(name_ext,"%s.%0.4d",name,i);
      f=fopen(name_ext,"r");
      if (f==NULL) { i-=2; break; }
      fclose(f);
    }
    auto_nb=i; if (auto_nb<-1) auto_nb=-1;

    /* open metafile */
    return (auto_fopen(name));
  }
  else
  {
    /* open metafile with highest number > auto_nb existing */
    auto_nb++;
    while(1)             /* waiting until auto_nb+1 exists */
    {
      sprintf(name_ext,"%s.%0.4d",name,auto_nb);
      f=fopen(name_ext,"r");
      if (f==NULL) { sleep(sleep_seconds); continue; }
      fclose(f); break;
    }
    while(1)             /* check if higher files exist */
    {
      sprintf(name_ext,"%s.%0.4d",name,auto_nb+1);
      f=fopen(name_ext,"r");
      if (f!=NULL) { fclose(f); auto_nb++; continue; }
      break;
    }
    if (auto_nb>9999) exit(0);

    while(1)             /* open file auto_nb when '.xugv_tail' pops up */
    {
      sleep(sleep_seconds);
      f=fopen(lock,"r"); if (f==NULL) continue;
      sprintf(lock,"rm %s/.xugv_tail",home); system(lock);
      sprintf(name_ext,"%s.%0.4d",name,auto_nb);
      f=fopen(name_ext,"r");
      if (f!=NULL)
      {
        printf("xugv: displaying '%s'\n",name_ext);
        return (f);
      }
    }
  }
}

static Boolean tail_film (void)
{
  XGCValues gcv;
  XPoint edges[4];
  short xshift, yshift;

  /* try next meta-file */
  mstream[0]=auto_fopen(file);
  GetFileScreen (mstream[0],mfx,mfy);

  /* create graphic context */
  gcv.foreground = WhitePixel(display, screen);
  gcv.background = WhitePixel(display, screen);
  gc = XCreateGC( XtDisplay(picture) , XtWindow(picture), GCForeground|GCBackground, &gcv);

  /* clear pixmap's area */
  edges[0].x = 0;
  edges[0].y = 0;
  edges[1].x = 0;
  edges[1].y = pheight;
  edges[2].x = pwidth;
  edges[2].y = pheight;
  edges[3].x = pwidth;
  edges[3].y = 0;

  XFillPolygon(display, pixmap, gc, edges, 4, Convex, CoordModeOrigin);

  /* draw picture into pixmap */
  xshift = 0; yshift = -fy_max;
  for (i_pic=0; i_pic<n_pic; i_pic++)
  {
    if (i_pic%nbreak == 0)
    {
      xshift = 0;
      yshift += fy_max;
    }
    else
      xshift += fx_max;

    stream = mstream[i_pic];
    RasterizePositionedFile(stream,xshift,yshift);
    fclose(stream);
  }

  XCopyArea(display, pixmap, XtWindow(picture), gc, 0, 0, pwidth, pheight , 0, 0);

  return(False);
}

/* main */
int main (int argc, char* argv[])
{
  Arg args[10];           /* argument list */
  short xshift, yshift;
  int i, nopt;
  XtWorkProcId film_work_id, tail_work_id;

  /* check for little/big endian storage type */
  littleEndian = !(*((char *) &littleEndian));
  if (VERBOSE) printf("littleEndian=%d",littleEndian);

  /* initialization and creating of application shell */

  applShell = XtAppInitialize (&kontext, APPL_KLASSE,
                               (XrmOptionDescRec*)NULL, 0,
                               &argc, argv,
                               (String*)NULL,
                               (Arg*)NULL, 0);

  /* initialize display variable */
  display = XtDisplay (applShell);

  /* enough colorcells */
  dispcells = XUGV_MIN(DisplayCells(display, DefaultScreen(display)), CSIZE);
  /* printf("%d color cells available\n",dispcells); */
  if (dispcells <= 0)
  {
    printf("%s: This program requires a color or greyscale monitor\n", argv[0]);
    exit(-1);
  }

  if (argc < 2) {
    printf("usage: xugv [<nb of files>] file [file2] [file3] ... [-v[n]] [-f first last] [-q increment] [-p|-g] [-c] [-s] [-N <nBreak>] [-n]\n");
    exit(-1);
  }

  /* for (i=0; i<argc; i++) printf("%d %s\n",i,argv[i]); */

  if (sscanf(argv[1],"%d",&n_pic)!=1)
  {
    i = 2;
    n_pic = 1;
    f_offset = 1;
  }
  else
  {
    i = 2+n_pic;
    f_offset = 2;
  }
  if (n_pic == 1) file = argv[1];else file = argv[2];
  film=0;
  option = "";
  count = 0;
  stoploop = 0;
  nopt = 0;
  while (i<argc)
  {
    if (argv[i][1]=='v') {option = argv[i]; i++; continue;}
    if (argv[i][1]=='V') {verbose = 1; i++; continue;}
    if (argv[i][1]=='c') {count = 1; i++; continue;}
    if (argv[i][1]=='w') {_wait = 1; i++; continue;}
    if (argv[i][1]=='s') {stoploop = 1; i++; continue;}
    if (argv[i][1]=='f')
    {
      if (i+2>=argc) {
        printf("not enough arguments for film option\n");
        exit(-1);
      }
      sscanf(argv[i+1],"%d",&first);
      sscanf(argv[i+2],"%d",&last);
      frame_number = first;
      film = 1;
      i+=3;
      continue;
    }
    if (argv[i][1]=='q')
    {
      if (i+1>=argc) {
        printf("not enough arguments for increment option\n");
        exit(-1);
      }
      sscanf(argv[i+1],"%d",&incr);
      i+=2;
      continue;
    }
    if (argv[i][1]=='p')
    {
      p_opt=1;
      i++;
      continue;
    }
    if (argv[i][1]=='g')
    {
      g_opt=1;
      i++;
      continue;
    }
    if (argv[i][1]=='N')
    {
      if (i+1>=argc) {
        printf("not enough arguments for N option\n");
        exit(-1);
      }
      sscanf(argv[i+1],"%d",&nbreak);
      if (nbreak<1) {
        printf("N must be > 0\n");
        exit(-1);
      }
      nopt=1;
      i+=2;
      continue;
    }

    if (argv[i][1]=='y')
    {
      if (i+1>=argc) {
        printf("not enough arguments for n option\n");
        exit(-1);
      }
      run_max = 0;
      sscanf(argv[i+1],"%d",&run_max);
      i+=2;
      continue;
    }
    i++;
  }
  if (!nopt) nbreak=n_pic;

  if (!film)
  {
    for (i_pic=0; i_pic<n_pic; i_pic++ )
    {
      mstream[i_pic]=fopen(argv[f_offset+i_pic],"r");
      if (mstream[i_pic]==NULL)
      {
        if (i_pic==0 && n_pic==1)
        {
          auto_nb=-2;
          mstream[0]=auto_fopen(argv[f_offset]);
          if (mstream[0]==NULL) { printf("Can't open file %s\n", file); exit(-1); }
        }
        else { printf("cannot open %d files\n", n_pic); exit(-1); }
      }
    }

    fx_max = 0; fy_max = 0;
    for (i_pic=0; i_pic<n_pic; i_pic++)
    {
      GetFileScreen (mstream[i_pic], mfx+i_pic, mfy+i_pic);
      fx_max = (fx_max > mfx[i_pic]) ? (fx_max) : (mfx[i_pic]);
      fy_max = (fy_max > mfy[i_pic]) ? (fy_max) : (mfy[i_pic]);
    }
    pwidth = nbreak * fx_max;
    pheight = ceil((float)n_pic/(float)nbreak) * fy_max;

    /* initialize viewport window */
    n = 0;
    XtSetArg(args[n], XtNwidth,  XUGV_MIN(pwidth,  WIN_WIDTH )); n++;
    XtSetArg(args[n], XtNheight, XUGV_MIN(pheight, WIN_HEIGHT)); n++;
    XtSetArg(args[n], XtNallowHoriz, True); n++;
    XtSetArg(args[n], XtNallowVert, True);  n++;

    /* create viewport widget */
    viewport = XtCreateManagedWidget ("viewport", viewportWidgetClass,applShell,args,n);

    /* set size of drawing widget */
    n = 0;
    XtSetArg(args[n], XtNwidth, pwidth); n++;
    XtSetArg(args[n], XtNheight, pheight); n++;

    /* create drawing widget */
    picture = XtCreateManagedWidget ("picture", simpleWidgetClass,viewport,args,n);

    /* establish exit dialog */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Do you really want\n to exit xugv?"); n++;
    XtSetArg(args[n], XtNx, 25); n++;
    XtSetArg(args[n], XtNy, 25); n++;
    xexit = XtCreateWidget ("exit", dialogWidgetClass,viewport,args,n);

    /* add some nice buttons */
    XawDialogAddButton(xexit, "confirm", exit_confirm, (XtPointer)NULL);
    XawDialogAddButton(xexit, "cancel", exit_cancel, (XtPointer)NULL);

    /* add event handler for expose events */
    XtAddEventHandler( picture, ExposureMask, FALSE, (XtEventHandler) exposeCB, (XtPointer)NULL);
    XtAddEventHandler(viewport, ButtonPressMask, FALSE, (XtEventHandler) exit_dial, (XtPointer)NULL);

    /* install work if film is tailed */
    if (auto_nb>=0)
    {
      tail_work_id = XtAppAddWorkProc(kontext,(XtWorkProc)tail_film,applShell);
    }

    /* realize widget tree */
    XtRealizeWidget (applShell);

    /* create colormap and graphic context */
    createGraphics();

    /* draw picture into pixmap */
    xshift = 0; yshift = -fy_max;
    for (i_pic=0; i_pic<n_pic; i_pic++)
    {
      if (i_pic%nbreak == 0)
      {
        xshift = 0;
        yshift += fy_max;
      }
      else
        xshift += fx_max;

      stream = mstream[i_pic];
      RasterizePositionedFile(stream,xshift,yshift);
      fclose(stream);
      ignore = 1;
    }

    /* main loop for event processing */
    XtAppMainLoop (kontext);
  }
  else
  {
    i_pic=0;
    while (1)
    {
      sprintf(frame,"%s.%04d",argv[f_offset+i_pic],first);
      mfile[i_pic] = argv[f_offset+i_pic];
      if (!_wait)
      {
        mstream[i_pic] = fopen(frame,"r");
        if (mstream[i_pic]==NULL) break;
        i_pic++;
      }
      else
      {
        mstream[i_pic] = fopen_with_wait(frame,"r");
        i_pic++;
        break;
      }
    }
    if (i_pic<1)
    {
      printf("Can't open file %s\n", file);
      exit(-1);
    }
    if (n_pic != i_pic)
    {
      printf("cannot open %d files\n", n_pic);
      exit(-1);
    }
    if (i_pic>1)
    {
      printf("opened %d files\n", n_pic);
    }

    fx_max = 0; fy_max = 0;
    for (i_pic=0; i_pic<n_pic; i_pic++)
    {
      GetFileScreen (mstream[i_pic], mfx+i_pic, mfy+i_pic);
      fx_max = (fx_max > mfx[i_pic]) ? (fx_max) : (mfx[i_pic]);
      fy_max = (fy_max > mfy[i_pic]) ? (fy_max) : (mfy[i_pic]);
    }
    pwidth = nbreak * fx_max;
    pheight = ceil((float)n_pic/(float)nbreak) * fy_max;

    /* initialize viewport window */
    n = 0;
    XtSetArg(args[n], XtNwidth,  XUGV_MIN(pwidth,  WIN_WIDTH )); n++;
    XtSetArg(args[n], XtNheight, XUGV_MIN(pheight, WIN_HEIGHT)); n++;
    XtSetArg(args[n], XtNallowHoriz, True); n++;
    XtSetArg(args[n], XtNallowVert, True); n++;


    /* create viewport widget */
    viewport = XtCreateManagedWidget ("viewport", viewportWidgetClass,
                                      applShell, args, n);

    /* set size of drawing widget */
    n = 0;
    XtSetArg(args[n], XtNwidth, pwidth); n++;
    XtSetArg(args[n], XtNheight, pheight); n++;


    /* create drawing widget */
    picture = XtCreateManagedWidget ("picture", simpleWidgetClass,
                                     viewport, args, n);

    /* establish exit dialog */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Do you really want\n to exit xugv?"); n++;
    XtSetArg(args[n], XtNx, 25); n++;
    XtSetArg(args[n], XtNy, 25); n++;

    xexit = XtCreateWidget ("exit", dialogWidgetClass,
                            viewport, args, n);

    /* add some nice buttons */
    XawDialogAddButton(xexit, "confirm", exit_confirm, (XtPointer)NULL);
    XawDialogAddButton(xexit, "cancel", exit_cancel, (XtPointer)NULL);

    /* add event handler for expose events */
    XtAddEventHandler( picture, ExposureMask, FALSE, (XtEventHandler) exposeCB, (XtPointer)NULL);

    /* add event handler for button press event */
    XtAddEventHandler(viewport, ButtonPressMask, FALSE, (XtEventHandler) exit_dial, (XtPointer)NULL);

    /* register a work procedure to display the film */
    film_work_id = XtAppAddWorkProc(kontext,(XtWorkProc) run_film,applShell);

    /* realize widget tree */
    XtRealizeWidget (applShell);

    /* use backing store for window
       {
            XSetWindowAttributes attr;
            unsigned long mask;

       attr.backing_store = Always;
       mask = CWBackingStore;
       XChangeWindowAttributes(display, XtWindow(picture), mask, &attr);
       } */

    /* create colormap and graphic context */
    createGraphics();

    /* draw picture into pixmap */
    xshift = 0; yshift = -fy_max;
    for (i_pic=0; i_pic<n_pic; i_pic++)
    {
      if (i_pic%nbreak == 0)
      {
        xshift = 0;
        yshift += fy_max;
      }
      else
        xshift += fx_max;

      stream = mstream[i_pic];
      RasterizePositionedFile(stream,xshift,yshift);
      fclose(stream);
      ignore = 1;
    }

    /* ignore color maps in subsequent frames */
    ignore = 1;

    /* main loop for event processing */
    run_count = 0;
    XtAppMainLoop (kontext);
  }
}
