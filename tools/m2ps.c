// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      m2ps.c                                                        */
/*                                                                          */
/* Purpose:   converts ug metafiles to postscript                           */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr.uni-heidelberg.de                       */
/*                                                                          */
/* History:   21 Mar 94 begin                                               */
/*                                                                          */
/* Remarks:                                                                 */
/*            This code uses some fragments from Juergen Bey's postscript   */
/*            output routines from ug3.                                     */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/


#include <string.h>
#include <strings.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define CREATOR "ug metafile converter"


/* postscript defaults */
#define  PS_DEFAULT_LINE_WIDTH     1
#define  PS_DEFAULT_FONT           "Helvetica"
#define  PS_FONT_FACTOR            10
#define  LW_FACTOR                 0.03
#define  LW_SCALE                  50.0
#define  PAPER_X                   600
#define  PAPER_Y                   850

/* metafile op codes for graphics primitives */
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

#define SIZE                    50
#define CSIZE                   256
#define COLORS          256

#define TRFMX(x,y) (((float)(x))*mxx + ((float)(y))*mxy + tx)
#define TRFMY(x,y) (((float)(x))*myx + ((float)(y))*myy + ty)

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/


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

static FILE *psfile;
static short IF_cx=0, IF_cy=0;
static float red[COLORS],green[COLORS],blue[COLORS];
static short IF_lw=-1, IF_ts=-1;
static unsigned char IF_cc=0;
static int ox,oy,sx,sy;
static float tx,ty,mxx,myy,mxy,myx;
static int landscape = 0;

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

void IF_LineWidth (short n);
void IF_TextSize (short n);


/***************************************************************************
*                                                                         *
* Function : OpenPSFile                                                   *
*                                                                         *
* Purpose  : Open postscript file                                         *
*                                                                         *
* Input    : STRING  filename   - name of the postscript file to create   *
*            STRING  header     - name of the file containing the header  *
*                                 to copy to the created file             *
*                                                                         *
* Output   : FILE * - the file handle                                     *
*                                                                         *
***************************************************************************/

FILE  *OpenPSFile(char *filename, char *title, int x, int y, int w, int h )
{
  FILE   *file;
  time_t tp;
  char date[64];

  if( ( file = fopen(filename, "w") ) == NULL ) return(NULL);

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


  return(file);
}

/***************************************************************************
*                                                                         *
* Function : ClosePSFile                                                  *
*                                                                         *
* Purpose  : Close postscript file                                        *
*                                                                         *
* Input    : FILE    *file      - handle of postscript file               *
*                                                                         *
***************************************************************************/

int ClosePSFile(FILE *file)
{
  if (file==NULL) return(0);

  /* write trailer */

  fprintf(file,"\nshowpage\n\n");
  fprintf(file,"%%%%Trailer\n");

  fclose(file);
  return(1);
}


/****************************************************************************/
/*                                                                          */
/* graphics primitives                                                      */
/*                                                                          */
/****************************************************************************/

void IF_MoveTo (short x, short y)
{
  IF_cx = x; IF_cy = y;
  return;
}

void IF_DrawTo (short x, short y)
{
  fprintf(psfile,"%g %g M %g %g S\n",
          TRFMX(IF_cx,IF_cy),TRFMY(IF_cx,IF_cy),TRFMX(x,y),TRFMY(x,y));
  IF_cx = x; IF_cy = y;
  return;
}

void IF_PolyLine (short n, short *x, short *y)
{
  int i;

  fprintf(psfile,"N\n");
  fprintf(psfile,"%g %g M\n",TRFMX(x[0],y[0]),TRFMY(x[0],y[0]));
  for (i=1; i<n; i++)
    fprintf(psfile,"%g %g L\n",TRFMX(x[i],y[i]),TRFMY(x[i],y[i]));
  fprintf(psfile,"stroke\n");

  return;
}

void IF_Polygon (short n, short *x, short *y)
{
  int i;

  fprintf(psfile,"N\n");
  fprintf(psfile,"%g %g M\n",TRFMX(x[0],y[0]),TRFMY(x[0],y[0]));
  for (i=1; i<n; i++)
    fprintf(psfile,"%g %g L\n",TRFMX(x[i],y[i]),TRFMY(x[i],y[i]));
  fprintf(psfile,"C\n");

  return;
}

static void IF_Marker (short n, short s, short x, short y)
{
  short top, left, bottom, right;
  short r,g,b;
  short xx[10],yy[10];

  s = s/2;

  top = y+s; bottom = y-s;
  left = x-s; right = x+s;

  n = n%12;

  switch (n)
  {
  case 0 :
    xx[0] = left; yy[0] = bottom;
    xx[1] = right; yy[1] = bottom;
    xx[2] = right; yy[2] = top;
    xx[3] = left; yy[3] = top;
    xx[4] = left; yy[4] = bottom;
    IF_PolyLine(5,xx,yy);
    break;
  case 1 :
    xx[0] = left; yy[0] = bottom;
    xx[1] = right; yy[1] = bottom;
    xx[2] = x; yy[2] = top;
    xx[3] = left; yy[3] = bottom;
    IF_PolyLine(4,xx,yy);
    break;
  case 2 :
    xx[0] = x; yy[0] = bottom;
    xx[1] = right; yy[1] = top;
    xx[2] = left; yy[2] = top;
    xx[3] = x; yy[3] = bottom;
    IF_PolyLine(4,xx,yy);
    break;
  case 3 :
    xx[0] = x; yy[0] = bottom;
    xx[1] = right; yy[1] = y;
    xx[2] = x; yy[2] = top;
    xx[3] = left; yy[3] = y;
    xx[4] = x; yy[4] = bottom;
    IF_PolyLine(5,xx,yy);
    break;
  case 4 :
    xx[0] = left; yy[0] = bottom;
    xx[1] = right; yy[1] = y;
    xx[2] = left; yy[2] = top;
    xx[3] = left; yy[3] = bottom;
    IF_PolyLine(4,xx,yy);
    break;
  case 5 :
    xx[0] = right; yy[0] = bottom;
    xx[1] = right; yy[1] = top;
    xx[2] = left; yy[2] = y;
    xx[3] = right; yy[3] = bottom;
    IF_PolyLine(4,xx,yy);
    break;
  case 6 :
    xx[0] = left; yy[0] = bottom;
    xx[1] = right; yy[1] = bottom;
    xx[2] = right; yy[2] = top;
    xx[3] = left; yy[3] = top;
    IF_Polygon(4,xx,yy);
    break;
  case 7 :
    xx[0] = left; yy[0] = bottom;
    xx[1] = right; yy[1] = bottom;
    xx[2] = x; yy[2] = top;
    IF_Polygon(3,xx,yy);
    break;
  case 8 :
    xx[0] = x; yy[0] = bottom;
    xx[1] = right; yy[1] = top;
    xx[2] = left; yy[2] = top;
    IF_Polygon(3,xx,yy);
    break;
  case 9 :
    xx[0] = x; yy[0] = bottom;
    xx[1] = right; yy[1] = y;
    xx[2] = x; yy[2] = top;
    xx[3] = left; yy[3] = y;
    IF_Polygon(4,xx,yy);
    break;
  case 10 :
    xx[0] = left; yy[0] = bottom;
    xx[1] = right; yy[1] = y;
    xx[2] = left; yy[2] = top;
    IF_Polygon(3,xx,yy);
    break;
  case 11 :
    xx[0] = right; yy[0] = bottom;
    xx[1] = right; yy[1] = top;
    xx[2] = left; yy[2] = y;
    IF_Polygon(3,xx,yy);
    break;
  }
}

static void PrintPSString(char *s)
{
  fputc('(',psfile);

  while( *s != '\0' )
  {
    switch( *s )
    {
    case '\\' : case '(' : case ')' : fputc('\\',psfile);
    default  : fputc(*s,psfile);
    }
    s++;
  }
  fputc(')',psfile);
}

void IF_Text (char *s)
{
  fprintf(psfile,"%g %g M\n",TRFMX(IF_cx,IF_cy),TRFMY(IF_cx,IF_cy));
  if (landscape) fprintf(psfile,"90 rotate\n");
  PrintPSString(s);
  fprintf(psfile,"show N\n");
  if (landscape) fprintf(psfile,"-90 rotate\n");
  return;
}

void IF_CenteredText (short x, short y, char *s)
{
  IF_MoveTo(x,y);
  IF_Text(s);
}

void IF_LineWidth (short n)
{
  if (n<1) n=1;
  if (n==IF_lw) return;

  fprintf(psfile,"%.3f W\n",LW_FACTOR+((float)(n-1))*LW_SCALE*LW_FACTOR);
  IF_lw = n;
  return;
}


void IF_TextSize (short n)
{
  if (n==IF_ts) return;

  fprintf(psfile,"/%s findfont %d scalefont setfont\n",PS_DEFAULT_FONT,(int)n);
  IF_ts = n;
  return;
}

void IF_PrintColor (float color)
{
  if (color==0.0) fprintf(psfile,"%d ",0);
  else if (color==1.0) fprintf(psfile,"%d ",1);
  else fprintf(psfile,"%.3f ",color);
  return;
}

void IF_Foreground (unsigned char n)
{
  if (IF_cc==n) return;

  IF_PrintColor(red[n]);
  IF_PrintColor(green[n]);
  IF_PrintColor(blue[n]);
  fprintf(psfile,"R\n");
  IF_cc = n;
  return;
}

void IF_SetEntry (unsigned char n, short r, short g, short b)
{
  red[n] = ((float)r)/255.0;
  green[n] = ((float)g)/255.0;
  blue[n] = ((float)b)/255.0;

  IF_PrintColor(red[n]);
  IF_PrintColor(green[n]);
  IF_PrintColor(blue[n]);
  fprintf(psfile,"R\n");
  IF_cc = n;
  return;
}

void IF_SetPalette (short x, short y, short *r, short *g, short *b)
{
  int i;

  for (i=x; i<=y; i++)
  {
    red[i] = ((float)r[i-x])/255.0;
    green[i] = ((float)g[i-x])/255.0;
    blue[i] = ((float)b[i-x])/255.0;
  }
  IF_PrintColor(red[x]);
  IF_PrintColor(green[x]);
  IF_PrintColor(blue[x]);
  fprintf(psfile,"R\n");
  IF_cc = (unsigned char) x;
  return;
}

/****************************************************************************/
/*                                                                          */
/* Function:  ConvertFile                                                   */
/*                                                                          */
/* Purpose:   parse an opened metafile                                      */
/*                                                                          */
/* Input:     FILE *stream open metafile                                    */
/*                                                                          */
/* Output:    0 OK, else error                                              */
/*                                                                          */
/****************************************************************************/

static void swap_short (void *data)
{
  char *s,c;

  s = (char *) data;

  c = s[0];
  s[0] = s[1];
  s[1] = c;
}

static void swap_long (void *data)
{
  char *s,c;

  s = (char *) data;

  c = s[0];
  s[0] = s[3];
  s[3] = c;
  c = s[1];
  s[1] = s[2];
  s[2] = c;
}

int ConvertFile (FILE *stream)
{
  char *buffer;                                                 /* input buffer						*/
  long blockSize;                                               /* METABUFFERSIZE					*/
  long blockUsed;                                               /* actual buffer size used			*/
  long itemCounter;                                             /* number of commands in buffer		*/
  char *data;                                                           /* data pointer in buffer			*/
  short fx,fy;                                                  /* file screen size					*/
  int i,error,j,size;
  char opCode;
  short x,y,r,g,b,n,lw,ts,m,ms,w,x1,y1,x2,y2;
  short xx[SIZE],yy[SIZE];
  short rr[COLORS],gg[COLORS],bb[COLORS];
  char s[CSIZE];
  unsigned char c;
  long l;
  short vector[2];

  /* get file parameters */
  rewind(stream);
  error = fread(&blockSize,4,1,stream); if (error!=1) return(1);       /* block size */
  error = fread(&fx,2,1,stream);            if (error!=1) return(1);       /* x size */
  error = fread(&fy,2,1,stream);        if (error!=1) return(1);       /* y size */
#ifdef __SWAPBYTES__
  swap_long(&blockSize);
  swap_short(&fx);
  swap_short(&fy);
#endif

  /* default values */
  lw = 1;
  ts = 12;
  m = 0;
  ms = 6;

  /* allocate input buffer */
  buffer = malloc(blockSize);
  if (buffer==NULL) return(1);

  /* loop through the blocks */
  while (!feof(stream))
  {
    /* read block parameters */
    error = fread(&blockUsed,4,1,stream);    if (error!=1) {free(buffer); return(1);}
    error = fread(&itemCounter,4,1,stream);  if (error!=1) {free(buffer); return(1);}
    error = fread(buffer,blockUsed,1,stream);if (error!=1) {free(buffer); return(1);}
#ifdef __SWAPBYTES__
    swap_long(&blockUsed);
    swap_long(&itemCounter);
#endif

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
        memcpy(&x,data,2);
        data += 2;
        memcpy(&y,data,2);
        data += 2;
#ifdef __SWAPBYTES__
        swap_short(&x);
        swap_short(&y);
#endif
        IF_MoveTo(x,y);
        break;

      case opDraw :
        memcpy(&x,data,2);
        data += 2;
        memcpy(&y,data,2);
        data += 2;
        IF_DrawTo(x,y);
        break;

      case opPolyline :
        memcpy(&n,data,2);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        IF_PolyLine(n,xx,yy);
        break;

      case opPolygon :
        memcpy(&n,data,2);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        IF_Polygon(n,xx,yy);
        break;

      case opPolymark :
        memcpy(&n,data,2);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        for (j=0; j<n; j++) IF_Marker(m,ms,xx[j],yy[j]);
        break;

      case opText :
        memcpy(&n,data,2);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        IF_Text(s);
        break;

      case opCenteredText :
        memcpy(&x,data,2);
        data += 2;
        memcpy(&y,data,2);
        data += 2;
        memcpy(&n,data,2);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        IF_CenteredText(x,y,s);
        break;

      case opSetLineWidth :
        memcpy(&n,data,2);
        data += 2;
        lw = n;
        IF_LineWidth(n);
        break;

      case opSetTextSize :
        memcpy(&n,data,2);
        data += 2;
        ts = n;
        IF_TextSize(n);
        break;

      case opSetMarker :
        memcpy(&n,data,2);
        data += 2;
        m = n;
        break;

      case opSetMarkerSize :
        memcpy(&n,data,2);
        data += 2;
        ms = n;
        break;

      case opSetColor :
        c = *((unsigned char *)data);
        data++;
        IF_Foreground(c);
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
        IF_SetEntry(c,r,g,b);
        break;

      case opSetPalette :
        x = (short) (*((unsigned char *)data));
        data++;
        y = (short) (*((unsigned char *)data));
        data++;
        if ((x<0)||(y>=COLORS)) break;
        for (j=x; j<=y; j++)
        {
          r = (short) (*((unsigned char *)data));
          data++;
          g = (short) (*((unsigned char *)data));
          data++;
          b = (short) (*((unsigned char *)data));
          data++;
          rr[j-x] = r;
          gg[j-x] = g;
          bb[j-x] = b;
        }
        IF_SetPalette(x,y,rr,gg,bb);
        break;

      case opNewLine :
        lw = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x1,data,2);
        data += 2;
        memcpy(&y1,data,2);
        data += 2;
        memcpy(&x2,data,2);
        data += 2;
        memcpy(&y2,data,2);
        data += 2;
        IF_LineWidth(lw);
        IF_Foreground(c);
        IF_MoveTo(x1,y1);
        IF_DrawTo(x2,y2);
        break;

      case opNewPolyline :
        lw = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,2);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        IF_LineWidth(lw);
        IF_Foreground(c);
        IF_PolyLine(n,xx,yy);
        break;

      case opNewPolygon :
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,2);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        if (n<3) break;
        IF_Foreground(c);
        IF_Polygon(n,xx,yy);
        break;

      case opNewPolymark :
        m = (short) *((unsigned char *)data);
        data++;
        ms = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&n,data,2);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        IF_Foreground(c);
        for (j=0; j<n; j++) IF_Marker(m,ms,xx[j],yy[j]);
        break;

      case opNewText :
        ts = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x,data,2);
        data += 2;
        memcpy(&y,data,2);
        data += 2;
        memcpy(&n,data,2);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        IF_TextSize(ts);
        IF_Foreground(c);
        IF_MoveTo(x,y);
        IF_Text(s);
        break;

      case opNewCenteredText :
        ts = (short) *((unsigned char *)data);
        data++;
        c = *((unsigned char *)data);
        data++;
        memcpy(&x,data,2);
        data += 2;
        memcpy(&y,data,2);
        data += 2;
        memcpy(&n,data,2);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        IF_TextSize(ts);
        IF_Foreground(c);
        IF_CenteredText(x,y,s);
        break;

      default :
        break;
      }
    }
  }

  return(0);
}

int GetFileScreen (FILE *stream, short *fx, short *fy)
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


int main (int argc, char **argv)
{
  int i,j,error;
  int sizeopt,origopt,outopt;
  FILE *stream;
  short fx,fy;
  char outname[80];

  if (argc==1)
  {
    printf("usage: m2ps <infile> [-o <outfile>] [-z x y] [-s x y] [-L]\n");
    return(0);
  }

  /* scan options */
  sizeopt = 0;
  origopt = 0;
  outopt = 0;
  i = 2;
  while (i<argc)
  {
    if (argv[i][1] == 'L') landscape=1;
    if (argv[i][1] == 'o')
    {
      if (i+1>=argc)
      {
        printf("not enough arguments for -o option\n");
        return(0);
      }
      strcpy(outname,argv[i+1]);
      outopt = 1;
      i += 1;
    }
    if (argv[i][1] == 'z')
    {
      if (i+2>=argc)
      {
        printf("not enough arguments for -z option\n");
        return(0);
      }
      sscanf(argv[i+1],"%d",&ox);
      sscanf(argv[i+2],"%d",&oy);
      origopt = 1;
      i += 2;
    }
    if (argv[i][1] == 's')
    {
      if (i+2>=argc)
      {
        printf("not enough arguments for -s option\n");
        return(0);
      }
      sscanf(argv[i+1],"%d",&sx);
      sscanf(argv[i+2],"%d",&sy);
      sizeopt = 1;
      i += 2;
    }
    i++;
  }

  stream = fopen(argv[1],"rb");
  if (stream==NULL)
  {
    printf("Could not open file %s\n", argv[1]);
    return(0);
  }

  GetFileScreen(stream,&fx,&fy);
  if (!sizeopt) { sx=fx; sy=fy; }
  if (!origopt) { ox=0;  oy=0;  }
  if (!outopt) sprintf(outname,"%s.ps",argv[1]);

  /* compute transformation */
  if (!landscape)
  {
    tx = (float) ox;
    ty = (float) oy;
    mxx = ((float)(sx))/((float)fx);
    myy = ((float)(sy))/((float)fy);
    mxy = myx = 0.0;
  }
  else
  {
    tx = (float) (PAPER_X-ox);
    ty = (float) oy;
    mxy = -((float)(sy))/((float)fy);
    myx = ((float)(sx))/((float)fx);
    mxx = myy = 0.0;
    ox = PAPER_X-ox-sy;
    j = sx; sx = sy;
    sy = j;
  }

  /* open post script file */
  psfile = OpenPSFile(outname,argv[1],ox,oy,sx+1,sy+1);
  if (psfile==NULL) exit(1);
  IF_LineWidth(1);
  IF_TextSize(10);

  /* convert file */
  ConvertFile(stream);
  fclose(stream);

  /* close postscript file */
  ClosePSFile(psfile);

  exit(0);
}
