// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugView.c                                                      */
/*                                                                          */
/* Purpose:   ug metafile display program for the IRIS Indigo               */
/*                                                                          */
/* Author:    Christoph Kindl, Peter Bastian                                */
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

#include <gl/gl.h>
#include <gl/device.h>

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

#define opNop                   0
#define opMove                  1
#define opDraw                  2
#define opPolyline              3
#define opPolygon               4
#define opPolymark              5
#define opText                  6
#define opCenteredText          7
#define opSetLineWidth          8
#define opSetMarker             9
#define opSetMarkerSize         10
#define opSetTextSize           11
#define opSetColor              12
#define opSetEntry              13
#define opSetPalette            14
#define opNewLine               15
#define opNewPolyline           16
#define opNewPolygon            17
#define opNewPolymark           18
#define opNewText               19
#define opNewCenteredText       20
#define opShadedPolygon     21

#define SIZE                    50
#define CSIZE                   256

#define TRFMX(x) l=x; l=l<<16; x = HiWrd(FixMul(l,sx))          /* use only in Rasterize ! */
#define TRFMY(y) l=y; l=l<<16; y = wy-1-HiWrd(FixMul(l,sy))

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

typedef struct {
  int winid;                    /* window identifier(gl) */
  short red[256];               /* color table	*/
  short green[256];
  short blue[256];
  int sx,sy;
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

static AWindowRecord myWindow;
static short cc;

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/


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

int CreateApplicationWindow (AWindowRecord *wr, char *name, short h, short v, short dh, short dv)
{
  /* read in resources */
  prefposition(h,h+dh-1,v,v+dv-1);
  if((wr->winid=winopen(name))==NULL) { fprintf(stderr,"Could not open window\n"); return(1); }
  winset(wr->winid);
  RGBmode();
  doublebuffer();
  gconfig();
  cpack(0xffffff);
  clear();
  swapbuffers();
  cpack(0xffffff);
  clear();
  wr->sx =dh; wr->sy = dv;

  return(0);
}

void DisposeApplicationWindow (AWindowRecord *wr)
{
  winclose(wr->winid);
  return;
}


static void Marker (short n, short s, short x, short y)
{
  short top, left, bottom, right;
  short r,g,b;
  short vector[2];

  gRGBcolor(&r,&g,&b);
  s = s/2;

  top = y+s; bottom = y-s;
  left = x-s; right = x+s;

  n = n%11;

  switch (n)
  {
  case 0 :
    rects(left,bottom,right,top);
    break;
  case 1 :
    RGBcolor(r/2,g/2,b/2);
    rectfs(left,bottom,right,top);
    RGBcolor(r,g,b);
    break;
  case 2 :
    rectfs(left,bottom,right,top);
    break;
  case 3 :
    circs(x,y,s);
    break;
  case 4 :
    RGBcolor(r/2,g/2,b/2);
    circfs(x,y,s);
    RGBcolor(r,g,b);
    break;
  case 5 :
    circfs(x,y,s);
    break;
  case 6 :
    move2s(x,y+s);
    draw2s(x+s,y);
    draw2s(x,y-s);
    draw2s(x-s,y);
    draw2s(x,y+s);
    break;
  case 7 :
    RGBcolor(r/2,g/2,b/2);
    bgnpolygon();
    vector[0]=x;vector[1]=y+s;v2s(vector);
    vector[0]=x+s;vector[1]=y;v2s(vector);
    vector[0]=x;vector[1]=y-s;v2s(vector);
    vector[0]=x-s;vector[1]=y;v2s(vector);
    endpolygon();
    RGBcolor(r,g,b);
    break;
  case 8 :
    bgnpolygon();
    vector[0]=x;vector[1]=y+s;v2s(vector);
    vector[0]=x+s;vector[1]=y;v2s(vector);
    vector[0]=x;vector[1]=y-s;v2s(vector);
    vector[0]=x-s;vector[1]=y;v2s(vector);
    endpolygon();
    break;
  case 9 :
    move2s(x,y+s);
    draw2s(x,y-s);
    move2s(x-s,y);
    draw2s(x+s,y-s);
    break;
  case 10 :
    move2s(x-s,y+s);
    draw2s(x+s,y-s);
    move2s(x-s,y-s);
    draw2s(x+s,y+s);
    break;
  }
}

int RasterizeFile (FILE *stream, AWindowRecord *myWindow)
{
  char *buffer;                                                 /* input buffer						*/
  long blockSize;                                               /* METABUFFERSIZE					*/
  long blockUsed;                                               /* actual buffer size used			*/
  long itemCounter;                                             /* number of commands in buffer		*/
  char *data;                                                           /* data pointer in buffer			*/
  short fx,fy;                                                  /* file screen size					*/
  int i,error,j,size;
  char opCode;
  short x,y,r,g,b,n,lw,ts,m,ms,si,w,x1,y1,x2,y2;
  short xx[SIZE],yy[SIZE];
  char s[CSIZE];
  unsigned char c;
  long l;
  short vector[2];

  /* get file parameters */
  rewind(stream);
  error = fread(&blockSize,4,1,stream); if (error!=1) return(1);       /* block size */
  error = fread(&fx,2,1,stream);            if (error!=1) return(1);       /* x size */
  error = fread(&fy,2,1,stream);        if (error!=1) return(1);       /* y size */

  /* set the window-coordinates */
  ortho2(0,fx-1,0,fy-1);

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
        move2s(x,y);
        break;

      case opDraw :
        memcpy(&x,data,2);
        data += 2;
        memcpy(&y,data,2);
        data += 2;
        draw2s(x,y);
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
        move2s(xx[0],yy[0]);
        for (j=1; j<n; j++) draw2s(xx[j],yy[j]);
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
        if (n<3) break;
        bgnpolygon();
        for (j=0; j<n; j++) { vector[0]=xx[j];vector[1]=yy[j];v2s(vector); }
        endpolygon();
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
        for (j=0; j<n; j++) Marker(m,ms,xx[j],yy[j]);
        break;

      case opText :
        memcpy(&n,data,2);
        data += 2;
        if (n>=CSIZE-1) {free(buffer); return(2);}
        memcpy(s,data,n);
        s[n] = 0;
        data += n;
        cmov2s(x,y);
        charstr(s);
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
        /*
                                                w = StringWidth((ConstStr255Param)s);
                                                MoveTo(x-w/2,y+ts/2);
                                                DrawString((ConstStr255Param)s);
         */
        cmov2s(x,y);
        charstr(s);
        break;

      case opSetLineWidth :
        memcpy(&n,data,2);
        data += 2;
        lw = n;
        linewidth(n);
        break;

      case opSetTextSize :
        memcpy(&n,data,2);
        data += 2;
        ts = n;
        /*
                                                TextSize(n);
         */
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
        c  = *((unsigned char *)data);
        cc = c;
        data++;
        RGBcolor(myWindow->red[c],myWindow->green[c],myWindow->blue[c]);
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
        myWindow->red[c]   = r;
        myWindow->green[c] = g;
        myWindow->blue[c]  = b;
        break;

      case opSetPalette :
        x = (short) (*((unsigned char *)data));
        data++;
        y = (short) (*((unsigned char *)data));
        data++;
        for (j=x; j<=y; j++)
        {
          r = (short) (*((unsigned char *)data));
          data++;
          g = (short) (*((unsigned char *)data));
          data++;
          b = (short) (*((unsigned char *)data));
          data++;
          myWindow->red[j]   = r;
          myWindow->green[j] = g;
          myWindow->blue[j]  = b;
        }
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
        linewidth(lw);
        RGBcolor(myWindow->red[c],myWindow->green[c],myWindow->blue[c]);
        move2s(x1,y1);
        draw2s(x2,y2);
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
        linewidth(lw);
        RGBcolor(myWindow->red[c],myWindow->green[c],myWindow->blue[c]);
        move2s(xx[0],yy[0]);
        for (j=1; j<n; j++) draw2s(xx[j],yy[j]);
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
        RGBcolor(myWindow->red[c],myWindow->green[c],myWindow->blue[c]);
        bgnpolygon();
        for (j=0; j<n; j++) { vector[0]=xx[j];vector[1]=yy[j];v2s(vector); }
        endpolygon();
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
        RGBcolor(myWindow->red[c],myWindow->green[c],myWindow->blue[c]);
        for (j=0; j<n; j++) Marker(m,ms,xx[j],yy[j]);
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
        RGBcolor(myWindow->red[c],myWindow->green[c],myWindow->blue[c]);
        cmov2s(x,y);
        charstr(s);
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
        RGBcolor(myWindow->red[c],myWindow->green[c],myWindow->blue[c]);
        cmov2s(x,y);
        charstr(s);
        break;

      case opShadedPolygon :
        memcpy(&n,data,2);
        data += 2;
        if (n>=SIZE) {free(buffer); return(2);}
        memcpy(&si,data,2);
        data += 2;
        r = ((int)myWindow->red[cc]  *(int)si)/1000;
        g = ((int)myWindow->green[cc]*(int)si)/1000;
        b = ((int)myWindow->blue[cc] *(int)si)/1000;
        if (n>=SIZE) {free(buffer); return(2);}
        size = n<<1;
        memcpy(xx,data,size);
        data += size;
        memcpy(yy,data,size);
        data += size;
        if (n<3) break;
        RGBcolor(r,g,b);
        bgnpolygon();
        for (j=0; j<n; j++) { vector[0]=xx[j];vector[1]=yy[j];v2s(vector); }
        endpolygon();
        RGBcolor(myWindow->red[cc],myWindow->green[cc],myWindow->blue[cc]);
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
  int error,i,j,sx,sy,from,to,dummy;
  int theEnd,sizeopt,pictureopt,filmopt,wait;
  FILE *stream;
  short fx,fy;
  long ox,oy,odx,ody;
  char name[80], msg[80];

  if (argc==1)
  {
    printf("usage: ugv <basename> [-f from to] [-s x y] [-p] [-w <n>]\n");
    return(0);
  }

  /* scan options */
  dummy = 0;
  sizeopt = 0;
  pictureopt = 0;
  filmopt = 0;
  i = 2;
  wait = 0;
  while (i<argc)
  {
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

    if (argv[i][1] == 'w')
    {
      if (i+1>=argc)
      {
        printf("not enough arguments for -w option\n");
        return(0);
      }
      sscanf(argv[i+1],"%d",&wait);
      i += 1;
    }

    if (argv[i][1] == 'f')
    {
      if (i+2>=argc)
      {
        printf("not enough arguments for -f option\n");
        return(0);
      }
      sscanf(argv[i+1],"%d",&from);
      sscanf(argv[i+2],"%d",&to);
      filmopt = 1;
      i += 2;
    }

    if (argv[i][1] == 'p')
    {
      pictureopt = 1;
    }

    i++;
  }

  if (!filmopt)
  {
    stream = fopen(argv[1],"rb");
    if (stream==NULL)
    {
      printf("Could not open file %s\n", argv[1]);
      return(0);
    }

    GetFileScreen(stream,&fx,&fy);
    if (!sizeopt) { sx=fx; sy=fy; }

    /* open a window */
    error = CreateApplicationWindow(&myWindow,argv[1],20,20,sx,sy);
    getorigin(&ox,&oy); getsize(&odx,&ody);

    /* draw something */
    RasterizeFile(stream,&myWindow);
    swapbuffers();
    fclose(stream);

    if (pictureopt)
    {
      sprintf(msg,"scrsave %s.rgb %d %d %d %d", argv[1], ox, ox+odx, oy, oy+ody);
      system(msg);
    }

    /* wait until escape pressed */
    while(!getvaluator(ESCKEY)) ;
  }
  else
  {
    sprintf(name,"%s.%04d",argv[1],from);
    stream = fopen(name,"rb");
    if (stream==NULL)
    {
      printf("Could not open file");
      return(0);
    }
    GetFileScreen(stream,&fx,&fy);
    fclose(stream);

    if (!sizeopt) { sx=fx; sy=fy; }

    /* open a window */
    error = CreateApplicationWindow(&myWindow,argv[1],20,20,sx,sy);
    if (error>0) return(0);
    getorigin(&ox,&oy);  getsize(&odx,&ody);

    while (1)
    {
      for (i=from; i<=to; i++)
      {
        /* construct file name */
        sprintf(name,"%s.%04d",argv[1],i);

        stream = fopen(name,"rb");
        if (stream==NULL)
        {
          printf("Could not open file %s\n", name);
          return(0);
        }

        /* draw something */
        cpack(0xffffff);
        clear();
        RasterizeFile(stream,&myWindow);
        swapbuffers();
        fclose(stream);

        for (j=0; j<wait*100; j++) dummy++;

        if (pictureopt)
        {
          sprintf(msg,"scrsave %s.rgb %d %d %d %d", name, ox, ox+odx, oy, oy+ody);
          system(msg);
        }
        if (getvaluator(ESCKEY)) goto exit;
      }
      if (pictureopt) goto exit;
    }

  }

exit:
  DisposeApplicationWindow(&myWindow);
  gexit();
  return(0);
}
