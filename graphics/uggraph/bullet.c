// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  bullet.c                                                                                                      */
/*																			*/
/* Purpose:   The bullet plotter rasterizes lines and polygons in a pixel   */
/*            buffer using---in 3D---the z buffer algorithm. Suitable       */
/*            output devices can display the pixel buffer in a single       */
/*            operation which may be faster than the standard plotting      */
/*            routine.                                                      */
/*																			*/
/* Author:	  Michael Lampe                                                                                                 */
/*			  Institut fuer Computeranwendungen                                                     */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de                                                    */
/*																			*/
/* History:   24.2.98 begin, ug3-version									*/
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

#include "bullet.h"
#include "ugdevices.h"
#include "commands.h"
#include "heaps.h"
#ifdef ModelP
#include "parallel.h"
#include "ppif.h"
#endif

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define ZTYP         FLOAT            /* type for z buffer                  */
#define ZEPS         (5*FLT_EPSILON)  /* eps for ZTYP                       */
#define FAR_AWAY     -1E30            /* a large negative number from ZTYP  */

#define DM_ROWS      8                /* dither matrix rows, power of 2     */
#define DM_COLS      8                /* dither matrix cols, power of 2     */
#define DM_ENTRIES   (DM_ROWS*DM_COLS)

/* macros for accessing pixel buffer lines */
#define Z_BUFFER(y)  ((ZTYP *)(ZBuffer)+(y)*Width)
#define P_BUFFER(y)  ((char *)(PBuffer)+(y)*Width)

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {
  INT x;
  INT y;
} POINT;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT Width;                    /* current pictures width              */
static INT Height;                   /* current pictures height             */
static INT NbPixels;                 /* number of pixels                    */
static INT Length;                   /* accumulated length of pixel & z buf.*/
static DOUBLE XShift;                /* to shift x coordinates              */
static DOUBLE YShift;                /* to shift y coordinates              */
static DOUBLE zOffsetFactor;         /* see comment in DrawTriangle         */
static INT MarkKey;                  /* for MarkTmpMem                      */
static OUTPUTDEVICE *OutputDevice;   /* current pictures output device      */
static void *PBuffer;                /* pixel buffer                        */
static void *ZBuffer;                /* z buffer                            */
static void *ABuffer;                /* auxiliary buffer                    */

static INT DitherMatrix[DM_ROWS][DM_COLS] =
{{ 0, 32,  8, 40,  2, 34, 10, 42},
 {48, 16, 56, 24, 50, 18, 58, 26},
 {12, 44,  4, 36, 14, 46,  6, 38},
 {60, 28, 52, 20, 62, 30, 54, 22},
 { 3, 35, 11, 43,  1, 33,  9, 41},
 {51, 19, 59, 27, 49, 17, 57, 25},
 {15, 47,  7, 39, 13, 45,  5, 37},
 {63, 31, 55, 23, 61, 29, 53, 21}};

/****************************************************************************/
/*D
   BulletOpen - initialize the bullet plotter

   SYNOPSIS:
   INT BulletOpen(PICTURE *picture, DOUBLE factor)

   PARAMETERS:
   .  picture  - the picture to plot
   .  factor   - factor for z offset

   DESCRIPTION:
   BulletOpen initializes the bullet plotter. It checks whether the pictures
   output device supports the bullet plotter and allocates the buffers for
   plotting.

   RETURN VALUE:
   INT
   .n     BULLET_OK:    ok
   .n     BULLET_CANT:  no suitable output device
   .n     BULLET_NOMEM: not enough memory
   D*/
/****************************************************************************/

INT BulletOpen(PICTURE *picture, DOUBLE factor)
{
  HEAP *heap;
  ZTYP *z;
  char *p;
  INT i, err;

  /* remember picture data */
  OutputDevice = PIC_OUTPUTDEV(picture);
  XShift = (DOUBLE)PIC_GLL(picture)[0];
  YShift = (DOUBLE)PIC_GUR(picture)[1];
  Width  = PIC_GUR(picture)[0] - PIC_GLL(picture)[0] + 1;
  Height = PIC_GLL(picture)[1] - PIC_GUR(picture)[1] + 1;

  zOffsetFactor = factor;

  /* check output device */
  if (OutputDevice->PlotPixelBuffer == NULL)
    return BULLET_CANT;

  /* allocate buffers */
  NbPixels = Width*Height;
#ifdef __THREEDIM__
  Length = NbPixels + NbPixels*sizeof(ZTYP);
#else
  Length = NbPixels;
#endif
  heap = GetCurrentMultigrid()->theHeap;
  MarkTmpMem(heap, &MarkKey);
  err =  ((ZBuffer = GetTmpMem(heap, Length, MarkKey)) == NULL);
#ifdef ModelP
  err = UG_GlobalMaxINT(err);
#endif
  if (err) {
    ReleaseTmpMem(heap, MarkKey);
    return BULLET_NOMEM;
  }
#ifdef ModelP
  err = ((ABuffer = GetTmpMem(heap, Length, MarkKey)) == NULL);
  if (UG_GlobalMaxINT(err)) {
    ReleaseTmpMem(heap, MarkKey);
    return BULLET_NOMEM;
  }
#endif
#ifdef __THREEDIM__

  /* init z buffer */
  z = (ZTYP *)ZBuffer;
  for (i = 0; i < NbPixels; i++)
    *z++ = FAR_AWAY;

  /* init pixel buffer */
  p = PBuffer = (char *)z;
#else
  p = PBuffer = (char *)ZBuffer;
#endif
  for (i = 0; i < NbPixels; i++)
    *p++ = (char)OutputDevice->white;

  return BULLET_OK;
}

/****************************************************************************/
/*D
   BulletClose - closes the bullet plotter

   SYNOPSIS:
   void BulletClose(void)

   PARAMETERS:
   .  none

   DESCRIPTION:
   BulletClose closes the bullet plotter by releasing the buffers allocated
   with BulletOpen.

   RETURN VALUE:
   none
   D*/
/****************************************************************************/

void BulletClose(void)
{
  HEAP *heap;

  heap = GetCurrentMultigrid()->theHeap;
  ReleaseTmpMem(heap, MarkKey);
}

/*****************************************************************************

   MergeBuffers - merge two buffers by plotting the 2nd on the 1st one.

*****************************************************************************/

static void MergeBuffers(void *buffer1, void *buffer2)
{
  INT i;
  char *p1, *p2;
#ifdef __THREEDIM__
  ZTYP *z1, *z2;

  z1 = (ZTYP *)buffer1;        z2 = (ZTYP *)buffer2;
  p1 = (char *)(z1+NbPixels);  p2 = (char *)(z2+NbPixels);

  for (i = 0; i < NbPixels; i++) {
    if (*z2 > *z1) {
      *p1 = *p2;
      *z1 = *z2;
    }
    p1++;  p2++;
    z1++;  z2++;
  }
#else
  p1 = (char *)buffer1;
  p2 = (char *)buffer2;
  for (i = 0; i < NbPixels; i++) {
    if (*p2 != (char)OutputDevice->white)
      *p1 = *p2;
    p1++;  p2++;
  }
#endif
}

/****************************************************************************/
/*D
   BulletPlot - plot pixel buffer

   SYNOPSIS:
   void BulletPlot(void)

   PARAMETERS:
   .  none

   DESCRIPTION:
   BulletPlot plots the pixel buffer via the output device's PlotPixelBuffer
   method. The parallel version first merges the buffers from all procs.

   RETURN VALUE:
   none
   D*/
/****************************************************************************/

void BulletPlot(void)
{
  INT son;
  void *data;

  /* reuse z buffer if possible */
#ifdef __THREEDIM__
  if (sizeof(ZTYP) >= 4)
    data = ZBuffer;
  else
#endif
  data = NULL;

  /* merge buffers */
#ifdef ModelP
  for (son = degree-1; son >= 0; son--) {
    GetConcentrate(son, ABuffer, Length);
    MergeBuffers(ZBuffer, ABuffer);
  }
  Concentrate(ZBuffer, Length);

  /* plot buffer */
  if (me == master)
#endif
  (*OutputDevice->PlotPixelBuffer)(PBuffer, data, NbPixels,
                                   (int)XShift, (int)YShift, (int)Width, (int)Height);
}

/*****************************************************************************

   DrawPoint - draw point in pixel buffer

*****************************************************************************/

static void DrawPoint(INT x, INT y, DOUBLE z, char c)
{
  ZTYP *zp;

  if (x < 0 || x >= Width || y < 0 || y >= Height) return;
#ifdef __THREEDIM__
  zp = Z_BUFFER(y)+x;
  if (z >= (*zp)-ZEPS*ABS(*zp)) {
    P_BUFFER(y)[x] = c;
    *zp = z;
  }
#else
  P_BUFFER(y)[x]=c;
#endif
}

/*****************************************************************************

   DrawLine - draw line in pixel buffer

*****************************************************************************/

static void DrawLine(POINT p1, DOUBLE z1, POINT p2, DOUBLE z2, char c)
{
  POINT pt;
  DOUBLE m1, m2, X, Y, Z, dz, zt;
  INT x, y, dx, dy;

  if (p1.x == p2.x && p1.y == p2.y) return;

  dx = p2.x - p1.x;
  dy = p2.y - p1.y;
  dz = z2 - z1;
  if (ABS(dx) >= ABS(dy)) {
    if (p1.x > p2.x) {
      SWAP(p1, p2, pt);
      SWAP(z1, z2, zt);
    }
    m1 = (DOUBLE)dy/(DOUBLE)dx;
    m2 = dz/(DOUBLE)dx;
    Y  = (DOUBLE)p1.y + 0.5;
    Z  = z1;
    for (x = p1.x; x <= p2.x; x++) {
      DrawPoint(x, (INT)Y, Z, c);
      Y += m1;
      Z += m2;
    }
  }
  else {
    if (p1.y > p2.y) {
      SWAP(p1, p2, pt);
      SWAP(z1, z2, zt);
    }
    m1 = (DOUBLE)dx/(DOUBLE)dy;
    m2 = dz/(DOUBLE)dy;
    X  = (DOUBLE)p1.x + 0.5;
    Z  = z1;
    for (y = p1.y; y <= p2.y; y++) {
      DrawPoint((INT)X, y, Z, c);
      X += m1;
      Z += m2;
    }
  }
}

/*****************************************************************************

   DrawSpan - draw horizontal line in pixel buffer

*****************************************************************************/

static void DrawSpan(INT x1, INT x2, INT y, DOUBLE z, DOUBLE dz, DOUBLE i, char c)
{
  char *pp;
  ZTYP *pz;
  INT *pd, x;

  if (y < 0 || y >= Height) return;

  pp = P_BUFFER(y)+x1;
  pz = Z_BUFFER(y)+x1;
  pd = DitherMatrix[y & (DM_ROWS-1)];

  if (x1 <= x2) {
    for (x = x1; x <= x2; x++) {
      if (x >= 0 && x < Width) {
#ifdef __THREEDIM__
        if (z >= *pz) {
          if (pd[x & (DM_COLS-1)] >= (INT)(DM_ENTRIES*i+0.5))
            *pp = (char)OutputDevice->black;
          else
            *pp = c;
          *pz = z;
        }
#else
        *pp = c;
#endif
      }
      pp++;
#ifdef __THREEDIM__
      pz++;
      z += dz;
#endif
    }
  }
  else {
    for (x = x1; x >=x2; x--) {
      if (x >= 0 && x < Width) {
#ifdef __THREEDIM__
        if (z >= *pz) {
          if (pd[x & (DM_COLS-1)] >= (INT)(DM_ENTRIES*i+0.5))
            *pp = (char)OutputDevice->black;
          else
            *pp = c;
          *pz = z;
        }
#else
        *pp = c;
#endif
      }
      pp--;
#ifdef __THREEDIM__
      pz--;
      z -= dz;
#endif
    }
  }
}

/*****************************************************************************

   DrawTriangle - draw filled triangle in pixel buffer

*****************************************************************************/

static void DrawTriangle(POINT p1, DOUBLE z1, POINT p2, DOUBLE z2, POINT p3, DOUBLE z3, DOUBLE i, char c)
{
  POINT pt;
  DOUBLE z31, z21, mx1, mx2, mz, dzNdx, dzNdy, x1, x2, z, zt, o;
  INT y31, y21, y32, x31, x32, x21, D, y;

  /* make p1.y <= p2.y <= p3.y */
  if (p1.y > p2.y) {
    SWAP(p1, p2, pt);
    SWAP(z1, z2, zt);
  }
  if (p1.y > p3.y) {
    SWAP(p1, p3, pt);
    SWAP(z1, z3, zt);
  }
  if (p2.y > p3.y) {
    SWAP(p2, p3, pt);
    SWAP(z2, z3, zt);
  }

  y31 = p3.y - p1.y;         /* height of long edge        */
  y21 = p2.y - p1.y;         /* height of lower short edge */
  y32 = p3.y - p2.y;         /* height of upper short edge */
  x31 = p3.x - p1.x;
  x32 = p3.x - p2.x;
  x21 = p2.x - p1.x;

  D   = y21*x31-y31*x21;

  /* degenerated triangle ? */
  if (D == 0) return;

  z31 = z3 - z1;
  z21 = z2 - z1;
  mx1 = (DOUBLE)x31 / (DOUBLE)y31;
  mz  = z31 / (DOUBLE)y31;

  /* grad z */
  dzNdx  = (DOUBLE)(y21*z31-z21*y31)/(DOUBLE)D;
  dzNdy  = (DOUBLE)(x31*z21-z31*x21)/(DOUBLE)D;

  /* offset for z coords proportional to max. depth slope     */
  /* should avoid line dropouts if sorrounding is drawn later */

  o = zOffsetFactor * sqrt(dzNdx*dzNdx+dzNdy*dzNdy);

  /* rasterize lower partial triangle */
  if (y21 != 0) {
    mx2 = (DOUBLE)x21 / (DOUBLE)y21;
    x1 = x2 = (DOUBLE)p1.x + 0.5;
    z = z1-o;
    for (y = p1.y; y <= p2.y; y++) {
      DrawSpan((INT)x1, (INT)x2, y, z, dzNdx, i,  c);
      x1 += mx1;
      x2 += mx2;
      z  += mz;
    }
  }

  /* rasterize upper partial triangle */
  if (y32 != 0) {
    mx2 = (DOUBLE)x32 / (DOUBLE)y32;
    x1 = x2 = (DOUBLE)p3.x + 0.5;
    z = z3-o;
    for (y = p3.y; y >= p2.y; y--) {
      DrawSpan((INT)x1, (INT)x2, y, z, dzNdx, i, c);
      x1 -= mx1;
      x2 -= mx2;
      z  -= mz;
    }
  }
}

/****************************************************************************/
/*D
   BulletLine - draw a line in pixel buffer

   SYNOPSIS:
   void BulletLine(DOUBLE *point1, DOUBLE *point2, long color)

   PARAMETERS:
   .  point1 - pointer to coordinates of 1st point
   .  point2 - pointer to coordinates of 2nd point
   .  color  - color

   DESCRIPTION:
   none

   RETURN VALUE:
   none
   D*/
/****************************************************************************/

void BulletLine(DOUBLE *point1, DOUBLE *point2, long color)
{
  POINT pp1, pp2;
  DOUBLE zz1, zz2;

  pp1.x = (INT)(point1[0] - XShift + 0.5);
  pp1.y = (INT)(point1[1] - YShift + 0.5);
  pp2.x = (INT)(point2[0] - XShift + 0.5);
  pp2.y = (INT)(point2[1] - YShift + 0.5);
#ifdef __THREEDIM__
  zz1   = point1[2];
  zz2   = point2[2];
#else
  zz1 = zz2 = 0.0;
#endif
  DrawLine(pp1, zz1, pp2, zz2, (char)color);
}

/****************************************************************************/
/*D
   BulletPolyLine - draw a polygon in pixel buffer

   SYNOPSIS:
   void BulletPolyLine(DOUBLE *points, INT nb, long color)

   PARAMETERS:
   .  points - pointer to coordinates of points
   .  nb     - number of points
   .  color  -

   DESCRIPTION:
   none

   RETURN VALUE:
   none
   D*/
/****************************************************************************/

void BulletPolyLine(DOUBLE *points, INT nb, long color)
{
  DOUBLE *p0, *p1;
  INT i, j;

  p0 = points;
  for (i = 0; i < nb-1; i++) {
    p1 = points+DIM;
    BulletLine(points, p1, color);
    points = p1;
  }
  BulletLine(p0, points, color);
}

/****************************************************************************/
/*D
   BulletPolygon - draw a filled (convex) polygon in pixel buffer

   SYNOPSIS:
   void BulletPolygon(DOUBLE *points, INT nb, DOUBLE intensity, long color)

   PARAMETERS:
   .  points    - pointer to coordinates of points
   .  nb        - number of points
   .  intensity -
   .  color     -

   DESCRIPTION:
   none

   RETURN VALUE:
   none
   D*/
/****************************************************************************/

void BulletPolygon(DOUBLE *points, INT nb, DOUBLE intensity, long color)
{
  DOUBLE z0, z1, z2;
  POINT p0, p1, p2;
  INT k;

  p0.x = (INT)((*points++) - XShift + 0.5);
  p0.y = (INT)((*points++) - YShift + 0.5);
#ifdef __THREEDIM__
  z0   = *points++;
#else
  z0 = 0.0;
#endif
  for (k = 1; k < nb-1; k++) {
    p1.x = (INT)((*points++) - XShift + 0.5);
    p1.y = (INT)((*points++) - YShift + 0.5);
#ifdef __THREEDIM__
    z1   = *points++;
#else
    z1   = 0.0;
#endif
    p2.x = (INT)((*points++) - XShift + 0.5);
    p2.y = (INT)((*points++) - YShift + 0.5);
#ifdef __THREEDIM__
    z2   = *points;
#else
    z2   = 0.0;
#endif
    points -= 2;
    DrawTriangle(p0, z0, p1, z1, p2, z2, intensity, (char)color);
  }
}
