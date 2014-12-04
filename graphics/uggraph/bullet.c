// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  bullet.c                                                                                                      */
/*																			*/
/* Purpose:   The bullet plotter rasterizes lines and polygons in a local   */
/*            pixel buffer using -- in 3D -- the z buffer algorithm. Main   */
/*            advantage is that rasterization can be done in parallel and   */
/*            that merging the buffers is a simple reduction operation.     */
/*            Besides, the output device can flush the pixels all at once.  */
/*            This speeds up a remote X connection considerably.            */
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

#include <config.h>
#include "architecture.h"
#include "bullet.h"
#include "ugdevices.h"
#include "commands.h"
#include "heaps.h"
#include "general.h"
#include "pixel.h"

#ifdef ModelP
#include "ppif.h"
USING_PPIF_NAMESPACE
#endif

USING_UG_NAMESPACES

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

/* macros for accessing pixel buffer lines */
#define Z_BUFFER(y)  ((ZTYP  *)(ZBuffer)+(y)*Width)
#define P_BUFFER(y)  ((PIXEL *)(PBuffer)+(y)*Width)

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
/* definition of exported global variables                                  */
/*																			*/
/****************************************************************************/

INT NS_DIM_PREFIX BulletDim;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT Width;                    /* current pictures width              */
static INT Height;                   /* current pictures height             */
static INT NbPixels;                 /* number of pixels                    */
static INT Length;                   /* accumulated length of pixel & z buf.*/
static DOUBLE XShift;                /* shift for x coordinates             */
static DOUBLE YShift;                /* shift for y coordinates             */
static DOUBLE zOffsetFactor;         /* see comment in DrawTriangle         */
static INT MarkKey;                  /* for MarkTmpMem                      */
static OUTPUTDEVICE *OutputDevice;   /* current pictures output device      */
static void *PBuffer;                /* pixel buffer                        */
static void *ZBuffer;                /* z buffer                            */
#ifdef ModelP
static void *ABuffer;                /* auxiliary buffer                    */
#endif

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static void DrawLine(POINT p1, DOUBLE z1, POINT p2, DOUBLE z2, char c);


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

INT NS_DIM_PREFIX BulletOpen(PICTURE *picture, DOUBLE factor)
{
  HEAP *heap;
  ZTYP *z;
  PIXEL *p;
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
  if (BulletDim == 3)
    Length = NbPixels * (sizeof(PIXEL)+sizeof(ZTYP));
  else
    Length = NbPixels * sizeof(PIXEL);

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
  if (BulletDim == 3)
  {
    /* init z buffer */
    z = (ZTYP *)ZBuffer;
    for (i = 0; i < NbPixels; i++)
      *z++ = FAR_AWAY;

    PBuffer = z;
  }
  else
    PBuffer = ZBuffer;

  /* init pixel buffer */
  p = (PIXEL *)PBuffer;
  for (i = 0; i < NbPixels; i++) {
    p->cindex = OutputDevice->white;
    p->intensity = 255;
    p++;
  }

  return BULLET_OK;
}

/****************************************************************************/
/*D
   BulletClose - closes the bullet plotter

   SYNOPSIS:
   void BulletClose(void)

   PARAMETERS:
   .  none -

   DESCRIPTION:
   BulletClose closes the bullet plotter by releasing the buffers allocated
   with BulletOpen.

   RETURN VALUE:
   none
   D*/
/****************************************************************************/

void NS_DIM_PREFIX BulletClose(void)
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
  ZTYP  *z1, *z2;
  PIXEL *p1, *p2;

  if (BulletDim == 3)
  {
    z1 = (ZTYP  *)buffer1;        z2 = (ZTYP  *)buffer2;
    p1 = (PIXEL *)(z1+NbPixels);  p2 = (PIXEL *)(z2+NbPixels);

    for (i = 0; i < NbPixels; i++) {
      if (*z2 > *z1) {
        *p1 = *p2;
        *z1 = *z2;
      }
      p1++;  p2++;
      z1++;  z2++;
    }
  }
  else
  {
    p1 = (PIXEL *)buffer1;
    p2 = (PIXEL *)buffer2;
    for (i = 0; i < NbPixels; i++) {
      if (p2->cindex != OutputDevice->white)
        *p1 = *p2;
      p1++;  p2++;
    }
  }

}

/****************************************************************************/
/*D
   BulletPlot - plot pixel buffer

   SYNOPSIS:
   void BulletPlot(void)

   PARAMETERS:
   .  none -

   DESCRIPTION:
   BulletPlot plots the pixel buffer via the output device''s PlotPixelBuffer
   method. The parallel version first merges the buffers from all procs.

   RETURN VALUE:
   none
   D*/
/****************************************************************************/

static void FramePicture(void)
{
  POINT p1, p2;

  p1.x = 0; p1.y = 0; p2.x = Width-1; p2.y = 0;
  DrawLine(p1, -FAR_AWAY, p2, -FAR_AWAY, OutputDevice->black);
  p1.x = Width-1; p1.y = Height-1;
  DrawLine(p1, -FAR_AWAY, p2, -FAR_AWAY, OutputDevice->black);
  p2.x = 0; p2.y = Height-1;
  DrawLine(p1, -FAR_AWAY, p2, -FAR_AWAY, OutputDevice->black);
  p1.x = 0; p1.y = 0;
  DrawLine(p1, -FAR_AWAY, p2, -FAR_AWAY, OutputDevice->black);
}

void NS_DIM_PREFIX BulletPlot(void)
{
#ifdef ModelP
  INT son;
#endif
  void *scratch;

  /* reuse z buffer if possible */
  if (BulletDim == 3 && sizeof(ZTYP) >= 4)
    scratch = ZBuffer;
  else
    scratch = NULL;

  /* merge buffers */
#ifdef ModelP
  for (son = degree-1; son >= 0; son--) {
    GetConcentrate(son, ABuffer, Length);
    MergeBuffers(ZBuffer, ABuffer);
  }
  Concentrate(ZBuffer, Length);

  if (me == master)
#endif
  {
    FramePicture();
    (*OutputDevice->PlotPixelBuffer)(PBuffer, scratch, XShift, YShift,
                                     Width, Height);
  }
}

/*****************************************************************************

   DrawPoint - draw point in pixel buffer

*****************************************************************************/

static void DrawPoint(INT x, INT y, DOUBLE z, char c)
{
  ZTYP *zp;

  if (x < 0 || x >= Width || y < 0 || y >= Height) return;

  if (BulletDim == 3)
  {
    zp = Z_BUFFER(y)+x;
    if (z >= (*zp)-ZEPS*ABS(*zp)) {
      P_BUFFER(y)[x].cindex = c;
      P_BUFFER(y)[x].intensity = 255;
      *zp = z;
    }
  }
  else {
    P_BUFFER(y)[x].cindex = c;
    P_BUFFER(y)[x].intensity = 255;
  }
}

/*****************************************************************************

   DrawLine - draw line in pixel buffer

*****************************************************************************/

static void DrawLine(POINT p1, DOUBLE z1, POINT p2, DOUBLE z2, char c)
{
  POINT pt;
  DOUBLE m1, m2, X, Y, Z, dz, zt;
  INT x, y, dx, dy;

  if (p1.x == p2.x && p1.y == p2.y) {
    DrawPoint(p1.x, p1.y, MAX(z1, z2), c);
    return;
  }

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
  PIXEL *pp;
  ZTYP *pz;
  INT x;

  if (y < 0 || y >= Height) return;

  pp = P_BUFFER(y)+x1;
  pz = Z_BUFFER(y)+x1;

  if (x1 <= x2) {
    for (x = x1; x <= x2; x++) {
      if (x >= 0 && x < Width) {
        if (BulletDim == 3)
        {
          if (z >= *pz) {
            pp->cindex = c;
            pp->intensity = i*255.0;
            *pz = z;
          }
        }
        else
        {
          pp->cindex = c;
          pp->intensity = i*255.0;
        }
      }
      pp++;
      if (BulletDim == 3)
      {
        pz++;
        z += dz;
      }
    }
  }
  else {
    for (x = x1; x >=x2; x--) {
      if (x >= 0 && x < Width) {
        if (BulletDim == 3)
        {
          if (z >= *pz) {
            pp->cindex = c;
            pp->intensity = i*255.0;
            *pz = z;
          }
        }
        else
        {
          pp->cindex = c;
          pp->intensity = i*255.0;
        }
      }
      pp--;
      if (BulletDim == 3)
      {
        pz--;
        z -= dz;
      }
    }
  }
}

/*****************************************************************************

   DrawTriangle - draw filled triangle in pixel buffer

*****************************************************************************/

static void DrawTriangle(POINT p1, DOUBLE z1, POINT p2, DOUBLE z2,
                         POINT p3, DOUBLE z3, DOUBLE i, char c)
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
  /* should avoid line dropouts if surrounding is drawn later */

  o = zOffsetFactor * sqrt(dzNdx*dzNdx+dzNdy*dzNdy);

  /* rasterize lower partial triangle */
  if (y21 != 0) {
    mx2 = (DOUBLE)x21 / (DOUBLE)y21;
    x1 = x2 = (DOUBLE)p1.x + 0.5;
    z = z1-o;
    for (y = p1.y; y <= p2.y; y++) {
      DrawSpan((INT)x1, (INT)x2, y, z, dzNdx, i, c);
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

void NS_DIM_PREFIX BulletLine(DOUBLE *point1, DOUBLE *point2, long color)
{
  POINT p1, p2;
  DOUBLE z1, z2;

  p1.x = (INT)(point1[0] - XShift + 0.5);
  p1.y = (INT)(point1[1] - YShift + 0.5);
  p2.x = (INT)(point2[0] - XShift + 0.5);
  p2.y = (INT)(point2[1] - YShift + 0.5);
  if (BulletDim == 3) {
    z1 = point1[2];
    z2 = point2[2];
  }
  else
    z1 = z2 = 0.0;
  DrawLine(p1, z1, p2, z2, color);
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

void NS_DIM_PREFIX BulletPolyLine(DOUBLE *points, INT nb, long color)
{
  DOUBLE *p0, *p1;
  INT i;

  p0 = points;
  for (i = 0; i < nb-1; i++) {
    p1 = points+BulletDim;
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

void NS_DIM_PREFIX BulletPolygon(DOUBLE *points, INT nb, DOUBLE intensity, long color)
{
  DOUBLE z0, z1, z2;
  POINT p0, p1, p2;
  INT k;

  p0.x = (INT)((*points++) - XShift + 0.5);
  p0.y = (INT)((*points++) - YShift + 0.5);
  if (BulletDim == 3)
    z0   = *points++;
  else
    z0 = 0.0;
  for (k = 1; k < nb-1; k++) {
    p1.x = (INT)((*points++) - XShift + 0.5);
    p1.y = (INT)((*points++) - YShift + 0.5);
    if (BulletDim == 3)
      z1 = *points++;
    else
      z1 = 0.0;
    p2.x = (INT)((*points++) - XShift + 0.5);
    p2.y = (INT)((*points++) - YShift + 0.5);
    if (BulletDim == 3)
      z2 = *points;
    else
      z2 = 0.0;
    points -= 2;
    DrawTriangle(p0, z0, p1, z1, p2, z2, intensity, color);
  }
}
