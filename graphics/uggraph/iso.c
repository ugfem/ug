// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  iso.c                                                                                                         */
/*																			*/
/* Purpose:   Extract isosurface from zoo element                           */
/*																			*/
/* Author:	  Michael Lampe                                                                                                 */
/*			  IWR - Technische Simulation                                   */
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69120 Heidelberg                                              */
/*			  Email: Michael.Lampe@iwr.uni-heidelberg.de                    */
/*																			*/
/* History:   30.10.04 begin, ug3-version                                   */
/*																			*/
/* Remarks:   Simple version that just decomposes hexahedra into 6 pyramids */
/*            and then into 12 tetrahedra.                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <limits.h>
#include <assert.h>
#include "general.h"
#include "iso.h"

USING_UG_NAMESPACES

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {
  double x[4][3];
  double v[4];
} TET;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static int Pyr[2][8] = {
  {0, 1, 2, 4, 0, 2, 3, 4},
  {0, 1, 3, 4, 1, 2, 3, 4}
};

static int Pri[8][12] = {
  {0, 4, 5, 3, 1, 4, 2, 0, 4, 5, 2, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 4, 5, 3, 1, 4, 5, 0, 1, 5, 2, 0},
  {1, 5, 3, 4, 0, 5, 3, 1, 0, 2, 5, 1},
  {0, 1, 2, 4, 2, 5, 3, 4, 0, 2, 3, 4},
  {2, 3, 4, 5, 0, 3, 1, 2, 1, 3, 4, 2},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {1, 5, 3, 4, 0, 2, 3, 1, 2, 5, 3, 1}
};

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


static void Interpolate(double *x, TET *tet, double lambda, int i, int j)
{
  double t = (lambda-tet->v[i])/(tet->v[j]-tet->v[i]);

  x[0] = tet->x[i][0]+t*(tet->x[j][0]-tet->x[i][0]);
  x[1] = tet->x[i][1]+t*(tet->x[j][1]-tet->x[i][1]);
  x[2] = tet->x[i][2]+t*(tet->x[j][2]-tet->x[i][2]);
}

static void ExtractTet(TET *tet, double lambda, POLY *poly)
{
  int index = 0;

  if (tet->v[0] >= lambda) index |= 1;
  if (tet->v[1] >= lambda) index |= 2;
  if (tet->v[2] >= lambda) index |= 4;
  if (tet->v[3] >= lambda) index |= 8;

  switch (index)
  {
  case 0x01 :
  case 0x0E :
    Interpolate(poly->x[0], tet, lambda, 0, 1);
    Interpolate(poly->x[1], tet, lambda, 0, 2);
    Interpolate(poly->x[2], tet, lambda, 0, 3);
    poly->n = 3;
    break;
  case 0x02 :
  case 0x0D :
    Interpolate(poly->x[0], tet, lambda, 1, 0);
    Interpolate(poly->x[1], tet, lambda, 1, 2);
    Interpolate(poly->x[2], tet, lambda, 1, 3);
    poly->n = 3;
    break;
  case 0x03 :
  case 0x0C :
    Interpolate(poly->x[0], tet, lambda, 0, 2);
    Interpolate(poly->x[1], tet, lambda, 1, 2);
    Interpolate(poly->x[2], tet, lambda, 1, 3);
    Interpolate(poly->x[3], tet, lambda, 0, 3);
    poly->n = 4;
    break;
  case 0x04 :
  case 0x0B :
    Interpolate(poly->x[0], tet, lambda, 0, 2);
    Interpolate(poly->x[1], tet, lambda, 1, 2);
    Interpolate(poly->x[2], tet, lambda, 2, 3);
    poly->n = 3;
    break;
  case 0x05 :
  case 0x0A :
    Interpolate(poly->x[0], tet, lambda, 0, 1);
    Interpolate(poly->x[1], tet, lambda, 1, 2);
    Interpolate(poly->x[2], tet, lambda, 2, 3);
    Interpolate(poly->x[3], tet, lambda, 0, 3);
    poly->n = 4;
    break;
  case 0x06 :
  case 0x09 :
    Interpolate(poly->x[0], tet, lambda, 0, 1);
    Interpolate(poly->x[1], tet, lambda, 1, 3);
    Interpolate(poly->x[2], tet, lambda, 2, 3);
    Interpolate(poly->x[3], tet, lambda, 0, 2);
    poly->n = 4;
    break;
  case 0x07 :
  case 0x08 :
    Interpolate(poly->x[0], tet, lambda, 0, 3);
    Interpolate(poly->x[1], tet, lambda, 1, 3);
    Interpolate(poly->x[2], tet, lambda, 2, 3);
    poly->n = 3;
    break;
  default :
    poly->n = 0;
  }
}

static int SplitSide(CELL *cell, int i0, int i1, int i2, int i3)
{
  int a[4], i, k, m;

  k = 0;
  m = cell->order[i0];
  a[0] = i0; a[1] = i1; a[2] = i2, a[3] = i3;
  for (i = 1; i < 4; i++)
    if (cell->order[a[i]] <= m) {
      m = cell->order[a[i]];
      k = i;
    }
  return a[k];
}

static int DecomposePyr(CELL *cell)
{
  return SplitSide(cell, 0, 1, 2, 3) & 1;
}

static int DecomposePri(CELL *cell)
{
  return
    ((SplitSide(cell, 0, 1, 4, 3) & 1) << 0) |
    ((SplitSide(cell, 1, 2, 5, 4) & 1) << 1) |
    ((SplitSide(cell, 2, 0, 3, 5) & 2) << 1);
}

static void CenterNode(CELL *cell, double *xc, double *vc)
{
  int i;

  xc[0] = xc[1] = xc[2] = *vc = 0.0;
  for (i = 0; i < cell->n; i++) {
    xc[0] += cell->x[i][0];
    xc[1] += cell->x[i][1];
    xc[2] += cell->x[i][2];
    *vc   += cell->v[i];
  }
  xc[0] /= cell->n;
  xc[1] /= cell->n;
  xc[2] /= cell->n;
  *vc   /= cell->n;
}

static void CopyNodes(TET *tet, CELL *cell, int *t)
{
  int i;

  for (i = 0; i < 4; i++) {
    memcpy(tet->x[i], cell->x[t[i]], 3*sizeof(double));
    memcpy(tet->v+i, &cell->v[t[i]],   sizeof(double));
  }
}

static void CopyNodes2(CELL *pyr, CELL *cell, int i0, int i1, int i2, int i3,
                       double *xc, double vc)
{
  memcpy(pyr->x[0], cell->x[i0], 3*sizeof(double));
  memcpy(pyr->x[1], cell->x[i1], 3*sizeof(double));
  memcpy(pyr->x[2], cell->x[i2], 3*sizeof(double));
  memcpy(pyr->x[3], cell->x[i3], 3*sizeof(double));
  memcpy(pyr->x[4], xc,          3*sizeof(double));

  pyr->v[0] = cell->v[i0];
  pyr->v[1] = cell->v[i1];
  pyr->v[2] = cell->v[i2];
  pyr->v[3] = cell->v[i3];
  pyr->v[4] = vc;

  pyr->order[0] = cell->order[i0];
  pyr->order[1] = cell->order[i1];
  pyr->order[2] = cell->order[i2];
  pyr->order[3] = cell->order[i3];

  pyr->n = 5;
}

void NS_DIM_PREFIX ExtractElement(CELL *cell, double lambda, POLY *poly, int *npoly)
{
  int k;
  double xc[3], vc;
  TET tet;
  CELL pyr;

  switch (cell->n)
  {
  case 4 :
    memcpy(tet.x, cell->x, 12*sizeof(double));
    memcpy(tet.v, cell->v,  4*sizeof(double));
    ExtractTet(&tet, lambda, poly);
    *npoly = 1;
    break;
  case 5 :
    k = DecomposePyr(cell);
    CopyNodes(&tet, cell, Pyr[k]+ 0);
    ExtractTet(&tet, lambda, poly+0);
    CopyNodes(&tet, cell, Pyr[k]+ 4);
    ExtractTet(&tet, lambda, poly+1);
    *npoly = 2;
    break;
  case 6 :
    k = DecomposePri(cell);
    CopyNodes(&tet, cell, Pri[k]+ 0);
    ExtractTet(&tet, lambda, poly+0);
    CopyNodes(&tet, cell, Pri[k]+ 4);
    ExtractTet(&tet, lambda, poly+1);
    CopyNodes(&tet, cell, Pri[k]+ 8);
    ExtractTet(&tet, lambda, poly+2);
    *npoly = 3;
    break;
  case 8 :
    CenterNode(cell, xc, &vc);
    CopyNodes2(&pyr, cell, 0, 4, 5, 1, xc, vc);
    ExtractElement(&pyr, lambda, poly+ 0, npoly);
    CopyNodes2(&pyr, cell, 1, 5, 6, 2, xc, vc);
    ExtractElement(&pyr, lambda, poly+ 2, npoly);
    CopyNodes2(&pyr, cell, 2, 6, 7, 3, xc, vc);
    ExtractElement(&pyr, lambda, poly+ 4, npoly);
    CopyNodes2(&pyr, cell, 0, 3, 7, 4, xc, vc);
    ExtractElement(&pyr, lambda, poly+ 6, npoly);
    CopyNodes2(&pyr, cell, 0, 1, 2, 3, xc, vc);
    ExtractElement(&pyr, lambda, poly+ 8, npoly);
    CopyNodes2(&pyr, cell, 4, 7, 6, 5, xc, vc);
    ExtractElement(&pyr, lambda, poly+10, npoly);
    *npoly = 12;
    break;
  default :
    assert(0);
  }
}
