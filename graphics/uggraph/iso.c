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
/* History:   19.10.04 begin, ug3-version                                   */
/*																			*/
/* Remarks:   - extracting an isosurface from a tetrahedron is easy, e.g.:  */
/*              "Tetra-Cubes: An algorithm to generate 3D isosurfaces based */
/*              upon tetrahedra" --  B. Carneiro, C. Silva, and A. Kaufman, */
/*              Proc. SIBGRAPI '96, pp. 205-210                             */
/*                                                                          */
/*            - consistent decomposition of zoo meshes into tetrahedra is   */
/*              explained here: "Approximate Volume Rendering for Curvi-    */
/*              linear and Unstructured Grids by Hardware-Assisted Poly-    */
/*              hedron Projection" --  N. Max, P. Williams, C. Silva,       */
/*              International Journal of Imaging Systems (2000), 11:53-61   */
/*                                                                          */
/*            - corners of elements are numbered as in UG                   */
/*																			*/
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

static int Hex[64][24] = {  /* computed with AH's procedural implementation */
  {7, 0, 3, 6, 0, 3, 2, 6, 0, 2, 1, 6, 7, 0, 4, 6, 0, 4, 5, 6, 0, 5, 1, 6},
  {7, 0, 3, 6, 3, 2, 1, 6, 0, 3, 1, 6, 7, 0, 4, 6, 0, 4, 5, 6, 0, 5, 1, 6},
  {7, 0, 3, 6, 0, 3, 2, 6, 0, 2, 1, 6, 7, 0, 4, 6, 4, 5, 1, 6, 0, 4, 1, 6},
  {0, 3, 7, 1, 3, 2, 6, 1, 3, 6, 7, 1, 0, 4, 7, 1, 4, 5, 6, 1, 4, 6, 7, 1},
  {3, 2, 6, 0, 2, 6, 5, 0, 1, 2, 5, 0, 3, 7, 6, 0, 4, 7, 6, 0, 4, 6, 5, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {4, 5, 7, 0, 1, 5, 7, 0, 1, 7, 3, 0, 3, 1, 2, 7, 1, 2, 6, 7, 1, 6, 5, 7},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {2, 7, 6, 1, 4, 7, 6, 1, 4, 6, 5, 1, 2, 3, 7, 1, 0, 3, 7, 1, 0, 7, 4, 1},
  {6, 4, 5, 2, 0, 1, 5, 2, 0, 5, 4, 2, 6, 4, 7, 2, 0, 3, 7, 2, 0, 7, 4, 2},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {1, 4, 5, 2, 4, 7, 6, 2, 4, 6, 5, 2, 1, 0, 4, 2, 0, 3, 7, 2, 0, 7, 4, 2},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {5, 0, 1, 6, 0, 1, 2, 6, 0, 2, 3, 6, 5, 0, 4, 6, 4, 7, 3, 6, 0, 4, 3, 6},
  {0, 1, 5, 3, 1, 2, 6, 3, 1, 6, 5, 3, 0, 4, 5, 3, 4, 7, 6, 3, 4, 6, 5, 3},
  {2, 0, 1, 6, 1, 5, 4, 6, 0, 1, 4, 6, 2, 0, 3, 6, 3, 7, 4, 6, 0, 3, 4, 6},
  {0, 1, 3, 4, 1, 2, 3, 6, 1, 4, 5, 6, 3, 4, 6, 7, 1, 3, 4, 6,-1,-1,-1,-1},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 1, 5, 3, 2, 6, 5, 3, 1, 2, 5, 3, 0, 4, 5, 3, 4, 7, 6, 3, 4, 6, 5, 3},
  {0, 1, 2, 4, 5, 6, 2, 4, 1, 5, 2, 4, 0, 2, 3, 4, 6, 7, 3, 4, 2, 6, 3, 4},
  {2, 1, 5, 3, 4, 5, 1, 3, 0, 4, 1, 3, 2, 6, 5, 3, 6, 5, 4, 3, 7, 6, 4, 3},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {7, 2, 6, 4, 1, 5, 6, 4, 1, 6, 2, 4, 7, 2, 3, 4, 0, 3, 2, 4, 1, 0, 2, 4},
  {2, 7, 6, 1, 4, 7, 6, 1, 4, 6, 5, 1, 2, 3, 7, 1, 3, 7, 4, 1, 0, 3, 4, 1},
  {3, 0, 4, 2, 0, 4, 5, 2, 0, 5, 1, 2, 3, 7, 4, 2, 6, 5, 4, 2, 7, 6, 4, 2},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {3, 0, 4, 2, 4, 5, 1, 2, 0, 4, 1, 2, 3, 7, 4, 2, 6, 5, 4, 2, 7, 6, 4, 2},
  {7, 2, 6, 4, 5, 6, 2, 4, 1, 5, 2, 4, 7, 2, 3, 4, 1, 0, 3, 4, 1, 3, 2, 4},
  {4, 5, 7, 0, 5, 7, 3, 0, 1, 5, 3, 0, 7, 5, 6, 3, 1, 2, 6, 3, 1, 6, 5, 3},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 3, 7, 1, 3, 2, 6, 1, 3, 6, 7, 1, 0, 4, 7, 1, 5, 6, 7, 1, 4, 5, 7, 1},
  {6, 3, 2, 5, 0, 1, 2, 5, 0, 2, 3, 5, 6, 3, 7, 5, 0, 4, 7, 5, 0, 7, 3, 5},
  {6, 3, 2, 5, 1, 2, 3, 5, 0, 1, 3, 5, 6, 3, 7, 5, 0, 4, 7, 5, 0, 7, 3, 5},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 2, 3, 7, 0, 1, 2, 5, 0, 4, 5, 7, 2, 5, 6, 7, 0, 2, 5, 7,-1,-1,-1,-1},       /* BUG!? 101000 */
  {3, 0, 1, 7, 0, 1, 5, 7, 0, 5, 4, 7, 3, 1, 2, 7, 1, 2, 6, 7, 1, 6, 5, 7},
  {4, 1, 5, 7, 1, 5, 6, 7, 1, 6, 2, 7, 4, 1, 0, 7, 0, 3, 2, 7, 1, 0, 2, 7},
  {2, 7, 6, 1, 7, 6, 5, 1, 4, 7, 5, 1, 2, 3, 7, 1, 0, 3, 7, 1, 0, 7, 4, 1},
  {0, 2, 3, 7, 0, 1, 2, 5, 0, 4, 5, 7, 2, 5, 6, 7, 0, 2, 5, 7,-1,-1,-1,-1},       /* BUG!? 101100 */
  {0, 2, 3, 7, 0, 1, 2, 5, 0, 4, 5, 7, 2, 5, 6, 7, 0, 2, 5, 7,-1,-1,-1,-1},       /* BUG!? 101101 */
  {1, 4, 5, 2, 7, 6, 5, 2, 4, 7, 5, 2, 1, 0, 4, 2, 0, 3, 7, 2, 0, 7, 4, 2},
  {4, 1, 5, 7, 5, 6, 2, 7, 1, 5, 2, 7, 4, 1, 0, 7, 1, 0, 3, 7, 1, 3, 2, 7},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 1, 5, 3, 1, 2, 6, 3, 1, 6, 5, 3, 0, 4, 5, 3, 7, 6, 5, 3, 4, 7, 5, 3},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {5, 4, 7, 1, 3, 7, 4, 1, 0, 3, 4, 1, 5, 6, 7, 1, 6, 7, 3, 1, 2, 6, 3, 1},
  {6, 3, 2, 5, 0, 1, 2, 5, 0, 2, 3, 5, 6, 3, 7, 5, 4, 7, 3, 5, 0, 4, 3, 5},
  {0, 1, 5, 3, 2, 6, 5, 3, 1, 2, 5, 3, 0, 4, 5, 3, 7, 6, 5, 3, 4, 7, 5, 3},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {2, 1, 5, 3, 4, 5, 1, 3, 0, 4, 1, 3, 2, 6, 5, 3, 7, 6, 5, 3, 7, 5, 4, 3},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {2, 7, 6, 1, 7, 6, 5, 1, 4, 7, 5, 1, 2, 3, 7, 1, 3, 7, 4, 1, 0, 3, 4, 1},
  {3, 0, 4, 2, 0, 4, 5, 2, 0, 5, 1, 2, 3, 7, 4, 2, 7, 6, 5, 2, 7, 5, 4, 2},
  {4, 0, 3, 5, 3, 2, 1, 5, 0, 3, 1, 5, 4, 3, 7, 5, 7, 6, 2, 5, 3, 7, 2, 5},
  {3, 0, 4, 2, 4, 5, 1, 2, 0, 4, 1, 2, 3, 7, 4, 2, 7, 6, 5, 2, 7, 5, 4, 2},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
};

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* Here goes. Someone please put in some documentation---thanks!            */
/*                                                                          */
/****************************************************************************/

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
  int a[4], i, k, m = INT_MAX;

  a[0] = i0; a[1] = i1; a[2] = i2, a[3] = i3;
  for (i = 0; i < 4; i++)
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

static int DecomposeHex(CELL *cell)
{
  int k, index = 0;

  k = SplitSide(cell, 0, 1, 2, 3);
  if (k == 1 || k == 3) index |= (1<<0);
  k = SplitSide(cell, 0, 1, 5, 4);
  if (k == 1 || k == 4) index |= (1<<1);
  k = SplitSide(cell, 1, 2, 6, 5);
  if (k == 2 || k == 5) index |= (1<<2);
  k = SplitSide(cell, 2, 3, 7, 6);
  if (k == 2 || k == 7) index |= (1<<3);
  k = SplitSide(cell, 3, 0, 4, 7);
  if (k == 3 || k == 4) index |= (1<<4);
  k = SplitSide(cell, 4, 5, 6, 7);
  if (k == 5 || k == 7) index |= (1<<5);
  return index;
}

static void CopyNodes(TET *tet, CELL *cell, int *t)
{
  int i;

  for (i = 0; i < 4; i++) {
    memcpy(tet->x[i], cell->x[t[i]], 3*sizeof(double));
    memcpy(tet->v+i, &cell->v[t[i]],   sizeof(double));
  }
}

void NS_DIM_PREFIX ExtractElement(CELL *cell, double lambda, POLY *poly, int *npoly)
{
  int k;
  TET tet;

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
    k = DecomposeHex(cell);
    CopyNodes(&tet, cell, Hex[k]+ 0);
    ExtractTet(&tet, lambda, poly+0);
    CopyNodes(&tet, cell, Hex[k]+ 4);
    ExtractTet(&tet, lambda, poly+1);
    CopyNodes(&tet, cell, Hex[k]+ 8);
    ExtractTet(&tet, lambda, poly+2);
    CopyNodes(&tet, cell, Hex[k]+12);
    ExtractTet(&tet, lambda, poly+3);
    CopyNodes(&tet, cell, Hex[k]+16);
    ExtractTet(&tet, lambda, poly+4);
    if (Hex[k][20] < 0)
      *npoly = 5;
    else {
      CopyNodes(&tet, cell, Hex[k]+20);
      ExtractTet(&tet, lambda, poly+5);
      *npoly = 6;
    }
    break;
  default :
    *npoly = 0;
  }
}
