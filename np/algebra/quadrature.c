// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      quadrature.c                                                  */
/*                                                                          */
/* Purpose:   quadrature formulas                                               */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   Sep 25 95 begin                                                                           */
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

#include <stdlib.h>

#include "evm.h"
#include "shapes.h"
#include "quadrature.h"
#include "general.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Quadrature formulas                                                      */
/*                                                                          */
/****************************************************************************/

static QUADRATURE Quadrature1D2 = {
  1,
  {{0.5,0,0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {1.0,0,0,0,0,0,0,0,0}
};

static QUADRATURE Quadrature1D4 = {
  2,
  {{ 0.211324865405187, 0, 0}, {0.788675134594813, 0, 0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {0.5,0.5,0,0,0,0,0,0,0}
};

static QUADRATURE Quadrature1D6 = {
  3,
  {{0.5,0,0},
   {0.112701665379258, 0, 0},
   {0.887298334620742, 0, 0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {0.2777777777777777,0.444444444444444,0.2777777777777777,0,0,0,0,0,0}
};

static QUADRATURE Quadrature2D31 = {
  1,
  {{0.333333333333333, 0.333333333333333, 0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {1.0,0,0,0,0,0,0,0,0}
};

static QUADRATURE Quadrature2D32 = {
  3,
  {{0.66666666666666666, 0.16666666666666666, 0},
   {0.16666666666666666, 0.66666666666666666, 0},
   {0.16666666666666666, 0.16666666666666666, 0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {0.333333333333333,     0.333333333333333,      0.333333333333333,0,0,0,0,0,0}
};

static QUADRATURE Quadrature2D33 = {
  4,
  {{0.333333333333333,        0.333333333333333,      0},
   {0.6,0.2,0},
   {0.2,0.6,0},
   {0.2,0.2,0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {-0.5625, 0.520833333333333, 0.520833333333333, 0.520833333333333,
   0,0,0,0,0}
};

static QUADRATURE Quadrature2D34 = {
  6,
  {{0.816847572980459, 0.091576213509771, 0},
   {0.091576213509771, 0.816847572980459, 0},
   {0.091576213509771, 0.091576213509771, 0},
   {0.108103018168070, 0.445948490915965, 0},
   {0.445948490915965, 0.108103018168070, 0},
   {0.445948490915965, 0.445948490915965, 0},
   {0,0,0},{0,0,0},{0,0,0}},
  { 0.109951743655322, 0.109951743655322, 0.109951743655322,
    0.223381589678011, 0.223381589678011, 0.223381589678011, 0,0,0}
};

static QUADRATURE Quadrature2D35 = {
  7,
  {{0.333333333333333, 0.333333333333333,      0},
   {0.059715871789770, 0.470142064105115,  0},
   {0.470142064105115, 0.059715871789770,  0},
   {0.470142064105115, 0.470142064105115,  0},
   {0.797426985368435,     0.101286507323456,      0},
   {0.101286507323456,     0.797426985368435,      0},
   {0.101286507323456,     0.101286507323456,      0},{0,0,0},{0,0,0}},
  { 0.2250300003,
    0.132394152788506,
    0.132394152788506,
    0.132394152788506,
    0.125939180544827,
    0.125939180544827,
    0.125939180544827,0,0}
};

static QUADRATURE Quadrature2D40 = {
  4,
  {{0.0, 0.0, 0.0},
   {1.0, 0.0, 0.0},
   {1.0, 1.0, 0.0},
   {0.0, 1.0, 0.0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0}},
  {0.25,0.25,0.25,0.25,0.0,0.0,0.0,0.0,0.0}
};

static QUADRATURE Quadrature2D42 = {
  4,
  {{0.211324865405187, 0.211324865405187, 0},
   {0.788675134594813, 0.211324865405187, 0},
   {0.211324865405187, 0.788675134594813, 0},
   {0.788675134594813, 0.788675134594813, 0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {0.25,0.25,0.25,0.25,0,0,0,0,0}
};

static QUADRATURE Quadrature3D40 = {
  7,
  {{0.25,0.25,0.25},
   {0.1,0.1,0.4},
   {0.1,0.4,0.1},
   {0.4,0.1,0.1},
   {0.1,0.4,0.4},
   {0.4,0.1,0.4},
   {0.4,0.4,0.1}},
  { 0.142857142857143,0.142857142857143,0.142857142857143,0.142857142857143,
    0.142857142857143,0.142857142857143,0.142857142857143,0,0}
};

static QUADRATURE Quadrature3D41 = {
  1,
  {{0.25,0.25,0.25},
   {0,0,0},{0,0,0},
   {0,0,0},{0,0,0},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {1.0,0,0,0,0,0,0,0,0}
};

static QUADRATURE Quadrature3D42 = {
  4,
  {{0.58541020, 0.13819660, 0.13819660},
   {0.13819660, 0.58541020, 0.13819660},
   {0.13819660, 0.13819660, 0.58541020},
   {0.13819660, 0.13819660, 0.13819660},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {0.25,0.25,0.25,0.25,0,0,0,0,0}
};

static QUADRATURE Quadrature3D43 = {
  5,
  {{0.25,0.25,0.25},
   {0.333333333333333,  0.1666666666666666, 0.1666666666666666},
   {0.1666666666666666, 0.333333333333333,  0.1666666666666666},
   {0.1666666666666666, 0.1666666666666666, 0.333333333333333 },
   {0.1666666666666666, 0.1666666666666666, 0.1666666666666666},
   {0,0,0},{0,0,0},{0,0,0},{0,0,0}},
  {-0.8,0.45,0.45,0.45,0.45,0,0,0,0}
};

static QUADRATURE Quadrature3D52 = {
  8,
  {{0.58541020, 0.72819660, 0.13819660},
   {0.13819660, 0.72819660, 0.13819660},
   {0.13819660, 0.27630920, 0.58541020},
   {0.13819660, 0.27630920, 0.13819660},
   {0.72819660, 0.13819660, 0.13819660},
   {0.72819660, 0.58541020, 0.13819660},
   {0.27630920, 0.13819660, 0.58541020},
   {0.27630920, 0.13819660, 0.13819660},
   {0,0,0}},
  {0.125,0.125,0.125,0.125,
   0.125,0.125,0.125,0.125,0}
};

static QUADRATURE Quadrature3D60 = {
  6,
  {{0.0, 0.0, 0.0},
   {1.0, 0.0, 0.0},
   {0.0, 1.0, 0.0},
   {0.0, 0.0, 1.0},
   {1.0, 0.0, 1.0},
   {0.0, 1.0, 1.0},
   {0,0,0},
   {0,0,0},
   {0,0,0}},
  {0.16666666666666666, 0.16666666666666666,
   0.16666666666666666, 0.16666666666666666,
   0.16666666666666666, 0.16666666666666666, 0, 0, 0}
};

static QUADRATURE Quadrature3D62 = {
  6,
  {{0.66666666666666666, 0.16666666666666666, 0.211324865405187},
   {0.16666666666666666, 0.66666666666666666, 0.211324865405187},
   {0.16666666666666666, 0.16666666666666666, 0.211324865405187},
   {0.66666666666666666, 0.16666666666666666, 0.788675134594813},
   {0.16666666666666666, 0.66666666666666666, 0.788675134594813},
   {0.16666666666666666, 0.16666666666666666, 0.788675134594813},
   {0,0,0},{0,0,0},{0,0,0}},
  {0.16666666666666666, 0.16666666666666666,
   0.16666666666666666, 0.16666666666666666,
   0.16666666666666666, 0.16666666666666666, 0, 0, 0}
};

static QUADRATURE Quadrature3D63 = {
  6,
  {{0.66666666666666666, 0.16666666666666666, 0.211324865405187},
   {0.16666666666666666, 0.66666666666666666, 0.211324865405187},
   {0.16666666666666666, 0.16666666666666666, 0.211324865405187},
   {0.66666666666666666, 0.16666666666666666, 0.788675134594813},
   {0.16666666666666666, 0.66666666666666666, 0.788675134594813},
   {0.16666666666666666, 0.16666666666666666, 0.788675134594813},
   {0,0,0},{0,0,0},{0,0,0}},
  {0.16666666666666666, 0.16666666666666666,
   0.16666666666666666, 0.16666666666666666,
   0.16666666666666666, 0.16666666666666666, 0, 0, 0}
};

static QUADRATURE Quadrature3D80 = {
  8,
  {{0.0, 0.0, 0.0},
   {0.0, 0.0, 1.0},
   {0.0, 1.0, 0.0},
   {0.0, 1.0, 1.0},
   {1.0, 0.0, 0.0},
   {1.0, 0.0, 1.0},
   {1.0, 1.0, 0.0},
   {1.0, 1.0, 1.0},
   {0,0,0}},
  {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.0}
};

static QUADRATURE Quadrature3D82 = {
  8,
  {{0.211324865405187, 0.211324865405187, 0.211324865405187},
   {0.788675134594813, 0.211324865405187, 0.211324865405187},
   {0.211324865405187, 0.788675134594813, 0.211324865405187},
   {0.788675134594813, 0.788675134594813, 0.211324865405187},
   {0.211324865405187, 0.211324865405187, 0.788675134594813},
   {0.788675134594813, 0.211324865405187, 0.788675134594813},
   {0.211324865405187, 0.788675134594813, 0.788675134594813},
   {0.788675134594813, 0.788675134594813, 0.788675134594813},
   {0,0,0}},
  {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0}
};

/****************************************************************************/
/*D
   GetQuadrature - providing a quadrature formula

   SYNOPSIS:
   QUADRATURE *GetQuadrature(INT dim, INT n, INT order);

   PARAMETERS:
   .  dim - dimension of the formular
   .  n - number of corners of the element
   .  order - order of approximation

   DESCRIPTION:
   This function returns a pointer to a quadrature formula.

   EXAMPLE:
   .vb
   QUADRATURE *quadrature;

   if ((quadrature = GetQuadrature(2,3,2)) == NULL)
          return(1);
   sum = 0.0;
   for (l=0; l<Q_NIP(quadrature); l++)
   {
    LOCAL_TO_GLOBAL(n,x,Q_LOCAL(quadrature,l),global);
    (*SourceFunc)(global,val);
    sum += val * Q_WEIGHT(quadrature,l);
   }
   sum = val * AreaOfTriangle;
   .ve

   RETURN VALUE:
   QUADRATURE *
   .n    pointer to quadrature
   .n    NULL if the quadrature formula cannot be found.
   D*/
/****************************************************************************/

QUADRATURE *GetQuadrature(INT dim, INT n, INT order)
{
  switch (dim)
  {
  case 1 :
    switch (order)
    {
    case 2 : return(&Quadrature1D2);
    case 4 : return(&Quadrature1D4);
    default : return(&Quadrature1D6);
    }
  case 2 :
    switch (n)
    {
    case 3 :
      switch (order)
      {
      case 1 : return(&Quadrature2D31);
      case 2 : return(&Quadrature2D32);
      case 3 : return(&Quadrature2D33);
      case 4 : return(&Quadrature2D34);
      default : return(&Quadrature2D35);
      }
    case 4 :
      switch (order)
      {
      case 0 : return(&Quadrature2D40);
      case 2 : return(&Quadrature2D42);
      default : return(&Quadrature2D42);
      }
    }
  case 3 :
    switch (n)
    {
    case 4 :
      switch (order)
      {
      case 0 :  return(&Quadrature3D40);
      case 1 :  return(&Quadrature3D41);
      case 2 :  return(&Quadrature3D42);
      default : return(&Quadrature3D43);
      }
    case 5 : return(&Quadrature3D52);
    case 6 :
      switch (order)
      {
      case 0 :  return(&Quadrature3D60);
      default : return(&Quadrature3D62);
      }
    case 8 :
      switch (order)
      {
      case 0 : return(&Quadrature3D80);
      case 2 : return(&Quadrature3D82);
      default : return(&Quadrature3D82);
      }
    }
  }
  return(NULL);
}

INT GaussPoints(INT dim, INT n, INT order, COORD_VECTOR *x, GAUSS_POINT *gp)
{
  DOUBLE Jdet,area;
  INT ip,nip;
  QUADRATURE *quadrature;

  AREA_OF_REF(n,area);

  if ((quadrature = GetQuadrature(dim,n,order)) == NULL)
    return(1);

  nip = Q_NIP(quadrature);
  for (ip=0; ip<nip; ip++, gp++)
  {
    V_DIM_COPY(Q_LOCAL(quadrature,ip),G_LOCAL(gp));
    LOCAL_TO_GLOBAL (n,x,G_LOCAL(gp),G_GLOBAL(gp));
    INVERSE_TRANSFORMATION(n,x,G_LOCAL(gp),G_JINV(gp),Jdet);
    G_WEIGHT(gp) = Q_WEIGHT(quadrature,ip) * area * ABS(Jdet);
  }

  return(nip);
}
