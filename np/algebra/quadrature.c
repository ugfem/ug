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
#include "devices.h"
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

static const DOUBLE_VECTOR_3D Quadrature1D2_l[] = {{0.5, 0.0}};
static const DOUBLE Quadrature1D2_w[] = {1.0};
static QUADRATURE Quadrature1D2 = {1, Quadrature1D2_l, Quadrature1D2_w};

static const DOUBLE_VECTOR_3D Quadrature1D4_l[] = {{ 0.211324865405187, 0.0},
                                                   {0.788675134594813, 0.0}};
static const DOUBLE Quadrature1D4_w[] = {0.5, 0.5};
static QUADRATURE Quadrature1D4 = {2, Quadrature1D4_l, Quadrature1D4_w};

static const DOUBLE_VECTOR_3D Quadrature1D6_l[] = {
  {0.112701665379258, 0.0},
  {0.5, 0.0},
  {0.887298334620742, 0.0}
};
static const DOUBLE Quadrature1D6_w[] = {0.2777777777777777, 0.444444444444444, 0.2777777777777777};
static QUADRATURE Quadrature1D6 = {3, Quadrature1D6_l, Quadrature1D6_w};

static const DOUBLE_VECTOR_3D Quadrature2D31_l[] = {{0.333333333333333, 0.333333333333333}};
static const DOUBLE Quadrature2D31_w[] = {1.0};
static QUADRATURE Quadrature2D31 = {1, Quadrature2D31_l, Quadrature2D31_w};

static const DOUBLE_VECTOR_3D Quadrature2D32_l[] = {{0.66666666666666666, 0.16666666666666666},
                                                    {0.16666666666666666, 0.66666666666666666},
                                                    {0.16666666666666666, 0.16666666666666666}};
static const DOUBLE Quadrature2D32_w[] = {0.333333333333333,    0.333333333333333,      0.333333333333333};
static QUADRATURE Quadrature2D32 = {3, Quadrature2D32_l, Quadrature2D32_w};

static const DOUBLE_VECTOR_3D Quadrature2D33_l[] = {{0.333333333333333, 0.3333333333333333},
                                                    {0.6,0.2},
                                                    {0.2,0.6},
                                                    {0.2,0.2}};
static const DOUBLE Quadrature2D33_w[] = {-0.5625, 0.520833333333333, 0.520833333333333, 0.520833333333333};
static QUADRATURE Quadrature2D33 = {4, Quadrature2D33_l, Quadrature2D33_w};

static const DOUBLE_VECTOR_3D Quadrature2D34_l[] = {{0.816847572980459, 0.091576213509771},
                                                    {0.091576213509771, 0.816847572980459},
                                                    {0.091576213509771, 0.091576213509771},
                                                    {0.108103018168070, 0.445948490915965},
                                                    {0.445948490915965, 0.108103018168070},
                                                    {0.445948490915965, 0.445948490915965}};
static const DOUBLE Quadrature2D34_w[] = {0.109951743655322, 0.109951743655322, 0.109951743655322,
                                          0.223381589678011, 0.223381589678011, 0.223381589678011};
static QUADRATURE Quadrature2D34 = {6, Quadrature2D34_l, Quadrature2D34_w};

static const DOUBLE_VECTOR_3D Quadrature2D35_l[] = {{0.333333333333333, 0.333333333333333},
                                                    {0.059715871789770, 0.470142064105115},
                                                    {0.470142064105115, 0.059715871789770},
                                                    {0.470142064105115, 0.470142064105115},
                                                    {0.797426985368435,   0.101286507323456},
                                                    {0.101286507323456,   0.797426985368435},
                                                    {0.101286507323456,   0.101286507323456}};
static const DOUBLE Quadrature2D35_w[] = {0.2250300003,
                                          0.132394152788506,
                                          0.132394152788506,
                                          0.132394152788506,
                                          0.125939180544827,
                                          0.125939180544827,
                                          0.125939180544827};
static QUADRATURE Quadrature2D35 = {7, Quadrature2D35_l, Quadrature2D35_w};

static const DOUBLE_VECTOR_3D Quadrature2D40_l[] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
static const DOUBLE Quadrature2D40_w[] = {0.25,0.25,0.25,0.25};
static QUADRATURE Quadrature2D40 = {4, Quadrature2D40_l, Quadrature2D40_w};

static const DOUBLE_VECTOR_3D Quadrature2D42_l[] = {{0.211324865405187, 0.211324865405187},
                                                    {0.788675134594813, 0.211324865405187},
                                                    {0.211324865405187, 0.788675134594813},
                                                    {0.788675134594813, 0.788675134594813}};
static const DOUBLE Quadrature2D42_w[] = {0.25,0.25,0.25,0.25};
static QUADRATURE Quadrature2D42 = {4, Quadrature2D42_l, Quadrature2D42_w};

static const DOUBLE_VECTOR_3D Quadrature2D44_l[] = {
  {0.112701665379258, 0.112701665379258},
  {0.5, 0.112701665379258},
  {0.887298334620742, 0.112701665379258},
  {0.112701665379258, 0.5},
  {0.5, 0.5},
  {0.887298334620742, 0.5},
  {0.112701665379258, 0.887298334620742},
  {0.5, 0.887298334620742},
  {0.887298334620742, 0.887298334620742}
};
static const DOUBLE Quadrature2D44_w[] = {
  0.07716049382716,
  0.12345679012346,
  0.07716049382716,
  0.12345679012346,
  0.19753086419753,
  0.12345679012346,
  0.07716049382716,
  0.12345679012346,
  0.07716049382716
};
static QUADRATURE Quadrature2D44 = {9, Quadrature2D44_l, Quadrature2D44_w};

static const DOUBLE_VECTOR_3D Quadrature3D40_l[] = {{0.25,0.25,0.25},
                                                    {0.1,0.1,0.4},
                                                    {0.1,0.4,0.1},
                                                    {0.4,0.1,0.1},
                                                    {0.1,0.4,0.4},
                                                    {0.4,0.1,0.4},
                                                    {0.4,0.4,0.1}};
static const DOUBLE Quadrature3D40_w[] = {0.142857142857143,0.142857142857143,0.142857142857143,0.142857142857143,
                                          0.142857142857143,0.142857142857143,0.142857142857143};
static QUADRATURE Quadrature3D40 = {7, Quadrature3D40_l, Quadrature3D40_w};

static const DOUBLE_VECTOR_3D Quadrature3D41_l[] = {{0.25,0.25,0.25}};
static const DOUBLE Quadrature3D41_w[] = {1.0};
static QUADRATURE Quadrature3D41 = {1, Quadrature3D41_l, Quadrature3D41_w};

static const DOUBLE_VECTOR_3D Quadrature3D42_l[] = {{0.58541020, 0.13819660, 0.13819660},
                                                    {0.13819660, 0.58541020, 0.13819660},
                                                    {0.13819660, 0.13819660, 0.58541020},
                                                    {0.13819660, 0.13819660, 0.13819660}};
static const DOUBLE Quadrature3D42_w[] = {0.25,0.25,0.25,0.25};
static QUADRATURE Quadrature3D42 = {4, Quadrature3D42_l, Quadrature3D42_w};

static const DOUBLE_VECTOR_3D Quadrature3D43_l[] = {
  {0.0,0.0,0.0},
  {1.0,0.0,0.0},
  {0.0,1.0,0.0},
  {1.0,0.0,1.0},
  {0.333333333333,0.333333333333,0.0},
  {0.333333333333,0.0,0.333333333333},
  {0.0,0.333333333333,0.333333333333},
  {0.333333333333,0.333333333333,0.333333333333}
};
static const DOUBLE Quadrature3D43_w[] = {
  0.025,0.025,0.025,0.025,0.225,0.225,0.225,0.225
};
static QUADRATURE Quadrature3D43 = {8, Quadrature3D43_l, Quadrature3D43_w};

/* Strout 1971 */

#define s_1 0.091971078   /* (7 - sq(15) ) / 34 */
#define s_2 0.31979363    /* (7 + sq(15) ) / 34 */
#define t_1 0.72408677    /* (13 + 3*sq(15) ) / 34 */
#define t_2 0.040619117   /* (13 - 3*sq(15) ) / 34 */
#define u   0.056350833   /* (10 - 2*sq(15) ) / 40 */
#define v   0.44364917    /* (10 + 2*sq(15) ) / 40 */
#define A   0.11851852    /* 16 / 135 */
#define B_1 0.071937084   /* (2665 - 14*sq(15) ) / 37800 */
#define B_2 0.069068207   /* (2665 - 14*sq(15) ) / 37800 */
#define C   0.052910053   /* 20 / 378 */

static const DOUBLE_VECTOR_3D Quadrature3D45_l[] = {
  {0.25,0.25,0.25},
  {s_1,s_1,s_1},
  {t_1,s_1,s_1},
  {s_1,t_1,s_1},
  {s_1,s_1,t_1},
  {s_2,s_2,s_2},
  {t_2,s_2,s_2},
  {s_2,t_2,s_2},
  {s_2,s_2,t_2},
  {v,u,u},
  {u,v,u},
  {u,u,v},
  {v,v,u},
  {v,u,v},
  {u,v,v}
};
static const DOUBLE Quadrature3D45_w[] = {
  A,B_1,B_1,B_1,B_1,B_2,B_2,B_2,B_2,C,C,C,C,C,C
};
static QUADRATURE Quadrature3D45 = {15, Quadrature3D45_l, Quadrature3D45_w};

#undef s_1
#undef s_2
#undef t_1
#undef t_2
#undef u
#undef v
#undef A
#undef B_1
#undef B_2
#undef C

static const DOUBLE_VECTOR_3D Quadrature3D52_l[] = {{0.58541020, 0.72819660, 0.13819660},
                                                    {0.13819660, 0.72819660, 0.13819660},
                                                    {0.13819660, 0.27630920, 0.58541020},
                                                    {0.13819660, 0.27630920, 0.13819660},
                                                    {0.72819660, 0.13819660, 0.13819660},
                                                    {0.72819660, 0.58541020, 0.13819660},
                                                    {0.27630920, 0.13819660, 0.58541020},
                                                    {0.27630920, 0.13819660, 0.13819660}};
static const DOUBLE Quadrature3D52_w[] = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};
static QUADRATURE Quadrature3D52 = {8, Quadrature3D52_l, Quadrature3D52_w};

static const DOUBLE_VECTOR_3D Quadrature3D60_l[] = {{0.0, 0.0, 0.0},
                                                    {1.0, 0.0, 0.0},
                                                    {0.0, 1.0, 0.0},
                                                    {0.0, 0.0, 1.0},
                                                    {1.0, 0.0, 1.0},
                                                    {0.0, 1.0, 1.0}};
static const DOUBLE Quadrature3D60_w[] = {0.16666666666666666, 0.16666666666666666,
                                          0.16666666666666666, 0.16666666666666666,
                                          0.16666666666666666, 0.16666666666666666};
static QUADRATURE Quadrature3D60 = {6, Quadrature3D60_l, Quadrature3D60_w};

static const DOUBLE_VECTOR_3D Quadrature3D62_l[] = {{0.66666666666666666, 0.16666666666666666, 0.211324865405187},
                                                    {0.16666666666666666, 0.66666666666666666, 0.211324865405187},
                                                    {0.16666666666666666, 0.16666666666666666, 0.211324865405187},
                                                    {0.66666666666666666, 0.16666666666666666, 0.788675134594813},
                                                    {0.16666666666666666, 0.66666666666666666, 0.788675134594813},
                                                    {0.16666666666666666, 0.16666666666666666, 0.788675134594813}};
static const DOUBLE Quadrature3D62_w[] = {0.16666666666666666, 0.16666666666666666,
                                          0.16666666666666666, 0.16666666666666666,
                                          0.16666666666666666, 0.16666666666666666};
static QUADRATURE Quadrature3D62 = {6, Quadrature3D62_l, Quadrature3D62_w};

static const DOUBLE_VECTOR_3D Quadrature3D63_l[] = {{0.333333333333333, 0.333333333333333, 0.211324865405187},
                                                    {0.6, 0.2, 0.211324865405187},
                                                    {0.2, 0.6, 0.211324865405187},
                                                    {0.2, 0.2, 0.211324865405187},
                                                    {0.333333333333333, 0.333333333333333, 0.788675134594813},
                                                    {0.6, 0.2, 0.788675134594813},
                                                    {0.2, 0.6, 0.788675134594813},
                                                    {0.2, 0.2, 0.788675134594813}};
static const DOUBLE Quadrature3D63_w[] = {-0.28125, 0.2604166666666666, 0.2604166666666666, 0.2604166666666666,
                                          -0.28125, 0.2604166666666666, 0.2604166666666666, 0.2604166666666666};
static QUADRATURE Quadrature3D63 = {8, Quadrature3D63_l, Quadrature3D63_w};

static const DOUBLE_VECTOR_3D Quadrature3D80_l[] = {{0.0, 0.0, 0.0},
                                                    {0.0, 0.0, 1.0},
                                                    {0.0, 1.0, 0.0},
                                                    {0.0, 1.0, 1.0},
                                                    {1.0, 0.0, 0.0},
                                                    {1.0, 0.0, 1.0},
                                                    {1.0, 1.0, 0.0},
                                                    {1.0, 1.0, 1.0}};
static const DOUBLE Quadrature3D80_w[] = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};
static QUADRATURE Quadrature3D80 = {8, Quadrature3D80_l, Quadrature3D80_w};

static const DOUBLE_VECTOR_3D Quadrature3D82_l[] = {{0.211324865405187, 0.211324865405187, 0.211324865405187},
                                                    {0.788675134594813, 0.211324865405187, 0.211324865405187},
                                                    {0.211324865405187, 0.788675134594813, 0.211324865405187},
                                                    {0.788675134594813, 0.788675134594813, 0.211324865405187},
                                                    {0.211324865405187, 0.211324865405187, 0.788675134594813},
                                                    {0.788675134594813, 0.211324865405187, 0.788675134594813},
                                                    {0.211324865405187, 0.788675134594813, 0.788675134594813},
                                                    {0.788675134594813, 0.788675134594813, 0.788675134594813}};
static const DOUBLE Quadrature3D82_w[] = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};
static QUADRATURE Quadrature3D82 = {8, Quadrature3D82_l, Quadrature3D82_w};

#define A_1 0.112701665379258
#define A_2 0.887298334620742
#define W_1 0.277777777777777
#define W_2 0.444444444444444

static const DOUBLE_VECTOR_3D Quadrature3D84_l[] = {
  {A_1,A_1,A_1},
  {0.5,A_1,A_1},
  {A_2,A_1,A_1},
  {A_1,0.5,A_1},
  {0.5,0.5,A_1},
  {A_2,0.5,A_1},
  {A_1,A_2,A_1},
  {0.5,A_2,A_1},
  {A_2,A_2,A_1},
  {A_1,A_1,0.5},
  {0.5,A_1,0.5},
  {A_2,A_1,0.5},
  {A_1,0.5,0.5},
  {0.5,0.5,0.5},
  {A_2,0.5,0.5},
  {A_1,A_2,0.5},
  {0.5,A_2,0.5},
  {A_2,A_2,0.5},
  {A_1,A_1,A_2},
  {0.5,A_1,A_2},
  {A_2,A_1,A_2},
  {A_1,0.5,A_2},
  {0.5,0.5,A_2},
  {A_2,0.5,A_2},
  {A_1,A_2,A_2},
  {0.5,A_2,A_2},
  {A_2,A_2,A_2}
};
static const DOUBLE Quadrature3D84_w[] =
{W_1*W_1*W_1,W_2*W_1*W_1,W_1*W_1*W_1,
 W_1*W_2*W_1,W_2*W_2*W_1,W_1*W_2*W_1,
 W_1*W_1*W_1,W_2*W_1*W_1,W_1*W_1*W_1,
 W_1*W_1*W_2,W_2*W_1*W_2,W_1*W_1*W_2,
 W_1*W_2*W_2,W_2*W_1*W_2,W_1*W_2*W_2,
 W_1*W_1*W_2,W_2*W_1*W_2,W_1*W_1*W_2,
 W_1*W_1*W_1,W_2*W_1*W_1,W_1*W_1*W_1,
 W_1*W_2*W_1,W_2*W_1*W_1,W_1*W_2*W_1,
 W_1*W_1*W_1,W_2*W_1*W_1,W_1*W_1*W_1};
static QUADRATURE Quadrature3D84 = {27, Quadrature3D84_l, Quadrature3D84_w};

#undef A_1
#undef A_2
#undef W_1
#undef W_2


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
      case 1 : return(&Quadrature2D40);
      case 2 : return(&Quadrature2D42);
      case 3 : return(&Quadrature2D42);
      case 4 : return(&Quadrature2D44);
      default : return(&Quadrature2D44);
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
      case 3 :  return(&Quadrature3D43);
      default : return(&Quadrature3D45);
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
      case 1 : return(&Quadrature3D82);
      case 2 : return(&Quadrature3D82);
      default : return(&Quadrature3D84);
      }
    }
  }
  return(NULL);
}

INT GaussPoints(INT dim, INT n, INT order, DOUBLE_VECTOR *x, GAUSS_POINT *gp)
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
