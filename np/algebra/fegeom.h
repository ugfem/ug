// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      fegeom.h                                                          */
/*                                                                          */
/* Purpose:   Finite--Element geometry related data							*/
/*			  dimension independent, general elmement						*/
/*                                                                          */
/* Author:	  Peter Bastian                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  fon: 0049-(0)711-685-7003										*/
/*			  fax: 0049-(0)711-685-7000										*/
/*																			*/
/* History:   06.05.96 begin, ug version 3.2								*/
/*            02.07.96 begin, adapted from FV to FE							*/
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#ifndef __FEGEOM__
#define __FEGEOM__

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "devices.h"
#include "enrol.h"
#include "compiler.h"
#include "misc.h"
#include "gm.h"
#include "ugenv.h"
#include "ugm.h"
#include "algebra.h"
#include "cmdint.h"
#include "commands.h"
#include "helpmsg.h"
#include "evm.h"
#include "domain.h"
#include "quadrature.h"

#ifndef MAXNC
#define MAXNC           MAX_CORNERS_OF_ELEM             /* just to make it more readable*/
#endif

#define MAXCON          (MAXNC*(MAXNC-1))/2             /* need also internal edges !	*/

#ifndef MAXS
#define MAXS            MAX_SIDES_OF_ELEM               /* just to make it more readable*/
#endif

#define MAXGP MAX_INT_POINTS

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

typedef struct {
  DOUBLE weight;                                        /* weight for that gauss point			*/
  DOUBLE local[DIM];                                    /* position of gp in local coordinates	*/
  DOUBLE N[MAXNC];                                      /* value of basis functions at gp		*/
  DOUBLE gradN[MAXNC][DIM];                     /* grad of basis functions at gp		*/
  DOUBLE Jinv[DIM][DIM];                        /* inverse of Jacobian of trafo at gp	*/
  DOUBLE AbsdetJ;                                       /* determinant of jacobian of trafo gp  */
} GaussPoint;

typedef struct {
  INT i,j;                                                      /* edge between node i and j			*/
  DOUBLE edge[DIM];                                     /* edge vector pointing from i to j		*/
  DOUBLE local[DIM];                                    /* mid point in local coords			*/
  DOUBLE Jinv[DIM][DIM];                        /* inverse of Jacobian of trafo at em	*/
  DOUBLE detJ;                                          /* determinant of jacobian of trafo em  */
} Edge;

typedef struct {
  DOUBLE weight;                                        /* weight for gauss points				*/
  DOUBLE local[DIM];                                    /* local coordinates of gauss point		*/
  DOUBLE param[DIM];                                    /* parameter for bgp on patch			*/
  DOUBLE surfel;                                        /* surface element at gauss point		*/
  DOUBLE N[MAXNC];                                      /* basis function of corner at bgp		*/
} BoundaryGaussPoint;

typedef struct {
  INT side;                                         /* side id			                    */
  INT nc;                                                       /* number of corners					*/
  INT corners[MAXNC];                                   /* at most four corners...				*/
  INT nbgp;                                                     /* number of boundary gauss points		*/
  BoundaryGaussPoint bgp[MAXGP];        /* boundary gauss points of this side	*/
} BoundarySide;

typedef struct {
  ELEMENT *e;                                                   /* data for this element				*/
  INT tag;                                                      /* element type							*/

  INT nc;                                                       /* number of corners					*/
  INT ngp;                                                      /* number of Gauss points in quadrature */
  INT ned;                                                      /* number of edges including internal	*/
  INT nbs;                                                      /* number of boundary sides				*/

  DOUBLE co_global[MAXNC][DIM];         /* points in global space, corners      */
  DOUBLE co_local[MAXNC][DIM];          /* points in local space, corners       */
  INT node_property[MAXNC];                     /* subdomain info for corner			*/

  GaussPoint gp[MAXGP];                         /* quantities at gauss points			*/
  Edge ed[MAXCON];                                      /* quantities on edges					*/
  BoundarySide bs[MAXS];                        /* element sides on domain boundary		*/

} FEElementGeometry;                            /* geometry data for a general element	*/

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

INT EvaluateFEGeometry (ELEMENT *e, FEElementGeometry *geo);

#endif
