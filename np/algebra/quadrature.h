// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      quadrature.h                                                  */
/*                                                                          */
/* Purpose:   quadrature formulas (header file)                                 */
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


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __QUADRATURE__
#define __QUADRATURE__

#include "gm.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#ifdef __TWODIM__
#define MAX_INT_POINTS 9
#else
#define MAX_INT_POINTS 27
#endif

#define Q_NIP(p)          ((p)->nip)
#define Q_LOCAL(p,i)      ((DOUBLE*)((const DOUBLE_VECTOR_3D *)((p)->local))[i])
#define Q_WEIGHT(p,i)     ((DOUBLE)((const DOUBLE *)((p)->weight))[i])

#define G_LOCAL(p)        ((p)->local)
#define G_GLOBAL(p)       ((p)->global)
#define G_WEIGHT(p)       ((p)->weight)
#define G_JINV(p)         ((p)->Jinv)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef const struct {
  INT nip;
  const DOUBLE_VECTOR_3D *local;        /* array[nip] for local coordinates */
  const DOUBLE *weight;                         /* array[nip] for weights */
} QUADRATURE;

typedef struct {
  DOUBLE_VECTOR local;
  DOUBLE_VECTOR global;
  DOUBLE weight;
  DOUBLE_VECTOR Jinv[DIM];
} GAUSS_POINT;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

QUADRATURE *GetQuadrature(INT dim, INT n, INT order);
INT GaussPoints(INT dim, INT n, INT order, DOUBLE_VECTOR *x, GAUSS_POINT *gp);

#endif
