// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gen.h	                                                                                                        */
/*																			*/
/* Purpose:   general domain declaration	                                                                */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   8 Dec 1999 begin                                                                  */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __GEN_DOMAIN__
#define __GEN_DOMAIN__

#include "compiler.h"
#include "domain.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define  BVP_GENERAL 2

#define MAX_LINES     16
#define MAX_SEGMENTS   8
#define MAX_CORNERS    8
#define MAX_FACES      6
#define MAX_BP         4
#define MAX_DIM        3

typedef struct {

  int id;
  int segment;
  int property;
  double x[MAX_DIM];

} BP;

typedef struct {

  int id;
  int segment;
  int property;
  int n;
  BP *bp[MAX_BP];

} BS;

typedef struct {

  int id;
  int bnd;
  int n;
  int subdomain;
  int property;
  int P[MAX_CORNERS];
  int F[MAX_FACES];
  int S[MAX_FACES];
  int C[MAX_FACES];
  BS *bs[MAX_FACES];

} CELL;

typedef struct {

  int id;
  int n;
  int C;
  int side;
  int S;
  int P[MAX_BP];

} FACE;

typedef struct {

  int id;
  int bnd;
  double x[MAX_DIM];
  BP *bp;

} POINT;

typedef struct {

  int id;
  int L[MAX_LINES];

} XPOINT;

typedef struct {

  int id;
  int S[MAX_SEGMENTS];

} LINE;

typedef struct {

  int id;

} SEGMENT;


typedef INT (*BndCondProcPtr)(int, int, DOUBLE *, DOUBLE *, DOUBLE *);
typedef INT (*BndGeomProcPtr)(int, int, DOUBLE *, DOUBLE *, DOUBLE *);

typedef struct {

  int nP;
  int nBP;
  POINT *P;
  int nC;
  int nBC;
  CELL *C;
  int nF;
  FACE *F;
  int sd;
  int nS;
  SEGMENT *S;
  int nL;
  LINE *L;
  int nX;
  XPOINT *X;

  BndCondProcPtr c;
  BndGeomProcPtr g;

} GEOMETRY;

INT InitGeometry (HEAP *Heap, GEOMETRY *g);

#endif
