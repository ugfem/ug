// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  iso.h                                                                                                         */
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
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __ISO__
#define __ISO__

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAXPOLY 6

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef struct {
  int n;
  int order[8];
  double x[8][3];
  double v[8];
} CELL;

typedef struct {
  int n;
  double x[4][3];
} POLY;

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

void ExtractElement(CELL *cell, double lambda, POLY *poly, int *npoly);

END_NAMESPACE

#endif
