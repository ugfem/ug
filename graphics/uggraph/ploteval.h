// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  plotproc.h													*/
/*																			*/
/* Purpose:   data structures for procedures which extract values from the	*/
/*			  multigrid structure											*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*/
/*																			*/
/* History:   09.05.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __PLOTPROC__
#define __PLOTPROC__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __GM__
#include "gm.h"
#endif

#ifndef __ALGEBRA__
#include "algebra.h"
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* data structures exported                                                                                             */
/*																			*/
/****************************************************************************/

/*----------- typedef for functions ----------------------------------------*/

typedef INT (*PreprocessingProcPtr)(MULTIGRID *);
typedef DOUBLE (*ElementPlotProcPtr)(const ELEMENT *,const COORD **,COORD *);
typedef void (*ElementVectorProcPtr)(const ELEMENT *,const COORD **,COORD *,DOUBLE *);
typedef DOUBLE (*MatrixPlotProcPtr)(const MATRIX *);

/*----------- definition of structs ----------------------------------------*/

struct elementvalues {

  /* fields for enironment list variable */
  ENVVAR v;

  PreprocessingProcPtr PreprocessProc;                  /* prepare plot values					*/
  ElementPlotProcPtr PlotProc;                                  /* pointer to corresponding function	*/
} ;

struct elementvector {

  /* fields for enironment list variable */
  ENVVAR v;

  PreprocessingProcPtr PreprocessProc;                  /* prepare plot values					*/
  ElementVectorProcPtr PlotProc;                                /* pointer to corresponding function	*/
  int dimension;                                                                /* dimension of result vector			*/
} ;

struct matrixvalues {

  /* fields for enironment list variable */
  ENVVAR v;

  PreprocessingProcPtr PreprocessProc;                  /* prepare plot values					*/
  MatrixPlotProcPtr PlotProc;                                   /* pointer to corresponding function	*/
} ;

/****************************************************************************/
/*																			*/
/*					typedef for structs                                                                     */
/*																			*/
/****************************************************************************/

typedef struct elementvalues EVALUES ;
typedef struct elementvector EVECTOR ;
typedef struct matrixvalues MVALUES ;

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT              InitPlotproc                                   ();

EVALUES         *CreateElementValuePlotProc     (const char *name, PreprocessingProcPtr PreProc, ElementPlotProcPtr PlotProc);
EVECTOR         *CreateElementVectorPlotProc    (const char *name, PreprocessingProcPtr PreProc, ElementVectorProcPtr PlotProc, INT d);
MVALUES         *CreateMatrixValuePlotProc              (const char *name, PreprocessingProcPtr PreProc, MatrixPlotProcPtr PlotProc);

EVALUES         *GetElementValuePlotProc                (const char *name);
EVECTOR         *GetElementVectorPlotProc               (const char *name);
MVALUES         *GetMatrixValuePlotProc                 (const char *name);

#endif
