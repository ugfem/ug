// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugm2d.h														*/
/*																			*/
/* Purpose:   dimension dependent utilities                                                             */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ica3.uni-stuttgart.de									*/
/*																			*/
/* History:   13.12.94 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UGM2D__
#define __UGM2D__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __GM__
#include "gm.h"
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
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/


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

/* miscellaneous */
INT              CheckOrientation               (INT n, VERTEX **vertices);
INT              QualityElement                 (int n, DOUBLE *xc, DOUBLE *yc, DOUBLE *lmin, DOUBLE *lmax);
EDGE            *CreateAuxEdge                  (GRID *theGrid, NODE *from, NODE *to);
INT              DisposeAuxEdges                (GRID *theGrid);

#endif
