// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_iter.h													*/
/*																			*/
/* Purpose:   iterative methods on the amg data structure					*/
/*																			*/
/* Author:	  Peter Bastian                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   04 FEB 1996 Begin												*/
/*			  01 OKT 1997 redesign											*/
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

#ifndef __AMG_ITER__
#define __AMG_ITER__

#include "amg_sp.h"

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

/* linear iteration kernels */
int     AMG_jac                 (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);
int     AMG_sorf                (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);
int     AMG_sorb                (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);

#endif
