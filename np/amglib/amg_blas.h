// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_blas.h													*/
/*																			*/
/* Purpose:   BLAS on amg data structure									*/
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
/* History:   01 FEB 1996 Begin												*/
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

#ifndef __AMG_BLAS__
#define __AMG_BLAS__

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

/* BLAS Level 1 */
int    AMG_dset      (AMG_VECTOR *x, double a                   );
int    AMG_randomize (AMG_VECTOR *x                                                             );
int    AMG_dcopy     (AMG_VECTOR *x, AMG_VECTOR *y              );
int    AMG_dscale    (AMG_VECTOR *x, double a                   );
int    AMG_daxpy     (AMG_VECTOR *x, double a,     AMG_VECTOR *y);
double AMG_ddot      (AMG_VECTOR *x, AMG_VECTOR *y              );

/* BLAS Level 2 */
int    AMG_dmatset   (AMG_MATRIX *A, double a                    );
int    AMG_dmatcopy  (AMG_MATRIX *A, AMG_MATRIX *B               );
int    AMG_dmatmul   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_VECTOR *y);
int    AMG_dmatminus (AMG_VECTOR *x, AMG_MATRIX *A, AMG_VECTOR *y);

#endif
