// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  coviseif.h                                                                                                    */
/*																			*/
/* Purpose:   header for Covise<->UG interface                                  */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de					                */
/*																			*/
/* History:   971216 begin                                                                      */
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

#ifndef __COVISEIF_H__
#define __COVISEIF_H__

#include <stdio.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif


/****************************************************************************/
/*																			*/
/*	constants																*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/*	macros																	*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/*	data types																*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* function exported by this module                                                                     */
/*																			*/
/****************************************************************************/


/* initialization and clean up */
INT InitCoviseIF                                (void);
INT ExitCoviseIF              (void);

/* connect */
INT ConnectCovise(MULTIGRID *, char *);

#endif
