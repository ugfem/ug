// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initparallel.h                                                                                                */
/*																			*/
/* Purpose:   call the init routines of the parallel modules                            */
/*			  (header file)                                                                                                 */
/*																			*/
/* Author:	  Stefan Lang                                                   */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   980604 start                                                  */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

#ifndef __INITPARALLEL_H__
#define __INITPARALLEL_H__

#include "compiler.h"

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* initialisation of the parallel modules */
INT InitParallel (int *argcp, char ***argvp);
INT ExitParallel (void);

#endif
