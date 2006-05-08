// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  enewton.h                                                                                                     */
/*																			*/
/* Purpose:   definition of the newton num proc with an nonlinear extension	*/
/*																			*/
/* Author:	  Klaus Johannsen                                                                                       */
/*			  IWR/TS								                                                */
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheiner Feld 368										*/
/*			  69120 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   August 14, 2000                                                                           */
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

#ifndef __ENEWTON__
#define __ENEWTON__

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

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* create extended newton num proc type */
INT InitENewtonSolver();

END_UGDIM_NAMESPACE

#endif
