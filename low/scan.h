// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  scan.h                                                                                                        */
/*																			*/
/* Purpose:   header file for scanning routines for npinit calls            */
/*																			*/
/* Author:	  Christian Wieners                                                                     */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   December 8, 1996                                                                  */
/*			  low part of former np/udm/scan.c, 15.5.97						*/
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

#ifndef __SCAN__
#define __SCAN__

#include "compiler.h"
#include "heaps.h"

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
/* general numerics defines													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* macros concerned with data descriptors and symbols						*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* macros concerned with solving											*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* structures concerned with symbolic user data management					*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* scanning argument lists                                                  */
INT ReadArgvDOUBLE (const char *name, DOUBLE *a, INT argc, char **argv);
INT ReadArgvINT (const char *name, INT *j, INT argc, char **argv);
INT ReadArgvDOUBLE_INT (const char *name, DOUBLE *a, INT *j, INT argc, char **argv);
INT ReadArgvChar (const char *name, char *buffer, INT argc, char **argv);
INT ReadArgvMEM (const char *name, MEM *mem_size, INT argc, char **argv);
INT ReadArgvOption (const char *name, INT argc, char **argv);
INT ReadArgvPosition (const char *name, INT argc, char **argv, DOUBLE *pos);

#endif
