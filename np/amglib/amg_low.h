// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_low.c														*/
/*																			*/
/* Purpose:   handlers, memory management for amg							*/
/*																			*/
/* Author:	  Peter Bastian					                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   28 Jan 1996 Begin												*/
/*            02 Apr 1996 new memory allocation strategy					*/
/*            30 Sep 1997 redesign											*/
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

#ifndef __AMGLOW__
#define __AMGLOW__

#include "amg_header.h"

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
/******** typedefs     ******************************************************/
/****************************************************************************/

typedef void (*AMG_PrintFuncPtr)(char *);
typedef void * (*AMG_MallocFuncPtr)(size_t n);

/****************************************************************************/
/******** management functions  *********************************************/
/****************************************************************************/

/* string i/o */
int       AMG_InstallPrintHandler (AMG_PrintFuncPtr print);
int       AMG_Print         (char *s);
int       AMG_RedirectToFile (char *name);
int       AMG_RedirectToScreen (void);

/* sp internal memory handling */
int       AMG_InstallMallocHandler (AMG_MallocFuncPtr mall);
void      *AMG_Malloc  (size_t nbytes);

#endif
