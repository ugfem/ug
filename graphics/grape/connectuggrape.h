// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*
 * File:	  connectuggrape.h
 *
 * Purpose:   header for connectuggrape.c
 *
 * Author:	  Monika Wierse
 *                           Universitaet Freiburg
 *                           Institut fuer Angewandte Mathematik
 *                           Hermann--Herder--Str. 10
 *                           D-79104 Freiburg
 *                Klaus Johannsen
 *			  Institut fuer Computeranwendungen III
 *			  Universitaet Stuttgart
 *			  Pfaffenwaldring 27
 *			  70550 Stuttgart
 *                         phone: 0049-(0)711-685-7003
 *			  fax  : 0049-(0)711-685-7000
 *
 * History:   9.10.95
 *
 ****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __GRAPE__
#define __GRAPE__

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

INT InitGrape(void);
int CallGrape (MULTIGRID *theMG);
void usleep (unsigned long time);

#endif
