// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_header.c													*/
/*																			*/
/* Purpose:   general header for common things (return values, misc..)		*/
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

#ifndef __AMGHEADER__
#define __AMGHEADER__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* general sizes */
#define AMG_NAME_SIZE                   32  /* for names of objects					*/

/* general return values*/
#define AMG_OK                                  0       /* operation succeded					*/
#define AMG_NULL                                NULL /* null pointer							*/
#define AMG_FATAL                               9999 /* fatal error							*/

/* misc macros */
#define AMG_MIN(x,y)            (((x)<(y)) ? (x) : (y))
#define AMG_MAX(x,y)            (((x)>(y)) ? (x) : (y))

#endif
