// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  bdf.h                                                                                                         */
/*																			*/
/* Purpose:   header file for BDF time solver								*/
/*																			*/
/* Author:	  Peter Bastian                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   JAN 10 1997                                                                                       */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __BDF__
#define __BDF__

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

/* create standard LinearSolver num proc type */
INT InitBDFSolver (void);

#endif
