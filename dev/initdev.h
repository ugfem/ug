// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initdev.h                                                                                                     */
/*																			*/
/* Purpose:   contains definitions of init functions of all devices             */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   16.12.94 restructured for ug version 3.0						*/
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

#ifndef __INITDEV__
#define __INITDEV__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DEVICESH__
#include "devices.h"
#endif

/****************************************************************************/
/*																			*/
/* init and exit functions called internally by devices.c							*/
/*																			*/
/****************************************************************************/


/* the screen device (same name for Mac, X11 and Remote!) */
OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error);
void ExitScreen (void);


/* metafile device */
INT InitMeta (void);


/* metafile device */
INT InitPostScript (void);

#endif
