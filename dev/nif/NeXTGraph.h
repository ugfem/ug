// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  NeXTGraph.h													*/
/*																			*/
/* Purpose:   handling windows for graphical output                                             */
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/*	History:  September 12, 1997 begin										*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __NEXTGRAPH__
#define __NEXTGRAPH__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DEVICESH__
#include "devices.h"
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

#define GRAPHWIN_MINSIZE                135
#define MAXTITLELENGTH                  128

/********************************************************************************/
/*																				*/
/* data structures exported by the corresponding source file					*/
/*																				*/
/********************************************************************************/

/********************************************************************************/
/*																				*/
/* macros for graph windows                                                                                                     */
/*																				*/
/********************************************************************************/

/********************************************************************************/
/*																				*/
/* data structures exported by the corresponding source file					*/
/*																				*/
/********************************************************************************/

/********************************************************************************/
/*																				*/
/* function declarations														*/
/*																				*/
/********************************************************************************/

OUTPUTDEVICE    *InitNeXTOutputDevice   (void);

INT                      NeXT_ActivateOutput    (WINDOWID win);
INT                      NeXT_UpdateOutput              (WINDOWID win, INT tool);

#endif
