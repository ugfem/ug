// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MInterface.h													*/
/*																			*/
/* Purpose:   OPENSTEP graphical user interface for ug 3.0                                      */
/*            All functions UG needs for user interaction are				*/
/*            defined here.													*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Interdisziplin"ares Zentrum f"ur Wissenschaftliches Rechnen	*/
/*			  Universit"at Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69210 Heidelberg												*/
/*			  email: Volker.Reichenberger@IWR.Uni-Heidelberg.DE		        */
/*																			*/
/*	History:  June 4, 1999 begin (based on OPENSTEP code)					*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef _MINTERFACE_H_
#define _MINTERFACE_H_

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DEVICESH__
#include "devices.h"
#endif

#import <NSWindow.h>

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

struct graphwindow {
  MGraphicWindow  *theGraphWindow;
  struct graphwindow *next;                             /* linked list of all graph windows     */

  /* save size of viewport */
  INT Global_LL[2];                                             /* global x-coord of plotting region	*/
  INT Global_UR[2];                                             /* global y-coord of plotting region	*/
  INT Local_LL[2];                                              /* local  x-coord of plotting region	*/
  INT Local_UR[2];                                              /* local  y-coord of plotting region	*/

  INT currTool;

  /* windows current point */
  int x;
  int y;

};

typedef struct graphwindow GRAPH_WINDOW;


OUTPUTDEVICE    *InitNeXTOutputDevice   (void);

INT                      NeXT_ActivateOutput    (WINDOWID win);
INT                      NeXT_UpdateOutput              (WINDOWID win, INT tool);

#endif /* _MINTERFACE_H_ */
