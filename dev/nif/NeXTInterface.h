// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  NeXTInterface.h												*/
/*																			*/
/* Purpose:   NEXTSTEP graphical user interface for ug 3.0                                      */
/*            All functions UG needs for its graphical output are			*/
/*            defined here.													*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/*	History:  January 16, 1998 begin										*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef _NEXTINTERFACE_H_
#define _NEXTINTERFACE_H_

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DEVICESH__
#include "ugdevices.h"
#endif

#import <appkit/Window.h>

/* RCS_ID
   $Header$
 */

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
  Window  *theWindow;
  View    *theView;
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

#define MENU_BAR                        39      /* the menu-bar has a hight of 39 pixels	*/
#define MARGIN_TO_SCREEN        2       /* leave a margin between window and screen */
#define TITLE_BAR                       19      /* height of the title bar of a window		*/
#define SCROLL_BAR                      15      /* width of h/v scrollbars					*/

#define NO_POS_CHANGE           0       /* drag and grow return code				*/
#define POS_CHANGE                      1       /* drag and grow return code				*/


#endif /* _NEXTINTERFACE_H_ */
