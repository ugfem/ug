// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  UGGraphWindow.h												*/
/*																			*/
/* Purpose:   Interface for UGGraphWindow object.							*/
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

#ifndef _UGGRAPHWINDOW_
#define _UGGRAPHWINDOW_

#import <appkit/appkit.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DEVICESH__
#include "ugdevices.h"
#endif

/* RCS_ID
   $Header$
 */



@interface UGGraphWindow : Window
{
  UGGraphWindow *next;                                  /* linked list of all graph windows     */

  /* drawing parameters */
  short marker_size;                                            /* size of markers in pixels			*/
  short marker_id;                                              /* number of marker                                     */
  short textSize;                                               /* text size							*/

  /* save size of viewport */
  INT Global_LL[2];                                             /* global x-coord of plotting region	*/
  INT Global_UR[2];                                             /* global y-coord of plotting region	*/
  INT Local_LL[2];                                              /* local  x-coord of plotting region	*/
  INT Local_UR[2];                                              /* local  y-coord of plotting region	*/

  INT currTool;

  /* windows current point */
  int x;
  int y;

  short nx,ny;                                                  /* discretization for scroll bars		*/
  short i,j;                                                            /* discrete midpoint coordinates		*/
}

- setUp;

@end

#endif
