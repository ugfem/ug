// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  graph.h														*/
/*																			*/
/* Purpose:   defines data structure for graph.c							*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Institut fuer Computeranwendungen                                                     */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: klaus@ica3.uni-stuttgart.de                                                 */
/*																			*/
/* History:   8.12.94 begin, ug3-version									*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __GRAPH__
#define __GRAPH__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __WPM__
#include "wpm.h"
#endif

#ifdef __MISC__
#include "misc.h"
#endif

#ifndef __EVM__
#include "evm.h"
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

#define MAX_POINTS_OF_POLY                      32

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* misc. functions */
INT     PrepareGraph                            (const PICTURE *thePicture);

/* low level drawing functions */
void    UgMove                                          (COORD_POINT in);
void    UgDraw                                          (COORD_POINT point);
void    UgLine                                          (COORD_POINT point1, COORD_POINT point2);
void    UgInverseLine                           (COORD_POINT point1, COORD_POINT point2);
void    UgPolyLine                                      (COORD_POINT *points, INT n);
void    UgPolygon                                       (COORD_POINT *points, INT n);
void    UgPolymark                                      (COORD_POINT *points, INT n);
void    UgText                                          (const char *s, INT mode);
void    UgCenteredText                          (COORD_POINT point, const char *s, INT mode);
void    UgClearViewPort                         (void);
void    UgErasePolygon                          (COORD_POINT *points, INT n);
void    UgInversePolygon                        (COORD_POINT *points, INT n);

/* set functions */
void    UgSetColor                                      (long colorIndex);
void    UgSetMarker                             (short index);
void    UgSetMarkerSize                         (short Index);
void    UgSetTextSize                           (short size);
void    UgSetLineWidth                          (short width);

/* miscellenious */
void    UgFlushCash                             (void);

#endif
