// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  bullet.h                                                                                                      */
/*																			*/
/* Purpose:   The bullet plotter rasterizes lines and polygons in a local   */
/*            pixel buffer using -- in 3D -- the z buffer algorithm. Main   */
/*            advantage is that rasterization can be done in parallel and   */
/*            that merging the buffers is a simple reduction operation.     */
/*            Besides, the output device can flush the pixels all at once.  */
/*            This speeds up a remote X connection considerably.            */
/*																			*/
/* Author:	  Michael Lampe                                                                                                 */
/*			  Institut fuer Computeranwendungen                                                     */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de                                                    */
/*																			*/
/* History:   24.2.98 begin, ug3-version									*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __BULLET__
#define __BULLET__

#include "compiler.h"
#include "wpm.h"

/****************************************************************************/
/*																			*/
/* constant declarations                                                                        */
/*																			*/
/****************************************************************************/

#define BULLET_OK             0
#define BULLET_CANT           1
#define BULLET_NOMEM          2

/****************************************************************************/
/*																			*/
/*	macros																	*/
/*																			*/
/****************************************************************************/

#define ZTYP         FLOAT             /* type for z buffer                 */
#define ZEPS         (5.0*FLT_EPSILON) /* eps for ZTYP                      */
#define FAR_AWAY     (-FLT_MAX)        /* a large negative number from ZTYP */

/****************************************************************************/
/*																			*/
/*  exported variables	                                                                                                */
/*																			*/
/****************************************************************************/

extern INT BulletDim;

/****************************************************************************/
/*																			*/
/*  exported functions	                                                                                                */
/*																			*/
/****************************************************************************/

INT BulletOpen(PICTURE *picture, DOUBLE factor);
void BulletClose(void);
void BulletPlot(void);
void BulletLine(DOUBLE *point1, DOUBLE *point2, long color);
void BulletPolyLine(DOUBLE *points, INT nb, long color);
void BulletPolygon(DOUBLE *points, INT nb, DOUBLE intensity, long color);

#endif
