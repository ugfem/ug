// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_coloring.h												*/
/*																			*/
/* Purpose:   parallel graph coloring functions for FAMG					*/
/*																			*/
/* Author:    Christian Wrobel												*/
/*			  Institut fuer Wissenschaftliches Rechnen						*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69120 Heidelberg												*/
/*			  internet: Christian.Wrobel@iwr.uni-heidelberg.de				*/
/*																			*/
/*																			*/
/* History:   February 99 begin, Stuttgart									*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/
#ifndef __FAMG_COLORING__
#define __FAMG_COLORING__

/* RCS_ID
   $Header$
 */

#ifdef ModelP

#define FAMGColorMaxProcs 512
#define FAMGColorMaxNb FAMGColorMaxProcs

typedef int FAMGColor;

extern FAMGColor FAMGMyColor;

int ConstructColoringGraph( DDD_ATTR grid_attr);
int ConstructColoring( int OrderingFunctionType );

#endif // ModelP

#endif
