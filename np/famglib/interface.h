// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      interface.C													*/
/*																			*/
/* Purpose:   cmg interface													*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   November 97 begin, Stuttgart									*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#ifndef CMG_INTERFACE
#define CMG_INTERFACE

/* RCS_ID
   $Header$
 */

// exported functions

int CMGDeconstructParameter();
int CMGConstructParameter(class CMGParameter *in_parameter);
int CMGConstruct(double *matrix, int *index, int *start, int n, int nl, double *tvA, double *tvB, void **extra);
int CMGPrepare(double *matrix, int *index, int *start, int n, int nl, void **extra);
int CMGSolve(double *rhs, double *unknown, double *defect);
int CMGDeconstruct();
int CMGRepair();


#endif
