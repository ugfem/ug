// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      interface.C													*/
/*																			*/
/* Purpose:   famg interface												*/
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

#ifndef FAMG_INTERFACE
#define FAMG_INTERFACE

/* RCS_ID
   $Header$
 */

// exported functions

int FAMGDeconstructParameter();
int FAMGConstructParameter(class FAMGParameter *in_parameter);
int FAMGConstruct(double *matrix, int *index, int *start, int n, int nl, double *tvA, double *tvB, void **extra);
int FAMGPrepare(double *matrix, int *index, int *start, int n, int nl, void **extra);
int FAMGSolve(double *rhs, double *unknown, double *defect);
int FAMGDeconstruct();
int FAMGRepair();


#endif
