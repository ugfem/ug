// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      misc.h														*/
/*																			*/
/* Purpose:   misc definitions												*/
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

#ifndef __CMG_MISC__
#define __CMG_MISC__

#include <iostream.h>
#include <strstream.h>

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* Functions                                                                */
/*                                                                          */
/****************************************************************************/
inline int Abs(int x) {
  return x > 0 ? x : -x;
}
inline double Abs(double x) {
  return x > 0 ? x : -x;
}
inline double Max(double a, double b)
{
  if (a > b) return a;
  else return b;
}
inline double Min(double a, double b)
{
  if (a > b) return b;
  else return a;
}

void CMGError(ostrstream &OutputString);
void CMGWarning(ostrstream &OutputString);
void CMGWrite(ostrstream &OutputString);
double CMGNorm(const int n, const double *v);
void CMGSetVector(const int n, double *v, const double val);
void CMGCopyVector(const int n, double *v1, const double *v2);
void CMGCopyScaledVector(const int n, double *v1, const double *v2, const double factor);
void CMGSubVector(const int n, double *v1, const double *v2);
void CMGAddVector(const int n, double *v1, const double *v2);
void CMGAddVector(const int n, double *v1, const double *v2, const double factor);
void CMGAddVector(const int n, double *v1, const double factor, const double *v2);
void CMGMultVector(const int n, double *v1, const double factor);
void CMGSetSubVector(const int n, double *v1, const double *v2, const double *v3);
double CMGSum(const int n, const double *v1);
double CMGScalProd(const int n, const double *v1, const double *v2);
void CMGEigenVector(int n, double *a, double *b, double *e);


#endif
