// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_misc.h													*/
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

#ifndef __FAMG_MISC__
#define __FAMG_MISC__

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

void FAMGError(ostrstream &OutputString);
void FAMGWarning(ostrstream &OutputString);
void FAMGWrite(ostrstream &OutputString);
double FAMGNorm(const int n, const double *v);
void FAMGSetVector(const int n, double *v, const double val);
void FAMGCopyVector(const int n, double *v1, const double *v2);
void FAMGCopyScaledVector(const int n, double *v1, const double *v2, const double factor);
void FAMGSubVector(const int n, double *v1, const double *v2);
void FAMGAddVector(const int n, double *v1, const double *v2);
void FAMGAddVector(const int n, double *v1, const double *v2, const double factor);
void FAMGAddVector(const int n, double *v1, const double factor, const double *v2);
void FAMGMultVector(const int n, double *v1, const double factor);
void FAMGSetSubVector(const int n, double *v1, const double *v2, const double *v3);
double FAMGSum(const int n, const double *v1);
double FAMGScalProd(const int n, const double *v1, const double *v2);
void FAMGEigenVector(int n, double *a, double *b, double *e);

// stuff for timing
#ifdef USE_UG_DS
extern "C"
{
#include "parallel.h"
}
extern double FAMGTimeVar;
void PrintTIME( double time, char *text );
inline void START_SYNC_TIME(void)
{
#ifdef ModelP
  Synchronize();
#endif
  FAMGTimeVar = CURRENT_TIME;
}

inline void END_SYNC_TIME( char *text )
{
  PrintTIME( CURRENT_TIME-FAMGTimeVar, text );
}
#endif
#endif
