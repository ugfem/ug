// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      multigrid.h													*/
/*																			*/
/* Purpose:   cmg multigrid classes											*/
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

#ifndef __CMG_MULTIGRID__
#define __CMG_MULTIGRID__

#include "grid.h"

/* RCS_ID
   $Header$
 */

const int CMGMAXGRIDS=32;

class CMGMultiGrid
{
public:
  int Step(int);
  int SGSStep(int);
  void Mult(double *vout, double *vin);
  CMGGrid *GetGrid(int) const;
  int GetN() const;
  int Construct();
  int Deconstruct();
  int Init(const class CMGSystem &);
  int Order();
  int Reorder();
private:
  int n;                    // grids
  CMGGrid *grid[CMGMAXGRIDS];
};

inline CMGGrid *CMGMultiGrid::GetGrid(int i) const {
  return grid[i];
}
inline int CMGMultiGrid::GetN() const {
  return n;
}

#endif
