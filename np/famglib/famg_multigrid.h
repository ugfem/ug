// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_multigrid.h												*/
/*																			*/
/* Purpose:   famg multigrid classes										*/
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

#ifndef __FAMG_MULTIGRID__
#define __FAMG_MULTIGRID__

#include "famg_grid.h"

/* RCS_ID
   $Header$
 */

const int FAMGMAXGRIDS=32;

class FAMGMultiGrid
{
public:
  int Step(int);
  int SGSStep(int);
  // weg void Mult(double *vout, double *vin);
  FAMGGrid *GetGrid(int) const;
  int GetN() const;
  int Construct();
  int Deconstruct();
  int Init(const class FAMGSystem &);
  int Order();
  int Reorder();
private:
  int n;                    // grids
  FAMGGrid *grid[FAMGMAXGRIDS];
};

inline FAMGGrid *FAMGMultiGrid::GetGrid(int i) const {
  return grid[i];
}
inline int FAMGMultiGrid::GetN() const {
  return n;
}

#endif
