// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef ModelP


#include <stdio.h>


#include "parallel.h"



void ddd_test (int mode, MULTIGRID *theMG)
{
  InitCurrMG(theMG);
  BalanceGrid(theMG);
}


#endif /* ModelP */
