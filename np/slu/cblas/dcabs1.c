// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*******************************************************/
/* This file was modified according to be integrated   */
/* in the ug-software-package. The changes are due to  */
/*                                                     */
/* Klaus Johannsen,                                    */
/* University of Heidelberg,                           */
/* Im Neuenheimer Feld 368,                            */
/* Germany.                                            */
/*******************************************************/

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
 */

#include "f2c.h"

doublereal dcabs1_(doublecomplex *z)
{
  /* >>Start of File<<

         System generated locals */
  doublereal ret_val;
  static doublecomplex equiv_0[1];

  /* Local variables */
#define t ((doublereal *)equiv_0)
#define zz (equiv_0)

  zz->r = z->r, zz->i = z->i;
  ret_val = abs(t[0]) + abs(t[1]);
  return ret_val;
} /* dcabs1_ */

#undef zz
#undef t
