// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   meshtype.cc                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/

#include <meshing/meshtype.hh>

ELEMENT & ELEMENT :: operator= (const ELEMENT & el2)
{
  pnum[0] = el2.pnum[0];
  pnum[1] = el2.pnum[1];
  pnum[2] = el2.pnum[2];
  pnum[3] = el2.pnum[3];
  surfaceind = el2.surfaceind;
  np = el2.np;
  return *this;
}
