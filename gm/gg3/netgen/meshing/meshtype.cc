// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <meshing/meshtype.hh>

Element & Element :: operator= (const Element & el2)
{
  pnum[0] = el2.pnum[0];
  pnum[1] = el2.pnum[1];
  pnum[2] = el2.pnum[2];
  pnum[3] = el2.pnum[3];
  pnum[4] = el2.pnum[4];
  pnum[5] = el2.pnum[5];
  surfaceind = el2.surfaceind;
  np = el2.np;
  return *this;
}
