// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   meshtype.hh                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#ifndef MESHTYPE
#define MESHTYPE

class ELEMENT
{
  int np;
  int pnum[4];
  int surfaceind;

public:
  ELEMENT () { np = 0; }
  ELEMENT (int anp) { np = anp; }
  ELEMENT & operator= (const ELEMENT & el2);

  void SetNP (int anp) { np = anp; }
  int NP () const { return np; }
  int & PNum (int i) { return pnum[i-1]; }
  const int & PNum (int i) const { return pnum[i-1]; }
  int & PNumMod (int i) { return pnum[(i-1) % np]; }
  const int & PNumMod (int i) const { return pnum[(i-1) % np]; }
  void SetSurfaceIndex (int si) { surfaceind = si; }
  int SurfaceIndex () const { return surfaceind; }
};


class SURFACE;
#endif
