// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef MESHTYPE
#define MESHTYPE

class Element
{
  int np;
  int pnum[6];
  int surfaceind;

public:
  Element () { np = 0; }
  Element (int anp) { np = anp; }
  Element & operator= (const Element & el2);

  void SetNP (int anp) { np = anp; }
  int NP () const { return np; }
  int & PNum (int i) { return pnum[i-1]; }
  const int & PNum (int i) const { return pnum[i-1]; }
  int & PNumMod (int i) { return pnum[(i-1) % np]; }
  const int & PNumMod (int i) const { return pnum[(i-1) % np]; }
  void SetSurfaceIndex (int si) { surfaceind = si; }
  int SurfaceIndex () const { return surfaceind; }
};


class Segment
{
public:
  Segment() { p1 = p2 = 0; }
  int p1, p2;         // point indices
  int s1, s2;         // surface indices
  int edgenr;         // edge index
  int invs1, invs2;   // inverse direction for surface 1 resp surface 2
};


class Surface;
#endif
