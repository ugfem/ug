// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   reftrans.hh                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/



#ifndef FILE_REFTRANS
#define FILE_REFTRANS

class referencetransform
{
  VEC3D ex, ey, ez;
  VEC3D exh, eyh, ezh;
  VEC3D ex_h, ey_h, ez_h;
  POINT3D rp;
  double h;

public:

  void Set (const POINT3D & p1, const POINT3D & p2,
            const POINT3D & p3, double ah);

  void ToPlain (const POINT3D & p, POINT3D & pp) const;
  void ToPlain (const ARRAY<POINT3D> & p, ARRAY<POINT3D> & pp) const;
  void FromPlain (const POINT3D & pp, POINT3D & p) const;
};


#endif
