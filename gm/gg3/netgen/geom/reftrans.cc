// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <new.h>

#include <template.hh>
#include <array.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>
// #include <linalg/linalg.hh>

#include <geom/reftrans.hh>



void referencetransform :: Set (const Point3d & p1, const Point3d & p2,
                                const Point3d & p3, double ah)
{
  ex = p2 - p1;
  ex /= ex.Length();
  ey = p3 - p1;
  ey -= (ex * ey) * ex;
  ey /= ey.Length();
  ez = Cross (ex, ey);
  rp = p1;
  h = ah;

  exh = ah * ex;
  eyh = ah * ey;
  ezh = ah * ez;
  ah = 1 / ah;
  ex_h = ah * ex;
  ey_h = ah * ey;
  ez_h = ah * ez;
}

void referencetransform :: ToPlain (const Point3d & p, Point3d & pp) const
{
  Vec3d v;
  v = p - rp;
  pp.X() = (ex_h * v);
  pp.Y() = (ey_h * v);
  pp.Z() = (ez_h * v);
}

void referencetransform :: ToPlain (const ARRAY<Point3d> & p,
                                    ARRAY<Point3d> & pp) const
{
  Vec3d v;
  int i;

  pp.SetSize (p.Size());
  for (i = 1; i <= p.Size(); i++)
  {
    v = p.Get(i) - rp;
    pp.Elem(i).X() = (ex_h * v);
    pp.Elem(i).Y() = (ey_h * v);
    pp.Elem(i).Z() = (ez_h * v);
  }
}

void referencetransform :: FromPlain (const Point3d & pp, Point3d & p) const
{
  Vec3d v;
  //  v = (h * pp.X()) * ex + (h * pp.Y()) * ey + (h * pp.Z()) * ez;
  //  p = rp + v;
  v.X() = pp.X() * exh.X() + pp.Y() * eyh.X() + pp.Z() * ezh.X();
  v.Y() = pp.X() * exh.Y() + pp.Y() * eyh.Y() + pp.Z() * ezh.Y();
  v.Z() = pp.X() * exh.Z() + pp.Y() * eyh.Z() + pp.Z() * ezh.Z();
  p.X() = rp.X() + v.X();
  p.Y() = rp.Y() + v.Y();
  p.Z() = rp.Z() + v.Z();
}
