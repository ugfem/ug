// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_REFTRANS
#define FILE_REFTRANS

class referencetransform
{
  Vec3d ex, ey, ez;
  Vec3d exh, eyh, ezh;
  Vec3d ex_h, ey_h, ez_h;
  Point3d rp;
  double h;

public:

  void Set (const Point3d & p1, const Point3d & p2,
            const Point3d & p3, double ah);

  void ToPlain (const Point3d & p, Point3d & pp) const;
  void ToPlain (const ARRAY<Point3d> & p, ARRAY<Point3d> & pp) const;
  void FromPlain (const Point3d & pp, Point3d & p) const;
};


#endif
