// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// find inner point

#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <new.h>
#include <ctype.h>
#include <string.h>

#include <myadt.hh>

#include <geom/geom2d.hh>
#include <geom/geom3d.hh>
#include <linalg/linalg.hh>

// #include <geom/solid.hh>

#ifdef MYGRAPH
#include <geom/rot3d.hh>
#include <graphics/mygraph.hh>
#endif

#include <meshing/global.hh>



int FindInnerPoint (const ARRAY<Point3d> & points,
                    const ARRAY<Element> & faces,
                    Point3d & p)
{
  int i;
  ARRAY<Vec3d> a;
  ARRAY<double> c;
  Point3d p1, pmin;
  int i1, i2, i3, i4;
  int nf;
  DenseMatrix m(3), inv(3);
  Vector rs(3), x(3);
  double f, fmin, hd,l;

  nf = faces.Size();

  //  testout << "#faces  = " << faces.Size() << endl;
  //  testout << "#points = " << points.Size() << endl;

  a.SetSize (nf);
  c.SetSize (nf);

  for (i = 1; i <= nf; i++)
  {
    p1 = points[faces[i].PNum(1)];
    a[i] = Cross (points[faces[i].PNum(2)] - points[faces[i].PNum(1)],
                  points[faces[i].PNum(3)] - points[faces[i].PNum(1)]);
    a[i] /= a[i].Length();
    c[i] = - (a[i].X() * p1.X() + a[i].Y() * p1.Y() + a[i].Z() * p1.Z());
  }

  fmin = 100;
  pmin = Point3d (0, 0, 0);
  for (i1 = 1; i1 <= nf; i1++)
    for (i2 = i1+1; i2 <= nf; i2++)
      for (i3 = i2+1; i3 <= nf; i3++)
        for (i4 = i3+1; i4 <= nf; i4++)
        {
          m.Elem(1, 1) = a.Get(i1).X() - a.Get(i2).X();
          m.Elem(1, 2) = a.Get(i1).Y() - a.Get(i2).Y();
          m.Elem(1, 3) = a.Get(i1).Z() - a.Get(i2).Z();
          rs.Elem(1) = c.Get(i2) - c.Get(i1);

          m.Elem(2, 1) = a.Get(i1).X() - a.Get(i3).X();
          m.Elem(2, 2) = a.Get(i1).Y() - a.Get(i3).Y();
          m.Elem(2, 3) = a.Get(i1).Z() - a.Get(i3).Z();
          rs.Elem(2) = c.Get(i3) - c.Get(i1);

          m.Elem(3, 1) = a.Get(i1).X() - a.Get(i4).X();
          m.Elem(3, 2) = a.Get(i1).Y() - a.Get(i4).Y();
          m.Elem(3, 3) = a.Get(i1).Z() - a.Get(i4).Z();
          rs.Elem(3) = c.Get(i4) - c.Get(i1);

          /*
                    testout << "i1, i2, i3, i4 = " << i1 << ", " << i2 << ", " << i3 << ", " << i4 << endl;
                    testout << "m = " << m << endl;
                    testout << "rs = " << rs << endl;
                    testout << "det = " << m.Det() << endl;
           */
          l = sqrt(rs.Elem(1)*rs.Elem(1) + rs.Elem(2)*rs.Elem(2) + rs.Elem(3)*rs.Elem(3));

          if ( (fabs (m.Det()) > 1e-10) && (l > 1e-10))
          {
            CalcInverse (m, inv);
            inv.Mult (rs, x);

            //            testout << "x = " << x << endl;

            f = -1e10;
            for (i = 1; i <= nf; i++)
            {
              hd = x.Elem(1) * a.Get(i).X()
                   + x.Elem(2) * a.Get(i).Y()
                   + x.Elem(3) * a.Get(i).Z()
                   + c.Get(i);
              if (hd > f) f = hd;
            }

            if (f < fmin)
            {
              fmin = f;
              pmin.X() = x.Elem(1);
              pmin.Y() = x.Elem(2);
              pmin.Z() = x.Elem(3);
            }
          }
        }

  //  testout << "fmin = " << fmin << endl;
  //  testout << "pmin = " << pmin << endl;

  p = pmin;
  return (fmin < 0) ? 1 : 0;
}
