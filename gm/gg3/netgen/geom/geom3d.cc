// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include <template.hh>
#include <array.hh>
#include <geom/geom3d.hh>

Vec3d & Vec3d :: operator/= (double s)
{
  if (s != 0)
  {
    vx /= s;
    vy /= s;
    vz /= s;
  }
  else
  {
    MyError ("Vec3d::operator /=: Divisioin by zero");
  }
  return *this;
}

ostream & operator<<(ostream  & s, const Point3d & p)
{
  return s << "(" << p.px << ", " << p.py << ", " << p.pz << ")";
}

ostream & operator<<(ostream  & s, const Vec3d & v)
{
  return s << "(" << v.vx << ", " << v.vy << ", " << v.vz << ")";
}

double Angle (const Vec3d & v1, const Vec3d & v2)
{
  double co = (v1 * v2) / (v1.Length() * v2.Length());
  if (co > 1) co = 1;
  if (co < -1) co = -1;
  return acos ( co );
}

void Vec3d :: GetNormal (Vec3d & n)
{
  if (fabs (X()) > fabs (Z()))
  {
    n.X() = -Y();
    n.Y() = X();
    n.Z() = 0;
  }
  else
  {
    n.X() = 0;
    n.Y() = Z();
    n.Z() = -n.Y();
  }
  double len = n.Length();
  if (len == 0)
  {
    n.X() = 1;
    n.Y() = n.Z() = 0;
  }
  else
    n /= len;
}

/*
   ostream & operator<<(ostream  & s, const ROTDenseMatrix3D & r)
   {
   return s << "{ (" << r.txx << ", " << r.txy << ", " << r.txz << ") , ("
                    << r.tyx << ", " << r.tyy << ", " << r.tyz << ") , ("
                    << r.tzx << ", " << r.tzy << ", " << r.tzz << ") }";
   }
 */

/*
   Vec3d operator- (const Point3d & p1, const Point3d & p2)
   {
   return Vec3d (p1.X() - p2.X(), p1.Y() - p2.Y(),p1.Z() - p2.Z());
   }

   Point3d operator- (const Point3d & p1, const Vec3d & v)
   {
   return Point3d (p1.X() - v.X(), p1.Y() - v.Y(),p1.Z() - v.Z());
   }

   Point3d operator+ (const Point3d & p1, const Vec3d & v)
   {
   return Point3d (p1.X() + v.X(), p1.Y() + v.Y(),p1.Z() + v.Z());
   }

   Vec3d operator- (const Vec3d & v1, const Vec3d & v2)
   {
   return Vec3d (v1.X() - v2.X(), v1.Y() - v2.Y(),v1.Z() - v2.Z());
   }

   Vec3d operator+ (const Vec3d & v1, const Vec3d & v2)
   {
   return Vec3d (v1.X() + v2.X(), v1.Y() + v2.Y(),v1.Z() + v2.Z());
   }

   Vec3d operator* (double scal, const Vec3d & v)
   {
   return Vec3d (scal * v.X(), scal * v.Y(), scal * v.Z());
   }
 */
/*
   double operator* (const Vec3d & v1, const Vec3d & v2)
   {
   return v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z();
   }

   double Cross (const Vec3d & v1, const Vec3d & v2)
   {
   return v1.X() * v2.Y() - v1.Y() * v2.X();
   }
 */

/*
   void ROTDenseMatrix3D :: CalcRotMat(double ag, double bg, double lg, double size2, Vec3d r)
   {
   size = size2;
   txx=size * ( cos(bg) * cos(lg) );
   txy=size * ( cos(bg) * sin(lg) );
   txz=size * (-sin(bg)           );

   tyx=size * ( sin(ag) * sin(bg) * cos(lg) - cos(ag) * sin(lg) );
   tyy=size * ( sin(ag) * sin(bg) * sin(lg) + cos(ag) * cos(lg) );
   tyz=size * ( sin(ag) * cos(bg)                               );

   tzx=size * ( cos(ag) * sin(bg) * cos(lg) + sin(ag) * sin(lg) );
   tzy=size * ( cos(ag) * sin(bg) * sin(lg) - sin(ag) * cos(lg) );
   tzz=size * ( cos(ag) * cos(bg)                               );

   deltaR=r;
   }
   ROTDenseMatrix3D :: ROTDenseMatrix3D(double ag, double bg, double lg, double size2, Vec3d r)
   {CalcRotMat(ag, bg, lg, size2, r); }

   ROTDenseMatrix3D :: ROTDenseMatrix3D(Vec3d rot2)
   {
   Vec3d r2(0,0,0);
   CalcRotMat(rot2.X(), rot2.Y(), rot2.Z(), 1, r2);
   }

   ROTDenseMatrix3D ROTDenseMatrix3D :: INV()
   {
   ROTDenseMatrix3D rinv(txx/sqr(size),tyx/sqr(size),tzx/sqr(size),
                   txy/sqr(size),tyy/sqr(size),tzy/sqr(size),
                   txz/sqr(size),tyz/sqr(size),tzz/sqr(size),
                   1/size,deltaR);
   return rinv;
   }

   Vec3d operator* (const ROTDenseMatrix3D & r, const Vec3d & v)
   {
   return Vec3d (r.XX() * v.X() + r.XY() * v.Y() + r.XZ() * v.Z(),
                r.YX() * v.X() + r.YY() * v.Y() + r.YZ() * v.Z(),
                r.ZX() * v.X() + r.ZY() * v.Y() + r.ZZ() * v.Z() );
   }

   Point3d operator* (const ROTDenseMatrix3D & r, const Point3d & p)
   {
   return Point3d (r.XX() * p.X() + r.XY() * p.Y() + r.XZ() * p.Z(),
                  r.YX() * p.X() + r.YY() * p.Y() + r.YZ() * p.Z(),
                  r.ZX() * p.X() + r.ZY() * p.Y() + r.ZZ() * p.Z() );
   }
 */







box3d :: box3d ( double aminx, double amaxx,
                 double aminy, double amaxy,
                 double aminz, double amaxz )
{
  minx = aminx; maxx = amaxx;
  miny = aminy; maxy = amaxy;
  minz = aminz; maxz = amaxz;
  CalcDiamCenter ();
}


void box3d :: CalcDiamCenter ()
{
  diam = sqrt( sqr (maxx - minx) + sqr (maxy - miny) + sqr (maxz - minz));

  c.X() = 0.5 * (minx + maxx);
  c.Y() = 0.5 * (miny + maxy);
  c.Z() = 0.5 * (minz + maxz);

  inner = min ( min (maxx - minx, maxy - miny), maxz - minz) / 2;
}


void box3d :: CalcSubBox (int i, box3d & sbox) const
{
  i--;
  if (i & 1)
  {
    sbox.minx = c.X();
    sbox.maxx = maxx;
  }
  else
  {
    sbox.minx = minx;
    sbox.maxx = c.X();
  }
  if (i & 2)
  {
    sbox.miny = c.Y();
    sbox.maxy = maxy;
  }
  else
  {
    sbox.miny = miny;
    sbox.maxy = c.Y();
  }
  if (i & 4)
  {
    sbox.minz = c.Z();
    sbox.maxz = maxz;
  }
  else
  {
    sbox.minz = minz;
    sbox.maxz = c.Z();
  }

  //  sbox.CalcDiamCenter ();

  sbox.c.X() = 0.5 * (sbox.minx + sbox.maxx);
  sbox.c.Y() = 0.5 * (sbox.miny + sbox.maxy);
  sbox.c.Z() = 0.5 * (sbox.minz + sbox.maxz);
  sbox.diam = 0.5 * diam;
  sbox.inner = 0.5 * inner;
}

void box3d :: GetPointNr (int i, Point3d & point) const
{
  i--;
  point.X() = (i & 1) ? maxx : minx;
  point.Y() = (i & 2) ? maxy : miny;
  point.Z() = (i & 4) ? maxz : minz;
}
