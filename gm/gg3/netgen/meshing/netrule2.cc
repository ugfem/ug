// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
// #include <conio.h>

#include <template.hh>
#include <array.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>
#include <meshing/ruler2.hh>


netrule :: netrule ()
{
  name = "";
  quality = 0;
}

/*
   void netrule :: GetFreeArea (ARRAY<Point2d> & afreearea)
   {
   int i;

   afreearea.SetSize (freearea.Size());
   for (i = 1; i <= freearea.Size(); i++)
    afreearea[i] = freearea[i];
   }
 */


void netrule :: SetFreeZoneTransformation (const Vector & u)
{
  int i;

  transfreezone.SetSize (freezone.Size());

  for (i = 1; i <= freezone.Size(); i++)
  {
    transfreezone.Elem(i).X() = freezone.Elem(i).X() + u.Get(2 * i - 1);
    transfreezone.Elem(i).Y() = freezone.Elem(i).Y() + u.Get(2 * i);

    if (i == 1)
    {
      fzmaxx = fzminx = transfreezone.Elem(1).X();
      fzmaxy = fzminy = transfreezone.Elem(1).Y();
    }
    else
    {
      if (transfreezone.Elem(i).X() > fzmaxx)
        fzmaxx = transfreezone.Elem(i).X();
      if (transfreezone.Elem(i).X() < fzminx)
        fzminx = transfreezone.Elem(i).X();
      if (transfreezone.Elem(i).Y() > fzmaxy)
        fzmaxy = transfreezone.Elem(i).Y();
      if (transfreezone.Elem(i).Y() < fzminy)
        fzminy = transfreezone.Elem(i).Y();
    }
  }



  for (i = 1; i <= transfreezone.Size(); i++)
  {
    const Point2d * p1, * p2;
    Vec2d vn;

    p1 = &transfreezone[i];
    p2 = &transfreezone[i % transfreezone.Size() + 1];

    vn.X() =  (p2->Y() - p1->Y());
    vn.Y() = -(p2->X() - p1->X());

    if (vn.Length2() < 1e-10)
    {
      freesetinequ.Set (i, 1, 0);
      freesetinequ.Set (i, 2, 0);
      freesetinequ.Set (i, 3, -1);
    }
    else

    {
      vn = (1/vn.Length()) * vn;

      freesetinequ.Set(i, 1, vn.X());
      freesetinequ.Set(i, 2, vn.Y());
      freesetinequ.Set(i, 3, -(p1->X() * vn.X() + p1->Y() * vn.Y()));
    }
  }
}


int netrule :: IsInFreeZone (const Point2d & p) const
{
  int i;

  if (p.X() < fzminx || p.X() > fzmaxx ||
      p.Y() < fzminy || p.Y() > fzmaxy) return 0;

  for (i = 1; i <= transfreezone.Size(); i++)
  {
    if (freesetinequ.Get(i, 1) * p.X() + freesetinequ.Get(i, 2) * p.Y() +
        freesetinequ.Get(i, 3) > 0) return 0;
  }
  return 1;
}

int netrule :: IsLineInFreeZone (const Point2d & p1, const Point2d & p2) const
{
  int i;
  int left, right, allleft, allright;
  double nx, ny, nl, c;

  if (p1.X() > fzmaxx && p2.X() > fzmaxx ||
      p1.X() < fzminx && p2.X() < fzminx ||
      p1.Y() > fzmaxy && p2.Y() > fzmaxy ||
      p1.Y() < fzminy && p2.Y() < fzminy) return 0;

  for (i = 1; i <= transfreezone.Size(); i++)
  {
    if (freesetinequ.Get(i, 1) * p1.X() + freesetinequ.Get(i, 2) * p1.Y() +
        freesetinequ.Get(i, 3) > -1e-5 &&
        freesetinequ.Get(i, 1) * p2.X() + freesetinequ.Get(i, 2) * p2.Y() +
        freesetinequ.Get(i, 3) > -1e-5
        ) return 0;
  }

  nx =  (p2.Y() - p1.Y());
  ny = -(p2.X() - p1.X());
  nl = sqrt (nx * nx + ny * ny);
  if (nl > 1e-8)
  {
    nx /= nl;
    ny /= nl;
    c = - (p1.X() * nx + p1.Y() * ny);

    allleft = 1;
    allright = 1;

    for (i = 1; i <= transfreezone.Size(); i++)
    {
      left  = transfreezone.Get(i).X() * nx + transfreezone.Get(i).Y() + c <  1e-7;
      right = transfreezone.Get(i).X() * nx + transfreezone.Get(i).Y() + c > -1e-7;

      if (!left) allleft = 0;
      if (!right) allright = 0;
    }
    if (allleft || allright) return 0;
  }

  return 1;
}

int netrule :: ConvexFreeZone () const
{
  int i, n;
  n = transfreezone.Size();
  for (i = 1; i <= n; i++)
  {
    if (! CCW (transfreezone.Get(i), transfreezone.Get(i % n + 1),
               transfreezone.Get( (i+1) % n + 1 ) ) )
      return 0;
  }
  return 1;
}



float netrule :: CalcPointDist (int pi, const Point2d & p) const
{
  float dx = p.X() - points.Get(pi).X();
  float dy = p.Y() - points.Get(pi).Y();
  const threefloat * tf = &tolerances.Get(pi);

  return tf->f1 * dx * dx + tf->f2 * dx * dy + tf->f3 * dy * dy;
}


float netrule :: CalcLineError (int li, const Vec2d & v) const
{
  float dx = v.X() - linevecs.Get(li).X();
  float dy = v.Y() - linevecs.Get(li).Y();

  const threefloat * tf = &linetolerances.Get(li);
  return tf->f1 * dx * dx + tf->f2 * dx * dy + tf->f3 * dy * dy;
}




/*
   int GetNRules ()
   {
   return rules.Size();
   }
 */
