// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   netrule3.cc                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#include <stdlib.h>
#include <fstream.h>
#include <iostream.h>
#include <math.h>

#include <template.hh>
#include <array.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>
#include <meshing/ruler3.hh>






vnetrule :: vnetrule ()
{
  name = "";
  quality = 0;
}


/*
   void vnetrule :: GetFreeZone (ARRAY<POINT3D> & afreezone)
   {
   int i;

   afreezone.SetSize (freezone.Size());
   for (i = 1; i <= freezone.Size(); i++)
    afreezone[i] = freezone[i];
   }
 */

int vnetrule :: TestFlag (char flag) const
{
  int i;
  for (i = 1; i <= flags.Size(); i++)
    if (flags.Get(i) == flag) return 1;
  return 0;
}


void vnetrule :: SetFreeZoneTransformation (const VECTOR & u)
{
  int i;
  double nx, ny, nz, v1x, v1y, v1z, v2x, v2y, v2z;
  double nl;
  const threeint * ti;
  int fs;


  transfreezone.SetSize (freezone.Size());

  for (i = 1; i <= freezone.Size(); i++)
  {
    transfreezone.Elem(i).X() = freezone.Elem(i).X() + u.Get(3 * i - 2);
    transfreezone.Elem(i).Y() = freezone.Elem(i).Y() + u.Get(3 * i - 1);
    transfreezone.Elem(i).Z() = freezone.Elem(i).Z() + u.Get(3 * i);

    if (i == 1)
    {
      fzmaxx = fzminx = transfreezone.Elem(1).X();
      fzmaxy = fzminy = transfreezone.Elem(1).Y();
      fzmaxz = fzminz = transfreezone.Elem(1).Z();
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
      if (transfreezone.Elem(i).Z() > fzmaxz)
        fzmaxz = transfreezone.Elem(i).Z();
      if (transfreezone.Elem(i).Z() < fzminz)
        fzminz = transfreezone.Elem(i).Z();
    }
  }


  for (fs = 1; fs <= freesets.Size(); fs++)
  {
    ARRAY<threeint> & freesetfaces = *freefaces.Get(fs);
    MATRIX & freesetinequ = *freefaceinequ.Get(fs);

    for (i = 1; i <= freesetfaces.Size(); i++)
    {
      ti = &freesetfaces.Get(i);
      const POINT3D & p1 = transfreezone.Get(ti->i1);
      const POINT3D & p2 = transfreezone.Get(ti->i2);
      const POINT3D & p3 = transfreezone.Get(ti->i3);

      v1x = p2.X() - p1.X();
      v1y = p2.Y() - p1.Y();
      v1z = p2.Z() - p1.Z();
      v2x = p3.X() - p1.X();
      v2y = p3.Y() - p1.Y();
      v2z = p3.Z() - p1.Z();

      nx = - v2y * v1z + v2z * v1y;
      ny = - v2z * v1x + v2x * v1z;
      nz = - v2x * v1y + v2y * v1x;

      nl = sqrt (nx * nx + ny * ny + nz * nz);

      if (nl < 1e-10)
      {
        freesetinequ.Set(1, 1, 0);
        freesetinequ.Set(1, 1, 0);
        freesetinequ.Set(1, 1, 0);
        freesetinequ.Set(1, 1, -1);
        (*testout) << "setfreezonetransformation: degenerated triangle" << endl << flush;
      }
      else
      {
        nx /= nl;
        ny /= nl;
        nz /= nl;

        freesetinequ.Set(i, 1, nx);
        freesetinequ.Set(i, 2, ny);
        freesetinequ.Set(i, 3, nz);
        freesetinequ.Set(i, 4,
                         -(p1.X() * nx + p1.Y() * ny + p1.Z() * nz));
      }
    }
  }
}


int vnetrule :: ConvexFreeZone () const
{
  int i, j, jj, fs;

  for (fs = 1; fs <= freesets.Size(); fs++)
  {
    ARRAY<threeint> & freesetfaces = *freefaces[fs];
    MATRIX & freesetinequ = *freefaceinequ[fs];
    ARRAY<int> & freeset = *freesets[fs];

    for (i = 1; i <= freesetfaces.Size(); i++)
    {
      for (jj = 1; jj <= freeset.Size(); jj++)
      {
        j = freeset.Get(jj);
        if ( j != freesetfaces.Get(i).i1 && j != freesetfaces.Get(i).i2 &&
             j != freesetfaces.Get(i).i3 &&
             freesetinequ.Get(i, 1) * transfreezone.Get(j).X() +
             freesetinequ.Get(i, 2) * transfreezone.Get(j).Y() +
             freesetinequ.Get(i, 3) * transfreezone.Get(j).Z() +
             freesetinequ.Get(i, 4) > 0 )
        {
          return 0;
        }
      }
    }
  }

  return 1;
}


int vnetrule :: IsInFreeZone (const POINT3D & p) const
{
  int i, fs;
  char inthis;

  if (p.X() < fzminx || p.X() > fzmaxx || p.Y() < fzminy ||
      p.Y() > fzmaxy || p.Z() < fzminz || p.Z() > fzmaxz) return 0;

  for (fs = 1; fs <= freesets.Size(); fs++)
  {
    inthis = 1;
    ARRAY<threeint> & freesetfaces = *freefaces[fs];
    MATRIX & freesetinequ = *freefaceinequ[fs];

    for (i = 1; i <= freesetfaces.Size() && inthis; i++)
    {
      if (freesetinequ.Get(i, 1) * p.X() + freesetinequ.Get(i, 2) * p.Y() +
          freesetinequ.Get(i, 3) * p.Z() + freesetinequ.Get(i, 4) > 0)
        inthis = 0;
    }

    if (inthis) return 1;
  }

  return 0;
}


int vnetrule :: IsTriangleInFreeZone (const POINT3D & p1, const POINT3D & p2, const POINT3D & p3)
{
  int fs;
  int infreeset, cannot = 0;

  for (fs = 1; fs <= freesets.Size(); fs++)
  {
    infreeset = IsTriangleInFreeSet(p1, p2, p3, fs);
    if (infreeset == 1) return 1;
    if (infreeset == -1) cannot = -1;
  }

  return cannot;
}


int vnetrule :: IsTriangleInFreeSet (const POINT3D & p1, const POINT3D & p2,
                                     const POINT3D & p3, int fs)
{
  int i, ii;
  VEC3D n;
  int allleft, allright;
  int hos1, hos2, hos3, os1, os2, os3;
  double hf, lam1, lam2, f, c1, c2, alpha;
  double v1n, v2n, h11, h12, h22, dflam1, dflam2;
  double lam1old, lam2old, fold;
  double hpx, hpy, hpz, v1x, v1y, v1z, v2x, v2y, v2z;
  int act1, act2, act3, it;
  int cntout;
  static ARRAY<int> activefaces;
  int isin;


  ARRAY<threeint> & freesetfaces = *freefaces.Get(fs);
  MATRIX & freesetinequ = *freefaceinequ[fs];


  if (p1.X() <= fzminx && p2.X() <= fzminx && p3.X() <= fzminx ||
      p1.X() >= fzmaxx && p2.X() >= fzmaxx && p3.X() >= fzmaxx ||
      p1.Y() <= fzminy && p2.Y() <= fzminy && p3.Y() <= fzminy ||
      p1.Y() >= fzmaxy && p2.Y() >= fzmaxy && p3.Y() >= fzmaxy ||
      p1.Z() <= fzminz && p2.Z() <= fzminz && p3.Z() <= fzminz ||
      p1.Z() >= fzmaxz && p2.Z() >= fzmaxz && p3.Z() >= fzmaxz) return 0;

  os1 = os2 = os3 = 0;

  activefaces.SetSize(0);

  for (i = 1; i <= freesetfaces.Size(); i++)
  {
    hos1 = freesetinequ.Get(i, 1) * p1.X() +
           freesetinequ.Get(i, 2) * p1.Y() +
           freesetinequ.Get(i, 3) * p1.Z() +
           freesetinequ.Get(i, 4) > -1E-5;

    hos2 = freesetinequ.Get(i, 1) * p2.X() +
           freesetinequ.Get(i, 2) * p2.Y() +
           freesetinequ.Get(i, 3) * p2.Z() +
           freesetinequ.Get(i, 4) > -1E-5;

    hos3 = freesetinequ.Get(i, 1) * p3.X() +
           freesetinequ.Get(i, 2) * p3.Y() +
           freesetinequ.Get(i, 3) * p3.Z() +
           freesetinequ.Get(i, 4) > -1E-5;

    if (hos1 && hos2 && hos3) return 0;

    if (hos1) os1 = 1;
    if (hos2) os2 = 1;
    if (hos3) os3 = 1;

    if (hos1 || hos2 || hos3) activefaces.Append (i);
  }

  if (!os1 || !os2 || !os3) return 1;


  v1x = p2.X() - p1.X();
  v1y = p2.Y() - p1.Y();
  v1z = p2.Z() - p1.Z();

  v2x = p3.X() - p1.X();
  v2y = p3.Y() - p1.Y();
  v2z = p3.Z() - p1.Z();

  n.X() = v1y * v2z - v1z * v2y;
  n.Y() = v1z * v2x - v1x * v2z;
  n.Z() = v1x * v2y - v1y * v2x;
  n /= n.Length();

  allleft = allright = 1;
  for (i = 1; i <= transfreezone.Size() && (allleft || allright); i++)
  {
    const POINT3D & p = transfreezone.Get(i);
    float scal = (p.X() - p1.X()) * n.X() +
                 (p.Y() - p1.Y()) * n.Y() +
                 (p.Z() - p1.Z()) * n.Z();

    if ( scal >  1E-8 ) allleft = 0;
    if ( scal < -1E-8 ) allright = 0;
  }

  if (allleft || allright) return 0;




  lam1old = lam2old = lam1 = lam2 = 1.0 / 3.0;


  //  testout << endl << endl << "Start minimizing" << endl;

  it = 0;
  int minit;
  minit = 1000;
  fold = 1E10;



  while (1)
  {
    it++;

    if (it > 1000) return -1;

    if (lam1 < 0) lam1 = 0;
    if (lam2 < 0) lam2 = 0;
    if (lam1 + lam2 > 1) lam1 = 1 - lam2;

    if (it > minit)
    {
      (*testout) << "it = " << it << endl;
      (*testout) << "lam1/2 = " << lam1 << "  " << lam2 << endl;
    }

    hpx = p1.X() + lam1 * v1x + lam2 * v2x;
    hpy = p1.Y() + lam1 * v1y + lam2 * v2y;
    hpz = p1.Z() + lam1 * v1z + lam2 * v2z;

    f = 0;

    h11 = h12 = h22 = dflam1 = dflam2 = 0;
    cntout = 0;

    isin = 1;

    for (i = 1; i <= activefaces.Size(); i++)
    {
      ii = activefaces[i];

      hf = freesetinequ.Get(ii, 1) * hpx +
           freesetinequ.Get(ii, 2) * hpy +
           freesetinequ.Get(ii, 3) * hpz +
           freesetinequ.Get(ii, 4);

      if (hf > -1E-7) isin = 0;

      hf += 1E-4;
      if (hf > 0)
      {
        f += hf * hf;

        v1n = freesetinequ.Get(ii, 1) * v1x +
              freesetinequ.Get(ii, 2) * v1y +
              freesetinequ.Get(ii, 3) * v1z;
        v2n = freesetinequ.Get(ii, 1) * v2x +
              freesetinequ.Get(ii, 2) * v2y +
              freesetinequ.Get(ii, 3) * v2z;

        h11 += 2 * v1n * v1n;
        h12 += 2 * v1n * v2n;
        h22 += 2 * v2n * v2n;
        dflam1 += 2 * hf * v1n;
        dflam2 += 2 * hf * v2n;
        cntout++;
      }
    }

    if (isin) return 1;

    if (it > minit)
    {
      (*testout) << "f = " << f
                 << "  dfdlam = " << dflam1 << "  " << dflam2 << endl;
      (*testout) << "h = " << h11 << "  " << h12 << "  " << h22 << endl;
      (*testout) << "active: " << cntout << endl;
      (*testout) << "lam1-lam1old = " << (lam1 - lam1old) << endl;
      (*testout) << "lam2-lam2old = " << (lam2 - lam2old) << endl;
    }


    if (f >= fold)
    {
      lam1 = 0.100000000000000 * lam1 + 0.9000000000000000 * lam1old;
      lam2 = 0.100000000000000 * lam2 + 0.9000000000000000 * lam2old;
    }
    else
    {
      lam1old = lam1;
      lam2old = lam2;
      fold = f;


      if (f < 1E-9) return 1;

      h11 += 1E-10;
      h22 += 1E-10;
      c1 = - ( h22 * dflam1 - h12 * dflam2) / (h11 * h22 - h12 * h12);
      c2 = - (-h12 * dflam1 + h11 * dflam2) / (h11 * h22 - h12 * h12);
      alpha = 1;


      if (it > minit)
        (*testout) << "c1/2 = " << c1 << "  " << c2 << endl;

      act1 = lam1 <= 1E-6 && c1 <= 0;
      act2 = lam2 <= 1E-6 && c2 <= 0;
      act3 = lam1 + lam2 >= 1 - 1E-6 && c1 + c2 >= 0;

      if (it > minit)
        (*testout) << "act1,2,3 = " << act1 << act2 << act3 << endl;

      if (act1 && act2 || act1 && act3 || act2 && act3) return 0;

      if (act1)
      {
        c1 = 0;
        c2 = - dflam2 / h22;
      }

      if (act2)
      {
        c1 = - dflam1 / h11;
        c2 = 0;
      }

      if (act3)
      {
        c1 = - (dflam1 - dflam2) / (h11 + h22 - 2 * h12);
        c2 = -c1;
      }

      if (it > minit)
        (*testout) << "c1/2 now = " << c1 << "  " << c2 << endl;


      if (f > 100 * sqrt (sqr (c1) + sqr (c2))) return 0;


      if (lam1 + alpha * c1 < 0 && !act1)
        alpha = -lam1 / c1;
      if (lam2 + alpha * c2 < 0 && !act2)
        alpha = -lam2 / c2;
      if (lam1 + lam2 + alpha * (c1 + c2) > 1 && !act3)
        alpha = (1 - lam1 - lam2) / (c1 + c2);

      if (it > minit)
        (*testout) << "alpha = " << alpha << endl;

      lam1 += alpha * c1;
      lam2 += alpha * c2;
    }
  }
}







float vnetrule :: CalcPointDist (int pi, const POINT3D & p) const
{
  float dx = p.X() - points.Get(pi).X();
  float dy = p.Y() - points.Get(pi).Y();
  float dz = p.Z() - points.Get(pi).Z();

  //  const threefloat * tf = &tolerances.Get(pi);
  //  return tf->f1 * dx * dx + tf->f2 * dx * dy + tf->f3 * dy * dy;
  return tolerances.Get(pi) * (dx * dx + dy * dy + dz * dz);
}


int vnetrule :: TestOk () const
{
  ARRAY<int> cntpused(points.Size());
  ARRAY<int> edge1, edge2;
  ARRAY<int> delf(faces.Size());
  int i, j, k;
  int pi1, pi2;
  int found;

  for (i = 1; i <= cntpused.Size(); i++)
    cntpused[i] = 0;
  for (i = 1; i <= faces.Size(); i++)
    delf[i] = 0;
  for (i = 1; i <= delfaces.Size(); i++)
    delf[delfaces[i]] = 1;


  for (i = 1; i <= faces.Size(); i++)
    if (delf[i] || i > noldf)
      for (j = 1; j <= faces[i].NP(); j++)
        cntpused[faces[i].PNum(j)]++;

  for (i = 1; i <= cntpused.Size(); i++)
    if (cntpused[i] > 0 && cntpused[i] < 2)
    {
      cout << "Fall 1";
      return 0;
    }


  for (i = 1; i <= faces.Size(); i++)
  {
    for (j = 1; j <= faces[i].NP(); j++)
    {
      pi1 = 0; pi2 = 0;
      if (delf[i])
      {
        pi1 = faces[i].PNumMod(j);
        pi2 = faces[i].PNumMod(j+1);
      }
      if (i > noldf)
      {
        pi1 = faces[i].PNumMod(j+1);
        pi2 = faces[i].PNumMod(j);
      }

      found = 0;
      if (pi1)
      {
        for (k = 1; k <= edge1.Size(); k++)
          if (edge1[k] == pi1 && edge2[k] == pi2)
          {
            found = 1;
            edge1.DeleteElement(k);
            edge2.DeleteElement(k);
            k--;
          }
        if (!found)
        {
          edge1.Append (pi2);
          edge2.Append (pi1);
        }
      }
    }
  }

  if (edge1.Size() > 0)
  {
    cout << "Fall 2" << endl;
    return 0;
  }

  /*
     cntpused.SetSize(freezone.Size());
     for (i = 1; i <= cntpused.Size(); i++)
      cntpused[i] = 0;

     for (i = 1; i <= freefaces.Size(); i++)
      {
      cntpused[freefaces[i].i1]++;
      cntpused[freefaces[i].i2]++;
      cntpused[freefaces[i].i3]++;
      }

     for (i = 1; i <= cntpused.Size(); i++)
      if (cntpused[i] < 3)
        {
        cout << "Fall 3" << endl;
        return 0;
        }



     for (i = 1; i <= freefaces.Size(); i++)
      {
      for (j = 1; j <= 3; j++)
        {
        if (j == 1)
          {
          pi1 = freefaces[i].i1;
          pi2 = freefaces[i].i2;
          }
        if (j == 2)
          {
          pi1 = freefaces[i].i2;
          pi2 = freefaces[i].i3;
          }
        if (j == 3)
          {
          pi1 = freefaces[i].i3;
          pi2 = freefaces[i].i1;
          }

        found = 0;
        for (k = 1; k <= edge1.Size(); k++)
          if (edge1[k] == pi1 && edge2[k] == pi2)
            {
            found = 1;
            edge1.DeleteElement(k);
            edge2.DeleteElement(k);
            k--;
            }

        if (!found)
          {
          edge1.Append (pi2);
          edge2.Append (pi1);
          }
        }
      }

     if (edge1.Size() > 0)
      {
      cout << "Fall 4" << endl;
      return 0;
      }
   */
  return 1;
}


int vnetrule :: IsDelFace (int fn) const
{
  int i;
  for (i = 1; i <= GetNDelF(); i++)
    if (GetDelFace(i) == fn) return 1;
  return 0;
}
