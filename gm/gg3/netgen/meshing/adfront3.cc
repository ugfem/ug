// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   adfront3.cc                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <math.h>
#include <limits.h>

#include <template.hh>
#include <array.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>


#include <meshing/adfront3.hh>



ADFRONT3 :: ADFRONT3 ()
{
  nff = 0;
  vol = 0;
}


/*
   void ADFRONT3 :: Load (char * filename, float & h)
   {
   cout << "Load File: " << filename << endl;

   ifstream listin(filename);

   INDEX linesize, pointsize;
   INDEX i;
   POINT3D p;
   ELEMENT el;

   if (!listin.good())
    {
    MyError ("Input file bad");
    exit (1);
    }

   listin >> h;
   listin >> pointsize;  //Anzahl der Punkte

   for (i = 1; i <= pointsize; i++)
    {
    listin >> p.X() >> p.Y() >> p.Z();
    AddPoint (p);
    }

   listin >> linesize;   //Anzahl der Linien

   for (i = 1; i <= linesize; i++)
    {
    el.SetNP(3);
    listin >> el.PNum(1) >> el.PNum(2) >> el.PNum(3);
    AddFace (el);
    }
   }


   void ADFRONT3 :: Save (char * filename, float h)
   {
   ofstream listout(filename);

   INDEX i;
   //  POINT3D p;
   //  ELEMENT el;
   ARRAY<INDEX> valnum(points.Size());
   INDEX vn;

   vn = 0;
   for (i = 1; i <= points.Size(); i++)
    {
    if (points[i].Valid())
      {
      vn++;
      valnum[i] = vn;
      }
    else
      {
      valnum[i] = 0;
      }
    }

   listout << h << endl;
   listout << vn << endl;

   for (i = 1; i <= points.Size(); i++)
    {
    if (points[i].Valid())
      {
      listout << points[i].p.X() << " "
              << points[i].p.Y() << " "
              << points[i].p.Z() << endl;
      }
    }

   listout << nff << endl;

   for (i = 1; i <= faces.Size(); i++)
    {
    if (faces[i].Valid())
      {
      listout << valnum[faces[i].f.PNum(1)] << " "
              << valnum[faces[i].f.PNum(2)] << " "
              << valnum[faces[i].f.PNum(3)] << endl;
      }
    }
   }
 */


void ADFRONT3 :: GetPoints (ARRAY<POINT3D> & apoints) const
{
  INDEX i;
  for (i = 1; i <= points.Size(); i++)
    apoints.Append (points[i].p);
}



INDEX ADFRONT3 :: AddPoint (const POINT3D & p, INDEX globind)
{
  if (delpointl.Size() != 0)
  {
    INDEX pi = delpointl[delpointl.Size()];
    delpointl.DeleteLast ();

    points[pi] = frontpoint3 (p, globind);
    return pi;
  }
  else
  {
    return points.Append (frontpoint3 (p, globind));
  }
}


INDEX ADFRONT3 :: AddFace (const ELEMENT & aface)
{
  int i, minfn;

  nff++;

  for (i = 1; i <= aface.NP(); i++)
    points[aface.PNum(i)].nlinetopoint++;

  const POINT3D & p1 = points[aface.PNum(1)].p;
  const POINT3D & p2 = points[aface.PNum(2)].p;
  const POINT3D & p3 = points[aface.PNum(3)].p;

  vol += 1.0/6.0 * (p1.X() + p2.X() + p3.X()) *
         ( (p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) -
           (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()) );


  minfn = 1000;
  for (i = 1; i <= aface.NP(); i++)
    if (i == 1 || points[aface.PNum(i)].frontnr < minfn)
      minfn = points[aface.PNum(i)].frontnr;

  for (i = 1; i <= aface.NP(); i++)
    if (points[aface.PNum(i)].frontnr > minfn+1)
      points[aface.PNum(i)].frontnr = minfn+1;

  return faces.Append(frontface (aface));
}


void ADFRONT3 :: DeleteFace (INDEX fi)
{
  int i;
  INDEX pi;

  nff--;

  for (i = 1; i <= faces[fi].f.NP(); i++)
  {
    pi = faces[fi].f.PNum(i);
    points[pi].nlinetopoint--;
    if (points[pi].nlinetopoint == 0)
    {
      points[pi].Invalidate();
      delpointl.Append (pi);
    }
  }

  const POINT3D & p1 = points[faces[fi].f.PNum(1)].p;
  const POINT3D & p2 = points[faces[fi].f.PNum(2)].p;
  const POINT3D & p3 = points[faces[fi].f.PNum(3)].p;

  vol -= 1.0/6.0 * (p1.X() + p2.X() + p3.X()) *
         ( (p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) -
           (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()) );

  faces[fi].Invalidate();
}



void ADFRONT3 :: IncrementClass (INDEX fi)
{
  faces[fi].qualclass++;
}


void ADFRONT3 :: ResetClass (INDEX fi)
{
  if (faces[fi].qualclass > 1)
  {
    faces[fi].qualclass = 1;
    faces[fi].oldfront = 0;
  }
}







int ADFRONT3 :: GetLocals (ARRAY<POINT3D> & locpoints,
                           ARRAY<ELEMENT> & locfaces,       // local index
                           ARRAY<INDEX> & pindex,
                           ARRAY<INDEX> & findex,
                           float xh)
{
  INDEX i, j, k;
  INDEX fstind, pstind;
  char found;
  int minval, hi;
  INDEX pi;
  POINT3D midp, p0;


  minval = INT_MAX;

  for (i = 1; i<= faces.Size(); i++)
  {
    hi = faces.Get(i).qualclass +
         2 * min (points.Get(faces.Get(i).f.PNum(1)).frontnr,
                  points.Get(faces.Get(i).f.PNum(2)).frontnr,
                  points.Get(faces.Get(i).f.PNum(3)).frontnr);

    if (hi < minval || i == 1)
    {
      minval = hi;
      fstind = i;
    }
  }

  pstind = faces[fstind].f.PNum(1);
  p0 = points[pstind].p;

  locfaces.Append(faces[fstind].f);
  findex.Append(fstind);

  for (i = 1; i <= faces.Size(); i++)
  {
    if (faces.Get(i).Valid() && i != fstind)
    {
      midp = Center (points.Get(faces.Get(i).f.PNum(1)).p,
                     points.Get(faces.Get(i).f.PNum(2)).p);

      // Center genauer !!!

      if (Dist (midp, p0) <= xh)
      {
        locfaces.Append(faces.Get(i).f);
        findex.Append(i);
      }
    }
  }

  for (i = 1; i <= locfaces.Size(); i++)
  {
    for (j = 1; j <= locfaces.Get(i).NP(); j++)
    {
      found = 0;
      pi = locfaces.Get(i).PNum(j);

      for (k = 1; k <= pindex.Size() && !found; k++)
      {
        if (pindex.Get(k) == pi)
        {
          locfaces.Elem(i).PNum(j) = k;
          found = 1;
        }
      }

      if (!found)
      {
        pindex.Append (pi);
        locfaces.Elem(i).PNum(j) = locpoints.Append (points.Get(pi).p);
      }
    }
  }

  /*
     for (i = 1; i <= points.Size(); i++)
      if (points.Elem(i).Valid() && Dist (points.Elem(i).p, p0) <= xh)
        {
        found = 0;
        for (k = 1; k <= pindex.Size(); k++)
          if (pindex[k] == i)
            found = 1;

        if (!found)
          {
          pindex.Append (i);
          locpoints.Append (points.Get(i).p);
          }
        }
   */
  return faces.Get(fstind).qualclass;
}



void ADFRONT3 :: GetGroup (int fi,
                           ARRAY<POINT3D> & grouppoints,
                           ARRAY<ELEMENT> & groupelements,
                           ARRAY<INDEX> & pindex,
                           ARRAY<INDEX> & findex
                           ) const
{
  ARRAY<char> pingroup(points.Size());
  INDEX i;
  int j, changed, fused;

  for (i = 1; i <= pingroup.Size(); i++)
    pingroup[i] = 0;
  pingroup[faces[fi].f.PNum(1)] = 1;

  do
  {
    changed = 0;

    for (i = 1; i <= faces.Size(); i++)
      if (faces[i].Valid())
      {
        fused = 0;
        for (j = 1; j <= faces[i].f.NP(); j++)
          if (pingroup[faces[i].f.PNum(j)])
            fused = 1;
        if (fused)
          for (j = 1; j <= faces[i].f.NP(); j++)
            if (!pingroup[faces[i].f.PNum(j)])
            {
              pingroup[faces[i].f.PNum(j)] = 1;
              changed = 1;
            }
      }

  }
  while (changed);


  for (i = 1; i <= points.Size(); i++)
  {
    grouppoints.Append (points[i].p);
    pindex.Append (i);
  }

  for (i = 1; i <= faces.Size(); i++)
    if (faces[i].Valid())
    {
      fused = 0;
      for (j = 1; j <= faces[i].f.NP(); j++)
        if (pingroup[faces[i].f.PNum(j)])
          fused = 1;
      if (fused)
      {
        groupelements.Append (faces[i].f);
        findex.Append (i);
      }
    }
}


void ADFRONT3 :: SetStartFront ()
{
  INDEX i;
  int j;

  for (i = 1; i <= faces.Size(); i++)
    if (faces.Get(i).Valid())
      for (j = 1; j <= faces[i].f.NP(); j++)
        points[faces[i].f.PNum(j)].frontnr = 0;
}

void ADFRONT3 :: Print () const
{
  /*
     INDEX i;

     testout << pointl.Size() << " Points: " << endl;
     for (i = 1; i <= pointl.Size(); i++)
      testout << pointl[i] << endl;
     testout << flush;

     testout << linel.Size() << " Lines: " << endl;
     for (i = 1; i <= linel.Size(); i++)
      testout << linel[i].I1() << " - " << linel[i].I2() << endl;
     testout << flush;
   */
}



int ADFRONT3 :: TestOk () const
{
  INDEX i;

  for (i = 1; i <= points.Size(); i++)
    if (points[i].nlinetopoint == 1 || points[i].nlinetopoint == 2)
      return 0;
  return 1;
}


void ADFRONT3 :: SaveSurface (char * filename, double h)
{
  INDEX i, np, nf;
  int j;
  ofstream outfile("nul");
  ARRAY<INDEX> pointnr;


  pointnr.SetSize (points.Size());
  for (i = 1; i <= pointnr.Size(); i++)
    pointnr[i] = 0;

  np = 0;
  nf = 0;
  for (i = 1; i <= faces.Size(); i++)
    if (faces[i].Valid())
    {
      nf++;
      for (j = 1; j <= faces[i].f.NP(); j++)
        if (pointnr[faces[i].f.PNum(j)] == 0)
        {
          np++;
          pointnr[faces[i].f.PNum(j)] = np;
        }
    }

  outfile << "surfacemesh" << endl;
  outfile << h << endl;
  outfile << np << endl;

  for (j = 1; j <= np; j++)
    for (i = 1; i <= points.Size(); i++)
      if (pointnr[i] == j)
        outfile << points[i].p.X() << " "
                << points[i].p.Y() << " "
                << points[i].p.Z() << endl;

  outfile << nf << endl;
  for (i = 1; i <= faces.Size(); i++)
    if (faces[i].Valid())
      outfile << pointnr[faces[i].f.PNum(1)] << " "
              << pointnr[faces[i].f.PNum(2)] << " "
              << pointnr[faces[i].f.PNum(3)] << endl;
}
