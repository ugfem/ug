// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <math.h>

#include <myadt.hh>

#include <geom/geom2d.hh>
#include <geom/geom3d.hh>


#include <meshing/adfront3.hh>

ADFRONT3 :: FrontPoint3 :: FrontPoint3 ()
{
  globalindex = 0;
  nfacetopoint = 0;
  frontnr = 1000;
}


ADFRONT3 :: FrontPoint3 :: FrontPoint3 (const Point3d & ap, INDEX agi)
{
  p = ap;
  globalindex = agi;
  nfacetopoint = 0;
  frontnr = 1000;
}


ADFRONT3 :: FrontFace :: FrontFace ()
{
  qualclass = 1;
  oldfront = 0;
}

ADFRONT3 :: FrontFace :: FrontFace (const Element & af)
{
  f = af;
  oldfront = 0;
  qualclass = 1;
}


void ADFRONT3 :: FrontFace :: Invalidate ()
{
  f.PNum(1) = 0;
  oldfront = 0;
  qualclass = 1000;
}


ADFRONT3 :: ADFRONT3 ()
{
  nff = 0;
  vol = 0;
}



void ADFRONT3 :: GetPoints (ARRAY<Point3d> & apoints) const
{
  INDEX i;
  for (i = 1; i <= points.Size(); i++)
    apoints.Append (points[i].P());
}



INDEX ADFRONT3 :: AddPoint (const Point3d & p, INDEX globind)
{
  if (delpointl.Size() != 0)
  {
    INDEX pi = delpointl[delpointl.Size()];
    delpointl.DeleteLast ();

    points[pi] = FrontPoint3 (p, globind);
    return pi;
  }
  else
  {
    return points.Append (FrontPoint3 (p, globind));
  }
}


INDEX ADFRONT3 :: AddFace (const Element & aface)
{
  int i, minfn;

  nff++;

  for (i = 1; i <= aface.NP(); i++)
    points[aface.PNum(i)].AddFace();

  const Point3d & p1 = points[aface.PNum(1)].P();
  const Point3d & p2 = points[aface.PNum(2)].P();
  const Point3d & p3 = points[aface.PNum(3)].P();

  vol += 1.0/6.0 * (p1.X() + p2.X() + p3.X()) *
         ( (p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) -
           (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()) );


  minfn = 1000;
  for (i = 1; i <= aface.NP(); i++)
    if (i == 1 || points[aface.PNum(i)].FrontNr() < minfn)
      minfn = points[aface.PNum(i)].FrontNr();

  for (i = 1; i <= aface.NP(); i++)
    points[aface.PNum(i)].DecFrontNr (minfn+1);

  return faces.Append(FrontFace (aface));
}


void ADFRONT3 :: DeleteFace (INDEX fi)
{
  int i;
  INDEX pi;

  nff--;

  for (i = 1; i <= faces[fi].Face().NP(); i++)
  {
    pi = faces[fi].Face().PNum(i);
    points[pi].RemoveFace();
    if (!points[pi].Valid())
      delpointl.Append (pi);
  }

  const Point3d & p1 = points[faces[fi].Face().PNum(1)].P();
  const Point3d & p2 = points[faces[fi].Face().PNum(2)].P();
  const Point3d & p3 = points[faces[fi].Face().PNum(3)].P();

  vol -= 1.0/6.0 * (p1.X() + p2.X() + p3.X()) *
         ( (p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) -
           (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()) );

  faces[fi].Invalidate();
}



void ADFRONT3 :: IncrementClass (INDEX fi)
{
  faces[fi].IncrementQualClass();
}


void ADFRONT3 :: ResetClass (INDEX fi)
{
  faces[fi].ResetQualClass();
}







int ADFRONT3 :: GetLocals (ARRAY<Point3d> & locpoints,
                           ARRAY<Element> & locfaces,       // local index
                           ARRAY<INDEX> & pindex,
                           ARRAY<INDEX> & findex,
                           float xh)
{
  INDEX i, j, k;
  INDEX fstind, pstind;
  char found;
  int minval, hi;
  INDEX pi;
  Point3d midp, p0;
  double min,f;

  minval = INT_MAX;
  min = INT_MAX;

  // Auswahl des Frontelements nach Schoeberl

  /*  for (i = 1; i<= faces.Size(); i++)
      {
     // ******************** bug **********************************
      hi = faces.Get(i).QualClass();
      hi = faces.Get(i).QualClass() +
           2 * min (points.Get(faces.Get(i).Face().PNum(1)).FrontNr(),
                    points.Get(faces.Get(i).Face().PNum(2)).FrontNr(),
                    points.Get(faces.Get(i).Face().PNum(3)).FrontNr() );

      if (hi < minval || i == 1)
        {
        minval = hi;
        fstind = i;
        }
      }*/

  // Neues Auswahlkriterium

  j = 1;
  fstind = -1;
  do
  {
    for (i = 1; i<= faces.Size(); i++)
    {
      if(faces.Get(i).QualClass()==j)
      {
        // Flaeche = Grundseite * Hoehe / 2
        f = 0.5 * Dist(points.Get(faces.Get(i).Face().PNum(1)).P(),
                       points.Get(faces.Get(i).Face().PNum(2)).P())
            * Dist(points.Get(faces.Get(i).Face().PNum(3)).P(),
                   Center(points.Get(faces.Get(i).Face().PNum(1)).P(),
                          points.Get(faces.Get(i).Face().PNum(2)).P()));
        f = f;
        if (f < min /* || i == 1*/)
        {
          min = f;
          fstind = i;
        }
      }
    }
    j++;
  }
  while(fstind<0);

  pstind = faces[fstind].Face().PNum(1);
  p0 = points[pstind].P();

  locfaces.Append(faces[fstind].Face());
  findex.Append(fstind);

  for (i = 1; i <= faces.Size(); i++)
  {
    if (faces.Get(i).Valid() && i != fstind)
    {
      midp = Center (points.Get(faces.Get(i).Face().PNum(1)).P(),
                     points.Get(faces.Get(i).Face().PNum(2)).P());

      // Center genauer !!!

      if (Dist (midp, p0) <= xh)
      {
        locfaces.Append(faces.Get(i).Face());
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
        locfaces.Elem(i).PNum(j) = locpoints.Append (points.Get(pi).P());
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
  return faces.Get(fstind).QualClass();
}



void ADFRONT3 :: GetGroup (int fi,
                           ARRAY<Point3d> & grouppoints,
                           ARRAY<Element> & groupelements,
                           ARRAY<INDEX> & pindex,
                           ARRAY<INDEX> & findex
                           ) const
{
  ARRAY<char> pingroup(points.Size());
  INDEX i;
  int j, changed, fused;

  for (i = 1; i <= pingroup.Size(); i++)
    pingroup[i] = 0;
  for (j = 1; j <= 3; j++)
    pingroup.Elem (faces.Get(fi).Face().PNum(j)) = 1;

  do
  {
    changed = 0;

    for (i = 1; i <= faces.Size(); i++)
      if (faces.Get(i).Valid())
      {
        fused = 0;
        for (j = 1; j <= faces.Get(i).Face().NP(); j++)
          if (pingroup[faces[i].Face().PNum(j)])
            fused++;

        if (fused >= 2)
          for (j = 1; j <= faces[i].Face().NP(); j++)
            if (!pingroup[faces[i].Face().PNum(j)])
            {
              pingroup[faces[i].Face().PNum(j)] = 1;
              changed = 1;
            }
      }

  }
  while (changed);


  for (i = 1; i <= points.Size(); i++)
  {
    grouppoints.Append (points[i].P());
    pindex.Append (i);
  }

  for (i = 1; i <= faces.Size(); i++)
    if (faces[i].Valid())
    {
      fused = 0;
      for (j = 1; j <= faces[i].Face().NP(); j++)
        if (pingroup[faces[i].Face().PNum(j)])
          fused++;

      if (fused >= 2)
      {
        groupelements.Append (faces[i].Face());
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
      for (j = 1; j <= faces[i].Face().NP(); j++)
        points[faces[i].Face().PNum(j)].DecFrontNr(0);
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
  //  INDEX i;

  //  for (i = 1; i <= points.Size(); i++)
  //    if (points[i].nlinetopoint == 1 || points[i].nlinetopoint == 2)
  //      return 0;
  return 1;
}


void ADFRONT3 :: SaveSurface (char * filename, double h)
{
  INDEX i, np, nf;
  int j;
  ofstream outfile(filename);
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
      for (j = 1; j <= faces[i].Face().NP(); j++)
        if (pointnr[faces[i].Face().PNum(j)] == 0)
        {
          np++;
          pointnr[faces[i].Face().PNum(j)] = np;
        }
    }

  outfile << "surfacemesh" << endl;
  outfile << h << endl;
  outfile << np << endl;

  for (j = 1; j <= np; j++)
    for (i = 1; i <= points.Size(); i++)
      if (pointnr[i] == j)
        outfile << points[i].P().X() << " "
                << points[i].P().Y() << " "
                << points[i].P().Z() << endl;

  outfile << nf << endl;
  for (i = 1; i <= faces.Size(); i++)
    if (faces[i].Valid())
      outfile << pointnr[faces[i].Face().PNum(1)] << " "
              << pointnr[faces[i].Face().PNum(2)] << " "
              << pointnr[faces[i].Face().PNum(3)] << endl;
}
