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

int pf;

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
  prism_flag = pf;
  pf = -1;
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


INDEX ADFRONT3 :: AddFace (const Element & aface, int prism_flag)
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

  if(aface.NP()==4)
    vol += 1.0/6.0 * (p1.X() + p3.X() + points[aface.PNum(4)].P().X()) *
           ( (p3.Y() - p1.Y()) * (points[aface.PNum(4)].P().Z() - p1.Z()) -
             (p3.Z() - p1.Z()) * (points[aface.PNum(4)].P().Y() - p1.Y()) );


  minfn = 1000;
  for (i = 1; i <= aface.NP(); i++)
    if (i == 1 || points[aface.PNum(i)].FrontNr() < minfn)
      minfn = points[aface.PNum(i)].FrontNr();

  for (i = 1; i <= aface.NP(); i++)
    points[aface.PNum(i)].DecFrontNr (minfn+1);

  pf = prism_flag;
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

  if(faces[fi].Face().NP()==4)
    vol -= 1.0/6.0 * (p1.X() + p3.X() + points[faces[fi].Face().PNum(4)].P().X()) *
           ( (p3.Y() - p1.Y()) * (points[faces[fi].Face().PNum(4)].P().Z() - p1.Z()) -
             (p3.Z() - p1.Z()) * (points[faces[fi].Face().PNum(4)].P().Y() - p1.Y()) );

  faces[fi].Invalidate();
}



void ADFRONT3 :: IncrementClass (INDEX fi)
{
  faces[fi].IncrementQualClass();
}

void ADFRONT3 :: ResetPrism (INDEX fi)
{
  faces[fi].ResetPrismFlag();
}


void ADFRONT3 :: ResetClass (INDEX fi)
{
  faces[fi].ResetQualClass();
}



int ADFRONT3 :: GetLocals_Tetrahedra_new (ARRAY<Point3d> & locpoints,
                                          ARRAY<Element> & locfaces, // local index
                                          ARRAY<INDEX> & pindex,
                                          ARRAY<INDEX> & findex,
                                          float xh,
                                          float relh,
                                          INDEX& facesplit,
                                          int dummy)
{
  INDEX i, j;
  INDEX fstind, pstind;
  int hi;
  INDEX pi;
  Point3d midp, p0;
  static ARRAY<int> invpindex;
  double min,f;

  static ARRAY<Element> locfaces2;           //all local faces in radius xh
  static ARRAY<int> locfaces3;           // all faces in outer radius relh
  static ARRAY<INDEX> findex2;

  locfaces2.SetSize(0);
  locfaces3.SetSize(0);
  findex2.SetSize(0);

  static int minval = -1;
  static int lasti = 0;

  if(dummy==0)
  {
    if (faces.Size() > 2 * nff)
    {
      // compress facelist

      hi = 0;
      for (i = 1; i <= faces.Size(); i++)
        if (faces.Get(i).Valid())
        {
          hi++;
          faces.Elem(hi) = faces.Get(i);
        }

      faces.SetSize (nff);
      lasti = 0;
    }

    fstind = 0;

    for (i = lasti+1; i <= faces.Size() && !fstind; i++)
      if (faces.Elem(i).Valid())
      {
        hi = faces.Get(i).QualClass() +
             points.Get(faces.Get(i).Face().PNum(1)).FrontNr() +
             points.Get(faces.Get(i).Face().PNum(2)).FrontNr() +
             points.Get(faces.Get(i).Face().PNum(3)).FrontNr();

        if (hi <= minval)
        {
          minval = hi;
          fstind = i;
          lasti = fstind;
        }
      }

    if (!fstind)
    {
      minval = INT_MAX;
      for (i = 1; i <= faces.Size(); i++)
        if (faces.Elem(i).Valid())
        {
          hi = faces.Get(i).QualClass() +
               points.Get(faces.Get(i).Face().PNum(1)).FrontNr() +
               points.Get(faces.Get(i).Face().PNum(2)).FrontNr() +
               points.Get(faces.Get(i).Face().PNum(3)).FrontNr();

          if (hi <= minval)
          {
            minval = hi;
            fstind = i;
            lasti = 0;
          }
        }
    }
  }
  else
  {
    min = INT_MAX;
    j = 1;
    fstind = -1;
    do
    {
      for (i = 1; i<= faces.Size(); i++)
      {
        if((faces.Get(i).QualClass()==j))
        {
          // Flaeche = Grundseite * Hoehe / 2
          f = 0.5 * Dist(points.Get(faces.Get(i).Face().PNum(1)).P(),
                         points.Get(faces.Get(i).Face().PNum(2)).P())
              * Dist(points.Get(faces.Get(i).Face().PNum(3)).P(),
                     Center(points.Get(faces.Get(i).Face().PNum(1)).P(),
                            points.Get(faces.Get(i).Face().PNum(2)).P()));
          f = f;
          if (f < min)
          {
            min = f;
            fstind = i;
          }
        }
      }
      j++;
    }
    while(fstind<0);
  }

  pstind = faces[fstind].Face().PNum(1);
  p0 = points[pstind].P();

  locfaces2.Append(faces[fstind].Face());
  findex2.Append(fstind);

  if (0)
  {}
  else
  {
    for (i = 1; i <= faces.Size(); i++)
    {
      const Element & face = faces.Get(i).Face();
      if (faces.Get(i).Valid() && i != fstind)
      {
        const Point3d & p1 = points.Get(face.PNum(1)).P();
        const Point3d & p2 = points.Get(face.PNum(2)).P();
        const Point3d & p3 = points.Get(face.PNum(3)).P();

        //	  midp = Center (p1, p2, p3);
        midp.X() = (1.0/3.0) * (p1.X() + p2.X() + p3.X());
        midp.Y() = (1.0/3.0) * (p1.Y() + p2.Y() + p3.Y());
        midp.Z() = (1.0/3.0) * (p1.Z() + p2.Z() + p3.Z());


        //	  if (Dist2 (midp, p0) <= xh*xh)
        if ( (midp.X()-p0.X()) * (midp.X()-p0.X()) +
             (midp.Y()-p0.Y()) * (midp.Y()-p0.Y()) +
             (midp.Z()-p0.Z()) * (midp.Z()-p0.Z())   <= xh * xh)
        {
          locfaces2.Append(faces.Get(i).Face());
          findex2.Append(i);
        }
      }
    }
  }

  //local faces for inner radius:
  for (i = 1; i <= locfaces2.Size(); i++)
  {
    /*
       midp = Center (points.Get(locfaces2.Get(i).PNum(1)).P(),
       points.Get(locfaces2.Get(i).PNum(2)).P());
       if (Dist (midp, p0) <= relh || i==1)
     */

    const Element & face = locfaces2.Get(i);
    const Point3d & p1 = points.Get(face.PNum(1)).P();
    const Point3d & p2 = points.Get(face.PNum(2)).P();
    const Point3d & p3 = points.Get(face.PNum(3)).P();

    midp.X() = (1.0/3.0) * (p1.X() + p2.X() + p3.X());
    midp.Y() = (1.0/3.0) * (p1.Y() + p2.Y() + p3.Y());
    midp.Z() = (1.0/3.0) * (p1.Z() + p2.Z() + p3.Z());

    if ( (midp.X()-p0.X()) * (midp.X()-p0.X()) +
         (midp.Y()-p0.Y()) * (midp.Y()-p0.Y()) +
         (midp.Z()-p0.Z()) * (midp.Z()-p0.Z())   <= relh * relh || i == 1)
    {
      locfaces.Append(locfaces2.Get(i));
      findex.Append(findex2.Get(i));
    }
    else
      locfaces3.Append (i);
  }

  facesplit=locfaces.Size();


  //local faces for outer radius:
  for (i = 1; i <= locfaces3.Size(); i++)
  {
    locfaces.Append (locfaces2.Get(locfaces3.Get(i)));
    findex.Append (findex2.Get(locfaces3.Get(i)));
  }


  invpindex.SetSize (points.Size());
  for (i = 1; i <= points.Size(); i++)
    invpindex.Elem(i) = 0;

  for (i = 1; i <= locfaces.Size(); i++)
  {
    for (j = 1; j <= 3; j++)
    {
      pi = locfaces.Get(i).PNum(j);
      if (invpindex.Get(pi) == 0)
      {
        pindex.Append (pi);
        invpindex.Elem(pi) = pindex.Size();
        locfaces.Elem(i).PNum(j) = locpoints.Append (points.Get(pi).P());
      }
      else
        locfaces.Elem(i).PNum(j) = invpindex.Get(pi);

    }
  }

  return faces.Get(fstind).QualClass();
}


int ADFRONT3 :: GetLocals_Tetrahedra (ARRAY<Point3d> & locpoints,
                                      ARRAY<Element> & locfaces, // local index
                                      ARRAY<INDEX> & pindex,
                                      ARRAY<INDEX> & findex,
                                      float xh)
{
  return(0);
}

int ADFRONT3 :: GetLocals_Prism(ARRAY<Point3d> & locpoints,
                                ARRAY<Element> & locfaces,
                                ARRAY<INDEX> & pindex,
                                ARRAY<INDEX> & findex,
                                float xh,
                                ARRAY<int> & prism_flags)
{
  INDEX i, j, k;
  INDEX fstind;
  int minval, hi,found;
  INDEX pi;
  Point3d midp,p0,p1,p2,p3,p4;

  minval = INT_MAX;
  fstind = -1;
  for (i = 1; i<= faces.Size(); i++)
    if(faces.Get(i).Valid())
      if(faces.Get(i).PrismFlag()!=-1)
        if(faces.Get(i).Face().NP()==3)
        {
          hi = faces.Get(i).QualClass();
          if (hi < minval || i == 1)
          {
            minval = hi;
            fstind = i;
          }
        }
  if(fstind==-1)
    return(-1);
  else
  {
    p1 =  points[faces[fstind].Face().PNum(1)].P();
    p2 =  points[faces[fstind].Face().PNum(2)].P();
    p3 =  points[faces[fstind].Face().PNum(3)].P();
    p0.X() = (p1.X()+p2.X()+p3.X())/3;
    p0.Y() = (p1.Y()+p2.Y()+p3.Y())/3;
    p0.Z() = (p1.Z()+p2.Z()+p3.Z())/3;

    locfaces.Append(faces[fstind].Face());
    findex.Append(fstind);
    prism_flags.Append(faces.Get(fstind).PrismFlag());

    for (i = 1; i <= faces.Size(); i++)
      if (faces.Get(i).Valid() && i != fstind)
      {
        if(faces.Get(i).Face().NP()==4)
        {
          p1 = points.Get(faces.Get(i).Face().PNum(1)).P();
          p2 = points.Get(faces.Get(i).Face().PNum(2)).P();
          p3 = points.Get(faces.Get(i).Face().PNum(3)).P();
          p4 = points.Get(faces.Get(i).Face().PNum(4)).P();

          midp.X() = (p1.X()+p2.X()+p3.X()+p4.X())/4;
          midp.Y() = (p1.Y()+p2.Y()+p3.Y()+p4.Y())/4;
          midp.Z() = (p1.Z()+p2.Z()+p3.Z()+p4.Z())/4;

          if (Dist (midp, p0) <= xh)
          {
            locfaces.Append(faces.Get(i).Face());
            findex.Append(i);
            prism_flags.Append(faces.Get(i).PrismFlag());
          }
        }
        else
        {
          p1 = points.Get(faces.Get(i).Face().PNum(1)).P();
          p2 = points.Get(faces.Get(i).Face().PNum(2)).P();
          p3 = points.Get(faces.Get(i).Face().PNum(3)).P();

          midp.X() = (p1.X()+p2.X()+p3.X())/3;
          midp.Y() = (p1.Y()+p2.Y()+p3.Y())/3;
          midp.Z() = (p1.Z()+p2.Z()+p3.Z())/3;

          if (Dist (midp, p0) <= xh)
          {
            locfaces.Append(faces.Get(i).Face());
            findex.Append(i);
            prism_flags.Append(faces.Get(i).PrismFlag());
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
          if (pindex.Get(k) == pi)
          {
            locfaces.Elem(i).PNum(j) = k;
            found = 1;
          }

        if (!found)
        {
          pindex.Append (pi);
          locfaces.Elem(i).PNum(j) = locpoints.Append (points.Get(pi).P());
        }
      }
    }
    return(minval);
  }
}

int ADFRONT3 :: GetLocals_Pyramid(      ARRAY<Point3d> & locpoints,
                                        ARRAY<Element> & locfaces,
                                        ARRAY<INDEX> & pindex,
                                        ARRAY<INDEX> & findex,
                                        float xh)
{
  INDEX i, j, k;
  INDEX fstind;
  char found;
  INDEX pi;
  Point3d midp,p0,p1,p2,p3,p4;

  fstind = -1;
  for (i = 1; i<= faces.Size(); i++)
    if(faces.Get(i).Valid())
      if(faces.Get(i).Face().NP()==4)
      {
        fstind = i;
        break;
      }

  if(fstind==-1)
    return(-1);
  else
  {
    p1 =  points[faces[fstind].Face().PNum(1)].P();
    p2 =  points[faces[fstind].Face().PNum(2)].P();
    p3 =  points[faces[fstind].Face().PNum(3)].P();
    p4 =  points[faces[fstind].Face().PNum(4)].P();
    p0.X() = (p1.X()+p2.X()+p3.X()+p4.X())/4;
    p0.Y() = (p1.Y()+p2.Y()+p3.Y()+p4.Y())/4;
    p0.Z() = (p1.Z()+p2.Z()+p3.Z()+p4.Z())/4;

    locfaces.Append(faces[fstind].Face());
    findex.Append(fstind);

    for (i = 1; i <= faces.Size(); i++)
      if (faces.Get(i).Valid() && i != fstind)
      {
        if(faces.Get(i).Face().NP()==4)
        {
          p1 = points.Get(faces.Get(i).Face().PNum(1)).P();
          p2 = points.Get(faces.Get(i).Face().PNum(2)).P();
          p3 = points.Get(faces.Get(i).Face().PNum(3)).P();
          p4 = points.Get(faces.Get(i).Face().PNum(4)).P();

          midp.X() = (p1.X()+p2.X()+p3.X()+p4.X())/4;
          midp.Y() = (p1.Y()+p2.Y()+p3.Y()+p4.Y())/4;
          midp.Z() = (p1.Z()+p2.Z()+p3.Z()+p4.Z())/4;

          if (Dist (midp, p0) <= xh)
          {
            locfaces.Append(faces.Get(i).Face());
            findex.Append(i);
          }
        }
        else
        {
          p1 = points.Get(faces.Get(i).Face().PNum(1)).P();
          p2 = points.Get(faces.Get(i).Face().PNum(2)).P();
          p3 = points.Get(faces.Get(i).Face().PNum(3)).P();

          midp.X() = (p1.X()+p2.X()+p3.X())/3;
          midp.Y() = (p1.Y()+p2.Y()+p3.Y())/3;
          midp.Z() = (p1.Z()+p2.Z()+p3.Z())/3;

          if (Dist (midp, p0) <= xh)
          {
            locfaces.Append(faces.Get(i).Face());
            findex.Append(i);
          }
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
        if (pindex.Get(k) == pi)
        {
          locfaces.Elem(i).PNum(j) = k;
          found = 1;
        }

      if (!found)
      {
        pindex.Append (pi);
        locfaces.Elem(i).PNum(j) = locpoints.Append (points.Get(pi).P());
      }
    }
  }

  return(1);
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

void ADFRONT3 :: SetClass (INDEX fi, int i)
{
  faces[fi].SetQualClass (i);
}

int ADFRONT3 :: Prism () const
{
  int i;
  for (i = 1; i <= faces.Size(); i++)
    if (faces.Get(i).PrismFlag()!=-1)
      return(1);
  return(0);
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
  INDEX i,j,n;

  cout << "Front:" << endl;
  n = 0;
  for (i=1; i<=points.Size(); i++)
    if (points[i].Valid())
      n++;

  cout << n << " Points: " << endl;
  for (i=1; i<=points.Size(); i++)
    if (points[i].Valid())
      cout << i <<  " (" << points[i].P().X() << "," << points[i].P().Y() << "," << points[i].P().Z() << ")" << endl;
  cout << endl;

  n = 0;
  for (i=1; i<=faces.Size(); i++)
    if (faces[i].Valid())
      n++;

  cout << n << " Faces: " << endl;
  for (i=1; i<=faces.Size(); i++)
  {
    if (faces[i].Valid())
    {
      for(j=1; j<faces[i].Face().NP(); j++)
        cout << faces[i].Face().PNum(j) << " - ";
      cout << faces[i].Face().PNum(faces[i].Face().NP()) << "      ";
      cout << faces[i].QualClass() << "  " << faces[i].Valid() << endl;
    }
  }
  cout << endl;
}

void ADFRONT3 :: Grape () const
{
  INDEX i,n;
  FILE *file;

  file = fopen("grape", "w+");
  n = 0;
  for (i=1; i<=points.Size(); i++)
    n++;

  fprintf(file, "%d\n", n);

  for (i=1; i<=points.Size(); i++)
    fprintf(file, "%lf %lf %lf\n", points[i].P().X(), points[i].P().Y(), points[i].P().Z());

  n = 0;
  for (i=1; i<=faces.Size(); i++)
    if (faces[i].Valid())
      n++;

  fprintf(file, "%d\n", n);

  for (i=1; i<=faces.Size(); i++)
    if (faces[i].Valid())
      fprintf(file, "%d %d %d\n", faces[i].Face().PNum(1), faces[i].Face().PNum(2), faces[i].Face().PNum(3));
  fclose(file);
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
