// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
   Advancing front class for surfaces
 */

#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include <limits.h>

#include <template.hh>
#include <array.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

extern int Test_Line(Point3d p1, Point3d p2, Point3d sp1, Point3d sp2, double xh);



#include <meshing/adfront2.hh>


ADFRONT2::frontpoint2 :: frontpoint2 ()
{
  globalindex = 0;
  nlinetopoint = 0;
  frontnr = INT_MAX-10;    // attention: overflow on calculating  INT_MAX + 1
}

ADFRONT2::frontpoint2 :: frontpoint2 (const Point3d & ap, INDEX agi)
{
  p = ap;
  globalindex = agi;
  nlinetopoint = 0;
  frontnr = INT_MAX-10;
}




ADFRONT2::frontline :: frontline ()
{
  lineclass = 1;
  surfaceindex = 0;
}

ADFRONT2::frontline :: frontline (const ILINE & al, int asi)
{
  l = al;
  surfaceindex = asi;
  lineclass = 1;
}






ADFRONT2 :: ADFRONT2 ()
{
  nfl = 0;
}



void ADFRONT2 :: GetPoints (ARRAY<Point3d> & apoints) const
{
  INDEX i;

  for (i = 1; i <= points.Size(); i++)
    apoints.Append (points.Get(i).P());
}





INDEX ADFRONT2 :: AddPoint (const Point3d & p, INDEX globind)
{
  // inserts at empty position or resizes array
  INDEX pi;

  if (delpointl.Size() != 0)
  {
    pi = delpointl.Last();
    delpointl.DeleteLast ();

    points.Elem(pi) = frontpoint2 (p, globind);
  }
  else
  {
    pi = points.Append (frontpoint2 (p, globind));
  }

  return pi;
}


INDEX ADFRONT2 :: AddLine (INDEX pi1, INDEX pi2, int asurfaceind)
{
  int minfn;
  INDEX li;

  frontpoint2 & p1 = points[pi1];
  frontpoint2 & p2 = points[pi2];

  nfl++;

  p1.AddLine();
  p2.AddLine();

  minfn = min (p1.FrontNr(), p2.FrontNr());
  p1.DecFrontNr (minfn+1);
  p2.DecFrontNr (minfn+1);

  if (dellinel.Size() != 0)
  {
    li = dellinel.Last();
    dellinel.DeleteLast ();

    lines.Elem(li) = frontline (ILINE(pi1, pi2), asurfaceind);
  }
  else
  {
    li = lines.Append(frontline (ILINE(pi1, pi2), asurfaceind));
  }

  return li;
}


void ADFRONT2 :: DeleteLine (INDEX li)
{
  int i;
  INDEX pi;

  nfl--;

  for (i = 1; i <= 2; i++)
  {
    pi = lines.Get(li).L().I(i);
    points.Elem(pi).RemoveLine();

    if (!points.Get(pi).Valid())
      delpointl.Append (pi);
  }

  lines.Elem(li).Invalidate();
  dellinel.Append (li);
}



void ADFRONT2 :: IncrementClass (INDEX li)
{
  lines[li].IncrementClass();
}


void ADFRONT2 :: ResetClass (INDEX li)
{
  lines[li].ResetClass();
}







int ADFRONT2 :: GetLocals (ARRAY<Point3d> & locpoints,
                           ARRAY<ILINE> & loclines,       // local index
                           ARRAY<INDEX> & pindex,
                           ARRAY<INDEX> & lindex,
                           int & asurfaceind,
                           double xh)
{
  INDEX i, j, k;
  INDEX lstind, pstind;
  char found;
  INDEX pi;
  Point3d midp, p0,p1,p2,sp1,sp2,p;
  int minval, hi, flag;
  double min,minn,l;

  minval = INT_MAX;
  min = INT_MAX;

  // Auswahl der Frontkante nach Schoeberl

  /*  for (i = 1; i<= lines.Size(); i++)
      if (lines.Get(i).Valid())
        {
        hi = lines.Get(i).LineClass() +
             2 * min (points.Get(lines.Get(i).L().I1()).FrontNr(),
                      points.Get(lines.Get(i).L().I2()).FrontNr());
     cout << lines.Get(i).LineClass() << "  "
       << points.Get(lines.Get(i).L().I1()).FrontNr() << "  "
       << points.Get(lines.Get(i).L().I2()).FrontNr() << endl;
        if (hi < minval)
          {
          minval = hi;
          lstind = i;
          }
        }
   */

  // Neues Auswahlkriterium

  // Neues Auswahlkriterium

  j = 1;
  lstind = -1;
  do
  {
    for (i = 1; i<= lines.Size(); i++)
    {
      if( (lines.Get(i).LineClass()==j) && (lines.Get(i).Valid()) )
      {
        l = Dist(points.Get(lines.Get(i).L().I1()).P(),points.Get(lines.Get(i).L().I2()).P());
        /*				cout << i << "  "
                                                 << j << "  "
                                                 << lines.Get(i).LineClass() << endl;*/
        if (l < min /*|| i == 1*/)
        {
          min = l;
          lstind = i;
        }
      }
    }
    j++;
  }
  while(lstind<0);


  asurfaceind = lines[lstind].SurfaceIndex();

  pstind = lines[lstind].L().I1();
  p0 = points[pstind].P();

  loclines.Append(lines[lstind].L());
  lindex.Append(lstind);

  /*  for (i = 1; i <= lines.Size(); i++)
      {
      if (lines.Get(i).Valid() && i != lstind &&
          lines.Get(i).SurfaceIndex() == asurfaceind)
        {
        midp = Center (points.Get(lines.Get(i).L().I1()).P(),
                       points.Get(lines.Get(i).L().I2()).P());

        if (Dist (midp, p0) <= xh)
          {
          loclines.Append(lines.Get(i).L());
          lindex.Append(i);
          }
        }
      }
   */
  sp1 = points.Get(lines.Get(lstind).L().I1()).P();
  sp2 = points.Get(lines.Get(lstind).L().I2()).P();

  for (i = 1; i <= lines.Size(); i++)
  {
    if (lines.Get(i).Valid() && i != lstind)
    {
      midp = Center (points.Get(lines.Get(i).L().I1()).P(),
                     points.Get(lines.Get(i).L().I2()).P());

      p1 = points.Get(lines.Get(i).L().I1()).P();
      p2 = points.Get(lines.Get(i).L().I2()).P();

      if (Test_Line(p1,p2,sp1,sp2,xh))
      {
        loclines.Append(lines.Get(i).L());
        lindex.Append(i);
      }
    }
  }

  for (i = 1; i <= loclines.Size(); i++)
  {
    for (j = 1; j <= 2; j++)
    {
      found = 0;
      pi = loclines.Get(i).I(j);

      for (k = 1; k <= pindex.Size() && !found; k++)
      {
        if (pindex.Get(k) == pi)
        {
          loclines.Elem(i).I(j) = k;
          found = 1;
        }
      }

      if (!found)
      {
        pindex.Append (pi);
        locpoints.Append (points.Get(pi).P());
        loclines.Elem(i).I(j) = locpoints.Size();
      }
    }
  }

  return lines.Get(lstind).LineClass();
}

INDEX ADFRONT2 :: GetGlobalIndex (INDEX pi) const
{
  return points[pi].GlobalIndex();
}


void ADFRONT2 :: SetStartFront ()
{
  INDEX i;
  int j;

  for (i = 1; i <= lines.Size(); i++)
    if (lines.Get(i).Valid())
      for (j = 1; j <= 2; j++)
        points[lines[i].L().I(j)].DecFrontNr(0);

  //  for (i = 1; i <= points.Size(); i++)
  //    points[i].DecFrontNr(0);
}




void ADFRONT2 :: Print (ostream & ost) const
{
  INDEX i;

  ost << points.Size() << " Points: " << endl;
  for (i = 1; i <= points.Size(); i++)
    if (points[i].Valid())
      ost << i << "  " << points[i].P() << endl;

  ost << nfl << " Lines: " << endl;
  for (i = 1; i <= lines.Size(); i++)
    if (lines[i].Valid())
      ost << lines[i].L().I1() << " - " << lines[i].L().I2() << endl;

  ost << flush;
}

void ADFRONT2 :: ugPrint (ostream & ost) const
{
  INDEX i;

  ost << points.Size() << endl;
  for (i = 1; i <= points.Size(); i++)
    if (points[i].Valid())
      ost << " " << points[i].P().X()
          << " " << points[i].P().Y()
          << " " << points[i].P().Z() << endl;

  ost << nfl << endl;
  for (i = 1; i <= lines.Size(); i++)
    if (lines[i].Valid())
      ost << lines[i].L().I1() << " " << lines[i].L().I2() << endl;

  ost << flush;
}
