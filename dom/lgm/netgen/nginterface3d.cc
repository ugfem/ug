// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   nginterface.cc                                               */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <strstream.h>
#include <math.h>
#include <new.h>
#include <ctype.h>
#include <string.h>
#include <time.h>



#include <myadt.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>
#include <linalg/linalg.hh>

// #include <meshing/meshtool.hh>
// #include <meshing/meshing2.hh>
// #include <meshing/meshsurf.hh>

#include <meshing/global.hh>
#include <meshing/meshing3.hh>
#include <meshing/meshing2.hh>
#include "lgm.hh"

extern "C"
{
#include "../../gginterface.h"
}

//int testmode;

extern "C"
void UserWriteF (char * ch, ...);

/*void MyError (char * ch)
   {
   UserWriteF (ch);
   }*/

class InputElement : public Element
{
  int neighbour[4];
public:
  int & Neighbour (int i) { return neighbour[i-1]; }

};

ARRAY<Surface*> surfaces;
static ARRAY<geompoint3d> geompoints;
static ARRAY<splinesegment3d*> splines;
static ARRAY<Point3d> points2;
static ARRAY<int> lp1, lp2;
//static ARRAY<Point3d> points;
static ARRAY<Element> elements;
static ARRAY<InputElement> geomelements;
static int nbp;
static int tt;
static int trid;
static Point3d tripoint;
extern int yyparse ();

static ARRAY<Point3d> points;
static ARRAY<Element> * volelements;
static int disp;
static double vol0;

static ARRAY<Point3d> locpoints;
static ARRAY<Point2d> plainpoints;
static ARRAY<ILINE> loclines;
static const char * rname;
static int cntelem;
static int oldnl;
static int qualclass, surfind;
static int LGM_DEBUG = 0;

class surfacemeshing
{
  ADFRONT2 * adfront;
  ARRAY<netrule*> rules;
  ARRAY<int> ruleused;
  double cxx, cyy, czz, cxy, cxz, cyz, cx, cy, cz, c1;
  Vec3d ex, ey, ez;
  Point3d globp1;

public:
  surfacemeshing (char * rulefilename);
  virtual ~surfacemeshing ();

  void LoadRules (char * filename);
  void Mesh (double gh);

  void ImproveMesh (ARRAY<Point3d> & points, const ARRAY<Element> & elements,
                    int improveedges, int numboundarypoints, double h, int steps, int err2);

  void AddPoint (const Point3d & p, INDEX globind);
  void AddBoundaryElement (INDEX i1, INDEX i2, int surfind);
  virtual void TestPoint (const Point3d & /* p */,int flag) { };

  virtual int SavePoint (const Point3d & p);
  virtual void SaveElement (const Element & elem);

  friend int StartNetgen (double h, int smooth, int display);

protected:
  virtual void StartMesh ();
  virtual void EndMesh ();
  virtual void DefineTransformation (INDEX surfind, Point3d & p1, Point3d & p2);
  virtual void TransformToPlain (INDEX ind, const Point3d & locpoint,
                                 Point2d & plainpoint, double h);
  virtual void TransformFromPlain (INDEX surfind, Point2d & plainpoint,
                                   Point3d & locpoint, double h);
  virtual void ProjectPoint (INDEX surfind, Point3d & p) const;
  virtual void ProjectPointold (INDEX surfind, Point3d & p) const;
  virtual void ProjectPoint2 (INDEX surfind, INDEX surfind2, Point3d & p) const;
  virtual void GetNormalVector(INDEX surfind, const Point3d & p, Vec3d & n) const;
  virtual void GetNormalVectorold(INDEX surfind, const Point3d & p, Vec3d & n) const;

  virtual double CalcLocalH (const Point3d & p, int surfind, double gh) const;
};

surfacemeshing :: surfacemeshing (char * rulefilename)
{
  LoadRules (rulefilename);
  adfront = new ADFRONT2();
}


surfacemeshing :: ~surfacemeshing ()
{}

void surfacemeshing :: AddPoint (const Point3d & p, INDEX globind)
{
  adfront ->AddPoint (p, globind);
}

void surfacemeshing :: AddBoundaryElement (INDEX i1, INDEX i2, int surfind)
{
  adfront ->AddLine (i1, i2, surfind);
}
void surfacemeshing :: LoadRules (char * filename)
{
  char buf[256];
  ifstream ist (filename);

  while (!ist.eof())
  {
    buf[0] = 0;
    ist >> buf;

    if (strcmp (buf, "rule") == 0)
    {
      netrule * rule = new netrule;
      rule -> LoadRule(ist);
      rules.Append (rule);
    }
  }
}

int GetTriangleId(const Point3d & pp, Point3d & pnew)
{
  Vec3d e1,e2,e3,n1,n2,n3,np;
  Point3d ep1,ep2,ep3,p1,pold,p;
  double min,test,dist,eps;
  int mi,i;
  p = pp;

  do
  {
    pold = p;
    min = 100000000.0;
    mi = 0;
    for(i=1; i<=geomelements.Size(); i++)
    {
      ep1 = geompoints[geomelements[i].PNum(1)].p;
      ep2 = geompoints[geomelements[i].PNum(2)].p;
      ep3 = geompoints[geomelements[i].PNum(3)].p;

      e1 = p - ep1;
      e2 = p - ep2;
      e3 = p - ep3;
      test = e1.Length() + e2.Length() + e3.Length();
      if(min > test)
      {
        min = test;
        mi = i;
      }
    }
    //    cout << mi << endl;

    if(mi==0)
      cout << "E R R O R" << endl;
    ep1 = geompoints[geomelements[mi].PNum(1)].p;
    ep2 = geompoints[geomelements[mi].PNum(2)].p;
    ep3 = geompoints[geomelements[mi].PNum(3)].p;

    eps = 0.0002;

    n1 = ep3 - ep1;
    if(n1.Length()>eps)
      n1 /= n1.Length();
    else
      cout << "ERROR n1" << endl;
    n2 = ep3 - ep2;
    if(n1.Length()>eps)
      n2 /= n2.Length();
    else
      cout << "ERROR n2" << endl;
    n3 = Cross(n1,n2);
    if(n1.Length()>eps)
      n3 /= n3.Length();
    else
      cout << "ERROR n3" << endl;

    p1.X() = 0.0;
    p1.Y() = 0.0;
    p1.Z() = 0.0;
    np = p - ep1;
    np = np - (n3 * np) * n3;
    p.X() = np.X() + ep1.X();
    p.Y() = np.Y() + ep1.Y();
    p.Z() = np.Z() + ep1.Z();
    dist = sqrt( (p.X() - pold.X()) * (p.X() - pold.X())
                 + (p.Y() - pold.Y()) * (p.Y() - pold.Y())
                 + (p.Z() - pold.Z()) * (p.Z() - pold.Z()));
  }
  while( dist >= 0.05 );
  pnew = p;
  p = pold;
  return(mi);
}


void surfacemeshing :: DefineTransformation (INDEX surfind, Point3d & p1, Point3d & p2)
{
  Vec3d n1,e1,e2,e3,ex1,ez1,ez2;
  Point3d ep11,ep12,ep13,p,pnew;
  Point3d ep21,ep22,ep23;
  Point3d ep1,ep2,ep3;
  double min,test;
  int mi,i;
  int m1,m2;

  p = Center(p1,p2);
  mi = GetTriangleId(p,pnew);

  ep1 = geompoints[geomelements[mi].PNum(1)].p;
  ep2 = geompoints[geomelements[mi].PNum(2)].p;
  ep3 = geompoints[geomelements[mi].PNum(3)].p;

  ex = p2 - p1;
  ex /= ex.Length();
  ez = Cross(ep3 - ep2, ep3 - ep1);
  ez = ez - (ez * ex) * ex;
  ez /= ez.Length();

  ey = Cross(ex,ez);
  //  ey = ey - (ey * ex) * ey - (ey * ez) * ey;
  ey /= ey.Length();

  globp1 = p1;
  trid = mi;
  tripoint = p;

  /*  cout << "input " <<  mi << endl;
     cout << p1.X() << "  " << p1.Y() << "  " << p1.Z() << endl;
     cout << p2.X() << "  " << p2.Y() << "  " << p2.Z() << endl;
     cout << ex.X() << "  " << ex.Y() << "  " << ex.Z() << endl;
     cout << ey.X() << "  " << ey.Y() << "  " << ey.Z() << endl;
     cout << ez.X() << "  " << ez.Y() << "  " << ez.Z() << endl;

     cout << ex*ex << "  " << ex*ey << "  " << ex*ez  << endl;
     cout << ey*ex << "  " << ey*ey << "  " << ey*ez  << endl;
     cout << ez*ex << "  " << ez*ey << "  " << ez*ez  << endl;*/
}

void surfacemeshing :: TransformToPlain (INDEX surfind, const Point3d & locpoint,
                                         Point2d & plainpoint, double h)
{
  Vec3d p1p,pp,dirvec,vec;
  Point3d lp,np1,np2,np3;
  int tidlp,tidgp,di;
  float d1,d2,d3;
  int MAXINT2;
  MAXINT2 = 100000;
  lp = locpoint;
  tidlp = GetTriangleId(lp,lp);
  //  if(tidlp==trid)
  {
    // Standard-Transformation (wie bisher)
    p1p = locpoint - globp1;
    pp =  (p1p * ez) * ez;
    p1p = p1p - pp;
    p1p /= h;
    plainpoint.X() = p1p * ex;
    plainpoint.Y() = p1p * ey;
  }
  /*  else
     {
      if(geomelements[trid].Neighbour(1)!=-1)
      {
        np1 = Center(geompoints[geomelements[geomelements[trid].Neighbour(1)].PNum(1)].p,
                     geompoints[geomelements[geomelements[trid].Neighbour(1)].PNum(2)].p,
                     geompoints[geomelements[geomelements[trid].Neighbour(1)].PNum(3)].p);
        d1 = (np1 - tripoint).Length();
      }
      else
        d1 = MAXINT2;
      if(geomelements[trid].Neighbour(2)!=-1)
      {
        np2 = Center(geompoints[geomelements[geomelements[trid].Neighbour(2)].PNum(1)].p,
                     geompoints[geomelements[geomelements[trid].Neighbour(2)].PNum(2)].p,
                     geompoints[geomelements[geomelements[trid].Neighbour(2)].PNum(3)].p);
        d2 = (np2 - tripoint).Length();
      }
      else
        d2 = MAXINT2;
      if(geomelements[trid].Neighbour(3)!=-1)
      {
        np3 = Center(geompoints[geomelements[geomelements[trid].Neighbour(3)].PNum(1)].p,
                     geompoints[geomelements[geomelements[trid].Neighbour(3)].PNum(2)].p,
                     geompoints[geomelements[geomelements[trid].Neighbour(3)].PNum(3)].p);
        d3 = (np3 - tripoint).Length();
      }
      else
        d3 = MAXINT2;
      p1p = locpoint - globp1;
      pp =  (p1p * ez) * ez;
      p1p = p1p - pp;
      p1p /= h;
      plainpoint.X() = 1.5*p1p * ex;
      plainpoint.Y() = 1.5*p1p * ey;

     }
   */
  /*  Vec3d p1p,pp;
     p1p = locpoint - globp1;
     pp =  (p1p * ez) * ez;
     p1p = p1p - pp;
     p1p /= h;
     plainpoint.X() = p1p * ex;
     plainpoint.Y() = p1p * ey;*/
  //  cout << "geo2plain" << endl;
  //  cout << plainpoint.X() << "  " << plainpoint.Y() << endl;
  //  cout << locpoint.X() << "  " << locpoint.Y() << "  " << locpoint.Z() << endl;

}

void ProjectPoint2geo (Point3d & p)
{
  Vec3d n1,n2,n3,e1,e2,e3,np;
  Point3d ep1,ep2,ep3,pnew,p1;
  double min,test;
  int mi,i;
  mi = GetTriangleId(p,pnew);
  p = pnew;

  ep1 = geompoints[geomelements[mi].PNum(1)].p;
  ep2 = geompoints[geomelements[mi].PNum(2)].p;
  ep3 = geompoints[geomelements[mi].PNum(3)].p;

  n1 = ep3 - ep1;
  n1 /= n1.Length();
  n2 = ep3 - ep2;
  n2 /= n2.Length();
  n3 = Cross(n1,n2);
  n3 /= n3.Length();
  p1.X() = 0.0;
  p1.Y() = 0.0;
  p1.Z() = 0.0;
  np = p - ep1;
  np = np - (n3 * np) * n3;
  p.X() = np.X() + ep1.X();
  p.Y() = np.Y() + ep1.Y();
  p.Z() = np.Z() + ep1.Z();
}

void ProjectPoint2Triangle (Point3d & p, int mi)
{
  Vec3d n1,n2,n3,e1,e2,e3,np;
  Point3d ep1,ep2,ep3,pnew,p1;
  double min,test;
  int i;

  ep1 = geompoints[geomelements[mi].PNum(1)].p;
  ep2 = geompoints[geomelements[mi].PNum(2)].p;
  ep3 = geompoints[geomelements[mi].PNum(3)].p;

  n1 = ep3 - ep1;
  n1 /= n1.Length();
  n2 = ep3 - ep2;
  n2 /= n2.Length();
  n3 = Cross(n1,n2);
  n3 /= n3.Length();
  p1.X() = 0.0;
  p1.Y() = 0.0;
  p1.Z() = 0.0;
  np = p - ep1;
  np = np - (n3 * np) * n3;
  p.X() = np.X() + ep1.X();
  p.Y() = np.Y() + ep1.Y();
  p.Z() = np.Z() + ep1.Z();
}

void surfacemeshing :: TransformFromPlain (INDEX surfind, Point2d & plainpoint,
                                           Point3d & locpoint, double h)
{

  Vec3d p1p;

  p1p = plainpoint.X() * ex + plainpoint.Y() * ey;
  p1p *= h;
  locpoint = globp1 + p1p;
  /*  cout << "plain2geo" << endl;
     cout << plainpoint.X() << "  " << plainpoint.Y() << endl;
     cout << p1p.X() << "  " << p1p.Y() << "  " << p1p.Z() << endl;
     cout << globp1.X() << "  " << globp1.Y() << "  " << globp1.Z() << endl;
     cout << locpoint.X() << "  " << locpoint.Y() << "  " << locpoint.Z() << endl;*/
  // fuer den Fall, dan p1 und p2 nicht auf einem Geometriedreieck liegen
  //  ProjectPoint2geo(locpoint);
  //  ProjectPoint2Triangle(locpoint,mi);
  //  cout << locpoint.X() << "  " << locpoint.Y() << "  " << locpoint.Z() << endl;
  // cout << endl;
}


void surfacemeshing :: ProjectPoint (INDEX surfind, Point3d & p) const
{
  Vec3d e1,e2,e3,n1,n2,n3,np;
  Point3d ep1,ep2,ep3,p1,pold,pnew;
  double min,test,dist;
  int mi,i;

  mi = GetTriangleId(p,pnew);
  p = pnew;

  /*do
     {
     pold = p;
     min = 100000000.0;
     mi = 0;
     for(i=1;i<=geomelements.Size();i++)
     {
      ep1 = geompoints[geomelements[i].PNum(1)].p;
      ep2 = geompoints[geomelements[i].PNum(2)].p;
      ep3 = geompoints[geomelements[i].PNum(3)].p;

      e1 = p - ep1;
      e2 = p - ep2;
      e3 = p - ep3;
      test = e1.Length() + e2.Length() + e3.Length();
      if(min > test)
      {
        min = test;
        mi = i;
      }
     }
   */
  /*  ep1 = geompoints[geomelements[mi].PNum(1)].p;
     ep2 = geompoints[geomelements[mi].PNum(2)].p;
     ep3 = geompoints[geomelements[mi].PNum(3)].p;

     n1 = ep3 - ep1;
     n1 /= n1.Length();
     n2 = ep3 - ep2;
     n2 /= n2.Length();
     n3 = Cross(n1,n2);
     n3 /= n3.Length();
     p1.X() = 0.0;
     p1.Y() = 0.0;
     p1.Z() = 0.0;
     np = p - ep1;
     np = np - (n3 * np) * n3;
     p.X() = np.X() + ep1.X();
     p.Y() = np.Y() + ep1.Y();
     p.Z() = np.Z() + ep1.Z();
     /*  cout << "projectpoint " << mi << endl;

     cout << p.X() << "  " << p.Y() << "  " << p.Z() << endl;
     cout << n3.X() << "  " << n3.Y() << "  " << n3.Z() << endl;*/
  dist = sqrt( (p.X() - pold.X()) * (p.X() - pold.X())
               + (p.Y() - pold.Y()) * (p.Y() - pold.Y())
               + (p.Z() - pold.Z()) * (p.Z() - pold.Z()));
  //  cout << " dist " << dist << endl;
  //}
  while( dist >= 0.5 ) ;
  //  cout << endl;
}

void surfacemeshing :: ProjectPointold (INDEX surfind, Point3d & p) const
{
  Vec3d e1,e2,e3,n1,n2,n3,np;
  Point3d ep1,ep2,ep3,p1,pnew;
  double min,test;
  int mi,i;

  /*  cout << "projectpoint " << endl;
     cout << p.X() << "  " << p.Y() << "  " << p.Z() << endl;*/

  mi = GetTriangleId(p,pnew);
  min = 100000000.0;


  /*  mi = 0;

     for(i=1;i<=geomelements.Size();i++)
     {
      ep1 = geompoints[geomelements[i].PNum(1)].p;
      ep2 = geompoints[geomelements[i].PNum(2)].p;
      ep3 = geompoints[geomelements[i].PNum(3)].p;
      e1 = p - ep1;
      e2 = p - ep2;
      e3 = p - ep3;
      test = e1.Length() + e2.Length() + e3.Length();
      if(min > test)
      {
        min = test;
        mi = i;
      }
     }*/

  ep1 = geompoints[geomelements[mi].PNum(1)].p;
  ep2 = geompoints[geomelements[mi].PNum(2)].p;
  ep3 = geompoints[geomelements[mi].PNum(3)].p;

  n1 = ep3 - ep1;
  n1 /= n1.Length();
  n2 = ep3 - ep2;
  n2 /= n2.Length();
  n3 = Cross(n1,n2);
  n3 /= n3.Length();
  p1.X() = 0.0;
  p1.Y() = 0.0;
  p1.Z() = 0.0;
  np = p - ep1;
  np = np - (n3 * np) * n3;
  p.X() = np.X() + ep1.X();
  p.Y() = np.Y() + ep1.Y();
  p.Z() = np.Z() + ep1.Z();
  /*  cout << "projectpoint " << mi << endl;

     cout << p.X() << "  " << p.Y() << "  " << p.Z() << endl;
     cout << n3.X() << "  " << n3.Y() << "  " << n3.Z() << endl;*/
}

void surfacemeshing :: ProjectPoint2 (INDEX surfind, INDEX surfind2, Point3d & p) const
{
  //  ProjectToEdge ( (*surfaces)[surfind], (*surfaces)[surfind2], p);
}

void surfacemeshing :: GetNormalVector(INDEX surfind, const Point3d & pp, Vec3d & n) const
{
  Vec3d e1,e2,e3,n1,n2,n3,np;
  Point3d ep1,ep2,ep3,pold,p;
  double min,test,dist;
  int mi,i;

  p = pp;
  do
  {
    pold = p;
    min = 100000000.0;
    mi = 0;
    for(i=1; i<=geomelements.Size(); i++)
    {
      ep1 = geompoints[geomelements[i].PNum(1)].p;
      ep2 = geompoints[geomelements[i].PNum(2)].p;
      ep3 = geompoints[geomelements[i].PNum(3)].p;

      e1 = p - ep1;
      e2 = p - ep2;
      e3 = p - ep3;
      test = e1.Length() + e2.Length() + e3.Length();
      if(min > test)
      {
        min = test;
        mi = i;
      }
    }

    ep1 = geompoints[geomelements[mi].PNum(1)].p;
    ep2 = geompoints[geomelements[mi].PNum(2)].p;
    ep3 = geompoints[geomelements[mi].PNum(3)].p;

    n1 = ep3 - ep1;
    n1 /= n1.Length();
    n2 = ep3 - ep2;
    n2 /= n2.Length();
    n3 = Cross(n1,n2);
    n3 /= n3.Length();
    np = p - ep1;
    np = np - (n3 * np) * n3;
    p.X() = np.X() + ep1.X();
    p.Y() = np.Y() + ep1.Y();
    p.Z() = np.Z() + ep1.Z();
    dist = sqrt( (p.X() - pold.X()) * (p.X() - pold.X())
                 + (p.Y() - pold.Y()) * (p.Y() - pold.Y())
                 + (p.Z() - pold.Z()) * (p.Z() - pold.Z()));
  }
  while( dist >= 0.5 );

  ep1 = geompoints[geomelements[mi].PNum(1)].p;
  ep2 = geompoints[geomelements[mi].PNum(2)].p;
  ep3 = geompoints[geomelements[mi].PNum(3)].p;

  n1 = ep3 - ep1;
  n1 /= n1.Length();
  n2 = ep3 - ep2;
  n2 /= n2.Length();
  n = Cross(n1,n2);
  n /= n.Length();
  /*  cout << "normalvector " <<  mi << endl;
     cout << n.X() << "  " << n.Y() << "  " << n.Z() << endl;
     cout << "im punkt " << p.X() << "  "
                        << p.Y() << "  "
                        << p.Z() << endl;*/
}

void surfacemeshing :: GetNormalVectorold(INDEX surfind, const Point3d & p, Vec3d & n) const
{
  Vec3d e1,e2,e3,n1,n2,n3;
  Point3d ep1,ep2,ep3;
  double min,test;
  int mi,i;

  min = 100000000.0;
  mi = 0;
  for(i=1; i<=geomelements.Size(); i++)
  {
    ep1 = geompoints[geomelements[i].PNum(1)].p;
    ep2 = geompoints[geomelements[i].PNum(2)].p;
    ep3 = geompoints[geomelements[i].PNum(3)].p;
    e1 = p - ep1;
    e2 = p - ep2;
    e3 = p - ep3;
    test = e1.Length() + e2.Length() + e3.Length();
    if(min > test)
    {
      min = test;
      mi = i;
    }
  }

  ep1 = geompoints[geomelements[mi].PNum(1)].p;
  ep2 = geompoints[geomelements[mi].PNum(2)].p;
  ep3 = geompoints[geomelements[mi].PNum(3)].p;

  n1 = ep3 - ep1;
  n1 /= n1.Length();
  n2 = ep3 - ep2;
  n2 /= n2.Length();
  n = Cross(n1,n2);
  n /= n.Length();
  /*  cout << "normalvector " <<  mi << endl;
     cout << n.X() << "  " << n.Y() << "  " << n.Z() << endl;
     cout << "im punkt " << p.X() << "  "
                        << p.Y() << "  "
                        << p.Z() << endl;*/
}

double surfacemeshing :: CalcLocalH (const Point3d & p, int surfind, double gh) const
{
  return gh;
}

int surfacemeshing :: SavePoint (const Point3d & p)
{
  return points.Append (p);
}

void surfacemeshing :: SaveElement (const Element & elem)
{
  elements.Append (elem);
}


class my_surfacemeshing : public surfacemeshing
{
public:
  //  my_meshinggeosurfaces (ARRAY<Surface*> & asurfaces);
  my_surfacemeshing (char * rulefilename);
  virtual int SavePoint (const Point3d & p);
  virtual void SaveElement (const Element & elem);
  virtual void TestPoint (const Point3d & ,int flag);
};

my_surfacemeshing :: my_surfacemeshing (char * rulefilename)
  : surfacemeshing (rulefilename)
{
  ;
}

int my_surfacemeshing :: SavePoint (const Point3d & p)
{
  return points.Append (p);
}

/*void my_surfacemeshing :: SaveElement (const Element & elem)
   {
   elements.Append (elem);
   }*/

void my_surfacemeshing :: TestPoint (const Point3d & ,int flag)
{
  int i;
  for(i=1; i<=elements.Size(); i++)
  {}
}

#define V3_ADD(A,B,C)                              {(C)[0] = (A)[0] + (B)[0];\
                                                    (C)[1] = (A)[1] + (B)[1];\
                                                    (C)[2] = (A)[2] + (B)[2];}

#define V3_SCALE(c,C)                              {(C)[0] = (c)*(C)[0];\
                                                    (C)[1] = (c)*(C)[1];\
                                                    (C)[2] = (c)*(C)[2];}

#define V3_EUKLIDNORM_OF_DIFF(A,B,b)    (b) = (sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])+((A)[2]-(B)[2])*((A)[2]-(B)[2]))));

#define V3_CLEAR(A)                                {(A)[0] = 0.0; (A)[1]= 0.0; (A)[2] = 0.0;}

void my_surfacemeshing :: SaveElement (const Element & elem)
{
  float x[4][3],diam,fac,global[3],inndiam,dist,percent;
  int i,n;
  float volume,percent1;

  elements.Append (elem);
}


static my_surfacemeshing * meshing;
static const surfacemeshing * meshthis;


//static surfacemeshing * meshing;

int AddGeomPoint (int id, double x, double y, double z)
{
  geompoints.Append (geompoint3d());
  geompoints.Last().p.X() = x;
  geompoints.Last().p.Y() = y;
  geompoints.Last().p.Z() = z;

  if (LGM_DEBUG) printf("%s %f %f %f \n","GEOMPOINT",x,y,z);

  return 0;
}

int AddGeomElement (int node0, int node1, int node2)
{
  geomelements.Append(InputElement());
  geomelements.Last().PNum(1) = node0;
  geomelements.Last().PNum(2) = node1;
  geomelements.Last().PNum(3) = node2;
  geomelements.Last().Neighbour(1) = 0;
  geomelements.Last().Neighbour(2) = 0;
  geomelements.Last().Neighbour(3) = 0;

  if (LGM_DEBUG) printf("%s %d %d %d \n","SURFACEELEMENT",node0,node1,node2);

  return 0;
}

int AddLinePoint (int id, double x, double y, double z)
{
  points.Append (Point3d(x,y,z));

  //  points->Append (Point3d());
  //  points->Last().X() = x;
  //  points->Last().Y() = y;
  //  points->Last().Z() = z;
  meshing->AddPoint (points.Last(), id);
  if (LGM_DEBUG) printf("%s %d %f %f %f \n","POINT",id,x,y,z);
  return 0;
}

int AddLineSegment (int i1,int i2)
{
  meshing->AddBoundaryElement (i1, i2, 1);
  if (LGM_DEBUG) printf("%s %d %d \n","LINESEGMENT",i1,i2);
  return 0;
}

void surfacemeshing :: StartMesh ()
{
  int i;

  ruleused.SetSize (rules.Size());
  for (i = 1; i <= ruleused.Size(); i++)
    ruleused[i] = 0;
  cntelem = 0;
}

void surfacemeshing :: EndMesh ()
{
  int i;

  for (i = 1; i <= ruleused.Size(); i++)
    (*testout) << setw(4) << ruleused[i]
               << " times used rule " << rules[i] -> Name() << endl;
}

void surfacemeshing :: Mesh (double gh)
{
  ARRAY<INDEX> pindex, lindex;
  ARRAY<int> delpoints, dellines;
  ARRAY<Element> locelements;
  int i, j, oldnp,flag,test;
  char ch, found;
  INDEX globind;
  int old;
  double in[5];
  static ARRAY<Point3d> locp;

  double h;
  /*  Point2d bemp, bemp1, bemp2;*/
  Point3d bemp, bemp1, bemp2;

  // cout << "Inputpoints for netgen" << endl;
  for (i = 1; i <= geompoints.Size(); i++)
  {
    /*cout << geompoints[i].p.X() << "  "
        << geompoints[i].p.Y() << "  "
        << geompoints[i].p.Z() << endl;*/
  }

  // cout << "Inputelements for netgen" << endl;
  for(i=1; i<=geomelements.Size(); i++)
  {
    /*cout << geomelements[i].PNum(1) << "  "
        << geomelements[i].PNum(2) << "  "
        << geomelements[i].PNum(3) << endl;*/
  }



  testmode = 0;
  meshthis = this;

  StartMesh();

  adfront ->SetStartFront ();
  ch = 0;
  //  cin >> test;
  while (ch != 'e' && !adfront ->Empty() /*&& cntelem<test*/)

  {

    h = gh;
    if (h>0.0)
      qualclass =
        adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                             surfind, 3 * h);
    else
      qualclass =
        adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                             surfind, -3 * h);
    oldnp = locpoints.Size();
    oldnl = loclines.Size();


    /*    bemp1.X() = locpoints[loclines[1].I1()].X();
        bemp1.Y() = locpoints[loclines[1].I1()].Y();
        bemp2.X() = locpoints[loclines[1].I2()].X();
        bemp2.Y() = locpoints[loclines[1].I2()].Y();

        bemp = Center (bemp1, bemp2);
        bemp.Y() += 0.1 * (bemp2.X() - bemp1.X());
        bemp.X() -= 0.1 * (bemp2.Y() - bemp1.Y());*/

    bemp1.X() = locpoints[loclines[1].I1()].X();
    bemp1.Y() = locpoints[loclines[1].I1()].Y();
    bemp1.Z() = locpoints[loclines[1].I1()].Z();
    bemp2.X() = locpoints[loclines[1].I2()].X();
    bemp2.Y() = locpoints[loclines[1].I2()].Y();
    bemp2.Z() = locpoints[loclines[1].I2()].Z();

    bemp = Center (bemp1, bemp2);

    if(gh<=0.0)
    {
      in[0] = bemp.X();
      in[1] = bemp.Y();
      in[2] = bemp.Z();
      in[3] = gh;
      Get_Local_h(in,&h);
    }

    DefineTransformation (surfind, locpoints[loclines[1].I1()],
                          locpoints[loclines[1].I2()]);

    plainpoints.SetSize (locpoints.Size());
    for (i = 1; i <= locpoints.Size(); i++)
    {
      /*      if (LGM_DEBUG) printf("%f %f %f \n",locpoints[i].X(),locpoints[i].Y(),locpoints[i].Z());*/
      TransformToPlain (surfind, locpoints[i], plainpoints[i], h);
    }
    /*    cout << " locpoints and plainpoints" << endl;
        for (i = 1; i <= locpoints.Size(); i++)
        {
          cout << locpoints[i].X() << "  " <<  locpoints[i].Y() << "  " <<  locpoints[i].Z() << endl;
          cout << plainpoints[i].X() << "  " <<  plainpoints[i].Y() << endl;
        }*/

    for (i = 2; i <= loclines.Size(); i++)  // don't remove first line
    {
      if (plainpoints[loclines[i].I1()].X() > 1e6 ||
          plainpoints[loclines[i].I2()].X() > 1e6)
      {
        loclines.DeleteElement(i);
        lindex.DeleteElement(i);
        oldnl--;
        i--;
      }
    }

    found = ApplyRules (rules, plainpoints, loclines, locelements,
                        dellines, qualclass);

    flag = 1;

    if (found)
    {
      //      adfront->Print(cout);
      ruleused[found]++;

      locpoints.SetSize (plainpoints.Size());
      for (i = oldnp+1; i <= plainpoints.Size(); i++)
        TransformFromPlain (surfind, plainpoints[i], locpoints[i], h);

      for (i = 1; i <= oldnl; i++)
        adfront -> ResetClass (lindex[i]);

      pindex.SetSize(locpoints.Size());

      for (i = oldnp+1; i <= locpoints.Size(); i++)
      {
        globind = SavePoint (locpoints[i]);
        pindex[i] = adfront -> AddPoint (locpoints[i], globind);
      }

      /*      for(i=1;i<=pindex.Size();i++)
                      cout << i << pindex[i] << endl;*/

      for (i = oldnl+1; i <= loclines.Size(); i++)
      {
        adfront -> AddLine (pindex[loclines[i].I1()],
                            pindex[loclines[i].I2()], surfind);
        /*        cout << " addline "
                     << loclines[i].I1() << "  "
                     << loclines[i].I2() << "  "
                     << pindex[loclines[i].I1()] << "  "
                     << pindex[loclines[i].I2()] << "  "
                     << plainpoints[pindex[loclines[i].I1()]].X() << "  "
                     << plainpoints[pindex[loclines[i].I1()]].Y() << "  "
                     << plainpoints[pindex[loclines[i].I2()]].X() << "  "
                     << plainpoints[pindex[loclines[i].I2()]].Y() << "  "
                     << endl;*/
      }

      /*      for(i=1;i<=pindex.Size();i++)
                      cout << i << pindex[i] << endl;*/

      for (i = 1; i <= locelements.Size(); i++)
      {
        for (j = 1; j <= locelements[i].NP(); j++)
        {
          locelements[i].PNum(j) =
            adfront -> GetGlobalIndex (pindex[locelements[i].PNum(j)]);

          //	  old = locelements[i].PNum(j);
          //          locelements[i].PNum(j) =
          //              adfront -> GetGlobalIndex (pindex[locelements[i].PNum(j)]);
          // das hier ging immer mit adfront -> GetGlobalIndex
          //          locelements[i].PNum(j) = pindex[locelements[i].PNum(j)];
          //	  cout << old << "  " << locelements[i].PNum(j)<< endl;
        }
        /*	cout << endl;
                cout << locelements[i].NP() << "  "
                     << locelements[i].PNum(1) << "  "
                     << locelements[i].PNum(2) << "  "
                     << locelements[i].PNum(3) << endl;*/

        locelements[i].SetSurfaceIndex (surfind);

        SaveElement (locelements[i]);
        cntelem++;
      }

      for (i = 1; i <= dellines.Size(); i++)
      {
        adfront -> DeleteLine (lindex[dellines[i]]);
      }

      rname = rules[found]->Name();
      adfront->GetPoints(locp);
    }
    else
    {
      adfront -> IncrementClass (lindex[1]);
    }

    locpoints.SetSize(0);
    loclines.SetSize(0);
    pindex.SetSize(0);
    lindex.SetSize(0);
    delpoints.SetSize(0);
    dellines.SetSize(0);
    locelements.SetSize(0);

    //    adfront->Print(cout);

  }

  EndMesh ();
  geompoints.SetSize(0);
  geomelements.SetSize(0);
}


int InitSurfaceNetgen (char * rulefilename)
{
  testout = new ofstream("test.out");
  meshing = new my_surfacemeshing(rulefilename);

  //  points = new ARRAY<Point3d>;
  volelements = new ARRAY<Element>;
  nbp = 0;

  return 0;
}

int StartSurfaceNetgen (double h, int smooth, int display)
{
  int i;

  nbp = points.Size();
  disp = display;
  //  vol0 = meshing->adfront->Volume();

  meshing -> Mesh (h);

  //  for (i=0; i<smooth; i++)
  //	meshing -> ImproveMesh (points, elements, 0, nbp, h,20,1);
  //points, elements, 0, nbp, h, 20, 1
  //  UserWriteF("\n");
  for (i = nbp + 1; i <= points.Size(); i++)
    AddInnerNode2ug (points.Get(i).X(),
                     points.Get(i).Y(),
                     points.Get(i).Z());

  for (i = 1; i <= elements.Size(); i++)
  {
    AddSurfaceTriangle2ug (elements.Get(i).PNum(1) - 1,
                           elements.Get(i).PNum(2) - 1,
                           elements.Get(i).PNum(3) - 1);
  }

  elements.SetSize(0);
  points.SetSize(0);
  return 0;
}
