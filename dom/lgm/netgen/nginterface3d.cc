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
static double SMALL = 0.0005;
static double Triangle_Angle2 = 60.0;
double hh;

int GetEdgeId(const Point3d & ep1, Point3d & ep2);
int Calc_Coord_Vectors(const Point3d p1,const Point3d p2,const int mi,Vec3d & nx,Vec3d & ny,Vec3d & nz);

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
  virtual int DefineTransformation (INDEX surfind, Point3d & p1, Point3d & p2);
  virtual int DefineTransformation_OLD (INDEX surfind, Point3d & p1, Point3d & p2);
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


int surfacemeshing :: DefineTransformation (INDEX surfind, Point3d & p1, Point3d & p2)
{
  int mi,i;

  mi = GetEdgeId(p1,p2);

  Calc_Coord_Vectors(p1,p2,mi,ex,ey,ez);

  //  globp1 = Center(p1,p2);
  globp1 = p1;

  return(mi);
}

int surfacemeshing :: DefineTransformation_OLD (INDEX surfind, Point3d & p1, Point3d & p2)
{
  /*

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
     tripoint = p;*/

  return (0);
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
     cout << "projectpoint " << mi << endl;

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

int AddGeomElement (int node0, int node1, int node2, int neigbor0, int neigbor1, int neigbor2)
{
  geomelements.Append(InputElement());
  geomelements.Last().PNum(1) = node0;
  geomelements.Last().PNum(2) = node1;
  geomelements.Last().PNum(3) = node2;
  geomelements.Last().Neighbour(1) = neigbor0;
  geomelements.Last().Neighbour(2) = neigbor1;
  geomelements.Last().Neighbour(3) = neigbor2;

  if (LGM_DEBUG)
    printf("%s %d %d %d %d %d %d \n","SURFACEELEMENT",node0,node1,node2,
           neigbor0,neigbor1,neigbor2);

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

  if(LGM_DEBUG)
  {
    // Ausgabe fuer Netgen
    printf("%s\n","ug");
    printf("%d\n",geompoints.Size());
    for (i = 1; i <= geompoints.Size(); i++)
    {
      cout << geompoints[i].p.X() << "  "
           << geompoints[i].p.Y() << "  "
           << geompoints[i].p.Z() << endl;
    }
    printf("%d\n",geomelements.Size());
    for(i=1; i<=geomelements.Size(); i++)
    {
      cout << geomelements[i].PNum(1) << "  "
           << geomelements[i].PNum(2) << "  "
           << geomelements[i].PNum(3) << "  "
           << geomelements[i].Neighbour(1) << "  "
           << geomelements[i].Neighbour(2) << "  "
           << geomelements[i].Neighbour(3)
           << endl;
    }

    // Ausgabe der Front
    adfront->ugPrint(cout);
  }

  testmode = 0;
  meshthis = this;

  StartMesh();

  adfront ->SetStartFront ();
  ch = 0;
  //  cin >> test;
  while (ch != 'e' && !adfront ->Empty() /*&& cntelem<test*/)

  {
    // adfront->Print(cout);
    locpoints.SetSize(0);
    loclines.SetSize(0);
    pindex.SetSize(0);
    lindex.SetSize(0);

    h = gh;
    if (h>0.0)
      qualclass =
        adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                             surfind, 3 * h);
    else
      qualclass =
        adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                             surfind, -3 * h);
    //	cout << "qualclass:   " << qualclass << endl;

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

    bemp.Y() += 0.1 * (bemp2.X() - bemp1.X());
    bemp.X() -= 0.1 * (bemp2.Y() - bemp1.Y());

    locpoints.SetSize(0);
    loclines.SetSize(0);
    pindex.SetSize(0);
    lindex.SetSize(0);
    if(gh<=0.0)
    {
      in[0] = bemp.X();
      in[1] = bemp.Y();
      in[2] = bemp.Z();
      in[3] = gh;
      Get_Local_h(in,&h);
    }
    qualclass =
      adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                           surfind, 3 * h);

    oldnp = locpoints.Size();
    oldnl = loclines.Size();

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

        locelements[i].SetSurfaceIndex (surfind);

        SaveElement (locelements[i]);
        cntelem++;
        if(LGM_DEBUG)
          cout << cntelem << "  "
               << locelements[i].PNum(1) << "  "
               << locelements[i].PNum(2) << "  "
               << locelements[i].PNum(3) << endl;
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
}

#define Det(a)                  ( - a[0] * a[4] * a[8]  \
                                  + a[0] * a[5] * a[7]  \
                                  + a[3] * a[1] * a[8]  \
                                  - a[3] * a[2] * a[7]  \
                                  - a[6] * a[1] * a[5]  \
                                  + a[6] * a[2] * a[4])

#define InvMatMult(b,c,a)               b[0] =    ( a[5] * a[7] - a[4] * a[8] ) * c[0]          \
                                               + ( a[1] * a[8] - a[2] * a[7] ) * c[1]          \
                                               + ( a[2] * a[4] - a[1] * a[5] ) * c[2];         \
  b[1] =    ( a[3] * a[8] - a[5] * a[6] ) * c[0]          \
         + ( a[2] * a[6] - a[0] * a[8] ) * c[1]          \
         + ( a[0] * a[5] - a[2] * a[3] ) * c[2];         \
  b[2] =    ( a[4] * a[6] - a[3] * a[7] ) * c[0]          \
         + ( a[0] * a[7] - a[1] * a[6] ) * c[1]          \
         + ( a[1] * a[3] - a[0] * a[4] ) * c[2];         \
  b[0] = b[0] / Det(a);                                   \
  b[1] = b[1] / Det(a);                                   \
  b[2] = b[2] / Det(a);

#define Det2d(a)                ( a[0] * a[3] - a[1] * a[2] )

#define InvMatMult2d(b,c,a)             b[0] =    a[3] * c[0]           \
                                               - a[1] * c[1];          \
  b[1] =  - a[2] * c[0]           \
         + a[0] * c[1];          \
  b[0] = b[0] / Det2d(a);         \
  b[1] = b[1] / Det2d(a);


int Calc_Vectors(const Point3d p0,
                 const Point3d p1,
                 const Point3d p2,
                 Vec3d & n0,
                 Vec3d & n1,
                 Vec3d & n2)
{
  // compute 2 vectors in the plane + normalvektor
  n0 = p2 - p0;
  n0 /= n0.Length();
  n1 = p2 - p1;
  n1 /= n1.Length();

  n2 = Cross(n0,n1);
  n2 /= n2.Length();
  return(0);
}


double ABS(double value)
{
  if(value<0.0)
    return(-value);
  else
    return(value);
}

double Calc_Local_Coordinates(const Point3d p0,
                              const Point3d p1,
                              const Point3d p2,
                              const Point3d p,
                              double & lam0,
                              double & lam1,
                              double & lam2)
{
  double a[9],aa[4],cc[2];
  double s1[3],s2[3],s3[3],det;
  Point3d np,pp1,pp2;
  Vec3d n0,n1,n2,nv;
  //double bb[3];
  nv = p - p0;

  Calc_Vectors(p0,p1,p2,n0,n1,n2);

  pp1.X() = p.X() - (n2*nv) * n2.X();
  pp1.Y() = p.Y() - (n2*nv) * n2.Y();
  pp1.Z() = p.Z() - (n2*nv) * n2.Z();

  // calculate local coordinates
  a[0] = p0.X();
  a[1] = p1.X();
  a[2] = p2.X();

  a[3] = p0.Y();
  a[4] = p1.Y();
  a[5] = p2.Y();

  a[6] = p0.Z();
  a[7] = p1.Z();
  a[8] = p2.Z();

  aa[0] = a[0] - a[2];
  aa[1] = a[1] - a[2];
  aa[2] = a[3] - a[5];
  aa[3] = a[4] - a[5];

  if(ABS(Det2d(aa))>0.0001)
  {
    det = ABS(Det2d(aa));
    cc[0] = pp1.X() - a[2];
    cc[1] = pp1.Y() - a[5];
    InvMatMult2d(s1,cc,aa);
    s1[2] = 1 - s1[0] - s1[1];
    lam0 = s1[0];
    lam1 = s1[1];
    lam2 = 1 - s1[0] - s1[1];
  }
  else
  {
    aa[0] = a[0] - a[2];
    aa[1] = a[1] - a[2];
    aa[2] = a[6] - a[8];
    aa[3] = a[7] - a[8];
    if(ABS(Det2d(aa))>0.001)
    {
      det = ABS(Det2d(aa));
      cc[0] = pp1.X() - a[2];
      cc[1] = pp1.Z() - a[8];
      InvMatMult2d(s2,cc,aa);
      s2[2] = 1 - s2[0] - s2[1];
      lam0 = s2[0];
      lam1 = s2[1];
      lam2 = 1 - s2[0] - s2[1];
    }
    else
    {
      aa[0] = a[3] - a[5];
      aa[1] = a[4] - a[5];
      aa[2] = a[6] - a[8];
      aa[3] = a[7] - a[8];
      if(ABS(Det2d(aa))>0.001)
      {
        det = ABS(Det2d(aa));
        cc[0] = pp1.Y() - a[5];
        cc[1] = pp1.Z() - a[8];
        InvMatMult2d(s3,cc,aa);
        s3[2] = 1 - s3[0] - s3[1];
        lam0 = s3[0];
        lam1 = s3[1];
        lam2 = 1 - s3[0] - s3[1];
      }
      else
      {
        printf("%s\n","hier geht was schief");
      }
    }
  }
  return(det);
}


int Calc_Coord_Vectors(const Point3d p1,
                       const Point3d p2,
                       const int mi,
                       Vec3d & nx,
                       Vec3d & ny,
                       Vec3d & nz)
{
  Point3d ep1,ep2,ep3;

  // Corner des Inputdreiecks
  ep1 = geompoints[geomelements[mi].PNum(1)].p;
  ep2 = geompoints[geomelements[mi].PNum(2)].p;
  ep3 = geompoints[geomelements[mi].PNum(3)].p;

  nz = Cross(ep2 - ep1, ep3 - ep1);
  nz /= nz.Length();

  nx = p2 - p1;
  nx = nx - (nx * nz) * nz;
  nx /= nx.Length();

  ny = Cross(nz,nx);
  ny /= ny.Length();

  return(0);
}

double Project2Plane(Point3d & p,
                     Vec3d & np,
                     Point3d & p0,
                     Vec3d & n0,
                     Vec3d & n1,
                     Vec3d & n2)
{
  double dist, dist1;
  np = p - p0;

  dist = (n2 * np);
  np = np - (n2 * np) * n2;

  if(dist<0.0)
    dist = -dist;
  return(dist);
}


// ***********************************************************************************
// Funktionen fuer die Projektion eines Punktes auf die Surface
// ***********************************************************************************

int c1(const Point3d & point, Point3d & returnpoint)
{
  int i,j,triang[10][4],trnb, mi, nn1, nn2, np1, np2;
  double lam[3],prod,help;
  Point3d p0,p1,p2, p;
  Vec3d edge,geomedge,n0,n1,n2,np;
  InputElement es;

  p = point;
  /* Punkt liegt auf der Surface
   * 0<=lamda_i<=1,  dist < eps
   * Falls der Punkt genau auf einer Kante zwischen 2 Dreiecken liegt,
   * ist es egal welches genommen wird
   */

  mi = 0;
  trnb = 0;
  for(i=1; i<=geomelements.Size(); i++)
  {
    p0 = geompoints[geomelements[i].PNum(1)].p;
    p1 = geompoints[geomelements[i].PNum(2)].p;
    p2 = geompoints[geomelements[i].PNum(3)].p;

    Calc_Local_Coordinates(p0,p1,p2,p,lam[0],lam[1],lam[2]);

    Calc_Vectors(p0,p1,p2,n0,n1,n2);
    help = Project2Plane(p,np,p0,n0,n1,n2);

    if( (lam[0]>=-SMALL) && (lam[1]>=-SMALL) && (lam[2]>=-SMALL) && (help<SMALL) )
    {
      mi = i;
      returnpoint.X() = np.X() + p0.X();
      returnpoint.Y() = np.Y() + p0.Y();
      returnpoint.Z() = np.Z() + p0.Z();
      break;
    }
  }

  return(mi);
}

int c2(Point3d & p, Point3d & returnpoint, double &min)
{
  int i,j, mi;
  double lam[3],help;
  Point3d p0,p1,p2;
  Vec3d n0,n1,n2,np;

  /* Punkt liegt nicht auf der Surface
   * Projeziere den Punkt in ein Dreieck mit  0<=lamda_i<=1 sodass der
   * Abstand zur Surface minimal wird
   */

  min = 10000000.0;
  mi = 0;
  for(i=1; i<=geomelements.Size(); i++)
  {
    p0 = geompoints[geomelements[i].PNum(1)].p;
    p1 = geompoints[geomelements[i].PNum(2)].p;
    p2 = geompoints[geomelements[i].PNum(3)].p;

    Calc_Local_Coordinates(p0,p1,p2,p,lam[0],lam[1],lam[2]);

    Calc_Vectors(p0,p1,p2,n0,n1,n2);
    help = Project2Plane(p,np,p0,n0,n1,n2);

    if( (lam[0]>=-SMALL) && (lam[1]>=-SMALL) && (lam[2]>=-SMALL) && (help<min) )
    {
      mi = i;
      min = help;
      returnpoint.X() = np.X() + p0.X();
      returnpoint.Y() = np.Y() + p0.Y();
      returnpoint.Z() = np.Z() + p0.Z();
    }
  }
  return(mi);
}

int c3(Point3d & p, Point3d & returnpoint, double &min)
{
  int i,j, mi;
  double lam[3],help, d0, d1, dist, m;
  Point3d p0,p1,p2, point;
  Vec3d n0,n1,n2,np, dist_vec;

  /* Punkt liegt nicht auf der Surface
   * Projeziere den Punkt auf die Kante,  sdass der
   *Abstand minimal wird */

  min = 10000000.0;

  for(i=1; i<=geomelements.Size(); i++)
  {
    for(j=0; j<3; j++)
    {
      p0 = geompoints[geomelements[i].PNum(j%3+1)].p;
      p1 = geompoints[geomelements[i].PNum((j+1)%3+1)].p;

      m = ( (p1 - p0) * (p - p0) ) / ( (p1 - p0) * (p1 - p0) );

      if((0.0<=m) && (m<=1.0))
      {
        point = p0 + m * (p1 - p0);

        dist_vec = p - point;
        help = dist_vec.Length();
      }
      else
      {
        dist_vec = p - p0;
        d0 = dist_vec.Length();

        dist_vec = p - p1;
        d1 = dist_vec.Length();

        if(d0>d1)
        {
          m = 0.0;
          help = d1;
          point = p1;
        }
        else
        {
          m = 1.0;
          help = d0;
          point = p0;
        }
      }

      if( (help<min) )
      {
        mi = i;
        min = help;
        returnpoint.X()=point.X();
        returnpoint.Y()=point.Y();
        returnpoint.Z()=point.Z();
      }
    }
  }
  return(mi);
}

int Project_Point2Surface(Point3d &inpoint, Point3d &outpoint)
{
  int i, j, mi,mi2,mi3;
  double min2, min3;
  Point3d p1, p2, p3;
  mi = 0;
  mi = c1(inpoint, p1);
  if(mi!=0)
    outpoint = p1;
  else
  {
    mi2 = c2(inpoint, p2, min2);
    mi3 = c3(inpoint, p3, min3);
    if((mi2==0)&&(mi3==0))
      printf("%s\n", "schotter");

    if(min2<min3)
    {
      mi = mi2;
      outpoint = p2;
    }
    else
    {
      mi = mi3;
      outpoint = p3;
    }
  }
  return(mi);
}

// ***********************************************************************************
//
// ***********************************************************************************

// ***********************************************************************************
//
// ***********************************************************************************



int case1(const Point3d & ep1, Point3d & ep2)
{
  int i,j,triang[10][4],trnb, mi, nn1, nn2, np1, np2;
  double lam[3],prod,help,det;
  Point3d p,p0,p1,p2,inp;
  Vec3d edge,geomedge,n0,n1,n2,np;
  InputElement es;

  // Punkt liegt genau auf der Kante zwischen 2 Dreiecken
  // passiert fuer erkannte Kanten auf der Surface

  mi = 0;
  trnb = 0;
  inp = Center(ep1,ep2);
  Project_Point2Surface(inp,p);
  for(i=1; i<=geomelements.Size(); i++)
  {

    edge = ep2 - ep1;
    edge /= edge.Length();
    p0 = geompoints[geomelements[i].PNum(1)].p;
    p1 = geompoints[geomelements[i].PNum(2)].p;
    p2 = geompoints[geomelements[i].PNum(3)].p;

    det = Calc_Local_Coordinates(p0,p1,p2,p,lam[0],lam[1],lam[2]);

    Calc_Vectors(p0,p1,p2,n0,n1,n2);
    help = Project2Plane(p,np,p0,n0,n1,n2);

    //cout << i << "  " << lam[0] << "  " << lam[1] << "  " << lam[2] << "  " << help << endl;

    if( (lam[0]>=-SMALL) && (lam[1]>=-SMALL) && (lam[2]>=-SMALL) && (help<SMALL) )
    {
      triang[trnb][0] = i;
      triang[trnb][1] = -1;
      triang[trnb][2] = -1;
      triang[trnb][3] = -1;
      Calc_Vectors(p0,p1,p2,n0,n1,n2);
      help = Project2Plane(p,np,p0,n0,n1,n2);
      for(j=0; j<3; j++)
        if(lam[j]<SMALL/det)
          triang[trnb][j+1] = j;
      trnb++;
    }
  }

  if(trnb==1)
  {
    mi = triang[0][0];
    return(mi);
  }
  else
  {
    for(i=0; i<trnb; i++)
    {
      if( (triang[i][1]==-1) && (triang[i][2]==-1) && (triang[i][3]==-1) )
        return(triang[i][0]);

      for(j=1; j<=3; j++)
      {
        if(triang[i][j]!=-1)
        {
          nn1 = (triang[i][j]+1)%3+1;
          nn2 = (triang[i][j]+2)%3+1;
          es = geomelements[triang[i][0]];
          np1 = es.PNum(nn1);
          np2 = es.PNum(nn2);
          geomedge = geompoints[np1].p - geompoints[np2].p;
          geomedge /= geomedge.Length();
          prod = edge*geomedge;
          if(prod<=0)
            mi = triang[i][0];
        }
      }
    }
  }

  return(mi);
}

int case2(const Point3d & ep1, Point3d & ep2)
{
  int i,j, mi;
  double lam[3],help,min;
  Point3d p,p0,p1,p2,inp;
  Vec3d edge,geomedge,n0,n1,n2,np;

  // Punkt liegt nicht auf der Surface, eps<lam[i]<1-eps
  // fuer mehrere Dreiecke
  // Waehle Dreieck mit kleinstem Abstand

  mi = 0;
  min = 1000000.0;
  inp = Center(ep1,ep2);
  Project_Point2Surface(inp,p);
  for(i=1; i<=geomelements.Size(); i++)
  {
    p0 = geompoints[geomelements[i].PNum(1)].p;
    p1 = geompoints[geomelements[i].PNum(2)].p;
    p2 = geompoints[geomelements[i].PNum(3)].p;

    Calc_Local_Coordinates(p0,p1,p2,p,lam[0],lam[1],lam[2]);

    Calc_Vectors(p0,p1,p2,n0,n1,n2);
    help = Project2Plane(p,np,p0,n0,n1,n2);

    if(help<min)
      if( (lam[0]>=0.0) && (lam[0]<=1.0) )
        if( (lam[1]>=0.0) && (lam[1]<=1.0) )
          if( (lam[2]>=0.0) && (lam[2]<=1.0) )
          {
            min = help;
            mi = i;
          }
  }
  if(min > 0.2*hh)
    return(0);
  else
    return(mi);
}

int case3(const Point3d & ep1, Point3d & ep2)
{
  int i,j, mi;
  double lam[3],help,min;
  Point3d p,p0,p1,p2,inp;
  Vec3d edge,geomedge,n0,n1,n2,np;

  // Punkt auf der Startfront. Aufgrund der neuen LineDisc
  // liegt die aktuelle Front ausserhalb der Dreiecke
  // Finde Dreieck, mit geringstem Abstand
  // Kante-Punkt

  mi = 0;
  min = 1000000.0;
  inp = Center(ep1,ep2);
  Project_Point2Surface(inp,p);
  for(i=1; i<=geomelements.Size(); i++)
  {
    for(j=0; j<3; j++)
    {
      p0 = geompoints[geomelements[i].PNum((j)%3+1)].p;
      p1 = geompoints[geomelements[i].PNum((j+1)%3+1)].p;

      help = Dist(Center(p0, p1), p);
      if(help<min)
      {
        min = help;
        mi =i;
      }
    }
  }
  if(min > 0.2*hh)
    return(0);
  else
    return(mi);
}

int case4(const Point3d & ep1, Point3d & ep2)
{
  int i,j, mi;
  double lam[3],help,min,m,dist;
  Point3d p,p0,p1,p2,cp,point,returnpoint,inp;
  Vec3d edge,geomedge,n0,n1,n2,np,dist_vec;

  // Punkt auf der Startfront. Aufgrund der neuen LineDisc
  // liegt die aktuelle Front ausserhalb der Dreiecke
  // Finde Dreieck, mit geringstem Abstand
  // Kante-Punkt

  mi = 0;
  min = 1000000.0;

  inp = Center(ep1,ep2);
  Project_Point2Surface(inp,p);
  for(i=1; i<=geomelements.Size(); i++)
  {
    for(j=0; j<3; j++)
    {
      p.X() = cp.X();
      p.Y() = cp.Y();
      p.Z() = cp.Z();
      p0 = geompoints[geomelements[i].PNum(j%3+1)].p;
      p1 = geompoints[geomelements[i].PNum((j+1)%3+1)].p;

      m = ( (p1 - p0) * (p - p0) ) / ( (p1 - p0) * (p1 - p0) );

      point = p0 + m * (p1 - p0);

      dist_vec = p - point;
      dist = dist_vec.Length();

      if( (dist<min) /*&& (m>=0.0) && (m<=1.0) */)
      {
        mi = i;
        min = dist;
        returnpoint.X()=point.X();
        returnpoint.Y()=point.Y();
        returnpoint.Z()=point.Z();
      }
    }
  }

  help = Dist(returnpoint, cp);

  if(min > 0.2*hh)
    return(0);
  else
    return(mi);
}

int GetEdgeId(const Point3d & ep1, Point3d & ep2)
{
  Vec3d e0,e1,e2,n0,n1,n2,hp,edge,geomedge,np,vec;
  Point3d p0,p1,p2,pold,p,em;
  double lam[3],prod;
  int mi ,i,j,triang[10][4],trnb, lam_i;
  double help, min, min_i,min_lam;
  InputElement es1,es2;
  int nn1,nn2,np1,np2;
  double fall1,fall2;
  int fall1_i,fall2_i;

  mi = 0;
  mi = case1(ep1, ep2);
  if(mi==0)
    mi = case2(ep1, ep2);
  if(mi==0)
    mi = case3(ep1, ep2);
  if(mi==0)
    mi = case4(ep1, ep2);
  if(mi==0)
  {
    cout << ep1 << endl;
    cout << ep2 << endl;

    printf("%s\n","schotter");
  }
  return(mi);
}

Vec3d NormalVector(InputElement e)
{
  Vec3d e1,e2,e3,n1,n2,n3,np;
  Point3d ep1,ep2,ep3,pold,p;
  double min,test,dist;
  int mi,i;

  ep1 = geompoints[e.PNum(1)].p;
  ep2 = geompoints[e.PNum(2)].p;
  ep3 = geompoints[e.PNum(3)].p;

  n1 = ep3 - ep1;
  n1 /= n1.Length();
  n2 = ep3 - ep2;
  n2 /= n2.Length();
  n3 = Cross(n1,n2);
  n3 /= n3.Length();

  return(n3);
}

double Calc_Angle(InputElement e1, InputElement e2)
{
  Vec3d n1, n2;
  float sp;

  n1 = NormalVector(e1);
  n2 = NormalVector(e2);

  sp = n1*n2;
  return(acos(sp));
}

int Test_Line(Point3d p1, Point3d p2, Point3d sp1, Point3d sp2, double xh)
{
  Point3d midp,midsp,ep0,ep1,ep2;
  Vec3d n0,n1,n2,np,front_vec, triang_direction,dist_vec;
  int front_id, line_id;
  double angle,help1,help2,nh,sp,help,L;

  midp = Center (p1,p2);
  midsp = Center (sp1,sp2);

  nh = Dist(sp1,sp2);
  nh = xh;
  if(Dist (midp, midsp) > nh*2/3)
    return(0);
  else
  {

    front_id = GetEdgeId(sp1,sp2);
    line_id = GetEdgeId(p1,p2);

    ep0 = geompoints[geomelements[front_id].PNum(1)].p;
    ep1 = geompoints[geomelements[front_id].PNum(2)].p;
    ep2 = geompoints[geomelements[front_id].PNum(3)].p;

    Calc_Vectors(ep0,ep1,ep2,n0,n1,n2);

    help1 = Project2Plane(p1,np,ep0,n0,n1,n2);
    help2 = Project2Plane(p2,np,ep0,n0,n1,n2);


    // Angle of Visibility

    front_vec = sp2 - sp1;
    triang_direction = Cross(n2,front_vec);
    dist_vec = midp - midsp;

    if(dist_vec*triang_direction<=0.0)
      return(0);





    help = Project2Plane(midp,np,ep0,n0,n1,n2);
    // vielleicht Dist(sp1,p1)
    L = Dist(p1,p2);

    if( (help1>0.5*nh/3) || (help2>0.5*nh/3) )
      return(0);

    angle = Calc_Angle(geomelements[front_id], geomelements[line_id]);
    if(angle>=3.1415*Triangle_Angle2/180)
      return(0);

    if(front_id==line_id)
    {
      if(help > 0.3*nh/3)
      {
        return(0);
      }
    }
    /*	  else
                    if(angle<help/L)
                    return(0);*/

    return(1);
  }
}





void Smooth_SurfaceMesh (ARRAY<Point3d> & points, const ARRAY<Element> & elements, int numboundarypoints, int steps)
{
  int i,j,k,l,m,num_neighbor_points;
  TABLE<INDEX> pointlist(points.Size());
  Point3d p_new;

  for (i = 1; i <= elements.Size(); i++)
  {
    for (j = 1; j <= elements[i].NP(); j++)
      pointlist.Add(elements[i].PNum(j),i);
  }

  if(LGM_DEBUG)
  {
    cout << points.Size() << " " << "points" << endl;
    for (i = 1; i <= points.Size(); i++)
      cout << points[i] << endl;
    cout << elements.Size() << " " << "elements" << endl;
    for (i = 1; i <= elements.Size(); i++)
      cout << elements[i].PNum(1) << "  " << elements[i].PNum(2) << "  " << elements[i].PNum(3) << endl;
    for (i = 1; i <= points.Size(); i++)
    {
      for (j = 1; j <= pointlist.EntrySize(i); j++)
        cout << pointlist.Get(i, j) << "  ";
      cout << endl;
    }
    cout << endl;
  }

  for (k = 1; k <= steps; k++)
  {
    for (i = numboundarypoints+1; i <= points.Size(); i++)
    {
      num_neighbor_points = 0;
      p_new.X() = 0.0;
      p_new.Y() = 0.0;
      p_new.Z() = 0.0;
      for (j = 1; j <= pointlist.EntrySize(i); j++)
      {
        for (l = 1; l <= elements[pointlist.Get(i,j)].NP(); l++)
        {
          for (m = 1; m <= elements[pointlist.Get(i,j)].NP(); m++)
            if(i!=elements[pointlist.Get(i,j)].PNum(m))
            {
              p_new.X() = p_new.X() + points[elements[pointlist.Get(i,j)].PNum(m)].X();
              p_new.Y() = p_new.Y() + points[elements[pointlist.Get(i,j)].PNum(m)].Y();
              p_new.Z() = p_new.Z() + points[elements[pointlist.Get(i,j)].PNum(m)].Z();
              num_neighbor_points++;
            }
        }
      }
      if(num_neighbor_points>0)
      {
        p_new.X() = p_new.X() / num_neighbor_points;
        p_new.Y() = p_new.Y() / num_neighbor_points;
        p_new.Z() = p_new.Z() / num_neighbor_points;
      }
      Project_Point2Surface(p_new, points[i]);
    }
  }
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

  Smooth_SurfaceMesh(points, elements, nbp, smooth);

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
  geompoints.SetSize(0);
  geomelements.SetSize(0);
  return 0;
}
