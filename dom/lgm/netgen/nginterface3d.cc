// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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
#include <meshing/global.hh>
#include <meshing/meshing3.hh>
#include <meshing/meshing2.hh>
#include "lgm.hh"

extern "C"
{
        #include "../../gginterface.h"
}

extern "C"
void UserWriteF (char * ch, ...);

static ARRAY<geompoint3d> geompoints;
static ARRAY<splinesegment3d*> splines;
static ARRAY<Point3d> points2;
static ARRAY<int> lp1, lp2;
static ARRAY<Element> elements;
static ARRAY<InputElement> geomelements;
static int nbp;
static Point3d tripoint;
extern int yyparse ();

static ARRAY<Point3d> points;
static ARRAY<Element> * volelements;
static int disp;
static double vol0;
static int write;
static int DD;

static ARRAY<Point3d> locpoints;
static ARRAY<Point2d> plainpoints;
static ARRAY<ILINE> loclines;
static const char * rname;
static int cntelem;
static int oldnl;
static int qualclass, surfind;
static int LGM_DEBUG = 1;
static double SMALL = 0.0005;
static double Triangle_Angle2 = 40.0;
static Point3d sp1;
static Vec3d n, t1, t2;
static ARRAY<INDEX> locelements(0);
Point3d gl_sp1,gl_sp2,gl_p1,gl_p2;


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

    ep1 = geompoints[geomelements[mi].PNum(1)].p;
    ep2 = geompoints[geomelements[mi].PNum(2)].p;
    ep3 = geompoints[geomelements[mi].PNum(3)].p;

    eps = 0.0002;

    n1 = ep3 - ep1;
    if(n1.Length()>eps)
      n1 /= n1.Length();
    n2 = ep3 - ep2;
    if(n1.Length()>eps)
      n2 /= n2.Length();
    n3 = Cross(n1,n2);
    if(n1.Length()>eps)
      n3 /= n3.Length();

    p1.X() = 0.0;
    p1.Y() = 0.0;
    p1.Z() = 0.0;
    np = p - ep1;
    np = np - (n3 * np) * n3;
    p.X() = np.X() + ep1.X();
    p.Y() = np.Y() + ep1.Y();
    p.Z() = np.Z() + ep1.Z();
    dist = sqrt(  (p.X() - pold.X()) * (p.X() - pold.X())
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
  int mi;

  mi = GetEdgeId(p1,p2);
  Calc_Coord_Vectors(p1,p2,mi,ex,ey,ez);
  globp1 = p1;

  return(mi);
}

void surfacemeshing :: TransformToPlain (INDEX surfind, const Point3d & locpoint,
                                         Point2d & plainpoint, double h)
{
  Vec3d p1p,pp;

  p1p = locpoint - globp1;
  pp =  (p1p * ez) * ez;
  p1p = p1p - pp;
  p1p /= h;
  plainpoint.X() = p1p * ex;
  plainpoint.Y() =  p1p * ey;

}

void ProjectPoint2geo (Point3d & p)
{
  Vec3d n1,n2,n3,e1,e2,e3,np;
  Point3d ep1,ep2,ep3,pnew,p1;
  int mi;

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

void surfacemeshing :: TransformFromPlain (INDEX surfind, Point2d & plainpoint,
                                           Point3d & locpoint, double h)
{
  Point3d locpoint1,savepoint;
  Vec3d p1p;

  p1p = plainpoint.X() * ex + plainpoint.Y() * ey;
  p1p *= h;
  locpoint1 = globp1 + p1p;
  savepoint = locpoint1;

  // Projeziere Punkt auf die Oberflaeche
  Project_Point2Surface_2(locpoint1, locpoint, ez);
}

void surfacemeshing :: ProjectPoint (INDEX surfind, Point3d & p) const
{
  Vec3d e1,e2,e3,n1,n2,n3,np;
  Point3d ep1,ep2,ep3,p1,pold,pnew;
  double dist;
  int mi;

  mi = GetTriangleId(p,pnew);
  p = pnew;

  dist = sqrt(  (p.X() - pold.X()) * (p.X() - pold.X())
                + (p.Y() - pold.Y()) * (p.Y() - pold.Y())
                + (p.Z() - pold.Z()) * (p.Z() - pold.Z()));
  while( dist >= 0.5 ) ;
}

void surfacemeshing :: ProjectPoint2 (INDEX surfind, INDEX surfind2, Point3d & p) const
{}

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
}

double surfacemeshing :: CalcLocalH (const Point3d & p, int surfind, double gh) const
{
  return gh;
}

int surfacemeshing :: Write_Surface_Grid ()
{
  int i,j;
  FILE *file;

  file = fopen("surface","w");
  fprintf(file,"%s\n","surfacemesh");
  fprintf(file,"%d\n",points.Size());
  for(i=1; i<=points.Size(); i++)
    fprintf(file,"%lf %lf %lf\n",points[i].X(),points[i].Y(),points[i].Z());
  fprintf(file,"%d\n",elements.Size());
  for(i=1; i<=elements.Size(); i++)
    fprintf(file,"%d %d %d\n",elements[i].PNum(1)-1,elements[i].PNum(2)-1,elements[i].PNum(3)-1);
  fclose(file);
  return(0);
}

int surfacemeshing :: SavePoint (const Point3d & p)
{
  return points.Append (p);
}

void surfacemeshing :: SaveElement (const Element & elem)
{
  elements.Append (elem);
}

void surfacemeshing :: Write2Shell (int n)
{
  if(n%DD==0)
  {
    if(write%10==0)
      UserWriteF("\n");
    UserWriteF("%s%d%s","[",n,"]");
    write++;
  }
}

class my_surfacemeshing : public surfacemeshing
{
public:
  my_surfacemeshing (char * rulefilename);
  virtual int SavePoint (const Point3d & p);
  virtual void SaveElement (const Element & elem);
};

my_surfacemeshing :: my_surfacemeshing (char * rulefilename)
  : surfacemeshing (rulefilename)
{}

int my_surfacemeshing :: SavePoint (const Point3d & p)
{
  return points.Append (p);
}

void my_surfacemeshing :: SaveElement (const Element & elem)
{
  elements.Append (elem);
}

static my_surfacemeshing * meshing;
static const surfacemeshing * meshthis;

int AddGeomPoint (int id, double x, double y, double z)
{
  geompoints.Append (geompoint3d());
  geompoints.Last().p.X() = x;
  geompoints.Last().p.Y() = y;
  geompoints.Last().p.Z() = z;

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

  return 0;
}

int AddLinePoint (int id, double x, double y, double z)
{
  points.Append (Point3d(x,y,z));

  meshing->AddPoint (points.Last(), id);
  return 0;
}

int AddLineSegment (int i1,int i2)
{
  meshing->AddBoundaryElement (i1, i2, 1);
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
{}

static double SMALLDOUBLE = 1.0e-6;

static int Cross_Check(int dummy)
{
  int i,j;
  double xRi,yRi,xQi,yQi,xRj,yRj,xQj,yQj,denominator,lambdai,lambdaj;
  for(i=1; i<=loclines.Size(); i++)
  {
    for(j=1; j<=loclines.Size(); j++)
    {
      if(i!=j)
      {
        xRi = plainpoints[loclines[i].I1()].X();
        yRi = plainpoints[loclines[i].I1()].Y();
        xQi = plainpoints[loclines[i].I2()].X();
        yQi = plainpoints[loclines[i].I2()].Y();
        xRj = plainpoints[loclines[j].I1()].X();
        yRj = plainpoints[loclines[j].I1()].Y();
        xQj = plainpoints[loclines[j].I2()].X();
        yQj = plainpoints[loclines[j].I2()].Y();

        denominator = (xRi-xQi)*(yQj-yRj) - (yRi-yQi)*(xQj-xRj);
        if(fabs(denominator)<SMALLDOUBLE)
          continue;

        lambdai = ( (yQj-yRj)*(xRi-xRj) + (xRj-xQj)*(yRi-yRj) ) / denominator;
        lambdaj = ( (yQi-yRi)*(xRi-xRj) + (xRi-xQi)*(yRi-yRj) ) / denominator;

        if( (lambdai<1.0-2*SMALLDOUBLE) && (lambdai>2*SMALLDOUBLE) && (lambdaj<1.0-2*SMALLDOUBLE) && (lambdaj>2*SMALLDOUBLE) )
        {
          return(1);
        }
      }
    }
  }
  return(0);
}


void surfacemeshing :: Mesh (double gh)
{
  ARRAY<INDEX> pindex, lindex;
  ARRAY<int> delpoints, dellines;
  ARRAY<Element> locelements;
  int i, j, oldnp;
  char ch, found;
  INDEX globind;
  double in[5];
  static ARRAY<Point3d> locp;
  FILE *file;
  double h,dist;
  Point3d bemp, bemp1, bemp2;

  if(LGM_DEBUG)
  {
    // Ausgabe fuer Netgen
    file = fopen("netgen","w");
    fprintf(file, "%s\n","ug");
    fprintf(file, "%d\n",geompoints.Size());
    for (i = 1; i <= geompoints.Size(); i++)
      fprintf(file, "%f %f %f\n", geompoints[i].p.X(),
              geompoints[i].p.Y(),
              geompoints[i].p.Z());
    fprintf(file, "%d\n",geomelements.Size());
    for(i=1; i<=geomelements.Size(); i++)
      fprintf(file, "%d %d %d %d %d %d\n",    geomelements[i].PNum(1),
              geomelements[i].PNum(2),
              geomelements[i].PNum(3),
              geomelements[i].Neighbour(1),
              geomelements[i].Neighbour(2),
              geomelements[i].Neighbour(3));

    fclose(file);

    //		 Ausgabe der Front
    adfront->ugPrint(cout);
  }

  testmode = 0;
  meshthis = this;

  StartMesh();

  adfront ->SetStartFront ();
  while (ch != 'e' && !adfront ->Empty())
  {
    //		adfront->Print(cout);

    dist = 1.5;
    do
    {
      locpoints.SetSize(0);
      loclines.SetSize(0);
      pindex.SetSize(0);
      lindex.SetSize(0);

      h = gh;
      if (h>0.0)
        qualclass =
          adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                               surfind, 3 * h, dist);
      else
        qualclass =
          adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                               surfind, -3 * h, dist);

      oldnp = locpoints.Size();
      oldnl = loclines.Size();

      bemp1.X() = locpoints[loclines[1].I1()].X();
      bemp1.Y() = locpoints[loclines[1].I1()].Y();
      bemp1.Z() = locpoints[loclines[1].I1()].Z();
      bemp2.X() = locpoints[loclines[1].I2()].X();
      bemp2.Y() = locpoints[loclines[1].I2()].Y();
      bemp2.Z() = locpoints[loclines[1].I2()].Z();

      bemp = Center (bemp1, bemp2);
      bemp.Y() += 0.1 * (bemp2.X() - bemp1.X());
      bemp.X() -= 0.1 * (bemp2.Y() - bemp1.Y());

      if(gh<=0.0)
      {
        locpoints.SetSize(0);
        loclines.SetSize(0);
        pindex.SetSize(0);
        lindex.SetSize(0);
        in[0] = bemp.X();
        in[1] = bemp.Y();
        in[2] = bemp.Z();
        in[3] = gh;
        Get_Local_h(in,&h);
        qualclass =
          adfront ->GetLocals (locpoints, loclines, pindex, lindex,
                               surfind, - 3 * gh, dist);
      }

      DefineTransformation (surfind, locpoints[loclines[1].I1()],
                            locpoints[loclines[1].I2()]);

      plainpoints.SetSize (locpoints.Size());
      for (i = 1; i <= locpoints.Size(); i++)
      {
        TransformToPlain (surfind, locpoints[i], plainpoints[i], h);
      }
      dist = dist / 2;
    }
    while(Cross_Check(0));

    for (i = 2; i <= loclines.Size(); i++)              // don't remove first line
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

    found = GenerateTriangle (plainpoints, loclines, locelements,
                              dellines, h);
    /*    found = ApplyRules (rules, plainpoints, loclines, locelements,
                            dellines, qualclass);*/

    if (found)
    {
      //			adfront->Print(cout);
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

      for (i = oldnl+1; i <= loclines.Size(); i++)
      {
        adfront -> AddLine (pindex[loclines[i].I1()],
                            pindex[loclines[i].I2()], surfind);
      }

      for (i = 1; i <= locelements.Size(); i++)
      {
        for (j = 1; j <= locelements[i].NP(); j++)
        {
          locelements[i].PNum(j) =
            adfront -> GetGlobalIndex (pindex[locelements[i].PNum(j)]);
        }

        locelements[i].SetSurfaceIndex (surfind);

        SaveElement (locelements[i]);
        Write2Shell(cntelem);
        //				Write_Surface_Grid();
        cntelem++;
        //				printf("%d %d %d %d\n",cntelem,locelements[i].PNum(1),locelements[i].PNum(2),
        //						locelements[i].PNum(3));
      }

      for (i = 1; i <= dellines.Size(); i++)
      {
        adfront -> DeleteLine (lindex[dellines[i]]);
      }

      //			rname = rules[found]->Name();
      //			adfront->GetPoints(locp);
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
    plainpoints.SetSize(0);
  }

  EndMesh ();
}

#define Det(a)                                  ( - a[0] * a[4] * a[8]  \
                                                  + a[0] * a[5] * a[7]    \
                                                  + a[3] * a[1] * a[8]    \
                                                  - a[3] * a[2] * a[7]    \
                                                  - a[6] * a[1] * a[5]    \
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

#define Det2d(a)                                ( a[0] * a[3] - a[1] * a[2] )

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

  det = ABS(Det2d(aa));
  if(ABS(det)>1e-8)
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

    det = ABS(Det2d(aa));
    if(ABS(det)>1e-8)
    {
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

      det = ABS(Det2d(aa));
      if(ABS(det)>1e-8)
      {
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

double Project2Plane(   Point3d & p,
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

int InitSurfaceNetgen (char * rulefilename)
{
  meshing = new my_surfacemeshing(rulefilename);

  volelements = new ARRAY<Element>;
  nbp = 0;

  return 0;
}

int StartSurfaceNetgen (double h, int smooth, int display,int D)
{
  int i;
  double hpx, hpy, hpz, v1x, v1y, v1z, v2x, v2y, v2z;
  Vec3d n;
  Point3d p1,p2,p3,p,p_in,p_out;

  nbp = points.Size();
  disp = display;
  write = 0;
  DD = D;

  meshing -> Mesh (h);

  Smooth_SurfaceMesh(points, elements, nbp, smooth);

  Allocate_Mem_Surfdisc(points.Size(), elements.Size());

  for (i = 1; i <= points.Size(); i++)
  {
    p_in = points.Get(i);
    p_out.X() = 0.0;
    p_out.Y() = 0.0;
    p_out.Z() = 0.0;
    n.X() = 0.0;
    n.Y() = 0.0;
    n.Z() = 0.0;
    Project_Point2Surface_2(p_in, p_out, n);
    AddInnerNode2ug (p_in.X(),p_in.Y(),p_in.Z());
    /*		AddInnerNode2ug (	points.Get(i).X(),
                                                            points.Get(i).Y(),
                                                            points.Get(i).Z());*/
  }
  for (i = 1; i <= elements.Size(); i++)
  {
    AddSurfaceTriangle2ug ( elements.Get(i).PNum(1) - 1,
                            elements.Get(i).PNum(2) - 1,
                            elements.Get(i).PNum(3) - 1);

    if(LGM_DEBUG)
    {
      p1 = points.Get(elements.Get(i).PNum(1));
      p2 = points.Get(elements.Get(i).PNum(2));
      p3 = points.Get(elements.Get(i).PNum(3));

      v1x = points.Get(elements.Get(i).PNum(2)).X() - points.Get(elements.Get(i).PNum(1)).X();
      v1y = points.Get(elements.Get(i).PNum(2)).Y() - points.Get(elements.Get(i).PNum(1)).Y();
      v1z = points.Get(elements.Get(i).PNum(2)).Z() - points.Get(elements.Get(i).PNum(1)).Z();

      v2x = points.Get(elements.Get(i).PNum(3)).X() - points.Get(elements.Get(i).PNum(1)).X();
      v2y = points.Get(elements.Get(i).PNum(3)).Y() - points.Get(elements.Get(i).PNum(1)).Y();
      v2z = points.Get(elements.Get(i).PNum(3)).Z() - points.Get(elements.Get(i).PNum(1)).Z();

      n.X() = v1y * v2z - v1z * v2y;
      n.Y() = v1z * v2x - v1x * v2z;
      n.Z() = v1x * v2y - v1y * v2x;
      n /= n.Length();
    }
  }

  elements.SetSize(0);
  points.SetSize(0);
  geompoints.SetSize(0);
  geomelements.SetSize(0);
  return 0;
}

int c11(const Point3d & point, Point3d & returnpoint, Vec3d n)
{
  int i,j,triang[10][4],trnb, mi, nn1, nn2, np1, np2;
  double lam[3],prod,help, angle,sp;
  Point3d p0,p1,p2, p;
  Vec3d edge,geomedge,n0,n1,n2,np, ne;
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
    ne = NormalVector(geomelements[i]);
    if(n.Length()>SMALL)
    {
      sp = ne*n;
      if(sp>1-SMALL)
        angle = 0.0;
      else
        angle = acos(sp);
    }
    else
      angle = 0;
    if(angle<3.1415*Triangle_Angle2/180)
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
  }

  return(mi);
}

int c22(Point3d & p, Point3d & returnpoint, double &min, Vec3d n)
{
  int i,j, mi;
  double lam[3],help, angle,sp;
  Point3d p0,p1,p2;
  Vec3d n0,n1,n2,np, ne;

  /* Punkt liegt nicht auf der Surface
   * Projeziere den Punkt in ein Dreieck mit  0<=lamda_i<=1 sodass der
   * Abstand zur Surface minimal wird unter der
   * Bedingung, dass der Normalenvektor mit der vorgegebenen Richtung
   * uebereinstimmt
   */

  min = 10000000.0;
  mi = 0;
  for(i=1; i<=geomelements.Size(); i++)
  {
    ne = NormalVector(geomelements[i]);
    if(n.Length()>SMALL)
    {
      sp = ne*n;
      if(sp>1-SMALL)
        angle = 0.0;
      else
        angle = acos(sp);
    }
    else
      angle = 0;
    if(angle<3.1415*Triangle_Angle2/180)
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
  }
  return(mi);
}

int c33(Point3d & p, Point3d & returnpoint, double &min, Vec3d n)
{
  int i,j, mi;
  double lam[3],help, d0, d1, dist, m, angle,sp;
  Point3d p0,p1,p2, point;
  Vec3d n0,n1,n2,np, dist_vec, ne;

  /* Punkt liegt nicht auf der Surface
   * Projeziere den Punkt auf die Kante,  sdass der
   *Abstand minimal wird */

  min = 10000000.0;

  for(i=1; i<=geomelements.Size(); i++)
  {
    ne = NormalVector(geomelements[i]);
    if(n.Length()>SMALL)
    {
      sp = ne*n;
      if(sp>1-SMALL)
        angle = 0.0;
      else
        angle = acos(sp);
    }
    else
      angle = 0;
    if(angle<3.1415*Triangle_Angle2/180)
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
  }
  return(mi);
}

int Project_Point2Surface_2(Point3d &inpoint, Point3d &outpoint, Vec3d n)
{
  int i, j, mi,mi2,mi3;
  double min2, min3;
  Point3d p1, p2, p3;
  mi = 0;
  mi = c11(inpoint, p1, n);
  if(mi!=0)
    outpoint = p1;
  else
  {
    mi2 = c22(inpoint, p2, min2, n);
    mi3 = c33(inpoint, p3, min3, n);
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

int Get_unique_Tangentialplane(const Point3d & point)
{
  int i,j,triang[10][4],trnb, mi, nn1, nn2, np1, np2,nb,list[10];
  double lam[3],prod,help,angle;
  Point3d p0,p1,p2, p;
  Vec3d edge,geomedge,n0,n1,n2,np,n_old,n_new;
  InputElement es;

  // Punkt liegt auf der Surface
  // liefere Inputdreieck zurueck, das eine sinnvolle Tangentialebene definiert

  p = point;
  mi = 0;
  nb = 0;
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
      list[nb] = i;
      nb++;
    }
  }
  if(nb==0)
    return(0);
  else
  {
    for(i=0; i<nb; i++)
      for(j=i+1; j<nb; j++)
      {
        angle = Calc_Angle(geomelements[list[i]], geomelements[list[j]]);
        if(angle>=3.1415*Triangle_Angle2/180)
          return(0);
      }
    return(list[0]);
  }
}

int case1(const Point3d & ep1, Point3d & ep2)
{
  int i,j,triang[10][4],trnb, mi, nn1, nn2, np1, np2;
  double lam[3],prod,help,det;
  Point3d p,p0,p1,p2,inp,e1,e2;
  Vec3d edge,geomedge,n0,n1,n2,np;
  InputElement es;
  int mi1,mi2;
  double angle1,angle2;

  // Punkt liegt genau auf der Kante zwischen 2 Dreiecken
  // passiert fuer erkannte Kanten auf der Surface

  mi = 0;
  trnb = 0;
  inp = Center(ep1,ep2);
  p = inp;
  // Projektion auf dei Surface faellt weg
  // Project_Point2Surface(inp,p);

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

    if( (lam[0]>=-SMALL) && (lam[1]>=-SMALL) && (lam[2]>=-SMALL) && (help<SMALL) )
    {
      triang[trnb][0] = i;
      triang[trnb][1] = -1;
      triang[trnb][2] = -1;
      triang[trnb][3] = -1;
      Calc_Vectors(p0,p1,p2,n0,n1,n2);
      help = Project2Plane(p,np,p0,n0,n1,n2);
      for(j=0; j<3; j++)
        if(lam[j]<10*SMALL)
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
  int i,j, mi, mi1, mi2;
  double lam[3],help,min, angle1, angle2;
  Point3d p,p0,p1,p2,inp;
  Vec3d edge,geomedge,n0,n1,n2,np;

  mi = 0;
  min = 1000000.0;
  inp = Center(ep1,ep2);
  p = inp;

  // Entscheide ob
  // sp1 und sp2 jeweils genau eine Tangetialebene definieren
  mi1 = Get_unique_Tangentialplane(ep1);
  mi2 = Get_unique_Tangentialplane(ep2);

  // Falls das moeglich ist, suche das Inputdreieck, dass genauso orientiert ist
  // und den kuerzesten Abstand zum Mittelpunkt von sp1 und sp2 hat
  if( (mi1!=0) || (mi2!=0) )
  {
    if(mi1==0)
      mi1 = mi2;
    for(i=1; i<=geomelements.Size(); i++)
    {
      angle1 = Calc_Angle(geomelements[mi1], geomelements[i]);
      angle2 = Calc_Angle(geomelements[mi1], geomelements[i]);
      if( (angle1<3.1415*Triangle_Angle2/180) && (angle2<3.1415*Triangle_Angle2/180) )
      {
        p0 = geompoints[geomelements[i].PNum(1)].p;
        p1 = geompoints[geomelements[i].PNum(2)].p;
        p2 = geompoints[geomelements[i].PNum(3)].p;

        Calc_Local_Coordinates(p0,p1,p2,p,lam[0],lam[1],lam[2]);

        Calc_Vectors(p0,p1,p2,n0,n1,n2);
        help = Project2Plane(p,np,p0,n0,n1,n2);

        if(help<min)
          if( (lam[0]>-SMALL) && (lam[0]<1.0+SMALL) )
            if( (lam[1]>-SMALL) && (lam[1]<1.0+SMALL) )
              if( (lam[2]>-SMALL) && (lam[2]<1.0+SMALL) )
              {
                min = help;
                mi = i;
              }

      }
    }
  }
  else
  {
    // Tangentialebenen nicht eindeutig
    // hier muesste aber case 1 greifen
    /*		assert(0);*/
  }
  return(mi);
}

int case3(const Point3d & ep1, Point3d & ep2)
{
  int i,j, mi;
  double lam[3],help,min;
  Point3d p,p0,p1,p2,inp;
  Vec3d edge,geomedge,n0,n1,n2,np,dummy;

  // Punkt auf der Startfront. Aufgrund der neuen LineDisc
  // liegt die aktuelle Front ausserhalb der Dreiecke
  // Finde Dreieck, mit geringstem Abstand
  // Kante-Punkt

  mi = 0;
  min = 1000000.0;
  inp = Center(ep1,ep2);
  dummy.X() = 0.0;
  dummy.Y() = 0.0;
  dummy.Z() = 0.0;

  Project_Point2Surface_2(inp,p,dummy);
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
  return(mi);
}

int case4(const Point3d & ep1, Point3d & ep2)
{
  int i,j, mi;
  double lam[3],help,min,m,dist;
  Point3d p,p0,p1,p2,cp,point,returnpoint,inp;
  Vec3d edge,geomedge,n0,n1,n2,np,dist_vec,dummy;

  // Punkt auf der Startfront. Aufgrund der neuen LineDisc
  // liegt die aktuelle Front ausserhalb der Dreiecke
  // Finde Dreieck, mit geringstem Abstand
  // Kante-Punkt

  mi = 0;
  min = 1000000.0;

  inp = Center(ep1,ep2);
  dummy.X() = 0.0;
  dummy.Y() = 0.0;
  dummy.Z() = 0.0;

  Project_Point2Surface_2(inp,p,dummy);
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
    //cout << ep1 << endl;
    //cout << ep2 << endl;

    printf("%s\n","schotter");
  }
  return(mi);
}

int Test_Line(Point3d p1, Point3d p2, Point3d sp1, Point3d sp2, double xh, double dist)
{
  Point3d midp,midsp,ep0,ep1,ep2;
  Vec3d n0,n1,n2,np,front_vec, triang_direction,dist_vec1,dist_vec2,front_n,edge_n;
  int front_id, line_id;
  double angle,help1,help2,nh,sp,help,L,front;

  gl_sp1 = sp1;
  gl_sp2 = sp2;
  gl_p1 = p1;
  gl_p2 = p2;

  //	Plot3D (PlotFrontLine, 1, 1);

  midp = Center (p1,p2);
  midsp = Center (sp1,sp2);

  nh = Dist(sp1,sp2);
  front = Dist(sp1,sp2);
  nh = xh;
  if(Dist (midp, midsp) > nh)
    return(0);
  else
  {
    front_id = GetEdgeId(sp1,sp2);
    assert(front_id);
    line_id = GetEdgeId(p1,p2);
    assert(line_id);

    front_n = NormalVector(geomelements[front_id]);
    edge_n = NormalVector(geomelements[line_id]);

    if(front_n*edge_n<0.0)
      return(0);

    // definiere lokale Ebene
    ep0 = geompoints[geomelements[front_id].PNum(1)].p;
    ep1 = geompoints[geomelements[front_id].PNum(2)].p;
    ep2 = geompoints[geomelements[front_id].PNum(3)].p;

    Calc_Vectors(ep0,ep1,ep2,n0,n1,n2);

    help1 = Project2Plane(p1,np,ep0,n0,n1,n2);
    help2 = Project2Plane(p2,np,ep0,n0,n1,n2);

    // Angle of Visibility
    front_vec = sp2 - sp1;
    triang_direction = Cross(n2,front_vec);
    dist_vec1 = midp - midsp;
    dist_vec2 = p2 - midsp;

    /*		if((dist_vec1*triang_direction<=0.0)&&(dist_vec2*triang_direction<=0.0))
                            return(0);*/
    if((dist_vec1*triang_direction<=0.0))
      return(0);

    help = Project2Plane(midp,np,ep0,n0,n1,n2);

    if( (help1>dist*xh) || (help2>dist*xh) )
      return(0);

    angle = Calc_Angle(geomelements[front_id], geomelements[line_id]);
    if(angle>=3.1415*Triangle_Angle2/180)
      return(0);
    if(front_id==line_id)
      if(help > 0.3*front)
        return(0);

  }
  return(1);
}
int Test_Line_OLD(Point3d p1, Point3d p2, Point3d sp1, Point3d sp2, double xh)
{
  Point3d midp,midsp,ep0,ep1,ep2;
  Vec3d n0,n1,n2,np,front_vec, triang_direction,dist_vec;
  int front_id, line_id;
  double angle,help1,help2,nh,sp,help,L,front;

  gl_sp1 = sp1;
  gl_sp2 = sp2;
  gl_p1 = p1;
  gl_p2 = p2;

  //	Plot3D (PlotFrontLine, 1, 1);

  midp = Center (p1,p2);
  midsp = Center (sp1,sp2);

  nh = Dist(sp1,sp2);
  front = Dist(sp1,sp2);
  nh = xh;
  if(Dist (midp, midsp) > nh)
    return(0);
  else
  {
    front_id = GetEdgeId(sp1,sp2);
    assert(front_id);
    line_id = GetEdgeId(p1,p2);
    assert(line_id);

    // definiere lokale Ebene
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

    if( (help1>1.5*front) || (help2>1.5*front) )
      return(0);

    angle = Calc_Angle(geomelements[front_id], geomelements[line_id]);
    if(angle>=3.1415*Triangle_Angle2/180)
      return(0);
    if(front_id==line_id)
      if(help > 0.3*front)
        return(0);

  }
  return(1);
}

int Test_Point(Point3d p, Point3d sp1, Point3d sp2, double xh)
{
  Point3d midp,midsp,ep0,ep1,ep2,p1;
  Vec3d n0,n1,n2,np;
  int front_id, line_id,point_id;
  double angle,help1,help2,nh,sp,help;

  midsp = Center (sp1,sp2);

  nh = Dist(sp1,sp2);
  front_id = GetEdgeId(sp1,sp2);

  ep0 = geompoints[geomelements[front_id].PNum(1)].p;
  ep1 = geompoints[geomelements[front_id].PNum(2)].p;
  ep2 = geompoints[geomelements[front_id].PNum(3)].p;

  if( (p.X()==ep0.X())&&(p.Y()==ep0.Y())&&(p.Z()==ep0.Z()) ||
      (p.X()==ep1.X())&&(p.Y()==ep1.Y())&&(p.Z()==ep1.Z()) ||
      (p.X()==ep2.X())&&(p.Y()==ep2.Y())&&(p.Z()==ep2.Z()) )
    return(1);
  else
    return(0);
}

void Smooth_SurfaceMesh (ARRAY<Point3d> & points, const ARRAY<Element> & elements, int numboundarypoints, int steps)
{
  int i,j,k,l,m,num_neighbor_points, mi;
  TABLE<INDEX> pointlist(points.Size());
  Point3d p_new;
  Vec3d n,dummy;

  for (i = 1; i <= elements.Size(); i++)
  {
    for (j = 1; j <= elements[i].NP(); j++)
      pointlist.Add(elements[i].PNum(j),i);
  }

  for (k = 1; k <= steps; k++)
  {
    for (i = numboundarypoints+1; i <= points.Size(); i++)
    {
      dummy.X() = 0.0;
      dummy.Y() = 0.0;
      dummy.Z() = 0.0;
      mi = Project_Point2Surface_2(points[i], points[i],dummy);
      n = NormalVector(geomelements[mi]);
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
      Project_Point2Surface_2(p_new, points[i], n);
    }
  }
}
