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


extern "C"
{
#include "../../gginterface.h"
}

int testmode;

extern "C"
void UserWriteF (char * ch, ...);


void MyError (char * ch)
{
  UserWriteF (ch);
}



static ARRAY<Point3d> * points;
static ARRAY<Element> * volelements;
static int nbp,disp;
static double vol0;

class my_meshing3 : public Meshing3
{
public:
  my_meshing3 (char * rulefilename);
  virtual int SavePoint (const Point3d & p);
  virtual void SaveElement (const Element & elem);

  friend int StartNetgen (double h, int smooth, int display);
};



my_meshing3 :: my_meshing3 (char * rulefilename)
  : Meshing3 (rulefilename)
{
  ;
}

int my_meshing3 :: SavePoint (const Point3d & p)
{
  return points -> Append (p);
}

#define V3_ADD(A,B,C)                              {(C)[0] = (A)[0] + (B)[0];\
                                                    (C)[1] = (A)[1] + (B)[1];\
                                                    (C)[2] = (A)[2] + (B)[2];}

#define V3_SCALE(c,C)                              {(C)[0] = (c)*(C)[0];\
                                                    (C)[1] = (c)*(C)[1];\
                                                    (C)[2] = (c)*(C)[2];}

#define V3_EUKLIDNORM_OF_DIFF(A,B,b)    (b) = (sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])+((A)[2]-(B)[2])*((A)[2]-(B)[2]))));

#define V3_CLEAR(A)                                {(A)[0] = 0.0; (A)[1]= 0.0; (A)[2] = 0.0;}

void my_meshing3 :: SaveElement (const Element & elem)
{
  float x[4][3],diam,fac,global[3],inndiam,dist,percent;
  int i,n;

  if (disp)
  {
    n = 4;
    x[0][0] = points -> Get(elem.PNum(1)).X();
    x[0][1] = points -> Get(elem.PNum(1)).Y();
    x[0][2] = points -> Get(elem.PNum(1)).Z();
    x[1][0] = points -> Get(elem.PNum(2)).X();
    x[1][1] = points -> Get(elem.PNum(2)).Y();
    x[1][2] = points -> Get(elem.PNum(2)).Z();
    x[2][0] = points -> Get(elem.PNum(3)).X();
    x[2][1] = points -> Get(elem.PNum(3)).Y();
    x[2][2] = points -> Get(elem.PNum(3)).Z();
    x[3][0] = points -> Get(elem.PNum(4)).X();
    x[3][1] = points -> Get(elem.PNum(4)).Y();
    x[3][2] = points -> Get(elem.PNum(4)).Z();

    V3_CLEAR(global);
    for (i=0; i<n; i++)
      V3_ADD(x[i],global,global);
    fac = 1.0 / n;
    V3_SCALE(fac,global);
    diam = 0.0;
    inndiam = 100000000.0;
    for (i=0; i<n; i++)
    {
      V3_EUKLIDNORM_OF_DIFF(x[i],global,fac);
      if (fac < inndiam)
        inndiam = fac;
      if (fac > diam)
        diam = fac;
    }
    dist = sqrt(global[0]*global[0]+global[1]*global[1]+global[2]*global[2]);

    volelements -> Append (elem);

    percent = 100.0 * adfront->Volume() / vol0;

    UserWriteF(" ID(Elem)=%4d midPoint %6.2f %6.2f %6.2f dist %6.2f diam %6.2f %6.2f vol %6.3f\%\n",
               volelements -> Size(),global[0],global[1],global[2],
               dist,inndiam,diam,percent);
  }
}


static my_meshing3 * meshing;


int AddSurfaceNode (int nodeid, double x, double y, double z)
{
  points -> Append (Point3d(x, y, z));
  meshing -> AddPoint (Point3d(x, y, z), nodeid+1);
  return 0;
}


int AddSurfaceTriangle (int node0, int node1, int node2)
{
  Element elem(3);

  elem.PNum(1) = node0 + 1;
  elem.PNum(2) = node1 + 1;
  elem.PNum(3) = node2 + 1;
  elem.SetSurfaceIndex (1);
  meshing -> AddBoundaryElement (elem, 0);

  return 0;
}


int InitNetgen (char * rulefilename)
{
  testout = new ofstream("test.out");
  meshing = new my_meshing3(rulefilename);

  points = new ARRAY<Point3d>;
  volelements = new ARRAY<Element>;
  nbp = 0;

  return 0;
}

int StartNetgen (double h, int smooth, int display)
{
  int i;

  nbp = points->Size();
  disp = display;
  vol0 = meshing->adfront->Volume();

  meshing -> Mesh (h);

  for (i=0; i<smooth; i++)
    meshing -> ImproveMesh (*points, *volelements, nbp, h);

  UserWriteF("\n");
  for (i = nbp + 1; i <= points -> Size(); i++)
    AddInnerNode (points -> Get(i).X(),
                  points -> Get(i).Y(),
                  points -> Get(i).Z());

  for (i = 1; i <= volelements -> Size(); i++)
  {
    AddTetrahedron (volelements -> Get(i).PNum(1) - 1,
                    volelements -> Get(i).PNum(2) - 1,
                    volelements -> Get(i).PNum(3) - 1,
                    volelements -> Get(i).PNum(4) - 1);
  }
  return 0;
}
