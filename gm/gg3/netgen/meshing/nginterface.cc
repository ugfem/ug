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



#include <template.hh>
#include <array.hh>
// #include <spbita2d.hh>
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



ARRAY<POINT3D> * points;
ARRAY<ELEMENT> * volelements;
int nbp;              // number of boundary points


class my_meshing3 : public meshing3
{
public:
  my_meshing3 (char * rulefilename);
  virtual int SavePoint (const POINT3D & p);
  virtual void SaveElement (const ELEMENT & elem);
};



my_meshing3 :: my_meshing3 (char * rulefilename)
  : meshing3 (rulefilename)
{
  ;
}

int my_meshing3 :: SavePoint (const POINT3D & p)
{
  return points -> Append (p);
}



void my_meshing3 :: SaveElement (const ELEMENT & elem)
{
  volelements -> Append (elem);
  UserWriteF("%4d",volelements -> Size());
  if (volelements -> Size()%10 == 0) UserWriteF ("\n");
}


static my_meshing3 * meshing;


int AddSurfaceNode (int nodeid, double x, double y, double z)
{
  points -> Append (POINT3D(x, y, z));
  meshing -> AddPoint (POINT3D(x, y, z), nodeid+1);
  return 0;
}


int AddSurfaceTriangle (int node0, int node1, int node2)
{
  ELEMENT elem(3);

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

  points = new ARRAY<POINT3D>;
  volelements = new ARRAY<ELEMENT>;
  nbp = 0;

  return 0;
}

int StartNetgen (double h, int smooth)
{
  int i;

  nbp = points->Size();

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
