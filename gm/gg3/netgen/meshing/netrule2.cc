// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
// #include <conio.h>

#include <template.hh>
#include <array.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>
#include <meshing/ruler2.hh>


netrule :: netrule ()
{}


void netrule :: SetFreeZoneTransformation (const Vector & u)
{}


int netrule :: IsInFreeZone (const Point2d & p) const
{
  return(0);
}

int netrule :: ConvexFreeZone () const
{
  return(0);
}



float netrule :: CalcPointDist (int pi, const Point2d & p) const
{
  return(0.0);
}


float netrule :: CalcLineError (int li, const Vec2d & v) const
{
  return(0.0);
}
