// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <fstream.h>
#include <math.h>

#include <template.hh>
#include <array.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>
#include <meshing/ruler2.hh>


int ApplyRules ( const ARRAY<netrule*> & rules,
                 ARRAY<Point2d> & lpoints, ARRAY<ILINE> & llines,
                 ARRAY<Element> & elements,
                 ARRAY<INDEX> & dellines, int tolerance)
{
  return(0);
}
