// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <new.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include <general/myadt.hh>

#include <geom/geom2d.hh>
#include <geom/geom3d.hh>
#include <geom/reftrans.hh>

#ifdef MYGRAPH
#include <geom/rot3d.hh>
#include <graphics/mygraph.hh>
#endif

#include <linalg/linalg.hh>

#include <meshing/global.hh>
#include <meshing/meshing3.hh>


// global variables for plotting:

static ARRAY<Point3d> locpoints;
static ARRAY<Element> locfaces;
static ARRAY<INDEX> pindex, findex;
static ARRAY<int> style;
static int drawrad = 100;
static int cntsucc = 0;
static long cnttrials = 0;
static int cntelem = 0;
static char buf[100];
static int qualclass;
static float vol0, err, h;
static int problemindex = 1;


static Meshing3 * meshing;
static int LGM_DEBUG = 0;

static ARRAY<Point3d> grouppoints;
static ARRAY<Element> groupfaces;
static ARRAY<INDEX> grouppindex, groupfindex;

extern int FindInnerPoint (const ARRAY<Point3d> & grouppoints,
                           const ARRAY<Element> & groupfaces,
                           Point3d & p);


void DrawPnF(const ARRAY<Point3d> & pointl, const ARRAY<Element> & facel,
             ARRAY<int> & style, float xh , const ROT3D & r)
{}




void PlotVolMesh (const ROT3D & r, char key)
{}



void Meshing3 :: Mesh (double ah,int prism)
{}
