// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   global.cc                                                    */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <strstream.h>
#include <math.h>
#include <new.h>
#include <ctype.h>
#include <string.h>
#include <time.h>


#include <template.hh>
#include <array.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#ifdef MYGRAPH
#include <geom/rot3d.hh>
ROT3D rot;
#endif

ostream * testout;
