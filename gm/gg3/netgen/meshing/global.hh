// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_GLOBAL
#define FILE_GLOBAL

#include <meshing/meshtype.hh>

/**************************************************************************/
/* File:   global.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
   global functions and variables
 */


class ROT3D;
extern ROT3D rot;

extern int GetTime ();
extern int testmode;

class ostream;
extern ostream * testout;

#endif
