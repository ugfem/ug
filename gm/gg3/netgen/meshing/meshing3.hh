// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   meshing3.hh                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/

#ifndef FILE_MESHING3
#define FILE_MESHING3

#include <meshing/meshtype.hh>
#include <meshing/adfront3.hh>
#include <meshing/ruler3.hh>

class meshing3
{
  ADFRONT3 * adfront;
  ARRAY<vnetrule*> rules;
  ARRAY<int> ruleused;
  ARRAY<char*> problems;
  double tolfak;

public:
  meshing3 (char * rulefilename);
  virtual ~meshing3 ();

  void LoadRules (char * filename);
  void Mesh (double gh);

  void ImproveMesh (ARRAY<POINT3D> & points, const ARRAY<ELEMENT> & surfelements,
                    const ARRAY<ELEMENT> & elements,
                    double h);

  void ImproveMesh (ARRAY<POINT3D> & points,
                    const ARRAY<ELEMENT> & elements,
                    int nboundnodes,
                    double h);

  void AddPoint (const POINT3D & p, INDEX globind);
  void AddBoundaryElement (const ELEMENT & elem, int inverse);

  virtual int SavePoint (const POINT3D & p) { return 0; }
  virtual void SaveElement (const ELEMENT & elem) { };

  friend void PlotVolMesh (const ROT3D & r, char key);
  friend void TestRules ();
};
#endif
