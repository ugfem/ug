// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   adfront3.hh                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/

#ifndef FILE_ADFRONT3
#define FILE_ADFRONT3

#include <meshing/meshtype.hh>

class frontpoint3
{
public:

  POINT3D p;
  INDEX globalindex;
  int nlinetopoint;
  int frontnr;

  frontpoint3 () { globalindex = nlinetopoint = 0; frontnr = 1000; }
  frontpoint3 (const POINT3D & ap, INDEX agi)
  { p = ap; globalindex = agi; nlinetopoint = 0; frontnr = 1000; }

  int Valid () const { return nlinetopoint >= 0; }
  void Invalidate () { nlinetopoint = -1; }
};

class frontface
{
public:

  ELEMENT f;
  int qualclass;
  char oldfront;

  frontface () { qualclass = 1; oldfront = 0; }
  frontface (const ELEMENT & af)
  { f = af; oldfront = 0; qualclass = 1; }

  int Valid () const { return f.PNum(1) != 0; }
  void Invalidate ()
  { f.PNum(1) = 0; oldfront = 0; qualclass = 1000; }
};




class ADFRONT3
{
  ARRAY<frontpoint3> points;
  ARRAY<frontface> faces;

  ARRAY<INDEX> delpointl;

  INDEX nff;  // number of front faces;

  float h;
  double vol;

public:

  ADFRONT3 ();
  //  void Load (char * filename, float & h);
  //  void Save (char * filename, float h);


  void GetPoints (ARRAY<POINT3D> & apoints) const;
  void Print () const;

  int Empty () const { return nff == 0; }

  int GetLocals (ARRAY<POINT3D> & locpoints,
                 ARRAY<ELEMENT> & locfaces,   // local index
                 ARRAY<INDEX> & pindex,
                 ARRAY<INDEX> & findex,
                 float xh);

  void GetGroup (int fi,
                 ARRAY<POINT3D> & grouppoints,
                 ARRAY<ELEMENT> & groupelements,
                 ARRAY<INDEX> & pindex,
                 ARRAY<INDEX> & findex
                 ) const;

  void DeleteFace (INDEX fi);
  INDEX AddPoint (const POINT3D & p, INDEX globind);
  INDEX AddFace (const ELEMENT & e);
  void IncrementClass (INDEX fi);
  void ResetClass (INDEX fi);
  int TestOk () const;
  void SetStartFront ();

  INDEX GetGlobalIndex (INDEX pi) const { return points[pi].globalindex; }
  float Volume () const { return float(vol); }

  void SaveSurface (char * filename, double h);
  void PrintSurface () const;
};

#endif
