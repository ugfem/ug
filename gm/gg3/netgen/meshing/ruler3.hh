// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   ruler3.hh                                                    */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/

#ifndef FILE_RULER3
#define FILE_RULER3


class vnetrule
{
private:

  int quality;
  char * name;
  ARRAY<POINT3D> points;
  ARRAY<ELEMENT> faces;
  ARRAY<twoint> edges;

  ARRAY<POINT3D> freezone;
  ARRAY<ARRAY<threeint>*> freefaces;
  ARRAY<ARRAY<int>*> freesets;
  ARRAY<POINT3D> transfreezone;

  ARRAY<int> delfaces;
  ARRAY<ELEMENT> elements;
  ARRAY<double> tolerances, linetolerances;
  SPARSE_MATRIX oldutonewu, oldutofreezone;
  ARRAY<MATRIX*> freefaceinequ;
  ARRAY<fourint> orientations;
  ARRAY<char> flags;
  ARRAY<int> fnearness;

  int noldp, noldf;
  float fzminx, fzmaxx, fzminy, fzmaxy, fzminz, fzmaxz;


public:

  vnetrule ();

  int GetNP () const { return points.Size(); }
  int GetNF () const { return faces.Size(); }
  int GetNE () const { return elements.Size(); }
  int GetNO () const { return orientations.Size(); }
  int GetNEd () const { return edges.Size(); }
  int GetNOldP () const { return noldp; }
  int GetNOldF () const { return noldf; }
  int GetNDelF () const { return delfaces.Size(); }
  int GetQuality () const { return quality; }
  int GetFNearness (int fi) const { return fnearness.Get(fi); }


  const POINT3D & GetPoint (int i) const { return points.Get(i); }
  const ELEMENT & GetFace (int i) const { return faces.Get(i); }
  const ELEMENT & GetElement (int i) const { return elements.Get(i); }
  const twoint & GetEdge (int i) const { return edges.Get(i); }
  int GetDelFace (int i) const { return delfaces.Get(i); }
  int IsDelFace (int fn) const;

  float CalcPointDist (int pi, const POINT3D & p) const;

  void SetFreeZoneTransformation (const VECTOR & u);
  int IsInFreeZone (const POINT3D & p) const;
  int IsTriangleInFreeZone (const POINT3D & p1, const POINT3D & p2,
                            const POINT3D & p3);
  int IsTriangleInFreeSet (const POINT3D & p1, const POINT3D & p2,
                           const POINT3D & p3, int fs);
  int ConvexFreeZone () const;
  const POINT3D & GetTransFreeZone (int i) { return transfreezone.Get(i); }

  int GetNP (int fn) const
  { return faces.Get(fn).NP(); }
  int GetPointNr (int fn, int endp) const
  { return faces.Get(fn).PNum(endp); }
  int GetPointNrMod (int fn, int endp) const
  { return faces.Get(fn).PNumMod(endp); }
  const fourint & GetOrientation (int i) { return orientations.Get(i); }

  int TestFlag (char flag) const;

  const SPARSE_MATRIX & GetOldUToNewU () const { return oldutonewu; }
  const SPARSE_MATRIX & GetOldUToFreeZone () const { return oldutofreezone; }
  const char * Name () const { return name; }
  void LoadRule (istream & ist);

  const ARRAY<POINT3D> & GetTransFreeZone () { return transfreezone; }
  int TestOk () const;

  friend void TestRules ();
  friend void Plot3DRule (const ROT3D & r, char key);
};



extern int ApplyVRules (const ARRAY<vnetrule*> & rules, double tolfak,
                        ARRAY<POINT3D> & lpoints, ARRAY<ELEMENT> & lfaces,
                        ARRAY<ELEMENT> & elements,
                        ARRAY<INDEX> & delfaces, int tolerance, int rotind1,
                        float & retminerr, ARRAY<char*> & problems);

#endif
