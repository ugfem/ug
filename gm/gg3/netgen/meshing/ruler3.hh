// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_RULER3
#define FILE_RULER3


class vnetrule
{
public:
  //  struct twoint { int i1, i2; };
  //  struct threeint { int i1, i2, i3; };
  //  struct fourint { int i1, i2, i3, i4; };

private:

  int quality;
  char * name;
  ARRAY<Point3d> points;
  ARRAY<Element> faces;
  ARRAY<twoint> edges;

  ARRAY<Point3d> freezone;
  ARRAY<ARRAY<threeint>*> freefaces;
  ARRAY<ARRAY<int>*> freesets;
  ARRAY<Point3d> transfreezone;

  ARRAY<int> delfaces;
  ARRAY<Element> elements;
  ARRAY<double> tolerances, linetolerances;
  SparseMatrix oldutonewu, oldutofreezone;
  ARRAY<DenseMatrix*> freefaceinequ;
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


  const Point3d & GetPoint (int i) const { return points.Get(i); }
  const Element & GetFace (int i) const { return faces.Get(i); }
  const Element & GetElement (int i) const { return elements.Get(i); }
  const twoint & GetEdge (int i) const { return edges.Get(i); }
  int GetDelFace (int i) const { return delfaces.Get(i); }
  int IsDelFace (int fn) const;

  float CalcPointDist (int pi, const Point3d & p) const;

  void SetFreeZoneTransformation (const Vector & u);
  int IsInFreeZone (const Point3d & p) const;
  int IsTriangleInFreeZone (const Point3d & p1, const Point3d & p2,
                            const Point3d & p3);
  int IsTriangleInFreeSet (const Point3d & p1, const Point3d & p2,
                           const Point3d & p3, int fs);
  int ConvexFreeZone () const;
  const Point3d & GetTransFreeZone (int i) { return transfreezone.Get(i); }

  int GetNP (int fn) const
  { return faces.Get(fn).NP(); }
  int GetPointNr (int fn, int endp) const
  { return faces.Get(fn).PNum(endp); }
  int GetPointNrMod (int fn, int endp) const
  { return faces.Get(fn).PNumMod(endp); }
  const fourint & GetOrientation (int i) { return orientations.Get(i); }

  int TestFlag (char flag) const;

  const SparseMatrix & GetOldUToNewU () const { return oldutonewu; }
  const SparseMatrix & GetOldUToFreeZone () const { return oldutofreezone; }
  const char * Name () const { return name; }
  void LoadRule (istream & ist);

  const ARRAY<Point3d> & GetTransFreeZone () { return transfreezone; }
  int TestOk () const;

  friend void TestRules ();
  friend void Plot3DRule (const ROT3D & r, char key);
};



extern int ApplyVRules (const ARRAY<vnetrule*> & rules, double tolfak,
                        ARRAY<Point3d> & lpoints, ARRAY<Element> & lfaces,
                        ARRAY<Element> & elements,
                        ARRAY<INDEX> & delfaces, int tolerance, int rotind1,
                        float & retminerr, ARRAY<char*> & problems);

#endif
