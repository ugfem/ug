// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_NETRULE
#define FILE_NETRULE

class netrule
{
  typedef struct tf { float f1, f2, f3; } threefloat;
  class threeint { public: int i1, i2, i3; threeint() { } };

  int quality;
  char * name;
  ARRAY<Point2d> points;
  ARRAY<ILINE> lines;
  ARRAY<Point2d> freezone;
  ARRAY<Point2d> transfreezone;

  ARRAY<int> dellines;
  ARRAY<Element> elements;
  ARRAY<threefloat> tolerances, linetolerances;
  ARRAY<threeint> orientations;
  //  SparseMatrix oldutonewu, oldutofreearea;
  //  SparseMatrix freesetinequ;

  ARRAY<Vec2d> linevecs;

  int noldp, noldl;
  float fzminx, fzmaxx, fzminy, fzmaxy;

public:

  netrule ();

  int GetNP () const { return points.Size(); }
  int GetNL () const { return lines.Size(); }
  int GetNE () const { return elements.Size(); }
  int GetNOldP () const { return noldp; }
  int GetNOldL () const { return noldl; }
  int GetNDelL () const { return dellines.Size(); }
  int GetNOrientations () const { return orientations.Size(); }
  int GetQuality () const { return quality; }

  const Point2d & GetPoint (int i) const { return points.Get(i); }
  const ILINE & GetLine (int i) const { return lines.Get(i); }
  const Element & GetElement (int i) const { return elements.Get(i); }
  const threeint & GetOrientation (int i) const { return orientations.Get(i); }
  int GetDelLine (int i) const { return dellines.Get(i); }

  void GetFreeZone (ARRAY<Point2d> & afreearea);
  float CalcPointDist (int pi, const Point2d & p) const;
  float CalcLineError (int li, const Vec2d & v) const;

  void SetFreeZoneTransformation (const Vector & u);
  int IsInFreeZone (const Point2d & p) const;
  int IsLineInFreeZone (const Point2d & p1, const Point2d & p2) const;
  int ConvexFreeZone () const;
  const ARRAY<Point2d> & GetTransFreeZone () { return transfreezone; }

  int GetPointNr (int ln, int endp) const { return lines.Get(ln).I(endp); }

  //  const SparseMatrix & GetOldUToNewU () const { return oldutonewu; }
  //  const SparseMatrix & GetOldUToFreeArea () const { return oldutofreearea; }
  const char * Name () const { return name; }

  void LoadRule (istream & ist);
};


extern int ApplyRules (const ARRAY<netrule*> & rules,
                       ARRAY<Point2d> & lpoints, ARRAY<ILINE> & llines,
                       ARRAY<Element> & elements, ARRAY<INDEX> & dellines,
                       int tolerance);

extern void DrawRules ();
#endif
