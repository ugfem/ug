// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_ADFRONT2
#define FILE_ADFRONT2

#include <meshing/meshtype.hh>

/*

    Advancing front class for surfaces

 */

class ADFRONT2
{

  class frontpoint2
  {
    Point3d p;             // coordinates
    INDEX globalindex;     // global node index
    int nlinetopoint;      // number of front lines connected to point
    int frontnr;           // distance to original boundary

  public:

    frontpoint2 ();
    frontpoint2 (const Point3d & ap, INDEX agi);

    const Point3d & P () const;
    INDEX GlobalIndex () const;

    void AddLine ();
    void RemoveLine ();
    int Valid () const;

    void DecFrontNr (int afrontnr);
    int FrontNr () const;
  };


  class frontline
  {
    ILINE l;             // Point Indizes
    int lineclass;       // quality class
    int surfaceindex;    // index of geometric surface

  public:

    frontline ();
    frontline (const ILINE & al, int asi);

    const ILINE & L () const;
    int SurfaceIndex () const;

    int LineClass () const;

    void IncrementClass ();
    void ResetClass ();

    int Valid () const;
    void Invalidate ();
  };



  ARRAY<frontpoint2> points;
  ARRAY<frontline> lines;

  ARRAY<INDEX> delpointl;
  ARRAY<INDEX> dellinel;

  INDEX nfl;  // number of front lines;


public:

  ADFRONT2 ();

  void GetPoints (ARRAY<Point3d> & apoints) const;
  void Print (ostream & ost) const;

  int Empty () const;

  int GetLocals (ARRAY<Point3d> & locpoints,
                 ARRAY<ILINE> & loclines,   // local index
                 ARRAY<INDEX> & pindex,
                 ARRAY<INDEX> & lindex,
                 int & surfind,
                 double xh);

  void DeleteLine (INDEX li);
  INDEX AddPoint (const Point3d & p, INDEX globind);
  INDEX AddLine (INDEX pi1, INDEX pi2, int asurfaceind);
  void IncrementClass (INDEX li);
  void ResetClass (INDEX li);

  INDEX GetGlobalIndex (INDEX pi) const;
  void SetStartFront ();
};








// inline functions

inline int ADFRONT2::frontpoint2 :: Valid () const
{
  return nlinetopoint >= 0;
}

inline const Point3d & ADFRONT2::frontpoint2 :: P () const
{
  return p;
}

inline INDEX ADFRONT2::frontpoint2 :: GlobalIndex () const
{
  return globalindex;
}

inline void ADFRONT2::frontpoint2 :: AddLine ()
{
  nlinetopoint++;
}

inline void ADFRONT2::frontpoint2 :: RemoveLine ()
{
  nlinetopoint--;
  if (nlinetopoint == 0)
    nlinetopoint = -1;
}

inline int ADFRONT2::frontpoint2 :: FrontNr () const
{
  return frontnr;
}

inline void ADFRONT2::frontpoint2 :: DecFrontNr (int afrontnr)
{
  if (frontnr > afrontnr) frontnr = afrontnr;
}










inline int ADFRONT2::frontline :: Valid () const
{
  return l.I1() != 0;
}

inline void ADFRONT2::frontline :: Invalidate ()
{
  l.I1() = 0;
  l.I2() = 0;
  lineclass = 1000;
}

inline const ILINE & ADFRONT2::frontline :: L () const
{
  return l;
}

inline int ADFRONT2::frontline :: SurfaceIndex () const
{
  return surfaceindex;
}


inline int ADFRONT2::frontline :: LineClass () const
{
  return lineclass;
}

inline void ADFRONT2::frontline :: IncrementClass ()
{
  lineclass++;
}

inline void ADFRONT2::frontline :: ResetClass ()
{
  lineclass = 1;
}


inline int ADFRONT2 :: Empty () const
{
  return nfl == 0;
}

/*inline INDEX ADFRONT2 :: GetGlobalIndex (INDEX pi) const
   {
   return points[pi].GlobalIndex();
   }*/

#endif
