// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_ADFRONT3
#define FILE_ADFRONT3

/**************************************************************************/
/* File:   adfront3.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
    Advancing front class for volume meshing
 */


#include <meshing/meshtype.hh>



class ADFRONT3
{

  class FrontPoint3
  {
    Point3d p;            // coordinates
    INDEX globalindex;    // global node index
    int nfacetopoint;     // number of faces connected to point
    int frontnr;          // distance to original boundary

  public:

    FrontPoint3 ();
    FrontPoint3 (const Point3d & ap, INDEX agi);

    const Point3d & P () const;
    INDEX GlobalIndex () const;

    void AddFace ();
    void RemoveFace();

    int Valid () const;
    //    void Invalidate ();

    void DecFrontNr (int afrontnr);
    int FrontNr () const;
  };

  class FrontFace
  {
    Element f;
    int qualclass;
    char oldfront;

  public:
    FrontFace ();
    FrontFace (const Element & af);
    const Element & Face () const;

    int QualClass () const;
    void IncrementQualClass ();
    void ResetQualClass ();

    int Valid () const;
    void Invalidate ();
  };




  ARRAY<FrontPoint3> points;
  ARRAY<FrontFace> faces;

  ARRAY<INDEX> delpointl;

  INDEX nff;  // number of front faces;

  double h;
  double vol;

public:

  ADFRONT3 ();
  //  void Load (char * filename, float & h);
  //  void Save (char * filename, float h);


  void GetPoints (ARRAY<Point3d> & apoints) const;
  void Print () const;

  int Empty () const;

  int GetLocals (ARRAY<Point3d> & locpoints,
                 ARRAY<Element> & locfaces,   // local index
                 ARRAY<INDEX> & pindex,
                 ARRAY<INDEX> & findex,
                 float xh);

  void GetGroup (int fi,
                 ARRAY<Point3d> & grouppoints,
                 ARRAY<Element> & groupelements,
                 ARRAY<INDEX> & pindex,
                 ARRAY<INDEX> & findex
                 ) const;

  void DeleteFace (INDEX fi);
  INDEX AddPoint (const Point3d & p, INDEX globind);
  INDEX AddFace (const Element & e);
  void IncrementClass (INDEX fi);
  void ResetClass (INDEX fi);
  int TestOk () const;
  void SetStartFront ();

  INDEX GetGlobalIndex (INDEX pi) const;
  double Volume () const;


  void SaveSurface (char * filename, double h);
  void PrintSurface () const;
};



inline const Point3d & ADFRONT3 :: FrontPoint3 :: P () const
{
  return p;
}

inline INDEX ADFRONT3 :: FrontPoint3 :: GlobalIndex () const
{
  return globalindex;
}

inline void ADFRONT3::FrontPoint3 :: AddFace ()
{
  nfacetopoint++;
}

inline void ADFRONT3::FrontPoint3 :: RemoveFace ()
{
  nfacetopoint--;
  if (nfacetopoint == 0)
    nfacetopoint = -1;
}



inline int ADFRONT3 :: FrontPoint3 :: Valid () const
{
  return nfacetopoint >= 0;
}

/*
   inline void ADFRONT3 :: FrontPoint3 :: Invalidate ()
   {
   nfacetopoint = -1;
   }
 */

inline int ADFRONT3 :: FrontPoint3 :: FrontNr () const
{
  return frontnr;
}

inline void ADFRONT3 :: FrontPoint3 :: DecFrontNr (int afrontnr)
{
  if (frontnr > afrontnr) frontnr = afrontnr;
}






inline int ADFRONT3 :: FrontFace :: Valid () const
{
  return f.PNum(1) != 0;
}

inline const Element & ADFRONT3 :: FrontFace :: Face() const
{
  return f;
}

inline int ADFRONT3 :: FrontFace :: QualClass () const
{
  return qualclass;
}

inline void ADFRONT3 :: FrontFace :: IncrementQualClass ()
{
  qualclass++;
}

inline void ADFRONT3 :: FrontFace :: ResetQualClass ()
{
  if (qualclass > 1)
  {
    qualclass = 1;
    oldfront = 0;
  }
}

inline int ADFRONT3 :: Empty () const
{
  return nff == 0;
}

inline INDEX ADFRONT3 :: GetGlobalIndex (INDEX pi) const
{
  return points[pi].GlobalIndex();
}

inline double ADFRONT3 :: Volume () const
{
  return vol;
}


#endif
