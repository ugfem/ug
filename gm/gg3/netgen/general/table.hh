// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_TABLE
#define FILE_TABLE

/**************************************************************************/
/* File:   table.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   Abstract data type TABLE

   To an integer i in the range from 1 to size a set of elements of the
   generic type T is associated.
 */


class BASE_TABLE
{
protected:

  class linestruct
  {
  public:
    int size, maxsize;
    void * col;
  };

  ARRAY<linestruct> data;

public:
  BASE_TABLE (int size);
  ~BASE_TABLE ();
  void SetSize (int size);
  void IncSize (int i, int elsize);
};


template <class T>
class TABLE : public BASE_TABLE
{
public:

  inline TABLE ();
  inline TABLE (int size);
  inline void SetSize (int size);
  inline void Add (INDEX i, const T & acont);
  inline void Set (INDEX i, int nr, const T & acont);
  inline const T & Get (INDEX i, int nr) const;
  inline int Size () const;
  inline int EntrySize (int i) const;
};

#include <table.icc>

#endif
