// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_ARRAY
#define FILE_ARRAY

/**************************************************************************/
/* File:   array.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   Abstract data type ARRAY

   ARRAY<T> is an automatically increasing array containing elements of the
   generic type T. The allocated size may be larger then the logical size
   of the array.
   The elements may be accessed by the brackets [ ] in the range from
   1 to size. In opposite to this possibility the elements can be accessed
   by the Methods Set, Get and Elem without range checking.
 */


class BASE_ARRAY
{
protected:
  void * data;
  INDEX actsize, allocsize, inc;

  BASE_ARRAY(INDEX asize, INDEX ainc, int elmentsize);
  void ReSize (INDEX minsize, int elementsize);
  int RangeOk (INDEX i) const;
  int CheckNonEmpty () const;
};



template <class T>
class ARRAY : private BASE_ARRAY
{
public:
  inline ARRAY(INDEX asize = 0, INDEX ainc = 0);
  inline ~ARRAY ();

  inline INDEX Size() const;
  inline void SetSize(INDEX nsize);
  inline void SetAllocSize (INDEX nallocsize);

  inline T & operator[] (INDEX i);
  inline const T & operator[] (INDEX i) const;
  inline T & Elem (INDEX i);
  inline const T & Get (INDEX i) const;
  inline void Set (INDEX i, const T & el);
  inline T & Last ();
  inline const T & Last () const;

  inline INDEX Append (const T & el);
  inline void DeleteElement (INDEX i);
  inline void DeleteLast ();
  inline void DeleteAll ();

private:
  ARRAY<T> & operator= (ARRAY<T> &);
  ARRAY (const ARRAY<T> &);
};


#include <array.icc>

#endif
