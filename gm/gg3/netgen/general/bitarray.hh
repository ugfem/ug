// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_BitArray
#define FILE_BitArray

/**************************************************************************/
/* File:   bitarray.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   data type BitArray

   BitArray is a compressed array of Boolean information. By Set and Clear
   the whole array or one bit can be set or reset, respectively.
   Test returns the state of the accoring bit.
   No range checking is done.
 */

#include <limits.h>

class BitArray
{
  INDEX size;
  unsigned char * data;

public:
  BitArray ();
  BitArray (INDEX asize);
  ~BitArray ();

  void SetSize (INDEX asize);
  inline INDEX Size () const;

  void Set ();
  inline void Set (INDEX i);
  void Clear ();
  inline void Clear (INDEX i);
  inline BOOL Test (INDEX i) const;

private:
  inline unsigned char Mask (INDEX i) const;
  inline INDEX Addr (INDEX i) const;
};


#include <bitarray.icc>

#endif
