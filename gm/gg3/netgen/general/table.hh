// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   table.hh                                                     */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#ifndef FILE_TABLE
#define FILE_TABLE


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

  void IncSize (int i, int elsize);

};





template <class T>

class TABLE : public BASE_TABLE

{

public:



  inline TABLE (int size) : BASE_TABLE (size)

  {

    ;

  }



  inline void Add (INDEX i, const T & acont)

  {

    IncSize (i, sizeof (T));

    ((T*)data.Elem(i).col)[data.Elem(i).size-1] = acont;

  }



  inline void Set (INDEX i, int nr, const T & acont)

  {

    ((T*)data.Get(i).col)[nr-1] = acont;

  }



  inline const T & Get (INDEX i, int nr) const

  {

    return ((T*)data.Get(i).col)[nr-1];

  }



  inline int Size () const

  {

    return data.Size();

  }



  inline int EntrySize (int i) const

  {

    return data.Get(i).size;

  }

};



#endif
