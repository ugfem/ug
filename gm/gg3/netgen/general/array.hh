// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   array.hh                                                     */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/

#ifndef FILE_ARRAY
#define FILE_ARRAY





extern void MyError (char *);



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



  inline ARRAY(INDEX asize = 0, INDEX ainc = 0)

    : BASE_ARRAY (asize, ainc, sizeof (T)) { };



  inline ~ARRAY ()

  {

    if (data) delete (void*)data;

  }



  inline INDEX Size() const

  {

    return actsize;

  }



  inline void SetSize(INDEX nsize)

  {

    if (nsize > allocsize)

      ReSize (nsize, sizeof(T));

    actsize = nsize;

  }



  inline void SetAllocSize (INDEX nallocsize)

  {

    if (nallocsize > allocsize)

      ReSize (nallocsize, sizeof(T));

  }



  inline T & operator[] (INDEX i)

  {

    RangeOk (i);

    return ((T*)data)[i-1];

  }



  inline const T & operator[] (INDEX i) const

  {

    RangeOk (i);

    return ((const T*)data)[i-1];

  }



  inline T & Elem (INDEX i)

  {

    return ((T*)data)[i-1];

  }



  inline const T & Get (INDEX i) const

  {

    return ((const T*)data)[i-1];

  }



  inline void Set (INDEX i, const T & el)

  {

    ((T*)data)[i-1] = el;

  }



  inline T & Last ()

  {

    CheckNonEmpty ();

    return ((T*)data)[actsize-1];

  }



  inline const T & Last () const

  {

    CheckNonEmpty ();

    return ((const T*)data)[actsize-1];

  }



  inline INDEX Append (const T & el)

  {

    if (actsize == allocsize) ReSize (actsize+1, sizeof (T));

    ((T*)data)[actsize] = el;

    actsize++;

    return actsize;

  }



  inline void DeleteElement (INDEX i)

  {

    RangeOk (i);

    ((T*)data)[i-1] = ((T*)data)[actsize-1];

    actsize--;

  }



  inline void DeleteLast ()

  {

    CheckNonEmpty ();

    actsize--;

  }



  inline void DeleteAll ()

  {

    if (data) delete (void*)data;

    data = NULL;

    actsize = allocsize = 0;

  }



private:



  ARRAY<T> & operator= (ARRAY<T> &);

  ARRAY (const ARRAY<T> &);

};


#endif
