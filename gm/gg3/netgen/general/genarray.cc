// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**************************************************************************/
/* File:   array.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   Abstract data type ARRAY
 */


#include <stdlib.h>
#include <string.h>

#include <template.hh>
#include <array.hh>


BASE_ARRAY :: BASE_ARRAY(INDEX asize, INDEX ainc, int elementsize)
{
  void * p;

  if (asize)
  {
    p = new char[asize * elementsize];

    if (p)
    {
      data = p;
      actsize = allocsize = asize; inc = ainc;
      return;
    }
    else
    {
      MyError ("Array not allocated");
    }
  }

  data = NULL;
  actsize = allocsize = 0;
  inc = ainc;
}



void BASE_ARRAY :: ReSize (INDEX minsize, int elementsize)
{
  void * p;
  INDEX nsize = (inc) ? allocsize + inc : 2 * allocsize;
  if (nsize < minsize) nsize = minsize;

  /*
     if (long(nsize) * long(sizeof(T)) > 32767)
      MyError ("Array too large");
   */

  if (data)
  {
    p = new char [nsize * elementsize];
    memcpy (p, data, (min (nsize, actsize)) * elementsize);

    delete data;
    data = p;
  }
  else
  {
    p = new char[nsize * elementsize];
    data = p;
  }

  if (!data) MyError ("Array not allocated");

  allocsize = nsize;
}



int BASE_ARRAY :: RangeOk (INDEX i) const
{
  if (i < 1 || i > actsize)
  {
    MyError ("Array out of Range");
    return 0;
  }
  return 1;
}

int BASE_ARRAY :: CheckNonEmpty () const
{
  if (!actsize)
  {
    MyError ("Array souldn't be empty");
    return 0;
  }
  return 1;
}
