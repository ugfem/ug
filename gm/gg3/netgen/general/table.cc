// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**************************************************************************/
/* File:   table.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   Abstract data type TABLE
 */

#include <stdlib.h>
#include <string.h>

#include <template.hh>
#include <array.hh>
#include <table.hh>


BASE_TABLE :: BASE_TABLE (int size)
  : data(size)
{
  int i;
  for (i = 1; i <= size; i++)
  {
    data.Elem(i).maxsize = 0;
    data.Elem(i).size = 0;
    data.Elem(i).col = NULL;
  }
}

BASE_TABLE :: ~BASE_TABLE ()
{
  int i;
  for (i = 1; i <= data.Size(); i++)
    if (data.Get(i).col)
      delete data[i].col;
}

void BASE_TABLE :: SetSize (int size)
{
  int i;
  data.SetSize(size);
  for (i = 1; i <= size; i++)
  {
    data.Elem(i).maxsize = 0;
    data.Elem(i).size = 0;
    data.Elem(i).col = NULL;
  }
}

void BASE_TABLE :: IncSize (int i, int elsize)
{
  if (i < 1 || i > data.Size())
  {
    MyError ("BASE_TABLE::Inc: Out of range");
    return;
  }

  linestruct & line = data.Elem (i);

  if (line.size == line.maxsize)
  {
    void * p = new char [(line.maxsize+5) * elsize];

    if (!p)
    {
      MyError ("BASE_TABLE::Inc: Out of memory");
      return;
    }

    memcpy (p, line.col, line.maxsize * elsize);
    delete line.col;
    line.col = p;
    line.maxsize += 5;
  }

  line.size++;
}
