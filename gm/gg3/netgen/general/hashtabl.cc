// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/**************************************************************************/
/* File:   hashtabl.cc                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   Abstract data type HASHTABLE
 */

#include <stdlib.h>
#include <string.h>

#include <template.hh>
#include <array.hh>
#include <table.hh>
#include <hashtabl.hh>


void INDEX_2 :: Sort ()
{
  if (i1 > i2) swap (i1, i2);
}

void INDEX_3 :: Sort ()
{
  if (i1 > i2) swap (i1, i2);
  if (i2 > i3) swap (i2, i3);
  if (i1 > i2) swap (i1, i2);
}

int BASE_INDEX_2_HASHTABLE :: Position (int bnr, const INDEX_2 & ind) const
{
  int i;
  for (i = 1; i <= hash.EntrySize (bnr); i++)
    if (hash.Get(bnr, i) == ind)
      return i;
  return 0;
}

int BASE_INDEX_3_HASHTABLE :: Position (int bnr, const INDEX_3 & ind) const
{
  int i;
  for (i = 1; i <= hash.EntrySize (bnr); i++)
    if (hash.Get(bnr, i) == ind)
      return i;
  return 0;
}
