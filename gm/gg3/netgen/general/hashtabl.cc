// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   tashtabl.cc                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/

#include <stdlib.h>

#include <string.h>



#include <template.hh>

#include <array.hh>

#include <table.hh>

#include <hashtabl.hh>





int BASE_INDEX_2_HASHTABLE :: Position (int bnr, const INDEX_2 & ind) const

{

  int i;

  for (i = 1; i <= hash.EntrySize (bnr); i++)

    if (hash.Get(bnr, i) == ind)

      return i;

  return 0;

}
