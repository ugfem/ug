// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**************************************************************************/
/* File:   bitarray.cc                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   data type BitArray
 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <template.hh>

#include "bitarray.hh"

BitArray :: BitArray ()
{
  size = 0;
  data = NULL;
}

BitArray :: BitArray (INDEX asize)
{
  size = 0;
  data = NULL;
  SetSize (asize);
}

BitArray :: ~BitArray ()
{
  if (data) delete data;
}

void BitArray :: SetSize (INDEX asize)
{
  if (size == asize) return;
  if (data) delete data;

  size = asize;
  data = new unsigned char [Addr (size)+1];
}

void BitArray :: Set ()
{
  INDEX i;
  if (!size) return;
  for (i = 0; i <= Addr (size); i++)
    data[i] = UCHAR_MAX;
}

void BitArray :: Clear ()
{
  INDEX i;
  if (!size) return;
  for (i = 0; i <= Addr (size); i++)
    data[i] = 0;
}
