// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>

#include <limits.h>

#include <template.hh>



#include "bitarray.hh"



BITARRAY :: BITARRAY ()

{

  size = 0;

  data = NULL;

}



BITARRAY :: BITARRAY (INDEX asize)

{

  size = 0;

  data = NULL;

  SetSize (asize);

}



BITARRAY :: ~BITARRAY ()

{

  if (data) delete data;

}



void BITARRAY :: SetSize (INDEX asize)

{

  if (size == asize) return;

  if (data) delete data;



  size = asize;

  data = new unsigned char [Addr (size)+1];

}



void BITARRAY :: Set ()

{

  INDEX i;

  for (i = 0; i <= Addr (size); i++)

    data[i] = UCHAR_MAX;

}



void BITARRAY :: Clear ()

{

  INDEX i;

  for (i = 0; i <= Addr (size); i++)

    data[i] = 0;

}
