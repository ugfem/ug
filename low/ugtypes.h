// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* $Id$ */

/*

   globally set data types

 */

#ifndef UGTYPES_H
#define UGTYPES_H

/* standard types */

/* !!! maybe insert values via ugtypes.h.in? */

typedef short SHORT;

/* these types are used for several bitfields. I'd guess that it needs
   at least 32 bits... */
typedef int INT;
typedef unsigned int UINT;

typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif
