// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
    \brief Globally set data types

    The basic types (normally 'short', 'int', 'float' and 'double') are
    replaced by 'SHORT', 'INT', 'FLOAT' and 'DOUBLE'. The type 'DOUBLE'
    is used for all Cartesian coordinates of the (x,y[,z])-directions of
    the grids and 'SCREEN_DOUBLE' is used for all transformed coordinates
    of the graphical interface.
 */

#ifndef UGTYPES_H
#define UGTYPES_H

#include "namespace.h"

START_UG_NAMESPACE

#if defined __DEC__ || defined __LINUXAXP__ || defined __LINUXIA64__ || defined __AMD64__ || defined __PPC64__

typedef short SHORT;
/* With the old build system we want INT to have the same size as void* */
typedef long INT;
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#elif defined __NECSX4__ || defined __SR2201__ || defined __SR8K__

typedef short SHORT;
/* With the old build system we want INT to have the same size as void* */
typedef long INT;
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef double COORD;
typedef float SCREEN_COORD;

#else

/* standard types */

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

#ifdef __cplusplus

typedef bool BOOL;
#define TRUE true
#define FALSE false

#else

typedef _Bool BOOL;
#define TRUE 1
#define FALSE 0

#endif

END_UG_NAMESPACE

#endif
