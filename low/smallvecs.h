// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* $Id$ */

/*

   small vector types

 */

#ifndef UG_SMALLVECS_H
#define UG_SMALLVECS_H

/* convert _2 and _3 to DIM and others */
#include "dimension.h"

/* basic types */
#include "ugtypes.h"

typedef DOUBLE DOUBLE_VECTOR[DIM];
typedef DOUBLE DOUBLE_VECTOR_2D[2];
typedef DOUBLE DOUBLE_VECTOR_3D[3];


#endif
