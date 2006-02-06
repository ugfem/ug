// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* $Id$ */

/** \file
    \brief Provide the usual preprocessor-defines for the dimension and complain
    if the dimension was set incorrectly.
 */


#ifndef DIMENSION_H
#define DIMENSION_H


#ifdef _2
#ifdef _3
#error ****    define EITHER dimension _2 OR _3       ****
#endif
#define __TWODIM__
#define DIM 2
#define DIM_OF_BND 1
#endif

#ifdef _3
#define __THREEDIM__
#define DIM 3
#define DIM_OF_BND 2
#endif

#ifndef _2
#ifndef _3
#error ****    define at least dimension two OR three        ****
#endif
#endif

#endif
