// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  switch.h														*/
/*																			*/
/* Purpose:   standard header file template                                                             */
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de                                        */
/*																			*/
/* History:   29.10.93 begin, ug 3D-version                                                             */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifndef __SWITCH__
#define __SWITCH__


#ifndef __GENERAL__
#include "general.h"
#endif

#ifndef __COMPILER__
#include "compiler.h"
#endif

/*D
   switch - switch header file, extracts makefile switches and sets some constants

   DESCRIPTION:

   This file processes the switches defined in the 'makefile' (see page 'makefiles')
   and sets some constants depending on these switches. This file is included
   in all UG source files.

   The following constants defined with the '#define' preprocessor command
   can be used for conditional compilation.

   .vb
   __TWODIM__
   __THREEDIM__
   .ve

   One of these two constants is defined, indicating compilation
   of a 2D or 3D version.

   .vb
   DIM
   DIM_OF_BND
   CORNERS_OF_BND_SEG
   MAXVECTORS
   .ve

   These constants are dimension-dependent and are set to the space
   dimension, the dimension of the boundary, the number of corners
   of a boundary segment and the possible number of different 'VECTOR'
   sizes.

   .vb
   NODEVECTOR
   EDGEVECTOR
   ELEMVECTOR
   SIDEVECTOR
   .ve

   These constants encode the geometric position of a 'VECTOR'.

   D*/


/****************************************************************************/
/*																			*/
/* check for defines from makefile											*/
/*																			*/
/****************************************************************************/

#if (!defined _2) && (!defined _3)
#error ****	define dimension _2 or _3		****
#endif

#ifdef _2
        #ifdef _3
        #error ****	define EITHER dimension _2 OR _3	   ****
        #endif
#define two
#endif

#ifdef _3
#define three
#endif

#ifdef two
#ifdef three
#error ****	define at most dimension two OR three		****
#endif
#endif

#ifndef two
#ifndef three
#error ****	define at least dimension two OR three		****
#endif
#endif

#ifdef two
#ifdef Sideon
#error ****   two dimensional case cannot have side data	****
#endif
#endif

#ifdef ModelP
#define MODEL "PARALLEL"
#else
#define MODEL "SEQUENTIAL"
#endif

#undef __GRAPE_TRUE__
#ifdef _3
#ifdef GRAPET
#define __GRAPE_TRUE__
#endif
#endif

#ifdef GRAPEF
#define GRAPE "OFF"
#endif
#ifdef GRAPET
#define GRAPE "ON"
#endif

#ifdef NETGENF
#define NETGEN "OFF"
#endif
#ifdef NETGENT
#define NETGEN "ON"
#endif

#ifdef Debug
#define DEBUG_MODE "ON"
#else
#define DEBUG_MODE "OFF"
#endif

/****************************************************************************/
/*																			*/
/* basic switch-defines for:												*/
/*	  dimensions															*/
/*	  user data in geometric components (nodes,edges,sides,elems)			*/
/*																			*/
/*	  dependent of makefile-defines                                                                                 */
/*																			*/
/****************************************************************************/

#ifdef two
#define __TWODIM__
#endif

#ifdef three
#define __THREEDIM__
#endif

/****************************************************************************/
/*																			*/
/* switch-define dependent defines											*/
/*																			*/
/****************************************************************************/

#define DIM_MAX                                         3                               /* maximal space dimension						*/
#define DIM_OF_BND_MAX                          2                               /* maximal dimension of boundary surface		*/

#ifdef __TWODIM__
        #define DIM                                     2                               /* space dimension								*/
        #define DIM_OF_BND                              1                               /* dimension of boundary surface				*/
        #define CORNERS_OF_BND_SEG              2                               /* number of corners of a boundary side                 */
        #define MAXVECTORS                              3                               /* three different data types					*/
#endif

#ifdef __THREEDIM__
        #define DIM                                     3                               /* space dimension								*/
        #define DIM_OF_BND                              2                               /* dimension of boundary surface				*/
        #define CORNERS_OF_BND_SEG              4                               /* number of corners of a boundary side                 */
        #define MAXVECTORS                              4                               /* four different data types					*/
#endif

/* values for vectors */
#define NODEVECTOR                                      0                               /* vector associated to a node					*/
#define EDGEVECTOR                                      1                               /* vector associated to an edge                                 */
#define ELEMVECTOR                                      2                               /* vector associated to an element				*/
#define SIDEVECTOR                                      3                               /* vector associated to an elementside			*/

#define DATATYPE_FROM_VECTORTYPE(t) (1<<(t))    /* transforms VectorType into DataType			*/

#endif
