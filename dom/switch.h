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
/*			  internet: ug@ica3.uni-stuttgart.de                            */
/*																			*/
/* History:   29.10.93 begin, ug 3D-version                                                             */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifndef __SWITCH__
#define __SWITCH__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifdef __MWCW__
#include "MWCW.cmdlinedefs"
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
   __NODEDATA__
   __EDGEDATA__
   __SIDEDATA__
   __ELEMDATA__
   .ve

   If defined, these constants indicate that degrees of freedom are allowed
   in the respective object.

   .vb
   __TWODIM__
   __THREEDIM__
   .ve

   One of these two constants is defined, indicating compilation
   of a 2D or 3D version.

   .vb
   __version23__
   __version3__
   .ve

   Indicates which version is compiled. '__version23__' is for a version
   that is compatible with UG version 2.3 (a version without the
   matrix vector structure).

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

#ifdef NodeT
#define Nodeon
#else
#define Nodeoff
#endif

#ifdef EdgeT
#define Edgeon
#else
#define Edgeoff
#endif

#ifdef SideT
#define Sideon
#else
#define Sideoff
#endif

#ifdef ElemT
#define Elemon
#else
#define Elemoff
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

#if defined Nodeoff && defined Edgeoff && defined Sideoff && defined Elemoff
        #ifdef _3
        #error ****	define at least one of NodeT, ElemT, SideT, EdgeT for 3d version	   ****
        #endif
#define __version23__
#endif

#ifdef two
#ifdef Sideon
#error ****   two dimensional case cannot have side data	****
#endif
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

#ifdef Nodeon
#define __NODEDATA__
#define NODE_DATA       1
#endif
#ifndef Nodeon
#define NODE_DATA       0
#endif

#ifdef Edgeon
#define __EDGEDATA__
#define EDGE_DATA       2
#endif
#ifndef Edgeon
#define EDGE_DATA       0
#endif

#ifdef Elemon
#define __ELEMDATA__
#define ELEM_DATA       4
#endif
#ifndef Elemon
#define ELEM_DATA       0
#endif

#ifdef Sideon
#define __SIDEDATA__
#define SIDE_DATA       8
#endif
#ifndef Sideon
#define SIDE_DATA       0
#endif

#ifndef __version23__
#define __version3__
#endif

#define DATA_TYPES                      (NODE_DATA | EDGE_DATA | SIDE_DATA | ELEM_DATA)
#define ALL_TYPES               (NODE_DATA | EDGE_DATA | SIDE_DATA | ELEM_DATA)
#define TYPES_EXISTING(t)       ((((t)|DATA_TYPES)==DATA_TYPES)&&(t))

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
        #define __MIDNODE__
#endif

/* values for vectors */
#define NODEVECTOR                                      0                               /* vector associated to a node					*/
#define EDGEVECTOR                                      1                               /* vector associated to an edge                                 */
#define ELEMVECTOR                                      2                               /* vector associated to an element				*/
#define SIDEVECTOR                                      3                               /* vector associated to an elementside			*/

#define DATATYPE_FROM_VECTORTYPE(t) (1<<(t))    /* transforms VectorType into DataType			*/

#endif
