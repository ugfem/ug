// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      compiler.h                                                    */
/*                                                                          */
/* Purpose:   define simple data types and standard include files           */
/*                                                                          */
/* Author:      Peter Bastian                                               */
/*              Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen */
/*              Universitaet Heidelberg                                     */
/*              Im Neuenheimer Feld 368                                     */
/*              6900 Heidelberg                                             */
/*                                                                          */
/* History:   29.01.92 begin, ug version 2.0                                */
/*                                                                          */
/* Revision:  04.09.95                                                      */
/*                                                                          */
/****************************************************************************/

#ifndef __COMPILER__
#define __COMPILER__

#include <limits.h>
#include <float.h>

#define __MWCW__  /* this is the default */

/****************************************************************************/
/*                                                                              */
/* #define exactly one of the following constants:                          */
/*                                                                          */
/*            __MPW32__     Apple MacIntosh Programmers Workshop version 3.2 */
/*            __SUN4GCC__  SUN Workstation SunOS version >= 4.0, gnu c comp! */
/*            __IRIS__     Silicon Graphics Workstations                    */
/*	                                                                        */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Definitions for Apple MacIntosh                                          */
/*                                                                          */
/****************************************************************************/

/* SMALL..: least number s.t. 1 + SMALL../SMALL_FAC != 1 */
#define SMALL_FAC            10

#ifdef __MPW32__
#undef __MWCW__

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#include <Memory.h>
#define malloc(n) ((void *) NewPtr((Size) n))
#define free(p) DisposPtr((Ptr) p)
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment */

/* Diese bloeden const pointer gehen immer noch nicht ! P.B. 27.7.95 */
#define const

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for Sun Version                                              */
/*                                                                          */
/****************************************************************************/

#ifdef __SUN4GCC__
#undef __MWCW__
#include <stddef.h>

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment */

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for Silicon Graphics Workstations (old 32 bit version)       */
/*                                                                          */
/****************************************************************************/

#ifdef __SGI__
#undef __MWCW__

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment */

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for HP Workstations                                          */
/*                                                                          */
/****************************************************************************/

#ifdef __GCC__
#undef __MWCW__

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD float

/* memory */
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment */

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for HP Workstations                                          */
/*                                                                          */
/****************************************************************************/

#ifdef __HP__
#undef __MWCW__

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment */

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for Apple Power Macintosh                                    */
/*                                                                          */
/****************************************************************************/

#ifdef __MWCW__

#define USES_UNIVERSAL_HEADERS

/*
    for older versions of the CW C-compiler you may have to
   #define __MWCW_oldVersion__
    (see MacShell.c)
 */

/* basic types */
#define SHORT            short
#define INT             int
#define FLOAT            float
#define DOUBLE            double
#define COORD            float
#define SCREEN_COORD    float

/* memory */
#include <Memory.h>
#define malloc(n) ((void *) NewPtr((Size) n))
#define free(p) DisposPtr((Ptr) p)
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment */

#endif

/****************************************************************************/
/*                                                                          */
/* some general definitions                                                 */
/*                                                                          */
/****************************************************************************/

/* limits of the basic types */
#define MAX_S            SHRT_MAX
#define MAX_I            INT_MAX
#define MAX_F            FLT_MAX
#define SMALL_F         (FLT_EPSILON*SMALL_FAC)
#define MAX_D            DBL_MAX
#define SMALL_D         (DBL_EPSILON*SMALL_FAC)
#define MAX_C            FLT_MAX
#define SMALL_C         (FLT_EPSILON*SMALL_FAC)

#endif
