// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      compiler.h                                                    */
/*                                                                          */
/* Purpose:   define simple data types and standard include files           */
/*                                                                          */
/* Author:    Peter Bastian, Klaus Birken, Stefan Lang                      */
/*            Institut fuer Computeranwendungen 3 ( ICA3 )                  */
/*            Pfaffenwaldring 27                                            */
/*            Universitaet Stuttgart                                        */
/*            70569 Stuttgart                                               */
/*                                                                          */
/* History:   29.01.92 begin, ug version 2.0                                */
/*            19.08.92 begin PARIX version                                  */
/*            02.06.93 begin PARAGON version                                */
/*            02.12.94 begin CRAY_T3D version                               */
/*            05.10.95 added GC (Explorer PowerPC) version                  */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */


#ifndef __COMPILER__
#define __COMPILER__

#include <limits.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif

#define __MWCW__  /* this is the default */

/****************************************************************************/
/*                                                                          */
/* #define exactly one of the following constants: (in Makefile)            */
/*                                                                          */
/*          __MPW32__    Apple MacIntosh Programmers Workshop version 3.2   */
/*          __SGI__      IRIS Indigo version                                */
/*          __PARIX__    PARIX Transputer version                           */
/*          __AIX__      AIX version (IBM)                                  */
/*          __PARAGON__  Intel Paragon version                              */
/*          __SUN4GCC__  Sun station 4 version                              */
/*          __HP__       HP Workstations                                    */
/*          __PC__       IBM compatible PC                                  */
/*          __T3D__      CRAY T3D version                                   */
/*          __T3E__      CRAY T3E version                                   */
/*          __POWERGC__  XPLORER (PowerPC)                                  */
/*          __CC__       CC-Parsytec (PowerPC)                              */
/*          __C90__      CRAY 90 version                                    */
/*          __YMP__      CRAY YMP version                                   */
/*          __NECSX4__   NEC SX4                                            */
/*          __SR2201__   Hitachi SR2201                                     */
/*          __MWCW__     Apple Power Macintosh                              */
/*          __MKLINUX__  Mikrokernel Linux (PowerPC)                        */
/*          __NEXTSTEP__ NEXTSTEP operating system                          */
/*                                                                          */
/* #define this if you are using NXLib                                      */
/*          __NXLIB__    NXLIB Paragon Library                              */
/*                                                                          */
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

#define ARCHNAME        "MPW32"

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

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)


/* Diese bloeden const pointer gehen immer noch nicht ! P.B. 27.7.95 */
#define const

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Silicon Graphics Workstations (old 32 bit version)       */
/*                                                                          */
/****************************************************************************/

#ifdef __SGI__
#undef __MWCW__

#define ARCHNAME        "SGI"

#include <time.h>

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

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym ## _

/* current time as DOUBLE value */
#define CURRENT_TIME   (((DOUBLE)clock())/((DOUBLE)CLOCKS_PER_SEC))

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for PARIX Transputer version                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __PARIX__
#undef __MWCW__

#define ARCHNAME        "PARIX"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD float
#define __SWAPBYTES__ 1

/* memory */
#define ALIGNMENT 4                                         /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFFC                    /* compatible to alignment			*/

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   (((DOUBLE)TimeNow())/((DOUBLE)CLOCK_TICK))

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for AIX                                                      */
/*                                                                          */
/****************************************************************************/


#ifdef __AIX__
#undef __MWCW__

#define ARCHNAME        "AIX"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float
#define __SWAPBYTES__ 1

/* memory */
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Intel Paragon version                                    */
/*                                                                          */
/****************************************************************************/

#ifdef __PARAGON__
#undef __MWCW__

#define ARCHNAME        "PARAGON"

#include <nx.h>

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float
#define __SWAPBYTES__ 1

/* memory */
#define ALIGNMENT 8                                         /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8                    /* compatible to alignment			*/

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym ## _

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)dclock())

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Solaris Version                                          */
/*                                                                          */
/****************************************************************************/

#ifdef __SOLARIS__
#undef __MWCW__

#define ARCHNAME        "SOLARIS"

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

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for Sun station 4 version                                    */
/*                                                                          */
/****************************************************************************/

#ifdef __SUN4GCC__
#undef __MWCW__

#define ARCHNAME        "SUN4GCC"

#include <stddef.h>

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Hewlett Packard HP700 series                             */
/*                                                                          */
/****************************************************************************/

#ifdef __HP__
#undef __MWCW__

#define ARCHNAME        "HP"

#include <time.h>

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

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   (((DOUBLE)clock())/((DOUBLE)CLOCKS_PER_SEC))

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Sun                                                      */
/*                                                                          */
/****************************************************************************/

#ifdef __SUN__
#undef __MWCW__

#define ARCHNAME        "SUN"

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

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for DEC                                                                                  */
/*                                                                          */
/****************************************************************************/

#ifdef __DEC__
#undef __MWCW__

#define ARCHNAME        "DEC"

/* basic types */
#define SHORT  short
#define INT    long
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD float

/* memory */
#define ALIGNMENT   8                   /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for IBM compatible Personal Computer                         */
/*                                                                          */
/****************************************************************************/

#ifdef __PC__
#undef __MWCW__

#define ARCHNAME        "PC"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float
#define __SWAPBYTES__ 1

/* memory */
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Microkernel Linux (on Apple PowerMacs)                   */
/*                                                                          */
/****************************************************************************/

#ifdef __MKLINUX__
#undef __MWCW__

#define ARCHNAME        "MKLINUX"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float
#define __SWAPBYTES__ 1

/* memory */
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY T3D                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __T3D__
#undef __MWCW__

#define ARCHNAME        "T3D"

#include <time.h>       /* for rtclock()   */
#include <unistd.h>     /* for _SC_CLK_TCK */

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  usym

/* current time as DOUBLE value */
#define CURRENT_TIME ((DOUBLE)(((double)rtclock())/((double)sysconf(_SC_CLK_TCK))))

#endif



/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY T3E                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __T3E__
#undef __MWCW__

#define ARCHNAME        "T3E"

#include <time.h>       /* for rtclock()   */
#include <unistd.h>     /* for _SC_CLK_TCK */

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  usym

/* current time as DOUBLE value */
#define CURRENT_TIME ((DOUBLE)(((double)rtclock())/((double)sysconf(_SC_CLK_TCK))))

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for XPLORER (PowerPC)  version                               */
/*                                                                          */
/****************************************************************************/

#if defined(__POWERGC__)
#undef __MWCW__

#define ARCHNAME        "POWERGC"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8               /* power of 2 and >=sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8     /*  compatible to alignment */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif




/****************************************************************************/
/*                                                                          */
/* Definitions for CC-Parsytec                                              */
/*                                                                          */
/****************************************************************************/

#ifdef __CC__
#undef __MWCW__

#define ARCHNAME        "CC"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 4               /* power of 2 and >=sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFFC     /*  compatible to alignment */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY 90                                                  */
/*                                                                          */
/****************************************************************************/

#ifdef __C90__
#undef __MWCW__

#define ARCHNAME        "C90"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif



/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY YMP                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __YMP__
#undef __MWCW__

#define ARCHNAME        "YMP"

/* basic types */
#define SHORT  short
#define INT    int
#define FLOAT  float
#define DOUBLE double
#define COORD  float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for NEC SX4                                                  */
/*                                                                          */
/****************************************************************************/

#ifdef __NECSX4__
#undef __MWCW__

#define ARCHNAME        "NECSX4"

#include <time.h>


/* basic types */
#define SHORT  short
#define INT    long    /* sizeof(int) != sizeof(void *) !! */
#define FLOAT  float
#define DOUBLE double
#define COORD  double
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym ## _

/* current time as DOUBLE value */
#define CURRENT_TIME   (((DOUBLE)clock())/((DOUBLE)CLOCKS_PER_SEC))

#endif



/****************************************************************************/
/*                                                                          */
/* Definitions for Hitachi SR2201                                           */
/*                                                                          */
/****************************************************************************/

#ifdef __SR2201__
#undef __MWCW__

#define ARCHNAME        "SR2201"

/* basic types */
#define SHORT  short
#define INT    long    /* sizeof(int) != sizeof(void *) !! */
#define FLOAT  float
#define DOUBLE double
#define COORD  double
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT 8                     /* power of 2 and >= sizeof(int) !  */
#define ALIGNMASK 0xFFFFFFF8            /* compatible to alignment          */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym ## _

/* current time as DOUBLE value */
#define CURRENT_TIME   (((DOUBLE)clock())/((DOUBLE)CLOCKS_PER_SEC))

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for NEXTSTEP                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __NEXTSTEP__
#undef __MWCW__

#define ARCHNAME        "NEXTSTEP"

/* basic types */
#define SHORT         short
#define INT           int
#define FLOAT         float
#define DOUBLE        double
#define COORD         float
#define SCREEN_COORD  float

/* memory */
#define ALIGNMENT     4             /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK     0xFFFFFFFC    /* compatible to alignment */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000
#endif
#define CURRENT_TIME   (((DOUBLE)times(NULL))/((DOUBLE)CLK_TCK))

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Apple Power Macintosh                                    */
/*                                                                          */
/****************************************************************************/

#ifdef __MWCW__

#define ARCHNAME        "MWCW"

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
/*#define malloc(n) ((void *) NewPtr((Size) n))
 #define free(p) DisposPtr((Ptr) p)*/
#define ALIGNMENT 4                     /* power of 2 and >= sizeof(int) ! */
#define ALIGNMASK 0xFFFFFFFC            /* compatible to alignment */

/* fortran interfacing */
#define F77SYM(lsym,usym)  lsym

/* current time as DOUBLE value */
#define CURRENT_TIME   ((DOUBLE)0.0)

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for NXLIB 1.1 version                                        */
/*                                                                          */
/****************************************************************************/

#ifdef __NXLIB__

#include <nxmalloc.h>                   /* redefine malloc and related calls*/

#endif


/****************************************************************************/
/*                                                                          */
/* some general definitions						                            */
/*                                                                          */
/****************************************************************************/

#define ARCH_VERSION "ARCH_1_0"
static char compilerrcs_id[] = "$Id$";

/* limits of the basic types */
#define MAX_S            SHRT_MAX
#define MAX_I            INT_MAX
#define MAX_F            FLT_MAX
#define SMALL_F         (FLT_EPSILON*SMALL_FAC)
#define MAX_D            DBL_MAX
#define SMALL_D         (DBL_EPSILON*SMALL_FAC)
#define MAX_C            FLT_MAX
#define SMALL_C         (FLT_EPSILON*SMALL_FAC)


#ifdef __cplusplus
}
#endif
#endif
