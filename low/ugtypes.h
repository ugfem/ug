// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* $Id$ */

/*

   globally set data types

 */

#ifndef UGTYPES_H
#define UGTYPES_H

#include "namespace.h"

START_UG_NAMESPACE

#ifdef AUTOTOOLS_BUILD
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

#else  /* AUTOTOOLS_BUILD */

/****************************************************************************/
/*                                                                          */
/* Definitions for Apple MacIntosh                                          */
/*                                                                          */
/****************************************************************************/

#ifdef __MPW32__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Silicon Graphics Workstations (old 32 bit version)       */
/*                                                                          */
/****************************************************************************/

#ifdef __SGI__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for Silicon Graphics Workstations (new 64 bit version)       */
/*                                                                          */
/****************************************************************************/

#ifdef __SGI10__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for PARIX Transputer version                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __PARIX__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for AIX                                                      */
/*                                                                          */
/****************************************************************************/


#ifdef __AIX__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Intel Paragon version                                    */
/*                                                                          */
/****************************************************************************/

#ifdef __PARAGON__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Solaris Version                                          */
/*                                                                          */
/****************************************************************************/

#ifdef __SOLARIS__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for Sun station 4 version                                    */
/*                                                                          */
/****************************************************************************/

#ifdef __SUN4GCC__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Hewlett Packard HP700 series                             */
/*                                                                          */
/****************************************************************************/

#ifdef __HP__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Sun                                                      */
/*                                                                          */
/****************************************************************************/

#ifdef __SUN__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for DEC                                                      */
/*                                                                          */
/****************************************************************************/

#ifdef __DEC__

typedef short SHORT;
typedef long INT;
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for IBM compatible Personal Computer                         */
/*                                                                          */
/****************************************************************************/

#ifdef __PC__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Microkernel Linux (on Apple PowerMacs)                   */
/*                                                                          */
/****************************************************************************/

#ifdef __MKLINUX__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for LINUXAXP                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __LINUXAXP__

typedef short SHORT;
typedef long INT;                    /* sizeof(int) != sizeof(void *) !! */
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for LINUXIA64                                                */
/*                                                                          */
/****************************************************************************/

#ifdef __LINUXIA64__

typedef short SHORT;
typedef long INT;                    /* sizeof(int) != sizeof(void *) !! */
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for Cygwin                                                   */
/*                                                                          */
/****************************************************************************/

#ifdef __CYGWIN__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for AMD64                                                    */
/*                                                                          */
/****************************************************************************/

#ifdef __AMD64__

typedef short SHORT;
typedef long INT;                   /* sizeof(int) != sizeof(void *) !! */
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for LINUXPPC                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __LINUXPPC__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY T3D                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __T3D__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif



/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY T3E                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __T3E__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for XPLORER (PowerPC)  version                               */
/*                                                                          */
/****************************************************************************/

#ifdef __POWERGC__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif




/****************************************************************************/
/*                                                                          */
/* Definitions for CC-Parsytec                                              */
/*                                                                          */
/****************************************************************************/

#ifdef __CC__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY 90                                                  */
/*                                                                          */
/****************************************************************************/

#ifdef __C90__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif



/****************************************************************************/
/*                                                                          */
/* Definitions for CRAY YMP                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __YMP__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for NEC SX4                                                  */
/*                                                                          */
/****************************************************************************/

#ifdef __NECSX4__

typedef short SHORT;
typedef long INT;   /* sizeof(int) != sizeof(void *) !! */
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef double COORD;
typedef float SCREEN_COORD;

#endif



/****************************************************************************/
/*                                                                          */
/* Definitions for Hitachi SR2201                                           */
/*                                                                          */
/****************************************************************************/

#ifdef __SR2201__

typedef short SHORT;
typedef long INT;    /* sizeof(int) != sizeof(void *) !! */
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef double COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Hitachi SR8000                                           */
/*                                                                          */
/****************************************************************************/

#ifdef __SR8K__

typedef short SHORT;
typedef long INT;    /* sizeof(int) != sizeof(void *) !! */
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef double COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for MACOSXSERVER                                             */
/*                                                                          */
/****************************************************************************/

#ifdef __MACOSXSERVER__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Mac OS X                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __MACOSX__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif

#ifdef __PPC64__

typedef short SHORT;
typedef long INT;                     /* sizeof(int) != sizeof(void *) !! */
typedef unsigned long UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


/****************************************************************************/
/*                                                                          */
/* Definitions for Apple Macintosh with CodeWarrior compiler                */
/*                                                                          */
/****************************************************************************/

#ifdef __MWCW__

typedef short SHORT;
typedef int INT;
typedef unsigned int UINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef float COORD;
typedef float SCREEN_COORD;

#endif


#endif  /* AUTOTOOLS_BUILD */

END_NAMESPACE

#endif
