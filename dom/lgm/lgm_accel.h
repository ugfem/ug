// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*****************************************************************************
* File:      lgm_accel.h                                                    *
* Purpose:   Accelerate geometry access by using bounding box trees         *
*                                                                           *
* Author:	  O. Sterz                                                       *
*                                                                           *
* History:   Nov 2002                                                       *
* Remarks:                                                                  *
*****************************************************************************/

/*****************************************************************************
* auto include mechanism and other include files                            *
*****************************************************************************/
#ifndef __LGM_ACCEL__
#define __LGM_ACCEL__

#include "namespace.h"

START_NAMESPACE


/*****************************************************************************
* defines in the following order:                                           *
*        compile time constants defining static data size (i.e. arrays)     *
*        other constants                                                    *
*        macros                                                             *
*****************************************************************************/
#define DebugLGMAccel 0

/*****************************************************************************
* exported data structures                                                  *
*****************************************************************************/

/*****************************************************************************
* exported global variables                                                 *
*****************************************************************************/

/*****************************************************************************
* public function declarations                                              *
*****************************************************************************/
INT LGM_InitAcceleration(HEAP *theHeap, LGM_SURFACE **sf, INT nsf);

END_NAMESPACE

#endif
