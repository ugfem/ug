// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File: block.h                                                            */
/*                                                                          */
/* Purpose: block num proc header                                           */
/*                                                                          */
/* Author: Klaus Johannsen                                                  */
/*         Sit                                                              */
/*         Universitaet Heidelberg                                          */
/*         INF 368                                                          */
/*         69120 Heidelberg                                                 */
/*         email: ug@ica3.uni-stuttgart.de                                  */
/*                                                                          */
/* History:   Sep 27, 2004 begin                                            */
/*                                                                          */
/* Remarks: not finished!                                                   */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __BLOCKING__
#define __BLOCKING__

#include "np.h"
#include "heaps.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*     compile time constants defining static data size (i.e. arrays)       */
/*     other constants                                                      */
/*     macros                                                               */
/*                                                                          */
/****************************************************************************/

#define BLOCKING_CLASS_NAME "blocking"

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef void *(*GetMemProcPtr)(NS_PREFIX MEM n);

typedef struct
{
  INT n;                                 /* number of blocks                */
  INT *nb;                               /* number of block entries         */
  VECTOR ***vb;                          /* VECTOR-blocks                   */
} BLOCKING_STRUCTUR;

struct np_blocking {
  NP_BASE base;                          /* inherits base class             */

  /* data (optional, necessary for calling the generic execute routine)   */
  MATDATA_DESC *A;                       /* matrix symbol                   */

  /* functions */
  INT (*PreProcess)
    (struct np_blocking *,               /* pointer to (derived) object     */
    INT level,                           /* level to block                  */
    INT *);                              /* result                          */
  INT (*Blocking)
    (struct np_blocking *,               /* pointer to (derived) object     */
    GetMemProcPtr GetMem,                /* memory allocation               */
    INT level,                           /* level to block                  */
    MATDATA_DESC *A,                     /* matrix                          */
    BLOCKING_STRUCTUR *bs,               /* block structur                  */
    INT *);                              /* result                          */
  INT (*PostProcess)
    (struct np_blocking *,               /* pointer to (derived) object     */
    INT level,                           /* level to block                  */
    INT *);                              /* result                          */
};
typedef struct np_blocking NP_BLOCKING;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT InitBlocking (void);

END_UGDIM_NAMESPACE

#endif
