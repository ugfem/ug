// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      dlmgr.h                                                       */
/*                                                                          */
/* Purpose:   defines for dynamic linked list management                    */
/*                                                                          */
/* Author:    Stefan Lang                                                   */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                                   */
/*                                                                          */
/* History:   960915 sl  start of dynamic list management					*/
/*                                                                          */
/* Remarks:                                                                 */
/*            Management of dynamic linked lists, which consist of          */
/*            several parts. An object can be mapped to a part of the list  */
/*            by its priority using PRIO2LISTPART().                        */
/*            The formulation on object basis allows for management of      */
/*            elements, nodes, vectors and vertices.                        */
/*            A list has the form:                                          */
/*                p0first-p0last->p1first-p1last->...->pnfirst..pnlast      */
/*            where each part p0,p1,..,pn is limited by two pointers,       */
/*            pxfirst and pxlast, x=0..n. The part numbers are ordered      */
/*            increasingly and connected in a manner that one can run       */
/*            through the whole list in increasing order (SUCC) but only    */
/*            through one listpart in decreasing order (PRED).              */
/*            Linking/Unlinking of objects in a list part is done in a      */
/*            way that preserves these conventions.                         */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

#ifndef __DLMGR_H__

#include "gm.h"
#include "misc.h"
#include "debug.h"

#ifdef ModelP
#include "parallel.h"
#endif

#ifdef ModelP

#define FIRSTPART_OF_LIST               0
#define LASTPART_OF_LIST(OTYPE) ((CAT(OTYPE,PRIOS)) -1)

/* define DDD_HDR macros */
#define HDRELEMENT      PARHDRE
#define HDRNODE         PARHDR
#define HDRVERTEX       PARHDRV
#define HDRVECTOR       PARHDR
#define HDR(OTYPE)      CAT(HDR,OTYPE)

/* define macros for formatted output */
#define ELEMENTFMT      EID_
#define NODEFMT         ID_
#define VERTEXFMT       VID_
#define VECTORFMT       VINDEX_
#define FORMAT(t)       CAT(t,FMT)

#endif

/* define Object COUNTER macros */
#define COUNTELEMENT    NT
#define COUNTNODE               NN
#define COUNTVERTEX             NV
#define COUNTVECTOR             NVEC

/* define header prototypes */
#define UNLINK(OTYPE)    void CAT(GRID_UNLINK_, OTYPE ) (GRID *Grid, OTYPE *Object)
#define LINK(OTYPE)      void CAT(GRID_LINK_,OTYPE) (GRID *Grid, OTYPE *Object, INT Prio)
#define LINKX(OTYPE)     void CAT(GRID_LINKX_,OTYPE) (GRID *Grid, OTYPE *Object, INT Prio, OTYPE *After)
#define INIT(OTYPE)      void CAT3(GRID_INIT_,OTYPE,_LIST(GRID *Grid))
#define CHECK(OTYPE)     void CAT3(GRID_CHECK_,OTYPE,_LIST(GRID *Grid))
#ifdef ModelP
#define PRINT_LIST(OTYPE) void CAT(PRINT_LIST_STARTS_,OTYPE) (GRID *Grid, INT prios)
#endif

LINK(ELEMENT);
LINKX(ELEMENT);
UNLINK(ELEMENT);
INIT(ELEMENT);
CHECK(ELEMENT);

LINK(NODE);
LINKX(NODE);
UNLINK(NODE);
INIT(NODE);
CHECK(NODE);

LINK(VERTEX);
LINKX(VERTEX);
UNLINK(VERTEX);
INIT(VERTEX);
CHECK(VERTEX);

LINK(VECTOR);
LINKX(VECTOR);
UNLINK(VECTOR);
INIT(VECTOR);
CHECK(VECTOR);

#ifdef ModelP
PRINT_LIST(ELEMENT);
PRINT_LIST(NODE);
PRINT_LIST(VERTEX);
PRINT_LIST(VECTOR);
#endif

#endif /* __DLMGR_H__ */
