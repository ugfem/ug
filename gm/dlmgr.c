// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      dlmgr.c                                                       */
/*                                                                          */
/* Purpose:   functions for dynamic linked list management                  */
/*                                                                          */
/* Author:    Stefan Lang                                                   */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            email: stefan@ica3.uni-stuttgart.de                           */
/*            phone: 0049-(0)711-685-7003                                   */
/*            fax  : 0049-(0)711-685-7000                                   */
/*                                                                          */
/* History:   961216 sl  start of dynamic list management					*/
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

#include "dlmgr.h"
#include "gm.h"

#ifdef ModelP
#include "parallel.h"
#endif

/* element */
#define OTYPE  ELEMENT
#define PRED   PREDE
#define SUCC   SUCCE
#include "dlmgr.t"

/* node */
#define OTYPE  NODE
#define PRED   PREDN
#define SUCC   SUCCN
#include "dlmgr.t"

/* vertex */
#define OTYPE  VERTEX
#define PRED   PREDV
#define SUCC   SUCCV
#include "dlmgr.t"

/* vertex */
#define OTYPE  VECTOR
#define PRED   PREDVC
#define SUCC   SUCCVC
#include "dlmgr.t"
#undef OTYPE
#undef PRED
#undef SUCC
