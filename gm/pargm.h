// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  pargm.h														*/
/*																			*/
/* Purpose:   defines for parallel grid manager                                                         */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: stefan@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   960410 kb  created from parallel.h					        */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __PARGM_H__
#define __PARGM_H__


#ifdef ModelP
#include "ppif.h"
#include "ddd.h"
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* object priorties */
enum Priorities
{
  PrioNone     = 0,
  PrioHGhost   = 1,
  PrioVGhost   = 2,
  PrioVHGhost  = 3,
  PrioBorder   = 4,
  PrioMaster   = 5
};

/* define dynamic lists ids */
#define ELEMENT_LIST     0
#define NODE_LIST        1
#define VECTOR_LIST      2
#define VERTEX_LIST      3

#ifdef ModelP

/* define number of priorities for objects */
#define ELEMENT_PRIOS                   4
#define NODE_PRIOS                              5
#define VECTOR_PRIOS                    5
#define VERTEX_PRIOS                    5

/* define number of listparts for objects */
#define MAX_LISTPARTS                   8
#define ELEMENT_LISTPARTS               2
#define NODE_LISTPARTS                  3
#define VECTOR_LISTPARTS                3
#define VERTEX_LISTPARTS                3

/* define mapping from object priority to position in linked list */
#define PRIO2LISTPART(listtype,prio)                                         \
  ((listtype == ELEMENT_LIST) ? ((prio == PrioHGhost) ? 0 :                \
                                 (prio == PrioVGhost) ? 0 : (prio == PrioVHGhost) ? 0 :               \
                                 (prio == PrioMaster) ? 1 : -1) :                                 \
   ((prio == PrioHGhost) ? 0 : (prio ==PrioVGhost) ? 0 :            \
      (prio == PrioVHGhost) ? 0 :                                  \
      (prio == PrioBorder) ? 2 : (prio == PrioMaster) ? 2 : -1))

/* define mapping from position in linked list to object priority */
#define LISTPART2PRIO(listtype,listpart,prios)                               \
  {                                                                        \
    INT Entry;                                                           \
    for (Entry=0; Entry<MAX_LISTPARTS; Entry++) prios[Entry] = -1;       \
    Entry = 0;                                                           \
    if (listtype == ELEMENT_LIST)                                        \
    {                                                                    \
      if (listpart == 0)                                               \
      {                                                                \
        prios[Entry++] = PrioHGhost;                                 \
        prios[Entry++] = PrioVGhost;                                 \
        prios[Entry++] = PrioVHGhost;                                \
      }                                                                \
      else if(listpart == 1) prios[Entry++] = PrioMaster;             \
    }                                                                    \
    else                                                                 \
    {                                                                    \
      if (listpart == 0)                                               \
      {                                                                \
        prios[Entry++] = PrioHGhost;                                 \
        prios[Entry++] = PrioVGhost;                                 \
        prios[Entry++] = PrioVHGhost;                                \
      }                                                                \
      else if (listpart == 2)                                          \
      {                                                                \
        prios[Entry++] = PrioBorder;                                 \
        prios[Entry++] = PrioMaster;                                 \
      }                                                                \
    }                                                                    \
  }

/* define mapping from element priority to index in son array of father */
#define PRIO2INDEX(prio)                                                     \
  ((prio==PrioHGhost || prio==PrioVGhost || prio==PrioVHGhost) ? 1 :       \
   (prio == PrioMaster) ? 0 : -1)

/* map pointer to structure onto a pointer to its DDD_HDR */
#define PARHDR(obj)    (&((obj)->ddd))


#else   /* not ModelP */

/* define number of priorities for objects */
#define ELEMENT_PRIOS                   4
#define NODE_PRIOS                              5
#define VECTOR_PRIOS                    5
#define VERTEX_PRIOS                    5

/* define number of listparts for objects */
#define MAX_LISTPARTS                   1
#define ELEMENT_LISTPARTS               1
#define NODE_LISTPARTS                  1
#define VECTOR_LISTPARTS                1
#define VERTEX_LISTPARTS                1

/* define mapping from object priority to position in linked list */

/* define mapping from object priority to position in linked list */
#define PRIO2LISTPART(listtype,prio) 0

/* define mapping from position in linked list to object priority */
#define LISTPART2PRIO(listtype,listpart,prios) 0

/* define mapping from position in linked list to object priority */
#define PRIO2INDEX(prio)  0

#endif



/*
        printing of IDs via printf()

        use ID_FMT as format string, and ID_PRT as macro for printing.
        example of usage:

                printf("NodeId " ID_FMT "\n", ID_PRT(theNode));

        ID_FFMT (fixed format) is a version of ID_FMT with fixed width.
        ID_FFMTE (extended) is a version of ID_FFMT with additional information.

        in ModelP, additionaly the DDD_GlobalID is printed for each object.

        NOTE: for vertices and elements, one must use the VID_ and EID_ macros,
              respectively. this is due to differences in data structures
              (unions -> PARHDRV/PARHDRE ->VID_/EID_)
 */

#ifdef ModelP

#define ID_FMT      "%ld/%08x"
#define ID_FFMT     "%9ld/%08x"
#define ID_PRT(x)   ((long)ID(x)),GID(x)
#define ID_FMTE     "%ld/%08x/%d"
#define ID_FFMTE    "%9ld/%08x/%02d"
#define ID_PRTE(x)  ((long)ID(x)),GID(x),PRIO(x)
#define ID_FMTX     "%d/%ld/%08x/%d"
#define ID_FFMTX    "%x/%9ld/%08x/%02d"
#define ID_PRTX(x)  KeyForObject((KEY_OBJECT *)x),((long)ID(x)),GID(x),PRIO(x)

#define VID_FMT     ID_FMT
#define VID_FFMT    ID_FFMT
#define VID_PRT(x)  ((long)ID(x)),VXGID(x)
#define VID_FMTE    ID_FMTE
#define VID_FFMTE   ID_FFMTE
#define VID_PRTE(x) ((long)ID(x)),VXGID(x),VXPRIO(x)
#define VID_FMTX    ID_FMTX "/%d"
#define VID_FFMTX   ID_FFMTX "/%d"
#define VID_PRTX(x)     KeyForObject((KEY_OBJECT *)x),((long)ID(x)),VXGID(x),VXPRIO(x),LEVEL(x)

#define EID_FMT     ID_FMT
#define EID_FFMT    ID_FFMT
#define EID_PRT(x)  ((long)ID(x)),EGID(x)
#define EID_FMTE    ID_FMTE
#define EID_FFMTE   ID_FFMTE
#define EID_PRTE(x) ((long)ID(x)),EGID(x),EPRIO(x)
#define EID_FMTX    ID_FMTX "/%d/%d/%d/%d"
#define EID_FFMTX   ID_FFMTX "/%d"
#define EID_PRTX(x) KeyForObject((KEY_OBJECT *)x),((long)ID(x)),EGID(x),EPRIO(x),TAG(x),\
  LEVEL(x),ECLASS(x),REFINECLASS(x)

#define VINDEX_FMT     ID_FMT
#define VINDEX_FFMT    ID_FFMT
#define VINDEX_PRT(x)  ((long)VINDEX(x)),GID(x)
#define VINDEX_FMTE    ID_FMTE
#define VINDEX_FFMTE   ID_FFMTE
#define VINDEX_PRTE(x) ((long)VINDEX(x)),GID(x),PRIO(x)
#define VINDEX_FMTX    ID_FMTX
#define VINDEX_FFMTX   ID_FFMTX
#define VINDEX_PRTX(x) KeyForObject((KEY_OBJECT *)x),((long)VINDEX(x)),GID(x),PRIO(x)

#ifdef __TWODIM__
#define EDID_FMT     "%08x"
#define EDID_FFMT    EDID_FMT
#define EDID_PRT(x)  (x)
#define EDID_FMTE    "%08x"
#define EDID_FFMTE   EDID_FMTE
#define EDID_PRTE(x) (x)
#define EDID_FMTX    "%08x"
#define EDID_FFMTX   EDID_FMTX
#define EDID_PRTX(x) (x)
#endif

#ifdef __THREEDIM__
#define EDID_FMT     "%08x"
#define EDID_FFMT    EDID_FMT
#define EDID_PRT(x)  GID(x)
#define EDID_FMTE    "%08x/%d"
#define EDID_FFMTE   EDID_FMTE
#define EDID_PRTE(x) GID(x),PRIO(x)
#define EDID_FMTX    "%x/%08x/%d"
#define EDID_FFMTX   EDID_FMTX
#define EDID_PRTX(x) x,GID(x),PRIO(x)
#endif


#define PFMT            "%3d:"

/* PAR/ENDPAR parallel preprocessor statements   */
/* to select code only valid for ModelP          */
#define PAR(x)          x
#define ENDPAR

#else

#define ID_FMT      "%ld"
#define ID_FFMT     "%9ld"
#define ID_PRT(x)   ((long)ID(x))
#define ID_FMTE     "%ld"
#define ID_FFMTE    "%9ld"
#define ID_PRTE(x)  ID_PRT(x)
#define ID_FMTX     "%ld"
#define ID_FFMTX    "%9ld"
#define ID_PRTX(x)  ID_PRT(x)

#define VID_FMT     ID_FMT
#define VID_FFMT    ID_FFMT
#define VID_PRT(x)  ID_PRT(x)
#define VID_FMTE    ID_FMTE
#define VID_FFMTE   ID_FFMTE
#define VID_PRTE(x) VID_PRT(x)
#define VID_FMTX    ID_FMTX
#define VID_FFMTX   ID_FFMTX
#define VID_PRTX(x) VID_PRT(x)

#define EID_FMT     ID_FMT
#define EID_FFMT    ID_FFMT
#define EID_PRT(x)  ID_PRT(x)
#define EID_FMTE    ID_FMTE
#define EID_FFMTE   ID_FFMTE
#define EID_PRTE(x) EID_PRT(x)
#define EID_FMTX    ID_FMTX
#define EID_FFMTX   ID_FFMTX
#define EID_PRTX(x) EID_PRT(x)

#define VINDEX_FMT     ID_FMT
#define VINDEX_FFMT    ID_FFMT
#define VINDEX_PRT(x)  ((long)VINDEX(x))
#define VINDEX_FMTE    ID_FMTE
#define VINDEX_FFMTE   ID_FFMTE
#define VINDEX_PRTE(x) VINDEX_PRT(x)
#define VINDEX_FMTX    ID_FMTX
#define VINDEX_FFMTX   ID_FFMTX
#define VINDEX_PRTX(x) VINDEX_PRT(x)

#define EDID_FMT     "%08x"
#define EDID_FFMT    EDID_FMT
#define EDID_PRT(x)  (x)
#define EDID_FMTE    "%08x"
#define EDID_FFMTE   EDID_FMTE
#define EDID_PRTE(x) (x)
#define EDID_FMTX    "%08x"
#define EDID_FFMTX   EDID_FMTX
#define EDID_PRTX(x) (x)

#define PFMT            "%1d:"
#define GIDFMT          "%1d"

/* dummies for global id */
#define EGID(e)     ID(e)
#define GID(e)      ((OBJT(e)==VEOBJ) ? VINDEX((VECTOR *)e) : ID(e))
#define VGID(e)     ID(e)

/* PAR/ENDPAR parallel preprocessor statements   */
/* for serial case expanded to code x is ignored */
#define PAR(x)
#define ENDPAR

#define GetAllSons(e,s)         GetSons(e,s)

/* dummy defines for serial case according to parallel defines in parallel.h */
/* dummies for elements */
#define EMASTER(p)              1
#define EGHOST(p)               0
#define EHGHOST(p)              0
#define EVGHOST(p)              0
#define EPRIO(p)                0
#define SETEPRIO(p,i)   ;
#define EMASTERPRIO(p)  1
#define EPROCLIST(p)    (&_proclist_)
#define ENCOPIES(p)             1
#define PARTITION(p)    _partition_

/* dummies for nodes, vectors, edges */
#define MASTER(p)               1
#define GHOST(p)                0
#define HGHOST(p)               0
#define VGHOST(p)               0
#define PRIO(p)                 0
#define SETPRIO(p,i)    ;
#define PROCLIST(p)             (&_proclist_)
#define NCOPIES(p)              1

/* dummies for vertices */
#define SETVXPRIO(e,p)  ;

/* ddd dummies */
#define DDD_OBJ                 void *
#define DDD_IdentifyBegin()
#define DDD_IdentifyEnd()
#define DDD_IdentifyNumber(o,p,n)
#define DDD_IFAOneway(p1,p2,p3,p4,p5,p6)
#define DDD_PrioritySet(e,p)

/* ppif dummies */
#define Broadcast(p,n)  ((int)0)

/* dummys for reduction functions implemented in parallel/dddif/support.c */
#define UG_GlobalSumINT(x)              x
#define UG_GlobalMaxINT(x)              x
#define UG_GlobalMinINT(x)              x
#define UG_GlobalSumNINT(x,y)
#define UG_GlobalMaxNINT(x,y)
#define UG_GlobalMinNINT(x,y)
#define UG_GlobalSumDOUBLE(x)   x
#define UG_GlobalMaxDOUBLE(x)   x
#define UG_GlobalMinDOUBLE(x)   x
#define UG_GlobalSumNDOUBLE(x,y)
#define UG_GlobalMaxNDOUBLE(x,y)
#define UG_GlobalMinNDOUBLE(x,y)
#endif


/****************************************************************************/
/*                                                                          */
/* exported global variables                                                */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP
/* DDD Interfaces */
extern DDD_IF ElementIF, ElementSymmIF, ElementVIF, ElementSymmVIF,
              ElementVHIF, ElementSymmVHIF;
extern DDD_IF BorderNodeIF, BorderNodeSymmIF, OuterNodeIF, NodeVIF,
              NodeIF, NodeAllIF;
extern DDD_IF BorderVectorIF, BorderVectorSymmIF,
              OuterVectorIF, OuterVectorSymmIF,
              VectorVIF, VectorVAllIF, VectorAllIF;
#ifdef __THREEDIM__
extern DDD_IF EdgeIF, BorderEdgeSymmIF, EdgeHIF, EdgeAllIF;
#endif
#endif


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* functions implemented in parallel/dddif/support.c */
#ifdef ModelP
INT    UG_GlobalSumINT     (INT x);
INT    UG_GlobalMaxINT     (INT x);
INT    UG_GlobalMinINT     (INT x);
void   UG_GlobalSumNINT    (INT n, INT *x);
void   UG_GlobalMaxNINT    (INT n, INT *x);
void   UG_GlobalMinNINT    (INT n, INT *x);
DOUBLE UG_GlobalSumDOUBLE  (DOUBLE i);
DOUBLE UG_GlobalMaxDOUBLE  (DOUBLE i);
DOUBLE UG_GlobalMinDOUBLE  (DOUBLE i);
void   UG_GlobalSumNDOUBLE (INT n, DOUBLE *x);
void   UG_GlobalMaxNDOUBLE (INT n, DOUBLE *x);
void   UG_GlobalMinNDOUBLE (INT n, DOUBLE *x);
#endif

#endif /* __PARGM_H__ */
