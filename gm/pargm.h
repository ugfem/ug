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
  PrioGhost    = 1,
  PrioVGhost   = 2,
  PrioBorder   = 3,
  PrioMaster   = 4
};

/* define dynamic lists ids */
#define ELEMENT_LIST     0
#define NODE_LIST        1
#define VECTOR_LIST      2
#define VERTEX_LIST      3

#ifdef ModelP

/* define the number of priorities for objects */
#define MAXPRIOS                        8
#define ELEMENTPRIOS            2
#define NODEPRIOS                       3
#define VECTORPRIOS                     3
#define VERTEXPRIOS                     3

/* define mapping from object priority to position in linked list */
#define PRIO2LISTPART(listtype,prio) \
  ((listtype == ELEMENT_LIST) ? ((prio == PrioGhost) ? 0 :\
                                 (prio == PrioVGhost) ? 0 :\
                                 (prio == PrioMaster) ? 1 : -1) :\
   ((prio == PrioGhost) ? 0 : (prio ==PrioVGhost) ? 0 :\
      (prio == PrioBorder) ? 2 : (prio == PrioMaster) ? 2 : -1))

/* define mapping from position in linked list to object priority */
#define LISTPART2PRIO(listtype,listpart,prios)                               \
  {                                                                        \
    INT Entry;                                                           \
    for (Entry=0; Entry<MAXPRIOS; Entry++) prios[Entry] = -1;            \
    Entry = 0;                                                           \
    if (listtype == ELEMENT_LIST)                                        \
    {                                                                    \
      if (listpart == 0)                                               \
      {                                                                \
        prios[Entry++] = PrioGhost;                                  \
        prios[Entry++] = PrioVGhost;                                 \
      }                                                                \
      else if(listpart == 1) prios[Entry++] = PrioMaster;             \
    }                                                                    \
    else                                                                 \
    {                                                                    \
      if (listpart == 0)                                               \
      {                                                                \
        prios[Entry++] = PrioGhost;                                  \
        prios[Entry++] = PrioVGhost;                                 \
      }                                                                \
      else if (listpart == 2)                                          \
      {                                                                \
        prios[Entry++] = PrioBorder;                                 \
        prios[Entry++] = PrioMaster;                                 \
      }                                                                \
    }                                                                    \
  }

/* define mapping from element priority to index in son array of father */
#define PRIO2INDEX(prio) \
  ((prio==PrioGhost || prio==PrioVGhost) ? 1 : \
   (prio == PrioMaster) ? 0 : -1)

/* map pointer to structure onto a pointer to its DDD_HDR */
#define PARHDR(obj)    (&((obj)->ddd))


#else   /* not ModelP */

/* define the number of prioities for objects */
#define MAXPRIOS                        1
#define ELEMENTPRIOS            1
#define NODEPRIOS                       1
#define VECTORPRIOS                     1
#define VERTEXPRIOS                     1

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
#define ID_FMTX     "%x/%ld/%08x/%d"
#define ID_FFMTX    "%x/%9ld/%08x/%02d"
#define ID_PRTX(x)  x,((long)ID(x)),GID(x),PRIO(x)

#define VID_FMT     ID_FMT
#define VID_FFMT    ID_FFMT
#define VID_PRT(x)  ((long)ID(x)),VXGID(x)
#define VID_FMTE    ID_FMTE
#define VID_FFMTE   ID_FFMTE
#define VID_PRTE(x) ((long)ID(x)),VXGID(x),VXPRIO(x)
#define VID_FMTX    ID_FMTX "/%d"
#define VID_FFMTX   ID_FFMTX "/%d"
#define VID_PRTX(x) x,((long)ID(x)),VXGID(x),VXPRIO(x),LEVEL(x)

#define EID_FMT     ID_FMT
#define EID_FFMT    ID_FFMT
#define EID_PRT(x)  ((long)ID(x)),EGID(x)
#define EID_FMTE    ID_FMTE
#define EID_FFMTE   ID_FFMTE
#define EID_PRTE(x) ((long)ID(x)),EGID(x),EPRIO(x)
#define EID_FMTX    ID_FMTX "/%d"
#define EID_FFMTX   ID_FFMTX "/%d"
#define EID_PRTX(x) x,((long)ID(x)),EGID(x),EPRIO(x),TAG(x)

#define VINDEX_FMT     ID_FMT
#define VINDEX_FFMT    ID_FFMT
#define VINDEX_PRT(x)  ((long)VINDEX(x)),GID(x)
#define VINDEX_FMTE    ID_FMTE
#define VINDEX_FFMTE   ID_FFMTE
#define VINDEX_PRTE(x) ((long)VINDEX(x)),GID(x),PRIO(x)
#define VINDEX_FMTX    ID_FMTX
#define VINDEX_FFMTX   ID_FFMTX
#define VINDEX_PRTX(x) x,((long)VINDEX(x)),GID(x),PRIO(x)

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

#define MASTER(p)       1
#define GHOST(p)        0
#define SETEPRIO(p,i)   ;


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
extern DDD_IF BorderNodeIF, BorderNodeSymmIF, OuterNodeIF, NodeVIF;
extern DDD_IF BorderVectorIF, BorderVectorSymmIF, OuterVectorIF,
              VectorVIF, VectorVAllIF;
#ifdef __THREEDIM__
extern DDD_IF EdgeIF, BorderEdgeSymmIF, EdgeHIF;
#endif
#endif


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* functions implemented in parallel/dddif/support.c */
DOUBLE UG_GlobalSumDOUBLE (DOUBLE);
INT UG_GlobalSumINT (INT);
INT UG_GlobalMaxINT (INT);
INT UG_GlobalMinINT (INT);
void UG_GlobalSumNDOUBLE (INT, DOUBLE *);
DOUBLE UG_GlobalMaxDOUBLE (DOUBLE i);
DOUBLE UG_GlobalMinDOUBLE (DOUBLE i);

#endif /* __PARGM_H__ */
