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
  PrioGhost    = 2,
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
#define ELEMENTPRIOS            2
#define NODEPRIOS                       3
#define VECTORPRIOS                     3
#define VERTEXPRIOS                     3

/* define mapping from object priority to position in linked list */
#define PRIO2LISTPART(listtype,prio) \
  ((listtype == ELEMENT_LIST) ? ((prio == PrioGhost) ? 0 : (prio == PrioMaster) ? 1 : -1) :\
   ((prio == PrioGhost) ? 0 : (prio == PrioBorder) ? 1 : \
      (prio == PrioMaster) ? 2 : -1))


/* map pointer to structure onto a pointer to its DDD_HDR */
#define PARHDR(obj)    (&((obj)->ddd))


#else   /* not ModelP */

/* define the number of prioities for objects */
#define ELEMENTPRIOS            1
#define NODEPRIOS                       1
#define VECTORPRIOS                     1
#define VERTEXPRIOS                     1

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
#define ID_PRT(x)   ((long)ID(x)),DDD_InfoGlobalId(PARHDR(x))
#define ID_FMTE     "%ld/%08x/%d"
#define ID_FFMTE    "%9ld/%08x/%02d"
#define ID_PRTE(x)  ((long)ID(x)),DDD_InfoGlobalId(PARHDR(x)),DDD_InfoPriority(PARHDR(x))

#define VID_FMT     ID_FMT
#define VID_FFMT    ID_FFMT
#define VID_PRT(x)  ((long)ID(x)),DDD_InfoGlobalId(PARHDRV(x))
#define VID_FMTE    ID_FMTE
#define VID_FFMTE   ID_FFMTE
#define VID_PRTE(x) ((long)ID(x)),DDD_InfoGlobalId(PARHDRV(x)),DDD_InfoPriority(PARHDRV(x))

#define EID_FMT     ID_FMT
#define EID_FFMT    ID_FFMT
#define EID_PRT(x)  ((long)ID(x)),DDD_InfoGlobalId(PARHDRE(x))
#define EID_FMTE    ID_FMTE
#define EID_FFMTE   ID_FFMTE
#define EID_PRTE(x) ((long)ID(x)),DDD_InfoGlobalId(PARHDRE(x)),DDD_InfoPriority(PARHDRE(x))

#define VINDEX_FMT     ID_FMT
#define VINDEX_FFMT    ID_FFMT
#define VINDEX_PRT(x)  ((long)VINDEX(x)),DDD_InfoGlobalId(PARHDR(x))
#define VINDEX_FMTE    ID_FMTE
#define VINDEX_FFMTE   ID_FFMTE
#define VINDEX_PRTE(x) ((long)VINDEX(x)),DDD_InfoGlobalId(PARHDR(x)),DDD_InfoPriority(PARHDR(x))

#else

#define ID_FMT      "%ld"
#define ID_FFMT     "%9ld"
#define ID_PRT(x)   ((long)ID(x))
#define ID_FMTE     "%ld"
#define ID_FFMTE    "%9ld"
#define ID_PRTE(x)  ID_PTR(x)

#define VID_FMT     ID_FMT
#define VID_FFMT    ID_FFMT
#define VID_PRT(x)  ID_PRT(x)
#define VID_FMTE    ID_FMTE
#define VID_FFMTE   ID_FFMTE
#define VID_PRTE(x) VID_PRT(x)

#define EID_FMT     ID_FMT
#define EID_FFMT    ID_FFMT
#define EID_PRT(x)  ID_PRT(x)
#define EID_FMTE    ID_FMTE
#define EID_FFMTE   ID_FFMTE
#define EID_PRTE(x) EID_PRT(x)

#define VINDEX_FMT     ID_FMT
#define VINDEX_FFMT    ID_FFMT
#define VINDEX_PRT(x)  ((long)VINDEX(x))
#define VINDEX_FMTE    ID_FMTE
#define VINDEX_FFMTE   ID_FFMTE
#define VINDEX_PRTE(x) VINDEX_PRT(x)

#endif


/****************************************************************************/
/*                                                                          */
/* exported global variables                                                */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP
/* DDD Interfaces */
extern DDD_IF ElementIF, ElementSymmIF;
extern DDD_IF BorderNodeIF, BorderNodeSymmIF, OuterNodeIF;
extern DDD_IF BorderVectorIF, BorderVectorSymmIF, OuterVectorIF;
#endif


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* functions implemented in parallel/dddif/support.c */
DOUBLE UG_GlobalSumDOUBLE (DOUBLE);
INT UG_GlobalMaxINT (INT);
void UG_GlobalSumNDOUBLE (INT, DOUBLE *);

#endif /* __PARGM_H__ */
