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

/* define dynamic lists */
#define ELEMENT_LIST     0
#define NODE_LIST        1
#define VECTOR_LIST      2
#define VERTEX_LIST      3

#ifdef ModelP

/* define the number of prioities for objects */
#define ELEMENTPRIOS            2
#define NODEPRIOS                       3
#define VECTORPRIOS                     3
#define VERTEXPRIOS                     2

/* define mapping from object priority to position in linked list */
#define PRIO2LISTPART(listtype,prio) \
  ((listtype == ELEMENT_LIST) ? ((prio == PrioGhost) ? 0 : (prio == PrioMaster) ? 1 : -1) :\
   ((prio == PrioGhost) ? 0 : (prio == PrioBorder) ? 1 : \
      (prio == PrioMaster) ? 2 : -1))

#else

/* define the number of prioities for objects */
#define ELEMENTPRIOS            1
#define NODEPRIOS                       1
#define VECTORPRIOS                     1
#define VERTEXPRIOS                     1

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

DOUBLE UG_GlobalSumDOUBLE (DOUBLE x);

#endif /* __PARGM_H__ */
