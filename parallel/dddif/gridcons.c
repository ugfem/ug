// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gridcons.c													*/
/*																			*/
/* Purpose:   basic functions for managing consistency of distributed grids */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: birken@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   960906 kb  begin                                                                                          */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>

#include "debug.h"
#include "parallel.h"
#include "general.h"
#include "gm.h"


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* macros for setting object priorities with related objects */
#define NODE_PRIORITY_SET(theNode,prio)                                      \
  /* set priorities of node */                                     \
  DDD_PrioritySet(PARHDR(theNode),prio);                           \
                                                                             \
  /* set priorities of vertex and node vector */                   \
  DDD_PrioritySet(PARHDRV(MYVERTEX(theNode)),prio);                \
                                                                             \
  if (TYPE_DEF_IN_GRID(theGrid,NODEVECTOR))                        \
    DDD_PrioritySet(PARHDR(NVECTOR(theNode)),prio);

#ifdef __TWODIM__
#define PRIO_SET_EDGE(theEdge,PrioGhost)
#endif
#ifdef __THREEDIM__
#define PRIO_SET_EDGE(theEdge,prio)                                          \
  DDD_PrioritySet(PARHDR(theEdge),prio);
#endif

#define EDGE_PRIORITY_SET(theEdge,prio)                                      \
  /* set priorities of node for 3D */                              \
  PRIO_SET_EDGE(theEdge,PrioGhost)                                 \
                                                                             \
  /* set priority of edge vector */                                \
  if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR))                        \
    DDD_PrioritySet(PARHDR(EDVECTOR(theEdge)),PrioGhost);


/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/*
        for all PrioMaster-nodes with remote copies, set exactly one
        to PrioMaster, the other copies to PrioBorder in order to establish
        the BorderNodeIF. this is done for one grid.
 */


static int dddif_ComputeNodeBorderPrios (DDD_OBJ obj)
{
  NODE    *node  = (NODE *)obj;
  int     *plist = DDD_InfoProcList(PARHDR(node));
  int i, min_proc = procs;

  /*
          minimum processor number will get Master-node,
          all others get Border-nodes
   */
  for(i=0; plist[i]>=0; i+=2)
  {
    if (plist[i+1]==PrioMaster && plist[i]<min_proc)
      min_proc = plist[i];
  }

  if (min_proc==procs)
    return;

  if (me!=min_proc)
    DDD_PrioritySet(PARHDR(node), PrioBorder);
}

static int dddif_ComputeVectorBorderPrios (DDD_OBJ obj)
{
  VECTOR  *vector  = (VECTOR *)obj;
  int     *plist = DDD_InfoProcList(PARHDR(vector));
  int i, min_proc = procs;

  /*
          minimum processor number will get Master-node,
          all others get Border-nodes
   */
  for(i=0; plist[i]>=0; i+=2)
  {
    if (plist[i+1]==PrioMaster && plist[i]<min_proc)
      min_proc = plist[i];
  }

  if (min_proc==procs)
    return;

  if (me!=min_proc)
    DDD_PrioritySet(PARHDR(vector), PrioBorder);
}


void dddif_SetGhostObjectPriorities (GRID *theGrid)
{
  ELEMENT *theElement;
  NODE    *theNode;
  INT i,prio,*proclist;

  /* reset USED flag for nodes of ghostelements */
  for (theElement=PRIO_LASTELEMENT(theGrid,PrioGhost);
       theElement!=NULL;
       theElement=PREDE(theElement))
  {
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      SETUSED(theNode,0);
    }
  }

  /* set USED flag for nodes of master elements */
  for (theElement=PRIO_FIRSTELEMENT(theGrid,PrioMaster);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      SETUSED(theNode,1);
    }
  }

  /* set node priorities for ghostelements */
  for (theElement=PRIO_LASTELEMENT(theGrid,PrioGhost);
       theElement!=NULL;
       theElement=PREDE(theElement))
  {
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);

      prio = DDD_InfoPriority(PARHDR(theNode));

      /* check if its a master node */
      if (prio != PrioGhost)
      {
        /* is it a ghost node */
        if (USED(theNode)==0)
        {
          PRINTDEBUG(dddif,3,(PFMT " dddif_SetGhostObjectPriorities():"
                              " downgrade n=" ID_FMTX " from=%d to PrioGhost\n",
                              me,ID_PRTX(theNode),prio));

          /* set node priorities of node to ghost */
          NODE_PRIORITY_SET(theNode,PrioGhost)
        }
      }
    }


    if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR) || DIM==3)
    {
      /* set edge priorities */
      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {
        NODE *theNode0, *theNode1;

        theNode0 = CORNER(theElement,CORNER_OF_EDGE(theElement,i,0));
        theNode1 = CORNER(theElement,CORNER_OF_EDGE(theElement,i,1));

        /* is one of the nodes a ghost node, */
        /* then this is a ghost edge         */
        if (USED(theNode0)==0 || USED(theNode1)==0)
        {
          EDGE *theEdge;

          theEdge = GetEdge(theNode0,theNode1);
          ASSERT(theEdge != NULL);

          EDGE_PRIORITY_SET(theEdge,PrioGhost);
        }
      }

                        #ifdef __THREEDIM__
      /* if one of the side nodes is a ghost node */
      /* then its a ghost side vector             */
      if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
        {
          INT j;

          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            if (USED(theNode) == 0)
            {
              /* set side vector priority */
              DDD_PrioritySet(PARHDR(SVECTOR(theElement,i)),
                              PrioGhost);
              break;
            }
          }
        }
                        #endif
    }
  }
}


void dddif_SetOverlapPriorities (GRID *theGrid)
{
  DDD_XferBegin();

  DDD_IFAExecLocal(BorderNodeSymmIF, GLEVEL(theGrid),
                   dddif_ComputeNodeBorderPrios);

  DDD_IFAExecLocal(BorderVectorSymmIF, GLEVEL(theGrid),
                   dddif_ComputeVectorBorderPrios);

  dddif_SetGhostObjectPriorities(theGrid);

  DDD_XferEnd();
}


INT CheckInterfaces(GRID *theGrid)
{
  /* check ddd interface consistency */
  DDD_ConsCheck();

  return(GM_OK);
}

/****************************************************************************/

#endif
