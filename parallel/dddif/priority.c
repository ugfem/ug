// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  priority.c													*/
/*																			*/
/* Purpose:   functions for managing priorities of distributed objects      */
/*																			*/
/* Author:	  Stefan Lang                                                                       */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: birken@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   980201 sl  begin                                                                                          */
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
#include "refine.h"
#include "ugm.h"
#include "evm.h"
#include "shapes.h"
#include "ugdevices.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* macros for merge new priority with objects existing one */
/* valid only for all types of ghost priorities            */
#define PRIO_CALC(e) ((USED(e) && THEFLAG(e)) ? PrioVHGhost :                \
                      (THEFLAG(e)) ? PrioVGhost : (USED(e)) ?              \
                      PrioHGhost : (assert(0),0))

/* macros for setting object priorities with related objects */
/* macros for setting object priorities with related objects */
#define NODE_PRIORITY_SET(g,n,prio)                                          \
  {                                                                    \
    /* set priorities of node */                                     \
    SETPRIOX(n,prio);                                                \
                                                                             \
    if (VEC_DEF_IN_OBJ_OF_GRID(g,NODEVEC))                           \
      if (NVECTOR(n) != NULL)                                      \
        SETPRIOX(NVECTOR(n),prio);                               \
  }

#ifdef __TWODIM__
#define PRIO_SET_EDGE(e,prio)
#endif
#ifdef __THREEDIM__
#define PRIO_SET_EDGE(e,prio)  SETPRIOX(e,prio);
#endif

#define EDGE_PRIORITY_SET(g,e,prio)                                          \
  {                                                                    \
    /* set priorities of node for 3D */                              \
    PRIO_SET_EDGE(e,prio)                                            \
                                                                             \
    /* set priority of edge vector */                                \
    if (VEC_DEF_IN_OBJ_OF_GRID(g,EDGEVEC))                           \
      if (EDVECTOR(e) != NULL)                                     \
        SETPRIOX(EDVECTOR(e),prio);                              \
  }


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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


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


/****************************************************************************/
/*
   ComputeNodeBorderPrios -

   SYNOPSIS:
   static int ComputeNodeBorderPrios (DDD_OBJ obj);

   PARAMETERS:
   .  obj -

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int ComputeNodeBorderPrios (DDD_OBJ obj)
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
    return(0);

  if (me!=min_proc)
    SETPRIO(node, PrioBorder);
}


/****************************************************************************/
/*
   ComputeVectorBorderPrios -

   SYNOPSIS:
   static int ComputeVectorBorderPrios (DDD_OBJ obj);

   PARAMETERS:
   .  obj -

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int ComputeVectorBorderPrios (DDD_OBJ obj)
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
    return(0);

  if (me!=min_proc)
    SETPRIO(vector, PrioBorder);
}

#ifdef __THREEDIM__


/****************************************************************************/
/*
   ComputeEdgeBorderPrios -

   SYNOPSIS:
   static int ComputeEdgeBorderPrios (DDD_OBJ obj);

   PARAMETERS:
   .  obj -

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static int ComputeEdgeBorderPrios (DDD_OBJ obj)
{
  EDGE    *edge  =        (EDGE *)obj;
  int             *plist =        DDD_InfoProcList(PARHDR(edge));
  int i, min_proc     = procs;

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
    return(0);

  if (me!=min_proc)
    SETPRIO(edge, PrioBorder);
}
#endif


/****************************************************************************/
/*
   SetGhostObjectPriorities -

   SYNOPSIS:
   void SetGhostObjectPriorities (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void SetGhostObjectPriorities (GRID *theGrid)
{
  ELEMENT *theElement,*theNeighbor,*SonList[MAX_SONS];
  NODE    *theNode;
  EDGE    *theEdge;
  VECTOR  *theVector;
  INT i,prio,*proclist,hghost,vghost;

  /* reset USED flag for objects of ghostelements */
  for (theElement=PFIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    SETUSED(theElement,0); SETTHEFLAG(theElement,0);
    for (i=0; i<EDGES_OF_ELEM(theElement); i++)
    {
      theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                        CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
      ASSERT(theEdge != NULL);
      SETUSED(theEdge,0); SETTHEFLAG(theEdge,0);
    }
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        theVector = SVECTOR(theElement,i);
        if (theVector != NULL) {
          SETUSED(theVector,0);
          SETTHEFLAG(theVector,0);
        }
      }
  }
  /* to reset also nodes which are at corners of the boundary */
  /* reset of nodes need to be done through the node list     */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    SETUSED(theNode,0); SETTHEFLAG(theNode,0);
    SETMODIFIED(theNode,0);
  }

  /* set FLAG for objects of horizontal and vertical overlap */
  for (theElement=PFIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if (PARTITION(theElement) == me) continue;

    /* check for horizontal ghost */
    hghost = 0;
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      theNeighbor = NBELEM(theElement,i);
      if (theNeighbor == NULL) continue;

      if (PARTITION(theNeighbor) == me)
      {
        hghost = 1;
        break;
      }
    }

    /* check for vertical ghost */
    vghost = 0;
    GetAllSons(theElement,SonList);
    for (i=0; SonList[i]!=NULL; i++)
    {
      if (PARTITION(SonList[i]) == me)
      {
        vghost = 1;
        break;
      }
    }

    /* one or both of vghost and hghost should be true here   */
    /* except for elements which will be disposed during Xfer */

    if (vghost) SETTHEFLAG(theElement,1);
    if (hghost) SETUSED(theElement,1);
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      if (vghost) SETTHEFLAG(theNode,1);
      if (hghost) SETUSED(theNode,1);
    }
    for (i=0; i<EDGES_OF_ELEM(theElement); i++)
    {
      theEdge = GetEdge(CORNER_OF_EDGE_PTR(theElement,i,0),
                        CORNER_OF_EDGE_PTR(theElement,i,1));
      ASSERT(theEdge != NULL);
      if (vghost) SETTHEFLAG(theEdge,1);
      if (hghost) SETUSED(theEdge,1);
    }
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        theVector = SVECTOR(theElement,i);
        if (theVector != NULL) {
          if (vghost) SETTHEFLAG(theVector,1);
          if (hghost) SETUSED(theVector,1);
        }
      }
  }

  DEBUG_TIME(0);

  /* set USED flag for objects of master elements */
  /* reset FLAG for objects of master elements  */
  for (theElement=PFIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if (PARTITION(theElement) != me) continue;

    SETUSED(theElement,0); SETTHEFLAG(theElement,0);
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      SETUSED(theNode,0); SETTHEFLAG(theNode,0);
      SETMODIFIED(theNode,1);
    }
    for (i=0; i<EDGES_OF_ELEM(theElement); i++)
    {
      theEdge = GetEdge(CORNER_OF_EDGE_PTR(theElement,i,0),
                        CORNER_OF_EDGE_PTR(theElement,i,1));
      ASSERT(theEdge != NULL);
      SETUSED(theEdge,0); SETTHEFLAG(theEdge,0);
    }
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        theVector = SVECTOR(theElement,i);
        if (theVector != NULL) {
          SETUSED(theVector,0);
          SETTHEFLAG(theVector,0);
        }
      }
  }

  DEBUG_TIME(0);

  /* set object priorities for ghostelements */
  for (theElement=PFIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if (PARTITION(theElement) == me) continue;

    if (USED(theElement) || THEFLAG(theElement))
    {
      prio = PRIO_CALC(theElement);
      PRINTDEBUG(gm,1,("SetGhostObjectPriorities(): e=" EID_FMTX " new prio=%d\n",
                       EID_PRTX(theElement),prio))
      SETEPRIOX(theElement,prio);

      if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
      {
        theVector = EVECTOR(theElement);
        if (theVector != NULL)
          SETPRIOX(theVector,prio);
      }
    }

    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC) || DIM==3)
    {
      /* set edge priorities */
      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {

        theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                          CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
        ASSERT(theEdge != NULL);

        if (USED(theEdge) || THEFLAG(theEdge))
        {
          PRINTDEBUG(dddif,3,(PFMT " dddif_SetGhostObjectPriorities():"
                              " downgrade edge=" EDID_FMTX " from=%d to PrioHGhost\n",
                              me,EDID_PRTX(theEdge),prio));

          EDGE_PRIORITY_SET(theGrid,theEdge,PRIO_CALC(theEdge));
        }
        else
          EDGE_PRIORITY_SET(theGrid,theEdge,PrioMaster);
      }

                        #ifdef __THREEDIM__
      /* if one(all) of the side nodes is (are) a hghost (vghost) node   */
      /* then its a hghost (vghost) side vector                          */
      if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
        {
          if (USED(theVector) || THEFLAG(theVector))
            SETPRIOX(theVector,PRIO_CALC(theVector));
        }
                        #endif
    }
  }
  /* to set also nodes which are at corners of the boundary   */
  /* set them through the node list                           */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    /* check if its a master node */
    if (USED(theNode) || THEFLAG(theNode))
    {
      PRINTDEBUG(dddif,3,(PFMT " dddif_SetGhostObjectPriorities():"
                          " downgrade node=" ID_FMTX " from=%d to PrioHGhost\n",
                          me,ID_PRTX(theNode),prio));

      /* set node priorities of node to ghost */
      NODE_PRIORITY_SET(theGrid,theNode,PRIO_CALC(theNode))
    }
    else if (MODIFIED(theNode) == 0)
    {
      /* this is a node of the boundary without connection to master elements */
      NODE_PRIORITY_SET(theGrid,theNode,PrioHGhost)
    }
    /* this is needed only for consistency after refinement */
    /* ghost nodes which belong after refinement to master  */
    /* elements have to be upgraded explicitly (980126 s.l.)*/
    else if (MODIFIED(theNode) == 1)
    {
      NODE_PRIORITY_SET(theGrid,theNode,PrioMaster)
    }
  }
}


/****************************************************************************/
/*
   SetBorderPriorities -

   SYNOPSIS:
   INT SetBorderPriorities (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT SetBorderPriorities (GRID *theGrid)
{
  DDD_IFAExecLocal(BorderNodeSymmIF,GRID_ATTR(theGrid),
                   ComputeNodeBorderPrios);

  /* TODO: distinguish two cases:
     1. only nodevectors then setting of vector prios can
          be done in ComputeNodeBorderPrios without extra communiction
     2. with other vectortypes (side and/or edgevectors) use
          ComputeVectorBorderPrios
   */
  DDD_IFAExecLocal(BorderVectorSymmIF,GRID_ATTR(theGrid),
                   ComputeVectorBorderPrios);

#ifdef __THREEDIM__
  DDD_IFAExecLocal(BorderEdgeSymmIF,GRID_ATTR(theGrid),
                   ComputeEdgeBorderPrios);
#endif

  return(GM_OK);
}


/****************************************************************************/
/*
   SetGridBorderPriorities -

   SYNOPSIS:
   INT SetGridBorderPriorities (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/


INT SetGridBorderPriorities (GRID *theGrid)
{
  /* set border priorities on next higher level */
  if (SetBorderPriorities(UPGRID(theGrid)) != GM_OK) return(GM_FATAL);

  return(GM_OK);
}

#endif
