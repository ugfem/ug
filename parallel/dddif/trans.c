// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  trans.c                                                                                                       */
/*																			*/
/* Purpose:   create new grid distribution according to lb-marks of master  */
/*            elements.                                                                                                 */
/*																			*/
/* Author:	  Klaus Birken, Stefan Lang                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   961216 kb  begin                                                                          */
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

#include <assert.h>

#include "parallel.h"
#include "evm.h"
#include "general.h"
#include "ugm.h"
#include "debug.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

enum GhostCmds { GC_Keep, GC_ToMaster, GC_Delete };


#define XferElement(elem,dest,prio) \
  { PRINTDEBUG(dddif,1,("%4d: XferElement(): XferCopy elem=%08x dest=%d prio=%d\n", \
                        me,DDD_InfoGlobalId(PARHDRE(elem)), dest, prio)); \
    DDD_XferCopyObjX(PARHDRE(elem), dest, prio, \
                     (OBJT(elem)==BEOBJ) ?   \
                     BND_SIZE_TAG(TAG(elem)) :   \
                     INNER_SIZE_TAG(TAG(elem)));  }


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

static int Gather_ElemDest (DDD_OBJ obj, void *data)
{
  ELEMENT *e = (ELEMENT *)obj;

  *(DDD_PROC *)data = PARTITION(e);
}

static int Scatter_ElemDest (DDD_OBJ obj, void *data)
{
  ELEMENT *e = (ELEMENT *)obj;

  PARTITION(e) = *(DDD_PROC *)data;
}


static int UpdateGhostDests (MULTIGRID *theMG)
{
  DDD_IFOneway(ElementIF, IF_FORWARD, sizeof(DDD_PROC),
               Gather_ElemDest, Scatter_ElemDest);

  return 0;
}



/****************************************************************************/


static int Gather_GhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *elem = (ELEMENT *)obj;

  if (PARTITION(elem) == proc)
  {
    *((int *)data) = GC_ToMaster;
  }
  else
  {
    int j;

    *((int *)data) = GC_Delete;

    for(j=0; j<SIDES_OF_ELEM(elem); j++)
    {
      ELEMENT *nb = NBELEM(elem,j);

      if (nb!=NULL)
      {
        if (PARTITION(nb)==proc)
        {
          *((int *)data) = GC_Keep;
          break;
        }
      }
    }
  }

  return 0;
}


static int Scatter_GhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *elem = (ELEMENT *)obj;

  switch (*(int *)data)
  {
  case GC_Keep :
    /* do nothing */
    break;

  case GC_ToMaster :
    /* not needed anymore. kb 9070108 */
    /*DDD_PrioritySet(PARHDRE(elem), PrioMaster);*/
    break;

  case GC_Delete :
    DDD_XferDeleteObj(PARHDRE(elem));
    break;

  default :
    assert(1);
  }

  return 0;
}


static int SendGhostCmds (MULTIGRID *theMG)
{
  DDD_IFOnewayX(ElementIF, IF_FORWARD, sizeof(int),
                Gather_GhostCmd, Scatter_GhostCmd);

  return 0;
}



/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* Function:  XferGridWithOverlap                                           */
/*                                                                          */
/* Purpose:   send elements to other procs, keep overlapping region of one  */
/*            element, maintain correct priorities at interfaces.           */
/*                                                                          */
/*            the destination procs have been computed by theRCB function   */
/*            and put into the elements' PARTITION-entries.                 */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void XferGridWithOverlap (GRID *theGrid)
{
  ELEMENT *elem;
  NODE *node;


  for(elem=FIRSTELEMENT(theGrid); elem!=NULL; elem=SUCCE(elem))
  {
    int has_local_nb = FALSE;
    int j;

    /* create Master copy */
    XferElement(elem, PARTITION(elem), PrioMaster);



    /* create 1-overlapping of elements */
    for(j=0; j<SIDES_OF_ELEM(elem); j++)
    {
      ELEMENT *nb = NBELEM(elem,j);

      if (nb!=NULL)
      {
        if (PARTITION(elem)!=PARTITION(nb))
        {
          /* create Ghost copy */
          XferElement(elem, PARTITION(nb), PrioGhost);
        }

        /* remember any local neighbour element */
        if (PARTITION(nb)==me)
          has_local_nb = TRUE;
      }
    }


    /* consider elements on master-proc */
    if (PARTITION(elem)!=me)
    {
      if (!has_local_nb)
      {
        /* element isn't needed */
        PRINTDEBUG(dddif,1,("%d: XferGridWithOverlap(): XferDel elem=%d to p=%d prio=%d\n",
                            me,DDD_InfoGlobalId(PARHDRE(elem)),PARTITION(elem),PrioGhost));
        DDD_XferDeleteObj(PARHDRE(elem));
      }
    }
  }
}



/****************************************************************************/

static void InheritPartitionBottomTop (ELEMENT *e)
{
  int i;
  ELEMENT *SonList[MAX_SONS];

  if (GetSons(e,SonList) != GM_OK) assert(0);

  for(i=0; i<SONS_OF_ELEM(e); i++)
  {
    ELEMENT *son = SonList[i];
    if (son==NULL) break;

    PARTITION(son) = PARTITION(e);
    InheritPartitionBottomTop(son);
  }
}


/****************************************************************************/


int TransferGridFromCoarse (MULTIGRID *theMG)
{
  HEAP *theHeap = theMG->theHeap;
  GRID *theGrid = GRID_ON_LEVEL(theMG,0);       /* balance coarse grid */
  ELEMENT *e;
  int i, son;
  int g;


  /* send son elements to destination of father element */
  for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
  {
    InheritPartitionBottomTop(e);
  }

  /* send new destination to ghost elements */
  UpdateGhostDests(theMG);


  /* init transfer */
  ddd_HandlerInit(HSET_XFER);

  /* start physical transfer */
  DDD_XferBegin();

  {
    /* send 'commands' to ghosts in old partitioning */
    SendGhostCmds(theMG);

    /* send all grids */
    for(g=TOPLEVEL(theMG); g>=0; g--)
    {
      GRID *grid = GRID_ON_LEVEL(theMG,g);
      if (NT(grid)>0) XferGridWithOverlap(grid);
    }
  }

  DDD_XferEnd();



  /* remove all connections for vectors with PrioGhost */
  if (0)
  {
    int g;
    for(g=TOPLEVEL(theMG); g>=0; g--)
    {
      GRID *grid = GRID_ON_LEVEL(theMG,g);
      VECTOR *vec;

      for(vec=PRIO_LASTVECTOR(grid,PrioGhost); vec!=NULL; vec=PREDVC(vec))
      {
        DisposeConnectionFromVector(grid, vec);
      }
    }
  }


  /* set priorities of border nodes */
  /* TODO this is an extra communication. eventually integrate this
              with grid distribution phase. */
  {
    int g;
    for(g=TOPLEVEL(theMG); g>=0; g--)
    {
      GRID *grid = GRID_ON_LEVEL(theMG,g);
      dddif_SetOverlapPriorities(grid);
    }
  }


  DDD_ConsCheck();

  return 0;
}


/****************************************************************************/

#endif  /* ModelP */
