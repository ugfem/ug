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
#include "gm.h"
#include "evm.h"
#include "general.h"
#include "ugm.h"
#include "algebra.h"
#include "debug.h"
#include "devices.h"

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
  { PRINTDEBUG(dddif,1,("%4d: XferElement(): XferCopy elem=" EID_FMTX " dest=%d prio=%d\n", \
                        me,EID_PRTX(elem), dest, prio)); \
    XFERECOPYX((elem), dest, prio, \
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
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/

void AMGAgglomerate(MULTIGRID *theMG)
{
  INT level,Size;
  GRID    *theGrid;
  VECTOR  *theVector;

  level = BOTTOMLEVEL(theMG);
  if (level >= 0)
  {
    UserWriteF("AMGAgglomerate(): no amg level found, current bottom level is %d\n", level);
    return;
  }
  theGrid = GRID_ON_LEVEL(theMG,level);

  DDD_XferBegin();
  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    Size = sizeof(VECTOR)-sizeof(DOUBLE)
           +FMT_S_VEC_TP(MGFORMAT(theMG),VTYPE(theVector));
    XFERCOPYX(theVector,master,PrioMaster,Size);
    SETPRIO(theVector,PrioVGhost);
  }
  DDD_XferEnd();

  return;
}

/****************************************************************************/

static int Gather_ElemDest (DDD_OBJ obj, void *data)
{
  ELEMENT *theElement = (ELEMENT *)obj;

  *(DDD_PROC *)data = PARTITION(theElement);
}

static int Scatter_ElemDest (DDD_OBJ obj, void *data)
{
  ELEMENT *theElement = (ELEMENT *)obj;

  PARTITION(theElement) = *(DDD_PROC *)data;
}


static int UpdateGhostDests (MULTIGRID *theMG)
{
  DDD_IFOneway(ElementIF, IF_FORWARD, sizeof(DDD_PROC),
               Gather_ElemDest, Scatter_ElemDest);

  DDD_IFOneway(ElementVIF, IF_FORWARD, sizeof(DDD_PROC),
               Gather_ElemDest, Scatter_ElemDest);

  return 0;
}



/****************************************************************************/


static int Gather_GhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *theElement = (ELEMENT *)obj;
  INT j;

  /* TODO: not needed anymore. kb 9070108 */
  if (PARTITION(theElement) == proc)
  {
    *((int *)data) = GC_ToMaster;

    return(0);
  }

  if (PARTITION(theElement) != proc)
  {
    *((int *)data) = GC_Delete;

    for(j=0; j<SIDES_OF_ELEM(theElement); j++)
    {
      ELEMENT *nb = NBELEM(theElement,j);

      if (nb!=NULL)
      {
        if (PARTITION(nb)==proc)
        {
          *((int *)data) = GC_Keep;
          break;
        }
      }
    }
    return(0);
  }

  return(1);
}


static int Scatter_GhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *theElement = (ELEMENT *)obj;
  ELEMENT *SonList[MAX_SONS];
  INT i;

  switch (*(int *)data)
  {
  case GC_Keep :
    /* do nothing */
    break;

  case GC_ToMaster :
    /* TODO: not needed anymore. kb 9070108
                            SETEPRIO(theElement, PrioMaster);*/
    break;

  case GC_Delete :
    if (NSONS(theElement) > 0)
    {
      if (GetAllSons(theElement,SonList) != 0) return(1);
      i = 0;
      while (SonList[i] != NULL)
      {
        if (PARTITION(SonList[i]) == me) return(0);
        i++;
      }
    }
    XFEREDELETE(theElement);
    break;

  default :
    assert(1);
  }

  return(0);
}

static int Gather_VHGhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *theElement = (ELEMENT *)obj;
  ELEMENT *theFather      = EFATHER(theElement);
  ELEMENT *theNeighbor;
  INT j;

  if (PARTITION(theElement) != proc)
  {
    *((int *)data) = GC_Delete;

    for(j=0; j<SIDES_OF_ELEM(theElement); j++)
    {
      theNeighbor = NBELEM(theElement,j);

      if (theNeighbor != NULL)
      {
        if (PARTITION(theNeighbor) == proc)
        {
          *((int *)data) = GC_Keep;
          return (0);
        }
      }
    }

    if (LEVEL(theElement) > 0)
    {
      ASSERT(theFather != NULL);

      if (PARTITION(theFather) == proc)
      {
        *((int *)data) = GC_Keep;
        return (0);
      }
    }
    return (0);
  }

  return(1);
}

static int Scatter_VHGhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *theElement = (ELEMENT *)obj;
  ELEMENT *SonList[MAX_SONS];
  INT i;

  /* if element is needed after transfer here */
  if ((*(int *)data) == GC_Keep) return(0);

  /* element becomes master copy */
  if (PARTITION(theElement) == me) return(0);

  /* if a son resides as master keep element as vghost */
  if (GetAllSons(theElement,SonList) != GM_OK) return(0);
  i = 0;
  while (SonList[i] != NULL)
  {
    if (PARTITION(SonList[i]) == me) return(0);
    i++;
  }

  /* element is not needed on me any more */
  if ((*(int *)data) == GC_Delete)
  {
    XFEREDELETE(theElement);
    return(0);
  }

  return(1);
}

static int ComputeGhostCmds (MULTIGRID *theMG)
{
  DDD_IFOnewayX(ElementVHIF, IF_FORWARD, sizeof(int),
                Gather_VHGhostCmd, Scatter_VHGhostCmd);

  return(0);
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
  ELEMENT *theElement, *theFather, *theNeighbor;
  ELEMENT *SonList[MAX_SONS];
  NODE    *theNode;
  INT i,j,overlap_elem;

  for(theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    overlap_elem = 0;

    /* create Master copy */
    XferElement(theElement, PARTITION(theElement), PrioMaster);

    /* create 1-overlapping of horizontal elements */
    for(j=0; j<SIDES_OF_ELEM(theElement); j++)
    {
      theNeighbor = NBELEM(theElement,j);

      if (theNeighbor != NULL)
      {
        if (PARTITION(theElement)!=PARTITION(theNeighbor))
        {
          /* create Ghost copy */
          XferElement(theElement, PARTITION(theNeighbor), PrioGhost);
        }

        /* remember any local neighbour element */
        if (PARTITION(theNeighbor)==me)
          overlap_elem = 1;
      }
    }

    /* create 1-overlapping of vertical elements */
    theFather = EFATHER(theElement);
    if (theFather != NULL)
    {
      if (PARTITION(theFather) != PARTITION(theElement) ||
          !EMASTER(theFather))
        /* create VGhost copy */
        XferElement(theFather, PARTITION(theElement), PrioVGhost);
    }
    else
    {
      ASSERT(LEVEL(theElement) == 0);
    }

    /* consider elements on master-proc */
    if (PARTITION(theElement)!=me)
    {
      if (NSONS(theElement) > 0)
      {
        if (GetAllSons(theElement,SonList) != 0) return;
        i = 0;
        while (SonList[i] != NULL)
        {
          if (PARTITION(SonList[i]) == me)
          {
            overlap_elem += 2;
            break;
          }
          i++;
        }
      }

      PRINTDEBUG(dddif,1,("%d: XferGridWithOverlap(): elem=" EID_FMTX " p=%d new prio=%d\n",
                          me,EGID(theElement),PARTITION(theElement),overlap_elem));

      if (overlap_elem > 0)
      {
        /* element is needed, set new prio */
        switch (overlap_elem)
        {
        case (1) :
          SETEPRIO(theElement,PrioGhost);
          break;

        case (2) :
          SETEPRIO(theElement,PrioVGhost);
          break;

        case (3) :
          SETEPRIO(theElement,PrioVGhost);
          break;

        default :
          assert(0);
        }
      }
      else
      {
        /* element isn't needed */
        PRINTDEBUG(dddif,2,("%d: XferGridWithOverlap(): XferDel elem=%d to p=%d prio=%d\n",
                            me,EGID(theElement),PARTITION(theElement),PrioGhost));
        XFEREDELETE(theElement);
      }
    }
  }

  /* set prios on ghost objects */
  SetGhostObjectPriorities(theGrid);
}



/****************************************************************************/

static void InheritPartitionBottomTop (ELEMENT *e)
{
  int i;
  ELEMENT *SonList[MAX_SONS];

  if (GetSons(e,SonList) != GM_OK) assert(0);

  for(i=0; i<MAX_SONS; i++)
  {
    ELEMENT *son = SonList[i];
    if (son==NULL) break;

    PARTITION(son) = PARTITION(e);
    InheritPartitionBottomTop(son);
  }
}


/****************************************************************************/


int TransferGridFromLevel (MULTIGRID *theMG, INT level)
{
  HEAP *theHeap = theMG->theHeap;
  GRID *theGrid = GRID_ON_LEVEL(theMG,level);       /* transfer grid starting at*/
  ELEMENT *e;
  int i, son;
  int g;

  /* dispose negative levels */
  if (level < 1)
    if (DisposeAMGLevels(theMG) != 0)
      return 1;

  /* send new destination to ghost elements */
  UpdateGhostDests(theMG);


  /* init transfer */
  ddd_HandlerInit(HSET_XFER);

  /* start physical transfer */
  DDD_XferBegin();

  {
    /* send 'commands' to ghosts in old partitioning */
    ComputeGhostCmds(theMG);

    /* send all grids */
    for(g=TOPLEVEL(theMG); g>=0; g--)
    {
      GRID *grid = GRID_ON_LEVEL(theMG,g);
      if (NT(grid)>0) XferGridWithOverlap(grid);
    }
  }

  DDD_XferEnd();

  /* set priorities of border nodes */
  /* TODO this is an extra communication. eventually integrate this
              with grid distribution phase. */
  {
    int g;
    for(g=TOPLEVEL(theMG); g>=0; g--)
    {
      GRID *grid = GRID_ON_LEVEL(theMG,g);
      ConstructConsistentGrid(grid);
    }
  }

    #ifndef __EXCHANGE_CONNECTIONS__
  MGCreateConnection(theMG);
        #endif

  DDD_ConsCheck();

  return 0;
}

int TransferGrid (MULTIGRID *theMG)
{
  TransferGridFromLevel(theMG,0);
}

/****************************************************************************/

#endif  /* ModelP */
