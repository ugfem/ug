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
#include "refine.h"
#include "algebra.h"
#include "debug.h"
#include "ugdevices.h"
#ifdef DYNAMIC_MEMORY_ALLOCMODEL
#include "mgheapmgr.h"
#endif
#ifdef __DLB__
#include "dlb.h"
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


/****************************************************************************/
/*
   AMGAgglomerate -

   SYNOPSIS:
   void AMGAgglomerate(MULTIGRID *theMG);

   PARAMETERS:
   .  theMG

   DESCRIPTION:

   RETURN VALUE:
   void
 */
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


/****************************************************************************/
/*
   Gather_ElemDest -

   SYNOPSIS:
   static int Gather_ElemDest (DDD_OBJ obj, void *data);

   PARAMETERS:
   .  obj
   .  data

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int Gather_ElemDest (DDD_OBJ obj, void *data)
{
  ELEMENT *theElement = (ELEMENT *)obj;

  *(DDD_PROC *)data = PARTITION(theElement);
}


/****************************************************************************/
/*
   Scatter_ElemDest -

   SYNOPSIS:
   static int Scatter_ElemDest (DDD_OBJ obj, void *data);

   PARAMETERS:
   .  obj
   .  data

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int Scatter_ElemDest (DDD_OBJ obj, void *data)
{
  ELEMENT *theElement = (ELEMENT *)obj;

  PARTITION(theElement) = *(DDD_PROC *)data;
}


/****************************************************************************/
/*
   UpdateGhostDests -

   SYNOPSIS:
   static int UpdateGhostDests (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int UpdateGhostDests (MULTIGRID *theMG)
{
  DDD_IFOneway(ElementIF, IF_FORWARD, sizeof(DDD_PROC),
               Gather_ElemDest, Scatter_ElemDest);

  DDD_IFOneway(ElementVIF, IF_FORWARD, sizeof(DDD_PROC),
               Gather_ElemDest, Scatter_ElemDest);

  return 0;
}



/****************************************************************************/


/****************************************************************************/
/*
   Gather_GhostCmd -

   SYNOPSIS:
   static int Gather_GhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

   PARAMETERS:
   .  obj
   .  data
   .  proc
   .  prio

   DESCRIPTION:

   RETURN VALUE:
   int
 */
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


/****************************************************************************/
/*
   Scatter_GhostCmd -

   SYNOPSIS:
   static int Scatter_GhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

   PARAMETERS:
   .  obj
   .  data
   .  proc
   .  prio

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int Scatter_GhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *theElement = (ELEMENT *)obj;
  ELEMENT *SonList[MAX_SONS];
  INT i,j;

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


/****************************************************************************/
/*
   Gather_VHGhostCmd -

   SYNOPSIS:
   static int Gather_VHGhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

   PARAMETERS:
   .  obj
   .  data
   .  proc
   .  prio

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

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

    /* wrong:		if (LEVEL(theElement) > 0)
                    {
                            ASSERT(theFather != NULL);

                            if (PARTITION(theFather) == proc)
                            {
       *((int *)data) = GC_Keep;
                                    return (0);
                            }
                    } */
    return (0);
  }

  return(1);
}


/****************************************************************************/
/*
   Scatter_VHGhostCmd -

   SYNOPSIS:
   static int Scatter_VHGhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

   PARAMETERS:
   .  obj
   .  data
   .  proc
   .  prio

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int Scatter_VHGhostCmd (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  ELEMENT *theElement = (ELEMENT *)obj;
  ELEMENT *SonList[MAX_SONS];
  INT i,j;

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


/****************************************************************************/
/*
   ComputeGhostCmds -

   SYNOPSIS:
   static int ComputeGhostCmds (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int ComputeGhostCmds (MULTIGRID *theMG)
{
  DDD_IFOnewayX(ElementVHIF, IF_FORWARD, sizeof(int),
                Gather_VHGhostCmd, Scatter_VHGhostCmd);

  return(0);
}

/****************************************************************************/

#ifdef __OVERLAP2__
static int XferNodesForOverlap2 (GRID *theGrid)
{
  ELEMENT *theElement;
  NODE    *theNode;
  INT i,part;
  MATRIX  *mat,*mat2;
  VECTOR  *vec,*dest;
  INT Size;

  for(theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    /* ensuring the overlap is for all elements necessary, because even for
       elements remaining on this processor their overlap may have been deleted
       some lines above */

    part = PARTITION(theElement);

    /* traverse all corner vectors */
    for(i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      vec = NVECTOR(theNode);

      PRINTDEBUG(dddif,2,(PFMT " XferGridWithOverlap():  e=" EID_FMTX
                          " Xfer n=" ID_FMTX " i=%d\n",
                          me,EID_PRTX(theElement),ID_PRTX(theNode),i))

      /* for master vectors all matrix neighbors within link depth 2
          must be copied too; the vec itself is automatically copied
              by the node-copy-handler */
      for(mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
      {
        dest = MDEST(mat);
        if (dest != NULL)
        {
          Size = sizeof(VECTOR)-sizeof(DOUBLE)
                 +FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(dest));

          PRINTDEBUG(dddif,2,(PFMT " XferGridWithOverlap(): n=" ID_FMTX
                              " Xfer NODEVEC=" VINDEX_FMTX " 1. NBvec  "
                              VINDEX_FMTX " size=%d\n",
                              me,ID_PRTX(theNode),VINDEX_PRTX(vec),VINDEX_PRTX(dest),Size))
          /* TODO: only vectors are necessary; only for debugging: send also the corresponding node to have geometric information */
          SETNO_DELETE_OVERLAP2((NODE*)VOBJECT(dest),1);
          DDD_XferCopyObj(PARHDR((NODE*)VOBJECT(dest)), part, PrioHGhost);
          /*DDD_XferCopyObjX(PARHDR(dest), part, PrioHGhost, Size);*/
        }

        for(mat2=MNEXT(VSTART(dest)); mat2!=NULL; mat2=MNEXT(mat2))
        {
          dest = MDEST(mat2);
          if (dest != NULL)
          {
            Size = sizeof(VECTOR)-sizeof(DOUBLE)
                   +FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(dest));

            PRINTDEBUG(dddif,2,(PFMT " XferGridWithOverlap(): n=" ID_FMTX
                                " Xfer NODEVEC=" VINDEX_FMTX " 2. NBvec  "
                                VINDEX_FMTX " size=%d\n",
                                me,ID_PRTX(theNode),VINDEX_PRTX(vec),VINDEX_PRTX(dest),Size))

            /* TODO: only vectors are necessary; only for debugging: send also the corresponding node to have geometric information */
            SETNO_DELETE_OVERLAP2((NODE*)VOBJECT(dest),1);
            DDD_XferCopyObj(PARHDR((NODE*)VOBJECT(dest)), part, PrioHGhost);
            /*DDD_XferCopyObjX(PARHDR(dest), part, PrioHGhost, Size);*/
          }
        }
      }
    }
  }
  return (0);
}
#endif

/****************************************************************************/
/*
   XferGridWithOverlap - send elements to other procs, keep overlapping region of one element, maintain correct priorities at interfaces.

   SYNOPSIS:
   static void XferGridWithOverlap (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:
   This function sends elements to other procs, keeps overlapping region of one element and maintains correct priorities at interfaces. The destination procs have been computed by theRCB function and put into the elements' PARTITION-entries.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static int XferGridWithOverlap (GRID *theGrid)
{
  ELEMENT *theElement, *theFather, *theNeighbor;
  ELEMENT *SonList[MAX_SONS];
  NODE    *theNode;
  INT i,j,overlap_elem,part;
  INT migrated = 0;

  for(theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    /* goal processor */
    part = PARTITION(theElement);

    /* create Master copy */
    XferElement(theElement, part, PrioMaster);

                #ifdef STAT_OUT
    /* count elems */
    if (part != me) migrated++;
                #endif
  }

  /* create grid overlap */
  for(theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    overlap_elem = 0;

    /* create 1-overlapping of horizontal elements */
    for(j=0; j<SIDES_OF_ELEM(theElement); j++)
    {
      theNeighbor = NBELEM(theElement,j);

      if (theNeighbor != NULL)
      {
        if (PARTITION(theElement)!=PARTITION(theNeighbor))
        {
          /* create Ghost copy */
          XferElement(theElement, PARTITION(theNeighbor), PrioHGhost);
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
          SETEPRIO(theElement,PrioHGhost);
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
                            me,EGID(theElement),PARTITION(theElement),PrioHGhost));

        XFEREDELETE(theElement);
      }
    }
  }

#ifdef __OVERLAP2__
  if (XferNodesForOverlap2(theGrid)) assert(0);
#endif

  return(migrated);
}



/****************************************************************************/


/****************************************************************************/
/*
   InheritPartitionBottomTop -

   SYNOPSIS:
   static void InheritPartitionBottomTop (ELEMENT *e);

   PARAMETERS:
   .  e

   DESCRIPTION:

   RETURN VALUE:
   void
 */
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


/****************************************************************************/
/*
   TransferGridFromLevel -

   SYNOPSIS:
   int TransferGridFromLevel (MULTIGRID *theMG, INT level);

   PARAMETERS:
   .  theMG
   .  level

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

int TransferGridFromLevel (MULTIGRID *theMG, INT level)
{
  GRID *theGrid = GRID_ON_LEVEL(theMG,level);       /* transfer grid starting at*/
  INT g,alreadydisposed;
  INT migrated = 0;       /* number of elements moved */
  DOUBLE trans_begin, trans_end, cons_end;

#ifdef __OVERLAP2__
  NODE *node;

  ASSERT(AllocateControlEntry(NODE_CW,NO_DELETE_OVERLAP2_LEN,&ce_NO_DELETE_OVERLAP2) == GM_OK);
  for (g=0; g<=TOPLEVEL(theMG); g++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,g);
    for( node=PFIRSTNODE(theGrid); node!= NULL; node=SUCCN(node) )
      SETNO_DELETE_OVERLAP2(node,0);                    /* reset flag */
  }
#endif

        #ifdef DYNAMIC_MEMORY_ALLOCMODEL
  if (theMG->bottomtmpmem == 1)
  {
    if (DisposeBottomHeapTmpMemory(theMG)) REP_ERR_RETURN(1);
    alreadydisposed = 0;
  }
  else
    alreadydisposed = 1;
        #endif

  trans_begin = CURRENT_TIME;

  /* dispose negative levels */
  if (level < 1)
    if (DisposeAMGLevels(theMG) != 0)
      return 1;

    #ifdef __DLB__
  /* optimize mapping */
  {
    DLB_CONFIG config;

    /* configuration of balancing */
    if (DLB_ConfigGet(&config)) REP_ERR_RETURN(1);
    if (config.map == 1) DLB_Mapping(theMG);
  }
    #endif

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
    /*		for (g=TOPLEVEL(theMG); g>=0; g--) */
    for (g=0; g<=TOPLEVEL(theMG); g++)
    {
      GRID *theGrid = GRID_ON_LEVEL(theMG,g);
      if (NT(theGrid)>0) migrated += XferGridWithOverlap(theGrid);
    }
  }

  DDD_XferEnd();

  trans_end = CURRENT_TIME;

#ifdef __OVERLAP2__
  for (g=0; g<=TOPLEVEL(theMG); g++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,g);
    for( node=PFIRSTNODE(theGrid); node!= NULL; node=SUCCN(node) )
      SETNO_DELETE_OVERLAP2(node,0);                    /* reset flag */
  }
  FreeControlEntry(ce_NO_DELETE_OVERLAP2);
  ce_NO_DELETE_OVERLAP2 = -1;           /* don't use further NO_DELETE_OVERLAP2 */
#endif

  /* set priorities of border nodes */
  /* TODO this is an extra communication. eventually integrate this
              with grid distribution phase. */
        #ifdef NEW_GRIDCONS_STYLE
  ConstructConsistentMultiGrid(theMG);
        #else
  {
    for(g=0; g<=TOPLEVEL(theMG); g++)
    {
      GRID *theGrid = GRID_ON_LEVEL(theMG,g);
      ConstructConsistentGrid(theGrid);
    }
  }
        #endif

    #ifndef __EXCHANGE_CONNECTIONS__
        #ifdef DYNAMIC_MEMORY_ALLOCMODEL
  if (alreadydisposed==0)
        #endif
  MGCreateConnection(theMG);
        #endif

  /* the grid has changed at least on one processor, thus reset MGSTATUS on all processors */
  RESETMGSTATUS(theMG);

  cons_end = CURRENT_TIME;

        #ifdef STAT_OUT
  /* sum up moved elements */
  migrated = UG_GlobalSumINT(migrated);

  UserWriteF("MIGRATION: migrated=%d t_migrate=%.2f t_cons=%.2f\n",
             migrated,trans_end-trans_begin,cons_end-trans_end);
        #endif

        #ifdef Debug
  if (0)
    DDD_ConsCheck();
        #endif

  return 0;
}


/****************************************************************************/
/*
   TransferGrid -

   SYNOPSIS:
   int TransferGrid (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

int TransferGrid (MULTIGRID *theMG)
{
  TransferGridFromLevel(theMG,0);
}

/****************************************************************************/

#endif  /* ModelP */
