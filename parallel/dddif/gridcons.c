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
#include "evm.h"
#include "shapes.h"

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
#define NODE_PRIORITY_SET(g,n,prio)                                          \
  /* set priorities of node */                                     \
  SETPRIO(n,prio);                                                 \
                                                                             \
  if (TYPE_DEF_IN_GRID(g,NODEVECTOR))                              \
    SETPRIO(NVECTOR(n),prio);

#ifdef __TWODIM__
#define PRIO_SET_EDGE(e,prio)
#endif
#ifdef __THREEDIM__
#define PRIO_SET_EDGE(e,prio)  SETPRIO(e,prio);
#endif

#define EDGE_PRIORITY_SET(g,e,prio)                                          \
  /* set priorities of node for 3D */                              \
  PRIO_SET_EDGE(e,prio)                                            \
                                                                             \
  /* set priority of edge vector */                                \
  if (TYPE_DEF_IN_GRID(g,EDGEVECTOR))                              \
    SETPRIO(EDVECTOR(e),prio);

#define CHECK_OBJECT_PRIO(o,prio,master,ghost,id,s)                          \
  if (USED(o)==1 && ! master (o))                                          \
    UserWriteF("MASTER %s=" id ## _FMTX " has WRONG prio=%d\n",      \
               s, id ## _PRTX(o),prio(o));                                  \
  if (USED( o )==0 && ! ghost ( o ))                                       \
    UserWriteF("GHOST %s=" id ## _FMTX " has WRONG prio=%d\n",       \
               s, id ## _PRTX( o ),prio(o));


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
    return;

  if (me!=min_proc)
    DDD_PrioritySet(PARHDR(node), PrioBorder);
}

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
    return;

  if (me!=min_proc)
    DDD_PrioritySet(PARHDR(vector), PrioBorder);
}

#ifdef __THREEDIM__
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
    return;

  if (me!=min_proc)
    DDD_PrioritySet(PARHDR(edge), PrioBorder);
}
#endif

void SetGhostObjectPriorities (GRID *theGrid)
{
  ELEMENT *theElement;
  NODE    *theNode;
  EDGE    *theEdge;
  INT i,prio,*proclist;

  /* reset USED flag for nodes of ghostelements */
  for (theElement=PFIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    SETUSED(theElement,0);
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      SETUSED(theNode,0);
    }
    for (i=0; i<EDGES_OF_ELEM(theElement); i++)
    {
      theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                        CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
      ASSERT(theEdge != NULL);
      SETUSED(theEdge,0);
    }
  }

  /* reset USED flag for nodes of master elements */
  for (theElement=FIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if (PARTITION(theElement) != me) continue;

    SETUSED(theElement,1);
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      SETUSED(theNode,1);
    }
    for (i=0; i<EDGES_OF_ELEM(theElement); i++)
    {
      theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                        CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
      ASSERT(theEdge != NULL);
      SETUSED(theEdge,1);
    }
  }

  /* set node priorities for ghostelements */
  for (theElement=PFIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if (TYPE_DEF_IN_GRID(theGrid,ELEMVECTOR))
    {
      if (USED(theElement) == 0)
        SETPRIO(EVECTOR(theElement),PrioGhost);
    }

    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);

      /* check if its a master node */
      if (USED(theNode) == 0)
      {
        PRINTDEBUG(dddif,3,(PFMT " dddif_SetGhostObjectPriorities():"
                            " downgrade node=" ID_FMTX " from=%d to PrioGhost\n",
                            me,ID_PRTX(theNode),prio));

        /* set node priorities of node to ghost */
        NODE_PRIORITY_SET(theGrid,theNode,PrioGhost)
      }
    }


    if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR) || DIM==3)
    {
      /* set edge priorities */
      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {

        theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                          CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
        ASSERT(theEdge != NULL);

        if (USED(theEdge) == 0)
        {
          PRINTDEBUG(dddif,3,(PFMT " dddif_SetGhostObjectPriorities():"
                              " downgrade edge=" EDID_FMTX " from=%d to PrioGhost\n",
                              me,EDID_PRTX(theEdge),prio));

          EDGE_PRIORITY_SET(theGrid,theEdge,PrioGhost);
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
              SETPRIO(SVECTOR(theElement,i),PrioGhost);
              break;
            }
          }
        }
                        #endif
    }
  }
}


void SetOverlapPriorities (GRID *theGrid)
{
  DDD_XferBegin();

  DDD_IFAExecLocal(BorderNodeSymmIF, GLEVEL(theGrid),
                   ComputeNodeBorderPrios);

  DDD_IFAExecLocal(BorderVectorSymmIF, GLEVEL(theGrid),
                   ComputeVectorBorderPrios);

#ifdef __THREEDIM__
  DDD_IFAExecLocal(BorderEdgeSymmIF, GLEVEL(theGrid),
                   ComputeEdgeBorderPrios);
#endif

  DDD_XferEnd();
}

void ConstructConsistentGrid (GRID *theGrid)
{
  INT i,j,k,l,m,o;
  DOUBLE fac,*local;
  ELEMENT *theElement,*theFather,*theNb;
  NODE    *theNode;
  EDGE    *theEdge;
  VERTEX  *theVertex;

  SetOverlapPriorities(theGrid);

    #ifdef __TWODIM__
  for (theVertex = PFIRSTVERTEX(theGrid); theVertex != NULL;
       theVertex = SUCCV(theVertex))
    if (OBJT(theVertex) == BVOBJ)
      if (MOVED(theVertex)) {
        INT n;
        DOUBLE *x[MAX_CORNERS_OF_ELEM];

        theElement = VFATHER(theVertex);
        if (theElement == NULL) continue;
        CORNER_COORDINATES(theElement,n,x);
        UG_GlobalToLocal(n,(const DOUBLE **)x,
                         CVECT(theVertex),LCVECT(theVertex));
      }
        #endif

        #ifdef __THREEDIM__
  /* reconstruct VFATHER pointers */
  for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      if (CORNERTYPE(theNode)) continue;

      theVertex = MYVERTEX(theNode);
      theFather = EFATHER(theElement);

      /*			if (VFATHER(theVertex)==NULL && theFather!=NULL) */
      if (VFATHER(theVertex)==NULL)
      {
        switch (NTYPE(theNode))
        {
        case (MID_NODE) :
        {
          INT co0,co1;

          for (j=0; j<EDGES_OF_ELEM(theFather); j++)
          {
            theEdge = GetEdge(CORNER(theFather,CORNER_OF_EDGE(theFather,j,0)),
                              CORNER(theFather,CORNER_OF_EDGE(theFather,j,1)));
            if (MIDNODE(theEdge) == theNode) break;
          }
          ASSERT(j<EDGES_OF_ELEM(theFather));

          /* reconstruct local coordinates of vertex */
          co0 = CORNER_OF_EDGE(theFather,j,0);
          co1 = CORNER_OF_EDGE(theFather,j,1);

          /* local coordinates have to be local towards pe */
          V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theFather,co0),
                        0.5, LOCAL_COORD_OF_ELEM(theFather,co1),
                        LCVECT(theVertex));
          SETONEDGE(theVertex,j);
          break;
        }

        case (SIDE_NODE) :
          /* always compute new coords for this case! */
          k =  GetSideIDFromScratch(theElement,theNode);
          ASSERT(k < SIDES_OF_ELEM(theFather));

          SETONSIDE(theVertex,k);

          m = CORNERS_OF_SIDE(theFather,k);
          local = LCVECT(theVertex);
          fac = 1.0 / m;
          V_DIM_CLEAR(local);
          for (o=0; o<m; o++)
          {
            l = CORNER_OF_SIDE(theFather,k,o);
            V_DIM_LINCOMB(1.0,local,1.0,
                          LOCAL_COORD_OF_ELEM(theFather,l),local);
          }
          V_DIM_SCALE(fac,local);

          theNb = NBELEM(theFather,k);
          if (theNb != NULL)
          {
            for (j=0; j<SIDES_OF_ELEM(theNb); j++)
            {
              if (NBELEM(theNb,j) == theFather) break;
            }
            ASSERT(j < SIDES_OF_ELEM(theNb));
            SETONNBSIDE(theVertex,j);
          }
          else SETONNBSIDE(theVertex,MAX_SIDES_OF_ELEM);
          break;

        case (CENTER_NODE) :
          /* nothing to do */
          break;

        case (CORNER_NODE) :
        default :
          assert(0);
          break;
        }
        VFATHER(theVertex) = theFather;

        if (OBJT(theVertex) == BVOBJ)
          if (MOVED(theVertex)) {
            INT n;
            DOUBLE *x[MAX_CORNERS_OF_ELEM];

            CORNER_COORDINATES(theFather,n,x);
            UG_GlobalToLocal(n,(const DOUBLE **)x,
                             CVECT(theVertex),LCVECT(theVertex));
          }
      }
    }
  }
        #endif

}

INT CheckProcListCons (int *proclist, int uniqueTag)
{
  int nunique = 0;

  /* check uniqueness */
  while (*proclist != -1)
  {
    if (*(proclist+1) == uniqueTag) nunique++;
    proclist += 2;
  }

  /* nunique must be 1 for master elements   */
  /* nunique can  be 0/1 for (inner) nodes   */
  /*   with PrioBorder/PrioMaster            */
  return (nunique);
}

INT ListProcList (int *proclist, int uniqueTag)
{
  while (*proclist != -1)
  {
    if (*(proclist+1) == uniqueTag)
      UserWriteF(" proc=%d",*proclist);
    proclist += 2;
  }
  return(0);
}

INT CheckVectorPrio (ELEMENT *theElement, VECTOR *theVector)
{
  INT nmaster;

  /* check vector prio */
  CHECK_OBJECT_PRIO(theVector,PRIO,MASTER,GHOST,ID,"Vector")

  /* master copy has to be unique */
  if ((nmaster = CheckProcListCons(PROCLIST(theVector),PrioMaster)) > 1)
  {
    UserWriteF("NODE=" ID_FMTX " ERROR: master copy not unique, nmaster=%d:",
               ID_PRTX(theVector),nmaster);
    ListProcList(PROCLIST(theVector),PrioMaster);
    UserWriteF("\n");
  }

  return(0);
}

INT CheckNodePrio (ELEMENT *theElement, NODE *theNode)
{
  INT nmaster;

  /* check node prio */
  CHECK_OBJECT_PRIO(theNode,PRIO,MASTER,GHOST,ID,"NODE")

  /* master copy has to be unique */
  if ((nmaster = CheckProcListCons(PROCLIST(theNode),PrioMaster)) > 1)
  {
    UserWriteF("NODE=" ID_FMTX " ERROR: master copy not unique, nmaster=%d:",
               ID_PRTX(theNode),nmaster);
    ListProcList(PROCLIST(theNode),PrioMaster);
    UserWriteF("\n");
  }

  if (dddctrl.nodeData)
    CheckVectorPrio(theElement,NVECTOR(theNode));

  return(0);
}


INT CheckEdgePrio (ELEMENT *theElement, EDGE *theEdge)
{
  INT nmaster;

        #ifdef __THREEDIM__
  /* check edge prio */
  CHECK_OBJECT_PRIO(theEdge,PRIO,MASTER,GHOST,ID,"EDGE")

  /* master copy has to be unique */
  if ((nmaster = CheckProcListCons(PROCLIST(theEdge),PrioMaster)) > 1)
  {
    UserWriteF("EDGE=" EDID_FMTX " ERROR: master copy not unique, nmaster=%d:",
               EDID_PRTX(theEdge),nmaster);
    ListProcList(PROCLIST(theEdge),PrioMaster);
    UserWriteF("\n");
  }
        #endif

  if (dddctrl.edgeData)
    CheckVectorPrio(theElement,EDVECTOR(theEdge));

  return(0);
}

INT CheckElementPrio (ELEMENT *theElement)
{
  INT i,nmaster,prio;
  NODE    *theNode;
  EDGE    *theEdge;
  ELEMENT *SonList[MAX_SONS];

  if (PARTITION(theElement)==me && !EMASTER(theElement))
  {
    UserWriteF(PFMT "#FATAL# MASTER ELEM=" EID_FMTX " has WRONG part=%d prio=%d\n",
               me,EID_PRTX(theElement),PARTITION(theElement),EPRIO(theElement));
  }
  if (PARTITION(theElement)!=me && !EGHOST(theElement))
  {
    UserWriteF(PFMT "#FATAL# GHOST ELEM=" EID_FMTX " has WRONG part=%d prio=%d\n",
               me,EID_PRTX(theElement),PARTITION(theElement),EPRIO(theElement));

    /* test ghost prio */
    prio = 0;
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      if (EMASTER(NBELEM(theElement,i))) prio = PrioGhost;
    }
    if (GetSons(theElement,SonList) != 0) RETURN(1);
    if (SonList[0] != NULL) prio += PrioVGhost;

    if (EPRIO(theElement) != prio)
    {
      UserWriteF(PFMT "ERROR GHOST ELEM=" EID_FMTX
                 " has WRONG prio=%d should be prio=%d\n",
                 me,EID_PRTX(theElement),EPRIO(theElement),prio);
    }
  }

  /* check element prio */
  CHECK_OBJECT_PRIO(theElement,EPRIO,EMASTER,EGHOST,EID,"ELEM")

  if (dddctrl.elemData)
    CheckVectorPrio(theElement,EVECTOR(theElement));

  if (dddctrl.sideData)
  {
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      CheckVectorPrio(theElement,SVECTOR(theElement,i));
  }

  /* master copy has to be unique */
  if ((nmaster = CheckProcListCons(EPROCLIST(theElement),PrioMaster)) != 1)
  {
    UserWriteF("ELEM=" EID_FMTX " ERROR: master copy not unique, nmaster:",
               EID_PRTX(theElement),nmaster);
    ListProcList(EPROCLIST(theElement),PrioMaster);
    UserWriteF("\n");
  }

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    theNode = CORNER(theElement,i);
    CheckNodePrio(theElement,theNode);
  }

  for (i=0; i<EDGES_OF_ELEM(theElement); i++)
  {
    theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                      CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
    ASSERT(theEdge != NULL);
    CheckEdgePrio(theElement,theEdge);
  }

  return (0);
}

INT CheckInterfaces(GRID *theGrid)
{
  INT i,j;
  ELEMENT *theElement;
  NODE    *theNode;
  EDGE    *theEdge;
  VECTOR  *theVector;
  int errors = 0;

  /* reset USED flag of all grid objects  */
  /* set USED flag of master grid objects */
  for (j=0; j<2; j++)
  {
    for (theElement =(j==0 ? PFIRSTELEMENT(theGrid) : FIRSTELEMENT(theGrid));
         theElement!=NULL;
         theElement=SUCCE(theElement))
    {
      SETUSED(theElement,j);
      if (dddctrl.elemData)
        SETUSED(EVECTOR(theElement),j);
      if (dddctrl.sideData)
      {
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
          SETUSED(SVECTOR(theElement,i),j);
      }

      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        theNode = CORNER(theElement,i);
        SETUSED(theNode,j);
        if (dddctrl.nodeData)
          SETUSED(NVECTOR(theNode),j);
        SETUSED(MYVERTEX(theNode),j);
      }

      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {
        theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                          CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
        ASSERT(theEdge != NULL);
        SETUSED(theEdge,j);
        if (dddctrl.edgeData)
          SETUSED(EDVECTOR(theEdge),j);
      }
    }
  }

  /* check validity of priorities */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    errors += CheckElementPrio(theElement);
  }

  /* check ddd interface consistency */
  errors += DDD_ConsCheck();

  return(errors);
}

/****************************************************************************/

#endif
