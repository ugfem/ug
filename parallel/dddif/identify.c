// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  identify.c													*/
/*																			*/
/* Purpose:   identification of distributed ug objects                                  */
/*																			*/
/* Author:	  Stefan Lang                                                                                   */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   26.11.96 begin, first version extracted from refine.c			*/
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

#include "general.h"
#include "compiler.h"
#include "debug.h"
#include "gm.h"
#include "rm.h"
#include "refine.h"
#include "ddd.h"
#include "parallel.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* flags for identification */
#define CLEAR 0
#define IDENT 1

/* maximum count of objects for identification */
#define MAX_OBJECT      3

/* maximum count of tokens for identification */
#define MAX_TOKEN       10

/* determine the ddd header for identification of a node */
#define NIDENT_HDR(node) ( (CORNERTYPE(node)) ? \
                           PARHDR((NODE *)NFATHER(node)) : PARHDR(node) )

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

#define NIDENT(p)                     THEFLAG(p)
#define SETNIDENT(p,n)                SETTHEFLAG(p,n)

#define EDIDENT(p)                    THEFLAG(p)
#define SETEDIDENT(p,n)               SETTHEFLAG(p,n)


/* this function is called for low level identification */
static INT (*Ident_FctPtr)(DDD_HDR *IdentObjectHdr, INT nobject,
                           int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident) = NULL;

#ifdef Debug
static INT debug = 0;
static INT identlevel = 0;
#endif

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/


/****************************************************************************/
/*
   compare_gid -

   SYNOPSIS:
   INT compare_gid (const void *e0, const void *e1);

   PARAMETERS:
   .  e0
   .  e1

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT compare_gid (const void *e0, const void *e1)
{
  INT num0, num1;

  num0 = DDD_InfoGlobalId(*(DDD_HDR *)e0);
  num1 = DDD_InfoGlobalId(*(DDD_HDR *)e1);

  if (num0 < num1) return(1);
  if (num0 > num1) return(-1);
  return(0);
}


/****************************************************************************/
/*
   ResetIdentFlags -

   SYNOPSIS:
   static void ResetIdentFlags (GRID *UpGrid);

   PARAMETERS:
   .  UpGrid

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void ResetIdentFlags (GRID *Grid)
{
  NODE *theNode;
  EDGE *theEdge;
  LINK *theLink;

  /* clear all IDENT flags */
  for (theNode=FIRSTNODE(Grid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    SETNIDENT(theNode,CLEAR);
    SETUSED(theNode,0);

    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      theEdge = MYEDGE(theLink);
      SETEDIDENT(theEdge,CLEAR);
    }
  }

}

#ifdef Debug


/****************************************************************************/
/*
   Print_Identify_ObjectList -

   SYNOPSIS:
   static INT Print_Identify_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject,
                                int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident);

   PARAMETERS:
   .  IdentObjectHdr
   .  nobject
   .  proclist
   .  skiptag
   .  IdentHdr
   .  nident

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT Print_Identify_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject,
                                      int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
{
  INT i;

  ASSERT(nobject>0);
  ASSERT(nident>0);
  ASSERT(*proclist!=-1);

  /* print the interesting call parameters */
  PrintDebug ("%d:    Print_Identify_ObjectList(): nobject=%d nident=%d"
              " skiptag=%d\n",me,nobject,nident,skiptag);

  /* print line prefix */
  PrintDebug ("%d: l=%d",me,identlevel);

  /* print the objects used for identify */
  PrintDebug ("    IdentHdr:");
  for (i=0; i<nident; i++)
    PrintDebug (" %d",DDD_InfoGlobalId(IdentHdr[i]));

  /* print type of objects to identify */
  PrintDebug ("    IdentObjectType:");
  for (i=0; i<nobject; i++)
    PrintDebug (" %d",DDD_InfoType(IdentObjectHdr[i]));

  /* print the proclist to identify to */
  PrintDebug ("    ProcList: %d",me);
  while (*proclist != -1)
  {
    if (*(proclist+1) == skiptag)
    {
      proclist += 2;
      continue;
    }

    PrintDebug (" %d",*proclist);
    proclist += 2;
  }

  /* print my processor number */
  PrintDebug ("    me:%d",me);

  /* print the objects to identify */
  PrintDebug ("    IdentObjectHdr:");
  for (i=0; i<nobject; i++)
    PrintDebug (" %d",DDD_InfoGlobalId(IdentObjectHdr[i]));

  PrintDebug ("\n");

  return(0);

}
#endif

#ifdef Debug


/****************************************************************************/
/*
   Print_Identified_ObjectList -

   SYNOPSIS:
   static INT Print_Identified_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject, int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident);

   PARAMETERS:
   .  IdentObjectHdr
   .  nobject
   .  proclist
   .  skiptag
   .  IdentHdr

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT Print_Identified_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject,
                                        int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
{
  INT i;

  ASSERT(nobject>0);
  ASSERT(nident>0);
  ASSERT(*proclist!=-1);

  /* print the interesting call parameters */
  PrintDebug ("%d:    Print_Identified_ObjectList(): nobject=%d nident=%d"
              " skiptag=%d\n",me,nobject,nident,skiptag);

  /* print the objects to identify */
  PrintDebug ("%d: l=%d   IdentObjectHdr:",me,identlevel);
  for (i=0; i<nobject; i++)
    PrintDebug (" %d",DDD_InfoGlobalId(IdentObjectHdr[i]));

  /* print the objects used for identify */
  PrintDebug ("    IdentHdr:");
  for (i=0; i<nident; i++)
    PrintDebug (" %d",DDD_InfoGlobalId(IdentHdr[i]));

  /* print the proclist to identify to */
  PrintDebug ("    ProcList: %d",me);
  while (*proclist != -1)
  {
    if (*(proclist+1) == skiptag)
    {
      proclist += 2;
      continue;
    }

    PrintDebug (" %d",*proclist);
    proclist += 2;
  }

  /* print my processor number */
  PrintDebug ("    me:%d",me);

  /* print type of objects to identify */
  PrintDebug ("    IdentObjectType:");
  for (i=0; i<nobject; i++)
    PrintDebug (" %d",DDD_InfoType(IdentObjectHdr[i]));

  PrintDebug ("\n");

  return(0);
}
#endif



/****************************************************************************/
/*
   Identify_by_ObjectList -

   SYNOPSIS:
   static INT Identify_by_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject,
                                int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident);

   PARAMETERS:
   .  IdentObjectHdr
   .  nobject
   .  proclist
   .  skiptag
   .  IdentHdr
   .  nident

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT Identify_by_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject,
                                   int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
{
  INT i,j,n;

  ASSERT(nobject>0);
  ASSERT(nident>0);
  ASSERT(*proclist!=-1);

  IFDEBUG(graph,1)
  Print_Identify_ObjectList(IdentObjectHdr,nobject,proclist,skiptag,
                            IdentHdr,nident);
  ENDDEBUG

    n = 0;
  while (*proclist != -1)
  {
    ASSERT(n<procs);

    if (*(proclist+1) == skiptag)
    {
      proclist += 2;
      continue;
    }

    /* identify the object */
    for (j=0; j<nobject; j++)
    {
      for (i=0; i<nident; i++)
      {
        PRINTDEBUG(dddif,5,("%d: Identify_by_ObjectList(): Type=%d"
                            " IdentObjectHdr=%08x proclist=%d IdentHdr=%08x me=%d\n",
                            me,DDD_InfoType(IdentObjectHdr[j]),
                            DDD_InfoGlobalId(IdentObjectHdr[j]),
                            *proclist,DDD_InfoGlobalId(IdentHdr[i]),me));

        /* hand identification hdr to ddd */
        DDD_IdentifyObject(IdentObjectHdr[j], *proclist, IdentHdr[i]);
      }
    }

    n++;
    assert(n<procs);
    proclist += 2;
  }

  /* identification should occur to at least one other proc */
  ASSERT(n>0);
}

#ifdef __THREEDIM__


/****************************************************************************/
/*
   IdentifySideVector -

   SYNOPSIS:
   static INT IdentifySideVector (ELEMENT* theElement, ELEMENT *theNeighbor,
                                ELEMENT *Son, INT SonSide) ;

   PARAMETERS:
   .  theElement
   .  theNeighbor
   .  theSon
   .  SonSide

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT IdentifySideVector (ELEMENT* theElement, ELEMENT *theNeighbor,
                               ELEMENT *Son, INT SonSide)
{
  INT k,nident;
  DDD_HDR IdentObjectHdr[MAX_OBJECT];
  DDD_HDR IdentHdr[MAX_TOKEN];
  int *proclist;
  NODE *theNode;

  nident = 0;

  IdentObjectHdr[0] = PARHDR(SVECTOR(Son,SonSide));

  proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

  /* identify using corner nodes */
  for (k=0; k<CORNERS_OF_SIDE(Son,SonSide); k++)
  {
    theNode = CORNER(Son,CORNER_OF_SIDE(Son,SonSide,k));
    if (CORNERTYPE(theNode))
      IdentHdr[nident++] = PARHDR((NODE *)NFATHER(theNode));
    else
      IdentHdr[nident++] = PARHDR(theNode);
  }

  proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

  Ident_FctPtr(IdentObjectHdr,1,proclist+2,PrioHGhost,IdentHdr,nident);

}
#endif


/****************************************************************************/
/*
   IdentifyNode -

   SYNOPSIS:
   static void IdentifyNode (ELEMENT *theNeighbor, NODE *theNode,
                                NODE *Nodes[MAX_SIDE_NODES], INT node, INT ncorners, INT Vec);

   PARAMETERS:
   .  theNeighbor
   .  theNode
   .  Nodes[MAX_SIDE_NODES]
   .  node
   .  ncorners
   .  Vec

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void IdentifyNode (GRID *theGrid, ELEMENT *theNeighbor, NODE *theNode,
                          NODE *Nodes[MAX_SIDE_NODES], INT node, INT ncorners, INT Vec)
{
  INT nobject,nident;
  DDD_HDR IdentObjectHdr[MAX_OBJECT];
  DDD_HDR IdentHdr[MAX_TOKEN];

  nobject = nident = 0;

  /* is this node identified? */
        #ifdef Debug
  if (debug == 1) {
    if (NIDENT(theNode) == CLEAR) return;
  }
  else
        #endif

  /* return if not needed any more */
  if (!USED(theNode)) return;

  /* return if already identified */
  if (NIDENT(theNode) == IDENT) return;

  /* only new created nodes are identified */
  if (!NEW_NIDENT(theNode)) return;

  if (Vec)
    if (GetVectorSize(theGrid,NODEVEC,(GEOM_OBJECT *)theNode) == 0)
      Vec = 0;

  switch (NTYPE(theNode))
  {
    int *proclist;

  case (CORNER_NODE) :

    /* identification of cornernodes is done */
    /* in Identify_SonNodes()     */
    return;

    PRINTDEBUG(dddif,1,("%d: Identify CORNERNODE gid=%08x node=%d "
                        "vec=%d\n",
                        me, DDD_InfoGlobalId(PARHDR(theNode)), node, Vec));

    IdentObjectHdr[nobject++] = PARHDR(theNode);
    if (Vec)
      IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

    /* identify to proclist of node */
    proclist = DDD_InfoProcList(PARHDR((NODE *)NFATHER(theNode)));

    /* identify using father node */
    IdentHdr[nident++] = PARHDR((NODE *)NFATHER(theNode));

    Ident_FctPtr(IdentObjectHdr, nobject,
                 proclist+2, PrioHGhost, IdentHdr, nident);

    break;

  case (MID_NODE) :
  {
                        #ifdef __TWODIM__
    NODE **EdgeNodes;
    EdgeNodes = Nodes;
                        #endif

                        #ifdef __THREEDIM__
    NODE *EdgeNodes[MAX_SIDE_NODES];
    EDGE *theEdge;

    /* identification of cornernodes is done */
    /* in Identify_SonEdges()     */
    return;

    EdgeNodes[0] = Nodes[node-ncorners];
    EdgeNodes[1] = Nodes[(node-ncorners+1)%ncorners];
    EdgeNodes[2] = theNode;
                        #endif

    ASSERT(EdgeNodes[0]!=NULL);
    ASSERT(EdgeNodes[1]!=NULL);
    ASSERT(EdgeNodes[2]!=NULL);

    PRINTDEBUG(dddif,1,("%d: Identify MIDNODE gid=%08x node=%d Vec=%d\n",
                        me, DDD_InfoGlobalId(PARHDR(theNode)), node, Vec));

    /* identify midnode, vertex, vector */
    IdentObjectHdr[nobject++] = PARHDR(theNode);
    IdentObjectHdr[nobject++] = PARHDRV(MYVERTEX(theNode));
    if (Vec)
      IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

                        #ifdef __TWODIM__
    if (!NEW_NIDENT(theNode)) break;

    /* 2D: identify to proclist of neighbor element */
    proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
                        #endif

                        #ifdef __THREEDIM__
    /* 3D: identify to proclist of edge */
    theEdge = GetEdge((NODE *)NFATHER(EdgeNodes[0]),
                      (NODE *)NFATHER(EdgeNodes[1]));
    ASSERT(theEdge!=NULL);

    proclist = DDD_InfoProcList(PARHDR(theEdge));
                        #endif

    /* identify using edge nodes */
    IdentHdr[nident++] = PARHDR((NODE *)NFATHER(EdgeNodes[0]));
    IdentHdr[nident++] = PARHDR((NODE *)NFATHER(EdgeNodes[1]));

    /* this is the buggy case
                            IdentHdr[nident++] = PARHDR(EdgeNodes[0]);
                            IdentHdr[nident++] = PARHDR(EdgeNodes[1]);
     */

    Ident_FctPtr(IdentObjectHdr, nobject,
                 proclist+2, PrioHGhost, IdentHdr, nident);

    break;
  }

                #ifdef __THREEDIM__
  case (SIDE_NODE) :
  {
    INT i;

    PRINTDEBUG(dddif,1,("%d: Identify SIDENODE gid=%08x node=%d "
                        "Vec=%d\n",me,DDD_InfoGlobalId(PARHDR(theNode)),node,Vec));

    /* identify sidenode, vertex and vector */
    IdentObjectHdr[nobject++] = PARHDR(theNode);
    IdentObjectHdr[nobject++] = PARHDRV(MYVERTEX(theNode));
    if (Vec)
      IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

    /* identify to proclist of neighbor element */
    proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

    /* identify using corner nodes of side */
    for (i=0; i<ncorners; i++)
      IdentHdr[nident++] = PARHDR((NODE *)NFATHER(Nodes[i]));

    /* identify side node */
    Ident_FctPtr(IdentObjectHdr, nobject,
                 proclist+2, PrioHGhost, IdentHdr, nident);

    break;
  }
                #endif
  default :
    ASSERT(0);
    break;
  }

        #ifdef Debug
  if (debug == 1) {
    SETNIDENT(theNode,CLEAR);
  }
  else
        #endif
  /* lock this node for identification */
  SETNIDENT(theNode,IDENT);

  return;
}



/****************************************************************************/
/*
   IdentifyEdge -

   SYNOPSIS:
   static INT IdentifyEdge (ELEMENT *theElement, ELEMENT *theNeighbor,
                        NODE **SideNodes, INT ncorners, ELEMENT *Son, INT SonSide, INT edgeofside, INT Vec);

   PARAMETERS:
   .  theElement
   .  theNeighbor
   .  SideNodes
   .  ncorners
   .  Son
   .  SonSide
   .  edgeofside
   .  Vec

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT IdentifyEdge (GRID *theGrid,
                         ELEMENT *theElement, ELEMENT *theNeighbor,
                         NODE **SideNodes, INT ncorners, ELEMENT *Son,
                         INT SonSide, INT edgeofside, INT Vec)
{
  NODE *Nodes[2];
  EDGE *theEdge;
  VECTOR *theVector;
  INT nobject,nident;
        #ifdef __THREEDIM__
  INT edge,corner0,corner1;
        #endif
  INT *proclist;
  DDD_HDR IdentObjectHdr[MAX_OBJECT];
  DDD_HDR IdentHdr[MAX_TOKEN];

  nobject = nident = 0;

        #ifdef __TWODIM__
  Nodes[0] = CORNER(Son,CORNER_OF_EDGE(Son,SonSide,0));
  Nodes[1] = CORNER(Son,CORNER_OF_EDGE(Son,SonSide,1));
        #endif

        #ifdef __THREEDIM__
  edge = EDGE_OF_SIDE(Son,SonSide,edgeofside);
  corner0 = CORNER_OF_EDGE(Son,edge,0);
  corner1 = CORNER_OF_EDGE(Son,edge,1);

  Nodes[0] = CORNER(Son,corner0);
  Nodes[1] = CORNER(Son,corner1);
  PRINTDEBUG(dddif,5,("%4d: edge=%d corner0=%d corner1=%d Nodes[0]=%d "
                      "Nodes[1]=%d\n",
                      me,edge, corner0, corner1, ID(Nodes[0]), ID(Nodes[1])));
        #endif

  ASSERT(Nodes[0]!=NULL);
  ASSERT(Nodes[1]!=NULL);

  theEdge = GetEdge(Nodes[0],Nodes[1]);
  ASSERT(theEdge!=NULL);

        #ifdef __TWODIM__
  /* no identfication to nonrefined neighbors */
  if (MARK(theNeighbor) == NO_REFINEMENT) return(0);
        #endif

        #ifdef __THREEDIM__
  /* identification of sonedges is done in Identify_SonEdges() */
  if (0)
    if (CORNERTYPE(Nodes[0]) && CORNERTYPE(Nodes[1]))
    {
      EDGE *FatherEdge;
      FatherEdge = GetEdge((NODE *)NFATHER(Nodes[0]),(NODE *)NFATHER(Nodes[1]));
      ASSERT(FatherEdge != NULL);
      return(0);
      /*
              if (FatherEdge != NULL) return(0);
       */
    }
  {
    EDGE *FatherEdge;
    FatherEdge = GetFatherEdge(theEdge);
    if (FatherEdge != NULL) return(0);
  }
        #endif

  /* only newly created edges are identified */
  if (!NEW_EDIDENT(theEdge)) return(0);

  /* edge unlocked -> no debugging occurs */
        #ifdef Debug
  if (debug == 1) {
    if (EDIDENT(theEdge) == CLEAR) return(0);
  }
  else
        #endif

  /* edge locked -> already identified */
  if (EDIDENT(theEdge) == IDENT) return(0);

  PRINTDEBUG(dddif,1,("%d: Identify EDGE edgeofside=%d pe=%08x/%x eID=%d"
                      " ntype0=%d  ntype1=%d Vec=%d\n",me,edgeofside,
                      DDD_InfoGlobalId(PARHDRE(Son)),Son,ID(Son),
                      NTYPE(Nodes[0]), NTYPE(Nodes[1]), Vec))

        #ifdef __THREEDIM__
  IdentObjectHdr[nobject++] = PARHDR(theEdge);
        #endif
  if (Vec)
    if (GetVectorSize(theGrid,EDGEVEC,(GEOM_OBJECT *)theEdge) > 0)
      if (EDVECTOR(theEdge) != NULL)
        IdentObjectHdr[nobject++] = PARHDR(EDVECTOR(theEdge));

        #ifdef __TWODIM__
  /* identify to proclist of neighbor */
  proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
        #endif

  /* identify to proclist of father edge or neighbor*/
        #ifdef __THREEDIM__
  {
    EDGE *fatherEdge = NULL;

    /* check whether edge inside the side of the element */
    fatherEdge = FatherEdge(SideNodes,ncorners,Nodes,theEdge);

    if (fatherEdge != NULL)
      proclist = DDD_InfoProcList(PARHDR(fatherEdge));
    else
      proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
  }
        #endif

  if (CORNERTYPE(Nodes[0]))
    IdentHdr[nident++] = PARHDR((NODE *)NFATHER(Nodes[0]));
        #ifdef __THREEDIM__
  /* since midnodes are identified later in Debug case */
  /* choose fatheredge here (s.l. 980227)              */
  else if (MIDTYPE(Nodes[0]))
    IdentHdr[nident++] = PARHDR((EDGE *)NFATHER(Nodes[0]));
        #endif
  else
    IdentHdr[nident++] = PARHDR(Nodes[0]);

  if (CORNERTYPE(Nodes[1]))
    IdentHdr[nident++] = PARHDR((NODE *)NFATHER(Nodes[1]));
        #ifdef __THREEDIM__
  /* since midnodes are identified later in Debug case */
  /* choose fatheredge here (s.l. 980227)              */
  else if (MIDTYPE(Nodes[1]))
    IdentHdr[nident++] = PARHDR((EDGE *)NFATHER(Nodes[1]));
        #endif
  else
    IdentHdr[nident++] = PARHDR(Nodes[1]);

  if (nobject > 0)
    Ident_FctPtr(IdentObjectHdr, nobject,
                 proclist+2, PrioHGhost, IdentHdr, nident);

  /* debugging unlocks the edge */
        #ifdef Debug
  if (debug == 1) {
    SETEDIDENT(theEdge,CLEAR);
  }
  else
        #endif
  /* lock this edge for identification */
  SETEDIDENT(theEdge,IDENT);

  return(0);
}


/****************************************************************************/
/*
   IdentifyObjectsOfElementSide -

   SYNOPSIS:
   static INT IdentifyObjectsOfElementSide(GRID *theGrid, ELEMENT *theElement,
                                                INT i, ELEMENT *theNeighbor);

   PARAMETERS:
   .  theGrid
   .  theElement
   .  i
   .  theNeighbor

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT IdentifyObjectsOfElementSide(GRID *theGrid, ELEMENT *theElement,
                                        INT i, ELEMENT *theNeighbor)
{
  INT nodes,j,n;
  NODE *SideNodes[MAX_SIDE_NODES];
  INT ncorners;
  NODE *theNode;

  GetSonSideNodes(theElement,i,&nodes,SideNodes,0);
  ncorners = CORNERS_OF_SIDE(theElement,i);
  n = 0;

  PRINTDEBUG(dddif,1,("%d: IdentifyObjectsOfElementSide():identify NODES "
                      "ncorners=%d nodes=%d\n",me,ncorners,nodes));

  /* identify nodes, vertices and node vectors of son elements */
  for (j=0; j<MAX_SIDE_NODES; j++)
  {
    theNode = SideNodes[j];
    if (theNode == NULL) continue;

    /* identify new node including its vector and vertex        */
    IdentifyNode(theGrid,theNeighbor, theNode, SideNodes, j, ncorners,
                 VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC));
    n++;
  }
  ASSERT(n == nodes);

  /* identify edge vectors (2D); edges, edge and side vectors (3D) */
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC) || DIM==3)
  {
    ELEMENT *SonList[MAX_SONS];
    INT SonsOfSide,SonSides[MAX_SONS];
    INT j;

    PRINTDEBUG(dddif,1,("%d: IdentifyObjectsOfElementSide(): identify "
                        "EDGES and VECTORS\n",me));

    if (Get_Sons_of_ElementSide(theElement,i,&SonsOfSide,
                                SonList,SonSides,1,0)!=GM_OK)
      RETURN(GM_FATAL);

    for (j=0; j<SonsOfSide; j++) {

      if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC) || DIM==3)
      {
        INT edgeofside;
        INT nedges = EDGES_OF_SIDE(SonList[j],SonSides[j]);

        /* identify the edge and vector */
        for (edgeofside=0; edgeofside<nedges; edgeofside++) {
          IdentifyEdge(theGrid,
                       theElement,theNeighbor,SideNodes,ncorners,
                       SonList[j],SonSides[j],edgeofside,
                       VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC));
        }
      }

                        #ifdef __THREEDIM__
      if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
      {
        IdentifySideVector(theElement,theNeighbor,SonList[j],SonSides[j]);
        /* this is not debugged */
        assert(0);
      }
                        #endif
    }
  }

  return(GM_OK);
}



/****************************************************************************/
/*
   IdentifyDistributedObjects -

   SYNOPSIS:
   INT	IdentifyDistributedObjects (MULTIGRID *theMG, INT FromLevel, INT ToLevel);

   PARAMETERS:
   .  theMG
   .  FromLevel
   .  ToLevel

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT     IdentifyDistributedObjects (MULTIGRID *theMG, INT FromLevel, INT ToLevel)
{
  INT l,i,j,prio;
  ELEMENT *theElement,*theNeighbor;
  NODE *theNode;
  GRID *theGrid;

  PRINTDEBUG(dddif,1,("%d: IdentifyDistributedObjects(): FromLevel=%d "
                      "ToLevel=%d\n",me,FromLevel,ToLevel));

  /* identify distributed objects */
  for (l=FromLevel; l<ToLevel; l++)
  {
    PRINTDEBUG(dddif,1,("%d: IdentifyDistributedObjects(): identification "
                        "level=%d\n",me,l));

    theGrid = GRID_ON_LEVEL(theMG,l);

                #ifdef Debug
    identlevel = l;
                #endif

    /* check control word flags for ident on upper level */
    if (debug != 1)
      ResetIdentFlags(GRID_ON_LEVEL(theMG,l+1));

    for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement))
    {
      prio = EPRIO(theElement);

      if (!IS_REFINED(theElement) || EGHOSTPRIO(prio)) continue;

      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        theNeighbor = NBELEM(theElement,i);
        if (theNeighbor == NULL) continue;

        /* TODO: change for full dynamic element distribution */
        prio = EPRIO(theNeighbor);
        if (!HGHOSTPRIO(prio) || NSONS(theNeighbor)!=0) continue;

        PRINTDEBUG(dddif,1,("%d: Identify element: pe=%08x/%x eID=%d "
                            "side=%d\n",me,DDD_InfoGlobalId(PARHDRE(theElement)),
                            theElement,ID(theElement),i));

        IdentifyObjectsOfElementSide(theGrid,theElement,i,theNeighbor);

      }
    }
  }

  return(GM_OK);
}


/****************************************************************************/
/*
   IdentifyGridLevels -

   SYNOPSIS:
   INT	IdentifyGridLevels (MULTIGRID *theMG, INT FromLevel, INT ToLevel);

   PARAMETERS:
   .  theMG
   .  FromLevel
   .  ToLevel

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT     IdentifyGridLevels (MULTIGRID *theMG, INT FromLevel, INT ToLevel)
{
        #ifdef Debug
  debug = 0;
        #endif

  /* allocate a control word entry to lock nodes */
  if (AllocateControlEntry(NODE_CW,NEW_NIDENT_LEN,&ce_NEW_NIDENT) != GM_OK)
    assert(0);

  /* allocate a control word entry to lock edges */
  if (AllocateControlEntry(EDGE_CW,NEW_EDIDENT_LEN,&ce_NEW_EDIDENT) != GM_OK)
    assert(0);

  /* set Ident_FctPtr to identification mode */
  Ident_FctPtr = Identify_by_ObjectList;

  /* identify new created objects */
  DDD_IdentifyBegin();

  /* give identification calls to DDD */
  IdentifyDistributedObjects(theMG,FromLevel,ToLevel);

  /* start identification process */
  DDD_IdentifyEnd();

  /* print formated info for all identified objects */
  IFDEBUG(dddif,1)

  debug = 1;

  /* set Ident_FctPtr to print mode */
  Ident_FctPtr = Print_Identified_ObjectList;

  PrintDebug("AFTER Identify\n");
  IdentifyDistributedObjects(theMG,FromLevel,ToLevel);

  ENDDEBUG

  FreeControlEntry(ce_NEW_NIDENT);
  FreeControlEntry(ce_NEW_EDIDENT);
}

#ifdef IDENT_ONLY_NEW
static int Gather_NewNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE *theNode = (NODE *)obj;
  NODE *SonNode = SONNODE(theNode);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theNode));

  if (SonNode!=NULL && NEW_NIDENT(SonNode))
  {
    IFDEBUG(dddif,1)
    UserWriteF(PFMT "new son node=" ID_FMTX  "node=" ID_FMTX "\n",
               me,ID_PRTX(SonNode),ID_PRTX(theNode));
    ENDDEBUG
    *((int *)data) = 1;
  }
  else
    *((int *)data) = 0;

  return(0);
}

static int Scatter_NewNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE    *theNode        = (NODE *)obj;
  NODE    *SonNode        = SONNODE(theNode);
  int has_newsonnode  = *((int *)data);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theNode));

  if (SonNode!=NULL && has_newsonnode) SETNEW_NIDENT(SonNode,1);

  return(0);
}

/*************************/

static int Gather_NodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE *theNode = (NODE *)obj;

  ASSERT(identlevel == LEVEL(theNode));

  if (!CORNERTYPE(theNode))
  {
    *((int *)data) = 0;
    return(0);
  }

  *((int *)data) = NEW_NIDENT(theNode);

  return(0);
}

#define NIDENTASSERT 1

static int Scatter_NodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE    *theNode = (NODE *)obj;
  int nprop    = *((int *)data);

  ASSERT(identlevel == LEVEL(theNode));

  if (!CORNERTYPE(theNode)) return(0);

  if (NIDENTASSERT) if (NEW_NIDENT(theNode)) assert(NFATHER(theNode) != NULL);

  if (nprop)
  {
    SETNEW_NIDENT(theNode,1);
    if (NFATHER(theNode) == NULL)
    {
      UserWriteF(PFMT "isolated node=" ID_FMTX "\n",
                 me,ID_PRTX(theNode));
      if (NIDENTASSERT) assert(0);
    }
    if (NIDENTASSERT) assert(NFATHER(theNode) != NULL);
  }

  return(0);
}

/*************************/

static int Gather_TestNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE *theNode = (NODE *)obj;

  ASSERT(identlevel == LEVEL(theNode));

  ((int *)data)[0] = NEW_NIDENT(theNode);
  if (NEW_NIDENT(theNode)) assert(NFATHER(theNode) != NULL);

  return(0);
}

static int Scatter_TestNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE    *theNode        = (NODE *)obj;
  int nprop   = *((int *)data);

  ASSERT(identlevel == LEVEL(theNode));

  if (NEW_NIDENT(theNode) != nprop)
  {
    UserWriteF(PFMT "nprop wrong mynprop=%d hisnprop=%d theNode=" ID_FMTX " LEVEL=%d PROC=%d PRIO=%d\n",
               me,NEW_NIDENT(theNode),nprop,ID_PRTX(theNode),LEVEL(theNode),proc,prio);
    fflush(stdout);
    assert(0);
  }

  return(0);
}

/*************************/

static int Gather_IdentSonNode (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE *theNode = (NODE *)obj;
  NODE *SonNode = SONNODE(theNode);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theNode));

  ((int *)data)[0] = 0;
  ((int *)data)[1] = 0;

  if (SonNode != NULL)
  {
    ((int *)data)[0] = 1;
    ((int *)data)[1] = NEW_NIDENT(SonNode);
  }

  return(0);
}

static int Scatter_IdentSonNode (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE    *theNode        = (NODE *)obj;
  NODE    *SonNode        = SONNODE(theNode);
  int sonnode         = ((int *)data)[0];
  int newsonnode      = ((int *)data)[1];

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theNode));

  if (SonNode!=NULL)
  {
    /*
            if (1 || NEW_NIDENT(SonNode))
     */
    if (NEW_NIDENT(SonNode))
    {
      if(sonnode)
      {
        if (!newsonnode)
        {
          UserWriteF(PFMT "theNode=" ID_FMTX " LEVEL=%d PROC=%d PRIO=%d sonnprop=%d\n",
                     me,ID_PRTX(theNode),LEVEL(theNode),proc,prio,NEW_NIDENT(SonNode));
          fflush(stdout);
          assert(0);
        }

        DDD_IdentifyObject(PARHDR(SonNode),proc,PARHDR(theNode));
        if (dddctrl.nodeData && NVECTOR(SonNode)!=NULL)
          DDD_IdentifyObject(PARHDR(NVECTOR(SonNode)),proc,PARHDR(theNode));
      }
    }
    else
    {
      if (newsonnode)
      {
        UserWriteF(PFMT "theNode=" ID_FMTX " LEVEL=%d PROC=%d PRIO=%d sonnprop=%d\n",
                   me,ID_PRTX(theNode),LEVEL(theNode),proc,prio,NEW_NIDENT(SonNode));
        fflush(stdout);
        assert(0);
      }
    }
  }

  return(0);
}

/* callback functions for edge identification */
#ifdef __THREEDIM__
static int Gather_NewObjectInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  INT nedges, ident_needed;
  EDGE *theEdge   = (EDGE *)obj;
  EDGE *SonEdges[MAX_SON_EDGES];
  NODE *MidNode   = MIDNODE(theEdge);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  nedges = GetSonEdges(theEdge,SonEdges);

  /* identification is done if one objects of MidNode and the one */
  /* or two sonedge have NEW_XXIDENT flags set.                   */
  ident_needed = ((MidNode!=NULL && NEW_NIDENT(MidNode)) ||
                  (SonEdges[0]!=NULL && NEW_EDIDENT(SonEdges[0])) ||
                  (SonEdges[1]!=NULL && NEW_EDIDENT(SonEdges[1])));

  if (ident_needed)
  {
    /* send number of objects that need identification */
    /* must be equal on all procs                      */
    *((int *)data) = (MidNode!=NULL) + nedges;
    ASSERT(*((int *)data)==1 || *((int *)data)==3);
  }
  else
    *((int *)data) = 0;

  return(0);
}

static int Scatter_NewObjectInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  INT nedges;
  int newsonobjects   = *((int *)data);
  EDGE    *theEdge                = (EDGE *)obj;
  EDGE    *SonEdges[MAX_SON_EDGES];
  NODE    *MidNode        = MIDNODE(theEdge);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  if (GHOST(theEdge)) return(0);

  nedges = GetSonEdges(theEdge,SonEdges);

  if (newsonobjects)
  {
    ASSERT(((newsonobjects==0 || newsonobjects==1) && (MidNode==NULL && nedges!=2))  ||
           (newsonobjects==3 && (MidNode!=NULL && nedges==2)));

    if (MidNode == NULL)
    {
      ASSERT(newsonobjects<=1 && nedges<=1 && SonEdges[1]==NULL);
      if (SonEdges[0] != NULL) SETNEW_EDIDENT(SonEdges[0],1);
    }
    else
    {
      ASSERT(MidNode!=NULL && SonEdges[0]!=NULL && SonEdges[1]!=NULL);
      ASSERT(nedges==2);

      SETNEW_NIDENT(MidNode,1);
      SETNEW_EDIDENT(SonEdges[0],1);
      SETNEW_EDIDENT(SonEdges[1],1);
    }
  }

  return(0);
}

/*************************/

static int Gather_EdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  EDGE *theEdge = (EDGE *)obj;
  NODE *theNode0 = NBNODE(LINK0(theEdge));
  NODE *theNode1 = NBNODE(LINK1(theEdge));

  ASSERT(identlevel == LEVEL(theEdge));

  if (!CORNERTYPE(theNode0) || !CORNERTYPE(theNode1))
  {
    *((int *)data) = 0;
    return(0);
  }

  *((int *)data) = NEW_EDIDENT(theEdge);

  return(0);
}

#define EDIDENTASSERT 1

static int Scatter_EdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  EDGE    *theEdge = (EDGE *)obj;
  int nprop    = *((int *)data);
  NODE *theNode0 = NBNODE(LINK0(theEdge));
  NODE *theNode1 = NBNODE(LINK1(theEdge));

  ASSERT(identlevel == LEVEL(theEdge));

  if (!CORNERTYPE(theNode0) || !CORNERTYPE(theNode1)) return(0);

  if (EDIDENTASSERT)
    if (NEW_EDIDENT(theEdge)) assert(GetFatherEdge(theEdge) != NULL);

  if (nprop)
  {
    if (EDIDENTASSERT) SETNEW_EDIDENT(theEdge,1);
    if (GetFatherEdge(theEdge) == NULL)
    {
      UserWriteF(PFMT "isolated edge=" ID_FMTX "\n",
                 me,ID_PRTX(theEdge));
      if (EDIDENTASSERT) assert(0);
    }
    if (EDIDENTASSERT) assert(GetFatherEdge(theEdge) != NULL);
  }

  return(0);
}

/*************************/

static int Gather_TestEdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  EDGE *theEdge = (EDGE *)obj;

  ASSERT(identlevel == LEVEL(theEdge));

  ((int *)data)[0] = NEW_EDIDENT(theEdge);
  if (NEW_EDIDENT(theEdge)) assert(GetFatherEdge(theEdge) != NULL);

  return(0);
}

static int Scatter_TestEdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  EDGE    *theEdge        = (EDGE *)obj;
  int nprop   = *((int *)data);

  ASSERT(identlevel == LEVEL(theEdge));

  if (NEW_EDIDENT(theEdge) != nprop)
  {
    UserWriteF(PFMT "nprop wrong mynprop=%d hisnprop=%d theEdge=" ID_FMTX " LEVEL=%d PROC=%d PRIO=%d\n",
               me,NEW_EDIDENT(theEdge),nprop,ID_PRTX(theEdge),LEVEL(theEdge),proc,prio);
    fflush(stdout);
    assert(0);
  }

  return(0);
}

/*************************/

static int Gather_IdentSonEdge (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  EDGE *theEdge = (EDGE *)obj;
  EDGE *SonEdge = GetSonEdge(theEdge);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  ((int *)data)[0] = 0;
  ((int *)data)[1] = 0;

  if (SonEdge != NULL)
  {
    ((int *)data)[0] = 1;
    ((int *)data)[1] = NEW_EDIDENT(SonEdge);
  }

  return(0);
}

static int Scatter_IdentSonEdge (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  INT nedges;
  EDGE    *theEdge        = (EDGE *)obj;
  EDGE    *SonEdge;
  int sonedge         = ((int *)data)[0];
  int newsonedge      = ((int *)data)[1];

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  if (SonEdge!=NULL)
  {
    /*
            if (1 || NEW_EDIDENT(SonEdge))
     */
    if (NEW_EDIDENT(SonEdge))
    {
      if(sonedge)
      {
        if (!newsonedge)
        {
          UserWriteF(PFMT "theEdge=" ID_FMTX " LEVEL=%d PROC=%d PRIO=%d sonnprop=%d\n",
                     me,ID_PRTX(theEdge),LEVEL(theEdge),proc,prio,NEW_EDIDENT(SonEdge));
          fflush(stdout);
          assert(0);
        }

        DDD_IdentifyObject(PARHDR(SonEdge),proc,PARHDR(theEdge));
        if (dddctrl.edgeData && EDVECTOR(SonEdge)!=NULL)
          DDD_IdentifyObject(PARHDR(EDVECTOR(SonEdge)),proc,PARHDR(theEdge));
      }
    }
    else
    {
      if (newsonedge)
      {
        UserWriteF(PFMT "theEdge=" ID_FMTX " LEVEL=%d PROC=%d PRIO=%d sonnprop=%d\n",
                   me,ID_PRTX(theEdge),LEVEL(theEdge),proc,prio,NEW_EDIDENT(SonEdge));
        fflush(stdout);
        assert(0);
      }
    }
  }

  return(0);
}

/*************************/

static int Gather_IdentSonObjects (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  INT nedges, ident_needed;
  EDGE *theEdge   = (EDGE *)obj;
  EDGE *SonEdges[2];
  NODE *MidNode   = MIDNODE(theEdge);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  nedges = GetSonEdges(theEdge,SonEdges);

  /* identification is done if one objects of MidNode and the one */
  /* or two sonedge have NEW_XXIDENT flags set.                   */
  ident_needed = ((MidNode!=NULL && NEW_NIDENT(MidNode)) ||
                  (SonEdges[0]!=NULL && NEW_EDIDENT(SonEdges[0])) ||
                  (SonEdges[1]!=NULL && NEW_EDIDENT(SonEdges[1])));

  if (ident_needed)
  {
    /* send number of objects that need identification */
    /* must be equal on all procs                      */
    *((int *)data) = (MidNode!=NULL) + nedges;
    ASSERT(*((int *)data)==1 || *((int *)data)==3);
  }
  else
    *((int *)data) = 0;

  return(0);
}

static int Scatter_IdentSonObjects (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  INT nedges;
  int newsonobjects   = *((int *)data);
  EDGE    *theEdge                = (EDGE *)obj;
  EDGE    *SonEdges[2];
  NODE    *MidNode        = MIDNODE(theEdge);
  NODE    *SonNode0,*SonNode1,*IdentNode;
  NODE    *Node0,*Node1;

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  nedges = GetSonEdges(theEdge,SonEdges);

  if (newsonobjects)
  {
    ASSERT(((newsonobjects==0 || newsonobjects==1) && (MidNode==NULL && nedges!=2))  ||
           (newsonobjects==3 && (MidNode!=NULL && nedges==2)));

    if (MidNode == NULL)
    {
      ASSERT(newsonobjects<=1 && nedges<=1 && SonEdges[1]==NULL);
      if (SonEdges[0] != NULL)
      {
        DDD_IdentifyObject(PARHDR(SonEdges[0]),proc,PARHDR(theEdge));
        if (dddctrl.edgeData && EDVECTOR(SonEdges[0])!=NULL)
          DDD_IdentifyObject(PARHDR(EDVECTOR(SonEdges[0])),proc,PARHDR(theEdge));
      }
    }
    else
    {
      ASSERT(MidNode!=NULL && SonEdges[0]!=NULL && SonEdges[1]!=NULL);
      ASSERT(nedges==2);

      /* identify midnode */
      if (1)
      {
        DDD_IdentifyObject(PARHDR(MidNode),proc,PARHDR(theEdge));
        DDD_IdentifyObject(PARHDRV(MYVERTEX(MidNode)),proc,PARHDR(theEdge));
        if (dddctrl.nodeData && NVECTOR(MidNode)!=NULL)
          DDD_IdentifyObject(PARHDR(NVECTOR(MidNode)),proc,PARHDR(theEdge));
      }
      else
      {
        Node0 = NBNODE(LINK0(theEdge));
        Node1 = NBNODE(LINK1(theEdge));
        DDD_IdentifyObject(PARHDR(MidNode),proc,PARHDR(Node0));
        DDD_IdentifyObject(PARHDR(MidNode),proc,PARHDR(Node1));
        DDD_IdentifyObject(PARHDRV(MYVERTEX(MidNode)),proc,PARHDR(Node0));
        DDD_IdentifyObject(PARHDRV(MYVERTEX(MidNode)),proc,PARHDR(Node1));
        if (dddctrl.nodeData && NVECTOR(MidNode)!=NULL)
        {
          DDD_IdentifyObject(PARHDR(NVECTOR(MidNode)),proc,PARHDR(Node0));
          DDD_IdentifyObject(PARHDR(NVECTOR(MidNode)),proc,PARHDR(Node1));
        }
      }

      /* identify edge0 */
      SonNode0 = NBNODE(LINK0(SonEdges[0]));
      SonNode1 = NBNODE(LINK1(SonEdges[0]));
      if (CORNERTYPE(SonNode0))
      {
        ASSERT(NFATHER(SonNode0)!=NULL);
        IdentNode = SonNode0;
      }
      else
      {
        ASSERT(CORNERTYPE(SonNode1));
        ASSERT(NFATHER(SonNode1)!=NULL);
        IdentNode = SonNode1;
      }
      DDD_IdentifyObject(PARHDR(SonEdges[0]),proc,PARHDR(theEdge));
      DDD_IdentifyObject(PARHDR(SonEdges[0]),proc,PARHDR((NODE *)NFATHER(IdentNode)));
      if (dddctrl.edgeData && EDVECTOR(SonEdges[0])!=NULL)
      {
        DDD_IdentifyObject(PARHDR(EDVECTOR(SonEdges[0])),proc,PARHDR(theEdge));
        DDD_IdentifyObject(PARHDR(EDVECTOR(SonEdges[0])),proc,PARHDR((NODE *)NFATHER(IdentNode)));
      }

      /* identify edge1 */
      SonNode0 = NBNODE(LINK0(SonEdges[1]));
      SonNode1 = NBNODE(LINK1(SonEdges[1]));
      if (CORNERTYPE(SonNode0))
      {
        ASSERT(NFATHER(SonNode0)!=NULL);
        IdentNode = SonNode0;
      }
      else
      {
        ASSERT(CORNERTYPE(SonNode1));
        ASSERT(NFATHER(SonNode1)!=NULL);
        IdentNode = SonNode1;
      }
      DDD_IdentifyObject(PARHDR(SonEdges[1]),proc,PARHDR(theEdge));
      DDD_IdentifyObject(PARHDR(SonEdges[1]),proc,PARHDR((NODE *)NFATHER(IdentNode)));
      if (dddctrl.edgeData && EDVECTOR(SonEdges[1])!=NULL)
      {
        DDD_IdentifyObject(PARHDR(EDVECTOR(SonEdges[1])),proc,PARHDR(theEdge));
        DDD_IdentifyObject(PARHDR(EDVECTOR(SonEdges[1])),proc,PARHDR((NODE *)NFATHER(IdentNode)));
      }
    }
  }

  return(0);
}

#endif
#endif



#ifndef IDENT_ONLY_NEW

/****************************************************************************/
/*
   Gather_SonNodeInfo -

   SYNOPSIS:
   static int Gather_SonNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

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

static int Gather_SonNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE *theNode = (NODE *)obj;

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theNode));

  if (SONNODE(theNode) != NULL)
    *((int *)data) = 1;
  else
    *((int *)data) = 0;

  return(0);
}

/****************************************************************************/
/*
   Scatter_SonNodeInfo -

   SYNOPSIS:
   static int Scatter_SonNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

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

static int Scatter_SonNodeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  NODE    *theNode        = (NODE *)obj;
  NODE    *SonNode        = SONNODE(theNode);
  INT has_sonnode     = *((int *)data);

  /* identification is only done between master objects */
  ASSERT(identlevel-1 == LEVEL(theNode));

  if (SonNode != NULL)
  {
    if (has_sonnode)
    {
      DDD_IdentifyObject(PARHDR(SonNode),proc,PARHDR(theNode));
      if (dddctrl.nodeData && NVECTOR(SonNode)!=NULL)
        DDD_IdentifyObject(PARHDR(NVECTOR(SonNode)),proc,PARHDR(theNode));
      IFDEBUG(dddif,1)
      if (dddctrl.nodeData && NVECTOR(SonNode)!=NULL)
        PrintDebug ("l=%d IdentHdr: %d Proc: %d me:%d IdentObjectHdr: %d %d\n",
                    identlevel,GID(theNode),proc,me,GID(SonNode),GID(EDVECTOR(SonNode)));
      else
        PrintDebug ("l=%d IdentHdr: %d Proc: %d me:%d IdentObjectHdr: %d\n",
                    identlevel,GID(theNode),proc,me,GID(SonNode));
      ENDDEBUG
    }
  }

  return(0);
}


#ifdef __THREEDIM__

/****************************************************************************/
/*
   Gather_SonEdgeInfo  -

   SYNOPSIS:
   static int Gather_SonEdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

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

static int Gather_SonEdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  EDGE *theEdge = (EDGE *)obj;
  EDGE *SonEdge;

  /* identification has to be done between all copies of an objects */
  /* otherwise this can result in unsymmetric interfaces            */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  SonEdge = GetSonEdge(theEdge);
  if (SonEdge != NULL)
    *((int *)data) = 1;
  else
    *((int *)data) = 0;

  return(0);
}


/****************************************************************************/
/*
   Scatter_SonEdgeInfo -

   SYNOPSIS:
   static int Scatter_SonEdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio);

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

static int Scatter_SonEdgeInfo (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  EDGE    *theEdge        = (EDGE *)obj;
  EDGE    *SonEdge;
  INT has_sonedge     = *((int *)data);

  /* identification has to be done between all copies of an objects */
  /* otherwise this can result in unsymmetric interfaces            */
  ASSERT(identlevel-1 == LEVEL(theEdge));

  SonEdge = GetSonEdge(theEdge);
  if (SonEdge != NULL)
  {
    if (has_sonedge)
    {
      DDD_IdentifyObject(PARHDR(SonEdge),proc,PARHDR(theEdge));
      if (dddctrl.edgeData && EDVECTOR(SonEdge)!=NULL)
        DDD_IdentifyObject(PARHDR(EDVECTOR(SonEdge)),proc,PARHDR(theEdge));
      IFDEBUG(dddif,1)
      if (dddctrl.edgeData && EDVECTOR(SonEdge)!=NULL)
        PrintDebug ("l=%d IdentHdr: %d Proc: %d me:%d IdentObjectHdr: %d %d\n",
                    identlevel,GID(theEdge),proc,me,GID(SonEdge),GID(EDVECTOR(SonEdge)));
      else
        PrintDebug ("l=%d IdentHdr: %d Proc: %d me:%d IdentObjectHdr: %d\n",
                    identlevel,GID(theEdge),proc,me,GID(SonEdge));
      ENDDEBUG
    }
  }

  return(0);
}
#endif
#endif

/****************************************************************************/
/*
   Identify_SonNodes - identify son nodes (type CORNERNODE)

   SYNOPSIS:
   INT Identify_SonNodes (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT Identify_SonNodes (GRID *theGrid)
{

#ifdef IDENT_ONLY_NEW
  DDD_IFAOnewayX(NodeAllIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(int),
                 Gather_NewNodeInfo,Scatter_NewNodeInfo);

  if (UPGRID(theGrid) != NULL)
  {
    DDD_IFAOnewayX(NodeAllIF,GRID_ATTR(UPGRID(theGrid)),IF_FORWARD,sizeof(int),
                   Gather_NodeInfo,Scatter_NodeInfo);
    if (0)
      DDD_IFAOnewayX(NodeAllIF,GRID_ATTR(UPGRID(theGrid)),IF_FORWARD,sizeof(int),
                     Gather_TestNodeInfo,Scatter_TestNodeInfo);
  }

  DDD_IFAOnewayX(NodeAllIF,GRID_ATTR(theGrid),IF_FORWARD,2*sizeof(int),
                 Gather_IdentSonNode,Scatter_IdentSonNode);

#else

  DDD_IFAOnewayX(NodeAllIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(int),
                 Gather_SonNodeInfo,Scatter_SonNodeInfo);

#endif
  return(GM_OK);
}

#ifdef __THREEDIM__

/****************************************************************************/
/*
   Identify_SonEdges - identify son edges

   SYNOPSIS:
   INT Identify_SonEdges (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT Identify_SonEdges (GRID *theGrid)
{

#ifdef IDENT_ONLY_NEW

  DDD_IFAOnewayX(BorderEdgeSymmIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(int),
                 Gather_NewObjectInfo,Scatter_NewObjectInfo);

  if (0)
    if (UPGRID(theGrid) != NULL)
    {
      DDD_IFAOnewayX(EdgeSymmVHIF,GRID_ATTR(UPGRID(theGrid)),IF_FORWARD,sizeof(int),
                     Gather_EdgeInfo,Scatter_EdgeInfo);
      DDD_IFAOnewayX(EdgeSymmVHIF,GRID_ATTR(UPGRID(theGrid)),IF_FORWARD,sizeof(int),
                     Gather_TestEdgeInfo,Scatter_TestEdgeInfo);
    }

  DDD_IFAOnewayX(BorderEdgeSymmIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(int),
                 Gather_IdentSonObjects,Scatter_IdentSonObjects);

#else

  /* identify the sonedges */
  DDD_IFAOnewayX(EdgeSymmVHIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(int),
                 Gather_SonEdgeInfo,Scatter_SonEdgeInfo);

#endif

  return(GM_OK);
}
#endif


/****************************************************************************/
/*
   Identify_SonObjects - identify son objects

   SYNOPSIS:
   INT Identify_SonObjects (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:
   This function identifies all objects which are not symmetrically created
   during grid adaption. These are edges and nodes of the type used by
   yellow elements, son nodes of type CORNERNODE and son edges.

   RETURN VALUE:
   INT
 */
/****************************************************************************/

#define NODESFIRST 1

INT Identify_SonObjects (GRID *theGrid)
{
    #ifdef Debug
  identlevel = GLEVEL(theGrid)+1;
        #endif

  if (NODESFIRST)
  {
    if (Identify_SonNodes (theGrid) != GM_OK) RETURN(GM_ERROR);
  }
  else
  {
        #ifdef __THREEDIM__
    if (Identify_SonEdges (theGrid) != GM_OK) RETURN(GM_ERROR);
        #endif
  }

  if (1)
  {
    printf(PFMT " 2. DDD_IdentifyEnd() in Identify_SonObjects()\n",me);
    fflush(stdout);
    Synchronize();
    DDD_IdentifyEnd();
    DDD_IdentifyBegin();
  }

  if (NODESFIRST)
  {
        #ifdef __THREEDIM__
    if (Identify_SonEdges (theGrid) != GM_OK) RETURN(GM_ERROR);
        #endif
  }
  else
  {
    if (Identify_SonNodes (theGrid) != GM_OK) RETURN(GM_ERROR);
  }

  return (GM_OK);
}

/****************************************************************************/
/*
   Identify_Objects_of_ElementSide -

   SYNOPSIS:
   INT Identify_Objects_of_ElementSide(GRID *theGrid, ELEMENT *theElement, INT i);

   PARAMETERS:
   .  theGrid
   .  theElement
   .  i

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT Identify_Objects_of_ElementSide(GRID *theGrid, ELEMENT *theElement, INT i)
{
  INT prio;
  ELEMENT *theNeighbor;

  theNeighbor = NBELEM(theElement,i);
  if (theNeighbor == NULL) return(GM_OK);

  prio = EPRIO(theNeighbor);
  /* identification is only needed if theNeighbor removed his refinement  */
  /* or was not refined before, thus has NSONS==0, if NSONS>0 the objects */
  /* shared between the element sides are already identified and no new   */
  /* objects are created for this element side which need identification  */
  /* (980217 s.l.)                                                        */
  /*
          if (!EHGHOSTPRIO(prio) || NSONS(theNeighbor)!=0) return(GM_OK);
   */
  if (!EHGHOSTPRIO(prio) || !MARKED(theNeighbor)) return(GM_OK);

        #ifdef Debug
  identlevel = GLEVEL(theGrid);
        #endif
  if (IdentifyObjectsOfElementSide(theGrid,theElement,i,theNeighbor)) RETURN(GM_FATAL);

  return(GM_OK);
}


/****************************************************************************/
/*
    IdentifyInit-

   SYNOPSIS:
   void IdentifyInit (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void IdentifyInit (MULTIGRID *theMG)
{
  INT i;

        #ifdef Debug
  debug = 0;
        #endif

  /* allocate a control word entry to lock nodes */
  if (AllocateControlEntry(NODE_CW,NEW_NIDENT_LEN,&ce_NEW_NIDENT) != GM_OK)
    assert(0);

  /* allocate a control word entry to lock edges */
  if (AllocateControlEntry(EDGE_CW,NEW_EDIDENT_LEN,&ce_NEW_EDIDENT) != GM_OK)
    assert(0);

  for (i=0; i<=TOPLEVEL(theMG); i++)
    ResetIdentFlags(GRID_ON_LEVEL(theMG,i));

  /* set Ident_FctPtr to identification mode */
  Ident_FctPtr = Identify_by_ObjectList;

}


/****************************************************************************/
/*
   IdentifyExit -

   SYNOPSIS:
   void IdentifyExit (void);

   PARAMETERS:
   .  void

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void IdentifyExit (void)
{
  FreeControlEntry(ce_NEW_NIDENT);
  FreeControlEntry(ce_NEW_EDIDENT);
}

#endif /* end ModelP */
