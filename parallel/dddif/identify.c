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
#define GET_IDENT_HDR(node) ( (NTYPE(node) == CORNER_NODE) ? \
                              PARHDR(NFATHER(node)) : PARHDR(node) )

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

/* temp node flag for Identification */
static INT ce_NIDENT;
#define NIDENT_LEN                    1
#define NIDENT(p)                     CW_READ(p,ce_NIDENT)
#define SETNIDENT(p,n)                CW_WRITE(p,ce_NIDENT,n)

#ifdef __THREEDIM__
/* temp edge flag for Identification */
static INT ce_EDIDENT;
#define EDIDENT_LEN                   1
#define EDIDENT(p)                    CW_READ(p,ce_EDIDENT)
#define SETEDIDENT(p,n)               CW_WRITE(p,ce_EDIDENT,n)
#endif

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

INT compare_gid (const void *e0, const void *e1)
{
  INT num0, num1;

  num0 = DDD_InfoGlobalId(*(DDD_HDR *)e0);
  num1 = DDD_InfoGlobalId(*(DDD_HDR *)e1);

  if (num0 < num1) return(1);
  if (num0 > num1) return(-1);
  return(0);
}

#ifdef Debug
static void ResetIdentFlags (GRID *UpGrid)
{
  NODE *theNode;
  EDGE *theEdge;
  LINK *theLink;

  /* clear all IDENT flags */
  for (theNode=FIRSTNODE(UpGrid); theNode!=NULL; theNode=SUCCN(theNode)) {

    SETNIDENT(theNode,CLEAR);

                #ifdef __THREEDIM__
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink)) {
      theEdge = MYEDGE(theLink);
      SETEDIDENT(theEdge,CLEAR);
    }
                #endif
  }

}
#endif

#ifdef Debug
static INT Print_Identify_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject, int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
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
  for (i=0; i<nident; i++) {
    PrintDebug (" %d",DDD_InfoGlobalId(IdentHdr[i]));
  }

  /* print type of objects to identify */
  PrintDebug ("    IdentObjectType:");
  for (i=0; i<nobject; i++) {
    PrintDebug (" %d",DDD_InfoType(IdentObjectHdr[i]));
  }

  /* print the proclist to identify to */
  PrintDebug ("    ProcList: %d",me);
  while (*proclist != -1) {
    if (*(proclist+1) == skiptag) {
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
  for (i=0; i<nobject; i++) {
    PrintDebug (" %d",DDD_InfoGlobalId(IdentObjectHdr[i]));
  }

  PrintDebug ("\n");

  return;

}
#endif

#ifdef Debug
static INT Print_Identified_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject, int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
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
  for (i=0; i<nobject; i++) {
    PrintDebug (" %d",DDD_InfoGlobalId(IdentObjectHdr[i]));
  }

  /* print the objects used for identify */
  PrintDebug ("    IdentHdr:");
  for (i=0; i<nident; i++) {
    PrintDebug (" %d",DDD_InfoGlobalId(IdentHdr[i]));
  }

  /* print the proclist to identify to */
  PrintDebug ("    ProcList: %d",me);
  while (*proclist != -1) {
    if (*(proclist+1) == skiptag) {
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
  for (i=0; i<nobject; i++) {
    PrintDebug (" %d",DDD_InfoType(IdentObjectHdr[i]));
  }
  PrintDebug ("\n");

  return;

}
#endif

static INT Identify_by_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject, int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
{
  INT i,j,n;

  ASSERT(nobject>0);
  ASSERT(nident>0);
  ASSERT(*proclist!=-1);

  IFDEBUG(dddif,5)
  Print_Identify_ObjectList(IdentObjectHdr,nobject,proclist,skiptag,IdentHdr,nident);
  ENDDEBUG

    n = 0;
  while (*proclist != -1) {
    ASSERT(n<procs);

    if (*(proclist+1) == skiptag) {
      proclist += 2;
      continue;
    }

    /* identify the object */
    for (j=0; j<nobject; j++) {
      for (i=0; i<nident; i++) {

        PRINTDEBUG(dddif,5,("%d: Identify_by_ObjectList(): Type=%d" \
                            " IdentObjectHdr=%08x proclist=%d IdentHdr=%08x me=%d\n",
                            me,
                            DDD_InfoType(IdentObjectHdr[j]),
                            DDD_InfoGlobalId(IdentObjectHdr[j]),
                            *proclist,
                            DDD_InfoGlobalId(IdentHdr[i]),
                            me));

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
static INT IdentifySideVector (ELEMENT* theElement, ELEMENT *theNeighbor, ELEMENT *Son, INT SonSide)
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
  for (k=0; k<CORNERS_OF_SIDE(Son,SonSide); k++) {
    theNode = CORNER(Son,CORNER_OF_SIDE(Son,SonSide,k));
    if (NFATHER(theNode) != NULL)
      IdentHdr[nident++] = PARHDR(NFATHER(theNode));
    else
      IdentHdr[nident++] = PARHDR(theNode);
  }

  proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

  Ident_FctPtr(IdentObjectHdr, 1,
               proclist+2, PrioGhost, IdentHdr, nident);

}
#endif

static void IdentifyNode (ELEMENT *theNeighbor, NODE *theNode, NODE *Nodes[MAX_SIDE_NODES], INT node, INT ncorners, INT Vec)
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
  /* return if already identified */
  if (NIDENT(theNode) == IDENT) return;

  switch (NTYPE(theNode)) {
    int *proclist;

  case (CORNER_NODE) :

    PRINTDEBUG(dddif,1,("%d: Identify CORNERNODE gid=%08x node=%d vec=%d\n",
                        me, DDD_InfoGlobalId(PARHDR(theNode)), node, Vec));

    IdentObjectHdr[nobject++] = PARHDR(theNode);
    if (Vec)
      IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

    /* identify to proclist of node */
    proclist = DDD_InfoProcList(PARHDR(NFATHER(theNode)));

    /* identify using father node */
    IdentHdr[nident++] = PARHDR(NFATHER(theNode));

    Ident_FctPtr(IdentObjectHdr, nobject,
                 proclist+2, PrioGhost, IdentHdr, nident);

    break;

  case (MID_NODE) : {

                        #ifdef __TWODIM__
    NODE **EdgeNodes;
    EdgeNodes = Nodes;
                        #endif

                        #ifdef __THREEDIM__
    NODE *EdgeNodes[MAX_SIDE_NODES];
    EDGE *theEdge;

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
    /* 2D: identify to proclist of neighbor element */
    proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
                        #endif

                        #ifdef __THREEDIM__
    /* 3D: identify to proclist of edge */
    theEdge = GetEdge(NFATHER(EdgeNodes[0]),NFATHER(EdgeNodes[1]));
    ASSERT(theEdge!=NULL);

    proclist = DDD_InfoProcList(PARHDR(theEdge));
                        #endif

    /* identify using edge nodes */
    IdentHdr[nident++] = PARHDR(NFATHER(EdgeNodes[0]));
    IdentHdr[nident++] = PARHDR(NFATHER(EdgeNodes[1]));

    Ident_FctPtr(IdentObjectHdr, nobject,
                 proclist+2, PrioGhost, IdentHdr, nident);

    break;
  }

                #ifdef __THREEDIM__
  case (SIDE_NODE) : {

    INT i;

    PRINTDEBUG(dddif,1,("%d: Identify SIDENODE gid=%08x node=%d Vec=%d\n",
                        me, DDD_InfoGlobalId(PARHDR(theNode)), node, Vec));

    /* identify sidenode, vertex and vector */
    IdentObjectHdr[nobject++] = PARHDR(theNode);
    IdentObjectHdr[nobject++] = PARHDRV(MYVERTEX(theNode));
    if (Vec)
      IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

    /* identify to proclist of neighbor element */
    proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

    /* identify using corner nodes of side */
    for (i=0; i<ncorners; i++)
      IdentHdr[nident++] = PARHDR(NFATHER(Nodes[i]));

    /* identify side node */
    Ident_FctPtr(IdentObjectHdr, nobject,
                 proclist+2, PrioGhost, IdentHdr, nident);

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

#ifdef __THREEDIM__
EDGE *FatherEdge (NODE **SideNodes, INT ncorners, NODE **Nodes, EDGE *theEdge)
{
  INT sonedge,pos0,pos1;
  EDGE *fatherEdge = NULL;

  ASSERT(Nodes[0]!=NULL);
  ASSERT(Nodes[1]!=NULL);

  /* one node is side node -> no father edge */
  if (NTYPE(Nodes[0])==SIDE_NODE || NTYPE(Nodes[1])==SIDE_NODE) return(NULL);

  /* both nodes are side nodes -> no father edge */
  if (NTYPE(Nodes[0])==MID_NODE && NTYPE(Nodes[1])==MID_NODE) return(NULL);

  for (pos0=0; pos0<MAX_SIDE_NODES; pos0++) {
    if (SideNodes[pos0] == Nodes[0])
      break;
  }
  ASSERT(pos0<MAX_SIDE_NODES);

  for (pos1=0; pos1<MAX_SIDE_NODES; pos1++) {
    if (SideNodes[pos1] == Nodes[1])
      break;
  }
  ASSERT(pos1<MAX_SIDE_NODES);

  switch (NTYPE(Nodes[0])) {

  case (CORNER_NODE) :

    ASSERT(pos0<ncorners);
    if ( (pos0+1 == pos1) ||
         (pos0+ncorners == pos1) ) {

      fatherEdge = GetEdge(NFATHER(Nodes[0]),NFATHER(SideNodes[pos0+1]));
      ASSERT(fatherEdge!=NULL);
    }

    if ( ((pos0-1+ncorners)%ncorners == pos1) ||
         ((pos0-1+ncorners)%ncorners+ncorners == pos1) ) {

      fatherEdge = GetEdge(NFATHER(Nodes[0]),
                           NFATHER(SideNodes[(pos0-1+ncorners)%ncorners]));
      ASSERT(fatherEdge!=NULL);
    }

    break;

  case (MID_NODE) :

    ASSERT(pos0>=ncorners);
    ASSERT(pos0<2*ncorners);

    if ((pos0+1)%ncorners == pos1) {

      fatherEdge = GetEdge(NFATHER(SideNodes[pos0%ncorners]),
                           NFATHER(Nodes[1]));
      ASSERT(fatherEdge!=NULL);
    }

    if (pos0%ncorners == pos1) {

      fatherEdge = GetEdge(NFATHER(SideNodes[(pos0+1)%ncorners]),
                           NFATHER(Nodes[1]));
      ASSERT(fatherEdge!=NULL);
    }

    break;

  case (SIDE_NODE) :

    /* this edge has no father edge */
    fatherEdge = NULL;
    break;

  default :
    assert(0);
    break;
  }

  IFDEBUG(dddif,0)
  INT i;
  EDGE* edge0, *edge1;

  edge0 = edge1 = NULL;

  /* test whether theEdge lies above fatherEdge */
  if (MIDNODE(fatherEdge) != NULL) {
    edge0 = GetEdge(MIDNODE(fatherEdge),SONNODE(NBNODE(LINK0(fatherEdge))));
    edge1 = GetEdge(MIDNODE(fatherEdge),SONNODE(NBNODE(LINK1(fatherEdge))));
  }
  else {
    edge0 = GetEdge(SONNODE(NBNODE(LINK0(fatherEdge))),SONNODE(NBNODE(LINK1(fatherEdge))));
  }

  IFDEBUG(dddif,5)
  UserWriteF("%4d: fatherEdge=%x theEdge=%x edge0=%x edge1=%x\n",me,fatherEdge,theEdge,edge0,edge1);
  UserWriteF("%4d: Nodes[0]=%d Nodes[1]=%d\n",me,ID(Nodes[0]),ID(Nodes[1]));

  UserWriteF("SideNodes\n");
  for (i=0; i<MAX_SIDE_NODES; i++) UserWriteF(" %5d",i);
  UserWriteF("\n");
  for (i=0; i<MAX_SIDE_NODES; i++)
    if (SideNodes[i]!=NULL) UserWriteF(" %5d",ID(SideNodes[i]));
  UserWriteF("\n");
  ENDDEBUG

  assert(edge0==theEdge || edge1==theEdge);
  ENDDEBUG

  return(fatherEdge);
}
#endif

static INT IdentifyEdge (ELEMENT *theElement, ELEMENT *theNeighbor, NODE **SideNodes, INT ncorners, ELEMENT *Son, INT SonSide, INT edgeofside, INT Vec)
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
  PRINTDEBUG(dddif,0,("%4d: edge=%d corner0=%d corner1=%d Nodes[0]=%d Nodes[1]=%d\n",
                      me,edge, corner0, corner1, ID(Nodes[0]), ID(Nodes[1])));
        #endif

  ASSERT(Nodes[0]!=NULL);
  ASSERT(Nodes[1]!=NULL);

  theEdge = GetEdge(Nodes[0],Nodes[1]);
  ASSERT(theEdge!=NULL);

  /* edge unlocked -> no debugging occurs */
        #ifdef Debug
  if (debug == 1) {
    if (EDIDENT(theEdge) == CLEAR) return;
  }
  else
        #endif
  /* edge locked -> already identified */
  if (EDIDENT(theEdge) == IDENT) return;

  PRINTDEBUG(dddif,1,("%d: Identify EDGE edgeofside=%d pe=%08x/%x eID=%d"
                      " ntype0=%d  ntype1=%d\n",me,edgeofside,
                      DDD_InfoGlobalId(PARHDRE(Son)),Son,ID(Son),
                      NTYPE(Nodes[0]), NTYPE(Nodes[1])))

        #ifdef __THREEDIM__
  IdentObjectHdr[nobject++] = PARHDR(theEdge);
        #endif
  if (Vec)
    IdentObjectHdr[nobject++] = PARHDR(EDVECTOR(theEdge));

        #ifdef __TWODIM__
  /* identify to proclist of neighbor */
  proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
        #endif

        #ifdef __THREEDIM__
  /* identify to proclist of father edge or neighbor*/
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

  if (NFATHER(Nodes[0]) != NULL)
    IdentHdr[nident++] = PARHDR(NFATHER(Nodes[0]));
  else
    IdentHdr[nident++] = PARHDR(Nodes[0]);

  if (NFATHER(Nodes[1]) != NULL)
    IdentHdr[nident++] = PARHDR(NFATHER(Nodes[1]));
  else
    IdentHdr[nident++] = PARHDR(Nodes[1]);

  Ident_FctPtr(IdentObjectHdr, nobject,
               proclist+2, PrioGhost, IdentHdr, nident);

  /* debugging unlocks the edge */
        #ifdef Debug
  if (debug == 1) {
    SETEDIDENT(theEdge,CLEAR);
  }
  else
        #endif
  /* lock this edge for identification */
  SETEDIDENT(theEdge,IDENT);

  return;
}


static INT IdentifyObjectsOfElementSide(GRID *theGrid, ELEMENT *theElement, INT i, ELEMENT *theNeighbor)
{
  INT nodes,j;
  NODE *SideNodes[MAX_SIDE_NODES];
  INT ncorners;
  NODE *theNode;

  GetSonSideNodes(theElement,i,&nodes,SideNodes);
  ncorners = CORNERS_OF_SIDE(theElement,i);

  PRINTDEBUG(dddif,1,("%d: IdentifyObjectsOfElementSide():identify NODES ncorners=%d nodes=%d\n",me,
                      ncorners,nodes));

  /* identify nodes, vertices and node vectors of son elements */
  for (j=0; j<nodes; j++) {
    INT prio;

    theNode = SideNodes[j];
    ASSERT(theNode != NULL);

    /* identify new node including its vector and vertex        */
    IdentifyNode(theNeighbor, theNode, SideNodes, j, ncorners,
                 TYPE_DEF_IN_GRID(theGrid,NODEVECTOR));
  }

  /* identify edge vectors (2D); edges, edge and side vectors (3D) */
  if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR) || DIM==3) {

    ELEMENT *SonList[MAX_SONS];
    INT SonsOfSide,SonSides[MAX_SONS];
    INT j;

    PRINTDEBUG(dddif,1,("%d: IdentifyObjectsOfElementSide(): identify EDGES and VECTORS\n",me));

    if (Get_Sons_of_ElementSide(theElement,i,&SonsOfSide,
                                SonList,SonSides,1)!=GM_OK)
      RETURN(GM_FATAL);

    for (j=0; j<SonsOfSide; j++) {

      if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR) || DIM==3) {

        INT edgeofside;
        INT nedges = EDGES_OF_SIDE(SonList[j],SonSides[j]);

        /* identify the edge and vector */
        for (edgeofside=0; edgeofside<nedges; edgeofside++) {
          IdentifyEdge(theElement,theNeighbor,SideNodes,ncorners,
                       SonList[j],SonSides[j],edgeofside,
                       TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR));
        }
      }

                        #ifdef __THREEDIM__
      if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR)) {

        IdentifySideVector(theElement,theNeighbor,SonList[j],SonSides[j]);
      }
                        #endif
    }
  }
}

INT     IdentifyDistributedObjects (MULTIGRID *theMG, INT FromLevel, INT ToLevel)
{
  INT l,i,j,prio;
  ELEMENT *theElement,*theNeighbor;
  NODE *theNode;
  GRID *theGrid;

  PRINTDEBUG(dddif,1,("%d: IdentifyDistributedObjects(): FromLevel=%d ToLevel=%d\n",
                      me,FromLevel,ToLevel));

  /* identify distributed objects */
  for (l=FromLevel; l<ToLevel; l++) {

    PRINTDEBUG(dddif,1,("%d: IdentifyDistributedObjects(): identification level=%d\n",me,l));

    theGrid = GRID_ON_LEVEL(theMG,l);

                #ifdef Debug
    identlevel = l;
                #endif

    /* check control word flags for ident on upper level */
                #ifdef Debug
    if (debug != 1)
      ResetIdentFlags(GRID_ON_LEVEL(theMG,l+1));
                #endif

    for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement)) {

      if (!IS_REFINED(theElement) ||
          (prio = DDD_InfoPriority(PARHDRE(theElement))) == PrioGhost) {
        continue;
      }

      for (i=0; i<SIDES_OF_ELEM(theElement); i++) {

        theNeighbor = NBELEM(theElement,i);
        if (theNeighbor == NULL) continue;

        if ((prio = DDD_InfoPriority(PARHDRE(theNeighbor))) != PrioGhost
            || NSONS(theNeighbor)!=0)
          continue;

        PRINTDEBUG(dddif,1,("%d: Identify element: pe=%08x/%x eID=%d side=%d\n",me,
                            DDD_InfoGlobalId(PARHDRE(theElement)),theElement,
                            ID(theElement),i));

        IdentifyObjectsOfElementSide(theGrid,theElement,i,theNeighbor);
      }
    }
  }

  return(GM_OK);
}


INT     IdentifyGridLevels (MULTIGRID *theMG, INT FromLevel, INT ToLevel)
{
        #ifdef Debug
  debug = 0;
        #endif

  /* allocate a control word entry to lock nodes */
  if (AllocateControlEntry(NODE_CW,1,&ce_NIDENT) != GM_OK)
    assert(0);

        #ifdef __THREEDIM__
  /* allocate a control word entry to lock edges */
  if (AllocateControlEntry(EDGE_CW,1,&ce_EDIDENT) != GM_OK)
    assert(0);
        #endif

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

  FreeControlEntry(ce_NIDENT);

        #ifdef __THREEDIM__
  FreeControlEntry(ce_EDIDENT);
        #endif
}

#endif /* end ModelP */
