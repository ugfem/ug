// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugm.c                                                                                                                 */
/*																			*/
/* Purpose:   unstructured grid manager                                                                         */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de	                                        */
/*																			*/
/* History:   09.03.92 begin, ug version 2.0								*/
/*			  Aug 28 1996, ug version 3.4                                                                   */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment ugm
#endif

/****************************************************************************/
/*																			*/
/*		defines to exclude functions										*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>

#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "debug.h"
#include "general.h"

#include "devices.h"

#include "switch.h"
#include "evm.h"
#include "gm.h"
#include "misc.h"
#include "dlmgr.h"
#include "algebra.h"
#include "ugm.h"
#include "elements.h"
#include "shapes.h"
#include "refine.h"
#include "domain.h"

#ifdef ModelP
#include "parallel.h"
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

#define RESOLUTION       20     /* resolution for creating boundary midnode */
#define SMALL1 0.001

#define ORDERRES                1e-3    /* resolution for OrderNodesInGrid			*/
#define LINKTABLESIZE   32              /* max number of inks per node for ordering	*/

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

static char buffer[4*256];                      /* general purpose text buffer			*/

static VIRT_HEAP_MGMT *theGenMGUDM; /* general user data space management	*/

static INT theMGDirID;                          /* env var ID for the multigrids		*/
static INT theMGRootDirID;                      /* env dir ID for the multigrids		*/

static INT UsedOBJT;                            /* for the dynamic OBJECT management	*/

/* used by OrderNodesInGrid */
static const INT *Order,*Sign;
static DOUBLE InvMeshSize;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static NODE *CreateNode (GRID *theGrid);
static VERTEX *CreateBoundaryVertex     (GRID *theGrid);
static VERTEX *CreateInnerVertex (GRID *theGrid);

static INT DisposeNode (GRID *theGrid, NODE *theNode);
static INT DisposeVertex (GRID *theGrid, VERTEX *theVertex);
static INT DisposeEdge (GRID *theGrid, EDGE *theEdge);

/****************************************************************************/
/*D
   GetFreeOBJT - Get an object type id not occupied in theMG

   SYNOPSIS:
   INT GetFreeOBJT ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function gets an object type id not occupied in theMG.

   RETURN VALUE:
   INT
   .n    id of object type if ok
   .n   -1 when error occured.
   D*/
/****************************************************************************/

INT GetFreeOBJT ()
{
  INT i;

  for (i=0; i<MAXOBJECTS; i++)
    if (!READ_FLAG(UsedOBJT,1<<i))
      break;

  if (i<MAXOBJECTS)
  {
    SET_FLAG(UsedOBJT,1<<i);
    return (i);
  }
  else
    return (-1);
}

/****************************************************************************/
/*D
   ReleaseOBJT - Release an object type id not needed anymore

   SYNOPSIS:
   INT ReleaseOBJT (INT type);

   PARAMETERS:
   .  type - object type

   DESCRIPTION:
   This function releases an object type id not needed anymore.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT ReleaseOBJT (INT type)
{
  if (type>=MAXOBJECTS)
    RETURN (GM_ERROR);

  CLEAR_FLAG(UsedOBJT,1<<type);

  return (GM_OK);
}

/****************************************************************************/
/*D
   GetMemoryForObject - Get an object from free list if possible

   SYNOPSIS:
   void *GetMemoryForObject (MULTIGRID *theMG, INT size, INT type);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  size - size of the object
   .  type - type of the requested object

   DESCRIPTION:
   This function gets an object of type `type` from free list if possible,
   otherwise it allocates memory from the multigrid heap using 'GetMem'.

   RETURN VALUE:
   void *
   .n   pointer to an object of the requested type
   .n   NULL if object of requested type is not available
   D*/
/****************************************************************************/

#ifdef ModelP
void *GetMemoryForObject_par (HEAP *theHeap, INT size, INT type)
{
  void *obj = GetFreelistMemory(theHeap, size);

  if (obj!=NULL && type!=NOOBJ)
  {
    memset(obj,0,size);
    /* link this object to DDD management */
    if (HAS_DDDHDR(type))
    {
      DDD_TYPE dddtype = DDDTYPE(type);
      DDD_HDR dddhdr = (DDD_HDR)(((char *)obj) + DDD_InfoHdrOffset(dddtype));
      DDD_HdrConstructor(dddhdr, dddtype, PrioNone, 0);
    }
  }

  return obj;
}
#endif

/****************************************************************************/
/*D
   PutFreeObject - Put an object in the free list

   SYNOPSIS:
   INT PutFreeObject (MULTIGRID *theMG, void *object, INT size, INT type);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  object - object to insert in free list
   .  size - size of the object
   .  type - type of the requested object

   DESCRIPTION:
   This function puts an object in the free list.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 when error occured.
   D*/
/****************************************************************************/

#ifdef ModelP
INT PutFreeObject_par (HEAP *theHeap, void *object, INT size, INT type)
{
  if (type!=NOOBJ)
  {
    /* unlink object from DDD management */
    if (HAS_DDDHDR(type))
    {
      DDD_HDR dddhdr = (DDD_HDR)
                       (((char *)object)+DDD_InfoHdrOffset(DDDTYPE(type)));
      DDD_HdrDestructor(dddhdr);
    }
  }

  return (PutFreelistMemory(theHeap, object, size));
}
#endif

/****************************************************************************/
/*D
   CreateBoundaryVertex - Return pointer to a new boundary vertex structure

   SYNOPSIS:
   static VERTEX *CreateBoundaryVertex (GRID *theGrid);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  after - vertex after which to insert

   DESCRIPTION:
   This function creates and initializes a new boundary vertex structure
        and returns a pointer to it.

   RETURN VALUE:
   VERTEX *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

static VERTEX *CreateBoundaryVertex (GRID *theGrid)
{
  VERTEX *pv;
  INT ds;
  INT i;

  pv = GetMemoryForObject(MYMG(theGrid),sizeof(struct bvertex),BVOBJ);
  if (pv==NULL) return(NULL);
  if ((ds=theGrid->mg->theFormat->sVertex)>0)
  {
    VDATA(pv) = GetMemoryForObject(MYMG(theGrid),ds,NOOBJ);
    if (VDATA(pv)==NULL) return(NULL);
  }
  else
    VDATA(pv) = NULL;

  /* initialize data */
  CTRL(pv) = 0;
  SETOBJT(pv,BVOBJ);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
  TOPNODE(pv) = NULL;
  for (i=0; i<DIM; i++) LCVECT(pv)[i] = 0.0;
  SETONEDGE(pv,0);
  SETMOVE(pv,DIM_OF_BND);
        #ifdef ModelP
  DDD_AttrSet(PARHDRV(pv),theGrid->level);
  DDD_PrioritySet(PARHDRV(pv),PrioMaster);
        #endif

  /* insert in vertex list */
  GRID_LINK_VERTEX(theGrid,pv,PrioMaster);

  return(pv);
}

/****************************************************************************/
/*D
   CreateInnerVertex - Return pointer to a new inner vertex structure

   SYNOPSIS:
   static VERTEX *CreateInnerVertex (GRID *theGrid);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  after - vertex after which to insert

   DESCRIPTION:
   This function creates and initializes a new inner vertex structure
   and returns a pointer to it.

   RETURN VALUE:
   VERTEX *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

static VERTEX *CreateInnerVertex (GRID *theGrid)
{
  VERTEX *pv;
  INT ds;
  INT i;

  pv = GetMemoryForObject(MYMG(theGrid),sizeof(struct ivertex),IVOBJ);
  if (pv==NULL) return(NULL);
  if ((ds=theGrid->mg->theFormat->sVertex)>0)
  {
    VDATA(pv) = GetMemoryForObject(MYMG(theGrid),ds,NOOBJ);
    if (VDATA(pv)==NULL) return(NULL);
  }
  else
    VDATA(pv) = NULL;

  /* initialize data */
  CTRL(pv) = 0;
  SETOBJT(pv,IVOBJ);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
  TOPNODE(pv) = NULL;
  SETMOVE(pv,DIM);
        #ifdef ModelP
  DDD_AttrSet(PARHDRV(pv),theGrid->level);
  DDD_PrioritySet(PARHDRV(pv),PrioMaster);
        #endif
  for (i=0; i<DIM; i++) LCVECT(pv)[i] = 0.0;

  /* insert in vertex list */
  GRID_LINK_VERTEX(theGrid,pv,PrioMaster);

  return(pv);
}

/****************************************************************************/
/*D
   CreateNode - Return pointer to a new node structure

   SYNOPSIS:
   static NODE *CreateNode (GRID *theGrid);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  after - node after which to insert new node

   DESCRIPTION:
   This function creates and initializes a new node structure
   and returns a pointer to it.

   RETURN VALUE:
   NODE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

static NODE *CreateNode (GRID *theGrid)
{
  NODE *pn;
  VECTOR *pv;
  INT size;

  size = sizeof(NODE);
  if (!TYPE_DEF_IN_GRID(theGrid,NODEVECTOR)) size -= sizeof(VECTOR *);
  if (NDATA_DEF_IN_GRID(theGrid)) size += sizeof(void *);
  if (NELIST_DEF_IN_GRID(theGrid)) size += sizeof(void *);

  pn = GetMemoryForObject(MYMG(theGrid),size,NDOBJ);
  if (pn==NULL) return(NULL);

  if (TYPE_DEF_IN_GRID(theGrid,NODEVECTOR)) {
    pv = CreateVector (theGrid,NODEVECTOR,(GEOM_OBJECT *)pn);
    if (pv == NULL)
    {
      DisposeNode (theGrid,pn);
      return (NULL);
    }
    NVECTOR(pn) = pv;
  }

  if (NDATA_DEF_IN_GRID(theGrid)) {
    NDATA(pn) = (void *) GetMemoryForObject(theGrid->mg,
                                            NDATA_DEF_IN_GRID(theGrid),-1);
    if (NDATA(pn) == NULL) {
      DisposeNode (theGrid,pn);
      return (NULL);
    }
  }

  /* initialize data */
  SETOBJT(pn,NDOBJ);
  SETCLASS(pn,4);
  SETLEVEL(pn,theGrid->level);
        #ifdef ModelP
  DDD_AttrSet(PARHDR(pn),theGrid->level);
  DDD_PrioritySet(PARHDR(pn),PrioMaster);
        #endif
  ID(pn) = (theGrid->mg->nodeIdCounter)++;
  START(pn) = NULL;
  SONNODE(pn) = NULL;
  if (NELIST_DEF_IN_GRID(theGrid)) NDATA(pn) = NULL;

  theGrid->status |= 1;          /* recalculate stiffness matrix */

  /* insert in vertex list */
  GRID_LINK_NODE(theGrid,pn,PrioMaster);

  return(pn);
}

/****************************************************************************/
/*D
   CreateSonNode - Return pointer to a new node structure on an edge

   SYNOPSIS:
   NODE *CreateSonNode (GRID *theGrid, NODE *FatherNode);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  FatherNode - node father

   DESCRIPTION:
   This function creates and initializes a new node structure
   at the midpoint of an element edge and returns a pointer to it.

   RETURN VALUE:
   NODE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

NODE *CreateSonNode (GRID *theGrid, NODE *FatherNode)
{
  NODE *pn;
  VERTEX *theVertex;

  pn = CreateNode(theGrid);
  if (pn == NULL)
    return(NULL);
  SONNODE(FatherNode) = pn;
  NFATHER(pn) = FatherNode;
  theVertex = MYVERTEX(FatherNode);
  MYVERTEX(pn) = theVertex;
  TOPNODE(theVertex) = pn;
  SETNTYPE(pn,CORNER_NODE);

  return(pn);
}

/****************************************************************************/
/*D
   CreateMidNode - Return pointer to a new node structure on an edge

   SYNOPSIS:
   NODE *CreateMidNode (GRID *theGrid, ELEMENT *theElement, INT edge);

   PARAMETERS:
   .  theGrid - grid where node should be inserted
   .  theElement - pointer to an element
   .  edge - id of an element edge

   DESCRIPTION:
   This function creates and initializes a new node structure
   at the midpoint of an element edge and returns a pointer to it.

   RETURN VALUE:
   NODE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

NODE *CreateMidNode (GRID *theGrid, ELEMENT *theElement, INT edge)
{
  NODE *theNode;
  EDGE *theEdge;
  VERTEX *theVertex,*v0,*v1;
  BNDP *bndp;
  DOUBLE *local,*x[MAX_CORNERS_OF_ELEM];
  DOUBLE_VECTOR bnd_global,global;
  DOUBLE diff;
  INT n,co0,co1;

  co0 = CORNER_OF_EDGE(theElement,edge,0);
  co1 = CORNER_OF_EDGE(theElement,edge,1);
  v0 = MYVERTEX(CORNER(theElement,co0));
  v1 = MYVERTEX(CORNER(theElement,co1));
  V_DIM_LINCOMB(0.5, CVECT(v0), 0.5, CVECT(v1), global);
  theVertex = NULL;

  /* allocate vertex */
  if ((OBJT(v0) == BVOBJ) && (OBJT(v1) == BVOBJ))
          #ifdef __TWODIM__
    if (OBJT(theElement) == BEOBJ)
      if (SIDE_ON_BND(theElement,edge))
      #endif
  {
    bndp = BNDP_CreateBndP(MGHEAP(MYMG(theGrid)),
                           V_BNDP(v0),V_BNDP(v1),0.5);
    if (bndp != NULL)
    {
      theVertex = CreateBoundaryVertex(theGrid);
      if (theVertex == NULL)
        return(NULL);
      if (BNDP_Global(bndp,bnd_global))
        return(NULL);
      V_BNDP(theVertex) = bndp;
      V_DIM_COPY(bnd_global,CVECT(theVertex));
      local = LCVECT(theVertex);
      V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
      if (diff > MAX_PAR_DIST)
      {
        SETMOVED(theVertex,1);
        CORNER_COORDINATES(theElement,n,x);
        UG_GlobalToLocal(n,(const DOUBLE **)x,global,local);
      }
      else
        V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theElement,co0),
                      0.5, LOCAL_COORD_OF_ELEM(theElement,co1),local);
    }
  }

  if (theVertex == NULL)
  {
    /* we need an inner vertex */
    theVertex = CreateInnerVertex(theGrid);
    if (theVertex==NULL) return(NULL);
    V_DIM_COPY(global,CVECT(theVertex));
    V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theElement,co0),
                  0.5, LOCAL_COORD_OF_ELEM(theElement,co1),
                  LCVECT(theVertex));
  }
  VFATHER(theVertex) = theElement;
  SETONEDGE(theVertex,edge);

  /* allocate node */
  theNode = CreateNode(theGrid);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }

  /* set MIDNODE pointer */
  theEdge = GetEdge(CORNER(theElement,co0),CORNER(theElement,co1));
  ASSERT(theEdge!=NULL);

  MIDNODE(theEdge) = theNode;
  MYVERTEX(theNode) = theVertex;
  NFATHER(theNode) = NULL;
  TOPNODE(theVertex) = theNode;
  SETNTYPE(theNode,MID_NODE);

  if (OBJT(theVertex) == BVOBJ)
    PRINTDEBUG(dom,1,(" MidPoint %d %f %f %f\n",ID(theNode),
                      bnd_global[0],
                      bnd_global[1],
                      bnd_global[2]));

  return(theNode);
}

/****************************************************************************/
/*D
   CreateSideNode - Return pointer to a new node structure on a side (3d)

   SYNOPSIS:
   NODE *CreateSideNode (GRID *theGrid, ELEMENT *theElement, INT side);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  theElement - pointer to an element
   .  side - id of an element side

   DESCRIPTION:
   This function creates and initializes a new node structure
   at the midpoint of an element side and returns a pointer to it.

   RETURN VALUE:
   NODE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

static INT SideOfNbElement(ELEMENT *theElement, INT side)
{
  ELEMENT *nb;
  NODE *nd[MAX_CORNERS_OF_SIDE];
  INT i,j,m,n,num;

  nb = NBELEM(theElement,side);
  if (nb == NULL) return(MAX_SIDES_OF_ELEM);

  n = CORNERS_OF_SIDE(theElement,side);
  for (i=0; i<n; i++)
    nd[i] = CORNER(theElement,CORNER_OF_SIDE(theElement,side,i));

  for (j=0; j<SIDES_OF_ELEM(nb); j++) {
    num = 0;
    for (i=0; i<n; i++)
      for (m=0; m<CORNERS_OF_SIDE(nb,j); m++)
        if (nd[i] == CORNER(nb,CORNER_OF_SIDE(nb,j,m))) num++;
    if (num == n) return(j);
  }

  return(MAX_SIDES_OF_ELEM);
}

#ifdef __THREEDIM__
NODE *CreateSideNode (GRID *theGrid, ELEMENT *theElement, INT side)
{
  DOUBLE_VECTOR bnd_global,global,local,bnd_local;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  VERTEX *theVertex;
  NODE *theNode;
  BNDP *bndp;
  BNDS *bnds;
  DOUBLE fac, diff;
  INT n,j,k;

  n = CORNERS_OF_SIDE(theElement,side);
  fac = 1.0 / n;
  V_DIM_CLEAR(local);
  V_DIM_CLEAR(global);
  for (j=0; j<n; j++)
  {
    k = CORNER_OF_SIDE(theElement,side,j);
    V_DIM_LINCOMB(1.0,local,1.0,
                  LOCAL_COORD_OF_ELEM(theElement,k),local);
    V_DIM_LINCOMB(1.0,global,1.0,
                  CVECT(MYVERTEX(CORNER(theElement,k))),global);
  }
  V_DIM_SCALE(fac,local);
  V_DIM_SCALE(fac,global);
  theVertex = NULL;

  /* check if boundary vertex */
  if (OBJT(theElement) == BEOBJ)
  {
    bnds = ELEM_BNDS(theElement,side);
    if (bnds != NULL)
    {
      if (n == 3)
        bnd_local[0] = bnd_local[1] = 0.33333333333333;
      else if (n == 4)
        bnd_local[0] = bnd_local[1] = 0.5;
      bndp = BNDS_CreateBndP(MGHEAP(MYMG(theGrid)),bnds,bnd_local);
      if (bndp != NULL)
      {
        theVertex = CreateBoundaryVertex(theGrid);
        if (theVertex == NULL)
          return(NULL);
        if (BNDP_Global(bndp,bnd_global))
          return(NULL);
        V_BNDP(theVertex) = bndp;
        V_DIM_COPY(bnd_global,CVECT(theVertex));
        V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
        if (diff > MAX_PAR_DIST)
        {
          SETMOVED(theVertex,1);
          CORNER_COORDINATES(theElement,n,x);
          UG_GlobalToLocal(n,(const DOUBLE **)x,global,local);
        }
        else
          V_DIM_COPY(local,LCVECT(theVertex));
      }
    }
  }

  if (theVertex == NULL)
  {
    theVertex = CreateInnerVertex(theGrid);
    if (theVertex == NULL) return(NULL);
    V_DIM_COPY(global,CVECT(theVertex));
  }
  VFATHER(theVertex) = theElement;
  SETONSIDE(theVertex,side);
  SETONNBSIDE(theVertex,SideOfNbElement(theElement,side));
  V_DIM_COPY(local,LCVECT(theVertex));

  /* create node */
  theNode = CreateNode(theGrid);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  MYVERTEX(theNode) = theVertex;
  TOPNODE(theVertex) = theNode;
  NFATHER(theNode) = NULL;
  SETNTYPE(theNode,SIDE_NODE);
  theGrid->status |= 1;

  return(theNode);
}

NODE *GetSideNode (ELEMENT *theElement, NODE *theNode0, NODE *theNode1, INT side)
{
  NODE *theNode=NULL;
  LINK *theLink0,*theLink1;
  INT l=0;

  ASSERT(theNode0!=NULL || theNode1!=NULL);

  /* search over one link list */
  if (theNode0!=NULL || theNode1!=NULL) {
    if (theNode0!=NULL)
      theLink0 = START(theNode0);
    if (theNode1!=NULL)
      theLink0 = START(theNode1);

    for ( ; theLink0!=NULL; theLink0=NEXT(theLink0)) {
      if (NTYPE(NBNODE(theLink0)) == SIDE_NODE) {
        theNode = NBNODE(theLink0);
        if (VFATHER(MYVERTEX(theNode)) == theElement) {
          if (ONSIDE(MYVERTEX(theNode)) == side)
            break;
        }
        if (NBELEM(theElement,side)==VFATHER(MYVERTEX(theNode))) {
          if (ONNBSIDE(MYVERTEX(theNode)) == side)
            break;
        }
        theNode = NULL;
      }
    }
  }

  /* search over both link lists */
  if (theNode0!=NULL && theNode1!=NULL)
    for (theLink0=START(theNode0); theLink0!=NULL; theLink0=NEXT(theLink0)) {
      for (theLink1=START(theNode1); theLink1!=NULL; theLink1=NEXT(theLink1))
        if (NBNODE(theLink0) == NBNODE(theLink1)) {
          if (NTYPE(NBNODE(theLink0)) == SIDE_NODE) {
            theNode = NBNODE(theLink0);
            if (VFATHER(MYVERTEX(theNode)) == theElement) {
              if (ONSIDE(MYVERTEX(theNode)) == side)
                break;
            }
            if (NBELEM(theElement,side)==VFATHER(MYVERTEX(theNode))) {
              if (ONNBSIDE(MYVERTEX(theNode)) == side)
                break;
            }
            theNode = NULL;
          }
        }
      if (theNode != NULL) break;
    }

  ASSERT(theNode==NULL || NTYPE(theNode)==SIDE_NODE &&
         (ONSIDE(MYVERTEX(theNode)) == side ||
          ONNBSIDE(MYVERTEX(theNode)) == side));

  return(theNode);
}
#endif /* __THREEDIM__ */

/****************************************************************************/
/*																			*/
/* Function:  CreateCenterNode	                                                                                        */
/*																			*/
/* Purpose:   allocate a new node on an side of an element. Includes vertex */
/*			  best fit boundary coordinates and local coordinates.			*/
/*																			*/
/* return:	  NODE* : pointer to new node									*/
/*			  NULL	: could not allocate									*/
/*																			*/
/****************************************************************************/

NODE *CreateCenterNode (GRID *theGrid, ELEMENT *theElement)
{
  DOUBLE *global,*local;
  INT n,m,j,moved;
  VERTEX *theVertex,*VertexOnEdge[MAX_EDGES_OF_ELEM];
  NODE *theNode;
  EDGE *theEdge;
  DOUBLE fac;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];

  /* check if moved side nodes exist */
  moved = 0;
  if (OBJT(theElement) == BEOBJ)
    for (j=0; j<EDGES_OF_ELEM(theElement); j++)
    {
      theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),
                      CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
      ASSERT(theEdge != NULL);
      theNode = MIDNODE(theEdge);
      if (theNode == NULL)
        VertexOnEdge[j] = NULL;
      else
      {
        VertexOnEdge[j] = MYVERTEX(theNode);
        moved += MOVED(VertexOnEdge[j]);
      }
    }

  theVertex = CreateInnerVertex(theGrid);
  if (theVertex==NULL)
    return(NULL);
  theNode = CreateNode(theGrid);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  CORNER_COORDINATES(theElement,n,x);
  global = CVECT(theVertex);
  local = LCVECT(theVertex);
  if (moved == 0)
  {
    V_DIM_CLEAR(local);
    fac = 1.0 / n;
    for (j=0; j<n; j++)
      V_DIM_LINCOMB(1.0,local,
                    fac,LOCAL_COORD_OF_ELEM(theElement,j),local);
    LOCAL_TO_GLOBAL(n,x,local,global);
  }
  else
  {
    V_DIM_CLEAR(global);
    m = EDGES_OF_ELEM(theElement);
    fac = 1.0 / m;
    for (j=0; j<m; j++)
      if (VertexOnEdge[j] == NULL)
      {
        V_DIM_LINCOMB(1.0,global,0.5*fac,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)))),global);
        V_DIM_LINCOMB(1.0,global,0.5*fac,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)))),global);
      }
      else
        V_DIM_LINCOMB(1.0,global,fac,CVECT(VertexOnEdge[j]),global);
    UG_GlobalToLocal(n,(const DOUBLE **)x,global,local);
  }
  VFATHER(theVertex) = theElement;
  NFATHER(theNode) = NULL;
  MYVERTEX(theNode) = theVertex;
  TOPNODE(theVertex) = theNode;
  SETNTYPE(theNode,CENTER_NODE);
  theGrid->status |= 1;

  return(theNode);
}

/****************************************************************************/
/*D
   GetEdge - Return pointer to edge if it exists

   SYNOPSIS:
   EDGE *GetEdge (NODE *from, NODE *to);

   PARAMETERS:
   .  from - starting node of edge
   .  to - end node of edge

   DESCRIPTION:
   This function returns the pointer to the specified edge if it exists.

   RETURN VALUE:
   EDGE *
   .n   pointer to specified object
   .n   NULL if not found
   D*/
/****************************************************************************/

EDGE *GetEdge (NODE *from, NODE *to)
{
  LINK *pl;

  /* run through neighbor list */
  for (pl=START(from); pl!=NULL; pl = NEXT(pl))
    if (NBNODE(pl)==to)
      return(MYEDGE(pl));

  /* return not found */
  return(NULL);
}

/****************************************************************************/
/*D
   CreateEdge - Return pointer to a new edge structure

   SYNOPSIS:
   EDGE *CreateEdge (GRID *theGrid, NODE *from, NODE *to, INT with_vector);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  from - starting node of new edge
   .  to - end node of new edge
   .  with_vector - also create vector for edge (TRUE/FALSE)

   DESCRIPTION:
   This function returns a pointer to a new edge structure.

   RETURN VALUE:
   EDGE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

#ifndef ModelP
static
#endif
EDGE *CreateEdge (GRID *theGrid, NODE *from, NODE *to, INT with_vector)
{
  EDGE *pe;
  LINK *link0,*link1;
  VECTOR *pv;

  /* check if edge exists already */
  if( (pe = GetEdge(from, to)) != NULL )
  {
    if (NO_OF_ELEM(pe)<NO_OF_ELEM_MAX-1)
      INC_NO_OF_ELEM(pe);
    else
      return (NULL);

    return(pe);
  }

  if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR))
    pe = GetMemoryForObject(theGrid->mg,sizeof(EDGE),EDOBJ);
  else
    pe = GetMemoryForObject(theGrid->mg,sizeof(EDGE)-sizeof(VECTOR*),EDOBJ);
  if (pe==NULL) return(NULL);

  /* initialize data */
  link0 = LINK0(pe);
  link1 = LINK1(pe);
  SETOBJT(pe,EDOBJ);
  SETLOFFSET(link0,0);
  SETLOFFSET(link1,1);
  SETLEVEL(pe,theGrid->level);
        #if (defined ModelP) && (defined __THREEDIM__)
  DDD_AttrSet(PARHDR(pe),theGrid->level);
  DDD_PrioritySet(PARHDR(pe),PrioMaster);
        #endif
  NBNODE(link0) = to;
  NBNODE(link1) = from;
  SET_NO_OF_ELEM(pe,1);
  SETTAG(pe,0);
  SETEDGENEW(pe,1);

  /* create vector if */
  if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR) && with_vector)
  {
    pv = CreateVector (theGrid,EDGEVECTOR,(GEOM_OBJECT *)pe);
    if (pv == NULL)
    {
      DisposeEdge (theGrid,pe);
      return (NULL);
    }
    EDVECTOR(pe) = pv;
  }

  /* put in neighbor lists */
  NEXT(link0) = START(from);
  START(from) = link0;
  NEXT(link1) = START(to);
  START(to) = link1;

  /* counters */
  theGrid->nEdge++;

  /* return ok */
  return(pe);
}

/****************************************************************************/
/*D
   GetLink - Return pointer to link if it exists

   SYNOPSIS:
   LINK *GetLink (NODE *from, NODE *to);

   PARAMETERS:
   .  from - starting node of link
   .  to - end node of link

   DESCRIPTION:
   This function returns the pointer to the specified link if it exists.

   RETURN VALUE:
   INT
   .n   pointer to specified link
   .n   NULL if not found.
   D*/
/****************************************************************************/

LINK *GetLink (NODE *from, NODE *to)
{
  LINK *pl;

  /* run through neighbor list */
  for (pl=START(from); pl!=NULL; pl = NEXT(pl))
    if (NBNODE(pl)==to)
      return(pl);

  /* return not found */
  return(NULL);
}

/****************************************************************************/
/*D
   CreateElement - Return a pointer to  a new element structure

   SYNOPSIS:
   ELEMENT *CreateElement (GRID *theGrid, INT tag, INT objtype,
   NODE **nodes, ELEMENT *Father);

   PARAMETERS:
   .  theGrid - grid structure to extend
   .  tag - the element type
   .  objtype - inner element (IEOBJ) or boundary element (BEOBJ)
   .  nodes - list of corner nodes in reference numbering
   .  Father - pointer to father element (NULL on base level)

   DESCRIPTION:
   This function creates and initializes a new element and returns a pointer to it.

   RETURN VALUE:
   ELEMENT *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

ELEMENT *CreateElement (GRID *theGrid, INT tag, INT objtype,
                        NODE **nodes, ELEMENT *Father)
{
  ELEMENT *pe;
  INT i;
  VECTOR *pv;
  void *q;

  if (objtype == IEOBJ)
    pe = GetMemoryForObject(MYMG(theGrid),INNER_SIZE_TAG(tag),
                            MAPPED_INNER_OBJT_TAG(tag));
  else if (objtype == BEOBJ)
    pe = GetMemoryForObject(MYMG(theGrid),BND_SIZE_TAG(tag),
                            MAPPED_BND_OBJT_TAG(tag));

  if (pe==NULL) return(NULL);

  /* initialize data */
  SETOBJT(pe,objtype);
  SETTAG(pe,tag);
  SETLEVEL(pe,theGrid->level);
        #ifdef ModelP
  DDD_AttrSet(PARHDRE(pe),theGrid->level);
  DDD_PrioritySet(PARHDRE(pe),PrioMaster);
        #endif
  SETEBUILDCON(pe,1);
  ID(pe) = (theGrid->mg->elemIdCounter)++;

  /* set corner nodes */
  for (i=0; i<CORNERS_OF_ELEM(pe); i++)
    SET_CORNER(pe,i,nodes[i]);

  /* create edges */
  for (i=0; i<EDGES_OF_ELEM(pe); i++)
    if (CreateEdge(theGrid,
                   nodes[CORNER_OF_EDGE(pe,i,0)],
                   nodes[CORNER_OF_EDGE(pe,i,1)],
                   TRUE) == NULL)
    {
      DisposeElement(theGrid,pe,TRUE);
      return(NULL);
    }

  /* create element vector if */
  if (TYPE_DEF_IN_GRID(theGrid,ELEMVECTOR))
  {
    pv = CreateVector (theGrid,ELEMVECTOR,(GEOM_OBJECT *)pe);
    if (pv == NULL)
    {
      DisposeElement (theGrid,pe,TRUE);
      return (NULL);
    }
    SET_EVECTOR(pe,pv);
  }

  if (EDATA_DEF_IN_GRID(theGrid)) {
    q = (void *) GetMemoryForObject(theGrid->mg,EDATA_DEF_IN_GRID(theGrid),-1);
    if (q == NULL) {
      DisposeElement (theGrid,pe,TRUE);
      return (NULL);
    }
    SET_EDATA(pe,q);
  }

  /* create side vectors if */
        #ifdef __THREEDIM__
  if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
    for (i=0; i<SIDES_OF_ELEM(pe); i++)
    {
      pv = CreateSideVector (theGrid,i,(GEOM_OBJECT *)pe);
      if (pv == NULL)
      {
        DisposeElement (theGrid,pe,TRUE);
        return (NULL);
      }
      SET_SVECTOR(pe,i,pv);
    }
        #endif

  /* insert in element list */
  GRID_LINK_ELEMENT(theGrid,pe,PrioMaster);

  SET_EFATHER(pe,Father);
        #ifndef ModelP
  if (Father != NULL)
    SETSUBDOMAIN(pe,SUBDOMAIN(Father));
        #else
  SETSUBDOMAIN(pe,me+1);
        #endif

  /* return ok */
  return(pe);
}

/****************************************************************************/
/*D
   CreateSonElementSide - creates the element sides of son elements

   SYNOPSIS:
   INT CreateSonElementSide (GRID *theGrid, ELEMENT *theElement, INT side,
   ELEMENT *theSon, INT son_side);

   PARAMETERS:
   .  theGrid - grid for which to create
   .  theElement - pointer to a boundary element
   .  side - side id of a side of the element
   .  theSon - pointer to a son element
   .  son_side - side id of a side of the son

   DESCRIPTION:
   This function creates and initializes an element side of a son element.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT CreateSonElementSide (GRID *theGrid, ELEMENT *theElement, INT side,
                          ELEMENT *theSon, INT son_side)
{
  INT n,i;
  BNDS *bnds;
  BNDP *bndp[MAX_CORNERS_OF_ELEM];

  ASSERT (OBJT(theElement) == BEOBJ);

  if (ELEM_BNDS(theElement,side) == NULL)
    return(GM_OK);

  n = CORNERS_OF_SIDE(theSon,son_side);
  for (i=0; i<n; i++)
    bndp[i] = V_BNDP(MYVERTEX(CORNER(theSon,
                                     CORNER_OF_SIDE(theSon,son_side,i))));
  bnds = BNDP_CreateBndS(MGHEAP(MYMG(theGrid)),bndp,n);
  if (bnds == NULL)
    RETURN(GM_ERROR);
  SET_BNDS(theSon,son_side,bnds);

  return(GM_OK);
}

/****************************************************************************/
/*D
   CreateNewLevel - Return pointer to new grid structure

   SYNOPSIS:
   GRID *CreateNewLevel (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid structure

   DESCRIPTION:
   This function creates and initialized a new grid structure for top level + 1
   and returns a pointer to it.

   RETURN VALUE:
   GRID *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

GRID *CreateNewLevel (MULTIGRID *theMG)
{
  GRID *theGrid;
  INT l;

  if (theMG->topLevel+1>=MAXLEVEL) return(NULL);

  l = theMG->topLevel+1;

  /* allocate grid object */
  theGrid = GetMemoryForObject(theMG,sizeof(GRID),GROBJ);
  if (theGrid==NULL) return(NULL);

  /* fill in data */
  CTRL(theGrid) = 0;
  SETOBJT(theGrid,GROBJ);
  theGrid->level = l;
  theGrid->nEdge = 0;
  theGrid->nCon = 0;
  /* other counters are init in INIT fcts below */

#ifdef __INTERPOLATION_MATRIX__
  theGrid->nIMat = 0;
#endif

  theGrid->status       = 0;
  GRID_INIT_ELEMENT_LIST(theGrid);
  GRID_INIT_NODE_LIST(theGrid);
  GRID_INIT_VERTEX_LIST(theGrid);
  GRID_INIT_VECTOR_LIST(theGrid);
  GFIRSTBV(theGrid) = NULL;
  GLASTBV(theGrid) = NULL;
  if (l>0)
  {
    theGrid->coarser = theMG->grids[l-1];
    theMG->grids[l-1]->finer = theGrid;
  }
  else
    theGrid->coarser = NULL;
  theGrid->finer = NULL;
  theGrid->mg = theMG;
  theMG->grids[l] = theGrid;
  theMG->topLevel = l;
  theMG->currentLevel = l;

  return(theGrid);
}

/****************************************************************************/
/*D
   MakeMGItem - Create a multigrid environment item

   SYNOPSIS:
   MULTIGRID *MakeMGItem (const char *name);

   PARAMETERS:
   .  name - name of the multigrid

   DESCRIPTION:
   This function creates a multigrid environment directory.

   RETURN VALUE:
   MULTIGRID *
   .n   pointer to new MULTIGRID
   .n   NULL if error occured
   D*/
/****************************************************************************/

MULTIGRID *MakeMGItem (const char *name)
{
  MULTIGRID *theMG;

  if (ChangeEnvDir("/Multigrids") == NULL) return (NULL);
  if (strlen(name)>=NAMESIZE || strlen(name)<=1) return (NULL);
  theMG = (MULTIGRID *) MakeEnvItem(name,theMGDirID,sizeof(MULTIGRID));
  if (theMG == NULL) return(NULL);

  return (theMG);
}

/****************************************************************************/
/*D
   GetMultigrid - Find the multigrid environment item with name

   SYNOPSIS:
   MULTIGRID *GetMultigrid (const char *name);

   PARAMETERS:
   .  name - name of the multigrid to find

   DESCRIPTION:
   This function find the multigrid environment item with `name` and
   returns a pointer to the multigrid structure.

   RETURN VALUE:
   MULTIGRID *
   .n   pointer to MULTIGRID
   .n   NULL if not found.
   D*/
/****************************************************************************/

MULTIGRID *GetMultigrid (const char *name)
{
  return ((MULTIGRID *) SearchEnv(name,"/Multigrids",
                                  theMGDirID,theMGRootDirID));
}

/****************************************************************************/
/*D
   GetFirstMultigrid - Return a pointer to the first multigrid

   SYNOPSIS:
   MULTIGRID *GetFirstMultigrid ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function returns a pointer to the first multigrid in the /Multigrids
   directory.

   RETURN VALUE:
   MULTIGRID *
   .n   pointer to MULTIGRID
   .n   NULL if not found.
   D*/
/****************************************************************************/

MULTIGRID *GetFirstMultigrid ()
{
  ENVDIR *theMGRootDir;
  MULTIGRID *theMG;

  theMGRootDir = ChangeEnvDir("/Multigrids");

  assert (theMGRootDir!=NULL);

  theMG = (MULTIGRID *) ENVDIR_DOWN(theMGRootDir);

  if (theMG != NULL)
    if (InitElementTypes(theMG) != GM_OK) {
      PrintErrorMessage('E',"GetFirstMultigrid",
                        "error in InitElementTypes");
      return(NULL);
    }

  return (theMG);
}

/****************************************************************************/
/*D
   GetNextMultigrid - Return a pointer to the next multigrid

   SYNOPSIS:
   MULTIGRID *GetNextMultigrid (const MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid structure

   DESCRIPTION:
   This function returns a pointer to the next multigrid in the /Multigrids
   directory.

   RETURN VALUE:
   MULTIGRID *
   .n   pointer to MULTIGRID
   .n   NULL if not found.
   D*/
/****************************************************************************/

MULTIGRID *GetNextMultigrid (const MULTIGRID *theMG)
{
  MULTIGRID *MG;

  MG = (MULTIGRID *) NEXT_ENVITEM(theMG);

  if (MG != NULL)
    if (InitElementTypes(MG)!=GM_OK) {
      PrintErrorMessage('E',"GetNextMultigrid",
                        "error in InitElementTypes");
      return(NULL);
    }

  return (MG);
}

/****************************************************************************/
/*D
   CreateMultiGrid - Return a pointer to new multigrid structure

   SYNOPSIS:
   MULTIGRID *CreateMultiGrid (char *MultigridName, char *domain, char *problem,
   char *format, unsigned long heapSize);

   PARAMETERS:
   .  MultigridName - name of multigrid
   .  domain - name of domain description from environment
   .  problem - name of problem description from environment
   .  format - name of format description from environment
   .  heapSize - size of heap to allocate for that multigrid in bytes

   DESCRIPTION:
   This function creates and initializes a new multigrid structure including
   allocation of heap, combining the domain and the boundary conditions
   and creation of the fixed corners of the domain.

   RETURN VALUE:
   MULTIGRID *
   .n   pointer to new object
   .n   NULL if an error occured.
   D*/
/****************************************************************************/

MULTIGRID *CreateMultiGrid (char *MultigridName, char *BndValProblem,
                            char *format, MEM heapSize)
{
  HEAP *theHeap,*theUserHeap;
  MULTIGRID *theMG;
  GRID *theGrid;
  INT i,ds;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  MESH mesh;
  FORMAT *theFormat;

  theFormat = GetFormat(format);
  if (theFormat==NULL)
  {
    PrintErrorMessage('E',"CreateMultiGrid","format not found");
    return(NULL);
  }

  /* allocate multigrid envitem */
  theMG = MakeMGItem(MultigridName);
  if (theMG==NULL) return(NULL);
  theMG->theFormat = theFormat;
  if (InitElementTypes(theMG)!=GM_OK)
  {
    PrintErrorMessage('E',"CreateMultiGrid","error in InitElementTypes");
    return(NULL);
  }

  /* allocate the heap */
  theHeap = NewHeap(SIMPLE_HEAP, heapSize, malloc(heapSize));
  if (theHeap==NULL)
  {
    PRINTDEBUG(gm,0,("CreateMultiGrid: cannot allocate %ld bytes\n",
                     heapSize));
    DisposeMultiGrid(theMG);
    return(NULL);
  }

  MarkTmpMem(theHeap);
  theBVP = BVP_Init(BndValProblem,theHeap,&mesh);
  if (theBVP==NULL)
  {
    PrintErrorMessage('E',"CreateMultiGrid","BVP not found");
    return(NULL);
  }
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc))
  {
    PrintErrorMessage('E',"CreateMultiGrid","BVP not evaluated");
    return(NULL);
  }

  /* 1: general user data space */
  if (!theGenMGUDM->locked)
    CalcAndFixTotalSize(theGenMGUDM);
  ds = theGenMGUDM->TotalSize;
  if (ds!=0)
  {
    GEN_MGUD(theMG) = GetMem(theHeap,ds,FROM_BOTTOM);
    if (GEN_MGUD(theMG)==NULL)
    {
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    /* clearing this heap provides the possibility of checking the
       initialization */
    memset(GEN_MGUD(theMG),0,ds);
  }
  else
    GEN_MGUD(theMG) = NULL;

  /* 2: user heap */
  ds = theFormat->sMultiGrid;
  if (ds!=0)
  {
    theUserHeap = NewHeap(SIMPLE_HEAP, ds, GetMem(theHeap,ds,FROM_BOTTOM));
    if (theUserHeap==NULL)
    {
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    MG_USER_HEAP(theMG) = theUserHeap;
  }
  else
    MG_USER_HEAP(theMG) = NULL;

  /* fill multigrid structure */
  theMG->status = 0;
  theMG->vertIdCounter = 0;
  theMG->nodeIdCounter = 0;
  theMG->elemIdCounter = 0;
  theMG->topLevel = -1;
  MG_BVP(theMG) = theBVP;
        #if defined(CAD) && defined(__THREEDIM__)
  MG_NPROPERTY(theMG) = 1;
        #else
  MG_NPROPERTY(theMG) = BVPD_NSUBDOM(theBVPDesc);
        #endif
  RESETMGSTATUS(theMG);
  MG_GENPURP(theMG) = NULL;

  theMG->theHeap = theHeap;
  SELECTIONSIZE(theMG) = 0;
  for (i=0; i<MAXLEVEL; i++)
    theMG->grids[i] = NULL;

  /* allocate level 0 grid */
  theGrid = CreateNewLevel(theMG);
  if (theGrid==NULL)
  {
    DisposeMultiGrid(theMG);
    return(NULL);
  }

  /* allocate predefined mesh, e. g. corner vertices pointers */
        #ifndef CAD
  if (InsertMesh(theMG,&mesh))
  {
    DisposeMultiGrid(theMG);
    return(NULL);
  }
        #endif

  ReleaseTmpMem(theHeap);

  /* return ok */
  return(theMG);
}

/****************************************************************************/
/*D
   DisposeEdge - Remove edge from the data structure

   SYNOPSIS:
   static INT DisposeEdge (GRID *theGrid, EDGE *theEdge);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  theEdge - edge to remove

   DESCRIPTION:
   This function remove an edge from the data structure including its
   vector (if one) and inserts them into the free list.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if an error occured.
   D*/
/****************************************************************************/

static INT DisposeEdge (GRID *theGrid, EDGE *theEdge)
{
  LINK *link0,*link1,*pl;
  NODE *from,*to;
  INT found;

  HEAPFAULT(theEdge);

  /* reconstruct data */
  link0 = LINK0(theEdge);
  link1 = LINK1(theEdge);
  from  = NBNODE(link1);
  to        = NBNODE(link0);
  found = 0;

  /* delete link0 in from vertex */
  if (START(from)==link0)
  {
    START(from) = NEXT(link0);
    found++;
  }
  else
  {
    for (pl=START(from); pl!=NULL; pl = NEXT(pl))
    {
      if (NEXT(pl)==link0)
      {
        NEXT(pl) = NEXT(link0);
        found++;
        break;
      }
    }
  }

  /* delete link1 in to vertex */
  if (START(to)==link1)
  {
    START(to) = NEXT(link1);
    found++;
  }
  else
  {
    for (pl=START(to); pl!=NULL; pl = NEXT(pl))
    {
      if (NEXT(pl)==link1)
      {
        NEXT(pl) = NEXT(link1);
        found++;
        break;
      }
    }
  }

  /* dispose vector and its matrices from edge-vector if */
  if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR))
  {
    if (DisposeVector (theGrid,EDVECTOR(theEdge)))
      RETURN(1);
    PutFreeObject(theGrid->mg,theEdge,sizeof(EDGE),EDOBJ);
  }
  else
    PutFreeObject(theGrid->mg,theEdge,sizeof(EDGE)-sizeof(VECTOR*),EDOBJ);

  /* check error condition */
  if (found!=2) RETURN(1);

  /* return ok */
  theGrid->nEdge--;
  return(0);
}

/****************************************************************************/
/*D
   DisposeNode - Remove node including its edges from the data structure

   SYNOPSIS:
   static INT DisposeNode (GRID *theGrid, NODE *theNode);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  theNode - node to remove

   DESCRIPTION:
   This function removes node including its edges and vector (if one)
   from the data structure and inserts all objects into the free list.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 when error occured.
   D*/
/****************************************************************************/

static INT DisposeNode (GRID *theGrid, NODE *theNode)
{
  VERTEX *theVertex;
  NODE *father;
  INT size;

  HEAPFAULT(theNode);

  /* call DisposeElement first! */
  assert(START(theNode) == NULL);
  assert(SONNODE(theNode) == NULL);

  /* remove node from node list */
  GRID_UNLINK_NODE(theGrid,theNode);

  theVertex = MYVERTEX(theNode);
  father = NFATHER(theNode);
  if (father != NULL)
  {
    SONNODE(father) = NULL;
    if (theVertex != NULL)
      TOPNODE(theVertex) = father;
  }
  else if (theVertex != NULL)
    DisposeVertex(theGrid,theVertex);

  /* dispose vector and its matrices from node-vector */
  size = sizeof(NODE);
  if (NDATA_DEF_IN_GRID(theGrid)) {
    size += sizeof(void *);
    PutFreeObject(theGrid->mg,NDATA(theNode),NDATA_DEF_IN_GRID(theGrid),-1);
  }
  if (NELIST_DEF_IN_GRID(theGrid)) {
    DisposeElementList(theGrid,theNode);
    size += sizeof(void *);
  }
  if (TYPE_DEF_IN_GRID(theGrid,NODEVECTOR))
  {
    if (DisposeVector (theGrid,NVECTOR(theNode)))
      RETURN(1);
  }
  else
    size -= sizeof(VECTOR *);
  PutFreeObject(theGrid->mg,theNode,size,NDOBJ);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*D
   DisposeVertex - Remove vertex from the data structure

   SYNOPSIS:
   static INT DisposeVertex (GRID *theGrid, VERTEX *theVertex);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  theVertex - vertex to remove

   DESCRIPTION:
   This function removes a vertex from the data structure
   and puts it into the free list.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 no valid object number
   D*/
/****************************************************************************/

static INT DisposeVertex (GRID *theGrid, VERTEX *theVertex)
{
  HEAPFAULT(theVertex);

  /* remove vertex from vertex list */
  GRID_UNLINK_VERTEX(theGrid,theVertex);

  if( OBJT(theVertex) == BVOBJ )
  {
    BNDP_Dispose(MGHEAP(MYMG(theGrid)),V_BNDP(theVertex));
    PutFreeObject(theGrid->mg,theVertex,sizeof(struct bvertex),BVOBJ);
  }
  else
    PutFreeObject(theGrid->mg,theVertex,sizeof(struct ivertex),IVOBJ);

  return(0);
}

/****************************************************************************/
/*D
   DisposeElement - Remove element from the data structure

   SYNOPSIS:
   INT DisposeElement (GRID *theGrid, ELEMENT *theElement, INT dispose_connections);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  theElement - element to remove
   .  dispose_connections - also dispose connections (TRUE/FALSE)

   DESCRIPTION:
   This function removes an element from the data structure and inserts it
   into the free list. This includes all elementsides, sidevectors and the
   elementvector if they exist.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 no valid object number.
   D*/
/****************************************************************************/

INT DisposeElement (GRID *theGrid, ELEMENT *theElement, INT dispose_connections)
{
  INT i,j,edge,tag;
  NODE *theNode;
  VERTEX *theVertex;
  EDGE *theEdge;
  BNDS *bnds;
  ELEMENT *theFather;
        #ifdef __THREEDIM__
  VECTOR *theVector;
  ELEMENT *theNeighbor;
        #endif

  HEAPFAULT(theElement);

  GRID_UNLINK_ELEMENT(theGrid,theElement);

  /* remove element sides if it's a boundary element */
  if (OBJT(theElement)==BEOBJ)
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      bnds = ELEM_BNDS(theElement,i);
      if (bnds != NULL)
      {
        HEAPFAULT(bnds);
        BNDS_Dispose(MGHEAP(MYMG(theGrid)),bnds);
      }
    }

  for (j=0; j<EDGES_OF_ELEM(theElement); j++)
  {
    theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),
                    CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
    ASSERT(theEdge!=NULL);

    if (NO_OF_ELEM(theEdge)<1)
      RETURN(GM_ERROR);
    if (NO_OF_ELEM(theEdge)==1)
      DisposeEdge(theGrid,theEdge);
    else
      DEC_NO_OF_ELEM(theEdge);
  }

  if (NELIST_DEF_IN_GRID(theGrid))
    for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      DisposeElementFromElementList(theGrid,
                                    CORNER(theElement,j),theElement);

  for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
  {
    theNode = CORNER(theElement,j);
    if (START(theNode) == NULL)
    {
      if (NTYPE(theNode)==MID_NODE)
      {
        theVertex = MYVERTEX(theNode);
        theFather = VFATHER(theVertex);
        edge = ONEDGE(theVertex);
        theEdge = GetEdge(CORNER(theFather,
                                 CORNER_OF_EDGE(theFather,edge,0)),
                          CORNER(theFather,
                                 CORNER_OF_EDGE(theFather,edge,1)));
        ASSERT(theEdge!=NULL);
        MIDNODE(theEdge) = NULL;
      }
      DisposeNode(theGrid,theNode);
    }
  }

  /* dispose matrices from element-vector */
  if (dispose_connections)
    if (DisposeConnectionFromElement(theGrid,theElement))
      RETURN(1);

  /* reset neighbor pointers referencing element and dispose vectors in sides if */
        #ifdef __TWODIM__
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    if (NBELEM(theElement,i)!=NULL) {
      for (j=0; j<SIDES_OF_ELEM(NBELEM(theElement,i)); j++)
        if (NBELEM(NBELEM(theElement,i),j)==theElement) {
          SET_NBELEM(NBELEM(theElement,i),j,NULL);
          break;
        }
      assert(j<SIDES_OF_ELEM(NBELEM(theElement,i)));
    }
  }
        #endif

        #ifdef __THREEDIM__
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR)) {
      theVector = SVECTOR(theElement,i);
      assert (VCOUNT(theVector) != 0);
      assert (VCOUNT(theVector) != 3);
      if (VCOUNT(theVector) == 1)
      {
        if (DisposeVector (theGrid,theVector))
          RETURN (1);
      }
      else
      {
        if (!FindNeighborElement (theElement,i,&theNeighbor,&j))
          RETURN (1);
        VOBJECT(theVector) = (void*)theNeighbor;
        SETVECTORSIDE(theVector,j);
        SETVCOUNT(SVECTOR(theElement,i),1);
      }
    }
    if (NBELEM(theElement,i)!=NULL) {
      for (j=0; j<SIDES_OF_ELEM(NBELEM(theElement,i)); j++)
        if (NBELEM(NBELEM(theElement,i),j)==theElement) {
          SET_NBELEM(NBELEM(theElement,i),j,NULL);
          break;
        }
      assert(j<SIDES_OF_ELEM(NBELEM(theElement,i)));
    }
  }
    #endif

  /* dispose vector in center of element */
  if (TYPE_DEF_IN_GRID(theGrid,ELEMVECTOR))
    if (DisposeVector (theGrid,EVECTOR(theElement)))
      RETURN(1);

  if (EDATA_DEF_IN_GRID(theGrid))
    PutFreeObject(theGrid->mg,EDATA(theElement),
                  EDATA_DEF_IN_GRID(theGrid),-1);

  /* dispose element */
  /* give it a new tag ! (I know this is somewhat ugly) */
  tag = TAG(theElement);
  if (OBJT(theElement)==BEOBJ)
  {
    SETOBJT(theElement,MAPPED_BND_OBJT_TAG(tag));
    PutFreeObject(theGrid->mg,theElement,
                  BND_SIZE_TAG(tag),MAPPED_BND_OBJT_TAG(tag));
  }
  else
  {
    SETOBJT(theElement,MAPPED_INNER_OBJT_TAG(tag));
    PutFreeObject(theGrid->mg,theElement,INNER_SIZE_TAG(tag),
                  MAPPED_INNER_OBJT_TAG(tag));
  }

  return(0);
}

/****************************************************************************/
/*D
   DisposeTopLevel - Remove top level grid from multigrid  structure

   SYNOPSIS:
   INT DisposeTopLevel (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid to remove from

   DESCRIPTION:
   This function removes the top level grid from multigrid structure.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 no valid object number
   .n   2 grid structure not empty or level 0
   D*/
/****************************************************************************/

INT DisposeTopLevel (MULTIGRID *theMG)
{
  int l;
  GRID *theGrid;

  /* level 0 can not be deleted */
  l = theMG->topLevel;
  if (l<=0) return(2);
  theGrid = theMG->grids[l];

  /* is level empty */
  if (FIRSTELEMENT(theGrid)!=NULL)
    return(2);
  if (FIRSTVERTEX(theGrid)!=NULL)
    return(2);
  if (FIRSTNODE(theGrid)!=NULL)
    return(2);

  /* remove from grids array */
  theMG->grids[l] = NULL;
  theMG->grids[l-1]->finer = NULL;
  (theMG->topLevel)--;
  if (theMG->currentLevel>theMG->topLevel)
    theMG->currentLevel = theMG->topLevel;

  PutFreeObject(theMG,theGrid,sizeof(GRID),GROBJ);

  return(0);
}

/****************************************************************************/
/*D
   DisposeGrid - dispose top level grid

   SYNOPSIS:
   INT DisposeGrid (GRID *theGrid);

   PARAMETERS:
   .  theMG - multigrid to remove from

   DESCRIPTION:
   This function removes the top level grid from multigrid structure.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 no valid object number
   .n   2 grid structure not empty or level 0
   D*/
/****************************************************************************/

INT DisposeGrid (GRID *theGrid)
{
  MULTIGRID *theMG;

  if (theGrid == NULL)
    return(0);

  if (theGrid->finer != NULL)
    return(1);

  /* clear level */
  while (FIRSTELEMENT(theGrid)!=NULL)
    if (DisposeElement(theGrid,FIRSTELEMENT(theGrid),1))
      return(2);

  while (FIRSTNODE(theGrid)!=NULL)
    if (DisposeNode(theGrid,FIRSTNODE(theGrid)))
      return(2);

  while (FIRSTVERTEX(theGrid)!=NULL)
    if (DisposeVertex(theGrid,FIRSTVERTEX(theGrid)))
      return(4);

  theMG = MYMG(theGrid);

  /* level 0 can not be deleted */
  if (GLEVEL(theGrid) > 0)
    return(DisposeTopLevel(theMG));

  /* remove from grids array */
  theMG->grids[0] = NULL;
  theMG->currentLevel = theMG->topLevel = -1;
  theMG->nodeIdCounter = 0;
  theMG->vertIdCounter = 0;
  theMG->elemIdCounter = 0;

  PutFreeObject(theMG,theGrid,sizeof(GRID),GROBJ);

  return(0);
}

/****************************************************************************/
/*D
   DisposeMultiGrid - Release memory for the whole multigrid  structure

   SYNOPSIS:
   INT DisposeMultiGrid (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid to remove

   DESCRIPTION:
   This function releases the memory for the whole multigrid  structure.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT DisposeMultiGrid (MULTIGRID *theMG)
{
  if (MGHEAP(theMG)!=NULL)
    free(MGHEAP(theMG));

  /* dispose BVP */
  if (BVP_Dispose(MG_BVP(theMG))) return (GM_ERROR);

  /* first unlock the mg */
  ((ENVITEM*) theMG)->v.locked = FALSE;

  /* delete mg */
  if (ChangeEnvDir("/Multigrids")==NULL) RETURN (GM_ERROR);
  if (RemoveEnvDir ((ENVITEM *)theMG)) RETURN (GM_ERROR);

  return(GM_OK);
}


/****************************************************************************/
/*D
   RenumberMultiGrid - Recalculate ids in the current order

   SYNOPSIS:
   INT RenumberMultiGrid (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - structure to renumber

   DESCRIPTION:
   This function recalculates ids in the current order.

   RETURN VALUE:
   INT
   .n   0 if ok
   D*/
/****************************************************************************/

INT RenumberMultiGrid (MULTIGRID *theMG)
{
  GRID *theGrid;
  VERTEX *theVertex;
  NODE *theNode;
  ELEMENT *theElement;
  INT nv,nn,ne;
  int k,j;

  nv = ne = nn = 0;

  j = theMG->topLevel;
  for (k=0; k<=j; k++)
  {
    theGrid = theMG->grids[k];

    /* vertices */
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
      ID(theVertex) = nv++;

    /* nodes */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      ID(theNode) = nn++;

    /* elements */
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      ID(theElement) = ne++;
  }

  theMG->vertIdCounter = nv;
  theMG->nodeIdCounter = nn;
  theMG->elemIdCounter = ne;

  return(0);
}


/****************************************************************************/
/*
   LexCompare - Define relation for lexicographic ordering

   SYNOPSIS:
   static int LexCompare (NODE **pnode1, NODE **pnode2);

   PARAMETERS:
   .  pnode1 - first node to compare
   .  pnode2 - second node to compare

   DESCRIPTION:
   This function defines a relation for lexicographic ordering

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT LexCompare (NODE **pnode1, NODE **pnode2)
{
  VERTEX *pv1,*pv2;
  DOUBLE diff[DIM];

  pv1 = MYVERTEX(*pnode1);
  pv2 = MYVERTEX(*pnode2);

  V_DIM_SUBTRACT(CVECT(pv2),CVECT(pv1),diff);
  V_DIM_SCALE(InvMeshSize,diff);

  if (fabs(diff[Order[DIM-1]])<ORDERRES)
  {
                #ifdef __THREEDIM__
    if (fabs(diff[Order[DIM-2]])<ORDERRES)
    {
      if (diff[Order[DIM-3]]>0.0) return (-Sign[DIM-3]);
      else return ( Sign[DIM-3]);
    }
    else
                #endif
    if (diff[Order[DIM-2]]>0.0) return (-Sign[DIM-2]);
    else return ( Sign[DIM-2]);
  }
  else
  {
    if (diff[Order[DIM-1]]>0.0) return (-Sign[DIM-1]);
    else return ( Sign[DIM-1]);
  }
}


/****************************************************************************/
/*
   LinkCompare - Define relation for lexicographic ordering of links

   SYNOPSIS:
   static int LinkCompare (LINK **LinkHandle1, LINK **LinkHandle2)

   PARAMETERS:
   .  LinkHandle1 - first link to compare
   .  LinkHandle2 - second link to compare

   DESCRIPTION:
   This function defines a relation for lexicographic ordering of links

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static int LinkCompare (LINK **LinkHandle1, LINK **LinkHandle2)
{
  INT ID1,ID2;;

  ID1 = ID(NBNODE(*LinkHandle1));
  ID2 = ID(NBNODE(*LinkHandle2));

  if (ID1>ID2)
    return ( 1);
  else
    return (-1);
}

/****************************************************************************/
/*D
   OrderNodesInGrid - reorder double linked 'NODE' list

   SYNOPSIS:
   INT OrderNodesInGrid (GRID *theGrid, const INT *order, const INT *sign);

   PARAMETERS:
   .  theGrid - grid to order
   .  order - precedence of coordinate directions
   .  sign - respective ordering direction
   .  AlsoOrderLinks - if 'TRUE' also order links

   DESCRIPTION:
   This function reorders the double linked 'NODE' list of the grid with
   qsort and order criteria LexCompare(). If specified the 'LINK's are ordered
   corresponding to the 'NODE' order.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   >0 when error occured.
   D*/
/****************************************************************************/

INT OrderNodesInGrid (GRID *theGrid, const INT *order, const INT *sign, INT AlsoOrderLinks)
{
  MULTIGRID *theMG;
  NODE **table,*theNode;
  LINK *theLink,*LinkTable[LINKTABLESIZE];
  INT i,entries,firstID,nl;
  HEAP *theHeap;
  BVP *theBVP;
  BVP_DESC theBVPDesc;

  theMG   = MYMG(theGrid);
  firstID = ID(FIRSTNODE(theGrid));
  entries = NN(theGrid);
  theBVP = MG_BVP(theMG);
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc)) RETURN (1);

  /* calculate the diameter of the bounding rectangle of the domain */
  InvMeshSize = POW2(GLEVEL(theGrid)) * pow(NN(GRID_ON_LEVEL(theMG,0)),1.0/DIM) / BVPD_RADIUS(theBVPDesc);

  /* allocate memory for the node list */
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap);
  if ((table=GetTmpMem(theHeap,entries*sizeof(NODE *)))==NULL)
  {
    ReleaseTmpMem(theHeap);
    PrintErrorMessage('E',"OrderNodesInGrid","ERROR: could not allocate memory from the MGHeap");
    RETURN (2);
  }

  /* fill array of pointers to nodes */
  entries = 0;
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    table[entries++] = theNode;

  /* sort array of pointers */
  Order = order;
  Sign  = sign;
  qsort(table,entries,sizeof(*table),(int (*)(const void *, const void *))LexCompare);

  /* reorder double linked list */
  for (i=0; i<entries-1; i++)
    SUCCN(table[i]) = table[i+1];

  for (i=1; i<entries; i++)
  {
    ID(table[i]) = i+firstID;
    PREDN(table[i]) = table[i-1];
  }
  ID(table[0]) = firstID;

  SUCCN(table[entries-1])  = NULL;
  PREDN(table[0]) = NULL;

        #ifdef ModelP
  LISTPART_FIRSTNODE(theGrid,2) = table[0];
        #else
  FIRSTNODE(theGrid) = table[0];
        #endif
  LASTNODE(theGrid)  = table[entries-1];


  ReleaseTmpMem(theHeap);

  if (!AlsoOrderLinks)
    return (0);

  /* now we also order the links of each node the same way (just using the new IDs) */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    /* fill array for qsort */
    for (nl=0, theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      if (nl>=LINKTABLESIZE)
        RETURN (1);

      LinkTable[nl++] = theLink;
    }
    qsort(LinkTable,nl,sizeof(LINK*),(int (*)(const void *, const void *))LinkCompare);

    /* establish pointer connections */
    NEXT(LinkTable[--nl]) = NULL;
    while (nl>0)
    {
      NEXT(LinkTable[nl-1]) = LinkTable[nl];
      --nl;
    }
    START(theNode) = LinkTable[0];
  }

  return (0);
}

/****************************************************************************/
/*D
   PutAtEndOfList - reorder a given set of elements and put them first in the list

   SYNOPSIS:
   INT PutAtEndOfList (GRID *theGrid, INT cnt, ELEMENT **elemList);

   PARAMETERS:
   .  theGrid - elements are part of that level (not checked)
   .  cnt - number of elements in list
   .  elemList - list of elements to reorder

   DESCRIPTION:
   This function reorders a given set of elements and put them last in the list.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   D*/
/****************************************************************************/

INT PutAtEndOfList (GRID *theGrid, INT cnt, ELEMENT **elemList)
{
  ELEMENT *theElement;
  INT i;

  /* remove all elements from list */
  for (i=0; i<cnt; i++)
  {
    theElement = elemList[i];
    GRID_UNLINK_ELEMENT(theGrid,theElement);
  }

  /* reorder elements locally */
  PREDE(elemList[0]) = NULL;
  SUCCE(elemList[cnt-1]) = NULL;
  for (i=0; i<cnt-1; i++)
  {
    SUCCE(elemList[i]) = elemList[i+1];
    PREDE(elemList[i+1]) = elemList[i];
  }

  /* and insert them at the end of the element list */
  for (i=0; i<cnt; i++)
    GRID_LINK_ELEMENT(theGrid,elemList[i],PrioMaster);

  return(GM_OK);
}

/****************************************************************************/
/*D
   FindNeighborElement - Determine neighbor and side of neighbor that goes back to elem

   SYNOPSIS:
   INT FindNeighborElement (const ELEMENT *theElement, INT Side,
   ELEMENT **theNeighbor, INT *NeighborSide);

   PARAMETERS:
   .  theElement - considered element
   .  Side - side of that element
   .  theNeighbor - handle to neighbor
   .  NeighborSide - number of side of neighbor that goes back to elem

   DESCRIPTION:
   This function determines the neighbor and side of the neighbor that goes back to elem.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 when error occured.
   D*/
/****************************************************************************/

INT FindNeighborElement (const ELEMENT *theElement, INT Side, ELEMENT **theNeighbor, INT *NeighborSide)
{
  INT i;

  /* find neighbor */
  *theNeighbor = NBELEM(theElement,Side);
  if (*theNeighbor == NULL) return (0);

  /* search the side */
  for (i=0; i<SIDES_OF_ELEM(*theNeighbor); i++)
    if (NBELEM(*theNeighbor,i) == theElement)
      break;

  /* found ? */
  if (i<SIDES_OF_ELEM(*theNeighbor))
  {
    *NeighborSide = i;
    return (1);
  }
  return (0);
}


/****************************************************************************/
/*D
   InsertInnerNode - Insert a inner node

   SYNOPSIS:
   NODE *InnerNode (MULTIGRID *theMG, DOUBLE *pos);

   PARAMETERS:
   .  theMG - multigrid structure
   .  pos - array containing position

   DESCRIPTION:
   This function inserts a inner node into level 0.

   RETURN VALUE:
   INT
   .n   pointer to a node if ok
   .n   NULL when error occured.
   D*/
/****************************************************************************/

NODE *InsertInnerNode (MULTIGRID *theMG, DOUBLE *pos)
{
  GRID *theGrid;
  VERTEX *theVertex;
  NODE *theNode;
  INT i;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertInnerNode","only a multigrid with exactly one level can be edited");
    return(NULL);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* create objects */
  theVertex = CreateInnerVertex(theGrid);
  if (theVertex==NULL)
  {
    PrintErrorMessage('E',"InsertInnerNode","cannot create vertex");
    return(NULL);
  }
  theNode = CreateNode(theGrid);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertInnerNode","cannot create node");
    return(NULL);
  }

  /* fill data */
  for (i=0; i<DIM; i++) CVECT(theVertex)[i] = pos[i];
  SETMOVE(theVertex,DIM);

  INDEX(theNode) = 0;
  MYVERTEX(theNode) = theVertex;

  return(theNode);
}

/****************************************************************************/
/*D
   InsertBoundaryNode - Insert a boundary node

   SYNOPSIS:
   NODE *InsertBoundaryNode (MULTIGRID *theMG, BNDP *bndp);

   PARAMETERS:
   .  theMG - multigrid structure
   .  bndp - boundary point descriptor

   DESCRIPTION:
   This function inserts a boundary node into level 0.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

NODE *InsertBoundaryNode (MULTIGRID *theMG, BNDP *bndp)
{
  GRID *theGrid;
  NODE *theNode;
  VERTEX *theVertex;
  INT move;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertBoundaryNode",
                      "only a multigrid with exactly one level can be edited");
    return(NULL);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* create objects */
  theVertex = CreateBoundaryVertex(theGrid);
  if (theVertex==NULL)
  {
    BNDP_Dispose(MGHEAP(MYMG(theGrid)),bndp);
    PrintErrorMessage('E',"InsertBoundaryNode","cannot create vertex");
    return(NULL);
  }
  if (BNDP_Global(bndp,CVECT(theVertex)))
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }

  if (BNDP_BndPDesc(bndp,&move))
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  SETMOVE(theVertex,move);
  V_BNDP(theVertex) = bndp;

  theNode = CreateNode(theGrid);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertBoundaryNode","cannot create node");
    return(NULL);
  }
  MYVERTEX(theNode) = theVertex;
  TOPNODE(theVertex) = theNode;

  /* fill data into node/vertex */
  INDEX(theNode) = 0;

  PRINTDEBUG(dom,1,("  ipn %ld nd %x bndp %x \n",
                    ID(theNode),theNode,V_BNDP(theVertex)));

  return(theNode);
}

/****************************************************************************/
/*D
   DeleteNode - Delete a node

   SYNOPSIS:
   INT DeleteNode (MULTIGRID *theMG, NODE *theNode);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theNode - node to delete

   DESCRIPTION:
   This function deletes a node from level 0.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT DeleteNode (MULTIGRID *theMG, NODE *theNode)
{
  GRID *theGrid;
  VERTEX *theVertex;
  ELEMENT *theElement;
  INT i;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"DeleteNode","only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  if (theNode==NULL)
  {
    PrintErrorMessage('E',"DeleteNode","node not found");
    RETURN(GM_ERROR);
  }

  /* check corner */
  theVertex = MYVERTEX(theNode);
  if (MOVE(theVertex)==0)
  {
    PrintErrorMessage('E',"DeleteNode","corners cannot be deleted");
    RETURN(GM_ERROR);
  }

  /* check if some element needs that node */
  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      if (CORNER(theElement,i)==theNode)
      {
        PrintErrorMessage('E',"DeleteNode","there is an element needing that node");
        RETURN(GM_ERROR);
      }

  /* now allowed to delete */
  DisposeNode(theGrid,theNode);

  return(GM_OK);
}

/****************************************************************************/
/*D
   DeleteNodeWithID - Delete the node with id

   SYNOPSIS:
   INT DeleteNodeWithID (MULTIGRID *theMG, INT id);

   PARAMETERS:
   .  theMG - multigrid structure
   .  id - id of node to delete

   DESCRIPTION:
   This function deletes the node with id `id`.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT DeleteNodeWithID (MULTIGRID *theMG, INT id)
{
  GRID *theGrid;
  NODE *theNode;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"DeleteNodeWithID","only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* find node */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    if (ID(theNode)==id) break;
  if (theNode==NULL)
  {
    PrintErrorMessage('E',"DeleteNodeWithID","node not found");
    RETURN(GM_ERROR);
  }
  return (DeleteNode(theMG,theNode));
}

/****************************************************************************/
/*D
   FindFather - Find the new father element

   SYNOPSIS:
   ELEMENT *FindFather(VERTEX *vptr);

   PARAMETERS:
   .  vptr - Pointer to 'VERTEX' whose father element is to be found.

   DESCRIPTION:
   This function finds the new father element of the given vertex.
   It assumes that the  new father is one of the neighbors of the
   old father element. The current father of 'vptr' is not changed.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to an element
   .n     NULL if none or no correct father is found or vertex is level 0
   D*/
/****************************************************************************/

ELEMENT *FindFather (VERTEX *theVertex)
{
  ELEMENT *theElement;
  INT i;

  theElement = VFATHER(theVertex);

  if (theElement == NULL)
    return(NULL);

  if (PointInElement(CVECT(theVertex),theElement))
    return(theElement);

  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    if (PointInElement(CVECT(theVertex),NBELEM(theElement,i)))
      return(NBELEM(theElement,i));

  if (i == SIDES_OF_ELEM(theElement))
    if (OBJT(theVertex) == BVOBJ)
      return(theElement);

  return(NULL);
}

/****************************************************************************/
/*D
   MoveMidNode - set new position for a midnode

   SYNOPSIS:
   INT MoveMidNode (MULTIGRID *theMG, NODE *theNode, DOUBLE lambda);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  lambda - parameter on the edge

   DESCRIPTION:
   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT MoveMidNode (MULTIGRID *theMG, NODE *theNode, DOUBLE lambda)
{
  ELEMENT *theElement;
  NODE *Node0,*Node1;
  VERTEX *theVertex;
  BNDP *bndp;
  DOUBLE *x[MAX_CORNERS_OF_ELEM],*global,*local;
  DOUBLE_VECTOR bnd_global;
  DOUBLE diff;
  INT n,k,co0,co1,edge;

  if ((lambda<0) || (lambda>1)) {
    PrintErrorMessage('E',"MoveMidNode","lambda not in range (0,1)");
    return(GM_ERROR);
  }
  if (NFATHER(theNode) != NULL) {
    PrintErrorMessage('E',"MoveMidNode","node not a midnode");
    return(GM_ERROR);
  }
  theVertex = MYVERTEX(theNode);
  theElement = VFATHER(theVertex);
  edge = ONEDGE(theVertex);
  co0 = CORNER_OF_EDGE(theElement,edge,0);
  co1 = CORNER_OF_EDGE(theElement,edge,1);
  Node0 = CORNER(theElement,co0);
  Node1 = CORNER(theElement,co1);
  global = CVECT(theVertex);
  local = LCVECT(theVertex);
  V_DIM_LINCOMB((1.0-lambda), CVECT(MYVERTEX(Node0)),
                lambda      , CVECT(MYVERTEX(Node1)), global);
  V_DIM_LINCOMB((1.0-lambda), LOCAL_COORD_OF_ELEM(theElement,co0),
                lambda      , LOCAL_COORD_OF_ELEM(theElement,co1),
                local);
  if (OBJT(theVertex) == BVOBJ) {
    if (BNDP_Dispose(MGHEAP(theMG),V_BNDP(theVertex)))
      return(GM_ERROR);
    bndp = BNDP_CreateBndP(MGHEAP(theMG),V_BNDP(MYVERTEX(Node0)),
                           V_BNDP(MYVERTEX(Node1)),lambda);
    if (bndp == NULL)
      return(GM_ERROR);
    V_BNDP(theVertex) = bndp;
    if (BNDP_Global(bndp,bnd_global))
      RETURN(GM_ERROR);
    V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
    if (diff > MAX_PAR_DIST) {
      SETMOVED(theVertex,1);
      CORNER_COORDINATES(theElement,n,x);
      UG_GlobalToLocal(n,(const DOUBLE **)x,global,local);
    }
  }

  /* Warning: O(n) Operation! */
  for(k=LEVEL(theNode)+1; k<=TOPLEVEL(theMG); k++)
    for (theVertex=FIRSTVERTEX(GRID_ON_LEVEL(theMG,k));
         theVertex!=NULL; theVertex=SUCCV(theVertex))
      if ((OBJT(theVertex) != BVOBJ)) {
        CORNER_COORDINATES(VFATHER(theVertex),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(theVertex),CVECT(theVertex));
      }

  return(GM_OK);
}

/****************************************************************************/
/*D
   MoveCenterNode - set new position for a centernode

   SYNOPSIS:
   INT MoveCenterNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  lambda - local coordinate in the father element

   DESCRIPTION:
   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT MoveCenterNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda)
{
  VERTEX *theVertex;
  ELEMENT *theElement;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE_VECTOR oldPos,newPos;
  INT n,k;

  if (NFATHER(theNode) != NULL) {
    PrintErrorMessage('E',"MoveCenterNode","node not a sidenode");
    return(GM_ERROR);
  }
  theVertex = MYVERTEX(theNode);
  theElement = VFATHER(theVertex);
  ASSERT(theElement != NULL);
  if (OBJT(theVertex) == BVOBJ) {
    PrintErrorMessage('E',"MoveCenterNode","no inner node");
    return(GM_ERROR);
  }
  CORNER_COORDINATES(theElement,n,x);
  LOCAL_TO_GLOBAL(n,x,lambda,newPos);
  V_DIM_COPY(CVECT(theVertex),oldPos);
  V_DIM_COPY(newPos,CVECT(theVertex));
  theElement = FindFather(theVertex);
  if (theElement == NULL) {
    PrintErrorMessage('W',"MoveCenterNode","cannot find father element");
    V_DIM_COPY(oldPos,CVECT(theVertex));
    return(GM_ERROR);
  }
  else {
    V_DIM_COPY(lambda,LCVECT(theVertex));
    VFATHER(theVertex) = theElement;
  }

  /* Warning: O(n) Operation! */
  for(k=LEVEL(theNode)+1; k<=TOPLEVEL(theMG); k++)
    for (theVertex=FIRSTVERTEX(GRID_ON_LEVEL(theMG,k));
         theVertex!=NULL; theVertex=SUCCV(theVertex))
      if ((OBJT(theVertex) != BVOBJ))
      {
        CORNER_COORDINATES(VFATHER(theVertex),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(theVertex),CVECT(theVertex));
      }

  return(GM_OK);
}

/****************************************************************************/
/*D
   MoveSideNode - set new position for a sidenode

   SYNOPSIS:
   INT MoveCenterNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  lambda - local coordinate on the father element side

   DESCRIPTION:
   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

#ifdef __THREEDIM__
INT MoveSideNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda)
{
  ELEMENT *theElement;
  NODE *Node[MAX_CORNERS_OF_SIDE];
  VERTEX *theVertex;
  BNDP *bndp;
  DOUBLE *x[MAX_CORNERS_OF_ELEM],*global,*local;
  DOUBLE_VECTOR bnd_global;
  DOUBLE diff;
  INT n,m,k,i,co[MAX_CORNERS_OF_SIDE],side;

  if ((lambda[0]<0) || (lambda[0]>1) ||(lambda[1]<0) || (lambda[1]>1)) {
    PrintErrorMessage('E',"MoveSideNode","lambda not in range (0,1)^2");
    return(GM_ERROR);
  }
  if (NFATHER(theNode) != NULL) {
    PrintErrorMessage('E',"MoveSideNode","node not a sidenode");
    return(GM_ERROR);
  }
  theVertex = MYVERTEX(theNode);
  theElement = VFATHER(theVertex);
  side = ONSIDE(theVertex);
  m = CORNERS_OF_SIDE(theElement,side);
  if (m  != 4) {
    PrintErrorMessage('E',"MoveSideNode","node not a sidenode");
    return(GM_ERROR);
  }
  global = CVECT(theVertex);
  local = LCVECT(theVertex);
  V_DIM_CLEAR(global);
  V_DIM_CLEAR(local);
  for (i=0; i<m; i++) {
    co[i] = CORNER_OF_SIDE(theElement,side,i);
    Node[i] = CORNER(theElement,co[i]);

  }
  V_DIM_LINCOMB((1.0-lambda[0])*(1.0-lambda[1]),
                CVECT(MYVERTEX(Node[0])),1.0,global,global);
  V_DIM_LINCOMB(lambda[0]*(1.0-lambda[1]),
                CVECT(MYVERTEX(Node[1])),1.0,global,global);
  V_DIM_LINCOMB(lambda[0]*lambda[1],
                CVECT(MYVERTEX(Node[2])),1.0,global,global);
  V_DIM_LINCOMB((1.0-lambda[0])*lambda[1],
                CVECT(MYVERTEX(Node[3])),1.0,global,global);
  V_DIM_LINCOMB((1.0-lambda[0])*(1.0-lambda[1]),
                LOCAL_COORD_OF_ELEM(theElement,co[0]),1.0,local,local);
  V_DIM_LINCOMB(lambda[0]*(1.0-lambda[1]),
                LOCAL_COORD_OF_ELEM(theElement,co[1]),1.0,local,local);
  V_DIM_LINCOMB(lambda[0]*lambda[1],
                LOCAL_COORD_OF_ELEM(theElement,co[2]),1.0,local,local);
  V_DIM_LINCOMB((1.0-lambda[0])*lambda[1],
                LOCAL_COORD_OF_ELEM(theElement,co[3]),1.0,local,local);

  if (OBJT(theVertex) == BVOBJ) {
    if (BNDP_Dispose(MGHEAP(theMG),V_BNDP(theVertex)))
      return(GM_ERROR);
    bndp = BNDS_CreateBndP(MGHEAP(theMG),ELEM_BNDS(theElement,side),lambda);
    if (bndp == NULL)
      return(GM_ERROR);
    V_BNDP(theVertex) = bndp;
    if (BNDP_Global(bndp,bnd_global))
      RETURN(GM_ERROR);
    V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
    if (diff > MAX_PAR_DIST) {
      SETMOVED(theVertex,1);
      CORNER_COORDINATES(theElement,n,x);
      UG_GlobalToLocal(n,(const DOUBLE **)x,global,local);
    }
  }

  /* Warning: O(n) Operation! */
  for(k=LEVEL(theNode)+1; k<=TOPLEVEL(theMG); k++)
    for (theVertex=FIRSTVERTEX(GRID_ON_LEVEL(theMG,k));
         theVertex!=NULL; theVertex=SUCCV(theVertex))
      if ((OBJT(theVertex) != BVOBJ)) {
        CORNER_COORDINATES(VFATHER(theVertex),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(theVertex),CVECT(theVertex));
      }

  return(GM_OK);
}
#endif

/****************************************************************************/
/*D
   MoveNode - Let user enter a new position for an inner node

   SYNOPSIS:
   INT MoveNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *newPos)

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  newPos - global coordinate for new position

   DESCRIPTION:
   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT MoveNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *newPos)
{
  VERTEX *theVertex;
  ELEMENT *theElement;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE_VECTOR oldPos;
  INT n,k;

  /* set k (and theNode) to the level where the node
     appears the first time */
  while (NFATHER(theNode)!=NULL)
    theNode = NFATHER(theNode);
  theVertex  = MYVERTEX(theNode);
  if (OBJT(theVertex) == BVOBJ)
  {
    PrintErrorMessage('E',"MoveNode","no inner node passed");
    return(GM_ERROR);
  }

  if (LEVEL(theNode) > 0)
  {
    V_DIM_COPY(CVECT(theVertex),oldPos);
    V_DIM_COPY(newPos,CVECT(theVertex));
    theElement = FindFather(theVertex);
    if (theElement == NULL)
    {
      PrintErrorMessage('W',"MoveNode",
                        "cannot find father element");
      V_DIM_COPY(oldPos,CVECT(theVertex));
      return(GM_ERROR);
    }
    else
    {
      CORNER_COORDINATES(theElement,n,x);
      UG_GlobalToLocal(n,(const DOUBLE **)x,newPos,LCVECT(theVertex));
      VFATHER(theVertex) = theElement;
    }
  }
  else
    V_DIM_COPY(newPos,CVECT(theVertex));

  /* Warning: O(n) Operation! */
  for(k=LEVEL(theNode)+1; k<=TOPLEVEL(theMG); k++)
    for (theVertex=FIRSTVERTEX(GRID_ON_LEVEL(theMG,k));
         theVertex!=NULL; theVertex=SUCCV(theVertex))
      if ((OBJT(theVertex) != BVOBJ))
      {
        CORNER_COORDINATES(VFATHER(theVertex),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(theVertex),CVECT(theVertex));
      }

  return(GM_OK);
}

INT GetMidNodeParam (NODE * theNode, DOUBLE *lambda)
{
  *lambda = 0.5;
  return(0);
}

/****************************************************************************/
/*D
   InsertElement - Insert an element

   SYNOPSIS:
   ELEMENT *InsertElement (MULTIGRID *theMG, INT n, NODE **Node);

   PARAMETERS:
   .  theMG - multigrid structure
   .  tag - type of element to insert
   .  idList - array of ids of corner nodes

   DESCRIPTION:
   This function inserts an element into level 0.

   RETURN VALUE:
   INT
   .n   pointer to an element if ok
   .n   NULL when error occured.
   D*/
/****************************************************************************/

#ifdef __TWODIM__

static INT CheckOrientation (INT n, VERTEX **vertices)
{
  int i;
  DOUBLE x1,x2,y1,y2;

  for (i=0; i<n; i++)
  {
    x1 = XC(vertices[(i+1)%n])-XC(vertices[i]);
    x2 = XC(vertices[(i+n-1)%n])-XC(vertices[i]);
    y1 = YC(vertices[(i+1)%n])-YC(vertices[i]);
    y2 = YC(vertices[(i+n-1)%n])-YC(vertices[i]);
    if (vp(x1,y1,x2,y2)<SMALL_C)
    {
      return(0);
    }
  }
  return(1);
}

#define SWAP_IJ(a,i,j,t)                        {t = a[i]; a[i] = a[j]; a[j] = t;}
#endif

#ifdef __THREEDIM__
static INT CheckOrientation (INT n, VERTEX **vertices)
{
  DOUBLE_VECTOR diff[3],rot;
  DOUBLE det;
  INT i;

  /* TODO: this case */
  if (n == 8 || n==6)
    return(0);

  for (i=1; i<n; i++)
    V3_SUBTRACT(CVECT(vertices[i]),CVECT(vertices[0]),diff[i-1]);
  V3_VECTOR_PRODUCT(diff[0],diff[1],rot);
  V3_SCALAR_PRODUCT(rot,diff[2],det);

  if (det < 0.0)
    return(1);

  return(0);
}
#endif

ELEMENT *InsertElement (MULTIGRID *theMG, INT n, NODE **Node, ELEMENT **ElemList, INT *NbgSdList)
{
  GRID             *theGrid;
  INT i,j,k,m,found,num,tag,ElementType;
  INT NeighborSide[MAX_SIDES_OF_ELEM];
  NODE             *sideNode[MAX_CORNERS_OF_SIDE],*NeighborNode;
  VERTEX           *Vertex[MAX_CORNERS_OF_ELEM],*sideVertex[MAX_CORNERS_OF_SIDE];
  ELEMENT          *theElement,*Neighbor[MAX_SIDES_OF_ELEM];
  BNDS         *bnds[MAX_SIDES_OF_ELEM];
  BNDP         *bndp[MAX_CORNERS_OF_ELEM];
        #ifdef __TWODIM__
  VERTEX           *theVertex;
  NODE             *theNode;
        #endif

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertElement",
                      "only a multigrid with exactly one level can be edited");
    return(NULL);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* check parameters */
    #ifdef __TWODIM__
  if ((n!=TRIANGLE)&&(n!=QUADRILATERAL))
  {
    PrintErrorMessage('E',"InsertElement",
                      "only triangles and quadrilaterals allowed in 2D");
    return(NULL);
  }
  tag = n;
        #endif

    #ifdef __THREEDIM__
  if (n == 4)
    tag = TETRAHEDRON;
  else if ( n == 5)
    tag = PYRAMID;
  else if ( n == 6)
    tag = PRISM;
  else if ( n == 8)
    tag = HEXAHEDRON;
  else
  {
    PrintErrorMessage('E',"InsertElement",
                      "only tetrahedrons, pyramids and hexahedrons are allowed in the 3D coarse grid");
    return(NULL);
  }
    #endif

  /* init vertices */
  for (i=0; i<n; i++)
    Vertex[i] = MYVERTEX(Node[i]);

    #ifdef __TWODIM__
  /* find orientation */
  if (!CheckOrientation(n,Vertex))
  {
    /* flip order */
    SWAP_IJ(Node,   0,n/2,theNode);
    SWAP_IJ(Vertex,0,n/2,theVertex);

    if (!CheckOrientation(n,Vertex))
    {
      /* this was the only possibility for a triangle: so is a nonconvex quadrilateral */
      /* interchange first two nodes and try again */
      SWAP_IJ(Node,   0,1,theNode);
      SWAP_IJ(Vertex,0,1,theVertex);
      if (!CheckOrientation(n,Vertex))
      {
        /* flip order */
        SWAP_IJ(Node,   0,n/2,theNode);
        SWAP_IJ(Vertex,0,n/2,theVertex);
        if (!CheckOrientation(n,Vertex))
        {
          /* flip order back */
          SWAP_IJ(Node,   0,n/2,theNode);
          SWAP_IJ(Vertex,0,n/2,theVertex);
          /* interchange second two nodes and try again */
          SWAP_IJ(Node,   1,2,theNode);
          SWAP_IJ(Vertex,1,2,theVertex);
          if (!CheckOrientation(n,Vertex))
          {
            /* flip order */
            SWAP_IJ(Node,   0,n/2,theNode);
            SWAP_IJ(Vertex,0,n/2,theVertex);
            if (!CheckOrientation(n,Vertex))
            {
              PrintErrorMessage('E',"InsertElement",
                                "cannot find orientation");
              return(NULL);
            }
          }
        }
      }
    }
  }
        #endif

    #ifdef __THREEDIM__
  if (CheckOrientation (n,Vertex))
  {
    sideNode[0] = Node[0];
    sideVertex[0] = Vertex[0];
    Node[0] = Node[1];
    Vertex[0] = Vertex[1];
    Node[1] = sideNode[0];
    Vertex[1] = sideVertex[0];
  }
        #endif

  /* init pointers */
  for (i=0; i<SIDES_OF_REF(n); i++)
  {
    Neighbor[i] = NULL;
    bnds[i] = NULL;
  }

  /* compute side information (theSeg[i]==NULL) means inner side */
  ElementType = IEOBJ;
  for (i=0; i<SIDES_OF_REF(n); i++)
  {
    m = CORNERS_OF_SIDE_REF(n,i);
    for(j=0; j<m; j++ )
    {
      k = CORNER_OF_SIDE_REF(n,i,j);
      sideNode[j] = Node[k];
      sideVertex[j] = Vertex[k];
    }
    found = 0;
    for(j=0; j<m; j++ )
    {
      if( OBJT(sideVertex[j]) == IVOBJ ) found = 1;
    }
    if( found ) continue;

    /* all vertices of side[i] are on the boundary now */

    /* We now assume, that side[i] is on the boundary if and only if */
    /* there is one boundary segment containing the three nodes.        */
    /* That means, one should not define elements at the boundary	*/
    /* with a boundary side covering more than one segment.			*/

    for (j=0; j<m; j++)
      bndp[j] = V_BNDP(sideVertex[j]);

    bnds[i] = BNDP_CreateBndS(MGHEAP(theMG),bndp,m);

    if (bnds[i] != NULL)
      ElementType = BEOBJ;
  }

  /* find neighboring elements */
  if (ElemList == NULL)
    /* for all sides of the element to be created */
    for (i=0; i<SIDES_OF_REF(n); i++)
    {
      for(j=0; j<CORNERS_OF_SIDE_REF(n,i); j++ )
        sideNode[j] = Node[CORNER_OF_SIDE_REF(n,i,j)];

      /* for all neighbouring elements allready inserted */
      for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
           theElement=SUCCE(theElement))
      {
        /* for all sides of the neighbour element */
        for (j=0; j<SIDES_OF_ELEM(theElement); j++)
        {
          num = 0;
          /* for all corners of the side of the neighbour */
          for (m=0; m<CORNERS_OF_SIDE(theElement,j); m++)
          {
            NeighborNode = CORNER(theElement,
                                  CORNER_OF_SIDE(theElement,j,m));
            /* for all corners of the side of the
                   element to be created */
            for (k=0; k<CORNERS_OF_SIDE_REF(n,i); k++)
              if(NeighborNode==sideNode[k])
              {
                num++;
                break;
              }
          }
          if(num==CORNERS_OF_SIDE_REF(n,i))
          {
            if (NBELEM(theElement,j)!=NULL)
            {
              PrintErrorMessage('E',"InsertElement",
                                "neighbor relation inconsistent");
              return(NULL);
            }
            Neighbor[i] = theElement;
            NeighborSide[i] = j;
          }
        }
      }
    }
  else {
    for (i=0; i<SIDES_OF_REF(n); i++) Neighbor[i] = ElemList[i];
    if (NbgSdList != NULL)
      for (i=0; i<SIDES_OF_REF(n); i++)
        NeighborSide[i] = NbgSdList[i];
  }

  /* create element */
  theElement = CreateElement(theGrid,tag,ElementType,Node,NULL);
  if (theElement==NULL)
  {
    PrintErrorMessage('E',"InsertElement","cannot allocate element");
    return(NULL);
  }

  IFDEBUG(dom,1)
  {
    DOUBLE *x[MAX_CORNERS_OF_ELEM],global[DIM],fac,diam;

    V_DIM_CLEAR(global);
    for (i=0; i<n; i++)
    {
      x[i] = CVECT(Vertex[i]);
      V_DIM_ADD(x[i],global,global);
    }
    fac = 1.0 / n;
    V_DIM_SCALE(fac,global);
    diam = 0.0;
    for (i=0; i<n; i++)
    {
      V_DIM_EUKLIDNORM_OF_DIFF(x[i],global,fac);
      diam = MAX(fac,diam);
    }
    printf(" created ID(Elem)=%d midPoint %5.2f %5.2f %5.2f diam %5.2f\n",
           ID(theElement),global[0],global[1],global[2],diam);
  }
  ENDDEBUG

  /* create element sides if necessary */
  if (OBJT(theElement)==BEOBJ)
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      SET_BNDS(theElement,i,bnds[i]);

  /* fill element data */
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    SET_NBELEM(theElement,i,Neighbor[i]);
    if (Neighbor[i]!=NULL)
    {
      if (NbgSdList == NULL)
        NeighborSide[i] = SideOfNbElement(theElement,i);
      if (NeighborSide[i] >= MAX_SIDES_OF_ELEM) {
        PrintErrorMessage('E',"InsertElement",
                          "neighbor relation inconsistent");
        return(NULL);
      }
      SET_NBELEM(Neighbor[i],NeighborSide[i],theElement);
            #ifdef __THREEDIM__
      if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
        if (DisposeDoubledSideVector(theGrid,Neighbor[i],
                                     NeighborSide[i],theElement,i))
          return(NULL);
            #endif
    }
  }

  SETNSONS(theElement,0);
  SET_SON(theElement,0,NULL);
  SET_EFATHER(theElement,NULL);
  SETECLASS(theElement,RED_CLASS);

  /* create connection to other elements. ATTENTION: this function is O(n)*/
  if (InsertedElementCreateConnection(theGrid,theElement))
  {
    PrintErrorMessage('E',"InsertElement",
                      "could not create algebra connections");
    DisposeElement (theGrid,theElement,TRUE);
    return(NULL);
  }

  return(theElement);
}

/****************************************************************************/
/*D
   InsertElementFromIDs - Insert element with node ids

   SYNOPSIS:
   ELEMENT *InsertElementFromIDs (MULTIGRID *theMG, INT n, INT *idList);

   PARAMETERS:
   .  theMG - multigrid structure
   .  n - number of nodes in node id list
   .  idList - ids of the nodes

   DESCRIPTION:
   This function inserts an element with nodes that have the ids
   given in `idList`,  on level 0.

   RETURN VALUE:
   INT
   .n   pointer to an element if ok
   .n   NULL when error occured.
   D*/
/****************************************************************************/

ELEMENT *InsertElementFromIDs (MULTIGRID *theMG, INT n, INT *idList)
{
  GRID *theGrid;
  NODE *Node[MAX_CORNERS_OF_ELEM],*theNode;
  INT i,j,found;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertElementFromIDs","only a multigrid with exactly one level can be edited");
    return(NULL);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* check data */
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++)
      if (idList[i]==idList[j])
      {
        PrintErrorMessage('E',"InsertElementFromIDs",
                          "nodes must be pairwise different");
        return(NULL);
      }

  /* init data */
  for (i=0; i<n; i++)
    Node[i] = NULL;

  /* find nodes */
  found = 0;
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (i=0; i<n; i++)
    {
      if ((Node[i]==NULL)&&(ID(theNode)==idList[i]))
      {
        Node[i] = theNode;
        found++;
      }
    }
    if (found==n) break;
  }
  if (found!=n)
  {
    PrintErrorMessage('E',"InsertElementFromIDs",
                      "could not find all nodes");
    return(NULL);
  }

  return (InsertElement(theMG,n,Node,NULL,NULL));
}

/****************************************************************************/
/*D
   DeleteElement - Delete an element

   SYNOPSIS:
   INT DeleteElement (MULTIGRID *theMG, ELEMENT *theElement);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theElement - element to delete

   DESCRIPTION:
   This function deletes an element from level 0.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT DeleteElement (MULTIGRID *theMG, ELEMENT *theElement) /* 3D VERSION */
{
  GRID *theGrid;
  ELEMENT *theNeighbor;
  INT i,j,found;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"DeleteElement",
                      "only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* delete pointers in neighbors */
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    theNeighbor = NBELEM(theElement,i);
    if (theNeighbor!=NULL)
    {
      found = 0;
      for (j=0; j<SIDES_OF_ELEM(theNeighbor); j++)
        if (NBELEM(theNeighbor,j)==theElement)
        {
          found++;
          SET_NBELEM(theNeighbor,j,NULL);
        }
      if (found!=1) RETURN(GM_ERROR);
    }
  }

  /* delete element now */
  DisposeElement(theGrid,theElement,TRUE);

  return(GM_OK);
}

INT DeleteElementWithID (MULTIGRID *theMG, INT id)
{
  GRID *theGrid;
  ELEMENT *theElement;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"DeleteElementWithId","only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* find element */
  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    if (ID(theElement)==id) break;
  if (theElement==NULL)
  {
    PrintErrorMessage('E',"DeleteElementWithId","element not found");
    RETURN(GM_ERROR);
  }

  return (DeleteElement(theMG,theElement));
}

/****************************************************************************/
/*D
   InsertMesh - Insert a mesh described by the domain

   SYNOPSIS:
   INT InsertMesh (MULTIGRID *theMG, MESH *theMesh);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theMesh - mesh structure

   DESCRIPTION:
   This function inserts all nodes and elements given by the mesh.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT InsertMesh (MULTIGRID *theMG, MESH *theMesh)
{
  GRID *theGrid;
  ELEMENT **EList,*nbList[MAX_SIDES_OF_ELEM],*theElement;
  NODE **NList,*Nodes[MAX_CORNERS_OF_ELEM],*theNode;
  INT sid,i,nel,nnd,k,n,nid,eid;

        #ifdef ModelP
  if (me!=master)
    return(GM_OK);
        #endif

  if (theMesh == NULL) return(GM_OK);
  if (theMesh->nElements == NULL)
  {
    for (i=0; i<theMesh->nBndP; i++)
      if (InsertBoundaryNode(theMG,theMesh->theBndPs[i]) == NULL)
        return(GM_ERROR);

    for (i=0; i<theMesh->nInnP; i++)
      if (InsertInnerNode(theMG,theMesh->Position[i]) == NULL)
        return(GM_ERROR);
    return(GM_OK);
  }

  theGrid = GRID_ON_LEVEL(theMG,0);
  nnd = theMesh->nBndP + theMesh->nInnP;
  NList = (NODE **) GetTmpMem(MGHEAP(theMG),nnd*sizeof(NODE *));
  if (NList == NULL) return(GM_ERROR);

  nid=0;
  for (i=0; i<theMesh->nBndP; i++)
  {
    theNode = InsertBoundaryNode(theMG,theMesh->theBndPs[i]);
    ID(theNode) = nid++;
    if (theNode == NULL) return(GM_ERROR);
    NList[i] = theNode;
  }
  for (i=0; i<theMesh->nInnP; i++)
  {
    theNode = InsertInnerNode(theMG,theMesh->Position[i]);
    ID(theNode) = nid++;
    if (theNode == NULL) return(GM_ERROR);
    NList[i+theMesh->nBndP] = theNode;
  }

  EList = NULL;
  if (theMesh->nbElements != NULL)
  {
    theGrid = GRID_ON_LEVEL(theMG,0);
    nel = 0;
    for (sid=1; sid<=theMesh->nSubDomains; sid++)
      nel += theMesh->nElements[sid];
    EList = (ELEMENT **) GetTmpMem(MGHEAP(theMG),nel*sizeof(ELEMENT *));
    if (EList != NULL)
      for (i=0; i<nel; i++) EList[i] = NULL;
  }

  nel = 0,eid=0;
  for (sid=1; sid<=theMesh->nSubDomains; sid++)
    for (i=0; i<theMesh->nElements[sid]; i++)
    {
      n = theMesh->Element_corners[sid][i];
      for (k=0; k<n; k++)
        Nodes[k] = NList[theMesh->Element_corner_ids[sid][i][k]];
      if (EList != NULL)
        for (k=0; k<SIDES_OF_REF(n); k++)
          if (theMesh->nbElements[sid][i][k]>=0)
            nbList[k] = EList[theMesh->nbElements[sid][i][k]];
          else
            nbList[k] = NULL;

      theElement = InsertElement (theMG,n,Nodes,nbList,NULL);
      ID(theElement) = eid++;
      if (theElement == NULL) return(GM_ERROR);
      if (EList != NULL) EList[nel++] = theElement;
    }

  return(GM_OK);
}

/****************************************************************************/
/*D
   FindNodeFromId - Find a node with given id

   SYNOPSIS:
   NODE *FindNodeFromId (GRID *theGrid, INT id);

   PARAMETERS:
   .  theGrid - grid level to search.

   DESCRIPTION:
   This function finds a node with given id.

   RETURN VALUE:
   NODE *
   .n   pointer to that NODE
   .n   NULL if not found.
   D*/
/****************************************************************************/

NODE *FindNodeFromId (GRID *theGrid, INT id)
{
  NODE *theNode;

  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    if (ID(theNode)==id) return(theNode);

  return(NULL);
}


/****************************************************************************/
/*D
   FindNodeFromPosition - Find node from position

   SYNOPSIS:
   NODE *FindNodeFromPosition (GRID *theGrid, DOUBLE *pos, DOUBLE *tol);

   PARAMETERS:
   .  theGrid - grid level to search
   .  pos - given position
   .  tol - tolerance to accept

   DESCRIPTION:
   This function finds the first node within `tol` from `pos` in 1-norm.

   RETURN VALUE:
   NODE *
   .n   pointer to NODE
   .n   NULL if not found.
   D*/
/****************************************************************************/

NODE *FindNodeFromPosition (GRID *theGrid, DOUBLE *pos, DOUBLE *tol)
{
  NODE *theNode;
  int i,found;

  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    found = 1;
    for (i=0; i<DIM; i++)
      if (fabs(pos[i]-CVECT(MYVERTEX(theNode))[i])>=tol[i]) {found=0; break;}
    if (found) return(theNode);
  }

  return(NULL);
}

/****************************************************************************/
/*D
   FindVectorFromPosition - Find vector from position

   SYNOPSIS:
   VECTOR *FindVectorFromPosition (GRID *theGrid, DOUBLE *pos, DOUBLE *tol);

   PARAMETERS:
   .  theGrid - grid level to search
   .  pos - given position
   .  tol - tolerance to accept

   DESCRIPTION:
   This function finds the first vector within `tol` from `pos` in 1-norm.

   RETURN VALUE:
   NODE *
   .n   pointer to NODE
   .n   NULL if not found.
   D*/
/****************************************************************************/

VECTOR *FindVectorFromPosition (GRID *theGrid, DOUBLE *pos, DOUBLE *tol)
{
  VECTOR *theVector;
  DOUBLE_VECTOR vpos;
  int i,found;

  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    VectorPosition(theVector,vpos);
    found = 1;
    for (i=0; i<DIM; i++)
      if (fabs(pos[i]-vpos[i])>=tol[i]) {found=0; break;}
    if (found) return(theVector);
  }

  return(NULL);
}

/****************************************************************************/
/*D
   FindElementFromId - Find element with id

   SYNOPSIS:
   ELEMENT *FindElementFromId (GRID *theGrid, INT id);

   PARAMETERS:
   .  theGrid - grid level to search
   .  id - id to search

   DESCRIPTION:
   This function finds an element with the identification `id`.

   RETURN VALUE:
   ELEMENT *
   .n   pointer to that ELEMENT
   .n   NULL if not found.
   D*/
/****************************************************************************/

ELEMENT *FindElementFromId (GRID *theGrid, INT id)
{
  ELEMENT *theElement;

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    if (ID(theElement)==id) return(theElement);

  return(NULL);
}

/****************************************************************************/
/*D
   PointInElement - Determine whether point is contained in element

   SYNOPSIS:
   INT PointInElement (const DOUBLE *x, const ELEMENT *theElement);

   PARAMETERS:
   .  x - coordinates of given point
   .  theElement - element to scan

   DESCRIPTION:
   This function determines whether a given point specified by coordinates `x`
   is contained in an element.

   RETURN VALUE:
   INT
   .n   0 an error occurred
   .n   1 point is contained in the element
   .n   2 point is nearly on one side of the the element
   .n   3 point is nearly on one edge of the the element
   .n   4 point is nearly one of the corners of the the element
   .n   5 point is not contained in the element
   D*/
/****************************************************************************/

#ifdef __TWODIM__
INT PointInElement (const DOUBLE *x, const ELEMENT *theElement) /* 2D version */
{
  COORD_POINT point[MAX_CORNERS_OF_ELEM],thePoint;
  int n,i;

  /* check element */
  if (theElement==NULL) return(0);

  /* load geometrical data of the corners */
  n = CORNERS_OF_ELEM(theElement);
  for (i=0; i<n; i++)
  {
    point[i].x = XC(MYVERTEX(CORNER(theElement,i)));
    point[i].y = YC(MYVERTEX(CORNER(theElement,i)));
  }
  thePoint.x = (DOUBLE)x[0];
  thePoint.y = (DOUBLE)x[1];

  return(PointInPolygon(point,n,thePoint));
}
#endif

#ifdef __THREEDIM__
INT PointInElement (const DOUBLE *global, const ELEMENT *theElement)
{
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE_VECTOR a,b,rot;
  DOUBLE det;
  INT n,i;

  CORNER_COORDINATES(theElement,n,x);

  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    V3_SUBTRACT(x[CORNER_OF_SIDE(theElement,i,1)],
                x[CORNER_OF_SIDE(theElement,i,0)],a);
    V3_SUBTRACT(x[CORNER_OF_SIDE(theElement,i,2)],
                x[CORNER_OF_SIDE(theElement,i,0)],b);
    V3_VECTOR_PRODUCT(a,b,rot);
    V3_SUBTRACT(global,x[CORNER_OF_SIDE(theElement,i,0)],b);
    V3_SCALAR_PRODUCT(rot,b,det);
    if (det > SMALL_C)
      return(0);
  }

  return(1);
}
#endif

/****************************************************************************/
/*D
   FindElementFromPosition - Find element containing position

   SYNOPSIS:
   ELEMENT *FindElementFromPosition (GRID *theGrid, DOUBLE *pos)

   PARAMETERS:
   .  theGrid - grid level to search
   .  pos - given position

   DESCRIPTION:
   This function finds the first element containing the position `pos`.

   RETURN VALUE:
   ELEMENT *
   .n   pointer to ELEMENT
   .n   NULL if not found.
   D*/
/****************************************************************************/

ELEMENT *FindElementFromPosition (GRID *theGrid, DOUBLE *pos)
{
  ELEMENT *theElement;

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
    if (PointInElement(pos,theElement) == 1)
      return(theElement);

  return(NULL);
}

ELEMENT *FindElementOnSurface (MULTIGRID *theMG, DOUBLE *global)
{
  ELEMENT *t;
  INT k;

  for (k=0; k<=TOPLEVEL(theMG); k++)
    for (t=FIRSTELEMENT(GRID_ON_LEVEL(theMG,k)); t!=NULL; t=SUCCE(t))
      if (EstimateHere(t))
        if (PointInElement(global,t)) return(t);

  return(NULL);
}

/****************************************************************************/
/*D
   NeighbourElement - get the neighbouring element

   SYNOPSIS:
   ELEMENT *NeighbourElement (ELEMENT *t, INT side);

   PARAMETERS:
   .  theElement - pointer to an element
   .  side - number of an element side

   DESCRIPTION:
   This function returns a pointer to the element on the given side.

   RETURN VALUE:
   ELEMENT *
   .n    pointer to an element
   .n    NULL if error occured.
   D*/
/****************************************************************************/

INT InnerBoundary (ELEMENT *t, INT side)
{
  INT left,right;

  ASSERT(OBJT(t) == BEOBJ);
  ASSERT(SIDE_ON_BND(t,side));

  BNDS_BndSDesc(ELEM_BNDS(t,side),&left,&right);

  return((left != 0) && (right != 0));
}

ELEMENT *NeighbourElement (ELEMENT *t, INT side)
{
  ELEMENT *e, *nb;

  nb = NBELEM(t,side);

  if (nb==NULL)
  {
    if (OBJT(t)==BEOBJ)
      if (SIDE_ON_BND(t,side))
        if (!INNER_BOUNDARY(t,side))
          return(NULL);

    /* It may happen, that a neighbor is not present on that level. */
    /* In that case go to father until a neighbor is found.			*/
    for (e=t; e!=NULL; e=EFATHER(e))
    {
      /* t must be a copy of e in order to have a correct side */
      if (NSONS(e)>1) return(NULL);

      /* now let's see if we have a neighbour */
      nb = NBELEM(e,side);
      if (nb != NULL) break;
    }
  }
  else if (NSONS(nb) == 1)
  {
    nb = SON(nb,0);
    if (NSONS(nb) == 1)
      nb = SON(nb,0);
  }

  return(nb);
}

/****************************************************************************/
/*D
   ListMultiGrid - List general information about multigrid structure

   SYNOPSIS:
   void ListMultiGrid (MULTIGRID *theMG, const INT isCurrent, const INT longformat);

   PARAMETERS:
   .  theMG - structure to list
   .  isCurrent - is `theMG` current multigrid
   .  longformat - print all information or only name of `theMG`

   DESCRIPTION:
   This function lists general information about a multigrid structure.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListMultiGridHeader (const INT longformat)
{
  if (longformat)
    sprintf(buffer,"   %-20.20s %-20.20s %-20.20s %10.10s %10.10s\n","mg name","domain name","problem name","heap size","heap used");
  else
    sprintf(buffer,"   %-20.20s\n","mg name");
}

void ListMultiGrid (MULTIGRID *theMG, const INT isCurrent, const INT longformat)
{
  char c;
  BVP *theBVP;
  BVP_DESC theBVPDesc;

  /* get BVP description */
  theBVP = MG_BVP(theMG);
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc))
  {
    PrintErrorMessage('E',"ListMultiGrid","cannot evaluate BVP");
    return;
  }

  c = isCurrent ? '*' : ' ';

  if (longformat)
    sprintf(buffer," %c %-20.20s %-20.20s %10lu %10lu\n",c,ENVITEM_NAME(theMG),
            BVPD_NAME(theBVPDesc), HeapSize(theMG->theHeap),HeapUsed(theMG->theHeap));
  else
    sprintf(buffer," %c %-20.20s\n",c,ENVITEM_NAME(theMG));

  UserWrite(buffer);

  return;
}


/****************************************************************************/
/*D
   ListGrids - list general information about grids of multigrid

   SYNOPSIS:
   void ListGrids (const MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid structure

   DESCRIPTION:
   This function lists general information about the grids of a multigrid.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListGrids (const MULTIGRID *theMG)
{
  GRID *theGrid;
  ELEMENT *theElement,*NBElem;
  VERTEX *myVertex,*nbVertex,*v0,*v1;
  NODE *theNode,*n0,*n1;
  EDGE *theEdge;
  LINK *theLink;
  VECTOR *vec;
  MATRIX *mat;
  char c;
  DOUBLE hmin,hmax,h;
  INT l,cl,minl,i,soe,eos,coe,side,e;
  INT nn,ne,nt,ns,nvec,nc,free,used,heap;

  cl = CURRENTLEVEL(theMG);

  UserWriteF("grids of '%s':\n",ENVITEM_NAME(theMG));

  UserWrite("level maxlevel    #vert    #node    #edge    #elem    #side    #vect    #conn  minedge  maxedge\n");
  for (l=0; l<=TOPLEVEL(theMG); l++)
  {
    theGrid = GRID_ON_LEVEL(theMG,l);

    c = (l==cl) ? '*' : ' ';

    /* calculate minimal and maximal edge */
    hmin = MAX_C;
    hmax = 0.0;
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      myVertex = MYVERTEX(theNode);
      for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
      {
        nbVertex = MYVERTEX(NBNODE(theLink));
        V_DIM_EUKLIDNORM_OF_DIFF(CVECT(myVertex),CVECT(nbVertex),h);
        hmin = MIN(hmin,h);
        hmax = MAX(hmax,h);
      }
    }
    ns = 0;
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement))
      if (OBJT(theElement) == BEOBJ)
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
          if (SIDE_ON_BND(theElement,i))
            ns++;

    UserWriteF("%c %3d %8d %8ld %8ld %8ld %8ld %8ld %8ld %8ld %9.3e %9.3e\n",c,l,(int)TOPLEVEL(theMG),
               (long)NV(theGrid),(long)NN(theGrid),(long)NE(theGrid),(long)NT(theGrid),
               (long)ns,(long)NVEC(theGrid),(long)NC(theGrid),(float)hmin,(float)hmax);
  }

  /* surface grid up to current level */
  minl = cl;
  hmin = MAX_C;
  hmax = 0.0;
  nn = ne = nt = ns = nvec = nc = 0;
  for (l=0; l<=cl; l++)
  {
    theGrid = GRID_ON_LEVEL(theMG,l);

    /* reset USED flags in all objects to be counted */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      SETUSED(theNode,0);
      for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
        SETUSED(MYEDGE(theLink),0);
    }
    for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
      for (mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
        SETCUSED(MMYCON(mat),0);

    /* count vectors and connections */
    for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
      if ((l==cl) || (VNCLASS(vec)<1))
      {
        nvec++;
        for (mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
        {
          if (MUSED(mat)) continue;
          SETCUSED(MMYCON(mat),1);

          if ((l==cl) || (VNCLASS(MDEST(mat))<1))
            nc++;
        }
      }

    /* count other objects */
    for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((NSONS(theElement)==0) || (l==cl))
      {
        nt++;

        minl = MIN(minl,l);

        coe = CORNERS_OF_ELEM(theElement);
        for (i=0; i<coe; i++)
        {
          theNode = CORNER(theElement,i);
          if (USED(theNode)) continue;
          SETUSED(theNode,1);

          if ((SONNODE(theNode)==NULL) || (l==cl))
            nn++;
        }

        soe = SIDES_OF_ELEM(theElement);
        for (side=0; side<soe; side++)
        {
          if (OBJT(theElement)==BEOBJ)
            if (ELEM_BNDS(theElement,side)!=NULL) ns++;

          /* check neighbour element */
          if (l<cl)
            if ((NBElem=NBELEM(theElement,side))!=NULL)
              if (NSONS(NBElem)>0)
                continue;                                                       /* objects of this side will be counted by the neighbour */

          eos = EDGES_OF_SIDE(theElement,side);
          for (i=0; i<eos; i++)
          {
            e  = EDGE_OF_SIDE(theElement,side,i);
            n0 = CORNER(theElement,CORNER_OF_EDGE(theElement,e,0));
            v0 = MYVERTEX(n0);
            n1 = CORNER(theElement,CORNER_OF_EDGE(theElement,e,1));
            v1 = MYVERTEX(n1);

            if ((theEdge=GetEdge(n0,n1))==NULL) continue;

            if (USED(theEdge)) continue;
            SETUSED(theEdge,1);

            ne++;

            V_DIM_EUKLIDNORM_OF_DIFF(CVECT(v0),CVECT(v1),h);
            hmin = MIN(hmin,h);
            hmax = MAX(hmax,h);
          }
        }
      }
  }

  UserWrite("\nsurface grid up to current level:\n");
  UserWriteF("%c %3d %8d %8s %8ld %8ld %8ld %8ld %8ld %8ld %9.3e %9.3e\n",' ',minl,(int)cl,
             "---",(long)nn,(long)ne,(long)nt,
             (long)ns,(long)nvec,(long)nc,(float)hmin,(float)hmax);

  /* storage */
  heap = HeapFreelistUsed(MGHEAP(theMG));
  used = HeapUsed(MGHEAP(theMG)) - heap;
  free = HeapSize(MGHEAP(theMG)) - used;
#ifdef MEM_SIZE_ULL
  if (heap == 0)
    UserWriteF("\n%llu bytes used out of %d allocated\n",used,used+free);
  else
    UserWriteF("\n%llu ( %llu + %llu ) bytes used out of %llu allocated\n",used+heap,used,heap,used+free);
#else
  if (heap == 0)
    UserWriteF("\n%lu bytes used out of %d allocated\n",used,used+free);
  else
    UserWriteF("\n%lu ( %lu + %lu ) bytes used out of %lu allocated\n",used+heap,used,heap,used+free);
#endif
}

/****************************************************************************/
/*D
   ListNode - List information about node in multigrid

   SYNOPSIS:
   void ListNode (MULTIGRID *theMG, NODE *theNode, INT dataopt,
   INT bopt, INT nbopt, INT vopt);

   PARAMETERS:
   .  theMG - structure containing the node
   .  theNode - node to list
   .  dataopt - list user data if true
   .  bopt - list boundary info if true
   .  nbopt - list info about neighbors if true
   .  vopt - list more information

   DESCRIPTION:
   This function lists information about a node in a multigrid.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListNode (MULTIGRID *theMG, NODE *theNode, INT dataopt, INT bopt, INT nbopt, INT vopt)
{
  VERTEX *theVertex;
  LINK *theLink;
  INT i;

  theVertex = MYVERTEX(theNode);

  /******************************/
  /* print standard information */
  /******************************/
  /* line 1 */ sprintf(buffer,"NODEID=" ID_FFMTE " CTRL=%8lx IX=%8ld VEID=" VID_FFMTX " LEVEL=%2d",
                       ID_PRTE(theNode),(long)CTRL(theNode),
                       (long)INDEX(theNode),VID_PRTX(theVertex),LEVEL(theNode));
  UserWrite(buffer);

  /* print coordinates of that node */
  for(i=0; i<DIM; i++)
  {
    sprintf(buffer," x%1d=%11.4E",i, (float)(CVECT(theVertex)[i]) );
    UserWrite(buffer);
  }
  UserWrite("\n");

  if (vopt)       /* verbose: print all information */
  {
    /* line 2 */ sprintf(buffer,"   CLASS=%d ",CLASS(theNode));
    UserWrite(buffer);

    /* print nfather information */
    if (NFATHER(theNode)!=NULL)
    {
      sprintf(buffer," NODEFATHER=%ld",(long)ID(NFATHER(theNode)));
      UserWrite(buffer);
    }
    UserWrite("\n");

    /* line 3 */	/* print vertex father information */
    if (VFATHER(theVertex)!=NULL)
    {
      sprintf(buffer,"   VERTEXFATHER=%ld ",(long)ID(VFATHER(theVertex)));
      UserWrite(buffer);
      for(i=0; i<DIM; i++)
      {
        sprintf(buffer,"XI[%d]=%11.4E ",i, (float)(LCVECT(theVertex)[i]) );
        UserWrite(buffer);
      }
    }
    UserWrite("\n");
  }

  /******************************/
  /* print boundary information */
  /******************************/
  if (bopt)
  {
    if (OBJT(theVertex) == BVOBJ)
    {
      if (BNDP_BndPDesc(V_BNDP(theVertex),&i))
        UserWrite("Error in boundary point\n");
      else
        UserWriteF("boundary point: move %d moved %d\n",i,
                   MOVED(theVertex));
    }
  }

  if (nbopt)
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
                        #if defined __THREEDIM__ && defined ModelP
      sprintf(buffer,"   EDGE=%x/%08x ",MYEDGE(theLink),
              DDD_InfoGlobalId(PARHDR(MYEDGE(theLink))));
                        #else
      sprintf(buffer,"   ");
                        #endif
      UserWrite(buffer);
      sprintf(buffer,"NB=" ID_FMTX " CTRL=%8lx NO_OF_ELEM=%3d ",
              ID_PRTX(NBNODE(theLink)),(long)CTRL(theLink),
              NO_OF_ELEM(MYEDGE(theLink)));
      UserWrite(buffer);

      theVertex=MYVERTEX(NBNODE(theLink));
      for(i=0; i<DIM; i++)
      {
        sprintf(buffer,"x%1d=%11.4 ",i, (float)(CVECT(theVertex)[i]) );
        UserWrite(buffer);
      }
      UserWrite("\n");
    }
  }
  return;
}


/****************************************************************************/
/*D
   ListNodeSelection - List information about all nodes in selection

   SYNOPSIS:
   void ListNodeSelection (MULTIGRID *theMG, INT dataopt, INT bopt, INT nbopt, INT vopt);

   PARAMETERS:
   .  theMG - structure containing the nodes
   .  dataopt - list user data if true
   .  bopt - list boundary info if true
   .  nbopt - list info about neighbors if true
   .  vopt - list more information

   DESCRIPTION:
   This function lists information about all nodes in the selection.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListNodeSelection (MULTIGRID *theMG, INT dataopt, INT bopt, INT nbopt, INT vopt)
{
  int j;
  NODE *theNode;

  if (SELECTIONMODE(theMG) != nodeSelection)
  {
    PrintErrorMessage('E',"ListNodeSelection","wrong selection type");
    return;
  }
  for(j=0; j<SELECTIONSIZE(theMG); j++)
  {
    theNode = (NODE *) SELECTIONOBJECT(theMG,j);
    ListNode(theMG,theNode,dataopt,bopt,nbopt,vopt);
  }
}


/****************************************************************************/
/*D
   IsNodeSelected - Check if element is in selection list

   SYNOPSIS:
   INT IsNodeSelected (MULTIGRID *theMG, NODE *theNode);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theNode - node to check

   DESCRIPTION:
   This function checks if an element is in the selection list.

   RETURN VALUE:
   INT
   .n   0 if NOT in list
   .n   1 if in list.
   D*/
/****************************************************************************/

INT IsNodeSelected (MULTIGRID *theMG, NODE *theNode)
{
  int j;

  if (SELECTIONMODE(theMG) != nodeSelection) return (0);
  for(j=0; j<SELECTIONSIZE(theMG); j++)
    if (theNode == (NODE *) SELECTIONOBJECT(theMG,j))
      return (1);
  return (0);
}


/****************************************************************************/
/*D
   ListNodeRange - List information about nodes in given range of ids

   SYNOPSIS:
   void ListNodeRange (MULTIGRID *theMG, INT from, INT to, INT dataopt,
   INT bopt, INT nbopt, INT vopt)

   PARAMETERS:
   .  theMG - structure to list
   .  from - first id
   .  to - last id
   .  dataopt - list user data if true
   .  bopt - list boundary info if true
   .  nbopt - list info about neighbors if true
   .  vopt - list more information

   DESCRIPTION:
   This function list information about all nodes in a given range of ids.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListNodeRange (MULTIGRID *theMG, INT from, INT to, INT dataopt, INT bopt, INT nbopt, INT vopt)
{
  int level;
  NODE *theNode;

  for (level=0; level<=TOPLEVEL(theMG); level++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,level)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      if ( (ID(theNode)>=from)&&(ID(theNode)<=to) )
        ListNode(theMG,theNode,dataopt,bopt,nbopt,vopt);
    }
}


/****************************************************************************/
/*D
   ListElement - List information about element

   SYNOPSIS:
   void ListElement (MULTIGRID *theMG, ELEMENT *theElement, INT dataopt,
   INT bopt, INT nbopt, INT vopt)

   PARAMETERS:
   .  theMG -  structure to list
   .  theElement - element to list
   .  dataopt - list user data if true
   .  vopt - list more information

   DESCRIPTION:
   This function lists information about an element

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListElement (MULTIGRID *theMG, ELEMENT *theElement, INT dataopt, INT bopt, INT nbopt, INT vopt)
{
  char etype[10];
  char ekind[8];
  int i,j;
        #ifdef __THREEDIM__
  ELEMENT *SonList[MAX_SONS];
        #endif

  if (DIM==2)
    switch (TAG(theElement))
    {
    case TRIANGLE :                  strcpy(etype,"TRI"); break;
    case QUADRILATERAL :             strcpy(etype,"QUA"); break;
    default :                strcpy(etype,"???"); break;
    }
  else
    switch (TAG(theElement))
    {
    case TETRAHEDRON :               strcpy(etype,"TET"); break;
    case PYRAMID :                   strcpy(etype,"PYR"); break;
    case PRISM :                             strcpy(etype,"PRI"); break;
    case HEXAHEDRON :                strcpy(etype,"HEX"); break;
    default :                strcpy(etype,"???"); break;
    }
  switch (ECLASS(theElement))
  {
  case YELLOW_CLASS :              strcpy(ekind,"YELLOW "); break;
  case GREEN_CLASS :               strcpy(ekind,"GREEN  "); break;
  case RED_CLASS :                 strcpy(ekind,"RED    "); break;
  default :                strcpy(ekind,"???    "); break;
  }
  sprintf(buffer,"ELEMID=" EID_FFMTE " %5s %5s CTRL=%8lx CTRL2=%8lx REFINE=%2d MARK=%2d LEVEL=%2d",
          EID_PRTE(theElement),ekind,etype,
          (long)CTRL(theElement),(long)FLAG(theElement),REFINE(theElement),MARK(theElement),LEVEL(theElement));
  UserWrite(buffer);
  if (COARSEN(theElement)) UserWrite(" COARSEN");
  UserWrite("\n");

  if (vopt)
  {
    UserWrite("   ");
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      sprintf(buffer," N%d=" ID_FMTX,i,ID_PRTX(CORNER(theElement,i)));
      UserWrite(buffer);
    }
    if (EFATHER(theElement)!=NULL)
    {
      sprintf(buffer," FA=%ld ",(long)ID(EFATHER(theElement)));
      UserWrite(buffer);
    }
    UserWrite("\n");
    UserWrite("   ");
                #ifdef __TWODIM__
    for (i=0; i<SONS_OF_ELEM(theElement); i++)
      if (SON(theElement,i)!=NULL)
      {
        sprintf(buffer,"S%d=%ld ",i,(long)ID(SON(theElement,i)));
        UserWrite(buffer);
      }
                #endif
                #ifdef __THREEDIM__
    if (GetSons(theElement,SonList)!=0) return;
    for (i=0; i<NSONS(theElement); i++)
    {
      sprintf(buffer,"S%d=%ld ",i,(long)ID(SonList[i]));
      UserWrite(buffer);
    }
                #endif
    UserWrite("\n");
  }
  if (nbopt)
  {
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if (NBELEM(theElement,i)!=NULL)
      {
        sprintf(buffer,"NB%d=%ld ",i,(long)ID(NBELEM(theElement,i)));
        UserWrite(buffer);
      }
    UserWrite("\n");
  }
  if (bopt)
  {
    UserWrite("   ");
    if (OBJT(theElement)==BEOBJ)
    {
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        for(j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
        {
                                                #ifdef __THREEDIM__
                        #ifdef ModelP
          sprintf(buffer,"  NODE[ID=%ld]: ",(long)(ID(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j)))));
          UserWrite(buffer);
                                                #endif
                                                #endif
          UserWrite("\n");
        }
      }
    }
    UserWrite("\n");
  }

  return;
}


/****************************************************************************/
/*D
    ListElementSelection - list information about elements in selection

   SYNOPSIS:
   void ListElementSelection (MULTIGRID *theMG, INT dataopt,
   INT bopt, INT nbopt, INT vopt);

   PARAMETERS:
   .  theMG: multigrid structure to list
   .  dataopt - list user data if true
   .  vopt - list more information

   DESCRIPTION:
   This function lists information about all elements in the selection.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListElementSelection (MULTIGRID *theMG, INT dataopt, INT bopt, INT nbopt, INT vopt)
{
  int j;
  ELEMENT *theElement;

  if (SELECTIONSIZE(theMG) <= 0) return;
  if (SELECTIONMODE(theMG) != elementSelection)
  {
    PrintErrorMessage('E',"ListElementSelection","wrong selection type");
    return;
  }
  for(j=0; j<SELECTIONSIZE(theMG); j++)
  {
    theElement = (ELEMENT *) SELECTIONOBJECT(theMG,j);
    ListElement(theMG,theElement,dataopt,bopt,nbopt,vopt);
  }
}


/****************************************************************************/
/*D
   IsElementSelected - Check whether element is in selection list

   SYNOPSIS:
   INT IsElementSelected (MULTIGRID *theMG, ELEMENT *theElement);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theElement - element to check

   DESCRIPTION:
   This function checks whether an element is in the selection list.

   RETURN VALUE:
   INT
   .n   0 if NOT in list
   .n   1 if in list.
   D*/
/****************************************************************************/

INT IsElementSelected (MULTIGRID *theMG, ELEMENT *theElement)
{
  int j;

  if (SELECTIONMODE(theMG) != elementSelection) return (0);
  for(j=0; j<SELECTIONSIZE(theMG); j++)
    if (theElement == (ELEMENT *) SELECTIONOBJECT(theMG,j))
      return (1);
  return (0);
}


/****************************************************************************/
/*D
   ListElementRange - List information about elements in range of ids

   SYNOPSIS:
   void ListElementRange (MULTIGRID *theMG, INT from, INT to, INT dataopt,
   INT bopt, INT nbopt, INT vopt, INT lopt);

   PARAMETERS:
   .  theMG - multigrid structure to list
   .  from - first id
   .  to - last id
   .  dataopt - list user data if true
   .  vopt - list more information

   DESCRIPTION:
   This function lists information about all elements in a range of ids.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListElementRange (MULTIGRID *theMG, INT from, INT to, INT dataopt, INT bopt, INT nbopt, INT vopt, INT lopt)
{
  int level;
  ELEMENT *theElement;

  if (lopt==FALSE)
    for (level=0; level<=TOPLEVEL(theMG); level++)
      for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,level)); theElement!=NULL; theElement=SUCCE(theElement))
      {
        if ( (ID(theElement)>=from)&&(ID(theElement)<=to) )
          ListElement(theMG,theElement,dataopt,bopt,nbopt,vopt);
      }
  else
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG))); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if ( (ID(theElement)>=from)&&(ID(theElement)<=to) )
        ListElement(theMG,theElement,dataopt,bopt,nbopt,vopt);
    }
}


/****************************************************************************/
/*D
   ListVector - List information about vector

   SYNOPSIS:
   void ListVector (MULTIGRID *theMG, VECTOR *theVector, INT matrixopt, INT dataopt);

   PARAMETERS:
   .  theMG - multigrid structure to list
   .  theVector - vector to list
   .  matrixopt - list line of matrix corresponding to theVector
   .  dataopt - list user data if true

   DESCRIPTION:
   This function lists information about a vector.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListVector (MULTIGRID *theMG, VECTOR *theVector, INT matrixopt, INT dataopt)
{
  FORMAT *theFormat;
  NODE *theNode;
  EDGE *theEdge;
  ELEMENT *theElement;
  MATRIX *theMatrix;
  void *Data;

  theFormat = MGFORMAT(theMG);

  /* print type of vector */
  if (VTYPE(theVector)==NODEVECTOR)
  {
    theNode = (NODE*)VOBJECT(theVector);
    sprintf(buffer,"NODE-V IND=" VINDEX_FFMTE " nodeID=" ID_FMTX
            "                VCLASS=%1d VNCLASS=%1d\n",
            VINDEX_PRTE(theVector),ID_PRTX(theNode),VCLASS(theVector),VNCLASS(theVector));
    UserWrite(buffer);
  }
  if (VTYPE(theVector)==EDGEVECTOR)
  {
    theEdge = (EDGE*)VOBJECT(theVector);
    sprintf(buffer,"EDGE-V IND=" VINDEX_FFMTE " fromID=" ID_FFMT
            " to__ID=%7ld VCLASS=%1d VNCLASS=%1d\n",
            VINDEX_PRTE(theVector),ID_PRT(NBNODE(LINK0(theEdge))),
            ID(NBNODE(LINK1(theEdge))),VCLASS(theVector),VNCLASS(theVector));
    UserWrite(buffer);
  }
        #ifdef __THREEDIM__
  if (VTYPE(theVector)==SIDEVECTOR)
  {
    theElement = (ELEMENT*)VOBJECT(theVector);
    sprintf(buffer,"SIDE-V IND=" VINDEX_FFMTE " elemID=" EID_FFMT
            "                VCLASS=%1d VNCLASS=%1d\n",
            VINDEX_PRTE(theVector),EID_PRT(theElement),
            VCLASS(theVector),VNCLASS(theVector));
    UserWrite(buffer);
  }
        #endif
  if (VTYPE(theVector)==ELEMVECTOR)
  {
    theElement = (ELEMENT*)VOBJECT(theVector);
    sprintf(buffer,"ELEM-V IND=" VINDEX_FFMTE " elemID=" EID_FFMT
            "                VCLASS=%1d VNCLASS=%1d\n",
            VINDEX_PRTE(theVector),EID_PRT(theElement),
            VCLASS(theVector),VNCLASS(theVector));
    UserWrite(buffer);
  }


  /* print vector data if */
  if (dataopt && theFormat->PrintVector[VTYPE(theVector)]!=NULL)
  {
    Data = (void*)(&VVALUE(theVector,0));
    if ((*(theFormat->PrintVector[VTYPE(theVector)]))(Data,"   ",buffer))
      return;
    UserWrite(buffer);
    UserWrite("\n");
  }

  /* print matrix list if */
  if (matrixopt)
    for (theMatrix = VSTART(theVector); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      UserWrite("    DEST(MATRIX): ");
      ListVector(theMG,MDEST(theMatrix),0,0);

      /* print matrix data if */
      if (dataopt && theFormat->PrintMatrix[MROOTTYPE(theMatrix)][MDESTTYPE(theMatrix)]!=NULL)
      {
        Data = (void*)(&MVALUE(theMatrix,0));
        if ((*(theFormat->PrintMatrix[MROOTTYPE(theMatrix)][MDESTTYPE(theMatrix)]))(Data,"       ",buffer))
          return;
        UserWrite(buffer);
        UserWrite("\n");
      }
    }
  return;
}

/****************************************************************************/
/*D
   ListVectorOfElementSelection - List info about vectors of elements in selection

   SYNOPSIS:
   void ListVectorOfElementSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt);

   PARAMETERS:
   .  theMG -  structure to list
   .  matrixopt - list line of matrix corresponding to theVector
   .  dataopt - list user data if true

   DESCRIPTION:
   This function lists info about all vectors of elements in the selection.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListVectorOfElementSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt)
{
  int i,j;
  ELEMENT *theElement;
  VECTOR *vList[20];
  INT cnt;

  if (SELECTIONMODE(theMG) != elementSelection)
  {
    PrintErrorMessage('E',"ListVectorOfElementSelection","wrong selection type");
    return;
  }
  for (j=0; j<SELECTIONSIZE(theMG); j++)
  {
    theElement = (ELEMENT *) SELECTIONOBJECT(theMG,j);
    sprintf(buffer,"ELEM(ID=%d):\n",ID(theElement));
    UserWrite(buffer);

    if (TYPE_DEF_IN_MG(theMG,NODEVECTOR))
    {
      GetVectorsOfNodes(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    }
    if (TYPE_DEF_IN_MG(theMG,EDGEVECTOR))
    {
      GetVectorsOfEdges(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    }
                #ifdef __THREEDIM__
    if (TYPE_DEF_IN_MG(theMG,SIDEVECTOR))
    {
      GetVectorsOfSides(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    }
                #endif
    if (TYPE_DEF_IN_MG(theMG,ELEMVECTOR))
    {
      GetVectorsOfElement(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    }
  }
}

/****************************************************************************/
/*D
   ListVectorSelection - list information about vectors in selection

   SYNOPSIS:
   void ListVectorSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt)

   PARAMETERS:
   .  theMG: multigrid structure to list
   .  matrixopt - list matrices of this vector
   .  dataopt - list user data if true

   DESCRIPTION:
   This function lists information about all elements in the selection.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListVectorSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt)
{
  int j;
  VECTOR *theVector;

  if (SELECTIONSIZE(theMG) <= 0) return;
  if (SELECTIONMODE(theMG) != vectorSelection)
  {
    PrintErrorMessage('E',"ListVectorSelection","wrong selection type");
    return;
  }
  for(j=0; j<SELECTIONSIZE(theMG); j++)
  {
    theVector = (VECTOR *) SELECTIONOBJECT(theMG,j);
    ListVector(theMG,theVector,matrixopt,dataopt);
  }
}

/****************************************************************************/
/*D
   IsVectorSelected - Check whether vector is in selection list

   SYNOPSIS:
   INT IsVectorSelected (MULTIGRID *theMG, VECTOR *theVector);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theVector - vector to check

   DESCRIPTION:
   This function checks whether an element is in the selection list.

   RETURN VALUE:
   INT
   .n   0 if NOT in list
   .n   1 if in list.
   D*/
/****************************************************************************/

INT IsVectorSelected (MULTIGRID *theMG, VECTOR *theVector)
{
  int j;

  if (SELECTIONMODE(theMG) != vectorSelection) return (0);
  for(j=0; j<SELECTIONSIZE(theMG); j++)
    if (theVector == (VECTOR *) SELECTIONOBJECT(theMG,j))
      return (1);
  return (0);
}

/****************************************************************************/
/*D
   ListVectorRange - list information about vectors in range of ids

   SYNOPSIS:
   void ListVectorRange (MULTIGRID *theMG, INT fl, INT tl,
   INT from, INT to, INT matrixopt, INT dataopt)

   PARAMETERS:
   .  theMG - structure to list
   .  from - first index
   .  to - last index
   .  matrixopt - list line of matrix corresponding to theVector
   .  dataopt - list user data if true

   DESCRIPTION:
   This function lists information about all vectors in a given range of indices.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListVectorRange (MULTIGRID *theMG, INT fl, INT tl, INT from, INT to, INT matrixopt, INT dataopt)
{
  int level;
  VECTOR *theVector;

  for (level=fl; level<=tl; level++)
    for (theVector=PFIRSTVECTOR(GRID_ON_LEVEL(theMG,level)); theVector!=NULL; theVector=SUCCVC(theVector))
    {
      if (VINDEX(theVector)>=from && VINDEX(theVector)<=to)
        ListVector(theMG,theVector,matrixopt,dataopt);
    }
}

static INT CheckElementold (ELEMENT *theElement, INT *SideError, INT *EdgeError, INT *NodeError)
{
  int i,j,k,n;
  EDGE *theEdge;
  ELEMENT *NbElement;
  VERTEX *theVertex;

  *SideError = 0;
  *NodeError = 0;
  n = SIDES_OF_ELEM(theElement);
  if (ECLASS(theElement)!=YELLOW_CLASS)
    for (i=0; i<n; i++)
    {
      NbElement = NBELEM(theElement,i);
      if (NbElement != NULL)
      {
        /* lets see if NbElement has the neighbor theElement */
        for (j=0; j<SIDES_OF_ELEM(NbElement); j++)
          if (NBELEM(NbElement,j) == theElement)
            break;
        if (j == SIDES_OF_ELEM(NbElement))
          *SideError |= (1<<i);
      }
      else
      {
        if (OBJT(theElement) == IEOBJ)
          *SideError |= (1<<(i+n));
        else
        {
          if (SIDE_ON_BND(theElement,i))
          {
            for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
            {
              theVertex = MYVERTEX(CORNER(theElement,(k=(i+j)%n)));
              if (OBJT(theVertex) == IVOBJ)
                *NodeError |= (1<<(k+CORNERS_OF_ELEM(theElement)));
            }
          }
          else
            *SideError |= (1<<(i+2*n));
        }
      }
    }

  *EdgeError = 0;
  for (i=0; i<n; i++)
  {
    SETUSED(CORNER(theElement,i),1);
    for (j=i+1; j<n; j++)
    {
      /* there are no diagonal edges in quadrilaterals */
      if ((n==4) && (j==i+2))
        continue;
      theEdge = GetEdge(CORNER(theElement,i),CORNER(theElement,j));
      if (theEdge==NULL)
        *EdgeError |= 1<<i;
      else
        SETUSED(theEdge,1);
    }
  }

  if (*SideError || *EdgeError || *NodeError)
    return (1);
  return (0);
}

static INT CheckGeometryold (GRID *theGrid) /* 2D VERSION */
{
  NODE *theNode;
  ELEMENT *theElement;
  LINK *theLink;
  int i,j;
  INT SideError, EdgeError, NodeError,count,n;

  /* check neighbors and edges of elements */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    n = SIDES_OF_ELEM(theElement);
    if (CheckElementold(theElement, &SideError, &EdgeError, &NodeError)==0) continue;
    sprintf(buffer,"ELEM " EID_FMTE ":\n", EID_PRTE(theElement));
    UserWrite(buffer);

    if (SideError)
      for (i=0; i<n; i++)
      {
        if (SideError & 1<<i)
        {
          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            sprintf(buffer,ID_FMT, ID_PRT(CORNER(theElement,(i+j)%n)));
            UserWrite(buffer);
            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has neighbour but a backPtr does not exist\n");
        }
        if (SideError & 1<<(i+n))
        {
          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            sprintf(buffer, ID_FMT, ID_PRT(CORNER(theElement,(i+j)%n)));
            UserWrite(buffer);
            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has no neighbour but element is IEOBJ\n");
        }
        if (SideError & 1<<(i+2*n))
        {
          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            sprintf(buffer, ID_FMT, ID_PRT(CORNER(theElement,(i+j)%n)));
            UserWrite(buffer);
            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has no neighbour, element is BEOBJ, "
                    "but there is no SIDE\n");
        }
      }
    if (EdgeError)
      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {
        if (!(EdgeError & 1<<i)) continue;
        sprintf(buffer,"   EDGE(" ID_FMT "," ID_FMT ") is missing\n",
                ID_PRT(CORNER(theElement,i)), ID_PRT(CORNER(theElement,(i+1)%n)));
        UserWrite(buffer);
      }
    if (NodeError)
      for (i=0; i<n; i++)
      {
        if (NodeError & (1<<i))
        {
          sprintf(buffer,"   CORNER " ID_FMT " is BVOBJ, "
                  "ids from elementside and vertexsegment are "
                  "not consistent\n", ID_PRT(CORNER(theElement,i)));
          UserWrite(buffer);
        }
        if (NodeError & (1<<(i+n)))
        {
          sprintf(buffer,"   CORNER " ID_FMT " is IVOBJ, but lies on elementside\n",
                  ID_PRT(CORNER(theElement,i)));
          UserWrite(buffer);
        }
      }
  }

  /* look for dead edges */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      SETUSED(theLink,1);
    }
  }
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      if (USED(theLink)!=1 || USED(REVERSE(theLink))!=1)
      {
                #ifdef ModelP
        UserWriteF("edge between " ID_FMT " and " ID_FMT
                   " dead: USED=%d USEDREV=%d \n",
                   ID_PRT(theNode), ID_PRT(NBNODE(theLink)),
                   USED(theLink),USED(REVERSE(theLink)));
                                #endif
      }
    }
  }
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      SETUSED(theLink,0);
    }
  }

  /* look for dead nodes */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    if (!USED(theNode))
    {
            #ifdef ModelP
      sprintf(buffer,"node " ID_FMTX " is dead ",ID_PRTX(theNode));
      UserWrite(buffer);
      UserWrite("\n");
                        #endif
    }
    else
      SETUSED(theNode,0);
  }

  /* check number of elem and their pointers */
  count = 0;
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (SUCCE(theElement)!=NULL)
    {
      if (OBJT(SUCCE(theElement))!=IEOBJ && OBJT(SUCCE(theElement))!=BEOBJ)
      {
        sprintf(buffer,"pointer of ELEM(" EID_FMT ") (number %ld) to next element is no pointer to an element\n", EID_PRT(theElement),(long)count);
        UserWrite(buffer);
        break;
      }
      if (PREDE(SUCCE(theElement))!=NULL)
      {
        if (PREDE(SUCCE(theElement))!=theElement)
        {
          sprintf(buffer,"pointer of ELEM(" EID_FMT ") (number %ld) to previous element is not the previous element\n",EID_PRT(SUCCE(theElement)),(long)(count+1));
          UserWrite(buffer);
        }
      }
      else
      {
        sprintf(buffer,"pointer of ELEM(" EID_FMT ") (number %ld) to previous element is NULL\n",EID_PRT(SUCCE(theElement)),(long)(count+1));
        UserWrite(buffer);
      }
    }
    count++;
  }
  if (FIRSTELEMENT(theGrid)!=NULL)
    if (PREDE(FIRSTELEMENT(theGrid)) != NULL)
    {
      sprintf(buffer,"first element of the grid has a previous 'element'\n");
      UserWrite(buffer);
    }
  if (LASTELEMENT(theGrid)!=NULL)
    if (SUCCE(LASTELEMENT(theGrid)) != NULL)
    {
      sprintf(buffer,"last element of the grid has a following 'element'\n");
      UserWrite(buffer);
    }
  if (count != theGrid->nElem)
  {
    sprintf(buffer,"there are %ld elements but %ld expected\n",(long)(count),(long)theGrid->nElem);
    UserWrite(buffer);
  }
  return(GM_OK);
}

static INT CheckElement (ELEMENT *theElement, INT *SideError, INT *EdgeError,
                         INT *NodeError, INT *ESonError, INT *NSonError)
{
  INT i,j,k,l,n;
  NODE *theNode;
  EDGE *theEdge;
  ELEMENT *NbElement;
  ELEMENT *SonList[MAX_SONS];
  VERTEX *theVertex;

  *SideError = 0;
  *NodeError = 0;
  *EdgeError = 0;
  *ESonError = 0;
  *NSonError = 0;

  /* check side information */
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    NbElement = NBELEM(theElement,i);
    if (NbElement != NULL)
    {
      /* lets see if NbElement has the neighbor theElement */
      for (j=0; j<SIDES_OF_ELEM(NbElement); j++)
        if (NBELEM(NbElement,j) == theElement)
          break;
      if (j == SIDES_OF_ELEM(NbElement))
        *SideError |= (1<<i);


      if (ECLASS(theElement)!=YELLOW_CLASS)
      {
        n = CORNERS_OF_SIDE(theElement,i);
        for (k=0; k<n; k++)
          if (CORNER(theElement,CORNER_OF_SIDE(theElement,i,k))
              == CORNER(NbElement,CORNER_OF_SIDE(NbElement,j,0)))
            break;
        if (k == n)
          *SideError |= (1<<i);
        for (l=1; l<n; l++)
          if (CORNER(theElement,
                     CORNER_OF_SIDE(theElement,i,(n+k-l)%n))
              != CORNER(NbElement,CORNER_OF_SIDE(NbElement,j,l)))
          {
            *SideError |= (1<<i);
          }
      }
    }
    else
    {
      if (ECLASS(theElement)!=YELLOW_CLASS)
        if (OBJT(theElement) == IEOBJ)
                                #ifdef ModelP
          if (DDD_InfoPriority(PARHDRE(theElement)) == PrioMaster)
                                #endif
          *SideError |= (1<<(i+MAX_SIDES_OF_ELEM));

      if (OBJT(theElement) == BEOBJ)
      {
        if (SIDE_ON_BND(theElement,i))
        {
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            theVertex = MYVERTEX(CORNER(theElement,(k=CORNER_OF_SIDE(theElement,i,j))));
            if (OBJT(theVertex) == IVOBJ)
              *NodeError |= (1<<(k+MAX_CORNERS_OF_ELEM));
          }
        }
        else if (ECLASS(theElement)!=YELLOW_CLASS)
                                        #ifdef ModelP
          if (DDD_InfoPriority(PARHDRE(theElement)) == PrioMaster)
                                        #endif
          *SideError |= (1<<(i+2*MAX_SIDES_OF_ELEM));
      }
    }
  }

  /* check edge information */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    SETUSED(CORNER(theElement,i),1);
  }

  for (i=0; i<EDGES_OF_ELEM(theElement); i++)
  {
    assert(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0))!=NULL);
    assert(CORNER(theElement,CORNER_OF_EDGE(theElement,i,1))!=NULL);

    theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                      CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
    if (theEdge==NULL)
    {
      *EdgeError |= 1<<i;
      continue;
    }
    else
      SETUSED(theEdge,1);

    theNode = MIDNODE(theEdge);
    if (theNode == NULL)
      continue;

    theVertex = MYVERTEX(theNode);
    if (VFATHER(theVertex) != theElement)
      continue;

    if (i != ONEDGE(theVertex))
    {
      UserWriteF(PFMT "EID=%d VID=%d edgenumber of vertex wrong\n",
                 me,ID(theElement),ID(theVertex));
      *EdgeError |= 1<<i;
    }
  }

  /* check son information */
  if (NSONS(theElement)!=0)
  {
    if (GetSons(theElement,SonList))
    {
      UserWrite("cannot get sons\n");
      return (1);
    }
    for (i=0; i<NSONS(theElement); i++)
    {
      if (EFATHER(SonList[i])!=theElement)
        *ESonError |= (1<<i);

      if (REFINE(theElement)==0)
      {
        UserWriteF("ELEM(" EID_FMTX "): element is not refined but "
                   "has NSONS=%d\n",EID_PRTX(theElement),NSONS(theElement));
      }
    }
  }

  if (*SideError || *EdgeError || *NodeError || *ESonError || *NSonError)
    return (1);

  return (0);
}

INT CheckGeometry (GRID *theGrid)
{
  NODE *theNode;
  ELEMENT *theElement;
  EDGE *theEdge;
  LINK *theLink;
  int i,j;
  INT SideError, EdgeError, NodeError, ESonError, NSonError, count;
  INT errors = 0;

  /* reset used flags */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    SETUSED(theNode,0);
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
      SETUSED(MYEDGE(theLink),0);
  }

  /* check elements */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if (CheckElement(theElement,&SideError,&EdgeError,
                     &NodeError,&ESonError,&NSonError)==0) continue;

    UserWriteF("ELEM=" EID_FMTX "\n",EID_PRTX(theElement));

    /* evaluate side information */
    if (SideError)
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        /* back pointer failure */
        if (SideError & 1<<i)
        {
          errors++;

          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            UserWriteF(ID_FMTX,ID_PRTX(CORNER(theElement,
                                              CORNER_OF_SIDE(theElement,i,j))));

            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has neighbour but a backPtr does not exist\n");
        }

        /* neighbor pointer failure */
        if (SideError & 1<<(i+MAX_SIDES_OF_ELEM))
        {
          errors++;

          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            UserWriteF(ID_FMTX,ID_PRTX(CORNER(theElement,
                                              CORNER_OF_SIDE(theElement,i,j))));

            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has no neighbour but element is IEOBJ\n");
        }

        /* boundary failure */
        if (SideError & 1<<(i+2*MAX_SIDES_OF_ELEM))
        {
          errors++;

          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            UserWriteF(ID_FMTX,ID_PRTX(CORNER(theElement,
                                              CORNER_OF_SIDE(theElement,i,j))));

            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has no neighbour, element is BEOBJ "
                    "but there is no SIDE\n");
        }
      }

    /* evaluate edge information */
    if (EdgeError)
      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {
        if (!(EdgeError & 1<<i)) continue;

        errors++;
        UserWriteF("   EDGE(" ID_FMTX " , " ID_FMTX ") is missing\n",
                   ID_PRTX(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0))),
                   ID_PRTX(CORNER(theElement,CORNER_OF_EDGE(theElement,i,1))));
      }

    /* evaluate node information */
    if (NodeError)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        if (NodeError & (1<<i))
        {
          errors++;
          UserWriteF("   CORNER=" ID_FMTX " is BVOBJ,"
                     " ids from elementside "
                     "and vertexsegment are not consistent\n",
                     ID_PRTX(CORNER(theElement,i)));
        }
        if (NodeError & (1<<(i+MAX_CORNERS_OF_ELEM)))
        {
          errors++;
          UserWriteF("   CORNER " ID_FMTX " is IVOBJ, but lies on "
                     "elementside\n",ID_PRTX(CORNER(theElement,i)));
        }
      }

    /* evaluate son information */
    if (ESonError)
    {
      for (i=0; i<NSONS(theElement); i++)
      {
        if ((ESonError & 1<<i))
        {
          errors++;
          UserWriteF("   ESON(%d) has wrong EFATHER "
                     "pointer\n",i);
        }
      }
    }

    if (NSonError)
    {
      for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
      {
        if (NSonError & (1<<i))
        {
          errors++;
          UserWriteF("   SONNODE(CORNER %d) != CORNER(ESON)\n",i);
        }
        if (NSonError & (1<<(i+MAX_CORNERS_OF_ELEM)))
        {
          errors++;
          UserWriteF("   CORNER %d != EFATHER(CORNER(ESON))\n",i);
        }
      }

      for (i=0; i<MAX_EDGES_OF_ELEM; i++)
      {

        if (NSonError & (1<<(i+MAX_CORNERS_OF_ELEM)))
        {
          errors++;
          UserWriteF("   MIDNODE(edge %d) != CORNER(ESON)\n",i);
        }
      }

      if (NSonError & (1<<(MAX_EDGES_OF_ELEM+2*MAX_CORNERS_OF_ELEM)))
      {
        errors++;
        UserWriteF("   NFATHER(CENTERNODE(ESON)) != NULL\n");
      }
    }
  }

  /* look for dead edges */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      theEdge = MYEDGE(theLink);
      if (!USED(theEdge))
      {
        errors++;
        UserWriteF("edge between " ID_FMTX " and " ID_FMTX
                   " has no element, NO_OF_ELEM=%d \n",
                   ID_PRTX(theNode),ID_PRTX(NBNODE(theLink)),
                   NO_OF_ELEM(theEdge));

                                #ifdef Debug
        {
          NODE *nb;
          LINK *theLink1;

          nb = NBNODE(theLink);
          UserWriteF("linklist of nbnode %d:",ID(nb));

          for (theLink1=START(nb); theLink1!=NULL;
               theLink1=NEXT(theLink1))
          {
            UserWriteF(" %d-%d",ID(NBNODE(theLink1)),
                       ID(NBNODE(REVERSE(theLink1))));
          }
          UserWrite("\n");
        }
                                #endif
      }
    }
  }

  /* look for dead nodes */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    if (!USED(theNode))
    {
      errors++;
      UserWriteF("node=" ID_FMTX " is dead\n",ID_PRTX(theNode));
    }
    else
      SETUSED(theNode,0);
  }

  /* check number of elem and their pointers */
  count = 0;
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if (SUCCE(theElement)!=NULL)
    {
      if (OBJT(SUCCE(theElement))!=IEOBJ &&
          OBJT(SUCCE(theElement))!=BEOBJ)
      {
        errors++;
        UserWriteF("pointer of ELEM(" EID_FMTX ") (number %ld) "
                   "to next element is no pointer to an element\n",
                   EID_PRTX(theElement),(long)count);
        break;
      }
      if (PREDE(SUCCE(theElement))!=NULL)
      {
        if (PREDE(SUCCE(theElement))!=theElement)
        {
          errors++;
          UserWriteF("pointer of ELEM(" EID_FMTX ") (number %ld) "
                     "to previous element is not the previous element\n",
                     EID_PRTX(SUCCE(theElement)),(long)(count+1));
        }
      }
                        #ifndef ModelP
      else
      {
        errors++;
        UserWriteF("pointer of ELEM(" EID_FMTX ") (number %ld) "
                   "to previous element is NULL\n",
                   EID_PRTX(SUCCE(theElement)),(long)(count+1));
      }
                        #endif
    }
    count++;
  }

  if (FIRSTELEMENT(theGrid) != NULL)
    if (PREDE(FIRSTELEMENT(theGrid)) != NULL)
    {
      errors++;
      UserWriteF("first element of the grid has a previous 'element'\n");
    }

  if (LASTELEMENT(theGrid) != NULL)
    if (SUCCE(LASTELEMENT(theGrid)) != NULL)
    {
      errors++;
      UserWriteF("last element of the grid has a following 'element'\n");
    }

  if (count != NT(theGrid))
  {
    errors++;
    UserWriteF("there are %ld elements but %ld expected\n",(long)(count),
               (long)NT(theGrid));
  }

  return(errors);
}

INT CheckLists (GRID *theGrid)
{
  int objs = 0;

  GRID_CHECK_ELEMENT_LIST(theGrid);
  GRID_CHECK_NODE_LIST(theGrid);
  GRID_CHECK_VERTEX_LIST(theGrid);
  GRID_CHECK_VECTOR_LIST(theGrid);

  return(GM_OK);
}

/****************************************************************************/
/*D
   CheckGrid - Check consistency of data structure

   SYNOPSIS:
   INT CheckGrid (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists, INT checkif)

   PARAMETERS:
   .  theGrid - grid to check
   .  checkgeom - check geomtry
   .  checkalgebra - check algebra
   .  checklists - checklists
   .  checkif - check the processor interfaces

   DESCRIPTION:
   This function checks the consistency of data structure.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an error occured.
   D*/
/****************************************************************************/

#ifndef ModelP
INT CheckGrid (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists)
#else
INT CheckGrid (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists,
               INT checkif)
#endif
{
  INT error               = 0;
  INT errors              = 0;
  INT totalerrors = 0;

  /* check geometrical data structures */
  if (checkgeom)
  {
    UserWrite(" geometry:");
    fflush(stdout);

    if ((errors = CheckGeometry(theGrid)) != GM_OK)
    {
      totalerrors += errors;
      error++;
      UserWriteF(" geometry BAD: %d errors",errors);
      fflush(stdout);
    }
    else
    {
      UserWrite(" ok");
    }
  }

  /* check algebraic data structures */
  if (checkalgebra)
  {
    UserWrite(", algebra:");
    fflush(stdout);

    if (errors = CheckAlgebra(theGrid) != GM_OK)
    {
      totalerrors += errors;
      error++;
      UserWriteF(" algebra BAD: %d errors",errors);
      fflush(stdout);
    }
    else
    {
      UserWrite(" ok");
    }
  }

  /* check lists and counters */
  if (checklists)
  {
    UserWrite(", lists:");
    fflush(stdout);

    if (errors = CheckLists(theGrid) != GM_OK)
    {
      totalerrors += errors;
      error++;
      UserWriteF(" lists BAD: %d errors",errors);
      fflush(stdout);
    }
    else
    {
      UserWrite(" ok");
    }
  }

        #ifdef ModelP
  /* check interfaces to other procs */
  if (checkif)
  {
    UserWrite(", interface:");
    fflush(stdout);

    if (errors = CheckInterfaces(theGrid) != GM_OK)
    {
      totalerrors += errors;
      error++;
      UserWriteF(" interfaces BAD: %d errors",errors);
      fflush(stdout);
    }
    else
    {
      UserWrite(" ok");
    }
  }
        #endif

  if (totalerrors)
    UserWriteF(", grid BAD: %d checks has %d totalerrors",error,totalerrors);
  else
    UserWrite(", grid ok");

  return(error);
}

/****************************************************************************/
/*D
   ClearSelection - Clear selection buffer

   SYNOPSIS:
   void ClearSelection (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid structure

   DESCRIPTION:
   This function clears the selection buffer of a multigrid.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/


void ClearSelection (MULTIGRID *theMG)
{
  SELECTIONSIZE(theMG) = 0;
}


/****************************************************************************/
/*D
   AddNodeToSelection - Add node to selection buffer

   SYNOPSIS:
   INT AddNodeToSelection (MULTIGRID *theMG, NODE *theNode);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theNode - node to add

   DESCRIPTION:
   This function adds an node to the selection buffer.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an  error occured.
   D*/
/****************************************************************************/

INT AddNodeToSelection (MULTIGRID *theMG, NODE *theNode)
{
  int i;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=nodeSelection) RETURN(GM_ERROR);
  }
  else SELECTIONMODE(theMG) = nodeSelection;

  g = (SELECTION_OBJECT *) theNode;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g) RETURN(GM_ERROR);

  if (SELECTIONSIZE(theMG)>=MAXSELECTION) RETURN(GM_ERROR);

  SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG)) = g;
  SELECTIONSIZE(theMG)++;
  return(GM_OK);
}


/****************************************************************************/
/*D
   AddElementToSelection - Add element to selection buffer

   SYNOPSIS:
   INT AddElementToSelection (MULTIGRID *theMG, ELEMENT *theElement);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theElement - element to add

   DESCRIPTION:
   This function adds an element to the selection buffer.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an  error occured.
   D*/
/****************************************************************************/

INT AddElementToSelection (MULTIGRID *theMG, ELEMENT *theElement)
{
  int i;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=elementSelection) RETURN(GM_ERROR);
  }
  else SELECTIONMODE(theMG) = elementSelection;

  g = (SELECTION_OBJECT *) theElement;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      SELECTIONSIZE(theMG)--;
      SELECTIONOBJECT(theMG,i) =
        SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG));
      return(GM_OK);
    }

  if (SELECTIONSIZE(theMG)>=MAXSELECTION)
    RETURN(GM_ERROR);

  SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG)) = g;
  SELECTIONSIZE(theMG)++;
  return(GM_OK);
}

/****************************************************************************/
/*D
   AddVectorToSelection - Add vector to selection buffer

   SYNOPSIS:
   INT AddVectorToSelection (MULTIGRID *theMG, VECTOR *theVector);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theVector - vector to add

   DESCRIPTION:
   This function adds a vector to the selection buffer.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an  error occured.
   D*/
/****************************************************************************/

INT AddVectorToSelection (MULTIGRID *theMG, VECTOR *theVector)
{
  int i;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=vectorSelection) RETURN(GM_ERROR);
  }
  else SELECTIONMODE(theMG) = vectorSelection;

  g = (SELECTION_OBJECT *) theVector;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g) return(GM_ERROR);

  if (SELECTIONSIZE(theMG)>=MAXSELECTION) RETURN(GM_ERROR);

  SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG)) = g;
  SELECTIONSIZE(theMG)++;
  return(GM_OK);
}

/****************************************************************************/
/*
   RemoveNodeFromSelection - Remove node from selection buffer

   SYNOPSIS:
   INT RemoveNodeFromSelection (MULTIGRID *theMG, NODE *theNode);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theNode - node to remove

   DESCRIPTION:
   This function removes an node from the selection buffer.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an error occured.
 */
/****************************************************************************/

INT RemoveNodeFromSelection (MULTIGRID *theMG, NODE *theNode)
{
  int i,j,found;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=nodeSelection) RETURN(GM_ERROR);
  }
  else RETURN(GM_ERROR);

  g = (SELECTION_OBJECT *) theNode;
  found = 0;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      found = 1;
      break;
    }

  if (!found) RETURN(GM_ERROR);

  for (j=i+1; j<SELECTIONSIZE(theMG); j++)
    SELECTIONOBJECT(theMG,j-1) = SELECTIONOBJECT(theMG,j);

  SELECTIONSIZE(theMG)--;
  return(GM_OK);
}


/****************************************************************************/
/*
   RemoveElementFromSelection - Remove element from selection buffer

   SYNOPSIS:
   INT RemoveElementFromSelection (MULTIGRID *theMG, ELEMENT *theElement);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theElement - element to remove

   DESCRIPTION:
   This function removes an element from the selection buffer.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an error occured.
 */
/****************************************************************************/

INT RemoveElementFromSelection (MULTIGRID *theMG, ELEMENT *theElement)
{
  int i,j,found;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=elementSelection) RETURN(GM_ERROR);
  }
  else RETURN(GM_ERROR);

  g = (SELECTION_OBJECT *) theElement;
  found = 0;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      found = 1;
      break;
    }

  if (!found) RETURN(GM_ERROR);

  for (j=i+1; j<SELECTIONSIZE(theMG); j++)
    SELECTIONOBJECT(theMG,j-1) = SELECTIONOBJECT(theMG,j);

  SELECTIONSIZE(theMG)--;
  return(GM_OK);
}

/****************************************************************************/
/*
   RemoveVectorFromSelection - Remove vector from selection buffer

   SYNOPSIS:
   INT RemoveVectorFromSelection (MULTIGRID *theMG, VECTOR *theVector);

   PARAMETERS:
   .  theMG - multigrid structure
   .  theVector - vector to remove

   DESCRIPTION:
   This function removes a vector from the selection buffer.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an error occured.
 */
/****************************************************************************/

INT RemoveVectorFromSelection (MULTIGRID *theMG, VECTOR *theVector)
{
  int i,j,found;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=vectorSelection) RETURN(GM_ERROR);
  }
  else RETURN(GM_ERROR);

  g = (SELECTION_OBJECT *) theVector;
  found = 0;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      found = 1;
      break;
    }

  if (!found) RETURN(GM_ERROR);

  for (j=i+1; j<SELECTIONSIZE(theMG); j++)
    SELECTIONOBJECT(theMG,j-1) = SELECTIONOBJECT(theMG,j);

  SELECTIONSIZE(theMG)--;
  return(GM_OK);
}

/****************************************************************************/
/*D
   MinMaxAngle - Determine min and max angle in degrees

   SYNOPSIS:
   INT MinMaxAngle (ELEMENT *theElement, DOUBLE *amin, DOUBLE *amax);

   PARAMETERS:
   .  theElement - element to check
   .  amin - minimum angle
   .  amax - maximum angle

   DESCRIPTION:
   This function determines min and max angle in degrees.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an error occured.
   D*/
/****************************************************************************/

#ifdef __TWODIM__
static INT QualityElement (int n, DOUBLE *xc, DOUBLE *yc, DOUBLE *lmin, DOUBLE *lmax)
{
  int i;
  DOUBLE x1,x2,y1,y2,w,c,s,l,x,y;

  for (i=0; i<n; i++)
  {
    x1 = xc[(i+1)%n]-xc[i];
    x2 = xc[(i+n-1)%n]-xc[i];
    y1 = yc[(i+1)%n]-yc[i];
    y2 = yc[(i+n-1)%n]-yc[i];
    l = sqrt(x1*x1+y1*y1);
    if (l<SMALL_D) RETURN(1);
    c = x1/l; s = -y1/l;
    x = c*x2-s*y2;
    y = s*x2+c*y2;
    if (fabs(x)<=SMALL_D)
    {
      if (y>=0) w=90.0;else w = 270.0;
    }
    else
    {
      if (x>=0)
      {
        w = atan(y/x)/PI*180.0;
        if (w<0) w = w+360.0;
      }
      else
      {
        x = -x;
        w = 180.0-(atan(y/x)/PI*180.0);
      }
    }
    if (w>*lmax) *lmax = w;
    if (w<*lmin) *lmin = w;
  }
  return (0);
}

INT MinMaxAngle (ELEMENT *theElement, DOUBLE *amin, DOUBLE *amax)
{
  int i,n;
  DOUBLE x[4],y[4];

  n = TAG(theElement);
  for (i=0; i<n; i++)
  {
    x[i] = XC(MYVERTEX(CORNER(theElement,i)));
    y[i] = YC(MYVERTEX(CORNER(theElement,i)));
  }
  if (QualityElement(n,x,y,amin,amax)) RETURN(GM_ERROR);
  return(GM_OK);
}
#endif

#ifdef __THREEDIM__
static INT QualityElement (INT type, ELEMENT *element, DOUBLE *angle)
{
  INT i,j,k,errorcode;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE delta[3][DIM],s[DIM],t[DIM];
  DOUBLE help,Scalarprdst,Scalarprd01,Scalarprd02,Scalarprd12;
  DOUBLE_VECTOR theNormal[MAX_CORNERS_OF_ELEM];

  if (type != TETRAHEDRON)
    return (1);

  /* load geometrical data of the corners */
  for (i=0; i<CORNERS_OF_ELEM(element); i++)
    x[i] = CVECT(MYVERTEX(CORNER(element,i)));

  /* calculate all corner-angles */
  for (i=0; i<CORNERS_OF_ELEM(element); i++)
  {
    /* calculate the normalized differencevectors d[] */
    for (j=0; j<(CORNERS_OF_ELEM(element)-1); j++)
      for(k=0; k<DIM; k++)
        delta[j][k] = x[(i+j+1)%CORNERS_OF_ELEM(element)][k] - x[i][k];
    if ((errorcode=V3_Normalize(delta[0]))!=0) return(errorcode);
    if ((errorcode=V3_Normalize(delta[1]))!=0) return(errorcode);
    if ((errorcode=V3_Normalize(delta[2]))!=0) return(errorcode);

    /* calculate necessary scalarproducts */
    V3_SCALAR_PRODUCT(delta[0],delta[1],Scalarprd01);
    V3_SCALAR_PRODUCT(delta[0],delta[2],Scalarprd02);
    V3_SCALAR_PRODUCT(delta[1],delta[2],Scalarprd12);

    /* calculate angle */
    V3_LINCOMB(1.0,delta[1],-1.0*Scalarprd01,delta[0],s);
    V3_LINCOMB(1.0,delta[2],-1.0*Scalarprd02,delta[0],t);
    if ((errorcode=V3_Normalize(s))!=0) return(errorcode);
    if ((errorcode=V3_Normalize(t))!=0) return(errorcode);
    V3_SCALAR_PRODUCT(s,t,Scalarprdst);
    if (Scalarprdst<-1.0) Scalarprdst=-1.0;
    if (Scalarprdst> 1.0) Scalarprdst= 1.0;
    angle[i] = (DOUBLE)acos((float)(Scalarprdst));

    V3_LINCOMB(1.0,delta[2],-1.0*Scalarprd12,delta[1],s);
    V3_LINCOMB(1.0,delta[0],-1.0*Scalarprd01,delta[1],t);
    if ((errorcode=V3_Normalize(s))!=0) return(errorcode);
    if ((errorcode=V3_Normalize(t))!=0) return(errorcode);
    V3_SCALAR_PRODUCT(s,t,Scalarprdst);
    if (Scalarprdst<-1.0) Scalarprdst=-1.0;
    if (Scalarprdst> 1.0) Scalarprdst= 1.0;
    angle[i] += (DOUBLE)acos((float)(Scalarprdst));

    V3_LINCOMB(1.0,delta[0],-1.0*Scalarprd02,delta[2],s);
    V3_LINCOMB(1.0,delta[1],-1.0*Scalarprd12,delta[2],t);
    if ((errorcode=V3_Normalize(s))!=0) return(errorcode);
    if ((errorcode=V3_Normalize(t))!=0) return(errorcode);
    V3_SCALAR_PRODUCT(s,t,Scalarprdst);
    if (Scalarprdst<-1.0) Scalarprdst=-1.0;
    if (Scalarprdst> 1.0) Scalarprdst= 1.0;
    angle[i] += (DOUBLE)acos((float)(Scalarprdst));
    angle[i] -= PI;

    /* normalize sperical angle to unit spere */
    angle[i] /= 4*PI;
  }
  if (TetraSideNormals (element,x,theNormal)) return (2);
  for (i=0; i<EDGES_OF_ELEM(element); i++)
  {
    V3_SCALAR_PRODUCT(theNormal[SIDE_WITH_EDGE(element,i,0)],theNormal[SIDE_WITH_EDGE(element,i,1)],help);
    help = MAX(help,-1.0);
    help = MIN(help, 1.0);
    angle[CORNERS_OF_ELEM(element)+i] = 180.0/PI*acos(-help);
  }

  return(0);
}

INT MinMaxAngle (ELEMENT *theElement, DOUBLE *amin, DOUBLE *amax)
{
  DOUBLE angle[MAX_CORNERS_OF_ELEM+MAX_EDGES_OF_ELEM];
  INT i;

  if (QualityElement(TETRAHEDRON,theElement,angle) != 0) return(GM_ERROR);

  *amin = MAX_D; *amax = -MAX_D;
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    *amax = MAX(*amax,angle[i]);
    *amin = MIN(*amin,angle[i]);
  }

  return(GM_OK);
}
#endif

/****************************************************************************/
/*D
   DefineMGUDBlock - Define block in general MG user data space

   SYNOPSIS:
   INT DefineMGUDBlock (BLOCK_ID id, MEM size);

   PARAMETERS:
   .  id - the id of the block to be allocated
   .  size - size of the data block

   DESCRIPTION:
   This function defines a block in the general MG user data space.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an  error occured.
   D*/
/****************************************************************************/

INT DefineMGUDBlock (BLOCK_ID id, MEM size)
{
  if (DefineBlock(theGenMGUDM,id,size)!=0)
    return (GM_ERROR);

  return (GM_OK);
}


/****************************************************************************/
/*D
   FreeMGUDBlock - Free block in general MG user data space

   SYNOPSIS:
   INT FreeMGUDBlock (BLOCK_ID id);

   PARAMETERS:
   .  id: the id of the block to be allocated

   DESCRIPTION:
   This function frees a block in the general MG user data space.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an error occured.
   D*/
/****************************************************************************/

INT FreeMGUDBlock (BLOCK_ID id)
{
  if (FreeBlock(theGenMGUDM,id)!=0)
    return (GM_ERROR);

  return (GM_OK);
}


/****************************************************************************/
/*D
   GetMGUDBlockDescriptor - Return pointer to block descriptor with id

   SYNOPSIS:
   BLOCK_DESC  *GetMGUDBlockDescriptor (INT where, MULTIGRID *theMG, BLOCK_ID id);

   PARAMETERS:
   .  id - the id of the block to be allocated

   DESCRIPTION:
   This function returns a pointer to the block descriptor with id.

   RETURN VALUE:
   BLOCK_DESC *
   .n   pointer to BLOCK_DESC
   .n   NULL if an error occured.
   D*/
/****************************************************************************/

BLOCK_DESC      *GetMGUDBlockDescriptor (BLOCK_ID id)
{
  return (GetBlockDesc(theGenMGUDM,id));
}

VIRT_HEAP_MGMT *GetGenMGUDM()
{
  return (theGenMGUDM);
}


/****************************************************************************/
/*D
   RenumberNodeElem - Init what is neccessary

   SYNOPSIS:
   INT RenumberNodeElem (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - ptr to multigrid

   DESCRIPTION:
   This function renumbers nodes and elements according to their double-linked
   lists.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   > 0 line in which error occured.
   D*/
/****************************************************************************/

INT RenumberNodeElem (MULTIGRID *theMG)
{
  INT i,nid,eid;
  GRID *theGrid;
  NODE *theNode;
  ELEMENT *theElement;

  nid=eid=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      ID(theElement) = eid++;
    if (i==0)
    {
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        if (OBJT(MYVERTEX(theNode))==BVOBJ)
          ID(theNode) = nid++;
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        if (OBJT(MYVERTEX(theNode))==IVOBJ)
          ID(theNode) = nid++;
    }
    else
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        ID(theNode) = nid++;
  }

  return (0);
}

/****************************************************************************/
/*D
   CheckEnumerationNodeElem - check if numeration is according to RenumberNodeElem

   SYNOPSIS:
   INT CheckEnumerationNodeElem (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - ptr to multigrid

   DESCRIPTION:
   This function checks numbers of nodes and elements;

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   > 0 line in which error occured.
   D*/
/****************************************************************************/

INT CheckEnumerationNodeElem (MULTIGRID *theMG)
{
  INT i,nid,eid;
  GRID *theGrid;
  NODE *theNode;
  ELEMENT *theElement;

  nid=eid=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (ID(theElement)!=eid) return (1);
      eid++;
    }
    if (i==0)
    {
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        if (OBJT(MYVERTEX(theNode))==BVOBJ)
        {
          if (ID(theNode) != nid) return (1);
          nid++;
        }
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        if (OBJT(MYVERTEX(theNode))==IVOBJ)
        {
          if (ID(theNode) != nid) return (1);
          nid++;
        }
    }
    else
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      {
        if (ID(theNode) != nid) return (1);
        nid++;
      }
  }

  return (0);
}

/****************************************************************************/
/*D
   InitUGManager - Init what is neccessary

   SYNOPSIS:
   INT InitUGManager ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes the grid manager.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   > 0 line in which error occured.
   D*/
/****************************************************************************/

INT InitUGManager ()
{
  INT i;

  theGenMGUDM = malloc(SIZEOF_VHM);
  if (theGenMGUDM==NULL)
    return (__LINE__);

  InitVirtualHeapManagement(theGenMGUDM,SIZE_UNKNOWN);

  /* install the /Multigrids directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitUGManager","could not changedir to root");
    return(__LINE__);
  }
  theMGRootDirID = GetNewEnvDirID();
  if (MakeEnvItem("Multigrids",theMGRootDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitUGManager","could not install /Multigrids dir");
    return(__LINE__);
  }
  theMGDirID = GetNewEnvDirID();

  /* init the OBJT management */
  UsedOBJT = 0;
  for (i=0; i<NPREDEFOBJ; i++)
    SET_FLAG(UsedOBJT,1<<i);

  return (GM_OK);
}
