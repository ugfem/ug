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
#include "fifo.h"

#include "devices.h"

#include "evm.h"
#include "gm.h"
#include "rm.h"
#include "misc.h"
#include "dlmgr.h"
#include "algebra.h"
#include "ugm.h"
#include "elements.h"
#include "shapes.h"
#include "refine.h"
#include "domain.h"
#include "pargm.h"

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

REP_ERR_FILE;

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static NODE *CreateNode (GRID *theGrid, VERTEX *vertex, GEOM_OBJECT *Father, INT NodeType, INT with_vector);
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

  /* skip predefined object types, they cannot be re-allocated */
  for (i=NPREDEFOBJ; i<MAXOBJECTS; i++)
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

  /* we cannot release predefined object types! */
  if (type<NPREDEFOBJ)
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
  if ((ds=FMT_S_VERTEX(MGFORMAT(MYMG(theGrid))))>0)
  {
    VDATA(pv) = GetMemoryForObject(MYMG(theGrid),ds,NOOBJ);
    if (VDATA(pv)==NULL) return(NULL);
  }
  else
    VDATA(pv) = NULL;

  /* initialize data */
  CTRL(pv) = 0;
  SETOBJT(pv,BVOBJ);
  SETNOOFNODE(pv,1);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
        #ifdef TOPNODE
  TOPNODE(pv) = NULL;
        #endif
  for (i=0; i<DIM; i++) LCVECT(pv)[i] = 0.0;
  SETONEDGE(pv,0);
  SETMOVE(pv,DIM_OF_BND);
        #ifdef ModelP
  DDD_AttrSet(PARHDRV(pv),GRID_ATTR(theGrid));
  SETVXPRIO(pv,PrioMaster);
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
  if ((ds=FMT_S_VERTEX(MGFORMAT(MYMG(theGrid))))>0)
  {
    VDATA(pv) = GetMemoryForObject(MYMG(theGrid),ds,NOOBJ);
    if (VDATA(pv)==NULL) return(NULL);
  }
  else
    VDATA(pv) = NULL;

  /* initialize data */
  CTRL(pv) = 0;
  SETOBJT(pv,IVOBJ);
  SETNOOFNODE(pv,1);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
        #ifdef TOPNODE
  TOPNODE(pv) = NULL;
        #endif
  SETMOVE(pv,DIM);
        #ifdef ModelP
  DDD_AttrSet(PARHDRV(pv),GRID_ATTR(theGrid));
  SETVXPRIO(pv,PrioMaster);
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
   static NODE *CreateNode (GRID *theGrid, VERTEX *vertex,
   GEOM_OBJECT *Father, INT NodeType);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  vertex  - vertex of the node
   .  FatherNode - father node (may be NULL)
   .  NodeType - node type (CORNER_NODE..., cf. gm.h)

   DESCRIPTION:
   This function creates and initializes a new node structure
   and returns a pointer to it.

   RETURN VALUE:
   NODE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

static NODE *CreateNode (GRID *theGrid, VERTEX *vertex,
                         GEOM_OBJECT *Father, INT NodeType, INT with_vector)
{
  NODE *pn;
  VECTOR *pv;
  INT size;

  size = sizeof(NODE);
  if (!VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC)) size -= sizeof(VECTOR *);
  if (NDATA_DEF_IN_GRID(theGrid)) size += sizeof(void *);
  if (NELIST_DEF_IN_GRID(theGrid)) size += sizeof(void *);

  pn = (NODE *)GetMemoryForObject(MYMG(theGrid),size,NDOBJ);
  if (pn==NULL) return(NULL);

  /* initialize data */
  SETOBJT(pn,NDOBJ);
  SETLEVEL(pn,theGrid->level);
        #ifdef ModelP
  DDD_AttrSet(PARHDR(pn),GRID_ATTR(theGrid));
  SETPRIO(pn,PrioMaster);
        #endif
  ID(pn) = (theGrid->mg->nodeIdCounter)++;
  START(pn) = NULL;
  SONNODE(pn) = NULL;
  if (NELIST_DEF_IN_GRID(theGrid)) NDATA(pn) = NULL;
  MYVERTEX(pn) = vertex;
  /* priliminary */
  if (Father != NULL)
    if ((OBJT(Father) == IEOBJ) || (OBJT(Father) == BEOBJ))
      Father = NULL;
  SETNFATHER(pn,Father);
  SETNTYPE(pn,NodeType);
  if (VFATHER(vertex) != NULL)
    SETNSUBDOM(pn,SUBDOMAIN(VFATHER(vertex)));
  else if (Father != NULL) {
    if (OBJT(Father) == NDOBJ)
      SETNSUBDOM(pn,NSUBDOM((NODE *)Father));
    else if (OBJT(Father) == EDOBJ)
      SETNSUBDOM(pn,EDSUBDOM((EDGE *)Father));
  }
  else
    SETNSUBDOM(pn,0);

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
    if (with_vector)
    {
      if (CreateVector (theGrid,NODEVEC,(GEOM_OBJECT *)pn,&pv))
      {
        DisposeNode (theGrid,pn);
        return (NULL);
      }
      NVECTOR(pn) = pv;
    }
    else
      NVECTOR(pn) = NULL;

  if (NDATA_DEF_IN_GRID(theGrid)) {
    NDATA(pn) = (void *) GetMemoryForObject(theGrid->mg,
                                            NDATA_DEF_IN_GRID(theGrid),-1);
    if (NDATA(pn) == NULL) {
      DisposeNode (theGrid,pn);
      return (NULL);
    }
  }

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

  theVertex = MYVERTEX(FatherNode);

  pn = CreateNode(theGrid,theVertex,(GEOM_OBJECT *)FatherNode,CORNER_NODE,1);
  if (pn == NULL)
    return(NULL);
  SONNODE(FatherNode) = pn;
  if (NOOFNODE(theVertex)<NOOFNODEMAX-1)
    INCNOOFNODE(theVertex);
  else
    ASSERT(0);

        #ifdef TOPNODE
  TOPNODE(theVertex) = pn;
        #endif

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
  INT n,co0,co1,move,part;

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
      if (BNDP_BndPDesc(bndp,&move,&part))
        return(NULL);
      SETMOVE(theVertex,move);
      V_BNDP(theVertex) = bndp;
      V_DIM_COPY(bnd_global,CVECT(theVertex));
      local = LCVECT(theVertex);
      V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
      if (diff > MAX_PAR_DIST)
      {
        SETMOVED(theVertex,1);
        CORNER_COORDINATES(theElement,n,x);
        UG_GlobalToLocal(n,(const DOUBLE **)x,bnd_global,local);
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

  /* set MIDNODE pointer */
  theEdge = GetEdge(CORNER(theElement,co0),CORNER(theElement,co1));
  ASSERT(theEdge!=NULL);

  /* allocate node */
  theNode = CreateNode(theGrid,theVertex,(GEOM_OBJECT *)theEdge,MID_NODE,1);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }

  MIDNODE(theEdge) = theNode;
        #ifdef TOPNODE
  TOPNODE(theVertex) = theNode;
        #endif

  if (OBJT(theVertex) == BVOBJ)
    PRINTDEBUG(dom,1,(" MidPoint %d %f %f %f\n",ID(theNode),
                      bnd_global[0],
                      bnd_global[1],
                      bnd_global[2]));

  PRINTDEBUG(dddif,1,(PFMT " CreateMidNode(): n=" ID_FMTX
                      " NTYPE=%d OBJT=%d father " ID_FMTX " \n",
                      me,ID_PRTX(theNode),NTYPE(theNode),
                      OBJT(NFATHER(theNode)),NFATHER(theNode)));

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

NODE *GetMidNode (ELEMENT *theElement, INT edge)
{
  EDGE *theEdge;
  NODE *theNode;
  VERTEX *theVertex;

  HEAPFAULT(theElement);
  theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,edge,0)),
                    CORNER(theElement,CORNER_OF_EDGE(theElement,edge,1)));
  if (theEdge == NULL) return(NULL);
  theNode = MIDNODE(theEdge);
  if (theNode == NULL) return(NULL);
  theVertex = MYVERTEX(theNode);
  if (VFATHER(theVertex) == NULL) {
    VFATHER(theVertex) = theElement;
    SETONEDGE(theVertex,edge);
    V_DIM_LINCOMB(0.5,
                  LOCAL_COORD_OF_ELEM(theElement,
                                      CORNER_OF_EDGE(theElement,edge,0)),
                  0.5,
                  LOCAL_COORD_OF_ELEM(theElement,
                                      CORNER_OF_EDGE(theElement,edge,1)),
                  LCVECT(theVertex));
  }
  return(theNode);
}

static INT SideOfNbElement(ELEMENT *theElement, INT side)
{
  ELEMENT *nb;
  NODE *nd[MAX_CORNERS_OF_SIDE];
  INT i,j,m,n,num;

  nb = NBELEM(theElement,side);
  if (nb == NULL) return(MAX_SIDES_OF_ELEM);

  for (j=0; j<SIDES_OF_ELEM(nb); j++)
    if (NBELEM(nb,j) == theElement) return(j);

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
  INT n,j,k,move,part;

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
  if (OBJT(theElement) == BEOBJ) {
    bnds = ELEM_BNDS(theElement,side);
    if (bnds != NULL) {
      if (n == 3)
        bnd_local[0] = bnd_local[1] = 0.33333333333333;
      else if (n == 4)
        bnd_local[0] = bnd_local[1] = 0.5;
      bndp = BNDS_CreateBndP(MGHEAP(MYMG(theGrid)),bnds,bnd_local);
      if (bndp != NULL) {
        theVertex = CreateBoundaryVertex(theGrid);
        if (theVertex == NULL)
          return(NULL);
        if (BNDP_BndPDesc(bndp,&move,&part))
          return(NULL);
        SETMOVE(theVertex,move);
        if (BNDP_Global(bndp,bnd_global))
          return(NULL);
        V_BNDP(theVertex) = bndp;
        V_DIM_COPY(bnd_global,CVECT(theVertex));
        V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
        if (diff > MAX_PAR_DIST) {
          SETMOVED(theVertex,1);
          CORNER_COORDINATES(theElement,k,x);
          UG_GlobalToLocal(k,(const DOUBLE **)x,bnd_global,local);
        }
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
  theNode = CreateNode(theGrid,theVertex,
                       (GEOM_OBJECT *)theElement,SIDE_NODE,1);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
        #ifdef TOPNODE
  TOPNODE(theVertex) = theNode;
        #endif
  theGrid->status |= 1;

  return(theNode);
}

NODE *GetSideNode (ELEMENT *theElement, INT side)
{
  ELEMENT *theFather;
  NODE *theNode;
  NODE *MidNodes[MAX_EDGES_OF_SIDE];
  VERTEX *theVertex;
  LINK *theLink0,*theLink1,*theLink2,*theLink3;
  DOUBLE fac,*local;
  INT i,n,m;

  m = EDGES_OF_SIDE(theElement,side);
  n = 0;
  for (i=0; i<m; i++) {
    theNode = GetMidNode(theElement,EDGE_OF_SIDE(theElement,side,i));
    if (theNode != NULL)
      MidNodes[n++] = theNode;
                #ifndef ModelP
    else return(NULL);
                #endif
  }
  if (n == 4) {
    for (theLink0=START(MidNodes[0]); theLink0!=NULL;
         theLink0=NEXT(theLink0)) {
      theNode = NBNODE(theLink0);
      if (NTYPE(theNode) != SIDE_NODE)
        continue;
      for (theLink1=START(MidNodes[1]); theLink1!=NULL;
           theLink1=NEXT(theLink1)) {
        if (theNode != NBNODE(theLink1))
          continue;
        for (theLink2=START(MidNodes[2]); theLink2!=NULL;
             theLink2=NEXT(theLink2)) {
          if (theNode != NBNODE(theLink2))
            continue;
          for (theLink3=START(MidNodes[3]); theLink3!=NULL;
               theLink3=NEXT(theLink3)) {
            if (theNode != NBNODE(theLink3))
              continue;
            theVertex = MYVERTEX(theNode);
            theFather = VFATHER(theVertex);
            if (theFather == theElement)
              assert(ONSIDE(theVertex) == side);
            else if (theFather == NBELEM(theElement,side))
              SETONNBSIDE(theVertex,side);
            else if (theFather == NULL) {
              VFATHER(theVertex) = theElement;
              SETONSIDE(theVertex,side);
              SETONNBSIDE(theVertex,
                          SideOfNbElement(theElement,side));
              fac = 1.0 / n;
              local = LCVECT(theVertex);
              V_DIM_CLEAR(local);
              for (i=0; i<n; i++) {
                V_DIM_LINCOMB(1.0,local,fac,
                              LOCAL_COORD_OF_ELEM(theElement,
                                                  CORNER_OF_SIDE(theElement,side,i)),
                              local);
              }
            }
            else
              assert(0);
            return(theNode);
          }
        }
      }
    }
  }
  else if (n == 3) {
    for (theLink0=START(MidNodes[0]); theLink0!=NULL;
         theLink0=NEXT(theLink0)) {
      theNode = NBNODE(theLink0);
      if (NTYPE(theNode) != SIDE_NODE)
        continue;
      for (theLink1=START(MidNodes[1]); theLink1!=NULL;
           theLink1=NEXT(theLink1)) {
        if (theNode != NBNODE(theLink1))
          continue;
        for (theLink2=START(MidNodes[2]); theLink2!=NULL;
             theLink2=NEXT(theLink2)) {
          if (theNode != NBNODE(theLink2))
            continue;
          theVertex = MYVERTEX(theNode);
          theFather = VFATHER(theVertex);
          if (theFather == theElement) {
            if (ONSIDE(theVertex) == side)
              return(theNode);
          }
          else if (theFather == NBELEM(theElement,side)) {
            SETONNBSIDE(theVertex,side);
            return(theNode);
          }
          else if (theFather == NULL) {
            VFATHER(theVertex) = theElement;
            SETONSIDE(theVertex,side);
            SETONNBSIDE(theVertex,
                        SideOfNbElement(theElement,side));
            fac = 1.0 / n;
            local = LCVECT(theVertex);
            V_DIM_CLEAR(local);
            for (i=0; i<n; i++) {
              V_DIM_LINCOMB(1.0,local,fac,
                            LOCAL_COORD_OF_ELEM(theElement,
                                                CORNER_OF_SIDE(theElement,side,i)),
                            local);
            }
            return(theNode);
          }
        }
      }
    }
  }
    #ifdef ModelP
  else if (n == 2) {
    for (theLink0=START(MidNodes[0]); theLink0!=NULL;
         theLink0=NEXT(theLink0)) {
      theNode = NBNODE(theLink0);
      if (NTYPE(theNode) != SIDE_NODE)
        continue;
      for (theLink1=START(MidNodes[1]); theLink1!=NULL;
           theLink1=NEXT(theLink1)) {
        if (theNode != NBNODE(theLink1))
          continue;
        theVertex = MYVERTEX(theNode);
        theFather = VFATHER(theVertex);
        if (theFather == theElement) {
          if (ONSIDE(theVertex) == side)
            return(theNode);
        }
        else if (theFather == NBELEM(theElement,side)) {
          SETONNBSIDE(theVertex,side);
          return(theNode);
        }
      }
    }
  }
  else if (n == 1) {
    for (theLink0=START(MidNodes[0]); theLink0!=NULL;
         theLink0=NEXT(theLink0)) {
      theNode = NBNODE(theLink0);
      if (NTYPE(theNode) != SIDE_NODE)
        continue;
      theVertex = MYVERTEX(theNode);
      theFather = VFATHER(theVertex);
      if (theFather == theElement) {
        if (ONSIDE(theVertex) == side)
          return(theNode);
      }
      else if (theFather == NBELEM(theElement,side)) {
        SETONNBSIDE(theVertex,side);
        return(theNode);
      }
    }
  }
    #endif

  return(NULL);
}

INT GetSideIDFromScratch (ELEMENT *theElement, NODE *theNode)
{
  ELEMENT *theFather;
  NODE *nd[MAX_EDGES_OF_ELEM];
  INT i,j,k,l,cnt;

  ASSERT(NTYPE(theNode) == SIDE_NODE);

  theFather = EFATHER(theElement);
  for (i=0; i<EDGES_OF_ELEM(theFather); i++)
    nd[i] = MIDNODE(
      GetEdge(CORNER(theFather,CORNER_OF_EDGE(theFather,i,0)),
              CORNER(theFather,CORNER_OF_EDGE(theFather,i,1))));
  for (j=0; j<SIDES_OF_ELEM(theElement); j++) {
    for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++)
      if (theNode == CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
        break;
    if (l == CORNERS_OF_SIDE(theElement,j)) continue;
    for (i=0; i<SIDES_OF_ELEM(theFather); i++) {
      cnt = 0;
      for (k=0; k<EDGES_OF_SIDE(theFather,i); k++)
        for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++) {
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
            cnt++;
          if (cnt == 2)
            return(i);
        }
    }
  }
  return(SIDES_OF_ELEM(theFather));
}

#endif /* __THREEDIM__ */

/****************************************************************************/
/*																			*/
/* Function:  GetCenterNode	                                                                                        */
/*																			*/
/* Purpose:   get the center node of an element of next finer level         */
/*																			*/
/* return:	  NODE* : pointer to center node								*/
/*			  NULL	: no node found                                                                         */
/*																			*/
/****************************************************************************/

NODE *GetCenterNode (ELEMENT *theElement)
{
  INT i,j;
  NODE    *theNode;
  ELEMENT *SonList[MAX_SONS],*theSon;

  theNode = NULL;
  if (GetAllSons(theElement,SonList) != GM_OK) assert(0);

  for (i=0; SonList[i]!=NULL; i++)
  {
    theSon = SonList[i];
    for (j=0; j<CORNERS_OF_ELEM(theSon); j++)
    {
      theNode = CORNER(theSon,j);
      if (NTYPE(theNode) == CENTER_NODE)
      {
        assert(VFATHER(MYVERTEX(theNode)) == theElement);
        return (theNode);
      }
    }
  }
  return (NULL);
}

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
  DOUBLE_VECTOR diff;
  INT n,j,moved;
  VERTEX *theVertex,*VertexOnEdge[MAX_EDGES_OF_ELEM];
  NODE *theNode;
  EDGE *theEdge;
  DOUBLE fac;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE len_opp, len_bnd;

  /* check if moved side nodes exist */
  CORNER_COORDINATES(theElement,n,x);
  moved = 0;
  if (OBJT(theElement) == BEOBJ) {
    for (j=0; j<EDGES_OF_ELEM(theElement); j++) {
      theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),
                      CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
      ASSERT(theEdge != NULL);
      theNode = MIDNODE(theEdge);
      if (theNode == NULL)
        VertexOnEdge[j] = NULL;
      else {
        VertexOnEdge[j] = MYVERTEX(theNode);
        moved += MOVED(VertexOnEdge[j]);
      }
    }
                #ifndef ModelP
    if (moved == 1) {
      for (j=0; j<EDGES_OF_ELEM(theElement); j++)
        if (VertexOnEdge[j] != NULL)
          if (MOVED(VertexOnEdge[j])) break;
      theVertex = VertexOnEdge[OPPOSITE_EDGE(theElement,j)];
      if (theVertex != NULL) {
        V_DIM_LINCOMB(0.5,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)))),
                      0.5,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)))),
                      diff);
        V_DIM_LINCOMB(1.0,CVECT(VertexOnEdge[j]),-1.0,diff,diff);
        /* scale diff according to length of edges */
        V_DIM_EUKLIDNORM_OF_DIFF(
          CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)))),
          CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)))),len_bnd);
        V_DIM_EUKLIDNORM_OF_DIFF(
          CVECT(MYVERTEX(CORNER(theElement,
                                CORNER_OF_EDGE(theElement,OPPOSITE_EDGE(theElement,j),0)))),
          CVECT(MYVERTEX(CORNER(theElement,
                                CORNER_OF_EDGE(theElement,OPPOSITE_EDGE(theElement,j),1)))),
          len_opp);
        V_DIM_SCALE(len_opp/len_bnd,diff);
        global = CVECT(theVertex);
        V_DIM_LINCOMB(1.0,global,0.5,diff,global);
        SETMOVED(VertexOnEdge[OPPOSITE_EDGE(theElement,j)],1);
        UG_GlobalToLocal(n,(const DOUBLE **)x,global,LCVECT(theVertex));
        SETONEDGE(theVertex,OPPOSITE_EDGE(theElement,j));
        VFATHER(theVertex) = theElement;
      }
    }
            #endif
  }

  theVertex = CreateInnerVertex(theGrid);
  if (theVertex==NULL)
    return(NULL);
  VFATHER(theVertex) = theElement;

  theNode = CreateNode(theGrid,theVertex,
                       (GEOM_OBJECT *)theElement,CENTER_NODE,1);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  global = CVECT(theVertex);
  local = LCVECT(theVertex);
  V_DIM_CLEAR(local);
  fac = 1.0 / n;
  for (j=0; j<n; j++)
    V_DIM_LINCOMB(1.0,local,
                  fac,LOCAL_COORD_OF_ELEM(theElement,j),local);
  LOCAL_TO_GLOBAL(n,x,local,global);
  if (moved) {
    V_DIM_CLEAR(diff);
    for (j=0; j<EDGES_OF_ELEM(theElement); j++)
      if (VertexOnEdge[j] != NULL) {
        V_DIM_COPY(CVECT(VertexOnEdge[j]),diff);
        V_DIM_LINCOMB(1.0,diff,-0.5,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)))),diff);
        V_DIM_LINCOMB(1.0,diff,-0.5,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)))),diff);
        V_DIM_LINCOMB(0.5,diff,1.0,global,global);
      }
    UG_GlobalToLocal(n,(const DOUBLE **)x,global,local);
    LOCAL_TO_GLOBAL(n,x,local,diff);
    SETMOVED(theVertex,1);
  }
        #ifdef TOPNODE
  TOPNODE(theVertex) = theNode;
        #endif
  theGrid->status |= 1;

  return(theNode);
}


/****************************************************************************/
/*D
   GetNodeContext - get all nodes related to this element on next level

   SYNOPSIS:
   INT GetNodeContext (ELEMENT *theElement, NODE **theElementContext)

   PARAMETERS:
   .  theElement - element for context
   .  theElementContext - node context of this element

   DESCRIPTION:
   This function returns the nodes related to the element on the next
   finer level. The ordering is according to the reference numbering.

   RETURN VALUE:
   INT
   .n   GM_OK    if ok
   .n   != GM_OK if not ok
   D*/
/****************************************************************************/

INT GetNodeContext (ELEMENT *theElement, NODE **theElementContext)
{
  NODE *theNode, **MidNodes, **CenterNode;
  EDGE *theEdge;
  INT i,Corner0, Corner1;
        #ifdef __THREEDIM__
  NODE **SideNodes;
        #endif

  /* reset context */
  for(i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++)
    theElementContext[i] = NULL;

  /* is element to refine */
  if (!IS_REFINED(theElement)) return(GM_OK);

  /* get corner nodes */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    theNode = CORNER(theElement,i);
    theElementContext[i] = SONNODE(theNode);
  }

  /* check for midpoint nodes */
  MidNodes = theElementContext+CORNERS_OF_ELEM(theElement);
  for (i=0; i<EDGES_OF_ELEM(theElement); i++)
  {
    Corner0 = CORNER_OF_EDGE(theElement,i,0);
    Corner1 = CORNER_OF_EDGE(theElement,i,1);

    theEdge = GetEdge(CORNER(theElement,Corner0),
                      CORNER(theElement,Corner1));
    ASSERT(theEdge != NULL);

    MidNodes[i] = MIDNODE(theEdge);
  }

        #ifdef __THREEDIM__
  SideNodes = theElementContext+CORNERS_OF_ELEM(theElement)+
              EDGES_OF_ELEM(theElement);
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
#ifdef TET_RULESET
    /* no side nodes for triangular sides yet */
    if (CORNERS_OF_SIDE(theElement,i) == 3) continue;
#endif
    /* check for side node */
    SideNodes[i] = GetSideNode(theElement,i);
  }
        #endif

  /* check for center node */
  CenterNode = MidNodes+CENTER_NODE_INDEX(theElement);
  CenterNode[0] = GetCenterNode(theElement);

  return(GM_OK);
}

/****************************************************************************/
/*D
   GetNbSideByNodes - return matching side of the neighboring element

   SYNOPSIS:
   void GetNbSideByNodes (ELEMENT *theNeighbor, INT *nbside, ELEMENT *theElement, INT side);

   PARAMETERS:
   .  theNeighbor - element to test for matching side
   .  nbside - the matching side
   .  theElement - element with side to match
   .  side - side of element to match

   DESCRIPTION:
   This function computes the matching side of neighboring element.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void GetNbSideByNodes (ELEMENT *theNeighbor, INT *nbside, ELEMENT *theElement, INT side)
{
  INT i,k,l,ec,nc;

  ec = CORNERS_OF_SIDE(theElement,side);

  for (i=0; i<SIDES_OF_ELEM(theNeighbor); i++)
  {
    nc = CORNERS_OF_SIDE(theNeighbor,i);
    if (ec != nc) continue;

    for (k=0; k<nc; k++)
      if (CORNER_OF_SIDE_PTR(theElement,side,0) ==
          CORNER_OF_SIDE_PTR(theNeighbor,i,k))
        break;
    if (k == nc) continue;

    for (l=1; l<ec; l++)
    {
      if (CORNER_OF_SIDE_PTR(theElement,side,l) !=
          CORNER_OF_SIDE_PTR(theNeighbor,i,(nc+k-l)%nc))
        break;
    }
    if (l == ec)
    {
      *nbside = i;
      return;
    }
  }

  /* no side of the neighbor matches */
  assert(0);
}

/****************************************************************************/
/*D
   SonEdge - Return pointer to son edge if it exists

   SYNOPSIS:
   EDGE *SonEdge (EDGE *theEdge);

   PARAMETERS:
   .  theEdge - edge for which son is searched

   DESCRIPTION:
   This function returns the pointer to the son edge if it exists.

   RETURN VALUE:
   EDGE *
   .n   pointer to specified object
   .n   NULL if not found
   D*/
/****************************************************************************/

EDGE *GetSonEdge (EDGE *theEdge)
{
  EDGE *SonEdge=NULL;
  NODE *Node0,*Node1,*SonNode0,*SonNode1;

  Node0 = NBNODE(LINK0(theEdge));
  Node1 = NBNODE(LINK1(theEdge));

  SonNode0 = SONNODE(Node0);
  SonNode1 = SONNODE(Node1);

  if (SonNode0!=NULL && SonNode1!=NULL)
    SonEdge = GetEdge(SonNode0,SonNode1);

  return(SonEdge);
}


/****************************************************************************/
/*D
   FatherEdge - Return pointer to father edge if it exists

   SYNOPSIS:
   EDGE *FatherEdge (NODE **SideNodes, INT ncorners, NODE **Nodes, EDGE *theEdge)

   PARAMETERS:
   .  SideNodes - nodes of the side
   .  ncorners - number of sidenodes
   .  Nodes - corners of edge for which father is searched
   .  theEdge - edge for which father is searched

   DESCRIPTION:
   This function returns the pointer to the father edge if it exists.

   RETURN VALUE:
   EDGE *
   .n   pointer to specified object
   .n   NULL if not found
   D*/
/****************************************************************************/

#ifdef __THREEDIM__
EDGE *FatherEdge (NODE **SideNodes, INT ncorners, NODE **Nodes, EDGE *theEdge)
{
  INT pos0,pos1;
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

  for (pos1=0; pos1<MAX_SIDE_NODES; pos1++)
    if (SideNodes[pos1] == Nodes[1])
      break;
  ASSERT(pos1<MAX_SIDE_NODES);

  switch (NTYPE(Nodes[0]))
  {
  case (CORNER_NODE) :

    ASSERT(pos0<ncorners);
    if ( ((pos0+1)%ncorners == pos1) ||
         (pos0+ncorners == pos1) )
    {
      ASSERT(OBJT(NFATHER(SideNodes[(pos0+1)%ncorners])) == NDOBJ);
      fatherEdge = GetEdge((NODE *)NFATHER(Nodes[0]),
                           (NODE *)NFATHER(SideNodes[(pos0+1)%ncorners]));
      ASSERT(fatherEdge!=NULL);
    }

    if ( ((pos0-1+ncorners)%ncorners == pos1) ||
         ((pos0-1+ncorners)%ncorners+ncorners == pos1) )
    {
      ASSERT(OBJT(NFATHER(SideNodes[(pos0-1+ncorners)%ncorners])) == NDOBJ);
      fatherEdge = GetEdge((NODE *)NFATHER(Nodes[0]),
                           (NODE *)NFATHER(SideNodes[(pos0-1+ncorners)%ncorners]));
      ASSERT(fatherEdge!=NULL);
    }

    break;

  case (MID_NODE) :

    ASSERT(pos0>=ncorners);
    ASSERT(pos0<2*ncorners);

    if ((pos0+1)%ncorners == pos1)
    {
      ASSERT(OBJT(NFATHER(SideNodes[pos0%ncorners])) == NDOBJ);
      fatherEdge = GetEdge((NODE *)NFATHER(SideNodes[pos0%ncorners]),
                           (NODE *)NFATHER(Nodes[1]));
      ASSERT(fatherEdge!=NULL);
    }

    if (pos0%ncorners == pos1)
    {
      ASSERT(OBJT(NFATHER(SideNodes[(pos0+1)%ncorners])) == NDOBJ);
      fatherEdge = GetEdge((NODE *)NFATHER(SideNodes[(pos0+1)%ncorners]),
                           (NODE *)NFATHER(Nodes[1]));
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

  IFDEBUG(dddif,1)
  INT i;
  EDGE* edge0, *edge1;

  edge0 = edge1 = NULL;

  /* test whether theEdge lies above fatherEdge */
  if (fatherEdge!=NULL)
  {
    if (MIDNODE(fatherEdge)!=NULL)
    {
      edge0 = GetEdge(MIDNODE(fatherEdge),SONNODE(NBNODE(LINK0(fatherEdge))));
      edge1 = GetEdge(MIDNODE(fatherEdge),SONNODE(NBNODE(LINK1(fatherEdge))));
    }
    else
      edge0 = GetEdge(SONNODE(NBNODE(LINK0(fatherEdge))),
                      SONNODE(NBNODE(LINK1(fatherEdge))));

    IFDEBUG(dddif,5)
    UserWriteF("%4d: fatherEdge=%x theEdge=%x edge0=%x edge1=%x\n",me,
               fatherEdge,theEdge,edge0,edge1);
    UserWriteF("%4d: Nodes[0]=%d Nodes[1]=%d\n",me,ID(Nodes[0]),ID(Nodes[1]));

    UserWriteF("SideNodes\n");
    for (i=0; i<MAX_SIDE_NODES; i++) UserWriteF(" %5d",i);
    UserWriteF("\n");
    for (i=0; i<MAX_SIDE_NODES; i++)
      if (SideNodes[i]!=NULL) UserWriteF(" %5d",ID(SideNodes[i]));
    UserWriteF("\n");
    ENDDEBUG

    assert(edge0==theEdge || edge1==theEdge);
  }
  ENDDEBUG

  return(fatherEdge);
}
#endif

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
   .  subdom_id - id for corresponding subdomain

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
  INT part;

  /* check if edge exists already */
  if( (pe = GetEdge(from, to)) != NULL )
  {
    if (NO_OF_ELEM(pe)<NO_OF_ELEM_MAX-1)
      INC_NO_OF_ELEM(pe);
    else
      ASSERT(0);

    return(pe);
  }

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    pe = GetMemoryForObject(theGrid->mg,sizeof(EDGE),EDOBJ);
  else
    pe = GetMemoryForObject(theGrid->mg,sizeof(EDGE)-sizeof(VECTOR*),EDOBJ);
  if (pe==NULL) return(NULL);

  /* initialize data */
  link0 = LINK0(pe);
  link1 = LINK1(pe);
  SETOBJT(pe,EDOBJ);
  SETLOFFSET(link0,0);
        #ifdef _DEBUG_CW_
  SETOBJT(link1,LIOBJ);
        #endif
  SETLOFFSET(link1,1);

  SETLEVEL(pe,GLEVEL(theGrid));
        #if (defined ModelP) && (defined __THREEDIM__)
  DDD_AttrSet(PARHDR(pe), GRID_ATTR(theGrid));
  SETPRIO(pe,PrioMaster);
        #endif
  NBNODE(link0) = to;
  NBNODE(link1) = from;
  SET_NO_OF_ELEM(pe,1);
  SETEDGENEW(pe,1);
  SETEDSUBDOM(pe,NSUBDOM(from));
  if (NSUBDOM(from) != NSUBDOM(to))
    SETEDSUBDOM(pe,0);
  else {
    if ((OBJT(MYVERTEX(from)) == BVOBJ) && (OBJT(MYVERTEX(to)) == BVOBJ)) {
      /* in parallel, it cannot be guaranteed that the boundary
         information is available */
                    #ifdef ModelP
      if (V_BNDP(MYVERTEX(from)) == NULL)
        SETEDSUBDOM(pe,0);
      else if (V_BNDP(MYVERTEX(to)) == NULL)
        SETEDSUBDOM(pe,0);
      else
                        #endif
      if (BNDP_BndEDesc(V_BNDP(MYVERTEX(from)),
                        V_BNDP(MYVERTEX(to)),&part) == 0)
        SETEDSUBDOM(pe,0);
    }
  }

  /* create vector if */
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    if (with_vector)
    {
      if (CreateVector (theGrid,EDGEVEC,(GEOM_OBJECT *)pe,&pv))
      {
        DisposeEdge (theGrid,pe);
        return (NULL);
      }
      EDVECTOR(pe) = pv;
    }
    else
      EDVECTOR(pe) = NULL;

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

ELEMENT *CreateElement (GRID *theGrid, INT tag, INT objtype, NODE **nodes,
                        ELEMENT *Father, INT with_vector)
{
  ELEMENT *pe;
  EDGE *ed;
  INT i,s_id;
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
  DDD_AttrSet(PARHDRE(pe),GRID_ATTR(theGrid));
  SETEPRIO(pe,PrioMaster);
  PARTITION(pe) = me;
        #endif
  SETEBUILDCON(pe,1);
  ID(pe) = (theGrid->mg->elemIdCounter)++;

  /* subdomain id */
  s_id = (Father != NULL) ? SUBDOMAIN(Father) : 0;
  SETSUBDOMAIN(pe,s_id);

  /* set corner nodes */
  for (i=0; i<CORNERS_OF_ELEM(pe); i++)
    SET_CORNER(pe,i,nodes[i]);

  /* create edges */
  for (i=0; i<EDGES_OF_ELEM(pe); i++)
    if ((ed=CreateEdge(theGrid,
                       nodes[CORNER_OF_EDGE(pe,i,0)],
                       nodes[CORNER_OF_EDGE(pe,i,1)],
                       with_vector)) == NULL)
    {
      DisposeElement(theGrid,pe,TRUE);
      return(NULL);
    }

  /* create element vector if */
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
    if (with_vector)
    {
      if (CreateVector (theGrid,ELEMVEC,(GEOM_OBJECT *)pe,&pv))
      {
        DisposeElement(theGrid,pe,TRUE);
        return (NULL);
      }
      SET_EVECTOR(pe,pv);
    }
    else
      SET_EVECTOR(pe,NULL);

  if (EDATA_DEF_IN_GRID(theGrid)) {
    q = (void *) GetMemoryForObject(theGrid->mg,EDATA_DEF_IN_GRID(theGrid),-1);
    if (q == NULL) {
      DisposeElement(theGrid,pe,TRUE);
      return (NULL);
    }
    SET_EDATA(pe,q);
  }

  /* create side vectors if */
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
    for (i=0; i<SIDES_OF_ELEM(pe); i++)
      if (with_vector)
      {
        if (CreateSideVector (theGrid,i,(GEOM_OBJECT *)pe,&pv))
        {
          DisposeElement(theGrid,pe,TRUE);
          return (NULL);
        }
        SET_SVECTOR(pe,i,pv);
      }
      else
        SET_SVECTOR(pe,i,NULL);

  /* insert in element list */
  GRID_LINK_ELEMENT(theGrid,pe,PrioMaster);

  SET_EFATHER(pe,Father);
  if (theGrid->level>0)
  {
    INT where = PRIO2INDEX(PrioMaster);

#ifndef ModelP
    ASSERT(Father != NULL);
#endif
    if (Father != NULL)
    {
      if (SON(Father,where) == NULL)
        SET_SON(Father,where,pe);
      SETNSONS(Father,NSONS(Father)+1);
    }
  }

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
   Here also the side vector (iff at all) is inspected in 'ReinspectSonSideVector'.
   The latter function eventually reallocates the vector if its size has changed and
   sets the VBUILDCON flag in the vector. The connections of the old vector are
   thereby disposed. The refine-module which is calling
   'CreateSonElementSide' will finally call 'GridCreateConnection' to reinstall
   the connections of the side-vector.

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
  VECTOR *vec;

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

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
  {
    vec = SVECTOR(theSon,son_side);
    ReinspectSonSideVector(theGrid,theSon,son_side,&vec);
    SET_SVECTOR(theSon,son_side,vec);
  }

  return(GM_OK);
}

/****************************************************************************/
/*D
   CreateNewLevel - Return pointer to new grid structure

   SYNOPSIS:
   GRID *CreateNewLevel (MULTIGRID *theMG, INT algebraic);

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

GRID *CreateNewLevel (MULTIGRID *theMG, INT algebraic)
{
  GRID *theGrid;
  INT l;

  if (theMG->bottomLevel>theMG->topLevel && algebraic) return (NULL);
  if (theMG->topLevel+1>=MAXLEVEL) return(NULL);
  if (algebraic) l = theMG->bottomLevel-1;
  else l = theMG->topLevel+1;

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
    DOWNGRID(theGrid) = GRID_ON_LEVEL(theMG,l-1);
    UPGRID(GRID_ON_LEVEL(theMG,l-1)) = theGrid;
    UPGRID(theGrid) = NULL;
  }
  else if (l==0)
  {
    DOWNGRID(theGrid) = NULL;
    UPGRID(theGrid) = NULL;
  }
  else
  {
    UPGRID(theGrid) = GRID_ON_LEVEL(theMG,l+1);
    DOWNGRID(theGrid) = NULL;
    DOWNGRID(GRID_ON_LEVEL(theMG,l+1)) = theGrid;
  }
  theGrid->mg = theMG;
  GRID_ON_LEVEL(theMG,l) = theGrid;
  if (algebraic) theMG->bottomLevel = l;
  else
  {
    theMG->topLevel = l;
    theMG->currentLevel = l;
  }

  return(theGrid);
}


/****************************************************************************/
/*D
   CreateNewLevelAMG - Create new amg level

   SYNOPSIS:
   GRID *CreateNewLevelAMG (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid structure

   DESCRIPTION:
   This function creates and initialized a new grid structure for bottomLevel - 1
   and returns a pointer to it.

   RETURN VALUE:
   GRID *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

GRID *CreateNewLevelAMG (MULTIGRID *theMG)
{
  GRID *theGrid;
  int l;

  if (theMG->bottomLevel-1<=-MAXLEVEL) return(NULL);

  l = theMG->bottomLevel-1;

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

  /* fill in further data */
  theGrid->mg = theMG;
  theGrid->level = l;
  theGrid->finer = theMG->grids[l+1];
  theMG->grids[l+1]->coarser = theGrid;

  theMG->grids[l] = theGrid;
  theMG->bottomLevel = l;

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


INT ClearMultiGridUsedFlags (MULTIGRID *theMG, INT FromLevel, INT ToLevel, INT mask)
{
  int i,level,elem,node,edge,vertex,vector,matrix;
  GRID *theGrid;
  ELEMENT *theElement;
  NODE *theNode;
  EDGE *theEdge;
  VECTOR *theVector;
  MATRIX *theMatrix;

  elem = mask & MG_ELEMUSED;
  node = mask & MG_NODEUSED;
  edge = mask & MG_EDGEUSED;
  vertex = mask & MG_VERTEXUSED;
  vector = mask & MG_VECTORUSED;
  matrix = mask & MG_MATRIXUSED;

  for (level=FromLevel; level<=ToLevel; level++)
  {
    theGrid = GRID_ON_LEVEL(theMG,level);
    if (elem || edge)
      for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      {
        if (elem) SETUSED(theElement,0);
        if (edge)
        {
          for (i=0; i<EDGES_OF_ELEM(theElement); i++)
          {
            theEdge = GetEdge(CORNER_OF_EDGE_PTR(theElement,i,0),
                              CORNER_OF_EDGE_PTR(theElement,i,1));
            SETUSED(theEdge,0);
          }
        }
      }
    if (node || vertex)
    {
      for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      {
        if (node) SETUSED(theNode,0);
        if (vertex) SETUSED(MYVERTEX(theNode),0);
      }
    }
    if (vector || matrix)
    {
      for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
      {
        if (vector) SETUSED(theVector,0);
        if (matrix)
          for (theMatrix=VSTART(theVector); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
            SETUSED(theMatrix,0);
      }
    }

  }

  return(GM_OK);
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
   MULTIGRID *CreateMultiGrid (char *MultigridName, char *BndValProblem,
   char *format, MEM heapSize, INT optimizedIE);

   PARAMETERS:
   .  MultigridName - name of multigrid
   .  domain - name of domain description from environment
   .  problem - name of problem description from environment
   .  format - name of format description from environment
   .  heapSize - size of heap to allocate for that multigrid in bytes
   .  optimizedIE - alloccate NodeElementList

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
                            char *format, MEM heapSize, INT optimizedIE)
{
  HEAP *theHeap,*theUserHeap;
  MULTIGRID *theMG;
  GRID *theGrid;
  INT i,ds;
  BVP *theBVP;
  BVP_DESC *theBVPDesc;
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
  if (BVP_SetBVPDesc(theBVP,&theMG->theBVPD))
  {
    PrintErrorMessage('E',"CreateMultiGrid","BVP not evaluated");
    return(NULL);
  }
  theBVPDesc = MG_BVPD(theMG);

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
  ds = FMT_S_MG(theFormat);
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
  MG_COARSE_FIXED(theMG) = 0;
  theMG->vertIdCounter = 0;
  theMG->nodeIdCounter = 0;
  theMG->elemIdCounter = 0;
  theMG->topLevel = -1;
  theMG->bottomLevel = 0;
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
  {
    GRID_ON_LEVEL(theMG,i) = NULL;
    GRID_ON_LEVEL(theMG,-i-1) = NULL;
  }

  if(optimizedIE == TRUE)
  {
    if ((MGNDELEMPTRARRAY(theMG)=
           GetTmpMem(theHeap,NDELEM_BLKS_MAX*sizeof(ELEMENT**)))==NULL)
    {
      ReleaseTmpMem(theHeap);
      PrintErrorMessage('E',"CreateMultiGrid",
                        "ERROR: could not allocate memory from the MGHeap");
      return (NULL);
    }
    for (i=0; i<NDELEM_BLKS_MAX; i++)
      MGNDELEMBLK(theMG,i) = NULL;
  }
  else
  {
    MGNDELEMPTRARRAY(theMG)=NULL;
  }

  /* allocate level 0 grid */
  theGrid = CreateNewLevel(theMG,0);
  if (theGrid==NULL)
  {
    DisposeMultiGrid(theMG);
    return(NULL);
  }

  /* allocate predefined mesh, e. g. corner vertices pointers */
        #ifdef ModelP
  if (me==master)
        #endif
  if (InsertMesh(theMG,&mesh))
  {
    DisposeMultiGrid(theMG);
    return(NULL);
  }
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

  /* reset pointer of midnode to edge */
  if (MIDNODE(theEdge) != NULL)
    SETNFATHER(MIDNODE(theEdge),NULL);

  /* dispose vector and its matrices from edge-vector if */
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
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
  GEOM_OBJECT *father;
  INT size;

  HEAPFAULT(theNode);

  /* call DisposeElement first! */
  assert(START(theNode) == NULL);
        #ifdef ModelP
  if (SONNODE(theNode) != NULL)
  {
    SETNFATHER(SONNODE(theNode),NULL);
  }
        #else
  assert(SONNODE(theNode) == NULL);
        #endif

  /* remove node from node list */
  GRID_UNLINK_NODE(theGrid,theNode);

  theVertex = MYVERTEX(theNode);
  father = NFATHER(theNode);
  if (father != NULL)
  {
    switch (NTYPE(theNode))
    {

    case (CORNER_NODE) :
      ASSERT(OBJT(father) == NDOBJ);
      SONNODE((NODE *)father) = NULL;
                                #ifdef TOPNODE
      if (theVertex != NULL)
        TOPNODE(theVertex) = (NODE *)father;
                                #endif
      break;

    case (MID_NODE) :
      ASSERT(OBJT(father) == EDOBJ);
      MIDNODE((EDGE *)father) = NULL;
      break;

    default :
      ASSERT(0);
      break;
    }
  }

  /* TODO delete old vertex handling */
  if (0)
    if (theVertex != NULL)
    {
                #ifdef ModelP
      /* vertices have to be linked and unlinked    */
      /* relative to the level they are created for */
      INT levelofvertex      = LEVEL(theVertex);
      MULTIGRID *MG           = MYMG(theGrid);
      GRID *GridOfVertex      = GRID_ON_LEVEL(MG,levelofvertex);

      if (SONNODE(theNode) == NULL)
        DisposeVertex(GridOfVertex,theVertex);
                #else
      DisposeVertex(theGrid,theVertex);
                #endif
    }

  if (NOOFNODE(theVertex)<1)
    RETURN(GM_ERROR);
  if (NOOFNODE(theVertex)==1)
    DisposeVertex(theGrid,theVertex);
  else
    DECNOOFNODE(theVertex);

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
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
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

  PRINTDEBUG(gm,1,(PFMT "DisposeVertex(): Gridlevel=%d theVertex=" VID_FMTX "\n",
                   me,GLEVEL(theGrid),VID_PRTX(theVertex)));

  theGrid = GRID_ON_LEVEL(MYMG(theGrid),LEVEL(theVertex));

  /* remove vertex from vertex list */
  GRID_UNLINK_VERTEX(theGrid,theVertex);

  if( OBJT(theVertex) == BVOBJ )
  {
    BNDP_Dispose(MGHEAP(MYMG(theGrid)),V_BNDP(theVertex));
    PutFreeObject(MYMG(theGrid),theVertex,sizeof(struct bvertex),BVOBJ);
  }
  else
    PutFreeObject(MYMG(theGrid),theVertex,sizeof(struct ivertex),IVOBJ);

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
  NODE    *theNode;
  VERTEX  *theVertex;
  EDGE    *theEdge;
  BNDS    *bnds;
  ELEMENT *theFather;
  ELEMENT *succe = SUCCE(theElement);
        #ifdef __THREEDIM__
  VECTOR  *theVector;
  DOUBLE *local,fac;
  INT k,m,o,l;
        #endif

  HEAPFAULT(theElement);

  GRID_UNLINK_ELEMENT(theGrid,theElement);

  theFather = EFATHER(theElement);

  if (LEVEL(theElement)>0)
  {
                #ifndef ModelP
    ASSERT(theFather != NULL);
                #endif

    /* check intergrid pointer from father */
    if (theFather != NULL)
    {
                        #ifdef ModelP
      int index = PRIO2INDEX(EPRIO(theElement));
                        #else
      int index = 0;
                        #endif
      ELEMENT *Next = NULL;
      ASSERT(index!=-1 && index<2);

      if (SON(theFather,index) == theElement)
      {
        if (succe != NULL)
        {
          if (EFATHER(succe)==theFather)
                                        #ifdef ModelP
            if (PRIO2INDEX(EPRIO(succe)) == PRIO2INDEX(EPRIO(theElement)))
                                        #endif
          {
            Next = succe;
          }
        }
        SET_SON(theFather,index,Next);
      }

      SETNSONS(theFather,NSONS(theFather)-1);

      PRINTDEBUG(gm,2,(PFMT "DisposeElement(): elem=" EID_FMTX
                       " father=" EID_FMTX " son0=%x son1=%x\n",
                       me,EID_PRTX(theElement),EID_PRTX(theFather),
                       SON(theFather,0),SON(theFather,1)));
    }
  }

        #ifdef ModelP
  /* reset father pointers of sons */
  /* TODO: possibly some son cannot be reached by GetAllSons, */
  /* because their father has not been on this proc and       */
  /* they lost their father pointers                          */
  if (NSONS(theElement)>0)
  {
    INT i,j,k,l,m,o;
    DOUBLE fac,*local;
    ELEMENT *SonList[MAX_SONS];

    if (GetAllSons(theElement,SonList)) RETURN(GM_FATAL);

    i = 0;
    while (SonList[i] != NULL)
    {
      PRINTDEBUG(gm,2,(PFMT "DisposeElement(): elem=" EID_FMTX
                       " deleting fatherpointer of son=" EID_FMTX "\n",
                       me,EID_PRTX(theElement),EID_PRTX(SonList[i])));
      SET_EFATHER(SonList[i],NULL);

                        #ifdef __THREEDIM__
      /* reset VFATHER of centernode vertex */
      for (j=0; j<CORNERS_OF_ELEM(SonList[i]); j++)
      {
        theNode = CORNER(SonList[i],j);
        if (NTYPE(theNode) != CENTER_NODE) continue;

        theVertex = MYVERTEX(theNode);
        if (VFATHER(theVertex) != NULL && VFATHER(theVertex) == theElement)
          VFATHER(theVertex) = NULL;
      }
                        #endif

      i++;
    }
  }
        #endif

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

        #ifdef __THREEDIM__
  /* reset VFATHER of sidenodes */
  for (j=0; j<SIDES_OF_ELEM(theElement); j++) {
    theNode = GetSideNode(theElement,j);
    if (theNode == NULL) continue;
    theVertex = MYVERTEX(theNode);
    if (VFATHER(MYVERTEX(theNode)) == theElement) {
      ELEMENT *theNb = NBELEM(theElement,j);

      VFATHER(theVertex) = theNb;
      if (theNb != NULL) {
        /* calculate new local coords */
        k = ONNBSIDE(theVertex);
        SETONSIDE(theVertex,k);
        m = CORNERS_OF_SIDE(theNb,k);
        local = LCVECT(theVertex);
        fac = 1.0 / m;
        V_DIM_CLEAR(local);
        for (o=0; o<m; o++) {
          l = CORNER_OF_SIDE(theNb,k,o);
          V_DIM_LINCOMB(1.0,local,1.0,
                        LOCAL_COORD_OF_ELEM(theNb,l),local);
        }
        V_DIM_SCALE(fac,local);
      }
    }
    SETONNBSIDE(theVertex,MAX_SIDES_OF_ELEM);
  }
        #endif

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

    /* edit VFATHER of midnode */
    if (MIDNODE(theEdge) != NULL)
    {
      theVertex = MYVERTEX(MIDNODE(theEdge));
      if (VFATHER(theVertex) == theElement) {
                #ifdef __TWODIM__
        theFather = NBELEM(theElement,j);
        VFATHER(theVertex) = theFather;
        if (theFather != NULL)
        {
          /* calculate new local coords */
          int co0,co1;

          /* reconstruct local coordinates of vertex */
          co0 = CORNER_OF_EDGE(theFather,j,0);
          co1 = CORNER_OF_EDGE(theFather,j,1);

          /* local coordinates have to be local towards pe */
          V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theFather,co0),
                        0.5, LOCAL_COORD_OF_ELEM(theFather,co1),
                        LCVECT(theVertex));
          SETONEDGE(theVertex,j);
        }
                            #endif
                #ifdef __THREEDIM__
        VFATHER(theVertex) = NULL;
                            #endif
      }
    }
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

        /* in serial case VFATHER pointer must always be present  */
        /* in parallel case only this MIDNODE pointer of          */
        /* edges have to be searched over EFATHER pointer,        */
        /* since VFATHER pointer can be NULL                      */
        /* overall MIDNODE pointer may be already NULL for ghosts */
#ifdef ModelP
        /* TODO: delete this */
        if (0)
          if (theFather != NULL)
          {
            INT i;

            for (i=0; i<EDGES_OF_ELEM(theFather); i++)
            {
              theEdge = GetEdge(CORNER(theFather,
                                       CORNER_OF_EDGE(theFather,i,0)),
                                CORNER(theFather,
                                       CORNER_OF_EDGE(theFather,i,1)));
              ASSERT(theEdge!=NULL);

              if (MIDNODE(theEdge) == theNode)
              {
                MIDNODE(theEdge) = NULL;
                break;
              }
            }
            ASSERT(i<EDGES_OF_ELEM(theFather) ||
                   EGHOST(theElement));
          }
#else
        theFather = VFATHER(theVertex);
        edge = ONEDGE(theVertex);
        theEdge = GetEdge(CORNER(theFather,
                                 CORNER_OF_EDGE(theFather,edge,0)),
                          CORNER(theFather,
                                 CORNER_OF_EDGE(theFather,edge,1)));
        ASSERT(theEdge!=NULL);
        MIDNODE(theEdge) = NULL;
#endif
      }
      DisposeNode(theGrid,theNode);
    }
  }

  /* dispose matrices from element-vector */
  if (dispose_connections)
    if (DisposeConnectionFromElement(theGrid,theElement))
      RETURN(1);

  /* reset neighbor pointers referencing element and dispose vectors in sides if */
  /* TODO: delete */
  if (0)
  {
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
  }

  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    ELEMENT *theNeighbor = NBELEM(theElement,i);

                #ifdef __THREEDIM__
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
    {
      theVector = SVECTOR(theElement,i);
      if (theVector!=NULL)
      {
        assert (VCOUNT(theVector) != 0);
        assert (VCOUNT(theVector) != 3);
        if (VCOUNT(theVector) == 1)
        {
          if (DisposeVector (theGrid,theVector))
            RETURN (1);
        }
        else
        {
                                        #ifdef ModelP
          /* TODO: sidevector case must be tested */
          assert(0);
                                        #endif
          if (!FindNeighborElement (theElement,i,&theNeighbor,&j))
            RETURN (1);
          VOBJECT(theVector) = (void*)theNeighbor;
          SETVECTORSIDE(theVector,j);
          SETVCOUNT(SVECTOR(theElement,i),1);
        }
      }
    }
                #endif
    if (theNeighbor!=NULL)
    {
      for (j=0; j<SIDES_OF_ELEM(theNeighbor); j++)
        if (NBELEM(theNeighbor,j)==theElement)
        {
          SET_NBELEM(theNeighbor,j,NULL);
          break;
        }
            #ifdef ModelP
      ASSERT(j<SIDES_OF_ELEM(theNeighbor) || EGHOST(theElement));
                        #else
      ASSERT(j<SIDES_OF_ELEM(theNeighbor));
            #endif
    }
  }

  /* dispose vector in center of element */
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
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

#ifndef ModelP
#define DO_NOT_DISPOSE  return (2)
#else
#define DO_NOT_DISPOSE  dispose=0
#endif

INT DisposeTopLevel (MULTIGRID *theMG)
{
  int l;
  GRID *theGrid;
        #ifdef ModelP
  int dispose = 1;
        #endif

  /* level 0 can not be deleted */
  l = theMG->topLevel;
  if (l<=0) DO_NOT_DISPOSE;
  if (theMG->bottomLevel<0) DO_NOT_DISPOSE;
  theGrid = GRID_ON_LEVEL(theMG,l);

  /* is level empty */
  if (FIRSTELEMENT(theGrid)!=NULL) DO_NOT_DISPOSE;
  if (FIRSTVERTEX(theGrid)!=NULL) DO_NOT_DISPOSE;
  if (FIRSTNODE(theGrid)!=NULL) DO_NOT_DISPOSE;

        #ifdef ModelP
  dispose = UG_GlobalMinINT(dispose);
  if (!dispose) return(2);
        #endif

  /* remove from grids array */
  GRID_ON_LEVEL(theMG,l) = NULL;
  GRID_ON_LEVEL(theMG,l-1)->finer = NULL;
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
   .  theGrid - grid to be removed

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

  theMG = MYMG(theGrid);

  if (GLEVEL(theGrid)<0)
    return (1);

  if (theGrid->finer != NULL)
    return(1);

  if (GLEVEL(theGrid)==0 && theMG->bottomLevel<0) return (1);

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

  /* level 0 can not be deleted */
  if (GLEVEL(theGrid) > 0)
    return(DisposeTopLevel(theMG));

  /* remove from grids array */
  GRID_ON_LEVEL(theMG,0) = NULL;
  theMG->currentLevel = theMG->topLevel = -1;
  theMG->nodeIdCounter = 0;
  theMG->vertIdCounter = 0;
  theMG->elemIdCounter = 0;

  PutFreeObject(theMG,theGrid,sizeof(GRID),GROBJ);

  return(0);
}

/****************************************************************************/
/*D
   DisposeAMGLevel - dispose bottom AMG level

   SYNOPSIS:
   INT DisposeAMGLevel (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid to remove from

   DESCRIPTION:
   This function removes the bottom AMG level from multigrid structure.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 no valid object number
   .n   2 no AMG levels
   D*/
/****************************************************************************/

INT DisposeAMGLevel (MULTIGRID *theMG)
{
  int l;
  GRID *theGrid;
  GRID *fineGrid;

  /* level 0 can not be deleted */
  l = theMG->bottomLevel;
  if (l>=0) return(2);
  theGrid = theMG->grids[l];
  fineGrid = theMG->grids[l+1];

  assert((FIRSTELEMENT(theGrid)==NULL)&&(FIRSTVERTEX(theGrid)==NULL)
         &&(FIRSTNODE(theGrid)==NULL));

  /* clear interpolation matrices from higher level!! */
  if (DisposeIMatricesInGrid(fineGrid))
    return(1);

  /* clear level */
  while (FIRSTVECTOR(theGrid)!=NULL)
    if (DisposeVector(theGrid,FIRSTVECTOR(theGrid)))
      return(1);

  /* remove from grids array */
  theMG->grids[l] = NULL;
  theMG->grids[l+1]->coarser = NULL;
  (theMG->bottomLevel)++;
  if (theMG->currentLevel<theMG->bottomLevel)
    theMG->currentLevel = theMG->bottomLevel;

  PutFreeObject(theMG,theGrid,sizeof(GRID),GROBJ);

  return(0);
}

/****************************************************************************/
/*D
   DisposeAMGLevels - dispose all AMG level

   SYNOPSIS:
   INT DisposeAMGLevels (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - multigrid to remove from

   DESCRIPTION:
   This function removes all AMG level from multigrid structure.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured
   D*/
/****************************************************************************/

INT DisposeAMGLevels (MULTIGRID *theMG)
{
  INT err;

  while ((err=DisposeAMGLevel(theMG))!=2)
    if (err==1)
    {
      PrintErrorMessage('E',"AMGTransferPreProcess","could not dispose AMG levels");
      REP_ERR_RETURN(1);
    }

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

/************************************************************************************/
/*                                                                                                                                                              */
/*    Enumeration of nodes from left to right corresponding to 0...n-1                          */
/*    ----------------------------------------------------------------				*/
/*                                                                                                                                                              */
/*                                         orphan nodes												*/
/*                                      |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|							*/
/*   N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N                      */
/*  |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|             */
/*                      ghost nodes						  master nodes						*/
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/*   Mapping should exist for all orphan nodes to vertices. 'foid' is id of first       */
/*   orphan node, 'non' number of orphan nodes. (MGIO_PAR) it looks (consequently): */
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/*                                                                      orphan nodes								*/
/*                                                      |~~~~~~~~~~~~~~~~~~~~|							*/
/*                                   N N N N N N N N N N N N N N N N N N N                      */
/*                                                          |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|            */
/*                                                                                        master nodes						*/
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/************************************************************************************/

static INT RenumberNodes (MULTIGRID *theMG, INT *foid, INT *non)
{
  INT i,nid;
  NODE *theNode;

  nid=0;
  if (procs==1)
  {
    /* ids for all nodes */
    for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,0)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      ID(theNode) = ID(MYVERTEX(theNode));
      nid = MAX(nid,ID(theNode));
    }
    nid++;
    non[0] = nid;
    for (i=1; i<=TOPLEVEL(theMG); i++)
      for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        ID(theNode) = nid++;
    foid[0] = 0;
  }
  else
  {
    /* ids for other nodes */
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (GHOST(theNode) && !USED(theNode))
          ID(theNode) = nid++;
    foid[0] = nid;

    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (GHOST(theNode) && USED(theNode))
          ID(theNode) = nid++;

    /* ids for master nodes */
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (MASTER(theNode) && USED(theNode))
          ID(theNode) = nid++;
    non[0] = nid-foid[0];

    /* ids for master nodes */
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (MASTER(theNode) && !USED(theNode))
          ID(theNode) = nid++;

  }

  return (0);
}

INT RenumberMultiGrid (MULTIGRID *theMG, INT *nboe, INT *nioe, INT *nbov, INT *niov, INT *foid, INT *non)
{
  NODE *theNode;
  ELEMENT *theElement;
  INT i,n_ioe,n_boe,vid,n_iov,n_bov,eid,lfoid,lnon,j;

  /* init used-flags */
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      SETUSED(theNode,0);
      SETUSED(MYVERTEX(theNode),0);
      SETTHEFLAG(MYVERTEX(theNode),0);
    }

  /* renumber elements and set orphan-flags for nodes and vertices */
  eid=n_ioe=n_boe=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
      if (EFATHER(theElement)==NULL || THEFLAG(theElement)==1)
      {
        ID(theElement) = eid++;
        if (OBJT(theElement)==BEOBJ) n_boe++;
        else n_ioe++;
        for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
        {
          SETUSED(CORNER(theElement,j),1);
          SETUSED(MYVERTEX(CORNER(theElement,j)),1);
        }
                                #ifdef ModelP
        assert(i==0 || EGHOST(theElement));
                                #endif
      }

  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
      if (EFATHER(theElement)!=NULL && THEFLAG(theElement)==0)
        ID(theElement) = eid++;
  if (nboe!=NULL) *nboe = n_boe;
  if (nioe!=NULL) *nioe = n_ioe;

  /* renumber vertices */
  vid=n_iov=n_bov=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
      if (THEFLAG(MYVERTEX(theNode))==0 && USED(MYVERTEX(theNode)) && OBJT(MYVERTEX(theNode))==BVOBJ)
      {
        ID(MYVERTEX(theNode)) = vid++;
        n_bov++;
        SETTHEFLAG(MYVERTEX(theNode),1);
      }
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
      if (THEFLAG(MYVERTEX(theNode))==0 && USED(MYVERTEX(theNode)) && OBJT(MYVERTEX(theNode))==IVOBJ)
      {
        ID(MYVERTEX(theNode)) = vid++;
        n_iov++;
        SETTHEFLAG(MYVERTEX(theNode),1);
      }
  for (i=0; i<=TOPLEVEL(theMG); i++)                                            /* not neccessary for i/o */
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
      if (THEFLAG(MYVERTEX(theNode))==0 && !USED(MYVERTEX(theNode)))
      {
        ID(MYVERTEX(theNode)) = vid++;
        SETTHEFLAG(MYVERTEX(theNode),1);
      }
  if (nbov!=NULL) *nbov = n_bov;
  if (niov!=NULL) *niov = n_iov;

  /* renumber nodes */
  if (RenumberNodes(theMG,&lfoid,&lnon)) return (1);
  if (foid!=NULL) *foid = lfoid;
  if (non!=NULL) *non = lnon;

  return (0);
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
  BVP_DESC *theBVPDesc;

  theMG   = MYMG(theGrid);
  entries = NN(theGrid);
  if (entries == 0) return (0);
  firstID = ID(FIRSTNODE(theGrid));
  theBVP = MG_BVP(theMG);
  theBVPDesc = MG_BVPD(theMG);

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
  ELEMENT *theElement,*After,*theFather;
  INT i,prio;

  if (cnt == 0) return (GM_OK);

  theFather = EFATHER(elemList[0]);
        #ifdef ModelP
  prio            = EPRIO(elemList[0]);
        #else
  prio            = 0;
        #endif

        #ifdef Debug
  /* check elements in list */
  for (i=0; i<cnt; i++)
  {
    theElement = elemList[i];

    ASSERT(theFather == EFATHER(theElement));
                #ifdef ModelP
    ASSERT(prio == EPRIO(theElement));
                #endif
  }
        #endif

  /* remove all elements from list */
  for (i=0; i<cnt; i++)
  {
    theElement = elemList[i];
    GRID_UNLINK_ELEMENT(theGrid,theElement);
  }

  /* and insert them at the end of the element list */
  After = NULL;
  for (i=0; i<cnt; i++)
  {
    GRID_LINKX_ELEMENT(theGrid,elemList[i],prio,After);
    After = elemList[i];
  }

  /* set son pointer of father to first son in list */
  if (EFATHER(elemList[0]) != NULL)
    SET_SON(EFATHER(elemList[0]),PRIO2INDEX(prio),elemList[0]);

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

NODE *InsertInnerNode (GRID *theGrid, DOUBLE *pos)
{
  VERTEX *theVertex;
  NODE *theNode;
  INT i;

  /* create objects */
  theVertex = CreateInnerVertex(theGrid);
  if (theVertex==NULL)
  {
    PrintErrorMessage('E',"InsertInnerNode","cannot create vertex");
    return(NULL);
  }
  theNode = CreateNode(theGrid,theVertex,NULL,LEVEL_0_NODE,0);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertInnerNode","cannot create node");
    return(NULL);
  }

  /* fill data */
  for (i=0; i<DIM; i++) CVECT(theVertex)[i] = pos[i];
  SETMOVE(theVertex,DIM);

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

NODE *InsertBoundaryNode (GRID *theGrid, BNDP *bndp)
{
  NODE *theNode;
  VERTEX *theVertex;
  INT move,part;

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

  if (BNDP_BndPDesc(bndp,&move,&part))
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  SETMOVE(theVertex,move);
  V_BNDP(theVertex) = bndp;

  theNode = CreateNode(theGrid,theVertex,NULL,LEVEL_0_NODE,0);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertBoundaryNode","cannot create node");
    return(NULL);
  }
        #ifdef TOPNODE
  TOPNODE(theVertex) = theNode;
        #endif

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

INT DeleteNode (GRID *theGrid, NODE *theNode)
{
  VERTEX *theVertex;
  ELEMENT *theElement;
  INT i;

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

INT DeleteNodeWithID (GRID *theGrid, INT id)
{
  NODE *theNode;

  /* find node */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    if (ID(theNode)==id) break;
  if (theNode==NULL)
  {
    PrintErrorMessage('E',"DeleteNodeWithID","node not found");
    RETURN(GM_ERROR);
  }
  return (DeleteNode(theGrid,theNode));
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

  if (OBJT(theElement)==BEOBJ && MOVED(theVertex)) return(theElement);

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
   RecreateBNDSofNode

   SYNOPSIS:
   static INT RecreateBNDSofNode(MULTIGRID *theMG, NODE *theNode)

   PARAMETERS:
   .  theMG - multigrid structure
   .  theNode - node with new BNDP

   DESCRIPTION:
   This function searches the boundary sides located at 'theNode' and recreate
   the corresponding BNDSs of these sides.
   It assumes that 'theNode' and the neighbour boundary nodes are in the same patch.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

static INT RecreateBNDSofNode (MULTIGRID *theMG, NODE *theNode)
{
  ELEMENT *theElement, *sonElem, *NBElem;
  ELEMENT *SonList[MAX_SONS];
  BNDS *bnds;
  BNDP *sidebndp[MAX_CORNERS_OF_SIDE];
  INT m,i,j,k,l;

  /* first scan father element of theNode */
  theElement = VFATHER(MYVERTEX(theNode));
  GetAllSons(theElement, SonList);
  for (i=0; i<NSONS(theElement); i++)
  {
    /* search side in son element with theNode as a corner */
    sonElem = SonList[i];
    if (OBJT(sonElem)!=BEOBJ) continue;
    for (j=0; j<SIDES_OF_ELEM(sonElem); j++)
      for (k=0; k<CORNERS_OF_SIDE(sonElem,j); k++)
        if (CORNER(sonElem,CORNER_OF_SIDE(sonElem,j,k))==theNode)
        {
          bnds = ELEM_BNDS(sonElem,j);
          if (bnds==NULL) continue;
          if (BNDS_Dispose(MGHEAP(theMG),bnds))
            return(GM_ERROR);
          /* create BNDS from BNDPs of this element side */
          for (l=0; l<CORNERS_OF_SIDE(sonElem,j); l++)
            sidebndp[l] = V_BNDP(MYVERTEX(CORNER(sonElem,CORNER_OF_SIDE(sonElem,j,l))));
          bnds = BNDP_CreateBndS(MGHEAP(theMG),sidebndp,CORNERS_OF_SIDE(sonElem,j));
          SET_BNDS(sonElem,j,bnds);
        }
  }

  if (NTYPE(theNode)==MID_NODE) return(GM_OK);

  /* scan all neighbour elements of the father element */
  for (m=0; m<SIDES_OF_ELEM(theElement); m++)
  {
    NBElem = NBELEM(theElement,m);
    if (NBElem==NULL) continue;
    if (OBJT(NBElem)!=BEOBJ) continue;
    GetAllSons(NBElem, SonList);
    for (i=0; i<NSONS(NBElem); i++)
    {
      /* search side in son element with theNode as a corner */
      sonElem = SonList[i];
      if (OBJT(sonElem)!=BEOBJ) continue;
      for (j=0; j<SIDES_OF_ELEM(sonElem); j++)
        for (k=0; k<CORNERS_OF_SIDE(sonElem,j); k++)
          if (CORNER(sonElem,CORNER_OF_SIDE(sonElem,j,k))==theNode)
          {
            bnds = ELEM_BNDS(sonElem,j);
            if (bnds==NULL) continue;
            if (BNDS_Dispose(MGHEAP(theMG),bnds))
              return(GM_ERROR);
            /* create BNDS from BNDPs of this element side */
            for (l=0; l<CORNERS_OF_SIDE(sonElem,j); l++)
              sidebndp[l] = V_BNDP(MYVERTEX(CORNER(sonElem,CORNER_OF_SIDE(sonElem,j,l))));
            bnds = BNDP_CreateBndS(MGHEAP(theMG),sidebndp,CORNERS_OF_SIDE(sonElem,j));
            SET_BNDS(sonElem,j,bnds);
          }
    }
  }
  return(GM_OK);
}

/****************************************************************************/
/*D
   MoveBndMidNode - set new position for a midnode on a boundary

   SYNOPSIS:
   INT MoveBndMidNode (MULTIGRID *theMG, VERTEX *theVertex);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theVertex - vertex to move

   DESCRIPTION:
   This function moves a given boundary vertex according to ist actual local
   coordinates. This function should only be called by MoveMidNode.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT MoveBndMidNode (MULTIGRID *theMG, VERTEX *theVertex)
{
  ELEMENT *theElement;
  NODE *Node0,*Node1,*sonNode, *theNode;
  EDGE *theEdge;
  BNDP *bndp;
  BNDS *bnds;
  DOUBLE *global,*local, diffmin, lambda_min;
  DOUBLE_VECTOR bnd_global, bnd_local, VecA, VecB, BndPoint;
  DOUBLE diff, lambda[DIM_OF_BND], *CornerPtrs[MAX_CORNERS_OF_ELEM], VecLen;
  INT co0,co1,edge,coe,n, ndiff;

  theElement = VFATHER(theVertex);
  edge = ONEDGE(theVertex);
  bnds = ELEM_BNDS(theElement,edge);
  if (bnds==NULL) return(GM_OK);

  co0 = CORNER_OF_EDGE(theElement,edge,0);
  co1 = CORNER_OF_EDGE(theElement,edge,1);
  theEdge=GetEdge(CORNER(theElement,co0),CORNER(theElement,co1));
  if (theEdge==NULL) return(GM_OK);    /* this should not happen */
  theNode = MIDNODE(theEdge);
  if (theNode==NULL) return(GM_OK);

  /* check actual local and global coordinates */
  local = LCVECT(theVertex);
  global = CVECT(theVertex);
  CORNER_COORDINATES(theElement,coe,CornerPtrs);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,global,bnd_local);
  if (V_DIM_ISEQUAL(bnd_local,local)) return(GM_OK);    /* nothing to do */

  /* calculate new lambda according to old local coordinates */
  Node0 = CORNER(theElement,co0);
  Node1 = CORNER(theElement,co1);
  /* global coordinates w.r.t. actual local coordinates */
  LOCAL_TO_GLOBAL(coe,CornerPtrs,local,bnd_global);

  /* find lambda with minimum distance to bnd_global */
  diffmin = 1E30;
  for (n=1; n<=100; n++)
  {
    lambda[0] = n/100.;
    BNDS_Global(bnds,lambda,BndPoint);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint,bnd_global,diff);
    if (diff < diffmin)
    {
      diffmin = diff;
      lambda_min = lambda[0];
      ndiff = n;
    }
  }
  /* more accurate */
  for (n=1; n<=100; n++)
  {
    lambda[0] = n/100./100. + ndiff/100.;
    BNDS_Global(bnds,lambda,BndPoint);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint,bnd_global,diff);
    if (diff < diffmin)
    {
      diffmin = diff;
      lambda_min = lambda[0];
    }
  }

  /* delete bndp ... */
  if (BNDP_Dispose(MGHEAP(theMG),V_BNDP(theVertex)))
    return(GM_ERROR);
  /* and create a new one */
  bndp = BNDP_CreateBndP(MGHEAP(theMG),V_BNDP(MYVERTEX(Node0)),
                         V_BNDP(MYVERTEX(Node1)),lambda_min);
  if (bndp == NULL)
    return(GM_ERROR);
  V_BNDP(theVertex) = bndp;
  /* new global coordinates */
  if (BNDP_Global(bndp,global))
    RETURN(GM_ERROR);

  /* global coordinates on a straight boundary */
  LOCAL_TO_GLOBAL(coe,CornerPtrs,local,bnd_global);

  V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
  if (diff > MAX_PAR_DIST) {
    SETMOVED(theVertex,1);
    UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,global,local);
  }
  /* update neighboring bnds */
  RecreateBNDSofNode(theMG,theNode);
  for (sonNode=SONNODE(theNode); sonNode!=0; sonNode=SONNODE(sonNode))
    RecreateBNDSofNode(theMG,sonNode);

  return(GM_OK);
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

INT MoveMidNode (MULTIGRID *theMG, NODE *theNode, DOUBLE lambda, INT update)
{
  ELEMENT *theElement;
  NODE *Node0,*Node1,*sonNode;
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
  if (NTYPE(theNode) != MID_NODE) {
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
      V_DIM_COPY(bnd_global,global);
      UG_GlobalToLocal(n,(const DOUBLE **)x,global,local);
    }
    RecreateBNDSofNode(theMG,theNode);
    for (sonNode=SONNODE(theNode); sonNode!=0; sonNode=SONNODE(sonNode))
      RecreateBNDSofNode(theMG,sonNode);
  }

  if (update==FALSE) return(GM_OK);
  /* Warning: O(n) Operation! */
  for(k=LEVEL(theNode)+1; k<=TOPLEVEL(theMG); k++)
    for (theVertex=FIRSTVERTEX(GRID_ON_LEVEL(theMG,k));
         theVertex!=NULL; theVertex=SUCCV(theVertex))
      if ((OBJT(theVertex) != BVOBJ)) {
        CORNER_COORDINATES(VFATHER(theVertex),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(theVertex),CVECT(theVertex));
      }
      else
      {
        if(MoveBndMidNode(theMG,theVertex)) return(GM_ERROR);
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

  if (NTYPE(theNode) != CENTER_NODE) {
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
  if (NTYPE(theNode) != SIDE_NODE) {
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
      V_DIM_COPY(bnd_global,global);
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

INT MoveNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *newPos, INT update)
{
  VERTEX *theVertex;
  EDGE *theEdge;
  ELEMENT *theElement;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE_VECTOR oldPos;
  INT n,k;

  /* set k (and theNode) to the level where the node
     appears the first time */
  while (CORNERTYPE(theNode))
    theNode = (NODE *)NFATHER(theNode);
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
      PrintErrorMessageF('W',"MoveNode",
                         "cannot find father element for Node %d",ID(theNode));
      V_DIM_COPY(oldPos,CVECT(theVertex));
      return(GM_ERROR);
    }
    else
    {
      CORNER_COORDINATES(theElement,n,x);
      UG_GlobalToLocal(n,(const DOUBLE **)x,newPos,LCVECT(theVertex));
      for (k=0; k<EDGES_OF_ELEM(theElement); k++) {
        theEdge =
          GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,k,0)),
                  CORNER(theElement,CORNER_OF_EDGE(theElement,k,1)));
        if (MIDNODE(theEdge) == theNode) {
          SETONEDGE(theVertex,k);
          break;
        }
      }
      VFATHER(theVertex) = theElement;
    }
  }
  else
    V_DIM_COPY(newPos,CVECT(theVertex));

  if (update==FALSE) return(GM_OK);
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
   MoveFreeBoundaryVertex - move a vertex on a free boundary

   SYNOPSIS:
   INT MoveFreeBoundaryVertex (MULTIGRID *theMG, VERTEX *vert, DOUBLE *newPos)

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  vert - vertex to move
   .  newPos - global coordinate for new position

   DESCRIPTION:
   This function moves a given vertex to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.
   To this end the function 'FinishMovingFreeBoundaryVertices' has to be
   called after having done the last call of 'MoveFreeBoundaryVertex'.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.

   SEE ALSO:
   FinishMovingFreeBoundaryVertices
   D*/
/****************************************************************************/

INT MoveFreeBoundaryVertex (MULTIGRID *theMG, VERTEX *vert, const DOUBLE *newPos)
{
  ELEMENT *theElement;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  INT n;

  if (OBJT(vert) != BVOBJ)
    REP_ERR_RETURN(GM_ERROR);
  if (MOVE(vert)!=DIM)
    REP_ERR_RETURN(GM_ERROR);

  if (BNDP_Move(V_BNDP(vert),newPos))
    REP_ERR_RETURN(GM_ERROR);
  V_DIM_COPY(newPos,CVECT(vert));
  if (LEVEL(vert) > 0)
  {
    theElement = VFATHER(vert);
    if (theElement == NULL)
      REP_ERR_RETURN(GM_ERROR)
      else
      {
        CORNER_COORDINATES(theElement,n,x);
        UG_GlobalToLocal(n,(const DOUBLE **)x,newPos,LCVECT(vert));
      }
  }

  return(GM_OK);
}

/****************************************************************************/
/*D
   FinishMovingFreeBoundaryVertices - finish moving of free boundary vertices

   SYNOPSIS:
   INT FinishMovingFreeBoundaryVertices (MULTIGRID *theMG)

   PARAMETERS:
   .  theMG - pointer to multigrid

   DESCRIPTION:
   The complete multigrid structure is moved hierachically such that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.

   SEE ALSO:
   MoveFreeBoundaryVertex
   D*/
/****************************************************************************/

INT FinishMovingFreeBoundaryVertices (MULTIGRID *theMG)
{
  VERTEX *vert;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  INT n,lev;

  for(lev=1; lev<=TOPLEVEL(theMG); lev++)
    for (vert=FIRSTVERTEX(GRID_ON_LEVEL(theMG,lev)); vert!=NULL; vert=SUCCV(vert))
      if ((OBJT(vert) != BVOBJ))
      {
        CORNER_COORDINATES(VFATHER(vert),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(vert),CVECT(vert));
      }
  return(GM_OK);
}

/****************************************************************************/
/*D
   GetMidNodeParam - Get local position of a midnode on an edge

   SYNOPSIS:
   INT GetMidNodeParam (NODE * theNode, DOUBLE *lambda)

   PARAMETERS:
   .  theNode - midnode
   .  lambda  - local coordinate of midnode w.r.t. the edge

   DESCRIPTION:
   This function gives the local coordinate of a midnode with respect to its edge
   (0 < lambda < 1). The function is called by SmoothGrid and can only be applied
   in 2D.

   RETURN VALUE:
   INT
   .n   GM_OK when ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT GetMidNodeParam (NODE * theNode, DOUBLE *lambda)
{

  ELEMENT *theElement;
  NODE *Node0,*Node1;
  VERTEX *theVertex;
  BNDS *bnds;
  DOUBLE *global;
  DOUBLE_VECTOR BndPoint0, BndPoint1;
  DOUBLE len0,len1, bndLambda[DIM_OF_BND], Lambda0, Lambda1, midLambda;
  INT co0,co1,edge,iter,MaxIter;

#ifdef __THREEDIM__
  PrintErrorMessage('E',"GetMidNodeParam","3D not implemented yet");
  return(GM_ERROR);
#endif

  if (!MIDTYPE(theNode)) {
    PrintErrorMessage('E',"GetMidNodeParam","node not a midnode");
    return(GM_ERROR);
  }
  theVertex = MYVERTEX(theNode);
  theElement = VFATHER(theVertex);

  edge = ONEDGE(theVertex);
  co0 = CORNER_OF_EDGE(theElement,edge,0);
  co1 = CORNER_OF_EDGE(theElement,edge,1);
  Node0 = CORNER(theElement,co0);
  Node1 = CORNER(theElement,co1);

  V_DIM_EUKLIDNORM_OF_DIFF(CVECT(theVertex),CVECT(MYVERTEX(Node0)),len0);
  V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(Node1)),CVECT(MYVERTEX(Node0)),len1);
  *lambda = len0/len1;

  if (OBJT(theVertex)!=BVOBJ) return(GM_OK);

  if (!MOVED(theVertex)) return(GM_OK);

  bnds = ELEM_BNDS(theElement,edge);
  Lambda0 = 0.;
  Lambda1 = 1.;
  iter=0;
  MaxIter = 40;

  global = CVECT(theVertex);
  do
  {
    iter++;
    midLambda = 0.5*(Lambda0+Lambda1);
    bndLambda[0] = Lambda0;
    BNDS_Global(bnds,bndLambda,BndPoint0);
    bndLambda[0]= midLambda;
    BNDS_Global(bnds,bndLambda,BndPoint1);
    V_DIM_EUKLIDNORM_OF_DIFF(global,BndPoint0,len0);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint1,BndPoint0,len1);
    if (len0<len1)
      Lambda1 = midLambda;
    else
      Lambda0 = midLambda;
  }
  while (!V_DIM_ISEQUAL(BndPoint0,global) && iter<MaxIter);
  *lambda = Lambda0;

  if (iter > (MaxIter-2))
    PrintErrorMessage('W',"GetMidNodeParam","could not determine lambda");

  return(GM_OK);
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

INT CheckOrientation (INT n, VERTEX **vertices)
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
INT CheckOrientation (INT n, VERTEX **vertices)
{
  DOUBLE_VECTOR diff[3],rot;
  DOUBLE det;
  INT i;

  /* TODO: this case */
  if (n == 8 || n==6)
    return(1);

  for (i=1; i<n; i++)
    V3_SUBTRACT(CVECT(vertices[i]),CVECT(vertices[0]),diff[i-1]);
  V3_VECTOR_PRODUCT(diff[0],diff[1],rot);
  V3_SCALAR_PRODUCT(rot,diff[2],det);

  if (det < 0.0)
    return(0);

  return(1);
}
#endif

INT CheckOrientationInGrid (GRID *theGrid)
{
  ELEMENT *theElement;
  NODE *theNode;
  VERTEX *vertices[MAX_CORNERS_OF_ELEM];
  INT i;

  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      if (theNode==NULL) return (1);
      vertices[i] = MYVERTEX(theNode);
      if (vertices[i]==NULL) return (1);
    }
    if (!CheckOrientation (CORNERS_OF_ELEM(theElement),vertices)) return (1);
  }

  return (0);
}

static INT NeighborSearch_O_n(INT n, NODE **Node, MULTIGRID *theMG, INT *NbrS, ELEMENT **Nbr)
{
  INT i,j,k,l,m,index,fnd,num,IndexOfDivPart,IndexOfModPart;
  NODE            *sideNode[MAX_CORNERS_OF_SIDE];
  ELEMENT         *interestingElements[ELEMS_OF_NODE_MAX];
  ELEMENT         *newinterestingElements[ELEMS_OF_NODE_MAX];
  ELEMENT         *theElement;
  NODE            *NeighborNode;

  /* for all sides of the element to be created */
  for (i=0; i<SIDES_OF_REF(n); i++)
  {
    for(j=0; j<CORNERS_OF_SIDE_REF(n,i); j++)            /* for the 3 corners*/
    {
      sideNode[j] = Node[CORNER_OF_SIDE_REF(n,i,j)];                  /* into sideNode*/
    }

    /*CAD*/
    j=0;

    /*fill interestingElements with elements of first node:*/
    m = 0;index=0;
    IndexOfDivPart = ID(sideNode[0]) / NO_NODES_OF_BLK;
    IndexOfModPart = ID(sideNode[0]) % NO_NODES_OF_BLK;
    index = MGNDELEMOFFS(IndexOfModPart,0);
    while(MGNDELEMBLKENTRY(theMG,IndexOfDivPart,index) != NULL)
    {
      interestingElements[m] = MGNDELEMBLKENTRY(theMG,IndexOfDivPart,index);
      m++; index++;
    }
    interestingElements[m] = NULL;


    j=1;             /*because the first node was already considered*/

    while( (j < CORNERS_OF_SIDE_REF(n,i) ) && (interestingElements[0] != NULL) )
    {
      k=0;
      l=0;
      IndexOfDivPart = ID(sideNode[j]) / NO_NODES_OF_BLK;
      IndexOfModPart = ID(sideNode[j]) % NO_NODES_OF_BLK;

      while (interestingElements[k] != NULL)
      {
        index = MGNDELEMOFFS(IndexOfModPart,0);

        while(MGNDELEMBLKENTRY(theMG,IndexOfDivPart,index) != NULL)
        {
          if(interestingElements[k] == MGNDELEMBLKENTRY(theMG,IndexOfDivPart,index))
          {
            newinterestingElements[l] = interestingElements[k];
            l++;
          }
          index++;
        }                         /* of while(MGNDELEMBLKENTRY(theMG,IndexOfDivPart,index) != NULL) */
        k++;
      }                  /* of while (interestingElements[k] != NULL)*/


      m = 0;
      newinterestingElements[l] = NULL;

      while (newinterestingElements[m] != NULL)
      {
        interestingElements[m] = newinterestingElements[m];
        m++;
      }
      interestingElements[m] = NULL;

      j++;                   /*next node of sideNode[j]*/

    }             /* while( (j < CORNERS_OF_SIDE_REF(n,i) ) && (interestingElements[0] != NULL) ) */

    theElement = interestingElements[0];

    /*CAD*/

    if (theElement != NULL)
    {
      /*for all sides of the already existing neighbour elememt*/
      fnd = 0;



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
            PrintErrorMessage('E',"InsertElement -> NeighborSearch_O_n",
                              "neighbor relation inconsistent");
            return(1);
          }
          Nbr[i] = theElement;
          NbrS[i] = j;
          fnd = 1;
        }

        if (fnd == 1)
          j = SIDES_OF_ELEM(theElement);               /*go on with next side*/

      }
    }
  }

  return(0);
  /*... O(n)InsertElement ...*/


} /*of static INT NeighborSearch_O_n()*/




static INT NeighborSearch_O_nn(INT n, NODE **Node, GRID *theGrid, INT *NghbrSide, ELEMENT **Nghbr)
{
  INT i,jj,k,m,num;
  NODE            *sideNode[MAX_CORNERS_OF_SIDE];
  ELEMENT         *theElement;
  NODE            *NeighborNode;

  /*O(n*n)InsertElement ...*/
  /* for all sides of the element to be created */
  for (i=0; i<SIDES_OF_REF(n); i++)
  {
    for(jj=0; jj<CORNERS_OF_SIDE_REF(n,i); jj++ )
      sideNode[jj] = Node[CORNER_OF_SIDE_REF(n,i,jj)];

    /* for all neighbouring elements allready inserted */
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement))
    {
      /* for all sides of the neighbour element */
      for (jj=0; jj<SIDES_OF_ELEM(theElement); jj++)
      {
        num = 0;
        /* for all corners of the side of the neighbour */
        for (m=0; m<CORNERS_OF_SIDE(theElement,jj); m++)
        {
          NeighborNode = CORNER(theElement,
                                CORNER_OF_SIDE(theElement,jj,m));
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
          if (NBELEM(theElement,jj)!=NULL)
          {
            PrintErrorMessage('E',"InsertElement -> NeighborSearch_O_nn",
                              "neighbor relation inconsistent");
            return(1);
          }
          Nghbr[i] = theElement;
          NghbrSide[i] = jj;
        }
      }
    }
  }
  /* ... O(n*n)InsertElement */


  return(0);
} /*of "static INT NeighborSearch_O_nn()"*/


static INT NdElPtrArray_evalIndexes(INT n, INT *cornerID, MULTIGRID *theMG, INT *MIndex, INT *MBlock, NODE **Node, GRID *theGrid, INT* NbrS, ELEMENT** Nbr)
{
  INT retval, j, IndexOfDivPart, IndexOfModPart, Index, merkeIndex, helpIndex;

  /*CAD*/
  /*evaluate indexes for "Insert in NodeElementMatrix"*/
  /*check if ELEMS_OF_NODE_MAX is too small */
  /* O(n)InsertElement  . . . */
  for(j=0; j<MAX_CORNERS_OF_ELEM; j++)
  {
    MIndex[j] = -1;
    MBlock[j] = -1;
    /* theElement is the new one !!! created in CreateElement !!*/
  }

  for(j=0; j<CORNERS_OF_REF(n); j++)      /* for the 4 corners*//*CORNERS_OF_ELEM_REF gibtsowas?*/
  {
    IndexOfDivPart = 2 *j;
    IndexOfModPart = IndexOfDivPart + 1;
    if (cornerID[IndexOfDivPart] < NDELEM_BLKS_MAX)
    {
      Index =   MGNDELEMOFFS(cornerID[IndexOfModPart],0);
      merkeIndex = Index;
      if(MGNDELEMBLK(theMG,cornerID[IndexOfDivPart]) != NULL)
      {
        /*only in this case  problems like Attention !!! Node has more than ELEMS_OF_NODE_MAX !!iare possible*/

        while(  (MGNDELEMBLKENTRY(theMG,cornerID[IndexOfDivPart],Index) != NULL) && (MGNDELEMPTRARRAY(theMG) != NULL) )
        {
          Index++;
          helpIndex = Index - merkeIndex;
          if (helpIndex == ELEMS_OF_NODE_MAX)
          {
            /* change to O(n*n)InsertElement ...*/
            /*IE_MEM_PROB  2*/
            MGNDELEMPTRARRAY(theMG) = NULL;                             /* set the O(n)InsertAlgFlag to zero*/
            UserWrite("Warning concerning InsertElement :\n");
            UserWrite("\n");
            UserWrite("The O(n) version of InsertElement cannot be continued \n");
            UserWrite("because some nodes have more than ELEMS_OF_NODE_MAX elements!\n");
            UserWrite("\n");
            UserWrite("Now the O(n*n) version will finish your job.\n");
            UserWrite("\n");
            UserWrite("Either you drink a cup of coffee and wait till the end of\n");
            UserWrite("the slow version or kill the job and do the following:\n");
            UserWrite("\n");
            UserWrite("1. increase ELEMS_OF_NODE_MAX in gm.h\n");
            UserWrite("2. remove ugm.o\n");
            UserWrite("3. ugmake\n");
            UserWrite("\n");
            UserWrite("\n");
            UserWrite("    .  .  .  .  .  .  .  .  .  .  .  .  please wait \n");
            /* O(n*n)InsertElement ...*/
            if ( (retval = NeighborSearch_O_nn(n, Node, theGrid, NbrS, Nbr)) == 1 )
            {
              PrintErrorMessage('E',"InsertElement -> NdElPtrArray_evalIndexes",
                                "neighbor relation inconsistent");
              return (1);
            }

          }

        }
      }
      MIndex[j] = Index;
      MBlock[j] = cornerID[IndexOfDivPart];
    }
                #ifdef Debug
    else
    {
      printf("FATAL1!!!\n");
      fflush(stdout);
    }
                #endif

    if (MGNDELEMPTRARRAY(theMG) == NULL)
      j = CORNERS_OF_REF(n);
  }
  /*CAD*/
  return(0);
} /* of static INT NdElPtrArray_evalIndexes() */


static INT NdElPtrArray_GetMemAndCheckIDs(INT n, MULTIGRID *theMG, INT *h_ID, NODE **c_Node, INT *c_ID, NODE **Node)
{
  INT i,j,maxi;
  INT IndexOfDivPart;
  INT IndexOfModPart;

  *h_ID = -1;

  for (i=0; i<CORNERS_OF_REF(n); i++)
  {
    IndexOfDivPart = 2 *i;
    IndexOfModPart = IndexOfDivPart + 1;


    c_Node[i] = Node[i];
    c_ID[IndexOfDivPart] = ID(c_Node[i]);
    c_ID[IndexOfModPart] = c_ID[IndexOfDivPart] % NO_NODES_OF_BLK;
    c_ID[IndexOfDivPart] = c_ID[IndexOfDivPart] / NO_NODES_OF_BLK;
    /*now you can find the cornea ID i in the two components c_ID[IndexOfDivPart] and c_ID[IndexOfModPart]
       where c_ID[IndexOfDivPart] includes the ID DIV NO_NODES_OF_BLK and c_ID[IndexOfModPart] the ID MOD NO_NODES_OF_BLK        */

    if (*h_ID < c_ID[IndexOfDivPart])
      *h_ID = c_ID[IndexOfDivPart];

    /*Speicher bereitstellen !!!*/
    if (c_ID[IndexOfDivPart] < NDELEM_BLKS_MAX)
    {
      if (MGNDELEMBLK(theMG,c_ID[IndexOfDivPart]) == NULL)                           /* wenn es noch gar keinen Speicherblock dafuer gibt ! */
      {                               /* Speicher bereitstellen */
        j=0;
        while (MGNDELEMBLK(theMG,j) != NULL)
        {
          j++;
        }

        maxi = NO_NODES_OF_BLK*ELEMS_OF_NODE_MAX*sizeof(ELEMENT*);

        while (j <= c_ID[IndexOfDivPart])
        {
          /*
                                  theMG		MGNDELEMBLK(theMG,j) = GetMem(theMG->theHeap,maxi,FROM_TOP);
           */
          /*IE_MEM_PROB*/
          /*
             MGNDELEMBLK(theMG,j) = malloc(maxi);
           */
          if ((MGNDELEMBLK(theMG,j)=GetTmpMem(MGHEAP(theMG),maxi))==NULL)
          {
            PrintErrorMessage('E',"InsertElement","  ==> NdElPtrArray_GetMemAndCheckIDs( ) ERROR: No memory for MGNDELEMBLK(theMG,j)");
            return(1);
          }
          /* Mustermemset(nodeflag_array,0,(statistik[0]+1)*sizeof(INT)); */
          memset(MGNDELEMBLK(theMG,j),0,maxi);
          j++;
        }
      }
    }

  }
  return(0);

} /* of static INT NdElPtrArray_GetMemAndCheckIDs()*/


static INT Neighbor_Direct_Insert(INT n, ELEMENT **ElemList, INT *NbgSdList, INT* NbrS, ELEMENT **Nbr)
{
  INT i;

  for (i=0; i<SIDES_OF_REF(n); i++)
    Nbr[i] = ElemList[i];
  if (NbgSdList != NULL)
    for (i=0; i<SIDES_OF_REF(n); i++)
      NbrS[i] = NbgSdList[i];

  return(0);
}

static INT NdElPtrArray_Update(INT *MIndex, INT *MBlock, ELEMENT *theElement, MULTIGRID *theMG)
{
  INT j;

  /*CAD*/
  /* "Insert in NodeElementMatrix"*/
  /* evaluated indexes now are used */
  for(j=0; j<MAX_CORNERS_OF_ELEM; j++)
  {
    if (MIndex[j] != -1)
      MGNDELEMBLKENTRY(theMG,MBlock[j],MIndex[j]) = theElement;
    /* theElement is the new one !!! created in CreateElement !!*/
  }

  return(0);
}

ELEMENT *InsertElement (GRID *theGrid, INT n, NODE **Node, ELEMENT **ElemList, INT *NbgSdList)
{
  MULTIGRID *theMG;
  INT i,j,k,m,rv,found,tag,ElementType;
  INT NeighborSide[MAX_SIDES_OF_ELEM];
  NODE             *sideNode[MAX_CORNERS_OF_SIDE];
  VERTEX           *Vertex[MAX_CORNERS_OF_ELEM],*sideVertex[MAX_CORNERS_OF_SIDE];
  ELEMENT          *theElement,*Neighbor[MAX_SIDES_OF_ELEM];
  BNDS         *bnds[MAX_SIDES_OF_ELEM];
  BNDP         *bndp[MAX_CORNERS_OF_ELEM];
  NODE            *cornerNode[MAX_CORNERS_OF_ELEM];
  INT cornerID[2*MAX_CORNERS_OF_ELEM];
  INT MIndex[MAX_CORNERS_OF_ELEM];
  INT MBlock[MAX_CORNERS_OF_ELEM];
  INT highestid;
        #ifdef __TWODIM__
  VERTEX           *theVertex;
  NODE             *theNode;
        #endif

  theMG = MYMG(theGrid);
  for (i=0; i<CORNERS_OF_REF(n); i++) cornerNode[i] = NULL;
  for (i=0; i<2*MAX_CORNERS_OF_ELEM; i++) cornerID [i] = 0;

  /* check parameters */
    #ifdef __TWODIM__
  switch (n)
  {
  case 3 :
    tag = TRIANGLE;
    break;
  case 4 :
    tag = QUADRILATERAL;
    break;
  default :
    PrintErrorMessage('E',"InsertElement","only triangles and quadrilaterals allowed in 2D");
    return(NULL);
  }
    #endif

        #ifdef __THREEDIM__
  switch (n)
  {
  case 4 :
    tag = TETRAHEDRON;
    break;
  case 5 :
    tag = PYRAMID;
    break;
  case 6 :
    tag = PRISM;
    break;
  case 8 :
    tag = HEXAHEDRON;
    break;
  default :
    PrintErrorMessage('E',"InsertElement","only tetrahedrons, pyramids and hexahedrons are allowed in the 3D coarse grid");
    return(NULL);
  }
    #endif

  /* init vertices */
  for (i=0; i<n; i++)
  {
    PRINTDEBUG(gm,1,("InsertElement(): node[%d]=" ID_FMTX "vertex[%d]=" VID_FMTX "\n",
                     i,ID_PRTX(Node[i]),i,VID_PRTX(MYVERTEX(Node[i]))))
    Vertex[i] = MYVERTEX(Node[i]);
  }

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
  if (!CheckOrientation (n,Vertex))
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

  /**********************************************************************/
  /* here begins the revised(03/97) part of InsertElement ...	        */
  /* documentation see ftp://ftp.ica3.uni-stuttgart.de/pub/text/dirk.   */
  /*                                                                    */
  if (ElemList == NULL)
  {
    /* find neighboring elements */
    if(MGNDELEMPTRARRAY(theMG) != NULL)
    {
      /* initialisations for using the node-element-pointer-array */
      if ( (rv = NdElPtrArray_GetMemAndCheckIDs(n, theMG, &highestid,
                                                cornerNode, cornerID,
                                                Node)) == 1 )
      {
        PrintErrorMessage('E',"InsertElement",
                          " ERROR by calling NdElPtrArray_GetMemAndCheckIDs()");
        return(NULL);
      }
      IFDEBUG(gm,2)
      {
        INT i;

        for (i=0; i<NDELEM_BLKS_MAX; i++)
        {
          printf("i=%d blk=%08x\n",i,MGNDELEMBLK(theMG,i));
          fflush(stdout);
        }
      }
      ENDDEBUG
      if ( (rv = NdElPtrArray_evalIndexes(n, cornerID, theMG, MIndex,
                                          MBlock, Node, theGrid,
                                          NeighborSide, Neighbor))  == 1)
      {
        PrintErrorMessage('E',"InsertElement",
                          " ERROR by calling NdElPtrArray_evalIndexes()");
        return(NULL);
      }
    }

    if ( (highestid < NDELEM_BLKS_MAX) && (MGNDELEMPTRARRAY(theMG) != NULL) )
    {
      /* using the fast O(n) algorithm */
      if( (rv = NeighborSearch_O_n(n, Node, theMG, NeighborSide, Neighbor)) == 1)
      {
        PrintErrorMessage('E',"InsertElement"," ERROR by calling NeighborSearch_O_n()");
        return(NULL);
      }
    }
    else
    /* using the slow O(n*n) algorithm */
    if ( (rv = NeighborSearch_O_nn(n, Node, theGrid, NeighborSide, Neighbor)) == 1 )
    {
      PrintErrorMessage('E',"InsertElement"," ERROR by calling NeighborSearch_O_nn()");
      return(NULL);
    }
  }
  else
  {
    /* use given neighboring elements */
    if ( (rv = Neighbor_Direct_Insert(n, ElemList, NbgSdList, NeighborSide, Neighbor)) == 1)
    {
      PrintErrorMessage('E',"InsertElement"," ERROR by calling NeighborSearch_O_nn()");
      return(NULL);
    }
  }

  /* create the new Element */
  theElement = CreateElement(theGrid,tag,ElementType,Node,NULL,0);
  if (theElement==NULL)
  {
    PrintErrorMessage('E',"InsertElement","cannot allocate element");
    return(NULL);
  }

  if (MGNDELEMPTRARRAY(theMG) != NULL && ElemList==NULL)
  {
    /* update of the node-element-pointer-array */
    if ( (rv = NdElPtrArray_Update(MIndex, MBlock, theElement, theMG))==1)
    {
      PrintErrorMessage('E',"InsertElement",
                        " ERROR by calling NdElPtrArray_Update()");
      return(NULL);
    }
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
      if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
        if (DisposeDoubledSideVector(theGrid,Neighbor[i],
                                     NeighborSide[i],theElement,i))
          return(NULL);
            #endif
    }
  }

  SET_EFATHER(theElement,NULL);
  SETECLASS(theElement,RED_CLASS);

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

ELEMENT *InsertElementFromIDs (GRID *theGrid, INT n, INT *idList)
{
  MULTIGRID *theMG;
  NODE *Node[MAX_CORNERS_OF_ELEM],*theNode;
  INT i,j,found;

  /* check level */
  theMG = MYMG(theGrid);
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertElementFromIDs","only a multigrid with exactly one level can be edited");
    return(NULL);
  }

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

  return (InsertElement(GRID_ON_LEVEL(theMG,0),n,Node,NULL,NULL));
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
  ELEMENT *theElement;
  NODE **NList,*Nodes[MAX_CORNERS_OF_ELEM],*ListNode;
  VERTEX **VList;
  INT i,k,n,nv,j,maxlevel,l,move,part;

  if (theMesh == NULL) return(GM_OK);
  if (theMesh->nElements == NULL)
  {
    assert(theMesh->VertexLevel==NULL);
    theGrid = GRID_ON_LEVEL(theMG,0);
    for (i=0; i<theMesh->nBndP; i++)
      if (InsertBoundaryNode(theGrid,theMesh->theBndPs[i]) == NULL)
        return(GM_ERROR);

    for (i=0; i<theMesh->nInnP; i++)
      if (InsertInnerNode(theGrid,theMesh->Position[i]) == NULL)
        return(GM_ERROR);
    return(GM_OK);
  }

  /* prepare */
  nv = theMesh->nBndP + theMesh->nInnP;
  VList = (VERTEX **) GetTmpMem(MGHEAP(theMG),nv*sizeof(VERTEX *));
  if (VList == NULL) return(GM_ERROR);
  NList = (NODE **) GetTmpMem(MGHEAP(theMG),nv*sizeof(NODE *));
  if (NList == NULL) return(GM_ERROR);
  for (j=0; j<nv; j++) NList[j] = NULL;

  maxlevel = 0;
  if (theMesh->VertexLevel!=NULL)
  {
    for (i=0; i<theMesh->nBndP; i++)
    {
      theGrid = GRID_ON_LEVEL(theMG,theMesh->VertexLevel[i]);
      VList[i] = CreateBoundaryVertex(theGrid);
      assert(VList[i]!=NULL);
      if (BNDP_Global(theMesh->theBndPs[i],CVECT(VList[i]))) assert(0);
      if (BNDP_BndPDesc(theMesh->theBndPs[i],&move,&part))
        return(NULL);
      SETMOVE(VList[i],move);
      V_BNDP(VList[i]) = theMesh->theBndPs[i];
      maxlevel = MAX(maxlevel,theMesh->VertexLevel[i]);
    }
    for (i=theMesh->nBndP; i<nv; i++)
    {
      theGrid = GRID_ON_LEVEL(theMG,theMesh->VertexLevel[i]);
      VList[i] = CreateInnerVertex(theGrid);
      V_DIM_COPY(theMesh->Position[i-theMesh->nBndP],CVECT(VList[i]));
      maxlevel = MAX(maxlevel,theMesh->VertexLevel[i]);
    }
  }
  else
  {
    theGrid = GRID_ON_LEVEL(theMG,0);
    for (i=0; i<theMesh->nBndP; i++)
    {
      VList[i] = CreateBoundaryVertex(theGrid);
      assert(VList[i]!=NULL);
      if (BNDP_Global(theMesh->theBndPs[i],CVECT(VList[i]))) assert(0);
      if (BNDP_BndPDesc(theMesh->theBndPs[i],&move,&part))
        return(NULL);
      SETMOVE(VList[i],move);
      V_BNDP(VList[i]) = theMesh->theBndPs[i];
    }
    for (i=theMesh->nBndP; i<nv; i++)
    {
      VList[i] = CreateInnerVertex(theGrid);
      V_DIM_COPY(theMesh->Position[i-theMesh->nBndP],CVECT(VList[i]));
    }
  }

  for (j=1; j<=theMesh->nSubDomains; j++)
    for (k=0; k<theMesh->nElements[j]; k++)
    {
      i = theMesh->ElementLevel[j][k];
      theGrid = GRID_ON_LEVEL(theMG,i);
      n = theMesh->Element_corners[j][k];
      for (l=0; l<n; l++)
      {
        ListNode = NList[theMesh->Element_corner_ids[j][k][l]];
        if (ListNode==NULL || LEVEL(ListNode)<i)
        {
          Nodes[l] = CreateNode(theGrid,VList[theMesh->Element_corner_ids[j][k][l]],NULL,LEVEL_0_NODE,0);
          if (Nodes[l]==NULL) assert(0);
          NList[theMesh->Element_corner_ids[j][k][l]] = Nodes[l];
          if (ListNode==NULL || LEVEL(ListNode)<i-1)
          {
            SETNFATHER(Nodes[l],NULL);
          }
          else
          {
            SETNFATHER(Nodes[l],(GEOM_OBJECT *)ListNode);
            SONNODE(ListNode) = Nodes[l];
          }
        }
        else
        {
          Nodes[l] = ListNode;
        }
      }
      theElement = InsertElement (theGrid,n,Nodes,NULL,NULL);
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

  /* check element */
  if (theElement==NULL) return(0);

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
  ELEMENT *theElement,*theFather,*Sons[MAX_SONS];
  INT i;

  if (GLEVEL(theGrid) == 0) {
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement))
      if (PointInElement(pos,theElement) == 1)
        return(theElement);
    return(NULL);
  }
  theFather = FindElementFromPosition(DOWNGRID(theGrid),pos);
  if (theFather == NULL) {
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement))
      if (PointInElement(pos,theElement) == 1)
        return(theElement);
    return(NULL);
  }
  if (GetSons(theFather,Sons)) {
    ASSERT(0);
    return(NULL);
  }
  for (i=0; i<NSONS(theFather); i++)
    if (PointInElement(pos,Sons[i]) == 1)
      return(Sons[i]);

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
  INT left,right,part;

  ASSERT(OBJT(t) == BEOBJ);
  ASSERT(SIDE_ON_BND(t,side));

  BNDS_BndSDesc(ELEM_BNDS(t,side),&left,&right,&part);

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
  else if (NSONS(nb) == 1) {
                #ifdef ModelP
    if (SON(nb,0) == NULL) return(nb);
                #endif
    nb = SON(nb,0);
    if (NSONS(nb) == 1) {
            #ifdef ModelP
      if (SON(nb,0) == NULL) return(nb);
                #endif
      nb = SON(nb,0);
    }
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
  BVP_DESC *theBVPDesc;

  /* get BVP description */
  theBVP = MG_BVP(theMG);
  theBVPDesc = MG_BVPD(theMG);

  c = isCurrent ? '*' : ' ';

  if (longformat)
    UserWriteF(" %c %-20.20s %-20.20s %10lu %10lu\n",c,ENVITEM_NAME(theMG),
               BVPD_NAME(theBVPDesc), HeapSize(theMG->theHeap),HeapUsed(theMG->theHeap));
  else
    UserWriteF(" %c %-20.20s\n",c,ENVITEM_NAME(theMG));

  return;
}

/****************************************************************************/
/*D
   MultiGridStatus - List information about refinement type distribution

   SYNOPSIS:
   INT MultiGridStatus (MULTIGRID *theMG)

   PARAMETERS:
   .  theMG - structure to list
   .  green - statistic about green elements
   .  load  - statistic about load balancing

   DESCRIPTION:
   This function lists information about multigrid's element types

   RETURN VALUE:
   INT
   D*/
/****************************************************************************/

INT MultiGridStatus (MULTIGRID *theMG, INT gridflag, INT greenflag, INT lbflag, INT verbose)
{
  INT i,j,sons,maxsons,heap,used,free_bytes;
  INT red, green, yellow;
  INT mg_red,mg_green,mg_yellow;
  INT mg_greenrulesons[MAXLEVEL+1][MAX_SONS+1],mg_greenrules[MAXLEVEL+1];
  INT markcount[MAXLEVEL+1],closuresides[MAXLEVEL+1];
  INT mg_red_size,mg_green_size,mg_yellow_size,mg_sum_size;
  FLOAT sum,sum_div_red,redplusgreen_div_red;
  FLOAT mg_sum,mg_sum_div_red,mg_redplusgreen_div_red;
  FLOAT mg_sum_size_div_red,mg_redplusgreen_size_div_red;
  ELEMENT *theElement;
  GRID    *theGrid;
        #ifdef ModelP
  INT             *infobuffer;
  INT             **lbinfo;
  INT total_elements,sum_elements;
  INT master_elements,hghost_elements,vghost_elements,vhghost_elements;
  VChannelPtr mych;
  void *ptr;
        #endif

  mg_red = mg_green = mg_yellow = mg_sum = 0;
  mg_sum_div_red = mg_redplusgreen_div_red = 0.0;

  for (i=0; i<MAXLEVEL+1; i++)
  {
    /* for greenflag */
    mg_greenrules[i] = 0;
    for (j=0; j<MAX_SONS+1; j++) mg_greenrulesons[i][j] = 0;
    maxsons = 0;

    /* for gridflag */
    markcount[i] = 0;
    closuresides[i] = 0;
  }

        #ifdef ModelP
  infobuffer      = (INT *) malloc((procs+1)*(MAXLEVEL+1)*ELEMENT_PRIOS*sizeof(INT));
  if (infobuffer == NULL) assert(0);

  lbinfo          = (INT **) malloc((procs+1)*sizeof(INT));
  if (lbinfo == NULL) assert(0);

  memset((void *)infobuffer,0,(procs+1)*(MAXLEVEL+1)*ELEMENT_PRIOS*sizeof(INT));
  for (i=0; i<procs+1; i++)
    lbinfo[i] = infobuffer+(i*(MAXLEVEL+1)*ELEMENT_PRIOS);
  total_elements = sum_elements = 0;
  master_elements = hghost_elements = vghost_elements = vhghost_elements = 0;
        #endif

  if (verbose && gridflag)
  {
    UserWriteF("\nMULTIGRID STATISTICS:\n");
    UserWriteF("LEVEL      RED     GREEN    YELLOW        SUM     SUM/RED (RED+GREEN)/RED\n");
  }

  /* compute multi grid infos */
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    red = green = yellow = 0;
    sum = sum_div_red = redplusgreen_div_red = 0.0;

    for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement))
    {
      SETUSED(theElement,0);
      /* count eclasses */
      switch (ECLASS(theElement))
      {
      case RED_CLASS :         red++;          break;
      case GREEN_CLASS :       green++;        break;
      case YELLOW_CLASS :      yellow++;       break;
      default :                        assert(0);
      }
      /* count marks and closuresides */
      if (EstimateHere(theElement))
      {
        ELEMENT *MarkElement = ELEMENT_TO_MARK(theElement);
        INT marktype = GetRefinementMarkType(theElement);

        if (marktype==1 &&
            USED(MarkElement)==0)
        {
          markcount[LEVEL(MarkElement)]++;
          markcount[MAXLEVEL]++;
          for (j=0; j<SIDES_OF_ELEM(MarkElement); j++)
          {
            ELEMENT *NbElement = NBELEM(MarkElement,j);
            if (NbElement!=NULL && MARKCLASS(NbElement)==RED_CLASS)
            {
              closuresides[LEVEL(MarkElement)]++;
              closuresides[MAXLEVEL]++;
            }
          }
          SETUSED(MarkElement,1);
        }
      }
      /* green refinement statistics */
      switch (REFINECLASS(theElement))
      {
      case GREEN_CLASS :
        sons = NSONS(theElement);
        mg_greenrulesons[i][sons]++;
        mg_greenrulesons[i][MAX_SONS]+=sons;
        mg_greenrules[i]++;
        mg_greenrulesons[MAXLEVEL][sons]++;
        mg_greenrulesons[MAXLEVEL][MAX_SONS]+=sons;
        mg_greenrules[MAXLEVEL]++;
        if (maxsons < sons) maxsons = sons;
        break;
      default :
        break;
      }
                        #ifdef ModelP
      /* count master, hghost, vghost and vhghost elements */
      switch (EPRIO(theElement))
      {
      case PrioMaster :
        lbinfo[me][ELEMENT_PRIOS*i]++;
        lbinfo[me][ELEMENT_PRIOS*MAXLEVEL]++;
        break;
      case PrioHGhost :
        lbinfo[me][ELEMENT_PRIOS*i+1]++;
        lbinfo[me][ELEMENT_PRIOS*MAXLEVEL+1]++;
        break;
      case PrioVGhost :
        lbinfo[me][ELEMENT_PRIOS*i+2]++;
        lbinfo[me][ELEMENT_PRIOS*MAXLEVEL+2]++;
        break;
      case PrioVHGhost :
        lbinfo[me][ELEMENT_PRIOS*i+3]++;
        lbinfo[me][ELEMENT_PRIOS*MAXLEVEL+3]++;
        break;
      default :
        assert(0);
        break;
      }
                        #endif
    }
    sum = red + green + yellow;
    if (red > 0)
    {
      sum_div_red = sum / red;
      redplusgreen_div_red = ((float)(red+green)) / red;
    }
    else
    {
      sum_div_red = 0.0;
      redplusgreen_div_red = 0.0;
    }

    if (verbose && gridflag)
      UserWriteF("   %2d  %9d %9d %9d  %9.0f    %2.3f      %2.3f\n",
                 i,red,green,yellow,sum,sum_div_red,redplusgreen_div_red);

    mg_red          += red;
    mg_green        += green;
    mg_yellow       += yellow;
    mg_sum          += sum;
  }
  if (mg_red > 0)
  {
    mg_sum_div_red                  = mg_sum / mg_red;
    mg_redplusgreen_div_red = ((float)(mg_red + mg_green)) / mg_red;
  }
  else
  {
    mg_sum_div_red                  = 0.0;
    mg_redplusgreen_div_red = 0.0;
  }

  if (verbose && gridflag)
    UserWriteF("  ALL  %9d %9d %9d  %9.0f    %2.3f      %2.3f\n",
               mg_red,mg_green,mg_yellow,mg_sum,mg_sum_div_red,mg_redplusgreen_div_red);

  /* compute heap info */
  if (gridflag)
  {
    heap = HeapFreelistUsed(MGHEAP(theMG));
    used = HeapUsed(MGHEAP(theMG))-heap;
    free_bytes = (HeapSize(MGHEAP(theMG))-used)>>10;
    mg_sum_size = used>>10;
    if (mg_sum > 0)
    {
      mg_red_size = mg_sum_size*mg_red/mg_sum;
      mg_green_size = mg_sum_size*mg_green/mg_sum;
      mg_yellow_size = (float)mg_sum_size*mg_yellow/mg_sum;
    }
    else
    {
      mg_red_size = 0.0;
      mg_green_size = 0.0;
      mg_yellow_size = 0.0;
    }
    if (mg_red > 0)
    {
      mg_sum_size_div_red = ((float)mg_sum_size)/mg_red;
      mg_redplusgreen_size_div_red = ((float)(mg_red_size+mg_green_size))/mg_red;
    }
    else
    {
      mg_sum_size_div_red = 0.0;
      mg_redplusgreen_size_div_red = 0.0;
    }
  }

  /* set heap info in refine info */
  if (gridflag)
  {
    float new;
    float newpergreen;
    float predmax;

    SETMARKCOUNT(REFINEINFO(theMG),markcount[MAXLEVEL]);

    new = markcount[MAXLEVEL]*(2<<(DIM-1))*mg_sum_div_red;
    SETPREDNEW0(REFINEINFO(theMG),new);

    if (mg_greenrules[MAXLEVEL] > 0)
      newpergreen = ((float)mg_greenrulesons[MAXLEVEL][MAX_SONS])/mg_greenrules[MAXLEVEL];
    else
      newpergreen = 0;
    new = markcount[MAXLEVEL]*(2<<(DIM-1))+newpergreen*closuresides[MAXLEVEL];
    SETPREDNEW1(REFINEINFO(theMG),new);

    SETREAL(REFINEINFO(theMG),mg_sum);
    if (mg_sum_size_div_red > 0.0)
      predmax = free_bytes/mg_sum_size_div_red;
    else
      predmax = free_bytes/
#ifdef __TWODIM__
                sizeof(struct quadrilateral);
#else
                sizeof(struct hexahedron);
#endif
    SETPREDMAX(REFINEINFO(theMG),predmax);
  }

  /* list heap info */
  if (verbose && gridflag)
  {
    float predmax0,predmax1;

    UserWriteF(" HEAP  %7dKB %7dKB %7dKB  %7dKB    %2.3fKB    %2.3fKB\n",
               mg_red_size,mg_green_size,mg_yellow_size,mg_sum_size,
               mg_sum_size_div_red,mg_redplusgreen_size_div_red);

    if (mg_sum_size_div_red > 0.0)
      predmax0 = free_bytes/mg_sum_size_div_red;
    else
      predmax0 = free_bytes/
#ifdef __TWODIMM__
                 sizeof(struct quadrilateral);
#else
                 sizeof(struct hexahedron);
#endif
    if (mg_redplusgreen_size_div_red > 0.0)
      predmax1 = free_bytes/mg_redplusgreen_size_div_red;
    else
      predmax1 = free_bytes/
#ifdef __TWODIM__
                 sizeof(struct quadrilateral);
#else
                 sizeof(struct hexahedron);
#endif

    UserWriteF(" EST  FREE=%7dKB  MAXNEWELEMENTS free_bytes/(SUM/RED)=%9.0f "
               "FREE/((RED+GREEN)/RED)=%9.0f\n",
               free_bytes,predmax0,predmax1);
    UserWriteF(" EST %2d  ELEMS=%9.0f MARKCOUNT=%9.0f PRED_NEW0=%9.0f PRED_NEW1=%9.0f PRED_MAX=%9.0f\n",
               REFINESTEP(REFINEINFO(theMG)),REAL(REFINEINFO(theMG)),MARKCOUNT(REFINEINFO(theMG)),
               PREDNEW0(REFINEINFO(theMG)),PREDNEW1(REFINEINFO(theMG)),PREDMAX(REFINEINFO(theMG)));
    UserWriteF(" EST TRACE step=%d\n",refine_info.step);
    for (i=0; i<refine_info.step; i++)
      UserWriteF(" EST  %2d  ELEMS=%9.0f MARKS=%9.0f REAL=%9.0f PRED0=%9.0f PRED1=%9.0f PRED_MAX=%9.0f\n",
                 i,refine_info.real[i],refine_info.markcount[i],
                 ((i<refine_info.step) ? refine_info.real[i+1]-refine_info.real[i] : 0),
                 refine_info.predicted_new[i][0],
                 refine_info.predicted_new[i][1],refine_info.predicted_max[i]);
  }

  /* compute and list green rule info */
  if (verbose && greenflag)
  {
    UserWriteF("\nGREEN RULE STATISTICS:\n");
    UserWriteF("  LEVEL GREENSONS     RULES GREENSONS/RUL");
    for (j=0; j<8 && j<maxsons; j++) UserWriteF("  %1d/%2d/...",j,j+8);
    UserWriteF("\n");

    for (i=0; i<=TOPLEVEL(theMG); i++)
    {
      UserWriteF("     %2d %9d %9d         %2.3f",i,mg_greenrulesons[i][MAX_SONS],
                 mg_greenrules[i],
                 (mg_greenrules[i]!=0) ? ((float)mg_greenrulesons[i][MAX_SONS])/mg_greenrules[i] : 0);
      for (j=0; j<maxsons; j++)
      {
        UserWriteF(" %9d",mg_greenrulesons[i][j]);
        if ((j+1)%8 == 0) UserWriteF("\n%41s"," ");
      }
      UserWriteF("\n");
    }
    UserWriteF("    ALL %9d %9d         %2.3f",mg_greenrulesons[MAXLEVEL][MAX_SONS],
               mg_greenrules[MAXLEVEL],
               (mg_greenrules[MAXLEVEL]) ? ((float)mg_greenrulesons[MAXLEVEL][MAX_SONS])/
               mg_greenrules[MAXLEVEL] : 0.0);
    for (j=0; j<maxsons; j++)
    {
      UserWriteF(" %9d",mg_greenrulesons[MAXLEVEL][j]);
      if ((j+1)%8 == 0) UserWriteF("\n%41s"," ");
    }
    UserWriteF("\n");
  }

        #ifdef ModelP
  /* compute and list loadbalancing info */
  if (verbose && lbflag)
  {
    UserWriteF("\nLB INFO:\n");
    /* now collect lb info on master */
    if (me == master)
      for (i=1; i<procs; i++)
      {
        mych =ConnSync(i,3917);
        RecvSync(mych,(void *)lbinfo[i],(MAXLEVEL+1)*ELEMENT_PRIOS*sizeof(INT));
        DiscSync(mych);
      }
    else
    {
      mych = ConnSync(master,3917);
      SendSync(mych,(void *)lbinfo[me],(MAXLEVEL+1)*ELEMENT_PRIOS*sizeof(INT));
      DiscSync(mych);
      return(GM_OK);
    }

    /* sum levels over procs */
    for (i=0; i<procs; i++)
    {
      for (j=0; j<TOPLEVEL(theMG)+1; j++)
      {
        lbinfo[procs][ELEMENT_PRIOS*j] += lbinfo[i][ELEMENT_PRIOS*j];
        lbinfo[procs][ELEMENT_PRIOS*j+1] += lbinfo[i][ELEMENT_PRIOS*j+1];
        lbinfo[procs][ELEMENT_PRIOS*j+2] += lbinfo[i][ELEMENT_PRIOS*j+2];
        lbinfo[procs][ELEMENT_PRIOS*j+3] += lbinfo[i][ELEMENT_PRIOS*j+3];
      }
    }

    /* only master */
    if (lbflag >= 3)
    {
      UserWriteF(" LEVEL");
      for (i=0; i<ELEMENT_PRIOS*(TOPLEVEL(theMG)+1); i++)
      {
        UserWriteF(" %9d",i/ELEMENT_PRIOS);
      }
      UserWrite("\n");
      UserWriteF("PROC  ");
      for (i=0; i<ELEMENT_PRIOS*(TOPLEVEL(theMG)+1); i++)
      {
        UserWriteF(" %9s",(i%ELEMENT_PRIOS==0) ? "MASTER" : (i%ELEMENT_PRIOS==1) ? "HGHOST" :
                   (i%ELEMENT_PRIOS==2) ? "VGHOST" : "VHGHOST");
      }
      UserWrite("\n");
      for (i=0; i<procs; i++)
      {
        UserWriteF("%4d  ",i);
        for (j=0; j<ELEMENT_PRIOS*(TOPLEVEL(theMG)+1); j++)
        {
          UserWriteF(" %9d",lbinfo[i][j]);
        }
        UserWrite("\n");
      }
      UserWriteF("\n");
    }

    if (lbflag >= 2)
    {
      float memeff;

      UserWriteF("%5s %9s %9s %9s %9s %9s %6s\n",
                 "LEVEL","SUM","MASTER","HGHOST","VGHOST","VHGHOST","MEMEFF");
      for (i=0; i<=TOPLEVEL(theMG); i++)
      {
        sum_elements = lbinfo[procs][ELEMENT_PRIOS*i]+lbinfo[procs][ELEMENT_PRIOS*i+1]+
                       lbinfo[procs][ELEMENT_PRIOS*i+2]+lbinfo[procs][ELEMENT_PRIOS*i+3];
        if (sum_elements > 0)
          memeff = ((float)lbinfo[procs][ELEMENT_PRIOS*i])/sum_elements*100;
        else
          memeff = 0.0;
        UserWriteF("%4d %9d %9d %9d %9d %9d  %3.2f\n",i,sum_elements,
                   lbinfo[procs][ELEMENT_PRIOS*i],lbinfo[procs][ELEMENT_PRIOS*i+1],
                   lbinfo[procs][ELEMENT_PRIOS*i+2],lbinfo[procs][ELEMENT_PRIOS*i+3],memeff);
      }
      UserWrite("\n");

      UserWriteF("%4s %9s %9s %9s %9s %9s %6s\n",
                 "PROC","SUM","MASTER","HGHOST","VGHOST","VHGHOST","MEMEFF");
      for (i=0; i<procs; i++)
      {
        sum_elements = lbinfo[i][ELEMENT_PRIOS*MAXLEVEL]+lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+1]+
                       lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+2]+lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+3];
        if (sum_elements > 0)
          memeff = ((float)lbinfo[i][ELEMENT_PRIOS*MAXLEVEL])/sum_elements*100;
        else
          memeff = 0.0;
        UserWriteF("%4d %9d %9d %9d %9d %9d  %3.2f\n",i,sum_elements,
                   lbinfo[i][ELEMENT_PRIOS*MAXLEVEL],lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+1],
                   lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+2],lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+3],memeff);
      }
      UserWrite("\n");
    }

    if (lbflag >= 1)
    {
      float memeff;

      for (i=0; i<procs; i++)
      {
        master_elements += lbinfo[i][ELEMENT_PRIOS*MAXLEVEL];
        hghost_elements += lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+1];
        vghost_elements += lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+2];
        vhghost_elements += lbinfo[i][ELEMENT_PRIOS*MAXLEVEL+3];
      }
      total_elements = master_elements + hghost_elements + vghost_elements;
      if (total_elements > 0)
        memeff = ((float)master_elements)/total_elements*100;
      else
        memeff = 0.0;
      UserWriteF("%9s %9s %9s %9s %9s %6s\n","TOTAL","MASTER","HGHOST","VGHOST","VHGHOST","MEMEFF");
      UserWriteF("%9d %9d %9d %9d %9d  %3.2f\n",total_elements,master_elements,hghost_elements,
                 vghost_elements,vhghost_elements,memeff);
    }

  }
  free(infobuffer);
  free(lbinfo);
        #endif

  return (GM_OK);
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

  UserWrite("level maxlevel    #vert    #node    #edge    #elem    #side    #vect    #conn");
#ifdef __INTERPOLATION_MATRIX__
  UserWrite("    #imat");
#endif
  UserWrite("  minedge  maxedge\n");

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
    for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement))
      if (OBJT(theElement) == BEOBJ)
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
          if (SIDE_ON_BND(theElement,i))
            ns++;

#ifdef __INTERPOLATION_MATRIX__
    UserWriteF("%c %3d %8d %8ld %8ld %8ld %8ld %8ld %8ld %8ld %8ld %9.3e %9.3e\n",c,l,(int)TOPLEVEL(theMG),
               (long)NV(theGrid),(long)NN(theGrid),(long)NE(theGrid),(long)NT(theGrid),
               (long)ns,(long)NVEC(theGrid),(long)NC(theGrid),(long)NIMAT(theGrid),(float)hmin,(float)hmax);
#else
    UserWriteF("%c %3d %8d %8ld %8ld %8ld %8ld %8ld %8ld %8ld %9.3e %9.3e\n",c,l,(int)TOPLEVEL(theMG),
               (long)NV(theGrid),(long)NN(theGrid),(long)NE(theGrid),(long)NT(theGrid),
               (long)ns,(long)NVEC(theGrid),(long)NC(theGrid),(float)hmin,(float)hmax);
#endif
  }

#ifdef __INTERPOLATION_MATRIX__
  if (BOTTOMLEVEL(theMG)<0)
  {
    UserWrite("AMG levels:\n");
    for (l=-1; l>=BOTTOMLEVEL(theMG); l--)
    {
      theGrid = GRID_ON_LEVEL(theMG,l);

      c = (l==cl) ? '*' : ' ';

      UserWriteF("%c %3d %8d %8ld %8ld %8ld %8ld %8ld %8ld %8ld %8ld\n",c,l,(int)TOPLEVEL(theMG),
                 (long)NV(theGrid),(long)NN(theGrid),(long)NE(theGrid),(long)NT(theGrid),
                 (long)0,(long)NVEC(theGrid),(long)NC(theGrid),(long)NIMAT(theGrid));
    }
  }
#endif

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

    #ifdef ModelP
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
    /* count vectors and connections */
    for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
      if ((l==cl) || (VNCLASS(vec)<1))
        if (PRIO(vec) == PrioMaster) nvec++;

    /* count other objects */
    for (theElement=FIRSTELEMENT(theGrid);
         theElement!=NULL; theElement=SUCCE(theElement))
      if (EstimateHere(theElement))
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
            if (PRIO(theNode) == PrioMaster) nn++;
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
            V_DIM_EUKLIDNORM_OF_DIFF(CVECT(v0),CVECT(v1),h);
            hmin = MIN(hmin,h);
            hmax = MAX(hmax,h);
          }
        }
      }
  }
  nn = UG_GlobalSumINT(nn);
  ne = UG_GlobalSumINT(ne);
  nt = UG_GlobalSumINT(nt);
  ns = UG_GlobalSumINT(ns);
  nvec = UG_GlobalSumINT(nvec);
  nc = UG_GlobalSumINT(nc);
  hmin = UG_GlobalMinDOUBLE(hmin);
  hmax = UG_GlobalMaxDOUBLE(hmax);
  UserWrite("\nsurface of all processors up to current level:\n");
  UserWriteF("%c %3d %8d %8s %8ld %8s %8ld %8ld %8ld %8s %9.3e %9.3e\n",
             ' ',minl,(int)cl,
             "---",(long)nn,"        ",(long)nt,
             (long)ns,(long)nvec,"        ",(float)hmin,(float)hmax);
        #endif

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
  INT i,part;

  theVertex = MYVERTEX(theNode);

  /******************************/
  /* print standard information */
  /******************************/
  /* line 1 */ UserWriteF("NODEID=" ID_FFMTE " CTRL=%8lx VEID="
                          VID_FMTX " LEVEL=%2d",
                          ID_PRTE(theNode),(long)CTRL(theNode),
                          VID_PRTX(theVertex),LEVEL(theNode));

  /* print coordinates of that node */
  for(i=0; i<DIM; i++)
  {
    UserWriteF(" x%1d=%11.4E",i, (float)(CVECT(theVertex)[i]) );
  }
  UserWrite("\n");

  if (vopt)       /* verbose: print all information */
  {
    /* print nfather information */
    if (NFATHER(theNode)!=NULL)
    {
      switch (NTYPE(theNode))
      {
      case (CORNER_NODE) :
        UserWriteF(" NFATHER(Node)=" ID_FMTX "\n",
                   ID_PRTX((NODE *)NFATHER(theNode)));
        break;
      case (MID_NODE) :
        UserWriteF(" NFATHER(Edge)=" EDID_FMTX "\n",
                   EDID_PRTX((EDGE *)NFATHER(theNode)));
        break;
      default :
        break;
      }
    }
    /* print nfather information */
    if (SONNODE(theNode)!=NULL)
    {
      UserWriteF(" SONNODE=" ID_FMTX "\n",ID_PRTX(SONNODE(theNode)));
    }

    /* line 3 */	/* print vertex father information */
    if (VFATHER(theVertex)!=NULL)
    {
      UserWriteF("   VERTEXFATHER=" EID_FMTX " ",
                 EID_PRTX(VFATHER(theVertex)));
      for(i=0; i<DIM; i++)
      {
        UserWriteF("XI[%d]=%11.4E ",i, (float)(LCVECT(theVertex)[i]) );
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
      if (BNDP_BndPDesc(V_BNDP(theVertex),&i,&part))
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
      UserWriteF("   EDGE=%x/%08x ",MYEDGE(theLink),
                 DDD_InfoGlobalId(PARHDR(MYEDGE(theLink))));
                        #else
      UserWrite("   ");
                        #endif
      UserWriteF("NB=" ID_FMTX " CTRL=%8lx NO_OF_ELEM=%3d MIDNODE=" ID_FMTX,
                 ID_PRTX(NBNODE(theLink)),(long)CTRL(theLink),
                 NO_OF_ELEM(MYEDGE(theLink)),ID_PRTX(MIDNODE(MYEDGE(theLink))));

      theVertex=MYVERTEX(NBNODE(theLink));
      for(i=0; i<DIM; i++)
      {
        UserWriteF("x%1d=%11.4 ",i, (float)(CVECT(theVertex)[i]) );
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
  ELEMENT *SonList[MAX_SONS];

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
  UserWriteF("ELEMID=" EID_FFMTE " %5s %5s CTRL=%8lx CTRL2=%8lx REFINE=%2d MARK=%2d LEVEL=%2d",
             EID_PRTE(theElement),ekind,etype,
             (long)CTRL(theElement),(long)FLAG(theElement),REFINE(theElement),MARK(theElement),LEVEL(theElement));
  if (COARSEN(theElement)) UserWrite(" COARSEN");
  UserWrite("\n");

  if (vopt)
  {
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      UserWriteF("    N%d=" ID_FMTX,i,ID_PRTX(CORNER(theElement,i)));
    }
    UserWriteF("\n");
    if (EFATHER(theElement))
      UserWriteF("    FA=" EID_FMTX ,EID_PRTX(EFATHER(theElement)));
    else
      UserWriteF("    FA=NULL");

    UserWriteF("  NSONS=%d\n",NSONS(theElement));
    /* TODO: delete this
       #ifdef __TWODIM__
                    for (i=0; i<SONS_OF_ELEM(theElement); i++)
                            if (SON(theElement,i)!=NULL)
                            {
                                    UserWriteF("    S%d=" EID_FMTX ,
                                            i,EID_PRTX(SON(theElement,i)));
                            }
       #endif
       #ifdef __THREEDIM__
     */
    if (GetAllSons(theElement,SonList)!=0) return;
    for (i=0; SonList[i] != NULL; i++)
    {
      UserWriteF("    S%d=" EID_FMTX ,i,EID_PRTX(SonList[i]));
      if ((i+1)%4 == 0) UserWrite("\n");
    }
    /*
       #endif
     */
    UserWrite("\n");
  }
  if (nbopt)
  {
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if (NBELEM(theElement,i)!=NULL)
      {
        UserWriteF("    NB%d=%ld ",i,(long)ID(NBELEM(theElement,i)));
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
                                                #if defined(ModelP) && defined(__THREEDIM__)
          UserWriteF("    NODE[ID=%ld]: ",
                     (long)(ID(CORNER(theElement,
                                      CORNER_OF_SIDE(theElement,i,j)))));
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
  UserWriteF("VTYPE=%d(%c) ",VTYPE(theVector),FMT_T2N(theFormat,VTYPE(theVector)));

  /* print object type of vector */
  if (VOTYPE(theVector)==NODEVEC)
  {
    theNode = (NODE*)VOBJECT(theVector);
    UserWriteF("NODE-V IND=" VINDEX_FFMTE " nodeID=" ID_FMTX
               "                VCLASS=%1d VNCLASS=%1d\n",
               VINDEX_PRTE(theVector),ID_PRTX(theNode),VCLASS(theVector),VNCLASS(theVector));
  }
  if (VOTYPE(theVector)==EDGEVEC)
  {
    theEdge = (EDGE*)VOBJECT(theVector);
    UserWriteF("EDGE-V IND=" VINDEX_FFMTE " fromID=" ID_FFMT
               " to__ID=%7ld VCLASS=%1d VNCLASS=%1d\n",
               VINDEX_PRTE(theVector),ID_PRT(NBNODE(LINK0(theEdge))),
               ID(NBNODE(LINK1(theEdge))),VCLASS(theVector),VNCLASS(theVector));
  }
        #ifdef __THREEDIM__
  if (VOTYPE(theVector)==SIDEVEC)
  {
    theElement = (ELEMENT*)VOBJECT(theVector);
    UserWriteF("SIDE-V IND=" VINDEX_FFMTE " elemID=" EID_FFMT
               "                VCLASS=%1d VNCLASS=%1d\n",
               VINDEX_PRTE(theVector),EID_PRT(theElement),
               VCLASS(theVector),VNCLASS(theVector));
  }
        #endif
  if (VOTYPE(theVector)==ELEMVEC)
  {
    theElement = (ELEMENT*)VOBJECT(theVector);
    UserWriteF("ELEM-V IND=" VINDEX_FFMTE " elemID=" EID_FFMT
               "                VCLASS=%1d VNCLASS=%1d\n",
               VINDEX_PRTE(theVector),EID_PRT(theElement),
               VCLASS(theVector),VNCLASS(theVector));
  }

  /* print vector data if */
  if (dataopt && FMT_PR_VEC(theFormat)!=NULL)
  {
    /* print skip flags */
    INT_2_bitpattern(VECSKIP(theVector),buffer);
    UserWriteF("  skip=%s\n",buffer);

    /* print data */
    Data = (void*)(&VVALUE(theVector,0));
    if ((*(FMT_PR_VEC(theFormat)))(VTYPE(theVector),Data,"   ",buffer))
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
      if (dataopt && theFormat->PrintMatrix!=NULL)
      {
        Data = (void*)(&MVALUE(theMatrix,0));
        if ((*(FMT_PR_MAT(theFormat)))(MROOTTYPE(theMatrix)*MAXVECTORS+MDESTTYPE(theMatrix),Data,"       ",buffer))
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
    UserWriteF("ELEM(ID=%d):\n",ID(theElement));

    if (VEC_DEF_IN_OBJ_OF_MG(theMG,NODEVEC))
    {
      GetVectorsOfNodes(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    }
    if (VEC_DEF_IN_OBJ_OF_MG(theMG,EDGEVEC))
    {
      GetVectorsOfEdges(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    }
                #ifdef __THREEDIM__
    if (VEC_DEF_IN_OBJ_OF_MG(theMG,SIDEVEC))
    {
      GetVectorsOfSides(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    }
                #endif
    if (VEC_DEF_IN_OBJ_OF_MG(theMG,ELEMVEC))
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

/*
   void ListConnections (GRID *theGrid)
   {
        VECTOR *v;
        MATRIX *m;
        INT len;
    buffer[256];

        for (v=PFIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
                len = sprintf(buffer,"%d: prio=%d  %8d ->",me,PRIO(v),GID(v));
                for (m=START(v); m!=NULL; m=MNEXT(m))
                    len += sprintf(buffer+len," %8d",GID(MDEST(m)));
                printf("%s\n",buffer);
        }
   }
 */

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
   This function adds an node to the selection buffer or removes it if it is already
   in the list.

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
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      /* remove g from list (replace it with last object) */
      SELECTIONSIZE(theMG)--;
      SELECTIONOBJECT(theMG,i) = SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG));
      return(GM_OK);
    }

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
   This function adds an element to the selection buffer or removes it if it is already
   in the list.

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
      /* remove g from list (replace it with last object) */
      SELECTIONSIZE(theMG)--;
      SELECTIONOBJECT(theMG,i) = SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG));
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
   This function adds a vector to the selection buffer or removes it if it is already
   in the list.

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
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      /* remove g from list (replace it with last object) */
      SELECTIONSIZE(theMG)--;
      SELECTIONOBJECT(theMG,i) = SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG));
      return(GM_OK);
    }

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

static INT GetAngle(DOUBLE *angle,DOUBLE *n1, DOUBLE *n2)
{
  DOUBLE norm1,norm2,s;

  V_DIM_EUKLIDNORM(n1,norm1);
  V_DIM_EUKLIDNORM(n2,norm2);

  if ((norm1<SMALL_D)||(norm2<SMALL_D))
    return(1);

  V_DIM_SCALE(1.0/norm1,n1);
  V_DIM_SCALE(1.0/norm2,n2);
  V_DIM_SCALAR_PRODUCT(n1,n2,s);

  s=MIN(1,s); s=MAX(-1,s);

  *angle=acos(s);
  return(0);
}

static INT SetNormal(DOUBLE *n, DOUBLE **x, INT nc)
{
  DOUBLE v[DIM];

        #ifdef __TWODIM__
  if (nc!=2) return(1);
  V_DIM_SUBTRACT(x[1],x[0],v);
  n[0] = v[1];
  n[1] = -v[0];
        #endif

        #ifdef __THREEDIM__
  DOUBLE w[DIM];

  /* this is the primitive form, it would be better to take a fit */
  if (nc<3) return(1);
  V_DIM_SUBTRACT(x[1],x[0],v);
  V_DIM_SUBTRACT(x[nc-1],x[0],w);
  V3_VECTOR_PRODUCT(v,w,n);
        #endif

  return(0);
}

INT MinMaxAngle (ELEMENT *theElement, DOUBLE *amin, DOUBLE *amax)
{
  INT error,i,s1,s2,tag;
  DOUBLE angle,*x[MAX_CORNERS_OF_SIDE],n1[DIM],n2[DIM];

  error=GM_OK;
  tag=TAG(theElement);
  for (s1=0; s1<SIDES_OF_TAG(tag); s1++)
  {
    /* get corner coordinates of side and evaluate normal */
    for (i=0; i<CORNERS_OF_SIDE_TAG(tag,s1); i++)
      x[i]=CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE_TAG(tag,s1,i))));

    if (SetNormal(n1,x,CORNERS_OF_SIDE_TAG(tag,s1))!=0) {
      error=GM_ERROR; continue;
    }

    for (s2=s1+1; s2<SIDES_OF_TAG(tag); s2++)
#ifdef __TWODIM__
      if ((CORNER_OF_SIDE_TAG(tag,s1,0)==CORNER_OF_SIDE_TAG(tag,s2,1)) ||
          (CORNER_OF_SIDE_TAG(tag,s1,1)==CORNER_OF_SIDE_TAG(tag,s2,0)))
#endif
#ifdef __THREEDIM__
      if (EDGE_OF_TWO_SIDES_TAG(tag,s1,s2)!=-1)
#endif
    {
      /* get corner coordinates of side and evaluate normal */
      for (i=0; i<CORNERS_OF_SIDE_TAG(tag,s2); i++)
        x[i]=CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE_TAG(tag,s2,i))));

      if (SetNormal(n2,x,CORNERS_OF_SIDE_TAG(tag,s2))!=0) {
        error=GM_ERROR; continue;
      }

      if (GetAngle(&angle,n1,n2)!=0) {
        error=GM_ERROR; continue;
      }

      *amax = MAX(*amax,PI-angle);
      *amin = MIN(*amin,PI-angle);
    }
  }

  return(error);
}


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
   SetEdgeAndNodeSubdomainFromElements - set subdomain id on level 0 edges

   SYNOPSIS:
   INT SetEdgeAndNodeSubdomainFromElements (GRID *theGrid)

   PARAMETERS:
   .  id - the id of the block to be allocated

   DESCRIPTION:
   This function sets the subdomain id taken from the elements for level 0 edges.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if error occured
   D*/
/****************************************************************************/

static INT SetEdgeAndNodeSubdomainFromElements (GRID *theGrid)
{
  ELEMENT *theElement;
  NODE *n0,*n1;
  EDGE *ed;
  INT s_id,s,i,k;

  /* first set subdomain id for all edges */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    /* all edges of the element acquire the subdomain id of the element */
    s_id = SUBDOMAIN(theElement);
    for (k=0; k<EDGES_OF_ELEM(theElement); k++)
    {
      n0 = CORNER(theElement,CORNER_OF_EDGE(theElement,k,0));
      n1 = CORNER(theElement,CORNER_OF_EDGE(theElement,k,1));
      ed = GetEdge(n0,n1);
      ASSERT(ed!=NULL);
      SETEDSUBDOM(ed,s_id);
    }
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      SETNSUBDOM(CORNER(theElement,i),s_id);
  }

  /* now change subdomain id for boundary edges to 0 */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
    if (OBJT(theElement)==BEOBJ)
      for (s=0; s<SIDES_OF_ELEM(theElement); s++)
      {
        if (ELEM_BNDS(theElement,s)==NULL)
          continue;

        for (i=0; i<EDGES_OF_SIDE(theElement,s); i++)
        {
          k  = EDGE_OF_SIDE(theElement,s,i);
          n0 = CORNER(theElement,CORNER_OF_EDGE(theElement,k,0));
          n1 = CORNER(theElement,CORNER_OF_EDGE(theElement,k,1));
          ed = GetEdge(n0,n1);
          ASSERT(ed!=NULL);
          SETEDSUBDOM(ed,0);
        }
      }
  IFDEBUG(gm,1)
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    PRINTDEBUG(gm,1,("el(%d)-sd=%d\n",ID(theElement),
                     SUBDOMAIN(theElement)));
    for (k=0; k<EDGES_OF_ELEM(theElement); k++)
    {
      n0 = CORNER(theElement,CORNER_OF_EDGE(theElement,k,0));
      n1 = CORNER(theElement,CORNER_OF_EDGE(theElement,k,1));
      ed = GetEdge(n0,n1);
      PRINTDEBUG(gm,1,("  ed(%d,%d)-sd=%d\n",ID(n0),ID(n1),
                       EDSUBDOM(ed)));
    }
  }
  ENDDEBUG

  return (GM_OK);
}

/****************************************************************************/
/*D
   SetSubdomainIDfromBndInfo - set subdomain id on level 0 elements and edges

   SYNOPSIS:
   INT SetSubdomainIDfromBndInfo (MULTIGRID *theMG)

   PARAMETERS:
   .  id - the id of the block to be allocated

   DESCRIPTION:
   This function sets the subdomain for level 0 elements and edges.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if error occured
   D*/
/****************************************************************************/

INT SetSubdomainIDfromBndInfo (MULTIGRID *theMG)
{
  HEAP *theHeap;
  GRID *theGrid;
  ELEMENT *theElement, *theNeighbor;
  NODE *theNode;
  void *buffer;
  INT i,n,id,nbid,part,j;
  FIFO myfifo;

  /* prepare */
  if (TOPLEVEL(theMG)<0) REP_ERR_RETURN (GM_ERROR);
  theGrid = GRID_ON_LEVEL(theMG,0);
  n = NT(theGrid);        if (n==0) return(0);

  /* allocate fifo and init */
  theHeap = MYMG(theGrid)->theHeap;
  buffer=(void *)GetTmpMem(theHeap,sizeof(ELEMENT*)*n);
  fifo_init(&myfifo,buffer,sizeof(ELEMENT*)*n);
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
    SETUSED(theElement,0);

  /* insert all boundary elements */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
    if (OBJT(theElement)==BEOBJ && !USED(theElement))
    {
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
        if (ELEM_BNDS(theElement,i)!=NULL)
          break;
      assert(i<SIDES_OF_ELEM(theElement));

      /* set id from BNDS */
      if (BNDS_BndSDesc(ELEM_BNDS(theElement,i),&id,&nbid,&part))
        REP_ERR_RETURN (GM_ERROR);
      assert(id>0);
      SETSUBDOMAIN(theElement,id);
      SETUSED(theElement,1);
      fifo_in(&myfifo,(void *)theElement);
      PRINTDEBUG(gm,1,("elem %3d sid %d\n",ID(theElement),SUBDOMAIN(theElement)));
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        theNode = CORNER(theElement,i);
        if (OBJT(MYVERTEX(theNode))==IVOBJ)
          SETNSUBDOM(theNode,id);
      }
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        if (NBELEM(theElement,i)==NULL || SIDE_ON_BND(theElement,i)) continue;
        theNeighbor = NBELEM(theElement,i);
        if (USED(theNeighbor))
          assert(SUBDOMAIN(theElement)==SUBDOMAIN(theNeighbor));
      }
    }

  /* set subdomain id for all elements */
  while(!fifo_empty(&myfifo))
  {
    theElement = (ELEMENT*)fifo_out(&myfifo);
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      if (NBELEM(theElement,i)==NULL) continue;
      theNeighbor = NBELEM(theElement,i);
      if (USED(theNeighbor))
      {
        if (INNER_SIDE(theElement,i))
          assert(SUBDOMAIN(theElement)==SUBDOMAIN(theNeighbor));
        continue;
      }
      SETSUBDOMAIN(theNeighbor,SUBDOMAIN(theElement));
      SETUSED(theNeighbor,1);
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      {
        theNode = CORNER(theElement,j);
        if (OBJT(MYVERTEX(theNode))==IVOBJ)
          SETNSUBDOM(theNode,SUBDOMAIN(theElement));
      }
      fifo_in(&myfifo,(void *)theNeighbor);
    }
  }
  if (SetEdgeAndNodeSubdomainFromElements(theGrid))
    REP_ERR_RETURN (GM_ERROR);

  return (GM_OK);
}

/****************************************************************************/
/*D
   FixCoarseGrid - do all that is necessary to complete the coarse grid

   SYNOPSIS:
   INT FixCoarseGrid (MULTIGRID *theMG)

   PARAMETERS:
   .  id - the id of the block to be allocated

   DESCRIPTION:
   This function does all that is necessary to complete the coarse grid.
   Finally the MG_COARSE_FIXED flag is set.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if error occured
   D*/
/****************************************************************************/

INT FixCoarseGrid (MULTIGRID *theMG)
{
  if (MG_COARSE_FIXED(theMG))
    return (GM_OK);

  if (SetSubdomainIDfromBndInfo(theMG))
    REP_ERR_RETURN (GM_ERROR);

  /* set this flag here because it is checked by CreateAlgebra */
  if (CreateAlgebra(theMG) != GM_OK)
    REP_ERR_RETURN (GM_ERROR);

  ReleaseTmpMem(MGHEAP(theMG));

  return (GM_OK);
}

/****************************************************************************/
/*D
   InitUGManager - Init what is necessary

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
