// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugm.c                                                         */
/*                                                                          */
/* Purpose:   unstructured grid manager                                     */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   09.03.92 begin, ug version 2.0                                */
/*            Aug 28 1996, ug version 3.4                                   */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef __MPW32__
#pragma segment ugm
#endif

/****************************************************************************/
/*                                                                          */
/*        defines to exclude functions                                      */
/*                                                                          */
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

#include "ugdevices.h"

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
#include "ugstruct.h"
#ifdef DYNAMIC_MEMORY_ALLOCMODEL
#include "mgheapmgr.h"
#endif

#ifdef ModelP
#include "identify.h"
#endif

#include "cw.h"

USING_UG_NAMESPACE
  USING_UGDIM_NAMESPACE

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

#ifdef __MWCW__
#define printf                                          PrintDebug
#endif

/* local refinement hack */
#undef _SCHALE_X_

/* macro for controlling debugging output by conditions on objects */
#define UGM_CDBG(x,y)

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*		  in the corresponding include file!)                         */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

#if defined ModelP && defined __OVERLAP2__
INT ce_NO_DELETE_OVERLAP2 = -1;
#endif


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static char buffer[4*256];                      /* general purpose text buffer			*/

static VIRT_HEAP_MGMT *theGenMGUDM; /* general user data space management	*/

static INT theMGDirID;                          /* env var ID for the multigrids		*/
static INT theMGRootDirID;                      /* env dir ID for the multigrids		*/

static UINT UsedOBJT;           /* for the dynamic OBJECT management	*/

/* used by OrderNodesInGrid */
static const INT *Order,*Sign;
static DOUBLE InvMeshSize;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

REP_ERR_FILE;

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

static NODE *CreateNode (GRID *theGrid, VERTEX *vertex, GEOM_OBJECT *Father, INT NodeType, INT with_vector);
static VERTEX *CreateBoundaryVertex     (GRID *theGrid);
static VERTEX *CreateInnerVertex (GRID *theGrid);

static INT DisposeVertex (GRID *theGrid, VERTEX *theVertex);
static INT DisposeEdge (GRID *theGrid, EDGE *theEdge);


#ifdef __PERIODIC_BOUNDARY__

#define SMALL_DOUBLE 1e-4
#define MAX_PERIODIC_PROCS      128

/* vid is unique id of vector */
#ifndef ModelP
#define PVID VINDEX
#else
#define PVID GID
#endif

#ifdef __TWODIM__
#define MIN_COORD(a,b) (((*a)[0]< (*b)[0]) ? (a) : (((*a)[0] > (*b)[0]) ? (b) : (((*a)[1] < (*b)[1]) ? (a) : (b))))
#else
#define MIN_COORD(a,b) (((*a)[0]< (*b)[0]) ? (a) : (((*a)[0] > (*b)[0]) ? (b) : (((*a)[1] < (*b)[1]) ? (a) : (((*a)[1] > (*b)[1]) ? (b) : (((*a)[2] < (*b)[2]) ? (a) : (b))))))
#endif

#ifdef __TWODIM__
#define NPERVEC 2
#else
#define NPERVEC 4
#endif

typedef struct id_tupel
{
  INT tpl[NPERVEC+1];
}
IDTPL;

typedef struct periodic_entries
{
  NODE *node;
  DOUBLE_VECTOR coord;
  INT periodic_id;
  INT n;
  VECTOR *vp[NPERVEC];
}
PERIODIC_ENTRIES;

static PeriodicBoundaryInfoProcPtr PeriodicBoundaryInfo = NULL;
#endif


/****************************************************************************/
/** \brief Get an object type id not occupied in theMG
 *
 * This function gets an object type id not occupied in theMG.
 *
 * @return <ul>
 *   <li> id of object type if ok </li>
 *   <li> -1 when error occurred </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetFreeOBJT ()
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
/** \brief Release an object type id not needed anymore
 *
 * @param  type - object type
 *
 * This function releases an object type id not needed anymore.
 *
 * @return <ul>
 * <li>   GM_OK if ok </li>
 * <li>   GM_ERROR when error occured. </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX ReleaseOBJT (INT type)
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
/** \brief Get an object from free list if possible
   \fn GetMemoryForObject

 * @param  theMG - pointer to multigrid
 * @param  size - size of the object
 * @param  type - type of the requested object

   This function gets an object of type `type` from free list if possible,
   otherwise it allocates memory from the multigrid heap using 'GetMem'.

   @return <ul>
   <li>   pointer to an object of the requested type </li>
   <li>   NULL if object of requested type is not available </li>
 * </ul>
 */
/****************************************************************************/

#ifdef ModelP
static void ConstructDDDObject (void *obj, INT size, INT type)
{
  if (obj!=NULL && type!=NOOBJ)
  {
    memset(obj,0,size);
    /* link this object to DDD management */
    if (HAS_DDDHDR(type))
    {
      DDD_TYPE dddtype = DDDTYPE(type);
      DDD_HDR dddhdr = (DDD_HDR)(((char *)obj) + DDD_InfoHdrOffset(dddtype));
      DDD_HdrConstructor(dddhdr, dddtype, PrioMaster, 0);
    }
  }
  return;
}
#endif

#ifndef DYNAMIC_MEMORY_ALLOCMODEL
void * NS_DIM_PREFIX GetMemoryForObject_par (HEAP *theHeap, INT size, INT type)
{
  void *obj = GetFreelistMemory(theHeap, size);

        #ifdef ModelP
  if (type!=MAOBJ && type!=COOBJ)
    ConstructDDDObject(obj,size,type);
        #endif

  return obj;
}
#else
void * NS_DIM_PREFIX GetMemoryForObjectNew (HEAP *theHeap, INT size, INT type)
{
  void                    *obj;

        #ifdef Debug
  check_of_getcallstack = 1;
        #endif

  if (usefreelistmemory == 1)
    obj = GetFreelistMemory(theHeap, size);
  else
  {
                #ifdef Debug
    switch (type)
    {
    case MAOBJ :
    case VEOBJ :
    case GROBJ :
    case BLOCKVOBJ :
      break;
    default : assert(0);
    }
                #endif
    obj = GetMem(theHeap,size,FROM_BOTTOM);
    if (obj != NULL)
      memset(obj,0,size);
  }

        #ifdef Debug
  check_of_getcallstack = 0;
        #endif

        #ifdef ModelP
  if (type!=MAOBJ && type!=COOBJ)
    ConstructDDDObject(obj,size,type);
        #endif

  return obj;
}
#endif

/****************************************************************************/
/** \brief  Put an object in the free list
   \fn PutFreeObject

 * @param  theMG - pointer to multigrid
 * @param  object - object to insert in free list
 * @param  size - size of the object
 * @param  type - type of the requested object

   This function puts an object in the free list.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 when error occured. </li>
 * </ul> */
/****************************************************************************/

#ifdef ModelP
static void DestructDDDObject(void *object, INT type)
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
  return;
}
#endif

#ifndef DYNAMIC_MEMORY_ALLOCMODEL
INT NS_DIM_PREFIX PutFreeObject_par (HEAP *theHeap, void *object, INT size, INT type)
{
        #ifdef ModelP
  if (type!=MAOBJ && type!=COOBJ)
    DestructDDDObject(object,type);
        #endif

  return (PutFreelistMemory(theHeap, object, size));
}
#else
INT NS_DIM_PREFIX PutFreeObjectNew (HEAP *theHeap, void *object, INT size, INT type)
{
  INT err;

        #ifdef ModelP
  if (type!=MAOBJ && type!=COOBJ)
    DestructDDDObject(object,type);
        #endif

        #ifdef Debug
  check_of_putcallstack = 1;
        #endif

  if (usefreelistmemory == 1)
  {
    err = PutFreelistMemory(theHeap, object, size);
                #ifdef Debug
    check_of_putcallstack = 0;
                #endif
    return (err);
  }

        #ifdef Debug
  switch (type)
  {
  case MAOBJ :
  case VEOBJ :
  case GROBJ :
  case BLOCKVOBJ :
    break;
  default : assert(0);
  }
        #endif

        #ifdef Debug
  check_of_putcallstack = 0;
        #endif

  /* memory is freed by release */
  return(0);
}
#endif

/****************************************************************************/
/** \brief Return pointer to a new boundary vertex structure
 *
 * @param theGrid grid where vertex should be inserted
 *
 * This function creates and initializes a new boundary vertex structure
 * and returns a pointer to it.
 *
 * @return <ul>
 *    <li> pointer to requested object </li>
 *    <li> NULL if out of memory </li>
 * </ul>
 */
/****************************************************************************/

static VERTEX *CreateBoundaryVertex (GRID *theGrid)
{
  VERTEX *pv;
  INT ds;
  INT i;

  pv = (VERTEX*)GetMemoryForObject(MYMG(theGrid),sizeof(struct bvertex),BVOBJ);
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
  SETNOOFNODE(pv,0);
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
  /* SETVXPRIO(pv,PrioMaster); */
        #endif

  /* insert in vertex list */
  GRID_LINK_VERTEX(theGrid,pv,PrioMaster);

  return(pv);
}

/****************************************************************************/
/** \brief Return pointer to a new inner vertex structure
 *
 * @param theGrid grid where vertex should be inserted
 *
 * This function creates and initializes a new inner vertex structure
 * and returns a pointer to it.
 *
 * @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

static VERTEX *CreateInnerVertex (GRID *theGrid)
{
  VERTEX *pv;
  INT ds;
  INT i;

  pv = (VERTEX*)GetMemoryForObject(MYMG(theGrid),sizeof(struct ivertex),IVOBJ);
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
  SETNOOFNODE(pv,0);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
        #ifdef TOPNODE
  TOPNODE(pv) = NULL;
        #endif
  SETMOVE(pv,DIM);
        #ifdef ModelP
  DDD_AttrSet(PARHDRV(pv),GRID_ATTR(theGrid));
  /* SETVXPRIO(pv,PrioMaster); */
        #endif
  for (i=0; i<DIM; i++) LCVECT(pv)[i] = 0.0;

  /* insert in vertex list */
  GRID_LINK_VERTEX(theGrid,pv,PrioMaster);

  return(pv);
}

/****************************************************************************/
/** \brief Return pointer to a new node structure

 * @param  theGrid - grid where vertex should be inserted
 * @param  vertex  - vertex of the node
 * @param  FatherNode - father node (may be NULL)
 * @param  NodeType - node type (CORNER_NODE..., cf. gm.h)
 * @param   with_vector

   This function creates and initializes a new node structure
   and returns a pointer to it.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
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
  /* SETPRIO(pn,PrioMaster); */
        #endif
  ID(pn) = (theGrid->mg->nodeIdCounter)++;
  START(pn) = NULL;
  SONNODE(pn) = NULL;
  if (NELIST_DEF_IN_GRID(theGrid)) NDATA(pn) = NULL;
  MYVERTEX(pn) = vertex;
  if (NOOFNODE(vertex)<NOOFNODEMAX)
    INCNOOFNODE(vertex);
  else
    ASSERT(0);
  /* priliminary */
  if (Father != NULL)
    if ((OBJT(Father) == IEOBJ) || (OBJT(Father) == BEOBJ))
      Father = NULL;
  SETNFATHER(pn,Father);
  SETNTYPE(pn,NodeType);
  SETNCLASS(pn,3),
  SETNNCLASS(pn,0);
  if (OBJT(vertex) == BVOBJ)
    SETNSUBDOM(pn,0);
  else if (VFATHER(vertex) != NULL)
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
/** \brief Return pointer to a new node structure on an edge

 * @param   theGrid - grid where vertex should be inserted
 * @param   FatherNode - node father

   This function creates and initializes a new node structure
   at the midpoint of an element edge and returns a pointer to it.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

NODE *NS_DIM_PREFIX CreateSonNode (GRID *theGrid, NODE *FatherNode)
{
  NODE *pn;
  VERTEX *theVertex;

  theVertex = MYVERTEX(FatherNode);

  pn = CreateNode(theGrid,theVertex,(GEOM_OBJECT *)FatherNode,CORNER_NODE,1);
  if (pn == NULL)
    return(NULL);
  SONNODE(FatherNode) = pn;

        #ifdef TOPNODE
  TOPNODE(theVertex) = pn;
        #endif

  return(pn);
}

/****************************************************************************/
/** \brief Return pointer to a new node structure on an edge

 * @param   theGrid - grid where node should be inserted
 * @param   theElement - pointer to an element
 * @param   theVertex - pointer to vertex if already existing
 * @param   edge - id of an element edge

   This function creates and initializes a new node structure
   at the midpoint of an element edge and returns a pointer to it.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

NODE *NS_DIM_PREFIX CreateMidNode (GRID *theGrid, ELEMENT *theElement, VERTEX *theVertex, INT edge)
{
  NODE *theNode;
  EDGE *theEdge;
  VERTEX *v0,*v1;
  BNDP *bndp;
  DOUBLE *local,*x[MAX_CORNERS_OF_ELEM];
  DOUBLE_VECTOR bnd_global,global;
  DOUBLE diff;
  INT n,co0,co1,move,part,vertex_null;

  co0 = CORNER_OF_EDGE(theElement,edge,0);
  co1 = CORNER_OF_EDGE(theElement,edge,1);
  v0 = MYVERTEX(CORNER(theElement,co0));
  v1 = MYVERTEX(CORNER(theElement,co1));
  V_DIM_LINCOMB(0.5, CVECT(v0), 0.5, CVECT(v1), global);

  /* set MIDNODE pointer */
  theEdge = GetEdge(CORNER(theElement,co0),CORNER(theElement,co1));
  ASSERT(theEdge!=NULL);

  /* allocate vertex */
  vertex_null = (theVertex==NULL);
  if (theVertex==NULL)
  {
    if ((OBJT(v0) == BVOBJ) && (OBJT(v1) == BVOBJ))
#ifdef __TWODIM__
      if (OBJT(theElement) == BEOBJ)
        if (SIDE_ON_BND(theElement,edge))
#endif
#ifdef __THREEDIM__
      if (EDSUBDOM(theEdge) == 0)
#endif
    {
      bndp = BNDP_CreateBndP(MGHEAP(MYMG(theGrid)),V_BNDP(v0),V_BNDP(v1),0.5);
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
        PRINTDEBUG(gm,1,("local = %f %f %f\n",local[0],local[1],local[2]));
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
  }
  /* allocate node */
  theNode = CreateNode(theGrid,theVertex,(GEOM_OBJECT *)theEdge,MID_NODE,1);
  if (theNode==NULL && vertex_null)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }

  MIDNODE(theEdge) = theNode;
        #ifdef TOPNODE
  if (TOPNODE(theVertex)==NULL || LEVEL(TOPNODE(theVertex))<LEVEL(theNode))
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


NODE * NS_DIM_PREFIX GetMidNode (ELEMENT *theElement, INT edge)
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

  /* this is a bad place for the following code (s.l. 981015) */
  theVertex = MYVERTEX(theNode);
  if (theVertex!=NULL && VFATHER(theVertex) == NULL) {
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


/****************************************************************************/
/** \brief ???
 */
/****************************************************************************/

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


/****************************************************************************/
/** \brief Return pointer to a new node structure on a side (3d)

 * @param   theGrid - grid where vertex should be inserted
 * @param   theElement - pointer to an element
 * @param   theVertex - pointer vertex
 * @param   side - id of an element side

   This function creates and initializes a new node structure
   at the midpoint of an element side and returns a pointer to it.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

NODE *NS_DIM_PREFIX CreateSideNode (GRID *theGrid, ELEMENT *theElement, VERTEX *theVertex, INT side)
{
  DOUBLE_VECTOR bnd_global,global,local,bnd_local;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  NODE *theNode;
  BNDP *bndp;
  BNDS *bnds;
  DOUBLE fac, diff;
  INT n,j,k,move,part,vertex_null;

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

  /* check if boundary vertex */
  vertex_null = (theVertex==NULL);
  if (theVertex==NULL)
  {
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
            PRINTDEBUG(gm,1,("local = %f %f %f\n",local[0],local[1],local[2]));
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
  }
  /* create node */
  theNode = CreateNode(theGrid,theVertex,
                       (GEOM_OBJECT *)theElement,SIDE_NODE,1);
  if (theNode==NULL && vertex_null)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
        #ifdef TOPNODE
  if (TOPNODE(theVertex) == NULL || LEVEL(TOPNODE(theVertex))<LEVEL(theNode))
    TOPNODE(theVertex) = theNode;
        #endif
  theGrid->status |= 1;

  return(theNode);
}


/****************************************************************************/
/*
   GetSideNode -

   SYNOPSIS:
   NODE *GetSideNode (ELEMENT *theElement, INT side);

   PARAMETERS:
 * @param   theElement
 * @param   side

   DESCRIPTION:

   @return
 */
/****************************************************************************/

static NODE *GetSideNodeX (ELEMENT *theElement, INT side, INT n,
                           NODE **MidNodes)
{
  ELEMENT *theFather;
  NODE *theNode;
  VERTEX *theVertex;
  LINK *theLink0,*theLink1,*theLink2,*theLink3;
  DOUBLE fac,*local;
  INT i,m;

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
            if (theFather == theElement) {
                        #ifndef ModelP
              /* HEAPFAULT in theFather possible,
                 if in a previous call of DisposeElement
                 some son is not reached by GetAllSons */
              assert(ONSIDE(theVertex) == side);
                        #endif
              SETONSIDE(theVertex,side);
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
                        #ifndef ModelP
            /* HEAPFAULT in theFather possible,
               if in a previous call of DisposeElement
               some son is not reached by GetAllSons */
            else
              assert(0);
                        #endif
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
                        #ifdef ModelP
            SETONSIDE(theVertex,side);
            return(theNode);
                        #endif
          }
          else if (theFather == NBELEM(theElement,side))
          {
            INT nbside = SideOfNbElement(theElement,side);
            if (nbside==ONSIDE(theVertex))
            {
              SETONNBSIDE(theVertex,side);
              return(theNode);
            }
                        #ifdef ModelP
            VFATHER(theVertex) = theElement;
            SETONSIDE(theVertex,side);
            SETONNBSIDE(theVertex,nbside);
            return(theNode);
                        #endif
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
                    #ifdef ModelP
          else {
            return(theNode);
          }
                    #endif
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
          SETONSIDE(theVertex,side);
          return(theNode);
        }
        else if (theFather == NBELEM(theElement,side)) {
          SETONNBSIDE(theVertex,side);
          return(theNode);
        }
        return(theNode);
      }
    }
  }
    #endif

  return(NULL);
}

NODE * NS_DIM_PREFIX GetSideNode (ELEMENT *theElement, INT side)
{
  ELEMENT *theFather;
  NODE *theNode;
  NODE *MidNodes[MAX_EDGES_OF_SIDE];
  VERTEX *theVertex;
  LINK *theLink0,*theLink1,*theLink2,*theLink3;
  DOUBLE fac,*local;
  INT i,k,n;

  n = 0;
  for (i=0; i<EDGES_OF_SIDE(theElement,side); i++) {
    theNode = GetMidNode(theElement,EDGE_OF_SIDE(theElement,side,i));
    if (theNode != NULL)
      MidNodes[n++] = theNode;
                #ifndef ModelP
    else return(NULL);
                #endif
  }
  PRINTDEBUG(gm,2,(PFMT " GetSideNode(): elem=" EID_FMTX
                   " side=%d nb. of midnodes=%d\n",
                   me,EID_PRTX(theElement),side,n));
#ifdef ModelP
  if (TAG(theElement)==PYRAMID && side!=0) return(NULL);
#endif
  theNode = GetSideNodeX(theElement,side,n,MidNodes);
    #ifdef ModelP
  if (theNode != NULL)
    return(theNode);
  if (n < 3)
    return(NULL);
  for (i=0; i<n; i++)
  {
    NODE *MidNodes1[MAX_EDGES_OF_SIDE-1];
    INT j,m;

    m = 0;
    for (j=0; j<n; j++) {
      if (i == j) continue;
      MidNodes1[m++] = MidNodes[j];
    }
    theNode = GetSideNodeX(theElement,side,n-1,MidNodes1);
    if (theNode != NULL)
      return(theNode);
  }
  if (n < 4)
    return(NULL);
  for (i=1; i<n; i++)
    for (k=0; k<i; k++)
    {
      NODE *MidNodes1[MAX_EDGES_OF_SIDE-2];
      INT j,m;

      m = 0;
      for (j=0; j<n; j++) {
        if (i == j) continue;
        if (k == j) continue;
        MidNodes1[m++] = MidNodes[j];
      }
      theNode = GetSideNodeX(theElement,side,n-2,MidNodes1);
      if (theNode != NULL)
        return(theNode);
    }
    #endif

  return(theNode);
}

/****************************************************************************/
/** \brief ???
   PARAMETERS:
 * @param   theElement
 * @param   theNode

   DESCRIPTION:

   @return ???
 */
/****************************************************************************/

static int CountSideNodes (ELEMENT *e)
{
  int i,side;
  NODE *n;

  side = 0;
  for (i=0; i<CORNERS_OF_ELEM(e); i++)
  {
    n = CORNER(e,i);
    if (SIDETYPE(n)) side++;
  }
  return(side);
}

int GetSideIDFromScratchSpecialRule17Pyr (ELEMENT *theElement, NODE *theNode)
{
  int i,k,l,nodes,cnodes,snodes;
  ELEMENT *f = EFATHER(theElement);
  NODE *fnode,*enode;
  int side = SIDES_OF_ELEM(f);

        #ifdef Debug
  assert(TAG(theElement)==PYRAMID);
  snodes = cnodes = 0;
  for (l=0; l<CORNERS_OF_ELEM(theElement); l++)
  {
    enode = CORNER(theElement,l);
    if (CORNERTYPE(enode)) cnodes++;
    if (SIDETYPE(enode)) snodes++;
  }
  assert(snodes == 1);
  assert(cnodes == 4);
        #endif

  for (i=0; i<SIDES_OF_ELEM(f); i++)
  {
    nodes = 0;
    for (k=0; k<CORNERS_OF_SIDE(f,i); k++)
    {
      fnode = CORNER(f,CORNER_OF_SIDE(f,i,k));
      for (l=0; l<CORNERS_OF_ELEM(theElement); l++)
      {
        enode = CORNER(theElement,l);
        if (enode == SONNODE(fnode)) nodes++;
      }
    }
    assert(nodes==0 || nodes==2 || nodes==4);
                #ifdef Debug
    if (nodes == 0) side = i;
                #else
    if (nodes == 0) return(i);
                #endif
  }

  assert(side<SIDES_OF_ELEM(f));
  return(side);
}


int GetSideIDFromScratchSpecialRule22Tet (ELEMENT *theElement, NODE *theNode)
{
  int i,k,l,nodes,cnodes,mnodes,snodes,midnodes;
  ELEMENT *f = EFATHER(theElement);
  NODE *fnode,*enode;
  EDGE *edge;
  int side = SIDES_OF_ELEM(f);

        #ifdef Debug
  assert(TAG(theElement)==TETRAHEDRON);
  snodes = cnodes = mnodes = 0;
  for (l=0; l<CORNERS_OF_ELEM(theElement); l++)
  {
    enode = CORNER(theElement,l);
    if (CORNERTYPE(enode)) cnodes++;
    if (MIDTYPE(enode)) mnodes++;
    if (SIDETYPE(enode)) snodes++;
  }
  assert(cnodes == 2);
  assert(mnodes == 1);
  assert(snodes == 1);
        #endif

  for (i=0; i<SIDES_OF_ELEM(f); i++)
  {
    nodes = 0;
    midnodes = 0;
    for (k=0; k<CORNERS_OF_SIDE(f,i); k++)
    {
      fnode = CORNER(f,CORNER_OF_SIDE(f,i,k));

      edge = GetEdge(CORNER_OF_SIDE_PTR(f,i,k),
                     CORNER_OF_SIDE_PTR(f,i,(k+1)%CORNERS_OF_SIDE(f,i)));
      assert(edge != NULL);

      for (l=0; l<CORNERS_OF_ELEM(theElement); l++)
      {
        enode = CORNER(theElement,l);
        if (enode == SONNODE(fnode)) nodes++;
        if (enode == MIDNODE(edge)) midnodes++;
      }
    }
    assert(nodes==0 || nodes==1 || nodes==2 || nodes==4);
                #ifdef Debug
    if (nodes==0 && midnodes==1) side = i;
                #else
    if (nodes==0 && midnodes==1) return(i);
                #endif
  }

  assert(side<SIDES_OF_ELEM(f));
  return(side);
}


INT GetSideIDFromScratchSpecialRule (ELEMENT *theElement, NODE *theNode)
{
  int j,l,side;
  ELEMENT *f = EFATHER(theElement);

  assert(TAG(f)==HEXAHEDRON);
  assert(ECLASS(theElement)==GREEN_CLASS);
  assert(NSONS(f)==9 || NSONS(f)==11 || EHGHOST(theElement));

  if (TAG(theElement)==PYRAMID)
  {
    return(GetSideIDFromScratchSpecialRule17Pyr(theElement,theNode));
  }

  assert(TAG(theElement)==TETRAHEDRON);
  /* centroid tetrahedron of special rule 22 */
  if (CountSideNodes(theElement) == 2)
  {
    /* if side not found search over neighbor */
    for (j=0; j<SIDES_OF_ELEM(theElement); j++)
    {
      ELEMENT *nb = NBELEM(theElement,j);

      if (nb == NULL) continue;

      for (l=0; l<CORNERS_OF_ELEM(nb); l++)
        if (theNode == CORNER(nb,l))
          return(GetSideIDFromScratch(nb,theNode));
    }
  }

  assert(CountSideNodes(theElement)==1);

  return(GetSideIDFromScratchSpecialRule22Tet(theElement,theNode));
}

INT NS_DIM_PREFIX GetSideIDFromScratch (ELEMENT *theElement, NODE *theNode)
{
  ELEMENT *theFather;
  NODE *nd[MAX_EDGES_OF_ELEM];
  EDGE *edge;
  INT i,j,k,l,cnt;

  ASSERT(NTYPE(theNode) == SIDE_NODE);

  theFather = EFATHER(theElement);

  /* determine midnodes of father */
  for (i=0; i<EDGES_OF_ELEM(theFather); i++)
  {
    edge = GetEdge(CORNER_OF_EDGE_PTR(theFather,i,0),
                   CORNER_OF_EDGE_PTR(theFather,i,1));
    nd[i] = MIDNODE(edge);
  }

  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
  {
    if (3 == CORNERS_OF_SIDE(theElement,j)) continue;

    for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++)
      if (theNode == CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
        break;
    if (l == CORNERS_OF_SIDE(theElement,j)) continue;

    for (i=0; i<SIDES_OF_ELEM(theFather); i++)
    {
                        #ifdef TET_RULESET
      if (3 == CORNERS_OF_SIDE(theFather,i)) continue;
                        #endif

      cnt = 0;
      for (k=0; k<EDGES_OF_SIDE(theFather,i); k++)
        for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++)
        {
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
            cnt++;
          if (cnt == 2)
            return(i);
        }
    }
  }


  /* if side not found search over neighbor */
  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
  {
    ELEMENT *nb = NBELEM(theElement,j);

    if (3 == CORNERS_OF_SIDE(theElement,j))
      continue;

    if (nb == NULL) continue;

    for (l=0; l<CORNERS_OF_ELEM(nb); l++)
      if (theNode == CORNER(nb,l))
        return(GetSideIDFromScratch(nb,theNode));
  }


  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
  {
    if (4 != CORNERS_OF_SIDE(theElement,j)) continue;
    for (l=0; l<4; l++)
      if (theNode == CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
        break;
    if (l < 4)
    {
      INT l1 = (l+1) % 4;
      INT l2 = (l+3) % 4;

      for (i=0; i<SIDES_OF_ELEM(theFather); i++) {
        if (3 == CORNERS_OF_SIDE(theFather,i)) continue;
        for (k=0; k<EDGES_OF_SIDE(theFather,i); k++) {
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l1)))
            return(i);
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l1)))
            return(i);
        }
      }
    }
  }

  return(GetSideIDFromScratchSpecialRule(theElement,theNode));

  return(SIDES_OF_ELEM(theFather));
}

INT GetSideIDFromScratchOld (ELEMENT *theElement, NODE *theNode)
{
  ELEMENT *theFather;
  NODE *nd[MAX_EDGES_OF_ELEM];
  EDGE *edge;
  INT i,j,k,l,cnt;

  ASSERT(NTYPE(theNode) == SIDE_NODE);

  theFather = EFATHER(theElement);

  /* determine midnodes of father */
  for (i=0; i<EDGES_OF_ELEM(theFather); i++)
  {
    edge = GetEdge(CORNER_OF_EDGE_PTR(theFather,i,0),
                   CORNER_OF_EDGE_PTR(theFather,i,1));
    nd[i] = MIDNODE(edge);
  }

  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
  {
    if (3 == CORNERS_OF_SIDE(theElement,j)) continue;

    for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++)
      if (theNode == CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
        break;
    if (l == CORNERS_OF_SIDE(theElement,j)) continue;

    for (i=0; i<SIDES_OF_ELEM(theFather); i++)
    {
      if (3 == CORNERS_OF_SIDE(theFather,i)) continue;

      cnt = 0;
      for (k=0; k<EDGES_OF_SIDE(theFather,i); k++)
        for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++)
        {
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
            cnt++;
          if (cnt == 2)
            return(i);
        }
    }
  }


  /* if side not found search over neighbor */
  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
  {
    ELEMENT *nb = NBELEM(theElement,j);

    if (3 == CORNERS_OF_SIDE(theElement,j))

      /* treatment of special green rule 17 and 22 */
      if ((((TAG(theElement)==PYRAMID && NSONS(theFather)==9) ||
            (TAG(theElement)==TETRAHEDRON && NSONS(theFather)==11)
            && 2==CountSideNodes(theElement))) &&
          TAG(theFather)==HEXAHEDRON &&
          ECLASS(theElement)==GREEN_CLASS)
        /* not continue */;
      else
        continue;

    if (nb == NULL) continue;

    for (l=0; l<CORNERS_OF_ELEM(nb); l++)
      if (theNode == CORNER(nb,l))
        return(GetSideIDFromScratch(nb,theNode));
  }


  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
  {
    if (4 != CORNERS_OF_SIDE(theElement,j)) continue;
    for (l=0; l<4; l++)
      if (theNode == CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
        break;
    if (l < 4)
    {
      INT l1 = (l+1) % 4;
      INT l2 = (l+3) % 4;

      for (i=0; i<SIDES_OF_ELEM(theFather); i++) {
        if (3 == CORNERS_OF_SIDE(theFather,i)) continue;
        for (k=0; k<EDGES_OF_SIDE(theFather,i); k++) {
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l1)))
            return(i);
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l1)))
            return(i);
        }
      }
    }
  }

  /* treatment of special green rule 17 and 22 */
  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
  {
    for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++)
      if (theNode == CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
        break;
    if (l == CORNERS_OF_SIDE(theElement,j)) continue;

    for (i=0; i<SIDES_OF_ELEM(theFather); i++)
    {
      if (3 == CORNERS_OF_SIDE(theFather,i)) continue;

      cnt = 0;
      for (k=0; k<EDGES_OF_SIDE(theFather,i); k++)
        for (l=0; l<CORNERS_OF_SIDE(theElement,j); l++)
        {
          if (nd[EDGE_OF_SIDE(theFather,i,k)] ==
              CORNER(theElement,CORNER_OF_SIDE(theElement,j,l)))
            cnt++;
          if (cnt==1 && ECLASS(theElement)==GREEN_CLASS &&
              TAG(theElement)==TETRAHEDRON &&
              TAG(theFather)==HEXAHEDRON &&
              (NSONS(theFather)==9 || NSONS(theFather)==11))
          {
            return(i);
          }

        }
    }
  }

  UserWriteF("GetSideIDFromScratch(): e=" EID_FMTX " f=" EID_FMTX "\n",
             EID_PRTX(theElement),EID_PRTX(theFather));
  return(0);
  return(SIDES_OF_ELEM(theFather));
}

#endif /* __THREEDIM__ */


/****************************************************************************/
/** \brief Get the center node of an element of next finer level
 *
 * This function gets the center node of an element of next finer level

   @return <ul>
   <li>   pointer to center node </li>
   <li>   NULL  no node found </li>
   </ul> */
/****************************************************************************/

NODE * NS_DIM_PREFIX GetCenterNode (ELEMENT *theElement)
{
  INT i,j;
  NODE    *theNode;
  ELEMENT *SonList[MAX_SONS],*theSon;

        #ifdef __CENTERNODE__
  return(CENTERNODE(theElement));
        #endif

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
        if (EMASTER(theElement))
          assert(VFATHER(MYVERTEX(theNode)) == theElement);
        return (theNode);
      }
    }
  }
  return (NULL);
}

/****************************************************************************/
/** \brief Allocate a new node on a side of an element
 *
 * Includes vertex
 * best fit boundary coordinates and local coordinates.
 *
 * @return <ul>
 *    <li> pointer to new node </li>
 *    <li> NULL: could not allocate </li>
   </ul> */
/****************************************************************************/
/* #define MOVE_MIDNODE */
NODE * NS_DIM_PREFIX CreateCenterNode (GRID *theGrid, ELEMENT *theElement, VERTEX *theVertex)
{
  DOUBLE *global,*local;
  DOUBLE_VECTOR diff;
  INT n,j,moved,vertex_null;
  VERTEX *VertexOnEdge[MAX_EDGES_OF_ELEM];
  NODE *theNode;
  EDGE *theEdge;
  DOUBLE fac;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
        #ifdef MOVE_MIDNODE
        #ifndef ModelP
  DOUBLE len_opp,len_bnd;
    #endif
    #endif

  /* check if moved side nodes exist */
  CORNER_COORDINATES(theElement,n,x);
  moved = 0;
  vertex_null = (theVertex==NULL);
  if (theVertex==NULL && OBJT(theElement) == BEOBJ) {
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
                #ifdef MOVE_MIDNODE
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
        PRINTDEBUG(gm,1,("CreateCenterNode: global_orig = %f %f %f\n",global[0],global[1],global[2]));
        PRINTDEBUG(gm,1,("CreateCenterNode: diff = %f %f %f\n",diff[0],diff[1],diff[2]));
        V_DIM_LINCOMB(1.0,global,0.5,diff,global);
        PRINTDEBUG(gm,1,("CreateCenterNode: global_mod = %f %f %f\n",global[0],global[1],global[2]));
        SETMOVED(VertexOnEdge[OPPOSITE_EDGE(theElement,j)],1);
        UG_GlobalToLocal(n,(const DOUBLE **)x,global,LCVECT(theVertex));
        PRINTDEBUG(gm,1,("CreateCenterNode: local = %f %f %f\n",LCVECT(theVertex)[0],LCVECT(theVertex)[1],LCVECT(theVertex)[2]));
        SETONEDGE(theVertex,OPPOSITE_EDGE(theElement,j));
        VFATHER(theVertex) = theElement;
      }
    }
            #endif
            #endif
  }

  if (vertex_null)
  {
    theVertex = CreateInnerVertex(theGrid);
    if (theVertex==NULL)
      return(NULL);
    VFATHER(theVertex) = theElement;
  }

  theNode = CreateNode(theGrid,theVertex,(GEOM_OBJECT *)theElement,CENTER_NODE,1);
  if (theNode==NULL && vertex_null)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }

        #ifdef TOPNODE
  if (TOPNODE(theVertex) == NULL || LEVEL(TOPNODE(theVertex))<LEVEL(theNode))
    TOPNODE(theVertex) = theNode;
        #endif
  theGrid->status |= 1;

  if (!vertex_null) return(theNode);

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
  return(theNode);
}


/****************************************************************************/
/** \brief Get all nodes related to this element on next level

 * @param   theElement - element for context
 * @param   theElementContext - node context of this element

   This function returns the nodes related to the element on the next
   finer level. The ordering is according to the reference numbering.

   @return <ul>
   <li>   GM_OK    if ok </li>
   <li>   != GM_OK if not ok </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX GetNodeContext (ELEMENT *theElement, NODE **theElementContext)
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
/** \brief Return matching side of the neighboring element

 * @param   theNeighbor - element to test for matching side
 * @param   nbside - the matching side
 * @param   theElement - element with side to match
 * @param   side - side of element to match

   This function computes the matching side of neighboring element.

 */
/****************************************************************************/

void NS_DIM_PREFIX GetNbSideByNodes (ELEMENT *theNeighbor, INT *nbside, ELEMENT *theElement, INT side)
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
/** \brief Return pointer to son edge if it exists
 *
 * @param theEdge edge for which son is searched
 *
 * This function returns the pointer to the son edge if it exists.
 *
 * @return <ul>
   <li>   pointer to specified object </li>
   <li>   NULL if not found </li>
   </ul> */
/****************************************************************************/

EDGE * NS_DIM_PREFIX GetSonEdge (EDGE *theEdge)
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
/** \brief
   GetSonEdges - Return pointer to son edges if it exists

 * @param   theEdge - edge for which son is searched
 * @param   SonEdges - array of pointers will be filled with son edges

   This function returns the pointer to the son edges if existing.

   @return <ul>
   <li>   number of found edges (0,1,2) </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX GetSonEdges (EDGE *theEdge, EDGE *SonEdges[MAX_SON_EDGES])
{
  INT nedges;
  NODE *Node0,*Node1,*SonNode0,*SonNode1,*MidNode;

  nedges = 0;
  SonEdges[0] = NULL;
  SonEdges[1] = NULL;

  Node0 = NBNODE(LINK0(theEdge));
  Node1 = NBNODE(LINK1(theEdge));

  if (GID(Node0)<GID(Node1))
  {
    SonNode0 = SONNODE(Node0);
    SonNode1 = SONNODE(Node1);
  }
  else
  {
    SonNode0 = SONNODE(Node1);
    SonNode1 = SONNODE(Node0);
  }
  MidNode = MIDNODE(theEdge);

  /* parallel note:                                                */
  /* since existance of MidNode decides whether for one SonEdge or */
  /* two half SonEdges is searched, the data structure must be     */
  /* consistent in a way that if the MidNode exists also the       */
  /* MIDNODE pointer is set to MidNode. (s.l. 980227)              */
  if (MidNode == NULL)
  {
    if (SonNode0!=NULL && SonNode1!=NULL)
      SonEdges[0] = GetEdge(SonNode0,SonNode1);
  }
  else
  {
    if (SonNode0!=NULL)
      SonEdges[0] = GetEdge(SonNode0,MidNode);

    if (SonNode1!=NULL)
      SonEdges[1] = GetEdge(MidNode,SonNode1);
  }

  if (SonEdges[0] != NULL) nedges++;
  if (SonEdges[1] != NULL) nedges++;

  return(nedges);
}

/****************************************************************************/
/** \brief
   GetFatherEdge - Return pointer to father edge if it exists

 * @param   theEdge - edge for which father is searched

   This function returns the pointer to the father edge if it exists.

   @return <ul>
   <li>   pointer to specified object </li>
   <li>   NULL if not found </li>
   </ul> */
/****************************************************************************/

EDGE * NS_DIM_PREFIX GetFatherEdge (EDGE *theEdge)
{
  NODE *theNode0 = NBNODE(LINK0(theEdge));
  NODE *theNode1 = NBNODE(LINK1(theEdge));
  EDGE *FatherEdge = NULL;

  /* one node is center node -> no father edge */
  if (CENTERTYPE(theNode0) || CENTERTYPE(theNode1)) return(NULL);

        #ifdef __THREEDIM__
  /* one node is side node -> no father edge */
  if (SIDETYPE(theNode0) || SIDETYPE(theNode1)) return(NULL);
        #endif

  /* both nodes are mid nodes -> no father edge */
  if (MIDTYPE(theNode0) && MIDTYPE(theNode1)) return(NULL);

  /* one node is mid node -> no father edge */
  if (MIDTYPE(theNode0) || MIDTYPE(theNode1))
  {
    NODE *FatherNode0,*FatherNode1, *theNode;

    if (MIDTYPE(theNode1))
    {
      theNode = theNode0; theNode0 = theNode1; theNode1 = theNode;
    }
    FatherEdge = (EDGE *) NFATHER(theNode0);
    if (FatherEdge == NULL) return(NULL);

    FatherNode0 = NBNODE(LINK0(FatherEdge));
    FatherNode1 = NBNODE(LINK1(FatherEdge));
    if (SONNODE(FatherNode0)==theNode1 || SONNODE(FatherNode1)==theNode1)
      return(FatherEdge);
    else
      return(NULL);
  }

  /* both nodes are corner nodes -> try to get the edge */
  if (CORNERTYPE(theNode0) && CORNERTYPE(theNode1))
  {
    if (NFATHER(theNode0)!=NULL && NFATHER(theNode1)!=NULL)
      return(GetEdge(NFATHER(theNode0),NFATHER(theNode1)));
    else
      return(NULL);
  }

  /* one case not considered */
  assert(0);

  return NULL;          /* in case NDEBUG defined */
}

#ifdef __THREEDIM__

/****************************************************************************/
/** \brief
   FatherEdge - Return pointer to father edge if it exists

 * @param   SideNodes - nodes of the side
 * @param   ncorners - number of sidenodes
 * @param   Nodes - corners of edge for which father is searched
 * @param   theEdge - edge for which father is searched

   This function returns the pointer to the father edge if it exists.

   @return <ul>
   <li>   pointer to specified object </li>
   <li>   NULL if not found </li>
   </ul> */
/****************************************************************************/

EDGE * NS_DIM_PREFIX FatherEdge (NODE **SideNodes, INT ncorners, NODE **Nodes, EDGE *theEdge)
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
/** \brief
   GetEdge - Return pointer to edge if it exists

 * @param   from - starting node of edge
 * @param   to - end node of edge

   This function returns the pointer to the specified edge if it exists.

   @return <ul>
   <li>   pointer to specified object </li>
   <li>   NULL if not found </li>
   </ul> */
/****************************************************************************/

EDGE * NS_DIM_PREFIX GetEdge (NODE *from, NODE *to)
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
/** \brief
   CreateEdge - Return pointer to a new edge structure

 * @param   theGrid - grid where vertex should be inserted
 * @param   theElement - pointer to element
 * @param   edge - number of edge
 * @param   with_vector - also create vector for edge (TRUE/FALSE)

   This function returns a pointer to a new edge structure.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

#ifndef ModelP
static
#endif
EDGE *
#ifdef ModelP
NS_DIM_PREFIX
#endif
CreateEdge (GRID *theGrid, ELEMENT *theElement, INT edge, INT with_vector)
{
  ELEMENT *theFather;
  EDGE *pe,*father_edge;
  NODE *from,*to,*n1,*n2;
  LINK *link0,*link1;
  VECTOR *pv;
#ifdef __THREEDIM__
  VERTEX *theVertex;
  NODE *nbn1,*nbn2,*nbn3,*nbn4;
  INT sc,found,side,k,j;
#endif

  from = CORNER(theElement,CORNER_OF_EDGE(theElement,edge,0));
  to = CORNER(theElement,CORNER_OF_EDGE(theElement,edge,1));

  /* check if edge exists already */
  if( (pe = GetEdge(from, to)) != NULL ) {
    if (NO_OF_ELEM(pe)<NO_OF_ELEM_MAX-1)
      INC_NO_OF_ELEM(pe);
    else
      ASSERT(0);

    return(pe);
  }

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    pe = (EDGE*)GetMemoryForObject(theGrid->mg,sizeof(EDGE),EDOBJ);
  else
    pe = (EDGE*)GetMemoryForObject(theGrid->mg,sizeof(EDGE)-sizeof(VECTOR*),EDOBJ);
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
  /* SETPRIO(pe,PrioMaster); */
        #endif
        #ifdef IDENT_ONLY_NEW
  if (GET_IDENT_MODE() == IDENT_ON)
    SETNEW_EDIDENT(pe,1);
        #endif

  UGM_CDBG(pe,
           UserWriteF(PFMT "create edge=" EDID_FMTX " from=" ID_FMTX "tf=%d to=" ID_FMTX "tt=%d"
                      "elem=" EID_FMTX "edge=%d\n",
                      me,EDID_PRTX(pe),ID_PRTX(from),NTYPE(from),ID_PRTX(to),NTYPE(to),
                      EID_PRTX(theElement),edge);
           if (0)
             UserWriteF(PFMT "nfatherf=" ID_FMTX "nfathert=" ID_FMTX " fatheredge=" EDID_FMTX "\n",
                        me,ID_PRTX((NODE*)NFATHER(from)),ID_PRTX((NODE*)NFATHER(to)),
                        EDID_PRTX(GetEdge((NODE*)NFATHER(from),(NODE*)NFATHER(to))));)

  NBNODE(link0) = to;
  NBNODE(link1) = from;
  SET_NO_OF_ELEM(pe,1);
  SETEDGENEW(pe,1);

  /* set edge-subdomain from topological information with respect to father-element */
  SETEDSUBDOM(pe,SUBDOMAIN(theElement));
  theFather = EFATHER(theElement);
  if (theFather!=NULL)
  {
    SETEDSUBDOM(pe,SUBDOMAIN(theFather));
    if (NTYPE(from)<NTYPE(to))
    {
      n1 = from;
      n2 = to;
    }
    else
    {
      n1 = to;
      n2 = from;
    }
    switch(NTYPE(n1)|(NTYPE(n2)<<4))
    {
#ifdef __TWODIM__
    case (CORNER_NODE | (CORNER_NODE<<4)) :
      father_edge = GetEdge(NFATHER(n1),NFATHER(n2));
      if (father_edge!=NULL) SETEDSUBDOM(pe,EDSUBDOM(father_edge));
      break;
    case (CORNER_NODE | (MID_NODE<<4)) :
      father_edge = NFATHEREDGE(n2);
#ifdef ModelP
      if (father_edge==NULL)
      {
        /* TODO: check this after priority set:
           assert( GHOST(n1) || GHOST(n2) ); */
        break;
      }
#endif
      assert(father_edge!=NULL);
      if (NBNODE(LINK0(father_edge))==NFATHER(n1) || NBNODE(LINK1(father_edge))==NFATHER(n1)) SETEDSUBDOM(pe,EDSUBDOM(father_edge));
      break;
#endif
#ifdef __THREEDIM__
    case (CORNER_NODE | (CORNER_NODE<<4)) :
      father_edge = GetEdge(NFATHER(n1),NFATHER(n2));
      if (father_edge!=NULL) SETEDSUBDOM(pe,EDSUBDOM(father_edge));
      else
      {
        /* do fathers of n1, n2 lies on a side (of the father) which has BNDS? */
        for (j=0; j<SIDES_OF_ELEM(theFather); j++)
        {
          found=0;
          for (k=0; k<CORNERS_OF_SIDE(theFather,j); k++)
          {
            sc = CORNER_OF_SIDE(theFather,j,k);
            if (CORNER(theFather,sc)==NFATHER(n1) || CORNER(theFather,sc)==NFATHER(n2)) found++;
          }
          if (found==2 && (OBJT(theFather)==BEOBJ) && SIDE_ON_BND(theFather,j))
          {
            SETEDSUBDOM(pe,0);
            break;
          }
        }
      }
      break;

    case (CORNER_NODE | (MID_NODE<<4)) :
      father_edge = NFATHEREDGE(n2);
      assert(father_edge!=NULL);
      nbn1 = NBNODE(LINK0(father_edge));
      nbn2 = NBNODE(LINK1(father_edge));
      if (nbn1==NFATHER(n1) || nbn2==NFATHER(n1)) SETEDSUBDOM(pe,EDSUBDOM(father_edge));
      else
      {
        /* do all nodes n1, nbn1, nbn2 ly on the same side of father? */
        side=-1;
        for (j=0; j<SIDES_OF_ELEM(theFather); j++)
        {
          found=0;
          for (k=0; k<CORNERS_OF_SIDE(theFather,j); k++)
          {
            sc = CORNER_OF_SIDE(theFather,j,k);
            if (CORNER(theFather,sc)==NFATHER(n1) || CORNER(theFather,sc)==nbn1 || CORNER(theFather,sc)==nbn2) found++;
          }
          if (found==3)
          {
            side = j;
            break;
          }
        }
        if (side>=0  && (OBJT(theFather)==BEOBJ) && SIDE_ON_BND(theFather,side)) SETEDSUBDOM(pe,0);
      }
      break;

    case (MID_NODE | (MID_NODE<<4)) :
      father_edge = NFATHEREDGE(n1);
      assert(father_edge!=NULL);
      nbn1 = NBNODE(LINK0(father_edge));
      nbn2 = NBNODE(LINK1(father_edge));
      father_edge = NFATHEREDGE(n2);
      assert(father_edge!=NULL);
      nbn3 = NBNODE(LINK0(father_edge));
      nbn4 = NBNODE(LINK1(father_edge));

      /* do all nodes nbn1, nbn2, nbn3, nbn4 ly on the same side of father? */
      side=-1;
      for (j=0; j<SIDES_OF_ELEM(theFather); j++)
      {
        found=0;
        for (k=0; k<CORNERS_OF_SIDE(theFather,j); k++)
        {
          sc = CORNER_OF_SIDE(theFather,j,k);
          if (CORNER(theFather,sc)==nbn1) found++;
          if (CORNER(theFather,sc)==nbn2) found++;
          if (CORNER(theFather,sc)==nbn3) found++;
          if (CORNER(theFather,sc)==nbn4) found++;
        }
        if (found==4)
        {
          side = j;
          break;
        }
      }
      if (side>=0 && (OBJT(theFather)==BEOBJ) && SIDE_ON_BND(theFather,side)) SETEDSUBDOM(pe,0);
      break;

    case (CORNER_NODE | (SIDE_NODE<<4)) :
      theVertex = MYVERTEX(n2);
      if (VFATHER(theVertex) == theFather)
        side = ONSIDE(theVertex);
      else
        side = ONNBSIDE(theVertex);
      if ((OBJT(theFather)==BEOBJ) && SIDE_ON_BND(theFather,side))
        for (k=0; k<CORNERS_OF_SIDE(theFather,side); k++)
          if (CORNER(theFather,CORNER_OF_SIDE(theFather,side,k))==NFATHER(n1))
          {
            SETEDSUBDOM(pe,0);
            break;
          }
      break;

    case (MID_NODE | (SIDE_NODE<<4)) :
      theVertex = MYVERTEX(n2);
      if (VFATHER(theVertex) == theFather)
        side = ONSIDE(theVertex);
      else
        side = ONNBSIDE(theVertex);
      if ((OBJT(theFather)==BEOBJ) && SIDE_ON_BND(theFather,side))
      {
        found=0;
        father_edge = NFATHEREDGE(n1);
        assert(father_edge!=NULL);
        nbn1 = NBNODE(LINK0(father_edge));
        nbn2 = NBNODE(LINK1(father_edge));
        for (k=0; k<CORNERS_OF_SIDE(theFather,side); k++)
        {
          if (CORNER(theFather,CORNER_OF_SIDE(theFather,side,k))==nbn1 || CORNER(theFather,CORNER_OF_SIDE(theFather,side,k))==nbn2)
            found++;
        }
        if (found==2) SETEDSUBDOM(pe,0);
      }
      break;
#endif
    }     /* end switch */
  }   /* end (theFather!=NULL) */

  /* create vector if */
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    if (with_vector) {
      if (CreateVector (theGrid,EDGEVEC,(GEOM_OBJECT *)pe,&pv)) {
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
  NE(theGrid)++;

  /* return ok */
  return(pe);
}

/****************************************************************************/
/** \brief
   GetLink - Return pointer to link if it exists

 * @param   from - starting node of link
 * @param   to - end node of link

   This function returns the pointer to the specified link if it exists.

   @return <ul>
   <li>   pointer to specified link </li>
   <li>   NULL if not found. </li>
   </ul> */
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
/** \brief
   CreateElement - Return a pointer to  a new element structure

 * @param   theGrid - grid structure to extend
 * @param   tag - the element type
 * @param   objtype - inner element (IEOBJ) or boundary element (BEOBJ)
 * @param   nodes - list of corner nodes in reference numbering
 * @param   Father - pointer to father element (NULL on base level)
 * @param   with_vector -

   This function creates and initializes a new element and returns a pointer to it.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX CreateElement (GRID *theGrid, INT tag, INT objtype, NODE **nodes,
                                       ELEMENT *Father, INT with_vector)
{
  ELEMENT *pe;
  INT i,s_id;
  VECTOR *pv;
  void *q;

  if (objtype == IEOBJ)
    pe = (ELEMENT*)GetMemoryForObject(MYMG(theGrid),INNER_SIZE_TAG(tag),
                                      MAPPED_INNER_OBJT_TAG(tag));
  else if (objtype == BEOBJ)
    pe = (ELEMENT*)GetMemoryForObject(MYMG(theGrid),BND_SIZE_TAG(tag),
                                      MAPPED_BND_OBJT_TAG(tag));

  if (pe==NULL) return(NULL);

  /* initialize data */
  SETNEWEL(pe,1);
  SETOBJT(pe,objtype);
  SETTAG(pe,tag);
  SETLEVEL(pe,theGrid->level);
        #ifdef ModelP
  DDD_AttrSet(PARHDRE(pe),GRID_ATTR(theGrid));
  /* SETEPRIO(pe,PrioMaster); */
  PARTITION(pe) = me;
        #endif
  SETEBUILDCON(pe,1);
  ID(pe) = (theGrid->mg->elemIdCounter)++;

  /* subdomain id */
  s_id = (Father != NULL) ? SUBDOMAIN(Father) : 0;
  SETSUBDOMAIN(pe,s_id);

        #ifdef __CENTERNODE__
  SET_CENTERNODE(pe,NULL);
        #endif

  SET_EFATHER(pe,Father);

  /* set corner nodes */
  for (i=0; i<CORNERS_OF_ELEM(pe); i++)
    SET_CORNER(pe,i,nodes[i]);

  /* create edges */
  for (i=0; i<EDGES_OF_ELEM(pe); i++)
    if (CreateEdge (theGrid,pe,i,with_vector) == NULL) {
      DisposeElement(theGrid,pe,TRUE);
      return(NULL);
    }

  UGM_CDBG(pe,
           UserWriteF(PFMT "create elem=" EID_FMTX,
                      me,EID_PRTX(pe));
           for (i=0; i<CORNERS_OF_ELEM(pe); i++)
             UserWriteF(" n%d=" ID_FMTX, i,ID_PRTX(CORNER(pe,i)));
           UserWriteF("\n");
           for (i=0; i<EDGES_OF_ELEM(pe); i++)
           {
             EDGE *theEdge;

             theEdge = GetEdge(CORNER_OF_EDGE_PTR(pe,i,0),
                               CORNER_OF_EDGE_PTR(pe,i,1));
             UserWriteF(" e%d=" EDID_FMTX, i,EDID_PRTX(theEdge));
           }
           UserWriteF("\n");)


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

  if (me == -1) {
    assert(KeyForObject((KEY_OBJECT *)pe) != -66529);
  }

  /* return ok */
  return(pe);
}

/****************************************************************************/
/** \brief
   CreateSonElementSide - creates the element sides of son elements

 * @param   theGrid - grid for which to create
 * @param   theElement - pointer to a boundary element
 * @param   side - side id of a side of the element
 * @param   theSon - pointer to a son element
 * @param   son_side - side id of a side of the son

   This function creates and initializes an element side of a son element.
   Here also the side vector (iff at all) is inspected in 'ReinspectSonSideVector'.
   The latter function eventually reallocates the vector if its size has changed and
   sets the VBUILDCON flag in the vector. The connections of the old vector are
   thereby disposed. The refine-module which is calling
   'CreateSonElementSide' will finally call 'GridCreateConnection' to reinstall
   the connections of the side-vector.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX CreateSonElementSide (GRID *theGrid, ELEMENT *theElement, INT side,
                                        ELEMENT *theSon, INT son_side)
{
  INT n,i;
  BNDS *bnds;
  BNDP *bndp[MAX_CORNERS_OF_ELEM];
  VECTOR *vec;
  EDGE *theEdge;

  ASSERT (OBJT(theElement) == BEOBJ);

  ASSERT (ELEM_BNDS(theElement,side) != NULL);

  /* check if Edges of theElement, which are on the side 'side' have EDSUBDOM 0 */
  n = CORNERS_OF_SIDE(theElement,side);
  for (i=0; i<n; i++)
  {
    theEdge = GetEdge(CORNER(theElement,CORNER_OF_SIDE(theElement,side,i)),CORNER(theElement,CORNER_OF_SIDE(theElement,side,(i+1)%n)));
    assert(EDSUBDOM(theEdge)==0);
  }

  n = CORNERS_OF_SIDE(theSon,son_side);
  for (i=0; i<n; i++)
  {
    /* check if vertices of Son ly on boundary */
    if (OBJT(MYVERTEX(CORNER(theSon,CORNER_OF_SIDE(theSon,son_side,i))))!=BVOBJ)
    {
      NODE *theNode,*NFather;
      EDGE *theFatherEdge;
      INT t1,t2;

      theNode = CORNER(theSon,CORNER_OF_SIDE(theSon,son_side,i));
      printf("ID=%d\n",(int)ID(theNode));
      switch (NTYPE(theNode))
      {
      case CORNER_NODE :
        printf("NTYPE = CORNER_NODE");
        NFather = NFATHER(theNode);
        break;
      case MID_NODE :
        printf(PFMT "el "EID_FMTX " son "EID_FMTX " vertex "VID_FMTX "\n",me,EID_PRTX(theElement),EID_PRTX(theSon),VID_PRTX(MYVERTEX(CORNER(theSon,CORNER_OF_SIDE(theSon,son_side,i)))));
        printf(PFMT "NTYPE = MID_NODE\n",me);
        theFatherEdge = NFATHEREDGE(theNode);
        printf(PFMT "EDSUBDOM = %d\n",me,(int)EDSUBDOM(theFatherEdge));
        t1 = (OBJT(MYVERTEX(NBNODE(LINK0(theFatherEdge))))==BVOBJ);
        t2 = (OBJT(MYVERTEX(NBNODE(LINK1(theFatherEdge))))==BVOBJ);
        printf(PFMT "BVOBJ(theFatherEdge): %d %d\n",me,(int)t1,(int)t2);
        break;
      case SIDE_NODE :
        printf("NTYPE = SIDE_NODE");
        break;
      case CENTER_NODE :
        printf("NTYPE = CENTER_NODE");
        break;
      }
      ASSERT(0);
    }
    bndp[i] = V_BNDP(MYVERTEX(CORNER(theSon,CORNER_OF_SIDE(theSon,son_side,i))));
  }
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
    #ifdef __TWODIM__
  theEdge = GetEdge(CORNER(theSon,CORNER_OF_EDGE(theSon,son_side,0)),
                    CORNER(theSon,CORNER_OF_EDGE(theSon,son_side,1)));
  ASSERT(theEdge != NULL);
  SETEDSUBDOM(theEdge,0);
        #endif

    #ifdef __THREEDIM__
  /** \todo is this necessary?
     for (i=0; i<EDGES_OF_SIDE(theSon,son_side); i++) {
          int k  = EDGE_OF_SIDE(theSon,son_side,i);
          theEdge = GetEdge(CORNER(theSon,CORNER_OF_EDGE(theSon,k,0)),
                                            CORNER(theSon,CORNER_OF_EDGE(theSon,k,1)));
          ASSERT(theEdge != NULL);
          SETEDSUBDOM(theEdge,0);
     } */
        #endif

  return(GM_OK);
}

/****************************************************************************/
/** \brief
   CreateNewLevel - Return pointer to new grid structure

 * @param   theMG - multigrid structure

   This function creates and initialized a new grid structure for top level + 1
   and returns a pointer to it.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

GRID * NS_DIM_PREFIX CreateNewLevel (MULTIGRID *theMG, INT algebraic)
{
  GRID *theGrid;
  INT l;

  if (BOTTOMLEVEL(theMG)>TOPLEVEL(theMG) && algebraic) return (NULL);
  if (TOPLEVEL(theMG)+1>=MAXLEVEL) return(NULL);
  if (algebraic) l = BOTTOMLEVEL(theMG)-1;
  else l = TOPLEVEL(theMG)+1;

  /* allocate grid object */
  theGrid = (GRID*)GetMemoryForObject(theMG,sizeof(GRID),GROBJ);
  if (theGrid==NULL) return(NULL);

  /* fill in data */
  CTRL(theGrid) = 0;
  SETOBJT(theGrid,GROBJ);
  GLEVEL(theGrid) = l;
  GATTR(theGrid) = GRID_ATTR(theGrid);
  NE(theGrid) = 0;
  NC(theGrid) = 0;
  /* other counters are init in INIT fcts below */

#ifdef __INTERPOLATION_MATRIX__
  NIMAT(theGrid) = 0;
#endif

  GSTATUS(theGrid,0);
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
  MYMG(theGrid) = theMG;
  GRID_ON_LEVEL(theMG,l) = theGrid;
  if (algebraic) BOTTOMLEVEL(theMG) = l;
  else
  {
    TOPLEVEL(theMG) = l;
    CURRENTLEVEL(theMG) = l;
  }

  return(theGrid);
}


/****************************************************************************/
/** \brief
   CreateNewLevelAMG - Create new amg level

 * @param   theMG - multigrid structure

   This function creates and initialized a new grid structure for bottomLevel - 1
   and returns a pointer to it.

   @return <ul>
   <li>   pointer to requested object </li>
   <li>   NULL if out of memory </li>
   </ul> */
/****************************************************************************/

GRID * NS_DIM_PREFIX CreateNewLevelAMG (MULTIGRID *theMG)
{
  GRID *theGrid;
  int l;

  if (theMG->bottomLevel-1<=-MAXLEVEL) return(NULL);

  l = theMG->bottomLevel-1;

  /* allocate grid object */
  theGrid = (GRID*)GetMemoryForObject(theMG,sizeof(GRID),GROBJ);
  if (theGrid==NULL) return(NULL);

  /* fill in data */
  CTRL(theGrid) = 0;
  SETOBJT(theGrid,GROBJ);
  theGrid->level = l;
  NE(theGrid) = 0;
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
/** \brief
   MakeMGItem - Create a multigrid environment item

 * @param   name - name of the multigrid

   This function creates a multigrid environment directory.

   @return <ul>
   <li>   pointer to new MULTIGRID </li>
   <li>   NULL if error occured </li>
   </ul> */
/****************************************************************************/

MULTIGRID * NS_DIM_PREFIX MakeMGItem (const char *name)
{
  MULTIGRID *theMG;

  if (ChangeEnvDir("/Multigrids") == NULL) return (NULL);
  if (strlen(name)>=NAMESIZE || strlen(name)<=1) return (NULL);
  theMG = (MULTIGRID *) MakeEnvItem(name,theMGDirID,sizeof(MULTIGRID));
  if (theMG == NULL) return(NULL);

  return (theMG);
}

/****************************************************************************/
/** \todo Please doc me!

 * @param   theMG
 * @param   FromLevel
 * @param   ToLevel
 * @param   mask

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

INT NS_DIM_PREFIX ClearMultiGridUsedFlags (MULTIGRID *theMG, INT FromLevel, INT ToLevel, INT mask)
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
/** \brief
   GetMultigrid - Find the multigrid environment item with name

 * @param   name - name of the multigrid to find

   This function find the multigrid environment item with `name` and
   returns a pointer to the multigrid structure.

   @return <ul>
   <li>   pointer to MULTIGRID  </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

MULTIGRID * NS_DIM_PREFIX GetMultigrid (const char *name)
{
  return ((MULTIGRID *) SearchEnv(name,"/Multigrids",
                                  theMGDirID,theMGRootDirID));
}

/****************************************************************************/
/** \brief
   GetFirstMultigrid - Return a pointer to the first multigrid

 * @param   void

   This function returns a pointer to the first multigrid in the /Multigrids
   directory.

   @return <ul>
   <li>   pointer to MULTIGRID </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

MULTIGRID * NS_DIM_PREFIX GetFirstMultigrid ()
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
/** \brief Return a pointer to the next multigrid

 * @param   theMG - multigrid structure

   This function returns a pointer to the next multigrid in the /Multigrids
   directory.

   @return <ul>
   <li>   pointer to MULTIGRID </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

MULTIGRID * NS_DIM_PREFIX GetNextMultigrid (const MULTIGRID *theMG)
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
/** \brief
   CreateMultiGrid - Return a pointer to new multigrid structure

 * @param   MultigridName - name of multigrid
 * @param   domain - name of domain description from environment
 * @param   problem - name of problem description from environment
 * @param   format - name of format description from environment
 * @param   heapSize - size of heap to allocate for that multigrid in bytes
 * @param   optimizedIE - alloccate NodeElementList

   This function creates and initializes a new multigrid structure including
   allocation of heap, combining the domain and the boundary conditions
   and creation of the fixed corners of the domain.

   @return <ul>
   <li>   pointer to new object </li>
   <li>   NULL if an error occured. </li>
   </ul> */
/****************************************************************************/

MULTIGRID * NS_DIM_PREFIX CreateMultiGrid (char *MultigridName, char *BndValProblem,
                                           char *format, MEM heapSize, INT optimizedIE, INT insertMesh)
{
  HEAP *theHeap,*theUserHeap;
  MULTIGRID *theMG;
  GRID *theGrid;
  INT i,ds;
  BVP *theBVP;
  BVP_DESC *theBVPDesc;
  MESH mesh;
  FORMAT *theFormat;
  INT MarkKey;

  theFormat = GetFormat(format);
  if (theFormat==NULL)
  {
    PrintErrorMessage('E',"CreateMultiGrid","format not found");
    return(NULL);
  }


  /* allocate multigrid envitem */
  theMG = MakeMGItem(MultigridName);
  if (theMG==NULL) return(NULL);
  MGFORMAT(theMG) = theFormat;
  if (InitElementTypes(theMG)!=GM_OK)
  {
    PrintErrorMessage('E',"CreateMultiGrid","error in InitElementTypes");
    return(NULL);
  }

  /* allocate the heap */
  theHeap = NewHeap(SIMPLE_HEAP, heapSize, malloc(heapSize));
  if (theHeap==NULL)
  {
    UserWriteF("CreateMultiGrid: cannot allocate %ld bytes\n", heapSize);
    PrintErrorMessage('E', "CreateMultiGrid","Cannot allocate heap!");

    DisposeMultiGrid(theMG);
    return(NULL);
  }

  /* mark temp memory here, release it after coarse grid construction in FixCoarseGrid */
  MarkTmpMem(theHeap,&MarkKey);
  MG_MARK_KEY(theMG) = MarkKey;

  if (insertMesh)
    theBVP = BVP_Init(BndValProblem,theHeap,&mesh,MarkKey);
  else
    theBVP = BVP_Init(BndValProblem,theHeap,NULL,MarkKey);
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
        #ifdef DYNAMIC_MEMORY_ALLOCMODEL
  theMG->bottomtmpmem = 0;
        #endif
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
           (ELEMENT***)GetTmpMem(theHeap,NDELEM_BLKS_MAX*sizeof(ELEMENT**),MarkKey))==NULL)
    {
      ReleaseTmpMem(theHeap,MarkKey);
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
  if (insertMesh)
  {
                #ifdef ModelP
    if (me==master)
    {
                #endif
    if (InsertMesh(theMG,&mesh))
    {
      DisposeMultiGrid(theMG);
      return(NULL);
    }
                #ifdef ModelP
  }
                #endif

    ASSERT(mesh.mesh_status!=MESHSTAT_NOTINIT);
    if (mesh.mesh_status==MESHSTAT_MESH)
      if (FixCoarseGrid(theMG))
      {
        DisposeMultiGrid(theMG);
        return(NULL);
      }
  }
  /* return ok */
  return(theMG);
}

/****************************************************************************/
/** \brief Remove edge from the data structure

 * @param   theGrid - grid to remove from
 * @param   theEdge - edge to remove

   This function remove an edge from the data structure including its
   vector (if one) and inserts them into the free list.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 if an error occured. </li>
   </ul> */
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
  NE(theGrid)--;
  return(0);
}

/****************************************************************************/
/** \brief Remove node including its edges from the data structure

 * @param   theGrid - grid to remove from
 * @param   theNode - node to remove

   This function removes node including its edges and vector (if one)
   from the data structure and inserts all objects into the free list.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DisposeNode (GRID *theGrid, NODE *theNode)
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
  father = (GEOM_OBJECT *)NFATHER(theNode);
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

                        #ifdef __CENTERNODE__
    case (CENTER_NODE) :
      ASSERT(OBJT(father)==IEOBJ || OBJT(father)==BEOBJ);
      SET_CENTERNODE((ELEMENT *)father,NULL);
      break;
                        #endif

    default :
      ASSERT(0);
      break;
    }
  }

  /** \todo delete old vertex handling */
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
#ifdef __PERIODIC_BOUNDARY__
    if (theNode == (NODE *)VOBJECT(NVECTOR(theNode)))
      VOBJECT(NVECTOR(theNode)) = NULL;
#endif
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
/** \brief
   DisposeVertex - Remove vertex from the data structure

 * @param   theGrid - grid to remove from
 * @param   theVertex - vertex to remove

   This function removes a vertex from the data structure
   and puts it into the free list.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 no valid object number </li>
   </ul> */
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
/** \brief
   DisposeElement - Remove element from the data structure

 * @param   theGrid - grid to remove from
 * @param   theElement - element to remove
 * @param   dispose_connections - also dispose connections (TRUE/FALSE)

   This function removes an element from the data structure and inserts it
   into the free list. This includes all elementsides, sidevectors and the
   elementvector if they exist.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 no valid object number. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DisposeElement (GRID *theGrid, ELEMENT *theElement, INT dispose_connections)
{
  INT i,j,tag;
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
        #ifndef ModelP
  INT edge;
        #endif

  HEAPFAULT(theElement);

  GRID_UNLINK_ELEMENT(theGrid,theElement);

        #ifdef __CENTERNODE__
  {
    theNode = CENTERNODE(theElement);

    if (theNode != NULL) SETNFATHER(theNode,NULL);
  }
        #endif

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
  /** \todo possibly some son cannot be reached by GetAllSons, */
  /* because their father has not been on this proc and       */
  /* they lost their father pointers                          */
  if (NSONS(theElement)>0)
  {
    INT i,j;
    ELEMENT *SonList[MAX_SONS];

    if (GetAllSons(theElement,SonList)) RETURN(GM_FATAL);

    i = 0;
    while (SonList[i] != NULL)
    {
      PRINTDEBUG(gm,2,(PFMT "DisposeElement(): elem=" EID_FMTX
                       " deleting fatherpointer of son=" EID_FMTX "\n",
                       me,EID_PRTX(theElement),EID_PRTX(SonList[i])));
      SET_EFATHER(SonList[i],NULL);

      /* reset VFATHER of centernode vertex */
      for (j=0; j<CORNERS_OF_ELEM(SonList[i]); j++)
      {
        theNode = CORNER(SonList[i],j);
                                #ifndef __CENTERNODE__
        if (CENTERTYPE(theNode) && NFATHER(theNode)!=NULL)
          SETNFATHER(theNode,NULL);
                                #endif
        theVertex = MYVERTEX(theNode);
        if (VFATHER(theVertex) != NULL && VFATHER(theVertex) == theElement)
          VFATHER(theVertex) = NULL;
      }
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
#ifndef __SWAPBYTES__  /* Don't check for bnds HEAPFAULTs on little endian machines!! */
        HEAPFAULT(bnds);
#endif
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

#ifdef __OVERLAP2__
    if( ce_NO_DELETE_OVERLAP2 != -1 && NO_DELETE_OVERLAP2(theNode) )
    {
      continue;
    }
#endif

    if (START(theNode) == NULL)
    {
      if (NTYPE(theNode)==MID_NODE)
      {
        if (NFATHER(theNode)!=NULL)
        {
          MIDNODE((EDGE *)NFATHER(theNode)) = NULL;
        }
                #ifndef ModelP
        /* HEAPFAULT in theFather possible, if in a previous call
           some son is not reached by GetAllSons */
        else
        {
          theVertex = MYVERTEX(theNode);
          theFather = VFATHER(theVertex);
          if (theFather != NULL)
          {
            INT edge = ONEDGE(theVertex);
            theEdge = GetEdge(CORNER(theFather,
                                     CORNER_OF_EDGE(theFather,edge,0)),
                              CORNER(theFather,
                                     CORNER_OF_EDGE(theFather,edge,1)));
            ASSERT(theEdge!=NULL);
            MIDNODE(theEdge) = NULL;
          }
        }
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
  /** \todo delete */
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
          VOBJECT(theVector) = (GEOM_OBJECT *)theNeighbor;
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

#ifndef ModelP
#define DO_NOT_DISPOSE  return (2)
#else
#define DO_NOT_DISPOSE  dispose=0
#endif

/****************************************************************************/
/** \brief
   Collapse - construct coarse grid from surface

 * @param   theMG - multigrid to collapse

   This function constructs coarse grid from surface. ATTENTION: Use refine $g
   to cover always the whole domain with the grid on each level.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 no valid object number </li>
   <li>   2 grid structure not empty or level 0 </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX Collapse (MULTIGRID *theMG)
{
  GRID *theGrid;
  ELEMENT *theElement;
  NODE *theNode;
  EDGE *theEdge;
  VERTEX *theVertex;
  VECTOR *vec;
  INT tl = TOPLEVEL(theMG);
  INT l,i;

        #ifdef DYNAMIC_MEMORY_ALLOCMODEL
  if (MG_COARSE_FIXED(theMG))
    if (DisposeBottomHeapTmpMemory(theMG))
      REP_ERR_RETURN(1);
        #endif

  if( DisposeAMGLevels(theMG) )
    REP_ERR_RETURN(1);

#ifdef ModelP
  DDD_XferBegin();
    #ifdef DDDOBJMGR
  DDD_ObjMgrBegin();
    #endif
#endif
  for (l=tl-1; l>=0; l--) {
    theGrid = GRID_ON_LEVEL(theMG,l);
#ifdef __PERIODIC_BOUNDARY__
    GridSetPerVecCount(theGrid);
#endif
    for (theNode=PFIRSTNODE(theGrid); theNode != NULL;
         theNode = SUCCN(theNode)) {
      SONNODE(theNode) = NULL;
      SETNFATHER(theNode,NULL);
    }
    for (theElement=PFIRSTELEMENT(theGrid); theElement != NULL;
         theElement = SUCCE(theElement)) {
      SETNSONS(theElement,0);
      SET_SON(theElement,0,NULL);
                #ifdef ModelP
      SET_SON(theElement,1,NULL);
                #endif
      for (i=0; i<EDGES_OF_ELEM(theElement); i++) {
        theEdge = GetEdge(CORNER(theElement,
                                 CORNER_OF_EDGE(theElement,i,0)),
                          CORNER(theElement,
                                 CORNER_OF_EDGE(theElement,i,1)));
        MIDNODE(theEdge) = NULL;
      }
    }
    while (PFIRSTELEMENT(theGrid)!=NULL)
      if (DisposeElement(theGrid,PFIRSTELEMENT(theGrid),1))
        return(1);
    while (PFIRSTNODE(theGrid)!=NULL)
    {
      if (DisposeNode(theGrid,PFIRSTNODE(theGrid)))
        return(1);
    }
    while (PFIRSTVERTEX(theGrid)!=NULL) {
      theVertex = PFIRSTVERTEX(theGrid);
      GRID_UNLINK_VERTEX(theGrid,theVertex);
      GRID_LINK_VERTEX(GRID_ON_LEVEL(theMG,tl),
                       theVertex,VXPRIO(theVertex));
    }
    GRID_ON_LEVEL(theMG,l) = NULL;
  }

#ifdef ModelP
    #ifdef DDDOBJMGR
  DDD_ObjMgrEnd();
    #endif
  DDD_XferEnd();
#endif

  /* move top level grid to bottom (level 0) */
  theGrid = GRID_ON_LEVEL(theMG,tl);
  theGrid->finer = NULL;
  theGrid->coarser = NULL;
  theGrid->level = 0;
  GATTR(theGrid) = GRID_ATTR(theGrid);
  GRID_ON_LEVEL(theMG,tl) = NULL;
  GRID_ON_LEVEL(theMG,0) = theGrid;
  theMG->topLevel = 0;
  theMG->fullrefineLevel = 0;
  theMG->currentLevel = 0;

  for (theNode=PFIRSTNODE(theGrid); theNode != NULL;
       theNode = SUCCN(theNode)) {
    SETNFATHER(theNode,NULL);
    SETNTYPE(theNode,LEVEL_0_NODE);
    SETNCLASS(theNode,3);
    SETNNCLASS(theNode,0);
    SETLEVEL(theNode,0);
    VFATHER(MYVERTEX(theNode)) = NULL;
                        #ifdef ModelP
    DDD_AttrSet(PARHDR(theNode),GRID_ATTR(theGrid));
                        #endif
  }
  for (theElement=PFIRSTELEMENT(theGrid); theElement != NULL;
       theElement = SUCCE(theElement)) {
    SETECLASS(theElement,RED_CLASS);
    SET_EFATHER(theElement,NULL);
    SETLEVEL(theElement,0);
                #ifdef ModelP
    DDD_AttrSet(PARHDRE(theElement),GRID_ATTR(theGrid));
                #endif
    for (i=0; i<EDGES_OF_ELEM(theElement); i++) {
      theEdge = GetEdge(CORNER(theElement,
                               CORNER_OF_EDGE(theElement,i,0)),
                        CORNER(theElement,
                               CORNER_OF_EDGE(theElement,i,1)));
      SETLEVEL(theEdge,0);
                        #if (defined ModelP) && (defined __THREEDIM__)
      DDD_AttrSet(PARHDR(theEdge),GRID_ATTR(theGrid));
                        #endif
    }
  }
  for (theVertex=PFIRSTVERTEX(theGrid); theVertex != NULL;
       theVertex = SUCCV(theVertex)) {
    SETLEVEL(theVertex,0);
                #ifdef ModelP
    DDD_AttrSet(PARHDRV(theVertex),GRID_ATTR(theGrid));
                #endif
    ASSERT(NOOFNODE(theVertex)==1);
  }

        #ifdef ModelP
  for (vec=PFIRSTVECTOR(theGrid); vec != NULL; vec = SUCCVC(vec))
    DDD_AttrSet(PARHDR(vec),GRID_ATTR(theGrid));

  /* rebiuld all DDD interfaces due to removed objects and changed attributes */
  DDD_IFRefreshAll();
        #endif

        #ifdef DYNAMIC_MEMORY_ALLOCMODEL
  if (MG_COARSE_FIXED(theMG))
    if (CreateAlgebra(theMG))
      REP_ERR_RETURN(1);
        #endif

  return(0);
}

/****************************************************************************/
/** \brief
   DisposeTopLevel - Remove top level grid from multigrid  structure

 * @param   theMG - multigrid to remove from

   This function removes the top level grid from multigrid structure.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 no valid object number </li>
   <li>   2 grid structure not empty or level 0 </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DisposeTopLevel (MULTIGRID *theMG)
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
  if (PFIRSTELEMENT(theGrid)!=NULL) DO_NOT_DISPOSE;
  if (PFIRSTVERTEX(theGrid)!=NULL) DO_NOT_DISPOSE;
  if (PFIRSTNODE(theGrid)!=NULL) DO_NOT_DISPOSE;

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
/** \brief
   DisposeGrid - dispose top level grid

 * @param   theGrid - grid to be removed

   This function removes the top level grid from multigrid structure.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 no valid object number </li>
   <li>   2 grid structure not empty or level 0 </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DisposeGrid (GRID *theGrid)
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
  while (PFIRSTELEMENT(theGrid)!=NULL)
    if (DisposeElement(theGrid,PFIRSTELEMENT(theGrid),1))
      return(2);

  while (PFIRSTNODE(theGrid)!=NULL)
    if (DisposeNode(theGrid,PFIRSTNODE(theGrid)))
      return(2);

  while (PFIRSTVERTEX(theGrid)!=NULL)
    if (DisposeVertex(theGrid,PFIRSTVERTEX(theGrid)))
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
/** \brief
   DisposeAMGLevel - dispose bottom AMG level

 * @param   theMG - multigrid to remove from

   This function removes the bottom AMG level from multigrid structure.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 no valid object number </li>
   <li>   2 no AMG levels </li>
   </ul> */
/****************************************************************************/

static INT DisposeAMGLevel (MULTIGRID *theMG)
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
  while (PFIRSTVECTOR(theGrid)!=NULL)
  {
    /* In ModelP, the DisposeVector is done on all procs which
           own copies. We do it without Xfer-communication. */
    if (DisposeVector(theGrid,PFIRSTVECTOR(theGrid)))
      return(1);
  }

        #ifdef ModelP
  /* stop dangerous mode. from now on DDD will issue warnings again. */
  DDD_SetOption(OPT_WARNING_DESTRUCT_HDR, OPT_ON);
        #endif

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
/** \brief
   DisposeAMGLevels - dispose all AMG level

 * @param   theMG - multigrid to remove from

   This function removes all AMG level from multigrid structure.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 if error occured </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DisposeAMGLevels (MULTIGRID *theMG)
{
  INT err;

        #ifdef ModelP
  /* tell DDD that we will 'inconsistently' delete objects.
     this is a dangerous mode as it switches DDD warnings off. */
  /** \briefDD_SetOption(OPT_WARNING_DESTRUCT_HDR, OPT_OFF);*/
  DDD_XferBegin();
    #ifdef DDDOBJMGR
  DDD_ObjMgrBegin();
        #endif
        #endif

  while ((err=DisposeAMGLevel(theMG))!=2)
    if (err==1)
    {
      PrintErrorMessage('E',"AMGTransferPreProcess","could not dispose AMG levels");
      REP_ERR_RETURN(1);
    }

        #ifdef ModelP
  /* stop dangerous mode. from now on DDD will issue warnings again. */
  /*DDD_SetOption(OPT_WARNING_DESTRUCT_HDR, OPT_ON);*/

  /* rebuild DDD-interfaces because distributed vectors have been
     deleted without communication */
  /*DDD_IFRefreshAll();*/
    #ifdef DDDOBJMGR
  DDD_ObjMgrEnd();
        #endif
  DDD_XferEnd();
        #endif

  return(0);
}


/****************************************************************************/
/** \brief Release memory for the whole multigrid  structure

 * @param   theMG - multigrid to remove

   This function releases the memory for the whole multigrid  structure.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DisposeMultiGrid (MULTIGRID *theMG)
{
  INT level;

        #ifdef DYNAMIC_MEMORY_ALLOCMODEL
  if (DisposeBottomHeapTmpMemory(theMG)) REP_ERR_RETURN(1);
        #else
  if (DisposeAMGLevels(theMG))
    RETURN(1);
        #endif

        #ifdef ModelP
  /* tell DDD that we will 'inconsistently' delete objects.
     this is a dangerous mode as it switches DDD warnings off. */
  DDD_SetOption(OPT_WARNING_DESTRUCT_HDR, OPT_OFF);
        #endif

  for (level = TOPLEVEL(theMG); level >= 0; level --)
    if (DisposeGrid(GRID_ON_LEVEL(theMG,level)))
      RETURN(1);

        #ifdef ModelP
  /* stop dangerous mode. from now on DDD will issue warnings again. */
  DDD_SetOption(OPT_WARNING_DESTRUCT_HDR, OPT_ON);

  /* rebuild DDD-interfaces because distributed vectors have been
     deleted without communication */
  DDD_IFRefreshAll();
        #endif

  /** \todo Normally the MG-heap should be cleaned-up before freeing.
           DDD depends on storage in the heap, even if no DDD objects
                   are allocated!! (due to free-lists, DDD type definitions
                   etc.) therefore, repeated new/close commands are inhibited
                   explicitly in parallel/dddif/initddd.c(InitCurrMG()). */
  if (MGHEAP(theMG)!=NULL)
    free(MGHEAP(theMG));

  /* dispose BVP */
  if (MG_BVP(theMG)!=NULL)
    if (BVP_Dispose(MG_BVP(theMG))) return (GM_ERROR);

  /* first unlock the mg */
  ((ENVITEM*) theMG)->v.locked = FALSE;

  /* delete mg */
  if (ChangeEnvDir("/Multigrids")==NULL) RETURN (GM_ERROR);
  if (RemoveEnvDir ((ENVITEM *)theMG)) RETURN (GM_ERROR);

  return(GM_OK);
}

/****************************************************************************/
/*
   LexCompare - Define relation for lexicographic ordering

 * @param   pnode1 - first node to compare
 * @param   pnode2 - second node to compare

   This function defines a relation for lexicographic ordering

   @return
   \todo Doc return value!
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

 * @param   LinkHandle1 - first link to compare
 * @param   LinkHandle2 - second link to compare

   This function defines a relation for lexicographic ordering of links

   \todo Doc return value!
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
/** \brief
   OrderNodesInGrid - reorder double linked 'NODE' list

 * @param   theGrid - grid to order
 * @param   order - precedence of coordinate directions
 * @param   sign - respective ordering direction
 * @param   AlsoOrderLinks - if 'TRUE' also order links

   This function reorders the double linked 'NODE' list of the grid with
   qsort and order criteria LexCompare(). If specified the 'LINK's are ordered
   corresponding to the 'NODE' order.

   @return <ul>
   <li>   0 if ok </li>
   <li>   >0 when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX OrderNodesInGrid (GRID *theGrid, const INT *order, const INT *sign, INT AlsoOrderLinks)
{
  MULTIGRID *theMG;
  NODE **table,*theNode;
  LINK *theLink,*LinkTable[LINKTABLESIZE];
  INT i,entries,firstID,nl;
  HEAP *theHeap;
  BVP *theBVP;
  BVP_DESC *theBVPDesc;
  INT MarkKey;

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
  MarkTmpMem(theHeap,&MarkKey);
  if ((table=(NODE**)GetTmpMem(theHeap,entries*sizeof(NODE *),MarkKey))==NULL)
  {
    ReleaseTmpMem(theHeap,MarkKey);
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


  ReleaseTmpMem(theHeap,MarkKey);

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
/** \brief
   PutAtEndOfList - reorder a given set of elements and put them first in the list

 * @param   theGrid - elements are part of that level (not checked)
 * @param   cnt - number of elements in list
 * @param   elemList - list of elements to reorder

   This function reorders a given set of elements and put them last in the list.

   @return <ul>
   <li>   GM_OK if ok </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX PutAtEndOfList (GRID *theGrid, INT cnt, ELEMENT **elemList)
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
/** \brief
 *     Determine neighbor and side of neighbor that goes back to element
 *
 * @param   theElement - considered element
 * @param   Side - side of that element
 * @param   theNeighbor - handle to neighbor
 * @param   NeighborSide - number of side of neighbor that goes back to elem
 *
   This function determines the neighbor and side of the neighbor that goes back to elem.

   @return <ul>
   <li>   0 if ok </li>
   <li>   1 when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX FindNeighborElement (const ELEMENT *theElement, INT Side, ELEMENT **theNeighbor, INT *NeighborSide)
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
/** \brief Insert an inner node
 *
 * @param theGrid grid structure
 * @param pos array containing position
 *
 * This function inserts a inner node into level 0.
 *
 * @return <ul>
 *    <li> pointer to new node if ok </li>
 *    <li> NULL when error occured </li>
 * </ul>
 */
/****************************************************************************/

NODE * NS_DIM_PREFIX InsertInnerNode (GRID *theGrid, DOUBLE *pos)
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
/** \brief
   InsertBoundaryNode - Insert a boundary node

 * @param   theGrid - grid structure
 * @param   bndp - boundary point descriptor

   This function inserts a boundary node into level 0.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

NODE * NS_DIM_PREFIX InsertBoundaryNode (GRID *theGrid, BNDP *bndp)
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
    REP_ERR_RETURN(NULL);
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
    REP_ERR_RETURN(NULL);
  }
        #ifdef TOPNODE
  TOPNODE(theVertex) = theNode;
        #endif

  PRINTDEBUG(dom,1,("  ipn %ld nd %x bndp %x \n",
                    ID(theNode),theNode,V_BNDP(theVertex)));

  SetStringValue(":bndp0",XC(theVertex));
  SetStringValue(":bndp1",YC(theVertex));
        #ifdef __THREEDIM__
  SetStringValue(":bndp2",ZC(theVertex));
        #endif

  return(theNode);
}

/****************************************************************************/
/** \brief
   DeleteNode - Delete a node

 * @param   theGrid - grid structure
 * @param   theNode - node to delete

   This function deletes a node from level 0.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DeleteNode (GRID *theGrid, NODE *theNode)
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
/** \brief
   DeleteNodeWithID - Delete the node with id

 * @param   theGrid - grid structure
 * @param   id - id of node to delete

   This function deletes the node with id `id`.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DeleteNodeWithID (GRID *theGrid, INT id)
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
/** \brief
   FindFather - Find the new father element

 * @param   theVertex -

   This function finds the new father element of the given vertex.
   It assumes that the  new father is one of the neighbors of the
   old father element.

   @return <ul>
   <li>     pointer to an element </li>
   <li>     NULL if none or no correct father is found or vertex is level 0 </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX FindFather (VERTEX *theVertex)
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
/** \brief
   RecreateBNDSofNode - searche the boundary sides and recreate the corresponding BNDS

 * @param   theMG - multigrid structure
 * @param   theNode - node with new BNDP

   This function searches the boundary sides located at 'theNode' and recreate
   the corresponding BNDSs of these sides.
   It assumes that 'theNode' and the neighbour boundary nodes are in the same patch.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
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
/** \brief
   MoveBndMidNode - set new position for a midnode on a boundary

 * @param   theMG - pointer to multigrid
 * @param   theVertex - vertex to move

   This function moves a given boundary vertex according to ist actual local
   coordinates. This function should only be called by MoveMidNode.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured.	 </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MoveBndMidNode (MULTIGRID *theMG, VERTEX *theVertex)
{
  ELEMENT *theElement;
  NODE *Node0,*Node1,*sonNode, *theNode;
  EDGE *theEdge;
  BNDP *bndp;
  BNDS *bnds;
  DOUBLE *global,*local, diffmin, lambda_min;
  DOUBLE_VECTOR bnd_global, bnd_local, BndPoint;
  DOUBLE diff, lambda[DIM_OF_BND], *CornerPtrs[MAX_CORNERS_OF_ELEM];
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
/** \brief
   MoveMidNode - set new position for a midnode

 * @param   theMG - pointer to multigrid
 * @param   theNode - node to move
 * @param   lambda - parameter on the edge
 * @param   update -

   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MoveMidNode (MULTIGRID *theMG, NODE *theNode, DOUBLE lambda, INT update)
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
/** \brief
   MoveCenterNode - set new position for a centernode

 * @param   theMG - pointer to multigrid
 * @param   theNode - node to move
 * @param   lambda - local coordinate in the father element

   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>
 */
/****************************************************************************/

INT NS_DIM_PREFIX MoveCenterNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda)
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
  /*    theElement = FindFather(theVertex); */
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
/** \brief
   MoveSideNode - set new position for a sidenode

 * @param   theMG - pointer to multigrid
 * @param   theNode - node to move
 * @param   lambda - local coordinate on the father element side

   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

#ifdef __THREEDIM__
INT NS_DIM_PREFIX MoveSideNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda)
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
/** \brief
   MoveNode - Let user enter a new position for an inner node

 * @param   theMG - pointer to multigrid
 * @param   theNode - node to move
 * @param   newPos - global coordinate for new position
 * @param   update -

   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MoveNode (MULTIGRID *theMG, NODE *theNode, DOUBLE *newPos, INT update)
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
    if (NTYPE(theNode) == CENTER_NODE)
      theElement = VFATHER(theVertex);
    else
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
/** \brief
 *    Assign new local and global coords to a vertex

 * @param   theMG - pointer to multigrid
 * @param   vert - vertex to move
 * @param   newPos - global coordinate for new position

   This function assigns new global and local coordinates to a vertex (MOVE==DIM).
   It is meant for restoring former consistent positions.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX SetVertexGlobalAndLocal (VERTEX *vert, const DOUBLE *global, const DOUBLE *local)
{
  if (MOVE(vert)!=DIM)
    REP_ERR_RETURN(GM_ERROR);

  if (OBJT(vert)==BVOBJ)
    if (BNDP_Move(V_BNDP(vert),global))
      REP_ERR_RETURN(GM_ERROR);
  V_DIM_COPY(global,CVECT(vert));
  V_DIM_COPY(local,LCVECT(vert));

  return (GM_OK);
}

/****************************************************************************/
/** \brief
   MoveFreeBoundaryVertex - move a vertex on a free boundary

 * @param   theMG - pointer to multigrid
 * @param   vert - vertex to move
 * @param   newPos - global coordinate for new position

   This function moves a given vertex to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.
   To this end the function 'FinishMovingFreeBoundaryVertices' has to be
   called after having done the last call of 'MoveFreeBoundaryVertex'.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul>
 *
   \sa FinishMovingFreeBoundaryVertices
 */
/****************************************************************************/

INT NS_DIM_PREFIX MoveFreeBoundaryVertex (MULTIGRID *theMG, VERTEX *vert, const DOUBLE *newPos)
{
        #ifdef ModelP
  /* TODO: parallel version */
  PrintErrorMessage('E',"MoveFreeBoundaryVertex","parallel not implemented");
  ASSERT(FALSE);
        #endif

  if (OBJT(vert) != BVOBJ)
    REP_ERR_RETURN(GM_ERROR);
  if (MOVE(vert)!=DIM)
    REP_ERR_RETURN(GM_ERROR);

  if (BNDP_Move(V_BNDP(vert),newPos))
    REP_ERR_RETURN(GM_ERROR);
  V_DIM_COPY(newPos,CVECT(vert));

  /* leave local coords for 'FinishMovingFreeBoundaryVertices' */

  return(GM_OK);
}

/****************************************************************************/
/** \brief
   FinishMovingFreeBoundaryVertices - finish moving of free boundary vertices

 * @param   theMG - pointer to multigrid

   The complete multigrid structure is moved hierachically such that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>

 * \sa   MoveFreeBoundaryVertex
 */
/****************************************************************************/

INT NS_DIM_PREFIX FinishMovingFreeBoundaryVertices (MULTIGRID *theMG)
{
  ELEMENT *theElement;
  VERTEX *vert;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  INT n,lev;

  /* adjust global coords of inner vertices */
  for(lev=1; lev<=TOPLEVEL(theMG); lev++)
    for (vert=FIRSTVERTEX(GRID_ON_LEVEL(theMG,lev)); vert!=NULL; vert=SUCCV(vert))
      if ((OBJT(vert) != BVOBJ))
      {
        CORNER_COORDINATES(VFATHER(vert),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(vert),CVECT(vert));
      }
  /* adjust local coords of boundary vertices */
  for(lev=1; lev<=TOPLEVEL(theMG); lev++)
    for (vert=FIRSTVERTEX(GRID_ON_LEVEL(theMG,lev)); vert!=NULL; vert=SUCCV(vert))
      if ((OBJT(vert) == BVOBJ))
      {
        theElement = VFATHER(vert);
        if (theElement == NULL)
          REP_ERR_RETURN(GM_ERROR)
          else
          {
            CORNER_COORDINATES(theElement,n,x);
            UG_GlobalToLocal(n,(const DOUBLE **)x,CVECT(vert),LCVECT(vert));
          }
      }

  RESETMGSTATUS(theMG);

  return(GM_OK);
}

/****************************************************************************/
/** \brief
   GetMidNodeParam - Get local position of a midnode on an edge

 * @param   theNode - midnode
 * @param   lambda  - local coordinate of midnode w.r.t. the edge

   This function gives the local coordinate of a midnode with respect to its edge
   (0 < lambda < 1). The function is called by SmoothGrid and can only be applied
   in 2D.

   @return <ul>
   <li>   GM_OK when ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX GetMidNodeParam (NODE * theNode, DOUBLE *lambda)
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
  {
    PrintErrorMessageF('W',"GetMidNodeParam","could not determine lambda for node %ld",ID(theNode));
  }

  return(GM_OK);
}

#ifdef __TWODIM__


/****************************************************************************/
/** \todo Please doc me!

   CheckOrientation -

   SYNOPSIS:
   INT CheckOrientation (INT n, VERTEX **vertices);


 * @param   n
 * @param   vertices

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

INT NS_DIM_PREFIX CheckOrientation (INT n, VERTEX **vertices)
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


/****************************************************************************/
/** \todo Please doc me!
   CheckOrientation -

   SYNOPSIS:
   INT CheckOrientation (INT n, VERTEX **vertices);


 * @param   n
 * @param   vertices

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

INT NS_DIM_PREFIX CheckOrientation (INT n, VERTEX **vertices)
{
  DOUBLE_VECTOR diff[3],rot;
  DOUBLE det;
  INT i;

  /* TODO: this case */
  if (n == 8 || n==6 || n==5)
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


/****************************************************************************/
/** \todo Please doc me!
   CheckOrientationInGrid -

   SYNOPSIS:
   INT CheckOrientationInGrid (GRID *theGrid);


 * @param   theGrid

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

INT NS_DIM_PREFIX CheckOrientationInGrid (GRID *theGrid)
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


/****************************************************************************/
/** \todo Please doc me!
   NeighborSearch_O_n -

   SYNOPSIS:
   static INT NeighborSearch_O_n(INT n, NODE **Node, MULTIGRID *theMG, INT *NbrS, ELEMENT **Nbr);


 * @param   n
 * @param   Node
 * @param   theMG
 * @param   NbrS
 * @param   Nbr

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

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

    /*CA*/
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

    /*CA*/

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



/****************************************************************************/
/** \todo Please doc me!
   NeighborSearch_O_nn -

   SYNOPSIS:
   static INT NeighborSearch_O_nn(INT n, NODE **Node, GRID *theGrid, INT *NghbrSide, ELEMENT **Nghbr);


 * @param   n
 * @param   Node
 * @param   theGrid
 * @param   NghbrSide
 * @param   Nghbr

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

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


/****************************************************************************/
/** \todo Please doc me!
   NdElPtrArray_evalIndexes -

   SYNOPSIS:
   static INT NdElPtrArray_evalIndexes(INT n, INT *cornerID, MULTIGRID *theMG, INT *MIndex, INT *MBlock, NODE **Node, GRID *theGrid, INT* NbrS, ELEMENT** Nbr);


 * @param   n
 * @param   cornerID
 * @param   theMG
 * @param   MIndex
 * @param   MBlock
 * @param   Node
 * @param   theGrid
 * @param   NbrS
 * @param   Nbr

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

static INT NdElPtrArray_evalIndexes(INT n, INT *cornerID, MULTIGRID *theMG, INT *MIndex, INT *MBlock, NODE **Node, GRID *theGrid, INT* NbrS, ELEMENT** Nbr)
{
  INT retval, j, IndexOfDivPart, IndexOfModPart, Index, merkeIndex, helpIndex;

  /*CA*/
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
  /*CA*/
  return(0);
} /* of static INT NdElPtrArray_evalIndexes() */


/****************************************************************************/
/** \todo Please doc me!
   NdElPtrArray_GetMemAndCheckIDs -

   SYNOPSIS:
   static INT NdElPtrArray_GetMemAndCheckIDs(INT n, MULTIGRID *theMG, INT *h_ID, NODE **c_Node, INT *c_ID, NODE **Node);


 * @param   n
 * @param   theMG
 * @param   h_ID
 * @param   c_Node
 * @param   c_ID
 * @param   NODE

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

static INT NdElPtrArray_GetMemAndCheckIDs(INT n, MULTIGRID *theMG, INT *h_ID, NODE **c_Node, INT *c_ID, NODE **Node)
{
  INT i,j,maxi;
  INT IndexOfDivPart;
  INT IndexOfModPart;
  INT MarkKey = MG_MARK_KEY(theMG);

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

    /* Speicher bereitstellen !!! */
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
                                  theMG		MGNDELEMBLK(theMG,j) = GetTmpMem(theMG->theHeap,maxi,MarkKey);
           */
          /*IE_MEM_PROB*/
          /*
             MGNDELEMBLK(theMG,j) = malloc(maxi);
           */
          if ((MGNDELEMBLK(theMG,j)=(ELEMENT**)GetTmpMem(MGHEAP(theMG),maxi,MarkKey))==NULL)
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



/****************************************************************************/
/** \todo Please doc me!
   Neighbor_Direct_Insert -

   SYNOPSIS:
   static INT Neighbor_Direct_Insert(INT n, ELEMENT **ElemList, INT *NbgSdList, INT* NbrS, ELEMENT **Nbr);


 * @param   n
 * @param   ElemList
 * @param   NbgSdList
 * @param   NbrS
 * @param   Nbr

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

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


/****************************************************************************/
/** \todo Please doc me!
   NdElPtrArray_Update -

   SYNOPSIS:
   static INT NdElPtrArray_Update(INT *MIndex, INT *MBlock, ELEMENT *theElement, MULTIGRID *theMG);


 * @param   MIndex
 * @param   MBlock
 * @param   theElement
 * @param   theMG

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

static INT NdElPtrArray_Update(INT *MIndex, INT *MBlock, ELEMENT *theElement, MULTIGRID *theMG)
{
  INT j;

  /*CA*/
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


/****************************************************************************/
/** \brief
   InsertElement - Insert an element

 * @param   theGrid - grid structure
 * @param   n
 * @param   Node
 * @param   ElemList
 * @param   NbgSdList
 * @param   bnds_flag

   This function inserts an element

   \todo Please doc the return value!

 */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX InsertElement (GRID *theGrid, INT n, NODE **Node, ELEMENT **ElemList, INT *NbgSdList, INT *bnds_flag)
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

    /* We now assume, that:                                         */
    /* if bnds_flag!=NULL && bnds_flag[i]!=0 there has to be a bnds */
    /* so, if not -->error                                          */
    /* or: if bnds_flag==NULL, the domain decides weather there     */
    /* should be a bnds or not (never an error)                     */

    for (j=0; j<m; j++)
      bndp[j] = V_BNDP(sideVertex[j]);

    if (bnds_flag==NULL)
    {
      bnds[i] = BNDP_CreateBndS(MGHEAP(theMG),bndp,m);
      if (bnds[i] != NULL) ElementType = BEOBJ;
    }
    else if (bnds_flag[i]!=0)
    {
      bnds[i] = BNDP_CreateBndS(MGHEAP(theMG),bndp,m);
      assert(bnds[i]!=NULL);
      ElementType = BEOBJ;
    }
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
/** \brief
   InsertElementFromIDs - Insert element with node ids

 * @param   theGrid - grid structure
 * @param   n - number of nodes in node id list
 * @param   idList - ids of the nodes

   This function inserts an element with nodes that have the ids
   given in `idList`,  on level 0.

   @return <ul>
   <li>   pointer to an element if ok </li>
   <li>   NULL when error occured. </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX InsertElementFromIDs (GRID *theGrid, INT n, INT *idList, INT *bnds_flag)
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

  return (InsertElement(GRID_ON_LEVEL(theMG,0),n,Node,NULL,NULL,bnds_flag));
}

/****************************************************************************/
/** \brief
   DeleteElement - Delete an element

 * @param   theMG - multigrid structure
 * @param   theElement - element to delete

   This function deletes an element from level 0.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DeleteElement (MULTIGRID *theMG, ELEMENT *theElement) /* 3D VERSION */
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


/****************************************************************************/
/** \todo Please doc me!
   DeleteElementWithID -

   SYNOPSIS:
   INT DeleteElementWithID (MULTIGRID *theMG, INT id);


 * @param   theMG
 * @param   id

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

INT NS_DIM_PREFIX DeleteElementWithID (MULTIGRID *theMG, INT id)
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
/** \brief
   InsertMesh - Insert a mesh described by the domain

 * @param   theMG - multigrid structure
 * @param   theMesh - mesh structure

   This function inserts all nodes and elements given by the mesh.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR when error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX InsertMesh (MULTIGRID *theMG, MESH *theMesh)
{
  GRID *theGrid;
  ELEMENT *theElement;
  NODE **NList,*Nodes[MAX_CORNERS_OF_ELEM],*ListNode;
  VERTEX **VList;
  INT i,k,n,nv,j,maxlevel,l,move,part;
  INT ElemSideOnBnd[MAX_SIDES_OF_ELEM];
  INT MarkKey = MG_MARK_KEY(theMG);

  if (theMesh == NULL) return(GM_OK);
  if (theMesh->nElements == NULL)
  {
    assert(theMesh->VertexLevel==NULL);
    theGrid = GRID_ON_LEVEL(theMG,0);
    for (i=0; i<theMesh->nBndP; i++)
      if (InsertBoundaryNode(theGrid,theMesh->theBndPs[i]) == NULL)
        REP_ERR_RETURN(GM_ERROR);

    for (i=0; i<theMesh->nInnP; i++)
      if (InsertInnerNode(theGrid,theMesh->Position[i]) == NULL)
        REP_ERR_RETURN(GM_ERROR);
    return(GM_OK);
  }

  /* prepare */
  nv = theMesh->nBndP + theMesh->nInnP;
  VList = (VERTEX **) GetTmpMem(MGHEAP(theMG),nv*sizeof(VERTEX *),MarkKey);
  if (VList == NULL) return(GM_ERROR);
  NList = (NODE **) GetTmpMem(MGHEAP(theMG),nv*sizeof(NODE *),MarkKey);
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
        return(GM_OK);
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
        return(GM_OK);
      SETMOVE(VList[i],move);
      V_BNDP(VList[i]) = theMesh->theBndPs[i];
    }
    for (i=theMesh->nBndP; i<nv; i++)
    {
      VList[i] = CreateInnerVertex(theGrid);
      V_DIM_COPY(theMesh->Position[i-theMesh->nBndP],CVECT(VList[i]));
    }
  }
  if (theMesh->nElements == NULL)
    return(GM_OK);
  for (j=1; j<=theMesh->nSubDomains; j++)
    for (k=0; k<theMesh->nElements[j]; k++)
    {
      if (theMesh->ElementLevel!=NULL) i = theMesh->ElementLevel[j][k];
      else i=0;
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
      if (theMesh->ElemSideOnBnd==NULL)
        theElement = InsertElement (theGrid,n,Nodes,NULL,NULL,NULL);
      else
      {
        for (l=0; l<SIDES_OF_REF(n); l++) ElemSideOnBnd[l] = (theMesh->ElemSideOnBnd[j][k]&(1<<l));
        theElement = InsertElement (theGrid,n,Nodes,NULL,NULL,ElemSideOnBnd);
      }
      SETSUBDOMAIN(theElement,j);
    }

  return(GM_OK);
}

/****************************************************************************/
/** \brief
   FindNodeFromId - Find a node with given id

 * @param   theGrid - grid level to search.

   This function finds a node with given id.

   @return <ul>
   <li>   pointer to that NODE </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

NODE * NS_DIM_PREFIX FindNodeFromId (GRID *theGrid, INT id)
{
  NODE *theNode;

  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    if (ID(theNode)==id) return(theNode);

  return(NULL);
}


/****************************************************************************/
/** \brief
   FindNodeFromPosition - Find node from position

 * @param   theGrid - grid level to search
 * @param   pos - given position
 * @param   tol - tolerance to accept

   This function finds the first node within `tol` from `pos` in 1-norm.

   @return <ul>
   <li>   pointer to NODE  </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

NODE * NS_DIM_PREFIX FindNodeFromPosition (GRID *theGrid, DOUBLE *pos, DOUBLE *tol)
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
/** \brief
   FindVectorFromPosition - Find vector from position

 * @param   theGrid - grid level to search
 * @param   pos - given position
 * @param   tol - tolerance to accept

   This function finds the first vector within `tol` from `pos` in 1-norm.

   @return <ul>
   <li>   pointer to NODE  </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

VECTOR * NS_DIM_PREFIX FindVectorFromPosition (GRID *theGrid, DOUBLE *pos, DOUBLE *tol)
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
/** \brief
   FindVectorFromIndex - Find vector from Index

 * @param   theGrid - grid level to search
 * @param   index - given index

   This function finds the first vector with index.

   @return <ul>
   <li>   pointer to VECTOR  </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

VECTOR * NS_DIM_PREFIX FindVectorFromIndex (GRID *theGrid, INT index)
{
  VECTOR *theVector;

  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
    if (VINDEX(theVector)==index)
      return(theVector);

  return(NULL);
}

/****************************************************************************/
/** \brief
   FindElementFromId - Find element with id

 * @param   theGrid - grid level to search
 * @param   id - id to search

   This function finds an element with the identification `id`. In parallel
   also ghost elements are searched.

   @return <ul>
   <li>   pointer to that ELEMENT </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX FindElementFromId (GRID *theGrid, INT id)
{
  ELEMENT *theElement;

  /* use PFIRSTELEMENT instead of FIRSTELEMENT to search also for ghost elements */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    if (ID(theElement)==id) return(theElement);

  return(NULL);
}

/****************************************************************************/
/** \brief
   PointInElement - Determine whether point is contained in element

 * @param   x - coordinates of given point
 * @param   theElement - element to scan

   This function determines whether a given point specified by coordinates `x`
   is contained in an element.

   @return <ul>
   <li>   0 an error occurred </li>
   <li>   1 point is contained in the element </li>
   <li>   2 point is nearly on one side of the the element </li>
   <li>   3 point is nearly on one edge of the the element </li>
   <li>   4 point is nearly one of the corners of the the element </li>
   <li>   5 point is not contained in the element </li>
   </ul> */
/****************************************************************************/

#ifdef __TWODIM__
INT NS_DIM_PREFIX PointInElement (const DOUBLE *x, const ELEMENT *theElement)  /* 2D version */
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
INT NS_DIM_PREFIX PointInElement (const DOUBLE *global, const ELEMENT *theElement)
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
/** \brief
   PointOnSide - Determine whether point is on an element side

 * @param   x - coordinates of given point
 * @param   theElement - element to scan
 * @param   side - the element side

   This function determines whether a given point specified by coordinates `x`
   is contained in an element side.

   Beware:  The function only tests if the Point is in the plane spawned by the element side.
   The point could be outside the element side area.

   @return <ul>
   <li>   0 not on side </li>
   <li>   1 x is on side </li>
   </ul> */
/****************************************************************************/

#ifdef __TWODIM__
INT NS_DIM_PREFIX PointOnSide(const DOUBLE *global, const ELEMENT *theElement, INT side)
{
  INT n;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE M[DIM+DIM];
  DOUBLE *a, *b;
  DOUBLE det;

  a = &M[0];
  b = &M[DIM];

  CORNER_COORDINATES(theElement,n,x);

  V2_SUBTRACT(x[CORNER_OF_SIDE(theElement,side,1)], x[CORNER_OF_SIDE(theElement,side,0)], a);
  V2_SUBTRACT(global, x[CORNER_OF_SIDE(theElement,side,0)], b);
  det = M2_DET(M);
  if (fabs(det) < SMALL_C)
    return 1;

  return 0;
}
#else
INT NS_DIM_PREFIX PointOnSide(const DOUBLE *global, const ELEMENT *theElement, INT side)
{
  INT n;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE M[DIM*DIM];
  DOUBLE *a, *b, *c;
  DOUBLE det;

  a = &M[0];
  b = &M[DIM];
  c = &M[2*DIM];

  CORNER_COORDINATES(theElement,n,x);

  V3_SUBTRACT(x[CORNER_OF_SIDE(theElement,side,1)], x[CORNER_OF_SIDE(theElement,side,0)], a);
  V3_SUBTRACT(x[CORNER_OF_SIDE(theElement,side,2)], x[CORNER_OF_SIDE(theElement,side,0)], b);
  V3_SUBTRACT(global, x[CORNER_OF_SIDE(theElement,side,0)], c);
  det = M3_DET(M);
  if (fabs(det) < SMALL_C)
    return 1;

  return 0;
}
#endif

/****************************************************************************/
/** \brief
   DOUBLESide - Determine distance of a point to an element side

 * @param   x - coordinates of given point
 * @param   theElement - element to scan
 * @param   side - the element side

   This function determines the distance of a given point specified by coordinates `x`
   from an element side.

   Beware:  The function only tests if the Point is in the plane spawned by the element side.
   The point could be outside the element side area.

   @return <ul>
   <li>   0 not on side </li>
   <li>   1 x is on side </li>
   </ul> */
/****************************************************************************/

#ifdef __TWODIM__
DOUBLE NS_DIM_PREFIX DistanceFromSide(const DOUBLE *global, const ELEMENT *theElement, INT side)
{
  INT n;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE M[DIM+DIM];
  DOUBLE *a, *b;
  DOUBLE det;

  a = &M[0];
  b = &M[DIM];

  CORNER_COORDINATES(theElement,n,x);

  V2_SUBTRACT(x[CORNER_OF_SIDE(theElement,side,1)], x[CORNER_OF_SIDE(theElement,side,0)], a);
  V2_SUBTRACT(global, x[CORNER_OF_SIDE(theElement,side,0)], b);
  det = M2_DET(M);

  return det;
}
#else
DOUBLE NS_DIM_PREFIX DistanceFromSide(const DOUBLE *global, const ELEMENT *theElement, INT side)
{
  INT n;
  DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE M[DIM*DIM];
  DOUBLE *a, *b, *c;
  DOUBLE det;

  a = &M[0];
  b = &M[DIM];
  c = &M[2*DIM];

  CORNER_COORDINATES(theElement,n,x);

  V3_SUBTRACT(x[CORNER_OF_SIDE(theElement,side,1)], x[CORNER_OF_SIDE(theElement,side,0)], a);
  V3_SUBTRACT(x[CORNER_OF_SIDE(theElement,side,2)], x[CORNER_OF_SIDE(theElement,side,0)], b);
  V3_SUBTRACT(global, x[CORNER_OF_SIDE(theElement,side,0)], c);
  det = M3_DET(M);

  return det;
}
#endif

/****************************************************************************/
/** \brief
   FindFlippedElements - Determine whether elements are flipped

 * @param   theMG - multigrid
 * @param   verbose - verbose mode

   This function checks, whether an element is flipped, i.e. posseses a
   negative volume. Works only for tetrahedra.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   > 0 error occured  </li>
   </ul> */
/****************************************************************************/
#ifdef __THREEDIM__
INT NS_DIM_PREFIX FindFlippedElements(MULTIGRID *theMG, INT verbose)
{
  GRID *theGrid;
  ELEMENT *e;
  INT l,n,i,j,fn,found,bfound,bfather;
  DOUBLE *x[MAX_CORNERS_OF_ELEM],*fx[MAX_CORNERS_OF_ELEM],a,b,c,vol;
  DOUBLE_VECTOR diff[3],cp;

  found=bfound=bfather=0;
  /* loop over all elements per grid level */
  for (l=0; l<=TOPLEVEL(theMG); l++)
    for (e=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l)); e!=NULL; e=SUCCE(e))
    {
      if(TAG(e)!=TETRAHEDRON) {
        UserWriteF("Command only for tetras implemented !\n");
        continue;
      }
      CORNER_COORDINATES(e,n,x);
      if(EFATHER(e)!=NULL)
        CORNER_COORDINATES(EFATHER(e),fn,fx);

      for(i=1; i<n; i++)
        V3_SUBTRACT(x[i],x[0],diff[i-1]);
      V3_VECTOR_PRODUCT(diff[0],diff[1],cp);
      V3_SCALAR_PRODUCT(cp,diff[2],vol);
      V_DIM_EUKLIDNORM(diff[0],a);
      V_DIM_EUKLIDNORM(diff[1],b);
      V_DIM_EUKLIDNORM(diff[2],c);
      vol /= (a*b*c);

      /* element is flipped */
      if(vol < FLT_EPSILON) {
        if(verbose) {
          if(EFATHER(e) == NULL)
            UserWriteF("No Father for element defined !\n");
          else
          {
            if(OBJT(EFATHER(e)) == BEOBJ)
              bfather++;
            UserWriteF("Father Element ID %d (SD %d): \n",ID(EFATHER(e)),SUBDOMAIN(EFATHER(e)));
            for(i=0; i<CORNERS_OF_ELEM(EFATHER(e)); i++) {
              UserWriteF("Vertex %d: ",i);
              for(j=0; j<DIM; j++)
                UserWriteF("%f ",fx[i][j]);
              UserWriteF("\n");
            }
          }
          if(OBJT(e) == BEOBJ)
            UserWriteF("Flipped boundary El %d (SD %d): \n",ID(e),SUBDOMAIN(e));
          else
            UserWriteF("Flipped inner El %d (SD %d): \n",ID(e),SUBDOMAIN(e));
          for(i=0; i<CORNERS_OF_ELEM(e); i++) {
            UserWriteF("Vertex %d: ",i);
            for(j=0; j<DIM; j++)
              UserWriteF("%f ",x[i][j]);
            UserWriteF("\n");
          }
        }
        if(OBJT(e) == BEOBJ)
          bfound++;
        else
          found++;
      }
    }
  UserWriteF("-> found %d flipped boundary father elements.\n",bfather);
  UserWriteF("-> found %d flipped boundary sons.\n",found);
  UserWriteF("-> found %d flipped inner sons.\n",bfound);

  return (0);
}
#endif
/****************************************************************************/
/** \brief
   FindElementFromPosition - Find element containing position

 * @param   theGrid - grid level to search
 * @param   pos - given position

   This function finds the first element containing the position `pos`.

   @return <ul>
   <li>   pointer to ELEMENT </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX FindElementFromPosition (GRID *theGrid, DOUBLE *pos)
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
  for (i=0; Sons[i]!=NULL; i++)
    if (PointInElement(pos,Sons[i]) == 1)
      return(Sons[i]);

  return(NULL);
}


/****************************************************************************/
/** \brief
   FindElementFromPosition - Find element containing position

 * @param   theMG - multigrid level to search
 * @param   global - given position

   This function finds the first element containing the position `pos`.

   @return <ul>
   <li>   pointer to ELEMENT </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX FindElementOnSurface (MULTIGRID *theMG, DOUBLE *global)
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
/** \brief
   FindElementOnSurfaceCached - Find element containing position

 * @param   theMG - multigrid level to search
 * @param   global - given position

   This function finds the first element containing the position `pos`.

   @return <ul>
   <li>   pointer to ELEMENT </li>
   <li>   NULL if not found. </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX FindElementOnSurfaceCached (MULTIGRID *theMG, DOUBLE *global)
{
  ELEMENT *t;
  INT k;
  static ELEMENT *e = NULL;

  if ( e!=NULL && EstimateHere(e) )
  {
    /* First try the cached element */
    if (PointInElement(global,e)) {
      return e;
    }

    /* Then try the neighbours */
    for (k=0; k<SIDES_OF_ELEM(e); k++) {
      t = NBELEM(e,k);
      if ( t!=NULL )
        if (PointInElement(global,t))
        {
          e = t;
          return t;
        }
    }
  }

  /* No luck? Do it the hard way. */
  e = FindElementOnSurface(theMG, global);
  return e;
}


/****************************************************************************/
/** \todo Please doc me!
   InnerBoundary -

   SYNOPSIS:
   INT InnerBoundary (ELEMENT *t, INT side);


 * @param   t
 * @param   side

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

INT NS_DIM_PREFIX InnerBoundary (ELEMENT *t, INT side)
{
  INT left,right,part;

  ASSERT(OBJT(t) == BEOBJ);
  ASSERT(SIDE_ON_BND(t,side));

  BNDS_BndSDesc(ELEM_BNDS(t,side),&left,&right,&part);

  return((left != 0) && (right != 0));
}


/****************************************************************************/
/** \brief
   NeighbourElement - get the neighbouring element

 * @param   theElement - pointer to an element
 * @param   side - number of an element side

   This function returns a pointer to the element on the given side.

   @return <ul>
   <li>    pointer to an element </li>
   <li>    NULL if error occured. </li>
   </ul> */
/****************************************************************************/

ELEMENT * NS_DIM_PREFIX NeighbourElement (ELEMENT *t, INT side)
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
/** \brief Calculate the center of mass for an element
 *
 * @param theElement the element
 * @param center_of_mass center of mass as the result
 *
 * This function calculates the center of mass for an arbitrary element.
 * DOUBLE_VECTOR is an array for a 2D resp. 3D coordinate.
 *
 * \sa DOUBLE_VECTOR, ELEMENT
 */
/****************************************************************************/

void NS_DIM_PREFIX CalculateCenterOfMass(ELEMENT *theElement, DOUBLE_VECTOR center_of_mass)
{
  DOUBLE *corner;
  INT i, nr_corners;

  nr_corners = CORNERS_OF_ELEM(theElement);
  V_DIM_CLEAR(center_of_mass);

  for (i=0; i<nr_corners; i++)
  {
    corner = CVECT(MYVERTEX(CORNER(theElement,i)));
    V_DIM_ADD(center_of_mass,corner,center_of_mass);
  }

  V_DIM_SCALE(1.0/nr_corners,center_of_mass);
}

/****************************************************************************/
/** \brief Calculate the center of mass for an element side
 *
 * @param theElement the element
 * @param side index of element side
 * @param center_of_mass center of mass as the result
 *
 * This function calculates the center of mass for an arbitrary element
 * side. DOUBLE_VECTOR is an array for a 2D resp. 3D coordinate.
 * The function calculates the center of mass in global and local coordinats.
 *
 * \sa DOUBLE_VECTOR, ELEMENT
 */
/****************************************************************************/

static void CalculateCenterOfMassOfSide(ELEMENT *theElement, int side, DOUBLE_VECTOR global, DOUBLE_VECTOR local)
{
  DOUBLE *corner;
  DOUBLE *l_corner;
  INT i, nr_corners;

  nr_corners = CORNERS_OF_SIDE(theElement,side);
  V_DIM_CLEAR(global);
  V_DIM_CLEAR(local);

  for (i=0; i<nr_corners; i++)
  {
    corner   = CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,side,i))));
    l_corner = LCVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,side,i))));
    V_DIM_ADD(global,corner,global);
    V_DIM_ADD(local,corner,local);
  }

  V_DIM_SCALE(1.0/nr_corners,global);
  V_DIM_SCALE(1.0/nr_corners,local);
}

/****************************************************************************/
/** \brief
   KeyForObject - calculate an (hopefully) unique key for the geometric object

 * @param   obj - geometric object which from the key is needed (can be one of VERTEX, ELEMENT, NODE or VECTOR)

   This function calculates an (hopefully) unique key for VERTEX,
   ELEMENT, NODE, EDGE and VECTOR typed objects.

   The heuristic is: calculate a 2D/3D position for the geometric object and
   transform this position to a single number by a weighted summation of the
   leading digits of the 2 resp. 3 coordinates and taking from this again
   the sigificant digits and adding the level number.

   APPLICATION:
   Use always an explicit cast to avoid compiler warnings, e.g.
        NODE *theNode;
                KeyForObject((KEY_OBJECT *)theNode);

 * \sa   VERTEX, ELEMENT, NODE, EDGE, VECTOR

   @return
 *   the resulting key
 */
/****************************************************************************/

INT NS_DIM_PREFIX KeyForObject( KEY_OBJECT *obj )
{
  int dummy,i;          /* dummy variable */
  DOUBLE_VECTOR coord;

  if (obj==NULL) return (-1);
  switch( OBJT(obj) )
  {
  /* vertex */
  case BVOBJ :
  case IVOBJ :                  /* both together cover all vertex types */
    return LEVEL(obj)+COORDINATE_TO_KEY(CVECT((VERTEX*)obj),&dummy);

  /* element */
  case BEOBJ :
  case IEOBJ :     for (i=0; i<CORNERS_OF_ELEM((ELEMENT*)obj); i++)
    {
      if(CORNER((ELEMENT*)obj,i)==NULL)
        return (-1);
      if(MYVERTEX(CORNER((ELEMENT*)obj,i))==NULL)
        return (-1);
    }
    /* both together cover all element types */
    CalculateCenterOfMass( (ELEMENT*)obj, coord );
    return LEVEL(obj)+COORDINATE_TO_KEY(coord,&dummy);

  /* node */
  case NDOBJ :     if( MYVERTEX((NODE*)obj) == NULL )
      return (-1);
    return LEVEL(obj)+COORDINATE_TO_KEY(CVECT(MYVERTEX((NODE*)obj)),&dummy);

  /* vector */
  case VEOBJ :     if( VOBJECT((VECTOR*)obj) == NULL )
      return (-1);
    VectorPosition( (VECTOR*)obj, coord );
    return LEVEL(obj)+COORDINATE_TO_KEY(coord,&dummy);

  /* edge */
  case EDOBJ :     if( NBNODE(LINK0((EDGE*)obj)) == NULL )
      return (-1);
    if( MYVERTEX(NBNODE(LINK0((EDGE*)obj))) == NULL )
      return (-1);
    if( NBNODE(LINK1((EDGE*)obj)) == NULL )
      return (-1);
    if( MYVERTEX(NBNODE(LINK1((EDGE*)obj))) == NULL )
      return (-1);
    V_DIM_CLEAR(coord);
    /* sum of the coordinates of the 2 edge corners */
    V_DIM_ADD(coord,CVECT(MYVERTEX(NBNODE(LINK0((EDGE*)obj)))),coord);
    V_DIM_ADD(coord,CVECT(MYVERTEX(NBNODE(LINK1((EDGE*)obj)))),coord);
    /* the midpoint of the line is half of the sum */
    V_DIM_SCALE(0.5,coord);
    /* return the key of the midpoint as the key for the edge */
    return LEVEL(obj)+COORDINATE_TO_KEY(coord,&dummy);

  default :        sprintf( buffer, "unrecognized object type %d", OBJT(obj) );
    PrintErrorMessage('E',"KeyForObject",buffer);
    return(0);
    assert(0);
  }
  return (GM_ERROR);
}

/****************************************************************************/
/** \brief
   ListMultiGrid - List general information about multigrid structure

 * @param   theMG - structure to list
 * @param   isCurrent - is `theMG` current multigrid
 * @param   longformat - print all information or only name of `theMG`

   This function lists general information about a multigrid structure.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListMultiGridHeader (const INT longformat)
{
  if (longformat)
    sprintf(buffer,"   %-20.20s %-20.20s %-20.20s %10.10s %10.10s\n","mg name","domain name","problem name","heap size","heap used");
  else
    sprintf(buffer,"   %-20.20s\n","mg name");
}

void NS_DIM_PREFIX ListMultiGrid (MULTIGRID *theMG, const INT isCurrent, const INT longformat)
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
/** \brief
   MultiGridStatus - List information about refinement type distribution

 * @param   theMG - structure to list
 * @param   gridflag -
 * @param   greenflag
 * @param   lbflag
 * @param   verbose

   This function lists information about multigrids element types.

 * \todo Please return value!
 */
/****************************************************************************/

INT NS_DIM_PREFIX MultiGridStatus (MULTIGRID *theMG, INT gridflag, INT greenflag, INT lbflag, INT verbose)
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
  INT MarkKey;
  INT             *infobuffer;
  INT             **lbinfo;
  INT total_elements,sum_elements;
  INT master_elements,hghost_elements,vghost_elements,vhghost_elements;
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
  MarkTmpMem(MGHEAP(theMG),&MarkKey);
  infobuffer      = (INT *) GetTmpMem(MGHEAP(theMG),(procs+1)*(MAXLEVEL+1)*ELEMENT_PRIOS*sizeof(INT),MarkKey);
  if (infobuffer == NULL) assert(0);

  lbinfo          = (INT **) GetTmpMem(MGHEAP(theMG),(procs+1)*sizeof(INT*),MarkKey);
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
        printf( PFMT "MultiGridStatus: wrong element prio %d\n",me,EPRIO(theElement));
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
    float New;
    float newpergreen;
    float predmax;

    SETMARKCOUNT(REFINEINFO(theMG),markcount[MAXLEVEL]);

    New = markcount[MAXLEVEL]*(2<<(DIM-1))*mg_sum_div_red;
    SETPREDNEW0(REFINEINFO(theMG),New);

    if (mg_greenrules[MAXLEVEL] > 0)
      newpergreen = ((float)mg_greenrulesons[MAXLEVEL][MAX_SONS])/mg_greenrules[MAXLEVEL];
    else
      newpergreen = 0;
    New = markcount[MAXLEVEL]*(2<<(DIM-1))+newpergreen*closuresides[MAXLEVEL];
    SETPREDNEW1(REFINEINFO(theMG),New);

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
    {
      VChannelPtr     *mych;

      mych = (VChannelPtr*)malloc(procs*sizeof(VChannelPtr));

      for (i=1; i<procs; i++)
      {
        mych[i] =ConnSync(i,3917);
        RecvSync(mych[i],(void *)lbinfo[i],(MAXLEVEL+1)*ELEMENT_PRIOS*sizeof(INT));
      }
      Synchronize();
      for (i=1; i<procs; i++)
      {
        DiscSync(mych[i]);
      }
      free(mych);
    }
    else
    {
      VChannelPtr mych;

      mych = ConnSync(master,3917);
      SendSync(mych,(void *)lbinfo[me],(MAXLEVEL+1)*ELEMENT_PRIOS*sizeof(INT));
      Synchronize();
      DiscSync(mych);
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
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
  ReleaseTmpMem(MGHEAP(theMG),MarkKey);
        #endif

  return (GM_OK);
}

/****************************************************************************/
/** \brief
   ListGrids - list general information about grids of multigrid

 * @param   theMG - multigrid structure

   This function lists general information about the grids of a multigrid.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListGrids (const MULTIGRID *theMG)
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
                #if defined(ModelP) && defined(Debug)
    /* output also the object priority counters on each level */
    if (0)
    {
      INT i;
      for (i=1; i<MAX_PRIOS; i++)
        UserWriteF("%c %3d %8d %8ld %8ld %8ld %8ld %8ld %8ld %8ld %8ld %9.3e %9.3e\n",c,l,(int)TOPLEVEL(theMG),
                   (long)NV_PRIO(theGrid,i),(long)NN_PRIO(theGrid,i),(long)NE(theGrid),(long)NT_PRIO(theGrid,i),
                   (long)ns,(long)NVEC_PRIO(theGrid,i),(long)NC(theGrid),(long)NIMAT(theGrid),(float)hmin,(float)hmax);
    }
                #endif
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

            /* any sons ? */
            if (SONNODE(n0)!=NULL && SONNODE(n1)!=NULL)
              if (GetEdge(SONNODE(n0),SONNODE(n1))!=NULL) continue;
            if (MIDNODE(theEdge) != NULL)
            {
              if (SONNODE(n0)!=NULL)
                if (GetEdge(MIDNODE(theEdge),SONNODE(n0))!=NULL) continue;
              if (SONNODE(n1)!=NULL)
                if (GetEdge(MIDNODE(theEdge),SONNODE(n1))!=NULL) continue;
            }
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
    UserWriteF("\n%lu bytes used out of %lu allocated\n",used,used+free);
  else
    UserWriteF("\n%lu ( %lu + %lu ) bytes used out of %lu allocated\n",
               used+heap,used,heap,used+free);
    #endif

    #ifdef ModelP
  used = used + heap;
  used = UG_GlobalMaxINT(used);
  UserWriteF("%lu bytes used on some processor %lu bytes used on all\n",used,UG_GlobalSumINT(used));
    #endif
}

/****************************************************************************/
/** \brief
   ListNode - List information about node in multigrid

 * @param   theMG - structure containing the node
 * @param   theNode - node to list
 * @param   dataopt - list user data if true
 * @param   bopt - list boundary info if true
 * @param   nbopt - list info about neighbors if true
 * @param   vopt - list more information

   This function lists information about a node in a multigrid.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListNode (MULTIGRID *theMG, NODE *theNode, INT dataopt, INT bopt, INT nbopt, INT vopt)
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

    UserWriteF(" key=%d\n", KeyForObject((KEY_OBJECT *)theNode) );

    if (NVECTOR(theNode) != NULL)
      UserWriteF(" vec=" VINDEX_FMTX "\n",
                 VINDEX_PRTX(NVECTOR(theNode)));

    UserWriteF(" classes: NCLASS = %d  NNCLASS = %d\n",NCLASS(theNode),NNCLASS(theNode));
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
      UserWriteF("NB=" ID_FMTX " CTRL=%8lx NO_OF_ELEM=%3d",
                 ID_PRTX(NBNODE(theLink)),(long)CTRL(theLink),NO_OF_ELEM(MYEDGE(theLink)));
      if (MIDNODE(MYEDGE(theLink))!=NULL)
        UserWriteF(" MIDNODE=" ID_FMTX, ID_PRTX(MIDNODE(MYEDGE(theLink))));
      theVertex=MYVERTEX(NBNODE(theLink));
      for(i=0; i<DIM; i++)
      {
        UserWriteF(" x%1d=%11.4E",i, (float)(CVECT(theVertex)[i]) );
      }
      UserWrite("\n");
    }
  }
  return;
}


/****************************************************************************/
/** \brief
   ListNodeSelection - List information about all nodes in selection

 * @param   theMG - structure containing the nodes
 * @param   dataopt - list user data if true
 * @param   bopt - list boundary info if true
 * @param   nbopt - list info about neighbors if true
 * @param   vopt - list more information

   This function lists information about all nodes in the selection.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListNodeSelection (MULTIGRID *theMG, INT dataopt, INT bopt, INT nbopt, INT vopt)
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
/** \brief Check if element is in selection list
 *
 * @param   theMG - multigrid structure
 * @param   theNode - node to check
 *
 * This function checks if an element is in the selection list.
 *
 * @return <ul>
 * <li>   0 if NOT in list
 * <li>   1 if in list.
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX IsNodeSelected (MULTIGRID *theMG, NODE *theNode)
{
  int j;

  if (SELECTIONMODE(theMG) != nodeSelection) return (0);
  for(j=0; j<SELECTIONSIZE(theMG); j++)
    if (theNode == (NODE *) SELECTIONOBJECT(theMG,j))
      return (1);
  return (0);
}


/****************************************************************************/
/** \brief
   ListNodeRange - List information about nodes in given range of ids

 * @param   theMG - structure to list
 * @param   from - first id
 * @param   to - last id
 * @param   idopt - determines the meaning of from/to
 * @param   dataopt - list user data if true
 * @param   bopt - list boundary info if true
 * @param   nbopt - list info about neighbors if true
 * @param   vopt - list more information

   This function list information about all nodes in a given range of ids.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListNodeRange (MULTIGRID *theMG, INT from, INT to, INT idopt, INT dataopt, INT bopt, INT nbopt, INT vopt)
{
  int level;
  NODE *theNode;

  for (level=0; level<=TOPLEVEL(theMG); level++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,level)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      switch( idopt )
      {
      case 0 :                          /* $i option */
        if ( (ID(theNode)>=from)&&(ID(theNode)<=to) )
          ListNode(theMG,theNode,dataopt,bopt,nbopt,vopt);
        break;
#ifdef ModelP
      case 1 :                          /* $g option */
        if (GID(theNode) == from)
          ListNode(theMG,theNode,dataopt,bopt,nbopt,vopt);
        break;
#endif
      case 2 :                          /* $k option */
        if ( KeyForObject((KEY_OBJECT *)theNode) == from)
          ListNode(theMG,theNode,dataopt,bopt,nbopt,vopt);
        break;

      default : PrintErrorMessage( 'E', "ListNodeRange", "unrecognized idopt" );
        assert(0);
      }
    }
}


/****************************************************************************/
/** \brief
   ListElement - List information about element

 * @param   theMG -  structure to list
 * @param   theElement - element to list
 * @param   dataopt - list user data if true
 * @param   vopt - list more information

   This function lists information about an element

 */
/****************************************************************************/

void NS_DIM_PREFIX ListElement (MULTIGRID *theMG, ELEMENT *theElement, INT dataopt, INT bopt, INT nbopt, INT vopt)
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
    UserWriteF("subdomain=%d \n", SUBDOMAIN(theElement));
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
    /** \todo delete this
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
  }
  if (nbopt)
  {
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if (NBELEM(theElement,i)!=NULL)
      {
        UserWriteF("    NB%d=" EID_FMTX ,
                   i,EID_PRTX(NBELEM(theElement,i)));
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
/** \brief
    ListElementSelection - list information about elements in selection

 * @param  theMG multigrid structure to list
 * @param   dataopt - list user data if true
 * @param   vopt - list more information

   This function lists information about all elements in the selection.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListElementSelection (MULTIGRID *theMG, INT dataopt, INT bopt, INT nbopt, INT vopt)
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
/** \brief
   IsElementSelected - Check whether element is in selection list

 * @param   theMG - multigrid structure
 * @param   theElement - element to check

   This function checks whether an element is in the selection list.

   @return <ul>
   <li>   0 if NOT in list </li>
   <li>   1 if in list.  </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX IsElementSelected (MULTIGRID *theMG, ELEMENT *theElement)
{
  int j;

  if (SELECTIONMODE(theMG) != elementSelection) return (0);
  for(j=0; j<SELECTIONSIZE(theMG); j++)
    if (theElement == (ELEMENT *) SELECTIONOBJECT(theMG,j))
      return (1);
  return (0);
}


/****************************************************************************/
/** \brief
   ListElementRange - List information about elements in range of ids

 * @param   theMG - multigrid structure to list
 * @param   from - first id
 * @param   to - last id
 * @param   idopt - determines the meaning of from/to
 * @param   dataopt - list user data if true
 * @param   vopt - list more information

   This function lists information about all elements in a range of ids.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListElementRange (MULTIGRID *theMG, INT from, INT to, INT idopt, INT dataopt, INT bopt, INT nbopt, INT vopt, INT lopt)
{
  int level,fromlevel,tolevel;
  ELEMENT *theElement;

  if (lopt==FALSE)
  {
    fromlevel = 0;
    tolevel = TOPLEVEL(theMG);
  }
  else
    fromlevel = tolevel = CURRENTLEVEL(theMG);

  for (level=fromlevel; level<=tolevel; level++)
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,level)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      switch( idopt )
      {
      case 0 :                          /* $i option */
        if ( (ID(theElement)>=from)&&(ID(theElement)<=to) )
          ListElement(theMG,theElement,dataopt,bopt,nbopt,vopt);
        break;
#ifdef ModelP
      case 1 :                          /* $g option */
        if (EGID(theElement) == from)
          ListElement(theMG,theElement,dataopt,bopt,nbopt,vopt);
        break;
#endif
      case 2 :                          /* $k option */
        if ( KeyForObject((KEY_OBJECT *)theElement) == from)
          ListElement(theMG,theElement,dataopt,bopt,nbopt,vopt);
        break;

      default : PrintErrorMessage( 'E', "ListElementRange", "unrecognized idopt" );
        assert(0);
      }
    }
}


/****************************************************************************/
/** \brief
   ListVector - List information about vector

 * @param   theMG - multigrid structure to list
 * @param   theVector - vector to list
 * @param   matrixopt - list line of matrix corresponding to theVector
 * @param   dataopt - list user data if true
 * @param   modifiers - flags modifying output style and verbose level

   This function lists information about a vector.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListVector (MULTIGRID *theMG, VECTOR *theVector, INT matrixopt, INT dataopt, INT modifiers)
{
  FORMAT *theFormat;
  NODE *theNode;
  EDGE *theEdge;
  ELEMENT *theElement;
  MATRIX *theMatrix;
  DOUBLE_VECTOR pos;
  void *Data;

  theFormat = MGFORMAT(theMG);

  /* print index and type of vector */
  UserWriteF("IND=" VINDEX_FFMTE " VTYPE=%d(%c) ",
             VINDEX_PRTE(theVector),
             VTYPE(theVector),
             FMT_T2N(theFormat,VTYPE(theVector)));

  if (READ_FLAG(modifiers,LV_POS))
  {
    if (VectorPosition(theVector,pos))
      return;
                #ifdef __TWODIM__
    UserWriteF("POS=(%10.2e,%10.2e)",pos[_X_],pos[_Y_]);
                #endif
                #ifdef __THREEDIM__
    UserWriteF("POS=(%10.2e,%10.2e,%10.2e)",pos[_X_],pos[_Y_],pos[_Z_]);
                #endif
  }

  /* print object type of vector */
  if (READ_FLAG(modifiers,LV_VO_INFO))
    switch (VOTYPE(theVector))
    {
    case NODEVEC :
    {
      theNode = (NODE*)VOBJECT(theVector);
                                #if defined __OVERLAP2__ || defined USE_FAMG
      if ( theNode == NULL )
        UserWriteF("NODE-V NULL                ");
      else
                                #endif
      UserWriteF("NODE-V nodeID=" ID_FMTX
                 "                ",
                 ID_PRTX(theNode));
    }
    break;
    case EDGEVEC :
    {
      theEdge = (EDGE*)VOBJECT(theVector);
      UserWriteF("EDGE-V fromID=" ID_FFMT
                 " to__ID=%7ld ",
                 ID_PRT(NBNODE(LINK0(theEdge))),
                 ID(NBNODE(LINK1(theEdge))));
    }
    break;
                #ifdef __THREEDIM__
    case SIDEVEC :
    {
      theElement = (ELEMENT*)VOBJECT(theVector);
      UserWriteF("SIDE-V elemID=" EID_FFMT
                 "                ",
                 EID_PRT(theElement));
    }
    break;
                #endif
    case ELEMVEC :
    {
      theElement = (ELEMENT*)VOBJECT(theVector);
      UserWriteF("ELEM-V elemID=" EID_FFMT
                 "                ",
                 EID_PRT(theElement));
    }
    break;

    default : PrintErrorMessage( 'E', "ListVector", "unrecognized VECTOR type" );
      assert(0);
    }

  UserWriteF("VCLASS=%1d VNCLASS=%1d",VCLASS(theVector),VNCLASS(theVector));
  UserWriteF(" key=%d\n", KeyForObject((KEY_OBJECT *)theVector) );

  /* print vector data if */
  if (dataopt && FMT_PR_VEC(theFormat)!=NULL)
  {
    /* print skip flags */
    if (READ_FLAG(modifiers,LV_SKIP))
    {
      INT_2_bitpattern(VECSKIP(theVector),buffer);
      UserWriteF("  skip=%s\n",buffer);
    }

    /* print data */
    Data = (void*)(&VVALUE(theVector,0));
    if ((*(FMT_PR_VEC(theFormat)))(VTYPE(theVector),Data,"   ",buffer))
      return;
    UserWrite(buffer);
  }

  /* print matrix list if */
  if (matrixopt > 0)
    for (theMatrix = VSTART(theVector); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      UserWrite("    DEST(MATRIX): ");
      ListVector(theMG,MDEST(theMatrix),0,0,modifiers);

      /* print matrix data if */
      if (dataopt && theFormat->PrintMatrix!=NULL)
      {
        Data = (void*)(&MVALUE(theMatrix,0));
        if ((*(FMT_PR_MAT(theFormat)))(MTYPE(theMatrix), Data, "       ", buffer))
          return;
        UserWrite(buffer);
      }
    }
  if (matrixopt < 0)
    for (theMatrix = VISTART(theVector); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      UserWrite("    DEST(MATRIX): ");
      ListVector(theMG,MDEST(theMatrix),0,0,modifiers);

      /* print matrix data if */
      if (dataopt)
      {
        UserWriteF("  P = %8.6lf, ", MVALUE(theMatrix,0));
        UserWriteF("  R = %8.6lf \n", MVALUE(theMatrix,1));
      }
    }
  return;
}

/****************************************************************************/
/** \brief
   ListVectorOfElementSelection - List info about vectors of elements in selection

 * @param   theMG -  structure to list
 * @param   matrixopt - list line of matrix corresponding to theVector
 * @param   dataopt - list user data if true
 * @param   modifiers - flags modifying output style and verbose level

   This function lists info about all vectors of elements in the selection.

 */
/****************************************************************************/

void NS_DIM_PREFIX  ListVectorOfElementSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt, INT modifiers)
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
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt,modifiers);
    }
    if (VEC_DEF_IN_OBJ_OF_MG(theMG,EDGEVEC))
    {
      GetVectorsOfEdges(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt,modifiers);
    }
                #ifdef __THREEDIM__
    if (VEC_DEF_IN_OBJ_OF_MG(theMG,SIDEVEC))
    {
      GetVectorsOfSides(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt,modifiers);
    }
                #endif
    if (VEC_DEF_IN_OBJ_OF_MG(theMG,ELEMVEC))
    {
      GetVectorsOfElement(theElement,&cnt,vList);
      for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt,modifiers);
    }
  }
}

/****************************************************************************/
/** \brief
   ListVectorSelection - list information about vectors in selection

 * @param   theMG: multigrid structure to list
 * @param   matrixopt - list matrices of this vector
 * @param   dataopt - list user data if true
 * @param   modifiers - flags modifying output style and verbose level

   This function lists information about all elements in the selection.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListVectorSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt, INT modifiers)
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
    ListVector(theMG,theVector,matrixopt,dataopt,modifiers);
  }
}

/****************************************************************************/
/** \brief
   IsVectorSelected - Check whether vector is in selection list

 * @param   theMG - multigrid structure
 * @param   theVector - vector to check

   This function checks whether an element is in the selection list.

   @return <ul>
   <li>   0 if NOT in list </li>
   <li>   1 if in list.  </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX IsVectorSelected (MULTIGRID *theMG, VECTOR *theVector)
{
  int j;

  if (SELECTIONMODE(theMG) != vectorSelection) return (0);
  for(j=0; j<SELECTIONSIZE(theMG); j++)
    if (theVector == (VECTOR *) SELECTIONOBJECT(theMG,j))
      return (1);
  return (0);
}

/****************************************************************************/
/** \brief
   ListVectorRange - list information about vectors in range of ids

 * @param   theMG - structure to list
 * @param   from - first index
 * @param   to - last index
 * @param   idopt - determines the meaning of from/to
 * @param   matrixopt - list line of matrix corresponding to theVector
 * @param   dataopt - list user data if true
 * @param   datatypes - list vectors with type in datatypes
 * @param   modifiers - flags modifying output style and verbose level

   This function lists information about all vectors in a given range of indices.

 */
/****************************************************************************/

void NS_DIM_PREFIX ListVectorRange (MULTIGRID *theMG, INT fl, INT tl, INT from, INT to, INT idopt, INT matrixopt, INT dataopt, INT datatypes, INT modifiers)
{
  int level;
  VECTOR *theVector;

  for (level=fl; level<=tl; level++)
    for (theVector=PFIRSTVECTOR(GRID_ON_LEVEL(theMG,level)); theVector!=NULL; theVector=SUCCVC(theVector))
      if (datatypes & VDATATYPE(theVector))
      {
        switch( idopt )
        {
        case LV_ID :                                            /* $i option */
          if (VINDEX(theVector)>=from && VINDEX(theVector)<=to)
            ListVector(theMG,theVector,matrixopt,dataopt,modifiers);
          break;
#ifdef ModelP
        case LV_GID :                                   /* $g option */
          if (GID(theVector) == from)
            ListVector(theMG,theVector,matrixopt,dataopt,modifiers);
          break;
#endif
        case LV_KEY :                                   /* $k option */
          if ( KeyForObject((KEY_OBJECT *)theVector) == from)
            ListVector(theMG,theVector,matrixopt,dataopt,modifiers);
          break;

        default : PrintErrorMessage( 'E', "ListVectorRange", "unrecognized idopt" );
          assert(0);
        }
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
/** \brief
   ClearSelection - Clear selection buffer

 * @param   theMG - multigrid structure

   This function clears the selection buffer of a multigrid.

 */
/****************************************************************************/


void NS_DIM_PREFIX ClearSelection (MULTIGRID *theMG)
{
  SELECTIONSIZE(theMG) = 0;
}


/****************************************************************************/
/** \brief Add node to selection buffer
 *
 * @param theMG multigrid structure
 * @param theNode node to add
 *
 * This function adds an node to the selection buffer or removes it if it is already
 * in the list.
 *
 * @return <ul>
 * <li>   GM_OK if ok </li>
 * <li>   GM_ERROR if an  error occured. </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX AddNodeToSelection (MULTIGRID *theMG, NODE *theNode)
{
  int i;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=nodeSelection) return(GM_ERROR);
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

  if (SELECTIONSIZE(theMG)>=MAXSELECTION) return(GM_ERROR);

  SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG)) = g;
  SELECTIONSIZE(theMG)++;
  return(GM_OK);
}


/****************************************************************************/
/** \brief Add element to selection buffer
 *
 * @param theMG multigrid structure
 * @param theElement element to add
 *
 * This function adds an element to the selection buffer or removes it if it is already
 * in the list.
 *
 * @return <ul>
 *    <li> GM_OK if ok </li>
 *    <li> GM_ERROR if an  error occured </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX AddElementToSelection (MULTIGRID *theMG, ELEMENT *theElement)
{
  int i;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=elementSelection) return(GM_ERROR);
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
    return(GM_ERROR);

  SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG)) = g;
  SELECTIONSIZE(theMG)++;
  return(GM_OK);
}

/**************************************************************/
/** \brief Add vector to selection buffer
 *
 * @param theMG multigrid structure
 * @param theVector vector to add
 *
 * This function adds a vector to the selection buffer or removes it if it is already
 * in the list.
 *
 * @return <ul>
 *   <li>   GM_OK if ok </li>
 *   <li>   GM_ERROR if an  error occured. </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX AddVectorToSelection (MULTIGRID *theMG, VECTOR *theVector)
{
  int i;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=vectorSelection) return(GM_ERROR);
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

  if (SELECTIONSIZE(theMG)>=MAXSELECTION) return(GM_ERROR);

  SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG)) = g;
  SELECTIONSIZE(theMG)++;
  return(GM_OK);
}

/****************************************************************************/
/*
   RemoveNodeFromSelection - Remove node from selection buffer

 * @param   theMG - multigrid structure
 * @param   theNode - node to remove

   This function removes an node from the selection buffer.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if an error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX RemoveNodeFromSelection (MULTIGRID *theMG, NODE *theNode)
{
  int i,j,found;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=nodeSelection) return(GM_ERROR);
  }
  else return(GM_ERROR);

  g = (SELECTION_OBJECT *) theNode;
  found = 0;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      found = 1;
      break;
    }

  if (!found) return(GM_ERROR);

  for (j=i+1; j<SELECTIONSIZE(theMG); j++)
    SELECTIONOBJECT(theMG,j-1) = SELECTIONOBJECT(theMG,j);

  SELECTIONSIZE(theMG)--;
  return(GM_OK);
}


/****************************************************************************/
/*
   RemoveElementFromSelection - Remove element from selection buffer

 * @param   theMG - multigrid structure
 * @param   theElement - element to remove

   This function removes an element from the selection buffer.

   @return <ul>
   INT
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if an error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX RemoveElementFromSelection (MULTIGRID *theMG, ELEMENT *theElement)
{
  int i,j,found;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=elementSelection) return(GM_ERROR);
  }
  else return(GM_ERROR);

  g = (SELECTION_OBJECT *) theElement;
  found = 0;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      found = 1;
      break;
    }

  if (!found) return(GM_ERROR);

  for (j=i+1; j<SELECTIONSIZE(theMG); j++)
    SELECTIONOBJECT(theMG,j-1) = SELECTIONOBJECT(theMG,j);

  SELECTIONSIZE(theMG)--;
  return(GM_OK);
}

/****************************************************************************/
/*
   RemoveVectorFromSelection - Remove vector from selection buffer

 * @param   theMG - multigrid structure
 * @param   theVector - vector to remove

   This function removes a vector from the selection buffer.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if an error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX RemoveVectorFromSelection (MULTIGRID *theMG, VECTOR *theVector)
{
  int i,j,found;
  SELECTION_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=vectorSelection) return(GM_ERROR);
  }
  else return(GM_ERROR);

  g = (SELECTION_OBJECT *) theVector;
  found = 0;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g)
    {
      found = 1;
      break;
    }

  if (!found) return(GM_ERROR);

  for (j=i+1; j<SELECTIONSIZE(theMG); j++)
    SELECTIONOBJECT(theMG,j-1) = SELECTIONOBJECT(theMG,j);

  SELECTIONSIZE(theMG)--;
  return(GM_OK);
}
/****************************************************************************/
/** \todo Please doc me!
   GetAngle -

   SYNOPSIS:
   static INT GetAngle(DOUBLE *angle,DOUBLE *n1, DOUBLE *n2);


 * @param   angle
 * @param   n1
 * @param   n2

   DESCRIPTION:

   @return
   INT
 */
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


/****************************************************************************/
/** \todo Please doc me!
   SetNormal-

   SYNOPSIS:
   static INT SetNormal(DOUBLE *n, DOUBLE **x, INT nc);


 * @param   n
 * @param   x
 * @param   nc

   DESCRIPTION:

   @return
   INT
 */
/****************************************************************************/

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



/****************************************************************************/
/** \brief
   MinMaxAngle - Determine min and max angle in degrees

 * @param   theElement - element to check
 * @param   amin - minimum angle
 * @param   amax - maximum angle

   This function determines min and max angle in degrees.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if an error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MinMaxAngle (ELEMENT *theElement, DOUBLE *amin, DOUBLE *amax)
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

static INT MinMaxEdge (ELEMENT *theElement, DOUBLE *amin, DOUBLE *amax)
{
  INT error,i,s1,s2,tag;
  DOUBLE angle,*x[MAX_CORNERS_OF_SIDE],l,n1[DIM],n2[DIM],min,max;

  error=GM_OK;
  tag=TAG(theElement);

  min = DBL_MAX;
  max = DBL_MIN;
  for (s1=0; s1<EDGES_OF_TAG(tag); s1++)
  {

    /* get corner coordinates of side and evaluate normal */
    for (i=0; i<CORNERS_OF_EDGE; i++)
      x[i]=CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE_TAG(tag,s1,i))));

    V_DIM_EUKLIDNORM_OF_DIFF(x[0],x[1],l);

    max = MAX(max,l);
    min = MIN(min,l);

  }
  l = max/min;
  /* UserWriteF("id=%d length=%lf max=%lf min=%lf\n",ID(theElement),l,max,min); */
  *amax = MAX(*amax,l);
  *amin = MIN(*amin,l);
  return(error);
}

/****************************************************************************/
/** \brief
   DefineMGUDBlock - Define block in general MG user data space

 * @param   id - the id of the block to be allocated
 * @param   size - size of the data block

   This function defines a block in the general MG user data space.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if an  error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX DefineMGUDBlock (BLOCK_ID id, MEM size)
{
  if (DefineBlock(theGenMGUDM,id,size)!=0)
    return (GM_ERROR);

  return (GM_OK);
}


/****************************************************************************/
/** \brief
   FreeMGUDBlock - Free block in general MG user data space

 * @param   id: the id of the block to be allocated

   This function frees a block in the general MG user data space.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if an error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX FreeMGUDBlock (BLOCK_ID id)
{
  if (FreeBlock(theGenMGUDM,id)!=0)
    return (GM_ERROR);

  return (GM_OK);
}

/****************************************************************************/
/** \brief
   GetMGUDBlockDescriptor - Return pointer to block descriptor with id

 * @param   id - the id of the block to be allocated

   This function returns a pointer to the block descriptor with id.

   @return <ul>
   <li>   pointer to BLOCK_DESC </li>
   <li>   NULL if an error occured. </li>
   </ul> */
/****************************************************************************/

BLOCK_DESC      * NS_DIM_PREFIX GetMGUDBlockDescriptor (BLOCK_ID id)
{
  return (GetBlockDesc(theGenMGUDM,id));
}

VIRT_HEAP_MGMT * NS_DIM_PREFIX GetGenMGUDM()
{
  return (theGenMGUDM);
}

/****************************************************************************/
/** \brief
   MaxNodeClass - Returns highest Node class of a dof on next level

 * @param   theElement - pointer to a element

   This function returns highest 'NCLASS' of a Node associated with the
   element.

   @return <ul>
   <li>    0 if ok </li>
   <li>    1 if error occured.		 </li>
   </ul> */
/****************************************************************************/

static INT MaxNodeClass (ELEMENT *theElement)
{
  INT m = 0;
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++) {
    INT c = NCLASS(CORNER(theElement,i));

    m = MAX(m,c);
  }

  return (m);
}

/****************************************************************************/
/** \brief
   MaxNextNodeClass - Returns highest Node class of a dof on next level

 * @param   theElement - pointer to a element

   This function returns highest 'NNCLASS' of a Node associated with the
   element.

   @return <ul>
   <li>    0 if ok </li>
   <li>    1 if error occured.  </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MaxNextNodeClass (ELEMENT *theElement)
{
  INT m = 0;
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++) {
    INT c = NNCLASS(CORNER(theElement,i));

    m = MAX(m,c);
  }

  return (m);
}

#ifdef __PERIODIC_BOUNDARY__
INT NS_DIM_PREFIX MaxPeriodicNextNodeClass (ELEMENT *theElement)
{
  INT m = 0;
  INT i;

  if (PeriodicBoundaryInfo==NULL)
    return (0);

  /* check if it is a periodic element */
  if (OBJT(theElement)==BEOBJ) {
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++) {
      INT c = NNCLASS(CORNER(theElement,i));
      INT nbnd,per_ids[MAX_PERIODIC_OBJ];
      DOUBLE_VECTOR coord,percoord[MAX_PERIODIC_OBJ];
      VERTEX *vtx = MYVERTEX(CORNER(theElement,i));

      if (OBJT(vtx)==BVOBJ)
        if ((PeriodicBoundaryInfo)(vtx,&nbnd,per_ids,coord,percoord))
          m = MAX(m,c);
    }
    return (m);
  }
  else
    return (0);
}
#endif

/****************************************************************************/
/** \brief
   MinNodeClass - Returns minimal Node class of a dof on next level

 * @param   theElement - pointer to a element

   This function returns highest 'NNCLASS' of a Node associated with the
   element.

   @return <ul>
   <li>    0 if ok </li>
   <li>    1 if error occured.  </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MinNodeClass (ELEMENT *theElement)
{
  INT m = 3;
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++) {
    INT c = NCLASS(CORNER(theElement,i));

    m = MIN(m,c);
  }

  return (m);
}

/****************************************************************************/
/** \brief
   MinNextNodeClass - Returns minimal Node class of a dof on next level

 * @param   theElement - pointer to a element

   This function returns highest 'NNCLASS' of a Node associated with the
   element.

   @return <ul>
   <li>    0 if ok </li>
   <li>    1 if error occured.  </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MinNextNodeClass (ELEMENT *theElement)
{
  INT m = 3;
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++) {
    INT c = NNCLASS(CORNER(theElement,i));

    m = MIN(m,c);
  }

  return (m);
}

/****************************************************************************/
/** \brief
   SeedNodeClasses - Initialize node classes

 * @param   theGrid - given grid
 * @param   theElement - given element

   Initialize Node class in all nodes associated with given element with 3.

   @return <ul>
   <li>    0 if ok </li>
   <li>    1 if error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX SeedNodeClasses (ELEMENT *theElement)
{
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    SETNCLASS(CORNER(theElement,i),3);

  return (0);
}

/****************************************************************************/
/** \brief
   ClearNodeClasses - Reset node classes

 * @param   theGrid - pointer to grid

   Reset all node classes in all nodes of given grid to 0.

   @return <ul>
   <li>     0 if ok </li>
   <li>     1 if error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX ClearNodeClasses (GRID *theGrid)
{
  NODE *theNode;

  /* reset class of each Node to 0 */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    SETNCLASS(theNode,0);

  return(0);
}
/****************************************************************************/
/** \brief
   PropagateNodeClasses - Compute Node classes after initialization

 * @param   theGrid - pointer to grid

   After Node classes have been reset and initialized, this function
   now computes the class 2 and class 1 Nodes.

   @return <ul>
   <li>      0 if ok </li>
   <li>      1 if error occured </li>
   </ul> */
/****************************************************************************/

#ifdef ModelP
static int Gather_NodeClass (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  ((INT *)data)[0] = NCLASS(theNode);

  return(0);
}

static int Scatter_NodeClass (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  SETNCLASS(theNode,MAX(NCLASS(theNode),((INT *)data)[0]));

  return(0);
}

static int Scatter_GhostNodeClass (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  SETNCLASS(theNode,((INT *)data)[0]);

  return(0);
}
#endif

static INT PropagateNodeClass (GRID *theGrid, INT nclass)
{
  ELEMENT *theElement;
  INT i;

  for (theElement=FIRSTELEMENT(theGrid);
       theElement!= NULL; theElement = SUCCE(theElement))
    if (MaxNodeClass(theElement) == nclass)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        NODE *theNode = CORNER(theElement,i);

        if (NCLASS(theNode) < nclass)
          SETNCLASS(theNode,nclass-1);
      }

        #ifdef __PERIODIC_BOUNDARY__
  {
    /* due to periodicity node class for periodic nodes has to be made consistent */
    VECTOR *vec;
    NODE *node;

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETVCLASS(vec,MAX(NCLASS(node),VCLASS(vec)));
    }

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETNCLASS(node,MAX(VCLASS(vec),NCLASS(node)));
    }

  }
        #endif

  /* only for this values valid */
  ASSERT(nclass==3 || nclass==2);

  return(0);
}

#ifdef __PERIODIC_BOUNDARY__
static INT PropagatePeriodicNodeClass (GRID *theGrid, INT nclass)
{
  ELEMENT *theElement;
  INT i;

  for (theElement=FIRSTELEMENT(theGrid);
       theElement!= NULL; theElement = SUCCE(theElement))
    if (MaxNodeClass(theElement) == nclass)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        NODE *theNode = CORNER(theElement,i);

        if (NCLASS(theNode) < nclass)
          SETNCLASS(theNode,nclass);
      }

  /* only for this values valid */
  ASSERT(nclass==1);

  return(0);
}
#endif

#ifdef _SCHALE_X_
static INT PropagateNodeClassX (GRID *theGrid, INT nclass)
{
  ELEMENT *theElement;
  INT i;

  /* first: set node class to nclass-1 */
  for (theElement=FIRSTELEMENT(theGrid);
       theElement!= NULL; theElement = SUCCE(theElement))
    if (MaxNodeClass(theElement) == nclass)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        NODE *theNode = CORNER(theElement,i);

        if (NCLASS(theNode) < nclass)
          SETNCLASS(theNode,nclass-1);
      }

  {
    NODE *node;

    /* second: set node class to nclass where class==nclass-1
           otherwise dependence on order not connectivity! */
    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node))
      if (NCLASS(node)==(nclass-1))
        SETNCLASS(node,nclass);
  }

        #ifdef __PERIODIC_BOUNDARY__
  {
    /* due to periodicity node class for periodic nodes has to be made consistent */
    VECTOR *vec;
    NODE *node;

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETVCLASS(vec,MAX(NCLASS(node),VCLASS(vec)));
    }

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETNCLASS(node,MAX(VCLASS(vec),NCLASS(node)));
    }

  }
        #endif

  return(0);
}
#endif

INT NS_DIM_PREFIX PropagateNodeClasses (GRID *theGrid)
{
  NODE *theNode;
  MATRIX *theMatrix;

#ifdef _SCHALE_X_
    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNodeClasses():"
                   " 0. communication on level %d\n",me,GLEVEL(theGrid)))
  /* exchange NCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NodeClass, Scatter_NodeClass);
    #endif

  /* set Node classes in the algebraic neighborhood to 3 */
  if (PropagateNodeClassX(theGrid,3)) REP_ERR_RETURN(1);
#endif

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNodeClasses():"
                   " 1. communication on level %d\n",me,GLEVEL(theGrid)))
  /* exchange NCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NodeClass, Scatter_NodeClass);
    #endif

  /* set Node classes in the algebraic neighborhood to 2 */
  if (PropagateNodeClass(theGrid,3)) REP_ERR_RETURN(1);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNodeClasses(): 2. communication\n",
                   me))
  /* exchange NCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NodeClass, Scatter_NodeClass);
    #endif

  /* set Node classes in the algebraic neighborhood to 1 */
  if (PropagateNodeClass(theGrid,2)) REP_ERR_RETURN(1);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNodeClasses(): 3. communication\n",
                   me))
  /* exchange NCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NodeClass, Scatter_NodeClass);
        #endif

#ifdef __PERIODIC_BOUNDARY__
  if (1) {
    /* set Node classes in periodic neighbourhood to 1 */
    if (PropagatePeriodicNodeClass(theGrid,1)) REP_ERR_RETURN(1);

    #ifdef ModelP
    PRINTDEBUG(gm,1,("\n" PFMT "PropagateNodeClasses(): 4. communication\n",
                     me))
    /* exchange NCLASS of Nodes */
    DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                    Gather_NodeClass, Scatter_NodeClass);
        #endif
  }
#endif

        #ifdef ModelP
  /* send NCLASS to ghosts */
  DDD_IFAOneway(NodeIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(INT),
                Gather_NodeClass, Scatter_GhostNodeClass);
    #endif

  return(0);
}

/****************************************************************************/
/** \brief Reset class of the Nodes on the next level

 * @param   theGrid - pointer to grid

   This function clears NNCLASS flag in all Nodes. This is the first step to
   compute the class of the dofs on the *NEXT* level, which
   is also the basis for determining copies.

   @return <ul>
   <li>     0 if ok </li>
   <li>     1 if error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX ClearNextNodeClasses (GRID *theGrid)
{
  NODE *theNode;

  /* reset class of each Node to 0 */
  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    SETNNCLASS(theNode,0);

  /* now the refinement algorithm will initialize the class 3 Nodes   */
  /* on the *NEXT* level.                                                                               */
  return(0);
}

/****************************************************************************/
/** \brief
   SeedNextNodeClasses - Set 'NNCLASS' in all Nodes associated with element

 * @param   theElement - pointer to element

   Set 'NNCLASS' in all nodes associated with the element to 3.

   @return <ul>
   <li>     0 if ok  </li>
   <li>     1 if error occured. </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX SeedNextNodeClasses (ELEMENT *theElement)
{
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    SETNNCLASS(CORNER(theElement,i),3);

  return (0);
}

/****************************************************************************/
/** \brief
   PropagateNextNodeClasses - Compute 'NNCLASS' in all Nodes of a grid level

 * @param   theGrid - pointer to grid

   Computes values of 'NNCLASS' field in all nodes after seed.

   @return <ul>
   <li>    0 if ok  </li>
   <li>    1 if error occured </li>
   </ul> */
/****************************************************************************/

#ifdef ModelP
static int Gather_NextNodeClass (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  ((INT *)data)[0] = NNCLASS(theNode);

  return(GM_OK);
}

static int Scatter_NextNodeClass (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  SETNNCLASS(theNode,MAX(NNCLASS(theNode),((INT *)data)[0]));

  return(GM_OK);
}

static int Scatter_GhostNextNodeClass (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  SETNNCLASS(theNode,((INT *)data)[0]);

  return(GM_OK);
}
#endif

static INT PropagateNextNodeClass (GRID *theGrid, INT nnclass)
{
  ELEMENT *theElement;
  INT i;

  for (theElement=FIRSTELEMENT(theGrid);
       theElement!= NULL; theElement = SUCCE(theElement))
    if (MaxNextNodeClass(theElement) == nnclass)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        NODE *theNode = CORNER(theElement,i);

        if (NNCLASS(theNode) < nnclass)
          SETNNCLASS(theNode,nnclass-1);
      }

        #ifdef __PERIODIC_BOUNDARY__
  {
    /* due to periodicity periodic next node class for periodic nodes has to be made consistent */
    VECTOR *vec;
    NODE *node;

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETVNCLASS(vec,MAX(NNCLASS(node),VNCLASS(vec)));
    }

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETNNCLASS(node,MAX(VNCLASS(vec),NNCLASS(node)));
    }

  }
        #endif

  /* only for this values valid */
  ASSERT(nnclass==3 || nnclass==2);

  return(0);
}

#ifdef __PERIODIC_BOUNDARY__
static INT PropagatePeriodicNextNodeClass (GRID *theGrid, INT nnclass)
{
  /* set NNCLASS=1 in periodic neighbourhood to set class consistently */
  ELEMENT *theElement;
  INT i;

  for (theElement=FIRSTELEMENT(theGrid);
       theElement!= NULL; theElement = SUCCE(theElement))
    if (MaxPeriodicNextNodeClass(theElement) == nnclass)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        NODE *theNode = CORNER(theElement,i);

        if (NNCLASS(theNode) < nnclass)
          SETNNCLASS(theNode,nnclass);
      }

  ASSERT(nnclass==1);

  return (0);
}
#endif

#ifdef _SCHALE_X_
static INT PropagateNextNodeClassX (GRID *theGrid, INT nnclass)
{
  ELEMENT *theElement;
  INT i;

  /* first: set next node class to nnclass-1 */
  for (theElement=FIRSTELEMENT(theGrid);
       theElement!= NULL; theElement = SUCCE(theElement))
    if (MaxNextNodeClass(theElement) == nnclass)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        NODE *theNode = CORNER(theElement,i);

        if (NNCLASS(theNode) < nnclass)
          SETNNCLASS(theNode,nnclass-1);
      }

  {
    NODE *node;

    /* second: set next node class to nnclass where class==nnclass-1
       otherwise dependence on order not connectivity! */
    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      if (NNCLASS(node)==nnclass-1)
        SETNNCLASS(node,nnclass);
    }
  }

#ifdef __PERIODIC_BOUNDARY__
  {
    /* due to periodicity periodic next node class for periodic nodes has to be made consistent */
    VECTOR *vec;
    NODE *node;

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETVNCLASS(vec,MAX(NNCLASS(node),VNCLASS(vec)));
    }

    for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node)) {
      vec=NVECTOR(node);

      SETNNCLASS(node,MAX(VNCLASS(vec),NNCLASS(node)));
    }

  }
#endif

  return (0);
}
#endif

INT NS_DIM_PREFIX PropagateNextNodeClasses (GRID *theGrid)
{
  NODE *theNode;
  MATRIX *theMatrix;

#ifdef _SCHALE_X_
    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextNodeClasses(): 0. communication\n",me))
  /* exchange NNCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NextNodeClass, Scatter_NextNodeClass);
    #endif

  if (PropagateNextNodeClassX(theGrid,3)) REP_ERR_RETURN(1);
#endif

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextNodeClasses(): 1. communication\n",me))
  /* exchange NNCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NextNodeClass, Scatter_NextNodeClass);
    #endif

  if (PropagateNextNodeClass(theGrid,3)) REP_ERR_RETURN(1);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextNodeClasses(): 2. communication\n",me))
  /* exchange NNCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NextNodeClass, Scatter_NextNodeClass);
    #endif

  if (PropagateNextNodeClass(theGrid,2)) REP_ERR_RETURN(1);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextNodeClasses(): 3. communication\n",me))
  /* exchange NNCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NextNodeClass, Scatter_NextNodeClass);
        #endif

#ifdef __PERIODIC_BOUNDARY__
  /* set next node class to 1 in periodic neighbourhood */
  if (PropagatePeriodicNextNodeClass(theGrid,1)) REP_ERR_RETURN (1);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextNodeClasses(): 4. communication\n",me))
  /* exchange NNCLASS of Nodes */
  DDD_IFAExchange(BorderNodeSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_NextNodeClass, Scatter_NextNodeClass);
        #endif
#endif

        #ifdef ModelP
  /* send NNCLASSn to ghosts */
  DDD_IFAOneway(NodeIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(INT),
                Gather_NextNodeClass, Scatter_GhostNextNodeClass);
    #endif

  return(0);
}

/****************************************************************************/
/** \brief
   SetEdgeAndNodeSubdomainFromElements - set subdomain id on level 0 edges

 * @param   id - the id of the block to be allocated

   This function sets the subdomain id taken from the elements for level 0 edges.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
   </ul> */
/****************************************************************************/

#ifndef __cplusplus
static
#endif
INT NS_DIM_PREFIX SetEdgeAndNodeSubdomainFromElements (GRID *theGrid)
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

  /* now change subdomain id for boundary edges and nodes to 0 */
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
          SETNSUBDOM(n0,0);
          ASSERT(OBJT(MYVERTEX(n0)) == BVOBJ);
          SETNSUBDOM(n1,0);
          ASSERT(OBJT(MYVERTEX(n1)) == BVOBJ);
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
      PRINTDEBUG(gm,1,("  ed(%d,%d)-sd=%d nsub %d %d\n",ID(n0),ID(n1),
                       EDSUBDOM(ed),NSUBDOM(n0),NSUBDOM(n1)));
    }
  }
  ENDDEBUG

  return (GM_OK);
}

/****************************************************************************/
/** \brief
   RemoveSpuriousBoundarySides - remove boundary side of element and neighbour

 * @param   elem - element
 * @param   side - boundary side to remove

   This function removes the boundary side of element and neighbour.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
   </ul> */
/****************************************************************************/

static INT RemoveSpuriousBoundarySides (HEAP *heap, ELEMENT *elem, INT side)
{
  ELEMENT *nb=NBELEM(elem,side);
  BNDS *nbbside,*bside=ELEM_BNDS(elem,side);
  INT nbside;

  ASSERT(bside!=NULL);
  ASSERT(OBJT(elem)==BEOBJ);
  ASSERT(nb!=NULL);
  ASSERT(OBJT(nb)==BEOBJ);

  /* search nbside */
  for (nbside=0; nbside<SIDES_OF_ELEM(nb); nbside++)
    if (NBELEM(nb,nbside)==elem)
      break;
  ASSERT(nbside<SIDES_OF_ELEM(nb));
  nbbside = ELEM_BNDS(nb,nbside);
  ASSERT(nbbside!=NULL);

  PRINTDEBUG(gm,1,("spurious bsides between elem %ld and elem %ld removed",(long)ID(elem),(long)ID(nb)));

  HEAPFAULT(bside);
  if (BNDS_Dispose(heap,bside))
    REP_ERR_RETURN(1);
  SET_BNDS(elem,side,NULL);

  HEAPFAULT(nbbside);
  if (BNDS_Dispose(heap,nbbside))
    REP_ERR_RETURN(2);
  SET_BNDS(nb,nbside,NULL);

  return (0);
}

/****************************************************************************/
/** \brief Replace boundary by inner element
 *
 * @param grid grid where elem resides
 * @param elemH handle to element to be replaced: will point to new element
 *
 * This function replaces a boundary by an inner element.
 *
 * @return <ul>
 *   <li>   GM_OK if ok </li>
 *   <li>   GM_ERROR if error occured </li>
 **</ul>
 */
/****************************************************************************/

static INT BElem2IElem (GRID *grid, ELEMENT **elemH)
{
  ELEMENT *nb[MAX_SIDES_OF_ELEM],*elem=*elemH,*ielem;
  NODE *nodes[MAX_CORNERS_OF_ELEM];
  INT i,j,nbside[MAX_SIDES_OF_ELEM],s_id;

  ASSERT(GLEVEL(grid)==0);

  /* save context */
  for (i=0; i<CORNERS_OF_ELEM(elem); i++)
    nodes[i] = CORNER(elem,i);

  for (i=0; i<SIDES_OF_ELEM(elem); i++)
  {
    nb[i] = NBELEM(elem,i);
    for (j=0; j<SIDES_OF_ELEM(nb[i]); j++)
      if (NBELEM(nb[i],j)==elem)
        break;
    ASSERT(j<SIDES_OF_ELEM(nb[i]));
    nbside[i] = j;
  }

  s_id = SUBDOMAIN(elem);

  /* create/dispose */
  ielem = CreateElement(grid,TAG(elem),IEOBJ,nodes,EFATHER(elem),NO);
  if (ielem==NULL)
    REP_ERR_RETURN(1);

  if (DisposeElement(grid,elem,NO))
    REP_ERR_RETURN(1);

  *elemH = ielem;

  /* set context */
  for (i=0; i<SIDES_OF_ELEM(ielem); i++)
  {
    SET_NBELEM(ielem,i,nb[i]);

    SET_NBELEM(nb[i],nbside[i],ielem);
  }
  SETSUBDOMAIN(ielem,s_id);
  SETECLASS(ielem,RED_CLASS);

  return (0);
}

/****************************************************************************/
/** \brief
   FinishGrid - remove erroneously introduced bsides and propagate sub domain IDs

 * @param   mg - multigrid

   This function removes erroneously introduced bsides and propagates sub domain IDs.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
   </ul> */
/****************************************************************************/

static INT FinishGrid (MULTIGRID *mg)
{
  GRID *grid;
  ELEMENT *elem,*nb,*succ;
  HEAP *heap=MGHEAP(mg);
  FIFO unused,shell;
  INT MarkKey = MG_MARK_KEY(mg);
  INT i,side,id,nbid,part,nsd,found,s_id;
  INT *sd_table;
  void *buffer;

  /* prepare */
  if (TOPLEVEL(mg)<0)
    REP_ERR_RETURN (GM_ERROR);
  grid = GRID_ON_LEVEL(mg,0);
  if (!NT(grid))
    return (GM_OK);

  for (elem=PFIRSTELEMENT(grid); elem!=NULL; elem=SUCCE(elem))
  {
    SETUSED(elem,FALSE);
    SETTHEFLAG(elem,FALSE);
  }

  /* table for subdomain ids */
  nsd = 1 + BVPD_NSUBDOM(MG_BVPD(mg));
  sd_table = (INT*)GetTmpMem(heap,nsd*sizeof(INT),MarkKey);
  if (sd_table==NULL)
    REP_ERR_RETURN (GM_ERROR);

  /* init two fifos */
  buffer=(void *)GetTmpMem(heap,sizeof(ELEMENT*)*NT(grid),MarkKey);
  if (buffer==NULL)
    REP_ERR_RETURN (GM_ERROR);
  fifo_init(&unused,buffer,sizeof(ELEMENT*)*NT(grid));
  buffer=(void *)GetTmpMem(heap,sizeof(ELEMENT*)*NT(grid),MarkKey);
  if (buffer==NULL)
    REP_ERR_RETURN (GM_ERROR);
  fifo_init(&shell,buffer,sizeof(ELEMENT*)*NT(grid));

  /* outermost loop handles nonconnected domains */
  while (TRUE)
  {
    for (elem=PFIRSTELEMENT(grid); elem!=NULL; elem=SUCCE(elem))
      if (!USED(elem))
        break;
    if (elem!=NULL)
      fifo_in(&unused,elem);
    else
      break;

    while (!fifo_empty(&unused))
    {
      /* grab next !USED element */
      do
        elem = (ELEMENT*) fifo_out(&unused);
      while (USED(elem) && !fifo_empty(&unused));
      if (USED(elem))
        /* we are done */
        break;

      /* shell algo (using FLAG): neighbours, but not across bside */
      fifo_clear(&shell);
      fifo_in(&shell,elem);
      SETTHEFLAG(elem,TRUE);
      for (i=0; i<=nsd; i++) sd_table[i] = 0;
      found = FALSE;
      while (!fifo_empty(&shell))
      {
        elem = (ELEMENT*) fifo_out(&shell);

        if (OBJT(elem)==BEOBJ)
          for (side=0; side<SIDES_OF_ELEM(elem); side++)
            if (SIDE_ON_BND(elem,side))
            {
              if (BNDS_BndSDesc(ELEM_BNDS(elem,side),&id,&nbid,&part))
                REP_ERR_RETURN (GM_ERROR);

              if ((nb=NBELEM(elem,side))==NULL)
              {
                /* this bside must be ok (outer boundary) */
                /* TODO (HRR 971012): parallel? */
                ASSERT(nbid==0);
                s_id = id;
                found = TRUE;
                break;
              }
              else
              if (USED(nb))
              {
                /* he must know! */
                if (nbid==SUBDOMAIN(nb))
                  s_id = id;
                else if (id==SUBDOMAIN(nb))
                  s_id = nbid;
                else
                  ASSERT(FALSE);
              }

              /* handle outer boundary cases */
              if (id==0)
              {
                ASSERT(nbid>0);
                s_id = nbid;
                found = TRUE;
                break;
              }
              if (nbid==0)
              {
                ASSERT(id>0);
                s_id = id;
                found = TRUE;
                break;
              }

              ++sd_table[id];
              if (sd_table[id]>1)
              {
                s_id = id;
                found = TRUE;
                break;
              }
            }
        if (found)
          break;

        /* push neighbours not across boundary */
        if (OBJT(elem)==BEOBJ)
        {
          for (side=0; side<SIDES_OF_ELEM(elem); side++)
            if (!SIDE_ON_BND(elem,side))
              if ((nb=NBELEM(elem,side))!=NULL)
                if (!USED(nb) && !THEFLAG(nb))
                {
                  fifo_in(&shell,nb);
                  SETTHEFLAG(nb,TRUE);
                }
        }
        else
        {
          for (side=0; side<SIDES_OF_ELEM(elem); side++)
            if ((nb=NBELEM(elem,side))!=NULL)
              if (!USED(nb) && !THEFLAG(nb))
              {
                fifo_in(&shell,nb);
                SETTHEFLAG(nb,TRUE);
              }
        }
      }

      /* count occurences of subdom ids (max 2 different) */
      for (found=0, i=0; i<=nsd; i++)
        if (sd_table[i])
          found++;
      if (found>2)
        /* FATAL: algorithm relies on assumptions obviously not fulfilled! */
        ASSERT(FALSE);

      /* again shell algo starting from last element */
      /* set USED, propagate subdomain ids and remove spurious bsides (has to have NB!) */
      /* use PRINTDEBUG */
      fifo_clear(&shell);
      fifo_in(&shell,elem);
      SETUSED(elem,TRUE);
      SETSUBDOMAIN(elem,s_id);
      while (!fifo_empty(&shell))
      {
        elem = (ELEMENT*) fifo_out(&shell);

        if (OBJT(elem)==BEOBJ)
        {
          for (side=0; side<SIDES_OF_ELEM(elem); side++)
            if (SIDE_ON_BND(elem,side))
            {
              if ((nb=NBELEM(elem,side))==NULL)
                continue;
              if (!USED(nb))
                /* push unused neighbour across boundary to unused fifo */
                fifo_in(&unused,nb);

              if (BNDS_BndSDesc(ELEM_BNDS(elem,side),&id,&nbid,&part))
                REP_ERR_RETURN (GM_ERROR);

              if (id!=s_id || nbid==0)
              {
                /* remove spurious bside of both elements */
                if (RemoveSpuriousBoundarySides(heap,elem,side))
                  REP_ERR_RETURN(1);
              }
            }
        }

        /* push neighbours not across boundary */
        if (OBJT(elem)==BEOBJ)
        {
          for (side=0; side<SIDES_OF_ELEM(elem); side++)
            if (!SIDE_ON_BND(elem,side))
              if ((nb=NBELEM(elem,side))!=NULL)
              {
                if (!USED(nb))
                {
                  fifo_in(&shell,nb);
                  SETUSED(nb,TRUE);
                  SETSUBDOMAIN(nb,s_id);
                }
              }
              else
                /* TODO (HRR 971012): ModelP: no error if EGHOST? */
                /* grid not closed */
                REP_ERR_RETURN(1);
        }
        else
        {
          for (side=0; side<SIDES_OF_ELEM(elem); side++)
            if ((nb=NBELEM(elem,side))!=NULL)
            {
              if (!USED(nb))
              {
                fifo_in(&shell,nb);
                SETUSED(nb,TRUE);
                SETSUBDOMAIN(nb,s_id);
              }
            }
            else
              /* TODO (HRR 971012): ModelP: no error if EGHOST? */
              /* grid not closed */
              REP_ERR_RETURN(1);
        }
      }
    }
  }

  for (elem=PFIRSTELEMENT(grid); elem!=NULL; elem=succ)
  {
    succ = SUCCE(elem);

    if (OBJT(elem)!=BEOBJ) continue;

    /* check whether element still has bsides */
    for (side=0; side<SIDES_OF_ELEM(elem); side++)
      if (ELEM_BNDS(elem,side)!=NULL)
        break;
    if (side>=SIDES_OF_ELEM(elem))
      if (BElem2IElem(grid,&elem))
        REP_ERR_RETURN(1);
  }

  if (SetEdgeAndNodeSubdomainFromElements(grid))
    REP_ERR_RETURN (GM_ERROR);

  return (GM_OK);
}

/****************************************************************************/
/** \brief
   SetSubdomainIDfromBndInfo - set subdomain id on level 0 elements and edges

 * @param   id - the id of the block to be allocated

   This function sets the subdomain for level 0 elements and edges.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX SetSubdomainIDfromBndInfo (MULTIGRID *theMG)
{
  HEAP *theHeap;
  GRID *theGrid;
  ELEMENT *theElement, *theNeighbor;
  NODE *theNode;
  void *buffer;
  INT i,n,id,nbid,part,j;
  FIFO myfifo;
  INT MarkKey = MG_MARK_KEY(theMG);

  /* prepare */
  if (TOPLEVEL(theMG)<0) REP_ERR_RETURN (GM_ERROR);
  theGrid = GRID_ON_LEVEL(theMG,0);
  n = NT(theGrid);        if (n==0) return(0);

  /* allocate fifo and init */
  theHeap = MYMG(theGrid)->theHeap;
  buffer=(void *)GetTmpMem(theHeap,sizeof(ELEMENT*)*n,MarkKey);
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

  IFDEBUG(gm,1)
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    assert(USED(theElement));
  ENDDEBUG

  if (SetEdgeAndNodeSubdomainFromElements(theGrid))
    REP_ERR_RETURN (GM_ERROR);

  return (GM_OK);
}



#ifdef __PERIODIC_BOUNDARY__

static int sort_entries_old (const void *e1, const void *e2)
{
  PERIODIC_ENTRIES *v1 = (PERIODIC_ENTRIES *)e1;
  PERIODIC_ENTRIES *v2 = (PERIODIC_ENTRIES *)e2;

  /* periodic ids */
  if (v1->periodic_id < v2->periodic_id) return (-1);
  if (v1->periodic_id > v2->periodic_id) return (1);

  if (v1->coord[0] < v2->coord[0] - SMALL_DOUBLE) return (-1);
  if (v1->coord[0] > v2->coord[0] + SMALL_DOUBLE) return (1);

  if (v1->coord[1] < v2->coord[1] - SMALL_DOUBLE) return (-1);
  if (v1->coord[1] > v2->coord[1] + SMALL_DOUBLE) return (1);
        #ifdef __THREEDIM__
  /* y coordinates almost equal compare z now */
  if (v1->coord[2] < v2->coord[2] - SMALL_DOUBLE) return (-1);
  if (v1->coord[2] > v2->coord[2] + SMALL_DOUBLE) return (1);
        #endif

  return (0);
}

static int sort_entries (const void *e1, const void *e2)
{
  INT i;
  PERIODIC_ENTRIES *v1 = (PERIODIC_ENTRIES *)e1;
  PERIODIC_ENTRIES *v2 = (PERIODIC_ENTRIES *)e2;

  if (v1->n < v2->n) return (-1);
  if (v1->n > v2->n) return (1);

  for (i=0; i<v1->n; i++)
  {
    if (PVID(v1->vp[i]) < PVID(v2->vp[i])) return (-1);
    if (PVID(v1->vp[i]) > PVID(v2->vp[i])) return (1);
  }

  /* periodic ids */
  if (v1->periodic_id < v2->periodic_id) return (-1);
  if (v1->periodic_id > v2->periodic_id) return (1);

  for (i=0; i<DIM; i++)
  {
    if (v1->coord[i] < v2->coord[i] - SMALL_DOUBLE) return (-1);
    if (v1->coord[i] > v2->coord[i] + SMALL_DOUBLE) return (1);
  }

  return(0);
}

/*
   static INT PositionsMatch(DOUBLE_VECTOR pos1, DOUBLE_VECTOR pos2)
   {
   if (ABS(pos1[0]-pos2[0])>SMALL_DOUBLE)
        return (1);
   if (ABS(pos1[1]-pos2[1])>SMALL_DOUBLE)
        return (1);
   #ifdef __THREEDIM__
   if (ABS(pos1[2]-pos2[2])>SMALL_DOUBLE)
        return (1);
   #endif

   return (0);
   }
 */

INT NS_DIM_PREFIX GridSetPerVecCount (GRID *g)
{
  NODE *n;

  if (PeriodicBoundaryInfo == NULL)
  {
    UserWriteF("Grid_GeometricToPeriodic: no function *PeriodicBoundaryInfo\n");
    return(GM_OK);
  }

  for (n=PFIRSTNODE(g); n!=NULL; n=SUCCN(n))
    SETPVCOUNT(NVECTOR(n),0);

  for (n=PFIRSTNODE(g); n!=NULL; n=SUCCN(n))
    SETPVCOUNT(NVECTOR(n),(PVCOUNT(NVECTOR(n))+1));

  IFDEBUG(gm,1)
  {
    INT max_nodes_per_vec=0;
    VECTOR *vec;

    for (n=PFIRSTNODE(g); n!=NULL; n=SUCCN(n)) {
      max_nodes_per_vec=MAX(max_nodes_per_vec,PVCOUNT(NVECTOR(n)));
    }
    PrintDebug(" at most %d nodes per vector on level %d !!\n",max_nodes_per_vec,GLEVEL(g));

    if (max_nodes_per_vec>1) {
      for (vec=PFIRSTVECTOR(g); vec!=NULL; vec=SUCCVC(vec))
        if (PVCOUNT(vec)==max_nodes_per_vec)
          PrintDebug(" vec(%08x)= " VINDEX_FMTX " with PVCOUNT=%d\n",vec,VINDEX_PRTX(vec),PVCOUNT(vec));
    }
  }
  ENDDEBUG

  return(0);
}

INT NS_DIM_PREFIX MGSetPerVecCount (MULTIGRID *mg)
{
  INT level;

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    GRID *g = GRID_ON_LEVEL(mg,level);

    if (GridSetPerVecCount(g)) return(GM_ERROR);
  }

  return (GM_OK);
}

static INT ModifyVectorPointer (GRID *g, PERIODIC_ENTRIES *list, INT i, INT j)
{
  MATRIX *m,*n;
  VECTOR *d,*w,*v;

  v = NVECTOR(list[i].node);
  w = NVECTOR(list[j].node);

  /* only node vectors implemented */
  if (VOTYPE(v) != NODEVEC) assert(0);
  if (VOTYPE(w) != NODEVEC) assert(0);

  /* same node or vectors already identified */
  if (v == w) return(0);

  m = START(w);
  n = START(v);

  if (m != NULL) {
    PRINTDEBUG(gm,1,("ModifyVectorPointer list[%d]:\n",j));

    for (m=NEXT(m); m != NULL; m=NEXT(m))
    {
      INT periodic_ids[MAX_PERIODIC_OBJ];
      DOUBLE_VECTOR own_coord, periodic_coords[MAX_PERIODIC_OBJ];
      INT nbnd,IsOnSamePerBnd;
      VERTEX *vtx;

      IsOnSamePerBnd=FALSE;
      d = MDEST(m);

      vtx = MYVERTEX((NODE*)VOBJECT(d));

      /* is d on periodic boundary */
      if (OBJT(vtx) == BVOBJ)
      {
        nbnd = 0;
        if ((*PeriodicBoundaryInfo)(vtx,&nbnd,periodic_ids,own_coord,periodic_coords))
        {
          INT k;

          for (k=0; k<nbnd; k++)
            if (periodic_ids[k]==list[i].periodic_id)
            {
              IsOnSamePerBnd=TRUE;
              break;
            }
        }

        if (IsOnSamePerBnd)
        {
          PRINTDEBUG(gm,1,("%d perid=%d: do not create connection to pos %f %f %f\n",i,list[i].periodic_id,XC(vtx),YC(vtx),ZC(vtx)));
          continue;
        }
      }
      PRINTDEBUG(gm,1,("%d perid=%d: create connection to pos %f %f %f\n",i,list[i].periodic_id,XC(vtx),YC(vtx),ZC(vtx)));

      n = GetMatrix(v,d);
      if (n == NULL)
      {
        /* create connection only to vectors not belonging to bnd with same periodic_id */
        n = CreateConnection(g,v,d);
        if (n == NULL) assert(0);
      }
    }
  }

  return(0);
}


static INT DisposeAndModVector(GRID *g, PERIODIC_ENTRIES *list, INT i, INT j)
{
  VECTOR *vec;
  MATRIX *m;

  /* dispose this vector*/
  vec = NVECTOR(list[j].node);

  if (vec == NVECTOR(list[i].node))
    /* nothing to do */
    return (0);

  /* modify vector pointer */
  NVECTOR(list[j].node) = NVECTOR(list[i].node);

  /* dispose vector and increase node count */
  SETPVCOUNT(NVECTOR(list[i].node),(PVCOUNT(NVECTOR(list[i].node))+1));

  PRINTDEBUG(gm,1,("increased node count for vec(%08x) = "VINDEX_FMTX ": %d\n",NVECTOR(list[i].node),VINDEX_PRTX(NVECTOR(list[i].node)),PVCOUNT(NVECTOR(list[i].node))));

  PRINTDEBUG(gm,1,(PFMT "DisposeAndModVector perid=%d vtx=%d node="
                   ID_FMTX " vec(%08x)=" VINDEX_FMTX " partner vec(%08x)=" VINDEX_FMTX " pos %lf %lf %lf level %d\n",
                   me,list[j].periodic_id,ID(MYVERTEX(list[j].node)),
                   ID_PRTX(list[j].node),vec,VINDEX_PRTX(vec),(NVECTOR(list[i].node)),VINDEX_PRTX(NVECTOR(list[i].node)),
                   XC(MYVERTEX(list[j].node)),
                   YC(MYVERTEX(list[j].node)),
                   ZC(MYVERTEX(list[j].node)),GLEVEL(g)));

  if (DisposeVector(g,vec))
    return(1);

  return (0);
}

INT NS_DIM_PREFIX SetPeriodicBoundaryInfoProcPtr (PeriodicBoundaryInfoProcPtr PBI)
{
  PeriodicBoundaryInfo = PBI;

  return(0);
}

INT NS_DIM_PREFIX GetPeriodicBoundaryInfoProcPtr (PeriodicBoundaryInfoProcPtr *PBI)
{
  *PBI = PeriodicBoundaryInfo;

  return(0);
}


static INT InsertVIDs (INT nn, PERIODIC_ENTRIES *coordlist, NODE *node)
{
  switch (NTYPE(node))
  {
  case (CORNER_NODE) :
  {
    NODE *fnode = NFATHER(node);
    coordlist[nn].vp[0] = NVECTOR(fnode);
    coordlist[nn].n = 1;
    break;
  }
  case (MID_NODE) :
  {
    EDGE *fedge = (EDGE *)NFATHER(node);
    NODE *n0 = NBNODE(LINK0(fedge));
    NODE *n1 = NBNODE(LINK1(fedge));
    assert(n0!=NULL && n1!=NULL);

    if (PVID(NVECTOR(n0)) < PVID(NVECTOR(n1)))
    {
      coordlist[nn].vp[0] = NVECTOR(n0);
      coordlist[nn].vp[1] = NVECTOR(n1);
    }
    else
    {
      coordlist[nn].vp[0] = NVECTOR(n1);
      coordlist[nn].vp[1] = NVECTOR(n0);
    }
    coordlist[nn].n = 2;
    break;
  }
#ifdef __THREEDIM__
  case (SIDE_NODE) :
  {
    INT side,i,k;
    VERTEX  *v              = MYVERTEX(node);
    ELEMENT *felem  = VFATHER(v);
    VECTOR *vp[MAX_CORNERS_OF_SIDE];
    assert(felem != NULL);
                        #ifdef ModelP
    assert(EMASTER(felem));
                        #endif

    for (i=0; i<MAX_CORNERS_OF_SIDE; i++) vp[i] = NULL;

    side = ONSIDE(v);
    assert(side<MAX_SIDES_OF_ELEM);

    /* extract vps */
    for (i=0; i<CORNERS_OF_SIDE(felem,side); i++)
      vp[i] = NVECTOR(CORNER_OF_SIDE_PTR(felem,side,i));

    /* sort vps */
    do
    {
      k = 0;
      for (i=1; i<CORNERS_OF_SIDE(felem,side); i++)
      {
        if (PVID(vp[i-1]) > PVID(vp[i]))
        {
          VECTOR *tmp = vp[i-1];
          vp[i-1] = vp[i];
          vp[i] = tmp;
          k = 1;
        }
      }
    }
    while (k!=0);

    /* store vps */
    for (i=0; i<CORNERS_OF_SIDE(felem,side); i++)
    {
      coordlist[nn].vp[i] = vp[i];
    }
    coordlist[nn].n = CORNERS_OF_SIDE(felem,side);
    break;
  }
#endif

  default :
    assert(0);
  }
  return(0);
}

static INT MatchingPeriodicEntries (PERIODIC_ENTRIES *coordlist, INT i, INT j)
{
  if (coordlist[i].n==0 && coordlist[j].n==0)
  {
    DOUBLE diff;

    V_DIM_EUKLIDNORM_OF_DIFF(coordlist[i].coord,coordlist[j].coord,diff);
    if (diff < SMALL_DOUBLE) return(1);
    else return(0);
  }
  else
  {
    INT k;
    if (coordlist[i].n != coordlist[j].n) return(0);

    for (k=0; k<coordlist[i].n; k++)
    {
      if (PVID(coordlist[i].vp[k]) != PVID(coordlist[j].vp[k])) return(0);
    }
    return (1);
  }

  return(0);
}

#ifdef ModelP
static int GetMatchingProcs (PERIODIC_ENTRIES *coordlist, INT i, int *np, int *theprocs);
#endif

INT NS_DIM_PREFIX SetPerVecVOBJECT(GRID *g)
{
  NODE *n;

  for (n=PFIRSTNODE(g); n!=NULL; n=SUCCN(n))
  {
    if (VOBJECT(NVECTOR(n)) == NULL)
    {
      VOBJECT(NVECTOR(n)) = (GEOM_OBJECT *)n;
    } else {
      if (PRIO((NODE *)VOBJECT(NVECTOR(n)))<PRIO(n))
        VOBJECT(NVECTOR(n)) = (void*)n;
    }
  }

        #ifdef Debug
  for (n=PFIRSTNODE(g); n!=NULL; n=SUCCN(n))
  {
    assert(VOBJECT(NVECTOR(n)) != NULL);
  }
        #endif

  return(0);
}

static void PrintListEntry (int i, PERIODIC_ENTRIES *coordlist)
{
  INT j;

  UserWriteF(PFMT " i=%d n=%d",me,i,coordlist[i].n);
  for (j=0; j<coordlist[i].n; j++)
    UserWriteF(" %08x",PVID(coordlist[i].vp[j]));
  UserWriteF("\n");

  return;
}

#ifdef ModelP
static void CountNTpls (GRID *g, INT nn, PERIODIC_ENTRIES *coordlist, int *ntpls)
{
  INT i;
  int j,theprocs[MAX_PERIODIC_PROCS],np;

  for (i=0; i<nn; i++)
  {
    assert(VOBJECT(NVECTOR(coordlist[i].node))!=NULL);

    if (coordlist[i].node != (NODE *)VOBJECT(NVECTOR(coordlist[i].node))) continue;

    np = 0;
    if (GetMatchingProcs(coordlist,i,&np,theprocs)) assert(0);

    for (j=0; j<np; j++)
    {
      ntpls[theprocs[j]]++;
    }
  }

  return;
}

static void MergeNTplsMax (int *gntpls, int *lntpls)
{
  int i;
  for (i=0; i<procs*procs; i++) gntpls[i] = MAX(gntpls[i],lntpls[i]);
  return;
}

static void AddMyNTpls (int *gntpls, int *ntpls)
{
  int i;
  for (i=0; i<procs; i++) gntpls[me*procs+i] = ntpls[i];
        #ifdef Debug
  if (0)
  {
    UserWriteF("AddMyNTpls:");
    for (i=0; i<procs; i++) UserWriteF(" %4d",ntpls[i]);
    UserWriteF("\n");
  }
        #endif
  return;
}

static void CommNTpls(GRID *g, int *send_ntpls, int *recv_ntpls)
{
  int l;
  INT n,MarkKey;
  int *gntpls,*lntpls;

  MarkTmpMem(MGHEAP(MYMG(g)),&MarkKey);

  gntpls = (int *)GetTmpMem(MGHEAP(MYMG(g)),procs*procs*sizeof(int),MarkKey);
  assert(gntpls!=NULL);
  memset(gntpls,0,procs*procs*sizeof(int));
  lntpls = (int *)GetTmpMem(MGHEAP(MYMG(g)),procs*procs*sizeof(int),MarkKey);
  assert(lntpls!=NULL);
  memset(lntpls,0,procs*procs*sizeof(int));

  /* construct global tupel array */
  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,lntpls,procs*procs*sizeof(int));
    MergeNTplsMax(gntpls,lntpls);
  }
  AddMyNTpls(gntpls,send_ntpls);
  Concentrate(gntpls,procs*procs*sizeof(int));
  Broadcast(gntpls,procs*procs*sizeof(int));
        #ifdef Debug
  if (1)
  {
    INT i,j;

    UserWriteF("CommNTplsMatrix:\n");

    for (i=0; i<procs; i++)
    {
      UserWriteF(PFMT,me);
      for (j=0; j<procs; j++) UserWriteF(" %4d",gntpls[i*procs+j]);
      UserWriteF("\n");
    }
  }
        #endif

  /* fill recv_ntpls array with recv counts */
  for (l=0; l<procs; l++) recv_ntpls[l] = gntpls[l*procs+me];

  ReleaseTmpMem(MGHEAP(MYMG(g)),MarkKey);

  return;
}

static void PrintIdTpl (int i, int j, IDTPL *idtpl)
{
  int k;

  UserWriteF(PFMT " j=%d n=%08x",i,j,(*idtpl).tpl[0]);
  for (k=0; k<idtpl->tpl[0]; k++)
    UserWriteF(" %08x",idtpl->tpl[k+1]);
  UserWriteF("\n");

  return;
}

static void CpyIDTpl (IDTPL **send_tpls, INT *send_tplcur, int proc, PERIODIC_ENTRIES *coordlist, int i)
{
  int j;

  send_tpls[proc][send_tplcur[proc]].tpl[0] = coordlist[i].n;

  for (j=0; j<coordlist[i].n; j++)
  {
    send_tpls[proc][send_tplcur[proc]].tpl[j+1] = PVID(coordlist[i].vp[j]);
    if (0)
      UserWriteF("CpyIDTpl copying: from=%08x to=%08x\n",
                 PVID(coordlist[i].vp[j]),send_tpls[proc][send_tplcur[proc]].tpl[j+1]);
  }

  /*
          UserWriteF("CpyIDTpl: ");
          PrintIdTpl(proc,send_tplcur[proc],&(send_tpls[proc][send_tplcur[proc]]));
   */

  send_tplcur[proc]++;
}

static void CommTpls (GRID *g, INT nn, PERIODIC_ENTRIES *coordlist, int *send_ntpls, IDTPL **send_tpls, int *recv_ntpls, IDTPL **recv_tpls, INT MarkKey)
{
  int i,j;

        #ifdef Debug
  if (0)
  {
    UserWriteF("Send_NTpls:");
    for (i=0; i<procs; i++) UserWriteF(" %4d",send_ntpls[i]);
    UserWriteF("\n");
    UserWriteF("Recv_NTpls:");
    for (i=0; i<procs; i++) UserWriteF(" %4d",recv_ntpls[i]);
    UserWriteF("\n");
  }
        #endif

  /* allocate send/recv tpls */
  for (i=0; i<procs; i++)
  {
    if (send_ntpls[i]>0 && recv_ntpls[i]>0)
    {
      send_tpls[i] = (IDTPL *)GetTmpMem(MGHEAP(MYMG(g)),send_ntpls[i]*sizeof(IDTPL),MarkKey);
      assert(send_tpls[i] != NULL);

      recv_tpls[i] = (IDTPL *)GetTmpMem(MGHEAP(MYMG(g)),recv_ntpls[i]*sizeof(IDTPL),MarkKey);
      assert(recv_tpls[i] != NULL);
    }
  }

  /* fill send_tpls */
  {
    INT *send_tplcur;
    int np,theprocs[MAX_PERIODIC_PROCS];

    send_tplcur = (INT *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(INT),MarkKey);
    assert(send_tplcur!=NULL);
    memset(send_tplcur,0,procs*sizeof(INT));

    for (i=0; i<nn; i++)
    {
      assert(VOBJECT(NVECTOR(coordlist[i].node))!=NULL);

      if (coordlist[i].node != (NODE *)VOBJECT(NVECTOR(coordlist[i].node))) continue;

      np = 0;
      if (GetMatchingProcs(coordlist,i,&np,theprocs)) assert(0);

      for (j=0; j<np; j++)
      {
        int p = theprocs[j];
        if (recv_ntpls[p] <= 0) continue;

        CpyIDTpl(send_tpls,send_tplcur,p,coordlist,i);
      }
    }
        #ifdef Debug
    /* check tuple sizes */
    for (i=0; i<procs; i++)
      if (send_ntpls[i]>0 && recv_ntpls[i]>0)
        assert(send_ntpls[i]==send_tplcur[i]);
        #endif
  }

  /* communicate IDTPLs */
  {
    VChannelPtr *mych;
    msgid *recv_msg, *send_msg;
    char *com_stat;
    int nc,rc,error;

    /* establish channels */
    mych = (VChannelPtr *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(VChannelPtr),MarkKey);
    assert(mych!=NULL);
    recv_msg = (msgid *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(msgid),MarkKey);
    assert(recv_msg!=NULL);
    send_msg = (msgid *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(msgid),MarkKey);
    assert(send_msg!=NULL);
    com_stat = (char *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(char),MarkKey);
    assert(com_stat!=NULL);
    memset(com_stat,0,procs*sizeof(char));

    nc = 0;
    for (i=0; i<procs; i++)
    {
      if (send_ntpls[i]>0 && recv_ntpls[i]>0)
      {
        mych[i] = ConnASync(i,3917);
        nc++;
      }
    }

    /* poll connect */
    rc = 0;
    while (rc < nc)
    {
      for (i=0; i<procs; i++)
      {
        if (send_ntpls[i]>0 && recv_ntpls[i]>0)
        {
          if ((com_stat[i] & 1) == 0)
            if (InfoAConn(mych[i]) > 0)
            {
              com_stat[i] |= 1;
              rc++;
            }
        }
      }
    }

    /* send/recv send/recv_tpls */
    for (i=0; i<procs; i++)
    {
      if (send_ntpls[i]>0 && recv_ntpls[i]>0)
      {
        recv_msg[i] = RecvASync(mych[i],(void *)recv_tpls[i],recv_ntpls[i]*sizeof(IDTPL),&error);
        send_msg[i] = SendASync(mych[i],(void *)send_tpls[i],send_ntpls[i]*sizeof(IDTPL),&error);
      }
    }

    /* poll send/recv */
    rc = 0;
    while (rc < 2*nc)
    {
      for (i=0; i<procs; i++)
      {
        if (send_ntpls[i]>0 && recv_ntpls[i]>0)
        {
          if ((com_stat[i] & 2) == 0)
            if (InfoASend(mych[i],send_msg[i]) > 0)
            {
              com_stat[i] |= 2;
              rc++;
            }

          if ((com_stat[i] & 4) == 0)
            if (InfoARecv(mych[i],recv_msg[i]) > 0)
            {
              com_stat[i] |= 4;
              rc++;
            }
        }
      }
    }

    /* disconnect channels */
    for (i=0; i<procs; i++)
    {
      if (send_ntpls[i]>0 && recv_ntpls[i]>0)
      {
        DiscASync(mych[i]);
      }
    }

    /* poll disconnect */
    rc = 0;
    while (rc < nc)
    {
      for (i=0; i<procs; i++)
      {
        if (send_ntpls[i]>0 && recv_ntpls[i]>0)
        {
          if ((com_stat[i] & 8) == 0)
            if (InfoADisc(mych[i]) > 0)
            {
              com_stat[i] |= 8;
              rc++;
            }
        }
      }
    }

  }
        #ifdef Debug
  /*
          Synchronize();
   */
        #endif

  return;
}

static INT SelectCorProc (VECTOR **vp, INT pos, INT p, int *np, int *theprocs)
{
  int *proclist = PROCLIST(vp[pos]);
  int j,Proclist[MAX_PERIODIC_PROCS];

  /* copy proclist */
  j = 0;
  while (*proclist != -1)
  {
    Proclist[j] = proclist[0];
    Proclist[j+1] = proclist[1];

    proclist += 2;
    j+=2;
    assert(j < MAX_PERIODIC_PROCS);
  }
  Proclist[j] = -1;
  proclist = Proclist;

  proclist += 2;
  while (*proclist != -1)
  {
    PRINTDEBUG(gm,1,("pos%d " VINDEX_FMTX " proc=%d prio=%d\n",
                     pos,VINDEX_PRTX(vp[pos]),*proclist,*(proclist+1)));

    if (*proclist==p && MASTERPRIO(*(proclist+1)))
    {
      if (pos == 0)
      {
        PRINTDEBUG(gm,1,("FOUND np=%d newproc=%d\n", *np,*proclist))
        theprocs[*np] = p;
        (*np)++;
      }
      else
      if (SelectCorProc(vp,pos-1,p,np,theprocs))
        assert(0);
    }
    proclist += 2;
  }

  return(0);
}

static int GetMatchingProcs (PERIODIC_ENTRIES *coordlist, INT i, int *np, int *theprocs)
{
  INT pos = coordlist[i].n-1;
  int *proclist = PROCLIST(coordlist[i].vp[pos]);
  int j,Proclist[MAX_PERIODIC_PROCS];

  /* copy proclist */
  j = 0;
  while (*proclist != -1)
  {
    Proclist[j] = proclist[0];
    Proclist[j+1] = proclist[1];

    proclist += 2;
    j+=2;
    assert(j < MAX_PERIODIC_PROCS);
  }
  Proclist[j] = -1;
  proclist = Proclist;

  proclist += 2;
  while (*proclist != -1)
  {
    PRINTDEBUG(gm,1,("pos%d " VINDEX_FMTX " proc=%d prio=%d\n",
                     pos,VINDEX_PRTX(coordlist[i].vp[pos]),*proclist,*(proclist+1)));

    if (MASTERPRIO(*(proclist+1)))
    {
      if (pos == 0)
      {
        PRINTDEBUG(gm,1,("FOUND np=%d newproc=%d\n", *np,*proclist))
        theprocs[*np] = *proclist;
        (*np)++;
        assert(*np < MAX_PERIODIC_PROCS);
      }
      else
      if (SelectCorProc(coordlist[i].vp,pos-1,*proclist,np,theprocs))
        assert(0);
    }

    proclist += 2;
  }

  return(0);
}

static INT Identify_PeriodicVectorX (PERIODIC_ENTRIES *coordlist, INT i, INT p)
{
  int j;

  PRINTDEBUG(gm,1,("IPVX:" VINDEX_FMTX " p=%d\n",
                   VINDEX_PRTX(NVECTOR(coordlist[i].node)),p))

  for (j=0; j<coordlist[i].n; j++)
  {
    ASSERT(NVECTOR(coordlist[i].node)!=NULL);
    ASSERT(coordlist[i].vp[j]!=NULL);
    PRINTDEBUG(gm,1,(" GID=%08x",
                     GID(coordlist[i].vp[j])))

    DDD_IdentifyNumber(PARHDR(NVECTOR(coordlist[i].node)),p,
                       GID(coordlist[i].vp[j]));
  }
        #ifdef Debug
  {
    DOUBLE *cv = CVECT(MYVERTEX(coordlist[i].node));
    PRINTDEBUG(gm,1,(" c %lf %lf %lf \n",cv[0],cv[1],cv[2]))
  }
        #endif

  return(0);
}

static INT Identify_PeriodicVector (PERIODIC_ENTRIES *coordlist, INT i)
{
  int p,j,theprocs[MAX_PERIODIC_PROCS],np;

  PRINTDEBUG(gm,1,("%d IPV: GetMatchingProcs\n",i));

  np = 0;
  if (GetMatchingProcs(coordlist,i,&np,theprocs)) return(1);

        #ifdef Debug
  if(0)
    if (coordlist[i].n != 2)
      return(0);
        #endif

  PRINTDEBUG(gm,1,("IPV:" VINDEX_FMTX " np=%d\n",
                   VINDEX_PRTX(NVECTOR(coordlist[i].node)),np))

  for (p=0; p<np; p++)
  {
    PRINTDEBUG(gm,1,("p=%d",theprocs[p]))

    for (j=0; j<coordlist[i].n; j++)
    {
      ASSERT(NVECTOR(coordlist[i].node)!=NULL);
      ASSERT(coordlist[i].vp[j]!=NULL);
      PRINTDEBUG(gm,1,(" GID=%08x",
                       GID(coordlist[i].vp[j])))

      if (0)
        DDD_IdentifyObject(PARHDR(NVECTOR(coordlist[i].node)),theprocs[p],
                           PARHDR(coordlist[i].vp[j]));
      else
        DDD_IdentifyNumber(PARHDR(NVECTOR(coordlist[i].node)),theprocs[p],
                           GID(coordlist[i].vp[j]));
    }
                #ifdef Debug
    {
      DOUBLE *cv = CVECT(MYVERTEX(coordlist[i].node));
      PRINTDEBUG(gm,1,(" c %lf %lf %lf \n",cv[0],cv[1],cv[2]))
    }
                #endif
  }

  return(0);
}

static int CompTpls (IDTPL *tpl, PERIODIC_ENTRIES *cle)
{
  int i;

  if (tpl->tpl[0] > cle->n) return(1);
  if (tpl->tpl[0] < cle->n) return(-1);

  for (i=0; i<cle->n; i++)
  {
    if (tpl->tpl[i+1] > PVID(cle->vp[i])) return(1);
    if (tpl->tpl[i+1] < PVID(cle->vp[i])) return(-1);
  }

  return(0);
}


static void PrintIdTpls (int *idntpls, IDTPL **idtpls)
{
  int i,j,k;

  for (i=0; i<procs; i++)
  {
    if (idntpls[i] <= 0) continue;

    for (j=0; j<idntpls[i]; j++)
    {
      PrintIdTpl(i,j,idtpls[i]+j);
    }
  }

  return;
}

static void IdentListX (GRID *g, INT nn, PERIODIC_ENTRIES *coordlist, int *recv_ntpls, IDTPL **recv_tpls, INT MarkKey)
{
  INT i,j;
  int theprocs[MAX_PERIODIC_PROCS],np;
  INT *recv_tplscur,*nidv,*njpv;

  recv_tplscur = (INT *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(INT),MarkKey);
  assert(recv_tplscur!=NULL);
  memset(recv_tplscur,0,procs*sizeof(INT));
  nidv = (INT *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(INT),MarkKey);
  assert(nidv!=NULL);
  memset(nidv,0,procs*sizeof(INT));
  njpv = (INT *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(INT),MarkKey);
  assert(njpv!=NULL);
  memset(njpv,0,procs*sizeof(INT));

  PRINTDEBUG(gm,1,("IPV: GetMatchingProcs\n"));

  DDD_IdentifyBegin();

  if (0) PrintIdTpls(recv_ntpls,recv_tpls);

  for (i=0; i<nn; i++)
  {
    assert(VOBJECT(NVECTOR(coordlist[i].node))!=NULL);

    if (coordlist[i].node != (NODE *)VOBJECT(NVECTOR(coordlist[i].node))) continue;

    np = 0;
    if (GetMatchingProcs(coordlist,i,&np,theprocs)) assert(0);

    for (j=0; j<np; j++)
    {
      int p = theprocs[j];

      if (recv_ntpls[p] <= 0) continue;

      if (0) PrintListEntry(i,coordlist);

      while (recv_tplscur[p]<recv_ntpls[p] &&
             CompTpls(&(recv_tpls[p][recv_tplscur[p]]),coordlist+i)<0)
      {
        UserWriteF(PFMT " skipping entry:",me);
        PrintIdTpl(p,recv_tplscur[p],recv_tpls[p]+recv_tplscur[p]);
        PrintListEntry(i,coordlist);

        recv_tplscur[p]++;
        njpv[p]++;
      }

      if (recv_tplscur[p]<recv_ntpls[p] &&
          CompTpls(&(recv_tpls[p][recv_tplscur[p]]),coordlist+i)==0)
      {
        /*
                                        UserWriteF(PFMT " ident entry:",me);
                                        PrintIdTpl(p,recv_tplscur[p],recv_tpls[p]+recv_tplscur[p]);
                                        PrintListEntry(i,coordlist);
         */

        Identify_PeriodicVectorX(coordlist,i,p);
        recv_tplscur[p]++;
        nidv[p]++;
      }
    }
  }

  if (1)
  {
    DDD_IdentifyEnd();
    DDD_IFRefreshAll();
  }

#ifdef Debug
  for (j=0; j<procs; j++)
  {
    if (nidv[j] != recv_ntpls[j])
    {
      UserWriteF(PFMT "SKIPS p=%d nrecv=%d nidv=%d njpv=%d\n",
                 me,j,recv_ntpls[j],nidv[j],njpv[j]);
    }
  }
  /*
     Synchronize();
     Synchronize();
   */
#endif

  return;
}

static void IdentList (INT nn, PERIODIC_ENTRIES *coordlist)
{
  INT i;

  DDD_IdentifyBegin();

  for (i=0; i<nn; i++)
  {
    assert(VOBJECT(NVECTOR(coordlist[i].node))!=NULL);

    if (coordlist[i].node != (NODE *)VOBJECT(NVECTOR(coordlist[i].node))) continue;

    Identify_PeriodicVector(coordlist,i);
  }

  if (1)
  {
    DDD_IdentifyEnd();
    DDD_IFRefreshAll();
  }

  return;
}
#endif

static INT CheckPerVectors(GRID *theGrid)
{
  VECTOR *theVector;

  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
  {
    if (GID(theVector) == 0x00122801)
      PRINTDEBUG(gm,0,(PFMT " FOUND disposed vec=" VINDEX_FMTX,
                       me,VINDEX_PRTX(theVector)));
  }

  return(0);
}

INT NS_DIM_PREFIX Grid_GeometricToPeriodic (GRID *g)
{
  NODE *node;
  PERIODIC_ENTRIES *coordlist = NULL;
  INT MarkKey,nn,i,level;
  INT heapmem = 0;

  level = GLEVEL(g);

  PRINTDEBUG(gm,1,("Grid_GeometricToPeriodic\n"))

  if (PeriodicBoundaryInfo == NULL)
  {
    UserWriteF("Grid_GeometricToPeriodic: no function *PeriodicBoundaryInfo\n");
    return(GM_OK);
  }

        #ifndef ModelP
  /* set PVID of vectors */
  i = 0;
  for (node=FIRSTNODE(g); node!=NULL; node=SUCCN(node))
  {
    if (node == (NODE *)VOBJECT(NVECTOR(node)))
    {
      PVID(NVECTOR(node)) = i;
      i++;
    }
  }
        #endif

  /* count vector-node pairs to identify */
  nn = 0;
  for (node=FIRSTNODE(g); node!=NULL; node=SUCCN(node))
  {
    VERTEX *vtx;
    DOUBLE_VECTOR own_coord, periodic_coords[MAX_PERIODIC_OBJ];
    INT n,periodic_ids[MAX_PERIODIC_OBJ];

    vtx = MYVERTEX(node);

    if (OBJT(vtx)!=BVOBJ) continue;

    /* only boundary vertices */
    n = 0;
    if ((*PeriodicBoundaryInfo)(vtx,&n,periodic_ids,own_coord,periodic_coords)) nn += n;
  }

  /* allocate arrays */
  if (heapmem)
  {
    MarkTmpMem(MGHEAP(MYMG(g)),&MarkKey);

    coordlist = (PERIODIC_ENTRIES *)GetTmpMem(MGHEAP(MYMG(g)),nn*sizeof(PERIODIC_ENTRIES),MarkKey);
  }
  else
    coordlist = (PERIODIC_ENTRIES *)malloc(nn*sizeof(PERIODIC_ENTRIES));
  if (coordlist == NULL) assert(0);

  nn=0;
  for (node=FIRSTNODE(g); node!=NULL; node=SUCCN(node))
  {
    VERTEX *vtx;
    DOUBLE_VECTOR own_coord, periodic_coords[MAX_PERIODIC_OBJ];
    DOUBLE diff;
    INT n,periodic_ids[MAX_PERIODIC_OBJ];

    vtx=MYVERTEX(node);

    if (OBJT(vtx)!=BVOBJ) continue;

    /* only boundary vertices */
    n = 0;
    if ((*PeriodicBoundaryInfo)(vtx,&n,periodic_ids,own_coord,periodic_coords))
    {
      INT i;

      for (i=0; i<n; i++)
      {
        DOUBLE_VECTOR *coord;

        coord = MIN_COORD(&own_coord,&periodic_coords[i]);

        V_DIM_EUKLIDNORM_OF_DIFF(own_coord,periodic_coords[i],diff);

        /* equal coordinates nothing to do */
        if (diff < SMALL_DOUBLE) continue;

        PRINTDEBUG(gm,1,("coordlist identify v=%d c0 %lf %lf %lf c1 %lf %lf %lf\n",
                         ID(vtx),own_coord[0],own_coord[1],own_coord[2],
                         periodic_coords[i][0],periodic_coords[i][1],periodic_coords[i][2]))

        V_DIM_COPY(*coord,coordlist[nn].coord);

        PRINTDEBUG(gm,1,("                    coord %lf %lf %lf\n",
                         coordlist[nn].coord[0],coordlist[nn].coord[1],coordlist[nn].coord[2]));

        coordlist[nn].node = node;
        coordlist[nn].periodic_id = periodic_ids[i];

        if (level > 0)
        {
          InsertVIDs(nn,coordlist,node);
        }
        else
          coordlist[nn].n = 0;
        nn++;
      }
    }
  }

  PRINTDEBUG(gm,1,("Grid_GeometricToPeriodic identify nn=%d, level %d\n",nn,level))

  /* sort list */
  qsort(coordlist,nn,sizeof(PERIODIC_ENTRIES),sort_entries);

  IFDEBUG(gm,1)
  for (i=0; i<nn; i++)
  {
    UserWriteF("%d perid=%d vtx=%08d node=" ID_FMTX " v="
               VINDEX_FMTX " c %lf %lf %lf %d\n",
               i,coordlist[i].periodic_id,
               ID(MYVERTEX(coordlist[i].node)),
               ID_PRTX(coordlist[i].node),VINDEX_PRTX(NVECTOR(coordlist[i].node)),
               coordlist[i].coord[0],
               coordlist[i].coord[1],
               coordlist[i].coord[2],
               coordlist[i].n
               );
    if (coordlist[i].n>0)
    {
      PrintListEntry(i,coordlist);
    }
  }
  ENDDEBUG

        #ifndef ModelP
  assert(nn%2==0);
        #endif

  /* initialize PVCOUNT */
  GridSetPerVecCount(g);

  /* identify */
  for (i=0; i<nn-1; i++)
  {
    if (MatchingPeriodicEntries(coordlist,i,i+1))
    {
      INT mode = 0;
      DOUBLE diff;

      V_DIM_EUKLIDNORM_OF_DIFF(coordlist[i].coord,CVECT(MYVERTEX(coordlist[i].node)),diff);
      if (diff < SMALL_DOUBLE) mode = 1;

      if (mode)
      {
        if (ModifyVectorPointer(g,coordlist,i,i+1))
          return (1);
      }
      else {
        if (ModifyVectorPointer(g,coordlist,i+1,i))
          return (1);
      }
      /* increase counter since matching pair already found -> step to next possible pair */
      i++;
    }
  }

  /* dispose vectors */
        #ifdef ModelP
  DDD_XferBegin();
    #ifdef DDDOBJMGR
  DDD_ObjMgrBegin();
    #endif
    #endif

  for (i=0; i<nn-1; i++)
  {
    if (MatchingPeriodicEntries(coordlist,i,i+1))
    {
      INT mode = 0;
      DOUBLE diff;

      V_DIM_EUKLIDNORM_OF_DIFF(coordlist[i].coord,CVECT(MYVERTEX(coordlist[i].node)),diff);
      if (diff < SMALL_DOUBLE) mode = 1;

      if (mode)
      {
        if (DisposeAndModVector(g,coordlist,i,i+1))
          return (1);
      }
      else {
        if (DisposeAndModVector(g,coordlist,i+1,i))
          return (1);
      }

      /* increase counter since matching pair already found -> step to next possible pair */
      i++;
    }
  }
        #ifdef ModelP
    #ifdef DDDOBJMGR
  DDD_ObjMgrEnd();
        #endif
  DDD_XferEnd();
  DDD_IFRefreshAll();
        #endif

  IFDEBUG(gm,1)
  UserWriteF("After Periodic Identification:\n");
  for (i=0; i<nn-1; i++)
  {
    if (MatchingPeriodicEntries(coordlist,i,i+1))
    {
      UserWriteF("%d v0=" VINDEX_FMTX " v1=" VINDEX_FMTX
                 " n=%d\n",i,
                 VINDEX_PRTX(NVECTOR(coordlist[i].node)),
                 VINDEX_PRTX(NVECTOR(coordlist[i+1].node)),
                 coordlist[i].n);
      i++;
    }
  }
  ENDDEBUG

        #ifdef ModelP
  /* identify periodic vectors between processors */
  if (GLEVEL(g)>0 && procs>1)
  {
    INT MarkKey,i;
    int *send_ntpls,*recv_ntpls;
    IDTPL **send_tpls,**recv_tpls;

    if (1) {
      /* count identify tupels per proc */
      MarkTmpMem(MGHEAP(MYMG(g)),&MarkKey);

      send_ntpls = (int *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(int),MarkKey);
      assert(send_ntpls!=NULL);
      memset(send_ntpls,0,procs*sizeof(int));

      recv_ntpls = (int *)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(int),MarkKey);
      assert(recv_ntpls!=NULL);
      memset(recv_ntpls,0,procs*sizeof(int));

      CountNTpls(g,nn,coordlist,send_ntpls);

      /* communicate local tupel list count */
      CommNTpls(g,send_ntpls,recv_ntpls);

      /* prepare, send and recv tupel lists */
      send_tpls = (IDTPL **)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(IDTPL *),MarkKey);
      assert(send_tpls!=NULL);
      for (i=0; i<procs; i++) send_tpls[i] = NULL;
      recv_tpls = (IDTPL **)GetTmpMem(MGHEAP(MYMG(g)),procs*sizeof(IDTPL *),MarkKey);
      assert(recv_tpls!=NULL);
      for (i=0; i<procs; i++) recv_tpls[i] = NULL;
      CommTpls(g,nn,coordlist,send_ntpls,send_tpls,recv_ntpls,recv_tpls,MarkKey);

      IdentListX(g,nn,coordlist,recv_ntpls,recv_tpls,MarkKey);
    }
    else {
      /* identify entries */
      IdentList(nn,coordlist);
    }

    ReleaseTmpMem(MGHEAP(MYMG(g)),MarkKey);
  }
        #endif

  if (heapmem)
    ReleaseTmpMem(MGHEAP(MYMG(g)),MarkKey);
  else
    free(coordlist);

  IFDEBUG(gm,1);
  GridSetPerVecCount(g);
  ENDDEBUG

  return (GM_OK);
}


/****************************************************************************/
/** \brief
   MG_GeometricToPeriodic - identify vectors on periodic boundaries

   SYNOPSIS:
   INT MG_GeometricToPeriodic (MULTIGRID *mg)


   .  mg - multigrid to work on

   DESCRIPTION:
   This function does all that is necessary to identify vectors on
   periodic boundaries. The node-vector pointers are unsymmetric.
   The vector points to one of the nodes, which points to the single
   periodic vector.

   @return <ul>
   INT
   .n   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX MG_GeometricToPeriodic (MULTIGRID *mg, INT fl, INT tl)
{
  INT level;

  for (level=fl; level<=tl; level++)
  {
    GRID *g = GRID_ON_LEVEL(mg,level);

    if (Grid_GeometricToPeriodic(g)) return(GM_ERROR);
  }

  return (GM_OK);
}

/****************************************************************************/
/** \brief
   Grid_CheckPeriodicity - check if periodicity is correct

   SYNOPSIS:
   INT Grid_CheckPeriodicity (GRID *grid)


 * @param   grid - grid to work on

   DESCRIPTION:
   This function checks all connections at periodic boundaries.

   @return <ul>
   INT
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
   </ul> */
/****************************************************************************/

INT NS_DIM_PREFIX Grid_CheckPeriodicity (GRID *grid)
{
  NODE *node;
  INT nn;

  if (PeriodicBoundaryInfo == NULL)
  {
    UserWriteF("Grid_CheckPeriodicity: no function *PeriodicBoundaryInfo\n");
    return(GM_OK);
  }

  PRINTDEBUG(gm,1,("Grid_CheckPeriodicity p=%d: level=%d\n",me,GLEVEL(grid)));

  nn=0;
  for (node=PFIRSTNODE(grid); node!=NULL; node=SUCCN(node)) {
    VERTEX *vtx;
    DOUBLE_VECTOR own_coord, periodic_coords[MAX_PERIODIC_OBJ];
    INT n,periodic_ids[MAX_PERIODIC_OBJ];

    vtx = MYVERTEX(node);

    if (vtx==NULL) continue;

    if (OBJT(vtx)!=BVOBJ) continue;

    /* only boundary vertices */
    n=0;
    if ((*PeriodicBoundaryInfo)(vtx,&n,periodic_ids,own_coord,periodic_coords)) {
      MATRIX *mat;

      /* print all information available */
      PRINTDEBUG(gm,0,("#%8d: per_node GID= %08x VGID= %08x PRIO=%d coord: ( %f %f %f )\n",nn,GID(node),GID(NVECTOR(node)), PRIO(NVECTOR(node)), own_coord[0],own_coord[1],own_coord[2]));
      nn++;

      if (((NODE *)VOBJECT(NVECTOR(node)))==NULL) {
        PRINTDEBUG(gm,0,("nodevector has no node!!!!!!\n"));
        return (GM_ERROR);
      }

      if (VSTART(NVECTOR(node))==NULL) {
        PRINTDEBUG(gm,0,("no matrix! -> continue\n"));
        continue;
      }

      PRINTDEBUG(gm,0,("connected to: \n"));
      for (mat=MNEXT(VSTART(NVECTOR(node))); mat!=NULL; mat=MNEXT(mat)) {
        if (MDEST(mat)==NULL) {
          PRINTDEBUG(gm,0,("no destination vector -> continue\n"));
          continue;
        }
        PRINTDEBUG(gm,0,("\tGID=%08x PRIO=%d\n",GID(MDEST(mat)),PRIO(MDEST(mat))));
      }
    }
  }

  return (GM_OK);
}

char pbuf[256];

static INT ListProclist (int *proclist)
{
  while (*proclist != -1)
  {
    sprintf(pbuf+strlen(pbuf),"%4d-%d ",proclist[0],proclist[1]);
    proclist += 2;
  }

  return(GM_OK);
}

static INT ListPeriodicNodeAndVec (GRID *g, INT vgid)
{
  VECTOR *v;
  NODE *n,*nref;
  int *proclist;
  INT found = 0;

  sprintf(pbuf," ");
  for (v=PFIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
  {
    if (vgid != GID(v)) continue;
    nref = (NODE*) VOBJECT(v);
    found++;

    /*                  sprintf(pbuf+strlen(pbuf),"v=" VINDEX_FMTX " ", VINDEX_PRTX(v)); */
    /*                  sprintf(pbuf+strlen(pbuf),"LEVEL %d v=%08x/%d  ",GLEVEL(g), GID(v),PRIO(v)); */

    proclist = PROCLIST(v);
    /*                  ListProclist(proclist); */
    /*                  sprintf(pbuf+strlen(pbuf),"  VNEW=%d  ",VNEW(v)); */
  }

  for (n=PFIRSTNODE(g); n!=NULL; n=SUCCN(n))
  {
    DOUBLE *cv;

    if (vgid != GID(NVECTOR(n))) continue;

    cv = CVECT(MYVERTEX(n));

                #ifdef __TWODIM__
    sprintf(pbuf+strlen(pbuf),"LEVEL %d c %g %g ",GLEVEL(g),cv[0],cv[1]);
                #else
    sprintf(pbuf+strlen(pbuf),"c %g %g %g ",cv[0],cv[1],cv[2]);
                #endif
    sprintf(pbuf+strlen(pbuf),"v=%08x/%d  ",GID(NVECTOR(n)),PRIO(NVECTOR(n)));
    proclist = PROCLIST(NVECTOR(n));
    ListProclist(proclist);

    if (n == nref) sprintf(pbuf+strlen(pbuf),"	X ");
    else sprintf(pbuf+strlen(pbuf),"	  ");

    /*                  sprintf(pbuf+strlen(pbuf),"n=" ID_FMTX " ",ID_PRTX(n)); */
    sprintf(pbuf+strlen(pbuf),"n=%08x/%d  ",GID(n),PRIO(n));

    proclist = PROCLIST(n);
    ListProclist(proclist);
    sprintf(pbuf+strlen(pbuf),"\n");
    found++;
  }
  if (found == 0) sprintf(pbuf+strlen(pbuf),"NOT FOUND\n");

  UserWriteF("%s\n",pbuf);

  return(found);
}

static INT Grid_ListPeriodicPos (GRID *g, DOUBLE_VECTOR pos)
{
  DOUBLE tol[3] = {SMALL_DOUBLE,SMALL_DOUBLE,SMALL_DOUBLE};
  NODE *n = FindNodeFromPosition(g,pos,tol);
  INT i,found;
  INT vgid;

  if (n != NULL)
    vgid = GID(NVECTOR(n));
  else
    vgid = -1;
  ASSERT(vgid==UG_GlobalMaxINT(vgid) || vgid==-1);
  vgid = UG_GlobalMaxINT(vgid);

  if (vgid == -1)
  {
    UserWriteF("NOT FOUND\n");
    return (GM_OK);
  }

  found = 0;
  for (i=0; i<procs; i++)
  {
                #ifdef ModelP
    /*
                    Synchronize();
                    fflush(stdout);
     */
                #endif
    if (me!=i) continue;

    found += ListPeriodicNodeAndVec(g,vgid);
  }
        #ifdef ModelP
  fflush(stdout);
        #endif
  found = UG_GlobalSumINT(found);

  if (me == master)
    UserWriteF("FOUND %d periodic nodes for this position\n",found);

  return (GM_OK);
}

/****************************************************************************/
/** \brief
   MG_ListPeriodicPos - list periodic positions and node-vector(-proc) info

 * @param   mg - multigrid to work on

   This function list periodic positions including vector and information
   of all periodic boundaries. In parallel also proclist are printed.
   This should give a clear view onto the periodic data structure.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
 */
/****************************************************************************/

INT NS_DIM_PREFIX MG_ListPeriodicPos (MULTIGRID *mg, INT fl, INT tl, DOUBLE_VECTOR pos)
{
  INT level;

  if (me == master)
  {
#ifdef __THREEDIM__
    UserWriteF("position is %lf %lf %lf\n",pos[0],pos[1],pos[2]);
#else
    UserWriteF("position is %lf %lf\n",pos[0],pos[1]);
#endif
  }
  for (level=fl; level<=tl; level++)
  {
    GRID *g = GRID_ON_LEVEL(mg,level);

    /*                  if (me == master) */
    UserWriteF("LEVEL %4d\n",level);
    if (Grid_ListPeriodicPos(g,pos)) return(GM_ERROR);
  }

  return (GM_OK);
}
#endif

/****************************************************************************/
/** \brief Do all that is necessary to complete the coarse grid

 * @param   id - the id of the block to be allocated

   This function does all that is necessary to complete the coarse grid.
   Finally the MG_COARSE_FIXED flag is set.

   @return <ul>
   <li>   GM_OK if ok </li>
   <li>   GM_ERROR if error occured </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX FixCoarseGrid (MULTIGRID *theMG)
{
  if (MG_COARSE_FIXED(theMG)) return (GM_OK);

  /** \todo (HRR 971031): check that before check-in!
     if (FinishGrid(theMG))
          REP_ERR_RETURN (GM_ERROR);*/

  /** \todo (HRR 971031): remove if above works */
  if (SetSubdomainIDfromBndInfo(theMG))
    REP_ERR_RETURN (GM_ERROR);

  /* set this flag here because it is checked by CreateAlgebra */
  if (CreateAlgebra(theMG) != GM_OK)
    REP_ERR_RETURN (GM_ERROR);

#ifdef __PERIODIC_BOUNDARY__
  if (MG_GeometricToPeriodic(theMG,0,0))
    REP_ERR_RETURN (GM_ERROR);
#endif

  /* here all temp memory since CreateMultiGrid is released */
  ReleaseTmpMem(MGHEAP(theMG),MG_MARK_KEY(theMG));
  MG_MARK_KEY(theMG) = 0;

  return (GM_OK);
}

/****************************************************************************/
/** \brief Init what is necessary
 *
 *  This function initializes the grid manager.
 *
 *  @return <ul>
 *     <li> GM_OK if ok </li>
 *     <li> > 0 line in which error occured </li>
 *  </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitUGManager ()
{
  INT i;

  theGenMGUDM = (VIRT_HEAP_MGMT*)malloc(SIZEOF_VHM);
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

INT NS_DIM_PREFIX ExitUGManager ()
{
  free(theGenMGUDM);

  return 0;
}

/* nur temporaer zum Debuggen drin (Christian Wrobel): */
/* TODO: entfernen nach Debuggphase */
char *PrintElementInfo (ELEMENT *theElement,INT full)
{
  static char out[2000];
  char tmp[200];
  char etype[10];
  char ekind[8];
  int i,j;
  ELEMENT *SonList[MAX_SONS];

  if (theElement==NULL)
  {
    printf( "PrintElementInfo: element == NULL\n");
    return (NULL);
  }

  if (DIM==2)
    switch (TAG(theElement))
    {
    case TRIANGLE :          strcpy(etype,"TRI"); break;
    case QUADRILATERAL :     strcpy(etype,"QUA"); break;
    default :                strcpy(etype,"???"); break;
    }
  else
    switch (TAG(theElement))
    {
    case TETRAHEDRON :       strcpy(etype,"TET"); break;
    case PYRAMID :           strcpy(etype,"PYR"); break;
    case PRISM :             strcpy(etype,"PRI"); break;
    case HEXAHEDRON :        strcpy(etype,"HEX"); break;
    default :                strcpy(etype,"???"); break;
    }
  switch (ECLASS(theElement))
  {
  case YELLOW_CLASS :      strcpy(ekind,"YELLOW "); break;
  case GREEN_CLASS :       strcpy(ekind,"GREEN  "); break;
  case RED_CLASS :         strcpy(ekind,"RED    "); break;
  default :                strcpy(ekind,"???    "); break;
  }
  if(full)
    sprintf(out,"ELEMID=" EID_FFMTE " %5s %5s CTRL=%8lx CTRL2=%8lx REFINE=%2d MARK=%2d LEVEL=%2d",
            EID_PRTE(theElement),ekind,etype,
            (long)CTRL(theElement),(long)FLAG(theElement),REFINE(theElement),MARK(theElement),LEVEL(theElement));
  else
    sprintf(out,"ELEMID=" EID_FFMTE, EID_PRTE(theElement));

  if (COARSEN(theElement)) strcat(out," COARSEN");
  strcat(out,"\n");
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
                #ifdef __TWODIM__
    sprintf(tmp,"    N%d=" ID_FMTX " x=%g  y=%g\n",
            i,
            ID_PRTX(CORNER(theElement,i)),
            CVECT(MYVERTEX(CORNER(theElement,i)))[0],
            CVECT(MYVERTEX(CORNER(theElement,i)))[1]
            );
                #endif
                #ifdef __THREEDIM__
    sprintf(tmp,"    N%d=" ID_FMTX " x=%g  y=%g z=%g\n",
            i,
            ID_PRTX(CORNER(theElement,i)),
            CVECT(MYVERTEX(CORNER(theElement,i)))[0],
            CVECT(MYVERTEX(CORNER(theElement,i)))[1],
            CVECT(MYVERTEX(CORNER(theElement,i)))[2]
            );
                #endif
    strcat( out, tmp );
  }

  if (EFATHER(theElement))
  {
    sprintf(tmp,"    FA=" EID_FMTX "\n" ,EID_PRTX(EFATHER(theElement)));
    strcat( out, tmp );
  }
  else
    strcat( out,"    FA=NULL\n");

  if( full)
  {
    UserWriteF("  NSONS=%d\n",NSONS(theElement));
    if (GetAllSons(theElement,SonList)==0)
    {
      for (i=0; SonList[i] != NULL; i++)
      {
        sprintf(tmp,"    SON%d "EID_FMTX "\n" ,i,EID_PRTX(SonList[i]));
        strcat( out, tmp );

        for (j=0; j<CORNERS_OF_ELEM(SonList[i]); j++)
        {
                                        #ifdef __TWODIM__
          sprintf(tmp,"        N%d= "ID_FMTX " x=%g  y=%g\n",
                  j,
                  ID_PRTX(CORNER(SonList[i],j)),
                  CVECT(MYVERTEX(CORNER(SonList[i],j)))[0],
                  CVECT(MYVERTEX(CORNER(SonList[i],j)))[1]
                  );
                                        #endif
                                        #ifdef __THREEDIM__
          sprintf(tmp,"        N%d= "ID_FMTX " x=%g  y=%g z=%g\n",
                  j,
                  ID_PRTX(CORNER(SonList[i],j)),
                  CVECT(MYVERTEX(CORNER(SonList[i],j)))[0],
                  CVECT(MYVERTEX(CORNER(SonList[i],j)))[1],
                  CVECT(MYVERTEX(CORNER(SonList[i],j)))[2]
                  );
                                        #endif
          strcat( out, tmp );
        }
      }
    }
  }
  sprintf(tmp," key=%d\n", KeyForObject((KEY_OBJECT *)theElement) );
  strcat( out, tmp );

  if(full)
  {
    if (OBJT(theElement)==BEOBJ)
      strcat( out," boundary element\n" );
    else
      strcat( out," no boundary element\n" );

    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      for(j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
      {
                                #ifdef __TWODIM__
        sprintf(tmp,"    NODE[ID=%ld]: x=%g y=%g",
                (long)(ID(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j)))),
                CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))))[0],
                CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))))[1]);
                                #endif
                                #ifdef __THREEDIM__
        sprintf(tmp,"    NODE[ID=%ld]: x=%g y=%g z=%g",
                (long)(ID(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j)))),
                CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))))[0],
                CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))))[1],
                CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))))[2]);
                                #endif
        strcat( out, tmp );
      }
      strcat( out,"\n");
    }
  }
#ifdef ModelP
  /*UserWriteF(PFMT"%s", me,out );*/
  printf(PFMT "%s",me,out);
#else
  UserWrite(out);
#endif

  return out;
}
