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
#include "algebra.h"
#include "ugm.h"
#include "elements.h"
#include "shapes.h"

/* include refine because of macros accessed in list functions */
#include "refine.h"

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

static char buffer[256];                        /* general purpose text buffer			*/

static VIRT_HEAP_MGMT *theGenMGUDM; /* general user data space management	*/

static INT theMGVarID;                          /* env var ID for the multigrids		*/
static INT theDocDirID;                         /* env dir ID for the multigrids		*/

static INT UsedOBJT;                            /* for the dynamic OBJECT management	*/

/* used by OrderNodesInGrid */
static const INT *Order,*Sign;
static DOUBLE InvMeshSize;

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

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
  MULTIGRID *theMG;

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
void *GetMemoryForObject (MULTIGRID *theMG, INT size, INT type)
{
  void *obj = GetMemoryLocal(theMG, size, type);

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


/* in ModelP, this function creates local memory. */
void *GetMemoryLocal (MULTIGRID *theMG, INT size, INT type)
#else

void *GetMemoryForObject (MULTIGRID *theMG, INT size, INT type)
#endif
{
  void **ptr, *obj;
  INT i,j,k,l;
  FORMAT *theFormat;

  if (size == 0)
    return(NULL);
  obj = NULL;

  /* 'ptr' will be set equal to 'theMG->freeObjects[k]' but with	        */
  /* different interpretation: void ** instead of void *. 'ptr'			*/
  /* points to the first two bytes of the object (i.e. unsigned INT ctrl	*/
  /* and INT id) but will be interpreted as a void * pointer, witch points*/
  /* to the next free object.                                                                                   */

  i = (size / ALIGNMENT);
  for (j=0; j<MAXFREEOBJECTS; j++)
  {
    k = (i + j) % MAXFREEOBJECTS;
    l = theMG->SizeOfFreeObjects[k];
    if (l == size)
    {
      if (theMG->freeObjects[k] != NULL)
      {
        ptr = (void **) theMG->freeObjects[k];
        theMG->freeObjects[k] = ptr[0];
        obj = (void *) ptr;
      }
      break;
    }
    if(l == -1)
      break;
  }

  if (obj == NULL)
    obj = GetMem(MGHEAP(theMG),(MEM)size,FROM_BOTTOM);

  if (obj != NULL)
    memset(obj,0,size);

  if (obj == NULL)
    PRINTDEBUG(gm,0,("GetMemoryForObject: cannot allocate %d bytes\n",size));

  return(obj);
}


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
INT PutFreeObject (MULTIGRID *theMG, void *object, INT size, INT type)
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

  return (PutFreeObjectLocal(theMG, object, size, type));
}


INT PutFreeObjectLocal (MULTIGRID *theMG, void *object, INT size, INT type)
#else

INT PutFreeObject (MULTIGRID *theMG, void *object, INT size, INT type)
#endif
{
  void **ptr;
  INT i,j,k,l;
  FORMAT *theFormat;

  memset(object,0,size);
  ((int *)object)[1] = -1;
  ptr = (void **) object;

  /* 'ptr' will be set equal to 'object' but with different inter-		*/
  /* pretation: void ** instead of void *. 'ptr' points to the first		*/
  /* two bytes of the object (i.e. unsigned INT ctrl	and INT id) but         */
  /* will be interpreted as a void * pointer, witch will be set equal   */
  /* to 'theMG->freeObjects[k]' i.e. the first free object.			    */

  i = (size / ALIGNMENT);
  for (j=0; j<MAXFREEOBJECTS; j++)
  {
    k = (i + j) % MAXFREEOBJECTS;
    l = theMG->SizeOfFreeObjects[k];
    if (l == size)
    {
      ptr[0] = theMG->freeObjects[k];
      theMG->freeObjects[k] = object;
      return(0);
    }
    if(l == -1)
    {
      theMG->SizeOfFreeObjects[k] = size;
      ptr[0] = theMG->freeObjects[k];
      theMG->freeObjects[k] = object;
      return(0);
    }
  }

  UserWrite("PutFreeObject: increase MAXFREEOBJECTS\n");

  RETURN(1);
}

/****************************************************************************/
/*D
   CreateVertexSegment - Return pointer to a new vertex segment structure

   SYNOPSIS:
   VSEGMENT *CreateVertexSegment (GRID *theGrid, VERTEX *vertex);

   PARAMETERS:
   .  theGrid - grid of corresponding vertex
   .  vertex - vertex the segment belongs to

   DESCRIPTION:
   This function return a pointer to vertex segment structure.

   RETURN VALUE:
   VSEGMENT *
   .n   pointer to new created VSEGMENT
   .n   NULL if out of memory.
   D*/
/****************************************************************************/

VSEGMENT *CreateVertexSegment (GRID *theGrid, VERTEX *vertex)
{
  VSEGMENT *vs;
  INT i;

  vs = (VSEGMENT *)GetMemoryForObject(MYMG(theGrid),sizeof(struct vsegment),
                                      VSOBJ);
  if (vs==NULL) return (NULL);

  /* initialize data */
  CTRL(vs) = 0;
  PARSETOBJT(vs,VSOBJ);
  VS_PATCH(vs) = NULL;
  for (i=0; i<DIM; i++) LAMBDA(vs,i) = 0.0;

  /* insert in vsegment list */
  NEXTSEG(vs) = VSEG(vertex);
  VSEG(vertex) = vs;

  return(vs);
}


/****************************************************************************/
/*D
   CreateBoundaryVertex - Return pointer to a new boundary vertex structure

   SYNOPSIS:
   VERTEX *CreateBoundaryVertex (GRID *theGrid, VERTEX *after);

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

VERTEX *CreateBoundaryVertex (GRID *theGrid, VERTEX *after)
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
  PARSETOBJT(pv,BVOBJ);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
  TOPNODE(pv) = NULL;
  for (i=0; i<DIM; i++) LCVECT(pv)[i] = 0.0;
  VSEG(pv) = NULL;
  SETONEDGE(pv,0);
  SETMOVE(pv,DIM_OF_BND);
        #ifdef ModelP
  DDD_AttrSet(PARHDRV(pv),theGrid->level);
  DDD_PrioritySet(PARHDRV(pv),PrioVertex);
        #endif

  /* insert in vertex list */
  if (after==NULL)
  {
    SUCCV(pv) = FIRSTVERTEX(theGrid);
    PREDV(pv) = NULL;
    if (SUCCV(pv)!=NULL) PREDV(SUCCV(pv)) = pv;
    else LASTVERTEX(theGrid) = pv;
    FIRSTVERTEX(theGrid) = pv;
  }
  else
  {
    SUCCV(pv) = SUCCV(after);
    PREDV(pv) = after;
    if (SUCCV(pv)!=NULL) PREDV(SUCCV(pv)) = pv;
    else LASTVERTEX(theGrid) = pv;
    SUCCV(after) = pv;
  }

  /* counters */
  theGrid->nVert++;

  return(pv);
}

/****************************************************************************/
/*D
   CreateInnerVertex - Return pointer to a new inner vertex structure

   SYNOPSIS:
   VERTEX *CreateInnerVertex (GRID *theGrid, VERTEX *after);

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

VERTEX *CreateInnerVertex (GRID *theGrid, VERTEX *after)
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
  PARSETOBJT(pv,IVOBJ);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
  TOPNODE(pv) = NULL;
  SETMOVE(pv,DIM);
        #ifdef ModelP
  DDD_AttrSet(PARHDRV(pv),theGrid->level);
  DDD_PrioritySet(PARHDRV(pv),PrioVertex);
        #endif
  for (i=0; i<DIM; i++) LCVECT(pv)[i] = 0.0;

  /* insert in vertex list */
  if (after==NULL)
  {
    SUCCV(pv) = FIRSTVERTEX(theGrid);
    PREDV(pv) = NULL;
    if (SUCCV(pv)!=NULL) PREDV(SUCCV(pv)) = pv;
    else LASTVERTEX(theGrid) = pv;
    FIRSTVERTEX(theGrid) = pv;
  }
  else
  {
    SUCCV(pv) = SUCCV(after);
    PREDV(pv) = after;
    if (SUCCV(pv)!=NULL) PREDV(SUCCV(pv)) = pv;
    else LASTVERTEX(theGrid) = pv;
    SUCCV(after) = pv;
  }

  /* counters */
  theGrid->nVert++;

  return(pv);
}

/****************************************************************************/
/*D
   CreateNode - Return pointer to a new node structure

   SYNOPSIS:
   NODE *CreateNode (GRID *theGrid, NODE *after)

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

NODE *CreateNode (GRID *theGrid, NODE *after)
{
  NODE *pn;
  VECTOR *pv;

  if (TYPE_DEF_IN_GRID(theGrid,NODEVECTOR))
  {
    /* create node and vector */
    pn = GetMemoryForObject(MYMG(theGrid),sizeof(NODE),NDOBJ);
    if (pn==NULL) return(NULL);
    if (CreateVector (theGrid,NULL,NODEVECTOR,&pv))
    {
      DisposeNode (theGrid,pn);
      return (NULL);
    }
    assert (pv != NULL);
    VOBJECT(pv) = (void*)pn;
    NVECTOR(pn) = (void*)pv;
  }
  else
  {
    /* create node */
    pn = GetMemoryForObject(MYMG(theGrid),
                            sizeof(NODE)-sizeof(VECTOR*),NDOBJ);
    if (pn==NULL) return(NULL);
  }

  /* initialize data */
  CTRL(pn) = 0;
  PARSETOBJT(pn,NDOBJ);
  SETCLASS(pn,4);
  SETLEVEL(pn,theGrid->level);
        #ifdef ModelP
  DDD_AttrSet(PARHDR(pn),theGrid->level);
  DDD_PrioritySet(PARHDR(pn),PrioNode);
        #endif
  ID(pn) = (theGrid->mg->nodeIdCounter)++;
  INDEX(pn) = 0;
  START(pn) = NULL;
  NFATHER(pn) = NULL;
  MYVERTEX(pn) = NULL;
  SONNODE(pn) = NULL;
  theGrid->status |= 1;          /* recalculate stiffness matrix */

  /* insert in vertex list */
  if (after==NULL)
  {
    SUCCN(pn) = theGrid->firstNode;
    PREDN(pn) = NULL;
    if (SUCCN(pn)!=NULL) PREDN(SUCCN(pn)) = pn;
    theGrid->firstNode = pn;
    if (theGrid->lastNode==NULL) theGrid->lastNode = pn;
  }
  else
  {
    SUCCN(pn) = SUCCN(after);
    PREDN(pn) = after;
    if (SUCCN(pn)!=NULL) PREDN(SUCCN(pn)) = pn;else theGrid->lastNode = pn;
    SUCCN(after) = pn;
  }

  /* counters */
  theGrid->nNode++;

  return(pn);
}

/****************************************************************************/
/*D
   CreateMidNode - Return pointer to a new node structure on an edge

   SYNOPSIS:
   NODE *CreateMidNode (GRID *theGrid, ELEMENT *theElement, INT edge, NODE *after);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  theElement - pointer to an element
   .  edge - id of an element edge
   .  after - node where the new node should be appended in the node list

   DESCRIPTION:
   This function creates and initializes a new node structure
   at the midpoint of an element edge and returns a pointer to it.

   RETURN VALUE:
   NODE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

NODE *CreateMidNode (GRID *theGrid, ELEMENT *theElement, INT edge, NODE *after)
{
  NODE *theNode;
  VERTEX *theVertex,*v0,*v1;
  VSEGMENT *vs0,*vs1,*vs;
  PATCH *thePatch;
  COORD *global,*local,*lambda,*lambda0,*lambda1,*x[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR bnd_global;
  DOUBLE diff;
  INT n,i,co0,co1,moved;

  co0 = CORNER_OF_EDGE(theElement,edge,0);
  co1 = CORNER_OF_EDGE(theElement,edge,1);
  v0 = MYVERTEX(CORNER(theElement,co0));
  v1 = MYVERTEX(CORNER(theElement,co1));
  theVertex = NULL;

  /* allocate vertex */
  if (OBJT(v0) == BVOBJ && OBJT(v1) == BVOBJ)
  {
    moved = 0;
    /* check if boundary vertex */
    for (vs0=VSEG(v0); vs0!=NULL; vs0=NEXTSEG(vs0))
    {
      thePatch = VS_PATCH(vs0);
      for (vs1=VSEG(v1); vs1!=NULL; vs1=NEXTSEG(vs1))
        if (VS_PATCH(vs1) == thePatch)
        {
          if (theVertex == NULL)
          {
            theVertex = CreateBoundaryVertex(theGrid,NULL);
            if (theVertex == NULL) return(NULL);
            global = CVECT(theVertex);
            local = LCVECT(theVertex);
            V_DIM_LINCOMB(0.5, CVECT(v0), 0.5, CVECT(v1), global);
            V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theElement,co0),
                          0.5, LOCAL_COORD_OF_ELEM(theElement,co1),
                          local);
          }
          if ((vs = CreateVertexSegment(theGrid,theVertex)) == NULL)
          {
            DisposeVertex(theGrid, theVertex);
            return(NULL);
          }
          VS_PATCH(vs) = thePatch;
          lambda = PVECT(vs);
          lambda0 = PVECT(vs0);
          lambda1 = PVECT(vs1);
          for (i=0; i<DIM-1; i++)
            lambda[i] = 0.5 * lambda0[i] + 0.5 * lambda1[i];
          if (Patch_local2global(thePatch,lambda,bnd_global))
            return (NULL);

          /* check if moved */
          V_DIM_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
          if (diff > MAX_PAR_DIST)
          {
            if (moved)
              PrintErrorMessage('W',"CreateMidNode",
                                "inconsistent boundary parametrization");
            else
            {
              SETMOVED(theVertex,1);
              moved = 1;
              CORNER_COORDINATES(theElement,n,x);
            }
            V_DIM_COPY(bnd_global,global);
            GlobalToLocal(n,(const COORD **)x,global,local);
          }
        }
    }
  }

  if (theVertex == NULL)
  {
    /* we need an inner vertex */
    theVertex = CreateInnerVertex(theGrid,NULL);
    if (theVertex==NULL) return(NULL);
    V_DIM_LINCOMB(0.5, CVECT(v0), 0.5, CVECT(v1), CVECT(theVertex));
    V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theElement,co0),
                  0.5, LOCAL_COORD_OF_ELEM(theElement,co1),
                  LCVECT(theVertex));
  }
  VFATHER(theVertex) = theElement;
  SETONEDGE(theVertex,edge);

  /* allocate node */
  theNode = CreateNode(theGrid,after);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  MYVERTEX(theNode) = theVertex;
  NFATHER(theNode) = NULL;
  TOPNODE(theVertex) = theNode;
  SETNTYPE(theNode,MID_NODE);

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

#ifdef __THREEDIM__
NODE *CreateSideNode (GRID *theGrid, ELEMENT *theElement, INT side)
{
  ELEMENTSIDE *theSide;
  VSEGMENT *vs;
  PATCH *thePatch;
  COORD *global,*local,*lambda,*x[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR bnd_global;
  INT n,m,j,k;
  VERTEX *theVertex;
  NODE *theNode;
  DOUBLE fac, diff;

  n = CORNERS_OF_SIDE(theElement,side);
  fac = 1.0 / n;
  theVertex = NULL;

  /* check if boundary vertex */
  if (OBJT(theElement) == BEOBJ)
  {
    theSide = SIDE(theElement,side);
    if (theSide != NULL)
    {
      theVertex = CreateBoundaryVertex(theGrid,NULL);
      if (theVertex == NULL) return(NULL);

      if ((vs = CreateVertexSegment(theGrid,theVertex)) == NULL)
      {
        DisposeVertex(theGrid, theVertex);
        return(NULL);
      }

      global = CVECT(theVertex);
      local = LCVECT(theVertex);
      lambda  = PVECT(vs);
      V_DIM_CLEAR(local);
      V_DIM_CLEAR(global);
      V2_CLEAR(lambda);
      for (j=0; j<n; j++)
      {
        k = CORNER_OF_SIDE(theElement,side,j);
        V_DIM_LINCOMB(1.0,local,1.0,
                      LOCAL_COORD_OF_ELEM(theElement,k),local);
        V_DIM_LINCOMB(1.0,global,1.0,
                      CVECT(MYVERTEX(CORNER(theElement,k))),global);
        V2_LINCOMB(1.0,lambda,1.0,
                   PARAMPTR(theSide,j),lambda);
      }
      V_DIM_SCALE(fac,local);
      V_DIM_SCALE(fac,global);
      V_DIM_SCALE(fac,lambda);

      thePatch = ES_PATCH(theSide);
      VS_PATCH(vs) = thePatch;

      if (Patch_local2global(thePatch,lambda,bnd_global))
        return (NULL);

      /* check if moved */
      V3_EUKLIDNORM_OF_DIFF(bnd_global,global,diff);
      if (diff>MAX_PAR_DIST)
      {
        SETMOVED(theVertex,1);
        V3_COPY(bnd_global,global);
        CORNER_COORDINATES(theElement,m,x);
        GlobalToLocal(m,(const COORD **)x,global,local);
      }
    }
  }

  if (theVertex == NULL)
  {
    theVertex = CreateInnerVertex(theGrid,NULL);
    if (theVertex == NULL) return(NULL);
    global = CVECT(theVertex);
    local = LCVECT(theVertex);
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
  }
  VFATHER(theVertex) = theElement;

  /* create node */
  theNode = CreateNode(theGrid,NULL);
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
  COORD *global,*local;
  INT n,j;
  VERTEX *theVertex;
  NODE *theNode;
  DOUBLE fac;
  COORD *x[MAX_CORNERS_OF_ELEM];

  theVertex = CreateInnerVertex(theGrid,NULL);
  if (theVertex==NULL)
    return(NULL);
  theNode = CreateNode(theGrid,NULL);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  theGrid->status |= 1;
  CORNER_COORDINATES(theElement,n,x);
  global = CVECT(theVertex);
  local = LCVECT(theVertex);
  V_DIM_CLEAR(local);
  for (j=0; j<n; j++)
    V_DIM_LINCOMB(1.0,local,1.0,LOCAL_COORD_OF_ELEM(theElement,j),local);
  fac = 1.0 / n;
  V_DIM_SCALE(fac,local);
  LOCAL_TO_GLOBAL(n,x,local,global);
  VFATHER(theVertex) = theElement;
  NFATHER(theNode) = NULL;
  MYVERTEX(theNode) = theVertex;
  TOPNODE(theVertex) = theNode;
  SETNTYPE(theNode,CENTER_NODE);

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
   static EDGE *CreateEdge (GRID *theGrid, NODE *from, NODE *to);

   PARAMETERS:
   .  theGrid - grid where vertex should be inserted
   .  from - starting node of new edge
   .  to - end node of new edge

   DESCRIPTION:
   This function returns a pointer to a new edge structure.

   RETURN VALUE:
   EDGE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

static EDGE *CreateEdge (GRID *theGrid, NODE *from, NODE *to)
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
  CTRL(pe) = 0;
  CTRL(link1) = 0;
  PARSETOBJT(pe,EDOBJ);
  SETLOFFSET(link0,0);
  SETLOFFSET(link1,1);
  PARSETLEVEL(pe,theGrid->level);
  NBNODE(link0) = to;
  NBNODE(link1) = from;
  SET_NO_OF_ELEM(pe,1);
  PARSETTAG(pe,0);
  SETEDGENEW(pe,1);
  MIDNODE(pe) = NULL;

  /* create vector if */
  if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR))
  {
    if (CreateVector (theGrid,NULL,EDGEVECTOR,&pv))
    {
      DisposeEdge (theGrid,pe);
      return (NULL);
    }
    assert (pv != NULL);
    VOBJECT(pv) = (void*)pe;
    EDVECTOR(pe) = (void*)pv;
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
   NODE **nodes, ELEMENT *after);

   PARAMETERS:
   .  theGrid - grid structure to extend
   .  tag - the element type
   .  objtype - inner element (IEOBJ) or boundary element (BEOBJ)
   .  nodes - list of corner nodes in reference numbering
   .  after - insert after that element (NULL if b.o.l.)

   DESCRIPTION:
   This function creates and initializes a new element and returns a pointer to it.

   RETURN VALUE:
   ELEMENT *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

ELEMENT *CreateElement (GRID *theGrid, INT tag, INT objtype,
                        NODE **nodes, ELEMENT *after)
{
  ELEMENT *pe;
  INT i,j;
  VECTOR *pv;

  if (objtype == IEOBJ)
    pe = GetMemoryForObject(MYMG(theGrid),INNER_SIZE(tag),
                            MAPPED_INNER_OBJT(tag));
  else if (objtype == BEOBJ)
    pe = GetMemoryForObject(MYMG(theGrid),BND_SIZE(tag),
                            MAPPED_BND_OBJT(tag));

  if (pe==NULL) return(NULL);

  /* initialize data */
  PARSETOBJT(pe,objtype);
  PARSETTAG(pe,tag);
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
                   nodes[CORNER_OF_EDGE(pe,i,1)]) == NULL)
    {
      DisposeElement(theGrid,pe);
      return(NULL);
    }

  /* create element vector if */
  if (TYPE_DEF_IN_GRID(theGrid,ELEMVECTOR))
  {
    if (CreateVector (theGrid,NULL,ELEMVECTOR,&pv))
    {
      DisposeElement (theGrid,pe);
      return (NULL);
    }
    assert (pv != NULL);
    VOBJECT(pv) = (void*)pe;
    SET_EVECTOR(pe,(void*)pv);
  }

  /* create side vectors if */
        #ifdef __THREEDIM__
  if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
    for (i=0; i<SIDES_OF_ELEM(pe); i++)
    {
      if (CreateVector (theGrid,NULL,SIDEVECTOR,&pv))
      {
        DisposeElement (theGrid,pe);
        return (NULL);
      }
      assert (pv != NULL);
      VOBJECT(pv) = (void*)pe;
      SET_SVECTOR(pe,i,(void*)pv);
      SETVECTORSIDE(pv,i);
      SETVCOUNT(pv,1);
    }
        #endif

  /* insert in element list */
  if (after==NULL)
  {
    SUCCE(pe) = theGrid->elements;
    PREDE(pe) = NULL;
    if (SUCCE(pe)!=NULL)
      PREDE(SUCCE(pe)) = pe;
    else
      theGrid->lastelement = pe;
    theGrid->elements = pe;
  }
  else
  {
    SUCCE(pe) = SUCCE(after);
    PREDE(pe) = after;
    if (SUCCE(pe)!=NULL)
      PREDE(SUCCE(pe)) = pe;
    else
      theGrid->lastelement = pe;
    SUCCE(after) = pe;
  }

  /* counters */
  theGrid->nElem++;

  /* return ok */
  return(pe);
}

/****************************************************************************/
/*D
   CreateElementSide - Return pointer to a new element side structure

   SYNOPSIS:
   ELEMENTSIDE *CreateElementSide (GRID *theGrid);

   PARAMETERS:
   .  theGrid - grid for which to create

   DESCRIPTION:
   This function creates and initializes a new element side structure and
   returns a pointer to it.

   RETURN VALUE:
   ELEMENTSIDE *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

ELEMENTSIDE *CreateElementSide (GRID *theGrid)
{
  ELEMENTSIDE *ps;

  ps = GetMemoryForObject(MYMG(theGrid),sizeof(ELEMENTSIDE),ESOBJ);
  if (ps==NULL) return(NULL);

  /* initialize data */
  CTRL(ps) = 0;
  PARSETOBJT(ps,ESOBJ);
  ES_PATCH(ps) = NULL;

  /* insert in side list */
  SUCCS(ps) = FIRSTELEMSIDE(theGrid);
  PREDS(ps) = NULL;
  if (SUCCS(ps)!=NULL) PREDS(SUCCS(ps)) = ps;
  FIRSTELEMSIDE(theGrid) = ps;

  /* counters */
  theGrid->nSide++;

  /* return ok */
  return(ps);
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
  ELEMENTSIDE *oldSide,*newSide;
  VSEGMENT *vs;
  PATCH *thePatch;
  COORD *lambda;
  INT i,j;

  IFDEBUG(gm,0)
  assert (OBJT(theElement) == BEOBJ);
  assert (SIDE(theElement,side) != NULL);
  ENDDEBUG

    oldSide = SIDE(theElement,side);
  newSide = CreateElementSide(theGrid);
  if (newSide == NULL)
    RETURN(GM_ERROR);
  SET_SIDE(theSon,son_side,newSide);
  thePatch = ES_PATCH(oldSide);
  ES_PATCH(newSide) = thePatch;

  for (i=0; i<CORNERS_OF_SIDE(theSon,son_side); i++)
  {
    for(vs=VSEG(MYVERTEX(CORNER(theSon,CORNER_OF_SIDE(theSon,son_side,i)))); vs!=NULL; vs = NEXTSEG(vs) )
      if (VS_PATCH(vs) == thePatch)
        break;
    assert(vs!=NULL);
    lambda = PARAMPTR(newSide,i);
    for (j=0; j<DIM-1; j++)
      lambda[j] = LAMBDA(vs,j);
  }

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
  int l;

  if (theMG->topLevel+1>=MAXLEVEL) return(NULL);

  l = theMG->topLevel+1;

  /* allocate grid object */
  theGrid = GetMemoryForObject(theMG,sizeof(GRID),GROBJ);
  if (theGrid==NULL) return(NULL);

  /* fill in data */
  CTRL(theGrid) = 0;
  PARSETOBJT(theGrid,GROBJ);
  theGrid->level = l;
  theGrid->nVert = 0;
  theGrid->nNode = 0;
  theGrid->nEdge = 0;
  theGrid->nElem = 0;
  theGrid->nSide = 0;
  theGrid->nVector = 0;
  theGrid->nCon = 0;

#ifdef __INTERPOLATION_MATRIX__
  theGrid->nIMat = 0;
#endif

  theGrid->status       = 0;
  theGrid->elements = NULL;
  theGrid->lastelement = NULL;
  theGrid->vertices = NULL;
  theGrid->sides = NULL;
  theGrid->firstNode = theGrid->lastNode = NULL;
  theGrid->firstVector = theGrid->lastVector = NULL;
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
   This function creates a multigrid environment item.

   RETURN VALUE:
   MULTIGRID *
   .n   pointer to new MULTIGRID
   .n   NULL if error occured
   D*/
/****************************************************************************/

MULTIGRID *MakeMGItem (const char *name)
{
  if (ChangeEnvDir("/Documents") == NULL) return (NULL);
  if (strlen(name)>=NAMESIZE || strlen(name)<=1) return (NULL);

  return ((MULTIGRID *) MakeEnvItem(name,theMGVarID,sizeof(MULTIGRID)));
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
  return ((MULTIGRID *) SearchEnv(name,"/Documents",theMGVarID,theDocDirID));
}

/****************************************************************************/
/*D
   GetFirstMultigrid - Return a pointer to the first multigrid

   SYNOPSIS:
   MULTIGRID *GetFirstMultigrid ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function returns a pointer to the first multigrid in the /Documents
   directory.

   RETURN VALUE:
   MULTIGRID *
   .n   pointer to MULTIGRID
   .n   NULL if not found.
   D*/
/****************************************************************************/

MULTIGRID *GetFirstMultigrid ()
{
  ENVDIR *theDocDir;
  MULTIGRID *theMG;

  theDocDir = ChangeEnvDir("/Documents");

  assert (theDocDir!=NULL);

  theMG = (MULTIGRID *) ENVDIR_DOWN(theDocDir);

  if (theMG != NULL)
    if (InitElementTypes(theMG)!=GM_OK)
    {
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
   This function returns a pointer to the next multigrid in the /Documents
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
    if (InitElementTypes(MG)!=GM_OK)
    {
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

MULTIGRID *CreateMultiGrid (char *MultigridName, char *BndValProblem, char *format, unsigned long heapSize)
{
  HEAP *theHeap,*theUserHeap;
  MULTIGRID *theMG;
  GRID *theGrid;
  VERTEX **pv;
  NODE *pn;
  INT i,j,k,n,ds,l,FatalError;
  COORD cvect[DIM];
  DOUBLE pardist;
  VSEGMENT *vs;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  PATCH *thePatch, **PatchList;
  PATCH_DESC thePatchDesc;
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
    PRINTDEBUG(gm,0,("CreateMultiGrid: cannot allocate %d bytes\n",
                     heapSize));
    DisposeMultiGrid(theMG);
    return(NULL);
  }

  theBVP = GetBVP(BndValProblem);
  if (theBVP==NULL)
  {
    PrintErrorMessage('E',"CreateMultiGrid","BVP not found");
    return(NULL);
  }
  if (BVP_GetBVPDesc(theBVP,&theBVPDesc))
  {
    PrintErrorMessage('E',"CreateMultiGrid","BVP not evaluated");
    return(NULL);
  }

  /* create a ID-orientated list of patches */
  Mark(theHeap,FROM_TOP);
  if ((PatchList=(PATCH**)GetMem(theHeap,
                                 BVPD_NPATCHES(theBVPDesc)*sizeof(PATCH*),
                                 FROM_TOP))==NULL)
  {
    Release(theHeap,FROM_TOP);
    UserWrite("ERROR: could not allocate memory for PatchList\n");
    return (NULL);
  }
  for (i=0; i<BVPD_NPATCHES(theBVPDesc); i++)
    PatchList[i] = NULL;
  for (thePatch=BVP_GetFirstPatch(theBVP); thePatch!=NULL;
       thePatch=BVP_GetNextPatch(theBVP,thePatch))
  {
    if (Patch_GetPatchDesc(thePatch,&thePatchDesc))
    {
      Release(theHeap,FROM_TOP);
      return (NULL);
    }
    i = PATCH_ID(thePatchDesc);
    if (i<0 || i>=BVPD_NPATCHES(theBVPDesc) || PatchList[i]!=NULL)
    {
      Release(theHeap,FROM_TOP);
      return (NULL);
    }
    PatchList[i] = thePatch;
  }

  /* allocate user data from that heap */

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
    /* clearing this heap provides the possibility of checking the initialization */
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
  theMG->theHeap = theHeap;
  SELECTIONSIZE(theMG) = 0;
  for (i=0; i<MAXLEVEL; i++) theMG->grids[i] = NULL;
  for (i=0; i<MAXFREEOBJECTS; i++)
  {
    theMG->freeObjects[i] = NULL;
    theMG->SizeOfFreeObjects[i] = -1;;
  }

  /* allocate level 0 grid */
  theGrid = CreateNewLevel(theMG);
  if (theGrid==NULL)
  {
    DisposeMultiGrid(theMG);
    return(NULL);
  }

  /* allocate corner vertices pointers */
  n = BVPD_NCORNERS(theBVPDesc);
  theMG->numOfCorners = n;
  pv = (VERTEX **) GetMem(theHeap,n*sizeof(VERTEX *),FROM_BOTTOM);
  if (pv==NULL) { DisposeMultiGrid(theMG); return(NULL); }
  theMG->corners = pv;

        #ifdef ModelP
  memset(pv,0,n*sizeof(VERTEX *));
  if (me!=master)
  {
    Release(theHeap,FROM_TOP);
    return(theMG);
  }
        #endif

  /* create nodes and vertices */
  for (i=0; i<n; i++)
  {
    pv[i] = CreateBoundaryVertex(theGrid,NULL);
    if (pv[i]==NULL) { DisposeMultiGrid(theMG); return(NULL); }
    SETMOVE(pv[i],0);
    SETUSED(pv[i],0);
    pn = CreateNode(theGrid,NULL);
    if (pn==NULL) { DisposeMultiGrid(theMG); return(NULL); }
    MYVERTEX(pn) = pv[i];
    TOPNODE(pv[i]) = pn;
  }

        #ifdef ModelP
  /* TODO: create the corner nodes or delete this */
  if (0) {
    /* identify boundary nodes and vertices on different procs */
    DDD_IdentifyBegin();
    for (i=0; i<n; i++)
    {
      int proc;
      for (proc=0; proc<procs; proc++)
        if (proc!=me)
        {
          NODE *node = TOPNODE(pv[i]);

          /* take number of vertex i as DDD identifier (id=i*3),
             (id+0 is vertex, id+1 is node, id+2 is vector) */

          DDD_IdentifyNumber(PARHDRV(pv[i]), proc, i*3);
          DDD_IdentifyNumber(PARHDR(node), proc, i*3+1);

          /* identify vectors of nodes */
          if (TYPE_DEF_IN_GRID(theGrid,NODEVECTOR))
          {
            DDD_IdentifyNumber(PARHDR(NVECTOR(node)), proc, i+3+2);
          }
        }
    }
    DDD_IdentifyEnd();

    if (me!=master)
    {
      Release(theHeap,FROM_TOP);
      return(theMG);
    }
  }
        #endif

  /* create and fill segment data of vertex */
  for (i=0; i<BVPD_NPATCHES(theBVPDesc); i++)
  {
    thePatch = PatchList[i];
    if (Patch_GetPatchDesc(thePatch,&thePatchDesc))
    {
      Release(theHeap,FROM_TOP);
      return (NULL);
    }

    for( k=0; k<PATCH_N(thePatchDesc); k++ )
    {
      j = PATCH_CID(thePatchDesc,k);
      SETUSED(pv[j],1);
      vs = CreateVertexSegment(theGrid,pv[j]);
      if (vs==NULL) { DisposeMultiGrid(theMG); return(NULL); }

      VS_PATCH(vs) = thePatch;

      FatalError = 0;
      if (DIM_OF_BND == 1)
        switch(k)
        {
        case 0 : LAMBDA(vs,0) = PATCH_LCVECT(thePatchDesc,0)[0];        break;
        case 1 : LAMBDA(vs,0) = PATCH_LCVECT(thePatchDesc,1)[0];        break;
        default : FatalError = 1;                                                break;
        }
      else if (DIM_OF_BND == 2)
        switch(k)
        {
        case 0 : LAMBDA(vs,0) = PATCH_LCVECT(thePatchDesc,0)[0];
          LAMBDA(vs,1) = PATCH_LCVECT(thePatchDesc,0)[1];        break;
        case 1 : LAMBDA(vs,0) = PATCH_LCVECT(thePatchDesc,1)[0];
          LAMBDA(vs,1) = PATCH_LCVECT(thePatchDesc,1)[1];        break;
        case 2 : LAMBDA(vs,0) = PATCH_LCVECT(thePatchDesc,2)[0];
          LAMBDA(vs,1) = PATCH_LCVECT(thePatchDesc,2)[1];        break;
        case 3 : LAMBDA(vs,0) = PATCH_LCVECT(thePatchDesc,3)[0];
          LAMBDA(vs,1) = PATCH_LCVECT(thePatchDesc,3)[1];        break;
        default : FatalError = 1;                                                                        break;
        }
      else
        FatalError = 1;
      if (FatalError)
      {
        PrintErrorMessage('E',"CreateMultiGrid","fatal error: DIM_OF_BND is wrong");
        DisposeMultiGrid(theMG);
        return(NULL);
      }

      /* if the vs is the first of that vertex: set vertex coordinates */
      if( NEXTSEG(vs) == NULL )
      {
        Patch_local2global(thePatch,PVECT(vs),CVECT(pv[j]));
        continue;                          /* next point (k) */
      }

      /* if the vs is not the first one compare the new geometrical data with the old ones */
      Patch_local2global(thePatch,PVECT(vs),cvect);
      for( l = 0, pardist = 0.0  ; l < DIM ; l++)
        pardist += ((pv[j]->iv.x[l]) - cvect[l])*((pv[j]->iv.x[l]) - cvect[l]);
      if ( pardist > (MAX_PAR_DIST*MAX_PAR_DIST) )
      {
        PrintErrorMessage('E',"CreateMultiGrid","multiple vertex segments incompatible");
        DisposeMultiGrid(theMG);
        return(NULL);
      }
    }
  }

  /* check if all corners are defined */
  for (i=0; i<n; i++)
  {
    if (!USED(pv[i]))
    {
      PrintErrorMessage('E',"CreateMultiGrid","not all corners are defined");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
  }
  for (i=0; i<n; i++) SETUSED(pv[i],0);

  /*release heap */
  Release(theHeap,FROM_TOP);

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
  LINK *link0,*link1,*pl;
  EDGE *pe;
  VERTEX *theVertex;
  NODE *to,*father;

  HEAPFAULT(theNode);

  /* call DisposeElement first! */
  assert(START(theNode) == NULL);
  assert(SONNODE(theNode) == NULL);

  /* remove node from node list */
  if (PREDN(theNode)!=NULL)
    SUCCN(PREDN(theNode)) = SUCCN(theNode);
  else
    theGrid->firstNode = SUCCN(theNode);
  if (SUCCN(theNode)!=NULL)
    PREDN(SUCCN(theNode)) = PREDN(theNode);
  else
    theGrid->lastNode = PREDN(theNode);
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
  if (TYPE_DEF_IN_GRID(theGrid,NODEVECTOR))
  {
    if (DisposeVector (theGrid,NVECTOR(theNode)))
      RETURN(1);
    PutFreeObject(theGrid->mg,theNode,sizeof(NODE),NDOBJ);
  }
  else
    PutFreeObject(theGrid->mg,theNode,sizeof(NODE)-sizeof(VECTOR*),NDOBJ);

  /* return ok */
  (theGrid->nNode)--;
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
  VSEGMENT *theVSeg;

  HEAPFAULT(theVertex);

  /* remove vertex from vertex list */
  if (PREDV(theVertex)!=NULL)
    SUCCV(PREDV(theVertex)) = SUCCV(theVertex);
  else
    FIRSTVERTEX(theGrid) = SUCCV(theVertex);

  if (SUCCV(theVertex)!=NULL)
    PREDV(SUCCV(theVertex)) = PREDV(theVertex);
  else
    LASTVERTEX(theGrid) = PREDV(theVertex);

  if( OBJT(theVertex) == BVOBJ )
  {
    theVSeg = VSEG(theVertex);
    while( theVSeg != NULL )
    {
      PutFreeObject(theGrid->mg,theVSeg,sizeof(struct vsegment),VSOBJ);
      theVSeg = NEXTSEG(theVSeg);
    }
    PutFreeObject(theGrid->mg,theVertex,sizeof(struct bvertex),BVOBJ);
  }
  else
    PutFreeObject(theGrid->mg,theVertex,sizeof(struct ivertex),IVOBJ);

  theGrid->nVert--;
  return(0);
}

/****************************************************************************/
/*D
   DisposeElement - Remove element from the data structure

   SYNOPSIS:
   INT DisposeElement (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  theElement - element to remove

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

INT DisposeElement (GRID *theGrid, ELEMENT *theElement)
{
  INT i,j,tag;
  VECTOR *theVector;
  NODE *theNode;
  EDGE *theEdge;
  ELEMENTSIDE *theElementSide;
  ELEMENT *theNeighbor;

  HEAPFAULT(theElement);

  /* remove element from element list */
  if (PREDE(theElement)!=NULL)
    SUCCE(PREDE(theElement)) = SUCCE(theElement);
  else
    theGrid->elements = SUCCE(theElement);
  if (SUCCE(theElement)!=NULL)
    PREDE(SUCCE(theElement)) = PREDE(theElement);
  else
    theGrid->lastelement = PREDE(theElement);

  /* remove element sides if it's a boundary element */
  if (OBJT(theElement)==BEOBJ)
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      theElementSide = SIDE(theElement,i);
      if (theElementSide != NULL)
      {
        if (PREDS(theElementSide)!=NULL)
          SUCCS(PREDS(theElementSide)) = SUCCS(theElementSide);
        else
          theGrid->sides = SUCCS(theElementSide);
        if (SUCCS(theElementSide)!=NULL)
          PREDS(SUCCS(theElementSide)) = PREDS(theElementSide);

        PutFreeObject(theGrid->mg,
                      theElementSide,sizeof(ELEMENTSIDE),ESOBJ);
        theGrid->nSide--;
      }
    }

  for (j=0; j<EDGES_OF_ELEM(theElement); j++)
  {
    theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),
                    CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
    if (NO_OF_ELEM(theEdge)<1)
      RETURN(GM_ERROR);
    if (NO_OF_ELEM(theEdge)==1)
      DisposeEdge(theGrid,theEdge);
    else
      DEC_NO_OF_ELEM(theEdge);
  }

  for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
  {
    theNode = CORNER(theElement,j);
    if (START(theNode) == NULL)
      DisposeNode(theGrid,theNode);
  }

  /* dispose matrices from element-vector */
  if (DisposeConnectionFromElement(theGrid,theElement))
    RETURN(1);

  /* reset neighbor pointers referencing element and dispose vectors in sides if */
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

  /* dispose element */
  /* give it a new tag ! (I know this is somewhat ugly) */
  tag = TAG(theElement);
  if (OBJT(theElement)==BEOBJ)
  {
    PARSETOBJT(theElement,MAPPED_BND_OBJT(tag));
    PutFreeObject(theGrid->mg,theElement,
                  BND_SIZE(tag),MAPPED_BND_OBJT(tag));
  }
  else
  {
    PARSETOBJT(theElement,MAPPED_INNER_OBJT(tag));
    PutFreeObject(theGrid->mg,theElement,INNER_SIZE(tag),
                  MAPPED_INNER_OBJT(tag));
  }

  theGrid->nElem--;
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

  PutFreeObject(theGrid->mg,theGrid,sizeof(GRID),GROBJ);

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

  /* first unlock the mg */
  ((ENVITEM*) theMG)->v.locked = FALSE;

  /* delete mg */
  if (ChangeEnvDir("/Documents")==NULL) RETURN (GM_ERROR);
  if (RemoveEnvItem ((ENVITEM *)theMG)) RETURN (GM_ERROR);

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
    for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
      ID(theNode) = nn++;

    /* elements */
    for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
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
  COORD diff[DIM];

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
   OrderNodesInGrid - reorder double linked ÔNODEÔ list

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
  if (BVP_GetBVPDesc(theBVP,&theBVPDesc)) RETURN (1);

  /* calculate the diameter of the bounding rectangle of the domain */
  InvMeshSize = POW2(GLEVEL(theGrid)) * pow(NN(GRID_ON_LEVEL(theMG,0)),1.0/DIM) / BVPD_RADIUS(theBVPDesc);

  /* allocate memory for the node list */
  theHeap = MGHEAP(theMG);
  Mark(theHeap,FROM_TOP);
  if ((table=GetMem(theHeap,entries*sizeof(NODE *),FROM_TOP))==NULL)
  {
    Release(theHeap,FROM_TOP);
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

  theGrid->firstNode = table[0];
  theGrid->lastNode  = table[entries-1];


  Release(theHeap,FROM_TOP);

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
   PutAtStartOfList - reorder a given set of elements and put them first in the list

   SYNOPSIS:
   INT PutAtStartOfList (GRID *theGrid, INT cnt, ELEMENT **elemList);

   PARAMETERS:
   .  theGrid - elements are part of that level (not checked)
   .  cnt - number of elements in list
   .  elemList - list of elements to reorder

   DESCRIPTION:
   This function reorders a given set of elements and put them first in the list.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   D*/
/****************************************************************************/

INT PutAtStartOfList (GRID *theGrid, INT cnt, ELEMENT **elemList)
{
  ELEMENT *theElement;
  INT i;

  /* remove all elements from list */
  for (i=0; i<cnt; i++)
  {
    theElement = elemList[i];
    if (PREDE(theElement)!=NULL)
      SUCCE(PREDE(theElement)) = SUCCE(theElement);
    else
      theGrid->elements = SUCCE(theElement);
    if (SUCCE(theElement)!=NULL)
      PREDE(SUCCE(theElement)) = PREDE(theElement);
    else
      theGrid->lastelement = PREDE(theElement);
  }

  /* reorder elements locally */
  PREDE(elemList[0]) = NULL;
  SUCCE(elemList[cnt-1]) = NULL;
  for (i=0; i<cnt-1; i++)
  {
    SUCCE(elemList[i]) = elemList[i+1];
    PREDE(elemList[i+1]) = elemList[i];
  }

  /* and insert them at the start of the element list */
  if (FIRSTELEMENT(theGrid)!=NULL) PREDE(FIRSTELEMENT(theGrid)) = elemList[cnt-1];
  SUCCE(elemList[cnt-1]) = FIRSTELEMENT(theGrid);
  FIRSTELEMENT(theGrid) = elemList[0];
  if (LASTELEMENT(theGrid)==NULL) LASTELEMENT(theGrid)=elemList[cnt-1];

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
   INT InsertInnerNode (MULTIGRID *theMG, COORD *pos);

   PARAMETERS:
   .  theMG - multigrid structure
   .  pos - array containing position

   DESCRIPTION:
   This function inserts a inner node into level 0.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT InsertInnerNode (MULTIGRID *theMG, COORD *pos)
{
  GRID *theGrid;
  VERTEX *theVertex;
  NODE *theNode;
  INT i;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertInnerNode","only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* create objects */
  theVertex = CreateInnerVertex(theGrid,NULL);
  if (theVertex==NULL)
  {
    PrintErrorMessage('E',"InsertInnerNode","cannot create vertex");
    RETURN(GM_ERROR);
  }
  theNode = CreateNode(theGrid,NULL);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertInnerNode","cannot create node");
    RETURN(GM_ERROR);
  }

  /* fill data */
  for (i=0; i<DIM; i++) CVECT(theVertex)[i] = pos[i];
  SETMOVE(theVertex,DIM);

  INDEX(theNode) = 0;
  MYVERTEX(theNode) = theVertex;

  return(GM_OK);
}

/****************************************************************************/
/*D
   InsertBoundaryNode - Insert a boundary node

   SYNOPSIS:
   INT InsertBoundaryNode (MULTIGRID *theMG, INT bnd_seg_id, COORD *pos);

   PARAMETERS:
   .  theMG - multigrid structure
   .  bnd_seg_id - number of boundary segment
   .  pos - array containing position

   DESCRIPTION:
   This function inserts a boundary node into level 0.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT InsertBoundaryNodeFromPatch (MULTIGRID *theMG, PATCH *thePatch, COORD *pos)
{
  GRID *theGrid;
  NODE *theNode;
  VERTEX *theVertex;
  VSEGMENT *vs,*vs1;
  PATCH *Patch;
  PATCH_DESC thePatchDesc;
  INT i,j,k,from,to,npc,node_on_edge;
  COORD *from_local,*to_local;
  COORD diff1[DIM_OF_BND],diff2[DIM_OF_BND],local[DIM_OF_BND];
  COORD_VECTOR global;
  DOUBLE val,lambda;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertBoundaryNodeFromPatch","only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  if (thePatch==NULL)
  {
    PrintErrorMessage('E',"InsertBoundaryNodeFromPatch","segment not found");
    RETURN(GM_ERROR);
  }
  if(Patch_GetPatchDesc(thePatch,&thePatchDesc)) RETURN(GM_ERROR);

  /* check distance from corner */
  if (DIM==2)
  {
    if ( (pos[0]<PATCH_LCVECT(thePatchDesc,0)[0])||(pos[0]> PATCH_LCVECT(thePatchDesc,1)[0]) )
    {
      PrintErrorMessage('E',"InsertBoundaryNodeFromPatch","parameter not in range of segment");
      RETURN(GM_ERROR);
    }
    if (  (fabs(pos[0]-PATCH_LCVECT(thePatchDesc,0)[0])<SMALL_C) || (fabs(pos[0]-PATCH_LCVECT(thePatchDesc,1)[0])<SMALL_C) )
    {
      PrintErrorMessage('E',"InsertBoundaryNodeFromPatch","parameter describes one of the corners of the segment");
      RETURN(GM_ERROR);
    }
  }
  if (DIM==3)
  {
    /* TODO: check range of parameters
       for(i=0; i<DIM_OF_BND; i++)
          if ((pos[0]<PATCH_LCVECT(thePatchDesc,0)[i])||(pos[0]> PATCH_LCVECT(thePatchDesc,1)[i]) )
                {
                  PrintErrorMessage('E',"InsertBoundaryNodeFromPatch",
                                                        "parameter not in range of segment");
                  RETURN(GM_ERROR);
                }*/

    i = 0;
    if (fabs(pos[0]-PATCH_LCVECT(thePatchDesc,0)[0])+fabs(pos[1]-PATCH_LCVECT(thePatchDesc,0)[1])<SMALL_C) i = 1;
    if (fabs(pos[0]-PATCH_LCVECT(thePatchDesc,1)[0])+fabs(pos[1]-PATCH_LCVECT(thePatchDesc,1)[1])<SMALL_C) i = 1;
    if (fabs(pos[0]-PATCH_LCVECT(thePatchDesc,2)[0])+fabs(pos[1]-PATCH_LCVECT(thePatchDesc,2)[1])<SMALL_C) i = 1;
    if (fabs(pos[0]-PATCH_LCVECT(thePatchDesc,3)[0])+fabs(pos[1]-PATCH_LCVECT(thePatchDesc,3)[1])<SMALL_C) i = 1;
    if (i)
    {
      PrintErrorMessage('E',"InsertBoundaryNodeFromPatch","parameter describes one of the corners of the segment");
      RETURN(GM_ERROR);
    }
  }

  /* create objects */
  theVertex = CreateBoundaryVertex(theGrid,NULL);
  if (theVertex==NULL)
  {
    PrintErrorMessage('E',"InsertBoundaryNodeFromPatch",
                      "cannot create vertex");
    RETURN(GM_ERROR);
  }
  Patch_local2global(thePatch,pos,CVECT(theVertex));
  SETMOVE(theVertex,DIM_OF_BND);

  vs = CreateVertexSegment(theGrid,theVertex);
  if (vs==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertBoundaryNodeFromPatch",
                      "cannot create vertexsegment");
    RETURN(GM_ERROR);
  }
  /* fill data into the first vertexsegment */
  for(i=0; i<DIM_OF_BND; i++)
    LAMBDA(vs,i) =        pos[i];
  VS_PATCH(vs) = thePatch;

  theNode = CreateNode(theGrid,NULL);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertBoundaryNodeFromPatch",
                      "cannot create node");
    RETURN(GM_ERROR);
  }
  MYVERTEX(theNode) = theVertex;

  /* check if the vertex is on an edge of the patch */
        #ifdef __THREEDIM__
  node_on_edge = FALSE;
  npc = PATCH_N(thePatchDesc);
  from = PATCH_CID(thePatchDesc,npc-1);
  from_local =  PATCH_LCVECT(thePatchDesc,npc-1);
  for (k=0; k<npc; k++ )
  {
    to = PATCH_CID(thePatchDesc,k);
    to_local = PATCH_LCVECT(thePatchDesc,k);
    V2_SUBTRACT (from_local,to_local,diff1);
    V2_SUBTRACT (from_local,pos,diff2);
    V2_VECTOR_PRODUCT(diff1,diff2,val);
    if (ABS(val) < SMALL_C * 1000)
    {
      node_on_edge = TRUE;
      break;
    }
    from = to;
    from_local =  to_local;
  }

  /* fill vertex segment list */
  if (node_on_edge)
  {
    V2_EUKLIDNORM(diff1,val);
    V2_EUKLIDNORM(diff2,lambda);
    lambda /= val;
    for (vs1 = VSEG(theMG->corners[from]); vs1 != NULL; vs1 = NEXTSEG(vs1))
    {
      Patch = VS_PATCH(vs1);
      if (Patch == thePatch)
        continue;
      if(Patch_GetPatchDesc(Patch,&thePatchDesc)) RETURN(GM_ERROR);
      npc = PATCH_N(thePatchDesc);
      for (i=0; i < npc; i++)
        if (to == PATCH_CID(thePatchDesc,i))
        {
          for (k=0; k < npc; k++)
            if (from == PATCH_CID(thePatchDesc,k))
            {
              vs = CreateVertexSegment(theGrid,theVertex);
              if (vs==NULL)
              {
                DisposeNode(theGrid, theNode);
                PrintErrorMessage('E',
                                  "InsertBoundaryNodeFromPatch",
                                  "cannot create vertexsegment");
                RETURN(GM_ERROR);
              }
              V2_LINCOMB(lambda,PATCH_LCVECT(thePatchDesc,k),
                         (1.0-lambda),PATCH_LCVECT(thePatchDesc,i),
                         local);
              Patch_local2global(Patch,local,global);
              V3_EUKLIDNORM_OF_DIFF(CVECT(theVertex),global,val);

              /* fill data into the vertex segment */
              if (ABS(val) < SMALL_C * 1000)
              {
                for(i=0; i<DIM_OF_BND; i++)
                  LAMBDA(vs,i) = local[i];
                VS_PATCH(vs) = Patch;
              }
              else
              {
                Patch_global2local(Patch,CVECT(theVertex),local);
                for(i=0; i<DIM_OF_BND; i++)
                  LAMBDA(vs,i) = local[i];
                VS_PATCH(vs) = Patch;
              }
              break;
            }
          break;
        }
    }
  }
        #endif

  /* fill data into node/vertex */
  INDEX(theNode) = 0;

  return(GM_OK);
}

INT InsertBoundaryNode (MULTIGRID *theMG, INT patch_id, COORD *pos)
{
  NODE *theNode;
  VERTEX *theVertex;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;
  int i;

  /* get BVP description */
  theBVP = MG_BVP(theMG);
  if (BVP_GetBVPDesc(theBVP,&theBVPDesc)) RETURN (GM_ERROR);

  /* scan input */
  if ((patch_id<0)||(patch_id>=BVPD_NPATCHES(theBVPDesc)))
  {
    PrintErrorMessage('E',"InsertBoundaryNode","segment id out of range");
    RETURN(GM_ERROR);
  }
  for (thePatch=BVP_GetFirstPatch(theBVP); thePatch!=NULL; thePatch=BVP_GetNextPatch(theBVP,thePatch))
  {
    if(Patch_GetPatchDesc(thePatch,&thePatchDesc)) RETURN(GM_ERROR);
    if (PATCH_ID(thePatchDesc)==patch_id) break;
  }

  return(InsertBoundaryNodeFromPatch(theMG,thePatch,pos));
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
  int i;

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
   static ELEMENT *FindFather(VERTEX *vptr);

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

#ifdef __TWODIM__
static ELEMENT *FindFather(VERTEX *vptr)
{
  short i;
  COORD x[2],lambda,lambda1,lambda2;
  ELEMENT *eptr,*eptr1;
  ELEMENTSIDE *eside;
  VSEGMENT *vseg;

  if ((eptr=VFATHER(vptr))==NULL)
    return(NULL);

  if (OBJT(vptr)==BVOBJ)
  {
    /* for boundary vertices only a consistency check is made */
    eside=SIDE(eptr,ONEDGE(vptr));
    vseg = VSEG(vptr);

    if (VS_PATCH(vseg)!=ES_PATCH(eside))
      return(NULL);

    /* for higher dimensions (3d) the following must be generalized */
    lambda=LAMBDA(vseg,0);
    lambda1=PARAM(eside,0,0);
    lambda2=PARAM(eside,1,0);
    if ( ((lambda1<=lambda)&&(lambda2>=lambda)) || ((lambda2<=lambda)&&(lambda1>=lambda)) )
      return(eptr);
    else
      return(NULL);
  }

  x[0]=XC(vptr); x[1]=YC(vptr);

  eptr=VFATHER(vptr);
  if (PointInElement(x,eptr))
    return(eptr);

  for (i=0; i<SIDES_OF_ELEM(eptr); i++)
    if (PointInElement(x,(eptr1=NBELEM(eptr,i))))
      return(eptr1);

  return(NULL);
}
#endif

/****************************************************************************/
/*D
   Local2Global - Updates the global coordinates of a vertex

   SYNOPSIS:
   static INT Local2Global (MULTIGRID *theMG, VERTEX *vptr);

   PARAMETERS:
   .  theMG - pointer to 'MULTIGRID' structure.
   .  vptr - pointer to a 'VERTEX'

   DESCRIPTION:
   This function updates the global coordinates of a vertex by
   evaluating the local coordinates in the father element of the vertex.
   This function overwrites the current coordinates, it works for interior
   and boundary vertices.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   >0 if an error occured or vertex is on level 0
   D*/
/****************************************************************************/

#ifdef __TWODIM__
static INT Local2Global (MULTIGRID *theMG, VERTEX *vptr)
{
  ELEMENTSIDE *theSide;
  SHORT i,n,side;
  COORD x[DIM],xi,eta;
  COORD lambda1,lambda2,lambdaa,lambdae;
  COORD lambda,z,lm0,lm1,lmopt,smin,dlambda,xc,yc,r[DIM],lambdaopt,s;
  ELEMENT *eptr;
  VERTEX *vptr1,*vptr2,*vptra,*vptre;
  VSEGMENT *vseg;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;

  if ((eptr=VFATHER(vptr))==NULL)
    return(8040);

  n=TAG(eptr);
  if (OBJT(vptr)==BVOBJ)
  {
    vseg = VSEG(vptr);
    if ((thePatch=VS_PATCH(vseg))==NULL) return(8041);
    if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return(8041);

    side=ONEDGE(vptr);
    vptr1=MYVERTEX(CORNER(eptr,side));
    vptr2=MYVERTEX(CORNER(eptr,(side+1)%n));

    vptra=theMG->corners[PATCH_CID(thePatchDesc,0)]; lambdaa=PATCH_LCVECT(thePatchDesc,0)[0];
    vptre=theMG->corners[PATCH_CID(thePatchDesc,1)]; lambdae=PATCH_LCVECT(thePatchDesc,1)[0];

    if (vptr1==vptra)
      lambda1=lambdaa;
    else
    if (vptr1==vptre)
      lambda1=lambdae;
    else
      lambda1=LAMBDA(VSEG(vptr1),0);

    if (vptr2==vptra)
      lambda2=lambdaa;
    else
    if (vptr2==vptre)
      lambda2=lambdae;
    else
      lambda2=LAMBDA(VSEG(vptr2),0);

    theSide = SIDE(eptr,side);
    smin = 1.0E10;
    lm0 = PARAM(theSide,0,0);
    lm1 = PARAM(theSide,1,0);
    xc = 0.5*(XC(vptr1)+XC(vptr2));
    yc = 0.5*(YC(vptr1)+YC(vptr2));
    dlambda = (lm1-lm0)/((COORD) RESOLUTION);
    lambda = lm0;
    for (i=1; i<RESOLUTION; i++)
    {
      lambda += dlambda;
      if (Patch_local2global(thePatch,&lambda,r)) return (NULL);
      s = (r[0]-xc)*(r[0]-xc)+(r[1]-yc)*(r[1]-yc);
      if (s<smin)
      {
        smin = s;
        lambdaopt = lambda;
      }
    }
    z = (lambdaopt-lm0)/(lm1-lm0);
    LAMBDA(vseg,0)=(1-z)*lambda1+z*lambda2;
    if (Patch_local2global(thePatch,PVECT(vseg),CVECT(vptr))) return(8041);
  }
  else
  {
    xi=XI(vptr); eta=ETA(vptr);
    x[0]=x[1]=0;
    for (i=0; i<n; i++)
    {
      vptr1=MYVERTEX(CORNER(eptr,i));
      x[0]+=XC(vptr1)*N(n,i,xi,eta);
      x[1]+=YC(vptr1)*N(n,i,xi,eta);
    }

    XC(vptr)=x[0]; YC(vptr)=x[1];
  }

  return(0);
}
#endif

/****************************************************************************/
/*D
   Global2Local - Updates the local coordinates of a vertex

   SYNOPSIS:
   static INT Global2Local (MULTIGRID *theMG, VERTEX *vptr);

   PARAMETERS:
   .  theMG - pointer to 'MULTIGRID'
   .  vptr - pointer to a 'VERTEX'

   DESCRIPTION:
   This function updates the local coordinates of a vertex in the
   local coordinate system of the father element. It is assumed that
   the vertex is inside the father element. If not an error is returned.

   RETURN VALUE:
   INT
   .n    0 if ok vertex is level 0
   .n    >0 if error occured (invalid trafo, not inside)
   D*/
/***************************************************************************/

#ifdef __TWODIM__
static INT Global2Local (MULTIGRID *theMG, VERTEX *vptr)
{
  short n,side;
  DOUBLE t1x,t2x,t3x,t1y,t2y,t3y,a,b,c,D,xi1,xi2,eta1,eta2;
  DOUBLE x1,x2,y1,y2,x3,y3,x4,y4,x,y;
  DOUBLE a1,a2,a3,a4,b1,b2,b3,b4;
  COORD lambda1,lambda2,lambdaa,lambdae;
  ELEMENT *eptr;
  VERTEX *vptr1,*vptr2,*vptra,*vptre;
  VSEGMENT *vseg;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;

  if ((eptr=VFATHER(vptr))==NULL)
    return(0);

  if (OBJT(vptr)==BVOBJ)
  {
    n=TAG(eptr);
    vseg = VSEG(vptr);
    thePatch=VS_PATCH(vseg);
    if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return (1);
    side=ONEDGE(vptr);
    vptr1=MYVERTEX(CORNER(eptr,side));
    vptr2=MYVERTEX(CORNER(eptr,(side+1)%n));

    vptra=theMG->corners[PATCH_CID(thePatchDesc,0)];
    lambdaa=PATCH_LCVECT(thePatchDesc,0)[0];
    vptre=theMG->corners[PATCH_CID(thePatchDesc,1)];
    lambdae=PATCH_LCVECT(thePatchDesc,1)[0];

    if (vptr1==vptra)
      lambda1=lambdaa;
    else
    if (vptr1==vptre)
      lambda1=lambdae;
    else
      lambda1=LAMBDA(VSEG(vptr1),0);

    if (vptr2==vptra)
      lambda2=lambdaa;
    else
    if (vptr2==vptre)
      lambda2=lambdae;
    else
      lambda2=LAMBDA(VSEG(vptr2),0);

    /* set local boundary coordinate */
    c=(LAMBDA(vseg,0)-lambda1)/(lambda2-lambda1);

    /* set local coordinates */
    switch(n)
    {
    case TRIANGLE :
      XI(vptr) = ETA(vptr) = 0.0;
      switch (side)
      {
      case 0 : XI(vptr)=c; break;
      case 1 : XI(vptr)=1-c; ETA(vptr)=c; break;
      case 2 : ETA(vptr)=1-c; break;
      }
      break;

    case QUADRILATERAL :
      switch (side)
      {
      case 0 : XI(vptr)=2*c-1; ETA(vptr) = -1;  break;
      case 1 : XI(vptr) = 1;   ETA(vptr)=2*c-1; break;
      case 2 : XI(vptr)=1-2*c; ETA(vptr) = 1;   break;
      case 3 : XI(vptr) = -1;  ETA(vptr)=1-2*c; break;
      }
      break;
    }
  }
  else
  {
    /* compute local coordinates */
    vptr1=MYVERTEX(CORNER(eptr,0)); x1 = XC(vptr1); y1 = YC(vptr1);
    vptr1=MYVERTEX(CORNER(eptr,1)); x2 = XC(vptr1); y2 = YC(vptr1);
    vptr1=MYVERTEX(CORNER(eptr,2)); x3 = XC(vptr1); y3 = YC(vptr1);
    x = XC(vptr);
    y = YC(vptr);

    switch (TAG(eptr))
    {
    case TRIANGLE :
      t1x = x2-x1; t2x = x3-x1; t3x = x1;
      t1y = y2-y1; t2y = y3-y1; t3y = y1;
      D = t1x*t2y-t2x*t1y;
      if (D<0.0)
        return(8051);

      XI(vptr) = (t2y*(x-t3x)-t2x*(y-t3y))/D;
      ETA(vptr) = (-t1y*(x-t3x)+t1x*(y-t3y))/D;
      break;

    case QUADRILATERAL :
      vptr1=MYVERTEX(CORNER(eptr,3)); x4 = XC(vptr1); y4 = YC(vptr1);

      a1=-x+.25*(x1+x2+x3+x4); b1=-y+.25*(y1+y2+y3+y4);
      a2=.25*(-x1+x2+x3-x4); b2=.25*(-y1+y2+y3-y4);
      a3=.25*(-x1-x2+x3+x4); b3=.25*(-y1-y2+y3+y4);
      a4=.25*(x1-x2+x3-x4); b4=.25*(y1-y2+y3-y4);

      c=a2*b1-a1*b2;
      b=a4*b1-a3*b2+a2*b3-a1*b4;
      a=a4*b3-a3*b4;

      if (ABS(a)<SMALL_D)
      {
        if (ABS(b)<SMALL_D)
          return(8051);

        eta1 = eta2 = -c/b;
      }
      else
      {
        D=b*b-4*a*c;
        if (D<0) return(8051);

        eta1=(-b+sqrt(D))/(2*a); eta2=(-b-sqrt(D))/(2*a);
      }

      c=a3*b1-a1*b3;
      b=a4*b1+a3*b2-a2*b3-a1*b4;
      a=a4*b2-a2*b4;

      if (ABS(a)<SMALL_D)
      {
        if (ABS(b)<SMALL_D)
          return(8051);

        xi1 = xi2 = -c/b;
      }
      else
      {
        D=b*b-4*a*c;
        if (D<0) return(8051);

        xi1=(-b+sqrt(D))/(2*a); xi2=(-b-sqrt(D))/(2*a);
      }

      if ((xi1>=-1-SMALL1)&&(xi1<=1+SMALL1)&&(eta1>=-1-SMALL1)&&(eta1<=1+SMALL1))
      {
        XI(vptr) = xi1;
        ETA(vptr) = eta1;
        break;
      }

      if ((xi2>=-1-SMALL1)&&(xi2<=1+SMALL1)&&(eta1>=-1-SMALL1)&&(eta1<=1+SMALL1))
      {
        XI(vptr) = xi2;
        ETA(vptr) = eta1;
        break;
      }

      if ((xi1>=-1-SMALL1)&&(xi1<=1+SMALL1)&&(eta2>=-1-SMALL1)&&(eta2<=1+SMALL1))
      {
        XI(vptr) = xi1;
        ETA(vptr) = eta2;
        break;
      }

      if ((xi2>=-1-SMALL1)&&(xi2<=1+SMALL1)&&(eta2>=-1-SMALL1)&&(eta2<=1+SMALL1))
      {
        XI(vptr) = xi2;
        ETA(vptr) = eta2;
        break;
      }

      return(8051);
      break;
    }
  }

  return(0);
}
#endif

/****************************************************************************/
/*D
   MoveInnerNode - Let user enter a new position for an inner node

   SYNOPSIS:
   INT MoveInnerNode (MULTIGRID *theMG, NODE *theNode, COORD *newPos);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  newPos - new position (x,y)

   DESCRIPTION:
   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   `Function only implemented in 2D version !`

   RETURN VALUE:
   INT
   .n   0 when ok
   .n   >0 when error occured.
   D*/
/****************************************************************************/

#ifdef __TWODIM__
INT MoveInnerNode (MULTIGRID *theMG, NODE *theNode, COORD *newPos)
{
  GRID *theGrid2;
  int k,k2;
  NODE *theNode2;
  VERTEX *theVertex,*theVertex2;
  ELEMENT *theElement,*oldElement;
  double x,y,oldx,oldy;

  k = LEVEL(theNode);

  /* set k (and theNode) to the level where the node appears the first time */
  while ((theNode2=NFATHER(theNode))!=NULL)
  {
    theNode=theNode2;
    k=k-1;
  }

  theVertex = MYVERTEX(theNode);
  oldx = XC(theVertex); oldy = YC(theVertex);

  if (OBJT(theVertex)!=IVOBJ)
  {
    PrintErrorMessage('E',"MoveInnerNode","no inner node passed");
    return(GM_ERROR);
  }

  x = newPos[0]; y = newPos[1];

  /* set values */
  XC(theVertex) = x; YC(theVertex) = y;
  if (VFATHER(theVertex)!=NULL)
  {
    oldElement=VFATHER(theVertex);
    if ((theElement=FindFather(theVertex))==NULL)
    {
      PrintErrorMessage('E',"MoveInnerNode","No father element! Probably you tried to move the vertex too far!");
      XC(theVertex) = oldx; YC(theVertex) = oldy;
      return(GM_ERROR);
    }
    else
      VFATHER(theVertex)=theElement;

    if (Global2Local(theMG,theVertex)!=0)
    {
      PrintErrorMessage('E',"MoveInnerNode","Error in Global2Local");
      VFATHER(theVertex)=oldElement;
      XC(theVertex) = oldx; YC(theVertex) = oldy;
      return(GM_ERROR);
    }
  }

  /*	now we correct the global coordinates for all levels above, since it is not
          easy to find exactly the vertices whose global coordinates have changed */

  for(k2=k+1; k2<=theMG->topLevel; k2++)
  {
    theGrid2 = theMG->grids[k2];
    for (theVertex2=FIRSTVERTEX(theGrid2); theVertex2!=NULL; theVertex2=SUCCV(theVertex2))
      if (Local2Global(theMG,theVertex2)!=0)
      {
        PrintErrorMessage('E',"MoveInnerNode","Fatal error in correcting global coordinates for higher levels. Grid may be inconsistent from now on.");
        return(GM_ERROR);
      }
  }

  /* OK, done */
  return(GM_OK);
}
#endif

#ifdef __THREEDIM__
INT MoveInnerNode (MULTIGRID *theMG, NODE *theNode, COORD *newPos)
{
  UserWrite("MoveInnerNode not implemented in 3D version jet\n");
  return (0);
}
#endif



/****************************************************************************/
/*D
   MoveBoundaryNode - Let user enter a new position

   SYNOPSIS:
   INT MoveBoundaryNode (MULTIGRID *theMG, NODE *theNode, INT segid,
   COORD *newPos);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  segid - new boundary segment id
   .  newPos - new parameter lambda

   DESCRIPTION:
   This function moves a boundary node to a new position. The position of
   all nodes on finer grids is updated to reflect these changes.

   `Function only implemented in 2D version!`

   RETURN VALUE:
   INT
   .n    0 when ok
   .n    >0 when error occured.
   D*/
/****************************************************************************/

#ifdef __TWODIM__
INT MoveBoundaryNode (MULTIGRID *theMG, NODE *theNode, INT patchid, COORD *newPos)
{
  GRID *theGrid,*theGrid2;
  int i,n,k,k2;
  NODE *theNode2;
  VERTEX *theVertex,*theVertex2;
  ELEMENT *theElement,*oldElement;
  ELEMENTSIDE *theSide;
  double oldx,oldy,l,oldl;
  BVP             *theBVP;
  BVP_DESC theBVPDesc;
  PATCH *thePatch, *oldPatch;
  PATCH_DESC thePatchDesc;

  /* get BVP description */
  theBVP = MG_BVP(theMG);
  if (BVP_GetBVPDesc(theBVP,&theBVPDesc))
  {
    PrintErrorMessage('E',"MoveBoundaryNode","cannot evaluate BVP");
    return(1);
  }

  k = LEVEL(theNode);
  theGrid = GRID_ON_LEVEL(theMG,k);

  /* set k (and theNode) to the level where the node appears the first time */
  while ((theNode2=NFATHER(theNode))!=NULL)
  {
    theNode=theNode2;
    k=k-1;
  }

  theVertex = MYVERTEX(theNode);
  oldx = XC(theVertex); oldy = YC(theVertex);

  if (OBJT(theVertex)!=BVOBJ)
  {
    PrintErrorMessage('E',"MoveBoundaryNode","no boundary node passed");
    return(GM_ERROR);
  }

  if (MOVE(theVertex)==0)
  {
    PrintErrorMessage('W',"MoveBoundaryNode","corners cannot be moved");
    return(GM_ERROR);
  }

  l = oldl = LAMBDA(VSEG(theVertex),0);
  oldPatch = VS_PATCH(VSEG(theVertex));
  if (START(theNode)==NULL)
  {
    if(patchid >= BVPD_NPATCHES(theBVPDesc))
    {
      PrintErrorMessage('E',"MoveBoundaryNode","patchid out of range");
      return(GM_ERROR);
    }
    thePatch = Patch_GetPatchByID(theBVP,patchid);
  }
  else thePatch = VS_PATCH(VSEG(theVertex));

  if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return (1);
  patchid = PATCH_ID(thePatchDesc);

  l = newPos[0];

  if ((l<PATCH_LCVECT(thePatchDesc,0)[0]) || (l>PATCH_LCVECT(thePatchDesc,1)[0]))
  {
    PrintErrorMessage('E',"MoveBoundaryNode","parameter out of range");
    return(GM_ERROR);
  }

  LAMBDA(VSEG(theVertex),0) = l;
  if (Patch_local2global(thePatch,PVECT(VSEG(theVertex)),CVECT(theVertex))) return (1);
  VS_PATCH(VSEG(theVertex)) = thePatch;

  if (VFATHER(theVertex)!=NULL)
  {
    oldElement=VFATHER(theVertex);
    if ((theElement=FindFather(theVertex))==NULL)
    {
      PrintErrorMessage('E',"MoveBoundaryNode","No father element! Probably you have tried to move the vertex too far");
      XC(theVertex) = oldx; YC(theVertex) = oldy; LAMBDA(VSEG(theVertex),0)=oldl; VS_PATCH(VSEG(theVertex)) = oldPatch;
      return(GM_ERROR);
    }
    else
      VFATHER(theVertex)=theElement;

    if (Global2Local(theMG,theVertex)!=0)
    {
      PrintErrorMessage('E',"MoveBoundaryNode","Error in Global2Local!");
      VFATHER(theVertex)=oldElement;
      XC(theVertex) = oldx; YC(theVertex) = oldy; LAMBDA(VSEG(theVertex),0)=oldl;
      return(GM_ERROR);
    }
  }

  /*	now we correct the global coordinates for all new vertices in the levels above, since it is not
          easy to find exactly the vertices whose global coordinates have changed */
  for(k2=k+1; k2<=theMG->topLevel; k2++)
  {
    theGrid2 = theMG->grids[k2];
    for (theVertex2=FIRSTVERTEX(theGrid2); theVertex2!=NULL; theVertex2=SUCCV(theVertex2))
      if (Local2Global(theMG,theVertex2)!=0)
      {
        PrintErrorMessage('E',"MoveBoundaryNode","Fatal error in correcting global coordinates for higher levels. Grid may be inconsistent from now on");
        return(GM_ERROR);
      }
  }

  /* update element sides on this level */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if (OBJT(theElement)==BEOBJ)
    {
      n = SIDES_OF_ELEM(theElement);
      for (i=0; i<n; i++)
      {
        theSide = SIDE(theElement,i);
        if (theSide!=NULL)
        {
          if (MYVERTEX(CORNER(theElement,i))==theVertex)
            PARAM(theSide,0,0) = l;
          if (MYVERTEX(CORNER(theElement,(i+1)%n))==theVertex)
            PARAM(theSide,1,0) = l;
        }
      }
    }

  /* update (for simplicity) all sides on that segment for all levels above */
  for(k2=k+1; k2<=theMG->topLevel; k2++)
  {
    theGrid2 = theMG->grids[k2];
    for (theElement=theGrid2->elements; theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (OBJT(theElement)==BEOBJ)
      {
        n = SIDES_OF_ELEM(theElement);
        for (i=0; i<n; i++)
        {
          theSide = SIDE(theElement,i);
          if (theSide!=NULL)
          {
            theNode2=CORNER(theElement,i);
            if (NFATHER(theNode2)==NULL)
            {
              theVertex2=MYVERTEX(theNode2);
              if (VS_PATCH(VSEG(theVertex2))==thePatch)
                PARAM(theSide,0,0) = LAMBDA(VSEG(theVertex2),0);
            }
            theNode2=CORNER(theElement,(i+1)%n);
            if (NFATHER(theNode2)==NULL)
            {
              theVertex2=MYVERTEX(theNode2);
              if (VS_PATCH(VSEG(theVertex2))==VS_PATCH(VSEG(theVertex)))
                PARAM(theSide,1,0) = LAMBDA(VSEG(theVertex2),0);
            }
          }
        }
      }
    }
  }

  return(GM_OK);
}
#endif

#ifdef __THREEDIM__
INT MoveBoundaryNode (MULTIGRID *theMG, NODE *theNode, INT segid, COORD *newPos)
{
  UserWrite("MoveBoundaryNode not implemented in 3D version jet\n");
  return (0);
}
#endif

/****************************************************************************/
/*D
   SmoothMultiGrid - Interprete and execute a smooth multigrid command

   SYNOPSIS:
   INT SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  niter - number of iterations to do
   .  bdryFlag - see 'smoothmg' command.

   DESCRIPTION:
   This function smoothes all grid levels of a multigrid hierarchy.
   It processes all grid levels from bottom to top. Within each grid
   level all nodes that are not already in coarser levels are set to
   a new position obtained by the center of gravity of their neighboring
   nodes. The processing from bottom to top level ensures fast convergence
   of the algorithm. Caution! The algorithm may produce undesirable
   results for non-convex domains.

   `This function is available in 2D version only!`

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 if error occured.
   D*/
/****************************************************************************/

#ifdef __TWODIM__
INT SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag)
{
  INT l,i,n,sides;
  DOUBLE ratio;
  DOUBLE N;
  COORD lambda,lambda0,lambda1,lambda2,lambdaa,lambdae;
  GRID *theGrid;
  NODE *node,*node2;
  ELEMENT *eptr;
  ELEMENTSIDE *eside;
  VERTEX *vptr0,*vptr1,*vptr2,*vptra,*vptre,*vptr;
  LINK *lptr;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;

  COORD x[2];

  n = niter;
  if (n<=0) n = 1;
  if (n>50) n = 50;

  ratio=.5;

  for (i=0; i<n; i++)
  {
    for (l=0; l<=theMG->topLevel; l++)
    {
      theGrid=theMG->grids[l];

      /* update global coordinates of new nodes */
      if (l!=0)
        for (node=theGrid->firstNode; node!=NULL; node=SUCCN(node))
          if (NFATHER(node)==NULL)
          {
            vptr=MYVERTEX(node);
            if ((OBJT(vptr)!=BVOBJ)||(bdryFlag!=0))
              if (Local2Global(theMG,vptr)!=0)
                return(GM_ERROR);
          }

      for (node=theGrid->firstNode; node!=NULL; node=SUCCN(node))
      {
        /* skip node if it is a copy from a lower level */
        if (NFATHER(node)!=NULL) continue;

        vptr0=MYVERTEX(node);
        if (OBJT(vptr0)==BVOBJ)
        {
          if (bdryFlag==0) continue;                                            /* boundary is not allowed to move */
          if (MOVE(vptr0)==0) continue;                                 /* corner: do not move it */

          /* test if free boundary: since that type is determined by the grid it may only
             be smoothed in certain cases and in a special way */
          thePatch=VS_PATCH(VSEG(vptr0));
          if(Patch_GetPatchDesc(thePatch,&thePatchDesc)) return(1);
          if (PATCH_TYPE(thePatchDesc)==FREE)
            continue;                                           /*free boundary is not allowed to be smoothed */

          /* first find endpoints of vptr0's boundary segment */
          lambda0=LAMBDA(VSEG(vptr0),0);
          vptra=theMG->corners[PATCH_CID(thePatchDesc,0)]; lambdaa=PATCH_LCVECT(thePatchDesc,0)[0];
          vptre=theMG->corners[PATCH_CID(thePatchDesc,1)]; lambdae=PATCH_LCVECT(thePatchDesc,1)[0];

          /* search for the two nearest neighbors on the same segment */
          vptr1=vptr2=NULL;
          for (lptr=START(node); lptr!=NULL; lptr=NEXT(lptr))
            if (EXTRA(MYEDGE(lptr))==0)
            {
              node2=NBNODE(lptr);
              vptr=MYVERTEX(node2);
              if (OBJT(vptr)!=BVOBJ)
                continue;

              /* if on same segment get lambda value */
              if (vptr==vptra)
                lambda=lambdaa;
              else
              if (vptr==vptre)
                lambda=lambdae;
              else
              if (thePatch==VS_PATCH(VSEG(vptr)))
                lambda=LAMBDA(VSEG(vptr),0);
              else
                continue;

              if (lambda>lambda0)
              {
                if (vptr2==NULL)
                {
                  vptr2=vptr;
                  lambda2=lambda;
                }
                else
                if (lambda<lambda2)
                {
                  vptr2=vptr;
                  lambda2=lambda;
                }
              }
              else
              {
                if (vptr1==NULL)
                {
                  vptr1=vptr;
                  lambda1=lambda;
                }
                else
                if (lambda>lambda1)
                {
                  vptr1=vptr;
                  lambda1=lambda;
                }
              }
            }

          if ((vptr1==NULL)||(vptr2==NULL))
            return(GM_ERROR);

          if (PATCH_TYPE(thePatchDesc)==FREE)
          {
            /*	This is only sensible if the free boundary is a line and the endpoints are moved.
                    A more general smoothing would be to use e.g. a quadratic interpolant. */
            XC(vptr0)=.5*(XC(vptr1)+XC(vptr2));
            YC(vptr0)=.5*(YC(vptr1)+YC(vptr2));
            /*	since local and global boundary coordinates are given by the grid itself
                    they are assumed to be correct. */
          }
          else
          {
            LAMBDA(VSEG(vptr0),0)=.5*(lambda1+lambda2);
            /* set global coordinates */
            if (Patch_local2global(thePatch,PVECT(VSEG(vptr0)),CVECT(vptr0))) RETURN (GM_ERROR);

            /* set local boundary coordinates */
            if (Global2Local(theMG,vptr0)!=0)
              RETURN(GM_ERROR);
          }
        }
        else
        {
          x[0]=x[1]=0; N=0;

          for (lptr=START(node); lptr!=NULL; lptr=NEXT(lptr))
          {
            node2=NBNODE(lptr);
            vptr=MYVERTEX(node2);

            if (EXTRA(MYEDGE(lptr))==0)
            {
              x[0]+=XC(vptr);
              x[1]+=YC(vptr);
              N+=1;
            }
            else
            {
              x[0]+=ratio*XC(vptr);
              x[1]+=ratio*YC(vptr);
              N+= ratio;
            }
          }

          XC(vptr0)=x[0]/N; YC(vptr0)=x[1]/N;

          /* if there is a father element, change local variables */
          if (l!=0)
            if ((eptr=FindFather(vptr0))!=NULL)
            {
              VFATHER(vptr0)=eptr;
              if (Global2Local(theMG,vptr0)!=0)
                RETURN(GM_ERROR);
            }
            else
              RETURN(GM_ERROR);
        }
      }

      /* at last, the boundary element sides of this level must be updated! */
      if (bdryFlag==0) continue;                        /* boundary is not allowed to move */
      for (eptr=theGrid->elements; eptr!=NULL; eptr=SUCCE(eptr))
        if (OBJT(eptr)==BEOBJ)
        {
          sides=SIDES_OF_ELEM(eptr);
          for (i=0; i<sides; i++)
            if ((eside=SIDE(eptr,i))!=NULL)
            {
              /* for higher dimensions (3d) the following must be generalized */
              vptr=MYVERTEX(CORNER(eptr,i));
              if (MOVE(vptr)!=0)
                PARAM(eside,0,0)=LAMBDA(VSEG(vptr),0);
              vptr=MYVERTEX(CORNER(eptr,(i+1)%sides));
              if (MOVE(vptr)!=0)
                PARAM(eside,1,0)=LAMBDA(VSEG(vptr),0);
            }
        }
    }
  }

  return(GM_OK);
}
#endif

#ifdef __THREEDIM__
INT  SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag)
{
  INT l,i,n,sides;
  DOUBLE ratio;
  DOUBLE N;
  COORD lambda,lambda0,lambda1,lambda2,lambdaa,lambdae;
  GRID *theGrid;
  NODE *node,*node2;
  ELEMENT *eptr;
  VERTEX *vptr0,*vptr1,*vptr2,*vptra,*vptre,*vptr;
  LINK *lptr;
  COORD_VECTOR x;

  n = niter;
  if (n<=0) n = 1;
  if (n>50) n = 50;

  ratio=.5;

  if (theMG->topLevel > 0)
  {
    UserWrite("SmoothMultiGrid not implemented in 3D version jet\n");
    return (0);
  }

  theGrid=theMG->grids[0];
  for (i=0; i<n; i++)
  {
    for (node=theGrid->firstNode; node!=NULL; node=SUCCN(node))
    {
      /* skip node if it is a copy from a lower level */
      if (NFATHER(node)!=NULL) continue;

      vptr0=MYVERTEX(node);
      {
        x[0]=x[1]=0; N=0;

        for (lptr=START(node); lptr!=NULL; lptr=NEXT(lptr))
        {
          node2=NBNODE(lptr);
          vptr=MYVERTEX(node2);

          if (EXTRA(MYEDGE(lptr))==0)
          {
            x[0]+=XC(vptr);
            x[1]+=YC(vptr);
            x[2]+=ZC(vptr);
            N+=1;
          }
          else
          {
            x[0]+=ratio*XC(vptr);
            x[1]+=ratio*YC(vptr);
            x[2]+=ratio*ZC(vptr);
            N+= ratio;
          }
        }

        XC(vptr0)=x[0]/N; YC(vptr0)=x[1]/N; ZC(vptr0)=x[2]/N;
      }
    }
  }

  return(GM_OK);
}
#endif

/****************************************************************************/
/*D
   InsertElement - Insert an element

   SYNOPSIS:
   InsertElement (MULTIGRID *theMG, INT n, NODE **Node);

   PARAMETERS:
   .  theMG - multigrid structure
   .  tag - type of element to insert
   .  idList - array of ids of corner nodes

   DESCRIPTION:
   This function inserts an element into level 0.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
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
  COORD_VECTOR diff[3],rot;
  COORD det;
  INT i;

  /* TODO: this case */
  if (n == 8)
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

INT InsertElement (MULTIGRID *theMG, INT n, NODE **Node)
{
  GRID             *theGrid;
  INT i,j,k,l,m,found,num,tag,ElementType;
  INT NeighborSide[MAX_SIDES_OF_ELEM];
  NODE             *sideNode[MAX_CORNERS_OF_SIDE],*NeighborNode,*theNode;
  VERTEX           *Vertex[MAX_CORNERS_OF_ELEM],*sideVertex[MAX_CORNERS_OF_SIDE],*theVertex;
  ELEMENT          *theElement,*Neighbor[MAX_SIDES_OF_ELEM];
  EDGE             *theEdge;
  COORD param[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE][DIM_OF_BND];
  COORD        *plambda[MAX_CORNERS_OF_SIDE];
  VSEGMENT         *vs,*vs1;
  PATCH            *thePatch[MAX_SIDES_OF_ELEM], *Patch;
  PATCH_DESC thePatchDesc;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertElement",
                      "only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* check parameters */
    #ifdef __TWODIM__
  if ((n!=TRIANGLE)&&(n!=QUADRILATERAL))
  {
    PrintErrorMessage('E',"InsertElement",
                      "only triangles and quadrilaterals allowed in 2D");
    RETURN(GM_ERROR);
  }
  tag = n;
        #endif

    #ifdef __THREEDIM__
  if (n == 4)
    tag = TETRAHEDRON;
  else if ( n == 8)
    tag = HEXAHEDRON;
  else
  {
    PrintErrorMessage('E',"InsertElement",
                      "only tetrahedrons and hexahedrons are allowed in the 3D coarse grid");
    RETURN(GM_ERROR);
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
              RETURN(GM_ERROR);
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
    thePatch[i] = NULL;
  }

  /* compute side information (theSeg[i]==NULL) means inner side */
  ElementType = IEOBJ;
  for (i=0; i<SIDES_OF_REF(n); i++)
  {
    for(j=0; j<CORNERS_OF_SIDE_REF(n,i); j++ )
    {
      k = CORNER_OF_SIDE_REF(n,i,j);
      sideNode[j] = Node[k];
      sideVertex[j] = Vertex[k];
    }
    found = 0;
    for(j=0; j<CORNERS_OF_SIDE_REF(n,i); j++ )
    {
      if( OBJT(sideVertex[j]) == IVOBJ ) found = 1;
    }
    if( found ) continue;

    /* all vertices of side[i] are on the boundary now */

    /* We now assume, that side[i] is on the boundary if and only if */
    /* there is one boundary segment containing the three nodes.        */
    /* That means, one should not define elements at the boundary	*/
    /* with a boundary side covering more than one segment.			*/

    for (vs1 = VSEG(sideVertex[0]); vs1 != NULL; vs1 = NEXTSEG(vs1))
    {
      Patch = VS_PATCH(vs1);
      plambda[0] = PVECT(vs1);
      for( k=1; k<CORNERS_OF_SIDE_REF(n,i); k++ )
      {
        found = 0;
        /* the segments of that corner */
        for( vs = VSEG(sideVertex[k]); vs!=NULL; vs = NEXTSEG(vs) )
        {
          if (Patch == VS_PATCH(vs))
          /* corner belongs to segment */
          {
            found = 1;
            break;
          }
        }
        if( !found      )
          break;
        plambda[k] = PVECT(vs);
      }
      if (found)
        break;
    }

    if (found)             /* a commen patch for all vertices on a side was found
                              i.e. the side is a boundary side */
    {
      /*	set boundary parameters for vertices of that side */
      ElementType = BEOBJ;
      thePatch[i] = Patch;
      for( k=0; k<CORNERS_OF_SIDE_REF(n,i); k++)
        for (j=0; j<DIM-1; j++)
          param[i][k][j] = plambda[k][j];
    }
  }

  /* find neighboring elements */

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
            PrintErrorMessage('E',"InsertElement","neighbor relation inconsistent");
            RETURN(GM_ERROR);
          }
          Neighbor[i] = theElement;
          NeighborSide[i] = j;
        }
      }
    }
  }

  /* create element */
  theElement = CreateElement(theGrid,tag,ElementType,Node,NULL);
  if (theElement==NULL)
  {
    PrintErrorMessage('E',"InsertElement","cannot allocate element");
    RETURN(GM_ERROR);
  }

  /* create element sides if necessary */
  if (OBJT(theElement)==BEOBJ)
  {
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      SET_SIDE(theElement,i,NULL);
      if (thePatch[i]!=NULL)
      {
        SET_SIDE(theElement,i,CreateElementSide(theGrid));
        if (SIDE(theElement,i)==NULL)
        {
          DisposeElement(theGrid,theElement);
          PrintErrorMessage('E',"InsertElement","cannot allocate element side");
          RETURN(GM_ERROR);
        }
        ES_PATCH(SIDE(theElement,i)) = thePatch[i];
        for(k=0; k<CORNERS_OF_SIDE(theElement,i); k++)
          for(l=0; l<DIM_OF_BND; l++)
            PARAM(SIDE(theElement,i),k,l) = param[i][k][l];
      }
    }
  }

  /* fill element data */
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    SET_NBELEM(theElement,i,Neighbor[i]);
    if (Neighbor[i]!=NULL)
    {
      SET_NBELEM(Neighbor[i],NeighborSide[i],theElement);
            #ifdef __THREEDIM__
      if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
        if (DisposeDoubledSideVector(theGrid,Neighbor[i],NeighborSide[i],theElement,i))
          RETURN(GM_ERROR);
            #endif
    }
  }
  SETNSONS(theElement,0);
  SET_SON(theElement,0,NULL);
  SET_EFATHER(theElement,NULL);
  SETECLASS(theElement,RED);

  /* create connection to other elements. ATTENTION: this function is O(n)*/
  if (InsertedElementCreateConnection(theGrid,theElement))
  {
    PrintErrorMessage('E',"InsertElement","could not create algebra connections");
    DisposeElement (theGrid,theElement);
    RETURN(GM_ERROR);
  }

  return(GM_OK);
}

/****************************************************************************/
/*D
   InsertElementFromIDs - Insert element with node ids

   SYNOPSIS:
   INT InsertElementFromIDs (MULTIGRID *theMG, INT n, INT *idList);

   PARAMETERS:
   .  theMG - multigrid structure
   .  n - number of nodes in node id list
   .  idList - ids of the nodes

   DESCRIPTION:
   This function inserts an element with nodes that have the ids
   given in `idList`,  on level 0.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR when error occured.
   D*/
/****************************************************************************/

INT InsertElementFromIDs (MULTIGRID *theMG, INT n, INT *idList)
{
  GRID *theGrid;
  NODE *Node[MAX_CORNERS_OF_ELEM],*theNode;
  INT i,j,found;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertElementFromIDs","only a multigrid with exactly one level can be edited");
    RETURN(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* check data */
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++)
      if (idList[i]==idList[j])
      {
        PrintErrorMessage('E',"InsertElementFromIDs","nodes must be pairwise different");
        RETURN(GM_ERROR);
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
    PrintErrorMessage('E',"InsertElementFromIDs","could not find all nodes");
    RETURN(GM_ERROR);
  }

  return (InsertElement(theMG,n,Node));
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
  LINK *theLink;
  EDGE *theEdge;
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
  DisposeElement(theGrid,theElement);

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
   NODE *FindNodeFromPosition (GRID *theGrid, COORD *pos, COORD *tol);

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

NODE *FindNodeFromPosition (GRID *theGrid, COORD *pos, COORD *tol)
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
   VECTOR *FindVectorFromPosition (GRID *theGrid, COORD *pos, COORD *tol);

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

VECTOR *FindVectorFromPosition (GRID *theGrid, COORD *pos, COORD *tol)
{
  VECTOR *theVector;
  COORD_VECTOR vpos;
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
   INT PointInElement (const COORD *x, const ELEMENT *theElement);

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
INT PointInElement (const COORD *x, const ELEMENT *theElement) /* 2D version */
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
  thePoint.x = (COORD)x[0];
  thePoint.y = (COORD)x[1];

  return(PointInPolygon(point,n,thePoint));
}
#endif

#ifdef __THREEDIM__
INT PointInElement (const COORD *global, const ELEMENT *theElement)
{
  COORD *x[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR a,b,rot;
  COORD det;
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
   ELEMENT *FindElementFromPosition (GRID *theGrid, COORD *pos)

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

ELEMENT *FindElementFromPosition (GRID *theGrid, COORD *pos)
{
  ELEMENT *theElement;

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
    if (PointInElement(pos,theElement) == 1)
      return(theElement);

  return(NULL);
}

ELEMENT *FindElementOnSurface (MULTIGRID *theMG, COORD *global)
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
  if (BVP_GetBVPDesc(theBVP,&theBVPDesc))
  {
    PrintErrorMessage('E',"InsertElement","cannot evaluate BVP");
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
  COORD hmin,hmax,h;
  INT l,cl,minl,i,soe,eos,coe,side,e;
  INT nn,ne,nt,ns,nvec,nc;

  cl = CURRENTLEVEL(theMG);

  sprintf(buffer,"grids of '%s':\n",ENVITEM_NAME(theMG));

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

    sprintf(buffer,"%c %3d %8d %8ld %8ld %8ld %8ld %8ld %8ld %8ld %9.3e %9.3e\n",c,l,(int)TOPLEVEL(theMG),
            (long)NV(theGrid),(long)NN(theGrid),(long)NE(theGrid),(long)NT(theGrid),
            (long)NS(theGrid),(long)NVEC(theGrid),(long)NC(theGrid),(float)hmin,(float)hmax);

    UserWrite(buffer);
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
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
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
            if (SIDE(theElement,side)!=NULL) ns++;

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
  sprintf(buffer,"\n%lu bytes used out of %lu allocated\n",HeapUsed(MGHEAP(theMG)),HeapSize(MGHEAP(theMG)));
  UserWrite(buffer);
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
  VSEGMENT *vs;
  LINK *theLink;
  int i;
  PATCH_DESC thePatchDesc;

  theVertex = MYVERTEX(theNode);

  /******************************/
  /* print standard information */
  /******************************/
  /* line 1 */ sprintf(buffer,"NODEID=%9ld CTRL=%8lx IX=%8ld VEID=%9ld LEVEL=%2d",(long)ID(theNode),(long)CTRL(theNode),
                       (long)INDEX(theNode),(long)ID(theVertex),LEVEL(theNode));
  UserWrite(buffer);

        #ifdef ModelP
  sprintf(buffer," NGID=%08x VGID=%08x",
          DDD_InfoGlobalId(PARHDR(theNode)),
          DDD_InfoGlobalId(PARHDRV(theVertex))
          );
  UserWrite(buffer);
        #endif

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
    if (OBJT(theVertex)==BVOBJ)
    {
      /* print vertex boundary information */
      for(vs=VSEG(theVertex); vs!=NULL; vs=NEXTSEG(vs))
      {
        if (Patch_GetPatchDesc(VS_PATCH(vs),&thePatchDesc)) return;
        sprintf(buffer,"   PATID=%9ld ",(long)PATCH_ID(thePatchDesc));
        UserWrite(buffer);
        for(i=0; i<DIM_OF_BND; i++)
        {
          sprintf(buffer,"LAMBDA[%d]=%11.4E ",i,(float)LAMBDA(vs,i) );
          UserWrite(buffer);
        }
        UserWrite("\n");
      }
    }
  }

  if (nbopt)
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      sprintf(buffer,"   NB=%5ld CTRL=%8lx ",(long)ID(NBNODE(theLink)),(long)CTRL(theLink));
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
    for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,level)); theNode!=NULL; theNode=SUCCN(theNode))
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
  int i,j,k;
  ELEMENTSIDE *theSide;
        #ifdef __THREEDIM__
  ELEMENT *SonList[MAX_SONS];
        #endif
  PATCH_DESC thePatchDesc;

  if (ECLASS(theElement)==COPY_CLASS) strcpy(etype,"COPY ");
  if (ECLASS(theElement)==IRREGULAR_CLASS) strcpy(etype,"IRREG");
  if (ECLASS(theElement)==REGULAR_CLASS) strcpy(etype,"REGUL");
  sprintf(buffer,"ELEMID=%9ld %5s CTRL=%8lx CTRL2=%8lx REFINE=%2d MARK=%2d LEVEL=%2d",(long)ID(theElement),etype,
          (long)CTRL(theElement),(long)FLAG(theElement),REFINE(theElement),MARK(theElement),LEVEL(theElement));
  UserWrite(buffer);
  if (COARSEN(theElement)) UserWrite(" COARSEN");

        #ifdef ModelP
  sprintf(buffer," EGID=%08x",DDD_InfoGlobalId(PARHDRE(theElement)));
  UserWrite(buffer);
        #endif

  UserWrite("\n");

  if (vopt)
  {
    UserWrite("   ");
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      sprintf(buffer,"N%d=%ld ",i,(long)ID(CORNER(theElement,i)));
      UserWrite(buffer);
    }
    if (EFATHER(theElement)!=NULL)
    {
      sprintf(buffer,"FA=%ld ",(long)ID(EFATHER(theElement)));
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
        if ((theSide=SIDE(theElement,i))!=NULL)
        {
          if (Patch_GetPatchDesc(ES_PATCH(theSide),&thePatchDesc)) return;
          sprintf(buffer,"PATID%d=%ld ",i,(long)PATCH_ID(thePatchDesc));
          UserWrite(buffer);
          for(j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
                                                #ifdef __THREEDIM__
                        #ifdef ModelP
            sprintf(buffer,"  NODE[ID=%ld]: ",(long)(ID(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j)))));
            UserWrite(buffer);
                                                #endif
                                                #endif
            for(k=0; k<DIM_OF_BND; k++)
            {
              sprintf(buffer,"LAMBDA[%d]=%.10g ",k,(float)(PARAM(theSide,j,k)) );
              UserWrite(buffer);
            }
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
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,level)); theElement!=NULL; theElement=SUCCE(theElement))
      {
        if ( (ID(theElement)>=from)&&(ID(theElement)<=to) )
          ListElement(theMG,theElement,dataopt,bopt,nbopt,vopt);
      }
  else
    for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG))); theElement!=NULL; theElement=SUCCE(theElement))
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
    sprintf(buffer,"NODE-V IND=%7ld nodeID=%7ld                VCLASS=%1d VNCLASS=%1d\n",VINDEX(theVector),ID(theNode),VCLASS(theVector),VNCLASS(theVector));
    UserWrite(buffer);
  }
  if (VTYPE(theVector)==EDGEVECTOR)
  {
    theEdge = (EDGE*)VOBJECT(theVector);
    sprintf(buffer,"EDGE-V IND=%7ld fromID=%7ld to__ID=%7ld VCLASS=%1d VNCLASS=%1d\n",VINDEX(theVector),ID(NBNODE(LINK0(theEdge))),ID(NBNODE(LINK1(theEdge))),VCLASS(theVector),VNCLASS(theVector));
    UserWrite(buffer);
  }
        #ifdef __THREEDIM__
  if (VTYPE(theVector)==SIDEVECTOR)
  {
    theElement = (ELEMENT*)VOBJECT(theVector);
    sprintf(buffer,"SIDE-V IND=%7ld elemID=%ld                VCLASS=%1d VNCLASS=%1d\n",VINDEX(theVector),ID(theElement),VCLASS(theVector),VNCLASS(theVector));
    UserWrite(buffer);
  }
        #endif
  if (VTYPE(theVector)==ELEMVECTOR)
  {
    theElement = (ELEMENT*)VOBJECT(theVector);
    sprintf(buffer,"ELEM-V IND=%7ld elemID=%ld                VCLASS=%1d VNCLASS=%1d\n",VINDEX(theVector),ID(theElement),VCLASS(theVector),VNCLASS(theVector));
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
    for (theVector=FIRSTVECTOR(GRID_ON_LEVEL(theMG,level)); theVector!=NULL; theVector=SUCCVC(theVector))
    {
      if (VINDEX(theVector)>=from && VINDEX(theVector)<=to)
        ListVector(theMG,theVector,matrixopt,dataopt);
    }
}

/****************************************************************************/
/*D
   CheckGrid - Check consistency of data structure

   SYNOPSIS:
   INT CheckGrid (GRID *theGrid);

   PARAMETERS:
   .  theGrid - grid to check

   DESCRIPTION:
   This function checks the consistency of data structure.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if an error occured.
   D*/
/****************************************************************************/

#ifdef __TWODIM__
static INT CheckElement (ELEMENT *theElement, INT *SideError, INT *EdgeError, INT *NodeError)
{
  int i,j,k,n;
  EDGE *theEdge;
  ELEMENT *NbElement;
  ELEMENTSIDE *theSide;
  VERTEX *theVertex;
  VSEGMENT *vs;
  PATCH *thePatch;

  *SideError = 0;
  *NodeError = 0;
  n = SIDES_OF_ELEM(theElement);
  if (ECLASS(theElement)!=YELLOW)
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
          if ((theSide=SIDE(theElement,i)) != NULL)
          {
            for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
            {
              theVertex = MYVERTEX(CORNER(theElement,(k=(i+j)%n)));
              if (OBJT(theVertex) == BVOBJ)
              {
                thePatch = ES_PATCH(theSide);
                for (vs=VSEG(theVertex); vs!=NULL; vs=NEXTSEG(vs))
                  if (thePatch == VS_PATCH(vs))
                    break;
                if (vs == NULL)
                  *NodeError |= (1<<k);
              }
              else
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

INT CheckGrid (GRID *theGrid) /* 2D VERSION */
{
  NODE *theNode;
  ELEMENT *theElement;
  LINK *theLink;
  int i,j;
  EDGE *theEdge;
  INT SideError, EdgeError, NodeError,count,n;

  /* reset used flags
     for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
     {
     #ifdef Debug
     #ifdef ModelP
          printf("%2d: CheckGrid:       n=%x\n",me,theNode);
     #endif
     #endif

          SETUSED(theNode,0);
          for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
            PARSETUSED(MYEDGE(theLink),0);
          }
     }

     #ifdef Debug
     #ifdef ModelP
     printf("%2d: CheckGrid: all nodes and links are SETUSED\n",me);
     #endif
   #endif */

  /* check neighbors and edges of elements */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    n = SIDES_OF_ELEM(theElement);
    if (CheckElement(theElement, &SideError, &EdgeError, &NodeError)==0) continue;
    sprintf(buffer,"ELEM %ld:\n",(long)ID(theElement));
    UserWrite(buffer);
    if (SideError)
      for (i=0; i<n; i++)
      {
        if (SideError & 1<<i)
        {
          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,(i+j)%n)));
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
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,(i+j)%n)));
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
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,(i+j)%n)));
            UserWrite(buffer);
            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has no neighbour, element is BEOBJ but there is no SIDE\n");
        }
      }
    if (EdgeError)
      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {
        if (!(EdgeError & 1<<i)) continue;
        sprintf(buffer,"   EDGE(%ld,%ld) is missing\n",(long)ID(CORNER(theElement,i)),(long)ID(CORNER(theElement,(i+1)%n)));
        UserWrite(buffer);
      }
    if (NodeError)
      for (i=0; i<n; i++)
      {
        if (NodeError & (1<<i))
        {
          sprintf(buffer,"   CORNER %ld is BVOBJ, ids from elementside and vertexsegment are not consistent\n",(long)ID(CORNER(theElement,i)));
          UserWrite(buffer);
        }
        if (NodeError & (1<<(i+n)))
        {
          sprintf(buffer,"   CORNER %ld is IVOBJ, but lies on elementside\n",(long)ID(CORNER(theElement,i)));
          UserWrite(buffer);
        }
      }
  }

  /* look for dead edges */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      PARSETUSED(theLink,1);
    }
  }
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
      if ((PARUSED(theLink)&&!PARUSED(REVERSE(theLink))) ||
          (PARUSED(REVERSE(theLink))&&!PARUSED(theLink)))
      {
                #ifdef ModelP
        sprintf(buffer,"edge between %ld and %ld dead ",(long)ID(theNode),(long)ID(NBNODE(theLink)));
                                #endif
        UserWrite(buffer);
        UserWrite("\n");
      }
      else
      {
        PARSETUSED(theLink,0);
        PARSETUSED(REVERSE(theLink),0);
      }
    }
  }

        #ifdef Debug
        #ifdef ModelP
  printf("%2d: CheckGrid(): no dead edges found!\n",me);
        #endif
        #endif

  /* look for dead nodes */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    if (!USED(theNode))
    {
            #ifdef ModelP
      sprintf(buffer,"node %ld is dead ",(long)ID(theNode));
      UserWrite(buffer);
      UserWrite("\n");
                        #endif
    }
    else
      SETUSED(theNode,0);
  }

        #ifdef Debug
        #ifdef ModelP
  printf("%2d: CheckGrid(): no dead nodes found!\n",me);
        #endif
        #endif

  /* check number of elem and their pointers */
  count = 0;
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (SUCCE(theElement)!=NULL)
    {
      if (OBJT(SUCCE(theElement))!=IEOBJ && OBJT(SUCCE(theElement))!=BEOBJ)
      {
        sprintf(buffer,"pointer of ELEM(%ld) (number %ld) to next element is no pointer to an element\n",(long)ID(theElement),(long)count);
        UserWrite(buffer);
        break;
      }
      if (PREDE(SUCCE(theElement))!=NULL)
      {
        if (PREDE(SUCCE(theElement))!=theElement)
        {
          sprintf(buffer,"pointer of ELEM(%ld) (number %ld) to previous element is not the previous element\n",(long)ID(SUCCE(theElement)),(long)(count+1));
          UserWrite(buffer);
        }
      }
      else
      {
        sprintf(buffer,"pointer of ELEM(%ld) (number %ld) to previous element is NULL\n",(long)ID(SUCCE(theElement)),(long)(count+1));
        UserWrite(buffer);
      }
    }
    count++;
  }
  if (PREDE(theGrid->elements) != NULL)
  {
    sprintf(buffer,"first element of the grid has a previous 'element'\n");
    UserWrite(buffer);
  }
  if (SUCCE(theGrid->lastelement) != NULL)
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
#endif

#ifdef __THREEDIM__
static INT CheckElement (ELEMENT *theElement, INT *SideError, INT *EdgeError, INT *NodeError, INT *ESonError, INT *NSonError)
{
  int i,j,k;
  EDGE *theEdge;
  ELEMENT *NbElement;
  ELEMENT *SonList[MAX_SONS];
  ELEMENTSIDE *theSide;
  VERTEX *theVertex;
  VSEGMENT *vs;
  PATCH *thePatch;

  *SideError = 0;
  *NodeError = 0;
  *EdgeError = 0;
  *ESonError = 0;
  *NSonError = 0;
  if (ECLASS(theElement)!=YELLOW)
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
      }
      else
      {
        if (OBJT(theElement) == IEOBJ)
          *SideError |= (1<<(i+MAX_SIDES_OF_ELEM));
        else
        {
          if ((theSide=SIDE(theElement,i)) != NULL)
          {
            for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
            {
              theVertex = MYVERTEX(CORNER(theElement,(k=CORNER_OF_SIDE(theElement,i,j))));
              if (OBJT(theVertex) == BVOBJ)
              {
                thePatch = ES_PATCH(theSide);
                for (vs=VSEG(theVertex); vs!=NULL; vs=NEXTSEG(vs))
                  if (thePatch == VS_PATCH(vs))
                    break;
                if (vs == NULL)
                  *NodeError |= (1<<k);
              }
              else
                *NodeError |= (1<<(k+MAX_CORNERS_OF_ELEM));
            }
          }
          else
            *SideError |= (1<<(i+2*MAX_SIDES_OF_ELEM));
        }
      }
    }
  else
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
      }
      else
      {
        if (OBJT(theElement) == BEOBJ)
        {
          if ((theSide=SIDE(theElement,i)) != NULL)
          {
            for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
            {
              theVertex = MYVERTEX(CORNER(theElement,(k=CORNER_OF_SIDE(theElement,i,j))));
              if (OBJT(theVertex) == BVOBJ)
              {
                thePatch = ES_PATCH(theSide);
                for (vs=VSEG(theVertex); vs!=NULL; vs=NEXTSEG(vs))
                  if (thePatch == VS_PATCH(vs))
                    break;
                if (vs == NULL)
                  *NodeError |= (i<<k);
              }
              else
                *NodeError |= (i<<(k+MAX_CORNERS_OF_ELEM));
            }
          }
        }
      }
    }

  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    SETUSED(CORNER(theElement,i),1);
  }

  for (i=0; i<EDGES_OF_ELEM(theElement); i++)
  {
    theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
    if (theEdge==NULL)
      *EdgeError |= 1<<i;
    else
      PARSETUSED(theEdge,1);
  }

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
        sprintf(buffer,"ELEM(ID=%d): element is not refined but has NSONS=%d\n",ID(theElement),NSONS(theElement));
        UserWrite(buffer);
      }
    }
  }

  if (*SideError || *EdgeError || *NodeError || *ESonError || *NSonError)
    return (1);

  return (0);
}


INT CheckGrid (GRID *theGrid) /* 3D VERSION */
{
  NODE *theNode;
  ELEMENT *theElement;
  LINK *theLink;
  int i,j;
  EDGE *theEdge;
  INT SideError, EdgeError, NodeError, ESonError, NSonError, count;

  /* reset used flags */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    SETUSED(theNode,0);
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
      PARSETUSED(MYEDGE(theLink),0);
  }

  /* check neighbors and edges of elements */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (CheckElement(theElement,&SideError,&EdgeError,&NodeError,&ESonError,&NSonError)==0) continue;
    sprintf(buffer,"ELEM %ld:\n",(long)ID(theElement));
    UserWrite(buffer);
    if (SideError)
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        if (SideError & 1<<i)
        {
          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))));
            UserWrite(buffer);
            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has neighbour but a backPtr does not exist\n");
        }
        if (SideError & 1<<(i+MAX_SIDES_OF_ELEM))
        {
          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))));
            UserWrite(buffer);
            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has no neighbour but element is IEOBJ\n");
        }
        if (SideError & 1<<(i+2*MAX_SIDES_OF_ELEM))
        {
          UserWrite("   SIDE(");
          for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,CORNER_OF_SIDE(theElement,i,j))));
            UserWrite(buffer);
            if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
          }
          UserWrite(") has no neighbour, element is BEOBJ but there is no SIDE\n");
        }
      }
    if (EdgeError)
      for (i=0; i<EDGES_OF_ELEM(theElement); i++)
      {
        if (!(EdgeError & 1<<i)) continue;
                                #ifdef ModelP
        sprintf(buffer,"   EDGE(%ld,%ld) is missing\n",(long)ID(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0))),(long)ID(CORNER(theElement,CORNER_OF_EDGE(theElement,i,1))));
        UserWrite(buffer);
                                #endif
      }

    if (NodeError)
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        if (NodeError & (1<<i))
        {
          sprintf(buffer,"   CORNER %ld is BVOBJ, ids from elementside and vertexsegment are not consistent\n",(long)ID(CORNER(theElement,i)));
          UserWrite(buffer);
        }
        if (NodeError & (1<<(i+MAX_CORNERS_OF_ELEM)))
        {
          sprintf(buffer,"   CORNER %ld is IVOBJ, but lies on elementside\n",(long)ID(CORNER(theElement,i)));
          UserWrite(buffer);
        }
      }
    if (ESonError)
    {
      for (i=0; i<NSONS(theElement); i++)
      {
        if ((ESonError & 1<<i))
        {
          sprintf(buffer,"   ESON(%d) has wrong EFATHER pointer\n",i);
          UserWrite(buffer);
        }
      }
    }
    if (NSonError)
    {
      for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
      {
        if (NSonError & (1<<i))
        {
          sprintf(buffer,"   SONNODE(CORNER %d) != CORNER(ESON)\n",i);
          UserWrite(buffer);
        }
        if (NSonError & (1<<(i+MAX_CORNERS_OF_ELEM)))
        {
          sprintf(buffer,"   CORNER %d != EFATHER(CORNER(ESON))\n",i);
          UserWrite(buffer);
        }
      }

      for (i=0; i<MAX_EDGES_OF_ELEM; i++)
      {

        if (NSonError & (1<<(i+2*MAX_CORNERS_OF_ELEM)))
        {
          sprintf(buffer,"   MIDNODE(edge %d) != CORNER(ESON)\n",i);
          UserWrite(buffer);
        }
      }
      if (NSonError & (1<<(6+2*MAX_CORNERS_OF_ELEM)))
      {
        sprintf(buffer,"   NFATHER(CENTERNODE(ESON)) != NULL\n");
        UserWrite(buffer);
      }
    }
  }

  /* look for dead edges */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
    {
                    #ifdef ModelP
      if (ID(NBNODE(theLink))>ID(theNode))
      {
        theEdge = MYEDGE(theLink);
        if (!PARUSED(theEdge))
        {
          sprintf(buffer,"edge between %ld and %ld dead, NO_OF_ELEM=%d ",(long)ID(theNode),(long)ID(NBNODE(theLink)),NO_OF_ELEM(theEdge));
          UserWrite(buffer);
          UserWrite("\n");
        }
        else
          PARSETUSED(theEdge,0);
      }
                        #endif
    }
  }

  /* look for dead nodes */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    if (!USED(theNode))
    {
                    #ifdef ModelP
      sprintf(buffer,"node %ld is dead ",(long)ID(theNode));
      UserWrite(buffer);
      UserWrite("\n");
                        #endif
    }
    else
      SETUSED(theNode,0);
  }

  /* check number of elem and their pointers */
  count = 0;
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (SUCCE(theElement)!=NULL)
    {
      if (OBJT(SUCCE(theElement))!=IEOBJ && OBJT(SUCCE(theElement))!=BEOBJ)
      {
        sprintf(buffer,"pointer of ELEM(%ld) (number %ld) to next element is no pointer to an element\n",(long)ID(theElement),(long)count);
        UserWrite(buffer);
        break;
      }
      if (PREDE(SUCCE(theElement))!=NULL)
      {
        if (PREDE(SUCCE(theElement))!=theElement)
        {
          sprintf(buffer,"pointer of ELEM(%ld) (number %ld) to previous element is not the previous element\n",(long)ID(SUCCE(theElement)),(long)(count+1));
          UserWrite(buffer);
        }
      }
      else
      {
        sprintf(buffer,"pointer of ELEM(%ld) (number %ld) to previous element is NULL\n",(long)ID(SUCCE(theElement)),(long)(count+1));
        UserWrite(buffer);
      }
    }
    count++;
  }
  if (PREDE(theGrid->elements) != NULL)
  {
    sprintf(buffer,"first element of the grid has a previous 'element'\n");
    UserWrite(buffer);
  }
  if (SUCCE(theGrid->lastelement) != NULL)
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
#endif


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
    if (SELECTIONOBJECT(theMG,i)==g) RETURN(GM_ERROR);

  if (SELECTIONSIZE(theMG)>=MAXSELECTION) RETURN(GM_ERROR);

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
    if (SELECTIONOBJECT(theMG,i)==g) RETURN(GM_ERROR);

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
  COORD *x[MAX_CORNERS_OF_ELEM];
  COORD delta[3][DIM],s[DIM],t[DIM];
  COORD help,Scalarprdst,Scalarprd01,Scalarprd02,Scalarprd12;
  COORD_VECTOR theNormal[MAX_CORNERS_OF_ELEM];

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

  /* install the /Documents directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitUGManager","could not changedir to root");
    return(__LINE__);
  }
  theDocDirID = GetNewEnvDirID();
  if (MakeEnvItem("Documents",theDocDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitUGManager","could not install /Documents dir");
    return(__LINE__);
  }
  theMGVarID = GetNewEnvVarID();

  /* init the OBJT management */
  UsedOBJT = 0;
  for (i=0; i<NPREDEFOBJ; i++)
    SET_FLAG(UsedOBJT,1<<i);

  return (GM_OK);
}
