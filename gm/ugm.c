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
/*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*/
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

#include "devices.h"

#include "switch.h"
#include "evm.h"
#include "gm.h"
#include "misc.h"
#include "simplex.h"
#include "algebra.h"
#include "ugm.h"

#ifdef __TWODIM__
#include "shapes2d.h"
#endif

#ifdef __THREEDIM__
#include "shapes3d.h"
#include "ugm3d.h"
#endif

/* include refine because of macros accessed in list functions */
#include "ugrefine.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

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

/* data for CVS */
static char rcsid[] = "$Header$";

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
    return (GM_ERROR);

  for (theMG=GetFirstMultigrid(); theMG!=NULL; theMG=GetNextMultigrid(theMG))
    if (theMG->freeObjects[type]!=NULL)
    {
      PrintErrorMessage('E',"ReleaseOBJT","first clean freeObjects list to this type FOR ALL mgs");
      return (GM_ERROR);
    }

  CLEAR_FLAG(UsedOBJT,1<<type);

  return (GM_OK);
}


/****************************************************************************/
/*D
   GetFreeObject - Get an object from free list if possible

   SYNOPSIS:
   void *GetFreeObject (MULTIGRID *theMG, int n);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  n - type of the requested object

   DESCRIPTION:
   This function gets an object of type `n` from free list if possible.

   RETURN VALUE:
   void *
   .n   pointer to an object of the requested type
   .n   NULL if object of requested type is not available
   D*/
/****************************************************************************/

void *GetFreeObject (MULTIGRID *theMG, int n)
{
  void **ptr;

  if ((n<0)||(n>=MAXOBJECTS)) return(NULL);
  if (theMG->freeObjects[n]==NULL) return(NULL);

  /* 'ptr' will be set equal to 'theMG->freeObjects[n]' but with			*/
  /* different interpretation: void ** instead of void *. 'ptr'			*/
  /* points to the first two bytes of the object (i.e. unsigned INT ctrl	*/
  /* and INT id) but will be interpreted as a void * pointer, witch points*/
  /* to the next free object.                                                                                   */
  ptr = (void **) theMG->freeObjects[n];
  theMG->freeObjects[n] = ptr[0];
  return((void *)ptr);
}


/****************************************************************************/
/*D
   PutFreeObject - Put an object in the free list

   SYNOPSIS:
   INT PutFreeObject (MULTIGRID *theMG, void *object);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  object - object to insert in free list

   DESCRIPTION:
   This function puts an object in the free list.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 when error occured.
   D*/
/****************************************************************************/

INT PutFreeObject (MULTIGRID *theMG, void *object)
{
  void **ptr;
  int n;

  n = OBJT(object);
  if ((n<0)||(n>=MAXOBJECTS))
  {
    UserWrite("wrong object given to PutFreeObject\n");
    return(1);
  }

  /* 'ptr' will be set equal to 'object' but with different inter-		*/
  /* pretation: void ** instead of void *. 'ptr' points to the first		*/
  /* two bytes of the object (i.e. unsigned INT ctrl	and INT id) but         */
  /* will be interpreted as a void * pointer, witch will be set equal   */
  /* to 'theMG->freeObjects[n]' i.e. the first free object.				*/
  ptr = (void **) object;
  ptr[0] = theMG->freeObjects[n];
  theMG->freeObjects[n] = object;
  return(0);
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
  int i;

  vs = (VSEGMENT *)GetFreeObject(theGrid->mg,VSOBJ);
  if (vs==NULL)
  {
    vs = (VSEGMENT *)GetMem(theGrid->mg->theHeap,sizeof(struct vsegment),FROM_BOTTOM);
    if (vs==NULL) return (NULL);
  }

  /* initialize data */
  CTRL(vs) = 0;
  SETOBJT(vs,VSOBJ);
  BSEGDESC(vs) = NULL;
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
  int i;

  pv = GetFreeObject(theGrid->mg,BVOBJ);
  if (pv==NULL)
  {
    pv = GetMem(theGrid->mg->theHeap,sizeof(struct bvertex),FROM_BOTTOM);
    if (pv==NULL) return(NULL);
    if ((ds=theGrid->mg->theFormat->sVertex)>0)
    {
      VDATA(pv) = GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM);
      if (VDATA(pv)==NULL) return(NULL);
      memset(VDATA(pv),0,ds);
    }
    else
      VDATA(pv) = NULL;
  }

  /* initialize data */
  CTRL(pv) = 0;
  SETOBJT(pv,BVOBJ);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
  TOPNODE(pv) = NULL;
  for (i=0; i<DIM; i++) LCVECT(pv)[i] = 0.0;
  VSEG(pv) = NULL;
  SETONEDGE(pv,0);
  SETMOVE(pv,DIM_OF_BND);

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
  int i;

  pv = GetFreeObject(theGrid->mg,IVOBJ);
  if (pv==NULL)
  {
    pv = GetMem(theGrid->mg->theHeap,sizeof(struct ivertex),FROM_BOTTOM);
    if (pv==NULL) return(NULL);
    if ((ds=theGrid->mg->theFormat->sVertex)>0)
    {
      VDATA(pv) = GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM);
      if (VDATA(pv)==NULL) return(NULL);
      memset(VDATA(pv),0,ds);
    }
    else
      VDATA(pv) = NULL;
  }

  /* initialize data */
  CTRL(pv) = 0;
  SETOBJT(pv,IVOBJ);
  SETLEVEL(pv,theGrid->level);
  ID(pv) = (theGrid->mg->vertIdCounter)++;
  VFATHER(pv) = NULL;
  TOPNODE(pv) = NULL;
  SETMOVE(pv,DIM);
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
   .n   poiinter to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

NODE *CreateNode (GRID *theGrid, NODE *after)
{
  NODE *pn;

        #ifdef __version23__
  INT ds;
        #endif
        #ifdef __NODEDATA__
  VECTOR *pv;
        #endif

  pn = GetFreeObject(theGrid->mg,NDOBJ);
  if (pn==NULL)
  {
    pn = GetMem(theGrid->mg->theHeap,sizeof(NODE),FROM_BOTTOM);
    if (pn==NULL) return(NULL);
                #ifdef __version23__
    if ((ds=theGrid->mg->theFormat->sNode)>0)
    {
      NDATA(pn) = GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM);
      if (NDATA(pn)==NULL) return(NULL);
      memset(NDATA(pn),0,ds);
    }
    else
      NDATA(pn) = NULL;
    if ((ds=theGrid->mg->theFormat->sDiag)>0)
    {
      NDIAG(pn) = GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM);
      if (NDIAG(pn)==NULL) return(NULL);
      memset(NDIAG(pn),0,ds);
    }
    else
      NDIAG(pn) = NULL;
                #endif
  }

  /* create vector */
        #ifdef __NODEDATA__
  if (CreateVector (theGrid,NULL,NODEVECTOR,&pv))
  {
    DisposeNode (theGrid,pn);
    return (NULL);
  }
  assert (pv != NULL);
  VOBJECT(pv) = (void*)pn;
  NVECTOR(pn) = (void*)pv;
        #endif

  /* initialize data */
  CTRL(pn) = 0;
  SETOBJT(pn,NDOBJ);
  SETCLASS(pn,4);
  SETLEVEL(pn,theGrid->level);
  ID(pn) = (theGrid->mg->nodeIdCounter)++;
  INDEX(pn) = 0;
        #ifdef __version23__
  VSKIP(pn) = 0;
  NSKIP(pn) = 0;
        #endif
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
  EDGE *pe;
  LINK *pl;

  /* run through neighbor list */
  for (pl=START(from); pl!=NULL; pl = NEXT(pl))
    if (NBNODE(pl)==to)
    {
      pe = (EDGE *) (pl-LOFFSET(pl));
      return(pe);
    }

  /* return not found */
  return(NULL);
}


/****************************************************************************/
/*D
   CreateEdge - Return pointer to a new edge structure

   SYNOPSIS:
   EDGE *CreateEdge (GRID *theGrid, NODE *from, NODE *to);

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

EDGE *CreateEdge (GRID *theGrid, NODE *from, NODE *to)
{
  EDGE *pe;
  LINK *link0,*link1;

        #ifdef __version23__
  INT ds;
        #endif
        #ifdef __EDGEDATA__
  VECTOR *pv;
        #endif

  /* check if edge exists already */
  if( (pe = GetEdge(from, to)) != NULL ) return(pe);

  pe = GetFreeObject(theGrid->mg,EDOBJ);
  if (pe==NULL)
  {
    pe = GetMem(theGrid->mg->theHeap,sizeof(EDGE),FROM_BOTTOM);
    if (pe==NULL) return(NULL);
                #ifdef __version23__
    if ((ds=theGrid->mg->theFormat->sLink)>0)
    {
      LDATA(LINK0(pe)) = GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM);
      if (LDATA(LINK0(pe))==NULL) return(NULL);
      memset(LDATA(LINK0(pe)),0,ds);
      LDATA(LINK1(pe)) = GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM);
      if (LDATA(LINK1(pe))==NULL) return(NULL);
      memset(LDATA(LINK1(pe)),0,ds);
    }
    else
    {
      LDATA(LINK0(pe)) = NULL;
      LDATA(LINK1(pe)) = NULL;
    }
    if ((ds=theGrid->mg->theFormat->sEdge)>0)
    {
      EDDATA(pe) = GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM);
      if (EDDATA(pe)==NULL) return(NULL);
      memset(EDDATA(pe),0,ds);
    }
    else
      EDDATA(pe) = NULL;
                #endif
  }

  /* initialize data */
  link0 = LINK0(pe);
  link1 = LINK1(pe);
  CTRL(pe) = 0;
  CTRL(link1) = 0;
  SETOBJT(pe,EDOBJ);
  SETLOFFSET(link0,0);
  SETLOFFSET(link1,1);
  SETLEVEL(pe,theGrid->level);
  NBNODE(link0) = to;
  NBNODE(link1) = from;
        #ifdef __THREEDIM__
  SET_NO_OF_ELEM(pe,0);
  SETTAG(pe,0);
  SETEDGENEW(pe,1);
        #endif
        #ifdef __MIDNODE__
  MIDNODE(pe)=NULL;
        #endif

  /* create vector if */
        #ifdef __EDGEDATA__
  if (CreateVector (theGrid,NULL,EDGEVECTOR,&pv))
  {
    DisposeEdge (theGrid,pe);
    return (NULL);
  }
  assert (pv != NULL);
  VOBJECT(pv) = (void*)pe;
  EDVECTOR(pe) = (void*)pv;
        #endif

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
   CreateBoundaryElement - Allocate a new element structure

   SYNOPSIS:
   ELEMENT *CreateBoundaryElement (GRID *theGrid, ELEMENT *after, INT tag);

   PARAMETERS:
   .  theGrid - grid structure to extend
   .  after - insert after that element (NULL if b.o.l.)
   .  tag - the element type

   DESCRIPTION:
   This function creates and initializes a new element of type `tag` and
   returns a pointer to it.

   RETURN VALUE:
   ELEMENT *
   .n   pointer to requested element
   .n   NULL if out of memory
   D*/
/****************************************************************************/

ELEMENT *CreateBoundaryElement (GRID *theGrid, ELEMENT *after, INT tag)
{
  ELEMENT *pe;
  INT i;

        #if defined __ELEMDATA__ || defined __SIDEDATA__
  VECTOR *pv;
        #endif
        #ifdef __version23__
  INT ds;
        #endif

  pe = GetFreeObject(theGrid->mg,MAPPED_BND_OBJT(tag));
  if (pe==NULL)
  {
    pe = GetMem(theGrid->mg->theHeap,BND_SIZE(tag),FROM_BOTTOM);
    if (pe==NULL) return(NULL);
                #ifdef __version23__
    if ((ds=theGrid->mg->theFormat->sElement)>0)
    {
      SET_EDATA(pe,GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM));
      if (EDATA(pe)==NULL) return(NULL);
      memset(EDATA(pe),0,ds);
    }
    else
      SET_EDATA(pe,NULL);
                #endif
  }

  /* initialize data */
  CTRL(pe) = 0;
  CTRL2(pe) = 0;
  SETOBJT(pe,BEOBJ);
  SETTAG(pe,tag);
  SETLEVEL(pe,theGrid->level);
        #ifdef __version3__
  SETEBUILDCON(pe,1);
        #endif
        #ifdef __THREEDIM__
  SETNEWEL(pe,TRUE);
        #endif
  ID(pe) = (theGrid->mg->elemIdCounter)++;
  SET_EFATHER(pe,NULL);
  for (i=0; i<CORNERS_OF_ELEM(pe); i++)
  {
    SET_CORNER(pe,i,NULL);
  }
  for (i=0; i<SIDES_OF_ELEM(pe); i++)
  {
    SET_NBELEM(pe,i,NULL);
    SET_SIDE(pe,i,NULL);
  }
        #ifdef __TWODIM__
  for (i=0; i<SONS_OF_ELEM(pe); i++) SET_SON(pe,i,NULL);
        #endif
        #ifdef __THREEDIM__
  SET_SON(pe,0,NULL);
        #endif

  /* create element vector if */
        #ifdef __ELEMDATA__
  if (CreateVector (theGrid,NULL,ELEMVECTOR,&pv))
  {
    DisposeElement (theGrid,pe);
    return (NULL);
  }
  assert (pv != NULL);
  VOBJECT(pv) = (void*)pe;
  SET_EVECTOR(pe,(void*)pv);
        #endif

  /* create side vectors if */
        #ifdef __SIDEDATA__
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
   CreateInnerElement - Return a pointer to  a new element structure

   SYNOPSIS:
   CreateInnerElement

   PARAMETERS:
   .  theGrid - grid structure to extend
   .  after - insert after that element (NULL if b.o.l.)
   .  tag - the element type

   DESCRIPTION:
   This function creates and initializes a new element and returns a pointer to it.

   RETURN VALUE:
   ELEMENT *
   .n   pointer to requested object
   .n   NULL if out of memory
   D*/
/****************************************************************************/

ELEMENT *CreateInnerElement (GRID *theGrid, ELEMENT *after, INT tag)
{
  ELEMENT *pe;
  INT i;
        #if defined __ELEMDATA__ || defined __SIDEDATA__
  VECTOR *pv;
        #endif
        #ifdef __version23__
  INT ds;
        #endif

  pe = GetFreeObject(theGrid->mg,MAPPED_INNER_OBJT(tag));
  if (pe==NULL)
  {
    pe = GetMem(theGrid->mg->theHeap,INNER_SIZE(tag),FROM_BOTTOM);
    if (pe==NULL) return(NULL);
                #ifdef __version23__
    if ((ds=theGrid->mg->theFormat->sElement)>0)
    {
      SET_EDATA(pe,GetMem(theGrid->mg->theHeap,ds,FROM_BOTTOM));
      if (EDATA(pe)==NULL) return(NULL);
      memset(EDATA(pe),0,ds);
    }
    else
      SET_EDATA(pe,NULL);
                #endif
  }

  /* initialize data */
  CTRL(pe) = 0;
  CTRL2(pe) = 0;
  SETOBJT(pe,IEOBJ);
  SETTAG(pe,tag);
  SETLEVEL(pe,theGrid->level);
        #ifdef __version3__
  SETEBUILDCON(pe,1);
        #endif
        #ifdef __THREEDIM__
  SETNEWEL(pe,TRUE);
        #endif
  ID(pe) = (theGrid->mg->elemIdCounter)++;
  SET_EFATHER(pe,NULL);
  for (i=0; i<CORNERS_OF_ELEM(pe); i++)
  {
    SET_CORNER(pe,i,NULL);
  }
  for (i=0; i<SIDES_OF_ELEM(pe); i++)
  {
    SET_NBELEM(pe,i,NULL);
  }
        #ifdef __TWODIM__
  for (i=0; i<SONS_OF_ELEM(pe); i++) SET_SON(pe,i,NULL);
        #endif
        #ifdef __THREEDIM__
  SET_SON(pe,0,NULL);
        #endif

  /* create element vector if */
        #ifdef __ELEMDATA__
  if (CreateVector (theGrid,NULL,ELEMVECTOR,&pv))
  {
    DisposeElement (theGrid,pe);
    return (NULL);
  }
  assert (pv != NULL);
  VOBJECT(pv) = (void*)pe;
  SET_EVECTOR(pe,(void*)pv);
        #endif

  /* create side vectors if */
        #ifdef __SIDEDATA__
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

  ps = GetFreeObject(theGrid->mg,ESOBJ);
  if (ps==NULL)
  {
    ps = GetMem(theGrid->mg->theHeap,sizeof(ELEMENTSIDE),FROM_BOTTOM);
    if (ps==NULL) return(NULL);
  }

  /* initialize data */
  CTRL(ps) = 0;
  SETOBJT(ps,ESOBJ);
  SEGDESC(ps) = NULL;

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
  theGrid = GetFreeObject(theMG, GROBJ);
  if (theGrid==NULL)
  {
    theGrid = GetMem(theMG->theHeap,sizeof(GRID),FROM_BOTTOM);
    if (theGrid==NULL) return(NULL);
  }
  memset(theGrid,0,sizeof(GRID));

  /* fill in data */
  CTRL(theGrid) = 0;
  SETOBJT(theGrid,GROBJ);
  theGrid->level = l;
  theGrid->nVert = 0;
  theGrid->nNode = 0;
  theGrid->nEdge = 0;
  theGrid->nElem = 0;
  theGrid->nSide = 0;
  theGrid->nVector = 0;
  theGrid->nCon = 0;
  theGrid->status = 0;
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

  theDocDir = ChangeEnvDir("/Documents");

  assert (theDocDir!=NULL);

  return ((MULTIGRID *) ENVDIR_DOWN(theDocDir));
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
  return ((MULTIGRID *) NEXT_ENVITEM(theMG));
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

MULTIGRID *CreateMultiGrid (char *MultigridName, char *domain, char *problem, char *format, unsigned long heapSize)
{
  HEAP *theHeap,*theUserHeap;
  MULTIGRID *theMG;
  GRID *theGrid;
  BOUNDARY_SEGMENT *theSegment;
  BOUNDARY_CONDITION *theBndCond;
  VERTEX **pv;
  NODE *pn;
  int i,j,k,n,ds,l,FatalError;
  COORD cvect[DIM];
  DOUBLE pardist;
  VSEGMENT *vs;
  DOMAIN *theDomain;
  PROBLEM *theProblem;
  FORMAT *theFormat;

  theDomain = GetDomain(domain);
  if (theDomain==NULL)
  {
    PrintErrorMessage('E',"CreateMultiGrid","domain not found");
    return(NULL);
  }

  theProblem = GetProblem(domain,problem);
  if (theProblem==NULL)
  {
    PrintErrorMessage('E',"CreateMultiGrid","problem not found");
    return(NULL);
  }

  theFormat = GetFormat(format);
  if (theFormat==NULL)
  {
    PrintErrorMessage('E',"CreateMultiGrid","format not found");
    return(NULL);
  }

  /* allocate multigrid envitem */
  theMG = MakeMGItem(MultigridName);
  if (theMG==NULL) return(NULL);

  /* allocate the heap */
  theHeap = NewHeap(SIMPLE_HEAP, heapSize, malloc(heapSize));
  if (theHeap==NULL)
  {
    DisposeMultiGrid(theMG);
    return(NULL);
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
  theMG->theDomain = theDomain;
  theMG->theFormat = theFormat;
  theMG->theProblem = theProblem;
  theMG->theHeap = theHeap;
  SELECTIONSIZE(theMG) = 0;
  for (i=0; i<MAXLEVEL; i++) theMG->grids[i] = NULL;
  for (i=0; i<MAXOBJECTS; i++) theMG->freeObjects[i] = NULL;
  for (i=0; i<MAXVECTORS; i++) theMG->freeVectors[i] = NULL;
  for (i=0; i<MAXCONNECTIONS; i++) theMG->freeConnections[i] = NULL;

  /* allocate level 0 grid */
  theGrid = CreateNewLevel(theMG);
  if (theGrid==NULL)
  {
    DisposeMultiGrid(theMG);
    return(NULL);
  }

  /* allocate boundary descriptors */
  n = theDomain->numOfSegments;
  theMG->numOfSegments = n;
  theMG->segments = (BNDSEGDESC *) GetMem(theHeap,n*sizeof(BNDSEGDESC),FROM_BOTTOM);
  if (theMG->segments==NULL) { DisposeMultiGrid(theMG); return(NULL); }
  for (i=0; i<n; i++)
  {
    theMG->segments[i].theSegment = NULL;
    theMG->segments[i].theBoundaryCondition = NULL;
  }

  /* combine boundary coordinates and boundary conditions */
  for (theSegment=GetFirstBoundarySegment(theDomain); theSegment!=NULL; theSegment = GetNextBoundarySegment(theSegment))
  {
    i = theSegment->id;
    if ((i<0)||(i>=n))
    {
      PrintErrorMessage('E',"CreateMultiGrid","segment id out of range");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    if (theMG->segments[i].theSegment!=NULL)
    {
      PrintErrorMessage('E',"CreateMultiGrid","segment id multiply defined");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    theMG->segments[i].theSegment = theSegment;
    theMG->segments[i].id = i;
  }
  for (theBndCond=GetFirstBoundaryCondition(theProblem); theBndCond!=NULL; theBndCond = GetNextBoundaryCondition(theBndCond))
  {
    i = theBndCond->id;
    if ((i<0)||(i>=n))
    {
      PrintErrorMessage('E',"CreateMultiGrid","boundary condition id out of range");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    if (theMG->segments[i].theBoundaryCondition!=NULL)
    {
      PrintErrorMessage('E',"CreateMultiGrid","boundary condition id multiply defined");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    theMG->segments[i].theBoundaryCondition = theBndCond;
  }

  /* check if all pointers are correctly defined */
  for (i=0; i<n; i++)
  {
    if (theMG->segments[i].theSegment==NULL)
    {
      PrintErrorMessage('E',"CreateMultiGrid","boundary segment not found");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    if (theMG->segments[i].theBoundaryCondition==NULL)
    {
      PrintErrorMessage('E',"CreateMultiGrid","boundary conditionnot found");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
  }

  /* allocate corner vertices pointers */
  n = theDomain->numOfCorners;
  theMG->numOfCorners = n;
  pv = (VERTEX **) GetMem(theHeap,n*sizeof(VERTEX *),FROM_BOTTOM);
  if (pv==NULL) { DisposeMultiGrid(theMG); return(NULL); }
  theMG->corners = pv;

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
  }

  /* create and fill segment data of vertex */
  for (i=0; i<theMG->numOfSegments; i++)
  {
    theSegment = theMG->segments[i].theSegment;

    for( k=0; k<CORNERS_OF_BND_SEG; k++ )
    {
      j = theSegment->points[k];
      SETUSED(pv[j],1);
      vs = CreateVertexSegment(theGrid,pv[j]);
      if (vs==NULL) { DisposeMultiGrid(theMG); return(NULL); }

      BSEGDESC(vs) = &(theMG->segments[i]);

      FatalError = 0;
      if (DIM_OF_BND == 1)
        switch(k)
        {
        case 0 : LAMBDA(vs,0) = theSegment->alpha[0];   break;
        case 1 : LAMBDA(vs,0) = theSegment->beta[0];    break;
        default : FatalError = 1;                                                break;
        }
      else if (DIM_OF_BND == 2)
        switch(k)
        {
        case 0 : LAMBDA(vs,0) = theSegment->alpha[0];
          LAMBDA(vs,1) = theSegment->alpha[1];   break;
        case 1 : LAMBDA(vs,0) = theSegment->beta[0];
          LAMBDA(vs,1) = theSegment->alpha[1];   break;
        case 2 : LAMBDA(vs,0) = theSegment->beta[0];
          LAMBDA(vs,1) = theSegment->beta[1];    break;
        case 3 : LAMBDA(vs,0) = theSegment->alpha[0];
          LAMBDA(vs,1) = theSegment->beta[1];    break;
        default : FatalError = 1;                                                break;
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
        (*theSegment->BndSegFunc)(theSegment->data,PVECT(vs),CVECT(pv[j]));
        continue;                          /* next point (k) */
      }

      /* if the vs is not the first one compare the new geometrical data with the old ones */
      (*theSegment->BndSegFunc)(theSegment->data,PVECT(vs),cvect);
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

  /* return ok */
  return(theMG);
}


/****************************************************************************/
/*D
   DisposeEdge - Remove edge from the data structure

   SYNOPSIS:
   INT DisposeEdge (GRID *theGrid, EDGE *theEdge);

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

INT DisposeEdge (GRID *theGrid, EDGE *theEdge)
{
  LINK *link0,*link1,*pl;
  NODE *from,*to;
  int found;

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

        #ifdef __EDGEDATA__
  /* dispose vector and its matrices from edge-vector if */
  if (DisposeVector (theGrid,EDVECTOR(theEdge)))
    return(1);
        #endif

  /* free edge object */
  if (PutFreeObject(theGrid->mg,theEdge)>0) return(1);

  /* check error condition */
  if (found!=2) return(1);

  /* return ok */
  theGrid->nEdge--;
  return(0);
}


/****************************************************************************/
/*D
   DisposeNode - Remove node including its edges from the data structure

   SYNOPSIS:
   INT DisposeNode (GRID *theGrid, NODE *theNode);

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

INT DisposeNode (GRID *theGrid, NODE *theNode)
{
  LINK *link0,*link1,*pl;
  EDGE *pe;
  NODE *to;
  int found;

  /* remove links in all neighbors lists */
  found = 0;
  for (link0=START(theNode); link0!=NULL; link0=NEXT(link0))
  {
    found--;
    link1 = REVERSE(link0);
    to = NBNODE(link0);
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
  }

  /* now delete the edges */
  pl = START(theNode);
  while (pl!=NULL)
  {
    pe = MYEDGE(pl);
    pl = NEXT(pl);
                #ifdef __EDGEDATA__
    /* dispose vector and its matrices from edge-vector if */
    if (DisposeVector (theGrid,EDVECTOR(pe)))
      return(1);
                #endif
    if (PutFreeObject(theGrid->mg,pe)==0) theGrid->nEdge--;
  }

  /* remove node from node list */
  if (PREDN(theNode)!=NULL)
    SUCCN(PREDN(theNode)) = SUCCN(theNode);
  else
    theGrid->firstNode = SUCCN(theNode);
  if (SUCCN(theNode)!=NULL)
    PREDN(SUCCN(theNode)) = PREDN(theNode);
  else
    theGrid->lastNode = PREDN(theNode);

        #ifdef __NODEDATA__
  /* dispose vector and its matrices from node-vector */
  if (DisposeVector (theGrid,NVECTOR(theNode)))
    return(1);
        #endif

  /* delete the node itself */
  if (PutFreeObject(theGrid->mg,theNode)>0) return(1);

  /* check error condition */
  if (found!=0) return(1);

  /* return ok */
  (theGrid->nNode)--;
  return(0);
}


/****************************************************************************/
/*D
   DisposeVertex - Remove vertex from the data structure

   SYNOPSIS:
   INT DisposeVertex (GRID *theGrid, VERTEX *theVertex);

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

INT DisposeVertex (GRID *theGrid, VERTEX *theVertex)
{
  VSEGMENT *theVSeg;

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
      PutFreeObject(theGrid->mg,theVSeg);
      theVSeg = NEXTSEG(theVSeg);
    }
  }

  if (PutFreeObject(theGrid->mg,theVertex)>0) return(1);

  theGrid->nVert--;
  return(0);
}


/****************************************************************************/
/*D
   DisposeElementSide - Put element side in free list

   SYNOPSIS:
   INT DisposeElementSide  (GRID *theGrid, ELEMENTSIDE *theElementSide);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  theElementSide - element side to remove

   DESCRIPTION:
   This function puts an element side into the free list.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 no valid object number
   D*/
/****************************************************************************/


INT DisposeElementSide  (GRID *theGrid, ELEMENTSIDE *theElementSide)
{
  /* remove elementside from elementside list */
  if (PREDS(theElementSide)!=NULL)
    SUCCS(PREDS(theElementSide)) = SUCCS(theElementSide);
  else
    theGrid->sides = SUCCS(theElementSide);
  if (SUCCS(theElementSide)!=NULL)
    PREDS(SUCCS(theElementSide)) = PREDS(theElementSide);

  if (PutFreeObject(theGrid->mg,theElementSide)>0) return(1);

  theGrid->nSide--;
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
  int i;
        #ifdef __SIDEDATA__
  int j;
  VECTOR *theVector;
  ELEMENT *theNeighbor;
        #endif

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
      if (SIDE(theElement,i)!=NULL) DisposeElementSide(theGrid,SIDE(theElement,i));

  /* dispose matrices from element-vector */
        #ifdef __version3__
  if (DisposeConnectionFromElement(theGrid,theElement))
    return(1);
        #endif

  /* dispose vectors in sides if */
        #ifdef __SIDEDATA__
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    theVector = SVECTOR(theElement,i);
    assert (VCOUNT(theVector) != 0);
    assert (VCOUNT(theVector) != 3);
    if (VCOUNT(theVector) == 1)
    {
      if (DisposeVector (theGrid,theVector))
        return (1);
    }
    else
    {
      if (!FindNeighborElement (theElement,i,&theNeighbor,&j))
        return (1);
      VOBJECT(theVector) = (void*)theNeighbor;
      SETVECTORSIDE(theVector,j);
      SETVCOUNT(SVECTOR(theElement,i),1);
    }
  }
        #endif

  /* dispose vector in center of element */
        #ifdef __ELEMDATA__
  if (DisposeVector (theGrid,EVECTOR(theElement)))
    return(1);
        #endif


  /* dispose element */
  /* give it a new tag ! (I know this is somewhat ugly) */
  if (OBJT(theElement)==BEOBJ)
    SETOBJT(theElement,MAPPED_BND_OBJT(TAG(theElement)));
  else
    SETOBJT(theElement,MAPPED_INNER_OBJT(TAG(theElement)));
  if (PutFreeObject(theGrid->mg,theElement)>0) return(1);

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
  if (FIRSTELEMENT(theGrid)!=NULL) return(2);
  if (FIRSTVERTEX(theGrid)!=NULL) return(2);
  if (FIRSTNODE(theGrid)!=NULL) return(2);

  /* remove from grids array */
  theMG->grids[l] = NULL;
  theMG->grids[l-1]->finer = NULL;
  (theMG->topLevel)--;
  if (theMG->currentLevel>theMG->topLevel) theMG->currentLevel = theMG->topLevel;

  if (PutFreeObject(theGrid->mg,theGrid)>0) return(1);

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
  if (ChangeEnvDir("/Documents")==NULL) return (GM_ERROR);
  if (RemoveEnvItem ((ENVITEM *)theMG)) return (GM_ERROR);

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

static int LexCompare (NODE **pnode1, NODE **pnode2)
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
  DOMAIN *theDomain;
  INT i,entries,firstID,nl;
  HEAP *theHeap;

  theMG   = MYMG(theGrid);
  firstID = ID(FIRSTNODE(theGrid));
  entries = NN(theGrid);

  /* calculate the diameter of the bounding rectangle of the domain */
  theDomain = MGDOMAIN(theMG);
  if (theDomain==NULL) return (1);
  InvMeshSize = POW2(GLEVEL(theGrid)) * pow(NN(GRID_ON_LEVEL(theMG,0)),1.0/DIM) / theDomain->radius;

  /* allocate memory for the node list */
  theHeap = MGHEAP(theMG);
  Mark(theHeap,FROM_TOP);
  if ((table=GetMem(theHeap,entries*sizeof(NODE *),FROM_TOP))==NULL)
  {
    Release(theHeap,FROM_TOP);
    PrintErrorMessage('E',"OrderNodesInGrid","ERROR: could not allocate memory from the MGHeap");
    return (2);
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
        return (1);

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
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* create objects */
  theVertex = CreateInnerVertex(theGrid,NULL);
  if (theVertex==NULL)
  {
    PrintErrorMessage('E',"InsertInnerNode","cannot create vertex");
    return(GM_ERROR);
  }
  theNode = CreateNode(theGrid,NULL);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertInnerNode","cannot create node");
    return(GM_ERROR);
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

INT InsertBoundaryNode (MULTIGRID *theMG, INT bnd_seg_id, COORD *pos)
{
  GRID *theGrid;
  NODE *theNode;
  VERTEX *theVertex;
  VSEGMENT *vsnew1;
  BNDSEGDESC *theSegment1;
  int i;


  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertBoundaryNode","only a multigrid with exactly one level can be edited");
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* scan input */
  if ((bnd_seg_id<0)||(bnd_seg_id>=MGNOOFSEG(theMG)))
  {
    PrintErrorMessage('E',"InsertBoundaryNode","segment id out of range");
    return(GM_ERROR);
  }
  theSegment1 = MGBNDSEGDESC(theMG,bnd_seg_id);

  /* check range of parameters */
  for(i=0; i<DIM_OF_BND; i++)
    if ( (pos[0]<ALPHA(theSegment1,i))||(pos[0]> BETA(theSegment1,i)) )
    {
      PrintErrorMessage('E',"InsertBoundaryNode","parameter not in range of segment");
      return(GM_ERROR);
    }

  /* check distance from corner */
  if (DIM==2)
  {
    if (  (fabs(pos[0]-ALPHA(theSegment1,0))<SMALL_C) || (fabs(pos[0]-BETA(theSegment1,0))<SMALL_C) )
    {
      PrintErrorMessage('E',"InsertBoundaryNode","parameter describes one of the corners of the segment");
      return(GM_ERROR);
    }
  }
  if (DIM==3)
  {
    i = 0;
    if (fabs(pos[0]-ALPHA(theSegment1,0))+fabs(pos[1]-ALPHA(theSegment1,1))<SMALL_C) i = 1;
    if (fabs(pos[0]-ALPHA(theSegment1,0))+fabs(pos[1]-BETA(theSegment1,1))<SMALL_C) i = 1;
    if (fabs(pos[0]-BETA(theSegment1,0))+fabs(pos[1]-ALPHA(theSegment1,1))<SMALL_C) i = 1;
    if (fabs(pos[0]-BETA(theSegment1,0))+fabs(pos[1]-BETA(theSegment1,1))<SMALL_C) i = 1;
    if (i)
    {
      PrintErrorMessage('E',"InsertBoundaryNode","parameter describes one of the corners of the segment");
      return(GM_ERROR);
    }
  }

  /* create objects */
  theVertex = CreateBoundaryVertex(theGrid,NULL);
  if (theVertex==NULL)
  {
    PrintErrorMessage('E',"InsertBoundaryNode","cannot create vertex");
    return(GM_ERROR);
  }
  theNode = CreateNode(theGrid,NULL);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    PrintErrorMessage('E',"InsertBoundaryNode","cannot create node");
    return(GM_ERROR);
  }
  vsnew1 = CreateVertexSegment(theGrid, theVertex);
  if (vsnew1==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    DisposeNode(theGrid, theNode);
    PrintErrorMessage('E',"InsertBoundaryNode","cannot create vertexsegment");
    return(GM_ERROR);
  }

  /* fill data into the first vertexsegment */
  for(i=0; i<DIM_OF_BND; i++) LAMBDA(vsnew1,i) =  pos[i];
  (*BNDSEGFUNC (theSegment1))(BNDDATA(theSegment1),pos,CVECT(theVertex));
  BSEGDESC(vsnew1) = theSegment1;

  /* fill data into node/vertex */
  INDEX(theNode) = 0;
  MYVERTEX(theNode) = theVertex;
  SETMOVE(theVertex,DIM_OF_BND);

  return(GM_OK);
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
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  if (theNode==NULL)
  {
    PrintErrorMessage('E',"DeleteNode","node not found");
    return(GM_ERROR);
  }

  /* check corner */
  theVertex = MYVERTEX(theNode);
  if (MOVE(theVertex)==0)
  {
    PrintErrorMessage('E',"DeleteNode","corners cannot be deleted");
    return(GM_ERROR);
  }

  /* check if some element needs that node */
  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      if (CORNER(theElement,i)==theNode)
      {
        PrintErrorMessage('E',"DeleteNode","there is an element needing that node");
        return(GM_ERROR);
      }

  /* now allowed to delete */
  DisposeVertex(theGrid,theVertex);
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
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* find node */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    if (ID(theNode)==id) break;
  if (theNode==NULL)
  {
    PrintErrorMessage('E',"DeleteNodeWithID","node not found");
    return(GM_ERROR);
  }
  return (DeleteNode(theMG,theNode));
}


/****************************************************************************/
/*D
   InsertElement - Insert an element

   SYNOPSIS:
   InsertElement (MULTIGRID *theMG, INT n, NODE *Node[4]); // 2D Version
   INT InsertElement (MULTIGRID *theMG, INT n,
   NODE *Node[MAX_CORNERS_OF_ELEM]); // 3D Version

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

INT InsertElement (MULTIGRID *theMG, INT n, NODE *Node[4]) /* 2D VERSION */
{
  GRID *theGrid;
  int i,j,m,found,NeighborSide[4],bndEdges;
  BNDSEGDESC *theSeg[4];
  NODE *theNode,*aNode,*bNode;
  VERTEX *Vertex[4],*theVertex,*aVertex,*bVertex;
  VSEGMENT *aVSeg, *bVSeg, *aUniqueSeg, *bUniqueSeg;
  ELEMENT *theElement,*Neighbor[4];
  COORD from[4],to[4];
  EDGE *theEdge;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertElement","only a multigrid with exactly one level can be edited");
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* check tag */
  if ((n<3)&&(n>4))
  {
    PrintErrorMessage('E',"InsertElement","only triangles and quadrilaterals allowed in 2D");
    return(GM_ERROR);
  }

  /* init data */
  for (i=0; i<n; i++)
  {
    Vertex[i]       = MYVERTEX(Node[i]);
    Neighbor[i] = NULL;
    theSeg[i]       = NULL;
  }

  /* find orientation */
  if (!CheckOrientation(n,Vertex))
  {
    /* flip order */
    for (i=0; i<n/2; i++)
    {
      j = n-i-1;
      theNode = Node[i]; Node[i] = Node[j]; Node[j] = theNode;
      theVertex = Vertex[i]; Vertex[i] = Vertex[j]; Vertex[j] = theVertex;
    }

    if (!CheckOrientation(n,Vertex))
    {
      PrintErrorMessage('E',"InsertElement","no convex polygon");
      return(GM_ERROR);
    }
  }

  /* compute side information (theSeg[i]==NULL) means inner side */
  bndEdges = 0;
  for (i=0; i<n; i++)
  {
    aNode = Node[i]; bNode = Node[(i+1)%n];
    aVertex = Vertex[i]; bVertex = Vertex[(i+1)%n];

    /* if at least one vertex is in the interior then skip */
    if ((OBJT(aVertex)==IVOBJ)||(OBJT(bVertex)==IVOBJ)) continue;

    /* aVertex and bVertex must be on the same boundary segment with unique parameters */
    found = 0;
    for (aVSeg = VSEG(aVertex); aVSeg != NULL; aVSeg = NEXTSEG(aVSeg))
      for (bVSeg = VSEG(bVertex); bVSeg != NULL; bVSeg = NEXTSEG(bVSeg))
        if( BSEGDESC(aVSeg) == BSEGDESC(bVSeg) )
        {
          found++;
          aUniqueSeg = aVSeg;
          bUniqueSeg = bVSeg;
        }

    /* no common boundary segment -> this will be an interior edge */
    if (found==0) continue;

    /* common boundary segment, but non unique parameters -> refuse to insert (no closed boundary segments!) */
    if (found>1)
    {
      PrintErrorMessage('E',"InsertElement","non unique parameter encountered, do not define closed boundary segments");
      return(GM_ERROR);
    }

    /* both nodes are on same boundary segment with unique parameters */
    theSeg[i] = BSEGDESC(aUniqueSeg);
    from[i] = LAMBDA(aUniqueSeg,0); to[i] = LAMBDA(bUniqueSeg,0);
    if (from[i]<to[i])
    {
      if (LEFT(theSeg[i])<=0)
      {
        PrintErrorMessage('E',"InsertElement","element outside of domain");
        return(GM_ERROR);
      }
    }
    else
    {
      if (RIGHT(theSeg[i])<=0)
      {
        PrintErrorMessage('E',"InsertElement","element outside of domain");
        return(GM_ERROR);
      }
    }
    bndEdges++;
  }

  /* find neighboring elements */
  for (i=0; i<n; i++)
  {
    aNode = Node[i]; bNode = Node[(i+1)%n];
    found = 0;
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      m = CORNERS_OF_ELEM(theElement);
      for (j=0; j<m; j++)
        if ((CORNER(theElement,j)==bNode)&&(CORNER(theElement,(j+1)%m)==aNode))
        {
          if (NBELEM(theElement,j)!=NULL)
          {
            PrintErrorMessage('E',"InsertElement","neighbor relation inconsistent");
            return(GM_ERROR);
          }
          found++;
          Neighbor[i] = theElement;
          NeighborSide[i] = j;
        }
    }
    if (found>1)
    {
      PrintErrorMessage('E',"InsertElement","ooops, found more than one neighbor");
      return(GM_ERROR);
    }
  }

  /* create element */
  if (bndEdges>0)
    theElement = CreateBoundaryElement(theGrid,NULL,n);
  else
    theElement = CreateInnerElement(theGrid,NULL,n);
  if (theElement==NULL)
  {
    PrintErrorMessage('E',"InsertElement","could not create element");
    return(GM_ERROR);
  }

  /* create element sides if necessary */
  if (bndEdges>0)
  {
    for (i=0; i<n; i++)
    {
      SET_SIDE(theElement,i,NULL);
      if (theSeg[i]!=NULL)
      {
        SET_SIDE(theElement,i,CreateElementSide(theGrid));
        if (SIDE(theElement,i)==NULL)
        {
          PrintErrorMessage('E',"InsertElement","could not create element side");
          DisposeElement(theGrid,theElement);
          return(GM_ERROR);
        }
        SEGDESC(SIDE(theElement,i)) = theSeg[i];
        PARAM(SIDE(theElement,i),0,0) = from[i];
        PARAM(SIDE(theElement,i),1,0) = to[i];
      }
    }
  }

  /* create edges */
  for (i=0; i<n; i++)
  {
    aNode = Node[i]; bNode = Node[(i+1)%n];
    theEdge = GetEdge(aNode,bNode);
    if (theEdge==NULL)
    {
      theEdge = CreateEdge(theGrid,aNode,bNode);
      if (theEdge==NULL)
      {
        PrintErrorMessage('E',"InsertElement","could not create edge");
        DisposeElement(theGrid,theElement);
        return(GM_ERROR);
      }
    }
  }

  /* create diagonals with extra flag */
        #ifdef __version23__
  if (n==4)
  {
    aNode = Node[0]; bNode = Node[2];
    theEdge = GetEdge(aNode,bNode);
    if (theEdge==NULL)
    {
      theEdge = CreateEdge(theGrid,aNode,bNode);
      if (theEdge==NULL)
      {
        PrintErrorMessage('E',"InsertElement","could not create edge");
        DisposeElement(theGrid,theElement);
        return(GM_ERROR);
      }
      SETEXTRA(theEdge,1);
    }
    aNode = Node[1]; bNode = Node[3];
    theEdge = GetEdge(aNode,bNode);
    if (theEdge==NULL)
    {
      theEdge = CreateEdge(theGrid,aNode,bNode);
      if (theEdge==NULL)
      {
        PrintErrorMessage('E',"InsertElement","could not create edge");
        DisposeElement(theGrid,theElement);
        return(GM_ERROR);
      }
      SETEXTRA(theEdge,1);
    }
  }
        #endif

  /* fill element data */
  SETECLASS(theElement,RED);       /* this is a coarse grid element */
  SETTAG(theElement,n);
  SET_EFATHER(theElement,NULL);
  for (i=0; i<n; i++)
  {
    SET_CORNER(theElement,i,Node[i]);
    SET_NBELEM(theElement,i,Neighbor[i]);
    if (Neighbor[i]!=NULL)
      SET_NBELEM(Neighbor[i],NeighborSide[i],theElement);
  }
  for (i=0; i<SONS_OF_ELEM(theElement); i++) SET_SON(theElement,i,NULL);


        #ifdef __version3__
  /* create algebra connections */
  if (InsertedElementCreateConnection(theGrid,theElement))
  {
    PrintErrorMessage('E',"InsertElement","could not create algebra connections");
    DisposeElement(theGrid,theElement);
    return(GM_ERROR);
  }
        #endif

  return(GM_OK);
}
#endif


#ifdef __THREEDIM__

INT InsertElement (MULTIGRID *theMG, INT n, NODE *Node[MAX_CORNERS_OF_ELEM]) /* 3D VERSION */
{
  GRID             *theGrid;
  int i,j,k,l,m,found,num;
  int NeighborSide[MAX_SIDES_OF_ELEM];
  BNDSEGDESC       *theSeg[MAX_SIDES_OF_ELEM], *Seg;
  NODE             *sideNode[MAX_CORNERS_OF_SIDE];
  VERTEX           *Vertex[MAX_CORNERS_OF_ELEM],*sideVertex[MAX_CORNERS_OF_SIDE];
  ELEMENT          *theElement,*Neighbor[MAX_SIDES_OF_ELEM];
  EDGE             *theEdge;
  COORD param[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE][DIM_OF_BND], *plambda[MAX_CORNERS_OF_SIDE];
  VSEGMENT         *vs;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"InsertElement","only a multigrid with exactly one level can be edited");
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* check parameters */
  if ( n != 4 )
  {
    PrintErrorMessage('E',"InsertElement","exactly four ID's of nodes are nesessary in 3D");
    return(GM_ERROR);
  }

  /* init vertices */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Vertex[i] = MYVERTEX(Node[i]);

  /* init pointers */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
  {
    Neighbor[i] = NULL;
    theSeg[i] = NULL;
  }

  /* compute side information (theSeg[i]==NULL) means inner side */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
  {
    for( j=0; j<MAX_CORNERS_OF_SIDE; j++ )
    {
      sideNode[j] = Node[CornerOfSide[i][j]];
      sideVertex[j] = Vertex[CornerOfSide[i][j]];
    }

    found = 0;
    for( j=0; j<MAX_CORNERS_OF_SIDE; j++ )
    {
      if( OBJT(sideVertex[j]) == IVOBJ ) found = 1;
    }
    if( found ) continue;

    /* all vertices of side[i] are on the boundary now */

    /* We now assume, that side[i] is on the boundary if and only if */
    /* there is one boundary segment containing the three nodes.        */
    /* That means, one should not define elements at the boundary	*/
    /* with a boundary side covering more than one segment.			*/

    for (j=0; j<theMG->numOfSegments; j++)
    {
      Seg = &(theMG->segments[j]);
      for( k=0; k<MAX_CORNERS_OF_SIDE; k++ )
      {
        found = 0;
        for( vs = VSEG(sideVertex[k]); vs!=NULL; vs = NEXTSEG(vs) )                           /* the segments of that corner */
        {
          if( BSEGDESC(vs) == Seg )                               /* corner belongs to segment */
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

    if (found)             /* a segment with three corners was found i.e. the side is a boundary side */
    {
      /*	set boundary parameters for vertices of that side */
      theSeg[i] = Seg;
      for( k=0; k<MAX_CORNERS_OF_SIDE; k++)                    /* vertices of side i */
      {
        param[i][k][0] = plambda[k][0];
        param[i][k][1] = plambda[k][1];
      }
    }
  }

  /* find neighboring elements */

  /* for all sides of the element to be created */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
  {
    for (j=0; j<MAX_CORNERS_OF_SIDE; j++)
    {
      sideNode[j] = Node[CornerOfSide[i][j]];
    }
    /* for all neighbouring elements allready inserted */
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      /* for all sides of the neighbour element */
      for (j=0; j<MAX_SIDES_OF_ELEM; j++)
      {
        num = 0;
        /* for all corners of the side of the neighbour */
        for (m=0; m<MAX_CORNERS_OF_SIDE; m++)
          /* for all corners of the side of the element to be created */
          for (k=0; k<MAX_CORNERS_OF_SIDE; k++)
          {
            if(CORNER(theElement,CornerOfSide[j][m])==sideNode[k]) num++;
          }
        if(num==MAX_CORNERS_OF_SIDE)
        {
          if (NBELEM(theElement,j)!=NULL)
          {
            PrintErrorMessage('E',"InsertElement","neighbor relation inconsistent");
            return(GM_ERROR);
          }
          Neighbor[i] = theElement;
          NeighborSide[i] = j;
        }
      }
    }
  }

  /* check type of element */

  found = 0;
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    if (theSeg[i]!=NULL)
      found++;

  /* create element */
  if ( found )
    theElement = CreateBoundaryElement(theGrid,NULL,4);
  else
    theElement = CreateInnerElement(theGrid,NULL,4);
  if (theElement==NULL)
  {
    PrintErrorMessage('E',"InsertElement","cannot allocate element");
    return(GM_ERROR);
  }

  /* create element sides if necessary */
  if (OBJT(theElement)==BEOBJ)
  {
    for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    {
      SET_SIDE(theElement,i,NULL);
      if (theSeg[i]!=NULL)
      {
        SET_SIDE(theElement,i,CreateElementSide(theGrid));
        if (SIDE(theElement,i)==NULL)
        {
          DisposeElement(theGrid,theElement);
          PrintErrorMessage('E',"InsertElement","cannot allocate element side");
          return(GM_ERROR);
        }
        SEGDESC(SIDE(theElement,i)) = theSeg[i];
        for(k=0; k<MAX_CORNERS_OF_SIDE; k++)
          for(l=0; l<DIM_OF_BND; l++)
            PARAM(SIDE(theElement,i),k,l) = param[i][k][l];
      }
    }
  }

  /* create edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
  {
    theEdge = GetEdge(Node[CornerOfEdge[i][0]],Node[CornerOfEdge[i][1]]);
    if (theEdge==NULL)
    {
      theEdge = CreateEdge(theGrid,Node[CornerOfEdge[i][0]],Node[CornerOfEdge[i][1]]);
      if (theEdge==NULL)
      {
        DisposeElement(theGrid,theElement);                         /* including element sides */
        PrintErrorMessage('E',"InsertElement","cannot allocate edge");
        return(GM_ERROR);
      }
      SET_NO_OF_ELEM(theEdge,1);
    }
    else
    {
      if (NO_OF_ELEM(theEdge)<NO_OF_ELEM_MAX-1)
        INC_NO_OF_ELEM(theEdge);
      else
        return(GM_ERROR);
    }
  }

  /* fill element data */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
  {
    SET_NBELEM(theElement,i,Neighbor[i]);
    if (Neighbor[i]!=NULL)
    {
      SET_NBELEM(Neighbor[i],NeighborSide[i],theElement);
                        #ifdef __SIDEDATA__
      if (DisposeDoubledSideVector(theGrid,Neighbor[i],NeighborSide[i],theElement,i))
        return(GM_ERROR);
                        #endif
    }
  }
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++) SET_CORNER(theElement,i,Node[i]);
  SETNSONS(theElement,0);
  SET_SON(theElement,0,NULL);
  SET_EFATHER(theElement,NULL);
  SETECLASS(theElement,RED);

  /* create connection to other elements. ATTENTION: this function is O(n) */
  if (InsertedElementCreateConnection(theGrid,theElement))
  {
    PrintErrorMessage('E',"InsertElement","could not create algebra connections");
    DisposeElement (theGrid,theElement);
    return(GM_ERROR);
  }

  /* set parity of element */
  if (SetParityOfElement(theElement)) return (1);


  return(GM_OK);
}
#endif


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
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* check data */
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++)
      if (idList[i]==idList[j])
      {
        PrintErrorMessage('E',"InsertElementFromIDs","nodes must be pairwise different");
        return(GM_ERROR);
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
    return(GM_ERROR);
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

#ifdef __TWODIM__
INT DeleteElement (MULTIGRID *theMG, ELEMENT *theElement) /* 2D VERSION */
{
  GRID *theGrid;
  NODE *aNode,*bNode;
  ELEMENT *theNeighbor;
  int i,j,n,m;
  EDGE *theEdge;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"DeleteElement","only a multigrid with exactly one level can be edited");
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* delete edges if possible */
  n = TAG(theElement);
  for (i=0; i<n; i++)
  {
    if (NBELEM(theElement,i)==NULL)
    {
      aNode = CORNER(theElement,i);
      bNode = CORNER(theElement,(i+1)%n);
      theEdge = GetEdge(aNode,bNode);
      if (theEdge!=NULL) DisposeEdge(theGrid,theEdge);
    }
  }
  if (n==4)
  {
    theEdge = GetEdge(CORNER(theElement,0),CORNER(theElement,2));
    if (theEdge!=NULL) DisposeEdge(theGrid,theEdge);
    theEdge = GetEdge(CORNER(theElement,1),CORNER(theElement,3));
    if (theEdge!=NULL) DisposeEdge(theGrid,theEdge);
  }

  /* delete pointers in neighbors */
  for (i=0; i<n; i++)
  {
    theNeighbor = NBELEM(theElement,i);
    if (theNeighbor!=NULL)
    {
      aNode = CORNER(theElement,i);
      bNode = CORNER(theElement,(i+1)%n);
      m = SIDES_OF_ELEM(theNeighbor);
      for (j=0; j<m; j++)
        if ((CORNER(theNeighbor,j)==bNode)&&(CORNER(theNeighbor,(j+1)%m)==aNode)) break;
      SET_NBELEM(theNeighbor,j,NULL);
    }
  }

  /* delete element now */
  DisposeElement(theGrid,theElement);

  return(GM_OK);
}
#endif

#ifdef __THREEDIM__
INT DeleteElement (MULTIGRID *theMG, ELEMENT *theElement) /* 3D VERSION */
{
  GRID *theGrid;
  ELEMENT *theNeighbor;
  LINK *theLink;
  EDGE *theEdge;
  int i,j,found;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"DeleteElement","only a multigrid with exactly one level can be edited");
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* delete edges if possible */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    for (j=i+1; j<CORNERS_OF_ELEM(theElement); j++)
    {
      theEdge = GetEdge(CORNER(theElement,i),CORNER(theElement,j));
      if (theEdge==NULL) return(GM_ERROR);
      if (NO_OF_ELEM(theEdge)<=0) return(GM_ERROR);
      if (NO_OF_ELEM(theEdge)==1)
        DisposeEdge(theGrid,theEdge);
      else
        DEC_NO_OF_ELEM(theEdge);
    }

    /* delete all EXTRA edges */
    for (theLink=START(CORNER(theElement,i)); theLink!=NULL; theLink=NEXT(theLink))
    {
      theEdge=MYEDGE(theLink);
      if (EXTRA(theEdge)) DisposeEdge(theGrid,theEdge);
    }
  }

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
      if (found!=1) return(GM_ERROR);
    }
  }

  /* delete element now */
  DisposeElement(theGrid,theElement);

  return(GM_OK);
}
#endif

INT DeleteElementWithID (MULTIGRID *theMG, INT id)
{
  GRID *theGrid;
  ELEMENT *theElement;

  /* check level */
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"DeleteElementWithId","only a multigrid with exactly one level can be edited");
    return(GM_ERROR);
  }
  theGrid = GRID_ON_LEVEL(theMG,0);

  /* find element */
  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    if (ID(theElement)==id) break;
  if (theElement==NULL)
  {
    PrintErrorMessage('E',"DeleteElementWithId","element not found");
    return(GM_ERROR);
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
INT PointInElement (const COORD *x, const ELEMENT *theElement) /* 3D version */
{
  COORD *point[MAX_CORNERS_OF_ELEM];
  COORD delta[MAX_CORNERS_OF_ELEM-1][DIM];
  COORD n[DIM];
  COORD a, b, norm;
  int i,j,k,rv;

  /* check element */
  if (theElement==NULL) return(0);

  /* load geometrical data of the corners */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    point[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  rv = 1;

  /* check for each side represented by a corner */
  for(i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    /* use this corner as reference point */
    for(j=0; j<CORNERS_OF_ELEM(theElement)-1; j++)
      for(k=0; k<DIM; k++)
        delta[j][k] = point[(i+j)%CORNERS_OF_ELEM(theElement)][k] - point[i][k];

    /* calculate a vector n[.] normal to that side */
    V3_VECTOR_PRODUCT(delta[0],delta[1], n);
    V3_EUKLIDNORM(n,norm);
    if(norm < SMALL_C) return(0);

    /* check if last point of tetrahedron is on the same side as x */
    V3_SCALAR_PRODUCT(n,x,a);
    V3_SCALAR_PRODUCT(n,delta[2],b);
    if (a*b < -SMALL_C*b*b) return(5);
    if (a*b <  SMALL_C*b*b) rv++;
  }

  /* x seems to lie on four sides: impossible error */
  if (rv == 5) rv = 0;

  /* return result */
  return(rv);
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

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    if (PointInElement(pos,theElement)==1) return(theElement);
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

  c = isCurrent ? '*' : ' ';

  if (longformat)
    sprintf(buffer," %c %-20.20s %-20.20s %-20.20s %10lu %10lu\n",c,ENVITEM_NAME(theMG),
            ENVITEM_NAME(theMG->theDomain),ENVITEM_NAME(theMG->theProblem),
            HeapSize(theMG->theHeap),HeapUsed(theMG->theHeap));
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
  int l;
  const GRID *theGrid;
  char c;

  sprintf(buffer,"grids of '%s':\n",ENVITEM_NAME(theMG));

  UserWrite("level maxlevel    #vert    #node    #edge    #elem    #side    #vect    #conn\n");
  for (l=0; l<=TOPLEVEL(theMG); l++)
  {
    theGrid = GRID_ON_LEVEL(theMG,l);

    c = (l==CURRENTLEVEL(theMG)) ? '*' : ' ';

    sprintf(buffer,"%c %3d %8d %8ld %8ld %8ld %8ld %8ld %8ld %8ld\n",c,l,(int)TOPLEVEL(theMG),
            (long)NV(theGrid),(long)NN(theGrid),(long)NE(theGrid),(long)NT(theGrid),
            (long)NS(theGrid),(long)NVEC(theGrid),(long)NC(theGrid));

    UserWrite(buffer);
  }
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
  FORMAT *theFormat;
  VERTEX *theVertex;
  VSEGMENT *vs;
  LINK *theLink;
  int i;

  theFormat = MGFORMAT(theMG);
  theVertex = MYVERTEX(theNode);

  /******************************/
  /* print standard information */
  /******************************/
  /* line 1 */ sprintf(buffer,"NODEID=%9ld CTRL=%8lx IX=%8ld VEID=%9ld LEVEL=%2d",(long)ID(theNode),(long)CTRL(theNode),
                       (long)INDEX(theNode),(long)ID(theVertex),LEVEL(theNode));
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
    /* line 2 */ sprintf(buffer,"   CLASS=%d NCLASS=%d",CLASS(theNode),NCLASS(theNode));
    UserWrite(buffer);

                #ifdef __version23__
    sprintf(buffer," NSKIP=%x VSKIP=%x ",(unsigned int)NSKIP(theNode),(unsigned int)VSKIP(theNode));
    UserWrite(buffer);
                #endif

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
        sprintf(buffer,"   SEGID=%9ld ",(long)SEGID(BSEGDESC(vs)) );
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

  /* handle dataoption */
        #ifdef __version23__
  if (dataopt)
  {
    if (theFormat->PrintNode!=NULL)
    {
      UserWrite("  ndata : ");
      (*theFormat->PrintNode)(NDATA(theNode),buffer);
      UserWrite(buffer);
      UserWrite("\n");
    }
    if (theFormat->PrintDiag!=NULL)
    {
      UserWrite("  ndiag : ");
      (*theFormat->PrintDiag)(NDIAG(theNode),buffer);
      UserWrite(buffer);
      UserWrite("\n");
    }
    if (theFormat->PrintVertex!=NULL)
    {
      UserWrite("  vdata : ");
      (*theFormat->PrintVertex)(VDATA(theVertex),buffer);
      UserWrite(buffer);
      UserWrite("\n");
    }
  }
        #endif

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
                        #ifdef __version23__
      if (dataopt)
      {
        if (theFormat->PrintLink!=NULL)
        {
          (*theFormat->PrintLink)(LDATA(theLink),buffer);
          UserWrite(buffer);
          UserWrite(" ");
        }
        if (theFormat->PrintEdge!=NULL)
        {
          (*theFormat->PrintEdge)(EDDATA(MYEDGE(theLink)),buffer);
          UserWrite(buffer);
        }
      }
      UserWrite("\n");
                        #endif
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
  FORMAT *theFormat;
  char etype[10];
  int i,j,k;
  ELEMENTSIDE *theSide;
        #ifdef __THREEDIM__
  ELEMENT *SonList[MAX_SONS];
        #endif

  theFormat = MGFORMAT(theMG);

  if (ECLASS(theElement)==COPY_CLASS) strcpy(etype,"COPY ");
  if (ECLASS(theElement)==IRREGULAR_CLASS) strcpy(etype,"IRREG");
  if (ECLASS(theElement)==REGULAR_CLASS) strcpy(etype,"REGUL");
  sprintf(buffer,"ELEMID=%9ld %5s CTRL=%8lx CTRL2=%8lx REFINE=%2d MARK=%2d LEVEL=%2d",(long)ID(theElement),etype,
          (long)CTRL(theElement),(long)FLAG(theElement),REFINE(theElement),MARK(theElement),LEVEL(theElement));
  UserWrite(buffer);
  if (COARSEN(theElement)) UserWrite(" COARSEN");
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
          sprintf(buffer,"SEGID%d=%ld ",i,(long)SEGID(SEGDESC(SIDE(theElement,i))));
          UserWrite(buffer);
          for(j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
          {
                                                #ifdef __THREEDIM__
            sprintf(buffer,"  NODE[ID=%ld]: ",(long)(ID(CORNER(theElement,CornerOfSide[i][j]))));
            UserWrite(buffer);
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

        #ifdef __version23__
  if (dataopt)
  {
    if (theFormat->PrintElement!=NULL)
    {
      UserWrite("  edata : ");
      (*theFormat->PrintElement)(EDATA(theElement),buffer);
      UserWrite(buffer);
      UserWrite("\n");
    }
  }
        #endif

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

#ifdef __version3__

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
#endif


/****************************************************************************/
/*D
   ListVectorSelection - List info about vectors of elements in selection

   SYNOPSIS:
   void ListVectorSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt);

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

#ifdef __version3__

void ListVectorSelection (MULTIGRID *theMG, INT matrixopt, INT dataopt)
{
  int i,j;
  ELEMENT *theElement;
  VECTOR *vList[20];
  INT cnt;

  if (SELECTIONMODE(theMG) != elementSelection)
  {
    PrintErrorMessage('E',"ListElementSelection","wrong selection type");
    return;
  }
  for (j=0; j<SELECTIONSIZE(theMG); j++)
  {
    theElement = (ELEMENT *) SELECTIONOBJECT(theMG,j);
    sprintf(buffer,"ELEM(ID=%d):\n",ID(theElement));
    UserWrite(buffer);
    GetVectorsOfNodes(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    GetVectorsOfEdges(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    GetVectorsOfSides(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
    GetVectorsOfElement(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) ListVector(theMG,vList[i],matrixopt,dataopt);
  }
}
#endif


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

#ifdef __version3__

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
#endif


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
  BNDSEGDESC *theSegment;

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
                theSegment = SEGDESC(theSide);
                for (vs=VSEG(theVertex); vs!=NULL; vs=NEXTSEG(vs))
                  if (theSegment == BSEGDESC(vs))
                    break;
                if (vs == NULL)
                  *NodeError |= (i<<k);
              }
              else
                *NodeError |= (i<<(k+CORNERS_OF_ELEM(theElement)));
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
                        #ifdef __version3__
      /* there are no diagonal edges in quadrilaterals */
      if ((n==4) && (j==i+2))
        continue;
                        #endif
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

  /* reset used flags */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    SETUSED(theNode,0);
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
      SETUSED(MYEDGE(theLink),0);
  }

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
      if (ID(NBNODE(theLink))>ID(theNode))
      {
        theEdge = MYEDGE(theLink);
        if (!USED(theEdge))
        {
          sprintf(buffer,"edge between %ld and %ld dead ",(long)ID(theNode),(long)ID(NBNODE(theLink)));
          UserWrite(buffer);
          UserWrite("\n");
        }
        else
          SETUSED(theEdge,0);
      }
    }
  }

  /* look for dead nodes */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    if (!USED(theNode))
    {
      sprintf(buffer,"node %ld is dead ",(long)ID(theNode));
      UserWrite(buffer);
      UserWrite("\n");
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
  BNDSEGDESC *theSegment;

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
              theVertex = MYVERTEX(CORNER(theElement,(k=CornerOfSide[i][j])));
              if (OBJT(theVertex) == BVOBJ)
              {
                theSegment = SEGDESC(theSide);
                for (vs=VSEG(theVertex); vs!=NULL; vs=NEXTSEG(vs))
                  if (theSegment == BSEGDESC(vs))
                    break;
                if (vs == NULL)
                  *NodeError |= (i<<k);
              }
              else
                *NodeError |= (i<<(k+MAX_CORNERS_OF_ELEM));
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
              theVertex = MYVERTEX(CORNER(theElement,(k=CornerOfSide[i][j])));
              if (OBJT(theVertex) == BVOBJ)
              {
                theSegment = SEGDESC(theSide);
                for (vs=VSEG(theVertex); vs!=NULL; vs=NEXTSEG(vs))
                  if (theSegment == BSEGDESC(vs))
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
    for (j=i+1; j<CORNERS_OF_ELEM(theElement); j++)
    {
      theEdge = GetEdge(CORNER(theElement,i),CORNER(theElement,j));
      if (theEdge==NULL)
        *EdgeError |= 1<<EdgeWithCorners[i][j];
      else
        SETUSED(theEdge,1);
    }
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
      SETUSED(MYEDGE(theLink),0);
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
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,CornerOfSide[i][j])));
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
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,CornerOfSide[i][j])));
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
            sprintf(buffer,"%ld",(long)ID(CORNER(theElement,CornerOfSide[i][j])));
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
        sprintf(buffer,"   EDGE(%ld,%ld) is missing\n",(long)ID(CORNER(theElement,CornerOfEdge[i][0])),(long)ID(CORNER(theElement,CornerOfEdge[i][1])));
        UserWrite(buffer);
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
      if (ID(NBNODE(theLink))>ID(theNode))
      {
        theEdge = MYEDGE(theLink);
        if (!USED(theEdge))
        {
          sprintf(buffer,"edge between %ld and %ld dead, NO_OF_ELEM=%d ",(long)ID(theNode),(long)ID(NBNODE(theLink)),NO_OF_ELEM(theEdge));
          UserWrite(buffer);
          UserWrite("\n");
        }
        else
          SETUSED(theEdge,0);
      }
    }
  }

  /* look for dead nodes */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    if (!USED(theNode))
    {
      sprintf(buffer,"node %ld is dead ",(long)ID(theNode));
      UserWrite(buffer);
      UserWrite("\n");
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
  GEOM_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=nodeSelection) return(GM_ERROR);
  }
  else SELECTIONMODE(theMG) = nodeSelection;

  g = (GEOM_OBJECT *) theNode;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g) return(GM_ERROR);

  if (SELECTIONSIZE(theMG)>=MAXSELECTION) return(GM_ERROR);

  SELECTIONOBJECT(theMG,SELECTIONSIZE(theMG)) = g;
  SELECTIONSIZE(theMG)++;
  return(GM_OK);
}


/****************************************************************************/
/*D
   AddElementToSelection - Add element to selection buffer

   SYNOPSIS:
   INT AddElementToSelection (MULTIGRID *theMG, NODE *theElement);

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
  GEOM_OBJECT *g;

  if (SELECTIONSIZE(theMG)!=0)
  {
    if (SELECTIONMODE(theMG)!=elementSelection) return(GM_ERROR);
  }
  else SELECTIONMODE(theMG) = elementSelection;

  g = (GEOM_OBJECT *) theElement;
  for (i=0; i<SELECTIONSIZE(theMG); i++)
    if (SELECTIONOBJECT(theMG,i)==g) return(GM_ERROR);

  if (SELECTIONSIZE(theMG)>=MAXSELECTION) return(GM_ERROR);

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
  GEOM_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=nodeSelection) return(GM_ERROR);
  }
  else return(GM_ERROR);

  g = (GEOM_OBJECT *) theNode;
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

   SYNOPSIS:
   INT RemoveElementFromSelection (MULTIGRID *theMG, NODE *theElement);

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
  GEOM_OBJECT *g;

  if (SELECTIONSIZE(theMG)>0)
  {
    if (SELECTIONMODE(theMG)!=elementSelection) return(GM_ERROR);
  }
  else return(GM_ERROR);

  g = (GEOM_OBJECT *) theElement;
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
    if (l<SMALL_D) return(1);
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
  if (QualityElement(n,x,y,amin,amax)) return(GM_ERROR);
  return(GM_OK);
}
#endif

#ifdef __THREEDIM__
static INT QualityElement (int type, ELEMENT *element, DOUBLE *angle)
{
  int i,j,k,errorcode;
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


  if (TetraSideNormals (x,theNormal)) return (2);
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
  {
    V3_SCALAR_PRODUCT(theNormal[SideWithEdge[i][0]],theNormal[SideWithEdge[i][1]],help);
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
