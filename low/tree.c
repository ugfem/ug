// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                    */
/* File:      tree.c                                                        */
/*                                                                          */
/* Purpose:   quadtree and octtree                                          */
/*                                                                          */
/* Author:      Carsten Schwarz                                             */
/*              Institut fuer Hydromechanik und Wasserwirtschaft            */
/*              ETH Hoenggerberg                                            */
/*              8093 Zuerich                                                */
/*                                                                          */
/* History:   07.04.97 begin, ug version 3.4                                */
/*                                                                          */
/* Revision:  07.04.97                                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*              system include files                                        */
/*              application include files                                   */
/*                                                                          */
/****************************************************************************/

#include <limits.h>
#include <stdlib.h>
#ifndef __MWCW__
#include <sys/types.h>
#endif

#include "compiler.h"
#include "general.h"
#include "fifo.h"
#include "heaps.h"

#include "tree.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define EPSILON 1e-10

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*D
   CreateTree - Create root of search tree

   SYNOPSIS:
   TREE *CreateTree(HEAP *theHeap, INT dim, DOUBLE *posrange)

   PARAMETERS:
   .  HEAP - reference to Heap to get memory from
   .  dim - dimension of underlaying space
   .  posrange - array of ranges of positions in order x1min,x1max,x2min,...

   DESCRIPTION:
   This function creates a new search tree.

   For dim equal to one the tree is a binary tree, for dim equal to two
   a quarttree...

   RETURN VALUE:
   TREE *
   .n     pointer to the tree
   .n     NULL if an error occurred
   D*/
/****************************************************************************/
TREE *CreateTree (HEAP *theHeap, INT dim, DOUBLE *posrange)
{
  TREE *newTree;
  INT i;

  if ((newTree = (TREE *) GetFreelistMemory(theHeap,sizeof(TREE)+sizeof(DOUBLE)*(4*dim-1)))
      == NULL)
    return(NULL);

  newTree->status = TREE_CHANGED;
  newTree->fifo_max_mem = (size_t) 0;
  newTree->heap = theHeap;
  newTree->fifo = NULL;
  newTree->dim = dim;
  TREEROOT(newTree) = NULL;
  for (i=0; i<dim; i++)
  {
    TREEPOS(newTree,i,0) = posrange[i];
    TREEPOS(newTree,i,1) = posrange[i+dim];
  }

  return(newTree);
}

/****************************************************************************/
/*D
   DeleteTree - Removes a search tree

   SYNOPSIS:
   INT DeleteTree(TREE *theTree);

   PARAMETERS:
   .  theTree - reference to the tree

   DESCRIPTION:
   This function deletes a search tree and frees all the memory
   occupied by the tree and its objects.

   RETURN VALUE:
   INT
   .n     0 if okay
   .n     1 if an error occurred
   D*/
/****************************************************************************/
INT DeleteTree(TREE *theTree)
{
  FIFO fifo;
  void *buffer;
  TREE_ENTRY *entry, *nextEntry;

  if (theTree == NULL)
  {
    return(1);
  }

  if (TREEROOT(theTree) != NULL)
  {
    if (theTree->status != TREE_SEARCH)
    {
      if ((buffer = (void *) GetFreelistMemory(theTree->heap,theTree->fifo_max_mem)) == NULL)
      {
        /* Cannot allocate fifo to dispose all leafs of tree */
        /* So I only dispose major entries */
        PutFreelistMemory(theTree->heap,TREEROOT(theTree),sizeof(TREE_ENTRY)+ (theTree->dim*2-1)*sizeof(DOUBLE));
        PutFreelistMemory(theTree->heap,theTree,sizeof(TREE)+sizeof(DOUBLE)*(4*theTree->dim-1));
        return(1);
      }
      fifo_init(&fifo,buffer,theTree->fifo_max_mem);
    }
    else
    {
      fifo = *theTree->fifo;
      fifo_clear(&fifo);
    }

    fifo_in(&fifo,(void *)theTree->root);
    while(!fifo_empty(&fifo))
    {
      entry = (TREE_ENTRY *)fifo_out(&fifo);

      switch (TNODETYPE(entry))
      {
      case TREELEAF :
        PutFreelistMemory(theTree->heap, entry,sizeof(TREE_LEAF) + (theTree->dim-1)*sizeof(DOUBLE));
        break;

      case TREENODE :
        if ((nextEntry = TNODESON(entry)) != NULL)
          fifo_in(&fifo,(void *)nextEntry);
        if ((nextEntry = TNODENEXT(entry)) != NULL)
          fifo_in(&fifo,(void *)nextEntry);
        PutFreelistMemory(theTree->heap,entry,sizeof(TREE_ENTRY)+(2*theTree->dim-1)*sizeof(DOUBLE));
        break;
      }
    }
  }
  if (theTree->fifo != NULL)
    PutFreelistMemory(theTree->heap, theTree->fifo, theTree->fifo_max_mem);
  PutFreelistMemory(theTree->heap, theTree, sizeof(TREE)+(4*theTree->dim-1)*sizeof(DOUBLE));

  return(0);
}

static INT PointInTNODE(TREE_ENTRY *entry, DOUBLE *pos, INT dim)
{
  INT i;

  for (i=0; i<dim; i++)
  {
    if (TNODEPOS(entry,i,0,dim) >= pos[i]) return(0);
    if (TNODEPOS(entry,i,1,dim) <  pos[i]) return(0);
  }
  return(1);
}

/****************************************************************************/
/*D
   SearchInsertPoint - Gets position of insert point

   SYNOPSIS:
   INT SearchInsertPoint(TREE *theTree, DOUBLE *position,
                                  TREE_ENTRY **insertPoint)

   PARAMETERS:
   .  theTree - reference to the tree
   .  Position - position of the object
   . insertPoint - Pointer to result

   DESCRIPTION:
   This function searches the tree for the point for insertation of a new
   object. *insertPoint points to the insertation place and the return value
   indicates if the insertation point is at a free end of the tree or already
   occupied.

   RETURN VALUE:
   INT
   .n 0 - tree not correct
   .n 1 - insertation point is a free end
   .n 2 - place is occupied
   D*/

/****************************************************************************/
static INT SearchInsertPoint(TREE *theTree, DOUBLE *position, TREE_ENTRY **insertPoint)
{
  TREE_ENTRY *act, *father;

  *insertPoint = NULL;
  if (theTree == NULL) return(0);
  if (theTree->root == NULL) return(1);

  act = theTree->root;
  father = NULL;
  while (act != NULL)
  {
    if (TNODETYPE(act) == TREELEAF)
    {
      *insertPoint = father;
      return(2);
    }
    while (!PointInTNODE(act, position, theTree->dim))
    {
      if (TNODENEXT(act) == NULL)
      {
        *insertPoint = act;
        return(1);
      }
      act = TNODENEXT(act);
    }
    father = act;
    act = TNODESON(act);
  }
  return(0);
}


/****************************************************************************/
/*D
   InsertinTree - Inserts object in search tree

   SYNOPSIS:
   void *InsertinTree(TREE *theTree, DOUBLE *Position);

   PARAMETERS:
   .  theTree - reference to the tree
   .  Position - position of the object
   .  obj - pointer to object

   DESCRIPTION:
   This function changes the search tree for a new object at Position.
   The tree leaf will hold a pointer to the object which has to exist
   before.

   RETURN VALUE:
   INT
   .n     error code
   .n     0 if okay
   D*/
/****************************************************************************/
INT InsertinTree (TREE *theTree, DOUBLE *Position, void *obj)
{
  TREE_ENTRY *entry, *pentry, *newentry;
  INT i,dim, collide;
  DOUBLE middle;

  if (theTree->status == TREE_INVALID)
    return(0);
  dim = theTree->dim;
  for (i=0; i < dim; i++)
  {
    if (theTree->posrange[i] >= Position[i])
      return(NULL);
    if (theTree->posrange[i+dim] < Position[i])
      return(NULL);
  }

  switch(SearchInsertPoint(theTree,Position,&entry))
  {
  case 1 :
    if ((newentry = (TREE_ENTRY *) GetFreelistMemory(theTree->heap, sizeof(TREE_ENTRY)
                                                     + (dim*2-1)*sizeof(DOUBLE)))
        == NULL)
      return (1);
    TNODENEXT(newentry) = NULL;
    TNODESON(newentry) = NULL;
    TNODETYPE(newentry) = TREENODE;
    theTree->status = TREE_CHANGED;
    theTree->fifo_max_mem += sizeof (void *);
    if (entry == NULL)
    {
      TNODEFATHER(newentry) = NULL;
      TREEROOT(theTree) = newentry;
      for (i=0; i<dim; i++)
      {
        TNODEPOS(newentry,i,0,dim) = theTree->posrange[i];
        TNODEPOS(newentry,i,1,dim) = theTree->posrange[i+dim];
      }
    }
    else
    {
      TNODEFATHER(newentry) = TNODEFATHER(entry);
      TNODENEXT(entry) = newentry;
      for (i=0; i<dim; i++)
        if ((middle = (TNODEPOS(TNODEFATHER(newentry),i,0,dim)+TNODEPOS(TNODEFATHER(newentry),i,1,dim))/2.0)
            >= Position[i])
        {
          TNODEPOS(newentry,i,0,dim) = TNODEPOS(TNODEFATHER(newentry),i,0,dim);
          TNODEPOS(newentry,i,1,dim) = middle;
        }
        else
        {
          TNODEPOS(newentry,i,0,dim) = middle;
          TNODEPOS(newentry,i,1,dim) = TNODEPOS(TNODEFATHER(newentry),i,1,dim);
        }
    }

    break;

  case 2 :
    /* refine until points are apart */
    pentry = entry;
    entry = TNODESON(pentry);
    do {
      if ((newentry = (TREE_ENTRY *) GetFreelistMemory(theTree->heap, sizeof(TREE_ENTRY)
                                                       + (dim*2-1)*sizeof(DOUBLE)))
          == NULL)
        return (1);
      TNODESON(pentry) = newentry;
      TNODENEXT(newentry) = NULL;
      TNODESON(newentry) = NULL;
      TNODETYPE(newentry) = TREENODE;
      TNODEFATHER(newentry) = pentry;
      theTree->status = TREE_CHANGED;
      theTree->fifo_max_mem += sizeof(void *);
      collide = 1;
      for (i=0; i<dim; i++)
        if ((middle = (TNODEPOS(pentry,i,0,dim)+TNODEPOS(pentry,i,1,dim))/2.0)
            >= TLEAFPOS(entry,i))
        {
          TNODEPOS(newentry,i,0,dim) = TNODEPOS(pentry,i,0,dim);
          TNODEPOS(newentry,i,1,dim) = middle;
          if (Position[i] > middle)
            collide = 0;
        }
        else
        {
          TNODEPOS(newentry,i,0,dim) = middle;
          TNODEPOS(newentry,i,1,dim) = TNODEPOS(pentry,i,1,dim);
          if (Position[i] <= middle)
            collide = 0;
        }
      pentry = newentry;
    } while (collide);
    if ((newentry = (TREE_ENTRY *) GetFreelistMemory(theTree->heap, sizeof(TREE_ENTRY)
                                                     + (dim*2-1)*sizeof(DOUBLE)))
        == NULL)
      return (1);
    TNODENEXT(pentry) = newentry;
    TNODESON(newentry) = NULL;
    TNODENEXT(newentry) = NULL;
    TNODEFATHER(newentry) = TNODEFATHER(pentry);
    TNODETYPE(newentry) = TREENODE;
    theTree->fifo_max_mem += sizeof(void *);
    for (i=0; i<dim; i++)
      if ((middle = (TNODEPOS(TNODEFATHER(newentry),i,0,dim)+TNODEPOS(TNODEFATHER(newentry),i,1,dim))/2.0)
          >= Position[i])
      {
        TNODEPOS(newentry,i,0,dim) = TNODEPOS(TNODEFATHER(newentry),i,0,dim);
        TNODEPOS(newentry,i,1,dim) = middle;
      }
      else
      {
        TNODEPOS(newentry,i,0,dim) = middle;
        TNODEPOS(newentry,i,1,dim) = TNODEPOS(TNODEFATHER(newentry),i,1,dim);
      }

    TNODESON(pentry) = entry;
    TNODEFATHER(entry) = pentry;
    break;

  default :
    /* not a legal tree, should not happen! */
    return(1);
  }

  /* append the object to the leaf node */
  if ((TNODESON(newentry) = (TREE_ENTRY *) GetFreelistMemory(theTree->heap,sizeof(TREE_LEAF)
                                                             + (theTree->dim-1)*sizeof(DOUBLE)))
      == NULL)
  {
    theTree->status = TREE_INVALID;
    return(1);
  }
  TNODEFATHER(TNODESON(newentry)) = newentry;
  entry = TNODESON(newentry);

  TNODETYPE(entry) = TREELEAF;
  for (i=0; i<theTree->dim; i++)
    entry->tleaf.pos[i] = Position[i];
  entry->tleaf.obj = obj;
  return (0);
}

/****************************************************************************/
/*D
   DeleteObjinTree - Deletes object in search tree

   SYNOPSIS:
   void *DeleteObjinTree(TREE *theTree, DOUBLE *Position);

   PARAMETERS:
   .  theTree - reference to the tree
   .  Position - position of the object

   DESCRIPTION:
   This function deletes the tree leaf at Position and all
   unnecassary tree nodes from the search tree.

   RETURN VALUE:
   void *
   .n     pointer to object at delete position
   .n     NULL if none or error occured
   D*/
/****************************************************************************/
void *DeleteObjinTree (TREE *theTree, DOUBLE *Position)
{
  TREE_ENTRY *entry, *pentry;
  INT i,dim;
  HEAP *heap;
  void *obj;

  if (theTree->status == TREE_INVALID)
    return(NULL);
  dim = theTree->dim;
  for (i=0; i < dim; i++)
  {
    if (theTree->posrange[i] >= Position[i])
      return(NULL);
    if (theTree->posrange[i+dim] <= Position[i])
      return(NULL);
  }

  switch(SearchInsertPoint(theTree,Position,&entry))
  {
  case 1 :
    /* object does not exist */
    return(NULL);
  case 2 :
    /* object exists? */
    for (i=0; i<dim; i++)
    {
      if (Position[i] < TLEAFPOS(TNODESON(entry),i)-EPSILON)
        return(NULL);                         /* Object does not exist */
      if (Position[i] > TLEAFPOS(TNODESON(entry),i)+EPSILON)
        return(NULL);                         /* Object does not exist */
    }

    /* delete it */
    heap =theTree->heap;
    obj = TNODEOBJ(TNODESON(entry));
    PutFreelistMemory(heap,TNODESON(entry),sizeof(TREE_ENTRY)+(dim-1)*sizeof(DOUBLE));
    if (entry == TREEROOT(theTree))
    {
      PutFreelistMemory(heap,entry,sizeof(TREE_ENTRY)+(2*dim-1)*sizeof(DOUBLE));
      theTree->fifo_max_mem = 0;
      TREEROOT(theTree) = NULL;
      theTree->status = TREE_CHANGED;
      return(obj);
    }
    pentry = TNODEFATHER(entry);
    if (entry==TNODESON(pentry))
    {
      TNODESON(pentry)=TNODENEXT(entry);
      PutFreelistMemory(heap,entry,sizeof(TREE_ENTRY)+(2*dim-1)*sizeof(DOUBLE));
      theTree->fifo_max_mem -= sizeof(void *);
      entry = TNODESON(pentry);
    }
    else
    {
      for (pentry=TNODESON(pentry); TNODENEXT(pentry)!=entry; pentry=TNODENEXT(pentry))
      {}
      TNODENEXT(pentry)=TNODENEXT(entry);
      PutFreelistMemory(heap,entry,sizeof(TREE_ENTRY)+(2*dim-1)*sizeof(DOUBLE));
      theTree->fifo_max_mem -= sizeof(void *);
      entry=TNODESON(TNODEFATHER(pentry));
    }
    theTree->status = TREE_CHANGED;

    /* try to move up the tree */
    while (TNODENEXT(entry) == NULL)
    {
      if ((pentry = TNODEFATHER(entry)) == NULL)
        return(0);
      else if (TNODESON(pentry) == entry)
      {
        TNODESON(pentry) = TNODESON(entry);
        PutFreelistMemory(heap,entry,sizeof(TREE_ENTRY)+(2*dim-1)*sizeof(DOUBLE));
        theTree->fifo_max_mem -= sizeof(void *);
        entry = pentry;
      }
      else
        return(obj);
    }
    return(obj);

  default :
    return(NULL);
  }
}

/****************************************************************************/
/*D
   GetFirstLeafinQuader,GetNextLeafinQuader - Gives objects whitch are inside a hyper-rectangular

   SYNOPSIS:
   void *GetFirstLeafinQuader(TREE *theTree, DOUBLE *ll, DOUBLE *ur);
   .n void *GetNextLeafinQuader(TREE *theTree);

   PARAMETERS:
   .  theTree - reference to the tree
   .  ll: lower left ... corner coordinates
   .  ur: upper right ... corner coordinates

   DESCRIPTION:
   GetFirstLeafinQuader initializes the fifo for a search and searches for
   a first leaf entry in the tree theTree which is lying inside the quader
   with lower left ... corner ll and upper right ... corner ur.
   .n GetNextLeafinQuader searches for the next object in the quader.

   RETURN VALUE:
   INT
   .n     Pointer to the TREE_ENTRY
   .n     NULL if no leaf entry inside exists or an error occured
   D*/
/****************************************************************************/
TREE_ENTRY *GetFirstLeafinQuader(TREE *theTree, DOUBLE *ll, DOUBLE *ur)
{
  INT i, dim, inside;
  void *buffer;
  TREE_ENTRY *entry;

  if (TREEROOT(theTree) == NULL)
    return(NULL);             /*Wo nix ist, brauch ich auch nicht zu suchen */

  /* Get the fifo and its memory */
  switch (theTree->status)
  {
  case TREE_CHANGED :
    if ((buffer = (void *) GetFreelistMemory(theTree->heap,theTree->fifo_max_mem)) == NULL)
      /* Cannot allocate fifo */
      return(NULL);
    if ((theTree->fifo = (FIFO *) GetFreelistMemory(theTree->heap,sizeof(FIFO))) == NULL)
      return(NULL);
    fifo_init(theTree->fifo,buffer,theTree->fifo_max_mem);
    theTree->status = TREE_SEARCH;
    break;

  case TREE_SEARCH :
    fifo_clear(theTree->fifo);
    break;

  default :
    return(NULL);
  }

  /* Save the question */
  dim = theTree->dim;
  for (i=0; i<dim; i++)
  {
    TREESEARCHLL(theTree,i) = ll[i];
    TREESEARCHUR(theTree,i) = ur[i];
  }

  /* Search! */
  fifo_in(theTree->fifo,TREEROOT(theTree));
  while(!fifo_empty(theTree->fifo))
  {
    entry = (TREE_ENTRY*)fifo_out(theTree->fifo);
    inside = 1;
    if (TNODETYPE(entry) == TREELEAF)
    {
      for(i=0; i<dim; i++)
      {
        if (TLEAFPOS(entry,i)<=ll[i])
        {
          inside = 0;
          break;
        }
        if (TLEAFPOS(entry,i)>ur[i])
        {
          inside = 0;
          break;
        }
      }
      if (inside == 1)
        return(entry);
    }
    else
    {
      for (i=0; i<dim; i++)
      {
        if(TNODEPOS(entry,i,1,dim)<=ll[i])
        {
          inside = 0;
          break;
        }
        if(TNODEPOS(entry,i,0,dim)>ur[i])
        {
          inside = 0;
          break;
        }
      }
      if (inside == 1)
        fifo_in(theTree->fifo,TNODESON(entry));
      if (TNODENEXT(entry)!=NULL)
        fifo_in(theTree->fifo,TNODENEXT(entry));
    }
  }
  return(NULL);
}

TREE_ENTRY *GetNextLeafinQuader(TREE *theTree)
{
  INT i, dim, inside;
  TREE_ENTRY *entry;

  if (theTree->status != TREE_SEARCH)
    return(NULL);             /* Not initialized */
  dim = theTree->dim;

  /* Search! */
  while(!fifo_empty(theTree->fifo))
  {
    entry = (TREE_ENTRY*)fifo_out(theTree->fifo);
    inside = 1;
    if (TNODETYPE(entry) == TREELEAF)
    {
      for(i=0; i<dim; i++)
      {
        if (TLEAFPOS(entry,i)<=TREESEARCHLL(theTree,i))
        {
          inside = 0;
          break;
        }
        if (TLEAFPOS(entry,i)>TREESEARCHUR(theTree,i))
        {
          inside = 0;
          break;
        }
      }
      if (inside == 1)
        return(entry);
    }
    else
    {
      for (i=0; i<dim; i++)
      {
        if(TNODEPOS(entry,i,1,dim)<=TREESEARCHLL(theTree,i))
        {
          inside = 0;
          break;
        }
        if(TNODEPOS(entry,i,0,dim)>TREESEARCHUR(theTree,i))
        {
          inside = 0;
          break;
        }
      }
      if (inside == 1)
        fifo_in(theTree->fifo,TNODESON(entry));
      if (TNODENEXT(entry)!=NULL)
        fifo_in(theTree->fifo,TNODENEXT(entry));
    }
  }
  return(NULL);
}
