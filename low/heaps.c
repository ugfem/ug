// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      heaps.c                                                       */
/*                                                                          */
/* Purpose:   low-level memory management for ug 2.0                        */
/*                                                                          */
/* Author:      Peter Bastian/Henrik Rentz-Reichert                         */
/*              Institut fuer Computeranwendungen III                       */
/*              Universitaet Stuttgart                                      */
/*              Pfaffenwaldring 27                                          */
/*              70569 Stuttgart                                             */
/*              email: ug@ica3.uni-stuttgart.de                             */
/*                                                                          */
/* History:   29.01.92 begin, ug version 2.0                                */
/*              02.02.95 begin, ug version 3.0                              */
/*                                                                          */
/* Revision:  04.09.95                                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files in the order                                               */
/*              system include files                                        */
/*              application include files                                   */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "compiler.h"
#include "heaps.h"
#include "misc.h"
#include "general.h"
#include "debug.h"
#include "devices.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

#define FLOOR(n)    ((n)&ALIGNMASK)     /* lower next multiple of four */


/* defines and macros for the virtual heap management                        */

#define B_OFFSET(bhm,i)         ((bhm)->BlockDesc[i].offset)
#define B_SIZE(bhm,i)            ((bhm)->BlockDesc[i].size)
#define B_ID(bhm,i)             ((bhm)->BlockDesc[i].id)

#define CALC_B_OFFSET(bhm,i)    (((i)==0) ? 0 : (B_OFFSET(theVHM,(i)-1)+B_SIZE(theVHM,(i)-1)))

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

REP_ERR_FILE;

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*D
   HeapStat - Get information on heap

   SYNOPSIS:
   MEM HeapStat (const HEAP *theHeap)

   PARAMETERS:
   .  theHeap - heap to get information

   DESCRIPTION:
   This function gets information on heap (objects and memory in freelists).

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void HeapStat (const HEAP *theHeap)
{
  INT i;
  INT usedfreelistentries,size,found;
  void ** ptr;

  usedfreelistentries = 0;

  UserWriteF("HeapStat: heap=%08x type=%d\n",theHeap,theHeap->type);
        #ifdef Debug
  UserWriteF("FreelistInfo:\n");
        #endif
  for (i=0; i<MAXFREEOBJECTS; i++)
  {
    size = theHeap->SizeOfFreeObjects[i];
    if (size != -1)
    {
                        #ifdef Debug
      /* inspect linked list */
      ptr = (void **) theHeap->freeObjects[i];
      found = 0;
      while (ptr != NULL)
      {
        ptr = (void **)ptr[0];
        found++;
      }
      UserWriteF("Entry %4d: objsize=%d objcount=%d entrymem=%d found=%d\n",
                 i,size,theHeap->objcount[i],size*theHeap->objcount[i],found);
                        #endif
      usedfreelistentries++;
    }
                #ifdef Debug
    else
      assert(theHeap->freeObjects[i] == NULL);
                #endif
  }
  UserWriteF("          size           =%d\n",theHeap->size);
  UserWriteF("          used           =%d\n",theHeap->used);
  UserWriteF("          freelistmem    =%d\n",theHeap->freelistmem);
  UserWriteF("          MAXFREEOBJECTS =%d\n",MAXFREEOBJECTS);
  UserWriteF("          usedfreelistent=%d\n",usedfreelistentries);
}


/****************************************************************************/
/*D
   InitHeaps - Initialize memory management module

   SYNOPSIS:
   INT InitHeaps ();

   PARAMETERS:
   void

   DESCRIPTION:
   This function initializes memory management module.

   CAUTION: code may be machine dependent.

   RETURN VALUE:
   INT
   .n     0 if module initialized correctly
   D*/
/****************************************************************************/

INT InitHeaps ()
{
  return(0);
}


/****************************************************************************/
/*D
   NewHeap - Install a new heap structure

   SYNOPSIS:
   HEAP *NewHeap (INT type, MEM size, void *buffer);

   PARAMETERS:
   .  type - type of heap
   .  size - size of new heap in bytes
   .  buffer - 4-aligned memory for the heap

   DESCRIPTION:
   This function installs a new heap structure.
   Valid 'type' is either 'SIMPLE_HEAP' or 'GENERAL_HEAP'.
   The allocation of memory starts at the address given by
   '*buffer' and is of 'size' bytes.

   RETURN VALUE:
   HEAP *
   .n      pointer to HEAP
   .n      NULL if not enough space available
   D*/
/****************************************************************************/

HEAP *NewHeap (INT type, MEM size, void *buffer)
{
  HEAP *theHeap;
  INT i;

  /* check size */
  if (buffer==NULL) return(NULL);
  if (size<MIN_HEAP_SIZE) return(NULL);

  /* initialize heap structure */
  theHeap = (HEAP *) buffer;
  theHeap->type = type;
  theHeap->size = size;
  theHeap->freelistmem = 0;
  theHeap->topStackPtr = theHeap->bottomStackPtr = 0;
  theHeap->heapptr = (BLOCK *) CEIL(((MEM)theHeap)+sizeof(HEAP));
  theHeap->used = ((MEM)theHeap->heapptr)-((MEM)theHeap);

  /* initialize first block */
  theHeap->heapptr->size = ((MEM)theHeap)+size-((MEM)theHeap->heapptr);
  theHeap->heapptr->next = theHeap->heapptr;
  theHeap->heapptr->previous = theHeap->heapptr;

  /* initialize free lists */
  for (i=0; i<MAXFREEOBJECTS; i++)
  {
    theHeap->SizeOfFreeObjects[i] = -1;;
    theHeap->freeObjects[i] = NULL;
  }
  /* return heap structure */
  return(theHeap);
}

/****************************************************************************/
/*D
   GetMem - Allocate memory from heap, depending on heap type

   SYNOPSIS:
   void *GetMem (HEAP *theHeap, MEM n, INT mode);

   PARAMETERS:
   .  theHeap - heap to allocate from
   .  n - number of bytes to allocate
   .  mode - allocation position for mark/release heap

   DESCRIPTION:
   This function allocates memory from 'HEAP', depending on heap type.

   If the heap type (theHeap->type) is 'SIMPLE_HEAP' new blocks in the
   allocated memory can be taken either from the top (mode = 'FROM_TOP')
   or from the bottom (mode = 'FROM_BOTTOM') of this block. In other words
   the total memory block to be provided for ug is used from both sides by
   introducing recursively new, smaller blocks.

   .vb
      --------------------------------
 | 1 | 3 |  4  |           |  2 |   total allocated memory for ug
      --------------------------------
      used blocks are #1,#2,#3 and #4

      -----
 ||5|6|   block #4 is separated also in block #5 and #6
      -----
   .ve

   If the heap type (theHeap->type) is 'GENERAL_HEAP' new blocks to be
   introduced in the total allocated memory for ug are laid at the position
   where enough memory is free. The search for this position is done by
   running round the memory and looking for the equivalent space.

   .vb
      --------------------------------
 |   |  1  |               |  2 |   total allocated memory for ug
      --------------------------------
      used blocks are #1 and #2

      block to be introduced #3:  |  3  |
      start position at the beginning

      --------------------------------
 |   |  1  |  3  |         |  2 |   total allocated memory for ug
      --------------------------------
   .ve

   RETURN VALUE:
   void *
   .n        NULL pointer                   if error occurs
   .n        'theBlock'                     if OK by type of 'SIMPLE_HEAP'
   .n        '((char *)newBlock)+ALIGNMENT' or
   .n        '((char *)theBlock)+ALIGNMENT' if OK by type of 'GENERAL_HEAP'
   D*/
/****************************************************************************/

void *GetMem (HEAP *theHeap, MEM n, INT mode)
{
  BLOCK *theBlock,*newBlock;
  long newsize,allocated;

  if (theHeap==NULL) return(NULL);              /* there is no heap         */
  if (theHeap->heapptr==NULL) return(NULL);      /* heap is full             */
  if (n==0) return(NULL);

  if (theHeap->type==SIMPLE_HEAP)               /* it's a mark/release heap */
  {
    if (mode==FROM_TOP)
    {
      theBlock = theHeap->heapptr;
      if (theBlock->size<n)
        return (NULL);
      newsize = FLOOR(theBlock->size-n);
      if (newsize<sizeof(BLOCK)) return(NULL);

      theHeap->used += theBlock->size-newsize;
      theBlock->size = newsize;
      return((void *)(((MEM)theBlock)+newsize));
    }
    if (mode!=FROM_BOTTOM) return(NULL);
    theBlock = theHeap->heapptr;
    allocated = CEIL(n);
    /*
        should it be
            if (allocated>theBlock->size) return(NULL);
        since the BLOCK struct is only shifted by allocated?
     */
    if (allocated>theBlock->size-sizeof(BLOCK)) return(NULL);
    newBlock = (BLOCK *) (((MEM)theBlock)+allocated);
    newBlock->size = theBlock->size-allocated;
    theHeap->heapptr = newBlock;
    theHeap->used += allocated;
    return((void *)(theBlock));
  }
  else                                          /* it's a general heap      */
  {
    /* add ALIGNMENT bytes at the beginning of the storage block
       in order to remember size of allocated memory chunk there */
    n += ALIGNMENT;

    /* if n is smaller than sizeof(BLOCK) then DisposeMem will fail:
       there won't be enough room for the next and previous pointer */
    if (n<sizeof(BLOCK)) n=sizeof(BLOCK);

    /* find first BLOCK that fits */
    for (theBlock=theHeap->heapptr; theBlock->next!=theHeap->heapptr;
         theBlock=theBlock->next)
      if (theBlock->size>=n) break;
    if (theBlock->size>=n)
    {
      /* if n is not an aligned size, we have to add some empty
         space at the end of the memory chunk.                  */
      MEM real_n = CEIL(n);

      /* BLOCK found */
      if (theBlock->size >= (CEIL(sizeof(BLOCK))+real_n))
      {
        /* decrease BLOCK size only */
        newsize = FLOOR(theBlock->size-real_n);
        allocated = theBlock->size-newsize;
        theHeap->used += allocated;
        theBlock->size = newsize;
        newBlock = (BLOCK *) (((MEM)theBlock)+newsize);
        newBlock->size = allocated;
        theHeap->heapptr = theBlock;
        return((void *)(((char *)newBlock)+ALIGNMENT));
      }
      else
      {
        /* BLOCK is allocated entirely */
        theHeap->used += theBlock->size;
        if (theBlock->next==theBlock)
        {
          /* this is the last BLOCK */
          theHeap->heapptr = NULL;
        }
        else
        {
          /* remove BLOCK from list of free blocks */
          theBlock->previous->next = theBlock->next;
          theBlock->next->previous = theBlock->previous;
          theHeap->heapptr = theBlock->next;
        }
        return((void *)(((char *)theBlock)+ALIGNMENT));
      }
    }
  }

  return(NULL);
}

void *GetMemUsingKey (HEAP *theHeap, MEM n, INT mode, INT key)
{
  if (theHeap->type==SIMPLE_HEAP)
  {
    if (mode==FROM_TOP)
    {
      if (theHeap->topStackPtr>0)
      {
        if (key>theHeap->topStackPtr)
        {
          /* Mark/Release calls not balanced */
          ASSERT(FALSE);
          return(NULL);
        }
        if (key<theHeap->topStackPtr)
        {
          /* stack pos already released */
          ASSERT(FALSE);
          return(NULL);
        }
        return(GetMem(theHeap,n,mode));
      }
      /* not marked */
      ASSERT(FALSE);
      return(NULL);
    }
    if (mode==FROM_BOTTOM)
    {
      if (theHeap->bottomStackPtr>0)
      {
        if (key>theHeap->bottomStackPtr)
        {
          /* Mark/Release calls not balanced */
          ASSERT(FALSE);
          return(NULL);
        }
        if (key<theHeap->bottomStackPtr)
        {
          /* stack pos already released */
          ASSERT(FALSE);
          return(NULL);
        }
        return(GetMem(theHeap,n,mode));
      }
      /* not marked */
      ASSERT(FALSE);
      return(NULL);
    }
    /* wrong mode */
    ASSERT(FALSE);
    return(NULL);
  }
  /* no key for GENERAL_HEAP */
  return (GetMem(theHeap,n,mode));
}

/****************************************************************************/
/*D
   DisposeMem - Free memory previously allocated from that heap

   SYNOPSIS:
   void DisposeMem (HEAP *theHeap, void *buffer);

   PARAMETERS:
   .  theHeap - heap from which memory has been allocated
   .  buffer - memory area previously allocated

   DESCRIPTION:
   This function creates free memory previously allocated from that heap.
   This function is only valid for a heap of type GENERAL_HEAP.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void DisposeMem (HEAP *theHeap, void *buffer)
{
  BLOCK *newBlock,*theBlock,*nextBlock;
  MEM b,n,p;

  if (theHeap->type!=GENERAL_HEAP) REP_ERR_RETURN_VOID;

  /* reconstruct BLOCK */
  newBlock = (BLOCK *) (((char *)buffer)-ALIGNMENT);
  theHeap->used -= newBlock->size;

  /* simplest case is when free list is empty */
  if (theHeap->heapptr==NULL)
  {
    theHeap->heapptr = newBlock;
    newBlock->next = newBlock;
    newBlock->previous = newBlock;
    return;
  }

  /* find next and previous blocks */
  b = (MEM) newBlock;

  for (theBlock=theHeap->heapptr; theBlock->next!=theHeap->heapptr;
       theBlock=theBlock->next)
  {
    p = (MEM) theBlock;
    n = (MEM) theBlock->next;
    if (n<=p)
    {
      /* wrap around position or single BLOCK */
      if ((b<n)||(b>p)) break;
    }
    else
    {
      /* n>p */
      if ((b>p)&&(b<n)) break;
    }
  }
  p = (MEM) theBlock;
  n = (MEM) theBlock->next;
  if (n<=p)
  {
    /* wrap around position or single BLOCK */
    if (b<n)
    {
      /* new block below lowest free block */
      if (b+newBlock->size==n)
      {
        /* merge newBlock and theBlock->next */
        theBlock = theBlock->next;
        newBlock->size += theBlock->size;
        theHeap->heapptr = newBlock;
        if (theBlock==theBlock->next)
          newBlock->next = newBlock->previous = newBlock;
        else
        {
          newBlock->next = theBlock->next;
          newBlock->previous = theBlock->previous;
          theBlock->previous->next = newBlock;
          theBlock->next->previous = newBlock;
        }
        return;
      }
      else
      {
        /* cannot merge blocks */
        newBlock->next = theBlock->next;
        newBlock->previous = theBlock;
        theBlock->next = newBlock;
        newBlock->next->previous = newBlock;
        return;
      }
    }
    if (b>p)
    {
      /* new block above highest free block */
      theHeap->heapptr = theBlock;
      if (p+theBlock->size==b)
      {
        /* merge the blocks */
        theBlock->size += newBlock->size;
        return;
      }
      else
      {
        /* cannot merge blocks */
        newBlock->next = theBlock->next;
        newBlock->previous = theBlock;
        theBlock->next = newBlock;
        newBlock->next->previous = newBlock;
        return;
      }
    }
  }
  else
  {
    /* n>p */
    if ((b>p)&&(b<n))
    {
      /* try to merge with both blocks */
      nextBlock = theBlock->next;
      if (p+theBlock->size==b)
      {
        /* merge with block below */
        theBlock->size += newBlock->size;
        newBlock = theBlock;
        b = p;
      }
      else
      {
        newBlock->previous = theBlock;
        newBlock->next = theBlock->next;
        theBlock->next = newBlock;
        newBlock->next->previous = newBlock;
      }
      if (b+newBlock->size==n)
      {
        /* merge newBlock and nextBlock */
        newBlock->size += nextBlock->size;
        theHeap->heapptr = newBlock;

        /* remove nextBlock from block list */
        nextBlock->previous->next = nextBlock->next;
        nextBlock->next->previous = nextBlock->previous;
      }
      return;
    }
  }

  /* here must be something wrong */
  assert(FALSE);
  return;
}

/****************************************************************************/
/*D
   GetFreelistMemory - Get an object from free list if possible

   SYNOPSIS:
   void *GetFreelistMemory (HEAP *theHeap, INT size);

   PARAMETERS:
   .  theHeap - pointer to Heap
   .  size - size of the object
   .  type - type of the requested object

   DESCRIPTION:
   This function gets an object of type `type` from free list if possible,
   otherwise it allocates memory from the heap using 'GetMem'.

   RETURN VALUE:
   void *
   .n   pointer to an object of the requested type
   .n   NULL if object of requested type is not available
   D*/
/****************************************************************************/

void *GetFreelistMemory (HEAP *theHeap, INT size)
{
  void **ptr, *obj;
  INT i,j,k,l;

  if (size == 0)
    return(NULL);
  obj = NULL;

  /* 'ptr' will be set equal to 'theHeap->freeObjects[k]' but with	        */
  /* different interpretation: void ** instead of void *. 'ptr'			*/
  /* points to the first two bytes of the object (i.e. unsigned INT ctrl	*/
  /* and INT id) but will be interpreted as a void * pointer, witch points*/
  /* to the next free object.                                                                                   */

  i = (size / ALIGNMENT);
  for (j=0; j<MAXFREEOBJECTS; j++)
  {
    k = (i + j) % MAXFREEOBJECTS;
    l = theHeap->SizeOfFreeObjects[k];
    if (l == size)
    {
      if (theHeap->freeObjects[k] != NULL)
      {
        ptr = (void **) theHeap->freeObjects[k];
        theHeap->freeObjects[k] = ptr[0];
        obj = (void *) ptr;
        theHeap->freelistmem -= size;
      }
      break;
    }
    if(l == -1)
      break;
  }

  if (obj == NULL)
    obj = GetMem(theHeap,size,FROM_BOTTOM);

  if (obj == NULL)
  {
    printf( "ERROR in low/heaps.c/GetFreelistMemory: not enough memory for %d bytes.\n", size );
    fprintf( stderr, "ERROR in low/heaps.c/GetFreelistMemory: not enough memory for %d bytes.\n", size );
    assert(0);
  }

  if (obj != NULL)
    memset(obj,0,size);

  return(obj);
}

/****************************************************************************/
/*D
   PutFreelistMemory - Put an object in the free list

   SYNOPSIS:
   INT PutFreelistMemory (HEAP *theHeap, void *object, INT size);

   PARAMETERS:
   .  theHeap - pointer to Heap
   .  object - object to insert in free list
   .  size - size of the object

   DESCRIPTION:
   This function puts an object in the free list.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 when error occured.
   D*/
/****************************************************************************/

INT PutFreelistMemory (HEAP *theHeap, void *object, INT size)
{
  void **ptr;
  INT i,j,k,l;

  memset(object,0,size);
#ifdef Debug
  SETHEAPFAULT(object,-1);              /* this is to check heap faults */
#endif
  ptr = (void **) object;

  /* 'ptr' will be set equal to 'object' but with different inter-		*/
  /* pretation: void ** instead of void *. 'ptr' points to the first		*/
  /* two bytes of the object (i.e. unsigned INT ctrl	and INT id) but         */
  /* will be interpreted as a void * pointer, witch will be set equal   */
  /* to 'theHeap->freeObjects[k]' i.e. the first free object.			    */

  i = (size / ALIGNMENT);
  for (j=0; j<MAXFREEOBJECTS; j++)
  {
    k = (i + j) % MAXFREEOBJECTS;
    l = theHeap->SizeOfFreeObjects[k];
    if (l == size)
    {
      ptr[0] = theHeap->freeObjects[k];
      theHeap->freeObjects[k] = object;
      theHeap->freelistmem += size;
      return(0);
    }
    if(l == -1)
    {
      theHeap->SizeOfFreeObjects[k] = size;
      ptr[0] = theHeap->freeObjects[k];
      theHeap->freeObjects[k] = object;
      theHeap->freelistmem += size;
      return(0);
    }
  }

  /* MAXFREEOBJECTS to small! */

  RETURN(1);
}

/****************************************************************************/
/*D
   Mark - Mark heap position for future release

   SYNOPSIS:
   INT Mark (HEAP *theHeap, INT mode, INT *key);

   PARAMETERS:
   .  theHeap - heap to mark
   .  mode - 'FROM_TOP' or 'FROM_BOTTOM' of the block

   DESCRIPTION:
   This function marks heap position for future release. Only valid in
   the 'SIMPLE_HEAP' type.

   RETURN VALUE:
   INT
   .n     0 if OK
   .n     1 if mark stack full or wrong heap type
   D*/
/****************************************************************************/

INT Mark (HEAP *theHeap, INT mode, INT *key)
{
  if (theHeap->type!=SIMPLE_HEAP) return(1);

  if (mode==FROM_TOP)
  {
    if (theHeap->topStackPtr<MARK_STACK_SIZE)
    {
      theHeap->topStack[theHeap->topStackPtr++] =
        ((MEM)theHeap->heapptr) + ((MEM)theHeap->heapptr->size);
      *key = theHeap->topStackPtr;
      return(0);
    }
  }
  if (mode==FROM_BOTTOM)
  {
    if (theHeap->bottomStackPtr<MARK_STACK_SIZE)
    {
      theHeap->bottomStack[theHeap->bottomStackPtr++] =
        ((MEM)theHeap->heapptr);
      *key = theHeap->bottomStackPtr;
      return(0);
    }
  }
  return(1);
}

/****************************************************************************/
/*D
   Release - Release to next stack position

   SYNOPSIS:
   INT Release (HEAP *theHeap, INT mode, INT key);

   PARAMETERS:
   .  theHeap - heap to release
   .  mode - 'FROM_TOP' or 'FROM_BOTTOM' of the block

   DESCRIPTION:
   This function releases to the next stack position. Only valid in the
   'SIMPLE_HEAP' type.

   RETURN VALUE:
   INT
   .n     0 if OK
   .n     1 if mark stack empty or wrong heap type.
   D*/
/****************************************************************************/

INT Release (HEAP *theHeap, INT mode, INT key)
{
  MEM oldsize;
  MEM newsize;

  if (theHeap->type!=SIMPLE_HEAP) return(1);

  if (mode==FROM_TOP)
  {
    if (theHeap->topStackPtr>0)
    {
      if (key>theHeap->topStackPtr)
      {
        /* Mark/Release calls not balanced */
        ASSERT(FALSE);
        return(1);
      }
      if (key<theHeap->topStackPtr)
      {
        /* stack pos already released */
        ASSERT(FALSE);
        return(2);
      }
      oldsize = theHeap->heapptr->size;
      newsize = theHeap->topStack[--theHeap->topStackPtr]-((MEM)theHeap->heapptr);
      theHeap->heapptr->size = newsize;
      theHeap->used -= newsize-oldsize;
      return(0);
    }
  }
  if (mode==FROM_BOTTOM)
  {
    if (theHeap->bottomStackPtr>0)
    {
      if (key>theHeap->bottomStackPtr)
      {
        /* Mark/Release calls not balanced */
        ASSERT(FALSE);
        return(3);
      }
      if (key<theHeap->bottomStackPtr)
      {
        /* stack pos already released */
        ASSERT(FALSE);
        return(4);
      }
      oldsize = theHeap->heapptr->size;
      newsize = (((MEM)theHeap->heapptr)+((MEM)theHeap->heapptr->size))
                -theHeap->bottomStack[--theHeap->bottomStackPtr];
      theHeap->heapptr = (BLOCK *) theHeap->bottomStack[theHeap->bottomStackPtr];
      theHeap->heapptr->size = newsize;
      theHeap->used -= newsize-oldsize;
      return(0);
    }
  }
  return(5);
}

/****************************************************************************/
/*D
   HeapSize - Get heap size

   SYNOPSIS:
   MEM HeapSize (const HEAP *theHeap)

   PARAMETERS:
   .  theHeap - heap to get heap size

   DESCRIPTION:
   This function gets the heap size.

   RETURN VALUE:
   MEM
   .n    theHeap->size
   D*/
/****************************************************************************/

MEM HeapSize (const HEAP *theHeap)
{
  return(theHeap->size);
}

/****************************************************************************/
/*D
   HeapUsed - Get used memory of heap

   SYNOPSIS:
   MEM HeapUsed (const HEAP *theHeap)

   PARAMETERS:
   .  theHeap - heap to get used memory of heap

   DESCRIPTION:
   This function gets the used memory of heap.

   RETURN VALUE:
   MEM
   .n    theHeap->used
   D*/
/****************************************************************************/

MEM HeapUsed (const HEAP *theHeap)
{
  return(theHeap->used);
}

/****************************************************************************/
/*D
   HeapFree - Get free memory of heap (without free lists)

   SYNOPSIS:
   MEM HeapFree (const HEAP *theHeap)

   PARAMETERS:
   .  theHeap - heap to get free memory of heap

   DESCRIPTION:
   This function gets the free memory of heap. The free momory in the free lists is
   not taken into account

   RETURN VALUE:
   MEM
   .n    theHeap->size-theHeap->used
   D*/
/****************************************************************************/

MEM HeapFree (const HEAP *theHeap)
{
  return(theHeap->size-theHeap->used);
}

/****************************************************************************/
/*D
   HeapFreelistUsed - Get memory of heap in freelists

   SYNOPSIS:
   MEM HeapFreelistUsed (const HEAP *theHeap)

   PARAMETERS:
   .  theHeap - heap to get used memory of heap

   DESCRIPTION:
   This function gets the used memory of heap which is available in the freelists.

   RETURN VALUE:
   MEM
   .n    theHeap->freelistmem
   D*/
/****************************************************************************/

MEM HeapFreelistUsed (const HEAP *theHeap)
{
  return(theHeap->freelistmem);
}

/****************************************************************************/
/*D
   HeapTotalFree - Get memory of heap in freelists

   SYNOPSIS:
   MEM HeapTotalFree (const HEAP *theHeap)

   PARAMETERS:
   .  theHeap - heap to get used memory of heap

   DESCRIPTION:
   This function gets the used memory of heap which is available in the freelists.

   RETURN VALUE:
   MEM
   .n    theHeap->size-theHeap->used+theHeap->freelistmem
   D*/
/****************************************************************************/

MEM HeapTotalFree (const HEAP *theHeap)
{
  return(theHeap->size-theHeap->used+theHeap->freelistmem);
}

/****************************************************************************/
/*D
   InitVirtualHeapManagement - Initialize the VIRT_HEAP_MGMT data structure

   SYNOPSIS:
   INT InitVirtualHeapManagement (VIRT_HEAP_MGMT *theVHM, MEM TotalSize);

   PARAMETERS:
   .  theVHM - pointer to the storage to initialize
   .  TotalSize - the total size of the heap to manage

   DESCRIPTION:
   This function initializes the VIRT_HEAP_MGMT data structure that provides
   additional memory independently of the 'SIMPLE_HEAP' and 'GENERAL_HEAP' for
   further use by handling virtual heaps.

   RETURN VALUE:
   INT
   .n    'BHM_OK' if OK
   .n          99 if error occurred.
   D*/
/****************************************************************************/

INT InitVirtualHeapManagement (VIRT_HEAP_MGMT *theVHM, MEM TotalSize)
{
  if (theVHM==NULL)
    return (99);

  /* first clear everything */
  memset(theVHM,0,sizeof(VIRT_HEAP_MGMT));

  /* now init what is necessary */
  if (TotalSize==SIZE_UNKNOWN)
    theVHM->locked    = FALSE;
  else
    theVHM->locked    = TRUE;

  theVHM->TotalSize    = TotalSize;
  theVHM->TotalUsed    = 0;
  theVHM->UsedBlocks   = 0;
  theVHM->LargestGap   = 0;
  theVHM->nGaps        = 0;

  return (BHM_OK);
}

/****************************************************************************/
/*D
   CalcAndFixTotalSize - Sum up the sizes of the blocks, set 'TotalSize' and
   lock it

   SYNOPSIS:
   MEM CalcAndFixTotalSize (VIRT_HEAP_MGMT *theVHM)

   PARAMETERS:
   .  theVHM - pointer to the storage to init

   DESCRIPTION:
   This function should be called while initializing virtual heaps. It
   sums up the sizes of all the heaps and returns the 'TotalSize' of all.

   RETURN VALUE:
   MEM
   .n    'theVHM->TotalSize' if OK
   .n                     99 if error occurred.
   D*/
/****************************************************************************/

MEM CalcAndFixTotalSize (VIRT_HEAP_MGMT *theVHM)
{
  if (theVHM==NULL)
    return (99);

  assert(theVHM->locked!=TRUE);

  theVHM->locked        = TRUE;
  theVHM->TotalSize    = theVHM->TotalUsed;
  theVHM->LargestGap    = 0;
  theVHM->nGaps        = 0;

  return (theVHM->TotalSize);
}

/****************************************************************************/
/*D
   GetNewBlockID - Return a unique block ID starting with FIRST_BLOCK_ID

   SYNOPSIS:
   BLOCK_ID GetNewBlockID ()

   PARAMETERS:
   void

   DESCRIPTION:
   This function returns a unique block ID starting with the FIRST_BLOCK_ID.

   RETURN VALUE:
   BLOCK_ID
   .n    'newID'
   D*/
/****************************************************************************/

BLOCK_ID GetNewBlockID ()
{
  static BLOCK_ID newID = 0;

  return (++newID);
}

/****************************************************************************/
/*D
   GetBlockDesc - Return a pointer to the block descriptor with 'id'

   SYNOPSIS:
   BLOCK_DESC *GetBlockDesc (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id)

   PARAMETERS:
   .  theVHM - pointer to the virtual heap management
   .  id - id of the desired block

   DESCRIPTION:
   As the location of the block descriptors is not fixed in the heap
   management 'GetBlockDesc' returns the address of the block descriptor.

   RETURN VALUE:
   BLOCK_DESC *
   .n    NULL if not defined in theVHM
   .n    'theVHM->BlockDesc' pointer to block descriptor

   D*/
/****************************************************************************/

BLOCK_DESC *GetBlockDesc (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id)
{
  INT i;

  if (theVHM==NULL)
    return (NULL);

  for (i=0; i<theVHM->UsedBlocks; i++)
    if (B_ID(theVHM,i)==id)
      break;

  if (i<theVHM->UsedBlocks)
    return (&(theVHM->BlockDesc[i]));
  else
    return (NULL);
}

/****************************************************************************/
/*D
   DefineBlock - Set size and offset of a new block

   SYNOPSIS:
   INT DefineBlock (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id, MEM size);

   PARAMETERS:
   .  theVHM - pointer to the virtual heap management
   .  id - id of the block to define
   .  size - size to be allocated

   DESCRIPTION:
   This function sets size and offset of a new block. It tries to fill gaps.

   RETURN VALUE:
   INT
   .n    'BHM_OK'        if OK
   .n    'HEAP_FULL'     if heap is full
   .n    'BLOCK_DEFINED' if block is already defined
   .n    'NO_FREE_BLOCK' if number of 'MAXNBLOCKS' reached
   .n     99             if error occurred.
   D*/
/****************************************************************************/

INT DefineBlock (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id, MEM size)
{
  BLOCK_DESC *theBlock;
  MEM Gap,BestFitGap,LargestGap;
  INT i,BestFitGapPos;

  if (theVHM==NULL)
    return (99);

  size = CEIL(size);

  /* check size */
  if (theVHM->TotalSize!=SIZE_UNKNOWN)
    if (size>theVHM->TotalSize-theVHM->TotalUsed)
      return (HEAP_FULL);

  theBlock = GetBlockDesc(theVHM,id);

  if (theBlock!=NULL)
    /* block already is defined */
    return (BLOCK_DEFINED);

  if (theVHM->UsedBlocks>=MAXNBLOCKS)
    return (NO_FREE_BLOCK);

  if (theVHM->TotalSize==SIZE_UNKNOWN)
  {
    /* the size is unbounded: take first new block */
    i = theVHM->UsedBlocks;

    theVHM->TotalUsed    += size;
    theVHM->UsedBlocks    ++;

    B_ID(theVHM,i)        = id;
    B_SIZE(theVHM,i)    = size;
    B_OFFSET(theVHM,i)    = CALC_B_OFFSET(theVHM,i);

    return (BHM_OK);
  }

  /* the TotalSize is fixed */

  /* is there a gap large enough */
  if (theVHM->nGaps>0)
    if (size<theVHM->LargestGap)
    {
      /*    find the minimal gap large enough */
      BestFitGap = theVHM->LargestGap;
      Gap = B_OFFSET(theVHM,0);
      if ((Gap>=size) && (Gap<BestFitGap))
      {
        BestFitGap      = Gap;
        BestFitGapPos = 0;
      }
      for (i=1; i<theVHM->UsedBlocks; i++)
      {
        Gap = B_OFFSET(theVHM,i) - (B_OFFSET(theVHM,i-1) + B_SIZE(theVHM,i-1));
        if ((Gap>=size) && (Gap<BestFitGap))
        {
          BestFitGap      = Gap;
          BestFitGapPos = i;
        }
      }

      /* shift the descriptors one up */
      for (i=theVHM->UsedBlocks-1; i>BestFitGapPos; i--)
        theVHM->BlockDesc[i] = theVHM->BlockDesc[i-1];

      theVHM->nGaps--;

      theVHM->TotalUsed    += size;
      theVHM->UsedBlocks    ++;

      B_ID(theVHM,BestFitGapPos)        = id;
      B_SIZE(theVHM,BestFitGapPos)    = size;
      B_OFFSET(theVHM,BestFitGapPos)    = CALC_B_OFFSET(theVHM,BestFitGapPos);

      /* recalculate LargestGap? */
      if (BestFitGap==theVHM->LargestGap)
      {
        LargestGap = 0;
        for (i=0; i<theVHM->TotalUsed; i++)
          if (LargestGap<B_SIZE(theVHM,i))
            LargestGap = theVHM->BlockDesc[i].size;

        theVHM->LargestGap = LargestGap;
      }

      return (BHM_OK);
    }

  /* there is no gap large enough: take the next new block */
  i = theVHM->UsedBlocks;

  theVHM->TotalUsed    += size;
  theVHM->UsedBlocks    ++;

  B_ID(theVHM,i)        = id;
  B_SIZE(theVHM,i)    = size;
  B_OFFSET(theVHM,i)    = CALC_B_OFFSET(theVHM,i);

  return (BHM_OK);
}

/****************************************************************************/
/*D
   FreeBlock - Free a block in the bhm defined before

   SYNOPSIS:
   INT FreeBlock (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id);

   PARAMETERS:
   .  theVHM - pointer to the virtual heap management
   .  id - id of the block to free

   DESCRIPTION:
   This function frees a block in the bhm defined before.

   RETURN VALUE:
   INT
   .n    'BHM_OK'            if OK
   .n    'BLOCK_NOT_DEFINED' if block is not defined, nothing to free
   .n     99                 if error occurred.
   D*/
/****************************************************************************/

INT FreeBlock (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id)
{
  MEM NewGap;
  INT i,i_free;

  if (theVHM==NULL)
    return (99);

  for (i_free=0; i_free<theVHM->UsedBlocks; i_free++)
    if (B_ID(theVHM,i_free)==id)
      break;

  if (i_free>=theVHM->UsedBlocks)
    /* block is not defined, nothing to free */
    return (BLOCK_NOT_DEFINED);

  assert(theVHM->TotalUsed > B_SIZE(theVHM,i_free));

  theVHM->UsedBlocks--;
  theVHM->TotalUsed -= B_SIZE(theVHM,i_free);

  if (theVHM->TotalSize==SIZE_UNKNOWN)
  {
    /* shift the blocks one down and recalculate the offset */
    for (i=i_free; i<theVHM->UsedBlocks; i++)
    {
      theVHM->BlockDesc[i] = theVHM->BlockDesc[i+1];
      B_OFFSET(theVHM,i)     = CALC_B_OFFSET(theVHM,i);
    }

    return (BHM_OK);
  }

  /* the TotalSize is fixed */

  /* shift the blocks one down (don't change the offset!) */
  for (i=i_free; i<theVHM->UsedBlocks; i++)
    theVHM->BlockDesc[i] = theVHM->BlockDesc[i+1];

  /* new gap? */
  if (i_free<theVHM->UsedBlocks)
  {
    theVHM->nGaps++;

    NewGap = B_OFFSET(theVHM,i_free) - (B_OFFSET(theVHM,i_free-1) + B_SIZE(theVHM,i_free-1));
    if (theVHM->LargestGap<NewGap)
      theVHM->LargestGap = NewGap;
  }

  return (BHM_OK);
}
