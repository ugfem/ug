// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      heaps.h                                                       */
/*                                                                          */
/* Purpose:   low-level memory management for ug                            */
/*                                                                          */
/* Author:      Peter Bastian                                               */
/*              Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen */
/*              Universitaet Heidelberg                                     */
/*              Im Neuenheimer Feld 368                                     */
/*              6900 Heidelberg                                             */
/*                                                                          */
/* History:   29.01.92 begin, ug version 2.0                                */
/*                                                                          */
/* Revision:  04.09.95                                                      */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __HEAPS__
#define __HEAPS__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#if defined __NECSX4__ && defined _MALLOC64
#define MEM_SIZE_ULL
#include "stdlib.h" /* for the patched malloc */
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/* defines for the simple and general heap management                       */
/****************************************************************************/

#define MIN_HEAP_SIZE   256              /* smallest heap to allocate       */
#define MARK_STACK_SIZE  20              /* max depth of mark/release calls */

#define GENERAL_HEAP      0              /* heap with alloc/free mechanism  */
#define SIMPLE_HEAP       1              /* heap with mark/release mechanism*/

#define FROM_TOP          1              /* allocate from top of stack      */
#define FROM_BOTTOM       2              /* allocate from bottom of stack   */

#define MAXFREEOBJECTS  128                     /* number of free object pionters  */

/* by convention, tempory memory on a simple heap should allocated FROM_TOP */
/* the Freelist memory is allocated FROM_BOTTOM                             */

#define MarkTmpMem(p,kp)     Mark(p,FROM_TOP,kp)
#define GetTmpMem(p,n,k)         GetMemUsingKey(p,n,FROM_TOP,k)
#define ReleaseTmpMem(p,k)       Release(p,FROM_TOP,k)

/****************************************************************************/
/* defines and macros for the virtual heap management                       */
/****************************************************************************/

#define MAXNBLOCKS         50        /* that many blocks can be allocated   */
#define SIZE_UNKNOWN        0        /* pass to init routine if no heap yet */
#define SIZEOF_VHM            sizeof(VIRT_HEAP_MGMT)
/* the memory sized neded for the vhm  */

#define BHM_OK              0        /* ok return code for virtual heap mgmt*/

/* return codes of DefineBlock */
#define HEAP_FULL            1       /* return code if storage exhausted    */
#define BLOCK_DEFINED        2       /* return code if block already defined*/
#define NO_FREE_BLOCK        3       /* return code if no free block found  */

/* return codes of FreeBlock */
#define BLOCK_NOT_DEFINED    1       /* return code if the block is not def */

/* some useful macros */
#define OFFSET_IN_HEAP(vhm,id)  (GetBlockDesc((VIRT_HEAP_MGMT*)vhm,id).offset)
#define TOTUSED_IN_HEAP(vhm)    ((vhm).TotalUsed)
#define IS_BLOCK_DEFINED(vhm,id) (GetBlockDesc((VIRT_HEAP_MGMT*)vhm,id)!=NULL)

#define CEIL(n)          ((n)+((ALIGNMENT-((n)&(ALIGNMENT-1)))&(ALIGNMENT-1)))

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

#ifdef MEM_SIZE_ULL
typedef unsigned long long MEM;
#else
typedef unsigned long MEM;
#endif

/****************************************************************************/
/* structs and typedefs for the simple and general heap management          */
/****************************************************************************/

struct block {
  MEM size;
  struct block *next,*previous;
};

typedef struct {
  INT type;
  MEM size;
  MEM used;
  MEM freelistmem;
  struct block *heapptr;
  INT topStackPtr,bottomStackPtr;
  MEM topStack[MARK_STACK_SIZE];
  MEM bottomStack[MARK_STACK_SIZE];
  INT SizeOfFreeObjects[MAXFREEOBJECTS];
  void *freeObjects[MAXFREEOBJECTS];
        #ifdef Debug
  INT objcount[MAXFREEOBJECTS];
        #endif
} HEAP;

/****************************************************************************/
/* structs and typedefs for the block virtual management                    */
/****************************************************************************/

typedef struct {

  INT id;                           /* id for this block                    */
  MEM offset;                       /* offset of the data in the heap       */
  MEM size;                         /* size of the data in the heap         */

} BLOCK_DESC;

typedef struct {

  INT locked;                       /* if TRUE the TotalSize is fixed        */
  MEM TotalSize;                    /* total size of the associated heap     */
  MEM TotalUsed;                    /* total size used                       */
  INT UsedBlocks;                   /* number of blocks initialized          */
  INT nGaps;                        /* TRUE if a gap between exist. blocks   */
  MEM LargestGap;                   /* largest free gap between blocks       */
  BLOCK_DESC BlockDesc[MAXNBLOCKS];
  /* the different block descriptors       */
} VIRT_HEAP_MGMT;

/****************************************************************************/
/* typedefs for the block virtual management                                */
/****************************************************************************/

typedef INT BLOCK_ID;
typedef struct block BLOCK;

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT          InitHeaps                (void);

/* functions for the simple and general heap management */
HEAP        *NewHeap                (INT type, MEM size, void *buffer);
void        *GetMem                 (HEAP *theHeap, MEM n, INT mode);
void            *GetMemUsingKey                 (HEAP *theHeap, MEM n, INT mode, INT key);
void         DisposeMem             (HEAP *theHeap, void *buffer);

void        *GetFreelistMemory      (HEAP *theHeap, INT size);
INT          PutFreelistMemory      (HEAP *theHeap, void *object, INT size);

INT          Mark                   (HEAP *theHeap, INT mode, INT *key);
INT          Release                (HEAP *theHeap, INT mode, INT key);

MEM          HeapSize               (const HEAP *theHeap);
MEM          HeapUsed               (const HEAP *theHeap);
MEM                      HeapFree                               (const HEAP *theHeap);
MEM          HeapFreelistUsed       (const HEAP *theHeap);
MEM                      HeapTotalFree                  (const HEAP *theHeap);

/* functions for the virtual heap management */
INT          InitVirtualHeapManagement(VIRT_HEAP_MGMT *theVHM, MEM TotalSize);
MEM          CalcAndFixTotalSize    (VIRT_HEAP_MGMT *theVHM);
BLOCK_ID     GetNewBlockID            (void);
BLOCK_DESC  *GetBlockDesc            (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id);
INT          DefineBlock            (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id, MEM size);
INT          FreeBlock                (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id);

#endif
