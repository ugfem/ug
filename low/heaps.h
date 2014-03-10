// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file heaps.h
 * \ingroup low
 */

/** \addtogroup low
 *
 * @{
 */

/****************************************************************************/
/*                                                                          */
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

/* Only for UG_USE_SYSTEM_HEAP */
#include <vector>

#include "ugtypes.h"
#include "namespace.h"

START_UG_NAMESPACE

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

/** \brief Smallest heap to allocate       */
#define MIN_HEAP_SIZE   256
/** \brief Max depth of mark/release calls */
#define MARK_STACK_SIZE 128

enum HeapType {GENERAL_HEAP,                  /**< Heap with alloc/free mechanism  */
               SIMPLE_HEAP         /**< Heap with mark/release mechanism*/
};

enum HeapAllocMode
{FROM_TOP=1,                       /**< Allocate from top of stack      */
 FROM_BOTTOM=2                       /**< Allocate from bottom of stack   */
};

/** \brief Number of free object pointers  */
#define MAXFREEOBJECTS  128

/* by convention, tempory memory on a simple heap should allocated FROM_TOP */
/* the Freelist memory is allocated FROM_BOTTOM                             */

#define MarkTmpMem(p,kp)     Mark(p,FROM_TOP,kp)
#define GetTmpMem(p,n,k)         GetMemUsingKey(p,n,FROM_TOP,k)
#define ReleaseTmpMem(p,k)       Release(p,FROM_TOP,k)

/****************************************************************************/
/****************************************************************************/
/** @name Defines and macros for the virtual heap management                 */

/** \brief That many blocks can be allocated   */
#define MAXNBLOCKS         50

/** \brief Pass to init routine if no heap yet */
#define SIZE_UNKNOWN        0

/** \brief The memory sized neded for the vhm  */
#define SIZEOF_VHM            sizeof(VIRT_HEAP_MGMT)

/** \brief Ok return code for virtual heap mgmt*/
#define BHM_OK              0

/** \brief Return codes of DefineBlock */
enum {HEAP_FULL =           1,           /**< Return code if storage exhausted    */
      BLOCK_DEFINED =       2,           /**< Return code if block already defined*/
      NO_FREE_BLOCK =       3           /**< Return code if no free block found  */
};

/* return codes of FreeBlock */
/** \brief Return code if the block is not defined */
#define BLOCK_NOT_DEFINED    1

/* some useful macros */
#define OFFSET_IN_HEAP(vhm,id)  (GetBlockDesc((VIRT_HEAP_MGMT*)vhm,id).offset)
#define TOTUSED_IN_HEAP(vhm)    ((vhm).TotalUsed)
#define IS_BLOCK_DEFINED(vhm,id) (GetBlockDesc((VIRT_HEAP_MGMT*)vhm,id)!=NULL)

/* @} */
/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef unsigned long MEM;

/****************************************************************************/
/* structs and typedefs for the simple and general heap management          */
/****************************************************************************/

struct block {
  MEM size;
  struct block *next,*previous;
};

typedef struct {
  enum HeapType type;
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


  /* This is used only if UG_USE_SYSTEM_HEAP is set, but I don't want the
   * #ifdef in an installed header, hence the data member is there all the time. */
  std::vector<void*> markedMemory[MARK_STACK_SIZE];
} HEAP;

/****************************************************************************/
/* structs and typedefs for the block virtual management                    */
/****************************************************************************/

typedef struct {

  INT id;                           /*!< Id for this block                    */
  MEM offset;                       /*!< Offset of the data in the heap       */
  MEM size;                         /*!< Size of the data in the heap         */

} BLOCK_DESC;

typedef struct {

  INT locked;                       /**< If TRUE the TotalSize is fixed        */
  MEM TotalSize;                    /**< Total size of the associated heap     */
  MEM TotalUsed;                    /**< Total size used                       */
  INT UsedBlocks;                   /**< Number of blocks initialized          */
  INT nGaps;                        /**< TRUE if a gap between exist. blocks   */
  MEM LargestGap;                   /**< Largest free gap between blocks       */
  BLOCK_DESC BlockDesc[MAXNBLOCKS];
  /**< The different block descriptors       */
} VIRT_HEAP_MGMT;

/****************************************************************************/
/* typedefs for the block virtual management                                */
/****************************************************************************/

typedef INT BLOCK_ID;
typedef struct block BLOCK;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

#if defined(DYNAMIC_MEMORY_ALLOCMODEL) && defined(Debug)
extern INT check_of_getcallstack;
extern INT check_of_putcallstack;
#endif

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT          InitHeaps                (void);

/** @name Functions for the simple and general heap management */
/* @{ */
HEAP        *NewHeap                (enum HeapType type, MEM size, void *buffer);
void         DisposeHeap            (HEAP *theHeap);
void        *GetMem                 (HEAP *theHeap, MEM n, HeapAllocMode mode);
void            *GetMemUsingKey                 (HEAP *theHeap, MEM n, HeapAllocMode mode, INT key);
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
void             HeapStat                               (const HEAP *theHeap);
/* @} */

/** @name Functions for the virtual heap management */
/* @{ */
INT          InitVirtualHeapManagement(VIRT_HEAP_MGMT *theVHM, MEM TotalSize);
MEM          CalcAndFixTotalSize    (VIRT_HEAP_MGMT *theVHM);
BLOCK_ID     GetNewBlockID            (void);
BLOCK_DESC  *GetBlockDesc            (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id);
INT          DefineBlock            (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id, MEM size);
INT          FreeBlock                (VIRT_HEAP_MGMT *theVHM, BLOCK_ID id);
/* @} */

END_UG_NAMESPACE

/** @} */

#endif
