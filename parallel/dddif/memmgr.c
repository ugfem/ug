// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      memmgr.c                                                      */
/*                                                                          */
/* Purpose:   memory management module                                      */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/04/27 kb  begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

/* standard C library */
#include <stdlib.h>
#include <stdio.h>

#include "compiler.h"
#include "heaps.h"
#include "misc.h"

#include "parallel.h"
#include "general.h"

#ifdef DYNAMIC_MEMORY_ALLOCMODEL
#include "ugm.h"
#endif

/****************************************************************************/

/* define this to protocol all alloc/free requests via hashtable */
/*
   #define WITH_HASH_CONTROL
 */

/* define this to map all PMEM, AMEM and TMEM requests to a UG general heap */
/*
   #define WITH_GENERAL_HEAP
 */


/* define this to detect size of allocatable memory */
/*
   #define DETECT_MEMORY_SIZE
 */




#define HARD_EXIT assert(0)
/*#define HARD_EXIT exit(1)*/


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#ifdef WITH_GENERAL_HEAP
/* constants for UG general heap */
#define HEAP_SIZE     3*1024*1024
#endif


#ifdef WITH_HASH_CONTROL
/* constants for hashing of alloc/free requests (for debugging) */
#define HASHTAB_SIZE  15731    /* prime number not near 2^n */
#define HASH_FUNC(k)   ((k)%HASHTAB_SIZE)
#endif



/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/


#ifdef WITH_HASH_CONTROL
typedef struct _HASH_ENTRY
{
  void   *ptr;           /* hashed key:  pointer to memory block */
  size_t size;           /* hashed data: size of memory block    */
  char info;             /* info character, one of { P, A, T }   */

  int flags;

  struct _HASH_ENTRY *next;

} HASH_ENTRY;
#endif



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


#ifdef GENERAL_HEAP
static HEAP *myheap;
#endif

static INT allocated=0;
static size_t pmem=0;
static size_t amem=0;
static size_t tmem=0;

static size_t mem_from_ug_freelists=0;


#ifdef WITH_HASH_CONTROL
/* hashing of alloc/free requests: hashtable and allocated entries */
/* (from ddd/memmgrs/memmgr_ctrl.c)                   */
HASH_ENTRY *htab[HASHTAB_SIZE];
int nHashEntries;
#endif


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


#ifdef WITH_HASH_CONTROL

/* auxiliary routines for hashing alloc/free requests */
/* (from ddd/memmgrs/memmgr_ctrl.c)                   */


/****************************************************************************/
/*
   NewHashEntry -

   SYNOPSIS:
   static HASH_ENTRY *NewHashEntry (void *ptr, size_t size, char info);

   PARAMETERS:
   .  ptr
   .  size
   .  info

   DESCRIPTION:

   RETURN VALUE:
   HASH_ENTRY
 */
/****************************************************************************/

static HASH_ENTRY *NewHashEntry (void *ptr, size_t size, char info)
{
  HASH_ENTRY *he;

  /*
     printf("%4d: alloc %c %08x %ld\n", me,info,ptr,(unsigned long)size);
   */

  he = (HASH_ENTRY *) malloc(sizeof(HASH_ENTRY));
  he->ptr = ptr;
  he->size = size;
  he->info = info;
  he->flags = 0;
  he->next = NULL;

  nHashEntries++;

  return(he);
}


/****************************************************************************/
/*
   FreeHashEntry -

   SYNOPSIS:
   static void FreeHashEntry (HASH_ENTRY *he);

   PARAMETERS:
   .  he

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void FreeHashEntry (HASH_ENTRY *he)
{
  /*
     printf("%4d: free  %c %08x %ld\n", me,he->info,he->ptr,(unsigned long)he->size);
   */

  nHashEntries--;
  free(he);
}

/****************************************************************************/
/*
   PushHash -

   SYNOPSIS:
   static void PushHash (void *ptr, size_t size, char info);

   PARAMETERS:
   .  ptr
   .  size
   .  info

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void PushHash (void *ptr, size_t size, char info)
{
  unsigned int idx = HASH_FUNC(((unsigned long)ptr));

  if (htab[idx] == NULL)
  {
    /* no collision */
    htab[idx] = NewHashEntry(ptr, size, info);
  }
  else
  {
    /* collision, find entry or none */
    HASH_ENTRY *he;

    for(he=htab[idx]; he->next!=NULL && he->ptr!=ptr; he=he->next)
      ;

    if (he->ptr==ptr)
    {
      UserWriteF("%4d: MEMMGR-ERROR, double alloc at %08x, size %ld (%c, %c)\n",
                 me, ptr, (unsigned long)size,
                 he->info, info);
      HARD_EXIT;
    }

    he->next = NewHashEntry(ptr, size, info);
  }
}


/****************************************************************************/
/*
   PopHash -

   SYNOPSIS:
   static size_t PopHash (void *ptr, char info);

   PARAMETERS:
   .  ptr
   .  info

   DESCRIPTION:

   RETURN VALUE:
   size_t
 */
/****************************************************************************/

static size_t PopHash (void *ptr, char info)
{
  unsigned int idx = HASH_FUNC(((unsigned long)ptr));
  HASH_ENTRY    *he, *helast;

  /* look for entry */
  if (htab[idx] != NULL)
  {
    helast = NULL;
    for(he=htab[idx]; he->next!=NULL && he->ptr!=ptr; he=he->next)
      helast = he;

    if (he->ptr==ptr)
    {
      /* found entry */
      size_t s = he->size;
      if (helast==NULL)
        htab[idx] = he->next;
      else
        helast->next = he->next;

      if (he->info!=info)
      {
        UserWriteF("%4d: MEMMGR-ERROR, wrong free-type %c for alloc %c\n",
                   me, info, he->info);
        HARD_EXIT;
      }

      FreeHashEntry(he);
      return(s);
    }
  }

  UserWriteF("%4d: MEMMGR-ERROR, no alloc for free %c at %08x\n", me, info, ptr);
  HARD_EXIT;

  return(0);       /* never reached */
}


/****************************************************************************/
/*
   HashMarkAll -

   SYNOPSIS:
   static void HashMarkAll (void);

   PARAMETERS:
   .  void

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void HashMarkAll (void)
{
  int i;

  for(i=0; i<HASHTAB_SIZE; i++)
  {
    HASH_ENTRY *he;
    for(he=htab[i]; he!=NULL; he=he->next)
      he->flags = 1;
  }
}


/****************************************************************************/
/*
   HashShowMarks -

   SYNOPSIS:
   static void HashShowMarks (char info);

   PARAMETERS:
   .  info

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void HashShowMarks (char info)
{
  int i;

  for(i=0; i<HASHTAB_SIZE; i++)
  {
    HASH_ENTRY *he;
    for(he=htab[i]; he!=NULL; he=he->next)
    {
      if (he->flags==0 && he->info==info)
      {
        UserWriteF("%4d: MALLOC %c adr=%08x size=%ld\n",
                   me, he->info, he->ptr, (unsigned long) he->size);
      }
    }
  }
}


#endif

/****************************************************************************/


/****************************************************************************/
/*
   memmgr_Report -

   SYNOPSIS:
   void memmgr_Report (void);

   PARAMETERS:
   .  void

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void memmgr_Report (void)
{
        #ifdef WITH_HASH_CONTROL
  UserWriteF("%04d memmgr_Report.  P=%9ld   A=%9ld   T=%9ld    SUM=%9ld\n",
             me, (long)pmem, (long)amem, (long)tmem, (long)allocated);
        #endif

        #ifdef WITH_HASH_CONTROL
  /* HashShowMarks('P'); */
  /* HashShowMarks('A'); */
  HashShowMarks('T');
        #endif

  UserWriteF("%04d memmgr_Report.  Memory from UG's freelists: %9ld\n",
             me, mem_from_ug_freelists);

  fflush(stdout);
}

/****************************************************************************/


/****************************************************************************/
/*
   memmgr_AllocOMEM -

   SYNOPSIS:
   void *memmgr_AllocOMEM (size_t size, int ddd_type, int prio, int attr);

   PARAMETERS:
   .  size
   .  ddd_type
   .  prio
   .  attr

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void *memmgr_AllocOMEM (size_t size, int ddd_type, int prio, int attr)
{
  void   *buffer;

        #ifndef DYNAMIC_MEMORY_ALLOCMODEL
  buffer = GetFreelistMemory(MGHEAP(dddctrl.currMG), size);
        #else
  buffer = GetMemoryForObject(dddctrl.currMG,size,MAOBJ);
        #endif

  /*
     printf("%4d: memmgr_AllocOMem: size=%05d ddd_type=%02d prio=%d attr=%d\n",
     me,size,ddd_type,prio,attr);
   */

  return(buffer);
}


/****************************************************************************/
/*
   memmgr_FreeOMEM -

   SYNOPSIS:
   void memmgr_FreeOMEM (void *buffer, size_t size, int ddd_type);

   PARAMETERS:
   .  buffer
   .  size
   .  ddd_type

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void memmgr_FreeOMEM (void *buffer, size_t size, int ddd_type)
{
  /*
     printf("%d: memmgr_FreeOMEM(): buffer=%x, ddd_type=%d\n", me, buffer, ddd_type);
   */

        #ifndef DYNAMIC_MEMORY_ALLOCMODEL
  PutFreelistMemory(MGHEAP(dddctrl.currMG), buffer, size);
        #else
  PutFreeObject(dddctrl.currMG,buffer,size,MAOBJ);
        #endif
}


/****************************************************************************/
/*
   memmgr_AllocPMEM -

   SYNOPSIS:
   void *memmgr_AllocPMEM (unsigned long size);

   PARAMETERS:
   .  size

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void *memmgr_AllocPMEM (unsigned long size)
{
  void   *buffer;

        #ifdef WITH_GENERAL_HEAP
  buffer = GetMem(myheap,size,0);
        #else
  buffer = malloc(size);
        #endif

  allocated += size;
  pmem      +=size;

        #ifdef WITH_HASH_CONTROL
  PushHash(buffer, size, 'P');
        #endif

  return(buffer);
}


/****************************************************************************/
/*
   memmgr_FreePMEM -

   SYNOPSIS:
   void memmgr_FreePMEM (void *buffer);

   PARAMETERS:
   .  buffer

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void memmgr_FreePMEM (void *buffer)
{
        #ifdef WITH_HASH_CONTROL
  {
    size_t hsize = PopHash(buffer,'P');
    allocated -= hsize;
    pmem -= hsize;
  }
        #endif

        #ifdef WITH_GENERAL_HEAP
  DisposeMem(myheap,buffer);
        #else
  free(buffer);
        #endif
}



/****************************************************************************/
/*
   memmgr_AllocAMEM -

   SYNOPSIS:
   void *memmgr_AllocAMEM (unsigned long size);

   PARAMETERS:
   .  size

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void *memmgr_AllocAMEM (unsigned long size)
{
  void   *buffer;

        #ifdef WITH_GENERAL_HEAP
  buffer = GetMem(myheap,size,0);
        #else
  buffer = malloc(size);
        #endif

  allocated += size;
  amem      += size;

        #ifdef WITH_HASH_CONTROL
  PushHash(buffer, size, 'A');
        #endif

  return(buffer);
}


/****************************************************************************/
/*
   memmgr_FreeAMEM -

   SYNOPSIS:
   void memmgr_FreeAMEM (void *buffer);

   PARAMETERS:
   .  buffer

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void memmgr_FreeAMEM (void *buffer)
{
        #ifdef WITH_HASH_CONTROL
  {
    size_t hsize = PopHash(buffer,'A');
    allocated -= hsize;
    amem -= hsize;
  }
        #endif

        #ifdef WITH_GENERAL_HEAP
  DisposeMem(myheap,buffer);
        #else
  free(buffer);
        #endif
}


/****************************************************************************/
/*
   memmgr_AllocTMEM -

   SYNOPSIS:
   void *memmgr_AllocTMEM (unsigned long size);

   PARAMETERS:
   .  size

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void *memmgr_AllocTMEM (unsigned long size, int kind)
{
  void   *buffer;


  if (kind==TMEM_XFER || kind==TMEM_CPL ||
      kind==TMEM_LOWCOMM || kind==TMEM_CONS || kind==TMEM_IDENT)
  {
    size_t real_size = size+sizeof(size_t);

                #ifndef DYNAMIC_MEMORY_ALLOCMODEL
    buffer = GetFreelistMemory(MGHEAP(dddctrl.currMG), real_size);
                #else
    buffer = GetMemoryForObject(dddctrl.currMG,real_size,MAOBJ);
                #endif
    if (buffer!=NULL)
    {
      /* store size at the beginning of memory chunk */
      *(size_t *)buffer = real_size;

      /* hide this information */
      buffer = (void *)(((char *)buffer) + sizeof(size_t));

      mem_from_ug_freelists += real_size;

      /*
         printf("%4d:    X MEMM adr=%08x kind=%d size=%ld\n", me,
                      buffer, kind, size);
       */
    }
  }
  else
  {
                #ifdef WITH_GENERAL_HEAP
    buffer = GetMem(myheap,size,0);
                #else
    buffer = malloc(size);
                #endif

    allocated += size;
    tmem      += size;

    /*
       printf("%4d:    O MEMM adr=%08x kind=%d size=%ld\n", me,
                    buffer, kind, size);
     */
  }

        #ifdef WITH_HASH_CONTROL
  PushHash(buffer, size, 'T');
        #endif

  return(buffer);
}


/****************************************************************************/
/*
   memmgr_FreeTMEM -

   SYNOPSIS:
   void memmgr_FreeTMEM (void *buffer);

   PARAMETERS:
   .  buffer

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void memmgr_FreeTMEM (void *buffer, int kind)
{
        #ifdef WITH_HASH_CONTROL
  {
    size_t hsize = PopHash(buffer,'T');
    allocated -= hsize;
    tmem -= hsize;
  }
        #endif


  if (kind==TMEM_XFER || kind==TMEM_CPL ||
      kind==TMEM_LOWCOMM || kind==TMEM_CONS || kind==TMEM_IDENT)
  {
    size_t real_size;

    /*
       printf("%4d:    X MEMF adr=%08x kind=%d\n", me, buffer, kind);
     */

    /* get real_size from beginning of buffer */
    buffer = (void *)(((char *)buffer) - sizeof(size_t));
    real_size = *(size_t *)buffer;

                #ifndef DYNAMIC_MEMORY_ALLOCMODEL
    PutFreelistMemory(MGHEAP(dddctrl.currMG), buffer, real_size);
                #else
    PutFreeObject(dddctrl.currMG,buffer,real_size,MAOBJ);
                #endif

    /*
       mem_from_ug_freelists -= real_size;
     */
  }
  else
  {
    /*
       printf("%4d:    O MEMF adr=%08x kind=%d\n", me, buffer, kind);
     */

                #ifdef WITH_GENERAL_HEAP
    DisposeMem(myheap,buffer);
                #else
    free(buffer);
                #endif
  }
}


/****************************************************************************/

void memmgr_MarkHMEM (long *theMarkKey)
{
  INT myMarkKey;
  MarkTmpMem(MGHEAP(dddctrl.currMG), &myMarkKey);
  *theMarkKey = (long)myMarkKey;
}

void *memmgr_AllocHMEM (size_t size, long theMarkKey)
{
  void *buffer;
  buffer = GetTmpMem(MGHEAP(dddctrl.currMG), size, (INT)theMarkKey);

  /*
     printf("%4d:    H MEMM adr=%08x           size=%ld\n", me, buffer, size);
   */

  return(buffer);
}

void memmgr_ReleaseHMEM (long theMarkKey)
{
  ReleaseTmpMem(MGHEAP(dddctrl.currMG), (INT)theMarkKey);
}


/****************************************************************************/

#define MAX_MALLOCS 256


/*
        DetectAllocatableMemory()

        This function tries to find out the amount of memory
        which can be allocated from the local heap. This is
        done by allocating a sequence of memory blocks,
        starting with big ones. NOTE: this function shouldn't
        be used on systems where virtual memory is available
        (e.g., workstations); the page swapping mechanism
        will break down if you try to do DetectAllocatableMemory()
        on such systems.
 */

#define MAX_MALLOCS 256

static size_t DetectAllocatableMemory (void)
{
  void *buffers[MAX_MALLOCS];
  size_t s = 1024*1024*1024;
  size_t all = 0;
  int i = 0;

  do {
    buffers[i] =(void *) malloc(s);
    if (buffers[i]==NULL)
    {
      /* couldnt get memory of size s, try with half size */
      s = s/2;
    }
    else
    {
      /* could allocate mem, continue */
      all += s;
      i++;
    }
  } while (i<MAX_MALLOCS && s>32);

  /* free memory */
  while (i>0)
  {
    i--;
    free(buffers[i]);
  };

  return(all);
}


/****************************************************************************/
/*
   memmgr_Init -

   SYNOPSIS:
   void memmgr_Init (void);

   PARAMETERS:
   .  void

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void memmgr_Init (void)
{
        #ifdef WITH_GENERAL_HEAP
  {
    void *buffer;

    buffer = malloc(HEAP_SIZE);
    if (buffer==NULL) {
      printf("not enough memory for DDD heap\n");
      return;
    }

    myheap = NewHeap(GENERAL_HEAP,HEAP_SIZE,buffer);
  }

        #else

  /* detect size of allocatable memory */
                #ifdef DETECT_MEMORY_SIZE
  printf("%4d: MemMgr. detecting size of allocatable memory ...\n", me);
  printf("%4d: MemMgr. size of allocatable memory: %ld\n", me,
         (unsigned long)DetectAllocatableMemory());
                #endif
        #endif


        #ifdef WITH_HASH_CONTROL
  {
    int i;

    /* init hash table */
    for(i=0; i<HASHTAB_SIZE; i++)
      htab[i] = NULL;
    nHashEntries = 0;
  }
        #endif
}


/****************************************************************************/

#endif /* ModelP */
