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


/****************************************************************************/

/* define this to protocol all alloc/free requests via hashtable */
/*
   #define WITH_HASH_CONTROL
 */

/* define this to map all PMEM, AMEM and TMEM requests to a UG general heap */
/*
   #define WITH_GENERAL_HEAP
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
#define HASHTAB_SIZE  15731
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

static void FreeHashEntry (HASH_ENTRY *he)
{
  /*
     printf("%4d: free  %c %08x %ld\n", me,he->info,he->ptr,(unsigned long)he->size);
   */

  nHashEntries--;
  free(he);
}


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

  fflush(stdout);
}

/****************************************************************************/



void *memmgr_AllocOMEM (size_t size, int ddd_type, int prio, int attr)
{
  void   *buffer;

  buffer = GetFreelistMemory(MGHEAP(dddctrl.currMG), size);

  /*
     printf("%4d: memmgr_AllocOMem: size=%05d ddd_type=%02d prio=%d attr=%d\n",
     me,size,ddd_type,prio,attr);
   */

  return(buffer);
}


void memmgr_FreeOMEM (void *buffer, size_t size, int ddd_type)
{
  /*
     printf("%d: memmgr_FreeOMEM(): buffer=%x, ddd_type=%d\n", me, buffer, ddd_type);
   */

  PutFreelistMemory(MGHEAP(dddctrl.currMG), buffer, size);
}




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


void *memmgr_AllocTMEM (unsigned long size)
{
  void   *buffer;

        #ifdef WITH_GENERAL_HEAP
  buffer = GetMem(myheap,size,0);
        #else
  buffer = malloc(size);
        #endif

  allocated += size;
  tmem      += size;

        #ifdef WITH_HASH_CONTROL
  PushHash(buffer, size, 'T');
        #endif

  return(buffer);
}


void memmgr_FreeTMEM (void *buffer)
{
        #ifdef WITH_HASH_CONTROL
  {
    size_t hsize = PopHash(buffer,'T');
    allocated -= hsize;
    tmem -= hsize;
  }
        #endif

        #ifdef WITH_GENERAL_HEAP
  DisposeMem(myheap,buffer);
        #else
  free(buffer);
        #endif
}


/****************************************************************************/

void *memmgr_AllocHMEM (unsigned long size)
{
  return(NULL);
}

void memmgr_FreeHMEM (void *buffer)
{}

void memmgr_MarkHMEM (void)
{}

void memmgr_ReleaseHMEM (void)
{}


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
