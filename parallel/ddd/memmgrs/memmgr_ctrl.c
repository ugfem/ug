// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      memmgr_ctrl.c                                                 */
/*                                                                          */
/* Purpose:   basic memory management module                                */
/*            (with standard malloc() calls)                                */
/*            with control feature, counts allocated bytes!                 */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/04/27 kb  begin                                            */
/*            96/01/20 kb  updated to DDD V1.5                              */
/*            96/09/05 kb  implemented for dddic                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

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
#include <assert.h>

#include "ppif.h"


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


#define HASHTAB_SIZE  15731

#define HASH_FUNC(k)   ((k)%HASHTAB_SIZE)


/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/


typedef struct _HASH_ENTRY
{
  void   *ptr;           /* hashed key:  pointer to memory block */
  size_t size;           /* hashed data: size of memory block    */

  int flags;

  struct _HASH_ENTRY *next;

} HASH_ENTRY;


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

unsigned long memmgr_AllocatedPMEM;
unsigned long memmgr_AllocatedOMEM;
unsigned long memmgr_AllocatedAMEM;
unsigned long memmgr_AllocatedTMEM;


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

HASH_ENTRY *htab[HASHTAB_SIZE];
int nHashEntries;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


HASH_ENTRY *NewHashEntry (void *ptr, size_t size)
{
  HASH_ENTRY *he;
  /*
     printf("%4d: alloc %08x %d\n", me,ptr,(unsigned long)size);
   */

  he = (HASH_ENTRY *) malloc(sizeof(HASH_ENTRY));
  he->ptr = ptr;
  he->size = size;
  he->flags = 0;
  he->next = NULL;

  nHashEntries++;

  return(he);
}

void FreeHashEntry (HASH_ENTRY *he)
{
  /*
     printf("%4d: free  %08x %d\n", me,he->ptr,(unsigned long)he->size);
   */
  nHashEntries--;
  free(he);
}


static void PushHash (void *ptr, size_t size)
{
  unsigned int idx = HASH_FUNC(((unsigned long)ptr));

  if (htab[idx] == NULL)
  {
    /* no collision */
    htab[idx] = NewHashEntry(ptr, size);
  }
  else
  {
    /* collision, find entry or none */
    HASH_ENTRY *he;

    for(he=htab[idx]; he->next!=NULL && he->ptr!=ptr; he=he->next)
      ;

    if (he->ptr==ptr)
    {
      printf("%4d: MEMMGR-ERROR, double alloc at %08x, size %ld\n",
             me, ptr, (unsigned long)size);
      exit(1);
    }

    he->next = NewHashEntry(ptr, size);
  }
}


static size_t PopHash (void *ptr)
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

      FreeHashEntry(he);
      return(s);
    }
  }

  printf("%4d: MEMMGR-ERROR, no alloc for free at %08x\n", me, ptr);
  assert(0);

  return(0);       /* never reached */
}



void HashMarkAll (void)
{
  int i;

  for(i=0; i<HASHTAB_SIZE; i++)
  {
    HASH_ENTRY *he;
    for(he=htab[i]; he!=NULL; he=he->next)
      he->flags = 1;
  }
}


void HashShowMarks (void)
{
  int i;

  for(i=0; i<HASHTAB_SIZE; i++)
  {
    HASH_ENTRY *he;
    for(he=htab[i]; he!=NULL; he=he->next)
    {
      if (he->flags==0)
      {
        printf("%4d: MALLOC adr=%08x size=%ld\n",
               me, he->ptr, (unsigned long) he->size);
      }
    }
  }
}


/****************************************************************************/


void *memmgr_AllocPMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreePMEM (void *buffer)
{
  free(buffer);
}




void *memmgr_AllocOMEM (size_t size, int ddd_typ, int proc, int attr)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeOMEM (void *buffer, size_t size, int ddd_typ)
{
  free(buffer);
}




void *memmgr_AllocAMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeAMEM (void *buffer)
{
  free(buffer);
}


void *memmgr_AllocTMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);

  PushHash(buffer, size);
  memmgr_AllocatedTMEM += size;

  return(buffer);
}


void memmgr_FreeTMEM (void *buffer)
{
  memmgr_AllocatedTMEM -= PopHash(buffer);

  free(buffer);
}


void *memmgr_AllocHMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeHMEM (void *buffer)
{
  free(buffer);
}


void memmgr_MarkHMEM (void)
{}


void memmgr_ReleaseHMEM (void)
{}



void memmgr_Init (void)
{
  int i;

  /* init hash table */
  for(i=0; i<HASHTAB_SIZE; i++)
    htab[i] = NULL;
  nHashEntries = 0;

  memmgr_AllocatedPMEM = 0;
  memmgr_AllocatedOMEM = 0;
  memmgr_AllocatedAMEM = 0;
  memmgr_AllocatedTMEM = 0;
}
