// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>

static int nmalloc=0;   /* number of calls to malloc */
static int nfree=0;     /* number of calls to free */

static struct smalloc_debug_data {
  int order;                            /* which smalloc call is it? */
  unsigned size;                        /* size of malloc invocation */
  double *ptr;                          /* memory location returned */
  struct smalloc_debug_data *next;      /* pointer to next element */
} *top = NULL;


void smalloc_stats()
{
  extern int DEBUG_MEMORY;              /* use debug memory allocator? */
  struct smalloc_debug_data *dbptr;     /* loops through debug list */

  if (DEBUG_MEMORY > 0) {
    {char buf[150]; sprintf(buf,"Calls to smalloc = %d,  Calls to sfree = %d\n", nmalloc, nfree);UserWrite(buf);}
  }
  if (DEBUG_MEMORY > 1) {
    if (top != NULL) {
      {char buf[150]; sprintf(buf,"Remaining allocations:\n");UserWrite(buf);}
      for (dbptr=top; dbptr!=NULL; dbptr=dbptr->next) {
        {char buf[150]; sprintf(buf," order=%d, size=%u, location=%d\n", dbptr->order,
                                dbptr->size, dbptr->ptr);UserWrite(buf);}
      }
    }
  }
}
