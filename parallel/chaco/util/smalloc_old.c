// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <malloc.h>

static int nmalloc=0;   /* number of calls to malloc */
static int nfree=0;     /* number of calls to free */

static struct smalloc_debug_data {
  int order;                            /* which smalloc call is it? */
  unsigned size;                        /* size of malloc invocation */
  double *ptr;                          /* memory location returned */
  struct smalloc_debug_data *next;      /* pointer to next element */
} *top = NULL;

/* Safe version of malloc.  Does not initialize memory .*/
double *c_malloc(n)
unsigned int n;                 /* number of bytes to be allocated */
{
  extern int DEBUG_MEMORY;      /* use debug memory allocator? */
  double *pntr;                 /* return value */
  struct smalloc_debug_data *new;       /* data structure for malloc data */
  void exit();

  nmalloc++;
  if (n==0) {
    {char buf[150]; sprintf(buf,"smalloc: Non-positive argument. (%u)\n", n);UserWrite(buf);}
    exit(0);
  }

  pntr = (double *) malloc(n);

  if (pntr == NULL) {
    {char buf[150]; sprintf(buf,"Program out of space while attempting to allocate (%u).  Sorry!\n",n);UserWrite(buf);}
    exit(0);
  }

  if (DEBUG_MEMORY > 1) {
    new = (struct smalloc_debug_data *)
          malloc(sizeof(struct smalloc_debug_data));

    if (new == NULL) {
      {char buf[150]; sprintf(buf,"Program out of space while attempting to allocate (%u).  Sorry!\n",n);UserWrite(buf);}
      exit(0);
    }

    new->order = nmalloc;
    new->size = n;
    new->ptr = pntr;
    new->next = top;
    top = new;
  }

  return(pntr);
}


/* Safe version of free. */
int c_free(ptr)
char *ptr;
{
  extern int DEBUG_MEMORY;              /* use debug memory allocator? */
  struct smalloc_debug_data *dbptr;     /* loops through debug list */
  struct smalloc_debug_data **prev;     /* holds previous pointer */

  if (DEBUG_MEMORY > 1) {
    if (ptr != NULL) {          /* search through debug list for it */
      prev = &top;
      for (dbptr=top; dbptr!=NULL && (char *) dbptr->ptr!=ptr; dbptr=dbptr->next) {
        prev = &(dbptr->next);
      }
      if (dbptr == NULL) {
        {char buf[150]; sprintf(buf,"Memory error: In sfree, address not found in debug list (%d)\n", ptr);UserWrite(buf);}
      }
      else {
        *prev = dbptr->next;
        free((char *) dbptr);
      }
    }
  }


  if (ptr != NULL) {
    nfree++;
    free(ptr);
    ptr = NULL;
  }

  return(0);
}


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
        {char buf[150]; sprintf(buf," order=%d, size=%u, location=%d\n", dbptr->order, dbptr->size, dbptr->ptr);UserWrite(buf);}
      }
    }
  }
}
