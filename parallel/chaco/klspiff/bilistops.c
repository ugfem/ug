// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        "../main/structs.h"
#include        "../main/defs.h"


/* Note: bi-directional lists aren't assumed to be sorted. */

void add2bilist(lptr, list)             /* add val to unsorted list */
struct bilist *lptr;                    /* element to add */
struct bilist **list;                   /* list added to */
{
  lptr->next = *list;
  if (*list != NULL) (*list)->prev = lptr;
  lptr->prev = NULL;
  *list = lptr;
}


void removebilist(lptr, list)
struct bilist *lptr;            /* ptr to element to remove */
struct bilist **list;           /* head of list to remove it from */
/* Remove an element from a bidirectional list. */
{
  if (lptr->next != NULL) lptr->next->prev = lptr->prev;
  if (lptr->prev != NULL) lptr->prev->next = lptr->next;
  else *list = lptr->next;
}


void movebilist(lptr, oldlist, newlist)
struct bilist *lptr;            /* ptr to element to move */
struct bilist **oldlist;        /* head of list to remove it from */
struct bilist **newlist;        /* head of list to add it to */
/* Move an element from a old bidirectional list to new one. */
{
  void removebilist(), add2bilist();

  removebilist(lptr, oldlist);

  add2bilist(lptr, newlist);
}
