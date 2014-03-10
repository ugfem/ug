// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      boxtree.c                                                     */
/*                                                                          */
/* Purpose:   boxtree data structure                                        */
/*                                                                          */
/* Author:    Michael Lampe                                                 */
/*            IWR - Technische Simulation                                   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            Email: Michael.Lampe@iwr.uni-heidelberg.de                    */
/*                                                                          */
/* History:   20010927 begin, ug3.8                                         */
/*                                                                          */
/* Remarks:   see Visual Computer (1987) 3:236-249                          */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include "config.h"

#include "ugtypes.h"
#include "architecture.h"
#include "general.h"
#include "boxtree.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define NCUT 10

#ifndef MIN
#define MIN(a, b)   ((a) <= (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)   ((a) >= (b) ? (a) : (b))
#endif

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/*    Sort - build a balanced boxtree in time O(n*log(n))                   */
/*																			*/
/****************************************************************************/

static void Sort(BT_OBJECT **box, int dim, int left, int right,
                 BT_OBJECT **root, double range[DIM][2])
{
  int i, j, k, l, r, m;
  double key, left_range[DIM][2], right_range[DIM][2];
  BT_OBJECT *p, *t;

  /* find median, 1st pass quicksort */
  l = left;
  r = right;
  m = (l+r)/2;
  while (r-l >= NCUT) {
    key = box[m]->range[dim][0];
    i = l;
    j = r;
    do {
      while (box[i]->range[dim][0] < key) i++;
      while (key < box[j]->range[dim][0]) j--;
      if (i <= j) {
        t = box[i];
        box[i] = box[j];
        box[j] = t;
        i++;
        j--;
      }
    } while (i <= j);
    if (j < m) l = i;
    if (m < i) r = j;
  }
  /* 2nd pass straight selection */
  for (i = l; i <= m; i++) {
    k = i;
    t = box[i];
    for (j = i+1; j <= r; j++)
      if (box[j]->range[dim][0] < t->range[dim][0]) {
        k = j;
        t = box[j];
      }
    box[k] = box[i];
    box[i] = t;
  }

  /* median becomes root */
  *root = p = box[m];

  /* determine sons */
  if (left < m) {
    if (left < m-1) {
      /* entire left subtree */
      Sort(box, (dim+1)%DIM, left, m-1, &(p->left_son), left_range);
      p->left_range[0] = left_range[dim][0];
      p->left_range[1] = left_range[dim][1];
    }
    else {
      /* left son is leave */
      p->left_son = t = box[left];
      p->left_range[0] = t->range[dim][0];
      p->left_range[1] = t->range[dim][1];
      t->left_range[1] = -MAX_D;
      t->right_range[1] = -MAX_D;
      for (i = 0; i < DIM; i++) {
        left_range[i][0] = t->range[i][0];
        left_range[i][1] = t->range[i][1];
      }
    }
    if (m+1 < right) {
      /* entire right subtree */
      Sort(box, (dim+1)%DIM, m+1, right, &(p->right_son), right_range);
      p->right_range[0] = right_range[dim][0];
      p->right_range[1] = right_range[dim][1];
    }
    else {
      /* right son is leave */
      p->right_son = t = box[right];
      p->right_range[0] = t->range[dim][0];
      p->right_range[1] = t->range[dim][1];
      t->left_range[1] = -MAX_D;
      t->right_range[1] = -MAX_D;
      for (i = 0; i < DIM; i++) {
        right_range[i][0] = t->range[i][0];
        right_range[i][1] = t->range[i][1];
      }
    }
  }
  else {
    /* no left son, right son is leave, or no sons(!) */
    p->right_son = t = box[right];
    p->left_range[1] = -MAX_D;
    p->right_range[0] = t->range[dim][0];
    p->right_range[1] = t->range[dim][1];
    t->left_range[1] = -MAX_D;
    t->right_range[1] = -MAX_D;
    for (i = 0; i < DIM; i++) {
      right_range[i][0] = t->range[i][0];
      right_range[i][1] = t->range[i][1];
      left_range[i][0] = MAX_D;
      left_range[i][1] = -MAX_D;
    }
  }

  /* update ranges */
  for (i = 0; i < DIM; i++) {
    range[i][0] = MIN(p->range[i][0], MIN(left_range[i][0], right_range[i][0]));
    range[i][1] = MAX(p->range[i][1], MAX(left_range[i][1], right_range[i][1]));
  }
}

/****************************************************************************/
/*D
   BT_Init - Initialize boxtree data structure

   SYNOPSIS:
   void BT_Init(BT_OBJECT **box, int n, BOXTREE *tree)

   PARAMETERS:
   .  box  - array of boxes to build tree from
   .  n    - #boxes
   .  tree - output

   DESCRIPTION:
   .  none

   RETURN VALUE:
   .  none

   D*/
/****************************************************************************/

void BT_Init(BT_OBJECT **box, int n, BOXTREE *tree)
{
  Sort(box, 0, 0, n-1, &(tree->root), tree->range);
}

/****************************************************************************/
/*																			*/
/*   Search - find all boxes that overlap with query box in time O(log n)   */
/*																			*/
/****************************************************************************/

static void Search(BT_OBJECT *root, int dim, double range[DIM][2], BT_FUNC func, void *data)
{
  int i;

  if (root->range[dim][0] <= range[dim][1]) {
    for (i = 0; i < DIM; i++)
      if (root->range[i][0] > range[i][1] || range[i][0] > root->range[i][1]) break;
    if (i == DIM) func(root, data);
    if (root->left_range[1] >= range[dim][0])
      Search(root->left_son, (dim+1)%DIM, range, func, data);
    if (root->right_range[1] >= range[dim][0] && root->right_range[0] <= range[dim][1])
      Search(root->right_son, (dim+1)%DIM, range, func, data);
  }
  else if (root->left_range[1] >= range[dim][0] && root->left_range[0] <= range[dim][1])
    Search(root->left_son, (dim+1)%DIM, range, func, data);
}

/****************************************************************************/
/*																			*/
/*    SearchP - find all boxes that contain query point in time O(log n)    */
/*																			*/
/****************************************************************************/

static void SearchP(BT_OBJECT *root, int dim, double *x, BT_FUNC func, void *data)
{
  int i;

  if (root->range[dim][0] <= x[dim]) {
    for (i = 0; i < DIM; i++)
      if (root->range[i][0] > x[i] || x[i] > root->range[i][1]) break;
    if (i == DIM) func(root, data);
    if (root->left_range[1] >= x[dim])
      SearchP(root->left_son, (dim+1)%DIM, x, func, data);
    if (root->right_range[1] >= x[dim] && root->right_range[0] <= x[dim])
      SearchP(root->right_son, (dim+1)%DIM, x, func, data);
  }
  else if (root->left_range[1] >= x[dim] && root->left_range[0] <= x[dim])
    SearchP(root->left_son, (dim+1)%DIM, x, func, data);
}

/****************************************************************************/
/*D
   BT_Search - call function for all objects in target boxes

   SYNOPSIS:
   void BT_Search(BOXTREE *tree, double range[DIM][2], BT_FUNC func, void *data);

   PARAMETERS:
   .  tree  - boxtree
   .  range - query box
   .  func  - callback function
   .  data  - user data

   DESCRIPTION:
   .  none

   RETURN VALUE:
   .  none

   D*/
/****************************************************************************/

void BT_Search(BOXTREE *tree, double range[DIM][2], BT_FUNC func, void *data)
{
  Search(tree->root, 0, range, func, data);
}

/****************************************************************************/
/*D
   BT_SearchP - call function for all objects in target boxes

   SYNOPSIS:
   void BT_SearchP(BOXTREE *tree, double *x, BT_FUNC func, void *data);

   PARAMETERS:
   .  tree  - boxtree
   .  x     - query point
   .  func  - callback function
   .  data  - user data

   DESCRIPTION:
   .  none

   RETURN VALUE:
   .  none

   D*/
/****************************************************************************/

void BT_SearchP(BOXTREE *tree, double *x, BT_FUNC func, void *data)
{
  SearchP(tree->root, 0, x, func, data);
}
