// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*****************************************************************************
* File:      bbtree.c                                                       *
* Purpose:   bounding box tree data structure                               *
*                                                                                                                                        *
* Author:	  O. Sterz                                                       *
*                                                                           *
* History:   Nov 2002                                                       *
* Remarks:   similar to boxtree.c, but more natural to compute distances    *
*****************************************************************************/

/*****************************************************************************
* include files                                                             *
*            system include files                                           *
*            application include files                                      *
*****************************************************************************/
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include "general.h"
#include "compiler.h"
#include "bbtree.h"

USING_UG_NAMESPACE

/*****************************************************************************
* defines in the following order                                            *
*        compile time constants defining static data size (i.e. arrays)     *
*        other constants                                                    *
*        macros                                                             *
*****************************************************************************/
#define MIN(a, b)   ((a) <= (b) ? (a) : (b))
#define MAX(a, b)   ((a) >= (b) ? (a) : (b))
#define SQR(x)      ((x)*(x))

/*****************************************************************************
* PRIVATE data structures (global to this source file)                      *
*****************************************************************************/
typedef struct tree_point_dist_data
{
  BBT_POINT_DIST_FUNC *dist;
  DOUBLE *x;
  DOUBLE min;
  void *obj;
} TREE_POINT_DIST_DATA;


/*****************************************************************************
* PRIVATE global variables (global to this source file)                     *
*****************************************************************************/
/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

INT theBBTDim = 0;
HEAP *theBBTHeap = NULL;
DOUBLE tmp[MAX_DIM];

/*****************************************************************************
* PRIVATE function declarations/definitions (used in this source file)      *
*****************************************************************************/
/*---------------------------------------------------------------------------*
* BBoxPointDistance2 - computes min and max squared distance of point to bb *
*---------------------------------------------------------------------------*/
static void BBoxPointDistance2(BBT_BBOX *bb, DOUBLE *x, DOUBLE *min,
                               DOUBLE *max)
#ifdef TIGHT_BOXES
{
  DOUBLE dmax, dmin, dll_i, dur_i, d, m, M;
  INT i, j;

  assert(theBBTDim <= MAX_DIM);
  for (i=0; i<theBBTDim; i++)
    tmp[i] = 0.0;

  dmin = 0.0;
  for (i=0; i<theBBTDim; i++)
  {
    d = x[i] - bb->ll[i];
    dll_i = SQR(d);
    d = x[i] - bb->ur[i];
    dur_i = SQR(d);
    dmin += x[i] < bb->ll[i] ? dll_i : (x[i] > bb->ur[i] ? dur_i : 0.0);

    m = MIN(dll_i, dur_i);
    M = MAX(dll_i, dur_i);
    for (j=0; j<theBBTDim; j++)
      tmp[j] += j == i ? m : M;
  }

  dmax  = tmp[0];
  for (i=1; i<theBBTDim; i++)
    dmax = MIN(dmax, tmp[i]);

  *min = dmin;
  *max = dmax;
}
#else
{
  DOUBLE dmax, dmin, dll_i, dur_i, d;
  INT i;

  dmax = dmin = 0.0;
  for (i=0; i<theBBTDim; i++)
  {
    d = x[i] - bb->ll[i];
    dll_i = SQR(d);
    d = x[i] - bb->ur[i];
    dur_i = SQR(d);
    dmin += x[i] < bb->ll[i] ? dll_i : (x[i] > bb->ur[i] ? dur_i : 0.0);
    dmax += MAX(dll_i, dur_i);
  }
  *min = dmin;
  *max = dmax;
}
#endif

/*---------------------------------------------------------------------------*
* MinMaxBBoxPointDist2 -                                                    *
*---------------------------------------------------------------------------*/
static DOUBLE MinMaxBBoxPointDist2(BBT_NODE *node, DOUBLE *x, DOUBLE minmax)
{
  DOUBLE lmin, lmax, rmin, rmax;

  if (node->left == NULL)
  {
    assert(node->right == NULL);
    return minmax;
  }

  BBoxPointDistance2(node->left->bb, x, &lmin, &lmax);
  if (lmax < minmax) minmax = lmax;
  BBoxPointDistance2(node->right->bb, x, &rmin, &rmax);
  if (rmax < minmax) minmax = rmax;
  if (lmin < rmin)
  {
    if (lmin < minmax)
    {
      minmax = MinMaxBBoxPointDist2(node->left, x, minmax);
      if (rmin < minmax)
        minmax = MinMaxBBoxPointDist2(node->right, x, minmax);
    }
  }
  else
  {
    if (rmin < minmax)
    {
      minmax = MinMaxBBoxPointDist2(node->right, x, minmax);
      if (lmin < minmax)
        minmax = MinMaxBBoxPointDist2(node->left, x, minmax);
    }
  }

  return minmax;
}

/*---------------------------------------------------------------------------*
* ClosestBoxesToPoint - find a minimum set of boxes, where one of them will *
*	                     include the nearest object to a point               *
*---------------------------------------------------------------------------*/
static void ClosestBBoxesToPoint(BBT_NODE *node, DOUBLE *x, BBT_FUNC func,
                                 void *bypass, DOUBLE *minmax)
{
  DOUBLE min, max;

  if (node->left == NULL)       /* leave */
  {
    assert(node->right == NULL);
    func(node->bb->object, bypass);
    return;
  }

  BBoxPointDistance2(node->left->bb, x, &min, &max);
  if (min < *minmax)
    ClosestBBoxesToPoint(node->left, x, func, bypass, minmax);
  BBoxPointDistance2(node->right->bb, x, &min, &max);
  if (min < *minmax)
    ClosestBBoxesToPoint(node->right, x, func, bypass, minmax);
}

/*---------------------------------------------------------------------------*
* BBoxForBBoxes - creates a new bbox containing all n given bboxes          *
*---------------------------------------------------------------------------*/
static BBT_BBOX *BBoxForBBoxes(BBT_BBOX **bboxes, INT n)
{
  BBT_BBOX *bbox;
  INT i, j;

  if (n < 1) return NULL;
  bbox = BBT_NewBBox(theBBTHeap, theBBTDim, bboxes[0]->ll, bboxes[0]->ur,
                     NULL);
  for (i=1; i<n; i++)
    for (j=0; j<theBBTDim; j++)
    {
      if (bboxes[i]->ll[j] < bbox->ll[j]) bbox->ll[j] = bboxes[i]->ll[j];
      if (bboxes[i]->ur[j] > bbox->ur[j]) bbox->ur[j] = bboxes[i]->ur[j];
    }
  return bbox;;
}

/*---------------------------------------------------------------------------*
* NewNode - creates a new bouding box tree node                             *
*---------------------------------------------------------------------------*/
static BBT_NODE *NewNode(BBT_BBOX *bbox)
{
  BBT_NODE *newNode;

  newNode = (BBT_NODE *) GetFreelistMemory(theBBTHeap, sizeof(BBT_NODE));
  if (newNode == NULL) return NULL;
  newNode->left = newNode->right = NULL;
  newNode->bb = bbox;
  return newNode;
}

/*---------------------------------------------------------------------------*
* BuildTree - builds a new bounding box tree                                *
*---------------------------------------------------------------------------*/
static BBT_NODE *BuildTree(BBT_BBOX **bboxes, INT n)
{
  BBT_NODE *node;
  BBT_BBOX *bbox, **lbboxes, **rbboxes, **new_bboxes;
  DOUBLE ext, maxext, cut;
  INT i, dir, nl, nr, MarkKey;

  if (n < 1) return NULL;
  assert(bboxes != NULL);
  assert(bboxes[0] != NULL);

  /* leave node ? */
  if (n == 1)
    return NewNode(bboxes[0]);

  bbox = BBoxForBBoxes(bboxes, n);
  node = NewNode(bbox);

  /* split set of bounding boxes by location of the center points */
  maxext = 0.0;
  for (i=0; i<theBBTDim; i++)
  {
    ext = bbox->ur[i] - bbox->ll[i];
    if (ext > maxext)
    {
      maxext = ext;
      dir = i;
    }
  }
  assert(maxext > 0.0);
  cut = (bbox->ur[dir] + bbox->ll[dir])/2.0;

  nr = 0;
  for (i=0; i<n; i++)
    if ((bboxes[i]->ll[dir] + bboxes[i]->ur[dir])/2.0 > cut) nr++;
  nl = n - nr;
  if (MarkTmpMem(theBBTHeap, &MarkKey))
  {
    printf("ERROR in low/bbtree.c/BuildTree(): Out of memory. Enlarge UG stack size.\n");
    fprintf( stderr, "ERROR in low/bbtree.c/BuildTree(): Out of memory. Enlarge UG stack size.\n");
    assert(0);
  }
  new_bboxes = (BBT_BBOX **) GetTmpMem(theBBTHeap, n*sizeof(BBT_BBOX **),
                                       MarkKey);
  lbboxes = new_bboxes;
  rbboxes = new_bboxes + nl;

  nr = nl = 0;
  for (i=0; i<n; i++)
    if ((bboxes[i]->ll[dir] + bboxes[i]->ur[dir])/2.0 > cut)
      rbboxes[nr++] = bboxes[i];
    else
      lbboxes[nl++] = bboxes[i];
  if (nr == 0)
  {
    nl /= 2;
    nr = n - nl;
    rbboxes = &lbboxes[nl];
  }
  else if (nl == 0)
  {
    nr /= 2;
    nl = n - nr;
    lbboxes = &rbboxes[nr];
  }

  /* extend the tree */
  node->left = BuildTree(lbboxes, nl);
  node->right = BuildTree(rbboxes, nr);
  ReleaseTmpMem(theBBTHeap, MarkKey);

  return node;
}

static void TreePointDistanceCallback(void *obj, void *bypass)
{
  TREE_POINT_DIST_DATA *data;
  DOUBLE d;

  data = (TREE_POINT_DIST_DATA *) bypass;
  d = data->dist(data->x, obj);
  if (d < data->min)
  {
    data->min = d;
    data->obj = obj;
  }
}


/*****************************************************************************
* PUBLIC global variables definitions (exported by the include file)        *
*****************************************************************************/


/*****************************************************************************
* PUBLIC function definitions (declarations are in the include file)        *
*****************************************************************************/
/*---------------------------------------------------------------------------*/
/*D
   BBT_NewTree - builds a new bounding box tree from an array of pointer to
                 bounding boxes

   SYNOPSIS:
   BBT_TREE *BBT_NewTree(HEAP *theHeap, BBT_BBOX **bboxes, INT n, INT dim)

   PARAMETERS:
   .  theHeap - pointer to UG heap
   .  bboxes  - array of pointers to bounding boxes to build tree from
   .  n       - #bounding boxes
   .  dim     - dimension

   DESCRIPTION:
   .  none

   RETURN VALUE:
   .  pointer to new bounding box tree
   D*/
/*---------------------------------------------------------------------------*/
BBT_TREE* NS_PREFIX BBT_NewTree(HEAP *theHeap, BBT_BBOX **bboxes, INT n, INT dim)
{
  BBT_TREE *newTree;

  newTree = (BBT_TREE *) GetFreelistMemory(theHeap, sizeof(BBT_TREE));
  if (newTree == NULL) return NULL;
  newTree->dim = theBBTDim = dim;
  newTree->heap = theBBTHeap = theHeap;
  newTree->n = 0;
  newTree->root = BuildTree(bboxes, n);

  return newTree;
}

/*---------------------------------------------------------------------------*/
/*D
   BBT_NewBBox - creates a new bounding box of dimension dim

   SYNOPSIS:
   BBT_BBOX *BBT_NewBBox(HEAP *theHeap, INT dim, DOUBLE *ll, DOUBLE *ur,
                         void *obj)

   PARAMETERS:
   .  theHeap - pointer to UG heap
   .  dim  - dimension
   .  ll   - lower left corner
   .  ur   - upper right corner
   .  obj  - object contained

   DESCRIPTION:
   .  none

   RETURN VALUE:
   .  pointer to new empty bounding box
   D*/
/*---------------------------------------------------------------------------*/
BBT_BBOX* NS_PREFIX BBT_NewBBox(HEAP *theHeap, INT dim, DOUBLE *ll, DOUBLE *ur, void *obj)
{
  BBT_BBOX *newBB;
  INT i;

  newBB = (BBT_BBOX *) GetFreelistMemory(theHeap, sizeof(BBT_BBOX)
                                         +sizeof(DOUBLE)*2*dim);
  if (newBB == NULL) return NULL;
  newBB->object = obj;
  newBB->ll = (DOUBLE *) (newBB + 1);
  newBB->ur = newBB->ll + dim;
  for (i=0; i<dim; i++)
  {
    newBB->ll[i] = ll[i];
    newBB->ur[i] = ur[i];
  }
  return newBB;
}

/*---------------------------------------------------------------------------*/
/*D
   BBT_ClosestBBoxesToPoint - call function func for a small number boxes,
                              one of the boxes is assured to contain the
                                                          nearest object to the given point x

   SYNOPSIS:
   void BBT_ClosestBBoxesToPoint(BBT_TREE *tree, DOUBLE *x, BBT_FUNC func,
                                 void *bypass)

   PARAMETERS:
   .  tree   - bounding box tree
   .  x      - point
   .  func   - callback function
   .  bypass - pointer to some data that will be bypassed to func
   DESCRIPTION:
   .  none

   RETURN VALUE:
   .  none
   D*/
/*---------------------------------------------------------------------------*/
void NS_PREFIX BBT_ClosestBBoxesToPoint(BBT_TREE *tree, DOUBLE *x, BBT_FUNC func,
                                        void *bypass)
{
  DOUBLE minmax;

  if (tree == NULL) return;
  assert(x != NULL);
  assert(func != NULL);
  theBBTDim = tree->dim;
  theBBTHeap = tree->heap;
  minmax = MinMaxBBoxPointDist2(tree->root, x, DBL_MAX);
  ClosestBBoxesToPoint(tree->root, x, func, bypass, &minmax);
}

/*---------------------------------------------------------------------------*/
/*D
   BBT_TreePointDistance - compute the distance of a point x to the closest
                           object in the tree; a call back function dist() is
                                                   used to compute the distance

   SYNOPSIS:
   DOUBLE BBT_TreePointDistance(BBT_TREE *tree, DOUBLE *x, void **obj,
                                                         BBT_POINT_DIST_FUNC dist)

   PARAMETERS:
   .  tree   - bounding box tree
   .  x      - point
   .  obj    - pointer to pointer to closest object in tree found
   .  dist   - callback function that is able to compute a distance between a
            point and an object that is contained tn the tree bounding boxes
                    DOUBLE dist(DOUBLE *x, *obj)

   DESCRIPTION:
   .  none

   RETURN VALUE:
   .  distance to closest object
   D*/
/*---------------------------------------------------------------------------*/
DOUBLE NS_PREFIX BBT_TreePointDistance(BBT_TREE *tree, DOUBLE *x, void **obj,
                                       BBT_POINT_DIST_FUNC dist)
{
  TREE_POINT_DIST_DATA bypass;
  DOUBLE minmax;

  if (tree == NULL) return DBL_MAX;;
  assert(x != NULL);
  theBBTDim = tree->dim;
  theBBTHeap = tree->heap;
  minmax = MinMaxBBoxPointDist2(tree->root, x, DBL_MAX);
  bypass.min = DBL_MAX;
  bypass.dist = dist;
  bypass.x = x;
  bypass.obj = NULL;
  ClosestBBoxesToPoint(tree->root, x, TreePointDistanceCallback,
                       &bypass, &minmax);
  *obj = bypass.obj;
  return bypass.min;
}
