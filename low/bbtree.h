// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*****************************************************************************
* File:      bbtree.h                                                       *
* Purpose:   bounding box tree data structure                               *
*                                                                           *
* Author:        O. Sterz                                                       *
*                                                                           *
* History:   Nov 2002 begin, ug3.8                                          *
* Remarks:   similar to boxtree.c, but more natural to compute distances    *
*****************************************************************************/

/*****************************************************************************
* auto include mechanism and other include files                            *
*****************************************************************************/
#ifndef __BBTREE__
#define __BBTREE__

#include "heaps.h"

#include "namespace.h"

START_NAMESPACE

/*****************************************************************************
* defines in the following order:                                           *
*        compile time constants defining static data size (i.e. arrays)     *
*        other constants                                                    *
*        macros                                                             *
*****************************************************************************/
#define MAX_DIM 3

/*
 * SWITCH:
 * Undefine the following line if your bboxes are not tight, i.e. if there
 * are bboxes containing objects that do not touch every side of the bounding
 * box. Note that here side means the side surface including its boundary.
 * This will result in somewhat slowered distance computations.
 */
#define TIGHT_BOXES

/*****************************************************************************
* exported data structures                                                  *
*****************************************************************************/
typedef struct bbt_bbox
{
  void *object;
  DOUBLE *ll;
  DOUBLE *ur;
} BBT_BBOX;

typedef struct bbt_node
{
  BBT_BBOX *bb;
  struct bbt_node *left;
  struct bbt_node *right;
} BBT_NODE;

typedef struct bbt_tree
{
  BBT_NODE *root;
  HEAP *heap;
  INT dim;
  INT n;
} BBT_TREE;

typedef void BBT_FUNC (void *obj, void *bypass);
typedef DOUBLE BBT_POINT_DIST_FUNC (DOUBLE *x, void *obj);


/*****************************************************************************
* exported global variables                                                 *
*****************************************************************************/

/*****************************************************************************
* public function declarations                                              *
*****************************************************************************/
BBT_TREE *BBT_NewTree(HEAP *theHeap, BBT_BBOX **bboxes, INT n, INT dim);
BBT_BBOX *BBT_NewBBox(HEAP *theHeap, INT dim, DOUBLE *ll, DOUBLE *ur,
                      void *obj);
void BBT_ClosestBBoxesToPoint(BBT_TREE *tree, DOUBLE *x, BBT_FUNC func,
                              void *bypass);
DOUBLE BBT_TreePointDistance(BBT_TREE *tree, DOUBLE *x, void **obj,
                             BBT_POINT_DIST_FUNC dist);

END_NAMESPACE

#endif
