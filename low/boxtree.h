// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      boxtree.h                                                     */
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

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __BOXTREE__
#define __BOXTREE__

#include "domain.h"
#include "dimension.h"

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef struct bt_object {
  double range[DIM][2];
  double left_range[2];
  double right_range[2];
  struct bt_object *left_son;
  struct bt_object *right_son;
} BT_OBJECT;

typedef struct {
  BT_OBJECT *root;
  double range[DIM][2];
} BOXTREE;

typedef void BT_FUNC (BT_OBJECT *, void *);

/****************************************************************************/
/*                                                                          */
/*  exported functions                                                      */
/*                                                                          */
/****************************************************************************/

void BT_Init(BT_OBJECT **box, int n, BOXTREE *tree);
void BT_Search(BOXTREE *tree, double range[DIM][2], BT_FUNC func, void *data);
void BT_SearchP(BOXTREE *tree, double *x, BT_FUNC func, void *data);

#endif
