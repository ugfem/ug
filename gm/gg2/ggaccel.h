// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ggaccel.h                                                         */
/*                                                                          */
/* Purpose:   header file for gg_accelerator                                */
/*                                                                          */
/* Author:    Dirk Feuchter													*/
/*                                                                          */
/* History:   17.05.95														*/
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __QUADTREE__
#define __QUADTREE__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include  "misc.h"
#include  "gm.h"
#include  "ugm.h"
#include "ggm.h"


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define SMALLDOUBLE     1e-6
#define MAXDOUBLE 999.0
#define EPSI 0.0003
#define MAXWIDTH 8.0
#define MAXDOUBLE 999.0
#define MAXNPOINTS 100
#define max(A,B)  ((A) > (B) ? (A) : (B))
#define min(A,B)  ((A) < (B) ? (A) : (B))

#define FROC(qfclp)             ((qfclp)->froco)
#define NXT(qfclp)              ((qfclp)->next)

#define NORMALCASE 0
#define LEFTNEIGHBOURCASE 1
#define RIGHTNEIGHBOURCASE 2
#define ININTERCASE 3
#define FINALCASE 4


/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

struct qfc_list_structure
{
  unsigned int control;
  struct qfc_list_structure *next;
  FRONTCOMP *froco;
};
typedef struct qfc_list_structure QFCLISTTYP;


typedef struct quadtree_structure
{
  unsigned int control;
  void* q_array[4];
  unsigned char q_flag;
}                                       QUADTREETYP;


typedef struct sourcetyp_structure
{
  unsigned int control;
  DOUBLE x;
  DOUBLE y;
}                                       SOURCETYP;

struct edgetree_structure
{
  struct edgetree_structure *left;
  struct edgetree_structure *right;
  int bal;
  FRONTCOMP *eFC;
  float key;
};
typedef struct edgetree_structure BALTREETYP;

INT TerminateAccel (MULTIGRID *theMG, INT flag);

void AccelUpdate(FRONTCOMP* theFC, FRONTCOMP* thenewFC, FRONTCOMP* the_old_succ,
                 int cas, int anglecrit, int edgecrit);

int AccelInit(GRID *the_Grid, int anglecrit, int edgecrit, GG_PARAM *params);

int AccelFCTreeSearch(INDEPFRONTLIST *theIFL,
                      FRONTCOMP* thefoundPoints[MAXNPOINTS],
                      FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS],
                      DOUBLE xt[3], DOUBLE yt[3], DOUBLE searchradis);

FRONTCOMP* AccelBaseTreeSearch(FRONTLIST** myList);
#endif
