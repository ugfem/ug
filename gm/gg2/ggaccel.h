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

#define SMALLCOORD      1e-6
#define MAXCOORD 999.0
#define EPSI 0.0003
#define MAXWIDTH 8.0
#define MAXCOORD 999.0
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
  COORD x;
  COORD y;
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










static INT InitAccelObjs (MULTIGRID *theMG);

INT TerminateAccel (MULTIGRID *theMG, INT flag);

static float length_of_edge(FRONTCOMP* edgefc, FRONTCOMP* edgefc_aim);

static float length_of_angle(FRONTCOMP* anglefc_ori, FRONTCOMP* anglefc,
                             FRONTCOMP* anglefc_aim);

static INT PointInTriangle (COORD pt[DIM], COORD x[3], COORD y[3]);

static INT PointInCircle (COORD pt[DIM], COORD x, COORD y, COORD searchrad2);

static COORD distance(FRONTCOMP *p1, FRONTCOMP *p2);

static QUADTREETYP* search(QUADTREETYP *q_pointer, SOURCETYP *so, COORD *wi,
                           FRONTCOMP *p_new);

static void environment_search(INDEPFRONTLIST *theIFL,
                               QUADTREETYP *q_pointer, SOURCETYP* so,
                               FRONTCOMP *thefoundPoints[MAXNPOINTS],
                               FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS],
                               COORD wi,
                               SOURCETYP *search_sq_ld, SOURCETYP *search_sq_ru,
                               SOURCETYP *big_search_sq_ld,
                               SOURCETYP *big_search_sq_ru,
                               COORD xt[3], COORD yt[3], COORD searchradis,
                               int *foundpoints, int *ii);

static void insert(QFCLISTTYP *p_new, QUADTREETYP *q_place,
                   SOURCETYP *src, COORD wi);

static void delete_node(QUADTREETYP *q_pointer, FRONTCOMP *p_del, COORD width,
                        SOURCETYP *so, int *aufhoeren, QFCLISTTYP **nd_mem);

static void delete_quadtree(QUADTREETYP *q_pointer);

static void FCTreeUpdate(FRONTCOMP *fc, int n);

static void InsertQuadtree(FRONTCOMP *pon, int ncomp);

static void DELETE_ND (FRONTCOMP *delete_p);

static void btree_ins(FRONTCOMP *edgefc, float x, BALTREETYP **p, int *h);

static void balance1(BALTREETYP **p, int *h);

static void balance2(BALTREETYP **p, int *h);

static void del(BALTREETYP **r, int *h);

static void delete(FRONTCOMP* edgefc, float x, BALTREETYP **p, int *h);

static void BaseTreeUpdate(FRONTCOMP* P, FRONTCOMP* Q, FRONTCOMP* S,
                           int ch, int anglecrit, int edgecrit);

void AccelUpdate(FRONTCOMP* theFC, FRONTCOMP* thenewFC, FRONTCOMP* the_old_succ,
                 int cas, int anglecrit, int edgecrit);

int AccelInit(GRID *the_Grid, int anglecrit, int edgecrit, GG_PARAM *params);

int AccelFCTreeSearch(INDEPFRONTLIST *theIFL,
                      FRONTCOMP* thefoundPoints[MAXNPOINTS],
                      FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS],
                      COORD xt[3], COORD yt[3], COORD searchradis);

FRONTCOMP* AccelBaseTreeSearch(FRONTLIST** myList);
#endif
