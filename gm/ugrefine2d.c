// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugrefine2d.c													*/
/*																			*/
/* Purpose:   unstructured grid refinement (tree version)					*/
/*																			*/
/* Author:	  Peter Bastian, Nicolas Neuss									*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de                                               */
/*																			*/
/* History:   09.03.92 begin, ug version 2.0								*/
/*			  17.11.92 new rules, more efficient							*/
/*			  01.12.92 tree structure locally abandoned (N. Neuss)			*/
/*			  01.01.93 new and more flexible rules, better					*/
/*					   concept for irregular elements (N.Neuss)                     */
/*			  09.02.93 grid rebuilding included (N.Neuss)					*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment ugrefine
#endif

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "compiler.h"

#include "devices.h"

#include "switch.h"
#include "cw.h"
#include "gm.h"
#include "misc.h"
#include "evm.h"
#include "shapes2d.h"
#include "ugm.h"
#include "ugrefine.h"
#include "ugrefine2d.h"
#include "algebra.h"


#ifdef __THREEDIM__
#error this source file is for 2D ONLY
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAXBASICRULES   14                      /* number of basic refinement rules     */

#define MAXRULES                100             /* maximum number of refinement rules	*/

#define MINANGLE                10.0            /* min angle for switching from T4 to T3*/
#define MAXANGLE                170.0           /* max angle for switching from T4 to T3*/

#define BMNRESOLUTION   10                      /* resolution for creating boundary midn*/

#define MINVNCLASS              2                       /* determines copies, dep. on discr. !	*/

/* numbers of basic rules */
#define NO_REF                  0
#define T_COPY                  1
#define T_RED                   2
#define T_BISECT_1              3
#define T_BISECT_2_T1   4
#define T_BISECT_2_T2   5
#define T_BISECT_2_Q    6
#define T_BISECT_3              7

#define Q_COPY                  8
#define Q_RED                   9
#define Q_CLOSE_1               10
#define Q_BLUE                  11
#define Q_CLOSE_2               12
#define Q_CLOSE_3               13

/* swap rules */
#define Q0 0
#define Q1 1
#define Q2 2
#define TT0 3
#define TT1 4
#define TT2 5
#define QT0 6
#define QT1 7
#define QT2 8
#define QT3 9
#define QT4 10
#define QT5 11
#define QT6 12
#define QQ0 13
#define QQ1 14
#define QQ2 15
#define QQ3 16
#define QQ4 17
#define QQ5 18
#define QQ6 19
#define QQ7 20
#define QQ8 21

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {
  int tag;
  int corners[4];
  int outer[4];
  int nb[4];
} SON_DATA ;

typedef struct {
  int tag;
  int rule;
  int class;
  int variant;
  int nsons;
  int pattern[4];
  int condensedPattern;
  SON_DATA sons[4];
  int midSon[5],midCorner[5],nbSon[8],nbSide[8];
  int followRule[16];
} REFRULE;

typedef struct {
  int rule;
  int nsons;
  int class;
  SON_DATA sons[4];
} CHANGE_RULE;

typedef struct {
  INT tag;                                                      /* for which element type				*/
  INT offset;                                           /* offset of rule... is computed		*/
  INT classOffset[3];                           /* offset of class... is computed		*/
  INT variantOffset[4];                         /* offset of variant... is computed     */
  INT classes[3];                                       /* 1 if class available                                 */
  INT variants[4];                                      /* 1 if variant available				*/
  INT pattern[4];                                       /* 1 if side is refined                                 */
  INT nsons;                                                    /* number of sons						*/
  SON_DATA sons[4];
  INT followRule[16];           /* is equal to T/QStandardFollowRule if ==-1	*/
  INT followVariant[16];        /* is equal to T/QStandardFollowVariant if ==-1 */
} basicRule;

typedef struct {
  NODE *Nodes[4];                                       /* corner nodes of an element to be ref.*/
  NODE *MidNodes[4];                                    /* nodes on refined edges				*/
  NODE *Interior;                                       /* the interior node					*/
  ELEMENT *Neighbors[8];                        /* neighbors on next higher level		*/
  ELEMENT **BackPtrs[8];                        /* locations where neighbs store my adr */
} ElementContext ;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT rFlag=GM_REFINE_TRULY_LOCAL; /* type of refine					*/
static INT MaxRules=0;  /* actual number of rules generated from BasicRules */
static REFRULE Rules[MAXRULES];                 /* the generated rules				*/

static basicRule BasicRules[MAXBASICRULES] = {
  {3,0,{0,0,0},{0,0,0,0},{1,0,0},{1,0,0,0},{0,0,0,0},0,{                                        /* NO_REF=0  */
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {3,0,{0,0,0},{0,0,0,0},{1,1,1},{1,0,0,0},{0,0,0,0},1,{                                        /* T_COPY=1  */
     {3,{0,1,2,0},{1,1,1,0},{0,2,4,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {3,0,{0,0,0},{0,0,0,0},{0,0,1},{1,0,0,0},{1,1,1,0},4,{                                        /* T_RED=2	*/
     {3,{0,4,6,0},{1,0,1,0},{0,3,5,0}},
     {3,{4,1,5,0},{1,1,0,0},{1,2,3,0}},
     {3,{6,5,2,0},{0,1,1,0},{3,3,4,0}},
     {3,{5,6,4,0},{0,0,0,0},{2,0,1,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {3,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,0},{1,0,0,0},2,{                                        /* T_BISECT_1=3  */
     {3,{4,2,0,0},{0,1,1,0},{1,4,0,0}},
     {3,{4,1,2,0},{1,1,0,0},{1,2,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,T_BISECT_2_T1,-1,T_BISECT_2_T2,-1,T_BISECT_3,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,                        0,-1,                    2,-1,                 0,-1,-1,-1,-1,-1,-1,-1,-1}},
  {3,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,0},{1,1,0,0},3,{                                        /* T_BISECT_2_T1=4	*/
     {3,{4,2,0,0},{0,1,1,0},{1,4,0,0}},
     {3,{4,5,2,0},{0,1,0,0},{2,3,0,0}},
     {3,{4,1,5,0},{1,1,0,0},{1,2,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {3,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,0},{1,1,0,0},3,{                                        /* T_BISECT_2_T2=5	*/
     {3,{5,4,1,0},{0,1,1,0},{1,1,2,0}},
     {3,{4,5,0,0},{0,0,1,0},{0,2,0,0}},
     {3,{0,5,2,0},{0,1,1,0},{1,3,4,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {3,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,0},{1,1,0,0},2,{                                        /* T_BISECT_2_Q=6  */
     {3,{5,4,1,0},{0,1,1,0},{1,1,2,0}},
     {4,{4,5,2,0},{0,1,1,1},{0,3,4,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {3,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,0},{1,1,1,0},4,{                                        /* T_BISECT_3=7  */
     {3,{4,1,5,0},{1,1,0,0},{1,2,1,0}},
     {3,{4,5,2,0},{0,1,0,0},{0,3,2,0}},
     {3,{4,2,6,0},{0,1,0,0},{1,4,3,0}},
     {3,{4,6,0,0},{0,1,1,0},{2,5,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {4,0,{0,0,0},{0,0,0,0},{1,1,1},{1,0,0,0},{0,0,0,0},1,{                                        /* Q_COPY=8  */
     {4,{0,1,2,3},{1,1,1,1},{0,2,4,6}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {4,0,{0,0,0},{0,0,0,0},{0,0,1},{1,0,0,0},{1,1,1,1},4,{                                        /* Q_RED=9 */
     {4,{0,4,8,7},{1,0,0,1},{0,1,3,7}},
     {4,{4,1,5,8},{1,1,0,0},{1,2,2,0}},
     {4,{8,5,2,6},{0,1,1,0},{1,3,4,3}},
     {4,{7,8,6,3},{0,0,1,1},{0,2,5,6}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {4,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,1},{1,0,0,0},3,{                                        /* Q_CLOSE_1=10 */
     {3,{4,3,0,0},{0,1,1,0},{1,6,0,0}},
     {3,{4,2,3,0},{0,1,0,0},{2,4,0,0}},
     {3,{4,1,2,0},{1,1,0,0},{1,2,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {4,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,0,0},{1,0,1,0},2,{                                        /* Q_BLUE=11 */
     {4,{0,4,6,3},{1,0,1,1},{0,1,5,6}},
     {4,{4,1,2,6},{1,1,1,0},{1,2,4,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {4,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,1},{1,1,0,0},3,{                                        /* Q_CLOSE_2=12 */
     {4,{4,8,3,0},{0,0,1,1},{1,2,6,0}},
     {4,{4,1,5,8},{1,1,0,0},{1,2,2,0}},
     {4,{5,2,3,8},{1,1,0,0},{3,4,0,1}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}},
  {4,0,{0,0,0},{0,0,0,0},{0,1,1},{1,1,1,1},{1,1,1,0},4,{                                        /* Q_CLOSE_3=13 */
     {3,{5,4,1,0},{0,1,1,0},{3,1,2,0}},
     {3,{5,2,6,0},{1,1,0,0},{3,4,3,0}},
     {4,{4,6,3,0},{0,1,1,1},{3,5,6,0}},
     {3,{5,6,4,0},{0,0,0,0},{1,2,0,0}}
   },
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
   {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}}
} ;

/*static INT StandardFollowRule_T[8]	=	{T_COPY,T_BISECT_1,T_BISECT_1,T_BISECT_2_Q,T_BISECT_1,T_BISECT_2_Q,T_BISECT_2_Q,T_RED};*/
static INT StandardFollowRule_T[8]      =       {T_COPY,T_BISECT_1,T_BISECT_1,T_BISECT_2_T1,T_BISECT_1,T_BISECT_2_T1,T_BISECT_2_T1,T_RED};
static INT StandardFollowVariant_T[8]=  {         0,             0,             1,                        0,             2,                2,                    1,    0};
static INT StandardFollowRule_Q[16] =   {Q_COPY,Q_CLOSE_1,Q_CLOSE_1,Q_CLOSE_2,Q_CLOSE_1,Q_BLUE,Q_CLOSE_2,Q_CLOSE_3,Q_CLOSE_1,Q_CLOSE_2,Q_BLUE,Q_CLOSE_3,Q_CLOSE_2,Q_CLOSE_3,Q_CLOSE_3,Q_RED};
static INT StandardFollowVariant_Q[16]= {         0,            0,                1,            0,                2,     0,        1,            0,        3,            3,     1,                3,            2,                2,            1,        0};

static CHANGE_RULE CRules[]=
{
  {Q0,1,REGULAR_CLASS,
   {
     {4,{0,1,2,3},{1,1,1,1},{0,1,2,3}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {Q1,2,REGULAR_CLASS,
   {
     {3,{0,1,2,0},{1,1,0,0},{0,1,1,0}},
     {3,{0,2,3,0},{0,1,1,0},{0,2,3,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {Q2,2,REGULAR_CLASS,
   {
     {3,{0,1,3,0},{1,0,1,0},{0,1,3,0}},
     {3,{1,2,3,0},{1,1,0,0},{1,2,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {TT0,2,REGULAR_CLASS,
   {
     {3,{0,1,2,0},{0,1,1,0},{1,1,2,0}},
     {3,{4,5,6,0},{0,1,1,0},{0,5,6,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {TT1,1,REGULAR_CLASS,
   {
     {4,{5,6,1,2},{1,1,1,1},{5,6,1,2}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {TT2,2,REGULAR_CLASS,
   {
     {3,{2,5,6,0},{1,1,0,0},{2,5,1,0}},
     {3,{2,6,1,0},{0,1,1,0},{0,6,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QT0,2,REGULAR_CLASS,
   {
     {4,{0,1,2,3},{0,1,1,1},{1,1,2,3}},
     {3,{4,5,6,0},{0,1,1,0},{0,5,6,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QT1,2,REGULAR_CLASS,
   {
     {3,{0,2,3,0},{0,1,1,0},{1,2,3,0}},
     {4,{5,6,1,2},{1,1,1,0},{5,6,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QT2,2,REGULAR_CLASS,
   {
     {4,{5,6,2,3},{1,0,1,1},{5,1,2,3}},
     {3,{6,1,2,0},{1,1,0,0},{6,1,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QT3,3,REGULAR_CLASS,
   {
     {3,{0,2,3,0},{0,1,1,0},{1,2,3,0}},
     {3,{5,6,2,0},{1,0,0,0},{5,2,0,0}},
     {3,{6,1,2,0},{1,1,0,0},{6,1,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QT4,2,REGULAR_CLASS,
   {
     {3,{3,1,2,0},{0,1,1,0},{1,1,2,0}},
     {4,{4,3,5,6},{0,1,1,1},{0,3,5,6}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QT5,2,REGULAR_CLASS,
   {
     {4,{3,6,1,2},{0,1,1,1},{1,6,1,2}},
     {3,{3,5,6,0},{1,1,0,0},{3,5,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QT6,3,REGULAR_CLASS,
   {
     {3,{3,1,2,0},{0,1,1,0},{1,1,2,0}},
     {3,{3,6,4,0},{0,1,0,0},{2,6,0,0}},
     {3,{3,5,6,0},{1,1,0,0},{3,5,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QQ0,2,REGULAR_CLASS,
   {
     {4,{0,1,2,3},{0,1,1,1},{1,1,2,3}},
     {4,{4,5,6,7},{0,1,1,1},{0,5,6,7}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QQ1,3,REGULAR_CLASS,
   {
     {3,{0,2,3,0},{0,1,1,0},{1,2,3,0}},
     {4,{5,6,1,2},{1,0,1,0},{5,2,1,0}},
     {3,{6,7,4,0},{1,1,0,0},{6,7,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QQ2,4,REGULAR_CLASS,
   {
     {3,{0,2,3,0},{0,1,1,0},{1,2,3,0}},
     {3,{5,6,2,1},{1,0,0,0},{5,2,0,0}},
     {3,{6,1,2,0},{0,1,0,0},{3,1,1,0}},
     {3,{6,7,4,0},{1,1,0,0},{6,7,2,0}}
   }},
  {QQ3,3,REGULAR_CLASS,
   {
     {3,{0,2,3,0},{0,1,1,0},{1,2,3,0}},
     {4,{5,7,1,2},{0,1,1,0},{2,7,1,0}},
     {3,{5,6,7,0},{1,1,0,0},{5,6,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QQ4,4,REGULAR_CLASS,
   {
     {3,{0,2,3,0},{0,1,1,0},{1,2,3,0}},
     {3,{5,7,2,0},{0,0,0,0},{2,3,0,0}},
     {3,{5,6,7,0},{1,1,0,0},{5,6,1,0}},
     {3,{1,2,7,0},{1,0,1,0},{1,1,7,0}}
   }},
  {QQ5,3,REGULAR_CLASS,
   {
     {3,{3,1,2,0},{0,1,1,0},{1,1,2,0}},
     {4,{5,7,1,3},{0,1,0,1},{2,7,0,3}},
     {3,{5,6,7,0},{1,1,0,0},{5,6,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QQ6,4,REGULAR_CLASS,
   {
     {3,{3,1,2,0},{0,1,1,0},{1,1,2,0}},
     {3,{3,7,4,0},{0,1,0,0},{2,7,0,0}},
     {3,{3,0,7,0},{1,0,0,0},{3,3,1,0}},
     {3,{5,6,7,0},{1,1,0,0},{5,6,2,0}}
   }},
  {QQ7,3,REGULAR_CLASS,
   {
     {3,{3,1,2,0},{0,1,1,0},{1,1,2,0}},
     {4,{5,6,1,3},{1,0,0,1},{5,2,0,3}},
     {3,{6,7,4,0},{1,1,0,0},{6,7,1,0}},
     {0,{0,0,0,0},{0,0,0,0},{0,0,0,0}}
   }},
  {QQ8,4,REGULAR_CLASS,
   {
     {3,{3,1,2,0},{0,1,1,0},{1,1,2,0}},
     {3,{3,6,4,0},{0,0,0,0},{3,2,0,0}},
     {3,{6,7,4,0},{1,1,0,0},{6,7,1,0}},
     {3,{3,5,6,0},{1,1,0,0},{3,5,1,0}}
   }}
};

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* Function:  GetRulePtr													*/
/*																			*/
/* Purpose:   permits access to the refinement rules						*/
/*																			*/
/* Input:	  INT rule														*/
/*																			*/
/* Output:	  pointer to rule												*/
/*																			*/
/****************************************************************************/

static REFRULE *GetRulePtr(INT rule)
{
  return(&(Rules[rule]));
}

/****************************************************************************/
/*																			*/
/* Function:  GetCRulePtr													*/
/*																			*/
/* Purpose:   permits access to the swap rules								*/
/*																			*/
/* Input:	  INT rule														*/
/*																			*/
/* Output:	  pointer to rule												*/
/*																			*/
/****************************************************************************/

static CHANGE_RULE *GetCRulePtr(INT rule)
{
  return(&(CRules[rule]));
}

/****************************************************************************/
/*																			*/
/*							   Grid Refinement								*/
/*																			*/
/*	The purpose of the following routines is to compute a new sequence of	*/
/*	nested triangulations from an existing one. Since the currently finest	*/
/*	triangulation of the domain omega is stored distributed over several	*/
/*	levels the input of the algorithm must be provided on the elements,     */
/*	which have no regular sons. (the elements where EstimateHere(e)==TRUE). */
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* Function:  CloseElement													*/
/*																			*/
/* Purpose:   makes the new refinement marks consistent                                         */
/*																			*/
/* Param:	  ELEMENT *theElement: pointer to element						*/
/*																			*/
/* return:	  0: element will not be refined								*/
/*			  1: this element will be refined in some way					*/
/*			  -1: error occurred											*/
/****************************************************************************/

static INT CloseElement (ELEMENT *theElement)
{
  ELEMENT *theNeighbor;
  INT i,l,n,myRule,myBasicRule,nbRule,myClass,myVariant;
  INT myPattern,newPattern,nbPattern0,nbPattern1,regFlag;

  n = TAG(theElement);

  /* get my pattern from current refinement rule */
  myRule = MARK(theElement);
  myPattern = Rules[myRule].condensedPattern;
  myClass = Rules[myRule].class;

  /* check if rule can be applied */
  if (MARK(theElement)!=0)
    if (Rules[myRule].tag!=n)
    {
      PrintErrorMessage('E',"CloseElement","rule does not fit");
      return(-1);
    }

  if (ECLASS(theElement)==COPY_CLASS) return(0);        /* nothing should be done for copies */

  /* compute pattern the neighbors with regular refinement rule imply (here could be some optimizing!) */
  nbPattern0 = 0; nbPattern1=0; regFlag=0;
  for (i=n-1; i>=0; i--)
  {
    nbPattern0 = nbPattern0<<1;
    nbPattern1 = nbPattern1<<1;
    theNeighbor = NBELEM(theElement,i);
    if (theNeighbor!=NULL)
      if (ECLASS(theNeighbor)==REGULAR_CLASS)
        if (Rules[nbRule=MARK(theNeighbor)].rule!=NO_REF)
          if (Rules[nbRule].class==REGULAR_CLASS)
          {
            regFlag=1;
            for (l=0; l<TAG(theNeighbor); l++)
              if (NBELEM(theNeighbor,l)==theElement)
                break;

            if (Rules[nbRule].pattern[l]==0)
              nbPattern0++;
            else
              nbPattern1++;
          }
  }

  newPattern=(myPattern&(~nbPattern0))|nbPattern1;

  /* if patterns are not consistent change to new rule */
  if ((myPattern!=newPattern)||regFlag)
  {
    if (myRule==NO_REF)
    {
      if (n==TRIANGLE)
      {
        myBasicRule=StandardFollowRule_T[newPattern];
        myVariant=StandardFollowVariant_T[newPattern];
      }
      else
      {
        myBasicRule=StandardFollowRule_Q[newPattern];
        myVariant=StandardFollowVariant_Q[newPattern];
      }
    }
    else
    {
      myRule=Rules[myRule].followRule[newPattern];
      if (myRule==-1)
      {
        PrintErrorMessage('E',"CloseElement","follow rule has not been set");
        return(-1);
      }

      myBasicRule=Rules[myRule].rule;
      myVariant=Rules[myRule].variant;
    }

    /* if previous rule was not regular the rule is irregular */
    if (myClass!=REGULAR_CLASS)
      myRule=BasicRules[myBasicRule].classOffset[IRREGULAR_CLASS]+BasicRules[myBasicRule].variantOffset[myVariant];
    else
      myRule=BasicRules[myBasicRule].classOffset[REGULAR_CLASS]+BasicRules[myBasicRule].variantOffset[myVariant];

    SETMARK(theElement,myRule);
  }

  if ((ECLASS(theElement)==REGULAR_CLASS)&&(Rules[myRule].class==REGULAR_CLASS))
    return(1);
  else
    return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  CloseGrid                                                                                                         */
/*																			*/
/* Purpose:   compute closure for next level								*/
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  INT >0: elements will be refined								*/
/*			  INT 0: no elements will be refined							*/
/*																			*/
/****************************************************************************/

static INT CloseGrid (GRID *theGrid)
{
  ELEMENT *theElement;
  INT status,cnt;

  cnt = 0;

  /* go once through all elements */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if ((status=CloseElement(theElement))>=0)
      cnt+=status;
    else
      return(-1);

  return(cnt);
}


/****************************************************************************/
/*																			*/
/* Function:  ComputeNodeClasses											*/
/*																			*/
/* Purpose:   compute node classes from refinement							*/
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  0: ok                                                                                                                 */
/*																			*/
/****************************************************************************/

static INT ComputeDistance (NODE *theNode)
{
  LINK *link1,*link2;

  /* distance 0 */
  if (NCLASS(theNode)==3) return(0);

  /* search distance 1 */
  for (link1=START(theNode); link1!=NULL; link1=NEXT(link1))
    if (NCLASS(NBNODE(link1))==3) return(1);

  /* search distance 2 */
  for (link1=START(theNode); link1!=NULL; link1=NEXT(link1))
    for (link2=START(NBNODE(link1)); link2!=NULL; link2=NEXT(link2))
      if (NCLASS(NBNODE(link2))==3) return(2);

  /* no class 3 node in distance <= 2 */
  return(999);
}

static INT ComputeNodeClasses (GRID *theGrid)
{
  NODE *theNode;
  ELEMENT *theElement;
  INT i,flag;

  /* reset classes on next level */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
    SETNCLASS(theNode,0);

  /* set corners of irregularly and regularly refined elements to class 3 */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if (Rules[MARK(theElement)].class>=IRREGULAR_CLASS)
    {
      flag=1;
      for (i=0; i<TAG(theElement); i++)
        SETNCLASS(CORNER(theElement,i),3);
    }

  if (rFlag==GM_COPY_ALL)
  {
    if (flag)
      for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
        SETNCLASS(theNode,3);
  }
  else
  {
    /* find class 2 nodes, this implementation better for parallel version */
    for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
      SETNCLASS(theNode,MAX(0,3-ComputeDistance(theNode)));
  }

  return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  ComputeCopies                                                                                             */
/*																			*/
/* Purpose:   determine copy elements from node classes                                         */
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  0: ok                                                                                                                 */
/*																			*/
/****************************************************************************/

static INT ComputeCopies (GRID *theGrid)
{
  ELEMENT *theElement;
  INT i;

  /* determine node classes in any version, since it doesn't hurt */
  ComputeNodeClasses(theGrid);

  /* now this is only correct in the old version				*/
  /* all unrefined elements with a corner >= 2 are set to copy*/
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if (MARK(theElement)==0)
    {
      for (i=0; i<TAG(theElement); i++)
        if (NCLASS(CORNER(theElement,i))>=2)
        {
          if (TAG(theElement)==TRIANGLE)
            SETMARK(theElement,BasicRules[T_COPY].classOffset[COPY_CLASS]);
          else
            SETMARK(theElement,BasicRules[Q_COPY].classOffset[COPY_CLASS]);
          break;
        }
    }

  return(0);
}

static INT Version3_ComputeCopies (GRID *theGrid)
{
  ELEMENT *theElement;
  INT flag;

  /* determine node classes in any version, since it doesn't hurt */
  ComputeNodeClasses(theGrid);

  /* set class of all dofs on next level to 0 */
  ClearNextVectorClasses(theGrid);

  /* seed dofs of regularly and irregularly refined elements to 3 */
  flag = 0;
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if (Rules[MARK(theElement)].class>=IRREGULAR_CLASS)
    {
      SeedNextVectorClasses(theElement);
      flag=1;                   /* there is at least one element to be refined */
    }

  /* copy all option or neighborhood */
  if (rFlag==GM_COPY_ALL)
  {
    if (flag)
      for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
        SeedNextVectorClasses(theElement);
  }
  else
  {
    PropagateNextVectorClasses(theGrid);
  }

  /* an element is copied if it has a dof of class 2 and higher */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if ((MARK(theElement)==0)&&(MaxNextVectorClass(theElement)>=MINVNCLASS))
    {
      if (TAG(theElement)==TRIANGLE)
        SETMARK(theElement,BasicRules[T_COPY].classOffset[COPY_CLASS]);
      else
        SETMARK(theElement,BasicRules[Q_COPY].classOffset[COPY_CLASS]);
    }

  return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  RestrictMarks                                                                                             */
/*																			*/
/* Purpose:   restrict refinement marks when going down                                         */
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void RestrictMarks (GRID *theGrid)
{
  ELEMENT *theElement,*theNeighbor;
  INT myClass,sonClass,myRule,sonRule,nbRule,newRule,newVariant;
  INT newPattern,nbPattern0,nbPattern1,orPattern;
  INT i,l,m,flag;

  /* this loop handles COARSEN commands and sets MARK=OLDREF for regularly refined elements */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    myRule = REFINE(theElement);
    if (myRule==NO_REF) continue;               /* leaf elements are marked by estimator */

    sonClass = Rules[myRule].class;
    myClass = ECLASS(theElement);

    /* copies that have a copy son take the tag of the son */
    if (myClass==COPY_CLASS)
    {
      SETMARK(theElement,MARK(SON(theElement,0)));
      SETCOARSEN(theElement,COARSEN(SON(theElement,0)));
      continue;
    }

    /* irregular elements that have a copy son take the tag of the son but cannot be coarsened */
    if (myClass==IRREGULAR_CLASS)
    {
      SETMARK(theElement,MARK(SON(theElement,0)));
      SETCOARSEN(theElement,0);
      continue;
    }

    /* now we are left with regular elements that are refined */

    /* son is a copy, take tag of the son */
    if (sonClass==COPY_CLASS)
    {
      SETMARK(theElement,MARK(SON(theElement,0)));
      SETCOARSEN(theElement,COARSEN(SON(theElement,0)));
      continue;
    }

    /* regular elements with irregular refinement are handled below */
    if (sonClass==IRREGULAR_CLASS) continue;

    /* regular elements with regular refinement are the only ones to coarsen */
    SETMARK(theElement,REFINE(theElement));
    if (DECOUPLED(theElement)) continue;                /* tree has been broken up here: no coarsening possible! */

    flag = 0;
    for (i=0; i<NSONS(theElement); i++)
    {
      /* if not all sons are marked no unrefinement is possible */
      if (!COARSEN(SON(theElement,i)))
      {
        flag = 1;
        break;
      }
    }

    if (!flag)
    {
      /* remove refinement */
      SETMARK(theElement,0);
      SETCOARSEN(theElement,0);
    }
  }

  /* now a second loop over all elements handles the regular elements with irregular sons */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    myRule = REFINE(theElement);
    if (myRule==NO_REF) continue;

    sonClass = Rules[myRule].class;
    if (sonClass!=IRREGULAR_CLASS) continue;                    /* only regular elements are allowed to have irregular sons */

    /* is there a refinement mark on the regular father? */
    if ((newRule=MARK(theElement))!=0)
      flag=1;
    else
    {
      /* no: check if its sons have a refinement flag set */
      flag=0;
      for (i=0; i<NSONS(theElement); i++)
        if ((sonRule=MARK(SON(theElement,i)))!=0)
        {
          flag=1;
          break;
        }
    }

    if (flag==0) continue;

    /* get the enforced pattern */
    nbPattern0=0; nbPattern1=0; orPattern=0;
    for (i=TAG(theElement)-1; i>=0; i--)
    {
      nbPattern0 = nbPattern0<<1;
      nbPattern1 = nbPattern1<<1;
      orPattern = orPattern<<1;
      theNeighbor = NBELEM(theElement,i);
      if (theNeighbor!=NULL)
        if (ECLASS(theNeighbor)==REGULAR_CLASS)
          if ((nbRule=MARK(theNeighbor))!=NO_REF)
            if (Rules[nbRule].class==REGULAR_CLASS)
            {
              if (Rules[nbRule].rule==Q_BLUE)
                orPattern++;
              m = TAG(theNeighbor);
              for (l=0; l<m; l++)
                if (NBELEM(theNeighbor,l)==theElement)
                  break;

              if (Rules[nbRule].pattern[l]==0)
                nbPattern0++;
              else
                nbPattern1++;
            }
    }

    if (newRule!=0)
    {
      /* there is a refinement order given: try to fulfill it as good as possible */
      newPattern=(Rules[newRule].condensedPattern&(~nbPattern0)) | nbPattern1;
      newRule=Rules[newRule].followRule[newPattern];
      SETMARK(theElement,newRule);
      continue;
    }

    /*	no refinement given on the father, yet we must change the refinement to regular:
            try to make it as good as possible (i.e. copy, red or blue) */
    if (TAG(theElement)==TRIANGLE)
    {
      if (nbPattern0==0)
      {
        newRule=T_RED;
        newVariant=0;
      }
      else
      {
        if (nbPattern1==0)
        {
          /* order copy refinement */
          newRule=T_COPY;
          newVariant=0;
        }
        else
        {
          /* no chance to make it more regular */
          newRule=StandardFollowRule_T[nbPattern1];
          newVariant=StandardFollowVariant_T[nbPattern1];
        }
      }
    }
    else
    {
      newRule=NO_REF;
      if ((orPattern==(nbPattern0|nbPattern1)) && (nbPattern0!=0))
      {
        /* all regular neighbors are blue refined and we have the possibility of locality by blue refinement */
        if (((5&(~nbPattern0))|nbPattern1)==5)
        {
          newRule=Q_BLUE;
          newVariant=0;
        }
        else
        {
          if (((10&(~nbPattern0))|nbPattern1)==10)
          {
            newRule=Q_BLUE;
            newVariant=1;
          }
        }
      }

      if (newRule==NO_REF)
      {
        if (nbPattern0==0)
        {
          newRule=Q_RED;
          newVariant=0;
        }
        else
        {
          if (nbPattern1==0)
          {
            newRule=Q_COPY;
            newVariant=0;
          }
          else
          {
            newRule=StandardFollowRule_Q[nbPattern1];
            newVariant=StandardFollowVariant_Q[nbPattern1];
          }
        }
      }

    }

    newRule=BasicRules[newRule].classOffset[REGULAR_CLASS]+BasicRules[newRule].variantOffset[newVariant];
    SETMARK(theElement,newRule);
    continue;

  }             /* endfor */

}


/****************************************************************************/
/*																			*/
/* Function:  CheckMemoryRequirements										*/
/*																			*/
/* Purpose:   check if there is enough memory for the following refinement	*/
/*																			*/
/* Param:	  MULTIGRID *theMG: pointer to multigrid structure with                 */
/*								computed refinement                                             */
/*																			*/
/* return:	  INT 1: enough memory for the operation						*/
/*			  INT 0: not enough memory for the operation					*/
/*																			*/
/****************************************************************************/

static INT CheckMemoryRequirements (MULTIGRID *theMG)
{
  if (theMG->theHeap->size-theMG->theHeap->used<=0)
    return(0);
  else
    return(1);
}


/****************************************************************************/
/*																			*/
/* Function:  CreateMidNode                                                                                             */
/*																			*/
/* Purpose:   allocate a new node on a side of an element. Includes vertex	*/
/*			  best fit boundary coordinates and local coordinates			*/
/*			  insert also links to endpoints of the refined edge			*/
/*																			*/
/* Param:	  ELEMENT *theElement: element to refine						*/
/*			  INT side: side to refine										*/
/*			  NODE *after: insert new node after that node					*/
/*																			*/
/* return:	  NODE* : pointer to new node									*/
/*			  NULL	: could not allocate									*/
/*																			*/
/****************************************************************************/

static NODE *CreateMidNode (GRID *theGrid,ELEMENT *theElement,INT side,NODE *after)
{
  ELEMENTSIDE *theSide;
  COORD x,y;
  COORD r[2],lambda1,lambda0,z;
  COORD lambda,dlambda,s,lambdaopt,smin;
  INT i,n;
  VERTEX *theVertex;
  VSEGMENT *vsnew;
  NODE *theNode;
  BNDSEGDESC *theSeg;

  n = TAG(theElement);

  /* calculate midpoint of edge */
  x = 0.5*(XC(MYVERTEX(CORNER(theElement,side)))+XC(MYVERTEX(CORNER(theElement,(side+1)%n))));
  y = 0.5*(YC(MYVERTEX(CORNER(theElement,side)))+YC(MYVERTEX(CORNER(theElement,(side+1)%n))));

  /* allocate vertex */
  if ((OBJT(theElement)==BEOBJ)&&(SIDE(theElement,side)!=NULL))
  {
    /* find optimal boundary coordinate for boundary vertex */
    theSide = SIDE(theElement,side);
    theSeg = SEGDESC(theSide);
    smin = 1.0E30;
    lambda0 = PARAM(theSide,0,0);
    lambda1 = PARAM(theSide,1,0);
    dlambda = (lambda1-lambda0)/((COORD) BMNRESOLUTION);
    lambda = lambda0;
    for (i=1; i<BMNRESOLUTION; i++)
    {
      lambda += dlambda;
      BNDSEGFUNC(theSeg) (BNDDATA(theSeg),&lambda,r);
      s = (r[0]-x)*(r[0]-x)+(r[1]-y)*(r[1]-y);
      if (s<smin)
      {
        smin = s;
        lambdaopt = lambda;
      }
    }

    /* create vertex */
    theVertex = CreateBoundaryVertex(theGrid,NULL);
    if (theVertex==NULL) return(NULL);
    vsnew = CreateVertexSegment(theGrid, theVertex);
    if (vsnew==NULL)
    {
      DisposeVertex(theGrid,theVertex);
      UserWrite("cannot create vertexsegment\n");
      return(NULL);
    }
    LAMBDA(vsnew,0) = lambdaopt;
    BNDSEGFUNC(theSeg) (BNDDATA(theSeg),PVECT(vsnew),CVECT(theVertex));
    z = (lambdaopt-lambda0)/(lambda1-lambda0);
    ZETA(vsnew) = z;
    SETONEDGE(theVertex,side);
    BSEGDESC(vsnew) = theSeg;
    if (n==3)
    {
      if (side==0) { XI(theVertex)= z  ; ETA(theVertex)= 0.0; }
      if (side==1) { XI(theVertex)= 1-z; ETA(theVertex)= z  ; }
      if (side==2) { XI(theVertex)= 0.0; ETA(theVertex)= 1-z; }
    }
    else
    {
      if (side==0) { XI(theVertex)= 2*z-1; ETA(theVertex)= -1.0; }
      if (side==1) { XI(theVertex)= 1.0;       ETA(theVertex)= 2*z-1; }
      if (side==2) { XI(theVertex)= 1-2*z; ETA(theVertex)= 1.0; }
      if (side==3) { XI(theVertex)= -1.0;  ETA(theVertex)= 1-2*z; }
    }
    VFATHER(theVertex) = theElement;
    SETMOVE(theVertex,1);
  }
  else
  {
    /* we need an inner vertex */
    theVertex = CreateInnerVertex(theGrid,NULL);
    if (theVertex==NULL) return(NULL);
    XC(theVertex) = x;
    YC(theVertex) = y;
    if (n==3)
    {
      if (side==0) { XI(theVertex)= 0.5; ETA(theVertex)= 0.0; }
      if (side==1) { XI(theVertex)= 0.5; ETA(theVertex)= 0.5; }
      if (side==2) { XI(theVertex)= 0.0; ETA(theVertex)= 0.5; }
    }
    else
    {
      if (side==0) { XI(theVertex)= 0.0; ETA(theVertex)=-1.0; }
      if (side==1) { XI(theVertex)= 1.0; ETA(theVertex)= 0.0; }
      if (side==2) { XI(theVertex)= 0.0; ETA(theVertex)= 1.0; }
      if (side==3) { XI(theVertex)=-1.0; ETA(theVertex)= 0.0; }
    }
    VFATHER(theVertex) = theElement;
    SETMOVE(theVertex,2);
  }

  /* allocate node */
  theNode = CreateNode(theGrid,after);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  MYVERTEX(theNode) = theVertex;
  NFATHER(theNode) = NULL;
  TOPNODE(theVertex) = theNode;

  return(theNode);
}


/****************************************************************************/
/*																			*/
/* Function:  GetCurrentContext                                                                                         */
/*																			*/
/* Purpose:   assemble references to objects which interact with the sons	*/
/*			  of the given element, as indicated by REFINE.                                 */
/*			  (i)	 corner nodes											*/
/*			  (ii)	 nodes at midpoints of edges							*/
/*			  (iii)  neighboring elements									*/
/*			  (iv)	 locations where neighbors store references to our sons */
/*																			*/
/* Param:	  ELEMENT *theElement: element to refine						*/
/*			  ElementContext *theContext: context structure to fill                 */
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void GetCurrentContext (ELEMENT *theElement, ElementContext *theContext)
{
  NODE *theNode;                                        /* corner node							*/
  ELEMENT *theNeighbor,*theSon;         /* neighbor and a son of current elem.	*/
  ELEMENT *theNeighborOfSon;                    /* a neighbor of theSon                                 */
  INT i,l,m,n;                                          /* some integer variables				*/
  NODE **Nodes;                                         /* corner nodes of an element to be ref.*/
  NODE **MidNodes;                                      /* nodes on refined edges				*/
  ELEMENT **Neighbors;                          /* neighbors on next higher level		*/
  ELEMENT ***BackPtrs;                                  /* locations where neighbs store my adr */
  REFRULE *rule;                                        /* current refinement rule of theElement*/
  REFRULE *nbrule;                                      /* current refinement rule of neighbor	*/

  n = TAG(theElement);
  rule = &(Rules[REFINE(theElement)]);

  Nodes = theContext->Nodes;
  MidNodes = theContext->MidNodes;
  Neighbors = theContext->Neighbors;
  BackPtrs = theContext->BackPtrs;
  theContext->Interior = NULL;

  for (i=0; i<4; i++)                           /* reset pointers						*/
  {
    Nodes[i] = MidNodes[i] = NULL;
    Neighbors[2*i] = Neighbors[2*i+1] = NULL;
    BackPtrs[2*i] = BackPtrs[2*i+1] = NULL;
  }

  for (i=0; i<n; i++)
  {                                                                     /* get corner nodes                                     */

    theNode = CORNER(theElement,i);
    Nodes[i] = SONNODE(theNode);
  }

  for (i=0; i<n; i++)
  {                                                                     /* get midpoint nodes					*/
    if (rule->midSon[i]>=0)
      MidNodes[i] = CORNER(SON(theElement,rule->midSon[i]),rule->midCorner[i]);
  }

  /* get interior node */
  if ((n==4)&&(rule->midSon[4]>=0))
    theContext->Interior = CORNER(SON(theElement,rule->midSon[4]),rule->midCorner[4]);

  for (i=0; i<n; i++)
  {                                                                     /* get pointers to neighbors and locatio*/
                                                                        /* of their pointers to us				*/

    theNeighbor = NBELEM(theElement,i);

    if (theNeighbor!=NULL)
    {
      /*	If neighbor is decoupled then get information from the own sons via son->nb.
              This is essential when abondoning locally the tree structure.	*/

      if (DECOUPLED(theNeighbor))
      {
        if (rule->nbSon[2*i]>=0)
        {
          theSon=SON(theElement,rule->nbSon[2*i]);
          Neighbors[2*i]=theNeighborOfSon=NBELEM(theSon,rule->nbSide[2*i]);
          if (theNeighborOfSon!=NULL)
          {
            /* now we have to search for sonEl */
            for (l=0; l<TAG(theNeighborOfSon); l++)
              if (NBELEM(theNeighborOfSon,l)==theSon)
                break;
            BackPtrs[2*i]  = (ELEMENT **) &(VOID_NBELEM(theNeighborOfSon,l));
          }

          if (rule->nbSon[2*i+1]>=0)
          {
            theSon=SON(theElement,rule->nbSon[2*i+1]);
            Neighbors[2*i+1]=theNeighborOfSon=NBELEM(theSon,rule->nbSide[2*i+1]);
            if (theNeighborOfSon!=NULL)
            {
              /* now we have to search for sonEl */
              for (l=0; l<TAG(theNeighborOfSon); l++)
                if (NBELEM(theNeighborOfSon,l)==theSon)
                  break;
              BackPtrs[2*i+1]  = (ELEMENT **) &(VOID_NBELEM(theNeighborOfSon,l));
            }
          }
        }
        continue;                               /* next side of theElement */
      }

      /* else get the neighbors of the sons from the tree structure via nb->son */

      m = TAG(theNeighbor);
      for (l=0; l<m; l++)
        if (NBELEM(theNeighbor,l)==theElement)
          break;
      nbrule = &(Rules[REFINE(theNeighbor)]);

      if (nbrule->nbSon[2*l]>=0)
      {
        if (nbrule->nbSon[2*l+1]>=0)
        {
          Neighbors[2*i] = SON(theNeighbor,nbrule->nbSon[2*l+1]);
          if (Neighbors[2*i]!=NULL)
            BackPtrs[2*i]  = (ELEMENT **) &(VOID_NBELEM(Neighbors[2*i],nbrule->nbSide[2*l+1]));
          Neighbors[2*i+1] = SON(theNeighbor,nbrule->nbSon[2*l]);
          if (Neighbors[2*i+1]!=NULL)
            BackPtrs[2*i+1]  = (ELEMENT **) &(VOID_NBELEM(Neighbors[2*i+1],nbrule->nbSide[2*l]));
        }
        else
        {
          Neighbors[2*i] = SON(theNeighbor,nbrule->nbSon[2*l]);
          if (Neighbors[2*i]!=NULL)
            BackPtrs[2*i]  = (ELEMENT **) &(VOID_NBELEM(Neighbors[2*i],nbrule->nbSide[2*l]));
        }
      }
    }
  }

  return;
}

/****************************************************************************/
/*																			*/
/* Function:  UpdateContext                                                                                             */
/*																			*/
/* Purpose:   assemble references to objects which interact with the sons	*/
/*			  of the given element, i.e.									*/
/*			  objects are allocated, kept or deleted as indicated by MARK	*/
/*			  (i)	 corner nodes											*/
/*			  (ii)	 nodes at midpoints of edges							*/
/*			  (iii)  neighboring elements									*/
/*			  (iv)	 locations where neighbors store references to our sons */
/*																			*/
/* Param:	  GRID *theGrid: grid level of the sons of theElement			*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ElementContext *theContext: context structure to update		*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*			  INT 1: fatal memory error                                                                     */
/*																			*/
/****************************************************************************/

static INT UpdateContext (GRID *theGrid, ELEMENT *theElement, ElementContext *theContext)
{
  NODE *theNode;                                        /* corner node							*/
  ELEMENT *theNeighbor;                         /* neighbor and a son of current elem.	*/
  EDGE *theEdge;                                        /* temporary storage for an edge		*/
  INT i,l,m,n;                                          /* some integer variables				*/
  NODE **Nodes;                                         /* corner nodes of an element to be ref.*/
  NODE **MidNodes;                                      /* nodes on refined edges				*/
  ELEMENT **Neighbors;                          /* neighbors on next higher level		*/
  ELEMENT ***BackPtrs;                          /* locations where neighbs store my adr */
  LINK *theLink;                                        /* scan through nodes neighbor list     */
  INT candelete;
  REFRULE *rule;                                        /* new refinement rule of theElement	*/
  REFRULE *nbrule;                                      /* refinement rule of neighbor			*/

  n = TAG(theElement);
  Nodes = theContext->Nodes;
  MidNodes = theContext->MidNodes;
  Neighbors = theContext->Neighbors;
  BackPtrs = theContext->BackPtrs;
  rule = &(Rules[MARK(theElement)]);

  /* theContext holds current context */

  /* allocate corner nodes if necessary */
  if (MARK(theElement)>0)
  {
    /* we need corner nodes */
    for (i=0; i<n; i++)
    {
      theNode = CORNER(theElement,i);
      if (SONNODE(theNode)==NULL)
      {
        SONNODE(theNode) = CreateNode(theGrid,NULL);
        if (SONNODE(theNode)==NULL) return(1);
        SETUSED(theNode,1);
        theGrid->status |= 1;
        MYVERTEX(SONNODE(theNode)) = MYVERTEX(theNode);
        NFATHER(SONNODE(theNode)) = theNode;
        Nodes[i] = SONNODE(theNode);
        SETCLASS(Nodes[i],NCLASS(theNode));
      }
    }
  }

  /* allocate,keep, or delete midpoint nodes */
  /* allocate,keep, or delete corner/corner edges */
  for (i=0; i<n; i++)
  {
    if (rule->pattern[i])
    {
      /* if a corner corner edge exists then delete it */
      if ((theEdge = GetEdge(Nodes[i],Nodes[(i+1)%n]))!=NULL)
        DisposeEdge(theGrid,theEdge);

      /* we need a midpoint node */
      if (MidNodes[i]!=NULL) continue;

      /* lets see if neighbor has already one */
      theNeighbor = NBELEM(theElement,i);
      if (theNeighbor!=NULL)
      {
        m = TAG(theNeighbor);
        for (l=0; l<m; l++)
          if (NBELEM(theNeighbor,l)==theElement)
            break;
        nbrule = &(Rules[REFINE(theNeighbor)]);
        if (nbrule->midSon[l]>=0)
        {
          MidNodes[i] = CORNER(SON(theNeighbor,nbrule->midSon[l]),nbrule->midCorner[l]);
          continue;
        }
      }

      /* now we must allocate a new one */
      MidNodes[i] = CreateMidNode(theGrid,theElement,i,Nodes[i]);
      if (MidNodes[i]==NULL) return(1);
      SETCLASS(MidNodes[i],4);
      SETUSED(MidNodes[i],1);
      theGrid->status |= 1;
      if (CreateEdge(theGrid,Nodes[i],MidNodes[i])==NULL) return(1);
      if (CreateEdge(theGrid,Nodes[(i+1)%n],MidNodes[i])==NULL) return(1);
    }
    else
    {
      /* if we need a corner corner edge then allocate it */
      if (MARK(theElement)>0)
      {
        if ((theEdge = GetEdge(Nodes[i],Nodes[(i+1)%n]))==NULL)
          theEdge = CreateEdge(theGrid,Nodes[i],Nodes[(i+1)%n]);
      }
      else
      {
        /* delete corner corner edge */
        if (Neighbors[2*i]==NULL)
        {
          theEdge = GetEdge(Nodes[i],Nodes[(i+1)%n]);
          if (theEdge!=NULL) DisposeEdge(theGrid,theEdge);
        }
      }

      /* we don't need a midpoint node on that edge, lets see if it can be deleted */
      if (MidNodes[i]==NULL) continue;

      candelete = 1;
      for (theLink=START(MidNodes[i]); theLink!=NULL; theLink=NEXT(theLink))
      {
        if (!EXTRA(theLink))
          if ((NBNODE(theLink)!=Nodes[i])&&(NBNODE(theLink)!=Nodes[(i+1)%n]))
          {
            candelete = 0;
            break;
          }
      }

      if (candelete)
      {
        DisposeVertex(theGrid,MYVERTEX(MidNodes[i]));
        DisposeNode(theGrid,MidNodes[i]);
        MidNodes[i] = NULL;
      }
    }
  }

  /* delete corner nodes if possible */
  if (MARK(theElement)==0)
  {
    /* we need no corner nodes */
    for (i=0; i<n; i++)
    {
      theNode = CORNER(theElement,i);
      if (Nodes[i]!=NULL)
      {
        candelete = 1;
        for (theLink=START(Nodes[i]); theLink!=NULL; theLink=NEXT(theLink))
          if (!EXTRA(theLink))
          {
            candelete = 0;
            break;
          }
        if (candelete)
        {
          DisposeNode(theGrid,Nodes[i]);
          SONNODE(theNode) = NULL;
          Nodes[i] = NULL;
        }
      }
    }
  }

  /* get neighbor refernces for following refine */
  for (i=0; i<n; i++)
  {
    /* get pointers to neighbors and location of their pointers to us */
    theNeighbor = NBELEM(theElement,i);
    if (theNeighbor!=NULL)
    {
      /* in the following case it is still useless to set the pointers */
      if (REFINE(theNeighbor)!=MARK(theNeighbor)) continue;

      /*	if theNeighbor is still in the same state nothing should have changed
              for the context */
      if (THEFLAG(theNeighbor)==0)
        if (ECLASS(theNeighbor)==REGULAR_CLASS)
          if (Rules[REFINE(theNeighbor)].class==REGULAR_CLASS)
            continue;

      /* get information about new neighborhood from tree */
      m = TAG(theNeighbor);
      for (l=0; l<m; l++)
        if (NBELEM(theNeighbor,l)==theElement)
          break;
      nbrule = &(Rules[REFINE(theNeighbor)]);
      if (nbrule->nbSon[2*l]>=0)
      {
        if (nbrule->nbSon[2*l+1]>=0)
        {
          Neighbors[2*i] = SON(theNeighbor,nbrule->nbSon[2*l+1]);
          if (Neighbors[2*i]!=NULL)
            BackPtrs[2*i]  = (ELEMENT **) &(VOID_NBELEM(Neighbors[2*i],nbrule->nbSide[2*l+1]));
          else
            BackPtrs[2*i]=NULL;

          Neighbors[2*i+1] = SON(theNeighbor,nbrule->nbSon[2*l]);
          if (Neighbors[2*i+1]!=NULL)
            BackPtrs[2*i+1]  = (ELEMENT **) &(VOID_NBELEM(Neighbors[2*i+1],nbrule->nbSide[2*l]));
          else
            BackPtrs[2*i+1]=NULL;
        }
        else
        {
          Neighbors[2*i] = SON(theNeighbor,nbrule->nbSon[2*l]);
          if (Neighbors[2*i]!=NULL)
            BackPtrs[2*i]  = (ELEMENT **) &(VOID_NBELEM(Neighbors[2*i],nbrule->nbSide[2*l]));
          else
            BackPtrs[2*i]=NULL;

          Neighbors[2*i+1]=NULL;
          BackPtrs[2*i+1]=NULL;
        }
      }
      else
      {
        Neighbors[2*i]=Neighbors[2*i+1]=NULL;
        BackPtrs[2*i]=BackPtrs[2*i+1]=NULL;
      }
    }
  }

  return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  UnrefineElement												*/
/*																			*/
/* Purpose:   remove previous refinement of an element						*/
/*			  (i)	 all interior nodes and edges are deletes				*/
/*			  (ii)	 sons are deleted and references to sons reset to NULL	*/
/*																			*/
/* Param:	  GRID *theGrid: grid level of sons of theElement				*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ElementContext *theContext: current context of element		*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void UnrefineElement (GRID *theGrid, ELEMENT *theElement, ElementContext *theContext)
{
  INT i,n,m,s;
  EDGE *theEdge;
  NODE *Nodes[9];                                       /* corner nodes of an element to be ref.*/
  ELEMENT **Neighbors;                          /* neighbors on next higher level		*/
  ELEMENT ***BackPtrs;                          /* locations where neighbs store my adr */
  REFRULE *rule;                                        /* current refinement rule of theElement*/
  REFRULE *newrule;                                     /* new refinement rule of theElement	*/
  ElementContext sonContext;
  ELEMENT *theSon;

  /* something to do ? */
  if ((REFINE(theElement)==0)||(theGrid==NULL)) return;

  /* remove elements above my sons (recursively!) */
  for (s=0; s<NSONS(theElement); s++)
  {
    theSon = SON(theElement,s);
    SETMARK(theSon,0);
    if (REFINE(theSon)>0)
    {
      GetCurrentContext(theSon,&sonContext);
      UnrefineElement(theGrid->finer,theSon,&sonContext);
      UpdateContext(theGrid->finer,theSon,&sonContext);
    }
    SETREFINE(theSon,MARK(theSon));
  }

  n = TAG(theElement);

#ifdef __version3__
  /* remove connections in neighborhood of sons */
  for (i=0; i<NSONS(theElement); i++)
    DisposeConnectionsInNeighborhood(theGrid,SON(theElement,i));
#endif

  /* remove my sons */
  Neighbors = theContext->Neighbors;
  BackPtrs = theContext->BackPtrs;
  for (i=0; i<9; i++) Nodes[i] = NULL;
  for (i=0; i<n; i++)
  {
    Nodes[i] = theContext->Nodes[i];
    Nodes[i+4] = theContext->MidNodes[i];
  }
  if (n==4) Nodes[8] = theContext->Interior;

  rule = &(Rules[REFINE(theElement)]);

  /* remove interior edges */
  for (s=0; s<rule->nsons; s++)
  {
    m = rule->sons[s].tag;
    for (i=0; i<m; i++)
      if ((!rule->sons[s].outer[i])&&(rule->sons[s].nb[i]>s))
      {
        theEdge = GetEdge(Nodes[rule->sons[s].corners[i]],Nodes[rule->sons[s].corners[(i+1)%m]]);
        if (theEdge!=NULL) DisposeEdge(theGrid,theEdge);
      }
    if (m==4)
    {
      theEdge = GetEdge(Nodes[rule->sons[s].corners[0]],Nodes[rule->sons[s].corners[2]]);
      if (theEdge!=NULL) DisposeEdge(theGrid,theEdge);
      theEdge = GetEdge(Nodes[rule->sons[s].corners[1]],Nodes[rule->sons[s].corners[3]]);
      if (theEdge!=NULL) DisposeEdge(theGrid,theEdge);
    }
  }

  /* remove son elements */
  for (i=0; i<NSONS(theElement); i++)
    DisposeElement(theGrid,SON(theElement,i));
  SETNSONS(theElement,0);
  for (i=0; i<SONS_OF_ELEM(theElement); i++) SET_SON(theElement,i,NULL);

  /* remove interior node, only if not needed in new refinement */
  if ((Nodes[8]!=NULL)&&(n==4))
  {
    newrule = &(Rules[MARK(theElement)]);
    if (newrule->midSon[4]<0)
    {
      DisposeVertex(theGrid,MYVERTEX(Nodes[8]));
      DisposeNode(theGrid,Nodes[8]);
      Nodes[8] = NULL;
      theContext->Interior = NULL;
    }
  }

  /* set pointers in neighbors to NULL */
  for (i=0; i<2*n; i++)
    if (BackPtrs[i]!=NULL)
      *(BackPtrs[i]) = NULL;

  SETTHEFLAG(theElement,1);
}


/****************************************************************************/
/*																			*/
/* Function:  RefineElement                                                                                             */
/*																			*/
/* Purpose:   refine an element in the given context						*/
/*			  (i)	 corner and midnodes are already allocated				*/
/*			  (ii)	 edges between corner and midnodes are ok				*/
/*			  (iii)  create interior nodes and edges						*/
/*			  (iv)	 create sons and set references to sons                                 */
/*																			*/
/* Param:	  GRID *theGrid: grid level of sons of theElement				*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ElementContext *theContext: current context of element		*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*			  INT 1: fatal memory error                                                                     */
/*																			*/
/****************************************************************************/

static INT RefineElement (GRID *theGrid, ELEMENT *theElement, ElementContext *theContext)
{
  INT i,j,n,m,s;
  NODE *theNode;
  ELEMENT *theSon;
  VERTEX *theVertex;
  COORD sx,sy;
  NODE *Nodes[9];                                       /* corner nodes of an element to be ref.*/
  NODE **MidNodes;                                      /* nodes on refined edges				*/
  ELEMENT **Neighbors;                          /* neighbors on next higher level		*/
  ELEMENT ***BackPtrs;                          /* locations where neighbs store my adr */
  ELEMENTSIDE *Sides[8];
  ELEMENTSIDE *oldSide,*newSide0,*newSide1;
  COORD_VECTOR pos,diff,corr;
  COORD *corners[MAX_CORNERS_OF_ELEM];
  INT boundaryelement,found;
  REFRULE *rule;
        #ifdef __version23__
  EDGE *theEdge;
        #endif

  if (MARK(theElement)==0) return(0);

  n = TAG(theElement);
  MidNodes = theContext->MidNodes;
  Neighbors = theContext->Neighbors;
  BackPtrs = theContext->BackPtrs;
  for (i=0; i<9; i++) Nodes[i] = NULL;
  for (i=0; i<n; i++)
  {
    Nodes[i] = theContext->Nodes[i];
    Nodes[i+4] = theContext->MidNodes[i];
  }
  if (n==4) Nodes[8] = theContext->Interior;            /* only possible for quadrilaterals */


  /* check if T4 on the boundary should be switched to T3 */
  /*	I don't think this is a good way...
     DOUBLE x[3],y[3],wmin,wmax;
     if ((n==3)&&(Rules[MARK(theElement)].rule==T_RED)&&(OBJT(theElement)==BEOBJ))
     {
          for (i=0; i<n; i++)
          {
                  x[i] = (DOUBLE) XC(MYVERTEX(MidNodes[i]));
                  y[i] = (DOUBLE) YC(MYVERTEX(MidNodes[i]));
          }
          wmin = 1000.0; wmax = -1000.0;
          QualityElement(n,x,y,&wmin,&wmax);
          if ((wmin<MINANGLE)||(wmax>MAXANGLE))
          {
                  for (i=0; i<n; i++)
                          if (SIDE(theElement,i)!=NULL) break;
                  SETMARK(theElement,T_BISECT_3+i);
          }
     }
   */
  /* get rule */
  rule = &(Rules[MARK(theElement)]);

  /* get interior node (8) */
  if ((Nodes[8]==NULL)&&(rule->midSon[4]>=0))
  {
    /* we need an interior node */
    theVertex = CreateInnerVertex(theGrid,NULL);
    theNode = CreateNode(theGrid,Nodes[0]);
    if ((theNode==NULL)||(theVertex==NULL)) return(1);
    SETUSED(theNode,1);
    theGrid->status |= 1;
    sx = sy = 0.0;
    for (j=0; j<n; j++)
    {
      corners[j] = CVECT(MYVERTEX(Nodes[j]));
      sx += corners[j][_X_];
      sy += corners[j][_Y_];
    }
    /* standard coordinates of the midnode */
    XC(theVertex) = 0.25*sx;
    YC(theVertex) = 0.25*sy;
    XI(theVertex) = 0.0;
    ETA(theVertex) = 0.0;
    /* adjust the midnode for quadrilaterals with boundary sides */
    found = 0;
    V2_CLEAR(corr);
    if (OBJT(theElement)==BEOBJ)
      for (j=0; j<n; j++)
        if (Nodes[j+4]!=NULL)
          if (SIDE(theElement,j)!=NULL)
          {
            /* the midnode lies on a boundary side */
            found++;
            /* calculate difference to sidemid */
            V2_LINCOMB(0.5,corners[j],0.5,corners[(j+1)%n],pos);
            V2_LINCOMB(0.5,CVECT(MYVERTEX(Nodes[j+4])),-0.5,pos,diff);
            V2_ADD(corr,diff,corr);
          }
    if (found)
    {
      /* now shift the midnode */
      V2_ADD(CVECT(theVertex),corr,CVECT(theVertex));

      /* calc local coordinates in the father element */
      if (GlobalToLocal2d(n,(const COORD **) corners,CVECT(theVertex),LCVECT(theVertex))!=0) return (1);
    }
    VFATHER(theVertex) = theElement;
    NFATHER(theNode) = NULL;
    MYVERTEX(theNode) = theVertex;
    TOPNODE(theVertex) = theNode;
    Nodes[8] = theNode;
  }

  /* create interior edges */
  for (s=0; s<rule->nsons; s++)
  {
    m = rule->sons[s].tag;
    for (i=0; i<m; i++)
      if ((!rule->sons[s].outer[i])&&(rule->sons[s].nb[i]>s))
      {
        if (CreateEdge(theGrid,Nodes[rule->sons[s].corners[i]],Nodes[rule->sons[s].corners[(i+1)%m]])==NULL) return(1);
      }
                #ifdef __version23__
    /* introduce extra diagonal edges in quadrilaterals */
    if (m==4)
    {
      theEdge = CreateEdge(theGrid,Nodes[rule->sons[s].corners[0]],Nodes[rule->sons[s].corners[2]]);
      if (theEdge==NULL) return(1);
      SETEXTRA(theEdge,1);
      theEdge = CreateEdge(theGrid,Nodes[rule->sons[s].corners[1]],Nodes[rule->sons[s].corners[3]]);
      if (theEdge==NULL) return(1);
      SETEXTRA(theEdge,1);
    }
                #endif
  }

  /* create element sides */
  for (i=0; i<8; i++) Sides[i] = NULL;
  if (OBJT(theElement)==BEOBJ)
  {
    for (i=0; i<n; i++)
    {
      oldSide = SIDE(theElement,i);
      if (oldSide!=NULL)
      {
        if (MidNodes[i]!=NULL)
        {
          newSide0 = CreateElementSide(theGrid);
          newSide1 = CreateElementSide(theGrid);
          if ((newSide0==NULL)||(newSide1==NULL)) return(1);
          Sides[2*i] = newSide0;
          Sides[2*i+1] = newSide1;
          SEGDESC(newSide0) = SEGDESC(newSide1) = SEGDESC(oldSide);
          PARAM(newSide0,0,0) = PARAM(oldSide,0,0);
          PARAM(newSide1,1,0) = PARAM(oldSide,1,0);
          PARAM(newSide0,1,0) = PARAM(newSide1,0,0) = LAMBDA(VSEG(MYVERTEX(MidNodes[i])),0);
        }
        else
        {
          newSide0 = CreateElementSide(theGrid);
          if (newSide0==NULL) return(1);
          Sides[2*i] = newSide0;
          SEGDESC(newSide0) = SEGDESC(oldSide);
          PARAM(newSide0,0,0) = PARAM(oldSide,0,0);
          PARAM(newSide0,1,0) = PARAM(oldSide,1,0);
        }
      }
    }
  }

  /* create elements */
  for (s=0; s<rule->nsons; s++)
  {
    m = rule->sons[s].tag;
    boundaryelement = 0;
    for (i=0; i<m; i++)
      if ((rule->sons[s].outer[i])&&(Sides[rule->sons[s].nb[i]]!=NULL))
      {
        boundaryelement = 1;
        break;
      }
    if (boundaryelement)
      theSon = CreateBoundaryElement(theGrid,NULL,m);
    else
      theSon = CreateInnerElement(theGrid,NULL,m);
    if (theSon==NULL) return(1);
    SET_SON(theElement,s,theSon);
    SETECLASS(theSon,rule->class);
    /*if (rule->class==copy) SETLEVEL(theSon,LEVEL(theElement));*/
    SETTAG(theSon,m);
    SET_EFATHER(theSon,theElement);
  }
  SETNSONS(theElement,rule->nsons);

  /* connect elements */
  for (s=0; s<rule->nsons; s++)
  {
    m = rule->sons[s].tag;
    theSon = SON(theElement,s);
    for (i=0; i<m; i++)
    {
      SET_CORNER(theSon,i,Nodes[rule->sons[s].corners[i]]);
      if (rule->sons[s].outer[i])
      {
        if (Neighbors[rule->sons[s].nb[i]]!=NULL)
        {
          SET_NBELEM(theSon,i,Neighbors[rule->sons[s].nb[i]]);
          *(BackPtrs[rule->sons[s].nb[i]]) = theSon;
        }
        if (OBJT(theSon)==BEOBJ)
          SET_SIDE(theSon,i,Sides[rule->sons[s].nb[i]]);
      }
      else
        SET_NBELEM(theSon,i,SON(theElement,rule->sons[s].nb[i]));
    }
  }

  /* set modified flags of nodes */
  for (i=0; i<9; i++)
    if (Nodes[i]!=NULL)
      SETMODIFIED(Nodes[i],1);

  /* mark element as changed and set modified flag in grid record */
  SETTHEFLAG(theElement,1);
  SETMODIFIED(theGrid,1);

  return(0);
}




/****************************************************************************/
/*																			*/
/* Function:  RefineGrid													*/
/*																			*/
/* Purpose:   refine one level of the grid									*/
/*																			*/
/* Param:	  GRID *theGrid: grid level to refine							*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*			  INT 1: fatal memory error                                                                     */
/*																			*/
/****************************************************************************/

static INT exitFlag=0;

static INT RefineGrid (GRID *theGrid)
{
  ELEMENT *theElement;
  ElementContext theContext;
  GRID *fineGrid;
  NODE *theNode;

  fineGrid = theGrid->finer;
  if (fineGrid==NULL) return(1);

  /* refine elements */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (REFINE(theElement)!=MARK(theElement))
    {
      GetCurrentContext(theElement,&theContext);
      UnrefineElement(fineGrid,theElement,&theContext);
      if (UpdateContext(fineGrid,theElement,&theContext)!=0) return(1);
      if (RefineElement(fineGrid,theElement,&theContext)!=0) return(1);
      SETREFINE(theElement,MARK(theElement));
    }
  }

  /* reset MARK flags */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    SETMARK(theElement,0);

  /* set node class on next level */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
    if (SONNODE(theNode)!=NULL)
    {
      SETCLASS(SONNODE(theNode),NCLASS(theNode));
      if (NCLASS(theNode)>=2) TOPNODE(MYVERTEX(theNode)) = SONNODE(theNode);
    }

  return(0);
}


/****************************************************************************/
/*D
   RefineMultiGrid - Refine whole multigrid structure

   SYNOPSIS:
   INT RefineMultiGrid (MULTIGRID *theMG, INT flag);

   PARAMETERS:
   .  theMG - multigrid to refine
   .  flag - refinement mode, see below
   .  direction - element evaluation direction

   DESCRIPTION:
   This function refines whole multigrid structure. It reads the
   refinement tags set with the 'MarkForRefinement' function and
   produces a consistent triangulation where at least all the
   selected elements are refined.

   The flag parameter can have two values. 'GM_REFINE_TRULY_LOCAL' refines
   in such a way that a grid level is not required to cover the whole
   domain. In this case the memory and run-time requirements are proportional
   to the number of elements. If 'flag' is set to 'GM_COPY_ALL' then
   elements not selected for refinement are copied to the next finer level,
   i.e. each grid level covers the whole domain but memory and run-time
   requirements are only optimal if there is a geometric growth in the
   number of refined elements.

   The direction is not considered.

   RETURN VALUE:
   INT
   .n    GM_OK when ok
   .n    GM_ERROR when out of memory, but data structure as before
   .n    GM_FATAL when fatal memory error, data structure corrupted
   D*/
/****************************************************************************/

INT RefineMultiGrid (MULTIGRID *theMG, INT flag, EVECTOR *direction)
{
  INT j,k,r;
  INT newlevel;
  NODE *theNode;
  GRID *theGrid,*theFineGrid;
  ELEMENT *theElement;

  rFlag=flag;           /* set global variable */

  j = theMG->topLevel;

  for (k=0; k<=j; k++)
    for (theElement=theMG->grids[k]->elements; theElement!=NULL; theElement=SUCCE(theElement))
      SETTHEFLAG(theElement,0);

#ifdef __version3__
  PrepareAlgebraModification(theMG);
#endif

  /* compute modification of coarser levels from above */
  for (k=j; k>0; k--)
  {
    theGrid = theMG->grids[k];

    /* make refinement rules consistent */
    if (CloseGrid(theGrid)<0)
    {
      PrintErrorMessage('E',"RefineMultiGrid","error in occured in CloseGrid");
      return(GM_ERROR);
    }

    /* restrict marks on irregular elements and handle COARSEN flags */
    RestrictMarks(theMG->grids[k-1]);
  }

  newlevel = 0;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    if (k<j) theFineGrid = GRID_ON_LEVEL(theMG,k+1);else theFineGrid = NULL;

    /* reset some flags */
    SETMODIFIED(theGrid,0);
    for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode)) SETMODIFIED(theNode,0);

    /* leave only regular marks */
    for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    {
      if ((ECLASS(theElement)==REGULAR_CLASS)&&(Rules[MARK(theElement)].class==REGULAR_CLASS)) continue;
      SETMARK(theElement,0);
    }

    /* determine regular and irregular elements on next level */
    if ((r = CloseGrid(theGrid))<0)
    {
      PrintErrorMessage('E',"RefineMultiGrid","error occured in CloseGrid");
      return(GM_FATAL);
    }

    /* determine copies on next level		*/
    /* this involves algebra in version 3 ! */
#ifdef __version23__
    ComputeCopies(theGrid);
#endif
#ifdef __version3__
    Version3_ComputeCopies(theGrid);
#endif

#ifdef __version3__
    /* dispose connections that may be changed on next level, this is determined */
    /* by the neighborhood of elements were MARK != REFINE.                                      */
    /* This will leave some flags where to rebuild connections later			 */
    if (k<j)
    {
      for (theElement=FIRSTELEMENT(theFineGrid); theElement!=NULL; theElement=SUCCE(theElement))
      {
        assert (EFATHER(theElement) != NULL);
        if (REFINE(EFATHER(theElement))!=MARK(EFATHER(theElement)))
          if (DisposeConnectionsInNeighborhood(theFineGrid,theElement)!=GM_OK)
            return(GM_FATAL);
      }
    }
#endif

    /* create a new grid level, if at least one element is refined on finest level */
    if ((r>0)&&(k==j))
    {
      newlevel = 1;
      if (CreateNewLevel(theMG)==NULL)
      {
        PrintErrorMessage('E',"RefineMultiGrid","could not create new level");
        return(GM_FATAL);
      }
      theFineGrid = GRID_ON_LEVEL(theMG,j+1);
    }

    /* now really manipulate the next finer level */
    if ((k<j)||(newlevel))
      if (RefineGrid(theGrid)!=0) return(GM_FATAL);

#ifdef __version3__
    if ((k<j)||(newlevel))
    {
      /* now rebuild connections in neighborhood of elements which have EBUILDCON set */
      /* This flag has been set either by GridDisposeConnection or by CreateElement	*/
      if (GridCreateConnection(theFineGrid)) return (GM_FATAL);
      /* and compute the vector classes on the new (or changed) level */

      if (rFlag==GM_COPY_ALL)
        for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
          SeedVectorClasses(theElement);
      else
      {
        ClearVectorClasses(theFineGrid);
        for (theElement=theFineGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
          if (ECLASS(theElement)>=IRREGULAR_CLASS) SeedVectorClasses(theElement);
        PropagateVectorClasses(theFineGrid);
      }
    }
#endif
  }

  DisposeTopLevel(theMG);       /* is only done when highest level is empty !*/

  return(GM_OK);
}


/****************************************************************************/
/*D
   MarkForRefinement - Mark an element for refinement

   SYNOPSIS:
   INT MarkForRefinement (ELEMENT *theElement, INT rule, INT side);

   PARAMETERS:
   .  theElement - element to be tagged for refinement
   .  rule - refinement rule
   .  side - side, if rule is not invariant under rotation

   DESCRIPTION:
   This function marks an element for refinement. Only elements where
   'EstimateHere' is true may be marked for refinement. Rule numbers
   are declared in the 'REFINEMENT' page.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 if wrong rule.
   D*/
/****************************************************************************/

INT MarkForRefinement (ELEMENT *theElement, INT rule, INT side)
{
  INT n;

  if (!EstimateHere(theElement))
    return (GM_ERROR);

  n = TAG(theElement);
  SETCOARSEN(theElement,0);

  if (n==TRIANGLE)
  {
    switch (rule)
    {
    case UNREFINE :
      SETCOARSEN(theElement,1);
      SETMARK(theElement,BasicRules[NO_REF].offset);
      break;

    case NO_REFINEMENT :
      SETMARK(theElement,BasicRules[NO_REF].offset);
      break;

    case COPY :
      SETMARK(theElement,BasicRules[T_COPY].classOffset[REGULAR_CLASS]);
      break;

    case RED :
      SETMARK(theElement,BasicRules[T_RED].classOffset[REGULAR_CLASS]);
      break;

    case BISECTION_3 :
      if (side<0) return (GM_ERROR);
      SETMARK(theElement,BasicRules[T_BISECT_3].classOffset[REGULAR_CLASS]+BasicRules[T_BISECT_3].variantOffset[side%3]);
      break;

    case BISECTION_1 :
      if (side<0) return (GM_ERROR);
      SETMARK(theElement,BasicRules[T_BISECT_1].classOffset[REGULAR_CLASS]+BasicRules[T_BISECT_1].variantOffset[side%3]);
      break;

    case BISECTION_2_Q :
      if (side<0) return (GM_ERROR);
      SETMARK(theElement,BasicRules[T_BISECT_2_Q].classOffset[REGULAR_CLASS]+BasicRules[T_BISECT_2_Q].variantOffset[side%3]);
      break;

    case BISECTION_2_T1 :
      if (side<0) return (GM_ERROR);
      SETMARK(theElement,BasicRules[T_BISECT_2_T1].classOffset[REGULAR_CLASS]+BasicRules[T_BISECT_2_T1].variantOffset[side%3]);
      break;

    case BISECTION_2_T2 :
      if (side<0) return (GM_ERROR);
      SETMARK(theElement,BasicRules[T_BISECT_2_T2].classOffset[REGULAR_CLASS]+BasicRules[T_BISECT_2_T2].variantOffset[side%3]);
      break;

    default :
      return(GM_ERROR);
    }
  }
  else
  {
    switch (rule)
    {
    case UNREFINE :
      SETCOARSEN(theElement,1);
      SETMARK(theElement,BasicRules[NO_REF].offset);
      break;

    case NO_REFINEMENT :
      SETMARK(theElement,BasicRules[NO_REF].offset);
      break;

    case COPY :
      SETMARK(theElement,BasicRules[Q_COPY].classOffset[REGULAR_CLASS]);
      break;

    case RED :
      SETMARK(theElement,BasicRules[Q_RED].classOffset[REGULAR_CLASS]);
      break;

    case BLUE :
      if (side<0) return (GM_ERROR);
      SETMARK(theElement,BasicRules[Q_BLUE].classOffset[REGULAR_CLASS]+BasicRules[Q_BLUE].variantOffset[side%2]);
      break;

    default :
      return(GM_ERROR);
    }
  }

  return(GM_OK);
}

/****************************************************************************/
/*D
   EstimateHere	- Determines elements allowed for refinement

   SYNOPSIS:
   INT EstimateHere (ELEMENT *theElement);

   PARAMETERS:
   .  theElement - element to be checked

   DESCRIPTION:
   This function returns true (1) when element can be tagged for refinement.

   RETURN VALUE:
   INT
   .n      0, do not tag element.
   .n      1, element can be tagged for refinement
   D*/
/****************************************************************************/

INT EstimateHere (ELEMENT *theElement)
{
  return(REFINE(theElement)==NO_REF);
}

/****************************************************************************/
/*D
   GetRefinementMark - Get rule and variant of refinement

   SYNOPSIS:
   INT GetRefinementMark (const ELEMENT *theElement, INT *rule, INT *side);

   PARAMETERS:
   .  theElement - element to be questioned
   .  rule - to be filled with current refinement rule
   .  side - to be filled with side, if rule is not invariant under rotation

   DESCRIPTION:
   This function retieves rule and side of refinement.

   RETURN VALUE:
   INT
   .n    GM_RULE_WITH_ORIENTATION if side information valid
   .n    GM_RULE_WITHOUT_ORIENTATION if rule without orientation.
   D*/
/****************************************************************************/

INT GetRefinementMark (const ELEMENT *theElement, INT *rule, INT *side)
{
  switch (Rules[MARK(theElement)].rule)
  {
  case NO_REF :
    *rule=NO_REFINEMENT;
    if (COARSEN(theElement)) *rule = UNREFINE;
    break;
  case T_COPY :                    *rule=COPY; break;
  case T_RED :                     *rule=RED; break;
  case T_BISECT_1 :                *rule=BISECTION_1; break;
  case T_BISECT_2_T1 :             *rule=BISECTION_2_T1;  break;
  case T_BISECT_2_T2 :     *rule=BISECTION_2_T2;  break;
  case T_BISECT_2_Q :              *rule=BISECTION_2_Q; break;
  case T_BISECT_3 :                *rule=BISECTION_3; break;

  case Q_COPY :                    *rule=COPY; break;
  case Q_RED :                     *rule=RED; break;
  case Q_CLOSE_1 :                 *rule=NO_REFINEMENT;  break;
  case Q_BLUE :                    *rule=BLUE; break;
  case Q_CLOSE_2 :                 *rule=NO_REFINEMENT;   break;
  case Q_CLOSE_3 :                 *rule=NO_REFINEMENT;   break;
  default :                                *rule=NO_REFINEMENT;  break;
  }
  if(Rules[MARK(theElement)].class == REGULAR_CLASS)
  {
    *side=Rules[MARK(theElement)].variant;
    return(GM_RULE_WITH_ORIENTATION);
  }
  *side=0;
  return(GM_RULE_WITHOUT_ORIENTATION);
}

/****************************************************************************/
/*D
   GetSons - Get pointers to all sons of an element

   SYNOPSIS:
   INT GetSons (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS]);

   PARAMETERS:
   .  theElement - pointer to an element
   .  SonList - array of pointers to elements where result is placed

   DESCRIPTION:
   This function fills a user supplied array of pointers with all
   sons of a given element. Use this function for access to sons,
   since it is available in 2D and 3D. In 3D the sons are not directly
   accessible!

   RETURN VALUE:
   INT
   .n     0 when ok
   .n     1 when error occured.
   D*/
/****************************************************************************/

INT GetSons (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS])
{
  INT i;

  for (i=0; i<NSONS(theElement); i++) SonList[i] = SON(theElement,i);
  return(GM_OK);
}

/****************************************************************************/
/*D
   GetFineNodeOnEdge - Get node in fine grid generated on edge

   SYNOPSIS:
   NODE *GetFineNodeOnEdge (const ELEMENT *theElement, INT side);

   PARAMETERS:
   .  theElement - pointer to an element
   .  side - side number of element

   DESCRIPTION:
   This function returns the node in the fine grid that is generated on the
   edge given by 'side'.

   `This function is only available in 2D!`

   RETURN VALUE:
   NODE*
   .n      pointer to node
   .n      NULL if none.
   D*/
/****************************************************************************/

NODE *GetFineNodeOnEdge (const ELEMENT *theElement, INT side)
{
  INT ref;

  ref = REFINE(theElement);

  if (!Rules[ref].pattern[side])
    return (NULL);

  return (CORNER(SON(theElement,Rules[ref].midSon[side]),Rules[ref].midCorner[side]));
}

/****************************************************************************/
/*
   InitUgRefine	- Initialize rules from basic rules

   SYNOPSIS:
   INT InitRefine2d (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initialize rules from basic rules.

   RETURN VALUE:
   void
 */
/****************************************************************************/

INT InitRefine2d (void)
{
  REFRULE *rule;
  basicRule *bRule;
  SON_DATA *bSond,*sond;
  INT i,j,k,l,m,n,s,BRule,Rule,Side;
  INT class,followRule,followVariant,cPattern,cPattern0,nrVar;

  Rule=0;

  for (BRule=0; BRule<MAXBASICRULES; BRule++)
  {
    bRule=&(BasicRules[BRule]);

    /* compute offsets */
    bRule->offset=Rule;

    nrVar=0;
    for (k=0; k<4; k++)
    {
      if (bRule->variants[k])
      {
        if (k>nrVar)
        {
          PrintErrorMessage('E',"InitGM","variants are not consecutively ordered");
          return(__LINE__);
        }

        bRule->variantOffset[k]=k;
        nrVar++;
      }
      else
        bRule->variantOffset[k]=k%nrVar;
    }

    /* go through all classes for which bRule exists */
    for (class=0; class<3; class++)
    {
      bRule->classOffset[class]=bRule->offset;                                  /* smallest possible damage if not existent */
      if (bRule->classes[class])
      {
        bRule->classOffset[class]=Rule;
        /* go through all variants that exist */
        for (k=0; k<4; k++)
          if (bRule->variants[k])
          {
            /* generate corresponding rule */
            if (Rule==MAXRULES)
            {
              PrintErrorMessage('E',"InitGM","too many rules to be generated");
              return(__LINE__);
            }
            else
              rule=&(Rules[Rule++]);

            n=rule->tag=bRule->tag;
            rule->rule=BRule;
            rule->class=class;
            rule->variant=k;
            rule->nsons=bRule->nsons;

            if (n==TRIANGLE)
              rule->pattern[3]=0;

            cPattern=0;
            for (Side=n-1; Side>=0; Side--)
              cPattern=2*cPattern+(rule->pattern[Side]=bRule->pattern[(Side-k+n)%n]);

            rule->condensedPattern=cPattern;

            /* now loop through the sons */
            for (l=0; l<4; l++)
            {
              bSond=&(bRule->sons[l]);
              sond=&(rule->sons[l]);
              if (l<rule->nsons)
              {
                sond->tag=bSond->tag;

                if (sond->tag==3)
                  sond->corners[3]=sond->outer[3]=sond->nb[3]=0;
                else
                if (sond->tag!=4)
                {
                  PrintErrorMessage('E',"InitGM","false tag for son");
                  return(__LINE__);
                }

                for (Side=0; Side<sond->tag; Side++)
                {
                  if (bSond->corners[Side]==8)
                    sond->corners[Side]=8;
                  else
                    sond->corners[Side]=(bSond->corners[Side] & 4) + ((bSond->corners[Side] & 3) + k) %n;

                  if ((sond->outer[Side]=bSond->outer[Side])==0)
                  {
                    /* inner side for refinement */
                    sond->nb[Side]=bSond->nb[Side];
                  }
                  else
                  {
                    /* outer side for refinement */
                    sond->nb[Side]=(bSond->nb[Side]+2*k)%(2*n);
                  }
                }
              }
              else
              {
                if (bSond->tag!=0)
                {
                  PrintErrorMessage('E',"InitGM","inconsistency with nsons");
                  return(__LINE__);
                }
                sond->tag=0;
                for (Side=0; Side<4; Side++)
                  sond->corners[Side]=sond->outer[Side]=sond->nb[Side]=0;
              }
            }

            /* compute data for access of the centerpoints of the outer edges */
            for (Side=0; Side<5; Side++)
            {
              m = Side+4;
              rule->midSon[Side] = -1;
              rule->midCorner[Side] = -1;
              for (s=0; s<rule->nsons; s++)
              {
                for (i=0; i<rule->sons[s].tag; i++)
                  if (rule->sons[s].corners[i]==m)
                  {
                    rule->midSon[Side] = s;
                    rule->midCorner[Side] = i;
                    break;
                  }
                if (rule->midSon[Side]>=0) break;
              }
            }

            /* compute data for access to the neighbors on the next level */
            for (Side=0; Side<8; Side++)
            {
              rule->nbSon[Side] = -1;
              rule->nbSide[Side] = -1;
              for (s=0; s<rule->nsons; s++)
              {
                for (j=0; j<rule->sons[s].tag; j++)
                  if ((rule->sons[s].outer[j])&&(rule->sons[s].nb[j]==Side))
                  {
                    rule->nbSon[Side]  = s;
                    rule->nbSide[Side] = j;
                    break;
                  }
                if (rule->nbSon[Side]>=0) break;
              }
            }

          }
      }
    }
  }

  /* set size of rules field */
  MaxRules=Rule;

  /* generate the follow rules */
  for (Rule=0; Rule<MaxRules; Rule++)
  {
    rule=&(Rules[Rule]);
    BRule=rule->rule;
    bRule=&(BasicRules[BRule]);
    n=rule->tag;
    class=rule->class;
    k=rule->variant;

    if (n==TRIANGLE)
    {
      for (cPattern=0; cPattern<8; cPattern++)
      {
        /* transform to 0-orientation */
        cPattern0=(cPattern>>k) + ((cPattern<<(n-k))&0x00000007);

        if (bRule->followRule[cPattern0]==-1)
        {
          followRule=StandardFollowRule_T[cPattern0];
          followVariant=StandardFollowVariant_T[cPattern0];
        }
        else
        {
          followRule=bRule->followRule[cPattern0];
          followVariant=bRule->followVariant[cPattern0];
        }

        /* transform back to Rules */
        rule->followRule[cPattern]=BasicRules[followRule].classOffset[class]+BasicRules[followRule].variantOffset[(followVariant+k)%n];
      }

      for (cPattern=8; cPattern<16; cPattern++)
        rule->followRule[cPattern]=0;                           /* should never be needed */
    }
    else
    {
      for (cPattern=0; cPattern<16; cPattern++)
      {
        /* transform to 0-orientation */
        cPattern0=(cPattern>>k) + ((cPattern<<(n-k))&0x0000000f);

        if (bRule->followRule[cPattern0]==-1)
        {
          followRule=StandardFollowRule_Q[cPattern0];
          followVariant=StandardFollowVariant_Q[cPattern0];
        }
        else
        {
          followRule=bRule->followRule[cPattern0];
          followVariant=bRule->followVariant[cPattern0];
        }

        /* transform back to Rules */
        rule->followRule[cPattern]=BasicRules[followRule].classOffset[class]+BasicRules[followRule].variantOffset[(followVariant+k)%n];
      }
    }
  }

  return(GM_OK);
}
