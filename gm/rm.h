// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  rule manager header file										*/
/*																			*/
/* Purpose:   defines data structures for refinement rules					*/
/*																			*/
/* Author:	  Stefan Lang                                                                           */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*																			*/
/* History:   20.11.95 begin, ugp version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __RULEMANAGER__
#define __RULEMANAGER__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __GM__
#include "gm.h"
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

/* defines for edge types */
#define INNER_EDGE          1
#define SIDE_EDGE           2
#define HALF_FATHER_EDGE    3
#define FATHER_EDGE         4

#define NEXTSIDEMASKHEX         0x00000007
#define NEXTSIDEHEX(i,n)        (((i) & (NEXTSIDEMASKHEX<<(3*(n))))>>(3*(n)))

#define FATHER_SIDE_OFFSET 20 /* greater values indicate outside faces */

#define MAX_NEW_CORNERS(tag) MaxNewCorners[tag]       /* midpoints on edges and sides plus center           */
#define MAX_NEW_EDGES(tag)   MaxNewEdges[tag]         /* maximal new edges of type 1/2                      */
#define MAX_RULES(tag)       MaxRules[tag]            /* number of rules for one element type               */
#define CENTER_NODE_INDEX(e) CenterNodeIndex[TAG(e)] /* index of centernode in sonandnode and elemcontext */

#define MARK2RULE(e,m)          (m)
#define MARK2RULEADR(e,m)       (&(RefRules[TAG(e)][m]))
#define RULE2PATTERN(e,r)       (RefRules[TAG(e)][r].pattern)
#define RULE2PAT(e,r)           (RefRules[TAG(e)][r].pat)
#define MARK2PAT(e,r)           (RefRules[TAG(e)][r].pat)
#define MARK2PATTERN(e,m)       (RefRules[TAG(e)][m].pattern)

/* TODO: delete this
 #define PATTERNOFELEM(e)	((dim==2)?(EDGEPATTERN(e)):((SIDEPATTERN(e)<<6) | EDGEPATTERN(e))) */
/* #define PATTERN2RULE(e,p)	(Pattern2Rule[TAG(e)][(p)]) */
#define PATTERN2RULE(e,p)       (Patterns2Rules((e),(p)))
#define RULE2MARK(e,r)          (RefRules[TAG(e)][r].mark)
/* TODO: delete this */
/* #define PATTERN2MARK(e,p)	(((Patterns2Rules((e),(p)))>=0)?(RefRules[TAG(e)][Pattern2Rule[TAG(e)][p]].mark):-1) */
#define PATTERN2MARK(e,p)       (((PATTERN2RULE(e,p))>=0) ? (RefRules[TAG(e)][PATTERN2RULE(e,p)].mark) : -1)

/* dimension dependent MAX_CORNERS_OF_ELEM */
#define MAX_CORNERS_OF_ELEM_2D  4
#define MAX_CORNERS_OF_ELEM_3D  8
#define MAX_CORNERS_OF_ELEM_DIM CONCAT3(MAX_CORNERS_OF_ELEM_,DIM,D)

/* dimension dependent MAX_EDGES_OF_ELEM */
#define MAX_EDGES_OF_ELEM_2D    4
#define MAX_EDGES_OF_ELEM_3D    12
#define MAX_EDGES_OF_ELEM_DIM CONCAT3(MAX_EDGES_OF_ELEM_,DIM,D)

/* dimension dependent MAX_NEW_CORNERS */
#define MAX_NEW_CORNERS_2D    5
#define MAX_NEW_CORNERS_3D   19
#define MAX_NEW_CORNERS_DIM CONCAT3(MAX_NEW_CORNERS_,DIM,D)

/* dimension dependent MAX_NEW_EDGES */
#define MAX_NEW_EDGES_2D   12
#define MAX_NEW_EDGES_3D   54
#define MAX_NEW_EDGES_DIM CONCAT3(MAX_NEW_EDGES_,DIM,D)

/* dimension dependent MAX_SIDES_OF_ELEM */
#define MAX_SIDES_OF_ELEM_2D    4
#define MAX_SIDES_OF_ELEM_3D    6
#define MAX_SIDES_OF_ELEM_DIM CONCAT3(MAX_SIDES_OF_ELEM_,DIM,D)

/* dimension dependent MAX_SONS */
#define MAX_SONS_2D    4
#define MAX_SONS_3D    12
#define MAX_SONS_DIM CONCAT3(MAX_SONS_,DIM,D)

#define LEAFELEM(e)             (NSONS(e)==0)

/* indices of rules in rule array */
#define NO_REF          0

#define T_NOREF                         0
#define T_COPY                  1
#define T_RED                           2
#define T_BISECT_1_0            3
#define T_BISECT_1_1            4
#define T_BISECT_1_2            5
#define T_BISECT_2_T1_0         7
#define T_BISECT_2_T1_1         8
#define T_BISECT_2_T1_2         6
#define T_BISECT_2_T2_0         11
#define T_BISECT_2_T2_1         9
#define T_BISECT_2_T2_2         10
#define T_BISECT_2_Q_0          12
#define T_BISECT_2_Q_1          13
#define T_BISECT_2_Q_2          14
#define T_BISECT_3_0            15
#define T_BISECT_3_1            16
#define T_BISECT_3_2            17

#define Q_NOREF                         0
#define Q_COPY                          1
#define Q_RED                           2
#define Q_CLOSE_1_0                     3
#define Q_CLOSE_1_1                     4
#define Q_CLOSE_1_2                     5
#define Q_CLOSE_1_3                     6
#define Q_BLUE_0                        7
#define Q_BLUE_1                        8
#define Q_CLOSE_2_0                     9
#define Q_CLOSE_2_1                     10
#define Q_CLOSE_2_2                     11
#define Q_CLOSE_2_3                     12
#define Q_CLOSE_3_0                     13
#define Q_CLOSE_3_1                     14
#define Q_CLOSE_3_2                     15
#define Q_CLOSE_3_3                     16

#define TET_COPY                        1
#define FULL_REFRULE        Pattern2Rule[TETRAHEDRON][0x3F]
#define FULL_REFRULE_0_5    (Pattern2Rule[TETRAHEDRON][0x3F]+1)
#define FULL_REFRULE_1_3    (Pattern2Rule[TETRAHEDRON][0x3F]+2)
#define FULL_REFRULE_2_4    (Pattern2Rule[TETRAHEDRON][0x3F]+0)

#define PYR_COPY                        1
#define PYR_RED                         2

#define HEXA_COPY                       1
#define HEXA_RED                        2

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef INT (*FULLREFRULEPTR)(ELEMENT *);

typedef struct {
  INT type;                                    /* one of the variable types above          */
  INT locked;                                  /* may not be changed or deleted            */
  union envitem *next;
  union envitem *previous;                     /* double linked list of environment items  */
  char name[NAMESIZE];                         /* name of that item (view)                 */

  /* full ref rule spezific stuff */
  FULLREFRULEPTR theFullRefRule;               /* the best full refrule                    */

} FULLREFRULE;

/*
   edge type  : 1 = inner edge of the father
                2 = inner edge of one side (3D only)
                3 = half an edge of the father element
                4 = edge of the father itself
 */

struct edgedata {
  SHORT type;                      /* interior edge of the tetra or only of one side        */
  SHORT from;                      /* indices of the first endpoint (0..9 for t)            */
  SHORT to;                        /* indices of the second endpoint (4..9 for t)           */
  SHORT side;                      /* side, for which edge is interior, if type 2 (3D only) */
};

struct sondata {
  SHORT tag;                                     /* which element type is the son     */
  SHORT corners[MAX_CORNERS_OF_ELEM_DIM];        /* corners of the son                */
  SHORT nb[MAX_SIDES_OF_ELEM_DIM];               /* neighbors of this son             */
  /* < 20 if neighbor has same father  */
  /* >= 20 if neighbor has other father*/
  INT path;                                      /* path used in GetSons() for tetras */
};

struct refrule {
  SHORT tag;                                                        /* which element type this rule can refine */
  SHORT mark;                                                       /* the mark of this rule                   */
  SHORT class;                                                      /* class of rule:3bits for COPY, IREG, REG */
  SHORT nsons;                                              /* number of sons rule creates             */
  SHORT pattern[MAX_NEW_CORNERS_DIM];                       /* stores which edges are refined          */
  SHORT pat;                                                /* bitwise format of pattern               */
  SHORT sonandnode[MAX_NEW_CORNERS_DIM][2];                 /* for each new node the number of the son */
  /* and the local node number of the node   */
  struct edgedata edges[MAX_NEW_EDGES_DIM];
  struct sondata sons[MAX_SONS_DIM];
};

typedef struct sondata SONDATA;
typedef struct edgedata EDGEDATA;
typedef struct refrule REFRULE;


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

extern INT MaxRules[TAGS];
extern INT MaxNewCorners[TAGS];
extern INT MaxNewEdges[TAGS];
extern INT CenterNodeIndex[TAGS];
extern REFRULE                  *RefRules[TAGS];
extern SHORT                    *Pattern2Rule[TAGS];
extern FULLREFRULEPTR theFullRefRule;


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT ShowRefRule (INT tag, INT nb);

INT InitRuleManager (void);
INT Patterns2Rules(ELEMENT *theElement,INT pattern);

#ifdef __THREEDIM__
INT SetAlignementPtr (MULTIGRID *theMG, EVECTOR *direction);
#endif

#endif
