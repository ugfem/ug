// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file rm.h
 * \ingroup gm
 */


/****************************************************************************/
/*                                                                                                                                                      */
/* File:          rule manager header file                                                                              */
/*                                                                                                                                                      */
/* Purpose:   defines data structures for refinement rules                                      */
/*                                                                                                                                                      */
/* Author:        Stefan Lang                                                                           */
/*                        Institut fuer Computeranwendungen III                                                 */
/*                        Universitaet Stuttgart                                                                                */
/*                        Pfaffenwaldring 27                                                                                    */
/*                        70550 Stuttgart                                                                                               */
/*                                                                                                                                                      */
/* History:   20.11.95 begin, ugp version 3.0                                                           */
/*                                                                                                                                                      */
/* Remarks:                                                                                                                             */
/*                                                                                                                                                      */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                                                                                                      */
/* auto include mechanism and other include files                                                       */
/*                                                                                                                                                      */
/****************************************************************************/

#ifndef __RULEMANAGER__
#define __RULEMANAGER__

#include "gm.h"
#include "refine.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                                                                                                      */
/* defines in the following order                                                                                       */
/*                                                                                                                                                      */
/*                compile time constants defining static data size (i.e. arrays)        */
/*                other constants                                                                                                       */
/*                macros                                                                                                                        */
/*                                                                                                                                                      */
/****************************************************************************/

/* Declaration of TET_RULESET has been moved to config.h */
/* uncomment this if you want to use the full rule set for tetrahedra */

/** \brief Defines for edge types */
enum {INNER_EDGE        = 1,
      SIDE_EDGE         = 2,
      HALF_FATHER_EDGE  = 3,
      FATHER_EDGE       = 4};

#define NEXTSIDEMASKHEX         0x00000007
#define NEXTSIDEHEX(i,n)        (((i) & (NEXTSIDEMASKHEX<<(3*(n))))>>(3*(n)))

#define FATHER_SIDE_OFFSET 100 /* greater values indicate outside faces */

#define MAX_NEW_CORNERS(tag) MaxNewCorners[(tag)]       /* midpoints on edges and sides plus center           */
#define MAX_NEW_EDGES(tag)   MaxNewEdges[(tag)]         /* maximal new edges of type 1/2                      */
#define MAX_RULES(tag)       MaxRules[(tag)]            /* number of rules for one element type               */
#define CENTER_NODE_INDEX(e) CenterNodeIndex[TAG(e)]  /* index of centernode in sonandnode and elemcontext */
#define CENTER_NODE_INDEX_TAG(t) CenterNodeIndex[(t)]     /* index of centernode in sonandnode and elemcontext */

#define MARK2RULE(e,m)          (m)
#define MARK2RULEADR(e,m)       (&(RefRules[TAG(e)][(m)]))
#define RULE2PATTERN(e,r)       (RefRules[TAG(e)][(r)].pattern)
#define RULE2PAT(e,r)           (RefRules[TAG(e)][(r)].pat)
#define MARK2PAT(e,r)           (RefRules[TAG(e)][(r)].pat)
#define MARK2PATTERN(e,m)       (RefRules[TAG(e)][(m)].pattern)

#define PATTERN2RULE(e,p)       (Patterns2Rules((e),(p)))
#define RULE2MARK(e,r)          (RefRules[TAG(e)][(r)].mark)
#define PATTERN2MARK(e,p)       (((PATTERN2RULE((e),(p)))>=0) ? (RefRules[TAG(e)][PATTERN2RULE((e),(p))].mark) : -1)

#ifdef __SR2201__
#define NODE_OF_RULE(e,m,i) ((*MARK2RULEADR((e),(m))).sonandnode[(i)][0]!=-1)
#else
#define NODE_OF_RULE(e,m,i)     (MARK2RULEADR((e),(m))->sonandnode[(i)][0]!=-1)
#endif

#define CONCAT(a,b,c)            CONCAT_AUX(a,b,c)
#define CONCAT_AUX(a,b,c)        a ## b ## c

/* dimension dependent MAX_CORNERS_OF_ELEM */
#define MAX_CORNERS_OF_ELEM_2D  4
#define MAX_CORNERS_OF_ELEM_3D  8
#define MAX_CORNERS_OF_ELEM_DIM CONCAT(MAX_CORNERS_OF_ELEM_,DIM,D)

/* dimension dependent MAX_EDGES_OF_ELEM */
#define MAX_EDGES_OF_ELEM_2D    4
#define MAX_EDGES_OF_ELEM_3D    12
#define MAX_EDGES_OF_ELEM_DIM CONCAT(MAX_EDGES_OF_ELEM_,DIM,D)

/* dimension dependent MAX_NEW_CORNERS */
#define MAX_NEW_CORNERS_2D    5
#define MAX_NEW_CORNERS_3D   19
#define MAX_NEW_CORNERS_DIM CONCAT(MAX_NEW_CORNERS_,DIM,D)

/* dimension dependent MAX_NEW_EDGES */
#define MAX_NEW_EDGES_2D   12
#define MAX_NEW_EDGES_3D   54
#define MAX_NEW_EDGES_DIM CONCAT(MAX_NEW_EDGES_,DIM,D)

/* dimension dependent MAX_SIDES_OF_ELEM */
#define MAX_SIDES_OF_ELEM_2D    4
#define MAX_SIDES_OF_ELEM_3D    6
#define MAX_SIDES_OF_ELEM_DIM CONCAT(MAX_SIDES_OF_ELEM_,DIM,D)

/* dimension dependent MAX_SONS */
#define MAX_SONS_2D    4
#define MAX_SONS_3D    12
#define MAX_SONS_DIM CONCAT(MAX_SONS_,DIM,D)

#if FATHER_SIDE_OFFSET<MAX_SONS
#error  /* increase FATHER_SIDE_OFFSET correspondigly */
#endif

/* max fine grid nodes of an element */
#define MAX_REFINED_CORNERS_DIM (MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM)

#define IS_REFINED(e)           (REFINE(e)!=NO_REFINEMENT)
#define MARKED(e)               (MARK(e)!=NO_REFINEMENT)
#define LEAFELEM(e)                     (!IS_REFINED(e))

/** \brief Indices of rules in rule array */
enum {T_NOREF               = 0,
      T_COPY                = 1,
      T_RED                 = 2,
      T_BISECT_1_0          = 3,
      T_BISECT_1_1          = 4,
      T_BISECT_1_2          = 5,
      T_BISECT_2_T1_0       = 7,
      T_BISECT_2_T1_1       = 8,
      T_BISECT_2_T1_2       = 6,
      T_BISECT_2_T2_0       = 11,
      T_BISECT_2_T2_1       = 9,
      T_BISECT_2_T2_2       = 10,
      T_BISECT_2_Q_0        = 12,
      T_BISECT_2_Q_1        = 13,
      T_BISECT_2_Q_2        = 14,
      T_BISECT_3_0          = 15,
      T_BISECT_3_1          = 16,
      T_BISECT_3_2          = 17};

enum {Q_NOREF,
      Q_COPY,
      Q_RED,
      Q_CLOSE_1_0,
      Q_CLOSE_1_1,
      Q_CLOSE_1_2,
      Q_CLOSE_1_3,
      Q_BLUE_0,
      Q_BLUE_1,
      Q_CLOSE_2_0,
      Q_CLOSE_2_1,
      Q_CLOSE_2_2,
      Q_CLOSE_2_3,
      Q_CLOSE_3_0,
      Q_CLOSE_3_1,
      Q_CLOSE_3_2,
      Q_CLOSE_3_3};

#define TET_COPY                        1
#ifdef TET_RULESET
#define FULL_REFRULE        Pattern2Rule[TETRAHEDRON][0x3F]
#define FULL_REFRULE_0_5    (Pattern2Rule[TETRAHEDRON][0x3F]+1)
#define FULL_REFRULE_1_3    (Pattern2Rule[TETRAHEDRON][0x3F]+2)
#define FULL_REFRULE_2_4    (Pattern2Rule[TETRAHEDRON][0x3F]+0)
#define TET_RED                         FULL_REFRULE
#define TET_RED_0_5                     FULL_REFRULE_0_5
#define TET_RED_1_3                     FULL_REFRULE_1_3
#define TET_RED_2_4                     FULL_REFRULE_2_4
#else
#define TET_RED                         2
#define TET_RED_0_5                     3
#define TET_RED_1_3                     4
#define TET_RED_2_4                     2
#define FULL_REFRULE        TET_RED
#define FULL_REFRULE_0_5    TET_RED_0_5
#define FULL_REFRULE_1_3    TET_RED_1_3
#define FULL_REFRULE_2_4    TET_RED_2_4
#define TET_RED_HEX                     5
#endif

enum {PYR_COPY       = 1,
      PYR_RED        = 2,
      PYR_BISECT_0_1 = 3,
      PYR_BISECT_0_2 = 4};

enum {PRI_COPY             = 1,
      PRI_RED              = 2,
      PRI_QUADSECT         = 3,
      PRI_BISECT_0_1       = 4,
      PRI_BISECT_0_2       = 5,
      PRI_BISECT_0_3       = 6,
      PRI_BISECT_1_2       = 7,
      PRI_BISECT_HEX0      = 8,
      PRI_BISECT_HEX1      = 9,
      PRI_BISECT_HEX2      = 10,
      PRI_RED_HEX          = 11,
      PRI_ROT_L            = 12,
      PRI_ROT_R            = 13,
      PRI_QUADSECT_HEXPRI0 = 14};

enum {HEXA_COPY           = 1,
      HEXA_RED            = 2,
      HEXA_BISECT_0_1     = 3,
      HEXA_BISECT_0_2     = 4,
      HEXA_BISECT_0_3     = 5,
      HEXA_QUADSECT_0     = 6,
      HEXA_QUADSECT_1     = 7,
      HEXA_QUADSECT_2     = 8,
      HEXA_TRISECT_0      = 9,
      HEXA_TRISECT_5      = 10,
      HEXA_BISECT_HEXPRI0 = 11,
      HEXA_BISECT_HEXPRI1 = 12};


/****************************************************************************/
/*                                                                                                                                                      */
/* macros for rules                                                                     */
/*                                                                                                                                                      */
/****************************************************************************/

#define TAG_OF_RULE(r)              ((r)->tag)
#define MARK_OF_RULE(r)             ((r)->mark)
#define CLASS_OF_RULE(r)            ((r)->rclass)
#define NSONS_OF_RULE(r)            ((r)->nsons)
#define SON_OF_RULE(r,s)            (&((r)->sons[(s)]))
#define PATTERN_OF_RULE(r,i)            ((r)->pattern[(i)])
#define PAT_OF_RULE(r)                          ((r)->pat)
#define SON_OF_NODE_OF_RULE(r,n)        ((r)->sonandnode[(n)][0])
#define SONNODE_OF_NODE_OF_RULE(r,n)((r)->sonandnode[(n)][0])
#define SON_TAG_OF_RULE(r,s)            ((r)->sons[(s)].tag)
#define SON_TAG(s)                                      ((s)->tag)
#define SON_CORNER_OF_RULE(r,s,n)       ((r)->sons[(s)].corners[(n)])
#define SON_CORNER(s,n)                         ((s)->corners[(n)])
#define SON_NB_OF_RULE(r,s,n)           ((r)->sons[(s)].nb[(n)])
#define SON_NB(s,n)                                     ((s)->nb[(n)])
#define SON_PATH_OF_RULE(r,s)           ((r)->sons[(s)].path)
#define SON_PATH(s)                                     ((s)->path)

/* macros for referencing of sons paths */
/* 4 high bits for no of neighbours to be passed */
#define PATHDEPTHMASK 0xF0000000
#define PATHDEPTHSHIFT 28
#define PATHDEPTH(i)                            (((i) & PATHDEPTHMASK)>>PATHDEPTHSHIFT)
#define SETPATHDEPTH(i,val)             (i) = ((i)&(~PATHDEPTHMASK))|(((val)<<PATHDEPTHSHIFT)&PATHDEPTHMASK)

/* 3 bits at position n for element side */
#define NEXTSIDEMASK 0x00000007
#define NEXTSIDE(i,n)                           (((i) & (NEXTSIDEMASK<<(3*(n))))>>(3*(n)))
#define SETNEXTSIDE(i,n,val)            (i) = ((i)&(~(NEXTSIDEMASK<<(3*(n)))))|(((val)&NEXTSIDEMASK)<<(3*(n)))

#define MAX_PATH_DEPTH                          (PATHDEPTHSHIFT/3)

#define NOCLASS(c)                                      ((c) == 0)
#define YELLOWCLASS(c)                          ((c) & 1)
#define GREENCLASS(c)                           ((c) & 2)
#define REDCLASS(c)                                     ((c) & 3)
#define SWITCHCLASS(c)                          ((c) & 4)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef INT (*FULLREFRULEPTR)(ELEMENT *);

typedef struct {
  INT type;                                    /* one of the variable types above          */
  INT locked;                                  /* may not be changed or deleted            */
  union envitem *next;
  union envitem *previous;                     /* double linked list of environment items  */
  char name[NS_PREFIX NAMESIZE];                         /* name of that item (view)                 */

  /* full ref rule spezific stuff */
  FULLREFRULEPTR theFullRefRule;               /* the best full refrule                    */

} FULLREFRULE;

struct sondata {
  SHORT tag;                                     /* which element type is the son                           */
  SHORT corners[MAX_CORNERS_OF_ELEM_DIM];        /* corners of the son                                  */
  SHORT nb[MAX_SIDES_OF_ELEM_DIM];               /* neighbors of this son                               */
  /* < FATHER_SIDE_OFFSET if nb has same father  */
  /* >= FATHER_SIDE_OFFSET if nb has other father*/
  INT path;                                      /* path used in GetSons() for tetras                   */
};

struct refrule {
  SHORT tag;                                                        /* which element type this rule can refine */
  SHORT mark;                                                       /* the mark of this rule                   */
  SHORT rclass;                                                     /* class of rule:3bits for COPY, IREG, REG */
  SHORT nsons;                                              /* number of sons rule creates             */
  SHORT pattern[MAX_NEW_CORNERS_DIM];                       /* stores which edges are refined          */
  INT pat;                                                      /* bitwise format of pattern               */
  SHORT sonandnode[MAX_NEW_CORNERS_DIM][2];                 /* for each new node the number of the son */
  /* and the local node number of the node   */
  struct sondata sons[MAX_SONS_DIM];
};

typedef struct sondata SONDATA;
typedef struct refrule REFRULE;


/****************************************************************************/
/*                                                                                                                                                      */
/* definition of exported global variables                                                                      */
/*                                                                                                                                                      */
/****************************************************************************/

extern INT MaxRules[TAGS];
extern INT MaxNewCorners[TAGS];
extern INT MaxNewEdges[TAGS];
extern INT CenterNodeIndex[TAGS];
extern REFRULE                  *RefRules[TAGS];
extern SHORT                    *Pattern2Rule[TAGS];
extern FULLREFRULEPTR theFullRefRule;


/****************************************************************************/
/*                                                                                                                                                      */
/* function declarations                                                                                                        */
/*                                                                                                                                                      */
/****************************************************************************/

INT ShowRefRule         (INT tag, INT nb);
INT ShowRefRuleX        (INT tag, INT nb, PrintfProcPtr Printf);

INT                     InitRuleManager                 (void);
INT                     Patterns2Rules                  (ELEMENT *theElement,INT pattern);
ELEMENT         *ELEMENT_TO_MARK                (ELEMENT *theElement);

#ifdef __THREEDIM__
INT             SetAlignmentPtr                 (MULTIGRID *theMG, EVECTOR *direction);
INT             GetRule_AnisotropicRed  (ELEMENT *theElement, INT *Rule);
#endif

END_UGDIM_NAMESPACE

#endif
