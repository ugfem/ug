// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugrefine3d.h													*/
/*																			*/
/* Purpose:   unstructured grid refinement header file						*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: klaus@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7007									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   09.03.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UGREFINE3__
#define __UGREFINE3__

#ifndef __gm__
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

#define LEAFELEM(e)     (REFINE(e)==0)
#define LEAFNODE(n)     (SONNODE(n)==NULL)

/* refinement rule numbers */
#define  NOREFRULE                       0
#define  COPY_REFRULE            1
#define  FULL_REFRULE            PatternToRefrule[0x3F]
#define  FULL_REFRULE_0_5        (PatternToRefrule[0x3F]+1)
#define  FULL_REFRULE_1_3        (PatternToRefrule[0x3F]+2)
#define  FULL_REFRULE_2_4        (PatternToRefrule[0x3F]+0)


#define MAXEDGES        16                      /* max. no of type 1/2 edges				*/
#define NEWCORNERS      7                       /* midpoints on edges plus center			*/
#define NOTUSED         -1                      /* SHORT has to be signed!					*/
#define NO_CENTER_NODE  NOTUSED

/* enumeration type for element class, now in gm.h ! */

/* other macros */
#define FATHER_SIDE_OFFSET      20 /* greater values indicate outside faces     */

/* macros specifying mark used in MarkForRefinement, this is now in gm.h ! */
#define NO_MARK                  0
#define FULL_REF_MARK    1
#define COPY_REF_MARK    2

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

#define FULLREFRULE_VAR          82
typedef INT (*FULLREFRULEPTR)(ELEMENT *);

typedef struct {
  INT type;                                                                                     /* one of the variable types above								*/
  INT locked;                                                                           /* may not be changed or deleted								*/
  union envitem *next;
  union envitem *previous;                                                      /* double linked list of environment items						*/
  char name[NAMESIZE];                                                          /* name of that item (view)                                                                     */

  /* full ref rule spezific stuff */
  FULLREFRULEPTR theFullRefRule;                        /* the best full refrule						*/

} FULLREFRULE;

struct sondata {
  SHORT corners[MAX_SIDES_OF_ELEM];                /* corners of the son */
  SHORT nb[MAX_SIDES_OF_ELEM];
  INT path;
};

/* edge type  : 1 = inner edge of the father
                                2 = inner edge of one side
                                3 = half an edge of the father element
                                4 = edge of the father itself
 */

struct edgedata {
  SHORT type;                              /* interior edge of the tetra or only of one side */
  SHORT from;                              /* indices of the first endpoint (0..9)			 */
  SHORT to;                                        /* indices of the second endpoint (4..9)              */
  SHORT side;                              /* side, for which edge is interior, if type 2	 */
};

struct refrule {
  SHORT nsons;
  SHORT pattern[MAX_EDGES_OF_ELEM];
  SHORT pat;
  SHORT sonandnode[NEWCORNERS][2];
  struct edgedata edges[MAXEDGES];
  struct sondata sons[MAX_SONS];
};


typedef struct sondata SONDATA;
typedef struct edgedata EDGEDATA;
typedef struct refrule REFRULE;


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

extern REFRULE *Rules;

/****************************************************************************/
/*																			*/
/* control word entries                                                                                                         */
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT InitRefine3d (void);
INT ShowRefRule (INT nb);

#endif
