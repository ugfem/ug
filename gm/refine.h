// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  refine.h														*/
/*																			*/
/* Purpose:   definitions for two AND three dimensional refinement			*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*/
/*																			*/
/* History:   09.03.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __REFINE__
#define __REFINE__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __GM__
#include "gm.h"
#endif

#ifndef __ALGEBRA__
#include "algebra.h"            /* just for ALGEBRA_N_CE */
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define NOTUSED    -1                      /* SHORT has to be signed!       */
#define NO_CENTER_NODE     NOTUSED


/****************************************************************************/
/*																			*/
/* control word definitions                                                                                             */
/*																			*/
/****************************************************************************/

enum REFINE_CE {

  PATTERN_CE = ALGEBRA_N_CE,            /* continue after algebra.h entries		*/
  ADDPATTERN_CE,
  REFINE_CE,
  MARK_CE,
  COARSEN_CE,
  DECOUPLED_CE,
  REFINECLASS_CE,
  UPDATE_GREEN_CE,
  SIDEPATTERN_CE,
  MARKCLASS_CE,

  REFINE_N_CE
};

/* edges */
#define PATTERN_SHIFT                           10
#define PATTERN_LEN                             1
#define PATTERN(p)                                      CW_READ(p,PATTERN_CE)
#define SETPATTERN(p,n)                         CW_WRITE(p,PATTERN_CE,n)

#define ADDPATTERN_SHIFT                        11
#define ADDPATTERN_LEN                          1
#define ADDPATTERN(p)                           CW_READ(p,ADDPATTERN_CE)
#define SETADDPATTERN(p,n)                      CW_WRITE(p,ADDPATTERN_CE,n)


/* element */
#define REFINE_SHIFT                                    0
#define REFINE_LEN                                              8
#define REFINE(p)                                               CW_READ(p,REFINE_CE)
#define SETREFINE(p,n)                                  CW_WRITE(p,REFINE_CE,n)

#define MARK_SHIFT                                              0
#define MARK_LEN                                                8
#define MARK(p)                                                 CW_READ(p,MARK_CE)
#define SETMARK(p,n)                                    CW_WRITE(p,MARK_CE,n)

#define COARSEN_SHIFT                                   10
#define COARSEN_LEN                                     1
#define COARSEN(p)                                              CW_READ(p,COARSEN_CE)
#define SETCOARSEN(p,n)                                 CW_WRITE(p,COARSEN_CE,n)

#define DECOUPLED_SHIFT                                 12
#define DECOUPLED_LEN                                   1
#define DECOUPLED(p)                                    CW_READ(p,DECOUPLED_CE)
#define SETDECOUPLED(p,n)                               CW_WRITE(p,DECOUPLED_SHIFT,n)

#define REFINECLASS_SHIFT                               15
#define REFINECLASS_LEN                                 2
#define REFINECLASS(p)                                  CW_READ(p,REFINECLASS_CE)
#define SETREFINECLASS(p,n)                     CW_WRITE(p,REFINECLASS_CE,n)

#define UPDATE_GREEN_SHIFT                              8
#define UPDATE_GREEN_LEN                                1
#define UPDATE_GREEN(p)                                 CW_READ(p,UPDATE_GREEN_CE)
#define SETUPDATE_GREEN(p,n)                    CW_WRITE(p,UPDATE_GREEN_CE,n)

#define SIDEPATTERN_SHIFT                               0
#define SIDEPATTERN_LEN                                 6
#define SIDEPATTERN(p)                                  CW_READ(p,SIDEPATTERN_CE)
#define SETSIDEPATTERN(p,n)                     CW_WRITE(p,SIDEPATTERN_CE,n)

#define MARKCLASS_SHIFT                                 13
#define MARKCLASS_LEN                                   2
#define MARKCLASS(p)                                    CW_READ(p,MARKCLASS_CE)
#define SETMARKCLASS(p,n)                               CW_WRITE(p,MARKCLASS_CE,n)

/* macros for refineinfo */
#define REFINEINFO(mg)                                  refine_info
#define REFINESTEP(r)                                   r.step
#define SETREFINESTEP(r,s)                              r.step = s
#define MARKCOUNT(r)                                    r.markcount[r.step]
#define SETMARKCOUNT(r,n)                               r.markcount[r.step] = n
#define PREDNEW0(r)                                             r.predicted_new[r.step][0]
#define SETPREDNEW0(r,n)                                r.predicted_new[r.step][0] = n
#define PREDNEW1(r)                                             r.predicted_new[r.step][1]
#define SETPREDNEW1(r,n)                                r.predicted_new[r.step][1] = n
#define PREDNEW2(r)                                             r.predicted_new[r.step][2]
#define SETPREDNEW2(r,n)                                r.predicted_new[r.step][2] = n
#define REAL(r)                                                 r.real[r.step]
#define SETREAL(r,n)                                    r.real[r.step] = n
#define PREDMAX(r)                                              r.predicted_max[r.step]
#define SETPREDMAX(r,n)                                 r.predicted_max[r.step] = n

/****************************************************************************/
/*																			*/
/* typedefs																	*/
/*																			*/
/****************************************************************************/

typedef struct refineinfo
{
  INT step;                                       /* count of calls to RefineMultiGrid      */
  float markcount[100];                   /* count of currently marked elements     */
  float predicted_new[100][3];        /* count of elements, would be created    */
  float real[100];                                /* count of elements before refinement    */
  float predicted_max[100];               /* count of elements which can be created */
} REFINEINFO;

typedef INT (*Get_Sons_of_ElementSideProcPtr)(ELEMENT *theElement, INT side, INT *Sons_of_Side,ELEMENT *SonList[MAX_SONS], INT *SonSides, INT NeedSons);

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

extern REFINEINFO refine_info;

/****************************************************************************/
/*																			*/
/* functions exported														*/
/*																			*/
/****************************************************************************/

INT     Set_Get_Sons_of_ElementSideProc                 (Get_Sons_of_ElementSideProcPtr Proc);
INT     GetSonSideNodes                                                 (ELEMENT *theElement, INT side, INT *nodes, NODE *SideNodes[MAX_SIDE_NODES]);
INT             Get_Sons_of_ElementSide                                 (ELEMENT *theElement, INT side, INT *Sons_of_Side, ELEMENT *SonList[MAX_SONS], INT *SonSides, INT NeedSons);
INT     Connect_Sons_of_ElementSide                     (GRID *theGrid, ELEMENT *theElement, INT side, INT Sons_of_Side, ELEMENT **Sons_of_Side_List, INT *SonSides, INT notHanging);

#endif
