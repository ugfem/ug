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

/* edges */
#define PATTERN_CE                                      39
#define PATTERN_SHIFT                           10
#define PATTERN_LEN                             1
#define PATTERN(p)                                      CW_READ(p,PATTERN_CE)
#define SETPATTERN(p,n)                         CW_WRITE(p,PATTERN_CE,n)

#define ADDPATTERN_CE                           40
#define ADDPATTERN_SHIFT                        11
#define ADDPATTERN_LEN                          1
#define ADDPATTERN(p)                           CW_READ(p,ADDPATTERN_CE)
#define SETADDPATTERN(p,n)                      CW_WRITE(p,ADDPATTERN_CE,n)


/* element */
#define REFINE_CE                                               45
#define REFINE_SHIFT                                    0
#define REFINE_LEN                                              8
#define REFINE(p)                                               CW_READ(p,REFINE_CE)
#define SETREFINE(p,n)                                  CW_WRITE(p,REFINE_CE,n)

#define MARK_CE                                                 50
#define MARK_SHIFT                                              0
#define MARK_LEN                                                8
#define MARK(p)                                                 CW_READ(p,MARK_CE)
#define SETMARK(p,n)                                    CW_WRITE(p,MARK_CE,n)

#define COARSEN_CE                                              51
#define COARSEN_SHIFT                                   10
#define COARSEN_LEN                                     1
#define COARSEN(p)                                              CW_READ(p,COARSEN_CE)
#define SETCOARSEN(p,n)                                 CW_WRITE(p,COARSEN_CE,n)

#define DECOUPLED_CE                                    53
#define DECOUPLED_SHIFT                                 12
#define DECOUPLED_LEN                                   1
#define DECOUPLED(p)                                    CW_READ(p,DECOUPLED_CE)
#define SETDECOUPLED(p,n)                               CW_WRITE(p,DECOUPLED_SHIFT,n)

#define REFINECLASS_CE                                  49
#define REFINECLASS_SHIFT                               15
#define REFINECLASS_LEN                                 2
#define REFINECLASS(p)                                  CW_READ(p,REFINECLASS_CE)
#define SETREFINECLASS(p,n)                     CW_WRITE(p,REFINECLASS_CE,n)

#define UPDATE_GREEN_CE                                 54
#define UPDATE_GREEN_SHIFT                              8
#define UPDATE_GREEN_LEN                                1
#define UPDATE_GREEN(p)                                 CW_READ(p,UPDATE_GREEN_CE)
#define SETUPDATE_GREEN(p,n)                    CW_WRITE(p,UPDATE_GREEN_CE,n)

#define SIDEPATTERN_CE                                  55
#define SIDEPATTERN_SHIFT                               0
#define SIDEPATTERN_LEN                                 6
#define SIDEPATTERN(p)                                  CW_READ(p,SIDEPATTERN_CE)
#define SETSIDEPATTERN(p,n)                     CW_WRITE(p,SIDEPATTERN_CE,n)

#define MARKCLASS_CE                                    58
#define MARKCLASS_SHIFT                                 13
#define MARKCLASS_LEN                                   2
#define MARKCLASS(p)                                    CW_READ(p,MARKCLASS_CE)
#define SETMARKCLASS(p,n)                               CW_WRITE(p,MARKCLASS_CE,n)

/* macros for refineinfo */
#define REFINEINFO(mg)                                  refine_info
#define MARKCOUNT(r)                                    r.markcount
#define SETMARKCOUNT(r,n)                               r.markcount = n
#define PREDNEW(r)                                              r.predicted_new
#define SETPREDNEW(r,n)                                 r.predicted_new = n
#define PREDMAX(r)                                              r.predicted_max
#define SETPREDMAX(r,n)                                 r.predicted_max = n

/****************************************************************************/
/*																			*/
/* typedefs																	*/
/*																			*/
/****************************************************************************/

typedef struct refineinfo
{
  float markcount;                      /* count of currently marked elements           */
  float predicted_new;          /* count of elements, which would be created    */
  float predicted_max;          /* count of elements which can be created       */
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
#ifdef ModelP
INT             Get_Sons_of_ElementSide                                 (ELEMENT *theElement, INT side, INT *Sons_of_Side, ELEMENT *SonList[MAX_SONS], INT *SonSides, INT NeedSons);
#endif
INT     Connect_Sons_of_ElementSide                     (GRID *theGrid, ELEMENT *theElement, INT side, INT Sons_of_Side, ELEMENT **Sons_of_Side_List, INT *SonSides);

#endif
