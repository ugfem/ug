// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  parallel.h													*/
/*																			*/
/* Purpose:   defines for parallel ugp version 3							*/
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: stefan@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   09.05.95 begin, ugp version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#ifndef __HEAPS__
#include "heaps.h"
#endif

#ifdef ModelP
#include "ppif.h"
#include "ddd.h"
#endif

#include "gm.h"
#include "pargm.h"


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/


#define MAXDDDTYPES   32

enum HandlerSets
{
  HSET_XFER = 0,
  HSET_REFINE
};


#ifdef ModelP

/* defines for control word entries */
#define TOUCHED     2 /* object is eventually to be copied  */
#define COPY        1 /* object is definitely to be copied   */
#define CLEAR       0 /* clear xfer flag					*/

/* CE for nodes */
#define KEEP_VECTOR  0 /* this is a node with vector */
#define DEL_VECTOR  1  /* this is a node without vector */




/* macros for processor-synchronized output */
#define SYNC_ALL   { int _p; for(_p=0; _p<procs; _p++) { \
                       Synchronize(); if(_p==me) {
#define SYNC_CONTEXT   { int _p; for(_p=0; _p<procs; _p++) { \
                           Synchronize(); if((_p==me)&&CONTEXT(me)) {
#define SYNC_END  }}}


#define UGTYPE(t)          (dddctrl.ugtypes[(t)])
#define DDDTYPE(t)         (dddctrl.types[(t)])
#define HAS_DDDHDR(t)      (dddctrl.dddObj[(t)])
#define CONTEXT(p)         (dddctrl.context[(p)])


#endif /* ModelP */


/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

#ifdef ModelP
/* DDD objects */
extern DDD_TYPE TypeVector;
extern DDD_TYPE TypeIVertex, TypeBVertex;
extern DDD_TYPE TypeNode;

#ifdef __TWODIM__
extern DDD_TYPE TypeTrElem, TypeTrBElem,
                TypeQuElem, TypeQuBElem;
#endif

#ifdef __THREEDIM__
extern DDD_TYPE TypeTeElem, TypeTeBElem;
extern DDD_TYPE TypePyElem, TypePyBElem;
extern DDD_TYPE TypePrElem, TypePrBElem;
extern DDD_TYPE TypeHeElem, TypeHeBElem;
#endif

/* DDD data objects */
extern DDD_TYPE TypeMatrix;
extern DDD_TYPE TypeBndP;
extern DDD_TYPE TypeEdge;
extern DDD_TYPE TypeBndS;



/* DDD Global Controls */
typedef struct
{
  /* data from ug */
  MULTIGRID *currMG;
  int nodeData;
  int edgeData;
  int elemData;
  int sideData;

  /* data for memmgr */
  MULTIGRID *memmgrMG;

  INT ugtypes[MAXDDDTYPES];                  /* dddtype -> ugtype */
  DDD_TYPE types[MAXOBJECTS];                /* ugtype -> dddtype */
  int dddObj[MAXOBJECTS];

  /* status of DDDIF */
  int allTypesDefined;

  /* context information, will be allocated with size=procs */
  int  *context;
} DDD_CTRL;

extern DDD_CTRL dddctrl;

#endif

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

#ifdef ModelP

/* from initddd.c */
int  InitParallel (int *, char ***);
void ExitParallel (void);
void InitDDDTypes (void);
void InitCurrMG (MULTIGRID *);

/* from debugger.c */
void ddd_pstat (int);
void ddd_DisplayContext (void);

/* from test.c */
void ddd_test ();

/* from handler.c */
void ddd_HandlerInit (INT);

/* from lbrcb.c */
int BalanceGrid (MULTIGRID *);

/* from gridcons.c */
void dddif_SetBorderPriorities (GRID *);


#endif /* ModelP */
#endif /* __PARALLEL_H__ */
