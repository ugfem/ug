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
/* not used, kb 961216
 #define KEEP_VECTOR  0*/ /* this is a node with vector */
/*#define DEL_VECTOR  1*/  /* this is a node without vector */



/* macros for processor-synchronized output */
#define SYNC_ALL   { int _p; for(_p=0; _p<procs; _p++) { \
                       Synchronize(); if(_p==me) {
#define SYNC_CONTEXT   { int _p; for(_p=0; _p<procs; _p++) { \
                           Synchronize(); if((_p==me)&&CONTEXT(me)) {
#define SYNC_END  }}}


#define UGTYPE(t)          (dddctrl.ugtypes[(t)])
#define DDDTYPE(t)         (dddctrl.types[(t)])
#define HAS_DDDHDR(t)      (dddctrl.dddObj[(t)])
#define CONTEXT(p)         (dddctrl._context[(p)])

#define DDD_DOMAIN_DATA    DDD_USER_DATA+1
#define DDD_EXTRA_DATA     DDD_USER_DATA+2

/* macros for ddd object info */
/* for elements */
#define EPRIO(e)                                                DDD_InfoPriority(PARHDRE(e))
#define SETEPRIO(e,p)                                   DDD_PrioritySet(PARHDRE(e),p)
#define EMASTER(e)                                              (DDD_InfoPriority(PARHDRE(e)) == PrioMaster)
#define EGHOST(e)                                               (DDD_InfoPriority(PARHDRE(e))==PrioGhost || \
                                                                 DDD_InfoPriority(PARHDRE(e))==PrioVGhost)
#define EVGHOST(e)                                              (DDD_InfoPriority(PARHDRE(e))==PrioVGhost)
#define EHGHOST(e)                                              (DDD_InfoPriority(PARHDRE(e))==PrioGhost)
#define EGID(e)                                                 DDD_InfoGlobalId(PARHDRE(e))
#define EPROCLIST(e)                                    DDD_InfoProcList(PARHDRE(e))
#define EATTR(e)                                                DDD_InfoAttr(PARHDRE(e))
#define XFEREDELETE(e)                                  DDD_XferDeleteObj(PARHDRE(e))
#define XFERECOPY(e,dest,prio)                  DDD_XferCopyObj(PARHDRE(e),dest,prio)
#define XFERECOPYX(e,dest,prio,size)    DDD_XferCopyObjX(PARHDRE(e),dest,prio,size)

/* for nodes, vectors, edges (edges only 3D) */
#define PRIO(e)                                                 DDD_InfoPriority(PARHDR(e))
#define SETPRIO(e,p)                                    DDD_PrioritySet(PARHDR(e),p)
#define MASTER(e)                                               (DDD_InfoPriority(PARHDR(e))==PrioMaster || \
                                                                 DDD_InfoPriority(PARHDR(e))==PrioBorder)
#define GHOST(e)                                                (DDD_InfoPriority(PARHDR(e))==PrioGhost || \
                                                                 DDD_InfoPriority(PARHDR(e))==PrioVGhost)
#define VGHOST(e)                                               (DDD_InfoPriority(PARHDR(e))==PrioVGhost)
#define HGHOST(e)                                               (DDD_InfoPriority(PARHDR(e))==PrioGhost)
#define GID(e)                                                  DDD_InfoGlobalId(PARHDR(e))
#define PROCLIST(e)                                             DDD_InfoProcList(PARHDR(e))
#define ATTR(e)                                                 DDD_InfoAttr(PARHDR(e))
#define XFERDELETE(e)                                   DDD_XferDeleteObj(PARHDR(e))
#define XFERCOPY(e,dest,prio)                   DDD_XferCopyObj(PARHDR(e),dest,prio)
#define XFERCOPYX(e,dest,prio,size)             DDD_XferCopyObjX(PARHDR(e),dest,prio,size)

/* for vertices */
#define VXPRIO(e)                                               DDD_InfoPriority(PARHDRV(e))
#define SETVXPRIO(e,p)                                  DDD_PrioritySet(PARHDRV(e),p)
#define VXMASTER(e)                                             (DDD_InfoPriority(PARHDRV(e))==PrioMaster || \
                                                                 DDD_InfoPriority(PARHDRV(e))==PrioBorder)
#define VXGHOST(e)                                              (DDD_InfoPriority(PARHDRV(e))==PrioGhost || \
                                                                 DDD_InfoPriority(PARHDRV(e))==PrioVGhost)
#define VXVGHOST(e)                                             (DDD_InfoPriority(PARHDRV(e))==PrioVGhost)
#define VXHGHOST(e)                                             (DDD_InfoPriority(PARHDRV(e))==PrioGhost)
#define VXGID(e)                                                DDD_InfoGlobalId(PARHDRV(e))
#define VXPROCLIST(e)                                   DDD_InfoProcList(PARHDRV(e))
#define VXATTR(e)                                               DDD_InfoAttr(PARHDRV(e))
#define XFERVXDELETE(e)                                 DDD_XferDeleteObj(PARHDRV(e))
#define XFERVXCOPY(e,dest,prio)                 DDD_XferCopyObj(PARHDRV(e),dest,prio,size)
#define XFERVXCOPYX(e,dest,prio,size)   DDD_XferCopyObjX(PARHDRV(e),dest,prio,size)

/* macros for priorities */
/* for elements */
#define EMASTERPRIO(p)                                  (p==PrioMaster)
#define EGHOSTPRIO(p)                                   (p==PrioGhost || p==PrioVGhost)
#define EVGHOSTPRIO(p)                                  (p==PrioVGhost)
#define EHGHOSTPRIO(p)                                  (p==PrioGhost)

/* for nodes, vertices, vectors, edges (edges only 3D) */
#define MASTERPRIO(p)                                   (p==PrioMaster || p==PrioBorder)
#define GHOSTPRIO(p)                                    (p==PrioGhost || p==PrioVGhost)
#define VGHOSTPRIO(p)                                   (p==PrioVGhost)
#define HGHOSTPRIO(p)                                   (p==PrioGhost)


#define GIDFMT                                                  "%08x"

#define __EXCHANGE_CONNECTIONS__

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

extern DDD_TYPE TypeUnknown;

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
  INT  *_context;
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
int  ExitParallel (void);
void InitDDDTypes (void);
void InitCurrMG (MULTIGRID *);

/* from debugger.c */
void ddd_pstat (char *);
void ddd_DisplayContext (void);

/* from test.c */
void ddd_test ();

/* from handler.c */
void            ddd_HandlerInit (INT);
DDD_TYPE        NFatherObjType  (DDD_OBJ obj, DDD_OBJ ref);

/* from lbrcb.c */
int BalanceGridRCB (MULTIGRID *, int);

/* from gridcons.c */
void    SetOverlapPriorities            (GRID *theGrid);
void    ConstructConsistentGrid         (GRID *theGrid);
INT             CheckInterfaces                         (GRID *theGrid);

/* from transfer.c */
int TransferGridFromLevel (MULTIGRID *theMG, INT level);

#endif /* ModelP */
#endif /* __PARALLEL_H__ */
