// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  parallel.h													*/
/*																			*/
/* Purpose:   defines for parallel ugp version 3							*/
/*																			*/
/* Author:	  Stefan Lang                                                                           */
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

#ifndef __GM__
#include "gm.h"
#endif

#ifdef ModelP
#ifndef __PPIF__
#include "ppif.h"
#endif
#endif

#ifndef __DDD__
#include "ddd.h"
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

#ifdef __TWODIM__
#define TAG2TYPE(t)
#define TYPE2TAG(t)
#endif

#ifdef __THREEDIM__
#define TAG2TYPE(t)
#define TYPE2TAG(t)
#endif

#ifdef ModelP
#define GetMemVector(h,s,f)       GetObjMem(TypeVector,(h),(s),(f))
#define GetMemIVertex(h,s,f)      GetObjMem(TypeIVertex,(h),(s),(f))
#define GetMemBVertex(h,s,f)      GetObjMem(TypeBVertex,(h),(s),(f))
#define GetMemNode(h,s,f)         GetObjMem(TypeNode,(h),(s),(f))
#ifdef __TWODIM__
#define GetMemIElement(h,s,f)     GetObjMem(((tag==3) ? 4 : 6),(h),(s),(f))
#define GetMemBElement(h,s,f)     GetObjMem(((tag==3) ? 5 : 7),(h),(s),(f))
#endif /* __TWODIM__ */
#ifdef __THREEDIM__
#define GetMemIElement(h,s,f)     GetObjMem(((tag==4) ? 4 : -1),(h),(s),(f))
#define GetMemBElement(h,s,f)     GetObjMem(((tag==4) ? 5 : -1),(h),(s),(f))
#endif /* __THREEDIM__ */

/* defines for control word entries */
#define DDD_OFFSET  4 /* size of DDD_HEADER in unsigned INT */
#define TOUCHED     2 /* object is eventually to be copied  */
#define COPY        1 /* object is definetly to be copied   */
#define CLEAR       0 /* clear xfer flag					*/

#define VECTOR_CTRL(p)          (p->control)
#define VERTEX_CTRL(p)          (p->iv.control)
#define NODE_CTRL(p)            (p->control)
#define ELEMENT_CTRL(p)         (p->ge.control)

#define VECTOR_ID(p)            (p->id)
#define VERTEX_ID(p)            (p->iv.id)
#define NODE_ID(p)                      (p->id)
#define ELEMENT_ID(p)           (p->ge.id)

#else   /* ModelS */

#define GetMemVector(h,s,f)       GetMem(h,s,f)
#define GetMemIVertex(h,s,f)      GetMem(h,s,f)
#define GetMemBVertex(h,s,f)      GetMem(h,s,f)
#define GetMemNode(h,s,f)         GetMem(h,s,f)
#define GetMemIElement(h,s,f)     GetMem(h,s,f)
#define GetMemBElement(h,s,f)     GetMem(h,s,f)

#define DDD_OFFSET      0 /* size of DDD_HEADER in unsigned INT */

#define VECTOR_CTRL(p)          CTRL(p)
#define VERTEX_CTRL(p)          CTRL(p)
#define NODE_CTRL(p)            CTRL(p)
#define ELEMENT_CTRL(p)         CTRL(p)

#define VECTOR_ID(p)            ID(p)
#define VERTEX_ID(p)            ID(p)
#define NODE_ID(p)                      ID(p)
#define ELEMENT_ID(p)           ID(p)
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
#endif

/* DDD data objects */
extern DDD_TYPE TypeConnection;
extern DDD_TYPE TypeVSegment;
extern DDD_TYPE TypeEdge;
extern DDD_TYPE TypeElementSide;

/* debug processor */
extern int theDebugProc;
#endif

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

#ifdef ModelP
/* from initddd.c */
int  InitParallel (int, char **);
void ExitParallel (void);

/* from memmgr.c */
void * GetObjMem (DDD_TYPE, HEAP *, MEM, int);
void InitMemMgr(HEAP *);

/* from debug.c */
void ddd_pstat (int cmd);

/* from test.c */
void ddd_test ();

/* from handler.c */
void ddd_handlerInit (void);
void InitCurrMG(MULTIGRID *MG)
#endif /* ModelP */

#endif /* PARALLEL */
