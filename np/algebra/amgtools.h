// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amgtools.h													*/
/*																			*/
/* Purpose:   header for amgtools.c			                                                                */
/*																			*/
/* Author:	  Nicolas Neuss                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de				                        */
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __AMGTOOLS__
#define __AMGTOOLS__

#include "np.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAXNEIGHBORS 128

#define ASBEFORE 40
#define COARSEFINE 41
#define FINECOARSE 42

/* we use the MUSED-Bit for Matrices */
#define STRONG(p) MUSED(p)
#define SETSTRONG(p,n) SETMUSED(p,n)

#define AVCOARSEMASK 0x00000001
#define AVCOARSESHIFT 0
#define AVCOARSE(p) (((*((unsigned int *)(p))) & AVCOARSEMASK)>>AVCOARSESHIFT)
#define SETAVCOARSE(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~AVCOARSEMASK))|(((n)<<AVCOARSESHIFT)&AVCOARSEMASK)

#define AVFINEMASK 0x00000002
#define AVFINESHIFT 1
#define AVFINE(p) (((*((unsigned int *)(p))) & AVFINEMASK)>>AVFINESHIFT)
#define SETAVFINE(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~AVFINEMASK))|(((n)<<AVFINESHIFT)&AVFINEMASK)

#define AVTESTEDMASK 0x00000004
#define AVTESTEDSHIFT 2
#define AVTESTED(p) (((*((unsigned int *)(p))) & AVTESTEDMASK)>>AVTESTEDSHIFT)
#define SETAVTESTED(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~AVTESTEDMASK))|(((n)<<AVTESTEDSHIFT)&AVTESTEDMASK)

#define STRONG_IN(p) ((p)->StronglyInfluencing)
#define STRONG_OUT(p) ((p)->StronglyInfluenced)
#define VECT(avect) (avect->vect)

/* some useful macros for avects (see below) */
#ifndef ModelP
        #define ELIMINATE_LIST1(ls,p) {if ((p)->pred!=NULL) (p)->pred->succ=(p)->succ;else ls=(p)->succ;if ((p)->succ!=NULL) (p)->succ->pred=(p)->pred;}
        #define ELIMINATE_LIST2(ls,le,p) {if ((p)->pred!=NULL) (p)->pred->succ=(p)->succ;else ls=(p)->succ;if ((p)->succ!=NULL) (p)->succ->pred=(p)->pred;else le=(p)->pred;}

        #define ADDATSTART_LIST1(ls,p) {(p)->succ=ls; (p)->pred=NULL; if (ls!=NULL) ls->pred=(p);ls=(p);}
        #define ADDATSTART_LIST2(ls,le,p) {(p)->succ=ls; (p)->pred=NULL; if (ls!=NULL) ls->pred=(p);else le=(p);ls=(p);}

        #define ADDATEND_LIST2(ls,le,p) {(p)->pred=le; (p)->succ=NULL; if (le!=NULL) (le)->succ=(p);else ls=(p);le=(p);}

        #define ADDBEFORE_LIST2(ls,le,pa,p) {(p)->succ=(pa); if (((p)->pred=(pa)->pred)==NULL) ls=(p);else (p)->pred->succ=(p);(pa)->pred=p;}
        #define APPEND_LIST2(la,le,aa,ae) { if ((aa)!=NULL) {if ((la)==NULL) {la=(aa); le=(ae);} else {(le)->succ=(aa); (aa)->pred=(le); le=(ae);}} }
#else
/* does not yet run in parallel version */
        #define ELIMINATE_LIST1(ls,p) {(p)=NULL;}
        #define ELIMINATE_LIST2(ls,le,p) {(p)=NULL;}
        #define ADDATSTART_LIST1(ls,p) {(p)=NULL;}
        #define ADDATSTART_LIST2(ls,le,p) {(p)=NULL;}
        #define ADDATEND_LIST2(ls,le,p) {(p)=NULL;}
        #define ADDBEFORE_LIST2(ls,le,pa,p) {(p)=NULL;}
        #define APPEND_LIST2(la,le,aa,ae) {(aa)=NULL;}
#endif

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/* the following structure is used to be able to do arbitrarily complex
   list operations without destroying the previous ordering */
struct avector {
  INT control;
  INT StronglyInfluencing;
  INT StronglyInfluenced;
  struct avector *pred;
  struct avector *succ;
  /*	struct avector *pred2;
          struct avector *succ2; */
  VECTOR *vect;
};
typedef struct avector AVECTOR;


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT UnmarkAll             (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta);
INT MarkAll               (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta);
INT MarkAbsolute          (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta);
INT MarkRelative          (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta);
INT SetupInitialList      (GRID *theGrid, HEAP *theHeap, AVECTOR **initialSH, AVECTOR **initialEH);
INT DistributeInitialList (AVECTOR **La, AVECTOR **Le, AVECTOR **Ta, AVECTOR **Te, AVECTOR **Ua, AVECTOR **Ue);
INT CountStrongNeighbors  (AVECTOR *initialS, DOUBLE *avNrOfStrongNbsHnd, INT *maxNeighbors);
INT CoarsenRugeStueben    (GRID *theGrid);
INT CoarsenVanek          (GRID *theGrid);
INT IpRugeStueben         (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT IpVanek               (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT GalerkinCGMatrixFromInterpolation(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT SparsenCGMatrix       (GRID *theGrid, MATDATA_DESC *A, INT lumpFlag);
INT ReorderFineGrid       (GRID *theGrid, INT orderType);

#endif
