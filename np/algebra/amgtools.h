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


/* RCS_ID
   $Header$
 */

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

#define AVSKIPMASK 0x00000008
#define AVSKIPSHIFT 3
#define AVSKIP(p) (((*((unsigned int *)(p))) & AVSKIPMASK)>>AVSKIPSHIFT)
#define SETAVSKIP(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~AVSKIPMASK))|(((n)<<AVSKIPSHIFT)&AVSKIPMASK)

#define STRONG_IN(p) ((p)->StronglyInfluencing)
#define STRONG_OUT(p) ((p)->StronglyInfluenced)
#define VECT(avect) (avect->vect)

/****************************************************************************/
/* some useful macros for avects (see below)                                            */
/****************************************************************************/
#define ELIMINATE_LIST1(ls,p) {if ((p)->pred!=NULL) (p)->pred->succ=(p)->succ;else ls=(p)->succ;if ((p)->succ!=NULL) (p)->succ->pred=(p)->pred;}
#define ELIMINATE_LIST2(ls,le,p) {if ((p)->pred!=NULL) (p)->pred->succ=(p)->succ;else ls=(p)->succ;if ((p)->succ!=NULL) (p)->succ->pred=(p)->pred;else le=(p)->pred;}

#define ADDATSTART_LIST1(ls,p) {(p)->succ=ls; (p)->pred=NULL; if (ls!=NULL) ls->pred=(p);ls=(p);}
#define ADDATSTART_LIST2(ls,le,p) {(p)->succ=ls; (p)->pred=NULL; if (ls!=NULL) ls->pred=(p);else le=(p);ls=(p);}

#define ADDATEND_LIST2(ls,le,p) {(p)->pred=le; (p)->succ=NULL; if (le!=NULL) (le)->succ=(p);else ls=(p);le=(p);}

#define ADDBEFORE_LIST2(ls,le,pa,p) {(p)->succ=(pa); if (((p)->pred=(pa)->pred)==NULL) ls=(p);else (p)->pred->succ=(p);(pa)->pred=p;}
#define APPEND_LIST2(la,le,aa,ae) { if ((aa)!=NULL) {if ((la)==NULL) {la=(aa); le=(ae);} else {(le)->succ=(aa); (aa)->pred=(le); le=(ae);}} }

/****************************************************************************/
/* Useful macros for block arithmetic

   DETAILS:
   - pointers are copied (hopefully to registers) and access then is
     to subsequent components (thus not allowing arbitrary sparse matrices)
   - loops are down to zero to make the end criterion easier
   - Mainly the variable "scalar" should be a register, to make the scalar case
     faster. Therefore use the BLOCK_SETUP-Macro for initializing.

   IMPORTANT:
   Be careful with terminating semicolons in if-else constructs!            */
/****************************************************************************/

/* This macro checks that only mtype=0x0 occurs and that the matrix
   entries are subsequent in memory. It further defines and sets the
   (register) variables error, scalar, blockN, blockNN. */
#define BLOCK_SETUP(md) INT error;\
  register int scalar,blockN,blockNN;\
  blockN=MD_ROWS_IN_MTYPE(md,0);\
  {register int mt; for (mt=1; mt<NMATTYPES; mt++) if (MD_ROWS_IN_MTYPE(md,mt)!=0) break;\
   if ((blockN==0) || (mt!=NMATTYPES)) error=1;\
   else if ((MD_SUCC_COMP(md))==0) error=2;\
   else {error=0; blockNN = blockN*blockN;\
         if (blockN==1) scalar=1;else scalar=0;} }

#define BLOCK_CLEAR(A) {if (scalar) *A=0.0;\
                        else {register int i; register DOUBLE *A_=A;\
                              for (i=blockNN; i>0; i--) *A_++ = 0.0;}}

#define BLOCK_IDENTITY(A) {if (scalar) *A=1.0;\
                           else {register int i,j; register DOUBLE *A_=A;\
                                 *A_++ = 1.0; for (i=blockN-1; i>0; i--) {\
                                   for (j=blockN; j>0; j--) *A_++ = 0.0;\
                                   *A_++ = 1.0; }}}

#define BLOCK_SCALIDENTITY(a,A) {if (scalar) *A=a;\
                                 else {register int i,j; register DOUBLE *A_=A,a_=a;\
                                       *A_++ = a_; for (i=blockN-1; i>0; i--) {\
                                         for (j=blockN; j>0; j--) *A_++ = 0.0;\
                                         *A_++ = a_; }}}

#define BLOCK_COPY(A,B) {if (scalar) *B=*A;\
                         else {register int i; register DOUBLE *A_=A, *B_=B;\
                               for (i=blockNN; i>0; i--) *B_++ = *A_++;}}

#define BLOCK_ADD(A,B,C) {if (scalar) *C = *A + *B;\
                          else {register int i;  register DOUBLE *A_=A, *B_=B, *C_=C;\
                                for (i=blockNN; i>0; i--) *C_++ = (*A_++)+(*B_++);}}

#define BLOCK_ADD1(A,C) {if (scalar) *C += *A;\
                         else {register int i; register DOUBLE *A_=A,*C_=C;\
                               for (i=blockNN; i>0; i--) (*C_++) += (*A_++);}}

#define BLOCK_SCALE(a,A,C) {if (scalar) *C = a * (*A);\
                            else {register int i; register DOUBLE a_=a,*A_=A,*C_=C;\
                                  for (i=blockNN; i>0; i--) *C_++ = a_ * (*A_++);}}

#define BLOCK_SCALE1(a,A) {if (scalar) *A *= a;\
                           else {register int i; register DOUBLE a_=a,*A_=A;\
                                 for (i=blockNN; i>0; i--) *A_++ *= a_;}}

#define BLOCK_NORM(A,a) {register DOUBLE *A_=A;\
                         if (scalar) a = sqrt((*A_)*(*A_));\
                         else {register int i; register DOUBLE s,t;\
                               s = 0.0; for (i=blockNN; i>0; i--) {t = *A_++; s += t*t;}\
                               a = sqrt(s);} }

#define BLOCK_MUL_NNN_TMPLT(A,B,C,OPERAND)  {if (scalar) *C OPERAND (*A) * (*B);\
                                             else {register int i,j,k; register DOUBLE s;\
                                                   register DOUBLE *A_=A, *B_=B, *C_=C, *A_2, *B_2;\
                                                   for (i=blockN; i>0; i--) {\
                                                     for (j=blockN; j>0; j--) {\
                                                       A_2=A_; B_2=B_; s=0.0;\
                                                       for (k=blockN; k>0; k--) {s += (*A_2++) * (*B_2); B_2+=blockN;\
                                                       } *C_++ OPERAND s; B_++;\
                                                     } A_ += blockN; B_ -=blockN;\
                                                   } } }

#define BLOCK_MUL_NTN_TMPLT(A,B,C,OPERAND) {if (scalar) *C OPERAND (*A) * (*B);\
                                            else {register int i,j,k; register DOUBLE s;\
                                                  register DOUBLE *A_=A, *B_=B, *C_=C;\
                                                  for (i=blockN; i>0; i--) {\
                                                    for (j=blockN; j>0; j--) {\
                                                      s=0.0; for (k=blockN; k>0; k--) {s += (*A_++) * (*B_++);}\
                                                      *C_++ OPERAND s; A_ -= blockN; \
                                                    } A_ += blockN; B_ -= blockNN;\
                                                  } } }

#define BLOCK_MUL_TNN_TMPLT(A,B,C,OPERAND) {if (scalar) *C OPERAND (*A) * (*B);\
                                            else {register int i,j,k; register DOUBLE s;\
                                                  register DOUBLE *A_=A, *B_=B, *C_=C, *A_2, *B_2;\
                                                  for (i=blockN; i>0; i--) {\
                                                    for (j=blockN; j>0; j--) {\
                                                      A_2=A_; B_2=B_; s=0.0;\
                                                      for (k=blockN; k>0; k--) {s += (*A_2) * (*B_2); A_2+=blockN; B_2+=blockN;}\
                                                      *C_++ OPERAND s; B_++;\
                                                    } A_++;\
                                                  } } }

#define BLOCK_MUL_TTN_TMPLT(A,B,C,OPERAND) {if (scalar) *C OPERAND (*A) * (*B);\
                                            else {register int i,j,k; register DOUBLE s;\
                                                  register DOUBLE *A_=A, *B_=B, *C_=C, *A_2, *B_2;\
                                                  for (i=blockN; i>0; i--) {\
                                                    for (j=blockN; j>0; j--) {\
                                                      A_2=A_; B_2=B_; s=0.0;\
                                                      for (k=blockN; k>0; k--) {s += (*A_2) * (*B_2++); A_2+=blockN;}\
                                                      *C_++ OPERAND s; A_++;\
                                                    } B_ += blockN; A_ -=blockN;\
                                                  } } }

#define BLOCK_MUL(A,B,C) BLOCK_MUL_NNN_TMPLT(A,B,C,=)
#define BLOCK_MUL_NNN(A,B,C) BLOCK_MUL(A,B,C)
#define BLOCK_MUL_TNN(A,B,C) BLOCK_MUL_TNN_TMPLT(A,B,C,=)
#define BLOCK_MUL_NTN(A,B,C) BLOCK_MUL_NTN_TMPLT(A,B,C,=)
#define BLOCK_MUL_TTN(A,B,C) BLOCK_MUL_TTN_TMPLT(A,B,C,=)
#define BLOCK_MUL_NNT(A,B,C) BLOCK_MUL_TTN(B,A,C)
#define BLOCK_MUL_TNT(A,B,C) BLOCK_MUL_TNN(B,A,C)
#define BLOCK_MUL_NTT(A,B,C) BLOCK_MUL_NTN(B,A,C)
#define BLOCK_MUL_TTT(A,B,C) BLOCK_MUL_NNN(B,A,C)

#define BLOCK_MUL_ADD(A,B,C) BLOCK_MUL_NNN_TMPLT(A,B,C,+=)
#define BLOCK_MUL_ADD_NNN(A,B,C) BLOCK_MUL_ADD(A,B,C)
#define BLOCK_MUL_ADD_TNN(A,B,C) BLOCK_MUL_TNN_TMPLT(A,B,C,+=)
#define BLOCK_MUL_ADD_NTN(A,B,C) BLOCK_MUL_NTN_TMPLT(A,B,C,+=)
#define BLOCK_MUL_ADD_TTN(A,B,C) BLOCK_MUL_TTN_TMPLT(A,B,C,+=)
#define BLOCK_MUL_ADD_NNT(A,B,C) BLOCK_MUL_ADD_TTN(B,A,C)
#define BLOCK_MUL_ADD_TNT(A,B,C) BLOCK_MUL_ADD_TNN(B,A,C)
#define BLOCK_MUL_ADD_NTT(A,B,C) BLOCK_MUL_ADD_NTN(B,A,C)
#define BLOCK_MUL_ADD_TTT(A,B,C) BLOCK_MUL_ADD_NNN(B,A,C)

#define BLOCK_INVERT(A,C) {register DOUBLE *A_=A,*C_=C;\
                           if (scalar) {if ((error=(*A_==0.0))==0) *C_=1.0/(*A_);}\
                           else if (blockN==2) {\
                             register DOUBLE det; det=A_[0]*A_[3]-A_[1]*A_[2];\
                             if ((error=(det==0.0))==0)\
                             {det=1.0/det; *C_++=A_[3]*det; *C_++=-A_[1]*det;\
                              *C_++=-A_[2]*det; *C_=A_[0]*det;}\
                           } else error=InvertFullMatrix_piv(blockN,A_,C_);}

#define BLOCK_TRANSPOSE(A,B) {register DOUBLE *A_=A,*B_=B;\
                              if (scalar) *B_=*A_;\
                              else {register int i,j; register DOUBLE s;\
                                    for (i=blockN; i>0; i--) for (j=i+1; j>0; j--)\
                                      {s=A_[i*N+j]; B_[i*N+j]=A_[j*N+i]; B_[j*N+i]=s;}} }

#define BLOCK_TRANSPOSE1(A) {if (!scalar)\
                             {register int i,j; register DOUBLE s; register DOUBLE *A_=A;\
                              for (i=blockN; i>0; i--) for (j=i+1; j>0; j--)\
                                {s=A_[i*N+j]; A_[i*N+j]=A_[j*N+i]; A_[j*N+i]=s;}} }

#define BLOCK_SKIP_N(skip,A) {register INT skip_=skip; register DOUBLE *A_=A;\
                              if (scalar) {if (skip_&1) *A_=0.0;}\
                              else {register int i,j;\
                                    for (i=blockN; i>0; i--) {\
                                      if (skip_&1) {for (j=blockN; j>0; j--) *A_++=0.0;} else A_+=blockN;\
                                      skip_=skip_>>1;} } }
#define BLOCK_SKIP(skip,A) BLOCK_SKIP_N(skip,A)

#define BLOCK_SKIP_T(skip,A) {register INT skip_=skip; register DOUBLE *A_=A;\
                              if (scalar) {if (skip_&1) *A_=0.0;}\
                              else {register int i,j;\
                                    for (i=blockN; i>0; i--) {\
                                      if (skip_&1) {\
                                        for (j=blockN; j>0; j--) {*A_=0.0; A_+=blockN;}\
                                        A_ -= blockNN;\
                                      } skip_=skip_>>1; A_++;} } }

/*ih*/
#define BLOCK_VECCLEAR(a) {if (scalar) *a=0.0;\
                           else {register int i; register DOUBLE *a_=a;\
                                 for (i=blockN; i>0; i--) *a_++ = 0.0;}}

#define BLOCK_VECCOPY(a,b) {if (scalar) *b = *a;\
                            else {register int i; register DOUBLE *a_=a, *b_=b;\
                                  for (i=blockN; i>0; i--) *b_++ = *a_++;}}

#define BLOCK_VECSCALE(s,a,b) {if (scalar) *b = *a * s;\
                               else {register int i; register DOUBLE s_=s,*a_=a,*b_=b;\
                                     for (i=blockN; i>0; i--) *b_++ = *a_++ * s_;}}

#define BLOCK_VECSCALE1(s,a) {if (scalar) *a *= s;\
                              else {register int i; register DOUBLE s_=s,*a_=a;\
                                    for (i=blockN; i>0; i--) *a_++ *= s_;}}


#define BLOCK_VECADD1(a,b) {if (scalar) *b += *a;\
                            else {register int i; register DOUBLE *a_=a, *b_=b;\
                                  for (i=blockN; i>0; i--) *b_++ += *a_++;}}

#define BLOCK_VECSUB1(a,b) {if (scalar) *b -= *a;\
                            else {register int i; register DOUBLE *a_=a, *b_=b;\
                                  for (i=blockN; i>0; i--) *b_++ -= *a_++;}}

#define BLOCK_MATVEC(A,b,c) {if (scalar) *c = (*A) * (*b);\
                             else {register int i,j; register DOUBLE s, *A_=A, *b_=b, *c_=c, *b_2;\
                                   for (i=blockN; i>0; i--) {\
                                     b_2=b_; s=0.0;\
                                     for (j=blockN; j>0; j--) {s += (*A_++) * (*b_2++);}\
                                     *c_++ = s;};\
                             } }

#define BLOCK_MATVECADD(A,b,c) {if (scalar) *c += (*A) * (*b);\
                                else {register int i,j; register DOUBLE s, *A_=A, *b_=b, *c_=c, *b_2;\
                                      for (i=blockN; i>0; i--) {\
                                        b_2=b_; s=0.0;\
                                        for (j=blockN; j>0; j--) {s += (*A_++) * (*b_2++);}\
                                        *c_++ += s;};\
                                } }

#define BLOCK_WRITEOUT(A) {if (scalar) UserWriteF("A = %g\n",*A);\
                           else {register int i,j,m=0; register DOUBLE *A_=A;\
                                 for (i=blockN; i>0; i--) {\
                                   for (j=blockN; j>0; j--) UserWriteF("A[%d] = %g\n",m++,*A_++);\
                                   UserWrite("\n");\
                                 } } }

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

INT UnmarkAll             (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp);
INT MarkAll               (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp);
INT MarkOffDiagWithoutDirichlet (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp);
INT MarkAbsolute          (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp);
INT MarkRelative          (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp);
INT MarkVanek             (GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp);
INT SetupInitialList      (GRID *theGrid, HEAP *theHeap, AVECTOR **initialSH, AVECTOR **initialEH, INT MarkKey);
INT DistributeInitialList (AVECTOR **La, AVECTOR **Le, AVECTOR **Ta, AVECTOR **Te, AVECTOR **Ua, AVECTOR **Ue);
INT CountStrongNeighbors  (AVECTOR *initialS, DOUBLE *avNrOfStrongNbsHnd, INT *maxNeighbors);
INT GeometricCoarsening   (GRID *theGrid);
INT CoarsenGreedy         (GRID *theGrid);
INT CoarsenGreedyWithBndLoop(GRID *theGrid);
INT CoarsenBreadthFirst   (GRID *theGrid);
INT CoarsenRugeStueben    (GRID *theGrid);
INT CoarsenVanek          (GRID *theGrid);
INT CoarsenAverage        (GRID *theGrid);
INT IpAverage             (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT IpRugeStueben         (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT IpReusken             (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT IpWagner              (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT IpPiecewiseConstant   (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT IpVanek               (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I);
INT FastGalerkinFromInterpolation(GRID *theGrid, MATDATA_DESC *A,
                                  MATDATA_DESC *I, INT type);
INT AssembleGalerkinFromInterpolation(GRID *theGrid, MATDATA_DESC *A,
                                      MATDATA_DESC *I, INT symmetric);
INT SparsenCGMatrix       (GRID *theGrid, MATDATA_DESC *A, INT lumpFlag);
INT ReorderFineGrid       (GRID *theGrid, INT orderType);

INT NBTransformDefect     (GRID *theGrid, const VECDATA_DESC *to,
                           const VECDATA_DESC *from,
                           const MATDATA_DESC *Mat);
INT NBFineGridCorrection  (GRID *theGrid, const VECDATA_DESC *to,
                           const VECDATA_DESC *from,
                           const MATDATA_DESC *Mat);
#endif
