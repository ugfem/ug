// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugblas.h														*/
/*																			*/
/* Purpose:   basic linear algebra routines                                                             */
/*			  working on the matrix-vector and								*/
/*			  matrix-blockvector structure									*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*																			*/
/*			  blockvector routines from:									*/
/*			  Christian Wrobel                                                                              */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*																			*/
/*			  email: ug@ica3.uni-stuttgart.de					                */
/*																			*/
/* History:   06.03.95 begin, ug version 3.0								*/
/*			  28.09.95 blockvector routines implemented (Christian Wrobel)	*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UGBLAS__
#define __UGBLAS__

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

#define UPPER_TRIANGLE          1
#define LOWER_TRIANGLE          -1

/* kinds of matrices														*/
#define R1C1                            RCKIND(1,1)
#define R1C2                            RCKIND(1,2)
#define R1C3                            RCKIND(1,3)
#define R2C1                            RCKIND(2,1)
#define R2C2                            RCKIND(2,2)
#define R2C3                            RCKIND(2,3)
#define R3C1                            RCKIND(3,1)
#define R3C2                            RCKIND(3,2)
#define R3C3                            RCKIND(3,3)

/* calc kind of matrix from nrows and ncols (nr,nc <=3)                                         */
#define RCKIND(nr,nc)           ((nr) << 4 | (nc))
#ifdef OLD_NUM_COMPATIBLE
#define MAT_RCKIND(M,rt,ct)     (((rt)<4 && (ct)<4) ? RCKIND(MDT_NROWCOMP(M,rt,ct),MDT_NCOLCOMP(M,rt,ct)) : -1)
#else
#define MAT_RCKIND(M,rt,ct)     (((rt)<4 && (ct)<4) ? RCKIND(MD_ROWS_IN_RT_CT(M,rt,ct),MD_COLS_IN_RT_CT(M,rt,ct)) : -1)
#endif

/* the vector loops used													*/
#define L_VLOOP__CLASS(v,first_v,c)                                                                                     \
  for (v=first_v; v!= NULL; v=SUCCVC(v))                          \
    if (VCLASS(v)>=c)

#define L_VLOOP__TYPE_CLASS(v,first_v,t,c)                                                                      \
  for (v=first_v; v!= NULL; v=SUCCVC(v))                          \
    if ((VTYPE(v)==t) && (VCLASS(v)>=c))

#define L_VLOOP__TYPE_CLASS2(v,first_v,end_v,t,c)                                                                       \
  for (v=first_v; v!= end_v; v=SUCCVC(v))                         \
    if ((VTYPE(v)==t) && (VCLASS(v)>=c))

#define L_REVERSE_VLOOP__TYPE_CLASS(v,last_v,t,c)                                                       \
  for (v=last_v; v!= NULL; v=PREDVC(v))                           \
    if ((VTYPE(v)==t) && (VCLASS(v)>=c))

#define L_REVERSE_VLOOP__CLASS(v,last_v,c)                                                                      \
  for (v=last_v; v!= NULL; v=PREDVC(v))                           \
    if (VCLASS(v)>=c)

#define L_REVERSE_VLOOP__CLASS(v,last_v,c)                                                              \
  for (v=last_v; v!= NULL; v=PREDVC(v))                       \
    if (VCLASS(v)>=c)

#define A_VLOOP__TYPE_CLASS(l,fl,tl,v,mg,t,c)                                                           \
  for (l=fl; l<=tl; l++)                                                          \
    for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,l)); v!= NULL; v=SUCCVC(v)) \
      if ((VTYPE(v)==t) && (VCLASS(v)>=c))

#define S_BELOW_VLOOP__TYPE(l,fl,tl,v,mg,t)                                                                     \
  for (l=fl; l<tl; l++)                                                           \
    for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,l)); v!= NULL; v=SUCCVC(v)) \
      if ((VTYPE(v)==t) && (FINE_GRID_DOF(v)))

#define S_FINE_VLOOP__TYPE(tl,v,mg,t)                                                                           \
  for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,tl)); v!= NULL; v=SUCCVC(v)) \
    if ((VTYPE(v)==t) && (NEW_DEFECT(v)))


/* the matrix loops used													*/
#define L_MLOOP__RCTYPE(v,first_v,m,rt,ct)                                                                      \
  for (v=first_v; v!= NULL; v=SUCCVC(v))          \
    if (VTYPE(v)==rt)                                                               \
      for (m=VSTART(v); m!=NULL; m=MNEXT(m))          \
        if (VTYPE(MDEST(m))==ct)

#define L_MLOOP__RCTYPE2(v,first_v,end_v,m,rt,ct)                                                                       \
  for (v=first_v; v!= end_v; v=SUCCVC(v))         \
    if (VTYPE(v)==rt)                                                               \
      for (m=VSTART(v); m!=NULL; m=MNEXT(m))          \
        if (VTYPE(MDEST(m))==ct)

#define S_BELOW_MLOOP__RCTYPE(l,fl,tl,v,mg,m,rt,ct)                                                     \
  for (l=fl; l<tl; l++)                                                           \
    for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,l)); v!= NULL; v=SUCCVC(v)) \
      if ((VTYPE(v)==rt) && (FINE_GRID_DOF(v)))       \
        for (m=VSTART(v); m!=NULL; m=MNEXT(m))  \
          if (VTYPE(MDEST(m))==ct)

#define S_FINE_MLOOP__RCTYPE(tl,v,mg,m,rt,ct)                                                           \
  for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,tl)); v!= NULL; v=SUCCVC(v)) \
    if ((VTYPE(v)==rt) && (NEW_DEFECT(v)))          \
      for (m=VSTART(v); m!=NULL; m=MNEXT(m))  \
        if (VTYPE(MDEST(m))==ct)


/* macros for matrix operations                                                                                         */

#define MATMUL_11(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0);}

#define MATMUL_12(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 01) * VVALUE(vec,c ## 1);}

#define MATMUL_13(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 01) * VVALUE(vec,c ## 1)       \
             + MVALUE(mat,m ## 02) * VVALUE(vec,c ## 2);}

#define MATMUL_21(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0);  \
   s ## 1 += MVALUE(mat,m ## 10) * VVALUE(vec,c ## 0);}

#define MATMUL_22(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 01) * VVALUE(vec,c ## 1);      \
   s ## 1 += MVALUE(mat,m ## 10) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 11) * VVALUE(vec,c ## 1);}

#define MATMUL_23(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 01) * VVALUE(vec,c ## 1)       \
             + MVALUE(mat,m ## 02) * VVALUE(vec,c ## 2);      \
   s ## 1 += MVALUE(mat,m ## 10) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 11) * VVALUE(vec,c ## 1)       \
             + MVALUE(mat,m ## 12) * VVALUE(vec,c ## 2);}

#define MATMUL_31(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0);  \
   s ## 1 += MVALUE(mat,m ## 10) * VVALUE(vec,c ## 0);  \
   s ## 2 += MVALUE(mat,m ## 20) * VVALUE(vec,c ## 0);}

#define MATMUL_32(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 01) * VVALUE(vec,c ## 1);      \
   s ## 1 += MVALUE(mat,m ## 10) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 11) * VVALUE(vec,c ## 1);      \
   s ## 2 += MVALUE(mat,m ## 20) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 21) * VVALUE(vec,c ## 1);}

#define MATMUL_33(s,mat,m,vec,c)                                                                                        \
  {s ## 0 += MVALUE(mat,m ## 00) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 01) * VVALUE(vec,c ## 1)       \
             + MVALUE(mat,m ## 02) * VVALUE(vec,c ## 2);      \
   s ## 1 += MVALUE(mat,m ## 10) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 11) * VVALUE(vec,c ## 1)       \
             + MVALUE(mat,m ## 12) * VVALUE(vec,c ## 2);      \
   s ## 2 += MVALUE(mat,m ## 20) * VVALUE(vec,c ## 0)   \
             + MVALUE(mat,m ## 21) * VVALUE(vec,c ## 1)       \
             + MVALUE(mat,m ## 22) * VVALUE(vec,c ## 2);}

/* the blockvector loops used												*/
#define BLOCK_L_VLOOP(v,first_v,end_v)  \
  for (v=first_v; v!= end_v; v=SUCCVC(v))

#define BLOCK_L_VLOOP__TYPE_CLASS(v,first_v,end_v,t,c)  \
  for (v=first_v; v!= end_v; v=SUCCVC(v))                         \
    if ((VTYPE(v)==t) && (VCLASS(v)>=c))


/* the matrix-blockvector loops used										*/
#define BLOCK_L_MLOOP(v,first_v,end_v,bvd_col,bvdf,m)   \
  for (v=first_v; v!= end_v; v=SUCCVC(v))                         \
    for (m=VSTART(v); m!=NULL; m=MNEXT(m))          \
      if (VMATCH(MDEST(m),bvd_col,bvdf))

#define BLOCK_L_MLOOP__RCTYPE(v,first_v,end_v,bvd_col,bvdf,m,rt,ct)     \
  for (v=first_v; v!= end_v; v=SUCCVC(v))                         \
    if (VTYPE(v)==rt)                                                               \
      for (m=VSTART(v); m!=NULL; m=MNEXT(m))          \
        if (VMATCH(MDEST(m),bvd_col,bvdf)&&(VTYPE(MDEST(m))==ct))


/* ptr to begin of vector values (not to be confused with VVALUEPTR(v,n))	*/
#define VVALPTR(v)                              ((v)->value)

/* ptr to begin of matrix values (not to be confused with MVALUEPTR(m,n))	*/
#define MVALPTR(m)                              ((m)->value)

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT  VecCheckConsistency                        (const VECDATA_DESC *x, const VECDATA_DESC *y);
INT  MatmulCheckConsistency             (const VECDATA_DESC *x, const MATDATA_DESC *M, const VECDATA_DESC *y);

#endif
