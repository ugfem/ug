// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugiter.c														*/
/*																			*/
/* Purpose:   iterative schemes and decompositions							*/
/*			  working on the matrix-vector structure						*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   25.03.95 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <math.h>

#include "algebra.h"
#include "compiler.h"
#include "disctools.h"
#include "debug.h"
#include "devices.h"
#include "general.h"
#include "gm.h"
#include "misc.h"
#include "pargm.h"

#include "np.h"
#include "ugblas.h"
#include "block.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define SMALL_DET                       1e-15
#define MAX_DEPTH           8
#define V_BVNUMBER(v,n)         (VINDEX(v)/n)

/* macros to define VEC_SCALAR, VECDATA_DESC and MATDATA_DESC components */
#define DEFINE_VS_CMPS(a)                               register DOUBLE a ## 0,a ## 1,a ## 2
#define DEFINE_VD_CMPS(x)                               register INT x ## 0,x ## 1,x ## 2
#define DEFINE_MD_CMPS(m)                               register INT m ## 00,m ## 01,m ## 02,m ## 10,m ## 11,m ## 12,m ## 20,m ## 21,m ## 22

/* macros to set VECDATA_DESC and MATDATA_DESC components */
#define SET_YCMP_1(y,v,tp,cp)                   {y ## 0 = (VD_CMPPTR_OF_TYPE(v,tp))[0];}
#define SET_YCMP_2(y,v,tp,cp)                   {cp=VD_CMPPTR_OF_TYPE(v,tp); y ## 0 = (cp)[0]; y ## 1 = (cp)[1];}
#define SET_YCMP_3(y,v,tp,cp)                   {cp=VD_CMPPTR_OF_TYPE(v,tp); y ## 0 = (cp)[0]; y ## 1 = (cp)[1]; y ## 2 = (cp)[2];}

#define SET_MCMP_11(m,M,rt,ct,cp)               {m ## 00 = MD_MCMPPTR_OF_RT_CT(M,rt,ct)[0];}
#define SET_MCMP_12(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; m ## 01 = (cp)[1];}
#define SET_MCMP_13(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; m ## 01 = (cp)[1]; m ## 02 = (cp)[2];}
#define SET_MCMP_21(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; \
                                                 m ## 10 = (cp)[1];}
#define SET_MCMP_22(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; m ## 01 = (cp)[1]; \
                                                 m ## 10 = (cp)[2]; m ## 11 = (cp)[3];}
#define SET_MCMP_23(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; m ## 01 = (cp)[1]; m ## 02 = (cp)[2]; \
                                                 m ## 10 = (cp)[3]; m ## 11 = (cp)[4]; m ## 12 = (cp)[5];}
#define SET_MCMP_31(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; \
                                                 m ## 10 = (cp)[1]; \
                                                 m ## 20 = (cp)[2];}
#define SET_MCMP_32(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; m ## 01 = (cp)[1]; \
                                                 m ## 10 = (cp)[2]; m ## 11 = (cp)[3]; \
                                                 m ## 20 = (cp)[4]; m ## 21 = (cp)[5];}
#define SET_MCMP_33(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; m ## 01 = (cp)[1]; m ## 02 = (cp)[2]; \
                                                 m ## 10 = (cp)[3]; m ## 11 = (cp)[4]; m ## 12 = (cp)[5]; \
                                                 m ## 20 = (cp)[6]; m ## 21 = (cp)[7]; m ## 22 = (cp)[8];}

#define SET_CMPS_11(y,v,m,M,rt,ct,cp)   SET_MCMP_11(m,M,rt,ct,cp); SET_YCMP_1(y,v,ct,cp);
#define SET_CMPS_12(y,v,m,M,rt,ct,cp)   SET_MCMP_12(m,M,rt,ct,cp); SET_YCMP_2(y,v,ct,cp);
#define SET_CMPS_13(y,v,m,M,rt,ct,cp)   SET_MCMP_13(m,M,rt,ct,cp); SET_YCMP_3(y,v,ct,cp);
#define SET_CMPS_21(y,v,m,M,rt,ct,cp)   SET_MCMP_21(m,M,rt,ct,cp); SET_YCMP_1(y,v,ct,cp);
#define SET_CMPS_22(y,v,m,M,rt,ct,cp)   SET_MCMP_22(m,M,rt,ct,cp); SET_YCMP_2(y,v,ct,cp);
#define SET_CMPS_23(y,v,m,M,rt,ct,cp)   SET_MCMP_23(m,M,rt,ct,cp); SET_YCMP_3(y,v,ct,cp);
#define SET_CMPS_31(y,v,m,M,rt,ct,cp)   SET_MCMP_31(m,M,rt,ct,cp); SET_YCMP_1(y,v,ct,cp);
#define SET_CMPS_32(y,v,m,M,rt,ct,cp)   SET_MCMP_32(m,M,rt,ct,cp); SET_YCMP_2(y,v,ct,cp);
#define SET_CMPS_33(y,v,m,M,rt,ct,cp)   SET_MCMP_33(m,M,rt,ct,cp); SET_YCMP_3(y,v,ct,cp);

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT StoreInverse=TRUE;

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* some functions to assign values to the global variables of this file		*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  IterCheckConsistencyEtc										*/
/*																			*/
/* Purpose:   check wether the VECDATA_DESC is consistent with the format	*/
/*			  and wether DESCRIPTORS match									*/
/*			  and wether the problem is scalar								*/
/*																			*/
/* Input:	  FORMAT		*theFormat										*/
/*			  VECDATA_DESC *theVec											*/
/*            MATDATA_DESC *M                   the matrix to precondition          */
/*            VECDATA_DESC *d           right hand side (the defect)            */
/*                                                                          */
/*            INT			*isscal     treat as scalar if TRUE             */
/*            INT			*vc         v component				            */
/*            INT			*mc         M component				            */
/*            INT			*dc         d component				            */
/*            INT			*mask       data type maks			            */
/*                                                                          */
/* Output:    INT 0: ok                                                         */
/*            INT 1: error	                                                */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   l_setindex - set index field in vectors of the specified types

   SYNOPSIS:
   INT l_setindex (GRID *g);

   PARAMETERS:
   .  g - pointer to grid

   DESCRIPTION:
   This function sets index field in vectors of the specified types.

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

INT l_setindex (GRID *g)
{
  VECTOR *v,*first_v;
  INT i;

  first_v = FIRSTVECTOR(g);

  /* enumerate VECTORS starting with 1 (so return value if decomp failed will be negative) */
  i = 1;
  L_VLOOP__CLASS(v,first_v,EVERY_CLASS)
  VINDEX(v) = i++;

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_ordervtypes - Rearrange doubly linked vector list according to TypeOrder

   SYNOPSIS:
   INT l_ordervtypes (GRID *g, const SHORT TypeOrder[NVECTYPES]);

   PARAMETERS:
   .  g - pointer to grid
   .  TypeOrder[NVECTYPES] - describes the order

   DESCRIPTION:
   This function rearranges doubly linked vector list according to TypeOrder.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if an error occured.
   D*/
/****************************************************************************/

INT l_ordervtypes (GRID *g, const SHORT TypeOrder[NVECTYPES])
{
  VECTOR *v,*first_v,*last_v_of_type[NVECTYPES];
  INT type,order,Types[NVECTYPES];

  /* consistency check before whole list is messed up */
  for (type=0; type<NVECTYPES; type++) Types[type] = FALSE;
  for (order=0; order<NVECTYPES; order++) Types[TypeOrder[order]] = TRUE;
  for (type=0; type<NVECTYPES; type++)
    if (!Types[type])
      REP_ERR_RETURN (NUM_ERROR);

  for (type=0; type<NVECTYPES; type++) last_v_of_type[type] = NULL;

  first_v = FIRSTVECTOR(g);

  /* run over succ list and sort in NVECTYPES lists using the pred ptrs */
  for (v=first_v; v!=NULL; v=SUCCVC(v))
  {
    type = VTYPE(v);
    PREDVC(v) = last_v_of_type[type];
    last_v_of_type[type] = v;
  }

  /* rebuild the complete list according to TypeOrder */
  SFIRSTVECTOR(g) = NULL;
  for (order=NVECTYPES-1; order>=0; order--)
    for (v=last_v_of_type[TypeOrder[order]]; v!=NULL; v=PREDVC(v))
    {
      SUCCVC(v) = FIRSTVECTOR(g);
      SFIRSTVECTOR(g) = v;
      if (SUCCVC(v)!=NULL)
        PREDVC(SUCCVC(v)) = v;
    }

  for (order=NVECTYPES-1; order>=0; order--)
    if (last_v_of_type[TypeOrder[order]]!=NULL)
    {
      LASTVECTOR(g)  = last_v_of_type[TypeOrder[order]];
      break;
    }

  PREDVC(FIRSTVECTOR(g)) = NULL;
  SUCCVC(LASTVECTOR(g))  = NULL;

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_jac - solve Diag(M)v=d

   SYNOPSIS:
   INT l_jac (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
   const VECDATA_DESC *d);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `Diag(M) v = d`, where only the diagonal of 'M'
   is used.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    NUM_BLOCK_TOO_LARGE if the blocks are larger as MAX_SINGLE_VEC_COMP
   .n    NUM_SMALL_DIAG if a diagonal block system is (nearly) singular.
   D*/
/****************************************************************************/

INT l_jac (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*first_vec;
  INT rtype,err;
  register SHORT vc,dc,mc,mask;
  register SHORT *vcomp,*dcomp;
  register SHORT i;
  register SHORT n;
  DOUBLE s[MAX_SINGLE_VEC_COMP];

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = FIRSTVECTOR(g);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve Diag(M)v=d */
    for (vec=first_vec; vec!= NULL; vec=SUCCVC(vec))
      if ((VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
        VVALUE(vec,vc) = VVALUE(vec,dc)/MVALUE(VSTART(vec),mc);

    return (NUM_OK);
  }

  for (rtype=0; rtype<NVECTYPES; rtype++)
    if (VD_NCMPS_IN_TYPE(v,rtype)>0)
    {
      n     = VD_NCMPS_IN_TYPE(v,rtype);
      vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
      dcomp = VD_CMPPTR_OF_TYPE(d,rtype);
      L_VLOOP__TYPE_CLASS(vec,first_vec,rtype,ACTIVE_CLASS)
      {
        /* rhs */
        for (i=0; i<n; i++)
          s[i] = VVALUE(vec,dcomp[i]);

        /* solve */
        if (SolveSmallBlock(n,vcomp,VVALPTR(vec),
                            MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                            MVALPTR(VSTART(vec)),s)!=0)
          REP_ERR_RETURN (NUM_SMALL_DIAG);
      }
    }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   jacBS - iterates Diag(K)u=f for a blockmatrix K (Jacobi method)

   SYNOPSIS:
   INT jacBS ( const BLOCKVECTOR *bv, const BV_DESC *bvd,
            const BV_DESC_FORMAT *bvdf, INT K_comp, INT u_comp, INT f_comp );

   PARAMETERS:
   .  bv - blockvector of the matrix
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  K_comp - position of the scalar in the MATRIXs of the blockmatrix
   .  u_comp - position of the solution scalar in the VECTORs of the blockvector
   .  f_comp - position of the right hand side scalar in the VECTORs of the blockvector

   DESCRIPTION:
   This function solves Ku = f by solving `Diag(K) u = f` (Jacobi method),
   where only the diagonal of 'K' is used. The blockmatrix K must be on the
   diagonal of the global matrix.

   u_comp and f_comp must be different!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   D*/
/****************************************************************************/

INT jacBS ( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT K_comp, INT u_comp, INT f_comp )
{
  VECTOR *v, *end_v, *first_v;

  ASSERT( (u_comp >= 0) && (K_comp >= 0) && (f_comp >= 0) );
  ASSERT( u_comp != f_comp );

  first_v = BVFIRSTVECTOR( bv );
  end_v   = BVENDVECTOR( bv );

  /* solve Diag(K)u=d */
  BLOCK_L_VLOOP(v,first_v,end_v)
  VVALUE(v,u_comp) = VVALUE(v,f_comp)/MVALUE(VSTART(v),K_comp);

  REP_ERR_RETURN (NUM_OK);
}


/****************************************************************************/
/*D
   l_lgs - solve LowerTriangle

   SYNOPSIS:
   INT l_lgs (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
   const VECDATA_DESC *d);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `LowerTriangle(M) v = d`,
   where only the lower part of 'M' is used.

   `Remark.` Index field must be set!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_lgs (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec;
  INT rtype,ctype,myindex,err;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = FIRSTVECTOR(g);

  if (MD_IS_SCALAR(M) && VD_IS_SCALAR(v) && VD_IS_SCALAR(d))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= NULL; vec=SUCCVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum)/MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  L_VLOOP__CLASS(vec,first_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n     = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    dcomp = VD_CMPPTR_OF_TYPE(d,rtype);
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++)
      s[i] = VVALUE(vec,dcomp[i]);
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          wmat  = VVALPTR(w);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
        }

    /* solve */
    if (SolveSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),VVALPTR(vec),
                        MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                        MVALPTR(VSTART(vec)),s)!=0)
      REP_ERR_RETURN (__LINE__);
  }

  return (NUM_OK);
}

INT l_lgs_SB (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*end_vec;
  INT myindex,err,first_index;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register DOUBLE sum;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = BVFIRSTVECTOR(theBV);
  first_index = VINDEX(first_vec);
  end_vec = BVENDVECTOR(theBV);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= end_vec; vec=SUCCVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) && (VINDEX(w)>=first_index))
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum)/MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  REP_ERR_RETURN (__LINE__);
}

INT l_tplgs_SB (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*end_vec;
  INT myindex,err,last_index;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register DOUBLE sum;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = BVLASTVECTOR(theBV);
  last_index = VINDEX(first_vec);
  end_vec = PREDVC(BVFIRSTVECTOR(theBV));

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= end_vec; vec=PREDVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) && (VINDEX(w)<=last_index))
            sum += MVALUE(MADJ(mat),mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum)/MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  REP_ERR_RETURN (__LINE__);
}

/****************************************************************************/
/*D
   l_ugs - solve UpperTriangle

   SYNOPSIS:
   INT l_ugs (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
   const VECDATA_DESC *d);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `UpperTriangle(M) v = d`,
   where only the upper part of 'M' is used.

   `Remark.` Index field must be set!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_ugs (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*last_vec;
  INT rtype,ctype,myindex,err;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  last_vec = LASTVECTOR(g);

  if (MD_IS_SCALAR(M) && VD_IS_SCALAR(v) && VD_IS_SCALAR(d))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve UpperTriangle(v)=d */
    for (vec=last_vec; vec!= NULL; vec=PREDVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum)/MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  L_REVERSE_VLOOP__CLASS(vec,last_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n               = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    dcomp   = VD_CMPPTR_OF_TYPE(d,rtype);
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,dcomp[i]);
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          wmat  = VVALPTR(w);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
        }

    /* solve */
    if (SolveSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),VVALPTR(vec),
                        MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                        MVALPTR(VSTART(vec)),s)!=0)
      REP_ERR_RETURN (__LINE__);
  }

  return (NUM_OK);
}

INT l_ugs_SB (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*last_vec,*end_vec;
  INT myindex,err,last_index;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register DOUBLE sum;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  last_vec = BVLASTVECTOR(theBV);
  end_vec = PREDVC(BVFIRSTVECTOR(theBV));
  last_index = VINDEX(last_vec);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve UpperTriangle(v)=d */
    for (vec=last_vec; vec!= end_vec; vec=PREDVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) && (VINDEX(w)<=last_index))
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum)/MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  REP_ERR_RETURN (__LINE__);
}

/****************************************************************************/
/*D
   l_lgsB - solve LowerBLOCKVECTORTriangle

   SYNOPSIS:
   INT l_lgsB (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const MATDATA_DESC *LU,
   const VECDATA_DESC *d);

   PARAMETERS:
   .  g  - pointer to grid
   .  v  - type vector descriptor to store correction
   .  M  - type matrix descriptor for precondition
   .  d  - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `LowerBLOCKVECTORTriangle(M) v = d`,
   where only the upper part of 'M' is used.

   `Remark.` Index field must be set!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_lgsB (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w;
  BLOCKVECTOR *theBV;
  INT rtype,ctype,myindex,err,maxBVmembers,bvn,index;
  register MATRIX *mat;
  register SHORT vc,mc,mask;
  register SHORT *mcomp,*wcomp,*vcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  maxBVmembers = NVEC(g);

  /* set index field again, because other prepares may have called l_setindex,
     other iterative schemes will still work since they only need the index field
     in ascending order */
  for (theBV=GFIRSTBV(g); theBV!=NULL; theBV=BVSUCC(theBV))
  {
    index = BVNUMBER(theBV)*maxBVmembers;
    for (vec=BVFIRSTVECTOR(theBV); vec!=BVENDVECTOR(theBV); vec=SUCCVC(vec))
      VINDEX(vec) = index++;
  }

  /* first copy rhs to v, then solve for v in place */
  if (l_dcopy(g,v,ACTIVE_CLASS,d)!=0)
    REP_ERR_RETURN (__LINE__);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerBLOCKVECTORTriangle(v)=d */
    for (theBV=GFIRSTBV(g); theBV!=NULL; theBV=BVSUCC(theBV))
    {
      bvn = BVNUMBER(theBV);
      /* lower triangle to rhs */
      for (vec=BVFIRSTVECTOR(theBV); vec!=BVENDVECTOR(theBV); vec=SUCCVC(vec))
      {
        if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
        {
          sum = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
          {
            w = MDEST(mat);
            if ((V_BVNUMBER(w,maxBVmembers)<bvn) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
              sum += MVALUE(mat,mc)*VVALUE(w,vc);
          }
          VVALUE(vec,vc) -= sum;
        }
      }
      /* solve */
      if (l_luiterB(g,theBV,v,M,v)!=0)
        REP_ERR_RETURN (-bvn);
    }

    return (NUM_OK);
  }

  for (theBV=GFIRSTBV(g); theBV!=NULL; theBV=BVSUCC(theBV))
  {
    bvn = BVNUMBER(theBV);
    /* lower triangle to rhs */
    for (vec=BVFIRSTVECTOR(theBV); vec!=BVENDVECTOR(theBV); vec=SUCCVC(vec))
    {
      if (VCLASS(vec)<ACTIVE_CLASS) continue;

      rtype = VTYPE(vec);

      n     = VD_NCMPS_IN_TYPE(v,rtype);
      if (n == 0) continue;
      myindex = VINDEX(vec);

      /* rhs */
      for (i=0; i<n; i++)
        s[i] = 0.0;
      for (ctype=0; ctype<NVECTYPES; ctype++)
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
          switch (MAT_RCKIND(M,rtype,ctype))
          {
          case R1C1 :
            SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_11(s,mat,m,w,cy)
                s[0] += s0;
            break;

          case R1C2 :
            SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_12(s,mat,m,w,cy)
                s[0] += s0;
            break;

          case R1C3 :
            SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_13(s,mat,m,w,cy)
                s[0] += s0;
            break;

          case R2C1 :
            SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = s1 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_21(s,mat,m,w,cy)
                s[0] += s0;
            s[1] += s1;
            break;

          case R2C2 :
            SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = s1 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_22(s,mat,m,w,cy)
                s[0] += s0;
            s[1] += s1;
            break;

          case R2C3 :
            SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = s1 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_23(s,mat,m,w,cy)
                s[0] += s0;
            s[1] += s1;
            break;

          case R3C1 :
            SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = s1 = s2 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_31(s,mat,m,w,cy)
                s[0] += s0;
            s[1] += s1;
            s[2] += s2;
            break;

          case R3C2 :
            SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = s1 = s2 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_32(s,mat,m,w,cy)
                s[0] += s0;
            s[1] += s1;
            s[2] += s2;
            break;

          case R3C3 :
            SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
            s0 = s1 = s2 = 0.0;
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                MATMUL_33(s,mat,m,w,cy)
                s[0] += s0;
            s[1] += s1;
            s[2] += s2;
            break;

          default :
            mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
            wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
            wmat  = VVALPTR(w);
            nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
            for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
              if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (V_BVNUMBER(w,maxBVmembers)<bvn))
                for (i=0; i<n; i++)
                  for (j=0; j<nc; j++)
                    s[i] += MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
          }
      vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
      for (i=0; i<n; i++)
        VVALUE(vec,vcomp[i]) -= s[i];
    }
    /* solve */
    if (l_luiterB(g,theBV,v,M,v)!=0)
      REP_ERR_RETURN (-bvn);
  }

  return (NUM_OK);
}


#ifdef __BLOCK_VECTOR_DESC__

/****************************************************************************/
/*D
   gs_solveBS - iterates LowerTriangle(K)u=f for a blockmatrix K (Gauss-Seidel method)

   SYNOPSIS:
   INT gs_solveBS ( const BLOCKVECTOR *bv, const BV_DESC *bvd,
   const BV_DESC_FORMAT *bvdf, DOUBLE eps, INT max_it, INT K_comp,
   INT u_comp, INT f_comp, INT aux_comp, INT verbose, INT eps_relative );

   PARAMETERS:
   .  bv - blockvector of the matrix
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  eps - final accuracy for the defect
   .  max_it - maximum of iterations to be performed
   .  K_comp - position of the scalar in the MATRIXs of the blockmatrix
   .  u_comp - position of the solution scalar in the VECTORs of the blockvector
   .  f_comp - position of the right hand side scalar in the VECTORs of the blockvector
   .  aux_comp - if nonnegative denotes auxiliary component for defect calculation
   .  verbose - if TRUE display the convergence rate of the iteration
   .  eps_relative - error reduction relative to the start defect

   DESCRIPTION:
   This function solve Ku = f by iterative solvings of `LowerTriangle(K) v = f`.
   The blockmatrix K must be on the diagonal of the
   global matrix.

   If there is no auxiliary component (dummy value < 0) given the given
   number of Gauss-Seidel-steps are performed. If the auxiliary component is
   given, the defect of the solution can be calculated and the iteration is
   stopped premature if the prescribed accuracy 'eps' is reached; this
   stopping-accuracy will be multiplied by the start defect norm (resulting
   in a relative criterion) if 'eps_relative' is TRUE.

   u_comp, f_comp (and aux_comp if given) must be different!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    1 if not converged
   D*/
/****************************************************************************/

INT gs_solveBS ( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, DOUBLE eps, INT max_it, INT K_comp, INT u_comp, INT f_comp, INT aux_comp, INT verbose, INT eps_relative )
{
  VECTOR *v, *end_v, *first_v, *w;
  register SHORT it;
  register MATRIX *mat;
  register DOUBLE sum;
  DOUBLE start_norm, new_norm;
  INT ret = 0;

  ASSERT( (u_comp >= 0) && (K_comp >= 0) && (f_comp >= 0) );
  ASSERT( u_comp != f_comp );

  first_v = BVFIRSTVECTOR( bv );
  end_v   = BVENDVECTOR( bv );

  if ( aux_comp > -1 )
  {
    ASSERT( u_comp != aux_comp );
    ASSERT( f_comp != aux_comp );
    start_norm = new_norm = CalculateDefectAndNormBS( bv, bvd, bvdf, aux_comp, f_comp, K_comp, u_comp );
    if ( eps_relative )
      eps *= start_norm;
  }
  else
    new_norm = eps + 1.0;               /* avoid premature terminating */

  it = 0;

  while ( ( new_norm > eps ) && ( it < max_it ) )
  {
    it++;

    /* other form of Gauss-Seidel: ui_new=(fi - sum(kij*uj_new){j<i} - sum(kij*uj_old){j>i} )/kii */
    BLOCK_L_VLOOP(v,first_v,end_v)
    {
      sum = 0.0;
      for (mat=MNEXT(VSTART(v)); mat!=NULL; mat = MNEXT(mat))
      /* skip the diagonal element */
      {
        w = MDEST(mat);
        if ( VMATCH(w, bvd, bvdf) )
          sum += MVALUE(mat,K_comp) * VVALUE(w,u_comp);
      }
      VVALUE(v,u_comp) = (VVALUE(v,f_comp) - sum) / MVALUE(VSTART(v),K_comp);
    }

    if ( aux_comp > -1 )
    {
      /* new defect */
      new_norm = CalculateDefectAndNormBS( bv, bvd, bvdf, aux_comp, f_comp, K_comp, u_comp );
    }
  }

  if ( aux_comp > -1 )
  {
    if ( it >= max_it )
    {
      UserWrite( "gauss seidel max. iteration not sufficient++++++++++\n" );
      printf( "gauss seidel max. iteration not sufficient++++++++++\n" );
      ret = 1;
    }

    if ( verbose )
    {
      UserWriteF( "gauss seidel avarage of convergency rate ( %d iterations) = %12lg, end defect = %g\n", it, pow( new_norm / start_norm, 1.0 / (DOUBLE)it ), new_norm );
      printf( "gauss seidel avarage of convergency rate ( %d iterations) = %12lg, end defect = %g\n", it, pow( new_norm / start_norm, 1.0 / (DOUBLE)it ), new_norm );
    }
  }

  if ( ret != NUM_OK )
    REP_ERR_RETURN(ret);

  return (ret);
}
#endif

/****************************************************************************/
/*D
   l_lsor - solve LowerTriangle

   SYNOPSIS:
   INT l_lsor (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
   const VECDATA_DESC *d, const DOUBLE *damp);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `LowerTriangle(M) v = d`,
   where only the lower part of 'M' is used,
   using successive overrelaxation.

   `Remark.` Index field must be set!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_lsor (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d, const DOUBLE *damp)
{
  VECTOR *vec,*w,*first_vec;
  INT rtype,ctype,myindex,err;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*wcomp,*dcomp,*vcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum,dmp;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat;
  const DOUBLE *tdmp;
  const SHORT *offset;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = FIRSTVECTOR(g);
  offset = VD_OFFSETPTR(v);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);
    dmp  = damp[0];

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= NULL; vec=SUCCVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = dmp*(VVALUE(vec,dc)-sum)/MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  L_VLOOP__CLASS(vec,first_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n               = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    vcomp   = VD_CMPPTR_OF_TYPE(v,rtype);
    dcomp   = VD_CMPPTR_OF_TYPE(d,rtype);
    myindex = VINDEX(vec);
    tdmp     = damp+offset[rtype];

    /* rhs */
    for (i=0; i<n; i++)
      s[i] = VVALUE(vec,dcomp[i]);
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy)
              s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy)
              s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          wmat  = VVALPTR(w);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
        }

    /* solve */
    if (SolveSmallBlock(n,vcomp,VVALPTR(vec),
                        MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                        MVALPTR(VSTART(vec)),s)!=0)
      REP_ERR_RETURN (__LINE__);

    /* damp */
    for (i=0; i<n; i++)
      VVALUE(vec,vcomp[i]) *= tdmp[i];
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_ilubthdecomp	- compute incomplete decomposition

   SYNOPSIS:
   INT l_ilubthdecomp (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta,
   const VEC_SCALAR threshold, const VECDATA_DESC *rest,
   const VEC_SCALAR oldrestthresh);

   PARAMETERS:
   .  g - pointer to grid
   .  M - type matrix descriptor
   .  beta - modification parameters
   .  threshold - introduce new connection if entry > threshold
   .  rest - if !=NULL store normalized rest matrix row sums here
   .  oldrestthresh - if !=NULL use rest field of last step to flag vectors where
                                                connections will be introduced

   DESCRIPTION:
   This function computes an incomplete decomposition of order 0
   with modification of the diagonal.
   A new connection is introduced if an entry is larger than threshold.

   The parameters beta and threshold are not considered if they are 'NULL'.

   The matrix M is overwritten!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    i<0 if decomposition failed
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_ilubthdecomp (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VEC_SCALAR threshold, const VECDATA_DESC *VD_rest, const VEC_SCALAR oldrestthresh)
{
  VECTOR *vi,*vj,*vk;
  MATRIX *Mij,*Mji,*Mjk,*Mik;
  DOUBLE InvMat[MAX_SINGLE_MAT_COMP],PivMat[MAX_SINGLE_MAT_COMP];
  DOUBLE CorMat[MAX_SINGLE_MAT_COMP];
  DOUBLE sum;
  INT offset[NVECTYPES+1];
  register DOUBLE *Diag,*Piv,*Elm,*Mat,*Djj,*Dkk;
  register SHORT *DiagComp,*PivComp,*ElmComp,*MatComp,*DjjComp,*DkkComp;
  register INT i0,j0,k0,l,m;
  INT type,ctype,rtype,PivIsZero,CorIsZero;
  INT n,n2,nr,nnr,nc,nrnc;
  INT i,mc,mask;
  DOUBLE RowSum[MAX_SINGLE_VEC_COMP],Damp[MAX_SINGLE_VEC_COMP];
  VEC_SCALAR j_Normalization,k_Normalization;
  DOUBLE diag,invdiag,pivot,AbsDjj;
  const DOUBLE *TypeBeta,*TypeThresh;
  const DOUBLE *TypeORT;
  DOUBLE *Rest;
  SHORT *RestComp;
  INT flag;

  /* consistency check: diagonal blocks are supposed to be square matrices */
  for (type=0; type<NVECTYPES; type++)
    if (MD_ROWS_IN_RT_CT(M,type,type)>0)
    {
      nr = MD_ROWS_IN_RT_CT(M,type,type);

      ASSERT (nr*nr <= MAX_SINGLE_MAT_COMP);
      /* if too little: increase MAX_SINGLE_VEC_COMP and recompile	*/
                        #ifdef NDEBUG
      if (nr*nr > MAX_SINGLE_MAT_COMP)
        /* check also in case NDEBUG is defined (assert off)	*/
        REP_ERR_RETURN (__LINE__);
                        #endif
      if (nr != MD_COLS_IN_RT_CT(M,type,type))
        REP_ERR_RETURN (__LINE__);
    }
  /* check VD_rest iff */
  if (VD_rest!=NULL)
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
      {
        if (VD_NCMPS_IN_TYPE(VD_rest,type)==0)
          REP_ERR_RETURN (__LINE__);
        if (VD_NCMPS_IN_TYPE(VD_rest,type)!=MD_ROWS_IN_RT_CT(M,type,type))
          REP_ERR_RETURN (__LINE__);
      }

  /* consistency check:
     the transpose block-matrices (iff) must have the same format */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype+1; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        if (MD_ROWS_IN_RT_CT(M,rtype,rtype)!=MD_ROWS_IN_RT_CT(M,rtype,ctype))
          REP_ERR_RETURN (__LINE__);
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)!=MD_COLS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
        if (MD_COLS_IN_RT_CT(M,rtype,ctype)!=MD_ROWS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
      }

  /* calculate offsets for beta and threshold for different types */
  offset[0] = 0;
  for (type=1; type<NVECTYPES; type++)
    offset[type] = offset[type-1] + MD_ROWS_IN_RT_CT(M,type-1,type-1);

  /* flag vectors where connections will be introduced */
  if ((VD_rest!=NULL) && (oldrestthresh!=NULL))
  {
    L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
    {
      type     = VTYPE(vi);
      n        = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      TypeORT  = oldrestthresh+offset[type];
      Rest     = VVALUEPTR(vi,0);
      RestComp = VD_CMPPTR_OF_TYPE(VD_rest,type);

      flag = FALSE;
      for (i=0; i<n; i++)
        if (Rest[RestComp[i]]>TypeORT[i])
          flag = TRUE;

      SETVCUSED(vi,flag);
    }
  }
  else
    L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
    SETVCUSED(vi,0);

  /* clear VD_rest iff */
  if (VD_rest!=NULL)
    l_dset(g,VD_rest,EVERY_CLASS,0.0);

  /* decompose the matrix */
  if (MD_IS_SCALAR(M))
  {
    mc = MD_SCALCMP(M);
    mask = 0;
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
        mask |= 1<<type;

    /* loop over all lines */
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      /* check class and mask */
      if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;
      i = VINDEX(vi);

      /* now we are at line i */
      diag = MVALUE(VSTART(vi),mc);                                     /* diagonal element */
      if (fabs(diag)<SMALL_D) REP_ERR_RETURN(-i);                               /* decomposition failed */

      /* store inverse back to diag */
      if (StoreInverse)
        MVALUE(VSTART(vi),mc) = invdiag = 1.0/diag;
      else
        invdiag = 1.0/diag;

      /* eliminate all entries (j,i) with j>i */
      for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
      {
        /* check class and mask */
        vj = MDEST(Mij);
        if ( !((VDATATYPE(vj)&mask)&&(VCLASS(vj)>=ACTIVE_CLASS)&&(VINDEX(vj)>i)) ) continue;
        Mji = MADJ(Mij);
        AbsDjj = fabs(MVALUE(VSTART(vj),mc));
        pivot = MVALUE(Mji,mc)*invdiag;

        /* the pivot becomes the entry of the lower triangular part */
        MVALUE(Mji,mc) = pivot;

        if (pivot==0.0)
          continue;                                             /* nothing to eliminate */

        /* do gaussian elimination on the pattern (all ujk, k>i) */
        for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
        {
          vk = MDEST(Mik);
          if ( !((VDATATYPE(vk)&mask)&&(VCLASS(vk)>=ACTIVE_CLASS)&&(VINDEX(vk)>i)) ) continue;
          Mjk = GetMatrix(vj,vk);
          if (threshold!=NULL)                                  /* only if threshold is defined */
            if (Mjk==NULL)
              /* the threshold is taken relative to the diagonal of row j */
              if (fabs(pivot*MVALUE(Mik,mc))>threshold[0]*AbsDjj)
              {
                /* introduce new connection */
                Mjk = CreateExtraConnection(g,vj,vk);
                if (Mjk==NULL)
                  REP_ERR_RETURN (NUM_OUT_OF_MEM);
              }
          if (Mjk==NULL)
            if (VCUSED(vj))
            {
              /* introduce new connection */
              Mjk = CreateExtraConnection(g,vj,vk);
              if (Mjk==NULL)
                REP_ERR_RETURN (NUM_OUT_OF_MEM);
            }
          if (Mjk!=NULL)
            MVALUE(Mjk,mc) -= pivot*MVALUE(Mik,mc);                                                                       /* entry is in pattern */
          else
          {
            if (beta!=NULL)                                     /* only if beta is defined */
              MVALUE(VSTART(vj),mc) += beta[0]*fabs(pivot*MVALUE(Mik,mc));                                           /* entry not in pattern */
            if (VD_rest!=NULL)
              VVALUE(vj,VD_CMP_OF_TYPE(VD_rest,VTYPE(vj),0)) += fabs(pivot*MVALUE(Mik,mc));
          }
        }
      }
    }
    return (NUM_OK);
  }

  /* loop over all lines */
  L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
  {
    type     = VTYPE(vi);
    n        = MD_ROWS_IN_RT_CT(M,type,type);
    if (n == 0) continue;
    DiagComp = MD_MCMPPTR_OF_RT_CT(M,type,type);
    n2       = n*n;

    i = VINDEX(vi);

    Diag = MVALUEPTR(VSTART(vi),0);

    if (InvertSmallBlock(n,DiagComp,Diag,InvMat)!=0)
      REP_ERR_RETURN (-i);                              /* decompostion failed */

    /* write inverse back to diagonal block */
    if (StoreInverse)
      for (l=0; l<n2; l++)
        Diag[DiagComp[l]] = InvMat[l];

    /* eliminate all entries (j,i) with j>i */
    for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
    {
      rtype   = VTYPE(vj=MDEST(Mij));

      if (!((MD_ROWS_IN_RT_CT(M,rtype,type)>0)
            && (VCLASS(vj)>=ACTIVE_CLASS)
            && (i<VINDEX(vj))))
        continue;

      DjjComp = MD_MCMPPTR_OF_RT_CT(M,rtype,rtype);
      Djj             = MVALUEPTR(VSTART(vj),0);
      PivComp = MD_MCMPPTR_OF_RT_CT(M,rtype,type);
      nr              = MD_ROWS_IN_RT_CT(M,rtype,type);
      nnr             = nr*n;
      if (threshold!=NULL)                      /* only if threshold is defined */
        TypeThresh = threshold+offset[rtype];

      /* the normalization factors of the j diagonal */
      for (l=0; l<nr; l++)
        j_Normalization[l] = 1.0 / sqrt(fabs(Djj[DjjComp[l*nr+l]]));

      /* and we use a further vector to store the row sums of the normalized rest matrix */
      for (l=0; l<nr; l++) RowSum[l] = 0.0;

      Mji = MADJ(Mij);
      Piv = MVALUEPTR(Mji,0);

      /* matrix multiplication */
      PivIsZero = TRUE;
      for (i0=0; i0<nr; i0++)
        for (j0=0; j0<n; j0++)
        {
          sum = 0.0;
          for (k0=0; k0<n; k0++)
            sum += Piv[PivComp[i0*n+k0]] * InvMat[k0*n+j0];
          PivMat[i0*n+j0] = sum;
          if (sum!=0.0)
            PivIsZero = FALSE;
        }

      /* the pivot becomes the entry of the lower triangular part */
      for (l=0; l<nnr; l++)
        Piv[PivComp[l]] = PivMat[l];

      if (PivIsZero)
        continue;                                       /* nothing to eliminate */

      /* do gaussian elimination on the pattern (all Mjk, k>i) */
      for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
      {
        ctype   = VTYPE(vk=MDEST(Mik));

        if (!((MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
              && (VCLASS(vk)>=ACTIVE_CLASS)
              && (i<VINDEX(vk))))
          continue;

        DkkComp = MD_MCMPPTR_OF_RT_CT(M,ctype,ctype);
        Dkk             = MVALUEPTR(VSTART(vk),0);
        ElmComp = MD_MCMPPTR_OF_RT_CT(M,type,ctype);
        nc              = MD_COLS_IN_RT_CT(M,type,ctype);
        nrnc    = nr*nc;

        MatComp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);

        Elm = MVALUEPTR(Mik,0);

        /* matrix multiplication */
        CorIsZero = TRUE;
        for (i0=0; i0<nr; i0++)
          for (j0=0; j0<nc; j0++)
          {
            sum = 0.0;
            for (k0=0; k0<n; k0++)
              sum += PivMat[i0*n+k0] * Elm[ElmComp[k0*nc+j0]];
            CorMat[i0*nc+j0] = sum;
            if (sum!=0.0)
              CorIsZero = FALSE;
          }

        if (CorIsZero)
          continue;                                             /* nothing to correct */

        /* the normalization factors of the k diagonal */
        if (ctype==rtype)
          /* normalize with diag of vj */
          for (l=0; l<nr; l++)
            k_Normalization[l] = j_Normalization[l];                                            /* this choice led to good results -*/
        /* but why not else in all cases?	*/
        else
          /* normalize with diag of vk */
          for (l=0; l<nc; l++)
            k_Normalization[l] = 1.0 / sqrt(fabs(Dkk[DkkComp[l*nc+l]]));
        Mjk = GetMatrix(vj,vk);
        if (threshold!=NULL)                            /* only if threshold is defined */
          if (Mjk==NULL)
          {
            /* check threshold vs row sums of rest matrix */
            for (l=0; l<nr; l++)
            {
              sum = 0.0;
              for (m=0; m<nc; m++)
                sum += fabs(CorMat[l*nc+m]*j_Normalization[l]*k_Normalization[m]);
              if (sum>TypeThresh[l])
                if (CreateExtraConnection(g,vj,vk)==NULL)
                  REP_ERR_RETURN (NUM_OUT_OF_MEM)
                  else
                    break;
            }
          }
        if (Mjk==NULL)
          if (VCUSED(vj))
          {
            /* introduce new connection */
            Mjk = CreateExtraConnection(g,vj,vk);
            if (Mjk==NULL)
              REP_ERR_RETURN (NUM_OUT_OF_MEM);
          }
        if (Mjk!=NULL)
        {
          /* we are on the pattern: subtract the entry from Mjk */
          Mat = MVALUEPTR(Mjk,0);
          for (l=0; l<nrnc; l++)
            Mat[MatComp[l]] -= CorMat[l];
        }
        else
        {
          /* we are off the pattern: add rowsum */
          for (l=0; l<nr; l++)
            for (m=0; m<nc; m++)
              RowSum[l] += fabs(CorMat[l*nc+m]*j_Normalization[l]*k_Normalization[m]);
        }
      }
      if (beta==NULL) continue;
      /* only if beta is defined */

      /* add to VD_rest iff */
      if (VD_rest!=NULL)
        for (m=0; m<n; m++)
          VVALUE(vj,VD_CMP_OF_TYPE(VD_rest,type,m)) += RowSum[m];

      /* finally we modify the diagonal Djj */
      /* NB: the diag of row j is inverted later, so we don't have to care about StoreInverse */

      /* the diagonal damp matrix */
      TypeBeta= beta+offset[rtype];
      for (m=0; m<nr; m++)
        Damp[m] = 1.0 + TypeBeta[m]*RowSum[m];

      /* Djj = Djj*Damp */
      for (m=0; m<nr; m++)
        for (l=0; l<nr; l++)
          Djj[DjjComp[m*nr+l]] *= Damp[l];
    }
  }

  return (NUM_OK);
}

INT l_ilubdecomp_SB (BLOCKVECTOR *theBV, const MATDATA_DESC *M, const VEC_SCALAR beta)
{
  VECTOR *vi,*vj,*vk,*first_vec,*last_vec;
  MATRIX *Mij,*Mji,*Mjk,*Mik;
  DOUBLE InvMat[MAX_SINGLE_MAT_COMP],PivMat[MAX_SINGLE_MAT_COMP];
  DOUBLE CorMat[MAX_SINGLE_MAT_COMP];
  DOUBLE sum;
  INT offset[NVECTYPES+1],last_index;
  register DOUBLE *Diag,*Piv,*Elm,*Mat,*Djj,*Dkk;
  register SHORT *DiagComp,*PivComp,*ElmComp,*MatComp,*DjjComp,*DkkComp;
  register INT i0,j0,k0,l,m;
  INT type,ctype,rtype,PivIsZero,CorIsZero;
  INT n,n2,nr,nnr,nc,nrnc;
  INT i,mc,mask;
  DOUBLE RowSum[MAX_SINGLE_VEC_COMP],Damp[MAX_SINGLE_VEC_COMP];
  VEC_SCALAR j_Normalization,k_Normalization;
  DOUBLE diag,invdiag,pivot,AbsDjj;
  const DOUBLE *TypeBeta,*TypeThresh;
  const DOUBLE *TypeORT;
  DOUBLE *Rest;
  SHORT *RestComp;
  INT flag;

  /* consistency check: diagonal blocks are supposed to be square matrices */
  for (type=0; type<NVECTYPES; type++)
    if (MD_ROWS_IN_RT_CT(M,type,type)>0)
    {
      nr = MD_ROWS_IN_RT_CT(M,type,type);

      ASSERT (nr*nr <= MAX_SINGLE_MAT_COMP);
      /* if too little: increase MAX_SINGLE_VEC_COMP and recompile	*/
                        #ifdef NDEBUG
      if (nr*nr > MAX_SINGLE_MAT_COMP)
        /* check also in case NDEBUG is defined (assert off)	*/
        REP_ERR_RETURN (__LINE__);
                        #endif
      if (nr != MD_COLS_IN_RT_CT(M,type,type))
        REP_ERR_RETURN (__LINE__);
    }

  /* consistency check:
     the transpose block-matrices (iff) must have the same format */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype+1; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        if (MD_ROWS_IN_RT_CT(M,rtype,rtype)!=MD_ROWS_IN_RT_CT(M,rtype,ctype))
          REP_ERR_RETURN (__LINE__);
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)!=MD_COLS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
        if (MD_COLS_IN_RT_CT(M,rtype,ctype)!=MD_ROWS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
      }

  /* calculate offsets for beta and threshold for different types */
  offset[0] = 0;
  for (type=1; type<NVECTYPES; type++)
    offset[type] = offset[type-1] + MD_ROWS_IN_RT_CT(M,type-1,type-1);

  /* decompose the matrix */
  first_vec = BVFIRSTVECTOR(theBV);
  last_vec = BVLASTVECTOR(theBV);
  last_index = VINDEX(last_vec);

  if (MD_IS_SCALAR(M))
  {
    mc = MD_SCALCMP(M);
    mask = 0;
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
        mask |= 1<<type;

    /* loop over all lines */
    for (vi=first_vec; vi!=SUCCVC(last_vec); vi=SUCCVC(vi))
    {
      /* check class and mask */
      if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;
      i = VINDEX(vi);

      /* now we are at line i */
      diag = MVALUE(VSTART(vi),mc);                                                     /* diagonal element */
      if (fabs(diag)<SMALL_D) REP_ERR_RETURN(-i);                               /* decomposition failed */
      invdiag = 1.0/diag;

      /* eliminate all entries (j,i) with j>i */
      for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
      {
        /* check class and mask */
        vj = MDEST(Mij);
        if ( !((VDATATYPE(vj)&mask)&&(VCLASS(vj)>=ACTIVE_CLASS)&&(VINDEX(vj)>i)&&(VINDEX(vj)<=last_index)) ) continue;
        Mji = MADJ(Mij);
        AbsDjj = fabs(MVALUE(VSTART(vj),mc));
        pivot = MVALUE(Mji,mc)*invdiag;

        /* the pivot becomes the entry of the lower triangular part */
        MVALUE(Mji,mc) = pivot;

        if (pivot==0.0)
          continue;                                             /* nothing to eliminate */

        /* do gaussian elimination on the pattern (all ujk, k>i) */
        for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
        {
          vk = MDEST(Mik);
          if ( !((VDATATYPE(vk)&mask)&&(VCLASS(vk)>=ACTIVE_CLASS)&&(VINDEX(vk)>i)&&(VINDEX(vk)<=last_index)) ) continue;
          Mjk = GetMatrix(vj,vk);
          if (Mjk!=NULL)
            MVALUE(Mjk,mc) -= pivot*MVALUE(Mik,mc);                                                                       /* entry is in pattern */
          else
          {
            if (beta!=NULL)                                     /* only if beta is defined */
              MVALUE(VSTART(vj),mc) += beta[0]*fabs(pivot*MVALUE(Mik,mc));                                           /* entry not in pattern */
          }
        }
      }
    }
    return (NUM_OK);
  }

  REP_ERR_RETURN (1);
}

/****************************************************************************/
/*D
   l_icdecomp	- compute incomplete decomposition

   SYNOPSIS:
   INT l_icdecomp (GRID *g, const MATDATA_DESC *M);

   PARAMETERS:
   .  g - pointer to grid
   .  M - type matrix descriptor

   DESCRIPTION:
   This function computes an incomplete decomposition of order 0.
   The matrix M is overwritten!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    i<0 if decomposition failed
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

#define SMALL_ROOT 0.001

static INT SquareRootOfSmallBlock(SHORT n, const SHORT *mcomp,
                                  DOUBLE *mat, DOUBLE *SqRoot)
{
  DOUBLE mat1[MAX_SINGLE_MAT_COMP],sum;
  DOUBLE mat2[MAX_SINGLE_MAT_COMP];
  SHORT comp[MAX_SINGLE_MAT_COMP];
  INT i,j,k,l,iter;

  if (n==1)
  {
    SqRoot[0] = 1.0 / sqrt(mat[mcomp[0]]);
    return(0);
  }

  if (n > MAX_SINGLE_MAT_COMP)
    REP_ERR_RETURN(1);

  for (i=0; i<n*n; i++)
    comp[i] = i;

  /* Newton */

  if (InvertSmallBlock(n,mcomp,mat,SqRoot))
    REP_ERR_RETURN(1);

  for (i=0; i<n; i++)
    SqRoot[i*n+i] += 1.0;

  for (iter=0; iter<2*n*n; iter++)
  {
    for (i=0; i<n; i++)
    {
      for (j=0; j<i; j++)
      {
        sum = 0.0;
        for (k=0; k<n; k++)
          for (l=0; l<n; l++)
            sum += SqRoot[i*n+k] * mat[mcomp[k*n+l]] * SqRoot[l*n+j];
        mat2[i*n+j] = sum;
        mat2[j*n+i] = sum;
      }
      sum = 0.0;
      for (k=0; k<n; k++)
        for (l=0; l<n; l++)
          sum += SqRoot[i*n+k] * mat[mcomp[k*n+l]] * SqRoot[l*n+i];
      mat2[i*n+i] = sum;
    }

    sum = 0.0;
    for (i=0; i<n; i++)
    {
      for (j=0; j<i; j++)
        sum += ABS(mat2[i*n+j]);
      sum += ABS(mat2[i*n+i]-1.0);
    }

    if (sum < SMALL_ROOT)
      return(0);

    if (InvertSmallBlock(n,comp,mat2,mat1))
      REP_ERR_RETURN(1);

    for (i=0; i<n; i++)
      mat1[i*n+i] += 1.0;

    if (iter%2 == 0)
      for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
          sum = 0.0;
          for (k=0; k<n; k++)
            sum += SqRoot[i*n+k] * mat1[k*n+j];
          mat2[i*n+j] = sum;
        }
    else
      for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
          sum = 0.0;
          for (k=0; k<n; k++)
            sum += mat1[i*n+k] * SqRoot[k*n+j];
          mat2[i*n+j] = sum;
        }

    for (i=0; i<n; i++)
    {
      for (j=0; j<i; j++)
      {
        sum = 0.25 * (mat2[i*n+j] + mat2[j*n+i]);
        SqRoot[i*n+j] = sum;
        SqRoot[j*n+i] = sum;
      }
      SqRoot[i*n+i] = 0.5 * mat2[i*n+i];
    }
  }

  REP_ERR_RETURN(1);
}

INT l_icdecomp (GRID *g, const MATDATA_DESC *M)
{
  VECTOR *vi,*vj,*vk;
  MATRIX *Mij,*Mjk,*Mik;
  DOUBLE InvMat[MAX_SINGLE_MAT_COMP];
  DOUBLE CorMat[MAX_SINGLE_MAT_COMP];
  DOUBLE sum;
  register DOUBLE *Dii,*Dij,*Dik,*Djk,*Dji;
  register SHORT *DiiComp,*DijComp,*DikComp,*DjkComp,*DjiComp;
  register INT i0,j0,k0,l;
  INT type,ctype,rtype,CorIsZero;
  SHORT n,n2,nr,nnr,nc;
  INT i,mc,mask;
  DOUBLE matdiag,invdiag;

  /* consistency check: diagonal blocks are supposed to be square matrices */
  for (type=0; type<NVECTYPES; type++)
    if (MD_ROWS_IN_RT_CT(M,type,type)>0)
    {
      nr = MD_ROWS_IN_RT_CT(M,type,type);

      ASSERT (nr*nr <= MAX_SINGLE_MAT_COMP);
      /* if too little: increase MAX_SINGLE_VEC_COMP and recompile	*/
                        #ifdef NDEBUG
      if (nr*nr > MAX_SINGLE_MAT_COMP)
        /* check also in case NDEBUG is defined (assert off)	*/
        REP_ERR_RETURN (__LINE__);
                        #endif
      if (nr != MD_COLS_IN_RT_CT(M,type,type))
        REP_ERR_RETURN (__LINE__);
    }

  /* consistency check:
     the transpose block-matrices (iff) must have the same format */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype+1; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        if (MD_ROWS_IN_RT_CT(M,rtype,rtype)
            !=MD_ROWS_IN_RT_CT(M,rtype,ctype))
          REP_ERR_RETURN (__LINE__);
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)
            !=MD_COLS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
        if (MD_COLS_IN_RT_CT(M,rtype,ctype)
            !=MD_ROWS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
      }

  /* decompose the matrix */
  if (MD_IS_SCALAR(M))
  {
    mc = MD_SCALCMP(M);
    mask = 0;
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
        mask |= 1<<type;

    /* loop over all lines */
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      if (VECSKIP(vi))
        continue;

      /* check class and mask */
      if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) )
        continue;
      i = VINDEX(vi);

      /* now we are at line i */
      matdiag = MVALUE(VSTART(vi),mc);                                  /* diagonal element */

      /* k < i */
      for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
      {
        vk = MDEST(Mik);
        if (VECSKIP(vk))
          continue;
        if ( !((VDATATYPE(vk)&mask)
               &&(VCLASS(vk)>=ACTIVE_CLASS)
               &&(VINDEX(vk)<i)) )
          continue;
        sum = MVALUE(Mik,mc);
        matdiag -= sum * sum;
      }

      if (matdiag<SMALL_D) REP_ERR_RETURN(-i);                                  /* decomposition failed */

      /* store sqrt of inverse back to diag */
      MVALUE(VSTART(vi),mc) = invdiag = 1.0/sqrt(matdiag);

      /* j > i */
      for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
      {
        /* check class and mask */
        vj = MDEST(Mij);
        if (VECSKIP(vj))
          continue;
        if ( !((VDATATYPE(vj)&mask)
               &&(VCLASS(vj)>=ACTIVE_CLASS)
               &&(VINDEX(vj)>i)) )
          continue;
        sum = MVALUE(Mij,mc);

        /* k < i */
        for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
        {
          vk = MDEST(Mik);
          if (VECSKIP(vk))
            continue;
          if ( !((VDATATYPE(vk)&mask)
                 &&(VCLASS(vk)>=ACTIVE_CLASS)
                 &&(VINDEX(vk)<i)) )
            continue;
          Mjk = GetMatrix(vj,vk);
          if (Mjk==NULL)
            continue;
          sum -= MVALUE(Mik,mc) * MVALUE(Mjk,mc);
        }
        MVALUE(MADJ(Mij),mc) = sum * invdiag;
      }
    }

    return (NUM_OK);
  }

  /* loop over all lines */
  L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
  {
    if (VECSKIP(vi))
      continue;
    type     = VTYPE(vi);
    n        = MD_ROWS_IN_RT_CT(M,type,type);
    if (n == 0) continue;
    n2       = n*n;
    i = VINDEX(vi);
    Dii = MVALUEPTR(VSTART(vi),0);
    DiiComp = MD_MCMPPTR_OF_RT_CT(M,type,type);

    /* k < i */
    for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
    {
      vk = MDEST(Mik);
      if (VECSKIP(vk))
        continue;
      ctype   = VTYPE(vk);
      if (!((MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
            &&(VCLASS(vk)>=ACTIVE_CLASS)
            &&(VINDEX(vk)<i)) )
        continue;
      DikComp = MD_MCMPPTR_OF_RT_CT(M,type,ctype);
      nc              = MD_COLS_IN_RT_CT(M,type,ctype);
      Dik = MVALUEPTR(Mik,0);
      CorIsZero = TRUE;

      /* matrix multiplication */
      for (i0=0; i0<n; i0++)
        for (j0=0; j0<n; j0++)
        {
          sum = 0.0;
          for (k0=0; k0<nc; k0++)
            sum += Dik[DikComp[i0*nc+k0]] * Dik[DikComp[j0*nc+k0]];
          CorMat[i0*n+j0] = sum;
          if (sum!=0.0)
            CorIsZero = FALSE;
        }

      if (CorIsZero)
        continue;                               /* nothing to correct */

      /* subtract correction */
      for (l=0; l<n2; l++)
        Dii[DiiComp[l]] -= CorMat[l];
    }

    if (SquareRootOfSmallBlock(n,DiiComp,Dii,InvMat)!=0)
      REP_ERR_RETURN (-i);                              /* decompostion failed */

    /* write inverse back to diagonal block */
    for (l=0; l<n2; l++)
      Dii[DiiComp[l]] = InvMat[l];

    /* eliminate all entries (j,i) with j>i */
    for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
    {
      vj=MDEST(Mij);
      if (VECSKIP(vj))
        continue;
      rtype   = VTYPE(vj);

      if (!((MD_ROWS_IN_RT_CT(M,rtype,type)>0)
            && (VCLASS(vj)>=ACTIVE_CLASS)
            && (i<VINDEX(vj))))
        continue;

      Dij = MVALUEPTR(Mij,0);
      DijComp = MD_MCMPPTR_OF_RT_CT(M,type,rtype);
      nr              = MD_ROWS_IN_RT_CT(M,type,rtype);
      nnr             = nr*n;

      /* k < i */
      for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
      {
        vk = MDEST(Mik);
        if (VECSKIP(vk))
          continue;
        ctype   = VTYPE(vk);
        if (!((MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
              &&(VCLASS(vk)>=ACTIVE_CLASS)
              &&(VINDEX(vk)<i)) )
          continue;
        Mjk = GetMatrix(vj,vk);
        if (Mjk==NULL)
          continue;
        DikComp = MD_MCMPPTR_OF_RT_CT(M,type,ctype);
        nc = MD_COLS_IN_RT_CT(M,type,ctype);
        Dik = MVALUEPTR(Mik,0);
        DjkComp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
        Djk = MVALUEPTR(Mjk,0);
        CorIsZero = TRUE;

        /* matrix multiplication */
        for (i0=0; i0<n; i0++)
          for (j0=0; j0<nr; j0++)
          {
            sum = 0.0;
            for (k0=0; k0<nc; k0++)
              sum += Dik[DikComp[i0*nc+k0]] * Djk[DjkComp[j0*nc+k0]];
            CorMat[i0*nr+j0] = sum;
            if (sum!=0.0)
              CorIsZero = FALSE;
          }

        if (CorIsZero)
          continue;                                     /* nothing to correct */

        /* subtract correction */
        for (l=0; l<nnr; l++)
          Dij[DijComp[l]] -= CorMat[l];
      }

      Dji = MVALUEPTR(MADJ(Mij),0);
      DjiComp = MD_MCMPPTR_OF_RT_CT(M,rtype,type);

      for (i0=0; i0<n; i0++)
        for (j0=0; j0<nr; j0++)
        {
          sum = 0.0;
          for (k0=0; k0<n; k0++)
            sum += InvMat[i0*n+k0] * Dij[DijComp[k0*nr+j0]];
          Dji[DjiComp[j0*n+i0]] = sum;
        }
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
    l_iluspdecomp - compute spectrally shifted incomplete decomposition of order 0

        SYNOPSIS:
        INT l_iluspdecomp (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta,
        const VECDATA_DESC *t, INT mode);

        PARAMETERS:
   .   g - grid level
   .   M - matrix to decompose
   .   t - temp vector to store shifting parameters
   .   mode - shift 'SP_GLOBAL' or 'SP_LOCAL'
   .   oldrestthresh - if !=NULL use t field of last step (it contains the normalized row sums of the
                                                rest natrix of the decomposition) to flag vectors where
                                                connections will be introduced (CAUTION: do not overwrite t before!)

        DESCRIPTION:
        The row sums of the rest matrix are used to modify the diagonal of the
        decomposition. This is done either globally or locally.

        RETURN VALUE:
        INT
   .n   NUM_OK when o.k.

   D*/
/****************************************************************************/

INT l_iluspdecomp (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VECDATA_DESC *t, INT mode, const VEC_SCALAR oldrestthresh)
{
  VECTOR *vi,*vj,*vk;
  MATRIX *Mij,*Mji,*Mjk,*Mik;
  DOUBLE InvMat[MAX_SINGLE_MAT_COMP],PivMat[MAX_SINGLE_MAT_COMP],CorMat[MAX_SINGLE_MAT_COMP];
  DOUBLE sum;
  INT offset[NVECTYPES+1];
  register DOUBLE *Diag,*Piv,*Elm,*Mat,*Tmp,*Djj,*Dkk;
  register SHORT *DiagComp,*PivComp,*ElmComp,*MatComp,*TmpComp,*DjjComp,*DkkComp;
  register INT i0,j0,k0;
  INT type,ctype,rtype,PivIsZero,CorIsZero;
  INT n,n2,nr,nnr,nc,nrnc;
  INT i,j,l,m,mc,tc,mask,found;
  DOUBLE diag,invdiag,pivot,min,max,Min[MAX_SINGLE_VEC_COMP],Max[MAX_SINGLE_VEC_COMP],*TypeMin,*TypeMax;
  DOUBLE RowSum[MAX_SINGLE_VEC_COMP],Damp[MAX_SINGLE_VEC_COMP];
  VEC_SCALAR j_Normalization,k_Normalization;
  const DOUBLE *TypeBeta;
  const DOUBLE *TypeORT;
  DOUBLE *Rest;
  SHORT *RestComp;
  INT flag;

  if (t==NULL)
    REP_ERR_RETURN (1);

  /* consistency check: diagonal blocks are supposed to be square matrices,
     t should be at least of same size */
  found = 0;
  for (type=0; type<NVECTYPES; type++)
    if (MD_ROWS_IN_RT_CT(M,type,type)>0)
    {
      found ++;
      nr = MD_ROWS_IN_RT_CT(M,type,type);

      if (VD_NCMPS_IN_TYPE(t,type)==0)
        REP_ERR_RETURN (__LINE__);
      if (VD_NCMPS_IN_TYPE(t,type)<nr)
        REP_ERR_RETURN (__LINE__);

      ASSERT (nr*nr <= MAX_SINGLE_MAT_COMP);
      /* if too little: increase MAX_SINGLE_VEC_COMP and recompile	*/
                        #ifdef NDEBUG
      if (nr*nr > MAX_SINGLE_MAT_COMP)
        /* check also in case NDEBUG is defined (assert off)	*/
        REP_ERR_RETURN (__LINE__);
                        #endif
      if (nr != MD_COLS_IN_RT_CT(M,type,type))
        REP_ERR_RETURN (__LINE__);
    }

  ASSERT(found==1);                     /* see the CAUTION above */

  /* consistency check: the transpose block-matrices
     (iff) must have the same format */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype+1; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        if (MD_ROWS_IN_RT_CT(M,rtype,rtype)!=MD_ROWS_IN_RT_CT(M,rtype,ctype))
          REP_ERR_RETURN (__LINE__);
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)!=MD_COLS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
        if (MD_COLS_IN_RT_CT(M,rtype,ctype)!=MD_ROWS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
      }

  /* calculate offsets for beta for different types */
  offset[0] = 0;
  for (type=1; type<NVECTYPES; type++)
    offset[type] = offset[type-1] + MD_ROWS_IN_RT_CT(M,type-1,type-1);

  /* flag vectors where connections will be introduced */
  if (oldrestthresh!=NULL)
  {
    /* t is supposed to contain the normalized rest matrix row sums from the last decomposition */
    L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
    {
      type     = VTYPE(vi);
      n        = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      TypeORT  = oldrestthresh+offset[type];
      Rest     = VVALUEPTR(vi,0);
      RestComp = VD_CMPPTR_OF_TYPE(t,type);

      flag = FALSE;
      for (i=0; i<n; i++)
        if (Rest[RestComp[i]]>TypeORT[i])
          flag = TRUE;

      SETVCUSED(vi,flag);
    }
  }
  else
    L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
    SETVCUSED(vi,0);

  /* clear temp (here we acumulate the rowsums of the rest matrix) */
  if (l_dset(g,t,EVERY_CLASS,0.0)) REP_ERR_RETURN (__LINE__);


  /* decompose the matrix */
  if (MD_IS_SCALAR(M))
  {
    mc = MD_SCALCMP(M);
    tc = VD_SCALCMP(t);
    mask = 0;
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
        mask |= 1<<type;

    /* loop over all lines */
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      /* check class and mask */
      if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;
      i = VINDEX(vi);

      /* now we are at line i */
      diag = MVALUE(VSTART(vi),mc);                                     /* diagonal element */
      if (fabs(diag)<SMALL_D) REP_ERR_RETURN(-i);                               /* decomposition failed */

      invdiag = 1.0/diag;

      /* eliminate all entries (j,i) with j>i */
      for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
      {
        /* check class and mask */
        vj = MDEST(Mij);
        if ( !((VDATATYPE(vj)&mask)&&(VCLASS(vj)>=ACTIVE_CLASS)&&(VINDEX(vj)>i)) ) continue;
        Mji = MADJ(Mij);
        pivot = MVALUE(Mji,mc)*invdiag;

        /* the pivot becomes the entry of the lower triangular part */
        MVALUE(Mji,mc) = pivot;

        if (pivot==0.0)
          continue;                                             /* nothing to eliminate */

        /* do gaussian elimination on the pattern (all ujk, k>i) */
        for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
        {
          vk = MDEST(Mik);
          if ( !((VDATATYPE(vk)&mask)&&(VCLASS(vk)>=ACTIVE_CLASS)&&(VINDEX(vk)>i)) ) continue;
          Mjk = GetMatrix(vj,vk);
          if (Mjk==NULL)
            if (VCUSED(vj))
            {
              /* introduce new connection */
              Mjk = CreateExtraConnection(g,vj,vk);
              if (Mjk==NULL)
                REP_ERR_RETURN (NUM_OUT_OF_MEM);
            }
          if (Mjk!=NULL)
            MVALUE(Mjk,mc) -= pivot*MVALUE(Mik,mc);                                                                     /* entry is in pattern */
          else
            VVALUE(vj,tc) += fabs(pivot*MVALUE(Mik,mc));                                                        /* entry not in pattern*/
        }
      }
    }

    /* after we are done with the decomposition we modify the diagonal */
    /* temp contains the row sums of the rest matrix */
    if (mode==SP_GLOBAL)
    {
      /* find the maximum of the rowsums */
      max = 0.0;
      for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
      {
        /* we don't have to check class and mask since temp was cleared
           if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue; */

        max = MAX(max,VVALUE(vi,tc));
      }
      /* NB: max is the row sum norm of the rest matrix of the decomposition */

      /* globally modify the diagonal */
      if (StoreInverse)
        for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
        {
          /* check class and mask */
          if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;

          MVALUE(VSTART(vi),mc) = 1.0/(MVALUE(VSTART(vi),mc) + beta[0]*max);
        }
      else
        for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
        {
          /* check class and mask */
          if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;

          MVALUE(VSTART(vi),mc) += beta[0]*max;
        }
    }
    else
    {
      /* find the minimum>0 of the rowsums */
      min = MAX_D;
      for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
      {
        /* we don't have to check class and mask since temp was cleared
           if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue; */

        if (VVALUE(vi,tc)>SMALL_D)
          min = MIN(min,VVALUE(vi,tc));
      }
      /* locally modify the diagonal */
      if (StoreInverse)
        for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
        {
          /* check class and mask */
          if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;

          if (VVALUE(vi,tc)>min)
            MVALUE(VSTART(vi),mc) = 1.0/(MVALUE(VSTART(vi),mc) + beta[0]*VVALUE(vi,tc));
          else
            MVALUE(VSTART(vi),mc) = 1.0/(MVALUE(VSTART(vi),mc) + beta[0]*min);
        }
      else
        for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
        {
          /* check class and mask */
          if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;

          if (VVALUE(vi,tc)>min)
            MVALUE(VSTART(vi),mc) += beta[0]*VVALUE(vi,tc);
          else
            MVALUE(VSTART(vi),mc) += beta[0]*min;
        }
    }

    return (NUM_OK);
  }

  /* loop over all lines */
  L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
  {
    type     = VTYPE(vi);
    n        = MD_ROWS_IN_RT_CT(M,type,type);
    if (n == 0) continue;
    DiagComp = MD_MCMPPTR_OF_RT_CT(M,type,type);
    n2       = n*n;

    i = VINDEX(vi);

    Diag = MVALUEPTR(VSTART(vi),0);

    if (InvertSmallBlock(n,DiagComp,Diag,InvMat)!=0)
      REP_ERR_RETURN (-i);                              /* decompostion failed */

    /* write inverse back to diagonal block */
    if (StoreInverse)
      for (l=0; l<n2; l++)
        Diag[DiagComp[l]] = InvMat[l];

    /* eliminate all entries (j,i) with j>i */
    for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
    {
      rtype   = VTYPE(vj=MDEST(Mij));

      if (!((MD_ROWS_IN_RT_CT(M,rtype,type)>0)
            && (VCLASS(vj)>=ACTIVE_CLASS)
            && (i<(j=VINDEX(vj)))))
        continue;

      DjjComp = MD_MCMPPTR_OF_RT_CT(M,rtype,rtype);
      Djj             = MVALUEPTR(VSTART(vj),0);
      TmpComp = VD_CMPPTR_OF_TYPE(t,rtype);
      Tmp             = VVALUEPTR(vj,0);
      PivComp = MD_MCMPPTR_OF_RT_CT(M,rtype,type);
      nr              = MD_ROWS_IN_RT_CT(M,rtype,type);
      nnr             = nr*n;

      /* the normalization factors of the j diagonal */
      for (l=0; l<nr; l++)
        j_Normalization[l] = 1.0 / sqrt(fabs(Djj[DjjComp[l*nr+l]]));

      /* and we use a further vector to store the row sums of the normalized rest matrix */
      for (l=0; l<nr; l++) RowSum[l] = 0.0;

      Mji = MADJ(Mij);
      Piv = MVALUEPTR(Mji,0);

      /* matrix multiplication */
      PivIsZero = TRUE;
      for (i0=0; i0<nr; i0++)
        for (j0=0; j0<n; j0++)
        {
          sum = 0.0;
          for (k0=0; k0<n; k0++)
            sum += Piv[PivComp[i0*n+k0]] * InvMat[k0*n+j0];
          PivMat[i0*n+j0] = sum;
          if (sum!=0.0)
            PivIsZero = FALSE;
        }

      /* the pivot becomes the entry of the lower triangular part */
      for (l=0; l<nnr; l++)
        Piv[PivComp[l]] = PivMat[l];

      if (PivIsZero)
        continue;                                       /* nothing to eliminate */

      /* do gaussian elimination on the pattern (all Mjk, k>i) */
      for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
      {
        ctype   = VTYPE(vk=MDEST(Mik));

        if (!((MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
              && (VCLASS(vk)>=ACTIVE_CLASS)
              && (i<VINDEX(vk))))
          continue;

        DkkComp = MD_MCMPPTR_OF_RT_CT(M,ctype,ctype);
        Dkk             = MVALUEPTR(VSTART(vk),0);
        ElmComp = MD_MCMPPTR_OF_RT_CT(M,type,ctype);
        nc              = MD_ROWS_IN_RT_CT(M,type,ctype);
        nrnc    = nr*nc;

        MatComp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);

        Elm = MVALUEPTR(Mik,0);

        /* matrix multiplication */
        CorIsZero = TRUE;
        for (i0=0; i0<nr; i0++)
          for (j0=0; j0<nc; j0++)
          {
            sum = 0.0;
            for (k0=0; k0<n; k0++)
              sum += PivMat[i0*n+k0] * Elm[ElmComp[k0*nc+j0]];
            CorMat[i0*nc+j0] = sum;
            if (sum!=0.0)
              CorIsZero = FALSE;
          }

        if (CorIsZero)
          continue;                                             /* nothing to correct */

        /* the normalization factors of the k diagonal */
        if (ctype==rtype)
          /* normalize with diag of vj */
          for (l=0; l<nr; l++)
            k_Normalization[l] = j_Normalization[l];                                            /* this choice led to good results -*/
        /* but why not else in all cases?	*/
        else
          /* normalize with diag of vk */
          for (l=0; l<nc; l++)
            k_Normalization[l] = 1.0 / sqrt(fabs(Dkk[DkkComp[l*nc+l]]));

        Mjk = GetMatrix(vj,vk);
        if (Mjk==NULL)
          if (VCUSED(vj))
          {
            /* introduce a new connection */
            Mjk = CreateExtraConnection(g,vj,vk);
            if (Mjk==NULL)
              REP_ERR_RETURN (NUM_OUT_OF_MEM);
          }
        if (Mjk!=NULL)
        {
          /* we are on the pattern: subtract the entry from Mjk */
          Mat = MVALUEPTR(Mjk,0);
          for (l=0; l<nrnc; l++)
            Mat[MatComp[l]] -= CorMat[l];
        }
        else
        {
          /* we are off the pattern: add rowsum */
          for (l=0; l<nr; l++)
            for (m=0; m<nc; m++)
              RowSum[l] += fabs(CorMat[l*nc+m]*j_Normalization[l]*k_Normalization[m]);
        }
      }

      /* finally we add the row sum of the normalized rest matrix to temp */
      for (m=0; m<nr; m++)
        Tmp[TmpComp[m]] += RowSum[m];
    }
  }

  /* after we are done with the decomposition we modify the diagonal */
  /* temp contains the row sums of the normalized rest matrix */
  if (mode==SP_GLOBAL)
  {
    /* find the maximum of the rowsums */
    for (i=0; i<MAX_SINGLE_VEC_COMP; i++) Max[i] = 0.0;
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      type     = VTYPE(vi);
      n       = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      TmpComp = VD_CMPPTR_OF_TYPE(t,type);
      Tmp             = VVALUEPTR(vi,0);

      TypeMax = Max+offset[type];
      for (i=0; i<n; i++)
        TypeMax[i] = MAX(TypeMax[i],Tmp[TmpComp[i]]);
    }

    /* globally modify the diagonal */
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      type    = VTYPE(vi);
      n       = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      DiagComp= MD_MCMPPTR_OF_RT_CT(M,type,type);
      Diag    = MVALUEPTR(VSTART(vi),0);
      TmpComp = VD_CMPPTR_OF_TYPE(t,type);
      Tmp             = VVALUEPTR(vi,0);

      TypeMax = Max+offset[type];
      TypeBeta= beta+offset[type];

      if (StoreInverse)
      {
        /* the inverse diagonal damp matrix */
        for (i=0; i<n; i++)
          Damp[i] = 1.0 / (1.0 + TypeBeta[i]*TypeMax[i]);

        /* InvDiag = InvDamp*InvDiag */
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
            Diag[DiagComp[i*n+j]] *= Damp[i];
      }
      else
      {
        /* the diagonal damp matrix */
        for (i=0; i<n; i++)
          Damp[i] = 1.0 + TypeBeta[i]*TypeMax[i];

        /* Diag = Diag*Damp */
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
            Diag[DiagComp[i*n+j]] *= Damp[j];
      }
    }
  }
  else
  {
    /* find the minimum>0 of the rowsums */
    for (i=0; i<MAX_SINGLE_VEC_COMP; i++) Min[i] = MAX_D;
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      type    = VTYPE(vi);
      n       = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      TmpComp = VD_CMPPTR_OF_TYPE(t,type);
      Tmp             = VVALUEPTR(vi,0);

      TypeMin = Min+offset[type];
      for (i=0; i<n; i++)
        if (Tmp[TmpComp[i]]>SMALL_D)
          TypeMin[i] = MIN(TypeMin[i],Tmp[TmpComp[i]]);
    }
    /* locally  modify the diagonal */
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      type     = VTYPE(vi);
      n       = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      DiagComp= MD_MCMPPTR_OF_RT_CT(M,type,type);
      Diag    = MVALUEPTR(VSTART(vi),0);
      TmpComp = VD_CMPPTR_OF_TYPE(t,type);
      Tmp             = VVALUEPTR(vi,0);

      TypeMin = Min+offset[type];
      TypeBeta= beta+offset[type];

      if (StoreInverse)
      {
        /* the inverse diagonal damp matrix */
        for (i=0; i<n; i++)
          if (Tmp[TmpComp[i]]>TypeMin[i])
            Damp[i] = 1.0 / (1.0 + TypeBeta[i]*Tmp[TmpComp[i]]);
          else
            Damp[i] = 1.0 / (1.0 + TypeBeta[i]*TypeMin[i]);

        /* InvDiag = InvDamp*InvDiag */
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
            Diag[DiagComp[i*n+j]] *= Damp[i];
      }
      else
      {
        /* the diagonal damp matrix */
        for (i=0; i<n; i++)
          if (Tmp[TmpComp[i]]>TypeMin[i])
            Damp[i] = 1.0 + TypeBeta[i]*Tmp[TmpComp[i]];
          else
            Damp[i] = 1.0 + TypeBeta[i]*TypeMin[i];

        /* Diag = Diag*Damp */
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
            Diag[DiagComp[i*n+j]] *= Damp[j];
      }
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_lrdecomp - compute left-right decomposition of M

   SYNOPSIS:
   INT l_lrdecomp (GRID *g, const MATDATA_DESC *M, const DOUBLE *threshold);

   PARAMETERS:
   .  g - pointer to grid
   .  M - type matrix descriptor for precondition

   DESCRIPTION:
   This function computes left-right decomposition of M.
   The matrix M is overwritten!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_OUT_OF_MEM if there is not enough memory to store the decomposition
   .n    i<0 if the decomposition fails at the 'VECTOR' with 'INDEX' '-i'
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_lrdecomp (GRID *g, const MATDATA_DESC *M)
{
  VECTOR *vi,*vj,*vk;
  MATRIX *Mij,*Mji,*Mjk,*Mik;
  DOUBLE InvMat[MAX_SINGLE_MAT_COMP],PivMat[MAX_SINGLE_MAT_COMP];
  DOUBLE CorMat[MAX_SINGLE_MAT_COMP];
  DOUBLE sum;
  register DOUBLE *Diag,*Piv,*Elm,*Mat;
  register SHORT *DiagComp,*PivComp,*ElmComp,*MatComp;
  register INT i0,j0,k0;
  INT type,ctype,rtype,PivIsZero,CorIsZero;
  INT n,n2,nr,nnr,nc,nrnc;
  INT i,l,mc,mask;
  DOUBLE diag,invdiag,pivot;

  /* consistency check: diagonal blocks are supposed to be square matrices */
  for (type=0; type<NVECTYPES; type++)
    if (MD_ROWS_IN_RT_CT(M,type,type)>0)
    {
      nr = MD_ROWS_IN_RT_CT(M,type,type);
      ASSERT (nr*nr <= MAX_SINGLE_MAT_COMP);
      /* if too little: increase MAX_SINGLE_VEC_COMP and recompile */
                        #ifdef NDEBUG
      if (nr*nr > MAX_SINGLE_MAT_COMP)
        /* check also in case NDEBUG is defined (assert off)	*/
        REP_ERR_RETURN (__LINE__);
                        #endif
      if (nr != MD_COLS_IN_RT_CT(M,type,type))
        REP_ERR_RETURN (__LINE__);
    }

  /* consistency check: the transpose block-matrices (iff) must have
     the same format */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype+1; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        if (MD_ROWS_IN_RT_CT(M,rtype,rtype)!=MD_ROWS_IN_RT_CT(M,rtype,ctype))
          REP_ERR_RETURN (__LINE__);
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)!=MD_COLS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
        if (MD_COLS_IN_RT_CT(M,rtype,ctype)!=MD_ROWS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
      }

  /* consistency check: the block pattern should allow for elimination */
  for (type=0; type<NVECTYPES; type++)
    for (rtype=type+1; rtype<NVECTYPES; rtype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,type)>0)
        for (ctype=type+1; ctype<NVECTYPES; ctype++)
          if ((MD_ROWS_IN_RT_CT(M,type,ctype)>0) &&
              (MD_ROWS_IN_RT_CT(M,rtype,ctype)==0))
            REP_ERR_RETURN (__LINE__);

  /* decompose the matrix */
  if (MD_IS_SCALAR(M))
  {
    mc    = MD_SCALCMP(M);
    mask = 0;
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
        mask |= 1<<type;

    /* loop over all lines */
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      /* check class and mask */
      if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;
      i = VINDEX(vi);

      /* now we are at line i */
      diag = MVALUE(VSTART(vi),mc);                                     /* diagonal element */
      if (fabs(diag)<SMALL_D) REP_ERR_RETURN(-i);                               /* decomposition failed */

      /* store inverse back to diag */
      if (StoreInverse)
        MVALUE(VSTART(vi),mc) = invdiag = 1.0/diag;
      else
        invdiag = 1.0/diag;

      /* eliminate all entries (j,i) with j>i */
      for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
      {
        /* check class and mask */
        vj = MDEST(Mij);
        if ( !((VDATATYPE(vj)&mask)&&(VCLASS(vj)>=ACTIVE_CLASS)&&(VINDEX(vj)>i)) ) continue;
        Mji = MADJ(Mij);
        pivot = MVALUE(Mji,mc)*invdiag;

        /* the pivot becomes the entry of the lower triangular part */
        MVALUE(Mji,mc) = pivot;

        if (pivot==0.0)
          continue;                                             /* nothing to eliminate */

        /* do gaussian elimination on the pattern (all ujk, k>i) */
        for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
        {
          vk = MDEST(Mik);
          if ( !((VDATATYPE(vk)&mask)&&(VCLASS(vk)>=ACTIVE_CLASS)&&(VINDEX(vk)>i)) ) continue;
          Mjk = GetMatrix(vj,vk);
          if (Mjk==NULL)
          {
            Mjk = CreateExtraConnection(g,vj,vk);
            if (Mjk==NULL)
              REP_ERR_RETURN (NUM_OUT_OF_MEM);
          }
          MVALUE(Mjk,mc) -= pivot*MVALUE(Mik,mc);                                                                 /* entry is in pattern */
        }
      }
    }
    return (NUM_OK);
  }

  /* loop over all lines */
  L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
  {
    type = VTYPE(vi);
    n        = MD_ROWS_IN_RT_CT(M,type,type);
    if (n == 0) continue;
    DiagComp = MD_MCMPPTR_OF_RT_CT(M,type,type);
    n2       = n*n;

    i = VINDEX(vi);

    Diag = MVALUEPTR(VSTART(vi),0);

    if (InvertSmallBlock(n,DiagComp,Diag,InvMat)!=0)
      REP_ERR_RETURN (-i);                              /* decompostion failed */

    /* write inverse back to diagonal block */
    if (StoreInverse)
      for (l=0; l<n2; l++)
        Diag[DiagComp[l]] = InvMat[l];

    /* eliminate all entries (j,i) with j>i */
    for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
    {
      rtype = VTYPE(vj=MDEST(Mij));

      if (!((MD_ROWS_IN_RT_CT(M,rtype,type)>0)
            && (VCLASS(vj)>=ACTIVE_CLASS)
            && (i<VINDEX(vj))))
        continue;

      PivComp = MD_MCMPPTR_OF_RT_CT(M,rtype,type);
      nr              = MD_ROWS_IN_RT_CT(M,rtype,type);
      nnr             = nr*n;

      Mji = MADJ(Mij);
      Piv = MVALUEPTR(Mji,0);

      /* matrix multiplication */
      PivIsZero = TRUE;
      for (i0=0; i0<nr; i0++)
        for (j0=0; j0<n; j0++)
        {
          sum = 0.0;
          for (k0=0; k0<n; k0++)
            sum += Piv[PivComp[i0*n+k0]] * InvMat[k0*n+j0];
          PivMat[i0*n+j0] = sum;
          if (sum!=0.0)
            PivIsZero = FALSE;
        }

      /* the pivot becomes the entry of the lower triangular part */
      for (l=0; l<nnr; l++)
        Piv[PivComp[l]] = PivMat[l];

      if (PivIsZero)
        continue;                                       /* nothing to eliminate */

      /* do gaussian elimination on the pattern (all Mjk, k>i) */
      for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
      {
        ctype   = VTYPE(vk=MDEST(Mik));

        if (!((MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
              && (VCLASS(vk)>=ACTIVE_CLASS)
              && (i<VINDEX(vk))))
          continue;

        ElmComp = MD_MCMPPTR_OF_RT_CT(M,type,ctype);
        nc              = MD_COLS_IN_RT_CT(M,type,ctype);
        nrnc    = nr*nc;

        MatComp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);

        Elm = MVALUEPTR(Mik,0);

        /* matrix multiplication */
        CorIsZero = TRUE;
        for (i0=0; i0<nr; i0++)
          for (j0=0; j0<nc; j0++)
          {
            sum = 0.0;
            for (k0=0; k0<n; k0++)
              sum += PivMat[i0*n+k0] * Elm[ElmComp[k0*nc+j0]];
            CorMat[i0*nc+j0] = sum;
            if (sum!=0.0)
              CorIsZero = FALSE;
          }

        if (CorIsZero)
          continue;                                             /* nothing to correct */

        Mjk = GetMatrix(vj,vk);
        if (Mjk==NULL)
        {
          Mjk = CreateExtraConnection(g,vj,vk);
          if (Mjk==NULL)
            REP_ERR_RETURN (NUM_OUT_OF_MEM);
        }
        /* subtract the entry from Mjk */
        Mat = MVALUEPTR(Mjk,0);
        for (l=0; l<nrnc; l++)
          Mat[MatComp[l]] -= CorMat[l];
      }
    }
  }

  return (NUM_OK);
}

#ifdef __BLOCK_VECTOR_DESC__

/****************************************************************************/
/*D
   LUDecomposeDiagBS - calculate the LU decomposition of a diagonal block of a matrix

   SYNOPSIS:
   INT LUDecomposeDiagBS( const BLOCKVECTOR *bv, const BV_DESC *bvd,
                          const BV_DESC_FORMAT *bvdf, INT A_comp, GRID *grid );

   PARAMETERS:
   .  bv - blockvector of the blockmatrix
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  A_comp - position of the matrix in the MATRIX-data of the blockmatrix

   DESCRIPTION:
   Calculate the LU decomposition of the given block of matrix 'A'. The result
   is stored again in this block. The matrix must be a digonalblock of the
   global matrix 'A'. For the fill-in extra connections are allocated.
   For 0-entries no fill-in is allocated.

   The VINDEX-field must be ascendend in the VECTOR-chain.

   Diagonalelements < SMALL_D are considered as 0 and the corresponding
   division is refused with an error message.

   To solve for such a decomposition use 'solveLUMatBS'.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_SMALL_DIAG if a diagonal element was < SMALL_D (i.e. division by nearly 0)
   .n    NUM_OUT_OF_MEM not enough memory to allocate the fill-in entries

   SEE ALSO:
   BLOCKVECTOR, blas_routines, solveLUMatBS
   D*/
/****************************************************************************/

INT LUDecomposeDiagBS( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT A_comp, GRID *grid )
{
  register VECTOR *vi, *end_vi, *vj, *vk;
  register DOUBLE aii, aji, ajk_corr;
  register MATRIX *mij, *mik, *mjk;
  register CONNECTION *con;
  register INT extra_con = 0;

  ASSERT( bv != NULL );
  ASSERT( A_comp != -1 );

  /* loop over all lines */
  end_vi = BVENDVECTOR( bv );
  for ( vi = BVFIRSTVECTOR( bv ); vi != end_vi; vi = SUCCVC( vi ) )
  {
    /* search aii */
    /*for ( mij = VSTART( vi ); mij != NULL; mij = MNEXT( mij ) )
            if ( MDEST( mij ) == vi )
            {
                    aii = MVALUE( mij, A_comp );
                    break;
            }*/
    /* there must be a diagonal element; assert it */
    /*ASSERT( mij != NULL );*/
    aii = MVALUE(VSTART(vi),A_comp);

    if ( fabs(aii) < SMALL_D )
    {
      PrintErrorMessage( 'E', "LUDecomposeDiagBS", "Diagonal element too small in LUDecompDiagBS!\n" );
      REP_ERR_RETURN (NUM_SMALL_DIAG);
    }

    /* eliminate all entries aji in column i; to do this, walk through the
       matrixlist for the upper diagonal (corresponding to row i)
       and go from there to the adjoint */
    for ( mij = VSTART( vi ); mij != NULL; mij = MNEXT( mij ) )
    {
      vj = MDEST( mij );
      if ( (VINDEX(vj) > VINDEX(vi)) && VMATCH( vj, bvd, bvdf ) )
      {
        /* store the pivot element in the lower triangular part */
        aji = MVALUE( MADJ(mij), A_comp ) /= aii ;

        if ( aji == 0.0 )
          continue;                                             /* nothing to eliminate */

        /* do gaussian elimination for row j on the pattern */
        for ( mik = VSTART( vi ); mik != NULL; mik = MNEXT( mik ) )
        {
          vk = MDEST( mik );
          if ( (VINDEX(vk) > VINDEX(vi)) && VMATCH( vk, bvd, bvdf ) )
          {
            ajk_corr = aji * MVALUE( mik, A_comp );
            if ( fabs( ajk_corr ) >= SMALL_D )
            {
              if ( (mjk = GetMatrix( vj, vk )) == NULL )
              {
                if ( (con = CreateExtraConnection( grid, vj, vk )) == NULL )
                {
                  PrintErrorMessage( 'E', "LUDecomposeDiagBS", "Not enough memory" );
                  REP_ERR_RETURN (NUM_OUT_OF_MEM);
                }
                mjk = CMATRIX0( con );
                extra_con++;
              }
              MVALUE( mjk, A_comp) -= ajk_corr;
            }
            /* else neglect the entry */
          }
        }
      }
    }
  }

  if (  (extra_con > 0) && (GetMuteLevel() >= 100) )
    UserWriteF( "%d extra connection allocated in LUDecompDiagBS.\n", extra_con );

  return NUM_OK;
}

#endif

/****************************************************************************/
/*D
   l_lrregularize - regularize last diagonal block

   SYNOPSIS:
   INT l_lrregularize (GRID *theGrid, const MATDATA_DESC *M)

   PARAMETERS:
   .  theGrid - pointer to grid
   .  M - type matrix descriptor for decomposition

   DESCRIPTION:
   This function regularizes the last diagonal block if it was found to be singular.
   It checks wether it is singular exactly for one component and replaces the zero by 1.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    1: more than one diagonal entry is zero
   .n    2: InvertSmallBlock failed after regularizing
   D*/
/****************************************************************************/

INT l_lrregularize (GRID *theGrid, const MATDATA_DESC *M)
{
  INT type,found,ncmp,i,l,cmp,singComp,matComp;
  DOUBLE value, min, InvMat[MAX_SINGLE_MAT_COMP];

  /* find singular component and regularize */
  type = VTYPE(LASTVECTOR(theGrid));
  ncmp = MD_ROWS_IN_RT_CT(M,type,type);
  found = 0;
  min = MAX_D;
  for (i=0; i<ncmp; i++)
  {
    cmp   = MD_MCMP_OF_RT_CT(M,type,type,i*ncmp+i);
    value = MVALUE(VSTART(LASTVECTOR(theGrid)),cmp);
    if (fabs(value)<10*SMALL_D)
    {
      found++;
      singComp = i;
      matComp  = cmp;
    }
    if (fabs(value)<min)
    {
      min = fabs(value);
      singComp = i;
      matComp  = cmp;
    }
  }
  /* hrr 31.3.96: check problematic because of round-up errors */
  if (found>1)
  {
    PrintErrorMessage('E',"l_lrregularize","more than one singular component in last block");
    REP_ERR_RETURN (1);
  }
  MVALUE(VSTART(LASTVECTOR(theGrid)),matComp) = 1.0;

  PRINTDEBUG(np,1,(" - decomposition regularized on level %d, component %d\n",(int)GLEVEL(theGrid),singComp));

  if (StoreInverse)
  {
    if (InvertSmallBlock(ncmp,MD_MCMPPTR_OF_RT_CT(M,type,type),MVALUEPTR(VSTART(LASTVECTOR(theGrid)),0),InvMat)!=0)
      REP_ERR_RETURN (2);                               /* decompostion failed */

    /* write inverse back to diagonal block */
    for (l=0; l<ncmp*ncmp; l++)
      MVALUEPTR(VSTART(LASTVECTOR(theGrid)),0)[MD_MCMP_OF_RT_CT(M,type,type,l)] = InvMat[l];
  }

  return (0);
}

/****************************************************************************/
/*D
   l_lrregularizeB - regularize last diagonal block of BLOCKVECTOR

   SYNOPSIS:
   INT l_lrregularizeB (GRID *theGrid, const MATDATA_DESC *M)

   PARAMETERS:
   .  theGrid - pointer to grid
   .  M - type matrix descriptor for decomposition

   DESCRIPTION:
   This function regularizes the last diagonal block if it was found to be singular.
   It checks wether it is singular exactly for one component and replaces the zero by 1.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    1: more than one diagonal entry is zero
   .n    2: InvertSmallBlock failed after regularizing
   D*/
/****************************************************************************/

static INT l_lrregularizeB (GRID *theGrid, VECTOR *vec, const MATDATA_DESC *M)
{
  INT type,found,ncmp,i,l,cmp,singComp,matComp;
  DOUBLE value, InvMat[MAX_SINGLE_MAT_COMP],min;

  /* find singular component and regularize */
  type = VTYPE(vec);
  ncmp = MD_ROWS_IN_RT_CT(M,type,type);
  found = 0;
  for (i=0; i<ncmp; i++)
  {
    cmp   = MD_MCMP_OF_RT_CT(M,type,type,i*ncmp+i);
    value = MVALUE(VSTART(vec),cmp);
    if (fabs(value)<SMALL_DET)
    {
      found++;
      singComp = i;
      matComp  = cmp;
    }
  }
  min = MAX_D;
  if (found!=1)
    /* just take smallest diag component */
    for (i=0; i<ncmp; i++)
    {
      cmp   = MD_MCMP_OF_RT_CT(M,type,type,i*ncmp+i);
      value = MVALUE(VSTART(vec),cmp);
      if (fabs(value)<min)
      {
        min = fabs(value);
        singComp = i;
        matComp  = cmp;
      }
    }
  MVALUE(VSTART(vec),matComp) = 1.0;

  if (TRUE /*mutelevel<=VERBOSE_SMOOTH*/)
    UserWriteF(" - BLOCKVECTOR decomposition regularized on level %d, component %d\n",(int)GLEVEL(theGrid),singComp);

  if (StoreInverse)
  {
    if (InvertSmallBlock(ncmp,MD_MCMPPTR_OF_RT_CT(M,type,type),MVALUEPTR(VSTART(vec),0),InvMat)!=0)
      REP_ERR_RETURN (2);                               /* decompostion failed */

    /* write inverse back to diagonal block */
    for (l=0; l<ncmp*ncmp; l++)
      MVALUEPTR(VSTART(vec),0)[MD_MCMP_OF_RT_CT(M,type,type,l)] = InvMat[l];
  }

  return (0);
}

/****************************************************************************/
/*D
   l_lrdecompB - compute left-right decomposition of M BLOCKVECTOR-wise

   SYNOPSIS:
   INT l_lrdecompB (GRID *g, const MATDATA_DESC *M, const DOUBLE *threshold);

   PARAMETERS:
   .  g - pointer to grid
   .  M - type matrix descriptor for precondition

   DESCRIPTION:
   This function computes left-right decomposition of M.
   The matrix M is overwritten!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_OUT_OF_MEM if there is not enough memory to store the decomposition
   .n    i<0 if the decomposition fails at the 'VECTOR' with 'INDEX' '-i'
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_lrdecompB (GRID *g, const MATDATA_DESC *M)
{
  VECTOR *vi,*vj,*vk,*vec;
  BLOCKVECTOR *theBV;
  MATRIX *Mij,*Mji,*Mjk,*Mik;
  DOUBLE InvMat[MAX_SINGLE_MAT_COMP],PivMat[MAX_SINGLE_MAT_COMP];
  DOUBLE CorMat[MAX_SINGLE_MAT_COMP];
  DOUBLE sum;
  register DOUBLE *Diag,*Piv,*Elm,*Mat;
  register SHORT *DiagComp,*PivComp,*ElmComp,*MatComp;
  register INT i0,j0,k0;
  INT type,ctype,rtype,PivIsZero,CorIsZero;
  INT n,n2,nr,nnr,nc,nrnc;
  INT i,l,mc,mask,index,maxBVmembers,bvn;
  DOUBLE diag,invdiag,pivot;

  /* consistency check: diagonal blocks are supposed to be square matrices */
  for (type=0; type<NVECTYPES; type++)
    if (MD_ROWS_IN_RT_CT(M,type,type)>0)
    {
      nr = MD_ROWS_IN_RT_CT(M,type,type);
      ASSERT (nr*nr <= MAX_SINGLE_MAT_COMP);
      /* if too little: increase MAX_SINGLE_VEC_COMP and recompile */
                        #ifdef NDEBUG
      if (nr*nr > MAX_SINGLE_MAT_COMP)
        /* check also in case NDEBUG is defined (assert off)	*/
        REP_ERR_RETURN (__LINE__);
                        #endif
      if (nr != MD_COLS_IN_RT_CT(M,type,type))
        REP_ERR_RETURN (__LINE__);
    }

  /* consistency check: the transpose block-matrices (iff) must have
     the same format */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype+1; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        if (MD_ROWS_IN_RT_CT(M,rtype,rtype)!=MD_ROWS_IN_RT_CT(M,rtype,ctype))
          REP_ERR_RETURN (__LINE__);
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)!=MD_COLS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
        if (MD_COLS_IN_RT_CT(M,rtype,ctype)!=MD_ROWS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
      }

  /* consistency check: the block pattern should allow for elimination */
  for (type=0; type<NVECTYPES; type++)
    for (rtype=type+1; rtype<NVECTYPES; rtype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,type)>0)
        for (ctype=type+1; ctype<NVECTYPES; ctype++)
          if ((MD_ROWS_IN_RT_CT(M,type,ctype)>0) &&
              (MD_ROWS_IN_RT_CT(M,rtype,ctype)==0))
            REP_ERR_RETURN (__LINE__);

  /* consistency check: BLOCKVECTOR-decomposition should consist of only one level */
  maxBVmembers = NVEC(g);
  for (theBV=GFIRSTBV(g); theBV!=NULL; theBV=BVSUCC(theBV))
    if (BVDOWNTYPE(theBV)!=BVDOWNTYPEVECTOR)
      REP_ERR_RETURN (__LINE__)
      else
      {
        /* encode associated BLOCKVECTOR number in vector index */
        /* this index numbering is consistent with what other iterative schemes need */
        index = BVNUMBER(theBV)*maxBVmembers;
        for (vec=BVFIRSTVECTOR(theBV); vec!=BVENDVECTOR(theBV); vec=SUCCVC(vec))
          VINDEX(vec) = index++;
      }

  /* decompose the matrix */
  if (MD_IS_SCALAR(M))
  {
    mc   = MD_SCALCMP(M);
    mask = 0;
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
        mask |= 1<<type;

    /* loop over all lines */
    for (theBV=GFIRSTBV(g); theBV!=NULL; theBV=BVSUCC(theBV))
    {
      bvn = BVNUMBER(theBV);
      for (vi=BVFIRSTVECTOR(theBV); vi!=BVENDVECTOR(theBV); vi=SUCCVC(vi))
      {
        /* check class and mask */
        if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;
        i = VINDEX(vi);

        /* now we are at line i */
        diag = MVALUE(VSTART(vi),mc);                                           /* diagonal element */
        /* decomposition failed */
        if (fabs(diag)<SMALL_DET)
        {
          if (vi!=((BVENDVECTOR(theBV)==NULL) ? LASTVECTOR(g) : PREDVC(BVENDVECTOR(theBV))))
            REP_ERR_RETURN (-i);
          diag = MVALUE(VSTART(vi),mc) = 1.0;
          UserWriteF("block %d regularized, vector %d, component %d\n",bvn,i,mc);
        }

        /* store inverse back to diag */
        if (StoreInverse)
          MVALUE(VSTART(vi),mc) = invdiag = 1.0/diag;
        else
          invdiag = 1.0/diag;

        /* eliminate all entries (j,i) with j>i */
        for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
        {
          /* check class and mask */
          vj = MDEST(Mij);
          if (V_BVNUMBER(vj,maxBVmembers)!=bvn) continue;
          if ( !((VDATATYPE(vj)&mask)&&(VCLASS(vj)>=ACTIVE_CLASS)&&(VINDEX(vj)>i)) ) continue;
          Mji = MADJ(Mij);
          pivot = MVALUE(Mji,mc)*invdiag;

          /* the pivot becomes the entry of the lower triangular part */
          MVALUE(Mji,mc) = pivot;

          if (pivot==0.0)
            continue;                                                   /* nothing to eliminate */

          /* do gaussian elimination on the pattern (all ujk, k>i) */
          for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
          {
            vk = MDEST(Mik);
            if (V_BVNUMBER(vk,maxBVmembers)!=bvn) continue;
            if ( !((VDATATYPE(vk)&mask)&&(VCLASS(vk)>=ACTIVE_CLASS)&&(VINDEX(vk)>i)) ) continue;
            Mjk = GetMatrix(vj,vk);
            if (Mjk==NULL)
            {
              Mjk = CreateExtraConnection(g,vj,vk);
              if (Mjk==NULL)
                REP_ERR_RETURN (NUM_OUT_OF_MEM);
            }
            MVALUE(Mjk,mc) -= pivot*MVALUE(Mik,mc);                                                                       /* entry is in pattern */
          }
        }
      }
    }
    return (NUM_OK);
  }

  /* loop over all lines */
  for (theBV=GFIRSTBV(g); theBV!=NULL; theBV=BVSUCC(theBV))
  {
    bvn = BVNUMBER(theBV);
    for (vi=BVFIRSTVECTOR(theBV); vi!=BVENDVECTOR(theBV); vi=SUCCVC(vi))
    {
      if (VCLASS(vi)<ACTIVE_CLASS) continue;
      type = VTYPE(vi);
      n        = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      DiagComp = MD_MCMPPTR_OF_RT_CT(M,type,type);
      n2       = n*n;

      i = VINDEX(vi);

      Diag = MVALUEPTR(VSTART(vi),0);

      /* decompostion failed */
      if (InvertSmallBlock(n,DiagComp,Diag,InvMat)!=0)
        if (l_lrregularizeB(g,(BVENDVECTOR(theBV)==NULL) ? LASTVECTOR(g) : PREDVC(BVENDVECTOR(theBV)),M)!=0)
          REP_ERR_RETURN (-i);
      /*if (InvertSmallBlock(n,DiagComp,Diag,InvMat)!=0)
         REP_ERR_RETURN (-i); */

      /* write inverse back to diagonal block */
      if (StoreInverse)
        for (l=0; l<n2; l++)
          Diag[DiagComp[l]] = InvMat[l];

      /* eliminate all entries (j,i) with j>i */
      for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
      {
        rtype = VTYPE(vj=MDEST(Mij));

        if (V_BVNUMBER(vj,maxBVmembers)!=bvn) continue;
        if (!((MD_ROWS_IN_RT_CT(M,rtype,type)>0)
              && (VCLASS(vj)>=ACTIVE_CLASS)
              && (i<VINDEX(vj))))
          continue;

        PivComp = MD_MCMPPTR_OF_RT_CT(M,rtype,type);
        nr              = MD_ROWS_IN_RT_CT(M,rtype,type);
        nnr             = nr*n;

        Mji = MADJ(Mij);
        Piv = MVALUEPTR(Mji,0);

        /* matrix multiplication */
        PivIsZero = TRUE;
        for (i0=0; i0<nr; i0++)
          for (j0=0; j0<n; j0++)
          {
            sum = 0.0;
            for (k0=0; k0<n; k0++)
              sum += Piv[PivComp[i0*n+k0]] * InvMat[k0*n+j0];
            PivMat[i0*n+j0] = sum;
            if (sum!=0.0)
              PivIsZero = FALSE;
          }

        /* the pivot becomes the entry of the lower triangular part */
        for (l=0; l<nnr; l++)
          Piv[PivComp[l]] = PivMat[l];

        if (PivIsZero)
          continue;                                             /* nothing to eliminate */

        /* do gaussian elimination on the pattern (all Mjk, k>i) */
        for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
        {
          ctype   = VTYPE(vk=MDEST(Mik));

          if (V_BVNUMBER(vk,maxBVmembers)!=bvn) continue;
          if (!((MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
                && (VCLASS(vk)>=ACTIVE_CLASS)
                && (i<VINDEX(vk))))
            continue;

          ElmComp = MD_MCMPPTR_OF_RT_CT(M,type,ctype);
          nc              = MD_COLS_IN_RT_CT(M,type,ctype);
          nrnc    = nr*nc;

          MatComp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);

          Elm = MVALUEPTR(Mik,0);

          /* matrix multiplication */
          CorIsZero = TRUE;
          for (i0=0; i0<nr; i0++)
            for (j0=0; j0<nc; j0++)
            {
              sum = 0.0;
              for (k0=0; k0<n; k0++)
                sum += PivMat[i0*n+k0] * Elm[ElmComp[k0*nc+j0]];
              CorMat[i0*nc+j0] = sum;
              if (sum!=0.0)
                CorIsZero = FALSE;
            }

          if (CorIsZero)
            continue;                                                   /* nothing to correct */

          Mjk = GetMatrix(vj,vk);
          if (Mjk==NULL)
          {
            Mjk = CreateExtraConnection(g,vj,vk);
            if (Mjk==NULL)
              REP_ERR_RETURN (NUM_OUT_OF_MEM);
          }
          /* subtract the entry from Mjk */
          Mat = MVALUEPTR(Mjk,0);
          for (l=0; l<nrnc; l++)
            Mat[MatComp[l]] -= CorMat[l];
        }
      }
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_luiter - solve L*U*v=d

   SYNOPSIS:
   INT l_luiter (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
   const VECDATA_DESC *d);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `L*U*v=d`, where `L`, `U` are the factors from a
   (in-)complete decomposition, stored in 'M'.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    NUM_BLOCK_TOO_LARGE if the blocks are larger as MAX_SINGLE_VEC_COMP
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_luiter (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*last_vec;
  INT rtype,ctype,myindex,err;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*vcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat,*vmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = FIRSTVECTOR(g);
  last_vec  = LASTVECTOR(g);

  if (MD_IS_SCALAR(M) && VD_IS_SCALAR(v) && VD_IS_SCALAR(d))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= NULL; vec=SUCCVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum);                           /* since Diag(L)=I per convention */
      }
    }

    /* solve UpperTriangle(v)=d */
    for (vec=last_vec; vec!= NULL; vec=PREDVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        /* solve (we stored the inverse of the diagonal) */
        if (StoreInverse)
          VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) * MVALUE(VSTART(vec),mc);
        else
          VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) / MVALUE(VSTART(vec),mc);
      }
    }
    return (NUM_OK);
  }

  /* solve lower traingle */
  L_VLOOP__CLASS(vec,first_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n           = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    dcomp   = VD_CMPPTR_OF_TYPE(d,rtype);
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,dcomp[i]);
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
            {
              wmat  = VVALPTR(w);
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
            }
        }

    /* solve (Diag(L)=I per convention) */
    vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
    vmat  = VVALPTR(vec);
    for (i=0; i<n; i++)
      vmat[vcomp[i]] = s[i];
  }

  last_vec = LASTVECTOR(g);

  /* solve upper triangle */
  L_REVERSE_VLOOP__CLASS(vec,last_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,VD_CMP_OF_TYPE(v,rtype,i));
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
            {
              wmat  = VVALPTR(w);
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j])
                          * wmat[wcomp[j]];
            }
        }

    /* solve */
    if (StoreInverse)
    {
      if (SolveInverseSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                                 VVALPTR(vec),
                                 MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                                 MVALPTR(VSTART(vec)),s)!=0)
        REP_ERR_RETURN (__LINE__);
    }
    else
    {
      if (SolveSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                          VVALPTR(vec),
                          MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                          MVALPTR(VSTART(vec)),s)!=0)
        REP_ERR_RETURN (__LINE__);
    }
  }

  return (NUM_OK);
}

INT l_luiter_SB (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*last_vec;
  INT rtype,ctype,myindex,err,first_index,last_index;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*vcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat,*vmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = BVFIRSTVECTOR(theBV);
  last_vec = BVLASTVECTOR(theBV);
  first_index = VINDEX(first_vec);
  last_index = VINDEX(last_vec);

  if (MD_IS_SCALAR(M) && VD_IS_SCALAR(v) && VD_IS_SCALAR(d))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= SUCCVC(last_vec); vec=SUCCVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>=first_index) && (VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum);                           /* since Diag(L)=I per convention */
      }
    }

    /* solve UpperTriangle(v)=d */
    for (vec=last_vec; vec!= PREDVC(first_vec); vec=PREDVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex) && (VINDEX(w)<=last_index) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) / MVALUE(VSTART(vec),mc);
      }
    }
    return (NUM_OK);
  }

  REP_ERR_RETURN (1);
}

INT l_tpluiter_SB (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*last_vec;
  INT rtype,ctype,myindex,err,first_index,last_index;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*vcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat,*vmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = BVFIRSTVECTOR(theBV);
  last_vec = BVLASTVECTOR(theBV);
  first_index = VINDEX(first_vec);
  last_index = VINDEX(last_vec);

  if (MD_IS_SCALAR(M) && VD_IS_SCALAR(v) && VD_IS_SCALAR(d))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= SUCCVC(last_vec); vec=SUCCVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>=first_index) && (VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(MADJ(mat),mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum) / MVALUE(VSTART(vec),mc);
      }
    }

    /* solve UpperTriangle(v)=d */
    for (vec=last_vec; vec!= PREDVC(first_vec); vec=PREDVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex) && (VINDEX(w)<=last_index) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(MADJ(mat),mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,vc)-sum);                          /* since Diag(L)=I per convention */
      }
    }
    return (NUM_OK);
  }

  REP_ERR_RETURN (1);
}

#ifdef __BLOCK_VECTOR_DESC__

/****************************************************************************/
/*D
   solveLUMatBS - solve for a LU-decomposed diagonal block of a matrix

   SYNOPSIS:
   INT solveLUMatBS( const BLOCKVECTOR *bv, const BV_DESC *bvd,
       const BV_DESC_FORMAT *bvdf, INT dest_comp, INT LU_comp, INT source_comp );

   PARAMETERS:
   .  bv - blockvector of the blockmatrix
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  dest_comp - position of the result in the VECTOR-data
   .  LU_comp - position of the LU decomposition in the MATRIX-data of the blockmatrix
   .  source_comp - position of the given right hand side in the VECTOR-data

   DESCRIPTION:
   Given a LU-decposition of a matrixblock 'A' and a right hand side vector
   'source_comp', this function solves 'A * dest_comp = source_comp'. The
   matrixblock must be on the diagonal of the global matrix.

   The VINDEX-field must be ascendend in the VECTOR-chain.

   Diagonalelements < SMALL_D are considered as 0 and the corresponding
   division is refused with an error message.

   The decomposition can be calculated by 'LUDecomposeDiagBS'.

   REMARK:
   dest_comp == source_comp is possible!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_SMALL_DIAG if a diagonal element was < SMALL_D (i.e. division by nearly 0)

   SEE ALSO:
   BLOCKVECTOR, blas_routines, LUDecomposeDiagBS
   D*/
/****************************************************************************/

INT solveLUMatBS( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT dest_comp, INT LU_comp, INT source_comp )
{
  register VECTOR *vi, *end_v, *vj;
  register MATRIX *m;
  register DOUBLE val, diag;
  register INT index_vi;

  ASSERT( bv != NULL );
  ASSERT( dest_comp != -1 );
  ASSERT( LU_comp != -1 );
  ASSERT( source_comp != -1 );

  /* solve lower triangular matrix */
  end_v = BVENDVECTOR( bv );
  vi = BVFIRSTVECTOR( bv );
  ASSERT( vi != NULL );
  /* the first VECTOR needs only copy */
  VVALUE( vi, dest_comp) = VVALUE( vi, source_comp );
  if ( vi != end_v )
    vi = SUCCVC( vi );                  /* be sure not to skip over the end */
  for ( ; vi != end_v; vi = SUCCVC( vi ) )
  {
    val = VVALUE( vi, source_comp );
    index_vi = VINDEX( vi );
    for ( m = VSTART( vi ); m != NULL; m = MNEXT( m ) )
    {
      vj = MDEST( m );
      if ( (VINDEX(vj) < index_vi) && VMATCH( vj, bvd, bvdf ) )
        val -= MVALUE( m, LU_comp ) * VVALUE( vj, dest_comp );
    }
    VVALUE( vi, dest_comp ) = val;
  }

  /* solve upper triangular matrix */
  end_v = PREDVC( BVFIRSTVECTOR( bv ) );
  vi = BVLASTVECTOR( bv );
  m = VSTART( vi );
  ASSERT( MDEST(m) == vi );
  diag = MVALUE( m, LU_comp );
  if ( fabs( diag ) < SMALL_D )
  {
    PrintErrorMessage( 'E', "solveLUMatBS", "Very small diagonal for division" );
    REP_ERR_RETURN (NUM_SMALL_DIAG);
  }
  VVALUE( vi, dest_comp ) = VVALUE( vi, dest_comp ) / diag;
  vi = PREDVC( vi );
  for ( ; vi != end_v; vi = PREDVC( vi ) )
  {
    val = VVALUE( vi, dest_comp );
    index_vi = VINDEX( vi );
    diag = 0.0;
    for ( m = VSTART( vi ); m != NULL; m = MNEXT( m ) )
    {
      vj = MDEST( m );
      if ( (VINDEX(vj) >= index_vi) && VMATCH( vj, bvd, bvdf ) )
        if ( VINDEX(vj) == index_vi )
          diag = MVALUE( m,  LU_comp );
        else
          val -= MVALUE( m, LU_comp ) * VVALUE( vj, dest_comp );
    }

    if ( fabs( diag ) < SMALL_D )
    {
      PrintErrorMessage( 'E', "solveLUMatBS", "Very small diagonal for division" );
      REP_ERR_RETURN (NUM_SMALL_DIAG);
    }

    VVALUE( vi, dest_comp ) = val / diag;
  }

  return NUM_OK;
}
#endif

/****************************************************************************/
/*D
   l_luiterB - solve L*U*v=d

   SYNOPSIS:
   INT l_luiterB (GRID *g, const BLOCKVECTOR *bv,
                const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);

   PARAMETERS:
   .  g  - pointer to grid
   .  bv - block vector to invert
   .  v  - type vector descriptor to store correction
   .  M  - type matrix descriptor for precondition
   .  d  - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `L*U*v=d`, where `L`, `U` are the factors from a
   (in-)complete decomposition, stored in 'M'.
   Only the block 'bv' is solved. Without updating any defect.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    NUM_BLOCK_TOO_LARGE if the blocks are larger as MAX_SINGLE_VEC_COMP
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_luiterB (GRID *g, const BLOCKVECTOR *bv, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*last_vec;
  INT rtype,ctype,myindex,err,bvn,maxBVmembers;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*vcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat,*vmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  maxBVmembers = NVEC(g);
  bvn = BVNUMBER(bv);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=BVFIRSTVECTOR(bv); vec!=BVENDVECTOR(bv); vec=SUCCVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if (V_BVNUMBER(w,maxBVmembers)!=bvn) continue;
          if ((VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum);                           /* since Diag(L)=I per convention */
      }
    }

    /* solve UpperTriangle(v)=d */
    first_vec = (BVFIRSTVECTOR(bv)==FIRSTVECTOR(g)) ? NULL : PREDVC(BVFIRSTVECTOR(bv));
    last_vec  = (BVENDVECTOR(bv)==NULL) ? LASTVECTOR(g) : PREDVC(BVENDVECTOR(bv));
    for (vec=last_vec; vec!=first_vec; vec=PREDVC(vec))
    {
      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if (V_BVNUMBER(w,maxBVmembers)!=bvn) continue;
          if ((VINDEX(w)>myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        /* solve (we stored the inverse of the diagonal) */
        if (StoreInverse)
          VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) * MVALUE(VSTART(vec),mc);
        else
          VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) / MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  /* solve lower traingle */
  for (vec=BVFIRSTVECTOR(bv); vec!=BVENDVECTOR(bv); vec=SUCCVC(vec))
  {
    if (VCLASS(vec)<ACTIVE_CLASS) continue;

    rtype = VTYPE(vec);

    n           = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    dcomp   = VD_CMPPTR_OF_TYPE(d,rtype);
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,dcomp[i]);
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
            {
              wmat  = VVALPTR(w);
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
            }
        }

    /* solve (Diag(L)=I per convention) */
    vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
    vmat  = VVALPTR(vec);
    for (i=0; i<n; i++)
      vmat[vcomp[i]] = s[i];
  }

  /* solve upper triangle */
  first_vec = (BVFIRSTVECTOR(bv)==FIRSTVECTOR(g)) ? NULL : PREDVC(BVFIRSTVECTOR(bv));
  last_vec  = (BVENDVECTOR(bv)==NULL) ? LASTVECTOR(g) : PREDVC(BVENDVECTOR(bv));
  for (vec=last_vec; vec!=first_vec; vec=PREDVC(vec))
  {
    if (VCLASS(vec)<ACTIVE_CLASS) continue;

    rtype = VTYPE(vec);

    n = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,VD_CMP_OF_TYPE(v,rtype,i));
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (V_BVNUMBER(w,maxBVmembers)==bvn) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
            {
              wmat  = VVALPTR(w);
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j])
                          * wmat[wcomp[j]];
            }
        }


    /* solve */
    if (StoreInverse)
    {
      if (SolveInverseSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                                 VVALPTR(vec),
                                 MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                                 MVALPTR(VSTART(vec)),s)!=0)
        REP_ERR_RETURN (__LINE__);
    }
    else
    {
      if (SolveSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                          VVALPTR(vec),
                          MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                          MVALPTR(VSTART(vec)),s)!=0)
        REP_ERR_RETURN (__LINE__);
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_lltiter - solve L * L^T * v = d

   SYNOPSIS:
   INT l_lltiter (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
   const VECDATA_DESC *d);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `L * L^T * v = d`, where `L` is a factor from a
   (in-)complete decomposition, stored in 'M'.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    NUM_BLOCK_TOO_LARGE if the blocks are larger as MAX_SINGLE_VEC_COMP
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_lltiter (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*last_vec;
  INT rtype,ctype,myindex,err;
  register MATRIX *mat,*dmat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc,nr;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = FIRSTVECTOR(g);
  last_vec  = LASTVECTOR(g);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= NULL; vec=SUCCVC(vec))
      if ( (VDATATYPE(vec)&mask)
           && (VCLASS(vec)>=ACTIVE_CLASS)
           && (VECSKIP(vec)==0) )
      {
        myindex = VINDEX(vec);
        sum = 0.0;
        dmat = VSTART(vec);
        for (mat=MNEXT(dmat); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)<myindex)
              && (VDATATYPE(w)&mask)
              && (VCLASS(w)>=ACTIVE_CLASS)
              && (VECSKIP(w)==0) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum) * MVALUE(dmat,mc);
      }

    /* solve LowerTriangleTransposed(v)=d */
    for (vec=last_vec; vec!= NULL; vec=PREDVC(vec))
      if ( (VDATATYPE(vec)&mask)
           && (VCLASS(vec)>=ACTIVE_CLASS)
           && (VECSKIP(vec)==0) )
      {
        myindex = VINDEX(vec);
        sum = 0.0;
        dmat = VSTART(vec);
        for (mat=MNEXT(dmat); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex)
              && (VDATATYPE(w)&mask)
              && (VCLASS(w)>=ACTIVE_CLASS)
              && (VECSKIP(w)==0) )
            sum += MVALUE(MADJ(mat),mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) * MVALUE(dmat,mc);
      }

    return (NUM_OK);
  }

  /* solve lower traingle */
  L_VLOOP__CLASS(vec,first_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n           = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0)
      continue;
    if (VECSKIP(vec))
      continue;
    dcomp   = VD_CMPPTR_OF_TYPE(d,rtype);
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,dcomp[i]);
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype)
                 && (VCLASS(w)>=ACTIVE_CLASS))
                && (VECSKIP(w) == 0)
                && (myindex>VINDEX(w)))
            {
              wmat  = VVALPTR(w);
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
            }
        }

    /* solve */
    if (SolveInverseSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                               VVALPTR(vec),
                               MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                               MVALPTR(VSTART(vec)),s)!=0)
      REP_ERR_RETURN (__LINE__);
  }

  last_vec = LASTVECTOR(g);

  /* solve upper triangle */
  L_REVERSE_VLOOP__CLASS(vec,last_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0)
      continue;
    if (VECSKIP(vec))
      continue;
    myindex = VINDEX(vec);

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,VD_CMP_OF_TYPE(v,rtype,i));
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        mcomp = MD_MCMPPTR_OF_RT_CT(M,ctype,rtype);
        wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
        nr    = MD_ROWS_IN_RT_CT(M,ctype,rtype);
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
          if (((VTYPE(w=MDEST(mat))==ctype)
               && (VCLASS(w)>=ACTIVE_CLASS))
              && (VECSKIP(w) == 0)
              && (myindex<VINDEX(w)))
          {
            wmat  = VVALPTR(w);
            for (i=0; i<n; i++)
              for (j=0; j<nr; j++)
                s[i] -= MVALUE(MADJ(mat),mcomp[j*n+i])
                        * wmat[wcomp[j]];
          }
      }

    /* solve */
    if (SolveInverseSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                               VVALPTR(vec),
                               MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                               MVALPTR(VSTART(vec)),s)!=0)
      REP_ERR_RETURN (__LINE__);
  }

  return (NUM_OK);
}



/****************************************************************************/
/*D
   l_ilubthdecomp_fine	- compute incomplete decomposition on fine grid nodes

   SYNOPSIS:
   INT l_ilubthdecomp_fine (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta,
   const VEC_SCALAR threshold, const VECDATA_DESC *rest,
   const VEC_SCALAR oldrestthresh);

   PARAMETERS:
   .  g - pointer to grid
   .  M - type matrix descriptor
   .  beta - modification parameters
   .  threshold - introduce new connection if entry > threshold
   .  rest - if !=NULL store normalized rest matrix row sums here
   .  oldrestthresh - if !=NULL use rest field of last step to flag vectors where
                                                connections will be introduced

   DESCRIPTION:
   This function computes an incomplete decomposition of order 0
   with modification of the diagonal.
   A new connection is introduced if an entry is larger than threshold.
   THIS VERSION SKIPS VECTORS THAT BELONG TO COARSE GRID NODES.

   The parameters beta and threshold are not considered if they are 'NULL'.

   The matrix M is overwritten!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    i<0 if decomposition failed
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_ilubthdecomp_fine (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VEC_SCALAR threshold, const VECDATA_DESC *VD_rest, const VEC_SCALAR oldrestthresh)
{
  VECTOR *vi,*vj,*vk;
  MATRIX *Mij,*Mji,*Mjk,*Mik;
  DOUBLE InvMat[MAX_SINGLE_MAT_COMP],PivMat[MAX_SINGLE_MAT_COMP];
  DOUBLE CorMat[MAX_SINGLE_MAT_COMP];
  DOUBLE sum;
  INT offset[NVECTYPES+1];
  register DOUBLE *Diag,*Piv,*Elm,*Mat,*Djj,*Dkk;
  register SHORT *DiagComp,*PivComp,*ElmComp,*MatComp,*DjjComp,*DkkComp;
  register INT i0,j0,k0,l,m;
  INT type,ctype,rtype,PivIsZero,CorIsZero;
  INT n,n2,nr,nnr,nc,nrnc;
  INT i,mc,mask;
  DOUBLE RowSum[MAX_SINGLE_VEC_COMP],Damp[MAX_SINGLE_VEC_COMP];
  VEC_SCALAR j_Normalization,k_Normalization;
  DOUBLE diag,invdiag,pivot,AbsDjj;
  const DOUBLE *TypeBeta,*TypeThresh;
  const DOUBLE *TypeORT;
  DOUBLE *Rest;
  SHORT *RestComp;
  INT flag;

  /* consistency check: diagonal blocks are supposed to be square matrices */
  for (type=0; type<NVECTYPES; type++)
    if (MD_ROWS_IN_RT_CT(M,type,type)>0)
    {
      nr = MD_ROWS_IN_RT_CT(M,type,type);

      ASSERT (nr*nr <= MAX_SINGLE_MAT_COMP);
      /* if too little: increase MAX_SINGLE_VEC_COMP and recompile	*/
                        #ifdef NDEBUG
      if (nr*nr > MAX_SINGLE_MAT_COMP)
        /* check also in case NDEBUG is defined (assert off)	*/
        REP_ERR_RETURN (__LINE__);
                        #endif
      if (nr != MD_COLS_IN_RT_CT(M,type,type))
        REP_ERR_RETURN (__LINE__);
    }
  /* check VD_rest iff */
  if (VD_rest!=NULL)
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
      {
        if (VD_NCMPS_IN_TYPE(VD_rest,type)==0)
          REP_ERR_RETURN (__LINE__);
        if (VD_NCMPS_IN_TYPE(VD_rest,type)!=MD_ROWS_IN_RT_CT(M,type,type))
          REP_ERR_RETURN (__LINE__);
      }

  /* consistency check:
     the transpose block-matrices (iff) must have the same format */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype+1; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
      {
        if (MD_ROWS_IN_RT_CT(M,rtype,rtype)!=MD_ROWS_IN_RT_CT(M,rtype,ctype))
          REP_ERR_RETURN (__LINE__);
        if (MD_ROWS_IN_RT_CT(M,rtype,ctype)!=MD_COLS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
        if (MD_COLS_IN_RT_CT(M,rtype,ctype)!=MD_ROWS_IN_RT_CT(M,ctype,rtype))
          REP_ERR_RETURN (__LINE__);
      }

  /* calculate offsets for beta and threshold for different types */
  offset[0] = 0;
  for (type=1; type<NVECTYPES; type++)
    offset[type] = offset[type-1] + MD_ROWS_IN_RT_CT(M,type-1,type-1);

  /* flag vectors where connections will be introduced */
  if ((VD_rest!=NULL) && (oldrestthresh!=NULL))
  {
    L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
    {
      type     = VTYPE(vi);
      n        = MD_ROWS_IN_RT_CT(M,type,type);
      if (n == 0) continue;
      TypeORT  = oldrestthresh+offset[type];
      Rest     = VVALUEPTR(vi,0);
      RestComp = VD_CMPPTR_OF_TYPE(VD_rest,type);

      flag = FALSE;
      for (i=0; i<n; i++)
        if (Rest[RestComp[i]]>TypeORT[i])
          flag = TRUE;

      SETVCUSED(vi,flag);
    }
  }
  else
    L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
    SETVCUSED(vi,0);

  /* clear VD_rest iff */
  if (VD_rest!=NULL)
    l_dset(g,VD_rest,EVERY_CLASS,0.0);

  /* decompose the matrix */
  if (MD_IS_SCALAR(M))
  {
    mc = MD_SCALCMP(M);
    mask = 0;
    for (type=0; type<NVECTYPES; type++)
      if (MD_ROWS_IN_RT_CT(M,type,type)>0)
        mask |= 1<<type;

    /* loop over all lines */
    for (vi=FIRSTVECTOR(g); vi!=NULL; vi=SUCCVC(vi))
    {
      /* check class and mask */
      if ( !((VDATATYPE(vi)&mask)&&(VCLASS(vi)>=ACTIVE_CLASS)) ) continue;
      i = VINDEX(vi);

      /* check coarse grid position */
      if (VTYPE(vi)==NODEVECTOR)
        if (CORNERTYPE(VMYNODE(vi))) continue;                         /* skip coarse grid node */

      /* now we are at line i */
      diag = MVALUE(VSTART(vi),mc);                                     /* diagonal element */
      if (fabs(diag)<SMALL_D) REP_ERR_RETURN(-i);                               /* decomposition failed */

      /* store inverse back to diag */
      if (StoreInverse)
        MVALUE(VSTART(vi),mc) = invdiag = 1.0/diag;
      else
        invdiag = 1.0/diag;

      /* eliminate all entries (j,i) with j>i */
      for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
      {
        /* check class and mask */
        vj = MDEST(Mij);
        if ( !((VDATATYPE(vj)&mask)&&(VCLASS(vj)>=ACTIVE_CLASS)&&(VINDEX(vj)>i)) ) continue;
        Mji = MADJ(Mij);
        AbsDjj = fabs(MVALUE(VSTART(vj),mc));
        pivot = MVALUE(Mji,mc)*invdiag;

        /* the pivot becomes the entry of the lower triangular part */
        MVALUE(Mji,mc) = pivot;

        if (pivot==0.0)
          continue;                                             /* nothing to eliminate */

        /* do gaussian elimination on the pattern (all ujk, k>i) */
        for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
        {
          vk = MDEST(Mik);
          if ( !((VDATATYPE(vk)&mask)&&(VCLASS(vk)>=ACTIVE_CLASS)&&(VINDEX(vk)>i)) ) continue;
          Mjk = GetMatrix(vj,vk);
          if (threshold!=NULL)                                  /* only if threshold is defined */
            if (Mjk==NULL)
              /* the threshold is taken relative to the diagonal of row j */
              if (fabs(pivot*MVALUE(Mik,mc))>threshold[0]*AbsDjj)
              {
                /* introduce new connection */
                Mjk = CreateExtraConnection(g,vj,vk);
                if (Mjk==NULL)
                  REP_ERR_RETURN (NUM_OUT_OF_MEM);
              }
          if (Mjk==NULL)
            if (VCUSED(vj))
            {
              /* introduce new connection */
              Mjk = CreateExtraConnection(g,vj,vk);
              if (Mjk==NULL)
                REP_ERR_RETURN (NUM_OUT_OF_MEM);
            }
          if (Mjk!=NULL)
            MVALUE(Mjk,mc) -= pivot*MVALUE(Mik,mc);                                                                       /* entry is in pattern */
          else
          {
            if (beta!=NULL)                                     /* only if beta is defined */
              MVALUE(VSTART(vj),mc) += beta[0]*fabs(pivot*MVALUE(Mik,mc));                                           /* entry not in pattern */
            if (VD_rest!=NULL)
              VVALUE(vj,VD_CMP_OF_TYPE(VD_rest,VTYPE(vj),0)) += fabs(pivot*MVALUE(Mik,mc));
          }
        }
      }
    }
    return (NUM_OK);
  }

  /* loop over all lines */
  L_VLOOP__CLASS(vi,FIRSTVECTOR(g),ACTIVE_CLASS)
  {
    type     = VTYPE(vi);
    n        = MD_ROWS_IN_RT_CT(M,type,type);
    if (n == 0) continue;
    DiagComp = MD_MCMPPTR_OF_RT_CT(M,type,type);
    n2       = n*n;

    /* check coarse grid position */
    if (VTYPE(vi)==NODEVECTOR)
      if (CORNERTYPE(VMYNODE(vi))) continue;                   /* skip coarse grid node */

    i = VINDEX(vi);

    Diag = MVALUEPTR(VSTART(vi),0);

    if (InvertSmallBlock(n,DiagComp,Diag,InvMat)!=0)
      REP_ERR_RETURN (-i);                              /* decompostion failed */

    /* write inverse back to diagonal block */
    if (StoreInverse)
      for (l=0; l<n2; l++)
        Diag[DiagComp[l]] = InvMat[l];

    /* eliminate all entries (j,i) with j>i */
    for (Mij=MNEXT(VSTART(vi)); Mij!=NULL; Mij=MNEXT(Mij))
    {
      rtype   = VTYPE(vj=MDEST(Mij));

      if (!((MD_ROWS_IN_RT_CT(M,rtype,type)>0)
            && (VCLASS(vj)>=ACTIVE_CLASS)
            && (i<VINDEX(vj))))
        continue;

      DjjComp = MD_MCMPPTR_OF_RT_CT(M,rtype,rtype);
      Djj             = MVALUEPTR(VSTART(vj),0);
      PivComp = MD_MCMPPTR_OF_RT_CT(M,rtype,type);
      nr              = MD_ROWS_IN_RT_CT(M,rtype,type);
      nnr             = nr*n;
      if (threshold!=NULL)                      /* only if threshold is defined */
        TypeThresh = threshold+offset[rtype];

      /* the normalization factors of the j diagonal */
      for (l=0; l<nr; l++)
        j_Normalization[l] = 1.0 / sqrt(fabs(Djj[DjjComp[l*nr+l]]));

      /* and we use a further vector to store the row sums of the normalized rest matrix */
      for (l=0; l<nr; l++) RowSum[l] = 0.0;

      Mji = MADJ(Mij);
      Piv = MVALUEPTR(Mji,0);

      /* matrix multiplication */
      PivIsZero = TRUE;
      for (i0=0; i0<nr; i0++)
        for (j0=0; j0<n; j0++)
        {
          sum = 0.0;
          for (k0=0; k0<n; k0++)
            sum += Piv[PivComp[i0*n+k0]] * InvMat[k0*n+j0];
          PivMat[i0*n+j0] = sum;
          if (sum!=0.0)
            PivIsZero = FALSE;
        }

      /* the pivot becomes the entry of the lower triangular part */
      for (l=0; l<nnr; l++)
        Piv[PivComp[l]] = PivMat[l];

      if (PivIsZero)
        continue;                                       /* nothing to eliminate */

      /* do gaussian elimination on the pattern (all Mjk, k>i) */
      for (Mik=MNEXT(VSTART(vi)); Mik!=NULL; Mik=MNEXT(Mik))
      {
        ctype   = VTYPE(vk=MDEST(Mik));

        if (!((MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
              && (VCLASS(vk)>=ACTIVE_CLASS)
              && (i<VINDEX(vk))))
          continue;

        DkkComp = MD_MCMPPTR_OF_RT_CT(M,ctype,ctype);
        Dkk             = MVALUEPTR(VSTART(vk),0);
        ElmComp = MD_MCMPPTR_OF_RT_CT(M,type,ctype);
        nc              = MD_COLS_IN_RT_CT(M,type,ctype);
        nrnc    = nr*nc;

        MatComp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);

        Elm = MVALUEPTR(Mik,0);

        /* matrix multiplication */
        CorIsZero = TRUE;
        for (i0=0; i0<nr; i0++)
          for (j0=0; j0<nc; j0++)
          {
            sum = 0.0;
            for (k0=0; k0<n; k0++)
              sum += PivMat[i0*n+k0] * Elm[ElmComp[k0*nc+j0]];
            CorMat[i0*nc+j0] = sum;
            if (sum!=0.0)
              CorIsZero = FALSE;
          }

        if (CorIsZero)
          continue;                                             /* nothing to correct */

        /* the normalization factors of the k diagonal */
        if (ctype==rtype)
          /* normalize with diag of vj */
          for (l=0; l<nr; l++)
            k_Normalization[l] = j_Normalization[l];                                            /* this choice led to good results -*/
        /* but why not else in all cases?	*/
        else
          /* normalize with diag of vk */
          for (l=0; l<nc; l++)
            k_Normalization[l] = 1.0 / sqrt(fabs(Dkk[DkkComp[l*nc+l]]));


        Mjk = GetMatrix(vj,vk);
        if (threshold!=NULL)                            /* only if threshold is defined */
          if (Mjk==NULL)
          {
            /* check threshold vs row sums of rest matrix */
            for (l=0; l<nr; l++)
            {
              sum = 0.0;
              for (m=0; m<nc; m++)
                sum += fabs(CorMat[l*nc+m]*j_Normalization[l]*k_Normalization[m]);
              if (sum>TypeThresh[l])
                if (CreateExtraConnection(g,vj,vk)==NULL)
                  REP_ERR_RETURN (NUM_OUT_OF_MEM)
                  else
                    break;
            }
          }
        if (Mjk==NULL)
          if (VCUSED(vj))
          {
            /* introduce new connection */
            Mjk = CreateExtraConnection(g,vj,vk);
            if (Mjk==NULL)
              REP_ERR_RETURN (NUM_OUT_OF_MEM);
          }
        if (Mjk!=NULL)
        {
          /* we are on the pattern: subtract the entry from Mjk */
          Mat = MVALUEPTR(Mjk,0);
          for (l=0; l<nrnc; l++)
            Mat[MatComp[l]] -= CorMat[l];
        }
        else
        {
          /* we are off the pattern: add rowsum */
          for (l=0; l<nr; l++)
            for (m=0; m<nc; m++)
              RowSum[l] += fabs(CorMat[l*nc+m]*j_Normalization[l]*k_Normalization[m]);
        }
      }
      if (beta==NULL) continue;
      /* only if beta is defined */

      /* add to VD_rest iff */
      if (VD_rest!=NULL)
        for (m=0; m<n; m++)
          VVALUE(vj,VD_CMP_OF_TYPE(VD_rest,type,m)) += RowSum[m];

      /* finally we modify the diagonal Djj */
      /* NB: the diag of row j is inverted later, so we don't have to care about StoreInverse */

      /* the diagonal damp matrix */
      TypeBeta= beta+offset[type];
      for (m=0; m<n; m++)
        Damp[m] = 1.0 + TypeBeta[m]*RowSum[m];

      /* Djj = Djj*Damp */
      for (m=0; m<n; m++)
        for (l=0; l<n; l++)
          Djj[DjjComp[m*n+l]] *= Damp[l];
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_luiter_fine - solve L*U*v=d

   SYNOPSIS:
   INT l_luiter_fine (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
   const VECDATA_DESC *d);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)

   DESCRIPTION:
   This function solves `L*U*v=d`, where `L`, `U` are the factors from a
   (in-)complete decomposition, stored in 'M'.
   THIS VERSION SKIPS VECTORS IN COARSE GRID NODES!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_DESC_MISMATCH if the type descriptors not match
   .n    NUM_BLOCK_TOO_LARGE if the blocks are larger as MAX_SINGLE_VEC_COMP
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_luiter_fine (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec,*last_vec;
  INT rtype,ctype,myindex,err;
  register MATRIX *mat;
  register SHORT vc,dc,mc,mask;
  register SHORT *mcomp,*vcomp,*wcomp,*dcomp;
  register SHORT i,j;
  register SHORT n,nc;
  register DOUBLE sum;
  DOUBLE s[MAX_SINGLE_VEC_COMP],*wmat,*vmat;
  DEFINE_VS_CMPS(s);
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  register SHORT *tmpptr;

#ifndef NDEBUG
  if ( (err = MatmulCheckConsistency(v,M,d)) != NUM_OK )
    REP_ERR_RETURN (err);
#endif

  first_vec = FIRSTVECTOR(g);
  last_vec  = LASTVECTOR(g);

  if (MD_IS_SCALAR(M))
  {
    vc    = VD_SCALCMP(v);
    mc    = MD_SCALCMP(M);
    dc    = VD_SCALCMP(d);
    mask  = VD_SCALTYPEMASK(v);

    /* solve LowerTriangle(v)=d */
    for (vec=first_vec; vec!= NULL; vec=SUCCVC(vec))
    {
      /* check coarse grid position */
      if (VTYPE(vec)==NODEVECTOR)
        if (CORNERTYPE(VMYNODE(vec))) {
          VVALUE(vec,vc) = 0.0;                                 /* no correction */
          continue;                                                             /* skip coarse grid node */
        }

      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)<myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        VVALUE(vec,vc) = (VVALUE(vec,dc)-sum);                           /* since Diag(L)=I per convention */
      }
    }

    /* solve UpperTriangle(v)=d */
    for (vec=last_vec; vec!= NULL; vec=PREDVC(vec))
    {
      /* check coarse grid position */
      if (VTYPE(vec)==NODEVECTOR)
        if (CORNERTYPE(VMYNODE(vec))) {
          VVALUE(vec,vc) = 0.0;                                 /* no correction */
          continue;                                                             /* skip coarse grid node */
        }

      myindex = VINDEX(vec);
      if ( (VDATATYPE(vec)&mask) && (VCLASS(vec)>=ACTIVE_CLASS) )
      {
        sum = 0.0;
        for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        {
          w = MDEST(mat);
          if ((VINDEX(w)>myindex) && (VDATATYPE(w)&mask) && (VCLASS(w)>=ACTIVE_CLASS) )
            sum += MVALUE(mat,mc)*VVALUE(w,vc);
        }
        /* solve (we stored the inverse of the diagonal) */
        if (StoreInverse)
          VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) * MVALUE(VSTART(vec),mc);
        else
          VVALUE(vec,vc) = (VVALUE(vec,vc)-sum) / MVALUE(VSTART(vec),mc);
      }
    }

    return (NUM_OK);
  }

  /* solve lower traingle */
  L_VLOOP__CLASS(vec,first_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n           = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    dcomp   = VD_CMPPTR_OF_TYPE(d,rtype);
    myindex = VINDEX(vec);

    /* check coarse grid position */
    if (rtype==NODEVECTOR)
      if (CORNERTYPE(VMYNODE(vec))) {
        vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
        vmat  = VVALPTR(vec);
        for (i=0; i<n; i++) vmat[vcomp[i]] = 0.0;
        continue;                         /* skip coarse grid node */
      }

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,dcomp[i]);
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
            {
              wmat  = VVALPTR(w);
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j]) * wmat[wcomp[j]];
            }
        }

    /* solve (Diag(L)=I per convention) */
    vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
    vmat  = VVALPTR(vec);
    for (i=0; i<n; i++)
      vmat[vcomp[i]] = s[i];
  }

  last_vec = LASTVECTOR(g);

  /* solve upper triangle */
  L_REVERSE_VLOOP__CLASS(vec,last_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    n = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    myindex = VINDEX(vec);

    /* check coarse grid position */
    if (rtype==NODEVECTOR)
      if (CORNERTYPE(VMYNODE(vec))) {
        vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
        vmat  = VVALPTR(vec);
        for (i=0; i<n; i++) vmat[vcomp[i]] = 0.0;
        continue;                         /* skip coarse grid node */
      }

    /* rhs */
    for (i=0; i<n; i++) s[i] = VVALUE(vec,VD_CMP_OF_TYPE(v,rtype,i));
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
        switch (MAT_RCKIND(M,rtype,ctype))
        {
        case R1C1 :
          SET_CMPS_11(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_11(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C2 :
          SET_CMPS_12(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_12(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R1C3 :
          SET_CMPS_13(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_13(s,mat,m,w,cy);
          s[0] -= s0;
          break;

        case R2C1 :
          SET_CMPS_21(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_21(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C2 :
          SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_22(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R2C3 :
          SET_CMPS_23(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_23(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          break;

        case R3C1 :
          SET_CMPS_31(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_31(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C2 :
          SET_CMPS_32(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_32(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        case R3C3 :
          SET_CMPS_33(cy,v,m,M,rtype,ctype,tmpptr);
          s0 = s1 = s2 = 0.0;
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
              MATMUL_33(s,mat,m,w,cy);
          s[0] -= s0;
          s[1] -= s1;
          s[2] -= s2;
          break;

        default :
          mcomp = MD_MCMPPTR_OF_RT_CT(M,rtype,ctype);
          wcomp = VD_CMPPTR_OF_TYPE(v,ctype);
          nc    = MD_COLS_IN_RT_CT(M,rtype,ctype);
          for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex<VINDEX(w)))
            {
              wmat  = VVALPTR(w);
              for (i=0; i<n; i++)
                for (j=0; j<nc; j++)
                  s[i] -= MVALUE(mat,mcomp[i*nc+j])
                          * wmat[wcomp[j]];
            }
        }


    /* solve */
    if (StoreInverse)
    {
      if (SolveInverseSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                                 VVALPTR(vec),
                                 MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                                 MVALPTR(VSTART(vec)),s)!=0)
        REP_ERR_RETURN (__LINE__);
    }
    else
    {
      if (SolveSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),
                          VVALPTR(vec),
                          MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                          MVALPTR(VSTART(vec)),s)!=0)
        REP_ERR_RETURN (__LINE__);
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   l_pgs - patch Gauss Seidel step

   SYNOPSIS:
   INT l_pgs (GRID *g, const VECDATA_DESC *v,
   const MATDATA_DESC *M, const VECDATA_DESC *d,
   INT depth, INT mode);

   PARAMETERS:
   .  g - pointer to grid
   .  v - type vector descriptor to store correction
   .  M - type matrix descriptor for precondition
   .  d - type vector descriptor for right hand side (the defect)
   .  depth - number of patch vectors
   .  mode - overlapping mode

   DESCRIPTION:
   This function performs a patch Gauss Seidel step.

   `Remark.` Index field must be set!

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    __LINE__ line where an error occured.
   D*/
/****************************************************************************/

INT l_pgs (GRID *g, const VECDATA_DESC *v,
           const MATDATA_DESC *M, const VECDATA_DESC *d,
           INT depth, INT mode)
{
  VECTOR *vec,*vlist[MAX_DEPTH],*w;
  MATRIX *mat;
  DOUBLE Mval[LOCAL_DIM*LOCAL_DIM],vval[LOCAL_DIM],dval[LOCAL_DIM];
  INT cnt,m,i,j,k,l,ncomp,vcnt,vtype,wtype,wncomp;
  INT myindex,index;
  const SHORT *Comp,*VComp;

  if (depth > MAX_DEPTH) {
    UserWriteF("l_pgs: MAX_DEPTH too small\n");
    REP_ERR_RETURN (__LINE__);
  }
  for (vec=FIRSTVECTOR(g); vec!= NULL; vec=SUCCVC(vec))
    SETVCUSED(vec,0);
  for (vec=FIRSTVECTOR(g); vec!= NULL; vec=SUCCVC(vec)) {
    if (VCUSED(vec)) continue;
    if (START(vec) == NULL) continue;
    if (VCLASS(vec) < ACTIVE_CLASS) continue;
    cnt = 1;
    vlist[0] = vec;
    if (depth > 1)
      for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat)) {
        w = MDEST(mat);
        if (VCUSED(w)) continue;
        vlist[cnt++] = w;
        if (cnt >= depth) break;
      }
    m = GetVlistMValues(cnt,vlist,M,Mval);
    if (m != GetVlistVValues(cnt,vlist,d,dval)) {
      UserWriteF("l_pgs: wrong dimension %d in local system %d\n",
                 m,GetVlistVValues(cnt,vlist,d,dval));
      REP_ERR_RETURN (__LINE__);
    }
    vcnt = 0;
    for (i=0; i<cnt; i++) {
      vtype = VTYPE(vlist[i]);
      ncomp = VD_NCMPS_IN_TYPE(d,vtype);
      for (mat=MNEXT(VSTART(vlist[i])); mat!=NULL; mat=MNEXT(mat)) {
        w = MDEST(mat);
        if (!VCUSED(w)) continue;
        if (VCLASS(w) < ACTIVE_CLASS) continue;
        wtype = VTYPE(w);
        Comp = MD_MCMPPTR_OF_MTYPE(M,MTP(vtype,wtype));
        wncomp = VD_NCMPS_IN_TYPE(d,wtype);
        VComp = VD_CMPPTR_OF_TYPE(v,wtype);
        for (k=0; k<ncomp; k++)
          for (l=0; l<wncomp; l++)
            dval[vcnt+k] -= MVALUE(mat,Comp[k*wncomp+l])
                            * VVALUE(w,VComp[l]);
      }
      vcnt += ncomp;
    }
    if (SolveFullMatrix(m,vval,Mval,dval)) {
      UserWriteF("l_pgs: solving on local patch failed\n");
      REP_ERR_RETURN (__LINE__);
    }
    SetVlistVValues(cnt,vlist,v,vval);
    for (i=0; i<cnt; i++)
      SETVCUSED(vec,1);
  }

  return (NUM_OK);
}
