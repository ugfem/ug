// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  np.h															*/
/*																			*/
/* Purpose:   numerics subsystem header file								*/
/*																			*/
/* Author:	  Klaus Johannsen/Henrik Rentz-Reichert							*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   25.03.95 begin, ug version 3.0								*/
/*			  09.12.95 transition to new descriptor formats (HRR)			*/
/*			  December 2, 1996 redesign of numerics                                                 */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __NP__
#define __NP__

#include "compiler.h"
#include "gm.h"
#include "algebra.h"
#include "ugenv.h"
#include "udm.h"
#include "numproc.h"
#include "scan.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* return codes of the numerics routines									*/
#define NUM_OK                                  0       /* everything ok						*/
#define NUM_OUT_OF_MEM                  1       /* could not allocate mem (for connect.)*/
/*#define NUM_DECOMP_FAILED	   -n	   any neg value: diag block singular	*/
#define NUM_DESC_MISMATCH               3       /* descriptors passed are inconsistent	*/
#define NUM_BLOCK_TOO_LARGE             4       /* block too large
                                                                                        (increase MAX_SINGLE_VEC_COMP)	*/
#define NUM_FORMAT_MISMATCH             5       /* user data size exceeded				*/
#define NUM_SMALL_DIAG                  6       /* diag entry too small to invert		*/
#define NUM_NO_COARSER_GRID             7       /* restrict called on grid level 0		*/
#define NUM_TYPE_MISSING                8       /* indicates one float for VEC_SCALAR	*/
#define NUM_ERROR                               9       /* other error							*/

/* modes for l_iluspdecomp */
#define SP_LOCAL                                0       /* modify locally						*/
#define SP_GLOBAL                               1       /* modify globally						*/

#define OPTIONLEN                       32
#define OPTIONLENSTR            "31"
#define VALUELEN                        64
#define VALUELENSTR                     "63"

/* matrix consitency modes */
#define MAT_DIAG_CONS         0
#define MAT_CONS              1
#define MAT_MASTER_CONS       2

/* special REP_ERR_RETURN macro */
#define NP_RETURN(err,intvar)           {intvar = __LINE__; REP_ERR_RETURN(err);}

/****************************************************************************/
/*																			*/
/* macros concerned with solving											*/
/*																			*/
/****************************************************************************/

/* formats for display routines */
#define DISPLAY_WIDTH                                   50
#define DISPLAY_NP_LI_FORMAT_SSSSS              "%-2s %-15.12s %-15.12s %-15.12s %-15.12s\n"
#define DISPLAY_NP_LI_FORMAT_SSSSSI             "%-2s %-15.12s %-15.12s %-15.12s %-15.12s %-2d\n"
#define DISPLAY_NP_LI_FORMAT_SSSSSS             "%-2s %-15.12s %-15.12s %-15.12s %-15.12s %-15.12s\n"
#define DISPLAY_NP_FORMAT_S                     "%-16.13s = "
#define DISPLAY_NP_FORMAT_SS                    "%-16.13s = %-35.32s\n"
#define DISPLAY_NP_FORMAT_SF                    "%-16.13s = %-7.4g\n"
#define DISPLAY_NP_FORMAT_SFF                   "%-16.13s = %-7.4g  %-7.4g\n"
#define DISPLAY_NP_FORMAT_SFFF                  "%-16.13s = %-7.4g  %-7.4g  %-7.4g\n"
#define DISPLAY_NP_FORMAT_SI                    "%-16.13s = %-2d\n"
#define DISPLAY_NP_FORMAT_SII                   "%-16.13s = %-2d  %-2d\n"
#define DISPLAY_NP_FORMAT_SIII                  "%-16.13s = %-2d  %-2d  %-2d\n"
#define DISPLAY_NP_FORMAT_FF                    "%-7.4g  %-7.4g\n"

#define CLEAR_VECTOR_OF_MG(m)                                 \
  { INT level;                                       \
    for (level=0; level<=TOPLEVEL((m)); level++)     \
      ClearIVector (GRID_ON_LEVEL((m),level));}

#define SCALE_VECTOR_OF_MG(m,v)                               \
  { INT level;                                       \
    for (level=0; level<=TOPLEVEL((m)); level++)     \
      ScaleIVector (GRID_ON_LEVEL((m),level),(v));}

#define CLEAR_VECSKIP_OF_GRID(g)                                \
  { VECTOR *theVector;                                 \
    for (theVector=FIRSTVECTOR((g)); theVector!= NULL; \
         theVector=SUCCVC(theVector))                  \
      VECSKIP(theVector) = 0;}

/****************************************************************************/
/*																			*/
/* structures concerned with symbolic user data management					*/
/*																			*/
/****************************************************************************/

typedef INT (*SetFuncProcPtr)(const DOUBLE_VECTOR, INT, DOUBLE *);
typedef INT (*TransGridProcPtr)(GRID *, const VECDATA_DESC *, const VECDATA_DESC *, const DOUBLE *);
typedef INT (*InterpolateNewVectorsProcPtr)(GRID *, const VECDATA_DESC *);

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* blas routines and iterative methods										*/
/*																			*/
/****************************************************************************/

#ifdef ModelP
INT l_vector_consistent (GRID *g, const VECDATA_DESC *x);
INT a_vector_consistent (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x);
INT l_ghostvector_consistent (GRID *g, const VECDATA_DESC *x);
INT l_vector_collect (GRID *g, const VECDATA_DESC *x);
INT a_vector_collect (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x);
INT l_matrix_consistent (GRID *g, const MATDATA_DESC *M, INT mode);
INT l_ghostvector_collect (GRID *g, const VECDATA_DESC *x);
INT l_vector_meanvalue (GRID *g, const VECDATA_DESC *x);
#endif

/* blas level 1 (vector operations) */
INT l_dset                      (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT a_dset                      (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT s_dset                      (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     DOUBLE a);

INT l_dcopy             (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
INT a_dcopy             (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
INT s_dcopy             (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const VECDATA_DESC *y);

INT l_dsetnonskip       (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT a_dsetnonskip       (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT s_dsetnonskip       (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, DOUBLE a);

INT l_dsetskip          (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE a);

INT l_dsetrandom        (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE a);

INT l_dsetfunc          (GRID *g,                                               const VECDATA_DESC *x, INT xclass, SetFuncProcPtr SetFunc);

INT l_dscale            (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
INT a_dscale        (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
INT s_dscale        (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, DOUBLE *a);

INT l_daxpy             (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT a_daxpy             (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT s_daxpy             (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const DOUBLE *a, const VECDATA_DESC *y);

INT l_dxdy                      (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT a_dxdy                      (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT s_dxdy                      (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const DOUBLE *a, const VECDATA_DESC *y);

INT l_ddot                      (const GRID *g,                                           const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y, DOUBLE *sp);
INT a_ddot                      (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y, DOUBLE *sp);
INT s_ddot                      (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const VECDATA_DESC *y, DOUBLE *sp);
INT s_ddot_sv           (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const VECDATA_DESC *y, DOUBLE *weight, DOUBLE *sv);

INT l_mean                      (const GRID *g, const VECDATA_DESC *x, INT xclass, DOUBLE *sp);


INT l_eunorm            (const GRID *g,                                           const VECDATA_DESC *x, INT xclass, DOUBLE *eu);
INT a_eunorm            (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, DOUBLE *eu);
INT s_eunorm            (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                      DOUBLE *eu);

/* blas level 1 (BLOCKVECTOR operations) on one gridlevel */
INT dsetB                       (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT dsetG                       (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT dsetBS                      (                                  const BLOCKVECTOR *bv,                                                         INT xcomp,                                              DOUBLE a);
INT dsetGS                      (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  DOUBLE a);

INT dsetfuncB           (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, SetFuncProcPtr SetFunc);
INT dsetfuncG           (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, SetFuncProcPtr SetFunc);
INT dsetfuncBS          (                                  const BLOCKVECTOR *bv,                                                 INT xcomp,                                              SetFuncProcPtr SetFunc);
INT dsetfuncGS          (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  SetFuncProcPtr SetFunc);

INT dcopyB                      (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
INT dcopyG                      (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
INT dcopyBS             (                                  const BLOCKVECTOR *bv,                                                         INT xcomp,                                              INT ycomp);
INT dcopyGS             (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  INT ycomp);

INT dscaleB             (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
INT dscaleG             (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
INT dscaleBS            (                                  const BLOCKVECTOR *bv,                                                         INT xcomp,                                              DOUBLE a);
INT dscaleGS            (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  DOUBLE a);

INT daddBS                      (                                       const BLOCKVECTOR *bv,                                                    INT xcomp,                                                                               INT ycomp);
INT dsubBS                      (                                       const BLOCKVECTOR *bv,                                                    INT xcomp,                                                                               INT ycomp);
INT dminusaddBS         (                                       const BLOCKVECTOR *bv,                                                    INT xcomp,                                                                               INT ycomp);

INT daxpyB                      (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT daxpyG                      (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT daxpyBS             (                                  const BLOCKVECTOR *bv,                                                         INT xcomp,                                              DOUBLE a,                INT ycomp);
INT daxpyGS             (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  DOUBLE a,                INT ycomp);

INT dxdyB                       (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT dxdyG                       (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
INT dxdyBS                      (                                  const BLOCKVECTOR *bv,                                                         INT xcomp,                                              DOUBLE a,        INT ycomp);
INT dxdyGS                      (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  DOUBLE a,        INT ycomp);

INT ddotB                       (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y, DOUBLE *sp);
INT ddotG                       (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y, DOUBLE *sp);
INT ddotBS                      (                                  const BLOCKVECTOR *bv,                                                         INT xcomp,                                              INT ycomp,                      DOUBLE *sp);
INT ddotGS                      (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  INT ycomp,                      DOUBLE *sp);

INT eunormB             (                                  const BLOCKVECTOR *bv,                                                         const VECDATA_DESC *x, INT xclass, DOUBLE *eu);
INT eunormG             (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, DOUBLE *eu);
INT eunormBS            (                                  const BLOCKVECTOR *bv,                                                         INT xcomp,                                              DOUBLE *eu);
INT eunormGS            (const GRID *grid, const BV_DESC *bvd,            const BV_DESC_FORMAT *bvdf, INT xcomp,                                                  DOUBLE *eu);


/* blas level 1 (Simple BLOCKVECTOR operations) on one gridlevel */
INT l_dcopy_SB          (BLOCKVECTOR *bv,                               const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
INT l_dscale_SB         (BLOCKVECTOR *bv,                               const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
INT l_daxpy_SB          (BLOCKVECTOR *theBV,                    const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);


/* blas level 2 (matrix (vector) operations) */
INT l_dmatset           (GRID *g,                                               const MATDATA_DESC *M, DOUBLE a);
INT s_dmatset           (MULTIGRID *mg, INT fl, INT tl, const MATDATA_DESC *M, DOUBLE a);

INT l_dmatcopy          (GRID *g,                                               const MATDATA_DESC *M1, const MATDATA_DESC *M2);
INT s_dmatcopy          (MULTIGRID *mg, INT fl, INT tl, const MATDATA_DESC *M1, const MATDATA_DESC *M2);

INT l_dmattranspose (GRID *g,                                           const MATDATA_DESC *M1, const MATDATA_DESC *M2);

INT l_dmatadd           (GRID *g, const MATDATA_DESC *M1, const MATDATA_DESC *M2);

INT l_dmatmul           (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT l_dmatmul_minus     (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT s_dmatmul           (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT s_dmatmul_minus     (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT s_dmatmul_set       (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT s_dtpmatmul_set     (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);

INT l_dtpmatmul         (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);

/* blas level 2 (matrix (BLOCKVECTOR) operations) on one gridlevel */
INT dmatsetB            (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const MATDATA_DESC *M, DOUBLE a);
INT dmatsetG            (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const MATDATA_DESC *M, DOUBLE a);
INT dmatsetBS           (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT mcomp,                            DOUBLE a);
INT dmatsetGS           (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT mcomp,                                DOUBLE a);

INT dmatcopyB           (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const MATDATA_DESC *M1, const MATDATA_DESC *M2);
INT dmatcopyG           (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const MATDATA_DESC *M1, const MATDATA_DESC *M2);
INT dmatcopyBS          (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT m1comp,                    INT m2comp);
INT dmatcopyGS          (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT m1comp,                        INT m2comp);

INT dmatcopyTransBS (const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT dest_comp, INT source_comp);

INT dmatmulB            (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT dmatmulG            (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT dmatmulBS           (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT xcomp,                                                    INT mcomp,                              INT ycomp);
INT dmatmulGS           (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT xcomp,                                                        INT mcomp,                              INT ycomp);

INT dmatmul_minusB      (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT dmatmul_minusG      (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
INT dmatmul_minusBS     (                                  const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT xcomp,                                                    INT mcomp,                              INT ycomp);
INT dmatmul_minusGS     (const GRID *grid, const BV_DESC *bvd_row,        const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT xcomp,                                                        INT mcomp,                              INT ycomp);

DOUBLE CalculateDefectAndNormBS( const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT d_comp, INT f_comp, INT K_comp, INT u_comp );

INT d2matmulBS          (const BLOCKVECTOR *bv_row1, const BV_DESC *bvd_col1, const BV_DESC *bvd_col2, const BV_DESC_FORMAT *bvdf, INT M_res_comp, INT M1comp, INT M2comp, GRID *grid );
INT d2matmul_minusBS(const BLOCKVECTOR *bv_row1, const BV_DESC *bvd_col1, const BV_DESC *bvd_col2, const BV_DESC_FORMAT *bvdf, INT M_res_comp, INT M1comp, INT M2comp, GRID *grid );

INT d3matmulBS          (const BLOCKVECTOR *bv_row1, const BV_DESC *bvd_col1, const BV_DESC *bvd_col2, const BV_DESC *bvd_col3, const BV_DESC_FORMAT *bvdf, INT M_res_comp, INT M1comp, INT M2comp, INT M3comp, GRID *grid );

/* blas level 2 (matrix (Simple BLOCKVECTOR) operations) on one gridlevel */
INT l_dmatset_SB        (BLOCKVECTOR *dest,     BLOCKVECTOR *source,const MATDATA_DESC *M, DOUBLE a);
INT l_dmatmul_set_SB(BLOCKVECTOR *theBVX,                       const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, BLOCKVECTOR *theBVY, const VECDATA_DESC *y, INT yclass);
INT l_dtpmatmul_set_SB(BLOCKVECTOR *theBVX,             const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, BLOCKVECTOR *theBVY, const VECDATA_DESC *y, INT yclass);
INT l_dmatmul_SB        (BLOCKVECTOR *theBVX,                   const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, BLOCKVECTOR *theBVY, const VECDATA_DESC *y, INT yclass);
INT l_dtpmatmul_SB      (BLOCKVECTOR *theBVX,                   const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, BLOCKVECTOR *theBVY, const VECDATA_DESC *y, INT yclass);
INT l_dmatmul_minus_SB (BLOCKVECTOR *theBVX,            const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, BLOCKVECTOR *theBVY, const VECDATA_DESC *y, INT yclass);

/* iterative methods */
INT l_ordervtypes       (GRID *g, const SHORT TypeOrder[NVECTYPES]);
INT l_setindex          (GRID *g);

INT l_jac                       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT jacBS                       (const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT K_comp, INT u_comp, INT f_comp );

INT l_lgs                       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_ugs                       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);

INT l_lgsB                      (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT gs_solveBS          (const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, DOUBLE eps, INT max_it, INT K_comp, INT u_comp, INT f_comp, INT aux_comp, INT verbose, INT eps_relative );

INT l_lsor                      (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d, const DOUBLE *damp);

/* iterative methods for Simple BLOCKVECTOR */
INT l_lgs_SB            (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_tplgs_SB          (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_ugs_SB            (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);


INT l_ilubthdecomp      (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VEC_SCALAR threshold, const VECDATA_DESC *rest, const VEC_SCALAR oldrestthresh);
INT l_ilubthdecomp_fine (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VEC_SCALAR threshold, const VECDATA_DESC *rest, const VEC_SCALAR oldrestthresh);
INT l_icdecomp      (GRID *g, const MATDATA_DESC *M);
INT l_iluspdecomp       (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VECDATA_DESC *t, INT mode, const VEC_SCALAR oldrestthresh);
INT l_lrdecomp          (GRID *g, const MATDATA_DESC *M);
INT LUDecomposeDiagBS(const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT A_comp, GRID *grid );
INT l_lrregularize      (GRID *theGrid, const MATDATA_DESC *M);
INT l_lrdecompB         (GRID *g, const MATDATA_DESC *M);
INT l_luiter            (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT solveLUMatBS        (const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT dest_comp, INT LU_comp, INT source_comp );
INT l_luiter_fine       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_luiterB           (GRID *g, const BLOCKVECTOR *bv, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_lltiter           (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);

/* intergrid transfer */

INT StandardRestrict (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp);
INT StandardInterpolateCorrection (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp);
INT StandardInterpolateNewVectors (GRID *FineGrid, const VECDATA_DESC *Sol);
INT StandardProject (GRID *CoarseGrid, const VECDATA_DESC *to,
                     const VECDATA_DESC *from);

INT MatDepRestrict                              (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const MATDATA_DESC *Mat, const DOUBLE *damp);
INT MatDepInterpolateCorrection (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const MATDATA_DESC *Mat, const DOUBLE *damp);
INT MatDepGalerkin                              (GRID *FineGrid, const MATDATA_DESC *Mat, const VECDATA_DESC *temp1, const VECDATA_DESC *temp2);
INT InstallInterpolationMatrix  (GRID *FineGrid, const MATDATA_DESC *Mat);

INT ScaledMGRestrict                       (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp);
INT InstallScaledRestrictionMatrix (GRID *FineGrid, const MATDATA_DESC *Mat, DOUBLE cut);
INT DiagonalScaleSystem                    (GRID *FineGrid, const MATDATA_DESC *Mat, const MATDATA_DESC *ConsMat, const VECDATA_DESC *rhs);

/* miscellaneous */
INT l_matflset (GRID *g, INT f);

/****************************************************************************/
/*																			*/
/* symbols and numprocs														*/
/*																			*/
/****************************************************************************/

/* create */
NP_BASE *       GetNumProcFromName                      (char *name);

/* step */
NP_BASE    *GetFirstNumProcType                 (void);
NP_BASE    *GetNextNumProcType                  (NP_BASE *);
NP_BASE   *GetFirstNumProc                              (void);
NP_BASE   *GetNextNumProc                               (NP_BASE *);

/* miscellaneous */
INT             ExecuteNumProc                          (NP_BASE *theNumProc, MULTIGRID *theMG, INT argc, char **argv);
INT             DisplayNumProc                          (NP_BASE *theNumProc);
INT             ListNumProc                             (NP_BASE *currNumProc);
INT                     SetNumProc                                      (NP_BASE *, INT, char **);
INT             InitNum                                         (void);
INT             GetVectorCompNames                      (VECDATA_DESC *theVDT, char *compNames, INT *nComp);
INT             WriteVEC_SCALAR                         (VECDATA_DESC *theVDT, VEC_SCALAR Scalar, char *structdir);

#ifdef __INTERPOLATION_MATRIX__
/* interpolation matrix functions */
INT AddInterpolationMatrix (GRID *theGrid,
                            ELEMENT *theElement, ELEMENT *theFather,
                            INT me, DOUBLE *IntMat, VECDATA_DESC *theVD);
INT ScaleIMatrix (GRID *g, VECDATA_DESC *theVD);
INT ClearIMatrix (GRID *g, VECDATA_DESC *theVD);
INT InterpolateCorrectionByMatrix (GRID *FineGrid, const VECDATA_DESC *to,
                                   const VECDATA_DESC *from,
                                   const DOUBLE *damp);
INT InterpolateCorrectionByMatrix_NoSkip (GRID *FineGrid, const VECDATA_DESC *to,
                                          const VECDATA_DESC *from,
                                          const DOUBLE *damp);
INT RestrictByMatrix              (GRID *FineGrid, const VECDATA_DESC *to,
                                   const VECDATA_DESC *from,
                                   const DOUBLE *damp);
INT RestrictByMatrix_s                    (GRID *FineGrid, const VECDATA_DESC *to,
                                           const VECDATA_DESC *from, const DOUBLE *damp);
INT InterpolateNewVectorsByMatrix (GRID *FineGrid, const VECDATA_DESC *sol);
INT AssembleGalerkinByMatrix (GRID *FineGrid, MATDATA_DESC *Mat);
#endif

INT ScaleIVector (GRID *g, VECDATA_DESC *theVD);
INT ClearIVector (GRID *g);

#endif
