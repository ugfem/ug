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


/* RCS_ID
   $Header$
 */

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
#include "npscan.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* if FF_PARALLEL_SIMULATION is defined, special functions from fe/ff are linked */
/*#define FF_PARALLEL_SIMULATION*/
/*#define FF_ModelP*/		/* um temp. Verweise von np nach fe/ff fuer die Allgemeinheit auszublenden */

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
#define MAT_GHOST_DIAG_CONS   3
#define MAT_DIAG_VEC_CONS     4

/* special REP_ERR_RETURN macro */
#define NP_RETURN(err,intvar)           {intvar = __LINE__; REP_ERR_RETURN(err);}

/****************************************************************************/
/*																			*/
/* macros concerned with solving											*/
/*																			*/
/****************************************************************************/

/* formats for display routines */
#define DISPLAY_WIDTH                                   50
#define DISPLAY_NP_BAR                                  "--------------------------------------------------\n"
#define DISPLAY_NP_LI_FORMAT_SSSSS              "%-2s %-15.12s %-15.12s %-15.12s %-15.12s\n"
#define DISPLAY_NP_LI_FORMAT_SSSSSI             "%-2s %-15.12s %-15.12s %-15.12s %-15.12s %-2d\n"
#define DISPLAY_NP_LI_FORMAT_SSSSSS             "%-2s %-15.12s %-15.12s %-15.12s %-15.12s %-15.12s\n"
#define DISPLAY_NP_FORMAT_S                     "%-16.13s = "
#define DISPLAY_NP_FORMAT_SS                    "%-16.13s = %-35.32s\n"
#define DISPLAY_NP_FORMAT_SSS                   "%-16.13s = %-15.12s %-15.12s\n"
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
INT a_outervector_consistent (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x);
INT l_ghostvector_project (GRID *g, const VECDATA_DESC *x);
INT l_vector_collect (GRID *g, const VECDATA_DESC *x);
INT a_vector_collect (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x);
INT l_matrix_consistent (GRID *g, const MATDATA_DESC *M, INT mode);
INT l_ghostvector_collect (GRID *g, const VECDATA_DESC *x);
INT l_vector_meanvalue (GRID *g, const VECDATA_DESC *x);
INT a_vector_meanvalue (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x);
INT l_ghostmatrix_collect (GRID *g, const MATDATA_DESC *A);
INT a_vector_vecskip (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x);
INT l_amgmatrix_collect (GRID *g, const MATDATA_DESC *A);
int DDD_InfoPrioCopies (DDD_HDR hdr);
INT a_elementdata_consistent (MULTIGRID *mg, INT fl, INT tl);
INT l_vector_consistent_noskip (GRID *g, const VECDATA_DESC *x);
INT a_vector_consistent_noskip (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x);

#ifdef __BLOCK_VECTOR_DESC__
INT l_vector_consistentBS (GRID *g, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT x);
#endif

#endif


/* modus for blas routines                                                  */
#define ON_SURFACE      -1      /* class on surface                                     */
#define ALL_VECTORS      0      /* all vectors                                          */

/* blas level 1 (vector operations) */

INT dset           (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    DOUBLE a);
INT dcopy          (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VECDATA_DESC *y);
INT dscal          (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    DOUBLE a);
INT dscalx         (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VEC_SCALAR a);
INT dadd           (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VECDATA_DESC *y);
INT dsub           (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VECDATA_DESC *y);
INT dminusadd      (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VECDATA_DESC *y);
INT daxpy          (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    DOUBLE a, const VECDATA_DESC *y);
INT daxpyx         (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VEC_SCALAR a, const VECDATA_DESC *y);
INT ddot           (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VECDATA_DESC *y, DOUBLE *a);
INT ddotx          (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VECDATA_DESC *y, VEC_SCALAR a);
INT ddotw          (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const VECDATA_DESC *y, const VEC_SCALAR w, DOUBLE *a);
INT dnrm2          (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    DOUBLE *a);
INT dnrm2x         (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    VEC_SCALAR a);


/* blas level 2 (matrix operations) */

INT dmatclear      (MULTIGRID *mg, INT fl, INT tl, INT mode, const MATDATA_DESC *M);

INT dmatset        (MULTIGRID *mg, INT fl, INT tl, INT mode, const MATDATA_DESC *M,
                    DOUBLE a);
INT dmatcopy       (MULTIGRID *mg, INT fl, INT tl, INT mode,
                    const MATDATA_DESC *M, const MATDATA_DESC *N);
INT dmatadd        (MULTIGRID *mg, INT fl, INT tl, INT mode,
                    const MATDATA_DESC *M, const MATDATA_DESC *N);
INT dmatmul        (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const MATDATA_DESC *M, const VECDATA_DESC *y);
INT dmatmul_add    (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const MATDATA_DESC *M, const VECDATA_DESC *y);
INT dmatmul_minus  (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x,
                    const MATDATA_DESC *M, const VECDATA_DESC *y);

/* for compatibility only */

#define l_dset(g,x,xclass,a)               dset(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,a)
#define a_dset(mg,fl,tl,x,xclass,a)        dset(mg,fl,tl,ALL_VECTORS,x,a)
#define s_dset(mg,fl,tl,x,a)               dset(mg,fl,tl,ON_SURFACE,x,a)

#define l_dcopy(g,x,xclass,y)              dcopy(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,y)
#define a_dcopy(mg,fl,tl,x,xclass,y)       dcopy(mg,fl,tl,ALL_VECTORS,x,y)
#define s_dcopy(mg,fl,tl,x,y)              dcopy(mg,fl,tl,ON_SURFACE,x,y)

#define l_dscale(g,x,xclass,a)             dscalx(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,a)
#define a_dscale(mg,fl,tl,x,xclass,a)      dscalx(mg,fl,tl,ALL_VECTORS,x,a)
#define s_dscale(mg,fl,tl,x,a)             dscalx(mg,fl,tl,ON_SURFACE,x,a)

#define l_daxpy(g,x,xclass,a,y)            daxpyx(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,a,y)
#define a_daxpy(mg,fl,tl,x,xclass,a,y)     daxpyx(mg,fl,tl,ALL_VECTORS,x,a,y)
#define s_daxpy(mg,fl,tl,x,a,y)            daxpyx(mg,fl,tl,ON_SURFACE,x,a,y)

#define l_ddot(g,x,xclass,y,a)             ddotx(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,y,a)
#define a_ddot(mg,fl,tl,x,xclass,y,a)      ddotx(mg,fl,tl,ALL_VECTORS,x,y,a)
#define s_ddot(mg,fl,tl,x,y,a)             ddotx(mg,fl,tl,ON_SURFACE,x,y,a)

#define l_ddot_sv(g,x,xclass,y,b,a)        ddotw(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,y,b,a)
#define a_ddot_sv(mg,fl,tl,x,xclass,y,b,a) ddotw(mg,fl,tl,ALL_VECTORS,x,y,b,a)
#define s_ddot_sv(mg,fl,tl,x,y,b,a)        ddotw(mg,fl,tl,ON_SURFACE,x,y,b,a)

#define l_eunorm(g,x,xclass,a)             dnrm2x(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,a)
#define a_eunorm(mg,fl,tl,x,xclass,a)      dnrm2x(mg,fl,tl,ALL_VECTORS,x,a)
#define s_eunorm(mg,fl,tl,x,a)             dnrm2x(mg,fl,tl,ON_SURFACE,x,a)

#define l_dmatset(g,M,a)                   dmatset(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,M,a)
#define s_dmatset(mg,fl,tl,M,a)            dmatset(mg,fl,tl,ON_SURFACE,M,a)

#define l_dmatcopy(g,M,N)                  dmatcopy(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,M,N)
#define s_dmatcopy(mg,fl,tl,M,N)           dmatcopy(mg,fl,tl,ON_SURFACE,M,N)

#define l_dmatadd(g,M,N)                   dmatadd(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,M,N)
#define s_dmatadd(mg,fl,tl,M,N)            dmatadd(mg,fl,tl,ON_SURFACE,M,N)

#define l_dmatmul_set(g,x,xc,M,y,yc)       dmatmul(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,M,y)
#define s_dmatmul_set(mg,fl,tl,x,M,y,yc)   dmatmul(mg,fl,tl,ON_SURFACE,x,M,y)

#define l_dmatmul(g,x,xc,M,y,yc)           dmatmul_add(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,M,y)
#define s_dmatmul(mg,fl,tl,x,M,y,yc)       dmatmul_add(mg,fl,tl,ON_SURFACE,x,M,y)

#define l_dmatmul_minus(g,x,xc,M,y,yc)     dmatmul_minus(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,x,M,y)
#define s_dmatmul_minus(mg,fl,tl,x,M,y,yc) dmatmul_minus(mg,fl,tl,ON_SURFACE,x,M,y)

/* old style **********************

   INT l_dset			(GRID *g,						const VECDATA_DESC *x, INT xclass, DOUBLE a);
   INT a_dset			(MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, DOUBLE a);
   INT s_dset			(MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,			   DOUBLE a);

   INT l_dscale		(GRID *g,						const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
   INT a_dscale        (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
   INT s_dscale        (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, DOUBLE *a);

   INT l_ddot			(const GRID *g,						  const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y, DOUBLE *sp);
   INT a_ddot			(const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y, DOUBLE *sp);
   INT s_ddot			(const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,			 const VECDATA_DESC *y, DOUBLE *sp);

   INT l_dcopy          (GRID *g,						const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
   INT a_dcopy          (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
   INT s_dcopy          (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,			   const VECDATA_DESC *y);

   INT l_daxpy          (GRID *g,						const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
   INT a_daxpy          (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);
   INT s_daxpy          (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,			   const DOUBLE *a, const VECDATA_DESC *y);

   INT l_eunorm                 (const GRID *g,                                           const VECDATA_DESC *x, INT xclass, DOUBLE *eu);
   INT a_eunorm                 (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, DOUBLE *eu);
   INT s_eunorm                 (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                      DOUBLE *eu);

   INT l_ddot_sv                (const GRID *g,						  const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y, DOUBLE *weight, DOUBLE *sv);
   INT s_ddot_sv                (const MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const VECDATA_DESC *y, DOUBLE *weight, DOUBLE *sv);

   INT l_dmatset		(GRID *g,						const MATDATA_DESC *M, DOUBLE a);
   INT s_dmatset		(MULTIGRID *mg, INT fl, INT tl, const MATDATA_DESC *M, DOUBLE a);

   INT l_dmatcopy		(GRID *g,						const MATDATA_DESC *M1, const MATDATA_DESC *M2);
   INT s_dmatcopy		(MULTIGRID *mg, INT fl, INT tl, const MATDATA_DESC *M1, const MATDATA_DESC *M2);

   INT l_dmatadd		(GRID *g, const MATDATA_DESC *M1, const MATDATA_DESC *M2);

   INT l_dmatmul		(GRID *g,						const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
   INT s_dmatmul		(MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,			   const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);

   INT l_dmatmul_set	(GRID *g,						const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
   INT s_dmatmul_set	(MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,			   const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);

   INT l_dmatmul_minus	(GRID *g,						const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);
   INT s_dmatmul_minus	(MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,			   const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);

 **************************** old style */



INT l_dsetnonskip       (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT a_dsetnonskip       (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT s_dsetnonskip       (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, DOUBLE a);

INT l_dsetskip          (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE a);

INT l_dsetrandom        (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE a);
INT l_dsetrandom2       (GRID *g,                                               const VECDATA_DESC *x, INT xclass, DOUBLE from, DOUBLE to, INT skip);

INT l_dsetfunc          (GRID *g,                                               const VECDATA_DESC *x, INT xclass, SetFuncProcPtr SetFunc);

INT l_mean                      (const GRID *g, const VECDATA_DESC *x, INT xclass, DOUBLE *sp);

/* blas level 1 (BLOCKVECTOR operations) on one gridlevel */
INT dsetBS                      (const BLOCKVECTOR *bv, INT xc, DOUBLE a);
INT dcopyBS             (const BLOCKVECTOR *bv, INT xc, INT yc);
INT dscalBS             (const BLOCKVECTOR *bv, INT xc, DOUBLE a);
INT daddBS                      (const BLOCKVECTOR *bv, INT xc, INT yc);
INT dsubBS                      (const BLOCKVECTOR *bv, INT xc, INT yc);
INT dminusaddBS         (const BLOCKVECTOR *bv, INT xc, INT yc);
INT daxpyBS             (const BLOCKVECTOR *bv, INT xc, DOUBLE a, INT yc);
INT ddotBS                      (const BLOCKVECTOR *bv, INT xc, INT yc,   DOUBLE *a);
INT dnrm2BS             (const BLOCKVECTOR *bv, INT xc, DOUBLE *a);

/* blas level 2 (matrix (BLOCKVECTOR) operations) on one gridlevel */
INT dmatsetBS           (const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT mc, DOUBLE a);
INT dmatcopyBS          (const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT mc, INT nc);
INT dmataddBS           (const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT mc, INT nc);
INT dmatmulBS           (const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT xc, INT mc, INT yc);
INT dmatmul_addBS       (const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT xc, INT mc, INT yc);
INT dmatmul_minusBS     (const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT xc, INT mc, INT yc);

INT d2matmulBS          (const BLOCKVECTOR *bv_row1, const BV_DESC *bvd_col1, const BV_DESC *bvd_col2, const BV_DESC_FORMAT *bvdf, INT M_res_comp, INT M1comp, INT M2comp, GRID *grid );
INT d2matmul_minusBS(const BLOCKVECTOR *bv_row1, const BV_DESC *bvd_col1, const BV_DESC *bvd_col2, const BV_DESC_FORMAT *bvdf, INT M_res_comp, INT M1comp, INT M2comp, GRID *grid );
INT d3matmulBS          (const BLOCKVECTOR *bv_row1, const BV_DESC *bvd_col1, const BV_DESC *bvd_col2, const BV_DESC *bvd_col3, const BV_DESC_FORMAT *bvdf, INT M_res_comp, INT M1comp, INT M2comp, INT M3comp, GRID *grid );
DOUBLE CalculateDefectAndNormBS( const BLOCKVECTOR *bv_row, const BV_DESC *bvd_col, const BV_DESC_FORMAT *bvdf, INT d_comp, INT f_comp, INT K_comp, INT u_comp );

/* blas level 1 (Simple BLOCKVECTOR operations) on one gridlevel */
INT l_dcopy_SB          (BLOCKVECTOR *bv,                               const VECDATA_DESC *x, INT xclass, const VECDATA_DESC *y);
INT l_dscale_SB         (BLOCKVECTOR *bv,                               const VECDATA_DESC *x, INT xclass, const DOUBLE *a);
INT l_daxpy_SB          (BLOCKVECTOR *theBV,                    const VECDATA_DESC *x, INT xclass, const DOUBLE *a, const VECDATA_DESC *y);

/* blas level 2 (matrix (vector) operations) */
INT l_dmattranspose (GRID *g,                                           const MATDATA_DESC *M1, const MATDATA_DESC *M2);
INT s_dtpmatmul_set     (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x,                     const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);

INT l_dtpmatmul         (GRID *g,                                               const VECDATA_DESC *x, INT xclass, const MATDATA_DESC *M, const VECDATA_DESC *y, INT yclass);

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

INT l_lgs                       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d, VECDATA_DESC *diag);
INT l_ugs                       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);

INT l_lgsB                      (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT gs_solveBS          (const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, DOUBLE eps, INT max_it, INT K_comp, INT u_comp, INT f_comp, INT aux_comp, INT verbose, INT eps_relative );

INT l_lsor                      (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d, const DOUBLE *damp);
INT l_usor (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
            const VECDATA_DESC *d, const DOUBLE *omega);
INT l_lsor_ld       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d, const VECDATA_DESC *damp);
INT l_usor_ld (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M,
               const VECDATA_DESC *d, VECDATA_DESC *omega);


/* iterative methods for Simple BLOCKVECTOR */
INT l_lgs_SB            (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_tplgs_SB          (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_ugs_SB            (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_luiter_SB         (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_tpluiter_SB       (BLOCKVECTOR *theBV, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_ilubdecomp_SB     (BLOCKVECTOR *theBV, const MATDATA_DESC *M, const VEC_SCALAR beta);


INT l_ilubthdecomp      (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VEC_SCALAR threshold, const VECDATA_DESC *rest, const VEC_SCALAR oldrestthresh);
INT l_ilubthdecomp_fine (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VEC_SCALAR threshold, const VECDATA_DESC *rest, const VEC_SCALAR oldrestthresh);
INT l_icdecomp      (GRID *g, const MATDATA_DESC *M);
INT l_iluspdecomp       (GRID *g, const MATDATA_DESC *M, const VEC_SCALAR beta, const VECDATA_DESC *t, INT mode, const VEC_SCALAR oldrestthresh);
INT l_lrdecomp          (GRID *g, const MATDATA_DESC *M);
INT LUDecomposeDiagBS(const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT A_comp, GRID *grid );
INT l_lrregularize      (GRID *theGrid, const MATDATA_DESC *M, INT restore);
INT l_lrdecompB         (GRID *g, const MATDATA_DESC *M);
INT l_luiter            (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT solveLUMatBS        (const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT dest_comp, INT LU_comp, INT source_comp );
INT l_luiter_fine       (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_luiterB           (GRID *g, const BLOCKVECTOR *bv, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_lltiter           (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d);
INT l_pgs           (GRID *g, const VECDATA_DESC *v, const MATDATA_DESC *M, const VECDATA_DESC *d, INT depth, INT mode, DOUBLE vdamp);

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
INT CreateStandardNodeRestProl  (GRID *fineGrid, INT ncomp);

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
INT             WriteVEC_SCALAR                         (const VECDATA_DESC *theVDT, const VEC_SCALAR Scalar, const char *structdir);

#ifdef __INTERPOLATION_MATRIX__
/* interpolation matrix functions */
INT GetInterpolationMatrix (ELEMENT *theElement, ELEMENT *theFather,
                            INT me, DOUBLE *IntMat, VECDATA_DESC *theVD);
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
INT AssembleGalerkinByMatrix (GRID *FineGrid, MATDATA_DESC *Mat, INT symmetric);
#endif

INT ScaleIVector (GRID *g, VECDATA_DESC *theVD);
INT ClearIVector (GRID *g);

#endif
