// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      blasm.h                                                       */
/*                                                                          */
/* Purpose:   interfaces to sparse matrix blas routines                     */
/*                                                                          */
/* Author:    Nicolas Neuss                                                 */
/*            email: Nicolas.Neuss@IWR.Uni-Heidelberg.De                    */
/*                                                                          */
/* History:   02.01.98 begin sparse matrix routines                         */
/*            20.01.98 end of implementation phase                          */
/*            28.01.98 scalar matrix operations work                        */
/*            xx.02.98 matrix operations work                               */
/*                                                                          */
/* Note:      The files sm.[ch], blasm.[ch] may be obtained also in a       */
/*            standalone form under the GNU General Public License.         */
/*            The use inside of UG under the actual UG license              */
/*            is allowed.                                                   */
/*                                  HD, 13.7.99,  Nicolas Neuss.            */
/*                                                                          */
/****************************************************************************/

#ifndef __BLASM__
#define __BLASM__

#ifdef _2
#define __UG__
#endif

#ifdef _3
#define __UG__
#endif

#ifdef __UG__

        #ifndef __COMPILER__
                #include "compiler.h"
        #endif

        #ifndef __UDM__
        #include "udm.h"
        #endif

        #ifndef __NP__
        #include "np.h"
        #endif


#else /* not __UG__ */

        #ifndef __COMPILER__
                #define SHORT short
                #define INT int
                #define DOUBLE double
                #define __COMPILER__
        #endif

        #ifndef __SMALG__
        #include "algebra.h"
        #endif

#endif /* not __UG__ */

/* the following macros are defined also in blasv.h */
#ifndef __BLAS__
  #define BLAS_OP_MASK    0x000f
  #define BLAS_OP_SHIFT   0
  #define BLAS_LOOP_MASK  0x00f0
  #define BLAS_LOOP_SHIFT 4
  #define BLAS_MODE_MASK  0x0f00
  #define BLAS_MODE_SHIFT 8

  #ifndef __UG__
    #define ALL_VECTORS     0x0000
    #define ON_SURFACE      0x0001
  #endif /* not __UG__ */

  #define __BLAS__
#endif

/* loop type */
#define BLAS_LOOP_M     0x0000
#define BLAS_LOOP_MN    0x0001
#define BLAS_LOOP_Mxy   0x0002

/* operations for BLAS_LOOP_M */
#define BLAS_M_CLEAR    0x0000
#define BLAS_M_SET      0x0001

/* operations for BLAS_LOOP_MN */
#define BLAS_M_COPY     0x0000
#define BLAS_M_ADD1     0x0001
#define BLAS_M_MINUS1   0x0002
#define BLAS_M_SCALMUL  0x0003

/* operations for BLAS_LOOP_Mxy */
#define BLAS_MV_MUL      0x0000
#define BLAS_MV_MULADD   0x0001
#define BLAS_MV_MULMINUS 0x0002
#define BLAS_MV_LGS      0x0003
#define BLAS_MV_BILFORM  0x0004

/* mode types */
#define BLAS_SURFACE     0x0001
#define BLAS_VACTIVE     0x0002
#define BLAS_MACTIVE     0x0004

/* further modes for matrix blas */
#define MBLAS_MTYPE_MASK  0xf000
#define MBLAS_MTYPE_SHIFT 12
#define MBLAS_NONE       0x0000
#define MBLAS_DIAG       0x0001
#define MBLAS_LOWER      0x0002
#define MBLAS_NOT_UPPER  0x0003
#define MBLAS_UPPER      0x0004
#define MBLAS_NOT_LOWER  0x0005
#define MBLAS_NOT_DIAG   0x0006
#define MBLAS_ALL        0x0007

#ifdef __UG__
INT Mark_and_Sort_Matrix (GRID *grid, int operation);
#endif

INT MG_Matrix_Loop(MULTIGRID *mg, INT fl, INT tl, INT mode,
                   const MATDATA_DESC *M, const MATDATA_DESC *N,
                   const VECDATA_DESC *x, const VECDATA_DESC *y,
                   int N_vals, const DOUBLE *value, DOUBLE *result);

#endif /* __BLASM__ */
