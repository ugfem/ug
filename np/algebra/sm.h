// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      sm.h                                                          */
/*                                                                          */
/* Purpose:   interfaces to sparse matrix handling routines                 */
/*                                                                          */
/* Author:    Nicolas Neuss                                                 */
/*            email: Nicolas.Neuss@IWR.Uni-Heidelberg.De                    */
/*                                                                          */
/* History:   01.98 begin sparse matrix routines                            */
/*                                                                          */
/* Note:      The files sm.[ch], blasm.[ch] may be obtained also in a       */
/*            standalone form under the GNU General Public License.         */
/*            The use inside of UG under the actual UG license              */
/*            is allowed.                                                   */
/*                                  HD, 13.7.99,  Nicolas Neuss.            */
/*                                                                          */
/****************************************************************************/

#ifndef __SM__
#define __SM__

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

    #include <stddef.h>

#else /* not __UG__ */

    #ifndef __GENERAL__
        #include "general.h"
    #endif

#endif /* not __UG__ */

/****************************************************************************/
/*																			*/
/* The sparse matrix structure                                                                          */
/*																			*/
/****************************************************************************/

typedef struct {

  SHORT nrows;                         /* number of rows                        */
  SHORT ncols;                         /* number of columns                     */
  SHORT N;                             /* total number of nonzero elements      */

  SHORT *row_start;                    /* pointer to nrows+1 row starts         */
  SHORT *col_ind;                      /* pointer to N column indices           */
  SHORT *offset;                       /* pointer to N value offsets            */

  SHORT components[1];

} SPARSE_MATRIX;
/* usually there will be nrows+2*N SHORTs allocated after this structure  */

/****************************************************************************/
/*																			*/
/* Routines for sparse matrix handling                                                                  */
/*																			*/
/****************************************************************************/

#ifdef __UG__
INT ComputeSMSizeOfArray (SHORT nr, SHORT nc, const SHORT *comps,
                          SHORT *NPtr, SHORT *NredPtr);
INT SM2Array             (const SPARSE_MATRIX *sm, SHORT *comps);
INT Array2SM             (SHORT nr, SHORT nc, const SHORT *comps,
                          SPARSE_MATRIX *sm);
INT String2SMArray       (SHORT n, char *str, SHORT *comps);
#endif

INT SM_Compute_Reduced_Size    (SPARSE_MATRIX *sm);
INT SM_Compute_Reduced_Offsets (SPARSE_MATRIX *sm, SHORT *reduced_offsets);
INT SM_Compare                 (SPARSE_MATRIX *sm1, SPARSE_MATRIX *sm2);

/* for the sparse BLAS routines */
INT SM_Compute_Diff_From_Offset  (INT N, SHORT *offset, ptrdiff_t *Diff);
INT SM_Compute_yDiff_From_Offset (INT N, SHORT *col_ind, SHORT *cmp_off,
                                  ptrdiff_t *Diff);

INT Decompose_LR_pivot    (int n, DOUBLE *mat, int *pivot);
INT Solve_LR              (int n, const DOUBLE *LR, const int *pivot,
                           DOUBLE *x, const DOUBLE *b);
INT SM_Decompose_LR_pivot (const SPARSE_MATRIX *sm, DOUBLE *values,
                           DOUBLE *LR, int *pivot);

#endif /* __SM__ */
