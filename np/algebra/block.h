// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      block.h                                                       */
/*                                                                          */
/* Purpose:   block solver (header file)                                    */
/*                                                                          */
/* Author:    Christian Wieners                                             */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/* email:     ug@ica3.uni-stuttgart.de                                      */
/*                                                                          */
/* History:   Nov 27 95                                                     */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __BLOCK__
#define __BLOCK__

#include "np.h"
#include "disctools.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define LOCAL_DIM MAX_NODAL_VALUES

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/** \brief Perform a matrix multiplication on the small blocks
 *
 * \param nr                    Number of rows
 * \param nc                    Number of columns
 * \param mcomp1                Components of matrix1
 * \param mat1                  Array of matrix1
 * \param mat2                  Array of matrix2
 * \param resmat                Store mat1*mat2 here
 *
 * \return 0: ok  1: error
 */
/****************************************************************************/
INT MatMulSmallBlock (SHORT nr, SHORT nc, SHORT n,
                      const SHORT *mcomp1, const DOUBLE *mat1,
                      const DOUBLE *mat2, DOUBLE *resmat);

/****************************************************************************/
/** \brief Invert a small system of equations
 *
 * \param n                     Size of the small system (n*n)
 * \param mcomp                 Components of matrix to invert
 * \param mat                   Array of this matrix
 * \param invmat                Store inverse matrix here
 *
 * \return 0: ok  1: error
 */
/****************************************************************************/
INT InvertSmallBlock (SHORT n, const SHORT *mcomp,
                      const DOUBLE *mat, DOUBLE *invmat);

/****************************************************************************/
/** \brief ???
 *
 * \param n                     Size of the small system (n*n)
 * \param scomp                 Components of the solution
 * \param sol                   Array of the solution
 * \param invcomp               Components of inverse matrix
 * \param inv                   Array of this matrix
 * \param rhs                   Find right hand side here
 *
 * \return 0: ok  1: error
 */
/****************************************************************************/
INT SolveInverseSmallBlock (SHORT n, const SHORT *scomp, DOUBLE *sol,
                            const SHORT *invcomp, const DOUBLE *inv,
                            const DOUBLE *rhs);

/****************************************************************************/
/** \brief Solve a small system of equations
 *
 * \param n                     Size of the small system (n*n)
 * \param scomp                 Components of the solution
 * \param sol                   Array of the solution
 * \param mcomp                 Components of matrix to invert
 * \param mat                   Array of this matrix
 * \param rhs                   Find right hand side here
 *
 * \return 0: ok  1: error
 */
/****************************************************************************/
INT SolveSmallBlock (SHORT n, const SHORT *scomp, DOUBLE *sol,
                     const SHORT *mcomp, const DOUBLE *mat, DOUBLE *rhs);


INT InvertFullMatrix (INT n, DOUBLE mat[LOCAL_DIM][LOCAL_DIM],
                      DOUBLE invmat[LOCAL_DIM][LOCAL_DIM]);
INT InvertSpdMatrix (INT n, DOUBLE mat[LOCAL_DIM][LOCAL_DIM],
                     DOUBLE invmat[LOCAL_DIM][LOCAL_DIM]);
INT Choleskydecomposition (INT n, DOUBLE *mat, DOUBLE *chol);
INT SolveFullMatrix (INT n, DOUBLE *sol, DOUBLE *mat, DOUBLE *rhs);
INT InvertFullMatrix_piv (INT n, DOUBLE *mat, DOUBLE *inv);
INT SolveFullMatrix2 (INT n, DOUBLE *sol, DOUBLE *mat, DOUBLE *rhs);
INT InvertFullMatrix_gen (INT n, DOUBLE *mat, DOUBLE *inv, DOUBLE *rhs,
                          INT *ipv);

#endif
