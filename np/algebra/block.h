// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      block.h                                                       */
/*                                                                          */
/* Purpose:   block solver (header file)                                                    */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de				                        */
/*																			*/
/* History:   Nov 27 95                                                                                 */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __BLOCK__
#define __BLOCK__

#include "np.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#ifdef __TWODIM__
#define LOCAL_DIM 21
#endif

#ifdef __THREEDIM__
#define LOCAL_DIM 42
#endif
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

INT MatMulSmallBlock (SHORT nr, SHORT nc, SHORT n,
                      const SHORT *mcomp1, const DOUBLE *mat1,
                      const DOUBLE *mat2, DOUBLE *resmat);
INT InvertSmallBlock (SHORT n, const SHORT *mcomp,
                      const DOUBLE *mat, DOUBLE *invmat);
INT SolveInverseSmallBlock (SHORT n, const SHORT *scomp, DOUBLE *sol,
                            const SHORT *invcomp, const DOUBLE *inv,
                            const DOUBLE *rhs);
INT SolveSmallBlock (SHORT n, const SHORT *scomp, DOUBLE *sol,
                     const SHORT *mcomp, const DOUBLE *mat, const DOUBLE *rhs);

INT InvertFullMatrix (INT n, DOUBLE mat[LOCAL_DIM][LOCAL_DIM],
                      DOUBLE invmat[LOCAL_DIM][LOCAL_DIM]);
INT InvertSpdMatrix (INT n, DOUBLE mat[LOCAL_DIM][LOCAL_DIM],
                     DOUBLE invmat[LOCAL_DIM][LOCAL_DIM]);
INT SolveFullMatrix (INT n, DOUBLE *sol, DOUBLE *mat, DOUBLE *rhs);
INT InvertFullMatrix_piv (INT n, DOUBLE *mat, DOUBLE *inv);

#endif
