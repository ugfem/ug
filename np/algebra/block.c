// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      block.c                                                       */
/*                                                                          */
/* Purpose:   block solver                                                                      */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
/*																			*/
/* History:   Nov 27 95 begin                                                                           */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <math.h>

#include "devices.h"
#include "compiler.h"
#include "gm.h"
#include "np.h"
#include "debug.h"
#include "general.h"

#include "block.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define SMALL_DET 1e-25

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* large matrices for InvertFullMatrix, InvertSpdMatrix and SolveFullMatrix2 */
/* (Macintosh cannot handle local data>32k) */
static DOUBLE BL_full_lrmat[LOCAL_DIM][LOCAL_DIM];
static DOUBLE BL_chol[LOCAL_DIM][LOCAL_DIM];
static DOUBLE BL_imat[LOCAL_DIM*LOCAL_DIM];
static DOUBLE BL_mat1[LOCAL_DIM*LOCAL_DIM];

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* Function:  SolveSmallBlock												*/
/*																			*/
/* Purpose:   solve a small system of equations								*/
/*																			*/
/* Input:	  SHORT   n				size of the small system (n*n)			*/
/*			  SHORT  *scomp			components of the solution				*/
/*			  DOUBLE *sol			DOUBLE array of the solution			*/
/*			  SHORT  *mcomp			components of matrix to invert			*/
/*			  DOUBLE *mat			DOUBLE array of this matrix				*/
/*			  DOUBLE *rhs			find right hand side here				*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT SolveSmallBlock (SHORT n, const SHORT *scomp, DOUBLE *sol,
                     const SHORT *mcomp, const DOUBLE *mat, const DOUBLE *rhs)
{
  DOUBLE BlockMat[MAX_SINGLE_MAT_COMP],BlockSol[MAX_SINGLE_VEC_COMP], det;
  DOUBLE aux,M3div0,M6div0;
  register DOUBLE dinv,piv,sum;
  register i,j,k;

  if (n>=MAX_SINGLE_VEC_COMP)
    return (1);

  switch (n)
  {
  case 1 :
    sol[scomp[0]] = rhs[0] / mat[mcomp[0]];
    return (NUM_OK);

  case 2 :
    det = mat[mcomp[0]]*mat[mcomp[3]] - mat[mcomp[1]]*mat[mcomp[2]];
    if (det==0.0) return (1);
    det = 1.0/det;
    sol[scomp[0]] = (rhs[0]*mat[mcomp[3]]-rhs[1]*mat[mcomp[1]])*det;
    sol[scomp[1]] = (rhs[1]*mat[mcomp[0]]-rhs[0]*mat[mcomp[2]])*det;
    return (NUM_OK);

  case 3 :
    M3div0 = mat[mcomp[3]]/mat[mcomp[0]];
    M6div0 = mat[mcomp[6]]/mat[mcomp[0]];
    aux = (mat[mcomp[7]]-M6div0*mat[mcomp[1]]) / (mat[mcomp[4]]-M3div0*mat[mcomp[1]]);

    sol[scomp[2]] = (rhs[2] - M6div0*rhs[0] - aux*(rhs[1]-M3div0*rhs[0]))
                    / (mat[mcomp[8]]-M6div0*mat[mcomp[2]]-aux*(mat[mcomp[5]]-M3div0*mat[mcomp[2]]));
    sol[scomp[1]] = (rhs[1] - mat[mcomp[3]]/mat[mcomp[0]]*rhs[0] - (mat[mcomp[5]]-M3div0*mat[mcomp[2]])*sol[scomp[2]])
                    / (mat[mcomp[4]]-M3div0*mat[mcomp[1]]);
    sol[scomp[0]] = (rhs[0] - mat[mcomp[1]] * sol[scomp[1]] - mat[mcomp[2]] * sol[scomp[2]])
                    /  mat[mcomp[0]];
    return (NUM_OK);

  default :
    /* copy matrix */
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        BlockMat[i*n+j] = mat[mcomp[i*n+j]];

    /* lr factorize mat */
    for (i=0; i<n; i++)
    {
      dinv = BlockMat[i*n+i];
      if (ABS(dinv)<SMALL_D)
        return (NUM_SMALL_DIAG);
      dinv = BlockMat[i*n+i] = 1.0/dinv;

      for (j=i+1; j<n; j++)
      {
        piv = (BlockMat[j*n+i] *= dinv);
        for (k=i+1; k<n; k++)
          BlockMat[j*n+k] -= BlockMat[i*n+k] * piv;
      }
    }

    /* solve */
    for (i=0; i<n; i++)
    {
      for (sum=rhs[i], j=0; j<i; j++)
        sum -= BlockMat[i*n+j] * BlockSol[j];
      BlockSol[i] = sum;                                /* Lii = 1 */
    }
    for (i=n-1; i>=0; i--)
    {
      for (sum=BlockSol[i], j=i+1; j<n; j++)
        sum -= BlockMat[i*n+j] * BlockSol[j];
      BlockSol[i] = sum * BlockMat[i*n+i];                              /* Uii = Inv(Mii) */
    }

    /* copy BlockSol to sol */
    for (i=0; i<n; i++)
      sol[scomp[i]] = BlockSol[i];

    return (NUM_OK);
  }
}

/****************************************************************************/
/*																			*/
/* Function:  SolveInverseSmallBlock										*/
/*																			*/
/* Purpose:   solve a small system of equations								*/
/*																			*/
/* Input:	  SHORT   n				size of the small system (n*n)			*/
/*			  SHORT  *scomp			components of the solution				*/
/*			  DOUBLE *sol			DOUBLE array of the solution			*/
/*			  SHORT  *invcomp		components of inverse matrix			*/
/*			  DOUBLE *inv			DOUBLE array of this matrix				*/
/*			  DOUBLE *rhs			find right hand side here				*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT SolveInverseSmallBlock (SHORT n, const SHORT *scomp, DOUBLE *sol,
                            const SHORT *invcomp, const DOUBLE *inv,
                            const DOUBLE *rhs)
{
  register DOUBLE sum;
  register i,j;

  if (n>=MAX_SINGLE_VEC_COMP)
    return (1);

  switch (n)
  {
  case 1 :
    sol[scomp[0]] = inv[invcomp[0]] * rhs[0];
    return (NUM_OK);

  default :

    /* sol = matrix * rhs */
    for (i=0; i<n; i++)
    {
      sum = 0.0;
      for (j=0; j<n; j++)
        sum += inv[invcomp[i*n+j]] * rhs[j];
      sol[scomp[i]] = sum;
    }

    return (NUM_OK);
  }
}

/****************************************************************************/
/*																			*/
/* Function:  InvertSmallBlock												*/
/*																			*/
/* Purpose:   solve a small system of equations								*/
/*																			*/
/* Input:	  SHORT   n				size of the small system (n*n)			*/
/*			  SHORT  *mcomp			components of matrix to invert			*/
/*			  DOUBLE *mat			DOUBLE array of this matrix				*/
/*			  DOUBLE *invmat		store inverse matrix here				*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT InvertSmallBlock (SHORT n, const SHORT *mcomp,
                      const DOUBLE *mat, DOUBLE *invmat)
{
  DOUBLE det,invdet,lrmat[MAX_SINGLE_MAT_COMP],sum,piv;
  INT i,j,k;

  switch (n)
  {
  case 1 :
    if (fabs(mat[mcomp[0]])<SMALL_DET)
    {
      UserWriteF("n=%d, c0=%d, m[c0]=%lf\n",(int)n,(int)(mcomp[0]),mat[mcomp[0]]);
      break;                            /* singular */
    }
    invmat[0] = 1.0 / mat[mcomp[0]];
    return (NUM_OK);

  case 2 :
    det = mat[mcomp[0]]*mat[mcomp[3]]-mat[mcomp[1]]*mat[mcomp[2]];
    if (ABS(det)<SMALL_DET)
      break;                            /* singular */
    invdet = 1.0/det;
    invmat[0] =  mat[mcomp[3]]*invdet;
    invmat[1] = -mat[mcomp[1]]*invdet;
    invmat[2] = -mat[mcomp[2]]*invdet;
    invmat[3] =  mat[mcomp[0]]*invdet;
    return (NUM_OK);

  case 3 :
    det = mat[mcomp[0]]*mat[mcomp[4]]*mat[mcomp[8]]
          + mat[mcomp[1]]*mat[mcomp[5]]*mat[mcomp[6]]
          + mat[mcomp[2]]*mat[mcomp[3]]*mat[mcomp[7]]
          - mat[mcomp[2]]*mat[mcomp[4]]*mat[mcomp[6]]
          - mat[mcomp[0]]*mat[mcomp[5]]*mat[mcomp[7]]
          - mat[mcomp[1]]*mat[mcomp[3]]*mat[mcomp[8]];
    if (ABS(det)<SMALL_DET)
      break;                            /* singular */
    invdet = 1.0/det;
    invmat[0] = ( mat[mcomp[4]]*mat[mcomp[8]] - mat[mcomp[5]]*mat[mcomp[7]]) * invdet;
    invmat[3] = (-mat[mcomp[3]]*mat[mcomp[8]] + mat[mcomp[5]]*mat[mcomp[6]]) * invdet;
    invmat[6] = ( mat[mcomp[3]]*mat[mcomp[7]] - mat[mcomp[4]]*mat[mcomp[6]]) * invdet;
    invmat[1] = (-mat[mcomp[1]]*mat[mcomp[8]] + mat[mcomp[2]]*mat[mcomp[7]]) * invdet;
    invmat[4] = ( mat[mcomp[0]]*mat[mcomp[8]] - mat[mcomp[2]]*mat[mcomp[6]]) * invdet;
    invmat[7] = (-mat[mcomp[0]]*mat[mcomp[7]] + mat[mcomp[1]]*mat[mcomp[6]]) * invdet;
    invmat[2] = ( mat[mcomp[1]]*mat[mcomp[5]] - mat[mcomp[2]]*mat[mcomp[4]]) * invdet;
    invmat[5] = (-mat[mcomp[0]]*mat[mcomp[5]] + mat[mcomp[2]]*mat[mcomp[3]]) * invdet;
    invmat[8] = ( mat[mcomp[0]]*mat[mcomp[4]] - mat[mcomp[1]]*mat[mcomp[3]]) * invdet;
    return (NUM_OK);

  default :
    if (n*n > MAX_SINGLE_MAT_COMP)
    {
      PrintErrorMessage('E',"InvertSmallMatrix","n too large");
      return (NUM_ERROR);
    }

    /* copy matrix */
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        lrmat[i*n+j] = mat[mcomp[i*n+j]];

    /* lr factorize mat */
    for (i=0; i<n; i++)
    {
      invdet = lrmat[i*n+i];
      if (ABS(invdet)<SMALL_DET)
        break;                          /* singular */
      invdet = lrmat[i*n+i] = 1.0/invdet;

      for (j=i+1; j<n; j++)
      {
        piv = (lrmat[j*n+i] *= invdet);
        for (k=i+1; k<n; k++)
          lrmat[j*n+k] -= lrmat[i*n+k] * piv;
      }
    }

    /* solve */
    for (k=0; k<n; k++)
    {
      for (i=0; i<k; i++)
        invmat[i*n+k] = 0.0;
      sum = 1.0;
      for (j=0; j<k; j++)
        sum -= lrmat[k*n+j] * invmat[j*n+k];
      invmat[k*n+k] = sum;                      /* Lii = 1 */
      for (i=k+1; i<n; i++)
      {
        sum = 0.0;
        for (j=0; j<i; j++)
          sum -= lrmat[i*n+j] * invmat[j*n+k];
        invmat[i*n+k] = sum;                            /* Lii = 1 */
      }
      for (i=n-1; i>=0; i--)
      {
        for (sum=invmat[i*n+k], j=i+1; j<n; j++)
          sum -= lrmat[i*n+j] * invmat[j*n+k];
        invmat[i*n+k] = sum * lrmat[i*n+i];                             /* Uii = Inv(Mii) */
      }
    }

    return (NUM_OK);
  }

  /*PrintErrorMessage('E',"InvertSmallBlock","singular block");*/
  return (1);
}

/****************************************************************************/
/*																			*/
/* Function:  MatMulSmallBlock												*/
/*																			*/
/* Purpose:   perform a matrix multiplication on the small blocks			*/
/*																			*/
/* Input:	  SHORT   nr			number of rows							*/
/*			  SHORT   nc			number of columns						*/
/*			  SHORT  *mcomp1		components of matrix1					*/
/*			  DOUBLE *mat1			DOUBLE array of matrix1					*/
/*			  DOUBLE *mat2			DOUBLE array of matrix2					*/
/*			  DOUBLE *resmat		store mat1*mat2 here					*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT MatMulSmallBlock (SHORT nr, SHORT nc, SHORT n,
                      const SHORT *mcomp1, const DOUBLE *mat1,
                      const DOUBLE *mat2, DOUBLE *resmat)
{
  register INT i,j,k;
  register DOUBLE sum;

  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
    {
      sum = 0.0;
      for (k=0; k<n; k++)
        sum += mat1[mcomp1[i*n+k]] * mat2[k*nc+j];
      resmat[i*nc+j] = sum;
    }

  return (NUM_OK);
}

INT InvertFullMatrix (INT n, DOUBLE mat[LOCAL_DIM][LOCAL_DIM],
                      DOUBLE invmat[LOCAL_DIM][LOCAL_DIM])
{
  DOUBLE det,invdet,piv,sum;
  INT i,j,k;

  switch (n)
  {
  case 1 :
    if (ABS(mat[0][0])<SMALL_DET)
      break;                    /* singular */
    invmat[0][0] = 1.0 / mat[0][0];
    return (0);

  case 2 :
    det = mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1];
    if (ABS(det)<SMALL_DET)
      break;                    /* singular */
    invdet = 1.0/det;
    invmat[0][0] =  mat[1][1]*invdet;
    invmat[0][1] = -mat[0][1]*invdet;
    invmat[1][0] = -mat[1][0]*invdet;
    invmat[1][1] =  mat[0][0]*invdet;
    return (0);

  case 3 :
    det = mat[0][0]*mat[1][1]*mat[2][2]
          + mat[0][1]*mat[1][2]*mat[2][0]
          + mat[0][2]*mat[1][0]*mat[2][1]
          - mat[0][2]*mat[1][1]*mat[2][0]
          - mat[0][0]*mat[1][2]*mat[2][1]
          - mat[0][1]*mat[1][0]*mat[2][2];
    if (ABS(det)<SMALL_DET)
      break;                    /* singular */
    invdet = 1.0/det;
    invmat[0][0] = ( mat[1][1]*mat[2][2]
                     - mat[1][2]*mat[2][1]) * invdet;
    invmat[1][0] = (-mat[1][0]*mat[2][2]
                    + mat[1][2]*mat[2][0]) * invdet;
    invmat[2][0] = ( mat[1][0]*mat[2][1]
                     - mat[1][1]*mat[2][0]) * invdet;
    invmat[0][1] = (-mat[0][1]*mat[2][2]
                    + mat[0][2]*mat[2][1]) * invdet;
    invmat[1][1] = ( mat[0][0]*mat[2][2]
                     - mat[0][2]*mat[2][0]) * invdet;
    invmat[2][1] = (-mat[0][0]*mat[2][1]
                    + mat[0][1]*mat[2][0]) * invdet;
    invmat[0][2] = ( mat[0][1]*mat[1][2]
                     - mat[0][2]*mat[1][1]) * invdet;
    invmat[1][2] = (-mat[0][0]*mat[1][2]
                    + mat[0][2]*mat[1][0]) * invdet;
    invmat[2][2] = ( mat[0][0]*mat[1][1]
                     - mat[0][1]*mat[1][0]) * invdet;
    return (0);

  default :
    if (n > LOCAL_DIM)
    {
      PrintErrorMessage('E',"InvertFullMatrix","n too large");
      return (1);
    }

    /* copy matrix */
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        BL_full_lrmat[i][j] = mat[i][j];

    /* lr factorize mat */
    for (i=0; i<n; i++)
    {
      invdet = BL_full_lrmat[i][i];
      if (ABS(invdet)<SMALL_DET)
        break;                          /* singular */
      invdet = BL_full_lrmat[i][i] = 1.0/invdet;

      for (j=i+1; j<n; j++)
      {
        piv = (BL_full_lrmat[j][i] *= invdet);
        for (k=i+1; k<n; k++)
          BL_full_lrmat[j][k] -= BL_full_lrmat[i][k] * piv;
      }
    }

    /* solve */
    for (k=0; k<n; k++)
    {
      for (i=0; i<k; i++)
        invmat[i][k] = 0.0;
      sum = 1.0;
      for (j=0; j<k; j++)
        sum -= BL_full_lrmat[k][j] * invmat[j][k];
      invmat[k][k] = sum;                       /* Lii = 1 */
      for (i=k+1; i<n; i++)
      {
        sum = 0.0;
        for (j=0; j<i; j++)
          sum -= BL_full_lrmat[i][j] * invmat[j][k];
        invmat[i][k] = sum;                             /* Lii = 1 */
      }
      for (i=n-1; i>=0; i--)
      {
        for (sum=invmat[i][k], j=i+1; j<n; j++)
          sum -= BL_full_lrmat[i][j] * invmat[j][k];
        invmat[i][k] = sum * BL_full_lrmat[i][i];                               /* Uii = Inv(Mii) */
      }
    }

    return (0);
  }

  PrintErrorMessage('E',"InvertFullMatrix","singular block");

  return (1);
}

static INT CholeskyDecomposition (INT n,
                                  DOUBLE mat[LOCAL_DIM][LOCAL_DIM],
                                  DOUBLE chol[LOCAL_DIM][LOCAL_DIM])
{
  DOUBLE piv,sum;
  INT i,j,k;

  for (i=0; i<n; i++)
  {
    sum = mat[i][i];
    for (k=0; k<i; k++)
      sum -= chol[i][k] * chol[i][k];
    if (sum < 0.0)
    {
      PrintErrorMessage('E',"CholeskyDecomposition","not spd");
      return (1);
    }
    chol[i][i] = piv = 1.0 / sqrt(sum);
    for (j=i+1; j<n; j++)
    {
      sum = mat[i][j];
      for (k=0; k<i; k++)
        sum -= chol[j][k] * chol[i][k];
      chol[j][i] = piv * sum;
    }
  }

  return(0);
}

INT InvertSpdMatrix (INT n, DOUBLE mat[LOCAL_DIM][LOCAL_DIM],
                     DOUBLE invmat[LOCAL_DIM][LOCAL_DIM])
{
  DOUBLE sum;
  INT i,j,k;

  if (n<4)
    return(InvertFullMatrix(n,mat,invmat));

  if (n > LOCAL_DIM)
  {
    PrintErrorMessage('E',"InvertSpdMatrix","n too large");
    return (1);
  }

  if (CholeskyDecomposition(n,mat,BL_chol))
    return(1);

  /* solve */
  for (k=0; k<n; k++)
  {
    for (i=0; i<k; i++)
      invmat[i][k] = 0.0;
    sum = 1.0;
    for (j=0; j<k; j++)
      sum -= BL_chol[k][j] * invmat[j][k];
    invmat[k][k] = sum * BL_chol[k][k];
    for (i=k+1; i<n; i++)
    {
      sum = 0.0;
      for (j=0; j<i; j++)
        sum -= BL_chol[i][j] * invmat[j][k];
      invmat[i][k] = sum * BL_chol[i][i];
    }
    for (i=n-1; i>=0; i--)
    {
      sum=invmat[i][k];
      for (j=i+1; j<n; j++)
        sum -= BL_chol[j][i] * invmat[j][k];
      invmat[i][k] = sum * BL_chol[i][i];
    }
  }

  return (0);
}

INT SolveFullMatrix (INT n, DOUBLE *sol, DOUBLE *mat, DOUBLE *rhs)
{
  register DOUBLE dinv,piv,sum;
  register i,j,k;
  INT ipv[LOCAL_DIM];

  if (n > LOCAL_DIM)
    return (1);

  for (i=0; i<n; i++)
    ipv[i] = i;

  /* lr factorize mat */
  for (i=0; i<n; i++)
  {
    k = i;
    piv = ABS(mat[i*n+i]);
    for (j=i+1; j<n; j++)
    {
      sum = ABS(mat[j*n+i]);
      if (sum > piv)
      {
        k = j;
        piv = sum;
      }
    }
    if (k != i)
    {
      j = ipv[i];
      ipv[i] = ipv[k];
      ipv[k] = j;
      for (j=0; j<n; j++)
      {
        sum = mat[k*n+j];
        mat[k*n+j] = mat[i*n+j];
        mat[i*n+j] = sum;
      }
    }

    dinv = mat[i*n+i];
    if (ABS(dinv)<SMALL_DET)
      RETURN (NUM_SMALL_DIAG);
    dinv = mat[i*n+i] = 1.0/dinv;
    for (j=i+1; j<n; j++)
    {
      piv = (mat[j*n+i] *= dinv);
      for (k=i+1; k<n; k++)
        mat[j*n+k] -= mat[i*n+k] * piv;
    }
  }

  /* solve */
  for (i=0; i<n; i++)
  {
    for (sum=rhs[ipv[i]], j=0; j<i; j++)
      sum -= mat[i*n+j] * sol[j];
    sol[i] = sum;               /* Lii = 1 */
  }
  for (i=n-1; i>=0; i--)
  {
    for (sum=sol[i], j=i+1; j<n; j++)
      sum -= mat[i*n+j] * sol[j];
    sol[i] = sum * mat[i*n+i];                  /* Uii = Inv(Mii) */
  }

  return (NUM_OK);
}

INT InvertFullMatrix_piv (INT n, DOUBLE *mat, DOUBLE *inv)
{
  register DOUBLE dinv,piv,sum;
  DOUBLE rhs[LOCAL_DIM];
  register i,j,k;
  INT ipv[LOCAL_DIM];

  if (n > LOCAL_DIM)
  {
    PrintErrorMessage('E',"InvertFullMatrix_piv","n too large");
    return (1);
  }

  for (i=0; i<n; i++)
    ipv[i] = i;

  /* lr factorize mat */
  for (i=0; i<n; i++)
  {
    k = i;
    piv = ABS(mat[i*n+i]);
    for (j=i+1; j<n; j++)
    {
      sum = ABS(mat[j*n+i]);
      if (sum > piv)
      {
        k = j;
        piv = sum;
      }
    }
    if (k != i)
    {
      j = ipv[i];
      ipv[i] = ipv[k];
      ipv[k] = j;
      for (j=0; j<n; j++)
      {
        sum = mat[k*n+j];
        mat[k*n+j] = mat[i*n+j];
        mat[i*n+j] = sum;
      }
    }

    dinv = mat[i*n+i];
    if (ABS(dinv)<SMALL_DET)
      RETURN (NUM_SMALL_DIAG);
    dinv = mat[i*n+i] = 1.0/dinv;
    for (j=i+1; j<n; j++)
    {
      piv = (mat[j*n+i] *= dinv);
      for (k=i+1; k<n; k++)
        mat[j*n+k] -= mat[i*n+k] * piv;
    }
  }

  /* solve */
  for (k=0; k<n; k++)
  {
    for (i=0; i<n; i++)
      rhs[i] = 0;
    rhs[k] = 1.0;
    for (i=0; i<n; i++)
    {
      for (sum=rhs[ipv[i]], j=0; j<i; j++)
        sum -= mat[i*n+j] * inv[j*n+k];
      inv[i*n+k] = sum;                         /* Lii = 1 */
    }
    for (i=n-1; i>=0; i--)
    {
      for (sum=inv[i*n+k], j=i+1; j<n; j++)
        sum -= mat[i*n+j] * inv[j*n+k];
      inv[i*n+k] = sum * mat[i*n+i];                    /* Uii = Inv(Mii) */
    }
  }

  return (NUM_OK);
}

INT SolveFullMatrix2 (INT n, DOUBLE *sol, DOUBLE *mat, DOUBLE *rhs)
{
  DOUBLE sum;
  INT i,j;

  for (i=0; i<n*n; i++)
    BL_mat1[i] = mat[i];
  if (InvertFullMatrix_piv(n,mat,BL_imat))
    return(NUM_ERROR);
  for (i=0; i<n; i++) {
    sum = 0.0;
    for (j=0; j<n; j++)
      sum += BL_imat[i*n+j] * rhs[j];
    sol[i] = sum;
  }
  for (i=0; i<n; i++) {
    sum = rhs[i];
    for (j=0; j<n; j++)
      sum -= BL_mat1[i*n+j] * sol[j];
    rhs[i] = sum;
  }
  for (i=0; i<n; i++) {
    sum = 0.0;
    for (j=0; j<n; j++)
      sum += BL_imat[i*n+j] * rhs[j];
    sol[i] += sum;
  }
  return(NUM_OK);
}
