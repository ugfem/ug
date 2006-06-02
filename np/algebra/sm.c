// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      sm.c                                                          */
/*                                                                          */
/* Purpose:   sparse matrix handling routines                               */
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

#include "config.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

#ifdef _2
#define __UG__
#endif

#ifdef _3
#define __UG__
#endif

#ifndef __UG__
   #include "sm.h"
   #define ERR_FILE_ID 6
#define NS_PREFIX
#else /* __UG__ */

#include "sm.h"
#include "debug.h"

#define ERR_RETURN(x) REP_ERR_RETURN(x)

REP_ERR_FILE;

#include "udm.h"  /* for MAX_MAT_COMP */

/****************************************************************************/
/** \brief Computes the size of a sparse matrix array

   \param nr - number of rows
   \param nc - number of columns
   \param comps - pointer to integer array
   \param NPtr - here N is stored
   \param NredPtr - here Nred is stored

   Positive numbers in the comps-array are transfered,
   negative ones mean non-existing entries in the sparse matrix. Equal
   positive numbers mean identified fields.

   \return <ul>
   <li> 0 ok </li>
   <li> 1 offset too large (increase MAX_MAT_COMP and recompile) </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX ComputeSMSizeOfArray (SHORT nr, SHORT nc, const SHORT *comps,
                                              SHORT *NPtr, SHORT *NredPtr)
{
  SHORT off,N,Nred;
  SHORT flag[MAX_NDOF];
  int i,j;

  /* reset flag field */
  for (i=0; i<MAX_NDOF; i++)
    flag[i] = 0;

  N = Nred = 0;
  for (i=0; i<nr; i++)
  {
    for (j=0; j<nc; j++)
    {
      off = comps[i*nc+j];
      if (off>=0)
      {
        if (off>=MAX_NDOF) return(1);

        N++;
        if (flag[off]==0)
        {
          Nred++;
          flag[off] = 1;
        }
      }
    }
  }

  *NPtr = N;
  *NredPtr = Nred;

  return (0);
}

/****************************************************************************/
/** \brief Computes the array form of a sparse matrix

   \param sm - sparse matrix
   \param comps - pointer to integer array

   Computes the array form of a sparse matrix. Positive numbers in the comps-array
   are taken from the sparse matrix format. -1 means a non-existing entry in the
   sparse matrix.
   Later these values will be interpreted as offsets in some value field.

   \return <ul>
   <li> 0 ok </li>
   <li> -1 if size of sparse matrix too large </li>
   <li> -2 if sparse matrix is not consistent </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX SM2Array (const SPARSE_MATRIX *sm, SHORT *comps)
{
  int i,j,nr,nc,off,posc;

  nr = sm->nrows;
  nc = sm->ncols;
  if (nr*nc>MAX_MAT_COMP)
    ERR_RETURN (-1);

  posc = sm->row_start[0];
  for (i=0; i<nr; i++)
  {
    for (j=0; j<nc; j++)
    {
      if (posc>=sm->row_start[i+1])
        off = -1;
      else
      {
        if (sm->col_ind[posc]==j)
        {
          off = sm->offset[posc];
          posc++;
        }
        else
          off = -1;
      }
      *comps++ = off;
    }
    if (posc!=sm->row_start[i+1])
      ERR_RETURN (-2);
  }

  return (0);
}


/****************************************************************************/
/** \brief Computes the sparse matrix form of an array

   \param nr - number of rows
   \param nc - number of columns
   \param comps - pointer to integer array
   \param sm - sparse matrix

   Computes the sparse matrix form of an offset
   array. Positive numbers in the comps-array are transfered, negative
   ones mean non-existing entries in the sparse matrix. Equal positive
   numbers mean identified fields. It is assumed that there is enough
   space in the arrays of the sparse matrix.

   \return <ul>
   <li> 0 ok </li>
   <li> 1 offset too large (increase MAX_MAT_COMP and recompile) </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX Array2SM (SHORT nr, SHORT nc, const SHORT *comps, SPARSE_MATRIX *sm)
{
  INT error;
  SHORT off,posc,N,Nred;
  SHORT flag[MAX_NDOF];
  int i,j;

  if (error=ComputeSMSizeOfArray(nr, nc, comps, &N, &Nred))
    return (error);

  /* reset flag field */
  for (i=0; i<MAX_NDOF; i++)
    flag[i] = 0;

  sm->nrows = nr;
  sm->ncols = nc;
  sm->N = N;

  sm->row_start = &(sm->components[0]);
  sm->col_ind = sm->row_start+nr+1;
  sm->offset = sm->col_ind+N;

  sm->row_start[0] = posc = 0;
  for (i=0; i<nr; i++)
  {
    for (j=0; j<nc; j++)
    {
      off = comps[i*nc+j];
      if (off>=0)
      {
        if (off>=MAX_NDOF)
          return(1);

        sm->col_ind[posc] = j;
        sm->offset[posc++] = off;

        if (flag[off]==0)
          flag[off] = 1;
      }
    }
    sm->row_start[i+1] = posc;
  }

  return (0);
}

/****************************************************************************/
/** \brief Transforms a string to a SM array

   \param n - size (i.e. rows*cols)
   \param str - pointer to char array consisting of [*0a-z]
   \param comps - pointer to integer array

   Transforms a string to a SM array. * means a non-zero entry, 0
   means a zero entry. a-z are used for identification of positions
   (they get the same offset). It is assumed that there is enough
   space in the sparse matrix array.

   \return <ul>
   <li> 0 ok </li>
   <li> 1 wrong format </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX String2SMArray (SHORT n, char *str, SHORT *comps)
{
  SHORT off;
  int i;
  char c;
  SHORT pos[26];         /* offset belonging to character */

  /* reset identification field */
  for (i=0; i<26; i++)
    pos[i] = -1;

  off=0;
  for (i=0; i<n; i++)
  {
    /* skip blanks */
    while ((c=*str++) != '\0')
    {
      if (c!=' ' && c!='\t' && c!= '\n')
        break;
    }
    if (c=='\0') return(1);

    /* zero entry? */
    if (c=='0')
    {
      *comps++ = -1;
      continue;
    }

    /* standard nonzero entry? */
    if (c=='*')
    {
      *comps++ = off++;
      continue;
    }

    if ((c<'a') || (c>'z')) return(-1);

    /* identified nonzero entry */
    if (pos[c-'a']>=0)
      *comps++ = pos[c-'a'];
    else
    {
      *comps++ = off;
      pos[c-'a'] = off++;
    }
  }

  return(0);
}

#endif /* __UG__ */

/****************************************************************************/
/** \brief Computes the reduced size of a sparse matrix

   \param sm - sparse matrix

   Computes the reduced size of a sparse matrix, i.e. takes possible
   identification into account.
   Warning: This is an O(N^2)-algorithm, but N should be moderate
   in comparison with the number of MATRIX structures anyhow.

   \return <ul>
   <li> positive: reduced size </li>
   <li> negative: an error occured </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX SM_Compute_Reduced_Size (SPARSE_MATRIX *sm)
{
  register INT i, j, off;
  register INT ident_count;

  if (sm->N<0)
    ERR_RETURN (-1);              /* Error: sparse matrix has negative size */

  ident_count = 0;
  for (i=0; i<sm->N; i++)
  {
    off = sm->offset[i];
    for (j=i+1; j<sm->N; j++)
    {
      if (sm->offset[j] == off)
      {
        /* identification: increase the counter and break the loop. */
        ident_count++;
        break;
      }
    }
  }

  return(sm->N-ident_count);
}

/****************************************************************************/
/** \brief Computes the reduced offset field of a sparse matrix

   \param sm - sparse matrix

   Computes the reduced offset field of a sparse matrix, i.e. counting
   identified only once. The field to take the reduced offsets should
   be large enough (i.e. apply SM_Compute_Reduced_Size before)!

   Warning: This is an O(N^2)-algorithm, but N should be moderate
   in comparison with the number of MATRIX structures anyhow.

   \return <ul>
   <li> positive: number of reduced offsets </li>
   <li> negative: an error occured </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX SM_Compute_Reduced_Offsets (SPARSE_MATRIX *sm, SHORT *reduced_offsets)
{
  register INT i, j, k, off;

  if (sm->N < 0)
    ERR_RETURN (-1);              /* Error: sparse matrix has negative size */

  k = 0;
  for (i=0; i<sm->N; i++)
  {
    off = sm->offset[i];

    /* did it already occur? */
    for (j=0; j<i; j++)
      if (sm->offset[j] == off)
        break;
    if (j<i)
      break;                   /* yes: no need to take it again */

    /* no: take it into the reduced_offsets */
    reduced_offsets[k++] = off;
  }

  return (k);
}

/****************************************************************************/
/** \brief Compares two sparse matrices

   \param sm1 - sparse matrix
   \param sm2 - sparse matrix

   Compares two sparse matrices

   \return <ul>
   <li> 0 ok </li>
   <li> 1 nrows not equal </li>
   <li> 2 ncols not equal </li>
   <li> 3 total number of nonzeros not equal </li>
   <li> 5 rows not equally long </li>
   <li> 6 column indices do not agree </li>
   <li> 7 offset structure is not compatible </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX SM_Compare (SPARSE_MATRIX *sm1, SPARSE_MATRIX *sm2)
{
  register INT i, j, off1, off2;

  if (sm1->nrows!=sm2->nrows) return(1);
  if (sm1->ncols!=sm2->ncols) return(2);
  if (sm1->N!=sm2->N) return(3);

  for (i=0; i<=sm1->ncols; i++)
    if (sm1->row_start[i]!=sm2->row_start[i])
      return(5);

  for (i=0; i<sm1->N; i++)
    if (sm1->col_ind[i]!=sm2->col_ind[i])
      return(6);

  /* now check, if the identification coincides */
  /* Warning: This is an O(N^2)-algorithm, but N should be moderate */
  /* in comparison with the number of MATRIX structures anyhow */
  for (i=0; i<sm1->N; i++)
  {
    off1 = sm1->offset[i];
    off2 = sm2->offset[i];
    for (j=i+1; j<sm1->N; j++)
    {
      if (sm1->offset[j] == off1)
      {
        if  (sm2->offset[j] != off2)
          return (7);
      }
      else
      {
        if (sm2->offset[j] == off2)
          return (7);
      }
    }
  }

  return(0);
}

/****************************************************************************/
/** \brief Computes a  ptrdiff_t-form of the sparse matrix offsets

   \param N - size of offset-field
   \param offset - pointer to offset field
   \param Diff - pointer to ptrdiff_t field

   Writes out the offset array in a ptrdiff_t-form, which should
   allow faster access for the (sparse) blas routines. The Diff-field
   should have place for N entries! Also the start_offset is returned,
   even if it might be easily obtained.

   \return <ul>
   <li> 0 ok </li>
   <li> -1 error </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX SM_Compute_Diff_From_Offset (INT N, SHORT *offset, ptrdiff_t *Diff)
{
  register int i;

  if (N<0)
    ERR_RETURN (-1);              /* Error: sparse matrix has negative size */

  if (N==0)
    return (0);              /* might be an error, but we prefer to do nothing */

  for (i=0; i<N; i++)
    *Diff++ = (ptrdiff_t) ((offset[(i+1)%N]-offset[i])*sizeof(DOUBLE));


  return(0);
}

/****************************************************************************/
/** \brief Computes a ptrdiff_t-form for a vector.

   \param N - size of col_ind-field
   \param col_ind - pointer to col_ind field
   \param cmp_off - pointer to y''s offset field
   \param Diff - pointer to ptrdiff_t field
   \param sm - sparse matrix
   \param compoff - component offsets
   \param startOff - start offset
   \param Diff - beginning of ptrdiff_t array

   Writes out a ptrdiff_t-field of for fast access of the components
   of a vector inside a sparse matrix-vector multiplication.  The
   Diff-field should have place for N entries!  The cmp_off field
   should be a VECDATA_DESC-offset-field for some type.

   \return <ul>
   <li> 0 ok </li>
   <li> -1 an error ocurred </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX SM_Compute_yDiff_From_Offset (INT N, SHORT *col_ind, SHORT *cmp_off,
                                                      ptrdiff_t *Diff)
{
  register int i;

  if (N<0)
    ERR_RETURN (-1);              /* Error: sparse matrix has negative size */

  if (N==0)
    return (0);              /* probably an error, but we prefer to do nothing */

  for (i=0; i<N; i++)
    *Diff++ = (ptrdiff_t) ((cmp_off[col_ind[(i+1)%N]]-cmp_off[col_ind[i]])
                           *sizeof(DOUBLE));

  return(0);
}

NS_PREFIX INT NS_PREFIX Decompose_LR_pivot (int n, DOUBLE *mat, int *pivot)
{
  register DOUBLE dinv, piv, sum, factor;
  register int i, j, k, off_i, off_j;

  for (i=0; i<n; i++)
    pivot[i] = i;

  /* LR factorization */
  for (i=0; i<n; i++)
  {
    /* pivot search */
    k = i;
    piv = fabs(mat[pivot[i]*n+i]);
    for (j=i+1; j<n; j++)
    {
      sum = fabs(mat[pivot[j]*n+i]);
      if (sum > piv)
      {
        k = j;
        piv = sum;
      }
    }

    /* if necessary, change pivot array */
    if (k != i)
    {
      j = pivot[k];
      pivot[k] = pivot[i];
      pivot[i] = j;
    }

    off_i = pivot[i]*n;
    dinv = mat[off_i+i];
    if (fabs(dinv)<DBL_EPSILON)
      return(1);

    dinv = mat[off_i+i] = 1.0/dinv;
    for (j=i+1; j<n; j++)
    {
      off_j = pivot[j]*n;
      factor = (mat[off_j+i] *= dinv);
      for (k=i+1; k<n; k++)
        mat[off_j+k] -= mat[off_i+k] * factor;
    }
  }

  return(0);
}

NS_PREFIX INT NS_PREFIX Solve_LR (int n, const DOUBLE *LR, const int *pivot, DOUBLE *x, const DOUBLE *b)
{
  register int i, j, off_i;
  register DOUBLE sum;

  /* solve L x' = b (note that b/LR must be accessed via pivot) */
  for (i=0; i<n; i++)
  {
    sum = b[pivot[i]];
    for (j=0; j<i; j++)
    {
      off_i = pivot[i]*n;
      sum -= LR[off_i+j] * x[j];
    }
    x[i] = sum;
  }

  /* now solve R x = x' */
  for (i=n-1; i>=0; i--)
  {
    off_i = pivot[i]*n;
    sum = x[i];
    for (j=i+1; j<n; j++)
      sum -= LR[off_i+j] * x[j];
    x[i] = sum * LR[off_i+i];                   /* in diagonal R_ii^-1 is stored */
  }

  return (0);
}

/****************************************************************************/
/** \brief LR decomposes a sparse matrix

   \sa
   Decompose_LR_pivot()
   Solve_LR()

   \param sm - sparse matrix
   \param LR - DOUBLE field of size nr^2
   \param pivot - SHORT field of size nr

   Computes an LR decomposition with pivoting from the sparse matrix sm.
   sm must be quadratic, LR and pivot field must be long enough!

   \return <ul>
   <li> 0 ok </li>
   <li> -1 error </li>
   </ul>
 */
/****************************************************************************/

NS_PREFIX INT NS_PREFIX SM_Decompose_LR_pivot (const SPARSE_MATRIX *sm, DOUBLE *values,
                                               DOUBLE *LR, int *pivot)
{
  register int i,j,k,n;
  register DOUBLE *Row_Ptr;

  n = sm->nrows;
  if (n!=sm->ncols)
    ERR_RETURN(-1);              /* Error: sparse matrix not symmetric */

  /* clear LR */
  for (i=0; i<n*n; i++)
    LR[i] = 0.0;

  /* then copy sm into LR */
  for (i=0; i<n; i++)
  {
    Row_Ptr = LR+i*n;
    for (j=sm->row_start[i]; j<sm->row_start[i+1]; j++)
    {
      k = sm->col_ind[j];
      if (k>=n)
        ERR_RETURN(-1);                          /* Error: column index too large */
      Row_Ptr[k] = values[sm->offset[j]];
    }
  }

  /* then call the block decomposition routine */
  return (Decompose_LR_pivot (n, LR, pivot));

}
