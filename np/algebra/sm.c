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
/*                                                                              */
/****************************************************************************/

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
#else /* __UG__ */

#include "sm.h"
#include "debug.h"

#define ERR_RETURN(x) REP_ERR_RETURN(x)

REP_ERR_FILE;

#include "udm.h"  /* for MAX_MAT_COMP */

/****************************************************************************/
/*D
   ComputeSMSizeOfArray - Computes the size of a sparse matrix array

   SYNOPSIS:
   static INT ComputeSMSizeOfArray (SHORT nr, SHORT nc, const SHORT *comps,
                                    SHORT *NPtr, SHORT *NredPtr)

   PARAMETERS:
   .  nr - number of rows
   .  nc - number of columns
   .  comps - pointer to integer array
   .  NPtr - here N is stored
   .  NredPtr - here Nred is stored

   DESCRIPTION:
   Positive numbers in the comps-array are transfered,
   negative ones mean non-existing entries in the sparse matrix. Equal
   positive numbers mean identified fields.

   RETURN VALUE:
   INT
   .n 0 ok
   .n 1 offset too large (increase MAX_MAT_COMP and recompile)

   SEE ALSO:
   SM2Array
   D*/
/****************************************************************************/

INT ComputeSMSizeOfArray (SHORT nr, SHORT nc, const SHORT *comps,
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
/*D
   SM2Array - Computes the array form of a sparse matrix

   SYNOPSIS:
   INT SM2Array (const SPARSE_MATRIX *sm, SHORT *comps)

   PARAMETERS:
   .  sm - sparse matrix
   .  comps - pointer to integer array

   DESCRIPTION:
   Computes the array form of a sparse matrix. Positive numbers in the comps-array
   are taken from the sparse matrix format. -1 means a non-existing entry in the
   sparse matrix.
   Later these values will be interpreted as offsets in some value field.

   RETURN VALUE:
   INT
   .n 0 ok
   .n -1 if size of sparse matrix too large
   .n -2 if sparse matrix is not consistent

   SEE ALSO:
   Array2SM
   D*/
/****************************************************************************/

INT SM2Array (const SPARSE_MATRIX *sm, SHORT *comps)
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
/*D
   Array2SM - Computes the sparse matrix form of an array

   SYNOPSIS:
   INT Array2SM (SHORT nr, SHORT nc, const SHORT *comps, SPARSE_MATRIX *sm)

   PARAMETERS:
   .  nr - number of rows
   .  nc - number of columns
   .  comps - pointer to integer array
   .  sm - sparse matrix

   DESCRIPTION:
   Computes the sparse matrix form of an offset
   array. Positive numbers in the comps-array are transfered, negative
   ones mean non-existing entries in the sparse matrix. Equal positive
   numbers mean identified fields. It is assumed that there is enough
   space in the arrays of the sparse matrix.

   RETURN VALUE:
   INT
   .n 0 ok
   .n 1 offset too large (increase MAX_MAT_COMP and recompile)

   SEE ALSO:
   SM2Array
   D*/
/****************************************************************************/

INT Array2SM (SHORT nr, SHORT nc, const SHORT *comps, SPARSE_MATRIX *sm)
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
/*D
   String2SMArray - Transforms a string to a SM array

   SYNOPSIS:
   INT String2SMArray (SHORT n, char *str, SHORT *comps)

   PARAMETERS:
   .  n - size (i.e. rows*cols)
   .  str - pointer to char array consisting of [*0a-z]
   .  comps - pointer to integer array

   DESCRIPTION:
   Transforms a string to a SM array. * means a non-zero entry, 0
   means a zero entry. a-z are used for identification of positions
   (they get the same offset). It is assumed that there is enough
   space in the sparse matrix array.

   RETURN VALUE:
   INT
   .n 0 ok
   .n 1 wrong format

   SEE ALSO:
   SM2Array
   D*/
/****************************************************************************/

INT String2SMArray (SHORT n, char *str, SHORT *comps)
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
/*D
   SM_Compute_Reduced_Size - Computes the reduced size of a sparse matrix

   SYNOPSIS:
   INT SM_Compute_Reduced_Size (SPARSE_MATRIX *sm)

   PARAMETERS:
   .  sm - sparse matrix

   DESCRIPTION:
   Computes the reduced size of a sparse matrix, i.e. takes possible
   identification into account.
   Warning: This is an O(N^2)-algorithm, but N should be moderate
   in comparison with the number of MATRIX structures anyhow.

   RETURN VALUE:
   INT
   .n positive: reduced size
   .n negative: an error occured
   D*/
/****************************************************************************/

INT SM_Compute_Reduced_Size (SPARSE_MATRIX *sm)
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
/*D
   SM_Compute_Reduced_Offsets - Computes the reduced offset field of a sparse matrix

   SYNOPSIS:
   INT SM_Compute_Reduced_Offsets (SPARSE_MATRIX *sm)

   PARAMETERS:
   .  sm - sparse matrix

   DESCRIPTION:
   Computes the reduced offset field of a sparse matrix, i.e. counting
   identified only once. The field to take the reduced offsets should
   be large enough (i.e. apply SM_Compute_Reduced_Size before)!

   Warning: This is an O(N^2)-algorithm, but N should be moderate
   in comparison with the number of MATRIX structures anyhow.

   RETURN VALUE:
   INT
   .n positive: number of reduced offsets
   .n negative: an error occured
   D*/
/****************************************************************************/

INT SM_Compute_Reduced_Offsets (SPARSE_MATRIX *sm, SHORT *reduced_offsets)
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
/*D
   SM_Compare - Compares two sparse matrices

   SYNOPSIS:
   INT SM_Compare (SPARSE_MATRIX *sm1, SPARSE_MATRIX *sm2)

   PARAMETERS:
   .  sm1 - sparse matrix
   .  sm2 - sparse matrix

   DESCRIPTION:
   Compares two sparse matrices

   RETURN VALUE:
   INT
   .n 0 ok
   .n 1 nrows not equal
   .n 2 ncols not equal
   .n 3 total number of nonzeros not equal
   .n 5 rows not equally long
   .n 6 column indices do not agree
   .n 7 offset structure is not compatible
   D*/
/****************************************************************************/

INT SM_Compare (SPARSE_MATRIX *sm1, SPARSE_MATRIX *sm2)
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
/*D
   SM_Compute_Diff_From_Offset - Computes a  ptrdiff_t-form of the sparse matrix offsets

   SYNOPSIS:
   INT SM_Compute_Diff_From_Offset (INT N, SHORT *offset,
                                 SHORT *start_off, ptrdiff_t Diff)

   PARAMETERS:
   .  N - size of offset-field
   .  offset - pointer to offset field
   .  Diff - pointer to ptrdiff_t field

   DESCRIPTION:
   Writes out the offset array in a ptrdiff_t-form, which should
   allow faster access for the (sparse) blas routines. The Diff-field
   should have place for N entries! Also the start_offset is returned,
   even if it might be easily obtained.

   RETURN VALUE:
   INT
   .n 0 ok
   .n -1 error
   D*/
/****************************************************************************/

INT SM_Compute_Diff_From_Offset (INT N, SHORT *offset, ptrdiff_t *Diff)
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
/*D
   SM_Compute_yDiff_From_Offset  - Computes a ptrdiff_t-form for a vector.

   SYNOPSIS:
   INT SM_Compute_yDiff_From_Offset (INT N, SHORT *col_ind, SHORT *cmp_off,
                                     ptrdiff_t *Diff)

   PARAMETERS:
   .  N - size of col_ind-field
   .  col_ind - pointer to col_ind field
   .  cmp_off - pointer to y's offset field
   .  Diff - pointer to ptrdiff_t field

   .  sm - sparse matrix
   .  compoff - component offsets
   .  startOff - start offset
   .  Diff - beginning of ptrdiff_t array

   DESCRIPTION:
   Writes out a ptrdiff_t-field of for fast access of the components
   of a vector inside a sparse matrix-vector multiplication.  The
   Diff-field should have place for N entries!  The cmp_off field
   should be a VECDATA_DESC-offset-field for some type.

   RETURN VALUE:
   INT
   .n 0 ok
   .n -1 an error ocurred
   D*/
/****************************************************************************/

INT SM_Compute_yDiff_From_Offset (INT N, SHORT *col_ind, SHORT *cmp_off,
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

/****************************************************************************/
/*D
   SM_Decompose_LR_pivot - LR decomposes a sparse matrix

   SYNOPSIS:
   INT SM_Decompose_LR_pivot (const SPARSE_MATRIX sm, DOUBLE *LR, int *pivot)

   associated functions:
   INT Decompose_LR_pivot (INT n, DOUBLE *mat, int *pivot)
   INT Solve_LR (const DOUBLE *LR, const int *pivot, DOUBLE *x, const DOUBLE *b)

   PARAMETERS:
   .  sm - sparse matrix
   .  LR - DOUBLE field of size nr^2
   .  pivot - SHORT field of size nr

   DESCRIPTION:
   Computes an LR decomposition with pivoting from the sparse matrix sm.
   sm must be quadratic, LR and pivot field must be long enough!

   RETURN VALUE:
   INT
   .n 0 ok
   .n -1 error
   D*/
/****************************************************************************/

INT Decompose_LR_pivot (int n, DOUBLE *mat, int *pivot)
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

INT Solve_LR (int n, const DOUBLE *LR, const int *pivot, DOUBLE *x, const DOUBLE *b)
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

INT SM_Decompose_LR_pivot (const SPARSE_MATRIX *sm, DOUBLE *values,
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
