// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_sp.c														*/
/*																			*/
/* Purpose:   support for sparse matrices in block compressed row storage	*/
/*																			*/
/* Author:	  Peter Bastian					                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   28 Jan 1996 Begin												*/
/*            02 Apr 1996 new memory allocation strategy					*/
/*            01 Oct 1997 redesign											*/
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
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "amg_header.h"
#include "amg_low.h"
#include "amg_sp.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#undef DEBUG

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

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*D
   AMG_NewVector - allocate a vector

   SYNOPSIS:
   AMG_VECTOR *AMG_NewVector (AMG_NODES *nd, int blocksize, char *name);

   PARAMETERS:
   .  n - number of blocks in this vector
   .  b - number of doubles per block, ie size in doubles is n*b
   .  name - name of this vector

   DESCRIPTION:
   This function allocates data (a double array) for a vector.

   RETURN VALUE:
   .n AMG_NULL memory overflow
   .n valid memory address else

   D*/
/****************************************************************************/

AMG_VECTOR *AMG_NewVector (int n, int b, char *name)
{
  AMG_VECTOR *new;
  double *values;

  /* allocate vector structure */
  new = AMG_Malloc(sizeof(AMG_VECTOR));
  if (new==NULL) return(AMG_NULL);

  /* allocate data */
  values = AMG_Malloc(n*b*sizeof(double));
  if (values==NULL) return(AMG_NULL);

  /* fill data structure */
  AMG_VECTOR_N(new) = n;
  AMG_VECTOR_B(new) = b;
  strncpy(AMG_VECTOR_NAME(new),name,AMG_NAME_SIZE-1);
  AMG_VECTOR_X(new) = values;

  return(new);
}

/****************************************************************************/
/*D
   AMG_NewMatrix - allocate a sparse matrix

   SYNOPSIS:
   AMG_MATRIX *AMG_NewMatrix (int n, int m, int b, int nonzeros, char *name);

   PARAMETERS:
   .  n,m - dimension of the matrix on blocks
   .  b - size of small blocks (are assumed to be quadratic)
   .  nonzeros - maximum number of nonzero blocks in this matrix
   .  system_as_scalar - system to be treated as scalar, see below
   .  name - name of this matrix

   DESCRIPTION:
   This function allocates data (a double array) for a matrix
   and corresponding integer arrays for block compressed row storage.

   Systems can be treated in two ways. First in a point block sense by setting
   b=system size and system_as_scalar=1. Each entry of the matrix is then a block
   matrix of dimension b*b. Or a system can be entered as individual scalar entries
   by setting b=1 and system_as_scalar=system size. Then IT IS ASSUMED THAT ROW NUMBER
   i BELONGS To COMPONENT i%system_as_scalar !

   EITHER b or system_as_scalar have to be 1 !

   RETURN VALUE:
   .n AMG_NULL memory overflow
   .n valid memory address else

   D*/
/****************************************************************************/

AMG_MATRIX *AMG_NewMatrix (int n, int b, int nonzeros, int system_as_scalar, char *name)
{
  AMG_MATRIX *new;
  double *a;
  int *r,*j;
  int i,k;

  if (b!=1 && system_as_scalar!=1)
  {
    AMG_Print("b or system_as_scalar must be 1\n");
    return(AMG_NULL);
  }

  /* allocate matrix structure */
  new = AMG_Malloc(sizeof(AMG_MATRIX));
  if (new==NULL) return(AMG_NULL);

  /* allocate data */
  a = AMG_Malloc(nonzeros*b*b*sizeof(double));
  if (a==NULL) return(AMG_NULL);

  /* col index array */
  j = AMG_Malloc(nonzeros*sizeof(int));
  if (j==NULL) return(AMG_NULL);

  /* row entry array */
  r = AMG_Malloc(n*sizeof(int));
  if (r==NULL) return(AMG_NULL);

  /* fill data structure */
  strncpy(AMG_MATRIX_NAME(new),name,AMG_NAME_SIZE-1);
  AMG_MATRIX_N(new) = n;
  AMG_MATRIX_B(new) = b;
  AMG_MATRIX_BB(new) = b*b;
  AMG_MATRIX_SAS(new) = system_as_scalar;
  AMG_MATRIX_NONZEROS(new) = nonzeros;
  AMG_MATRIX_CONNECTIONS(new) = 0;
  AMG_MATRIX_RA(new) = r;
  AMG_MATRIX_JA(new) = j;
  AMG_MATRIX_A(new) = a;

  /* initialize data structure */
  for (i=0; i<n; i++) r[i] = -1;
  for (k=0; k<nonzeros; k++) j[k] = -1;
  for (k=0; k<nonzeros*b*b; k++) a[k] = 0.0;

  return(new);
}

/****************************************************************************/
/*D
   AMG_CopyMatrix - make an identical copy of a matrix

   SYNOPSIS:
   AMG_MATRIX *AMG_CopyMatrix (AMG_MATRIX *A)

   PARAMETERS:
   .  A - matrix to copy

   DESCRIPTION:
   This function makes an identical copy of a matrix.

   RETURN VALUE:
   .n AMG_NULL error, no memory
   .n valid pointer else

   D*/
/****************************************************************************/

AMG_MATRIX *AMG_CopyMatrix (AMG_MATRIX *A, char *name)
{
  AMG_MATRIX *B;
  int i,k,n,*ra,*ja,*rb,*jb,nonzeros,bb;
  double *a, *b;

  /* allocate a matrix with same size */
  B = AMG_NewMatrix(A->n,A->b,A->nonzeros,A->system_as_scalar,name);
  if (B==AMG_NULL) return(B);

  /* copy data structure */
  ra = A->ra; ja = A->ja; a = A->a;
  rb = B->ra; jb = B->ja; b = B->a;
  nonzeros = A->nonzeros; bb = A->bb; n = A->n;
  for (i=0; i<n; i++) rb[i] = ra[i];
  for (k=0; k<nonzeros; k++) jb[k] = ja[k];
  for (k=0; k<nonzeros*bb; k++) b[k] = a[k];

  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_SetRowLength - initialize the length of a row

   SYNOPSIS:
   int AMG_SetRowLength (AMG_MATRIX *A, int i, int l);

   PARAMETERS:
   .  A - pointer to matrix
   .  i - row number
   .  l - number of nonzeros in row i, includes diagonal

   DESCRIPTION:
   Initialize the number of nonzeros needed in row i of the matrix.
   This function must be called once for each row in consecutive order
   (ie starting with row 0 up to row n-1). All entries in a row must
   be initialized subsequently with AMG_InsertEntry or AMG_InsertValues
   otherwise errors will occur when doing computations with the matrix.

   RETURN VALUE:
   .n AMG_FATAL memory overflow
   .n AMG_OK operation was ok

   D*/
/****************************************************************************/

int AMG_SetRowLength (AMG_MATRIX *A, int i, int l)
{
  if (i==0)
  {
    A->ra[0] = 0;                                               /* this is the first row                        */
    A->ja[0] = l;
    A->ra[1] = l;
    A->connections+=l;
  }
  else
  {
    if (A->ra[i]<0)                                                     /* check if i-1 has been set already*/
      return(AMG_FATAL);
    A->ja[A->ra[i]] = l;                                        /* store length of row i                        */
    A->connections+=l;
    if (i+1<A->n) {                                                     /* check if one more row                        */
      A->ra[i+1] = A->ra[i]+l;                                  /* set entry point for next row         */
      if (A->ra[i+1]>=A->nonzeros)                      /* check overflow					*/
        return(AMG_FATAL);
    }
  }

  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_FindEntry - retrieve a link from pattern

   SYNOPSIS:
   int AMG_FindEntry (AMG_MATRIX *A, int i, int j);

   PARAMETERS:
   .  A - matrix to search
   .  i,j - row and col indices

   DESCRIPTION:
   This function returns the index of the block i,j in the a array. The function
   can only be applied after the length of the given row has been initialized.

   RETURN VALUE:
   .n -1 entry not found
   .n valid index else

   D*/
/****************************************************************************/

int AMG_FindEntry (AMG_MATRIX *A, int i, int j)
{
  int k,*ra,*ja;

  /* check range */
  if ((i<0)||(i>=A->n)) return(-1);
  if ((j<0)||(j>=A->n)) return(-1);

  /* fast access */
  ra = A->ra; ja = A->ja;

  /* check initialization */
  if (ra[i]<0) return(-1);

  /* diagonal element ? */
  if (i==j) return(ra[i]);

  /* search row wise */
  for (k=ra[i]+1; k<ra[i]+ja[ra[i]]; k++)
    if (ja[k]==j) return(k);

  /* not found exit */
  return(-1);
}

/****************************************************************************/
/*D
   AMG_Insert - fill sparse matrix

   SYNOPSIS:
   int AMG_InsertEntry (AMG_MATRIX *A, int i, int j);
   int AMG_InsertValues (AMG_MATRIX *A, int i, int j, double *aij);
   int AMG_AddValues (AMG_MATRIX *A, int i, int j, double *aij);

   PARAMETERS:
   .  A - matrix to insert link
   .  i,j - row and col indices
   .  aij - array of values for block aij

   DESCRIPTION:
   These functions insert matrix elements into the sparse matrix. The length
   of the row must have been set before by the AMG_SetRowLength function.
   Then AMG_InsertEntry extends only the structure of the matrix,
   AMG_InsertValues extends the structure if necessary and copies the
   given values of a block, AMG_AddValues adds the values to the matrix block.

   RETURN VALUE:
   .n -1 if not enough memory or other error
   .n index of the block in the a array

   D*/
/****************************************************************************/

int AMG_InsertEntry (AMG_MATRIX *A, int i, int j)
{
  int k,*ra,*ja,found;

  /* check range */
  if ((i<0)||(i>=A->n)) return(-1);
  if ((j<0)||(j>=A->n)) return(-1);

  /* fast access */
  ra = A->ra; ja = A->ja;

  /* check initialization */
  if (ra[i]<0) return(-1);

  /* find or create new entry */
  found = -1;
  if (i==j)
    found = ra[i];
  else
  {
    for (k=ra[i]+1; k<ra[i]+ja[ra[i]]; k++)
    {
      if (ja[k]==j) {                                   /* already there */
        found=k;
        break;
      }
      if (ja[k]<0) {                                    /* assign new entry */
        ja[k]=j;
        found=k;
        break;
      }
    }
  }

  return(found);
}

int AMG_InsertValues (AMG_MATRIX *A, int i, int j, double *aij)
{
  int k,found;
  double *a;

  /* get entry */
  found = AMG_InsertEntry(A,i,j);
  if (found<0) return(found);

  /* copy values */
  a = (A->a)+(found*A->bb);
  for (k=0; k<A->bb; k++) {
    *a = aij[k];
    a++;
  }

  return(found);
}

int AMG_AddValues (AMG_MATRIX *A, int i, int j, double *aij)
{
  int k,found;
  double *a;

  /* get entry */
  found = AMG_InsertEntry(A,i,j);
  if (found<0) return(found);

  /* copy values */
  a = (A->a)+(found*A->bb);
  for (k=0; k<A->bb; k++) {
    *a += aij[k];
    a++;
  }

  return(found);
}

/****************************************************************************/
/*D
   AMG_PrintVector - print entries of a list of (block) vectors

   SYNOPSIS:
   int AMG_PrintVector (int k, AMG_VECTOR **vlist, char *text);

   PARAMETERS:
   .  k - number of vectors in vlist
   .  vlist - array of pointers to vectors
   .  text - title for printout

   DESCRIPTION:
   This function lists the values of a list of vectors. All vectors must
   share the same nodes and the same blocksize. One line of output
   is printed for one component in all vectors.

   RETURN VALUE:
   .n AMG_FATAL vectors not compatible, too many components
   .n AMG_OK else

   D*/
/****************************************************************************/

int AMG_PrintVector (int k, AMG_VECTOR **vlist, char *text)
{
  int block,blocksize,n;
  int component,i;
  char line[128];

  /* check if not too many */
  if (k>8) return(AMG_FATAL);
  n = AMG_VECTOR_N(vlist[0]);
  blocksize = AMG_VECTOR_B(vlist[0]);

  /* title */
  AMG_Print("------------------------------------------------------------------------\n");
  AMG_Print(text); AMG_Print("\n");
  AMG_Print("------------------------------------------------------------------------\n");

  /* now the loop */
  for (block=0; block<n; block++)
  {
    if (block%AMG_LINES_PER_PAGE==0)
    {
      sprintf(line,"%5s.%1s","BLOCK","C");
      AMG_Print(line);
      for (i=0; i<k; i++)
      {
        sprintf(line,"  %12s",AMG_VECTOR_NAME(vlist[i]));
        AMG_Print(line);
      }
      AMG_Print("\n");
    }
    for (component=0; component<blocksize; component++)
    {
      if (component==0) {
        sprintf(line,"%5d.",block); AMG_Print(line);
      } else {
        sprintf(line,"     ."); AMG_Print(line);
      }
      sprintf(line,"%1d",component); AMG_Print(line);
      for (i=0; i<k; i++)
      {
        sprintf(line,"  %12.4le",AMG_VECTOR_ENTRY(vlist[i],block,component));
        AMG_Print(line);
      }
      AMG_Print("\n");
    }
  }
  return(AMG_OK);
}


/****************************************************************************/
/*D
   AMG_PrintMatrix - print entries of a (block) matrix rowwise

   SYNOPSIS:
   int AMG_PrintMatrix (AMG_MATRIX *A, char *text)

   PARAMETERS:
   .  A - matrix to print
   .  text - title for printout

   DESCRIPTION:
   This function lists the contents of a matrix, i.e. values and structure.

   RETURN VALUE:
   .n AMG_OK

   D*/
/****************************************************************************/

int AMG_PrintMatrix (AMG_MATRIX *A, char *text)
{
  int n,b,i,k,c,*ra,*ja;
  double *a,*aa;
  char line[128];

  /* title */
  AMG_Print("------------------------------------------------------------------------\n");
  AMG_Print(AMG_MATRIX_NAME(A)); AMG_Print(": "); AMG_Print(text); AMG_Print("\n");
  AMG_Print("------------------------------------------------------------------------\n");

  n = A->n; b = A->b;
  ra = A->ra; ja = A->ja; a = A->a;

  if (b==1)       /* scalar case */
  {
    for (i=0; i<n; i++)
    {
      sprintf(line,"\nR %4d ",i); AMG_Print(line);
      sprintf(line,"[%4d:%12.4le] ",i,a[ra[i]]); AMG_Print(line);
      for (k=1; k<ja[ra[i]]; k++)
      {
        if (k%AMG_COLS_PER_LINE==0) AMG_Print("\n       ");
        sprintf(line,"[%4d:%12.4le] ",ja[k+ra[i]],a[k+ra[i]]);
        AMG_Print(line);
      }
      AMG_Print("\n");
    }
  }
  else
  {
    for (i=0; i<n; i++)
    {
      sprintf(line,"R %4d ",i); AMG_Print(line);

      sprintf(line,"[%4d:",i); AMG_Print(line);
      aa = a+(b*b*ra[i]);
      for (c=0; c<b*b; c++) {
        sprintf(line," %12.4le",i,aa[c]); AMG_Print(line);
      }
      AMG_Print("]\n");
      for (k=ra[i]; k<ra[i]+ja[ra[i]]; k++)
      {
        AMG_Print("       ");
        sprintf(line,"[%4d:",ja[k]); AMG_Print(line);
        aa = a+(b*b*k);
        for (c=0; c<b*b; c++) {
          sprintf(line," %12.4le",i,aa[c]); AMG_Print(line);
        }
        AMG_Print("]\n");
      }
    }
  }

  return(AMG_OK);
}
