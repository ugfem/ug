// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Dynamically allocate a 2 dimensional array. */

#include "../main/defs.h"

double *array_alloc_2D(dim1, dim2, size)
int dim1;		/* size of first dimension */
int dim2;		/* size of second dimension */
unsigned int size;	/* size of array elements */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   int total;		/* Total size of the array */
   int aligned_dim;	/* dim1 or dim1+1 to ensure data alignement */
   int offset;		/* offset of array elements */
   double *field;	/* The multi-dimensional array */
   char **ptr;		/* Pointer offset */
   char *data;		/* Data offset*/
   int	j;		/* loop counter */

   aligned_dim =  (dim1 % 2) ? dim1 + 1 : dim1;
   offset = aligned_dim * sizeof(void *);
   total = offset + dim1 * dim2 * size;
   field = (double *) (MEM_OK = smalloc((unsigned) total);
   if (!MEM_OK) return;

   ptr = (char **) (field);
   data = (char *) (field);
   data += offset;
   for (j=0; j<dim1; j++) {
      ptr[j] = data + j * size * dim2;
   }

   return((double *) field);
}
