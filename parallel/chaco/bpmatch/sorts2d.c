// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"../main/defs.h"
#include	"../main/params.h"
   

void sorts2d(vals, indices, nvtxs)
/* Sort the lists needed to find the splitter. */
double *vals[4][MAXSETS];	/* lists of values to sort */
int *indices[4][MAXSETS];	/* indices of sorted lists */
int nvtxs;			/* number of vertices */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   int *space;			/* space for mergesort routine */
   int *temp[4];		/* place holders for indices */
   int nlists = 4;		/* number of directions to sort */
   int i;			/* loop counter */
   
   void mergesort();

   space = (int *) (MEM_OK = smalloc((unsigned) nvtxs*sizeof(int));
   if (!MEM_OK) return;

   for (i=0; i<nlists; i++) {
      temp[i] = (int *) (MEM_OK = smalloc((unsigned) nvtxs*sizeof(int));
      if (!MEM_OK) return;
   }

   mergesort(vals[0][1], nvtxs, temp[0], space);
   mergesort(vals[0][2], nvtxs, temp[1], space);
   mergesort(vals[0][3], nvtxs, temp[2], space);
   mergesort(vals[1][2], nvtxs, temp[3], space);

   sfree((char *) space);

   indices[0][1] = indices[1][0] = indices[2][3] = indices[3][2] = temp[0];
   indices[0][2] = indices[2][0] = indices[1][3] = indices[3][1] = temp[1];
   indices[0][3] = indices[3][0] = temp[2];
   indices[1][2] = indices[2][1] = temp[3];
}
