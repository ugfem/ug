// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include "../main/defs.h"

/* Allocates a double vector with range [nl..nh]. */
/* After Numerical Recipies, p. 706 */
double *mkvec(nl,nh)
int nl, nh;
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */
	double *v;

	v=(double *) (MEM_OK = smalloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!MEM_OK) return;
	return(v-nl);
}

/* Free a double vector with range [nl..nh]. */
void frvec(v,nl)
double *v;
int nl;
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
	

	sfree((char *)(v+nl));
	v = NULL;
}
