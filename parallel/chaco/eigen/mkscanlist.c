// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

struct scanlink *mkscanlist(depth)
int depth;
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */
	struct scanlink *prevlnk;
	struct scanlink *newlnk;
	int i;


	prevlnk = (struct scanlink *) (MEM_OK = smalloc(sizeof(struct scanlink));
    if (!MEM_OK) return;
	prevlnk->pntr = NULL;
	newlnk = prevlnk; /* in case the list is one long */
	for (i=1; i<=(depth-1); i++) {
		newlnk = (struct scanlink *) (MEM_OK = smalloc(sizeof(struct scanlink));
        if (!MEM_OK) return; 
		newlnk->pntr = prevlnk;
		prevlnk = newlnk;
	}
	return(newlnk);
}
