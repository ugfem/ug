// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/structs.h"
#include "../main/defs.h"

/* Allocate space for new orthlink.
--------------------------------------*/
struct orthlink *makeorthlnk()
{
        extern Heap   *heap;     /* pointer to heap of multigrid */
        extern double *MEM_OK;   /* variable for memory overflow exeception */
        struct orthlink *newlnk;

        newlnk =  (struct orthlink *) (MEM_OK = smalloc(sizeof(struct orthlink));
        if (!MEM_OK) return;
        return(newlnk);
}
