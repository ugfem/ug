// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "heaps.h"
#include "cmdint.h"
/*#include "../util/smalloc.h"*/


/* RCS_ID
   $Header$
 */

#define Heap HEAP

/*#define __DEBUG__*/

#define smalloc(n)      GetMem(heap,(n),FROM_BOTTOM))
#define sfree(p)     DisposeMem(heap,(p))
/* #define smalloc(n)	cmalloc(n))
 #define sfree(p)   cfree(p) */

#define MED(x)  ((x)/2+(x)%2)
#define BI_1     0
#define BI_2     1
#define  QUAD_1   2
#define  QUAD_2   1
#define  QUAD_3   0
#define  QUAD_4   3

#ifndef max
#define max(A, B)       ((A) > (B) ? (A) : (B))
#endif
#ifndef min
#define min(A, B)       ((A) < (B) ? (A) : (B))
#endif
#ifndef sign
#define sign(A)         ((A) <  0  ? -1  :  1)
#endif
#define absval(A)       ((A) <  0  ? -(A) : (A))
#define TRUE            1
#define FALSE           0

/* Define constants that are needed in various places */
#define PI                       3.141592653589793238462643383279
#define TWOPI   6.283185307179586
#define HALFPI  1.570796326794896
