// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*****************************************************************************
* File:      lgm_accel.c                                                    *
* Purpose:   Accelerate geometry access by using bounding box trees         *
*                                                                           *
* Author:	  O. Sterz                                                       *
*                                                                           *
* History:   Nov 2002                                                       *
* Remarks:   Each surface gets a pointer to a tree with triangles           *
*****************************************************************************/

/*****************************************************************************
* include files                                                             *
*            system include files                                           *
*            application include files                                      *
*****************************************************************************/
#include <assert.h>
#include "compiler.h"   /* for INT etc.                              */
#include "debug.h"      /* for IFDEBUG() and  MIN(), MAX() by misc.h */
#include "ugdevices.h"  /* for UserWrite() and friends               */
#include "heaps.h"      /* for GetTmpMem() etc.                      */
#include "lgm_domain.h" /* for LGM_SURFACE etc.                      */
#include "bbtree.h"     /* for BBT_NewTree() etc.                    */
#include "lgm_accel.h"
#include "tree.h"
#include "fifo.h"

#include "namespace.h"

USING_UG_NAMESPACE
USING_UGDIM_NAMESPACE


#ifdef LGM_ACCELERATE
/*****************************************************************************
* defines in the following order                                            *
*        compile time constants defining static data size (i.e. arrays)     *
*        other constants                                                    *
*        macros                                                             *
*****************************************************************************/

/*****************************************************************************
* PRIVATE data structures (global to this source file)                      *
*****************************************************************************/


/*****************************************************************************
* PRIVATE global variables (global to this source file)                     *
*****************************************************************************/


/*****************************************************************************
* PRIVATE function declarations/definitions (used in this source file)      *
*****************************************************************************/
/*---------------------------------------------------------------------------*
* CreateTriangleBBox - creates a bounding box for an lgm triangle           *
*---------------------------------------------------------------------------*/
static BBT_BBOX *CreateTriangleBBox(HEAP *theHeap, LGM_TRIANGLE *triangle)
{
  DOUBLE *x, ll[DIM], ur[DIM];

  assert(DIM == 3);
  x = triangle->corner[0]->position;
  ll[0] = ur[0] = x[0];
  ll[1] = ur[1] = x[1];
  ll[2] = ur[2] = x[2];

  x = triangle->corner[1]->position;
  ll[0] = MIN(ll[0], x[0]);
  ur[0] = MAX(ur[0], x[0]);
  ll[1] = MIN(ll[1], x[1]);
  ur[1] = MAX(ur[1], x[1]);
  ll[2] = MIN(ll[2], x[2]);
  ur[2] = MAX(ur[2], x[2]);

  x = triangle->corner[2]->position;
  ll[0] = MIN(ll[0], x[0]);
  ur[0] = MAX(ur[0], x[0]);
  ll[1] = MIN(ll[1], x[1]);
  ur[1] = MAX(ur[1], x[1]);
  ll[2] = MIN(ll[2], x[2]);
  ur[2] = MAX(ur[2], x[2]);

  return BBT_NewBBox(theHeap, DIM, ll, ur, triangle);
}


/*****************************************************************************
* PUBLIC function definitions (declarations are in the include file)        *
*****************************************************************************/
/*---------------------------------------------------------------------------*/
/*D
   LGM_InitAcceleration - init the lgm geometry acces acceleration

   SYNOPSIS:
   LGM_InitAcceleration(LGM_SURFACE **sf, INT nsf)

   PARAMETERS:
   .  theHeap - pointer to heap
   .  sf      - array of pointer to lgm surfaces
   .  n       - #lgm surfaces

   DESCRIPTION:
   .  This fuction will be called from LGM_LoadDomain() in lgm_load.c to generate
   .  for each surface a tree with bounding boxes containing the triangles of the
   .  lgm surface.

   RETURN VALUE:
   .  INT value 0: o.k. 1: error
   D*/
/*---------------------------------------------------------------------------*/
INT NS_DIM_PREFIX LGM_InitAcceleration(HEAP *theHeap, LGM_SURFACE **sf, INT nsf)
{
  LGM_TRIANGLE *triangle;
  BBT_BBOX **bboxes;
  INT MarkKey, max_nsf, nt, ntsum, i, j;

  UserWriteF("Building %d trees to speed up geometry: ", nsf);
  IFDEBUG(LGMAccel,1)  UserWrite("\n"); ENDDEBUG

  /* Max number of triangles per surface */
    max_nsf = 0;
  for (i=0; i<nsf; i++)
    max_nsf = MAX(max_nsf, LGM_SURFACE_NTRIANGLE(sf[i]));

  /* Allocate array of pointers to bboxes */
  MarkTmpMem(theHeap, &MarkKey);
  bboxes = (BBT_BBOX **) GetTmpMem(theHeap, sizeof(BBT_BBOX *)*max_nsf, MarkKey);
  if (bboxes == NULL) return(1);

  ntsum = 0;
  for (i=0; i<nsf; i++)
  {
    IFDEBUG(LGMAccel,1)  UserWriteF("[%3d", i); ENDDEBUG
      nt = LGM_SURFACE_NTRIANGLE(sf[i]);
    for(j=0; j<nt; j++)
    {
      triangle = LGM_SURFACE_TRIANGLE(sf[i], j);
      bboxes[j] = CreateTriangleBBox(theHeap, triangle);
    }
    sf[i]->bbtree = BBT_NewTree(theHeap, bboxes, nt, DIM);
    if (sf[i]->bbtree == NULL) return 1;
    ntsum += nt;
    IFDEBUG(LGMAccel,1) UserWriteF(":%4d]\n", nt); ENDDEBUG
  }
  UserWriteF("%d triangles\n", ntsum);

  /* Free array of pointers to bboxes */
  ReleaseTmpMem(theHeap, MarkKey);
  return 0;
}
#endif
