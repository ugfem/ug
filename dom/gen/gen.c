// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gen.c                                                                                                         */
/*																			*/
/* Purpose:   general ug domain description                                                             */
/*																			*/
/* Author:	  Christian Wieners	                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Dez 8 1999 begin                                                                          */
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

/* standard C library */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/* low modules */
#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "bio.h"
#include "misc.h"
#include "fileopen.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"
#include "scan.h"
#include "evm.h"

/* dev modules */
#include "ugdevices.h"

/* domain module */
#include "gen.h"
#include "domain.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAX_LEN      128

#define Malloc(m)    GetTmpMem(Heap,(m),HKey)

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* functions for the domain interface                                           */
/*																			*/
/****************************************************************************/

BVP *BVP_GetNext (BVP *theBVP)
{
  return(NULL);
}

INT BVP_Save (BVP *theBVP, char *name,
              char *mgname, HEAP *theHeap, INT argc, char **argv)
{
  return(1);
}

BVP *BVP_Load (char *name, INT argc, char **argv)
{
  return(NULL);
}

BVP *BVP_GetByName (char *name)
{
  return(NULL);
}

BVP *BVP_Init (char *name, HEAP *Heap, MESH *Mesh, INT MarkKey)
{
  return(NULL);
}

INT BVP_Dispose (BVP *theBVP)
{
  return(1);
}

INT BVP_SetBVPDesc (BVP *theBVP, BVP_DESC *theBVPDesc)
{
  GEOMETRY *G = (GEOMETRY *) theBVP;
  INT i,k;

  theBVPDesc->convex = FALSE;
  theBVPDesc->nSubDomains = G->sd;
  theBVPDesc->nDomainParts = 1;
  for (i=0; i<=theBVPDesc->nSubDomains; i++)
    theBVPDesc->s2p[i] = 0;
  for (k=0; k<DIM; k++)
    theBVPDesc->midpoint[k] = 0.0;
  for (i=0; i<G->nP; i++)
    if (G->P[i].bnd == 1)
      for (k=0; k<DIM; k++)
        theBVPDesc->midpoint[k] += G->P[i].x[k];
  ASSERT(G->nBP > 0);
  for (k=0; k<DIM; k++)
    theBVPDesc->midpoint[k] *= 1.0 / G->nBP;
  theBVPDesc->radius = 0.0;
  for (i=0; i<G->nP; i++)
    if (G->P[i].bnd == 1)
      for (k=0; k<DIM; k++)
        theBVPDesc->radius =
          MAX(theBVPDesc->radius,
              ABS(theBVPDesc->midpoint[k]-G->P[i].x[k]));
  theBVPDesc->numOfCoeffFct = 0;
  theBVPDesc->numOfUserFct = 0;

  return (0);
}

INT BVP_SetCoeffFct (BVP *theBVP, INT n, CoeffProcPtr *CoeffFct)
{
  return(1);
}

INT BVP_SetUserFct (BVP *theBVP, INT n, UserProcPtr *UserFct)
{
  return(1);
}

INT BVP_Check (BVP *aBVP)
{
  return(1);
}

BNDP* BVP_InsertBndP (HEAP *Heap, BVP *theBVP, INT argc, char **argv)
{
  return(NULL);
}

INT BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size)
{
  return(1);
}

MESH *BVP_GenerateMesh (HEAP *Heap, BVP *aBVP, INT argc, char **argv,
                        INT MarkKey)
{
  return(NULL);
}

INT BNDP_Global (BNDP *theBndP, DOUBLE *global)
{
  BP *bp = (BP *)theBndP;
  int k;

  for (k=0; k<DIM; k++)
    global[k] = bp->x[k];

  return(0);
}

INT BNDP_Move (BNDP *aBndP, const DOUBLE global[])
{
  return(1);
}

INT BNDP_BndCond (BNDP *theBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value,
                  INT *type)
{
  return(1);
}

INT BNDP_BndPDesc (BNDP *theBndP, INT *move, INT *part)
{
  *move = 0;
  *part = 0;

  return(0);
}

INT BNDP_BndEDesc (BNDP *theBndP0, BNDP *theBndP1, INT *part)
{
  *part = 0;

  return(0);
}

BNDS* BNDP_CreateBndS (HEAP *Heap, BNDP **theBndP, INT n)
{
  BS *bs = (BS *) GetFreelistMemory(Heap,sizeof(BS));
  INT k;

  ASSERT(bs != NULL);

  bs->n = n;
  bs->id = 0;
  bs->segment = 0;
  bs->property = 0;
  for (k=0; k<n; k++) {
    memcpy(&(bs->bp[k]),theBndP,sizeof(BP));
    bs->segment = MAX(bs->segment,bs->bp[k].segment);
    bs->property = MAX(bs->property,bs->bp[k].property);
  }
  return((BNDS *)bs);
}

BNDP* BNDP_CreateBndP (HEAP *Heap, BNDP *theBndP0, BNDP *theBndP1,
                       DOUBLE lcoord)
{
  BP *bp = (BP *) GetFreelistMemory(Heap,sizeof(BP));
  BP *b0 = (BP *) theBndP0;
  BP *b1 = (BP *) theBndP1;
  INT k;

  ASSERT(bp != NULL);

  bp->id = 0;
  for (k=0; k<DIM; k++)
    bp->x[k] = (1.0 - lcoord) * b0->x[k] + lcoord * b1->x[k];
  bp->segment = MAX(b0->segment,b1->segment);
  bp->property = MAX(b0->property,b1->property);

  return((BNDP *)bp);
}

INT BNDP_Dispose (HEAP *Heap, BNDP *theBndP)
{
  return(PutFreelistMemory(Heap,theBndP,sizeof(BP)));
}

INT BNDP_SaveBndP (BNDP *theBndP)
{
  return(1);
}

INT BNDP_SaveBndP_Ext (BNDP *theBndP)
{
  return(1);
}

BNDP *BNDP_LoadBndP (BVP *theBVP, HEAP *Heap)
{
  return(NULL);
}

BNDP *BNDP_LoadBndP_Ext (void)
{
  return(NULL);
}

INT BNDS_Global (BNDS *theBndS, DOUBLE *local, DOUBLE *global)
{
  BS *bs = (BS *) theBndS;
  INT k;

  if (bs->n == 2)
  {
    BP *b0 = &(bs->bp[0]);
    BP *b1 = &(bs->bp[1]);

    for (k=0; k<DIM; k++)
      global[k] = (1.0 - local[0]) * b0->x[k] + local[0] * b1->x[k];
  }
  else if (bs->n ==3)
  {
    BP *b0 = &(bs->bp[0]);
    BP *b1 = &(bs->bp[1]);
    BP *b2 = &(bs->bp[2]);

    for (k=0; k<DIM; k++)
      global[k] = (1.0 - local[0] - local[1]) * b0->x[k]
                  + local[0] * b1->x[k] + local[1] * b2->x[k];
  }
  else if (bs->n == 4)
  {
    BP *b0 = &(bs->bp[0]);
    BP *b1 = &(bs->bp[1]);
    BP *b2 = &(bs->bp[2]);
    BP *b3 = &(bs->bp[3]);

    for (k=0; k<DIM; k++)
      global[k] = (1.0 - local[0]) * (1.0 - local[1]) * b0->x[k]
                  + local[0] * (1.0 - local[1]) * b1->x[k]
                  + (1.0 - local[0]) * local[1] * b2->x[k]
                  + local[0] * local[1] * b3->x[k];
  }

  return(0);
}

INT BNDS_BndCond (BNDS *theBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value,
                  INT *type)
{
  return(1);
}

INT BNDS_BndSDesc (BNDS *theBndS, INT *id, INT *nbid, INT *part)
{
  *id = 1;
  *nbid = 0;
  *part = 0;

  return(0);
}

BNDP* BNDS_CreateBndP (HEAP *Heap, BNDS *theBndS, DOUBLE *local)
{
  BS *bs = (BS *) theBndS;
  BP *bp = (BP *) GetFreelistMemory(Heap,sizeof(BP));

  ASSERT(bp != NULL);

  BNDS_Global(theBndS,local,bp->x);
  bp->id = bs->id;
  bp->segment = bs->segment;
  bp->property = bs->property;

  return((BNDP *)bp);
}

BVP *BVP_GetFirst (void)
{
  return(NULL);
}

INT BNDS_Dispose (HEAP *Heap, BNDS *theBndS)
{
  return(PutFreelistMemory(Heap,theBndS,sizeof(BS)));
}

#ifdef __TWODIM__
static int nfaces[5] = {0,0,0,3,4};
static int ncfaces[5][4] = {
  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0},
  {2,2,2,0},
  {2,2,2,2}
};
static int faces[5][4][2] = {
  {{0,0},{0,0},{0,0},{0,0}},
  {{0,0},{0,0},{0,0},{0,0}},
  {{0,0},{0,0},{0,0},{0,0}},
  {{0,1},{1,2},{2,0},{0,0}},
  {{0,1},{1,2},{2,3},{3,4}}
};
#else
static int nfaces[9] = {0,0,0,0,4,5,5,0,6};
static int ncfaces[9][6] = {
  {0,0,0,0,0,0},
  {0,0,0,0,0,0},
  {0,0,0,0,0,0},
  {0,0,0,0,0,0},
  {3,3,3,3,0,0},
  {4,3,3,3,3,0},
  {3,4,4,4,3,0},
  {4,4,4,4,4,4}
};
static int faces[9][6][4] = {
  {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
  {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
  {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
  {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
  {{0,2,1,9},{1,2,3,9},{0,3,2,9},{0,1,3,9},{0,0,0,0},{0,0,0,0}},
  {{0,3,2,1},{0,1,4,9},{1,2,4,9},{2,3,4,9},{3,0,4,9},{0,0,0,0}},
  {{0,2,1,9},{0,1,4,3},{1,2,5,4},{2,0,3,5},{3,4,5,9},{0,0,0,0}},
  {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
  {{0,3,2,1},{0,1,4,5},{1,2,5,6},{2,3,6,7},{3,0,7,4},{4,5,6,7}}
};
#endif

INT InitGeometry (HEAP *Heap, GEOMETRY *G)
{
  INT i;

  for (i=0; i<G->nP; i++)
    if (G->P[i].bnd == 1)
    {
      int k;

      G->P[i].bp = (BP *) GetFreelistMemory(Heap,sizeof(BP));
      ASSERT(G->P[i].bp != NULL);
      G->P[i].bp->id = G->P[i].id;
      G->P[i].bp->segment = 0;
      G->P[i].bp->property = 0;
      for (k=0; k<DIM; k++)
        G->P[i].bp->x[k] = G->P[i].x[k];
    }
  for (i=0; i<G->nC; i++)
    if (G->C[i].bnd == 1)
    {
      int j;

      for (j=0; j<nfaces[G->C[i].n]; j++) {
        if (G->C[i].C[j] < 0)
        {
          int k;

          G->C[i].bs[j] = (BS *) GetFreelistMemory(Heap,sizeof(BS));
          ASSERT(G->C[i].bs[j] != NULL);
          G->C[i].bs[j]->id = G->P[i].id;
          G->C[i].bs[j]->segment = G->C[i].S[j];
          G->C[i].bs[j]->property = G->S[G->C[i].S[j]].property;
          G->C[i].bs[j]->n = ncfaces[G->C[i].n][j];
          for (k=0; k<G->C[i].bs[j]->n; k++)
            memcpy(&(G->C[i].bs[j]->bp[k]),
                   G->P[G->C[i].P[faces[G->C[i].n][j][k]]].bp,
                   sizeof(BP));

          for (k=0; k<G->C[i].bs[j]->n; k++) {
            PRINTDEBUG(dom,5,("C %d bs[%d]= %d bp %d\n",
                              G->C[i].id,j,
                              G->C[i].bs[j]->id,
                              G->C[i].bs[j]->bp[k].id));
          }

        }
        else
          G->C[i].bs[j] = NULL;
      }
    }

  return (0);
}

INT InitDom (void)
{
  return (0);
}
