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

static BndPropProcPtr Prop;
static BndGeomProcPtr Geom;

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
  theBVPDesc = (BVP_DESC *) theBVP;

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

static int FindSegment (int n, BP **bp)
{
  int i,j,k;

  for (i=0; i<bp[0]->n; i++) {
    for (k=1; k<n; k++) {
      for (j=0; j<bp[k]->n; j++)
        if (bp[k]->segment[j] == bp[0]->segment[i])
          break;
      if (j == bp[k]->n)
        break;
    }
    if (k == n)
      return(bp[0]->segment[i]);
  }
  printf("gen.c: inconsistent lines in geometry\n");

  for (k=0; k<n; k++) {
    printf("k %d n %d: ",k,bp[k]->n);
    for (j=0; j<bp[k]->n; j++)
      printf("%d: ",bp[k]->segment[j]);
    printf("\n");
  }

  ASSERT(0);
}

static int CountCommonSegments (BP *b0, BP *b1)
{
  int i,j;
  int n = 0;

  for (i=0; i<b0->n; i++)
  {
    int s = b0->segment[i];

    for (j=0; j<b1->n; j++)
      if (b1->segment[j] == s)
        n++;
  }
  return(n);
}

static int CommonSegments (BP *b0, BP *b1, int *S)
{
  int i,j;
  int n = 0;

  for (i=0; i<b0->n; i++)
  {
    int s = b0->segment[i];

    for (j=0; j<b1->n; j++)
      if (b1->segment[j] == s)
        S[n++] = s;
  }
  return(n);
}

BNDS* BNDP_CreateBndS (HEAP *Heap, BNDP **theBndP, INT n)
{
  BS *bs = (BS *) GetFreelistMemory(Heap,sizeof(BS));
  BP **bp = (BP **) theBndP;
  INT j,k;

  ASSERT(bs != NULL);

  bs->n = n;
  bs->id = 0;
  bs->segment = 0;
  bs->segment = FindSegment(n,bp);
  if (Prop == NULL)
    bs->property = 0;
  else
    bs->property = Prop(bs->segment);
  for (j=0; j<n; j++)
    for (k=0; k<DIM; k++)
      bs->x[j][k] = bp[j]->x[k];

  return((BNDS *)bs);
}

BNDP* BNDP_CreateBndP (HEAP *Heap, BNDP *theBndP0, BNDP *theBndP1,
                       DOUBLE lcoord)
{
  BP *bp;
  BP *b0 = (BP *) theBndP0;
  BP *b1 = (BP *) theBndP1;
  INT k;
  INT n = CountCommonSegments(b0,b1);

  ASSERT(n>0);
  bp = (BP *) GetFreelistMemory(Heap,sizeof(BP)+(n-1)*sizeof(int));
  ASSERT(bp != NULL);
  bp->n = CommonSegments(b0,b1,bp->segment);
  bp->id = 0;
  for (k=0; k<DIM; k++)
    bp->x[k] = (1.0 - lcoord) * b0->x[k] + lcoord * b1->x[k];
  if (Geom != NULL)
    Geom(bp->segment[0],b0->x,b1->x,lcoord,bp->x);
  if (Prop == NULL)
    bp->property = 0;
  else
    bp->property = Prop(bp->segment[0]);

  return((BNDP *)bp);
}

INT BNDP_Dispose (HEAP *Heap, BNDP *theBndP)
{
  return(PutFreelistMemory(Heap,theBndP,sizeof(BP)));
}

INT BNDP_SaveBndP (BNDP *theBndP)
{
  BP *p = (BP *)theBndP;
  INT j;
  int iList[3+MAX_SEGMENTS];
  double dList[DIM];

  iList[0] = p->id;
  iList[1] = p->property;
  iList[2] = p->n;
  if (Bio_Write_mint(3,iList)) return (1);
  for (j=0; j<p->n; j++)
    iList[j] = p->segment[j];
  if (Bio_Write_mint(p->n,iList)) return (1);

  for (j=0; j<DIM; j++)
    dList[j] = p->x[j];
  if (Bio_Write_mdouble(DIM,dList)) return (1);

  return(0);
}

INT BNDP_SaveBndP_Ext (BNDP *theBndP)
{
  BP *p = (BP *)theBndP;
  INT j;
  int iList[3+MAX_SEGMENTS];
  double dList[DIM];

  iList[0] = p->id;
  iList[1] = p->property;
  iList[2] = p->n;
  if (Bio_Write_mint(3,iList)) return (1);
  for (j=0; j<p->n; j++)
    iList[j] = p->segment[j];
  if (Bio_Write_mint(p->n,iList)) return (1);

  for (j=0; j<DIM; j++)
    dList[j] = p->x[j];
  if (Bio_Write_mdouble(DIM,dList)) return (1);

  return(0);
}

BNDP *BNDP_LoadBndP (BVP *theBVP, HEAP *Heap)
{
  BP *p;
  INT j;
  int iList[3+MAX_SEGMENTS];
  double dList[DIM];

  if (Bio_Read_mint(3,iList)) return (NULL);
  p->id = iList[0];
  p->property = iList[1];
  p->n = iList[2];
  p = (BP *) GetFreelistMemory(Heap,sizeof(BP)+(p->n-1)*sizeof(int));
  if (Bio_Read_mint(p->n,iList)) return (NULL);
  for (j=0; j<p->n; j++)
    p->segment[j] = iList[j];
  if (Bio_Read_mdouble(DIM,dList)) return (NULL);
  for (j=0; j<DIM; j++)
    p->x[j] = dList[j];

  return((BNDP *)p);
}

BNDP *BNDP_LoadBndP_Ext (void)
{
  BP *p;
  INT j;
  int iList[3+MAX_SEGMENTS];
  double dList[DIM];

  if (Bio_Read_mint(3,iList)) return (NULL);
  p->id = iList[0];
  p->property = iList[1];
  p->n = iList[2];
  p  = (BP *)malloc(sizeof(BP)+(p->n-1)*sizeof(int));
  if (Bio_Read_mint(p->n,iList)) return (NULL);
  for (j=0; j<p->n; j++)
    p->segment[j] = iList[j];
  if (Bio_Read_mdouble(DIM,dList)) return (NULL);
  for (j=0; j<DIM; j++)
    p->x[j] = dList[j];

  return((BNDP *)p);
}

INT BNDS_Global (BNDS *theBndS, DOUBLE *local, DOUBLE *global)
{
  BS *bs = (BS *) theBndS;
  INT k;

  if (bs->n == 2)
  {
    for (k=0; k<DIM; k++)
      global[k] = (1.0 - local[0]) * bs->x[0][k]
                  + local[0] * bs->x[1][k];
    if (Geom != NULL)
      Geom(bs->segment,bs->x[0],bs->x[1],local[0],global);
  }
  else if (bs->n == 3)
  {
    if (Geom == NULL) {
      for (k=0; k<DIM; k++)
        global[k] = (1.0 - local[0] - local[1]) * bs->x[0][k]
                    + local[0] * bs->x[1][k] + local[1] * bs->x[2][k];
    }
    else
    {
      DOUBLE pos[DIM];
      DOUBLE a = local[0] + local[1];
      DOUBLE b;

      if (a == 0.0) {
        for (k=0; k<DIM; k++)
          global[k] = bs->x[0][k];
        return(0);
      }
      b = local[1] / a;

      for (k=0; k<DIM; k++)
        pos[k] = (1.0 - b) * bs->x[1][k] + b * bs->x[2][k];
      Geom(bs->segment,bs->x[1],bs->x[2],b,pos);
      for (k=0; k<DIM; k++)
        global[k] = (1.0 - a) * bs->x[0][k] + a * pos[k];
      Geom(bs->segment,bs->x[0],pos,a,global);
    }
  }
  else if (bs->n == 4)
  {
    if (Geom == NULL) {
      for (k=0; k<DIM; k++)
        global[k] = (1.0 - local[0]) * (1.0 - local[1]) * bs->x[0][k]
                    + local[0] * (1.0 - local[1]) * bs->x[1][k]
                    + (1.0 - local[0]) * local[1] * bs->x[2][k]
                    + local[0] * local[1] * bs->x[3][k];
    }
    else
    {
      DOUBLE pos0[DIM];
      DOUBLE pos1[DIM];

      for (k=0; k<DIM; k++)
        pos0[k] = (1.0 - local[0]) * bs->x[0][k]
                  + local[0] * bs->x[1][k];
      Geom(bs->segment,bs->x[0],bs->x[1],local[0],pos0);
      for (k=0; k<DIM; k++)
        pos1[k] = (1.0 - local[0]) * bs->x[3][k]
                  + local[0] * bs->x[2][k];
      Geom(bs->segment,bs->x[3],bs->x[2],local[0],pos1);
      for (k=0; k<DIM; k++)
        global[k] = (1.0 - local[1]) * pos0[k] + local[1] * pos1[k];
      Geom(bs->segment,pos0,pos1,local[1],global);
    }
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
  BS *bs = (BS *) theBndS;

  *id = bs->in;
  *nbid = bs->out;
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
  bp->n = 1;
  bp->segment[0] = bs->segment;
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
  {0,0,0,0,0,0},
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
      int k,n;

      if (G->P[i].segment <G->nS)
        n = 1;
      else
        n = G->L[G->P[i].segment-G->nS].n;
      ASSERT(n>0);
      G->P[i].bp = (BP *)
                   GetFreelistMemory(Heap,sizeof(BP)+(n-1)*sizeof(int));
      ASSERT(G->P[i].bp != NULL);
      G->P[i].bp->id = G->P[i].id;
      G->P[i].bp->n = n;
      if (n == 1)
        G->P[i].bp->segment[0] = G->P[i].segment;
      else
        for (k=0; k<n; k++)
          G->P[i].bp->segment[k] = G->L[G->P[i].segment-G->nS].S[k];
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
          int in = G->S[G->C[i].S[j]].in;
          int out = G->S[G->C[i].S[j]].out;

          G->C[i].bs[j] = (BS *) GetFreelistMemory(Heap,sizeof(BS));
          ASSERT(G->C[i].bs[j] != NULL);
          G->C[i].bs[j]->id = G->F[G->C[i].F[j]].id;
          G->C[i].bs[j]->in = G->C[i].subdomain;
          if (G->C[i].subdomain == in)
            G->C[i].bs[j]->out = out;
          else
            G->C[i].bs[j]->out = in;
          G->C[i].bs[j]->segment = G->C[i].S[j];
          G->C[i].bs[j]->property = G->S[G->C[i].S[j]].property;
          G->C[i].bs[j]->n = ncfaces[G->C[i].n][j];
          for (k=0; k<G->C[i].bs[j]->n; k++) {
            int l;
            for (l=0; l<DIM; l++)
              G->C[i].bs[j]->x[k][l] =
                G->P[G->C[i].P[faces[G->C[i].n][j][k]]].x[l];
          }
        }
        else
          G->C[i].bs[j] = NULL;
      }
    }

  Prop = G->Prop;
  Geom = G->Geom;

  return (0);
}

INT InitDom (void)
{
  Prop = NULL;
  Geom = NULL;

  return (0);
}
