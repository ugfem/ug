// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      domain.h                                                      */
/*                                                                          */
/* Purpose:   standard header file template                                 */
/*                                                                          */
/* Author:      klaus Johannsen                                             */
/*              Institut fuer Computeranwendungen III                       */
/*              Universitaet Stuttgart                                      */
/*              Pfaffenwaldring 27                                          */
/*              70550 Stuttgart                                             */
/*              email: ug@ica3.uni-stuttgart.de                             */
/*                                                                          */
/* History:   18.06.96 begin, ug version 3.2                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DOMAIN__
#define __DOMAIN__

#include <stdio.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

#ifndef __HEAPS__
#include "heaps.h"
#endif

#ifdef __MWCW__
#include "MWCW.cmdlinedefs"
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

#ifdef _2
#ifdef _3
#error ****    define EITHER dimension _2 OR _3       ****
#endif
#define __TWODIM__
#define DIM 2
#define DIM_OF_BND 1
#define DOM_PARAM_OFFSET DIM
#endif

#ifdef _3
#define __THREEDIM__
#define DIM 3
#define DIM_OF_BND 2
#define DOM_PARAM_OFFSET DIM
#endif

#ifndef _2
#ifndef _3
#error ****    define at least dimension two OR three        ****
#endif
#endif

/* object identification */
#define BSOBJ 8                                                 /* boundary surface                 */
#define BPOBJ 9                                                 /* boundary point object			*/

/* boundary types */
#define FIXED         0
#define FREE          1
#define PERIODIC      2
#define NON_PERIODIC  3

/* function formats */
typedef INT (*ConfigProcPtr)(INT argc, char **argv);
typedef INT (*CoeffProcPtr)(DOUBLE *, DOUBLE *);
typedef INT (*UserProcPtr)(DOUBLE *, DOUBLE *);

/* macros for BVPDescriptor */
#define BVPD_NAME(d)         ((d).name)
#define BVPD_MIDPOINT(d)     ((d).midpoint)
#define BVPD_RADIUS(d)       ((d).radius)
#define BVPD_CONVEX(d)       ((d).convex)
#define BVPD_NSUBDOM(d)      ((d).nSubDomains)
#define BVPD_CONFIG(d)       ((d).ConfigProc)
#define BVPD_NCOEFFF(d)      ((d).numOfCoeffFct)
#define BVPD_NUSERF(d)       ((d).numOfUserFct)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef void *BVP;                 /* structure handled by domain module    */
typedef void *BNDS;                /* structure handled by domain module    */
typedef void *BNDP;                /* structure handled by domain module    */

struct BVP_Descriptor
{
  /* general part */
  char name[NAMELEN];                /* name of the BVP                     */

  /* domain part */
  DOUBLE midpoint[DIM];               /* sphere in which the domain lies     */
  DOUBLE radius;
  INT convex;                        /* 1 if domain is convex, 0 if not     */
  INT nSubDomains;                   /* nb. of subdomains,
                                                                                exterior not counted                */

  /* problem part */
  ConfigProcPtr ConfigProc;          /* configuration function              */
  INT numOfCoeffFct;                 /* nb. of coefficient functions        */
  INT numOfUserFct;                  /* nb. of user functions               */
};
typedef struct BVP_Descriptor BVP_DESC;

struct mesh
{
  INT nBndP;                         /* nb. of boundary points              */
  BNDP **theBndPs;                                       /* list of boundary points	            */
  INT nInnP;                         /* nb. of inner nodes                  */
  DOUBLE **Position;                  /* positions of inner nodes            */
  INT nSubDomains;                   /* nb. of subdomains                   */
  INT *nSides;                       /* nb. of boundary sides per subdomain */
  INT **Side_corners;                /* nb. of side corners                 */
  INT ***Side_corner_ids;                /* corner ids                          */
  INT *nElements;                    /* nb. of element corners              */
  INT **Element_corners;             /* nb. of element corners              */
  INT ***Element_corner_ids;         /* nb. of side corners                 */
  INT ***nbElements;                 /* nb. of side corners                 */
};
typedef struct mesh MESH;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP

void DomInitParallel     (INT TypeBndP, INT TypeBndS);
void DomHandlerInit      (INT handlerSet);

void BElementXferBndS    (BNDS **bnds, int n, int proc, int prio);
void BElementGatherBndS  (BNDS **bnds, int n, int cnt, char *data);
void BElementScatterBndS (BNDS **bnds, int n, int cnt, char *data);

void BVertexXferBndP     (BNDP *bndp, int proc, int prio);
void BVertexGatherBndP   (BNDP *bndp, int cnt, char *data);
void BVertexScatterBndP  (BNDP **bndp, int cnt, char *data);

#endif

/* functions for BVP */
BVP        *BVP_GetFirst          (void);
BVP        *BVP_GetNext           (BVP *theBVP);
INT         BVP_Save              (BVP *theBVP, char *name, char *mgname, HEAP *theHeap, INT argc, char **argv);
BVP        *BVP_Load              (char *name, INT argc, char **argv);
BVP        *BVP_GetByName         (char *name);
BVP        *BVP_Init              (char *name, HEAP *Heap, MESH *Mesh);
INT         BVP_Dispose           (BVP *theBVP);
INT         BVP_SetBVPDesc        (BVP *theBVP, BVP_DESC *theBVPDesc);
INT         BVP_SetCoeffFct       (BVP *theBVP, INT n, CoeffProcPtr *CoeffFct);
INT         BVP_SetUserFct        (BVP *theBVP, INT n, UserProcPtr *UserFct);
INT             BVP_Check                         (BVP *aBVP);

/* functions called by script commands */
BNDP*           BVP_InsertBndP            (HEAP *Heap, BVP *theBVP, INT argc, char **argv);
INT         BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size);
MESH       *BVP_GenerateMesh      (HEAP *Heap, BVP *theBVP, INT argc, char **argv);

/* functions for BNDP */
INT         BNDP_Global           (BNDP *theBndP, DOUBLE *global);
INT         BNDP_BndCond          (BNDP *theBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value, INT *type);
INT         BNDP_BndPDesc         (BNDP *theBndP, INT *move);
BNDS*       BNDP_CreateBndS       (HEAP *Heap, BNDP **theBndP, INT n);
BNDP*       BNDP_CreateBndP       (HEAP *Heap, BNDP *theBndP0, BNDP *theBndP1, DOUBLE lcoord);
INT         BNDP_Dispose          (HEAP *Heap, BNDP *theBndP);
INT         BNDP_SaveBndP         (BNDP *theBndP);
BNDP       *BNDP_LoadBndP         (BVP *theBVP, HEAP *Heap);

/* functions for BNDS */
INT         BNDS_Global           (BNDS *theBndS, DOUBLE *local, DOUBLE *global);
INT         BNDS_BndCond          (BNDS *theBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value, INT *type);
INT         BNDS_BndSDesc         (BNDS *theBndS, INT *left, INT *right);
BNDP*       BNDS_CreateBndP       (HEAP *Heap, BNDS *theBndS, DOUBLE *local);
INT         BNDS_Dispose          (HEAP *Heap, BNDS *theBndS);

/* miscellanious */
INT         InitDom               (void);

#endif
