// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  domain.h														*/
/*																			*/
/* Purpose:   interface to domain description                               */
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __DOMAIN__
#define __DOMAIN__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __SWITCH__
#include "switch.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* compile time constants */
#define BVPD_MAX_COEFFFCTS              10
#define BVPD_MAX_USERFCTS               10

/* boundary segment types */
#define PERIODIC                                1
#define NON_PERIODIC                    2
#define FREE                                    3


/* some function formats */
typedef INT (*ConfigProcPtr)(INT argc, char **argv);
typedef INT (*CoeffProcPtr)(COORD *, DOUBLE *);
typedef INT (*UserProcPtr)(DOUBLE *, DOUBLE *);


/* macros for BVPDescriptor */
#define BVPD_NAME(d)                    ((d).name)
#define BVPD_MIDPOINT(d)                ((d).midpoint)
#define BVPD_RADIUS(d)                  ((d).radius)
#define BVPD_NCORNERS(d)                ((d).nCorners)
#define BVPD_CONVEX(d)                  ((d).convex)
#define BVPD_NPATCHES(d)                ((d).nPatches)
#define BVPD_NSUBDOM(d)                 ((d).nSubDomains)
#define BVPD_CONFIG(d)                  ((d).ConfigProc)
#define BVPD_NCOEFFF(d)                 ((d).numOfCoeffFct)
#define BVPD_NUSERF(d)                  ((d).numOfUserFct)

/* macros for PatchDescriptor */
#define PATCH_LEFT(p)                   ((p).left)
#define PATCH_RIGHT(p)                  ((p).right)
#define PATCH_TYPE(p)                   ((p).type)
#define PATCH_N(p)                              ((p).n)
#define PATCH_ID(p)                             ((p).id)
#define PATCH_RES(p)                    ((p).resolution)
#define PATCH_CID(p,i)                  ((p).pc[i].id)
#define PATCH_LCVECT(p,i)               ((p).pc[i].lcoord)

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct BVPDescriptor
{
  /* general part */
  char name[NAMELEN];

  /* domain part */
  COORD midpoint[DIM];
  COORD radius;
  INT nCorners;
  INT convex;
  INT nPatches;
  INT nSubDomains;

  /* problem part */
  ConfigProcPtr ConfigProc;
  INT numOfCoeffFct;
  INT numOfUserFct;
};

struct patch_corner
{
  INT id;
  COORD lcoord[DIM_OF_BND];
};

struct PatchDescriptor
{
  INT left, right;
  INT type;
  INT n;
  INT id;
  INT resolution;
  struct patch_corner pc[CORNERS_OF_BND_SEG];
};

typedef struct BVPDescriptor BVP_DESC;
typedef struct PatchDescriptor PATCH_DESC;

/* dummies for ug */
typedef char BVP;
typedef char PATCH;


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* function declarations: IF-functions defined in the module				*/
/*																			*/
/****************************************************************************/

/* functionallity for BVPs (BoundaryValueProblems) */
BVP               *GetBVP                                                       (char *name);
INT        BVP_GetBVPDesc                               (BVP *theBVP, BVP_DESC *theBVPDesc);
PATCH     *BVP_GetFirstPatch                            (BVP *theBVP);
PATCH     *BVP_GetNextPatch                             (BVP *theBVP, PATCH *thePatch);
INT                BVP_GetCoeffFct                                      (BVP *theBVP, INT n, CoeffProcPtr *CoeffFct);
INT                BVP_GetUserFct                                       (BVP *theBVP, INT n, UserProcPtr *UserFct);

/* functionallity for PATCHes */
INT                Patch_GetPatchDesc                           (PATCH *thePatch, PATCH_DESC *thePatchDesc);
INT        Patch_global2local                   (PATCH *thePatch, COORD *global, COORD *local);
INT        Patch_local2global                   (PATCH *thePatch, COORD *local, COORD *global);
INT                Patch_local2bndcond                          (PATCH *thePatch, COORD *local, DOUBLE *value, INT *type);

/* miscellanious */
INT                InitDom                                                      (void);

/****************************************************************************/
/*																			*/
/* function declarations: derived function									*/
/*																			*/
/****************************************************************************/

PATCH      *Patch_GetPatchByID                          (BVP *theBVP, INT id);
INT             Patch_GetPatchID                                (PATCH *thePatch);
#endif
