// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  std_domain.c													*/
/*																			*/
/* Purpose:   standard ug domain description                                                            */
/*																			*/
/* Author:	  Klaus Johannsen / Christian Wieners	                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Feb 18 1996 begin, ug version 3.1								*/
/*            Sep 12 1996 ug version 3.4                                                                */
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
#include "defaults.h"
#include "general.h"
#include "debug.h"

/* dev modules */
#include "devices.h"

/* domain module */
#include "std_domain.h"
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

#define SMALL_DIFF   SMALL_C*100
#define RESOLUTION   100

#define DEFAULTDOMMEMORY 50000

#define STD_BVP_NAME(p)                         ENVITEM_NAME(p)
#define STD_BVP_MIDPOINT(p)                     ((p)->MidPoint)
#define STD_BVP_RADIUS(p)                       ((p)->radius)
#define STD_BVP_CONVEX(p)                       ((p)->domConvex)
#define STD_BVP_NCORNER(p)                      ((p)->nCorner)
#define STD_BVP_NSUBDOM(p)                      ((p)->nSubDomain)
#define STD_BVP_NPATCH(p)                       ((p)->nPatch)
#define STD_BVP_NCOEFFPROC(p)           ((p)->nCoeffFct)
#define STD_BVP_NUSERPROC(p)            ((p)->nUserFct)
#define STD_BVP_CONFIGPROC(p)           ((p)->ConfigProblem)
#define STD_BVP_COEFFPROC(p,i)          (CoeffProcPtr)((p)->CoeffF[i])
#define STD_BVP_USERPROC(p,i)           (UserProcPtr)((p)->UserF[i])

#define GetSTD_BVP(p)              ((STD_BVP *)(p))

#define V2_LINCOMB(a,A,b,B,C)              {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
                                            (C)[1] = (a)*(A)[1] + (b)*(B)[1];}

#define V2_EUKLIDNORM_OF_DIFF(A,B,b)    (b) = sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])));

#define V3_EUKLIDNORM_OF_DIFF(A,B,b)    (b) = (sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])+((A)[2]-(B)[2])*((A)[2]-(B)[2]))));

#ifdef __TWODIM__
#define V_DIM_EUKLIDNORM_OF_DIFF(A,B,b) V2_EUKLIDNORM_OF_DIFF(A,B,b)
#else
#define V_DIM_EUKLIDNORM_OF_DIFF(A,B,b) V3_EUKLIDNORM_OF_DIFF(A,B,b)
#endif

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

static INT theProblemDirID;             /* env type for Problem dir                     */
static INT theBdryCondVarID;            /* env type for Problem vars			*/

static INT theDomainDirID;                      /* env type for Domain dir				*/
static INT theBdrySegVarID;             /* env type for bdry segment vars		*/

static INT theBVPDirID;                         /* env type for BVP dir					*/

static STD_BVP *currBVP;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   CreateProblem -  Create a new PROBLEM structure

   SYNOPSIS:
   PROBLEM *CreateProblem (char *domain,char *name, int id,
   ConfigProcPtr config, int numOfCoefficients, CoeffProcPtr coeffs[],
   int numOfUserFct, UserProcPtr userfct[]);

   PARAMETERS:
   .  domain - name of the domain
   .  name - name of the problem
   .  id - identification number for the problem
   .  config - pointer to the configuration function
   .  numOfCoefficients - number of coefficient functions
   .  coeffs[] - pointer to coefficient functions
   .  numOfUserFct - number of user coefficient functions
   .  userfct[] - pointer to user coefficient functions

   DESCRIPTION:
   This function allocates and initializes a new PROBLEM structure in the /domains/domain directory.

   RETURN VALUE:
   PROBLEM *
   .n    pointer to new PROBLEM
   .n    NULL if out of memory.
   D*/
/****************************************************************************/

PROBLEM *CreateProblem (char *domain,char *name, int id, ConfigProcPtr config,
                        int numOfCoefficients, CoeffProcPtr coeffs[],
                        int numOfUserFct,          UserProcPtr userfct[])

{
  PROBLEM *newProblem;
  int i;

  if (ChangeEnvDir("/Domains")==NULL) return(NULL);
  if (ChangeEnvDir(domain)==NULL) return(NULL);

  /* allocate new problem structure */
  newProblem = (PROBLEM *) MakeEnvItem (name,theProblemDirID,sizeof(PROBLEM)+(numOfCoefficients+numOfUserFct-1)*sizeof(void*));
  if (newProblem==NULL) return(NULL);

  /* fill in data */
  newProblem->problemID = id;
  newProblem->ConfigProblem = config;
  newProblem->numOfCoeffFct = numOfCoefficients;
  newProblem->numOfUserFct  = numOfUserFct;
  for (i=0; i<numOfCoefficients; i++) newProblem->CU_ProcPtr[i] = (void*)(coeffs[i]);
  for (i=0; i<numOfUserFct; i++) newProblem->CU_ProcPtr[i+numOfCoefficients] = (void*)(userfct[i]);

  if (ChangeEnvDir(name)==NULL) return(NULL);
  UserWrite("problem "); UserWrite(name); UserWrite(" installed\n");

  return(newProblem);
}

/****************************************************************************/
/*
   GetProblem -  Get pointer to a problem structure

   SYNOPSIS:
   PROBLEM *GetProblem (const char * domain, const char *name);

   PARAMETERS:
   .  domain - name of the domain
   .  name - name of the problem

   DESCRIPTION:
   This function searches a problem structure in the /domains/domain directory.

   RETURN VALUE:
   PROBLEM *
   .n      pointer to PROBLEM
   .n      NULL if not found or error.
 */
/****************************************************************************/

static PROBLEM *GetProblem (const char * domain, const char *name)
{
  if (ChangeEnvDir("/Domains")==NULL) return(NULL);
  return((PROBLEM *) SearchEnv(name,domain,theProblemDirID,theDomainDirID));
}

/****************************************************************************/
/*D
   CreateBoundaryCondition - Allocate a new BNDCOND structure

   SYNOPSIS:
   BOUNDARY_CONDITION *CreateBoundaryCondition (char *name, INT id,
   BndCondProcPtr theBndCond, void *Data);

   PARAMETERS:
   .  name - name of the boundary condition
   .  id - identification number of condition
   .  theBndCond - the boundary conditions
   .  Data - user defined data

   DESCRIPTION:
   This function allocates and initializes a new BNDCOND structure in the previously allocated
   PROBLEM structure

   THE RETURN VALUE:
   BOUNDARY_CONDITION *
   .n      pointer to
   .n      0 if  out of memory.
   D*/
/****************************************************************************/

BOUNDARY_CONDITION *CreateBoundaryCondition (char *name, INT id, BndCondProcPtr theBndCond, void *Data)
{
  BOUNDARY_CONDITION *newBndCond;

  /* allocate boundary condition */
  newBndCond = (BOUNDARY_CONDITION *) MakeEnvItem (name,theBdryCondVarID,sizeof(BOUNDARY_CONDITION));
  if (newBndCond==NULL) return(NULL);

  /* fill in data */
  newBndCond->id = id;
  newBndCond->BndCond = theBndCond;
  newBndCond->data = Data;

  return(newBndCond);
}

/****************************************************************************/
/*
   GetNextBoundaryCondition - Allocate a new BNDCOND structure

   SYNOPSIS:
   static BOUNDARY_CONDITION *GetNextBoundaryCondition (BOUNDARY_CONDITION *theBCond);

   PARAMETERS:
   .  theBCond - the boundary condition

   DESCRIPTION:
   This function gets next BNDCOND structure

   THE RETURN VALUE:
   BOUNDARY_CONDITION *
   .n      pointer to BOUNDARY_CONDITION
   .n      NULL if not found or error.
 */
/****************************************************************************/
static BOUNDARY_CONDITION *GetNextBoundaryCondition (BOUNDARY_CONDITION *theBCond)
{
  ENVITEM *theItem;

  theItem = (ENVITEM *) theBCond;

  do
    theItem = NEXT_ENVITEM(theItem);
  while ((theItem!=NULL) && (ENVITEM_TYPE(theItem)!=theBdryCondVarID));

  return ((BOUNDARY_CONDITION *) theItem);
}

/****************************************************************************/
/*
   GetFirstBoundaryCondition - Get first BNDCOND structure of `theProblem`

   SYNOPSIS:
   static BOUNDARY_CONDITION *GetFirstBoundaryCondition (PROBLEM *theProblem);

   PARAMETERS:
   .  theProblem - pointer to PROBLEM

   DESCRIPTION:
   This function gets the first BNDCOND structure of a problem

   THE RETURN VALUE:
   BOUNDARY_CONDITION *
   .n      pointer to BOUNDARY_CONDITION
   .n      NULL if not found or error.
 */
/****************************************************************************/
static BOUNDARY_CONDITION *GetFirstBoundaryCondition (PROBLEM *theProblem)
{
  ENVITEM *theItem;

  theItem = ENVITEM_DOWN(theProblem);

  if (ENVITEM_TYPE(theItem)==theBdryCondVarID)
    return ((BOUNDARY_CONDITION *) theItem);
  else
    return (GetNextBoundaryCondition((BOUNDARY_CONDITION *) theItem));
}

/****************************************************************************/
/*D
   CreateDomain	- Create a new DOMAIN data structure

   SYNOPSIS:
   DOMAIN *CreateDomain (char *name, COORD *MidPoint, COORD radius,
   INT segments, INT corners, INT Convex);

   PARAMETERS:
   .  name - name of the domain
   .  MidPoint - coordinates of some inner point
   .  radius - radius of a circle, containing the domain
   .  segments - number of the boundary segments
   .  corners - number of corners
   .  Convex - 0, if convex, 1 - else

   DESCRIPTION:
   This function allocates and initializes a new DOMAIN data structure in the
   /domains directory in the environment.

   RETURN VALUE:
   DOMAIN *
   .n     pointer to a DOMAIN
   .n     NULL if out of memory.
   D*/
/****************************************************************************/

DOMAIN *CreateDomain (char *name, COORD *MidPoint, COORD radius, INT segments, INT corners, INT Convex)
{
  DOMAIN *newDomain;
  INT i;

  /* change to /domains directory */
  if (ChangeEnvDir("/Domains")==NULL)
    return (NULL);

  /* allocate new domain structure */
  newDomain = (DOMAIN *) MakeEnvItem (name,theDomainDirID,sizeof(DOMAIN));
  if (newDomain==NULL) return(NULL);

  /* fill in data */
  for( i = 0 ; i < DIM ; i++)
    newDomain->MidPoint[i] = MidPoint[i];
  newDomain->radius = radius;
  newDomain->numOfSegments = segments;
  newDomain->numOfCorners = corners;
  newDomain->domConvex = Convex;

  if (ChangeEnvDir(name)==NULL) return(NULL);
  UserWrite("domain "); UserWrite(name); UserWrite(" installed\n");

  return(newDomain);
}

/****************************************************************************/
/*
   GetDomain  - Get a pointer to a domain structure

   SYNOPSIS:
   DOMAIN *GetDomain (char *name);

   PARAMETERS:
   .  name - name of the domain

   DESCRIPTION:
   This function searches the environment for a domain with the name `name`
   and return a pointer to the domain structure.

   RETURN VALUE:
   DOMAIN *
   .n     pointer to a DOMAIN
   .n     NULL if not found or error.
 */
/****************************************************************************/

DOMAIN *GetDomain (char *name)
{
  return((DOMAIN *) SearchEnv(name,"/Domains",theDomainDirID,theDomainDirID));
}

/****************************************************************************/
/*D
   CreateBoundarySegment - Create a new BOUNDARY_SEGMENT

   SYNOPSIS:
   BOUNDARY_SEGMENT *CreateBoundarySegment (char *name
   INT left, INT right,INT id,INT type,INT res,INT *point,
   COORD *alpha,COORD *beta,BndSegFuncPtr BndSegFunc, void *data);

   PARAMETERS:
   .  name - name of the boundary segment
   .  left - id of left subdomain
   .  right - id of right subdomain
   .  id - id of this boundary segment
   .  type - type of the boundary segment
   .  res  - resolution of the boundary segment
   .  point - the endpoints of the boundary segment
   .  alpha - list where the parameter interval begins
   .  beta - list where the parameter interval ends
   .  BndSegFunc - function mapping parameters
   .  data - user defined space

   DESCRIPTION:
   This function allocates and initializes a new BOUNDARY_SEGMENT for the previously
   allocated DOMAIN structure.

   RETURN VALUE:
   BOUNDARY_SEGMENT *
   .n     Pointer to a BOUNDARY_SEGMENT
   .n     NULL if out of memory.
   D*/
/****************************************************************************/

BOUNDARY_SEGMENT *CreateBoundarySegment (char *name,
                                         INT left, INT right,INT id,INT type,INT res,INT *point,
                                         COORD *alpha,COORD *beta,BndSegFuncPtr BndSegFunc, void *data)
{
  BOUNDARY_SEGMENT *newSegment;
  INT i;

  /* allocate the boundary segment */
  newSegment = (BOUNDARY_SEGMENT *) MakeEnvItem (name,theBdrySegVarID,sizeof(BOUNDARY_SEGMENT));
  if (newSegment==NULL) return(NULL);

  /* fill in data */
  newSegment->left = left;
  newSegment->right = right;
  newSegment->id = id;
  newSegment->segType = type;
  for(i = 0 ; i < CORNERS_OF_BND_SEG ; i++)
    newSegment->points[i] = point[i];
  newSegment->resolution = res;
  for(i = 0 ; i < DIM_OF_BND ; i++)
  {
    newSegment->alpha[i] = alpha[i];
    newSegment->beta[i] = beta[i];
  }
  newSegment->BndSegFunc = BndSegFunc;
  newSegment->data = data;

  return(newSegment);
}

/****************************************************************************/
/*
   GetNextBoundarySegment - Get next boundary segment of a domain

   SYNOPSIS:
   static BOUNDARY_SEGMENT *GetNextBoundarySegment (BOUNDARY_SEGMENT *theBSeg);

   PARAMETERS:
   .  theBSeg - pointer to a boundary segment

   DESCRIPTION:
   This function gets the next boundary segment of a domain.

   RETURN VALUE:
   BOUNDARY_SEGMENT *
   .n     pointer to next BOUNDARY_SEGMENT
   .n     NULL if no more found.
 */
/****************************************************************************/

static BOUNDARY_SEGMENT *GetNextBoundarySegment (BOUNDARY_SEGMENT *theBSeg)
{
  ENVITEM *theItem;

  theItem = (ENVITEM *) theBSeg;

  do
    theItem = NEXT_ENVITEM(theItem);
  while ((theItem!=NULL) && (ENVITEM_TYPE(theItem)!=theBdrySegVarID));

  return ((BOUNDARY_SEGMENT*) theItem);
}

/****************************************************************************/
/*
   GetFirstBoundarySegment - Get first boundary segment of a domain

   SYNOPSIS:
   static BOUNDARY_SEGMENT *GetFirstBoundarySegment (DOMAIN *theDomain);

   PARAMETERS:
   .  theDomain - pointer to the domain structure

   DESCRIPTION:
   This function returns the first boundary segment of a domain.

   RETURN VALUE:
   BOUNDARY_SEGMENT *
   .n     pointer to a BOUNDARY_SEGMENT
   .n     NULL if not found or error.
 */
/****************************************************************************/

static BOUNDARY_SEGMENT *GetFirstBoundarySegment (DOMAIN *theDomain)
{
  ENVITEM *theItem;

  theItem = ENVITEM_DOWN(theDomain);

  if (ENVITEM_TYPE(theItem)==theBdrySegVarID)
    return ((BOUNDARY_SEGMENT *) theItem);
  else
    return (GetNextBoundarySegment((BOUNDARY_SEGMENT *) theItem));
}

/****************************************************************************/
/*

   (This function is here for upward compatibilty)

   SYNOPSIS:
   BOUNDARY_SEGMENT *CreateBoundarySegment2D (char *name, int left, int right,
        int id, int from, int to, int res, COORD alpha, COORD beta,
        BndSegFuncPtr BndSegFunc, void *data);


   PARAMETERS:
   .  name - name of the boundary segment
   .  left - id of the left subdomain
   .  right - id of the right subdomain
   .  id -
   .  from -
   .  to -
   .  res -
   .  alpha -
   .  beta -
   .  BndSegFunc -
   .  data -

   DESCRIPTION:
   This function defines CreateBoundarySegment2D for old style 2D definitions
   (they were not dimension independent)

   RETURN VALUE:
   BOUNDARY_SEGMENT *
   .n     pointer to
   .n     NULL if object not available.
 */
/****************************************************************************/

BOUNDARY_SEGMENT *CreateBoundarySegment2D (char *name, int left, int right,
                                           int id, int from, int to, int res, COORD alpha, COORD beta,
                                           BndSegFuncPtr BndSegFunc, void *data)
{
  INT pt[3];
  COORD alp[3],bet[3];

  pt[0] = from;
  pt[1] = to;
  alp[0] = alpha;
  bet[0] = beta;

  return(CreateBoundarySegment(name,left,right,id,NON_PERIODIC,
                               res,pt,alp,bet,BndSegFunc,data));
}

/****************************************************************************/
/*D
   CreateBVP - Create BoundaryValueProblem from Domain and Problem

   SYNOPSIS:
   BVP *CreateBVP (char *BVP, char *Domain, char *Problem);

   PARAMETERS:
   .  BVP - name of BVP to be created
   .  Domain - name of domain
   .  Problem - name of Problem

   DESCRIPTION:
   This function creates a BoundaryValueProblem from the already created
   structures domain and problem

   RETURN VALUE:
   BVP *
   .n     pointer BoundaryValueProblem
   .n     NULL if not creatable.
   D*/
/****************************************************************************/

BVP *CreateBoundaryValueProblem (char *BVPName, BndCondProcPtr theBndCond,
                                 int numOfCoeffFct, CoeffProcPtr coeffs[],
                                 int numOfUserFct, UserProcPtr userfct[])
{
  STD_BVP *theBVP;
  INT i,n;

  /* change to /BVP directory */
  if (ChangeEnvDir("/BVP")==NULL)
    return (NULL);

  /* allocate new domain structure */
  n = (numOfCoeffFct+numOfUserFct-1)*sizeof(void*);
  theBVP = (STD_BVP *) MakeEnvItem (BVPName,theBVPDirID,sizeof(STD_BVP)+n);
  if (theBVP==NULL)
    return(NULL);
  if (ChangeEnvDir(BVPName)==NULL)
    return(NULL);

  theBVP->numOfCoeffFct = numOfCoeffFct;
  theBVP->numOfUserFct = numOfUserFct;
  for (i=0; i<numOfCoeffFct; i++)
    theBVP->CU_ProcPtr[i] = (void*)(coeffs[i]);
  for (i=0; i<numOfUserFct; i++)
    theBVP->CU_ProcPtr[i+numOfCoeffFct] = (void*)(userfct[i]);

  theBVP->Domain = NULL;
  theBVP->Problem = NULL;
  theBVP->ConfigProc = STD_BVP_Configure;
  theBVP->GeneralBndCond = theBndCond;

  UserWriteF("BVP %s installed.\n",BVPName);

  return ((BVP*)theBVP);
}

BVP *CreateBVP (char *BVPName, char *DomainName, char *ProblemName)
{
  STD_BVP *theBVP;
  DOMAIN *theDomain;
  PROBLEM *theProblem;
  INT i,n;

  /* get domain and problem */
  theDomain = GetDomain(DomainName);
  if (theDomain==NULL)
    return (NULL);
  theProblem = GetProblem(DomainName,ProblemName);
  if (theProblem==NULL)
    return (NULL);

  /* change to /BVP directory */
  if (ChangeEnvDir("/BVP")==NULL)
    return (NULL);

  /* allocate new domain structure */
  n = (theProblem->numOfCoeffFct+theProblem->numOfUserFct-1)*sizeof(void*);
  theBVP = (STD_BVP *) MakeEnvItem (BVPName,theBVPDirID,sizeof(STD_BVP)+n);
  if (theBVP==NULL)
    return(NULL);
  if (ChangeEnvDir(BVPName)==NULL)
    return(NULL);
  for (i=0; i<theProblem->numOfCoeffFct; i++)
    theBVP->CU_ProcPtr[i] = theProblem->CU_ProcPtr[i];
  for (i=0; i<theProblem->numOfUserFct; i++)
    theBVP->CU_ProcPtr[i+theProblem->numOfCoeffFct] =
      theProblem->CU_ProcPtr[i+theProblem->numOfCoeffFct];
  theBVP->numOfCoeffFct = theProblem->numOfCoeffFct;
  theBVP->numOfUserFct = theProblem->numOfUserFct;

  theBVP->Domain = theDomain;
  theBVP->Problem = theProblem;

  /* fill in data of problem */
  theBVP->ConfigProc = theProblem->ConfigProblem;
  theBVP->GeneralBndCond = NULL;

  UserWriteF("BVP %s installed.\n",BVPName);

  return ((BVP*)theBVP);
}


/****************************************************************************/
/*D
   BVP_CreateMesh - sets a BNDP from command input

   SYNOPSIS:
   MESH *BVP_CreateMesh (HEAP *Heap, BVP *theBVP);

   PARAMETERS:
   .  theBVP - BVP structure
   .  i - number of boundary point
   .  theBndP - the BNDP to set

   DESCRIPTION:
   This function sets a BNDP from command input parameters

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

static INT CreateCornerPoints (HEAP *Heap, STD_BVP *theBVP, BNDP **bndp)
{
  BND_PS *ps;
  PATCH *p,*pp;
  INT i,j,m;

  for (i=0; i<theBVP->ncorners; i++)
  {
    p = theBVP->patches[i];
    m = POINT_PATCH_N(p);
    ps = (BND_PS *)GetFreelistMemory(Heap,sizeof(BND_PS)
                                     +(m-1)*sizeof(COORD_BND_VECTOR));
    ps->n = m;
    if (ps == NULL)
      return(1);
    ps->patch_id = PATCH_ID(p);


    for (j=0; j<m; j++)
    {
      pp = theBVP->patches[POINT_PATCH_PID(p,j)];


      PRINTDEBUG(dom,1,("cp i %d r %f %f %f %f\n",i,
                        PARAM_PATCH_RANGE(pp)[0][0],
                        PARAM_PATCH_RANGE(pp)[0][1],
                        PARAM_PATCH_RANGE(pp)[1][0],
                        PARAM_PATCH_RANGE(pp)[1][1]));

      switch (POINT_PATCH_CID(p,j))
      {
              #ifdef __TWODIM__
      case 0 :
        ps->local[j][0] = PARAM_PATCH_RANGE(pp)[0][0];
        break;
      case 1 :
        ps->local[j][0] = PARAM_PATCH_RANGE(pp)[1][0];
        break;
                          #endif
              #ifdef __THREEDIM__
      case 0 :
        ps->local[j][0] = PARAM_PATCH_RANGE(pp)[0][0];
        ps->local[j][1] = PARAM_PATCH_RANGE(pp)[0][1];
        break;
      case 1 :
        ps->local[j][0] = PARAM_PATCH_RANGE(pp)[1][0];
        ps->local[j][1] = PARAM_PATCH_RANGE(pp)[0][1];
        break;
      case 2 :
        ps->local[j][0] = PARAM_PATCH_RANGE(pp)[1][0];
        ps->local[j][1] = PARAM_PATCH_RANGE(pp)[1][1];
        break;
      case 3 :
        ps->local[j][0] = PARAM_PATCH_RANGE(pp)[0][0];
        ps->local[j][1] = PARAM_PATCH_RANGE(pp)[1][1];
        break;
                          #endif
      }
      PRINTDEBUG(dom,1,("mesh loc j %d pid %d cid %d loc %f %f\n",
                        j,
                        POINT_PATCH_PID(p,j),
                        POINT_PATCH_CID(p,j),
                        ps->local[j][0],ps->local[j][1]));
    }
    bndp[i] = (BNDP *)ps;
  }

  for (i=0; i<theBVP->ncorners; i++)
    PRINTDEBUG(dom,1,(" id %d\n",PATCH_ID(theBVP->patches[i])));

  return(0);
}

BVP *BVP_Init (char *name, HEAP *Heap, MESH *Mesh)
{
  STD_BVP *theBVP;
  DOMAIN *theDomain;
  PROBLEM *theProblem;
  BOUNDARY_SEGMENT *theSegment;
  BOUNDARY_CONDITION *theBndCond;
  PATCH **corners,**lines,**sides,*thePatch;
  INT i,j,k,n,m,maxSubDomains,ncorners,nlines,nsides;

  theBVP = (STD_BVP *)BVP_GetByName(name);
  if (theBVP == NULL)
    return(NULL);
  theDomain = theBVP->Domain;
  if (theDomain == NULL)
    return(NULL);
  theProblem = theBVP->Problem;

  /* fill in data of domain */
  for (i=0 ; i<DIM ; i++)
    theBVP->MidPoint[i] = theDomain->MidPoint[i];
  theBVP->radius = theDomain->radius;
  theBVP->domConvex = theDomain->domConvex;

  ncorners = theDomain->numOfCorners;
  nsides = theDomain->numOfSegments;

  /* create parameter patches */
  maxSubDomains = 0;
  sides = (PATCH **)GetTmpMem(Heap,nsides*sizeof(PATCH *));
  for (i=0; i<nsides; i++)
    sides[i] = NULL;

  if (sides == NULL)
    return(NULL);
  theBVP->nsides = nsides;
  for (theSegment=GetFirstBoundarySegment(theDomain); theSegment!=NULL;
       theSegment = GetNextBoundarySegment(theSegment))
  {
    if ((theSegment->id<0)||(theSegment->id>=nsides))
      return(NULL);
    thePatch=(PATCH *)GetFreelistMemory(Heap,sizeof(PARAMETER_PATCH));
    if (thePatch == NULL)
      return (NULL);
    PATCH_TYPE(thePatch) = PARAMETRIC_PATCH_TYPE;
    PATCH_ID(thePatch) = theSegment->id;
    PARAM_PATCH_LEFT(thePatch) = theSegment->left;
    PARAM_PATCH_RIGHT(thePatch) = theSegment->right;
    PARAM_PATCH_BC(thePatch) = NULL;
    PARAM_PATCH_BCD(thePatch) = NULL;
    for (i=0; i<2*DIM_OF_BND; i++)
      PARAM_PATCH_POINTS(thePatch,i) = theSegment->points[i];
    for (i=0; i<DIM_OF_BND; i++)
    {
      PARAM_PATCH_RANGE(thePatch)[0][i] = theSegment->alpha[i];
      PARAM_PATCH_RANGE(thePatch)[1][i] = theSegment->beta[i];
    }
    PARAM_PATCH_BS(thePatch) = theSegment->BndSegFunc;
    PARAM_PATCH_BSD(thePatch) = theSegment->data;
    maxSubDomains = MAX(maxSubDomains,theSegment->left);
    maxSubDomains = MAX(maxSubDomains,theSegment->right);
    sides[theSegment->id] = thePatch;
    PRINTDEBUG(dom,1,("sides id %d type %d left %d right %d\n",
                      PATCH_ID(thePatch),PATCH_TYPE(thePatch),
                      PARAM_PATCH_LEFT(thePatch),
                      PARAM_PATCH_RIGHT(thePatch)));
    for (i=0; i<2*DIM_OF_BND; i++)
      PRINTDEBUG(dom,1,("   corners %d",PARAM_PATCH_POINTS(thePatch,i)));
    PRINTDEBUG(dom,1,("\n"));
  }
  theBVP->numOfSubdomains = maxSubDomains;
  PRINTDEBUG(dom,1,(" bvp nsubcf %x\n",theBVP->numOfSubdomains));
  for (i=0; i<nsides; i++)
    if (sides[i] == NULL)
      return(NULL);

  if (theProblem != NULL)
    for (theBndCond=GetFirstBoundaryCondition(theProblem);
         theBndCond!=NULL; theBndCond = GetNextBoundaryCondition(theBndCond))
    {
      i = theBndCond->id;
      if ((i<0)||(i>=nsides))
        return(NULL);
      thePatch = sides[i];
      PARAM_PATCH_BC(thePatch) = theBndCond->BndCond;
      PARAM_PATCH_BCD(thePatch) = theBndCond->data;
    }

  /* create point patches */
  corners = (PATCH **)GetTmpMem(Heap,ncorners*sizeof(PATCH *));
  if (corners == NULL)
    return(NULL);
  theBVP->ncorners = ncorners;
  for (i=0; i<ncorners; i++)
  {
    m = 0;
    /* count parameter patchs */
    for (j=0; j<nsides; j++)
      for (n=0; n<2*DIM_OF_BND; n++)
        if (PARAM_PATCH_POINTS(sides[j],n) == i)
          m++;
    thePatch=
      (PATCH *)GetFreelistMemory(Heap,sizeof(POINT_PATCH)
                                 +(m-1)*sizeof(struct point_on_patch));
    if (thePatch == NULL)
      return (NULL);
    PATCH_TYPE(thePatch) = POINT_PATCH_TYPE;
    PATCH_ID(thePatch) = i;
    m = 0;
    for (j=0; j<nsides; j++)
      for (n=0; n<2*DIM_OF_BND; n++)
        if (PARAM_PATCH_POINTS(sides[j],n) == i)
        {
          POINT_PATCH_PID(thePatch,m) = j;
          POINT_PATCH_CID(thePatch,m++) = n;
        }
    POINT_PATCH_N(thePatch) = m;
    corners[i] = thePatch;
    PRINTDEBUG(dom,1,("corners id %d type %d n %d\n",
                      PATCH_ID(thePatch),PATCH_TYPE(thePatch),
                      POINT_PATCH_N(thePatch)));
    for (j=0; j<POINT_PATCH_N(thePatch); j++)
      PRINTDEBUG(dom,1,(" pid %d cid %d\n",
                        POINT_PATCH_PID(thePatch,j),
                        POINT_PATCH_CID(thePatch,j)));
  }

  /* create line patches */
  nlines = 0;
    #ifdef __THREEDIM__
  lines = (PATCH **)GetTmpMem(Heap,ncorners*ncorners*sizeof(PATCH *));
  if (lines == NULL)
    return(NULL);
  for (i=0; i<ncorners; i++)
    for (j=i+1; j<ncorners; j++)
    {
      k = 0;
      for (n=0; n<POINT_PATCH_N(corners[i]); n++)
        for (m=0; m<POINT_PATCH_N(corners[j]); m++)
          if (POINT_PATCH_PID(corners[i],n) ==
              POINT_PATCH_PID(corners[j],m))
            k++;
      if (k < 2)
        continue;
      thePatch=
        (PATCH *)GetFreelistMemory(Heap,sizeof(LINE_PATCH)
                                   +(k-1)*sizeof(struct line_on_patch));
      if (thePatch == NULL)
        return (NULL);
      PATCH_TYPE(thePatch) = LINE_PATCH_TYPE;
      PATCH_ID(thePatch) = nlines;
      k = 0;
      for (n=0; n<POINT_PATCH_N(corners[i]); n++)
        for (m=0; m<POINT_PATCH_N(corners[j]); m++)
          if (POINT_PATCH_PID(corners[i],n) ==
              POINT_PATCH_PID(corners[j],m))
          {
            LINE_PATCH_PID(thePatch,k) = POINT_PATCH_PID(corners[i],n);
            LINE_PATCH_CID0(thePatch,k) = POINT_PATCH_CID(corners[i],n);
            LINE_PATCH_CID1(thePatch,k) = POINT_PATCH_CID(corners[j],m);
            k++;
          }
      LINE_PATCH_N(thePatch) = k;
      lines[nlines++] = thePatch;
      PRINTDEBUG(dom,1,("lines id %d type %d n %d\n",
                        PATCH_ID(thePatch),PATCH_TYPE(thePatch),
                        LINE_PATCH_N(thePatch)));
      for (n=0; n<LINE_PATCH_N(thePatch); n++)
        PRINTDEBUG(dom,1,(" pid %d cid %d %d",
                          LINE_PATCH_PID(thePatch,n),
                          LINE_PATCH_CID0(thePatch,n),
                          LINE_PATCH_CID1(thePatch,n)));
      PRINTDEBUG(dom,1,("\n"));
    }
        #endif

  m = ncorners + nlines;
  theBVP->sideoffset = m;
  n = m + nsides;
  theBVP->patches=(PATCH **)GetFreelistMemory(Heap,n*sizeof(PATCH *));
  n = 0;
  for (i=0; i<ncorners; i++)
  {
    thePatch = corners[i];
    for (j=0; j<POINT_PATCH_N(thePatch); j++)
      POINT_PATCH_PID(thePatch,j) += m;
    theBVP->patches[n++] = thePatch;
  }
    #ifdef __THREEDIM__
  for (i=0; i<nlines; i++)
  {
    thePatch = lines[i];
    PATCH_ID(thePatch) = n;
    for (j=0; j<LINE_PATCH_N(thePatch); j++)
      LINE_PATCH_PID(thePatch,j) += m;
    theBVP->patches[n++] = thePatch;
  }
        #endif
  for (i=0; i<nsides; i++)
  {
    thePatch = sides[i];
    PATCH_ID(thePatch) = n;
    theBVP->patches[n++] = thePatch;
  }

  PRINTDEBUG(dom,1,("ncorners %d\n",theBVP->ncorners));
  for (i=0; i<theBVP->ncorners; i++)
    PRINTDEBUG(dom,1,("   id %d\n",PATCH_ID(theBVP->patches[i])));

  Mesh->nBndP = theBVP->ncorners;
  Mesh->nInnP = 0;
  Mesh->nElements = NULL;
  Mesh->theBndPs = (BNDP **) GetTmpMem(Heap,n*sizeof(BNDP *));
  if (Mesh->theBndPs == NULL)
    return(NULL);

  if (CreateCornerPoints(Heap,theBVP,Mesh->theBndPs))
    return(NULL);

  PRINTDEBUG(dom,1,("mesh n %d\n",Mesh->nBndP));
  for (i=0; i<theBVP->ncorners; i++)
    PRINTDEBUG(dom,1,(" id %d\n",
                      BND_PATCH_ID((BND_PS*)(Mesh->theBndPs[i]))));

  currBVP = theBVP;
  return ((BVP*)theBVP);
}

INT BVP_Dispose (BVP *theBVP)
{
  return (0);
}

/****************************************************************************/
/*D
   BVP_GetFirst - Return a pointer to the first STD_BVP

   SYNOPSIS:
   BVP *BVP_GetFirst (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function returns a pointer to the first BVP in the /STD_BVP
   directory.

   RETURN VALUE:
   BVP *
   .n   pointer to BVP
   .n   NULL if not found.
   D*/
/****************************************************************************/

BVP *BVP_GetFirst (void)
{
  ENVDIR *theSBVPDir;
  BVP *theBVP;

  theSBVPDir = ChangeEnvDir("/STD_BVP");
  assert (theSBVPDir!=NULL);
  theBVP = (BVP *) ENVDIR_DOWN(theSBVPDir);

  return (theBVP);
}

/****************************************************************************/
/*D
   BVP_GetNext - Return a pointer to the next multigrid

   SYNOPSIS:
   BVP *BVP_GetNext (BVP *theBVP);

   PARAMETERS:
   .  theBVP - BVP structure

   DESCRIPTION:
   This function returns a pointer to the next BVP in the /STD_BVP
   directory.

   RETURN VALUE:
   BVP *
   .n   pointer to BVP
   .n   NULL if not found.
   D*/
/****************************************************************************/

BVP *BVP_GetNext (BVP *theBVP)
{
  if (theBVP==NULL) return (NULL);
  return ((BVP *) NEXT_ENVITEM(theBVP));
}

/****************************************************************************/
/*D
   BVP_GetByName - get pointer to BVP by name

   SYNOPSIS:
   BVP *BVP_GetByName (char *name);

   PARAMETERS:
   .  name - name of file
   .  argc, argv - command parameters

   DESCRIPTION:
   This function gives the pointer to the BVP by its <name>.

   RETURN VALUE:
   BVP *
   .n   pointer to BVP
   .n   NULL if error.
   D*/
/****************************************************************************/

BVP *BVP_GetByName (char *name)
{
  return((BVP *) SearchEnv(name,"/BVP",theBVPDirID,theBVPDirID));
}

/****************************************************************************/
/*D
   BVP_SetBVPDesc - set BVP-descriptor

   SYNOPSIS:
   INT BVP_SetBVPDesc (BVP *theBVP, BVP_DESC *theBVPDesc);

   PARAMETERS:
   .  theBVP - BVP structure
   .  theBVPDesc - descriptor to set

   DESCRIPTION:
   This function sets the BVP descriptor according to the BVP.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BVP_SetBVPDesc (BVP *aBVP, BVP_DESC *theBVPDesc)
{
  STD_BVP *theBVP;
  INT i;

  /* cast */
  theBVP = GetSTD_BVP(aBVP);

  /* general part */
  strcpy(theBVPDesc->name,ENVITEM_NAME(theBVP));

  /* the domain part */
  for (i=0; i<DIM; i++)
    theBVPDesc->midpoint[i] = theBVP->MidPoint[i];
  theBVPDesc->radius              = theBVP->radius;
  theBVPDesc->convex              = theBVP->domConvex;
  theBVPDesc->nSubDomains         = theBVP->numOfSubdomains;
  theBVPDesc->numOfCoeffFct       = theBVP->numOfCoeffFct;
  theBVPDesc->numOfUserFct        = theBVP->numOfUserFct;
  theBVPDesc->ConfigProc          = theBVP->ConfigProc;

  currBVP = theBVP;

  return (0);
}

/****************************************************************************/
/*D
   BVP_SetCoeffFct - set coefficient function(s)

   SYNOPSIS:
   INT BVP_SetCoeffFct (BVP *theBVP, INT n, CoeffProcPtr *CoeffFct);

   PARAMETERS:
   .  theBVP - BVP structure
   .  n - nb. of coefficient function or -1 for all

   DESCRIPTION:
   This function one or all coefficient functions.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BVP_SetCoeffFct (BVP *aBVP, INT n, CoeffProcPtr *CoeffFct)
{
  STD_BVP *theBVP;
  INT i;

  theBVP = GetSTD_BVP(aBVP);

  /* check */
  if (n<-1 || n>=theBVP->numOfCoeffFct) return (1);

  if (n==-1)
    for (i=0; i<theBVP->numOfCoeffFct; i++)
      CoeffFct[i] = (CoeffProcPtr)theBVP->CU_ProcPtr[i];
  else
    CoeffFct[0] = (CoeffProcPtr)theBVP->CU_ProcPtr[n];

  return (0);
}

/****************************************************************************/
/*D
   BVP_SetUserFct - set coefficient function(s)

   SYNOPSIS:
   INT BVP_SetUserFct (BVP *theBVP, INT n, UserProcPtr *UserFct);

   PARAMETERS:
   .  theBVP - BVP structure
   .  n - nb. of user function or -1 for all

   DESCRIPTION:
   This function gives one or all user functions.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BVP_SetUserFct (BVP *aBVP, INT n, UserProcPtr *UserFct)
{
  STD_BVP *theBVP;
  INT i;

  theBVP = GetSTD_BVP(aBVP);

  /* check */
  if (n<-1 || n>=theBVP->numOfUserFct) return (1);

  if (n==-1)
    for (i=0; i<theBVP->numOfUserFct; i++)
      UserFct[i] = (UserProcPtr)theBVP->CU_ProcPtr[i+theBVP->numOfCoeffFct];
  else
    UserFct[0] = (UserProcPtr)theBVP->CU_ProcPtr[n+theBVP->numOfCoeffFct];

  return (0);
}

/****************************************************************************/
/*D
   BVP_InsertBndP - sets a BNDP from command input

   SYNOPSIS:
   BNDP *BVP_InsertBndP (HEAP *Heap, BVP *theBVP, INT argc, char **argv);

   PARAMETERS:
   .  theBVP - BVP structure
   .  argc, argv - command parameters
   .  theBndP - the BNDP to set

   DESCRIPTION:
   This function sets a BNDP from command input parameters.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

static INT GetNumberOfPatches(PATCH *p)
{
  switch (PATCH_TYPE(p))
  {
  case PARAMETRIC_PATCH_TYPE :
    return(1);
  case POINT_PATCH_TYPE :
    return(POINT_PATCH_N(p));
      #ifdef __THREEDIM__
  case LINE_PATCH_TYPE :
    return(LINE_PATCH_N(p));
      #endif
  }

  return(-1);
}

static INT GetPatchId (PATCH *p, INT i)
{
  switch (PATCH_TYPE(p))
  {
  case PARAMETRIC_PATCH_TYPE :
    return(PATCH_ID(p));
  case POINT_PATCH_TYPE :
    return(POINT_PATCH_PID(p,i));
      #ifdef __THREEDIM__
  case LINE_PATCH_TYPE :
    return(LINE_PATCH_PID(p,i));
      #endif
  }

  return(-1);
}

static INT GetNumberOfCommonPatches (PATCH *p0, PATCH *p1)
{
  INT i,j,cnt;

  cnt = 0;
  for (i=0; i<GetNumberOfPatches(p0); i++)
    for (j=0; j<GetNumberOfPatches(p1); j++)
      if (GetPatchId(p0,i) == GetPatchId(p1,j))
        cnt++;

  return(cnt);
}

static INT GetCommonPatchId (PATCH *p0, PATCH *p1, INT k)
{
  INT i,j,cnt;

  cnt = 0;
  for (i=0; i<GetNumberOfPatches(p0); i++)
    for (j=0; j<GetNumberOfPatches(p1); j++)
      if (GetPatchId(p0,i) == GetPatchId(p1,j))
        if (k == cnt++)
          return(GetPatchId(p1,j));

  return(-1);
}

#ifdef __THREEDIM__
static INT GetCommonLinePatchId (PATCH *p0, PATCH *p1)
{
  INT i,j,k,l,cnt,cnt1;
  PATCH *p;

  if (PATCH_TYPE(p0) == LINE_PATCH_TYPE)
    return(PATCH_ID(p0));
  else if (PATCH_TYPE(p1) == LINE_PATCH_TYPE)
    return(PATCH_ID(p1));

  cnt = GetNumberOfCommonPatches(p0,p1);

  if (cnt < 1)
    return(-1);

  for (k=currBVP->ncorners; k<currBVP->sideoffset; k++)
  {
    p = currBVP->patches[k];
    if (LINE_PATCH_N(p) != cnt)
      continue;
    cnt1 = 0;
    for (i=0; i<cnt; i++)
      for (l=0; l<LINE_PATCH_N(p); l++)
        if (GetCommonPatchId(p0,p1,i) == LINE_PATCH_PID(p,l))
          cnt1++;
    if (cnt == cnt1)
      return(k);
  }

  return(-1);
}

BNDP *CreateBndPOnLine (HEAP *Heap, PATCH *p0, PATCH *p1, COORD lcoord)
{
  BND_PS *bp;
  PATCH *p,*pp;
  COORD local0[DIM_OF_BND];
  COORD local1[DIM_OF_BND];
  INT k,l,cnt;

  if (PATCH_TYPE(p0) != POINT_PATCH_TYPE)
    return(NULL);
  if (PATCH_TYPE(p1) != POINT_PATCH_TYPE)
    return(NULL);

  PRINTDEBUG(dom,1,("    p0 p1 %d %d\n",PATCH_ID(p0),PATCH_ID(p1)));
  for (l=0; l<GetNumberOfPatches(p0); l++)
    PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p0,l)));
  for (l=0; l<GetNumberOfPatches(p1); l++)
    PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p1,l)));

  cnt = GetNumberOfCommonPatches(p0,p1);

  if (cnt < 2)
    return(NULL);

  bp = (BND_PS *)GetFreelistMemory(Heap,(cnt-1)*sizeof(COORD_BND_VECTOR)
                                   + sizeof(BND_PS));
  if (bp == NULL)
    return(NULL);
  bp->n = cnt;

  k = GetCommonLinePatchId (p0,p1);
  if ((k<currBVP->ncorners) || (k>=currBVP->sideoffset))
    return(NULL);
  p = currBVP->patches[k];
  bp->patch_id = k;

  PRINTDEBUG(dom,1,(" Create BNDP line %d cnt %d\n",k,cnt));
  for (l=0; l<GetNumberOfPatches(p0); l++)
    PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p0,l)));
  for (l=0; l<GetNumberOfPatches(p1); l++)
    PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p1,l)));

  for (l=0; l<LINE_PATCH_N(p); l++)
  {
    pp = currBVP->patches[LINE_PATCH_PID(p,l)];
    switch (LINE_PATCH_CID0(p,l))
    {
    case 0 :
      local0[0] = PARAM_PATCH_RANGE(pp)[0][0];
      local0[1] = PARAM_PATCH_RANGE(pp)[0][1];
      break;
    case 1 :
      local0[0] = PARAM_PATCH_RANGE(pp)[1][0];
      local0[1] = PARAM_PATCH_RANGE(pp)[0][1];
      break;
    case 2 :
      local0[0] = PARAM_PATCH_RANGE(pp)[1][0];
      local0[1] = PARAM_PATCH_RANGE(pp)[1][1];
      break;
    case 3 :
      local0[0] = PARAM_PATCH_RANGE(pp)[0][0];
      local0[1] = PARAM_PATCH_RANGE(pp)[1][1];
      break;
    }
    switch (LINE_PATCH_CID1(p,l))
    {
    case 0 :
      local1[0] = PARAM_PATCH_RANGE(pp)[0][0];
      local1[1] = PARAM_PATCH_RANGE(pp)[0][1];
      break;
    case 1 :
      local1[0] = PARAM_PATCH_RANGE(pp)[1][0];
      local1[1] = PARAM_PATCH_RANGE(pp)[0][1];
      break;
    case 2 :
      local1[0] = PARAM_PATCH_RANGE(pp)[1][0];
      local1[1] = PARAM_PATCH_RANGE(pp)[1][1];
      break;
    case 3 :
      local1[0] = PARAM_PATCH_RANGE(pp)[0][0];
      local1[1] = PARAM_PATCH_RANGE(pp)[1][1];
      break;
    }
    bp->local[l][0] = (1.0-lcoord)*local0[0]+lcoord*local1[0];
    bp->local[l][1] = (1.0-lcoord)*local0[1]+lcoord*local1[1];

    PRINTDEBUG(dom,1,(" Create bndp %d line %d l %d %f %f\n",
                      bp->patch_id,
                      LINE_PATCH_PID(p,l),l,
                      bp->local[l][0],bp->local[l][1]));
  }

  return((BNDP *)bp);
}
#endif

BNDP *BVP_InsertBndP (HEAP *Heap, BVP *aBVP, INT argc, char **argv)
{
  STD_BVP *theBVP;
  BND_PS *ps;
  PATCH *p;
  COORD lc;
  INT j,n,pid;
  int i;
  float pos[DIM_OF_BND];

  theBVP = GetSTD_BVP(aBVP);

    #ifdef __TWODIM__
  if (sscanf(argv[0],"bn %d %f",&i,pos) != 2)
    return(NULL);
    #endif
    #ifdef __THREEDIM__
  if (sscanf(argv[0],"bn %d %f %f",&i,pos,pos+1) != 3)
    return(NULL);
    #endif

  pid = i + theBVP->sideoffset;
  p = theBVP->patches[pid];

    #ifdef __THREEDIM__
  if (ABS(pos[0] - PARAM_PATCH_RANGE(p)[0][0]) < SMALL_DIFF)
  {
    lc =  (pos[1] - PARAM_PATCH_RANGE(p)[0][1])
         / (PARAM_PATCH_RANGE(p)[1][1] - PARAM_PATCH_RANGE(p)[0][1]);
    return(CreateBndPOnLine(Heap,
                            currBVP->patches[PARAM_PATCH_POINTS(p,0)],
                            currBVP->patches[PARAM_PATCH_POINTS(p,3)],lc));
  }
  else if (ABS(pos[0] - PARAM_PATCH_RANGE(p)[1][0]) < SMALL_DIFF)
  {
    lc =  (pos[1] - PARAM_PATCH_RANGE(p)[0][1])
         / (PARAM_PATCH_RANGE(p)[1][1] - PARAM_PATCH_RANGE(p)[0][1]);
    return(CreateBndPOnLine(Heap,
                            currBVP->patches[PARAM_PATCH_POINTS(p,1)],
                            currBVP->patches[PARAM_PATCH_POINTS(p,2)],lc));
  }
  else if (ABS(pos[1] - PARAM_PATCH_RANGE(p)[0][1]) < SMALL_DIFF)
  {
    lc =  (pos[0] - PARAM_PATCH_RANGE(p)[0][0])
         / (PARAM_PATCH_RANGE(p)[1][0] - PARAM_PATCH_RANGE(p)[0][0]);
    return(CreateBndPOnLine(Heap,
                            currBVP->patches[PARAM_PATCH_POINTS(p,0)],
                            currBVP->patches[PARAM_PATCH_POINTS(p,1)],lc));
  }
  else if (ABS(pos[1] - PARAM_PATCH_RANGE(p)[1][1]) < SMALL_DIFF)
  {
    lc =  (pos[0] - PARAM_PATCH_RANGE(p)[0][0])
         / (PARAM_PATCH_RANGE(p)[1][0] - PARAM_PATCH_RANGE(p)[0][0]);
    return(CreateBndPOnLine(Heap,
                            currBVP->patches[PARAM_PATCH_POINTS(p,3)],
                            currBVP->patches[PARAM_PATCH_POINTS(p,2)],lc));
  }
    #endif

  if (PATCH_TYPE(p) == PARAMETRIC_PATCH_TYPE)
  {
    PRINTDEBUG(dom,1,(" id %d i %d ns %d %s pos %f\n",
                      pid,i,currBVP->sideoffset,argv[0],pos[0]));
    ps = (BND_PS *)GetFreelistMemory(Heap,sizeof(BND_PS));
    if (ps == NULL)
      return(NULL);
    ps->patch_id = pid;
    ps->n = 1;
    for (j=0; j<DIM_OF_BND; j++)
      ps->local[0][j] = pos[j];
  }
  else
    return(NULL);

  return ((BNDP *)ps);
}

static DOUBLE LengthOfSide (PATCH *p)
{
  DOUBLE length,step,s;
  COORD lambda[DIM_OF_BND],x[2][DIM];
  INT i;

  length = 0.0;

  lambda[0] = PARAM_PATCH_RANGE(p)[0][0];
  step  = (PARAM_PATCH_RANGE(p)[1][0]
           - PARAM_PATCH_RANGE(p)[0][0]) / RESOLUTION;
  if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),lambda,x[0]))
    return(0.0);
  for (i=1; i<RESOLUTION; i++)
  {
    lambda[0] += step;
    if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),lambda,x[i%2]))
      return(0.0);
    V_DIM_EUKLIDNORM_OF_DIFF(x[0],x[1],s);
    length += s;
  }
  lambda[0] = PARAM_PATCH_RANGE(p)[1][0];
  if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),lambda,x[RESOLUTION%2]))
    return(0.0);
  V_DIM_EUKLIDNORM_OF_DIFF(x[0],x[1],s);
  length += s;

  return(length);
}

static DOUBLE MeshSize (CoeffProcPtr coeff, PATCH *p, COORD *lambda)
{
  DOUBLE step;
  COORD global[DIM];

  if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),lambda,global))
    return(0.0);
  (*coeff)(global,&step);

  PRINTDEBUG(dom,1,(" c:lambda %f x %f %f step %f\n",
                    lambda[0],global[0],global[1],step));


  return(step);
}

INT GenerateBnodes (HEAP *Heap, STD_BVP *theBVP, BNDP **bndp,
                    INT *sides, INT ***corners, CoeffProcPtr coeff)
{
  INT i,j,n,nside,left,right;
  DOUBLE length,plength,step,step1;
  COORD lambda[DIM_OF_BND],lambda1,x[DIM];
  PATCH *p;
  BND_PS *ps;

  n = theBVP->ncorners;

        #ifdef __THREEDIM__
  return(n);
    #endif

  for (i=0; i<=theBVP->numOfSubdomains; i++)
    sides[i] = 0;

  for (j=0; j<=theBVP->numOfSubdomains; j++)
    PRINTDEBUG(dom,1,("  g h       nside %d %d %d\n",nside,
                      theBVP->nsides,sides[j]));

  for (i=theBVP->sideoffset; i<theBVP->sideoffset+theBVP->nsides; i++)
  {
    nside = n;
    p = theBVP->patches[i];
    length = LengthOfSide(p);
    if (length == 0.0)
      return(-1);
    plength = ABS(PARAM_PATCH_RANGE(p)[1][0]
                  - PARAM_PATCH_RANGE(p)[0][0]);
    lambda[0] = MIN(PARAM_PATCH_RANGE(p)[0][0],
                    PARAM_PATCH_RANGE(p)[1][0]);
    lambda1 = MAX(PARAM_PATCH_RANGE(p)[0][0],
                  PARAM_PATCH_RANGE(p)[1][0]);

    PRINTDEBUG(dom,1,(" side i %d lenght %f pl %f\n",i,length,plength));

    step = MeshSize (coeff,p,lambda) * plength / length;
    if (step == 0.0)
      return(-1);
    lambda[0] += step;
    while (lambda[0] < lambda1)
    {
      if (lambda[0]+step > lambda1)
        lambda[0] = 0.5 * (lambda1 + lambda[0] - step);
      else
      {
        step1 = MeshSize (coeff,p,lambda) * plength / length;
        if (step1 == 0.0)
          return(-1);
        if (step1 < step)
        {
          lambda[0] += step1 - step;
          step = MeshSize (coeff,p,lambda) * plength / length;
        }
        else
          step = step1;
      }
      if (bndp != NULL)
      {
        ps = (BND_PS *)GetFreelistMemory(Heap,sizeof(BND_PS));
        if (ps == NULL)
          return(0);
        ps->n = 1;
        bndp[n] = (BNDP *)ps;
        ps->patch_id = i;
        ps->local[0][0] = lambda[0];
      }
      n++;
      lambda[0] += step;

      PRINTDEBUG(dom,1,("  lam %f sp %f\n",lambda[0],step));

    }
    if (n == nside)
    {
      if (bndp != NULL)
      {
        ps = (BND_PS *)bndp[n];
        ps->patch_id = i;
        ps->local[0][0] = 0.5 *(PARAM_PATCH_RANGE(p)[0][0]+
                                PARAM_PATCH_RANGE(p)[1][0]);
      }
      n++;
    }
    left = PARAM_PATCH_LEFT(p);
    right = PARAM_PATCH_RIGHT(p);
    if (corners != NULL)
    {
      if (left > 0)
      {
        if (PARAM_PATCH_RANGE(p)[1][0] > PARAM_PATCH_RANGE(p)[0][0])
        {
          corners[left][sides[left]][0] = PARAM_PATCH_POINTS(p,0);
          for (j=nside; j<n; j++)
          {
            corners[left][sides[left]++][1] = j;
            corners[left][sides[left]][0] = j;
          }
          corners[left][sides[left]++][1] = PARAM_PATCH_POINTS(p,1);
        }
        else
        {
          corners[left][sides[left]][1] = PARAM_PATCH_POINTS(p,0);
          for (j=n-1; j>=nside; j--)
          {
            corners[left][sides[left]++][0] = j;
            corners[left][sides[left]][1] = j;
          }
          corners[left][sides[left]++][0] = PARAM_PATCH_POINTS(p,1);
        }
      }
      if (right > 0)
      {
        if (PARAM_PATCH_RANGE(p)[1][0] > PARAM_PATCH_RANGE(p)[0][0])
        {
          corners[right][sides[right]][1] = PARAM_PATCH_POINTS(p,0);
          for (j=nside; j<n; j++)
          {
            corners[right][sides[right]++][0] = j;
            corners[right][sides[right]][1] = j;
          }
          corners[right][sides[right]++][0] = PARAM_PATCH_POINTS(p,1);
        }
        else
        {
          corners[right][sides[right]][0] = PARAM_PATCH_POINTS(p,0);
          for (j=n-1; j>=nside; j--)
          {
            corners[right][sides[right]++][1] = j;
            corners[right][sides[right]][0] = j;
          }
          corners[right][sides[right]++][1] = PARAM_PATCH_POINTS(p,1);
        }
      }
    }
    else
    {
      if (left > 0)
        sides[left] += n-nside+1;
      if (right > 0)
        sides[right] += n-nside+1;
    }
  }

  return (n);
}

#ifdef __TWODIM__
INT GenerateBnodes_h (HEAP *Heap, STD_BVP *theBVP, BNDP **bndp,
                      INT *sides, INT ***corners, DOUBLE h)
{
  INT i,j,m,n,nside,left,right;
  DOUBLE length,plength,step,step1;
  COORD lambda[DIM_OF_BND],lambda1,x[DIM];
  PATCH *p;
  BND_PS *ps;

  for (j=0; j<=theBVP->numOfSubdomains; j++)
    PRINTDEBUG(dom,1,("  g h       nside %d %d %d\n",nside,n,sides[j]));

  n = theBVP->ncorners;
  for (i=0; i<=theBVP->numOfSubdomains; i++)
    sides[i] = 0;

  for (i=theBVP->sideoffset; i<theBVP->sideoffset+theBVP->nsides; i++)
  {
    nside = n;
    p = theBVP->patches[i];
    length = LengthOfSide(p);
    if (length == 0.0)
      return(-1);
    plength = ABS(PARAM_PATCH_RANGE(p)[1][0]
                  - PARAM_PATCH_RANGE(p)[0][0]);
    lambda[0] = MIN(PARAM_PATCH_RANGE(p)[0][0],
                    PARAM_PATCH_RANGE(p)[1][0]);
    lambda1 = MAX(PARAM_PATCH_RANGE(p)[0][0],
                  PARAM_PATCH_RANGE(p)[1][0]);


    PRINTDEBUG(dom,1,(" h:side n %d lenght %f pl %f\n",i,length,plength));

    m = length / h;
    if (m < 2)
      m = 2;
    step  = plength / m;
    for (j=1; j<m; j++)
    {
      lambda[0] += step;
      if (bndp != NULL)
      {
        ps = (BND_PS *)GetFreelistMemory(Heap,sizeof(BND_PS));
        if (ps == NULL)
          return(0);
        ps->n = 1;
        bndp[n] = (BNDP *)ps;
        ps->patch_id = i;
        ps->local[0][0] = lambda[0];
      }
      n++;
    }
    left = PARAM_PATCH_LEFT(p);
    right = PARAM_PATCH_RIGHT(p);
    if (corners != NULL)
    {
      if (left > 0)
      {
        if (PARAM_PATCH_RANGE(p)[1][0] > PARAM_PATCH_RANGE(p)[0][0])
        {
          corners[left][sides[left]][0] = PARAM_PATCH_POINTS(p,0);
          for (j=nside; j<n; j++)
          {
            corners[left][sides[left]++][1] = j;
            corners[left][sides[left]][0] = j;
          }
          corners[left][sides[left]++][1] = PARAM_PATCH_POINTS(p,1);
        }
        else
        {
          corners[left][sides[left]][1] = PARAM_PATCH_POINTS(p,0);
          for (j=n-1; j>=nside; j--)
          {
            corners[left][sides[left]++][0] = j;
            corners[left][sides[left]][1] = j;
          }
          corners[left][sides[left]++][0] = PARAM_PATCH_POINTS(p,1);
        }
      }
      if (right > 0)
      {
        if (PARAM_PATCH_RANGE(p)[1][0] > PARAM_PATCH_RANGE(p)[0][0])
        {
          corners[right][sides[right]][1] = PARAM_PATCH_POINTS(p,0);
          for (j=nside; j<n; j++)
          {
            corners[right][sides[right]++][0] = j;
            corners[right][sides[right]][1] = j;
          }
          corners[right][sides[right]++][0] = PARAM_PATCH_POINTS(p,1);
        }
        else
        {
          corners[right][sides[right]][0] = PARAM_PATCH_POINTS(p,0);
          for (j=n-1; j>=nside; j--)
          {
            corners[right][sides[right]++][1] = j;
            corners[right][sides[right]][0] = j;
          }
          corners[right][sides[right]++][1] = PARAM_PATCH_POINTS(p,1);
        }
      }
    }
    else
    {
      if (left > 0)
        sides[left] += n-nside+1;
      if (right > 0)
        sides[right] += n-nside+1;
    }
  }

  return (n);
}
#endif

/****************************************************************************/
/*D
   AddBoundaryElements - decompose a strip into triangles

   SYNOPSIS:
   static INT AddBoundaryElements (INT n, INT m,
   INT c0, INT c1, INT c2, INT c3,
   INT s0, INT s1, INT s2, INT s3);

   PARAMETERS:
   .  n,m - stripe with n+1 nodes on the bottom and m+1 nodes on the top
   .  c0,c1,c2,c3 - corner node ids
   .  s0,s1,s2,s3 - side node ids

   DESCRIPTION:
   This function splits a stripe into triangles an calls 'AddBoundaryElements'.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

static INT AddBoundaryElement (INT n, INT *nodelist,
                               INT left, INT right,
                               INT ***corners, INT *sides)
{
  if (left > 0)
  {
    if (corners != NULL)
    {
      corners[left][sides[left]][0] = nodelist[0];
      corners[left][sides[left]][1] = nodelist[1];
      corners[left][sides[left]][2] = nodelist[2];
    }
    sides[left]++;
  }
  if (right > 0)
  {
    if (corners != NULL)
    {
      corners[right][sides[right]][0] = nodelist[0];
      corners[right][sides[right]][1] = nodelist[2];
      corners[right][sides[right]][2] = nodelist[1];
    }
    sides[right]++;
  }

  return(0);
}

static INT AddBoundaryElements (INT n, INT m,
                                INT c0, INT c1, INT c2, INT c3,
                                INT s0, INT s1, INT s2, INT s3,
                                INT left, INT right,
                                INT ***corners, INT *sides)
{
  INT nodelist[3];

  PRINTDEBUG(dom,1,("    add n %d m %d c0 %d c1 %d s0 %d s1 %d\n",
                    n,m,c0,c1,s0,s1));
  if (m < n)
  {
    if (n == 1)
    {
      nodelist[0] = c0;
      nodelist[1] = c2;
      nodelist[2] = c1;
      AddBoundaryElement(3,nodelist,left,right,corners,sides);
    }
    else
    {
      nodelist[0] = c0;
      nodelist[1] = c2;
      nodelist[2] = s0;
      AddBoundaryElement(3,nodelist,left,right,corners,sides);
      c0 = s0;
      if (s0<s1) s0++;else s0--;
      AddBoundaryElements (n-1,m,c0,c1,c2,c3,s0,s1,s2,s3,
                           left,right,corners,sides);
    }
  }
  else
  {
    if (m == 1)
    {
      nodelist[0] = c0;
      nodelist[1] = c2;
      nodelist[2] = c1;
      AddBoundaryElement(3,nodelist,left,right,corners,sides);
      nodelist[0] = c1;
      nodelist[1] = c2;
      nodelist[2] = c3;
      AddBoundaryElement(3,nodelist,left,right,corners,sides);
    }
    else
    {
      nodelist[0] = c2;
      nodelist[1] = s2;
      nodelist[2] = c0;
      AddBoundaryElement(3,nodelist,left,right,corners,sides);
      c2 = s2;
      if (s2<s3) s2++;else s2--;
      AddBoundaryElements (n,m-1,c0,c1,c2,c3,s0,s1,s2,s3,
                           left,right,corners,sides);
    }
  }

  return(0);
}

/****************************************************************************/
/*D
   TriangulatePatch - decompose a patch into stripes

   SYNOPSIS:
   static INT TriangulatePatch (DOUBLE h, PATCH *thePatch,
   INT npc, INT *cornerid, COORD local[CORNERS_OF_BND_SEG][DIM-1],
   INT sideid[CORNERS_OF_BND_SEG][2], INT *siden);

   PARAMETERS:
   .  h - maximal length of an edge of the boundary triangles
   .  thePatch - poiter to a patch
   .  npc - number of corner nodes of the patch
   .  cornerid - ids of the corner nodes
   .  local - local coordinates of the corners
   .  sideid - ids of the nodes on the edges of the patch
   .  siden - number of nodes on the edges of the patch

   DESCRIPTION:
   This function splits the patch into stripes and calls
   'AddBoundaryElements' to decompose the stripes into triangles.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

static INT nc,nodeid;

static INT TriangulatePatch (HEAP *Heap, PATCH *p, BNDP **bndp,
                             INT *sides, INT ***corners, DOUBLE h,
                             INT npc, INT *cornerid,
                             COORD local[CORNERS_OF_BND_SEG][DIM-1],
                             INT sideid[CORNERS_OF_BND_SEG][2],
                             INT *siden)
{
  BND_PS *ps;
  INT i,j,k,n,left,right,nside;
  COORD lvect[DIM-1],gvect[DIM],gvect1[DIM];
  COORD next_local[CORNERS_OF_BND_SEG][DIM-1];
  DOUBLE lambda,dist,step;
  INT next_cornerid[CORNERS_OF_BND_SEG];
  INT next_sideid[CORNERS_OF_BND_SEG][2];
  INT next_siden[CORNERS_OF_BND_SEG];

  left = PARAM_PATCH_LEFT(p);
  right = PARAM_PATCH_RIGHT(p);


  PRINTDEBUG(dom,1,("Triang  nid %d sn %d %d %d %d\n",nodeid,
                    siden[0],siden[1],siden[2],siden[3]));

  if ((siden[npc-1] > 1) && (siden[1] > 1))
  {
    for (k=2; k<npc; k++)
      for (i=0; i<DIM-1; i++)
        next_local[k][i] = local[k][i];
    next_siden[0] = n;
    next_siden[1] = siden[1] - 1;
    next_siden[2] = siden[2];
    next_siden[npc-1] = siden[npc-1] - 1;
    next_sideid[2][0] = sideid[2][0];
    next_sideid[2][1] = sideid[2][1];
    for (k=2; k<npc; k++)
      next_cornerid[k] = cornerid[k];
    next_cornerid[0] = sideid[npc-1][1];
    next_sideid[npc-1][0] = sideid[npc-1][0];
    if (sideid[npc-1][0] < sideid[npc-1][1])
      next_sideid[npc-1][1] = sideid[npc-1][1] - 1;
    else
      next_sideid[npc-1][1] = sideid[npc-1][1] + 1;
    lambda = (siden[npc-1] - 1.0) / siden[npc-1];
    V2_LINCOMB(lambda,local[0],(1.0-lambda),local[npc-1],
               next_local[0]);
    next_cornerid[1] = sideid[1][0];
    next_sideid[1][1] = sideid[1][1];
    if (sideid[1][0] < sideid[1][1])
      next_sideid[1][0] = sideid[1][0] + 1;
    else
      next_sideid[1][0] = sideid[1][0] - 1;
    lambda = (siden[1] - 1.0) / siden[1];
    V2_LINCOMB(lambda,local[1],(1.0-lambda),local[2],next_local[1]);
    next_siden[0] = siden[0];
    if (h > 0)
    {
      if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),next_local[0],gvect))
        return(-1);
      if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),next_local[1],gvect1))
        return(-1);
      V_DIM_EUKLIDNORM_OF_DIFF(gvect,gvect1,dist);
      if (((INT)(dist/h)) < siden[0])
        next_siden[0] -= 1;
      else if (((INT)(dist/h)) > siden[0])
        next_siden[0] += 1;
    }
    next_sideid[0][0] = nodeid;
    next_sideid[0][1] = nodeid + next_siden[0] - 2;
    step = 1.0 / next_siden[0];
    if (bndp != NULL)
      for (i=1; i<next_siden[0]; i++)
      {
        ps = (BND_PS *)GetFreelistMemory(Heap,sizeof(BND_PS));
        if (ps == NULL)
          return(0);
        ps->n = 1;
        ps->patch_id = PATCH_ID(p);
        lambda = i * step;
        V2_LINCOMB(lambda,next_local[1],(1.0-lambda),next_local[0],
                   ps->local[0]);
        bndp[nodeid++] = (BNDP *)ps;


        PRINTDEBUG(dom,1,("    lambda nid %d %f %f\n",
                          nodeid-1,ps->local[0][0],ps->local[0][1]));

      }
    else
      nodeid += next_siden[0] -1;

    AddBoundaryElements(siden[0],next_siden[0],
                        cornerid[0],cornerid[1],
                        next_cornerid[0],next_cornerid[1],
                        sideid[0][0],sideid[0][1],
                        next_sideid[0][0],next_sideid[0][1],
                        left,right,corners,sides);

    return(TriangulatePatch(Heap,p,bndp,sides,corners,h,npc,next_cornerid,
                            next_local,next_sideid,next_siden));
  }
  else if ((siden[npc-1] == 1) && (siden[1] == 1))
  {
    if (npc == 3)
      return(AddBoundaryElements(siden[0],0,
                                 cornerid[0],cornerid[1],cornerid[2],0,
                                 sideid[0][0],sideid[0][1],0,0,
                                 left,right,corners,sides));
    else
      return(AddBoundaryElements(siden[0],siden[2],
                                 cornerid[0],cornerid[1],
                                 cornerid[3],cornerid[2],
                                 sideid[0][0],sideid[0][1],
                                 sideid[2][1],sideid[2][0],
                                 left,right,corners,sides));
  }
  else if (((siden[npc-1] > 1) && (siden[1] == 1)) ||
           ((siden[npc-1] == 1) && (siden[1] > 1))   )
  {
    if ((siden[npc-2] > 1) && (siden[1] == 1) && (siden[0] == 1))
    {
      UserWrite("TriangulatePatch: this case is not implemented\n");
      return(1);
    }

    /* swap sides */
    for (k=0; k<npc; k++)
    {
      for (i=0; i<DIM-1; i++)
        next_local[k][i] = local[(k+1)%npc][i];
      next_siden[k] = siden[(k+1)%npc];
      next_sideid[k][0] = sideid[(k+1)%npc][0];
      next_sideid[k][1] = sideid[(k+1)%npc][1];
      next_cornerid[k] = cornerid[(k+1)%npc];
    }
    return(TriangulatePatch(Heap,p,bndp,sides,corners,h,npc,next_cornerid,
                            next_local,next_sideid,next_siden));
  }

  return(0);
}

#define INDEX_IN_LIST(from,to,nc)  (from < to ? from*nc+to : to*nc+from)

#ifdef __THREEDIM__
static INT GenerateBnodes_h (HEAP *Heap, STD_BVP *theBVP, BNDP **bndp,
                             INT *sides, INT ***corners, DOUBLE h)
{
  INT i,j,m,n,nside,left,right,from,to,k;
  DOUBLE length,plength,step,step1;
  COORD lambda;
  COORD dist, global[CORNERS_OF_BND_SEG][DIM];
  COORD local[CORNERS_OF_BND_SEG][DIM-1];
  PATCH *p;
  BND_PS *ps;
  BNDP *bp;
  INT *vlist;
  INT cornerid[CORNERS_OF_BND_SEG],npc;
  INT sideid[CORNERS_OF_BND_SEG][2],siden[CORNERS_OF_BND_SEG];

  for (j=0; j<=theBVP->numOfSubdomains; j++)
    PRINTDEBUG(dom,1,("  g h       nside %d %d %d\n",nside,n,sides[j]));

  nc = theBVP->ncorners;
  nodeid = nc;
  if (bndp == NULL)
  {
    vlist = (INT*) GetTmpMem(Heap,nc*nc*sizeof(INT));
    if (vlist == NULL)
      return(nc);
  }

  for (i=0; i<nc; i++)
    for (j=0; j<nc; j++)
      vlist[INDEX_IN_LIST(i,j,nc)] = 0;

  for (i=0; i<=theBVP->numOfSubdomains; i++)
    sides[i] = 0;

  for (i=theBVP->sideoffset; i<theBVP->sideoffset+theBVP->nsides; i++)
  {
    p = theBVP->patches[i];
    npc = 4;
    local[0][0] = PARAM_PATCH_RANGE(p)[0][0];
    local[0][1] = PARAM_PATCH_RANGE(p)[0][1];
    local[1][0] = PARAM_PATCH_RANGE(p)[1][0];
    local[1][1] = PARAM_PATCH_RANGE(p)[0][1];
    local[2][0] = PARAM_PATCH_RANGE(p)[1][0];
    local[2][1] = PARAM_PATCH_RANGE(p)[1][1];
    local[3][0] = PARAM_PATCH_RANGE(p)[0][0];
    local[3][1] = PARAM_PATCH_RANGE(p)[1][1];
    for( k=0; k<npc; k++ )
      if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),local[k],global[k]))
        return(-1);

    for( k=0; k<npc; k++ )
      PRINTDEBUG(dom,1,(" gl %d %f %f %f\n",
                        PARAM_PATCH_POINTS(p,k),
                        global[k][0],
                        global[k][1],
                        global[k][2]));

    /* sorry for that */
    local[1][0] = PARAM_PATCH_RANGE(p)[0][0];
    local[1][1] = PARAM_PATCH_RANGE(p)[0][1];
    local[2][0] = PARAM_PATCH_RANGE(p)[1][0];
    local[2][1] = PARAM_PATCH_RANGE(p)[0][1];
    local[3][0] = PARAM_PATCH_RANGE(p)[1][0];
    local[3][1] = PARAM_PATCH_RANGE(p)[1][1];
    local[0][0] = PARAM_PATCH_RANGE(p)[0][0];
    local[0][1] = PARAM_PATCH_RANGE(p)[1][1];

    from = PARAM_PATCH_POINTS(p,npc-1);
    for( k=0; k<npc; k++ )
    {
      to = PARAM_PATCH_POINTS(p,k);
      cornerid[k] = from;
      V_DIM_EUKLIDNORM_OF_DIFF(global[(k+npc-1)%npc],global[k],dist);
      if (h > 0)
        n = MAX(1,1.00001*dist/h);
      else
        n = MAX(1,-h);
      siden[k] = n;
      sideid[k][0] = vlist[INDEX_IN_LIST(from,to,nc)];
      if (sideid[k][0] > 0)
      {
        if (from < to)
          sideid[k][1] = sideid[k][0] + n - 2;
        else
        {
          sideid[k][1] = sideid[k][0];
          sideid[k][0] = sideid[k][1] + n - 2;
        }
      }
      else
      {
        vlist[INDEX_IN_LIST(from,to,nc)] = nodeid;
        if (from < to)
        {
          sideid[k][0] = nodeid;
          sideid[k][1] = nodeid + n - 2;
        }
        else
        {
          sideid[k][1] = nodeid;
          sideid[k][0] = nodeid + n - 2;
        }
        if (bndp != NULL)
        {
          step = 1.0 / n;
          lambda = 0.0;
          for (j=1; j<n; j++)
          {
            lambda += step;
            bp = CreateBndPOnLine (Heap,theBVP->patches[from],
                                   theBVP->patches[to],lambda);
            if (bp == NULL)
              return(-1);
            bndp[nodeid++] = bp;
          }
        }
        else
          nodeid += (n-1);
      }

      PRINTDEBUG(dom,1,(" VID %d %d n %d dist %f sideid %d %d\n",
                        from,to,n,dist,
                        sideid[k][0],sideid[k][1]));

      from = to;
    }
    if (TriangulatePatch(Heap,p,bndp,sides,corners,h,npc,
                         cornerid,local,sideid,siden))
      return(-1);
  }

  return(nodeid);
}
#endif

MESH *BVP_GenerateMesh (HEAP *Heap, BVP *aBVP, INT argc, char **argv)
{
  STD_BVP *theBVP;
  BND_PS *ps;
  INT i,j,m,n;
  MESH *mesh;
  CoeffProcPtr coeff;
  float h;
  int ic;

  theBVP = GetSTD_BVP(aBVP);

  PRINTDEBUG(dom,1,(" bvp nsubcf %x\n",theBVP->numOfSubdomains));
  for (i=currBVP->sideoffset; i<currBVP->sideoffset+theBVP->nsides; i++)
    PRINTDEBUG(dom,1,(" addr %x\n",PARAM_PATCH_BS(currBVP->patches[i])));

  mesh = (MESH *) GetTmpMem(Heap,sizeof(MESH));
  if (mesh == NULL)
    return(NULL);

  coeff = NULL;
  h = 0.0;
  for (i=1; i<argc; i++)
    if (argv[i][0] == 'h')
    {
      if (sscanf(argv[i],"h %f",&h) != 1)
        h = 0.0;
    }
    else if (argv[i][0] == 'm')
    {
      if (sscanf(argv[i],"m %d",&ic) == 1)
        if (BVP_SetCoeffFct((BVP *)theBVP,ic,&coeff))
          coeff = NULL;
    }

  mesh->nInnP = 0;
  mesh->nElements = NULL;
  mesh->Element_corners = NULL;
  mesh->Element_corner_ids = NULL;
  mesh->nSubDomains = theBVP->numOfSubdomains;
  mesh->nSides =
    (INT *) GetTmpMem(Heap,(theBVP->numOfSubdomains+1)*sizeof(INT));
  if (mesh->nSides == NULL)
    return(NULL);
  for (i=0; i<=mesh->nSubDomains; i++)
    mesh->nSides[i] = 0;
  mesh->Side_corners =
    (INT **) GetTmpMem(Heap,(theBVP->numOfSubdomains+1)*sizeof(INT*));
  if (mesh->Side_corners == NULL)
    return(NULL);
  mesh->Side_corner_ids =
    (INT ***) GetTmpMem(Heap,(theBVP->numOfSubdomains+1)*sizeof(INT**));
  if (mesh->Side_corner_ids == NULL)
    return(NULL);

  n = theBVP->ncorners;
  if (coeff != NULL)
    n = GenerateBnodes(Heap,theBVP,NULL,mesh->nSides,NULL,coeff);
  else if (h>0)
    n = GenerateBnodes_h(Heap,theBVP,NULL,mesh->nSides,NULL,(DOUBLE) h);

  if (n == -1)
    return(NULL);
  mesh->nBndP = n;
  mesh->theBndPs = (BNDP **) GetTmpMem(Heap,n*sizeof(BNDP *));
  if (mesh->theBndPs == NULL)
    return(NULL);

  for (i=0; i<=mesh->nSubDomains; i++)
    PRINTDEBUG(dom,1,("mesh sd i %d m %d\n",i,mesh->nSides[i]));
  PRINTDEBUG(dom,1,("mesh n %d\n",mesh->nBndP));
  PRINTDEBUG(dom,1,("ncorners %d\n",theBVP->ncorners));
  for (i=0; i<theBVP->ncorners; i++)
    PRINTDEBUG(dom,1,(" id %d\n",PATCH_ID(theBVP->patches[i])));

  if (CreateCornerPoints(Heap,theBVP,mesh->theBndPs))
    return(NULL);

  for (i=0; i<theBVP->ncorners; i++)
    PRINTDEBUG(dom,1,("   i %d  patch id %d\n",i,
                      ((BND_PS *)(mesh->theBndPs[i]))->patch_id));

  for (i=0; i<=mesh->nSubDomains; i++)
  {
    m = mesh->nSides[i];
    if (m == 0)
    {
      mesh->Side_corners[i] = NULL;
      mesh->Side_corner_ids[i] = NULL;
    }
    else
    {
      mesh->Side_corners[i] = (INT *) GetTmpMem(Heap,m*sizeof(INT));
      if (mesh->Side_corners[i] == NULL)
        return(NULL);
      mesh->Side_corner_ids[i] = (INT **) GetTmpMem(Heap,m*sizeof(INT*));
      if (mesh->Side_corner_ids[i] == NULL)
        return(NULL);
    }
    for (j=0; j<m; j++)
    {
      mesh->Side_corners[i][j] = DIM;

      PRINTDEBUG(dom,1,("  i %d m %d j %d size %d\n",i,m,j,
                        mesh->Side_corners[i][j]));

      mesh->Side_corner_ids[i][j] =
        (INT *)GetTmpMem(Heap,DIM*sizeof(INT));
      if (mesh->Side_corner_ids[i][j] == NULL)
        return(NULL);
    }
  }

  if (coeff != NULL)
    n = GenerateBnodes(Heap,theBVP,mesh->theBndPs,
                       mesh->nSides,mesh->Side_corner_ids,coeff);
  else if (h>0)
    n = GenerateBnodes_h(Heap,theBVP,mesh->theBndPs,
                         mesh->nSides,mesh->Side_corner_ids,(DOUBLE) h);
  if (n == -1)
    return(NULL);

  for (i=0; i<=mesh->nSubDomains; i++)
  {
    m = mesh->nSides[i];
    PRINTDEBUG(dom,1,("mesh sd i %d m %d\n",i,m));
    for (j=0; j<m; j++)
      PRINTDEBUG(dom,1,("mesh face j %d (%d,%d,%d)\n",j,
                        mesh->Side_corner_ids[i][j][0],
                        mesh->Side_corner_ids[i][j][1],
                        mesh->Side_corner_ids[i][j][2]));
  }

  for (i=0; i<mesh->nBndP; i++)
    PRINTDEBUG(dom,1,("   i %d  patch id %d\n",i,
                      ((BND_PS *)(mesh->theBndPs[i]))->patch_id));

  return (mesh);
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/* functions for BNDS														*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDS_Global - gets global coordinates of local position

   SYNOPSIS:
   INT BNDS_Local2Global (BNDS *aBndS, COORD *local, COORD *global);

   PARAMETERS:
   .  aBndS - BNDS structure
   .  local - local coordinate on BNDS
   .  global - global coordinate

   DESCRIPTION:
   This function gets global coordinates of local position on BNDS

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BNDS_Global (BNDS *aBndS, COORD *local, COORD *global)
{
  BND_PS *ps;
  PATCH *p;
  COORD lambda[DIM_OF_BND];
  INT i,j;

  ps = (BND_PS *)aBndS;
  p = currBVP->patches[ps->patch_id];
  if (p == NULL)
    return(1);

  PRINTDEBUG(dom,1,(" Bnds global loc %f pid %d\n",
                    local[0],PATCH_ID(p)));

  if (PATCH_TYPE(p) == PARAMETRIC_PATCH_TYPE)
      #ifdef __TWODIM__
    lambda[0] = (1.0-local[0])*ps->local[0][0]+local[0]*ps->local[1][0];
      #endif
      #ifdef __THREEDIM__
    switch (ps->n)
    {
    case 3 :
      lambda[0] = (1.0-local[0]-local[1]) * ps->local[0][0]
                  + local[0]*ps->local[1][0] + local[1]*ps->local[2][0];
      lambda[1] = (1.0-local[0]-local[1]) * ps->local[0][1]
                  + local[0]*ps->local[1][1] + local[1]*ps->local[2][1];
      break;
    case 4 :
      lambda[0] = (1.0-local[0])*(1.0-local[1]) * ps->local[0][0]
                  + local[0]*(1.0-local[1])*ps->local[1][0]
                  + local[0]*local[1]*ps->local[2][0]
                  + (1.0-local[0])*local[1]*ps->local[3][0];
      lambda[1] = (1.0-local[0])*(1.0-local[1]) * ps->local[0][1]
                  + local[0]*(1.0-local[1])*ps->local[1][1]
                  + local[0]*local[1]*ps->local[2][1]
                  + (1.0-local[0])*local[1]*ps->local[3][1];
      break;
    }
       #endif
    else
      return(1);



  (*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),lambda,global);
  PRINTDEBUG(dom,1,(" Bnds global loc %f lam %f glo %f %f\n",
                    local[0],lambda[0],global[0],global[1]));

  return((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),lambda,global));
}

/****************************************************************************/
/*D
   BNDS_BndCond - gets global coordinates of local position

   SYNOPSIS:
   INT BNDS_BndCond (BNDS *aBndS, COORD *local, COORD *in,
   INT *type, DOUBLE *value);

   PARAMETERS:
   .  aBndS - BNDS structure
   .  local - local coordinate on BNDS
   .  time - time parameter or NULL
   .  type - type of bnd cond
   .  value - values

   DESCRIPTION:
   This function gets bnd conditions of local position on BNDS

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BNDS_BndCond (BNDS *aBndS, COORD *local, COORD *in, DOUBLE *value, INT *type)
{
  BND_PS *ps;
  PATCH *p;
  COORD lambda[DIM_OF_BND],global[DIM];
  INT i,j;

  PRINTDEBUG(dom,1,(" BndCond loc %f\n",local[0]));

  ps = (BND_PS *)aBndS;
  if (ps == NULL)
    return(1);

  p = currBVP->patches[ps->patch_id];

  PRINTDEBUG(dom,1,(" BndCond %d %x\n",PATCH_TYPE(p),p));

  if (PATCH_TYPE(p) == PARAMETRIC_PATCH_TYPE)
      #ifdef __TWODIM__
    lambda[0] = (1.0-local[0])*ps->local[0][0] + local[0]*ps->local[1][0];
      #endif
      #ifdef __THREEDIM__
    switch (ps->n)
    {
    case 3 :
      lambda[0] = (1.0-local[0]-local[1]) * ps->local[0][0]
                  + local[0]*ps->local[1][0] + local[1]*ps->local[2][0];
      lambda[1] = (1.0-local[0]-local[1]) * ps->local[0][1]
                  + local[0]*ps->local[1][1] + local[1]*ps->local[2][1];
      break;
    case 4 :
      lambda[0] = (1.0-local[0])*(1.0-local[1]) * ps->local[0][0]
                  + local[0]*(1.0-local[1])*ps->local[1][0]
                  + local[0]*local[1]*ps->local[2][0]
                  + (1.0-local[0])*local[1]*ps->local[3][0];
      lambda[1] = (1.0-local[0])*(1.0-local[1]) * ps->local[0][1]
                  + local[0]*(1.0-local[1])*ps->local[1][1]
                  + local[0]*local[1]*ps->local[2][1]
                  + (1.0-local[0])*local[1]*ps->local[3][1];
      break;
    }
       #endif
    else
      return(1);

  if (currBVP->GeneralBndCond != NULL)
  {
    type[0] = PATCH_ID(p) - currBVP->sideoffset;
    if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),lambda,global))
      return(1);
    if (in == NULL)
      return((*(currBVP->GeneralBndCond))(NULL,NULL,global,value,type));
    for (i=0; i<DIM; i++)
      in[i] = global[i];
    return((*(currBVP->GeneralBndCond))(NULL,NULL,in,value,type));
  }

  if (in == NULL)
    return((*PARAM_PATCH_BC (p))(PARAM_PATCH_BCD(p),NULL,lambda,value,type));

  for (i=0; i<DIM_OF_BND; i++)
    in[i] = lambda[i];
  return((*PARAM_PATCH_BC (p))(PARAM_PATCH_BCD(p),NULL,in,value,type));
}

INT BNDS_BndSDesc (BNDS *theBndS, INT *left, INT *right)
{
  BND_PS *ps;
  PATCH *p;

  ps = (BND_PS *)theBndS;
  p = currBVP->patches[ps->patch_id];

  if (PATCH_TYPE(p) == PARAMETRIC_PATCH_TYPE)
  {
    *left = PARAM_PATCH_LEFT(p);
    *right = PARAM_PATCH_RIGHT(p);
    return(0);
  }

  return(1);
}

/****************************************************************************/
/*D
   BNDS_CreateBndP - set BNDP of BNDS

   SYNOPSIS:
   BNDP *BNDS_CreateBndP (HEAP *Heap, BNDS *aBndS, COORD *local);

   PARAMETERS:
   .  aBndS - BNDS structure
   .  local - local coordinate on BNDS
   .  size - size used for aBndP

   DESCRIPTION:
   This function sets a boundary point (BNDP) of the BNDS

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

BNDP *BNDS_CreateBndP (HEAP *Heap, BNDS *aBndS, COORD *local)
{
  BND_PS *ps,*pp;
  PATCH *p;
  COORD lambda[DIM_OF_BND];
  INT n,i,j;

  if (aBndS == NULL)
    return(NULL);
  ps = (BND_PS *)aBndS;
  p = currBVP->patches[ps->patch_id];

  pp = (BND_PS *)GetFreelistMemory(Heap,sizeof(BND_PS));
  if (pp == NULL)
    return(NULL);

  pp->patch_id = ps->patch_id;
  pp->n = 1;

  if (PATCH_TYPE(p) == PARAMETRIC_PATCH_TYPE)
      #ifdef __TWODIM__
    pp->local[0][0]=(1.0-local[0])*ps->local[0][0]+local[0]*ps->local[1][0];
      #endif
      #ifdef __THREEDIM__
    switch (ps->n)
    {
    case 3 :
      pp->local[0][0] = (1.0-local[0]-local[1]) * ps->local[0][0]
                        + local[0]*ps->local[1][0] + local[1]*ps->local[2][0];
      pp->local[0][1] = (1.0-local[0]-local[1]) * ps->local[0][1]
                        + local[0]*ps->local[1][1] + local[1]*ps->local[2][1];
      break;
    case 4 :
      pp->local[0][0] = (1.0-local[0])*(1.0-local[1]) * ps->local[0][0]
                        + local[0]*(1.0-local[1])*ps->local[1][0]
                        + local[0]*local[1]*ps->local[2][0]
                        + (1.0-local[0])*local[1]*ps->local[3][0];
      pp->local[0][1] = (1.0-local[0])*(1.0-local[1]) * ps->local[0][1]
                        + local[0]*(1.0-local[1])*ps->local[1][1]
                        + local[0]*local[1]*ps->local[2][1]
                        + (1.0-local[0])*local[1]*ps->local[3][1];
      break;
    }
       #endif
    else
      return(NULL);


  PRINTDEBUG(dom,1,(" BNDP s %d\n",pp->patch_id));

  return((BNDP *)pp);
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/* functions for BNDP														*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_Global - gets global coordinates of BNDP

   SYNOPSIS:
   INT BNDP_Global (BNDP *aBndP, COORD *global);

   PARAMETERS:
   .  aBndP - BNDP structure
   .  global - global coordinate

   DESCRIPTION:
   This function gets global coordinates of BNDP

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BNDP_Global (BNDP *aBndP, COORD *global)
{
  BND_PS *ps,*pp;
  PATCH *p,*s;
  COORD lambda[DIM_OF_BND];
  INT j,k;
  COORD pglobal[DIM];

  ps = (BND_PS *)aBndP;
  p = currBVP->patches[ps->patch_id];

  PRINTDEBUG(dom,1,(" bndp pid %d %d %d\n",ps->patch_id,
                    PATCH_ID(p),PATCH_TYPE(p)));

  switch (PATCH_TYPE(p))
  {
  case PARAMETRIC_PATCH_TYPE :
    (*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),ps->local[0],global);
    PRINTDEBUG(dom,1,(" bndp param  loc %f %f gl %f %f %f\n",
                      ps->local[0][0],
                      ps->local[0][1],
                      global[0],global[1],global[2]));
    return((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),ps->local[0],global));
  case POINT_PATCH_TYPE :
    s = currBVP->patches[POINT_PATCH_PID(p,0)];
    PRINTDEBUG(dom,1,(" bndp n %d %d loc %f %f gl \n",
                      POINT_PATCH_N(p),
                      POINT_PATCH_PID(p,0),
                      ps->local[0][0],
                      ps->local[0][1]));
    if ((*PARAM_PATCH_BS (s))(PARAM_PATCH_BSD(s),ps->local[0],global))
      return(1);
    PRINTDEBUG(dom,1,(" bndp n %d %d loc %f %f gl %f %f %f\n",
                      POINT_PATCH_N(p),
                      POINT_PATCH_PID(p,0),
                      ps->local[0][0],
                      ps->local[0][1],
                      global[0],global[1],global[2]));
    for (j=1; j<POINT_PATCH_N(p); j++)
    {
      s = currBVP->patches[POINT_PATCH_PID(p,j)];
      if ((*PARAM_PATCH_BS (s))(PARAM_PATCH_BSD(s),ps->local[j],pglobal))
        return(1);
      PRINTDEBUG(dom,1,(" bndp    j %d %d loc %f %f gl %f %f %f\n",j,
                        POINT_PATCH_PID(p,j),
                        ps->local[j][0],
                        ps->local[j][1],
                        pglobal[0],pglobal[1],pglobal[2]));
      for (k=0; k<DIM; k++)
        if (ABS(pglobal[k] - global[k]) > SMALL_DIFF)
          return(1);
    }
    return(0);
      #ifdef __THREEDIM__
  case LINE_PATCH_TYPE :
    s = currBVP->patches[LINE_PATCH_PID(p,0)];
    if ((*PARAM_PATCH_BS (s))(PARAM_PATCH_BSD(s),ps->local[0],global))
      return(1);
    PRINTDEBUG(dom,1,(" bndp    n %d %d loc %f %f gl %f %f %f\n",
                      POINT_PATCH_N(p),
                      LINE_PATCH_PID(p,0),
                      ps->local[0][0],
                      ps->local[0][1],
                      global[0],global[1],global[2]));
    for (j=1; j<LINE_PATCH_N(p); j++)
    {
      s = currBVP->patches[LINE_PATCH_PID(p,j)];
      if ((*PARAM_PATCH_BS (s))(PARAM_PATCH_BSD(s),ps->local[j],pglobal))
        return(1);
      PRINTDEBUG(dom,1,(" bndp    j %d %d loc %f %f gl %f %f %f\n",j,
                        LINE_PATCH_PID(p,j),
                        ps->local[j][0],
                        ps->local[j][1],
                        pglobal[0],pglobal[1],pglobal[2]));
      for (k=0; k<DIM; k++)
        if (ABS(pglobal[k] - global[k]) > SMALL_DIFF)
          return(1);
    }
    return(0);
      #endif
  }

  return(1);
}

/****************************************************************************/
/*D
   BNDP_BndPDesc - sets descriptor for BNDP

   SYNOPSIS:
   INT BNDP_BndPDesc (BNDP *theBndP, INT *move);

   PARAMETERS:
   .  aBndP - BNDP structure
   .  move

   DESCRIPTION:
   This function sets the descriptor for a BNDP

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BNDP_BndPDesc (BNDP *theBndP, INT *move)
{
  BND_PS *ps,*pp;
  PATCH *p;
  COORD lambda[DIM_OF_BND];
  INT n,i,j;

  ps = (BND_PS *)theBndP;
  p = currBVP->patches[ps->patch_id];

  switch (PATCH_TYPE(p))
  {
  case PARAMETRIC_PATCH_TYPE :
    *move = DIM_OF_BND;
    return(0);
  case POINT_PATCH_TYPE :
    *move = 0;
    return(0);
      #ifdef __THREEDIM__
  case LINE_PATCH_TYPE :
    *move = 1;
    return(0);
      #endif
  }

  return(1);
}

/****************************************************************************/
/*D
   BNDP_CreateBndS - sets BNDS from a nb of BNDPs

   SYNOPSIS:
   BNDS *BNDP_CreateBndS (HEAP *Heap, BNDP **aBndP, INT n);

   PARAMETERS:
   .  aBndP - ptr to list of BNDP structures
   .  n - nb of BNDPs

   DESCRIPTION:
   This function sets a BNDS from n BNDPs

   RETURN VALUE:
   BNDS *
   .n   pointer
   .n   NULL if the points describe an inner side
   D*/
/****************************************************************************/

BNDS *BNDP_CreateBndS (HEAP *Heap, BNDP **aBndP, INT n)
{
  BND_PS *bp[4],*bs;
  PATCH *p[4];
  COORD *lambda[4];
  INT i,j,k,l,pid;

  PRINTDEBUG(dom,1,("Create BNDS:\n"));
  for (i=0; i<n; i++)
  {
    bp[i] = (BND_PS *)aBndP[i];
    p[i] = currBVP->patches[bp[i]->patch_id];

    PRINTDEBUG(dom,1,(" bp %d p %d n %d\n",
                      bp[i]->patch_id,PATCH_ID(p[i]),n));
    for (l=0; l<GetNumberOfPatches(p[i]); l++)
      PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p[i],l)));
  }

  pid = -1;
    #ifdef __TWODIM__
  if (n!=2)
    return(NULL);
  for (i=0; i<GetNumberOfPatches(p[0]); i++)
    for (j=0; j<GetNumberOfPatches(p[1]); j++)
    {
      PRINTDEBUG(dom,1,(" pid i, j %d %d %d %d\n",i,j,
                        GetPatchId(p[0],i),GetPatchId(p[1],j)));
      if (GetPatchId(p[0],i) == GetPatchId(p[1],j))
      {
        pid = GetPatchId(p[0],i);
        lambda[0] = bp[0]->local[i];
        lambda[1] = bp[1]->local[j];
        break;
        PRINTDEBUG(dom,1,(" pid %d \n",pid));
      }
    }
    #endif
    #ifdef __THREEDIM__
  switch (n)
  {
  case 3 :
    for (i=0; i<GetNumberOfPatches(p[0]); i++)
      for (j=0; j<GetNumberOfPatches(p[1]); j++)
        if (GetPatchId(p[0],i) == GetPatchId(p[1],j))
          for (k=0; k<GetNumberOfPatches(p[2]); k++)
            if (GetPatchId(p[0],i) == GetPatchId(p[2],k))
            {
              pid = GetPatchId(p[0],i);
              lambda[0] = bp[0]->local[i];
              lambda[1] = bp[1]->local[j];
              lambda[2] = bp[2]->local[k];
              break;
            }
    break;
  case 4 :
    for (i=0; i<GetNumberOfPatches(p[0]); i++)
      for (j=0; j<GetNumberOfPatches(p[1]); j++)
        if (GetPatchId(p[0],i) == GetPatchId(p[1],j))
          for (k=0; k<GetNumberOfPatches(p[2]); k++)
            if (GetPatchId(p[0],i) == GetPatchId(p[2],k))
              for (l=0; l<GetNumberOfPatches(p[3]); l++)
                if (GetPatchId(p[0],i) == GetPatchId(p[3],l))
                {
                  pid = GetPatchId(p[0],i);
                  lambda[0] = bp[0]->local[i];
                  lambda[1] = bp[1]->local[j];
                  lambda[2] = bp[2]->local[k];
                  lambda[3] = bp[3]->local[l];
                  break;
                }
    break;
  }
    #endif

  if (pid == -1)
    return(NULL);

  bs = (BND_PS *)GetFreelistMemory(Heap,(n-1)*sizeof(COORD_BND_VECTOR)
                                   + sizeof(BND_PS));
  if (bs == NULL)
    return(NULL);
  bs->n = n;
  bs->patch_id = pid;
  for (i=0; i<n; i++)
    for (j=0; j<DIM_OF_BND; j++)
      bs->local[i][j] = lambda[i][j];

  for (i=0; i<n; i++)
    for (j=0; j<DIM_OF_BND; j++)
      PRINTDEBUG(dom,1,(" bnds i, j %d %d %f\n",i,j,lambda[i][j]));

  PRINTDEBUG(dom,1,(" Create BNDS %x %d\n",bs,pid));

  return((BNDS *)bs);
}

/****************************************************************************/
/*D
   BNDP_CreateBndP - sets BNDP from a two of BNDPs

   SYNOPSIS:
   BNDP *BNDP_CreateBndP (HEAP *Heap, BNDP *aBndP0, BNDP *aBndP1, COORD lcoord);

   PARAMETERS:
   .  aBndP0 - first BNDP
   .  aBndP1 - second BNDP

   DESCRIPTION:
   This function sets a BNDP from two BNDPs

   RETURN VALUE:
   BNDS *
   .n   pointer
   .n   NULL if the points describe an inner point
   D*/
/****************************************************************************/

BNDP *BNDP_CreateBndP (HEAP *Heap, BNDP *aBndP0, BNDP *aBndP1, COORD lcoord)
{
  BND_PS *bp0,*bp1,*bp;
  PATCH *p0,*p1,*p;
  COORD *lambda0,*lambda1;
  INT i,j,k,l,pid,cnt,cnt1;

  bp0 = (BND_PS *)aBndP0;
  bp1 = (BND_PS *)aBndP1;

  if ((bp0 == NULL) || (bp1 == NULL))
    return(NULL);

  p0 = currBVP->patches[bp0->patch_id];
  p1 = currBVP->patches[bp1->patch_id];

  PRINTDEBUG(dom,1,("   bp0 %d pid %d\n",
                    bp0->patch_id,PATCH_ID(p0)));
  for (l=0; l<GetNumberOfPatches(p0); l++)
    PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p0,l)));
  PRINTDEBUG(dom,1,("   bp1 %d pid %d\n",
                    bp1->patch_id,PATCH_ID(p1)));
  for (l=0; l<GetNumberOfPatches(p1); l++)
    PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p1,l)));

  cnt = GetNumberOfCommonPatches(p0,p1);
  if (cnt == 0)
    return(NULL);

  bp = (BND_PS *)GetFreelistMemory(Heap,(cnt-1)*sizeof(COORD_BND_VECTOR)
                                   + sizeof(BND_PS));
  if (bp == NULL)
    return(NULL);
  bp->n = cnt;

    #ifdef __THREEDIM__
  if (cnt > 1)
  {
    k = GetCommonLinePatchId (p0,p1);
    if ((k<currBVP->ncorners) || (k>=currBVP->sideoffset))
      return(NULL);
    p = currBVP->patches[k];
    bp->patch_id = k;

    PRINTDEBUG(dom,1,(" Create BNDP line %d cnt %d\n",k,cnt));
    PRINTDEBUG(dom,1,("   bp0 %d pid %d\n",
                      bp0->patch_id,PATCH_ID(p0)));
    for (l=0; l<GetNumberOfPatches(p0); l++)
      PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p0,l)));
    PRINTDEBUG(dom,1,("   bp1 %d pid %d\n",
                      bp1->patch_id,PATCH_ID(p1)));
    for (l=0; l<GetNumberOfPatches(p1); l++)
      PRINTDEBUG(dom,1,("    bp pid %d\n",GetPatchId(p1,l)));

    for (l=0; l<LINE_PATCH_N(p); l++)
    {
      for (i=0; i<GetNumberOfPatches(p0); i++)
        if (GetPatchId(p0,i) == LINE_PATCH_PID(p,l))
          for (j=0; j<GetNumberOfPatches(p1); j++)
            if (GetPatchId(p1,j) == LINE_PATCH_PID(p,l))
              for (k=0; k<DIM_OF_BND; k++)
                bp->local[l][k] = (1.0-lcoord)*bp0->local[i][k]
                                  +lcoord*bp1->local[j][k];
    }

    for (l=0; l<LINE_PATCH_N(p); l++)
      PRINTDEBUG(dom,1,(" Create BNDP line %d l %d %f %f\n",
                        LINE_PATCH_PID(p,l),l,
                        bp->local[l][0],bp->local[l][1]));

    return((BNDP *)bp);
  }
    #endif

  for (i=0; i<GetNumberOfPatches(p0); i++)
    for (j=0; j<GetNumberOfPatches(p1); j++)
      if (GetPatchId(p0,i) == GetPatchId(p1,j))
      {
        bp->patch_id = GetPatchId(p0,i);
        for (k=0; k<DIM_OF_BND; k++)
          bp->local[0][k] =(1.0-lcoord)*bp0->local[i][k]
                            +lcoord * bp1->local[j][k];
        PRINTDEBUG(dom,1,(" Create BNDP param %d \n",GetPatchId(p0,i)));

        break;
      }

  return((BNDP *)bp);
}

INT BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size)
{
  BND_PS *bp;
  PATCH *p;
  INT pid;

  bp = (BND_PS *)theBndP;

  if (bp == NULL)
    return(1);

  pid = bp->patch_id;
  p = currBVP->patches[pid];

  switch (PATCH_TYPE(p))
  {
  case PARAMETRIC_PATCH_TYPE :
    pid -= currBVP->sideoffset;
    break;
  case POINT_PATCH_TYPE :
    pid = POINT_PATCH_PID(p,0) - currBVP->sideoffset;
    break;
      #ifdef __THREEDIM__
  case LINE_PATCH_TYPE :
    pid = LINE_PATCH_PID(p,0) - currBVP->sideoffset;
    break;
      #endif
  }

  PRINTDEBUG(dom,1,(" Insert pid %d %d\n",bp->patch_id,pid));

    #ifdef __TWODIM__
  if (sprintf(data,"bn %d %f",pid,(float) bp->local[0][0])>max_data_size)
    return(1);
    #endif

    #ifdef __THREEDIM__
  if (sprintf(data,"bn %d %f %f",(int) pid,
              (float) bp->local[0][0],(float) bp->local[0][1])>max_data_size)
    return(1);
    #endif

  return(0);
}

/****************************************************************************/
/*D
   BNDP_BndCond - gets global coordinates of local position

   SYNOPSIS:
   INT BNDP_BndCond (BNDP *aBndP, INT *n, INT i,
   COORD *in, DOUBLE *value, INT *type);

   PARAMETERS:
   .  aBndP - BNDP structure
   .  n - number of BNDS
   .  in - input vector
   .  type - type of bnd cond
   .  value - values

   DESCRIPTION:
   This function gets bnd conditions of local position on BNDS

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BNDP_BndCond (BNDP *aBndP, INT *n, INT i, COORD *in, DOUBLE *value, INT *type)
{
  BND_PS *ps;
  PATCH *p;
  COORD global[DIM];
  COORD *local;
  INT j;

  if (i < 0)
    return(1);

  ps = (BND_PS *)aBndP;
  if (ps == NULL)
    return(1);

  p = currBVP->patches[ps->patch_id];

  switch (PATCH_TYPE(p))
  {
  case PARAMETRIC_PATCH_TYPE :
    n[0] = 1;
    local = ps->local[0];
    break;
  case POINT_PATCH_TYPE :
    n[0] = POINT_PATCH_N(p);
    if (i >= n[0])
      return(1);
    p = currBVP->patches[POINT_PATCH_PID(p,i)];
    local = ps->local[i];
    break;
      #ifdef __THREEDIM__
  case LINE_PATCH_TYPE :
    n[0] = LINE_PATCH_N(p);
    if (i >= n[0])
      return(1);
    p = currBVP->patches[LINE_PATCH_PID(p,i)];
    local = ps->local[i];
    break;
      #endif
  }

  PRINTDEBUG(dom,1,(" BndCond %d loc %f %f\n",
                    PATCH_ID(p),local[0],local[1]));

  if (PATCH_TYPE(p) != PARAMETRIC_PATCH_TYPE)
    return(1);

  if (currBVP->GeneralBndCond != NULL)
  {
    type[0] = PATCH_ID(p) - currBVP->sideoffset;
    if ((*PARAM_PATCH_BS (p))(PARAM_PATCH_BSD(p),local,global))
      return(1);
    if (in == NULL)
      return((*(currBVP->GeneralBndCond))(NULL,NULL,global,value,type));
    for (i=0; i<DIM; i++)
      in[i] = global[i];
    return((*(currBVP->GeneralBndCond))(NULL,NULL,in,value,type));
  }

  if (in == NULL)
    return((*PARAM_PATCH_BC (p))(PARAM_PATCH_BCD(p),NULL,local,value,type));

  for (i=0; i<DIM_OF_BND; i++)
    in[i] = local[i];
  return((*PARAM_PATCH_BC (p))(PARAM_PATCH_BCD(p),NULL,in,value,type));
}

INT BNDP_Dispose (HEAP *Heap, BNDP *theBndP)
{
  BND_PS *ps;

  if (theBndP == NULL)
    return(0);

  ps = (BND_PS *)theBndP;
  return(PutFreelistMemory(Heap,ps,BND_SIZE(ps)));
}

INT BNDS_Dispose (HEAP *Heap, BNDS *theBndS)
{
  BND_PS *ps;

  if (theBndS == NULL)
    return(0);

  ps = (BND_PS *)theBndS;
  return(PutFreelistMemory(Heap,ps,BND_SIZE(ps)));
}

/* the following interface functions are not available in std_domain.c */

INT BVP_Save (BVP *theBVP, char *name, INT argc, char **argv)
{
  return (1);
}

BVP *BVP_Load (char *name, INT argc, char **argv)
{
  return (NULL);
}

/****************************************************************************/
/*D
   BNDP_SaveBndP - save a BNDP

   SYNOPSIS:
   INT BNDP_SaveBndP (BNDP *theBndP, FILE *stream);

   PARAMETERS:
   .  theBndP - BNDP
   .  stream - file

   DESCRIPTION:
   This function saves a BNDP on a file.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT BNDP_SaveBndP (BNDP *BndP)
{
  BND_PS *bp;
  INT i,j;
  int iList[2];
  double dList[DIM-1];

  iList[0] = BND_PATCH_ID(BndP);
  iList[1] = BND_N(BndP);
  if (Bio_Write_mint(2,iList)) return (1);

  bp = (BND_PS *)BndP;
  for (i=0; i<BND_N(BndP); i++)
  {
    for (j=0; j<DIM-1; j++)
      dList[j] = bp->local[i][j];
    if (Bio_Write_mdouble(DIM-1,dList)) return (1);
  }

  return(0);
}

/****************************************************************************/
/*D
   BVP_LoadBndP - load a BNDP

   SYNOPSIS:
   BNDP *BNDP_LoadBndP (BVP *theBVP, HEAP *Heap, FILE *stream);

   PARAMETERS:
   .  theBndP - BNDP
   .  Heap - heap
   .  stream - file

   DESCRIPTION:
   This function loads a BNDP with the format given by BVP_SaveBndP.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

BNDP *BNDP_LoadBndP (BVP *theBVP, HEAP *Heap)
{
  BND_PS *bp;
  int i,j,pid,n;
  double local;
  int iList[2];
  double dList[DIM-1];

  if (Bio_Read_mint(2,iList)) return (NULL);
  pid = iList[0]; n = iList[1];
  bp = (BND_PS *)GetFreelistMemory(Heap,(n-1)*sizeof(COORD_BND_VECTOR) + sizeof(BND_PS));
  bp->n = n;
  bp->patch_id = pid;
  for (i=0; i<n; i++)
  {
    if (Bio_Read_mdouble(DIM-1,dList)) return (NULL);
    for (j=0; j<DIM-1; j++)
      bp->local[i][j] = dList[j];
  }

  return((BNDP *)bp);
}

/****************************************************************************/
/*D
   InitDom - Create and initialize the std_domain

   SYNOPSIS:
   INT InitDom ();

   PARAMETERS:
   .  void - no argument

   DESCRIPTION:
   This function creates the environments domain, problem and BVP.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 when error occured.
   D*/
/****************************************************************************/

INT InitDom ()
{
  int heapSize;
  char buffer[128];

  /* change to root directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitDom","could not changedir to root");
    return(__LINE__);
  }

  /* get env dir/var IDs for the problems */
  theProblemDirID  = GetNewEnvDirID();
  theBdryCondVarID = GetNewEnvVarID();

  /* install the /Domains directory */
  theDomainDirID = GetNewEnvDirID();
  if (MakeEnvItem("Domains",theProblemDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitDom","could not install '/Domains' dir");
    return(__LINE__);
  }
  theBdrySegVarID = GetNewEnvVarID();

  /* install the /BVP directory */
  theBVPDirID = GetNewEnvDirID();
  if (MakeEnvItem("BVP",theBVPDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitDom","could not install '/BVP' dir");
    return(__LINE__);
  }

  return (0);
}
