// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  std_domain.c													*/
/*																			*/
/* Purpose:   standard ug domain description                                                            */
/*																			*/
/* Author:	  Klaus Johannsen				                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   18.02.96 begin, ug version 3.1								*/
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
#include <math.h>

/* low modules */
#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "misc.h"

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



/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct bndsegcond
{
  BOUNDARY_CONDITION *theBC;
  BOUNDARY_SEGMENT *theBS;
};

struct tmp
{
  /* fields for environment directory */
  ENVVAR v;

  /* pointers to bndseg and bndcond */
  struct bndsegcond theBSC[1];
};

struct tmp2
{
  /* fields for environment directory */
  ENVVAR v;

  /* integer list */
  INT used[1];
};

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
static INT thePatchVarID;                       /* env type for patch vars				*/
static INT theTmpVarID;                         /* env type for temporary use           */

/* data for CVS */
static char rcsid[] = "$Header$";

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
/*D
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
   D*/
/****************************************************************************/

PROBLEM *GetProblem (const char * domain, const char *name)
{
  if (ChangeEnvDir("/Domains")==NULL) return(NULL);
  return((PROBLEM *) SearchEnv(name,domain,theProblemDirID,theDomainDirID));
}

/****************************************************************************/
/*
   CreateBoundaryCondition - Create a new BNDCOND structure

   SYNOPSIS:
   BOUNDARY_CONDITION *CreateBoundaryCondition (char *name, INT id,
   BndCondProcPtr theBndCond);

   PARAMETERS:
   .  name - name of the boundary condition
   .  id - identification number of condition
   .  theBndCond - the boundary conditions

   DESCRIPTION:
   This function allocates and initializes a new BNDCOND structure in the previously allocated
   PROBLEM structure.

   RETURN VALUE:
   BOUNDARY_CONDITION *
   .n      Pointer to BOUNDARY_CONDITION
   .n      NULL if if out of memory.
 */
/****************************************************************************/

#ifdef __version23__
BOUNDARY_CONDITION *CreateBoundaryCondition (char *name, INT id, BndCondProcPtr theBndCond)
{
  BOUNDARY_CONDITION *newBndCond;

  /* allocate boundary condition */
  newBndCond = (BOUNDARY_CONDITION *) MakeEnvItem (name,theBdryCondVarID,sizeof(BOUNDARY_CONDITION));
  if (newBndCond==NULL) return(NULL);

  /* fill in data */
  newBndCond->id = id;
  newBndCond->BndCond = theBndCond;

  return(newBndCond);
}
#endif

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

#ifdef __version3__
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
#endif

/****************************************************************************/
/*D
   GetFirstBoundaryCondition - Get first BNDCOND structure of `theProblem`

   SYNOPSIS:
   BOUNDARY_CONDITION *GetFirstBoundaryCondition (PROBLEM *theProblem);

   PARAMETERS:
   .  theProblem - pointer to PROBLEM

   DESCRIPTION:
   This function gets the first BNDCOND structure of a problem

   THE RETURN VALUE:
   BOUNDARY_CONDITION *
   .n      pointer to BOUNDARY_CONDITION
   .n      NULL if not found or error.
   D*/
/****************************************************************************/
BOUNDARY_CONDITION *GetFirstBoundaryCondition (PROBLEM *theProblem)
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
   GetNextBoundaryCondition - Allocate a new BNDCOND structure

   SYNOPSIS:
   BOUNDARY_CONDITION *GetNextBoundaryCondition (BOUNDARY_CONDITION *theBCond);

   PARAMETERS:
   .  theBCond - the boundary condition

   DESCRIPTION:
   This function gets next BNDCOND structure

   THE RETURN VALUE:
   BOUNDARY_CONDITION *
   .n      pointer to BOUNDARY_CONDITION
   .n      NULL if not found or error.
   D*/
/****************************************************************************/
BOUNDARY_CONDITION *GetNextBoundaryCondition (BOUNDARY_CONDITION *theBCond)
{
  ENVITEM *theItem;

  theItem = (ENVITEM *) theBCond;

  do
    theItem = NEXT_ENVITEM(theItem);
  while ((theItem!=NULL) && (ENVITEM_TYPE(theItem)!=theBdryCondVarID));

  return ((BOUNDARY_CONDITION *) theItem);
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
/*D
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
   D*/
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
/*D
   GetFirstBoundarySegment - Get first boundary segment of a domain

   SYNOPSIS:
   BOUNDARY_SEGMENT *GetFirstBoundarySegment (DOMAIN *theDomain);

   PARAMETERS:
   .  theDomain - pointer to the domain structure

   DESCRIPTION:
   This function returns the first boundary segment of a domain.

   RETURN VALUE:
   BOUNDARY_SEGMENT *
   .n     pointer to a BOUNDARY_SEGMENT
   .n     NULL if not found or error.
   D*/
/****************************************************************************/

BOUNDARY_SEGMENT *GetFirstBoundarySegment (DOMAIN *theDomain)
{
  ENVITEM *theItem;

  theItem = ENVITEM_DOWN(theDomain);

  if (ENVITEM_TYPE(theItem)==theBdrySegVarID)
    return ((BOUNDARY_SEGMENT *) theItem);
  else
    return (GetNextBoundarySegment((BOUNDARY_SEGMENT *) theItem));
}

/****************************************************************************/
/*D
   GetNextBoundarySegment - Get next boundary segment of a domain

   SYNOPSIS:
   BOUNDARY_SEGMENT *GetNextBoundarySegment (BOUNDARY_SEGMENT *theBSeg);

   PARAMETERS:
   .  theBSeg - pointer to a boundary segment

   DESCRIPTION:
   This function gets the next boundary segment of a domain.

   RETURN VALUE:
   BOUNDARY_SEGMENT *
   .n     pointer to next BOUNDARY_SEGMENT
   .n     NULL if no more found.
   D*/
/****************************************************************************/

BOUNDARY_SEGMENT *GetNextBoundarySegment (BOUNDARY_SEGMENT *theBSeg)
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

BVP *CreateBVP (char *BVPName, char *DomainName, char *ProblemName)
{
  STD_BVP *theBVP;
  STD_PATCH *thePatch;
  DOMAIN *theDomain;
  PROBLEM *theProblem;
  BOUNDARY_SEGMENT *theSegment;
  BOUNDARY_CONDITION *theBndCond;
  struct tmp *theTmp;
  struct tmp2 *theTmp2;
  INT i, j, n, maxSubDomains;
  char name[NAMESIZE];

  /* get domain and problem */
  theDomain = GetDomain(DomainName);                                      if (theDomain==NULL) return (NULL);
  theProblem = GetProblem(DomainName,ProblemName);        if (theProblem==NULL) return (NULL);

  /* change to /BVP directory */
  if (ChangeEnvDir("/BVP")==NULL) return (NULL);

  /* allocate new domain structure */
  n = (theProblem->numOfCoeffFct+theProblem->numOfUserFct-1)*sizeof(void*);
  theBVP = (BVP *) MakeEnvItem (BVPName,theBVPDirID,sizeof(STD_BVP)+n);
  if (theBVP==NULL) return(NULL);
  if (ChangeEnvDir(BVPName)==NULL) return(NULL);

  /* fill in data of domain */
  for( i = 0 ; i < DIM ; i++)
    theBVP->MidPoint[i] = theDomain->MidPoint[i];
  theBVP->radius                  = theDomain->radius;
  theBVP->numOfPatches    = theDomain->numOfSegments;
  theBVP->numOfCorners    = theDomain->numOfCorners;
  theBVP->domConvex               = theDomain->domConvex;

  /* fill in data of problem */
  theBVP->problemID               = theProblem->problemID;
  theBVP->ConfigProblem   = theProblem->ConfigProblem;
  theBVP->numOfCoeffFct   = theProblem->numOfCoeffFct;
  theBVP->numOfUserFct    = theProblem->numOfUserFct;
  for (i=0; i<theProblem->numOfCoeffFct; i++) theBVP->CU_ProcPtr[i] = theProblem->CU_ProcPtr[i];
  for (i=0; i<theProblem->numOfUserFct; i++) theBVP->CU_ProcPtr[i+theProblem->numOfCoeffFct] = theProblem->CU_ProcPtr[i+theProblem->numOfCoeffFct];

  /* create patches */
  theTmp = (struct tmp *) MakeEnvItem ("tmp",theTmpVarID,sizeof(struct tmp)+(theBVP->numOfPatches-1)*sizeof(struct bndsegcond));
  ENVITEM_LOCKED(theTmp) = 0;
  for (i=0; i<theBVP->numOfPatches; i++)
  {
    theTmp->theBSC[i].theBS = NULL;
    theTmp->theBSC[i].theBC = NULL;
  }
  maxSubDomains = 0;
  for (theSegment=GetFirstBoundarySegment(theDomain); theSegment!=NULL; theSegment = GetNextBoundarySegment(theSegment))
  {
    i = theSegment->id;
    if ((i<0)||(i>=theBVP->numOfPatches)) return(NULL);
    if (theTmp->theBSC[i].theBS!=NULL) return (NULL);
    theTmp->theBSC[i].theBS = theSegment;
    maxSubDomains = MAX(maxSubDomains,theSegment->left);
    maxSubDomains = MAX(maxSubDomains,theSegment->right);
  }
  for (theBndCond=GetFirstBoundaryCondition(theProblem); theBndCond!=NULL; theBndCond = GetNextBoundaryCondition(theBndCond))
  {
    i = theBndCond->id;
    if ((i<0)||(i>=theBVP->numOfPatches)) return(NULL);
    if (theTmp->theBSC[i].theBC!=NULL) return (NULL);
    theTmp->theBSC[i].theBC = theBndCond;
  }
  for (i=0; i<theBVP->numOfPatches; i++)
  {
    if (theTmp->theBSC[i].theBS==NULL) return (NULL);
    theSegment = theTmp->theBSC[i].theBS;
    if (theTmp->theBSC[i].theBC==NULL) return (NULL);
    theBndCond = theTmp->theBSC[i].theBC;
    strcpy(name,ENVITEM_NAME(theSegment)); strcat(name,".p");
    thePatch = (PATCH *) MakeEnvItem (name,thePatchVarID,sizeof(STD_PATCH));
    ENVITEM_LOCKED(thePatch) = 1;

    if (thePatch == NULL) return (NULL);

    /* set domain part */
    thePatch->left                  = theSegment->left;
    thePatch->right                 = theSegment->right;
    thePatch->id                    = theSegment->id;
    thePatch->segType               = theSegment->segType;
    if (DIM==2) thePatch->nPoints = 2;
    if (DIM==3) thePatch->nPoints = 4;
    for (j=0; j<thePatch->nPoints; j++) thePatch->points[j] = theSegment->points[j];
    thePatch->resolution    = theSegment->resolution;
    for (j=0; j<DIM_OF_BND; j++) thePatch->alpha[j] = theSegment->alpha[j];
    for (j=0; j<DIM_OF_BND; j++) thePatch->beta[j]  = theSegment->beta[j];
    thePatch->BndSegFunc    = theSegment->BndSegFunc;
    thePatch->bs_data               = theSegment->data;

    /* set problem part */
    thePatch->BndCond               = theBndCond->BndCond;
    thePatch->bc_data               = theBndCond->data;
  }
  if (RemoveEnvItem((ENVITEM *)theTmp)) return (NULL);
  theTmp2 = (struct tmp2 *) MakeEnvItem ("tmp2",theTmpVarID,sizeof(struct tmp2)+maxSubDomains*sizeof(INT));
  ENVITEM_LOCKED(theTmp2) = 0;
  for (i=1; i<=maxSubDomains; i++)
    theTmp2->used[i] = 0;
  for (theSegment=GetFirstBoundarySegment(theDomain); theSegment!=NULL; theSegment = GetNextBoundarySegment(theSegment))
  {
    theTmp2->used[theSegment->left] = 1;
    theTmp2->used[theSegment->right] = 1;
  }
  for (i=1; i<=maxSubDomains; i++)
    if (!theTmp2->used[i])
    {
      UserWrite("number of subdomains not consistent\n");
      if (RemoveEnvItem ((ENVITEM *)theTmp2)) return (NULL);
      return (NULL);
    }
  if (RemoveEnvItem ((ENVITEM *)theTmp2)) return (NULL);
  theBVP->numOfSubdomains = maxSubDomains;

  UserWrite("BVP "); UserWrite(BVPName); UserWrite(" installed\n");

  return ((BVP*)theBVP);
}

/****************************************************************************/
/*D
   GetBVP -  Get pointer to a BVP structure

   SYNOPSIS:
   BVP *GetBVP (char *name);

   PARAMETERS:
   .  name - name of the BoundaryValueProblem

   DESCRIPTION:
   This function searches a BoundaryValueProblem structure.

   RETURN VALUE:
   BVP *
   .n      pointer to BoundaryValueProblem
   .n      NULL if not found or error.
   D*/
/****************************************************************************/

BVP     *GetBVP (char *name)
{
  return((BVP *) SearchEnv(name,"/BVP",theBVPDirID,theBVPDirID));
}

/****************************************************************************/
/*D
   BVP_GetBVPDesc -  Get BVP description

   SYNOPSIS:
   INT BVP_GetBVPDesc (BVP *theBVP, BVP_DESC *theBVPDesc)

   PARAMETERS:
   .  theBVP - the BoundaryValueProblem
   .  theDomainDesc - domain description

   DESCRIPTION:
   This function fills the theDomainDesc structure with information about the
   domain.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT BVP_GetBVPDesc (BVP *aBVP, BVP_DESC *theBVPDesc)
{
  STD_BVP *theBVP;
  INT i;

  /* cast */
  theBVP = (STD_BVP*)aBVP;

  /* general part */
  strcpy(theBVPDesc->name,ENVITEM_NAME(theBVP));

  /* the domain part */
  for (i=0; i<DIM; i++) theBVPDesc->midpoint[i] = theBVP->MidPoint[i];
  theBVPDesc->radius              = theBVP->radius;
  theBVPDesc->nCorners            = theBVP->numOfCorners;
  theBVPDesc->convex              = theBVP->domConvex;
  theBVPDesc->nPatches            = theBVP->numOfPatches;
  theBVPDesc->nSubDomains         = theBVP->numOfSubdomains;
  theBVPDesc->numOfCoeffFct       = theBVP->numOfCoeffFct;
  theBVPDesc->numOfUserFct        = theBVP->numOfUserFct;
  theBVPDesc->ConfigProc          = theBVP->ConfigProblem;
  theBVPDesc->id                  = theBVP->problemID;

  return (0);
}

/****************************************************************************/
/*D
   BVP_GetCoeffFct -  Get BVP coefficient fct ptrs

   SYNOPSIS:
   INT BVP_GetCoeffFct (BVP *aBVP, INT n, CoeffProcPtr *CoeffFct);

   PARAMETERS:
   .  aBVP - the BoundaryValueProblem
   .  n - number of coeff function: >=0 gives the special coeff fct of this number
                                 ==-1 gives all coeff fcts

   DESCRIPTION:
   This function fills in 'CoeffFct[0], CoeffFct[1], ... one or more coeff fcts

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT BVP_GetCoeffFct (BVP *aBVP, INT n, CoeffProcPtr *CoeffFct)
{
  STD_BVP *theBVP;
  INT i;

  /* cast */
  theBVP = (STD_BVP*)aBVP;

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
   BVP_GetUserFct -  Get BVP user fct ptrs

   SYNOPSIS:
   INT BVP_GetUserFct (BVP *aBVP, INT n, UserProcPtr *UserFct);

   PARAMETERS:
   .  aBVP - the BoundaryValueProblem
   .  n - number of user function: >=0 gives the special user fct of this number
                                ==-1 gives all user fcts

   DESCRIPTION:
   This function fills in 'UserFct[0], UserFct[1], ... one or more user fcts

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT BVP_GetUserFct (BVP *aBVP, INT n, UserProcPtr *UserFct)
{
  STD_BVP *theBVP;
  INT i;

  /* cast */
  theBVP = (STD_BVP*)aBVP;

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
   BVP_GetNextPatch - Get next PATCH structure from previous

   SYNOPSIS:
   PATCH *BVP_GetNextPatch (BVP *theBVP, PATCH *thePatch);

   PARAMETERS:
   .  theBVP - the BoundaryValueProblem
   .  thePatch - the previous patch

   DESCRIPTION:
   This function gets next PATCH structure

   THE RETURN VALUE:
   PATCH *
   .n      pointer to PATCH
   .n      NULL if not found or error.
   D*/
/****************************************************************************/
PATCH *BVP_GetNextPatch (BVP *theBVP, PATCH *thePatch)
{
  ENVITEM *theItem;

  theItem = (ENVITEM *) thePatch;

  do
    theItem = NEXT_ENVITEM(theItem);
  while ((theItem!=NULL) && (ENVITEM_TYPE(theItem)!=thePatchVarID));

  return ((PATCH *) theItem);
}

/****************************************************************************/
/*D
   BVP_GetFirstPatch - Get first PATCH structure of BVP

   SYNOPSIS:
   PATCH *BVP_GetFirstPatch (BVP *theBVP);

   PARAMETERS:
   .  theBVP - the BoundaryValueProblem

   DESCRIPTION:
   This function gets the first Patch structure of a BVP.

   THE RETURN VALUE:
   PATCH *
   .n      pointer to PATCH
   .n      NULL if not found or error.
   D*/
/****************************************************************************/
PATCH *BVP_GetFirstPatch (BVP *theBVP)
{
  ENVITEM *theItem;

  theItem = ENVITEM_DOWN(theBVP);

  if (ENVITEM_TYPE(theItem)==thePatchVarID)
    return ((PATCH *) theItem);
  else
    return (BVP_GetNextPatch(theBVP, (PATCH *)theItem));
}

/****************************************************************************/
/*D
   Patch_global2local - Project global point onto patch

   SYNOPSIS:
   INT Patch_global2local (PATCH *thePatch, COORD *global, COORD *local);

   PARAMETERS:
   .  thePatch - the patch

   DESCRIPTION:
   This function projects a point in physical space onto the patch

   THE RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT Patch_global2local (PATCH *aPatch, COORD *global, COORD *local)
{
  BndSegFuncPtr BndSegFunc;
  COORD loc[DIM_OF_BND], step[DIM_OF_BND], search_global[DIM];
  COORD min, dist;
  STD_PATCH *thePatch;

  thePatch = (STD_PATCH *)aPatch;
  BndSegFunc = thePatch->BndSegFunc;

#ifdef __TWODIM__

  /* first try midpoint */
  local[0] = 0.5*(thePatch->beta[0] + thePatch->alpha[0]);
  if ((*BndSegFunc)(thePatch->bs_data,local,search_global)) return (1);
  min = dist = (search_global[0]-global[0])*(search_global[0]-global[0]) +
               (search_global[1]-global[1])*(search_global[1]-global[1]);

  /* step trough the patch */
  if (thePatch->resolution <= 0) return (1);
  step[0] = (thePatch->beta[0] - thePatch->alpha[0])/(COORD)thePatch->resolution;
  for (loc[0]=thePatch->alpha[0]; loc[0]<=thePatch->beta[0]; loc[0] +=step[0])
  {
    if ((*BndSegFunc)(thePatch->bs_data,loc,search_global)) return (1);
    dist = (search_global[0]-global[0])*(search_global[0]-global[0]) +
           (search_global[1]-global[1])*(search_global[1]-global[1]);
    if (dist < min)
    {
      min = dist;
      local[0] = loc[0];
    }
  }

#endif

#ifdef __THREEDIM__

  /* first try midpoint */
  local[0] = 0.5*(thePatch->beta[0] + thePatch->alpha[0]);
  local[1] = 0.5*(thePatch->beta[1] + thePatch->alpha[1]);
  if ((*BndSegFunc)(thePatch->bs_data,local,search_global)) return (1);
  min = dist = (search_global[0]-global[0])*(search_global[0]-global[0]) +
               (search_global[1]-global[1])*(search_global[1]-global[1]) +
               (search_global[2]-global[2])*(search_global[2]-global[2]);

  /* step trough the patch */
  if (thePatch->resolution <= 0) return (1);
  step[0] = (thePatch->beta[0] - thePatch->alpha[0])/(COORD)thePatch->resolution;
  step[1] = (thePatch->beta[1] - thePatch->alpha[1])/(COORD)thePatch->resolution;
  for (loc[0]=thePatch->alpha[0]; loc[0]<=thePatch->beta[0]; loc[0] +=step[0])
    for (loc[1]=thePatch->alpha[1]; loc[1]<=thePatch->beta[1]; loc[1] +=step[1])
    {
      if ((*BndSegFunc)(thePatch->bs_data,loc,search_global)) return (1);
      dist = (search_global[0]-global[0])*(search_global[0]-global[0]) +
             (search_global[1]-global[1])*(search_global[1]-global[1]) +
             (search_global[2]-global[2])*(search_global[2]-global[2]);
      if (dist < min)
      {
        min = dist;
        local[0] = loc[0];
        local[1] = loc[1];
      }
    }

#endif

  return (0);
}

/****************************************************************************/
/*D
   Patch_local2global - get global from local coordinates on the patch

   SYNOPSIS:
   INT Patch_local2global (PATCH *thePatch, COORD *local, COORD *global);

   PARAMETERS:
   .  thePatch - the patch

   DESCRIPTION:
   This function get the point in physical space from local coordinates

   THE RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT Patch_local2global (PATCH *aPatch, COORD *local, COORD *global)
{
  BndSegFuncPtr BndSegFunc;
  STD_PATCH *thePatch;

  thePatch = (STD_PATCH *)aPatch;
  BndSegFunc = thePatch->BndSegFunc;
  return ((*BndSegFunc)(thePatch->bs_data,local,global));
}

/****************************************************************************/
/*D
   Patch_local2global - get boundary condition from local coordinates on the patch

   SYNOPSIS:
   INT Patch_local2bndcond (PATCH *thePatch, COORD *local, INT *type, DOUBLE *value);

   PARAMETERS:
   .  thePatch - the patch

   DESCRIPTION:
   This function get the boundary conditions from local coordinates

   THE RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT Patch_local2bndcond (PATCH *aPatch, COORD *local, DOUBLE *value, INT *type)
{
  BndCondProcPtr BndCond;
  STD_PATCH *thePatch;

  thePatch = (STD_PATCH *)aPatch;
  BndCond = thePatch->BndCond;
  return ((*BndCond)(thePatch->bc_data,NULL,local,value,type));
}

/****************************************************************************/
/*D
   Patch_GetPatchDesc - get patch description

   SYNOPSIS:
   INT Patch_GetPatchDesc (PATCH *thePatch, PATCH_DESC *thePatchDesc);

   PARAMETERS:
   .  thePatch - the patch

   DESCRIPTION:
   This function sets the patch descriptor according to a patch

   THE RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

INT Patch_GetPatchDesc (PATCH *aPatch, PATCH_DESC *thePatchDesc)
{
  INT i;
  STD_PATCH *thePatch;

  thePatch = (STD_PATCH *)aPatch;

  thePatchDesc->left              = thePatch->left;
  thePatchDesc->right             = thePatch->right;
  thePatchDesc->type              = thePatch->segType;
  thePatchDesc->n                 = thePatch->nPoints;
  thePatchDesc->id                = thePatch->id;
  thePatchDesc->resolution= thePatch->resolution;
  for (i=0; i<thePatchDesc->n; i++)
    thePatchDesc->pc[i].id = thePatch->points[i];

#ifdef __TWODIM__

  thePatchDesc->pc[0].lcoord[0] = thePatch->alpha[0];
  thePatchDesc->pc[1].lcoord[0] = thePatch->beta[0];

#endif

#ifdef __THREEDIM__

  thePatchDesc->pc[0].lcoord[0] = thePatch->alpha[0];
  thePatchDesc->pc[0].lcoord[1] = thePatch->alpha[1];

  thePatchDesc->pc[1].lcoord[0] = thePatch->beta[0];
  thePatchDesc->pc[1].lcoord[1] = thePatch->alpha[1];

  thePatchDesc->pc[2].lcoord[0] = thePatch->beta[0];
  thePatchDesc->pc[2].lcoord[1] = thePatch->beta[1];

  thePatchDesc->pc[3].lcoord[0] = thePatch->alpha[0];
  thePatchDesc->pc[3].lcoord[1] = thePatch->beta[1];

#endif

  return (0);
}

/****************************************************************************/
/*D
   InitStd_Domain - Create and initialize the std_domain

   SYNOPSIS:
   INT InitEnrol ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function creates the environment

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 when error occured.
   D*/
/****************************************************************************/

INT InitDom ()
{
  /* change to root directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitStd_Domain","could not changedir to root");
    return(__LINE__);
  }

  /* get env dir/var IDs for the problems */
  theProblemDirID  = GetNewEnvDirID();
  theBdryCondVarID = GetNewEnvVarID();

  /* install the /Domains directory */
  theDomainDirID = GetNewEnvDirID();
  if (MakeEnvItem("Domains",theProblemDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitEnrol","could not install '/Domains' dir");
    return(__LINE__);
  }
  theBdrySegVarID = GetNewEnvVarID();

  /* install the /BVP directory */
  theBVPDirID = GetNewEnvDirID();
  if (MakeEnvItem("BVP",theBVPDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitStd_Domain","could not install '/BVP' dir");
    return(__LINE__);
  }
  theTmpVarID = GetNewEnvVarID();
  thePatchVarID = GetNewEnvVarID();

  return (0);
}
