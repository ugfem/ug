// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lgm_domain.c													*/
/*																			*/
/* Purpose:   source for lgm_domain                                                                             */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: klaus@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   08.07.96 begin												*/
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compiler.h"
#include "heaps.h"
#include "bio.h"
#include "devices.h"
#include "domain.h"
#include "lgm_load.h"
#include "lgm_domain.h"
#include "misc.h"
#include "general.h"

#define BVP2LGM(p)                                                      ((LGM_DOMAIN*)(p))
#define BNDP2LGM(p)                                                     ((LGM_BNDP*)(p))
#define BNDS2LGM(p)                                                     ((LGM_BNDS*)(p))

static INT theBVPDirID;
static INT theLGMDomainVarID;
static INT theProblemDirID;
static INT theProblemVarID;

/****************************************************************************/
/*D
   CreateProblem -  Create a new PROBLEM structure

   SYNOPSIS:
   PROBLEM *CreateProblem (char *name,
   ConfigProcPtr config, int numOfCoefficients, CoeffProcPtr coeffs[],
   int numOfUserFct, UserProcPtr userfct[]);

   PARAMETERS:
   .  name - name of the problem
   .  config - pointer to the configuration function
   .  numOfCoefficients - number of coefficient functions
   .  coeffs[] - pointer to coefficient functions
   .  numOfUserFct - number of user coefficient functions
   .  userfct[] - pointer to user coefficient functions

   DESCRIPTION:
   This function allocates and initializes a new PROBLEM structure in the /LGM_PROBLEM directory.

   RETURN VALUE:
   PROBLEM *
   .n    pointer to new PROBLEM
   .n    NULL if out of memory.
   D*/
/****************************************************************************/

LGM_PROBLEM *CreateProblem (char *name, ConfigProcPtr config, DomainSizeConfig domconfig, BndCondProcPtr BndCond, int numOfCoefficients, CoeffProcPtr coeffs[], int numOfUserFct, UserProcPtr userfct[])
{
  LGM_PROBLEM *newProblem;
  int i;

  if (ChangeEnvDir("/LGM_PROBLEM")==NULL) return(NULL);

  /* allocate new problem structure */
  newProblem = (LGM_PROBLEM *) MakeEnvItem (name,theProblemVarID,sizeof(LGM_PROBLEM)+(numOfCoefficients+numOfUserFct-1)*sizeof(void*));
  if (newProblem==NULL) return(NULL);

  /* fill in data */
  LGM_PROBLEM_CONFIG(newProblem)          = config;
  LGM_PROBLEM_DOMCONFIG(newProblem)       = domconfig;
  LGM_PROBLEM_BNDCOND(newProblem)         = BndCond;
  LGM_PROBLEM_NCOEFF(newProblem)          = numOfCoefficients;
  LGM_PROBLEM_NUSERF(newProblem)          = numOfUserFct;
  for (i=0; i<numOfCoefficients; i++) LGM_PROBLEM_SETCOEFF(newProblem,i,coeffs[i]);
  for (i=0; i<numOfUserFct; i++) LGM_PROBLEM_SETUSERF(newProblem,i,userfct[i]);

  UserWrite("lgm_problem "); UserWrite(name); UserWrite(" installed\n");

  return (newProblem);
}

/****************************************************************************/
/*
   Lgm_Problem_GetByName - get pointer to LGM_PROBLEM by name

   SYNOPSIS:
   LGM_PROBLEM *Lgm_Problem_GetByName (char *name);

   PARAMETERS:
   .  name - name

   DESCRIPTION:
   This function gives the pointer to the LGM_PROBLEM by its <name>.

   RETURN VALUE:
   LGM_PROBLEM *
   .n   pointer to LGM_PROBLEM
   .n   NULL if error.
 */
/****************************************************************************/

static LGM_PROBLEM *Lgm_Problem_GetByName (char *name)
{
  return((LGM_PROBLEM *) SearchEnv(name,"/LGM_PROBLEM",theProblemVarID,theProblemDirID));
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
  return((BVP *) SearchEnv(name,"/LGM_BVP",theLGMDomainVarID,theBVPDirID));
}

/****************************************************************************/
/*D
   BVP_GetFirst - Return a pointer to the first LGM_BVP

   SYNOPSIS:
   BVP *BVP_GetFirst (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function returns a pointer to the first BVP in the /LGM_BVP
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

  theSBVPDir = ChangeEnvDir("/LGM_BVP");
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
   This function returns a pointer to the next BVP in the /LGMD_BVP
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

/* domain interface function: for description see domain.h */
BVP *BVP_Init (char *name, HEAP *Heap, MESH *Mesh)
{
  LGM_DOMAIN *theDomain;
  LGM_PROBLEM *theProblem;
  BndCondProcPtr BndCond;
  INT i;

  if ((theDomain = (LGM_DOMAIN *)BVP_GetByName(name))==NULL)
  {
    if ((theDomain = LGM_LoadDomain(name,name,Heap,theLGMDomainVarID))==NULL)
    {
      UserWrite("ERROR in BVP_Init: cannot load domain\n");
      return (NULL);
    }

    /* set problem */
    theProblem = Lgm_Problem_GetByName(LGM_DOMAIN_PROBLEMNAME(theDomain));
    if (theProblem==NULL)
    {
      UserWrite("ERROR in BVP_Init: cannot find problem\n");
      return (NULL);
    }
    LGM_DOMAIN_PROBLEM(theDomain) = theProblem;

    /* set boundary conditions */
    BndCond = LGM_PROBLEM_BNDCOND(theProblem);
    if (SetBoundaryCondition(theDomain,BndCond)) return (NULL);
  }

  /* set bounding sphere */
  if (SetDomainSize(theDomain)) return (NULL);


  /* set mesh with nothing */
  Mesh->nBndP             = 0;
  Mesh->nInnP             = 0;
  Mesh->nSubDomains       = 0;
  Mesh->nbElements        = NULL;
  Mesh->nElements = NULL;
  Mesh->VertexLevel = NULL;
  Mesh->VertexPrio = NULL;
  Mesh->ElementLevel = NULL;
  Mesh->ElementPrio = NULL;

  /* allocate s2p table */
  LGM_DOMAIN_NPART(theDomain) = 1;
  LGM_DOMAIN_S2P_PTR(theDomain) = (INT*)GetFreelistMemory(Heap,LGM_DOMAIN_NSUBDOM(theDomain)*sizeof(INT));
  if (LGM_DOMAIN_S2P_PTR(theDomain)==NULL)
    return (NULL);
  /* HRR_TODO: fill number of parts */
  for (i=0; i<LGM_DOMAIN_NSUBDOM(theDomain); i++)
    LGM_DOMAIN_S2P(theDomain,i) = 0;

  return ((BVP *)theDomain);
}

/* domain interface function: for description see domain.h */
INT BVP_Dispose (BVP *theBVP)
{
  LGM_DOMAIN *theDomain;

  theDomain = (LGM_DOMAIN *)theBVP;
  theDomain->v.locked = 0;
  if (ChangeEnvDir("/LGM_BVP")==NULL) return (1);
  if (RemoveEnvItem((ENVITEM *)theBVP)) return (1);

  return (0);
}

/* domain interface function: for description see domain.h */
BVP *BVP_Load (char *name, INT argc, char **argv)
{
  return (NULL);
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
  LGM_DOMAIN *theDomain;
  LGM_PROBLEM *theProblem;
  INT i;

  /* cast */
  theDomain = BVP2LGM(aBVP);

  /* general part */
  strcpy(theBVPDesc->name,LGM_DOMAIN_NAME(theDomain));

  /* the domain part */
  for (i=0; i<DIM; i++)
    theBVPDesc->midpoint[i] = LGM_DOMAIN_MIDPOINT(theDomain)[i];
  BVPD_RADIUS(theBVPDesc)         = LGM_DOMAIN_RADIUS(theDomain);
  BVPD_CONVEX(theBVPDesc)         = LGM_DOMAIN_CONVEX(theDomain);
  BVPD_NSUBDOM(theBVPDesc)        = LGM_DOMAIN_NSUBDOM(theDomain);
  BVPD_NPARTS(theBVPDesc)         = LGM_DOMAIN_NPART(theDomain);
  BVPD_S2P_PTR(theBVPDesc)        = LGM_DOMAIN_S2P_PTR(theDomain);
  theProblem                                      = LGM_DOMAIN_PROBLEM(theDomain);
  if (theProblem==NULL) return (1);
  BVPD_NCOEFFF(theBVPDesc)        = LGM_PROBLEM_NCOEFF(theProblem);
  BVPD_NUSERF(theBVPDesc)         = LGM_PROBLEM_NUSERF(theProblem);
  BVPD_CONFIG(theBVPDesc)         = LGM_PROBLEM_CONFIG(theProblem);

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
  LGM_DOMAIN *theDomain;
  LGM_PROBLEM *theProblem;
  INT i;

  /* cast */
  theDomain       = BVP2LGM(aBVP);
  if (theDomain==NULL) return (1);
  theProblem      = LGM_DOMAIN_PROBLEM(theDomain);
  if (theProblem==NULL) return (1);

  /* check */
  if (n<-1 || n>=LGM_PROBLEM_NCOEFF(theProblem)) return (1);

  if (n==-1)
    for (i=0; i<LGM_PROBLEM_NCOEFF(theProblem); i++)
      CoeffFct[i] = LGM_PROBLEM_COEFF(theProblem,i);
  else
    CoeffFct[0] = LGM_PROBLEM_COEFF(theProblem,n);

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
  LGM_DOMAIN *theDomain;
  LGM_PROBLEM *theProblem;
  INT i;

  /* cast */
  theDomain       = BVP2LGM(aBVP);
  if (theDomain==NULL) return (1);
  theProblem      = LGM_DOMAIN_PROBLEM(theDomain);
  if (theProblem==NULL) return (1);

  /* check */
  if (n<-1 || n>=LGM_PROBLEM_NUSERF(theProblem)) return (1);

  if (n==-1)
    for (i=0; i<LGM_PROBLEM_NUSERF(theProblem); i++)
      UserFct[i] = LGM_PROBLEM_USERF(theProblem,i);
  else
    UserFct[0] = LGM_PROBLEM_USERF(theProblem,n);

  return (0);
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/* functions called by script commands										*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

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

BNDP *BVP_InsertBndP (HEAP *Heap, BVP *aBVP, INT argc, char **argv)
{
  return (NULL);
}

/* domain interface function: for description see domain.h */
INT BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size)
{
  return (1);
}

/* domain interface function: for description see domain.h */
INT BNDS_Dispose (HEAP *Heap, BNDS *aBndS)
{
  LGM_BNDS *theBndS;

  if (aBndS == NULL) return(0);
  theBndS = BNDS2LGM(aBndS);
  return (PutFreelistMemory(Heap,theBndS,sizeof(theBndS)));
}

/****************************************************************************/
/*D
   BNDP_BndEDesc - set part for boundary edge

   SYNOPSIS:
   INT BNDP_BndEDesc (BNDP *aBndP0, BNDP *aBndP1, INT *part);

   PARAMETERS:
   .  aBndP0, aBndP1 - endpoints of edge
   .  part - part-info for edge

   DESCRIPTION:
   This function set part of edge on boundary

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BNDP_BndEDesc (BNDP *aBndP0, BNDP *aBndP1, INT *part)
{
  return (1);
}

/****************************************************************************/
/*
   InitDom - init for lgm-domain

   SYNOPSIS:
   INT InitDom (void);

   PARAMETERS:

   DESCRIPTION:
   init the lgm-domain module

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.

   SEE ALSO:
 */
/****************************************************************************/

INT InitDom (void)
{
  /* change to root directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitLgm_Domain","could not changedir to root");
    return(__LINE__);
  }

  /* install the /LGM_BVP directory */
  theBVPDirID = GetNewEnvDirID();
  if (MakeEnvItem("LGM_BVP",theBVPDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitLgm_Domain","could not install '/LGM_BVP' dir");
    return(__LINE__);
  }
  theLGMDomainVarID = GetNewEnvVarID();

  /* change to root directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitLgm_Domain","could not changedir to root");
    return(__LINE__);
  }

  /* install the /LGM_PROBLEM directory */
  theProblemDirID = GetNewEnvDirID();
  if (MakeEnvItem("LGM_PROBLEM",theProblemDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitLgm_Domain","could not install '/LGM_PROBLEM' dir");
    return(__LINE__);
  }
  theProblemVarID = GetNewEnvVarID();

  /* init load procedures */
  if (InitLGMLoad()) return (1);

  return (0);
}
