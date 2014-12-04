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

#include <config.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ugtypes.h"
#include "heaps.h"
#include "bio.h"
#include "ugdevices.h"
#include "domain.h"
#include "lgm_load.h"
#include "lgm_domain.h"
#include "misc.h"
#include "general.h"

#include "namespace.h"

USING_UG_NAMESPACE
  USING_UGDIM_NAMESPACE


#define BVP2LGM(p)                                                      ((LGM_DOMAIN*)(p))
#define BNDP2LGM(p)                                                     ((LGM_BNDP*)(p))
#define BNDS2LGM(p)                                                     ((LGM_BNDS*)(p))

static INT theBVPDirID;
static INT theLGMDomainVarID;
static INT theProblemDirID;
static INT theProblemVarID;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/**
* \brief  Create a new PROBLEM structure.
*
* This function allocates and initializes a new PROBLEM structure in the /LGM_PROBLEM directory.
*
* @param name  name of the problem.
* @param config  pointer to the configuration function.
* @param numOfCoefficients  number of coefficient functions.
* @param coeffs  pointer to coefficient functions.
* @param numOfUserFct  number of user coefficient functions.
* @param userfct  pointer to user coefficient functions.
*
* @return  pointer to new PROBLEM, NULL if out of memory.
   */
/****************************************************************************/

LGM_PROBLEM* NS_DIM_PREFIX CreateProblem (const char *name, InitProcPtr init, DomainSizeConfig domconfig, BndCondProcPtr BndCond, int numOfCoefficients, CoeffProcPtr coeffs[], int numOfUserFct, UserProcPtr userfct[])
{
  LGM_PROBLEM *newProblem;
  int i;

  if (ChangeEnvDir("/LGM_PROBLEM")==NULL) return(NULL);

  /* allocate new problem structure */
  newProblem = (LGM_PROBLEM *) MakeEnvItem (name,theProblemVarID,sizeof(LGM_PROBLEM)+(numOfCoefficients+numOfUserFct-1)*sizeof(void*));
  if (newProblem==NULL) return(NULL);

  /* fill in data */
  LGM_PROBLEM_INIT(newProblem)            = init;
  LGM_PROBLEM_CONFIG(newProblem)          = NULL;
  LGM_PROBLEM_DOMCONFIG(newProblem)       = domconfig;
  LGM_PROBLEM_BNDCOND(newProblem)         = BndCond;
  LGM_PROBLEM_INNERBNDCOND(newProblem)= NULL;
  LGM_PROBLEM_NCOEFF(newProblem)          = numOfCoefficients;
  LGM_PROBLEM_NUSERF(newProblem)          = numOfUserFct;
  for (i=0; i<numOfCoefficients; i++) LGM_PROBLEM_SETCOEFF(newProblem,i,coeffs[i]);
  for (i=0; i<numOfUserFct; i++) LGM_PROBLEM_SETUSERF(newProblem,i,userfct[i]);

  UserWrite("lgm_problem "); UserWrite(name); UserWrite(" installed\n");

  return (newProblem);
}

LGM_PROBLEM* NS_DIM_PREFIX CreateProblemWithInnerBCs (const char *name, InitProcPtr init, DomainSizeConfig domconfig, BndCondProcPtr BndCond, BndCondProcPtr InnerBndCond, int numOfCoefficients, CoeffProcPtr coeffs[], int numOfUserFct, UserProcPtr userfct[])
{
  LGM_PROBLEM *newProblem;
  int i;

  if (ChangeEnvDir("/LGM_PROBLEM")==NULL) return(NULL);

  /* allocate new problem structure */
  newProblem = (LGM_PROBLEM *) MakeEnvItem (name,theProblemVarID,sizeof(LGM_PROBLEM)+(numOfCoefficients+numOfUserFct-1)*sizeof(void*));
  if (newProblem==NULL) return(NULL);

  /* fill in data */
  LGM_PROBLEM_INIT(newProblem)            = init;
  LGM_PROBLEM_CONFIG(newProblem)          = NULL;
  LGM_PROBLEM_DOMCONFIG(newProblem)       = domconfig;
  LGM_PROBLEM_BNDCOND(newProblem)         = BndCond;
  LGM_PROBLEM_INNERBNDCOND(newProblem)= InnerBndCond;
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

static LGM_PROBLEM *Lgm_Problem_GetByName (const char *name)
{
  return((LGM_PROBLEM *) SearchEnv(name,"/LGM_PROBLEM",theProblemVarID,theProblemDirID));
}

void NS_DIM_PREFIX
Set_Current_BVP(BVP* theBVP)
{
  return;
}

/* domain interface function: for description see domain.h */
BVP *NS_DIM_PREFIX
BVP_GetByName (const char *name)
{
  return((BVP *) SearchEnv(name,"/LGM_BVP",theLGMDomainVarID,theBVPDirID));
}

/* domain interface function: for description see domain.h */
BVP* NS_DIM_PREFIX BVP_GetFirst (void)
{
  ENVDIR *theSBVPDir;
  BVP *theBVP;

  theSBVPDir = ChangeEnvDir("/LGM_BVP");
  assert (theSBVPDir!=NULL);
  theBVP = (BVP *) ENVDIR_DOWN(theSBVPDir);

  return (theBVP);
}

/* domain interface function: for description see domain.h */
BVP* NS_DIM_PREFIX BVP_GetNext (BVP *theBVP)
{
  if (theBVP==NULL) return (NULL);
  return ((BVP *) NEXT_ENVITEM(theBVP));
}

/* domain interface function: for description see domain.h */
BVP* NS_DIM_PREFIX BVP_Init (const char *name, HEAP *Heap, MESH *Mesh, INT MarkKey)
{
  LGM_DOMAIN *theDomain;
  LGM_PROBLEM *theProblem;
  BndCondProcPtr BndCond,InnerBndCond;
  INT i,nSubDom,conf_df_problem;
  char **argv;

  if ((theDomain = (LGM_DOMAIN *)BVP_GetByName(name))==NULL)
  {
    if ((theDomain = LGM_LoadDomain(name,name,Heap,theLGMDomainVarID,MarkKey))==NULL)
    {
      UserWriteF("ERROR in BVP_Init: cannot load domain '%s'\n",name);
      return (NULL);
    }

    /* set problem */
    theProblem = Lgm_Problem_GetByName(LGM_DOMAIN_PROBLEMNAME(theDomain));
    conf_df_problem = 0;
    if (theProblem==NULL)
    {
      theProblem = Lgm_Problem_GetByName("configurable");
      if (theProblem==NULL)
      {
        UserWrite("ERROR in BVP_Init: cannot find problem\n");
        return (NULL);
      }
      conf_df_problem = 1;
    }
    LGM_DOMAIN_PROBLEM(theDomain) = theProblem;

    /* initialize problem */
    if (conf_df_problem)
    {
      INT maxLineId = 0;

      if (theProblem->InitProblem==NULL) return (NULL);
      nSubDom = LGM_DOMAIN_NSUBDOM(theDomain);
      argv = (char **) GetTmpMem(Heap, sizeof(char *)*(nSubDom+1),MarkKey);
      if (argv==NULL)
      {
        UserWrite("ERROR in BVP_Init: cannot allocate argv\n");
        return (NULL);
      }
      for(i=1; i<=nSubDom; i++)
      {
        LGM_SUBDOMAIN *subdom = LGM_DOMAIN_SUBDOM(theDomain,i);

        argv[i] = LGM_SUBDOMAIN_UNIT(subdom);
      }

      /* get maximum Line/Surface-Id */
      maxLineId = GetMaximumSurfaceID(theDomain);

      if ((*(theProblem->InitProblem))(nSubDom, argv, maxLineId+1, LGM_DOMAIN_PROBLEMNAME(theDomain), Heap))
      {
        UserWrite("ERROR in BVP_Init: cannot initialize problem\n");
        return (NULL);
      }
    }

    /* set boundary conditions */
    BndCond = LGM_PROBLEM_BNDCOND(theProblem);
    InnerBndCond = LGM_PROBLEM_INNERBNDCOND(theProblem);
    if (SetBoundaryCondition(theDomain,BndCond,InnerBndCond)) return (NULL);
  }

  /* set bounding sphere */
  if (SetDomainSize(theDomain)) return (NULL);

  /* set mesh with nothing */
  if (Mesh!=NULL && LGM_LoadMesh(name, Heap,Mesh,theDomain,MarkKey))
  {
    Mesh->mesh_status   = MESHSTAT_EMPTY;
    Mesh->nBndP             = 0;
    Mesh->nInnP             = 0;
    Mesh->nSubDomains       = 0;
    Mesh->nbElements        = NULL;
    Mesh->nElements     = NULL;
    Mesh->VertexLevel   = NULL;
    Mesh->VertexPrio    = NULL;
    Mesh->ElementLevel  = NULL;
    Mesh->ElementPrio   = NULL;
    Mesh->ElemSideOnBnd = NULL;
  }

  /* allocate s2p table */
  LGM_DOMAIN_NPART(theDomain) = 1;
  LGM_DOMAIN_S2P_PTR(theDomain) = (INT*)GetFreelistMemory(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT));
  if (LGM_DOMAIN_S2P_PTR(theDomain)==NULL)
    return (NULL);
  /* HRR_TODO: fill number of parts */
  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
    LGM_DOMAIN_S2P(theDomain,i) = 0;
  theDomain->theHeap = Heap;

  return ((BVP *)theDomain);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BVP_Dispose (BVP *theBVP)
{
  LGM_DOMAIN *theDomain;

  theDomain = (LGM_DOMAIN *)theBVP;
  theDomain->v.locked = 0;
  if (ChangeEnvDir("/LGM_BVP")==NULL) return (1);
  if (RemoveEnvItem((ENVITEM *)theBVP)) return (1);

  return (0);
}

/* domain interface function: for description see domain.h */
BVP* NS_DIM_PREFIX BVP_Load (const char *name, INT argc, char **argv)
{
  return (NULL);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BVP_SetBVPDesc (BVP *aBVP, BVP_DESC *theBVPDesc)
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

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BVP_SetCoeffFct (BVP *aBVP, INT n, CoeffProcPtr *CoeffFct)
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

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BVP_SetUserFct (BVP *aBVP, INT n, UserProcPtr *UserFct)
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

/* domain interface function: for description see domain.h */
BNDP* NS_DIM_PREFIX BVP_InsertBndP (HEAP *Heap, BVP *aBVP, INT argc, char **argv)
{

        #ifdef __THREEDIM__
  double global[3];

  /* scan global coordinates */
  if (sscanf(argv[0],"ngbn %lf %lf %lf",global,global+1,global+2) != DIM)
  {
    UserWriteF("could not scan");

    return(NULL);
  }

  return(BNDP_InsertBndP(Heap,aBVP,global));
        #endif

  return (NULL);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size)
{
  return (1);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDS_Dispose (HEAP *Heap, BNDS *aBndS)
{
  LGM_BNDS *theBndS;

  if (aBndS == NULL) return(0);
  theBndS = BNDS2LGM(aBndS);
  return (PutFreelistMemory(Heap,theBndS,sizeof(theBndS)));
}
/* domain interface function: for description see domain.h */
#ifdef __TWODIM__
INT NS_DIM_PREFIX BNDP_Move (BNDP *aBndP, const DOUBLE global[])
{
  return(1);
}
#endif
#ifdef __THREEDIM__
INT NS_DIM_PREFIX BNDP_Move (BNDP *aBndP, const DOUBLE global[])
{
#ifdef NO_PROJECT
  INT i,j;
  LGM_BNDP *theBndP;

  theBndP=(LGM_BNDP*)aBndP;
  for (i=0; i<theBndP->nsurf; i++)
    for (j=0; j<DIM; j++)
      theBndP->Surf[i].global[j]=global[j];

  return (0);
#else
  return (1);
#endif
}
#endif

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDP_BndEDesc (BNDP *aBndP0, BNDP *aBndP1, INT *part)
{
  /* HRR_TODO: assign part */
  *part = 0;

  return(0);
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

INT NS_DIM_PREFIX InitDom ()
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
