// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  enrol.c														*/
/*																			*/
/* Purpose:   contains functions to enrol user defineable structures to         */
/*			  ug's environment.                                             */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   12.11.94 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment ugm
#endif

/****************************************************************************/
/*																			*/
/*		defines to exclude functions										*/
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

/* devices module */
#include "devices.h"

/* grid manager module */
#include "switch.h"
#include "gm.h"
#include "algebra.h"
#include "misc.h"
#include "enrol.h"

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

static INT theFormatDirID;                      /* env type for Format dir				*/
static INT theSymbolVarID;                      /* env type for Format vars                     */

static INT theProblemDirID;             /* env type for Problem dir                     */
static INT theBdryCondVarID;            /* env type for Problem vars			*/

static INT theDomainDirID;                      /* env type for Domain dir				*/
static INT theBdrySegVarID;             /* env type for bdry segment vars		*/

/* data for CVS */
static char rcsid[] = "$Header$";


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
/*
   CreateFormat	- Create a new FORMAT structure in the environment

   SYNOPSIS:
   FORMAT *CreateFormat (char *name,int sVertex,int sNode,int sDiag,
   int sElement,int sLink,int sEdge, int sMultiGrid,
                                          ConversionProcPtr SaveVertex,
                                          ConversionProcPtr SaveNode,
                                          ConversionProcPtr SaveDiag,
                                          ConversionProcPtr SaveElement,
                                          ConversionProcPtr SaveLink,
                                          ConversionProcPtr SaveEdge,
                                          ConversionProcPtr SaveGrid,
                                          ConversionProcPtr SaveMultiGrid,
                                          ConversionProcPtr LoadVertex,
                                          ConversionProcPtr LoadNode,
                                          ConversionProcPtr LoadDiag,
                                          ConversionProcPtr LoadElement,
                                          ConversionProcPtr LoadLink,
                                          ConversionProcPtr LoadEdge,
                                          ConversionProcPtr LoadGrid,
                                          ConversionProcPtr LoadMultiGrid,
                                          ConversionProcPtr PrintVertex,
                                          ConversionProcPtr PrintNode,
                                          ConversionProcPtr PrintDiag,
                                          ConversionProcPtr PrintElement,
                                          ConversionProcPtr PrintLink,
                                          ConversionProcPtr PrintEdge,
                                          ConversionProcPtr PrintGrid,
                                          ConversionProcPtr PrintMultiGrid);

   PARAMETERS:
   .  name - name of new format structure
   .  sVertex -
   .  sNode -
   .  sDiag -
   .  sElement -
   .  sLink -
   .  sEdge -
   .  sMultigrid -
   .  SaveVertex -
   .  SaveNode -
   .  SaveDiag -
   .  SaveElement -
   .  SaveLink -
   .  SaveEdge -
   .  SaveGrid -
   .  SaveMultiGrid -
   .  ConversionProcPtr LoadVertex -
   .  ConversionProcPtr LoadNode -
   .  ConversionProcPtr LoadDiag -
   .  ConversionProcPtr LoadElement -
   .  ConversionProcPtr LoadLink -
   .  ConversionProcPtr LoadEdge -
   .  ConversionProcPtr LoadGrid -
   .  ConversionProcPtr LoadMultiGrid -
   .  ConversionProcPtr PrintVertex -
   .  ConversionProcPtr PrintNode -
   .  ConversionProcPtr PrintDiag -
   .  ConversionProcPtr PrintElement -
   .  ConversionProcPtr PrintLink -
   .  ConversionProcPtr PrintEdge -
   .  ConversionProcPtr PrintGrid -
   .  ConversionProcPtr PrintMultiGrid -

   DESCRIPTION:
   This function allocates and initializes a new FORMAT structure in the environment.

   RETURN VALUE:
   FORMAT *
   .n     pointer to
   .n     NULL if out of memory.
 */
/****************************************************************************/

#ifdef __version23__
FORMAT *CreateFormat (char *name,int sVertex,int sNode,int sDiag,int sElement,int sLink,int sEdge, int sMultiGrid,
                      ConversionProcPtr SaveVertex,
                      ConversionProcPtr SaveNode,
                      ConversionProcPtr SaveDiag,
                      ConversionProcPtr SaveElement,
                      ConversionProcPtr SaveLink,
                      ConversionProcPtr SaveEdge,
                      ConversionProcPtr SaveGrid,
                      ConversionProcPtr SaveMultiGrid,
                      ConversionProcPtr LoadVertex,
                      ConversionProcPtr LoadNode,
                      ConversionProcPtr LoadDiag,
                      ConversionProcPtr LoadElement,
                      ConversionProcPtr LoadLink,
                      ConversionProcPtr LoadEdge,
                      ConversionProcPtr LoadGrid,
                      ConversionProcPtr LoadMultiGrid,
                      ConversionProcPtr PrintVertex,
                      ConversionProcPtr PrintNode,
                      ConversionProcPtr PrintDiag,
                      ConversionProcPtr PrintElement,
                      ConversionProcPtr PrintLink,
                      ConversionProcPtr PrintEdge,
                      ConversionProcPtr PrintGrid,
                      ConversionProcPtr PrintMultiGrid)
{
  FORMAT *newFormat;

  /* change to /Formats directory */
  if (ChangeEnvDir("/Formats")==NULL)
    return(NULL);

  /* allocate new format structure */
  newFormat = (FORMAT *) MakeEnvItem (name,theFormatDirID,sizeof(FORMAT));
  if (newFormat==NULL) return(NULL);

  /* fill in data */
  newFormat->sVertex      = sVertex;
  newFormat->sNode        = sNode;
  newFormat->sDiag        = sDiag;
  newFormat->sElement = sElement;
  newFormat->sLink        = sLink;
  newFormat->sEdge        = sEdge;
  newFormat->sMultiGrid = sMultiGrid;

  newFormat->SaveVertex  = SaveVertex;
  newFormat->SaveNode    = SaveNode;
  newFormat->SaveDiag    = SaveDiag;
  newFormat->SaveElement = SaveElement;
  newFormat->SaveLink    = SaveLink;
  newFormat->SaveEdge    = SaveEdge;
  newFormat->SaveGrid    = SaveGrid;
  newFormat->SaveMultiGrid = SaveMultiGrid;

  newFormat->LoadVertex  = LoadVertex;
  newFormat->LoadNode    = LoadNode;
  newFormat->LoadDiag    = LoadDiag;
  newFormat->LoadElement = LoadElement;
  newFormat->LoadLink    = LoadLink;
  newFormat->LoadEdge    = LoadEdge;
  newFormat->LoadGrid    = LoadGrid;
  newFormat->LoadMultiGrid = LoadMultiGrid;

  newFormat->PrintVertex  = PrintVertex;
  newFormat->PrintNode    = PrintNode;
  newFormat->PrintDiag    = PrintDiag;
  newFormat->PrintElement = PrintElement;
  newFormat->PrintLink    = PrintLink;
  newFormat->PrintEdge    = PrintEdge;
  newFormat->PrintGrid    = PrintGrid;
  newFormat->PrintMultiGrid = PrintMultiGrid;

  if (ChangeEnvDir(name)==NULL) return(NULL);
  UserWrite("format "); UserWrite(name); UserWrite(" installed\n");

  return(newFormat);
}
#endif


#ifdef __version3__
FORMAT *Ugly_CreateFormat (char *name,
                           INT sVertex, INT sMultiGrid,
                           INT *VectorSizes,
                           INT *FromType, INT *ToType, INT *MatrixSizes, INT *ConnectionDepth,

                           ConversionProcPtr PrintVertex,
                           ConversionProcPtr PrintGrid,
                           ConversionProcPtr PrintMultigrid,
                           ConversionProcPtr PrintVector[MAXVECTORS],
                           ConversionProcPtr PrintMatrix[MAXVECTORS][MAXVECTORS] )
{
  FORMAT *newFormat;
  INT i, j, MaxDepth, NeighborhoodDepth;

  /* check if format definition matches switch defines */
  if ((NODE_DATA&&VectorSizes[NODEVECTOR]<=0) || (!NODE_DATA&&VectorSizes[NODEVECTOR]>0)) return (NULL);
  if ((EDGE_DATA&&VectorSizes[EDGEVECTOR]<=0) || (!EDGE_DATA&&VectorSizes[EDGEVECTOR]>0)) return (NULL);
#ifdef __THREEDIM__
  if ((SIDE_DATA&&VectorSizes[SIDEVECTOR]<=0) || (!SIDE_DATA&&VectorSizes[SIDEVECTOR]>0)) return (NULL);
#endif
  if ((ELEM_DATA&&VectorSizes[ELEMVECTOR]<=0) || (!ELEM_DATA&&VectorSizes[ELEMVECTOR]>0)) return (NULL);

  /* change to /Formats directory */
  if (ChangeEnvDir("/Formats")==NULL)
    return (NULL);

  /* allocate new format structure */
  newFormat = (FORMAT *) MakeEnvItem (name,theFormatDirID,sizeof(FORMAT));
  if (newFormat==NULL) return(NULL);

  /* fill in data */
  newFormat->sVertex         = sVertex;
  newFormat->sMultiGrid  = sMultiGrid;

  /* set vector sizes */
  newFormat->VectorSizes[NODEVECTOR] = VectorSizes[NODEVECTOR];
  newFormat->VectorSizes[EDGEVECTOR] = VectorSizes[EDGEVECTOR];
#ifdef __THREEDIM__
  newFormat->VectorSizes[SIDEVECTOR] = VectorSizes[SIDEVECTOR];
#endif
  newFormat->VectorSizes[ELEMVECTOR] = VectorSizes[ELEMVECTOR];

  /* set matrix sizes and matrix depths */
  MaxDepth = NeighborhoodDepth = 0;
  for (i=0; i<MAXMATRICES; i++)
  {
    if (newFormat->VectorSizes[FromType[i]]>0 && newFormat->VectorSizes[ToType[i]]>0 && MatrixSizes[i]>0 && ConnectionDepth[i]>=0)
    {
      newFormat->MatrixSizes[MatrixType[FromType[i]][ToType[i]]] = MatrixSizes[i];
      newFormat->ConnectionDepth[MatrixType[FromType[i]][ToType[i]]] = ConnectionDepth[i];
    }
    else
    {
      newFormat->MatrixSizes[MatrixType[FromType[i]][ToType[i]]] = 0;
      newFormat->ConnectionDepth[MatrixType[FromType[i]][ToType[i]]] = 0;
    }
    MaxDepth = MAX(MaxDepth,ConnectionDepth[i]);
    if (FromType[i] == ELEMVECTOR && ToType[i] == ELEMVECTOR)
      NeighborhoodDepth = MAX(NeighborhoodDepth,ConnectionDepth[i]);
    else
      NeighborhoodDepth = MAX(NeighborhoodDepth,ConnectionDepth[i]+1);
  }
  newFormat->MaxConnectionDepth = MaxDepth;
  newFormat->NeighborhoodDepth  = NeighborhoodDepth;

  newFormat->PrintVertex          = PrintVertex;
  newFormat->PrintGrid            = PrintGrid;
  newFormat->PrintMultigrid               = PrintMultigrid;
  for (i=0; i<MAXVECTORS; i++)
    newFormat->PrintVector[i] = PrintVector[i];
  for (i=0; i<MAXVECTORS; i++)
    for (j=0; j<MAXVECTORS; j++)
      newFormat->PrintMatrix[i][j] = PrintMatrix[i][j];

  if (ChangeEnvDir(name)==NULL) return(NULL);
  UserWrite("format "); UserWrite(name); UserWrite(" installed\n");

  return(newFormat);
}


/****************************************************************************/
/*D
   CreateFormat	- Create a new FORMAT structure in the environment

   SYNOPSIS:
   FORMAT *CreateFormat (char *name, INT sVertex, INT sMultiGrid,
                ConversionProcPtr PrintVertex, ConversionProcPtr PrintGrid,
                ConversionProcPtr PrintMultigrid, INT nvDesc, VectorDescriptor *vDesc,
                INT nmDesc, MatrixDescriptor *mDesc);

   PARAMETERS:
   .  name - name of new format structure
   .  sVertex - size of user data space in VERTEX counted in bytes
   .  sMultiGrid -  size of user data space in MULTIGRID counted in bytes
   .  PrintVertex - pointer to conversion procedure
   .  PrintGrid - pointer to conversion procedure
   .  PrintMultigrid - pointer to conversion procedure
   .  nvDesc - number of vector descriptors
   .  vDesc - pointer to vector descriptor
   .  nmDesc - number of matrix desciptors
   .  mDesc - pointer to matrix descriptor

   DESCRIPTION:
   This function allocates and initializes a new FORMAT structure in the environment.
   The parameters vDesc and mDesc are pointers to structures which describe the
   VECTOR or MATRIX types used.
   VectorDescriptor is defined as
   .vb
          typedef struct {
          int pos;
          int size;
          ConversionProcPtr print;
          } VectorDescriptor ;
   .ve
        The components have the following meaning

   .   pos - this is the position for which this description is valid, one of
                        NODEVECTOR, EDGEVECTOR, ELEMVECTOR, SIDEVECTOR
   .   size - the data size of a VECTOR structure on this position in bytes
   .   print - pointer to a function which is called for printing the contents of the data
                        fields.

        MatrixDescriptor has the definition
   .vb
           typedef struct {
           int from;
           int to;
           int size;
           int depth;
           ConversionProcPtr print;
           } MatrixDescriptor ;
   .ve
        The meaning of the components is

   .   from - this connection goes from position
   .   to - to this position, position is one of NODEVECTOR, EDGEVECTOR, ELEMVECTOR, SIDEVECTOR
   .   size - this defines the size in bytes per connection
   .   depth - this connection has the depth defined here
   .   print - function to print the data.

   EXAMPLES:
   A small example to create a format looks like the following. In this format only
   vectors in nodes are used and therfore all connections connect two nodevectors.
   .vb
      // we need dofs only in nodes
      vd[0].pos   = NODEVECTOR;
      vd[0].size  = 3*sizeof(DOUBLE);
      vd[0].print = Print_3_NodeVectorData;

      // and the following connection: node-node
      md[0].from  = NODEVECTOR;
      md[0].to    = NODEVECTOR;
      md[0].size  = sizeof(DOUBLE);
      md[0].depth = 0;
      md[0].print = Print_1_NodeNodeMatrixData;

      newFormat = CreateFormat("full scalar",0,0,
                  (ConversionProcPtr)NULL,(ConversionProcPtr)NULL,
                  (ConversionProcPtr)NULL,1,vd,1,md);
   .ve

   RETURN VALUE:
   FORMAT *
   .n     pointer to FORMAT
   .n     NULL if out of memory.
   D*/
/****************************************************************************/

FORMAT *CreateFormat (char *name, INT sVertex, INT sMultiGrid,
                      ConversionProcPtr PrintVertex, ConversionProcPtr PrintGrid, ConversionProcPtr PrintMultigrid,
                      INT nvDesc, VectorDescriptor *vDesc, INT nmDesc, MatrixDescriptor *mDesc)
{
  FORMAT *newFormat;
  INT i, j, MaxDepth, NeighborhoodDepth;


  /* change to /Formats directory */
  if (ChangeEnvDir("/Formats")==NULL)
    return (NULL);

  /* allocate new format structure */
  newFormat = (FORMAT *) MakeEnvItem (name,theFormatDirID,sizeof(FORMAT));
  if (newFormat==NULL) return(NULL);

  /* fill in data */
  newFormat->sVertex                      = sVertex;
  newFormat->sMultiGrid           = sMultiGrid;
  newFormat->PrintVertex          = PrintVertex;
  newFormat->PrintGrid            = PrintGrid;
  newFormat->PrintMultigrid       = PrintMultigrid;


  /* initialize with zero */
  for (i=0; i<MAXVECTORS; i++)
  {
    newFormat->VectorSizes[i] = 0;
    newFormat->PrintVector[i] = (ConversionProcPtr) NULL;
    for (j=0; j<MAXVECTORS; j++) newFormat->PrintMatrix[i][j] = (ConversionProcPtr) NULL;
  }
  for (i=0; i<MAXMATRICES; i++)
  {
    newFormat->MatrixSizes[i] = 0;
    newFormat->ConnectionDepth[i] = 0;
  }
  MaxDepth = NeighborhoodDepth = 0;

  /* set vector stuff */
  for (i=0; i<nvDesc; i++)
  {
    if ((vDesc[i].pos<0)||(vDesc[i].pos>=MAXVECTORS)||(vDesc[i].size<0)) return(NULL);
    newFormat->VectorSizes[vDesc[i].pos] = vDesc[i].size;
    newFormat->PrintVector[vDesc[i].pos] = vDesc[i].print;
  }

  /* set connection stuff */
  for (i=0; i<nmDesc; i++)
  {
    if ((mDesc[i].from<0)||(mDesc[i].from>=MAXVECTORS)) return(NULL);
    if ((mDesc[i].to<0)  ||(mDesc[i].to>=MAXVECTORS)) return(NULL);
    if ((mDesc[i].size<0)||(mDesc[i].depth<0)) return(NULL);
    if (newFormat->VectorSizes[mDesc[i].from]>0 &&
        newFormat->VectorSizes[mDesc[i].to]>0 &&
        mDesc[i].size>0 && mDesc[i].depth>=0)
    {
      newFormat->MatrixSizes[MatrixType[mDesc[i].from][mDesc[i].to]] = mDesc[i].size;
      newFormat->ConnectionDepth[MatrixType[mDesc[i].from][mDesc[i].to]] = mDesc[i].depth;
      MaxDepth = MAX(MaxDepth,mDesc[i].depth);
      if ((mDesc[i].from==ELEMVECTOR)&&(mDesc[i].to==ELEMVECTOR))
        NeighborhoodDepth = MAX(NeighborhoodDepth,mDesc[i].depth);
      else
        NeighborhoodDepth = MAX(NeighborhoodDepth,mDesc[i].depth+1);
      newFormat->PrintMatrix[mDesc[i].from][mDesc[i].to] = mDesc[i].print;
      newFormat->PrintMatrix[mDesc[i].to][mDesc[i].from] = mDesc[i].print;
    }
  }
  newFormat->MaxConnectionDepth = MaxDepth;
  newFormat->NeighborhoodDepth  = NeighborhoodDepth;

  if (ChangeEnvDir(name)==NULL) return(NULL);
  UserWrite("format "); UserWrite(name); UserWrite(" installed\n");

  return(newFormat);
}
#endif

/****************************************************/
/*D
   GetFormat - Get a format pointer from the environment

   SYNOPSIS:
   FORMAT *GetFormat (char *name)

   PARAMETERS:
   .  name - name of the format

   DESCRIPTION:
   This function searches the directory /Formats for a format

   RETURN VALUE:
   FORMAT *
   .n   pointer to FORMAT
   .n   NULL  if not found or error.
   D*/
/****************************************************/

FORMAT *GetFormat (const char *name)
{
  return((FORMAT *) SearchEnv(name,"/Formats",theFormatDirID,theFormatDirID));
}

/****************************************************/
/*D
   ChangeToFormatDir - change to format directory with name

   SYNOPSIS:
   INT ChangeToFormatDir (const char *name)

   PARAMETERS:
   .  name - name of the format

   DESCRIPTION:
   This function changes to the format directory with name.

   RETURN VALUE:
   INT
   .n   0: ok
   .n   1: could not change to /Formats dir
   .n   1: could not change to /Formats/<name> dir
   D*/
/****************************************************/

INT ChangeToFormatDir (const char *name)
{
  if (ChangeEnvDir("/Formats")==NULL)
    return (1);
  if (ChangeEnvDir(name)==NULL)
    return (2);

  return (0);
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
   InitEnrol - Create and initialize the environment

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

INT InitEnrol ()
{
  /* install the /Formats directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitEnrol","could not changedir to root");
    return(__LINE__);
  }
  theFormatDirID = GetNewEnvDirID();
  if (MakeEnvItem("Formats",theFormatDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitEnrol","could not install '/Formats' dir");
    return(__LINE__);
  }
  theSymbolVarID = GetNewEnvVarID();

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

  return (GM_OK);
}
