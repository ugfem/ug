// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  evalproc.c													*/
/*																			*/
/* Purpose:   evaluation functions											*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de			                        */
/*																			*/
/* History:   31.03.92 begin, ug version 2.0								*/
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"

#include "gm.h"
#include "ugenv.h"
#include "evm.h"
#include "ugdevices.h"
#include "shapes.h"
#include "elements.h"

USING_UG_NAMESPACES

/****************************************************************************/
/*                                                                          */
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAX_COEFFPROC_ELEMEVAL                  50

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct Coubling_CoeffProc_Name {

  INT nUsed;
  char ElemEvalName[MAX_COEFFPROC_ELEMEVAL][NAMESIZE];
  CoeffProcPtr theCoeffProc[MAX_COEFFPROC_ELEMEVAL];
  CoeffProcPtr theActualCoeffProc;
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

/* environment directory and item IDs used in this source file */
static INT theEEvalProcDirID;
static INT theElemValVarID;

static INT theMEvalProcDirID;
static INT theMatrixValVarID;

static INT theVEvalProcDirID;
static INT theElemVectorVarID;

/* variables used for CoeffProcElemEval */
static struct Coubling_CoeffProc_Name Couple_for_ElemValue;
static struct Coubling_CoeffProc_Name Couple_for_ElemVector;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*
   CreateElementValueEvalProc - Create element value plot proceedure in evironement

   SYNOPSIS:
   EVALUES *CreateElementValueEvalProc (const char *name,
   PreprocessingProcPtr PreProc, ElementEvalProcPtr EvalProc);

   PARAMETERS:
   .  name -
   .  PreProc -
   .  EvalProc -

   DESCRIPTION:
   This function creates element value plot proceedure in evironement.

   RETURN VALUE:
   EVALUES *
   .n     pointer to
   .n     Null if error occured.
 */
/****************************************************************************/

EVALUES * NS_DIM_PREFIX CreateElementValueEvalProc (const char *name, PreprocessingProcPtr PreProc, ElementEvalProcPtr EvalProc)
{
  EVALUES *newElementValues;

  /* change to directory */
  if (ChangeEnvDir("/ElementEvalProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newElementValues = (EVALUES*) MakeEnvItem (name,theElemValVarID,sizeof(EVALUES));
  if (newElementValues==NULL) return(NULL);

  /* fill in data */
  newElementValues->PreprocessProc = PreProc;
  newElementValues->EvalProc = EvalProc;

  UserWrite("ElementValueEvalProc "); UserWrite(name); UserWrite(" installed\n");

  return(newElementValues);
}

/****************************************************************************/
/*
   GetElementValueEvalProc - Get element value plot proceedure in evironement from name

   SYNOPSIS:
   EVALUES *GetElementValueEvalProc (const char *name);

   PARAMETERS:
   .  name -

   DESCRIPTION:
   This function gets element value plot proceedure in evironement from name.

   RETURN VALUE:
   EVALUES *
   .n     pointer to
   .n     NULL if error occured.
 */
/****************************************************************************/

EVALUES * NS_DIM_PREFIX GetElementValueEvalProc (const char *name)
{
  if (ChangeEnvDir("/ElementEvalProcs")==NULL) return(NULL);
  return((EVALUES*) SearchEnv(name,".",theElemValVarID,SEARCHALL));
}

/****************************************************************************/
/*
   GetFirstElementValueEvalProc - Get the first element eval proc

   SYNOPSIS:
   EVALUES *GetFirstElementValueEvalProc (void);

   PARAMETERS:

   DESCRIPTION:
   This function gets the first element eval proc.

   RETURN VALUE:
   EVALUES *
   .n      pointer to EVALUES
   .n      NULL if there is no
 */
/****************************************************************************/

EVALUES * NS_DIM_PREFIX GetFirstElementValueEvalProc (void)
{
  ENVITEM *theEvalProc;

  if ((theEvalProc=(ENVITEM*)ChangeEnvDir("/ElementEvalProcs")) == NULL) return (NULL);

  for (theEvalProc=ENVITEM_DOWN(theEvalProc); theEvalProc!=NULL; theEvalProc=NEXT_ENVITEM(theEvalProc))
    if (ENVITEM_TYPE(theEvalProc) == theElemValVarID)
      return ((EVALUES*)theEvalProc);
  return (NULL);
}

/****************************************************************************/
/*
   GetNextElementValueEvalProc - Get the next element eval proc

   SYNOPSIS:
   EVALUES *GetNextElementValueEvalProc (EVALUES *EvalProc);

   PARAMETERS:
   .  EvalProc - evaluation procedure

   DESCRIPTION:
   This function gets the next 'EVALUES' in the element-eval-proc-list

   RETURN VALUE:
   EVALUES *
   .n      pointer to EVALUES *
   .n      NULL if there is no
 */
/****************************************************************************/

EVALUES * NS_DIM_PREFIX GetNextElementValueEvalProc (EVALUES *EvalProc)
{
  ENVITEM *theEvalProc;

  for (theEvalProc=EvalProc->v.next; theEvalProc!=NULL; theEvalProc=theEvalProc->v.next )
    if (theEvalProc->d.type == theElemValVarID)
      return ((EVALUES*)theEvalProc);
  return (NULL);
}

/****************************************************************************/
/*
   CreateElementVectorEvalProc - Create element vector plot proceedure in evironement

   SYNOPSIS:
   EVECTOR *CreateElementVectorEvalProc (const char *name,
   PreprocessingProcPtr PreProc, ElementVectorProcPtr EvalProc, INT d);

   PARAMETERS:
   .  name -
   .  PreProc -
   .  EvalProc -
   .  d - dimension of vector

   DESCRIPTION:
   This function creates element vector plot proceedure in evironement.

   RETURN VALUE:
   EVECTOR *
   .n      pointer to EVECTOR *
   .n      Null if not found.
 */
/****************************************************************************/

EVECTOR * NS_DIM_PREFIX CreateElementVectorEvalProc (const char *name, PreprocessingProcPtr PreProc, ElementVectorProcPtr EvalProc, INT d)
{
  EVECTOR *newElemVector;

  /* change directory */
  if (ChangeEnvDir("/ElementVectorEvalProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newElemVector = (EVECTOR*) MakeEnvItem (name,theElemVectorVarID,sizeof(EVECTOR));
  if (newElemVector==NULL) return(NULL);

  /* fill in data */
  newElemVector->PreprocessProc = PreProc;
  newElemVector->EvalProc = EvalProc;
  newElemVector->dimension = d;

  UserWrite("ElementVectorEvalProc "); UserWrite(name); UserWrite(" installed\n");
  return(newElemVector);
}

/****************************************************************************/
/*
   GetElementVectorEvalProc -

   SYNOPSIS:
   EVECTOR *GetElementVectorEvalProc (const char *name);

   PARAMETERS:
   .  name -

   DESCRIPTION:
   This function gets element vector plot proceedure in evironement from name.

   RETURN VALUE:
   EVECTOR *
   .n      pointer to
   .n      NULL if not found.
 */
/****************************************************************************/

EVECTOR * NS_DIM_PREFIX GetElementVectorEvalProc (const char *name)
{
  if (ChangeEnvDir("/ElementVectorEvalProcs")==NULL) return(NULL);
  return((EVECTOR*) SearchEnv(name,".",theElemVectorVarID,SEARCHALL));
}

/****************************************************************************/
/*
   GetFirstElementVectorEvalProc - Get the first element evec proc

   SYNOPSIS:
   EVECTOR *GetFirstElementVectorEvalProc (void);

   PARAMETERS:

   DESCRIPTION:
   This function gets the first element evec proc.

   RETURN VALUE:
   EVECTOR *
   .n      pointer to EVECTOR
   .n      NULL if there is no
 */
/****************************************************************************/

EVECTOR * NS_DIM_PREFIX GetFirstElementVectorEvalProc (void)
{
  ENVITEM *theEvecProc;

  if ((theEvecProc=(ENVITEM*)ChangeEnvDir("/ElementVectorEvalProcs")) == NULL) return (NULL);

  for (theEvecProc=ENVITEM_DOWN(theEvecProc); theEvecProc!=NULL; theEvecProc=NEXT_ENVITEM(theEvecProc))
    if (ENVITEM_TYPE(theEvecProc) == theElemVectorVarID)
      return ((EVECTOR*)theEvecProc);
  return (NULL);
}

/****************************************************************************/
/*D
   GetNextElementValueEvalProc - Get the next element evec proc

   SYNOPSIS:
   EVALUES *GetNextElementVectorEvalProc (EVALUES *EvecProc);

   PARAMETERS:
   .  EvecProc -

   DESCRIPTION:
   This function gets the next 'EVALUES' in the element-evec-proc-list

   RETURN VALUE:
   EVECTOR *
   .n      pointer to EVECTOR *
   .n      NULL if there is no
   D*/
/****************************************************************************/

EVECTOR * NS_DIM_PREFIX GetNextElementVectorEvalProc (EVECTOR *EvecProc)
{
  ENVITEM *theEvecProc;

  for (theEvecProc=EvecProc->v.next; theEvecProc!=NULL; theEvecProc=theEvecProc->v.next )
    if (theEvecProc->d.type == theElemVectorVarID)
      return ((EVECTOR*)theEvecProc);
  return (NULL);
}

/****************************************************************************/
/*
   CreateMatrixValueEvalProc - Create matrix value plot proceedure in evironement

   SYNOPSIS:
   MVALUES *CreateMatrixValueEvalProc (const char *name, PreprocessingProcPtr PreProc,
   MatrixEvalProcPtr EvalProc);

   PARAMETERS:
   .  name -
   .  PreProc -
   .  EvalProc -

   DESCRIPTION:
   This function create matrix value plot proceedure in evironement.

   RETURN VALUE:
   MVALUES *
   .n      pointer to
   .n      NULL, if error occured.
 */
/****************************************************************************/

MVALUES * NS_DIM_PREFIX CreateMatrixValueEvalProc (const char *name, PreprocessingProcPtr PreProc, MatrixEvalProcPtr EvalProc)
{
  MVALUES *newMatrixValues;

  /* change to directory */
  if (ChangeEnvDir("/MatrixEvalProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newMatrixValues = (MVALUES*) MakeEnvItem (name,theMatrixValVarID,sizeof(MVALUES));
  if (newMatrixValues==NULL) return(NULL);

  /* fill in data */
  newMatrixValues->PreprocessProc = PreProc;
  newMatrixValues->EvalProc = EvalProc;

  UserWrite("MatrixValueEvalProc "); UserWrite(name); UserWrite(" installed\n");
  return(newMatrixValues);
}

/****************************************************************************/
/*
   GetMatrixValueEvalProc - Get matrix value plot proceedure in evironement from name

   SYNOPSIS:
   MVALUES *GetMatrixValueEvalProc (const char *name);

   PARAMETERS:
   .  name -

   DESCRIPTION:
   This function gets matrix value plot proceedure in evironement from name.

   RETURN VALUE:
   MVALUES *
   .n      pointer to
   .n      NULL if error occured.
 */
/****************************************************************************/

MVALUES * NS_DIM_PREFIX GetMatrixValueEvalProc (const char *name)
{
  if (ChangeEnvDir("/MatrixEvalProcs")==NULL) return(NULL);
  return((MVALUES*) SearchEnv(name,".",theMatrixValVarID,SEARCHALL));
}


/****************************************************************************/
/*
   CreateElementValueEvalProcFromCoeffProc - Create element value plot proceedure from CoeffProc

   SYNOPSIS:
   EVALUES *CreateElementValueEvalProcFromCoeffProc (const char *name,
   CoeffProcPtr CoeffProc)

   PARAMETERS:
   .  name -
   .  CoeffProc -

   DESCRIPTION:
   This function creates element value plot proceedure in evironement from a
   .n CoeffProc

   RETURN VALUE:
   EVALUES *
   .n     pointer to
   .n     Null if error occured.
 */
/****************************************************************************/

static INT CoeffProcElementValuePreProc (const char *name, MULTIGRID *theMG)
{
  INT i;

  for (i=0; i<Couple_for_ElemValue.nUsed; i++)
    if (strcmp(Couple_for_ElemValue.ElemEvalName[i],name)==0)
    {
      Couple_for_ElemValue.theActualCoeffProc = Couple_for_ElemValue.theCoeffProc[i];
      return (0);
    }

  return (1);
}

static DOUBLE CoeffProcElementValueEvalProc (const ELEMENT *theElement,const DOUBLE **CornersCoord, DOUBLE *LocalCoord)
{
  INT i,n;
  DOUBLE phi;
  DOUBLE_VECTOR EvalPoint;

  n = CORNERS_OF_ELEM(theElement);
  V_DIM_CLEAR(EvalPoint)
  for (i=0; i<n; i++)
  {
    phi = GN(n,i,LocalCoord);
    V_DIM_LINCOMB(1.0,EvalPoint,phi,CornersCoord[i],EvalPoint)
  }
  (*Couple_for_ElemValue.theActualCoeffProc)(EvalPoint,&phi);

  return(phi);
}

EVALUES * NS_DIM_PREFIX CreateElementValueEvalProcFromCoeffProc (const char *name, CoeffProcPtr CoeffProc)
{
  EVALUES *newElementValues;

  /* check if there is space */
  if (Couple_for_ElemValue.nUsed >= MAX_COEFFPROC_ELEMEVAL)
    return(NULL);

  /* change to directory */
  if (ChangeEnvDir("/ElementEvalProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newElementValues = (EVALUES*) MakeEnvItem (name,theElemValVarID,sizeof(EVALUES));
  if (newElementValues==NULL) return(NULL);

  /* fill in data */
  newElementValues->PreprocessProc = CoeffProcElementValuePreProc;
  newElementValues->EvalProc = CoeffProcElementValueEvalProc;

  /* store coupling */
  strcpy(Couple_for_ElemValue.ElemEvalName[Couple_for_ElemValue.nUsed],name);
  Couple_for_ElemValue.theCoeffProc[Couple_for_ElemValue.nUsed] = CoeffProc;
  Couple_for_ElemValue.nUsed++;

  UserWrite("ElementValueEvalProc "); UserWrite(name); UserWrite(" installed\n");

  return(newElementValues);
}

/****************************************************************************/
/*
   CreateElementVectorEvalProcFromCoeffProc - Create element vector plot proceedure

   SYNOPSIS:
   EVECTOR *CreateElementVectorEvalProcFromCoeffProc (const char *name,
   CoeffProcPtr CoeffProc, INT d);

   PARAMETERS:
   .  name -
   .  CoeffProc -
   .  d - dimension of vector

   DESCRIPTION:
   This function creates element vector plot proceedure in evironement.

   RETURN VALUE:
   EVECTOR *
   .n      pointer to EVECTOR *
   .n      Null if not found.
 */
/****************************************************************************/

static INT CoeffProcElementVectorPreProc (const char *name, MULTIGRID *theMG)
{
  INT i;

  for (i=0; i<Couple_for_ElemVector.nUsed; i++)
    if (strcmp(Couple_for_ElemVector.ElemEvalName[i],name)==0)
    {
      Couple_for_ElemVector.theActualCoeffProc = Couple_for_ElemVector.theCoeffProc[i];
      return (0);
    }

  return (1);
}

static void CoeffProcElementVectorEvalProc (const ELEMENT *theElement,const DOUBLE **CornersCoord, DOUBLE *LocalCoord, DOUBLE *values)
{
  INT i,n;
  DOUBLE phi;
  DOUBLE_VECTOR EvalPoint;

  n = CORNERS_OF_ELEM(theElement);
  V_DIM_CLEAR(EvalPoint)
  for (i=0; i<n; i++)
  {
    phi = GN(n,i,LocalCoord);
    V_DIM_LINCOMB(1.0,EvalPoint,phi,CornersCoord[i],EvalPoint)
  }
  (*Couple_for_ElemVector.theActualCoeffProc)(EvalPoint,values);

  return;
}

EVECTOR * NS_DIM_PREFIX CreateElementVectorEvalProcFromCoeffProc (const char *name, CoeffProcPtr CoeffProc, INT d)
{
  EVECTOR *newElemVector;

  /* check if there is space */
  if (Couple_for_ElemVector.nUsed >= MAX_COEFFPROC_ELEMEVAL)
    return(NULL);

  /* change directory */
  if (ChangeEnvDir("/ElementVectorEvalProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newElemVector = (EVECTOR*) MakeEnvItem (name,theElemVectorVarID,sizeof(EVECTOR));
  if (newElemVector==NULL) return(NULL);

  /* fill in data */
  newElemVector->PreprocessProc = CoeffProcElementVectorPreProc;
  newElemVector->EvalProc = CoeffProcElementVectorEvalProc;
  newElemVector->dimension = d;

  /* store coupling */
  strcpy(Couple_for_ElemVector.ElemEvalName[Couple_for_ElemVector.nUsed],name);
  Couple_for_ElemVector.theCoeffProc[Couple_for_ElemVector.nUsed] = CoeffProc;
  Couple_for_ElemVector.nUsed++;

  UserWrite("ElementVectorEvalProc "); UserWrite(name); UserWrite(" installed\n");
  return(newElemVector);
}

/****************************************************************************/
/*																			*/
/*				  ug plotfunctions											*/
/*																			*/
/****************************************************************************/

/**************************************************************************/
/*
   NodeIndex - plots vector index of the node vectors

   SYNOPSIS:
   static DOUBLE NodeIndex (const ELEMENT *theElement,
   const DOUBLE **CornersCoord, DOUBLE *LocalCoord);

   PARAMETERS:
   .  theElement - pointer to an element
   .  CornersCoord - corner coordinates
   .  LocalCoord - local coordinate

   DESCRIPTION:
   This function plots the vector index of the node vectors.

   RETURN VALUE:
   DOUBLE  value
 */
/*************************************************************************/

static INT PreprocessNodeIndex (const char *name, MULTIGRID *theMG)
{
  INT i, index;
  VECTOR *theVector;

  if (!VEC_DEF_IN_OBJ_OF_MG(theMG,NODEVEC))
    return (1);
  for (i=0; i<=CURRENTLEVEL(theMG); i++)
  {
    index = 0;
    for (theVector=FIRSTVECTOR(GRID_ON_LEVEL(theMG,i)); theVector!=NULL; theVector=SUCCVC(theVector))
      if (VOTYPE(theVector)==NODEVEC)
        VINDEX(theVector) = index++;
  }

  return (0);
}

static DOUBLE NodeIndex (const ELEMENT *theElement,const DOUBLE **CornersCoord, DOUBLE *LocalCoord)
{
  int i;
  DOUBLE phi;

  phi = 0.0;
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    phi += GN(CORNERS_OF_ELEM(theElement),i,LocalCoord)
           * VINDEX(NVECTOR(CORNER(theElement,i)));
  return(phi);
}

/**************************************************************************/
/*
   GradNodeIndex - plots the gradient multiplied with the vector index

   SYNOPSIS:
   static void GradNodeIndex (const ELEMENT *theElement,
   const DOUBLE **theCorners, DOUBLE *LocalCoord, DOUBLE *values);

   PARAMETERS:
   .  theElement - pointer to an element
   .  CornersCoord - corner coordinates
   .  LocalCoord - local coordinate
   .  values - plot vector

   DESCRIPTION:
   This function plots the gradient multiplied with the vector index.

   RETURN VALUE:
   void
 */
/*************************************************************************/

#ifdef __TWODIM__

static void GradNodeIndex (const ELEMENT *theElement, const DOUBLE **theCorners,  DOUBLE *LocalCoord, DOUBLE *values)
{
  INT i;
  DOUBLE_VECTOR theGradient[MAX_CORNERS_OF_ELEM];
  DOUBLE v,det;

  Gradients(CORNERS_OF_ELEM(theElement),theCorners,LocalCoord[0],LocalCoord[1],theGradient,&det);
  V2_CLEAR(values)
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    v = VINDEX(NVECTOR(CORNER(theElement,i)));
    V2_LINCOMB(v,theGradient[i],1.0,values,values)
  }
  return;
}

#else

static void GradNodeIndex (const ELEMENT *theElement, const DOUBLE **theCorners, DOUBLE *LocalCoord, DOUBLE *values)
{
  INT i;
  DOUBLE_VECTOR theGradient[MAX_CORNERS_OF_ELEM];
  DOUBLE v;

  TetraDerivative((ELEMENT *)theElement,theCorners,theGradient);
  V3_CLEAR(values)
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    v = VINDEX(NVECTOR(CORNER(theElement,i)));
    V3_LINCOMB(v,theGradient[i],1.0,values,values)
  }
  return;
}

#endif

/****************************************************************************/
/*
   InitEvalProc	- Init this file

   SYNOPSIS:
   INT InitEvalProc ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitEvalProc ()
{
  /* install the /ElementEvalProcs directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitEvalProc","could not changedir to root");
    return(__LINE__);
  }
  theEEvalProcDirID = GetNewEnvDirID();
  if (MakeEnvItem("ElementEvalProcs",theEEvalProcDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitEvalProc","could not install '/ElementEvalProcs' dir");
    return(__LINE__);
  }
  theElemValVarID = GetNewEnvVarID();

  /* install the /MatrixEvalProcs directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitEvalProc","could not changedir to root");
    return(__LINE__);
  }
  theMEvalProcDirID = GetNewEnvDirID();
  if (MakeEnvItem("MatrixEvalProcs",theMEvalProcDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitEvalProc","could not install '/MatrixEvalProcs' dir");
    return(__LINE__);
  }
  theMatrixValVarID = GetNewEnvVarID();

  /* install the /ElementVectorEvalProcs directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitEvalProc","could not changedir to root");
    return(__LINE__);
  }
  theVEvalProcDirID = GetNewEnvDirID();
  if (MakeEnvItem("ElementVectorEvalProcs",theVEvalProcDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitEvalProc","could not install '/ElementVectorEvalProcs' dir");
    return(__LINE__);
  }
  theElemVectorVarID = GetNewEnvVarID();

  /* install general plot procs */
  if (CreateElementValueEvalProc("nindex",PreprocessNodeIndex,NodeIndex)==NULL) return(1);
  if (CreateElementVectorEvalProc("gradnindex",PreprocessNodeIndex,GradNodeIndex,DIM)==NULL) return(1);

  /* init variables used for CoeffProcElemEval */
  Couple_for_ElemValue.nUsed  = 0;
  Couple_for_ElemVector.nUsed = 0;

  return (0);
}
