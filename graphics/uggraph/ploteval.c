// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  plotproc.c													*/
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "gm.h"
#include "plotproc.h"
#include "ugenv.h"
#include "evm.h"
#include "devices.h"

#ifdef __TWODIM__
#include "shapes2d.h"
#endif

#ifdef __THREEDIM__
#include "shapes3d.h"
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

/* environment directory and item IDs used in this source file */
static INT theEPlotProcDirID;
static INT theElemValVarID;

static INT theMPlotProcDirID;
static INT theMatrixValVarID;

static INT theVPlotProcDirID;
static INT theElemVectorVarID;

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*
   CreateElementValuePlotProc - Create element value plot proceedure in evironement

   SYNOPSIS:
   EVALUES *CreateElementValuePlotProc (const char *name,
   PreprocessingProcPtr PreProc, ElementPlotProcPtr PlotProc);

   PARAMETERS:
   .  name -
   .  PreProc -
   .  PlotProc -

   DESCRIPTION:
   This function creates element value plot proceedure in evironement.

   RETURN VALUE:
   EVALUES *
   .n     pointer to
   .n     Null if error occured.
 */
/****************************************************************************/

EVALUES *CreateElementValuePlotProc (const char *name, PreprocessingProcPtr PreProc, ElementPlotProcPtr PlotProc)
{
  EVALUES *newElementValues;

  /* change to directory */
  if (ChangeEnvDir("/ElementPlotProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newElementValues = (EVALUES*) MakeEnvItem (name,theElemValVarID,sizeof(EVALUES));
  if (newElementValues==NULL) return(NULL);

  /* fill in data */
  newElementValues->PreprocessProc = PreProc;
  newElementValues->PlotProc = PlotProc;

  UserWrite("ElementValuePlotProc "); UserWrite(name); UserWrite(" installed\n");

  return(newElementValues);
}

/****************************************************************************/
/*
   GetElementValuePlotProc - Get element value plot proceedure in evironement from name

   SYNOPSIS:
   EVALUES *GetElementValuePlotProc (const char *name);

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

EVALUES *GetElementValuePlotProc (const char *name)
{
  if (ChangeEnvDir("/ElementPlotProcs")==NULL) return(NULL);
  return((EVALUES*) SearchEnv(name,".",theElemValVarID,SEARCHALL));
}

/****************************************************************************/
/*
   CreateElementVectorPlotProc - Create element vector plot proceedure in evironement

   SYNOPSIS:
   EVECTOR *CreateElementVectorPlotProc (const char *name,
   PreprocessingProcPtr PreProc, ElementVectorProcPtr PlotProc, INT d);

   PARAMETERS:
   .  name -
   .  PreProc -
   .  PlotProc -
   .  d - dimension of vector

   DESCRIPTION:
   This function creates element vector plot proceedure in evironement.

   RETURN VALUE:
   EVECTOR *
   .n      pointer to EVECTOR *
   .n      Null if not found.
 */
/****************************************************************************/

EVECTOR *CreateElementVectorPlotProc (const char *name, PreprocessingProcPtr PreProc, ElementVectorProcPtr PlotProc, INT d)
{
  EVECTOR *newElemVector;

  /* change directory */
  if (ChangeEnvDir("/ElementVectorPlotProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newElemVector = (EVECTOR*) MakeEnvItem (name,theElemVectorVarID,sizeof(EVECTOR));
  if (newElemVector==NULL) return(NULL);

  /* fill in data */
  newElemVector->PreprocessProc = PreProc;
  newElemVector->PlotProc = PlotProc;
  newElemVector->dimension = d;

  UserWrite("ElementVectorPlotProc "); UserWrite(name); UserWrite(" installed\n");
  return(newElemVector);
}

/****************************************************************************/
/*
   GetElementVectorPlotProc -

   SYNOPSIS:
   EVECTOR *GetElementVectorPlotProc (const char *name);

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

EVECTOR *GetElementVectorPlotProc (const char *name)
{
  if (ChangeEnvDir("/ElementVectorPlotProcs")==NULL) return(NULL);
  return((EVECTOR*) SearchEnv(name,".",theElemVectorVarID,SEARCHALL));
}

/****************************************************************************/
/*
   CreateMatrixValuePlotProc - Create matrix value plot proceedure in evironement

   SYNOPSIS:
   MVALUES *CreateMatrixValuePlotProc (const char *name, PreprocessingProcPtr PreProc,
   MatrixPlotProcPtr PlotProc);

   PARAMETERS:
   .  name -
   .  PreProc -
   .  PlotProc -

   DESCRIPTION:
   This function create matrix value plot proceedure in evironement.

   RETURN VALUE:
   MVALUES *
   .n      pointer to
   .n      NULL, if error occured.
 */
/****************************************************************************/

MVALUES *CreateMatrixValuePlotProc (const char *name, PreprocessingProcPtr PreProc, MatrixPlotProcPtr PlotProc)
{
  MVALUES *newMatrixValues;

  /* change to directory */
  if (ChangeEnvDir("/MatrixPlotProcs")==NULL)
    return(NULL);

  /* allocate structure */
  newMatrixValues = (MVALUES*) MakeEnvItem (name,theMatrixValVarID,sizeof(MVALUES));
  if (newMatrixValues==NULL) return(NULL);

  /* fill in data */
  newMatrixValues->PreprocessProc = PreProc;
  newMatrixValues->PlotProc = PlotProc;

  UserWrite("MatrixValuePlotProc "); UserWrite(name); UserWrite(" installed\n");
  return(newMatrixValues);
}

/****************************************************************************/
/*
   GetMatrixValuePlotProc - Get matrix value plot proceedure in evironement from name

   SYNOPSIS:
   MVALUES *GetMatrixValuePlotProc (const char *name);

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

MVALUES *GetMatrixValuePlotProc (const char *name)
{
  if (ChangeEnvDir("/MatrixPlotProcs")==NULL) return(NULL);
  return((MVALUES*) SearchEnv(name,".",theMatrixValVarID,SEARCHALL));
}


/****************************************************************************/
/*																			*/
/*				  ug plotfunctions											*/
/*																			*/
/****************************************************************************/

#ifdef __NODEDATA__

/****************************************************************************/
/*D
   NodeVectorOrder - Plotfunction for NodeVector Order

   SYNOPSIS:
   int PreprocessNodeIndex (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - pointer to multigrid

   DESCRIPTION:
   This function is the plotfunction for NodeVector Order.

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT PreprocessNodeIndex (const char *name, MULTIGRID *theMG)
{
  INT i, index;
  VECTOR *theVector;

  for (i=0; i<=CURRENTLEVEL(theMG); i++)
  {
    index = 0;
    for (theVector=FIRSTVECTOR(GRID_ON_LEVEL(theMG,i)); theVector!=NULL; theVector=SUCCVC(theVector))
      if (VTYPE(theVector)==NODEVECTOR)
        VINDEX(theVector) = index++;
  }

  return (0);
}
/**************************************************************************/
/*
   NodeIndex -

   SYNOPSIS:
   static DOUBLE NodeIndex (const ELEMENT *theElement,const COORD **CornersCoord,
   COORD *LocalCoord);

   PARAMETERS:
   .  theElement -
   .  CornersCoord -
   .  LocalCoord -

   DESCRIPTION:
   This function

   RETURN VALUE:
   DOUBLE
 */
/*************************************************************************/
static DOUBLE NodeIndex (const ELEMENT *theElement,const COORD **CornersCoord, COORD *LocalCoord)
{
  int i;
  DOUBLE phi;

  phi = 0.0;
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
        #ifdef __TWODIM__
    phi += N(CORNERS_OF_ELEM(theElement),i,LocalCoord[0],LocalCoord[1])*VINDEX(NVECTOR(CORNER(theElement,i)));
        #else
    phi += N(i,LocalCoord)*VINDEX(NVECTOR(CORNER(theElement,i)));
        #endif
  return(phi);
}

/**************************************************************************/
/*
   GradNodeIndex -

   SYNOPSIS:
   static void GradNodeIndex (const ELEMENT *theElement,
   const COORD **theCorners, COORD *LocalCoord, DOUBLE *values);

   PARAMETERS:
   .  theElement -
   .  theCorners -
   .  LocalCoord -
   .  values -

   DESCRIPTION:
   This function

   RETURN VALUE:
   void
 */
/*************************************************************************/

#ifdef __TWODIM__

static void GradNodeIndex (const ELEMENT *theElement, const COORD **theCorners,  COORD *LocalCoord, DOUBLE *values)
{
  int i;
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

static void GradNodeIndex (const ELEMENT *theElement, const COORD **theCorners, COORD *LocalCoord, DOUBLE *values)
{
  int i;
  COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM];
  DOUBLE v;

  TetraDerivative(theCorners,theGradient);
  V3_CLEAR(values)
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    v = VINDEX(NVECTOR(CORNER(theElement,i)));
    V3_LINCOMB(v,theGradient[i],1.0,values,values)
  }
  return;
}

#endif

#endif

/****************************************************************************/
/*
   InitPlotproc	- Init this file

   SYNOPSIS:
   INT InitPlotproc ();

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

INT InitPlotproc ()
{
  /* install the /ElementPlotProcs directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitPlotproc","could not changedir to root");
    return(__LINE__);
  }
  theEPlotProcDirID = GetNewEnvDirID();
  if (MakeEnvItem("ElementPlotProcs",theEPlotProcDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitPlotproc","could not install '/ElementPlotProcs' dir");
    return(__LINE__);
  }
  theElemValVarID = GetNewEnvVarID();

  /* install the /MatrixPlotProcs directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitPlotproc","could not changedir to root");
    return(__LINE__);
  }
  theMPlotProcDirID = GetNewEnvDirID();
  if (MakeEnvItem("MatrixPlotProcs",theMPlotProcDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitPlotproc","could not install '/MatrixPlotProcs' dir");
    return(__LINE__);
  }
  theMatrixValVarID = GetNewEnvVarID();

  /* install the /ElementVectorPlotProcs directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitPlotproc","could not changedir to root");
    return(__LINE__);
  }
  theVPlotProcDirID = GetNewEnvDirID();
  if (MakeEnvItem("ElementVectorPlotProcs",theVPlotProcDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitPlotproc","could not install '/ElementVectorPlotProcs' dir");
    return(__LINE__);
  }
  theElemVectorVarID = GetNewEnvVarID();

  /* install general plot procs */
        #ifdef __NODEDATA__
  if (CreateElementValuePlotProc("nindex",PreprocessNodeIndex,NodeIndex)==NULL) return(1);
  if (CreateElementVectorPlotProc("gradnindex",PreprocessNodeIndex,GradNodeIndex,DIM)==NULL) return(1);
        #endif


  return (0);
}
