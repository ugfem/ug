// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  transfer.c													*/
/*																			*/
/* Purpose:   grid transfer num procs                                           */
/*																			*/
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   November 29, 1996                                                                         */
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

#include "general.h"
#include "gm.h"
#include "np.h"

#include "transfer.h"

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

typedef struct
{
  NP_TRANSFER transfer;

  TransGridProcPtr res;
  TransGridProcPtr intcor;
  InterpolateSolutionProcPtr intnew;

} NP_STANDARD_TRANSFER;

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

INT NPTransferInit (NP_TRANSFER *np, INT argc , char **argv)
{
  INT i;

  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->c = ReadArgvVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);

  if (sc_read(np->damp,np->x,"damp",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->damp[i] = 1.0;

  if ((np->A == NULL) && (np->b == NULL) && (np->x == NULL) && (np->c == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPTransferDisplay (NP_TRANSFER *np)
{
  if ((np->A == NULL) && (np->x == NULL) && (np->b == NULL) && (np->c == NULL))
    return(0);

  UserWrite("symbolic user data:\n");
  if (np->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  if (np->b != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->b));
  if (np->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
  if (np->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->x));
  UserWrite("\n");

  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->b,"damp"))
    return (1);

  return(0);
}

INT NPTransferExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_TRANSFER *np;
  INT result,level,bl;

  np = (NP_TRANSFER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (ReadOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no PreProcess");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if (np->b == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector b");
      return (1);
    }
    if (np->A == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no matrix A");
      return (1);
    }
    if ((*np->PreProcess)(np,level,np->x,np->b,np->A,&bl,&result)) {
      UserWriteF("NPTransferExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("R",argc,argv)) {
    if (np->RestrictDefect == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no RestrictDefect");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if (np->A == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no matrix A");
      return (1);
    }
    if ((*np->RestrictDefect)(np,level,np->b,np->b,np->A,np->damp,&result)) {
      UserWriteF("NPTransferExecute: RestrictDefect failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("I",argc,argv)) {
    if (np->InterpolateCorrection == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no InterpolateCorrection");
      return (1);
    }
    if (np->c == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector c");
      return (1);
    }
    if (np->A == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no matrix A");
      return (1);
    }
    if ((*np->InterpolateCorrection)(np,level,np->c,np->c,np->A,np->damp,&result)) {
      UserWriteF("NPTransferExecute: InterpolateCorrection failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("N",argc,argv)) {
    if (np->InterpolateNewVectors == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no InterpolateNewVectors");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if ((*np->InterpolateNewVectors)(np,level,np->x,&result)) {
      UserWriteF("NPTransferExecute: InterpolateNewVectors failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("P",argc,argv)) {
    if (np->ProjectSolution == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no ProjectSolution");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if ((*np->ProjectSolution)(np,level,np->x,&result)) {
      UserWriteF("NPTransferExecute: ProjectSolution failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no PostProcess");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if (np->b == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector b");
      return (1);
    }
    if (np->A == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no matrix A");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPTransferExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  return(0);
}

static INT TransferInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;

  np->res = StandardRestrict;
  np->intcor = StandardInterpolateCorrection;
  np->intnew = StandardInterpolateNewVectors;

  if (ReadArgvOption("M",argc,argv)) {
    np->res = RestrictByMatrix;
    np->intcor = InterpolateCorrectionByMatrix;
    np->intnew = InterpolateNewVectorsByMatrix;
  }

  return (TransferInitSymbols(&np->transfer,argc,argv));
}

static INT TransferDisplay (NP_BASE *theNP)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;

  TransferDisplaySymbols(&np->transfer);

  if (np->res == StandardRestrict)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Restrict","StandardRestrict");
  if (np->res == RestrictByMatrix)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Restrict","RestrictByMatrix");
  if (np->intcor = StandardInterpolateCorrection)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateCor","StandardInterpolateCorrection");
  if (np->intcor == InterpolateCorrectionByMatrix)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateCor","InterpolateCorrectionByMatrix");
  if (np->intnew == StandardInterpolateNewVectors)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateNew","StandardInterpolateNewVectors");
  if (np->intnew == InterpolateNewVectorsByMatrix)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateNew","StandardRestrict");

  return (0);
}

static INT TransferPreProcess (NP_TRANSFER *theNP, INT level,
                               VECDATA_DESC *from, VECDATA_DESC *to,
                               MATDATA_DESC *A, INT *baselevel, INT *result)
{
  return(0);
}

static INT RestrictDefect (NP_TRANSFER *theNP, INT level,
                           VECDATA_DESC *from, VECDATA_DESC *to,
                           MATDATA_DESC *A, VEC_SCALAR damp,
                           INT *result)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;
  result[0] = (*np->res)(GRID_ON_LEVEL(theNP->base.mg,level),from,to,damp);

  return(result[0]);
}

static INT InterpolateCorrection (NP_TRANSFER *theNP, INT level,
                                  VECDATA_DESC *from, VECDATA_DESC *to,
                                  MATDATA_DESC *A, VEC_SCALAR damp,
                                  INT *result)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;
  result[0] = (*np->intcor)(GRID_ON_LEVEL(theNP->base.mg,level),from,to,damp);

  return(result[0]);
}

static INT InterpolateNewVectors (NP_TRANSFER *theNP, INT level,
                                  VECDATA_DESC *x, INT *result)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;
  result[0] = (*np->intnew)(GRID_ON_LEVEL(theNP->base.mg,level),x);

  return(result[0]);
}

static INT TransferPostProcess (NP_TRANSFER *theNP, INT level,
                                VECDATA_DESC *from, VECDATA_DESC *to,
                                MATDATA_DESC *A, INT *result)
{
  return(0);
}

static INT TransferConstruct (NP_BASE *theNP)
{
  NP_TRANSFER *np;

  theNP->Init = TransferInit;
  theNP->Display = TransferDisplay;
  theNP->Execute = NPTransferExecute;

  np = (NP_TRANSFER *) theNP;
  np->PreProcess = TransferPreProcess;
  np->RestrictDefect = RestrictDefect;
  np->InterpolateCorrection = InterpolateCorrection;
  np->InterpolateNewVectors = InterpolateNewVectors;
  np->ProjectSolution = NULL;
  np->PostProcess = TransferPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitTransfer	- Init this file

   SYNOPSIS:
   INT InitPlotProc ();

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

INT InitTransfer ()
{
  if (CreateClass("transfer",sizeof(NP_STANDARD_TRANSFER),TransferConstruct))
    return (__LINE__);

  return (0);
}
