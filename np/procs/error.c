// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      error.c                                                       */
/*                                                                          */
/* Purpose:   most simple error indicators                                                      */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de				                        */
/*																			*/
/* History:   Sep 4, 1996, ug version 3.4                                                               */
/*            December 8, 1996, new np subsystem                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "gm.h"       /* for data structure               */
#include "evm.h"      /* for data structure               */
#include "shapes.h"   /* for data structure               */
#include "np.h"
#include "ugm.h"
#include "quadrature.h"
#include "disctools.h"
#include "general.h"
#include "devices.h"
#include "ugstruct.h"

#include "error.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

typedef struct {

  NP_ERROR error;

  INT from;
  INT to;
  DOUBLE refine;
  DOUBLE coarse;

  INT clear;
  INT update;
  INT interpolate;
  INT project;

} NP_INDICATOR;

/****************************************************************************/
/*D
   NP_ERROR - type definition for error indicators

   DESCRIPTION:
   This numproc type is used for the configuration of error indicators.
   PreProcessError has to be called once befor the first call.
   Error can be applied several times.
   Finally, PostProcessError has to be called once at the very end.

   Initializing the solution Symbol is optional; it can be done with

   'INT NPErrorInit (NP_ERROR *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPErrorDisplay (NP_ERROR *theNP);'
   'INT NPErrorExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb


   ..... fill in data structure here when the realizition is finished


   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPErrorInit (NP_ERROR *np, INT argc , char **argv)
{
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->o = ReadArgvVecDesc(np->base.mg,"o",argc,argv);

  if (np->x == NULL)
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPErrorDisplay (NP_ERROR *np)
{
  if ((np->x == NULL) && (np->o == NULL))
    return(0);
  UserWrite("symbolic user data:\n");
  if (np->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
  if (np->o != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"o",ENVITEM_NAME(np->o));
  UserWrite("\n");

  return(0);
}

INT NPErrorExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ERROR *np;
  DOUBLE Time,step;
  INT result,level;
  ERESULT eresult;

  np = (NP_ERROR *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"NPErrorExecute","no vector x");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,&result)) {
      UserWriteF("NPErrorExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("e",argc,argv)) {
    if (np->Error == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no PreProcess");
      return (1);
    }
    if ((*np->Error)(np,level,np->x,&eresult)) {
      UserWriteF("NPErrorExecute: Error failed, error code %d\n",
                 eresult.error_code);
      return (1);
    }
  }

  if (ReadArgvOption("t",argc,argv)) {
    if (np->TimeError == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no PreProcess");
      return (1);
    }
    if (np->o == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no vector o");
      return (1);
    }
    if (ReadArgvDOUBLE("t",&Time,argc,argv)) {
      PrintErrorMessage('E',"NPErrorExecute","no time");
      return (1);
    }
    if (ReadArgvDOUBLE("s",&step,argc,argv)) {
      PrintErrorMessage('E',"NPErrorExecute","no time step");
      return (1);
    }
    if ((*np->TimeError)(np,level,Time,&step,np->x,np->o,&eresult)) {
      UserWriteF("NPErrorExecute: PreProcess failed, error code %d\n",
                 eresult.error_code);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,&result)) {
      UserWriteF("NPErrorExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  return(0);
}

/****************************************************************************/
/*D
   indicator - numproc for the compuation of coarsening marks

   DESCRIPTION:
   This numproc computes a hierachical basis representation for the
   surface elements and marks them for coarsening or refinement.
   It requieres node values and uses the shapefunctions defined in shapes.c.

   'npinit $x <vec sym> [$coarse <val>] [$refine <val>]'

   .  $x~<vec~sym> - symbol for the (solution) vector
   .  $coarse~<val> - bound for coarsening
   .  $refine~<val> - bound for refinement

   'npexecute <name> [$p] [$r] [$i]'

   .  $p - calls StandardProject
   .  $r - calls RefineMultigrid
   .  $i - calls StandardInterpolateNewVectors

   The estimator runs the program 'SurfaceIndicator'.

   EXAMPLE:
   .vb
   npcreate indicator $t error;
   npinit $x sol $coarse 0.5;

   npexecute coarse $p;
   .ve
   D*/
/****************************************************************************/

/* gradient computation using the general shape functions of shape.c */

static INT GradientAtLMP (ELEMENT *theElement, INT ncomp, VECDATA_DESC *theVD,
                          DOUBLE_VECTOR *grad, DOUBLE *area)
{
  DOUBLE_VECTOR gradient[MAX_SINGLE_VEC_COMP];
  DOUBLE_VECTOR Jinv[DIM];
  DOUBLE Jdet;
  DOUBLE *vptr[MAX_NODAL_VALUES];
  DOUBLE_VECTOR derivative;
  DOUBLE *x[MAX_CORNERS_OF_ELEM],*lmp;
  INT i,j,m,n;

  CORNER_COORDINATES(theElement,n,x);
  GetElementVPtrs (theElement,theVD,vptr);

  for (j=0; j<ncomp; j++)
    V_DIM_CLEAR (gradient[j]);

  m = 0;
  lmp = LMP(n);
  for (i=0; i<n; i++)
  {
    D_GN (n,i,lmp,derivative);
    for (j=0; j<ncomp; j++)
    {
      V_DIM_LINCOMB (1.0,gradient[j],*vptr[m],derivative,gradient[j]);
      m++;
    }
  }
  INVERSE_TRANSFORMATION(n,x,lmp,Jinv,Jdet);
  for (j=0; j<ncomp; j++)
    MM_TIMES_V_DIM (Jinv,gradient[j],grad[j]);

  AREA_OF_REF(n,*area);
  *area *= ABS(Jdet);

  return(0);
}

static DOUBLE ElementIndicator (ELEMENT *t, INT ncomp, VECDATA_DESC *theVD)
{
  ELEMENT *f;
  DOUBLE_VECTOR egrad[MAX_SINGLE_VEC_COMP];
  DOUBLE_VECTOR fgrad[MAX_SINGLE_VEC_COMP];
  INT i;
  DOUBLE est,diff,area;

  f = EFATHER(t);
  if (f == NULL)
    return(0.0);

  GradientAtLMP (f,ncomp,theVD,fgrad,&area);
  GradientAtLMP (t,ncomp,theVD,egrad,&area);

  est = 0.0;
  for (i=0; i<ncomp; i++)
  {
    V_DIM_EUKLIDNORM_OF_DIFF(egrad[i],fgrad[i],diff);
    est += diff;
  }

  return(est * area);
}

INT SurfaceIndicator (MULTIGRID *theMG, VECDATA_DESC *theVD,
                      DOUBLE refine, DOUBLE coarse, INT project,
                      INT from, INT to, INT clear, ERESULT *eresult)
{
  ELEMENT *t;
  DOUBLE *List,min,max,est,rf,cr;
  INT k,toplevel,nel,mfr,mfc,ncomp;

  ncomp = VD_NCMPS_IN_TYPE(theVD,NODEVECTOR);
  if (ncomp == 0)
    return(1);

  /* toplevel */
  toplevel = TOPLEVEL(theMG);
  if (toplevel<0)
    return(1);

  if (project)
    for (k=toplevel-1; k>=0; k--)
      if (StandardProject (GRID_ON_LEVEL(theMG,k),
                           (const VECDATA_DESC *)theVD,
                           (const VECDATA_DESC *)theVD))
        return(1);

  /* count the number of red surface elements */
  nel = 0;
  for (k=0; k<=toplevel; k++)
    for (t=FIRSTELEMENT(GRID_ON_LEVEL(theMG,k)); t!=NULL; t=SUCCE(t))
      if (EstimateHere(t) && ECLASS(t) != YELLOW_CLASS)
      {
        nel++;
        if (clear)
          MarkForRefinement(t,NO_REFINEMENT,0);
      }

  Mark(MGHEAP(theMG),FROM_TOP);
  List = (DOUBLE*) GetMem(MGHEAP(theMG),nel*sizeof(DOUBLE),FROM_TOP);
  if (List == NULL)
    return(-1);

  max = 0.0;
  min = 1e20;
  nel = 0;
  mfr = 0;
  mfc = 0;
  for (k=0; k<=toplevel; k++)
    for (t=FIRSTELEMENT(GRID_ON_LEVEL(theMG,k)); t!=NULL; t=SUCCE(t))
      if (EstimateHere(t) && ECLASS(t) != YELLOW_CLASS)
      {
        est = ElementIndicator(t,ncomp,theVD);
        min = MIN(min,est);
        max = MAX(max,est);
        List[nel++] = est;
      }

  rf = max * refine;
  cr = max * coarse;
  nel = 0;
  for (k=0; k<=toplevel; k++)
    for (t=FIRSTELEMENT(GRID_ON_LEVEL(theMG,k)); t!=NULL; t=SUCCE(t))
      if (EstimateHere(t) && ECLASS(t) != YELLOW_CLASS)
      {
        est = List[nel++];
        if ((est > rf) && (k < to))
        {
          MarkForRefinement(t,RED,0);
          mfr++;
        }
        else if ((est < cr) && (k > from))
        {
          MarkForRefinement(t,COARSE,0);
          mfc++;
        }
      }
  Release(MGHEAP(theMG),FROM_TOP);

  if (SetStringValue("indicator:mfr",(DOUBLE)mfr))
    return (-1);
  if (SetStringValue("indicator:mfc",(DOUBLE)mfc))
    return (-1);

  /* print result */
  UserWrite("Indicator:");
  if (mfr > 0)
    UserWriteF(" %d elements marked for refinement",mfr);
  if (mfc > 0)
    UserWriteF("    %d elements marked for coarsening",mfc);
  UserWrite("\n");

  /* return number of flagged elements */
  eresult->nel = nel;
  eresult->refine = mfr;
  eresult->coarse = mfc;

  return(0);
}

static INT IndicatorInit (NP_BASE *theNumProc, INT argc, char **argv)
{
  NP_INDICATOR *theNP;

  theNP = (NP_INDICATOR*)theNumProc;

  if (ReadArgvINT("from",&(theNP->from),argc,argv))
    theNP->from = 0;
  if (ReadArgvINT("to",&(theNP->to),argc,argv))
    theNP->to = 32;

  if (ReadArgvDOUBLE("refine",&(theNP->refine),argc,argv))
    theNP->refine = 2.0;
  if (ReadArgvDOUBLE("coarse",&(theNP->coarse),argc,argv))
    theNP->coarse = 0.0;

  return (NPErrorInit(&theNP->error,argc,argv));
}

static INT IndicatorDisplay (NP_BASE *theNumProc)
{
  NP_INDICATOR *theNP;

  theNP = (NP_INDICATOR*)theNumProc;

  NPErrorDisplay(&theNP->error);

  UserWriteF(DISPLAY_NP_FORMAT_SI,"from level",(int)theNP->from);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"to level",(int)theNP->to);
  if (theNP->refine < 1.0)
    UserWriteF("%-16.13s = %-12.9lf\n","refine",theNP->refine);
  if (theNP->coarse > 0.0)
    UserWriteF("%-16.13s = %-12.9lf\n","refine",theNP->coarse);

  return(0);
}

static INT Indicator (NP_ERROR *theNP, INT level, VECDATA_DESC *x,
                      ERESULT *eresult)
{
  NP_INDICATOR *np;
  MULTIGRID *theMG;
  GRID *theGrid;
  INT i;

  np = (NP_INDICATOR*) theNP;
  theMG = theNP->base.mg;

  if (SurfaceIndicator(theMG,x,np->refine,np->coarse,
                       np->project,np->from,np->to,np->clear,eresult) == -1) {
    eresult->error_code = __LINE__;
    return(1);
  }
  i = 0;
  if (np->update) {
    i = 1;
    if (RefineMultiGrid(theMG,GM_REFINE_TRULY_LOCAL) != GM_OK) {
      eresult->error_code = __LINE__;
      return(1);
    }
    UserWrite("[r]");
  }
  if (np->interpolate) {
    for (i=1; i<=TOPLEVEL(theMG); i++) {
      theGrid = GRID_ON_LEVEL(theMG,i);
      if (GSTATUS(theGrid)&1) {
        GSTATUS(theGrid) &= 0xFFFFFFFE;
        if (StandardInterpolateNewVectors
              (theGrid,(const VECDATA_DESC *)x) != NUM_OK) {
          eresult->error_code = __LINE__;
          return(1);
        }
        UserWriteF(" [i%d]",i);
      }
    }
  }
  if (i > 0) UserWrite("\n");

  return(0);
}

static INT IndicatorExecute (NP_BASE *theNumProc, INT argc, char **argv)
{
  NP_INDICATOR *theNP;
  ERESULT eresult;

  theNP = (NP_INDICATOR*)theNumProc;

  theNP->clear = ReadArgvOption("c",argc,argv);
  theNP->project = ReadArgvOption("p",argc,argv);
  theNP->update = ReadArgvOption("r",argc,argv);
  theNP->interpolate = ReadArgvOption("i",argc,argv);

  if (Indicator(&theNP->error,CURRENTLEVEL(theNumProc->mg),
                theNP->error.x,&eresult)) {
    UserWriteF("Indicator failed, error code %d\n",eresult.error_code);
    return (1);
  }

  return (NUM_OK);
}

static INT IndicatorConstruct (NP_BASE *theNP)
{
  NP_ERROR *np;

  theNP->Init = IndicatorInit;
  theNP->Display = IndicatorDisplay;
  theNP->Execute = IndicatorExecute;

  np = (NP_ERROR *) theNP;
  np->PreProcess = NULL;
  np->Error = Indicator;
  np->TimeError = NULL;
  np->PostProcess = NULL;

  return(0);
}

/****************************************************************************/
/*
   InitError - Init this file

   SYNOPSIS:
   INT InitError ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

INT InitError (void)
{
  if (CreateClass (ERROR_CLASS_NAME ".indicator",
                   sizeof(NP_INDICATOR), IndicatorConstruct))
    return (__LINE__);

  if (MakeStruct(":indicator"))
    return (__LINE__);

  return(0);
}
