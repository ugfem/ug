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
#include "debug.h"
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
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

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
   typedef struct {
        INT error_code;                     // error code
        INT nel;                            // number of surface elements
        INT refine;                         // nb. el. marked for refinement
        INT coarse;                         // nb. el. marked for coarsening
        DOUBLE error;                       // error estimate
   } ERESULT;

   struct np_error {

        NP_BASE base;                        // inherits base class

        // data (optinal, necessary for calling the generic execute routine)
    VECDATA_DESC *x;                     // solution
    VECDATA_DESC *o;                     // old solution
    NP_T_SOLVER *ts;                     // reference to timesolver

        // functions
        INT (*PreProcess)
             (struct np_error *,             // pointer to (derived) object
                  INT,                           // level
                  INT *);                        // result
    INT (*Error)
             (struct np_error *,             // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  ERESULT *);                    // result
    INT (*TimeError)
             (struct np_error *,             // pointer to (derived) object
                  INT,                           // level
                  DOUBLE,                        // time
                  DOUBLE *,                      // time step
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // old solution vector
          NP_T_SOLVER *,                 // reference to timesolver
                  ERESULT *);                    // result
        INT (*PostProcess)
             (struct np_error *,             // pointer to (derived) object
                  INT,                           // level
                  INT *);                        // result
   };
   typedef struct np_error NP_ERROR;
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
    REP_ERR_RETURN(1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no PreProcess");
      REP_ERR_RETURN(1);
    }
    if ((*np->PreProcess)(np,level,&result)) {
      UserWriteF("NPErrorExecute: PreProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN(1);
    }
  }

  if (ReadArgvOption("e",argc,argv)) {
    if (np->Error == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no Error");
      REP_ERR_RETURN(1);
    }
    if ((*np->Error)(np,level,np->x,&eresult)) {
      UserWriteF("NPErrorExecute: Error failed, error code %d\n",
                 eresult.error_code);
      REP_ERR_RETURN(1);
    }
  }

  if (ReadArgvOption("t",argc,argv)) {
    if (np->TimeError == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no PreProcess");
      REP_ERR_RETURN(1);
    }
    if (np->o == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no vector o");
      REP_ERR_RETURN(1);
    }
    if (ReadArgvDOUBLE("t",&Time,argc,argv)) {
      PrintErrorMessage('E',"NPErrorExecute","no time");
      REP_ERR_RETURN(1);
    }
    if (ReadArgvDOUBLE("s",&step,argc,argv)) {
      PrintErrorMessage('E',"NPErrorExecute","no time step");
      REP_ERR_RETURN(1);
    }
    if ((*np->TimeError)(np,level,Time,&step,np->x,np->o,np->ts,&eresult)) {
      UserWriteF("NPErrorExecute: PreProcess failed, error code %d\n",
                 eresult.error_code);
      REP_ERR_RETURN(1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPErrorExecute","no PostProcess");
      REP_ERR_RETURN(1);
    }
    if ((*np->PostProcess)(np,level,&result)) {
      UserWriteF("NPErrorExecute: PostProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN(1);
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

   'npinit <name> $x <vec sym> [$coarse <val>] [$refine <val>];'

   .  $x~<vec~sym> - symbol for the (solution) vector
   .  $coarse~<val> - bound for coarsening
   .  $refine~<val> - bound for refinement

   'npexecute <name> [$p] [$r] [$i];'

   .  $p - calls StandardProject
   .  $r - calls RefineMultigrid
   .  $i - calls StandardInterpolateNewVectors

   The estimator runs the program 'SurfaceIndicator'.

   EXAMPLE:
   .vb
   npcreate coarse $c indicator;
   npinit coarse $x sol $coarse 0.5;

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

static DOUBLE ElementIndicator_christian (ELEMENT *t, INT ncomp, VECDATA_DESC *theVD)
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

static DOUBLE ElementIndicator (ELEMENT *t, INT ncomp, VECDATA_DESC *theVD)
{
  DOUBLE theMin=1.0E100, theMax=-1.0E100;
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(t); i++)
  {
    theMin = MIN(theMin,VVALUE(NVECTOR(CORNER(t,i)),VD_CMP_OF_TYPE(theVD,NODEVEC,1)));
    theMax = MAX(theMax,VVALUE(NVECTOR(CORNER(t,i)),VD_CMP_OF_TYPE(theVD,NODEVEC,1)));
  }

  return(theMax-theMin);
}

INT SurfaceIndicator (MULTIGRID *theMG, VECDATA_DESC *theVD,
                      DOUBLE refine, DOUBLE coarse, INT project,
                      INT from, INT to, INT clear, ERESULT *eresult)
{
  ELEMENT *t;
  DOUBLE *List,min,max,est,rf,cr;
  INT k,toplevel,nel,mfr,mfc,ncomp;

  ncomp = VD_ncmps_in_otype(theVD,NODEVEC);
  if (ncomp <= 0)
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
        if (est > rf)
        {
          if (MarkForRefinementX(t,from,to,RED,0) == GM_OK)
            mfr++;
        }
        if (est < cr)
        {
          if (MarkForRefinementX(t,from,to,COARSE,0) == GM_OK)
            mfc++;
        }
      }
  Release(MGHEAP(theMG),FROM_TOP);

        #ifdef ModelP
  mfr = UG_GlobalSumDOUBLE((DOUBLE)mfr);
  mfc = UG_GlobalSumDOUBLE((DOUBLE)mfc);
        #endif

  if (SetStringValue("indicator:mfr",(DOUBLE)mfr))
    return (-1);
  if (SetStringValue("indicator:mfc",(DOUBLE)mfc))
    return (-1);

  /* print result */
  if (mfr + mfc > 0)
    UserWrite("Indicator:");
  if (mfr > 0)
    UserWriteF(" %d elements marked for refinement",mfr);
  if (mfc > 0)
    UserWriteF("    %d elements marked for coarsening",mfc);
  if (mfr + mfc > 0)
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

  theNP->project = ReadArgvOption("p",argc,argv);
  theNP->update = ReadArgvOption("r",argc,argv);
  theNP->interpolate = ReadArgvOption("i",argc,argv);
  theNP->clear = ReadArgvOption("c",argc,argv);

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
    UserWriteF("%-16.13s = %-12.9lf\n","coarse",theNP->coarse);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"p",(int)theNP->project);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"r",(int)theNP->update);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"i",(int)theNP->interpolate);

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
  theMG = NP_MG(theNP);

  if (SurfaceIndicator(theMG,x,np->refine,np->coarse,
                       np->project,np->from,np->to,np->clear,eresult) == -1) NP_RETURN(1,eresult->error_code);
  i = 0;
  if (np->update) {
    i = 1;
    if (RefineMultiGrid(theMG,GM_REFINE_TRULY_LOCAL,GM_REFINE_PARALLEL,GM_REFINE_NOHEAPTEST) != GM_OK) NP_RETURN(1,eresult->error_code);
    UserWrite("[r]");
  }
  if (np->interpolate) {
    for (i=1; i<=TOPLEVEL(theMG); i++) {
      theGrid = GRID_ON_LEVEL(theMG,i);
      if (GSTATUS(theGrid,GSTATUS_INTERPOLATE)) {
        RESETGSTATUS(theGrid,GSTATUS_INTERPOLATE);
        if (StandardInterpolateNewVectors
              (theGrid,(const VECDATA_DESC *)x) != NUM_OK) NP_RETURN(1,eresult->error_code);
        UserWriteF(" [i%d]",i);
      }
    }
  }
  if (i > 0) UserWrite("\n");

  return(0);
}

static INT TimeIndicator (NP_ERROR *theNP, INT level, DOUBLE t,
                          DOUBLE *dt, VECDATA_DESC *x,
                          VECDATA_DESC *o, NP_T_SOLVER *ts, ERESULT *eresult)
{
  return(Indicator(theNP,level,x,eresult));
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
    REP_ERR_RETURN(1);
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
  np->TimeError = TimeIndicator;
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
