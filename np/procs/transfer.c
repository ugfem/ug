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
#include "debug.h"
#include "devices.h"
#include "gm.h"
#include "pcr.h"
#include "np.h"
#include "disctools.h"

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

#define STANDARD_MODE           1               /* transfer via shape functions			*/
#define IMAT_MODE                       2               /* $M option							*/
#define SCALEDMG_MODE           3               /* $S option							*/

typedef struct
{
  NP_TRANSFER transfer;

  NP_TRANSFER *amg;                            /* reference to algebraic mg     */
  TransGridProcPtr res;
  TransGridProcPtr intcor;
  InterpolateNewVectorsProcPtr intnew;

  MATDATA_DESC *L;
  VECDATA_DESC *t;
  INT mode;                                                                /* mode selected in init			*/
  DOUBLE cut;                                                              /* cut value for scaled mg		*/
  INT display;                                 /* display modus                 */
  INT level;                                   /* level optimization            */
  INT meanvalue;                               /* in parallel for nonconforming */
  /* interpolation                 */
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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   NP_TRANSFER - type definition for grid transfer

   DESCRIPTION:
   This numproc type is used for the description of grid transfers.
   It can be called by the given interface from a muligrid solver.
   Initializing the data is optional; it can be done with

   'INT NPTransferInit (NP_TRANSFER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPTransferDisplay (NP_TRANSFER *theNP);'
   'INT NPTransferExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb

   struct np_transfer {

        NP_BASE base;                        // inherits base class

        // data (optinal, necessary for calling the generic execute routine)
    VECDATA_DESC *x;                     // solution
    VECDATA_DESC *c;                     // correction
    VECDATA_DESC *b;                     // defect
    MATDATA_DESC *A;                     // matrix
        VEC_SCALAR damp;                     // damping factor

        // functions
        INT (*PreProcess)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*PreProcessSolution)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
        INT (*PreProcessProject)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  INT *);                        // result
    INT (*InterpolateCorrection)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // destination vector
                  VECDATA_DESC *,                // source vector
                  MATDATA_DESC *,                // matrix
                  VEC_SCALAR,                    // damping factor
                  INT *);                        // result
    INT (*RestrictDefect)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // destination vector
                  VECDATA_DESC *,                // source vector
                  MATDATA_DESC *,                // matrix
                  VEC_SCALAR,                    // damping factor
                  INT *);                        // result
    INT (*InterpolateNewVectors)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
    INT (*ProjectSolution)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
    INT (*AdaptCorrection)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // correction vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*PostProcess)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*PostProcessProject)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  INT *);                        // result
        INT (*PostProcessSolution)
             (struct np_transfer *,          // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
   };
   typedef struct np_transfer NP_TRANSFER;
   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPTransferInit (NP_TRANSFER *np, INT argc , char **argv)
{
  INT i;

  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->c = ReadArgvVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);
  np->baselevel = 0;
  ReadArgvINT("baselevel",&(np->baselevel),argc,argv);

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
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (sc_disp(np->damp,np->b,"damp"))
    return (1);

  return(0);
}

INT NPTransferExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_TRANSFER *np;
  INT result,level;

  np = (NP_TRANSFER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (ReadArgvOption("i",argc,argv)) {
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
    if ((*np->PreProcess)(np,&(np->baselevel),level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPTransferExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->PreProcessSolution == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no PreProcessSolution");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if ((*np->PreProcessSolution)(np,np->baselevel,level,np->x,&result)) {
      UserWriteF("NPTransferExecute: PreProcessSolution failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("R",argc,argv)) {
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

  if (ReadArgvOption("I",argc,argv)) {
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

  if (ReadArgvOption("N",argc,argv)) {
    if (np->InterpolateNewVectors == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no InterpolateNewVectors");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if ((*np->InterpolateNewVectors)(np,0,level,np->x,&result)) {
      UserWriteF("NPTransferExecute: InterpolateNewVectors failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("P",argc,argv)) {
    if (np->ProjectSolution == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no ProjectSolution");
      return (1);
    }
    if (np->x == NULL) {
      PrintErrorMessage('E',"NPTransferExecute","no vector x");
      return (1);
    }
    if ((*np->ProjectSolution)(np,0,level,np->x,&result)) {
      UserWriteF("NPTransferExecute: ProjectSolution failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
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
    if ((*np->PostProcess)(np,&(np->baselevel),level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPTransferExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  return(0);
}

INT MinimizeLevel (GRID *theGrid, VECDATA_DESC *c, VECDATA_DESC *b,
                   MATDATA_DESC *A, VECDATA_DESC *t, INT display)
{
  VEC_SCALAR scal;
  DOUBLE a0,a1;
  INT j,ncomp;

  ncomp = VD_NCOMP(c);
  if (l_dset(theGrid,t,EVERY_CLASS,0.0) != NUM_OK)
    return(1);
  if (l_dmatmul(theGrid,t,NEWDEF_CLASS,A,c,ACTIVE_CLASS) != NUM_OK)
    return(1);
    #ifdef ModelP
  if (l_vector_collect(theGrid,t) != NUM_OK)
    return (1);
    #endif
  if (l_ddot (theGrid,t,NEWDEF_CLASS,b,scal) != NUM_OK)
    return(1);
  a0 = 0.0;
  for (j=0; j<ncomp; j++)
    a0 += scal[j];
  if (l_ddot (theGrid,t,NEWDEF_CLASS,t,scal) != NUM_OK)
    return(1);
  a1 = 0.0;
  for (j=0; j<ncomp; j++)
    a1 += scal[j];
  if (a1 <= 0.0)
    return(1);
  if (display == PCR_FULL_DISPLAY)
    UserWriteF("       min  %7.4lf\n",1+a0/a1);
  for (j=0; j<ncomp; j++)
    scal[j] = 1 + a0 / a1;
  if (l_dscale (theGrid,c,ACTIVE_CLASS,scal))
    return(1);
  for (j=0; j<ncomp; j++)
    scal[j] = - a0 / a1;
  if (l_daxpy (theGrid,b,NEWDEF_CLASS,scal,t))
    return(1);

  return(0);
}

/****************************************************************************/
/*D
   transfer - num proc for grid transfer configuration

   DESCRIPTION:
   This numproc configures the grid transfer.
   It can be used in 'lmgc'.

   .vb
   npinit <name> [$x <sol>] [$c <cor>] [$b <rhs>] [$A <mat>]
              [$M] [$L] [$m] [$d {full|red|no}] [$S <cut>];
   .ve

   .  $x~<sol> - solution vector
   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $M - use iterpolation matrix for transfer (default: standard transfer)
   .  $m - parallel nonconforming interpolation
   .  $L - level optimization
   .  $d - display modus
   .  $S~<cut> - scaled restriction with cut (in combination with standard prologation)

   'npexecute <name> [$i] [$R] [$I] [$N] [$P] [$p];'

   .  $i - preprocess
   .  $R - restrict defect
   .  $I - interpolate correction
   .  $N - interpolate new vectors
   .  $P - project solution
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT TransferInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;

  np->mode = STANDARD_MODE;
  np->res = StandardRestrict;
  np->intcor = StandardInterpolateCorrection;
  np->intnew = StandardInterpolateNewVectors;

  if (ReadArgvOption("M",argc,argv)) {
    np->mode = IMAT_MODE;
    np->res = RestrictByMatrix;
    np->intcor = InterpolateCorrectionByMatrix;
    np->intnew = InterpolateNewVectorsByMatrix;
  }
  np->meanvalue = ReadArgvOption("m",argc,argv);
  np->level = ReadArgvOption("L",argc,argv);
  np->display = ReadArgvDisplay(argc,argv);

  if (ReadArgvOption("S",argc,argv))
    if (ReadArgvDOUBLE("S",&(np->cut),argc,argv))
      UserWrite("$S option not active!\n");
    else
    {
      np->mode = SCALEDMG_MODE;
      np->res = ScaledMGRestrict;
      np->intcor = StandardInterpolateCorrection;
      np->intnew = StandardInterpolateNewVectors;
    }
  np->L = ReadArgvMatDesc(theNP->mg,"B",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->amg = (NP_TRANSFER *)
            ReadArgvNumProc(theNP->mg,"amg",TRANSFER_CLASS_NAME,argc,argv);

  return (NPTransferInit(&np->transfer,argc,argv));
}

static INT TransferDisplay (NP_BASE *theNP)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;

  NPTransferDisplay(&np->transfer);

  if (np->res == StandardRestrict)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Restrict","StandardRestrict");
  if (np->res == RestrictByMatrix)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Restrict","RestrictByMatrix");
  if (np->intcor == StandardInterpolateCorrection)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateCor",
               "StandardInterpolateCorrection");
  if (np->intcor == InterpolateCorrectionByMatrix)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateCor",
               "InterpolateCorrectionByMatrix");
  if (np->intnew == StandardInterpolateNewVectors)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateNew",
               "StandardInterpolateNewVectors");
  if (np->intnew == InterpolateNewVectorsByMatrix)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"InterpolateNew",
               "InterpolateNewVectorsByMatrix");

  UserWriteF(DISPLAY_NP_FORMAT_SI,"meanvalue",(int)np->meanvalue);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"level",(int)np->level);

  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->L != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(np->L));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->amg != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"amg",ENVITEM_NAME(np->amg));

  return (0);
}

static INT TransferPreProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
                               VECDATA_DESC *x, VECDATA_DESC *b,
                               MATDATA_DESC *A, INT *result)
{
  NP_STANDARD_TRANSFER *np;
  MULTIGRID *theMG;
  INT i,err;
#       ifdef ModelP
  GRID *theGrid;
#       endif

  np = (NP_STANDARD_TRANSFER *) theNP;
  theMG = theNP->base.mg;

  if (np->amg != NULL) {
    if (*fl == 0)
      if (np->amg->PreProcess(np->amg,fl,0,x,b,A,result))
        REP_ERR_RETURN(1);
        #ifdef ModelP
    if (a_vector_vecskip(theMG,*fl,tl,x) != NUM_OK)
      NP_RETURN(1,result[0]);
        #endif
    for (i=*fl; i<=tl; i++)
      if (AssembleDirichletBoundary (GRID_ON_LEVEL(theMG,i),A,x,b))
        NP_RETURN(1,result[0]);
    if (np->display != PCR_NO_DISPLAY)
      UserWrite(" [d]\n");
  }

  if (np->mode == SCALEDMG_MODE)
  {
                #ifdef ModelP
    if (AllocMDFromMD(theNP->base.mg,*fl,tl,A,&np->L)) {
      result[0] = __LINE__;
      return (1);
    }
    for (i=tl; i>=*fl; i--)
    {
      theGrid = GRID_ON_LEVEL(theMG,i);
      if (l_dmatcopy(theGrid,np->L,A) != NUM_OK) {
        result[0] = __LINE__;
        return (1);
      }
      if (l_matrix_consistent(theGrid,np->L,MAT_GHOST_DIAG_CONS) != NUM_OK) {
        result[0] = __LINE__;
        return (1);
      }
    }
                #else
    np->L = A;
                #endif
    /* create restriction matrices */
    for (i=tl; i>*fl; i--)
    {
      err = InstallScaledRestrictionMatrix(GRID_ON_LEVEL(theMG,i),np->L,np->cut);
      if (err!=NUM_OK) {
        UserWriteF("InstallScaledRestrictionMatrix failed in %d\n",err);
        result[0] = __LINE__;
        return(1);
      }
    }
    /* scale equations */
    for (i=tl; i>=*fl; i--)
      if (DiagonalScaleSystem(GRID_ON_LEVEL(theMG,i),A,np->L,b)!=NUM_OK) {
        result[0] = __LINE__;
        return (1);
      }
                #ifdef ModelP
    FreeMD(theNP->base.mg,*fl,tl,np->L);
                #endif
  }

  return(0);
}

static INT RestrictDefect (NP_TRANSFER *theNP, INT level,
                           VECDATA_DESC *to, VECDATA_DESC *from,
                           MATDATA_DESC *A, VEC_SCALAR damp,
                           INT *result)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;
  if (level < 1)
    result[0] = RestrictByMatrix(GRID_ON_LEVEL(theNP->base.mg,level),
                                 to,from,damp);
  else
    result[0] = (*np->res)(GRID_ON_LEVEL(theNP->base.mg,level),
                           to,from,damp);

    #ifdef ModelP
  if (l_ghostvector_collect(GRID_ON_LEVEL(theNP->base.mg,level-1),to)
      != NUM_OK) NP_RETURN(1,result[0]);
        #endif

  return(result[0]);
}

static INT InterpolateCorrection (NP_TRANSFER *theNP, INT level,
                                  VECDATA_DESC *to, VECDATA_DESC *from,
                                  MATDATA_DESC *A, VEC_SCALAR damp,
                                  INT *result)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;
    #ifdef ModelP
  if (l_ghostvector_consistent(GRID_ON_LEVEL(theNP->base.mg,level-1),from)
      != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (level < 1)
    result[0] =
      InterpolateCorrectionByMatrix(GRID_ON_LEVEL(theNP->base.mg,level),
                                    to,from,damp);
  else
    result[0] =
      (*np->intcor)(GRID_ON_LEVEL(theNP->base.mg,level),to,from,damp);
    #ifdef ModelP
  if (np->meanvalue)
    if (l_vector_meanvalue(GRID_ON_LEVEL(theNP->base.mg,level),to)
        != NUM_OK) NP_RETURN(1,result[0]);
        #endif

  return(result[0]);
}

static INT InterpolateNewVectors (NP_TRANSFER *theNP,  INT fl, INT tl,
                                  VECDATA_DESC *x, INT *result)
{
  NP_STANDARD_TRANSFER *np;
  INT i;

  np = (NP_STANDARD_TRANSFER *) theNP;
  for (i=fl+1; i<=tl; i++) {
        #ifdef ModelP
    if (l_ghostvector_consistent(GRID_ON_LEVEL(theNP->base.mg,i-1),x)
        != NUM_OK) NP_RETURN(1,result[0]);
        #endif
    result[0] = (*np->intnew)(GRID_ON_LEVEL(theNP->base.mg,i),x);
    if (result[0]) NP_RETURN(1,result[0]);
  }

  return(0);
}

static INT ProjectSolution (NP_TRANSFER *theNP,  INT fl, INT tl,
                            VECDATA_DESC *x, INT *result)
{
  INT i;

  result[0] = 0;
  for (i=tl-1; i>=fl; i--) {
    result[0] = StandardProject(GRID_ON_LEVEL(theNP->base.mg,i),x,x);
    if (result[0]) NP_RETURN(1,result[0]);
        #ifdef ModelP
    if (l_ghostvector_project(GRID_ON_LEVEL(theNP->base.mg,i),x)
        != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  }

  return(0);
}

static INT AdaptCorrection (NP_TRANSFER *theNP, INT level,
                            VECDATA_DESC *c, VECDATA_DESC *b,
                            MATDATA_DESC *A, INT *result)
{
  NP_STANDARD_TRANSFER *np;
  GRID *theGrid;

  np = (NP_STANDARD_TRANSFER *) theNP;
  if (np->level) {
    theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
    if (AllocVDFromVD(theNP->base.mg,level,level,c,&np->t)) NP_RETURN(1,result[0]);
    if (MinimizeLevel(theGrid,c,b,A,np->t,np->display)) NP_RETURN(1,result[0]);
    if (FreeVD(theNP->base.mg,level,level,np->t)) NP_RETURN(1,result[0]);
  }

  return(0);
}

static INT TransferPostProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
                                VECDATA_DESC *x, VECDATA_DESC *b,
                                MATDATA_DESC *A, INT *result)
{
  NP_STANDARD_TRANSFER *np;

  np = (NP_STANDARD_TRANSFER *) theNP;
  if (np->amg != NULL)
    if (np->amg->PostProcess(np->amg,fl,0,x,b,A,result))
      REP_ERR_RETURN(1);

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
  np->PreProcessProject = NULL;
  np->PreProcessSolution = NULL;
  np->RestrictDefect = RestrictDefect;
  np->InterpolateCorrection = InterpolateCorrection;
  np->InterpolateNewVectors = InterpolateNewVectors;
  np->ProjectSolution = ProjectSolution;
  np->AdaptCorrection = AdaptCorrection;
  np->PostProcess = TransferPostProcess;
  np->PostProcessProject = NULL;
  np->PostProcessSolution = NULL;

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
  if (CreateClass(TRANSFER_CLASS_NAME ".transfer",
                  sizeof(NP_STANDARD_TRANSFER),TransferConstruct))
    return (__LINE__);

  return (0);
}
