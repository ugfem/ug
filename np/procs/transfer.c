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
#include "formats.h"
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

#define MAX_PT                          2
#define MAX_VD                          5
#define PT_NOT_INIT                     -1

enum {
  PT_PRE_x,
  PT_PRE_b,

  PT_PRESOL_x,

  PT_REST_f,
  PT_REST_t,

  PT_ICOR_f,
  PT_ICOR_t,

  PT_INEW_x,

  PT_PROJSOL_x,

  PT_ACOR_c,
  PT_ACOR_b,

  PT_POSTSOL_x,

  PT_POST_x,
  PT_POST_b,

  PT_N_IND
};

#define PT_MVT(pt)                              ((pt)->mvt)
#define PT_NTRANS(pt)                   ((pt)->ntrans)
#define PT_TRANS(pt,i)                  ((pt)->trans[i])
#define PT_SVT(pt,i)                    ((pt)->svt[i])
#define PT_SMD(pt,i)                    ((pt)->smd[i])
#define PT_SMDI(pt,i)                   ((pt)->smdi[i])
#define PT_NVD(pt)                              ((pt)->nvd)
#define PT_MVD(pt,i)                    ((pt)->mvd[i])
#define PT_IND(pt,i)                    ((pt)->ind[i])
#define PT_SVD(pt,k,i)                  ((pt)->svd[k][i])
#define PT_SVDI(pt,k,i)                 ((pt)->svdi[k][i])

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

typedef struct
{
  NP_TRANSFER transfer;

  /* additional data */
  const VEC_TEMPLATE *mvt;                      /* main vector template					*/
  INT ntrans;                                                   /* number of part transfer numprocs		*/
  NP_TRANSFER *trans[MAX_PT];                   /* pointers to part transfer numprocs	*/
  INT svt[MAX_PT];                                      /* indices of sub vector templates		*/

  /* sub XXXDATA_DESCs needed */
  MATDATA_DESC *smd[MAX_PT];                    /* sub descriptors of A					*/
  MATDATA_DESC *smdi[MAX_PT];                   /* interface descs of A					*/

  INT nvd;                                                      /* number of mains with constructed subs*/
  const VECDATA_DESC *mvd[MAX_VD];      /* main VECDATA_DESCs					*/
  INT ind[PT_N_IND];                                    /* indices to sub descriptors of mains	*/
  VECDATA_DESC *svd[MAX_VD][MAX_PT];       /* sub descriptors of mains			*/
  VECDATA_DESC *svdi[MAX_VD][MAX_PT];      /* interface descs of mains			*/

} NP_PART_TRANSFER;

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

  if (sc_read(np->damp,NP_FMT(np),np->x,"damp",argc,argv))
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
/*D
   parttransfer - num proc for grid transfer with domain parts

   DESCRIPTION:
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

static INT PartTransferInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  FORMAT *fmt;
  VEC_TEMPLATE *mvt;
  INT i,ntrans,nsub;
  char name[NAMESIZE],buffer[VALUELEN];

  fmt = NP_FMT(theNP);

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get name of main vector template */
  if (ReadArgvChar("m",buffer,argc,argv)!=NULL)
  {
    PrintErrorMessage('E',"PartTransferInit","m option with main vector template not found");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  mvt = GetVectorTemplate(fmt,buffer);
  if (mvt == NULL)
  {
    PrintErrorMessageF('E',"PartTransferInit",
                       "cannot find specified vector template '%s'",buffer);
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PT_MVT(thePT) = mvt;

  PT_NTRANS(thePT) = ntrans = nsub = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      /* scan name of subvector of main vector template */
      if (sscanf(argv[i],expandfmt(CONCAT3("s %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"PartTransferInit","specify a sub vector template with $s");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      for (i=0; i<VF_NSUB(mvt); i++)
        if (strcmp(SUBV_NAME(VF_SUB(mvt,i)),name)==0)
          break;
      if (i>=VF_NSUB(mvt))
      {
        PrintErrorMessageF('E',"PartTransferInit","name '%s' of sub template not found",name);
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PT_SVT(thePT,nsub++) = i;
      break;

    case 't' :
      /* scan transfer numproc */
      if (ntrans>=MAX_PT)
      {
        PrintErrorMessage('E',"PartTransferInit","max number of part transfer numprocs exceeded");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      if (sscanf(argv[i],expandfmt(CONCAT3("trans %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"PartTransferInit","specify a transfer numproc with $trans");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      trans = (NP_TRANSFER*) GetNumProcByName (NP_MG(theNP),name,TRANSFER_CLASS_NAME);
      if (trans == NULL)
      {
        PrintErrorMessage('E',"PartTransferInit",
                          "cannot find specified numerical procedure");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PT_TRANS(thePT,ntrans++) = trans;
      break;
    }

  if (ntrans==0)
  {
    PrintErrorMessage('E',"PartTransferInit","specify at least one transfer numproc with $trans");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  if (nsub!=ntrans)
  {
    PrintErrorMessage('E',"PartTransferInit","number of sub templates and transfer numprocs have to be equal");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PT_NTRANS(thePT) = ntrans;

  return (NPTransferInit(&thePT->transfer,argc,argv));
}

static INT PartTransferDisplay (NP_BASE *theNP)
{
  NP_PART_TRANSFER *thePT;
  const VEC_TEMPLATE *mvt;
  INT i;
  char text1[16],text2[64];

  thePT = (NP_PART_TRANSFER *) theNP;

  NPTransferDisplay(&thePT->transfer);

  mvt = PT_MVT(thePT);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"main tplt",ENVITEM_NAME(mvt));

  UserWrite("\npart transfer numprocs and sub templates:\n");
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    sprintf(text1,"trans+sub%d",i);
    sprintf(text2,"%20s %20s",ENVITEM_NAME(PT_TRANS(thePT,i)),SUBV_NAME(VF_SUB(mvt,PT_SVT(thePT,i))));
    UserWriteF(DISPLAY_NP_FORMAT_SS,text1,text2);
  }

  return (0);
}

static INT ConstructSubVD (NP_PART_TRANSFER *thePT, const VECDATA_DESC *x, INT *index)
{
  const VEC_TEMPLATE *mvt;
  INT i,k,n;

  if (*index==PT_NOT_INIT)
  {
    /* first check existing ones */
    for (k=0; k<PT_NVD(thePT); k++)
      if (VDequal(x,PT_MVD(thePT,k)))
        break;
    if (k<PT_NVD(thePT))
      *index = k;
  }

  if (*index!=PT_NOT_INIT)
    return (0);

  /* x is new: store pointer to x and decompose it into sub descriptors */
  n = PT_NTRANS(thePT);
  k = PT_NVD(thePT);
  ASSERT(k<MAX_VD);
  PT_MVD(thePT,k) = x;
  mvt = PT_MVT(thePT);
  for (i=0; i<n; i++)
  {
    if (VDsubDescFromVT(x,mvt,PT_SVT(thePT,i),&PT_SVD(thePT,k,i)))
      REP_ERR_RETURN(1);
    if (VDinterfaceDesc(x,PT_SVD(thePT,k,i),&PT_SVDI(thePT,k,i)))
      REP_ERR_RETURN(1);
  }
  PT_NVD(thePT)++;

  *index = k;

  return (0);
}

static INT ConstructSubMD (NP_PART_TRANSFER *thePT, const MATDATA_DESC *A)
{
  const VEC_TEMPLATE *mvt;
  INT i,n;

  n = PT_NTRANS(thePT);
  mvt = PT_MVT(thePT);
  for (i=0; i<n; i++)
  {
    if (MDsubDescFromVT(A,mvt,PT_SVT(thePT,i),&PT_SMD(thePT,i)))
      REP_ERR_RETURN(1);
    if (MDinterfaceDesc(A,PT_SMD(thePT,i),&PT_SMDI(thePT,i)))
      REP_ERR_RETURN(1);
  }
  return (0);
}

static INT PartTransferPreProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
                                   VECDATA_DESC *x, VECDATA_DESC *b,
                                   MATDATA_DESC *A, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* indicate that decomposition into sub descriptors has to be done */
  for (i=0; i<PT_N_IND; i++)
    PT_IND(thePT,i) = PT_NOT_INIT;
  PT_NVD(thePT) = 0;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,x,&PT_IND(thePT,PT_PRE_x)))
    REP_ERR_RETURN(1);
  if (ConstructSubVD(thePT,b,&PT_IND(thePT,PT_PRE_b)))
    REP_ERR_RETURN(1);

  /* get the MATDATA_DESC needed (once and for all) */
  if (ConstructSubMD(thePT,A))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 2;
  SPID_NMD(&spid) = 1;

  /* call pre of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_PRE(trans)!=NULL)
    {
      SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_PRE_x),i);
      SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_PRE_x),i);
      SPID_VD(&spid,1)  = PT_SVD(thePT,PT_IND(thePT,PT_PRE_b),i);
      SPID_VDI(&spid,1) = PT_SVDI(thePT,PT_IND(thePT,PT_PRE_b),i);
      SPID_MD(&spid,0)  = PT_SMD(thePT,i);
      SPID_MDI(&spid,0) = PT_SMDI(thePT,i);
      if (SwapPartInterfaceData(*fl,tl,&spid,SPID_FORTH))
        REP_ERR_RETURN(1);

      if (NPTR_PRE(trans) (trans,fl,tl,
                           PT_SVD(thePT,PT_IND(thePT,PT_PRE_x),i),
                           PT_SVD(thePT,PT_IND(thePT,PT_PRE_b),i),
                           PT_SMD(thePT,i),
                           result))
        REP_ERR_RETURN(1);
      if (SwapPartInterfaceData(*fl,tl,&spid,SPID_BACK))
        REP_ERR_RETURN(1);
    }
  }

  return (0);
}

static INT PartPreProcessProject (NP_TRANSFER *theNP,  INT fl, INT tl, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* call preprocess project of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_PREPROJ(trans)==NULL)
      continue;

    if (NPTR_PREPROJ(trans) (trans,fl,tl,result))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartPreProcessSolution (NP_TRANSFER *theNP,  INT fl, INT tl,
                                   VECDATA_DESC *x, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,x,&PT_IND(thePT,PT_PRESOL_x)))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 1;
  SPID_NMD(&spid) = 0;

  /* call preprocess solution of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_PRESOL(trans)==NULL)
      continue;

    /* swap interface data */
    SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_PRESOL_x),i);
    SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_PRESOL_x),i);
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_FORTH))
      REP_ERR_RETURN(1);

    if (NPTR_PRESOL(trans) (trans,fl,tl,
                            PT_SVD(thePT,PT_IND(thePT,PT_PRESOL_x),i),
                            result))
      REP_ERR_RETURN(1);

    /* swap interface data back */
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_BACK))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartRestrictDefect (NP_TRANSFER *theNP, INT level,
                               VECDATA_DESC *to, VECDATA_DESC *from,
                               MATDATA_DESC *A, VEC_SCALAR damp,
                               INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,from,&PT_IND(thePT,PT_REST_f)))
    REP_ERR_RETURN(1);
  if (ConstructSubVD(thePT,to,  &PT_IND(thePT,PT_REST_t)))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 2;
  SPID_NMD(&spid) = 1;

  /* call restrict defect of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    /* swap interface data */
    SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_REST_t),i);
    SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_REST_t),i);
    SPID_VD(&spid,1)  = PT_SVD(thePT,PT_IND(thePT,PT_REST_f),i);
    SPID_VDI(&spid,1) = PT_SVDI(thePT,PT_IND(thePT,PT_REST_f),i);
    SPID_MD(&spid,0)  = PT_SMD(thePT,i);
    SPID_MDI(&spid,0) = PT_SMDI(thePT,i);
    if (SwapPartInterfaceData(level,level,&spid,SPID_FORTH))
      REP_ERR_RETURN(1);

    trans = PT_TRANS(thePT,i);
    if (NPTR_RESTRICT(trans) (trans,level,
                              PT_SVD(thePT,PT_IND(thePT,PT_REST_t),i),
                              PT_SVD(thePT,PT_IND(thePT,PT_REST_f),i),
                              PT_SMD(thePT,i),
                              damp,result))
      REP_ERR_RETURN(1);

    /* swap interface data back */
    if (SwapPartInterfaceData(level,level,&spid,SPID_BACK))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartInterpolateCorrection (NP_TRANSFER *theNP, INT level,
                                      VECDATA_DESC *to, VECDATA_DESC *from,
                                      MATDATA_DESC *A, VEC_SCALAR damp,
                                      INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,from,&PT_IND(thePT,PT_ICOR_f)))
    REP_ERR_RETURN(1);
  if (ConstructSubVD(thePT,to,  &PT_IND(thePT,PT_ICOR_t)))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 2;
  SPID_NMD(&spid) = 1;

  /* call interpolate correction of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    /* swap interface data */
    SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_ICOR_t),i);
    SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_ICOR_t),i);
    SPID_VD(&spid,1)  = PT_SVD(thePT,PT_IND(thePT,PT_ICOR_f),i);
    SPID_VDI(&spid,1) = PT_SVDI(thePT,PT_IND(thePT,PT_ICOR_f),i);
    SPID_MD(&spid,0)  = PT_SMD(thePT,i);
    SPID_MDI(&spid,0) = PT_SMDI(thePT,i);
    if (SwapPartInterfaceData(level,level,&spid,SPID_FORTH))
      REP_ERR_RETURN(1);

    trans = PT_TRANS(thePT,i);
    if (NPTR_INTCOR(trans) (trans,level,
                            PT_SVD(thePT,PT_IND(thePT,PT_ICOR_t),i),
                            PT_SVD(thePT,PT_IND(thePT,PT_ICOR_f),i),
                            PT_SMD(thePT,i),
                            damp,result))
      REP_ERR_RETURN(1);

    /* swap interface data back */
    if (SwapPartInterfaceData(level,level,&spid,SPID_BACK))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartInterpolateNewVectors (NP_TRANSFER *theNP,  INT fl, INT tl,
                                      VECDATA_DESC *x, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,x,&PT_IND(thePT,PT_INEW_x)))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 1;
  SPID_NMD(&spid) = 0;

  /* call interpolate new of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    /* swap interface data */
    SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_INEW_x),i);
    SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_INEW_x),i);
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_FORTH))
      REP_ERR_RETURN(1);

    trans = PT_TRANS(thePT,i);
    if (NPTR_INTNEW(trans) (trans,fl,tl,
                            PT_SVD(thePT,PT_IND(thePT,PT_INEW_x),i),
                            result))
      REP_ERR_RETURN(1);

    /* swap interface data back */
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_BACK))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartProjectSolution (NP_TRANSFER *theNP,  INT fl, INT tl,
                                VECDATA_DESC *x, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,x,&PT_IND(thePT,PT_PROJSOL_x)))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 1;
  SPID_NMD(&spid) = 0;

  /* call project solution of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_PROJSOL(trans)==NULL)
      continue;

    /* swap interface data */
    SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_PROJSOL_x),i);
    SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_PROJSOL_x),i);
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_FORTH))
      REP_ERR_RETURN(1);

    if (NPTR_PROJSOL(trans) (trans,fl,tl,
                             PT_SVD(thePT,PT_IND(thePT,PT_PROJSOL_x),i),
                             result))
      REP_ERR_RETURN(1);

    /* swap interface data back */
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_BACK))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartAdaptCorrection (NP_TRANSFER *theNP, INT level,
                                VECDATA_DESC *c, VECDATA_DESC *b,
                                MATDATA_DESC *A, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,c,&PT_IND(thePT,PT_ACOR_c)))
    REP_ERR_RETURN(1);
  if (ConstructSubVD(thePT,b,&PT_IND(thePT,PT_ACOR_b)))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 2;
  SPID_NMD(&spid) = 1;

  /* call adapt correction of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_ADPTCOR(trans)==NULL)
      continue;

    /* swap interface data */
    SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_ACOR_c),i);
    SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_ACOR_c),i);
    SPID_VD(&spid,1)  = PT_SVD(thePT,PT_IND(thePT,PT_ACOR_b),i);
    SPID_VDI(&spid,1) = PT_SVDI(thePT,PT_IND(thePT,PT_ACOR_b),i);
    SPID_MD(&spid,0)  = PT_SMD(thePT,i);
    SPID_MDI(&spid,0) = PT_SMDI(thePT,i);
    if (SwapPartInterfaceData(level,level,&spid,SPID_FORTH))
      REP_ERR_RETURN(1);

    if (NPTR_ADPTCOR(trans) (trans,level,
                             PT_SVD(thePT,PT_IND(thePT,PT_ACOR_c),i),
                             PT_SVD(thePT,PT_IND(thePT,PT_ACOR_b),i),
                             PT_SMD(thePT,i),
                             result))
      REP_ERR_RETURN(1);

    /* swap interface data back */
    if (SwapPartInterfaceData(level,level,&spid,SPID_BACK))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartPostProcessProject (NP_TRANSFER *theNP,  INT *fl, INT tl, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* call postprocess project of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_POSTPROJ(trans)==NULL)
      continue;

    if (NPTR_POSTPROJ(trans) (trans,fl,tl,result))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartPostProcessSolution (NP_TRANSFER *theNP,  INT fl, INT tl,
                                    VECDATA_DESC *x, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  SPID_DESC spid;
  INT i;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,x,&PT_IND(thePT,PT_POSTSOL_x)))
    REP_ERR_RETURN(1);

  SPID_NVD(&spid) = 1;
  SPID_NMD(&spid) = 0;

  /* call postprocess solution of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_POSTSOL(trans)==NULL)
      continue;

    /* swap interface data */
    SPID_VD(&spid,0)  = PT_SVD(thePT,PT_IND(thePT,PT_POSTSOL_x),i);
    SPID_VDI(&spid,0) = PT_SVDI(thePT,PT_IND(thePT,PT_POSTSOL_x),i);
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_FORTH))
      REP_ERR_RETURN(1);

    if (NPTR_POSTSOL(trans) (trans,fl,tl,
                             PT_SVD(thePT,PT_IND(thePT,PT_POSTSOL_x),i),
                             result))
      REP_ERR_RETURN(1);

    /* swap interface data back */
    if (SwapPartInterfaceData(fl,tl,&spid,SPID_BACK))
      REP_ERR_RETURN(1);
  }

  return (0);
}

static INT PartTransferPostProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
                                    VECDATA_DESC *x, VECDATA_DESC *b,
                                    MATDATA_DESC *A, INT *result)
{
  NP_PART_TRANSFER *thePT;
  NP_TRANSFER *trans;
  INT i,k;

  thePT = (NP_PART_TRANSFER *) theNP;

  /* get the VECDATA_DESCs needed */
  if (ConstructSubVD(thePT,x,&PT_IND(thePT,PT_POST_x)))
    REP_ERR_RETURN(1);
  if (ConstructSubVD(thePT,b,&PT_IND(thePT,PT_POST_b)))
    REP_ERR_RETURN(1);

  /* call post of part transfer numprocs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    trans = PT_TRANS(thePT,i);
    if (NPTR_POST(trans)!=NULL)
      if (NPTR_POST(trans) (trans,fl,tl,
                            PT_SVD(thePT,PT_IND(thePT,PT_POST_x),i),
                            PT_SVD(thePT,PT_IND(thePT,PT_POST_b),i),
                            PT_SMD(thePT,i),
                            result))
        REP_ERR_RETURN(1);
  }

  /* dispose the auxiliary XXXDATA_DESCs */
  for (i=0; i<PT_NTRANS(thePT); i++)
  {
    if (DisposeMD(PT_SMD(thePT,i)))
      REP_ERR_RETURN(1);
    if (DisposeMD(PT_SMDI(thePT,i)))
      REP_ERR_RETURN(1);
    for (k=0; k<PT_NVD(thePT); k++)
    {
      if (DisposeVD(PT_SVD(thePT,k,i)))
        REP_ERR_RETURN(1);
      if (DisposeVD(PT_SVDI(thePT,k,i)))
        REP_ERR_RETURN(1);
    }
  }

  return (0);
}

static INT PartTransferConstruct (NP_BASE *theNP)
{
  NP_TRANSFER *np;

  theNP->Init                                     = PartTransferInit;
  theNP->Display                          = PartTransferDisplay;
  theNP->Execute                          = NPTransferExecute;

  np = (NP_TRANSFER *) theNP;

  np->PreProcess                          = PartTransferPreProcess;
  np->PreProcessProject           = PartPreProcessProject;
  np->PreProcessSolution          = PartPreProcessSolution;
  np->RestrictDefect                      = PartRestrictDefect;
  np->InterpolateCorrection       = PartInterpolateCorrection;
  np->InterpolateNewVectors       = PartInterpolateNewVectors;
  np->ProjectSolution                     = PartProjectSolution;
  np->AdaptCorrection                     = PartAdaptCorrection;
  np->PostProcessProject          = PartPostProcessProject;
  np->PostProcessSolution         = PartPostProcessSolution;
  np->PostProcess                         = PartTransferPostProcess;

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

  if (CreateClass(TRANSFER_CLASS_NAME ".parttransfer",
                  sizeof(NP_PART_TRANSFER),PartTransferConstruct))
    return (__LINE__);

  return (0);
}
