// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ts.c                                                              */
/*                                                                          */
/* Purpose:   time-stepping scheme for (non-)linear time-dependent problems */
/*                                                                          */
/* Author:    Peter Bastian                                                                                             */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   January  09, 1997  begin                                      */
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "devices.h"
#include "general.h"
#include "gm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"

#include "nls.h"
#include "ls.h"
#include "nls.h"
#include "assemble.h"
#include "transfer.h"
#include "ts.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/*  Implement default functions                                             */
/*                                                                          */
/****************************************************************************/

INT NPTSolverInit (NP_T_SOLVER *np, INT argc , char **argv)
{
  INT i,r;

  r = NP_EXECUTABLE;       /* highest state */

  /* The solution is required for execution */
  np->y = ReadArgvVecDesc(np->nlass.base.mg,"y",argc,argv);
  if (np->y == NULL)
  {
    r = NP_NOT_ACTIVE;
    UserWrite("Warning: solution y is required for execution !\n");
  }

  /* assemble numproc is required for execution */
  np->tass = (NP_T_ASSEMBLE *)
             ReadArgvNumProc(np->nlass.base.mg,"A",T_ASSEMBLE_CLASS_NAME,argc,argv);
  if (np->tass == NULL) r = NP_NOT_ACTIVE;

  /* solver numproc is required for execution */
  np->nlsolve = (NP_NL_SOLVER *)
                ReadArgvNumProc(np->nlass.base.mg,"S",NL_SOLVER_CLASS_NAME,argc,argv);
  if (np->nlsolve == NULL) r = NP_NOT_ACTIVE;

  return(r);
}

INT NPTSolverDisplay (NP_T_SOLVER *np)
{
  UserWrite("symbolic user data:\n");
  if (np->y != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"y",ENVITEM_NAME(np->y));
  UserWrite("\n");

  UserWrite("configuration parameters:\n");
  if (np->y != NULL)
  {
    if (sc_disp(np->reduction,np->y,"reduction")) return (1);
    if (sc_disp(np->abslimit,np->y,"abslimit")) return (1);
  }
  if (np->tass != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"tass",ENVITEM_NAME(np->tass));
  if (np->nlsolve != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"nlsolve",ENVITEM_NAME(np->nlsolve));

  return(0);
}

INT NPTSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_T_SOLVER *np;
  INT result,level;

  np = (NP_T_SOLVER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->y == NULL) {
    PrintErrorMessage('E',"NPTSolverExecute","no vector y");
    return (1);
  }
  if (np->tass == NULL) {
    PrintErrorMessage('E',"NPTSolverExecute","no assemble num proc");
    return (1);
  }
  if (np->nlsolve == NULL) {
    PrintErrorMessage('E',"NPTSolverExecute","no solver num proc");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (*np->TimePreProcess != NULL)
      if ((*np->TimePreProcess)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimePreProcess failed, error code %d\n",
                   result);
        return (1);
      }
  }

  if (ReadArgvOption("0",argc,argv)) {
    if (*np->TimeInit != NULL)
      if ((*np->TimeInit)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeInit failed, error code %d\n",
                   result);
        return (1);
      }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (*np->TimeStep != NULL)
      if ((*np->TimeStep)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimeStep failed, error code %d\n",
                   result);
        return (1);
      }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (*np->TimePostProcess != NULL)
      if ((*np->TimePostProcess)(np,level,&result)) {
        UserWriteF("NPTSolverExecute: TimePostProcess failed, error code %d\n",
                   result);
        return (1);
      }
  }

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  Init                                                                                        */
/*                                                                          */
/****************************************************************************/

INT InitTSolver (void)
{
  return(0);
}
