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
  INT r;

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
/*D
   NP_T_SOLVER - type definition for time solvers

   DESCRIPTION:
   This numproc type is used for the description of time solvers.
   It can be calls a nonlinear solver of type 'NP_NL_SOLVER'.
   Therefore, it constructs a nonlinear assemble numpro
   of type 'NP_NL_ASSEMBLE' using an assemble numproc of type
   'NP_T_ASSEMBLE'.
   Initializing the data is optional; it can be done with

   'INT NPTSolverInit (NP_T_SOLVER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPTSolverDisplay (NP_T_SOLVER *theNP);'
   'INT NPTSolverExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   struct np_t_solver {

        NP_NL_ASSEMBLE nlass;                // derived from nonlinear assemble

        // things to be initialized by generic init
    VECDATA_DESC *y;                     // solution vector
    NP_T_ASSEMBLE *tass;                 // time assemble numproc
    NP_NL_SOLVER *nlsolve;               // nonlinear solver numproc
        VEC_SCALAR reduction;                // reduction factor per time step
        VEC_SCALAR abslimit;                 // absolute limit for the defect

        // functions
        INT (*TimePreProcess)                            // called before first time step
             (struct np_t_solver *,          // pointer to (derived) object
                  INT,                           // level
                  INT *);                        // result
    INT (*TimeInit)                      // initialize, set initial values
             (struct np_t_solver *,          // pointer to (derived) object
                  INT,                           // level
                  INT *);                        // result
    INT (*TimeStep)                      // b := b - Ax
             (struct np_t_solver *,          // pointer to (derived) object
                  INT,                           // level
                  INT *);                        // result
        INT (*TimePostProcess)                           // to be called after last timestep
             (struct np_t_solver *,          // pointer to (derived) object
                  INT,                           // level
                  INT *);                        // result
   };
   typedef struct np_t_solver NP_T_SOLVER;
   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/*  Init                                                                                        */
/*                                                                          */
/****************************************************************************/

INT InitTSolver (void)
{
  return(0);
}
