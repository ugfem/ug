// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  nls.c	                                                                                                        */
/*																			*/
/* Purpose:   nonlinear solver num procs                                        */
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

#include "scan.h"
#include "devices.h"
#include "numproc.h"
#include "scan.h"
#include "np.h"
#include "general.h"

#include "assemble.h"
#include "transfer.h"
#include "ls.h"
#include "nls.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define ABS_LIMIT 1e-10

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   NP_NL_SOLVER - type definition for nonlinear solvers

   DESCRIPTION:
   This numproc type is used for the description of nonlinear solvers.
   It can be called by the given interface from a nonlinear time solver.
   Initializing the data is optional; it can be done with

   'INT NPNLSolverInit (NP_NL_SOLVER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPNLSolverDisplay (NP_NL_SOLVER *theNP);'
   'INT NPNLSolverExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   typedef struct {
        INT error_code;                     // error code
        INT converged;                      // error code
        VEC_SCALAR first_defect;            // first defect
        VEC_SCALAR last_defect;             // last defect
        INT number_of_nonlinear_iterations; // number of iterations
        INT number_of_line_searches;        // number of line search steps
        INT rho_first;                      // first rho
    INT total_linear_iterations;        // total number
    INT max_linear_iterations;          // max number of linear iterations
    DOUBLE exec_time;                   // for this nonlinear solve ...
   } NLRESULT;

   struct np_nl_solver {
        NP_BASE base;                        // inherits base class

        // data (optinal, necessary for calling the generic execute routine)
    VECDATA_DESC *x;                     // solution
    NP_NL_ASSEMBLE *Assemble;            // the assemble numproc
        VEC_SCALAR reduction;                // reduction factor
        VEC_SCALAR abslimit;                 // absolute limit for the defect

        // functions
        INT (*PreProcess)
             (struct np_nl_solver *,         // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
    INT (*Solver)                        // b := b - Ax
             (struct np_nl_solver *,         // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
          NP_NL_ASSEMBLE *,			     // the assemble numproc
                  VEC_SCALAR,                    // absolute limit for the defect
                  VEC_SCALAR,                    // reduction factor
                  NLRESULT *);                   // result structure
        INT (*PostProcess)
             (struct np_nl_solver *,         // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
   };
   typedef struct np_nl_solver NP_NL_SOLVER;
   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPNLSolverInit (NP_NL_SOLVER *np, INT argc , char **argv)
{
  INT i,r;

  r = NP_EXECUTABLE;       /* highest state */

  /* The solution is required for execution */
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  if (np->x == NULL) r = NP_ACTIVE;

  /* abslimit is required for execution */
  if (sc_read(np->abslimit,NP_FMT(np),np->x,"abslimit",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++) np->abslimit[i] = 1.0E-10;             /* default */

  /* reduction factor is required for execution */
  if (sc_read(np->reduction,NP_FMT(np),NULL,"red",argc,argv))
    r = NP_ACTIVE;

  /* assemble numproc is required for execution */
  np->Assemble = (NP_NL_ASSEMBLE *)
                 ReadArgvNumProc(np->base.mg,"A",NL_ASSEMBLE_CLASS_NAME,argc,argv);
  if (np->Assemble == NULL) r = NP_ACTIVE;

  return(r);
}

INT NPNLSolverDisplay (NP_NL_SOLVER *np)
{
  UserWrite("symbolic user data:\n");
  if (np->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
  UserWrite("\n");

  UserWrite("configuration parameters:\n");
  if (np->x != NULL)
    if (sc_disp(np->reduction,np->x,"red"))
      return (1);
  if (sc_disp(np->abslimit,np->x,"abslimit"))
    return (1);
  if (np->Assemble != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Assemble",ENVITEM_NAME(np->Assemble));

  return(0);
}

INT NPNLSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_NL_SOLVER *np;
  NLRESULT nlresult;
  INT result,level;

  np = (NP_NL_SOLVER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"NPNLSolverExecute","no vector x");
    return (1);
  }
  if (np->Assemble == NULL) {
    PrintErrorMessage('E',"NPNLSolverExecute","no assemble num proc");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess != NULL)
      if ((*np->PreProcess)(np,level,np->x,&result)) {
        UserWriteF("NPNLSolverExecute: PreProcess failed, error code %d\n",
                   result);
        return (1);
      }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->Solver == NULL) {
      PrintErrorMessage('E',"NPNLSolverExecute","no Solver");
      return (1);
    }
    if ((*np->Solver)(np,level,np->x,np->Assemble,
                      np->abslimit,np->reduction,&nlresult)) {
      UserWriteF("NPNLSolverExecute: Solver failed, error code %d\n",
                 nlresult.error_code);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess != NULL)
      if ((*np->PostProcess)(np,level,np->x,&result)) {
        UserWriteF("NPNLSolverExecute: PostProcess failed, error code %d\n",
                   result);
        return (1);
      }
  }

  return(0);
}


/****************************************************************************/
/*
   InitNonlinearSolver	- Init this file

   SYNOPSIS:
   INT InitNonlinearSolver ();

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

INT InitNonlinearSolver ()
{
  return (0);
}
