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

/* derive a final class from NP_NL_SOLVER */
typedef struct
{
  NP_NL_SOLVER nlsolver;                /* derived from abstract class NP_NL_SOLVER	*/

  /* parameters to be set via npinit */
  INT DisplayMode;                              /* for PCR									*/
  INT maxit;                                            /* maximum number of newton iterations		*/
  NP_LINEAR_SOLVER *Solve;              /* uses linear solver						*/
  NP_TRANSFER *Trans;                       /* uses transgrid	                                        */
  DOUBLE rho_reass;                             /* reassemble if nonlin conv worth than this*/
  INT linearRate;                               /* 1 if nonquadratic nonlin rate assumed        */
  DOUBLE lin_min_red;                           /* minimum reduction for linear solver		*/
  INT LineSearch;                               /* do line search                                                       */
  DOUBLE lambda;                                /* nonlinear damp factor in                 */
  /* $step and $nmg_step                      */
  INT max_line_search;                  /* maximum number of line search steps		*/

  /* and DATA_DESC (optional) */
  MATDATA_DESC *J;                              /* the Matrix to be solved					*/
  VECDATA_DESC *d;                              /* nonlinear defect							*/
  VECDATA_DESC *v;                              /* correction computed by newton step           */
  VECDATA_DESC *s;                              /* saved nonlinear solution					*/

} NP_NEWTON;

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

static VEC_SCALAR Factor_One;

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


   ..... fill in data structure here when the realizition is finished


   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPNLSolverInit (NP_NL_SOLVER *np, INT argc , char **argv)
{
  INT i;

  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  if (sc_read(np->abslimit,np->x,"abslimit",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->abslimit[i] = 1.0;
  if (sc_read(np->reduction,NULL,"red",argc,argv))
    return(NP_ACTIVE);
  np->Assemble = (NP_ASSEMBLE *)ReadArgvNumProc(np->base.mg,"A",argc,argv);

  if ((np->x == NULL) || (np->Assemble == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
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
  INT result,level,bl;

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
    if (*np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPNLSolverExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,&bl,&result)) {
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
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPNLSolverExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->x,&result)) {
      UserWriteF("NPNLSolverExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  return(0);
}


/****************************************************************************/
/*D
   newton - numproc for ...

   DESCRIPTION:
   This numproc executes ...

   .vb
   npinit [$x <sol>] [$A <assemble numproc>] [$red <sc double list>]
       [$abslimit <sc double list>] ....
   .ve

   .  $x~<sol> - the solution vector
   .  $abslimit~<sc~double~list> - absolute limit for the defect
   .  $reduction~<sc~double~list> - reduction factor

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate nl $t newton;
   npinit $x sol $A box $red 0.0001 ..............

   npexecute solver $i $d $r $s $p;
   .ve
   D*/
/****************************************************************************/

static INT NewtonConstruct (NP_BASE *theNP)
{
  NP_NL_SOLVER *np;

  /*
     theNP->Init = NewtonInit;
      theNP->Display = NewtonDisplay;
   */
  theNP->Execute = NPNLSolverExecute;
  np = (NP_NL_SOLVER *) theNP;
  /*
      np->PreProcess = NewtonPreProcess;
      np->Solver = NewtonSolver;
      np->PostProcess = NewtonPostProcess;
   */
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
  INT i;

  if (CreateClass ("newton", sizeof(NP_NEWTON), NewtonConstruct))
    return (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  return (0);
}
