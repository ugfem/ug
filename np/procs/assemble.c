// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  assemble.c		                                                                                */
/*																			*/
/* Purpose:   assemble num procs                                                */
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

#include "assemble.h"

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   NP_ASSEMBLE - type definition for assembling

   DESCRIPTION:
   This numproc type is used for the description of assembling.
   It can be called by the given interface from a nonlinear solver.
   Initializing the data is optional; it can be done with

   'INT NPAssembleInit (NP_ASSEMBLE *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPAssembleDisplay (NP_ASSEMBLE *theNP);'
   'INT NPAssembleExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb


   ..... fill in data structure here when the realizition is finished


   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPAssembleInit (NP_ASSEMBLE *np, INT argc , char **argv)
{
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->c = ReadArgvVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);

  if ((np->A == NULL) || (np->b == NULL) || (np->x == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPAssembleDisplay (NP_ASSEMBLE *np)
{
  if ((np->A == NULL) && (np->b == NULL) && (np->x == NULL))
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

  return(0);
}

INT NPAssmbleExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ASSEMBLE *np;
  INT result,level;

  np = (NP_ASSEMBLE *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"NPAssembleExecute","no vector x");
    return (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPAssembleExecute","no vector b");
    return (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPAssembleExecute","no matrix A");
    return (1);
  }

  if (ReadOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPAssembleExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPAssembleExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("a",argc,argv)) {
    if (np->Assemble == NULL) {
      PrintErrorMessage('E',"NPAssembleExecute","no Assemble");
      return (1);
    }
    if ((*np->Assemble)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPAssembleExecute: Assemble failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("s",argc,argv)) {
    if (np->AssembleSolution == NULL) {
      PrintErrorMessage('E',"NPAssembleExecute","no AssembleSolution");
      return (1);
    }
    if ((*np->AssembleSolution)(np,level,np->x,&result)) {
      UserWriteF("NPAssembleExecute: AssembleSolution failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("d",argc,argv)) {
    if (np->AssembleDefect == NULL) {
      PrintErrorMessage('E',"NPAssembleExecute","no AssembleDefect");
      return (1);
    }
    if ((*np->AssembleDefect)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPAssembleExecute: AssembleDefect failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("M",argc,argv)) {
    if (np->AssembleMatrix == NULL) {
      PrintErrorMessage('E',"NPAssembleExecute","no AssembleMatrix");
      return (1);
    }
    if ((*np->AssembleMatrix)(np,level,np->x,np->b,np->c,np->A,&result)) {
      UserWriteF("NPAssembleExecute: AssembleMatrix failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPAssembleExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPAssembleExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  return(0);
}
