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
#include "disctools.h"
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

static DOUBLE *mat,*sol,*def;
static INT *vecskip;

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

INT NPAssembleExecute (NP_BASE *theNP, INT argc , char **argv)
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

  if (ReadArgvOption("i",argc,argv)) {
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

  if (ReadArgvOption("a",argc,argv)) {
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

  if (ReadArgvOption("s",argc,argv)) {
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

  if (ReadArgvOption("d",argc,argv)) {
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

  if (ReadArgvOption("M",argc,argv)) {
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

  if (ReadArgvOption("p",argc,argv)) {
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

/****************************************************************************/
/*D
   NP_LOCAL_ASSEMBLE - type definition for local assembling

   DESCRIPTION:
   This numproc type is used for the description of local assembling.
   It can be called by the given interface from a nonlinear multigrid
   solver.
   Initializing the data is optional; it can with

   'INT NPLocalAssembleInit (NP_ASSEMBLE *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPLocalAssembleDisplay (NP_ASSEMBLE *theNP);'
   'INT NPAssembleExecute (NP_BASE *theNP, INT argc , char **argv);'

   The interface functions 'AssemblePreProcess', 'Assemble'
   'AssembleDefect', 'AssembleMatrix' and 'AssemblePostProcess'
   of NP_ASSEMBLE can be constructed by the interface of NP_LOCAL_ASSEMBLE
   by

   'INT LocalAssembleConstruct (NP_ASSEMBLE *theNP);'

   .vb


   ..... fill in data structure here when the realizition is finished


   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPLocalAssembleInit (NP_LOCAL_ASSEMBLE *np, INT argc , char **argv)
{
  if (ReadArgvINT("g",&np->galerkin,argc,argv))
    np->galerkin = 0;

  return(NPAssembleInit(&np->assemble,argc,argv));
}

INT NPLocalAssembleDisplay (NP_LOCAL_ASSEMBLE *np)
{
  NPAssembleDisplay(&np->assemble);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->galerkin);

  return(0);
}

static INT AssemblePreProcess (NP_ASSEMBLE *theNP, INT level, VECDATA_DESC *x,
                               VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_LOCAL_ASSEMBLE *np;

  np = (NP_LOCAL_ASSEMBLE *) theNP;
  if ((*np->PreProcess)(np,level,x,b,A,&sol,&mat,&def,&vecskip,result)) {
    UserWriteF("PreProcess failed, error code %d\n",result[0]);
    return (1);
  }

  return(0);
}

static INT Assemble (NP_ASSEMBLE *theNP, INT level, VECDATA_DESC *x,
                     VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_LOCAL_ASSEMBLE *np;
  MULTIGRID *theMG;
  GRID *theGrid;
  ELEMENT *theElement;
  DOUBLE *mptr[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
  DOUBLE *sptr[MAX_NODAL_VALUES];
  DOUBLE *rptr[MAX_NODAL_VALUES];
  INT i,j,l,m;

  np = (NP_LOCAL_ASSEMBLE *) theNP;
  theMG = theNP->base.mg;
  for (l=0; l<level; l++) {
    UserWriteF(" [%d:",l);
    theGrid = GRID_ON_LEVEL(theMG,l);
    if (l_dset(theGrid,b,EVERY_CLASS,0.0)!=NUM_OK) {
      result[0] = __LINE__;
      return(1);
    }
    if (l_dmatset(theGrid,A,0.0)!=NUM_OK) {
      result[0] = __LINE__;
      return(1);
    }
    CLEAR_VECSKIP_OF_GRID(theGrid);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
         theElement=SUCCE(theElement)) {
      if (np->galerkin)
        if (NSONS(theElement) > 2) continue;
      m = GetElementVVMPtrs(theElement,x,b,A,sptr,rptr,mptr,vecskip);
      for (i=0; i<m; i++) sol[i] = *sptr[i];
      for (i=0; i<m; i++) def[i] = 0.0;
      for (i=0; i<m*m; i++) mat[i] = 0.0;
      if ((*np->AssembleLocal)(theElement,result)) {
        UserWriteF("AssembleLocal failed for element, error code %d\n",
                   ID(theElement),result[0]);
        return (1);
      }
      for (i=0; i<m; i++) *rptr[i] += def[i];
      for (i=0; i<m*m; i++) *mptr[i] += mat[i];
      for (i=0; i<m; i++) *sptr[i] = sol[i];
      if (OBJT(theElement) == BEOBJ)
        SetElementDirichletFlags(theElement,x,vecskip);
      UserWrite("a]");
    }
  }
  if (theNP->AssembleSolution != NULL)
    if ((*theNP->AssembleSolution)(theNP,level,x,result)) {
      UserWriteF("(AssembleSolution failed, error code %d\n",result[0]);
      return (1);
    }
  if (np->PostMatrix != NULL)
    if ((*np->PostMatrix)(np,level,x,b,A,result)) {
      UserWriteF("(PostMatrix failed, error code %d\n",result[0]);
      return (1);
    }

  return(0);
}

static INT AssemblePostProcess (NP_ASSEMBLE *theNP, INT level, VECDATA_DESC *x,
                                VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_LOCAL_ASSEMBLE *np;

  np = (NP_LOCAL_ASSEMBLE *) theNP;
  if ((*np->PostProcess)(np,level,x,b,A,result)) {
    UserWriteF("PreProcess failed, error code %d\n",result[0]);
    return (1);
  }
  UserWrite("\n");
  return(0);
}

INT LocalAssembleConstruct (NP_ASSEMBLE *np)
{
  np->PreProcess = AssemblePreProcess;
  np->Assemble = Assemble;
  /* np->AssembleSolution  depends on the application */
  np->AssembleDefect = NULL;       /* TODO  */
  np->AssembleMatrix = NULL;       /* TODO */
  np->PostProcess = AssemblePostProcess;

  return(0);
}
