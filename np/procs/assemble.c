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
#include "debug.h"
#include "devices.h"
#include "gm.h"
#include "disctools.h"
#include "np.h"

#include "devices.h"
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

#define MAX_PA                          2

#define PA_VT(pa)                       ((pa)->vt)
#define PA_VD_g(pa)                     ((pa)->gridvel)
#define PA_NASS(pa)                     ((pa)->nass)
#define PA_SUB(pa,i)            ((pa)->sub[i])
#define PA_ASS(pa,i)            ((pa)->ass[i])
#define PA_DT(pa)                       ((pa)->dt)

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  NP_NL_ASSEMBLE pa;                                    /* class for nonlinear ass routines		*/

  /* additional data */
  VEC_TEMPLATE *vt;                                     /* vector template for part decomp		*/
  VECDATA_DESC *gridvel;                        /* velocity of moving grid (iff)		*/

  INT nass;                                                     /* number of part assembling numprocs	*/
  INT sub[MAX_PA];                                      /* sub vector of vt for part decomp		*/
  NP_NL_PARTASS *ass[MAX_PA];                   /* pointers to part assembling numprocs	*/

} NP_PA_NL;

typedef struct
{
  NP_T_ASSEMBLE pa;                                     /* class for nonlinear ass routines		*/

  /* additional data */
  VEC_TEMPLATE *vt;                                     /* vector template for part decomp		*/
  VECDATA_DESC *gridvel;                        /* velocity of moving grid (iff)		*/

  INT nass;                                                     /* number of part assembling numprocs	*/
  INT sub[MAX_PA];                                      /* sub vector of vt for part decomp		*/
  NP_T_PARTASS *ass[MAX_PA];                    /* pointers to part assembling numprocs	*/

  DOUBLE dt;                                                    /* time step							*/

} NP_PA_T;

/* TODO (HRR 971118): remove old part assemble structs (2 items) */
typedef struct
{
  NP_NL_ASSEMBLE pa;                                    /* class for nonlinear ass routines		*/

  /* additional data */
  INT nass;                                                     /* number of part assembling numprocs	*/
  NP_NL_ASSEMBLE *ass[MAX_PA];          /* pointers to part assembling numprocs	*/

} OLD_NP_NL_PARTASS;

typedef struct
{
  NP_T_ASSEMBLE pa;                                     /* class for nonlinear ass routines		*/

  /* additional data */
  INT nass;                                                     /* number of part assembling numprocs	*/
  NP_T_ASSEMBLE *ass[MAX_PA];                   /* pointers to part assembling numprocs	*/

} OLD_NP_T_PARTASS;

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
   assemble - ug offers several classes of assmbling 'num_proc's

   DESCRIPTION:
   The classes defined in ug are:~
   .  NP_ASSEMBLE - type definition for assembling
   .  NP_NL_ASSEMBLE - type definition for nonlinear assembling
   .  NP_LOCAL_ASSEMBLE - type definition for local assembling
   .  NP_T_ASSEMBLE - type definition for time dependent assembling

   For realizations of assembling num_procs (which are defined in the problem
   classes) try typing
   .n     help assemble $k

   SEE ALSO:
   'num_proc', 'NP_ASSEMBLE', 'NP_NL_ASSEMBLE', 'NP_LOCAL_ASSEMBLE', 'NP_T_ASSEMBLE'
   D*/
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
   The data can be displayed and the num proc can be executed by

   'INT NPAssembleDisplay (NP_ASSEMBLE *theNP);'
   'INT NPAssembleExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   struct np_assemble {

        NP_BASE base;                        // inherits base class

        // data (optinal, necessary for calling the generic execute routine)
    VECDATA_DESC *x;                     // solution
    VECDATA_DESC *b;                     // defect
    MATDATA_DESC *A;                     // matrix

        // functions
        INT (*PreProcess)
             (struct np_assemble *,          // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // rhs vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
    INT (*Assemble)
             (struct np_assemble *,          // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,			     // current solution	(initial)
                  VECDATA_DESC *,			     // right hand side
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*PostProcess)
             (struct np_assemble *,          // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
   };
   typedef struct np_assemble NP_ASSEMBLE;
   .ve

   SEE ALSO:
   'num_proc', 'NP_NL_ASSEMBLE', 'NP_LOCAL_ASSEMBLE', 'NP_T_ASSEMBLE'
   D*/
/****************************************************************************/

INT NPAssembleInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ASSEMBLE *np;

  np = (NP_ASSEMBLE *) theNP;
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);

  if ((np->A == NULL) || (np->b == NULL) || (np->x == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPAssembleDisplay (NP_BASE *theNP)
{
  NP_ASSEMBLE *np;

  np = (NP_ASSEMBLE *) theNP;
  if ((np->A == NULL) && (np->b == NULL) && (np->x == NULL))
    return(0);
  UserWrite("symbolic user data:\n");
  if (np->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  if (np->b != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->b));
  if (np->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
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
   NP_NL_ASSEMBLE - type definition for nonlinear assembling

   DESCRIPTION:
   This numproc type is used for the description of assembling.
   It can be called by the given interface from a nonlinear solver.
   Initializing the data is optional; it can be done with

   'INT NPNLAssembleInit (NP_BASE *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data can be displayed and the num proc can be executed by

   'INT NPNLAssembleDisplay (NP_BASE *theNP);'
   'INT NPNLAssembleExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   struct np_nl_assemble {

        NP_BASE base;                        // inherits base class

        // data (optinal, necessary for calling the generic execute routine)
    VECDATA_DESC *x;                     // solution
    VECDATA_DESC *c;                     // correction
    VECDATA_DESC *b;                     // defect
    MATDATA_DESC *A;                     // matrix

        // functions
        INT (*PreProcess)
             (struct np_nl_assemble *,       // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
    INT (*NLAssembleSolution)
             (struct np_nl_assemble *,       // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  INT *);                        // result
    INT (*NLAssembleDefect)
             (struct np_nl_assemble *,       // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
    INT (*NLAssembleMatrix)
             (struct np_nl_assemble *,       // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,			     // current solution	(initial)
                  VECDATA_DESC *,			     // defect for current solution
                  VECDATA_DESC *,			     // correction to be computed
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*PostProcess)
             (struct np_nl_assemble *,       // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
   };
   typedef struct np_nl_assemble NP_NL_ASSEMBLE;
   .ve

   SEE ALSO:
   'num_proc', 'NP_ASSEMBLE', 'NP_LOCAL_ASSEMBLE', 'NP_T_ASSEMBLE'
   D*/
/****************************************************************************/

INT NPNLAssembleInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_NL_ASSEMBLE *np;

  np = (NP_NL_ASSEMBLE *) theNP;
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->c = ReadArgvVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);

  if ((np->A == NULL) || (np->b == NULL) || (np->x == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPNLAssembleDisplay (NP_BASE *theNP)
{
  NP_NL_ASSEMBLE *np;

  np = (NP_NL_ASSEMBLE *) theNP;
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

INT NPNLAssembleExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_NL_ASSEMBLE *np;
  INT result,level;

  np = (NP_NL_ASSEMBLE *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"NPNLAssembleExecute","no vector x");
    return (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPNLAssembleExecute","no vector b");
    return (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPNLAssembleExecute","no matrix A");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,0,level,np->x,&result)) {
      UserWriteF("NPNLAssembleExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->NLAssembleSolution == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no NLAssembleSolution");
      return (1);
    }
    if ((*np->NLAssembleSolution)(np,0,level,np->x,&result)) {
      UserWriteF("NPNLAssembleExecute: NLAssembleSolution failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("d",argc,argv)) {
    if (np->NLAssembleDefect == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no NLAssembleDefect");
      return (1);
    }
    if ((*np->NLAssembleDefect)(np,0,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPNLAssembleExecute: NLAssembleDefect failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("M",argc,argv)) {
    if (np->NLAssembleMatrix == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no NLAssembleMatrix");
      return (1);
    }
    if ((*np->NLAssembleMatrix)(np,0,level,np->x,np->b,np->c,np->A,&result)) {
      UserWriteF("NPNLAssembleExecute: NLAssembleMatrix failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,0,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPNLAssembleExecute: PostProcess failed, error code %d\n",
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

   'INT NPLocalAssembleInit (NP_LOCAL_ASSEMBLE *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data can be displayed and the num proc can be executed by

   'INT NPLocalAssembleDisplay (NP_LOCAL_ASSEMBLE *theNP);'
   'INT NPAssembleExecute (NP_BASE *theNP, INT argc , char **argv);'

   The interface functions 'LocalAssemblePreProcess', 'LocalAssemble'
   'AssembleDefect', 'AssembleMatrix' and 'LocalAssemblePostProcess'
   of NP_ASSEMBLE can be constructed by the interface of NP_LOCAL_ASSEMBLE
   by

   'INT LocalAssembleConstruct (NP_ASSEMBLE *theNP);'

   .vb
   struct np_local_assemble {

        NP_ASSEMBLE assemble;                // inherits assemble class

        // data
        INT galerkin;                        // Galerkin assembling

        // functions
        INT (*PreProcess)
             (struct np_local_assemble *,    // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
          DOUBLE **,                     // local solution
          DOUBLE **,                     // local defect
          DOUBLE **,                     // local matrix
          INT **,                        // local vecskip
                  INT *);                        // result
    INT (*AssembleLocal)
             (ELEMENT *,                     // pointer to an element
                  INT *);                        // result
    INT (*AssembleLocalDefect)
             (ELEMENT *,                     // pointer to an element
                  INT *);                        // result
    INT (*AssembleLocalMatrix)
             (ELEMENT *,                     // pointer to an element
                  INT *);                        // result
        INT (*PostMatrix)
             (struct np_local_assemble *,    // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*PostProcess)
             (struct np_local_assemble *,    // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
   };
   typedef struct np_local_assemble NP_LOCAL_ASSEMBLE;
   .ve

   SEE ALSO:
   'num_proc', 'NP_ASSEMBLE', 'NP_NL_ASSEMBLE', 'NP_T_ASSEMBLE'
   D*/
/****************************************************************************/

INT NPLocalAssembleInit (NP_LOCAL_ASSEMBLE *np, INT argc , char **argv)
{
  if (ReadArgvINT("g",&np->galerkin,argc,argv))
    np->galerkin = 0;

  return(NPAssembleInit(&np->assemble.base,argc,argv));
}

INT NPLocalAssembleDisplay (NP_LOCAL_ASSEMBLE *np)
{
  NPAssembleDisplay(&np->assemble.base);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->galerkin);

  return(0);
}

INT NPLocalAssemblePostMatrix (NP_LOCAL_ASSEMBLE *theNP, INT level,
                               VECDATA_DESC *x,
                               VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  INT lev;

        #ifdef ModelP
  if (a_vector_vecskip(theNP->assemble.base.mg,0,level,x) != NUM_OK)
    return (1);
        #endif

  for (lev=0; lev<=level; lev++)
    AssembleDirichletBoundary(GRID_ON_LEVEL(theNP->assemble.base.mg,lev),
                              A,x,b);
  UserWrite(" [d]");

  return(0);
}

static INT LocalAssemblePreProcess (NP_ASSEMBLE *theNP, INT level, VECDATA_DESC *x,
                                    VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_LOCAL_ASSEMBLE *np;

  np = (NP_LOCAL_ASSEMBLE *) theNP;
  if ((*np->PreProcess)(np,level,x,b,A,&sol,&def,&mat,&vecskip,result)) {
    UserWriteF("PreProcess failed, error code %d\n",result[0]);
    return (1);
  }

  return(0);
}

static INT LocalAssemble (NP_ASSEMBLE *theNP, INT level, VECDATA_DESC *x,
                          VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_LOCAL_ASSEMBLE *np;
  MULTIGRID *theMG;
  GRID *theGrid;
  ELEMENT *theElement;
  DOUBLE *mptr[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
  DOUBLE *sptr[MAX_NODAL_VALUES];
  DOUBLE *rptr[MAX_NODAL_VALUES];
  INT i,l,m;

  np = (NP_LOCAL_ASSEMBLE *) theNP;
  theMG = NP_MG(theNP);
  for (l=0; l<=level; l++) {
    UserWriteF(" [%d:",l);
    theGrid = GRID_ON_LEVEL(theMG,l);
    if (dset(theMG,l,l,ALL_VECTORS,b,0.0)!=NUM_OK) NP_RETURN(1,result[0]);
    if (dmatset(theMG,l,l,ALL_VECTORS,A,0.0)!=NUM_OK) NP_RETURN(1,result[0]);
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
        UserWriteF("AssembleLocal failed for element %d, error code %d\n",
                   ID(theElement),result[0]);
        return (1);
      }
      for (i=0; i<m; i++) *rptr[i] += def[i];
      for (i=0; i<m*m; i++) *mptr[i] += mat[i];
      for (i=0; i<m; i++) *sptr[i] = sol[i];
      if (OBJT(theElement) == BEOBJ)
        SetElementDirichletFlags(theElement,x,vecskip);
    }
    UserWrite("a]");
  }
  if (np->PostMatrix != NULL)
    if ((*np->PostMatrix)(np,level,x,b,A,result)) {
      UserWriteF("(PostMatrix failed, error code %d\n",result[0]);
      return (1);
    }
  UserWrite("\n");

  return(0);
}

static INT LocalAssemblePostProcess (NP_ASSEMBLE *theNP, INT level, VECDATA_DESC *x,
                                     VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_LOCAL_ASSEMBLE *np;

  np = (NP_LOCAL_ASSEMBLE *) theNP;
  if (np->PostProcess != NULL)
    if ((*np->PostProcess)(np,level,x,b,A,result)) {
      UserWriteF("PostProcess failed, error code %d\n",result[0]);
      return (1);
    }
  UserWrite("\n");

  return(0);
}

INT NPLocalAssembleConstruct (NP_ASSEMBLE *np)
{
  np->PreProcess          = LocalAssemblePreProcess;
  np->Assemble            = LocalAssemble;
  np->PostProcess         = LocalAssemblePostProcess;

  return(0);
}

/****************************************************************************/
/*D
        partass - nonlinear assemble numproc calling several assembling numprocs

        DESCRIPTION:
   .vb
        npcreate pa $c partass
        npinit	$A <mat data desc>
                        $s <vec data desc>
                        $r <vec data desc>

                        {$ass <nl part ass>}+
   .ve
        Data descriptors:~
   .   A~<mat~data~desc>		- stiffnes matrix to assemble
   .   s~<vec~data~desc>		- current solution
   .   r~<vec~data~desc>                - source term

        Parameters:~
   .   ass~<nl~part~ass>		- specification of a nonlinear assembling numproc called by 'pa'

        DESCRIPTION:
        The numerical procedure class 'partass' can call several nonlinear
        assembling numprocs. 'partass' does not check any compatibility or
        consistency!

        This numproc is useful for coupled problems using different discretizations
        in different parts of the domain.

        KEYWORDS:
        assemble, part
   D*/
/****************************************************************************/

static INT OldNLPartAssInit (NP_BASE *theNP, INT argc, char **argv)
{
  OLD_NP_NL_PARTASS *thePA;
  NP_NL_ASSEMBLE *ass;
  INT r,i,nass;
  char name[NAMESIZE];

  r = NPNLAssembleInit(theNP,argc,argv);

  thePA = (OLD_NP_NL_PARTASS*)theNP;

  PA_NASS(thePA) = nass = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (nass>=MAX_PA)
      {
        PrintErrorMessage('E',"OldNLPartAssInit","max number of part assembling numprocs exceeded");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      if (sscanf(argv[i],expandfmt(CONCAT3("ass %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"OldNLPartAssInit","specify a nonlinear assembling numproc with $ass");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      ass = (NP_NL_ASSEMBLE*) GetNumProcByName (NP_MG(theNP),name,NL_ASSEMBLE_CLASS_NAME);
      if (ass == NULL)
      {
        PrintErrorMessage('E',"OldNLPartAssInit",
                          "cannot find specified numerical procedure");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PA_ASS(thePA,nass++) = ass;
      break;
    }
  if (nass==0)
  {
    PrintErrorMessage('E',"OldNLPartAssInit","specify at least one nonlinear assembling numproc with $ass");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PA_NASS(thePA) = nass;

  if ((r==NP_ACTIVE) || (r==NP_EXECUTABLE))
    return (r);
  else
    REP_ERR_RETURN (r);
}

/****************************************************************************/
/*
        the follwing functions do nothing but just calling the corresponding
        routines of the nonlinear assembling numprocs specified when initializing
        the numproc
 */
/****************************************************************************/

static INT OldNLPartAssDisplay (NP_BASE *theNP)
{
  OLD_NP_NL_PARTASS *thePA;
  INT i;
  char text[8];

  NPNLAssembleDisplay(theNP);

  thePA   = (OLD_NP_NL_PARTASS*)theNP;

  UserWrite("\npart assembling numprocs:\n");
  for (i=0; i<PA_NASS(thePA); i++)
  {
    sprintf(text,"ass%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SS,text,ENVITEM_NAME(PA_ASS(thePA,i)));
  }

  return (0);
}

static INT OldNLPartAssPreProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, INT *res)
{
  OLD_NP_NL_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_NL_PARTASS*)ass;

  /* call prep routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPANL_PRE(PA_ASS(thePA,i))!=NULL)
      if (NPANL_PRE(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,x,res))
        REP_ERR_RETURN(1);

  return(0);
}

static INT OldNLPartAssPostProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x,
                                    VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  OLD_NP_NL_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_NL_PARTASS*)ass;

  /* call post routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPANL_POST(PA_ASS(thePA,i))!=NULL)
      if (NPANL_POST(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,x,d,J,res))
        REP_ERR_RETURN(1);

  return(0);
}

static INT OldNLPartAssAssembleSolution (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, INT *res)
{
  OLD_NP_NL_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_NL_PARTASS*)ass;

  /* call assemble solution routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPANL_ASSSOL(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,u,res))
      REP_ERR_RETURN(1);

  return(0);
}

static INT OldNLPartAssAssembleDefect (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *s,
                                       VECDATA_DESC *r, MATDATA_DESC *A, INT *res)
{
  OLD_NP_NL_PARTASS *thePA;
  INT i,l;

  thePA   = (OLD_NP_NL_PARTASS*)ass;

  /* first clear skip flags of global problem */
  for (l=fl; l<=tl; l++)
    ClearVecskipFlags(GRID_ON_LEVEL(NP_MG(ass),l),s);

  /* call assemble defect routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPANL_ASSDEF(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,s,r,A,res))
      REP_ERR_RETURN(1);

  return(0);
}

static INT OldNLPartAssAssembleMatrix (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *s,
                                       VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *A, INT *res)
{
  OLD_NP_NL_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_NL_PARTASS*)ass;

  /* clear global matrix */
  if (dmatset(NP_MG(ass),fl,tl,ALL_VECTORS,A,0.0)!=NUM_OK) REP_ERR_RETURN (__LINE__);

  /* call assemble matrix routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPANL_ASSMAT(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,s,d,v,A,res))
      REP_ERR_RETURN(1);

  return(0);
}

static INT NLPartAssembleConstruct (NP_BASE *theNP)
{
  NP_NL_ASSEMBLE *np;

  np = (NP_NL_ASSEMBLE *) theNP;

  np->base.Init                       = OldNLPartAssInit;
  np->base.Display            = OldNLPartAssDisplay;
  np->base.Execute            = NPNLAssembleExecute;
  np->PreProcess                      = OldNLPartAssPreProcess;
  np->PostProcess             = OldNLPartAssPostProcess;
  np->NLAssembleSolution      = OldNLPartAssAssembleSolution;
  np->NLAssembleDefect        = OldNLPartAssAssembleDefect;
  np->NLAssembleMatrix        = OldNLPartAssAssembleMatrix;

  return(0);
}

/****************************************************************************/
/*D
        parttass - nonlinear time dependent assemble numproc calling several assembling numprocs

        DESCRIPTION:
   .vb
        npcreate pa $c parttass
        npinit	$A <mat data desc>
                        $s <vec data desc>
                        $r <vec data desc>

                        {$ass <tnl part ass>}+
   .ve
        Data descriptors:~
   .   A~<mat~data~desc>		- stiffnes matrix to assemble
   .   s~<vec~data~desc>		- current solution
   .   r~<vec~data~desc>                - source term

        Parameters:~
   .   ass~<tnl~part~ass>		- specification of a nonlinear time dependent assembling numproc
                                                                called by 'pa'

        DESCRIPTION:
        The numerical procedure class 'parttass' can call several nonlinear
        assembling numprocs. 'parttass' does not check any compatibility or
        consistency!

        This numproc is useful for coupled problems using different discretizations
        in different parts of the domain.

        KEYWORDS:
        assemble, part
   D*/
/****************************************************************************/

static INT OldTPartAssInit (NP_BASE *theNP, INT argc, char **argv)
{
  OLD_NP_T_PARTASS *thePA;
  NP_T_ASSEMBLE *ass;
  INT r,i,nass;
  char name[NAMESIZE];

  r = NPTAssembleInit(theNP,argc,argv);

  thePA = (OLD_NP_T_PARTASS*)theNP;

  PA_NASS(thePA) = nass = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (nass>=MAX_PA)
      {
        PrintErrorMessage('E',"TPartAssInit","max number of part assembling numprocs exceeded");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      if (sscanf(argv[i],expandfmt(CONCAT3("ass %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"TPartAssInit","specify a nonlinear assembling numproc with $ass");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      ass = (NP_T_ASSEMBLE*) GetNumProcByName (NP_MG(theNP),name,T_ASSEMBLE_CLASS_NAME);
      if (ass == NULL)
      {
        PrintErrorMessage('E',"TPartAssInit",
                          "cannot find specified numerical procedure");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PA_ASS(thePA,nass++) = ass;
      break;
    }
  if (nass==0)
  {
    PrintErrorMessage('E',"TPartAssInit","specify at least one nonlinear assembling numproc with $ass");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PA_NASS(thePA) = nass;

  if ((r==NP_ACTIVE) || (r==NP_EXECUTABLE))
    return (r);
  else
    REP_ERR_RETURN (r);
}

/****************************************************************************/
/*
        the follwing functions do nothing but just calling the corresponding
        routines of the nonlinear assembling numprocs specified when initializing
        the numproc
 */
/****************************************************************************/

static INT OldTPartAssDisplay (NP_BASE *theNP)
{
  OLD_NP_T_PARTASS *thePA;
  INT i;
  char text[8];

  NPTAssembleDisplay(theNP);

  thePA   = (OLD_NP_T_PARTASS*)theNP;

  UserWrite("\npart assembling numprocs:\n");
  for (i=0; i<PA_NASS(thePA); i++)
  {
    sprintf(text,"ass%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SS,text,ENVITEM_NAME(PA_ASS(thePA,i)));
  }

  return (0);
}

static INT OldTPartAssPreProcess (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE t_p1, DOUBLE t_0,
                                  DOUBLE t_m1, VECDATA_DESC *u_p1, VECDATA_DESC *u_0,
                                  VECDATA_DESC *u_m1, INT *res)
{
  OLD_NP_T_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_T_PARTASS*)ass;

  /* call prep routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPAT_PRE(PA_ASS(thePA,i))!=NULL)
      if (NPAT_PRE(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,t_p1,t_0,t_m1,u_p1,u_0,u_m1,res))
        REP_ERR_RETURN(1);

  return(0);
}

static INT OldTPartAssPostProcess (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE t_p1, DOUBLE t_0,
                                   DOUBLE t_m1, VECDATA_DESC *u_p1, VECDATA_DESC *u_0,
                                   VECDATA_DESC *u_m1, INT *res)
{
  OLD_NP_T_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_T_PARTASS*)ass;

  /* call post routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPAT_POST(PA_ASS(thePA,i))!=NULL)
      if (NPAT_POST(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,t_p1,t_0,t_m1,u_p1,u_0,u_m1,res))
        REP_ERR_RETURN(1);

  return(0);
}

static INT OldTPartAssAssembleInitial (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE time, VECDATA_DESC *u, INT *res)
{
  OLD_NP_T_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_T_PARTASS*)ass;

  /* call assemble solution routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPAT_INITIAL(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,time,u,res))
      REP_ERR_RETURN(1);

  return(0);
}

static INT OldTPartAssAssembleSolution (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE time, VECDATA_DESC *u, INT *res)
{
  OLD_NP_T_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_T_PARTASS*)ass;

  /* call assemble solution routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPAT_ASSSOL(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,time,u,res))
      REP_ERR_RETURN(1);

  return(0);
}

static INT OldTPartAssAssembleDefect (NP_T_ASSEMBLE *ass, INT fl, INT tl,
                                      DOUBLE time, DOUBLE s_m, DOUBLE s_a, VECDATA_DESC *s,
                                      VECDATA_DESC *r, MATDATA_DESC *A, INT *res)
{
  OLD_NP_T_PARTASS *thePA;
  INT i,l;

  thePA   = (OLD_NP_T_PARTASS*)ass;

  /* first clear skip flags of global problem */
  for (l=fl; l<=tl; l++)
    ClearVecskipFlags(GRID_ON_LEVEL(NP_MG(ass),l),s);

  /* call assemble defect routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPAT_ASSDEF(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,time,s_m,s_a,s,r,A,res))
      REP_ERR_RETURN(1);

  return(0);
}

static INT OldTPartAssAssembleMatrix (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE time, DOUBLE s_a,
                                      VECDATA_DESC *s, VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *A, INT *res)
{
  OLD_NP_T_PARTASS *thePA;
  INT i;

  thePA   = (OLD_NP_T_PARTASS*)ass;

  /* clear global matrix */
  if (dmatset(NP_MG(ass),fl,tl,ALL_VECTORS,A,0.0)!=NUM_OK) REP_ERR_RETURN (__LINE__);

  /* call assemble matrix routines */
  for (i=0; i<PA_NASS(thePA); i++)
    if (NPAT_ASSMAT(PA_ASS(thePA,i)) (PA_ASS(thePA,i),fl,tl,time,s_a,s,d,v,A,res))
      REP_ERR_RETURN(1);

  return(0);
}

static INT TPartAssembleConstruct (NP_BASE *theNP)
{
  NP_T_ASSEMBLE *np;

  np = (NP_T_ASSEMBLE *) theNP;

  np->base.Init                               = OldTPartAssInit;
  np->base.Display                    = OldTPartAssDisplay;
  np->base.Execute                    = NPTAssembleExecute;
  np->TAssemblePreProcess     = OldTPartAssPreProcess;
  np->TAssemblePostProcess    = OldTPartAssPostProcess;
  np->TAssembleInitial                = OldTPartAssAssembleInitial;
  np->TAssembleSolution               = OldTPartAssAssembleSolution;
  np->TAssembleDefect                 = OldTPartAssAssembleDefect;
  np->TAssembleMatrix                 = OldTPartAssAssembleMatrix;

  return(0);
}

/****************************************************************************/
/*D
   NP_T_ASSEMBLE - type definition for time dependent assembling

   DESCRIPTION:
   This is the interface for a time dependent problem as it is required
   by the tsolver. An NP_T_ASSEMBLE object is never executable, only its
   functional interface is used.

   .vb
   struct np_t_assemble {

        NP_BASE base;                        // inherits base class

        // functions
        INT (*TAssemblePreProcess)               // call at begin of timestep
             (struct np_t_assemble *,        // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  DOUBLE,						 // time t_k+1
                  DOUBLE,						 // time t_k
                  DOUBLE,						 // time t_k-1
                  VECDATA_DESC *,                // (unknown) solution at t_k+1
                  VECDATA_DESC *,                // solution vector at t_k
                  VECDATA_DESC *,                // solution vector at t_k-1
                  INT *);                        // result
    INT (*TAssembleInitial)                      // set initial values
             (struct np_t_assemble *,        // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  DOUBLE,						 // time value t
                  VECDATA_DESC *,                // solution vector at time t
                  INT *);                        // result
    INT (*TAssembleSolution)             // set dirichlet conditions in sol.
             (struct np_t_assemble *,        // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  DOUBLE,						 // time value t
                  VECDATA_DESC *,                // solution vector at time t
                  INT *);                        // result
    INT (*TAssembleDefect)                   // accumulate to defect vector
             (struct np_t_assemble *,        // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  DOUBLE,						 // time value t
                  DOUBLE,						 // scaling for m-term: s_m
                  DOUBLE,						 // scaling for a-term: s_a
                  VECDATA_DESC *,                // solution vector y
                  VECDATA_DESC *,                // accumulate s_m*m(t,y)+s_a*a(t,y)
                  MATDATA_DESC *,                // matrix may be handy for Picard
                  INT *);                        // result
    INT (*TAssembleMatrix)                       // compute linearization (Jacobian)
             (struct np_t_assemble *,        // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  DOUBLE,						 // time value t
                  DOUBLE,						 // scaling for a-term: s_a	(s_m=1!)
                  VECDATA_DESC *,			     // current sol (linearization pt)
                  VECDATA_DESC *,			     // defect for current solution
                  VECDATA_DESC *,			     // correction to be computed
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*TAssemblePostProcess)          // call after solution t_k+1 known
             (struct np_t_assemble *,        // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  DOUBLE,						 // time t_k+1
                  DOUBLE,						 // time t_k
                  DOUBLE,						 // time t_k-1
                  VECDATA_DESC *,                // solution t_k+1 (just computed!)
                  VECDATA_DESC *,                // solution vector at t_k
                  VECDATA_DESC *,                // solution vector at t_k-1
                  INT *);                        // result
   };
   typedef struct np_t_assemble NP_T_ASSEMBLE;
   .ve

   SEE ALSO:
   'num_proc', 'NP_ASSEMBLE', 'NP_NL_ASSEMBLE', 'NP_LOCAL_ASSEMBLE'
   D*/
/****************************************************************************/

INT NPTAssembleInit (NP_BASE *theNP, INT argc , char **argv)
{
  return(NP_ACTIVE);
}

INT NPTAssembleDisplay (NP_BASE *theNP)
{
  return(0);
}

INT NPTAssembleExecute (NP_BASE *theNP, INT argc , char **argv)
{
  REP_ERR_RETURN(1);   /* never executable */
}

/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

void DefaultPartassParams (PARTASS_PARAMS *pp)
{
  memset(pp,0.0,sizeof(PARTASS_PARAMS));

  PP_ACTION(pp)           = PARTASS_UNKNOWN;
  PP_SCALE_A(pp)          = 1.0;
  PP_SCALE_M(pp)          = 0.0;
  PP_TIME(pp)             = 0.0;
  PP_DELTA_T(pp)          = 0.0;
  PP_ASS_PART(pp)         = NO;
  PP_SKIP(pp)             = NULL;
  PP_CO_SKIP(pp)          = NULL;
  PP_MD_A(pp)             = NULL;
  PP_MD_A_glob(pp)        = NULL;
  PP_VD_s(pp)             = NULL;
  PP_VD_s_glob(pp)        = NULL;
  PP_VD_s_i(pp)           = NULL;
  PP_VD_s_co(pp)          = NULL;
  PP_VD_o(pp)             = NULL;
  PP_VD_o_glob(pp)        = NULL;
  PP_VD_c(pp)             = NULL;
  PP_VD_c_glob(pp)        = NULL;
  PP_VD_r(pp)             = NULL;
  PP_VD_r_glob(pp)        = NULL;
  PP_VD_gridvel(pp)       = NULL;

  return;
}

/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

INT SetPartassParams (PARTASS_PARAMS *pp, const VEC_TEMPLATE *vt, INT sub,
                      DOUBLE s_a, DOUBLE s_m, DOUBLE t, DOUBLE dt,
                      VECDATA_DESC *s, VECDATA_DESC *r, VECDATA_DESC *o,
                      VECDATA_DESC *c, VECDATA_DESC *g, MATDATA_DESC *A)
{
  /* checks */
  if (s==NULL)
    REP_ERR_RETURN (1);
  if (vt==NULL)
    REP_ERR_RETURN (1);
  if (sub<0 || sub>=VT_NSUB(vt))
    REP_ERR_RETURN (1);

  /* clear */
  memset(pp,0.0,sizeof(PARTASS_PARAMS));

  /* general */
  PP_ASS_PART(pp)         = TRUE;
  PP_ACTION(pp)           = PARTASS_UNKNOWN;
  PP_SCALE_A(pp)          = s_a;
  PP_SCALE_M(pp)          = s_m;
  PP_TIME(pp)                     = t;
  PP_DELTA_T(pp)          = dt;
  PP_MD_A_glob(pp)        = A;
  PP_VD_s_glob(pp)        = s;
  PP_VD_o_glob(pp)        = o;
  PP_VD_c_glob(pp)        = c;
  PP_VD_r_glob(pp)        = r;
  PP_VD_gridvel(pp)       = g;

  /* decompose descriptors */
  if (!VDmatchesVT(s,vt))
    REP_ERR_RETURN(1);
  if (VDsubDescFromVT(s,vt,sub,&PP_VD_s(pp)))
    REP_ERR_RETURN(1)
    if (VDinterfaceDesc(s,PP_VD_s(pp),&PP_VD_s_i(pp)))
      REP_ERR_RETURN(1)
      if (VDinterfaceCoDesc(s,PP_VD_s(pp),&PP_VD_s_co(pp)))
        REP_ERR_RETURN(1)
        if (ComputePartVecskip(s,PP_VD_s(pp),PP_SKIP(pp),PP_CO_SKIP(pp)))
          REP_ERR_RETURN(1)

          if (o!=NULL)
          {
            if (!VDmatchesVT(o,vt))
              REP_ERR_RETURN(1);
            if (VDsubDescFromVT(o,vt,sub,&PP_VD_o(pp)))
              REP_ERR_RETURN(1);
          }
  if (c!=NULL)
  {
    if (!VDmatchesVT(c,vt))
      REP_ERR_RETURN(1);
    if (VDsubDescFromVT(c,vt,sub,&PP_VD_c(pp)))
      REP_ERR_RETURN(1);
  }
  if (r!=NULL)
  {
    if (!VDmatchesVT(r,vt))
      REP_ERR_RETURN(1);
    if (VDsubDescFromVT(r,vt,sub,&PP_VD_r(pp)))
      REP_ERR_RETURN(1);
  }
  if (A!=NULL)
  {
    if (!MDmatchesVT(A,vt))
      REP_ERR_RETURN(1);
    if (MDsubDescFromVT(A,vt,sub,&PP_MD_A(pp)))
      REP_ERR_RETURN(1);
  }
  return (0);
}

/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

static INT NLPartAssInit (NP_BASE *theNP, INT argc, char **argv)
{
  NP_PA_NL *pa = (NP_PA_NL*) theNP;
  VEC_TEMPLATE *mvt;
  INT r,i,j,nass;
  char name[NAMESIZE],buffer[VALUELEN];

  r = NPNLAssembleInit(theNP,argc,argv);

  /* get name of main vector template */
  if (ReadArgvChar("m",buffer,argc,argv)!=NULL)
  {
    PrintErrorMessage('E',"NLPartAssInit","m option with main vector template not found");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  mvt = GetVectorTemplate(MGFORMAT(NP_MG(theNP)),buffer);
  if (mvt == NULL)
  {
    PrintErrorMessageF('E',"NLPartAssInit",
                       "cannot find specified vector template '%s'",buffer);
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PA_VT(pa) = mvt;

  PA_VD_g(pa) = ReadArgvVecDesc(NP_MG(theNP),"g",argc,argv);

  PA_NASS(pa) = nass = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (nass>=MAX_PA)
      {
        PrintErrorMessage('E',"NLPartAssInit","max number of part assembling numprocs exceeded");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      if (sscanf(argv[i],expandfmt(CONCAT3("ass %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"NLPartAssInit","specify a nonlinear part assembling numproc with $ass");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PA_ASS(pa,nass) = (NP_NL_PARTASS*) GetNumProcByName (NP_MG(theNP),name,NL_PARTASS_CLASS_NAME);
      if (PA_ASS(pa,nass) == NULL)
      {
        PrintErrorMessage('E',"NLPartAssInit",
                          "cannot find specified numerical procedure");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }

      /* next arg has to specify a valid sub vector */
      if (++i>=argc)
      {
        PrintErrorMessage('E',"NLPartAssInit",
                          "last ass option has no sub option");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      /* scan name of subvector of main vector template */
      if (sscanf(argv[i],expandfmt(CONCAT3("sub %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"NLPartAssInit","s option expected after ass option");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      for (j=0; j<VT_NSUB(mvt); j++)
        if (strcmp(SUBV_NAME(VT_SUB(mvt,j)),name)==0)
          break;
      if (j>=VT_NSUB(mvt))
      {
        PrintErrorMessageF('E',"NLPartAssInit","name '%s' of sub template not found",name);
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PA_SUB(pa,nass) = j;
      nass++;
      break;
    }
  if (nass==0)
  {
    PrintErrorMessage('E',"NLPartAssInit","specify at least one nonlinear assembling numproc with $ass");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PA_NASS(pa) = nass;

  if ((r==NP_ACTIVE) || (r==NP_EXECUTABLE))
    return (r);
  else
    REP_ERR_RETURN (r);
}

static INT NLPartAssDisplay (NP_BASE *theNP)
{
  NP_PA_NL *pa = (NP_PA_NL*) theNP;
  INT i;
  char text[8];

  NPNLAssembleDisplay(theNP);

  if (PA_VD_g(pa)!=NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"g",ENVITEM_NAME(PA_VD_g(pa)));

  UserWriteF(DISPLAY_NP_FORMAT_SS,"vec tmplt",ENVITEM_NAME(PA_VT(pa)));

  UserWrite("\npart assembling numprocs:\n");
  for (i=0; i<PA_NASS(pa); i++)
  {
    sprintf(text,"ass%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SSS,text,ENVITEM_NAME(PA_ASS(pa,i)),SUBV_NAME(VT_SUB(PA_VT(pa),PA_SUB(pa,i))));
  }

  return (0);
}

static INT NLPartAssPreProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, INT *res)
{
  NP_PA_NL *pa = (NP_PA_NL*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* call prep routines */
  for (i=0; i<PA_NASS(pa); i++)
    if (NPPNL_PRE(PA_ASS(pa,i))!=NULL)
    {
      if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,0.,0.,x,NULL,NULL,NULL,PA_VD_g(pa),NULL))
        REP_ERR_RETURN (1);
      if (NPPNL_PRE(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
        REP_ERR_RETURN(1);
    }

  return(0);
}

static INT NLPartAssPostProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x,
                                 VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  NP_PA_NL *pa = (NP_PA_NL*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* call post routines */
  for (i=0; i<PA_NASS(pa); i++)
    if (NPANL_POST(PA_ASS(pa,i))!=NULL)
    {
      if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,0.,0.,x,d,NULL,NULL,PA_VD_g(pa),J))
        REP_ERR_RETURN (1);
      if (NPPNL_POST(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
        REP_ERR_RETURN(1);
    }

  return(0);
}

static INT NLPartAssSolution (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, INT *res)
{
  NP_PA_NL *pa = (NP_PA_NL*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* call assemble solution routines */
  for (i=0; i<PA_NASS(pa); i++)
    if (NPPNL_ASSSOL(PA_ASS(pa,i))!=NULL)
    {
      if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,0.,0.,u,NULL,NULL,NULL,PA_VD_g(pa),NULL))
        REP_ERR_RETURN(1);
      if (NPPNL_ASSSOL(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
        REP_ERR_RETURN(1);
    }

  return(0);
}

static INT NLPartAssDefect (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *s,
                            VECDATA_DESC *r, MATDATA_DESC *A, INT *res)
{
  NP_PA_NL *pa = (NP_PA_NL*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i,l;

  /* first clear skip flags of global problem */
  for (l=fl; l<=tl; l++)
    ClearVecskipFlags(GRID_ON_LEVEL(NP_MG(ass),l),s);

  /* call assemble defect routines */
  for (i=0; i<PA_NASS(pa); i++)
  {
    if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,0.,0.,s,r,NULL,NULL,PA_VD_g(pa),A))
      REP_ERR_RETURN(1);
    PP_ACTION(pp) = PARTASS_DEFECT;
    if (NPPNL_ASS(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
      REP_ERR_RETURN(1);
  }

  return(0);
}

static INT NLPartAssMatrix (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *s,
                            VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *A, INT *res)
{
  NP_PA_NL *pa = (NP_PA_NL*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* clear global matrix */
  if (dmatset(NP_MG(ass),fl,tl,ALL_VECTORS,A,0.0)!=NUM_OK) REP_ERR_RETURN (__LINE__);

  /* call assemble defect routines */
  for (i=0; i<PA_NASS(pa); i++)
  {
    if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,0.,0.,s,d,NULL,v,PA_VD_g(pa),A))
      REP_ERR_RETURN(1);
    PP_ACTION(pp) = PARTASS_MATRIX;
    if (NPPNL_ASS(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
      REP_ERR_RETURN(1);
  }

  return(0);
}

static INT NLPartAssConstruct (NP_BASE *theNP)
{
  NP_NL_ASSEMBLE *np = (NP_NL_ASSEMBLE *) theNP;

  NP_INIT(np)             = NLPartAssInit;
  NP_EXECUTE(np)          = NPNLAssembleExecute;
  NP_DISPLAY(np)          = NLPartAssDisplay;

  NPANL_PRE(np)           = NLPartAssPreProcess;
  NPANL_POST(np)          = NLPartAssPostProcess;
  NPANL_ASSSOL(np)        = NLPartAssSolution;
  NPANL_ASSDEF(np)        = NLPartAssDefect;
  NPANL_ASSMAT(np)        = NLPartAssMatrix;
  return(0);
}

/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

static INT TPartAssInit (NP_BASE *theNP, INT argc, char **argv)
{
  NP_PA_T *pa = (NP_PA_T*) theNP;
  VEC_TEMPLATE *mvt;
  INT r,i,j,nass;
  char name[NAMESIZE],buffer[VALUELEN];

  r = NPTAssembleInit(theNP,argc,argv);

  /* get name of main vector template */
  if (ReadArgvChar("m",buffer,argc,argv)!=NULL)
  {
    PrintErrorMessage('E',"NLPartAssInit","m option with main vector template not found");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  mvt = GetVectorTemplate(MGFORMAT(NP_MG(theNP)),buffer);
  if (mvt == NULL)
  {
    PrintErrorMessageF('E',"NLPartAssInit",
                       "cannot find specified vector template '%s'",buffer);
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PA_VT(pa) = mvt;

  PA_VD_g(pa) = ReadArgvVecDesc(NP_MG(theNP),"g",argc,argv);

  PA_NASS(pa) = nass = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (nass>=MAX_PA)
      {
        PrintErrorMessage('E',"NLPartAssInit","max number of part assembling numprocs exceeded");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      if (sscanf(argv[i],expandfmt(CONCAT3("ass %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"NLPartAssInit","specify a nonlinear part assembling numproc with $ass");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PA_ASS(pa,nass) = (NP_T_PARTASS*) GetNumProcByName (NP_MG(theNP),name,T_PARTASS_CLASS_NAME);
      if (PA_ASS(pa,nass) == NULL)
      {
        PrintErrorMessage('E',"NLPartAssInit",
                          "cannot find specified numerical procedure");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }

      /* next arg has to specify a valid sub vector */
      if (++i>=argc)
      {
        PrintErrorMessage('E',"NLPartAssInit",
                          "last ass option has no sub option");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      /* scan name of subvector of main vector template */
      if (sscanf(argv[i],expandfmt(CONCAT3("sub %",NAMELENSTR,"[ -~]")),name)!=1)
      {
        PrintErrorMessage('E',"NLPartAssInit","s option expected after ass option");
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      for (j=0; j<VT_NSUB(mvt); j++)
        if (strcmp(SUBV_NAME(VT_SUB(mvt,j)),name)==0)
          break;
      if (j>=VT_NSUB(mvt))
      {
        PrintErrorMessageF('E',"NLPartAssInit","name '%s' of sub template not found",name);
        REP_ERR_RETURN (NP_NOT_ACTIVE);
      }
      PA_SUB(pa,nass) = j;
      nass++;
      break;
    }
  if (nass==0)
  {
    PrintErrorMessage('E',"NLPartAssInit","specify at least one nonlinear assembling numproc with $ass");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  PA_NASS(pa) = nass;

  if ((r==NP_ACTIVE) || (r==NP_EXECUTABLE))
    return (r);
  else
    REP_ERR_RETURN (r);
}

static INT TPartAssDisplay (NP_BASE *theNP)
{
  NP_PA_T *pa = (NP_PA_T*) theNP;
  INT i;
  char text[8];

  NPTAssembleDisplay(theNP);

  if (PA_VD_g(pa)!=NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"g",ENVITEM_NAME(PA_VD_g(pa)));

  UserWriteF(DISPLAY_NP_FORMAT_SS,"vec tmplt",ENVITEM_NAME(PA_VT(pa)));

  UserWrite("\npart assembling numprocs:\n");
  for (i=0; i<PA_NASS(pa); i++)
  {
    sprintf(text,"ass%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SSS,text,ENVITEM_NAME(PA_ASS(pa,i)),SUBV_NAME(VT_SUB(PA_VT(pa),PA_SUB(pa,i))));
  }

  return (0);
}

static INT TPartAssPreProcess (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE t_p1, DOUBLE t_0,
                               DOUBLE t_m1, VECDATA_DESC *u_p1, VECDATA_DESC *u_0,
                               VECDATA_DESC *u_m1, INT *res)
{
  NP_PA_T *pa = (NP_PA_T*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  PA_DT(pa) = t_p1-t_0;

  /* call prep routines */
  for (i=0; i<PA_NASS(pa); i++)
    if (NPPT_PRE(PA_ASS(pa,i))!=NULL)
    {
      if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,t_p1,PA_DT(pa),u_p1,NULL,u_0,NULL,PA_VD_g(pa),NULL))
        REP_ERR_RETURN(1);
      if (NPPT_PRE(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
        REP_ERR_RETURN(1);
    }

  return(0);
}

static INT TPartAssPostProcess (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE t_p1, DOUBLE t_0,
                                DOUBLE t_m1, VECDATA_DESC *u_p1, VECDATA_DESC *u_0,
                                VECDATA_DESC *u_m1, INT *res)
{
  NP_PA_T *pa = (NP_PA_T*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* call post routines */
  for (i=0; i<PA_NASS(pa); i++)
    if (NPPT_POST(PA_ASS(pa,i))!=NULL)
    {
      if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,t_p1,t_p1-t_0,u_p1,NULL,u_0,NULL,PA_VD_g(pa),NULL))
        REP_ERR_RETURN(1);
      if (NPPT_POST(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
        REP_ERR_RETURN(1);
    }

  return(0);
}

static INT TPartAssInitial (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE time, VECDATA_DESC *u, INT *res)
{
  NP_PA_T *pa = (NP_PA_T*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* call assemble solution routines */
  for (i=0; i<PA_NASS(pa); i++)
  {
    if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),1.,0.,time,PA_DT(pa),u,NULL,NULL,NULL,PA_VD_g(pa),NULL))
      REP_ERR_RETURN(1);
    if (NPPT_INITIAL(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
      REP_ERR_RETURN(1);
  }

  return(0);
}

static INT TPartAssSolution (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE time, VECDATA_DESC *u, INT *res)
{
  NP_PA_T *pa = (NP_PA_T*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* call assemble solution routines */
  for (i=0; i<PA_NASS(pa); i++)
  {
    if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),0.,1.,time,PA_DT(pa),u,NULL,NULL,NULL,PA_VD_g(pa),NULL))
      REP_ERR_RETURN(1);
    if (NPPT_ASSSOL(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
      REP_ERR_RETURN(1);
  }

  return(0);
}

static INT TPartAssDefect (NP_T_ASSEMBLE *ass, INT fl, INT tl,
                           DOUBLE time, DOUBLE s_m, DOUBLE s_a, VECDATA_DESC *s,
                           VECDATA_DESC *r, MATDATA_DESC *A, INT *res)
{
  NP_PA_T *pa = (NP_PA_T*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i,l;

  /* first clear skip flags of global problem */
  for (l=fl; l<=tl; l++)
    ClearVecskipFlags(GRID_ON_LEVEL(NP_MG(ass),l),s);

  /* call assemble defect routines */
  for (i=0; i<PA_NASS(pa); i++)
  {
    if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),s_a,s_m,time,PA_DT(pa),s,r,NULL,NULL,PA_VD_g(pa),A))
      REP_ERR_RETURN(1);
    PP_ACTION(pp) = PARTASS_DEFECT;
    if (NPPT_ASS(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
      REP_ERR_RETURN(1);
  }

  return(0);
}

static INT TPartAssMatrix (NP_T_ASSEMBLE *ass, INT fl, INT tl, DOUBLE time, DOUBLE s_a,
                           VECDATA_DESC *s, VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *A, INT *res)
{
  NP_PA_T *pa = (NP_PA_T*) ass;
  PARTASS_PARAMS papa,*pp=&papa;
  INT i;

  /* clear global matrix */
  if (dmatset(NP_MG(ass),fl,tl,ALL_VECTORS,A,0.0)!=NUM_OK) REP_ERR_RETURN (__LINE__);

  /* call assemble matrix routines */
  for (i=0; i<PA_NASS(pa); i++)
  {
    if (SetPartassParams(pp,PA_VT(pa),PA_SUB(pa,i),s_a,1.,time,PA_DT(pa),s,d,NULL,v,PA_VD_g(pa),A))
      REP_ERR_RETURN(1);
    PP_ACTION(pp) = PARTASS_DEFECT;
    if (NPPT_ASS(PA_ASS(pa,i)) (PA_ASS(pa,i),fl,tl,pp,res))
      REP_ERR_RETURN(1);
  }

  return(0);
}

static INT TPartAssConstruct (NP_BASE *theNP)
{
  NP_T_ASSEMBLE *np = (NP_T_ASSEMBLE *) theNP;

  NP_INIT(np)             = TPartAssInit;
  NP_DISPLAY(np)          = TPartAssDisplay;
  NP_EXECUTE(np)          = NPTAssembleExecute;

  NPAT_PRE(np)            = TPartAssPreProcess;
  NPAT_POST(np)           = TPartAssPostProcess;
  NPAT_INITIAL(np)        = TPartAssInitial;
  NPAT_ASSSOL(np)         = TPartAssSolution;
  NPAT_ASSDEF(np)         = TPartAssDefect;
  NPAT_ASSMAT(np)         = TPartAssMatrix;

  return(0);
}

/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

INT NPNLPartAssInit (NP_BASE *theNP, INT argc, char **argv)
{
  NP_NL_PARTASS *np       = (NP_NL_PARTASS *) theNP;
  MULTIGRID *mg           = NP_MG(theNP);

  NPPNL_A(np) = ReadArgvMatDesc(mg,"A",argc,argv);
  NPPNL_x(np) = ReadArgvVecDesc(mg,"x",argc,argv);
  NPPNL_c(np) = ReadArgvVecDesc(mg,"c",argc,argv);
  NPPNL_b(np) = ReadArgvVecDesc(mg,"b",argc,argv);
  NPPNL_g(np) = ReadArgvVecDesc(mg,"g",argc,argv);
  NPPNL_t(np) = ReadArgvVecTemplateSub(MGFORMAT(mg),"part",argc,argv,&NPPNL_s(np));

  if ((NPPNL_A(np) == NULL) || (NPPNL_b(np) == NULL) || (NPPNL_x(np) == NULL) || (NPPNL_t(np) == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPNLPartAssDisplay (NP_BASE *theNP)
{
  NP_NL_PARTASS *np       = (NP_NL_PARTASS *) theNP;

  UserWrite("part description:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SSS,"vt+sub",ENVITEM_NAME(NPPNL_t(np)),SUBV_NAME(VT_SUB(NPPNL_t(np),NPPNL_s(np))));

  UserWrite("\nsymbolic user data:\n");
  if (NPPNL_A(np)!=NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(NPPNL_A(np)));
  if (NPPNL_x(np)!=NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(NPPNL_x(np)));
  if (NPPNL_c(np)!=NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(NPPNL_c(np)));
  if (NPPNL_b(np)!=NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(NPPNL_b(np)));
  if (NPPNL_g(np)!=NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"g",ENVITEM_NAME(NPPNL_g(np)));
  UserWrite("\n");

  return (0);
}

INT NPNLPartAssExecute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_NL_PARTASS *np       = (NP_NL_PARTASS *) theNP;
  INT result                  = NUM_OK;
  int level                   = CURRENTLEVEL(NP_MG(theNP));
  PARTASS_PARAMS papa,*pp=&papa;

  if (NPPNL_x(np) == NULL) {
    PrintErrorMessage('E',"NPNLAssembleExecute","no vector x");
    REP_ERR_RETURN (1);
  }
  if (NPPNL_b(np) == NULL) {
    PrintErrorMessage('E',"NPNLAssembleExecute","no vector b");
    REP_ERR_RETURN (1);
  }
  if (NPPNL_A(np) == NULL) {
    PrintErrorMessage('E',"NPNLAssembleExecute","no matrix A");
    REP_ERR_RETURN (1);
  }

  if (NPPNL_t(np)!=NULL)
  {
    if (SetPartassParams(pp,NPPNL_t(np),NPPNL_s(np),1.,0.,0.,0.,NPPNL_x(np),NPPNL_b(np),NULL,NULL,NPPNL_g(np),NPPNL_A(np)))
      REP_ERR_RETURN (1);
  }
  else
  {
    DefaultPartassParams(pp);
    PP_MD_A(pp)                     = NPPNL_A(np);
    PP_VD_s(pp)                     = NPPNL_x(np);
    PP_VD_r(pp)                     = NPPNL_b(np);
    PP_VD_gridvel(pp)       = NPPNL_g(np);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (NPPNL_PRE(np) == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no PreProcess");
      REP_ERR_RETURN (1);
    }
    if (NPPNL_PRE(np) (np,0,level,pp,&result)) {
      PrintErrorMessageF('E',"NPNLAssembleExecute","PreProcess failed, error code %d\n",result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (NPPNL_ASSSOL(np) == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no NLAssembleSolution");
      REP_ERR_RETURN (1);
    }
    if (NPPNL_ASSSOL(np) (np,0,level,pp,&result)) {
      PrintErrorMessageF('E',"NPNLAssembleExecute","NLAssembleSolution failed, error code %d\n",result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("a",argc,argv)) {
    if (NPPNL_ASS(np) == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no NLAssembleDefect");
      REP_ERR_RETURN (1);
    }
    if (NPPNL_ASS(np) (np,0,level,pp,&result)) {
      PrintErrorMessageF('E',"NPNLAssembleExecute","NLPassemble failed, error code %d\n",result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (NPPNL_POST(np) == NULL) {
      PrintErrorMessage('E',"NPNLAssembleExecute","no PostProcess");
      REP_ERR_RETURN (1);
    }
    if (NPPNL_POST(np) (np,0,level,pp,&result)) {
      PrintErrorMessageF('E',"NPNLAssembleExecute","PostProcess failed, error code %d\n",result);
      REP_ERR_RETURN (1);
    }
  }

  return(0);
}

/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

INT NPTPartAssInit (NP_BASE *theNP, INT argc, char **argv)
{
  return(NP_ACTIVE);
}

INT NPTPartAssDisplay (NP_BASE *theNP)
{
  return(0);
}

INT NPTPartAssExecute (NP_BASE *theNP, INT argc, char **argv)
{
  /* never executable */
  REP_ERR_RETURN(1);
}

/****************************************************************************/
/*D
        InitAssemble - init the file assemble.c

        SYNOPSIS:
        INT InitAssemble (void)

    PARAMETERS:
   .   void -

        DESCRIPTION:
        Create a numerical procedure class 'partass' which can call several
        assembling numprocs. 'partass' does not check any compatibility or
        consistency!

        RETURN VALUE:
        INT
   .n   __LINE__: CreateClass failed
   .n   0: ok
   D*/
/****************************************************************************/

INT InitAssemble (void)
{
  /* create partass class */
  if (CreateClass(NL_ASSEMBLE_CLASS_NAME ".nlpass",
                  sizeof(NP_NL_PARTASS), NLPartAssConstruct))
    REP_ERR_RETURN (__LINE__);

  /* create time partass class */
  if (CreateClass(T_ASSEMBLE_CLASS_NAME ".tpass",
                  sizeof(NP_T_PARTASS), TPartAssConstruct))
    REP_ERR_RETURN (__LINE__);

  /* TODO (HRR 971118): remove old part assemble num proc (2 items) */
  /* create partass class */
  if (CreateClass(NL_ASSEMBLE_CLASS_NAME ".partass",
                  sizeof(OLD_NP_NL_PARTASS), NLPartAssembleConstruct))
    REP_ERR_RETURN (__LINE__);

  /* create time partass class */
  if (CreateClass(T_ASSEMBLE_CLASS_NAME ".parttass",
                  sizeof(OLD_NP_T_PARTASS), TPartAssembleConstruct))
    REP_ERR_RETURN (__LINE__);

  return(0);
}
