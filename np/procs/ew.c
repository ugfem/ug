// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ew.c	                                                                                                        */
/*																			*/
/* Purpose:   eigenvalue solver num procs                                       */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Januar 7, 1997                                                                            */
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
#include <math.h>

#include "general.h"
#include "debug.h"
#include "ugstruct.h"
#include "devices.h"
#include "debug.h"
#include "gm.h"
#include "ugblas.h"
#include "disctools.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "formats.h"
#include "np.h"
#include "ugstruct.h"

#include "assemble.h"
#include "transfer.h"
#include "ls.h"

#include "ew.h"

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
#define VERY_SMALL 1e-10

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  NP_EW_SOLVER ew;

  NP_LINEAR_SOLVER *LS;
  NP_TRANSFER *Transfer;

  INT maxiter;
  INT baselevel;
  INT display;
  INT Orthogonalize;
  INT Quadratic;
  INT Neumann;
  INT assemble;
  INT interpolate;
  INT reset;

  VEC_SCALAR damp;

  VECDATA_DESC *r;
  VECDATA_DESC *t;
  VECDATA_DESC *q;
  MATDATA_DESC *M;

} NP_EW;

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
static INT global;

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
   NP_EW_SOLVER - type definition for linear solvers

   DESCRIPTION:
   This numproc type is used for the description of eigenvalue solvers.
   It can be called by the given interface from.
   Initializing the data is optional; it can be done with

   'INT NPEWSolverInit (NP_LINEAR_SOLVER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPEWSolverDisplay (NP_LINEAR_SOLVER *theNP);'
   'INT NPEWSolverExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   typedef struct {
        INT error_code;                     // error code
    LRESULT lresult[MAX_NUMBER_EW];     // result of linear solver
    VEC_SCALAR defect[MAX_NUMBER_EW];   // defect
        INT converged[MAX_NUMBER_EW];       // convergence flag
        INT number_of_iterations[MAX_NUMBER_EW];  // number of iterations
   } EWRESULT;

   struct np_ew_solver {

        NP_BASE base;                        // inherits base class

        // data (optional, necessary for calling the generic execute routine)
    INT nev;                             // number of eigenvectors
    VECDATA_DESC *ev[MAX_NUMBER_EW];     // eigenvectors
    DOUBLE ew[MAX_NUMBER_EW];            // eigenvalues
    NP_NL_ASSEMBLE *Assemble;            // assembling stiffness matrix
                                         // and right hand side
        VEC_SCALAR reduction;                // reduction factor
        VEC_SCALAR abslimit;                 // absolute limit for the defect

        // functions
        INT (*PreProcess)
             (struct np_ew_solver *,         // pointer to (derived) object
                  INT,                           // level
                  INT,                           // number of eigenvectors
                  VECDATA_DESC **,               // eigenvectors
                  NP_NL_ASSEMBLE *,              // matrix and right hand side
                  INT *);                        // result
    INT (*Rayleigh)
             (struct np_ew_solver *,         // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // eigenvector
                  NP_NL_ASSEMBLE *,              // matrix and right hand side
                  DOUBLE *,                      // numerator and denominator
                  DOUBLE *,                      // quotient
                  INT *);                        // result
    INT (*Solver)
             (struct np_ew_solver *,         // pointer to (derived) object
                  INT,                           // level
                  INT,                           // number of eigenvectors
                  VECDATA_DESC **,               // eigenvectors
                  DOUBLE *,                      // eigenvalues
                  NP_NL_ASSEMBLE *,              // matrix and right hand side
                  VEC_SCALAR,                    // reduction factor
                  VEC_SCALAR,                    // absolut limit for the defect
                  EWRESULT *);                   // result structure
        INT (*PostProcess)
             (struct np_ew_solver *,         // pointer to (derived) object
                  INT,                           // level
                  INT,                           // number of eigenvectors
                  VECDATA_DESC **,               // eigenvectors
                  NP_NL_ASSEMBLE *,              // matrix and right hand side
                  INT *);                        // result
   };
   typedef struct np_ew_solver NP_EW_SOLVER;
   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPEWSolverInit (NP_EW_SOLVER *np, INT argc , char **argv)
{
  INT i;
  int n;
  char *token,*names,buffer[128];

  n = 0;
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e') {
      if (sscanf(argv[i],"e %s",buffer)!=1) {
        UserWrite("Missing symbol for eigenvector in init of ew\n");
        return(NP_NOT_ACTIVE);
      }
      names =  argv[i];
      names++;
      while ((*names==' ')||(*names=='\t')) names++;
      token = strtok(names," ");
      np->ev[n] = GetVecDataDescByName(np->base.mg,token);
      if (np->ev[n] == NULL)
        np->ev[n] = CreateVecDescOfTemplate(np->base.mg,token,NULL);
      if (np->ev[n++] == NULL)
        return(NP_NOT_ACTIVE);
      token = strtok(NULL," ");
      if (token != NULL)
        if (sscanf(token,"%d",&n) != 1) {
          n = 1;
          while (token!=NULL) {
            np->ev[n] = GetVecDataDescByName(np->base.mg,token);
            if (np->ev[n] == NULL)
              np->ev[n] = CreateVecDescOfTemplate(np->base.mg,
                                                  token,NULL);
            if (np->ev[n++] == NULL)
              return(NP_NOT_ACTIVE);
            token = strtok(NULL," ");
          }
        }
    }
  np->nev = n;
  if (sc_read(np->abslimit,NP_FMT(np),np->ev[0],"abslimit",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->abslimit[i] = ABS_LIMIT;
  if (sc_read(np->reduction,NP_FMT(np),np->ev[0],"red",argc,argv))
    return(NP_ACTIVE);
  np->Assemble = (NP_NL_ASSEMBLE *)
                 ReadArgvNumProc(np->base.mg,"A",NL_ASSEMBLE_CLASS_NAME,argc,argv);

  if ((np->Assemble == NULL) || (np->nev == 0))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPEWSolverDisplay (NP_EW_SOLVER *np)
{
  INT i;

  if (np->nev > 0) UserWrite("symbolic user data:\n");
  for (i=0; i<np->nev; i++)
    if (i<10)
      UserWriteF("ev[%d]            = %-35.32s\n",
                 i,ENVITEM_NAME(np->ev[i]));
    else
      UserWriteF("ev[%d]           = %-35.32s\n",
                 i,ENVITEM_NAME(np->ev[i]));
  UserWrite("\n");

  UserWrite("configuration parameters:\n");
  if (np->Assemble != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Assemble",ENVITEM_NAME(np->Assemble));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Assemble","---");

  if (sc_disp(np->reduction,np->ev[0],"red"))
    return (1);
  if (sc_disp(np->abslimit,np->ev[0],"abslimit"))
    return (1);

  return(0);
}

INT NPEWSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EW_SOLVER *np;
  DOUBLE a[2],q;
  EWRESULT ewresult;
  INT result,level;

  np = (NP_EW_SOLVER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->Assemble == NULL) {
    PrintErrorMessage('E',"NPEWSolverExecute","no assemble num proc");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (*np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPEWSolverExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,np->nev,np->ev,np->Assemble,&result)) {
      UserWriteF("NPEWSolverExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("q",argc,argv)) {
    if (*np->Rayleigh == NULL) {
      PrintErrorMessage('E',"NPEWSolverExecute","no Rayleigh");
      return (1);
    }
    if ((*np->Rayleigh)(np,level,np->ev[0],np->Assemble,a,&q,&result)) {
      UserWriteF("NPEWSolverExecute: Rayleigh failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (*np->Solver == NULL) {
      PrintErrorMessage('E',"NPEWSolverExecute","no Solver");
      return (1);
    }
    if ((*np->Solver)(np,level,np->nev,np->ev,np->ew,np->Assemble,
                      np->abslimit,np->reduction,&ewresult)) {
      UserWriteF("NPEWSolverExecute: Solver failed, error code %d\n",
                 ewresult.error_code);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (*np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPEWSolverExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->nev,np->ev,np->Assemble,&result)) {
      UserWriteF("NPEWSolverExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }
  return(0);
}

/* tools for eingenvalue computations */

static INT SetUnsymmetric (MULTIGRID *mg, INT fl, INT tl,
                           const VECDATA_DESC *x, INT xclass)
{
  register VECTOR *v;
  register SHORT i;
  register SHORT ncomp;
  register INT vtype;
  INT lev;

  for (lev=fl; lev<=tl; lev++)
    l_setindex(GRID_ON_LEVEL(mg,lev));

  for (vtype=0; vtype<NVECTYPES; vtype++)
    if (VD_ISDEF_IN_TYPE(x,vtype)) {
      ncomp = VD_NCMPS_IN_TYPE(x,vtype);
      A_VLOOP__TYPE_CLASS(lev,fl,tl,v,mg,vtype,xclass)
      for (i=0; i<ncomp; i++)
        VVALUE(v,VD_CMP_OF_TYPE(x,vtype,i)) =
          VINDEX(v) + i * 0.3 + 0.1;
    }

    #ifdef ModelP
  if (a_vector_consistent(mg,fl,tl,x))
    return(NUM_ERROR);
    #endif

  return (NUM_OK);
}

static INT RayleighQuotient (MULTIGRID *theMG,
                             MATDATA_DESC *M, VECDATA_DESC *x,
                             VECDATA_DESC *r, VECDATA_DESC *t, DOUBLE *a)
{
  VEC_SCALAR scal1,scal2;
  INT i,bl,tl;

  bl = 0;
  tl = CURRENTLEVEL(theMG);

  if (dset(theMG,bl,tl,ON_SURFACE,t,0.0))
    return(1);

  if (dmatmul(theMG,bl,tl,ON_SURFACE,t,M,x) != NUM_OK)
    return(1);

        #ifdef ModelP
  if (a_vector_collect(theMG,bl,tl,t))
    return(1);
        #endif

  if (ddot(theMG,bl,tl,ON_SURFACE,t,x,a) != NUM_OK)
    return(1);

  if (ddot(theMG,bl,tl,ON_SURFACE,r,x,a+1) != NUM_OK)
    return(1);

  return (0);
}

static INT RayleighQuotientQ (MULTIGRID *theMG,
                              MATDATA_DESC *M, VECDATA_DESC *x,
                              VECDATA_DESC *r, VECDATA_DESC *t,
                              VECDATA_DESC *q, DOUBLE *a)
{
  GRID *theGrid;
  VEC_SCALAR scal1,scal2;
  INT i,bl,tl;

  bl = 0;
  tl = CURRENTLEVEL(theMG);

  for (i=tl-1; i>=bl; i--)
    if (StandardProject(GRID_ON_LEVEL(theMG,i),r,r))
      return(1);
  if (dset(theMG,bl,tl,ALL_VECTORS,t,0.0) != NUM_OK)
    return(1);
  if (dset(theMG,bl,tl,ALL_VECTORS,q,0.0) != NUM_OK)
    return(1);
  if (dmatmul(theMG,bl,tl,ALL_VECTORS,q,M,x) != NUM_OK)
    return(1);

  for (i=tl-1; i>=bl; i--)
    if (StandardProject(GRID_ON_LEVEL(theMG,i),q,q))
      return(1);
  if (dmatmul(theMG,bl,tl,ALL_VECTORS,t,M,q) != NUM_OK)
    return(1);
  for (i=tl-1; i>=bl; i--)
    if (StandardProject(GRID_ON_LEVEL(theMG,i),t,t))
      return(1);


  if (ddot(theMG,bl,tl,ON_SURFACE,q,q,a) != NUM_OK)
    return(1);
  if (ddot(theMG,bl,tl,ON_SURFACE,r,x,a+1) != NUM_OK)
    return(1);

  return (0);
}

static INT RayleighDefect (MULTIGRID *theMG, VECDATA_DESC *r,
                           VECDATA_DESC *t, DOUBLE rq, DOUBLE *defect)
{
  INT i,bl,tl;

  bl = 0;
  tl = CURRENTLEVEL(theMG);

  if (daxpy(theMG,bl,tl,ON_SURFACE,t,-rq,r))
    return(1);
  if (dnrm2x(theMG,bl,tl,ON_SURFACE,t,defect))
    return(1);

  return (0);
}

static INT Orthogonalize (MULTIGRID *theMG, INT level, INT m,
                          VECDATA_DESC **ev, VECDATA_DESC *b, INT display)
{
  VEC_SCALAR scal;
  DOUBLE scalP;
  INT i,j;

  for (i=0; i<m; i++) {
    if (display == PCR_FULL_DISPLAY)
      UserWriteF("%s ",ENVITEM_NAME(ev[i]));
    if (ddot(theMG,0,level,ON_SURFACE,ev[i],b,&scalP) != NUM_OK)
      return(1);
    if (display == PCR_FULL_DISPLAY)
      UserWriteF(" %lf",scalP);
    if (daxpy(theMG,0,level,ALL_VECTORS,ev[m],-scalP,ev[i]) != NUM_OK)
      return(1);
  }
  if ((display == PCR_FULL_DISPLAY) && (m > 0)) UserWrite("\n");

  return(0);
}

/****************************************************************************/
/*D
   ew - numproc for nonlinear multigrid algorithm

   DESCRIPTION:
   This numproc executes an inverse iteration for the computation of
   eigenvalues and eigenvectors of symmetric elliptic problems.

   .vb
   npinit <name> $e <sym list> [$t <tmp sym>] [$r <rhs sym>] [$M <mat sym>]
              [$d {no|red|full}] $m <maxit> [$O] [$N]
                  $red <sc double list> [$damp <sc double list>]
                  $L <linear solver> $A <assemble> $T <transfer>;
   .ve

   .  $L~<linear~solver> - linear solver num proc
   .  $A~<assemble> - assemble num proc
   .  $T~<transfer> - transfer num proc
   .  $e~<sym~list> - list of symbols for the eigenvectors
   .  $r~<rhs~sym> - symbol for the right hnd side vector
   .  $t~<tmp~sym> - symbol for a tempory vector
   .  $m~<maxit> - maximal number of multigrid cylces
   .  $O - orthogonalize using right hand side
   .  $N - set the first eigenvector = 1 for Neumann problems
   .  $d - no, reduced or full display
   .  $red~<sc~double~list> - reduction factors for each component
   .  $damp~<sc~double~list> - damping factors for each component
   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$a] [$r];'

   .  $r - reset solution
   .  $i - interpolate by the standard interpolation
   .  $a - assemble the matrix

   EXAMPLE:
   .vb
   npcreate smooth $c ilu;
   npinit smooth;

   npcreate baseiter $c lu;
   npinit baseiter;

   npcreate basesolver $c ls;
   npinit basesolver $red 0.001 $m 1 $I baseiter;

   npcreate transfer $c transfer;
   npinit transfer;

   npcreate lmgc $c lmgc;
   npinit lmgc $S smooth smooth basesolver $T transfer $n1 1 $n2 1;

   npcreate mgs $c cg;
   npinit mgs $m 8 $red 0.00001 $I lmgc $display no;

   npcreate assemble $c ewassemble;
   npinit assemble $w weight $ew;

   npcreate ew $c ew;
   npinit ew $e ew0 ew1 $m 10 $red 0.000001
          $A assemble $T transfer $L mgs $display full $N;

   clear weight $v 1.0 $a;
   npexecute ew $r $a;
   .ve
   D*/
/****************************************************************************/

static INT EWPreProcess (NP_EW_SOLVER *theNP, INT level, INT nev,
                         VECDATA_DESC **ev, NP_NL_ASSEMBLE *Assemble,
                         INT *result)
{
  NP_EW *np;
  INT i,bl;

  np = (NP_EW *) theNP;

  bl = 0;
  for (i=1; i<nev; i++)
    if (AllocVDFromVD(theNP->base.mg,bl,level,ev[0],&ev[i]))
      NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theNP->base.mg,bl,level,ev[0],&np->r))
    NP_RETURN(1,result[0]);
  if (AllocMDFromVD(theNP->base.mg,bl,level,ev[0],ev[0],&np->M))
    NP_RETURN(1,result[0]);
  if (Assemble->PreProcess != NULL)
    if ((*Assemble->PreProcess)(Assemble,bl,level,ev[0],result))
      return(1);
  if (np->reset)
    for (i=0; i<nev; i++)
      if (SetUnsymmetric(theNP->base.mg,bl,level,ev[i],EVERY_CLASS))
        NP_RETURN(1,result[0]);
  np->reset = 0;
  if (np->interpolate) {
    if (np->Transfer->PreProcessSolution != NULL)
      if ((*np->Transfer->PreProcessSolution)
            (np->Transfer,bl,level,ev[0],result))
        return(1);
    for (i=0; i<nev; i++)
      if ((*np->Transfer->InterpolateNewVectors)
            (np->Transfer,bl,level,ev[i],result))
        return(1);
  }
  if (np->assemble) {
    if (AllocVDFromVD(theNP->base.mg,bl,level,ev[0],&np->t))
      NP_RETURN(1,result[0]);
    if ((*Assemble->NLAssembleMatrix)(Assemble,bl,level,
                                      ev[0],np->r,np->t,np->M,result))
      return(1);
    if (FreeVD(theNP->base.mg,bl,level,np->t)) NP_RETURN(1,result[0]);
    if (np->LS->PreProcess != NULL)
      if ((*np->LS->PreProcess)(np->LS,level,ev[0],np->r,np->M,
                                &np->baselevel,result))
        return(1);
    np->assemble = 0;
  }
  if (np->Quadratic)
    for (i=bl; i<=level; i++)
      AssembleTotalDirichletBoundary(GRID_ON_LEVEL(theNP->base.mg,i),
                                     np->M,ev[0],np->r);

  return (0);
}

static INT Rayleigh (NP_EW_SOLVER *theNP, INT level,
                     VECDATA_DESC *ev, NP_NL_ASSEMBLE *Assemble,
                     DOUBLE *a, DOUBLE *q, INT *result)
{
  NP_EW *np;

  np = (NP_EW *) theNP;
  if (np->M == NULL) NP_RETURN(1,result[0]);
  if (np->r == NULL) NP_RETURN(1,result[0]);
  if (np->t == NULL) NP_RETURN(1,result[0]);
  if ((*Assemble->NLAssembleDefect)(Assemble,0,level,ev,np->r,np->M,result))
    return(1);

        #ifdef ModelP
  if (a_vector_collect(theNP->base.mg,0,level,np->r))
    return(1);
        #endif

  IFDEBUG(np,5)
  UserWriteF("r\n");
  PrintVector(GRID_ON_LEVEL(theNP->base.mg,level),np->r,3,3);
  UserWriteF("ev\n");
  PrintVector(GRID_ON_LEVEL(theNP->base.mg,level),ev,3,3);
  ENDDEBUG

  if (np->Quadratic) {
    if (AllocVDFromVD(theNP->base.mg,0,level,ev,&np->q))
      NP_RETURN(1,result[0]);
    if (RayleighQuotientQ(theNP->base.mg,np->M,ev,np->r,np->t,np->q,a))
      NP_RETURN(1,result[0]);
    if (FreeVD(theNP->base.mg,0,level,np->q))
      NP_RETURN(1,result[0]);
  }
  else if (RayleighQuotient(theNP->base.mg,np->M,ev,np->r,np->t,a))
    NP_RETURN(1,result[0]);

  PRINTDEBUG(np,1,("a0 %f a1 %f ",a[0],a[1]));

  if (ABS(a[1]) <= ABS(a[0]) * VERY_SMALL)
    NP_RETURN(1,result[0]);
  *q = a[0] / a[1];

  return (0);
}

static INT EWSolver (NP_EW_SOLVER *theNP, INT level, INT nev,
                     VECDATA_DESC **ev, DOUBLE *ew, NP_NL_ASSEMBLE *Assemble,
                     VEC_SCALAR abslimit, VEC_SCALAR reduction,
                     EWRESULT *ewresult)
{
  NP_EW *np;
  MULTIGRID *theMG;
  INT i,PrintID,j,bl,iter;
  char text[DISPLAY_WIDTH+4];
  VEC_SCALAR defect, defect2reach,scal;
  DOUBLE a[2],rq,s;

  np = (NP_EW *) theNP;
  theMG = theNP->base.mg;
  bl = 0;

  if (Assemble->NLAssembleDefect == NULL)
    NP_RETURN(1,ewresult->error_code);
  ewresult->error_code = 0;
  i = 0;
  if (np->Neumann) {                   /* set ev[0] = 1 */
    if (dset(theMG,bl,level,ON_SURFACE,ev[0],1.0))
      NP_RETURN(1,ewresult->error_code);
    if (np->Neumann == 2)
      SetUnsymmetric(theMG,bl,level,ev[0],EVERY_CLASS);
    if ((*Assemble->NLAssembleDefect)(Assemble,bl,level,ev[0],np->r,np->M,
                                      &ewresult->error_code))
      return(1);
            #ifdef ModelP
    if (a_vector_collect(theMG,bl,level,np->r))
      NP_RETURN(1,ewresult->error_code);
            #endif
    if (ddot(theMG,0,level,ON_SURFACE,ev[0],np->r,a+1) != NUM_OK)
      NP_RETURN(1,ewresult->error_code);
    if (dscal(theMG,0,level,ALL_VECTORS,ev[0],1.0/sqrt(a[1])))
      NP_RETURN(1,ewresult->error_code);
    ew[0] = 0.0;
    i++;
  }
  while (i < nev) {
    if (np->display == PCR_FULL_DISPLAY)
      UserWriteF("%s:\n",ENVITEM_NAME(ev[i]));

    /* orthogonalize iteration vector */
    if (AllocVDFromVD(theMG,bl,level,ev[0],&np->t))
      NP_RETURN(1,ewresult->error_code);
    if (np->Orthogonalize) {
      if ((*Assemble->NLAssembleDefect)(Assemble,bl,level,
                                        ev[i],np->t,np->M,
                                        &ewresult->error_code))
        return(1);
      if (ew[i] < 0.0)
        if (dscal(theMG,bl,level,ALL_VECTORS,np->t,-1.0))
          NP_RETURN(1,ewresult->error_code);
    }
    else {
      if (dset(theMG,0,level,ON_SURFACE,np->t,0.0))
        NP_RETURN(1,ewresult->error_code);
      if (dmatmul (theMG,0,level,ON_SURFACE,np->t,np->M,ev[i]) != NUM_OK)
        NP_RETURN(1,ewresult->error_code);
    }
            #ifdef ModelP
    if (a_vector_collect(theMG,bl,level,np->t))
      NP_RETURN(1,ewresult->error_code);
            #endif

    if (Orthogonalize(theMG,level,i,ev,np->t,np->display))
      NP_RETURN(1,ewresult->error_code);
    if (Rayleigh(&(np->ew),level,ev[i],Assemble,a,&rq,
                 &ewresult->error_code))
      return(1);
    if (np->display == PCR_FULL_DISPLAY)
      UserWriteF("Rayleigh quotient %lf\n",rq);
    if (np->Orthogonalize) {
      if (ABS(a[1]) <= VERY_SMALL) NP_RETURN(1,ewresult->error_code);
      s = 1.0 / sqrt(ABS(a[1]));
    }
    else {
      if (a[0] <= 0.0) NP_RETURN(1,ewresult->error_code);
      s = 1.0 / sqrt(a[0]);
    }
    if (dscal(theMG,bl,level,ALL_VECTORS,ev[i],s) != NUM_OK)
      NP_RETURN(1,ewresult->error_code);
    if (dscal(theMG,bl,level,ALL_VECTORS,np->r,s) != NUM_OK)
      NP_RETURN(1,ewresult->error_code);
    if (dscal(theMG,bl,level,ALL_VECTORS,np->t,s) != NUM_OK)
      NP_RETURN(1,ewresult->error_code);
    CenterInPattern(text,DISPLAY_WIDTH," inverse iteration ",'%',"\n");
    if (PreparePCR(np->r,np->display,text,&PrintID))
      NP_RETURN(1,ewresult->error_code);
    if (RayleighDefect(theMG,np->r,np->t,rq,defect))
      NP_RETURN(1,ewresult->error_code);
    if (sc_mul(defect2reach,defect,reduction,np->t))
      NP_RETURN(1,ewresult->error_code);
    if (DoPCR(PrintID,defect,PCR_CRATE))
      NP_RETURN(1,ewresult->error_code);
    for (iter=0; iter<np->maxiter; iter++) {
      if (sc_cmp(defect,defect2reach,np->t))
        break;
      if (sc_cmp(defect,abslimit,np->t))
        break;
      /* orthogonalize iteration vector */
      if (np->Orthogonalize) {
        if ((*Assemble->NLAssembleDefect)(Assemble,bl,level,
                                          ev[i],np->t,
                                          np->M,&ewresult->error_code))
          return(1);
        if (ew[i] < 0.0)
          if (dscal(theMG,bl,level,ALL_VECTORS,np->t,-1.0))
            NP_RETURN(1,ewresult->error_code);
      }
      else {
        if (dset(theMG,0,level,ON_SURFACE,np->t,0.0))
          NP_RETURN(1,ewresult->error_code);
        if (dmatmul (theMG,0,level,ON_SURFACE,np->t,np->M,ev[i])
            != NUM_OK) NP_RETURN(1,ewresult->error_code);
      }

            #ifdef ModelP
      if (a_vector_collect(theMG,bl,level,np->t))
        NP_RETURN(1,ewresult->error_code);
            #endif
      if (Orthogonalize(theMG,level,i,ev,np->t,np->display))
        NP_RETURN(1,ewresult->error_code);
      if (Rayleigh(&(np->ew),level,ev[i],Assemble,a,&rq,
                   &ewresult->error_code))
        return(1);
      if (dscal(theMG,0,level,ALL_VECTORS,np->r,rq) != NUM_OK)
        NP_RETURN(1,ewresult->error_code);
      /* solve */
      if (np->Quadratic) {
        if (dcopy (theMG,bl,level,ALL_VECTORS,np->t,ev[i]) != NUM_OK)
          NP_RETURN(1,ewresult->error_code);
        if ((*np->Transfer->ProjectSolution)
              (np->Transfer,bl,level,ev[i],&ewresult->error_code))
          NP_RETURN(1,ewresult->error_code);
        if ((*np->Transfer->ProjectSolution)
              (np->Transfer,bl,level,np->r,&ewresult->error_code))
          NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Defect)(np->LS,level,np->t,np->r,np->M,
                              &ewresult->error_code))
          return (1);
        if ((*np->LS->Residuum)(np->LS,level,level,np->t,np->r,np->M,
                                &ewresult->lresult[i]))
          NP_RETURN(1,ewresult->error_code);

        PRINTDEBUG(np,2,("res1 %f\n",
                         ewresult->lresult[0].last_defect[0]));

        if ((*np->LS->Solver)(np->LS,level,np->t,np->r,np->M,
                              abslimit,reduction,
                              &ewresult->lresult[i]))
          return (1);
        if ((*np->Transfer->ProjectSolution)
              (np->Transfer,bl,level,np->t,&ewresult->error_code))
          NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Defect)(np->LS,level,ev[i],np->t,np->M,
                              &ewresult->error_code))
          return (1);
        if ((*np->LS->Residuum)(np->LS,level,level,ev[i],np->t,np->M,
                                &ewresult->lresult[i]))
          return (1);

        PRINTDEBUG(np,2,("res2 %f\n",
                         ewresult->lresult[0].last_defect[0]));

        if ((*np->LS->Solver)(np->LS,level,ev[i],np->t,np->M,
                              abslimit,reduction,
                              &ewresult->lresult[i]))
          return (1);
        if (FreeVD(theMG,bl,level,np->t))
          NP_RETURN(1,ewresult->error_code);
      }
      else {
        if (FreeVD(theMG,bl,level,np->t))
          NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Defect)(np->LS,level,ev[i],np->r,np->M,
                              &ewresult->error_code))
          NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Residuum)(np->LS,level,level,ev[i],np->r,np->M,
                                &ewresult->lresult[i]))
          NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Solver)(np->LS,level,ev[i],np->r,np->M,
                              abslimit,reduction,
                              &ewresult->lresult[i]))
          NP_RETURN(1,ewresult->error_code);
      }
      if (AllocVDFromVD(theMG,bl,level,ev[0],&np->t))
        NP_RETURN(1,ewresult->error_code);
      if (Rayleigh(&(np->ew),level,ev[i],Assemble,a,&rq,
                   &ewresult->error_code))
        return(1);
      if (np->display == PCR_FULL_DISPLAY)
        UserWriteF("Rayleigh quotient %lf\n",rq);
      if (np->Orthogonalize) {
        if (ABS(a[1]) <= VERY_SMALL) NP_RETURN(1,ewresult->error_code);
        s = 1.0 / sqrt(ABS(a[1]));
      }
      else {
        if (a[0] <= 0.0) NP_RETURN(1,ewresult->error_code);
        s = 1.0 / sqrt(a[0]);
      }
      if (dscal(theMG,bl,level,ALL_VECTORS,ev[i],s) != NUM_OK)
        NP_RETURN(1,ewresult->error_code);
      if (dscal(theMG,bl,level,ALL_VECTORS,np->r,s) != NUM_OK)
        NP_RETURN(1,ewresult->error_code);
      if (dscal(theMG,bl,level,ALL_VECTORS,np->t,s) != NUM_OK)
        NP_RETURN(1,ewresult->error_code);
      /* print defect */
      if (RayleighDefect(theMG,np->r,np->t,rq,defect))
        NP_RETURN(1,ewresult->error_code);
      if (FreeVD(theMG,bl,level,np->t))
        NP_RETURN(1,ewresult->error_code);
      if (DoPCR(PrintID,defect,PCR_CRATE))
        NP_RETURN(1,ewresult->error_code);
    }

    /* print average and finish */
    if (DoPCR(PrintID,defect,PCR_AVERAGE))
      NP_RETURN(1,ewresult->error_code);
    if (PostPCR(PrintID,":ew:avg"))
      NP_RETURN(1,ewresult->error_code);
    ewresult->number_of_iterations[i] = iter + 1;
    ewresult->converged[i] = !(iter == np->maxiter);
    ew[i++] = rq;
  }

  return (0);
}

static INT EWPostProcess (NP_EW_SOLVER *theNP, INT level, INT nev,
                          VECDATA_DESC **ev, NP_NL_ASSEMBLE *Assemble,
                          INT *result)
{
  NP_EW *np;
  INT i,bl;

  np = (NP_EW *) theNP;
  bl = 0;

  for (i=1; i<nev; i++)
    if (FreeVD(theNP->base.mg,bl,level,ev[i])) NP_RETURN(1,result[0]);
  if (FreeVD(theNP->base.mg,bl,level,np->r)) NP_RETURN(1,result[0]);
  if (FreeMD(theNP->base.mg,bl,level,np->M)) NP_RETURN(1,result[0]);
  if (Assemble->PostProcess != NULL)
    if ((*Assemble->PostProcess)(Assemble,bl,level,ev[0],
                                 np->r,np->M,result))
      return(1);
  for (i=0; i<nev; i++)
    if ((*np->Transfer->ProjectSolution)
          (np->Transfer,bl,level,ev[i],result))
      NP_RETURN(1,result[0]);
  if (np->Transfer->PostProcess != NULL)
    if ((*np->Transfer->PostProcess)(np->Transfer,0,level,
                                     ev[0],np->r,np->M,result))
      NP_RETURN(1,result[0]);
  if (np->LS->PostProcess != NULL)
    if ((*np->LS->PostProcess)(np->LS,level,ev[0],np->r,np->M,result))
      NP_RETURN(1,result[0]);

  return (0);
}

static INT EWInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EW *np;
  INT i;

  np = (NP_EW *) theNP;

  np->interpolate = 0;
  np->reset = 1;
  np->LS = (NP_LINEAR_SOLVER *)
           ReadArgvNumProc(theNP->mg,"L",LINEAR_SOLVER_CLASS_NAME,argc,argv);
  if (np->LS == NULL)
    return(NP_NOT_ACTIVE);
  np->Transfer = (NP_TRANSFER *)
                 ReadArgvNumProc(theNP->mg,"T",TRANSFER_CLASS_NAME,argc,argv);

  np->M = ReadArgvMatDesc(theNP->mg,"M",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  if (sc_read(np->damp,NP_FMT(np),np->r,"damp",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->damp[i] = 1.0;
  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    return(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->baselevel = 0;
  if (ReadArgvOption("O",argc,argv))
    np->Orthogonalize = 1;
  else
    np->Orthogonalize = 0;
  if (ReadArgvOption("Q",argc,argv))
    np->Quadratic = 1;
  else
    np->Quadratic = 0;
  if (ReadArgvOption("N",argc,argv)) {
    if (ReadArgvOption("S",argc,argv))
      np->Neumann = 2;
    else
      np->Neumann = 1;
  }
  else
    np->Neumann = 0;
  if (np->Neumann)
    np->Orthogonalize = 1;
  if (ReadArgvOption("na",argc,argv))
    np->assemble = 0;
  else
    np->assemble = 1;

  return(NPEWSolverInit(&np->ew,argc,argv));
}

static INT EWDisplay (NP_BASE *theNP)
{
  NP_EW *np;

  np = (NP_EW *) theNP;
  NPEWSolverDisplay(&np->ew);

  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  if (np->LS != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(np->LS));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L","---");
  if (np->Transfer != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(np->Transfer));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->r != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->q != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"q",ENVITEM_NAME(np->q));
  if (np->M != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"M",ENVITEM_NAME(np->M));
  if (sc_disp(np->damp,np->r,"damp"))
    return(1);
  if (np->Orthogonalize)
    UserWrite("\nuse right hand side for orthogolization\n");
  else
    UserWrite("\nuse left hand side for orthogolization\n");
  if (np->Quadratic)
    UserWrite("\nuse quadratic stiffness matrix\n");
  if (np->Neumann)
    UserWrite("\nNeumann boundary\n");

  return(0);
}

static INT EWExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EW *np;
  EWRESULT ewresult;
  INT i,result,level;

  np = (NP_EW *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->ew.Assemble == NULL) {
    PrintErrorMessage('E',"EWExecute","no assemble num proc");
    return (1);
  }
  np->assemble = ReadArgvOption("a",argc,argv);
  np->interpolate = ReadArgvOption("i",argc,argv);
  np->reset = ReadArgvOption("r",argc,argv);
  global = ReadArgvOption("g",argc,argv);
  if (np->reset && np->interpolate) {
    PrintErrorMessage('E',"EWExecute",
                      "Only one option $r or $i can be specified.\n");
    return(1);
  }
  if ((*np->ew.PreProcess)(&np->ew,level,np->ew.nev,np->ew.ev,
                           np->ew.Assemble,&result)) {
    UserWriteF("EWExecute: PreProcess failed, error code %d\n",result);
    return (1);
  }
  if ((*np->ew.Solver)(&np->ew,level,np->ew.nev,np->ew.ev,
                       np->ew.ew,np->ew.Assemble,
                       np->ew.abslimit,np->ew.reduction,&ewresult)) {
    UserWriteF("NPEWSolverExecute: Solver failed, error code %d\n",
               ewresult.error_code);
    return (1);
  }

  if ((*np->ew.PostProcess)(&np->ew,level,np->ew.nev,np->ew.ev,
                            np->ew.Assemble,&result)) {
    UserWriteF("EWExecute: PostProcess failed, error code %d\n",result);
    return (1);
  }
  if (ChangeStructDir(":ew")==NULL)
    return (1);
  for (i=0; i<np->ew.nev; i++) {
    if (np->display > PCR_NO_DISPLAY)
      UserWriteF("  ew%d = %10.5lf\n",i,np->ew.ew[i]);
    if (SetStringValue(ENVITEM_NAME(np->ew.ev[i]),np->ew.ew[i]))
      return (1);
  }
  if (ChangeStructDir(":")==NULL)
    return (1);

  return(0);
}

static INT EWConstruct (NP_BASE *theNP)
{
  NP_EW *np;

  theNP->Init = EWInit;
  theNP->Display = EWDisplay;
  theNP->Execute = EWExecute;

  np = (NP_EW *) theNP;
  np->ew.PreProcess = EWPreProcess;
  np->ew.Rayleigh = Rayleigh;
  np->ew.Solver = EWSolver;
  np->ew.PostProcess = EWPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitEW	- Init this file

   SYNOPSIS:
   INT InitEW ();

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

INT InitEW ()
{
  INT i;

  if (CreateClass(EW_SOLVER_CLASS_NAME ".ew",sizeof(NP_EW),EWConstruct))
    return (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;
  if (MakeStruct(":ew")) return(__LINE__);
  if (MakeStruct(":ew:avg")) return(__LINE__);

  return (0);
}
