// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ls.c	                                                                                                        */
/*																			*/
/* Purpose:   linear solver num procs                                           */
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
#include <time.h>
#include <math.h>

#include "devices.h"
#include "general.h"
#include "debug.h"
#include "ugstruct.h"
#include "gm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"

#include "iter.h"
#include "ls.h"

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
#define MAX_RESTART 20
#ifdef ModelP
#include "ppif.h"
#define CSTART()    clock_start=CurrentTime()
#define CSTOP(t,c)  t+=(CurrentTime()-clock_start);c++
#else
#define CSTART()    clock_start=clock()
#define CSTOP(t,c)  t+=((double)(clock()-clock_start))/((double)CLOCKS_PER_SEC);c++
#endif

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct np_ls
{
  NP_LINEAR_SOLVER ls;

  NP_ITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;

  VECDATA_DESC *c;

  INT (*Prepare)
    (struct np_ls *,                         /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*Update)
    (struct np_ls *,                         /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* correction vector               */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*Close)
    (struct np_ls *,                         /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_ls NP_LS;

typedef struct
{
  NP_LS ls;

  DOUBLE rho;
  VECDATA_DESC *p;
  VECDATA_DESC *t;

} NP_CG;

typedef struct
{
  NP_LINEAR_SOLVER ls;

  NP_ITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;
  INT restart;

  DOUBLE sp;
  VEC_SCALAR weight;
  VECDATA_DESC *p;
  VECDATA_DESC *pp;
  VECDATA_DESC *t;
  VECDATA_DESC *h1;
  VECDATA_DESC *h2;
  VECDATA_DESC *h3;

} NP_CR;

typedef struct
{
  NP_LINEAR_SOLVER ls;

  INT maxiter;
  INT baselevel;
  INT display;
  INT restart;

  DOUBLE rbr;
  VECDATA_DESC *rb;
  VECDATA_DESC *p;
  VECDATA_DESC *pb;
  VECDATA_DESC *h;

} NP_BCG;

typedef struct
{
  NP_LINEAR_SOLVER ls;

  NP_ITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;
  INT restart;

  DOUBLE rho, omega;
  VEC_SCALAR weight;
  VECDATA_DESC *r;
  VECDATA_DESC *p;
  VECDATA_DESC *v;
  VECDATA_DESC *s;
  VECDATA_DESC *t;
  VECDATA_DESC *q;

} NP_BCGS;

typedef struct
{
  NP_LINEAR_SOLVER ls;

  NP_ITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;
  INT restart;

  DOUBLE rho, omega;
  VEC_SCALAR weight;
  VECDATA_DESC *c;
  VECDATA_DESC *r;
  VECDATA_DESC *p;
  VECDATA_DESC *s;
  VECDATA_DESC *t;
  VECDATA_DESC *q;
  VECDATA_DESC *w;
  VECDATA_DESC *v[MAX_RESTART+1];
} NP_GMRES;

typedef struct
{
  NP_LINEAR_SOLVER ls;

  INT maxiter;
  INT baselevel;
  INT display;
  INT restart;

  VECDATA_DESC *r;
  VECDATA_DESC *p;
  VECDATA_DESC *h;
  VECDATA_DESC *d;

} NP_SQCG;

struct np_ldcs
{
  NP_LINEAR_SOLVER ls;

  NP_LINEAR_SOLVER *linsol;

  VECDATA_DESC *b;                       /* tmp defect                                  */
  VECDATA_DESC *c;                       /* tmp corr                                            */
  MATDATA_DESC *DC;                      /* defect correction matrix        */
  INT maxiter;
  INT display;
  INT baselevel;
};
typedef struct np_ldcs NP_LDCS;

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
   NP_LINEAR_SOLVER - type definition for linear solvers

   DESCRIPTION:
   This numproc type is used for the description of linear solvers.
   It can be called by the given interface from a nonlinearsolver.
   Initializing the data is optional; it can be done with

   'INT NPLinearSolverInit (NP_LINEAR_SOLVER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPLinearSolverDisplay (NP_LINEAR_SOLVER *theNP);'
   'INT NPLinearSolverExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   typedef struct {
        INT error_code;                     // error code
        INT converged;                      // error code
        VEC_SCALAR first_defect;            // first defect
        VEC_SCALAR last_defect;             // last defect
        INT number_of_linear_iterations;    // number of iterations
   } LRESULT;

   struct np_linear_solver {
        NP_BASE base;                        // inherits base class

        // data (optional, necessary for calling the generic execute routine)
    VECDATA_DESC *x;                     // solution
    VECDATA_DESC *b;                     // defect
    MATDATA_DESC *A;                     // matrix
        VEC_SCALAR reduction;                // reduction factor
        VEC_SCALAR abslimit;                 // absolute limit for the defect

        // functions
        INT (*PreProcess)
             (struct np_linear_solver *,     // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *,                         // baselevel used by the solver
                  INT *);                        // result
    INT (*Defect)                        // b := b - Ax
             (struct np_linear_solver *,     // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
    INT (*Residuum)                      // computes norm of the defect
             (struct np_linear_solver *,     // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  LRESULT *);                    // result structure
    INT (*Solver)                        // b := b - Ax
             (struct np_linear_solver *,     // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  VEC_SCALAR,                    // reduction factor
                  VEC_SCALAR,                    // absolute limit for the defect
                  LRESULT *);                    // result structure
        INT (*PostProcess)
             (struct np_linear_solver *,     // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
   };
   typedef struct np_linear_solver NP_LINEAR_SOLVER;
   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPLinearSolverInit (NP_LINEAR_SOLVER *np, INT argc , char **argv)
{
  INT i;

  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"b",argc,argv);
  if (sc_read(np->abslimit,np->x,"abslimit",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->abslimit[i] = ABS_LIMIT;
  if (sc_read(np->reduction,NULL,"red",argc,argv))
    return(NP_ACTIVE);

  if ((np->x == NULL) || (np->b == NULL) || (np->A == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPLinearSolverDisplay (NP_LINEAR_SOLVER *np)
{
  if ((np->x != NULL) || (np->b != NULL) || (np->A != NULL)) {
    UserWrite("symbolic user data:\n");
    if (np->A != NULL)
      UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
    if (np->x != NULL)
      UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
    if (np->b != NULL)
      UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->b));
    UserWrite("\n");
  }
  UserWrite("configuration parameters:\n");
  if (np->x != NULL)
    if (sc_disp(np->reduction,np->x,"red"))
      REP_ERR_RETURN (1);
  if (sc_disp(np->abslimit,np->x,"abslimit"))
    REP_ERR_RETURN (1);

  return(0);
}

INT NPLinearSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LINEAR_SOLVER *np;
  LRESULT lresult;
  INT result,level,bl;

  np = (NP_LINEAR_SOLVER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no vector x");
    REP_ERR_RETURN (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no vector b");
    REP_ERR_RETURN (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPLinearSolverExecute","no matrix A");
    REP_ERR_RETURN (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no PreProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PreProcess)(np,level,np->x,np->b,np->A,&bl,&result)) {
      UserWriteF("NPLinearSolverExecute: PreProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("d",argc,argv)) {
    if (np->Defect == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Defect");
      REP_ERR_RETURN (1);
    }
    if ((*np->Defect)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPLinearSolverExecute: Defect failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("r",argc,argv)) {
    if (np->Residuum == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Residuum");
      REP_ERR_RETURN (1);
    }
    if ((*np->Residuum)(np,bl,level,np->x,np->b,np->A,&lresult)) {
      UserWriteF("NPLinearSolverExecute: Residuum failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->Solver == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no Solver");
      REP_ERR_RETURN (1);
    }
    if ((*np->Solver)(np,level,np->x,np->b,np->A,
                      np->abslimit,np->reduction,&lresult)) {
      UserWriteF("NPLinearSolverExecute: Solver failed, error code %d\n",
                 lresult.error_code);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPLinearSolverExecute","no PostProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PostProcess)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("NPLinearSolverExecute: PostProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }
  return(0);
}

/* tools for linear solvers */

static INT LinearSolverInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    REP_ERR_RETURN(NP_NOT_ACTIVE);

  np->display = ReadArgvDisplay(argc,argv);
  np->Iter = (NP_ITER *)
             ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);
  if (np->Iter == NULL)
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->baselevel = 0;
  np->c = ReadArgvVecDesc(theNP->mg,"c",argc,argv);

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT LinearSolverDisplay (NP_BASE *theNP)
{
  NP_LS *np;

  np = (NP_LS *) theNP;
  NPLinearSolverDisplay(&np->ls);

  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->Iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");
  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->c));

  return (0);
}

static INT LinearSolverPreProcess (NP_LINEAR_SOLVER *theNP, INT level,
                                   VECDATA_DESC *x, VECDATA_DESC *b,
                                   MATDATA_DESC *A,
                                   INT *baselevel, INT *result)
{
  NP_LS *np;

  /* store passed XXXDATA_DESCs */
  NPLS_A(theNP) = A;
  NPLS_x(theNP) = x;
  NPLS_b(theNP) = b;

  np = (NP_LS *) theNP;
  if (np->Iter->PreProcess != NULL)
    if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result))
      REP_ERR_RETURN(1);
  np->baselevel = MIN(*baselevel,level);

  return(0);
}

static INT LinearDefect (NP_LINEAR_SOLVER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *result)
{
  NP_LS *np;

  np = (NP_LS *) theNP;
  if (s_dmatmul_minus(theNP->base.mg,np->baselevel,level,b,A,x,EVERY_CLASS)
      != NUM_OK) NP_RETURN(1,result[0]);
  return (0);
}

static INT LinearResiduum (NP_LINEAR_SOLVER *theNP, INT bl, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           LRESULT *lresult)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

        #ifdef ModelP
  if (a_vector_collect(theNP->base.mg,bl,level,b)) NP_RETURN(1,lresult->error_code);
        #endif
  if (s_eunorm(theNP->base.mg,bl,level,b,lresult->last_defect)) NP_RETURN(1,lresult->error_code);

  return(0);
}

static INT LinearSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_LS *np;
  VEC_SCALAR defect2reach;
  INT i,bl,PrintID;
  char text[DISPLAY_WIDTH+4];
  double ti;
  int ii;
#ifdef ModelP
  double clock_start;
#else
  clock_t clock_start;
#endif

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  np = (NP_LS *) theNP;
  bl = np->baselevel;
  if (np->Iter->Iter == NULL) NP_RETURN(1,lresult->error_code);
  if (np->Update == NULL) NP_RETURN(1,lresult->error_code);
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->c)) NP_RETURN(1,lresult->error_code);
  if (np->Prepare != NULL)
    if ((*np->Prepare)(np,level,x,&lresult->error_code)) REP_ERR_RETURN (1);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY)
    if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  CSTART(); ti=0; ii=0;
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY)
    if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;
    if (l_dset(GRID_ON_LEVEL(theNP->base.mg,level),np->c,EVERY_CLASS,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    if ((*np->Iter->Iter)(np->Iter,level,np->c,b,A,&lresult->error_code)) REP_ERR_RETURN (1);
    if ((*np->Update)(np,level,x,np->c,b,A,&lresult->error_code)) REP_ERR_RETURN (1);
    if (LinearResiduum(theNP,bl,level,x,b,A,lresult)) REP_ERR_RETURN(1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  if (!lresult->converged)
    lresult->number_of_linear_iterations=i;
  FreeVD(theNP->base.mg,bl,level,np->c);
  if (np->Close != NULL)
    if ((*np->Close)(np,level,&lresult->error_code))
      REP_ERR_RETURN (1);
  CSTOP(ti,ii);

  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
    UserWriteF("LS  : L=%2d N=%2d TSOLVE=%10.4lg TIT=%10.4lg\n",level,
               lresult->number_of_linear_iterations,ti,
               ti/lresult->number_of_linear_iterations);
  }

  return (0);
}

static INT LinearSolverPostProcess (NP_LINEAR_SOLVER *theNP,
                                    INT level,
                                    VECDATA_DESC *x, VECDATA_DESC *b,
                                    MATDATA_DESC *A,
                                    INT *result)
{
  NP_LS *np;

  np = (NP_LS *) theNP;

  if (np->Iter->PostProcess == NULL)
    return(0);
  return((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result));
}

/****************************************************************************/
/*D
   ls - numproc for linear solvers

   DESCRIPTION:
   This numproc executes a linear solver: it performs an iteration
   up to convergence.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n   nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata
   .n   if only a single value is specified, this will be used for all components

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate pre $c ilu;           npinit pre;
   npcreate post $c ilu;          npinit post;
   npcreate base $c ilu;          npinit base $n 3;
   npcreate basesolver $c ls;     npinit basesolver $red 0.001 $I base;
   npcreate transfer $c transfer; npinit transfer;
   npcreate lmgc $c lmgc;         npinit lmgc $S pre post basesolver $T transfer;
   npcreate mgs $c ls;            npinit mgs $A MAT $x sol $b rhs
                                          $red 0.00001 $I lmgc $d full;
   npexecute mgs $i $d $r $s $p;
   .ve
   D*/
/****************************************************************************/

static INT LSUpdate (NP_LS *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *c,
                     VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  if (a_daxpy (theNP->ls.base.mg,theNP->baselevel,level,
               x,EVERY_CLASS,Factor_One,c) != NUM_OK) NP_RETURN(1,result[0]);

  return(0);
}

static INT LSConstruct (NP_BASE *theNP)
{
  NP_LS *np;

  theNP->Init = LinearSolverInit;
  theNP->Display = LinearSolverDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_LS *) theNP;
  np->ls.PreProcess = LinearSolverPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = LinearSolver;
  np->ls.PostProcess = LinearSolverPostProcess;

  np->Prepare = NULL;
  np->Update = LSUpdate;
  np->Close = NULL;

  return(0);
}

/****************************************************************************/
/*D
   cg - numproc for the conjugate gradient method

   DESCRIPTION:
   This numproc executes a conjugate gradient step. It is preconditioned
   by an iteration numproc, e. g. a multi grid cycle or a smoother.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
       [$p <con>] [$t <tmp>];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus
   .  $p~<con> - conjugate vector
   .  $t~<tmp> - temporaty vector
   .  $r~<restart> - restart index

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT CGInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  np->p = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);

  return (LinearSolverInit(theNP,argc,argv));
}

static INT CGDisplay (NP_BASE *theNP)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  LinearSolverDisplay(theNP);
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->p != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));

  return (0);
}

static INT CGPrepare (NP_LS *theNP, INT level, VECDATA_DESC *x, INT *result)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  if (AllocVDFromVD(theNP->ls.base.mg,theNP->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (a_dset(theNP->ls.base.mg,theNP->baselevel,level,np->p,EVERY_CLASS,0.0)
      != NUM_OK) NP_RETURN(1,result[0]);
  np->rho = 1.0;

  return(0);
}

static INT CGUpdate (NP_LS *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *c,
                     VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_CG *np;
  MULTIGRID *theMG;
  VEC_SCALAR scal;
  DOUBLE lambda;
  INT ncomp,j;

  np = (NP_CG *) theNP;
  theMG = theNP->ls.base.mg;
  ncomp = VD_NCOMP(x);
  if (AllocVDFromVD(theMG,theNP->baselevel,level,x,&np->t))
    NP_RETURN(1,result[0]);
  if (a_dset(theMG,theNP->baselevel,level,np->t,EVERY_CLASS,0.0) != NUM_OK)
    NP_RETURN(1,result[0]);
  for (j=theNP->baselevel; j<=level; j++)
    if (l_dmatmul(GRID_ON_LEVEL(theMG,j),np->t,EVERY_CLASS,A,c,EVERY_CLASS)
        !=NUM_OK) NP_RETURN(1,result[0]);
  if (a_daxpy(theMG,theNP->baselevel,level,b,EVERY_CLASS,Factor_One,np->t))
    NP_RETURN(1,result[0]);
  if (s_ddot(theMG,theNP->baselevel,level,c,b,scal) !=NUM_OK)
    NP_RETURN(1,result[0]);
  lambda = 0.0;
  for (j=0; j<ncomp; j++) lambda += scal[j];
  for (j=0; j<ncomp; j++) scal[j] = lambda / np->rho;
  np->rho = lambda;
  if (a_dscale(theMG,theNP->baselevel,level,np->p,EVERY_CLASS,scal)
      != NUM_OK) NP_RETURN(1,result[0]);
  if (a_daxpy (theMG,theNP->baselevel,level,np->p,EVERY_CLASS,Factor_One,c)
      != NUM_OK) NP_RETURN(1,result[0]);
  if (a_dset(theMG,theNP->baselevel,level,np->t,EVERY_CLASS,0.0) != NUM_OK)
    NP_RETURN(1,result[0]);
  for (j=theNP->baselevel; j<=level; j++)
    if (l_dmatmul (GRID_ON_LEVEL(theMG,j),np->t,EVERY_CLASS,
                   A,np->p,EVERY_CLASS) != NUM_OK) NP_RETURN(1,result[0]);
  if (s_ddot (theMG,theNP->baselevel,level,np->t,np->p,scal) != NUM_OK)
    NP_RETURN(1,result[0]);
  lambda = 0.0;
  for (j=0; j<ncomp; j++) lambda += scal[j];
  if (lambda == 0.0) NP_RETURN(1,result[0]);
  for (j=0; j<ncomp; j++) scal[j] = np->rho / lambda;
  if (a_daxpy(theMG,theNP->baselevel,level,x,EVERY_CLASS,scal,np->p)
      != NUM_OK) NP_RETURN(1,result[0]);
  for (j=0; j<ncomp; j++) scal[j] = - np->rho / lambda;
  if (a_daxpy (theMG,theNP->baselevel,level,b,EVERY_CLASS,scal,np->t)
      != NUM_OK) NP_RETURN(1,result[0]);
  FreeVD(theNP->ls.base.mg,theNP->baselevel,level,np->t);
  if (theNP->display == PCR_FULL_DISPLAY)
    UserWriteF("      rho %-.4g \n",np->rho);

  return(0);
}

static INT CGClose (NP_LS *theNP, INT level, INT *result)
{
  NP_CG *np;

  np = (NP_CG *) theNP;
  FreeVD(theNP->ls.base.mg,theNP->baselevel,level,np->p);

  return(0);
}

static INT CGConstruct (NP_BASE *theNP)
{
  NP_LS *np;

  theNP->Init = CGInit;
  theNP->Display = CGDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_LS *) theNP;
  np->ls.PreProcess = LinearSolverPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = LinearSolver;
  np->ls.PostProcess = LinearSolverPostProcess;

  np->Prepare = CGPrepare;
  np->Update = CGUpdate;
  np->Close = CGClose;

  return(0);
}

/****************************************************************************/
/*D
   cr - numproc for the conjugate residuum method

   DESCRIPTION:
   This numproc executes a conjugate residuum step. It is preconditioned
   by an iteration numproc, e. g. a multi grid cycle or a smoother.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
       [$p <con>] [$t <tmp>]
           [$r <restart>] [$w <sc double list>];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $m~<maxit> - maximal number of iterations
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus
   .  $p~<con> - conjugate vector
   .  $t~<tmp> - temporaty vector
   .  $r~<restart> - restart index
   .  $w~<sc~double~list> - weighting factor

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT CRInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_CR *np;
  INT i;

  np = (NP_CR *) theNP;
  if (sc_read (np->weight,NULL,"weight",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] = 1.0;
  np->p  = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->pp = ReadArgvVecDesc(theNP->mg,"pp",argc,argv);
  np->h1 = ReadArgvVecDesc(theNP->mg,"h1",argc,argv);
  np->h2 = ReadArgvVecDesc(theNP->mg,"h2",argc,argv);
  np->h3 = ReadArgvVecDesc(theNP->mg,"h3",argc,argv);
  np->t  = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("r",&(np->restart),argc,argv))
    np->restart = 0;
  if (np->restart<0) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->Iter = (NP_ITER *) ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);
  np->baselevel = 0;

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT CRDisplay (NP_BASE *theNP)
{
  NP_CR *np;

  np = (NP_CR *) theNP;
  NPLinearSolverDisplay(&np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"r",(int)np->restart);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->Iter != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");
  if (np->display == PCR_NO_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->p != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  if (np->pp != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"pp",ENVITEM_NAME(np->pp));
  if (np->h1 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h1",ENVITEM_NAME(np->h1));
  if (np->h2 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h2",ENVITEM_NAME(np->h2));
  if (np->h3 != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h3",ENVITEM_NAME(np->h3));
  if (np->t != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h3",ENVITEM_NAME(np->t));
  if (np->p != NULL) if (sc_disp(np->weight,np->p,"weight")) REP_ERR_RETURN (1);

  return (0);
}

static INT CRPreProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_CR *np;

  np = (NP_CR *) theNP;
  if (np->Iter!=NULL)
    if (np->Iter->PreProcess != NULL)
      if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);
  np->baselevel = MIN(*baselevel,level);

  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->pp)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->t)) NP_RETURN(1,result[0]);

  return(0);
}

static INT CRPostProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_CR *np;

  np = (NP_CR *) theNP;
  FreeVD(np->ls.base.mg,np->baselevel,level,np->p);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->pp);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->t);

  if (np->Iter!=NULL)
  {
    if (np->Iter->PostProcess == NULL) return(0);
    return((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result));
  }
  else
    return (0);

  return(0);
}

static void PrintEunorm (MULTIGRID *theMG, VECDATA_DESC *v, char *name)
{
  DOUBLE eu;

  s_eunorm(theMG,0,0,v,&eu);
  UserWriteF("EUNORM(%s): %f\n",name,(float)eu);

  return;
}

static INT CRSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_CR *np;
  VEC_SCALAR defect2reach,scal,alpha0, alpha1;
  INT i,j,bl,PrintID,restart;
  char text[DISPLAY_WIDTH+4];
  DOUBLE s,t;

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  /* prepare */
  np = (NP_CR *) theNP;
  bl = np->baselevel;
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->h1)) NP_RETURN(1,lresult->error_code);
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->h2)) NP_RETURN(1,lresult->error_code);
  if (AllocVDFromVD(theNP->base.mg,bl,level,x,&np->h3)) NP_RETURN(1,lresult->error_code);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY) if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY) if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;

  if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  if (s_dset(theNP->base.mg,np->baselevel,level,np->p,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  if (np->Iter!=NULL)
  {
    if ((*np->Iter->Iter)(np->Iter,level,np->p,np->h1,A,&lresult->error_code)) REP_ERR_RETURN (1);
  }
  else
  {
    if (s_dcopy(theNP->base.mg,np->baselevel,level,np->p,np->h1)!= NUM_OK) REP_ERR_RETURN (1);
  }
  if (s_dset(theNP->base.mg,np->baselevel,level,np->pp,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  if (s_dset(theNP->base.mg,np->baselevel,level,np->t,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
  np->sp = 1.0; restart = 0;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;

    /* update x, b */
    if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->h2,A,np->p,EVERY_CLASS)) REP_ERR_RETURN (1);
    if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,np->h2)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    if (s_dset(theNP->base.mg,np->baselevel,level,np->h3,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    if (np->Iter!=NULL)
    {
      if ((*np->Iter->Iter)(np->Iter,level,np->h3,np->h2,A,&lresult->error_code)) REP_ERR_RETURN (1);
    }
    else
    {
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h3,np->h2)!= NUM_OK) REP_ERR_RETURN (1);
    }
    if (s_ddot (theNP->base.mg,np->baselevel,level,np->h1,np->h3,scal)!=NUM_OK) REP_ERR_RETURN (1);
    s=0.0; for (j=0; j<VD_NCOMP(x); j++) s += np->weight[j]*scal[j];
    if (s_ddot (theNP->base.mg,np->baselevel,level,b,np->h3,scal)!=NUM_OK) REP_ERR_RETURN (1);
    t=0.0; for (j=0; j<VD_NCOMP(x); j++) t += np->weight[j]*scal[j];
    if (ABS(s)>1e-3*ABS(t) || (np->restart>0 && s!=0.0))
    {
      for (j=0; j<VD_NCOMP(x); j++) scal[j] = t/s;
      if (s_daxpy (theNP->base.mg,np->baselevel,level,x,scal,np->p)!= NUM_OK) REP_ERR_RETURN (1);
      for (j=0; j<VD_NCOMP(x); j++) scal[j] = -scal[j];
      if (s_daxpy (theNP->base.mg,np->baselevel,level,b,scal,np->h1)!= NUM_OK) REP_ERR_RETURN (1);
    }
    else
    {
      restart = 1;
    }

    /* update p */
    if ((np->restart>0 && i%np->restart==0) || restart==1)
    {
      restart = 0;
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dset(theNP->base.mg,np->baselevel,level,np->p,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (np->Iter!=NULL)
      {
        if ((*np->Iter->Iter)(np->Iter,level,np->p,np->h1,A,&lresult->error_code)) REP_ERR_RETURN (1);
      }
      else
      {
        if (s_dcopy(theNP->base.mg,np->baselevel,level,np->p,np->h1)!= NUM_OK) REP_ERR_RETURN (1);
      }
      if (s_dset(theNP->base.mg,np->baselevel,level,np->pp,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dset(theNP->base.mg,np->baselevel,level,np->t,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      np->sp = 1.0;
    }
    else
    {
      if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->h2,A,np->h3,EVERY_CLASS)) REP_ERR_RETURN (1);
      if (s_ddot (theNP->base.mg,np->baselevel,level,np->h2,np->h3,scal)!=NUM_OK) REP_ERR_RETURN (1);
      t=0.0; for (j=0; j<VD_NCOMP(x); j++) t += np->weight[j]*scal[j];
      for (j=0; j<VD_NCOMP(x); j++) alpha0[j] = -t/s;
      if (s_ddot (theNP->base.mg,np->baselevel,level,np->h2,np->t,scal)!=NUM_OK) REP_ERR_RETURN (1);
      t=0.0; for (j=0; j<VD_NCOMP(x); j++) t += np->weight[j]*scal[j];
      for (j=0; j<VD_NCOMP(x); j++) alpha1[j] = -t/np->sp;
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->t,np->h3)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->h1,np->pp)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->pp,np->p)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->p,np->h3)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,alpha0,np->pp)!= NUM_OK) REP_ERR_RETURN (1);
      if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,alpha1,np->h1)!= NUM_OK) REP_ERR_RETURN (1);
      np->sp = s;
    }

    /* redisuum */
    if (LinearResiduum(theNP,bl,level,x,b,A,lresult)) REP_ERR_RETURN (1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  FreeVD(theNP->base.mg,bl,level,np->h1);
  FreeVD(theNP->base.mg,bl,level,np->h2);
  FreeVD(theNP->base.mg,bl,level,np->h3);
  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
  }

  return (0);
}

static INT CRConstruct (NP_BASE *theNP)
{
  NP_CR *np;

  theNP->Init = CRInit;
  theNP->Display = CRDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_CR *) theNP;
  np->ls.PreProcess = CRPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = CRSolver;
  np->ls.PostProcess = CRPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   bcg - numproc for the bi conjugate gradient method

   DESCRIPTION:
   This numproc executes the bi-conjugate gradient method.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
       [$p <con>] [$pb <p-bar>] [$rb <r-bar>] [$h <help>]
           [$r <restart>];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $m~<maxit> - maximal number of iterations
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus
   .  $p~<con> - conjugate vector
   .  $pb~<p-bar> - temporaty vector
   .  $rb~<r-bar> - temporaty vector
   .  $h~<tmp> - temporaty vector
   .  $r~<restart> - restart index

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT BCGInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_BCG *np;

  np = (NP_BCG *) theNP;
  np->p  = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->pb = ReadArgvVecDesc(theNP->mg,"pb",argc,argv);
  np->rb = ReadArgvVecDesc(theNP->mg,"rb",argc,argv);
  np->h = ReadArgvVecDesc(theNP->mg,"h",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("r",&(np->restart),argc,argv))
    np->restart = 0;
  if (np->restart<0) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->baselevel = 0;

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT BCGDisplay (NP_BASE *theNP)
{
  NP_BCG *np;

  np = (NP_BCG *) theNP;
  NPLinearSolverDisplay(&np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"r",(int)np->restart);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->display == PCR_NO_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->p != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  if (np->pb != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"pp",ENVITEM_NAME(np->pb));
  if (np->rb != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h1",ENVITEM_NAME(np->rb));
  if (np->h != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h",ENVITEM_NAME(np->h));

  return (0);
}

static INT BCGPreProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_BCG *np;

  np = (NP_BCG *) theNP;

  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->pb)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->rb)) NP_RETURN(1,result[0]);

  return(0);
}

static INT BCGPostProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_BCG *np;

  np = (NP_BCG *) theNP;
  FreeVD(np->ls.base.mg,np->baselevel,level,np->p);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->pb);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->rb);

  return(0);
}

static INT BCGSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_BCG *np;
  VEC_SCALAR defect2reach,alpha,alpha_m,beta;
  INT i,j,PrintID,restart;
  char text[DISPLAY_WIDTH+4];
  DOUBLE sigma,rbr_new;

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  /* prepare */
  np = (NP_BCG *) theNP;
  if (AllocVDFromVD(theNP->base.mg,level,level,x,&np->h)) NP_RETURN(1,lresult->error_code);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY) if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY) if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;

  /* go */
  restart = 1;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;

    /* restart ? */
    if ((np->restart>0 && i%np->restart==0) || restart)
    {
      if (s_dset(theNP->base.mg,np->baselevel,level,np->p,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dset(theNP->base.mg,np->baselevel,level,np->pb,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->rb,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      np->rbr = 1.0;
      restart = 0;
    }

    /* update x, b */
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,b,np->rb,Factor_One,&rbr_new)!=NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) beta[j]=rbr_new/np->rbr;np->rbr = rbr_new;
    if (s_dscale (theNP->base.mg,np->baselevel,level,np->p,beta)) REP_ERR_RETURN (1);
    if (s_dscale (theNP->base.mg,np->baselevel,level,np->pb,beta)) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,Factor_One,b)!= NUM_OK) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->pb,Factor_One,np->rb)!= NUM_OK) REP_ERR_RETURN (1);
    if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->h,A,np->p,EVERY_CLASS)) REP_ERR_RETURN (1);
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->h,np->pb,Factor_One,&sigma)!=NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) {alpha[j]=np->rbr/sigma; alpha_m[j]=-np->rbr/sigma;}
    if (s_daxpy (theNP->base.mg,np->baselevel,level,x,alpha,np->p)!= NUM_OK) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,b,alpha_m,np->h)!= NUM_OK) REP_ERR_RETURN (1);
    if (s_dtpmatmul_set(theNP->base.mg,np->baselevel,level,np->h,A,np->pb,EVERY_CLASS)) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->rb,alpha_m,np->h)!= NUM_OK) REP_ERR_RETURN (1);

    /* redisuum */
    if (LinearResiduum(theNP,np->baselevel,level,x,b,A,lresult)) REP_ERR_RETURN (1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  FreeVD(theNP->base.mg,level,level,np->h);
  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
  }

  return (0);
}

static INT BCGConstruct (NP_BASE *theNP)
{
  NP_BCG *np;

  theNP->Init = BCGInit;
  theNP->Display = BCGDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_BCG *) theNP;
  np->ls.PreProcess = BCGPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = BCGSolver;
  np->ls.PostProcess = BCGPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   bcgs - numproc for the bi cg stab method

   DESCRIPTION:
   This numproc executes the bi-conjugate gradient method.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
       [$p <con>] [$t <tmp>]
           [$R <restart>] [$w <sc double list>];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus
   .  $p~<con> - conjugate vector
   .  $t~<tmp> - temporaty vector
   .  $R~<restart> - restart index
   .  $w~<sc~double~list> - weighting factor

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT BCGSInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_BCGS *np;
  INT i;

  np = (NP_BCGS *) theNP;
  if (sc_read (np->weight,NULL,"weight",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] = 1.0;
  for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] *= np->weight[i];
  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  np->p = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->v = ReadArgvVecDesc(theNP->mg,"v",argc,argv);
  np->s = ReadArgvVecDesc(theNP->mg,"s",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->q = ReadArgvVecDesc(theNP->mg,"q",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("R",&(np->restart),argc,argv))
    np->restart = 0;
  if (np->restart<0)
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->baselevel = 0;
  np->Iter = (NP_ITER *) ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT BCGSDisplay (NP_BASE *theNP)
{
  NP_BCGS *np;

  np = (NP_BCGS *) theNP;
  NPLinearSolverDisplay(&np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"R",(int)np->restart);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->Iter != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");
  if (np->display == PCR_NO_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->r != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->p != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  if (np->v != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"v",ENVITEM_NAME(np->v));
  if (np->s != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(np->s));
  if (np->t != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->q != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"q",ENVITEM_NAME(np->q));
  if (np->p != NULL) if (sc_disp(np->weight,np->p,"weight")) REP_ERR_RETURN (1);

  return (0);
}

static INT BCGSPreProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_BCGS *np;

  np = (NP_BCGS *) theNP;

  np->baselevel = MIN(*baselevel,level);

  if (np->Iter!=NULL)
    if (np->Iter->PreProcess != NULL)
      if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);

  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->r)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->v)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->s)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->t)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->q)) NP_RETURN(1,result[0]);

  return(0);
}

static INT BCGSPostProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_BCGS *np;

  np = (NP_BCGS *) theNP;
  FreeVD(np->ls.base.mg,np->baselevel,level,np->r);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->p);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->v);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->s);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->t);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->q);

  if (np->Iter!=NULL)
  {
    if (np->Iter->PostProcess == NULL) return(0);
    return((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result));
  }
  else
    return (0);

  return(0);
}

static INT BCGSSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_BCGS *np;
  VEC_SCALAR defect2reach,scal;
  INT i,j,PrintID,restart;
  char text[DISPLAY_WIDTH+4];
  DOUBLE alpha,rho_new,beta,tt;
  double ti;
  int ii;
#ifdef ModelP
  double clock_start;
#else
  clock_t clock_start;
#endif

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  /* prepare */
  np = (NP_BCGS *) theNP;

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY) if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);

  CSTART(); ti=0; ii=0;

  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY) if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;

  /* go */
  restart = 1;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;

    /* restart ? */
    if ((np->restart>0 && i%np->restart==0) || restart)
    {
      if (s_dset(theNP->base.mg,np->baselevel,level,np->p,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dset(theNP->base.mg,np->baselevel,level,np->v,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->r,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      alpha = np->rho = np->omega = 1.0;
      restart = 0;
    }

    /* update x, b */
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,b,np->r,np->weight,&rho_new)!=NUM_OK) REP_ERR_RETURN (1);
    beta=rho_new*alpha/np->rho/np->omega;
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=beta;
    if (s_dscale (theNP->base.mg,np->baselevel,level,np->p,scal)) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,Factor_One,b)!= NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-beta*np->omega;
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,scal,np->v)!= NUM_OK) REP_ERR_RETURN (1);
    if (np->Iter!=NULL)
    {
      if (s_dset(theNP->base.mg,np->baselevel,level,np->q,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->s,np->p)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->p,A,&lresult->error_code)) REP_ERR_RETURN (1);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->p,np->s)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->v,A,np->q,EVERY_CLASS)) REP_ERR_RETURN (1);
            #ifdef ModelP
      if (a_vector_collect(theNP->base.mg,np->baselevel,level,np->v)
          != NUM_OK)
        NP_RETURN(1,lresult->error_code);
            #endif
      if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->v,np->r,np->weight,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      alpha = rho_new/alpha;
      for (j=0; j<VD_NCOMP(x); j++) scal[j]=alpha;
      if (s_daxpy (theNP->base.mg,np->baselevel,level,x,scal,np->q)!= NUM_OK) REP_ERR_RETURN (1);
    }
    else
    {
            #ifdef ModelP
      if (a_vector_consistent(theNP->base.mg,np->baselevel,level,np->p)
          != NUM_OK)
        NP_RETURN(1,lresult->error_code);
            #endif
      if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->v,A,np->p,EVERY_CLASS)) REP_ERR_RETURN (1);
            #ifdef ModelP
      if (a_vector_collect(theNP->base.mg,np->baselevel,level,np->v)
          != NUM_OK)
        NP_RETURN(1,lresult->error_code);
            #endif
      if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->v,np->r,np->weight,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      alpha = rho_new/alpha;
      for (j=0; j<VD_NCOMP(x); j++) scal[j]=alpha;
      if (s_daxpy (theNP->base.mg,np->baselevel,level,x,scal,np->p)!= NUM_OK) REP_ERR_RETURN (1);
    }
    lresult->number_of_linear_iterations++;
    if (s_dcopy(theNP->base.mg,np->baselevel,level,np->s,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-alpha;
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->s,scal,np->v)!= NUM_OK) REP_ERR_RETURN (1);
    if (LinearResiduum(theNP,np->baselevel,level,x,np->s,A,lresult))
      REP_ERR_RETURN (1);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      if (s_dcopy(theNP->base.mg,np->baselevel,level,b,np->s) != NUM_OK)
        NP_RETURN(1,lresult->error_code);
      lresult->converged = 1;
      if (np->display > PCR_NO_DISPLAY)
        if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE))
          NP_RETURN(1,lresult->error_code);
      break;
    }
    if (np->Iter!=NULL)
    {
      if (s_dset(theNP->base.mg,np->baselevel,level,np->q,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->t,np->s)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->s,A,&lresult->error_code)) REP_ERR_RETURN (1);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->s,np->t)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    }
    else
    {
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->q,np->s)!= NUM_OK) NP_RETURN(1,lresult->error_code);
            #ifdef ModelP
      if (a_vector_consistent(theNP->base.mg,np->baselevel,level,np->q)
          != NUM_OK)
        NP_RETURN(1,lresult->error_code);
            #endif
    }
    if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->t,A,np->q,EVERY_CLASS)) REP_ERR_RETURN (1);
        #ifdef ModelP
    if (a_vector_collect(theNP->base.mg,np->baselevel,level,np->t)
        != NUM_OK)
      NP_RETURN(1,lresult->error_code);
        #endif
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->t,np->t,np->weight,&tt)!=NUM_OK) REP_ERR_RETURN (1);
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->s,np->t,np->weight,&(np->omega))!=NUM_OK) REP_ERR_RETURN (1);
    PRINTDEBUG(np,2,("tt %f omega %f\n",tt,np->omega));
    if (tt!=0.0)
      np->omega /= tt;
    else
      np->omega /= 1.0E-20;
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=np->omega;
    if (s_daxpy (theNP->base.mg,np->baselevel,level,x,scal,np->q)!= NUM_OK) REP_ERR_RETURN (1);
    if (s_dcopy(theNP->base.mg,np->baselevel,level,b,np->s)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-np->omega;
    if (s_daxpy (theNP->base.mg,np->baselevel,level,b,scal,np->t)!= NUM_OK) REP_ERR_RETURN (1);
    np->rho = rho_new;

    /* redisuum */
    if (LinearResiduum(theNP,np->baselevel,level,x,b,A,lresult)) REP_ERR_RETURN (1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    lresult->number_of_linear_iterations++;
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      break;
    }
  }
  CSTOP(ti,ii);

  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
    UserWriteF("BCGS: L=%2d N=%2d TSOLVE=%10.4lg TIT=%10.4lg\n",level,
               lresult->number_of_linear_iterations,ti,
               ti/lresult->number_of_linear_iterations);
  }

  return (0);
}

static INT BCGSConstruct (NP_BASE *theNP)
{
  NP_BCGS *np;

  theNP->Init = BCGSInit;
  theNP->Display = BCGSDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_BCGS *) theNP;
  np->ls.PreProcess = BCGSPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = BCGSSolver;
  np->ls.PostProcess = BCGSPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   gmres - numproc for the gmres method

   DESCRIPTION:
   This numproc executes the generalized mimimum residual method.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
       [$p <con>] [$t <tmp>]
           [$R <restart>] [$w <sc double list>];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus
   .  $p~<con> - conjugate vector
   .  $t~<tmp> - temporaty vector
   .  $R~<restart> - restart index
   .  $w~<sc~double~list> - weighting factor

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT GMRESInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_GMRES *np;
  INT i;

  np = (NP_GMRES *) theNP;
  if (sc_read (np->weight,NULL,"weight",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] = 1.0;
  for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] *= np->weight[i];
  np->c = ReadArgvVecDesc(theNP->mg,"c",argc,argv);
  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  np->p = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->s = ReadArgvVecDesc(theNP->mg,"s",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->q = ReadArgvVecDesc(theNP->mg,"q",argc,argv);
  np->w = ReadArgvVecDesc(theNP->mg,"w",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("R",&(np->restart),argc,argv))
    np->restart = 0;
  if (np->restart<0)
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  for (i=0; i<=MAX_RESTART; i++) np->v[i] = NULL;
  np->display = ReadArgvDisplay(argc,argv);
  np->baselevel = 0;
  np->Iter = (NP_ITER *)
             ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT GMRESDisplay (NP_BASE *theNP)
{
  NP_GMRES *np;
  INT i;

  np = (NP_GMRES *) theNP;
  NPLinearSolverDisplay(&np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"R",(int)np->restart);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->Iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");
  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->c));
  if (np->r != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->p != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  for (i=0; i<=MAX_RESTART; i++) if (np->v[i] != NULL)
      if (i<10)
        UserWriteF("v[%d]            = %-35.32s\n",
                   i,ENVITEM_NAME(np->v[i]));
      else
        UserWriteF("v[%d]           = %-35.32s\n",
                   i,ENVITEM_NAME(np->v[i]));
  if (np->s != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(np->s));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->q != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"q",ENVITEM_NAME(np->q));
  if (np->w != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"w",ENVITEM_NAME(np->w));
  if (np->p != NULL)
    if (sc_disp(np->weight,np->p,"weight")) REP_ERR_RETURN (1);

  return (0);
}

static INT GMRESPreProcess (NP_LINEAR_SOLVER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b,
                            MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_GMRES *np;
  INT i;
  np = (NP_GMRES *) theNP;

  np->baselevel = MIN(*baselevel,level);

  if (np->Iter!=NULL)
    if (np->Iter->PreProcess != NULL)
      if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result))
        REP_ERR_RETURN(1);

  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->c))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->r))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->p))
    NP_RETURN(1,result[0]);
  for (i=0; i<=np->restart; i++)
    if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->v[i]))
      NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->s))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->t))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->q))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->w))
    NP_RETURN(1,result[0]);

  return(0);
}

static INT GMRESPostProcess (NP_LINEAR_SOLVER *theNP, INT level,
                             VECDATA_DESC *x, VECDATA_DESC *b,
                             MATDATA_DESC *A, INT *result)
{
  NP_GMRES *np;
  INT i;

  np = (NP_GMRES *) theNP;
  FreeVD(np->ls.base.mg,np->baselevel,level,np->c);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->r);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->p);
  for (i=0; i<=np->restart; i++)
    FreeVD(np->ls.base.mg,np->baselevel,level,np->v[i]);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->s);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->t);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->q);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->w);

  if (np->Iter!=NULL) {
    if (np->Iter->PostProcess == NULL) return(0);
    return((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result));
  }
  else
    return (0);

  return(0);
}

#ifdef ModelP
static INT l_vector_makeinconsistent (GRID *g, const VECDATA_DESC *x)
{
  VECTOR *v;
  INT vc,i,type,mask,n,m;
  const SHORT *Comp;

  if (VD_IS_SCALAR(x)) {
    mask = VD_SCALTYPEMASK(x);
    vc = VD_SCALCMP(x);
    for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
      if ((mask & VDATATYPE(v)) &&
          (DDD_InfoPriority(PARHDR(v)) != PrioMaster))
        VVALUE(v,vc) = 0.0;
  }
  else
    for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v)) {
      type = VTYPE(v);
      n = VD_NCMPS_IN_TYPE(x,type);
      if (n == 0) continue;
      if (DDD_InfoPriority(PARHDR(v)) != PrioMaster) continue;
      Comp = VD_CMPPTR_OF_TYPE(x,type);
      for (i=0; i<n; i++)
        VVALUE(v,Comp[i]) = 0.0;
    }

  return(NUM_OK);
}

static INT a_vector_makeinconsistent (MULTIGRID *mg, INT fl, INT tl,
                                      const VECDATA_DESC *x)
{
  INT level;

  for (level=fl; level<=tl; level++)
    if (l_vector_makeinconsistent(GRID_ON_LEVEL(mg,level),x))
      REP_ERR_RETURN(NUM_ERROR);

  return (NUM_OK);
}
#endif

static INT GMRESSolver (NP_LINEAR_SOLVER *theNP, INT level,
                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                        VEC_SCALAR abslimit, VEC_SCALAR reduction,
                        LRESULT *lresult)
{
  NP_GMRES *np;
  VEC_SCALAR defect2reach,rnorm,scal;
  INT i,bl,PrintID;
  char text[DISPLAY_WIDTH+4];
  double ti;
  int ii;
  DOUBLE H[MAX_RESTART+1][MAX_RESTART];
  DOUBLE Hsq[MAX_RESTART*MAX_RESTART];
  DOUBLE s[MAX_RESTART+1];
  DOUBLE cs[MAX_RESTART];
  DOUBLE sn[MAX_RESTART];
  DOUBLE y[MAX_RESTART];
  DOUBLE lambda;
  DOUBLE tol;
  INT k,j;
  INT it,i1,i2,ncomp,*result;
  MULTIGRID *theMG;
    #ifdef ModelP
  double clock_start;
    #else
  clock_t clock_start;
    #endif
  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++) {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }
  result = &lresult->error_code;
  np = (NP_GMRES *) theNP;
  bl = np->baselevel;
  if (np->Iter->Iter == NULL)
    NP_RETURN(1,lresult->error_code);
  theMG = theNP->base.mg;
  ncomp = VD_NCOMP(x);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY)
    if (PreparePCR(x,np->display,text,&PrintID))
      NP_RETURN(1,lresult->error_code);
  CSTART(); ti=0; ii=0;
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b))
    NP_RETURN(1,lresult->error_code);

  /* Compute the tolerance */
  tol = 0.0;
  for (j=0; j<ncomp; j++) tol += defect2reach[j]*defect2reach[j];
  tol = sqrt(tol);

  if (np->display > PCR_NO_DISPLAY)
    if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE))
      NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;

  /* GMRES outer loop */
  for (it=0; it<np->maxiter; it++) {
    if (lresult->converged) break;

    if (s_ddot(theMG,bl,level,b,b,rnorm) !=NUM_OK)
      NP_RETURN(1,result[0]);
    lambda = 0.0;
    for (j=0; j<ncomp; j++) lambda += rnorm[j];
    printf("res %12.8f\n",lambda);
    if (s_dcopy(theMG,bl,level,np->s,b)!= NUM_OK)
      NP_RETURN(1,lresult->error_code);

    /* Solve preconditioner for the initial residual */
    if (a_dset(theMG,bl,level,np->r,EVERY_CLASS,0.0) != NUM_OK)
      NP_RETURN(1,result[0]);
    if ((*np->Iter->Iter)(np->Iter,level,np->r,np->s,A,
                          &lresult->error_code))
      REP_ERR_RETURN (1);
    /* form the norm of the initial residual */
                #ifdef ModelP
    if (s_dcopy(theMG,bl,level,np->s,np->r)!= NUM_OK)
      NP_RETURN(1,lresult->error_code);
    if (a_vector_makeinconsistent(theMG,bl,level,np->s) != NUM_OK)
      NP_RETURN(1,lresult->error_code);
    if (s_ddot(theMG,bl,level,np->s,np->r,rnorm) !=NUM_OK)
      NP_RETURN(1,result[0]);
                #else
    if (s_ddot(theMG,bl,level,np->r,np->r,rnorm) !=NUM_OK)
      NP_RETURN(1,result[0]);
                #endif

    lambda = 0.0;
    for (j=0; j<ncomp; j++) lambda += rnorm[j];
    if (lambda <= 0.0) NP_RETURN(1,result[0]);
    lambda = 1.0/sqrt(lambda);

    /* copy the initial residual into v[0] */
    if (s_dcopy(theNP->base.mg,bl,level,np->v[0],np->r) != NUM_OK)
      NP_RETURN(1,lresult->error_code);

    /* scale v[0] = r/norm(r) */
    for (j=0; j<ncomp; j++) rnorm[j] = lambda;
    if (a_dscale(theMG,bl,level,np->v[0],EVERY_CLASS,rnorm) != NUM_OK)
      NP_RETURN(1,result[0]);

    /* form s = norm(r)*e1 */
    for (j=0; j<=MAX_RESTART; j++) s[j] = 0.0;
    s[0] = 1.0/lambda;

    /* GMRES inner loop */
    for (i=0; i<np->restart; i++) {
      printf("#####################################\n");
      printf("#### GMRES inner loop: i=%i\n",i);
      printf("####   s[0]: %12.8f\n",s[0]);
      lresult->number_of_linear_iterations++;
      /* Matrix-vector mutliply: A*v[i] */
      if (s_dmatmul_set(theMG,bl,level,np->s,A,np->v[i],EVERY_CLASS))
        NP_RETURN(1,lresult->error_code);
      /* preconditioner solve: Mw = A*v[i] */
      if (a_dset(theMG,bl,level,np->w,EVERY_CLASS,0.0) != NUM_OK)
        NP_RETURN(1,result[0]);
      if ((*np->Iter->Iter)(np->Iter,level,np->w,np->s,A,
                            &lresult->error_code))
        REP_ERR_RETURN (1);

      /* form a column of the upper Hessenberg matrix / Arnoldi
             process */
      for (k=0; k<=i; k++) {
        /* form an entry of the upper Hessenberg matrix */
                            #ifdef ModelP
        if (s_dcopy(theMG,bl,level,np->s,np->w)!= NUM_OK)
          NP_RETURN(1,lresult->error_code);
        if (a_vector_makeinconsistent(theMG,bl,level,np->s) != NUM_OK)
          NP_RETURN(1,lresult->error_code);
        if (s_ddot(theMG,bl,level,np->s,np->v[k],rnorm) !=NUM_OK)
          NP_RETURN(1,result[0]);
                                #else
        if (s_ddot(theMG,bl,level,np->w,np->v[k],rnorm) !=NUM_OK)
          NP_RETURN(1,result[0]);
                                #endif
        lambda = 0.0;
        for (j=0; j<ncomp; j++) lambda += rnorm[j];
        H[k][i] = lambda;
        printf("#### H[k,i]: %8.4f \n",H[k][i]);

        /* update the vector w  = w-h[k,i]*v */
        for (j=0; j<ncomp; j++) scal[j] = -lambda;
        if (a_daxpy(theMG,bl,level,np->w,EVERY_CLASS,scal,np->v[k])
            != NUM_OK)
          NP_RETURN(1,result[0]);
      }                   /* k */

      /* form H[i+1][i] = norm(w) */
                        #ifdef ModelP
      if (s_dcopy(theMG,bl,level,np->s,np->w)!= NUM_OK)
        NP_RETURN(1,lresult->error_code);
      if (a_vector_makeinconsistent(theMG,bl,level,np->s) != NUM_OK)
        NP_RETURN(1,lresult->error_code);
      if (s_ddot(theMG,bl,level,np->s,np->w,rnorm) !=NUM_OK)
        NP_RETURN(1,result[0]);
                        #else
      if (s_ddot(theMG,bl,level,np->w,np->w,rnorm) !=NUM_OK)
        NP_RETURN(1,result[0]);
                        #endif
      lambda = 0.0;
      for (j=0; j<ncomp; j++) lambda += rnorm[j];
      lambda = sqrt(lambda);
      H[i+1][i] = lambda;
      printf("####norm(w): %8.4f \n",H[i+1][i]);

      /* set v[i+1] = w/H[i+1][i] */ /* #### check scaling #### */
      if (s_dcopy(theNP->base.mg,bl,level,np->v[i+1],np->w)!= NUM_OK)
        NP_RETURN(1,lresult->error_code);
      for (j=0; j<ncomp; j++) rnorm[j] = 1.0/lambda;
      if (a_dscale(theMG,bl,level,np->v[i+1],EVERY_CLASS,rnorm)!=NUM_OK)
        NP_RETURN(1,result[0]);

      /* apply Givens rotations */
      for (k=0; k<i; k++) {
        lambda   =  cs[k]*H[k][i] + sn[k]*H[k+1][i];
        H[k+1][i] = -sn[k]*H[k][i] + cs[k]*H[k+1][i];
        H[k][i]   = lambda;
      }                   /* k */

      /* construct next Givens rotation in a numerically stable manner */
      if (H[i+1][i] == 0.0) {
        cs[i] = 1.0;
        sn[i] = 0.0;
      } else if (ABS(H[i+1][i]) > ABS(H[i][i])) {
        lambda = H[i][i] / H[i+1][i];
        sn[i] = 1.0 / sqrt(1.0 + lambda*lambda);
        cs[i] = lambda*sn[i];
      } else {
        lambda = H[i+1][i] / H[i][i];
        cs[i] = 1.0 / sqrt(1.0 + lambda*lambda);
        sn[i] = lambda*cs[i];
      }

      /* form the (recursively computed) residual norm */
      lambda   = cs[i]*s[i];
      s[i+1]   = -sn[i]*s[i];
      s[i]     = lambda;
      H[i][i]   = cs[i]*H[i][i] + sn[i]*H[i+1][i];
      H[i+1][i] = 0.0;
      printf("####      i: %12i\n",i);
      printf("#### MAXRES: %12i\n",MAX_RESTART);
      printf("####   s[i]: %12.8f\n",s[i]);
      printf("#### s[i+1]: %12.8f \n",s[i+1]);
      printf("####  cs[i]: %12.8f \n",cs[i]);
      printf("####  sn[i]: %12.8f\n",sn[i]);
      printf("#### H[i,i]: %12.8f\n",H[i][i]);
      printf("####    TOL: %12.8f\n", tol);
      /* if the error is sufficiently small,
         perform the update and exit */
      printf("#### abs(s[i+1]): %12.8f\n",fabs(s[i+1]));
      printf("####         TOL: %12.8f\n",tol);
      if (fabs(s[i+1]) < tol) break;
    }             /* end: GMRES inner loop:i */

    /* ensure i is the proper value */
    if (i>=np->restart) i=np->restart-1;
    printf("#####################################\n");
    /* solve the upper Hessenberg system */
    for (i1=0; i1<=i; i1++)
      for (i2=0; i2<=i; i2++)
        Hsq[i1*(i+1)+i2] = H[i1][i2];
    for (i1=0; i1<=i; i1++) {
      printf("% i  @@@@ ",i1);
      for (i2=0; i2<=i; i2++)
        printf("%12.8f ",Hsq[i1*(i+1)+i2]);
      printf("\n");
    }

    if (SolveFullMatrix(i+1,y,Hsq,s)) {
      UserWriteF("GMRESSolver: decompostion failed");
      NP_RETURN(1,lresult->error_code);
    }
    printf("y ");
    for (i2=0; i2<=i; i2++)
      printf("%12.8f ",y[i2]);
    printf("\n");
    printf("s ");
    for (i2=0; i2<=i; i2++)
      printf("%12.8f ",s[i2]);
    printf("\n");
    /* Perform a full matrix-vector multiply
       (NOT THE COEFFICIENT MATRIX!) */
    /*x = x + V(:,1:i)*y; */

    if (s_dset(theMG,bl,level,np->c,0.0)!= NUM_OK)
      NP_RETURN(1,lresult->error_code);

    for (i1=0; i1<=i; i1++) {
      for (j=0; j<ncomp; j++) scal[j] = y[i1];
      if (a_daxpy(theMG,bl,level,np->c,EVERY_CLASS,scal,np->v[i1])
          != NUM_OK)
        NP_RETURN(1,lresult->error_code);
    }
    if (a_daxpy(theMG,bl,level,x,EVERY_CLASS,Factor_One,np->c) != NUM_OK)
      NP_RETURN(1,lresult->error_code);
    if (s_dmatmul_minus(theMG,bl,level,b,A,np->c,EVERY_CLASS) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (LinearResiduum(theNP,bl,level,x,b,A,lresult))
      NP_RETURN(1,lresult->error_code);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE))
        NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) ||
        sc_cmp(lresult->last_defect,defect2reach,b))
      lresult->converged = 1;

  }       /* end: GMRES outer loop:it*/
  CSTOP(ti,ii);

  if (np->display > PCR_NO_DISPLAY) {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE))
      NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg"))
      NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1)))
      NP_RETURN(1,lresult->error_code);
    UserWriteF("LS  : L=%2d N=%2d TSOLVE=%10.4lg TIT=%10.4lg\n",level,
               lresult->number_of_linear_iterations,ti,
               ti/lresult->number_of_linear_iterations);
  }

  return (0);
}

static INT GMRESConstruct (NP_BASE *theNP)
{
  NP_GMRES *np;

  theNP->Init = GMRESInit;
  theNP->Display = GMRESDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_GMRES *) theNP;
  np->ls.PreProcess = GMRESPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = GMRESSolver;
  np->ls.PostProcess = GMRESPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   sqcg - numproc for the squared cg method

   DESCRIPTION:
   This numproc executes the squared conjugate gradient method.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
       [$p <con>] [$h <tmp>]
           [$R <restart>];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus
   .  $p~<con> - conjugate vector
   .  $h~<tmp> - temporaty vector
   .  $R~<restart> - restart index

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT SQCGInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SQCG *np;

  np = (NP_SQCG *) theNP;
  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  np->p = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->h = ReadArgvVecDesc(theNP->mg,"h",argc,argv);
  np->d = ReadArgvVecDesc(theNP->mg,"d",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("R",&(np->restart),argc,argv))
    np->restart = 0;
  if (np->restart<0) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->baselevel = 0;

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT SQCGDisplay (NP_BASE *theNP)
{
  NP_SQCG *np;

  np = (NP_SQCG *) theNP;
  NPLinearSolverDisplay(&np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"r",(int)np->restart);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->display == PCR_NO_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->r != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->p != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  if (np->h != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"h",ENVITEM_NAME(np->h));
  if (np->d != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"d",ENVITEM_NAME(np->d));

  return (0);
}

static INT SQCGPreProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_SQCG *np;

  np = (NP_SQCG *) theNP;

  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->r)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->h)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->ls.base.mg,np->baselevel,level,x,&np->d)) NP_RETURN(1,result[0]);

  return(0);
}

static INT SQCGPostProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SQCG *np;

  np = (NP_SQCG *) theNP;
  FreeVD(np->ls.base.mg,np->baselevel,level,np->r);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->p);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->h);
  FreeVD(np->ls.base.mg,np->baselevel,level,np->d);

  return(0);
}

static INT SQCGSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_SQCG *np;
  VEC_SCALAR defect2reach,scal;
  INT i,j,PrintID,restart;
  char text[DISPLAY_WIDTH+4];
  DOUBLE s,t,lambda;

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  /* prepare */
  np = (NP_SQCG *) theNP;

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY) if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY) if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;

  /* go */
  restart = 1;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;

    /* restart ? */
    if ((np->restart>0 && i%np->restart==0) || restart)
    {
      if (s_dtpmatmul_set(theNP->base.mg,np->baselevel,level,np->r,A,b,EVERY_CLASS)) REP_ERR_RETURN (1);
      if (s_dcopy(theNP->base.mg,np->baselevel,level,np->p,np->r)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      restart = 0;
    }

    /* update x, b */
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->p,np->r,Factor_One,&s)!=NUM_OK) REP_ERR_RETURN (1);
    if (s_dmatmul_set(theNP->base.mg,np->baselevel,level,np->h,A,np->p,EVERY_CLASS)) REP_ERR_RETURN (1);
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->h,np->h,Factor_One,&t)!=NUM_OK) REP_ERR_RETURN (1);
    lambda = s/t;
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=lambda;
    if (s_daxpy (theNP->base.mg,np->baselevel,level,x,scal,np->p)!= NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-lambda;
    if (s_daxpy (theNP->base.mg,np->baselevel,level,b,scal,np->h)!= NUM_OK) REP_ERR_RETURN (1);
    if (s_dtpmatmul_set(theNP->base.mg,np->baselevel,level,np->d,A,np->h,EVERY_CLASS)) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->r,scal,np->d)!= NUM_OK) REP_ERR_RETURN (1);
    if (s_ddot_sv (theNP->base.mg,np->baselevel,level,np->d,np->r,Factor_One,&s)!=NUM_OK) REP_ERR_RETURN (1);
    lambda = s/t;
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-lambda;
    if (s_dscale (theNP->base.mg,np->baselevel,level,np->p,scal)) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,np->p,Factor_One,np->r)!= NUM_OK) REP_ERR_RETURN (1);

    /* redisuum */
    if (LinearResiduum(theNP,np->baselevel,level,x,b,A,lresult)) REP_ERR_RETURN (1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
  }

  return (0);
}

static INT SQCGConstruct (NP_BASE *theNP)
{
  NP_SQCG *np;

  theNP->Init = SQCGInit;
  theNP->Display = SQCGDisplay;
  theNP->Execute = NPLinearSolverExecute;

  np = (NP_SQCG *) theNP;
  np->ls.PreProcess = SQCGPreProcess;
  np->ls.Defect = LinearDefect;
  np->ls.Residuum = LinearResiduum;
  np->ls.Solver = SQCGSolver;
  np->ls.PostProcess = SQCGPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   ldcs - numproc for defect correction linear solvers

   DESCRIPTION:
   This numproc executes a defect correction scheme for a linear solver.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $LS <linear solver [$d {full|red|no}]
       $DC <mat sym>;
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $LS~<linear~solver> - linear solver numproc
   .  $d~{full|red|no} - display modus
   .  $DC~<mat~sym> - defect correction matrix

   'npexecute <name> [$i] [$d] [$r] [$s] [$p];'

   .  $i - preprocess
   .  $d - replace right hand side by the defect
   .  $r - compute the residuum of the defect
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
   npcreate damp		$c smd;
   npcreate cdeq		$c cd;
   npcreate transfer    $c transfer;
   npcreate lu			$c lu;
   npcreate base		$c ls;
   npcreate sm			$c scgs;
   npcreate lmgc                $c lmgc;
   npcreate lin                 $c ls;
   npcreate dc			$c ldcs;

   npinit transfer;
   npinit damp $rp 0 $off 0;
   npinit lu;
   npinit base $red 1e-8 $m 50 $I lu $display no $abslimit 0;
   npinit sm $n 1 $mode ff $limit 0 $bl 0 $gamma 1 $display no;
   npinit lmgc $S sm sm base $T transfer $n1 1 $n2 1 $g 1 $b @:BL;
   npinit lin $red 1e-8 $m 100 $I lmgc $display no $absimit 0;
   npinit dc $A mat $DC dcmat $x x $b b $red 0 $m 20 $LS lin $display full;
   npinit cdeq $A mat $b b $x x $v @:VELO $p 1 $d 0 $m fu $a 1;
   npex cdeq;
   orderv $m CCFFLL $c @:CUT $a;
   npex dc $i $r $d $s $p;
   .ve
   D*/
/****************************************************************************/

static INT LDCSInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LDCS *np;

  np = (NP_LDCS *) theNP;

  if (ReadArgvINT("m",&(np->maxiter),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->linsol = (NP_LINEAR_SOLVER *)       ReadArgvNumProc(theNP->mg,"LS",LINEAR_SOLVER_CLASS_NAME,argc,argv);
  if (np->linsol == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->DC = ReadArgvMatDesc(theNP->mg,"DC",argc,argv);
  if (np->DC == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);

  return (NPLinearSolverInit(&np->ls,argc,argv));
}

static INT LDCSDisplay (NP_BASE *theNP)
{
  NP_LDCS *np;

  np = (NP_LDCS *) theNP;
  NPLinearSolverDisplay(&np->ls);

  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->linsol != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"LS",ENVITEM_NAME(np->linsol));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"LS","---");
  if (np->display == PCR_NO_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->DC != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"DC",ENVITEM_NAME(np->DC));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"DC","---");
  if (np->b != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->b));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"b","---");
  if (np->c != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->c));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"c","---");

  return (0);
}

static INT LDCSPreProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_LDCS *np;

  /* store passed XXXDATA_DESCs */
  NPLS_A(theNP) = A;
  NPLS_x(theNP) = x;
  NPLS_b(theNP) = b;

  np = (NP_LDCS *) theNP;
  if (np->linsol->PreProcess != NULL)
    if ((*np->linsol->PreProcess)(np->linsol,level,x,b,np->DC,baselevel,result))
      REP_ERR_RETURN(1);
  np->baselevel = MIN(*baselevel,level);

  return(0);
}

static INT LDCSSolver (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, VEC_SCALAR abslimit, VEC_SCALAR reduction, LRESULT *lresult)
{
  NP_LDCS *np;
  VEC_SCALAR defect2reach;
  INT i,bl,PrintID;
  char text[DISPLAY_WIDTH+4];
  NP_LINEAR_SOLVER * linsol;
  LRESULT linsolresult;

  /* store passed reduction and abslimit */
  for (i=0; i<VD_NCOMP(x); i++)
  {
    NPLS_red(theNP)[i] = reduction[i];
    NPLS_abs(theNP)[i] = abslimit[i];
  }

  np = (NP_LDCS *) theNP;
  linsol = np->linsol;
  if (linsol == NULL) NP_RETURN(1,lresult->error_code);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY)
    if (PreparePCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY)
    if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;
  if (AllocVDFromVD(theNP->base.mg,np->baselevel,level,x,&np->b)) REP_ERR_RETURN(1);
  if (AllocVDFromVD(theNP->base.mg,np->baselevel,level,x,&np->c)) REP_ERR_RETURN(1);
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;
    if (s_dcopy(theNP->base.mg,np->baselevel,level,np->b,b)) REP_ERR_RETURN(1);
    if (s_dset(theNP->base.mg,np->baselevel,level,np->c,0.0)) REP_ERR_RETURN(1);
    if ((*linsol->Residuum)(linsol,np->baselevel,level,np->c,np->b,np->DC,&linsolresult)) NP_RETURN(1,lresult->error_code);
    if ((*linsol->Solver)(linsol,level,np->c,np->b,np->DC,linsol->abslimit,linsol->reduction,&linsolresult)) REP_ERR_RETURN (1);
    if (s_daxpy (theNP->base.mg,np->baselevel,level,x,Factor_One,np->c)) REP_ERR_RETURN(1);
    if (s_dmatmul_minus(theNP->base.mg,np->baselevel,level,b,A,np->c,EVERY_CLASS)) REP_ERR_RETURN(1);
    if (LinearResiduum(theNP,bl,level,x,b,A,lresult)) REP_ERR_RETURN(1);
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (sc_cmp(lresult->last_defect,abslimit,b) || sc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  FreeVD(theNP->base.mg,np->baselevel,level,np->b);
  FreeVD(theNP->base.mg,np->baselevel,level,np->c);
  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
  }

  return (0);
}

static INT LDCSPostProcess (NP_LINEAR_SOLVER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_LDCS *np;

  np = (NP_LDCS *) theNP;

  if (np->linsol->PostProcess == NULL) return(0);
  return((*np->linsol->PostProcess)(np->linsol,level,x,b,np->DC,result));
}

static INT LDCSConstruct (NP_BASE *theNP)
{
  NP_LDCS *np;

  theNP->Init             = LDCSInit;
  theNP->Display          = LDCSDisplay;
  theNP->Execute          = NPLinearSolverExecute;

  np = (NP_LDCS *) theNP;
  np->ls.PreProcess       = LDCSPreProcess;
  np->ls.Defect           = LinearDefect;
  np->ls.Residuum         = LinearResiduum;
  np->ls.Solver           = LDCSSolver;
  np->ls.PostProcess      = LDCSPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitLinearSolver	- Init this file

   SYNOPSIS:
   INT InitLinearSolver ();

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

INT InitLinearSolver ()
{
  INT i;

  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".ls",sizeof(NP_LS),LSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".cg",sizeof(NP_CG),CGConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".cr",sizeof(NP_CR),CRConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".bcg",sizeof(NP_BCG),BCGConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".bcgs",sizeof(NP_BCGS),BCGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".gmres",sizeof(NP_GMRES),GMRESConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".sqcg",sizeof(NP_SQCG),SQCGConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".ldcs",sizeof(NP_LDCS),LDCSConstruct))
    REP_ERR_RETURN (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;
  if (MakeStruct(":ls")) REP_ERR_RETURN(__LINE__);
  if (MakeStruct(":ls:avg")) REP_ERR_RETURN(__LINE__);

  return (0);
}
