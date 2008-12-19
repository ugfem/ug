// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  els.c	                                                                                                        */
/*																			*/
/* Purpose:   extended linear solver num procs                                  */
/*																			*/
/* Author:    Klaus Johannsen                                               */
/*            IWR/Technische Simulation                                     */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: klaus.johannsen@iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   19.07.02 begin                                                */
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

#include "config.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "ugdevices.h"
#include "general.h"
#include "debug.h"
#include "ugstruct.h"
#include "gm.h"
#include "scan.h"
#include "block.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"
#include "transgrid.h"

#include "eiter.h"
#include "els.h"

USING_UG_NAMESPACES

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
#define MAX_RESTART 30

#define CSTART()    clock_start=CURRENT_TIME_LONG;
#define CSTOP(t,c)  t+=(CURRENT_TIME_LONG-clock_start);c++

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct np_els
{
  NP_ELINEAR_SOLVER ls;

  NP_EITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;

  EVECDATA_DESC *c;

  INT (*Prepare)
    (struct np_els *,                        /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*Update)
    (struct np_els *,                        /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector                 */
    EVECDATA_DESC *,                             /* correction vector               */
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    INT *);                                      /* result                          */
  INT (*Close)
    (struct np_els *,                        /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_els NP_ELS;

typedef struct
{
  NP_ELINEAR_SOLVER ls;

  NP_EITER *Iter;

  INT maxiter;
  INT baselevel;
  INT display;
  INT restart;

  DOUBLE rho, omega;
  EVEC_SCALAR weight;
  EVEC_SCALAR old_defect;
  EMATDATA_DESC *B;
  EVECDATA_DESC *r;
  EVECDATA_DESC *p;
  EVECDATA_DESC *v;
  EVECDATA_DESC *s;
  EVECDATA_DESC *t;
  EVECDATA_DESC *q;

} NP_EBCGS;

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

static DOUBLE clock_start;
static DOUBLE basetime;

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
   NP_ELINEAR_SOLVER - type definition for extended linear solvers

   DESCRIPTION:
   This numproc type is used for the description of extended linear solvers.
   It can be called by the given interface from a nonlinearsolver.
   Initializing the data is optional; it can be done with

   'INT NPELinearSolverInit (NP_ELINEAR_SOLVER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPELinearSolverDisplay (NP_ELINEAR_SOLVER *theNP);'
   'INT NPELinearSolverExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   typedef struct {
        INT error_code;                     // error code
        INT converged;                      // error code
        EVEC_SCALAR first_defect;           // first defect
        EVEC_SCALAR last_defect;            // last defect
        INT number_of_linear_iterations;    // number of iterations
   } ELRESULT;

   struct np_elinear_solver {
        NP_BASE base;                        // inherits base class

        // data (optional, necessary for calling the generic execute routine)
    EVECDATA_DESC *x;                    // solution
    EVECDATA_DESC *b;                    // defect
    EMATDATA_DESC *A;                    // matrix
        EVEC_SCALAR reduction;               // reduction factor
        EVEC_SCALAR abslimit;                // absolute limit for the defect

        // functions
        INT (*PreProcess)
             (struct np_elinear_solver *,    // pointer to (derived) object
                  INT,                           // level
                  EVECDATA_DESC *,               // solution vector
                  EVECDATA_DESC *,               // defect vector
                  EMATDATA_DESC *,               // matrix
                  INT *,                         // baselevel used by the solver
                  INT *);                        // result
    INT (*Defect)                        // b := b - Ax
             (struct np_elinear_solver *,    // pointer to (derived) object
                  INT,                           // level
                  EVECDATA_DESC *,               // solution vector
                  EVECDATA_DESC *,               // defect vector
                  EMATDATA_DESC *,               // matrix
                  INT *);                        // result
    INT (*Residuum)                      // computes norm of the defect
             (struct np_elinear_solver *,    // pointer to (derived) object
                  INT,                           // from level
                  INT,                           // to level
                  EVECDATA_DESC *,               // solution vector
                  EVECDATA_DESC *,               // defect vector
                  EMATDATA_DESC *,               // matrix
                  ELRESULT *);                   // result structure
    INT (*Solver)                        // b := b - Ax
             (struct np_elinear_solver *,    // pointer to (derived) object
                  INT,                           // level
                  EVECDATA_DESC *,               // solution vector
                  EVECDATA_DESC *,               // defect vector
                  EMATDATA_DESC *,               // matrix
                  EVEC_SCALAR,                   // reduction factor
                  EVEC_SCALAR,                   // absolute limit for the defect
                  ELRESULT *);                   // result structure
        INT (*PostProcess)
             (struct np_elinear_solver *,    // pointer to (derived) object
                  INT,                           // level
                  EVECDATA_DESC *,               // solution vector
                  EVECDATA_DESC *,               // defect vector
                  EMATDATA_DESC *,               // matrix
                  INT *);                        // result
   };
   typedef struct np_elinear_solver NP_ELINEAR_SOLVER;
   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NS_DIM_PREFIX NPELinearSolverInit (NP_ELINEAR_SOLVER *np, INT argc , char **argv)
{
  INT i;

  np->A = ReadArgvEMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvEVecDesc(np->base.mg,"x",argc,argv);
  np->b = ReadArgvEVecDesc(np->base.mg,"b",argc,argv);
  if (esc_read(np->abslimit,NP_FMT(np),np->x,"abslimit",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->abslimit[i] = ABS_LIMIT;

  if (ReadArgvINT("setbasetime",&NPELS_setbasetime(np),argc,argv))
  {
    NPELS_setbasetime(np) = 0;
  }

  if (ReadArgvINT("printbasetime",&NPELS_printbasetime(np),argc,argv))
  {
    NPELS_printbasetime(np) = 0;
  }

  if (esc_read(np->reduction,NP_FMT(np),np->x,"red",argc,argv)) return(NP_ACTIVE);
  if (esc_read(np->abslimit,NP_FMT(np),np->x,"abslimit",argc,argv)) return(NP_ACTIVE);
  if ((np->x == NULL) || (np->b == NULL) || (np->A == NULL)) return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NS_DIM_PREFIX NPELinearSolverDisplay (NP_ELINEAR_SOLVER *np)
{
  if ((np->x != NULL) || (np->b != NULL) || (np->A != NULL)) {
    UserWrite("symbolic user data:\n");
    if (np->A != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
    if (np->x != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
    if (np->b != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"b",ENVITEM_NAME(np->b));
    UserWrite("\n");
  }
  UserWrite("configuration parameters:\n");
  if (np->x != NULL)
  {
    if (esc_disp(np->reduction,np->x,"red")) REP_ERR_RETURN (1);
    if (esc_disp(np->abslimit,np->x,"abslimit")) REP_ERR_RETURN (1);
  }
  UserWriteF(DISPLAY_NP_FORMAT_SI,"setbasetime",(int)NPELS_setbasetime(np));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"printbasetime",(int)NPELS_printbasetime(np));

  return(0);
}

INT NS_DIM_PREFIX NPELinearSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ELINEAR_SOLVER *np;
  ELRESULT lresult;
  INT result,level,bl;

  np = (NP_ELINEAR_SOLVER *) theNP;

#ifdef USE_FAMG
  if (ReadArgvOption("t",argc,argv))
    level = -1;
  else
    level = CURRENTLEVEL(theNP->mg);
#else
  level = CURRENTLEVEL(theNP->mg);
#endif

  result = 0;
  bl = 0;

  if (np->x == NULL)
  {
    PrintErrorMessage('E',"NPELinearSolverExecute","no vector x");
    REP_ERR_RETURN (1);
  }
  if (np->b == NULL)
  {
    PrintErrorMessage('E',"NPELinearSolverExecute","no vector b");
    REP_ERR_RETURN (1);
  }
  if (np->A == NULL)
  {
    PrintErrorMessage('E',"NPELinearSolverExecute","no matrix A");
    REP_ERR_RETURN (1);
  }

  if (ReadArgvOption("i",argc,argv))
  {
    if (np->PreProcess == NULL)
    {
      PrintErrorMessage('E',"NPELinearSolverExecute","no PreProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PreProcess)(np,level,np->x,np->b,np->A,&bl,&result))
    {
      UserWriteF("NPELinearSolverExecute: PreProcess failed, error code %d\n", result);
            #ifndef ModelP
      REP_ERR_RETURN (1);
            #endif
      REP_ERR_INC;
                #ifndef Debug
      REP_ERR_RETURN (1);
            #endif
    }
  }

  if (ReadArgvOption("d",argc,argv))
  {
    if (np->Defect == NULL)
    {
      PrintErrorMessage('E',"NPELinearSolverExecute","no Defect");
      REP_ERR_RETURN (1);
    }
    if ((*np->Defect)(np,level,np->x,np->b,np->A,&result))
      UserWriteF("NPELinearSolverExecute: Defect failed, error code %d\n", result);
  }

  if (ReadArgvOption("r",argc,argv))
  {
    if (np->Residuum == NULL)
    {
      PrintErrorMessage('E',"NPELinearSolverExecute","no Residuum");
      REP_ERR_RETURN (1);
    }
    if ((*np->Residuum)(np,bl,level,np->x,np->b,np->A,&lresult)) {
      UserWriteF("NPELinearSolverExecute: Residuum failed, error code %d\n", result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("s",argc,argv))
  {
    if (np->Solver == NULL)
    {
      PrintErrorMessage('E',"NPELinearSolverExecute","no Solver");
      REP_ERR_RETURN (1);
    }
    if ((*np->Solver)(np,level,np->x,np->b,np->A,np->abslimit,np->reduction,&lresult))
    {
      UserWriteF("NPELinearSolverExecute: Solver failed, error code %d\n",lresult.error_code);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("p",argc,argv))
  {
    if (np->PostProcess == NULL)
    {
      PrintErrorMessage('E',"NPELinearSolverExecute","no PostProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PostProcess)(np,level,np->x,np->b,np->A,&result))
    {
      UserWriteF("NPELinearSolverExecute: PostProcess failed, error code %d\n",result);
      REP_ERR_RETURN (1);
    }
  }
  return(0);
}

/* tools for linear solvers */

static INT ELinearSolverPreProcess (NP_ELINEAR_SOLVER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_ELS *np;

  /* store passed XXXDATA_DESCs */
  NPELS_A(theNP) = A;
  NPELS_x(theNP) = x;
  NPELS_b(theNP) = b;

  np = (NP_ELS *) theNP;
  if (np->Iter==NULL) REP_ERR_RETURN(1);
  if (np->Iter->PreProcess != NULL)
    if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);
  np->baselevel = MIN(*baselevel,level);

  return(0);
}

static INT ELinearDefect (NP_ELINEAR_SOLVER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  NP_ELS *np;
  INT bl;

  np = (NP_ELS *) theNP;
  bl = MIN(FULLREFINELEVEL(NP_MG(theNP)),MAX(0,np->baselevel));
  if (dematmul_minus(NP_MG(theNP),bl,level,ON_SURFACE,b,A,x)) NP_RETURN(1,result[0]);

  return (*result);
}

static INT ELinearResiduum (NP_ELINEAR_SOLVER *theNP, INT bl, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, ELRESULT *lresult)
{
  NP_ELS *np;

  np = (NP_ELS *) theNP;

        #ifdef ModelP
  if (a_vector_collect(NP_MG(theNP),bl,level,b->vd)) NP_RETURN(1,lresult->error_code);
        #endif
  if (denrm2x(NP_MG(theNP),bl,level,ON_SURFACE,b,lresult->last_defect)) NP_RETURN(1,lresult->error_code);

  return(0);
}

static INT EnergyResiduum (NP_ELINEAR_SOLVER *theNP, INT bl, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, ELRESULT *lresult)
{
  EVEC_SCALAR err;
  EVECDATA_DESC *temp;
  INT i;

        #ifdef ModelP
  if (a_vector_collect(theNP->base.mg,bl,level,b->vd)) NP_RETURN(1,lresult->error_code);
        #endif

  temp = NULL;
  if(AllocEVDFromEVD(theNP->base.mg,bl,level,x,&temp)) NP_RETURN(1,lresult->error_code);
  if (deset(theNP->base.mg, bl,level,ON_SURFACE,temp, 0.0)) NP_RETURN(1,lresult->error_code);
  if(dematmul(theNP->base.mg,bl,level,ON_SURFACE,temp,A,x)) NP_RETURN(1,lresult->error_code);
  if (dedotx(theNP->base.mg,bl, level,ON_SURFACE, temp,x,err)) NP_RETURN(1,lresult->error_code);

  /*ABS is necessary, because residuum is computed for the solution.
     In this case, the value at the dirichlet nodes is not zero und thus
     A = MAT does not represent a(.,.) */

  lresult->last_defect[0] = 0.0;
  for (i=0; i<VD_NCOMP(x->vd)+x->n; i++) lresult->last_defect[0] += err[i];
  lresult->last_defect[0] = sqrt((DOUBLE) ABS(lresult->last_defect[0]));
  lresult->last_defect[0] /= sqrt((DOUBLE) (VD_NCOMP(x->vd)+x->n));
  for(i=1; i<VD_NCOMP(x->vd)+x->n; i++)
    lresult->last_defect[i] =  lresult->last_defect[0];

  FreeEVD(theNP->base.mg,bl,level,temp);

  return(0);
}

static INT ELinearSolver (NP_ELINEAR_SOLVER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, EVEC_SCALAR abslimit, EVEC_SCALAR reduction, ELRESULT *lresult)
{
  NP_ELS *np;
  EVEC_SCALAR defect2reach;
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
  for (i=0; i<VD_NCOMP(x->vd)+x->n; i++)
  {
    NPELS_red(theNP)[i] = reduction[i];
    NPELS_abs(theNP)[i] = abslimit[i];
  }

  np = (NP_ELS *) theNP;
  bl = np->baselevel;
  if (np->Iter->Iter == NULL) NP_RETURN(1,lresult->error_code);
  if (np->Update == NULL) NP_RETURN(1,lresult->error_code);
  if (AllocEVDFromEVD(NP_MG(theNP),bl,level,x,&np->c)) NP_RETURN(1,lresult->error_code);
  if (np->Prepare != NULL)
    if ((*np->Prepare)(np,level,x,&lresult->error_code)) REP_ERR_RETURN (1);

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display > PCR_NO_DISPLAY) if (PrepareEPCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  if(NPELS_printbasetime(theNP)) basetime = 0.0;
  CSTART(); ti=0; ii=0;
  for (i=0; i<VD_NCOMP(x->vd)+x->n; i++) lresult->first_defect[i] = lresult->last_defect[i];
  if (esc_mul(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display > PCR_NO_DISPLAY)
    if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (esc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;
  lresult->number_of_linear_iterations = 0;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;
    if (deset(NP_MG(theNP),level,level,ALL_VECTORS,np->c,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
    if ((*np->Iter->Iter)(np->Iter,level,np->c,b,A,&lresult->error_code)) REP_ERR_RETURN (1);
    if ((*np->Update)(np,level,x,np->c,b,A,&lresult->error_code)) REP_ERR_RETURN (1);
    if (np->ls.Residuum == EnergyResiduum)
    {
      if((np->ls.Residuum)(theNP,bl,level,np->c,b,A,lresult)) REP_ERR_RETURN(1);
    }
    else
    {
      if (ELinearResiduum(theNP,bl,level,x,b,A,lresult)) REP_ERR_RETURN(1);
    }
    if (np->display > PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    if (esc_cmp(lresult->last_defect,abslimit,b) || esc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      lresult->number_of_linear_iterations=i+1;
      break;
    }
  }
  if (!lresult->converged) lresult->number_of_linear_iterations=i;
  if (FreeEVD(NP_MG(theNP),bl,level,np->c)) REP_ERR_RETURN(1);
  if (np->Close != NULL) if ((*np->Close)(np,level,&lresult->error_code)) REP_ERR_RETURN (1);
  CSTOP(ti,ii);
  if(NPELS_setbasetime(theNP)) basetime += ti;

  if (np->display > PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (lresult->number_of_linear_iterations != 0)
      if(NPELS_printbasetime(theNP))
        UserWriteF("LS  : L=%2d N=%2d TSOLVE=%10.4g TIT=%10.4g TBASE=%g\n",level,lresult->number_of_linear_iterations,ti,ti/lresult->number_of_linear_iterations,basetime);
      else
        UserWriteF("LS  : L=%2d N=%2d TSOLVE=%10.4g TIT=%10.4g\n",level,lresult->number_of_linear_iterations,ti,ti/lresult->number_of_linear_iterations);
    else
      UserWriteF("LS  : L=%2d N=%2d TSOLVE=%10.4g\n",level,lresult->number_of_linear_iterations,ti);
  }

  return (0);
}

static INT ELinearSolverPostProcess (NP_ELINEAR_SOLVER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  NP_ELS *np;

  np = (NP_ELS *) theNP;

  if (np->Iter != NULL)
    if (np->Iter->PostProcess != NULL)
      if ((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result)) NP_RETURN(1,result[0]);

  np->baselevel = MAX(BOTTOMLEVEL(theNP->base.mg),np->baselevel);

  return(0);
}

static INT ELinearSolverInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ELS *np;

  np = (NP_ELS *) theNP;

  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    REP_ERR_RETURN(NP_NOT_ACTIVE);

  np->display = ReadArgvDisplay(argc,argv);
  np->Iter = (NP_EITER *) ReadArgvNumProc(theNP->mg,"I",EITER_CLASS_NAME,argc,argv);
  if (np->Iter == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->baselevel = 0;
  np->c = ReadArgvEVecDesc(theNP->mg,"c",argc,argv);
  if (ReadArgvOption("E",argc,argv)) np->ls.Residuum = EnergyResiduum;

  return (NPELinearSolverInit(&np->ls,argc,argv));
}

static INT ELinearSolverDisplay (NP_BASE *theNP)
{
  NP_ELS *np;

  np = (NP_ELS *) theNP;
  NPELinearSolverDisplay(&np->ls);

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

/****************************************************************************/
/*D
   els - numproc for extended linear solvers

   DESCRIPTION:
   This numproc executes a extended linear solver: it performs an iteration
   up to convergence.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <sc double list>] [$abslimit <sc double list>]
       $m <maxit> $I <ext.iteration> [$d {full|red|no}];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<sc~double~list> - reduction factor
   .  $abslimit~<sc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $I~<ext.iteration> - extended iteration numproc
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
   npcreate basesolver $c els;    npinit basesolver $red 0.001 $I base;
   npcreate transfer $c transfer; npinit transfer;
   npcreate lmgc $c elmgc;        npinit lmgc $S pre post basesolver $T transfer;
   npcreate mgs $c els;           npinit mgs $A MAT $x sol $b rhs
                                          $red 0.00001 $I lmgc $d full;
   npexecute mgs $i $d $r $s $p;
   .ve
   D*/
/****************************************************************************/

static INT ELSUpdate (NP_ELS *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *c, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  if (deadd(theNP->ls.base.mg,theNP->baselevel,level,ALL_VECTORS,x,c) != NUM_OK) NP_RETURN(1,result[0]);

  return(0);
}

static INT ELSConstruct (NP_BASE *theNP)
{
  NP_ELS *np;

  theNP->Init = ELinearSolverInit;
  theNP->Display = ELinearSolverDisplay;
  theNP->Execute = NPELinearSolverExecute;

  np = (NP_ELS *) theNP;
  np->ls.PreProcess = ELinearSolverPreProcess;
  np->ls.Defect = ELinearDefect;
  np->ls.Residuum = ELinearResiduum;
  np->ls.Solver = ELinearSolver;
  np->ls.PostProcess = ELinearSolverPostProcess;

  np->Prepare = NULL;
  np->Update = ELSUpdate;
  np->Close = NULL;

  return(0);
}

/****************************************************************************/
/*D
   ebcgs - numproc for the bi cg stab method for extended linear systems

   DESCRIPTION:
   This numproc executes the bi-conjugate gradient method for extended systems.

   .vb
   npinit <name> [$x <sol>] [$b <rhs>] [$A <mat sym>]
       [$red <esc double list>] [$abslimit <esc double list>]
       $m <maxit> $I <iteration> [$d {full|red|no}]
       [$p <con>] [$t <tmp>]
       [$R <restart>] [$w <esc double list>];
   .ve

   .  $x~<sol> - solution vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $red~<esc~double~list> - reduction factor
   .  $abslimit~<esc~double~list> - absolute limit for the defect (default 1E-10)
   .  $m~<maxit> - maximal number of iterations
   .  $I~<iteration> - iteration numproc
   .  $d~{full|red|no} - display modus
   .  $p~<con> - conjugate vector
   .  $t~<tmp> - temporaty vector
   .  $R~<restart> - restart index
   .  $w~<esc~double~list> - weighting factor

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

static INT EBCGSInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EBCGS *np=(NP_EBCGS *)theNP;
  INT i;

  if (esc_read (np->weight,NP_FMT(np),NULL,"weight",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] = 1.0;
  for (i=0; i<MAX_VEC_COMP; i++) np->weight[i] *= np->weight[i];
  np->B=ReadArgvEMatDesc(theNP->mg,"B",argc,argv);
  np->r=ReadArgvEVecDesc(theNP->mg,"r",argc,argv);
  np->p=ReadArgvEVecDesc(theNP->mg,"p",argc,argv);
  np->v=ReadArgvEVecDesc(theNP->mg,"v",argc,argv);
  np->s=ReadArgvEVecDesc(theNP->mg,"s",argc,argv);
  np->t=ReadArgvEVecDesc(theNP->mg,"t",argc,argv);
  np->q=ReadArgvEVecDesc(theNP->mg,"q",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("R",&(np->restart),argc,argv)) np->restart = 0;
  if (np->restart<0) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->display=ReadArgvDisplay(argc,argv);
  np->baselevel=0;
  np->Iter=(NP_EITER *)ReadArgvNumProc(theNP->mg,"I",EITER_CLASS_NAME,argc,argv);

  return (NPELinearSolverInit(&np->ls,argc,argv));
}

static INT EBCGSDisplay (NP_BASE *theNP)
{
  NP_EBCGS *np=(NP_EBCGS *)theNP;

  NPELinearSolverDisplay(&np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"R",(int)np->restart);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);
  if (np->Iter!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");
  if (np->display==PCR_NO_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display==PCR_RED_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display==PCR_FULL_DISPLAY) UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->B!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"B",ENVITEM_NAME(np->B));
  if (np->r!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->p!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  if (np->v!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"v",ENVITEM_NAME(np->v));
  if (np->s!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(np->s));
  if (np->t!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->q!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"q",ENVITEM_NAME(np->q));
  if (np->p!=NULL) if (esc_disp(np->weight,np->p,"weight")) REP_ERR_RETURN (1);

  return (0);
}

static INT EBCGSPreProcess (NP_ELINEAR_SOLVER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_EBCGS *np=(NP_EBCGS*)theNP;
  INT i;

  np->baselevel=MIN(*baselevel,level);
  if (np->Iter!=NULL)
    if (np->Iter->PreProcess != NULL)
    {
      if (np->B == NULL) if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);
      if (np->B != NULL) if ((*np->Iter->PreProcess)(np->Iter,level,x,b,np->B,baselevel,result)) REP_ERR_RETURN(1);
    }
  if (AllocEVDFromEVD(np->ls.base.mg,np->baselevel,level,x,&np->r)) NP_RETURN(1,result[0]);
  if (AllocEVDFromEVD(np->ls.base.mg,np->baselevel,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (AllocEVDFromEVD(np->ls.base.mg,np->baselevel,level,x,&np->v)) NP_RETURN(1,result[0]);
  if (AllocEVDFromEVD(np->ls.base.mg,np->baselevel,level,x,&np->s)) NP_RETURN(1,result[0]);
  if (AllocEVDFromEVD(np->ls.base.mg,np->baselevel,level,x,&np->t)) NP_RETURN(1,result[0]);
  if (AllocEVDFromEVD(np->ls.base.mg,np->baselevel,level,x,&np->q)) NP_RETURN(1,result[0]);

  for (i=0; i<EVD_NCOMP(x); i++) np->old_defect[i] = -1.0;

  return 0;
}

static INT EBCGSPostProcess (NP_ELINEAR_SOLVER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, INT *result)
{
  NP_EBCGS *np=(NP_EBCGS*)theNP;

  if (FreeEVD(np->ls.base.mg,np->baselevel,level,np->r)) REP_ERR_RETURN(1);
  if (FreeEVD(np->ls.base.mg,np->baselevel,level,np->p)) REP_ERR_RETURN(1);
  if (FreeEVD(np->ls.base.mg,np->baselevel,level,np->v)) REP_ERR_RETURN(1);
  if (FreeEVD(np->ls.base.mg,np->baselevel,level,np->s)) REP_ERR_RETURN(1);
  if (FreeEVD(np->ls.base.mg,np->baselevel,level,np->t)) REP_ERR_RETURN(1);
  if (FreeEVD(np->ls.base.mg,np->baselevel,level,np->q)) REP_ERR_RETURN(1);

  if (np->Iter!=NULL)
  {
    if (np->Iter->PostProcess != NULL)
    {
      if (np->B == NULL) if ((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result)) NP_RETURN(1,result[0]);
      if (np->B != NULL) if ((*np->Iter->PostProcess)(np->Iter,level,x,b,np->B,result)) NP_RETURN(1,result[0]);
    }
  }
  else
    return (0);

  np->baselevel=MAX(BOTTOMLEVEL(theNP->base.mg),np->baselevel);

  return 0;
}

static INT EBCGSSolver (NP_ELINEAR_SOLVER *theNP, INT level, EVECDATA_DESC *x, EVECDATA_DESC *b, EMATDATA_DESC *A, EVEC_SCALAR abslimit, EVEC_SCALAR reduction, ELRESULT *lresult)
{
  NP_EBCGS *np=(NP_EBCGS*)theNP;
  EVEC_SCALAR defect2reach;
  EMATDATA_DESC *DC;
  INT i,j,PrintID,restart,eq_count;
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
  for (i=0; i<EVD_NCOMP(x); i++)
  {
    NPELS_red(theNP)[i] = reduction[i];
    NPELS_abs(theNP)[i] = abslimit[i];
  }

  /* prepare */
  if (np->B==NULL) DC=A;
  else DC=np->B;
  alpha=rho_new=beta=tt=0.0;

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'*',"\n");
  if (np->display>PCR_NO_DISPLAY) if (PrepareEPCR(x,np->display,text,&PrintID)) NP_RETURN(1,lresult->error_code);
  if(NPELS_printbasetime(theNP)) basetime = 0.0;
  CSTART(); ti=0; ii=0;

  for (i=0; i<EVD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (esc_mul_check(defect2reach,lresult->first_defect,reduction,b)) NP_RETURN(1,lresult->error_code);
  if (np->display>PCR_NO_DISPLAY) if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
  if (esc_cmp(lresult->first_defect,abslimit,b)) lresult->converged=1;
  else lresult->converged=0;
  lresult->number_of_linear_iterations=0;

  /* go */
  restart=1;
  eq_count=0;
  for (i=0; i<np->maxiter; i++)
  {
    if (lresult->converged) break;

    /* restart ? */
    if ((np->restart>0 && i%np->restart==0) || restart)
    {
      if (deset(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->p,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (deset(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->v,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->r,b)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      alpha=np->rho=np->omega=1.0;
      restart = 0;
    }

    /* update x, b */
    if (dedotw(NP_MG(theNP),np->baselevel,level,ON_SURFACE,b,np->r,np->weight,&rho_new)!=NUM_OK) REP_ERR_RETURN (1);
    if ((np->rho!=0.0) && (np->omega!=0.0))
      beta = rho_new*alpha/np->rho/np->omega;
    if (descal(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->p,beta)) REP_ERR_RETURN (1);
    if (deadd(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->p,b)!=NUM_OK) REP_ERR_RETURN (1);
    if (deaxpy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->p,-beta*np->omega,np->v)!=NUM_OK) REP_ERR_RETURN (1);
    if (np->Iter!=NULL)
    {
      if (deset(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->q,0.0)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->s,np->p)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->p,DC,&lresult->error_code)) REP_ERR_RETURN (1);
      if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->p,np->s)!= NUM_OK) NP_RETURN(1,lresult->error_code);
      if (dematmul(NP_MG(theNP),np->baselevel,level,ON_SURFACE,np->v,A,np->q)) REP_ERR_RETURN (1);
            #ifdef ModelP
      if (a_vector_collect(NP_MG(theNP),np->baselevel,level,np->v->vd)!=NUM_OK) NP_RETURN(1,lresult->error_code);
            #endif
      if (dedotw(NP_MG(theNP),np->baselevel,level,ON_SURFACE,np->v,np->r,np->weight,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      if (alpha!=0.0)
        alpha=rho_new/alpha;
      if (deaxpy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,x,alpha,np->q)!= NUM_OK) REP_ERR_RETURN (1);
    }
    else
    {
            #ifdef ModelP
      if (a_vector_consistent(NP_MG(theNP),np->baselevel,level,np->p->vd)!=NUM_OK) NP_RETURN(1,lresult->error_code);
            #endif
      if (dematmul(NP_MG(theNP),np->baselevel,level,ON_SURFACE,np->v,A,np->p)) REP_ERR_RETURN (1);
            #ifdef ModelP
      if (a_vector_collect(NP_MG(theNP),np->baselevel,level,np->v->vd)!=NUM_OK) NP_RETURN(1,lresult->error_code);
            #endif
      if (dedotw(NP_MG(theNP),np->baselevel,level,ON_SURFACE,np->v,np->r,np->weight,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      if (alpha!=0.0)
        alpha=rho_new/alpha;
      if (deaxpy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,x,alpha,np->p)!= NUM_OK) REP_ERR_RETURN (1);
    }
    lresult->number_of_linear_iterations++;
    if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->s,b)!=NUM_OK) NP_RETURN(1,lresult->error_code);
    if (deaxpy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->s,-alpha,np->v)!=NUM_OK) REP_ERR_RETURN (1);
    if (ELinearResiduum(theNP,np->baselevel,level,x,np->s,A,lresult)) REP_ERR_RETURN (1);
    if (esc_cmp(lresult->last_defect,abslimit,b) || esc_cmp(lresult->last_defect,defect2reach,b))
    {
      if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,b,np->s) != NUM_OK) NP_RETURN(1,lresult->error_code);
      lresult->converged=1;
      if (np->display>PCR_NO_DISPLAY)
        if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
      break;
    }
    if (np->Iter!=NULL)
    {
      if (deset(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->q,0.0)!=NUM_OK) NP_RETURN(1,lresult->error_code);
      if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->t,np->s)!=NUM_OK) NP_RETURN(1,lresult->error_code);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->s,DC,&lresult->error_code)) REP_ERR_RETURN (1);
      if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->s,np->t)!=NUM_OK) NP_RETURN(1,lresult->error_code);
    }
    else
    {
      if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,np->q,np->s)!=NUM_OK) NP_RETURN(1,lresult->error_code);
            #ifdef ModelP
      if (a_vector_consistent(NP_MG(theNP),np->baselevel,level,np->q->vd)!=NUM_OK) NP_RETURN(1,lresult->error_code);
            #endif
    }
    if (dematmul(NP_MG(theNP),np->baselevel,level,ON_SURFACE,np->t,A,np->q)) REP_ERR_RETURN (1);
        #ifdef ModelP
    if (a_vector_collect(NP_MG(theNP),np->baselevel,level,np->t->vd)!=NUM_OK) NP_RETURN(1,lresult->error_code);
        #endif
    if (dedotw(NP_MG(theNP),np->baselevel,level,ON_SURFACE,np->t,np->t,np->weight,&tt)!=NUM_OK) REP_ERR_RETURN (1);
    if (dedotw(NP_MG(theNP),np->baselevel,level,ON_SURFACE,np->s,np->t,np->weight,&(np->omega))!=NUM_OK) REP_ERR_RETURN (1);
    PRINTDEBUG(np,2,("tt %f omega %f\n",tt,np->omega));
    if (tt!=0.0) np->omega /= tt;
    if (deaxpy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,x,np->omega,np->q)!=NUM_OK) REP_ERR_RETURN (1);
    if (decopy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,b,np->s)!=NUM_OK) NP_RETURN(1,lresult->error_code);
    if (deaxpy(NP_MG(theNP),np->baselevel,level,ALL_VECTORS,b,-np->omega,np->t)!=NUM_OK) REP_ERR_RETURN (1);
    np->rho=rho_new;

    /* residuum */
    if (ELinearResiduum(theNP,np->baselevel,level,x,b,A,lresult)) REP_ERR_RETURN (1);
    if (np->display>PCR_NO_DISPLAY)
      if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) NP_RETURN(1,lresult->error_code);
    lresult->number_of_linear_iterations++;
    if (esc_cmp(lresult->last_defect,abslimit,b) || esc_cmp(lresult->last_defect,defect2reach,b))
    {
      lresult->converged = 1;
      break;
    }
    if (esc_eq (lresult->last_defect,np->old_defect,1e-4,x)) eq_count++;
    else eq_count=0;
    for (j=0; j<EVD_NCOMP(x); j++) np->old_defect[j] = lresult->last_defect[j];
    if (eq_count>4)
    {
      lresult->converged = 0;
      break;
    }
  }
  CSTOP(ti,ii);
  if(NPELS_setbasetime(theNP)) basetime += ti;

  if (np->display>PCR_NO_DISPLAY)
  {
    if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) NP_RETURN(1,lresult->error_code);
    if (PostPCR(PrintID,":ls:avg")) NP_RETURN(1,lresult->error_code);
    if (SetStringValue(":ls:avg:iter",(DOUBLE) (i+1))) NP_RETURN(1,lresult->error_code);
    if (lresult->number_of_linear_iterations>0)
      if(NPELS_printbasetime(theNP))
        UserWriteF("BCGS: L=%2d N=%2d TSOLVE=%10.4g TIT=%10.4g TBASE=%g\n",level,lresult->number_of_linear_iterations,ti,ti/lresult->number_of_linear_iterations,basetime);
      else
        UserWriteF("BCGS: L=%2d N=%2d TSOLVE=%10.4g TIT=%10.4g\n",level,lresult->number_of_linear_iterations,ti,ti/lresult->number_of_linear_iterations);
  }

  return (0);
}

static INT EBCGSConstruct (NP_BASE *theNP)
{
  NP_EBCGS *np;

  theNP->Init = EBCGSInit;
  theNP->Display = EBCGSDisplay;
  theNP->Execute = NPELinearSolverExecute;

  np = (NP_EBCGS *) theNP;
  np->ls.PreProcess = EBCGSPreProcess;
  np->ls.Defect = ELinearDefect;
  np->ls.Residuum = ELinearResiduum;
  np->ls.Solver = EBCGSSolver;
  np->ls.PostProcess = EBCGSPostProcess;

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

INT NS_DIM_PREFIX InitELinearSolver ()
{
  if (CreateClass(ELINEAR_SOLVER_CLASS_NAME ".els",sizeof(NP_ELS),ELSConstruct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ELINEAR_SOLVER_CLASS_NAME ".ebcgs",sizeof(NP_EBCGS),EBCGSConstruct)) REP_ERR_RETURN (__LINE__);

  return (0);
}
