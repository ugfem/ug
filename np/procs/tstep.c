// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      tstep.c                                                           */
/*                                                                          */
/* Purpose:   time-step schemes                                                                 */
/*                                                                          */
/* Author:    Klaus Johannsen                                                                                           */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   Jan  31, 1998  begin                                          */
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
#include "debug.h"
#include "ugstruct.h"
#include "misc.h"
#include "gm.h"
#include "evm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"

#include "nls.h"
#include "ls.h"
#include "nls.h"
#include "assemble.h"
#include "transfer.h"

#include "tstep.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/*  Class Definition                                                                    */
/*                                                                          */
/****************************************************************************/

typedef struct
{
  NP_T_STEP tstep;                       /* derived from class NP_T_STEP    */

  /* nonlinear assemmly used */
  NP_NL_ASSEMBLE nlass;

  /* parameters (to be set with init function */
  INT baselevel;                         /* for nested iteration            */
  INT nested;                            /* use nested iteration            */
  INT displayMode;
  char scaleName[NAMELEN];
  DOUBLE scale;

  /* references to numprocs  */
  NP_TRANSFER *trans;                    /* uses transgrid for nested iter  */
  NP_T_ASSEMBLE *tass;                   /* time assemble numproc			*/
  NP_NL_SOLVER *nlsolve;                 /* nonlinear solver numproc		*/

} NP_BE;                                 /* backward euler                  */

/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* Nonlinear Assemble Interface provided to nonlinear solver				*/
/*                                                                          */
/****************************************************************************/

static NP_T_STEP *current_tstep;

INT BE_NLPreProcess (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x, INT *res)
{
  return(0);
}

INT BE_NLAssembleSolution (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, INT *res)
{
  NP_T_ASSEMBLE *tass;

  tass = ((NP_BE*)current_tstep)->tass;
  return((*tass->TAssembleSolution)(tass,fl,tl,current_tstep->t1,u,res));
}

INT BE_NLAssembleDefect (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  NP_BE *be;
  NP_T_ASSEMBLE *tass;
  DOUBLE s_a,s_m;

  /* get tass numproc */
  be = (NP_BE*)current_tstep;
  tass = be->tass;

  /* compute coefficients */
  s_a = -current_tstep->t1+current_tstep->t0;
  s_m = 1.0;

  /* assemble defect */
  dset(NP_MG(be),fl,tl,ALL_VECTORS,d,0.0);
  if ((*tass->TAssembleDefect)(tass,fl,tl,be->tstep.t0,-1.0,0.0,be->tstep.sol_t0,d,NULL,res)) NP_RETURN(1,*res);
  return ((*tass->TAssembleDefect)(tass,fl,tl,current_tstep->t1,s_m,s_a,u,d,J,res));
}

INT BE_NLAssembleMatrix (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *u, VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *J, INT *res)
{
  NP_BE *be;
  NP_T_ASSEMBLE *tass;
  DOUBLE s_a;

  /* get tass numproc */
  be = (NP_BE*)current_tstep;
  tass = be->tass;

  /* compute coefficients */
  s_a = -current_tstep->t1+current_tstep->t0;

  /* call function from time assemble interface */
  return ((*tass->TAssembleMatrix)(tass,fl,tl,current_tstep->t1,s_a,u,d,v,J,res));
}

INT BE_NLNAssembleMatrix (NP_NL_ASSEMBLE *ass, INT fl, INT tl, NODE *n, VECDATA_DESC *u, VECDATA_DESC *d, VECDATA_DESC *v, MATDATA_DESC *J, INT *res)
{
  NP_BE *be;
  NP_T_ASSEMBLE *tass;
  DOUBLE s_a;

  /* get tass numproc */
  be = (NP_BE*)current_tstep;
  tass = be->tass;

  /* compute coefficients */
  s_a = -current_tstep->t1+current_tstep->t0;

  /* call function from time assemble interface */
  return( (*tass->TNAssembleMatrix)(tass,fl,tl,n,current_tstep->t1,s_a,u,d,v,J,res) );
}

INT BE_NLPostProcess
  (NP_NL_ASSEMBLE *ass, INT fl, INT tl, VECDATA_DESC *x,
  VECDATA_DESC *d, MATDATA_DESC *J, INT *res)
{
  return(0);
}

/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* Implementation of NP_T_SOLVER functions                                                              */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/

INT BE_TimePreProcess (NP_T_STEP *tstep, INT level, INT *res)
{
  return(0);
}

INT BE_TimeInit (NP_T_STEP *tstep, INT level, DOUBLE t, VECDATA_DESC *sol, INT *res)
{
  NP_BE *be;
  NP_T_ASSEMBLE *tass;
  MULTIGRID *mg;
  char buffer[128];

  /* get numproc */
  be = (NP_BE *)tstep;
  tass = be->tass;
  mg = tass->base.mg;

  /* set parameters for debugging */
  be->tstep.t0 = t;
  be->tstep.sol_t0 = sol;

  /* set initial values and boundary conditions in sol */
  *res = 1;
  current_tstep = tstep;
  if (tass->TAssembleInitial!=NULL)
    if ((*tass->TAssembleInitial)(tass,0,level,t,sol,res)) return(1);
  if ((*tass->TAssembleSolution)(tass,0,level,t,sol,res)) return(1);

  /* write time to shell */
  sprintf(buffer,"%12.4lE",t);
  SetStringVar(":BE:T0",buffer);

  *res = 0;
  return(*res);
}

static INT BE_TimeStep (NP_T_STEP *tstep, INT level, DOUBLE t0, VECDATA_DESC *sol_t0, DOUBLE t1, VECDATA_DESC *sol_t1, TSTEP_RESULT *tstep_result)
{
  NP_BE *be;
  NP_T_ASSEMBLE *tass;
  NP_NL_SOLVER *nlsolve;
  NP_TRANSFER *trans;
  MULTIGRID *mg;
  INT i,ncomp,flevel,k,res;
  NLRESULT nlresult;
  DOUBLE Factor[MAX_VEC_COMP];
  char buffer[128],text[128];

  /* get numprocs */
  be = (NP_BE *)tstep;
  tass = be->tass;
  trans = be->trans;
  nlsolve = be->nlsolve;
  mg = nlsolve->base.mg;

  /* set parameters for debugging */
  be->tstep.t0 = t0;
  be->tstep.sol_t0 = sol_t0;
  be->tstep.t1 = t1;
  be->tstep.sol_t1 = sol_t1;

  /* get number of components */
  ncomp = VD_NCOMP(sol_t0);

  /* reset tstep_result */
  tstep_result->converged = 1;
  tstep_result->number_of_nonlinear_iterations = 0;
  tstep_result->number_of_linear_iterations = 0;
  tstep_result->max_linear_iterations = 0;
  tstep_result->exec_time = 0.0;

  /* determine level where to predict to new time step */
  if (be->nested) flevel = MIN(be->baselevel,level);
  else flevel = level;

  /* predict to new time step on level low */
  dcopy(mg,0,flevel,ALL_VECTORS,sol_t1,sol_t0);

  /* prepare transfer */
  if (trans->PreProcessProject!=NULL)
    if ((*trans->PreProcessProject)(trans,0,level,&res)) NP_RETURN(1,res);

  /* assemble */
  current_tstep = tstep;
  if ((*tass->TAssemblePreProcess)(tass,0,level,t1,t0,0.0,sol_t1,sol_t0,NULL,&res)) NP_RETURN(1,res);
  if ((*tass->TAssembleSolution)(tass,0,flevel,t1,sol_t1,&res)) NP_RETURN(1,res);

  /* do (nested) iteration on new time step */
  for (k=flevel; k<=level; k++)
  {
    if (be->nested) UserWriteF("Nested Iteration on level %d (%d)\n",k,level);

    /* prepare nonlinear solver on level k */
    if (nlsolve->PreProcess!=NULL)
      if ((*nlsolve->PreProcess)(nlsolve,k,sol_t1,&res)) NP_RETURN(1,res);

    /* solve nonlinear on level k */
    if ((*nlsolve->Solver)
          (nlsolve,k,sol_t1,&be->nlass,nlsolve->abslimit,nlsolve->reduction,&nlresult))
      return(__LINE__);

    /* postprocess nonlinear solver on level k */
    if (nlsolve->PostProcess!=NULL)
      if ((*nlsolve->PostProcess)(nlsolve,k,sol_t1,&res)) NP_RETURN(1,res);

    /* update result */
    tstep_result->number_of_nonlinear_iterations += nlresult.number_of_nonlinear_iterations;
    if (nlresult.number_of_nonlinear_iterations==0) tstep_result->jumped=1;
    else tstep_result->jumped=0;
    tstep_result->number_of_linear_iterations += nlresult.total_linear_iterations;
    tstep_result->max_linear_iterations = MAX(tstep_result->max_linear_iterations,nlresult.max_linear_iterations);
    tstep_result->exec_time += nlresult.exec_time;
    if (!nlresult.converged)
    {
      tstep_result->converged = 0;
      break;
    }

    /* are we ready ? */
    if (k==level) break;

    /* interpolate up */
    for (i=0; i<ncomp; i++) Factor[i] = 1.0;
    if ((*trans->InterpolateCorrection)(trans,k+1,sol_t1,sol_t1,NULL,Factor,&res)) NP_RETURN(1,res);

    /* set Dirichlet conditions in predicted solution */
    if ((*tass->TAssembleSolution)(tass,k+1,k+1,t1,sol_t1,&res)) NP_RETURN(1,res);
  }

  /* postprocess transfer */
  if (trans->PostProcessProject!=NULL)
    if ((*trans->PostProcessProject)(trans,0,level,&res)) NP_RETURN(1,res);

  /* postprocess assemble */
  if ((*tass->TAssemblePostProcess)(tass,0,level,t1,t0,0.0,sol_t1,sol_t0,NULL,&res)) NP_RETURN(1,res);

  /* write time to shell */
  sprintf(buffer,"%12.4lE",t0);
  SetStringVar(":BE:T0",buffer);
  sprintf(buffer,"%12.4lE",t1);
  SetStringVar(":BE:T1",buffer);

  if (be->displayMode == PCR_RED_DISPLAY || be->displayMode == PCR_FULL_DISPLAY )
  {
    UserWriteF("\n");
    CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(be),'%',"\n");  UserWrite(text);
    UserWriteF("Time:          t0: %e %s\n",(float)t0/be->scale,be->scaleName);
    UserWriteF("               t1: %e %s\n",(float)t1/be->scale,be->scaleName);
    UserWriteF("\n");
  }
  if (be->displayMode == PCR_FULL_DISPLAY)
  {
    UserWriteF("Exec:           t: %e sec.\n",(float)tstep_result->exec_time);
    UserWriteF("               nl: %4d it.\n",(int)tstep_result->number_of_nonlinear_iterations);
    UserWriteF("              lin: %4d it.\n",(int)tstep_result->number_of_linear_iterations);
    UserWriteF("       max lin/nl: %4d it.\n",(int)tstep_result->max_linear_iterations);
    UserWriteF("\n");
  }

  return(0);
}

INT BE_TimePostProcess (NP_T_STEP *tstep, INT level, INT *res)
{
  return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  Init															*/
/*																			*/
/* Purpose:   init solve cycle												*/
/*																			*/
/* Input:         NumProcType Init Function										*/
/*																			*/
/* Output: INT 0: ok														*/
/*			   1: error														*/
/*																			*/
/****************************************************************************/

INT BE_Init (NP_BASE *base, INT argc, char **argv)
{
  NP_BE *be;
  VECDATA_DESC *tmp;
  INT ret;
  MULTIGRID *mg;

  /* get numprocs ... */
  be = (NP_BE *)base;
  mg = base->mg;

  /* read other numprocs */
  be->tass = (NP_T_ASSEMBLE *) ReadArgvNumProc(mg,"A",T_ASSEMBLE_CLASS_NAME,argc,argv);
  if (be->tass==NULL) return (NP_NOT_ACTIVE);
  be->nlsolve = (NP_NL_SOLVER *) ReadArgvNumProc(mg,"S",NL_SOLVER_CLASS_NAME,argc,argv);
  if (be->nlsolve==NULL) return (NP_NOT_ACTIVE);
  be->trans = (NP_TRANSFER *) ReadArgvNumProc(mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  if (be->trans==NULL) return(NP_NOT_ACTIVE);

  /* set configuration parameters */
  if (ReadArgvINT("baselevel",&(be->baselevel),argc,argv)) be->baselevel=0;
  if ((be->baselevel<0)||(be->baselevel>32)) return(NP_NOT_ACTIVE);
  if (ReadArgvINT("nested",&(be->nested),argc,argv)) be->nested=0;
  if ((be->nested<0)||(be->nested>1)) return(NP_NOT_ACTIVE);
  if (ReadArgvChar("scale",be->scaleName,argc,argv))
  {
    be->scale=1.0;
    be->scaleName[0]='\0';
  }
  else
  {
    if (strcmp(be->scaleName,"second")==0) be->scale = 1.0;
    else if (strcmp(be->scaleName,"minute")==0) be->scale = 60.0;
    else if (strcmp(be->scaleName,"hour")==0) be->scale = 3600.0;
    else if (strcmp(be->scaleName,"day")==0) be->scale = 86400.0;
    else if (strcmp(be->scaleName,"week")==0) be->scale = 604800;
    else if (strcmp(be->scaleName,"month")==0) be->scale = 2628000;
    else if (strcmp(be->scaleName,"year")==0) be->scale = 31536000;
    else
    {
      UserWrite("ERROR: cannot read scale-option\n");
      return(NP_NOT_ACTIVE);
    }
  }
  be->displayMode=ReadArgvDisplay(argc,argv);
  ret=NP_EXECUTABLE;

  /* read configuration for executability */
  be->tstep.sol_t0 = ReadArgvVecDesc(base->mg,"sol",argc,argv);
  if (be->tstep.sol_t0==NULL) ret=NP_ACTIVE;
  if (ReadArgvDOUBLE("t0",&(be->tstep.t0),argc,argv)) ret=NP_ACTIVE;
  if (ReadArgvDOUBLE("t1",&(be->tstep.t1),argc,argv)) ret=NP_ACTIVE;

  return (ret);
}

/****************************************************************************/
/*																			*/
/* Function:  Display														*/
/*																			*/
/* Purpose:   display linear multigrid cycle								*/
/*																			*/
/* Input:     none															*/
/*																			*/
/* Output:    INT 0: ok														*/
/*			   1: error														*/
/*																			*/
/****************************************************************************/

INT BE_Display (NP_BASE *theNumProc)
{
  NP_BE *be;

  /* get numprocs */
  be = (NP_BE *)theNumProc;

  UserWrite("\nBE configuration:\n");
  if (be->tass!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(be->tass));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"A","---");
  if (be->nlsolve!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(be->nlsolve));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"S","---");
  if (be->trans!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(be->trans));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");

  UserWriteF(DISPLAY_NP_FORMAT_SF,"t0",(float)be->tstep.t0);
  if (be->tstep.sol_t0!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_t0",ENVITEM_NAME(be->tstep.sol_t0));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_t0","---");
  UserWriteF(DISPLAY_NP_FORMAT_SF,"t1",(float)be->tstep.t1);
  if (be->tstep.sol_t1!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_t1",ENVITEM_NAME(be->tstep.sol_t1));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"sol_t1","---");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"nested",(int)be->nested);
  if (be->displayMode==PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (be->displayMode==PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (be->displayMode==PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  Execute                                                       */
/*                                                                          */
/* Purpose:   init solve cycle                                              */
/*                                                                          */
/* Input:     NumProcType Init Function                                     */
/*                                                                          */
/* Output: INT 0: ok                                                        */
/*             1: error                                                     */
/*                                                                          */
/****************************************************************************/

static INT TSTEP_Execute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_T_STEP *tstep;
  INT level,result;
  TSTEP_RESULT tresult;
  MULTIGRID *mg;
  DOUBLE dt;

  /* get numprocs */
  tstep = (NP_T_STEP *)theNP;

  /* get level */
  mg = theNP->mg;
  level = CURRENTLEVEL(mg);

  /* call preprocess function */
  if (ReadArgvOption("pre",argc,argv))
    if (tstep->TimePreProcess!=NULL)
      if ((*tstep->TimePreProcess)(tstep,level,&result))
      {
        UserWriteF("TSTEP_Execute: TimePreProcess failed, error code %d\n",result);
        return (1);
      }

  /* set initial values */
  if (ReadArgvOption("init",argc,argv))
    if (tstep->TimeInit != NULL)
      if ((*tstep->TimeInit)(tstep,level,tstep->t0,tstep->sol_t0,&result))
      {
        UserWriteF("TSTEP_Execute: TimeInit failed, error code %d\n",result);
        return (1);
      }

  /* execute bdf1, nonnested */
  if (ReadArgvOption("step",argc,argv))
    if (tstep->TimeStep!=NULL)
    {
      /* allocate tmp vector */
      if (AllocVDFromVD(mg,0,level,tstep->sol_t0,&(tstep->sol_t1))) return (1);

      /* calculate solution */
      if ((*tstep->TimeStep)(tstep,level,tstep->t0,tstep->sol_t0,tstep->t1,tstep->sol_t1,&tresult))
      {
        UserWriteF("TSTEP_Execute: TimeStep failed, error code\n");
        return (1);
      }
      if (!tresult.converged)
      {
        UserWriteF("TSTEP_Execute: TimeInit failed, cannot calculate solution at t1\n");
        return (1);
      }
      else
      {
        dcopy(mg,0,level,ALL_VECTORS,tstep->sol_t0,tstep->sol_t1);
        dt = tstep->t1-tstep->t0;
        tstep->t0=tstep->t1;
        tstep->t1+=dt;
      }

      /* free tmp-vector */
      if (FreeVD(mg,0,level,tstep->sol_t1)) return (1);
    }

  /* postprocess */
  if (ReadArgvOption("post",argc,argv))
    if (tstep->TimePostProcess!=NULL)
      if ((*tstep->TimePostProcess)(tstep,level,&result))
      {
        UserWriteF("TSTEP_Execute: TimePostProcess failed, error code %d\n", result);
        return (1);
      }

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  BDFConstruct													*/
/*																			*/
/* Purpose:	  create an object of class bdf									*/
/*																			*/
/* Input:	  theNP	num proc data structure									*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*			  else : error													*/
/*																			*/
/****************************************************************************/

static INT BEConstruct (NP_BASE *theNP)
{
  NP_T_STEP *tstep;
  NP_NL_ASSEMBLE *na;
  NP_BE *be;

  /* get numprocs */
  tstep = (NP_T_STEP *)theNP;
  be = (NP_BE *)theNP;
  na = &be->nlass;

  /* set base functions */
  theNP->Init    = BE_Init;
  theNP->Display = BE_Display;
  theNP->Execute = TSTEP_Execute;

  /* nonlinear assemble functions */
  na->PreProcess                      = BE_NLPreProcess;
  na->PostProcess             = BE_NLPostProcess;
  na->NLAssembleSolution      = BE_NLAssembleSolution;
  na->NLAssembleDefect        = BE_NLAssembleDefect;
  na->NLAssembleMatrix        = BE_NLAssembleMatrix;
  na->NLNAssembleMatrix       = BE_NLNAssembleMatrix;

  /* and the time solver */
  tstep->TimePreProcess   = BE_TimePreProcess;
  tstep->TimeInit                 = BE_TimeInit;
  tstep->TimeStep                 = BE_TimeStep;
  tstep->TimePostProcess  = BE_TimePostProcess;

  be->tstep.t0 = 0.0;
  be->tstep.sol_t0 = NULL;
  be->tstep.t1 = 0.0;
  be->tstep.sol_t1 = NULL;

  return(0);
}

/****************************************************************************/
/*D
   be - backward-euler time stepping scheme

   DESCRIPTION:
   This numproc executes a BDF time stepping scheme of first and second order.
   In every time step it calls a nonlinear solver, where the linearization
   in the nonlinear solver is constructed from a time dependent assembling
   routine of type 'NP_T_ASSEMBLE'.

   .vb
   npinit <name> $y <sol> $A <tnlass> $S <nlsolver>
          $T <transfer> [$baselevel <bl>] $order {1|2}
                  $predictorder {0|1} $nested {0|1}
          $tstart <time> $dtstart <step> $dtmin <min step> $dtmax <max step>
          $dtscale <fac> $rhogood <rho> [$J <mat>];
   .ve

   .  $y~<sol>            - solution vector
   .  $A~<tnlass>         - assemble numproc of type 'NP_T_ASSEMBLE'
   .  $S~<nlsolver>       - nonlinear solver numproc of type 'NP_NL_SOLVER'
   .  $T~<transfer>       - transfer numproc
   .  $baselevel~<bl>     - baselevel for nested iteration
   .  $order~{1|2}        - order of the time stepping scheme
   .  $predictorder~{0|1} - order of the predictor for the next start vector
   .  $nested~{0|1}       - flag nested nonlinear solver in every time step
   .  $tstart~<time>      - start time
   .  $dtstart~<step>     - time step to begin with
   .  $dtmin~<min step>   - smallest time step allowed
   .  $dtmax~<max step>   - largest time step allowed
   .  $dtscale~<fac>      - scaling factor applied after step
   .  $rhogood~<rho>      - threshold for step doubling

   'npexecute <name> [$pre] [$init] [$bdf1|$bdf1n|$bdf2|$bdf2n] [$post];'

   .  $pre - preprocess
   .  $init - initialization
   .  $bdf1 - backward Euler of first order
   .  $bdf1n - nested backward Euler of first order
   .  $bdf1 - backward Euler of second order
   .  $bdf1n - nested backward Euler of second order
   .  $post - postprocess

   EXAMPLE:
   .vb
 # grid transfer numproc
   npcreate transfer $c transfer;
   npinit transfer;

 # assemble numproc
   npcreate tnlass $c boxtnl;
   npinit tnlass;

 # linear solver and iteration numprocs
   npcreate smooth $c gs;
   npinit smooth;

   npcreate basesolver $c ls;
   npinit basesolver $red 0.001 $m 10 $I smooth;

   npcreate lmgc $c lmgc;
   npinit lmgc $S smooth smooth basesolver $T transfer $n1 2 $n2 2;

   npcreate mgs $c ls;
   npinit mgs $m 40 $I lmgc $display full;

 # nonlinear solver numproc to be used by time solver
   npcreate newton $c newton;
   npinit newton $red 1.0E-8 $T transfer $S mgs
              $rhoreass 0.8 $lsteps 6 $maxit 50 $line 1 $linrate 0
              $lambda 1.0 $linminred 1.0E-2 $display full;

 # the time solver
   npcreate ts $c bdf;
   npinit ts $y sol $A tnlass $S newton
          $T transfer $baselevel 0 $order 1 $predictorder 0 $nested 0
          $dtstart 0.00125 $dtmin 0.001 $dtmax 0.00125 $dtscale 1.0
          $rhogood 0.001;

 # first time step
   TIMESTEPS = 100;
   step = 1;
   npexecute ts $pre $init;
   npexecute ts $bdf1;

 # time step loop
   repeat {
    if (step==TIMESTEPS) break;
    step = step+1;
    npexecute ts $bdf2;
   }
   npexecute ts $post;
   .ve
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   InitTStep - Init this file

   SYNOPSIS:
   INT InitTStep (void);

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

INT InitTStep (void)
{
  INT error=0;

  if (MakeStruct(":BE")!=0)
  {
    SetHiWrd(error,__LINE__);
    return (1);
  }

  if (CreateClass (T_STEP_CLASS_NAME ".be",sizeof(NP_BE), BEConstruct)) return (__LINE__);

  return (0);
}
