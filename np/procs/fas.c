// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  fas.c													        */
/*																			*/
/* Purpose:   nonlinear multigrid solver (nls type)			                            */
/*																			*/
/* Author:	  Gabriele Beddies                                                                                              */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   28.07.97 begin, ug version 3.8								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#define __COMPILE_NMGSOLVER__

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "ugdevices.h"
#include "general.h"
#include "debug.h"
#include "gm.h"
#include "scan.h"
#include "misc.h"
#include "numproc.h"
#include "pcr.h"
#include "shapes.h"
#include "np.h"

#include "nls.h"
#include "assemble.h"

#include "nliter.h"
#include "transfer.h"
#include "fas.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define CSTART()    clock_start=CURRENT_TIME_LONG;
#define CSTOP(t,c)  t+=(CURRENT_TIME_LONG-clock_start);c++

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static DOUBLE Factor_One[MAX_VEC_COMP];
static DOUBLE Factor_Minus_One[MAX_VEC_COMP];

/* variables for timing measurement		*/
static int fas_c;
static double fas_t;
static DOUBLE clock_start;

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
  NP_NL_SOLVER nlsolver;                /* derived from abstract class NP_NL_SOLVER	*/

  /* parameters to be set via npinit */
  NP_TRANSFER *trans;                           /* uses transgrid						        */
  NP_NL_ITER *nliter;           /* uses nonlinear iteration                 */

  INT displayMode;                                     /* for PCR					                */
  INT baselevel;                                       /* baselevel					                */
  INT gamma;                                                   /* gamma						                */
  INT nu1;                                 /* number of pre-smoothing steps     */
  INT nu2;                                 /* number of post-smoothing steps    */
  INT niter;                               /* number of basesolver steps        */
  INT maxit;                           /* maximal number of iterations      */
  DOUBLE damp[MAX_VEC_COMP];           /* damp factor for solution		    */
  DOUBLE restriction[MAX_VEC_COMP];        /* restrict factor for solution              */

  /* and XDATA_DESCs */
  MATDATA_DESC *J;                              /* the Matrix to be solved				        */
  VECDATA_DESC *l;                              /* solution of last iterate				        */
  VECDATA_DESC *v;                              /* nonlinear correction						*/
  VECDATA_DESC *d;                              /* nonlinear defect						        */

} NP_FAS;

INT FasStep (NP_FAS *fas, NP_NL_ASSEMBLE *ass, INT level, VECDATA_DESC *x);

/****************************************************************************/
/*D
   RestrictSolNodeVector - Restrict solution of fine node vectors

   SYNOPSIS:
   static INT RestrictSolNodeVector (GRID *FineGrid,
   const VECDATA_DESC *to,
   const VECDATA_DESC *from, const DOUBLE *damp);

   PARAMETERS:
   .  FineGrid - pointer to grid
   .  to - type vector descriptor
   .  from  - type vector descriptor
   .  damp - damping factor for every component

   DESCRIPTION:
   This function restricts solution of fine node vectors with NEWDEFECT_CLASS
   to the next coarser grid. It is the transposent operation to
   'StandardIntCorNodeVector'.
   First, resets all components to zero if vecskip==0.
   It considers the VECSKIP-flags.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
   D*/
/****************************************************************************/

static INT RestrictSolNodeVector (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  GRID *CoarseGrid;
  ELEMENT *theElement;
  VERTEX *theVertex;
  NODE *theNode;
  VECTOR *v,*vc;
  DOUBLE c[MAX_CORNERS_OF_ELEM],s[MAX_SINGLE_VEC_COMP];
  const SHORT *toComp,*fromComp;
  INT i,j,n,ncomp,vecskip,dt;

  CoarseGrid = DOWNGRID(FineGrid);

  toComp    = VD_ncmp_cmpptr_of_otype_mod(to,NODEVEC,&ncomp,NON_STRICT);
  fromComp  = VD_cmpptr_of_otype_mod(from,NODEVEC,NON_STRICT);
  if (ncomp <= 0)
    return(NUM_ERROR);
  if (ncomp>MAX_SINGLE_VEC_COMP)
    return (NUM_BLOCK_TOO_LARGE);

  /* reset coarser value at positions where a new value is restricted */
  dt = VD_DATA_TYPES(to);
  for (v=FIRSTVECTOR(CoarseGrid); v!= NULL; v=SUCCVC(v))
    if ((VOTYPE(v)==NODEVEC) && V_IN_DATATYPE(v,dt))
      if (VCLASS(v)>=NEWDEF_CLASS)
        for (i=0; i<ncomp; i++)
          VVALUE(v,toComp[i]) = 0.0;

  /* compute contributions to all coarse node vectors */
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    v = NVECTOR(theNode);
    if (VCLASS(v)<NEWDEF_CLASS) continue;
    if (!V_IN_DATATYPE(v,dt)) continue;

    if (CORNERTYPE(theNode))
    {
      vc = NVECTOR((NODE*)NFATHER(theNode));
      vecskip = VECSKIP(vc);
      for (i=0; i<ncomp; i++)
      {
        if (!(vecskip & (1<<i)))
          VVALUE(vc,toComp[i]) = damp[i] * VVALUE(v,fromComp[i]);
      }
      for (i=0; i<ncomp; i++)
        if (vecskip)
          VVALUE(vc,toComp[i]) = VVALUE(v,fromComp[i]);
    }
    /*else
       {
            theVertex  = MYVERTEX(theNode);
            theElement = VFATHER(theVertex);
            n = CORNERS_OF_ELEM(theElement);
            GNs(n,LCVECT(theVertex),c);
            for (i=0; i<ncomp; i++)
                    s[i] = damp[i] * VVALUE(v,fromComp[i]);
            for (i=0; i<n; i++)
            {
                    vc = NVECTOR(CORNER(theElement,i));
                    vecskip = VECSKIP(vc);
                    for (j=0; j<ncomp; j++)
                            if (!(vecskip & (1<<j)))
                                    VVALUE(vc,toComp[j]) += c[i] * s[j];
            }
            for (i=0; i<ncomp; i++)
                    if (vecskip)
                            VVALUE(vc,toComp[i]) = VVALUE(v,fromComp[i]);
       }*/
  }

  return (0);
}

/****************************************************************************/
/*D
   RestrictValue - Restrict solution of fine vectors with NEWDEFECT_CLASS

   SYNOPSIS:
   INT RestrictValue (GRID *FineGrid, const VECDATA_DESC *to,
   const VECDATA_DESC *from, const DOUBLE *damp);

   PARAMETERS:
   .  FineGrid - pointer to grid
   .  to - type vector descriptor
   .  from  - type vector descriptor
   .  damp - damping factor for every component

   DESCRIPTION:
   This function restricts solution of fine vectors with NEWDEFECT_CLASS,
   considers the VECSKIP-flags.
   It calls 'RestrictSolElemVector', 'RestrictSolNodeVector',
   'RestrictSolEdgeVector' and 'RestrictSolSideVector',
   depending on the type vector descriptor.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
   D*/
/****************************************************************************/

INT RestrictValue (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  FORMAT *fmt;
  INT vtype,rv,otype;
  const SHORT *offset;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  offset = VD_OFFSETPTR(to);
  fmt = MGFORMAT(MYMG(FineGrid));

  for (otype=0; otype<MAXVOBJECTS; otype++)
    if (VD_OBJ_USED(to) & BITWISE_TYPE(otype))
      switch (otype)
      {
      case ELEMVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      case NODEVEC :
        for (vtype=0; vtype<NVECTYPES; vtype++)
          if (VD_ISDEF_IN_TYPE(to,vtype))
            if (GetUniqueOTypeOfVType(fmt,vtype)<0)
              REP_ERR_RETURN(1)
              if (RestrictSolNodeVector(FineGrid,to,from,damp+offset[otype]))
                return (1);
        break;
      case EDGEVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      case SIDEVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      }

  return (NUM_OK);
}
/****************************************************************************/
/*D
        fas - nonlinear solver numproc

   DESCRIPTION:

   .vb
   npinit <name> $x <sol>
              $A <assemble> $T <transfer> $S <solver>
              $abslimit <sc double list> $red <sc double list>
              $maxit <m>
                          [$display {no|red|full}] [$J <mat>];
   .ve

   .  $x~<sol>           - solution vector
   .  $A~<assemble>      - assemble numproc of type 'NP_NL_ASSEMBLE'
   .  $S~<linear solver> - linear solver numproc of type 'NP_LINEAR_SOLVER'
   .  $T~<transfer>      - transfer numproc
   .  $abslimit~<sc~double~list>  - absolute limit for the defect (default 1.0E-10)
   .  $red~<sc~double~list> - reduction factor
   .  $maxit~<m>            - maximum number of nonlinear iterations
   .  $display~{no|red|full}] - display mode
   .  $J~<mat>                - Jacobi matrix (optional)

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n   nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata
   .n   if only a single value is specified, this will be used for all components

   'npexecute <name> [$i] [$s] [$p];

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   EXAMPLE:
   .vb
 # grid transfer numproc
   npcreate transfer $c transfer;
   npinit transfer;

 # assemble numproc
   npcreate nlass $c tldf;
   npinit nlass;

 # nonlinear iteration numproc
   npcreate smooth $c nlgs;
   npinit smooth $damp 0.92;

 # nonlinear solver numproc
   npcreate fas $c fas;
   npinit fas $x sol
              $A nlass $T transfer $S mgs
              $abslimit 1.0E-10 $red 1.0E-5
              $maxit 50
              $display full;

   npexecute fas $i $s $p;
   .ve
   D*/
/****************************************************************************/

static INT FasInitialStep (MULTIGRID *mg, INT level, INT init, VECDATA_DESC *x,
                           NP_FAS *fas, NP_NL_ASSEMBLE *ass, VEC_SCALAR defect)
{
  INT n_unk;
  INT error;

  /* get number of components */
  n_unk = VD_NCOMP(x);

  /*if (fas->trans->PreProcess != NULL)
          if ((*fas->trans->PreProcess)
                  (fas->trans,&(fas->baselevel),level,x,fas->d,ass->A,&error)) {
      error = __LINE__;
          REP_ERR_RETURN(error);
     }*/

  /* project solution to all grid levels */
  if (fas->trans->PreProcessProject!=NULL)
    if ((*fas->trans->PreProcessProject)
          (fas->trans,0,level,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }
  if ((*fas->trans->ProjectSolution)(fas->trans,0,level,x,&error)) {
    error = __LINE__;
    REP_ERR_RETURN(error);
  }
  if (fas->trans->PostProcessProject!=NULL)
    if ((*fas->trans->PostProcessProject)
          (fas->trans,0,level,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }

  /* set initial values */
  if (init)
  {
    /* preprocess assemble once before all calls */
    if (ass->PreProcess!=NULL)
      if ((*ass->PreProcess)(ass,0,level,x,&error)) {
        error = __LINE__;
        REP_ERR_RETURN(error);
      }

    /* set dirichlet conditions on all grid levels */
    if ((*ass->NLAssembleSolution)(ass,0,level,x,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }
  }

  /* compute new nonlinear defect */
  dset(mg,0,level,ALL_VECTORS,fas->d,0.0);
  if ((*ass->NLAssembleDefect)(ass,0,level,x,fas->d,ass->A,&error)) {
    error = __LINE__;
    REP_ERR_RETURN(error);
  }

  /* compute norm of defect */
        #ifdef ModelP
  if (a_vector_collect(mg,0,level,fas->d)) {
    error = __LINE__;
    REP_ERR_RETURN(error);
  }
        #endif
  if (dnrm2x(mg,0,level,ON_SURFACE,fas->d,defect)) {
    error = __LINE__;
    REP_ERR_RETURN(error);
  }

  return(0);
}

static INT FasPreProcess  (NP_NL_SOLVER *solve, INT level, VECDATA_DESC *x, INT *result)
{
  NP_FAS *fas;
  NP_NL_ASSEMBLE *ass;
  MULTIGRID *mg;
  GRID *g;
  NODE *theNode;
  ELEMENT *theElement;
  INT i,k;
  INT baselevel;

  fas = (NP_FAS *) solve;
  ass = (NP_NL_ASSEMBLE *) solve;

  /* read multigrid */
  mg = solve->base.mg;
  if (AllocVDFromVD(mg,0,level,x,&(fas->l)))
    NP_RETURN(1,result[0]);
  if (AllocMDFromVD(mg,0,level,x,x,&(fas->J)))
    NP_RETURN(1,result[0]);

  if (ass->A == NULL)
    ass->A = fas->J;

  if (fas->nliter->PreProcess != NULL)
    for (i = fas->baselevel+1; i <= level; i++)
      if ((*fas->nliter->PreProcess )
            (fas->nliter,i,x,fas->d,ass->A,&fas->baselevel,result))
        REP_ERR_RETURN (1);

  baselevel = MIN(fas->baselevel,level);
  if (fas->nliter->PreProcess != NULL)
    if ((*fas->nliter->PreProcess)
          (fas->nliter,baselevel,x,fas->d,ass->A,&fas->baselevel,result))
      REP_ERR_RETURN(1);

  /* create ElementList */
  for (k=0; k<=level; k++)
  {
    g = GRID_ON_LEVEL(mg,k);

    /*for (theNode=FIRSTNODE(g); theNode!= NULL; theNode=SUCCN(theNode))
            NODE_ELEMENT_LIST(theNode) = NULL;*/

    for (theElement=FIRSTELEMENT(g); theElement!=NULL; theElement=SUCCE(theElement))
    {
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        theNode = CORNER(theElement,i);
        if(CreateElementList(g,theNode,theElement)) return(__LINE__);
      }
    }
  }

  return(0);
}


static INT FasPostProcess (NP_NL_SOLVER *solve, INT level, VECDATA_DESC *x, INT *result)
{
  NP_FAS *fas;
  NP_NL_ASSEMBLE *ass;
  INT i;

  fas = (NP_FAS *) solve;
  ass = (NP_NL_ASSEMBLE *) solve;

  FreeVD(solve->base.mg,0,level,fas->l);
  FreeMD(solve->base.mg,0,level,fas->J);

  if (fas->nliter->PostProcess != NULL)
    for (i = fas->baselevel+1; i <= level; i++)
      if ((*fas->nliter->PostProcess)
            (fas->nliter,level,x,fas->d,ass->A,result))
        REP_ERR_RETURN (1);

  if (fas->nliter->PostProcess != NULL)
    if ((*fas->nliter->PostProcess)
          (fas->nliter,fas->baselevel,x,fas->d,ass->A,result))
      REP_ERR_RETURN(1);

  if (fas->trans->PostProcess != NULL)
    if ((*fas->trans->PostProcess)
          (fas->trans,&(fas->baselevel),level,x,fas->d,ass->A,result))
      REP_ERR_RETURN(1);

  return(0);
}

static INT FasSolver (NP_NL_SOLVER *nls, INT level, VECDATA_DESC *x,
                      NP_NL_ASSEMBLE *ass, VEC_SCALAR abslimit, VEC_SCALAR reduction,
                      NLRESULT *res)
{
  NP_FAS *fas;
  MULTIGRID *mg;
  INT i,j, PrintID;
  VEC_SCALAR defect, defect2reach;
  char text[DISPLAY_WIDTH+4];
  INT n_unk;

  /* get status */
  fas = (NP_FAS *) nls;      /* cast from abstract base class to final class*/
  mg = nls->base.mg;

  /* fill result variable with error condition */
  res->error_code = 0;
  res->converged = 0;
  res->number_of_nonlinear_iterations = 0;
  res->exec_time = 0.0;

  /* initialize timers and counters */
  fas_c = 0;
  fas_t = 0.0;

  /* check function pointers in numprocs */
  if (ass->NLAssembleSolution==NULL)
  {
    UserWrite("Fas: ass->NLAssembleSolution not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  if (ass->NLAssembleDefect==NULL)
  {
    UserWrite("Fas: ass->NLAssembleDefect not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  if (ass->NLAssembleMatrix==NULL)
  {
    UserWrite("Fas: ass->NLAssembleMatrix not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }
  if (ass->NLNAssembleMatrix==NULL)
  {
    UserWrite("Fas: ass->NLNAssembleMatrix not defined\n");
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }

  /* dynamic XDATA_DESC allocation */
  if (ass->A == NULL)
    ass->A = fas->J;
  if (AllocVDFromVD(mg,0,level,x,  &(fas->v)))
  {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}
  if (AllocVDFromVD(mg,0,level,x,  &(fas->d)))
  {res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);}

  /* get number of components */
  n_unk = VD_NCOMP(x);

  for (i=0; i<n_unk; i++) {
    Factor_One[i] = 1.0;
    Factor_Minus_One[i] = -1.0;
  }

  /* init ass once and compute nonlinear defect */
  if (FasInitialStep (mg,level,TRUE,x,fas,ass,defect)!=0)
  {
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }

  /* print norm of defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(fas),'#',"\n");
  if (PreparePCR(fas->d,fas->displayMode,text,&PrintID)) {
    res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);
  }
  if (sc_mul(defect2reach,defect,reduction,fas->d)) {
    res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);
  }
  if (DoPCR(PrintID,defect,PCR_CRATE)) {
    res->error_code = __LINE__; REP_ERR_RETURN(res->error_code);
  }
  for (i=0; i<n_unk; i++) res->first_defect[i] = defect[i];

  /* check if iteration is necessary */
  if (sc_cmp(defect,abslimit,fas->d)) {
    res->converged = 1;
    for (i=0; i<n_unk; i++) res->last_defect[i] = defect[i];
    res->error_code = 0;
    goto exit;
  }

  /* run cycles */
  for (i=0; i<fas->maxit; i++)
  {
    if (res->converged) break;

    /* set correction to zero */
    if (dset(mg,0,level,ALL_VECTORS,fas->v,0.0))
      return (1);

    CSTART();
    /* iterate */
    if (FasStep(fas,ass,level,x))
      return (1);
    CSTOP(fas_t,fas_c);

    /* compute nonlinear defect */
    if (FasInitialStep (mg,level,FALSE,x,fas,ass,defect)!=0)
    {
      res->error_code = __LINE__;
      REP_ERR_RETURN(res->error_code);
    }


    /* print norm of defect */
    if (DoPCR(PrintID,defect,PCR_CRATE)) {
      res->error_code = __LINE__;
      REP_ERR_RETURN(res->error_code);
    }

    /* check if limit reached */
    if (sc_cmp(defect,abslimit,fas->d)) {res->converged=1;  break;}
    if (sc_cmp(defect,defect2reach,fas->d)) {res->converged=1;      break;}
  }

  /* print norm of defect */
  if (DoPCR(PrintID,defect,PCR_AVERAGE)) {
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }

  /* if converged, then report results */
  if (res->converged) {
    res->error_code = 0;
    res->number_of_nonlinear_iterations = fas_c;
    res->exec_time = fas_t;
  }

exit:
  if (PostPCR(PrintID,NULL)) {
    res->error_code = __LINE__;
    REP_ERR_RETURN(res->error_code);
  }

  /* deallocate local XDATA_DESCs */
  FreeVD(mg,0,level,fas->v);
  FreeVD(mg,0,level,fas->d);
  if (res->error_code==0)
    return (0);
  else
    REP_ERR_RETURN (res->error_code);
}

INT FasStep (NP_FAS *fas, NP_NL_ASSEMBLE *ass, INT level,
             VECDATA_DESC *x)
{
  VECTOR *theVector;
  NODE *theNode;
  MULTIGRID *mg;
  GRID *g;
  DOUBLE damp_factor[MAX_VEC_COMP];
  INT i,error;
  INT n_unk;

  fas->nlsolver.Assemble = ass;

  /* get number of components */
  n_unk = VD_NCOMP(x);

  for (i=0; i<n_unk; i++) damp_factor[i] = -1.0*fas->damp[i];

  /* get multigrid and grid */
  mg = NP_MG(fas);
  g = GRID_ON_LEVEL(mg,level);

  if (level<=fas->baselevel)
  {
    /* keep value */
    if (dcopy(mg,level,level,ALL_VECTORS,fas->l,x))
      return (1);

    for (i=0; i<fas->niter; i++) {
      if ((*fas->nliter->NLIter)
            (fas->nliter,fas->baselevel,x,fas->d,ass->A,fas->nlsolver.Assemble,&error)) {
        error = __LINE__;
        REP_ERR_RETURN (error);
      }
    }
    return(0);
  }
  /* keep value */
  if (dcopy(mg,level,level,ALL_VECTORS,fas->l,x))
    return (1);

  /* presmooth */
  for (i=0; i<fas->nu1; i++) {
    if ((*fas->nliter->NLIter)
          (fas->nliter,level,x,fas->d,ass->A,fas->nlsolver.Assemble,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }
  }

  /* restrict value */
  if (RestrictValue(g,x,x,fas->restriction))
    return (1);

  /* nonlinear defect on current level */
  if ((*ass->NLAssembleDefect)(ass,level,level,x,fas->d,ass->A,&error)) {
    error = __LINE__;
    REP_ERR_RETURN(error);
  }

  /* restrict defect */
  if (StandardRestrict (g,fas->d,fas->d,Factor_One))
    return (1);

  if (dcopy (mg,level-1,level-1,ALL_VECTORS,fas->v,x))
    return (1);

  /* recursiv call */
  for (i=0; i<fas->gamma; i++)
    if (FasStep(fas,ass,level-1,x))
      return (1);

  /* interpolate correction */
  /*if (dcopy(mg,level-1,level-1,ALL_VECTORS,fas->v,x))
      return(1);*/
  if (daxpyx (mg,level-1,level-1,ALL_VECTORS,fas->v,Factor_Minus_One,fas->l))
    return (1);
  if (StandardInterpolateCorrection(g,fas->v,fas->v,Factor_One))
    return(1);

  /* update solution */
  if (daxpyx  (mg,level,level,ALL_VECTORS,x,damp_factor,fas->v))
    return(1) ;

  /* postsmooth */
  for (i=0; i<fas->nu2; i++) {
    if ((*fas->nliter->NLIter)
          (fas->nliter,level,x,fas->d,ass->A,fas->nlsolver.Assemble,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }
  }

  return (0);
}

static INT FasSolverInit (NP_BASE *base, INT argc, char **argv)
{
  NP_FAS *fas;
  INT i,j;

  fas = (NP_FAS *) base;

  /* read  data descs */
  fas->l = ReadArgvVecDesc(base->mg,"l",argc,argv);
  fas->v = ReadArgvVecDesc(base->mg,"v",argc,argv);
  fas->d = ReadArgvVecDesc(base->mg,"d",argc,argv);

  /* read other numprocs */
  fas->trans = (NP_TRANSFER *) ReadArgvNumProc(base->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  if (fas->trans == NULL) {
    PrintErrorMessage('E',"FasSolverInit","cannot read transfer num proc");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  fas->nliter = (NP_NL_ITER *) ReadArgvNumProc(base->mg,"S",NL_ITER_CLASS_NAME,argc,argv);
  if (fas->nliter == NULL) {
    PrintErrorMessage('E',"FasSolverInit","cannot read iter num proc");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }

  /* set configuration parameters */
  if(sc_read(fas->damp,NP_FMT(fas),fas->l,"damp",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      fas->damp[i] = 1.0;
  if(sc_read(fas->restriction,NP_FMT(fas),fas->l,"res",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      fas->restriction[i] = 1.0;
  if (ReadArgvINT("maxit",&(fas->maxit),argc,argv))
    fas->maxit = 50;
  if ((fas->maxit<0)||(fas->maxit>100)) {
    PrintErrorMessage('E',"FasSolverInit","maxit <= 100");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("gamma",&(fas->gamma),argc,argv))
    fas->gamma = 1;
  if ((fas->gamma<0)||(fas->gamma>2)) {
    PrintErrorMessage('E',"FasSolverInit","gamma <= 2");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("baselevel",&(fas->baselevel),argc,argv))
    fas->baselevel = 0;
  if ((fas->baselevel<0)||(fas->baselevel>10)) {
    PrintErrorMessage('E',"FasSolverInit","baselevel <= 10");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }
  if (ReadArgvINT("n1",&(fas->nu1),argc,argv))
    fas->nu1 = 1;
  if (ReadArgvINT("n2",&(fas->nu2),argc,argv))
    fas->nu2 = 1;
  if (ReadArgvINT("niter",&(fas->niter),argc,argv))
    fas->niter = 1;

  /* set display option */
  fas->displayMode = ReadArgvDisplay(argc,argv);

  /* call general nls init */
  return (NPNLSolverInit(&(fas->nlsolver),argc,argv));
}

static INT FasSolverDisplay (NP_BASE *theNumProc)
{
  NP_FAS *fas;

  fas     = (NP_FAS*) theNumProc;

  /* general nls display */
  NPNLSolverDisplay(&(fas->nlsolver));

  if (fas->l != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"l",ENVITEM_NAME(fas->l));
  if (fas->v != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"v",ENVITEM_NAME(fas->v));
  if (fas->d != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"d",ENVITEM_NAME(fas->d));

  if (fas->nliter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(fas->nliter));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S","---");
  if (fas->trans != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(fas->trans));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");

  if (fas->displayMode == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (fas->displayMode == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (fas->displayMode == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  UserWriteF(DISPLAY_NP_FORMAT_SI,"maxit",(int)fas->maxit);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"gamma",(int)fas->gamma);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)fas->nu1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n2",(int)fas->nu2);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"niter",(int)fas->niter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)fas->baselevel);
  if (sc_disp(fas->damp,fas->l,"damp")) REP_ERR_RETURN (1);
  if (sc_disp(fas->restriction,fas->l,"res")) REP_ERR_RETURN (1);

  return (0);
}

static INT FasConstruct (NP_BASE *theNP)
{
  NP_NL_SOLVER *np;

  /* set general functions */
  theNP->Init = FasSolverInit;
  theNP->Display = FasSolverDisplay;
  theNP->Execute = NPNLSolverExecute;

  np = (NP_NL_SOLVER *) theNP;
  np->PreProcess = FasPreProcess;
  np->Solver = FasSolver;
  np->PostProcess = FasPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitFasSolver  - Init this file

   SYNOPSIS:
   INT InitFasSolver (void);

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

INT InitFasSolver (void)
{
  if (CreateClass (NL_SOLVER_CLASS_NAME ".fas",
                   sizeof(NP_FAS), FasConstruct))
    REP_ERR_RETURN (__LINE__);

  return (0);
}
