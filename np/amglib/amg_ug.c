// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amgsolver.c													*/
/*																			*/
/* Purpose:   interface ug<->sp sparse matrix package with algebraic		*/
/*			  multigrid.													*/
/*																			*/
/* Author:	  Peter Bastian					                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   05 FEB 1996 Begin    This is a very simple interface: It		*/
/*                                 allocates the problem completely new     */
/*								   when LMGC_Solve is called.				*/
/*            21 JAN 1997 new numprocs										*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#define __COMPILE_AMGSOLVER__

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "compiler.h"
#include "devices.h"
#include "gm.h"
#include "misc.h"
#include "evm.h"
#include "np.h"
#include "ls.h"
#include "ugstruct.h"
#include "commands.h"    /* for GetCurrentMultigrid          */
#include "cmdint.h"      /* for CreateCommand                */
#include "uginterface.h" /* for UserInterrupt                */
#include "scan.h"
#include "pcr.h"

#include "amg_header.h"
#include "amg_low.h"
#include "amg_sp.h"
#include "amg_blas.h"
#include "amg_iter.h"
#include "amg_coarsen.h"
#include "amg_solve.h"

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

typedef struct
{
  NP_LINEAR_SOLVER ls;

  VECDATA_DESC *c;
  INT display;

  INT scale;
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;

  AMG_MATRIX *A;
  AMG_VECTOR *x,*b;

} NP_AMG;

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

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   amgc - numproc for algebraic multigrid cycle

   DESCRIPTION:
   This numproc executes linear multigrid cycles.

   .vb
   npinit $x <sol sym> $b <rhs sym> $A <mat sym> $s {0|1}

   .ve

   .  $x~<sol~sym> - symbol for the solution vector
   .  $b~<rhs~sym> - symbol for the right hand side vector
   .  $A~<mat~sym> - symbol for the stiffness matrix
   .  $s~<mode> - 1: store always as scalar, 0 store as small block system


   'npexecute <name> '

   EXAMPLE:

   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   AMGSolverPreProcess - Prepare algebraic multigrid solver

   SYNOPSIS:
   static INT AMGSolverPreProcess (NP_LINEAR_SOLVER *theNP, INT level,
                                                                   VECDATA_DESC *VD_x, VECDATA_DESC *VD_b,
                                                                   MATDATA_DESC *MD_A,
                                                                   INT *baselevel, INT *result);

   PARAMETERS:
   see ls.c

   DESCRIPTION:
   This function prepares the algebraic multigrid cycle. It marks the heap
   and produces the coarse grid hierarchy.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static int mark_counter=0;

static MULTIGRID *amgMG;

static void *amgmalloc (size_t n)
{
  return(GetMem(MGHEAP(amgMG),n,FROM_BOTTOM));
}

static INT AMGSolverPreProcess (NP_LINEAR_SOLVER *theNP, INT level,
                                VECDATA_DESC *VD_x, VECDATA_DESC *VD_b,
                                MATDATA_DESC *MD_A,
                                INT *baselevel, INT *result)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  void *buffer;
  int n,nonzeros,blocksize;
  int Acomp;
  MATRIX *theMatrix,*theNeighbor;
  VECTOR *theVector;
  int index;
  int rv,i,j,k,block_i,block_j;
  char buf[128];
  int nRows_A,nCols_A,nComp_x,nComp_b;
  NP_AMG *theAMGC;
  double ti;
  int ii;
        #ifdef ModelP
  double clock_start;
        #else
  clock_t clock_start;
        #endif

  theAMGC = (NP_AMG *) theNP;

  /* prepare solving */
  theMG = theAMGC->ls.base.mg;
  theGrid = GRID_ON_LEVEL(theMG,level);

  /* mark heap for use by amg */
  Mark(MGHEAP(theMG),FROM_BOTTOM);
  mark_counter++;

  /* initialize sp package */
  AMG_InstallPrintHandler((AMG_PrintFuncPtr)UserWrite);
  amgMG=theMG;       /* make it global for memory handler */
  AMG_InstallMallocHandler((AMG_MallocFuncPtr)amgmalloc);

  /* get access to components */
  nRows_A = MD_ROWS_IN_RT_CT(MD_A,NODEVEC,NODEVEC);
  nCols_A = MD_COLS_IN_RT_CT(MD_A,NODEVEC,NODEVEC);
  nComp_x = VD_NCMPS_IN_TYPE(VD_x,NODEVEC);
  nComp_b = VD_NCMPS_IN_TYPE(VD_b,NODEVEC);
  blocksize = nComp_x;
  if (blocksize==0) goto exit;
  if (nComp_b!=blocksize) goto exit;
  if (nCols_A!=blocksize) goto exit;
  if (nRows_A!=blocksize) goto exit;
  Acomp = MD_MCMP_OF_RT_CT(MD_A,NODEVEC,NODEVEC,0);

  CSTART(); ti=0; ii=0;
  /* diagonal scaling */
  if (theAMGC->scale)
    if (DiagonalScaleSystem(theGrid,MD_A,MD_A,VD_b)!=NUM_OK)
    {
      UserWrite("Error in scaling system\n");
      goto exit;
    }

  /* gather some data for the matrix */
  n = nonzeros = 0;

  /* loop through all vectors, we assume there are only node vectors ! */
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
  {
    VINDEX(theVector) = n++;             /* renumber vectors just to be sure ... */
    /* now speed through this row */
    for (theMatrix=VSTART(theVector); theMatrix!=NULL; theMatrix = MNEXT(theMatrix))
      nonzeros++;
  }

  /* now allocate fine grid vectors x and b */
  theAMGC->x = AMG_NewVector(n*blocksize,1,"x");
  if (theAMGC->x==NULL) {
    UserWrite("no memory for x\n");
    goto exit;
  }
  theAMGC->b = AMG_NewVector(n*blocksize,1,"b");
  if (theAMGC->b==NULL) {
    UserWrite("no memory for b\n");
    goto exit;
  }

  /* and a new matrix */
  theAMGC->A = AMG_NewMatrix(n*blocksize,1,nonzeros*blocksize*blocksize,blocksize,"fine grid A");
  if (theAMGC->A==NULL) {
    UserWrite("no memory for A\n");
    goto exit;
  }

  /* now fill matrix */
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
  {
    i = VINDEX(theVector);

    /* count row length */
    nonzeros=0;
    for (theMatrix=VSTART(theVector); theMatrix!=NULL; theMatrix = MNEXT(theMatrix))
      nonzeros++;

    /* for each row */
    for (block_i=0; block_i<blocksize; block_i++)
    {
      /* allocate row */
      if (AMG_SetRowLength(theAMGC->A,i*blocksize+block_i,nonzeros*blocksize)!=AMG_OK)
      {
        UserWrite("Error in AMG_SetRowLength\n");
        goto exit;
      }

      /* the diagonal block, be careful to allocate the main diagonal first */
      theMatrix=VSTART(theVector);
      if (AMG_InsertValues(theAMGC->A,i*blocksize+block_i,i*blocksize+block_i,
                           &(MVALUE(theMatrix,Acomp+block_i*blocksize+block_i)))<0)
      {
        UserWrite("Error in AMG_InsertValues\n");
        goto exit;
      }
      for (block_j=0; block_j<blocksize; block_j++)
      {
        if (block_j==block_i) continue;
        if (AMG_InsertValues(theAMGC->A,i*blocksize+block_i,i*blocksize+block_j,
                             &(MVALUE(theMatrix,Acomp+block_i*blocksize+block_j)))<0)
        {
          UserWrite("Error in AMG_InsertValues\n");
          goto exit;
        }
      }

      /* all the offdiagonal blocks */
      for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix = MNEXT(theMatrix))
      {
        j = VINDEX(MDEST(theMatrix));
        for (block_j=0; block_j<blocksize; block_j++)
        {
          if (AMG_InsertValues(theAMGC->A,i*blocksize+block_i,j*blocksize+block_j,
                               &(MVALUE(theMatrix,Acomp+block_i*blocksize+block_j)))<0)
          {
            UserWrite("Error in AMG_InsertValues\n");
            goto exit;
          }
        }
      }
    }
  }
  /*AMG_PrintMatrix(theAMGC->A,"Matrix");*/

  /* call algebraic multigrid solver */
  if (AMG_Build(&theAMGC->sc,&theAMGC->cc,theAMGC->A)!=AMG_OK) goto exit;
  CSTOP(ti,ii);
  UserWriteF("AMG : L=%2d BUILD=%10.4lg\n",level,ti);

  return(0);       /* ok, matrix is set up */

exit: /* error */
  if (mark_counter>0) {
    Release(MGHEAP(theMG),FROM_BOTTOM);
    mark_counter--;
  }
  return(1);
}

static INT AMGSolverDefect (NP_LINEAR_SOLVER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *result)
{
  if (dmatmul_minus(theNP->base.mg,0,level,ON_SURFACE,b,A,x)
      != NUM_OK) {
    result[0] = __LINE__;
    return(1);
  }
  return(0);
}

static INT AMGSolverResiduum (NP_LINEAR_SOLVER *theNP, INT bl, INT level,
                              VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                              LRESULT *lresult)
{
        #ifdef ModelP
  if (a_vector_collect(theNP->base.mg,bl,level,b)) {
    lresult->error_code = __LINE__;
    return(1);
  }
        #endif
  if (dnrm2x(theNP->base.mg,bl,level,ON_SURFACE,b,lresult->last_defect)) {
    lresult->error_code = __LINE__;
    return(1);
  }

  return(0);
}

/********************************************************/
/*D
   AMGSolver - algebraic multigrid solver

   SYNOPSIS:
   static INT AMGSolver (NP_LINEAR_SOLVER *theNP, INT level,
                                                 VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                                                 VEC_SCALAR abslimit, VEC_SCALAR reduction,
                                                 LRESULT *lresult)

   PARAMETERS:
   see ls.c

   DESCRIPTION:
   copies right hand side to amg package and calls the solver.

   RETURN VALUE:
   INT
   .n    0 if converged
   .n    1 if error occured
   and structure lresult is filled
   D*/
/****************************************************/

static INT AMGSolver (NP_LINEAR_SOLVER *theNP, INT level,
                      VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                      VEC_SCALAR abslimit, VEC_SCALAR reduction,
                      LRESULT *lresult)
{
  NP_AMG *theAMGC;
  VEC_SCALAR defect2reach;
  INT rv,i,bl,PrintID;
  char text[DISPLAY_WIDTH+4];
  VEC_SCALAR Factor_One;
  MULTIGRID *theMG;
  GRID *theGrid;
  INT converged=0;
  int xcomp,bcomp;
  VECTOR *theVector;
  int k;
  int nComp_x,nComp_b,blocksize;
  double ti;
  int ii;
        #ifdef ModelP
  double clock_start;
        #else
  clock_t clock_start;
        #endif

  /* prepare solving */
  theAMGC = (NP_AMG *) theNP;
  theMG = theAMGC->ls.base.mg;
  theGrid = GRID_ON_LEVEL(theMG,level);
  theAMGC->sc.red_factor=reduction[0];
  theAMGC->sc.dnorm_min=abslimit[0];

  bl = 0;
  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  /* allocate correction */
  if (AllocVDFromVD(theNP->base.mg,0,level,x,&theAMGC->c)) {
    lresult->error_code = __LINE__;
    return(1);
  }

  /* print defect */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(theAMGC),'*',"\n");
  if (PreparePCR(x,theAMGC->display,text,&PrintID)) {
    lresult->error_code = __LINE__;
    return(1);
  }
  for (i=0; i<VD_NCOMP(x); i++)
    lresult->first_defect[i] = lresult->last_defect[i];
  if (sc_mul_check(defect2reach,lresult->first_defect,reduction,b)) {
    lresult->error_code = __LINE__;
    return(1);
  }
  if (DoPCR(PrintID,lresult->first_defect,PCR_CRATE)) {
    lresult->error_code = __LINE__;
    return(1);
  }
  if (sc_cmp(lresult->first_defect,abslimit,b)) lresult->converged = 1;
  else lresult->converged = 0;

  CSTART(); ti=0; ii=0;

  /* fill values in x,b */
  xcomp = VD_ncmp_cmpptr_of_otype(theAMGC->c,NODEVEC,&nComp_x)[0];
  bcomp = VD_ncmp_cmpptr_of_otype(b,NODEVEC,&nComp_b)[0];
  blocksize = nComp_x;
  if (blocksize==0) goto exit;
  if (nComp_b!=blocksize) goto exit;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
  {
    i = VINDEX(theVector);
    for (k=0; k<blocksize; k++)
    {
      AMG_VECTOR_ENTRY(theAMGC->b,i*blocksize+k,0)=VVALUE(theVector,bcomp+k);
    }
  }
  AMG_dset(theAMGC->x,0.0);
  if ((rv=AMG_Solve(theAMGC->x,theAMGC->b))<0)
  {
    lresult->error_code = __LINE__;
    lresult->converged = 0;
    goto exit;
  }
  lresult->converged = 1;
  lresult->number_of_linear_iterations = rv;

  /* write back solution values */
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
  {
    i = VINDEX(theVector);
    for (k=0; k<blocksize; k++)
      VVALUE(theVector,xcomp+k)=AMG_VECTOR_ENTRY(theAMGC->x,i*blocksize+k,0);
  }

  if (dmatmul_minus(theNP->base.mg,0,level,ON_SURFACE,b,A,theAMGC->c)
      != NUM_OK) {
    lresult->error_code = __LINE__;
    return(1);
  }
  if (daxpyx(theAMGC->ls.base.mg,0,level,
             ON_SURFACE,x,Factor_One,theAMGC->c) != NUM_OK) {
    lresult->error_code = __LINE__;
    return(1);
  }
  if (AMGSolverResiduum(theNP,bl,level,x,b,A,lresult))
    return(1);
  if (DoPCR(PrintID, lresult->last_defect,PCR_CRATE)) {
    lresult->error_code = __LINE__;
    return (1);
  }

  if (DoPCR(PrintID,lresult->last_defect,PCR_AVERAGE)) {
    lresult->error_code = __LINE__;
    return (1);
  }
  FreeVD(theNP->base.mg,0,level,theAMGC->c);
  if (PostPCR(PrintID,NULL)) {
    lresult->error_code = __LINE__;
    return (1);
  }
  CSTOP(ti,ii);
  if (lresult->number_of_linear_iterations != 0)
    UserWriteF("AMG : L=%2d N=%2d TSOLVE=%10.4lg TIT=%10.4lg\n",level,
               lresult->number_of_linear_iterations,ti,
               ti/lresult->number_of_linear_iterations);
  else
    UserWriteF("AMG : L=%2d N=%2d TSOLVE=%10.4lg\n",level,
               lresult->number_of_linear_iterations,ti);

  return (0);

exit:
  return(1);
}

static INT AMGSolverPostProcess (NP_LINEAR_SOLVER *theNP,
                                 INT level,
                                 VECDATA_DESC *x, VECDATA_DESC *b,
                                 MATDATA_DESC *A,
                                 INT *result)
{
  if (mark_counter>0) {
    Release(MGHEAP(theNP->base.mg),FROM_BOTTOM);
    mark_counter--;
  }

  return(0);
}

/****************************************************************************/
/*
   AMGC_Init - Init algebraic multigrid

   SYNOPSIS:
   static INT AMGSolverInit (NP_BASE *theNP, INT argc , char **argv);

   PARAMETERS:
   .  theNP - pointer to numproc
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function sets all parameters of the algebraic multigrid method. All
   options are explained in detail below.

   CURRENT FUNCTIONALITY:
   AMG is implemented for scalar equations and with restrictions for
   systems. (1) Unknowns may be only in the nodes of the mesh (this can
   be extended easily e.g. for edge unknowns, only the amg-ug interface
   is affected. Systems are filled into AMG as a sparse matrix with individual
   entries. AMG makes use of the system character by assuming that row number
   modulo blocksize gives the system component. Example: If the system has
   components u,v,p then the row of sparse matrix A correspond to u1,v1,p1,u2,v2,p2,..
   and so on.

   HOW IT WORKS:
   This amg is an agglomeration multigrid variant. It uses the usual strong/weak
   coupling definition. Then the set of unknowns is partitioned into clusters
   (subsets) according to a number of rules (most important: unknowns in
   a cluster are strongly coupled). In the case of systems only unknowns of the
   same component are clustered. The error to be computed by the multigrid
   coarse grid correction is now assumed to be constant in each cluster. The coarse
   grid operator is obtained by Galerkin which is very simple because prolongation
   and restriction is constant on each cluster.

   AMG uses its own compressed row storage data structure for sparse matrices,
   it can be used stand alone without ug (see test.c in the amglib directory)
   and has also cg and bicgstab krylov solvers.

   The interface to ug is realized as a NP_LINEAR_SOLVER numproc with clas name "amg".
   It can be used on any level that is UNIFORMELY REFINED (!).
   It can also be used as a base solver in geometric multigrid.

   OPTIONS:
   .  $A <sym> - MATDATA_DESC for the matrix as in all linear solvers, only necessary if executed.
   .  $x <sym> - VECDATA_DESC for the solution.
   .  $b <sym> - VECDATA_DESC for right hand side.
   .  $abslimit <VEC_SCALAR> -  for absolute limit.
   .  $red <VEC_SCALAR> - for reduction factor.
   .  $display full|red|no - affects only first and last defect.

   .  $vc <value> - verbose value for coarsening, default is 1, higher means more output.
   .  $dependency <string> - use 'sym' for symmetric problems, 'unsym' for others, default is unsym.
   .  $alpha <value> - threshold for strong coupling default is 0.4
   .  $beta  <value> - threshold for diagonal dominant (isolated) nodes, default is 0.001
   .  $minc <value> - smallest desired cluster size, default is 4 in 2D, 8 in 3D.
   .  $maxc <value> - largest possible cluster size, default is 6 in 2D, 10 in 3D.
   .  $maxd <value> - maximum diameter of a cluster in strong connections, default 2.
   .  $maxcon <value> - stops clusters from having many neighbors, default is 30.
   .  $major <value> - clustering depends on one component only if >=0, default -1.
   .  $dt <value> - maximum depth allowed, n means n coarsening steps are allowed, default 20.
   .  $ct <value> - coarsen target, stop coarsening if number of unknowns is reached, default 10.
   .  $cr <value> - coarsen rate, stop coarsening if it becomes too slow, default 1.5.
   .  $scale 0|1 - scale system by point block diagonal before loading matrix to AMG.
                This should be used for systems, but NOT if AMG is used as basesolver in lmgc.
                (Then use the $S option in transfer if your system becomes zero on diagonal)

   .  $vs <value> - verbose level for solver, higher means more aoutput, default is 1.
   .  $solver <string> - select solver:  bcgs | cg | ls, default bcgs.
   .  $prec <string> - select preconditioner: ssor | sor | jac | mgc, default mgc.
   .  $it <value> - maximum number of iterations allowed, default 80.
   .  $sm <string> - select smoother for multigrid: ssor | sor | jac, default ssor.
   .  $oms <value> - damping factor for smoother, default 1.0
   .  $n1 <value> -  number of presmoothing steps in mg, default 2.
   .  $n2 <value> -  number of postsmoothing steps in mg, default 2.
   .  $g <value> -  gamma in mg, default 1 (V-cycle).
   .  $omp <value> - damping factor for prolongation, default 1.8, adapts automatically.
   .  $csm <string> - select smoother for coarse grid solver: ssor | sor | jac | ex default ssor.
   .  $cit <value> - maximum number of iterations on coarsest grid, default 100.
   .  $cred <value> - reduction factor on coarsest grid, default 0.001.

   EXAMPLE:
   npcreate amgs $c amg;
   npinit amgs $red 1.0E-8 $abslimit 1.0E-15 $display full
      $alpha 0.4 $beta 1.0E-3 $minc 4 $maxc 6 $maxd 2
      $maxcon 30 $vc 1 $dt 20 $ct 10 $cr 1.3
      $scale 1 $vs 1 $solver bcgs $prec mgc $it 80 $sm ssor $n1 2 $n2 2 $g 2 $omp 1.8 $oms 1.0
      $csm ssor $cit 100 $cred 1.0E-3;


   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT AMGSolverInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_AMG *theAMGC;
  INT i,ret;
  double w;
  char buf[128];

  theAMGC = (NP_AMG *)theNP;

  ret = NPLinearSolverInit(&theAMGC->ls,argc,argv);       /* standard linear solver init */

  theAMGC->display = ReadArgvDisplay(argc,argv);


  /* fill coarsen context */
  if (!ReadArgvDOUBLE("alpha",&w,argc,argv))
    theAMGC->cc.alpha = w;
  else
    theAMGC->cc.alpha = 0.4;
  if (!ReadArgvDOUBLE("beta",&w,argc,argv))
    theAMGC->cc.beta = w;
  else
    theAMGC->cc.beta = 1.0E-3;
  if (!ReadArgvINT("minc",&i,argc,argv))
    theAMGC->cc.mincluster=i;
  else
  {
    if (DIM==2) theAMGC->cc.mincluster=4;
    if (DIM==3) theAMGC->cc.mincluster=8;
  }
  if (!ReadArgvINT("maxc",&i,argc,argv))
    theAMGC->cc.maxcluster=i;
  else
  {
    if (DIM==2) theAMGC->cc.maxcluster=6;
    if (DIM==3) theAMGC->cc.maxcluster=10;
  }
  if (!ReadArgvINT("maxd",&i,argc,argv))
    theAMGC->cc.maxdistance=i;
  else
    theAMGC->cc.maxdistance=2;
  if (!ReadArgvINT("maxcon",&i,argc,argv))
    theAMGC->cc.maxconnectivity=i;
  else
    theAMGC->cc.maxconnectivity=30;
  if (!ReadArgvINT("vc",&i,argc,argv))
    theAMGC->cc.verbose=i;
  else
    theAMGC->cc.verbose=1;
  if (!ReadArgvINT("dt",&i,argc,argv))
    theAMGC->cc.depthtarget=i;
  else
    theAMGC->cc.depthtarget=20;
  if (!ReadArgvINT("ct",&i,argc,argv))
    theAMGC->cc.coarsentarget=i;
  else
    theAMGC->cc.coarsentarget=10;
  if (!ReadArgvDOUBLE("cr",&w,argc,argv))
    theAMGC->cc.coarsenrate=w;
  else
    theAMGC->cc.coarsenrate=1.5;
  if (!ReadArgvINT("major",&i,argc,argv))
    theAMGC->cc.major=i;
  else
    theAMGC->cc.major=-1;
  if (!ReadArgvChar("dependency",buf,argc,argv))
  {
    if (strcmp(buf,"sym")==0) theAMGC->cc.dependency=AMG_SYM;
    if (strcmp(buf,"unsym")==0) theAMGC->cc.dependency=AMG_UNSYM;
  }
  else
    theAMGC->cc.dependency=AMG_UNSYM;

  /* fill solver context */
  if (!ReadArgvINT("vs",&i,argc,argv))
    theAMGC->sc.verbose=i;
  else
    theAMGC->sc.verbose=1;
  if (!ReadArgvChar("solver",buf,argc,argv))
  {
    if (strcmp(buf,"bcgs")==0) theAMGC->sc.solver=AMG_BCGS;
    if (strcmp(buf,"cg")==0) theAMGC->sc.solver=AMG_CG;
    if (strcmp(buf,"ls")==0) theAMGC->sc.solver=AMG_LS;
  }
  else
    theAMGC->sc.solver=AMG_BCGS;
  if (!ReadArgvChar("prec",buf,argc,argv))
  {
    if (strcmp(buf,"ssor")==0) theAMGC->sc.preconditioner=AMG_SSOR;
    if (strcmp(buf,"sor")==0) theAMGC->sc.preconditioner=AMG_SOR;
    if (strcmp(buf,"jac")==0) theAMGC->sc.preconditioner=AMG_DJAC;
    if (strcmp(buf,"mgc")==0) theAMGC->sc.preconditioner=AMG_MGC;
  }
  else
    theAMGC->sc.preconditioner=AMG_MGC;
  if (!ReadArgvINT("it",&i,argc,argv))
    theAMGC->sc.maxit=i;
  else
    theAMGC->sc.maxit=80;
  theAMGC->sc.red_factor=theAMGC->ls.reduction[0];
  theAMGC->sc.dnorm_min=theAMGC->ls.abslimit[0];
  if (!ReadArgvChar("csm",buf,argc,argv))
  {
    if (strcmp(buf,"ssor")==0) theAMGC->sc.coarse_smoother=AMG_SSOR;
    if (strcmp(buf,"sor")==0) theAMGC->sc.coarse_smoother=AMG_SOR;
    if (strcmp(buf,"jac")==0) theAMGC->sc.coarse_smoother=AMG_DJAC;
    if (strcmp(buf,"ex")==0) theAMGC->sc.coarse_smoother=AMG_EX;
  }
  else
    theAMGC->sc.coarse_smoother=AMG_SSOR;
  if (!ReadArgvINT("cit",&i,argc,argv))
    theAMGC->sc.coarse_maxit=i;
  else
    theAMGC->sc.coarse_maxit=100;
  if (!ReadArgvDOUBLE("cred",&w,argc,argv))
    theAMGC->sc.coarse_red_factor=w;
  else
    theAMGC->sc.coarse_red_factor=1.0E-3;
  if (!ReadArgvINT("n1",&i,argc,argv))
    theAMGC->sc.n1=i;
  else
    theAMGC->sc.n1=2;
  if (!ReadArgvINT("n2",&i,argc,argv))
    theAMGC->sc.n2=i;
  else
    theAMGC->sc.n2=2;
  if (!ReadArgvINT("g",&i,argc,argv))
    theAMGC->sc.gamma=i;
  else
    theAMGC->sc.gamma=1;
  if (!ReadArgvChar("sm",buf,argc,argv))
  {
    if (strcmp(buf,"ssor")==0) theAMGC->sc.smoother=AMG_SSOR;
    if (strcmp(buf,"sor")==0) theAMGC->sc.smoother=AMG_SOR;
    if (strcmp(buf,"jac")==0) theAMGC->sc.smoother=AMG_DJAC;
  }
  else
    theAMGC->sc.smoother=AMG_SSOR;
  if (!ReadArgvDOUBLE("omp",&w,argc,argv))
  {
    for (i=0; i<AMG_MAX_COMP; i++)
      theAMGC->sc.omega_p[i]=w;
  }
  else
  {
    for (i=0; i<AMG_MAX_COMP; i++)
      theAMGC->sc.omega_p[i]=1.8;
  }
  if (!ReadArgvDOUBLE("oms",&w,argc,argv))
  {
    for (i=0; i<AMG_MAX_COMP; i++)
      theAMGC->sc.omega[i]=w;
  }
  else
  {
    for (i=0; i<AMG_MAX_COMP; i++)
      theAMGC->sc.omega[i]=1.8;
  }

  if (!ReadArgvINT("scale",&i,argc,argv))
    theAMGC->scale=i;
  else
    theAMGC->scale=0;

  return (ret);
}

/****************************************************************************/
/*
   AMGC_Display	- Display linear multigrid cycle

   SYNOPSIS:
   static INT AMGC_Display (NUM_PROC *theNumProc);

   PARAMETERS:
   .  theNumProc - pointer to numproc

   DESCRIPTION:
   This function displays linear multigrid cycle.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT AMGSolverDisplay (NP_BASE *theNP)
{
  NP_AMG *theAMGC;
  char buffer[128];

  theAMGC = (NP_AMG *)theNP;

  /* display symbols */
  NPLinearSolverDisplay(&theAMGC->ls);

  /* display configuration parameters */
  UserWrite("configuration parameters:\n");

  return (0);
}

/****************************************************************************/
/*
   AMGC_Execute	- Execute linear multigrid cycle

   SYNOPSIS:
   static INT AMGC_Execute (NUM_PROC *theNumProc, MULTIGRID *theMG,
   INT argc, char **argv);

   PARAMETERS:
   .  theNumProc - pointer to numproc
   .  theMG - pointer to multigrid
   .  argc - argument counter
   .  argv -argument vector

   DESCRIPTION:
   This function executes linear multigrid cycle.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT AMGSolverExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LINEAR_SOLVER *np;
  VECDATA_DESC *x,*b;
  MATDATA_DESC *A;
  LRESULT lresult;
  INT result,level,bl=0;

  np = (NP_LINEAR_SOLVER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->x == NULL) {
    PrintErrorMessage('E',"AMGSolverExecute","no vector x");
    return (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"AMGSolverExecute","no vector b");
    return (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"AMGSolverExecute","no matrix A");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"AMGSolverExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,np->x,np->b,np->A,&bl,&result)) {
      UserWriteF("AMGSolverExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("d",argc,argv)) {
    if (np->Defect == NULL) {
      PrintErrorMessage('E',"AMGSolverExecute","no Defect");
      return (1);
    }
    if ((*np->Defect)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("AMGSolverExecute: Defect failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("r",argc,argv)) {
    if (np->Residuum == NULL) {
      PrintErrorMessage('E',"AMGSolverExecute","no Residuum");
      return (1);
    }
    if ((*np->Residuum)(np,bl,level,np->x,np->b,np->A,&lresult)) {
      UserWriteF("AMGSolverExecute: Residuum failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->Solver == NULL) {
      PrintErrorMessage('E',"AMGSolverExecute","no Solver");
      return (1);
    }
    if ((*np->Solver)(np,level,np->x,np->b,np->A,
                      np->abslimit,np->reduction,&lresult)) {
      UserWriteF("AMGSolverExecute: Solver failed, error code %d\n",
                 lresult.error_code);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"AMGSolverExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->x,np->b,np->A,&result)) {
      UserWriteF("AMGSolverExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }
  return(0);
}

/****************************************************************************/
/*D
   InitMgSolver - Enrol mgsolvers

   SYNOPSIS:
   INT InitMgSolver (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function enrols mgsolvers: its creates the numproc 'AMGC'
   and sets global variables in 'mgsolver.c'.
   It is called in InitNumerics.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT AMGConstruct (NP_BASE *theNP)
{
  NP_AMG *np;

  theNP->Init                     = AMGSolverInit;
  theNP->Display          = AMGSolverDisplay;
  theNP->Execute          = AMGSolverExecute;

  np = (NP_AMG *) theNP;
  np->ls.PreProcess       = AMGSolverPreProcess;
  np->ls.Defect           = AMGSolverDefect;
  np->ls.Residuum         = AMGSolverResiduum;
  np->ls.Solver           = AMGSolver;
  np->ls.PostProcess      = AMGSolverPostProcess;

  return(0);
}

INT InitAMGSolver (void)
{
  if (CreateClass(LINEAR_SOLVER_CLASS_NAME ".amg",sizeof(NP_AMG),AMGConstruct))
    return (__LINE__);
  return (0);
}
