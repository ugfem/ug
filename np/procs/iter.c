// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  iter.c			                                                                                */
/*																			*/
/* Purpose:   iteration num procs                                               */
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
#include "gm.h"
#include "algebra.h"
#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "devices.h"
#include "udm.h"

#include "transfer.h"
#include "ls.h"
#include "iter.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* macros for the symmetric Gauss-Seidel smoother */
#define NP_SGS_t(p)             ((p)->t)

/* macros for the Block Gauss-Seidel smoother */
#define MAX_BLOCKS                  3
#define MAX_ORDER                   6
#define NP_BGS_NBLOCKS(p)             ((p)->nBlocks)
#define NP_BGS_BLOCKS(p)              ((p)->Block)
#define NP_BGS_BLOCKSTART(p,i)        ((p)->Block[i])
#define NP_BGS_BLOCKEND(p,i)          ((p)->Block[(i)+1])
#define NP_BGS_BLOCKITERS(p)          ((p)->BlockIter)
#define NP_BGS_BLOCKITER(p,i)         ((p)->BlockIter[i])
#define NP_BGS_BLOCKITNAME(p,i)       ENVITEM_NAME(NP_BGS_BLOCKITER(p,i))
#define NP_BGS_NBLOCKITER(p)          ((p)->nBlockIter)
#define NP_BGS_BLOCKORDER(p)          ((p)->BlockOrder)
#define NP_BGS_BLOCKORD(p,i)          ((p)->BlockOrder[i])
#define NP_BGS_MD_Ab(p,cb)            (&((p)->MD_Ab[cb]))
#define NP_BGS_VD_rb(p,b)             (&((p)->VD_rb[b]))
#define NP_BGS_VD_tb(p,b)             (&((p)->VD_tb[b]))
#define NP_BGS_COMPS_A(p)                       ((p)->COMP_A)
#define NP_BGS_COMPS_r(p)                       ((p)->COMP_r)

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct np_smoother {

  NP_ITER iter;

  VEC_SCALAR damp;
  MATDATA_DESC *L;

  INT (*Step)
    (struct np_smoother *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* correction vector               */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    MATDATA_DESC *,                              /* temporary matrix                */
    INT *);                                      /* result                          */
};
typedef struct np_smoother NP_SMOOTHER;

typedef struct
{
  NP_SMOOTHER smoother;

  VECDATA_DESC *t;

} NP_SGS;

typedef struct
{
  NP_SMOOTHER smoother;

  INT nBlocks;                      /* number of blocks                     */
  INT Block[MAX_BLOCKS+1];          /* block subdivision                    */
  NP_SMOOTHER *BlockIter[MAX_BLOCKS];  /* block iteration scheme            */
  INT nBlockIter;                   /* number of block iterations           */
  INT BlockOrder[MAX_ORDER];        /* iteration order for the blocks       */
  MATDATA_DESC MD_Ab[MAX_BLOCKS];   /* blocked lower stiffness matrix       */
  VECDATA_DESC VD_rb[MAX_BLOCKS];   /* blocked right hand side              */
  VECDATA_DESC VD_tb[MAX_BLOCKS];   /* blocked temporary                    */

  /* storage for components */
  SHORT COMP_A[MAX_BLOCKS*MAX_MAT_COMP];
  SHORT COMP_r[MAX_BLOCKS*MAX_VEC_COMP];

} NP_BGS;

typedef struct
{
  NP_SMOOTHER smoother;

  VEC_SCALAR beta;

} NP_ILU;

typedef struct
{
  NP_SMOOTHER smoother;

  VEC_SCALAR beta;
  INT mode;

} NP_SPILU;

typedef struct
{
  NP_SMOOTHER smoother;

  VEC_SCALAR beta;
  VEC_SCALAR thresh;

} NP_THILU;

typedef struct
{
  NP_ITER iter;

  INT gamma;
  INT nu1;
  INT nu2;
  INT baselevel;

  NP_TRANSFER *Transfer;
  NP_ITER *PreSmooth;
  NP_ITER *PostSmooth;
  NP_LINEAR_SOLVER *BaseSolver;

  VECDATA_DESC *t;

} NP_LMGC;

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
   NP_ITER - type definition for iterations

   DESCRIPTION:
   This numproc type is used for the description of linear iterations.
   It can be called by the given interface from a linear solver.
   Initializing the data is optional; it can be done with

   'INT NPIterInit (NP_ITER *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPIterDisplay (NP_ITER *theNP);'
   'INT NPIterExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb


   ..... fill in data structure here when the realizition is finished


   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

INT NPIterInit (NP_ITER *np, INT argc , char **argv)
{
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->c = ReadArgvVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"r",argc,argv);

  if ((np->A == NULL) || (np->b == NULL) || (np->c == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPIterDisplay (NP_ITER *np)
{
  if ((np->A == NULL) && (np->b == NULL) && (np->c == NULL))
    return(0);
  UserWrite("symbolic user data:\n");
  if (np->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  if (np->b != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->b));
  if (np->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(np->c));
  UserWrite("\n");

  return(0);
}

INT NPIterExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ITER *np;
  INT result,bl,level;

  np = (NP_ITER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->c == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no vector c");
    return (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no vector b");
    return (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no matrix A");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no PreProcess");
      return (1);
    }
    if ((*np->PreProcess)(np,level,np->c,np->b,np->A,&bl,&result)) {
      UserWriteF("NPIterExecute: PreProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->Iter == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no Iter");
      return (1);
    }
    if ((*np->Iter)(np,level,np->c,np->b,np->A,&result)) {
      UserWriteF("NPIterExecute: Iter failed, error code %d\n", result);
      return (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no PostProcess");
      return (1);
    }
    if ((*np->PostProcess)(np,level,np->c,np->b,np->A,&result)) {
      UserWriteF("NPIterExecute: PostProcess failed, error code %d\n",
                 result);
      return (1);
    }
  }

  return(0);
}

/* tools for all smoothers */

static INT SmootherInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SMOOTHER *np;
  INT i;

  np = (NP_SMOOTHER *) theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i] = 1.0;
  sc_read(np->damp,np->iter.b,"damp",argc,argv);
  np->L = ReadArgvMatDesc(theNP->mg,"L",argc,argv);

  return (NPIterInit(&np->iter,argc,argv));
}

static INT SmootherDisplay (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  np = (NP_SMOOTHER *) theNP;
  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) return (1);
  if (np->L != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(np->L));

  return (0);
}

static INT Smoother (NP_ITER *theNP, INT level,
                     VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                     INT *result)
{
  NP_SMOOTHER *np;
  GRID *theGrid;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  np = (NP_SMOOTHER *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if ((*np->Step)(np,level,x,b,A,np->L,result))
    return (1);
    #ifdef ModelP
  if (l_vector_consistent(theGrid,x) != NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_dscale(theGrid,x,ACTIVE_CLASS,np->damp) != NUM_OK) NP_RETURN(1,result[0]);
  if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,x,ACTIVE_CLASS)
      != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SmootherPostProcess (NP_ITER *theNP, INT level,
                                VECDATA_DESC *x, VECDATA_DESC *b,
                                MATDATA_DESC *A, INT *result)
{
  NP_SMOOTHER *np;

  np = (NP_SMOOTHER *) theNP;
  if (np->L != NULL)
    FreeMD(theNP->base.mg,level,level,np->L);

  return(0);
}

/****************************************************************************/
/*D
   jac - numproc for Jacobi smoother

   DESCRIPTION:
   This numproc executes a block Jacobi smoother, using the blas routine
   'l_jac'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT JacobiPreProcess  (NP_ITER *theNP, INT level,
                              VECDATA_DESC *x, VECDATA_DESC *b,
                              MATDATA_DESC *A, INT *baselevel, INT *result)
{
        #ifdef ModelP
  NP_SMOOTHER *np;
  GRID *theGrid;

  np = (NP_SMOOTHER *) theNP;
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->L)) NP_RETURN(1,result[0]);
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_dmatcopy(theGrid,np->L,A) != NUM_OK) NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->L,MAT_DIAG_CONS) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  *baselevel = level;

  return (0);
}


static INT JacobiStep (NP_SMOOTHER *theNP, INT level,
                       VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                       MATDATA_DESC *L,
                       INT *result)
{
    #ifdef ModelP
  if (l_jac(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (l_jac(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,A,b) != NUM_OK) NP_RETURN(1,result[0]);
    #endif

  return (0);
}

static INT JacobiConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SmootherInit;
  theNP->Display = SmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = JacobiPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = JacobiStep;

  return(0);
}

/****************************************************************************/
/*D
   gs - numproc for Gauss-Seidel smoother

   DESCRIPTION:
   This numproc executes a Gauss-Seidel smoother, using the blas routines
   'l_lgs'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT GSPreProcess  (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b,
                          MATDATA_DESC *A, INT *baselevel, INT *result)
{
        #ifdef ModelP
  NP_SMOOTHER *np;
  GRID *theGrid;

  np = (NP_SMOOTHER *) theNP;
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->L)) NP_RETURN(1,result[0]);
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_dmatcopy(theGrid,np->L,A) != NUM_OK) NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->L,MAT_MASTER_CONS) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  *baselevel = level;

  return (0);
}


static INT GSStep (NP_SMOOTHER *theNP, INT level,
                   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                   MATDATA_DESC *L,
                   INT *result)
{
    #ifdef ModelP
  if (l_vector_collect(GRID_ON_LEVEL(theNP->iter.base.mg,level),b)!=NUM_OK) NP_RETURN(1,result[0]);

  if (l_lgs(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (l_lgs(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,A,b) != NUM_OK) NP_RETURN(1,result[0]);
    #endif

  return (0);
}

static INT GSConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SmootherInit;
  theNP->Display = SmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = GSPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = GSStep;

  return(0);
}

/****************************************************************************/
/*D
   sgs - numproc for Gauss-Seidel smoother

   DESCRIPTION:
   This numproc executes a symmetric Gauss-Seidel smoother, using the blas routines
   'l_lgs' and 'l_ugs'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT SGSPreProcess  (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *baselevel, INT *result)
{
        #ifdef ModelP
  NP_SGS *np;
  GRID *theGrid;

  np = (NP_SGS *) theNP;
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L))
    NP_RETURN(1,result[0]);
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS) != NUM_OK)
    NP_RETURN(1,result[0]);
        #endif
  *baselevel = level;

  /* get storage for extra temp */
  if (AllocVDFromVD(theNP->base.mg,level,level,x,&NP_SGS_t((NP_SGS *)theNP)))
    NP_RETURN(1,result[0]);

  return (0);
}


static INT SGSStep (NP_SMOOTHER *theNP, INT level,
                    VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                    MATDATA_DESC *L,
                    INT *result)
{
  return (0);
}

static INT SGSSmoother (NP_ITER *theNP, INT level,
                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                        INT *result)
{
  NP_SGS *np;
  GRID *theGrid;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  np = (NP_SGS *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);

  /* iterate forward */
    #ifdef ModelP
  if (l_vector_collect(theGrid,b)!=NUM_OK) NP_RETURN(1,result[0]);
  if (l_lgs(theGrid,NP_SGS_t(np),(const MATDATA_DESC *)&np->smoother.L,b)
      != NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_vector_consistent(theGrid,NP_SGS_t(np)) != NUM_OK)
    NP_RETURN(1,result[0]);
    #else
  if (l_lgs(theGrid,NP_SGS_t(np),A,b)) NP_RETURN(1,result[0]);
    #endif

  /* damp */
  if (l_dscale(theGrid,NP_SGS_t(np),ACTIVE_CLASS,np->smoother.damp) != NUM_OK) NP_RETURN(1,result[0]);

  /* update defect */
  if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,NP_SGS_t(np),ACTIVE_CLASS)
      != NUM_OK) NP_RETURN(1,result[0]);

  /* iterate backward */
    #ifdef ModelP
  if (l_vector_collect(theGrid,b)!=NUM_OK) NP_RETURN(1,result[0]);
  if (l_ugs(theGrid,x,(const MATDATA_DESC *)&np->smoother.L,b))
    NP_RETURN(1,result[0]);
  if (l_vector_consistent(theGrid,x) != NUM_OK) NP_RETURN(1,result[0]);
        #else
  if (l_ugs(theGrid,x,A,b)) NP_RETURN(1,result[0]);
    #endif

  /* damp */
  if (l_dscale(theGrid,x,ACTIVE_CLASS,np->smoother.damp) != NUM_OK) NP_RETURN(1,result[0]);

  /* update defect */
  if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,x,ACTIVE_CLASS)
      != NUM_OK) NP_RETURN(1,result[0]);

  /* now add the two corrections */
  if (l_daxpy(theGrid,x,ACTIVE_CLASS,Factor_One,NP_SGS_t(np)) != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SGSPostProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *result)
{
  NP_SGS *np;

  np = (NP_SGS *) theNP;
  if (np->smoother.L != NULL)
    FreeMD(theNP->base.mg,level,level,np->smoother.L);
  FreeVD(NP_MG(theNP),level,level,NP_SGS_t(np));

  return(0);
}

static INT SGSConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SmootherInit;
  theNP->Display = SmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = SGSPreProcess;
  np->iter.Iter = SGSSmoother;
  np->iter.PostProcess = SGSPostProcess;
  np->Step = SGSStep;

  return(0);
}

/****************************************************************************/
/*D
   sor - numproc for SOR smoother

   DESCRIPTION:
   This numproc executes an SOR (successive over relaxation) smoother,
   using the blas routine
   'l_lsor'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT SORPreProcess  (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *baselevel, INT *result)
{
        #ifdef ModelP
  NP_SMOOTHER *np;
  GRID *theGrid;

  np = (NP_SMOOTHER *) theNP;
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->L)) NP_RETURN(1,result[0]);
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_dmatcopy(theGrid,np->L,A) != NUM_OK) NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->L,TRUE) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  *baselevel = level;

  return (0);
}


static INT SORStep (NP_SMOOTHER *theNP, INT level,
                    VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                    MATDATA_DESC *L,
                    INT *result)
{
    #ifdef ModelP
  if (l_lsor(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,L,b,theNP->damp) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (l_lsor(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,A,b,theNP->damp) != NUM_OK) NP_RETURN(1,result[0]);
    #endif

  return (0);
}

static INT SORSmoother (NP_ITER *theNP, INT level,
                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                        INT *result)
{
  NP_SMOOTHER *np;
  GRID *theGrid;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  np = (NP_SMOOTHER *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if ((*np->Step)(np,level,x,b,A,np->L,result))
    return (1);
    #ifdef ModelP
  if (l_vector_consistent(theGrid,x) != NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,x,ACTIVE_CLASS)
      != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SORConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SmootherInit;
  theNP->Display = SmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = SORPreProcess;
  np->iter.Iter = SORSmoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = SORStep;

  return(0);
}

/****************************************************************************/
/*D
   ebgs - numproc for equation block Gauss Seidel smoother

   DESCRIPTION:
   This numproc executes an equation block Gauss Seidel smoother.
   It can be used in 'lmgc'.

   .vb
   npinit $n <n-iter> $A <mat sym> $r <rhs sym> $x <sol sym> $t <temp sym>
   $D {no|red|full} $damp <sc double list>
   $Blocking <sc int list> $BlockOrder <ord list> $BlockIter <sc numproc list>
   .ve

   .  $n~<n-iter> - number of iterations
   .  $A~<mat~sym> - symbol for the stiffness matrix
   .  $r~<rhs~sym> - symbol for the right hand side vector
   .  $x~<sol~sym> - symbol for the solution vector
   .  $t~<temp~sym> - symbol for a tempory vector
   .  $D - no, reduced or full display
   .  $damp~<sc~double~list> - damping factors for each component
   .  $Blocking~<sc~int~list> - describes the block structure (in ascending order starting always with 0)
   .  $BlockOrder~<ord~list> - describes the order of the blocks
   .  $BlockIter~<numproc~list> - smoother per block
   .  <sc~double~list>  - [nd <double list>] | [ed <double list>] | [el <double list>] | [si <double list>]
   .  <sc~int~list>     - [nd <int list>] | [ed <int list>] | [el <int list>] | [si <int list>]
   . <sc~numproc~list>  - [nd <numproc list>] | [ed <numproc list>] | [el <numproc list>] | [si <numproc list>]
   . <ord~list>         - nd|ed|el|sd<number>{ nd|ed|el|sd<number>}+

   'npexecute <name>'

   EXAMPLE:
   For the Navier-Sokes problem class library we have three unknowns per node
   (the velocity vector (u,v) and the pressure p). We want to use the equation-block
   Gauss-Seidel in the following way: update all pressure values first, then update
   u and finally update p. The block inverses should be computed inexactly by a point-block
   ILU method (actually the point-blocks are then 1x1, i.e. scalar!).

   This can be done in the following way

   .vb
 # creation of smoother numproc realizations

 # inner (inexact) solvers for the ebgs: pbilu
   npcreate u_ilu $t pbilu $f ns;
   npcreate v_ilu $t pbilu $f ns;
   npcreate p_ilu $t pbilu $f ns;

 # equation block Gauss-Seidel
   npcreate ebgs  $t ebgs  $f ns;


 # initialization of the smoother num procs

 # inner (inexact) solvers for the ebgs: pbilu
   scnp u_ilu;
   npinit $A MATu $L LUu $t utmp $r urhs $x ucor $n 1 $damp nd 1.0 $beta nd 1.0 $D no;
   scnp v_ilu;
   npinit $A MATv $L LUv $t vtmp $r vrhs $x vcor $n 1 $damp nd 1.0 $beta nd 1.0 $D no;
   scnp p_ilu;
   npinit $A MATp $L LUp $t ptmp $r prhs $x pcor $n 1 $damp nd 0.5 $beta nd 5.0 $D no;

 # equation block Gauss-Seidel with inner solvers u_ilu v_ilu p_ilu
   scnp pre_ebgs;
   npinit $A MAT $t tmp $r rhs $x cor $n 1 $Blocking nd 0 1 2 3 $BlockOrder nd2 nd0 nd1 $BlockIter nd u_ilu v_ilu p_ilu $D no;
   .ve

   The 'Blocking nd 0 1 2 3' means that block 0 is defined in the node and consists
   of the component 0 of the corresponding 'VECDATA_DESC' only. Block number 1
   consists of component 1 and finally block number 2 consists of component 2.
   (If we had chosen 'Blocking nd 0 2 3' block 0 would include components 0 and 1
   while block 2 would consist of component 2). Blocking descriptions underly
   two constraints: the associated 'VECDATA_DESC' has to be exhausted by the description
   and blocks can never exceed single types.

   The 'BlockOrder nd2 nd0 nd1' means: solve first for block 2 of the node then
   for block 0 of the node and finally for block 2 of the node.

   'BlockIter nd u_ilu v_ilu p_ilu' means that the (in this case inexact) solvers
   for the blocks are u_ilu for block 0, v_ilu for block 1 and p_ilu for block 2.

   SEE ALSO:
   num_procs

   D*/
/****************************************************************************/
#ifdef __BGS__
static INT BGS_Init (NP_BASE *theNP, INT argc , char **argv)
{
  NP_BGS *theBGS;
  char option[OPTIONLEN],value[VALUELEN];
  NP_BASE *BlockIter[MAX_BLOCKS][NVECTYPES];
  INT i,bopt,boopt,biopt,nBlocks,nIter,LastBlockEnd,LastTypeBlock,type,nTypeBlocksPred;
  INT nTypeBlocks[NVECTYPES],TypeBlocks[MAX_BLOCKS+1][NVECTYPES];
  INT nBlockIter[NVECTYPES];

  theBGS = (NP_BGS*) theNP;

  /* set configuration parameters */
  bopt = boopt = biopt = FALSE;
  for (i=1; i<argc; i++)
    if (sscanf(argv[i],expandfmt(CONCAT5("%",OPTIONLENSTR,"[a-zA-Z0-9_] %",VALUELENSTR,"[ -~]")),option,value)==2)
    {
      /* Blocking */
      if (strstr(option,"Blocking")!=NULL)
        if (ReadVecTypeINTs(value,MAX_BLOCKS+1,nTypeBlocks,TypeBlocks)!=0)
          return (NP_NOT_ACTIVE);
        else {bopt = TRUE; continue;}

      /* BlockOrder */
      if (strstr(option,"BlockOrder")!=NULL)
        if (ReadVecTypeOrder(value,MAX_ORDER,MAX_BLOCKS,&NP_BGS_NBLOCKITER(theBGS),NP_BGS_BLOCKORDER(theBGS))!=0)
          return (NP_NOT_ACTIVE);
        else {boopt = TRUE; continue;}

      /* BlockIter */
      if (strstr(option,"BlockIter")!=NULL)
        if (ReadVecTypeNUMPROCs(NP_MG(theNP),value,ITER_CLASS_NAME,MAX_BLOCKS,nBlockIter,BlockIter)!=0)
          return (NP_NOT_ACTIVE);
        else {biopt = TRUE; continue;}
    }

  if (!(bopt && boopt && biopt))
  {
    PrintErrorMessage('E',"NP_BGS_Init","one or several options missing");
    return (NP_NOT_ACTIVE);
  }

  /* combine TypeBlocks to 'global' blocks */
  nBlocks = 0;
  LastBlockEnd = NP_BGS_BLOCKSTART(theBGS,nBlocks++) = 0;
  for (type=0; type<NVECTYPES; type++)
    for (LastTypeBlock=0, i=(TypeBlocks[0][type]==0) ? 1 : 0; i<nTypeBlocks[type]; i++)
    {
      LastBlockEnd  = NP_BGS_BLOCKSTART(theBGS,nBlocks++) = LastBlockEnd + (TypeBlocks[i][type]-LastTypeBlock);
      LastTypeBlock = TypeBlocks[i][type];
    }
  NP_BGS_NBLOCKS(theBGS) = --nBlocks;    /* don't count first 0 */

  /* condense BlockOrder to 'global' BlockOrder */
  for (type=0; type<NVECTYPES; type++)
    if ((nTypeBlocks[type]>0) && (TypeBlocks[0][type]==0))
      nTypeBlocks[type]--;
  for (i=0; i<NP_BGS_NBLOCKITER(theBGS); i++)
  {
    type = NP_BGS_BLOCKORD(theBGS,i) / MAX_BLOCKS;
    nTypeBlocksPred = (type>0) ? nTypeBlocks[type-1] : 0;
    NP_BGS_BLOCKORD(theBGS,i) = nTypeBlocksPred + NP_BGS_BLOCKORD(theBGS,i) % MAX_BLOCKS;
  }

  /* combine BlockIter to 'global' BlockIter */
  nIter = 0;
  for (type=0; type<NVECTYPES; type++)
    for (i=0; i<nBlockIter[type]; i++)
      NP_BGS_BLOCKITER(theBGS,nIter++) = (NP_SMOOTHER *)BlockIter[i][type];

  /* check blocks of iteration order */
  for (i=0; i<NP_BGS_NBLOCKITER(theBGS); i++)
    if (NP_BGS_BLOCKORD(theBGS,i)>=nBlocks)
    {
      PrintErrorMessage('E',"NP_BGS_Init","block id in BlockOrder too large");
      return (NP_NOT_ACTIVE);
    }

  /* check number of block iteration schemes */
  if (nIter!=nBlocks)
  {
    PrintErrorMessage('E',"NP_BGS_Init","number of specified block iteration schemes does not match number of blocks");
    return (NP_NOT_ACTIVE);
  }

  return (SmootherInit(theNP,argc,argv));
}

static INT BGS_Display (NP_BASE *theNP)
{
  NP_BGS *theBGS;
  char name[16];
  INT i;

  theBGS = (NP_BGS*) theNP;

  SmootherDisplay(theNP);

  /* now display additional stuff for BGS */
  UserWrite("Blocking:\n");
  for (i=0; i<NP_BGS_NBLOCKS(theBGS); i++)
  {
    sprintf(name," block%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SII,name,NP_BGS_BLOCKSTART(theBGS,i),NP_BGS_BLOCKEND(theBGS,i));
  }

  UserWrite("BlockOrder:\n");
  for (i=0; i<NP_BGS_NBLOCKITER(theBGS); i++)
  {
    sprintf(name," blockord%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SI,name,NP_BGS_BLOCKORD(theBGS,i));
  }

  UserWrite("BlockIterations:\n");
  for (i=0; i<NP_BGS_NBLOCKS(theBGS); i++)
  {
    sprintf(name," blockiter%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SS,name,NP_BGS_BLOCKITNAME(theBGS,i));
  }

  return (0);
}

static INT BGSPreProcess (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                          INT *baselevel, INT *result)
{
  GRID *theGrid;
  NP_BGS *theBGS;
  VECDATA_DESC *vd;
  MATDATA_DESC *md;
  SHORT *cmpptr_r,*cmpptr_A;
  INT offset[MAX_BLOCKS],rowoffset,coloffset,ncols_a,compR,rowcA;
  INT i,j,bl,rbl,cbl,rtype,ctype,currComp,currType,NextBlockStart,blComp,nrc,ncc,cmp;
  INT bopt,boopt,biopt,nBlocks,nIter,LastBlockEnd,LastTypeBlock,type,nTypeBlocksPred;

  theBGS = (NP_BGS*) theNP;
  theGrid = GRID_ON_LEVEL(NP_MG(theNP),level);

  /* construct VEC_DESCs and VECDATA_DESCs of tmp block vectors */
  currComp = currType = 0;
  for (bl=0; bl<nBlocks; bl++)
  {
    /* init DESCriptor */
    for (i=0; i<NVECTYPES; i++)
      VD_NCMPS_IN_TYPE(NP_BGS_VD_tb(theBGS,bl),i)=0;

    /* set currType to next used type */
    while ((VD_NCMPS_IN_TYPE(x,currType)==0)
           && currType<NVECTYPES)
      currType++;
    if (currType>=NVECTYPES)
    {
      PrintErrorMessage('E',"NP_BGS_Init","too many blocks specified");
      return (NP_NOT_ACTIVE);
    }
    if ((VD_NCMPS_IN_TYPE(NPIT_b(theBGS),currType)==0)
        ||(VD_NCMPS_IN_TYPE(x,currType)!=
           VD_NCMPS_IN_TYPE(NPIT_b(theBGS),currType)))
    {
      PrintErrorMessage('E',"NP_BGS_Init",
                        "rhs and tmp descriptors do not match");
      return (NP_NOT_ACTIVE);
    }

    NextBlockStart = NP_BGS_BLOCKSTART(theBGS,bl+1);
    if (NextBlockStart > VD_NCMPS_IN_TYPE(x,currType))
    {
      PrintErrorMessage('E',"NP_BGS_Init","block exceeds type-block");
      return (NP_NOT_ACTIVE);
    }

    offset[bl] = currComp;
    blComp = NextBlockStart-currComp;
    VD_CMPPTR_OF_TYPE(NP_BGS_VD_tb(theBGS,bl),currType) =
      VD_CMPPTR_OF_TYPE(x,currType) + currComp;
    VD_NCMPS_IN_TYPE(NP_BGS_VD_tb(theBGS,bl),currType) = blComp;
    currComp += blComp;

    if (NextBlockStart==VD_NCMPS_IN_TYPE(x,currType))
    {
      currType++;
      currComp = 0;
    }
  }
  /* check consistency of the decomposition into blocks */
  if (NextBlockStart!=VD_NCMPS_IN_TYPE(x,currType-1))
  {
    PrintErrorMessage('E',"NP_BGS_Init","descriptor is not exhausted by blocks");
    return (NP_NOT_ACTIVE);
  }
  while ((VD_NCMPS_IN_TYPE(x,currType)==0)
         && currType<NVECTYPES)
    currType++;
  if (currType<NVECTYPES)
  {
    PrintErrorMessage('E',"NP_BGS_Init","type descriptor is not exhausted by blocks");
    return (NP_NOT_ACTIVE);
  }

  /* construct block-rhs VEC_DESC and MAT_DESCs and MATDATA_DESCs needed */
  cmpptr_r = NP_BGS_COMPS_r(theBGS);
  cmpptr_A = NP_BGS_COMPS_A(theBGS);
  for (cbl=0; cbl<nBlocks; cbl++)
  {
    /* find column type of this block */
    for (ctype=0; ctype<NVECTYPES; ctype++)
      if (VD_NCMPS_IN_TYPE(NP_BGS_VD_tb(theBGS,cbl),ctype)>0)
        break;

    /* init DESCriptors */
    for (i=0; i<NVECTYPES; i++)
      VD_NCMPS_IN_TYPE(NP_BGS_VD_rb(theBGS,cbl),i) = 0;
    for (i=0; i<NVECTYPES; i++)
      for (j=0; j<NVECTYPES; j++)
      {
        MD_ROWS_IN_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),i,j) = 0;
        MD_COLS_IN_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),i,j) = 0;
      }

    ncc = VD_NCMPS_IN_TYPE(NP_BGS_VD_tb(theBGS,cbl),ctype);

    currType = NVECTYPES;

    for (rbl=0; rbl<nBlocks; rbl++)
    {
      if (rbl==cbl) continue;

      /* find row-type of this block */
      for (rtype=0; rtype<NVECTYPES; rtype++)
        if (VD_NCMPS_IN_TYPE(NP_BGS_VD_tb(theBGS,rbl),rtype)>0)
          break;

      if (MD_ROWS_IN_RT_CT(NPIT_A(theBGS),rtype,ctype)==0) continue;

      if (rtype!=currType)
      {
        /* the type has changed */

        if (currType!=NVECTYPES)
        {
          /* set number of comps for last type handled */
          VD_NCMPS_IN_TYPE(NP_BGS_VD_rb(theBGS,cbl),ctype) = compR;
          MD_ROWS_IN_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),rtype,ctype) = rowcA;
          MD_COLS_IN_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),rtype,ctype) = ncc;
        }

        /* init comp ptrs for new type */
        VD_CMPPTR_OF_TYPE(NP_BGS_VD_rb(theBGS,cbl),rtype) = cmpptr_r;
        MD_MCMPPTR_OF_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),rtype,ctype) = cmpptr_A;

        /* init counters for new type */
        currType = rtype;
        rowcA = compR = 0;
      }

      ncols_a = MD_COLS_IN_RT_CT(NPIT_A(theBGS),rtype,ctype);

      coloffset = offset[cbl];
      rowoffset = offset[rbl];
      nrc = VD_NCMPS_IN_TYPE(NP_BGS_VD_tb(theBGS,rbl),rtype);

      for (i=0; i<nrc; i++)
      {
        /* VD_COMP(rb,compR++) = VD_COMP(theVD_r,rowoffset+i);*/
        VD_CMP_OF_TYPE(NP_BGS_VD_rb(theBGS,cbl),ctype,compR++) = VD_CMP_OF_TYPE(NPIT_b(theBGS),rtype,rowoffset+i);
        for (j=0; j<ncc; j++)
        {
          cmp = (rowoffset+i) * ncols_a + (coloffset+j);
          MD_MCMP_OF_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),
                           rtype,ctype,rowcA*ncc+j)
            = MD_MCMP_OF_RT_CT(NPIT_A(theBGS),rtype,ctype,cmp);
        }
        rowcA++;
      }
    }
    VD_NCMPS_IN_TYPE(NP_BGS_VD_rb(theBGS,cbl),ctype) = compR;
    MD_ROWS_IN_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),rtype,ctype) = rowcA;
    MD_COLS_IN_RT_CT(NP_BGS_MD_Ab(theBGS,cbl),rtype,ctype) = ncc;

    cmpptr_r += compR;
    cmpptr_A += rowcA*ncc;
  }

  /* fill redundant information in DESCriptors */
  for (rbl=0; rbl<nBlocks; rbl++)
  {
    FillRedundantComponentsOfVD(NP_BGS_VD_tb(theBGS,rbl));
    FillRedundantComponentsOfVD(NP_BGS_VD_rb(theBGS,rbl));
    FillRedundantComponentsOfMD(NP_BGS_MD_Ab(theBGS,rbl));
  }

  /* call prepares of block iteration schemes */
  for (bl=0; bl<NP_BGS_NBLOCKS(theBGS); bl++)
    if ((*NP_BGS_BLOCKITER(theBGS,bl)->iter.PreProcess)((NP_ITER*)NP_BGS_BLOCKITER(theBGS,bl),theGrid)!=0)
      return (bl+1);

  return (0);
}

static INT BGSSmoother (NP_ITER *theNP, INT level,
                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                        INT *result)
{
  NP_BGS *theBGS;
  GRID *theGrid;
  INT i,blo,bl,BlockItShellUsed;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  theBGS = (NP_BGS*) theNP;
  theGrid = GRID_ON_LEVEL(NP_MG(theNP),level);

  BlockItShellUsed = FALSE;
  for (blo=0; blo<NP_BGS_NBLOCKITER(theBGS); blo++)
  {
    bl = NP_BGS_BLOCKORD(theBGS,blo);

    /* iterate */
    if ((*NP_BGS_BLOCKITER(theBGS,bl)->iter.Iter)((NP_ITER*)NP_BGS_BLOCKITER(theBGS,bl),theGrid,&BlockItShellUsed)!=0)
      return (bl+1);

    /* now temp contains the corresponding update of the corr-field
       the corr-field is updated already
       we have to update the remaining defects */

    if (l_dmatmul_minus(theGrid,NP_BGS_VD_rb(theBGS,bl),NEWDEF_CLASS,
                        NP_BGS_MD_Ab(theBGS,bl),
                        NP_BGS_VD_tb(theBGS,bl),ACTIVE_CLASS)) return (1);
  }

  return (0);
}

static INT BGSConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = BGS_Init;
  theNP->Display = BGS_Display;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = BGSPreProcess;
  np->iter.Iter = BGSSmoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = SGSStep;

  return(0);
}
#endif
/****************************************************************************/
/*D
   ilu - numproc for point block beta-modified ilu smoother

   DESCRIPTION:
   This numproc executes a point block ilu smoother, using the blas routines
   'l_ilubthdecomp' and 'l_luiter'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list> $beta <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT ILUInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ILU *np;
  INT i;

  np = (NP_ILU *) theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->beta[i] = 0.0;
  sc_read(np->beta,np->smoother.iter.b,"beta",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT ILUDisplay (NP_BASE *theNP)
{
  NP_ILU *np;

  SmootherDisplay(theNP);
  np = (NP_ILU *) theNP;
  if (sc_disp(np->beta,np->smoother.iter.b,"beta")) return (1);

  return (0);
}

static INT ILUPreProcess (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                          INT *baselevel, INT *result)
{
  NP_ILU *np;
  GRID *theGrid;

  np = (NP_ILU *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_ilubthdecomp(theGrid,np->smoother.L,np->beta,NULL,NULL,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"ILUPreProcess","decomposition failed");
    result[0] = __LINE__;
    return (1);
  }
  *baselevel = level;

  return (0);
}

static INT ILUStep (NP_SMOOTHER *theNP, INT level,
                    VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                    MATDATA_DESC *L,
                    INT *result)
{
    #ifdef ModelP
  if (l_vector_collect(GRID_ON_LEVEL(theNP->iter.base.mg,level),b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_luiter(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT ILUConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = ILUInit;
  theNP->Display = ILUDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = ILUPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = ILUStep;

  return(0);
}

/****************************************************************************/
/*D
   filu - numproc for point block beta-modified ilu smoother working on FINE NODES ONLY

   DESCRIPTION:
   This numproc executes a point block ilu smoother, using the blas routines
   'l_ilubthdecomp_fine' and 'l_luiter_fine'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list> $beta <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT FILUPreProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           INT *baselevel, INT *result)
{
  NP_ILU *np;
  GRID *theGrid;

  np = (NP_ILU *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_ilubthdecomp_fine(theGrid,np->smoother.L,np->beta,NULL,NULL,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"FILUPreProcess","decomposition failed");
    result[0] = __LINE__;
    return (1);
  }
  *baselevel = level;

  return (0);
}

static INT FILUStep (NP_SMOOTHER *theNP, INT level,
                     VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                     MATDATA_DESC *L,
                     INT *result)
{
    #ifdef ModelP
  if (l_vector_collect(GRID_ON_LEVEL(theNP->iter.base.mg,level),b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_luiter_fine(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT FILUConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = ILUInit;
  theNP->Display = ILUDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = FILUPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = FILUStep;

  return(0);
}

/****************************************************************************/
/*D
   thilu - numproc for point block beta-modified ilu smoother with threshold for extending
                        the sparsity pattern

   DESCRIPTION:
   This numproc executes a point block ilu smoother, using the blas routines
   'l_ilubthdecomp' and 'l_luiter'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list> $beta <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT THILUInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_THILU *np;
  INT i;

  np = (NP_THILU *) theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->beta[i] = np->thresh[i] = 0.0;
  sc_read(np->beta,np->smoother.iter.b,"beta",argc,argv);
  sc_read(np->thresh,np->smoother.iter.b,"thresh",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT THILUDisplay (NP_BASE *theNP)
{
  NP_THILU *np;

  SmootherDisplay(theNP);
  np = (NP_THILU *) theNP;
  if (sc_disp(np->beta,np->smoother.iter.b,"beta")) return (1);
  if (sc_disp(np->thresh,np->smoother.iter.b,"thresh")) return (1);

  return (0);
}

static INT THILUPreProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *baselevel, INT *result)
{
  NP_THILU *np;
  GRID *theGrid;

  np = (NP_THILU *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_ilubthdecomp(theGrid,np->smoother.L,np->beta,np->thresh,NULL,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"THILUPreProcess","decomposition failed");
    result[0] = __LINE__;
    return (1);
  }
  *baselevel = level;

  return (0);
}

static INT THILUConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = THILUInit;
  theNP->Display = THILUDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = THILUPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = ILUStep;

  return(0);
}

/****************************************************************************/
/*D
   spilu - numproc for point block spectrally shifted ilu smoother

   DESCRIPTION:
   This numproc executes a point block spectrally shifted ilu smoother, using the blas routines
   'l_iluspbthdecomp' and 'l_luiter'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list> $beta <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT SPILUInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SPILU *np;
  char buffer[32];
  INT i;

  np = (NP_SPILU *) theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->beta[i] = 0.0;
  sc_read(np->beta,np->smoother.iter.b,"beta",argc,argv);

  np->mode = SP_LOCAL;
  if (ReadArgvChar("mode",buffer,argc,argv))
  {
    PrintErrorMessage('E',"SPILUInit","specify mode");
    return (NP_NOT_ACTIVE);
  }
  if (strncmp(buffer,"global",3)==0)
    np->mode = SP_GLOBAL;
  else if (strncmp(buffer,"local",3)==0)
    np->mode = SP_LOCAL;
  else
  {
    PrintErrorMessage('E',"SPILUInit","specify local/global for mode");
    return (NP_NOT_ACTIVE);
  }

  return (SmootherInit(theNP,argc,argv));
}

static INT SPILUDisplay (NP_BASE *theNP)
{
  NP_SPILU *np;

  SmootherDisplay(theNP);
  np = (NP_SPILU *) theNP;
  if (sc_disp(np->beta,np->smoother.iter.b,"beta")) return (1);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"mode",
             (np->mode==SP_GLOBAL) ? "global" : "local");

  return (0);
}

static INT SPILUPreProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *baselevel, INT *result)
{
  NP_SPILU *np;
  GRID *theGrid;

  np = (NP_SPILU *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_iluspdecomp(theGrid,np->smoother.L,np->beta,NULL,np->mode,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"SPILUPreProcess","decomposition failed");
    result[0] = __LINE__;
    return (1);
  }
  *baselevel = level;

  return (0);
}

static INT SPILUConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SPILUInit;
  theNP->Display = SPILUDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = SPILUPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = ILUStep;

  return(0);
}

/****************************************************************************/
/*D
   ic - numproc for point block icncomplete Cholesky smoother

   DESCRIPTION:
   This numproc executes a point block icncomplete Cholesky smoother, using the blas routines
   'l_icdecomp' and 'l_lltiter'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list> $beta <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT ICPreProcess (NP_ITER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *baselevel, INT *result)
{
  NP_ILU *np;
  GRID *theGrid;

  np = (NP_ILU *) theNP;
  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_icdecomp(theGrid,np->smoother.L)
      !=NUM_OK) {
    PrintErrorMessage('E',"ICPreProcess","decomposition failed");
    result[0] = __LINE__;
    return (1);
  }
  *baselevel = level;

  return (0);
}

static INT ICStep (NP_SMOOTHER *theNP, INT level,
                   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                   MATDATA_DESC *L,
                   INT *result)
{
    #ifdef ModelP
  if (l_vector_collect(GRID_ON_LEVEL(theNP->iter.base.mg,level),b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_lltiter(GRID_ON_LEVEL(theNP->iter.base.mg,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT ICConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SmootherInit;
  theNP->Display = SmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = ICPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = ICStep;

  return(0);
}

/****************************************************************************/
/*D
   lu - numproc for lu smoother

   DESCRIPTION:
   This numproc executes lu smoother, using the blas routines
   'l_lrdecomp' and 'l_luiter'. It can be used in 'lmgc'.

   .vb
   npinit [$c <cor>] [$b <rhs>] [$A <mat>]
       $n <it> $damp <sc double list>
   .ve

   .  $c~<sol> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $n~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT LUPreProcess (NP_ITER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *baselevel, INT *result)
{
  NP_SMOOTHER *np;
  GRID *theGrid;
  INT err;
  char warn[255];

  np = (NP_SMOOTHER *) theNP;

  theGrid = GRID_ON_LEVEL(theNP->base.mg,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->L,MAT_MASTER_CONS) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  err = l_lrdecomp(theGrid,np->L);
  if (err != NUM_OK) {
    if (err>0) {
      switch (err) {
      case NUM_OUT_OF_MEM :
        PrintErrorMessage('E',"LUPreProcess","out of memory");
        result[0] = __LINE__;
        return (1);
      default :
        PrintErrorMessage('E',"LUPreProcess","err > 0");
        result[0] = __LINE__;
        return (1);
      }
    }
    if (err!=-VINDEX(LASTVECTOR(theGrid))) {
      sprintf(warn,"decomp failed: IDX %ld on level %d",
              -err,GLEVEL(theGrid));
      PrintErrorMessage('E',"LUPreProcess",warn);
      UserWriteF(" - LASTVECTOR has IDX %ld\n",
                 VINDEX(LASTVECTOR(theGrid)));
      result[0] = __LINE__;
      return (1);
    }
    if (l_lrregularize(theGrid,np->L) !=NUM_OK) {
      PrintErrorMessage('E',"LUPreProcess","cannot regularize");
      result[0] = __LINE__;
      return (1);
    }
  }
  *baselevel = level;

  return (0);
}

static INT LUConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SmootherInit;
  theNP->Display = SmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = LUPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = ILUStep;

  return(0);
}

/****************************************************************************/
/*D
   lmgc - numproc for linear multigrid cycle

   DESCRIPTION:
   This numproc executes

   .vb
   npinit [$c <cor>] [$r <rhs>] [$A <mat>]
       $S <pre post base> $T <transfer>
       [$b <baselevel>] [$g <gamma>] [$n1 <it>] [$n2 <it>]
   .ve

   .  $c~<sol> - correction vector
   .  $r~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $T~<transfer> - transfer numproc
   .  $S~<pre~post~base> - numprocs for pre- and postsmoother, base solver
   .  $b~<baselevel> - baselevel where the base solver is called
   .  $g~<gamma> - number of iterations of Lmgc per level (default gamma = 1)
   .  $n1~<it> - number of iterations of the presmoother (default n1 = 1)
   .  $n2~<it> - number of iteration of the postsmoother (default n2 = 1)

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT LmgcInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LMGC *np;
  INT i;
  char post[VALUELEN],pre[VALUELEN],base[VALUELEN];

  np = (NP_LMGC *) theNP;

  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->Transfer = (NP_TRANSFER *)
                 ReadArgvNumProc(theNP->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  for (i=1; i<argc; i++)
    if (argv[i][0]=='S') {
      if (sscanf(argv[i],"S %s %s %s",pre,post,base)!=3)
        continue;
      np->PreSmooth = (NP_ITER *)
                      GetNumProcByName(theNP->mg,pre,ITER_CLASS_NAME);
      np->PostSmooth = (NP_ITER *)
                       GetNumProcByName(theNP->mg,post,ITER_CLASS_NAME);
      np->BaseSolver = (NP_LINEAR_SOLVER *)
                       GetNumProcByName(theNP->mg,base,LINEAR_SOLVER_CLASS_NAME);
      break;
    }

  if (ReadArgvINT("g",&(np->gamma),argc,argv))
    np->gamma = 1;
  if (ReadArgvINT("n1",&(np->nu1),argc,argv))
    np->nu1 = 1;
  if (ReadArgvINT("n2",&(np->nu2),argc,argv))
    np->nu2 = 1;
  if (ReadArgvINT("b",&(np->baselevel),argc,argv))
    np->baselevel = 0;

  if (np->Transfer == NULL) return(NP_NOT_ACTIVE);
  if (np->PreSmooth == NULL) return(NP_NOT_ACTIVE);
  if (np->PostSmooth == NULL) return(NP_NOT_ACTIVE);
  if (np->BaseSolver == NULL) return(NP_NOT_ACTIVE);

  return (NPIterInit(&np->iter,argc,argv));
}

static INT LmgcDisplay (NP_BASE *theNP)
{
  NP_LMGC *np;

  np = (NP_LMGC *) theNP;

  NPIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->gamma);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)np->nu1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n2",(int)np->nu2);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"b",(int)np->baselevel);

  if (np->Transfer != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(np->Transfer));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (np->PreSmooth != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"pre",ENVITEM_NAME(np->PreSmooth));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"pre","---");
  if (np->PostSmooth != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"post",ENVITEM_NAME(np->PostSmooth));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"post","---");
  if (np->BaseSolver != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"base",ENVITEM_NAME(np->BaseSolver));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"base","---");

  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));

  return (0);
}

static INT LmgcPreProcess  (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *baselevel, INT *result)
{
  NP_LMGC *np;
  INT i;

  np = (NP_LMGC *) theNP;

  if (np->Transfer->PreProcess != NULL)
    if ((*np->Transfer->PreProcess)
          (np->Transfer,np->baselevel,level,x,b,A,result))
      return(1);

  if (np->PreSmooth->PreProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PreProcess)
            (np->PreSmooth,i,x,b,A,baselevel,result))
        return(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PreProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PreSmooth->PreProcess)
              (np->PostSmooth,i,x,b,A,baselevel,result))
          return(1);

  *baselevel = MIN(np->baselevel,level);
  if (np->BaseSolver->PreProcess != NULL)
    if ((*np->BaseSolver->PreProcess)
          (np->BaseSolver,*baselevel,x,b,A,baselevel,result))
      return(1);

  return (0);
}

static INT Lmgc (NP_ITER *theNP, INT level,
                 VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A,
                 INT *result)
{
  NP_LMGC *np;
  MULTIGRID *theMG;
  GRID *theGrid;
  LRESULT lresult;
  INT i;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = c;
  NPIT_b(theNP) = b;

  np = (NP_LMGC *) theNP;

  if (level <= np->baselevel) {
    if ((*np->BaseSolver->Residuum)
          (np->BaseSolver,np->baselevel,level,c,b,A,&lresult))
      return(1);
    if ((*np->BaseSolver->Solver)(np->BaseSolver,level,c,b,A,
                                  np->BaseSolver->abslimit,
                                  np->BaseSolver->reduction,&lresult)) NP_RETURN(1,result[0]);
    if (!lresult.converged)
      PrintErrorMessage('W',"Lmgc","no convergence of BaseSolver");
    return(0);
  }
  theMG = theNP->base.mg;
  theGrid = GRID_ON_LEVEL(theMG,level);
  if (AllocVDFromVD(theMG,level,level,c,&np->t)) NP_RETURN(1,result[0]);
  for (i=0; i<np->nu1; i++) {
    if ((*np->PreSmooth->Iter)(np->PreSmooth,level,np->t,b,A,result))
      return(1);
    if (l_daxpy(theGrid,c,ACTIVE_CLASS,Factor_One,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }
  if ((*np->Transfer->RestrictDefect)
        (np->Transfer,level,b,b,A,Factor_One,result))
    return(1);

  if (l_dset(DOWNGRID(theGrid),c,EVERY_CLASS,0.0) != NUM_OK) NP_RETURN(1,result[0]);
  for (i=0; i<np->gamma; i++)
    if (Lmgc(theNP,level-1,c,b,A,result))
      return(1);
  if ((*np->Transfer->InterpolateCorrection)
        (np->Transfer,level,np->t,c,A,Factor_One,result))
    return(1);
  if (l_daxpy(theGrid,c,EVERY_CLASS,Factor_One,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,np->t,EVERY_CLASS) != NUM_OK) NP_RETURN(1,result[0]);
  for (i=0; i<np->nu2; i++) {
    if ((*np->PostSmooth->Iter)(np->PostSmooth,level,np->t,b,A,result))
      return(1);
    if (l_daxpy(theGrid,c,ACTIVE_CLASS,Factor_One,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }
  FreeVD(theNP->base.mg,level,level,np->t);
  if (np->Transfer->AdaptCorrection != NULL)
    if ((*np->Transfer->AdaptCorrection)(np->Transfer,level,c,b,A,result))
      return(1);

  return (0);
}

static INT LmgcPostProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *result)
{
  NP_LMGC *np;
  INT i;

  np = (NP_LMGC *) theNP;

  if (np->Transfer->PostProcess != NULL)
    for (i = np->baselevel; i <= level; i++)
      if ((*np->Transfer->PostProcess)(np->Transfer,i,x,b,A,result))
        return(1);

  if (np->PreSmooth->PostProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PostProcess)
            (np->PreSmooth,i,x,b,A,result))
        return(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PostProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PreSmooth->PostProcess)
              (np->PostSmooth,i,x,b,A,result))
          return(1);

  if (np->BaseSolver->PostProcess != NULL)
    if ((*np->BaseSolver->PostProcess)
          (np->BaseSolver,np->baselevel,x,b,A,result))
      return(1);

  return (0);
}

static INT LmgcConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = LmgcInit;
  theNP->Display = LmgcDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = LmgcPreProcess;
  np->Iter = Lmgc;
  np->PostProcess = LmgcPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitIter	- Init this file

   SYNOPSIS:
   INT InitIter ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

INT InitIter ()
{
  INT i;

  if (CreateClass(ITER_CLASS_NAME ".jac",sizeof(NP_SMOOTHER),JacobiConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".gs",sizeof(NP_SMOOTHER),GSConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sgs",sizeof(NP_SGS),SGSConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sor",sizeof(NP_SMOOTHER),SORConstruct))
    return (__LINE__);
#ifdef __BGS__
  if (CreateClass(ITER_CLASS_NAME ".bgs",sizeof(NP_BGS),BGSConstruct))
    return (__LINE__);
#endif
  if (CreateClass(ITER_CLASS_NAME ".ilu",sizeof(NP_ILU),ILUConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".filu",sizeof(NP_ILU),FILUConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".thilu",sizeof(NP_THILU),THILUConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".spilu",sizeof(NP_ILU),SPILUConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ic",sizeof(NP_ILU),ICConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".lu",sizeof(NP_SMOOTHER),LUConstruct))
    return (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".lmgc",sizeof(NP_LMGC),LmgcConstruct))
    return (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  return (0);
}
