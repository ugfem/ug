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
#include <math.h>

#include "general.h"
#include "debug.h"
#include "gm.h"
#include "algebra.h"
#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "devices.h"
#include "udm.h"
#include "pcr.h"
#include "debug.h"

#include "transfer.h"
#include "ls.h"
#include "iter.h"

#include "ff_gen.h"
#include "ff.h"


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define TYPE_TFF 1
#define TYPE_FF 2

#define NPFF_tv(p)                              (((p)->tv))
#define NPFF_tv2(p)                             (((p)->tv2))
#define NPFF_t(p)                               (((p)->t))

#define NPFF_TYPE(p)                    ((p)->type)
#define NPFF_DO_TFF(p)                  (NPFF_TYPE(p)==TYPE_TFF)
#define NPFF_DO_FF(p)                   (NPFF_TYPE(p)==TYPE_FF)
#define NPFF_MESHWIDTH(p)               ((p)->meshwidth)
#define NPFF_WaveNrRel(p)               ((p)->wave_nr_rel)
#define NPFF_WaveNrRel3D(p)             ((p)->wave_nr_rel3D)
#define NPFF_ALLFREQ(p)                 ((p)->all_freq)
#define NPFF_DISPLAY(p)                 ((p)->display)
#define NPFF_BVDF(p)                    (&(p)->bvdf)

/* macros for the symmetric Gauss-Seidel smoother */
#define NP_SGS_t(p)                             ((p)->t)

/* macros for the Block Gauss-Seidel smoother */
#define MAX_BLOCKS                        3
#define MAX_ORDER                         6
#define SBGS_NBLOCKS(p)               ((p)->nBlocks)
#define SBGS_BLOCKS(p)                ((p)->Block)
#define SBGS_BLOCKDESC(p,i)                   ((p)->BlockDesc[i])
#define SBGS_BLOCKITERS(p)            ((p)->BlockIter)
#define SBGS_BLOCKITER(p,i)           ((p)->BlockIter[i])
#define SBGS_BLOCKITNAME(p,i)         ENVITEM_NAME(SBGS_BLOCKITER(p,i))
#define SBGS_NBLOCKITER(p)            ((p)->nBlockIter)
#define SBGS_BLOCKORDER(p)            ((p)->BlockOrder)
#define SBGS_BLOCKORD(p,i)            ((p)->BlockOrder[i])
#define SBGS_MD_Ad(p,cb)              (&((p)->MD_Ad[cb]))
#define SBGS_MD_Ao(p,cb)              (&((p)->MD_Ao[cb]))
#define SBGS_VD_cd(p)                 (&((p)->VD_cd))
#define SBGS_VD_rd(p,b)               (&((p)->VD_rd[b]))
#define SBGS_VD_ro(p,b)               (&((p)->VD_ro[b]))
#define SBGS_COMPS_Ad(p)                      ((p)->COMP_Ad)
#define SBGS_COMPS_Ao(p)                      ((p)->COMP_Ao)
#define SBGS_COMPS_cd(p)                      ((p)->COMP_cd)
#define SBGS_COMPS_rd(p)                      ((p)->COMP_rd)
#define SBGS_COMPS_ro(p)                      ((p)->COMP_ro)

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
  INT tp;                                                       /* block type							*/
  INT fc;                                                       /* first block comp in type				*/
  INT tc;                                                       /* last block comp +1 in type			*/

} SBLOCK_DESC;

typedef struct
{
  NP_ITER iter;

  INT nBlocks;                      /* number of blocks                     */
  INT Block[MAX_BLOCKS+1];          /* block subdivision                    */
  NP_SMOOTHER *BlockIter[MAX_BLOCKS];  /* block iteration scheme            */
  INT nBlockIter;                   /* number of block iterations           */
  INT BlockOrder[MAX_ORDER];        /* iteration order for the blocks       */
  SBLOCK_DESC BlockDesc[MAX_BLOCKS];  /* block descriptors					*/
  MATDATA_DESC MD_Ad[MAX_BLOCKS];   /* diagonal blocks of stiffness matrix  */
  MATDATA_DESC MD_Ao[MAX_BLOCKS];   /* off-diag blocks of stiffness matrix  */
  VECDATA_DESC VD_cd;                       /* diagonal blocks of right hand side   */
  VECDATA_DESC VD_rd[MAX_BLOCKS];   /* diagonal blocks of right hand side   */
  VECDATA_DESC VD_ro[MAX_BLOCKS];   /* off-diag blocks of right hand side   */

  /* storage for components */
  SHORT COMP_Ad[MAX_BLOCKS*MAX_MAT_COMP];
  SHORT COMP_Ao[MAX_BLOCKS*MAX_MAT_COMP];
  SHORT COMP_cd[MAX_VEC_COMP];
  SHORT COMP_rd[MAX_BLOCKS*MAX_VEC_COMP];
  SHORT COMP_ro[MAX_BLOCKS*MAX_VEC_COMP];

} NP_SBGS;

typedef struct
{
  NP_SMOOTHER smoother;

  VEC_SCALAR beta;
  VEC_SCALAR mindiag;

} NP_ILU;

typedef struct
{
  NP_SMOOTHER smoother;

  /* abstract memory description */
  VECDATA_DESC *tv;                     /* testvector */
  VECDATA_DESC *tv2;                    /* 2. testvector for FF */
  VECDATA_DESC *t;                      /* temp. vector for the update of the correction */

  /* configuration */
  INT type;                                     /* TYPE_TFF or TYPE_FF */
  DOUBLE meshwidth;                     /* meshwidth of the grid */
  DOUBLE wave_nr_rel;                   /* wavenumber for the testing frequency */
  DOUBLE wave_nr_rel3D;         /* wavenumber for the testing frequency; only for 3D */
  INT all_freq;                         /* flag; TRUE == smooth for all relevant frequencies */
  INT display;
#ifdef __BLOCK_VECTOR_DESC__
  BV_DESC_FORMAT bvdf;
#endif

} NP_FF;

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

struct np_bcgssmoother {

  NP_SMOOTHER np_smoother;

  NP_ITER *Iter;
  DOUBLE rho, omega;
  INT maxiter;
  INT restart;
  VECDATA_DESC *r;
  VECDATA_DESC *p;
  VECDATA_DESC *v;
  VECDATA_DESC *s;
  VECDATA_DESC *t;
  VECDATA_DESC *q;
};
typedef struct np_bcgssmoother NP_BCGSSMOOTHER;

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

static MATDATA_DESC *FF_MATDATA_DESC_ARRAY[FF_MAX_MATS];
static VECDATA_DESC *FF_VECDATA_DESC_ARRAY[FF_MAX_VECS];

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
   The data can be displayed and the num proc can be executed by

   'INT NPIterDisplay (NP_ITER *theNP);'
   'INT NPIterExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   struct np_iter {

        NP_BASE base;                        // inherits base class

        // data (optinal, necessary for calling the generic execute routine)
    VECDATA_DESC *c;                     // correction
    VECDATA_DESC *b;                     // defect
    MATDATA_DESC *A;                     // matrix

        // functions
        INT (*PreProcess)
             (struct np_iter *,              // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *,                         // baselevel used by iter
                  INT *);                        // result
    INT (*Iter)
             (struct np_iter *,              // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // correction vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
        INT (*PostProcess)
             (struct np_iter *,              // pointer to (derived) object
                  INT,                           // level
                  VECDATA_DESC *,                // solution vector
                  VECDATA_DESC *,                // defect vector
                  MATDATA_DESC *,                // matrix
                  INT *);                        // result
   };
   typedef struct np_iter NP_ITER;
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
    REP_ERR_RETURN (1);
  }
  if (np->b == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no vector b");
    REP_ERR_RETURN (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPIterExecute","no matrix A");
    REP_ERR_RETURN (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no PreProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PreProcess)(np,level,np->c,np->b,np->A,&bl,&result)) {
      UserWriteF("NPIterExecute: PreProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->Iter == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no Iter");
      REP_ERR_RETURN (1);
    }
    if ((*np->Iter)(np,level,np->c,np->b,np->A,&result)) {
      UserWriteF("NPIterExecute: Iter failed, error code %d\n", result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no PostProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PostProcess)(np,level,np->c,np->b,np->A,&result)) {
      UserWriteF("NPIterExecute: PostProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
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
  if (sc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);
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
  theGrid = NP_GRID(theNP,level);
  if ((*np->Step)(np,level,x,b,A,np->L,result))
    REP_ERR_RETURN (1);
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
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

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
  theGrid = NP_GRID(theNP,level);
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
  if (l_jac(NP_GRID(theNP,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (l_jac(NP_GRID(theNP,level),x,A,b) != NUM_OK) NP_RETURN(1,result[0]);
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
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

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
  theGrid = NP_GRID(theNP,level);
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
  if (l_vector_collect(NP_GRID(theNP,level),b)!=NUM_OK) NP_RETURN(1,result[0]);

  if (l_lgs(NP_GRID(theNP,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (l_lgs(NP_GRID(theNP,level),x,A,b) != NUM_OK) NP_RETURN(1,result[0]);
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
   bcgss - numproc for bi-cg-stab-smoother

   DESCRIPTION:
   This numproc executes a bi-cg-stab as smoother.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$I iter]
       $m <it> $damp <sc double list> [$I <iteration>] [$R <restart>];
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $m~<it> - number of iterations
   .  $damp~<sc~double~list> - damping factors for each component
       $m <maxit>  [$w <sc double list>];
   .  $I~<iteration> - iteration numproc
   .  $R~<restart> - restart index

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT BCGSSmootherInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_BCGSSMOOTHER *np;

  np = (NP_BCGSSMOOTHER *) theNP;

  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  np->p = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->v = ReadArgvVecDesc(theNP->mg,"v",argc,argv);
  np->s = ReadArgvVecDesc(theNP->mg,"s",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->q = ReadArgvVecDesc(theNP->mg,"q",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvINT("R",&(np->restart),argc,argv))
    np->restart = 0;
  if (np->restart<0) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->Iter = (NP_ITER *) ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT BCGSSmootherDisplay (NP_BASE *theNP)
{
  NP_BCGSSMOOTHER *np;

  SmootherDisplay(theNP);
  np = (NP_BCGSSMOOTHER *) theNP;
  if (np->r != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->p != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));
  if (np->v != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"v",ENVITEM_NAME(np->v));
  if (np->s != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(np->s));
  if (np->t != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->q != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"q",ENVITEM_NAME(np->q));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"R",(int)np->restart);
  if (np->Iter != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter",ENVITEM_NAME(np->Iter));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"Iter","---");

  return (0);
}

static INT BCGSSPreProcess  (NP_ITER *theNP, INT level,
                             VECDATA_DESC *x, VECDATA_DESC *b,
                             MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_BCGSSMOOTHER *np;

  *baselevel = level;
  np = (NP_BCGSSMOOTHER *) theNP;

  if (np->Iter!=NULL)
    if (np->Iter->PreProcess != NULL)
      if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);

  if (AllocVDFromVD(np->np_smoother.iter.base.mg,level,level,x,&np->r)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->np_smoother.iter.base.mg,level,level,x,&np->p)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->np_smoother.iter.base.mg,level,level,x,&np->v)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->np_smoother.iter.base.mg,level,level,x,&np->s)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->np_smoother.iter.base.mg,level,level,x,&np->t)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(np->np_smoother.iter.base.mg,level,level,x,&np->q)) NP_RETURN(1,result[0]);

  return (0);
}


static INT BCGSSStep (NP_SMOOTHER *theNP, INT level,
                      VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                      MATDATA_DESC *L,
                      INT *result)
{
  NP_BCGSSMOOTHER *np;
  VEC_SCALAR scal;
  INT i,j,restart;
  DOUBLE alpha,rho_new,beta,tt;

  np = (NP_BCGSSMOOTHER *) theNP;
  restart = 1;
  for (i=0; i<np->maxiter; i++)
  {
    /* restart ? */
    if ((np->restart>0 && i%np->restart==0) || restart)
    {
      if (l_dset(NP_GRID(theNP,level),np->p,EVERY_CLASS,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (l_dset(NP_GRID(theNP,level),np->v,EVERY_CLASS,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (l_dcopy(NP_GRID(theNP,level),np->r,EVERY_CLASS,b)!= NUM_OK) REP_ERR_RETURN (1);
      alpha = np->rho = np->omega = 1.0;
      restart = 0;
    }

    /* update x, b */
    if (l_ddot_sv (NP_GRID(theNP,level),b,EVERY_CLASS,np->r,Factor_One,&rho_new)!=NUM_OK) REP_ERR_RETURN (1);
    beta=rho_new*alpha/np->rho/np->omega;
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=beta;
    if (l_dscale (NP_GRID(theNP,level),np->p,EVERY_CLASS,scal)) REP_ERR_RETURN (1);
    if (l_daxpy (NP_GRID(theNP,level),np->p,EVERY_CLASS,Factor_One,b)!= NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-beta*np->omega;
    if (l_daxpy (NP_GRID(theNP,level),np->p,EVERY_CLASS,scal,np->v)!= NUM_OK) REP_ERR_RETURN (1);
    if (np->Iter!=NULL)
    {
      if (l_dset(NP_GRID(theNP,level),np->q,EVERY_CLASS,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (l_dcopy(NP_GRID(theNP,level),np->s,EVERY_CLASS,np->p)!= NUM_OK) REP_ERR_RETURN (1);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->p,A,result)) REP_ERR_RETURN (1);
      if (l_dcopy(NP_GRID(theNP,level),np->p,EVERY_CLASS,np->s)!= NUM_OK) REP_ERR_RETURN (1);
      if (l_dmatmul_set(NP_GRID(theNP,level),np->v,EVERY_CLASS,A,np->q,EVERY_CLASS)) REP_ERR_RETURN (1);
      if (l_ddot_sv (NP_GRID(theNP,level),np->v,EVERY_CLASS,np->r,Factor_One,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      alpha = rho_new/alpha;
      for (j=0; j<VD_NCOMP(x); j++) scal[j]=alpha;
      if (l_daxpy (NP_GRID(theNP,level),x,EVERY_CLASS,scal,np->q)!= NUM_OK) REP_ERR_RETURN (1);
    }
    else
    {
      if (l_dmatmul_set(NP_GRID(theNP,level),np->v,EVERY_CLASS,A,np->p,EVERY_CLASS)) REP_ERR_RETURN (1);
      if (l_ddot_sv (NP_GRID(theNP,level),np->v,EVERY_CLASS,np->r,Factor_One,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      alpha = rho_new/alpha;
      for (j=0; j<VD_NCOMP(x); j++) scal[j]=alpha;
      if (l_daxpy (NP_GRID(theNP,level),x,EVERY_CLASS,scal,np->p)!= NUM_OK) REP_ERR_RETURN (1);
    }
    if (l_dcopy(NP_GRID(theNP,level),np->s,EVERY_CLASS,b)!= NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-alpha;
    if (l_daxpy (NP_GRID(theNP,level),np->s,EVERY_CLASS,scal,np->v)!= NUM_OK) REP_ERR_RETURN (1);
    if (np->Iter!=NULL)
    {
      if (l_dset(NP_GRID(theNP,level),np->q,EVERY_CLASS,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (l_dcopy(NP_GRID(theNP,level),np->t,EVERY_CLASS,np->s)!= NUM_OK) REP_ERR_RETURN (1);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->s,A,result)) REP_ERR_RETURN (1);
      if (l_dcopy(NP_GRID(theNP,level),np->s,EVERY_CLASS,np->t)!= NUM_OK) REP_ERR_RETURN (1);
    }
    else
    {
      if (l_dcopy(NP_GRID(theNP,level),np->q,EVERY_CLASS,np->s)!= NUM_OK) REP_ERR_RETURN (1);
    }
    if (l_dmatmul_set(NP_GRID(theNP,level),np->t,EVERY_CLASS,A,np->q,EVERY_CLASS)) REP_ERR_RETURN (1);
    if (l_ddot_sv (NP_GRID(theNP,level),np->t,EVERY_CLASS,np->t,Factor_One,&tt)!=NUM_OK) REP_ERR_RETURN (1);
    if (l_ddot_sv (NP_GRID(theNP,level),np->s,EVERY_CLASS,np->t,Factor_One,&(np->omega))!=NUM_OK) REP_ERR_RETURN (1);
    np->omega /= tt;
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=np->omega;
    if (l_daxpy (NP_GRID(theNP,level),x,EVERY_CLASS,scal,np->q)!= NUM_OK) REP_ERR_RETURN (1);
    if (l_dcopy(NP_GRID(theNP,level),b,EVERY_CLASS,np->s)!= NUM_OK) REP_ERR_RETURN (1);
    for (j=0; j<VD_NCOMP(x); j++) scal[j]=-np->omega;
    if (l_daxpy (NP_GRID(theNP,level),b,EVERY_CLASS,scal,np->t)!= NUM_OK) REP_ERR_RETURN (1);
    np->rho = rho_new;
  }


  return (0);
}

static INT BCGSSmootherPostProcess (NP_ITER *theNP, INT level,
                                    VECDATA_DESC *x, VECDATA_DESC *b,
                                    MATDATA_DESC *A, INT *result)
{
  NP_SMOOTHER *nps;
  NP_BCGSSMOOTHER *np;

  nps = (NP_SMOOTHER *) theNP;
  if (nps->L != NULL) FreeMD(theNP->base.mg,level,level,nps->L);
  np = (NP_BCGSSMOOTHER *) theNP;
  FreeVD(nps->iter.base.mg,level,level,np->r);
  FreeVD(nps->iter.base.mg,level,level,np->p);
  FreeVD(nps->iter.base.mg,level,level,np->v);
  FreeVD(nps->iter.base.mg,level,level,np->s);
  FreeVD(nps->iter.base.mg,level,level,np->t);
  FreeVD(nps->iter.base.mg,level,level,np->q);

  if (np->Iter!=NULL)
  {
    if (np->Iter->PostProcess == NULL) return(0);
    return((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result));
  }
  else
    return (0);

  return(0);
}

static INT BCGSSConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = BCGSSmootherInit;
  theNP->Display = BCGSSmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = BCGSSPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = BCGSSmootherPostProcess;
  np->Step = BCGSSStep;

  return(0);
}

/****************************************************************************/
/*D
   sgs - numproc for symmetric Gauss-Seidel smoother

   DESCRIPTION:
   This numproc executes a symmetric Gauss-Seidel smoother,
   using the blas routines
   'l_lgs' and 'l_ugs'. It can be used in 'lmgc'.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT SGSInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SGS *np;

  np = (NP_SGS *) theNP;
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT SGSDisplay (NP_BASE *theNP)
{
  NP_SGS *np;

  SmootherDisplay(theNP);
  np = (NP_SGS *) theNP;
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));

  return (0);
}

static INT SGSPreProcess  (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *baselevel, INT *result)
{
        #ifdef ModelP
  NP_SMOOTHER *np;
  GRID *theGrid;

  np = (NP_SMOOTHER *) theNP;
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->L))
    NP_RETURN(1,result[0]);
  theGrid = NP_GRID(theNP,level);
  if (l_dmatcopy(theGrid,np->L,A) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->L,MAT_MASTER_CONS) != NUM_OK)
    NP_RETURN(1,result[0]);
        #endif
  *baselevel = level;

  /* get storage for extra temp */
  if (AllocVDFromVD(theNP->base.mg,level,level,x,&NP_SGS_t((NP_SGS *)theNP)))
    NP_RETURN(1,result[0]);

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
  theGrid = NP_GRID(theNP,level);

  /* iterate forward */
    #ifdef ModelP
  if (l_vector_collect(theGrid,b)!=NUM_OK) NP_RETURN(1,result[0]);
  if (l_lgs(theGrid,NP_SGS_t(np),np->smoother.L,b)
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
  if (l_ugs(theGrid,x,np->smoother.L,b))
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
  NP_ITER *np;

  theNP->Init = SGSInit;
  theNP->Display = SGSDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = SGSPreProcess;
  np->Iter = SGSSmoother;
  np->PostProcess = SGSPostProcess;

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
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

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
  theGrid = NP_GRID(theNP,level);
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
  if (l_lsor(NP_GRID(theNP,level),x,L,b,theNP->damp) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (l_lsor(NP_GRID(theNP,level),x,A,b,theNP->damp) != NUM_OK) NP_RETURN(1,result[0]);
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
  theGrid = NP_GRID(theNP,level);
  if ((*np->Step)(np,level,x,b,A,np->L,result))
    REP_ERR_RETURN (1);
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
   sbgs - numproc for equation block Gauss Seidel smoother

   DESCRIPTION:
   This numproc executes an equation block Gauss Seidel smoother.
   It can be used in 'lmgc'.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
               $damp <sc double list>
               $Blocking <sc int list> $BlockOrder <ord list>
               $BlockIter <sc numproc list>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $Blocking~<sc~int~list> - describes the block structure (in ascending order starting always with 0)
   .  $BlockOrder~<ord~list> - describes the order of the blocks
   .  $BlockIter~<numproc~list> - smoother per block
   .  <sc~double~list>  - [nd <double list>] | [ed <double list>] | [el <double list>] | [si <double list>]
   .  <sc~int~list>     - [nd <int list>] | [ed <int list>] | [el <int list>] | [si <int list>]
   . <sc~numproc~list>  - [nd <numproc list>] | [ed <numproc list>] | [el <numproc list>] | [si <numproc list>]
   . <ord~list>         - nd|ed|el|sd<number>{ nd|ed|el|sd<number>}+
   .  <double~list>  - <double> {: <double>}*

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess

   EXAMPLE:
   For the Navier-Sokes problem class library we have three unknowns per node
   (the velocity vector (u,v) and the pressure p). We want to use the equation-block
   Gauss-Seidel in the following way: update all pressure values first, then update
   u and finally update p. The block inverses should be computed inexactly by a point-block
   ILU method (actually the point-blocks are then 1x1, i.e. scalar!).

   This can be done in the following way

   .vb
 # creation of smoother numproc realizations

 # inner (inexact) solvers for the sbgs: pbilu
   npcreate u_ilu $c ilu;
   npcreate v_ilu $c ilu;
   npcreate p_ilu $c ilu;

 # equation block Gauss-Seidel
   npcreate sbgs  $c sbgs;


 # initialization of the smoother num procs

 # inner (inexact) solvers for the sbgs: pbilu
   scnp u_ilu;
   npinit $damp nd 1.0 $beta nd 1.0;
   scnp v_ilu;
   npinit $damp nd 1.0 $beta nd 1.0;
   scnp p_ilu;
   npinit $damp nd 0.5 $beta nd 5.0;

 # equation block Gauss-Seidel with inner solvers u_ilu v_ilu p_ilu
   scnp pre_sbgs;
   npinit $Blocking nd 0 1 2 3 $BlockOrder nd2 nd0 nd1 $BlockIter nd u_ilu v_ilu p_ilu;
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

static INT SBGS_Init (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SBGS *theSBGS;
  char option[OPTIONLEN],value[VALUELEN];
  NP_BASE *BlockIter[MAX_BLOCKS][NVECTYPES];
  INT i,bopt,boopt,biopt,nBlocks,nIter,type,nTypeBlocksPred;
  INT nTypeBlocks[NVECTYPES],TypeBlocks[MAX_BLOCKS+1][NVECTYPES];
  INT nBlockIter[NVECTYPES];

  theSBGS = (NP_SBGS*) theNP;

  /* set configuration parameters */
  bopt = boopt = biopt = FALSE;
  for (i=1; i<argc; i++)
    if (sscanf(argv[i],expandfmt(CONCAT5("%",OPTIONLENSTR,"[a-zA-Z0-9_] %",VALUELENSTR,"[ -~]")),option,value)==2)
    {
      /* Blocking */
      if (strstr(option,"Blocking")!=NULL)
        if (ReadVecTypeINTs(value,MAX_BLOCKS+1,nTypeBlocks,TypeBlocks)!=0)
          REP_ERR_RETURN (NP_NOT_ACTIVE)
          else {bopt = TRUE; continue;}

      /* BlockOrder */
      if (strstr(option,"BlockOrder")!=NULL)
        if (ReadVecTypeOrder(value,MAX_ORDER,MAX_BLOCKS,&SBGS_NBLOCKITER(theSBGS),SBGS_BLOCKORDER(theSBGS))!=0)
          REP_ERR_RETURN (NP_NOT_ACTIVE)
          else {boopt = TRUE; continue;}

      /* BlockIter */
      if (strstr(option,"BlockIter")!=NULL)
        if (ReadVecTypeNUMPROCs(NP_MG(theNP),value,ITER_CLASS_NAME,MAX_BLOCKS,nBlockIter,BlockIter)!=0)
          REP_ERR_RETURN (NP_NOT_ACTIVE)
          else {biopt = TRUE; continue;}
    }

  if (!(bopt && boopt && biopt))
  {
    PrintErrorMessage('E',"SBGS_Init","one or several options missing");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }

  /* fill SBLOCK_DESC data structure */
  nBlocks = 0;
  for (type=0; type<NVECTYPES; type++)
    for (i=0; i<nTypeBlocks[type]-1; i++)
    {
      SBGS_BLOCKDESC(theSBGS,nBlocks).tp = type;
      SBGS_BLOCKDESC(theSBGS,nBlocks).fc = TypeBlocks[i][type];
      SBGS_BLOCKDESC(theSBGS,nBlocks).tc = TypeBlocks[i+1][type];
      nBlocks++;
    }
  SBGS_NBLOCKS(theSBGS) = nBlocks;

  /* condense BlockOrder to 'global' BlockOrder */
  for (type=0; type<NVECTYPES; type++)
    if ((nTypeBlocks[type]>0) && (TypeBlocks[0][type]==0))
      nTypeBlocks[type]--;
  for (i=0; i<SBGS_NBLOCKITER(theSBGS); i++)
  {
    type = SBGS_BLOCKORD(theSBGS,i) / MAX_BLOCKS;
    nTypeBlocksPred = (type>0) ? nTypeBlocks[type-1] : 0;
    SBGS_BLOCKORD(theSBGS,i) = nTypeBlocksPred + SBGS_BLOCKORD(theSBGS,i) % MAX_BLOCKS;
  }

  /* check blocks of iteration order */
  for (i=0; i<SBGS_NBLOCKITER(theSBGS); i++)
    if (SBGS_BLOCKORD(theSBGS,i)>=nBlocks)
    {
      PrintErrorMessage('E',"SBGS_Init","block id in BlockOrder too large");
      REP_ERR_RETURN (NP_NOT_ACTIVE);
    }

  /* combine BlockIter to 'global' BlockIter */
  nIter = 0;
  for (type=0; type<NVECTYPES; type++)
    for (i=0; i<nBlockIter[type]; i++)
      SBGS_BLOCKITER(theSBGS,nIter++) = (NP_SMOOTHER *)BlockIter[i][type];

  /* check number of block iteration schemes */
  if (nIter!=nBlocks)
  {
    PrintErrorMessage('E',"SBGS_Init","number of specified block iteration schemes does not match number of blocks");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }

  return (NPIterInit(&theSBGS->iter,argc,argv));
}

static INT SBGS_Display (NP_BASE *theNP)
{
  NP_SBGS *theSBGS;
  char name[16];
  INT i;

  theSBGS = (NP_SBGS*) theNP;

  NPIterDisplay(&theSBGS->iter);

  /* now display additional stuff for SBGS */
  UserWrite("Blocking:\n");
  for (i=0; i<SBGS_NBLOCKS(theSBGS); i++)
  {
    sprintf(name," block%d(%s)",i,VecTypeName[SBGS_BLOCKDESC(theSBGS,i).tp]);
    UserWriteF(DISPLAY_NP_FORMAT_SII,name,SBGS_BLOCKDESC(theSBGS,i).fc,SBGS_BLOCKDESC(theSBGS,i).tc);
  }

  UserWrite("BlockOrder:\n");
  for (i=0; i<SBGS_NBLOCKITER(theSBGS); i++)
  {
    sprintf(name," blockord%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SI,name,SBGS_BLOCKORD(theSBGS,i));
  }

  UserWrite("BlockIterations:\n");
  for (i=0; i<SBGS_NBLOCKS(theSBGS); i++)
  {
    sprintf(name," blockiter%d",i);
    UserWriteF(DISPLAY_NP_FORMAT_SS,name,SBGS_BLOCKITNAME(theSBGS,i));
  }

  return (0);
}

static INT SetCorComps (NP_SBGS *theSBGS, const VECDATA_DESC *c, INT bl)
{
  INT i,type,kcd;

  type = SBGS_BLOCKDESC(theSBGS,bl).tp;

  for (i=0; i<NVECTYPES; i++)
    VD_NCMPS_IN_TYPE(SBGS_VD_cd(theSBGS),i) = 0;

  VD_NCMPS_IN_TYPE(SBGS_VD_cd(theSBGS),type) = SBGS_BLOCKDESC(theSBGS,bl).tc
                                               - SBGS_BLOCKDESC(theSBGS,bl).fc;
  VD_CMPPTR_OF_TYPE(SBGS_VD_cd(theSBGS),type) = SBGS_COMPS_cd(theSBGS);
  kcd = 0;
  for (i=0; i<VD_NCMPS_IN_TYPE(c,type); i++)
    if ((i>=SBGS_BLOCKDESC(theSBGS,bl).fc) && (i<SBGS_BLOCKDESC(theSBGS,bl).tc))
      VD_CMP_OF_TYPE(SBGS_VD_cd(theSBGS),type,kcd++) = VD_CMP_OF_TYPE(c,type,i);

  FillRedundantComponentsOfVD(SBGS_VD_cd(theSBGS));

  return (0);
}

static INT SBGSPreProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           INT *baselevel, INT *result)
{
  GRID *theGrid;
  NP_SBGS *theSBGS;
  SHORT *cmpptr_rd,*cmpptr_ro,*cmpptr_Ad,*cmpptr_Ao;
  INT i,j,bl,rtype,ctype;
  INT n,nd,no,krd,kro,kad,kao;
  INT nBlocks;

  theSBGS = (NP_SBGS*) theNP;
  theGrid = NP_GRID(theNP,level);
  nBlocks = SBGS_NBLOCKS(theSBGS);

  /* construct auxiliary block XXXDATA_DESCs */
  cmpptr_rd = SBGS_COMPS_rd(theSBGS);
  cmpptr_ro = SBGS_COMPS_ro(theSBGS);
  cmpptr_Ad = SBGS_COMPS_Ad(theSBGS);
  cmpptr_Ao = SBGS_COMPS_Ao(theSBGS);
  for (bl=0; bl<nBlocks; bl++)
  {
    ctype = SBGS_BLOCKDESC(theSBGS,bl).tp;

    n  = VD_NCMPS_IN_TYPE(b,ctype);
    nd = SBGS_BLOCKDESC(theSBGS,bl).tc - SBGS_BLOCKDESC(theSBGS,bl).fc;
    no = n - nd;

    for (i=0; i<NVECTYPES; i++)
      VD_NCMPS_IN_TYPE(SBGS_VD_rd(theSBGS,bl),i) =
        VD_NCMPS_IN_TYPE(SBGS_VD_ro(theSBGS,bl),i) = 0;
    for (i=0; i<NMATTYPES; i++)
      MD_ROWS_IN_MTYPE(SBGS_MD_Ad(theSBGS,bl),i) =
        MD_COLS_IN_MTYPE(SBGS_MD_Ad(theSBGS,bl),i) =
          MD_ROWS_IN_MTYPE(SBGS_MD_Ao(theSBGS,bl),i) =
            MD_COLS_IN_MTYPE(SBGS_MD_Ao(theSBGS,bl),i) = 0;

    VD_NCMPS_IN_TYPE(SBGS_VD_rd(theSBGS,bl),ctype) = nd;
    VD_CMPPTR_OF_TYPE(SBGS_VD_rd(theSBGS,bl),ctype) = cmpptr_rd;
    VD_NCMPS_IN_TYPE(SBGS_VD_ro(theSBGS,bl),ctype) = no;
    VD_CMPPTR_OF_TYPE(SBGS_VD_ro(theSBGS,bl),ctype) = cmpptr_ro;

    kro = krd = kao = kad = 0;
    for (rtype=0; rtype<NVECTYPES; rtype++)
    {
      MD_MCMPPTR_OF_RT_CT(SBGS_MD_Ad(theSBGS,bl),rtype,ctype) = cmpptr_Ad;
      MD_MCMPPTR_OF_RT_CT(SBGS_MD_Ao(theSBGS,bl),rtype,ctype) = cmpptr_Ao;

      if (rtype==ctype)
      {
        /* here the diagonal block is contained */
        MD_ROWS_IN_RT_CT(SBGS_MD_Ad(theSBGS,bl),rtype,ctype) = nd;
        MD_COLS_IN_RT_CT(SBGS_MD_Ad(theSBGS,bl),rtype,ctype) = nd;
        MD_ROWS_IN_RT_CT(SBGS_MD_Ao(theSBGS,bl),rtype,ctype) = no;
        MD_COLS_IN_RT_CT(SBGS_MD_Ao(theSBGS,bl),rtype,ctype) = nd;

        for (j=0; j<n; j++)
          if ((j<SBGS_BLOCKDESC(theSBGS,bl).fc) || (j>=SBGS_BLOCKDESC(theSBGS,bl).tc))
            VD_CMP_OF_TYPE(SBGS_VD_ro(theSBGS,bl),ctype,kro++) = VD_CMP_OF_TYPE(b,ctype,j);
          else
          {
            VD_CMP_OF_TYPE(SBGS_VD_rd(theSBGS,bl),ctype,krd++) = VD_CMP_OF_TYPE(b,ctype,j);
            for (i=0; i<n; i++)
              if ((i<SBGS_BLOCKDESC(theSBGS,bl).fc) || (i>=SBGS_BLOCKDESC(theSBGS,bl).tc))
                MD_MCMP_OF_RT_CT(SBGS_MD_Ao(theSBGS,bl),rtype,ctype,kao++) = MD_IJ_CMP_OF_RT_CT(A,rtype,ctype,i,j);
              else
                MD_MCMP_OF_RT_CT(SBGS_MD_Ad(theSBGS,bl),rtype,ctype,kad++) = MD_IJ_CMP_OF_RT_CT(A,rtype,ctype,i,j);
          }
        ASSERT(krd==nd);
        ASSERT(kro==no);
        ASSERT(kad==nd*nd);
        ASSERT(kao==nd*no);
      }
      else
      {
        /* this is in any case off diag */
        MD_ROWS_IN_RT_CT(SBGS_MD_Ao(theSBGS,bl),rtype,ctype) = MD_ROWS_IN_RT_CT(A,rtype,ctype);
        MD_COLS_IN_RT_CT(SBGS_MD_Ao(theSBGS,bl),rtype,ctype) = MD_COLS_IN_RT_CT(A,rtype,ctype);

        for (j=0; j<n; j++)
        {
          VD_CMP_OF_TYPE(SBGS_VD_ro(theSBGS,bl),ctype,kro++) = VD_CMP_OF_TYPE(b,ctype,j);
          if ((j>=SBGS_BLOCKDESC(theSBGS,bl).fc) && (j<SBGS_BLOCKDESC(theSBGS,bl).tc))
            for (i=0; i<MD_ROWS_IN_RT_CT(A,rtype,ctype); i++)
              MD_MCMP_OF_RT_CT(SBGS_MD_Ao(theSBGS,bl),rtype,ctype,kao++) = MD_IJ_CMP_OF_RT_CT(A,rtype,ctype,i,j);
        }
      }
    }
    cmpptr_rd += krd;
    cmpptr_ro += kro;
    cmpptr_Ad += kad;
    cmpptr_Ao += kao;
  }

  /* fill redundant information in DESCriptors */
  for (bl=0; bl<nBlocks; bl++)
  {
    FillRedundantComponentsOfVD(SBGS_VD_rd(theSBGS,bl));
    FillRedundantComponentsOfVD(SBGS_VD_ro(theSBGS,bl));
    FillRedundantComponentsOfMD(SBGS_MD_Ad(theSBGS,bl));
    FillRedundantComponentsOfMD(SBGS_MD_Ao(theSBGS,bl));
  }

  /* call prepares of block iteration schemes */
  for (bl=0; bl<SBGS_NBLOCKS(theSBGS); bl++)
  {
    SetCorComps(theSBGS,x,bl);

    if ((*SBGS_BLOCKITER(theSBGS,bl)->iter.PreProcess)
          ((NP_ITER*)SBGS_BLOCKITER(theSBGS,bl),
          level,
          SBGS_VD_cd(theSBGS),
          SBGS_VD_rd(theSBGS,bl),
          SBGS_MD_Ad(theSBGS,bl),
          baselevel,
          result
          )!=0)
      REP_ERR_RETURN (bl+1);
  }

  return (0);
}

static INT SBGSSmoother (NP_ITER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *result)
{
  NP_SBGS *theSBGS;
  GRID *theGrid;
  INT blo,bl;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  theSBGS = (NP_SBGS*) theNP;
  theGrid = NP_GRID(theNP,level);

  for (blo=0; blo<SBGS_NBLOCKITER(theSBGS); blo++)
  {
    bl = SBGS_BLOCKORD(theSBGS,blo);

    SetCorComps(theSBGS,x,bl);

    /* iterate */
    if ((*SBGS_BLOCKITER(theSBGS,bl)->iter.Iter)
          ((NP_ITER*)SBGS_BLOCKITER(theSBGS,bl),
          level,
          SBGS_VD_cd(theSBGS),
          SBGS_VD_rd(theSBGS,bl),
          SBGS_MD_Ad(theSBGS,bl),
          result
          )!=0)
      NP_RETURN (bl+1,result[0]);

    /* now temp contains the corresponding update of the corr-field
       the corr-field is updated already
       we have to update the remaining defects */

    if (l_dmatmul_minus(theGrid,SBGS_VD_ro(theSBGS,bl),NEWDEF_CLASS,
                        SBGS_MD_Ao(theSBGS,bl),
                        SBGS_VD_cd(theSBGS),ACTIVE_CLASS)) NP_RETURN (1,result[0]);
  }

  return (0);
}

static INT SBGSPostProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b,
                            MATDATA_DESC *A, INT *result)
{
  NP_SBGS *theSBGS;
  INT bl;

  theSBGS = (NP_SBGS*) theNP;

  for (bl=0; bl<SBGS_NBLOCKS(theSBGS); bl++)
    if ((*SBGS_BLOCKITER(theSBGS,bl)->iter.PostProcess!=NULL))
      if ((*SBGS_BLOCKITER(theSBGS,bl)->iter.PostProcess)
            (&SBGS_BLOCKITER(theSBGS,bl)->iter,level,x,b,A,result))
        NP_RETURN (1,result[0]);

  return(0);
}

static INT SBGSConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = SBGS_Init;
  theNP->Display = SBGS_Display;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = SBGSPreProcess;
  np->Iter = SBGSSmoother;
  np->PostProcess = SBGSPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   gbggs - numproc for grid (or geometric) block Gauss-Seidel smoother

   DESCRIPTION:
   This numproc executes an grid (or geometric) block Gauss-Seidel smoother,
   using the blas routine
   'l_lrdecompB' and 'l_lgsB'. It can be used in 'lmgc'.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
              [$damp <sc double list>] [$L <mat>]
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $L~<mat> - decomposition matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT GBGSPreProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           INT *baselevel, INT *result)
{
  NP_SMOOTHER *np;
  GRID *theGrid;

  np = (NP_SMOOTHER*) theNP;
  theGrid = NP_GRID(theNP,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_lrdecompB(theGrid,np->L)!=NUM_OK) {
    PrintErrorMessage('E',"GBGSPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
  }
  *baselevel = level;

  return (0);
}

static INT GBGSStep (NP_SMOOTHER *theNP, INT level,
                     VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                     MATDATA_DESC *L,
                     INT *result)
{
    #ifdef ModelP
  if (l_vector_collect(NP_GRID(theNP,level),b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_lgsB(NP_GRID(theNP,level),x,L,b)) NP_RETURN(1,result[0]);

  return (0);
}

static INT GBGSConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SmootherInit;
  theNP->Display = SmootherDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = GBGSPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = GBGSStep;

  return(0);
}

/****************************************************************************/
/*D
   ilu - numproc for point block beta-modified ilu smoother

   DESCRIPTION:
   This numproc executes a point block ilu smoother, using the blas routines
   'l_ilubthdecomp' and 'l_luiter'. It can be used in 'lmgc'.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$L <mat>]
       [$damp <sc double list>] [$beta <sc double list>]
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $L~<mat> - decomposition matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
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
  for (i=0; i<MAX_VEC_COMP; i++) np->mindiag[i] = 0.0;
  sc_read(np->mindiag,np->smoother.iter.b,"mindiag",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT ILUDisplay (NP_BASE *theNP)
{
  NP_ILU *np;

  SmootherDisplay(theNP);
  np = (NP_ILU *) theNP;
  if (sc_disp(np->beta,np->smoother.iter.b,"beta"))
    REP_ERR_RETURN (1);
  if (sc_disp(np->mindiag,np->smoother.iter.b,"mindiag"))
    REP_ERR_RETURN (1);

  return (0);
}

static INT ILUPreProcess (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                          INT *baselevel, INT *result)
{
  NP_ILU *np;
  GRID *theGrid;

  np = (NP_ILU *) theNP;
  theGrid = NP_GRID(theNP,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  /*	if (np->mindiag[0] > 0.0)
      if (l_shift_diagonal(theGrid,np->smoother.L,np->mindiag) != NUM_OK)
              NP_RETURN(1,result[0]); */
  if (l_ilubthdecomp(theGrid,np->smoother.L,np->beta,NULL,NULL,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"ILUPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
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
  if (l_vector_collect(NP_GRID(theNP,level),b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_luiter(NP_GRID(theNP,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

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
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$L <mat>]
       [$damp <sc double list>] [$beta <sc double list>]
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $L~<mat> - decomposition matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
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
  theGrid = NP_GRID(theNP,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_ilubthdecomp_fine(theGrid,np->smoother.L,np->beta,NULL,NULL,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"FILUPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
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
  if (l_vector_collect(NP_GRID(theNP,level),b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_luiter_fine(NP_GRID(theNP,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

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
   thilu - numproc for point block beta-modified ilu smoother with threshold for extending the sparsity pattern

   DESCRIPTION:
   This numproc executes a point block ilu smoother, using the blas routines
   'l_ilubthdecomp' and 'l_luiter'. It can be used in 'lmgc'.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$L <mat>]
       [$damp <sc double list>] [$beta <sc double list>]
           [$thresh <sc double list>];
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $L~<mat> - decomposition matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal
   .  $thresh~<sc~double~list> - parameter for the threshold

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
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
  if (sc_disp(np->beta,np->smoother.iter.b,"beta")) REP_ERR_RETURN (1);
  if (sc_disp(np->thresh,np->smoother.iter.b,"thresh")) REP_ERR_RETURN (1);

  return (0);
}

static INT THILUPreProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *baselevel, INT *result)
{
  NP_THILU *np;
  GRID *theGrid;

  np = (NP_THILU *) theNP;
  theGrid = NP_GRID(theNP,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_ilubthdecomp(theGrid,np->smoother.L,np->beta,np->thresh,NULL,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"THILUPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
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
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$L <mat>]
       [$damp <sc double list>] [$beta <sc double list>]
           $mode {local|global};
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $L~<mat> - decomposition matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  $beta~<sc~double~list> - parameter for modification of the diagonal
   .  $mode~{local|global} - global or local shift

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
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
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  if (strncmp(buffer,"global",3)==0)
    np->mode = SP_GLOBAL;
  else if (strncmp(buffer,"local",3)==0)
    np->mode = SP_LOCAL;
  else
  {
    PrintErrorMessage('E',"SPILUInit","specify local/global for mode");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }

  return (SmootherInit(theNP,argc,argv));
}

static INT SPILUDisplay (NP_BASE *theNP)
{
  NP_SPILU *np;

  SmootherDisplay(theNP);
  np = (NP_SPILU *) theNP;
  if (sc_disp(np->beta,np->smoother.iter.b,"beta")) REP_ERR_RETURN (1);
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
  VECDATA_DESC *tmp=NULL;

  np = (NP_SPILU *) theNP;
  theGrid = NP_GRID(theNP,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theNP->base.mg,level,level,x,&tmp)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_iluspdecomp(theGrid,np->smoother.L,np->beta,tmp,np->mode,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"SPILUPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
  }
  *baselevel = level;

  FreeVD(theNP->base.mg,level,level,tmp);

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
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$L <mat>]
       [$damp <sc double list>];
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $L~<mat> - decomposition matrix
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
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
  theGrid = NP_GRID(theNP,level);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if (l_dmatcopy(theGrid,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
        #endif
  if (l_icdecomp(theGrid,np->smoother.L)
      !=NUM_OK) {
    PrintErrorMessage('E',"ICPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
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
  if (l_vector_collect(NP_GRID(theNP,level),b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (l_lltiter(NP_GRID(theNP,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

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
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$L <mat>]
       [$damp <sc double list>];
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $L~<mat> - decomposition matrix
   .  $damp~<sc~double~list> - damping factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
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

  np = (NP_SMOOTHER *) theNP;

  theGrid = NP_GRID(theNP,level);
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
        NP_RETURN(1,result[0]);
      default :
        PrintErrorMessage('E',"LUPreProcess","err > 0");
        NP_RETURN(1,result[0]);
      }
    }
    if (err!=-VINDEX(LASTVECTOR(theGrid))) {
      PrintErrorMessageF('E',"LUPreProcess","decomp failed: IDX %ld on level %d",
                         -err,GLEVEL(theGrid));
      UserWriteF(" - LASTVECTOR has IDX %ld\n",
                 VINDEX(LASTVECTOR(theGrid)));
      NP_RETURN(1,result[0]);
    }
    if (l_lrregularize(theGrid,np->L) !=NUM_OK) {
      PrintErrorMessage('E',"LUPreProcess","cannot regularize");
      NP_RETURN(1,result[0]);
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
   ff - numproc for frequency filtering solvers

   DESCRIPTION:
   This numproc solves an equation with the tangential frequency filtering
   method.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>] [$damp <sc double list>]
       [$FF <FF-mat sym>] [$FF3D <3D FF-mat sym>] [$L <LU mat sym>]
       [$tv <testvector sym>] [$tv2 <2. testvector sym>]
           [$t <update for correction sym>]
           $display {no|red|full} $wr <"all"|number> $wr3D <number> $type <FF|TFF>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  $FF~<FF-mat~sym> - symbol for the frequency filtered matrix
   .  $FF3D~<3D~FF-mat~sym> - symbol for an additional frequency filtered matrix for 3D
   .  $L~<LU-mat~sym> - symbol for the LU decomposed matrix
   .  $tv~<testvector~sym> - symbol for the testvector
   .  $tv2~<2.~testvector~sym> - symbol for the second testvector if neccessary
   .  $t~<update~for~correction~sym> - temp. vector
   .  $type~<type of frequency filter> - "TFF" for Wagners or "FF" for Wittums
   .  $display - display mode: 'no', 'red'uced or 'full'
   .  $wr - relative frequency [0..1] for 2D OR 'all' for the whole logarithmic sequence of frequencies
   .  $wr3D - relative frequency [0..1] for 3D

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess

   EXAMPLE:
   FF as smoother:

   .vb
        npcreate smooth $c ff;
        npinit smooth $wr 0.5 $wr3D 0.5 $type TFF $display full;

        npcreate ls_iter $c lmgc;
        npinit ls_iter $S smooth smooth basesolver $T transfer
                        $b @BASELEVEL $n1 1 $n2 1 $g 1;
   .ve

   FF as solver:

   .vb
        npcreate ls_iter $c ff;
        npinit ls_iter $wr ALL $wr3D -1.0 $type TFF $display full;

        npcreate mgs $c ls;
        npinit mgs $A MAT $x sol $b rhs $m 8 $abslimit 1e-8 $red 1e-30
                $I ls_iter $display full;

   .ve
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   FFInit - Init tangential frequency filtering iterator numproc

   SYNOPSIS:
   static INT FFInit (NP_BASE *theNP, INT argc , char **argv);

   PARAMETERS:
   .  theNP - pointer to numproc
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function inits the numerical procedure for
   tangential frequency filtering iterator. The data descriptors,
   the display mode and the blockvector description format are set.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT FFInit (NP_BASE *theNP, INT argc , char **argv)
{
#ifdef __BLOCK_VECTOR_DESC__

  NP_FF *np;
  char buffer[128];
  MULTIGRID *theMG;
  INT i;

  np = (NP_FF *) theNP;
  theMG = np->smoother.iter.base.mg;

  TOS_FF_Vecs = 0;

  for( i=0; i < FF_MAX_VECS; i++ )
  {
    FF_Vecs[i] = DUMMY_COMP;
    FF_VECDATA_DESC_ARRAY[i] = NULL;
  }

  for( i=0; i < FF_MAX_MATS; i++ )
  {
    FF_Mats[i] = DUMMY_COMP;
    FF_MATDATA_DESC_ARRAY[i] = NULL;
  }

#ifdef __THREEDIM__
  if ( ReadArgvChar ( "wr3D", buffer, argc, argv) )
  {
    PrintErrorMessage('E',"FFInit", "Option $wr3D mandatory");
    REP_ERR_RETURN(1);
  }
  sscanf(buffer,"%lf", &NPFF_WaveNrRel3D(np) );
#else
  NPFF_WaveNrRel3D(np) = -1.0;
#endif

  NPFF_tv(np)  = ReadArgvVecDesc(theMG,"tv",argc,argv);
  NPFF_tv2(np)  = ReadArgvVecDesc(theMG,"tv2",argc,argv);
  NPFF_t(np)   = ReadArgvVecDesc(theMG,"t",argc,argv);

  NPFF_DISPLAY(np) = ReadArgvDisplay(argc,argv);
  NPFF_MESHWIDTH(np) = 0.0;

  if ( ReadArgvChar ( "wr", buffer, argc, argv) )
  {
    PrintErrorMessage('E',"FFInit", "Option $wr mandatory");
    REP_ERR_RETURN(1);
  }
  if( strcmp( buffer, "ALL") == 0 || strcmp( buffer, "all") == 0 )
  {
    NPFF_ALLFREQ(np) = TRUE;
    NPFF_WaveNrRel(np) = -1.0;
  }
  else
  {
    NPFF_ALLFREQ(np) = FALSE;
    sscanf(buffer,"%lf", &NPFF_WaveNrRel(np) );
  }

  if ( ReadArgvChar ( "type", buffer, argc, argv) )
  {
    PrintErrorMessage('W',"FFInit", "default type TFF set");
    NPFF_TYPE(np) = TYPE_TFF;
  }
  else
  {
    if( strcmp( buffer, "TFF") == 0 )
      NPFF_TYPE(np) = TYPE_TFF;
    else if( strcmp( buffer, "FF") == 0 )
      NPFF_TYPE(np) = TYPE_FF;
    else
    {
      PrintErrorMessage('E',"FFInit", "Option $type: wrong argument");
      REP_ERR_RETURN(1);
    }
  }

#ifdef __TWODIM__
  *NPFF_BVDF(np) = two_level_bvdf;
#else
  *NPFF_BVDF(np) = three_level_bvdf;
#endif

  /* reset other parameters */
  NPFF_MESHWIDTH(np) = 0.0;

  return (SmootherInit(theNP,argc,argv));

#else
  PrintErrorMessage( 'E', "FFInit", "__BLOCK_VECTOR_DESC__ must be defined in gm.h" );
  REP_ERR_RETURN(1);
#endif /* __BLOCK_VECTOR_DESC__ */
}

/****************************************************************************/
/*
   FFDisplay - Display tangential frequency filtering iterator numproc

   SYNOPSIS:
   static INT FFDisplay (NP_BASE *theNP);

   PARAMETERS:
   .  theNP - pointer to numproc

   DESCRIPTION:
   This function displays the parameters set for the tangential frequency
   filtering iterator numerical procedure.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT FFDisplay (NP_BASE *theNP)
{
  NP_FF *np;
  INT i;

  SmootherDisplay(theNP);
  np = (NP_FF *) theNP;

  UserWrite("FF specific data:\n");

  if (NPFF_tv(np) != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"tv",ENVITEM_NAME(NPFF_tv(np)));
  if (NPFF_tv2(np) != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"tv2",ENVITEM_NAME(NPFF_tv2(np)));

  UserWrite("matrix hierarchy:");
  i = 0;
  while ( (FF_Mats[i] != DUMMY_COMP) && (i < FF_MAX_MATS) )
    UserWriteF("  %d", FF_Mats[i++] );
  UserWrite("\naux vector list:");
  i = 0;
  while ( (FF_Vecs[i] != DUMMY_COMP) && (i < FF_MAX_VECS) )
    UserWriteF("  %d", FF_Vecs[i++] );
  UserWrite("\n");

  UserWriteF(DISPLAY_NP_FORMAT_SF,"meshwidth",(double)NPFF_MESHWIDTH(np));

  if ( NPFF_ALLFREQ(np) == TRUE )
    UserWriteF(DISPLAY_NP_FORMAT_SS,"frequency","ALL");
  else
  {
                #ifdef __THREEDIM__
    UserWriteF(DISPLAY_NP_FORMAT_SF,"frequency (2D)",(double)NPFF_WaveNrRel(np));
    UserWriteF(DISPLAY_NP_FORMAT_SF,"frequency (3D)",(double)NPFF_WaveNrRel3D(np));
                #else
    UserWriteF(DISPLAY_NP_FORMAT_SF,"frequency",(double)NPFF_WaveNrRel(np));
                #endif
  }

  if (NPFF_DO_TFF(np) )
    UserWriteF(DISPLAY_NP_FORMAT_SS,"type","tangential FF (Wagner)");
  if (NPFF_DO_FF(np) )
    UserWriteF(DISPLAY_NP_FORMAT_SS,"type","ordinary FF (Wittum)");

  if (NPFF_DISPLAY(np) == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (NPFF_DISPLAY(np) == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (NPFF_DISPLAY(np) == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  return (0);
}

/****************************************************************************/
/*D
   FFPreProcess - Prepare tangential frequency filtering solver

   SYNOPSIS:
   static INT FFPreProcess (NP_ITER *theNP, INT level,
   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
   INT *baselevel, INT *result);

   PARAMETERS:
   .  theNP - pointer to numproc
   .  level - gridlevel to be prepared
   .  x - solution vector
   .  b - defect vector
   .  A - stiffness matrix
   .  baselevel - output: baselevel used by iter (== level)
   .  result - return value of the function

   DESCRIPTION:
   This function prepares a tangential frequency filtering iteration:
   allocate temporarily the neccessary data descriptors,
   determine the meshwidth of the grid, construct the linewise (and in
   3D additional planewise) blockvector decomposition, puts the dirichlet
   values on the right hand side and disposes all connections
   consisting entirely of matrixvalues 0, if only one testfrequency should be
   considered calculate the FF decomposition of the stiffnes matrix in
   smoothers matrix L.

   Points must be ordered lexicographic, boundary nodes at the end of the
   list. The grid must be a square.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT FFPreProcess (NP_ITER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *baselevel, INT *result)
{
#ifdef __BLOCK_VECTOR_DESC__
  NP_FF *np;
  GRID *theGrid;
  DOUBLE wavenr, wavenr3D, meshwidth;
  INT i, n;
  BV_DESC bvd;

  np = (NP_FF *) theNP;
  theGrid = NP_GRID(theNP,level);

  BVD_INIT( &bvd );
  BVD_PUSH_ENTRY( &bvd, 0, NPFF_BVDF(np) );

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  if (AllocMDFromMD(theNP->base.mg,level,level,A,&np->smoother.L))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theNP->base.mg,level,level,x,&NPFF_tv(np)))
    NP_RETURN(1,result[0]);

  if (NPFF_DO_FF(np))
  {
    if (AllocVDFromVD(theNP->base.mg,level,level,x,&NPFF_tv2(np)))
      NP_RETURN(1,result[0]);
  }

  /* check if all objects are valid and scalar */
  if ( A == NULL )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol A is not defined" );
    NP_RETURN(1,result[0]);
  }
  if ( !MD_IS_SCALAR(A) )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol A is not scalar" );
    NP_RETURN(1,result[0]);
  }


  if ( np->smoother.L == NULL )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol L is not defined" );
    NP_RETURN(1,result[0]);
  }
  if ( !MD_IS_SCALAR( np->smoother.L ) )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol L is not scalar" );
    NP_RETURN(1,result[0]);
  }


  if ( x == NULL )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol x is not defined" );
    NP_RETURN(1,result[0]);
  }
  if ( !VD_IS_SCALAR( x ) )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol x is not scalar" );
    NP_RETURN(1,result[0]);
  }

  if ( b == NULL )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol b is not defined" );
    NP_RETURN(1,result[0]);
  }
  if ( !VD_IS_SCALAR( b ) )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol b is not scalar" );
    NP_RETURN(1,result[0]);
  }

  if ( NPFF_tv(np) == NULL )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol tv is not defined" );
    NP_RETURN(1,result[0]);
  }
  if ( !VD_IS_SCALAR( NPFF_tv(np) ) )
  {
    PrintErrorMessage( 'E', "FFPreProcess", "Symbol tv is not scalar" );
    NP_RETURN(1,result[0]);
  }

  if (NPFF_DO_FF(np))
  {
    if ( NPFF_tv2(np) == NULL )
    {
      PrintErrorMessage( 'E', "FFPreProcess", "Symbol tv2 is not defined" );
      NP_RETURN(1,result[0]);
    }
    if ( !VD_IS_SCALAR( NPFF_tv2(np) ) )
    {
      PrintErrorMessage( 'E', "FFPreProcess", "Symbol tv2 is not scalar" );
      NP_RETURN(1,result[0]);
    }
  }

  /*	if ( FF_Mats[0] == DUMMY_COMP )
     {*/
  /* init the global arrays once ! */

#ifdef __THREEDIM__
  n = 2;
#else
  n = 1;
#endif

  /* first component is always the original, global stiffness matrix */
  if( FF_Mats[0] == DUMMY_COMP )
    FF_Mats[0] = MD_SCALCMP( A );
  else
  {
    ASSERT(FF_Mats[0]==MD_SCALCMP( A ));
  }

  for( i=1; i<=n; i++ )
  {
    if (AllocMDFromMD(theNP->base.mg,level,level,A,FF_MATDATA_DESC_ARRAY+i))
      NP_RETURN(1,result[0]);
    if ( FF_Mats[i] == DUMMY_COMP )
      FF_Mats[i] = MD_SCALCMP( FF_MATDATA_DESC_ARRAY[i] );
    else
    {
      ASSERT(FF_Mats[i]==MD_SCALCMP( FF_MATDATA_DESC_ARRAY[i] ));
    }
  }

  /* in the last component there is stored the LU decomp of elementary blocks */
  if( FF_Mats[i] == DUMMY_COMP )
    FF_Mats[i] = MD_SCALCMP( np->smoother.L );
  else
  {
    ASSERT(FF_Mats[i]==MD_SCALCMP( np->smoother.L ));
  }

#ifdef __THREEDIM__
  n = 2;
#else
  n = 1;
#endif

  if ( NPFF_DO_FF(np) )
    n *= 2;

  for( i=0; i<n; i++ )
  {
    if (AllocVDFromVD(theNP->base.mg,level,level,x,FF_VECDATA_DESC_ARRAY+i))
      NP_RETURN(1,result[0]);
    FF_Vecs[i] = VD_SCALCMP( FF_VECDATA_DESC_ARRAY[i] );
  }
  /*	}*/


  if (FF_PrepareGrid( theGrid, &meshwidth, TRUE, MD_SCALCMP( A ), VD_SCALCMP( x ), VD_SCALCMP( b ), NPFF_BVDF(np) )!=NUM_OK)
  {
    PrintErrorMessage('E',"FFPreProcess","preparation of the grid failed");
    NP_RETURN(1,result[0]);
  }
  NPFF_MESHWIDTH(np) = meshwidth;

  if ( !NPFF_ALLFREQ(np) )
  {
    n = (INT)( log(1.0/meshwidth)/M_LN2 + 0.5 );
    wavenr = (DOUBLE)(1<<(INT)( (n-1) * NPFF_WaveNrRel(np) + 0.5 ));
    wavenr3D = (DOUBLE)(1<<(INT)( (n-1) * NPFF_WaveNrRel3D(np) + 0.5 ));

    if (NPFF_DO_TFF(np) && TFFDecomp( wavenr, wavenr3D, GFIRSTBV(theGrid), &bvd,
                                      NPFF_BVDF(np),
                                      VD_SCALCMP( NPFF_tv(np) ),
                                      theGrid ) != NUM_OK )
    {
      PrintErrorMessage('E',"FFPreProcess","decomposition failed");
      NP_RETURN(1,result[0]);
    }
    else if (NPFF_DO_FF(np) && FFDecomp( wavenr, wavenr3D, GFIRSTBV(theGrid), &bvd,
                                         NPFF_BVDF(np),
                                         VD_SCALCMP( NPFF_tv(np) ),
                                         VD_SCALCMP( NPFF_tv2(np) ),
                                         theGrid ) != NUM_OK )
    {
      PrintErrorMessage('E',"FFPreProcess","decomposition failed");
      NP_RETURN(1,result[0]);
    }
  }

  *baselevel = level;

  return (0);

#else
  PrintErrorMessage( 'E', "FFPreProcess", "__BLOCK_VECTOR_DESC__ must be defined in gm.h" );
  REP_ERR_RETURN (1);
#endif
}

/****************************************************************************/
/*D
   FFPostProcess - Prepare tangential frequency filtering solver

   SYNOPSIS:
   static INT FFPostProcess (NP_ITER *theNP, INT level,
   VECDATA_DESC *x, VECDATA_DESC *b,
   MATDATA_DESC *A, INT *result);

   PARAMETERS:
   .  theNP - pointer to numproc
   .  level - gridlevel to be postprocessed
   .  x - solution vector
   .  b - defect vector
   .  A - stiffness matrix
   .  result - return value of the function

   DESCRIPTION:
   This function postprocesses a tangential frequency filtering iteration:
   Free all temporarily allocated data descriptors, free all 'BLOCKVECTOR's
   in the grid and rebuild the 0-connections freed in the preprocess.
   Then proceed with SmootherPostProcess.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT FFPostProcess (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b,
                          MATDATA_DESC *A, INT *result)
{
  NP_FF *np;
  MULTIGRID *theMG;
  INT i;

  np = (NP_FF *) theNP;
  theMG = np->smoother.iter.base.mg;

  if (NPFF_tv(np) != NULL)
    FreeVD(theMG,level,level,NPFF_tv(np));
  if (NPFF_tv2(np) != NULL)
    FreeVD(theMG,level,level,NPFF_tv2(np));

  for( i=1; i<FF_MAX_MATS; i++ )
  {
    if ( FF_MATDATA_DESC_ARRAY[i] != NULL )
    {
      FreeMD(theMG,level,level,FF_MATDATA_DESC_ARRAY[i]);
      FF_MATDATA_DESC_ARRAY[i] = NULL;
    }
    FF_Mats[i] = DUMMY_COMP;
  }


  for( i=0; i<FF_MAX_VECS; i++ )
  {
    if ( FF_Vecs[i] != DUMMY_COMP )
    {
      FreeVD(theMG,level,level,FF_VECDATA_DESC_ARRAY[i]);
      FF_VECDATA_DESC_ARRAY[i] = NULL;
      FF_Vecs[i] = DUMMY_COMP;
    }
  }

  FreeAllBV( GRID_ON_LEVEL(theMG,level) );
  if (MGCreateConnection(theMG))        /* restore the disposed connections */
  {
    PrintErrorMessage('E',"FFPostProcess","MGCreateConnection failed");
    NP_RETURN(1,result[0]);
  }
  return (SmootherPostProcess (theNP, level, x, b, A, result));
}


/****************************************************************************/
/*D
   FFIter - Perform one tangential frequency filtering iteration

   SYNOPSIS:
   static INT FFIter (NP_ITER *theNP, INT level,
                                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                                        INT *result);

   PARAMETERS:
   .  theNP - pointer to numproc
   .  level - gridlevel to be postprocessed
   .  x - correction or correction-update vector
   .  b - defect vector
   .  A - stiffness matrix
   .  result - return value of the function

   DESCRIPTION:
   This function performs one tangential frequency filtering iteration.
   It makes a difference using it as solver or as smoothing iteration.

   As smoothing iteration:
   Give at the '$wr' (and $wr3D in 3D) a certain relatively wavenumber in the
   range 0..1. Then the defect 'b' will be updated and the correction-update
   is returned in 'x'.

   As solver iteration:
   Give 'ALL' as '$wr' argument. For each absolute wavenumber 1..(1/h)/2
   a FF step is performed and the defect 'b' will be updated and the
   correction is returned in 'x'. According to '$display' option after
   each FF step the defect and convergence rate is printed.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT FFIter (NP_ITER *theNP, INT level,
                   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                   INT *result)
{
#ifdef __BLOCK_VECTOR_DESC__
  NP_FF *np;
  BV_DESC bvd;
  GRID *theGrid;
  DOUBLE end_wave, wavenr, start_norm, new_norm;

  np = (NP_FF *) theNP;
  BVD_INIT( &bvd );
  BVD_PUSH_ENTRY( &bvd, 0, NPFF_BVDF(np) );

  theGrid = NP_GRID(theNP,level);

  /* make a copy for displaying */
  np->smoother.iter.c = x;

  if ( !NPFF_ALLFREQ(np) )
  {             /* smooth only for 1 testvector frequency */
                /* copy defect to tv because FFMultWithMInv destroys its defect */
    dcopyBS( GFIRSTBV(theGrid), VD_SCALCMP( NPFF_tv(np) ), VD_SCALCMP( b ) );
    if (FFMultWithMInv( GFIRSTBV(theGrid), &bvd, NPFF_BVDF(np),
                        VD_SCALCMP( x ),
                        VD_SCALCMP( NPFF_tv(np) ) ) != NUM_OK)
    {
      PrintErrorMessage('E',"FFStep","inversion failed");
      NP_RETURN(1,result[0]);
    }
    /* defect -= A * corr_update */
    dmatmul_minusBS( GFIRSTBV(theGrid), &bvd, NPFF_BVDF(np),
                     VD_SCALCMP( b ), MD_SCALCMP( A ), VD_SCALCMP( x ));
  }
  else
  {             /* smooth for all testvector frequencies */

    /* alloc temp. for correction update (in x!) */
    if (AllocVDFromVD(theNP->base.mg,level,level,x,&NPFF_t(np))) NP_RETURN(1,result[0]);

    if ( NPFF_DISPLAY(np) != PCR_NO_DISPLAY )
      if(eunormBS( GFIRSTBV(theGrid), VD_SCALCMP( b ), &new_norm ) ) NP_RETURN(1,result[0]);

    end_wave = 1.0 / NPFF_MESHWIDTH(np) - 0.5;             /* rounding */
    for ( wavenr = 1.0; wavenr < end_wave; wavenr *= 2.0 )
    {                   /* wave 1.0 ... (1/h)/2 */
      if (NPFF_DO_TFF(np) )
      {
        if ( TFFDecomp( wavenr, wavenr, GFIRSTBV(theGrid), &bvd,
                        NPFF_BVDF(np),
                        VD_SCALCMP( NPFF_tv(np) ),
                        theGrid ) != NUM_OK )
        {
          PrintErrorMessage('E',"FFStep","TFF decomposition failed");
          NP_RETURN(1,result[0]);
        }
      }
      else if (NPFF_DO_FF(np))
      {
        /*if (wavenr == 2.0) wavenr = 3.0;*/ /* wavenr==2 already in the last step */
        printf("wavenr %g\n", wavenr);
        if (FFDecomp( wavenr, wavenr, GFIRSTBV(theGrid), &bvd,
                      NPFF_BVDF(np),
                      VD_SCALCMP( NPFF_tv(np) ),
                      VD_SCALCMP( NPFF_tv2(np) ),
                      theGrid ) != NUM_OK )
        {
          PrintErrorMessage('E',"FFStep","FF decomposition failed");
          NP_RETURN(1,result[0]);
        }
      }

      /* copy defect to aux because FFMultWithMInv destroys its defect */
      dcopyBS( GFIRSTBV(theGrid), VD_SCALCMP( NPFF_t(np) ), VD_SCALCMP( b ) );
      if (FFMultWithMInv( GFIRSTBV(theGrid), &bvd, NPFF_BVDF(np),
                          VD_SCALCMP( NPFF_t(np) ),
                          VD_SCALCMP( NPFF_t(np) ) ) != NUM_OK)
      {
        PrintErrorMessage('E',"FFStep","inversion failed");
        NP_RETURN(1,result[0]);
      }

      /* corr += corr_update */
      daddBS( GFIRSTBV(theGrid), VD_SCALCMP( x ), VD_SCALCMP( NPFF_t(np) ) );

      /* defect -= A * corr_update */
      dmatmul_minusBS( GFIRSTBV(theGrid), &bvd, NPFF_BVDF(np),
                       VD_SCALCMP( b ), MD_SCALCMP( A ), VD_SCALCMP( NPFF_t(np) ));

      if ( NPFF_DISPLAY(np) != PCR_NO_DISPLAY )
      {
        start_norm = new_norm;
        if(eunormBS( GFIRSTBV(theGrid), VD_SCALCMP( b ), &new_norm ) ) NP_RETURN(1,result[0]);

        UserWriteF( "Wnr plane = %4g Wnr line = %4g new defect = %12lg "
                    "conv. rate = %12lg\n", wavenr, wavenr, new_norm,
                    new_norm/start_norm );
      }
    }

    FreeVD(theNP->base.mg,level,level,NPFF_t(np));
  }

  return (0);
#else
  PrintErrorMessage( 'E', "FFStep", "__BLOCK_VECTOR_DESC__ must be defined in gm.h" );
  REP_ERR_RETURN (1);
#endif
}

static INT FFConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = FFInit;
  theNP->Display = FFDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = FFPreProcess;
  np->iter.Iter = FFIter;
  np->iter.PostProcess = FFPostProcess;
  np->Step = NULL;

  return(0);
}

/****************************************************************************/
/*D
   lmgc - numproc for linear multigrid cycle

   DESCRIPTION:
   This numproc executes

   .vb
   npinit <name> [$c <cor>] [$r <rhs>] [$A <mat>]
       $S <pre post base> $T <transfer>
       [$b <baselevel>] [$g <gamma>] [$n1 <it>] [$n2 <it>]
   .ve

   .  $c~<cor> - correction vector
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

   SEE ALSO:
   ls
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

  if (np->Transfer == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PreSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PostSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->BaseSolver == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);

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
      REP_ERR_RETURN(1);

  if (np->PreSmooth->PreProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PreProcess)
            (np->PreSmooth,i,x,b,A,baselevel,result))
        REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PreProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PreSmooth->PreProcess)
              (np->PostSmooth,i,x,b,A,baselevel,result))
          REP_ERR_RETURN(1);

  *baselevel = MIN(np->baselevel,level);
  if (np->BaseSolver->PreProcess != NULL)
    if ((*np->BaseSolver->PreProcess)
          (np->BaseSolver,*baselevel,x,b,A,baselevel,result))
      REP_ERR_RETURN(1);

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
      REP_ERR_RETURN(1);
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
      REP_ERR_RETURN(1);
    if (l_daxpy(theGrid,c,ACTIVE_CLASS,Factor_One,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }
  if ((*np->Transfer->RestrictDefect)
        (np->Transfer,level,b,b,A,Factor_One,result))
    REP_ERR_RETURN(1);

  if (l_dset(DOWNGRID(theGrid),c,EVERY_CLASS,0.0) != NUM_OK) NP_RETURN(1,result[0]);
  for (i=0; i<np->gamma; i++)
    if (Lmgc(theNP,level-1,c,b,A,result))
      REP_ERR_RETURN(1);
  if ((*np->Transfer->InterpolateCorrection)
        (np->Transfer,level,np->t,c,A,Factor_One,result))
    REP_ERR_RETURN(1);
  if (l_daxpy(theGrid,c,EVERY_CLASS,Factor_One,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,np->t,EVERY_CLASS) != NUM_OK) NP_RETURN(1,result[0]);
  for (i=0; i<np->nu2; i++) {
    if ((*np->PostSmooth->Iter)(np->PostSmooth,level,np->t,b,A,result))
      REP_ERR_RETURN(1);
    if (l_daxpy(theGrid,c,ACTIVE_CLASS,Factor_One,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }
  FreeVD(theNP->base.mg,level,level,np->t);
  if (np->Transfer->AdaptCorrection != NULL)
    if ((*np->Transfer->AdaptCorrection)(np->Transfer,level,c,b,A,result))
      REP_ERR_RETURN(1);

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
    if ((*np->Transfer->PostProcess)
          (np->Transfer,np->baselevel,level,x,b,A,result))
      REP_ERR_RETURN(1);

  if (np->PreSmooth->PostProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PostProcess)
            (np->PreSmooth,i,x,b,A,result))
        REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PostProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PreSmooth->PostProcess)
              (np->PostSmooth,i,x,b,A,result))
          REP_ERR_RETURN(1);

  if (np->BaseSolver->PostProcess != NULL)
    if ((*np->BaseSolver->PostProcess)
          (np->BaseSolver,np->baselevel,x,b,A,result))
      REP_ERR_RETURN(1);

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
/*D
   addmgc - numproc for additive linear multigrid cycle

   DESCRIPTION:
   This numproc executes

   .vb
   npinit <name> [$c <cor>] [$r <rhs>] [$A <mat>]
       $S <smooth> $T <transfer>
       [$b <baselevel>] [$g <gamma>] [$n1 <it>] [$n2 <it>]
   .ve

   .  $c~<cor> - correction vector
   .  $r~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $T~<transfer> - transfer numproc
   .  $S~<smooth> - numproc for smoother
   .  $b~<baselevel> - baselevel where the base solver is called
   .  $g~<gamma> - number of iterations of Lmgc per level (default gamma = 1)
   .  $n1~<it> - number of iterations of the presmoother (default n1 = 1)
   .  $n2~<it> - number of iteration of the postsmoother (default n2 = 1)

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT AddmgcInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LMGC *np;
  INT i;
  char pre[VALUELEN];

  np = (NP_LMGC *) theNP;

  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->Transfer = (NP_TRANSFER *)
                 ReadArgvNumProc(theNP->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  for (i=1; i<argc; i++)
    if (argv[i][0]=='S') {
      if (sscanf(argv[i],"S %s",pre)!=1)
        continue;
      np->PreSmooth = (NP_ITER *)
                      GetNumProcByName(theNP->mg,pre,ITER_CLASS_NAME);
      break;
    }
  if (ReadArgvINT("n1",&(np->nu1),argc,argv))
    np->nu1 = 1;
  if (ReadArgvINT("n2",&(np->nu2),argc,argv))
    np->nu2 = 0;
  np->nu1 += np->nu2;
  if (ReadArgvINT("b",&(np->baselevel),argc,argv))
    np->baselevel = 0;

  if (np->Transfer == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PreSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);

  return (NPIterInit(&np->iter,argc,argv));
}

static INT AddmgcDisplay (NP_BASE *theNP)
{
  NP_LMGC *np;

  np = (NP_LMGC *) theNP;

  NPIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)np->nu1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"b",(int)np->baselevel);

  if (np->Transfer != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(np->Transfer));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (np->PreSmooth != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"smooth",ENVITEM_NAME(np->PreSmooth));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"smooth","---");

  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));

  return (0);
}

static INT AddmgcPreProcess  (NP_ITER *theNP, INT level,
                              VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                              INT *baselevel, INT *result)
{
  NP_LMGC *np;
  INT i;

  np = (NP_LMGC *) theNP;

  if (np->Transfer->PreProcess != NULL)
    if ((*np->Transfer->PreProcess)
          (np->Transfer,np->baselevel,level,x,b,A,result))
      REP_ERR_RETURN(1);

  if (np->PreSmooth->PreProcess != NULL)
    for (i = np->baselevel; i <= level; i++)
      if ((*np->PreSmooth->PreProcess)
            (np->PreSmooth,i,x,b,A,baselevel,result))
        REP_ERR_RETURN(1);

  *baselevel = MIN(np->baselevel,level);

  return (0);
}

static INT Addmgc (NP_ITER *theNP, INT level,
                   VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A,
                   INT *result)
{
  NP_LMGC *np;
  MULTIGRID *theMG;
  GRID *theGrid;
  LRESULT lresult;
  INT i;
  INT mylevel;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = c;
  NPIT_b(theNP) = b;

  np = (NP_LMGC *) theNP;
  theMG = theNP->base.mg;
  for (mylevel=level; mylevel > np->baselevel ; mylevel--)
    if ((*np->Transfer->RestrictDefect)
          (np->Transfer,mylevel,b,b,A,Factor_One,result))
      REP_ERR_RETURN(1);
    #ifdef ModelP
  if (a_vector_collect(theMG,np->baselevel,level,b)!=NUM_OK)
    NP_RETURN(1,result[0]);
  /* TODO: supress collect and consistent routines in smoother */
    #endif
  for (mylevel=np->baselevel; mylevel < level ; mylevel++) ;{
    theGrid = GRID_ON_LEVEL(theMG,mylevel);
    if (AllocVDFromVD(theMG,mylevel,mylevel,c,&np->t))
      NP_RETURN(1,result[0]);
    for (i=0; i<np->nu1; i++) {
      if ((*np->PreSmooth->Iter)
            (np->PreSmooth,mylevel,np->t,b,A,result))
        REP_ERR_RETURN(1);
      if (l_daxpy(theGrid,c,ACTIVE_CLASS,Factor_One,np->t) != NUM_OK)
        NP_RETURN(1,result[0]);
    }
    FreeVD(theNP->base.mg,mylevel,mylevel,np->t);
  }
    #ifdef ModelP
  if (a_vector_consistent(theMG,np->baselevel,level,np->t) != NUM_OK)
    NP_RETURN(1,result[0]);
    #endif
  for (mylevel=np->baselevel+1; mylevel<level ; mylevel++) {
    if (AllocVDFromVD(theMG,mylevel,mylevel,c,&np->t))
      NP_RETURN(1,result[0]);
    if ((*np->Transfer->InterpolateCorrection)
          (np->Transfer,mylevel,np->t,c,A,Factor_One,result))
      REP_ERR_RETURN(1);
    if (l_daxpy(theGrid,c,EVERY_CLASS,Factor_One,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (l_dmatmul_minus(theGrid,b,NEWDEF_CLASS,A,np->t,EVERY_CLASS)
        != NUM_OK) NP_RETURN(1,result[0]);
    FreeVD(theNP->base.mg,mylevel,mylevel,np->t);
  }
  return (0);
}

static INT AddmgcPostProcess (NP_ITER *theNP, INT level,
                              VECDATA_DESC *x, VECDATA_DESC *b,
                              MATDATA_DESC *A, INT *result)
{
  NP_LMGC *np;
  INT i;

  np = (NP_LMGC *) theNP;

  if (np->Transfer->PostProcess != NULL)
    if ((*np->Transfer->PostProcess)
          (np->Transfer,np->baselevel,level,x,b,A,result))
      REP_ERR_RETURN(1);

  if (np->PreSmooth->PostProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PostProcess)
            (np->PreSmooth,i,x,b,A,result))
        REP_ERR_RETURN(1);

  return (0);
}

static INT AddmgcConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = AddmgcInit;
  theNP->Display = AddmgcDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = AddmgcPreProcess;
  np->Iter = Addmgc;
  np->PostProcess = AddmgcPostProcess;

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
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".gs",sizeof(NP_SMOOTHER),GSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".bcgss",sizeof(NP_BCGSSMOOTHER),BCGSSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sgs",sizeof(NP_SGS),SGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sor",sizeof(NP_SMOOTHER),SORConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sbgs",sizeof(NP_SBGS),SBGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".gbgs",sizeof(NP_SBGS),GBGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ilu",sizeof(NP_ILU),ILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".filu",sizeof(NP_ILU),FILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".thilu",sizeof(NP_THILU),THILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".spilu",sizeof(NP_ILU),SPILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ic",sizeof(NP_ILU),ICConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ff",sizeof(NP_FF),FFConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".lu",sizeof(NP_SMOOTHER),LUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".lmgc",sizeof(NP_LMGC),LmgcConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".addmgc",sizeof(NP_LMGC),AddmgcConstruct))
    REP_ERR_RETURN (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  return (0);
}
