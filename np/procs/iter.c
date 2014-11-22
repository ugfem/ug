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

#include <config.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ugtypes.h"
#include "architecture.h"
#include "ugtime.h"
#include "general.h"
#include "debug.h"
#include "dlmgr.h"
#include "gm.h"
#include "ugm.h"
#include "algebra.h"
#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "ugdevices.h"
#include "udm.h"
#include "pcr.h"
#include "debug.h"
#include "fifo.h"
#include "evm.h"
#include "misc.h"
#include "ugstruct.h"

#include "transfer.h"
#include "ls.h"
#include "iter.h"
#include "project.h"
#include "disctools.h"
#include "block.h"

#include "ff_gen.h"
#include "ff.h"
#include "ugblas.h"
#include "blasm.h"
#include "order.h"

#ifdef USE_FAMG
#include "ug-famg.h"
#endif

#include "ppif_namespace.h"

USING_UG_NAMESPACES
  USING_PPIF_NAMESPACE

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#undef _DEBUG_ITER_

#define TYPE_TFF 1
#define TYPE_FF 2

enum LU_REGULARIZE {
  REG_IF_SING,
  REG_ALWAYS,
  REG_NEVER,
  N_REG
};

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
#define NPFF_ParSim(p)                  ((p)->par_sim)
#define NPFF_AssDirichlet(p)    ((p)->ass_dirichlet)
#define NPFF_SymmFrq(p)                 ((p)->symm_frq)
#define NPFF_CheckSymm(p)               ((p)->check_symm)

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
  NP_ORDER *Order;
  INT Liscomplete;
  INT AutoDamp;
  VECDATA_DESC *DampVector;

    #ifdef ModelP
  INT cons_mode;
  VECDATA_DESC *diag;
    #endif

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

  VEC_SCALAR omega;
  VECDATA_DESC *t;

} NP_SSOR;

typedef struct
{
  NP_SMOOTHER smoother;

  VECDATA_DESC *t;
  INT mode;
  INT depth;
  DOUBLE vdamp;

} NP_PGS;

typedef struct
{
  NP_SMOOTHER smoother;

  VEC_TEMPLATE *vt;
  INT u_sub;
  INT p_sub;
  MAT_TEMPLATE *mt;
  INT uu_sub;
  INT pu_sub;
  INT up_sub;
  INT pp_sub;
  VECDATA_DESC *t;
  VECDATA_DESC *ux;
  VECDATA_DESC *px;
  VECDATA_DESC *ub;
  VECDATA_DESC *pb;
  MATDATA_DESC *uuA;
  MATDATA_DESC *upA;
  MATDATA_DESC *puA;
  MATDATA_DESC *ppA;

  DOUBLE s1;
  DOUBLE s2;
  DOUBLE s3;
  DOUBLE s4;

} NP_BLOCK;

typedef struct
{
  NP_ITER iter;

  VEC_SCALAR damp;
  VECDATA_DESC *u;
  VECDATA_DESC *t;
  VECDATA_DESC *s;
  VECDATA_DESC *p;
  VECDATA_DESC *q;
  VECDATA_DESC *r;
  MATDATA_DESC *L;
  MATDATA_DESC *S;
  VECDATA_DESC *ux;
  VECDATA_DESC *px;
  VECDATA_DESC *ub;
  VECDATA_DESC *pb;
  MATDATA_DESC *uuA;
  MATDATA_DESC *upA;
  MATDATA_DESC *puA;
  MATDATA_DESC *ppA;

  VEC_TEMPLATE *vt;
  INT u_sub;
  INT p_sub;
  MAT_TEMPLATE *mt;
  INT uu_sub;
  INT pu_sub;
  INT up_sub;
  INT pp_sub;

  INT dc;
  INT dc_max;
  INT extra;
  INT display;
  INT ls;
  INT diag;

  DOUBLE thresh;

  NP_ITER *u_iter;
  NP_ITER *v_iter;
  NP_ITER *p_iter;
  NP_LINEAR_SOLVER *u_solver;
  NP_LINEAR_SOLVER *p_solver;

  VEC_SCALAR red;

    #ifdef ModelP
  INT cons_mode;
    #endif

} NP_TS;

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

  VEC_SCALAR beta;
  VEC_SCALAR mindiag;

} NP_SPBLILU;

typedef struct
{
  NP_SMOOTHER smoother;

  INT regularize;

} NP_LU;

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
  INT par_sim;                          /* temp for: simulating parallel algo on SEQ */
  INT ass_dirichlet;                    /* assemble Dirichlet boundary condition (only necessary if not done otherwise (fetransfer) */
  INT symm_frq;                         /* TRUE, if series of testfrequencies should be symmetric, i.e. 1,2,4,8,...,n-1,n,n-1,...,8,4,2,1 */
  INT check_symm;                       /* check, whether the preconditioner is symmetric.
                                           This is done by checking <M^-1*M^-1*d,d> == <M^-1*d,M^-1*d>*/
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

  VEC_SCALAR damp;

} NP_LMGC;

typedef struct
{
  NP_ITER iter;

  INT gamma;
  INT nu1;
  INT nu2;
  INT basenu;
  INT baselevel;

  NP_TRANSFER *Transfer;
  NP_ITER *PreSmooth;
  NP_ITER *PostSmooth;

  VECDATA_DESC *t;
  VECDATA_DESC *d;

  VEC_SCALAR damp;

} NP_II;

typedef struct
{
  NP_ITER iter;

  NP_ITER *Iter;
  NP_TRANSFER *Transfer;

  VECDATA_DESC *s;
  VECDATA_DESC *t;
  VECDATA_DESC *u;

  INT display;
  INT n;

  DOUBLE damp[2*MAXLEVEL];

} NP_CALIBRATE;

typedef struct
{
  NP_SMOOTHER smoother;

  INT nv;                                                       /* # vectors							*/
  INT bw;                                                       /* bandwidth							*/
  INT fmode;                                                    /* apply float-matrix					*/
  INT optimizeBand;                                     /* 1 for optimization of bandwidth		*/
  INT CopyBack;                                         /* 1 for copy decomposed mat back		*/
  INT MarkKey[MAXLEVEL];                        /* key for Mark/Release					*/
  INT count;                            /* counter for MarkKey                  */
  FLOAT *FMat[MAXLEVEL];                        /* float-matrix							*/
  DOUBLE *DMat[MAXLEVEL];                       /* double-matrix						*/
  INT mem;                                                      /* memory used temporary (bytes)		*/
  INT pp_failed;                                        /* 1 if preproc failed, used in smooth  */

  DOUBLE *Vec;                                  /* vector								*/
  INT printdecomptime;                          /* 1 for printing decopmposition time	*/

} NP_EX;

typedef struct
{
  NP_SMOOTHER smoother;

  NP_PROJECT *project;

  VECDATA_DESC *p;
  VECDATA_DESC *t;

  INT mem;                                                          /* memory used temporary (bytes)    */
  INT nv;                                                           /* nb. vectors					    */
  INT bw;                                                           /* bandwidth					    */
  INT fmode;                                                    /* apply float-matrix			        */
  INT optimizeBand;                                         /* 1 for optimization of bandwidth	*/
  FLOAT *FMat;                                              /* float-matrix						*/
  DOUBLE *DMat;                                             /* double-matrix					*/
  DOUBLE *Vec;                                              /* vector						    */

} NP_EXPRJ;

#define MI_MAXITER                                              10

typedef struct
{
  NP_ITER iter;

  VECDATA_DESC *t;                                              /* temporary field					*/
  INT n;                                                                /* number of sub iterations			*/
  NP_ITER *si[MI_MAXITER];                              /* sub iterations					*/

} NP_MI;

#define SP_NOTINIT                                              0
#define SP_OMS                                                  1
#define SP_A                                                    2
#define SP_O                                                    3

typedef struct
{
  NP_ITER iter;

  /* necessary componentes for 'activ' mode */
  VECDATA_DESC *t;                      /* temporary field                  */
  INT nu;                               /* number of sub iterations         */
  NP_ITER *S;                                       /* sub iteration                    */
  INT mode;                                                             /* SP_OMS: (1-S)*S^nu				*/
  /* SP_A:       A*S^nu				*/

  /* necessary componentes for 'executable' mode */
  INT maxiter;                                                  /* maximum number of iterations     */
  DOUBLE tol;                                                           /* tolerance, stopping criterion    */

} NP_SP;

typedef struct
{
  NP_BASE base;

  MATDATA_DESC *A;                      /* matrix to handle                 */
  VECDATA_DESC *x;                      /* in-field                         */
  INT xc;                                                               /* scalar component in x			*/
  NP_ITER *S;                           /* sub iteration                    */
  INT skip;                             /* 1 for skipping dirichlet flags   */
  char file[128];                                               /* file name						*/

} NP_IM;


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
static char LU_reg[N_REG][16];

REP_ERR_FILE;

static VECDATA_DESC *FF_VECDATA_DESC_ARRAY[FF_MAX_VECS];


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static INT SetAutoDamp (NP_ITER *theNP, GRID *g, INT mode, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv);

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

INT DPrintVector (MULTIGRID *mg, VECDATA_DESC *x)
{
  VECTOR *v;
  DOUBLE_VECTOR pos;
  INT vtype,comp,ncomp;
  FILE *f;

  f=fopen("logfiles/x","w");
  for (vtype=0; vtype<NVECTYPES; vtype++)
  {
    ncomp = VD_NCMPS_IN_TYPE(x,vtype);
    if (ncomp == 0) continue;
    comp = VD_CMP_OF_TYPE(x,vtype,0);
    S_FINE_VLOOP__TYPE(CURRENTLEVEL(mg),v,mg,vtype)
    {
      VectorPosition(v,pos);
      fprintf(f,"%15.8e\n",VVALUE(v,comp));
    }
  }
  fclose(f);

  return(NUM_OK);
}

INT NS_DIM_PREFIX NPIterInit (NP_ITER *np, INT argc , char **argv)
{
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->c = ReadArgvVecDesc(np->base.mg,"c",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"r",argc,argv);

  if ((np->A == NULL) || (np->b == NULL) || (np->c == NULL))
    return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NS_DIM_PREFIX NPIterDisplay (NP_ITER *np)
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

INT NS_DIM_PREFIX NPIterExecute (NP_BASE *theNP, INT argc , char **argv)
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
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  np->L = ReadArgvMatDesc(theNP->mg,"L",argc,argv);
  np->Order = (NP_ORDER*)ReadArgvNumProc(theNP->mg,"O",ORDER_CLASS_NAME,argc,argv);
    #ifdef ModelP
  if (ReadArgvINT("Liscomplete",&(np->Liscomplete),argc,argv))
    np->Liscomplete = 0;
  if( np->Liscomplete && (np->L==NULL) )
  {
    PrintErrorMessage('E',"SmootherInit","L is mandatory for option Liscomplete");
    REP_ERR_RETURN (1);
  }
  np->diag = ReadArgvVecDesc(theNP->mg,"diag",argc,argv);
  if (ReadArgvOption("M",argc,argv))
    np->cons_mode = MAT_CONS;
  else if (ReadArgvOption("D",argc,argv))
    np->cons_mode = MAT_DIAG_CONS;
  else if (np->diag != NULL)
    np->cons_mode = MAT_DIAG_VEC_CONS;
  else
    np->cons_mode = MAT_MASTER_CONS;
        #endif
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
  if (np->Order != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Order",ENVITEM_NAME(np->Order));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Order","---");
    #ifdef ModelP
  UserWriteF(DISPLAY_NP_FORMAT_SI,"cons_mode",(int)np->cons_mode);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"Liscomplete",(int)np->Liscomplete);
        #endif

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
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->damp) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK)
    NP_RETURN(1,result[0]);

  return (0);
}

static INT SmootherPostProcess (NP_ITER *theNP, INT level,
                                VECDATA_DESC *x, VECDATA_DESC *b,
                                MATDATA_DESC *A, INT *result)
{
  NP_SMOOTHER *np;

  np = (NP_SMOOTHER *) theNP;
  if (np->L != NULL)
    if (FreeMD(NP_MG(theNP),level,level,np->L))
      REP_ERR_RETURN(1);
  np->Liscomplete = 0;
        #ifdef ModelP
  if (np->cons_mode == MAT_DIAG_VEC_CONS)
    if (FreeVD(NP_MG(theNP),level,level,np->diag))
      REP_ERR_RETURN(1);
        #endif

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
  NP_SMOOTHER *np=(NP_SMOOTHER *) theNP;

        #ifdef ModelP
  GRID *theGrid;
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->L)) NP_RETURN(1,result[0]);
  theGrid = NP_GRID(theNP,level);
  if( !np->Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,A) != NUM_OK) NP_RETURN(1,result[0]);
    if (l_matrix_consistent(theGrid,np->L,MAT_DIAG_CONS) != NUM_OK) NP_RETURN(1,result[0]);
  }
        #else
  np->L = A;
        #endif
  *baselevel = level;

  return (0);
}

static INT JacobiStep (NP_SMOOTHER *theNP, INT level,
                       VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                       MATDATA_DESC *L,
                       INT *result)
{
  if (l_jac(NP_GRID(theNP,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

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

static INT CopyDiagMatrix (GRID *theGrid, MATDATA_DESC *A, VECDATA_DESC *diag)
{
  VECTOR *v;

  if (MD_IS_SCALAR(A)) {
    register SHORT mc = MD_SCALCMP(A);
    register SHORT vc = VD_SCALCMP(diag);
    register SHORT mask = MD_SCAL_RTYPEMASK(A);

    for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v))
      if ((VDATATYPE(v)&mask))
        VVALUE(v,vc) = MVALUE(START(v),mc);

    return (0);
  }
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v))
  {
    INT tp = VTYPE(v);
    INT i;
    INT n = VD_NCMPS_IN_TYPE(diag,tp);

    for (i=0; i<n; i++)
      VVALUE(v,VD_CMP_OF_TYPE(diag,tp,i)) =
        MVALUE(START(v),MD_MCMP_OF_RT_CT(A,tp,tp,i));
  }

  return(0);
}

static INT GSPreProcess  (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b,
                          MATDATA_DESC *A, INT *baselevel, INT *result)
{
  GRID *theGrid=NP_GRID(theNP,level);

  NP_SMOOTHER *np=(NP_SMOOTHER *) theNP;

        #ifdef ModelP
  if (np->cons_mode == MAT_DIAG_VEC_CONS) {
    CopyDiagMatrix(theGrid,A,np->diag);
    if (l_vector_consistent(theGrid,np->diag) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->L))
      NP_RETURN(1,result[0]);
    if( !np->Liscomplete )
    {
      if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,A) != NUM_OK)
        NP_RETURN(1,result[0]);
      if (l_matrix_consistent(theGrid,np->L,np->cons_mode) != NUM_OK)
        NP_RETURN(1,result[0]);
    }
  }
        #endif
  if (np->Order!=NULL)
    if ((*np->Order->Order)(np->Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid))
    NP_RETURN(1,result[0]);
  *baselevel = level;

  return (0);
}


static INT GSStep (NP_SMOOTHER *theNP, INT level,
                   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                   MATDATA_DESC *L,
                   INT *result)
{
  NP_SMOOTHER *np;
  GRID *theGrid;

  np = (NP_SMOOTHER *) theNP;
  theGrid = NP_GRID(theNP,level);
    #ifdef ModelP
  if (np->cons_mode == MAT_MASTER_CONS) {
    if (l_vector_collect(theGrid,b)!=NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  if (np->cons_mode == MAT_DIAG_VEC_CONS) {
    if (l_lgs(NP_GRID(theNP,level),x,A,b,np->diag) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else
  if (l_lgs(NP_GRID(theNP,level),x,L,b,NULL) != NUM_OK)
    NP_RETURN(1,result[0]);
    #else
  if (l_lgs(NP_GRID(theNP,level),x,A,b,NULL) != NUM_OK) NP_RETURN(1,result[0]);
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
  INT i,restart;
  DOUBLE alpha,rho_new,beta,tt;

  np = (NP_BCGSSMOOTHER *) theNP;
  restart = 1;
  for (i=0; i<np->maxiter; i++)
  {
    /* restart ? */
    if ((np->restart>0 && i%np->restart==0) || restart)
    {
      if (dset(NP_MG(theNP),level,level,ALL_VECTORS,np->p,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (dset(NP_MG(theNP),level,level,ALL_VECTORS,np->v,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->r,b)!= NUM_OK) REP_ERR_RETURN (1);
      alpha = np->rho = np->omega = 1.0;
      restart = 0;
    }

    /* update x, b */
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,b,np->r,&rho_new)!=NUM_OK) REP_ERR_RETURN (1);
    beta = rho_new*alpha/np->rho/np->omega;
    if (dscal(NP_MG(theNP),level,level,ALL_VECTORS,np->p,beta)) REP_ERR_RETURN (1);
    if (dadd(NP_MG(theNP),level,level,ALL_VECTORS,np->p,b)!= NUM_OK) REP_ERR_RETURN (1);
    if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,np->p,-beta*np->omega,np->v)!= NUM_OK) REP_ERR_RETURN (1);
    if (np->Iter!=NULL)
    {
      if (dset(NP_MG(theNP),level,level,ALL_VECTORS,np->q,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->s,np->p)!= NUM_OK) REP_ERR_RETURN (1);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->p,A,result)) REP_ERR_RETURN (1);
      if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->p,np->s)!= NUM_OK) REP_ERR_RETURN (1);
      if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,np->v,A,np->q)) REP_ERR_RETURN (1);
      if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,np->v,np->r,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      alpha = rho_new/alpha;
      if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,x,alpha,np->q)!= NUM_OK) REP_ERR_RETURN (1);
    }
    else
    {
      if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,np->v,A,np->p)) REP_ERR_RETURN (1);
      if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,np->v,np->r,&alpha)!=NUM_OK) REP_ERR_RETURN (1);
      alpha = rho_new/alpha;
      if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,x,alpha,np->p)!= NUM_OK) REP_ERR_RETURN (1);
    }
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->s,b)!= NUM_OK) REP_ERR_RETURN (1);
    if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,np->s,-alpha,np->v)!= NUM_OK) REP_ERR_RETURN (1);
    if (np->Iter!=NULL)
    {
      if (dset(NP_MG(theNP),level,level,ALL_VECTORS,np->q,0.0)!= NUM_OK) REP_ERR_RETURN (1);
      if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->t,np->s)!= NUM_OK) REP_ERR_RETURN (1);
      if ((*np->Iter->Iter)(np->Iter,level,np->q,np->s,A,result)) REP_ERR_RETURN (1);
      if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->s,np->t)!= NUM_OK) REP_ERR_RETURN (1);
    }
    else
    {
      if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->q,np->s)!= NUM_OK) REP_ERR_RETURN (1);
    }
    if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,np->t,A,np->q)) REP_ERR_RETURN (1);
    if (dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,np->t,&tt)!=NUM_OK) REP_ERR_RETURN (1);
    tt *= tt;
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,np->s,np->t,&(np->omega))!=NUM_OK) REP_ERR_RETURN (1);
    np->omega /= tt;
    if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,x,np->omega,np->q)!= NUM_OK) REP_ERR_RETURN (1);
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,b,np->s)!= NUM_OK) REP_ERR_RETURN (1);
    if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,b,-np->omega,np->t)!= NUM_OK) REP_ERR_RETURN (1);
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
  if (nps->L != NULL) if (FreeMD(NP_MG(theNP),level,level,nps->L)) REP_ERR_RETURN(1);
  nps->Liscomplete = 0;

  np = (NP_BCGSSMOOTHER *) theNP;
  if (FreeVD(nps->iter.base.mg,level,level,np->r)) REP_ERR_RETURN(1);
  if (FreeVD(nps->iter.base.mg,level,level,np->p)) REP_ERR_RETURN(1);
  if (FreeVD(nps->iter.base.mg,level,level,np->v)) REP_ERR_RETURN(1);
  if (FreeVD(nps->iter.base.mg,level,level,np->s)) REP_ERR_RETURN(1);
  if (FreeVD(nps->iter.base.mg,level,level,np->t)) REP_ERR_RETURN(1);
  if (FreeVD(nps->iter.base.mg,level,level,np->q)) REP_ERR_RETURN(1);

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
  GRID *theGrid=NP_GRID(theNP,level);
  NP_SMOOTHER *np=(NP_SMOOTHER *) theNP;

        #ifdef ModelP
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->L))
    NP_RETURN(1,result[0]);
  if( !np->Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,A) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (l_matrix_consistent(theGrid,np->L,np->cons_mode) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
        #endif
  if (np->Order!=NULL)
    if ((*np->Order->Order)(np->Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid))
    NP_RETURN(1,result[0]);
  *baselevel = level;

  /* get storage for extra temp */
  if (AllocVDFromVD(NP_MG(theNP),level,level,x,&NP_SGS_t((NP_SGS *)theNP)))
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
  if (np->smoother.cons_mode == MAT_MASTER_CONS) {
    if (l_vector_collect(theGrid,b)!=NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  if (l_lgs(theGrid,NP_SGS_t(np),np->smoother.L,b,NULL)
      != NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_vector_consistent(theGrid,NP_SGS_t(np)) != NUM_OK)
    NP_RETURN(1,result[0]);
    #else
  if (l_lgs(theGrid,NP_SGS_t(np),A,b,NULL)) NP_RETURN(1,result[0]);
    #endif

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,NP_SGS_t(np),np->smoother.damp)!= NUM_OK)
    NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,NP_SGS_t(np))
      != NUM_OK) NP_RETURN(1,result[0]);

  /* iterate backward */
    #ifdef ModelP
  if (np->smoother.cons_mode == MAT_MASTER_CONS) {
    if (l_vector_collect(theGrid,b)!=NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  if (l_ugs(theGrid,x,np->smoother.L,b))
    NP_RETURN(1,result[0]);
  if (l_vector_consistent(theGrid,x) != NUM_OK) NP_RETURN(1,result[0]);
        #else
  if (l_ugs(theGrid,x,A,b)) NP_RETURN(1,result[0]);
    #endif

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->smoother.damp) != NUM_OK) NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK) NP_RETURN(1,result[0]);

  /* now add the two corrections */
  if (dadd(NP_MG(theNP),level,level,ALL_VECTORS,x,NP_SGS_t(np)) != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SGSPostProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *result)
{
  NP_SGS *np;

  np = (NP_SGS *) theNP;
  if (FreeVD(NP_MG(theNP),level,level,NP_SGS_t(np)) ) REP_ERR_RETURN(1);

  return(SmootherPostProcess(theNP,level,x,b,A,result));
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
   ssor - numproc for symmetric SOR smoother

   DESCRIPTION:
   This numproc executes a symmetric SOR smoother,
   using the blas routines
   'l_lsor' and 'l_usor'. It can be used in 'lmgc'.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list> $omega <sc double list>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  $omega~<sc~double~list> - overrelaxation factors for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT SSORInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SSOR *np = (NP_SSOR *) theNP;
  INT i;

  if (sc_read (np->omega,NP_FMT(np),NULL,"omega",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->omega[i] = 1.0;
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->smoother.AutoDamp = ReadArgvOption("autodmp",argc,argv);
  np->smoother.DampVector = ReadArgvVecDesc(NP_MG(np),"dv",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT SSORDisplay (NP_BASE *theNP)
{
  NP_SSOR *np = (NP_SSOR *) theNP;

  SmootherDisplay(theNP);

  if (np->smoother.iter.b != NULL)
    if (sc_disp(np->omega,np->smoother.iter.b,"omega"))
      REP_ERR_RETURN (1);
  if (np->t != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"t","---");
  if (np->smoother.DampVector != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"dv",ENVITEM_NAME(np->smoother.DampVector));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"dv","---");


  return (0);
}

static INT cmp_comp;
static INT CompareADV (const void *p, const void *q)
{
  VECTOR *v1, *v2;

  v1 = *((VECTOR **)p);
  v2 = *((VECTOR **)q);
  if (VVALUE(v1,cmp_comp)<VVALUE(v2,cmp_comp)) {return (1);}
  if (VVALUE(v1,cmp_comp)>VVALUE(v2,cmp_comp)) {return (-1);}
  return (0);
}

static INT SSORPreProcess  (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b,
                            MATDATA_DESC *A, INT *baselevel, INT *result)
{
  GRID *theGrid=NP_GRID(theNP,level);
  NP_SMOOTHER *np=(NP_SMOOTHER *) theNP;
  NP_SSOR *ssor = (NP_SSOR *) theNP;
  MATDATA_DESC *myMat;

        #ifdef ModelP
  if (np->cons_mode == MAT_DIAG_VEC_CONS) {
    CopyDiagMatrix(theGrid,A,np->diag);
    if (l_vector_consistent(theGrid,np->diag) != NUM_OK)
      NP_RETURN(1,result[0]);
    np->L = A;
  }
  else {
    if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->L))
      NP_RETURN(1,result[0]);
    if( !np->Liscomplete )
    {
      if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,A) != NUM_OK)
        NP_RETURN(1,result[0]);
      if (l_matrix_consistent(theGrid,np->L,np->cons_mode) != NUM_OK)
        NP_RETURN(1,result[0]);
    }
    np->diag = NULL;
  }
  myMat = np->L;
        #else
  myMat = A;
        #endif

  /* get storage for extra temp */
  if (AllocVDFromVD(NP_MG(theNP),level,level,x,&ssor->t))
    NP_RETURN(1,result[0]);

  /* set auto-damping iff */
  if (np->AutoDamp)
  {
    if (AllocVDFromVD(NP_MG(theNP),level,level,x,&np->DampVector)) NP_RETURN(1,result[0]);
    if (SetAutoDamp(theNP,theGrid,np->AutoDamp,myMat,ssor->omega,np->DampVector)) NP_RETURN(1,result[0]);
    /* DPrintVector(MYMG(theGrid),np->DampVector); */
  }
  if (np->Order!=NULL)
    if ((*np->Order->Order)(np->Order,level,myMat,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid))
    NP_RETURN(1,result[0]);
  *baselevel = level;

  return (0);
}

static INT SSORSmoother (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SSOR *np = (NP_SSOR *) theNP;
  GRID *theGrid=NP_GRID(theNP,level);
  VECDATA_DESC *b0;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  /* iterate forward */
    #ifdef ModelP
  if (np->smoother.cons_mode == MAT_MASTER_CONS)
  {
    if (l_vector_collect(theGrid,b)!=NUM_OK) NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
  }
  if (np->smoother.AutoDamp)
  {
    if (l_lsor_ld(theGrid,np->t,np->smoother.L,b, np->smoother.DampVector, np->smoother.diag) != NUM_OK) NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_lsor(theGrid,np->t,np->smoother.L,b, np->omega,np->smoother.diag) != NUM_OK) NP_RETURN(1,result[0]);
  }
  if (l_vector_consistent(theGrid,np->t) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (np->smoother.AutoDamp)
  {
    if (l_lsor_ld(theGrid,np->t,A,b,np->smoother.DampVector,NULL) != NUM_OK) NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_lsor(theGrid,np->t,A,b,np->omega,NULL) != NUM_OK) NP_RETURN(1,result[0]);
  }
    #endif

  /* store unmodified b */
  b0=NULL; if (AllocVDFromVD(NP_MG(theNP),level,level,b,&b0)) NP_RETURN(1,result[0]);
  if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,b0,b) != NUM_OK) NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,np->t) != NUM_OK) NP_RETURN(1,result[0]);

  /* iterate backward */
    #ifdef ModelP
  if (np->smoother.cons_mode == MAT_MASTER_CONS)
  {
    if (l_vector_collect(theGrid,b)!=NUM_OK) NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
  }
  if (np->smoother.AutoDamp)
  {
    if (l_usor_ld(theGrid,x,np->smoother.L,b, np->smoother.DampVector, np->smoother.diag) != NUM_OK) NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_usor(theGrid,x,np->smoother.L,b,np->omega,np->smoother.diag) != NUM_OK) NP_RETURN(1,result[0]);
  }
  if (l_vector_consistent(theGrid,x) != NUM_OK) NP_RETURN(1,result[0]);
    #else
  if (np->smoother.AutoDamp)
  {
    if (l_usor_ld(theGrid,x,A,b,np->smoother.DampVector,NULL) != NUM_OK) NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_usor(theGrid,x,A,b,np->omega,NULL) != NUM_OK) NP_RETURN(1,result[0]);
  }
    #endif

  /* copy back unmodified b */
  if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,b,b0) != NUM_OK) NP_RETURN(1,result[0]);
  if (FreeVD(NP_MG(theNP),level,level,b0)) NP_RETURN(1,result[0]);

  /* add the two corrections */
  if (dadd(NP_MG(theNP),level,level,ALL_VECTORS,x,np->t) != NUM_OK) NP_RETURN(1,result[0]);

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->smoother.damp) != NUM_OK) NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SSORPostProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b,
                            MATDATA_DESC *A, INT *result)
{
  NP_SSOR *np = (NP_SSOR *) theNP;

  if (FreeVD(NP_MG(theNP),level,level,np->t))
    REP_ERR_RETURN(1);

  if (np->smoother.AutoDamp)
    if (FreeVD(NP_MG(theNP),level,level,np->smoother.DampVector))
      REP_ERR_RETURN(1);

  return(SmootherPostProcess(theNP,level,x,b,A,result));
}

static INT SSORConstruct (NP_BASE *theNP)
{
  NP_ITER *np = (NP_ITER *) theNP;

  theNP->Init = SSORInit;
  theNP->Display = SSORDisplay;
  theNP->Execute = NPIterExecute;
  np->PreProcess = SSORPreProcess;
  np->Iter = SSORSmoother;
  np->PostProcess = SSORPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   pgs - numproc for patch Gauss-Seidel smoother

   DESCRIPTION:
   This numproc executes a patch Gauss-Seidel smoother,
   using the blas routine 'l_pgs'.
   It depends strongly on the order in the matrix list.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list> $mode <mode> $depth <depth>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - damping factors for each component
   .  $depth~<depth> - number of patch vectors
   .  $mode~<mode> - overlapping mode

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT PGSInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_PGS *np;

  np = (NP_PGS *) theNP;
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  if (ReadArgvINT("mode",&(np->mode),argc,argv))
    np->mode = 0;
  if (ReadArgvINT("depth",&(np->depth),argc,argv))
    np->depth = 2;

  if (ReadArgvDOUBLE("vdamp",&(np->vdamp),argc,argv))
    np->vdamp = 1.0;

  return (SmootherInit(theNP,argc,argv));
}

static INT PGSDisplay (NP_BASE *theNP)
{
  NP_PGS *np;

  SmootherDisplay(theNP);
  np = (NP_PGS *) theNP;
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"mode",(int)np->mode);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"depth",(int)np->depth);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"vdamp",np->vdamp);

  return (0);
}

struct matrix_list {
  MATRIX *m;
  DOUBLE value;
};

static int sort_MatArray (const void *e1, const void *e2)
{
  DOUBLE v1 = (*((struct matrix_list *)e1)).value;
  DOUBLE v2 = (*((struct matrix_list *)e2)).value;

  if (v1 < v2) return(-1);
  if (v1 > v2) return(1);
  return (0);
}

static INT SortMatrices (GRID *theGrid, MATDATA_DESC *A)
{
  VECTOR *v;
  MATRIX *m;
  struct matrix_list *table;
  INT cnt,max,rtype,ctype,rcomp,ccomp,i,j;
  DOUBLE val;
  SHORT *Mcomp;
  INT MarkKey;

  max = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(A,rtype,rtype);
    if (rcomp == 0) continue;
    cnt = 0;
    m = VSTART(v);
    ASSERT(m != NULL);
    ASSERT(MDEST(m) == v);
    for (m=MNEXT(m); m!=NULL; m=MNEXT(m)) cnt ++;
    max = MAX(max,cnt);
  }
  MarkTmpMem(MGHEAP(MYMG(theGrid)),&MarkKey);
  table = (struct matrix_list *)GetTmpMem(MGHEAP(MYMG(theGrid)),
                                          max*sizeof(struct matrix_list),MarkKey);
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(A,rtype,rtype);
    if (rcomp == 0) continue;
    cnt = 0;
    m = VSTART(v);
    for (m=MNEXT(m); m!=NULL; m=MNEXT(m)) {
      ctype = MDESTTYPE(m);
      ccomp = MD_COLS_IN_RT_CT(A,rtype,ctype);
      table[cnt].m = m;
      if (ccomp > 0)
        Mcomp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
      val = 0;
      for (i=0; i<rcomp*ccomp; i++)
        val += ABS(MVALUE(m,Mcomp[i]));
      table[cnt].value = val;
      cnt++;
    }
    if (cnt < 2) continue;
    qsort(table,cnt,sizeof(struct matrix_list),sort_MatArray);
    m=VSTART(v);
    for(j=0; j<cnt; j++) {
      MNEXT(m) = table[j].m;
      m = MNEXT(m);
    }
    MNEXT(m)=NULL;
  }
  ReleaseTmpMem(MGHEAP(MYMG(theGrid)),MarkKey);

  return(0);
}

static INT PGSPreProcess  (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_PGS *np;
  GRID *theGrid;

  np = (NP_PGS *) theNP;
  theGrid = NP_GRID(theNP,level);
        #ifdef ModelP
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L))
    NP_RETURN(1,result[0]);
  if( !np->smoother.Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (l_matrix_consistent(theGrid,np->smoother.L,MAT_CONS) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
        #endif
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid))
    NP_RETURN(1,result[0]);
  *baselevel = level;

  /* get storage for extra temp */
  if (AllocVDFromVD(NP_MG(theNP),level,level,x,&np->t))
    NP_RETURN(1,result[0]);
  /*
     if (SortMatrices(theGrid,A))
          NP_RETURN(1,result[0]);
   */
  return (0);
}

static INT PGSSmoother (NP_ITER *theNP, INT level,
                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                        INT *result)
{
  NP_PGS *np;
  GRID *theGrid;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  np = (NP_PGS *) theNP;
  theGrid = NP_GRID(theNP,level);

  /* iterate forward */
    #ifdef ModelP
  if (l_vector_collect(theGrid,b)!=NUM_OK) NP_RETURN(1,result[0]);
  if (dset(NP_MG(theNP),level,level,ALL_VECTORS,x,0.0)!= NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_pgs(theGrid,x,np->smoother.L,b,np->depth,np->mode,np->vdamp)
      != NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_vector_consistent(theGrid,x) != NUM_OK)
    NP_RETURN(1,result[0]);
    #else
  if (l_pgs(theGrid,x,A,b,np->depth,np->mode,np->vdamp))
    NP_RETURN(1,result[0]);
    #endif

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->smoother.damp)!= NUM_OK)
    NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK)
    NP_RETURN(1,result[0]);

  return (0);
}

static INT PGSPostProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *result)
{
  NP_PGS *np;

  np = (NP_PGS *) theNP;

  if (FreeVD(NP_MG(theNP),level,level,np->t)) REP_ERR_RETURN(1);

  return(SmootherPostProcess(theNP,level,x,b,A,result));
}

static INT PGSConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = PGSInit;
  theNP->Display = PGSDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = PGSPreProcess;
  np->Iter = PGSSmoother;
  np->PostProcess = PGSPostProcess;

  return(0);
}

#define MAX_PATCH 100
#define MAX_VPATCH 100

static INT SolveFullMatrix3 (INT n, DOUBLE *sol, DOUBLE mat[MAX_PATCH][MAX_PATCH],
                             DOUBLE *rhs)
{
  register DOUBLE dinv,piv,sum;
  register int i,j,k;
  INT ipv[MAX_PATCH];

  /*
     UserWrite("m\n");
     for (i=0; i<n; i++) {
     for (j=0; j<n; j++)
          UserWriteF("%8.1e",mat[i][j]);
     UserWrite("\n");
     }
   */

  if (n > MAX_PATCH)
    return (1);

  for (i=0; i<n; i++)
    ipv[i] = i;

  /* lr factorize mat */
  for (i=0; i<n; i++)
  {
    k = i;
    piv = ABS(mat[i][i]);
    for (j=i+1; j<n; j++)
    {
      sum = ABS(mat[j][i]);
      if (sum > piv)
      {
        k = j;
        piv = sum;
      }
    }
    if (k != i)
    {
      j = ipv[i];
      ipv[i] = ipv[k];
      ipv[k] = j;
      for (j=0; j<n; j++)
      {
        sum = mat[k][j];
        mat[k][j] = mat[i][j];
        mat[i][j] = sum;
      }
    }

    dinv = mat[i][i];
    if (ABS(dinv)<1e-30)
      return(NUM_SMALL_DIAG);
    dinv = mat[i][i] = 1.0/dinv;
    for (j=i+1; j<n; j++)
    {
      piv = (mat[j][i] *= dinv);
      for (k=i+1; k<n; k++)
        mat[j][k] -= mat[i][k] * piv;
    }
  }

  /*
     UserWrite("im\n");
     for (i=0; i<n; i++) {
     for (j=0; j<n; j++)
          UserWriteF("%8.1e",mat[i][j]);
     UserWrite("\n");
     }
   */

  /* solve */
  for (i=0; i<n; i++)
  {
    for (sum=rhs[ipv[i]], j=0; j<i; j++)
      sum -= mat[i][j] * sol[j];
    sol[i] = sum;               /* Lii = 1 */
  }
  for (i=n-1; i>=0; i--)
  {
    for (sum=sol[i], j=i+1; j<n; j++)
      sum -= mat[i][j] * sol[j];
    sol[i] = sum * mat[i][i];                   /* Uii = Inv(Mii) */
  }

  return (NUM_OK);
}

static INT InvertFullMatrix3 (INT n, DOUBLE mat[MAX_PATCH][MAX_PATCH],
                              DOUBLE *inv)
{
  register DOUBLE dinv,piv,sum;
  DOUBLE rhs[MAX_PATCH];
  register INT i,j,k;
  INT ipv[MAX_PATCH];

  if (n > MAX_PATCH) {
    PrintErrorMessage('E',"InvertFullMatrix3","n too large");
    return (1);
  }
  for (i=0; i<n; i++)
    ipv[i] = i;

  /* lr factorize mat */
  for (i=0; i<n; i++) {
    k = i;
    piv = ABS(mat[i][i]);
    for (j=i+1; j<n; j++)
    {
      sum = ABS(mat[j][i]);
      if (sum > piv)
      {
        k = j;
        piv = sum;
      }
    }
    if (k != i)
    {
      j = ipv[i];
      ipv[i] = ipv[k];
      ipv[k] = j;
      for (j=0; j<n; j++)
      {
        sum = mat[k][j];
        mat[k][j] = mat[i][j];
        mat[i][j] = sum;
      }
    }

    dinv = mat[i][i];
    if (ABS(dinv)<1e-30) {
      UserWriteF("InvertFullMatrix3: dinv too small\n");

      printf("n %d\n",n);

      /*
         for (i=0; i<n; i++) {
         for (j=1; j<n; j++)
          printf("%8.1e",mat[i][j]);
         printf("\n");
         }
       */


      return (1);
    }
    dinv = mat[i][i] = 1.0/dinv;
    for (j=i+1; j<n; j++)
    {
      piv = (mat[j][i] *= dinv);
      for (k=i+1; k<n; k++)
        mat[j][k] -= mat[i][k] * piv;
    }
  }

  /* solve */
  for (k=0; k<n; k++)
  {
    for (i=0; i<n; i++)
      rhs[i] = 0;
    rhs[k] = 1.0;
    for (i=0; i<n; i++)
    {
      for (sum=rhs[ipv[i]], j=0; j<i; j++)
        sum -= mat[i][j] * inv[j*n+k];
      inv[i*n+k] = sum;                         /* Lii = 1 */
    }
    for (i=n-1; i>=0; i--)
    {
      for (sum=inv[i*n+k], j=i+1; j<n; j++)
        sum -= mat[i][j] * inv[j*n+k];
      inv[i*n+k] = sum * mat[i][i];                     /* Uii = Inv(Mii) */
    }
  }

  return (0);
}

INT l_block_collect (GRID *theGrid,
                     VECDATA_DESC *ux, VECDATA_DESC *px,
                     VECDATA_DESC *ub, VECDATA_DESC *pb,
                     MATDATA_DESC *uuA,
                     MATDATA_DESC *upA,
                     MATDATA_DESC *puA,
                     MATDATA_DESC *ppA,
                     double s1,
                     double s2,
                     double s3,
                     double s4)
{
  VECTOR *v;

  INT nmax = NDATA_DEF_IN_GRID(theGrid) / 8;
  DOUBLE *ndata;


  ASSERT(nmax >= MAX_PATCH*MAX_PATCH);

  /* s1 = 4;    s2 = 1;         s3 = 1;         s4 = 2; */

  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v))
  {
    INT ptype = VTYPE(v);
    INT pcomp = VD_NCMPS_IN_TYPE(pb,ptype);
    DOUBLE matr[MAX_PATCH][MAX_PATCH];
    VECTOR *w[MAX_VPATCH];
    INT ucomp[MAX_VPATCH];
    INT utype[MAX_VPATCH];
    MATRIX *m;
    INT n = pcomp;
    INT cnt = 0;
    INT i,j;

    if (pcomp == 0) continue;
    ASSERT(ptype==NODEVEC);
    ndata = (DOUBLE *)NDATA((NODE*)VOBJECT(v));

    if (ndata == NULL)
    {
      NODE *theNode = (NODE *) VOBJECT(v);
      NDATA(theNode) = (void *)GetMemoryForObject(theGrid->mg,nmax*8,-1);
      ndata = (DOUBLE *)NDATA(theNode);
    }

    for (i=0; i<MAX_PATCH; i++)
      for (j=0; j<MAX_PATCH; j++)
        matr[i][j] = 0.0;

    /*
       printf("A ");
       for (i=0; i<MAX_PATCH; i++)
        printf("%8.1e",matr[0][i]);
       printf("\n");
     */

    m = START(v);
    for (i=0; i<pcomp; i++) {
      for (j=0; j<pcomp; j++)
        matr[i][j] = s3 *
                     MVALUE(m,MD_IJ_CMP_OF_RT_CT(ppA,ptype,ptype,i,j));
    }
    for (; m!=NULL; m=NEXT(m))
    {
      w[cnt] = MDEST(m);
      utype[cnt] = VTYPE(w[cnt]);
      ucomp[cnt] = VD_NCMPS_IN_TYPE(ux,utype[cnt]);
      if (ucomp[cnt] == 0) continue;
      if (n+ucomp[cnt] > MAX_PATCH) break;
      for (i=0; i<ucomp[cnt]; i++) {
        for (j=0; j<ucomp[cnt]; j++)
          matr[n+i][n+j] = s1 *
                           MVALUE(START(w[cnt]),MD_IJ_CMP_OF_RT_CT(uuA,utype[cnt],
                                                                   utype[cnt],i,j));
        for (j=0; j<pcomp; j++) {
          matr[n+i][j] = s2 *
                         MVALUE(MADJ(m),MD_IJ_CMP_OF_RT_CT(upA,utype[cnt],ptype,
                                                           i,j));
          matr[j][n+i] = s2 *
                         MVALUE(m,MD_IJ_CMP_OF_RT_CT(puA,ptype,utype[cnt],
                                                     j,i));
        }
      }
      /*
             printf("%d:matr %d %d type %d mtp %d comp %d val %8.1e\n",
                                me,VINDEX(v),VINDEX(w[cnt]),utype[cnt],
                                MTP(ptype,utype[cnt]),
                                MD_IJ_CMP_OF_RT_CT(puA,ptype,utype[cnt],0,1),
                 MVALUE(m,MD_IJ_CMP_OF_RT_CT(puA,ptype,utype[cnt],0,1)));
                      printf("B ");
                      for (i=0; i<MAX_PATCH; i++)
                          printf("%8.1e",matr[0][i]);
                      printf("\n");
       */

      n += ucomp[cnt];
      cnt++;
      if (cnt == MAX_VPATCH) break;
    }
    n = pcomp;
    for (i=0; i<cnt; i++)
    {
      INT n1 = pcomp;

      for (j=0; j<cnt; j++) {
        if (i != j)
        {
          INT k,l;

          m = GetMatrix(w[i],w[j]);
          if (m != NULL)
            for (k=0; k<ucomp[i]; k++)
              for (l=0; l<ucomp[j]; l++)
                matr[n+k][n1+l] = s4 *
                                  MVALUE(m,MD_IJ_CMP_OF_RT_CT(uuA,utype[i],utype[j],k,l));
        }
        n1 += ucomp[j];
      }
      n += ucomp[i];
    }
    /*
       if (me == 1) {
       printf("b %d\n",n);
       printf("%d:matr %d %8.1e\n",me,VINDEX(v),matr[0][0]);
                    for (i=0; i<n; i++) {
       for (j=0; j<n; j++)
        printf("%8.1e",matr[i][j]);
       printf("\n");
       }
       }
       printf("%d: matr %d =",me,VINDEX(v));
       for (i=0; i<n; i++)
        printf("%8.1e",matr[0][i]);
       printf("\n");
     */

    if (InvertFullMatrix3(n,matr,ndata)) {
      for (i=0; i<n; i++)
        for (j=0; j<n; j++)
          if (i != j)
            matr[i][j] = 0.0;
          else
            matr[i][i] = 1.0;
      if (InvertFullMatrix3(n,matr,ndata)) {
        matr[0][0] = 1.0;
        RETURN(1);
      }
    }
  }
  return(0);
}

INT l_block (GRID *theGrid,
             VECDATA_DESC *ux, VECDATA_DESC *px,
             VECDATA_DESC *ub, VECDATA_DESC *pb,
             MATDATA_DESC *uuA,
             MATDATA_DESC *upA,
             MATDATA_DESC *puA,
             MATDATA_DESC *ppA)
{
  VECTOR *v;
  int v_cnt = 0;

  DOUBLE *ndata;

  ASSERT((NDATA_DEF_IN_GRID(theGrid) / 8) >= MAX_PATCH);

  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v))
  {
    INT ptype = VTYPE(v);
    INT pcomp = VD_NCMPS_IN_TYPE(pb,ptype);
    DOUBLE cor[MAX_PATCH];
    DOUBLE def[MAX_PATCH];
    VECTOR *w[MAX_VPATCH];
    MATRIX *m;
    INT ucomp[MAX_VPATCH];
    INT utype[MAX_VPATCH];
    INT n = pcomp;
    INT cnt = 0;
    INT i,j;

    if (pcomp == 0)
    {
      INT ncomp = VD_NCMPS_IN_TYPE(ub,ptype);
      DOUBLE mat[MAX_SINGLE_MAT_COMP];
      MATRIX *mm = START(v);

      for (i=0; i<ncomp; i++)
        def[i] = VVALUE(v,VD_CMP_OF_TYPE(ub,ptype,i));
      for (i=0; i<ncomp*ncomp; i++)
        mat[i] = MVALUE(mm,MD_MCMP_OF_RT_CT(uuA,ptype,ptype,i));
      for (mm = MNEXT(mm); mm!=NULL; mm=NEXT(mm)) {
        VECTOR *vv = MDEST(mm);
        for (j=0; j<VD_NCMPS_IN_TYPE(ux,VTYPE(vv)); j++)
          for (i=0; i<ncomp; i++)
            def[i] -=
              MVALUE(mm,MD_IJ_CMP_OF_RT_CT(uuA,ptype,
                                           VTYPE(vv),i,j))
              * VVALUE(vv,VD_CMP_OF_TYPE(ux,VTYPE(vv),j));
        for (j=0; j<VD_NCMPS_IN_TYPE(px,VTYPE(vv)); j++)
          for (i=0; i<ncomp; i++)
            def[i] -=
              MVALUE(mm,MD_IJ_CMP_OF_RT_CT(upA,ptype,
                                           VTYPE(vv),i,j))
              * VVALUE(vv,VD_CMP_OF_TYPE(px,VTYPE(vv),j));
      }
      SolveFullMatrix2(ncomp,cor,mat,def);
      for (i=0; i<ncomp; i++)
        VVALUE(v,VD_CMP_OF_TYPE(ux,ptype,i)) += cor[i];

      continue;
    }
    ASSERT(ptype==NODEVEC);
    ndata = (DOUBLE *)NDATA((NODE*)VOBJECT(v));

    /*
       if (me == 1) printf("ndata %d %x\n",VINDEX(v),ndata);
     */
    m = START(v);
    for (i=0; i<pcomp; i++) {
      def[i] = VVALUE(v,VD_CMP_OF_TYPE(pb,ptype,i));
      for (j=0; j<pcomp; j++)
        def[i] -=
          MVALUE(m,MD_IJ_CMP_OF_RT_CT(ppA,ptype,ptype,i,j))
          * VVALUE(v,VD_CMP_OF_TYPE(px,ptype,j));
    }
    for ( ; m!=NULL; m=NEXT(m))
    {
      MATRIX *mm;

      w[cnt] = MDEST(m);
      utype[cnt] = VTYPE(w[cnt]);
      ucomp[cnt] = VD_NCMPS_IN_TYPE(ux,utype[cnt]);

      if (ucomp[cnt] == 0) continue;
      if (n+ucomp[cnt] > MAX_PATCH) break;
      for (i=0; i<ucomp[cnt]; i++) {
        def[n+i] = VVALUE(w[cnt],VD_CMP_OF_TYPE(ub,utype[cnt],i));
        for (j=0; j<pcomp; j++)
          def[j] -=
            MVALUE(m,MD_IJ_CMP_OF_RT_CT(puA,ptype,utype[cnt],j,i))
            * VVALUE(w[cnt],VD_CMP_OF_TYPE(ux,utype[cnt],i));
      }
      for (i=0; i< VD_NCMPS_IN_TYPE(px,utype[cnt]); i++)
        for (j=0; j<pcomp; j++)
          def[j] -=
            MVALUE(m,MD_IJ_CMP_OF_RT_CT(ppA,ptype,utype[cnt],j,i))
            * VVALUE(w[cnt],VD_CMP_OF_TYPE(px,utype[cnt],i));
      for (mm=START(w[cnt]); mm!=NULL; mm=NEXT(mm)) {
        VECTOR *vv = MDEST(mm);
        for (j=0; j<VD_NCMPS_IN_TYPE(ux,VTYPE(vv)); j++)
          for (i=0; i<ucomp[cnt]; i++)
            def[n+i] -=
              MVALUE(mm,MD_IJ_CMP_OF_RT_CT(uuA,utype[cnt],
                                           VTYPE(vv),i,j))
              * VVALUE(vv,VD_CMP_OF_TYPE(ux,VTYPE(vv),j));
        for (j=0; j<VD_NCMPS_IN_TYPE(px,VTYPE(vv)); j++)
          for (i=0; i<ucomp[cnt]; i++)
            def[n+i] -=
              MVALUE(mm,MD_IJ_CMP_OF_RT_CT(upA,utype[cnt],
                                           VTYPE(vv),i,j))
              * VVALUE(vv,VD_CMP_OF_TYPE(px,VTYPE(vv),j));
      }
      n += ucomp[cnt];
      cnt++;
      if (cnt == MAX_VPATCH) break;
    }
    /*
       UserWrite("d ");
       for (i=0; i<n; i++)
        UserWriteF("%8.1e",def[i]);
       UserWriteF("s %d\n",VECSKIP(v));
     */

    /*mat[0][0] -= s;*/

    /*
       if (SolveFullMatrix3(n,cor,mat,def))
        RETURN(1);
     */
    /*
       if (InvertFullMatrix3(n,mat,ndata))
        RETURN(1);
     */
    for (i=0; i<n; i++) {
      DOUBLE s = 0.0;

      for (j=0; j<n; j++)
        s += ndata[n*i+j] * def[j];
      cor[i] = s;
    }
    for (n=0; n<pcomp; n++)
      VVALUE(v,VD_CMP_OF_TYPE(px,ptype,n)) += cor[n];
    for (i=0; i<cnt; i++)
      for (j=0; j<ucomp[i]; j++)
        VVALUE(w[i],VD_CMP_OF_TYPE(ux,utype[i],j)) += cor[n++];

    /*
       UserWrite("c ");
       for (i=0; i<n; i++)
        UserWriteF("%8.1e",cor[i]);
       UserWriteF("\n");
     */

    if (v_cnt++ == -2) return(0);
  }
  return(0);
}

static INT BLOCKInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_BLOCK *np = (NP_BLOCK *) theNP;

  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "u",argc,argv,&(np->u_sub));
  if (np->vt == NULL) {
    UserWriteF("BLOCKInit: no subtemplate u found\n");
    return(NP_NOT_ACTIVE);
  }
  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "p",argc,argv,&(np->p_sub));
  if (np->vt == NULL) {
    UserWriteF("BLOCKInit: no subtemplate p found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "uu",argc,argv,&(np->uu_sub));
  if (np->mt == NULL) {
    UserWriteF("BLOCKInit: no subtemplate uu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "up",argc,argv,&(np->up_sub));
  if (np->mt == NULL) {
    UserWriteF("BLOCKInit: no subtemplate up found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pu",argc,argv,&(np->pu_sub));
  if (np->mt == NULL) {
    UserWriteF("BLOCKInit: no subtemplate pu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pp",argc,argv,&(np->pp_sub));
  if (np->mt == NULL) {
    UserWriteF("BLOCKInit: no subtemplate pp found\n");
    return(NP_NOT_ACTIVE);
  }
  if (ReadArgvDOUBLE("s1",&(np->s1),argc,argv))
    np->s1 = 1.0;
  if (ReadArgvDOUBLE("s2",&(np->s2),argc,argv))
    np->s2 = 1.0;
  if (ReadArgvDOUBLE("s3",&(np->s3),argc,argv))
    np->s3 = 1.0;
  if (ReadArgvDOUBLE("s4",&(np->s4),argc,argv))
    np->s4 = 1.0;

  return (SmootherInit(theNP,argc,argv));
}

static INT BLOCKDisplay (NP_BASE *theNP)
{
  NP_BLOCK *np = (NP_BLOCK *) theNP;

  SmootherDisplay(theNP);
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->vt != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"vt",ENVITEM_NAME(np->vt));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"u_sub",(int)np->u_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"p_sub",(int)np->p_sub);
  if (np->ux != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ux",ENVITEM_NAME(np->ux));
  if (np->px != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"px",ENVITEM_NAME(np->px));
  if (np->ub != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ub",ENVITEM_NAME(np->ub));
  if (np->pb != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"pb",ENVITEM_NAME(np->pb));
  if (np->mt != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"vt",ENVITEM_NAME(np->mt));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"uu_sub",(int)np->uu_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"up_sub",(int)np->up_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"pu_sub",(int)np->pu_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"pp_sub",(int)np->pp_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"s1",np->s1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"s2",np->s2);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"s3",np->s3);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"s4",np->s4);
  if (np->uuA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"uuA",ENVITEM_NAME(np->uuA));
  if (np->upA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"upA",ENVITEM_NAME(np->upA));
  if (np->puA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"puA",ENVITEM_NAME(np->puA));
  if (np->ppA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ppA",ENVITEM_NAME(np->ppA));

  return (0);
}

static INT CountMultiplicity (GRID *g, VECDATA_DESC *x)
{
  ELEMENT *e;
  VECTOR *v;

  for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v))
    VINDEX(v) = 0;
  for (e = PFIRSTELEMENT(g); e != NULL; e = SUCCE(e))
  {
    VECTOR *vlist[MAX_NODAL_VECTORS];
    INT cnt = GetAllVectorsOfElementOfType(e,vlist,x);
    INT i;

    for (i=0; i<cnt; i++)
      VINDEX(vlist[i])++;
  }

  return(0);
}

static INT BLOCKPreProcess  (NP_ITER *theNP, INT level,
                             VECDATA_DESC *x, VECDATA_DESC *b,
                             MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_BLOCK *np = (NP_BLOCK *) theNP;
  GRID *theGrid = NP_GRID(theNP,level);

        #ifdef ModelP
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L))
    NP_RETURN(1,result[0]);
  if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A)
      != NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->smoother.L,MAT_CONS) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(np->smoother.L,np->mt,np->uu_sub,&np->uuA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(np->smoother.L,np->mt,np->up_sub,&np->upA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(np->smoother.L,np->mt,np->pu_sub,&np->puA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(np->smoother.L,np->mt,np->pp_sub,&np->ppA))
    NP_RETURN(1,result[0]);
        #else
  if (MDsubDescFromMT(A,np->mt,np->uu_sub,&np->uuA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->up_sub,&np->upA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pu_sub,&np->puA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pp_sub,&np->ppA))
    NP_RETURN(1,result[0]);
        #endif

  *baselevel = level;

  /* get storage for extra temp */
  if (AllocVDFromVD(NP_MG(theNP),level,level,x,&np->t))
    NP_RETURN(1,result[0]);

  /* CountMultiplicity(theGrid,x); */

  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->u_sub,&np->ub))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->p_sub,&np->pb))
    NP_RETURN(1,result[0]);

  /*
     if (l_block_check(theGrid,np->ux,np->px,np->ub,np->pb,
                                    np->uuA,np->upA,np->puA,np->ppA))
      NP_RETURN(1,result[0]);
   */

  if (l_block_collect(theGrid,np->ux,np->px,np->ub,np->pb,
                      np->uuA,np->upA,np->puA,np->ppA,
                      np->s1,np->s2,np->s3,np->s4))
    NP_RETURN(1,result[0]);

  return (0);
}

static INT BLOCKSmoother (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                          INT *result)
{
  NP_BLOCK *np = (NP_BLOCK *) theNP;
  GRID *theGrid = NP_GRID(theNP,level);

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->u_sub,&np->ub))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->p_sub,&np->pb))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->uu_sub,&np->uuA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->up_sub,&np->upA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pu_sub,&np->puA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pp_sub,&np->ppA))
    NP_RETURN(1,result[0]);

    #ifdef ModelP
  if (l_vector_collect(theGrid,b)!=NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (dset(NP_MG(theNP),level,level,ALL_VECTORS,x,0.0)!= NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_block(theGrid,np->ux,np->px,np->ub,np->pb,
              np->uuA,np->upA,np->puA,np->ppA))
    NP_RETURN(1,result[0]);
    #ifdef ModelP
  if (l_vector_consistent(theGrid,x) != NUM_OK)
    NP_RETURN(1,result[0]);
    #endif

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->smoother.damp)
      != NUM_OK)
    NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK)
    NP_RETURN(1,result[0]);

  /* PrintSVector(NP_MG(theNP),b); */

  return (0);
}

static INT BLOCKPostProcess (NP_ITER *theNP, INT level,
                             VECDATA_DESC *x, VECDATA_DESC *b,
                             MATDATA_DESC *A, INT *result)
{
  NP_BLOCK *np;

  np = (NP_BLOCK *) theNP;

  if (FreeVD(NP_MG(theNP),level,level,np->t)) REP_ERR_RETURN(1);

  return(SmootherPostProcess(theNP,level,x,b,A,result));
}

static INT BLOCKConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = BLOCKInit;
  theNP->Display = BLOCKDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = BLOCKPreProcess;
  np->Iter = BLOCKSmoother;
  np->PostProcess = BLOCKPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   ts - transforming smoother

   DESCRIPTION:
   This numproc executes a transforming smoother for
   a saddle point problem.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list> $u <ucomp> $p <vcomp>
       $UI <smoother1> $PI <smoother2>;
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

static INT II_Init (NP_BASE *theNP, INT argc , char **argv)
{
  NP_II *np;
  INT i;
  char post[VALUELEN],pre[VALUELEN],base[VALUELEN];

  np = (NP_II *) theNP;

  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->d = ReadArgvVecDesc(theNP->mg,"d",argc,argv);
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
      break;
    }

  if (ReadArgvINT("g",&(np->gamma),argc,argv))
    np->gamma = 1;
  if (ReadArgvINT("n1",&(np->nu1),argc,argv))
    np->nu1 = 1;
  if (ReadArgvINT("n2",&(np->nu2),argc,argv))
    np->nu2 = 1;
  if (ReadArgvINT("basen",&(np->basenu),argc,argv))
    np->basenu = 1;
  if (ReadArgvINT("b",&(np->baselevel),argc,argv))
    np->baselevel = 0;

  if (np->Transfer == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PreSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PostSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);

  if (sc_read(np->damp,NP_FMT(np),NULL,"damp",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->damp[i] = 1.0;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT II_Display (NP_BASE *theNP)
{
  NP_II *np;

  np = (NP_II *) theNP;

  NPIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->gamma);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)np->nu1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n2",(int)np->nu2);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"basen",(int)np->basenu);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);

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

  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->d != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"d",ENVITEM_NAME(np->d));

  return (0);
}

static INT II_PreProcess (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                          INT *baselevel, INT *result)
{
  NP_II *np;
  INT i;

  np = (NP_II *) theNP;

  if (np->Transfer->PreProcess != NULL)
    if ((*np->Transfer->PreProcess)
          (np->Transfer,&(np->baselevel),level,x,b,A,result))
      REP_ERR_RETURN(1);

  if (np->PreSmooth->PreProcess != NULL)
    for (i = level; i <= level; i++)
      if ((*np->PreSmooth->PreProcess)
            (np->PreSmooth,i,x,b,A,baselevel,result))
        REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PreProcess != NULL)
      for (i = level; i <= level; i++)
        if ((*np->PreSmooth->PreProcess)
              (np->PostSmooth,i,x,b,A,baselevel,result))
          REP_ERR_RETURN(1);

  *baselevel = MIN(np->baselevel,level);

  return (0);
}

static MATDATA_DESC *II_uuA,*II_upA,*II_puA,*II_ppA;
static VECDATA_DESC *II_t,*II_u,*II_p;
static NP_ITER *II_Smooth;

static INT II_MultiplySchurComplement (MULTIGRID *theMG, INT level,
                                       VECDATA_DESC *c, VECDATA_DESC *d,
                                       INT *result)
{
  if (dmatmul(theMG,level,level,ALL_VECTORS,II_t,II_upA,c) != NUM_OK)
    NP_RETURN(1,result[0]);
    #ifdef ModelP
  if (l_vector_meanvalue(GRID_ON_LEVEL(theMG,level),II_t)!=NUM_OK)
    NP_RETURN(1,result[0]);
    #endif
  if (dset(theMG,level,level,ALL_VECTORS,II_u,0.0)!= NUM_OK)
    NP_RETURN(1,result[0]);
  if ((*II_Smooth->Iter)(II_Smooth,level,II_u,II_t,II_uuA,result))
    REP_ERR_RETURN(1);
  if (dmatmul(theMG,level,level,ALL_VECTORS,d,II_puA,II_u)
      != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dmatmul_minus(theMG,level,level,ALL_VECTORS,d,II_ppA,c)
      != NUM_OK)
    NP_RETURN(1,result[0]);

  return(0);
}

static INT II_Iter (NP_ITER *theNP, INT level,
                    VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A,
                    INT *result)
{
  NP_II *np = (NP_II *) theNP;
  MULTIGRID *theMG = NP_MG(theNP);
  INT i;
  INT iter = np->nu1;

  if (AllocVDFromVD(theMG,level,level,c,&np->t)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,level,level,c,&np->d)) NP_RETURN(1,result[0]);
  if (dcopy(theMG,level,level,ALL_VECTORS,np->d,b) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dset(theMG,level,level,ALL_VECTORS,c,0.0) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (level == np->baselevel)
    iter = np->basenu - np->nu2;
  for (i=0; i<iter; i++) {
    if ((*np->PreSmooth->Iter)(np->PreSmooth,level,np->t,b,A,result))
      REP_ERR_RETURN(1);
    if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (II_MultiplySchurComplement(theMG,level,c,np->t,result))
      NP_RETURN(1,result[0]);
    if (dcopy(theMG,level,level,ALL_VECTORS,b,np->d) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (dadd(theMG,level,level,ALL_VECTORS,b,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  if (level > np->baselevel)
  {
    if ((*np->Transfer->RestrictDefect)
          (np->Transfer,level,b,b,A,Factor_One,result))
      REP_ERR_RETURN(1);

    if (dset(theMG,level-1,level-1,ALL_VECTORS,c,0.0) != NUM_OK)
      NP_RETURN(1,result[0]);
    for (i=0; i<np->gamma; i++)
      if (II_Iter(theNP,level-1,c,b,A,result))
        REP_ERR_RETURN(1);
    if ((*np->Transfer->InterpolateCorrection)
          (np->Transfer,level,np->t,c,A,np->damp,result))
      REP_ERR_RETURN(1);
    if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (II_MultiplySchurComplement(theMG,level,c,np->t,result))
      NP_RETURN(1,result[0]);
    if (dcopy(theMG,level,level,ALL_VECTORS,b,np->d) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (dadd(theMG,level,level,ALL_VECTORS,b,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  for (i=0; i<np->nu2; i++) {
    if ((*np->PostSmooth->Iter)(np->PostSmooth,level,np->t,b,A,result))
      REP_ERR_RETURN(1);
    if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
    /*	    if (i<np->nu2-1) {
        if (II_MultiplySchurComplement(theMG,level,c,np->t,result))
                NP_RETURN(1,result[0]);
            if (dcopy(theMG,level,level,ALL_VECTORS,b,np->d) != NUM_OK)
                NP_RETURN(1,result[0]);
            if (dadd(theMG,level,level,ALL_VECTORS,b,np->t) != NUM_OK)
                NP_RETURN(1,result[0]);
       } */
  }
  if (FreeVD(NP_MG(theNP),level,level,np->t)) REP_ERR_RETURN(1);
  if (FreeVD(NP_MG(theNP),level,level,np->d)) REP_ERR_RETURN(1);
  if (np->Transfer->AdaptCorrection != NULL)
    if ((*np->Transfer->AdaptCorrection)(np->Transfer,level,c,b,A,result))
      REP_ERR_RETURN(1);

  return (0);
}

static INT II_PostProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           INT *result)
{
  NP_II *np;
  INT i;

  np = (NP_II *) theNP;

  if (np->PreSmooth->PostProcess != NULL)
    for (i = level; i >= np->baselevel; i--)
      if ((*np->PreSmooth->PostProcess)
            (np->PreSmooth,i,x,b,A,result))
        REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PostProcess != NULL)
      for (i = level; i >= np->baselevel; i--)
        if ((*np->PreSmooth->PostProcess)
              (np->PostSmooth,i,x,b,A,result))
          REP_ERR_RETURN(1);

  if (np->Transfer->PostProcess != NULL)
    if ((*np->Transfer->PostProcess)
          (np->Transfer,&(np->baselevel),level,x,b,A,result))
      REP_ERR_RETURN(1);

  return (0);
}

static INT IIConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = II_Init;
  theNP->Display = II_Display;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = II_PreProcess;
  np->Iter = II_Iter;
  np->PostProcess = II_PostProcess;

  return(0);
}

static INT TSInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_TS *np;
  INT i;

  np = (NP_TS *) theNP;
  np->u = ReadArgvVecDesc(theNP->mg,"U",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->s = ReadArgvVecDesc(theNP->mg,"s",argc,argv);
  np->p = ReadArgvVecDesc(theNP->mg,"P",argc,argv);
  np->q = ReadArgvVecDesc(theNP->mg,"q",argc,argv);
  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  np->L = ReadArgvMatDesc(theNP->mg,"L",argc,argv);

  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "u",argc,argv,&(np->u_sub));
  if (np->vt == NULL) {
    UserWriteF("TSInit: no subtemplate u found\n");
    return(NP_NOT_ACTIVE);
  }
  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "p",argc,argv,&(np->p_sub));
  if (np->vt == NULL) {
    UserWriteF("TSInit: no subtemplate p found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "uu",argc,argv,&(np->uu_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate uu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "up",argc,argv,&(np->up_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate up found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pu",argc,argv,&(np->pu_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate pu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pp",argc,argv,&(np->pp_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate pp found\n");
    return(NP_NOT_ACTIVE);
  }
  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i] = 1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  for (i=0; i<MAX_VEC_COMP; i++) np->red[i] = 0.0;
  sc_read(np->red,NP_FMT(np),np->iter.b,"red",argc,argv);

  np->u_iter = (NP_ITER *)
               ReadArgvNumProc(theNP->mg,"UI",ITER_CLASS_NAME,argc,argv);
  if (np->u_iter == NULL) {
    np->u_solver = (NP_LINEAR_SOLVER *)
                   ReadArgvNumProc(theNP->mg,"UI",LINEAR_SOLVER_CLASS_NAME,argc,argv);
    if (np->u_solver == NULL) {
      UserWriteF("TSInit: no iter UI and no solver US found\n");
      return(NP_NOT_ACTIVE);
    }
  }
  else
    np->u_solver = NULL;
  np->v_iter = (NP_ITER *)
               ReadArgvNumProc(theNP->mg,"VI",ITER_CLASS_NAME,argc,argv);
  if (np->v_iter == NULL) np->v_iter = np->u_iter;
  np->p_iter = (NP_ITER *)
               ReadArgvNumProc(theNP->mg,"PI",ITER_CLASS_NAME,argc,argv);
  if (np->p_iter == NULL) {
    UserWriteF("TSInit: no iter PI found\n");
    return(NP_NOT_ACTIVE);
  }
    #ifdef ModelP
  if (ReadArgvOption("M",argc,argv))
    np->cons_mode = MAT_CONS;
  else if (ReadArgvOption("D",argc,argv))
    np->cons_mode = MAT_DIAG_CONS;
  else
    np->cons_mode = MAT_MASTER_CONS;
        #endif
  if (ReadArgvINT("dc",&np->dc,argc,argv))
    np->dc = 0;
  np->extra = ReadArgvOption("extra",argc,argv);
  np->ls = ReadArgvOption("ls",argc,argv);
  np->diag = ReadArgvOption("diag",argc,argv);
  np->display = ReadArgvDisplay(argc,argv);
  np->dc_max = 0;
  if (ReadArgvDOUBLE("thresh",&np->thresh,argc,argv))
    np->thresh = 0.0;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT TSDisplay (NP_BASE *theNP)
{
  NP_TS *np;

  np = (NP_TS *) theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);
  if (sc_disp(np->red,np->iter.b,"red")) REP_ERR_RETURN (1);
    #ifdef ModelP
  UserWriteF(DISPLAY_NP_FORMAT_SI,"cons_mode",(int)np->cons_mode);
        #endif
  if (np->u_iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"UI",ENVITEM_NAME(np->u_iter));
  if (np->u_solver != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"UI",ENVITEM_NAME(np->u_solver));
  if (np->v_iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"VI",ENVITEM_NAME(np->v_iter));
  if (np->p_iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"PI",ENVITEM_NAME(np->p_iter));
  if (np->u != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"U",ENVITEM_NAME(np->u));
  if (np->r != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->s != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(np->s));
  if (np->p != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"P",ENVITEM_NAME(np->p));
  if (np->q != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"q",ENVITEM_NAME(np->q));
  if (np->L != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(np->L));
  if (np->S != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(np->S));
  if (np->vt != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"vt",ENVITEM_NAME(np->vt));
  if (np->mt != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"mt",ENVITEM_NAME(np->mt));
  if (np->ux != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ux",ENVITEM_NAME(np->ux));
  if (np->px != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"px",ENVITEM_NAME(np->px));
  if (np->ub != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ub",ENVITEM_NAME(np->ub));
  if (np->pb != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"pb",ENVITEM_NAME(np->pb));
  if (np->uuA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"uuA",ENVITEM_NAME(np->uuA));
  if (np->puA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"puA",ENVITEM_NAME(np->puA));
  if (np->upA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"upA",ENVITEM_NAME(np->upA));
  if (np->ppA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ppA",ENVITEM_NAME(np->ppA));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"u_sub",(int)np->u_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"p_sub",(int)np->p_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"uu_sub",(int)np->uu_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"up_sub",(int)np->up_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"pu_sub",(int)np->pu_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"pp_sub",(int)np->pp_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"dc",(int)np->dc);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"dc_max",(int)np->dc_max);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"extra",(int)np->extra);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"ls",(int)np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"diag",(int)np->diag);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"thresh",(float)np->thresh);

  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  return (0);
}

static INT ConstructSchurComplement (GRID *theGrid,
                                     MATDATA_DESC *uuA,
                                     MATDATA_DESC *upA,
                                     MATDATA_DESC *puA,
                                     MATDATA_DESC *ppA,
                                     MATDATA_DESC *S, INT extra)
{
  VECTOR *v,*w,*z;
  MATRIX *m,*m0,*m1;
  DOUBLE *sval,InvMat[MAX_SINGLE_MAT_COMP],sum,*puval,*upval;
  INT i,j,k,l,vtype,wtype,ztype,vncomp,wncomp,zncomp;
  SHORT *pu,*up,*s;

  if (dmatcopy(MYMG(theGrid),GLEVEL(theGrid),GLEVEL(theGrid),
               ALL_VECTORS,S,ppA) != NUM_OK)
    REP_ERR_RETURN (1);

  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    vncomp = MD_ROWS_IN_RT_CT(puA,vtype,vtype);
    if (vncomp == 0) continue;
    if (extra) {
      w = FIRSTVECTOR(theGrid);
    }
    else {
      m = VSTART(v);
      if (m != NULL)
        w = MDEST(m);
      else
        w = NULL;
    }
    while (w != NULL) {
      wtype = VTYPE(w);
      wncomp = MD_COLS_IN_RT_CT(upA,vtype,wtype);
      if (wncomp > 0) {
        if (extra) {
          m = GetMatrix(v,w);
          if (m == NULL)
            m = CreateExtraConnection(theGrid,v,w);
          ASSERT(m != NULL);
        }
        s = MD_MCMPPTR_OF_RT_CT(S,vtype,wtype);
        sval = MVALUEPTR(m,0);
        for (m0=VSTART(v); m0!=NULL; m0=NEXT(m0)) {
          z = MDEST(m0);
          ztype = VTYPE(z);
          zncomp = MD_COLS_IN_RT_CT(puA,vtype,ztype);
          if (zncomp == 0) continue;
          for (m1=VSTART(w); m1!=NULL; m1=NEXT(m1))
            if (MDEST(m1) == z) break;
          if (m1 == NULL) continue;
          ASSERT(zncomp == MD_COLS_IN_RT_CT(uuA,ztype,ztype));
          if (InvertSmallBlock(zncomp,
                               MD_MCMPPTR_OF_RT_CT(uuA,ztype,ztype),
                               MVALUEPTR(VSTART(z),0),InvMat)) {
            /*   REP_ERR_RETURN (1);  */
            for (i=0; i<zncomp*zncomp; i++) InvMat[i] = 0.0;
            for (i=0; i<zncomp; i++) InvMat[i*zncomp+i] = 1.0;
          }
          pu = MD_MCMPPTR_OF_RT_CT(puA,vtype,ztype);
          up = MD_MCMPPTR_OF_RT_CT(puA,wtype,ztype);
          puval = MVALUEPTR(m0,0);
          upval = MVALUEPTR(m1,0);
          for (i=0; i<vncomp; i++)
            for (j=0; j<wncomp; j++) {
              sum = 0.0;
              for (k=0; k<zncomp; k++)
                for (l=0; l<zncomp; l++)
                  sum += puval[pu[i*zncomp+k]]
                         * InvMat[k*zncomp+l]
                         * upval[up[j*zncomp+l]];
              sval[s[i*wncomp+j]] -= sum;
            }
        }
      }
      if (extra)
        w = SUCCVC(w);
      else {
        m = NEXT(m);
        if (m != NULL)
          w = MDEST(m);
        else
          w = NULL;
      }
    }
  }

  return (0);
}

static INT ConstructDiagSchurComplement (GRID *theGrid,
                                         MATDATA_DESC *uuA,
                                         MATDATA_DESC *upA,
                                         MATDATA_DESC *puA,
                                         MATDATA_DESC *ppA,
                                         MATDATA_DESC *S)
{
  VECTOR *v;

  if (dmatcopy(MYMG(theGrid),GLEVEL(theGrid),GLEVEL(theGrid),
               ALL_VECTORS,S,ppA) != NUM_OK)
    REP_ERR_RETURN (1);

  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v))
  {
    INT vtype = VTYPE(v);
    INT vncomp = MD_ROWS_IN_RT_CT(ppA,vtype,vtype);

    if (vncomp > 0)
    {
      SHORT *s = MD_MCMPPTR_OF_RT_CT(S,vtype,vtype);
      DOUBLE *sval = MVALUEPTR(VSTART(v),0);
      MATRIX *m;

      for (m=VSTART(v); m!=NULL; m=NEXT(m))
      {
        INT wtype = MDESTTYPE(m);
        INT wncomp = MD_COLS_IN_RT_CT(puA,vtype,wtype);

        if (wncomp > 0)
        {
          VECTOR *w = MDEST(m);
          SHORT *pu = MD_MCMPPTR_OF_RT_CT(puA,vtype,wtype);
          DOUBLE *puval = MVALUEPTR(m,0);
          SHORT *up = MD_MCMPPTR_OF_RT_CT(upA,wtype,vtype);
          DOUBLE *upval = MVALUEPTR(MADJ(m),0);
          DOUBLE InvMat[MAX_SINGLE_MAT_COMP];
          INT i,j,k,l;

          if (InvertSmallBlock(wncomp,
                               MD_MCMPPTR_OF_RT_CT(uuA,wtype,wtype),
                               MVALUEPTR(VSTART(w),0),InvMat)) {
            for (i=0; i<wncomp*wncomp; i++) InvMat[i] = 0.0;
            for (i=0; i<wncomp; i++) InvMat[i*wncomp+i] = 1.0;
          }
          for (i=0; i<vncomp; i++)
            for (j=0; j<vncomp; j++)
            {
              DOUBLE sum = 0.0;

              for (k=0; k<wncomp; k++)
                for (l=0; l<wncomp; l++) {
                  sum += puval[pu[i*wncomp+k]]
                         * InvMat[k*wncomp+l]
                         * upval[up[l*vncomp+j]];
                }
              sval[s[i*vncomp+j]] -= sum;
            }
        }
      }
    }
  }

  return (0);
}

static INT ConstructExtraDiagSchurComplement (GRID *theGrid,
                                              MATDATA_DESC *uuA,
                                              MATDATA_DESC *upA,
                                              MATDATA_DESC *puA,
                                              MATDATA_DESC *ppA,
                                              MATDATA_DESC *S)
{
  VECTOR *w;

  if (dmatcopy(MYMG(theGrid),GLEVEL(theGrid),GLEVEL(theGrid),
               ALL_VECTORS,S,ppA) != NUM_OK)
    REP_ERR_RETURN (1);

  for (w=FIRSTVECTOR(theGrid); w!=NULL; w=SUCCVC(w))
  {
    INT wtype = VTYPE(w);
    INT wncomp = MD_ROWS_IN_RT_CT(uuA,wtype,wtype);
    DOUBLE InvMat[MAX_SINGLE_MAT_COMP];
    INT i,j,k,l;
    MATRIX *m0;

    if (wncomp == 0) continue;

    if (InvertSmallBlock(wncomp,
                         MD_MCMPPTR_OF_RT_CT(uuA,wtype,wtype),
                         MVALUEPTR(VSTART(w),0),InvMat)) {
      for (i=0; i<wncomp*wncomp; i++) InvMat[i] = 0.0;
      for (i=0; i<wncomp; i++) InvMat[i*wncomp+i] = 1.0;
    }
    for (m0=VSTART(w); m0!=NULL; m0=NEXT(m0))
    {
      VECTOR *v = MDEST(m0);
      INT vtype = VTYPE(v);
      INT vncomp = MD_ROWS_IN_RT_CT(puA,vtype,wtype);
      SHORT *pu = MD_MCMPPTR_OF_RT_CT(puA,vtype,wtype);
      DOUBLE *puval = MVALUEPTR(MADJ(m0),0);
      MATRIX *m1;

      if (vncomp == 0) continue;

      for (m1=VSTART(w); m1!=NULL; m1=NEXT(m1))
      {
        VECTOR *z = MDEST(m1);
        INT ztype = VTYPE(z);
        INT zncomp = MD_COLS_IN_RT_CT(upA,wtype,ztype);
        SHORT *up = MD_MCMPPTR_OF_RT_CT(upA,wtype,ztype);
        DOUBLE *upval = MVALUEPTR(m1,0);
        SHORT *s = MD_MCMPPTR_OF_RT_CT(S,vtype,ztype);
        DOUBLE *sval;
        MATRIX *m = GetMatrix(v,z);

        if (zncomp == 0) continue;
        if (m == NULL) {
          m = CreateExtraConnection(theGrid,v,z);
          ASSERT(m != NULL);
        }
        sval = MVALUEPTR(m,0);

        for (i=0; i<vncomp; i++)
          for (j=0; j<zncomp; j++)
          {
            DOUBLE sum = 0.0;

            for (k=0; k<wncomp; k++)
              for (l=0; l<wncomp; l++) {
                sum += puval[pu[i*wncomp+k]]
                       * InvMat[k*wncomp+l]
                       * upval[up[l*zncomp+j]];
              }
            sval[s[i*zncomp+j]] -= sum;
          }
      }
    }
  }

  return (0);
}

static INT TSPreProcess  (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b,
                          MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_TS *np;
  GRID *theGrid;

  np = (NP_TS *) theNP;
  theGrid = NP_GRID(theNP,level);
  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->u_sub,&np->ub))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->p_sub,&np->pb))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->uu_sub,&np->uuA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->up_sub,&np->upA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pu_sub,&np->puA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pp_sub,&np->ppA))
    NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (AllocMDFromMD(NP_MG(theNP),level,level,np->uuA,&np->L))
    NP_RETURN(1,result[0]);
  if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,np->uuA)
      != NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->L,np->cons_mode) != NUM_OK)
    NP_RETURN(1,result[0]);
        #else
  np->L = np->uuA;
        #endif
  if (AllocMDFromMD(NP_MG(theNP),level,level,np->ppA,&np->S))
    NP_RETURN(1,result[0]);
  /*	if (AssembleTotalDirichletBoundary(theGrid,A,x,b))
      NP_RETURN(1,result[0]); */

  if (np->diag) {
    if (np->extra) {
      if (ConstructExtraDiagSchurComplement(theGrid,np->L,
                                            np->upA,np->puA,np->ppA,
                                            np->S))
        NP_RETURN(1,result[0]);
    }
    else {
      if (ConstructDiagSchurComplement(theGrid,np->L,
                                       np->upA,np->puA,np->ppA,
                                       np->S))
        NP_RETURN(1,result[0]);
    }
  }
  else {
    if (ConstructSchurComplement(theGrid,np->L,np->upA,np->puA,np->ppA,
                                 np->S,np->extra))
      NP_RETURN(1,result[0]);
  }
        #ifdef ModelP
  FreeMD(NP_MG(theNP),level,level,np->L);
        #endif
  *baselevel = level;

  if (np->u_iter != NULL)
    if (np->u_iter->PreProcess != NULL)
      if ((*np->u_iter->PreProcess)
            (np->u_iter,level,np->ux,np->ub,np->uuA,baselevel,result))
        REP_ERR_RETURN(1);
  if (np->u_solver != NULL)
    if (np->u_solver->PreProcess != NULL)
      if ((*np->u_solver->PreProcess)
            (np->u_solver,*baselevel,np->ux,np->ub,np->uuA,
            baselevel,result))
        REP_ERR_RETURN(1);
  if (np->v_iter != np->u_iter)
    if (np->v_iter->PreProcess != NULL)
      if ((*np->v_iter->PreProcess)
            (np->v_iter,level,np->ux,np->ub,np->uuA,baselevel,result))
        REP_ERR_RETURN(1);
  if (np->p_iter->PreProcess != NULL)
    if ((*np->p_iter->PreProcess)
          (np->p_iter,level,np->px,np->pb,np->S,baselevel,result))
      REP_ERR_RETURN(1);

  np->dc_max = 0;

  return (0);
}

static INT TSSmoother (NP_ITER *theNP, INT level,
                       VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                       INT *result)
{
  NP_TS *np;
  MULTIGRID *theMG;
  VEC_SCALAR defect2reach,defect;
  DOUBLE rho,lambda;
  INT i,PrintID;
  char text[DISPLAY_WIDTH+4],text1[DISPLAY_WIDTH];
  LRESULT lresult;

  np = (NP_TS *) theNP;
  theMG = NP_MG(theNP);

  /* get storage for extra temp */
  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->u_sub,&np->ub))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->p_sub,&np->pb))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,level,level,np->ux,&np->u))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,level,level,np->ux,&np->t))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,level,level,np->px,&np->s))
    NP_RETURN(1,result[0]);
  if (np->dc) {
    if (AllocVDFromVD(theMG,level,level,np->px,&np->q))
      NP_RETURN(1,result[0]);
    if (AllocVDFromVD(theMG,level,level,np->px,&np->r))
      NP_RETURN(1,result[0]);
    if (AllocVDFromVD(theMG,level,level,np->px,&np->p))
      NP_RETURN(1,result[0]);
  }
  if (dcopy(theMG,level,level,ALL_VECTORS,np->t,np->ub) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dcopy(theMG,level,level,ALL_VECTORS,np->s,np->pb) != NUM_OK)
    NP_RETURN(1,result[0]);
    #ifdef ModelP
  if (l_vector_meanvalue(GRID_ON_LEVEL(theMG,level),np->t)!=NUM_OK)
    NP_RETURN(1,result[0]);
    #endif
  if (dset(theMG,level,level,ALL_VECTORS,x,0.0)!= NUM_OK)
    NP_RETURN(1,result[0]);

  II_uuA = np->uuA;
  II_puA = np->puA;
  II_upA = np->upA;
  II_ppA = np->ppA;
  II_t = np->t;
  II_u = np->u;
  II_p = np->p;
  II_Smooth = np->v_iter;

  if (np->u_iter != NULL) {
    if ((*np->u_iter->Iter)(np->u_iter,level,np->ux,np->t,np->uuA,result))
      REP_ERR_RETURN(1);
  }
  else {
    if ((*np->u_solver->Residuum)
          (np->u_solver,0,level,
          np->ux,np->t,np->uuA,&lresult))
      NP_RETURN(1,result[0]);
    if ((*np->u_solver->Solver)(np->u_solver,level,
                                np->ux,np->t,np->uuA,
                                np->u_solver->abslimit,
                                np->u_solver->reduction,&lresult))
      NP_RETURN(1,result[0]);
  }
  if (dmatmul_minus(theMG,level,level,ALL_VECTORS,np->s,np->puA,np->ux)
      != NUM_OK)
    NP_RETURN(1,result[0]);
    #ifdef ModelP
  if (l_vector_meanvalue(GRID_ON_LEVEL(theMG,level),np->s)!=NUM_OK)
    NP_RETURN(1,result[0]);
    #endif

  if (np->dc)
  {
    INT display = np->display;

    if (level < TOPLEVEL(theMG))
      display = PCR_NO_DISPLAY;
    if (dcopy(theMG,level,level,ALL_VECTORS,np->q,np->s) != NUM_OK)
      NP_RETURN(1,result[0]);
    sprintf(text1,"level %d: %s",level,ENVITEM_NAME(np));
    CenterInPattern(text,DISPLAY_WIDTH,text1,' ',"\n");
    if (PreparePCR(np->q,display,text,&PrintID))
      NP_RETURN(1,result[0]);
        #ifdef ModelP
    if (l_vector_collect(GRID_ON_LEVEL(theMG,level),np->q)!=NUM_OK)
      NP_RETURN(1,result[0]);
        #endif
    if (dnrm2x(theMG,level,level,ALL_VECTORS,np->q,defect))
      NP_RETURN(1,result[0]);
    if (sc_mul_check(defect2reach,defect,np->red,np->q))
      NP_RETURN(1,result[0]);
    if (DoPCR(PrintID,defect,PCR_CRATE_SD))
      NP_RETURN(1,result[0]);
    if (dset(theMG,level,level,ALL_VECTORS,np->p,0.0)!= NUM_OK)
      NP_RETURN(1,result[0]);
    rho = 1.0;
  }
  else {
    np->q = np->s;
    np->r = np->px;
  }
  if (dset(theMG,level,level,ALL_VECTORS,np->r,0.0)!= NUM_OK)
    NP_RETURN(1,result[0]);
  if ((*np->p_iter->Iter)(np->p_iter,level,np->r,np->q,np->S,result))
    REP_ERR_RETURN(1);
  /* defect correction for the Schur complement*/
  for (i=0; i<np->dc; i++) {

    /*
          if (ddot(theMG,level,level,ALL_VECTORS,np->s,np->s,&lambda))
              NP_RETURN(1,result[0]);
          printf("s %e\n",lambda);
          if (ddot(theMG,level,level,ALL_VECTORS,np->q,np->q,&lambda))
              NP_RETURN(1,result[0]);
          printf("q %e\n",lambda);
          if (ddot(theMG,level,level,ALL_VECTORS,np->r,np->r,&lambda))
              NP_RETURN(1,result[0]);
          printf("r %e\n",lambda);
     */

    /* compute q = S * p */
    if (ddot(theMG,level,level,ALL_VECTORS,np->s,np->r,&lambda))
      NP_RETURN(1,result[0]);
    rho = - lambda / rho;
    if (np->ls) rho = 0.0;
    PRINTDEBUG(np,1,("rho %e lambda %e\n",rho,lambda));
    if (daxpy(theMG,level,level,ALL_VECTORS,np->r,rho,np->p) != NUM_OK)
      NP_RETURN(1,result[0]);
    rho = - lambda;
    ASSERT(lambda != 0.0);
    if (dcopy(theMG,level,level,ALL_VECTORS,np->p,np->r) != NUM_OK)
      NP_RETURN(1,result[0]);

    /* Multiply Schur complement */
    if (dmatmul(theMG,level,level,ALL_VECTORS,np->t,np->upA,np->p)
        != NUM_OK)
      NP_RETURN(1,result[0]);
        #ifdef ModelP
    if (l_vector_meanvalue(GRID_ON_LEVEL(theMG,level),np->t)!=NUM_OK)
      NP_RETURN(1,result[0]);
        #endif
    if (dset(theMG,level,level,ALL_VECTORS,np->u,0.0)!= NUM_OK)
      NP_RETURN(1,result[0]);

    /*
       if ((*np->v_iter->Iter)(np->v_iter,level,np->u,np->t,np->uuA,result))
        REP_ERR_RETURN(1);
     */

    if (np->v_iter != NULL) {
      if ((*np->v_iter->Iter)(np->v_iter,level,
                              np->u,np->t,np->uuA,result))
        REP_ERR_RETURN(1);
    }
    else {
      if ((*np->u_solver->Residuum)
            (np->u_solver,0,level,
            np->u,np->t,np->uuA,&lresult))
        NP_RETURN(1,result[0]);
      if ((*np->u_solver->Solver)(np->u_solver,level,
                                  np->u,np->t,np->uuA,
                                  np->u_solver->abslimit,
                                  np->u_solver->reduction,&lresult))
        NP_RETURN(1,result[0]);
    }


    if (dmatmul(theMG,level,level,ALL_VECTORS,np->q,np->puA,np->u)
        != NUM_OK)
      NP_RETURN(1,result[0]);
    if (dmatmul_minus(theMG,level,level,ALL_VECTORS,np->q,np->ppA,np->p)
        != NUM_OK)
      NP_RETURN(1,result[0]);
    /* for testing ...
       II_MultiplySchurComplement(theMG,level,np->p,np->q,result);
     */
    if (ddot(theMG,level,level,ALL_VECTORS,np->q,np->p,&lambda))
      NP_RETURN(1,result[0]);
    ASSERT(lambda != 0);
    lambda = rho / lambda;
    if (np->ls) lambda = 1.0;
    /* update inner defect */
    if (daxpy(theMG,level,level,ALL_VECTORS,np->s,lambda,np->q) != NUM_OK)
      NP_RETURN(1,result[0]);
    /* update Lagrange corrector */
    if (daxpy(theMG,level,level,ALL_VECTORS,np->px,lambda,np->p) != NUM_OK)
      NP_RETURN(1,result[0]);
        #ifdef ModelP
    if (l_vector_collect(GRID_ON_LEVEL(theMG,level),np->s)!=NUM_OK)
      NP_RETURN(1,result[0]);
        #endif
    if (dnrm2x(theMG,level,level,ALL_VECTORS,np->s,defect))
      NP_RETURN(1,result[0]);
    if (DoPCR(PrintID,defect,PCR_CRATE_SD))
      NP_RETURN(1,result[0]);
    if (sc_cmp(defect,defect2reach,np->q))
      break;
        #ifdef ModelP
    if (l_vector_meanvalue(GRID_ON_LEVEL(theMG,level),np->q)!=NUM_OK)
      NP_RETURN(1,result[0]);
        #endif
    if (dset(theMG,level,level,ALL_VECTORS,np->r,0.0)!= NUM_OK)
      NP_RETURN(1,result[0]);
    if (dcopy(theMG,level,level,ALL_VECTORS,np->q,np->s) != NUM_OK)
      NP_RETURN(1,result[0]);
    if ((*np->p_iter->Iter)(np->p_iter,level,np->r,np->q,np->S,result))
      REP_ERR_RETURN(1);
  }
  if (np->dc) {
    np->dc_max = MAX(i+1,np->dc_max);
    if (DoPCR(PrintID,defect,PCR_AVERAGE_SD))
      NP_RETURN(1,result[0]);
    if (PostPCR(PrintID,NULL))
      NP_RETURN(1,result[0]);
    if (dscal(theMG,level,level,ALL_VECTORS,np->px,1.0/(1.0-np->red[0]))
        != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  if (dmatmul(theMG,level,level,ALL_VECTORS,np->t,np->upA,np->px)
      != NUM_OK)
    NP_RETURN(1,result[0]);

    #ifdef ModelP
  if (l_vector_meanvalue(GRID_ON_LEVEL(theMG,level),np->t)!=NUM_OK)
    NP_RETURN(1,result[0]);
    #endif
  if (dset(theMG,level,level,ALL_VECTORS,np->u,0.0)!= NUM_OK)
    NP_RETURN(1,result[0]);
  if (np->v_iter != NULL) {
    if ((*np->v_iter->Iter)(np->v_iter,level,
                            np->u,np->t,np->uuA,result))
      REP_ERR_RETURN(1);
  }
  else {
    if ((*np->u_solver->Residuum)
          (np->u_solver,0,level,
          np->u,np->t,np->uuA,&lresult))
      NP_RETURN(1,result[0]);
    if ((*np->u_solver->Solver)(np->u_solver,level,
                                np->u,np->t,np->uuA,
                                np->u_solver->abslimit,
                                np->u_solver->reduction,&lresult))
      NP_RETURN(1,result[0]);
  }
  if (dsub(theMG,level,level,ALL_VECTORS,np->ux,np->u) != NUM_OK)
    NP_RETURN(1,result[0]);

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->damp)!= NUM_OK)
    NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK)
    NP_RETURN(1,result[0]);

  FreeVD(NP_MG(theNP),level,level,np->u);
  FreeVD(NP_MG(theNP),level,level,np->t);
  FreeVD(NP_MG(theNP),level,level,np->s);
  if (np->dc) {
    FreeVD(NP_MG(theNP),level,level,np->p);
    FreeVD(NP_MG(theNP),level,level,np->q);
    FreeVD(NP_MG(theNP),level,level,np->r);
  }
  return (0);
}

static INT TSPostProcess (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b,
                          MATDATA_DESC *A, INT *result)
{
  NP_TS *np;

  np = (NP_TS *) theNP;
  FreeMD(NP_MG(theNP),level,level,np->S);
  if (np->u_iter != NULL)
    if (np->u_iter->PostProcess != NULL)
      if ((*np->u_iter->PostProcess)
            (np->u_iter,level,x,b,A,result))
        REP_ERR_RETURN(1);
  if (np->u_solver != NULL)
    if (np->u_solver->PostProcess != NULL)
      if ((*np->u_solver->PostProcess)
            (np->u_solver,level,x,b,A,result))
        REP_ERR_RETURN(1);
  if (np->v_iter != np->u_iter)
    if (np->v_iter->PostProcess != NULL)
      if ((*np->v_iter->PostProcess)
            (np->v_iter,level,x,b,A,result))
        REP_ERR_RETURN(1);
  if (np->p_iter->PostProcess != NULL)
    if ((*np->p_iter->PostProcess)
          (np->p_iter,level,x,b,A,result))
      REP_ERR_RETURN(1);

  if (TOPLEVEL(NP_MG(theNP)) == level) {
    if (np->display > PCR_NO_DISPLAY)
      UserWriteF("maximal number of inner iterations: %d\n",np->dc_max);
    if (SetStringValue(":iter:inner",(DOUBLE)np->dc_max))
      NP_RETURN(1,result[0]);
  }

  return(0);
}

static INT TSConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = TSInit;
  theNP->Display = TSDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = TSPreProcess;
  np->Iter = TSSmoother;
  np->PostProcess = TSPostProcess;

  return(0);
}

static INT APInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_TS *np;
  INT i;

  np = (NP_TS *) theNP;
  np->u = ReadArgvVecDesc(theNP->mg,"U",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->s = ReadArgvVecDesc(theNP->mg,"s",argc,argv);
  np->p = ReadArgvVecDesc(theNP->mg,"P",argc,argv);
  np->q = ReadArgvVecDesc(theNP->mg,"q",argc,argv);
  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  np->L = ReadArgvMatDesc(theNP->mg,"L",argc,argv);

  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "u",argc,argv,&(np->u_sub));
  if (np->vt == NULL) {
    UserWriteF("TSInit: no subtemplate u found\n");
    return(NP_NOT_ACTIVE);
  }
  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "p",argc,argv,&(np->p_sub));
  if (np->vt == NULL) {
    UserWriteF("TSInit: no subtemplate p found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "uu",argc,argv,&(np->uu_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate uu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "up",argc,argv,&(np->up_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate up found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pu",argc,argv,&(np->pu_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate pu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pp",argc,argv,&(np->pp_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate pp found\n");
    return(NP_NOT_ACTIVE);
  }
  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i] = 1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  for (i=0; i<MAX_VEC_COMP; i++) np->red[i] = 0.0;
  sc_read(np->red,NP_FMT(np),np->iter.b,"red",argc,argv);

  np->u_iter = (NP_ITER *)
               ReadArgvNumProc(theNP->mg,"UI",ITER_CLASS_NAME,argc,argv);
  if (np->u_iter == NULL) {
    np->u_solver = (NP_LINEAR_SOLVER *)
                   ReadArgvNumProc(theNP->mg,"UI",LINEAR_SOLVER_CLASS_NAME,argc,argv);
    if (np->u_solver == NULL) {
      UserWriteF("TSInit: no iter UI and no solver UI found\n");
      return(NP_NOT_ACTIVE);
    }
    np->p_solver = (NP_LINEAR_SOLVER *)
                   ReadArgvNumProc(theNP->mg,"PI",LINEAR_SOLVER_CLASS_NAME,argc,argv);
    if (np->p_solver == NULL) {
      UserWriteF("TSInit: no iter PI and no solver PI found\n");
      return(NP_NOT_ACTIVE);
    }
  }
  else
    np->u_solver = NULL;
  np->v_iter = (NP_ITER *)
               ReadArgvNumProc(theNP->mg,"VI",ITER_CLASS_NAME,argc,argv);
  if (np->v_iter == NULL) np->v_iter = np->u_iter;
  np->p_iter = (NP_ITER *)
               ReadArgvNumProc(theNP->mg,"PI",ITER_CLASS_NAME,argc,argv);
    #ifdef ModelP
  if (ReadArgvOption("M",argc,argv))
    np->cons_mode = MAT_CONS;
  else if (ReadArgvOption("D",argc,argv))
    np->cons_mode = MAT_DIAG_CONS;
  else
    np->cons_mode = MAT_MASTER_CONS;
        #endif
  if (ReadArgvINT("dc",&np->dc,argc,argv))
    np->dc = 0;
  np->extra = ReadArgvOption("extra",argc,argv);
  np->ls = ReadArgvOption("ls",argc,argv);
  np->diag = ReadArgvOption("diag",argc,argv);
  np->display = ReadArgvDisplay(argc,argv);
  np->dc_max = 0;
  if (ReadArgvDOUBLE("thresh",&np->thresh,argc,argv))
    np->thresh = 0.0;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT APDisplay (NP_BASE *theNP)
{
  NP_TS *np;

  np = (NP_TS *) theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);
  if (sc_disp(np->red,np->iter.b,"red")) REP_ERR_RETURN (1);
    #ifdef ModelP
  UserWriteF(DISPLAY_NP_FORMAT_SI,"cons_mode",(int)np->cons_mode);
        #endif
  if (np->u_iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"UI",ENVITEM_NAME(np->u_iter));
  if (np->u_solver != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"UI",ENVITEM_NAME(np->u_solver));
  if (np->p_solver != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"UI",ENVITEM_NAME(np->p_solver));
  if (np->v_iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"VI",ENVITEM_NAME(np->v_iter));
  if (np->p_iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"PI",ENVITEM_NAME(np->p_iter));
  if (np->u != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"U",ENVITEM_NAME(np->u));
  if (np->r != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->s != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(np->s));
  if (np->p != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"P",ENVITEM_NAME(np->p));
  if (np->q != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"q",ENVITEM_NAME(np->q));
  if (np->L != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(np->L));
  if (np->S != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"S",ENVITEM_NAME(np->S));
  if (np->vt != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"vt",ENVITEM_NAME(np->vt));
  if (np->mt != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"mt",ENVITEM_NAME(np->mt));
  if (np->ux != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ux",ENVITEM_NAME(np->ux));
  if (np->px != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"px",ENVITEM_NAME(np->px));
  if (np->ub != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ub",ENVITEM_NAME(np->ub));
  if (np->pb != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"pb",ENVITEM_NAME(np->pb));
  if (np->uuA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"uuA",ENVITEM_NAME(np->uuA));
  if (np->puA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"puA",ENVITEM_NAME(np->puA));
  if (np->upA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"upA",ENVITEM_NAME(np->upA));
  if (np->ppA != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"ppA",ENVITEM_NAME(np->ppA));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"u_sub",(int)np->u_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"p_sub",(int)np->p_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"uu_sub",(int)np->uu_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"up_sub",(int)np->up_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"pu_sub",(int)np->pu_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"pp_sub",(int)np->pp_sub);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"dc",(int)np->dc);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"dc_max",(int)np->dc_max);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"extra",(int)np->extra);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"ls",(int)np->ls);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"diag",(int)np->diag);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"thresh",(float)np->thresh);

  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  return (0);
}

static INT APPreProcess  (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b,
                          MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_TS *np = (NP_TS *) theNP;

  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->u_sub,&np->ub))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->p_sub,&np->pb))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->uu_sub,&np->uuA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->up_sub,&np->upA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pu_sub,&np->puA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pp_sub,&np->ppA))
    NP_RETURN(1,result[0]);

  if (np->u_solver != NULL)
    if (np->u_solver->PreProcess != NULL)
      if ((*np->u_solver->PreProcess)
            (np->u_solver,level,np->ux,np->ub,np->uuA,
            baselevel,result))
        REP_ERR_RETURN(1);

  if (np->p_solver != NULL)
    if (np->p_solver->PreProcess != NULL)
      if ((*np->p_solver->PreProcess)
            (np->p_solver,level,np->px,np->pb,np->ppA,
            baselevel,result))
        REP_ERR_RETURN(1);

  return (0);
}

static INT APSmoother (NP_ITER *theNP, INT level,
                       VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                       INT *result)
{
  NP_TS *np = (NP_TS *) theNP;
  MULTIGRID *theMG = NP_MG(theNP);
  LRESULT lresult;

  /* get storage for extra temp */
  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->u_sub,&np->ub))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->p_sub,&np->pb))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,0,level,np->ux,&np->t))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,0,level,np->px,&np->s))
    NP_RETURN(1,result[0]);
  if (dcopy(theMG,0,level,ALL_VECTORS,np->t,np->ub) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dcopy(theMG,0,level,ALL_VECTORS,np->s,np->pb) != NUM_OK)
    NP_RETURN(1,result[0]);

    #ifdef ModelP
  if (a_vector_collect(theMG,0,level,np->t)!=NUM_OK)
    NP_RETURN(1,result[0]);
    #endif
  if (dset(theMG,0,level,ALL_VECTORS,x,0.0)!= NUM_OK)
    NP_RETURN(1,result[0]);
  if ((*np->u_solver->Residuum)
        (np->u_solver,0,level,np->ux,np->t,np->uuA,&lresult))
    NP_RETURN(1,result[0]);
  if ((*np->u_solver->Solver)(np->u_solver,level,
                              np->ux,np->t,np->uuA,
                              np->u_solver->abslimit,
                              np->u_solver->reduction,&lresult))
    NP_RETURN(1,result[0]);

  /*
     UserWriteF("np->ux \n");
     PrintSVector(theMG,np->ux);
   */

  if (dmatmul_minus(theMG,0,level,ALL_VECTORS,np->s,np->puA,np->ux)
      != NUM_OK)
    NP_RETURN(1,result[0]);
    #ifdef ModelP
  if (a_vector_collect_noskip(theMG,0,level,np->s)!=NUM_OK)
    NP_RETURN(1,result[0]);
    #endif

  /*
     UserWriteF("np->s \n");
     PrintSVector(theMG,np->s);
   */

  if ((*np->p_solver->Residuum)
        (np->p_solver,0,level,np->px,np->s,np->ppA,&lresult))
    NP_RETURN(1,result[0]);
  if ((*np->p_solver->Solver)(np->p_solver,level,
                              np->px,np->s,np->ppA,
                              np->p_solver->abslimit,
                              np->p_solver->reduction,&lresult))
    NP_RETURN(1,result[0]);
  /*
     UserWriteF("np->px \n");
     PrintSVector(theMG,np->px);
   */

  /* damp */
  if (dscalx(NP_MG(theNP),0,level,ALL_VECTORS,x,np->damp)!= NUM_OK)
    NP_RETURN(1,result[0]);
  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),0,level,ALL_VECTORS,b,A,x)!= NUM_OK)
    NP_RETURN(1,result[0]);

  FreeVD(NP_MG(theNP),0,level,np->t);
  FreeVD(NP_MG(theNP),0,level,np->s);

  return (0);
}

static INT APPostProcess (NP_ITER *theNP, INT level,
                          VECDATA_DESC *x, VECDATA_DESC *b,
                          MATDATA_DESC *A, INT *result)
{
  NP_TS *np = (NP_TS *) theNP;

  if (np->p_solver != NULL)
    if (np->p_solver->PostProcess != NULL)
      if ((*np->p_solver->PostProcess)
            (np->p_solver,level,np->px,np->pb,np->ppA,result))
        REP_ERR_RETURN(1);

  if (np->u_solver != NULL)
    if (np->u_solver->PostProcess != NULL)
      if ((*np->u_solver->PostProcess)
            (np->u_solver,level,np->ux,np->ub,np->uuA,result))
        REP_ERR_RETURN(1);

  return(0);
}

static INT APConstruct (NP_BASE *theNP)
{
  NP_ITER *np = (NP_ITER *) theNP;

  theNP->Init     = APInit;
  theNP->Display  = APDisplay;
  theNP->Execute  = NPIterExecute;
  np->PreProcess  = APPreProcess;
  np->Iter        = APSmoother;
  np->PostProcess = APPostProcess;

  return(0);
}


static INT BHRInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_TS *np;

  np = (NP_TS *) theNP;

  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "u",argc,argv,&(np->u_sub));
  if (np->vt == NULL) {
    UserWriteF("TSInit: no subtemplate u found\n");
    return(NP_NOT_ACTIVE);
  }
  np->vt = ReadArgvVecTemplateSub(MGFORMAT(theNP->mg),
                                  "p",argc,argv,&(np->p_sub));
  if (np->vt == NULL) {
    UserWriteF("TSInit: no subtemplate p found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "uu",argc,argv,&(np->uu_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate uu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "up",argc,argv,&(np->up_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate up found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pu",argc,argv,&(np->pu_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate pu found\n");
    return(NP_NOT_ACTIVE);
  }
  np->mt = ReadArgvMatTemplateSub(MGFORMAT(theNP->mg),
                                  "pp",argc,argv,&(np->pp_sub));
  if (np->mt == NULL) {
    UserWriteF("TSInit: no subtemplate pp found\n");
    return(NP_NOT_ACTIVE);
  }
  np->u_iter = (NP_ITER *)
               ReadArgvNumProc(theNP->mg,"UI",ITER_CLASS_NAME,argc,argv);
  if (np->u_iter == NULL) {
    UserWriteF("TSInit: no iter UI found\n");
    return(NP_NOT_ACTIVE);
  }

  np->t = NULL;
  return (NPIterInit(&np->iter,argc,argv));
}

static INT BHRDisplay (NP_BASE *theNP)
{
  NP_TS *np;

  np = (NP_TS *) theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (np->u_iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"UI",ENVITEM_NAME(np->u_iter));
  return (0);
}

static INT BHRPreProcess  (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_TS *np;
  GRID *theGrid;

  np = (NP_TS *) theNP;
  theGrid = NP_GRID(theNP,level);
  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->u_sub,&np->ub))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(b,np->vt,np->p_sub,&np->pb))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->uu_sub,&np->uuA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->up_sub,&np->upA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pu_sub,&np->puA))
    NP_RETURN(1,result[0]);
  if (MDsubDescFromMT(A,np->mt,np->pp_sub,&np->ppA))
    NP_RETURN(1,result[0]);
  *baselevel = level;

  if (np->u_iter->PreProcess != NULL)
    if ((*np->u_iter->PreProcess)
          (np->u_iter,level,np->ux,np->ub,np->puA,baselevel,result))
      REP_ERR_RETURN(1);

  return (0);
}

static INT BHRSmoother (NP_ITER *theNP, INT level,
                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                        INT *result)
{
  NP_TS *np;
  MULTIGRID *theMG;

  np = (NP_TS *) theNP;
  theMG = NP_MG(theNP);

  /* get storage for extra temp */
  if (VDsubDescFromVT(x,np->vt,np->u_sub,&np->ux))
    NP_RETURN(1,result[0]);
  if (VDsubDescFromVT(x,np->vt,np->p_sub,&np->px))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,level,level,np->ux,&np->t))
    NP_RETURN(1,result[0]);
  if (dcopy(theMG,level,level,ALL_VECTORS,np->t,np->ub) != NUM_OK)
    NP_RETURN(1,result[0]);
  if ((*np->u_iter->Iter)(np->u_iter,level,np->px,np->t,np->upA,result))
    REP_ERR_RETURN(1);

  if (dcopy(theMG,level,level,ALL_VECTORS,np->t,np->pb) != NUM_OK)
    NP_RETURN(1,result[0]);
  if ((*np->u_iter->Iter)(np->u_iter,level,np->ux,np->t,np->puA,result))
    REP_ERR_RETURN(1);

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->damp)!= NUM_OK)
    NP_RETURN(1,result[0]);

  /* update defect */
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK)
    NP_RETURN(1,result[0]);

  FreeVD(NP_MG(theNP),level,level,np->t);

  return (0);
}

static INT BHRPostProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *result)
{
  NP_TS *np;

  np = (NP_TS *) theNP;
  if (np->u_iter->PostProcess != NULL)
    if ((*np->u_iter->PostProcess)
          (np->u_iter,level,x,b,A,result))
      REP_ERR_RETURN(1);

  return(0);
}

static INT BHRConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = BHRInit;
  theNP->Display = BHRDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = BHRPreProcess;
  np->Iter = BHRSmoother;
  np->PostProcess = BHRPostProcess;

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

static INT GetSpec (DOUBLE *m, DOUBLE *l)
{
  DOUBLE s,t;

  s=0.5*ABS(m[0]+m[3]);
  t=s*s-m[0]*m[3]+m[1]*m[2];
  if (t>=0.0)
  {
    t=sqrt(t);
    l[0]=ABS(s-t);
    l[1]=ABS(s+t);
    if (l[0]>l[1])
    {
      s=l[0]; l[0]=l[1]; l[1]=s;
    }
  }
  else
  {
    l[0]=l[1]=sqrt(s*s+t*t);
  }
  return(0);
}

static INT SetAutoDamp_SSpec (GRID *g, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv)
{
  INT n,comp;
  SHORT *advcomp;
  VECTOR *v;
  MATRIX *m;
  DOUBLE diag,sum,spec[2],d;

  advcomp = VD_ncmp_cmpptr_of_otype(adv,NODEVEC,&n);
  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
  {
    comp = MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,0);
    GetSpec(&MVALUE(VSTART(v),comp),spec);
    diag=spec[0];
    if (diag==0.0) return(1);
    sum=0.0;
    for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
    {
      GetSpec(&MVALUE(m,comp),spec);
      sum+=spec[1];
    }
    if (diag>=sum)
    {
      VVALUE(v,advcomp[0])=damp[0];
      VVALUE(v,advcomp[1])=damp[1];
    }
    else
    {
      d=diag/sum;
      VVALUE(v,advcomp[0])=damp[0]*d;
      VVALUE(v,advcomp[1])=damp[1]*d;
    }
  }

  return(0);
}

static INT SetAutoDamp_WSpec (GRID *g, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv)
{
  INT n,comp;
  SHORT *advcomp;
  VECTOR *v;
  MATRIX *m;
  DOUBLE diag,sum,spec[2];

  advcomp = VD_ncmp_cmpptr_of_otype(adv,NODEVEC,&n);
  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
  {
    comp = MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,0);
    GetSpec(&MVALUE(VSTART(v),comp),spec);
    diag=spec[1];
    if (diag==0.0) return(1);
    sum=0.0;
    for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
    {
      GetSpec(&MVALUE(m,comp),spec);
      sum+=spec[1];
    }
    if (diag>=sum)
    {
      VVALUE(v,advcomp[0])=damp[0];
      VVALUE(v,advcomp[1])=damp[1];
    }
    else
    {
      VVALUE(v,advcomp[0])=damp[0]*diag/sum;
      VVALUE(v,advcomp[1])=damp[1]*diag/sum;
    }
  }

  return(0);
}

static INT SetAutoDamp_Scalar (GRID *g, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv)
{
  INT i,n,comp;
  SHORT *advcomp;
  VECTOR *v;
  MATRIX *m;
  DOUBLE diag,sum;

  advcomp = VD_ncmp_cmpptr_of_otype(adv,NODEVEC,&n);
  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
    for (i=0; i<n; i++)
    {
      comp = MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,i*(n+1));
      diag=ABS(MVALUE(VSTART(v),comp));
      if (diag==0.0) return(1);
      sum=0.0;
      for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
        sum+=ABS(MVALUE(m,comp));
      if (diag>=sum) VVALUE(v,advcomp[i])=damp[i];
      else VVALUE(v,advcomp[i])=damp[i]*diag/sum;
    }

  return(0);
}

static INT SetAutoDamp_Test (GRID *g, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv)
{
  INT n,comp;
  SHORT *advcomp;
  VECTOR *v;
  MATRIX *m;
  DOUBLE diag,sum;

  advcomp = VD_ncmp_cmpptr_of_otype(adv,NODEVEC,&n);
  assert(n==2);
  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
  {
    /* c comp */
    comp = MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,0);
    diag=ABS(MVALUE(VSTART(v),comp));
    if (diag==0.0) return(1);
    sum=0.0;
    for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
      sum+=ABS(MVALUE(m,comp));
    if (diag>=sum) VVALUE(v,advcomp[0])=damp[0];
    else VVALUE(v,advcomp[0])=damp[0]*diag/sum;

    /* p comp */
    VVALUE(v,advcomp[1])=damp[1];

    /* check determinant */
    sum=ABS(MVALUE(VSTART(v),comp))*ABS(MVALUE(VSTART(v),comp+3))+ABS(MVALUE(VSTART(v),comp+1))*ABS(MVALUE(VSTART(v),comp+2));
    sum=ABS(MVALUE(VSTART(v),comp)*MVALUE(VSTART(v),comp+3)-MVALUE(VSTART(v),comp+1)*MVALUE(VSTART(v),comp+2))/sum;
    if (sum<0.5)
    {
      /* UserWriteF("regularizing %e: ",sum);
         if (sum<VVALUE(v,advcomp[0])) UserWriteF("c ");
         if (sum<VVALUE(v,advcomp[1])) UserWriteF("p ");
         UserWriteF("\n"); */
      VVALUE(v,advcomp[0])=MIN(sum,VVALUE(v,advcomp[0]));
      VVALUE(v,advcomp[1])=MIN(sum,VVALUE(v,advcomp[1]));
    }
  }

  return(0);
}

static INT SetAutoDamp_HZ (GRID *g, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv)
{
  INT i,n,comp;
  SHORT *advcomp;
  VECTOR *v;
  MATRIX *m;
  DOUBLE diag,sum;

  advcomp = VD_ncmp_cmpptr_of_otype(adv,NODEVEC,&n);
  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
    for (i=0; i<n; i++)
    {
      comp = MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,i*(n+1));
      diag=MVALUE(VSTART(v),comp);
      if (diag==0.0) return(1);
      sum=0.0;
      for (m=VSTART(v); m!=NULL; m=MNEXT(m))
        sum+=MVALUE(m,comp);
      if (diag*sum>=0.0) VVALUE(v,advcomp[i])=damp[i];
      else VVALUE(v,advcomp[i])=damp[i]*ABS(diag)/(ABS(diag)+ABS(sum));
    }

  return(0);
}

static INT SetAutoDamp_EV (NP_ITER *theNP, GRID *g, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv, INT n)
{
  NP_SMOOTHER *np = (NP_SMOOTHER *) theNP;
  INT i,r,level,nc,autodamp;
  VECDATA_DESC *t,*p;
  SHORT *advc,*pc;
  VECTOR *v;

  level=GLEVEL(g);
  t=p=NULL;
  if (AllocVDFromVD(MYMG(g),level,level,adv,&t)) return(1);
  if (AllocVDFromVD(MYMG(g),level,level,adv,&p)) return(1);
  l_dsetrandom(g,adv,EVERY_CLASS,1.0);
  l_dset(g,adv,EVERY_CLASS,1.0);
  l_dset(g,t,EVERY_CLASS,0.0);
  autodamp=np->AutoDamp; np->AutoDamp=0;
  for (i=0; i<n; i++)
  {
    if (dcopy(MYMG(g),level,level,ALL_VECTORS,p,adv)!= NUM_OK) return(1);
    if (SSORSmoother(theNP,level,t,adv,A,&r)) return(1);
  }
  np->AutoDamp=autodamp;
  advc = VD_ncmp_cmpptr_of_otype(adv,NODEVEC,&nc);
  pc = VD_ncmp_cmpptr_of_otype(p,NODEVEC,&nc);
  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
    for (i=0; i<nc; i++)
    {
      if (ABS(VVALUE(v,advc[i]))>ABS(VVALUE(v,pc[i]))) VVALUE(v,advc[i])=sqrt(ABS(VVALUE(v,pc[i])/VVALUE(v,advc[i])));
      else VVALUE(v,advc[i])=1.0;
    }

  if (FreeVD(MYMG(g),level,level,t)) return(1);
  if (FreeVD(MYMG(g),level,level,p)) return(1);

  return(0);
}

static INT SetAutoDamp (NP_ITER *theNP, GRID *g, INT mode, MATDATA_DESC *A, const DOUBLE *damp, VECDATA_DESC *adv)
{
  INT r,n;

  switch (mode)
  {
  case 1 :
    r=SetAutoDamp_Scalar(g,A,damp,adv);
    break;
  case 2 :
    r=SetAutoDamp_SSpec(g,A,damp,adv);
    break;
  case 3 :
    r=SetAutoDamp_WSpec(g,A,damp,adv);
    break;
  case 4 :
    r=SetAutoDamp_HZ(g,A,damp,adv);
    break;
  case 5 :
    r=SetAutoDamp_Test(g,A,damp,adv);
    break;
  default :
    r=1;
    if (mode<0)
    {
      n=-mode;
      r=SetAutoDamp_EV(theNP,g,A,damp,adv,n);
    }
  }

  return(r);
}

static INT SORInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SMOOTHER *np;

  np = (NP_SMOOTHER *) theNP;
  np->AutoDamp = ReadArgvOption("autodmp",argc,argv);
  np->DampVector = ReadArgvVecDesc(NP_MG(np),"dv",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT SORDisplay (NP_BASE *theNP)
{
  NP_SMOOTHER *np = (NP_SMOOTHER *) theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"autodmp",np->AutoDamp);
  if (np->DampVector != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"dv",ENVITEM_NAME(np->DampVector));

  return (0);
}

static INT SORPreProcess  (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_SMOOTHER *np;
  GRID *theGrid;
  MATDATA_DESC *myMat;

  np = (NP_SMOOTHER *) theNP;
  theGrid = NP_GRID(theNP,level);
        #ifdef ModelP
  if (np->cons_mode == MAT_DIAG_VEC_CONS) {
    CopyDiagMatrix(theGrid,A,np->diag);
    if (l_vector_consistent(theGrid,np->diag) != NUM_OK)
      NP_RETURN(1,result[0]);
    np->L = A;
  }
  else {
    if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->L))
      NP_RETURN(1,result[0]);
    if( !np->Liscomplete )
    {
      if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,A) != NUM_OK)
        NP_RETURN(1,result[0]);
      if (l_matrix_consistent(theGrid,np->L,np->cons_mode) != NUM_OK)
        NP_RETURN(1,result[0]);
    }
    np->diag = NULL;
  }
  myMat = np->L;
        #else
  myMat = A;
        #endif

  if (np->Order!=NULL)
    if ((*np->Order->Order)(np->Order,level,myMat,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (np->AutoDamp)
  {
    if (AllocVDFromVD(NP_MG(theNP),level,level,x,&np->DampVector)) NP_RETURN(1,result[0]);
    if (SetAutoDamp(theNP,theGrid,np->AutoDamp,myMat,np->damp,np->DampVector)) NP_RETURN(1,result[0]);
    DPrintVector(MYMG(theGrid),np->DampVector);
  }
  *baselevel = level;

  return (0);
}


static INT SORStep (NP_SMOOTHER *theNP, INT level,
                    VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                    MATDATA_DESC *L,
                    INT *result)
{
  NP_SMOOTHER *np = (NP_SMOOTHER *) theNP;

    #ifdef ModelP
  GRID *theGrid = NP_GRID(theNP,level);
  if (np->cons_mode == MAT_MASTER_CONS) {
    if (l_vector_collect(theGrid,b)!=NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  if (np->AutoDamp)
  {
    if (l_lsor_ld(NP_GRID(theNP,level),x,L,b,np->DampVector,
                  np->diag) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_lsor(NP_GRID(theNP,level),x,L,b,theNP->damp,np->diag) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
    #else
  if (np->AutoDamp)
  {
    if (l_lsor_ld(NP_GRID(theNP,level),x,A,b,np->DampVector,NULL) != NUM_OK) NP_RETURN(1,result[0]);
  }
  else
  {
    if (l_lsor(NP_GRID(theNP,level),x,A,b,theNP->damp,NULL) != NUM_OK) NP_RETURN(1,result[0]);
  }
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
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)
      != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SORPostProcess (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SMOOTHER *np;

  np = (NP_SMOOTHER *) theNP;
  if (np->AutoDamp)
    if (FreeVD(np->iter.base.mg,level,level,np->DampVector)) REP_ERR_RETURN(1);

  return(SmootherPostProcess(theNP,level,x,b,A,result));
}

static INT SORConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SORInit;
  theNP->Display = SORDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = SORPreProcess;
  np->iter.Iter = SORSmoother;
  np->iter.PostProcess = SORPostProcess;
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
        if (ReadVecTypeINTs(MGFORMAT(NP_MG(theNP)),value,MAX_BLOCKS+1,nTypeBlocks,TypeBlocks)!=0)
          REP_ERR_RETURN (NP_NOT_ACTIVE)
          else {bopt = TRUE; continue;}

      /* BlockOrder */
      if (strstr(option,"BlockOrder")!=NULL)
        if (ReadVecTypeOrder(MGFORMAT(NP_MG(theNP)),value,MAX_ORDER,MAX_BLOCKS,&SBGS_NBLOCKITER(theSBGS),SBGS_BLOCKORDER(theSBGS))!=0)
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
    sprintf(name," block%d(%s)",i,ObjTypeName[SBGS_BLOCKDESC(theSBGS,i).tp]);
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

    if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,
                      SBGS_VD_ro(theSBGS,bl),
                      SBGS_MD_Ao(theSBGS,bl),
                      SBGS_VD_cd(theSBGS)))
      NP_RETURN (1,result[0]);
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
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->L)) NP_RETURN(1,result[0]);
  if( !np->Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,A) != NUM_OK) NP_RETURN(1,result[0]);
                #ifdef ModelP
    if (l_matrix_consistent(theGrid,np->L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
                #endif
  }
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
  sc_read(np->beta,NP_FMT(np),np->smoother.iter.b,"beta",argc,argv);
  for (i=0; i<MAX_VEC_COMP; i++) np->mindiag[i] = 0.0;
  sc_read(np->mindiag,NP_FMT(np),np->smoother.iter.b,"mindiag",argc,argv);

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
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L))
    NP_RETURN(1,result[0]);
  if( !np->smoother.Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK)
      NP_RETURN(1,result[0]);
                #ifdef ModelP
    if (l_matrix_consistent(theGrid,np->smoother.L,np->smoother.cons_mode)
        != NUM_OK)
      NP_RETURN(1,result[0]);
                #endif
  }
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
  NP_ILU *np;
  GRID *theGrid;

  np = (NP_ILU *) theNP;
  theGrid = NP_GRID(theNP,level);
    #ifdef ModelP
  if (np->smoother.cons_mode == MAT_MASTER_CONS) {
    if (l_vector_collect(theGrid,b)!=NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
    #endif

  /*
     PrintVector(theGrid,x,3,3);
     PrintMatrix(theGrid,L,3,3);
     PrintVector(theGrid,b,3,3);
   */

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

static INT BDInit (NP_BASE *theNP, INT argc , char **argv)
{
  return (SmootherInit(theNP,argc,argv));
}

static INT BDDisplay (NP_BASE *theNP)
{
  SmootherDisplay(theNP);

  return (0);
}

INT l_bdpreprocess1 (GRID *g, MATDATA_DESC *A, MATDATA_DESC *L)
{
  VECTOR *v;

  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
  {
    MATRIX *m = VSTART(v);
    if (VECSKIP(v) == 0)
      MVALUE(m,MD_SCALCMP(L)) = 1.0 / MVALUE(m,MD_SCALCMP(L));
    else
      MVALUE(m,MD_SCALCMP(L)) = 0.0;
    for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
      MVALUE(m,MD_SCALCMP(L)) = 0.0;
  }

  return (0);
}

INT l_bdpreprocess2 (GRID *g, VECDATA_DESC *x,
                     MATDATA_DESC *A, MATDATA_DESC *L)
{
  ELEMENT *e;
  VECTOR *v;

  dmatset(MYMG(g),LEVEL(g),LEVEL(g),ALL_VECTORS,L,0.0);

  for (e=FIRSTELEMENT(g); e!=NULL; e=SUCCE(e)) {
    VECTOR *v[MAX_NODAL_VECTORS];
    INT cnt = GetAllVectorsOfElementOfType(e,v,x);
    DOUBLE mat[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
    DOUBLE imat[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
    DOUBLE rmat[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
    INT m = GetVlistMValues(cnt,v,A,mat);
    INT i;

    /*
       for (i=0; i<m; i++)
        for (j=0; j<m; j++)
              if (i!= j)
                mat[i*m+j] = 0.0;
     */
    /*
       for (i=0; i<m; i++)
        for (j=0; j<m; j++)
                rmat[i*m+j] = 0.0;

       for (j=0; j<m; j++)
        rmat[j*m+j] = 1.0;
     */

    if (InvertFullMatrix_piv(m,mat,imat))
      return(1);

    GetVlistMValues(cnt,v,L,rmat);

    for (i=0; i<m*m; i++)
      imat[i] -= rmat[i];

    AddVlistMValues(g,cnt,v,L,imat);
  }
  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    INT type = VTYPE(v);
    INT ncomp = VD_NCMPS_IN_TYPE (x,type);
    INT i,j;

    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(v) & (1<<j))
      {
        MATRIX *m = VSTART(v);

        for (i=j*ncomp; i<(j+1)*ncomp; i++)
          MVALUE(m,MD_MCMP_OF_RT_CT(L,type,type,i)) = 0.0;
        for (m=MNEXT(m); m!=NULL; m=MNEXT(m))
        {
          INT dtype = MDESTTYPE(m);
          INT dcomp = VD_NCMPS_IN_TYPE (x,dtype);

          if (dcomp == 0) continue;
          for (i=j*dcomp; i<(j+1)*dcomp; i++)
            MVALUE(m,MD_MCMP_OF_RT_CT(L,type,dtype,i)) = 0.0;
        }
      }
  }

  return (0);
}

INT l_bdpreprocess (GRID *g, VECDATA_DESC *x,
                    MATDATA_DESC *A, MATDATA_DESC *L)
{
  ELEMENT *e;
  VECTOR *v;

  dmatset(MYMG(g),GLEVEL(g),GLEVEL(g),ALL_VECTORS,L,0.0);

  for (e=FIRSTELEMENT(g); e!=NULL; e=SUCCE(e)) {
    VECTOR *v[MAX_NODAL_VECTORS];
    INT cnt = GetAllVectorsOfElementOfType(e,v,x);
    DOUBLE mat[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
    DOUBLE imat[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
    DOUBLE rmat[MAX_NODAL_VALUES*MAX_NODAL_VALUES];
    INT m = GetVlistMValues(cnt,v,A,mat);
    INT i,j,k,ci,cj;

    for (i=0; i<m*m; i++)
      rmat[i] = 0.0;
    for (j=0; j<m; j++)
      rmat[j*m+j] = 1.0;

    GetVlistMValues(cnt,v,L,imat);

    ci = 0;
    for (i=0; i<cnt; i++)
    {
      INT vitype = VTYPE(v[i]);
      INT vincomp = VD_NCMPS_IN_TYPE (x,vitype);

      cj = 0;
      for (j=0; j<cnt; j++)
      {
        MATRIX *m1,*m2;
        DOUBLE s[MAX_SINGLE_MAT_COMP];
        INT vjtype = VTYPE(v[j]);
        INT vjncomp = VD_NCMPS_IN_TYPE (x,vjtype);
        INT i1,i2,i3;

        for (i1=0; i1<vincomp; i1++)
          for (i2=0; i2<vjncomp; i2++)
            s[i1*vjncomp+i2] = 0;
        for (m1=VSTART(v[i]); m1!=NULL; m1=MNEXT(m1))
          for (m2=VSTART(v[j]); m2!=NULL; m2=MNEXT(m2)) {
            INT vktype;
            INT vkncomp = VD_NCMPS_IN_TYPE (x,vktype);

            if (MDEST(m1) != MDEST(m2)) continue;
            vktype = VTYPE(MDEST(m1));
            for (i1=0; i1<vincomp; i1++)
              for (i2=0; i2<vjncomp; i2++)
              {
                DOUBLE t = 0.0;

                for (i3=0; i3<vkncomp; i3++)
                  t +=
                    MVALUE(m1,MD_MCMP_OF_RT_CT(L,vitype,vktype,i1*vkncomp+i3)) *
                    MVALUE(MADJ(m2),MD_MCMP_OF_RT_CT(A,vktype,vjtype,i3*vjncomp+i2));
                s[i1*vjncomp+i2] = t;
              }
          }


        for (i1=0; i1<vincomp; i1++)
          for (i2=0; i2<vjncomp; i2++)
            rmat[(ci+i1)*m+cj+i2] -= s[i1*vjncomp+i2];
        cj += vjncomp;
      }
      ci += vincomp;
    }

    /*
       printf("mat\n");
       for (i=0; i<m; i++) {
        for (j=0; j<m; j++)
                printf("%6.1e",mat[i*m+j]);
            printf("\n");
       }
     */

    if (InvertFullMatrix_piv(m,mat,imat))
      return(1);

    /*
       printf("imat\n");
       for (i=0; i<m; i++) {
        for (j=0; j<m; j++)
                printf("%6.1e",imat[i*m+j]);
            printf("\n");
       }
       printf("rmat\n");
       for (i=0; i<m; i++) {
        for (j=0; j<m; j++)
                printf("%6.1e",rmat[i*m+j]);
            printf("\n");
       }
     */
    for (i=0; i<m; i++)
      for (j=0; j<m; j++)
      {
        DOUBLE s = 0.0;

        for (k=0; k<m; k++)
          s += rmat[i*m+k] * imat[k*m+j];
        mat[i*m+j] = s;
      }
    AddVlistMValues(g,cnt,v,L,mat);
  }
  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    INT type = VTYPE(v);
    INT ncomp = VD_NCMPS_IN_TYPE (x,type);
    INT i,j;

    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(v) & (1<<j))
      {
        MATRIX *m = VSTART(v);

        for (i=j*ncomp; i<(j+1)*ncomp; i++)
          MVALUE(m,MD_MCMP_OF_RT_CT(L,type,type,i)) = 0.0;
        for (m=MNEXT(m); m!=NULL; m=MNEXT(m))
        {
          INT dtype = MDESTTYPE(m);
          INT dcomp = VD_NCMPS_IN_TYPE (x,dtype);

          if (dcomp == 0) continue;
          for (i=j*dcomp; i<(j+1)*dcomp; i++)
            MVALUE(m,MD_MCMP_OF_RT_CT(L,type,dtype,i)) = 0.0;
        }
      }
  }

  return (0);
}

static INT BDPreProcess (NP_ITER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *baselevel, INT *result)
{
  NP_SMOOTHER *np = (NP_SMOOTHER *)theNP;
  GRID *theGrid = NP_GRID(theNP,level);

  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->L))
    NP_RETURN(1,result[0]);
  if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->L,A) != NUM_OK)
    NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->L,np->cons_mode)
      != NUM_OK)
    NP_RETURN(1,result[0]);
        #endif
  *baselevel = level;

  l_bdpreprocess(theGrid,x,A,np->L);

  return (0);
}

static INT BDStep (NP_SMOOTHER *theNP, INT level,
                   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                   MATDATA_DESC *L,
                   INT *result)
{
    #ifdef ModelP
  GRID *theGrid = NP_GRID(theNP,level);

  if (theNP->cons_mode == MAT_MASTER_CONS) {
    if (l_vector_collect(theGrid,b)!=NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
    #endif

  /* PrintVector(theGrid,b,3,3); */

  /*dmatset(MYMG(theGrid),level,level,ALL_VECTORS,L,0.0);*/

  if (dmatmul(NP_MG(theNP),level,level,ON_SURFACE,x,L,b) != NUM_OK)
    NP_RETURN(1,result[0]);

  return (0);
}

static INT BDConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = BDInit;
  theNP->Display = BDDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = BDPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = BDStep;

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
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if( !np->smoother.Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
                #ifdef ModelP
    if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
                #endif
  }
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
  sc_read(np->beta,NP_FMT(np),np->smoother.iter.b,"beta",argc,argv);
  sc_read(np->thresh,NP_FMT(np),np->smoother.iter.b,"thresh",argc,argv);

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
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if( !np->smoother.Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
                #ifdef ModelP
    if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
                #endif
  }
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
  sc_read(np->beta,NP_FMT(np),np->smoother.iter.b,"beta",argc,argv);

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
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(NP_MG(theNP),level,level,x,&tmp)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if( !np->smoother.Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
                #ifdef ModelP
    if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
                #endif
  }
  if (l_iluspdecomp(theGrid,np->smoother.L,np->beta,tmp,np->mode,NULL)
      !=NUM_OK) {
    PrintErrorMessage('E',"SPILUPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
  }
  *baselevel = level;

  if (FreeVD(NP_MG(theNP),level,level,tmp)) REP_ERR_RETURN(1);

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
   spblilu - numproc for point block ilu smoother for sparse block matrices

   DESCRIPTION:
   This numproc executes a point block ilu smoother, using the blas routines
   'l_iluspbldecomp' and 'l_iluspbliter'. It can be used in 'lmgc'.
   The Beta-correction and the threshhold mode are not implemented (yet).

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
   .n  ($beta is only used if the matrix blocks are not sparse. Else it
   .n   is ignored.)

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT SPBLILUInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SPBLILU *np;
  INT i;

  np = (NP_SPBLILU *) theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->beta[i] = 0.0;
  sc_read(np->beta,NP_FMT(np),np->smoother.iter.b,"beta",argc,argv);
  for (i=0; i<MAX_VEC_COMP; i++) np->mindiag[i] = 0.0;
  sc_read(np->mindiag,NP_FMT(np),np->smoother.iter.b,"mindiag",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT SPBLILUDisplay (NP_BASE *theNP)
{
  NP_SPBLILU *np;

  SmootherDisplay(theNP);
  np = (NP_SPBLILU *) theNP;
  if (sc_disp(np->beta,np->smoother.iter.b,"beta"))
    REP_ERR_RETURN (1);
  if (sc_disp(np->mindiag,np->smoother.iter.b,"mindiag"))
    REP_ERR_RETURN (1);

  return (0);
}

static INT SPBLILUPreProcess (NP_ITER *theNP, INT level,
                              VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                              INT *baselevel, INT *result)
{
  NP_SPBLILU *np;
  GRID *theGrid;

  np = (NP_SPBLILU *) theNP;
  theGrid = NP_GRID(theNP,level);
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  /* Natalia - np->smoother.L = A; */
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L))
    NP_RETURN(1,result[0]);
  if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK)
    NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (l_matrix_consistent(theGrid,np->smoother.L,np->smoother.cons_mode)
      != NUM_OK)
    NP_RETURN(1,result[0]);
        #endif
  if (l_iluspbldecomp(theGrid,np->smoother.L,np->beta)
      !=NUM_OK) {
    PrintErrorMessage('E',"SPBLILUPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
  }
  *baselevel = level;

  return (0);
}

static INT SPBLILUStep (NP_SMOOTHER *theNP, INT level,
                        VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                        MATDATA_DESC *L,
                        INT *result)
{
  NP_SPBLILU *np;
  GRID *theGrid;

  np = (NP_SPBLILU *) theNP;
  theGrid = NP_GRID(theNP,level);
    #ifdef ModelP
  if (np->smoother.cons_mode == MAT_MASTER_CONS) {
    if (l_vector_collect(theGrid,b)!=NUM_OK)
      NP_RETURN(1,result[0]);
  }
  else {
    if (l_vector_meanvalue(theGrid,b) != NUM_OK)
      NP_RETURN(1,result[0]);
  }
    #endif
  if (l_iluspbliter(NP_GRID(theNP,level),x,L,b) != NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SPBLILUConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = SPBLILUInit;
  theNP->Display = SPBLILUDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_SMOOTHER *) theNP;
  np->iter.PreProcess = SPBLILUPreProcess;
  np->iter.Iter = Smoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = SPBLILUStep;

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
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if( !np->smoother.Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
                #ifdef ModelP
    if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS)!=NUM_OK) NP_RETURN(1,result[0]);
                #endif
  }
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
       [$damp <sc double list>]
       [$regularize always|never|ifsing];
   .ve

   .  $c~<cor>							- correction vector
   .  $b~<rhs>							- right hand side vector
   .  $A~<mat>							- stiffness matrix
   .  $L~<mat>							- decomposition matrix
   .  $damp~<sc~double~list>			- damping factors for each component
   .  $regularize~always|never|ifsing	- regularize last block (default: ifsing)

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT LUInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_LU *np;
  char reg[32];

  np = (NP_LU *) theNP;
  if (ReadArgvChar("regularize",reg,argc,argv))
    np->regularize = REG_IF_SING;
  else if (strcmp(reg,LU_reg[REG_ALWAYS])==0)
    np->regularize = REG_ALWAYS;
  else if (strcmp(reg,LU_reg[REG_NEVER])==0)
    np->regularize = REG_NEVER;
  else if (strcmp(reg,LU_reg[REG_IF_SING])==0)
    np->regularize = REG_IF_SING;
  else
    REP_ERR_RETURN (NP_NOT_ACTIVE);

  return (SmootherInit(theNP,argc,argv));
}

static INT LUDisplay (NP_BASE *theNP)
{
  NP_LU *np;

  SmootherDisplay(theNP);
  np = (NP_LU *) theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SS,"regularize",LU_reg[np->regularize]);

  return (0);
}

static INT LUPreProcess (NP_ITER *theNP, INT level,
                         VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                         INT *baselevel, INT *result)
{
  NP_LU *np;
  GRID *theGrid;
  INT err;

  np = (NP_LU *) theNP;

  theGrid = NP_GRID(theNP,level);
  if (np->smoother.Order!=NULL)
    if ((*np->smoother.Order->Order)(np->smoother.Order,level,A,result)) NP_RETURN(1,result[0]);
  if (l_setindex(theGrid)) NP_RETURN(1,result[0]);
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L)) NP_RETURN(1,result[0]);
  if( !np->smoother.Liscomplete )
  {
    if (dmatcopy(NP_MG(theNP),level,level,ALL_VECTORS,np->smoother.L,A) != NUM_OK) NP_RETURN(1,result[0]);
                #ifdef ModelP
    if (l_matrix_consistent(theGrid,np->smoother.L,MAT_MASTER_CONS) != NUM_OK) NP_RETURN(1,result[0]);
                #endif
  }
  err = l_lrdecomp(theGrid,np->smoother.L);
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
    if ((err!=-VINDEX(LASTVECTOR(theGrid))) || (np->regularize==REG_NEVER))
    {
      PrintErrorMessageF('E',"LUPreProcess","decomp failed: IDX %ld on level %d",
                         -err,GLEVEL(theGrid));
      UserWriteF(" - LASTVECTOR has IDX %ld\n",
                 VINDEX(LASTVECTOR(theGrid)));
      NP_RETURN(1,result[0]);
    }
    if (l_lrregularize(theGrid,np->smoother.L,NO) !=NUM_OK) {
      PrintErrorMessage('E',"LUPreProcess","cannot regularize");
      NP_RETURN(1,result[0]);
    }
  }
  if (np->regularize==REG_ALWAYS)
    if (l_lrregularize(theGrid,np->smoother.L,YES) !=NUM_OK) {
      PrintErrorMessage('E',"LUPreProcess","cannot regularize");
      NP_RETURN(1,result[0]);
    }
  *baselevel = level;

  /*PrintMatrix(theGrid,np->smoother.L,3,3);*/

  return (0);
}

static INT LUConstruct (NP_BASE *theNP)
{
  NP_SMOOTHER *np;

  theNP->Init = LUInit;
  theNP->Display = LUDisplay;
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
   .  $tv2~<2.~testvector~sym> - symbol for the second testvector if necessary
   .  $t~<update~for~correction~sym> - temp. vector
   .  $type~<type of frequency filter> - "TFF" for Wagners or "FF" for Wittums
   .  $display - display mode: 'no', 'red'uced or 'full'
   .  $wr - relative frequency [0..1] for 2D OR 'all' for the whole logarithmic sequence of frequencies
   .  $wr3D - relative frequency [0..1] for 3D
   .  $AssDirichlet - assemble Dirichlet boundary conditions
   .  $parsim~[0|1] - perform simulation of the parallel algorithm on a sequential machine
   .  $SymmFrq~[0|1] - symmetric sequence of testfrequencies 1,2,4,...,n,...,4,2,1
   .  $CheckSymm~[0|1] - check if the preconditioner is symmetric

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
  if ( ReadArgvDOUBLE ( "wr3D", &NPFF_WaveNrRel3D(np), argc, argv) )
  {
    PrintErrorMessage('E',"FFInit", "Option $wr3D mandatory");
    REP_ERR_RETURN(1);
  }
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

  NPFF_ParSim(np) = FALSE;
  if ( ReadArgvINT ( "parsim", &NPFF_ParSim(np), argc, argv) )
  {
    NPFF_ParSim(np) = FALSE;
  }
  else
    NPFF_ParSim(np) = (NPFF_ParSim(np)==1);

  NPFF_AssDirichlet(np) = ReadArgvOption("AssDirichlet",argc,argv);
  NPFF_SymmFrq(np) = ReadArgvOption("SymmFrq",argc,argv);
  NPFF_CheckSymm(np) = ReadArgvOption("CheckSymm",argc,argv);


#ifdef __TWODIM__

#ifdef ModelP
  (void)InitBVDF( NPFF_BVDF(np), 256 );
#else
  if ( NPFF_ParSim(np) )
    (void)InitBVDF( NPFF_BVDF(np), 64 );
  else
    *NPFF_BVDF(np) = two_level_bvdf;
#endif

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

  UserWriteF(DISPLAY_NP_FORMAT_SI,"ParSim",(int)NPFF_ParSim(np));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"AssDirichlet",(int)NPFF_AssDirichlet(np));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"SymmFrq",(int)NPFF_SymmFrq(np));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"CheckSymm",(int)NPFF_CheckSymm(np));

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
   allocate temporarily the necessary data descriptors,
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

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(NP_MG(theNP),level,level,x,&NPFF_tv(np)))
    NP_RETURN(1,result[0]);

  if (NPFF_DO_FF(np))
  {
    if (AllocVDFromVD(NP_MG(theNP),level,level,x,&NPFF_tv2(np)))
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

#ifdef ModelP
  n = 3;                /* 2 for BPS, 3 for general */
#else
  n = 1;
#endif

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
    if (AllocMDFromMD(NP_MG(theNP),level,level,A,FF_MATDATA_DESC_ARRAY+i))
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

#ifdef __TWODIM__
  n+=10; /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! temp fuer CHECK_CALCS */
#endif

  for( i=0; i<n; i++ )
  {
    if (AllocVDFromVD(NP_MG(theNP),level,level,x,FF_VECDATA_DESC_ARRAY+i))
      NP_RETURN(1,result[0]);
    FF_Vecs[i] = VD_SCALCMP( FF_VECDATA_DESC_ARRAY[i] );
  }
  /*	}*/

  /* construction of Dirichlet boundary conditions only necessary if no
     fetransfer is done */
  if ( NPFF_AssDirichlet(np) )
  {
                #ifdef ModelP
    if (a_vector_vecskip(MYMG(theGrid), GLEVEL(theGrid), GLEVEL(theGrid), x)!= NUM_OK)
      NP_RETURN(1,result[0]);
                #endif

    if (AssembleDirichletBoundary (theGrid,A,x,b))
      NP_RETURN(1,result[0]);
    UserWrite(" [d]\n");
  }

#ifdef FF_ModelP
  /*    if (l_matrix_consistent(theGrid,np->smoother.L,np->smoother.cons_mode)
                  != NUM_OK)
              NP_RETURN(1,result[0]);
   */

  if ( PFFPreProcessIter( theGrid, &meshwidth, MD_SCALCMP( A ), VD_SCALCMP( x ), VD_SCALCMP( b ), NPFF_BVDF(np) )!=NUM_OK)
  {
    PrintErrorMessage('E',"FFPreProcess","preparation of the grid failed for ParSim");
    NP_RETURN(1,result[0]);
  }
#elif defined FF_PARALLEL_SIMULATION
  if ( NPFF_ParSim(np) )
  {
    if ( PFFPreProcessIter( theGrid, &meshwidth, MD_SCALCMP( A ), VD_SCALCMP( x ), VD_SCALCMP( b ), NPFF_BVDF(np) )!=NUM_OK)
    {
      PrintErrorMessage('E',"FFPreProcess","preparation of the grid failed for ParSim");
      NP_RETURN(1,result[0]);
    }
  }
  else
  {
    if ( FF_PrepareGrid( theGrid, &meshwidth, TRUE, MD_SCALCMP( A ), VD_SCALCMP( x ), VD_SCALCMP( b ), NPFF_BVDF(np) )!=NUM_OK)
    {
      PrintErrorMessage('E',"FFPreProcess","preparation of the grid failed");
      NP_RETURN(1,result[0]);
    }
  }
#else
  if (FF_PrepareGrid( theGrid, &meshwidth, TRUE, MD_SCALCMP( A ), VD_SCALCMP( x ), VD_SCALCMP( b ), NPFF_BVDF(np) )!=NUM_OK)
  {
    PrintErrorMessage('E',"FFPreProcess","preparation of the grid failed");
    NP_RETURN(1,result[0]);
  }
#endif /* FF_PARALLEL_SIMULATION */

  NPFF_MESHWIDTH(np) = meshwidth;

  BVD_INIT( &bvd );
  BVD_PUSH_ENTRY( &bvd, BVNUMBER(GFIRSTBV(theGrid)), NPFF_BVDF(np) );

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
    else
    if (NPFF_DO_FF(np) && FFDecomp( wavenr, wavenr3D, GFIRSTBV(theGrid), &bvd,
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

#ifdef FF_ModelP
  if ( PFFPostProcessIter(GRID_ON_LEVEL(theMG,level))!=NUM_OK)
  {
    PrintErrorMessage('E',"FFPostProcess","postprocessing of the grid failed in PFFPostProcessIter");
    NP_RETURN(1,result[0]);
  }
#endif

  if (NPFF_tv(np) != NULL)
    if (FreeVD(theMG,level,level,NPFF_tv(np))) REP_ERR_RETURN(1);
  if (NPFF_tv2(np) != NULL)
    if (FreeVD(theMG,level,level,NPFF_tv2(np))) REP_ERR_RETURN(1);

  for( i=1; i<FF_MAX_MATS; i++ )
  {
    if ( FF_MATDATA_DESC_ARRAY[i] != NULL )
    {
      if (FreeMD(theMG,level,level,FF_MATDATA_DESC_ARRAY[i])) REP_ERR_RETURN(1);
      FF_MATDATA_DESC_ARRAY[i] = NULL;
    }
    FF_Mats[i] = DUMMY_COMP;
  }

  for( i=0; i<FF_MAX_VECS; i++ )
  {
    if ( FF_Vecs[i] != DUMMY_COMP )
    {
      if (FreeVD(theMG,level,level,FF_VECDATA_DESC_ARRAY[i])) REP_ERR_RETURN(1);
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

static INT FFApplyPreconditioner( NP_FF *np, INT level,VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result, BV_DESC *bvd, GRID *theGrid)
{
#ifdef __BLOCK_VECTOR_DESC__
  DOUBLE end_wave, wavenr, start_norm, new_norm;
  INT ascending_frq;

  if ( !NPFF_ALLFREQ(np) )
  {             /* smooth only for 1 testvector frequency */
                /* copy defect to tv because FFMultWithMInv destroys its defect */
    dcopyBS( GFIRSTBV(theGrid), VD_SCALCMP( NPFF_tv(np) ), VD_SCALCMP( b ) );
#ifdef ModelP
    if (FFMultWithMInv( GFIRSTBV(theGrid), bvd, NPFF_BVDF(np),
                        VD_SCALCMP( x ),
                        VD_SCALCMP( NPFF_tv(np)),
                        x,
                        theGrid) != NUM_OK)
    {
      PrintErrorMessage('E',"FFStep","inversion failed");
      NP_RETURN(1,result[0]);
    }
#else
    if (FFMultWithMInv( GFIRSTBV(theGrid), bvd, NPFF_BVDF(np),
                        VD_SCALCMP( x ),
                        VD_SCALCMP( NPFF_tv(np) ) ) != NUM_OK)
    {
      PrintErrorMessage('E',"FFStep","inversion failed");
      NP_RETURN(1,result[0]);
    }
#endif
    /* defect -= A * corr_update */
    dmatmul_minusBS( GFIRSTBV(theGrid), bvd, NPFF_BVDF(np),
                     VD_SCALCMP( b ), MD_SCALCMP( A ), VD_SCALCMP( x ));
  }
  else
  {             /* smooth for all testvector frequencies */

    /* alloc temp. for correction update (in x!) */
    if (AllocVDFromVD(NP_MG((NP_ITER*)np),level,level,x,&NPFF_t(np))) NP_RETURN(1,result[0]);

    if ( NPFF_DISPLAY(np) != PCR_NO_DISPLAY )
      if(dnrm2BS( GFIRSTBV(theGrid), VD_SCALCMP( b ), &new_norm ) ) NP_RETURN(1,result[0]);

    ascending_frq = 1;
    end_wave = 1.0 / NPFF_MESHWIDTH(np) - 0.5;             /* rounding */
    /*for ( wavenr = 1.0; wavenr < end_wave; wavenr *= 2.0 )*/
    wavenr = 1.0;
    while (1)
    {                   /* wave 1.0 ... (1/h)/2 */
      if (NPFF_DO_TFF(np) )
      {
        if ( TFFDecomp( wavenr, wavenr, GFIRSTBV(theGrid), bvd,
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
        /*printf("wavenr %g\n", wavenr);*/
        if (FFDecomp( wavenr, wavenr, GFIRSTBV(theGrid), bvd,
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
#ifdef ModelP
      if (FFMultWithMInv( GFIRSTBV(theGrid), bvd, NPFF_BVDF(np),
                          VD_SCALCMP( NPFF_t(np) ),
                          VD_SCALCMP( NPFF_t(np) ),
                          NPFF_t(np),
                          theGrid ) != NUM_OK)
      {
        PrintErrorMessage('E',"FFStep","inversion failed");
        NP_RETURN(1,result[0]);
      }

      /* NOTE: the corr update is already consistent from FFMultWithMInv!
               don't try to make it consistent again! */
#else
      if (FFMultWithMInv( GFIRSTBV(theGrid), bvd, NPFF_BVDF(np),
                          VD_SCALCMP( NPFF_t(np) ),
                          VD_SCALCMP( NPFF_t(np) ) ) != NUM_OK)
      {
        PrintErrorMessage('E',"FFStep","inversion failed");
        NP_RETURN(1,result[0]);
      }
#endif
      /* corr += corr_update */
      daddBS( GFIRSTBV(theGrid), VD_SCALCMP( x ), VD_SCALCMP( NPFF_t(np) ) );

      /* defect -= A * corr_update */
      dmatmul_minusBS( GFIRSTBV(theGrid), bvd, NPFF_BVDF(np),
                       VD_SCALCMP( b ), MD_SCALCMP( A ), VD_SCALCMP( NPFF_t(np) ));

      if ( NPFF_DISPLAY(np) != PCR_NO_DISPLAY )
      {
        start_norm = new_norm;
        if(dnrm2BS( GFIRSTBV(theGrid), VD_SCALCMP( b ), &new_norm ) ) NP_RETURN(1,result[0]);

#ifdef ModelP
        UserWrite( "ONLY LOCAL:" );                             /* otherwise communication is necessary for the residuum */
#endif
        if( fabs(start_norm)<1e-20)
          printf(PFMT " start_norm == 0\n",me);
        else
          UserWriteF( "Wnr plane = %4g Wnr line = %4g new defect = %12g "
                      "conv. rate = %12g\n", wavenr, wavenr, new_norm,
                      new_norm/start_norm );
      }

      if ( ascending_frq )
      {
        wavenr *= 2.0;
        if( wavenr >= end_wave )
        {
          if( !NPFF_SymmFrq(np) )
            break;
          wavenr /= 4.0;
          ascending_frq = 0;
        }
      }
      else
      {
        wavenr /= 2.0;
        if ( wavenr < 0.999 )
          break;
      }
    }


    if (FreeVD(NP_MG((NP_ITER*)np),level,level,NPFF_t(np))!=NUM_OK) REP_ERR_RETURN(1);
  }
#else
  PrintErrorMessage( 'E', "FFApplyPreconditioner", "__BLOCK_VECTOR_DESC__ must be defined in gm.h" );
  REP_ERR_RETURN (1);
#endif

  return (NUM_OK);
}

static void FFGenerateCheckA( BLOCKVECTOR *bv, INT v_comp, INT nr_of_call, INT scaling_comp )
/* the result is inconsistent */
{
  register VECTOR *v;
  DOUBLE pos[DIM];

  BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), BVENDVECTOR( bv ) )
  {
    VectorPosition(v,pos);
    VVALUE( v, v_comp ) = exp(pos[0])*(1.0-pos[1]);
    if( nr_of_call > 1 )
      VVALUE( v, v_comp ) *= VVALUE( v, scaling_comp );
  }
}

static void FFGenerateCheckB( BLOCKVECTOR *bv, INT v_comp, INT nr_of_call, INT scaling_comp )
/* the result is inconsistent */
{
  register VECTOR *v;
  DOUBLE pos[DIM];

  BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), BVENDVECTOR( bv ) )
  {
    VectorPosition(v,pos);
    VVALUE( v, v_comp ) = sin(13.423*pos[0])*exp(1.0-pos[1]);
    if( nr_of_call > 1 )
      VVALUE( v, v_comp ) *= VVALUE( v, scaling_comp );
  }
}

void FFCopyVector( GRID *grid, INT dest, INT source )
/* copy the whole vector in the grid; vectors are identified by their component number */
{
  register VECTOR *v;

  for( v=PFIRSTVECTOR(grid); v!=NULL; v=SUCCVC(v) )
    VVALUE(v,dest) = VVALUE(v,source);
}


static INT FFIter (NP_ITER *theNP, INT level,
                   VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                   INT *result)
{
#ifdef __BLOCK_VECTOR_DESC__
  NP_FF *np;
  BV_DESC bvd;
  GRID *theGrid;
  INT err;
  INT x_save, b_start_save, b_end_save;

  np = (NP_FF *) theNP;
  theGrid = NP_GRID(theNP,level);

  BVD_INIT( &bvd );
  BVD_PUSH_ENTRY( &bvd, BVNUMBER(GFIRSTBV(theGrid)), NPFF_BVDF(np) );

  /* make a copy for displaying */
  np->smoother.iter.c = x;

  if ( NPFF_CheckSymm(np) )
  {
    x_save =  GET_AUX_VEC;
    b_start_save =  GET_AUX_VEC;
    b_end_save =  GET_AUX_VEC;

    FFCopyVector( theGrid, b_start_save, VD_SCALCMP( b ) );
  }

  if( (err = FFApplyPreconditioner( np, level, x, b, A, result, &bvd, theGrid)) != 0 )
    REP_ERR_RETURN (err);

  if ( NPFF_CheckSymm(np) )
  {
    /* the idea for the symmetry check: if M^-1 is symmetric, then for any
       2 arbitrary vectors (M^-1 * a, b) == ( a, M^-1 * b ) must hold,
       especially for a:=M^-1*d and b:=d. Thus we check:
       (M^-1 * M^-1*d, d) == ( M^-1*d, M^-1 * d ). */

    DOUBLE norm1, norm2;
    register VECTOR *v;
    static INT nr_of_call=0;

    nr_of_call++;

    /* save x and b, the actual results of this routine */
    FFCopyVector( theGrid, b_end_save, VD_SCALCMP( b ) );
    FFCopyVector( theGrid, x_save, VD_SCALCMP( x ) );

    /* norm2 := ( M^-1*d, M^-1 * d ) notice: x contains in this Moment M^-1*d, the approximate solution of Mx=d with initial guess x=0 */
    /* make x collect to can apply ddot; consider now x is consistent!  */
    for( v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v) )
      if(PRIO(v)==PrioBorder)
        VVALUE(v,VD_SCALCMP( x )) = 0.0;
    if( ddot(NP_MG((NP_ITER*)np),level,level,ALL_VECTORS,x,x,&norm2) != NUM_OK )
      REP_ERR_RETURN (1);

    /* set up for solving Mx=M^-1*d */
    /* note: x already has been made inconsistent, which is necessary for b */
    FFCopyVector( theGrid, VD_SCALCMP( b ), VD_SCALCMP( x ) );
    dsetBS( GFIRSTBV(theGrid), VD_SCALCMP( x ), 0.0 );

    /* solving Mx=M^-1*d */
    UserWrite("Solving with FF for symmetry check (A):\n");
    if( (err = FFApplyPreconditioner( np, level, x, b, A, result, &bvd, theGrid)) != 0 )
      REP_ERR_RETURN (err);

    /* norm1 := (M^-1 * M^-1*d, d), x containing M^-1 * M^-1*d, b containing d */
    FFCopyVector( theGrid, VD_SCALCMP( b ), b_start_save );
    if( ddot(NP_MG((NP_ITER*)np),level,level,ALL_VECTORS,b,x,&norm1) != NUM_OK )
      REP_ERR_RETURN (1);

    /* check whether the 2 norms are equal or not
       use a relative criterion */
    if ( fabs( (norm1-norm2) / (norm1+norm2) ) > 1e-5 )
    /* do not expect too accurate identity */
    {
      UserWriteF( "(A) FF preconditioner is NOT symmetric: (M^-1M^-1d,d)=%17.15g<>%17.15g=(M^-1d,M^-1d), difference=%17.15g\n", norm1, norm2, fabs(norm1-norm2) );
    }
    else
    {
      UserWriteF( "(A) FF preconditioner is symmetric: (M^-1M^-1d,d)=%17.15g==%17.15g=(M^-1d,M^-1d)\n", norm1, norm2 );
    }

    /* second kind of check: arbitrary, but fixed vectors; caled after the first call by x resp. b */
    FFGenerateCheckA( GFIRSTBV(theGrid), VD_SCALCMP( b ), nr_of_call, x_save );
    dsetBS( GFIRSTBV(theGrid), VD_SCALCMP( x ), 0.0 );
    UserWrite("Solving with FF for symmetry check (B):\n");
    if( (err = FFApplyPreconditioner( np, level, x, b, A, result, &bvd, theGrid)) != 0 )
      REP_ERR_RETURN (err);
    FFGenerateCheckB( GFIRSTBV(theGrid), VD_SCALCMP( b ), nr_of_call, b_start_save );
    if( ddot(NP_MG((NP_ITER*)np),level,level,ALL_VECTORS,b,x,&norm1) != NUM_OK )
      REP_ERR_RETURN (1);

    FFGenerateCheckB( GFIRSTBV(theGrid), VD_SCALCMP( b ), nr_of_call, b_start_save );
    dsetBS( GFIRSTBV(theGrid), VD_SCALCMP( x ), 0.0 );
    UserWrite("Solving with FF for symmetry check (B):\n");
    if( (err = FFApplyPreconditioner( np, level, x, b, A, result, &bvd, theGrid)) != 0 )
      REP_ERR_RETURN (err);
    FFGenerateCheckA( GFIRSTBV(theGrid), VD_SCALCMP( b ), nr_of_call, x_save );
    if( ddot(NP_MG((NP_ITER*)np),level,level,ALL_VECTORS,b,x,&norm2) != NUM_OK )
      REP_ERR_RETURN (1);

    /* check whether the 2 norms are equal or not
       use a relative criterion */
    if ( fabs( (norm1-norm2) / (norm1+norm2) ) > 1e-5 )
    /* do not expect too accurate identity */
    {
      UserWriteF( "(B) FF preconditioner is NOT symmetric: (M^-1a,b)=%17.15g<>%17.15g=(a,M^-1b), difference=%17.15g\n", norm1, norm2, fabs(norm1-norm2) );
    }
    else
    {
      UserWriteF( "(B) FF preconditioner is symmetric: (M^-1a,b)=%17.15g==%17.15g=(a,M^-1b)\n", norm1, norm2 );
    }

    /* restore the saved values */
    FFCopyVector( theGrid, VD_SCALCMP( b ), b_end_save );
    FFCopyVector( theGrid, VD_SCALCMP( x ), x_save );

    FREE_AUX_VEC(b_end_save);
    FREE_AUX_VEC(b_start_save);
    FREE_AUX_VEC(x_save);
  }

  /* set all vectors with VCLASS < ACTIVE_CLASS to 0.0
     since in the FF routines all inner vectors are calculated, only the
     dirichlet boundary vectors remain as < ACTIVE_CLASS
     BVSUCC(GFIRSTBV(theGrid)) are exactly the dirichlet boundary vectors */
  dsetBS( BVSUCC(GFIRSTBV(theGrid)), VD_SCALCMP( x ), 0.0 );

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
   This numproc executes a linear multigrid cycle.

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
  INT i,ret;
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
  if (np->baselevel<0)
  {
    for (i=FULLREFINELEVEL(NP_MG(theNP)); i>0; i--)
      if (NVEC(GRID_ON_LEVEL(NP_MG(theNP),i))<=-np->baselevel)
        break;
    np->baselevel=i;
  }

  if (np->Transfer == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PreSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->PostSmooth == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->BaseSolver == NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);

  ret = NPIterInit(&np->iter,argc,argv);

  if (sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      np->damp[i] = 1.0;

  return (ret);
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
  UserWriteF(DISPLAY_NP_FORMAT_SI,"baselevel",(int)np->baselevel);

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

  if (sc_disp(np->damp,np->iter.b,"damp"))
    REP_ERR_RETURN (1);

  return (0);
}

static INT LmgcPreProcess  (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *baselevel, INT *result)
{
  NP_LMGC *np = (NP_LMGC *) theNP;
  INT i;

  if (np->Transfer->PreProcess != NULL)
    if ((*np->Transfer->PreProcess)
          (np->Transfer,&(np->baselevel),level,x,b,A,result))
      REP_ERR_RETURN(1);

  if (np->PreSmooth->PreProcess != NULL)
    for (i = np->baselevel+1; i <= level; i++)
      if ((*np->PreSmooth->PreProcess)
            (np->PreSmooth,i,x,b,A,baselevel,result))
        REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PreProcess != NULL)
      for (i = np->baselevel+1; i <= level; i++)
        if ((*np->PostSmooth->PreProcess)
              (np->PostSmooth,i,x,b,A,baselevel,result))
          REP_ERR_RETURN(1);

  *baselevel = MIN(np->baselevel,level);
  if (np->gamma>0 && np->BaseSolver->PreProcess != NULL)
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
  DOUBLE eunorm;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = c;
  NPIT_b(theNP) = b;

  np = (NP_LMGC *) theNP;
  theMG = NP_MG(theNP);
  theGrid = GRID_ON_LEVEL(theMG,level);

  if (level <= np->baselevel) {
    if ((*np->BaseSolver->Residuum)
          (np->BaseSolver,MIN(level,np->baselevel),level,c,b,A,&lresult))
      REP_ERR_RETURN(1);

    IFDEBUG(np,4)
                #ifdef ModelP
    if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
                #endif
    dnrm2(theMG,level,level,ALL_VECTORS,b,&eunorm);
    UserWriteF("defect before base solver : %f\n",eunorm);
    dnrm2(theMG,level,level,ALL_VECTORS,c,&eunorm);
    UserWriteF("norm before base solver : %f\n",eunorm);
    ddot(theMG,level,level,ALL_VECTORS,b,c,&eunorm);
    UserWriteF("c*b before base solver : %f\n",eunorm);
    ENDDEBUG

    if ((*np->BaseSolver->Solver)(np->BaseSolver,level,c,b,A,
                                  np->BaseSolver->abslimit,
                                  np->BaseSolver->reduction,&lresult))
      NP_RETURN(1,result[0]);

    IFDEBUG(np,4)
                #ifdef ModelP
    if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
                #endif
    dnrm2(theMG,level,level,ALL_VECTORS,b,&eunorm);
    UserWriteF("defect after base solver : %f\n",eunorm);
    dnrm2(theMG,level,level,ALL_VECTORS,c,&eunorm);
    UserWriteF("norm after base solver : %f\n",eunorm);
    ddot(theMG,level,level,ALL_VECTORS,b,c,&eunorm);
    UserWriteF("c*b after base solver : %f\n",eunorm);
    ENDDEBUG

    /*
       if (!lresult.converged)
            PrintErrorMessage('W',"Lmgc","no convergence of BaseSolver");
     */
    return(0);
  }

  IFDEBUG(np,4)
        #ifdef ModelP
  if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  dnrm2(theMG,level,level,ALL_VECTORS,b,&eunorm);
  UserWriteF("defect before smoothing : %f\n",eunorm);
  ENDDEBUG

  if (AllocVDFromVD(theMG,level,level,c,&np->t)) NP_RETURN(1,result[0]);
  for (i=0; i<np->nu1; i++) {
    if ((*np->PreSmooth->Iter)(np->PreSmooth,level,np->t,b,A,result))
      REP_ERR_RETURN(1);
    if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
  }

  IFDEBUG(np,4)
        #ifdef ModelP
  if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,b,&eunorm);
  UserWriteF("defect after presmoothing : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,c,&eunorm);
  UserWriteF("correction of presmoothing : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,c,&eunorm);
  UserWriteF("c*b after presmoothing : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,np->t,&eunorm);
  UserWriteF("last correction update of presmoothing : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,np->t,&eunorm);
  UserWriteF("last t*b of presmoothing : %f\n",eunorm);
  ENDDEBUG

#ifdef USE_FAMG
  if((*np->Transfer->RestrictDefect) == FAMGRestrictDefect)
  {
    ((NP_FAMG_TRANSFER*)np->Transfer)->smooth_globsol = c;
    ((NP_FAMG_TRANSFER*)np->Transfer)->smooth_sol = np->t;
    ((NP_FAMG_TRANSFER*)np->Transfer)->smooth_def = b;
  }
#endif

  if ((*np->Transfer->RestrictDefect)
        (np->Transfer,level,b,b,A,Factor_One,result))
    REP_ERR_RETURN(1);

#ifdef USE_FAMG
  IFDEBUG(np,4)
        #ifdef ModelP
  if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  dnrm2(theMG,level,level,ALL_VECTORS,b,&eunorm);
  UserWriteF("defect on fine grid after restriction : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,np->t,&eunorm);
  UserWriteF("correction update on fine grid after restriction : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,c,&eunorm);
  UserWriteF("correction on fine grid after restriction : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,c,&eunorm);
  UserWriteF("c*b on fine grid after restriction : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,np->t,&eunorm);
  UserWriteF("t*b on fine grid after restriction : %f\n",eunorm);
  ENDDEBUG
#endif

  IFDEBUG(np,4)
        #ifdef ModelP
  if (l_vector_collect(DOWNGRID(theGrid),b) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  dnrm2(theMG,level-1,level-1,ALL_VECTORS,b,&eunorm);
  UserWriteF("defect on coarse grid after restriction : %f\n",eunorm);
  ENDDEBUG

  if (dset(theMG,level-1,level-1,ALL_VECTORS,c,0.0) != NUM_OK)
    NP_RETURN(1,result[0]);
  for (i=0; i<np->gamma; i++)
    if (Lmgc(theNP,level-1,c,b,A,result))
      REP_ERR_RETURN(1);

#ifdef USE_FAMG
  if((*np->Transfer->InterpolateCorrection) == FAMGInterpolateCorrection)
  {
    ((NP_FAMG_TRANSFER*)np->Transfer)->smooth_sol = np->t;
    ((NP_FAMG_TRANSFER*)np->Transfer)->smooth_def = b;
  }
#endif
  if ((*np->Transfer->InterpolateCorrection)
        (np->Transfer,level,np->t,c,A,np->damp,result))
    REP_ERR_RETURN(1);
  IFDEBUG(np,4)
        #ifdef ModelP
  if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  dnrm2(theMG,level,level,ALL_VECTORS,b,&eunorm);
  UserWriteF("defect after prolongation : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,np->t,&eunorm);
  UserWriteF("norm of interpolated correction update: %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,c,&eunorm);
  UserWriteF("c*b after prolongation : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,np->t,&eunorm);
  UserWriteF("t*b after prolongation : %f\n",eunorm);
  ENDDEBUG

  if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dmatmul_minus(theMG,level,level,ALL_VECTORS,b,A,np->t) != NUM_OK)
    NP_RETURN(1,result[0]);

  IFDEBUG(np,4)
        #ifdef ModelP
  if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  dnrm2(theMG,level,level,ALL_VECTORS,b,&eunorm);
  UserWriteF("defect after CG correction : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,c,&eunorm);
  UserWriteF("correction after CG correction : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,c,&eunorm);
  UserWriteF("c*b after CG correction : %f\n",eunorm);
  ENDDEBUG


  for (i=0; i<np->nu2; i++) {
    if ((*np->PostSmooth->Iter)(np->PostSmooth,level,np->t,b,A,result))
      REP_ERR_RETURN(1);
    if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }

  IFDEBUG(np,4)
        #ifdef ModelP
  if (l_vector_collect(theGrid,b) != NUM_OK) NP_RETURN(1,result[0]);
        #endif
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,b,&eunorm);
  UserWriteF("defect after postsmoothing : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,c,&eunorm);
  UserWriteF("correction after postsmoothing : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,c,&eunorm);
  UserWriteF("c*b after postsmoothing : %f\n",eunorm);
  dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,np->t,&eunorm);
  UserWriteF("last correction update of postsmoothing : %f\n",eunorm);
  ddot(theMG,level,level,ALL_VECTORS,b,np->t,&eunorm);
  UserWriteF("last t*b of presmoothing : %f\n",eunorm);
  ENDDEBUG

  if (FreeVD(NP_MG(theNP),level,level,np->t)) REP_ERR_RETURN(1);
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

  if (np->gamma>0 && np->BaseSolver->PostProcess != NULL)
    if ((*np->BaseSolver->PostProcess)
          (np->BaseSolver,np->baselevel,x,b,A,result))
      REP_ERR_RETURN(1);

  if (np->PreSmooth != np->PostSmooth)
    if (np->PostSmooth->PostProcess != NULL)
      for (i = level; i >= np->baselevel+1; i--)
        if ((*np->PostSmooth->PostProcess)
              (np->PostSmooth,i,x,b,A,result))
          REP_ERR_RETURN(1);

  if (np->PreSmooth->PostProcess != NULL)
    for (i = level; i >= np->baselevel+1; i--)
      if ((*np->PreSmooth->PostProcess)
            (np->PreSmooth,i,x,b,A,result))
        REP_ERR_RETURN(1);

  if (np->Transfer->PostProcess != NULL)
    if ((*np->Transfer->PostProcess)
          (np->Transfer,&(np->baselevel),level,x,b,A,result))
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
   This numproc executes an additive multigrid cycle.

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
          (np->Transfer,&(np->baselevel),level,x,b,A,result))
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
  INT i;
  INT mylevel;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = c;
  NPIT_b(theNP) = b;

  np = (NP_LMGC *) theNP;
  theMG = NP_MG(theNP);
  for (mylevel=level; mylevel > np->baselevel ; mylevel--)
    if ((*np->Transfer->RestrictDefect)
          (np->Transfer,mylevel,b,b,A,Factor_One,result))
      REP_ERR_RETURN(1);
    #ifdef ModelP
  if (a_vector_collect(theMG,np->baselevel,level,b)!=NUM_OK)
    NP_RETURN(1,result[0]);
  /* TODO: supress collect and consistent routines in smoother */
    #endif
  for (mylevel=np->baselevel; mylevel < level ; mylevel++) {
    theGrid = GRID_ON_LEVEL(theMG,mylevel);
    if (AllocVDFromVD(theMG,mylevel,mylevel,c,&np->t))
      NP_RETURN(1,result[0]);
    for (i=0; i<np->nu1; i++) {
      if ((*np->PreSmooth->Iter)
            (np->PreSmooth,mylevel,np->t,b,A,result))
        REP_ERR_RETURN(1);
      if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK)
        NP_RETURN(1,result[0]);
    }
    if (FreeVD(NP_MG(theNP),mylevel,mylevel,np->t)) REP_ERR_RETURN(1);
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
    if (dadd(theMG,level,level,ALL_VECTORS,c,np->t) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (dmatmul_minus(theMG,level,level,ALL_VECTORS,b,A,np->t)!= NUM_OK)
      NP_RETURN(1,result[0]);

    if (FreeVD(NP_MG(theNP),mylevel,mylevel,np->t)) REP_ERR_RETURN(1);
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
          (np->Transfer,&(np->baselevel),level,x,b,A,result))
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
/*D
   ex - numproc for exact solver

   DESCRIPTION:
   This numproc executes an exact band-lu decomposition.

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

   .  $p - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT EXInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EX *np;

  np = (NP_EX *) theNP;
  np->fmode = ReadArgvOption ("f",argc,argv);
  if (ReadArgvINT ("o",&np->optimizeBand,argc,argv)) np->optimizeBand=1;
  if (ReadArgvINT ("copyback",&np->CopyBack,argc,argv)) np->CopyBack=0;
  if (ReadArgvINT ("printdecomptime",&np->printdecomptime,argc,argv)) np->printdecomptime=0;
  np->nv = -1;
  np->count = -1;

  return (SmootherInit(theNP,argc,argv));
}

static INT EXDisplay (NP_BASE *theNP)
{
  NP_EX *np;
  DOUBLE mem;
  char name [32];

  np = (NP_EX *) theNP;
  SmootherDisplay(theNP);
  if (np->mem > MBYTE)
  {
    strcpy(name,"memory(MByte)");
    mem = (np->mem)/((DOUBLE)MBYTE);
  }
  else if (np->mem > KBYTE)
  {
    strcpy(name,"memory(KByte)");
    mem = (np->mem)/((DOUBLE)KBYTE);
  }
  else
  {
    strcpy(name,"memory(Byte)");
    mem = (np->mem);
  }
  UserWriteF(DISPLAY_NP_FORMAT_SF,name,(float)mem);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"optimize",BOOL_2_YN(np->optimizeBand));
  UserWriteF(DISPLAY_NP_FORMAT_SS,"copy back",BOOL_2_YN(np->CopyBack));
  UserWriteF(DISPLAY_NP_FORMAT_SI,"nv",(int)np->nv);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"count",(int)np->count);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"print decomp time",(int)np->printdecomptime);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"f-mode",(int)np->fmode);

  return (0);
}

/****************************************************************************/
/*D
   EXCopyMatrixFLOAT - copy ug MATRIX to band matrix (FLOAT numbers)

   SYNOPSIS:
   static INT EXCopyMatrixFLOAT (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, FLOAT *Mat);

   PARAMETERS:
   .  theGrid - grid which holds the matrix
   .  x - vector data descriptor for row vector
   .  A - matrix data descriptor
   .  bw - bandwidth
   .  Mat - pointer to FLOAT array to get the bandmatrix

   DESCRIPTION:
   This function copies an ug MATRIX to an array; it is stored as a band matrix.

   'Mat' must provide enough memory!
   At least (2*bw*number_rows + number_cols)*sizeof(FLOAT).

   RETURN VALUE:
   INT  0: o.k.

   SEE ALSO:
   EXCopyMatrixDOUBLE, EXDecomposeMatrixFLOAT, EXApplyLUFLOAT
   D*/
/****************************************************************************/

static INT EXCopyMatrixFLOAT (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, FLOAT *Mat)
{
  INT ment,index,rindex,rtype,rcomp,cindex,ctype,ccomp,i,j;
  VECTOR *theV,*theW;
  MATRIX *theM;
  SHORT *comp;

        #ifdef ModelP
  if (FIRSTVECTOR(theGrid) == NULL)
    return(0);
        #endif

        #ifdef _SPARSE_
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    SPARSE_MATRIX *sm;
    rindex = VINDEX(theV);
    rtype = VTYPE(theV);
    rcomp = VD_NCMPS_IN_TYPE(x,rtype);
    for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
    {
      theW = MDEST(theM);
      cindex = VINDEX(theW);
      ctype = VTYPE(theW);
      sm = (MDIAG(theM)) ? A->sm[DMTP(rtype)] : A->sm[MTP(rtype,ctype)];
      ccomp = VD_NCMPS_IN_TYPE(x,ctype);
      for (i=0; i<rcomp; i++)
        for (j=sm->row_start[i]; j<sm->row_start[i+1]; j++)
          EX_MAT(Mat,bw,rindex+i,cindex+sm->col_ind[j])
            = MVALUE(theM,sm->offset[j]);
    }
  }
        #else
  if (MD_IS_SCALAR(A))
  {
    ment = MD_SCALCMP(A);
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      index = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      if (rcomp == 0) continue;
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM)) {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        if (ccomp == 0) continue;
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        EX_MAT(Mat,bw,index,MDESTINDEX(theM)) = MVALUE(theM,ment);
      }
    }
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      rindex = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        for (i=0; i<rcomp; i++)
          for (j=0; j<ccomp; j++)
            EX_MAT(Mat,bw,rindex+i,cindex+j) = MVALUE(theM,comp[i*ccomp+j]);
      }
    }
  }
        #endif
  return (0);
}

/****************************************************************************/
/*D
   EXCopyMatrixFLOATback - copy ug MATRIX to band matrix (FLOAT numbers)

   SYNOPSIS:
   static INT EXCopyMatrixFLOATback (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, FLOAT *Mat);

   PARAMETERS:
   .  theGrid - grid which holds the matrix
   .  x - vector data descriptor for row vector
   .  A - matrix data descriptor
   .  bw - bandwidth
   .  Mat - pointer to FLOAT array to get the bandmatrix

   DESCRIPTION:
   This function copies an array to an ug MATRIX. It is meant for use
   with the Matrix plot object for visualization of the decomposed matrix.
   CAUTION: the entries not in the ug-MATRIX-pattern are dropped.

   RETURN VALUE:
   INT  0: o.k.

   SEE ALSO:
   EXCopyMatrixDOUBLE, EXDecomposeMatrixFLOAT, EXApplyLUFLOAT
   D*/
/****************************************************************************/

static INT EXCopyMatrixFLOATback (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, FLOAT *Mat)
{
#ifndef _SPARSE_
  INT ment,index,rindex,rtype,rcomp,cindex,ctype,ccomp,i,j;
  VECTOR *theV,*theW;
  MATRIX *theM;
  SHORT *comp;

        #ifdef ModelP
  if (FIRSTELEMENT(theGrid) == NULL)
    return(0);
        #endif

  if (MD_IS_SCALAR(A))
  {
    ment = MD_SCALCMP(A);
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      index = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      if (rcomp == 0) continue;
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM)) {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        if (ccomp == 0) continue;
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        MVALUE(theM,ment) = EX_MAT(Mat,bw,index,MDESTINDEX(theM));
      }
    }
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      rindex = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        for (i=0; i<rcomp; i++)
          for (j=0; j<ccomp; j++)
            MVALUE(theM,comp[i*ccomp+j]) = EX_MAT(Mat,bw,rindex+i,cindex+j);
      }
    }
  }
  return (0);
#else
  PrintErrorMessage ('E', "ex",
                     "No copy-back for sparse matrix blocks implemented");
  return(__LINE__);
#endif
}

/****************************************************************************/
/*D
   EXCopyMatrixDOUBLE - copy ug MATRIX to band matrix (DOUBLE numbers)

   SYNOPSIS:
   static INT EXCopyMatrixDOUBLE (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, DOUBLE *Mat);

   PARAMETERS:
   .  theGrid - grid which holds the matrix
   .  x - vector data descriptor for row vector
   .  A - matrix data descriptor
   .  bw - bandwidth
   .  Mat - pointer to DOUBLE array to get the bandmatrix

   DESCRIPTION:
   This function copies an ug MATRIX to an array; it is stored as a band matrix.

   'Mat' must provide enough memory!
   At least (2*bw*number_rows + number_cols)*sizeof(DOUBLE).

   RETURN VALUE:
   INT  0: o.k.

   SEE ALSO:
   EXCopyMatrixFLOAT, EXDecomposeMatrixDOUBLE, EXApplyLUDOUBLE
   D*/
/****************************************************************************/

static INT EXCopyMatrixDOUBLE (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, DOUBLE *Mat)
{
  INT ment,index,rindex,rtype,rcomp,cindex,ctype,ccomp,i,j;
  VECTOR *theV,*theW;
  MATRIX *theM;
  SHORT *comp;
#ifdef _SPARSE_
  SPARSE_MATRIX *sm;
#endif

        #ifdef ModelP
  if (FIRSTVECTOR(theGrid) == NULL)
    return(0);
        #endif

        #ifdef _SPARSE_
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    rindex = VINDEX(theV);
    rtype = VTYPE(theV);
    rcomp = VD_NCMPS_IN_TYPE(x,rtype);
    for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
    {
      theW = MDEST(theM);
      cindex = VINDEX(theW);
      ctype = VTYPE(theW);
      sm = (MDIAG(theM)) ? A->sm[DMTP(rtype)] : A->sm[MTP(rtype,ctype)];
      ccomp = VD_NCMPS_IN_TYPE(x,ctype);
      for (i=0; i<rcomp; i++)
        for (j=sm->row_start[i]; j<sm->row_start[i+1]; j++)
          EX_MAT(Mat,bw,rindex+i,cindex+sm->col_ind[j])
            = MVALUE(theM,sm->offset[j]);
    }
  }
        #else
  if (MD_IS_SCALAR(A))
  {
    ment = MD_SCALCMP(A);
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      index = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      if (rcomp == 0) continue;
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM)) {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        if (ccomp == 0) continue;
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        EX_MAT(Mat,bw,index,MDESTINDEX(theM)) = MVALUE(theM,ment);
      }
    }
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      rindex = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        for (i=0; i<rcomp; i++)
          for (j=0; j<ccomp; j++)
            EX_MAT(Mat,bw,rindex+i,cindex+j) = MVALUE(theM,comp[i*ccomp+j]);
      }
    }
  }
        #endif
  return (0);
}

/****************************************************************************/
/*D
   EXCopyMatrixDOUBLEback - copy ug MATRIX to band matrix (DOUBLE numbers)

   SYNOPSIS:
   static INT EXCopyMatrixDOUBLEback (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, const DOUBLE *Mat);

   PARAMETERS:
   .  theGrid - grid which holds the matrix
   .  x - vector data descriptor for row vector
   .  A - matrix data descriptor
   .  bw - bandwidth
   .  Mat - pointer to DOUBLE array to get the bandmatrix

   DESCRIPTION:
   This function copies an array to an ug MATRIX. It is meant for use
   with the Matrix plot object for visualization of the decomposed matrix.
   CAUTION: the entries not in the ug-MATRIX-pattern are dropped.

   RETURN VALUE:
   INT  0: o.k.

   SEE ALSO:
   EXCopyMatrixFLOAT, EXDecomposeMatrixDOUBLE, EXApplyLUDOUBLE
   D*/
/****************************************************************************/

static INT EXCopyMatrixDOUBLEback (GRID *theGrid, VECDATA_DESC *x, MATDATA_DESC *A, INT bw, const DOUBLE *Mat)
{
#ifndef _SPARSE_
  INT ment,index,rindex,rtype,rcomp,cindex,ctype,ccomp,i,j;
  VECTOR *theV,*theW;
  MATRIX *theM;
  SHORT *comp;

        #ifdef ModelP
  if (FIRSTELEMENT(theGrid) == NULL)
    return(0);
        #endif

  if (MD_IS_SCALAR(A))
  {
    ment = MD_SCALCMP(A);
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      index = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      if (rcomp == 0) continue;
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM)) {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        if (ccomp == 0) continue;
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        MVALUE(theM,ment) = EX_MAT(Mat,bw,index,MDESTINDEX(theM));
      }
    }
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      rindex = VINDEX(theV);
      rtype = VTYPE(theV);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      for (theM=VSTART(theV); theM!=NULL; theM=MNEXT(theM))
      {
        theW = MDEST(theM);
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        for (i=0; i<rcomp; i++)
          for (j=0; j<ccomp; j++)
            MVALUE(theM,comp[i*ccomp+j]) = EX_MAT(Mat,bw,rindex+i,cindex+j);
      }
    }
  }
  return (0);
#else
  PrintErrorMessage ('E', "ex",
                     "No copy-back for sparse matrix blocks implemented");
  return(__LINE__);
#endif
}

static INT EXPreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_EX *np = (NP_EX *) theNP;
  FIFO myfifo;
  void *buffer;
  VECTOR **vlist;
  VECTOR *theV;
  MATRIX *theM;
  HEAP *theHeap = MGHEAP(NP_MG(theNP));
  GRID *theGrid = NP_GRID(theNP,level);
  INT bw,i,k,index,max, outOfMem, globalOutOfMem;
  INT MarkKey;
  INT optimizeBand = np->optimizeBand;
  INT n = 0;
  DOUBLE decomptime = -1.0;

  if( np->printdecomptime )
    decomptime = CURRENT_TIME;

  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VD_NCMPS_IN_TYPE(x,VTYPE(theV)) > 0) n++;
  np->nv = n;
  np->pp_failed=0;
  if (n == 0)
  {
#ifdef ModelP
    /* exit premature; but nevertheless participate all communications in this function */
    outOfMem = 0;

    /* i counts the number of communications */
    i = 1;                      /* mem for matrix copy */
    if( (np->count+1)==0 )
      i++;
    if (np->optimizeBand)
      i++;
    for( ; i>0; i-- )
    {
      globalOutOfMem = UG_GlobalMaxINT( outOfMem );
      if( globalOutOfMem )
        REP_ERR_RETURN(1);
    }
#endif
    return(0);
  }

  *baselevel = level;
  if (np->count >= 0)
    optimizeBand = 0;
  if (np->optimizeBand)
  {
    IFDEBUG(np,1)
    INT ne=0;
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
        if (!CEXTRA(MMYCON(theM)))
          ne++;
    if (ne)
    {
#ifdef Debug
      PrintDebug("WARNING: %d extra connections found by ex!\n",(int)ne);
#endif
      PrintErrorMessageF('W',"EXPreProcess","%d extra connections found by ex!\n",(int)ne);
    }
    ENDDEBUG

    /* reorder vector-list */
    MarkTmpMem(theHeap,&MarkKey);
    buffer=(void *)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);
    vlist = (VECTOR**)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);

    outOfMem = ( buffer==NULL || vlist==NULL );
#ifdef ModelP
    globalOutOfMem = UG_GlobalMaxINT( outOfMem );
#else
    globalOutOfMem = outOfMem;
#endif
    if( globalOutOfMem )
    {
#ifdef ModelP
      if( outOfMem )
      {
        PrintErrorMessageF('W',"EXPreProcess","cannot allocate mem for %d vectors on processor %d for reorder vector-list\n", n, (int)me);
      }
#else
      PrintErrorMessageF('W',"EXPreProcess","cannot allocate mem for %d vectors for reorder vector-list\n", n);
#endif
      REP_ERR_RETURN(1);
    }

    fifo_init(&myfifo,buffer,sizeof(VECTOR*)*n);
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
      SETVCUSED(theV,0);
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
      if (VD_NCMPS_IN_TYPE(x,VTYPE(theV)) > 0)
        break;
    fifo_in(&myfifo,theV);
    SETVCUSED(FIRSTVECTOR(theGrid),1);
    while(!fifo_empty(&myfifo))
    {
      theV = (VECTOR *)fifo_out(&myfifo);
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
        if (!VCUSED(MDEST(theM)))
        {
          if (VD_NCMPS_IN_TYPE(x,VTYPE(MDEST(theM))) == 0)
            continue;
          fifo_in(&myfifo,(void *)MDEST(theM));
          SETVCUSED(MDEST(theM),1);
        }
    }
    fifo_in(&myfifo,(void *)theV);
    SETVCUSED(theV,0); i=0;
    while(!fifo_empty(&myfifo))
    {
      theV = (VECTOR *)fifo_out(&myfifo);
      vlist[i++] = theV;
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
        if (VCUSED(MDEST(theM)))
        {
          if (VD_NCMPS_IN_TYPE(x,VTYPE(MDEST(theM))) == 0)
            continue;
          fifo_in(&myfifo,(void *)MDEST(theM));
          SETVCUSED(MDEST(theM),0);
        }
    }
    assert(i==n);
    for (i=0; i<n; i++) GRID_UNLINK_VECTOR(theGrid,vlist[i]);
    for (i=0; i<n; i++) GRID_LINK_VECTOR(theGrid,vlist[i],PRIO(vlist[i]));
    ReleaseTmpMem(theHeap,MarkKey);
  }
  if (MD_IS_SCALAR(A))
  {
    k = 0;
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV),i++)
      if (VD_NCMPS_IN_TYPE(x,VTYPE(theV)) > 0)
        VINDEX(theV) = k++;
    bw = 0;
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      index = VINDEX(theV);
      if (VD_NCMPS_IN_TYPE(x,VTYPE(theV)) == 0) continue;
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      {
        if (VD_NCMPS_IN_TYPE(x,VTYPE(MDEST(theM))) == 0) continue;
        k = index-MDESTINDEX(theM);
        k = ABS(k);
        bw = MAX(bw,k);
      }
    }
    np->bw = bw;
  }
  else
  {
    max = 0;
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV))
    {
      VINDEX(theV) = i;
      k = VD_NCMPS_IN_TYPE(x,VTYPE(theV));
      i += k;
      max = MAX(max,k);
    }
    bw = 0;
    for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    {
      index = VINDEX(theV);
      if (VD_NCMPS_IN_TYPE(x,VTYPE(theV)) == 0) continue;
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      {
        if (VD_NCMPS_IN_TYPE(x,VTYPE(MDEST(theM))) == 0) continue;
        k = index-MDESTINDEX(theM);
        k = ABS(k);
        bw = MAX(bw,k);
      }
    }
    np->bw = bw + max - 1;
    np->nv = i;
  }
  if (np->CopyBack)
    /* alloc a decomposed matrix for copying back */
    /* TODO: possibly enlarge pattern to full bandwidth for storing the complete decomposed matrix */
    if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->smoother.L)) REP_ERR_RETURN(1);

  /* get storage for matrix */
  bw = np->bw;
  np->count++;
  if (MarkTmpMem(theHeap,&(np->MarkKey[np->count])))
    REP_ERR_RETURN(1);
  if (np->count == 0)
  {
    np->Vec = (DOUBLE*)GetTmpMem(theHeap,np->nv*sizeof(DOUBLE),np->MarkKey[np->count]);
    outOfMem = ( np->Vec==NULL );
#ifdef ModelP
    globalOutOfMem = UG_GlobalMaxINT( outOfMem );
#else
    globalOutOfMem = outOfMem;
#endif
    if( globalOutOfMem )
    {
#ifdef ModelP
      if( outOfMem )
      {
        PrintErrorMessageF('W',"EXPreProcess","cannot allocate mem for %d vectors on processor %d for np->Vec\n", np->nv, (int)me);
      }
#else
      PrintErrorMessageF('W',"EXPreProcess","cannot allocate mem for %d vectors for np->Vec\n", np->nv);
#endif
      REP_ERR_RETURN(1);
    }
  }


  if (np->fmode == 1)
  {
    np->mem = np->nv*(2*bw+1)*sizeof(FLOAT);
    np->FMat[np->count] = (FLOAT*)GetTmpMem(theHeap,np->mem,np->MarkKey[np->count]);
    outOfMem = ( np->FMat[np->count]==NULL );
#ifdef ModelP
    globalOutOfMem = UG_GlobalMaxINT( outOfMem );
#else
    globalOutOfMem = outOfMem;
#endif
    if( globalOutOfMem )
    {
#ifdef ModelP
      if( outOfMem )
      {
        PrintErrorMessageF('W',"EXPreProcess","cannot allocate %d bytes on processor %d for FMat\n", np->mem, (int)me);
      }
#else
      PrintErrorMessageF('W',"EXPreProcess","cannot allocate %d bytes for FMat\n", np->mem, (int)me);
#endif
      REP_ERR_RETURN(1);
    }

    memset((void*)(np->FMat[np->count]),0,np->mem);
    if (EXCopyMatrixFLOAT (theGrid,x,A,np->bw,np->FMat[np->count]))
      REP_ERR_RETURN(1);
    if (EXDecomposeMatrixFLOAT (np->FMat[np->count],np->bw,np->nv))
      np->pp_failed=1;
    if (np->CopyBack)
      if (EXCopyMatrixFLOATback(theGrid,x,np->smoother.L,np->bw,np->FMat[np->count]))
        REP_ERR_RETURN(1);
  }
  else
  {
    np->mem = np->nv*(2*bw+1)*sizeof(DOUBLE);
    np->DMat[np->count] = (DOUBLE*)GetTmpMem(theHeap,np->mem,np->MarkKey[np->count]);
    outOfMem = ( np->DMat[np->count]==NULL );
#ifdef ModelP
    globalOutOfMem = UG_GlobalMaxINT( outOfMem );
#else
    globalOutOfMem = outOfMem;
#endif
    if( globalOutOfMem )
    {
#ifdef ModelP
      if( outOfMem )
      {
        PrintErrorMessageF('W',"EXPreProcess","cannot allocate %d bytes on processor %d for DMat\n", np->mem, (int)me);
      }
#else
      PrintErrorMessageF('W',"EXPreProcess","cannot allocate %d bytes for DMat\n", np->mem, (int)me);
#endif
      REP_ERR_RETURN(1);
    }

    memset((void*)(np->DMat[np->count]),0,np->mem);
    if (EXCopyMatrixDOUBLE (theGrid,x,A,np->bw,np->DMat[np->count]))
      REP_ERR_RETURN(1);
    if (EXDecomposeMatrixDOUBLE (np->DMat[np->count],np->bw,np->nv))
      np->pp_failed=1;
    if (np->CopyBack)
      if (EXCopyMatrixDOUBLEback(theGrid,x,np->smoother.L,np->bw,np->DMat[np->count]))
        REP_ERR_RETURN(1);
  }

  if( np->printdecomptime )
  {
    decomptime = CURRENT_TIME - decomptime;
                #ifdef ModelP
    if( me == master )
                #endif
    UserWriteF("EX decomp %g sec\n", decomptime );
  }

  return (0);
}

static INT EXSmoother (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  INT i,j,n,bw,vent,type;
  NP_EX *np;
  GRID *theGrid;
  DOUBLE *Vec;
  VECTOR *theV;
  SHORT *comp;

  /* store passed XXXDATA_DESCs */
  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  np = (NP_EX *) theNP;
  theGrid = NP_GRID(theNP,level);

  /* check failure of preprocess */
  if (np->pp_failed)
  {
    if (dset(NP_MG(theNP),level,level,ALL_VECTORS,x,0.0)!= NUM_OK) NP_RETURN(1,result[0]);
    return(1);
  }

  /* init */
  n               = np->nv;
  if (n == 0)
    return(0);

  bw              = np->bw;
  Vec             = np->Vec;

  /* copy b to Vec */
  if (MD_IS_SCALAR(A))
  {
    vent=VD_SCALCMP(b);
    j = 0;
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV),i++)
      if (VD_NCMPS_IN_TYPE(b,VTYPE(theV)) > 0)
        Vec[j++] = VVALUE(theV,vent);
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV))
    {
      type = VTYPE(theV);
      comp = VD_CMPPTR_OF_TYPE(b,type);
      for (j=0; j<VD_NCMPS_IN_TYPE(b,type); j++)
        Vec[i++] = VVALUE(theV,comp[j]);
    }
  }

  /* solve for vector */
  if (np->fmode == 1)
  {
    if (EXApplyLUFLOAT (np->FMat[np->count],bw,n,Vec)) REP_ERR_RETURN(1);
  }
  else
  {
    if (EXApplyLUDOUBLE (np->DMat[np->count],bw,n,Vec)) REP_ERR_RETURN(1);
  }

  /* copy Vec to x */
  if (MD_IS_SCALAR(A))
  {
    vent=VD_SCALCMP(x);
    j = 0;
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV),i++)
      if (VD_NCMPS_IN_TYPE(x,VTYPE(theV)) > 0)
        VVALUE(theV,vent) = Vec[j++];
  }
  else
  {
    for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV))
    {
      type = VTYPE(theV);
      comp = VD_CMPPTR_OF_TYPE(x,type);
      for (j=0; j<VD_NCMPS_IN_TYPE(x,type); j++)
        VVALUE(theV,comp[j]) = Vec[i++];
    }
  }

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->smoother.damp) != NUM_OK) NP_RETURN(1,result[0]);

  /* update defect */
    #ifdef ModelP
  if (l_vector_consistent(theGrid,x) != NUM_OK) NP_RETURN(1,result[0]);
    #endif


  /*

     PrintVector(theGrid,x,3,3);	PrintVector(theGrid,b,3,3);
   */

  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT EXPostProcess (NP_ITER *theNP, INT level,VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_EX *np=(NP_EX*)theNP;
  HEAP *theHeap;

        #ifdef ModelP
  if (FIRSTVECTOR(NP_GRID(theNP,level)) == NULL)
    return(0);
        #endif

  theHeap = MGHEAP(NP_MG(theNP));

  if (np->smoother.L != NULL)
    if (FreeMD(NP_MG(theNP),level,level,np->smoother.L)) REP_ERR_RETURN(1);

  /* we need a better concept ...
     if (np->fmode) {
          ASSERT(np->FMat!=NULL);
          if (PutFreelistMemory(theHeap, np->FMat, np->mem))
              REP_ERR_RETURN(1);
     }
     else {
          ASSERT(np->DMat!=NULL);
          if (PutFreelistMemory(theHeap, np->DMat, np->mem))
              REP_ERR_RETURN(1);
     }

     UserWriteF("post: freelist %d used %d size %d\n",
                     HeapFreelistUsed(theHeap),HeapUsed(theHeap),
                     HeapSize(theHeap));

     ASSERT(np->Vec!=NULL);
     if (PutFreelistMemory(theHeap, np->Vec, np->nv*sizeof(DOUBLE)))
          REP_ERR_RETURN(1);
   */
  ReleaseTmpMem(theHeap,np->MarkKey[np->count]);
  np->FMat[np->count] = NULL;
  np->DMat[np->count] = NULL;
  if (np->count == 0)
    np->Vec = NULL;
  np->count--;

  return(0);
}

static INT EXConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init     = EXInit;
  theNP->Display  = EXDisplay;
  theNP->Execute  = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess  = EXPreProcess;
  np->Iter                = EXSmoother;
  np->PostProcess = EXPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   exact-solver - numproc for exact solver for Neumann-Boundary

   DESCRIPTION:
   This numproc executes exact solver for a Neumann_boundary-Probelm,
   which extend the global stiffness-matrix with the projection
   to the solution-space.

   A_neu = (A_stiff , A_project)^t ; b_neu = (b_stiff , 0_project)

   x_exact = (A_neu^t * A_neu)^(-1) * (A_neu^t * b_neu)

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $P projection ;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $P projection - connect the projection

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT EXPRJInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EXPRJ *np = (NP_EXPRJ *) theNP;

  np->project = (NP_PROJECT *)
                ReadArgvNumProc(theNP->mg,"P",PROJECT_CLASS_NAME,argc,argv);

  np->p = ReadArgvVecDesc(theNP->mg,"p",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);

  return (SmootherInit(theNP,argc,argv));
}

static INT EXPRJDisplay (NP_BASE *theNP)
{
  NP_EXPRJ *np = (NP_EXPRJ *) theNP;

  if (np->p != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"p",ENVITEM_NAME(np->p));

  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));

  if (np->project != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Project",ENVITEM_NAME(np->project));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"Project","---");

  return (0);
}

static INT EXPRJSmoother (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_EXPRJ *np = (NP_EXPRJ *) theNP;
  MULTIGRID *theMG = NP_MG(theNP);
  GRID *theGrid = NP_GRID(theNP,level);
  HEAP *theHeap = MGHEAP(NP_MG(theNP));
  DOUBLE *a, *r, *a_neu, *r_neu, *inv_a_neu, *x_neu;
  DOUBLE *ev[6];
  VECTOR *v;
  MATRIX *m;
  SHORT *Mcomp;
  INT rtype, ctype, rcomp, ccomp, i, j, k, n, vtype, comp, ncomp;
  INT n0,n1,d;
  DOUBLE a0, a1 ;
  INT l,MarkKey;
  /* store passed XXXDATA_DESCs */

  NPIT_A(theNP) = A;
  NPIT_c(theNP) = x;
  NPIT_b(theNP) = b;

  MarkTmpMem(theHeap,&MarkKey);
  if (AllocVDFromVD(theMG,level,level,x,&np->p))
    NP_RETURN(1,result[0]);
  d = np->project->dim;

  /* Counting the Rows */

  n = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    VINDEX(v) = n;
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(A,rtype,rtype);
    n += rcomp;
  }
  n0 = n ;
  n1 = n + d;

  a = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n0 * n1,MarkKey);
  for (i=0; i<n0*n1; i++)
    a[i] = 0.0;

  r = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n1,MarkKey);
  for (i=0; i<n1; i++)
    r[i] = 0.0;

  for (j=0; j<d; j++) {
    ev[j] = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n0,MarkKey);
    for (i=0; i<n0; i++)
      ev[j][i] = 0.0;

    if ((*np->project->ProjectionVector)
          (np->project,level,level,j,np->p,result))
      NP_RETURN(1,result[0]);
    n = 0;
    for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      ncomp = VD_NCMPS_IN_TYPE(np->p,vtype);
      if (ncomp == 0) continue;
      comp =  VD_CMP_OF_TYPE(np->p,vtype,0);
      for (l=0; l<ncomp; l++) {
        ev[j][n] = VVALUE(v,comp+l);
        n++;
      }
    }

    if (j!=0) {
      for (k=0 ; k<j ; k++) {
        a0=0.0;
        a1=0.0;
        for (l=0 ; l<n0 ; l++) {
          a1 += ev[k][l] * ev[j][l] ;
          a0 += ev[k][l] * ev[k][l] ;
        }
        for (l=0 ; l<n0 ; l++) {
          ev[j][l] = ev[j][l] - (a1 / a0) * ev[k][l] ;
        }
      }
    }
  }

  IFDEBUG(np,4)
  {
    for (j=0; j<d; j++) {
      UserWriteF("%d-ten Eigenvektor : \n",j);
      for (l = 0 ; l < n0 ; l++)
        UserWriteF("%f \n",ev[j][l]);
    }
  }
  ENDDEBUG

    n = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    ncomp = VD_NCMPS_IN_TYPE(b,vtype);
    if (ncomp == 0) continue;
    comp =  VD_CMP_OF_TYPE(b,vtype,0);
    for (j=0; j<ncomp; j++) {
      r[n] = VVALUE(v,comp+j);
      n++;
    }
  }

  for (i=0 ; i<d ; i++) {
    a0=0.0;
    a1=0.0;
    for (k=0 ; k<n0 ; k++) {
      a1 += ev[i][k] * r[k] ;
      a0 += ev[i][k] * ev[i][k] ;
    }
    for (k=0 ; k<n0 ; k++) {
      r[k] = r[k] - (a1 / a0 ) * ev[i][k] ;
    }
  }

  n = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    ncomp = VD_NCMPS_IN_TYPE(b,vtype);
    if (ncomp == 0) continue;
    comp =  VD_CMP_OF_TYPE(b,vtype,0);
    for (j=0; j<ncomp; j++) {
      VVALUE(v,comp+j) = r[n] ;
      n++;
    }
  }

  IFDEBUG(np,1)
  {
    DOUBLE norm,norm2;


    IFDEBUG(np,4)
    UserWrite("global matrix at the beginning of EXProjSmoother \n");
    PrintMatrix(theGrid,A,3,3);
    UserWrite("right hand side at the beginning of EXProjSmoother \n");
    PrintVector(theGrid,b,3,3);
    UserWrite("global solution at the beginning of EXProjSmoother \n");
    PrintVector(theGrid,x,3,3);
    ENDDEBUG

    /* calculate the test-vector p with the stiffness-matrix A :
       t = A * p  ; norm = (t,t) = t * t   */

    if (AllocVDFromVD(theMG,level,level,x,&np->t))
      NP_RETURN(1,result[0]);

    for (i=0; i<d; i++) {
      if ((*np->project->ProjectionVector)
            (np->project,level,level,i,np->p,result))
        NP_RETURN(1,result[0]);

      IFDEBUG(np,2)
      UserWriteF("Projektionsvektor %d : \n",i);
      PrintVector(theGrid,np->p,3,3);
      ENDDEBUG

      if (dmatmul(theMG,level,level,ALL_VECTORS,np->t,A,np->p)!= NUM_OK)
        NP_RETURN(1,result[0]);
      if (ddot(theMG,level,level,ON_SURFACE,np->t,np->t,&norm) != NUM_OK)
        NP_RETURN(1,result[0]);
      UserWriteF("Norm zu A vom %d-ten Eigenvektor : %f \n",i,(float)norm);
      if (ddot(theMG,level,level,ON_SURFACE,np->p,b,&norm2) != NUM_OK)
        NP_RETURN(1,result[0]);
      UserWriteF("Norm zu b vom %d-ten Eigenvektor : %f \n",i,(float)norm2);
    }
    FreeVD(theMG,level,level,np->t);
  }
  ENDDEBUG

    n = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    ncomp = VD_NCMPS_IN_TYPE(b,vtype);
    if (ncomp == 0) continue;
    comp =  VD_CMP_OF_TYPE(b,vtype,0);
    for (j=0; j<ncomp; j++) {
      r[n] = VVALUE(v,comp+j);
      n++;
    }
  }
  for (i=0; i<d; i++) {
    if ((*np->project->ProjectionVector)
          (np->project,level,level,i,np->p,result))
      NP_RETURN(1,result[0]);
    n = 0;
    for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      ncomp = VD_NCMPS_IN_TYPE(np->p,vtype);
      if (ncomp == 0) continue;
      comp =  VD_CMP_OF_TYPE(np->p,vtype,0);
      for (j=0; j<ncomp; j++) {
        a[(n0+i)*n0+n] = VVALUE(v,comp+j);
        n++;
      }
    }
  }

  /* stiffness matrix of the Neumannproblem */

  n = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(A,rtype,rtype);
    for (i=0; i<rcomp; i++) {
      for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
        ctype = MDESTTYPE(m);
        ccomp = MD_COLS_IN_RT_CT(A,rtype,ctype);
        if (ccomp == 0) continue;
        Mcomp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        k = VINDEX(MDEST(m));
        for (j=0; j<ccomp; j++) {
          a[n0*n+k+j] = MVALUE(m,Mcomp[i*ccomp+j]);
        }
      }
      n++;
    }
  }

  /* Printing of the modified matrix and the modified right hand side */

  IFDEBUG(np,4)
  UserWriteF("A \n");
  for (i=0; i<n0+d; i++)
  {
    for (j=0; j<n0; j++)
      UserWriteF("%8.3f",a[i*n0+j]);
    UserWriteF("\n");
  }
  UserWriteF("\n");

  UserWriteF("b \n");
  for (i=0; i<n0+d; i++)
  {
    UserWriteF("%8.3f",r[i]);
  }
  UserWriteF("\n");
  ENDDEBUG

    a_neu = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n0 * n0,MarkKey);

  for (i=0; i<n0; i++)
    for (j=0; j<n0; j++)
      a_neu[i*n0+j] = 0.0 ;

  r_neu = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n0,MarkKey);

  for (i=0; i<n0; i++)
    r_neu[i] = 0.0 ;

  /* Calculating of the new stiffness matrix */

  for (i=0; i<n0; i++)
    for (j=0; j<n0; j++)
      for (k=0; k<n1; k++)
        a_neu[i*n0+j] += a[j+n0*k] * a[i+n0*k] ;


  /* Calculating of the new right hand side */

  for (i=0; i<n0; i++)
    for (j=0; j<n1; j++)
      r_neu[i] +=  a[i+n0*j] * r[j] ;



  /* Printing of the new stiffness matrix and new right hand side */

  IFDEBUG(np,4)
  UserWriteF("A_neu \n");
  for (i=0; i<n0; i++)
  {
    for (j=0; j<n0; j++)
      UserWriteF("%8.3f",a_neu[i*n0+j]);
    UserWriteF("\n");
  }
  UserWriteF("\n");

  UserWriteF("b_neu \n");
  for (i=0; i<n0; i++)
  {
    UserWriteF("%8.3f",r_neu[i]);
  }
  UserWriteF("\n");
  ENDDEBUG


  /* Solving the equation-system :  a_neu * x = b_neu */

  /* Calculating of the inverse */

    inv_a_neu = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n0 * n0,MarkKey);

  for (i=0; i<n0; i++)
    for (j=0; j<n0; j++)
      inv_a_neu[i*n0+j] = 0.0 ;

  {
    DOUBLE * rhs_tmp = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n,MarkKey);
    INT    * ipv_tmp = (INT *)   GetTmpMem(theHeap,sizeof(DOUBLE) * n,MarkKey);
    InvertFullMatrix_gen(n0, a_neu, inv_a_neu,rhs_tmp,ipv_tmp);
  }

  IFDEBUG(np,4)
  UserWriteF("inv_a_neu \n");   /* Printing of the inverse matrix */
  for (i=0; i<n0; i++)
  {
    for (j=0; j<n0; j++)
      UserWriteF("%8.3f",inv_a_neu[i*n0+j]);
    UserWriteF("\n");
  }
  UserWriteF("\n");
  ENDDEBUG

    x_neu = (DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE) * n0,MarkKey);

  for (i=0; i<n0; i++)
    x_neu[i] = 0.0 ;

  /* Solving : x = inv_a_neu * r_neu */

  for (i=0; i<n0; i++)
    for (j=0; j<n0; j++)
      x_neu[i] +=  inv_a_neu[i*n0+j] * r_neu[j] ;


  IFDEBUG(np,4)
  UserWriteF("x_neu \n");
  for (i=0; i<n0; i++)
    UserWriteF("%8.3f",x_neu[i]);
  UserWriteF("\n");

  ENDDEBUG

  /* Saving of the solution x_neu into the VECDATA_DESC-Structur */

    n = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    ncomp = VD_NCMPS_IN_TYPE(x,vtype);
    if (ncomp == 0) continue;
    comp =  VD_CMP_OF_TYPE(x,vtype,0);
    for (j=0; j<ncomp; j++) {
      VVALUE(v,comp+j) = x_neu[n] ;
      n++;
    }
  }

  ReleaseTmpMem(theHeap,MarkKey);
  FreeVD(theMG,level,level,np->p);

  /* damp */
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->smoother.damp)
      != NUM_OK)
    NP_RETURN(1,result[0]);

  /* update defect */
    #ifdef ModelP
  if (l_vector_consistent(theGrid,x) != NUM_OK) NP_RETURN(1,result[0]);
    #endif
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)!= NUM_OK)
    NP_RETURN(1,result[0]);

  return (0);
}

static INT EXPRJConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init     = EXPRJInit;
  theNP->Display  = EXPRJDisplay;
  theNP->Execute  = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess  = NULL;
  np->Iter                = EXPRJSmoother;
  np->PostProcess = NULL;

  return(0);
}

/****************************************************************************/
/*D
   calibrate - numproc for additive linear multigrid cycle

   DESCRIPTION:
   This numproc calibrates the damping factors of an iteration.

   .vb
   npinit <name> [$c <cor>] [$r <rhs>] [$A <mat>]
       $I <iter> [$s <tmp1>] [$t <tmp2>]
   .ve

   .  $c~<cor> - correction vector
   .  $r~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $s~<tmp1> - tempory vector
   .  $t~<tmp2> - tempory vector
   .  $I~<iter> - iteration numproc

   'npexecute <name> [$i] [$s] [$p]'

   .  $i - preprocess
   .  $s - solve
   .  $p - postprocess

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT CalibrateInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_CALIBRATE *np = (NP_CALIBRATE *) theNP;
  INT i;

  np->s = ReadArgvVecDesc(theNP->mg,"s",argc,argv);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->Iter = (NP_ITER *)
             ReadArgvNumProc(theNP->mg,"I",ITER_CLASS_NAME,argc,argv);
  if (np->Iter == NULL)
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->Transfer = (NP_TRANSFER *)
                 ReadArgvNumProc(theNP->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  if (ReadArgvINT("n",&(np->n),argc,argv))
    np->n = 1;

  for (i=0; i<2*MAXLEVEL; i++) np->damp[i] = SMALL_D;
  np->display = ReadArgvDisplay(argc,argv);

  return (NPIterInit(&np->iter,argc,argv));
}

static INT CalibrateDisplay (NP_BASE *theNP)
{
  NP_CALIBRATE *np = (NP_CALIBRATE *) theNP;
  INT i;

  NPIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  if (np->Iter != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"I",ENVITEM_NAME(np->Iter));
  if (np->Transfer != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(np->Transfer));
  if (np->s != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(np->s));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->u != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->u));

  for (i=0; i<2*MAXLEVEL; i++)
    if (np->damp[i] != SMALL_D)
      UserWriteF("ev[%3d]         = %-7.4g\n",i-MAXLEVEL,np->damp[i]);

  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",np->n);

  return (0);
}

static INT CalibratePreProcess  (NP_ITER *theNP, INT level,
                                 VECDATA_DESC *x, VECDATA_DESC *b,
                                 MATDATA_DESC *A,
                                 INT *baselevel, INT *result)
{
  NP_CALIBRATE      *np = (NP_CALIBRATE *) theNP;
  MULTIGRID      *theMG = NP_MG(theNP);
  DOUBLE a0,a1;
  INT i;

  if (level == BOTTOMLEVEL(theMG)) return(0);

  if (np->Iter->PreProcess != NULL)
    if ((*np->Iter->PreProcess)(np->Iter,level,x,b,A,baselevel,result))
      REP_ERR_RETURN(1);

  if (AllocVDFromVD(theMG,level-1,level,x,&np->s))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,level,level,x,&np->t))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,level,level,x,&np->u))
    NP_RETURN(1,result[0]);
  np->damp[level+MAXLEVEL] = 0.0;
  for (i=0; i<np->n; i++) {
    if (np->Transfer != NULL) {          /* choose a high frequent vector */
      l_dsetrandom(GRID_ON_LEVEL(theMG,level),np->t,EVERY_CLASS,1.0);
      /* here we assume that the transfer proprocess is done before */
      if (dmatmul(theMG,level,level,ALL_VECTORS,np->s,A,np->t)
          != NUM_OK)
        NP_RETURN(1,result[0]);
      ClearDirichletValues(GRID_ON_LEVEL(theMG,level),np->s);
      if ((*np->Transfer->ProjectSolution)
            (np->Transfer,level-1,level,np->s,result))
        REP_ERR_RETURN(1);
      if ((*np->Transfer->InterpolateCorrection)
            (np->Transfer,level,np->t,np->s,A,Factor_One,result))
        REP_ERR_RETURN(1);
      if (dsub(theMG,level,level,ALL_VECTORS,np->s,np->t) != NUM_OK)
        NP_RETURN(1,result[0]);
    }
    else {                                    /* choose a low frequent vector */
      if (dset(theMG,level,level,ALL_VECTORS,np->s,1.0) != NUM_OK)
        NP_RETURN(1,result[0]);
      ClearDirichletValues(GRID_ON_LEVEL(theMG,level),np->s);
    }
    if (dmatmul(theMG,level,level,ALL_VECTORS,np->t,A,np->s) != NUM_OK)
      NP_RETURN(1,result[0]);
    if ((*np->Iter->Iter)(np->Iter,level,np->u,np->t,A,result))
      REP_ERR_RETURN(1);
    if (ddot(theMG,level,level,ALL_VECTORS,np->s,np->s,&a0) != NUM_OK)
      NP_RETURN(1,result[0]);
    if (ddot(theMG,level,level,ALL_VECTORS,np->u,np->s,&a1) != NUM_OK)
      NP_RETURN(1,result[0]);

    IFDEBUG(np,1)
    UserWriteF("test vector x\n");
    PrintVector(GRID_ON_LEVEL(theMG,level),np->s,3,3);
    UserWriteF("y = Iter(Ax)\n");
    PrintVector(GRID_ON_LEVEL(theMG,level),np->u,3,3);
    ENDDEBUG

    if (a1 == 0.0)
      a0 = a1 = 1.0;
    else
      a0 /= a1;
    np->damp[level+MAXLEVEL] += a0;
    if (np->display == PCR_FULL_DISPLAY)
      UserWriteF(" test %d: damping factor for %s on level %d = %f\n",
                 i,ENVITEM_NAME(np->Iter),level,a0);
  }
  if (np->n > 0)
    np->damp[level+MAXLEVEL] *= 1.0/np->n;
  else
    np->damp[level+MAXLEVEL] = 1.0;
  if (np->display > PCR_NO_DISPLAY)
    UserWriteF("calibrated damping factor for %s on level %d = %f\n",
               ENVITEM_NAME(np->Iter),level,np->damp[level+MAXLEVEL]);
  FreeVD(theMG,level-1,level,np->s);
  FreeVD(theMG,level,level,np->t);
  FreeVD(theMG,level,level,np->u);

  return (0);
}

static INT Calibrate (NP_ITER *theNP, INT level,
                      VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A,
                      INT *result)
{
  NP_CALIBRATE      *np = (NP_CALIBRATE *) theNP;
  MULTIGRID      *theMG = NP_MG(theNP);

  if ((*np->Iter->Iter)(np->Iter,level,c,b,A,result))
    REP_ERR_RETURN(1);
  if (np->display > PCR_RED_DISPLAY)
    UserWriteF("calibrated damping factor for %s on level %d = %f\n",
               ENVITEM_NAME(np->Iter),level,np->damp[level+MAXLEVEL]);
  if (ABS(np->damp[level+MAXLEVEL]-1.0) < SMALL_D)
    return(0);
  if (dscal(theMG,level,level,ALL_VECTORS,c,(np->damp[level+MAXLEVEL]-1.0))
      != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dmatmul_minus(theMG,level,level,ALL_VECTORS,b,A,c) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (dscal(theMG,level,level,ALL_VECTORS,c,
            np->damp[level+MAXLEVEL]/(np->damp[level+MAXLEVEL]-1.0))
      != NUM_OK)
    NP_RETURN(1,result[0]);

  return (0);
}

static INT CalibratePostProcess (NP_ITER *theNP, INT level,
                                 VECDATA_DESC *x, VECDATA_DESC *b,
                                 MATDATA_DESC *A, INT *result)
{
  NP_CALIBRATE *np = (NP_CALIBRATE *) theNP;

  if (np->Iter->PostProcess != NULL)
    if ((*np->Iter->PostProcess)(np->Iter,level,x,b,A,result))
      REP_ERR_RETURN(1);

  return (0);
}

static INT CalibrateConstruct (NP_BASE *theNP)
{
  NP_ITER *np = (NP_ITER *) theNP;

  theNP->Init = CalibrateInit;
  theNP->Display = CalibrateDisplay;
  theNP->Execute = NPIterExecute;
  np->PreProcess = CalibratePreProcess;
  np->Iter = Calibrate;
  np->PostProcess = CalibratePostProcess;

  return(0);
}

/****************************************************************************/
/*D
   mi - numproc for multiple iterations

   DESCRIPTION:
   This numproc executes a sequence of iterations

   .vb
   npinit <name>
   .ve

   .  $d~<dummy> - correction vector

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT MIInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_MI *np = (NP_MI *) theNP;
  INT i;
  char option[16],buffer[256];

  if (ReadArgvINT("n",&(np->n),argc,argv)) np->n = 0;
  if (np->n<1 || np->n>MI_MAXITER) REP_ERR_RETURN(NP_NOT_ACTIVE);
  for (i=0; i<np->n; i++)
  {
    sprintf(option,"i%d",i);
    if (ReadArgvChar(option,buffer,argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
    np->si[i]=(NP_ITER *)GetNumProcByName(theNP->mg,buffer,ITER_CLASS_NAME);
    if (np->si[i]==NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  }

  return (NPIterInit(&np->iter,argc,argv));
}

static INT MIDisplay (NP_BASE *theNP)
{
  INT i;
  char option[16];
  NP_MI *np = (NP_MI *) theNP;

  NPIterDisplay(&np->iter);

  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",(int)np->n);
  for (i=0; i<np->n; i++)
  {
    sprintf(option,"i%d",i);
    if (np->si[i] != NULL)
      UserWriteF(DISPLAY_NP_FORMAT_SS,option,ENVITEM_NAME(np->si[i]));
    else
      UserWriteF(DISPLAY_NP_FORMAT_SS,option,"---");
  }

  return (0);
}

static INT MIPreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  INT i,j;
  NP_MI *np = (NP_MI *) theNP;

  for (i=0; i<np->n; i++)
  {
    for (j=0; j<i; j++)
      if (np->si[j]==np->si[i]) break;
    if (j!=i) continue;
    if (np->si[i]->PreProcess!=NULL)
      if ((*np->si[i]->PreProcess)(np->si[i],level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);
  }

  return (0);
}

static INT MI (NP_ITER *theNP, INT level, VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  INT i;
  NP_MI *np = (NP_MI *) theNP;

  if (AllocVDFromVD(NP_MG(theNP),level,level,c,&np->t)) NP_RETURN(1,result[0]);
  if (dset(NP_MG(theNP),level,level,ALL_VECTORS,c,0.0)!= NUM_OK) NP_RETURN(1,result[0]);
  for (i=0; i<np->n; i++)
  {
    if ((*np->si[i]->Iter)(np->si[i],level,np->t,b,A,result)) NP_RETURN(1,result[0]);
    if (dadd(NP_MG(theNP),level,level,ALL_VECTORS,c,np->t) != NUM_OK) NP_RETURN(1,result[0]);
  }
  if (FreeVD(NP_MG(theNP),level,level,np->t)) NP_RETURN(1,result[0]);

  return (0);
}

static INT MIPostProcess (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  INT i,j;
  NP_MI *np = (NP_MI *) theNP;

  for (i=0; i<np->n; i++)
  {
    for (j=0; j<i; j++)
      if (np->si[j]==np->si[i]) break;
    if (j!=i) continue;
    if (np->si[i]->PostProcess!=NULL)
      if ((*np->si[i]->PostProcess)(np->si[i],level,x,b,A,result)) REP_ERR_RETURN(1);
  }

  return (0);
}

static INT MIConstruct (NP_BASE *theNP)
{
  NP_ITER *np = (NP_ITER *) theNP;

  theNP->Init = MIInit;
  theNP->Display = MIDisplay;
  theNP->Execute = NPIterExecute;
  np->PreProcess = MIPreProcess;
  np->Iter = MI;
  np->PostProcess = MIPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   sp - numproc for investigation  of  smoothing property

   DESCRIPTION:
   This numproc iterates x <-- A*S^nu, where A is the system matrix and S an iterator

   .vb
   npinit <name>
   .ve

   .  $d~<dummy> - correction vector

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT SPInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_SP *np = (NP_SP *) theNP;
  char buffer[256];
  INT ret;

  /* base class initialization */
  ret=NPIterInit(&np->iter,argc,argv);

  /* initialization for 'executable' */
  if (ReadArgvDOUBLE("tol",&(np->tol),argc,argv)) { np->tol=-1.0; ret=NP_ACTIVE; }
  if (ReadArgvINT("n",&(np->maxiter),argc,argv)) { np->maxiter=-1; ret=NP_ACTIVE; }
  if (np->maxiter<1) { np->maxiter=-1; ret=NP_ACTIVE; }

  /* initialization for 'active' */
  if (ReadArgvINT("nu",&(np->nu),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (np->nu<0) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvChar("I",buffer,argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->S=(NP_ITER *)GetNumProcByName(theNP->mg,buffer,ITER_CLASS_NAME);
  if (ReadArgvChar("mode",buffer,argc,argv))
  {
    PrintErrorMessage('E',"SPInit","specify mode");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }
  np->mode=SP_NOTINIT;
  if (strcmp(buffer,"oms")==0) np->mode=SP_OMS;
  if (strcmp(buffer,"a")==0) np->mode=SP_A;
  if (strcmp(buffer,"o")==0) np->mode=SP_O;
  if (np->mode==SP_NOTINIT)
  {
    PrintErrorMessage('E',"SPInit","wrong mode specification");
    REP_ERR_RETURN (NP_NOT_ACTIVE);
  }

  return (ret);
}

static INT SPDisplay (NP_BASE *theNP)
{
  NP_SP *np = (NP_SP *) theNP;

  NPIterDisplay(&np->iter);

  if (np->tol>=0.0) UserWriteF(DISPLAY_NP_FORMAT_SF,"tol",np->tol);
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"tol","---");
  if (np->maxiter>=0) UserWriteF(DISPLAY_NP_FORMAT_SI,"n",np->maxiter);
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"n","---");

  UserWriteF(DISPLAY_NP_FORMAT_SI,"nu",(int)np->nu);
  if (np->S != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"I",ENVITEM_NAME(np->S));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"I","---");
  if (np->mode==SP_NOTINIT) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","---");
  if (np->mode==SP_OMS) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","oms");
  if (np->mode==SP_A) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","a");
  if (np->mode==SP_O) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","o");

  return (0);
}

static INT SPPreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_SP *np = (NP_SP *) theNP;

  if (np->S->PreProcess!=NULL)
    if ((*np->S->PreProcess)(np->S,level,x,b,A,baselevel,result)) REP_ERR_RETURN(1);

  return (0);
}

static INT SP (NP_ITER *theNP, INT level, VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SP *np = (NP_SP *) theNP;
  INT i;

  if (AllocVDFromVD(NP_MG(theNP),level,level,c,&np->t)) NP_RETURN(1,result[0]);
  for (i=0; i<np->nu; i++)
  {
    if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,np->t,A,b)!=NUM_OK) NP_RETURN(1,result[0]);
    if ((*np->S->Iter)(np->S,level,c,np->t,A,result)) NP_RETURN(1,result[0]);
    if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,b,-1.0,c)!=NUM_OK) NP_RETURN(1,result[0]);
  }
  if (np->mode==SP_OMS)
  {
    if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,np->t,A,b)!=NUM_OK) NP_RETURN(1,result[0]);
    if ((*np->S->Iter)(np->S,level,c,np->t,A,result)) NP_RETURN(1,result[0]);
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,b,c)!= NUM_OK) REP_ERR_RETURN (1);
  }
  else if (np->mode==SP_A)
  {
    if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,np->t,A,b)!=NUM_OK) NP_RETURN(1,result[0]);
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,b,np->t)!= NUM_OK) REP_ERR_RETURN (1);
  }
  else if (np->mode==SP_O)
  {}
  else REP_ERR_RETURN (1);

  if (FreeVD(NP_MG(theNP),level,level,np->t)) NP_RETURN(1,result[0]);
  if (dset(NP_MG(theNP),level,level,ALL_VECTORS,c,0.0)!= NUM_OK) NP_RETURN(1,result[0]);

  return (0);
}

static INT SPPostProcess (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SP *np = (NP_SP *) theNP;

  if (np->S->PostProcess!=NULL)
    if ((*np->S->PostProcess)(np->S,level,x,b,A,result)) REP_ERR_RETURN(1);

  return (0);
}

INT SPExecute (NP_BASE *theNP, INT argc , char **argv)
{
  INT i,level,baselevel,result;
  NP_ITER *npiter=(NP_ITER *)theNP;
  NP_SP *np = (NP_SP *) theNP;
  DOUBLE sp,real,imag,norm_old,norm,s11,s12,s21,s22,a,b;
  char text[DISPLAY_WIDTH+4];
  VECDATA_DESC *x1,*x2,*s1,*s2;

  /* prepare */
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'@',"\n"); UserWriteF("\n%s",text);
  level = CURRENTLEVEL(theNP->mg);
  if (SPPreProcess(npiter,level,npiter->c,npiter->b,npiter->A,&baselevel,&result)) REP_ERR_RETURN (1);
  if (dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,npiter->b,&norm)!=NUM_OK) REP_ERR_RETURN (1);
  if (norm==0.0) REP_ERR_RETURN (1);
  if (dscal(NP_MG(theNP),level,level,ALL_VECTORS,npiter->b,1.0/norm)) REP_ERR_RETURN (1);

  /* allocate vectors and initialize */
  x1=NULL;
  if (AllocVDFromVD(NP_MG(theNP),level,level,npiter->c,&x1)) REP_ERR_RETURN (1);
  if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,x1,npiter->b)!= NUM_OK) REP_ERR_RETURN (1);
  x2=NULL;
  if (AllocVDFromVD(NP_MG(theNP),level,level,npiter->c,&x2)) REP_ERR_RETURN (1);
  if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,x2,x1)!= NUM_OK) REP_ERR_RETURN (1);
  if (SP(npiter,level,npiter->c,x2,npiter->A,&result)) REP_ERR_RETURN (1);
  if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,x1,x2,&sp)!=NUM_OK) REP_ERR_RETURN (1);
  if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,x2,-sp,x1)!=NUM_OK) REP_ERR_RETURN (1);
  if (dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,x2,&norm)!=NUM_OK) REP_ERR_RETURN (1);
  if (norm==0.0) REP_ERR_RETURN (1);
  if (dscal(NP_MG(theNP),level,level,ALL_VECTORS,x2,1.0/norm)) REP_ERR_RETURN (1);
  s1=NULL;
  if (AllocVDFromVD(NP_MG(theNP),level,level,npiter->c,&s1)) REP_ERR_RETURN (1);
  s2=NULL;
  if (AllocVDFromVD(NP_MG(theNP),level,level,npiter->c,&s2)) REP_ERR_RETURN (1);

  /* go */
  real=imag=0.0;
  for (i=1; i<=np->maxiter; i++)
  {
    /* iterate 2-dim-subspace */
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,s1,x1)!= NUM_OK) REP_ERR_RETURN (1);
    if (SP(npiter,level,npiter->c,s1,npiter->A,&result)) REP_ERR_RETURN (1);
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,s2,x2)!= NUM_OK) REP_ERR_RETURN (1);
    if (SP(npiter,level,npiter->c,s2,npiter->A,&result)) REP_ERR_RETURN (1);

    /* solve reduced eigen-value problem */
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,x1,s1,&s11)!=NUM_OK) REP_ERR_RETURN (1);
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,x1,s2,&s12)!=NUM_OK) REP_ERR_RETURN (1);
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,x2,s1,&s21)!=NUM_OK) REP_ERR_RETURN (1);
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,x2,s2,&s22)!=NUM_OK) REP_ERR_RETURN (1);

    norm_old=sqrt(real*real+imag*imag);
    a=0.5*(s11+s22);
    b=0.25*(s11-s22)*(s11-s22)+s12*s21;
    if (b>0.0)
    {
      if (a>=0.0) real=a+sqrt(b);
      else real=a-sqrt(b);
      imag=0.0;
    }
    else
    {
      real=a;
      imag=sqrt(-b);
    }
    norm=sqrt(real*real+imag*imag);
    UserWriteF(" %-3d  %c: %-12.7e   %-12.7e   %-12.7e\n",i,'e',real,imag,norm);
    if (ABS(norm-norm_old)<np->tol*norm) { i++; break; }

    /* update subspace */
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,x1,s1)!= NUM_OK) REP_ERR_RETURN (1);
    if (dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,x1,&norm)!=NUM_OK) REP_ERR_RETURN (1);
    if (norm==0.0) REP_ERR_RETURN (1);
    if (dscal(NP_MG(theNP),level,level,ALL_VECTORS,x1,1.0/norm)) REP_ERR_RETURN (1);
    if (dcopy(NP_MG(theNP),level,level,ALL_VECTORS,x2,s2)!= NUM_OK) REP_ERR_RETURN (1);
    if (ddot(NP_MG(theNP),level,level,ALL_VECTORS,x1,x2,&sp)!=NUM_OK) REP_ERR_RETURN (1);
    if (daxpy(NP_MG(theNP),level,level,ALL_VECTORS,x2,-sp,x1)!=NUM_OK) REP_ERR_RETURN (1);
    if (dnrm2(NP_MG(theNP),level,level,ALL_VECTORS,x2,&norm)!=NUM_OK) REP_ERR_RETURN (1);
    if (norm==0.0) REP_ERR_RETURN (1);
    if (dscal(NP_MG(theNP),level,level,ALL_VECTORS,x2,1.0/norm)) REP_ERR_RETURN (1);
  }

  if (FreeVD(NP_MG(theNP),level,level,x1)) REP_ERR_RETURN (1);
  if (FreeVD(NP_MG(theNP),level,level,x2)) REP_ERR_RETURN (1);
  if (FreeVD(NP_MG(theNP),level,level,s1)) REP_ERR_RETURN (1);
  if (FreeVD(NP_MG(theNP),level,level,s2)) REP_ERR_RETURN (1);

  norm=sqrt(real*real+imag*imag);
  UserWriteF("\n");
  UserWriteF(" %-3d  %c: %-12.7e   %-12.7e   %-12.7e\n\n",i-1,'r',real,imag,norm);

  if (SPPostProcess(npiter,level,npiter->c,npiter->b,npiter->A,&result)) return (1);

  return(0);
}

static INT SPConstruct (NP_BASE *theNP)
{
  NP_ITER *np = (NP_ITER *) theNP;

  theNP->Init = SPInit;
  theNP->Display = SPDisplay;
  theNP->Execute = SPExecute;
  np->PreProcess = SPPreProcess;
  np->Iter = SP;
  np->PostProcess = SPPostProcess;

  return(0);
}

/****************************************************************************/
/*D
   im - calculates iteration matrix of an iteration, usefull only for small matrices

   DESCRIPTION:
   This numproc calculates iteration matrix of an iteration and writes it to file, scilab format

   .vb
   npinit <name>
   .ve

   .  $d~<dummy> - correction vector

   SEE ALSO:
   ls
   D*/
/****************************************************************************/

static INT IMInit (NP_BASE *theNP, INT argc, char **argv)
{
  NP_IM *np=(NP_IM*)theNP;
  char buffer[256];
  INT n,vtype,ncomp;

  np->A = ReadArgvMatDesc(theNP->mg,"A",argc,argv);
  if (np->A==NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->x = ReadArgvVecDesc(theNP->mg,"x",argc,argv);
  if (np->x==NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  for (vtype=0,ncomp=0; vtype<NVECTYPES; vtype++)
  {
    n=VD_NCMPS_IN_TYPE(np->x,vtype);
    if (n==0) continue;
    np->xc=VD_CMP_OF_TYPE(np->x,vtype,0);
    ncomp+=n;
  }
  if (ncomp!=1) REP_ERR_RETURN(NP_NOT_ACTIVE);
  np->S=NULL;
  if (ReadArgvChar("I",buffer,argc,argv)==0)
    np->S=(NP_ITER *)GetNumProcByName(theNP->mg,buffer,ITER_CLASS_NAME);
  if (ReadArgvINT("s",&(np->skip),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvChar ("file",np->file,argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);

  return(NP_EXECUTABLE);
}

static INT IMDisplay (NP_BASE *theNP)
{
  NP_IM *np = (NP_IM*)theNP;

  if (np->S != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"I",ENVITEM_NAME(np->S));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"I","---");
  if (np->A != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"A","---");
  if (np->x != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"x","---");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"s",(int)np->skip);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"file",np->file);

  return (0);
}

INT IMExecute (NP_BASE *theNP, INT argc , char **argv)
{
  INT nv,vtype,bc,n,ncomp,i,j,level,baselevel,result,MarkKey;
  NP_IM *np = (NP_IM*)theNP;
  VECDATA_DESC *b,*Ax;
  VECTOR *v,*w;
  GRID *g;
  HEAP *Heap;
  DOUBLE *Mat;
  FILE *file;

  level = CURRENTLEVEL(theNP->mg);
  g = GRID_ON_LEVEL(theNP->mg,level);
  Heap=MGHEAP(theNP->mg);
  MarkTmpMem(Heap,&MarkKey);
  for (v=FIRSTVECTOR(g),i=0; v!=NULL; v=SUCCVC(v)) if (np->skip*VECSKIP(v)==0) i++;
  nv=i;
  Mat=(DOUBLE*)GetTmpMem(Heap,sizeof(DOUBLE)*nv*nv,MarkKey);
  b=NULL; if (AllocVDFromVD(NP_MG(theNP),0,level,np->x,&b)) REP_ERR_RETURN (1);
  Ax=NULL; if (AllocVDFromVD(NP_MG(theNP),0,level,np->x,&Ax)) REP_ERR_RETURN (1);
  for (vtype=0,ncomp=0; vtype<NVECTYPES; vtype++)
  {
    n=VD_NCMPS_IN_TYPE(b,vtype);
    if (n==0) continue;
    bc=VD_CMP_OF_TYPE(b,vtype,0);
    ncomp+=n;
  }
  if (ncomp!=1) REP_ERR_RETURN(1);
  if (np->S!=NULL) if (np->S->PreProcess!=NULL) if ((*np->S->PreProcess)(np->S,level,np->x,b,np->A,&baselevel,&result)) REP_ERR_RETURN(1);
  if (dset(NP_MG(theNP),level,level,ALL_VECTORS,np->x,0.0)!= NUM_OK) REP_ERR_RETURN(1);
  if (np->S!=NULL)
    for (v=FIRSTVECTOR(g),j=0; v!=NULL; v=SUCCVC(v))
    {
      if (np->skip*VECSKIP(v)!=0) continue;
      VVALUE(v,np->xc)=1.0;
      if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,Ax,np->A,np->x)!=NUM_OK) REP_ERR_RETURN(1);
      if ((*np->S->Iter)(np->S,level,b,Ax,np->A,&result)) REP_ERR_RETURN(1);
      for (w=FIRSTVECTOR(g),i=0; w!=NULL; w=SUCCVC(w))
        if (np->skip*VECSKIP(w)==0)
        {
          if (w==v) Mat[nv*i+j]=1.0-VVALUE(w,bc);
          else Mat[nv*i+j]=   -VVALUE(w,bc);
          i++;
        }
      VVALUE(v,np->xc)=0.0;
      j++;
    }
  else
    for (v=FIRSTVECTOR(g),j=0; v!=NULL; v=SUCCVC(v))
    {
      if (np->skip*VECSKIP(v)!=0) continue;
      VVALUE(v,np->xc)=1.0;
      if (dmatmul(NP_MG(theNP),level,level,ALL_VECTORS,b,np->A,np->x)!=NUM_OK) REP_ERR_RETURN(1);
      for (w=FIRSTVECTOR(g),i=0; w!=NULL; w=SUCCVC(w))
        if (np->skip*VECSKIP(w)==0)
        {
          Mat[nv*i+j]=VVALUE(w,bc);
          i++;
        }
      VVALUE(v,np->xc)=0.0;
      j++;
    }
  if (j!=nv) REP_ERR_RETURN(1);
  if (np->S!=NULL) if (np->S->PostProcess!=NULL) if ((*np->S->PostProcess)(np->S,level,np->x,b,np->A,&result)) REP_ERR_RETURN(1);
  if (FreeVD(NP_MG(theNP),0,level,b)) REP_ERR_RETURN (1);
  if (FreeVD(NP_MG(theNP),0,level,Ax)) REP_ERR_RETURN (1);
  file=fopen(np->file,"w"); if (file==NULL) REP_ERR_RETURN(1);
  for (i=0; i<nv; i++)
  {
    for (j=0; j<nv; j++)
      fprintf(file,"%e ",Mat[nv*i+j]);
    fprintf(file,"\n");
  }
  fclose(file);
  ReleaseTmpMem(Heap,MarkKey);

  return(0);
}

static INT IMConstruct (NP_BASE *theNP)
{
  theNP->Init = IMInit;
  theNP->Display = IMDisplay;
  theNP->Execute = IMExecute;

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

INT NS_DIM_PREFIX InitIter ()
{
  INT i;

  if (MakeStruct(":iter")) REP_ERR_RETURN(__LINE__);

  strcpy(LU_reg[REG_ALWAYS],      "always");
  strcpy(LU_reg[REG_NEVER],       "never");
  strcpy(LU_reg[REG_IF_SING],     "ifsing");

  if (CreateClass(ITER_CLASS_NAME ".jac",sizeof(NP_SMOOTHER),
                  JacobiConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".gs",sizeof(NP_SMOOTHER),GSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".bcgss",sizeof(NP_BCGSSMOOTHER),
                  BCGSSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sgs",sizeof(NP_SGS),SGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".pgs",sizeof(NP_PGS),PGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".block",sizeof(NP_BLOCK),BLOCKConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ts",sizeof(NP_TS),TSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ap",sizeof(NP_TS),APConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ii",sizeof(NP_II),IIConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".bhr",sizeof(NP_TS),BHRConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sor",sizeof(NP_SMOOTHER),
                  SORConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ssor",sizeof(NP_SSOR),SSORConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sbgs",sizeof(NP_SBGS),SBGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".gbgs",sizeof(NP_SBGS),GBGSConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ilu",sizeof(NP_ILU),ILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".bd",sizeof(NP_SMOOTHER),BDConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".filu",sizeof(NP_ILU),FILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".thilu",sizeof(NP_THILU),
                  THILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".spilu",sizeof(NP_ILU),SPILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".spblilu",sizeof(NP_SPBLILU),SPBLILUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ic",sizeof(NP_ILU),ICConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ff",sizeof(NP_FF),FFConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".lu",sizeof(NP_LU),LUConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".lmgc",sizeof(NP_LMGC),LmgcConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".addmgc",sizeof(NP_LMGC),
                  AddmgcConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ex",sizeof(NP_EX),EXConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".exprj",sizeof(NP_EXPRJ),
                  EXPRJConstruct))
    REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".calibrate",sizeof(NP_CALIBRATE),CalibrateConstruct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".mi",sizeof(NP_MI),MIConstruct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".sp",sizeof(NP_SP),SPConstruct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".im",sizeof(NP_SP),IMConstruct)) REP_ERR_RETURN (__LINE__);

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  return (0);
}
