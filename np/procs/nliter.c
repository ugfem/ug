// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      nliter.c                                                          */
/*                                                                          */
/* Purpose:   nonlinear iteration num procs                                 */
/*                                                                          */
/* Author:    Gabriele Beddies                                              */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   30.07.97 begin, ug version 3.8                                */
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

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "compiler.h"
#include "devices.h"
#include "gm.h"
#include "misc.h"
#include "algebra.h"
#include "assert.h"
#include "evm.h"
#include "general.h"
#include "ugstruct.h"
#include "ugblas.h"
#include "pcr.h"
#include "np.h"

#include "nls.h"
#include "assemble.h"

#include "fas.h"
#include "nliter.h"

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

#define DAMP_FACTOR     1.0

#define BUFFER_SIZE                 MAX(256,2*DISPLAY_WIDTH+4)

#define OPTIONLEN           32
#define OPTIONLENSTR        "31"
#define VALUELEN            64
#define VALUELENSTR         "63"

/* macros to define VEC_SCALAR, VECDATA_DESC and MATDATA_DESC components */
#define DEFINE_VS_CMPS(a)                               register DOUBLE a ## 0,a ## 1,a ## 2
#define DEFINE_VD_CMPS(x)                               register INT x ## 0,x ## 1,x ## 2
#define DEFINE_MD_CMPS(m)                               register INT m ## 00,m ## 01,m ## 02,m ## 10,m ## 11,m ## 12,m ## 20,m ## 21,m ## 22

#define SET_YCMP_2(y,v,tp,cp)                   {cp=VD_CMPPTR_OF_TYPE(v,tp); y ## 0 = (cp)[0]; y ## 1 = (cp)[1];}
#define SET_MCMP_22(m,M,rt,ct,cp)               {cp = MD_MCMPPTR_OF_RT_CT(M,rt,ct); \
                                                 m ## 00 = (cp)[0]; m ## 01 = (cp)[1]; \
                                                 m ## 10 = (cp)[2]; m ## 11 = (cp)[3];}
#define SET_CMPS_22(y,v,m,M,rt,ct,cp)   SET_MCMP_22(m,M,rt,ct,cp); SET_YCMP_2(y,v,ct,cp);

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* class definition                                                         */
/*                                                                          */
/****************************************************************************/

struct np_nl_smoother {

  NP_NL_ITER iter;

  VEC_SCALAR damp;
  VECDATA_DESC *c;                       /* correction                  */
  MATDATA_DESC *L;                           /* temporary matrix                */

    #ifdef ModelP
  INT cons_mode;
    #endif

  INT (*Step)
    (struct np_nl_smoother *,                /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    VECDATA_DESC *,                              /* correction vector               */
    MATDATA_DESC *,                              /* matrix                          */
    MATDATA_DESC *,                              /* temporary matrix                */
    INT *);                                      /* result                          */
};
typedef struct np_nl_smoother NP_NL_SMOOTHER;

typedef struct
{
  NP_NL_SMOOTHER smoother;

  VECDATA_DESC *t;

  INT displayMode;                       /* for PCR                         */
  INT niter;                             /* number of iterations            */

} NP_NLGS;

INT l_nlgs (NP_NLGS *nlgs, NP_NL_ASSEMBLE *ass, GRID *grid, const DOUBLE *damp,
            VECDATA_DESC *x, VECDATA_DESC *c, MATDATA_DESC *M, VECDATA_DESC *d);

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* general purpose text buffer */
static char buffer[BUFFER_SIZE];

/* for daxpy */
static DOUBLE Factor_One[MAX_VEC_COMP];

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

static INT cy0,cy1,cy2;
static INT m00,m01,m02,m10,m11,m12,m20,m21,m22;
static DOUBLE s0,s1,s2;
static DOUBLE a0,a1,a2;

static INT MySetYComp (const SHORT *cmp, SHORT n)
{
  switch (n)
  {
  case 1 : cy0 = cmp[0];                                              return (0);
  case 2 : cy0 = cmp[0]; cy1 = cmp[1];                        return (0);
  case 3 : cy0 = cmp[0]; cy1 = cmp[1]; cy2 = cmp[2]; return (0);
  }

  return (1);
}

static INT MySetMComp (const SHORT *cmp, SHORT kind)
{
  switch (kind)
  {
  case R1C1 : m00 = cmp[0];                                                         return (0);
  case R1C2 : m00 = cmp[0]; m01 = cmp[1];                           return (0);
  case R1C3 : m00 = cmp[0]; m01 = cmp[1]; m02 = cmp[2]; return (0);

  case R2C1 : m00 = cmp[0];
    m10 = cmp[1];                                                         return (0);
  case R2C2 : m00 = cmp[0]; m01 = cmp[1];
    m10 = cmp[2]; m11 = cmp[3];                           return (0);
  case R2C3 : m00 = cmp[0]; m01 = cmp[1]; m02 = cmp[2];
    m10 = cmp[3]; m11 = cmp[4]; m12 = cmp[5]; return (0);

  case R3C1 : m00 = cmp[0];
    m10 = cmp[1];
    m20 = cmp[2];                                                         return (0);
  case R3C2 : m00 = cmp[0]; m01 = cmp[1];
    m10 = cmp[2]; m11 = cmp[3];
    m20 = cmp[4]; m21 = cmp[5];                           return (0);
  case R3C3 : m00 = cmp[0]; m01 = cmp[1]; m02 = cmp[2];
    m10 = cmp[3]; m11 = cmp[4]; m12 = cmp[5];
    m20 = cmp[6]; m21 = cmp[7]; m22 = cmp[8]; return (0);
  }

  return (1);
}

static INT MySolveSmallBlock (SHORT n, const SHORT *scomp, DOUBLE *sol, const SHORT *mcomp, const DOUBLE *mat,
                              const DOUBLE *rhs)
{
  DOUBLE BlockMat[MAX_SINGLE_MAT_COMP],BlockSol[MAX_SINGLE_VEC_COMP], det;
  DOUBLE aux,M3div0,M6div0;
  register DOUBLE dinv,piv,sum;
  register i,j,k;

  if (n>=MAX_SINGLE_VEC_COMP)
    return (1);

  switch (n)
  {
  case 1 :
    sol[scomp[0]] = rhs[0] / mat[mcomp[0]];
    return (NUM_OK);

  case 2 :
    det = mat[mcomp[0]]*mat[mcomp[3]] - mat[mcomp[1]]*mat[mcomp[2]];
    if (det==0.0) return (1);
    det = 1.0/det;
    sol[scomp[0]] = (rhs[0]*mat[mcomp[3]]-rhs[1]*mat[mcomp[1]])*det;
    sol[scomp[1]] = (rhs[1]*mat[mcomp[0]]-rhs[0]*mat[mcomp[2]])*det;
    return (NUM_OK);

  case 3 :
    M3div0 = mat[mcomp[3]]/mat[mcomp[0]];
    M6div0 = mat[mcomp[6]]/mat[mcomp[0]];
    aux = (mat[mcomp[7]]-M6div0*mat[mcomp[1]]) / (mat[mcomp[4]]-M3div0*mat[mcomp[1]]);

    sol[scomp[2]] = (rhs[2] - M6div0*rhs[0] - aux*(rhs[1]-M3div0*rhs[0]))
                    / (mat[mcomp[8]]-M6div0*mat[mcomp[2]]-aux*(mat[mcomp[5]]-M3div0*mat[mcomp[2]]));
    sol[scomp[1]] = (rhs[1] - mat[mcomp[3]]/mat[mcomp[0]]*rhs[0] - (mat[mcomp[5]]-M3div0*mat[mcomp[2]])*sol[scomp[2]])
                    / (mat[mcomp[4]]-M3div0*mat[mcomp[1]]);
    sol[scomp[0]] = (rhs[0] - mat[mcomp[1]] * sol[scomp[1]] - mat[mcomp[2]] * sol[scomp[2]])
                    /  mat[mcomp[0]];
    return (NUM_OK);

  default :
    /* copy matrix */
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        BlockMat[i*n+j] = mat[mcomp[i*n+j]];

    /* lr factorize mat */
    for (i=0; i<n; i++)
    {
      dinv = BlockMat[i*n+i];
      if (ABS(dinv)<SMALL_D)
        return (NUM_SMALL_DIAG);
      dinv = BlockMat[i*n+i] = 1.0/dinv;

      for (j=i+1; j<n; j++)
      {
        piv = (BlockMat[j*n+i] *= dinv);
        for (k=i+1; k<n; k++)
          BlockMat[j*n+k] -= BlockMat[i*n+k] * piv;
      }
    }

    /* solve */
    for (i=0; i<n; i++)
    {
      for (sum=rhs[i], j=0; j<i; j++)
        sum -= BlockMat[i*n+j] * BlockSol[j];
      BlockSol[i] = sum;                                /* Lii = 1 */
    }
    for (i=n-1; i>=0; i--)
    {
      for (sum=BlockSol[i], j=i+1; j<n; j++)
        sum -= BlockMat[i*n+j] * BlockSol[j];
      BlockSol[i] = sum * BlockMat[i*n+i];                              /* Uii = Inv(Mii) */
    }

    /* copy BlockSol to sol */
    for (i=0; i<n; i++)
      sol[scomp[i]] = BlockSol[i];

    return (NUM_OK);
  }
}

/****************************************************************************/
/****************************************************************************/

INT NPNLIterInit (NP_NL_ITER *np, INT argc , char **argv)
{
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  np->x = ReadArgvVecDesc(np->base.mg,"x",argc,argv);
  np->b = ReadArgvVecDesc(np->base.mg,"r",argc,argv);

  if ((np->A == NULL) || (np->x == NULL) || (np->b == NULL))
    return(NP_ACTIVE);

  /* assemble numproc is required for execution */
  np->Assemble = (NP_NL_ASSEMBLE *)
                 ReadArgvNumProc(np->base.mg,"A",NL_ASSEMBLE_CLASS_NAME,argc,argv);
  if (np->Assemble == NULL) return(NP_ACTIVE);

  return(NP_EXECUTABLE);
}

INT NPNLIterDisplay (NP_NL_ITER *np)
{
  if ((np->A == NULL) && (np->x == NULL) && (np->b == NULL))
    return(0);
  UserWrite("symbolic user data:\n");
  if (np->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));
  if (np->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
  if (np->b != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->b));
  UserWrite("\n");

  return(0);
}

INT NPNLIterExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_NL_ITER *np;
  INT result,bl,level;

  np = (NP_NL_ITER *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  if (np->b == NULL) {
    PrintErrorMessage('E',"NPINLterExecute","no vector b");
    REP_ERR_RETURN (1);
  }
  if (np->x == NULL) {
    PrintErrorMessage('E',"NPINLterExecute","no vector x");
    REP_ERR_RETURN (1);
  }
  if (np->A == NULL) {
    PrintErrorMessage('E',"NPNLIterExecute","no matrix A");
    REP_ERR_RETURN (1);
  }
  if (np->Assemble == NULL) {
    PrintErrorMessage('E',"NPNLIterExecute","no assemble num proc");
    return (1);
  }

  if (ReadArgvOption("i",argc,argv)) {
    if (np->PreProcess == NULL) {
      PrintErrorMessage('E',"NPIterExecute","no PreProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PreProcess)(np,level,np->b,np->x,np->A,&bl,&result)) {
      UserWriteF("NPIterExecute: PreProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("s",argc,argv)) {
    if (np->NLIter == NULL) {
      PrintErrorMessage('E',"NPNLIterExecute","no Iter");
      REP_ERR_RETURN (1);
    }
    if ((*np->NLIter)(np,level,np->b,np->x,np->A,np->Assemble,&result)) {
      UserWriteF("NPIterExecute: Iter failed, error code %d\n", result);
      REP_ERR_RETURN (1);
    }
  }

  if (ReadArgvOption("p",argc,argv)) {
    if (np->PostProcess == NULL) {
      PrintErrorMessage('E',"NPNLIterExecute","no PostProcess");
      REP_ERR_RETURN (1);
    }
    if ((*np->PostProcess)(np,level,np->b,np->x,np->A,&result)) {
      UserWriteF("NPIterExecute: PostProcess failed, error code %d\n",
                 result);
      REP_ERR_RETURN (1);
    }
  }

  return(0);
}

static INT NLSmootherInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_NL_SMOOTHER *np;
  INT i;

  np = (NP_NL_SMOOTHER *) theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i] = DAMP_FACTOR;
  sc_read(np->damp,NP_FMT(np),np->c,"damp",argc,argv);

  np->c = ReadArgvVecDesc(theNP->mg,"c",argc,argv);
  np->L = ReadArgvMatDesc(theNP->mg,"L",argc,argv);
    #ifdef ModelP
  if (ReadArgvOption("M",argc,argv))
    np->cons_mode = MAT_CONS;
  else if (ReadArgvOption("D",argc,argv))
    np->cons_mode = MAT_DIAG_CONS;
  else
    np->cons_mode = MAT_MASTER_CONS;
        #endif
  return (NPNLIterInit(&np->iter,argc,argv));
}

static INT NLSmootherDisplay (NP_NL_SMOOTHER *theNP)
{
  NPNLIterDisplay(&theNP->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(theNP->damp,theNP->iter.b,"damp")) REP_ERR_RETURN (1);
  if (theNP->c != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"c",ENVITEM_NAME(theNP->c));
  if (theNP->L != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(theNP->L));
    #ifdef ModelP
  UserWriteF(DISPLAY_NP_FORMAT_SI,"cons_mode",(int)theNP->cons_mode);
        #endif

  return (0);
}

static INT NLSmoother (NP_NL_ITER *theNP, INT level,
                       VECDATA_DESC *x, VECDATA_DESC *b,MATDATA_DESC *A,
                       NP_NL_ASSEMBLE *ass, INT *result)
{
  NP_NL_SMOOTHER *np;
  GRID *theGrid;

  /* store passed XXXDATA_DESCs */
  NPINL_A(theNP) = A;
  NPINL_x(theNP) = x;
  NPINL_b(theNP) = b;

  np = (NP_NL_SMOOTHER *) theNP;
  theGrid = NP_GRID(theNP,level);
  /* check function pointers in numprocs */
  if (ass->NLNAssembleDefect==NULL)
  {
    UserWrite("NLGS: ass->NLNAssembleDefect not defined\n");
    REP_ERR_RETURN (1);
  }
  if (ass->NLAssembleMatrix==NULL)
  {
    UserWrite("NLGS: ass->NLAssembleMatrix not defined\n");
    REP_ERR_RETURN (1);
  }
  if (ass->NLNAssembleMatrix==NULL)
  {
    UserWrite("NLGS: ass->NLNAssembleMatrix not defined\n");
    REP_ERR_RETURN (1);
  }

  np->iter.Assemble = ass;

  if ((*np->Step)(np,level,x,b,np->c,A,np->L,result))
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

static INT SmootherPostProcess (NP_NL_ITER *theNP, INT level,
                                VECDATA_DESC *x, VECDATA_DESC *b,
                                MATDATA_DESC *A, INT *result)
{
  NP_NL_SMOOTHER *np;

  np = (NP_NL_SMOOTHER *) theNP;
  FreeVD(NP_MG(theNP),level,level,np->c);
  if (np->L != NULL)
    FreeMD(NP_MG(theNP),level,level,np->L);

  return(0);
}

/****************************************************************************/
/*D
   nlgs - numproc for nonlinear Gauss-Seidel smoother

   DESCRIPTION:
   This numproc executes a Gauss-Seidel smoother, using the blas routines
   'l_nlgs'.

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

static INT NLGSPreProcess  (NP_NL_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b,
                            MATDATA_DESC *A, INT *baselevel, INT *result)
{

  NP_NL_SMOOTHER *np =(NP_NL_SMOOTHER *) theNP;
  GRID *theGrid = NP_GRID(theNP,level);
  MULTIGRID *mg = theNP->base.mg;

  if (AllocVDFromVD(mg,level,level,x,&(np->c)))
    NP_RETURN(1,result[0]);
        #ifdef ModelP
  if (AllocMDFromMD(mg,level,level,A,&np->L))
    NP_RETURN(1,result[0]);
  if (dmatcopy(mg,level,level,ALL_VECTORS,np->L,A) != NUM_OK)
    NP_RETURN(1,result[0]);
  if (l_matrix_consistent(theGrid,np->L,np->cons_mode) != NUM_OK)
    NP_RETURN(1,result[0]);
        #endif
  if (l_setindex(theGrid))
    NP_RETURN(1,result[0]);
  *baselevel = level;

  return (0);
}

static INT NLGSStep (NP_NL_SMOOTHER *theNP, INT level,
                     VECDATA_DESC *x, VECDATA_DESC *b, VECDATA_DESC *c,
                     MATDATA_DESC *A, MATDATA_DESC *L, INT *result)
{
  NP_NLGS *nlgs;
  NP_NL_ASSEMBLE *ass;
  MULTIGRID *mg;
  INT i,error;

  nlgs = (NP_NLGS*) theNP;
  ass = nlgs->smoother.iter.Assemble;

  mg = nlgs->smoother.iter.base.mg;

  for (i=0; i<nlgs->niter; i++)
  {
    dmatset(mg,level,level,ALL_VECTORS,A,0.0);
    dset (mg,level,level,ALL_VECTORS,theNP->c,0.0);

    /* iterate */
    if (l_nlgs(nlgs,ass,NP_GRID(theNP,level),theNP->damp,x,theNP->c,A,b) != NUM_OK) NP_RETURN(1,result[0]);
  }

  return (0);
}

INT l_nlgs (NP_NLGS *nlgs, NP_NL_ASSEMBLE *ass, GRID *grid, const DOUBLE *damp,
            VECDATA_DESC *x, VECDATA_DESC *v, MATDATA_DESC *M,
            VECDATA_DESC *d)
{
  VECTOR *vec,*w,*first_vec;
  NODE *theNode;
  MULTIGRID *mg;
  INT level;
  INT rtype,ctype,myindex,error;
  register MATRIX *mat;
  register SHORT vc,dc,mc,xc,mask;
  register SHORT *tmpptr,*mcomp,*wcomp,*dcomp,*xcomp,*vcomp;
  register SHORT i,j,l;
  register SHORT n,nc;
  register DOUBLE sum;
  DEFINE_VD_CMPS(cy);
  DEFINE_MD_CMPS(m);
  DOUBLE r[MAX_SINGLE_VEC_COMP],*wmat;

  mg = nlgs->smoother.iter.base.mg;
  level = GLEVEL(grid);
  first_vec = FIRSTVECTOR(grid);

  L_VLOOP__CLASS(vec,first_vec,ACTIVE_CLASS)
  {
    rtype = VTYPE(vec);

    /* get node */
    theNode = (NODE*)VOBJECT(vec);

    n     = VD_NCMPS_IN_TYPE(v,rtype);
    if (n == 0) continue;
    dcomp = VD_CMPPTR_OF_TYPE(d,rtype);
    xcomp = VD_CMPPTR_OF_TYPE(x,rtype);
    vcomp = VD_CMPPTR_OF_TYPE(v,rtype);
    myindex = VINDEX(vec);

    /* Jacobi matrix */
    if ((*ass->NLNAssembleMatrix)(ass,level,level,theNode,x,d,v,M,&error)) {
      error = __LINE__;
      REP_ERR_RETURN(error);
    }

    /* get defect */
    for (i=0; i<n; i++)
      r[i] = VVALUE(vec,dcomp[i]);

    /* rhs
       for (ctype=0; ctype<=NVECTYPES; ctype++)
       if (MD_ROWS_IN_RT_CT(M,rtype,ctype)>0)
            {
                    SET_CMPS_22(cy,v,m,M,rtype,ctype,tmpptr);
                s0 = s1 = 0.0;
                    for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
                            if (((VTYPE(w=MDEST(mat))==ctype) && (VCLASS(w)>=ACTIVE_CLASS)) && (myindex>VINDEX(w)))
                                    MATMUL_22(s,mat,m,w,cy);
                                    r[0] -= s0;
                                    r[1] -= s1;
       }*/

    /* solve */
    if (MySolveSmallBlock(n,VD_CMPPTR_OF_TYPE(v,rtype),VVALPTR(vec),
                          MD_MCMPPTR_OF_RT_CT(M,rtype,rtype),
                          MVALPTR(VSTART(vec)),r)!=0)
      return (__LINE__);

    /* damp */
    for (i=0; i<n; i++)
      VVALUE(vec,vcomp[i]) *= damp[i];

    /* update solution */
    for (i=0; i<n; i++)
      VVALUE(vec,xcomp[i]) -= VVALUE(vec,vcomp[i]);
  }

  return (0);
}

static INT NLGS_Init (NP_BASE *base, INT argc , char **argv)
{
  NP_NLGS *nlgs;
  INT j,nopt,Dopt;
  char option[OPTIONLEN],value[VALUELEN];

  nlgs = (NP_NLGS*) base;

  /* set configuration parameters */
  if (ReadArgvINT("n",&(nlgs->niter),argc,argv))
    nlgs->niter = 1;
  if ((nlgs->niter<0)||(nlgs->niter>10)) {
    PrintErrorMessage('E',"NLGS_Init","n <= 10");
    REP_ERR_RETURN(NP_NOT_ACTIVE);
  }

  /* set display option */
  nlgs->displayMode = ReadArgvDisplay(argc,argv);

  /*return (NPNLIterInit(&(nlgs->smoother.iter),argc,argv));*/
  return (NLSmootherInit(base,argc,argv));
}

static INT NLGS_Display (NP_BASE *theNumProc)
{
  NP_NLGS *nlgs;

  nlgs = (NP_NLGS*) theNumProc;

  /* general nl_iter display */
  NLSmootherDisplay(&nlgs->smoother);

  if (nlgs->displayMode == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (nlgs->displayMode == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (nlgs->displayMode == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");

  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",(int)nlgs->niter);

  return (0);
}

static INT NLGSConstruct (NP_BASE *theNP)
{
  NP_NL_SMOOTHER *np;

  theNP->Init = NLGS_Init;
  theNP->Display = NLGS_Display;
  theNP->Execute = NPNLIterExecute;

  np = (NP_NL_SMOOTHER *) theNP;
  np->iter.PreProcess = NLGSPreProcess;
  np->iter.NLIter = NLSmoother;
  np->iter.PostProcess = SmootherPostProcess;
  np->Step = NLGSStep;

  return(0);
}

/****************************************************************************/
/*D
   InitNLIter - Initialization of nonlinear smoother numprocs

   SYNOPSIS:
   INT InitNLIter (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes the numprocs 'nlgs'.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT InitNLIter (void)
{
  INT i;

  /* init static variables */
  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  if (CreateClass (NL_ITER_CLASS_NAME ".nlgs",
                   sizeof(NP_NLGS), NLGSConstruct))
    REP_ERR_RETURN (__LINE__);

  return (0);
}
