// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_solve.c													*/
/*																			*/
/* Purpose:   solvers for AMG												*/
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
/* History:   05 FEB 1996 Begin												*/
/*			  02 OKT 1997 redesign											*/
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "amg_header.h"
#include "amg_low.h"
#include "amg_sp.h"
#include "amg_blas.h"
#include "amg_iter.h"
#include "amg_coarsen.h"
#include "amg_solve.h"

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

/* Iteration procedure, this may be multigrid, sor, jacobi, ilu, ...		*/
/* The aim is to solve Mx=b with the following constraints:					*/
/* M is an auxiliary matrix containing e.g. an ILU decomposition of A		*/
/* x contains an iterate on entry and holds the new iterate on exit			*/
/* b is the right hand side THAT MAY NOT BE CHANGED BY THE PROCEDURE !		*/
/* d contains the defect d=b-Ax on entry, arbitrary on exit					*/
typedef int (*IterProcPtr)(AMG_SolverContext *sc, int level, int depth,        \
                           AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],                    \
                           AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_VECTOR *x[AMG_MAX_LEVELS],                   \
                           AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

static IterProcPtr smoother,coarse_smoother;
static IterProcPtr preconditioner;

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


static int pc_restrict (AMG_GRAPH *g, AMG_VECTOR *fine, AMG_VECTOR *coarse)
{
  int *ca=AMG_GRAPH_CA(g);       /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i;

  /* consistency */
  if (AMG_VECTOR_B(fine)!=AMG_VECTOR_B(coarse)) return(AMG_FATAL);
  if (AMG_VECTOR_N(fine)!=AMG_GRAPH_N(g)) return(AMG_FATAL);

  if (b==1)
  {
    for (i=0; i<nc; i++) c[i] = 0.0;
    for (i=0; i<nf; i++) c[ca[i]] += f[i];
  }
  else
  {
    for (i=0; i<nc; i++) c[i] = 0.0;
    for (i=0; i<nf; i++) c[ca[i/b]*b+(i%b)] += f[i];
  }
  return(AMG_OK);
}


static int pc_prolongate (AMG_GRAPH *g, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp)
{
  int *ca=AMG_GRAPH_CA(g);       /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i;
  double om;

  /* consistency */
  if (AMG_VECTOR_B(fine)!=AMG_VECTOR_B(coarse)) return(AMG_FATAL);
  if (AMG_VECTOR_N(fine)!=AMG_GRAPH_N(g)) return(AMG_FATAL);

  if (b==1)
  {
    om=damp[0];
    for (i=0; i<nf; i++) f[i] += c[ca[i]]*om;
  }
  else
  {
    for (i=0; i<nf; i++) f[i] += c[ca[i/b]*b+(i%b)]*damp[i%b];
  }
  return(AMG_OK);
}


static int pc_prolongate_auto (AMG_GRAPH *g, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp)
{
  int *ca=AMG_GRAPH_CA(g);       /* cluster mapping information */
  float *da=AMG_GRAPH_DA(g);       /* auto damping factor */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i;
  double w1,w2;

  /* consistency */
  if (AMG_VECTOR_B(fine)!=AMG_VECTOR_B(coarse)) return(AMG_FATAL);
  if (AMG_VECTOR_N(fine)!=AMG_GRAPH_N(g)) return(AMG_FATAL);

  if (b==1)
  {
    w1=2.0-damp[0]; w2=damp[0]-1;
    for (i=0; i<nf; i++) f[i] += c[ca[i]]*(w1+w2*da[i]);
  }
  else
  {
    for (i=0; i<nf; i++) f[i] += c[ca[i/b]*b+(i%b)]*damp[i%b];
  }
  return(AMG_OK);
}


static int ssor (AMG_SolverContext *sc, int k, int depth,                               \
                 AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],            \
                 AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_VECTOR *x[AMG_MAX_LEVELS],           \
                 AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  /* defect d=b-Ax is valid on entry */
  AMG_sorf(A[k],d[k],d[k],sc->omega);       /* compute correction, overwrite d */
  AMG_daxpy(x[k],1.0,d[k]);                             /* update solution x */
  AMG_dcopy(d[k],b[k]);                                 /* recompute d=b-Ax */
  AMG_dmatminus(d[k],A[k],x[k]);
  AMG_sorb(A[k],d[k],d[k],sc->omega);       /* backward step */
  AMG_daxpy(x[k],1.0,d[k]);
  return(AMG_OK);
}

static int sor (AMG_SolverContext *sc, int k, int depth,                                \
                AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],            \
                AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_VECTOR *x[AMG_MAX_LEVELS],           \
                AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  /* defect d=b-Ax is valid on entry */
  AMG_sorf(A[k],d[k],d[k],sc->omega);       /* compute correction, overwrite d */
  AMG_daxpy(x[k],1.0,d[k]);                             /* update solution x */
  return(AMG_OK);
}

static int jac (AMG_SolverContext *sc, int k, int depth,                                \
                AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],            \
                AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_VECTOR *x[AMG_MAX_LEVELS],           \
                AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  /* defect d=b-Ax is valid on entry */
  AMG_jac(A[k],d[k],d[k],sc->omega);            /* compute correction, overwrite d */
  AMG_daxpy(x[k],1.0,d[k]);                             /* update solution x */
  return(AMG_OK);
}

static int ex (AMG_SolverContext *sc, int k, int depth,                         \
               AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],            \
               AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_VECTOR *x[AMG_MAX_LEVELS],           \
               AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  /* defect d=b-Ax is valid on entry */
  AMG_EXApplyLU(AMG_MATRIX_A(M[k]),AMG_MATRIX_BW(M[k]),AMG_MATRIX_N(M[k]),AMG_VECTOR_X(d[k]));
  AMG_daxpy(x[k],1.0,d[k]);                             /* update solution x */
  return(AMG_OK);
}

static AMG_MATRIX *prepare_ex (AMG_MATRIX *A)
{
  int bw;
  int i,n=A->n,k,start,end;
  int *ra=A->ra, *ja=A->ja;
  double *a=A->a,*lu;
  AMG_MATRIX *new;

  /* compute bandwith */
  bw=0;
  for (i=0; i<n; i++)
  {
    start = ra[i]; end = start+ja[start];
    for (k=start+1; k<end; k++)
      bw=AMG_MAX(bw,AMG_ABS(i-ja[k]));
  }

  /* allocate new matrix */
  new = AMG_NewMatrix(n,1,n*(2*bw+1),AMG_MATRIX_SAS(A),"ex matrix");
  if (new==NULL) return(new);
  lu=new->a;
  AMG_MATRIX_BW(new)=bw;

  /* insert & copy entries */
  for (i=0; i<n*(2*bw+1); i++) lu[i]=0.0;
  for (i=0; i<n; i++)
  {
    start = ra[i]; end = start+ja[start];
    AMG_EX_MAT(lu,bw,i,i) = a[start];
    for (k=start+1; k<end; k++)
      AMG_EX_MAT(lu,bw,i,ja[k]) = a[k];
  }

  /* decompose */
  if (AMG_EXDecomposeMatrix(lu,bw,n)) return(AMG_NULL);

  /* return matrix */
  return(new);
}


static int mgc (AMG_SolverContext *sc, int k, int depth,                                \
                AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],            \
                AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_VECTOR *x[AMG_MAX_LEVELS],           \
                AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  int i;
  int dnorm, dnorm0;

  /* defect d=b-Ax is valid on entry */
  if (k==depth)
  {
    /* coarse grid solve */
    dnorm=dnorm0=sqrt(AMG_ddot(d[k],d[k]));
    for (i=0; i<sc->coarse_maxit; i++)
    {
      coarse_smoother(sc,k,depth,A,G,M,x,b,d);
      AMG_dcopy(d[k],b[k]);
      AMG_dmatminus(d[k],A[k],x[k]);
      dnorm=sqrt(AMG_ddot(d[k],d[k]));
      if (dnorm<dnorm0*sc->coarse_red_factor) break;
      if (dnorm<sc->dnorm_min) break;
    }
    if (i==sc->coarse_maxit)
      AMG_Print("coarse grid solver not converged\n");
  }
  else
  {
    for (i=0; i<sc->n1; i++)
    {
      smoother(sc,k,depth,A,G,M,x,b,d);
      AMG_dcopy(d[k],b[k]);
      AMG_dmatminus(d[k],A[k],x[k]);
    }
    pc_restrict(G[k],d[k],b[k+1]);
    AMG_dcopy(d[k+1],b[k+1]);
    AMG_dset(x[k+1],0.0);
    for (i=0; i<AMG_MIN(sc->gamma,depth-k); i++)
    {
      mgc(sc,k+1,depth,A,G,M,x,b,d);
      if (i+1==AMG_MIN(sc->gamma,depth-k)) break;
      AMG_dcopy(d[k+1],b[k+1]);                   /* because d must cont defect on entry */
      AMG_dmatminus(d[k+1],A[k+1],x[k+1]);
    }
    pc_prolongate_auto(G[k],x[k],x[k+1],sc->omega_p);
    for (i=0; i<sc->n2; i++)
    {
      AMG_dcopy(d[k],b[k]);
      AMG_dmatminus(d[k],A[k],x[k]);
      smoother(sc,k,depth,A,G,M,x,b,d);
    }
  }

  return(AMG_OK);
}

/* global data for solvers */
static AMG_MATRIX *A[AMG_MAX_LEVELS];
static AMG_MATRIX *M[AMG_MAX_LEVELS];
static AMG_GRAPH *G[AMG_MAX_LEVELS];
static AMG_VECTOR *x[AMG_MAX_LEVELS];
static AMG_VECTOR *b[AMG_MAX_LEVELS];
static AMG_VECTOR *d[AMG_MAX_LEVELS];
static AMG_VECTOR *z[AMG_MAX_LEVELS];
static AMG_VECTOR *r[AMG_MAX_LEVELS];
static AMG_VECTOR *q;
static AMG_VECTOR *p[AMG_MAX_LEVELS];
static AMG_VECTOR *w;
static int depth;
static AMG_CoarsenContext *global_cc;
static AMG_SolverContext *global_sc;

/* linear solver */

static int ls_build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A_in)
{
  int k;

  /* set up coarse grid matrices */
  switch (sc->preconditioner)
  {
  case AMG_MGC :
    depth=AMG_BuildHierarchy(cc,A_in,A,G);
    if (depth<0) {
      AMG_Print("Could not set up coarse grid matrices\n");
      return(AMG_FATAL);
    }
    break;
  default :
    depth=0;
    A[0]=A_in;
  }

  /* allocate vectors, we need only 3 on each level */
  d[0] = AMG_NewVector(AMG_MATRIX_N(A[0]),AMG_MATRIX_B(A[0]),"d");
  if (d[0]==NULL) {AMG_Print("no memory for d\n"); return(AMG_FATAL);}
  for (k=1; k<=depth; k++)
  {
    x[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"x");
    if (x[k]==NULL) {AMG_Print("no memory for x\n"); return(AMG_FATAL);}
    b[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"b");
    if (b[k]==NULL) {AMG_Print("no memory for b\n"); return(AMG_FATAL);}
    d[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"d");
    if (d[k]==NULL) {AMG_Print("no memory for d\n"); return(AMG_FATAL);}
  }

  /* preconditioner preprocess */
  switch (sc->preconditioner)
  {
  case AMG_DJAC :
    preconditioner = jac;
    break;
  case AMG_SOR :
    preconditioner = sor;
    break;
  case AMG_SSOR :
    preconditioner = ssor;
    break;
  case AMG_MGC :
    preconditioner = mgc;
    break;
  default :
    AMG_Print("invalid preconditioner\n");
    return(AMG_FATAL);
  }

  /* smoother preprocess */
  for (k=0; k<=depth; k++) M[k]=A[k];       /* default, change only for ILU */
  if (sc->preconditioner==AMG_MGC)
    switch (sc->smoother)
    {
    case AMG_DJAC :
      smoother = jac;
      break;
    case AMG_SOR :
      smoother = sor;
      break;
    case AMG_SSOR :
      smoother = ssor;
      break;
    default :
      AMG_Print("invalid smoother\n");
      return(AMG_FATAL);
    }

  /* smoother preprocess */
  if (sc->preconditioner==AMG_MGC)
    switch (sc->coarse_smoother)
    {
    case AMG_DJAC :
      coarse_smoother = jac;
      break;
    case AMG_SOR :
      coarse_smoother = sor;
      break;
    case AMG_SSOR :
      coarse_smoother = ssor;
      break;
    case AMG_EX :
      M[depth]=prepare_ex(A[depth]);
      if (M[depth]==NULL) {
        AMG_Print("error in prepare_ex\n");
        return(AMG_FATAL);
      }
      coarse_smoother = ex;
      break;
    default :
      AMG_Print("invalid coarse smoother\n");
      return(AMG_FATAL);
    }

  return(AMG_OK);
}


static int ls_solve (AMG_VECTOR *x_in, AMG_VECTOR *b_in)
{
  int i;
  double dnorm, dnorm0, dnormlast;
  char buf[128];
  AMG_SolverContext *sc=global_sc;

  x[0] = x_in;
  b[0] = b_in;

  /* iterate */
  AMG_dcopy(d[0],b[0]);
  AMG_dmatminus(d[0],A[0],x[0]);
  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(d[0],d[0]));
  if (sc->verbose>0)
  {
    sprintf(buf,"%4d %12.4lE \n",0,dnorm);
    AMG_Print(buf);
  }
  for (i=0; i<sc->maxit; i++)
  {
    preconditioner(sc,0,depth,A,G,M,x,b,d);
    AMG_dcopy(d[0],b[0]);
    AMG_dmatminus(d[0],A[0],x[0]);
    dnorm=sqrt(AMG_ddot(d[0],d[0]));
    if (sc->verbose>0)
    {
      sprintf(buf,"%4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
      AMG_Print(buf);
    }
    dnormlast=dnorm;
    if (sc->ex_maxit) continue;
    if (dnorm<dnorm0*sc->red_factor) break;
    if (dnorm<sc->dnorm_min) break;
  }
  if (i==sc->maxit && !sc->ex_maxit)
  {
    AMG_Print("solver not converged\n");
    return(-1);
  }
  return(i+1);
}


/* cg solver */

static int cg_build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A_in)
{
  int k;

  /* set up coarse grid matrices */
  switch (sc->preconditioner)
  {
  case AMG_MGC :
    depth=AMG_BuildHierarchy(cc,A_in,A,G);
    if (depth<0) {
      AMG_Print("Could not set up coarse grid matrices\n");
      return(AMG_FATAL);
    }
    break;
  default :
    depth=0;
    A[0]=A_in;
  }

  /* allocate vectors, we need only 3 on each level */
  z[0] = AMG_NewVector(AMG_MATRIX_N(A[0]),AMG_MATRIX_B(A[0]),"z");
  if (z[0]==NULL) {AMG_Print("no memory for z\n"); return(AMG_FATAL);}
  d[0] = AMG_NewVector(AMG_MATRIX_N(A[0]),AMG_MATRIX_B(A[0]),"d");
  if (d[0]==NULL) {AMG_Print("no memory for d\n"); return(AMG_FATAL);}
  q = AMG_NewVector(AMG_MATRIX_N(A[0]),AMG_MATRIX_B(A[0]),"q");
  if (q==NULL) {AMG_Print("no memory for q\n"); return(AMG_FATAL);}
  for (k=1; k<=depth; k++)
  {
    z[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"z");
    if (z[k]==NULL) {AMG_Print("no memory for z\n"); return(AMG_FATAL);}
    r[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"r");
    if (r[k]==NULL) {AMG_Print("no memory for r\n"); return(AMG_FATAL);}
    d[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"d");
    if (d[k]==NULL) {AMG_Print("no memory for d\n"); return(AMG_FATAL);}
  }

  /* preconditioner preprocess */
  switch (sc->preconditioner)
  {
  case AMG_DJAC :
    preconditioner = jac;
    break;
  case AMG_SOR :
    preconditioner = sor;
    break;
  case AMG_SSOR :
    preconditioner = ssor;
    break;
  case AMG_MGC :
    preconditioner = mgc;
    break;
  default :
    AMG_Print("invalid preconditioner\n");
    return(AMG_FATAL);
  }

  /* smoother preprocess */
  for (k=0; k<=depth; k++) M[k]=A[k];       /* default, change only for ILU */
  if (sc->preconditioner==AMG_MGC)
    switch (sc->smoother)
    {
    case AMG_DJAC :
      smoother = jac;
      break;
    case AMG_SOR :
      smoother = sor;
      break;
    case AMG_SSOR :
      smoother = ssor;
      break;
    default :
      AMG_Print("invalid smoother\n");
      return(AMG_FATAL);
    }

  /* smoother preprocess */
  if (sc->preconditioner==AMG_MGC)
    switch (sc->coarse_smoother)
    {
    case AMG_DJAC :
      coarse_smoother = jac;
      break;
    case AMG_SOR :
      coarse_smoother = sor;
      break;
    case AMG_SSOR :
      coarse_smoother = ssor;
      break;
    case AMG_EX :
      M[depth]=prepare_ex(A[depth]);
      if (M[depth]==NULL) {
        AMG_Print("error in prepare_ex\n");
        return(AMG_FATAL);
      }
      coarse_smoother = ex;
      break;
    default :
      AMG_Print("invalid coarse smoother\n");
      return(AMG_FATAL);
    }

  return(AMG_OK);
}



static int cg_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i;
  double dnorm, dnorm0, dnormlast;
  double rho,rho_last,alpha;
  char buf[128];
  AMG_SolverContext *sc=global_sc;

  /* iterate */
  r[0] = b;       /* overwrite b with residual */
  AMG_dmatminus(r[0],A[0],x);
  rho_last = 1.0;

  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(r[0],r[0]));
  if (sc->verbose>0)
  {
    sprintf(buf,"%4d %12.4lE \n",0,dnorm);
    AMG_Print(buf);
  }
  if (dnorm<sc->dnorm_min) return(0);

  for (i=0; i<sc->maxit; i++)
  {
    AMG_dset(z[0],0.0);
    AMG_dcopy(d[0],r[0]);
    preconditioner(sc,0,depth,A,G,M,z,r,d);
    rho = AMG_ddot(r[0],z[0]);
    if (i>0) AMG_daxpy(z[0],rho/rho_last,q);
    rho_last=rho;
    AMG_dcopy(q,z[0]);
    AMG_dmatmul(d[0],A[0],q);
    alpha=rho/AMG_ddot(q,d[0]);
    AMG_daxpy(x,alpha,q);
    AMG_daxpy(r[0],-alpha,d[0]);
    dnorm=sqrt(AMG_ddot(r[0],r[0]));
    if (sc->verbose>0)
    {
      sprintf(buf,"%4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
      AMG_Print(buf);
    }
    dnormlast=dnorm;
    if (sc->ex_maxit) continue;
    if (dnorm<dnorm0*sc->red_factor) break;
    if (dnorm<sc->dnorm_min) break;
  }
  if (i==sc->maxit && !sc->ex_maxit)
  {
    AMG_Print("solver not converged\n");
    return(-1);
  }
  return(i+1);
}


static int bcgs_build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A_in)
{
  int k;

  /* set up coarse grid matrices */
  switch (sc->preconditioner)
  {
  case AMG_MGC :
    depth=AMG_BuildHierarchy(cc,A_in,A,G);
    if (depth<0) {
      AMG_Print("Could not set up coarse grid matrices\n");
      return(AMG_FATAL);
    }
    break;
  default :
    depth=0;
    A[0]=A_in;
  }

  /* allocate vectors, we need only 3 on each level */
  w = AMG_NewVector(AMG_MATRIX_N(A[0]),AMG_MATRIX_B(A[0]),"w");
  if (w==NULL) {AMG_Print("no memory for w\n"); return(AMG_FATAL);}
  for (k=0; k<=depth; k++)
  {
    z[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"z");
    if (z[k]==NULL) {AMG_Print("no memory for z\n"); return(AMG_FATAL);}
    r[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"r");
    if (r[k]==NULL) {AMG_Print("no memory for r\n"); return(AMG_FATAL);}
    p[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"p");
    if (p[k]==NULL) {AMG_Print("no memory for p\n"); return(AMG_FATAL);}
    d[k] = AMG_NewVector(AMG_MATRIX_N(A[k]),AMG_MATRIX_B(A[k]),"d");
    if (d[k]==NULL) {AMG_Print("no memory for d\n"); return(AMG_FATAL);}
  }

  /* preconditioner preprocess */
  switch (sc->preconditioner)
  {
  case AMG_DJAC :
    preconditioner = jac;
    break;
  case AMG_SOR :
    preconditioner = sor;
    break;
  case AMG_SSOR :
    preconditioner = ssor;
    break;
  case AMG_MGC :
    preconditioner = mgc;
    break;
  default :
    AMG_Print("invalid preconditioner\n");
    return(AMG_FATAL);
  }

  /* smoother preprocess */
  for (k=0; k<=depth; k++) M[k]=A[k];       /* default, change only for ILU */
  if (sc->preconditioner==AMG_MGC)
    switch (sc->smoother)
    {
    case AMG_DJAC :
      smoother = jac;
      break;
    case AMG_SOR :
      smoother = sor;
      break;
    case AMG_SSOR :
      smoother = ssor;
      break;
    default :
      AMG_Print("invalid smoother\n");
      return(AMG_FATAL);
    }

  /* smoother preprocess */
  if (sc->preconditioner==AMG_MGC)
    switch (sc->coarse_smoother)
    {
    case AMG_DJAC :
      coarse_smoother = jac;
      break;
    case AMG_SOR :
      coarse_smoother = sor;
      break;
    case AMG_SSOR :
      coarse_smoother = ssor;
      break;
    case AMG_EX :
      M[depth]=prepare_ex(A[depth]);
      if (M[depth]==NULL) {
        AMG_Print("error in prepare_ex\n");
        return(AMG_FATAL);
      }
      coarse_smoother = ex;
      break;
    default :
      AMG_Print("invalid coarse smoother\n");
      return(AMG_FATAL);
    }

  return(AMG_OK);
}



static int bcgs_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i;
  double dnorm, dnorm0, dnormlast;
  double rho,rho_last,alpha,beta,omega,delta;
  char buf[128];
  AMG_SolverContext *sc=global_sc;

  /* iterate */
  AMG_dcopy(r[0],b);
  AMG_dmatminus(r[0],A[0],x);       /* r[0] is initial residual */
  AMG_dcopy(b,r[0]);            /* b is initial residual, r^tilde in algorithm */

  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(r[0],r[0]));
  if (sc->verbose>0)
  {
    sprintf(buf,"%4d %12.4lE \n",0,dnorm);
    AMG_Print(buf);
  }
  if (dnorm<sc->dnorm_min) return(0);

  for (i=0; i<sc->maxit; i++)
  {
    rho=AMG_ddot(r[0],b);             /* b is r[0]^tilde */
    if (fabs(rho)<1.0E-50)
    {
      AMG_Print("BCGS break down\n");
      return(-1);
    }

    if (i==0)             /* compute p */
      AMG_dcopy(p[0],r[0]);
    else
    {
      beta=(rho/rho_last)*(alpha/omega);
      AMG_dscale(p[0],beta);
      AMG_daxpy(p[0],1.0,r[0]);
      AMG_daxpy(p[0],-beta*omega,w);
    }
    rho_last=rho;

    AMG_dset(z[0],0.0);
    AMG_dcopy(d[0],p[0]);
    preconditioner(sc,0,depth,A,G,M,z,p,d);
    AMG_dmatmul(w,A[0],z[0]);
    delta=AMG_ddot(b,w);
    if (fabs(delta)<1.0E-50) delta=1.0E-50;
    alpha=rho/delta;
    AMG_daxpy(r[0],-alpha,w);
    AMG_daxpy(x,alpha,z[0]);
    dnorm=sqrt(AMG_ddot(r[0],r[0]));
    if ( (dnorm<dnorm0*sc->red_factor || dnorm<sc->dnorm_min) && !sc->ex_maxit )
    {
      if (sc->verbose>0)
      {
        sprintf(buf,"%4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
        AMG_Print(buf);
      }
      break;
    }

    AMG_dset(z[0],0.0);
    AMG_dcopy(d[0],r[0]);
    preconditioner(sc,0,depth,A,G,M,z,r,d);
    AMG_dmatmul(d[0],A[0],z[0]);
    delta=AMG_ddot(d[0],d[0]);
    if (fabs(delta)<1.0E-50) delta=1.0E-50;
    omega=AMG_ddot(d[0],r[0])/delta;
    AMG_daxpy(r[0],-omega,d[0]);
    AMG_daxpy(x,omega,z[0]);
    dnorm=sqrt(AMG_ddot(r[0],r[0]));
    if (sc->verbose>0)
    {
      sprintf(buf,"%4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
      AMG_Print(buf);
    }
    dnormlast=dnorm;
    if (sc->ex_maxit) continue;
    if (dnorm<dnorm0*sc->red_factor) break;
    if (dnorm<sc->dnorm_min) break;
  }
  if (i==sc->maxit && !sc->ex_maxit)
  {
    AMG_Print("solver not converged\n");
    return(-1);
  }
  return(i+1);
}

/****************************************************************************/
/*D
   amg - solve system of linear equations with algebraic multigrid

   SYNOPSIS:
   int AMG_Build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A);
   int AMG_Solve (AMG_VECTOR *x, AMG_VECTOR *b);

   PARAMETERS:

   DESCRIPTION:

   RETURN VALUE:

   D*/
/****************************************************************************/

int AMG_Build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A)
{
  int rv;

  /* store context */
  global_cc=cc; global_sc=sc;

  switch (sc->solver)
  {
  case AMG_CG :
    rv=cg_build(sc,cc,A);
    break;

  case AMG_BCGS :
    rv=bcgs_build(sc,cc,A);
    break;

  case AMG_LS :
    rv=ls_build(sc,cc,A);
    break;

  default :
    AMG_Print("solver not implemented\n");
    return(AMG_FATAL);
  }

  return(rv);
}


int AMG_Solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int rv;

  switch (global_sc->solver)
  {
  case AMG_CG :
    rv=cg_solve(x,b);
    break;

  case AMG_BCGS :
    rv=bcgs_solve(x,b);
    break;

  case AMG_LS :
    rv=ls_solve(x,b);
    break;

  default :
    AMG_Print("solver not implemented\n");
    return(-1);
  }

  return(rv);
}
