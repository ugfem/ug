// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        <math.h>
#include        "../main/defs.h"
#include        "../main/structs.h"

/* Perform Rayleigh Quotient Iteration */

void RQI(A, u, n, r1, r2, v, w, x, y, work, tol, initshift, evalest,
         vwsqrt, orthlist)
struct vtx_data **A;                    /* matrix/graph being analyzed */
double *u;                              /* eigenvector being refined */
double *r1, *r2, *v, *w, *x, *y, *work; /* work space for symmlq */
int n;                                  /* number of rows/columns in matrix */
double tol;                             /* error tolerance in eigenpair */
double initshift;                       /* initial shift */
double *evalest;                        /* returned eigenvalue */
double *vwsqrt;                         /* square roots of vertex weights */
struct orthlink *orthlist;              /* lower evecs to orthogonalize against */
{
  extern Heap   *heap;       /* pointer to heap of multigrid */
  extern double *MEM_OK;     /* variable for memory overflow exeception */
  extern int DEBUG_EVECS;               /* debug flag for eigen computation */
  extern int WARNING_EVECS;             /* warning flag for eigen computation */
  int rqisteps;                         /* # rqi rqisteps */
  double res;                           /* convergence quant for rqi */
  double last_res;                      /* res on previous rqi step */
  double macheps;                       /* machine precision calculated by symmlq */
  double normxlim;                      /* a stopping criteria for symmlq */
  double normx;                         /* norm of the solution vector */
  int symmlqitns;                       /* # symmlq itns */
  int inv_it_steps;                     /* intial steps of inverse iteration */
  long itnmin;                          /* symmlq input */
  double shift, rtol;                   /* symmlq input */
  long precon, goodb, nout;             /* symmlq input */
  long checka, intlim;                  /* symmlq input */
  double anorm, acond;                  /* symmlq output */
  double rnorm, ynorm;                  /* symmlq output */
  long istop, itn;                      /* symmlq output */
  long long_n;                          /* copy of n for passing to symmlq */
  int warning;                          /* warning on possible misconvergence */
  double factor;                        /* ratio between previous res and new tol */
  double minfactor;                     /* minimum acceptable value of factor */
  double dot(), norm();
  int symmlq_();
  void splarax(), scadd(), vecscale(), doubleout();

  /* Initialize RQI loop */
  splarax(y, A, n, u, vwsqrt, r1);
  if (!MEM_OK) return;
  shift = dot(u,1,n,y);
  scadd(y,1,n,-shift,u);
  res = norm(y,1,n);                    /* eigen-residual */
  rqisteps = 0;                         /* a counter */
  symmlqitns = 0;                       /* a counter */

  /* Set invariant symmlq parameters */
  precon = FALSE;       /* FALSE until we figure out a good way */
  goodb = TRUE;         /* should be TRUE for this application */
  nout = 0;             /* set to 0 for no Symmlq output; 6 for lots */
  checka = FALSE;       /* if don't know by now, too bad */
  intlim = n;           /* set to enforce a maximum number of Symmlq itns */
  itnmin = 0;           /* set to enforce a minimum number of Symmlq itns */
  long_n = n;           /* type change for alint */

  if (DEBUG_EVECS > 0) {
    {char buf[150]; sprintf(buf,"Using RQI/Symmlq refinement on graph with %d vertices.\n",n);UserWrite(buf);}
  }
  if (DEBUG_EVECS > 1) {
    {char buf[150]; sprintf(buf,"  step      lambda est.            Ares          Symmlq its.   istop  factor\n");UserWrite(buf);}
    {char buf[150]; sprintf(buf,"    0");UserWrite(buf);}
    doubleout(shift,1);
    doubleout(res,1);
    {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
  }

  /* Perform RQI */
  inv_it_steps = 2;
  warning = FALSE;
  factor = 10;
  minfactor = factor/2;
  while(res >= tol) {
    if (res/tol < 1.2) {
      factor = max (factor/2,minfactor);
    }
    rtol = res/factor;

    /* exit Symmlq if iterate is this large */
    normxlim = 1.0/rtol;

    if (rqisteps < inv_it_steps) {
      shift = initshift;
    }

    symmlq_(&long_n,&u[1],&r1[1],&r2[1],&v[1],&w[1],&x[1],&y[1],
            work, &checka,&goodb,&precon,&shift,&nout,
            &intlim,&rtol,&istop,&itn,&anorm,&acond,
            &rnorm,&ynorm,(double *)A,vwsqrt,(double *)orthlist,
            &macheps,&normxlim,&itnmin);
    if (!MEM_OK) return;
    symmlqitns += itn;
    normx = norm(x,1,n);
    vecscale(u,1,n,1.0/normx,x);
    splarax(y, A, n, u, vwsqrt, r1);
    if (!MEM_OK) return;
    shift = dot(u,1,n,y);
    scadd(y,1,n,-shift,u);
    last_res = res;
    res = norm(y,1,n);
    if (res > last_res) {warning = TRUE;}
    rqisteps++;
    if (DEBUG_EVECS > 1) {
      {char buf[150]; sprintf(buf,"   %2d",rqisteps);UserWrite(buf);}
      doubleout(shift,1);
      doubleout(res,1);
      {char buf[150]; sprintf(buf,"     %3d",itn);UserWrite(buf);}
      {char buf[150]; sprintf(buf,"          %d",istop);UserWrite(buf);}
      {char buf[150]; sprintf(buf,"      %g\n",factor);UserWrite(buf);}
    }
  }
  *evalest = shift;

  if (WARNING_EVECS > 0 && warning) {
    {char buf[150]; sprintf(buf,"WARNING: Residual convergence not monotonic; RQI may have misconverged.\n");UserWrite(buf);}
  }

  if (DEBUG_EVECS > 0) {
    {char buf[150]; sprintf(buf,"Eval ");UserWrite(buf);}
    doubleout(*evalest,1);
    {char buf[150]; sprintf(buf,"   RQI steps %d,  Symmlq iterations. %d\n\n",rqisteps,symmlqitns);UserWrite(buf);}
  }
}
