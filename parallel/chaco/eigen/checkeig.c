// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

/* Check an eigenpair of A by direct multiplication.  */
double checkeig(err,A,y,n,lambda,vwsqrt,work)
double *err;
struct vtx_data **A;
double *y;
int n;
double lambda;
double *vwsqrt;
double *work;
{
  extern Heap   *heap;       /* pointer to heap of multigrid */
  extern double *MEM_OK;     /* variable for memory overflow exeception */
  double resid;
  double normy;
  double norm();
  void splarax(), scadd();

  splarax(err, A, n, y, vwsqrt, work);
  if (!MEM_OK) return;
  scadd(err,1,n,-lambda,y);
  normy = norm(y,1,n);
  resid = norm(err,1,n)/normy;
  return(resid);
}
