// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugeblas.c                                                     */
/*                                                                          */
/* Purpose:   source for extended blas routines                             */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            IWR/Technische Simulation                                     */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: klaus.johannsen@iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   19.07.02 begin                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include "general.h"
#include "ugblas.h"

USING_UG_NAMESPACES

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* blas level 1																*/
/*																			*/
/****************************************************************************/

INT NS_PREFIX deset (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, DOUBLE a)
{
  INT i,ret,level;

  ret=dset(mg,fl,tl,mode,x->vd,a); if (ret!=NUM_OK) return ret;
  for (level=fl; level<=tl; level++)
    for (i=0; i<x->n; i++) EVDD_E(x,level,i)=a;

  return NUM_OK;
}

INT NS_PREFIX deadd (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, const EVECDATA_DESC *y)
{
  INT i,ret,level;

  ret=dadd(mg,fl,tl,mode,x->vd,y->vd); if (ret!=NUM_OK) return ret;
  for (level=fl; level<=tl; level++)
    for (i=0; i<x->n; i++) EVDD_E(x,level,i)+=EVDD_E(y,level,i);

  return NUM_OK;
}

INT NS_PREFIX decopy (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, const EVECDATA_DESC *y)
{
  INT i,ret,level;

  ret=dcopy(mg,fl,tl,mode,x->vd,y->vd); if (ret!=NUM_OK) return ret;
  for (level=fl; level<=tl; level++)
    for (i=0; i<x->n; i++) EVDD_E(x,level,i)=EVDD_E(y,level,i);

  return NUM_OK;
}

INT NS_PREFIX dedotx (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, const EVECDATA_DESC *y, EVEC_SCALAR a)
{
  INT i,n,ret;

  if (x->n!=y->n) return NUM_ERROR;
  ret=ddotx(mg,fl,tl,mode,x->vd,y->vd,a); if (ret!=NUM_OK) return ret;
  n=x->n;
  for (i=0; i<x->n; i++) a[i+n]=EVDD_E(x,tl,i)*EVDD_E(y,tl,i);

  return NUM_OK;
}

INT NS_PREFIX dedotw (MULTIGRID *mg, INT fl, INT tl, INT mode, const EVECDATA_DESC *x, const EVECDATA_DESC *y, const EVEC_SCALAR w, DOUBLE *a)
{
  INT i,n,ret;

  if (x->n!=y->n) return NUM_ERROR;
  ret=ddotw(mg,fl,tl,mode,x->vd,y->vd,w,a); if (ret!=NUM_OK) return ret;
  n=VD_NCOMP(x->vd);
  for (i=0; i<x->n; i++) (*a)+=w[n+i]*EVDD_E(x,tl,i)*EVDD_E(y,tl,i);

  return NUM_OK;
}

INT NS_PREFIX denrm2x (MULTIGRID *mg, INT fl, INT tl, INT mode, const EVECDATA_DESC *x, EVEC_SCALAR a)
{
  INT i,n,ret;

  ret=dnrm2x(mg,fl,tl,mode,x->vd,a); if (ret!=NUM_OK) return ret;
  n=VD_NCOMP(x->vd);
  for (i=0; i<x->n; i++) a[i+n]=fabs(EVDD_E(x,tl,i));

  return NUM_OK;
}

INT NS_PREFIX descal (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, DOUBLE a)
{
  INT ret,level,i;

  ret=dscal(mg,fl,tl,mode,x->vd,a); if (ret!=NUM_OK) return ret;
  for (level=fl; level<=tl; level++)
    for (i=0; i<x->n; i++) EVDD_E(x,level,i)*=a;

  return NUM_OK;
}

INT NS_PREFIX deaxpy (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, DOUBLE a, const EVECDATA_DESC *y)
{
  INT i,level,ret;

  ret=daxpy(mg,fl,tl,mode,x->vd,a,y->vd); if (ret!=NUM_OK) return ret;
  for (level=fl; level<=tl; level++)
    for (i=0; i<x->n; i++) EVDD_E(x,level,i)+=a*EVDD_E(y,level,i);

  return NUM_OK;
}

/****************************************************************************/
/*																			*/
/* blas level 2																*/
/*																			*/
/****************************************************************************/

INT NS_PREFIX dematmul_minus (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, const EMATDATA_DESC *M, const EVECDATA_DESC *y)
{
  INT i,j,n,ret,level;
  DOUBLE a;

  if (x->n!=M->n || M->n!=y->n) return NUM_ERROR;
  n=M->n;
  ret=dmatmul_minus(mg,fl,tl,mode,x->vd,M->mm,y->vd); if (ret!=NUM_OK) return ret;
  for (i=0; i<n; i++)
  {
    ret=daxpy(mg,fl,tl,mode,x->vd,-EVDD_E(y,tl,i),M->me[i]); if (ret!=NUM_OK) return ret;
    ret=ddot(mg,fl,tl,mode,y->vd,M->em[i],&a); if (ret!=NUM_OK) return ret;EVDD_E(x,tl,i)-=a;
    for (level=fl; level<=tl; level++)
      for (j=0; j<n; j++)
        EVDD_E(x,tl,i)-=EMDD_EE(M,level,i*n+j)*EVDD_E(y,tl,j);
  }
  return NUM_OK;
}

INT NS_PREFIX dematmul (MULTIGRID *mg, INT fl, INT tl, INT mode, EVECDATA_DESC *x, const EMATDATA_DESC *M, const EVECDATA_DESC *y)
{
  INT i,j,n,ret,level;
  DOUBLE a;

  if (x->n!=M->n || M->n!=y->n) return NUM_ERROR;
  n=M->n;
  ret=dmatmul(mg,fl,tl,mode,x->vd,M->mm,y->vd); if (ret!=NUM_OK) return ret;
  for (i=0; i<n; i++)
  {
    ret=daxpy(mg,fl,tl,mode,x->vd,EVDD_E(y,tl,i),M->me[i]); if (ret!=NUM_OK) return ret;
    ret=ddot(mg,fl,tl,mode,y->vd,M->em[i],&a); if (ret!=NUM_OK) return ret;EVDD_E(x,tl,i)=a;
    for (level=fl; level<=tl; level++)
      for (j=0; j<n; j++)
        EVDD_E(x,tl,i)+=EMDD_EE(M,level,i*n+j)*EVDD_E(y,tl,j);
  }
  return NUM_OK;
}

INT NS_PREFIX dematset (MULTIGRID *mg, INT fl, INT tl, INT mode, EMATDATA_DESC *M, DOUBLE a)
{
  INT i,ret,level;

  ret=dmatset(mg,fl,tl,mode,M->mm,a); if (ret!=NUM_OK) return ret;
  for (i=0; i<M->n; i++)
  {
    ret=dset(mg,fl,tl,mode,M->me[i],a); if (ret!=NUM_OK) return ret;
    ret=dset(mg,fl,tl,mode,M->em[i],a); if (ret!=NUM_OK) return ret;
  }
  for (level=fl; level<=tl; level++)
    for (i=0; i<M->n*M->n; i++) EMDD_EE(M,level,i)=a;

  return NUM_OK;
}
