// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  test.c														*/
/*																			*/
/* Purpose:   example application for algebraic multigrid                                       */
/*			  This file shows how to allocate a system of linear eqs. in sp	*/
/*            and fill it with values										*/
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
/* History:   31 JAN 1996 Begin												*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#undef MWCW

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#ifdef MWCW
#include <console.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

/* our linear system */
static AMG_MATRIX *A0;
static AMG_VECTOR *x0,*b0;


/****************************************************************************/
/*																			*/
/* Discretization, set up of system	2D										*/
/*																			*/
/****************************************************************************/

static double u0 (double x, double y)
{
  return(x*x+y*y);
}

static double kx (double x, double y)
{
  double r;

  r = 0.0;
  return(1.0E-6);
}

static double ky (double x, double y)
{
  double r;

  r = 0.0;
  return(1.0E-6);
}

static double rx (double x, double y)
{
  return(1.0);
}

static double ry (double x, double y)
{
  return(0.0);
}

static double f (double x, double y)
{
  return(0.0);
}

static int SetupSystem2D (int N) /* N is the number of points in one direction */
{
  int i,j;
  int me,nonzeros;
  double x,y,h,w;
  double x_e,x_w,y_n,y_s;
  double kx_nw,kx_ne,kx_sw,kx_se,kx_n,kx_e,kx_s,kx_w;
  double ky_nw,ky_ne,ky_sw,ky_se,ky_n,ky_e,ky_s,ky_w;
  double rx_e,rx_w,ry_n,ry_s;

  /* allocate the matrix and the vectors */
  A0 = AMG_NewMatrix(N*N,1,N*N*5,1,"fine grid A0");
  if (A0==NULL) {AMG_Print("no memory for A0\n"); return(1);}
  x0 = AMG_NewVector(N*N,1,"fine grid x0");
  if (x0==NULL) {AMG_Print("no memory for x0\n"); return(1);}
  b0 = AMG_NewVector(N*N,1,"fine grid b0");
  if (b0==NULL) {AMG_Print("no memory for b0\n"); return(1);}

  AMG_randomize(x0);

  /* mesh size */
  h = 1.0/((double)(N-1));

  /* set up the matrix by running through the grid in lexicographic order */
  for (j=0; j<N; j++)
    for (i=0; i<N; i++)
    {
      /* center point */
      x = i*h;
      y = j*h;
      me = j*N+i;

      /* half points */
      x_w = x-0.5*h;
      x_e = x+0.5*h;
      y_n = y+0.5*h;
      y_s = y-0.5*h;

      /* permeability at center of elements */
      kx_nw = kx(x_w,y_n);
      kx_ne = kx(x_e,y_n);
      kx_sw = kx(x_w,y_s);
      kx_se = kx(x_e,y_s);
      ky_nw = ky(x_w,y_n);
      ky_ne = ky(x_e,y_n);
      ky_sw = ky(x_w,y_s);
      ky_se = ky(x_e,y_s);

      /* harmonic mean of permeabilities */
      kx_w = 2.0 / ( (1.0/kx_nw)+(1.0/kx_sw) );
      kx_e = 2.0 / ( (1.0/kx_ne)+(1.0/kx_se) );
      kx_n = 2.0 / ( (1.0/kx_nw)+(1.0/kx_ne) );
      kx_s = 2.0 / ( (1.0/kx_sw)+(1.0/kx_se) );
      ky_w = 2.0 / ( (1.0/ky_nw)+(1.0/ky_sw) );
      ky_e = 2.0 / ( (1.0/ky_ne)+(1.0/ky_se) );
      ky_n = 2.0 / ( (1.0/ky_nw)+(1.0/ky_ne) );
      ky_s = 2.0 / ( (1.0/ky_sw)+(1.0/ky_se) );

      /* velocities */
      rx_e = rx(x_e,y);
      rx_w = rx(x_w,y);
      ry_n = ry(x,y_n);
      ry_s = ry(x,y_s);

      /* count the number of nonzeros in row me and allocate row */
      nonzeros=1;
      if (i>0) nonzeros++;
      if (i<N-1) nonzeros++;
      if (j>0) nonzeros++;
      if (j<N-1) nonzeros++;
      if (AMG_SetRowLength(A0,me,nonzeros)!=AMG_OK) return(1);

      if ( (i>0) && (i<N-1) && (j>0) && (j<N-1) )
      {                         /* interior node */
                                /* diagonal coefficient and right hand side */
        w = kx_w + kx_e + ky_n + ky_s + h*(rx_e+ry_n);
        if (AMG_InsertValues(A0,me,me,&w)<0) return(1);
        AMG_VECTOR_ENTRY(b0,me,0) = f(x,y)*h*h;
        /*AMG_VECTOR_ENTRY(x0,me,0) = u0(x,y);*/

        /* west coefficient */
        w = -kx_w - h*rx_w;
        if (AMG_InsertValues(A0,me,(j  )*N+(i-1),&w)<0) return(1);

        /* south coefficient */
        w = -ky_s - h*ry_s;
        if (AMG_InsertValues(A0,me,(j-1)*N+(i  ),&w)<0) return(1);

        /* east coefficient */
        w = -kx_e;
        if (AMG_InsertValues(A0,me,(j  )*N+(i+1),&w)<0) return(1);

        /* north coefficient */
        w = -ky_n;
        if (AMG_InsertValues(A0,me,(j+1)*N+(i  ),&w)<0) return(1);
      }
      else
      {
        /* Dirichlet node */
        w = 1.0;
        if (AMG_InsertValues(A0,me,me,&w)<0) return(1);
        AMG_VECTOR_ENTRY(b0,me,0) = 0.0;
        AMG_VECTOR_ENTRY(x0,me,0) = 0.0;

        w = 0.0;
        if (i>0)
          if (AMG_InsertValues(A0,me,(j  )*N+(i-1),&w)<0) return(1);
        if (j>0)
          if (AMG_InsertValues(A0,me,(j-1)*N+(i  ),&w)<0) return(1);
        if (i<N-1)
          if (AMG_InsertValues(A0,me,(j  )*N+(i+1),&w)<0) return(1);
        if (j<N-1)
          if (AMG_InsertValues(A0,me,(j+1)*N+(i  ),&w)<0) return(1);
      }
    }

  return(AMG_OK);
}



/****************************************************************************/
/*																			*/
/* Discretization, set up of system	2D										*/
/*																			*/
/****************************************************************************/

static double u0_3d (double x, double y, double z)
{
  return(x*x+y*y);
}

static double kx_3d (double x, double y, double z)
{
  double r;

  r = 0.0;
  return(1.0);
}

static double ky_3d (double x, double y, double z)
{
  double r;

  r = 0.0;
  return(1.0);
}

static double kz_3d (double x, double y, double z)
{
  double r;

  r = 0.0;
  return(1.0E-6);
}

/* !!! rx,ry,rz positive !!! */
static double rx_3d (double x, double y, double z)
{
  return(0.0);
}

static double ry_3d (double x, double y, double z)
{
  return(0.0);
}

static double rz_3d (double x, double y, double z)
{
  return(0.0);
}

static double f_3d (double x, double y, double z)
{
  return(0.0);
}

static int SetupSystem3D (int N) /* N is the number of points in one direction */
{
  int i,j,k;
  int me,nonzeros;
  double x,y,z,h,w;
  double k_x,k_y,k_z,r_x,r_y,r_z;

  /* allocate the matrix and the vectors */
  A0 = AMG_NewMatrix(N*N*N,1,N*N*N*7,1,"fine grid A0");
  if (A0==NULL) {AMG_Print("no memory for A0\n"); return(1);}
  x0 = AMG_NewVector(N*N*N,1,"fine grid x0");
  if (x0==NULL) {AMG_Print("no memory for x0\n"); return(1);}
  b0 = AMG_NewVector(N*N*N,1,"fine grid b0");
  if (b0==NULL) {AMG_Print("no memory for b0\n"); return(1);}

  AMG_randomize(x0);

  /* mesh size */
  h = 1.0/((double)(N-1));

  /* set up the matrix by running through the grid in lexicographic order */
  for (k=0; k<N; k++)
    for (j=0; j<N; j++)
      for (i=0; i<N; i++)
      {
        /* center point */
        x = i*h;
        y = j*h;
        z = k*h;
        me = k*N*N+j*N+i;

        k_x=kx_3d(x,y,z);
        k_y=ky_3d(x,y,z);
        k_z=kz_3d(x,y,z);

        r_x = rx_3d(x,y,z);
        r_y = ry_3d(x,y,z);
        r_z = rz_3d(x,y,z);

        /* count the number of nonzeros in row me and allocate row */
        nonzeros=1;
        if (i>0) nonzeros++;
        if (i<N-1) nonzeros++;
        if (j>0) nonzeros++;
        if (j<N-1) nonzeros++;
        if (k>0) nonzeros++;
        if (k<N-1) nonzeros++;
        if (AMG_SetRowLength(A0,me,nonzeros)!=AMG_OK) return(1);

        if ( (i>0) && (i<N-1) && (j>0) && (j<N-1) && (k>0) && (k<N-1) )
        {                               /* interior node */
                                        /* diagonal coefficient and right hand side */
          w = 2.0*(k_x+k_y+k_z) + h*(r_x+r_y+r_z);
          if (AMG_InsertValues(A0,me,me,&w)<0) return(1);
          AMG_VECTOR_ENTRY(b0,me,0) = f_3d(x,y,z)*h*h;
          /*AMG_VECTOR_ENTRY(x0,me,0) = u0_3d(x,y,z);*/

          /* west coefficient */
          w = -k_x - h*r_x;
          if (AMG_InsertValues(A0,me,(k  )*N*N+(j  )*N+(i-1),&w)<0) return(1);

          /* east coefficient */
          w = -k_x;
          if (AMG_InsertValues(A0,me,(k  )*N*N+(j  )*N+(i+1),&w)<0) return(1);

          /* south coefficient */
          w = -k_y - h*r_y;
          if (AMG_InsertValues(A0,me,(k  )*N*N+(j-1)*N+(i  ),&w)<0) return(1);

          /* north coefficient */
          w = -k_y;
          if (AMG_InsertValues(A0,me,(k  )*N*N+(j+1)*N+(i  ),&w)<0) return(1);

          /* down coefficient */
          w = -k_z - h*r_z;
          if (AMG_InsertValues(A0,me,(k-1)*N*N+(j  )*N+(i  ),&w)<0) return(1);

          /* up coefficient */
          w = -k_z;
          if (AMG_InsertValues(A0,me,(k+1)*N*N+(j  )*N+(i  ),&w)<0) return(1);
        }
        else
        {
          /* Dirichlet node */
          w = 1.0;
          if (AMG_InsertValues(A0,me,me,&w)<0) return(1);
          AMG_VECTOR_ENTRY(b0,me,0) = 0.0;
          AMG_VECTOR_ENTRY(x0,me,0) = 0.0;

          w = 0.0;
          if (i>0)
            if (AMG_InsertValues(A0,me,(k  )*N*N+(j  )*N+(i-1),&w)<0) return(1);
          if (j>0)
            if (AMG_InsertValues(A0,me,(k  )*N*N+(j-1)*N+(i  ),&w)<0) return(1);
          if (k>0)
            if (AMG_InsertValues(A0,me,(k-1)*N*N+(j  )*N+(i  ),&w)<0) return(1);
          if (i<N-1)
            if (AMG_InsertValues(A0,me,(k  )*N*N+(j  )*N+(i+1),&w)<0) return(1);
          if (j<N-1)
            if (AMG_InsertValues(A0,me,(k  )*N*N+(j+1)*N+(i  ),&w)<0) return(1);
          if (k<N-1)
            if (AMG_InsertValues(A0,me,(k+1)*N*N+(j  )*N+(i  ),&w)<0) return(1);
        }
      }

  return(AMG_OK);
}




int main (int argc, char *argv[])
{
  int N;
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;

  /* read size parameter from command line */
#ifdef MWCW
  argc = ccommand(&argv);
#endif
  if (argc<=1) {printf("usage: test <unknowns per dimension>\n"); return(1);}
  sscanf(argv[1],"%ld",&N);
  if ((N<0)||(N>10000000)) return(1);

  AMG_Print("Setting up system\n");
  if (SetupSystem2D(N))
    AMG_Print("Error in setup\n");
  /*AMG_PrintMatrix(A0,"A initially");*/

  /* coarsen context */
  cc.alpha = 0.33333333;
  cc.beta = 1.0E-3;
  cc.mincluster=4;
  cc.maxcluster=6;
  cc.maxdistance=2;
  cc.maxconnectivity=15;
  cc.verbose=1;
  cc.depthtarget=20;
  cc.coarsentarget=10;
  cc.coarsenrate=1.2;
  cc.major=-1;

  /* solver context */
  sc.verbose=1;
  sc.solver=AMG_BCGS;
  sc.preconditioner=AMG_MGC;
  sc.maxit=800;
  sc.red_factor=1.0E-8;
  sc.dnorm_min=1.0E-15;
  sc.coarse_smoother=AMG_SSOR;
  sc.coarse_maxit=200;
  sc.coarse_red_factor=1.0E-3;
  sc.n1=3;
  sc.n2=3;
  sc.gamma=1;
  sc.omega_p[0]=1.8;
  sc.smoother=AMG_DJAC;
  sc.omega[0]=0.5;

  AMG_Build(&sc,&cc,A0);
  AMG_Solve(x0,b0);

  return(0);
}
