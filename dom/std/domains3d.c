// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      domains3d.c                                                   */
/*                                                                          */
/* Purpose:   domain definitions                                                */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de					                */
/*																			*/
/* History:   Sep 11 1996                                                                               */
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* low modules */
#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "misc.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"

/* dev modules */
#include "devices.h"

/* domain module */
#include "std_domain.h"

/* TODO: this violates the subsystem hierarchy */
#include "scan.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* needed again if scan.h removed: typedef DOUBLE DOUBLE_VECTOR[DIM];*/

#define V3_EUKLIDNORM(A,b)                              (b) = (sqrt((double)((A)[0]*(A)[0]+(A)[1]*(A)[1]+(A)[2]*(A)[2])));

#define V3_SCALE(c,C)                              {(C)[0] = (c)*(C)[0];\
                                                    (C)[1] = (c)*(C)[1];\
                                                    (C)[2] = (c)*(C)[2];}


/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

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

static DOUBLE_VECTOR x_hex[8];

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/****************************************************************************/
/*                             domain definitions                           */
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  define the unit cube                                                    */
/*                                                                          */
/****************************************************************************/

static INT southBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] =   x_hex[0][0]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[1][0]*lambda1*(1.0-lambda2)
              + x_hex[2][0]*lambda1*lambda2
              + x_hex[3][0]*(1.0-lambda1)*lambda2;
  result[1] =   x_hex[0][1]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[1][1]*lambda1*(1.0-lambda2)
              + x_hex[2][1]*lambda1*lambda2
              + x_hex[3][1]*(1.0-lambda1)*lambda2;
  result[2] =   x_hex[0][2]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[1][2]*lambda1*(1.0-lambda2)
              + x_hex[2][2]*lambda1*lambda2
              + x_hex[3][2]*(1.0-lambda1)*lambda2;

  /* return ok */
  return(0);
}

static INT eastBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] =   x_hex[1][0]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[2][0]*lambda1*(1.0-lambda2)
              + x_hex[6][0]*lambda1*lambda2
              + x_hex[5][0]*(1.0-lambda1)*lambda2;
  result[1] =   x_hex[1][1]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[2][1]*lambda1*(1.0-lambda2)
              + x_hex[6][1]*lambda1*lambda2
              + x_hex[5][1]*(1.0-lambda1)*lambda2;
  result[2] =   x_hex[1][2]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[2][2]*lambda1*(1.0-lambda2)
              + x_hex[6][2]*lambda1*lambda2
              + x_hex[5][2]*(1.0-lambda1)*lambda2;

  /* return ok */
  return(0);
}

static INT northBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] =   x_hex[4][0]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[5][0]*lambda1*(1.0-lambda2)
              + x_hex[6][0]*lambda1*lambda2
              + x_hex[7][0]*(1.0-lambda1)*lambda2;
  result[1] =   x_hex[4][1]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[5][1]*lambda1*(1.0-lambda2)
              + x_hex[6][1]*lambda1*lambda2
              + x_hex[7][1]*(1.0-lambda1)*lambda2;
  result[2] =   x_hex[4][2]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[5][2]*lambda1*(1.0-lambda2)
              + x_hex[6][2]*lambda1*lambda2
              + x_hex[7][2]*(1.0-lambda1)*lambda2;

  /* return ok */
  return(0);
}

static INT westBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] =   x_hex[0][0]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[3][0]*lambda1*(1.0-lambda2)
              + x_hex[7][0]*lambda1*lambda2
              + x_hex[4][0]*(1.0-lambda1)*lambda2;
  result[1] =   x_hex[0][1]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[3][1]*lambda1*(1.0-lambda2)
              + x_hex[7][1]*lambda1*lambda2
              + x_hex[4][1]*(1.0-lambda1)*lambda2;
  result[2] =   x_hex[0][2]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[3][2]*lambda1*(1.0-lambda2)
              + x_hex[7][2]*lambda1*lambda2
              + x_hex[4][2]*(1.0-lambda1)*lambda2;

  /* return ok */
  return(0);
}

static INT frontBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] =   x_hex[0][0]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[1][0]*lambda1*(1.0-lambda2)
              + x_hex[5][0]*lambda1*lambda2
              + x_hex[4][0]*(1.0-lambda1)*lambda2;
  result[1] =   x_hex[0][1]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[1][1]*lambda1*(1.0-lambda2)
              + x_hex[5][1]*lambda1*lambda2
              + x_hex[4][1]*(1.0-lambda1)*lambda2;
  result[2] =   x_hex[0][2]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[1][2]*lambda1*(1.0-lambda2)
              + x_hex[5][2]*lambda1*lambda2
              + x_hex[4][2]*(1.0-lambda1)*lambda2;

  /* return ok */
  return(0);
}

static INT backBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 1.0;
  result[2] = lambda2;
  result[0] =   x_hex[3][0]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[2][0]*lambda1*(1.0-lambda2)
              + x_hex[6][0]*lambda1*lambda2
              + x_hex[7][0]*(1.0-lambda1)*lambda2;
  result[1] =   x_hex[3][1]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[2][1]*lambda1*(1.0-lambda2)
              + x_hex[6][1]*lambda1*lambda2
              + x_hex[7][1]*(1.0-lambda1)*lambda2;
  result[2] =   x_hex[3][2]*(1.0-lambda1)*(1.0-lambda2)
              + x_hex[2][2]*lambda1*(1.0-lambda2)
              + x_hex[6][2]*lambda1*lambda2
              + x_hex[7][2]*(1.0-lambda1)*lambda2;

  /* return ok */
  return(0);
}

static INT InitHexahedron (void)
{
  INT point[CORNERS_OF_BND_SEG],i,j;
  DOUBLE radius,MidPoint[3], alpha[DIM_OF_BND], beta[DIM_OF_BND];

  /* allocate new domain structure */
  MidPoint[0] = 0.0;
  MidPoint[1] = 0.0;
  MidPoint[2] = 0.0;
  for (i=0; i<8; i++) {
    MidPoint[0] += x_hex[i][0];
    MidPoint[1] += x_hex[i][1];
    MidPoint[2] += x_hex[i][2];
  }
  MidPoint[0] /= 8;
  MidPoint[1] /= 8;
  MidPoint[2] /= 8;
  radius = 0.0;
  for (i=0; i<8; i++)
    for (j=0; j<3; j++)
      radius = MAX(radius,ABS(MidPoint[j]-x_hex[i][j]));

  if (CreateDomain("Hexahedron",MidPoint,radius,6,8,YES)==NULL) return(1);

  /* allocate the boundary segments */
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  point[0]=0; point[1]=1; point[2]=2; point[3]=3;
  if (CreateBoundarySegment("south",1,0,0,NON_PERIODIC,1,point,alpha,beta,southBoundary,NULL)==NULL) return(1);
  point[0]=0; point[1]=3; point[2]=7; point[3]=4;
  if (CreateBoundarySegment("west", 1,0,1,NON_PERIODIC,1,point,alpha,beta,westBoundary, NULL)==NULL) return(1);
  point[0]=0; point[1]=1; point[2]=5; point[3]=4;
  if (CreateBoundarySegment("front",0,1,2,NON_PERIODIC,1,point,alpha,beta,frontBoundary,NULL)==NULL) return(1);
  point[0]=4; point[1]=5; point[2]=6; point[3]=7;
  if (CreateBoundarySegment("north",0,1,3,NON_PERIODIC,1,point,alpha,beta,northBoundary,NULL)==NULL) return(1);
  point[0]=1; point[1]=2; point[2]=6; point[3]=5;
  if (CreateBoundarySegment("east", 0,1,4,NON_PERIODIC,1,point,alpha,beta,eastBoundary, NULL)==NULL) return(1);
  point[0]=3; point[1]=2; point[2]=6; point[3]=7;
  if (CreateBoundarySegment("back", 1,0,5,NON_PERIODIC,1,point,alpha,beta,backBoundary, NULL)==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the unit ball                                                    */
/*                                                                          */
/****************************************************************************/

# define RADIUS 1.0

static void ProjectOnBall (DOUBLE x, DOUBLE y, DOUBLE z,
                           DOUBLE radius, DOUBLE *result)
{
  DOUBLE d;

  result[0] = x - 0.5;
  result[1] = y - 0.5;
  result[2] = z - 0.5;

  V3_EUKLIDNORM(result,d);
  d = radius / d;
  V3_SCALE(d,result);

  result[0] += 0.5;
  result[1] += 0.5;
  result[2] += 0.5;
}

static INT southboundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,lambda2,0.0,RADIUS,result);

  /* return ok */
  return(0);
}

static INT eastboundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (1.0,lambda1,lambda2,RADIUS,result);

  /* return ok */
  return(0);
}

static INT northboundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,lambda2,1.0,RADIUS,result);

  /* return ok */
  return(0);
}

static INT westboundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (0.0,lambda1,lambda2,RADIUS,result);

  /* return ok */
  return(0);
}

static INT frontboundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,0.0,lambda2,RADIUS,result);

  /* return ok */
  return(0);
}

static INT backboundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,1.0,lambda2,RADIUS,result);

  /* return ok */
  return(0);
}

static INT InitBall (void)
{
  INT point[CORNERS_OF_BND_SEG];
  DOUBLE radius,MidPoint[3], alpha[DIM_OF_BND], beta[DIM_OF_BND];

  /* allocate new domain structure */
  MidPoint[0] = 0.5;
  MidPoint[1] = 0.5;
  MidPoint[2] = 0.5;
  radius = 1.0;
  if (CreateDomain("Ball",MidPoint,radius,6,8,YES)==NULL) return(1);

  /* allocate the boundary segments */
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  point[0]=0; point[1]=1; point[2]=2; point[3]=3;
  if (CreateBoundarySegment("south",1,0,0,NON_PERIODIC,1,point,alpha,beta,southboundary,NULL)==NULL) return(1);
  point[0]=0; point[1]=3; point[2]=7; point[3]=4;
  if (CreateBoundarySegment("west", 1,0,1,NON_PERIODIC,1,point,alpha,beta,westboundary, NULL)==NULL) return(1);
  point[0]=0; point[1]=1; point[2]=5; point[3]=4;
  if (CreateBoundarySegment("front",0,1,2,NON_PERIODIC,1,point,alpha,beta,frontboundary,NULL)==NULL) return(1);
  point[0]=4; point[1]=5; point[2]=6; point[3]=7;
  if (CreateBoundarySegment("north",0,1,3,NON_PERIODIC,1,point,alpha,beta,northboundary,NULL)==NULL) return(1);
  point[0]=1; point[1]=2; point[2]=6; point[3]=5;
  if (CreateBoundarySegment("east", 0,1,4,NON_PERIODIC,1,point,alpha,beta,eastboundary, NULL)==NULL) return(1);
  point[0]=3; point[1]=2; point[2]=6; point[3]=7;
  if (CreateBoundarySegment("back", 1,0,5,NON_PERIODIC,1,point,alpha,beta,backboundary, NULL)==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the torus                                                        */
/*                                                                          */
/****************************************************************************/

#define R0              10.0
#define R1               4.0

#define PI_2    (0.5*PI)
#define PI_23   (1.5*PI)

/**************************** first quarter *********************************/

static INT Bnd_11 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(0+lambda));      s1 = sin(PI_2*(0+lambda));
  c2 = cos(PI_2*mu);                      s2 = sin(PI_2*mu);

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_12 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(0+lambda));      s1 = sin(PI_2*(0+lambda));
  c2 = cos(PI_2*(1+mu));      s2 = sin(PI_2*(1+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_13 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(0+lambda));      s1 = sin(PI_2*(0+lambda));
  c2 = cos(PI_2*(2+mu));      s2 = sin(PI_2*(2+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_14 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(0+lambda));      s1 = sin(PI_2*(0+lambda));
  c2 = cos(PI_2*(3+mu));      s2 = sin(PI_2*(3+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

/**************************** second quarter *********************************/

static INT Bnd_21 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(1+lambda));      s1 = sin(PI_2*(1+lambda));
  c2 = cos(PI_2*mu);                      s2 = sin(PI_2*mu);

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_22 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(1+lambda));      s1 = sin(PI_2*(1+lambda));
  c2 = cos(PI_2*(1+mu));      s2 = sin(PI_2*(1+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_23 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(1+lambda));      s1 = sin(PI_2*(1+lambda));
  c2 = cos(PI_2*(2+mu));      s2 = sin(PI_2*(2+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_24 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(1+lambda));      s1 = sin(PI_2*(1+lambda));
  c2 = cos(PI_2*(3+mu));      s2 = sin(PI_2*(3+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

/**************************** third quarter *********************************/

static INT Bnd_31 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(2+lambda));      s1 = sin(PI_2*(2+lambda));
  c2 = cos(PI_2*mu);                      s2 = sin(PI_2*mu);

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_32 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(2+lambda));      s1 = sin(PI_2*(2+lambda));
  c2 = cos(PI_2*(1+mu));      s2 = sin(PI_2*(1+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_33 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(2+lambda));      s1 = sin(PI_2*(2+lambda));
  c2 = cos(PI_2*(2+mu));      s2 = sin(PI_2*(2+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_34 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(2+lambda));      s1 = sin(PI_2*(2+lambda));
  c2 = cos(PI_2*(3+mu));      s2 = sin(PI_2*(3+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

/**************************** fourth quarter *********************************/

static INT Bnd_41 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(3+lambda));      s1 = sin(PI_2*(3+lambda));
  c2 = cos(PI_2*mu);                      s2 = sin(PI_2*mu);

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_42 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(3+lambda));      s1 = sin(PI_2*(3+lambda));
  c2 = cos(PI_2*(1+mu));      s2 = sin(PI_2*(1+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_43 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(3+lambda));      s1 = sin(PI_2*(3+lambda));
  c2 = cos(PI_2*(2+mu));      s2 = sin(PI_2*(2+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}

static INT Bnd_44 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,mu,c1,s1,c2,s2;

  lambda = param[0];
  mu     = param[1];

  /* check range */
  if ( lambda<0.0 || lambda>1.0 || mu<0.0 || mu>1.0 ) return(1);

  c1 = cos(PI_2*(3+lambda));      s1 = sin(PI_2*(3+lambda));
  c2 = cos(PI_2*(3+mu));      s2 = sin(PI_2*(3+mu));

  /* fill result */
  result[0] = (R0+c2*R1)*c1;
  result[1] = (R0+c2*R1)*s1;
  result[2] = R1*s2;

  /* return ok */
  return(0);
}


static INT InitTorus (void)
{
  INT point[CORNERS_OF_BND_SEG];
  DOUBLE radius,MidPoint[3], alpha[DIM_OF_BND], beta[DIM_OF_BND];

  /* allocate new domain structure */
  MidPoint[0] = 0.0;
  MidPoint[1] = 0.0;
  MidPoint[2] = 0.0;
  radius = (R0+R1) * 1.1;
  if (CreateDomain("Torus",MidPoint,radius,16,16,NO)==NULL) return(1);

  /* allocate the boundary segments */
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;

  /* first quarter */
  point[0]=0; point[1]=4; point[2]=5; point[3]=1; if (CreateBoundarySegment("Bnd_11",0,1,0 ,NON_PERIODIC,10,point,alpha,beta,Bnd_11,NULL)==NULL) return(1);
  point[0]=1; point[1]=5; point[2]=6; point[3]=2; if (CreateBoundarySegment("Bnd_12",0,1,1 ,NON_PERIODIC,10,point,alpha,beta,Bnd_12,NULL)==NULL) return(1);
  point[0]=2; point[1]=6; point[2]=7; point[3]=3; if (CreateBoundarySegment("Bnd_13",0,1,2 ,NON_PERIODIC,10,point,alpha,beta,Bnd_13,NULL)==NULL) return(1);
  point[0]=3; point[1]=7; point[2]=4; point[3]=0; if (CreateBoundarySegment("Bnd_14",0,1,3 ,NON_PERIODIC,10,point,alpha,beta,Bnd_14,NULL)==NULL) return(1);
  /* second quarter */
  point[0]=4; point[1]=8; point[2]=9; point[3]=5; if (CreateBoundarySegment("Bnd_21",0,1,4 ,NON_PERIODIC,10,point,alpha,beta,Bnd_21,NULL)==NULL) return(1);
  point[0]=5; point[1]=9; point[2]=10;point[3]=6; if (CreateBoundarySegment("Bnd_22",0,1,5 ,NON_PERIODIC,10,point,alpha,beta,Bnd_22,NULL)==NULL) return(1);
  point[0]=6; point[1]=10;point[2]=11;point[3]=7; if (CreateBoundarySegment("Bnd_23",0,1,6 ,NON_PERIODIC,10,point,alpha,beta,Bnd_23,NULL)==NULL) return(1);
  point[0]=7; point[1]=11;point[2]=8; point[3]=4; if (CreateBoundarySegment("Bnd_24",0,1,7 ,NON_PERIODIC,10,point,alpha,beta,Bnd_24,NULL)==NULL) return(1);
  /* third quarter */
  point[0]=8; point[1]=12;point[2]=13;point[3]=9; if (CreateBoundarySegment("Bnd_31",0,1,8 ,NON_PERIODIC,10,point,alpha,beta,Bnd_31,NULL)==NULL) return(1);
  point[0]=9; point[1]=13;point[2]=14;point[3]=10;if (CreateBoundarySegment("Bnd_32",0,1,9 ,NON_PERIODIC,10,point,alpha,beta,Bnd_32,NULL)==NULL) return(1);
  point[0]=10;point[1]=14;point[2]=15;point[3]=11;if (CreateBoundarySegment("Bnd_33",0,1,10,NON_PERIODIC,10,point,alpha,beta,Bnd_33,NULL)==NULL) return(1);
  point[0]=11;point[1]=15;point[2]=12;point[3]=8; if (CreateBoundarySegment("Bnd_34",0,1,11,NON_PERIODIC,10,point,alpha,beta,Bnd_34,NULL)==NULL) return(1);
  /* fourth quarter */
  point[0]=12;point[1]=0; point[2]=1;point[3]=13; if (CreateBoundarySegment("Bnd_41",0,1,12,NON_PERIODIC,10,point,alpha,beta,Bnd_41,NULL)==NULL) return(1);
  point[0]=13;point[1]=1; point[2]=2;point[3]=14; if (CreateBoundarySegment("Bnd_42",0,1,13,NON_PERIODIC,10,point,alpha,beta,Bnd_42,NULL)==NULL) return(1);
  point[0]=14;point[1]=2; point[2]=3;point[3]=15; if (CreateBoundarySegment("Bnd_43",0,1,14,NON_PERIODIC,10,point,alpha,beta,Bnd_43,NULL)==NULL) return(1);
  point[0]=15;point[1]=3; point[2]=0;point[3]=12; if (CreateBoundarySegment("Bnd_44",0,1,15,NON_PERIODIC,10,point,alpha,beta,Bnd_44,NULL)==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*
   InitFEMDomains - Calls all inits of format definitions

   SYNOPSIS:
   INT InitFEMDomains (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initialize the domains 'unit square',
   'reentrant corner' and 'cook'.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
 */
/****************************************************************************/

INT STD_BVP_Configure (INT argc, char **argv)
{
  STD_BVP *theBVP;
  DOMAIN *theDomain;
  char BVPName[NAMESIZE];
  char DomainName[NAMESIZE];
  INT i;

  /* get BVP name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" configure %",NAMELENSTR,"[ -~]")),
              BVPName)!=1) || (strlen(BVPName)==0))
    return(1);

  theBVP = (STD_BVP *) BVP_GetByName(BVPName);
  if (theBVP == NULL)
    return(1);

  for (i=0; i<argc; i++)
    if (argv[i][0] == 'd')
      if ((sscanf(argv[i],expandfmt(CONCAT3("d %",NAMELENSTR,"[ -~]")),
                  DomainName)!=1) || (strlen(DomainName)==0))
        continue;

  theDomain = GetDomain(DomainName);

  if (strcmp(DomainName,"Hexahedron") == 0) {
    if (ReadArgvPosition("x0",argc,argv,x_hex[0])) {
      x_hex[0][0] = 0.0;
      x_hex[0][1] = 0.0;
      x_hex[0][2] = 0.0;
    }
    if (ReadArgvPosition("x1",argc,argv,x_hex[1])) {
      x_hex[1][0] = 1.0;
      x_hex[1][1] = 0.0;
      x_hex[1][2] = 0.0;
    }
    if (ReadArgvPosition("x2",argc,argv,x_hex[2])) {
      x_hex[2][0] = 1.0;
      x_hex[2][1] = 1.0;
      x_hex[2][2] = 0.0;
    }
    if (ReadArgvPosition("x3",argc,argv,x_hex[3])) {
      x_hex[3][0] = 0.0;
      x_hex[3][1] = 1.0;
      x_hex[3][2] = 0.0;
    }
    if (ReadArgvPosition("x4",argc,argv,x_hex[4])) {
      x_hex[4][0] = 0.0;
      x_hex[4][1] = 0.0;
      x_hex[4][2] = 1.0;
    }
    if (ReadArgvPosition("x5",argc,argv,x_hex[5])) {
      x_hex[5][0] = 1.0;
      x_hex[5][1] = 0.0;
      x_hex[5][2] = 1.0;
    }
    if (ReadArgvPosition("x6",argc,argv,x_hex[6])) {
      x_hex[6][0] = 1.0;
      x_hex[6][1] = 1.0;
      x_hex[6][2] = 1.0;
    }
    if (ReadArgvPosition("x7",argc,argv,x_hex[7])) {
      x_hex[7][0] = 0.0;
      x_hex[7][1] = 1.0;
      x_hex[7][2] = 1.0;
    }
  }

  if (theDomain == NULL)
  {
    if (strcmp(DomainName,"Ball") == 0)
    {
      if (InitBall())
        return(1);
    }
    else if (strcmp(DomainName,"Hexahedron") == 0)
    {
      if (InitHexahedron())
        return(1);
    }
    else if (strcmp(DomainName,"Torus") == 0)
    {
      if (InitTorus())
        return(1);
    }
    else
      return(1);

    theDomain = GetDomain(DomainName);

    if (theDomain == NULL)
      return(1);
  }

  theBVP->Domain = theDomain;

  return(0);
}
