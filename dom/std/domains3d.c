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

#define RADIUS 1.0

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
/*                                                                          */
/*  define the cylinder                                                     */
/*                                                                          */
/****************************************************************************/



static INT ZylFront1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.0 || lambda2>0.075 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.075 || lambda2>0.2 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.2 || lambda2>0.33 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.33 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront5 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.4 || lambda1>0.5 || lambda2<0.0 || lambda2>0.075 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront6 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4 + 0.1*lambda2;
  pz = 0.075;
  qx = 0.5 + 0.05*cos(5.0*PI*0.25 + 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(5.0*PI*0.25 + 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront7 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4;
  pz = 0.075 + 0.125*lambda2;
  qx = 0.5 + 0.05*cos(5.0*PI*0.25 - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(5.0*PI*0.25 - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront8 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4;
  pz = 0.2 + 0.13*lambda2;
  qx = 0.5 + 0.05*cos(PI - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(PI - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront9 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4 + 0.1*lambda2;
  pz = 0.33;
  qx = 0.5 + 0.05*cos(3.0*PI*0.25 - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(3.0*PI*0.25 - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront10 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.4 || lambda1>0.5 || lambda2<0.33 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront11 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.5 || lambda1>0.6 || lambda2<0.0 || lambda2>0.075 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront12 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.5 + 0.1*lambda2;
  pz = 0.075;
  qx = 0.5 + 0.05*cos(3.0*PI*0.5 + 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(3.0*PI*0.5 + 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront13 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.6;
  pz = 0.075 + 0.125*lambda2;
  qx = 0.5 + 0.05*cos(7.0*PI*0.25 + 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(7.0*PI*0.25 + 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront14 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.6;
  pz = 0.2 + 0.13*lambda2;
  qx = 0.5 + 0.05*cos(0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront15 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.5 + 0.1*lambda2;
  pz = 0.33;
  qx = 0.5 + 0.05*cos(PI*0.5 - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(PI*0.5 - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.0;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylFront16 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.5 || lambda1>0.6 || lambda2<0.33 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront17 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.0 || lambda2>0.075 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront18 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.075 || lambda2>0.2 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront19 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.2 || lambda2>0.33 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylFront20 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.33 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.0 || lambda2>0.075 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.075 || lambda2>0.2 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.2 || lambda2>0.33 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.33 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack5 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.4 || lambda1>0.5 || lambda2<0.0 || lambda2>0.075 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack6 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4 + 0.1*lambda2;
  pz = 0.075;
  qx = 0.5 + 0.05*cos(5.0*PI*0.25 + 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(5.0*PI*0.25 + 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack7 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4;
  pz = 0.075 + 0.125*lambda2;
  qx = 0.5 + 0.05*cos(5.0*PI*0.25 - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(5.0*PI*0.25 - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack8 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4;
  pz = 0.2 + 0.13*lambda2;
  qx = 0.5 + 0.05*cos(PI - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(PI - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack9 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.4 + 0.1*lambda2;
  pz = 0.33;
  qx = 0.5 + 0.05*cos(3.0*PI*0.25 - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(3.0*PI*0.25 - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack10 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.4 || lambda1>0.5 || lambda2<0.33 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack11 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.5 || lambda1>0.6 || lambda2<0.0 || lambda2>0.075 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack12 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.5 + 0.1*lambda2;
  pz = 0.075;
  qx = 0.5 + 0.05*cos(3.0*PI*0.5 + 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(3.0*PI*0.5 + 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack13 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.6;
  pz = 0.075 + 0.125*lambda2;
  qx = 0.5 + 0.05*cos(7.0*PI*0.25 + 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(7.0*PI*0.25 + 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack14 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.6;
  pz = 0.2 + 0.13*lambda2;
  qx = 0.5 + 0.05*cos(0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack15 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 0.5 + 0.1*lambda2;
  pz = 0.33;
  qx = 0.5 + 0.05*cos(PI*0.5 - 0.25*PI*lambda2);
  qz = 0.2 + 0.05*sin(PI*0.5 - 0.25*PI*lambda2);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = 0.41;
  result[2] = (1.0 - lambda1)*pz + lambda1*qz;


  /* return ok */
  return(0);
}

static INT ZylBack16 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.5 || lambda1>0.6 || lambda2<0.33 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack17 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.0 || lambda2>0.075 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack18 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.075 || lambda2>0.2 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack19 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.2 || lambda2>0.33 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylBack20 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.33 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.41;
  result[2] = lambda2;


  /* return ok */
  return(0);
}

static INT ZylSouth1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;


  /* return ok */
  return(0);
}

static INT ZylSouth2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.4 || lambda1>0.5 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;


  /* return ok */
  return(0);
}

static INT ZylSouth3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.5 || lambda1>0.6 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;


  /* return ok */
  return(0);
}

static INT ZylSouth4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;


  /* return ok */
  return(0);
}

static INT ZylNorth1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.4 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.41;


  /* return ok */
  return(0);
}

static INT ZylNorth2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.4 || lambda1>0.5 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.41;


  /* return ok */
  return(0);
}

static INT ZylNorth3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.5 || lambda1>0.6 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.41;


  /* return ok */
  return(0);
}

static INT ZylNorth4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.6 || lambda1>2.5 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.41;


  /* return ok */
  return(0);
}

static INT ZylWest1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.075 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 0.0;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT ZylWest2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.075 || lambda1>0.2 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 0.0;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT ZylWest3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.2 || lambda1>0.33 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 0.0;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT ZylWest4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.33 || lambda1>0.41 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 0.0;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT ZylEast1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>0.075 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 2.5;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT ZylEast2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.075 || lambda1>0.2 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 2.5;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT ZylEast3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.2 || lambda1>0.33 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 2.5;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT ZylEast4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2, px, pz, qx, qz;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.33 || lambda1>0.41 || lambda2<0.0 || lambda2>0.41 ) return(1);
  /* fill result */
  result[0] = 2.5;
  result[1] = lambda2;
  result[2] = lambda1;


  /* return ok */
  return(0);
}

static INT Zyl1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(PI*0.25*lambda1);

  /* return ok */
  return(0);
}

static INT Zyl2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(PI*0.25 + PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(PI*0.25 + PI*0.25*lambda1);

  /* return ok */
  return(0);
}

static INT Zyl3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(PI*0.5 + PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(PI*0.5 + PI*0.25*lambda1);

  /* return ok */
  return(0);
}

static INT Zyl4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(0.75*PI + PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(0.75*PI + PI*0.25*lambda1);

  /* return ok */
  return(0);
}

static INT Zyl5 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(PI + PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(PI + PI*0.25*lambda1);

  /* return ok */
  return(0);
}

static INT Zyl6 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(1.25*PI + PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(1.25*PI + PI*0.25*lambda1);

  /* return ok */
  return(0);
}

static INT Zyl7 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(1.5*PI + PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(1.5*PI + PI*0.25*lambda1);

  /* return ok */
  return(0);
}

static INT Zyl8 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>0.41 ) return(1);

  /* fill result */
  result[0] = 0.5 + 0.05*cos(1.75*PI + PI*0.25*lambda1);
  result[1] = lambda2;
  result[2] = 0.2 + 0.05*sin(1.75*PI + PI*0.25*lambda1);

  /* return ok */
  return(0);
}


static INT InitCylinder (void)
{
  INT point[CORNERS_OF_BND_SEG];
  DOUBLE radius,MidPoint[3], alpha[DIM_OF_BND], beta[DIM_OF_BND];

  /* allocate new domain structure */
  MidPoint[0] = 1.25;
  MidPoint[1] = 0.205;
  MidPoint[2] = 0.205;
  radius = 1.6;
  if (CreateDomain("Cylinder",MidPoint,radius,64,64,NO)==NULL) return(1);

  /* allocate the boundary segments */
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =0.4; beta[1] =0.075;
  point[0]=0; point[1]=5; point[2]=6; point[3]=1;
  if (CreateBoundarySegment("zylfront1",0,1,0,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront1,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.075;
  beta[0] =0.4; beta[1] =0.2;
  point[0]=1; point[1]=6; point[2]=7; point[3]=2;
  if (CreateBoundarySegment("zylfront2", 0,1,1,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront2, NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.2;
  beta[0] =0.4; beta[1] =0.33;
  point[0]=2; point[1]=7; point[2]=8; point[3]=3;
  if (CreateBoundarySegment("zylfront3",0,1,2,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront3,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.33;
  beta[0] =0.4; beta[1] =0.41;
  point[0]=3; point[1]=8; point[2]=9; point[3]=4;
  if (CreateBoundarySegment("zylfront4",0,1,3,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront4,NULL)==NULL) return(1);

  alpha[0]=0.4; alpha[1]=0.0;
  beta[0] =0.5; beta[1] =0.075;
  point[0]=5; point[1]=13; point[2]=14; point[3]=6;
  if (CreateBoundarySegment("zylfront5",0,1,4,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront5,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  point[0]=6; point[1]=11; point[2]=15; point[3]=14;
  if (CreateBoundarySegment("zylfront6",0,1,5,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront6,NULL)==NULL) return(1);

  point[0]=6; point[1]=11; point[2]=10; point[3]=7;
  if (CreateBoundarySegment("zylfront7",0,1,6,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront7,NULL)==NULL) return(1);

  point[0]=7; point[1]=10; point[2]=12; point[3]=8;
  if (CreateBoundarySegment("zylfront8",0,1,7,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront8,NULL)==NULL) return(1);

  point[0]=8; point[1]=12; point[2]=16; point[3]=17;
  if (CreateBoundarySegment("zylfront9",0,1,8,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront9,NULL)==NULL) return(1);

  alpha[0]=0.4; alpha[1]=0.33;
  beta[0] =0.5; beta[1] =0.41;
  point[0]=8; point[1]=17; point[2]=18; point[3]=9;
  if (CreateBoundarySegment("zylfront10",0,1,9,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront10,NULL)==NULL) return(1);

  alpha[0]=0.5; alpha[1]=0.0;
  beta[0] =0.6; beta[1] =0.075;
  point[0]=13; point[1]=22; point[2]=23; point[3]=14;
  if (CreateBoundarySegment("zylfront11",0,1,10,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront11,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  point[0]=14; point[1]=15; point[2]=19; point[3]=23;
  if (CreateBoundarySegment("zylfront12",0,1,11,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront12,NULL)==NULL) return(1);

  point[0]=23; point[1]=19; point[2]=21; point[3]=24;
  if (CreateBoundarySegment("zylfront13",0,1,12,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront13,NULL)==NULL) return(1);

  point[0]=24; point[1]=21; point[2]=20; point[3]=25;
  if (CreateBoundarySegment("zylfront14",0,1,13,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront14,NULL)==NULL) return(1);

  point[0]=17; point[1]=16; point[2]=20; point[3]=25;
  if (CreateBoundarySegment("zylfront15",0,1,14,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront15,NULL)==NULL) return(1);

  alpha[0]=0.5; alpha[1]=0.33;
  beta[0] =0.6; beta[1] =0.41;
  point[0]=17; point[1]=25; point[2]=26; point[3]=18;
  if (CreateBoundarySegment("zylfront16",0,1,15,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront16,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =0.075;
  point[0]=22; point[1]=27; point[2]=28; point[3]=23;
  if (CreateBoundarySegment("zylfront17",0,1,16,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront17,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.075;
  beta[0] =2.5; beta[1] =0.2;
  point[0]=23; point[1]=28; point[2]=29; point[3]=24;
  if (CreateBoundarySegment("zylfront18",0,1,17,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront18,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.2;
  beta[0] =2.5; beta[1] =0.33;
  point[0]=24; point[1]=29; point[2]=30; point[3]=25;
  if (CreateBoundarySegment("zylfront19",0,1,18,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront19,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.33;
  beta[0] =2.5; beta[1] =0.41;
  point[0]=25; point[1]=30; point[2]=31; point[3]=26;
  if (CreateBoundarySegment("zylfront20",0,1,19,NON_PERIODIC,1,point,
                            alpha,beta,ZylFront20,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =0.4; beta[1] =0.075;
  point[0]=32; point[1]=37; point[2]=38; point[3]=33;
  if (CreateBoundarySegment("zylback1",1,0,20,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack1,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.075;
  beta[0] =0.4; beta[1] =0.2;
  point[0]=33; point[1]=38; point[2]=39; point[3]=34;
  if (CreateBoundarySegment("zylback2", 1,0,21,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack2, NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.2;
  beta[0] =0.4; beta[1] =0.33;
  point[0]=34; point[1]=39; point[2]=40; point[3]=35;
  if (CreateBoundarySegment("zylback3",1,0,22,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack3,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.33;
  beta[0] =0.4; beta[1] =0.41;
  point[0]=35; point[1]=40; point[2]=41; point[3]=36;
  if (CreateBoundarySegment("zylback4",1,0,23,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack4,NULL)==NULL) return(1);

  alpha[0]=0.4; alpha[1]=0.0;
  beta[0] =0.5; beta[1] =0.075;
  point[0]=37; point[1]=45; point[2]=46; point[3]=38;
  if (CreateBoundarySegment("zylback5",1,0,24,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack5,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  point[0]=38; point[1]=43; point[2]=47; point[3]=46;
  if (CreateBoundarySegment("zylback6",1,0,25,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack6,NULL)==NULL) return(1);

  point[0]=38; point[1]=43; point[2]=42; point[3]=39;
  if (CreateBoundarySegment("zylback7",1,0,26,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack7,NULL)==NULL) return(1);

  point[0]=39; point[1]=42; point[2]=44; point[3]=40;
  if (CreateBoundarySegment("zylback8",1,0,27,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack8,NULL)==NULL) return(1);

  point[0]=40; point[1]=44; point[2]=48; point[3]=49;
  if (CreateBoundarySegment("zylback9",1,0,28,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack9,NULL)==NULL) return(1);

  alpha[0]=0.4; alpha[1]=0.33;
  beta[0] =0.5; beta[1] =0.41;
  point[0]=40; point[1]=49; point[2]=50; point[3]=41;
  if (CreateBoundarySegment("zylback10",1,0,29,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack10,NULL)==NULL) return(1);

  alpha[0]=0.5; alpha[1]=0.0;
  beta[0] =0.6; beta[1] =0.075;
  point[0]=45; point[1]=54; point[2]=55; point[3]=46;
  if (CreateBoundarySegment("zylback11",1,0,30,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack11,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  point[0]=46; point[1]=47; point[2]=51; point[3]=55;
  if (CreateBoundarySegment("zylback12",1,0,31,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack12,NULL)==NULL) return(1);

  point[0]=55; point[1]=51; point[2]=53; point[3]=56;
  if (CreateBoundarySegment("zylback13",1,0,32,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack13,NULL)==NULL) return(1);

  point[0]=56; point[1]=53; point[2]=52; point[3]=57;
  if (CreateBoundarySegment("zylback14",1,0,33,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack14,NULL)==NULL) return(1);

  point[0]=49; point[1]=48; point[2]=52; point[3]=57;
  if (CreateBoundarySegment("zylback15",1,0,34,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack15,NULL)==NULL) return(1);

  alpha[0]=0.5; alpha[1]=0.33;
  beta[0] =0.6; beta[1] =0.41;
  point[0]=49; point[1]=57; point[2]=58; point[3]=50;
  if (CreateBoundarySegment("zylback16",1,0,35,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack16,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =0.075;
  point[0]=54; point[1]=59; point[2]=60; point[3]=55;
  if (CreateBoundarySegment("zylback17",1,0,36,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack17,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.075;
  beta[0] =2.5; beta[1] =0.2;
  point[0]=55; point[1]=60; point[2]=61; point[3]=56;
  if (CreateBoundarySegment("zylback18",1,0,37,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack18,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.2;
  beta[0] =2.5; beta[1] =0.33;
  point[0]=56; point[1]=61; point[2]=62; point[3]=57;
  if (CreateBoundarySegment("zylback19",1,0,38,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack19,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.33;
  beta[0] =2.5; beta[1] =0.41;
  point[0]=57; point[1]=62; point[2]=63; point[3]=58;
  if (CreateBoundarySegment("zylback20",1,0,39,NON_PERIODIC,1,point,
                            alpha,beta,ZylBack20,NULL)==NULL) return(1);



  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =0.4; beta[1] =0.41;
  point[0]=0; point[1]=5; point[2]=37; point[3]=32;
  if (CreateBoundarySegment("zylsouth1",1,0,40,NON_PERIODIC,1,point,
                            alpha,beta,ZylSouth1,NULL)==NULL) return(1);

  alpha[0]=0.4; alpha[1]=0.0;
  beta[0] =0.5; beta[1] =0.41;
  point[0]=5; point[1]=13; point[2]=45; point[3]=37;
  if (CreateBoundarySegment("zylsouth2",1,0,41,NON_PERIODIC,1,point,
                            alpha,beta,ZylSouth2,NULL)==NULL) return(1);

  alpha[0]=0.5; alpha[1]=0.0;
  beta[0] =0.6; beta[1] =0.41;
  point[0]=13; point[1]=22; point[2]=54; point[3]=45;
  if (CreateBoundarySegment("zylsouth3",1,0,42,NON_PERIODIC,1,point,
                            alpha,beta,ZylSouth3,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =0.41;
  point[0]=22; point[1]=27; point[2]=59; point[3]=54;
  if (CreateBoundarySegment("zylsouth4",1,0,43,NON_PERIODIC,1,point,
                            alpha,beta,ZylSouth4,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =0.4; beta[1] =0.41;
  point[0]=4; point[1]=9; point[2]=41; point[3]=36;
  if (CreateBoundarySegment("zylnorth1",0,1,44,NON_PERIODIC,1,point,
                            alpha,beta,ZylNorth1,NULL)==NULL) return(1);

  alpha[0]=0.4; alpha[1]=0.0;
  beta[0] =0.5; beta[1] =0.41;
  point[0]=9; point[1]=18; point[2]=50; point[3]=41;
  if (CreateBoundarySegment("zylnorth2",0,1,45,NON_PERIODIC,1,point,
                            alpha,beta,ZylNorth2,NULL)==NULL) return(1);

  alpha[0]=0.5; alpha[1]=0.0;
  beta[0] =0.6; beta[1] =0.41;
  point[0]=18; point[1]=26; point[2]=58; point[3]=50;
  if (CreateBoundarySegment("zylnorth3",0,1,46,NON_PERIODIC,1,point,
                            alpha,beta,ZylNorth3,NULL)==NULL) return(1);

  alpha[0]=0.6; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =0.41;
  point[0]=26; point[1]=31; point[2]=63; point[3]=58;
  if (CreateBoundarySegment("zylnorth4",0,1,47,NON_PERIODIC,1,point,
                            alpha,beta,ZylNorth4,NULL)==NULL) return(1);


  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =0.075; beta[1] =0.41;
  point[0]=0; point[1]=1; point[2]=33; point[3]=32;
  if (CreateBoundarySegment("zylwest1",1,0,48,NON_PERIODIC,1,point,
                            alpha,beta,ZylWest1,NULL)==NULL) return(1);

  alpha[0]=0.075; alpha[1]=0.0;
  beta[0] =0.2; beta[1] =0.41;
  point[0]=1; point[1]=2; point[2]=34; point[3]=33;
  if (CreateBoundarySegment("zylwest2",1,0,49,NON_PERIODIC,1,point,
                            alpha,beta,ZylWest2,NULL)==NULL) return(1);

  alpha[0]=0.2; alpha[1]=0.0;
  beta[0] =0.33; beta[1] =0.41;
  point[0]=2; point[1]=3; point[2]=35; point[3]=34;
  if (CreateBoundarySegment("zylwest3",1,0,50,NON_PERIODIC,1,point,
                            alpha,beta,ZylWest3,NULL)==NULL) return(1);

  alpha[0]=0.33; alpha[1]=0.0;
  beta[0] =0.41; beta[1] =0.41;
  point[0]=3; point[1]=4; point[2]=36; point[3]=35;
  if (CreateBoundarySegment("zylwest4",1,0,51,NON_PERIODIC,1,point,
                            alpha,beta,ZylWest4,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =0.075; beta[1] =0.41;
  point[0]=27; point[1]=28; point[2]=60; point[3]=59;
  if (CreateBoundarySegment("zyleast1",0,1,52,NON_PERIODIC,1,point,
                            alpha,beta,ZylEast1,NULL)==NULL) return(1);

  alpha[0]=0.075; alpha[1]=0.0;
  beta[0] =0.2; beta[1] =0.41;
  point[0]=28; point[1]=29; point[2]=61; point[3]=60;
  if (CreateBoundarySegment("zyleast2",0,1,53,NON_PERIODIC,1,point,
                            alpha,beta,ZylEast2,NULL)==NULL) return(1);

  alpha[0]=0.2; alpha[1]=0.0;
  beta[0] =0.33; beta[1] =0.41;
  point[0]=29; point[1]=30; point[2]=62; point[3]=61;
  if (CreateBoundarySegment("zyleast3",0,1,54,NON_PERIODIC,1,point,
                            alpha,beta,ZylEast3,NULL)==NULL) return(1);

  alpha[0]=0.33; alpha[1]=0.0;
  beta[0] =0.41; beta[1] =0.41;
  point[0]=30; point[1]=31; point[2]=63; point[3]=62;
  if (CreateBoundarySegment("zyleast4",0,1,55,NON_PERIODIC,1,point,
                            alpha,beta,ZylEast4,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =0.41;
  point[0]=21; point[1]=20; point[2]=52; point[3]=53;
  if (CreateBoundarySegment("z1",0,1,56,NON_PERIODIC,20,point,
                            alpha,beta,Zyl1,NULL)==NULL) return(1);

  point[0]=20; point[1]=16; point[2]=48; point[3]=52;
  if (CreateBoundarySegment("z2",0,1,57,NON_PERIODIC,20,point,
                            alpha,beta,Zyl2,NULL)==NULL) return(1);

  point[0]=16; point[1]=12; point[2]=44; point[3]=48;
  if (CreateBoundarySegment("z3",0,1,58,NON_PERIODIC,20,point,
                            alpha,beta,Zyl3,NULL)==NULL) return(1);

  point[0]=12; point[1]=10; point[2]=42; point[3]=44;
  if (CreateBoundarySegment("z4",0,1,59,NON_PERIODIC,20,point,
                            alpha,beta,Zyl4,NULL)==NULL) return(1);

  point[0]=10; point[1]=11; point[2]=43; point[3]=42;
  if (CreateBoundarySegment("z5",0,1,60,NON_PERIODIC,20,point,
                            alpha,beta,Zyl5,NULL)==NULL) return(1);

  point[0]=11; point[1]=15; point[2]=47; point[3]=43;
  if (CreateBoundarySegment("z6",0,1,61,NON_PERIODIC,20,point,
                            alpha,beta,Zyl6,NULL)==NULL) return(1);

  point[0]=15; point[1]=19; point[2]=51; point[3]=47;
  if (CreateBoundarySegment("z7",0,1,62,NON_PERIODIC,20,point,
                            alpha,beta,Zyl7,NULL)==NULL) return(1);

  point[0]=19; point[1]=21; point[2]=53; point[3]=51;
  if (CreateBoundarySegment("z8",0,1,63,NON_PERIODIC,20,point,
                            alpha,beta,Zyl8,NULL)==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the benchmark domain                                             */
/*                                                                          */
/****************************************************************************/

#define BENCH_H 1.0


static INT bsouth1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<0.0 || lambda2>2.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth5Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<0.0 || lambda2>2.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth6Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth7Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth8Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth9Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<0.0 || lambda2>2.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth10Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth11Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth12Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth13Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 7.5;
  py = 2.5*lambda2;
  qx = 10.0 + cos(PI - lambda2*PI*0.25);
  qy = sin(PI - lambda2*PI*0.25);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = (1.0 - lambda1)*py + lambda1*qy;
  result[2] = 0.0;



  /* return ok */
  return(0);
}

static INT bsouth14Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 7.5 + 2.5*lambda2;
  py = 2.5;
  qx = 10.0 + cos(0.75*PI - lambda2*PI*0.25);
  qy = sin(0.75*PI - lambda2*PI*0.25);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = (1.0 - lambda1)*py + lambda1*qy;
  result[2] = 0.0;



  /* return ok */
  return(0);
}

static INT bsouth15Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth16Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT bsouth17Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}


static INT bnorth1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<0.0 || lambda2>2.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth5Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<0.0 || lambda2>2.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth6Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth7Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth8Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth9Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<0.0 || lambda2>2.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth10Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth11Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth12Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth13Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 7.5;
  py = 2.5*lambda2;
  qx = 10.0 + cos(PI - lambda2*PI*0.25);
  qy = sin(PI - lambda2*PI*0.25);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = (1.0 - lambda1)*py + lambda1*qy;
  result[2] = BENCH_H;



  /* return ok */
  return(0);
}

static INT bnorth14Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  px = 7.5 + 2.5*lambda2;
  py = 2.5;
  qx = 10.0 + cos(0.75*PI - lambda2*PI*0.25);
  qy = sin(0.75*PI - lambda2*PI*0.25);
  result[0] = (1.0 - lambda1)*px + lambda1*qx;
  result[1] = (1.0 - lambda1)*py + lambda1*qy;
  result[2] = BENCH_H;



  /* return ok */
  return(0);
}

static INT bnorth15Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<2.5 || lambda2>5.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth16Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<5.0 || lambda2>7.5 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bnorth17Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2, px, py, qx, qy;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<7.5 || lambda2>10.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = BENCH_H;

  /* return ok */
  return(0);
}

static INT bfront1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bfront2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bfront3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bfront4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>9.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT beast1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;
  double d;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 10.0 + cos(PI - 0.25*lambda1*PI);
  result[1] = sin(PI - 0.25*lambda1*PI);

  result[2] = lambda2;
  /* return ok */
  return(0);
}

static INT beast2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;
  double d;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 10.0 + cos(0.75*PI - 0.25*lambda1*PI);
  result[1] = sin(0.75*PI - 0.25*lambda1*PI);

  result[2] = lambda2;
  /* return ok */
  return(0);
}

static INT beast3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<1.0 || lambda1>2.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 10.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT beast4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 10.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT beast5Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 10.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT beast6Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 10.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bback1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 10.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bback2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 10.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bback3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 10.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bback4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 10.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bwest1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>2.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 0.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bwest2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<2.5 || lambda1>5.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 0.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bwest3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<5.0 || lambda1>7.5 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 0.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT bwest4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1, lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<7.5 || lambda1>10.0 || lambda2<0.0 || lambda2>BENCH_H ) return(1);

  /* fill result */
  result[0] = 0.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT InitBenchmark (void)
{
  INT point[CORNERS_OF_BND_SEG];
  DOUBLE radius,MidPoint[3], alpha[DIM_OF_BND], beta[DIM_OF_BND];

  /* allocate new domain structure */
  MidPoint[0] = 5.0;
  MidPoint[1] = 5.0;
  MidPoint[2] = 2.5;
  radius = 7.0;
  if (CreateDomain("Benchmark",MidPoint,radius,52,54,NO)==NULL) return(1);

  /* allocate the boundary segments */
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =2.5;
  point[0]=0; point[1]=5; point[2]=6; point[3]=1;
  if (CreateBoundarySegment("bsouth1",1,0,0,NON_PERIODIC,1,point,
                            alpha,beta,bsouth1Boundary,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=2.5;
  beta[0] =2.5; beta[1] =5.0;
  point[0]=1; point[1]=6; point[2]=7; point[3]=2;
  if (CreateBoundarySegment("bsouth2", 1,0,1,NON_PERIODIC,1,point,
                            alpha,beta,bsouth2Boundary, NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=5.0;
  beta[0] =2.5; beta[1] =7.5;
  point[0]=2; point[1]=7; point[2]=8; point[3]=3;
  if (CreateBoundarySegment("bsouth3",1,0,2,NON_PERIODIC,1,point,
                            alpha,beta,bsouth3Boundary,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=7.5;
  beta[0] =2.5; beta[1] =10.0;
  point[0]=3; point[1]=8; point[2]=9; point[3]=4;
  if (CreateBoundarySegment("bsouth4",1,0,3,NON_PERIODIC,1,point,
                            alpha,beta,bsouth4Boundary,NULL)==NULL) return(1);

  alpha[0]=2.5; alpha[1]=0.0;
  beta[0] =5.0; beta[1] =2.5;
  point[0]=5; point[1]=10; point[2]=11; point[3]=6;
  if (CreateBoundarySegment("bsouth5",1,0,4,NON_PERIODIC,1,point,
                            alpha,beta,bsouth5Boundary,NULL)==NULL) return(1);

  alpha[0]=2.5; alpha[1]=2.5;
  beta[0] =5.0; beta[1] =5.0;
  point[0]=6; point[1]=11; point[2]=12; point[3]=7;
  if (CreateBoundarySegment("bsouth6",1,0,5,NON_PERIODIC,1,point,
                            alpha,beta,bsouth6Boundary,NULL)==NULL) return(1);

  point[0]=7; point[1]=12; point[2]=13; point[3]=8;
  alpha[0]=2.5; alpha[1]=5.0;
  beta[0] =5.0; beta[1] =7.5;
  if (CreateBoundarySegment("bsouth7",1,0,6,NON_PERIODIC,1,point,
                            alpha,beta,bsouth7Boundary,NULL)==NULL) return(1);

  alpha[0]=2.5; alpha[1]=7.5;
  beta[0] =5.0; beta[1] =10.0;
  point[0]=8; point[1]=13; point[2]=14; point[3]=9;
  if (CreateBoundarySegment("bsouth8",1,0,7,NON_PERIODIC,1,point,
                            alpha,beta,bsouth8Boundary,NULL)==NULL) return(1);

  point[0]=10; point[1]=15; point[2]=16; point[3]=11;
  alpha[0]=5.0; alpha[1]=0.0;
  beta[0] =7.5; beta[1] =2.5;
  if (CreateBoundarySegment("bsouth9",1,0,8,NON_PERIODIC,1,point,
                            alpha,beta,bsouth9Boundary,NULL)==NULL) return(1);

  alpha[0]=5.0; alpha[1]=2.5;
  beta[0] =7.5; beta[1] =5.0;
  point[0]=11; point[1]=16; point[2]=17; point[3]=12;
  if (CreateBoundarySegment("bsouth10",1,0,9,NON_PERIODIC,1,point,
                            alpha,beta,bsouth10Boundary,NULL)==NULL) return(1);

  alpha[0]=5.0; alpha[1]=5.0;
  beta[0] =7.5; beta[1] =7.5;
  point[0]=12; point[1]=17; point[2]=18; point[3]=13;
  if (CreateBoundarySegment("bsouth11",1,0,10,NON_PERIODIC,1,point,
                            alpha,beta,bsouth11Boundary,NULL)==NULL) return(1);

  alpha[0]=5.0; alpha[1]=7.5;
  beta[0] =7.5; beta[1] =10.0;
  point[0]=13; point[1]=18; point[2]=19; point[3]=14;
  if (CreateBoundarySegment("bsouth12",1,0,11,NON_PERIODIC,1,point,
                            alpha,beta,bsouth12Boundary,NULL)==NULL) return(1);

  point[0]=15; point[1]=20; point[2]=21; point[3]=16;
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  if (CreateBoundarySegment("bsouth13",1,0,12,NON_PERIODIC,1,point,
                            alpha,beta,bsouth13Boundary,NULL)==NULL) return(1);

  point[0]=16; point[1]=21; point[2]=22; point[3]=23;
  if (CreateBoundarySegment("bsouth14",1,0,13,NON_PERIODIC,1,point,
                            alpha,beta,bsouth14Boundary,NULL)==NULL) return(1);

  point[0]=16; point[1]=23; point[2]=24; point[3]=17;
  alpha[0]=7.5; alpha[1]=2.5;
  beta[0] =10.0; beta[1] =5.0;
  if (CreateBoundarySegment("bsouth15",1,0,14,NON_PERIODIC,1,point,
                            alpha,beta,bsouth15Boundary,NULL)==NULL) return(1);

  alpha[0]=7.5; alpha[1]=5.0;
  beta[0] =10.0; beta[1] =7.5;
  point[0]=17; point[1]=24; point[2]=25; point[3]=18;
  if (CreateBoundarySegment("bsouth16",1,0,15,NON_PERIODIC,1,point,
                            alpha,beta,bsouth16Boundary,NULL)==NULL) return(1);

  alpha[0]=7.5; alpha[1]=7.5;
  beta[0] =10.0; beta[1] =10.0;
  point[0]=18; point[1]=25; point[2]=26; point[3]=19;
  if (CreateBoundarySegment("bsouth17",1,0,16,NON_PERIODIC,1,point,
                            alpha,beta,bsouth17Boundary,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =2.5;
  point[0]=27; point[1]=32; point[2]=33; point[3]=28;
  if (CreateBoundarySegment("bnorth1",0,1,17,NON_PERIODIC,1,point,
                            alpha,beta,bnorth1Boundary,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=2.5;
  beta[0] =2.5; beta[1] =5.0;
  point[0]=28; point[1]=33; point[2]=34; point[3]=29;
  if (CreateBoundarySegment("bnorth2", 0,1,18,NON_PERIODIC,1,point,
                            alpha,beta,bnorth2Boundary, NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=5.0;
  beta[0] =2.5; beta[1] =7.5;
  point[0]=29; point[1]=34; point[2]=35; point[3]=30;
  if (CreateBoundarySegment("bnorth3",0,1,19,NON_PERIODIC,1,point,
                            alpha,beta,bnorth3Boundary,NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=7.5;
  beta[0] =2.5; beta[1] =10.0;
  point[0]=30; point[1]=35; point[2]=36; point[3]=31;
  if (CreateBoundarySegment("bnorth4",0,1,20,NON_PERIODIC,1,point,
                            alpha,beta,bnorth4Boundary,NULL)==NULL) return(1);

  alpha[0]=2.5; alpha[1]=0.0;
  beta[0] =5.0; beta[1] =2.5;
  point[0]=32; point[1]=37; point[2]=38; point[3]=33;
  if (CreateBoundarySegment("bnorth5",0,1,21,NON_PERIODIC,1,point,
                            alpha,beta,bnorth5Boundary,NULL)==NULL) return(1);

  alpha[0]=2.5; alpha[1]=2.5;
  beta[0] =5.0; beta[1] =5.0;
  point[0]=33; point[1]=38; point[2]=39; point[3]=34;
  if (CreateBoundarySegment("bnorth6",0,1,22,NON_PERIODIC,1,point,
                            alpha,beta,bnorth6Boundary,NULL)==NULL) return(1);

  point[0]=34; point[1]=39; point[2]=40; point[3]=35;
  alpha[0]=2.5; alpha[1]=5.0;
  beta[0] =5.0; beta[1] =7.5;
  if (CreateBoundarySegment("bnorth7",0,1,23,NON_PERIODIC,1,point,
                            alpha,beta,bnorth7Boundary,NULL)==NULL) return(1);

  alpha[0]=2.5; alpha[1]=7.5;
  beta[0] =5.0; beta[1] =10.0;
  point[0]=35; point[1]=40; point[2]=41; point[3]=36;
  if (CreateBoundarySegment("bnorth8",0,1,24,NON_PERIODIC,1,point,
                            alpha,beta,bnorth8Boundary,NULL)==NULL) return(1);

  point[0]=37; point[1]=42; point[2]=43; point[3]=38;
  alpha[0]=5.0; alpha[1]=0.0;
  beta[0] =7.5; beta[1] =2.5;
  if (CreateBoundarySegment("bnorth9",0,1,25,NON_PERIODIC,1,point,
                            alpha,beta,bnorth9Boundary,NULL)==NULL) return(1);

  alpha[0]=5.0; alpha[1]=2.5;
  beta[0] =7.5; beta[1] =5.0;
  point[0]=38; point[1]=43; point[2]=44; point[3]=39;
  if (CreateBoundarySegment("bnorth10",0,1,26,NON_PERIODIC,1,point,
                            alpha,beta,bnorth10Boundary,NULL)==NULL) return(1);

  alpha[0]=5.0; alpha[1]=5.0;
  beta[0] =7.5; beta[1] =7.5;
  point[0]=39; point[1]=44; point[2]=45; point[3]=40;
  if (CreateBoundarySegment("bnorth11",0,1,27,NON_PERIODIC,1,point,
                            alpha,beta,bnorth11Boundary,NULL)==NULL) return(1);

  alpha[0]=5.0; alpha[1]=7.5;
  beta[0] =7.5; beta[1] =10.0;
  point[0]=40; point[1]=45; point[2]=46; point[3]=41;
  if (CreateBoundarySegment("bnorth12",0,1,28,NON_PERIODIC,1,point,
                            alpha,beta,bnorth12Boundary,NULL)==NULL) return(1);

  point[0]=42; point[1]=47; point[2]=48; point[3]=43;
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  if (CreateBoundarySegment("bnorth13",0,1,29,NON_PERIODIC,1,point,
                            alpha,beta,bnorth13Boundary,NULL)==NULL) return(1);

  point[0]=43; point[1]=48; point[2]=49; point[3]=50;
  if (CreateBoundarySegment("bnorth14",0,1,30,NON_PERIODIC,1,point,
                            alpha,beta,bnorth14Boundary,NULL)==NULL) return(1);

  point[0]=43; point[1]=50; point[2]=51; point[3]=44;
  alpha[0]=7.5; alpha[1]=2.5;
  beta[0] =10.0; beta[1] =5.0;
  if (CreateBoundarySegment("bnorth15",0,1,31,NON_PERIODIC,1,point,
                            alpha,beta,bnorth15Boundary,NULL)==NULL) return(1);

  alpha[0]=7.5; alpha[1]=5.0;
  beta[0] =10.0; beta[1] =7.5;
  point[0]=44; point[1]=51; point[2]=52; point[3]=45;
  if (CreateBoundarySegment("bnorth16",0,1,32,NON_PERIODIC,1,point,
                            alpha,beta,bnorth16Boundary,NULL)==NULL) return(1);

  alpha[0]=7.5; alpha[1]=7.5;
  beta[0] =10.0; beta[1] =10.0;
  point[0]=45; point[1]=52; point[2]=53; point[3]=46;
  if (CreateBoundarySegment("bnorth17",0,1,33,NON_PERIODIC,1,point,
                            alpha,beta,bnorth17Boundary,NULL)==NULL) return(1);


  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =BENCH_H;
  point[0]=0; point[1]=5; point[2]=32; point[3]=27;
  if (CreateBoundarySegment("bfront1", 0,1,34,NON_PERIODIC,1,point,alpha,beta,
                            bfront1Boundary, NULL)==NULL) return(1);

  point[0]=5; point[1]=10; point[2]=37; point[3]=32;
  alpha[0]=2.5; alpha[1]=0.0;
  beta[0] =5.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bfront2", 0,1,35,NON_PERIODIC,1,point,alpha,beta,
                            bfront2Boundary, NULL)==NULL) return(1);

  point[0]=10; point[1]=15; point[2]=42; point[3]=37;
  alpha[0]=5.0; alpha[1]=0.0;
  beta[0] =7.5; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bfront3", 0,1,36,NON_PERIODIC,1,point,alpha,beta,
                            bfront3Boundary, NULL)==NULL) return(1);

  alpha[0]=7.5; alpha[1]=0.0;
  beta[0] =9.0; beta[1] =BENCH_H;
  point[0]=15; point[1]=20; point[2]=47; point[3]=42;
  if (CreateBoundarySegment("bfront4", 0,1,37,NON_PERIODIC,1,point,alpha,beta,
                            bfront4Boundary, NULL)==NULL) return(1);

  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =BENCH_H;
  point[0]=20; point[1]=21; point[2]=48; point[3]=47;
  if (CreateBoundarySegment("beast1", 0,1,38,NON_PERIODIC,1,point,alpha,beta,
                            beast1Boundary, NULL)==NULL) return(1);

  point[0]=21; point[1]=22; point[2]=49; point[3]=48;
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("beast2", 0,1,39,NON_PERIODIC,1,point,alpha,beta,
                            beast2Boundary, NULL)==NULL) return(1);

  alpha[0]=1.0; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =BENCH_H;
  point[0]=22; point[1]=23; point[2]=50; point[3]=49;
  if (CreateBoundarySegment("beast3", 0,1,40,NON_PERIODIC,1,point,alpha,beta,
                            beast3Boundary, NULL)==NULL) return(1);

  point[0]=23; point[1]=24; point[2]=51; point[3]=50;
  alpha[0]=2.5; alpha[1]=0.0;
  beta[0] =5.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("beast4", 0,1,41,NON_PERIODIC,1,point,alpha,beta,
                            beast4Boundary, NULL)==NULL) return(1);

  point[0]=24; point[1]=25; point[2]=52; point[3]=51;
  alpha[0]=5.0; alpha[1]=0.0;
  beta[0] =7.5; beta[1] =BENCH_H;
  if (CreateBoundarySegment("beast5", 0,1,42,NON_PERIODIC,1,point,alpha,beta,
                            beast5Boundary, NULL)==NULL) return(1);

  point[0]=25; point[1]=26; point[2]=53; point[3]=52;
  alpha[0]=7.5; alpha[1]=0.0;
  beta[0] =10.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("beast6", 0,1,43,NON_PERIODIC,1,point,alpha,beta,
                            beast6Boundary, NULL)==NULL) return(1);

  point[0]=4; point[1]=9; point[2]=36; point[3]=31;
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bback1", 1,0,44,NON_PERIODIC,1,point,alpha,beta,
                            bback1Boundary, NULL)==NULL) return(1);

  point[0]=9; point[1]=14; point[2]=41; point[3]=36;
  alpha[0]=2.5; alpha[1]=0.0;
  beta[0] =5.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bback2", 1,0,45,NON_PERIODIC,1,point,alpha,beta,
                            bback2Boundary, NULL)==NULL) return(1);

  point[0]=14; point[1]=19; point[2]=46; point[3]=41;
  alpha[0]=5.0; alpha[1]=0.0;
  beta[0] =7.5; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bback3", 1,0,46,NON_PERIODIC,1,point,alpha,beta,
                            bback3Boundary, NULL)==NULL) return(1);

  point[0]=19; point[1]=26; point[2]=53; point[3]=46;
  alpha[0]=7.5; alpha[1]=0.0;
  beta[0] =10.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bback4", 1,0,47,NON_PERIODIC,1,point,alpha,beta,
                            bback4Boundary, NULL)==NULL) return(1);

  point[0]=0; point[1]=1; point[2]=28; point[3]=27;
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =2.5; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bwest1", 1,0,48,NON_PERIODIC,1,point,alpha,beta,
                            bwest1Boundary, NULL)==NULL) return(1);

  point[0]=1; point[1]=2; point[2]=29; point[3]=28;
  alpha[0]=2.5; alpha[1]=0.0;
  beta[0] =5.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bwest2", 1,0,49,NON_PERIODIC,1,point,alpha,beta,
                            bwest2Boundary, NULL)==NULL) return(1);

  point[0]=2; point[1]=3; point[2]=30; point[3]=29;
  alpha[0]=5.0; alpha[1]=0.0;
  beta[0] =7.5; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bwest3", 1,0,50,NON_PERIODIC,1,point,alpha,beta,
                            bwest3Boundary, NULL)==NULL) return(1);

  point[0]=3; point[1]=4; point[2]=31; point[3]=30;
  alpha[0]=7.5; alpha[1]=0.0;
  beta[0] =10.0; beta[1] =BENCH_H;
  if (CreateBoundarySegment("bwest4", 1,0,51,NON_PERIODIC,1,point,alpha,beta,
                            bwest4Boundary, NULL)==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the hole domain                                                  */
/*                                                                          */
/****************************************************************************/

#define RADIUS1 0.2

static INT south_Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 0.0;

  /* return ok */
  return(0);
}

static INT east_Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] = 1.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT north_Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = lambda2;
  result[2] = 1.0;

  /* return ok */
  return(0);
}

static INT west_Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] = 0.0;
  result[1] = lambda1;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT front_Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  result[0] = lambda1;
  result[1] = 0.0;
  result[2] = lambda2;

  /* return ok */
  return(0);
}

static INT back_Boundary (void *data, DOUBLE *param, DOUBLE *result)
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

  /* return ok */
  return(0);
}


static INT south_Boundary1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,lambda2,0.0,RADIUS1,result);

  /* return ok */
  return(0);
}

static INT east_Boundary1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (1.0,lambda1,lambda2,RADIUS1,result);

  /* return ok */
  return(0);
}

static INT north_Boundary1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,lambda2,1.0,RADIUS1,result);

  /* return ok */
  return(0);
}

static INT west_Boundary1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (0.0,lambda1,lambda2,RADIUS1,result);

  /* return ok */
  return(0);
}

static INT front_Boundary1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,0.0,lambda2,RADIUS1,result);

  /* return ok */
  return(0);
}

static INT back_Boundary1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda1,lambda2;

  lambda1 = param[0];
  lambda2 = param[1];

  /* check range */
  if ( lambda1<0.0 || lambda1>1.0 || lambda2<0.0 || lambda2>1.0 ) return(1);

  /* fill result */
  ProjectOnBall (lambda1,1.0,lambda2,RADIUS1,result);

  /* return ok */
  return(0);
}

static INT InitHole (void)
{
  INT point[CORNERS_OF_BND_SEG];
  DOUBLE radius,MidPoint[3], alpha[DIM_OF_BND], beta[DIM_OF_BND];

  /* allocate new domain structure */
  MidPoint[0] = 0.5;
  MidPoint[1] = 0.5;
  MidPoint[2] = 0.5;
  radius = 1.0;
  if (CreateDomain("Hole",MidPoint,radius,12,16,NO)==NULL) return(1);

  /* allocate the boundary segments */
  alpha[0]=0.0; alpha[1]=0.0;
  beta[0] =1.0; beta[1] =1.0;
  point[0]=0; point[1]=1; point[2]=2; point[3]=3;
  if (CreateBoundarySegment("loch0",1,0,0,NON_PERIODIC,1,point,alpha,beta,south_Boundary,NULL)==NULL) return(1);
  point[0]=0; point[1]=3; point[2]=7; point[3]=4;
  if (CreateBoundarySegment("loch1", 1,0,1,NON_PERIODIC,1,point,alpha,beta,west_Boundary, NULL)==NULL) return(1);
  point[0]=0; point[1]=1; point[2]=5; point[3]=4;
  if (CreateBoundarySegment("loch2",0,1,2,NON_PERIODIC,1,point,alpha,beta,front_Boundary,NULL)==NULL) return(1);
  point[0]=4; point[1]=5; point[2]=6; point[3]=7;
  if (CreateBoundarySegment("loch3",0,1,3,NON_PERIODIC,1,point,alpha,beta,north_Boundary,NULL)==NULL) return(1);
  point[0]=1; point[1]=2; point[2]=6; point[3]=5;
  if (CreateBoundarySegment("loch4", 0,1,4,NON_PERIODIC,1,point,alpha,beta,east_Boundary, NULL)==NULL) return(1);
  point[0]=3; point[1]=2; point[2]=6; point[3]=7;
  if (CreateBoundarySegment("loch5", 1,0,5,NON_PERIODIC,1,point,alpha,beta,back_Boundary, NULL)==NULL) return(1);

  point[0]=8; point[1]=9; point[2]=10; point[3]=11;
  if (CreateBoundarySegment("loch6",0,1,6,NON_PERIODIC,1,point,alpha,beta,south_Boundary1,NULL)==NULL) return(1);
  point[0]=8; point[1]=11; point[2]=15; point[3]=12;
  if (CreateBoundarySegment("loch7", 0,1,7,NON_PERIODIC,1,point,alpha,beta,west_Boundary1, NULL)==NULL) return(1);
  point[0]=8; point[1]=9; point[2]=13; point[3]=12;
  if (CreateBoundarySegment("loch8",0,1,8,NON_PERIODIC,1,point,alpha,beta,front_Boundary1,NULL)==NULL) return(1);
  point[0]=12; point[1]=13; point[2]=14; point[3]=15;
  if (CreateBoundarySegment("loch9",0,1,9,NON_PERIODIC,1,point,alpha,beta,north_Boundary1,NULL)==NULL) return(1);
  point[0]=9; point[1]=10; point[2]=14; point[3]=13;
  if (CreateBoundarySegment("loch10", 0,1,10,NON_PERIODIC,1,point,alpha,beta,east_Boundary1, NULL)==NULL) return(1);
  point[0]=11; point[1]=10; point[2]=14; point[3]=15;
  if (CreateBoundarySegment("loch11", 0,1,11,NON_PERIODIC,1,point,alpha,beta,back_Boundary1, NULL)==NULL) return(1);
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
    else if (strcmp(DomainName,"Cylinder") == 0)
    {
      if (InitCylinder())
        return(1);
    }
    else if (strcmp(DomainName,"Benchmark") == 0)
    {
      if (InitBenchmark())
        return(1);
    }
    else if (strcmp(DomainName,"Hole") == 0)
    {
      if (InitHole())
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
