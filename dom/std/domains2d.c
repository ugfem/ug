// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      domains2d.c                                                   */
/*                                                                          */
/* Purpose:   domain definitions                                                */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
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
#include "scan.h"

/* dev modules */
#include "ugdevices.h"

/* domain module */
#include "std_domain.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

typedef DOUBLE DOUBLE_VECTOR[DIM];

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

static DOUBLE_VECTOR x_quad[9];
static DOUBLE alpha,left,top,rad1,L,D,glob_h,form;

static DOUBLE Rand[54][2] = {
  {189,22.5},
  {196.915521148934573,18.006474101645949},
  {200.842627739465058,16.578435341453044},
  {208.333333333333343,19.666666666666668},
  {211.666666666666657,22.333333333333332},
  {209,27.5},
  {198.844284782490547,32.257891016875057},
  {198.373683304990749,36.963905791873074},
  {201.5,43},
  {192,49.5},
  {172.228192115054185,54.442175586914956},
  {161.906544346471549,57.844425265434694},
  {152.388503946133994,63.090244788751953},
  {141.1782032302755,66.468911086754474},
  {136.968911086754474,64.8217967697245},
  {132.5,65.5},
  {128.968911086754474,68.321796769724486},
  {127.3217967697245,72.531088913245526},
  {117,83},
  {99.25,95.75},
  {93.175138136185083,96.47385778767196},
  {89.113077615715099,95.796847700926975},
  {79.599336232946285,99.334042940108191},
  {75.46883503526341,104.152961004071557},
  {69,113.5},
  {58.670424705746314,117.890407897836639},
  {53.091134641144684,121.009361569373624},
  {49.017209870606372,125.934753593786596},
  {52.38112061936242,133.540635580472667},
  {53.853873546595345,138.940729646993418},
  {43,146},
  {31,154},
  {17.952268353586355,150.785944133371856},
  {13.613767716931292,140.136338006452434},
  {13.630509212754852,128.636928536610355},
  {26.620407208445954,103.338620868347988},
  {37.124876377588805,89.962388739829962},
  {49.324527587088951,78.111820447633519},
  {79.805639735655504,73.822280949530011},
  {94.581750807323075,75.015715327395981},
  {106.829739109061194,73.414390858952558},
  {116.643986820153827,70.320696897596306},
  {124.356858456271823,66.581768347010197},
  {130.355557147568618,62.711575606862581},
  {133.888822151335773,50.97118227233117},
  {133.854543533457388,42.865908970841083},
  {134.897772300798522,34.827979880927224},
  {141.5,23.5},
  {158.292290881968853,15.255852308602602},
  {167.64644902215224,15.310713580865084},
  {174.033287865701823,17.081497911805641},
  {178.241954993803603,19.153901560940003},
  {182.283095189484527,25.63238075793813},
  {189,22.5}
};

static DOUBLE kiel[83][2] = {
  {1.,.46},
  {1.,.7},
  {1.,1.},
  {.75,1.},
  {.5,1.},
  {.25,1.},
  {0.,1.},
  {0.,.87},
  {.03,.88},
  {.095,.875},
  {.13,.835},
  {.17,.805},
  {.24,.745},
  {.265,.740},
  {.270,.725},
  {.255,.710},
  {.245,.685},
  {.220,.675},
  {.221,.650},
  {.219,.635},
  {.218,.600},
  {.205,.550},
  {.230,.520},
  {.245,.480},
  {.240,.460},
  {.225,.400},
  {.220,.400},
  {.140,.370},
  {.145,.355},
  {.165,.340},
  {.145,.325},
  {.140,.305},
  {.100,.310},
  {.090,.270},
  {.120,.255},
  {.115,.245},
  {.095,.250},
  {.090,.205},
  {.140,.170},
  {.145,.145},
  {.125,.125},
  {.110,.090},
  {.080,.065},
  {.075,.045},
  {.070,.020},
  {.085,.000},
  {.105,.020},
  {.115,.040},
  {.125,.030},
  {.140,.045},
  {.135,.050},
  {.155,.045},
  {.170,.060},
  {.175,.080},
  {.215,.085},
  {.185,.120},
  {.190,.135},
  {.195,.150},
  {.200,.155},
  {.210,.165},
  {.205,.195},
  {.245,.205},
  {.255,.215},
  {.250,.265},
  {.255,.285},
  {.280,.295},
  {.330,.465},
  {.335,.485},
  {.365,.525},
  {.405,.535},
  {.440,.535},
  {.500,.535},
  {.515,.535},
  {.520,.560},
  {.495,.570},
  {.495,.585},
  {.520,.595},
  {.600,.630},
  {.630,.630},
  {.745,.605},
  {.820,.585},
  {.885,.525},
  {1.,.46}
};


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/****************************************************************************/
/*                             domain definitions                           */
/****************************************************************************/
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/*  define a quadrilateral                                                  */
/*                                                                          */
/****************************************************************************/

static INT southBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[0][0] + lambda*x_quad[1][0];
  result[1] = (1.0-lambda)*x_quad[0][1] + lambda*x_quad[1][1];
  if (alpha != 0.0)
    result[1] += alpha*result[0]*result[0]*(1-result[0])*(1-result[0]);

  return(0);
}

static INT eastBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[1][0] + lambda*x_quad[2][0];
  result[1] = (1.0-lambda)*x_quad[1][1] + lambda*x_quad[2][1];

  return(0);
}

static INT northBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[2][0] + lambda*x_quad[3][0];
  result[1] = (1.0-lambda)*x_quad[2][1] + lambda*x_quad[3][1];

  return(0);
}

static INT westBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[0][0] + lambda*x_quad[3][0];
  result[1] = (1.0-lambda)*x_quad[0][1] + lambda*x_quad[3][1];

  return(0);
}

static INT InitQuadrilateral (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = 0.25*(x_quad[0][0]+x_quad[1][0]+x_quad[2][0]+x_quad[3][0]);
  MidPoint[1] = 0.25*(x_quad[0][1]+x_quad[1][1]+x_quad[2][1]+x_quad[3][1]);
  radius =            ABS(x_quad[0][0]-MidPoint[0]);
  radius = MAX(radius,ABS(x_quad[1][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[2][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[3][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[0][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[1][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[2][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[3][1]-MidPoint[1]));
  radius *= 1.42;
  if (CreateDomain("Quadrilateral",MidPoint,radius,4,4,YES)==NULL) return(1);
  if (CreateBoundarySegment2D("south",1,0,0,0,1,1,0.0,1.0,
                              southBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east", 1,0,1,1,2,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north",1,0,2,2,3,1,0.0,1.0,
                              northBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west", 0,1,3,0,3,1,0.0,1.0,
                              westBoundary, NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define two quadrilaterals                                               */
/*                                                                          */
/****************************************************************************/

static INT south2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[1][0] + lambda*x_quad[4][0];
  result[1] = (1.0-lambda)*x_quad[1][1] + lambda*x_quad[4][1];

  return(0);
}

static INT east2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[4][0] + lambda*x_quad[5][0];
  result[1] = (1.0-lambda)*x_quad[4][1] + lambda*x_quad[5][1];

  return(0);
}

static INT north2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[2][0] + lambda*x_quad[5][0];
  result[1] = (1.0-lambda)*x_quad[2][1] + lambda*x_quad[5][1];

  return(0);
}

static const INT two_sd2p[3] = {0,0,3};
static const INT two_sg2p[8] = {0,1,0,0,3,3,3,2};
static const INT two_pt2p[8] = {0,1,1,0,3,3,2,2};
static const DOMAIN_PART_INFO two_dpi = {two_sd2p,two_sg2p,two_pt2p};

static INT InitTwo (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = 0.16666666666*(x_quad[0][0]+x_quad[1][0]
                               +x_quad[2][0]+x_quad[3][0]
                               +x_quad[4][0]+x_quad[5][0]);
  MidPoint[1] = 0.16666666666*(x_quad[0][1]+x_quad[1][1]
                               +x_quad[2][1]+x_quad[3][1]
                               +x_quad[4][1]+x_quad[5][1]);
  radius =            ABS(x_quad[0][0]-MidPoint[0]);
  radius = MAX(radius,ABS(x_quad[1][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[2][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[3][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[4][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[5][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[0][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[1][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[2][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[3][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[4][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[5][1]-MidPoint[1]));

  if (CreateDomainWithParts("Two",MidPoint,radius,8,8,YES,3,&two_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("south", 1,0,0,0,1,1,0.0,1.0,
                              southBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east",  1,2,1,1,2,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north", 1,0,2,2,3,1,0.0,1.0,
                              northBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west",  0,1,3,0,3,1,0.0,1.0,
                              westBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south2",2,0,4,6,4,1,0.0,1.0,
                              south2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east2", 2,0,5,4,5,1,0.0,1.0,
                              east2Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north2",0,2,6,7,5,1,0.0,1.0,
                              north2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east1", 1,2,7,6,7,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define a triangle                                                       */
/*                                                                          */
/****************************************************************************/

static INT diagonal (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[0][0] + lambda*x_quad[2][0];
  result[1] = (1.0-lambda)*x_quad[0][1] + lambda*x_quad[2][1];

  return(0);
}

static INT InitTriangle (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = 0.333333*(x_quad[0][0]+x_quad[1][0]+x_quad[2][0]);
  MidPoint[1] = 0.333333*(x_quad[0][1]+x_quad[1][1]+x_quad[2][1]);
  radius =            ABS(x_quad[0][0]-MidPoint[0]);
  radius = MAX(radius,ABS(x_quad[1][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[2][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[0][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[1][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[2][1]-MidPoint[1]));
  if (CreateDomain("Triangle",MidPoint,radius,3,3,YES)==NULL) return(1);
  if (CreateBoundarySegment2D("south",   1,0,0,0,1,1,0.0,1.0,
                              southBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east",    1,0,1,1,2,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("diagonal",0,1,2,0,2,1,0.0,1.0,
                              diagonal,NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the puctured disc                                                */
/*                                                                          */
/****************************************************************************/

static INT BottomBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[0][0] + lambda*9.0;
  result[1] = (1.0-lambda)*x_quad[0][1] + lambda*0.0;

  return(0);
}

static INT CircleBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>0.5)) return(1);

  /* fill result */
  result[0] = 10.0 + cos(PI*(1.0-lambda));       /* x */ /* PI defined in misc.h */
  result[1] = sin(PI*lambda);       /* y */

  /* return ok */
  return(0);
}

static INT RightBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<1.0)||(lambda>10.0)) return(1);

  /* fill result */
  result[0] = 10.0;
  result[1] = lambda;

  /* return ok */
  return(0);
}

static INT TopBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>10.0)) return(1);

  /* fill result */
  result[0] = 10.0-lambda;
  result[1] = 10.0;

  /* return ok */
  return(0);
}

static INT LeftBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[0][0] + lambda* 0.0;
  result[1] = (1.0-lambda)*x_quad[0][1] + lambda*10.0;

  return(0);
}

static INT InitPuncturedDisc (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = MidPoint[1] = 5.0;
  radius = 7.5;
  if (CreateDomain("Punctured Disc",
                   MidPoint,radius,5,5,NO)==NULL) return(1);

  /* allocate the boundary segments */
  if (CreateBoundarySegment2D("bottom",1,0,
                              0,0,1,1,
                              0.0,1.0,
                              BottomBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("circle", 1,0,
                              1,1,2,20,
                              0.0,0.5,
                              CircleBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right",1,0,
                              2,2,3,1,
                              1.0,10.0,
                              RightBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("top", 1,0,
                              3,3,4,1,
                              0.0,10.0,
                              TopBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left", 0,1,
                              4,0,4,1,
                              0.0,1.0,
                              LeftBoundary, NULL)==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define a variable disc                                                  */
/*                                                                          */
/****************************************************************************/

static INT BVar1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (top - rad1)*lambda;
  result[1] = 0.0;

  return(0);
}

static INT BVar2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>0.5)) return(1);
  result[0] = top  + rad1*cos(PI*(1.0-lambda));
  result[1] = rad1*sin(PI*lambda);

  /* straight for error bounds */
  result[0] = top  - rad1;
  result[1] = rad1*2*lambda;

  return(0);
}

static INT BVar2aBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>0.5)) return(1);

  /* straight for error bounds */
  result[0] = top  + rad1*(2*lambda-1.0);
  result[1] = rad1;

  return(0);
}

static INT BVar3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = top;
  result[1] = rad1 + (left - rad1)*lambda;

  return(0);
}

static INT BVar4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = top - top*lambda;
  result[1] = left;

  return(0);
}

static INT BVar5Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0.0;
  result[1] = left - left*lambda;

  return(0);
}

static INT InitVDisc (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 0.5*top;
  MidPoint[1] = 0.5*left;
  radius = 0.5*sqrt((top*top)+(left*left)) + 0.1;
  if (CreateDomain("Variable Disc",
                   MidPoint,radius,6,6,NO)==NULL) return(1);

  /* allocate the boundary segments */
  if (CreateBoundarySegment2D("bvar1",1,0,
                              0,0,1,1,
                              0.0,1.0,
                              BVar1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bvar2", 1,0,
                              1,1,2,20,
                              0.0,0.5,
                              BVar2Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bvar3",1,0,
                              2,3,4,1,
                              0.0,1.0,
                              BVar3Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bvar4", 1,0,
                              3,4,5,1,
                              0.0,1.0,
                              BVar4Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bvar5", 1,0,
                              4,5,0,1,
                              0.0,1.0,
                              BVar5Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bvar2a", 1,0,
                              5,2,3,20,
                              0.0,0.5,
                              BVar2aBoundary, NULL)==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the CT disc                                                      */
/*                                                                          */
/****************************************************************************/

static INT Bottom1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>10.0)) return(1);
  result[0] = lambda;
  result[1] = 0.5;

  return(0);
}

static INT Bottom2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>0.5)) return(1);
  result[0] = 10.0 + 0.5*sin(PI*lambda);       /* x */ /* PI defined in misc.h */
  result[1] = 0.5*cos(PI*lambda);       /* y */

  return(0);
}

static INT Bottom3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>9.5)) return(1);
  result[0] = 10.5 + lambda;
  result[1] = 0.0;

  return(0);
}

static INT Top1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>20.0)) return(1);
  result[0] = 20.0-lambda;
  result[1] = 12.0;

  return(0);
}

static INT Right1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>12.0)) return(1);
  result[0] = 20.0;
  result[1] = lambda;

  return(0);
}

static INT Left1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>11.5)) return(1);
  result[0] = 0.0;
  result[1] = 12.0 - lambda;

  return(0);
}

static INT CircleUpperBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 6.0 + 0.5*sin(PI*lambda);       /* x */ /* PI defined in misc.h */
  result[1] = 6.5 + 0.5*cos(PI*lambda);       /* y */

  return(0);
}

static INT CircleLowerBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<1.0)||(lambda>2.0)) return(1);
  result[0] = 6.0 + 0.5*sin(PI*lambda);       /* x */ /* PI defined in misc.h */
  result[1] = 6.5 + 0.5*cos(PI*lambda);       /* y */

  return(0);
}

static INT InitCTDisc (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = 10; MidPoint[1] = 6.0;
  radius = 10;
  if (CreateDomain("CT Disc",
                   MidPoint,radius,8,8,NO)==NULL) return(1);
  if (CreateBoundarySegment2D("bottom1",1,0,
                              0,0,1,1,
                              0.0,10.0,
                              Bottom1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bottom2", 1,0,
                              1,1,2,20,
                              0.0,0.5,
                              Bottom2Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bottom3",1,0,
                              2,2,3,1,
                              0.0,9.5,
                              Bottom3Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right1",1,0,
                              3,3,4,1,
                              0.0,12.0,
                              Right1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("top1", 1,0,
                              4,4,5,1,
                              0.0,20.0,
                              Top1Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left1", 1,0,
                              5,5,0,1,
                              0.0,11.5,
                              Left1Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("upper", 1,0,
                              6,6,7,20,
                              0.0,1.0,
                              CircleUpperBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("lower", 1,0,
                              7,7,6,20,
                              1.0,2.0,
                              CircleLowerBoundary, NULL)==NULL) return(1);
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define a circle segment                                                 */
/*                                                                          */
/****************************************************************************/

/*
 * A domain is defined by its boundary, the boundary is composed of
 * any number of boundary segments. Each segment is a mapping from some
 * interval [a,b] to R^2 in two dimensions. In 2D also internal boundaries
 * are allowed to create different subdomains. The points where two boundary
 * segments are joined are called corners. All corners and subdomains are
 * numbered, the exterior of the domain has always number 0.
 */

/* circleBoundary maps [0.0,1.0] to the upper half circle */
static INT circleBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  /* retrieve parameter value */
  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.5)) return(1);

  /* fill result */
  result[0] = cos(PI*lambda);       /* x */ /* PI defined in misc.h */
  result[1] = sin(PI*lambda);       /* y */

  /* return ok */
  return(0);
}

/* circleBoundary maps [0.0,1.0] to the lower half circle */
static INT circleBoundaryHorizontal (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  /* retrieve parameter value */
  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] = lambda;
  result[1] = 0;

  /* return ok */
  return(0);
}

/* circleBoundary maps [0.0,1.0] to the lower half circle */
static INT circleBoundaryVertical (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  /* retrieve parameter value */
  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] = 0;
  result[1] = lambda-1;

  /* return ok */
  return(0);
}

static INT InitSegment (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = MidPoint[1] = 0.0;
  radius = 1.05;
  if (CreateDomain(
        "Segment",                              /* name of the new domain                               */
        MidPoint,radius,                /* circle containing the domain			*/
        3,                                                              /* number of boundary segments                  */
        3,                                                              /* number of corners					*/
        NO                                                              /* true if domain is convex				*/
        )==NULL) return(1);

  /* allocate the boundary segments, segment allocation must                            */
  /* immediately follow the domain definition.							*/
  if (CreateBoundarySegment2D(
        "circle bnd",                                   /* name of the boundary segment                 */
        1,                                                              /* number of left subdomain				*/
        0,                                                              /* number of right subdomain			*/
        0,                                                              /* number of segment, starting with 0	*/
        0,                                                              /* number of corner where segment starts*/
        1,                                                              /* number of corner where segment ends  */
        20,                                                     /* resolution, use 1 for straight line  */
        0.0,                                                    /* begin of parameter interval			*/
        1.5,                                                    /* end of parameter interval			*/
        circleBoundary,                                 /* function mapping parameter to world  */
        NULL                                                    /* user defined pointer to be supplied  */
        )==NULL) return(1);
  if (CreateBoundarySegment2D(
        "circle bnd vert",                              /* name of the boundary segment                 */
        1,                                                              /* number of left subdomain				*/
        0,                                                              /* number of right subdomain			*/
        1,                                                              /* number of segment, starting with 0	*/
        1,                                                              /* number of corner where segment starts*/
        2,                                                              /* number of corner where segment ends  */
        1,                                                              /* resolution, use 1 for straight line  */
        0.0,                                                    /* begin of parameter interval			*/
        1.0,                                                    /* end of parameter interval			*/
        circleBoundaryVertical,                 /* function mapping parameter to world  */
        NULL                                                    /* user defined pointer to be supplied  */
        )==NULL) return(1);
  if (CreateBoundarySegment2D(
        "circle bnd hor",                               /* name of the boundary segment                 */
        1,                                                              /* number of left subdomain				*/
        0,                                                              /* number of right subdomain			*/
        2,                                                              /* number of segment, starting with 0	*/
        2,                                                              /* number of corner where segment starts*/
        0,                                                              /* number of corner where segment ends  */
        1,                                                              /* resolution, use 1 for straight line  */
        0.0,                                                    /* begin of parameter interval			*/
        1.0,                                                    /* end of parameter interval			*/
        circleBoundaryHorizontal,                /* function mapping parameter to world */
        NULL                                                    /* user defined pointer to be supplied  */
        )==NULL) return(1);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the unit circle                                                  */
/*                                                                          */
/****************************************************************************/

static INT circleBoundaryUpper (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = cos(PI*lambda);
  result[1] = sin(PI*lambda);

  return(0);
}

static INT circleBoundaryLower (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = cos(PI+PI*lambda);
  result[1] = sin(PI+PI*lambda);

  return(0);
}

static INT InitCircle (void)
{
  DOUBLE_VECTOR MidPoint;
  DOUBLE radius;

  MidPoint[0] = MidPoint[1] = 0.0;
  radius = 1.05;
  if (CreateDomain("Circle",MidPoint,radius,2,2,YES)==NULL)
    return(1);

  if (CreateBoundarySegment2D("circle bnd upper",
                              1,0,0,0,1,20,0.0,1.0,
                              circleBoundaryUpper,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("circle bnd lower",
                              1,0,1,1,0,20,0.0,1.0,
                              circleBoundaryLower,NULL)==NULL)
    return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the ring                                                         */
/*                                                                          */
/****************************************************************************/

#define INNER_RADIUS 0.3

static INT circleBoundaryUpper1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = INNER_RADIUS * cos(PI*lambda);
  result[1] = INNER_RADIUS * sin(PI*lambda);

  return(0);
}

static INT circleBoundaryLower1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = INNER_RADIUS * cos(PI+PI*lambda);
  result[1] = INNER_RADIUS * sin(PI+PI*lambda);

  return(0);
}

static INT InitRing (void)
{
  DOUBLE_VECTOR MidPoint;
  DOUBLE radius;

  MidPoint[0] = MidPoint[1] = 0.0;
  radius = 1.05;
  if (CreateDomain("Ring",MidPoint,radius,4,4,YES)==NULL)
    return(1);

  if (CreateBoundarySegment2D("ring bnd upper",
                              1,0,0,0,1,20,0.0,1.0,
                              circleBoundaryUpper,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring bnd lower",
                              1,0,1,1,0,20,0.0,1.0,
                              circleBoundaryLower,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring inner bnd upper",
                              0,1,2,2,3,20,0.0,1.0,
                              circleBoundaryUpper1,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring inner bnd lower",
                              0,1,3,3,2,20,0.0,1.0,
                              circleBoundaryLower1,NULL)==NULL)
    return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define Hexagon                                                          */
/*                                                                          */
/****************************************************************************/

static INT southBound (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] = lambda;
  result[1] = 0.0;

  /* return ok */
  return(0);
}
static INT east1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] = 1+lambda*0.5;
  result[1] = lambda*sqrt(3.0)/2.0;

  /* return ok */
  return(0);
}

static INT east2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] = 3.0/2.0-lambda/2.0;
  result[1] = sqrt(3.0)/2.0+sqrt(3.0)/2.0*lambda;

  /* return ok */
  return(0);
}

static INT northBound (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] = 1.0-lambda;
  result[1] = sqrt(3.0);

  /* return ok */
  return(0);
}
static INT west1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] =-lambda*0.5;
  result[1] = sqrt(3.0)-lambda*sqrt(3.0)/2.0;

  /* return ok */
  return(0);
}

static INT west2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>1.0)) return(1);

  /* fill result */
  result[0] = -0.5 +lambda*0.5;
  result[1] = sqrt(3.0)/2.0-lambda*sqrt(3.0)/2.0;

  /* return ok */
  return(0);
}

static INT InitHexagon (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 0.5;
  MidPoint[1] = sqrt(3.0)/2.0;
  radius = 1.0;
  if (CreateDomain("Hexagon",MidPoint,radius,6,6,YES)==NULL) return(1);

  /* allocate the boundary segments */
  if (CreateBoundarySegment2D("south", 1,0,0,0,1,1,0.0,1.0,
                              southBound,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east1", 1,0,1,1,2,1,0.0,1.0,
                              east1, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east2", 1,0,2,2,3,1,0.0,1.0,
                              east2, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north", 1,0,3,3,4,1,0.0,1.0,
                              northBound,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west1", 1,0,4,4,5,1,0.0,1.0,
                              west1, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west2" ,1,0,5,5,0,1,0.0,1.0,
                              west2, NULL)==NULL) return(1);
  /* return ok */
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the Wolfgangsee domain                                           */
/*                                                                          */
/****************************************************************************/

#define MAXSEG 6

static INT unteresWolfgangseeUfer (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;
  DOUBLE c;
  INT i,j;

  lambda = param[0];

  /* check range */
  if ((lambda<0.0)||(lambda>20.0)) return(1);

  /* fill result */
  c=lambda-floor(lambda);
  i=(INT)ceil(lambda);
  j=(INT)floor(lambda);

  result[0] = (1.0-c)*Rand[j][0] + c*Rand[i][0];
  result[1] = (1.0-c)*Rand[j][1] + c*Rand[i][1];

  /* return ok */
  return(0);
}

static INT rechtesWolfgangseeUfer (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;
  DOUBLE c;
  INT i,j;

  lambda = param[0];

  /* check range */
  if ((lambda<20.0)||(lambda>40.0)) return(1);

  /* fill result */
  c=lambda-floor(lambda);
  i=(INT)ceil(lambda);
  j=(INT)floor(lambda);

  result[0] = (1.0-c)*Rand[j][0] + c*Rand[i][0];
  result[1] = (1.0-c)*Rand[j][1] + c*Rand[i][1];

  /* return ok */
  return(0);
}

static INT oberesWolfgangseeUfer (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;
  DOUBLE c;
  INT i,j;

  lambda = param[0];

  /* check range */
  if ((lambda<40.0)||(lambda>53.0)) return(1);

  /* fill result */
  c=lambda-floor(lambda);
  i=(INT)ceil(lambda);
  j=(INT)floor(lambda);

  result[0] = (1.0-c)*Rand[j][0] + c*Rand[i][0];
  result[1] = (1.0-c)*Rand[j][1] + c*Rand[i][1];

  /* return ok */
  return(0);
}

static INT InitWolfgangsee (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] =115.0;
  MidPoint[1] =70.0;
  radius =150.0;
  if (CreateDomain(
        "Wolfgangsee",                                          /* name of the new domain           */
        MidPoint,radius,                                        /* circle containing the domain	*/
        3,                                                                              /* number of boundary segments  */
        3,                                                                              /* number of corners			*/
        NO                                                                              /* true if domain is convex	    */
        )==NULL)
    return(1);

  if (CreateBoundarySegment2D("WolfgangseeUfer0",1,0,0,0,1,20,0.0,20.0,
                              unteresWolfgangseeUfer,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("WolfgangseeUfer1",1,0,1,1,2,20,20.0,40.0,
                              rechtesWolfgangseeUfer,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("WolfgangseeUfer2",1,0,2,2,0,20,40.0,53.0,
                              oberesWolfgangseeUfer,NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define the Kiel domain                                                  */
/*                                                                          */
/****************************************************************************/

static INT kiel_lower (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda = param[0];
  DOUBLE c;
  INT i,j;

  if ((lambda<0.0)||(lambda>20.0)) return(1);
  c=lambda-floor(lambda);
  i=(INT)ceil(lambda);
  j=(INT)floor(lambda);
  result[0] = (1.0-c)*kiel[j][0] + c*kiel[i][0];
  result[1] = (1.0-c)*kiel[j][1] + c*kiel[i][1];

  return(0);
}

static INT kiel_right (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda = param[0];
  DOUBLE c;
  INT i,j;

  if ((lambda<20.0)||(lambda>60.0)) return(1);
  c=lambda-floor(lambda);
  i=(INT)ceil(lambda);
  j=(INT)floor(lambda);
  result[0] = (1.0-c)*kiel[j][0] + c*kiel[i][0];
  result[1] = (1.0-c)*kiel[j][1] + c*kiel[i][1];

  return(0);
}

static INT kiel_upper (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda = param[0];
  DOUBLE c;
  INT i,j;

  if ((lambda<40.0)||(lambda>78.0)) return(1);
  c=lambda-floor(lambda);
  i=(INT)ceil(lambda);
  j=(INT)floor(lambda);
  result[0] = (1.0-c)*kiel[j][0] + c*kiel[i][0];
  result[1] = (1.0-c)*kiel[j][1] + c*kiel[i][1];

  return(0);
}

static INT kiel_left (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda = param[0];
  DOUBLE c;
  INT i,j;

  if ((lambda<60.0)||(lambda>82.0)) return(1);
  c=lambda-floor(lambda);
  i=(INT)ceil(lambda);
  j=(INT)floor(lambda);
  result[0] = (1.0-c)*kiel[j][0] + c*kiel[i][0];
  result[1] = (1.0-c)*kiel[j][1] + c*kiel[i][1];

  return(0);
}

static INT InitKiel (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] =0.5;
  MidPoint[1] =0.5;
  radius =1;
  if (CreateDomain("Kiel",MidPoint,radius,4,4,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("kiel_lower",1,0,0,0,1,20,0.0,20.0,
                              kiel_lower,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("kiel_right",1,0,1,1,2,20,20.0,40.0,
                              kiel_right,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("kiel_upper",1,0,2,2,3,20,40.0,60.0,
                              kiel_upper,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("kiel_left" ,1,0,3,3,0,20,60.0,82.0,
                              kiel_left,NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define two holes                                                        */
/*                                                                          */
/****************************************************************************/

static INT Start1Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = lambda;
  result[1] = 0.0;

  return(0);
}

static INT Start2Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0;
  result[1] = lambda;

  return(0);
}

static INT Start3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0 - lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0.0;
  result[1] = 1.0 - lambda;

  return(0);
}

static INT Start5Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start6Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0;
  result[1] = 1.0 + lambda;

  return(0);
}

static INT Start7Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0 - lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start8Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0.0;
  result[1] = 2.0 - lambda;

  return(0);
}

static INT Start9Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start10Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0;
  result[1] = 2.0 + lambda;

  return(0);
}

static INT Start11Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0 - lambda;
  result[1] = 3.0;

  return(0);
}

static INT Start12Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0.0;
  result[1] = 3.0 - lambda;

  return(0);
}

static INT Start13Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0 + lambda;
  result[1] = 0.0;

  return(0);
}

static INT Start14Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0;
  result[1] = lambda;

  return(0);
}

static INT Start15Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0 - lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start16Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0;
  result[1] = 1.0 - lambda;

  return(0);
}

static INT Start17Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0 + lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start18Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0;
  result[1] = 2.0 + lambda;

  return(0);
}

static INT Start19Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0 - lambda;
  result[1] = 3.0;

  return(0);
}

static INT Start20Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 1.0;
  result[1] = 3.0 - lambda;

  return(0);
}

static INT Start21Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0 + lambda;
  result[1] = 0.0;

  return(0);
}

static INT Start22Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0;
  result[1] = lambda;

  return(0);
}

static INT Start23Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0 - lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start24Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0;
  result[1] = 1.0 - lambda;

  return(0);
}

static INT Start25Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0 + lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start26Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0;
  result[1] = 1.0 + lambda;

  return(0);
}

static INT Start27Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0 - lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start28Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0;
  result[1] = 2.0 - lambda;

  return(0);
}

static INT Start29Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0 + lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start30Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0;
  result[1] = 2.0 + lambda;

  return(0);
}

static INT Start31Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0 - lambda;
  result[1] = 3.0;

  return(0);
}

static INT Start32Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 2.0;
  result[1] = 3.0 - lambda;

  return(0);
}

static INT Start33Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0 + lambda;
  result[1] = 0.0;

  return(0);
}

static INT Start34Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0;
  result[1] = lambda;

  return(0);
}

static INT Start35Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0 -lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start36Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0;
  result[1] = 1.0 - lambda;

  return(0);
}

static INT Start37Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0 + lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start38Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0;
  result[1] = 2.0 + lambda;

  return(0);
}

static INT Start39Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0 - lambda;
  result[1] = 3.0;

  return(0);
}

static INT Start40Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3.0;
  result[1] = 3.0 - lambda;

  return(0);
}

static INT Start41Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0 + lambda;
  result[1] = 0.0;

  return(0);
}

static INT Start42Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 5.0;
  result[1] = lambda;

  return(0);
}

static INT Start43Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 5.0 - lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start44Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0;
  result[1] = 1.0 - lambda;

  return(0);
}

static INT Start45Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0 + lambda;
  result[1] = 1.0;

  return(0);
}

static INT Start46Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 5.0;
  result[1] = 1.0 + lambda;

  return(0);
}

static INT Start47Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 5.0 - lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start48Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0;
  result[1] = 2.0 - lambda;

  return(0);
}

static INT Start49Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0 + lambda;
  result[1] = 2.0;

  return(0);
}

static INT Start50Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 5.0;
  result[1] = 2.0 + lambda;

  return(0);
}

static INT Start51Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 5.0 - lambda;
  result[1] = 3.0;

  return(0);
}

static INT Start52Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 4.0;
  result[1] = 3.0 - lambda;

  return(0);
}

static const INT h_sd2p[2] = {0,0};
static const INT h_sg2p[52] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const INT h_pt2p[52] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const DOMAIN_PART_INFO h_dpi = {h_sd2p,h_sg2p,h_pt2p};

static INT InitHoles (void)
{
  DOUBLE radius,MidPoint[2];
  INT vert[52],i;

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomainWithParts("Holes",MidPoint,radius,52,52,NO,2,&h_dpi)
      ==NULL)
    return(1);

  for (i=0; i<52; i++) vert[i] = i;

  if (CreateBoundarySegment2D("start1" ,1,0,0,vert[0],vert[1],1.0,0.0,1.0,
                              Start1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2" ,1,0,1,vert[1],vert[2],1.0,0.0,1.0,
                              Start2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start3" ,1,0,2,vert[2],vert[3],1.0,0.0,1.0,
                              Start3Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start4" ,1,0,3,vert[3],vert[0],1.0,0.0,1.0,
                              Start4Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5" ,1,0,4,vert[4],vert[5],1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start6" ,1,0,5,vert[5],vert[6],1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start7" ,1,0,6,vert[6],vert[7],1.0,0.0,1.0,
                              Start7Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start8" ,1,0,7,vert[7],vert[4],1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start9" ,1,0,8,vert[8],vert[9],1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start10",1,0,9,vert[9],vert[10],1.0,0.0,1.0,
                              Start10Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start11",1,0,10,vert[10],vert[11],1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start12",1,0,11,vert[11],vert[8],1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start13",1,0,12,vert[12],vert[13],1.0,0.0,1.0,
                              Start13Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start14",1,0,13,vert[13],vert[14],1.0,0.0,1.0,
                              Start14Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start15",1,0,14,vert[14],vert[15],1.0,0.0,1.0,
                              Start15Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start16",1,0,15,vert[15],vert[12],1.0,0.0,1.0,
                              Start16Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start17",1,0,16,vert[16],vert[17],1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start18",1,0,17,vert[17],vert[18],1.0,0.0,1.0,
                              Start18Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start19",1,0,18,vert[18],vert[19],1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start20",1,0,19,vert[19],vert[16],1.0,0.0,1.0,
                              Start20Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start21",1,0,20,vert[20],vert[21],1.0,0.0,1.0,
                              Start21Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start22",1,0,21,vert[21],vert[22],1.0,0.0,1.0,
                              Start22Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start23",1,0,22,vert[22],vert[23],1.0,0.0,1.0,
                              Start23Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start24",1,0,23,vert[23],vert[20],1.0,0.0,1.0,
                              Start24Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start25",1,0,24,vert[24],vert[25],1.0,0.0,1.0,
                              Start25Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start26",1,0,25,vert[25],vert[26],1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start27",1,0,26,vert[26],vert[27],1.0,0.0,1.0,
                              Start27Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start28",1,0,27,vert[27],vert[24],1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start29",1,0,28,vert[28],vert[29],1.0,0.0,1.0,
                              Start29Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start30",1,0,29,vert[29],vert[30],1.0,0.0,1.0,
                              Start30Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start31",1,0,30,vert[30],vert[31],1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start32",1,0,31,vert[31],vert[28],1.0,0.0,1.0,
                              Start32Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start33",1,0,32,vert[32],vert[33],1.0,0.0,1.0,
                              Start33Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start34",1,0,33,vert[33],vert[34],1.0,0.0,1.0,
                              Start34Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start35",1,0,34,vert[34],vert[35],1.0,0.0,1.0,
                              Start35Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start36",1,0,35,vert[35],vert[32],1.0,0.0,1.0,
                              Start36Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start37",1,0,36,vert[36],vert[37],1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start38",1,0,37,vert[37],vert[38],1.0,0.0,1.0,
                              Start38Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start39",1,0,38,vert[38],vert[39],1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start40",1,0,39,vert[39],vert[36],1.0,0.0,1.0,
                              Start40Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start41",1,0,40,vert[40],vert[41],1.0,0.0,1.0,
                              Start41Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start42",1,0,41,vert[41],vert[42],1.0,0.0,1.0,
                              Start42Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start43",1,0,42,vert[42],vert[43],1.0,0.0,1.0,
                              Start43Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start44",1,0,43,vert[43],vert[40],1.0,0.0,1.0,
                              Start44Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start45",1,0,44,vert[44],vert[45],1.0,0.0,1.0,
                              Start45Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start46",1,0,45,vert[45],vert[46],1.0,0.0,1.0,
                              Start46Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start47",1,0,46,vert[46],vert[47],1.0,0.0,1.0,
                              Start47Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start48",1,0,47,vert[47],vert[44],1.0,0.0,1.0,
                              Start48Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start49",1,0,48,vert[48],vert[49],1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start50",1,0,49,vert[49],vert[50],1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start51",1,0,50,vert[50],vert[51],1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start52",1,0,51,vert[51],vert[48],1.0,0.0,1.0,
                              Start52Boundary,NULL)==NULL) return(1);

  return(0);
}

static INT InitHoles2 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Holes2",MidPoint,radius,24,24,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("start2_1" ,1,0,0,0,1,1.0,0.0,1.0,
                              Start1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_2" ,1,0,1,1,2,1.0,0.0,1.0,
                              Start13Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_3" ,1,0,2,2,3,1.0,0.0,1.0,
                              Start21Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_4" ,1,0,3,3,4,1.0,0.0,1.0,
                              Start33Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_5" ,1,0,4,4,5,1.0,0.0,1.0,
                              Start41Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_6" ,1,0,5,5,6,1.0,0.0,1.0,
                              Start42Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_7" ,1,0,6,6,7,1.0,0.0,1.0,
                              Start46Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_8" ,1,0,7,7,8,1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_9" ,1,0,8,8,9,1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_10",1,0,9,9,10,1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_11",1,0,10,10,11,1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_12",1,0,11,11,12,1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_13",1,0,12,12,13,1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_14",1,0,13,13,14,1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_15",1,0,14,14,15,1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_16",1,0,15,15,0,1.0,0.0,1.0,
                              Start4Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_17",1,0,16,16,17,1.0,0.0,1.0,
                              Start15Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_18",1,0,17,17,18,1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_19",1,0,18,18,19,1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_20",1,0,19,19,16,1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_21",1,0,20,20,21,1.0,0.0,1.0,
                              Start35Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_22",1,0,21,21,22,1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_23",1,0,22,22,23,1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_24",1,0,23,23,20,1.0,0.0,1.0,
                              Start48Boundary,NULL)==NULL) return(1);
  return(0);
}

static INT InitHoles3 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Holes3",MidPoint,radius,24,24,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("start2_1" ,1,0,0,0,1,1.0,0.0,1.0,
                              Start1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_2" ,1,0,1,1,2,1.0,0.0,1.0,
                              Start13Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_3" ,1,0,2,2,3,1.0,0.0,1.0,
                              Start21Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_4" ,1,0,3,3,4,1.0,0.0,1.0,
                              Start33Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_5" ,1,0,4,4,5,1.0,0.0,1.0,
                              Start41Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_6" ,1,0,5,5,6,1.0,0.0,1.0,
                              Start42Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_7" ,1,0,6,6,20,1.0,0.0,1.0,
                              Start43Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_8" ,1,0,7,7,8,1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_9" ,1,0,8,8,9,1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_10",1,0,9,9,10,1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_11",1,0,10,10,11,1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_12",1,0,11,11,12,1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_13",1,0,12,12,13,1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_14",1,0,13,13,14,1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_15",1,0,14,14,15,1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_16",1,0,15,15,0,1.0,0.0,1.0,
                              Start4Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_17",1,0,16,16,17,1.0,0.0,1.0,
                              Start15Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_18",1,0,17,17,18,1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_19",1,0,18,18,19,1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_20",1,0,19,19,16,1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_21",1,0,20,20,21,1.0,0.0,1.0,
                              Start35Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_22",1,0,21,21,22,1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_23",1,0,22,22,23,1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_24",0,1,23,7,23,1.0,0.0,1.0,
                              Start47Boundary,NULL)==NULL) return(1);

  return(0);
}

static INT InitHoles4 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Holes4",MidPoint,radius,24,24,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("start2_1" ,1,0,0,0,1,1.0,0.0,1.0,
                              Start1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_2" ,1,0,1,1,2,1.0,0.0,1.0,
                              Start13Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_3" ,1,0,2,2,3,1.0,0.0,1.0,
                              Start21Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_4" ,1,0,3,3,4,1.0,0.0,1.0,
                              Start33Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_5" ,1,0,4,4,5,1.0,0.0,1.0,
                              Start41Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_6" ,1,0,5,5,6,1.0,0.0,1.0,
                              Start42Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_7" ,1,0,6,6,20,1.0,0.0,1.0,
                              Start43Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_8" ,1,0,7,7,8,1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_9" ,1,0,8,8,9,1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_10",1,0,9,9,10,1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_11",1,0,10,10,11,1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_12",1,0,11,11,12,1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_13",1,0,12,12,13,1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_14",1,0,13,13,14,1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_16",1,0,15,15,0,1.0,0.0,1.0,
                              Start4Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_17",1,0,16,16,17,1.0,0.0,1.0,
                              Start15Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_19",1,0,18,18,19,1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_20",1,0,19,19,16,1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_21",1,0,20,20,21,1.0,0.0,1.0,
                              Start35Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_22",1,0,21,21,22,1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_23",1,0,22,22,23,1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_24",0,1,23,7,23,1.0,0.0,1.0,
                              Start47Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_18",1,0,17,17,15,1.0,0.0,1.0,
                              Start3Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2_15",1,0,14,14,18,1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);

  return(0);
}


static INT InitHoles5 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Holes5",MidPoint,radius,20,20,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("start5_1" ,1,0,0,0,1,1.0,0.0,1.0,
                              Start1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_2" ,1,0,1,1,2,1.0,0.0,1.0,
                              Start13Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_3" ,1,0,2,2,3,1.0,0.0,1.0,
                              Start21Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_8" ,1,0,7,7,8,1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_9" ,1,0,8,8,9,1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_10",1,0,9,9,10,1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_11",1,0,10,10,11,1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_12",1,0,11,11,12,1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_13",1,0,12,12,13,1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_14",1,0,13,13,14,1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_15",1,0,14,14,18,1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_16",1,0,15,15,0,1.0,0.0,1.0,
                              Start4Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_17",1,0,16,16,17,1.0,0.0,1.0,
                              Start15Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_18",1,0,17,17,15,1.0,0.0,1.0,
                              Start3Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_19",1,0,18,18,19,1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_20",1,0,19,19,16,1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("start5_4" ,1,0,3,3,4,1.0,0.0,1.0,
                              Start22Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_5" ,1,0,4,4,5,1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_6" ,1,0,5,5,6,1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_7" ,1,0,6,6,7,1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);

  return(0);
}

static INT InitHoles6 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Holes6",MidPoint,radius,16,16,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("start5_3" ,1,0,2,2,3,1.0,0.0,1.0,
                              Start21Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_8" ,1,0,7,7,8,1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_9" ,1,0,8,8,9,1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_10",1,0,9,9,10,1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_11",1,0,10,10,11,1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_12",1,0,11,11,12,1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_13",1,0,12,12,13,1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_14",1,0,13,13,14,1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_15",1,0,14,14,15,1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_16",1,0,15,15,0,1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("start5_4" ,1,0,3,3,4,1.0,0.0,1.0,
                              Start22Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_5" ,1,0,4,4,5,1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_6" ,1,0,5,5,6,1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_7" ,1,0,6,6,7,1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("start5_1" ,1,0,0,0,1,1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5_2" ,1,0,1,1,2,1.0,0.0,1.0,
                              Start24Boundary,NULL)==NULL) return(1);
  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define a variable ring                                                  */
/*                                                                          */
/****************************************************************************/

static INT kreisBoundaryUpper (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = cos(PI*lambda);
  result[1] = sin(PI*lambda);

  return(0);
}

static INT kreisBoundaryLower (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = cos(PI+PI*lambda);
  result[1] = sin(PI+PI*lambda);

  return(0);
}

#define INNER_RADIUS2 0.8

static INT kreisBoundaryUpper1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = rad1 * cos(PI*lambda);
  result[1] = rad1 * sin(PI*lambda);

  return(0);
}

static INT kreisBoundaryLower1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = rad1 * cos(PI+PI*lambda);
  result[1] = rad1 * sin(PI+PI*lambda);

  return(0);
}

static INT kreisBoundaryUpper2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = rad1 * cos(alpha + PI*lambda);
  result[1] = rad1 * sin(alpha + PI*lambda);

  return(0);
}

static INT kreisBoundaryLower2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = rad1 * cos(alpha + PI + PI*lambda);
  result[1] = rad1 * sin(alpha + PI + PI*lambda);

  return(0);
}

static INT kreisBoundaryUpper3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;
  DOUBLE_VECTOR x,y,z;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  x[0] = rad1 * cos(alpha-0.0) * 0.8;
  x[1] = rad1 * sin(alpha-0.0) * 0.8;
  y[0] = rad1 * cos(alpha + PI-0.0) * 0.8;
  y[1] = rad1 * sin(alpha + PI-0.0) * 0.8;
  z[0] = rad1 * cos(alpha + PI * 0.5-0.0) * 0.3;
  z[1] = rad1 * sin(alpha + PI * 0.5-0.0) * 0.3;
  if (lambda <= 0.5) {
    result[0] = x[0] * (1 - lambda*2) + z[0] * lambda*2
                * (lambda - 0.4)/0.1;
    result[1] = x[1] * (1 - lambda*2) + z[1] * lambda*2
                * (lambda - 0.4)/0.1;
  }
  else {
    result[0] = z[0] * (1 - (lambda-0.5)*2)
                * (lambda - 0.4)/0.1
                + y[0] * (lambda-0.5)*2;
    result[1] = z[1] * (1 - (lambda-0.5)*2)
                * (lambda - 0.4)/0.1
                + y[1] * (lambda-0.5)*2;
  }

  return(0);
}

static INT kreisBoundaryLower3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;
  DOUBLE_VECTOR x,y,z;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  x[0] = rad1 * cos(alpha + PI-0.0) * 0.8;
  x[1] = rad1 * sin(alpha + PI-0.0) * 0.8;
  y[0] = rad1 * cos(alpha-0.0) * 0.8;
  y[1] = rad1 * sin(alpha-0.0) * 0.8;
  z[0] = rad1 * cos(alpha + PI * 1.5-0.0) * 0.3;
  z[1] = rad1 * sin(alpha + PI * 1.5-0.0) * 0.3;
  if (lambda <= 0.5) {
    result[0] = x[0] * (1 - lambda*2) + z[0] * lambda*2
                * (lambda - 0.4)/0.1;
    result[1] = x[1] * (1 - lambda*2) + z[1] * lambda*2
                * (lambda - 0.4)/0.1;
  }
  else {
    result[0] = z[0] * (1 - (lambda-0.5)*2)
                * (lambda - 0.4)/0.1
                + y[0] * (lambda-0.5)*2;
    result[1] = z[1] * (1 - (lambda-0.5)*2)
                * (lambda - 0.4)/0.1
                + y[1] * (lambda-0.5)*2;
  }

  return(0);
}

static const INT ring1_sd2p[3] = {0,0,3};
static const INT ring1_sg2p[8] = {0,0,1,1,2,2,3,3};
static const INT ring1_pt2p[8] = {0,0,1,1,2,2,3,3};
static const DOMAIN_PART_INFO ring1_dpi = {ring1_sd2p,ring1_sg2p,ring1_pt2p};

static INT InitRings1 (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = MidPoint[1] = 0.0;
  radius = 1.05;

  if (CreateDomainWithParts("Rings1",MidPoint,radius,8,8,NO,3,&ring1_dpi)
      ==NULL) return(1);

  if (CreateBoundarySegment2D("ring2 bnd upper",
                              1,0,0,0,1,20,0.0,1.0,
                              kreisBoundaryUpper,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 bnd lower",
                              1,0,1,1,0,20,0.0,1.0,
                              kreisBoundaryLower,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner bnd upper",
                              2,1,2,2,3,20,0.0,1.0,
                              kreisBoundaryUpper1,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner bnd lower",
                              2,1,3,3,2,20,0.0,1.0,
                              kreisBoundaryLower1,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner2 bnd upper",
                              2,1,4,4,5,20,0.0,1.0,
                              kreisBoundaryUpper2,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner2 bnd lower",
                              2,1,5,5,4,20,0.0,1.0,
                              kreisBoundaryLower2,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner3 bnd upper",
                              0,2,6,6,7,20,0.0,1.0,
                              kreisBoundaryUpper3,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner3 bnd lower",
                              0,2,7,7,6,20,0.0,1.0,
                              kreisBoundaryLower3,NULL)==NULL)
    return(1);

  return(0);
}

static const INT ring_sd2p[3] = {0,0,3};
static const INT ring_sg2p[8] = {0,0,1,1,2,2};
static const INT ring_pt2p[8] = {0,0,1,1,2,2};
static const DOMAIN_PART_INFO ring_dpi = {ring_sd2p,ring_sg2p,ring_pt2p};

static INT InitRings (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = MidPoint[1] = 0.0;
  radius = 1.05;

  if (CreateDomainWithParts("Rings",MidPoint,radius,6,6,YES,3,&ring_dpi)
      ==NULL) return(1);

  if (CreateBoundarySegment2D("ring2 bnd upper",
                              1,0,0,0,1,20,0.0,1.0,
                              kreisBoundaryUpper,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 bnd lower",
                              1,0,1,1,0,20,0.0,1.0,
                              kreisBoundaryLower,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner bnd upper",
                              2,1,2,2,3,20,0.0,1.0,
                              kreisBoundaryUpper1,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner bnd lower",
                              2,1,3,3,2,20,0.0,1.0,
                              kreisBoundaryLower1,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner2 bnd upper",
                              2,1,4,4,5,20,0.0,1.0,
                              kreisBoundaryUpper2,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner2 bnd lower",
                              2,1,5,5,4,20,0.0,1.0,
                              kreisBoundaryLower2,NULL)==NULL)
    return(1);

  return(0);
}

static INT InitRings2 (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = MidPoint[1] = 0.0;
  radius = 1.05;

  if (CreateDomain("Rings2",MidPoint,radius,4,4,YES)
      ==NULL) return(1);

  if (CreateBoundarySegment2D("ring2 bnd upper",
                              1,0,0,0,1,20,0.0,1.0,
                              kreisBoundaryUpper,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 bnd lower",
                              1,0,1,1,0,20,0.0,1.0,
                              kreisBoundaryLower,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner bnd upper",
                              2,1,2,2,3,20,0.0,1.0,
                              kreisBoundaryUpper1,NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("ring2 inner bnd lower",
                              2,1,3,3,2,20,0.0,1.0,
                              kreisBoundaryLower1,NULL)==NULL)
    return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define four quadrilaterals                                              */
/*                                                                          */
/****************************************************************************/

static INT south3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[6][0] + lambda*x_quad[7][0];
  result[1] = (1.0-lambda)*x_quad[6][1] + lambda*x_quad[7][1];

  return(0);
}

static INT east3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[7][0] + lambda*x_quad[1][0];
  result[1] = (1.0-lambda)*x_quad[7][1] + lambda*x_quad[1][1];

  return(0);
}

static INT west3Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[6][0] + lambda*x_quad[0][0];
  result[1] = (1.0-lambda)*x_quad[6][1] + lambda*x_quad[0][1];

  return(0);
}

static INT south4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[7][0] + lambda*x_quad[8][0];
  result[1] = (1.0-lambda)*x_quad[7][1] + lambda*x_quad[8][1];

  return(0);
}

static INT east4Boundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (1.0-lambda)*x_quad[8][0] + lambda*x_quad[4][0];
  result[1] = (1.0-lambda)*x_quad[8][1] + lambda*x_quad[4][1];

  return(0);
}

static const INT four_sd2p[ 5] = {0,0,0,0,0};
static const INT four_sg2p[16] = {1,1,0,0,2,0,0,2,0,2,2,0,0,0,1,1};
static const INT four_pt2p[16] = {1,1,1,0,3,0,0,3,0,3,3,3,3,1,1,1};
static const DOMAIN_PART_INFO four_dpi = {four_sd2p,four_sg2p,four_pt2p};

static INT InitFour (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = (1.0/9.0)*(x_quad[0][0]+x_quad[1][0]
                           +x_quad[2][0]+x_quad[3][0]
                           +x_quad[4][0]+x_quad[5][0]
                           +x_quad[6][0]+x_quad[7][0]
                           +x_quad[8][0]);
  MidPoint[1] = (1.0/9.0)*(x_quad[0][1]+x_quad[1][1]
                           +x_quad[2][1]+x_quad[3][1]
                           +x_quad[4][1]+x_quad[5][1]
                           +x_quad[6][1]+x_quad[7][1]
                           +x_quad[8][1]);
  radius =            ABS(x_quad[0][0]-MidPoint[0]);
  radius = MAX(radius,ABS(x_quad[1][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[2][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[3][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[4][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[5][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[6][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[7][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[8][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[0][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[1][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[2][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[3][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[4][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[5][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[6][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[7][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[8][1]-MidPoint[1]));

  if (CreateDomainWithParts("Four",MidPoint,radius,16,16,YES,3,&four_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("south", 1,3, 0, 0, 1,1,0.0,1.0,
                              southBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east",  1,2, 1, 1, 2,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north", 1,0, 2, 2, 3,1,0.0,1.0,
                              northBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west",  0,1, 3, 0, 3,1,0.0,1.0,
                              westBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south2",2,4, 4, 9, 4,1,0.0,1.0,
                              south2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east2", 2,0, 5, 4, 5,1,0.0,1.0,
                              east2Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north2",0,2, 6,10, 5,1,0.0,1.0,
                              north2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east1", 1,2, 7, 9,10,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south3",3,0, 8, 6, 7,1,0.0,1.0,
                              south3Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east3", 3,4, 9, 7,11,1,0.0,1.0,
                              east3Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north3",1,3,10,12,11,1,0.0,1.0,
                              southBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west3", 0,3,11, 6,12,1,0.0,1.0,
                              west3Boundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south4",4,0,12,13, 8,1,0.0,1.0,
                              south4Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east4", 4,0,13, 8,14,1,0.0,1.0,
                              east4Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north4",2,4,14,15,14,1,0.0,1.0,
                              south2Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west4", 3,4,15,13,15,1,0.0,1.0,
                              east3Boundary, NULL)==NULL) return(1);

  return(0);
}

static const INT fourp1cr_sd2p[ 5] = {0,0,3,3,0};
static const INT fourp1cr_sg2p[16] = {1,1,0,0,2,3,3,2,3,2,2,3,0,0,1,1};
static const INT fourp1cr_pt2p[16] = {1,1,1,0,2,3,3,2,0,2,2,2,2,1,1,1};
static const DOMAIN_PART_INFO fourp1cr_dpi =
{fourp1cr_sd2p,fourp1cr_sg2p,fourp1cr_pt2p};

static INT InitFour_cr (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = (1.0/9.0)*(x_quad[0][0]+x_quad[1][0]
                           +x_quad[2][0]+x_quad[3][0]
                           +x_quad[4][0]+x_quad[5][0]
                           +x_quad[6][0]+x_quad[7][0]
                           +x_quad[8][0]);
  MidPoint[1] = (1.0/9.0)*(x_quad[0][1]+x_quad[1][1]
                           +x_quad[2][1]+x_quad[3][1]
                           +x_quad[4][1]+x_quad[5][1]
                           +x_quad[6][1]+x_quad[7][1]
                           +x_quad[8][1]);
  radius =            ABS(x_quad[0][0]-MidPoint[0]);
  radius = MAX(radius,ABS(x_quad[1][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[2][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[3][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[4][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[5][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[6][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[7][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[8][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[0][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[1][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[2][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[3][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[4][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[5][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[6][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[7][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[8][1]-MidPoint[1]));

  if (CreateDomainWithParts("Four_cr",MidPoint,radius,16,16,YES,3,
                            &fourp1cr_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("south", 1,3, 0, 0, 1,1,0.0,1.0,
                              southBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east",  1,2, 1, 1, 2,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north", 1,0, 2, 2, 3,1,0.0,1.0,
                              northBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west",  0,1, 3, 0, 3,1,0.0,1.0,
                              westBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south2",2,4, 4, 9, 4,1,0.0,1.0,
                              south2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east2", 2,0, 5, 4, 5,1,0.0,1.0,
                              east2Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north2",0,2, 6,10, 5,1,0.0,1.0,
                              north2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east1", 1,2, 7, 9,10,1,0.0,1.0,
                              eastBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south3",3,0, 8, 6, 7,1,0.0,1.0,
                              south3Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east3", 3,4, 9, 7,11,1,0.0,1.0,
                              east3Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north3",1,3,10,12,11,1,0.0,1.0,
                              southBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west3", 0,3,11, 6,12,1,0.0,1.0,
                              west3Boundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south4",4,0,12,13, 8,1,0.0,1.0,
                              south4Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east4", 4,0,13, 8,14,1,0.0,1.0,
                              east4Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north4",2,4,14,15,14,1,0.0,1.0,
                              south2Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west4", 3,4,15,13,15,1,0.0,1.0,
                              east3Boundary, NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define skin-domain                                                      */
/*                                                                          */
/****************************************************************************/

static INT south0skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * lambda;
  result[1] = 2 * glob_h + D * 1.5;

  return(0);
}

static INT east0skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L;
  result[1] = 2 * glob_h + (1.5*D) + (D *0.5)* lambda;

  return(0);
}

static INT north0skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * (1.0-lambda);
  result[1] = 2 * (D + glob_h);

  return(0);
}

static INT west0skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0;
  result[1] = 2 * glob_h + (1.5 * D) + (D * 0.5) * (1 - lambda);

  return(0);
}

static INT south1skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h + (L*0.5 * lambda);
  result[1] = 2 * glob_h + 1.5 * D;

  return(0);
}

static INT east1skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 3*0.5 * L + glob_h;
  result[1] = 2 * glob_h + (1.5 * D) + (D * 0.5) * lambda;

  return(0);
}

static INT north1skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h + (L*0.5 * (1 - lambda));
  result[1] = 2 *(D + glob_h);

  return(0);
}

static INT west1skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h;
  result[1] = 2 * glob_h + (1.5 * D) + (D * 0.5) * (1 - lambda);

  return(0);
}


static INT south2skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L*0.5 * lambda;
  result[1] = D * 0.5 + glob_h;

  return(0);
}

static INT east2skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L*0.5;
  result[1] = (D * 0.5 + glob_h) + (D * lambda);

  return(0);
}

static INT north2skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L*0.5 * (1.0-lambda);
  result[1] = (D * 0.5 + glob_h) + D;

  return(0);
}

static INT west2skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0;
  result[1] = (D * 0.5 + glob_h) + D * (1.0-lambda);

  return(0);
}


static INT south3skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L*0.5 + glob_h) + L * lambda;
  result[1] = (D * 0.5 + glob_h);

  return(0);
}

static INT east3skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L*0.5 + glob_h) + L;
  result[1] = (D * 0.5 + glob_h) + D * lambda;

  return(0);
}

static INT north3skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L*0.5 + glob_h) + L * (1 - lambda);
  result[1] = (D * 0.5 + glob_h) + D;

  return(0);
}

static INT west3skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L*0.5 + glob_h);
  result[1] = (D * 0.5 + glob_h) + D * (1 - lambda);

  return(0);
}


static INT south4skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L *lambda;
  result[1] = 0.0;

  return(0);
}

static INT east4skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L ;
  result[1] = (D * 0.5) * lambda;

  return(0);
}

static INT north4skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * (1 - lambda);
  result[1] = D * 0.5;

  return(0);
}

static INT west4skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0.0;
  result[1] = (D * 0.5) * (1 - lambda);

  return(0);
}


static INT south5skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L + glob_h) + L*0.5 *lambda;
  result[1] = 0.0;

  return(0);
}

static INT east5skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L + glob_h) + L*0.5 ;
  result[1] = (D * 0.5) * lambda;

  return(0);
}

static INT north5skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L + glob_h) + L*0.5 * (1 - lambda);
  result[1] = D * 0.5;

  return(0);
}

static INT west5skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = (L + glob_h) ;
  result[1] = (D * 0.5) * (1 - lambda);

  return(0);
}

static INT north6_11skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h  + L * 0.5 * (1 - lambda);
  result[1] = 0.5 * D + glob_h ;

  return(0);
}

static INT north6_12skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h * (1 - lambda);
  result[1] = 0.5 * D + glob_h ;

  return(0);
}

static INT north6_13skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 + glob_h + (0.5 * L - glob_h) * (1 - lambda);
  result[1] = 0.5 * D + glob_h ;

  return(0);
}

static INT south6_31skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 + glob_h  + (L * 0.5 - glob_h) * lambda ;
  result[1] = 0.5 * D + glob_h + D ;

  return(0);
}

static INT south6_32skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h * lambda ;
  result[1] = 0.5 * D + glob_h + D ;

  return(0);
}

static INT south6_33skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h + 0.5 * L * lambda ;
  result[1] = 0.5 * D + glob_h + D ;

  return(0);
}

static INT north6_41skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 + glob_h  + (L * 0.5 - glob_h ) * (1 - lambda);
  result[1] = 0.5 * D + glob_h + D + glob_h;

  return(0);
}

static INT north6_42skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 + glob_h * (1 - lambda);
  result[1] = 0.5 * D + glob_h + D + glob_h ;

  return(0);
}

static INT north6_43skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 * (1 - lambda);
  result[1] = 0.5 * D + glob_h + D + glob_h ;

  return(0);
}

static INT south6_51skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 * lambda ;
  result[1] = 0.5 * D ;

  return(0);
}

static INT south6_52skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 + glob_h * lambda ;
  result[1] = 0.5 * D ;

  return(0);
}

static INT south6_53skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L * 0.5 + glob_h + (0.5 * L - glob_h ) * lambda ;
  result[1] = 0.5 * D ;

  return(0);
}

static INT north6_3skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h * (1 - lambda);
  result[1] = 2 * (D + glob_h);

  return(0);
}

static INT west6_2skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0.0;
  result[1] = (D * 0.5 + glob_h) + D + glob_h * (1.0-lambda);

  return(0);
}

static INT east6_4skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h + L*0.5 ;
  result[1] = 1.5*D + glob_h + glob_h * lambda ;

  return(0);
}

static INT west6_4skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = 0 ;
  result[1] = 0.5 * D + (1 - lambda) * glob_h ;

  return(0);
}

static INT east6_2skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + glob_h + L*0.5 ;
  result[1] = 0.5 * D + glob_h * lambda ;

  return(0);
}

static INT south6_1skBoundary (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = L + lambda * glob_h ;
  result[1] = 0 ;

  return(0);
}

static const INT skin_sd2p[8] = {0,0,0,0,0,0,0,0};
static const INT skin_sg2p[52] = {0,0,0,1,1,0,1,1,1,1,
                                  0,1,1,0,1,1,1,1,0,1,
                                  1,0,0,0,0,2,2,2,2,0,
                                  0,2,2,2,2,2,2,0,0,2,
                                  2,2,2,0,2,2,2,2,2,2,
                                  2,2};
static const INT skin_pt2p[52] = {0,1,1,0,1,1,1,1,1,1,
                                  1,1,1,1,1,1,1,1,1,1,
                                  0,1,1,0,3,3,3,3,3,3,
                                  3,3,3,3,3,3,3,3,3,3,
                                  3,3,3,3,2,2,2,2,2,2,
                                  2,2};

static const DOMAIN_PART_INFO skin_dpi = {skin_sd2p,skin_sg2p,skin_pt2p};

static INT InitSkin (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = L*0.75 + glob_h*0.5 ;
  MidPoint[1] = glob_h + D ;
  radius = sqrt(MidPoint[0] * MidPoint[0] + MidPoint[1] * MidPoint[1]) ;

  if (CreateDomainWithParts("Skin",MidPoint,radius,52,52,YES,3,&skin_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("south0", 1, 7, 6, 4, 5,1,0.0,1.0,
                              south0skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east0",  1, 7, 3, 5, 1,1,0.0,1.0,
                              east0skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north0", 1, 0, 0, 1, 0,1,0.0,1.0,
                              north0skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west0",  1, 0, 2, 0, 4,1,0.0,1.0,
                              west0skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south1", 2, 7, 7, 6, 7,1,0.0,1.0,
                              south1skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east1",  2, 0, 5, 7, 3,1,0.0,1.0,
                              east1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north1", 2, 0, 1, 3, 2,1,0.0,1.0,
                              north1skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west1",  2, 7, 4, 2, 6,1,0.0,1.0,
                              west1skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south2", 3, 7,14,12,13,1,0.0,1.0,
                              south2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east2",  3, 7,11,13, 9,1,0.0,1.0,
                              east2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north2", 3, 7, 8, 9, 8,1,0.0,1.0,
                              north2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west2",  3, 0,10, 8,12,1,0.0,1.0,
                              west2skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south3", 4, 7,15,14,15,1,0.0,1.0,
                              south3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east3",  4, 0,13,15,11,1,0.0,1.0,
                              east3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north3", 4, 7, 9,11,10,1,0.0,1.0,
                              north3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west3",  4, 7,12,10,14,1,0.0,1.0,
                              west3skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south4", 5, 0,22,20,21,1,0.0,1.0,
                              south4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east4",  5, 7,19,21,17,1,0.0,1.0,
                              east4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north4", 5, 7,16,17,16,1,0.0,1.0,
                              north4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west4",  5, 0,18,16,20,1,0.0,1.0,
                              west4skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south5", 6, 0,23,22,23,1,0.0,1.0,
                              south5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east5",  6, 0,21,23,19,1,0.0,1.0,
                              east5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north5", 6, 7,17,19,18,1,0.0,1.0,
                              north5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west5",  6, 7,20,18,22,1,0.0,1.0,
                              west5skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south6_1", 7, 0,43,42,43,1,0.0,1.0,
                              south6_1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_1" , 6, 7,42,40,43,1,0.0,1.0,
                              west5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_2", 6, 7,40,41,40,1,0.0,1.0,
                              north5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_2" , 7, 0,38,41,37,1,0.0,1.0,
                              east6_2skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_11", 7, 4,49,37,49,1,0.0,1.0,
                              north6_11skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_12", 7, 4,48,49,48,1,0.0,1.0,
                              north6_12skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_13", 7, 4,36,48,36,1,0.0,1.0,
                              north6_13skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_3" , 4, 7,34,32,36,1,0.0,1.0,
                              west3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_31", 7, 4,32,32,46,1,0.0,1.0,
                              south6_31skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_32", 7, 4,46,46,47,1,0.0,1.0,
                              south6_32skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_33", 7, 4,47,47,33,1,0.0,1.0,
                              south6_33skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_4" , 7, 0,30,33,29,1,0.0,1.0,
                              east6_4skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_2", 2, 7,28,28,29,1,0.0,1.0,
                              south1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_5" , 2, 7,26,25,28,1,0.0,1.0,
                              west1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_3", 7, 0,24,25,24,1,0.0,1.0,
                              north6_3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_1" , 1, 7,25,27,24,1,0.0,1.0,
                              east0skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_41", 7, 1,45,27,45,1,0.0,1.0,
                              north6_41skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_42", 7, 1,44,45,44,1,0.0,1.0,
                              north6_42skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_43", 7, 1,27,44,26,1,0.0,1.0,
                              north6_43skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_2" , 7, 0,29,26,30,1,0.0,1.0,
                              west6_2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_4", 3, 7,31,31,30,1,0.0,1.0,
                              north2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_3" , 3, 7,33,35,31,1,0.0,1.0,
                              east2skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_5", 3, 7,35,34,35,1,0.0,1.0,
                              south2skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_4" , 7, 0,37,34,38,1,0.0,1.0,
                              west6_4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_51", 7, 5,39,38,50,1,0.0,1.0,
                              south6_51skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_52", 7, 5,50,50,51,1,0.0,1.0,
                              south6_52skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_53", 7, 5,51,51,39,1,0.0,1.0,
                              south6_53skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_5" , 5, 7,41,42,39,1,0.0,1.0,
                              east4skBoundary, NULL)==NULL) return(1);

  return(0);
}

static const INT skin1_sd2p[8] = {0,0,0,0,0,0,0,0};
static const INT skin1_sg2p[52] = {0,0,0,2,2,0,2,2,2,2,
                                   0,2,2,0,2,2,2,2,0,2,
                                   2,0,0,0,0,1,1,1,1,0,
                                   0,1,1,1,1,1,1,0,0,1,
                                   1,1,1,0,1,1,1,1,1,1,
                                   1,1};
static const INT skin1_pt2p[52] = {0,3,3,0,3,3,3,3,3,3,
                                   3,3,3,3,3,3,3,3,3,3,
                                   0,3,3,0,
                                   1,1,1,1,1,1,
                                   1,1,1,1,1,1,1,1,1,1,
                                   1,1,1,1,1,1,1,1,1,1,
                                   1,1};

static const DOMAIN_PART_INFO skin1_dpi = {skin1_sd2p,skin1_sg2p,skin1_pt2p};

static INT InitSkin1 (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = L*0.75 + glob_h*0.5 ;
  MidPoint[1] = glob_h + D ;
  radius = sqrt(MidPoint[0] * MidPoint[0] + MidPoint[1] * MidPoint[1]) ;

  if (CreateDomainWithParts("Skin1",MidPoint,radius,52,52,YES,3,&skin1_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("south0", 1, 7, 6, 4, 5,1,0.0,1.0,
                              south0skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east0",  1, 7, 3, 5, 1,1,0.0,1.0,
                              east0skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north0", 1, 0, 0, 1, 0,1,0.0,1.0,
                              north0skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west0",  1, 0, 2, 0, 4,1,0.0,1.0,
                              west0skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south1", 2, 7, 7, 6, 7,1,0.0,1.0,
                              south1skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east1",  2, 0, 5, 7, 3,1,0.0,1.0,
                              east1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north1", 2, 0, 1, 3, 2,1,0.0,1.0,
                              north1skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west1",  2, 7, 4, 2, 6,1,0.0,1.0,
                              west1skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south2", 3, 7,14,12,13,1,0.0,1.0,
                              south2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east2",  3, 7,11,13, 9,1,0.0,1.0,
                              east2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north2", 3, 7, 8, 9, 8,1,0.0,1.0,
                              north2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west2",  3, 0,10, 8,12,1,0.0,1.0,
                              west2skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south3", 4, 7,15,14,15,1,0.0,1.0,
                              south3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east3",  4, 0,13,15,11,1,0.0,1.0,
                              east3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north3", 4, 7, 9,11,10,1,0.0,1.0,
                              north3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west3",  4, 7,12,10,14,1,0.0,1.0,
                              west3skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south4", 5, 0,22,20,21,1,0.0,1.0,
                              south4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east4",  5, 7,19,21,17,1,0.0,1.0,
                              east4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north4", 5, 7,16,17,16,1,0.0,1.0,
                              north4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west4",  5, 0,18,16,20,1,0.0,1.0,
                              west4skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south5", 6, 0,23,22,23,1,0.0,1.0,
                              south5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east5",  6, 0,21,23,19,1,0.0,1.0,
                              east5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north5", 6, 7,17,19,18,1,0.0,1.0,
                              north5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west5",  6, 7,20,18,22,1,0.0,1.0,
                              west5skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("south6_1", 7, 0,43,42,43,1,0.0,1.0,
                              south6_1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_1" , 6, 7,42,40,43,1,0.0,1.0,
                              west5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_2", 6, 7,40,41,40,1,0.0,1.0,
                              north5skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_2" , 7, 0,38,41,37,1,0.0,1.0,
                              east6_2skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_11", 7, 4,49,37,49,1,0.0,1.0,
                              north6_11skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_12", 7, 4,48,49,48,1,0.0,1.0,
                              north6_12skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_13", 7, 4,36,48,36,1,0.0,1.0,
                              north6_13skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_3" , 4, 7,34,32,36,1,0.0,1.0,
                              west3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_31", 7, 4,32,32,46,1,0.0,1.0,
                              south6_31skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_32", 7, 4,46,46,47,1,0.0,1.0,
                              south6_32skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_33", 7, 4,47,47,33,1,0.0,1.0,
                              south6_33skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_4" , 7, 0,30,33,29,1,0.0,1.0,
                              east6_4skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_2", 2, 7,28,28,29,1,0.0,1.0,
                              south1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("east6_5" , 2, 7,26,25,28,1,0.0,1.0,
                              west1skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_3", 7, 0,24,25,24,1,0.0,1.0,
                              north6_3skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_1" , 1, 7,25,27,24,1,0.0,1.0,
                              east0skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_41", 7, 1,45,27,45,1,0.0,1.0,
                              north6_41skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_42", 7, 1,44,45,44,1,0.0,1.0,
                              north6_42skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("north6_43", 7, 1,27,44,26,1,0.0,1.0,
                              north6_43skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_2" , 7, 0,29,26,30,1,0.0,1.0,
                              west6_2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_4", 3, 7,31,31,30,1,0.0,1.0,
                              north2skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_3" , 3, 7,33,35,31,1,0.0,1.0,
                              east2skBoundary, NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("north6_5", 3, 7,35,34,35,1,0.0,1.0,
                              south2skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_4" , 7, 0,37,34,38,1,0.0,1.0,
                              west6_4skBoundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_51", 7, 5,39,38,50,1,0.0,1.0,
                              south6_51skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_52", 7, 5,50,50,51,1,0.0,1.0,
                              south6_52skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("south6_53", 7, 5,51,51,39,1,0.0,1.0,
                              south6_53skBoundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("west6_5" , 5, 7,41,42,39,1,0.0,1.0,
                              east4skBoundary, NULL)==NULL) return(1);

  return(0);
}


/* end of skin-domain-definition */


/****************************************************************************/
/*                                                                          */
/*  define Composed1 (A Part of the domain "Holes")                         */
/*    -----------------                                                     */
/*    |               |                                                     */
/*    -----------------                                                     */
/*    |  |   |  |  |  |                                                     */
/*    ----   ----  ----                                                     */
/*                                                                          */
/****************************************************************************/

static const INT c1_sd2p[2] = {0,0};
static const INT c1_sg2p[24] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const INT c1_pt2p[24] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const DOMAIN_PART_INFO c1_dpi = {c1_sd2p,c1_sg2p,c1_pt2p};

static INT InitComposed1 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomainWithParts("Composed1",MidPoint,radius,24,24,NO,2,&c1_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("left1"  ,1,0, 0, 0,1, 1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left2"  ,1,0, 1, 1,2, 1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left3"  ,1,0, 2, 2,3, 1.0,0.0,1.0,
                              Start7Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left4"  ,1,0, 3, 3,0, 1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle1",1,0, 4, 4,5, 1.0,0.0,1.0,
                              Start25Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle2",1,0, 5, 5,6, 1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle3",1,0, 6, 6,7, 1.0,0.0,1.0,
                              Start27Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle4",1,0, 7, 7,4, 1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right1" ,1,0, 8, 8,9, 1.0,0.0,1.0,
                              Start45Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right2" ,1,0, 9, 9,10, 1.0,0.0,1.0,
                              Start46Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right3" ,1,0, 10, 10,11, 1.0,0.0,1.0,
                              Start47Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right4" ,1,0, 11, 11,8, 1.0,0.0,1.0,
                              Start48Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth1",1,0, 12, 12,13, 1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth2",1,0, 13, 13,14, 1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth3",1,0, 14, 14,15, 1.0,0.0,1.0,
                              Start29Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth4",1,0, 15, 15,16, 1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth5",1,0, 16, 16,17, 1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppereast",1,0, 17, 17,18, 1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth1",1,0, 18, 18,19, 1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppernorth2",1,0, 19, 19,20, 1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth3",1,0, 20, 20,21, 1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth4",1,0, 21, 21,22, 1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth5",1,0, 22, 22,23, 1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("upperwest",1,0, 23, 23,12, 1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);

  return(0);
}

static INT InitComposed1a (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Composed1a",MidPoint,radius,18,18,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("left1"  ,1,0, 0, 0,1, 1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left2"  ,1,0, 1, 1,2, 1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth2",1,0, 2, 2,3, 1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("middle4",1,0, 3, 3,4, 1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle1",1,0, 4, 4,5, 1.0,0.0,1.0,
                              Start25Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle2",1,0, 5, 5,6, 1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth4",1,0, 6, 6,7, 1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("right4" ,1,0, 7, 7,8, 1.0,0.0,1.0,
                              Start48Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right1" ,1,0, 8, 8,9, 1.0,0.0,1.0,
                              Start45Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("right2" ,1,0, 9, 9,10, 1.0,0.0,1.0,
                              Start46Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppereast",1,0, 10, 10,11, 1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth1",1,0, 11, 11,12, 1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth2",1,0, 12, 12,13, 1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth3",1,0, 13, 13,14, 1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth4",1,0, 14, 14,15, 1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth5",1,0, 15, 15,16, 1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("upperwest",1,0, 16, 16,17, 1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left4"  ,1,0, 17, 17,0, 1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define Composed2 (A Part of the domain "Holes")                         */
/*    -----------------                                                     */
/*    |               |                                                     */
/*    -----------------                                                     */
/*    |  |   |  |                                                           */
/*    ----   ----                                                           */
/*                                                                          */
/****************************************************************************/

static const INT c2_sd2p[2] = {0,0};
static const INT c2_sg2p[20] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const INT c2_pt2p[20] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const DOMAIN_PART_INFO c2_dpi = {c2_sd2p,c2_sg2p,c2_pt2p};

static INT InitComposed2 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomainWithParts("Composed2",MidPoint,radius,20,20,NO,2,&c2_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("left1"  ,1,0, 0, 0,1, 1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left2"  ,1,0, 1, 1,2, 1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left3"  ,1,0, 2, 2,3, 1.0,0.0,1.0,
                              Start7Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left4"  ,1,0, 3, 3,0, 1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle1",1,0, 4, 4,5, 1.0,0.0,1.0,
                              Start25Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle2",1,0, 5, 5,6, 1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle3",1,0, 6, 6,7, 1.0,0.0,1.0,
                              Start27Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle4",1,0, 7, 7,4, 1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth1",1,0, 8, 8,9, 1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth2",1,0, 9, 9,10, 1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth3",1,0, 10, 10,11, 1.0,0.0,1.0,
                              Start29Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth4",1,0, 11, 11,12, 1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth5",1,0, 12, 12,13, 1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppereast",1,0, 13, 13,14, 1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth1",1,0, 14, 14,15, 1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppernorth2",1,0, 15, 15,16, 1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth3",1,0, 16, 16,17, 1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth4",1,0, 17, 17,18, 1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth5",1,0, 18, 18,19, 1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("upperwest",1,0, 19, 19,8, 1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);

  return(0);
}

static INT InitComposed2a (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Composed2a",MidPoint,radius,16,16,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("left1"  ,1,0, 0, 0,1, 1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left2"  ,1,0, 1, 1,2, 1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth2",1,0, 2, 2,3, 1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("middle4",1,0, 3, 3,4, 1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle1",1,0, 4, 4,5, 1.0,0.0,1.0,
                              Start25Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("middle2",1,0, 5, 5,6, 1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth4",1,0, 6, 6,7, 1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth5",1,0, 7, 7,8, 1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppereast",1,0, 8, 8,9, 1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth1",1,0, 9, 9,10, 1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth2",1,0, 10, 10,11, 1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth3",1,0, 11, 11,12, 1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth4",1,0, 12, 12,13, 1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth5",1,0, 13, 13,14, 1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("upperwest",1,0, 14, 14,15, 1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left4"  ,1,0, 15, 15,0, 1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define Composed3 (A Part of the domain "Holes")                         */
/*    -----------------                                                     */
/*    |               |                                                     */
/*    -----------------                                                     */
/*    |  |                                                                  */
/*    ----                                                                  */
/*                                                                          */
/****************************************************************************/

static const INT c3_sd2p[2] = {0,0};
static const INT c3_sg2p[16] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const INT c3_pt2p[16] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static const DOMAIN_PART_INFO c3_dpi = {c3_sd2p,c3_sg2p,c3_pt2p};

static INT InitComposed3 (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomainWithParts("Composed3",MidPoint,radius,16,16,NO,2,&c3_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("left1"  ,1,0, 0, 0,1, 1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left2"  ,1,0, 1, 1,2, 1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left3"  ,1,0, 2, 2,3, 1.0,0.0,1.0,
                              Start7Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left4"  ,1,0, 3, 3,0, 1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth1",1,0, 4, 4,5, 1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth2",1,0, 5, 5,6, 1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth3",1,0, 6, 6,7, 1.0,0.0,1.0,
                              Start29Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth4",1,0, 7, 7,8, 1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth5",1,0, 8, 8,9, 1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppereast",1,0, 9, 9,10, 1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth1",1,0, 10, 10,11, 1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppernorth2",1,0, 11, 11,12, 1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth3",1,0, 12, 12,13, 1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth4",1,0, 13, 13,14, 1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth5",1,0, 14, 14,15, 1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("upperwest",1,0, 15, 15,4, 1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);

  return(0);
}

static INT InitComposed3a (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomain("Composed3a",MidPoint,radius,14,14,NO)==NULL)
    return(1);

  if (CreateBoundarySegment2D("left1"  ,1,0, 0, 0,1, 1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left2"  ,1,0, 1, 1,2, 1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppersouth2",1,0, 2, 2,3, 1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth3",1,0, 3, 3,4, 1.0,0.0,1.0,
                              Start29Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth4",1,0, 4, 4,5, 1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppersouth5",1,0, 5, 5,6, 1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);

  if (CreateBoundarySegment2D("uppereast",1,0, 6, 6,7, 1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth1",1,0, 7, 7,8, 1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth2",1,0, 8, 8,9, 1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth3",1,0, 9, 9,10, 1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth4",1,0, 10, 10,11, 1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("uppernorth5",1,0, 11, 11,12, 1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("upperwest",1,0, 12, 12,13, 1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("left4"  ,1,0, 13, 13,0, 1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  define Beam-domain                                                    */
/*                                                                          */
/****************************************************************************/


static INT T_Beam2_0 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = x_quad[0][0]+(x_quad[1][0]-x_quad[0][0])*lambda;
  result[1] = x_quad[0][1]+(x_quad[1][1]-x_quad[0][1])*lambda;

  return(0);
}

static INT T_Beam2_1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = x_quad[1][0]+(x_quad[2][0]-x_quad[1][0])*lambda;
  result[1] = x_quad[1][1]+(x_quad[2][1]-x_quad[1][1])*lambda;

  return(0);
}
static INT T_Beam2_2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda,a,alpha;

  lambda = param[0];
  a = sqrt((x_quad[3][1]-x_quad[2][1])*(x_quad[3][1]-x_quad[2][1])
           +(x_quad[3][0]-x_quad[2][0])*(x_quad[3][0]-x_quad[2][0]));
  alpha = atan((x_quad[3][1]-x_quad[2][1])/(x_quad[3][0]-x_quad[2][0]));

  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] =  x_quad[2][0]+(x_quad[3][0]-x_quad[2][0])*lambda
              - cos(alpha)*(-a*lambda + a*lambda*lambda)*form;
  result[1] =  x_quad[2][1]+(x_quad[3][1]-x_quad[2][1])*lambda
              + sin(alpha)*(-a*lambda + a*lambda*lambda)*form;

  return(0);
}

static INT T_Beam2_3 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = x_quad[3][0]+(x_quad[4][0]-x_quad[3][0])*lambda;
  result[1] = x_quad[3][1]+(x_quad[4][1]-x_quad[3][1])*lambda;

  return(0);
}

static INT T_Beam2_4 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = x_quad[4][0]+(x_quad[5][0]-x_quad[4][0])*lambda;
  result[1] = x_quad[4][1]+(x_quad[5][1]-x_quad[4][1])*lambda;

  return(0);
}

static INT T_Beam2_5 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = x_quad[5][0]+(x_quad[6][0]-x_quad[5][0])*lambda;
  result[1] = x_quad[5][1]+(x_quad[6][1]-x_quad[5][1])*lambda;

  return(0);
}

static INT T_Beam2_6 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = x_quad[6][0]+(x_quad[7][0]-x_quad[6][0])*lambda;
  result[1] = x_quad[6][1]+(x_quad[7][1]-x_quad[6][1])*lambda;

  return(0);
}

static INT T_Beam2_7 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];
  if ((lambda<0.0)||(lambda>1.0)) return(1);
  result[0] = x_quad[7][0]+(x_quad[0][0]-x_quad[7][0])*lambda;
  result[1] = x_quad[7][1]+(x_quad[0][1]-x_quad[7][1])*lambda;

  return(0);
}

static INT InitBeam (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = 0.125*(x_quad[0][0]+x_quad[1][0]+x_quad[2][0]+
                       x_quad[3][0]+x_quad[4][0]+x_quad[5][0]+
                       x_quad[6][0]+x_quad[7][0]) ;
  MidPoint[1] = 0.125*(x_quad[0][1]+x_quad[1][1]+x_quad[2][1]+
                       x_quad[3][1]+x_quad[4][1]+x_quad[5][1]+
                       x_quad[6][1]+x_quad[7][1]) ;

  radius =            ABS(x_quad[0][0]-MidPoint[0]);
  radius = MAX(radius,ABS(x_quad[1][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[2][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[3][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[4][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[5][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[6][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[7][0]-MidPoint[0]));
  radius = MAX(radius,ABS(x_quad[0][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[1][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[2][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[3][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[4][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[5][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[6][1]-MidPoint[1]));
  radius = MAX(radius,ABS(x_quad[7][1]-MidPoint[1]));

  if (CreateDomain("Beam",MidPoint,radius,8,8,YES)==NULL)
    return(1);

  if (CreateBoundarySegment2D("T_Beam_2_0",1,0,0,0,1, 1,0.0,1.0,
                              T_Beam2_0,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("T_Beam_2_1",1,0,1,1,2, 1,0.0,1.0,
                              T_Beam2_1,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("T_Beam_2_2",1,0,2,2,3,20,0.0,1.0,
                              T_Beam2_2,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("T_Beam_2_3",1,0,3,3,4, 1,0.0,1.0,
                              T_Beam2_3,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("T_Beam_2_4",1,0,4,4,5, 1,0.0,1.0,
                              T_Beam2_4,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("T_Beam_2_5",1,0,5,5,6, 1,0.0,1.0,
                              T_Beam2_5,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("T_Beam_2_6",1,0,6,6,7, 1,0.0,1.0,
                              T_Beam2_6,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("T_Beam_2_7",1,0,7,7,0, 1,0.0,1.0,
                              T_Beam2_7,NULL)==NULL) return(1);

  return(0);
}
/* end of T_Beam-domain-definition */



/****************************************************************************/
/*                                                                          */
/*  Definition des Brenner/Scott Titelgebietes, Name : Channel
 */
/*                                                                          */
/****************************************************************************/

/* Zur vollstaendigen Beschreibung eines Randsteueckes benoetigen wir
 * lediglich die Koordinaten der Endpunkte, weil alle Randstuecke Teile
   von
 * Geraden sind. Die Endpunkte werden von (0,0) linksherum fortlaufend
 * durchnumeriert. Es gibt 26 solcher Endpunkte mit den folgenden
 * Koordinaten:
 */

#define NAME(n) Segment ## n
#define QUOTE(t) STR(t)
#define NO_OF_SEG 26

static DOUBLE SegCorner[NO_OF_SEG][2] = { {.00, .00},  {.19, .00},
                                          {.19, .18},  {.13, .29},
                                          {.13, .48},  {.19, .59},
                                          {.19, .82},  {.25, .82},
                                          {.25, .00},  {.69, .00},
                                          {.69, .53},  {.88, .59},
                                          {.88, .00},  {1.0, .00},
                                          {1.0, 1.0},  {.81, .65},
                                          {.75, .88},  {.56, .88},
                                          {.56, .12},  {.31, .12},
                                          {.31, .30},  {.50, .47},
                                          {.50, .71},  {.33, .82},
                                          {.33, 1.0},  {.00, 1.0} };



#define LINE_FROM_POINT_TO_POINT(number, startpt, endpt) \
  static INT Segment ## number  (void *data, DOUBLE *param, DOUBLE *result)\
  {\
    DOUBLE lambda;\
    ASSERT(startpt < NO_OF_SEG && endpt < NO_OF_SEG);\
    lambda = param[0];\
    if ((lambda<0.0)||(lambda>1.0)) return (1);\
    result[0] = (1.0 - lambda) * SegCorner[startpt][0] \
                +lambda * SegCorner[endpt]  [0];\
    result[1] = (1.0 - lambda) * SegCorner[startpt][1] \
                +lambda * SegCorner[endpt][1];\
    return (0); }

LINE_FROM_POINT_TO_POINT(0,0,1)
LINE_FROM_POINT_TO_POINT(1,1,2)
LINE_FROM_POINT_TO_POINT(2,2,3)
LINE_FROM_POINT_TO_POINT(3,3,4)
LINE_FROM_POINT_TO_POINT(4,4,5)
LINE_FROM_POINT_TO_POINT(5,5,6)
LINE_FROM_POINT_TO_POINT(6,6,7)
LINE_FROM_POINT_TO_POINT(7,7,8)
LINE_FROM_POINT_TO_POINT(8,8,9)
LINE_FROM_POINT_TO_POINT(9,9,10)
LINE_FROM_POINT_TO_POINT(10,10,11)
LINE_FROM_POINT_TO_POINT(11,11,12)
LINE_FROM_POINT_TO_POINT(12,12,13)
LINE_FROM_POINT_TO_POINT(13,13,14)
LINE_FROM_POINT_TO_POINT(14,14,15)
LINE_FROM_POINT_TO_POINT(15,15,16)
LINE_FROM_POINT_TO_POINT(16,16,17)
LINE_FROM_POINT_TO_POINT(17,17,18)
LINE_FROM_POINT_TO_POINT(18,18,19)
LINE_FROM_POINT_TO_POINT(19,19,20)
LINE_FROM_POINT_TO_POINT(20,20,21)
LINE_FROM_POINT_TO_POINT(21,21,22)
LINE_FROM_POINT_TO_POINT(22,22,23)
LINE_FROM_POINT_TO_POINT(23,23,24)
LINE_FROM_POINT_TO_POINT(24,24,25)
LINE_FROM_POINT_TO_POINT(25,25,0 )

/*
 *  so, die naechsten Punkte (das sind 26-33) und die dazugehhoerigen
 *  Segmente (26-33) definieren den Part zwei, das Gebiet zwischen den
   beiden
 *  linken Saeulen
 *
 */

LINE_FROM_POINT_TO_POINT(26,1,8)

/*
 *  so, die naechsten Punkte (das sind 34-37 und die dazugehhoerigen
 *  Segmente (26-33) definieren den Part drei, das Gebiet zwischen den
   beiden
 *  rechten Saeulen
 *
 */

LINE_FROM_POINT_TO_POINT(34,9,12)


/*
 *  so, die naechsten Punkte (das sind 38-48 und die dazugehhoerigen
 *  Segmente (26-33) definieren den Part vier, das Gebiet rechts oben
 *
 */

LINE_FROM_POINT_TO_POINT(38,14,24)
/*
    Numerierung: 0 - 25: Schlauch
                26 - 33: links unten
                34 - 37: rechts unten
                38 - 48: rechts oben

    Der Schlauch wurde von links unten (Koordinate (0,0)) bis in die
   rechte obere
    Ecke (1,1), Nummer 14, fortlaufend durchnumeriert, dann weiter
   linksherum bis
    Nummer 25 oben links (0,1)
    Weiter sind folgende Punkte am selben Ort

    1-26, 8-27,7,28,6-29,etc. bis 33-2
    9-34,12-35,11-36,10-37
    14-38,24-39,23-40,22-41, etc. bis 15-48
 */

static const INT channelp1cr_sd2p[5] = {0,3,0,0,0};
static const INT channelp1cr_sg2p[49] = {3,2,2,2,2,2,2,2,3,2, /* 10
                                                                 Zahlen pro Zeile */
                                         2,2,3,3,2,2,2,2,2,2,
                                         2,2,2,2,3,3,         /* das war
                                                                 der Schlauch */
                                         0,1,1,1,1,1,1,1,     /* links
                                                                 unten */
                                         0,1,1,1,             /* rechts
                                                                 unten */
                                         0,1,1,1,1,1,1,1,1,1, /* rechts
                                                                 oben */
                                         1};                  /* auch
                                                                 rechts oben */
static const INT channelp1cr_pt2p[49] = {3,2,2,2,2,2,2,2,2,2, /* 10
                                                                 Zahlen pro Zeile */
                                         2,2,2,3,2,2,2,2,2,2,
                                         2,2,2,2,2,3,         /* das war
                                                                 der Schlauch */
                                         1,1,1,1,1,1,1,1,     /* links
                                                                 unten */
                                         1,1,1,1,             /* rechts
                                                                 unten */
                                         1,1,1,1,1,1,1,1,1,1, /* rechts
                                                                 oben */
                                         1};
static const DOMAIN_PART_INFO channelp1cr_dpi =
{channelp1cr_sd2p,channelp1cr_sg2p,channelp1cr_pt2p};


static INT InitChannel (void)
{
  DOUBLE radius,MidPoint[2];

  MidPoint[0] = MidPoint[1] = 0.5;
  radius = 1;

  if (CreateDomainWithParts("Channel",MidPoint,radius,49,49,NO,3,
                            &channelp1cr_dpi)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(0)),1,0, 0,0,1,
                              1,0.0,1.0,NAME(0),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(1)),1,0, 1,1,2,
                              1,0.0,1.0,NAME(1),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(2)),1,0, 2,2,3,
                              1,0.0,1.0,NAME(2),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(3)),1,0, 3,3,4,
                              1,0.0,1.0,NAME(3),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(4)),1,0, 4,4,5,
                              1,0.0,1.0,NAME(4),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(5)),1,0, 5,5,6,
                              1,0.0,1.0,NAME(5),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(6)),1,0, 6,6,7,
                              1,0.0,1.0,NAME(6),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(7)),1,0, 7,7,8,
                              1,0.0,1.0,NAME(7),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(8)),1,0, 8,8,9,
                              1,0.0,1.0,NAME(8),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(9)),1,0, 9,9,10,
                              1,0.0,1.0,NAME(9),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(10)),1,0, 10,10,11,
                              1,0.0,1.0,NAME(10),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(11)),1,0, 11,11,12,
                              1,0.0,1.0,NAME(11),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(12)),1,0, 12,12,13,
                              1,0.0,1.0,NAME(12),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(13)),1,0, 13,13,14,
                              1,0.0,1.0,NAME(13),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(14)),1,0, 14,14,15,
                              1,0.0,1.0,NAME(14),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(15)),1,0, 15,15,16,
                              1,0.0,1.0,NAME(15),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(16)),1,0, 16,16,17,
                              1,0.0,1.0,NAME(16),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(17)),1,0, 17,17,18,
                              1,0.0,1.0,NAME(17),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(18)),1,0, 18,18,19,
                              1,0.0,1.0,NAME(18),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(19)),1,0, 19,19,20,
                              1,0.0,1.0,NAME(19),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(20)),1,0, 20,20,21,
                              1,0.0,1.0,NAME(20),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(21)),1,0, 21,21,22,
                              1,0.0,1.0,NAME(21),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(22)),1,0, 22,22,23,
                              1,0.0,1.0,NAME(22),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(23)),1,0, 23,23,24,
                              1,0.0,1.0,NAME(23),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(24)),1,0, 24,24,25,
                              1,0.0,1.0,NAME(24),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(25)),1,0, 25,25,0,
                              1,0.0,1.0,NAME(25),NULL)==NULL)
    return(1);

  /* Subdomain 2 */
  if (CreateBoundarySegment2D("Segment26", 2,0, 26,26,27,
                              1,0.0,1.0,NAME(26),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment27",0,2, 27,28,27,
                              1,0.0,1.0,NAME(7),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment28",0,2, 28,29,28,
                              1,0.0,1.0,NAME(6),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment29",0,2, 29,30,29,
                              1,0.0,1.0,NAME(5), NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment30",0,2, 30,31,30,
                              1,0.0,1.0,NAME(4),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment31",0,2, 31,32,31,
                              1,0.0,1.0,NAME(3),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment32",0,2, 32,33,32,
                              1,0.0,1.0,NAME(2),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment33",0,2, 33,26,33,
                              1,0.0,1.0,NAME(1),NULL)==NULL)
    return(1);

  /* Subdomain 3 */
  if (CreateBoundarySegment2D(QUOTE(NAME(34)),3,0, 34,34,35,
                              1,0.0,1.0,NAME(34),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment35",0,3, 35,36,35,
                              1,0.0,1.0,NAME(11),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment36",0,3, 36,37,36,
                              1,0.0,1.0,NAME(10),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment37",0,3, 37,34,37,
                              1,0.0,1.0,NAME(9),NULL)==NULL)
    return(1);

  /* Subdomain 4 */
  if (CreateBoundarySegment2D("Segment38",4,0, 38,38,39,
                              1,0.0,1.0,NAME(38),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment39",0,4, 39,40,39,
                              1,0.0,1.0,NAME(23),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment40",0,4, 40,41,40,
                              1,0.0,1.0,NAME(22),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment41",0,4, 41,42,41,
                              1,0.0,1.0,NAME(21),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment42",0,4, 42,43,42,
                              1,0.0,1.0,NAME(20),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(43)),0,4, 43,44,43,
                              1,0.0,1.0,NAME(19),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment44",0,4, 44,45,44,
                              1,0.0,1.0,NAME(18),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment45",0,4, 45,46,45,
                              1,0.0,1.0,NAME(17),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment46",0,4, 46,47,46,
                              1,0.0,1.0,NAME(16),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment47",0,4, 47,48,47,
                              1,0.0,1.0,NAME(15),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D("Segment48",0,4, 48,38,48,
                              1,0.0,1.0,NAME(14),NULL)==NULL)
    return(1);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/*  Definition des Brenner/Scott Titelgebietes, Name : Channel
 */
/*                                                                          */
/****************************************************************************/

/* Zur vollstaendigen Beschreibung eines Randsteueckes benoetigen wir
 * lediglich die Koordinaten der Endpunkte, weil alle Randstuecke Teile
   von
 * Geraden sind. Die Endpunkte werden von (0,0) linksherum fortlaufend
 * durchnumeriert. Es gibt 26 solcher Endpunkte mit den folgenden
 * Koordinaten:
 */


static INT InitChannelNoParts (void)
{
  DOUBLE radius,MidPoint[2];

  /* allocate new domain structure */
  MidPoint[0] = MidPoint[1] = 0.5;
  radius = 1.05;

  if (CreateDomain(
        "ChannelNoParts",               /* name of the new domain                               */
        MidPoint,radius,                /* circle containing the domain			*/
        NO_OF_SEG+3,                            /* number of boundary segments                  */
        NO_OF_SEG,                              /* number of corners					*/
        NO                                                      /* true if domain is convex				*/
        )==NULL) return(1);
  if (CreateBoundarySegment2D(
        QUOTE(NAME(0)),                         /* name of the boundary segment                 */
        1,                                                      /* number of left subdomain				*/
        0,                                                      /* number of right subdomain			*/
        0,                                                      /* number of segment, starting with 0	*/
        0,                                                      /* number of corner where segment starts*/
        1,                                                      /* number of corner where segment ends  */
        1,                                              /* resolution, use 1 for straight line  */
        0.0,                                            /* begin of parameter interval			*/
        1.0,                                            /* end of parameter interval			*/
        NAME(0),                                /* function mapping parameter to world  */
        NULL                                            /* user defined pointer to be supplied  */
        )==NULL) return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(1)),1,2, 1,1,2,
                              1,0.0,1.0,NAME(1),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(2)),1,2, 2,2,3,
                              1,0.0,1.0,NAME(2),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(3)),1,2, 3,3,4,
                              1,0.0,1.0,NAME(3),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(4)),1,2, 4,4,5,
                              1,0.0,1.0,NAME(4),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(5)),1,2, 5,5,6,
                              1,0.0,1.0,NAME(5),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(6)),1,2, 6,6,7,
                              1,0.0,1.0,NAME(6),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(7)),1,2, 7,7,8,
                              1,0.0,1.0,NAME(7),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(8)),1,0, 8,8,9,
                              1,0.0,1.0,NAME(8),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(9)),1,3, 9,9,10,
                              1,0.0,1.0,NAME(9),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(10)),1,3, 10,10,11,
                              1,0.0,1.0,NAME(10),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(11)),1,3, 11,11,12,
                              1,0.0,1.0,NAME(11),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(12)),1,0, 12,12,13,
                              1,0.0,1.0,NAME(12),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(13)),1,0, 13,13,14,
                              1,0.0,1.0,NAME(13),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(14)),1,4, 14,14,15,
                              1,0.0,1.0,NAME(14),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(15)),1,4, 15,15,16,
                              1,0.0,1.0,NAME(15),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(16)),1,4, 16,16,17,
                              1,0.0,1.0,NAME(16),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(17)),1,4, 17,17,18,
                              1,0.0,1.0,NAME(17),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(18)),1,4, 18,18,19,
                              1,0.0,1.0,NAME(18),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(19)),1,4, 19,19,20,
                              1,0.0,1.0,NAME(19),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(20)),1,4, 20,20,21,
                              1,0.0,1.0,NAME(20),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(21)),1,4, 21,21,22,
                              1,0.0,1.0,NAME(21),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(22)),1,4, 22,22,23,
                              1,0.0,1.0,NAME(22),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(23)),1,4, 23,23,24,
                              1,0.0,1.0,NAME(23),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(24)),1,0, 24,24,25,
                              1,0.0,1.0,NAME(24),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(25)),1,0, 25,25,0,
                              1,0.0,1.0,NAME(25),NULL)==NULL)
    return(1);

  /* die drei abschliessenden Kanten */
  if (CreateBoundarySegment2D(QUOTE(NAME(26)),2,0, 26,1,8,
                              1,0.0,1.0,NAME(26),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(34)),3,0, 27,9,12,
                              1,0.0,1.0,NAME(34),NULL)==NULL)
    return(1);
  if (CreateBoundarySegment2D(QUOTE(NAME(38)),4,0, 28,14,24,
                              1,0.0,1.0,NAME(38),NULL)==NULL)
    return(1);

  return(0);
}

/* configuring a domain */

INT STD_BVP_Configure (INT argc, char **argv)
{
  STD_BVP *theBVP;
  DOMAIN *theDomain;
  char BVPName[NAMESIZE];
  char DomainName[NAMESIZE];
  DOUBLE dalpha;
  INT i;

  /* get BVP name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" configure %",NAMELENSTR,"[ -~]")),
              BVPName)!=1) || (strlen(BVPName)==0)) {
    for (i=0; i<argc; i++)
      if (argv[i][0] == 'b')
        if ((sscanf(argv[i],expandfmt(CONCAT3("b %",NAMELENSTR,"[ -~]")),
                    BVPName)==1) && (strlen(BVPName)>0))
          break;
    if (i >= argc) RETURN(1);
  }

  theBVP = (STD_BVP *) BVP_GetByName(BVPName);
  if (theBVP == NULL)
    RETURN(1);

  for (i=0; i<argc; i++)
    if ((argv[i][0] == 'd') && (argv[i][1] == ' '))
      if ((sscanf(argv[i],expandfmt(CONCAT3("d %",NAMELENSTR,"[ -~]")),
                  DomainName)!=1) || (strlen(DomainName)==0))
        continue;

  if (strcmp(DomainName,"Quadrilateral") == 0) {
    if (ReadAndPrintArgvPosition("x0",argc,argv,x_quad[0])) {
      x_quad[0][0] = 0.0;
      x_quad[0][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x1",argc,argv,x_quad[1])) {
      x_quad[1][0] = 1.0;
      x_quad[1][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x2",argc,argv,x_quad[2])) {
      x_quad[2][0] = 1.0;
      x_quad[2][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x3",argc,argv,x_quad[3])) {
      x_quad[3][0] = 0.0;
      x_quad[3][1] = 1.0;
    }
    if (ReadArgvDOUBLE("alpha",&alpha,argc,argv))
      alpha = 0.0;
  }
  else if (strcmp(DomainName,"Two") == 0) {
    if (ReadAndPrintArgvPosition("x0",argc,argv,x_quad[0]))
    {
      x_quad[0][0] = 0.0;
      x_quad[0][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x1",argc,argv,x_quad[1]))
    {
      x_quad[1][0] = 1.0;
      x_quad[1][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x2",argc,argv,x_quad[2]))
    {
      x_quad[2][0] = 1.0;
      x_quad[2][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x3",argc,argv,x_quad[3]))
    {
      x_quad[3][0] = 0.0;
      x_quad[3][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x4",argc,argv,x_quad[4]))
    {
      x_quad[4][0] = 2.0;
      x_quad[4][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x5",argc,argv,x_quad[5]))
    {
      x_quad[5][0] = 2.0;
      x_quad[5][1] = 1.0;
    }
  }
  else if (strcmp(DomainName,"Four") == 0) {
    if (ReadAndPrintArgvPosition("x0",argc,argv,x_quad[0]))
    {
      x_quad[0][0] = 0.0;
      x_quad[0][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x1",argc,argv,x_quad[1]))
    {
      x_quad[1][0] = 1.0;
      x_quad[1][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x2",argc,argv,x_quad[2]))
    {
      x_quad[2][0] = 1.0;
      x_quad[2][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x3",argc,argv,x_quad[3]))
    {
      x_quad[3][0] = 0.0;
      x_quad[3][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x4",argc,argv,x_quad[4]))
    {
      x_quad[4][0] = 2.0;
      x_quad[4][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x5",argc,argv,x_quad[5]))
    {
      x_quad[5][0] = 2.0;
      x_quad[5][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x6",argc,argv,x_quad[6]))
    {
      x_quad[6][0] = 0.0;
      x_quad[6][1] = -1.0;
    }
    if (ReadAndPrintArgvPosition("x7",argc,argv,x_quad[7]))
    {
      x_quad[7][0] = 1.0;
      x_quad[7][1] = -1.0;
    }
    if (ReadAndPrintArgvPosition("x8",argc,argv,x_quad[8]))
    {
      x_quad[8][0] = 2.0;
      x_quad[8][1] = -1.0;
    }
  }
  else if (strcmp(DomainName,"Four_cr") == 0) {
    if (ReadAndPrintArgvPosition("x0",argc,argv,x_quad[0]))
    {
      x_quad[0][0] = 0.0;
      x_quad[0][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x1",argc,argv,x_quad[1]))
    {
      x_quad[1][0] = 1.0;
      x_quad[1][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x2",argc,argv,x_quad[2]))
    {
      x_quad[2][0] = 1.0;
      x_quad[2][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x3",argc,argv,x_quad[3]))
    {
      x_quad[3][0] = 0.0;
      x_quad[3][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x4",argc,argv,x_quad[4]))
    {
      x_quad[4][0] = 2.0;
      x_quad[4][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x5",argc,argv,x_quad[5]))
    {
      x_quad[5][0] = 2.0;
      x_quad[5][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x6",argc,argv,x_quad[6]))
    {
      x_quad[6][0] = 0.0;
      x_quad[6][1] = -1.0;
    }
    if (ReadAndPrintArgvPosition("x7",argc,argv,x_quad[7]))
    {
      x_quad[7][0] = 1.0;
      x_quad[7][1] = -1.0;
    }
    if (ReadAndPrintArgvPosition("x8",argc,argv,x_quad[8]))
    {
      x_quad[8][0] = 2.0;
      x_quad[8][1] = -1.0;
    }
  }
  else if (strcmp(DomainName,"Triangle") == 0)
  {
    if (ReadAndPrintArgvPosition("x0",argc,argv,x_quad[0]))
    {
      x_quad[0][0] = 0.0;
      x_quad[0][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x1",argc,argv,x_quad[1]))
    {
      x_quad[1][0] = 1.0;
      x_quad[1][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x2",argc,argv,x_quad[2]))
    {
      x_quad[2][0] = 1.0;
      x_quad[2][1] = 1.0;
    }
  }
  else if (strcmp(DomainName,"Punctured Disc") == 0)
  {
    if (ReadAndPrintArgvPosition("x0",argc,argv,x_quad[0]))
    {
      x_quad[0][0] = 0.0;
      x_quad[0][1] = 0.0;
    }
    x_quad[1][0] = 10.0;
    x_quad[1][1] = 0.0;
    x_quad[2][0] = 10.0;
    x_quad[2][1] = 10.0;
    x_quad[3][0] = 0.0;
    x_quad[3][1] = 10.0;
  }
  else if (strcmp(DomainName,"Variable Disc") == 0)
  {
    if (ReadArgvDOUBLE("l",&left,argc,argv))
    {
      left = 10.0;
    }
    if (ReadArgvDOUBLE("t",&top,argc,argv))
    {
      top = 10.0;
    }
    if (ReadArgvDOUBLE("r",&rad1,argc,argv))
    {
      rad1 = 1.0;
    }
  }
  else if (strcmp(DomainName,"Rings") == 0) {
    if (ReadArgvDOUBLE("r",&rad1,argc,argv)) {
      rad1 = INNER_RADIUS2;
    }
    if (ReadArgvDOUBLE("dalpha",&dalpha,argc,argv) == 0) {
      alpha += dalpha;
    }
    else if (ReadArgvDOUBLE("alpha",&alpha,argc,argv)) {
      alpha = 0.0;
    }
  }
  else if (strcmp(DomainName,"Rings1") == 0) {
    if (ReadArgvDOUBLE("r",&rad1,argc,argv)) {
      rad1 = INNER_RADIUS2;
    }
    if (ReadArgvDOUBLE("dalpha",&dalpha,argc,argv) == 0) {
      alpha += dalpha;
    }
    else if (ReadArgvDOUBLE("alpha",&alpha,argc,argv)) {
      alpha = 0.0;
    }
  }
  else if (strcmp(DomainName,"Rings2") == 0) {
    if (ReadArgvDOUBLE("r",&rad1,argc,argv)) {
      rad1 = INNER_RADIUS2;
    }
    if (ReadArgvDOUBLE("dalpha",&dalpha,argc,argv) == 0) {
      alpha += dalpha;
    }
    else if (ReadArgvDOUBLE("alpha",&alpha,argc,argv)) {
      alpha = 0.0;
    }
  }
  else if (strcmp(DomainName,"Skin") == 0) {
    if (ReadArgvDOUBLE("L",&L,argc,argv))
    {
      L = 10.0;
    }
    if (ReadArgvDOUBLE("D",&D,argc,argv))
    {
      D = 5.0 ;
    }
    if (ReadArgvDOUBLE("h",&glob_h,argc,argv))
    {
      glob_h = 0.5 ;
    }
  }
  else if (strcmp(DomainName,"Skin1") == 0) {
    if (ReadArgvDOUBLE("L",&L,argc,argv))
    {
      L = 10.0;
    }
    if (ReadArgvDOUBLE("D",&D,argc,argv))
    {
      D = 5.0 ;
    }
    if (ReadArgvDOUBLE("h",&glob_h,argc,argv))
    {
      glob_h = 0.5 ;
    }
  }
  else if (strcmp(DomainName,"Beam") == 0) {
    if (ReadAndPrintArgvPosition("x0",argc,argv,x_quad[0]))
    {
      x_quad[0][0] = 0.0;
      x_quad[0][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x1",argc,argv,x_quad[1]))
    {
      x_quad[1][0] = 1.0;
      x_quad[1][1] = 0.0;
    }
    if (ReadAndPrintArgvPosition("x2",argc,argv,x_quad[2]))
    {
      x_quad[2][0] = 1.0;
      x_quad[2][1] = 1.0;
    }
    if (ReadAndPrintArgvPosition("x3",argc,argv,x_quad[3]))
    {
      x_quad[3][0] = 2.0;
      x_quad[3][1] = 2.0;
    }
    if (ReadAndPrintArgvPosition("x4",argc,argv,x_quad[4]))
    {
      x_quad[4][0] = 2.0;
      x_quad[4][1] = 3.0;
    }
    if (ReadAndPrintArgvPosition("x5",argc,argv,x_quad[5]))
    {
      x_quad[5][0] = 1.0;
      x_quad[5][1] = 3.0;
    }
    if (ReadAndPrintArgvPosition("x6",argc,argv,x_quad[6]))
    {
      x_quad[6][0] = 1.0;
      x_quad[6][1] = 4.0;
    }
    if (ReadAndPrintArgvPosition("x7",argc,argv,x_quad[7]))
    {
      x_quad[7][0] = 0.0;
      x_quad[7][1] = 4.0;
    }
    if (ReadArgvDOUBLE("form",&form,argc,argv))
    {
      form = 1.0;
    }
  }

  theDomain = GetDomain(DomainName);

  if (theDomain == NULL)
  {
    if (strcmp(DomainName,"Circle") == 0)
    {
      if (InitCircle())
        return(1);
    }
    else if (strcmp(DomainName,"Ring") == 0)
    {
      if (InitRing())
        return(1);
    }
    else if (strcmp(DomainName,"Quadrilateral") == 0)
    {
      if (InitQuadrilateral())
        return(1);
    }
    else if (strcmp(DomainName,"Two") == 0)
    {
      if (InitTwo())
        return(1);
    }
    else if (strcmp(DomainName,"Four") == 0)
    {
      if (InitFour())
        return(1);
    }
    else if (strcmp(DomainName,"Four_cr") == 0)
    {
      if (InitFour_cr())
        return(1);
    }
    else if (strcmp(DomainName,"Triangle") == 0)
    {
      if (InitTriangle())
        return(1);
    }
    else if (strcmp(DomainName,"Hexagon") == 0)
    {
      if (InitHexagon())
        return(1);
    }
    else if (strcmp(DomainName,"Segment") == 0)
    {
      if (InitSegment())
        return(1);
    }
    else if (strcmp(DomainName,"Kiel") == 0)
    {
      if (InitKiel())
        return(1);
    }
    else if (strcmp(DomainName,"Wolfgangsee") == 0)
    {
      if (InitWolfgangsee())
        return(1);
    }
    else if (strcmp(DomainName,"Punctured Disc") == 0)
    {
      if (InitPuncturedDisc())
        return(1);
    }
    else if (strcmp(DomainName,"CT Disc") == 0)
    {
      if (InitCTDisc())
        return(1);
    }
    else if (strcmp(DomainName,"Variable Disc") == 0)
    {
      if (InitVDisc())
        return(1);
    }
    else if (strcmp(DomainName,"Rings") == 0)
    {
      if (InitRings())
        return(1);
    }
    else if (strcmp(DomainName,"Rings1") == 0)
    {
      if (InitRings1())
        return(1);
    }
    else if (strcmp(DomainName,"Rings2") == 0)
    {
      if (InitRings2())
        return(1);
    }
    else if (strcmp(DomainName,"Holes") == 0)
    {
      if (InitHoles())
        return(1);
    }
    else if (strcmp(DomainName,"Holes2") == 0)
    {
      if (InitHoles2())
        return(1);
    }
    else if (strcmp(DomainName,"Holes3") == 0)
    {
      if (InitHoles3())
        return(1);
    }
    else if (strcmp(DomainName,"Holes4") == 0)
    {
      if (InitHoles4())
        return(1);
    }
    else if (strcmp(DomainName,"Holes5") == 0)
    {
      if (InitHoles5())
        return(1);
    }
    else if (strcmp(DomainName,"Holes6") == 0)
    {
      if (InitHoles6())
        return(1);
    }
    else if (strcmp(DomainName,"Skin") == 0)
    {
      if (InitSkin())
        return(1);
    }
    else if (strcmp(DomainName,"Skin1") == 0)
    {
      if (InitSkin1())
        return(1);
    }
    else if (strcmp(DomainName,"Composed1") == 0)
    {
      if (InitComposed1())
        return(1);
    }
    else if (strcmp(DomainName,"Composed1a") == 0)
    {
      if (InitComposed1a())
        return(1);
    }
    else if (strcmp(DomainName,"Composed2") == 0)
    {
      if (InitComposed2())
        return(1);
    }
    else if (strcmp(DomainName,"Composed2a") == 0)
    {
      if (InitComposed2a())
        return(1);
    }
    else if (strcmp(DomainName,"Composed3") == 0)
    {
      if (InitComposed3())
        return(1);
    }
    else if (strcmp(DomainName,"Composed3a") == 0)
    {
      if (InitComposed3a())
        return(1);
    }
    else if (strcmp(DomainName,"Beam") == 0)
    {
      if (InitBeam())
        return(1);
    }
    else if (strcmp(DomainName,"Channel") == 0)
    {
      if (InitChannel())
        return(1);
    }
    else if (strcmp(DomainName,"ChannelNoParts") == 0)
    {
      if (InitChannelNoParts())
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
