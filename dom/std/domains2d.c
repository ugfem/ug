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
#include "devices.h"

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

static DOUBLE_VECTOR x_quad[6];
static DOUBLE alpha,left,top,rad1;

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
static const INT two_pt2p[8] = {0,1,1,0,3,3,3,3};
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
  radius = 5.0;
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
                   MidPoint,radius,5,5,NO)==NULL) return(1);

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
                              2,2,3,1,
                              0.0,1.0,
                              BVar3Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bvar4", 1,0,
                              3,3,4,1,
                              0.0,1.0,
                              BVar4Boundary, NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("bvar5", 1,0,
                              4,4,0,1,
                              0.0,1.0,
                              BVar5Boundary, NULL)==NULL) return(1);

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
  radius =110.0;
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

  /* allocate new domain structure */
  MidPoint[0] = 2.5;
  MidPoint[1] = 1.5;
  radius = 3.0;
  if (CreateDomainWithParts("Holes",MidPoint,radius,52,52,NO,2,&h_dpi)
      ==NULL)
    return(1);

  if (CreateBoundarySegment2D("start1" ,1,0,0,0,1,1.0,0.0,1.0,
                              Start1Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start2" ,1,0,1,1,2,1.0,0.0,1.0,
                              Start2Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start3" ,1,0,2,2,3,1.0,0.0,1.0,
                              Start3Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start4" ,1,0,3,3,0,1.0,0.0,1.0,
                              Start4Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start5" ,1,0,4,4,5,1.0,0.0,1.0,
                              Start5Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start6" ,1,0,5,5,6,1.0,0.0,1.0,
                              Start6Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start7" ,1,0,6,6,7,1.0,0.0,1.0,
                              Start7Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start8" ,1,0,7,7,4,1.0,0.0,1.0,
                              Start8Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start9" ,1,0,8,8,9,1.0,0.0,1.0,
                              Start9Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start10",1,0,9,9,10,1.0,0.0,1.0,
                              Start10Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start11",1,0,10,10,11,1.0,0.0,1.0,
                              Start11Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start12",1,0,11,11,8,1.0,0.0,1.0,
                              Start12Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start13",1,0,12,12,13,1.0,0.0,1.0,
                              Start13Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start14",1,0,13,13,14,1.0,0.0,1.0,
                              Start14Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start15",1,0,14,14,15,1.0,0.0,1.0,
                              Start15Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start16",1,0,15,15,12,1.0,0.0,1.0,
                              Start16Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start17",1,0,16,16,17,1.0,0.0,1.0,
                              Start17Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start18",1,0,17,17,18,1.0,0.0,1.0,
                              Start18Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start19",1,0,18,18,19,1.0,0.0,1.0,
                              Start19Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start20",1,0,19,19,16,1.0,0.0,1.0,
                              Start20Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start21",1,0,20,20,21,1.0,0.0,1.0,
                              Start21Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start22",1,0,21,21,22,1.0,0.0,1.0,
                              Start22Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start23",1,0,22,22,23,1.0,0.0,1.0,
                              Start23Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start24",1,0,23,23,20,1.0,0.0,1.0,
                              Start24Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start25",1,0,24,24,25,1.0,0.0,1.0,
                              Start25Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start26",1,0,25,25,26,1.0,0.0,1.0,
                              Start26Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start27",1,0,26,26,27,1.0,0.0,1.0,
                              Start27Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start28",1,0,27,27,24,1.0,0.0,1.0,
                              Start28Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start29",1,0,28,28,29,1.0,0.0,1.0,
                              Start29Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start30",1,0,29,29,30,1.0,0.0,1.0,
                              Start30Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start31",1,0,30,30,31,1.0,0.0,1.0,
                              Start31Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start32",1,0,31,31,28,1.0,0.0,1.0,
                              Start32Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start33",1,0,32,32,33,1.0,0.0,1.0,
                              Start33Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start34",1,0,33,33,34,1.0,0.0,1.0,
                              Start34Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start35",1,0,34,34,35,1.0,0.0,1.0,
                              Start35Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start36",1,0,35,35,32,1.0,0.0,1.0,
                              Start36Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start37",1,0,36,36,37,1.0,0.0,1.0,
                              Start37Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start38",1,0,37,37,38,1.0,0.0,1.0,
                              Start38Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start39",1,0,38,38,39,1.0,0.0,1.0,
                              Start39Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start40",1,0,39,39,36,1.0,0.0,1.0,
                              Start40Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start41",1,0,40,40,41,1.0,0.0,1.0,
                              Start41Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start42",1,0,41,41,42,1.0,0.0,1.0,
                              Start42Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start43",1,0,42,42,43,1.0,0.0,1.0,
                              Start43Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start44",1,0,43,43,40,1.0,0.0,1.0,
                              Start44Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start45",1,0,44,44,45,1.0,0.0,1.0,
                              Start45Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start46",1,0,45,45,46,1.0,0.0,1.0,
                              Start46Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start47",1,0,46,46,47,1.0,0.0,1.0,
                              Start47Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start48",1,0,47,47,44,1.0,0.0,1.0,
                              Start48Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start49",1,0,48,48,49,1.0,0.0,1.0,
                              Start49Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start50",1,0,49,49,50,1.0,0.0,1.0,
                              Start50Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start51",1,0,50,50,51,1.0,0.0,1.0,
                              Start51Boundary,NULL)==NULL) return(1);
  if (CreateBoundarySegment2D("start52",1,0,51,51,48,1.0,0.0,1.0,
                              Start52Boundary,NULL)==NULL) return(1);

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

#define INNEN_RADIUS 0.6

static INT kreisBoundaryUpper1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = INNEN_RADIUS * cos(PI*lambda);
  result[1] = INNEN_RADIUS * sin(PI*lambda);

  return(0);
}

static INT kreisBoundaryLower1 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = INNEN_RADIUS * cos(PI+PI*lambda);
  result[1] = INNEN_RADIUS * sin(PI+PI*lambda);

  return(0);
}

static INT kreisBoundaryUpper2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = INNEN_RADIUS * cos(alpha + PI*lambda);
  result[1] = INNEN_RADIUS * sin(alpha + PI*lambda);

  return(0);
}

static INT kreisBoundaryLower2 (void *data, DOUBLE *param, DOUBLE *result)
{
  DOUBLE lambda;

  lambda = param[0];

  if ((lambda<0.0)||(lambda>1.0)) return(1);

  result[0] = INNEN_RADIUS * cos(alpha + PI + PI*lambda);
  result[1] = INNEN_RADIUS * sin(alpha + PI + PI*lambda);

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

  if (CreateDomainWithParts("Rings",MidPoint,radius,6,6,YES,4,&ring_dpi)
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
  else if (strcmp(DomainName,"BVariabel") == 0)
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
    if (ReadArgvDOUBLE("dalpha",&dalpha,argc,argv) == 0) {
      alpha += dalpha;
    }
    else if (ReadArgvDOUBLE("alpha",&alpha,argc,argv)) {
      alpha = 0.0;
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
    else if (strcmp(DomainName,"Holes") == 0)
    {
      if (InitHoles())
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
