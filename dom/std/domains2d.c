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
static DOUBLE alpha;

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
  if (CreateDomain("Two",MidPoint,radius,8,8,YES)==NULL) return(1);
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
  if (CreateBoundarySegment2D("left", 1,0,
                              4,0,4,1,
                              0.0,1.0,
                              LeftBoundary, NULL)==NULL) return(1);

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
/*D
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
   D*/
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

  if (strcmp(DomainName,"Quadrilateral") == 0)
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
    if (ReadAndPrintArgvPosition("x3",argc,argv,x_quad[3]))
    {
      x_quad[3][0] = 0.0;
      x_quad[3][1] = 1.0;
    }
    if (ReadArgvDOUBLE("alpha",&alpha,argc,argv))
      alpha = 0.0;
  }
  else if (strcmp(DomainName,"Two") == 0)
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
    else
      return(1);

    theDomain = GetDomain(DomainName);

    if (theDomain == NULL)
      return(1);
  }

  theBVP->Domain = theDomain;

  return(0);
}
