// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lgm_domain3d.c												*/
/*																			*/
/* Purpose:   source for lgm_domain                                                                             */
/*																			*/
/* Author:	  Christoph Tapp                                                                                                */
/*			  Institut fuer Angewandte Mathematik	                                                */
/*			  Universitaet Erlangen											*/
/*			  Martenstr. 3													*/
/*			  91058 Erlangen												*/
/*			  email: tapp@am.uni-erlangen.de								*/
/*																			*/
/* History:   17.2.1997														*/
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

#include <config.h>

#undef OCC_GEOMETRY
#ifdef OCC_GEOMETRY
#include "occ/occ_geom.hh"
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "ugtypes.h"
#include "architecture.h"
#include "heaps.h"
#include "bio.h"
#include "domain.h"
#include "general.h"
#include "misc.h"
#include "evm.h"
#include "lgm_domain.h"
#include "lgm_load.h"
#include "lgm_macros.h"
#include "netgen/lgm_gginterface.h"

#include "ugdevices.h"

#include "namespace.h"

USING_UG_NAMESPACE
  USING_UGDIM_NAMESPACE


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define LGM_BUFFERLEN                                           128

#define BVP2LGM(p)                                                      ((LGM_DOMAIN*)(p))
#define BNDP2LGM(p)                                                     ((LGM_BNDP*)(p))
#define BNDS2LGM(p)                                                     ((LGM_BNDS*)(p))

#ifdef OCC_GEOMETRY
extern int *LineMap;
extern int *SurfMap;
extern OCC_GEOM occ_geom;
#endif


/* global mark key for temp mem */
static INT LGM_MarkKey;

static INT SurfaceInfoId;
static INT LineInfoId;
static INT LGM_DEBUG = 0;
static INT SAVE_SURFACE;
static INT SURFACE_ADAPT;
static DOUBLE LINE_DISTANCE = 0.05;

#define SMALL_FAK 1e-6
static DOUBLE SMALL;

//static INT VAR_H = 1; /* OS_CHANGED: not used */

//#define SMALL1 1e-6 OS_CHANGED
#define SMALL1 SMALL_FAK

/*static DOUBLE cosAngle = 0.99;		komplex3d */
/*static DOUBLE cosAngle = 0.99999;	*/
static DOUBLE cosAngle;         /* Winkel zwischen Inputdreiecken */
static DOUBLE Triangle_Angle2 = 40.0;
static DOUBLE EPS = 0.1;

/* INT NS_DIM_PREFIX Surface_Local2Global (LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local);*/
INT GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local, DOUBLE *n);

static INT Get_NBNDS_Per_Subdomain              (HEAP *Heap, LGM_DOMAIN *theDomain, INT **sizes, DOUBLE h);
static INT DiscretizeLine                               (HEAP *Heap, LGM_LINE *theLine, DOUBLE h, LGM_POINT *pointlist, INT norp, INT MarkKey);
static INT DiscretizeLineNew                    (HEAP *Heap, LGM_LINE *theLine, DOUBLE h, LGM_POINT *pointlist, INT norp, INT MarkKey);
static INT Get_NBNDP                                    (LGM_DOMAIN *theDomain, INT *nBND, DOUBLE h);
static INT DiscretizeDomain                     (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h, INT MarkKey);
static INT DiscretizeSurface                    (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h, LGM_POINT *pointlist,
                                                 INT norp, INT D, INT MarkKey);
static INT TransferSurfaces2Mesh                (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h);
static INT P_Dist(DOUBLE *p1, DOUBLE *p2);
static DOUBLE E_Distance(DOUBLE *p1, DOUBLE *p2);


INT TEST(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local);

extern INT GenerateSurfaceGrid (HEAP *theHeap, INT MarkKey, LGM_SURFACE *aSurface, DOUBLE h, INT smooth,INT display,INT D);
extern INT InitSurface(CoeffProcPtr Coeff);
extern int AddGeomElement (int node0, int node1, int node2, int neigbor0, int neigbor1, int neigbor2);

/* data for CVS */
/*static char RCS_ID("$Header$",UG_RCS_STRING);
 */
static INT currSubdom, currSurface, currLine;

static CoeffProcPtr Coefficients[8];
static CoeffProcPtr LOCAL_H;

static LGM_SURFACE *xy_Surface;
static LGM_LINE *xy_Line[4];
static LGM_POINT *xy_Point[4];
static INT xy_bpoints;
static DOUBLE xy_n[3];
static DOUBLE **xy_lcoord;
static INT *xy_flag;
static INT *xy_new_id;
static INT *xy_element_id;
static INT surface_error;
static INT surf_err;

static INT P_Dist(DOUBLE *p1, DOUBLE *p2)
{
  DOUBLE x, y, z, r;

  x = p1[0] - p2[0];
  y = p1[1] - p2[1];
  z = p1[2] - p2[2];

  r = sqrt(x*x+y*y+z*z);
  if(r<SMALL)
    return(1);
  else
    return(0);
}

static DOUBLE E_Distance(DOUBLE *p1, DOUBLE *p2)
{
  DOUBLE d;

  d = sqrt( (p1[0]-p2[0])*(p1[0]-p2[0])
            +     (p1[1]-p2[1])*(p1[1]-p2[1])
            +     (p1[2]-p2[2])*(p1[2]-p2[2]) );
  return(d);
}



static INT ResetSurfaceFlags (LGM_DOMAIN *theDomain)
{
  INT i,j;
  LGM_SUBDOMAIN *theSubdom;

  for (i=1; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    theSubdom = LGM_DOMAIN_SUBDOM(theDomain,i);
    for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSubdom); j++)
      LGM_SURFACE_FLAG(LGM_SUBDOMAIN_SURFACE(theSubdom,j)) = 0;
  }

  return (0);
}

LGM_SURFACE *NS_DIM_PREFIX FirstSurface (LGM_DOMAIN *theDomain)
{
  LGM_SURFACE *theSurface;

  if (ResetSurfaceFlags(theDomain))
    return (NULL);
  theSurface = LGM_SUBDOMAIN_SURFACE(LGM_DOMAIN_SUBDOM(theDomain,1),0);
  LGM_SURFACE_FLAG(theSurface) = 1;
  currSurface = 0;
  currSubdom = 1;
  return (theSurface);
}

static LGM_SURFACE *helpNextSurface (LGM_DOMAIN *theDomain)
{
  LGM_SUBDOMAIN *theSubdom;
  LGM_SURFACE *theSurface;

  theSubdom = LGM_DOMAIN_SUBDOM(theDomain,currSubdom);
  if (currSurface<LGM_SUBDOMAIN_NSURFACE(theSubdom)-1)
  {
    currSurface++;
    theSurface =  LGM_SUBDOMAIN_SURFACE(theSubdom,currSurface);
    return (theSurface);
  }
  else if (currSubdom<LGM_DOMAIN_NSUBDOM(theDomain))
  {
    currSubdom++;
    theSubdom = LGM_DOMAIN_SUBDOM(theDomain,currSubdom);
    currSurface = 0;
    theSurface =  LGM_SUBDOMAIN_SURFACE(theSubdom,currSurface);
    return (theSurface);
  }

  return (NULL);
}

LGM_SURFACE *NS_DIM_PREFIX NextSurface (LGM_DOMAIN *theDomain)
{
  LGM_SURFACE *theSurface;

  while(1)
  {
    theSurface = helpNextSurface(theDomain);
    if (theSurface==NULL) \
      break;
    if (LGM_SURFACE_FLAG(theSurface)==0)
    {
      LGM_SURFACE_FLAG(theSurface) = 1;
      break;
    }
  }
  return (theSurface);
}

static INT ResetLineFlags (LGM_DOMAIN *theDomain)
{
  INT i,j,k;
  LGM_SUBDOMAIN *theSubdom;
  LGM_SURFACE *theSurface;

  for (i=1; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    theSubdom = LGM_DOMAIN_SUBDOM(theDomain,i);
    for(j=0; j<LGM_SUBDOMAIN_NSURFACE(theSubdom); j++)
    {
      theSurface = LGM_SUBDOMAIN_SURFACE(theSubdom,j);
      for(k=0; k<LGM_SURFACE_NLINE(theSurface); k++)
      {
        LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface,k)) = 0;
      }
    }
  }

  return (0);
}

LGM_LINE *NS_DIM_PREFIX
FirstLine (LGM_DOMAIN *theDomain)
{
  LGM_LINE *theLine;

  if (ResetLineFlags(theDomain))
    return (NULL);
  theLine = LGM_SURFACE_LINE(LGM_SUBDOMAIN_SURFACE(LGM_DOMAIN_SUBDOM(theDomain,1),0),0);;
  LGM_LINE_FLAG(theLine) = 1;
  currLine = 0;
  currSurface = 0;
  currSubdom = 1;
  return (theLine);
}

static LGM_LINE *helpNextLine (LGM_DOMAIN *theDomain)
{
  LGM_SUBDOMAIN *theSubdom;
  LGM_SURFACE * theSurface;
  LGM_LINE *theLine;

  theSubdom = LGM_DOMAIN_SUBDOM(theDomain,currSubdom);
  theSurface = LGM_SUBDOMAIN_SURFACE(theSubdom,currSurface);
  if (currLine<LGM_SURFACE_NLINE(theSurface)-1)
  {
    currLine++;
    theLine = LGM_SURFACE_LINE(theSurface,currLine);
    return (theLine);
  }
  else if (currSurface<LGM_SUBDOMAIN_NSURFACE(theSubdom)-1)
  {
    currSurface++;
    theSurface = LGM_SUBDOMAIN_SURFACE(theSubdom,currSurface);
    currLine = 0;
    theLine =  LGM_SURFACE_LINE(theSurface,currLine);
    return (theLine);
  }
  else if (currSubdom<LGM_DOMAIN_NSUBDOM(theDomain))
  {
    currSubdom++;
    currSurface = 0;
    currLine = 0;
    theSurface = LGM_SUBDOMAIN_SURFACE(theSubdom,currSurface);
    theLine = LGM_SURFACE_LINE(theSurface,currLine);
    return (theLine);
  }
  return (NULL);
}

LGM_LINE *NS_DIM_PREFIX
NextLine (LGM_DOMAIN *theDomain)
{
  LGM_LINE *theLine;

  while(1)
  {
    theLine = helpNextLine(theDomain);
    if (theLine==NULL) break;
    if (LGM_LINE_FLAG(theLine)==0)
    {
      LGM_LINE_FLAG(theLine) = 1;
      break;
    }
  }
  return (theLine);
}

INT NS_DIM_PREFIX SetBoundaryCondition (LGM_DOMAIN *theDomain, BndCondProcPtr BndCond, BndCondProcPtr InnerBndCond)
{
  INT i,k;
  LGM_SUBDOMAIN *theSubdom;
  LGM_SURFACE *theSurface;

  for (i=1; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    theSubdom = LGM_DOMAIN_SUBDOM(theDomain,i);
    for (k=0; k<LGM_SUBDOMAIN_NSURFACE(theSubdom); k++)
    {
      theSurface = LGM_SUBDOMAIN_SURFACE(theSubdom,k);
      if (LGM_SURFACE_LEFT(theSurface)*LGM_SURFACE_RIGHT(theSurface)!=0)
        LGM_SURFACE_BNDCOND(theSurface) = InnerBndCond;
      else
        LGM_SURFACE_BNDCOND(theSurface) = BndCond;
    }
  }

  return (0);
}

INT NS_DIM_PREFIX SetDomainSize (LGM_DOMAIN *theDomain)
{
  LGM_PROBLEM *theProblem;
  LGM_LINE *theLine;
  DOUBLE min[3], max[3];
  INT i;

  min[0]=min[1]=min[2]=MAX_C;
  max[0]=max[1]=max[2]=-MAX_C;

  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
    for (i=0; i<LGM_LINE_NPOINT(theLine); i++)
    {
      min[0] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[0],min[0]);
      min[1] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[1],min[1]);
      min[2] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[2],min[2]);
      max[0] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[0],max[0]);
      max[1] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[1],max[1]);
      max[2] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[2],max[2]);
    }


  LGM_DOMAIN_MIDPOINT(theDomain)[0] = 0.5*(min[0]+max[0]);
  LGM_DOMAIN_MIDPOINT(theDomain)[1] = 0.5*(min[1]+max[1]);
  LGM_DOMAIN_MIDPOINT(theDomain)[2] = 0.5*(min[2]+max[2]);
  LGM_DOMAIN_RADIUS(theDomain) = 0.55*sqrt((max[0]-min[0])*(max[0]-min[0])
                                           +(max[1]-min[1])*(max[1]-min[1])
                                           +(max[2]-min[2])*(max[2]-min[2]));

  theProblem = LGM_DOMAIN_PROBLEM(theDomain);
  if (LGM_PROBLEM_DOMCONFIG(theProblem)!=NULL)
    if ((*LGM_PROBLEM_DOMCONFIG (theProblem))(min,max))
      return (1);

  SMALL = SMALL_FAK * LGM_DOMAIN_RADIUS(theDomain);
  return (0);
}

static INT PrintDomainInfo (LGM_DOMAIN *aDomain)
{

  printf("********* domain-info *********\n");

  printf("%s %s\n","Name: ",LGM_DOMAIN_PROBLEMNAME(aDomain));
  printf("%s %d\n","nSubdomain: ",(int)LGM_DOMAIN_NSUBDOM(aDomain));
  printf("%s %d\n","nPoint: ",(int)LGM_DOMAIN_NPOINT(aDomain));
  printf("%s %f\n","radius: ",LGM_DOMAIN_RADIUS(aDomain));
  printf("%s %f %f %f\n","midpoint: ",LGM_DOMAIN_MIDPOINT(aDomain)[0],
         LGM_DOMAIN_MIDPOINT(aDomain)[1],
         LGM_DOMAIN_MIDPOINT(aDomain)[2]);
  printf("%s %d\n","convex: ",(int)LGM_DOMAIN_CONVEX(aDomain));

  return (0);
}

static INT PrintSubdomainInfo (LGM_DOMAIN *aDomain)
{
  INT i;
  LGM_SUBDOMAIN *aSubdom;

  printf("********* subdomain-info *********\n");
  for(i=1; i<=LGM_DOMAIN_NSUBDOM(aDomain); i++)
  {
    aSubdom = LGM_DOMAIN_SUBDOM(aDomain,i);
    printf("%s %d\n","id: ",(int)LGM_SUBDOMAIN_ID(aSubdom));
    printf("%s %d\n","nSurface: ",(int)LGM_SUBDOMAIN_NSURFACE(aSubdom));
  }
  return (0);
}

static INT PrintSurfaceInfo (LGM_SURFACE *aSurface)
{
  INT i;

  printf("********* surface-info *********\n");
  printf("%s %d\n","SurfaceId: ",(int)SurfaceInfoId);
  SurfaceInfoId++;
  printf("%s %d\n","nPoint: ",(int)LGM_SURFACE_NPOINT(aSurface));
  printf("%s %d\n","nTriangle: ",(int)LGM_SURFACE_NTRIANGLE(aSurface));
  printf("%s %d\n","nLine: ",(int)LGM_SURFACE_NLINE(aSurface));
  printf("%s %d\n","left: ",(int)LGM_SURFACE_LEFT(aSurface));
  printf("%s %d\n","right: ",(int)LGM_SURFACE_RIGHT(aSurface));
  /*	printf("%s %d\n","BndCond: ",LGM_SURFACE_BNDCOND(aSurface));*/
  for(i=0; i<LGM_SURFACE_NPOINT(aSurface); i++)
    printf("%s %f %f %f\n","Point: ",LGM_SURFACE_POINT(aSurface,i)->position[0],
           LGM_SURFACE_POINT(aSurface,i)->position[1],
           LGM_SURFACE_POINT(aSurface,i)->position[2]);
  for(i=0; i<LGM_SURFACE_NTRIANGLE(aSurface); i++)
  {
    printf("%s %d %d %d\n","Triangle: ",LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(aSurface,i),0),
           LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(aSurface,i),1),
           LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(aSurface,i),2));
    printf("%f %f %f\n",LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),0))[0],
           LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),0))[1],
           LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),0))[2]);
    printf("%f %f %f\n",LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),1))[0],
           LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),1))[1],
           LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),1))[2]);
    printf("%f %f %f\n",LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),2))[0],
           LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),2))[1],
           LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(aSurface,i),2))[2]);
  }

  for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(aSurface)); i++)
  {
    printf("%d %d %d\n",(int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(aSurface),i,0),
           (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(aSurface),i,1),
           (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(aSurface),i,2));
  }

  return (0);
}

static INT Line_Local2Global (LGM_LINE *theLine, DOUBLE *global, INT i)
{
  DOUBLE slocal;
  INT ilocal;

  ilocal = (int)floor(LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),i));
  slocal = LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),i) - ilocal;

  assert(slocal>=0.0);
  assert(ilocal<LGM_LINE_NPOINT(theLine) && ilocal>=0);

  if (ilocal<LGM_LINE_NPOINT(theLine)-1)
  {
    global[0] = (1.0-slocal)*LGM_LINE_POINT(theLine,ilocal)->position[0]
                + slocal*LGM_LINE_POINT(theLine,ilocal+1)->position[0];
    global[1] = (1.0-slocal)*LGM_LINE_POINT(theLine,ilocal)->position[1]
                + slocal*LGM_LINE_POINT(theLine,ilocal+1)->position[1];
    global[2] = (1.0-slocal)*LGM_LINE_POINT(theLine,ilocal)->position[2]
                + slocal*LGM_LINE_POINT(theLine,ilocal+1)->position[2];
  }
  else
  {
    if (slocal!=0.0)
    {
      UserWrite("*\n");
      UserWriteF("slocal: %f\n",(float)slocal);
    }
    assert(slocal==0.0);
    global[0] = LGM_LINE_POINT(theLine,ilocal)->position[0];
    global[1] = LGM_LINE_POINT(theLine,ilocal)->position[1];
    global[2] = LGM_LINE_POINT(theLine,ilocal)->position[2];
  }
  return(0);
}

INT NS_DIM_PREFIX Line_Local2GlobalNew (LGM_LINE *theLine, DOUBLE *global, DOUBLE local)
{
  DOUBLE slocal;
  INT ilocal;

  ilocal = (int)floor(local);
  slocal = local - ilocal;

  if(slocal<0.0)
    assert(slocal>=0.0);
  assert(ilocal<LGM_LINE_NPOINT(theLine) && ilocal>=0);

  if (ilocal<LGM_LINE_NPOINT(theLine)-1)
  {
    global[0] = (1.0-slocal)*LGM_LINE_POINT(theLine,ilocal)->position[0]
                + slocal*LGM_LINE_POINT(theLine,ilocal+1)->position[0];
    global[1] = (1.0-slocal)*LGM_LINE_POINT(theLine,ilocal)->position[1]
                + slocal*LGM_LINE_POINT(theLine,ilocal+1)->position[1];
    global[2] = (1.0-slocal)*LGM_LINE_POINT(theLine,ilocal)->position[2]
                + slocal*LGM_LINE_POINT(theLine,ilocal+1)->position[2];
  }
  else
  {
    if (slocal!=0.0)
    {
      UserWrite("*\n");
      UserWriteF("slocal: %f\n",(float)slocal);
    }
    /*	        assert(slocal==0.0);*/
    global[0] = LGM_LINE_POINT(theLine,ilocal)->position[0];
    global[1] = LGM_LINE_POINT(theLine,ilocal)->position[1];
    global[2] = LGM_LINE_POINT(theLine,ilocal)->position[2];
  }
  return(0);
}

INT Line_Global2Local (LGM_LINE *theLine, DOUBLE *global, DOUBLE *local)
{
  DOUBLE start[3], end[3], lambda[3], globalnew[3], d;
  INT i, j, id[3], r, lambda_counter;

  *local = -1.0;
  for(i=0; i<LGM_LINE_NPOINT(theLine)-1; i++)
  {
    id[0] = id[1] = id[2] = 0;
    lambda[0] = lambda[1] = lambda[2] = 0.0;
    lambda_counter = 0;
    Line_Local2GlobalNew(theLine,start,(DOUBLE)(i));
    Line_Local2GlobalNew(theLine,end,(DOUBLE)(i+1));
    for(j=0; j<3; j++)
    {
      if(sqrt( (end[j] - start[j])*(end[j] - start[j]) ) <SMALL)
      {
        if(sqrt( (end[j] - global[j])*(end[j] - global[j]) ) < SMALL)
          id[j] = 1;
      }
      else
      {
        lambda[j] = (global[j] - start[j]) / (end[j] - start[j]);
        if( (lambda[j]>=0.0)&&(lambda[j]<=1.0) )
        {
          r = j;
          id[j] = 1;
          // OS_CHANGED		l = j;
          lambda_counter++;
        }
      }
    }
    if(id[0]+id[1]+id[2] == 3)
    {
      d =     ((lambda[0]+lambda[1]+lambda[2])/lambda_counter - lambda[r])
          * ((lambda[0]+lambda[1]+lambda[2])/lambda_counter - lambda[r]);
      if(sqrt(d)<SMALL)
      {
        /* lokale Koordinaten gefunden */
        *local = (DOUBLE)i + lambda[r];
        /* check */
        Line_Local2GlobalNew(theLine, globalnew, *local);
        if(sqrt((global[0]-globalnew[0])*(global[0]-globalnew[0])
                +    (global[1]-globalnew[1])*(global[1]-globalnew[1])
                +    (global[2]-globalnew[2])*(global[2]-globalnew[2]))>SMALL)
          printf("%s\n", "Line_Global2Local ist falsch");

      }
    }
  }

  return(0);
}

static INT PrintLineInfo (LGM_LINE *aLine)
{
  INT i;
  DOUBLE global[3];
  LINEPOINT *help;

  printf("********* line-info *********\n");
  printf("%s %d\n","LineId: ",(int)LGM_LINE_ID(aLine));
  LineInfoId++;
  printf("%s %d\n","nPoint: ",(int)LGM_LINE_NPOINT(aLine));
  printf("%s %d %d\n","firstPoint lastPoint: ",(int)LGM_LINE_BEGIN(aLine),(int)LGM_LINE_END(aLine));
  for(i=0; i<LGM_LINE_NPOINT(aLine); i++)
    printf("%s %f %f %f\n","Point: ",LGM_LINE_POINT(aLine,i)->position[0],
           LGM_LINE_POINT(aLine,i)->position[1],
           LGM_LINE_POINT(aLine,i)->position[2]);
  /*	printf("%s\n","linedisc");
          printf("%s %d\n","nPoint: ",LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(aLine)));
          for(i=0;i<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(aLine));i++)
          {
                  printf("%f\n",LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(aLine),i));
                  Line_Local2Global(aLine,global,i);
                  printf("%f %f %f\n",global[0],global[1],global[2]);
          }*/
  printf("%s\n","linediscnew");
  printf("%s %d\n","nPoint: ",(int)LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(aLine)));
  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(aLine));
  for(i=0; i<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(aLine)); i++)
  {
    printf("%f\n",help->local);
    Line_Local2GlobalNew(aLine,global,help->local);
    printf("%f %f %f\n",global[0],global[1],global[2]);
    help = help->next;
  }
  return (0);
}

static INT PrintMeshInfo (MESH *mesh)
{
  INT i,j,k;
  DOUBLE global[3];

  printf("********* mesh-info *********\n");
  printf("\n");

  printf("nBNPs = %d\n",(int)mesh->nBndP);
  for (i=0; i<mesh->nBndP; i++)
  {
    if (BNDP_Global(mesh->theBndPs[i],global))
      return (1);
    printf("    BNP(%d) = ",(int)i);;
    for(j=0; j<3; j++)
      printf("%f ",(float)global[j]);
    printf("%d \n",(int)LGM_BNDP_N(BNDP2LGM(mesh->theBndPs[i])));
  }
  printf("\n");

  /*	for (i=0; i<mesh->nBndP; i++)
     {
          printf("%s \n","globali");

          for(index=0;index<LGM_BNDP_N(BNDP2LGM(mesh->theBndPs[i]));index++)
          {
                  if (BNDP_Globali(mesh->theBndPs[i],global,index))
                          return (1);
                  printf("%d ",index);
                  printf("    BNP(%d) = ",(int)i);
                  for(j=0;j<3;j++)
                          printf("%f ",(float)global[j]);
                  printf("\n");
          }
     }
     printf("\n");*/

  printf("nSub = %d\n",(int)mesh->nSubDomains);
  printf("\n");

  for (i=1; i<=mesh->nSubDomains; i++)
  {
    printf("Subdomain %d\n",(int)i);
    printf("    nSides = %d\n",(int)mesh->nSides[i]);
    for (j=0; j<mesh->nSides[i]; j++)
    {
      printf("    Side_corner_ids = ");
      for (k=0; k<mesh->Side_corners[i][j]; k++)
        printf("%d ",(int)mesh->Side_corner_ids[i][j][k]);
      printf("\n");
    }
  }
  printf("\n");
  printf("********* mesh-info *********\n");

  return (0);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BVP_Save (BVP *theBVP, const char *name, const char *mgname, HEAP *theHeap, INT argc, char **argv)
{
  UserWrite("SORRY: not implemented yet\n");
  return (1);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BVP_Check (BVP *aBVP)
{
  INT i,i2,j,j2,p,ret,sbd,sfce,left,right,l,m, m2;
  INT commmonIDs_Counter, tr_index, ntr_index, richtung_TheTria, richtung_TheNgbTria;
  LGM_POINT **TemporaryPointArray;
  INT indexesTheTria[2];
  INT indexesTheNgbTria[2];
  INT yes[3];
  LGM_DOMAIN *theD,*lgmd;
  LGM_SUBDOMAIN *theSD, *theSD2, *lfSD;
  LGM_LINE *theL,*theL2;
  LGM_SURFACE *theSF, *theSF2, *lfSF;
  LGM_TRIANGLE *TheTria, *TheTria2, *TheNgbTria;
  LGM_POINT *ThePoint, *Point_TheTria, *Point_TheNgbTria;
  INT CHECKA,CHECKB,CHECKC,CHECKD,CHECKE,CHECKF,CHECKG;
  INT NachfolgerLine;INT VorgaengerLine;
  INT memsize_variable, CornerCounter, Lf_Var;
  INT NewCorners[3];


  CHECKA = 1;
  /*Check , ob ueberhaupt alles angelegt ist*/
  /* check subdomains */

  CHECKB = 1;
  /*Checke die Left/RightInformation aller Surfaces ...*/
  /*Laufe ueber Subdomains ...*/

  CHECKC = 1;
  /*laufe ueber die Subdomains und pruefe ob es stimmt dass wirklich nur die referenzierten Sbdms
     diese Surface besitzen! else {ErrorCase} */
  /*Laufe ueber Subdomains und dort ueber Surfaces ...*/

  CHECKD = 1;
  /*CHECK aller Lines, ob sie gar zyklisch sind und ob sie eine Vorgaenger- und eine Nachfolger-Line besitzen*/

  CHECKE = 1;
  /*laufe ueber Points der Surface und suche zu jedem Point, ob es mind. ein Triangle
     gibt, das es in dieser Surface gibt, und das diesen Point besitzt.
     d.h. Pruefe ob die Surface Points besitzt, die von keinem Triangle referenziert werden.
     =>{ErrorCase}*/
  /*Laufe ueber alle Subdomains:*/

  CHECKF = 0;       /*macht nur Sinn fuer den Ansyseinlesevorgang not tested yet, 11.7.97*/
  /*Gebe die 3 NachbarTriangles an und pruefe ob es ueberhaupt ein Nachbartriangle ist und ob
     die selbe Direction vorliegt! d.h dass die beiden gemeinsamen IDs in der umgekehrten
     Reihenfolge durchlaufen werden. und pruefe auch ob 3 unterschiedliche CornerIDs vorliegen*/
  /*NEU: macht nach C.Tapps Aenderugen auch wieder SInn*/


  CHECKG = 1;
  /*Pruefe zu jedem Triangle, ob es nur einmal vorkommt*/

  UserWrite("BVP_Check3D: \n");

  theD = BVP2LGM(aBVP);
  if (theD==NULL)
  {
    UserWriteF("theD = BVP2LGM(aBVP) returned NilPointer\n\n");
    return (1);
  }

  UserWrite("Test A:  Do all components exist:");


  if(CHECKA ==1)
  {
    /*CHECKA*/
    /*Check , ob ueberhaupt alles angelegt ist*/
    /* check subdomains */
    for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
    {
      /* CHECK[0] */
      theSD = LGM_DOMAIN_SUBDOM(theD,i);
      if (theSD==NULL)
      {
        if (!ret) UserWrite("\n");
        UserWriteF("Subdomain %d is not referenced from Domain\n\n",(int)i);
        ret = 1;
        return(1);
      }

      /*Laufe ueber Surfaces*/
      for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSD); j++)
      {
        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        /*UserWriteF("The %d . Surface of Sbd %d has got %d Points",(int)j,(int)i,LGM_SURFACE_NPOINT(theSF));*/
        /*UserWriteF("and has got the ID %d \n",LGM_SURFACE_ID(theSF));*/
        for (l=0; l<LGM_SURFACE_NLINE(theSF); l++)
        {
          theL = LGM_SURFACE_LINE(theSF, l);
          /*zyklische Line ?*/
          /*UserWriteF("	The %d . Line :\n Begin: %d , End: %d\n",(int)l,LGM_LINE_BEGIN(theL),LGM_LINE_END(theL));*/
          /*UserWriteF("	The %d . Line has got the following ID: %d\n",(int)l,LGM_LINE_ID(theL));*/
        }
        if (theSF==NULL)
        {
          /* CHECK[13] */
          if (!ret) UserWrite("\n");
          UserWriteF("Surfcae %d is not referenced from Subdomain %d\n\n",(int)j,(int)i);
          ret = 1;
          return(1);
        }

        /*Laufe ueber Lines ...*/
        if(LGM_SURFACE_NLINE(theSF) <2)
        {
          if (!ret) UserWrite("\n");
          UserWriteF("Warning: Surface %d has less than 2 Lines\n\n",(int)j);
          ret = 1;
        }
        for (l=0; l<LGM_SURFACE_NLINE(theSF); l++)
        {
          theL = LGM_SURFACE_LINE(theSF, l);
          if (theL==NULL)
          {
            /* CHECK[13] */
            if (!ret) UserWrite("\n");
            UserWriteF("Line %d is not referenced from Surface %d\n\n",(int)l,(int)j,(int)i);
            ret = 1;
            return(1);
          }
        }
        /*Laufe ueber Triangles ...*/
        for (m=0; m<LGM_SURFACE_NTRIANGLE(theSF); m++)
        {
          TheTria = LGM_SURFACE_TRIANGLE(theSF,m);
          /*UserWriteF("Triangle %d of Surface %d of Subdomain %d\n has got the following neighbour triangles: \n Neighbour_0: %d, Neighbour_1: %d, Neighbour_2: %d\n",(int)m,(int)j,(int)i,(int)(LGM_TRIANGLE_NEIGHBOR(TheTria, 0)),(int)(LGM_TRIANGLE_NEIGHBOR(TheTria, 2)),(int)(LGM_TRIANGLE_NEIGHBOR(TheTria, 2)));*/
          if (TheTria==NULL)
          {
            /* CHECK[5] */
            if (!ret) UserWrite("\n");
            UserWriteF("Triangle %d is not referenced from Surface %d\n\n",(int)m,(int)j);
            ret = 1;
            return(1);
          }

          /*laufe ueber die 3CornerIDs und pruefe, ob fuer diese CornerID ueberhaupt ein
             lokaler Point vorgesehen ist : */
          /*GOONHERE Wie gross ist eigentlich der rechte Wert?*/
          /* auskomentiert,da LGM_TRIANGLE_CORNERID die globalen IDs aus dem LGM-File liefert
                                          if(LGM_TRIANGLE_CORNERID(TheTria,0) >= LGM_SURFACE_NPOINT(theSF))
                                          {
                                                  if (!ret) UserWrite("\n");
                                                  UserWriteF("Triangle %d has the impossible CornerID %d, \n LGM_SURFACE_NPOINT-Value is %d\n\n",(int)m,(int)(LGM_TRIANGLE_CORNERID(TheTria,0)),(int)(LGM_SURFACE_NPOINT(theSF)));
                                                  ret = 1;
                                                  return(1);
                                          }
                                          if(LGM_TRIANGLE_CORNERID(TheTria,1) >= LGM_SURFACE_NPOINT(theSF))
                                          {
                                                  if (!ret) UserWrite("\n");
                                                  UserWriteF("Triangle %d has the impossible CornerID %d,\n LGM_SURFACE_NPOINT-Value is %d\n\n",(int)m,(int)(LGM_TRIANGLE_CORNERID(TheTria,1)),(int)(LGM_SURFACE_NPOINT(theSF)));
                                                  ret = 1;
                                                  return(1);
                                          }
                                          if(LGM_TRIANGLE_CORNERID(TheTria,2) >= LGM_SURFACE_NPOINT(theSF))
                                          {
                                                  if (!ret) UserWrite("\n");
                                                  UserWriteF("Triangle %d has the impossible CornerID %d,\n LGM_SURFACE_NPOINT-Value is %d\n\n",(int)m,(int)(LGM_TRIANGLE_CORNERID(TheTria,2)),(int)(LGM_SURFACE_NPOINT(theSF)));
                                                  ret = 1;
                                                  return(1);
                                          }
           */
        }
        /*Laufe ueber die Points*/
        for(p=0; p<LGM_SURFACE_NPOINT(theSF); p++)
        {
          ThePoint = LGM_SURFACE_POINT(theSF,p);
          if (ThePoint==NULL)
          {
            if (!ret) UserWrite("\n");
            UserWriteF("The %d . Point is not referenced from Surface %d\n\n",(int)p,(int)j);
            ret = 1;
            return(1);
          }
        }
      }
    }
  }
  UserWrite(" OK\n");
  UserWrite("Test B:  LeftRightInformations ofSurfaces:");


  if(CHECKB ==1)
  {
    /*CHECKB*/
    /*Checke die Left/RightInformation aller Surfaces ...*/
    /*Laufe ueber Subdomains ...*/
    for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
    {
      theSD = LGM_DOMAIN_SUBDOM(theD,i);
      /*Laufe ueber Surfaces ...*/
      for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSD); j++)
      {

        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        if ((LGM_SURFACE_LEFT(theSF) != i) && (LGM_SURFACE_RIGHT(theSF)!=i))
        {
          /* CHECK[1] */
          if (!ret) UserWrite("\n");
          UserWriteF("Surface %d does not reference Subdomain %d,\n [%d (left), %d (right)]\n\n",(int)LGM_SURFACE_ID(theSF),(int)i,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF));
          ret = 1;
          continue;
        }
        if ((LGM_SURFACE_LEFT(theSF)==i) && (LGM_SURFACE_RIGHT(theSF)==i))
        {
          /* CHECK[1] */
          if (!ret) UserWrite("\n");
          UserWriteF("Surface %d references Subdomain %d two times ,\n [%d (left), %d (right)]\n\n",(int)LGM_SURFACE_ID(theSF),(int)i,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF));
          ret = 1;
          continue;
        }
        /*UserWriteF("Surface %d does reference Subdomain %d,\n [%d (left), %d (right)]\n\n",(int)LGM_SURFACE_ID(theSF),(int)i,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF));*/

      }
    }
  }
  UserWrite(" OK\n");
  UserWrite("Test C:  Subdomains - Surfaces - Relations:");


  if(CHECKC ==1)
  {
    /*CHECKC*/
    /*laufe ueber die Subdomains und pruefe ob es stimmt dass wirklich nur die referenzierten Sbdms
       diese Surface besitzen! else {ErrorCase} */
    /*Laufe ueber Subdomains und dort ueber Surfaces ...*/
    /* CHECK[2]] */

    for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
    {
      theSD = LGM_DOMAIN_SUBDOM(theD,i);
      /*Laufe ueber Surfaces ...*/
      for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSD); j++)
      {
        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        /* CHECK[2] */
        /*laufe ueber die Subdomains und pruefe ob es stimmt dass wirklich nur die referenzierten Sbdms
           diese Surface besitzen! else {ErrorCase} */
        left = 0; right = 0;
        for (sbd=1; sbd<=LGM_DOMAIN_NSUBDOM(theD); sbd++)
        {
          lfSD = LGM_DOMAIN_SUBDOM(theD,sbd);

          for (sfce=0; sfce<LGM_SUBDOMAIN_NSURFACE(lfSD); sfce++)
          {
            lfSF = LGM_SUBDOMAIN_SURFACE(lfSD,sfce);

            if (LGM_SURFACE_ID(lfSF) == LGM_SURFACE_ID(theSF))                              /*Jede Surface wird ja nicht nur einmal abgespeichert ?!*/
            /*d.h. es genuegt nicht die Adressen zu vergleichen sondern die
               IDs muessen verglichen werden.?*/
            {
              /*eine gefunden*/
              if(LGM_SURFACE_RIGHT(theSF) == sbd)
              {
                /*UserWriteF("Surface %d is referenced by Subdomain %d  , [%d (left), %d (right)]\n\n",(int)LGM_SURFACE_ID(theSF),(int)sbd,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF));*/
                right ++;
                if(right >1)
                {
                  if (!ret) UserWrite("\n");
                  UserWriteF("Surface %d is referenced by Subdomain %d  more than one time !!, [%d (left), %d (right)]\n\n",(int)LGM_SURFACE_ID(theSF),(int)sbd,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF));
                  ret = 1;
                  continue;
                }
              }
              else if(LGM_SURFACE_LEFT(theSF) == sbd)
              {
                /*	UserWriteF("Surface %d is referenced by Subdomain %d  ,[%d (left), %d (right)]\n",(int)LGM_SURFACE_ID(theSF),(int)sbd,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF)); */
                left ++;
                if(left >1)
                {
                  if (!ret) UserWrite("\n");
                  UserWriteF("Surface %d is referenced by Subdomain %d  more than one time !!,[%d (left), %d (right)]\n",(int)LGM_SURFACE_ID(theSF),(int)sbd,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF));
                  ret = 1;
                  continue;
                }
              }
              else
              {
                if (!ret) UserWrite("\n");
                UserWriteF("Surface %d is referenced by an additional Subdomain %d, \n[%d (left), %d (right)]\n\n",(int)LGM_SURFACE_ID(theSF),(int)sbd,(int)LGM_SURFACE_LEFT(theSF),LGM_SURFACE_RIGHT(theSF));
                ret = 1;
                continue;
              }

            }
          }

        }
      }

    }
  }
  UserWrite(" OK\n");

  UserWrite("Test D:  PolylineTest:");


  if(CHECKD ==1)
  {
    /*CHECKD*/
    /*CHECK aller Lines, ob sie gar zyklisch sind und ob sie eine Vorgaenger- und eine Nachfolger-Line besitzen*/
    for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
    {
      theSD = LGM_DOMAIN_SUBDOM(theD,i);
      /*Laufe ueber Surfaces ...*/
      for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSD); j++)
      {
        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        for (l=0; l<LGM_SURFACE_NLINE(theSF); l++)
        {
          theL = LGM_SURFACE_LINE(theSF, l);
          NachfolgerLine =0;VorgaengerLine=0;
          /*zyklische Line ?*/
          if (LGM_LINE_BEGIN(theL)==LGM_LINE_END(theL))
          {
            if (!ret) UserWrite("\n");
            UserWriteF("Warning : Line %d of Surface %d of Subdomain %d is cyclic\n",(int)LGM_LINE_ID(theL),j,i);
            ret = 1;
            continue;
          }
          /* VorgaengerLine und NachfolgerLine ?*/
          for (m=0; m<LGM_SURFACE_NLINE(theSF); m++)
          {
            if (m != l)                             /*wenn nicht Line selbst...*/
            {
              theL2 = LGM_SURFACE_LINE(theSF, m);
              if (LGM_LINE_BEGIN(theL)==LGM_LINE_END(theL2))
              {
                /* NachfolgerLine gefunden*/
                VorgaengerLine++;
              }
              if (LGM_LINE_BEGIN(theL)==LGM_LINE_BEGIN(theL2))
              {
                /* NachfolgerLine gefunden*/
                VorgaengerLine++;
              }
              if (LGM_LINE_END(theL)==LGM_LINE_END(theL2))
              {
                /* VorgaengerLine gefunden*/
                NachfolgerLine++;
              }
              if (LGM_LINE_END(theL)==LGM_LINE_BEGIN(theL2))
              {
                /* VorgaengerLine gefunden*/
                NachfolgerLine++;
              }
            }
          }
          if((VorgaengerLine != 1)||(NachfolgerLine != 1))
          {
            if (!ret) UserWrite("\n");
            UserWriteF("Line %d of Surface %d of Subdomain %d \n has no or too much successor/predescessor-Lines\n\n",(int)LGM_LINE_ID(theL),j,i);
            ret = 1;
            continue;
          }
        }
      }
    }
  }
  UserWrite(" OK\n");
  UserWrite("Test E:  unused surface-Points:");


  if(CHECKE ==1)
  {
    /*CHECKE*/
    /*laufe ueber Points der Surface und suche zu jedem Point, ob es mind. ein Triangle
       gibt, das es in dieser Surface gibt, und das diesen Point besitzt.
       d.h. Pruefe ob die Surface Points besitzt, die von keinem Triangle referenziert werden.
       =>{ErrorCase}*/
    /*Laufe ueber alle Subdomains:*/

    for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
    {
      theSD = LGM_DOMAIN_SUBDOM(theD,i);
      /*Laufe ueber Surfaces ...*/
      for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSD); j++)
      {
        INT MarkKey;

        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        /*Laufe ueber Points*/

        lgmd = (LGM_DOMAIN*)aBVP;
        MarkTmpMem(lgmd->theHeap,&MarkKey);

        memsize_variable = LGM_SURFACE_NPOINT(theSF)*sizeof(LGM_POINT*);

        if ((TemporaryPointArray = (LGM_POINT **)GetTmpMem(lgmd->theHeap,memsize_variable,MarkKey))==NULL)
        {
          PrintErrorMessage('E',"BVP_Check","  ERROR: No memory for TemporaryPointArray");
          ReleaseTmpMem(lgmd->theHeap,MarkKey);
          return(1);
        }
        /* Mustermemset(nodeflag_array,0,(statistik[0]+1)*sizeof(INT)); */
        memset(TemporaryPointArray,0,memsize_variable);


        /*den allerersten Corner vor der Schleife eintargen*/
        TemporaryPointArray[0] = LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSF,0),0);
        CornerCounter = 1;
        /*laufe ueber alle Dreiecke der Surface*/
        for(m=0; m<LGM_SURFACE_NTRIANGLE(theSF); m++)
        {
          TheTria = LGM_SURFACE_TRIANGLE(theSF,m);

          /*Suche welche Corners bereits existieren*/
          Lf_Var = 0;
          NewCorners[0] = 1;
          NewCorners[1] = 1;
          NewCorners[2] = 1;
          while((Lf_Var != LGM_SURFACE_NPOINT(theSF)) &&(TemporaryPointArray[Lf_Var] != NULL))
          {
            if(TemporaryPointArray[Lf_Var] == LGM_TRIANGLE_CORNER(TheTria,0))
            {
              NewCorners[0] = 0;                                  /*gibt es bereits*/
            }
            if(TemporaryPointArray[Lf_Var] == LGM_TRIANGLE_CORNER(TheTria,1))
            {
              NewCorners[1] = 0;                                  /*gibt es bereits*/
            }
            if(TemporaryPointArray[Lf_Var] == LGM_TRIANGLE_CORNER(TheTria,2))
            {
              NewCorners[2] = 0;                                  /*gibt es bereits*/
            }
            Lf_Var++;
          }                      /*end of while*/
          if(NewCorners[0] == 1)
          {
            /*d.h. muss noch eingetragen werden*/
            if(Lf_Var == LGM_SURFACE_NPOINT(theSF))
            {
              if (!ret) UserWrite("\n");
              UserWriteF("Surface %d of Subdomain %d \n has got more different CornerIDs than available Points !!!\n\n",j,i);
              ret = 1;
              continue;
            }
            else if(CornerCounter == LGM_SURFACE_NPOINT(theSF))
            {
              if (!ret) UserWrite("\n");
              UserWriteF("CornerCounter too big - Surface %d of Subdomain %d \n has got more different CornerIDs than available Points !!!\n\n",j,i);
              ret = 1;
              continue;
            }
            else
            {
              TemporaryPointArray[CornerCounter++] = LGM_TRIANGLE_CORNER(TheTria,0);
            }
          }
          if(NewCorners[1] == 1)
          {
            /*d.h. muss noch eingetragen werden*/
            if(Lf_Var == LGM_SURFACE_NPOINT(theSF))
            {
              if (!ret) UserWrite("\n");
              UserWriteF("Surface %d of Subdomain %d \n has got more different CornerIDs than available Points !!!\n\n",j,i);
              ret = 1;
              continue;
            }
            else if(CornerCounter == LGM_SURFACE_NPOINT(theSF))
            {
              if (!ret) UserWrite("\n");
              UserWriteF("CornerCounter too big - Surface %d of Subdomain %d \n has got more different CornerIDs than available Points !!!\n\n",j,i);
              ret = 1;
              continue;
            }
            else
            {
              TemporaryPointArray[CornerCounter++] = LGM_TRIANGLE_CORNER(TheTria,1);
            }
          }
          if(NewCorners[2] == 1)
          {
            /*d.h. muss noch eingetragen werden*/
            if(Lf_Var == LGM_SURFACE_NPOINT(theSF))
            {
              if (!ret) UserWrite("\n");
              UserWriteF("Surface %d of Subdomain %d \n has got more different CornerIDs than available Points !!!\n\n",j,i);
              ret = 1;
              continue;
            }
            else if(CornerCounter == LGM_SURFACE_NPOINT(theSF))
            {
              if (!ret) UserWrite("\n");
              UserWriteF("CornerCounter too big - Surface %d of Subdomain %d \n has got more different CornerIDs than available Points !!!\n\n",j,i);
              ret = 1;
              continue;
            }
            else
            {
              TemporaryPointArray[CornerCounter++] = LGM_TRIANGLE_CORNER(TheTria,2);
            }
          }
        }

        /*Nun muessen alle CornerIDs gefunden worden sein:*/
        if(CornerCounter < LGM_SURFACE_NPOINT(theSF))
        {
          if (!ret) UserWrite("\n");
          UserWriteF("Surface %d of Subdomain %d \n has got unused Points !!!\n\n",j,i);
          ret = 1;
          continue;
        }
        /*Hat ein Point kein Flag erhalten und ist somit in keinem einzigen Dreieck benutzt?*/
        for(p=0; p<LGM_SURFACE_NPOINT(theSF); p++)
        {
          if(TemporaryPointArray[p] == NULL)
          {
            if (!ret) UserWrite("\n");
            UserWriteF("The %d . Point of Surface %d of Subdomain %d \n ist not (!) used in any Triangle of that Surface !!!\n\n",p,j,i);
            ret = 1;
            continue;
          }
        }

        ReleaseTmpMem(lgmd->theHeap,MarkKey);

      }
    }
  }
  UserWrite(" OK\n");


  /**********************/
  /*Pruefe die  Dreiecke*/
  /**********************/
  UserWrite("Test F: TriangleNeighbourhood and -orientation:");
  if(CHECKF ==1)
  {
    /*CHECKF*/
    /*Gebe die 3 NachbarTriangles an und pruefe ob es ueberhaupt ein Nachbartriangle ist und ob
       die selbe Direction vorliegt! d.h dass die beiden gemeinsamen IDs in der umgekehrten
       Reihenfolge durchlaufen werden. und pruefe auch ob 3 unterschiedliche CornerIDs vorliegen*/
    for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
    {
      theSD = LGM_DOMAIN_SUBDOM(theD,i);
      /*Laufe ueber Surfaces ...*/
      for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSD); j++)
      {
        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        /*Laufe ueber die Triangles*/
        for(m=0; m<LGM_SURFACE_NTRIANGLE(theSF); m++)
        {
          TheTria = LGM_SURFACE_TRIANGLE(theSF,m);
          /*Hat Tria 3 unterschiedliche Corner IDs?*/
          if((LGM_TRIANGLE_CORNER(TheTria,0) == LGM_TRIANGLE_CORNER(TheTria,1)) ||
             (LGM_TRIANGLE_CORNER(TheTria,1) == LGM_TRIANGLE_CORNER(TheTria,2)) ||
             (LGM_TRIANGLE_CORNER(TheTria,0) == LGM_TRIANGLE_CORNER(TheTria,2))   )
          {
            if (!ret) UserWrite("\n");
            UserWriteF("The %d . Triangle of Surface %d of Subdomain %d \n doesnt have 3 different CornerIDs !!!\n\n",m,j,i);
            ret = 1;
            continue;
          }

          /*Laufe ueber die Nachbarn*/
          for(p=0; p<3; p++)
          {
            /*wenn ueberhaupt ein Nachbar ...*/
            if(LGM_TRIANGLE_NEIGHBOR(TheTria, p) != -1)
            {
              TheNgbTria = LGM_SURFACE_TRIANGLE(theSF, (LGM_TRIANGLE_NEIGHBOR(TheTria, p)));
              /*Hat TheNgbTria 3 unterschiedliche Corner IDs?*/
              if((LGM_TRIANGLE_CORNER(TheNgbTria,0) == LGM_TRIANGLE_CORNER(TheNgbTria,1)) ||
                 (LGM_TRIANGLE_CORNER(TheNgbTria,1) == LGM_TRIANGLE_CORNER(TheNgbTria,2)) ||
                 (LGM_TRIANGLE_CORNER(TheNgbTria,0) == LGM_TRIANGLE_CORNER(TheNgbTria,2))   )
              {
                if (!ret) UserWrite("\n");
                UserWriteF("The %d .Triangle of Surface %d of Subdomain %d\n  doesnt have 3 different CornerIDs !!!\n\n",m,j,i);
                ret = 1;
                continue;
              }
              /*laufe ueber CornerIDs von Tria und suche die beiden gemeinsamen IDS */
              /*laufe ueber CornerIDs von TheTria*/
              commmonIDs_Counter = -1;
              for(tr_index = 0; tr_index<3; tr_index++)
              {
                /*laufe ueber CornerIDs von TheNgbTria*/
                for(ntr_index = 0; ntr_index<3; ntr_index++)
                {
                  /*Sind die CornerIDs identisch?*/
                  if(LGM_TRIANGLE_CORNER(TheTria,tr_index) == LGM_TRIANGLE_CORNER(TheNgbTria,ntr_index))
                  {
                    commmonIDs_Counter++;
                    if(commmonIDs_Counter ==2)
                    {
                      if (!ret) UserWrite("\n");
                      UserWriteF("The %d . Triangle and its negibour the %d . Triangle \n of Surface %d of Subdomain %d \n do have the same CornerIDs !!!\n\n",m,LGM_TRIANGLE_NEIGHBOR(TheTria, p),j,i);
                      ret = 1;
                      continue;
                    }
                    else
                    {
                      if(tr_index == p)
                      {
                        if (!ret) UserWrite("\n");
                        UserWriteF("The third corner of the %d. Triangle is also found in its negibour the %d . Triangle \n of Surface %d of Subdomain %d \n \n",m,LGM_TRIANGLE_NEIGHBOR(TheTria, p),j,i);
                        ret = 1;
                        continue;
                      }
                      /*die beiden Indizes des gemeinsamen COrners merken:*/
                      indexesTheTria[commmonIDs_Counter] = tr_index;
                      indexesTheNgbTria[commmonIDs_Counter] = ntr_index;
                    }
                  }
                }
              }
              if(commmonIDs_Counter != 1)
              {
                if (!ret) UserWrite("\n");
                UserWriteF("The %d . Triangle and its negibour the %d . Triangle \n of Surface %d of Subdomain %d \n have less than 2 common CornerIDs !!!\n\n",m,LGM_TRIANGLE_NEIGHBOR(TheTria, p),j,i);
                ret = 1;
                continue;
              }
              else
              {
                /*Pruefung der Orientierung der beiden benachbarten Dreiecke*/
                /* Die beiden gemeinsamen Corneer muessen in unterschiedlicher Reihenfolge vorkommen*/

                /*Reihenfolge der gemeinsamen CornerIDs bei TheTria*/
                if((indexesTheTria[0] == 0) && (indexesTheTria[1] == 2))
                {
                  richtung_TheTria = -1;
                }
                else
                {
                  richtung_TheTria = 1;
                }
                /*Reihenfolge der gemeinsamen CornerIDs bei TheNgbTria*/
                if( ((indexesTheNgbTria[0] == 0) && (indexesTheNgbTria[1] == 2)) ||
                    ((indexesTheNgbTria[0] == 1) && (indexesTheNgbTria[1] == 0)) ||
                    ((indexesTheNgbTria[0] == 2) && (indexesTheNgbTria[1] == 1))
                    )
                {
                  richtung_TheNgbTria = -1;
                }
                else
                {
                  richtung_TheNgbTria = 1;
                }

                /*Sind die beiden Durchlaufrichtungen falscherwiese identisch ?*/
                if(richtung_TheTria==richtung_TheNgbTria)
                {
                  if (!ret) UserWrite("\n");
                  UserWriteF("TriangleOrientationError: The 2 common Corners of the \n  %d . Triangle and its negibour the \n %d . Triangle of Surface %d of Subdomain %d \n have the same direction!!!\n\n",m,LGM_TRIANGLE_NEIGHBOR(TheTria, p),j,i);
                  ret = 1;
                  continue;
                }
              }
            }
          }
        }
      }
    }
  }
  UserWrite(" OK if ANSYS-INPUT \n");

  UserWrite("Test G: double triangles:");

  if(CHECKG ==1)
  {
    /*CHECKG*/
    /*Pruefe zu jedem Triangle, ob es nur einmal vorkommt*/
    for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
    {
      theSD = LGM_DOMAIN_SUBDOM(theD,i);
      /*Laufe ueber Surfaces ...*/
      for (j=0; j<LGM_SUBDOMAIN_NSURFACE(theSD); j++)
      {
        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        /*Laufe ueber die Triangles*/
        for(m=0; m<LGM_SURFACE_NTRIANGLE(theSF); m++)
        {
          TheTria = LGM_SURFACE_TRIANGLE(theSF,m);

          for (i2=1; i2<=LGM_DOMAIN_NSUBDOM(theD); i2++)
          {
            theSD2 = LGM_DOMAIN_SUBDOM(theD,i2);
            /*Laufe ueber Surfaces ...*/
            for (j2=0; j2<LGM_SUBDOMAIN_NSURFACE(theSD2); j2++)
            {
              theSF2 = LGM_SUBDOMAIN_SURFACE(theSD2,j2);
              /*Laufe ueber die Triangles*/
              for(m2=0; m2<LGM_SURFACE_NTRIANGLE(theSF2); m2++)
              {
                TheTria2 = LGM_SURFACE_TRIANGLE(theSF2,m2);
                if(TheTria != TheTria2)
                {
                  yes[0] = 0;
                  yes[1] = 0;
                  yes[2] = 0;
                  /*laufe ueber die 3 Corners von TheTria*/
                  for(tr_index=0; tr_index<3; tr_index++)
                  {
                    Point_TheTria = LGM_TRIANGLE_CORNER(TheTria,tr_index);
                    /*laufe ueber die 3 Corners von TheNgbTria*/
                    for(ntr_index=0; ntr_index<3; ntr_index++)
                    {
                      Point_TheNgbTria = LGM_TRIANGLE_CORNER(TheTria2,ntr_index);
                      if(LGM_POINT_DIST(Point_TheTria,Point_TheNgbTria) < 1E-6)
                      {
                        yes[tr_index] = 1;
                      }
                    }
                  }
                  /*UserWriteF("The %d .Triangle of Surface %d of Subdomain %d and The \n %d . Triangle of Surface %d of Subdomain %d \n have the following yes-array: [0]: %d, [1]: %d, [2]: %d\n\n",m,j,i,m2,j2,i2,yes[0],yes[1],yes[2]); */
                  if((yes[0] == 1)&&(yes[1] == 1)&&(yes[2] == 1))
                  {
                    if (!ret) UserWrite("\n");
                    UserWriteF("The %d .Triangle of Surface %d of Subdomain %d and The \n %d . Triangle of Surface %d of Subdomain %d \n have identical corners !!!!!\n\n",m,j,i,m2,j2,i2);
                    ret = 1;
                    continue;
                  }
                }
              }
            }
          }

        }
      }
    }
  }

  UserWrite(" OK\n\n");
  return (ret);
}

static INT Normal_Vector(DOUBLE *p0, DOUBLE *p1, DOUBLE *p2, DOUBLE *n)
{
  DOUBLE l1, l2, l, n1[3], n2[3];

  n1[0] = p2[0] - p0[0];
  n1[1] = p2[1] - p0[1];
  n1[2] = p2[2] - p0[2];

  V3_EUKLIDNORM(n1,l1);
  /*    Scale(n1, 1.0/l1); */
  V3_SCALE(1.0/l1,n1);

  n2[0] = p2[0] - p1[0];
  n2[1] = p2[1] - p1[1];
  n2[2] = p2[2] - p1[2];

  V3_EUKLIDNORM(n2,l2);
  /*    Scale(n2, 1.0/l2); */
  V3_SCALE(1.0/l2, n2);

  /*    Cross(n, n1, n2); */
  V3_VECTOR_PRODUCT(n1,n2,n);
  V3_EUKLIDNORM(n,l);
  /*    Scale(n, 1.0/l); */
  V3_SCALE(1.0/l, n);

  return(0);
}

static DOUBLE Calc_Triangle_Angle(LGM_SURFACE *theSurface, INT i,  INT j)
{
  DOUBLE n1[3], n2[3], *p0, *p1, *p2, scalarproduct;

  p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
  p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
  p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

  /* get Normalvector to the Triangle i */
  Normal_Vector(p0, p1, p2, &n1[0]);

  p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,j),0);
  p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,j),1);
  p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,j),2);

  /* get Normalvector to the Triangle j */
  Normal_Vector(p0, p1, p2, &n2[0]);

  /*    scalarproduct =  Mult(n1, n2); */
  V3_SCALAR_PRODUCT(n1, n2, scalarproduct);
  return(scalarproduct);
}

static INT Count_Folds_Surface(LGM_SURFACE *theSurface)
{
  INT i, j, folds;
  DOUBLE scalarproduct;

  folds = 0;

  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    for(j=0; j<3; j++)
    {
      if(LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j)!=-1)
      {
        scalarproduct = Calc_Triangle_Angle(theSurface, i, LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j));

        if(scalarproduct < cosAngle)
          folds++;
      }
    }
  }
  return(folds);
}

static INT POINT_DIST(LGM_POINT *p0, LGM_POINT *p1)
{
  DOUBLE x, y, z;

  x = LGM_POINT_POS(p0)[0] - LGM_POINT_POS(p1)[0];
  y = LGM_POINT_POS(p0)[1] - LGM_POINT_POS(p1)[1];
  z = LGM_POINT_POS(p0)[2] - LGM_POINT_POS(p1)[2];
  if( (x<SMALL) && (x>-SMALL) &&
      (y<SMALL) && (y>-SMALL) &&
      (z<SMALL) && (z>-SMALL) )
    return(1);
  else
    return(0);
}

static INT AddPoint2List(LGM_POINT *Point, LGM_POINT *pointlist, INT *norp)
{
  int i, flag;
  LGM_POINT p;

  flag = 0;
  for(i=0; i<(*norp); i++)
  {
    p.position[0] = pointlist[i].position[0];
    p.position[1] = pointlist[i].position[1];
    p.position[2] = pointlist[i].position[2];
    if(POINT_DIST(Point, &p)==1)
      flag = 1;
  }
  if(flag==0)
  {
    pointlist[(*norp)].position[0] = LGM_POINT_POS(Point)[0];
    pointlist[(*norp)].position[1] = LGM_POINT_POS(Point)[1];
    pointlist[(*norp)].position[2] = LGM_POINT_POS(Point)[2];
    (*norp)++;
  }
  return (0);
}

static INT ResolvePoints(HEAP *Heap, LGM_DOMAIN *theDomain, LGM_POINT *pointlist, INT *norp)
{
  INT i, j;
  DOUBLE scalarproduct;
  LGM_SURFACE *theSurface;

  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
    {
      for(j=0; j<3; j++)
      {
        if(LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j)!=-1)
        {
          scalarproduct = Calc_Triangle_Angle(theSurface, i, LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j));

          if(scalarproduct < cosAngle)
          {
            AddPoint2List(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3), pointlist, norp);
            AddPoint2List(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+2)%3), pointlist, norp);
          }
        }
      }
    }

  }
  return(0);
}

static INT Count_Points(LGM_DOMAIN *theDomain)
{
  INT i, j, npoints;
  DOUBLE scalarproduct;
  LGM_SURFACE *theSurface;

  npoints = 0;
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
    {
      for(j=0; j<3; j++)
      {
        if(LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j)!=-1)
        {
          scalarproduct = Calc_Triangle_Angle(theSurface, i, LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j));

          if(scalarproduct < cosAngle)
            npoints = npoints + 2;
        }
      }
    }

  }
  return(npoints);
}

static INT Write_Line(LGM_LINE *theLine)
{
  char buff[5], name[12];
  INT i, id;
  FILE *stream;
  DOUBLE global[3];
  LINEPOINT *help;

  id = LGM_LINE_ID(theLine);
  name[0] = 'l';
  name[1] = 'i';
  name[2] = 'n';
  name[3] = 'e';
  sprintf(buff,"%d",(int)id);
  name[4] = buff[0];
  name[5] = buff[1];
  name[6] = buff[2];
  name[7] = buff[3];
  name[8] = buff[4];

  stream = fopen(name,"w+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }

  fprintf(stream, "%d\n", (int)LGM_LINE_NPOINT(theLine));
  for(i=0; i<LGM_LINE_NPOINT(theLine); i++)
    fprintf(stream, "%f %f %f\n",LGM_LINE_POINT(theLine,i)->position[0],
            LGM_LINE_POINT(theLine,i)->position[1],
            LGM_LINE_POINT(theLine,i)->position[2]);

  fprintf(stream, "%d\n", (int)LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine)));
  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
  for(i=0; i<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine)); i++)
  {
    Line_Local2GlobalNew(theLine,global,help->local);
    fprintf(stream, "%f %f %f\n",global[0],global[1],global[2]);
    help = help->next;
  }

  fclose(stream);

  return(0);
}


static INT Write_Surface(LGM_SURFACE *theSurface, const char *name, const char *name1)
{
  INT i, id;
  FILE *stream;
  DOUBLE local[3], global[3];

  id = LGM_SURFACE_ID(theSurface);

  stream = fopen(name,"w+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }
  fprintf(stream, "%d\n", (int)LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)));
  for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); i++)
  {
    local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,0);
    local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,1);
    fprintf(stream, "%20.16f %20.16f\n", local[0], local[1]);
  }
  fprintf(stream, "%d\n", (int)LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)));
  for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)); i++)
  {
    fprintf(stream, "%d %d %d\n",   (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,0),
            (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,1),
            (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,2));
  }
  if(SURFACE_ADAPT)
    if(LGM_SURFACE_ID(theSurface)!=-1)
    {
      for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(xy_Surface)); i++)
        fprintf(stream,"%d\n",(int)xy_new_id[i]);
      for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(xy_Surface)); i++)
        fprintf(stream,"%d\n",(int)xy_element_id[i]);
    }
  fclose(stream);

  stream = fopen(name1,"w+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }
  fprintf(stream, "%s\n", "surfacemesh");
  fprintf(stream, "%d\n", (int)LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)));
  for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); i++)
  {
    local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,0);
    local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,1);
    Surface_Local2Global(theSurface, global, local);
    fprintf(stream, "%f %f %f\n", global[0], global[1], global[2]);
  }
  fprintf(stream, "%d\n", (int)LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)));
  for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)); i++)
  {
    fprintf(stream, "%d %d %d\n",   (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,0),
            (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,1),
            (int)LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,2));
  }
  fclose(stream);

  return(0);
}

static DOUBLE Calc_Angle_to_Normal(LGM_SURFACE *theSurface, INT i, DOUBLE *n)
{
  DOUBLE n1[3], *p0, *p1, *p2, scalarproduct;

  p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
  p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
  p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

  /* get Normalvector to the Triangle i */
  Normal_Vector(p0, p1, p2, &n1[0]);

  /*    scalarproduct =  Mult(n1, n); */
  V3_SCALAR_PRODUCT(n1, n, scalarproduct);
  return(scalarproduct);
}

static INT Surface_Plane(LGM_SURFACE *theSurface)
{
  INT i, plain1, plain2;
  DOUBLE sp;

  plain1 = 1;
  plain2 = 1;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    sp = Calc_Angle_to_Normal(theSurface, i, xy_n);
    if(sp<0.995)
      plain1 = 0;
    if(sp>-0.995)
      plain2 = 0;
  }
  if(plain1||plain2)
    return(1);
  else
    return(0);
}


static INT Read_Line(HEAP *Heap, LGM_LINE *theLine, INT MarkKey, char name[12])
{
  int d;
  INT i;
  FILE *stream;
  DOUBLE global[3];
  LINEPOINT *help;

  stream = fopen(name,"r+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }

  fscanf(stream, "%d\n",&d);
  LGM_LINE_NPOINT(theLine) = d;
  for(i=0; i<LGM_LINE_NPOINT(theLine); i++)
    fscanf(stream, "%lf %lf %lf\n",&LGM_LINE_POINT(theLine,i)->position[0],
           &LGM_LINE_POINT(theLine,i)->position[1],
           &LGM_LINE_POINT(theLine,i)->position[2]);

  fscanf(stream, "%d\n",&d);
  LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine)) = d;
  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
  for(i=0; i<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine)); i++)
  {
    help = (LINEPOINT*)GetTmpMem(Heap,sizeof(LINEPOINT), MarkKey);
    fscanf(stream, "%lf %lf %lf\n",&global[0],&global[1],&global[2]);
    Line_Local2GlobalNew(theLine,global,help->local);
    help = help->next;
  }

  fclose(stream);

  return(0);
}

static INT Read_Surface(HEAP *Heap, LGM_SURFACE *theSurface, INT MarkKey, char name[12])
{
  INT i;
  int d, d1, d2, d3;
  FILE *stream;
  double local0, local1;

  stream = fopen(name,"r+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }

  fscanf(stream,"%d\n",&d);
  LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)) = d;

  LGM_SURFACE_DISC(theSurface)->local = (DOUBLE **) GetTmpMem(Heap,(d+1)*sizeof(DOUBLE*), MarkKey);
  if(LGM_SURFACE_DISC(theSurface)->local==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<d; i++)
  {
    LGM_SURFACE_DISC(theSurface)->local[i] = (DOUBLE *)  GetTmpMem(Heap,3*sizeof(DOUBLE), MarkKey);
    if(LGM_SURFACE_DISC(theSurface)->local[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
    fscanf(stream,"%lg %lg\n",&local0, &local1);
    LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,0) = local0;
    LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,1) = local1;
  }

  fscanf(stream,"%d",&d);
  LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)) = d;

  LGM_SURFACE_DISC(theSurface)->triangle = (INT **) GetTmpMem(Heap,(d+1)*sizeof(INT*), MarkKey);
  if(LGM_SURFACE_DISC(theSurface)->triangle==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<d; i++)
  {
    LGM_SURFACE_DISC(theSurface)->triangle[i] = (INT *)  GetTmpMem(Heap,4*sizeof(INT), MarkKey);
    if(LGM_SURFACE_DISC(theSurface)->triangle[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
    fscanf(stream,"%d %d %d\n",&d1, &d2, &d3);
    LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,0) = d1;
    LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,1) = d2;
    LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,2) = d3;
  }
  if(SURFACE_ADAPT)
    if(LGM_SURFACE_ID(theSurface)!=-1)
    {
      for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(xy_Surface)); i++)
      {
        fscanf(stream,"%d\n",&d1);
        xy_new_id[i] = d1;
      }
      for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(xy_Surface)); i++)
      {
        fscanf(stream,"%d\n",&d1);
        xy_element_id[i] = d1;
      }
    }
  fclose(stream);

  return(0);
}

static INT Generate_Basis_Mesh(LGM_DOMAIN *theDomain, DOUBLE h, HEAP *Heap, INT D, INT MarkKey)
{
  INT i, j, size;
  DOUBLE min[3], max[3];
  char name[11], name1[15];
  LGM_LINE *theLine;
  FILE *stream;

  min[0]=min[1]=min[2]=MAX_C;
  max[0]=max[1]=max[2]=-MAX_C;

  xy_n[0] = 0.0;
  xy_n[1] = 0.0;
  xy_n[2] = 1.0;

  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
    for (i=0; i<LGM_LINE_NPOINT(theLine); i++)
    {
      min[0] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[0],min[0]);
      min[1] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[1],min[1]);
      min[2] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[2],min[2]);
      max[0] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[0],max[0]);
      max[1] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[1],max[1]);
      max[2] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[2],max[2]);
    }

  for (i=0; i<4; i++)
    if ((xy_Point[i] = (LGM_POINT*)GetTmpMem(Heap,sizeof(LGM_POINT), MarkKey)) == NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }

  LGM_POINT_POS(xy_Point[0])[0] = min[0] - 2 * h;
  LGM_POINT_POS(xy_Point[0])[1] = min[1] - 2 * h;
  LGM_POINT_POS(xy_Point[0])[2] = 0.0;

  LGM_POINT_POS(xy_Point[1])[0] = max[0] + 2 * h;
  LGM_POINT_POS(xy_Point[1])[1] = min[1] - 2 * h;
  LGM_POINT_POS(xy_Point[1])[2] = 0.0;

  LGM_POINT_POS(xy_Point[2])[0] = max[0] + 2 * h;
  LGM_POINT_POS(xy_Point[2])[1] = max[1] + 2 * h;
  LGM_POINT_POS(xy_Point[2])[2] = 0.0;

  LGM_POINT_POS(xy_Point[3])[0] = min[0] - 2 * h;
  LGM_POINT_POS(xy_Point[3])[1] = max[1] + 2 * h;
  LGM_POINT_POS(xy_Point[3])[2] = 0.0;

  for (i=0; i<4; i++)
  {
    size = sizeof(LGM_LINE) + (2)*sizeof(LGM_POINT);
    if ((xy_Line[i] = (LGM_LINE*)GetTmpMem(Heap,size, MarkKey)) == NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }

    LGM_LINE_NPOINT(xy_Line[i]) = 2;
    LGM_LINE_ID(xy_Line[i])  = i;
  }

  LGM_LINE_BEGIN(xy_Line[0]) = 0;
  LGM_LINE_BEGIN(xy_Line[1]) = 1;
  LGM_LINE_BEGIN(xy_Line[2]) = 2;
  LGM_LINE_BEGIN(xy_Line[3]) = 3;
  LGM_LINE_END(xy_Line[0]) = 1;
  LGM_LINE_END(xy_Line[1]) = 2;
  LGM_LINE_END(xy_Line[2]) = 3;
  LGM_LINE_END(xy_Line[3]) = 0;

  for (i=0; i<4; i++)
  {
    for(j=0; j<3; j++)
    {
      LGM_POINT_POS(LGM_LINE_POINT(xy_Line[i], 0))[j] = LGM_POINT_POS(xy_Point[LGM_LINE_BEGIN(xy_Line[i])])[j];
      LGM_POINT_POS(LGM_LINE_POINT(xy_Line[i], 1))[j] = LGM_POINT_POS(xy_Point[LGM_LINE_END(xy_Line[i])])[j];
    }
  }
  for (i=0; i<4; i++)
  {
    xy_Line[i]->ldiscnew = (LGM_LINEDISCNEW*)GetTmpMem(Heap,sizeof(LGM_LINEDISCNEW)*1, MarkKey);
    if( (xy_Line[i]->ldiscnew==NULL) )
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
    if (DiscretizeLineNew(Heap,xy_Line[i],h,NULL, 0, MarkKey))
      return(1);
    if(LGM_DEBUG) PrintLineInfo(xy_Line[i]);
  }
  xy_bpoints =    LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(xy_Line[0]))
               +       LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(xy_Line[1]))
               +       LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(xy_Line[2]))
               +       LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(xy_Line[3])) - 4;

  size = sizeof(LGM_SURFACE)+(4)*sizeof(LGM_POINT);
  if ((xy_Surface=(LGM_SURFACE*)GetTmpMem(Heap,size, MarkKey)) == NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  if ((LGM_SURFACE_FPOINT(xy_Surface) = (LGM_POINT*)GetTmpMem(Heap,sizeof(LGM_POINT)*4, MarkKey)) == NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  if ((LGM_SURFACE_FTRIANGLE(xy_Surface) = (LGM_TRIANGLE*)GetTmpMem(Heap,sizeof(LGM_TRIANGLE)*2, MarkKey)) == NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  LGM_SURFACE_NPOINT(xy_Surface)          = 4;
  LGM_SURFACE_ID(xy_Surface)                      = -1;
  LGM_SURFACE_NTRIANGLE(xy_Surface)       = 2;
  LGM_SURFACE_LEFT(xy_Surface)            = 0;
  LGM_SURFACE_RIGHT(xy_Surface)           = 0;
  LGM_SURFACE_NLINE(xy_Surface)           = 4;
  LGM_SURFACE_DISC(xy_Surface)            = NULL;
  LGM_SURFACE_BNDCOND(xy_Surface)         = NULL;

  for (j=0; j<LGM_SURFACE_NLINE(xy_Surface); j++)
    LGM_SURFACE_LINE(xy_Surface,j) = xy_Line[j];

  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,0))[0] = min[0] - 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,0))[1] = min[1] - 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,0))[2] = 0.0;

  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,1))[0] = max[0] + 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,1))[1] = min[1] - 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,1))[2] = 0.0;

  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,2))[0] = max[0] + 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,2))[1] = max[1] + 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,2))[2] = 0.0;

  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,3))[0] = min[0] - 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,3))[1] = max[1] + 2 * h;
  LGM_POINT_POS(LGM_SURFACE_POINT(xy_Surface,3))[2] = 0.0;

  LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(xy_Surface,0),0) = LGM_SURFACE_POINT(xy_Surface,0);
  LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(xy_Surface,0),1) = LGM_SURFACE_POINT(xy_Surface,1);
  LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(xy_Surface,0),2) = LGM_SURFACE_POINT(xy_Surface,2);
  LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(xy_Surface,1),0) = LGM_SURFACE_POINT(xy_Surface,0);
  LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(xy_Surface,1),1) = LGM_SURFACE_POINT(xy_Surface,2);
  LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(xy_Surface,1),2) = LGM_SURFACE_POINT(xy_Surface,3);

  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(xy_Surface,0),0) = 0;
  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(xy_Surface,0),1) = 1;
  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(xy_Surface,0),2) = 2;
  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(xy_Surface,1),0) = 0;
  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(xy_Surface,1),1) = 2;
  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(xy_Surface,1),2) = 3;

  for(i=0; i<2; i++)
    for(j=0; j<3; j++)
      LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(xy_Surface,i),j) = -1;

  LGM_SURFACE_DISC(xy_Surface) = (LGM_SURFDISC*)GetTmpMem(Heap,sizeof(LGM_SURFDISC)*1,MarkKey);
  if( (LGM_SURFACE_DISC(xy_Surface)==NULL) )
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  if(SAVE_SURFACE)
  {
    name[0] = 'x';
    name[1] = 'y';
    name[2] = '_';
    name[3] = 's';
    name[4] = 'u';
    name[5] = 'r';
    name[6] = 'f';
    name[7] = 'a';
    name[8] = 'c';
    name[9] = 'e';
    name[10] = '\0';

    name1[0] = 'x';
    name1[1] = 'y';
    name1[2] = '_';
    name1[3] = 's';
    name1[4] = 'u';
    name1[5] = 'r';
    name1[6] = 'f';
    name1[7] = 'a';
    name1[8] = 'c';
    name1[9] = 'e';
    name1[10] = 'm';
    name1[11] = 'e';
    name1[12] = 's';
    name1[13] = 'h';
    name1[14] = '\0';

    stream = fopen(name,"r+");
    if (stream==NULL)
    {
      /* fclose(stream); */
      printf("%s\n", "xy_Surface not triangulated, do now");

      if (DiscretizeSurface(Heap,xy_Surface,NULL,h,NULL,0,D,MarkKey))
        return(1);
      Write_Surface(xy_Surface, name, name1);
    }
    else
    {
      fclose(stream);
      printf("%s %d\n", "Read Surface ", (int)LGM_SURFACE_ID(xy_Surface));
      Read_Surface(Heap, xy_Surface, MarkKey, name);
    }
  }
  else
  {
    if (DiscretizeSurface(Heap,xy_Surface,NULL,h,NULL,0,D,MarkKey))
      return(1);
  }

  size = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(xy_Surface)) + 1;
  xy_lcoord = (DOUBLE**)GetTmpMem(Heap,sizeof(DOUBLE*)*size, MarkKey);
  if(xy_lcoord==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  for(i=0; i<size; i++)
  {
    xy_lcoord[i] = (DOUBLE*)GetTmpMem(Heap,sizeof(DOUBLE)*2, MarkKey);
    if(xy_lcoord[i]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
  }
  xy_flag = (INT*)GetTmpMem(Heap,sizeof(INT)*size, MarkKey);
  if(xy_flag==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  xy_new_id = (INT*)GetTmpMem(Heap,sizeof(INT)*size, MarkKey);
  if(xy_new_id==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  xy_element_id = (INT*)GetTmpMem(Heap,sizeof(INT)*(LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(xy_Surface))+1), MarkKey);
  if(xy_element_id==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  return(0);
}

#define MAX_T 30

INT Aux_Check(HEAP *Heap, LGM_DOMAIN *theDomain, INT MarkKey)
{
  INT i, j, k, l, ni, sid, ntriangle, npoint, **point_list, **nb1, **nb2, corner_id, max_tr;
  INT a, b, c, d, flag1, flag2, return_flag1, return_flag2, **triangle_list, ntr;
  LGM_SURFACE *theSurface;

  return_flag1 = return_flag2 = 0;
  npoint = LGM_DOMAIN_NPOINT(theDomain);
  point_list = (INT **) GetTmpMem(Heap,(npoint+1)*sizeof(INT*), LGM_MarkKey);
  if(point_list==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<npoint; i++)
  {
    point_list[i] = (INT *)  GetTmpMem(Heap,MAX_T*sizeof(INT),LGM_MarkKey);
    if(point_list[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }

  max_tr = 0;
  for(sid=1; sid<=LGM_DOMAIN_NSUBDOM(theDomain); sid++)
  {
    ntriangle = 0;
    for(theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
      if( (LGM_SURFACE_LEFT(theSurface)==sid)||(LGM_SURFACE_RIGHT(theSurface)==sid) )
        ntriangle = ntriangle + LGM_SURFACE_NTRIANGLE(theSurface);
    if(max_tr<ntriangle)
      max_tr = ntriangle;
  }

  nb1 = (INT **) GetTmpMem(Heap,(max_tr+1)*sizeof(INT*), LGM_MarkKey);
  if(nb1==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<max_tr; i++)
  {
    nb1[i] = (INT *)  GetTmpMem(Heap,3*sizeof(INT),LGM_MarkKey);
    if(nb1[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }
  nb2 = (INT **) GetTmpMem(Heap,(max_tr+1)*sizeof(INT*), LGM_MarkKey);
  if(nb2==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<max_tr; i++)
  {
    nb2[i] = (INT *)  GetTmpMem(Heap,3*sizeof(INT),LGM_MarkKey);
    if(nb2[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }
  triangle_list = (INT **) GetTmpMem(Heap,(max_tr+1)*sizeof(INT*), LGM_MarkKey);
  if(triangle_list==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<max_tr; i++)
  {
    triangle_list[i] = (INT *)  GetTmpMem(Heap,3*sizeof(INT),LGM_MarkKey);
    if(triangle_list[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }

  for(sid=1; sid<=LGM_DOMAIN_NSUBDOM(theDomain); sid++)
  {
    ntriangle = 0;
    ntr = 0;
    for(theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
      if( (LGM_SURFACE_LEFT(theSurface)==sid)||(LGM_SURFACE_RIGHT(theSurface)==sid) )
      {
        ntriangle = ntriangle + LGM_SURFACE_NTRIANGLE(theSurface);
        for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
        {
          if(LGM_SURFACE_LEFT(theSurface)==sid)
          {
            triangle_list[ntr][0] = LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface, i), 0);
            triangle_list[ntr][1] = LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface, i), 1);
            triangle_list[ntr][2] = LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface, i), 2);
          }
          if(LGM_SURFACE_RIGHT(theSurface)==sid)
          {
            triangle_list[ntr][0] = LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface, i), 1);
            triangle_list[ntr][1] = LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface, i), 0);
            triangle_list[ntr][2] = LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface, i), 2);
          }
          ntr++;
        }
      }
    assert(ntr=ntriangle);
    for(i=0; i<ntriangle; i++)
      for(j=0; j<3; j++)
      {
        nb1[i][j] = -1;
        nb2[i][j] = -1;
      }

    for(i=0; i<npoint; i++)
    {
      point_list[i][0] = 0;
      for(j=1; j<MAX_T; j++)
        point_list[i][j] = -1;
    }

    for(i=0; i<ntriangle; i++)
      for(j=0; j<3; j++)
      {
        corner_id = triangle_list[i][j];
        point_list[corner_id][++point_list[corner_id][0]] = i;
      }

    for(ni=0; ni<npoint; ni++)
      for(i=1; i<=point_list[ni][0]; i++)
        for(j=1; j<=point_list[ni][0]; j++)
          if(i!=j)
            for(k=0; k<3; k++)
              for(l=0; l<3; l++)
              {
                a = triangle_list[point_list[ni][i]][(k+1)%3];
                b = triangle_list[point_list[ni][i]][(k+2)%3];
                c = triangle_list[point_list[ni][j]][(l+2)%3];
                d = triangle_list[point_list[ni][j]][(l+1)%3];
                if( ((a==c)&&(b==d)) )
                  nb1[point_list[ni][i]][k] = point_list[ni][j];
                if( ((a==c)&&(b==d))||((a==d)&&(b==c)) )
                  nb2[point_list[ni][i]][k] = point_list[ni][j];

              }

    flag1 = flag2 = 0;
    for(i=0; i<ntriangle; i++)
      for(j=0; j<3; j++)
      {
        if(nb1[i][j]==-1)
          flag1 = 1;
        if(nb2[i][j]==-1)
          flag2 = 1;
      }

    if(flag2)
    {
      UserWriteF("%s %d %s\n","Subdomain", sid, "not closed");
      UserWriteF("%s\n","Check the geometry file");
      return_flag2 = 1;
    }
    if(!flag2&&flag1)
    {
      UserWriteF("%s %d %s\n","Orientation of surfaces in subdoamin", sid, "is wrong");
      UserWriteF("%s\n","Check the geometry file");
      return_flag1 = 1;
    }

  }

  if(return_flag1||return_flag2)
    return(1);

  return(0);
}

/* domain interface function: for description see domain.h */
MESH* NS_DIM_PREFIX BVP_GenerateMesh (HEAP *Heap, BVP *aBVP, INT argc, char **argv, INT MarkKey)
{
  LGM_DOMAIN *theDomain;
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  MESH *mesh;
  LGM_POINT *pointlist;
  float fValue;
  int iValue, D;
  DOUBLE h;
  INT i, j, k, *BNDS_Per_Subdom, *p, coeff, old, surf_id[3], xy_id[3], norp, np;
  char buff[5], name[12], name1[16];
  FILE *stream;

  /* make MarkKey global */
  LGM_MarkKey = MarkKey;

  surface_error = 0;
  old = 0;
  h = 0.0;
  D = 1;
  coeff = 0;
  SURFACE_ADAPT = 1;
  cosAngle = 0.94;
  SAVE_SURFACE = 0;
  for (i=1; i<argc; i++)
  {
    if (argv[i][0] == 'a')                      /* read option for surface adaption */
    {
      if (sscanf(argv[i],"a %d",&iValue) != 1)
        return (NULL);
      SURFACE_ADAPT = iValue;
    }
    if (argv[i][0] == 'h')                      /* read h-option */
    {
      if (sscanf(argv[i],"h %f",&fValue) != 1)
        return (NULL);
      h = fValue;
    }
    if (argv[i][0] == 'c')
    {
      if (sscanf(argv[i],"c %d",&iValue) != 1)
        return (NULL);
      coeff = iValue;
    }
    if (argv[i][0] == 'D')                      /* display all D surface triangles */
    {
      if (sscanf(argv[i],"D %d",&iValue) != 1)
        return (NULL);
      D = iValue;
    }
    if (argv[i][0] == 'w')                      /* resolve angle */
    {
      if (sscanf(argv[i],"w %f",&fValue) != 1)
        return (NULL);
      cosAngle = cos(PI*(180-fValue)/180);
    }
    if (argv[i][0] == 's')                      /* save surfaces and volumes */
    {
      if (sscanf(argv[i],"s %d",&iValue) != 1)
        return (NULL);
      SAVE_SURFACE = iValue;
    }

  }

  if(h<=0.0)
  {
    /* get Coefficientfunctions */
    if (BVP_SetCoeffFct(aBVP,-1,Coefficients))
      return (NULL);
    LOCAL_H = Coefficients[coeff];
  }

  if(h<0.0)
    SURFACE_ADAPT = 0;
  /* set LGM_BVP */
  theDomain = BVP2LGM(aBVP);
  if (theDomain==NULL)
    return (NULL);

  if (LGM_DEBUG)
  {
    PrintDomainInfo(theDomain);
    PrintSubdomainInfo(theDomain);
    /*	for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
       PrintSurfaceInfo(theSurface);	*/
    /*	for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
       PrintLineInfo(theLine);	*/
  }

  if(Aux_Check(Heap, theDomain, LGM_MarkKey))
    return(NULL);

  /* allocate mesh */
  mesh = (MESH *) GetTmpMem(Heap,sizeof(MESH),LGM_MarkKey);
  if (mesh == NULL)
    return(NULL);

  /* init mesh: only surface-mesh */
  mesh->nInnP = 0;
  mesh->Position = NULL;
  mesh->nElements = NULL;
  mesh->nbElements        = NULL;
  mesh->Element_corners = NULL;
  mesh->Element_corner_ids = NULL;
  mesh->nSubDomains = LGM_DOMAIN_NSUBDOM(theDomain);

  mesh->VertexLevel = NULL;
  mesh->VertexPrio = NULL;
  mesh->ElementLevel = NULL;
  mesh->ElementPrio = NULL;

  /* get heap for surface-mesh-substructures: the subdomain-dependence */
  mesh->nSides = (INT *) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT),LGM_MarkKey);
  if (mesh->nSides==NULL)
    return(NULL);
  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
    mesh->nSides[i] = 0;


  mesh->Side_corners = (INT **) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT*),LGM_MarkKey);
  if (mesh->Side_corners == NULL)
    return (NULL);
  for(i=0; i<LGM_DOMAIN_NSUBDOM(theDomain)+1; i++)
    mesh->Side_corners[i] = NULL;
  mesh->xy_Side = (INT **) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT*),LGM_MarkKey);
  if (mesh->xy_Side == NULL)
    return (NULL);
  for(i=0; i<LGM_DOMAIN_NSUBDOM(theDomain)+1; i++)
    mesh->xy_Side[i] = NULL;
  mesh->Side_corner_ids = (INT ***) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT**),LGM_MarkKey);
  if (mesh->Side_corner_ids == NULL)
    return (NULL);

  if((SURFACE_ADAPT)&&(h>0.0))
  {
    Generate_Basis_Mesh(theDomain, h, Heap, D, MarkKey);
  }
  /* search points, which must be resolved by DiscretizeLine and DiscretizeSurface */
  norp = 0;
  np = Count_Points(theDomain);
  pointlist = (LGM_POINT *)GetTmpMem(Heap,(np+1)*sizeof(LGM_POINT),LGM_MarkKey);
  if(pointlist==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  ResolvePoints(Heap, theDomain, pointlist, &norp);

  if(LGM_DEBUG)
    for(i=0; i<norp; i++)
      printf("%d %f %f %f\n", (int)i, pointlist[i].position[0], pointlist[i].position[1],  pointlist[i].position[2]);


  /* prepare for surface-mesh */
  mesh->nBndP = 0;

  /* discretize lines */
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    theLine->ldisc = (LGM_LINEDISC*)GetTmpMem(Heap,sizeof(LGM_LINEDISC)*1,LGM_MarkKey);
    if( (theLine->ldisc==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }
    theLine->ldiscnew = (LGM_LINEDISCNEW*)GetTmpMem(Heap,sizeof(LGM_LINEDISCNEW)*1, MarkKey);
    if( (theLine->ldiscnew==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }
    /*		if (DiscretizeLine(Heap,theLine,h,pointlist, norp, MarkKey))
                            return(NULL);*/
    if (DiscretizeLineNew(Heap,theLine,h,pointlist, norp, MarkKey))
      return(NULL);
    if (LGM_DEBUG)
      PrintLineInfo(theLine);
    if(SAVE_SURFACE)
      Write_Line(theLine);
  }

  /* discretize surfaces */
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    LGM_SURFACE_DISC(theSurface) = (LGM_SURFDISC*)GetTmpMem(Heap,sizeof(LGM_SURFDISC)*1,LGM_MarkKey);
    if( (LGM_SURFACE_DISC(theSurface)==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }

    if(SAVE_SURFACE)
    {
      /* Abspeichern der schon erstellten Oberflaechentriangulierungen */

      name[0] = 's';
      name[1] = 'u';
      name[2] = 'r';
      name[3] = 'f';
      name[4] = 'a';
      name[5] = 'c';
      name[6] = 'e';
      sprintf(buff,"%d",(int)LGM_SURFACE_ID(theSurface));
      name[7] = buff[0];
      name[8] = buff[1];
      name[9] = buff[2];
      name[10] = buff[3];
      name[11] = buff[4];

      name1[0] = 's';
      name1[1] = 'u';
      name1[2] = 'r';
      name1[3] = 'f';
      name1[4] = 'a';
      name1[5] = 'c';
      name1[6] = 'e';
      name1[7] = 'm';
      name1[8] = 'e';
      name1[9] = 's';
      name1[10] = 'h';
      sprintf(buff,"%d",(int)LGM_SURFACE_ID(theSurface));
      name1[11] = buff[0];
      name1[12] = buff[1];
      name1[13] = buff[2];
      name1[14] = buff[3];
      name1[15] = buff[4];
      stream = fopen(name,"r+");
      if (stream==NULL)
      {
        /* fclose(stream); */
        printf("%s %d %s\n", "Surface ", (int)LGM_SURFACE_ID(theSurface), "not triangulated, do now");

        if (DiscretizeSurface(Heap,theSurface,mesh,h,pointlist,norp, D,MarkKey))
          return(NULL);
        if(!surf_err)
          Write_Surface(theSurface, name, name1);
      }
      else
      {
        fclose(stream);
        printf("%s %d\n", "Read Surface ", (int)LGM_SURFACE_ID(theSurface));
        Read_Surface(Heap, theSurface, MarkKey, name);
      }
    }
    else
    {
      if (DiscretizeSurface(Heap,theSurface,mesh,h,pointlist,norp,D,MarkKey))
        return(NULL);
    }

    LGM_SURFDISC_FMESH_ID(LGM_SURFACE_DISC(theSurface)) =
      (INT*)GetTmpMem(Heap,(LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))+1)*sizeof(INT),LGM_MarkKey);
    if(LGM_SURFDISC_FMESH_ID(LGM_SURFACE_DISC(theSurface))==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }

    for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); i++)
      LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), i) = -1;

    if (LGM_DEBUG)
    {
      PrintSurfaceInfo(theSurface);
      printf("%s\n", "#####################################################");
    }

    if(SURFACE_ADAPT)
    {
      if(LGM_SURFACE_ID(theSurface)!=-1)
      {
        LGM_SURFACE_DISC(theSurface)->xy_id = (INT *) GetTmpMem(Heap,(LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface))+1)*sizeof(INT), LGM_MarkKey);
        if(LGM_SURFACE_DISC(theSurface)->xy_id==NULL)
        {
          printf("%s\n","Not enough memory");
          assert(0);
        }
        for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)); i++)
          LGM_SURFDISC_XY_ID(LGM_SURFACE_DISC(theSurface), i) = -1;

        for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(xy_Surface)); i++)
        {
          if(xy_element_id[i])
          {
            for(k=0; k<3; k++)
              xy_id[k] = xy_new_id[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, k)]-1;
            for(j=0; j<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)); j++)
            {
              for(k=0; k<3; k++)
                surf_id[k] = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), j, k);
              if( ( (xy_id[0]==surf_id[0]) && (xy_id[1]==surf_id[1]) && (xy_id[2]==surf_id[2]) )
                  ||  ( (xy_id[0]==surf_id[1]) && (xy_id[1]==surf_id[2]) && (xy_id[2]==surf_id[0]) )
                  ||  ( (xy_id[0]==surf_id[2]) && (xy_id[1]==surf_id[0]) && (xy_id[2]==surf_id[1]) )
                  ||  ( (xy_id[1]==surf_id[0]) && (xy_id[0]==surf_id[1]) && (xy_id[2]==surf_id[2]) )
                  ||  ( (xy_id[1]==surf_id[1]) && (xy_id[0]==surf_id[2]) && (xy_id[2]==surf_id[0]) )
                  ||  ( (xy_id[1]==surf_id[2]) && (xy_id[0]==surf_id[0]) && (xy_id[2]==surf_id[1]) ) )
                LGM_SURFDISC_XY_ID(LGM_SURFACE_DISC(theSurface), j) = i;
            }
          }
        }
      }
    }
  }

  if(surface_error)
  {
    UserWriteF("\n");
    UserWriteF("Error in Surface Triangulation !\n");
    UserWriteF("Check the Input Data !\n");
    return(NULL);
  }
  /* discretize domain */
  if (DiscretizeDomain(Heap,theDomain,mesh,h, MarkKey))
    return(NULL);

  if (Get_NBNDS_Per_Subdomain(Heap,theDomain,&BNDS_Per_Subdom,h))
    return (NULL);

  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    mesh->Side_corners[i]    = (INT *)  GetTmpMem(Heap,sizeof(INT)*(BNDS_Per_Subdom[i]+1),LGM_MarkKey);
    if (mesh->Side_corners[i]==NULL)
      return (NULL);
    mesh->xy_Side[i] = (INT *)  GetTmpMem(Heap,sizeof(INT)*(BNDS_Per_Subdom[i]+1),LGM_MarkKey);
    if (mesh->xy_Side[i]==NULL)
      return (NULL);
    for(j=0; j<BNDS_Per_Subdom[i]; j++)
      mesh->xy_Side[i][j] = -1;
    mesh->Side_corner_ids[i] = (INT **) GetTmpMem(Heap,sizeof(INT*)*(BNDS_Per_Subdom[i]+1),LGM_MarkKey);
    if (mesh->Side_corner_ids[i]==NULL)
      return (NULL);
    p = (INT *) GetTmpMem(Heap,3*sizeof(INT)*BNDS_Per_Subdom[i],LGM_MarkKey);
    if(p==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
    for (j=0; j<BNDS_Per_Subdom[i]; j++)
      mesh->Side_corner_ids[i][j] = p+3*j;
  }

  /* discretize surfaces */
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    if (TransferSurfaces2Mesh(Heap,theSurface,mesh,h))
      return(NULL);
  }

  mesh->VertexLevel = NULL;
  mesh->VertexPrio = NULL;
  mesh->ElemSideOnBnd = NULL;

  /* print mesh-info */
  if (LGM_DEBUG)
    if (PrintMeshInfo(mesh))
      return (NULL);

  if (LGM_DEBUG) printf("%s\n","und tschuessss");

  return (mesh);
}

static INT InnerPointsPerSurfaceSegment (LGM_SURFACE *theSurface, DOUBLE h, INT i, INT *n)
{
  INT j,npointsall,innerpoints;

  npointsall = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));
  innerpoints = npointsall;
  for(j=0; j<LGM_SURFACE_NLINE(theSurface); j++)
  {
    innerpoints = innerpoints - ( LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(LGM_SURFACE_LINE(theSurface,j))) - 1);
  }
  *n = innerpoints;
  return (0);
}

static INT Get_NBNDS_Per_Subdomain (HEAP *Heap, LGM_DOMAIN *theDomain, INT **sizes, DOUBLE h)
{
  INT i,*p,n;
  LGM_SURFACE *theSurface;

  /* get heap for information */
  p = (INT *) GetTmpMem(Heap,sizeof(INT)*(LGM_DOMAIN_NSUBDOM(theDomain)+1),LGM_MarkKey);
  if (p == NULL) return(1);
  *sizes = p;

  /* get information */
  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
    p[i]=0;
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    /*		n = 0;
                    if (InnerPointsPerSurfaceSegment(theSurface,h,i,&pn))
                            return (1);
                    n += pn;*/
    n = LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface));
    p[LGM_SURFACE_LEFT(theSurface)]  += LGM_SURFACE_NPOINT(theSurface)-1+n;
    p[LGM_SURFACE_RIGHT(theSurface)] += LGM_SURFACE_NPOINT(theSurface)-1+n;
  }

  return (0);
}

static INT Get_NBNDP (LGM_DOMAIN *theDomain, INT *nBND, DOUBLE h)
{
  INT n;

  n = LGM_DOMAIN_NPOINT(theDomain);
  *nBND = n;

  return (0);
}


/* OS_CHANGED **********************************************************************************/
#ifdef NEW_LGM
#define SPAT(x,y,z)                     (((x[1])*(y[2])-(x[2])*(y[1]))*z[0] + \
                                         ((x[2])*(y[0])-(x[0])*(y[2]))*z[1] + \
                                         ((x[0])*(y[1])-(x[1])*(y[0]))*z[2])

static DOUBLE Calc_Local_Coord(DOUBLE *p0, DOUBLE *p1, DOUBLE *p2, DOUBLE *x, DOUBLE *l)
{
  DOUBLE a0[3], a1[3], n[3], p[3], xmp2[3], mag, invmag, dist;

  /* two edges */
  V3_SUBTRACT(p0,p2,a0);
  V3_SUBTRACT(p1,p2,a1);

  /* normal and half triangle area (mag) */
  V3_VECTOR_PRODUCT(a0,a1,n);
  V3_EUKLIDNORM(n,mag);
  invmag = 1.0/mag;
  V3_SCALE(invmag,n);

  /* project x onto triangle plane */
  V3_SUBTRACT(x,p2,xmp2);
  V3_SCALAR_PRODUCT(n,xmp2,dist);
  V3_SCALESET(-dist,n,p);
  V3_ADD1(xmp2,p);

  /* local co-ordinates */
  l[0] = SPAT(p,a1,n)*invmag;
  l[1] = SPAT(a0,p,n)*invmag;
  l[2] = 1.0 - l[1] - l[0];

  return fabs(dist);
}

#ifdef LGM_ACCELERATE
static DOUBLE PointLGMTriangleDistance(DOUBLE *x, void *obj)
{
  LGM_TRIANGLE *theTriangle;
  DOUBLE d, dmin, d_vec[3],  alpha, lambda[3], a[3], b[3], new_x[3];
  DOUBLE *p0, *p1, *p2;
  INT i, found = 0;

  theTriangle = (LGM_TRIANGLE *) obj;
  dmin = DBL_MAX;

  /* Check if point is close to its orthog. projection onto this triangle */
  p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,0);
  p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,1);
  p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,2);
  d = Calc_Local_Coord(p0,p1,p2,x,lambda);
  if ((lambda[0]>=-SMALL1) && (lambda[1]>=-SMALL1) && (lambda[2]>=-SMALL1))
    dmin = d;
  else
  {
    /* Check if point is close to its orthog. prj. onto one of the edges */
    for (i=0; i<3; i++)
    {
      p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,i);
      p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,(i+1)%3);
      V3_SUBTRACT(x,p0,a);
      V3_SUBTRACT(p1,p0,b);
      alpha = V3_SCAL_PROD(a,b)/V3_SCAL_PROD(b,b);
      if (alpha >= 0.0 && alpha <= 1.0)
      {
        found = 1;
        V3_SCALESET(alpha,b,new_x);
        V3_ADD1(p0,new_x);
        V3_SUBTRACT(x,new_x,d_vec);
        V3_EUKLIDNORM(d_vec,d);
        dmin = MIN(dmin,d);
      }
    }

    if (!found)
    {
      /* Get distance to closest node */
      for (i=0; i<3; i++)
      {
        p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,i);
        V3_SUBTRACT(x,p0,d_vec);
        V3_EUKLIDNORM(d_vec,d);
        dmin = MIN(dmin,d);
      }
    }
  }

  return dmin;
}

INT GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local, DOUBLE *dummy)
{
  LGM_TRIANGLE *theTriangle;
  DOUBLE dist, min_dist, lambda[3], *p0, *p1, *p2, new_global[3];
  DOUBLE dist_vec[3], min_global[3], a[3], b[3], alpha;
  INT tri, i, found;

  dist = BBT_TreePointDistance(theSurface->bbtree,global,(void **)&theTriangle,
                               PointLGMTriangleDistance);
  if (dist == DBL_MAX) return -1;

  tri = theTriangle - LGM_SURFACE_TRIANGLE(theSurface,0);       /* Get triangle ID */
  assert(LGM_SURFACE_TRIANGLE(theSurface,tri)==theTriangle);
  found = 0;
  min_dist = DBL_MAX;

  /* Check if point is close to its orthog. projection onto this triangle */
  p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,0);
  p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,1);
  p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,2);
  dist = Calc_Local_Coord(p0,p1,p2,global,lambda);

  if (!((lambda[0]>=-SMALL1) && (lambda[1]>=-SMALL1) && (lambda[2]>=-SMALL1)))
  {
    /* Check if point is close to its orthog. prj. onto one of the edges */
    for (i=0; i<3; i++)
    {
      p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,i);
      p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,(i+1)%3);
      V3_SUBTRACT(global,p0,a);
      V3_SUBTRACT(p1,p0,b);
      alpha = V3_SCAL_PROD(a,b)/V3_SCAL_PROD(b,b);
      if (alpha >= 0.0 && alpha <= 1.0)
      {
        found = 1;
        V3_SCALESET(alpha,b,new_global);
        V3_ADD1(p0,new_global);
        V3_SUBTRACT(global,new_global,dist_vec);
        V3_EUKLIDNORM(dist_vec,dist);
        if (dist < min_dist)                         /* Currently the nearest */
        {
          min_dist = dist;
          V3_COPY(new_global,min_global);
        }
      }
    }

    if (!found)
    {
      /* Take closest node of the triangle as new point */
      for (i=0; i<3; i++)
      {
        p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,i);
        V3_SUBTRACT(global,p0,dist_vec);
        V3_EUKLIDNORM(dist_vec,dist);
        if (dist < min_dist)                         /* Currently the nearest */
        {
          min_dist = dist;
          V3_COPY(p0,min_global);
        }
      }
    }

    /* Compute local */
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,2);
    dist = Calc_Local_Coord(p0,p1,p2,min_global,lambda);

  }

  if (lambda[0] < 0.0) lambda[0] = 0.0;
  if (lambda[1] < 0.0) lambda[1] = 0.0;

  local[0] = lambda[0] + tri;
  local[1] = lambda[1] + tri;

#if 0 //OS_DEBUG
  {
    printf("# lambda =(%f,%f,%f) \n",lambda[0],lambda[1],lambda[1]);
    printf("  x      =(%f,%f,%f)\n",global[0],global[1],global[2]);
    dist = Calc_Local_Coord(p0,p1,p2,min_global,lambda);
    printf("  dist   = %f, tri=%d\n",dist,tri);
    printf("  lambda =(%f,%f,%f) \n",lambda[0],lambda[1],lambda[1]);
    printf("  xmin   =(%f,%f,%f)\n",min_global[0],min_global[1],min_global[2]);
    printf("  p0     =(%f,%f,%f)\n",p0[0],p0[1],p0[2]);
    printf("  p1     =(%f,%f,%f)\n",p1[0],p1[1],p1[2]);
    printf("  p2     =(%f,%f,%f)\n",p2[0],p2[1],p2[2]);
  }
#endif

  return(tri);
}

#else /* no LGM acceleration */
INT GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local, DOUBLE *dummy)
{
  LGM_TRIANGLE *theTriangle;
  DOUBLE dist, min_dist, lambda[3], min_lambda[3], *p0, *p1, *p2, new_global[3];
  DOUBLE dist_vec[3], min_global[3], a[3], b[3], alpha;
  INT tri, min_tri, i, found;
  enum {none,triangle,edge,node} obj;

  min_dist = DBL_MAX;       /* Largest value possible */
  min_tri = -1;
  obj = none;

  for(tri=0; tri<LGM_SURFACE_NTRIANGLE(theSurface); tri++)
  {
    found = 0;
    theTriangle = LGM_SURFACE_TRIANGLE(theSurface,tri);

    /* Check if point is close to its orthog. projection onto this triangle */
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,2);

    dist = Calc_Local_Coord(p0,p1,p2,global,lambda);
    if ((lambda[0]>=-SMALL1) && (lambda[1]>=-SMALL1) && (lambda[2]>=-SMALL1))
    {
      if (0 && dist < SMALL1)                   /* The point is already on the triangle OS_TODO: use global citerium, not SMALL1, which is local*/
      {
        if (lambda[0] < 0.0) lambda[0] = 0.0;
        if (lambda[1] < 0.0) lambda[1] = 0.0;
        local[0] = lambda[0] + tri;
        local[1] = lambda[1] + tri;
        return(tri);
      }
      found = 1;
      if (dist < min_dist)                   /* Currently the nearest */
      {
        min_dist = dist;
        min_tri = tri;
        obj = triangle;
        V3_COPY(lambda,min_lambda);
      }
    }

    if (!found)
    {
      /* Check if point is close to its orthog. prj. onto one of the edges */
      for (i=0; i<3; i++)
      {
        p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,i);
        p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,(i+1)%3);
        V3_SUBTRACT(global,p0,a);
        V3_SUBTRACT(p1,p0,b);
        alpha = V3_SCAL_PROD(a,b)/V3_SCAL_PROD(b,b);
        if (alpha >= 0.0 && alpha <= 1.0)
        {
          found = 1;
          V3_SCALESET(alpha,b,new_global);
          V3_ADD1(p0,new_global);
          V3_SUBTRACT(global,new_global,dist_vec);
          V3_EUKLIDNORM(dist_vec,dist);
          if (dist < min_dist)                               /* Currently the nearest */
          {
            min_dist = dist;
            min_tri = tri;
            obj = edge;
            V3_COPY(new_global,min_global);
          }
        }
      }
    }

    if (!found)
    {
      /* Check if point is close to one of the nodes */
      for (i=0; i<3; i++)
      {
        p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,i);
        V3_SUBTRACT(global,p0,dist_vec);
        V3_EUKLIDNORM(dist_vec,dist);
        if (dist < min_dist)                         /* Currently the nearest */
        {
          min_dist = dist;
          min_tri = tri;
          obj = node;
          V3_COPY(p0,min_global);
        }
      }
    }

  }

  /* Compute local */
  if (obj != triangle)
  {
    theTriangle = LGM_SURFACE_TRIANGLE(theSurface,min_tri);
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(theTriangle,2);
    Calc_Local_Coord(p0,p1,p2,min_global,min_lambda);
  }

  if (min_lambda[0] < 0.0) min_lambda[0] = 0.0;
  if (min_lambda[1] < 0.0) min_lambda[1] = 0.0;
  local[0] = min_lambda[0] + min_tri;
  local[1] = min_lambda[1] + min_tri;

#if 0   //OS_DEBUG
  {
    dist = Calc_Local_Coord(p0,p1,p2,min_global,min_lambda);
    if (obj == triangle) printf("triangle case\n");
    else if (obj == edge) printf("edge case\n");
    else if (obj == node) printf("node case\n");
    else printf("????? case\n");

    printf("# lambda =(%f,%f,%f) \n",,min_lambda[0],min_lambda[1],min_lambda[1]);
    printf("  x      =(%f,%f,%f)\n",global[0],global[1],global[2]);
    dist = Calc_Local_Coord(p0,p1,p2,min_global,lambda);
    printf("  dist   = %f, tri=%d\n",dist,tri);
    printf("  lambda =(%f,%f,%f) \n",lambda[0],lambda[1],lambda[1]);
    printf("  xmin   =(%f,%f,%f)\n",min_global[0],min_global[1],min_global[2]);
    printf("  p0     =(%f,%f,%f)\n",p0[0],p0[1],p0[2]);
    printf("  p1     =(%f,%f,%f)\n",p1[0],p1[1],p1[2]);
    printf("  p2     =(%f,%f,%f)\n",p2[0],p2[1],p2[2]);
    printf("  mindist = %f, tri=%d\n",min_dist,tri);
  }
#endif

  assert(min_tri != -1);
  return(min_tri);
}
#endif

/* OS_CHANGED **********************************************************************************/
#else /* use old LGM */

static DOUBLE Calc_Local_Coord(DOUBLE *p0, DOUBLE*p1, DOUBLE *p2, DOUBLE *global, DOUBLE *lam)
{
  INT i,j;
  DOUBLE e0[3],e1[3],e2[3];
  DOUBLE n0[3],n1[3],n2[3],p[3],np[3],hp[3],l;
  DOUBLE a[9], det;
  DOUBLE aa[4],bb[2],cc[2];

  Minus(e0,global,p0);
  Minus(e1,global,p1);
  Minus(e2,global,p2);

  /* compute 2 vectors in the plane + normalvector */
  Minus(n0,e2,e0);
  V3_EUKLIDNORM(n0,l);
  Scale(n0,1.0/l);

  Minus(n1,e2,e1);
  V3_EUKLIDNORM(n1,l);
  Scale(n1,1.0/l);

  Cross(n2,n0,n1);
  V3_EUKLIDNORM(n2,l);
  Scale(n2,1.0/l);

  p[0] = global[0];
  p[1] = global[1];
  p[2] = global[2];

  Minus(np,p,p0);

  hp[0] = n2[0];
  hp[1] = n2[1];
  hp[2] = n2[2];

  Scale(hp,Mult(n2,np));
  Minus(np,np,hp);
  Plus(p,np,p0);

  /* calculate local coordinates */

  a[0] = p0[0];
  a[1] = p1[0];
  a[2] = p2[0];

  a[3] = p0[1];
  a[4] = p1[1];
  a[5] = p2[1];

  a[6] = p0[2];
  a[7] = p1[2];
  a[8] = p2[2];

  aa[0] = a[0] - a[2];
  aa[1] = a[1] - a[2];
  aa[2] = a[3] - a[5];
  aa[3] = a[4] - a[5];

  det = Det2d(aa);
  if(Det2d(aa)!=0)
  {
    cc[0] = p[0] - a[2];
    cc[1] = p[1] - a[5];
    InvMatMult2d(bb,cc,aa);
  }
  else
  {
    aa[0] = a[0] - a[2];
    aa[1] = a[1] - a[2];
    aa[2] = a[6] - a[8];
    aa[3] = a[7] - a[8];
    det = Det2d(aa);
    if(Det2d(aa)!=0)
    {
      cc[0] = p[0] - a[2];
      cc[1] = p[2] - a[8];
      InvMatMult2d(bb,cc,aa);
    }
    else
    {
      aa[0] = a[3] - a[5];
      aa[1] = a[4] - a[5];
      aa[2] = a[6] - a[8];
      aa[3] = a[7] - a[8];
      det = Det2d(aa);
      if(Det2d(aa)!=0)
      {
        cc[0] = p[1] - a[5];
        cc[1] = p[2] - a[8];
        InvMatMult2d(bb,cc,aa);
      }
      else
        printf("%s\n","E R R O R 1");
    }
  }
  lam[0] = bb[0];
  lam[1] = bb[1];
  lam[2] = 1 - bb[0] - bb[1];

  V3_EUKLIDNORM(hp,l);
  return l;
}

static INT NormalVector(LGM_TRIANGLE *theSurfaceTriangle, DOUBLE *n)
{
  DOUBLE *p0,*p1,*p2,n0[3],n1[3];

  p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(theSurfaceTriangle,0);
  p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(theSurfaceTriangle,1);
  p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(theSurfaceTriangle,2);

  V_DIM_SUBTRACT(p2, p0, n0);
  V_DIM_SCALE(1.0/sqrt(V_DIM_SCAL_PROD(n0, n0)), n0);
  V_DIM_SUBTRACT(p2, p1, n1);
  V_DIM_SCALE(1.0/sqrt(V_DIM_SCAL_PROD(n1, n1)), n1);

  V_DIM_VECTOR_PRODUCT(n0, n1, n);
  V_DIM_SCALE(1.0/sqrt(V_DIM_SCAL_PROD(n, n)), n);

  return(0);
}

static INT CASE1(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *lam, DOUBLE *d, DOUBLE *n)
{
  DOUBLE dist;
  INT i, mi;
  DOUBLE *p0,*p1,*p2, ne[3], angle, sp;

  /* Punkt liegt auf der Surface
   * 0<=lamda_i<=1,  dist < eps
   * Falls der Punkt genau auf einer Kante zwischen 2 Dreiecken liegt,
   * ist es egal welches genommen wird
   */
  dist = DBL_MAX;
  mi = -1;
  d[0] = DBL_MAX;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    NormalVector(LGM_SURFACE_TRIANGLE(theSurface,i), ne);
    if(sqrt(V_DIM_SCAL_PROD(n, n))>SMALL)
    {
      sp = V_DIM_SCAL_PROD(ne, n);
      if(sp>1-SMALL)
        angle = 0.0;
      else
        angle = acos(sp);
    }
    else
      angle = 0.0;

    if(angle<PI*Triangle_Angle2/180)
    {
      p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
      p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
      p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

      dist = Calc_Local_Coord(p0, p1, p2, global, lam);
      if( (lam[0]>=-SMALL1) && (lam[1]>=-SMALL1) && (lam[2]>=-SMALL1) && (dist<SMALL1) )
      {
        mi = i;
        d[0] = dist;
        break;
      }
    }
  }
  return(mi);
}

static INT CASE2(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *lam, DOUBLE *n)
{
  DOUBLE dist1, dist2, min, d1, d2;
  INT i, j, mi1, mi2, mi;
  DOUBLE *p0,*p1,*p2, mm;
  DOUBLE A, B, m, point[3], dist_vec[3],pp0[3], pp1[3], new_global[3], local[3], angle, sp, ne[3];

  /* Punkt liegt nicht auf der Surface
   * Projeziere den Punkt in ein Dreieck mit  0<=lamda_i<=1 und dist minimal
   * oder bestimme einen Punkt auf einer Kante,  sodass der Abstand
   * noch kleiner wird
   */

  /* kuerzeste Abstand zur einem Dreieck */
  min = DBL_MAX;
  mi1 = -1;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    NormalVector(LGM_SURFACE_TRIANGLE(theSurface,i), ne);
    if(sqrt(V_DIM_SCAL_PROD(n, n))>SMALL)
    {
      sp = V_DIM_SCAL_PROD(ne, n);
      if(sp>1-SMALL)
        angle = 0.0;
      else
        angle = acos(sp);
    }
    else
      angle = 0.0;

    if(angle<PI*Triangle_Angle2/180)
    {
      p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
      p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
      p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

      dist1 = Calc_Local_Coord(p0, p1, p2, global, lam);
      if( (lam[0]>=-SMALL1) && (lam[1]>=-SMALL1) && (lam[2]>=-SMALL1) && (dist1<min) )
      {
        mi1 = i;
        min = dist1;
      }
    }
  }
  dist1 = min;

  /* kuerzeste Abstand zur eine Kante */
  /*	mi1 = -1;*/
  mi2 = -1;
  min = DBL_MAX;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    NormalVector(LGM_SURFACE_TRIANGLE(theSurface,i), ne);
    if(sqrt(V_DIM_SCAL_PROD(n, n))>SMALL)
    {
      sp = V_DIM_SCAL_PROD(ne, n);
      if(sp>1-SMALL)
        angle = 0.0;
      else
        angle = acos(sp);
    }
    else
      angle = 0.0;

    if(angle<PI*Triangle_Angle2/180)
    {
      for(j=0; j<3; j++)
      {
        p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),j%3);
        p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3);

        pp0[0] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j)%3))[0];
        pp0[1] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j)%3))[1];
        pp0[2] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j)%3))[2];
        pp1[0] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3))[0];
        pp1[1] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3))[1];
        pp1[2] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3))[2];

        A =   (pp1[0]-pp0[0]) * (global[0]-pp0[0])
            + (pp1[1]-pp0[1]) * (global[1]-pp0[1])
            + (pp1[2]-pp0[2]) * (global[2]-pp0[2]);
        B =   (pp1[0]-pp0[0]) * (pp1[0]-pp0[0])
            + (pp1[1]-pp0[1]) * (pp1[1]-pp0[1])
            + (pp1[2]-pp0[2]) * (pp1[2]-pp0[2]);

        m = A / B;

        if((0.0<=m) && (m<=1.0))
        {
          point[0] = pp0[0] + m * (pp1[0]-pp0[0]);
          point[1] = pp0[1] + m * (pp1[1]-pp0[1]);
          point[2] = pp0[2] + m * (pp1[2]-pp0[2]);

          dist_vec[0] = global[0] - point[0];
          dist_vec[1] = global[1] - point[1];
          dist_vec[2] = global[2] - point[2];

          dist2 =sqrt(dist_vec[0]*dist_vec[0]+dist_vec[1]*dist_vec[1]+dist_vec[2]*dist_vec[2]);
        }
        else
        {
          dist_vec[0] = global[0] - pp0[0];
          dist_vec[1] = global[1] - pp0[1];
          dist_vec[2] = global[2] - pp0[2];

          d1 =sqrt(dist_vec[0]*dist_vec[0]+dist_vec[1]*dist_vec[1]+dist_vec[2]*dist_vec[2]);

          dist_vec[0] = global[0] - pp1[0];
          dist_vec[1] = global[1] - pp1[1];
          dist_vec[2] = global[2] - pp1[2];

          d2 =sqrt(dist_vec[0]*dist_vec[0]+dist_vec[1]*dist_vec[1]+dist_vec[2]*dist_vec[2]);
          if(d1>d2)
          {
            m = 1.0;
            dist2 = d2;
            point[0] = pp1[0];
            point[1] = pp1[1];
            point[2] = pp1[2];
          }
          else
          {
            m = 0.0;
            dist2 = d1;
            point[0] = pp0[0];
            point[1] = pp0[1];
            point[2] = pp0[2];
          }
        }
        if( (dist2<min) )
        {
          mi2 = i;
          min = dist2;
          mm = m;
          new_global[0] = point[0];
          new_global[1] = point[1];
          new_global[2] = point[2];
        }
      }
    }
  }
  dist2 = min;

  if( (dist1<dist2)&&(mi1!=-1) )
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi1),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi1),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi1),2);

    Calc_Local_Coord(p0, p1, p2, global, lam);
    return(mi1);
  }
  else
  {
    if(mi2!=-1)
    {
      /*		mi = GetLocalKoord(theSurface, new_global, local, n);*/
      p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi2),0);
      p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi2),1);
      p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi2),2);

      Calc_Local_Coord(p0, p1, p2, new_global, lam);
      return(mi2);

    }
  }
  return(-1);
}

int GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local, DOUBLE *n)
{
  INT i,j,test,min,mi, dist_i;
  DOUBLE *p0,*p1,*p2,e0[3],e1[3],e2[3],eps;
  DOUBLE n0[3],n1[3],n2[3],p[3],np[3],hp[3],l, d;
  DOUBLE a[9],b[3],c[3];
  DOUBLE aa[4],bb[2],cc[2];
  DOUBLE lam[3], dist, min_dist, new_lam[3], A, B, m, point[3], dist_vec[3],pp0[3], pp1[3], new_global[3];

  min = INT_MAX;
  min_dist = DBL_MAX;
  dist_i = -1;
  mi = -1;

  mi = CASE1(theSurface, global, lam, &d, n);
  if(mi==-1)
    mi = CASE2(theSurface, global, lam, n);
  /*	if(mi==-1)
                  printf("%s\n", "schotter");	*/

  if( (lam[0]<0.0) || (lam[1]<0.0) )
  {
    if( (lam[0]<0.0) && (lam[0]>-SMALL) )
      lam[0] = 0.0;
    if( (lam[1]<0.0) && (lam[1]>-SMALL) )
      lam[1] = 0.0;
  }
  if(lam[0]<0.0) lam[0] = 0.0;
  if(lam[1]<0.0) lam[1] = 0.0;
  local[0] = lam[0] + mi;
  local[1] = lam[1] + mi;

  return(mi);
}
#endif

INT Project2Surface(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local, DOUBLE *n)
{
  INT i,mi;
  DOUBLE *p0,*p1,*p2;
  DOUBLE lam[3], dist;

  dist = DBL_MAX;
  mi = -1;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

    dist = Calc_Local_Coord(p0, p1, p2, global, lam);
    if( (lam[0]>=-SMALL1) && (lam[1]>=-SMALL1) && (lam[2]>=-SMALL1))
    {
      mi = i;
      break;
    }
  }

  if( (lam[0]<0.0) || (lam[1]<0.0) )
  {
    if( (lam[0]<0.0) && (lam[0]>-SMALL) )
      lam[0] = 0.0;
    if( (lam[1]<0.0) && (lam[1]>-SMALL) )
      lam[1] = 0.0;
  }

  local[0] = lam[0] + mi;
  local[1] = lam[1] + mi;

  return(mi);
}

DOUBLE Check_Surface(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  DOUBLE dist, d;
  INT i, mi;
  DOUBLE *p0,*p1,*p2, lam[3], lam1[3];

  dist = DBL_MAX;
  mi = -1;
  d = DBL_MAX;

  lam[0] = -1.0;
  lam[1] = -1.0;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

    dist = Calc_Local_Coord(p0, p1, p2, global, lam1);
    if( (lam1[0]>=-SMALL) && (lam1[1]>=-SMALL) && (lam1[2]>=-SMALL) && (d>dist))
    {
      mi = i;
      d = dist;
      lam[0] = lam1[0];
      lam[1] = lam1[1];
      lam[2] = lam1[2];
    }
  }

  /* if mi<0, no triangle has been found. hence, lam[] isn't valid! */
  if (mi==-1)
    return(d);

  if( (lam[0]<0.0) || (lam[1]<0.0) )
  {
    if( (lam[0]<0.0) && (lam[0]>-SMALL) )
      lam[0] = 0.0;
    if( (lam[1]<0.0) && (lam[1]>-SMALL) )
      lam[1] = 0.0;
  }

  local[0] = lam[0] + mi;
  local[1] = lam[1] + mi;

  return(d);
}

static INT Compare_Points(LGM_POINT *p, DOUBLE *global)
{
  DOUBLE d;

  d = sqrt( (LGM_POINT_POS(p)[0]-global[0])*(LGM_POINT_POS(p)[0]-global[0])
            +         (LGM_POINT_POS(p)[1]-global[1])*(LGM_POINT_POS(p)[1]-global[1])
            +         (LGM_POINT_POS(p)[2]-global[2])*(LGM_POINT_POS(p)[2]-global[2]) );
  if( d < SMALL )
    return(1);
  else
    return(0);
}

static INT DiscretizeDomain (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h, INT MarkKey)
{
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  INT i,j;
  INT *nRef, *newID, *nRefLines;
  INT flag, dummy;
  DOUBLE global[3],local[2], line_local, local_left, local_right, local1, local2;
#ifdef NO_PROJECT
  DOUBLE global_left[3], global_right[3];
#endif
  LGM_POINT *PointList, Point;
  INT npoints, npsurface, nPointList, nbndpoints, ID;
  LINEPOINT *help;

  /* allocate Memory for the mesh */
  npoints = 0;
  for(theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
    npoints = npoints + LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine));
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    npsurface = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      npsurface = npsurface - ( LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(LGM_SURFACE_LINE(theSurface,i))) - 1);
    npoints = npoints + npsurface;
  }

  PointList = (LGM_POINT *) GetTmpMem(Heap,(npoints+1)*sizeof(LGM_POINT),LGM_MarkKey);
  if (PointList==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  nRef = (INT *) GetTmpMem(Heap,(npoints+1)*sizeof(INT),LGM_MarkKey);
  if (nRef==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  newID = (INT *) GetTmpMem(Heap,(npoints+1)*sizeof(INT),LGM_MarkKey);
  if (newID==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  nRefLines = (INT *) GetTmpMem(Heap,(npoints+1)*sizeof(INT),LGM_MarkKey);
  if (nRefLines==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  for (i=0; i<=npoints; i++)
  {
    nRef[i] = 0;
    newID[i] = -1;
    nRefLines[i] = 0;
  }

  nPointList = 0;
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    nbndpoints = 0;
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      nbndpoints = nbndpoints + LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(LGM_SURFACE_LINE(theSurface, i)));

    for(j=0; j<MIN(nbndpoints, LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))); j++)
    {
      local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,0);
      local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,1);
      Surface_Local2Global(theSurface,global,local);
      flag = -1;
      for(i=0; i<nPointList; i++)
      {
        Point = PointList[i];
        if(Compare_Points(&Point, global)==1)
        {
          nRef[i]++;
          flag = 1;
          LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), j) = i;
        }
      }
      if(flag==-1)
      {
        PointList[nPointList].position[0] =  global[0];
        PointList[nPointList].position[1] =  global[1];
        PointList[nPointList].position[2] =  global[2];
        LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), j) = nPointList;
        nPointList++;
        nRef[i]++;
      }
    }
    for(j=MIN(nbndpoints, LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))); j<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); j++)
    {
      local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,0);
      local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,1);
      Surface_Local2Global(theSurface,global,local);
      PointList[nPointList].position[0] =  global[0];
      PointList[nPointList].position[1] =  global[1];
      PointList[nPointList].position[2] =  global[2];
      LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), j) = nPointList;
      nRef[nPointList]++;
      nPointList++;
    }
  }

  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    if(LGM_LINE_BEGIN(theLine)==LGM_LINE_END(theLine))
      dummy = 1;
    else
      dummy = 0;
    help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
    for(j=0; j<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-dummy; j++)
    {
      line_local = help->local;
      Line_Local2GlobalNew(theLine, global, line_local);
      help = help->next;
      flag = 0;
      for(i=0; i<nPointList; i++)
      {
        Point = PointList[i];
        if(Compare_Points(&Point, global)==1)
        {
          nRefLines[i]++;
          flag++;
        }
      }
      assert(flag<2);
    }
  }

  /*	for(i=0;i<nPointList;i++)
                  printf("%d %lf %lf %lf %d %d\n", i,
                                                                                  PointList[i].position[0],
                                                                                  PointList[i].position[1],
                                                                                  PointList[i].position[2],
                                                                                  nRef[i],
                                                                                  nRefLines[i]);*/

  theMesh->theBndPs = (BNDP**)GetFreelistMemory(Heap,sizeof(LGM_BNDP*)*(nPointList+1));
  if (theMesh->theBndPs == NULL)
    return (1);

  for(i=0; i<nPointList; i++)
  {
    theMesh->theBndPs[i] = (BNDP*)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
    if(theMesh->theBndPs[i]==NULL)
      return(1);
    if(nRefLines[i]>0)
    {
      BNDP2LGM(theMesh->theBndPs[i])->Line = (LGM_BNDP_PLINE*)GetFreelistMemory(Heap,(nRefLines[i]+1)*sizeof(LGM_BNDP_PLINE));
      if(BNDP2LGM(theMesh->theBndPs[i])->Line==NULL)
        return(1);
    }
    else
      BNDP2LGM(theMesh->theBndPs[i])->Line = NULL;
    LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i])) = 0;
    BNDP2LGM(theMesh->theBndPs[i])->Surf = (LGM_BNDP_PSURFACE*)GetFreelistMemory(Heap,(nRef[i]+1)*sizeof(LGM_BNDP_PSURFACE));
    if((nRef[i]>1)&&(nRefLines[i]==0))
      assert(0);
    if(BNDP2LGM(theMesh->theBndPs[i])->Surf==NULL)
      return(1);
    LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])) = 0;
  }
  theMesh->nBndP = nPointList;

  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    nbndpoints = 0;
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      nbndpoints = nbndpoints + LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(LGM_SURFACE_LINE(theSurface, i)));

    for(j=0; j<MIN(nbndpoints, LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))); j++)
    {
      local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,0);
      local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,1);
      Surface_Local2Global(theSurface,global,local);

      for(i=0; i<nPointList; i++)
      {
        Point = PointList[i];
        if(Compare_Points(&Point, global)==1)
        {
                                        #ifdef NO_PROJECT
          LGM_BNDP_GLOBAL(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])))[0] = global[0];
          LGM_BNDP_GLOBAL(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])))[1] = global[1];
          LGM_BNDP_GLOBAL(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])))[2] = global[2];
                                        #else
          LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])))[0] = local[0];
          LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])))[1] = local[1];
                                        #endif
          LGM_BNDP_SURFACE(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i]))) = theSurface;
          LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i]))++;
          break;
        }

      }
    }
    for(j=MIN(nbndpoints, LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))); j<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); j++)
    {
      local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,0);
      local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,1);
      Surface_Local2Global(theSurface,global,local);
      ID = LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), j);
                        #ifdef NO_PROJECT
      LGM_BNDP_GLOBAL(BNDP2LGM(theMesh->theBndPs[ID]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[ID])))[0] = global[0];
      LGM_BNDP_GLOBAL(BNDP2LGM(theMesh->theBndPs[ID]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[ID])))[1] = global[1];
      LGM_BNDP_GLOBAL(BNDP2LGM(theMesh->theBndPs[ID]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[ID])))[2] = global[2];
                        #else
      LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[ID]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[ID])))[0] = local[0];
      LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[ID]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[ID])))[1] = local[1];
                        #endif
      LGM_BNDP_SURFACE(BNDP2LGM(theMesh->theBndPs[ID]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[ID]))) = theSurface;
      LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[ID]))++;
    }
  }

  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    if(LGM_LINE_BEGIN(theLine)==LGM_LINE_END(theLine))
      dummy = 1;
    else
      dummy = 0;
    help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
    local1 = local2 = 0;
    for(j=0; j<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-dummy; j++)
    {
      line_local = help->local;
      Line_Local2GlobalNew(theLine, global, line_local);
      flag = 0;
      if(j>0)
        local_left = local1;                                                                                                            /* j-1 */
      else
        local_left = -1.0;
      if(j<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1)
        local_right = help->next->local;                                                                                        /* j+1 */
      else
        local_right = 12345677890.0;
      for(i=0; i<nPointList; i++)
      {
        Point = PointList[i];
        if(Compare_Points(&Point, global)==1)
        {
                                        #ifdef NO_PROJECT
          if(local_left>-SMALL)
            Line_Local2GlobalNew(theLine, global_left, local_left);
          else
            global_left[0] = global_left[1] = global_left[2] = -1e50;
          if(local_right<12345677890.0)
            Line_Local2GlobalNew(theLine, global_right, local_right);
          else
            global_right[0] = global_right[1] = global_right[2] = 1e50;
          LGM_BNDP_LINE_GLOBALLEFT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i])))[0] = global_left[0];
          LGM_BNDP_LINE_GLOBALLEFT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i])))[1] = global_left[1];
          LGM_BNDP_LINE_GLOBALLEFT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i])))[2] = global_left[2];
          LGM_BNDP_LINE_GLOBALRIGHT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i])))[0] = global_right[0];
          LGM_BNDP_LINE_GLOBALRIGHT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i])))[1] = global_right[1];
          LGM_BNDP_LINE_GLOBALRIGHT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i])))[2] = global_right[2];
                                        #else
          LGM_BNDP_LINE_LEFT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i]))) = local_left;
          LGM_BNDP_LINE_RIGHT(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i]))) = local_right;
                                        #endif
          LGM_BNDP_LINE(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i]))) = theLine;
          LGM_BNDP_NLINE(BNDP2LGM(theMesh->theBndPs[i]))++;
          break;
        }
      }
      local1 = help->local;
      help = help->next;
      assert(flag<2);
    }
  }

  /*	for(i=0;i<nPointList;i++)
          {
                  printf("%d ", i);
                  for(j=0;j<LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i]));j++)
                          printf("%d ", LGM_SURFACE_ID(LGM_BNDP_SURFACE(BNDP2LGM(theMesh->theBndPs[i]),j)));
                  printf("\n");
          }*/



  return (0);
}

static DOUBLE Calc_Line_Segment_Angle(LGM_LINE *theLine, INT i,  INT j,  INT k)
{
  DOUBLE l1, l2;
  DOUBLE s1[3], s2[3], scalarproduct;

  s1[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine, k))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine, j))[0];
  s1[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine, k))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine, j))[1];
  s1[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine, k))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine, j))[2];
  s2[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine, j))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine, i))[0];
  s2[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine, j))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine, i))[1];
  s2[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine, j))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine, i))[2];

  l1 = sqrt( s1[0]*s1[0] + s1[1]*s1[1] + s1[2]*s1[2] );
  l2 = sqrt( s2[0]*s2[0] + s2[1]*s2[1] + s2[2]*s2[2] );

  scalarproduct = ( s1[0]*s2[0] + s1[1]*s2[1] + s1[2]*s2[2] ) / (l1 * l2);
  return(scalarproduct);
}

static INT Point_Counter(LGM_LINE *theLine, DOUBLE h, INT StartIndex, INT EndIndex)
{
  int i,j,k,npoints;
  DOUBLE length_of_line,dist,dist1,dist2,dist3,slocal;
  DOUBLE lh, lh1, lh2, rest_length, in[5];

  if(h>0.0)
  {
    /* uniform dicretisation */
    length_of_line = 0;
    for(i=StartIndex; i<EndIndex-1; i++)
      length_of_line = length_of_line + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));

    /* calculate number of points on the line */
    npoints = (int)floor(length_of_line / h + 0.51) + 1;
    if(npoints<2)
      npoints = 2;
    /*		if(StartIndex!=0)
                            npoints--;*/
  }
  else
  {
    /* nonuniform dicretisation */
    length_of_line = 0;
    for(i=StartIndex; i<EndIndex-1; i++)
      length_of_line = length_of_line + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));
    rest_length = length_of_line;

    /* choose startpoint for linediscretization */
    /* the direction of discretisation is detremined by the information of
       the first and the last point in the line */
    in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[0];
    in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[1];
    in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[2];
    in[3] = h;
    (*LOCAL_H)(in, &lh1);
    in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[0];
    in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[1];
    in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[2];
    in[3] = h;
    (*LOCAL_H)(in, &lh2);

    if(lh1<=lh2)
    {
      /* start with first point */
      in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[0];
      in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[1];
      in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[2];
      in[3] = h;
      (*LOCAL_H)(in, &lh);

      if(StartIndex==0)
        npoints = 2;
      else
        npoints = 2;
      while(rest_length > 2 * lh)
      {
        rest_length = rest_length - lh;
        npoints++;
        dist = length_of_line - rest_length;
        dist1 = 0.0;
        for(j=StartIndex; j<EndIndex-1; j++)
        {
          dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if(dist1>dist)
            break;
        }

        dist2 = 0.0;
        for(k=StartIndex; k<j; k++)
          dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
        dist3 = dist - dist2;
        if(j==LGM_LINE_NPOINT(theLine)-1)
          slocal = 0.0;
        else
          slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if( (slocal<0.0) && (slocal>-SMALL) )
          slocal = 0.0;
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
        in[3] = h;
        (*LOCAL_H)(in, &lh);
      }
      if(rest_length > lh)
        npoints++;
    }
    else
    {
      /* start with last point */
      in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[0];
      in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[1];
      in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[2];
      in[3] = h;
      (*LOCAL_H)(in, &lh);

      if(StartIndex==0)
        npoints = 2;
      else
        npoints = 1;
      while(rest_length > 2 * lh)
      {
        rest_length = rest_length - lh;
        npoints++;
        dist = rest_length;
        dist1 = 0.0;
        for(j=StartIndex; j<EndIndex-1; j++)
        {
          dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if(dist1>dist)
            break;
        }

        dist2 = 0.0;
        for(k=StartIndex; k<j; k++)
          dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
        dist3 = dist - dist2;
        if(j==LGM_LINE_NPOINT(theLine)-1)
          slocal = 0.0;
        else
          slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if( (slocal<0.0) && (slocal>-SMALL) )
          slocal = 0.0;
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
        in[3] = h;
        (*LOCAL_H)(in, &lh);
      }
      if(rest_length > lh)
        npoints++;
    }
  }
  return(npoints);
}

static INT DiscretizeLineSegment(LGM_LINE *theLine, DOUBLE h, INT StartIndex, INT EndIndex, INT *Disc_nPoints)
{
  int i,j,k,npoints, sv, from_id, to_id;
  DOUBLE length_of_line,lenght_of_segment,dist,dist1,dist2,dist3,slocal;
  DOUBLE lh, lh1, lh2, rest_length, in[5], dummy;

  /* calculate number of points on the line */
  npoints = Point_Counter(theLine, h, StartIndex, EndIndex);

  from_id = (*Disc_nPoints);
  if(h>0.0)
  {
    /* uniform dicretisation */
    length_of_line = 0;
    for(i=StartIndex; i<EndIndex-1; i++)
      length_of_line = length_of_line + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));
    rest_length = length_of_line;

    lenght_of_segment = length_of_line / (npoints - 1);
    if(StartIndex==0)
      sv = 0;
    else
      sv = 1;
    for(i=sv; i<npoints; i++)
    {
      dist = i * lenght_of_segment;

      dist1 = 0.0;
      for(j=StartIndex; j<EndIndex-1; j++)
      {
        dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if(dist1>dist)
          break;
      }
      if(j>LGM_LINE_NPOINT(theLine)-1)
        j = LGM_LINE_NPOINT(theLine) - 1;

      dist2 = 0.0;
      for(k=StartIndex; k<j; k++)
        dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
      dist3 = dist - dist2;
      if(j==LGM_LINE_NPOINT(theLine)-1)
        slocal = 0.0;
      else
        slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
      if( (slocal<0.0) && (slocal>-SMALL) )
        slocal = 0.0;
      LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) = j + slocal;
      (*Disc_nPoints)++;
    }
  }
  else
  {
    /* nonuniform dicretisation */
    length_of_line = 0;
    for(i=StartIndex; i<EndIndex-1; i++)
      length_of_line = length_of_line + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));
    rest_length = length_of_line;

    /* choose startpoint for linediscretization */
    /* the direction of discretisation is detremined by the information of
       the first and the last point in the line */
    in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[0];
    in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[1];
    in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[2];
    in[3] = h;
    (*LOCAL_H)(in, &lh1);
    in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[0];
    in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[1];
    in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[2];
    in[3] = h;
    (*LOCAL_H)(in, &lh2);

    if(lh1<=lh2)
    {

      /* calculate the points -- start with first point */
      rest_length = length_of_line;
      in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[0];
      in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[1];
      in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,StartIndex))[2];
      in[3] = h;
      (*LOCAL_H)(in, &lh);
      if(StartIndex==0)
        npoints = 2;
      else
        npoints = 1;
      i = 0;
      LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) = (DOUBLE)StartIndex;
      (*Disc_nPoints)++;
      while(rest_length > 2 * lh)
      {
        rest_length = rest_length - lh;
        dist = length_of_line - rest_length;
        dist1 = 0.0;
        for(j=StartIndex; j<LGM_LINE_NPOINT(theLine)-1; j++)
        {
          dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if(dist1>dist)
            break;
        }
        if(j>LGM_LINE_NPOINT(theLine)-1)
          j = LGM_LINE_NPOINT(theLine) - 1;

        dist2 = 0.0;
        for(k=StartIndex; k<j; k++)
          dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
        dist3 = dist - dist2;
        if(j==LGM_LINE_NPOINT(theLine)-1)
          slocal = 0.0;
        else
          slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if( (slocal<0.0) && (slocal>-SMALL) )
          slocal = 0.0;
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
        in[3] = h;
        (*LOCAL_H)(in, &lh);
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) = j + slocal;
        (*Disc_nPoints)++;
        i++;
      }
      if(rest_length > lh)
      {
        lh = rest_length / 2;
        npoints++;
        rest_length = rest_length - lh;
        dist = length_of_line - rest_length;
        dist1 = 0.0;
        for(j=StartIndex; j<LGM_LINE_NPOINT(theLine)-1; j++)
        {
          dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if(dist1>dist)
            break;
        }
        if(j>LGM_LINE_NPOINT(theLine)-1)
          j = LGM_LINE_NPOINT(theLine) - 1;

        dist2 = 0.0;
        for(k=StartIndex; k<j; k++)
          dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
        dist3 = dist - dist2;
        if(j==LGM_LINE_NPOINT(theLine)-1)
          slocal = 0.0;
        else
          slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if( (slocal<0.0) && (slocal>-SMALL) )
          slocal = 0.0;
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
        in[3] = h;
        (*LOCAL_H)(in, &lh);
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) =
          j + slocal;
        (*Disc_nPoints)++;
        i++;
      }
      if(EndIndex==LGM_LINE_NPOINT(theLine))
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) =  LGM_LINE_NPOINT(theLine) - 1;
    }
    else
    {
      /* calculate the points -- start with last point */
      rest_length = length_of_line;
      in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[0];
      in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[1];
      in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,EndIndex-1))[2];
      in[3] = h;
      (*LOCAL_H)(in, &lh);
      npoints = 2;
      i = 0;
      LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) = (DOUBLE)StartIndex;
      (*Disc_nPoints)++;
      while(rest_length > 2 * lh)
      {
        rest_length = rest_length - lh;
        dist = rest_length;
        dist1 = 0.0;
        for(j=StartIndex; j<LGM_LINE_NPOINT(theLine)-1; j++)
        {
          dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if(dist1>dist)
            break;
        }

        if(j>LGM_LINE_NPOINT(theLine)-1)
          j = LGM_LINE_NPOINT(theLine) - 1;

        dist2 = 0.0;
        for(k=StartIndex; k<j; k++)
          dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
        dist3 = dist - dist2;
        if(j==LGM_LINE_NPOINT(theLine)-1)
          slocal = 0.0;
        else
          slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if( (slocal<0.0) && (slocal>-SMALL) )
          slocal = 0.0;
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
        in[3] = h;
        (*LOCAL_H)(in, &lh);
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) = j + slocal;
        (*Disc_nPoints)++;
      }
      if(rest_length > lh)
      {
        lh = rest_length / 2;
        npoints++;
        rest_length = rest_length - lh;
        dist = rest_length;
        dist1 = 0.0;
        for(j=StartIndex; j<LGM_LINE_NPOINT(theLine)-1; j++)
        {
          dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if(dist1>dist)
            break;
        }
        if(j>LGM_LINE_NPOINT(theLine)-1)
          j = LGM_LINE_NPOINT(theLine) - 1;

        dist2 = 0.0;
        for(k=StartIndex; k<j; k++)
          dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
        dist3 = dist - dist2;
        if(j==LGM_LINE_NPOINT(theLine)-1)
          slocal = 0.0;
        else
          slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if( (slocal<0.0) && (slocal>-SMALL) )
          slocal = 0.0;
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
        in[3] = h;
        (*LOCAL_H)(in, &lh);
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) = j + slocal;
        (*Disc_nPoints)++;
      }
      if(EndIndex==LGM_LINE_NPOINT(theLine))
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),*Disc_nPoints) =  LGM_LINE_NPOINT(theLine) - 1;

      /* change direction of points */
      to_id = from_id + (int)((*Disc_nPoints-from_id)/2);
      for(i=from_id+1; i<to_id+1; i++)
      {
        dummy = LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),i);
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),i) =
          LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),(*Disc_nPoints) - i + from_id);
        LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),(*Disc_nPoints) - i + from_id)
          = dummy;
      }
    }
  }
  return (0);
}

static INT DiscretizeLine (HEAP *Heap, LGM_LINE *theLine, DOUBLE h, LGM_POINT *pointlist, INT norp, INT MarkKey)
{
  INT i, j, np, npoints, folds, StartIndex, EndIndex, Disc_nPoints, resolve;
  DOUBLE scalarprodukt;
  LGM_POINT p;

  np = 0;
  npoints = 0;
  folds = 0;

  /* count the new points on the line for allocate Mem */
  if(LGM_LINE_NPOINT(theLine) == 2)
    npoints = Point_Counter(theLine, h, 0, 2);
  else
  {
    StartIndex = 0;
    for(i=1; i<LGM_LINE_NPOINT(theLine)-1; i++)
    {
      resolve = 0;
      for(j=0; j<norp; j++)
      {
        p.position[0] = pointlist[j].position[0];
        p.position[1] = pointlist[j].position[1];
        p.position[2] = pointlist[j].position[2];
        if(POINT_DIST(LGM_LINE_POINT(theLine, i), &p)==1)
          resolve = 1;
      }
      scalarprodukt = Calc_Line_Segment_Angle(theLine, i-1, i, i+1);
      if( (scalarprodukt<cosAngle)||(resolve==1) )
      {
        folds++;
        EndIndex = i;
        npoints = npoints + Point_Counter(theLine, h, StartIndex, EndIndex+1);
        if(StartIndex!=0)
          npoints--;
        StartIndex = EndIndex;
      }
    }
    npoints = npoints + Point_Counter(theLine, h, StartIndex, LGM_LINE_NPOINT(theLine)) /*- folds + 1*/;
    if(StartIndex!=0)
      npoints--;
  }
  if(LGM_DEBUG)
    printf("%s %d\n", "folds: ", (int)folds);
  /* store npoints in the line-structure */
  LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine)) = npoints;

  /* allocate memery for local coordinates */
  LGM_LINE_LINEDISC(theLine)->local = (DOUBLE*)GetTmpMem(Heap,sizeof(DOUBLE)*(npoints+5),LGM_MarkKey);

  Disc_nPoints = 0;
  /* generate the new points on the line */
  if(LGM_LINE_NPOINT(theLine) == 2)
    DiscretizeLineSegment(theLine, h, 0, 2, &Disc_nPoints);
  else
  {
    StartIndex = 0;
    for(i=1; i<LGM_LINE_NPOINT(theLine)-1; i++)
    {
      resolve = 0;
      for(j=0; j<norp; j++)
      {
        p.position[0] = pointlist[j].position[0];
        p.position[1] = pointlist[j].position[1];
        p.position[2] = pointlist[j].position[2];
        if(POINT_DIST(LGM_LINE_POINT(theLine, i), &p)==1)
          resolve = 1;
      }
      scalarprodukt = Calc_Line_Segment_Angle(theLine, i-1, i, i+1);
      if( (scalarprodukt<cosAngle)||(resolve==1) )
      {
        EndIndex = i;
        DiscretizeLineSegment(theLine, h, StartIndex, EndIndex+1, &Disc_nPoints);
        StartIndex = EndIndex;
      }
    }
    DiscretizeLineSegment(theLine, h, StartIndex, LGM_LINE_NPOINT(theLine), &Disc_nPoints);
  }
  /* fuer den Fall, dass die Surface nur von einer Line berandet wird */
  /*	if(LGM_LINE_BEGIN(theLine)==LGM_LINE_END(theLine))
          {
                  LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))--;
                  UserWriteF(" Line %4d: %4d Linesegments\n",
                                             LGM_LINE_ID(theLine),
                                             LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine)));

          }
          else*/
  UserWriteF(" Line %4d: %4d Linesegments\n",
             LGM_LINE_ID(theLine),
             LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1);
  return(0);
}

static INT AddPoint2Line(HEAP *Heap, LGM_LINE *theLine, DOUBLE local, INT k, INT MarkKey)
{
  LINEPOINT *help, *help1;

  if(k==0)
  {
    help1 = (LINEPOINT*)GetTmpMem(Heap,sizeof(LINEPOINT), MarkKey);
    if(help1==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
    help1->local = local;
    help1->next = NULL;
    LGM_LINEDISCNEW_POINT(LGM_LINE_LINEDISCNEW(theLine)) = help1;
    LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine)) = LGM_LINEDISCNEW_POINT(LGM_LINE_LINEDISCNEW(theLine));
    LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine)) = 1;
  }
  else
  {
    help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
    while((help->next!=NULL)&&(help->next)->local<local)
      help = help->next;

    help1 = (LINEPOINT*)GetTmpMem(Heap,sizeof(LINEPOINT), MarkKey);
    if(help1==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
    help1->local = local;
    help1->next = help->next;
    help->next = help1;
    LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))++;
  }

  return(0);
}

static DOUBLE LinePointDistance(LGM_LINE *theLine, INT Index)
{
  int i, start_i, end_i;
  DOUBLE length_of_segment, start_local, end_local, start_s, end_s;
  LINEPOINT *help;

  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
  for(i=0; i<Index; i++)
    help = help->next;

  start_local = help->local;
  end_local = (help->next)->local;

  start_i = (int)floor(start_local);
  start_s = start_local - start_i;
  end_i = (int)floor(end_local);
  end_s = end_local - end_i;

  length_of_segment = 0;

  for(i=start_i+1; i<end_i; i++)
    length_of_segment = length_of_segment + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));

  if(end_i<LGM_LINE_NPOINT(theLine)-1)
    length_of_segment = length_of_segment +
                        end_s * LGM_POINT_DIST(LGM_LINE_POINT(theLine,end_i),LGM_LINE_POINT(theLine,end_i+1));
  if(end_i>start_i)
    length_of_segment = length_of_segment +
                        (1.0-start_s) * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));
  else
    length_of_segment = length_of_segment -
                        start_s * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));

  return(length_of_segment);
}

static INT Check_Line_Linedisc_Distance(HEAP *Heap, LGM_LINE *theLine, INT Index, DOUBLE h, INT MarkKey)
{
  int i, start_i, end_i;
  DOUBLE length_of_segment, start_local, end_local, start_s, end_s;
  DOUBLE start_p[3], end_p[3], mid_p[3], len, line_p[3], mid_i, mid_local, local, dist, abs_h;
  LINEPOINT *help;

  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
  for(i=0; i<Index; i++)
    help = help->next;

  start_local = help->local;
  end_local = (help->next)->local;

  start_i = (int)floor(start_local);
  start_s = start_local - start_i;
  end_i = (int)floor(end_local);
  end_s = end_local - end_i;
  mid_i = -1;
  mid_local = 0.0;

  length_of_segment = LinePointDistance(theLine, Index);

  Line_Local2GlobalNew(theLine,start_p,start_local);
  Line_Local2GlobalNew(theLine,end_p,end_local);
  for(i=0; i<3; i++)
    mid_p[i] = 0.5 * (start_p[i] + end_p[i]);

  dist = length_of_segment / 2;

  len = 0;
  if(start_i==end_i)
  {
    len = (end_s - start_s) * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));
    mid_i = start_i;
    mid_local = (end_s + start_s) / 2;
  }
  else
  {
    len = (1.0 - start_s) * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));
    if(len>dist)
    {
      mid_i = start_i;
      /*			mid_local = (1.0 + start_s) / 2;*/
      /*			mid_local = dist/len;*/
      mid_local = dist/len * (1-start_s) + start_s;
    }
    else
    {
      for(i=start_i+1; i<end_i; i++)
        if(len + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1)) < dist)
          len = len + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));
        else
        {
          mid_i = i;
          mid_local = (dist - len) / LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));
          break;
        }
      if(mid_i==-1)
      {
        mid_i = end_i;
        mid_local = (dist - len) / LGM_POINT_DIST(LGM_LINE_POINT(theLine,end_i),LGM_LINE_POINT(theLine,end_i+1));
      }
    }
  }

  local = (DOUBLE)mid_i + mid_local;
  Line_Local2GlobalNew(theLine,line_p,local);
  if(h<0.0)
    abs_h = -h;
  else
    abs_h = h;
  if(sqrt( (line_p[0]-mid_p[0])*(line_p[0]-mid_p[0]) +
           (line_p[1]-mid_p[1])*(line_p[1]-mid_p[1]) +
           (line_p[2]-mid_p[2])*(line_p[2]-mid_p[2]) ) < LINE_DISTANCE * abs_h)
    return(0);
  else
  {
    local = (DOUBLE)mid_i + mid_local;
    AddPoint2Line(Heap, theLine, local, 1, MarkKey);
    if(LGM_DEBUG)
      printf("%s\n", "Line refined");
    return(1);
  }
}

static INT RefineLine(HEAP *Heap, LGM_LINE *theLine, INT Index, DOUBLE h, INT MarkKey)
{
  int i, j, k, start_i, end_i, npoints;
  DOUBLE length_of_segment, length_of_linesegment, start_local, end_local, start_s, end_s, slocal, local;
  DOUBLE dist, dist1, dist2, dist3, rest_length, in[4], lh, lh1, lh2;
  LINEPOINT *help;

  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
  for(i=0; i<Index; i++)
    help = help->next;

  start_local = help->local;
  end_local = (help->next)->local;

  start_i = (int)floor(start_local);
  start_s = start_local - start_i;
  end_i = (int)floor(end_local);
  end_s = end_local - end_i;

  length_of_linesegment = 0;

  for(i=start_i; i<end_i; i++)
    length_of_linesegment = length_of_linesegment + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));

  if(length_of_linesegment>h)                           /* refine line */
  {
    if(h>0.0)
    {
      npoints = (int)ceil(length_of_linesegment / h + 0.51);
      /* uniform dicretisation */
      length_of_segment = length_of_linesegment / npoints;

      for(i=1; i<npoints; i++)
      {
        dist = i * length_of_segment;

        dist1 = 0.0;
        for(j=start_i; j<end_i; j++)
        {
          dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if(dist1>dist)
            break;
        }
        dist2 = 0.0;
        for(k=start_i; k<j; k++)
          dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
        dist3 = dist - dist2;
        slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
        if( (slocal<0.0) && (slocal>-SMALL) )
          slocal = 0.0;

        local = (DOUBLE)j + slocal;
        AddPoint2Line(Heap, theLine, local, start_i+i, MarkKey);
      }
    }
    else
    {
      /* nonuniform dicretisation */
      rest_length = length_of_linesegment;

      /* choose startpoint for linediscretization */
      /* the direction of discretisation is detremined by the information of
              the first and the last point in the line */
      in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,start_i))[0];
      in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,start_i))[1];
      in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,start_i))[2];
      in[3] = h;
      (*LOCAL_H)(in, &lh1);
      in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,end_i))[0];
      in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,end_i))[1];
      in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,end_i))[2];
      in[3] = h;
      (*LOCAL_H)(in, &lh2);

      if(lh1<=lh2)
      {
        /* calculate the points -- start with first point */
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,start_i))[0];
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,start_i))[1];
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,start_i))[2];
        in[3] = h;
        (*LOCAL_H)(in, &lh);
        i = 0;
        while(rest_length > 2 * lh)
        {
          rest_length = rest_length - lh;
          dist = length_of_linesegment - rest_length;
          dist1 = 0.0;
          for(j=start_i; j<LGM_LINE_NPOINT(theLine)-1; j++)
          {
            dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
            if(dist1>dist)
              break;
          }
          if(j>LGM_LINE_NPOINT(theLine)-1)
            j = LGM_LINE_NPOINT(theLine) - 1;

          dist2 = 0.0;
          for(k=start_i; k<j; k++)
            dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
          dist3 = dist - dist2;
          if(j==LGM_LINE_NPOINT(theLine)-1)
            slocal = 0.0;
          else
            slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if( (slocal<0.0) && (slocal>-SMALL) )
            slocal = 0.0;
          in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
          in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
          in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
          in[3] = h;
          (*LOCAL_H)(in, &lh);
          local = (DOUBLE)j + slocal;
          AddPoint2Line(Heap, theLine, local, 1, MarkKey);
        }
        if(rest_length > lh)
        {
          lh = rest_length / 2;
          npoints++;
          rest_length = rest_length - lh;
          dist = length_of_linesegment - rest_length;
          dist1 = 0.0;
          for(j=start_i; j<LGM_LINE_NPOINT(theLine)-1; j++)
          {
            dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
            if(dist1>dist)
              break;
          }
          if(j>LGM_LINE_NPOINT(theLine)-1)
            j = LGM_LINE_NPOINT(theLine) - 1;

          dist2 = 0.0;
          for(k=start_i; k<j; k++)
            dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
          dist3 = dist - dist2;
          if(j==LGM_LINE_NPOINT(theLine)-1)
            slocal = 0.0;
          else
            slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if( (slocal<0.0) && (slocal>-SMALL) )
            slocal = 0.0;
          in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
          in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
          in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
          in[3] = h;
          (*LOCAL_H)(in, &lh);
          local = (DOUBLE)j + slocal;
          AddPoint2Line(Heap, theLine, local, 1, MarkKey);
        }
      }
      else
      {
        /* calculate the points -- start with last point */
        rest_length = length_of_linesegment;
        in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,end_i))[0];
        in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,end_i))[1];
        in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,end_i))[2];
        in[3] = h;
        (*LOCAL_H)(in, &lh);
        npoints = 2;
        while(rest_length > 2 * lh)
        {
          rest_length = rest_length - lh;
          dist = rest_length;
          dist1 = 0.0;
          for(j=start_i; j<LGM_LINE_NPOINT(theLine)-1; j++)
          {
            dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
            if(dist1>dist)
              break;
          }

          if(j>LGM_LINE_NPOINT(theLine)-1)
            j = LGM_LINE_NPOINT(theLine) - 1;

          dist2 = 0.0;
          for(k=start_i; k<j; k++)
            dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
          dist3 = dist - dist2;
          if(j==LGM_LINE_NPOINT(theLine)-1)
            slocal = 0.0;
          else
            slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if( (slocal<0.0) && (slocal>-SMALL) )
            slocal = 0.0;
          in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
          in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
          in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
          in[3] = h;
          (*LOCAL_H)(in, &lh);
          local = (DOUBLE)j + slocal;
          AddPoint2Line(Heap, theLine, local, 1, MarkKey);
        }
        if(rest_length > lh)
        {
          lh = rest_length / 2;
          npoints++;
          rest_length = rest_length - lh;
          dist = rest_length;
          dist1 = 0.0;
          for(j=start_i; j<LGM_LINE_NPOINT(theLine)-1; j++)
          {
            dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
            if(dist1>dist)
              break;
          }
          if(j>LGM_LINE_NPOINT(theLine)-1)
            j = LGM_LINE_NPOINT(theLine) - 1;

          dist2 = 0.0;
          for(k=start_i; k<j; k++)
            dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
          dist3 = dist - dist2;
          if(j==LGM_LINE_NPOINT(theLine)-1)
            slocal = 0.0;
          else
            slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
          if( (slocal<0.0) && (slocal>-SMALL) )
            slocal = 0.0;
          in[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[0] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0]);
          in[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[1] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1]);
          in[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]
                  + slocal * (LGM_POINT_POS(LGM_LINE_POINT(theLine,j+1))[2] - LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[2]);
          in[3] = h;
          (*LOCAL_H)(in, &lh);
          local = (DOUBLE)j + slocal;
          AddPoint2Line(Heap, theLine, local, 1, MarkKey);
        }
      }
    }
  }
  return(0);
}

static INT DiscretizeLineNew (HEAP *Heap, LGM_LINE *theLine, DOUBLE h, LGM_POINT *pointlist, INT norp, INT MarkKey)
{
  INT i, j, resolve, flag;
  DOUBLE scalarprodukt, local, len, in[4], lh, lh1, lh2;
  LGM_POINT p;
  LINEPOINT *help;

  LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine)) = 0;

  /* resolve start-, end- and essential points in linedisc */

  /* startpoint of the line */
  local = 0.0;
  AddPoint2Line(Heap, theLine, local, 0, MarkKey);

  /* endpoint of the line */
  local = (DOUBLE)LGM_LINE_NPOINT(theLine) - 1.0;
  AddPoint2Line(Heap, theLine, local, 1, MarkKey);

  /* essential points of the line */
  for(i=1; i<LGM_LINE_NPOINT(theLine)-1; i++)
  {
    resolve = 0;
    for(j=0; j<norp; j++)
    {
      p.position[0] = pointlist[j].position[0];
      p.position[1] = pointlist[j].position[1];
      p.position[2] = pointlist[j].position[2];
      if(POINT_DIST(LGM_LINE_POINT(theLine, i), &p)==1)
        resolve = 1;
    }
    scalarprodukt = Calc_Line_Segment_Angle(theLine, i-1, i, i+1);
    if( (scalarprodukt<cosAngle)||(resolve==1) )
    {
      local = (DOUBLE)i;
      AddPoint2Line(Heap, theLine, local, LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1, MarkKey);
    }
  }

  /* refine line */
  do
  {
    for(i=0; i<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1; i++)
    {
      len = LinePointDistance(theLine, i);
      if(h>0.0)
        lh = h;
      else
      {
        help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
        for(j=0; j<i; j++)
          help = help->next;
        Line_Local2GlobalNew(theLine,in,help->local);
        in[3] = h;
        (*LOCAL_H)(in, &lh1);
        help = help->next;
        Line_Local2GlobalNew(theLine,in,help->local);
        in[3] = h;
        (*LOCAL_H)(in, &lh2);
        if(lh1<lh2)
          lh = lh1;
        else
          lh = lh2;

      }
      if(len>(1.0+EPS)*lh)
      {
        RefineLine(Heap, theLine, i, h, MarkKey);
        break;
      }
    }
  }
  while(i<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1);

  /* check, if linedisc aproximates the line */
  do
  {
    for(i=0; i<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1; i++)
    {
      flag = Check_Line_Linedisc_Distance(Heap, theLine, i, h, MarkKey);
      if(flag)
        break;
    }
  }
  while(flag);

  UserWriteF(" Line %4d: %4d Linesegments\n",
             LGM_LINE_ID(theLine),
             LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1);
  return(0);
}

static INT TransferSurfaces2Mesh (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h)
{
  INT i,id;
  INT *newId;

  /* find id's of the trianglecorners */
  newId = (INT*)GetTmpMem(Heap,(LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))+1)*sizeof(INT),LGM_MarkKey);
  if(newId==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); i++)
    newId[i] = - 1;

  for (i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)); i++)
  {
    id = LGM_SURFACE_LEFT(theSurface);
    if (id!=0)
    {
      theMesh->Side_corners[id][theMesh->nSides[id]] = 3;
      if(SURFACE_ADAPT)
        theMesh->xy_Side[id][theMesh->nSides[id]] = LGM_SURFDISC_XY_ID(LGM_SURFACE_DISC(theSurface), i);
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] =
        LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,0));
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] =
        LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,2));
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][2] =
        LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,1));
      theMesh->nSides[id]++;
    }
    id = LGM_SURFACE_RIGHT(theSurface);
    if (id!=0)
    {
      theMesh->Side_corners[id][theMesh->nSides[id]] = 3;
      if(SURFACE_ADAPT)
        theMesh->xy_Side[id][theMesh->nSides[id]] = LGM_SURFDISC_XY_ID(LGM_SURFACE_DISC(theSurface), i);
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] =
        LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,0));
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] =
        LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,1));
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][2] =
        LGM_SURFDISC_MESH_ID(LGM_SURFACE_DISC(theSurface), LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,2));
      theMesh->nSides[id]++;
    }
  }
  return (0);
}

#ifdef _NETGEN
static INT Get_Folds_Surface(LGM_SURFACE *theSurface, LGM_LINE **lineptr)
{
  INT i, j, folds, a;
  DOUBLE scalarproduct;

  folds = 0;

  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    for(j=0; j<3; j++)
    {
      if(LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j)!=-1)
      {
        scalarproduct = Calc_Triangle_Angle(theSurface, i, LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j));
        if(scalarproduct < cosAngle)
        {
          LGM_LINE_ID((lineptr)[folds]) = -1;
          LGM_LINE_NPOINT((lineptr)[folds]) = 2;
          LGM_LINE_BEGIN((lineptr)[folds]) = (j+1)%3;
          LGM_LINE_END((lineptr)[folds]) = (j+2)%3;

          LGM_POINT_POS(LGM_LINE_POINT((lineptr)[folds], 0))[0] =
            LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3))[0];
          LGM_POINT_POS(LGM_LINE_POINT((lineptr)[folds], 0))[1] =
            LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3))[1];
          LGM_POINT_POS(LGM_LINE_POINT((lineptr)[folds], 0))[2] =
            LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+1)%3))[2];

          LGM_POINT_POS(LGM_LINE_POINT((lineptr)[folds], 1))[0] =
            LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+2)%3))[0];
          LGM_POINT_POS(LGM_LINE_POINT((lineptr)[folds], 1))[1] =
            LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+2)%3))[1];
          LGM_POINT_POS(LGM_LINE_POINT((lineptr)[folds], 1))[2] =
            LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(j+2)%3))[2];
          folds++;
        }
      }
    }
  }
  return(folds);
}

static INT Get_Id(LGM_LINE *theLine, LGM_SURFACE *theSurface, INT index1, INT index2,
                  DOUBLE **line_points, INT *nline_points, DOUBLE **lcoord_lpoints)
{
  INT i, j, k,  flag1, flag2, start_id, end_id;
  DOUBLE global[3], startpoint[3], endpoint[3], local, n[3];
  LINEPOINT *help;

  flag1 = 0;
  flag2 = 0;

  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
  for(i=0; i<index1; i++)
    help = help->next;

  Line_Local2GlobalNew(theLine,startpoint,help->local);

  help = LGM_LINEDISCNEW_START(LGM_LINE_LINEDISCNEW(theLine));
  for(i=0; i<index2; i++)
    help = help->next;

  Line_Local2GlobalNew(theLine,endpoint,help->local);

  start_id = *nline_points;
  for(k=0; k<*nline_points; k++)
  {
    if( P_Dist(startpoint, line_points[k]) )
    {
      flag1 = 1;
      start_id = k;
    }
  }
  if(flag1==0)
  {
    AddLinePoint(start_id+1, startpoint[0],startpoint[1],startpoint[2]);
    for(i=0; i<3; i++)
      line_points[*nline_points][i] = startpoint[i];
    Line_Global2Local(theLine, startpoint, &local);
    lcoord_lpoints[*nline_points][0] = - (DOUBLE) LGM_LINE_ID(theLine) - 2.0;
    lcoord_lpoints[*nline_points][1] = local;
    (*nline_points)++;
  }

  end_id = *nline_points;
  for(k=0; k<*nline_points; k++)
  {
    for(i=0; i<3; i++)
      global[i] = line_points[k][i];
    if( P_Dist(endpoint, line_points[k]) )
    {
      flag2 = 1;
      end_id = k;
    }
  }
  if(flag2==0)
  {
    AddLinePoint(end_id+1, endpoint[0],endpoint[1],endpoint[2]);
    for(i=0; i<3; i++)
      line_points[*nline_points][i] = endpoint[i];
    Line_Global2Local(theLine, endpoint, &local);
    lcoord_lpoints[*nline_points][0] = - (DOUBLE) LGM_LINE_ID(theLine) - 2.0;
    lcoord_lpoints[*nline_points][1] = local;
    (*nline_points)++;
  }

  AddLineSegment(start_id+1, end_id+1);

  return(0);
}

static INT TransferLine2Surface(LGM_LINE *theLine,
                                LGM_SURFACE *theSurface,
                                INT direction,
                                DOUBLE **line_points,
                                INT *nline_points,
                                DOUBLE **lcoord_lpoints)
{
  INT i, j, k, flag1, flag2, start_id, end_id;
  DOUBLE global[3], startpoint[3], endpoint[3], local[2];

  LGM_LINE_USED(theLine) = 1;
  if(direction==1)
  {
    for(j=0; j<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1; j++)
      Get_Id(theLine, theSurface, j, j+1, line_points, nline_points, lcoord_lpoints);
  }
  if(direction==-1)
  {
    for(j=LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1; j>0; j--)
      Get_Id(theLine, theSurface, j, j-1, line_points, nline_points, lcoord_lpoints);
  }
  if(direction==-2)
  {
    for(j=0; j<LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(theLine))-1; j++)
      Get_Id(theLine, theSurface, j, j+1, line_points, nline_points, lcoord_lpoints);
  }
  return(0);
}

static INT Get_Direction(LGM_SURFACE *theSurface, LGM_LINE *theLine)
{
  DOUBLE x1, x2, y1, y2, z1, z2;
  INT i, j, k, direction;

  /* find direction for this line */
  direction = 0;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    for(k=0; k<3; k++)
    {
      x1 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+0)%3))[0]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[0];
      y1 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+0)%3))[1]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[1];
      z1 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+0)%3))[2]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[2];
      x2 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+1)%3))[0]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,1))[0];
      y2 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+1)%3))[1]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,1))[1];
      z2 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+1)%3))[2]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,1))[2];

      if(    (x1<SMALL) && (x1>-SMALL)
             && (y1<SMALL) && (y1>-SMALL)
             && (z1<SMALL) && (z1>-SMALL)
             && (x2<SMALL) && (x2>-SMALL)
             && (y2<SMALL) && (y2>-SMALL)
             && (z2<SMALL) && (z2>-SMALL) )
        direction = 1;
    }
  }

  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    for(k=0; k<3; k++)
    {
      x1 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+0)%3))[0]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,1))[0];
      y1 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+0)%3))[1]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,1))[1];
      z1 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+0)%3))[2]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,1))[2];
      x2 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+1)%3))[0]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[0];
      y2 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+1)%3))[1]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[1];
      z2 = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),(k+1)%3))[2]
           -         LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[2];

      if(    (x1<SMALL) && (x1>-SMALL)
             && (y1<SMALL) && (y1>-SMALL)
             && (z1<SMALL) && (z1>-SMALL)
             && (x2<SMALL) && (x2>-SMALL)
             && (y2<SMALL) && (y2>-SMALL)
             && (z2<SMALL) && (z2>-SMALL) )
        direction = -1;
    }
  }
  return(direction);
}

static INT DiscretizeSurface (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h, LGM_POINT *pointlist, INT norp,INT D, INT MarkKey)
{
  INT i,n,ni,j,k,offset,nils,id,ls_offset, mi, dummy;
  LGM_BNDP *theBndPList, *theBndP;
  DOUBLE global[3], globalnew[3],globalbndp[3],local[2], x1, y1, z1, x2, y2, z2;
  INT ii,oldline,newline,direction,a,b;
  INT nnew[5],*pi;
  INT *oldId,nsp,aksp;
  INT *newId,nId[3], neighbor[3],oldnb,newpoints,oldpoints,size, folds;
  LGM_POINT **ptrlst;
  LGM_LINE **lineptr, *theLine;
  DOUBLE startpoint[3], endpoint[3];
  INT start_id, end_id, flag, all_used, nline_points, nlp;
  DOUBLE **line_points, **lcoord_lpoints;

  surf_err = 0;
  ptrlst = (LGM_POINT**)GetTmpMem(Heap,(LGM_SURFACE_NPOINT(theSurface)+1)*sizeof(LGM_POINT*),MarkKey);
  if(ptrlst==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  for(j=0; j<LGM_SURFACE_NPOINT(theSurface); j++)
    ptrlst[j] = LGM_SURFACE_POINT(theSurface,j);

  /* transfer data to gg */
  /* points from inputdata */
  for(i=0; i<LGM_SURFACE_NPOINT(theSurface); i++)
  {
    global[0] = LGM_POINT_POS(LGM_SURFACE_POINT(theSurface,i))[0];
    global[1] = LGM_POINT_POS(LGM_SURFACE_POINT(theSurface,i))[1];
    global[2] = LGM_POINT_POS(LGM_SURFACE_POINT(theSurface,i))[2];
    AddGeomPoint(i,global[0],global[1],global[2]);
  }

  /* triangles from input data */
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    for(j=0; j<LGM_SURFACE_NPOINT(theSurface); j++)
      for(k=0; k<3; k++)
        if(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),k)==ptrlst[j])
          nId[k] = j;
    for(k=0; k<3; k++)
      LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),k) = nId[k];
  }

  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    for(j=0; j<3; j++)
      if(LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j)!=-1)
        neighbor[j] = LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j) + 1;
      else
        neighbor[j] = -1;
    AddGeomElement( LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),0)+1,
                    LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),1)+1,
                    LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),2)+1,
                    neighbor[0], neighbor[1], neighbor[2]);
  }

  /* set used-flags to 0 */
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    LGM_LINE_USED(LGM_SURFACE_LINE(theSurface, i)) = 0;

  folds = 0;
  /* search folds on the Surface */
  folds = Count_Folds_Surface(theSurface);
  if(LGM_DEBUG)
    printf("%s %d\n", "FOLDS: ", folds);

  /* allocate lines */
  if(folds>0)
    if ((lineptr=(LGM_LINE**)GetTmpMem(Heap,sizeof(LGM_LINE*)*(folds+1),LGM_MarkKey)) == NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }

  for (i=0; i<folds; i++)
  {
    size = sizeof(LGM_LINE) + 2*sizeof(LGM_POINT);
    if ((lineptr[i] = (LGM_LINE*)GetTmpMem(Heap,size,LGM_MarkKey)) == NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
  }

  Get_Folds_Surface(theSurface, lineptr);

  for (i=0; i<folds; i++)
  {
    lineptr[i]->ldisc = (LGM_LINEDISC*)GetTmpMem(Heap,sizeof(LGM_LINEDISC),LGM_MarkKey);
    if( (lineptr[i]->ldisc==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }
    lineptr[i]->ldiscnew = (LGM_LINEDISCNEW*)GetTmpMem(Heap,sizeof(LGM_LINEDISCNEW),LGM_MarkKey);
    if( (lineptr[i]->ldiscnew==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }
    if (DiscretizeLine(Heap,lineptr[i],h, pointlist, norp, MarkKey))
      return(NULL);
    if (DiscretizeLineNew(Heap,lineptr[i],h, pointlist, norp, MarkKey))
      return(NULL);
    if(LGM_DEBUG)
      PrintLineInfo(lineptr[i]);
  }

  nsp = 0;
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    nsp = nsp + LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(LGM_SURFACE_LINE(theSurface,i)))-1;
  for (i=0; i<folds; i++)
    nsp = nsp + LGM_LINEDISCNEW_NPOINT(LGM_LINE_LINEDISCNEW(lineptr[i]))-1;

  line_points = (DOUBLE **) GetTmpMem(Heap,(nsp+1)*sizeof(DOUBLE*),LGM_MarkKey);
  if(line_points==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  for(i=0; i<nsp+1; i++)
  {
    line_points[i] = (DOUBLE *) GetTmpMem(Heap,3*sizeof(DOUBLE),LGM_MarkKey);
    if(line_points[i]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
  }
  nline_points = 0;
  lcoord_lpoints = (DOUBLE **) GetTmpMem(Heap,(nsp+1)*sizeof(DOUBLE*),LGM_MarkKey);
  if(lcoord_lpoints==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  for(i=0; i<nsp+1; i++)
  {
    lcoord_lpoints[i] = (DOUBLE *) GetTmpMem(Heap,2*sizeof(DOUBLE),LGM_MarkKey);
    if(lcoord_lpoints[i]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
  }
  InitSurface(LOCAL_H);

  LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)) = 0;
  LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)) = 0;

  /* input of the linediscretisation */

  do
  {
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      if(LGM_LINE_USED(LGM_SURFACE_LINE(theSurface, i)) == 0)
        break;
    theLine = LGM_SURFACE_LINE(theSurface, i);
    oldline = i;

    direction = Get_Direction(theSurface, theLine);

    if(direction==0)
      return(NULL);

    TransferLine2Surface(theLine, theSurface, direction, line_points, &nline_points, lcoord_lpoints);

    for(i=1; i<LGM_SURFACE_NLINE(theSurface); i++)
    {
      /* search next line */
      for(k=1; k<LGM_SURFACE_NLINE(theSurface); k++)
      {
        if((k!=oldline))
        {
          if(direction==1)
          {
            if (LGM_LINE_END(LGM_SURFACE_LINE(theSurface,oldline))
                ==LGM_LINE_BEGIN(LGM_SURFACE_LINE(theSurface,k)))
            {
              newline = k;
              direction = 1;
              theLine = LGM_SURFACE_LINE(theSurface,newline);
              if(LGM_LINE_USED(theLine)==1)
                break;
              TransferLine2Surface(theLine, theSurface, direction, line_points, &nline_points, lcoord_lpoints);
              oldline = newline;
              break;
            }
            if (LGM_LINE_END(LGM_SURFACE_LINE(theSurface,oldline))
                ==LGM_LINE_END(LGM_SURFACE_LINE(theSurface,k)))
            {
              newline = k;
              direction = -1;
              theLine = LGM_SURFACE_LINE(theSurface,newline);
              if(LGM_LINE_USED(theLine)==1)
                break;
              TransferLine2Surface(theLine, theSurface, direction, line_points, &nline_points, lcoord_lpoints);
              oldline = newline;
              break;
            }
          }
          if(direction==-1)
          {
            if (LGM_LINE_BEGIN(LGM_SURFACE_LINE(theSurface,oldline))
                ==LGM_LINE_BEGIN(LGM_SURFACE_LINE(theSurface,k)))
            {
              newline = k;
              direction = 1;
              theLine = LGM_SURFACE_LINE(theSurface,newline);
              if(LGM_LINE_USED(theLine)==1)
                break;
              TransferLine2Surface(theLine, theSurface, direction, line_points, &nline_points, lcoord_lpoints);
              oldline = newline;
              break;
            }
            if (LGM_LINE_BEGIN(LGM_SURFACE_LINE(theSurface,oldline))
                ==LGM_LINE_END(LGM_SURFACE_LINE(theSurface,k)))
            {
              newline = k;
              direction = -1;
              theLine = LGM_SURFACE_LINE(theSurface,newline);
              if(LGM_LINE_USED(theLine)==1)
                break;
              TransferLine2Surface(theLine, theSurface, direction, line_points, &nline_points, lcoord_lpoints);
              oldline = newline;
              break;
            }
          }
        }
      }
    }
    all_used = 1;
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      if(LGM_LINE_USED(LGM_SURFACE_LINE(theSurface,i))==0)
        all_used = 0;                                                   /* noch nicht alle Lines uebergeben */
  }
  while(all_used==0);

  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    for(j=0; j<3; j++)
      if(LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j)!=-1)
        neighbor[j] = LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(theSurface,i),j) + 1;
      else
        neighbor[j] = -1;
  }

  for (i=0; i<folds; i++)
  {
    theLine = lineptr[i];
    TransferLine2Surface(theLine, theSurface, -2, line_points, &nline_points, lcoord_lpoints);
  }

  /* check used-flags */
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    assert(LGM_LINE_USED(LGM_SURFACE_LINE(theSurface, i))==1);

  oldnb = LGM_SURFACE_NPOINT(theSurface);

  if(SURFACE_ADAPT)
    if(LGM_SURFACE_ID(theSurface)!=-1)
      if(Surface_Plane(theSurface))
      {
        if(LGM_DEBUG)
          printf("%s %d %s\n", "Surface", LGM_SURFACE_ID(theSurface), "is plain");
        nlp = nline_points + 1;
        id = 0;
        for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(xy_Surface)); i++)
        {
          xy_flag[i] = 0;
          xy_new_id[i] = -1;
        }
        for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(xy_Surface)); i++)
          xy_element_id[i] = 0;
        for(i=xy_bpoints; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(xy_Surface)); i++)
        {
          local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(xy_Surface),i,0);
          local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(xy_Surface),i,1);
          Surface_Local2Global(xy_Surface, global, local);
          mi = Project2Surface(theSurface, global, local, xy_n);
          /*mi = GetLocalKoord(theSurface, global, local, xy_n);*/
          if(mi!=-1)
          {
            Surface_Local2Global(theSurface, global, local);
            flag = 0;
            for(j=0; j<nline_points; j++)
            {
              /*local[0] = lcoord_lpoints[j][0];
                 local[1] = lcoord_lpoints[j][1];
                 Surface_Local2Global(theSurface, globalnew, local);*/
              for(k=0; k<3; k++)
                globalnew[k] = line_points[j][k];
              if(E_Distance(global, globalnew)<h)
                flag = 1;
            }
            if(flag==0)
            {
              xy_lcoord[i][0] = local[0];
              xy_lcoord[i][1] = local[1];
              xy_flag[i] = 1;
            }
            else
              xy_flag[i] = 0;
          }
          else
            xy_flag[i] = 0;
        }

        for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(xy_Surface)); i++)
        {
          for(j=0; j<3; j++)
            if( (xy_flag[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+0)%3)]>0)
                &&  (xy_flag[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+1)%3)]>0) )
            {
              xy_flag[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+0)%3)] = 2;
              xy_flag[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+1)%3)] = 2;
            }
        }

        for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(xy_Surface)); i++)
          if(xy_flag[i]==2)
          {
            local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(xy_Surface),i,0);
            local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(xy_Surface),i,1);
            Surface_Local2Global(xy_Surface, global, local);

            mi = Project2Surface(theSurface, global, local, xy_n);
            Surface_Local2Global(theSurface, global, local);

            if(LGM_DEBUG)
              printf("%d %d\n", i, nlp);
            xy_new_id[i] = nlp;
            AddLinePoint(nlp, global[0], global[1], global[2]);
            nlp++;
          }
        for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(xy_Surface)); i++)
        {
          dummy = 0;
          for(j=0; j<3; j++)
          {
            a = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+0)%3);
            b = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+1)%3);
            /*						printf("%d %d\n", a, b);*/
            if( (xy_flag[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+0)%3)]==2)
                &&  (xy_flag[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+1)%3)]==2) )
            {
              if(LGM_DEBUG)
                printf("%d %d\n", xy_new_id[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+0)%3)],
                       xy_new_id[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+1)%3)]);
              AddLineSegment( xy_new_id[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+0)%3)],
                              xy_new_id[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(xy_Surface), i, (j+1)%3)]);
              dummy++;
            }
          }
          if(dummy==3)
            xy_element_id[i] = 1;
        }
      }

  UserWriteF("Triangulate surface %4d: ", LGM_SURFACE_ID(theSurface));

  if(GenerateSurfaceGrid (Heap, MarkKey, theSurface, h, 5, 1, D))
  {
    surface_error = 1;
    surf_err = 1;
  }
  else
  {
    j = 0;
    for(i=0; i<nline_points; i++)
    {
      if((INT)lcoord_lpoints[i][0]!=-1)
      {
        LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,0) = lcoord_lpoints[i][0];
        LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,1) = lcoord_lpoints[i][1];
        j++;
      }
    }

    UserWriteF("\n");
    UserWriteF(" Surface %4d: %4d Triangles\n",
               LGM_SURFACE_ID(theSurface),
               LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)));
  }
  return (0);
}

#else

static INT DiscretizeSurface (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h, LGM_POINT *pointlist, INT norp, INT D, INT MarkKey)
{
  return (1);
}

#endif

INT TEST(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  INT i, ilocal,id;
  DOUBLE sloc;
  LGM_LINE *theLine;
  if(local[0]<0.0)
  {
    /* point is on the boundary of the surface */
    id = - (INT)floor(local[0]) - 2;
    ilocal = (INT)floor(local[1]);
    sloc = local[1]-ilocal;
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      if(LGM_LINE_ID(LGM_SURFACE_LINE(theSurface, i))==id)
        theLine = LGM_SURFACE_LINE(theSurface, i);
    Line_Local2GlobalNew(theLine,global,local[1]);
  }
  return(0);
}

INT NS_DIM_PREFIX Surface_Local2Global (LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  INT i, ilocal,ilocal1,id;
  DOUBLE slocal[2],sloc;
  LGM_LINE *theLine;

  /* OS_CHANGED: Changed by Andreas */
  for (i=0; i<2; i++)
    if(floor(local[i] + SMALL_FAK) != floor(local[i]))
    {
      local[i] += SMALL_FAK;
      local[i] = floor(local[i]);
    }

  if(local[0]<0.0)
  {
    /* point is on the boundary of the surface */
    id = - (INT)floor(local[0]) - 2;
    ilocal = (int)floor(local[1]);
    sloc = local[1]-ilocal;
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      if(LGM_LINE_ID(LGM_SURFACE_LINE(theSurface, i))==id)
        theLine = LGM_SURFACE_LINE(theSurface, i);
    Line_Local2GlobalNew(theLine,global,local[1]);
  }
  else
  {
    ilocal = (int)floor(local[0]);
    ilocal1 = (int)floor(local[1]);
    if(ilocal>ilocal1)
      ilocal = ilocal1;
    slocal[0] = local[0]-ilocal;
    slocal[1] = local[1]-ilocal;

    assert(slocal[0]>=0.0);
    assert(slocal[1]>=0.0);
    assert(ilocal<LGM_SURFACE_NTRIANGLE(theSurface) && ilocal>=0);
    if (ilocal<LGM_SURFACE_NTRIANGLE(theSurface))
    {
      global[0] = slocal[0]                           *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[0]
                  + slocal[1]                           *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[0]
                  + (1-slocal[0]-slocal[1])     *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[0];
      global[1] = slocal[0]                           *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[1]
                  + slocal[1]                           *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[1]
                  + (1-slocal[0]-slocal[1])     *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[1];
      global[2] = slocal[0]                           *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[2]
                  + slocal[1]                           *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[2]
                  + (1-slocal[0]-slocal[1])     *LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[2];
    }
  }

  return (0);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDP_Global (BNDP *aBndP, DOUBLE *global)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE *local;
#ifdef NO_PROJECT
  DOUBLE *g;
  INT k;
#endif

  theBndP = BNDP2LGM(aBndP);
  theSurface = LGM_BNDP_SURFACE(theBndP,0);

  assert(LGM_BNDP_N(theBndP));
  if(theSurface==NULL)
    assert(0);
        #ifdef NO_PROJECT
  g = LGM_BNDP_GLOBAL(theBndP,0);
  for(k=0; k<3; k++)
    global[k] = g[k];
        #else
  local = LGM_BNDP_LOCAL(theBndP,0);
  Surface_Local2Global (theSurface,global,local);
        #endif
  return (0);
}

/* domain interface function: for description see domain.h */
static INT BNDP_Globali (BNDP *aBndP, DOUBLE *global, INT i)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE *local;
#ifdef NO_PROJECT
  DOUBLE *g;
  INT k;
#endif

  theBndP = BNDP2LGM(aBndP);
  theSurface = LGM_BNDP_SURFACE(theBndP,i);

        #ifdef NO_PROJECT
  g = LGM_BNDP_GLOBAL(theBndP,i);
  for(k=0; k<3; k++)
    global[k] = g[k];
        #else
  local = LGM_BNDP_LOCAL(theBndP,i);
  Surface_Local2Global (theSurface,global,local);
        #endif

  return (0);
}

/* domain interface function: for description see domain.h */
INT BNDP_BndCond (BNDP *aBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE global[DOM_PARAM_OFFSET],*local;
  INT ilocal=0;
#ifdef NO_PROJECT
  DOUBLE *g;
  INT k;
#endif

  /* general */
  theBndP = BNDP2LGM(aBndP);
  *n = LGM_BNDP_N(theBndP);
  if (i<0 || i>=LGM_BNDP_N(theBndP))
    ilocal = ilocal;
  assert(i>=0 && i<LGM_BNDP_N(theBndP));
  theSurface = LGM_BNDP_SURFACE(theBndP,i);
  if (LGM_SURFACE_BNDCOND(theSurface)==NULL)
    return (2);

  /* global coordinates */
        #ifdef NO_PROJECT
  g = LGM_BNDP_GLOBAL(theBndP,i);
  for(k=0; k<3; k++)
    global[k] = g[k];
        #else
  local = LGM_BNDP_LOCAL(theBndP,i);
  Surface_Local2Global(theSurface,global,local);
        #endif

  /* get values */
  if (in!=NULL)
  {
    in[0] = global[0];
    in[1] = global[1];
    in[2] = global[2];
    in[DIM] = LGM_SURFACE_ID(theSurface);
    if ((*LGM_SURFACE_BNDCOND (theSurface))(in,value,type))
      return (1);
  }
  else
  {
    global[DIM] = LGM_SURFACE_ID(theSurface);
    if ((*LGM_SURFACE_BNDCOND (theSurface))(global,value,type))
      return (1);
  }

  return (0);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDP_SurfaceId (BNDP *aBndP, INT *n, INT i)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  INT ilocal=0;

  /* general */
  theBndP = BNDP2LGM(aBndP);
  *n = LGM_BNDP_N(theBndP);
  if (i<0 || i>=LGM_BNDP_N(theBndP))
    ilocal = ilocal;
  assert(i>=0 && i<LGM_BNDP_N(theBndP));
  theSurface = LGM_BNDP_SURFACE(theBndP,i);

  return LGM_SURFACE_ID(theSurface);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDP_BndPDesc (BNDP *aBndP, INT *move, INT *part)
{
  LGM_BNDP *theBndP;

  part[0] = 0;
  theBndP = BNDP2LGM(aBndP);
  if(LGM_BNDP_N(theBndP)==1)
    *move=1;                                            /* Point in the Interior of the surface */
  else
    *move=0;                                            /* Point on the Boundary of the surface */

  *move=0;

  /* HRR_TODO: assign part */
  *part = 0;

  return(0);
}

static INT Check_Local_Coord(LGM_SURFACE *theSurface, DOUBLE *local)
{
  INT ilocal, ilocal1;
  DOUBLE slocal[2];

  ilocal = (INT)floor(local[0]);
  ilocal1 = (INT)floor(local[1]);
  if(ilocal>ilocal1)
    ilocal = ilocal1;
  slocal[0] = local[0]-ilocal;
  slocal[1] = local[1]-ilocal;
  /* Punkt liegt im Innern des Dreiecks */
  if( (slocal[0]>-SMALL) && (slocal[0]<1.0+SMALL)
      && (slocal[1]>-SMALL) && (slocal[1]<1.0+SMALL)
      && (1-slocal[0]-slocal[1]>-SMALL) && (1-slocal[0]-slocal[1]<1.0+SMALL)
      && (ilocal<LGM_SURFACE_NTRIANGLE(theSurface)) )
    return(1);
  else
    assert(0);

  return(0);
}

BNDS* NS_DIM_PREFIX BNDP_CreateBndS (HEAP *Heap, BNDP **aBndP, INT n)
{
  INT i,j,k,q,i0,j0,k0,q0,count, ilocal,ilocal1;
  LGM_BNDP *theBndP1, *theBndP2, *theBndP3, *theBndP4;
  LGM_SURFACE *theSurface;
  LGM_BNDS *theBndS;
  DOUBLE local[2], nv[3];
  DOUBLE globalp0[3],globalp1[3],globalp2[3], globalp3[3], global[3], globalnew[3];
  double sp, d, min_d;
  DOUBLE A[3], B[3], BNDP_NV[3], Surface_NV[3];
  INT mi;

  if(n!=3 && n!=4)
  {
    return(NULL);
  }

  theBndP1 = BNDP2LGM(aBndP[0]);
  theBndP2 = BNDP2LGM(aBndP[1]);
  theBndP3 = BNDP2LGM(aBndP[2]);
  if (n==4) theBndP4 = BNDP2LGM(aBndP[3]);

  BNDP_Global(aBndP[0],globalp0);
  BNDP_Global(aBndP[1],globalp1);
  BNDP_Global(aBndP[2],globalp2);
  if (n==4) BNDP_Global(aBndP[3],globalp3);


  A[0] = globalp2[0] - globalp0[0];
  A[1] = globalp2[1] - globalp0[1];
  A[2] = globalp2[2] - globalp0[2];
  B[0] = globalp2[0] - globalp1[0];
  B[1] = globalp2[1] - globalp1[1];
  B[2] = globalp2[2] - globalp1[2];
  V3_VECTOR_PRODUCT(A, B, BNDP_NV);
  // This is wrong. OS_CHANGED
  //	V_DIM_SCALE(1.0/sqrt(V_DIM_SCAL_PROD(BNDP_NV, BNDP_NV)),BNDP_NV);

  // We don't need this. OS_CHANGED
  //	  /* check */
  //	  for(i=0;i<3;i++)
  //		  m1[i] = globalp2[i]-globalp0[i];
  //	  for(i=0;i<3;i++)
  //		  m2[i] = globalp1[i]-globalp0[i];
  //	  m3[0] = m1[1]*m2[2] - m1[2]*m2[1];
  //	  m3[1] = m1[2]*m2[0] - m1[0]*m2[2];
  //	  m3[2] = m1[0]*m2[1] - m1[1]*m2[0];
  //
  //	  area = 0.0;
  //	  for(i=0;i<3;i++)
  //		  area = area + m3[i]*m3[i];
  //	  area = sqrt(area);
  //	  if(area<0.0)
  //		  area = - area;

  if(E_Distance(globalp2, globalp0)<SMALL)
    assert(E_Distance(globalp2, globalp0)>SMALL);
  if(E_Distance(globalp2, globalp1)<SMALL)
    assert(E_Distance(globalp2, globalp1)>SMALL);
  if(E_Distance(globalp1, globalp0)<SMALL)
    assert(E_Distance(globalp1, globalp0)>SMALL);
  if (n==4)
  {
    if(E_Distance(globalp3, globalp0)<SMALL)
      assert(E_Distance(globalp1, globalp0)>SMALL);
    if(E_Distance(globalp3, globalp1)<SMALL)
      assert(E_Distance(globalp1, globalp0)>SMALL);
    if(E_Distance(globalp3, globalp2)<SMALL)
      assert(E_Distance(globalp1, globalp0)>SMALL);
  }

  if (n==3)
  {
    global[0] = ( globalp0[0] + globalp1[0] +  globalp2[0] ) / 3.0;
    global[1] = ( globalp0[1] + globalp1[1] +  globalp2[1] ) / 3.0;
    global[2] = ( globalp0[2] + globalp1[2] +  globalp2[2] ) / 3.0;
  }
  if (n==4)
  {
    global[0] = ( globalp0[0] + globalp1[0] +  globalp2[0]  +  globalp3[0]) / 4.0;
    global[1] = ( globalp0[1] + globalp1[1] +  globalp2[1]  +  globalp3[1]) / 4.0;
    global[2] = ( globalp0[2] + globalp1[2] +  globalp2[2]  +  globalp3[2]) / 4.0;
  }
  // We don't need this. OS_CHANGED
  //	  count = 0;
  //	  min_d = 1000000.0;
  //	  d =  1000000.0;

  count = 0;
  min_d = DBL_MAX;
  if (n==3)
    for (i=0; i<LGM_BNDP_N(theBndP1); i++)
      for (j=0; j<LGM_BNDP_N(theBndP2); j++)
        for (k=0; k<LGM_BNDP_N(theBndP3); k++)
        {
          if((LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j))
             && (LGM_BNDP_SURFACE(theBndP2,j)==LGM_BNDP_SURFACE(theBndP3,k)) )
          {
            count++;
            theSurface = LGM_BNDP_SURFACE(theBndP1,i);
            i0=i;
            j0=j;
            k0=k;
          }
        }

  if (n==4)
    for (i=0; i<LGM_BNDP_N(theBndP1); i++)
      for (j=0; j<LGM_BNDP_N(theBndP2); j++)
        for (k=0; k<LGM_BNDP_N(theBndP3); k++)
          for (q=0; q<LGM_BNDP_N(theBndP4); q++)
          {
            if ((LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j)) &&
                (LGM_BNDP_SURFACE(theBndP2,j)==LGM_BNDP_SURFACE(theBndP3,k)) &&
                (LGM_BNDP_SURFACE(theBndP3,k)==LGM_BNDP_SURFACE(theBndP4,q)) )
            {
              count++;
              theSurface = LGM_BNDP_SURFACE(theBndP1,i);
              i0=i;
              j0=j;
              k0=k;
              q0=q;
            }
          }


  if (count==0)
  {
    UserWrite("BNDP_CreateBndS(): No common surface found.\n");
    return(NULL);
  }
  if(count>1 && n==3)
  {
    count = 0;
    min_d = DBL_MAX;
    for (i=0; i<LGM_BNDP_N(theBndP1); i++)
      for (j=0; j<LGM_BNDP_N(theBndP2); j++)
        for (k=0; k<LGM_BNDP_N(theBndP3); k++)
        {
          if((LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j))
             && (LGM_BNDP_SURFACE(theBndP2,j)==LGM_BNDP_SURFACE(theBndP3,k)) )
          {
            /* Zusatzabfrage fuer den Fall, dass 3 Punkte auf zwei Surfaces liegen */
            theSurface = LGM_BNDP_SURFACE(theBndP1,i);
            V_DIM_CLEAR(nv);
            mi = GetLocalKoord(theSurface, global, local, nv);
            if(mi!=-1)
            {
              Surface_Local2Global(theSurface, globalnew, local);
              d = E_Distance(global, globalnew);
            }
            else
              d = DBL_MAX;
            if(d<min_d)
            {
              if( Check_Local_Coord(theSurface, local) )
              {
                min_d = d;
                theSurface = LGM_BNDP_SURFACE(theBndP1,i);
                count++;
                i0=i;
                j0=j;
                k0=k;
              }
            }
          }
        }
  }
  else if (count>1 && n==4)
  {
    count = 0;
    min_d = DBL_MAX;
    for (i=0; i<LGM_BNDP_N(theBndP1); i++)
      for (j=0; j<LGM_BNDP_N(theBndP2); j++)
        for (k=0; k<LGM_BNDP_N(theBndP3); k++)
          for (q=0; q<LGM_BNDP_N(theBndP4); q++)
          {
            if((LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j))
               && (LGM_BNDP_SURFACE(theBndP2,j)==LGM_BNDP_SURFACE(theBndP3,k))
               && (LGM_BNDP_SURFACE(theBndP3,k)==LGM_BNDP_SURFACE(theBndP4,q)) )
            {
              /* Zusatzabfrage fuer den Fall, dass 3 Punkte auf zwei Surfaces liegen */
              theSurface = LGM_BNDP_SURFACE(theBndP1,i);
              V_DIM_CLEAR(nv);
              mi = GetLocalKoord(theSurface, global, local, nv);
              if(mi!=-1)
              {
                Surface_Local2Global(theSurface, globalnew, local);
                d = E_Distance(global, globalnew);
              }
              else
                d = DBL_MAX;
              if(d<min_d)
              {
                if( Check_Local_Coord(theSurface, local) )
                {
                  min_d = d;
                  theSurface = LGM_BNDP_SURFACE(theBndP1,i);
                  count++;
                  i0=i;
                  j0=j;
                  k0=k;
                }
              }
            }
          }
  }
  if (count==0) return (NULL);

  theSurface = LGM_BNDP_SURFACE(theBndP1,i0);

  theBndS = (LGM_BNDS *)GetFreelistMemory(Heap,sizeof(LGM_BNDS));
  assert(theBndS!=NULL);
  LGM_BNDS_SURFACE(theBndS) = theSurface;

        #ifdef NO_PROJECT
  LGM_BNDS_GLOBAL(theBndS,0,0) = LGM_BNDP_GLOBAL(theBndP1,i0)[0];
  LGM_BNDS_GLOBAL(theBndS,0,1) = LGM_BNDP_GLOBAL(theBndP1,i0)[1];
  LGM_BNDS_GLOBAL(theBndS,0,2) = LGM_BNDP_GLOBAL(theBndP1,i0)[2];
  LGM_BNDS_GLOBAL(theBndS,1,0) = LGM_BNDP_GLOBAL(theBndP2,j0)[0];
  LGM_BNDS_GLOBAL(theBndS,1,1) = LGM_BNDP_GLOBAL(theBndP2,j0)[1];
  LGM_BNDS_GLOBAL(theBndS,1,2) = LGM_BNDP_GLOBAL(theBndP2,j0)[2];
  LGM_BNDS_GLOBAL(theBndS,2,0) = LGM_BNDP_GLOBAL(theBndP3,k0)[0];
  LGM_BNDS_GLOBAL(theBndS,2,1) = LGM_BNDP_GLOBAL(theBndP3,k0)[1];
  LGM_BNDS_GLOBAL(theBndS,2,2) = LGM_BNDP_GLOBAL(theBndP3,k0)[2];
  if (n==4)
  {
    LGM_BNDS_GLOBAL(theBndS,3,0) = LGM_BNDP_GLOBAL(theBndP4,q0)[0];
    LGM_BNDS_GLOBAL(theBndS,3,1) = LGM_BNDP_GLOBAL(theBndP4,q0)[1];
    LGM_BNDS_GLOBAL(theBndS,3,2) = LGM_BNDP_GLOBAL(theBndP4,q0)[2];
  }
        #else
  LGM_BNDS_LOCAL(theBndS,0,0) = LGM_BNDP_LOCAL(theBndP1,i0)[0];
  LGM_BNDS_LOCAL(theBndS,0,1) = LGM_BNDP_LOCAL(theBndP1,i0)[1];
  LGM_BNDS_LOCAL(theBndS,1,0) = LGM_BNDP_LOCAL(theBndP2,j0)[0];
  LGM_BNDS_LOCAL(theBndS,1,1) = LGM_BNDP_LOCAL(theBndP2,j0)[1];
  LGM_BNDS_LOCAL(theBndS,2,0) = LGM_BNDP_LOCAL(theBndP3,k0)[0];
  LGM_BNDS_LOCAL(theBndS,2,1) = LGM_BNDP_LOCAL(theBndP3,k0)[1];
  if (n==4) {
    LGM_BNDS_LOCAL(theBndS,3,0) = LGM_BNDP_LOCAL(theBndP4,q0)[0];
    LGM_BNDS_LOCAL(theBndS,3,1) = LGM_BNDP_LOCAL(theBndP4,q0)[1];
  }
        #endif
  /* lege Richtung fuer die BNDS fest */

  V_DIM_CLEAR(nv);
  GetLocalKoord(theSurface, global, local, nv);
  ilocal = (INT)floor(local[0]);
  ilocal1 = (INT)floor(local[1]);
  if(ilocal>ilocal1)
    ilocal = ilocal1;

  /* bestimme den Normalenvektor der des Surfacedreieck */
  A[0] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[0]
         - LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[0];
  A[1] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[1]
         - LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[1];
  A[2] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[2]
         - LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[2];
  B[0] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[0]
         - LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[0];
  B[1] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[1]
         - LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[1];
  B[2] = LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[2]
         - LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[2];

  V3_VECTOR_PRODUCT(A, B, Surface_NV);

  V3_SCALAR_PRODUCT(BNDP_NV, Surface_NV, sp);


  if(sp>0.0)
    LGM_BNDS_N(theBndS) = n;
  else
    LGM_BNDS_N(theBndS) = -n;

  return((BNDS *)theBndS);
}

static DOUBLE LinePointDistance_local(LGM_LINE *theLine, DOUBLE *globalp1, DOUBLE *globalp2, DOUBLE *start_local, DOUBLE *end_local)
{
  int i, start_i, end_i;
  DOUBLE length_of_segment, start_s, end_s;

  start_i = (int)floor(*start_local);
  start_s = *start_local - start_i;
  end_i = (int)floor(*end_local);
  end_s = *end_local - end_i;

  length_of_segment = 0;

  for(i=start_i+1; i<end_i; i++)
    length_of_segment = length_of_segment + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));

  if(end_i<LGM_LINE_NPOINT(theLine)-1)
    length_of_segment = length_of_segment +
                        end_s * LGM_POINT_DIST(LGM_LINE_POINT(theLine,end_i),LGM_LINE_POINT(theLine,end_i+1));
  if(end_i>start_i)
    length_of_segment = length_of_segment +
                        (1.0-start_s) * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));
  else
    length_of_segment = length_of_segment -
                        start_s * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));

  return(length_of_segment);
}


static int Find_Midpoint(LGM_LINE *theLine, DOUBLE *globalp1, DOUBLE *globalp2, DOUBLE *global, DOUBLE *start_local, DOUBLE *end_local)
{
  int i, start_i, end_i, help_i;
  DOUBLE length_of_segment, start_s, end_s, help_s;
  DOUBLE mid_p[3], len, mid_i, mid_local, local, dist, *help, *help_local;

  for(i=0; i<3; i++)
    mid_p[i] = 0.5 * (globalp1[i] + globalp2[i]);

  start_i = (int)floor(*start_local);
  start_s = *start_local - start_i;
  end_i = (int)floor(*end_local);
  end_s = *end_local - end_i;
  if(start_i>end_i)
  {
    help_i = start_i;
    help_s = start_s;
    help_local = start_local;
    help = globalp1;

    start_i = end_i;
    start_s = end_s;
    start_local = end_local;
    globalp1 = globalp2;

    end_i = help_i;
    end_s = help_s;
    end_local = help_local;
    globalp2 = help;
  }
  else
  if((start_s>end_s)&&(start_i==end_i))
  {
    help_s = start_s;
    help_local = start_local;
    help = globalp1;

    start_s = end_s;
    start_local = end_local;
    globalp1 = globalp2;

    end_s = help_s;
    end_local = help_local;
    globalp2 = help;
  }

  length_of_segment = LinePointDistance_local(theLine, globalp1, globalp2, start_local, end_local);
  dist = length_of_segment / 2;

  len = 0;
  mid_i = -1;
  if(start_i==end_i)
  {
    len = (end_s - start_s) * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));
    mid_i = start_i;
    mid_local = (end_s + start_s) / 2;
  }
  else
  {
    len = (1.0 - start_s) * LGM_POINT_DIST(LGM_LINE_POINT(theLine,start_i),LGM_LINE_POINT(theLine,start_i+1));
    if(len>dist)
    {
      mid_i = start_i;
      mid_local = dist/len * (1-start_s) + start_s;
    }
    else
    {
      if(start_i+1==end_i)
      {
        mid_i = end_i;
        mid_local = (dist - len) / LGM_POINT_DIST(LGM_LINE_POINT(theLine,end_i),LGM_LINE_POINT(theLine,end_i+1));
      }
      else
        for(i=start_i+1; i<end_i; i++)
          if(len + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1)) < dist)
            len = len + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));
          else
          {
            mid_i = i;
            mid_local = (dist - len) / LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));
            break;
          }
      if(mid_i==-1)
      {
        mid_i = end_i;
        mid_local = (dist - len) / LGM_POINT_DIST(LGM_LINE_POINT(theLine,end_i),LGM_LINE_POINT(theLine,end_i+1));
      }
    }
  }

  local = (DOUBLE)mid_i + mid_local;
  Line_Local2GlobalNew(theLine,global,local);
  return(0);
}

static INT Points_in_Line(LGM_LINE *theLine, DOUBLE *globalp1, DOUBLE *globalp2)
{
  DOUBLE local1,  local2;

  Line_Global2Local(theLine, globalp1, &local1);
  Line_Global2Local(theLine, globalp2, &local2);
  if( (local1>0.0-SMALL) && (local2>0.0-SMALL))
    return(1);
  else
    return(0);
}

#define MAX_SURFACES 20
#define MAX_LINES 20

static DOUBLE Nearest_Surface(LGM_SURFACE **Surfaces, LGM_SURFACE **theNewSurface, INT count, DOUBLE *global)
{
  DOUBLE min_d, d, nv[3], local[2], globalnew[3];
  INT i, mi;
  LGM_SURFACE *theSurface;

  min_d = DBL_MAX;
  for(i=0; i<count; i++)
  {
    theSurface = Surfaces[i];
    V_DIM_CLEAR(nv);
    mi = GetLocalKoord(theSurface,global,local, nv);
    if(mi!=-1)
    {
      Surface_Local2Global(theSurface, globalnew, local);
      d = sqrt((global[0]-globalnew[0])*(global[0]-globalnew[0])
               +    (global[1]-globalnew[1])*(global[1]-globalnew[1])
               +    (global[2]-globalnew[2])*(global[2]-globalnew[2]));
      if(min_d>d)
      {
        *theNewSurface = theSurface;
        min_d = d;
      }
    }
  }
  return(min_d);
}

static INT Count_Common_Surfaces(LGM_BNDP *theBndP1, LGM_BNDP *theBndP2, LGM_SURFACE **Surfaces)
{
  INT i, j, count;
  /* finde die Surfaces, die theBndP1 und theBndP2 enthalten */
  count = 0;
  for (i=0; i<LGM_BNDP_N(theBndP1); i++)
    for (j=0; j<LGM_BNDP_N(theBndP2); j++)
      if (LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j))
      {
        Surfaces[count] = LGM_BNDP_SURFACE(theBndP1,i);
        count++;
      }
  return(count);
}

static INT Count_Common_Lines(LGM_BNDP *theBndP1, LGM_BNDP *theBndP2, LGM_SURFACE **Surfaces, LGM_LINE **Lines, INT count)
{
  INT i, j, k, nlines, flag;
  LGM_LINE *theLine;
  DOUBLE globalp1[3],globalp2[3];

  BNDP_Global((BNDP *)theBndP1,globalp1);
  BNDP_Global((BNDP *)theBndP2,globalp2);

  nlines = 0;
  for(i=0; i<count; i++)
    for(j=0; j<LGM_SURFACE_NLINE(Surfaces[i]); j++)
    {
      theLine = LGM_SURFACE_LINE(Surfaces[i], j);
      if(Points_in_Line(theLine, globalp1, globalp2))
      {
        flag = 1;
        for(k=0; k<nlines; k++)
          if(theLine==Lines[k])
            flag = 0;
        if(flag)
        {
          Lines[nlines] = theLine;
          nlines++;
        }
      }
    }
  return(nlines);
}

/* domain interface function: for description see domain.h */
BNDP* NS_DIM_PREFIX BNDP_InsertBndP (HEAP *Heap, BVP *aBVP, double *global)
{
  LGM_DOMAIN *theDomain;
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  int nline;
  double local[2];

  theDomain = BVP2LGM(aBVP);

  UserWriteF("B %lf %lf %lf\n", global[0],global[1],global[2]);

  nline = 0;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    Line_Global2Local(theLine,global,local);

    if (*local >= 0.0)
    {
      UserWriteF("	L %d %f\n", LGM_LINE_ID(theLine),(FLOAT)*local);
      nline++;
    }
  }

  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    INT mi;
    DOUBLE d,nv[3],globalnew[3];

    V_DIM_CLEAR(nv);
    mi = GetLocalKoord(theSurface,global,local, nv);
    if(mi!=-1)
    {
      Surface_Local2Global(theSurface, globalnew, local);
      d = sqrt((global[0]-globalnew[0])*(global[0]-globalnew[0])
               +    (global[1]-globalnew[1])*(global[1]-globalnew[1])
               +    (global[2]-globalnew[2])*(global[2]-globalnew[2]));
      if(d<SMALL)
      {
        UserWriteF("	S %d %d %f %f\n",
                   LGM_SURFACE_ID(theSurface),
                   mi,
                   (FLOAT)(local[1]-(DOUBLE)mi),
                   (FLOAT)(1.0-(local[0]-(DOUBLE)mi)-(local[1]-(DOUBLE)mi)));
      }
    }
  }
  UserWriteF(";\n");

  return (NULL);
}

/* domain interface function: for description see domain.h */
BNDP* NS_DIM_PREFIX BNDP_CreateBndP (HEAP *Heap, BNDP *aBndP0, BNDP *aBndP1, DOUBLE lcoord)
{
  LGM_BNDP *theBndP1, *theBndP2, *theBndP;
  LGM_SURFACE *Surfaces[MAX_SURFACES], *theNewSurface;
  LGM_LINE *theLine;
  INT i,j, count, max, found, iold, jold, mi;
  DOUBLE globalp1[3],globalp2[3],global[3], local[2], nv[3], local1, local2, globalnew[3], midp[3], newlocal;
  DOUBLE localp1, localp2;
#ifdef NO_PROJECT
  DOUBLE len, value, n[3];
#endif

  theBndP1 = BNDP2LGM(aBndP0);
  theBndP2 = BNDP2LGM(aBndP1);
  max = LGM_BNDP_N(theBndP1);
  if(max<LGM_BNDP_N(theBndP2))
    max = LGM_BNDP_N(theBndP2);

  BNDP_Global(aBndP0,globalp1);
  BNDP_Global(aBndP1,globalp2);

  /* check */
  assert(E_Distance(globalp1, globalp2)>SMALL);

  midp[0] = global[0] = ( globalp1[0] + globalp2[0] ) / 2;
  midp[1] = global[1] = ( globalp1[1] + globalp2[1] ) / 2;
  midp[2] = global[2] = ( globalp1[2] + globalp2[2] ) / 2;

  assert(lcoord>0.0 && lcoord<1.0);
  theBndP1 = BNDP2LGM(aBndP0);
  theBndP2 = BNDP2LGM(aBndP1);

  count = Count_Common_Surfaces(theBndP1, theBndP2, Surfaces);
  if(count==0)
    return(NULL);

  found = 0;
  for (i=0; i<LGM_BNDP_NLINE(theBndP1); i++)
    for (j=0; j<LGM_BNDP_NLINE(theBndP2); j++)
      if (LGM_BNDP_LINE(theBndP1,i)==LGM_BNDP_LINE(theBndP2,j))
      {
                                #ifdef NO_PROJECT
        /* nothing to do */
                                #else
        Line_Global2Local(LGM_BNDP_LINE(theBndP1, i), globalp1, &localp1);
        Line_Global2Local(LGM_BNDP_LINE(theBndP2, j), globalp2, &localp2);
                                #endif
        /* teste, ob BndP's benachbart sind */
        if(LGM_LINE_BEGIN(LGM_BNDP_LINE(theBndP1, i))==LGM_LINE_END(LGM_BNDP_LINE(theBndP1, i)))
        {
          /* line cyclic */

        }
        else
        {
          /* line not cyclic */
                                        #ifdef NO_PROJECT

          if( (E_Distance(LGM_BNDP_LINE_GLOBALRIGHT(theBndP1, i), globalp2)<SMALL)
              ||  (E_Distance(LGM_BNDP_LINE_GLOBALLEFT (theBndP2, j), globalp1)<SMALL)
              ||  (E_Distance(LGM_BNDP_LINE_GLOBALLEFT (theBndP1, i), globalp2)<SMALL)
              ||  (E_Distance(LGM_BNDP_LINE_GLOBALRIGHT(theBndP2, j), globalp1)<SMALL) )
          {
            theLine = LGM_BNDP_LINE(theBndP1, i);
            found++;
          }
                                        #else
          if((    (ABS(LGM_BNDP_LINE_RIGHT(theBndP1, i) - localp2)<SMALL)
                  ||  (ABS(LGM_BNDP_LINE_LEFT (theBndP2, j) - localp1)<SMALL) )
             || ((ABS(LGM_BNDP_LINE_LEFT (theBndP1, i) - localp2)<SMALL)
                 ||  (ABS(LGM_BNDP_LINE_RIGHT(theBndP2, j) - localp1)<SMALL) ))
          {
            theLine = LGM_BNDP_LINE(theBndP1, i);
            found++;
          }
                                        #endif
        }
      }
  if(found>1)
    assert(found<2);

  if(found==1)                                  /* projeziere auf die Linie */
  {
                #ifdef NO_PROJECT
    global[0] = (globalp1[0] + globalp2[0]) / 2;
    global[1] = (globalp1[1] + globalp2[1]) / 2;
    global[2] = (globalp1[2] + globalp2[2]) / 2;
#ifdef OCC_GEOMETRY
    int i = LGM_LINE_ID(theLine);
    double ppt[3],dist,err;
    //cout << "P  " << global[0] << " " << global[1] << " " << global[2] << endl;
    err = occ_geom.ProjectPointOnEdge(LineMap[i],global,global,&dist);
    //cout << "PP " << global[0] << " " << global[1] << " " << global[2] << " " << dist << endl;
    if (err != 0)
    {
      assert(0);
    }
    midp[0] = global[0];
    midp[1] = global[1];
    midp[2] = global[2];
#endif
                #else
    Line_Global2Local(theLine, globalp1, &local1);
    Line_Global2Local(theLine, globalp2, &local2);
    global[0] = global[1] = global[2] = 0.0;
    Find_Midpoint(theLine, globalp1, globalp2, global, &local1, &local2);
    Line_Global2Local(theLine, global, &newlocal);
                #endif
    theBndP = (LGM_BNDP*)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
    if(theBndP==NULL)
    {
      assert(0);
      return(NULL);
    }

    theBndP->Line = (LGM_BNDP_PLINE*)GetFreelistMemory(Heap,sizeof(LGM_BNDP_PLINE));
    if(theBndP->Line==NULL)
    {
      assert(0);
      return(NULL);
    }
    LGM_BNDP_NLINE(theBndP) = 0;
    /* set local coord of left and right point */
    LGM_BNDP_LINE(theBndP, LGM_BNDP_NLINE(theBndP)) = theLine;
    LGM_BNDP_NLINE(theBndP)++;
    /* change left and right in bndp1,bndp2 */
    for(i=0; i<LGM_BNDP_NLINE(theBndP1); i++)
      if(LGM_BNDP_LINE(theBndP1, i)==theLine)
        iold = i;
    for(i=0; i<LGM_BNDP_NLINE(theBndP2); i++)
      if(LGM_BNDP_LINE(theBndP2, i)==theLine)
        jold = i;
                #ifdef NO_PROJECT
    V_DIM_SUBTRACT(globalp2, globalp1, n);
    V_DIM_EUKLIDNORM(n,len);
    V_DIM_SCALE(1/len, n);
    value = n[0] * (global[0] - globalp1[0]) +
            n[1] * (global[1] - globalp1[1]) +
            n[2] * (global[2] - globalp1[2]);
    if(value>0.0)
    {
      LGM_BNDP_LINE_GLOBALLEFT (theBndP, 0)[0] = globalp1[0];
      LGM_BNDP_LINE_GLOBALLEFT (theBndP, 0)[1] = globalp1[1];
      LGM_BNDP_LINE_GLOBALLEFT (theBndP, 0)[2] = globalp1[2];
      LGM_BNDP_LINE_GLOBALRIGHT(theBndP, 0)[0] = globalp2[0];
      LGM_BNDP_LINE_GLOBALRIGHT(theBndP, 0)[1] = globalp2[1];
      LGM_BNDP_LINE_GLOBALRIGHT(theBndP, 0)[2] = globalp2[2];
    }
    else
    {
      LGM_BNDP_LINE_GLOBALLEFT (theBndP, 0)[0] = globalp2[0];
      LGM_BNDP_LINE_GLOBALLEFT (theBndP, 0)[1] = globalp2[1];
      LGM_BNDP_LINE_GLOBALLEFT (theBndP, 0)[2] = globalp2[2];
      LGM_BNDP_LINE_GLOBALRIGHT(theBndP, 0)[0] = globalp1[0];
      LGM_BNDP_LINE_GLOBALRIGHT(theBndP, 0)[1] = globalp1[1];
      LGM_BNDP_LINE_GLOBALRIGHT(theBndP, 0)[2] = globalp1[2];
    }
                #else
    if( (local1<newlocal) && (newlocal<local2) )
    {
      LGM_BNDP_LINE_LEFT(theBndP, 0) = local1;
      LGM_BNDP_LINE_RIGHT(theBndP, 0) = local2;
    }
    if( (local1>newlocal) && (newlocal>local2) )
    {
      LGM_BNDP_LINE_LEFT(theBndP, 0) = local2;
      LGM_BNDP_LINE_RIGHT(theBndP, 0) = local1;
    }
                #endif

    max = 0;
    for(i=0; i<count; i++)
      for(j=0; j<LGM_SURFACE_NLINE(Surfaces[i]); j++)
        if(LGM_SURFACE_LINE(Surfaces[i], j)==theLine)
          max++;
    theBndP->Surf = (LGM_BNDP_PSURFACE*)GetFreelistMemory(Heap,max*sizeof(LGM_BNDP_PSURFACE));
    if(theBndP->Surf==NULL)
    {
      assert(0);
      return(NULL);
    }
    LGM_BNDP_N(theBndP) = 0;

    for(i=0; i<count; i++)
      for(j=0; j<LGM_SURFACE_NLINE(Surfaces[i]); j++)
        if(LGM_SURFACE_LINE(Surfaces[i], j)==theLine)
        {
          LGM_BNDP_SURFACE(theBndP,LGM_BNDP_N(theBndP)) = Surfaces[i];
                                        #ifdef NO_PROJECT
          LGM_BNDP_GLOBAL(theBndP,LGM_BNDP_N(theBndP))[0] = midp[0];
          LGM_BNDP_GLOBAL(theBndP,LGM_BNDP_N(theBndP))[1] = midp[1];
          LGM_BNDP_GLOBAL(theBndP,LGM_BNDP_N(theBndP))[2] = midp[2];
                                        #else
          Line_Global2Local(theLine, global, local);
          LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[0] = - (DOUBLE) LGM_LINE_ID(theLine) - 2.0;
          LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[1] = local[0];
                                        #endif
          LGM_BNDP_N(theBndP)++;
        }
    ASSERT(max==LGM_BNDP_N(theBndP));
  }
  else                                                  /* projeziere auf die Surface */
  {
    /* die BndP's liegen auf keiner Line, koennen aber noch auf mehreren Surfaces liegen */
    /* finde die richtige Surface */
    theBndP = (LGM_BNDP*)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
    if(theBndP==NULL)
    {
      assert(0);
      return(NULL);
    }
    theBndP->Line = NULL;
    LGM_BNDP_NLINE(theBndP) = 0;
    theBndP->Surf = (LGM_BNDP_PSURFACE*)GetFreelistMemory(Heap,sizeof(LGM_BNDP_PSURFACE));
    if(theBndP->Surf==NULL)
    {
      assert(0);
      return(NULL);
    }
    LGM_BNDP_N(theBndP) = 0;

    Nearest_Surface(Surfaces, &theNewSurface, count, global);
    V_DIM_CLEAR(nv);
    mi = GetLocalKoord(theNewSurface,global,local, nv);
    if(mi!=-1)
    {
      Surface_Local2Global(theNewSurface, globalnew, local);
      if( Check_Local_Coord(theNewSurface, local) )
      {
        LGM_BNDP_SURFACE(theBndP,LGM_BNDP_N(theBndP)) = theNewSurface;
                                #ifdef NO_PROJECT
#ifdef OCC_GEOMETRY
        int i = LGM_SURFACE_ID(theNewSurface);
        double ppt[3],dist,err;
        //cout << "P  " << global[0] << " " << global[1] << " " << global[2] << endl;
        err = occ_geom.ProjectPointOnFace(SurfMap[i],global,global,&dist);
        //cout << "PP " << global[0] << " " << global[1] << " " << global[2] << " " << dist << endl;
        if (err != 0)
        {
          assert(0);
        }
        midp[0] = global[0];
        midp[1] = global[1];
        midp[2] = global[2];
#endif
        LGM_BNDP_GLOBAL(theBndP,LGM_BNDP_N(theBndP))[0] = midp[0];
        LGM_BNDP_GLOBAL(theBndP,LGM_BNDP_N(theBndP))[1] = midp[1];
        LGM_BNDP_GLOBAL(theBndP,LGM_BNDP_N(theBndP))[2] = midp[2];
                                #else
        LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[0] = local[0];
        LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[1] = local[1];
                                #endif
        LGM_BNDP_N(theBndP)++;

        if( (local[0]<0.0) || (local[1]<0.0) )
          assert(0);
      }
      else
        assert(0);
    }
    else
      assert(0);
  }

  /* check */
#ifndef OCC_GEOMETRY
  if(E_Distance(globalp1, global)<SMALL)
    assert(E_Distance(globalp1, global)>SMALL);
  if(E_Distance(globalp2, global)<SMALL)
    assert(E_Distance(globalp2, global)>SMALL);
  if(E_Distance(globalp1, globalp2)<SMALL)
    assert(E_Distance(globalp1, globalp2)>SMALL);
#endif

  return((BNDP *)theBndP);

}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDP_Dispose (HEAP *Heap, BNDP *aBndP)
{
  LGM_BNDP *theBndP;

  if (aBndP == NULL) return(0);

  theBndP = BNDP2LGM(aBndP);
  if (LGM_BNDP_NLINE(theBndP)>0)
    PutFreelistMemory(Heap,theBndP->Line,LGM_BNDP_NLINE(theBndP)*sizeof(LGM_BNDP_PLINE));
  ASSERT(LGM_BNDP_N(theBndP)>0);
  PutFreelistMemory(Heap,theBndP->Surf,LGM_BNDP_N(theBndP)*sizeof(LGM_BNDP_PSURFACE));

  return (PutFreelistMemory(Heap,theBndP,sizeof(LGM_BNDP)));
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDP_SaveBndP (BNDP *aBndP)
{
  INT i;
  LGM_BNDP *theBndP;
  int n;
  double d[3], e;
#ifdef NO_PROJECT
  double *global;
#endif

  theBndP = BNDP2LGM(aBndP);
  /* the lines */
  n = LGM_BNDP_NLINE(theBndP);
  if (Bio_Write_mint(1,&n)) return (1);
  n = LGM_BNDP_N(theBndP);
  if (Bio_Write_mint(1,&n)) return (1);

  /* the lines */
  for (i=0; i<LGM_BNDP_NLINE(theBndP); i++)
  {
    n = LGM_LINE_ID(LGM_BNDP_LINE(theBndP,i));
    if (Bio_Write_mint(1,&n)) return (1);
                #ifdef NO_PROJECT
    global =  LGM_BNDP_LINE_GLOBALLEFT(theBndP,i);
    if (Bio_Write_mdouble(3,global)) return (1);
    global =  LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i);
    if (Bio_Write_mdouble(3,global)) return (1);
                #else
    e = LGM_BNDP_LINE_LEFT(theBndP,i);
    if (Bio_Write_mdouble(1,&e)) return (1);
    e = LGM_BNDP_LINE_RIGHT(theBndP,i);
    if (Bio_Write_mdouble(1,&e)) return (1);
                #endif
  }

  /* the surfaces */
  for (i=0; i<LGM_BNDP_N(theBndP); i++)
  {
    n = LGM_SURFACE_ID(LGM_BNDP_SURFACE(theBndP,i));
    if (Bio_Write_mint(1,&n)) return (1);
                #ifdef NO_PROJECT
    d[0] = LGM_BNDP_GLOBAL(theBndP,i)[0];
    d[1] = LGM_BNDP_GLOBAL(theBndP,i)[1];
    d[2] = LGM_BNDP_GLOBAL(theBndP,i)[2];
    if (Bio_Write_mdouble(3,d)) return (1);
                #else
    d[0] = LGM_BNDP_LOCAL(theBndP,i)[0];
    d[1] = LGM_BNDP_LOCAL(theBndP,i)[1];
    if (Bio_Write_mdouble(2,d)) return (1);
                #endif
  }
  return(0);
}

INT NS_DIM_PREFIX BNDP_SaveBndP_Ext (BNDP *aBndP)
{
  INT i;
  LGM_BNDP *theBndP;
  int n;
  double d[3], e;
#ifdef NO_PROJECT
  double *global;
#endif

  theBndP = BNDP2LGM(aBndP);
  /* the lines */
  n = LGM_BNDP_NLINE(theBndP);
  if (Bio_Write_mint(1,&n)) return (1);
  n = LGM_BNDP_N(theBndP);
  if (Bio_Write_mint(1,&n)) return (1);

  /* the lines */
  for (i=0; i<LGM_BNDP_NLINE(theBndP); i++)
  {
    /* store only an id, no line data */
    n = LGM_BNDP_ID(theBndP,i);
    if (Bio_Write_mint(1,&n)) return (1);
                #ifdef NO_PROJECT
    global =  LGM_BNDP_LINE_GLOBALLEFT(theBndP,i);
    if (Bio_Write_mdouble(3,global)) return (1);
    global =  LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i);
    if (Bio_Write_mdouble(3,global)) return (1);
                #else
    e = LGM_BNDP_LINE_LEFT(theBndP,i);
    if (Bio_Write_mdouble(1,&e)) return (1);
    e = LGM_BNDP_LINE_RIGHT(theBndP,i);
    if (Bio_Write_mdouble(1,&e)) return (1);
                #endif
  }

  /* the surfaces */
  for (i=0; i<LGM_BNDP_N(theBndP); i++)
  {
    n = LGM_BNDP_SURFACEID(theBndP,i);
    if (Bio_Write_mint(1,&n)) return (1);
                #ifdef NO_PROJECT
    d[0] = LGM_BNDP_GLOBAL(theBndP,i)[0];
    d[1] = LGM_BNDP_GLOBAL(theBndP,i)[1];
    d[2] = LGM_BNDP_GLOBAL(theBndP,i)[2];
    if (Bio_Write_mdouble(3,d)) return (1);
                #else
    d[0] = LGM_BNDP_LOCAL(theBndP,i)[0];
    d[1] = LGM_BNDP_LOCAL(theBndP,i)[1];
    if (Bio_Write_mdouble(2,d)) return (1);
                #endif
  }
  return(0);
}

/* domain interface function: for description see domain.h */
BNDP* NS_DIM_PREFIX BNDP_LoadBndP (BVP *theBVP, HEAP *Heap)
{
  LGM_DOMAIN *theDomain;
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  int i,id, nline, nsurf;
  double local[2], left, right;
  LGM_BNDP *theBndP;
#ifdef NO_PROJECT
  double global[3];
#endif

  theDomain = BVP2LGM(theBVP);

  if (Bio_Read_mint(1,&nline)) return (NULL);
  if (Bio_Read_mint(1,&nsurf)) return (NULL);

  theBndP = (LGM_BNDP*)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
  if(nline>0)
    theBndP->Line = (LGM_BNDP_PLINE*)GetFreelistMemory(Heap,nline*sizeof(LGM_BNDP_PLINE));
  else
    theBndP->Line = NULL;
  LGM_BNDP_NLINE(theBndP) = nline;

  theBndP->Surf = (LGM_BNDP_PSURFACE*)GetFreelistMemory(Heap,nsurf*sizeof(LGM_BNDP_PSURFACE));
  LGM_BNDP_N(theBndP) = nsurf;

  for (i=0; i<nline; i++)
  {
    if (Bio_Read_mint(1,&id)) return (NULL);
    for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
      if (LGM_LINE_ID(theLine)==id) break;
    if (theLine==NULL) return (NULL);

                #ifdef NO_PROJECT
    if (Bio_Read_mdouble(3,global)) return (NULL);
    LGM_BNDP_LINE_GLOBALLEFT(theBndP,i)[0] = global[0];
    LGM_BNDP_LINE_GLOBALLEFT(theBndP,i)[1] = global[1];
    LGM_BNDP_LINE_GLOBALLEFT(theBndP,i)[2] = global[2];
    if (Bio_Read_mdouble(3,global)) return (NULL);
    LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i)[0] = global[0];
    LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i)[1] = global[1];
    LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i)[2] = global[2];
    LGM_BNDP_LINE(theBndP,i) = theLine;
                #else
    if (Bio_Read_mdouble(1,&left)) return (NULL);
    if (Bio_Read_mdouble(1,&right)) return (NULL);
    LGM_BNDP_LINE(theBndP,i) = theLine;
    LGM_BNDP_LINE_LEFT(theBndP,i) = left;
    LGM_BNDP_LINE_RIGHT(theBndP,i) = right;
                #endif
  }

  for (i=0; i<nsurf; i++)
  {
    if (Bio_Read_mint(1,&id)) return (NULL);
    for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
      if (LGM_SURFACE_ID(theSurface)==id) break;
    if (theSurface==NULL) return (NULL);
                #ifdef NO_PROJECT
    if (Bio_Read_mdouble(3,global)) return (NULL);
    LGM_BNDP_SURFACE(theBndP,i) = theSurface;
    LGM_BNDP_GLOBAL(theBndP,i)[0] = global[0];
    LGM_BNDP_GLOBAL(theBndP,i)[1] = global[1];
    LGM_BNDP_GLOBAL(theBndP,i)[2] = global[2];
                #else
    if (Bio_Read_mdouble(2,local)) return (NULL);
    LGM_BNDP_SURFACE(theBndP,i) = theSurface;
    LGM_BNDP_LOCAL(theBndP,i)[0] = local[0];
    LGM_BNDP_LOCAL(theBndP,i)[1] = local[1];
                #endif
  }

  return((BNDP *)theBndP);
}

BNDP* NS_DIM_PREFIX BNDP_LoadBndP_Ext (void)
{
  int i,id, nline, nsurf;
  double local[2], left, right;
  LGM_BNDP *theBndP;
#ifdef NO_PROJECT
  double global[3];
#endif

  if (Bio_Read_mint(1,&nline)) return (NULL);
  if (Bio_Read_mint(1,&nsurf)) return (NULL);

  theBndP = (LGM_BNDP*)malloc(sizeof(LGM_BNDP));
  if(nline>0)
    theBndP->Line = (LGM_BNDP_PLINE*)malloc(nline*sizeof(LGM_BNDP_PLINE));
  else
    theBndP->Line = NULL;
  LGM_BNDP_NLINE(theBndP) = nline;

  theBndP->Surf = (LGM_BNDP_PSURFACE*)malloc(nsurf*sizeof(LGM_BNDP_PSURFACE));
  LGM_BNDP_N(theBndP) = nsurf;

  for (i=0; i<nline; i++)
  {
    if (Bio_Read_mint(1,&id)) return (NULL);
                #ifdef NO_PROJECT
    if (Bio_Read_mdouble(3,global)) return (NULL);
    LGM_BNDP_LINE_GLOBALLEFT(theBndP,i)[0] = global[0];
    LGM_BNDP_LINE_GLOBALLEFT(theBndP,i)[1] = global[1];
    LGM_BNDP_LINE_GLOBALLEFT(theBndP,i)[2] = global[2];
    if (Bio_Read_mdouble(3,global)) return (NULL);
    LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i)[0] = global[0];
    LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i)[1] = global[1];
    LGM_BNDP_LINE_GLOBALRIGHT(theBndP,i)[2] = global[2];
    LGM_BNDP_ID(theBndP,i) = id;
                #else
    if (Bio_Read_mdouble(1,&left)) return (NULL);
    if (Bio_Read_mdouble(1,&right)) return (NULL);
    LGM_BNDP_ID(theBndP,i) = id;
    LGM_BNDP_LINE_LEFT(theBndP,i) = left;
    LGM_BNDP_LINE_RIGHT(theBndP,i) = right;
                #endif
  }

  for (i=0; i<nsurf; i++)
  {
    if (Bio_Read_mint(1,&id)) return (NULL);
                #ifdef NO_PROJECT
    if (Bio_Read_mdouble(3,global)) return (NULL);
    LGM_BNDP_SURFACEID(theBndP,i) = id;
    LGM_BNDP_GLOBAL(theBndP,i)[0] = global[0];
    LGM_BNDP_GLOBAL(theBndP,i)[1] = global[1];
    LGM_BNDP_GLOBAL(theBndP,i)[2] = global[2];
                #else
    if (Bio_Read_mdouble(2,local)) return (NULL);
    LGM_BNDP_SURFACEID(theBndP,i) = id;
    LGM_BNDP_LOCAL(theBndP,i)[0] = local[0];
    LGM_BNDP_LOCAL(theBndP,i)[1] = local[1];
                #endif
  }

  return((BNDP *)theBndP);
}

INT NS_DIM_PREFIX LGM_BNDS_SurfId (BNDS *aBndS)
{
  int surf;
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;

  if (theBndS==NULL) return (-1);
  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);
  if (theSurface==NULL) return (-1);
  surf=LGM_SURFACE_ID(theSurface);

  return (surf);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDS_Global (BNDS *aBndS, DOUBLE *local, DOUBLE *global)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;
  INT mi;
  DOUBLE nv[3], loc0[2], loc1[2], loc2[2], loc3[2], global0[3], global1[3], global2[3], global3[3];

  /* global coordinates */
  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);

        #ifdef NO_PROJECT
  if (LGM_BNDS_N(theBndS)==3 || LGM_BNDS_N(theBndS)==-3)
  {
    global[0] =     (1.0 - local[0] - local[1] ) * LGM_BNDS_GLOBAL(theBndS,0,0)
                +   (      local[0]            ) * LGM_BNDS_GLOBAL(theBndS,1,0)
                +   (                 local[1] ) * LGM_BNDS_GLOBAL(theBndS,2,0);
    global[1] =     (1.0 - local[0] - local[1] ) * LGM_BNDS_GLOBAL(theBndS,0,1)
                +   (      local[0]            ) * LGM_BNDS_GLOBAL(theBndS,1,1)
                +   (                 local[1] ) * LGM_BNDS_GLOBAL(theBndS,2,1);
    global[2] =     (1.0 - local[0] - local[1] ) * LGM_BNDS_GLOBAL(theBndS,0,2)
                +   (      local[0]            ) * LGM_BNDS_GLOBAL(theBndS,1,2)
                +   (                 local[1] ) * LGM_BNDS_GLOBAL(theBndS,2,2);
  }
  else
  {
    global[0] =     (1.0 - local[0])*(1.0 - local[1])       * LGM_BNDS_GLOBAL(theBndS,0,0)
                +   local[0] * (1.0 - local[1])                 * LGM_BNDS_GLOBAL(theBndS,1,0)
                +   local[0] * local[1]                                 * LGM_BNDS_GLOBAL(theBndS,2,0)
                +   (1.0 - local[0])*local[1]                   * LGM_BNDS_GLOBAL(theBndS,3,0);
    global[1] =     (1.0 - local[0])*(1.0 - local[1])       * LGM_BNDS_GLOBAL(theBndS,0,1)
                +   local[0] * (1.0 - local[1])                 * LGM_BNDS_GLOBAL(theBndS,1,1)
                +   local[0] * local[1]                                 * LGM_BNDS_GLOBAL(theBndS,2,1)
                +   (1.0 - local[0])*local[1]                   * LGM_BNDS_GLOBAL(theBndS,3,1);
    global[2] =     (1.0 - local[0])*(1.0 - local[1])       * LGM_BNDS_GLOBAL(theBndS,0,2)
                +   local[0] * (1.0 - local[1])                 * LGM_BNDS_GLOBAL(theBndS,1,2)
                +   local[0] * local[1]                                 * LGM_BNDS_GLOBAL(theBndS,2,2)
                +   (1.0 - local[0])*local[1]                   * LGM_BNDS_GLOBAL(theBndS,3,2);
  }
        #else
  loc0[0] =  LGM_BNDS_LOCAL(theBndS,0,0);
  loc0[1] =  LGM_BNDS_LOCAL(theBndS,0,1);
  loc1[0] =  LGM_BNDS_LOCAL(theBndS,1,0);
  loc1[1] =  LGM_BNDS_LOCAL(theBndS,1,1);
  loc2[0] =  LGM_BNDS_LOCAL(theBndS,2,0);
  loc2[1] =  LGM_BNDS_LOCAL(theBndS,2,1);
  if (LGM_BNDS_N(theBndS)==4 || LGM_BNDS_N(theBndS)==-4) {
    loc3[0] =  LGM_BNDS_LOCAL(theBndS,3,0);
    loc3[1] =  LGM_BNDS_LOCAL(theBndS,3,1);
  }

  if (LGM_BNDS_N(theBndS)==3 || LGM_BNDS_N(theBndS)==-3) {
    Surface_Local2Global (theSurface,global0,loc0);
    Surface_Local2Global (theSurface,global1,loc1);
    Surface_Local2Global (theSurface,global2,loc2);
    global[0] =   (1.0 - local[0] - local[1] ) * global0[0]
                +   (      local[0]            ) * global1[0]
                +   (                 local[1] ) * global2[0];
    global[1] =   (1.0 - local[0] - local[1] ) * global0[1]
                +   (      local[0]            ) * global1[1]
                +   (                 local[1] ) * global2[1] ;
    global[2] =   (1.0 - local[0] - local[1] ) * global0[2]
                +   (      local[0]            ) * global1[2]
                +   (                 local[1] ) * global2[2];
  } else {
    Surface_Local2Global (theSurface,global0,loc0);
    Surface_Local2Global (theSurface,global1,loc1);
    Surface_Local2Global (theSurface,global2,loc2);
    Surface_Local2Global (theSurface,global3,loc3);
    global[0] =     (1.0 - local[0])*(1.0 - local[1]) * global0[0]
                +   local[0] * (1.0 - local[1])           * global1[0]
                +   local[0] * local[1]                           * global2[0]
                +   (1.0 - local[0])*local[1]                 * global3[0];
    global[1] =     (1.0 - local[0])*(1.0 - local[1]) * global0[1]
                +   local[0] * (1.0 - local[1])           * global1[1]
                +   local[0] * local[1]                           * global2[1]
                +   (1.0 - local[0])*local[1]                 * global3[1];
    global[2] =     (1.0 - local[0])*(1.0 - local[1]) * global0[2]
                +   local[0] * (1.0 - local[1])           * global1[2]
                +   local[0] * local[1]                           * global2[2]
                +   (1.0 - local[0])*local[1]             * global3[2];
  }

  V_DIM_CLEAR(nv);
  mi = GetLocalKoord(theSurface,global,local, nv);
  if(mi!=-1)
    Surface_Local2Global(theSurface,global,local);
  else
    assert(0);
        #endif

  return(0);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDS_BndCond (BNDS *aBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;
  DOUBLE global[DOM_PARAM_OFFSET];
  DOUBLE global0[3],global1[3],global2[3],global3[3],new_global[DIM+1];
  DOUBLE loc0[2],loc1[2],loc2[2],loc3[2],loc[2], nv[3];

  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);
  if (LGM_SURFACE_BNDCOND(theSurface)==NULL)
    return (2);

        #ifdef NO_PROJECT
  global0[0] = LGM_BNDS_GLOBAL(theBndS,0,0);
  global0[1] = LGM_BNDS_GLOBAL(theBndS,0,1);
  global0[2] = LGM_BNDS_GLOBAL(theBndS,0,2);
  global1[0] = LGM_BNDS_GLOBAL(theBndS,1,0);
  global1[1] = LGM_BNDS_GLOBAL(theBndS,1,1);
  global1[2] = LGM_BNDS_GLOBAL(theBndS,1,2);
  global2[0] = LGM_BNDS_GLOBAL(theBndS,2,0);
  global2[1] = LGM_BNDS_GLOBAL(theBndS,2,1);
  global2[2] = LGM_BNDS_GLOBAL(theBndS,2,2);
  if (LGM_BNDS_N(theBndS)==4 || LGM_BNDS_N(theBndS)==-4)
  {
    global3[0] = LGM_BNDS_GLOBAL(theBndS,3,0);
    global3[1] = LGM_BNDS_GLOBAL(theBndS,3,1);
    global3[2] = LGM_BNDS_GLOBAL(theBndS,3,2);
  }
        #else
  loc0[0] =  LGM_BNDS_LOCAL(theBndS,0,0);
  loc0[1] =  LGM_BNDS_LOCAL(theBndS,0,1);
  loc1[0] =  LGM_BNDS_LOCAL(theBndS,1,0);
  loc1[1] =  LGM_BNDS_LOCAL(theBndS,1,1);
  loc2[0] =  LGM_BNDS_LOCAL(theBndS,2,0);
  loc2[1] =  LGM_BNDS_LOCAL(theBndS,2,1);
  if (LGM_BNDS_N(theBndS)==4 || LGM_BNDS_N(theBndS)==-4) {
    loc3[0] =  LGM_BNDS_LOCAL(theBndS,3,0);
    loc3[1] =  LGM_BNDS_LOCAL(theBndS,3,1);
  }

  Surface_Local2Global (theSurface,global0,loc0);
  Surface_Local2Global (theSurface,global1,loc1);
  Surface_Local2Global (theSurface,global2,loc2);
  if (LGM_BNDS_N(theBndS)==4 || LGM_BNDS_N(theBndS)==-4)
    Surface_Local2Global (theSurface,global3,loc3);
        #endif
  if (LGM_BNDS_N(theBndS)==3 || LGM_BNDS_N(theBndS)==-3)
  {
    global[0] = (1-local[0]-local[1]) * global0[0] + local[0] * global1[0] + local[1] * global2[0];
    global[1] = (1-local[0]-local[1]) * global0[1] + local[0] * global1[1] + local[1] * global2[1];
    global[2] = (1-local[0]-local[1]) * global0[2] + local[0] * global1[2] + local[1] * global2[2];
  }
  else
  {
    global[0] = (1-local[0])*(1.0-local[1]) * global0[0] +
                local[0]*(1.0-local[1])         * global1[0] +
                local[0]*local[1]           * global2[0] +
                (1-local[0])*local[1]           * global3[0];
    global[1] = (1-local[0])*(1.0-local[1]) * global0[1] +
                local[0]*(1.0-local[1])         * global1[1] +
                local[0]*local[1]           * global2[1] +
                (1-local[0])*local[1]           * global3[1];
    global[2] = (1-local[0])*(1.0-local[1]) * global0[2] +
                local[0]*(1.0-local[1])         * global1[2] +
                local[0]*local[1]           * global2[2] +
                (1-local[0])*local[1]           * global3[2];
  }


        #ifdef NO_PROJECT
  new_global[0] = global[0];
  new_global[1] = global[1];
  new_global[2] = global[2];
        #else
  V_DIM_CLEAR(nv);
  GetLocalKoord(theSurface,global,loc, nv);
  Surface_Local2Global(theSurface, new_global, loc);
        #endif

  if (in!=NULL)
  {
    in[0] = new_global[0];
    in[1] = new_global[1];
    in[2] = new_global[2];
    in[DIM] = LGM_SURFACE_ID(theSurface);
    if ((*LGM_SURFACE_BNDCOND (theSurface))(in,value,type))
      return (1);
  }
  else
  {
    new_global[DIM] = LGM_SURFACE_ID(theSurface);
    if ((*LGM_SURFACE_BNDCOND (theSurface))(new_global,value,type))
      return (1);
  }

  return (0);
}

/* domain interface function: for description see domain.h */
INT NS_DIM_PREFIX BNDS_BndSDesc (BNDS *aBndS, INT *left, INT *right, INT *part)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;

  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);

  part[0] = 0;
  if(LGM_BNDS_N(theBndS)<0)
  {
    *left = LGM_SURFACE_LEFT(theSurface);
    *right = LGM_SURFACE_RIGHT(theSurface);
  }
  else
  {
    *right = LGM_SURFACE_LEFT(theSurface);
    *left = LGM_SURFACE_RIGHT(theSurface);
  }

  /* HRR_TODO: assign part */
  *part = 0;

  return(0);
}

/* domain interface function: for description see domain.h */
BNDP* NS_DIM_PREFIX BNDS_CreateBndP (HEAP *Heap, BNDS *aBndS, DOUBLE *local)
{
  LGM_BNDS *theBndS;
  LGM_BNDP *theBndP;
  LGM_SURFACE *theSurface;
  DOUBLE loc0[2],loc1[2],loc2[2],loc3[2],loc[2];
  DOUBLE global0[3],global1[3],global2[3],global3[3],global[3], nv[3];

  if (local[0]<=0.0 || local[0]>=1.0)
    return (NULL);
  if (local[1]<=0.0 || local[1]>=1.0)
    return (NULL);

  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);

        #ifdef NO_PROJECT
  global0[0] = LGM_BNDS_GLOBAL(theBndS,0,0);
  global0[1] = LGM_BNDS_GLOBAL(theBndS,0,1);
  global0[2] = LGM_BNDS_GLOBAL(theBndS,0,2);
  global1[0] = LGM_BNDS_GLOBAL(theBndS,1,0);
  global1[1] = LGM_BNDS_GLOBAL(theBndS,1,1);
  global1[2] = LGM_BNDS_GLOBAL(theBndS,1,2);
  global2[0] = LGM_BNDS_GLOBAL(theBndS,2,0);
  global2[1] = LGM_BNDS_GLOBAL(theBndS,2,1);
  global2[2] = LGM_BNDS_GLOBAL(theBndS,2,2);
  if (LGM_BNDS_N(theBndS)==4 || LGM_BNDS_N(theBndS)==-4)
  {
    global3[0] = LGM_BNDS_GLOBAL(theBndS,3,0);
    global3[1] = LGM_BNDS_GLOBAL(theBndS,3,1);
    global3[2] = LGM_BNDS_GLOBAL(theBndS,3,2);
  }
        #else
  loc0[0] =  LGM_BNDS_LOCAL(theBndS,0,0);
  loc0[1] =  LGM_BNDS_LOCAL(theBndS,0,1);
  loc1[0] =  LGM_BNDS_LOCAL(theBndS,1,0);
  loc1[1] =  LGM_BNDS_LOCAL(theBndS,1,1);
  loc2[0] =  LGM_BNDS_LOCAL(theBndS,2,0);
  loc2[1] =  LGM_BNDS_LOCAL(theBndS,2,1);
  if (LGM_BNDS_N(theBndS)==4 || LGM_BNDS_N(theBndS)==-4)
  {
    loc3[0] = LGM_BNDS_LOCAL(theBndS,3,0);
    loc3[1] = LGM_BNDS_LOCAL(theBndS,3,1);
  }
  Surface_Local2Global (theSurface,global0,loc0);
  Surface_Local2Global (theSurface,global1,loc1);
  Surface_Local2Global (theSurface,global2,loc2);
  if (LGM_BNDS_N(theBndS)==4 || LGM_BNDS_N(theBndS)==-4)
    Surface_Local2Global (theSurface,global3,loc3);
        #endif

  if (LGM_BNDS_N(theBndS)==3 || LGM_BNDS_N(theBndS)==-3)
  {
    global[0] = (1-local[0]-local[1]) * global0[0] + local[0] * global1[0] + local[1] * global2[0];
    global[1] = (1-local[0]-local[1]) * global0[1] + local[0] * global1[1] + local[1] * global2[1];
    global[2] = (1-local[0]-local[1]) * global0[2] + local[0] * global1[2] + local[1] * global2[2];
  }
  else
  {
    global[0] = (1-local[0])*(1.0-local[1]) * global0[0] +
                local[0]*(1.0-local[1])         * global1[0] +
                local[0]*local[1]           * global2[0] +
                (1-local[0])*local[1]           * global3[0];
    global[1] = (1-local[0])*(1.0-local[1]) * global0[1] +
                local[0]*(1.0-local[1])         * global1[1] +
                local[0]*local[1]           * global2[1] +
                (1-local[0])*local[1]           * global3[1];
    global[2] = (1-local[0])*(1.0-local[1]) * global0[2] +
                local[0]*(1.0-local[1])         * global1[2] +
                local[0]*local[1]           * global2[2] +
                (1-local[0])*local[1]           * global3[2];
  }

  theBndP = (LGM_BNDP *)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
  LGM_BNDP_NLINE(theBndP) = 0;
  theBndP->Surf = (LGM_BNDP_PSURFACE*)GetFreelistMemory(Heap,sizeof(LGM_BNDP_PSURFACE));
  LGM_BNDP_N(theBndP) = 1;
  LGM_BNDP_SURFACE(theBndP,0) = theSurface;

        #ifdef NO_PROJECT
  LGM_BNDP_GLOBAL(theBndP,0)[0] = global[0];
  LGM_BNDP_GLOBAL(theBndP,0)[1] = global[1];
  LGM_BNDP_GLOBAL(theBndP,0)[2] = global[2];
        #else
  V_DIM_CLEAR(nv);
  GetLocalKoord(theSurface,global,loc, nv);
  LGM_BNDP_LOCAL(theBndP,0)[0] = loc[0];                        /* TODO */
  LGM_BNDP_LOCAL(theBndP,0)[1] = loc[1];
        #endif

  return((BNDP *)theBndP);
}



/* auxiliary function for getting the maximum surface-ID */
INT NS_DIM_PREFIX GetMaximumSurfaceID (LGM_DOMAIN *theDomain)
{
  INT nSubDom, i, l, maxLineId=0;

  nSubDom = LGM_DOMAIN_NSUBDOM(theDomain);
  for(i=1; i<=nSubDom; i++)
  {
    LGM_SUBDOMAIN *subdom = LGM_DOMAIN_SUBDOM(theDomain,i);

    for(l=0; l<LGM_SUBDOMAIN_NSURFACE(subdom); l++)
    {
      INT id = LGM_SURFACE_ID(LGM_SUBDOMAIN_SURFACE(subdom,l));
      if (maxLineId < id)
        maxLineId = id;
    }
  }

  return(maxLineId);
}

/* OS_CHANGED: an new auxiliary function to mark surface_IDs of the domain boundary */
INT NS_DIM_PREFIX OuterBndSurfaceIDs (LGM_DOMAIN *theDomain, INT *sf)
{
  INT nSubDom, i, l;

  nSubDom = LGM_DOMAIN_NSUBDOM(theDomain);
  for(i=1; i<=nSubDom; i++)
  {
    LGM_SUBDOMAIN *subdom = LGM_DOMAIN_SUBDOM(theDomain,i);

    for(l=0; l<LGM_SUBDOMAIN_NSURFACE(subdom); l++)
    {
      INT id = LGM_SURFACE_ID(LGM_SUBDOMAIN_SURFACE(subdom,l));
      if ((LGM_SURFACE_LEFT(LGM_SUBDOMAIN_SURFACE(subdom,l)) == 0) ||
          (LGM_SURFACE_RIGHT(LGM_SUBDOMAIN_SURFACE(subdom,l)) == 0))
        sf[id] = 1;
      else
        sf[id] = 0;
    }
  }

  return 0;
}

/* OS_CHANGED: an new auxiliary function to get surface_IDs of a subdomain,
 * return value is the number of surfaces of the subdomain i
 */
INT SurfaceIDsOfSubdomain (LGM_DOMAIN *theDomain, INT *sf, INT i)
{
  INT nsf, l;

  LGM_SUBDOMAIN *subdom;
  if (i>LGM_DOMAIN_NSUBDOM(theDomain)) return 0;

  subdom = LGM_DOMAIN_SUBDOM(theDomain,i);
  nsf = LGM_SUBDOMAIN_NSURFACE(subdom);
  for(l=0; l<nsf; l++)
    sf[l] = LGM_SURFACE_ID(LGM_SUBDOMAIN_SURFACE(subdom,l));

  return nsf;
}
