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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compiler.h"
#include "heaps.h"
#include "bio.h"
#include "domain.h"
#include "general.h"
#include "misc.h"
#include "lgm_domain.h"
#include "lgm_load.h"
#include "lgm_macros.h"
#include "lgm_gginterface.h"

#include "devices.h"

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

static char buffer[LGM_BUFFERLEN];

static INT SurfaceInfoId;
static INT LineInfoId;
static INT LGM_DEBUG = 0;
static SMALL = 0.001;
static INT VAR_H = 1;
static double EPS = 0.00001;

static DOUBLE cosAngle = 0.94;          /* Winkel zwischen Inputdreiecken */

INT Surface_Local2Global (LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local);
INT GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local);

static INT Get_NBNDS_Per_Subdomain              (HEAP *Heap, LGM_DOMAIN *theDomain, INT **sizes, DOUBLE h);
static INT DiscretizeLine                               (HEAP *Heap, LGM_LINE *theLine, DOUBLE h, LGM_POINT *pointlist, INT norp);
static INT Get_NBNDP                                    (LGM_DOMAIN *theDomain, INT *nBND, DOUBLE h);
static INT DiscretizeDomain                     (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h);
static INT DiscretizeSurface                    (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h, LGM_POINT *pointlist, INT norp);
static INT TransferLines2Mesh                   (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h);
static INT TransferSurfaces2Mesh                (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h);

extern INT GenerateSurfaceGrid (LGM_SURFACE *aSurface, DOUBLE h, INT smooth,INT display);
extern INT InitSurface(CoeffProcPtr Coeff);
extern int AddGeomElement (int node0, int node1, int node2, int neigbor0, int neigbor1, int neigbor2);

/* data for CVS */
/*static char RCS_ID("$Header$",UG_RCS_STRING);
 */
static INT currSubdom, currSurface, currLine;

static CoeffProcPtr Coefficients[8];
static CoeffProcPtr LOCAL_H;

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

static LGM_SURFACE *FirstSurface (LGM_DOMAIN *theDomain)
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

static LGM_SURFACE *NextSurface (LGM_DOMAIN *theDomain)
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

static LGM_LINE *FirstLine (LGM_DOMAIN *theDomain)
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

/*static LGM_LINE *helpNextLine (LGM_DOMAIN *theDomain)
   {
        INT i,j,k;
        LGM_SUBDOMAIN *theSubdom;
        LGM_SURFACE *theSurface;

        for (i=1; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
        {
                theSubdom = LGM_DOMAIN_SUBDOM(theDomain,i);
                for(j=0;j<LGM_SUBDOMAIN_NSURFACE(theSubdom);j++)
                {
                        theSurface = LGM_SUBDOMAIN_SURFACE(theSubdom,j);
                        for(k=0;k<LGM_SURFACE_NLINE(theSurface);k++)
                        {
                                if(LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface,k))==0)
                                {
                                        LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface,k))==1;
                                        return(LGM_SURFACE_LINE(theSurface,k));
                                }
                        }
                }
        }

        return (NULL);
   }*/

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

static LGM_LINE *NextLine (LGM_DOMAIN *theDomain)
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

INT SetBoundaryCondition (LGM_DOMAIN *theDomain, BndCondProcPtr BndCond)
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
        LGM_SURFACE_BNDCOND(theSurface) = NULL;
      else
        LGM_SURFACE_BNDCOND(theSurface) = BndCond;
    }
  }

  return (0);
}

INT SetDomainSize (LGM_DOMAIN *theDomain)
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

  return (0);
}

static INT PrintDomainInfo (LGM_DOMAIN *aDomain)
{

  printf("********* domain-info *********\n");

  printf("%s %s\n","Name: ",LGM_DOMAIN_PROBLEMNAME(aDomain));
  printf("%s %d\n","nSubdomain: ",LGM_DOMAIN_NSUBDOM(aDomain));
  printf("%s %d\n","nPoint: ",LGM_DOMAIN_NPOINT(aDomain));
  printf("%s %f\n","radius: ",LGM_DOMAIN_RADIUS(aDomain));
  printf("%s %f %f %f\n","midpoint: ",LGM_DOMAIN_MIDPOINT(aDomain)[0],
         LGM_DOMAIN_MIDPOINT(aDomain)[1],
         LGM_DOMAIN_MIDPOINT(aDomain)[2]);
  printf("%s %d\n","convex: ",LGM_DOMAIN_CONVEX(aDomain));

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
    printf("%s %d\n","id: ",LGM_SUBDOMAIN_ID(aSubdom));
    printf("%s %d\n","nSurface: ",LGM_SUBDOMAIN_NSURFACE(aSubdom));
  }
  return (0);
}

static INT PrintSurfaceInfo (LGM_SURFACE *aSurface)
{
  INT i;

  printf("********* surface-info *********\n");
  printf("%s %d\n","SurfaceId: ",SurfaceInfoId);
  SurfaceInfoId++;
  printf("%s %d\n","nPoint: ",LGM_SURFACE_NPOINT(aSurface));
  printf("%s %d\n","nTriangle: ",LGM_SURFACE_NTRIANGLE(aSurface));
  printf("%s %d\n","nLine: ",LGM_SURFACE_NLINE(aSurface));
  printf("%s %d\n","left: ",LGM_SURFACE_LEFT(aSurface));
  printf("%s %d\n","right: ",LGM_SURFACE_RIGHT(aSurface));
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
  /*	for(i=0;i<LGM_SURFACE_NTRIANGLE(aSurface);i++)
                  printf("%s %d %d %d\n","Neighbor: ",LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(aSurface,i),0),
                  LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(aSurface,i),1),
                  LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(aSurface,i),2));*/

  for(i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(aSurface)); i++)
  {
    printf("%d %d %d\n",LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(aSurface),i,0),
           LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(aSurface),i,1),
           LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(aSurface),i,2));
  }
  /*	for(i=0;i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(aSurface));i++)
          {
                  local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(aSurface),i,0);
                  local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(aSurface),i,1);
                  printf("%f %f\n",local[0],local[1]);
                  Surface_Local2Global(aSurface,global,local);
                  printf("%f %f %f\n",global[0],global[1],global[2]);
          }*/


  return (0);
}

INT Line_Local2Global (LGM_LINE *theLine, DOUBLE *global, INT i)
{
  DOUBLE slocal;
  INT ilocal;

  ilocal = floor(LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),i));
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

static INT PrintLineInfo (LGM_LINE *aLine)
{
  INT i;
  DOUBLE global[3];

  printf("********* line-info *********\n");
  printf("%s %d\n","LineId: ",LGM_LINE_ID(aLine));
  LineInfoId++;
  printf("%s %d\n","nPoint: ",LGM_LINE_NPOINT(aLine));
  printf("%s %d %d\n","firstPoint lastPoint: ",LGM_LINE_BEGIN(aLine),LGM_LINE_END(aLine));
  for(i=0; i<LGM_LINE_NPOINT(aLine); i++)
    printf("%s %f %f %f\n","Point: ",LGM_LINE_POINT(aLine,i)->position[0],
           LGM_LINE_POINT(aLine,i)->position[1],
           LGM_LINE_POINT(aLine,i)->position[2]);
  printf("%s\n","linedisc");
  printf("%s %d\n","nPoint: ",LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(aLine)));
  for(i=0; i<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(aLine)); i++)
  {
    printf("%f\n",LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(aLine),i));
    Line_Local2Global(aLine,global,i);
    printf("%f %f %f\n",global[0],global[1],global[2]);
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
    printf("\n");
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
INT BVP_Save (BVP *theBVP, char *name, char *mgname, HEAP *theHeap, INT argc, char **argv)
{
  UserWrite("SORRY: not implemented yet\n");
  return (1);
}

/* domain interface function: for description see domain.h */
INT BVP_Check (BVP *aBVP)
{
  INT i,i2,j,j2,k,p,ret,at_begin,at_left,flags,sbd,sfce,left,right,l,found,m, m2;
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
        theSF = LGM_SUBDOMAIN_SURFACE(theSD,j);
        /*Laufe ueber Points*/

        lgmd = (LGM_DOMAIN*)aBVP;
        MarkTmpMem(lgmd->theHeap);

        memsize_variable = LGM_SURFACE_NPOINT(theSF)*sizeof(LGM_POINT*);

        if ((TemporaryPointArray = GetTmpMem(lgmd->theHeap,memsize_variable))==NULL)
        {
          PrintErrorMessage('E',"BVP_Check","  ERROR: No memory for TemporaryPointArray");
          return(1);
        }
        /* Mustermemset(nodeflag_array,0,(statistik[0]+1)*sizeof(INT)); */
        memset(TemporaryPointArray,NULL,memsize_variable);


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

        ReleaseTmpMem(lgmd->theHeap);

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
                /*Die beiden gemeinsamen Corneer muessen in unterschiedlicher Reihenfolge vorkommen*/

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

  l1 = Lenght(n1);
  Scale(n1, 1/l1);

  n2[0] = p2[0] - p1[0];
  n2[1] = p2[1] - p1[1];
  n2[2] = p2[2] - p1[2];

  l2 = Lenght(n2);
  Scale(n2, 1/l2);

  Cross(n, n1, n2);
  l = Lenght(n);
  Scale(n, 1/l);

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

  scalarproduct =  Mult(n1, n2);
  return(scalarproduct);
}

static INT Count_Folds_Surface(LGM_SURFACE *theSurface)
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
          folds++;
      }
    }
  }
  return(folds);
}

INT POINT_DIST(LGM_POINT *p0, LGM_POINT *p1)
{
  DOUBLE x, y, z;

  x = LGM_POINT_POS(p0)[0] - LGM_POINT_POS(p1)[0];
  y = LGM_POINT_POS(p0)[1] - LGM_POINT_POS(p1)[1];
  z = LGM_POINT_POS(p0)[2] - LGM_POINT_POS(p1)[2];
  if( (x<EPS) && (x>-EPS) &&
      (y<EPS) && (y>-EPS) &&
      (z<EPS) && (z>-EPS) )
    return(1);
  else
    return(0);
}

INT AddPoint2List(LGM_POINT *Point, LGM_POINT *pointlist, INT *norp)
{
  int i, flag;
  DOUBLE x1, y1, z1;
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
}

INT ResolvePoints(HEAP *Heap, LGM_DOMAIN *theDomain, LGM_POINT *pointlist, INT *norp)
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

/* domain interface function: for description see domain.h */
MESH *BVP_GenerateMesh (HEAP *Heap, BVP *aBVP, INT argc, char **argv)
{
  LGM_DOMAIN *theDomain;
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  MESH *mesh;
  LGM_POINT *pointlist;
  INT i,j;
  INT *BNDS_Per_Subdom, *p;
  float fValue;
  int iValue;
  DOUBLE h;
  INT coeff;
  INT old;
  INT npoints,npsurface, norp;

  old = 0;
  /* read h-option */
  h = 0.0;
  for (i=1; i<argc; i++)
    if (argv[i][0] == 'h')
    {
      if (sscanf(argv[i],"h %f",&fValue) != 1)
        return (NULL);
      h = fValue;
    }
  coeff = 0;
  for (i=1; i<argc; i++)
    if (argv[i][0] == 'c')
    {
      if (sscanf(argv[i],"c %d",&iValue) != 1)
        return (NULL);
      coeff = iValue;
    }
  /*	if (h<=0.0) return (NULL);*/

  if(h<=0.0)
  {
    /* get Coefficientfunctions */
    if (BVP_SetCoeffFct(aBVP,-1,Coefficients))
      return (NULL);
    LOCAL_H = Coefficients[coeff];
  }

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

  /* allocate mesh */
  mesh = (MESH *) GetTmpMem(Heap,sizeof(MESH));
  if (mesh == NULL)
    return(NULL);

  /* init mesh: only surface-mesh */
  mesh->nInnP = 0;
  mesh->Position = NULL;
  mesh->nElements = NULL;
  mesh->Element_corners = NULL;
  mesh->Element_corner_ids = NULL;
  mesh->nSubDomains = LGM_DOMAIN_NSUBDOM(theDomain);

  /* get heap for surface-mesh-substructures: the subdomain-dependence */
  mesh->nSides = (INT *) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT));
  if (mesh->nSides==NULL)
    return(NULL);
  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
    mesh->nSides[i] = 0;


  mesh->Side_corners = (INT **) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT*));
  if (mesh->Side_corners == NULL)
    return (NULL);
  mesh->Side_corner_ids = (INT ***) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT**));
  if (mesh->Side_corner_ids == NULL)
    return (NULL);

  /*	if (Get_NBNDS_Per_Subdomain(Heap,theDomain,&BNDS_Per_Subdom,h))
                  return (NULL);
          for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
          {
          mesh->Side_corners[i]    = (INT *)  GetTmpMem(Heap,sizeof(INT)*BNDS_Per_Subdom[i]);
          mesh->Side_corner_ids[i] = (INT **) GetTmpMem(Heap,sizeof(INT*)*BNDS_Per_Subdom[i]);
                  p = (INT *) GetTmpMem(Heap,3*sizeof(INT)*BNDS_Per_Subdom[i]);
                  for (j=0; j<BNDS_Per_Subdom[i]; j++)
                          mesh->Side_corner_ids[i][j] = p+3*j;
          }*/

  /*	if (Get_NBNDP(theDomain,&nBNDP_on_Domain,h))
                  return (NULL);
          mesh->theBndPs = (BNDP**)GetTmpMem(Heap,sizeof(LGM_BNDP*)*nBNDP_on_Domain);

          if (mesh->theBndPs == NULL)
                  return (NULL);*/

  /* search points, which must be resolved by DiscretizeLine and DiscretizeSurface */
  norp = 0;
  pointlist = (LGM_POINT *)GetTmpMem(Heap,100*sizeof(LGM_POINT));

  ResolvePoints(Heap, theDomain, pointlist, &norp);

  if(LGM_DEBUG)
    for(i=0; i<norp; i++)
      printf("%d %f %f %f\n", i, pointlist[i].position[0], pointlist[i].position[1],  pointlist[i].position[2]);


  /* prepare for surface-mesh */
  mesh->nBndP = 0;

  /* discretize lines */
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    theLine->ldisc = (LGM_LINEDISC*)GetTmpMem(Heap,sizeof(LGM_LINEDISC)*1);
    if( (theLine->ldisc==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }
    if (DiscretizeLine(Heap,theLine,h,pointlist, norp))
      return(NULL);
    if (LGM_DEBUG)
      PrintLineInfo(theLine);
  }

  /* discretize surfaces */
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    LGM_SURFACE_DISC(theSurface) = (LGM_SURFDISC*)GetTmpMem(Heap,sizeof(LGM_SURFDISC)*1);
    if( (LGM_SURFACE_DISC(theSurface)==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }
    if (LGM_DEBUG)
      for(j=0; j<LGM_SURFACE_NLINE(theSurface); j++)
      {
        PrintLineInfo(LGM_SURFACE_LINE(theSurface,j));
      }
    if (DiscretizeSurface(Heap,theSurface,mesh,h,pointlist,  norp))
      return(NULL);
    if (LGM_DEBUG)
    {
      PrintSurfaceInfo(theSurface);
      printf("%s\n", "#####################################################");
    }
  }

  /* discretize domain */
  if (DiscretizeDomain(Heap,theDomain,mesh,h))
    return(NULL);

  /*	if (LGM_DEBUG)
                  if (PrintMeshInfo(mesh))
                          return (NULL);*/

  /* discretize lines */
  /*	if (TransferLines2Mesh(Heap,theDomain,mesh,h))
                  return(NULL);

          if (LGM_DEBUG)
                  if (PrintMeshInfo(mesh))
                          return (NULL);*/

  if (Get_NBNDS_Per_Subdomain(Heap,theDomain,&BNDS_Per_Subdom,h))
    return (NULL);

  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    mesh->Side_corners[i]    = (INT *)  GetTmpMem(Heap,sizeof(INT)*BNDS_Per_Subdom[i]);
    if (mesh->Side_corners[i]==NULL)
      return (NULL);
    mesh->Side_corner_ids[i] = (INT **) GetTmpMem(Heap,sizeof(INT*)*BNDS_Per_Subdom[i]);
    if (mesh->Side_corner_ids[i]==NULL)
      return (NULL);
    p = (INT *) GetTmpMem(Heap,3*sizeof(INT)*BNDS_Per_Subdom[i]);
    for (j=0; j<BNDS_Per_Subdom[i]; j++)
      mesh->Side_corner_ids[i][j] = p+3*j;
  }

  /* discretize surfaces */
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    if (TransferSurfaces2Mesh(Heap,theSurface,mesh,h))
      return(NULL);
    /*		if (LGM_DEBUG)
                            if (PrintMeshInfo(mesh))
                                    return (NULL);*/
  }

  mesh->VertexLevel = NULL;
  mesh->VertexPrio = NULL;

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
    innerpoints = innerpoints - ( LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,j))) - 1);
  }
  /*	*n = LGM_SURFACE_NPOINT(theSurface);*/
  *n = innerpoints;
  return (0);
}

static INT Get_NBNDS_Per_Subdomain (HEAP *Heap, LGM_DOMAIN *theDomain, INT **sizes, DOUBLE h)
{
  INT i,*p,n;
  LGM_SURFACE *theSurface;

  /* get heap for information */
  p = (INT *) GetTmpMem(Heap,sizeof(INT)*(LGM_DOMAIN_NSUBDOM(theDomain)+1));
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

  /*	n = LGM_DOMAIN_NPOINT(theDomain);
          for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
          {
                  for (i=0; i<LGM_SURFACE_NPOINT(theSurface)-1; i++)
                  {
                          if (InnerPointsPerSurfaceSegment(theSurface,h,i,&pn))
                                  return (1);
                          n += pn;
                  }
          }*/

  n = LGM_DOMAIN_NPOINT(theDomain);
  /*	for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
          {
                  n += LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine)) - 2;
          }*/
  *nBND = n;

  return (0);
}

static DOUBLE Calc_Local_Coord(DOUBLE *p0, DOUBLE*p1, DOUBLE *p2, DOUBLE *global, DOUBLE *lam)
{
  INT i,j;
  DOUBLE small;
  DOUBLE e0[3],e1[3],e2[3];
  DOUBLE n0[3],n1[3],n2[3],p[3],np[3],hp[3],l;
  DOUBLE a[9],b[3],c[3];
  DOUBLE aa[4],bb[2],cc[2];

  small = 0.00001;

  Minus(e0,global,p0);
  Minus(e1,global,p1);
  Minus(e2,global,p2);

  /* compute 2 vectors in the plane + normalvector */
  Minus(n0,e2,e0);
  l = Lenght(n0);
  Scale(n0,1/l);

  Minus(n1,e2,e1);
  l = Lenght(n1);
  Scale(n1,1/l);

  Cross(n2,n0,n1);
  l = Lenght(n2);
  Scale(n2,1/l);

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

  if((Det2d(aa)>small)||(Det2d(aa)<-small))
  /*	if(Det2d(aa)!=0)*/
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
    if((Det2d(aa)>small)||(Det2d(aa)<-small))
    /*		if(Det2d(aa)!=0)*/
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
      if((Det2d(aa)>small)||(Det2d(aa)<-small))
      /*			if(Det2d(aa)!=0)*/
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

  return(Lenght(hp));
}

static INT CASE1(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *lam, DOUBLE *d)
{
  DOUBLE small, dist;
  INT i, mi;
  DOUBLE *p0,*p1,*p2;

  /* Punkt liegt auf der Surface
   * 0<=lamda_i<=1,  dist < eps
   * Falls der Punkt genau auf einer Kante zwischen 2 Dreiecken liegt,
   * ist es egal welches genommen wird
   */
  small = 0.000001;
  dist = 10000000.0;
  mi = -1;
  d[0] = 10000000.0;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

    dist = Calc_Local_Coord(p0, p1, p2, global, lam);
    if( (lam[0]>=-small) && (lam[1]>=-small) && (lam[2]>=-small) && (dist<small) )
    {
      mi = i;
      d[0] = dist;
      break;
    }
  }
  return(mi);
}

static INT CASE2(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *lam)
{
  DOUBLE small, dist1, dist2, min, d1, d2;
  INT i, j, mi1, mi2, mi;
  DOUBLE *p0,*p1,*p2, mm;
  DOUBLE A, B, m, point[3], dist_vec[3],pp0[3], pp1[3], new_global[3], local[3];

  /* Punkt liegt nicht auf der Surface
   * Projeziere den Punkt in ein Dreieck mit  0<=lamda_i<=1 und dist minimal
   * oder bestimme einen Punkt auf einer Kante,  sodass der Abstand
   * noch kleiner wird
   */

  /* kuerzeste Abstand zur einem Dreieck */
  min = 10000000.0;
  mi1 = -1;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

    dist1 = Calc_Local_Coord(p0, p1, p2, global, lam);
    if( (lam[0]>=-small) && (lam[1]>=-small) && (lam[2]>=-small) && (dist1<min) )
    {
      mi1 = i;
      min = dist1;
    }
  }
  dist1 = min;

  /* kuerzeste Abstand zur eine Kante */
  min = 10000000.0;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
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
  dist2 = min;

  if(dist1<=dist2)
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi1),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi1),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi1),2);

    Calc_Local_Coord(p0, p1, p2, global, lam);
    return(mi1);
  }
  else
  {
    mi = GetLocalKoord(theSurface, new_global, local);
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,mi),2);

    Calc_Local_Coord(p0, p1, p2, new_global, lam);
    return(mi);
  }
}
INT GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  INT i,j,test,min,mi, dist_i;
  DOUBLE *p0,*p1,*p2,e0[3],e1[3],e2[3],eps;
  DOUBLE n0[3],n1[3],n2[3],p[3],np[3],hp[3],l, d;
  DOUBLE a[9],b[3],c[3];
  DOUBLE aa[4],bb[2],cc[2];
  DOUBLE lam[3],small, dist, min_dist, new_lam[3], A, B, m, point[3], dist_vec[3],pp0[3], pp1[3], new_global[3];
  min = 10000000.0;
  min_dist = 10000000.0;

  dist_i = -1;
  mi = -1;
  small = 0.000001;

  mi = CASE1(theSurface, global, lam, &d);
  if(mi==-1)
    mi = CASE2(theSurface, global, lam);
  if(mi==-1)
    printf("%s\n", "schotter");

  if( (lam[0]<0.0) || (lam[1]<0.0) )
  {
    if( (lam[0]<0.0) && (lam[0]>-0.000001) )
      lam[0] = 0.0;
    if( (lam[1]<0.0) && (lam[1]>-0.000001) )
      lam[1] = 0.0;
  }

  local[0] = lam[0] + mi;
  local[1] = lam[1] + mi;

  /*	printf("%f %f %d %f %f %f\n", local[0], local[1], mi, global[0], global[1], global[2]);*/
  /*	Surface_Local2Global(theSurface, new_global, local);

          if(LGM_DEBUG)
          {
                  printf("%f %f %f\n", global[0], global[1], global[2]);
                  printf("%f %f %f\n", new_global[0], new_global[1], new_global[2]);
          }*/
  return(0);
}

DOUBLE Check_Surface(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local, DOUBLE *d)
{
  DOUBLE small, dist;
  INT i, mi;
  DOUBLE *p0,*p1,*p2, lam[3], lam1[3];

  small = 0.000001;
  dist = 10000000.0;
  mi = -1;
  d[0] = 10000000.0;

  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

    dist = Calc_Local_Coord(p0, p1, p2, global, lam1);
    if( (lam1[0]>=-small) && (lam1[1]>=-small) && (lam1[2]>=-small) && (d[0]>dist))
    {
      mi = i;
      d[0] = dist;
      lam[0] = lam1[0];
      lam[1] = lam1[1];
      lam[2] = lam1[2];
    }
  }

  /* if mi<0, no triangle has been found. hence, lam[] isn't valid! */
  if (mi==-1)
    return(1);

  if( (lam[0]<0.0) || (lam[1]<0.0) )
  {
    if( (lam[0]<0.0) && (lam[0]>-0.000001) )
      lam[0] = 0.0;
    if( (lam[1]<0.0) && (lam[1]>-0.000001) )
      lam[1] = 0.0;
  }

  local[0] = lam[0] + mi;
  local[1] = lam[1] + mi;

  return(0);
}

INT GetLocalKoord1(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  INT i,min,mi;
  DOUBLE *p0,*p1,*p2,e0[3],e1[3],e2[3];
  DOUBLE n0[3],n1[3],n2[3],p[3],np[3],hp[3],l;
  DOUBLE a[9];
  DOUBLE aa[4],bb[2],cc[2];
  DOUBLE lam[3],small;
  min = 10000000.0;
  mi = -1;
  small = 0.000001;

  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    p0 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0);
    p1 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1);
    p2 = (DOUBLE*)LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2);

    Minus(e0,global,p0);
    Minus(e1,global,p1);
    Minus(e2,global,p2);

    /* compute 2 vectors in the plane + normalvector */
    Minus(n0,e2,e0);
    l = Lenght(n0);
    Scale(n0,1/l);

    Minus(n1,e2,e1);
    l = Lenght(n1);
    Scale(n1,1/l);

    Cross(n2,n0,n1);
    l = Lenght(n2);
    Scale(n2,1/l);

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
        if(Det2d(aa)!=0)
        {
          cc[0] = p[1] - a[5];
          cc[1] = p[2] - a[8];
          InvMatMult2d(bb,cc,aa);
        }
        else
          printf("%s\n","E R R O R");
      }
    }
    lam[0] = bb[0];
    lam[1] = bb[1];
    lam[2] = 1 - bb[0] - bb[1];

    if( (lam[0]>=-small) && (lam[1]>=-small) && (lam[2]>=-small) )
      break;
  }

  if( (lam[0]<0.0) || (lam[1]<0.0) )
  {
    /*		printf("%s %lf %lf\n","WARNING1", local[0], local[1]);*/
    if( (lam[0]<0.0) && (lam[0]>-0.000001) )
      lam[0] = 0.0;
    if( (lam[1]<0.0) && (lam[1]>-0.000001) )
      lam[1] = 0.0;
  }
  local[0] = lam[0] + i;
  local[1] = lam[1] + i;

  /*	printf("%f %f %d %f %f %f\n", local[0], local[1], i, global[0], global[1], global[2]);*/
  return(0);
}

static INT Compare_Points(LGM_POINT *p, DOUBLE *global)
{
  DOUBLE small;

  small = 0.00001;

  if( sqrt( (LGM_POINT_POS(p)[0]-global[0])*(LGM_POINT_POS(p)[0]-global[0])
            +         (LGM_POINT_POS(p)[1]-global[1])*(LGM_POINT_POS(p)[1]-global[1])
            +         (LGM_POINT_POS(p)[2]-global[2])*(LGM_POINT_POS(p)[2]-global[2]) ) < small )
    return(1);
  else
    return(0);
}

static INT DiscretizeDomain (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h)
{
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  LGM_TRIANGLE *theTriangle;
  LGM_BNDP *theBndP,*testBndP;
  INT ii, i,j,n,m,id, l;
  INT *nRef, *newID;
  INT size,a,b,flag, counter;
  DOUBLE local1,local2;
  DOUBLE global[3],local[2];
  LGM_POINT *PointList, Point;
  INT npoints, npsurface, nPointList;

  /* allocate Memory for the mesh */
  npoints = 0;
  for(theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
    npoints = npoints + LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine));
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    npsurface = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
      npsurface = npsurface - ( LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,i))) - 1);
    npoints = npoints + npsurface;
  }

  PointList = (LGM_POINT *) GetTmpMem(Heap,(npoints+1)*sizeof(LGM_POINT));
  if (PointList==NULL)
    return (1);
  nRef = (INT *) GetTmpMem(Heap,(npoints+1)*sizeof(INT));
  if (nRef==NULL)
    return (1);
  newID = (INT *) GetTmpMem(Heap,(npoints+1)*sizeof(INT));
  if (newID==NULL)
    return (1);
  for (i=0; i<=npoints; i++)
  {
    nRef[i] = 0;
    newID[i] = -1;
  }

  nPointList = 0;
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    for(j=0; j<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); j++)
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
        }
      }
      if(flag==-1)
      {
        PointList[nPointList].position[0] =  global[0];
        PointList[nPointList].position[1] =  global[1];
        PointList[nPointList].position[2] =  global[2];
        nPointList++;
        nRef[i]++;
      }
    }
  }

  theMesh->theBndPs = (BNDP**)GetTmpMem(Heap,sizeof(LGM_BNDP*)*(nPointList+1));
  if (theMesh->theBndPs == NULL)
    return (NULL);

  for(i=0; i<nPointList; i++)
  {
    size = sizeof(LGM_BNDP)+(nRef[i])*sizeof(struct lgm_bndp_surf);
    theMesh->theBndPs[i] = (BNDP*)GetFreelistMemory(Heap,size);
    if(theMesh->theBndPs[i]==NULL)
      return(NULL);
    LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])) = 0;
    a = LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i]));
  }
  theMesh->nBndP = nPointList;

  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    for(j=0; j<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); j++)
    {
      local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,0);
      local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),j,1);
      Surface_Local2Global(theSurface,global,local);

      for(i=0; i<nPointList+1; i++)
      {
        Point = PointList[i];
        if(Compare_Points(&Point, global)==1)
        {
          a = LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i]));
          LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])))[0] = local[0];
          LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i])))[1] = local[1];
          LGM_BNDP_SURFACE(BNDP2LGM(theMesh->theBndPs[i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i]))) = theSurface;
          LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[i]))++;
          i = nPointList+1;
        }

      }
    }
  }

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
  DOUBLE length_of_line,lenght_of_segment,dist,dist1,dist2,dist3,slocal;
  DOUBLE lh, lh1, lh2, rest_length, x, y, z, in[5], dummy;

  if(h>0.0)
  {
    /* uniform dicretisation */
    length_of_line = 0;
    for(i=StartIndex; i<EndIndex-1; i++)
      length_of_line = length_of_line + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));

    /* calculate number of points on the line */
    npoints = floor(length_of_line / h + 0.51) + 1;
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
        if( (slocal<0.0) && (slocal>-0.0001) )
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
        if( (slocal<0.0) && (slocal>-0.0001) )
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
  DOUBLE lh, lh1, lh2, rest_length, x, y, z, in[5], dummy;

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
      if( (slocal<0.0) && (slocal>-0.0001) )
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
        if( (slocal<0.0) && (slocal>-0.0001) )
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
        if( (slocal<0.0) && (slocal>-0.0001) )
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
        if( (slocal<0.0) && (slocal>-0.0001) )
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
        if( (slocal<0.0) && (slocal>-0.0001) )
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
}

static INT DiscretizeLine (HEAP *Heap, LGM_LINE *theLine, DOUBLE h, LGM_POINT *pointlist, INT norp)
{
  INT i, j, k, np, npoints, folds, StartIndex, EndIndex, Disc_nPoints, resolve;
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
        folds = folds++;
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
    printf("%s %d\n", "folds: ", folds);
  /* store npoints in the line-structure */
  LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine)) = npoints;

  /* allocate memery for local coordinates */
  LGM_LINE_LINEDISC(theLine)->local = (DOUBLE*)GetTmpMem(Heap,sizeof(DOUBLE)*(npoints+5));

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
static INT TransferLines2Mesh (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h)
{
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine, *aLine;
  LGM_BNDP *theBndP;
  INT i,j,ni,nl,offset,count,count1,*nRef,size;
  DOUBLE local[2],global[3];

  offset = theMesh->nBndP;

  /* calculate number of lines */
  /* ************************************************** */
  nl = 0;
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
    nl = nl + LGM_SURFACE_NLINE(theSurface);
  /* ************************************************** */

  nRef = (INT *) GetTmpMem(Heap,(nl+1)*sizeof(INT));
  if (nRef==NULL)
    return (1);
  for (i=0; i<nl; i++)
    nRef[i] = 0;

  /* calculate number of Surfaces which contain theLine */
  count = 0;                    /* number of current line */
  count1 = 0;

  aLine = FirstLine(theDomain);

  do
  {
    aLine = FirstLine(theDomain);
    count1 = 0;
    while( (aLine!=NULL) && (count1<count) )
    {
      aLine = NextLine(theDomain);
      count1++;
    }
    for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
    {
      for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
        if(LGM_SURFACE_LINE(theSurface,i)==aLine)
          nRef[count]++;
    }
    count++;
  }
  while(aLine!=NULL);

  /* ni number of new inner points an all lines */
  /* (without begin- and end-points */
  ni = 0;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    ni = ni + LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine)) - 2;
  }

  /* put new points into the mesh-structure */

  /* insert the points off the line */
  for (i=0; i<ni; i++)
  {
    size = sizeof(LGM_BNDP)+(2)*sizeof(struct lgm_bndp_surf);
    theMesh->theBndPs[offset + i] = (BNDP*)GetFreelistMemory(Heap,size);
    LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[offset + i])) = 0;

    /*		theBndP = theBndPList+i;
                    theMesh->theBndPs[offset + i] = (BNDP*)theBndP;
                    theMesh->theBndPs[theMesh->nBndP++] = (BNDP*)theBndP;*/
  }

  count = 0;
  count1 = 0;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    for(i=1; i<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1; i++)
    {
      size = sizeof(LGM_BNDP)+(nRef[count]-1)*sizeof(struct lgm_bndp_surf);
      theMesh->theBndPs[offset + count1] = (BNDP*)GetFreelistMemory(Heap,size);
      count1++;
    }
    count++;
  }

  count = 0;                    /* number of current line */
  count1 = 0;

  aLine = FirstLine(theDomain);

  do
  {
    aLine = FirstLine(theDomain);
    count1 = 0;
    while( (aLine!=NULL) && (count1<count) )
    {
      aLine = NextLine(theDomain);
      count1++;
    }
    if(aLine!=NULL)
    {
      for(i=1; i<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(aLine))-1; i++)
      {
        theBndP = BNDP2LGM(theMesh->theBndPs[theMesh->nBndP]);
        theMesh->nBndP++;
        for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
        {
          for(j=0; j<LGM_SURFACE_NLINE(theSurface); j++)
          {
            if(LGM_SURFACE_LINE(theSurface,j)==aLine)
            {
              /* add this Surface to the List in BndP */
              LGM_BNDP_SURFACE(theBndP,LGM_BNDP_N(theBndP)) = theSurface;

              /* calculate glocal coordinates */
              Line_Local2Global(aLine,global,i);
              /* calculate local coordinates */
              GetLocalKoord(theSurface,global,local);
              LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[0] = local[0];
              LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[1] = local[1];
              LGM_BNDP_N(theBndP)++;
            }
          }

        }
      }
      count++;
    }
  }
  while(aLine!=NULL);

  return (0);
}

static INT TransferSurfaces2Mesh (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h)
{
  INT i,j,offset,id;
  LGM_BNDP *theBndPList;
  DOUBLE global[3],globalbndp[3],local[2];
  INT *newId,newpoints,oldpoints,size;
  DOUBLE small;

  /* find id's of the trianglecorners */
  newId = (INT*)GetTmpMem(Heap,LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))*sizeof(INT));
  for(i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); i++)
    newId[i] = - 1;

  for (i=0; i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); i++)
  {
    local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,0);
    local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,1);
    if(Surface_Local2Global(theSurface,global,local))
      return(1);

    for (j=0; j<theMesh->nBndP; j++)
    {
      if (BNDP_Global(theMesh->theBndPs[j],globalbndp))
        return (1);

      small = 0.00001;
      if( sqrt( (global[0]-globalbndp[0]) * (global[0]-globalbndp[0])
                + (global[1]-globalbndp[1]) * (global[1]-globalbndp[1])
                + (global[2]-globalbndp[2]) * (global[2]-globalbndp[2]) ) < 0.0000001 )
        newId[i] = j;
    }
  }

  for (i=0; i<LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)); i++)
  {
    id = LGM_SURFACE_LEFT(theSurface);
    if (id!=0)
    {
      theMesh->Side_corners[id][theMesh->nSides[id]] = 3;
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] =
        newId[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,0)];
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] =
        newId[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,2)];
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][2] =
        newId[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,1)];
      theMesh->nSides[id]++;
    }
    id = LGM_SURFACE_RIGHT(theSurface);
    if (id!=0)
    {
      theMesh->Side_corners[id][theMesh->nSides[id]] = 3;
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] =
        newId[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,0)];
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] =
        newId[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,1)];
      theMesh->Side_corner_ids[id][theMesh->nSides[id]][2] =
        newId[LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,2)];
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
          LGM_LINE_ID((lineptr)[folds]) = 42;
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

static INT Get_Id(LGM_LINE *theLine, LGM_SURFACE *theSurface, INT index1, INT index2)
{
  INT i, j, k,  flag1, flag2, start_id, end_id;
  DOUBLE global[3], startpoint[3], endpoint[3], local[2];

  flag1 = 0;
  flag2 = 0;

  Line_Local2Global(theLine,startpoint,index1);
  Line_Local2Global(theLine,endpoint,index2);

  start_id = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));
  for(k=0; k<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); k++)
  {
    local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),k,0);
    local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),k,1);
    Surface_Local2Global(theSurface,global,local);
    if(     (  ( startpoint[0]-global[0]<EPS ) && ( global[0]-startpoint[0]<EPS ) )
            &&  (  ( startpoint[1]-global[1]<EPS ) && ( global[1]-startpoint[1]<EPS ) )
            &&  (  ( startpoint[2]-global[2]<EPS ) && ( global[2]-startpoint[2]<EPS ) ) )
    {
      flag1 = 1;
      start_id = k;
    }
  }
  if(flag1==0)
  {
    AddLinePoint(start_id+1, startpoint[0],startpoint[1],startpoint[2]);
    GetLocalKoord(theSurface,startpoint,local);
    LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
    LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
    LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
  }

  end_id = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));
  for(k=0; k<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)); k++)
  {
    local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),k,0);
    local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),k,1);
    Surface_Local2Global(theSurface,global,local);
    if(     (  ( endpoint[0]-global[0]<EPS ) && ( global[0]-endpoint[0]<EPS ) )
            &&  (  ( endpoint[1]-global[1]<EPS ) && ( global[1]-endpoint[1]<EPS ) )
            &&  (  ( endpoint[2]-global[2]<EPS ) && ( global[2]-endpoint[2]<EPS ) ) )
    {
      flag2 = 1;
      end_id = k;
    }
  }
  if(flag2==0)
  {
    AddLinePoint(end_id+1, endpoint[0],endpoint[1],endpoint[2]);
    GetLocalKoord(theSurface,endpoint,local);
    LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
    LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
    LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
  }

  AddLineSegment(start_id+1, end_id+1);

  return(0);
}

static INT TransferLine2Surface(LGM_LINE *theLine, LGM_SURFACE *theSurface, INT direction)
{
  INT i, j, k, flag1, flag2, start_id, end_id;
  DOUBLE global[3], startpoint[3], endpoint[3], local[2];

  LGM_LINE_USED(theLine) = 1;
  if(direction==1)
  {
    for(j=0; j<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1; j++)
      Get_Id(theLine, theSurface, j, j+1);
    /*		if(LGM_SURFACE_NLINE(theSurface)==1)
                            Get_Id(theLine, theSurface, LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1, 0);	*/
  }
  if(direction==-1)
  {
    for(j=LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1; j>0; j--)
      Get_Id(theLine, theSurface, j, j-1);
    /*		if(LGM_SURFACE_NLINE(theSurface)==1)
                            Get_Id(theLine, theSurface, 0, LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1);*/
  }
  if(direction==-2)
  {
    for(j=0; j<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1; j++)
      Get_Id(theLine, theSurface, j, j+1);
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

      if(    (x1<EPS) && (x1>-EPS)
             && (y1<EPS) && (y1>-EPS)
             && (z1<EPS) && (z1>-EPS)
             && (x2<EPS) && (x2>-EPS)
             && (y2<EPS) && (y2>-EPS)
             && (z2<EPS) && (z2>-EPS) )
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

      if(    (x1<EPS) && (x1>-EPS)
             && (y1<EPS) && (y1>-EPS)
             && (z1<EPS) && (z1>-EPS)
             && (x2<EPS) && (x2>-EPS)
             && (y2<EPS) && (y2>-EPS)
             && (z2<EPS) && (z2>-EPS) )
        direction = -1;
    }
  }
  return(direction);
}

static INT DiscretizeSurface (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h, LGM_POINT *pointlist, INT norp)
{
  INT i,n,ni,j,k,offset,nils,id,ls_offset;
  LGM_BNDP *theBndPList, *theBndP;
  DOUBLE global[3],globalbndp[3],local[2], x1, y1, z1, x2, y2, z2;
  INT ii,oldline,newline,direction,a,b;
  INT nnew[5],*pi;
  INT *oldId,nsp,aksp;
  INT *newId,nId[3], neighbor[3],oldnb,newpoints,oldpoints,size, folds;
  LGM_POINT **ptrlst;
  LGM_LINE **lineptr, *theLine;
  DOUBLE startpoint[3], endpoint[3];
  INT start_id, end_id, flag, all_used;

  /* set used-flags to 0 */
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    LGM_LINE_USED(LGM_SURFACE_LINE(theSurface, i)) = 0;

  folds = 0;
  nsp = 0;
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    nsp = nsp + LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,i)))-1;

  LGM_SURFACE_DISC(theSurface)->local = (DOUBLE **) GetTmpMem(Heap,10*nsp*nsp*sizeof(DOUBLE*));
  for(i=0; i<10*nsp*nsp; i++)
    LGM_SURFACE_DISC(theSurface)->local[i] = (DOUBLE *)  GetTmpMem(Heap,3*sizeof(DOUBLE));

  LGM_SURFACE_DISC(theSurface)->triangle = (INT **) GetTmpMem(Heap,10*nsp*nsp*sizeof(INT*));
  for(i=0; i<10*nsp*nsp; i++)
    LGM_SURFACE_DISC(theSurface)->triangle[i] = (INT *)  GetTmpMem(Heap,4*sizeof(INT));

  InitSurface(LOCAL_H);

  LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)) = 0;
  LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)) = 0;

  ptrlst = (LGM_POINT**)GetTmpMem(Heap,LGM_SURFACE_NPOINT(theSurface)*sizeof(LGM_POINT*));

  if(ptrlst==NULL)
    return(1);

  for(j=0; j<LGM_SURFACE_NPOINT(theSurface); j++)
    ptrlst[j] = LGM_SURFACE_POINT(theSurface,j);

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
  if(LGM_DEBUG)
    for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
      printf("%d %d %d %d\n", i,  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),0),
             LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),1),
             LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),2));
  /*	Search_Neighbors(theSurface, Heap);*/
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

  /* search folds on the Surface */
  folds = Count_Folds_Surface(theSurface);
  if(LGM_DEBUG)
    printf("%s %d\n", "FOLDS: ", folds);

  /* allocate lines */
  if(folds>0)
    if ((lineptr=(LGM_LINE**)GetTmpMem(Heap,sizeof(LGM_LINE*)*folds)) == NULL)
      return (NULL);

  for (i=0; i<folds; i++)
  {
    size = sizeof(LGM_LINE) + 2*sizeof(LGM_POINT);
    if ((lineptr[i] = (LGM_LINE*)GetFreelistMemory(Heap,size)) == NULL)
      return (NULL);
  }

  Get_Folds_Surface(theSurface, lineptr);

  for (i=0; i<folds; i++)
  {
    lineptr[i]->ldisc = (LGM_LINEDISC*)GetTmpMem(Heap,sizeof(LGM_LINEDISC)*1);
    if( (lineptr[i]->ldisc==NULL) )
    {
      printf("%s\n", "Not enough memory");
      return(NULL);
    }
    if (DiscretizeLine(Heap,lineptr[i],h, pointlist, norp))
      return(NULL);
    PrintLineInfo(lineptr[i]);
  }

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

    TransferLine2Surface(theLine, theSurface, direction);

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
              TransferLine2Surface(theLine, theSurface, direction);
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
              TransferLine2Surface(theLine, theSurface, direction);
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
              TransferLine2Surface(theLine, theSurface, direction);
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
              TransferLine2Surface(theLine, theSurface, direction);
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
    if(LGM_DEBUG)
      printf("%d %d %d %d %d %d \n",  LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),0)+1,
             LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),1)+1,
             LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(theSurface,i),2)+1,
             neighbor[0], neighbor[1], neighbor[2]);
  }

  for (i=0; i<folds; i++)
  {
    theLine = lineptr[i];
    TransferLine2Surface(theLine, theSurface, -2);
  }

  /* check used-flags */
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    assert(LGM_LINE_USED(LGM_SURFACE_LINE(theSurface, i))==1);

  oldnb = LGM_SURFACE_NPOINT(theSurface);

  GenerateSurfaceGrid (theSurface, h, 5, 1);

  UserWriteF(" Surface %4d: %4d Triangles\n",
             LGM_SURFACE_ID(theSurface),
             LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)));

  /*	for (i=0; i<folds; i++)
                  DisposeMem(Heap, lineptr[i]);
          DisposeMem(Heap, lineptr);*/

  return (0);
}

#else

static INT DiscretizeSurface (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h, LGM_POINT *pointlist, INT norp)
{
  return (1);
}

#endif

INT Surface_Local2Global (LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  INT ilocal,ilocal1;
  DOUBLE slocal[2];
  ilocal = floor(local[0]);
  ilocal1 = floor(local[1]);
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

  return (0);
}

/* domain interface function: for description see domain.h */
INT BNDP_Global (BNDP *aBndP, DOUBLE *global)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE *local;
  theBndP = BNDP2LGM(aBndP);
  theSurface = LGM_BNDP_SURFACE(theBndP,0);

  local = LGM_BNDP_LOCAL(theBndP,0);
  Surface_Local2Global (theSurface,global,local);

  return (0);
}

/* domain interface function: for description see domain.h */
static INT BNDP_Globali (BNDP *aBndP, DOUBLE *global, INT i)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE *local;

  theBndP = BNDP2LGM(aBndP);
  theSurface = LGM_BNDP_SURFACE(theBndP,i);

  local = LGM_BNDP_LOCAL(theBndP,i);
  Surface_Local2Global (theSurface,global,local);

  return (0);
}

/* domain interface function: for description see domain.h */
INT BNDP_BndCond (BNDP *aBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE global[DOM_PARAM_OFFSET],*local;
  INT ilocal=0;

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
  local = LGM_BNDP_LOCAL(theBndP,i);
  Surface_Local2Global (theSurface,global,local);

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
INT BNDP_BndPDesc (BNDP *aBndP, INT *move, INT *part)
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

/* domain interface function: for description see domain.h */
BNDS *BNDP_CreateBndS (HEAP *Heap, BNDP **aBndP, INT n)
{
  INT i,j,k,i0,j0,k0,count, ilocal,ilocal1, found, direction;
  LGM_BNDP *theBndP1, *theBndP2, *theBndP3;
  LGM_SURFACE *theSurface;
  LGM_BNDS *theBndS;
  DOUBLE loc1[2], loc2[2], loc3[2], local[2], slocal[2];
  DOUBLE globalp0[3],globalp1[3],globalp2[3], global[3];
  DOUBLE small, sp, d, min_d, l, l1;
  DOUBLE A[3], B[3], BNDP_NV[3], Surface_NV[3];

  small = 0.000001;
  assert(n==3);
  theBndP1 = BNDP2LGM(aBndP[0]);
  theBndP2 = BNDP2LGM(aBndP[1]);
  theBndP3 = BNDP2LGM(aBndP[2]);

  BNDP_Global(aBndP[0],globalp0);
  BNDP_Global(aBndP[1],globalp1);
  BNDP_Global(aBndP[2],globalp2);

  global[0] = ( globalp0[0] + globalp1[0] +  globalp2[0] ) / 3;
  global[1] = ( globalp0[1] + globalp1[1] +  globalp2[1] ) / 3;
  global[2] = ( globalp0[2] + globalp1[2] +  globalp2[2] ) / 3;

  l1 = sqrt(      (globalp0[0] - globalp1[0])*(globalp0[0] - globalp1[0]) +
                  (globalp0[1] - globalp1[1])*(globalp0[1] - globalp1[1]) +
                  (globalp0[2] - globalp1[2])*(globalp0[2] - globalp1[2]) );
  l = l1;
  l1 = sqrt(      (globalp0[0] - globalp2[0])*(globalp0[0] - globalp2[0]) +
                  (globalp0[1] - globalp2[1])*(globalp0[1] - globalp2[1]) +
                  (globalp0[2] - globalp2[2])*(globalp0[2] - globalp2[2]) );
  if(l<l1)
    l = l1;
  l1 = sqrt(      (globalp1[0] - globalp2[0])*(globalp1[0] - globalp2[0]) +
                  (globalp1[1] - globalp2[1])*(globalp1[1] - globalp2[1]) +
                  (globalp1[2] - globalp2[2])*(globalp1[2] - globalp2[2]) );
  if(l<l1)
    l = l1;

  count = 0;
  min_d = 1000000.0;
  for (i=0; i<LGM_BNDP_N(theBndP1); i++)
    for (j=0; j<LGM_BNDP_N(theBndP2); j++)
      for (k=0; k<LGM_BNDP_N(theBndP3); k++)
      {
        if((LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j))
           && (LGM_BNDP_SURFACE(theBndP2,j)==LGM_BNDP_SURFACE(theBndP3,k)) )
        {
          /* Zusatzabfrage fuer den Fall, dass 3 Punkte auf zwei Surfaces liegen */
          theSurface = LGM_BNDP_SURFACE(theBndP1,i);
          Check_Surface(theSurface,global,local,&d);
          if( (d<min_d) && (2*d/l<0.16) )
          {
            ilocal = floor(local[0]);
            ilocal1 = floor(local[1]);
            if(ilocal>ilocal1)
              ilocal = ilocal1;
            slocal[0] = local[0]-ilocal;
            slocal[1] = local[1]-ilocal;
            /* Punkt liegt im Innern des Dreiecks */
            if( (slocal[0]>-small) && (slocal[0]<1.0+small)
                && (slocal[1]>-small) && (slocal[1]<1.0+small)
                && (1-slocal[0]-slocal[1]>-small) && (1-slocal[0]-slocal[1]<1.0+small)
                )
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

  if (count==0)
    return (NULL);

  theSurface = LGM_BNDP_SURFACE(theBndP1,i0);

  loc1[0] = LGM_BNDP_LOCAL(theBndP1,i0)[0];
  loc1[1] = LGM_BNDP_LOCAL(theBndP1,i0)[1];
  loc2[0] = LGM_BNDP_LOCAL(theBndP2,j0)[0];
  loc2[1] = LGM_BNDP_LOCAL(theBndP2,j0)[1];
  loc3[0] = LGM_BNDP_LOCAL(theBndP3,k0)[0];
  loc3[1] = LGM_BNDP_LOCAL(theBndP3,k0)[1];
  theBndS = (LGM_BNDS *)GetFreelistMemory(Heap,sizeof(LGM_BNDS));
  LGM_BNDS_SURFACE(theBndS) = theSurface;

  LGM_BNDS_LOCAL(theBndS,0,0) = LGM_BNDP_LOCAL(theBndP1,i0)[0];
  LGM_BNDS_LOCAL(theBndS,0,1) = LGM_BNDP_LOCAL(theBndP1,i0)[1];
  LGM_BNDS_LOCAL(theBndS,1,0) = LGM_BNDP_LOCAL(theBndP2,j0)[0];
  LGM_BNDS_LOCAL(theBndS,1,1) = LGM_BNDP_LOCAL(theBndP2,j0)[1];
  LGM_BNDS_LOCAL(theBndS,2,0) = LGM_BNDP_LOCAL(theBndP3,k0)[0];
  LGM_BNDS_LOCAL(theBndS,2,1) = LGM_BNDP_LOCAL(theBndP3,k0)[1];

  /* lege Richtung fuer die BNDS fest */

  /* bestimme den Normalenvektor der Ebene, die durch
     die 3 BndP's festgelegt wird */
  found = 0;
  direction = 0;

  A[0] = globalp2[0] - globalp0[0];
  A[1] = globalp2[1] - globalp0[1];
  A[2] = globalp2[2] - globalp0[2];
  B[0] = globalp2[0] - globalp1[0];
  B[1] = globalp2[1] - globalp1[1];
  B[2] = globalp2[2] - globalp1[2];

  LGM_VECTOR_PRODUCT(A, B, BNDP_NV);

  /* global ist der Schwerpunkt der Dreiecks,
     dass durch die 3 BNDP's aufgespannt wird.
     Projeziere diesen Punkt auf die Surface und
     finde das zugehoerige Inputdreieck */

  GetLocalKoord(theSurface, global, local);
  ilocal = floor(local[0]);
  ilocal1 = floor(local[1]);
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

  LGM_VECTOR_PRODUCT(A, B, Surface_NV);

  LGM_SCALAR_PRODUCT(BNDP_NV, Surface_NV, sp);


  if(sp>0.0)
    LGM_BNDS_N(theBndS) = n;
  else
    LGM_BNDS_N(theBndS) = -n;

  if(LGM_DEBUG)
    printf("%d %f %f %f %f\n", LGM_BNDS_N(theBndS), global[0], global[1], global[2], min_d);
  return((BNDS *)theBndS);
}

/* domain interface function: for description see domain.h */
BNDP *BNDP_CreateBndP (HEAP *Heap, BNDP *aBndP0, BNDP *aBndP1, DOUBLE lcoord)
{
  LGM_BNDP *theBndP1, *theBndP2, *theBndP;
  LGM_SURFACE *theSurface,*s;
  INT i,j,count,size, max, ilocal, ilocal1;
  DOUBLE globalp1[3],globalp2[3],global[3],local[2], slocal[2];
  DOUBLE small;

  small = 0.000001;
  theBndP1 = BNDP2LGM(aBndP0);
  theBndP2 = BNDP2LGM(aBndP1);
  max = LGM_BNDP_N(theBndP1);
  if(max<LGM_BNDP_N(theBndP2))
    max = LGM_BNDP_N(theBndP2);

  BNDP_Global(aBndP0,globalp1);
  BNDP_Global(aBndP1,globalp2);
  global[0] = ( globalp1[0] + globalp2[0] ) / 2;
  global[1] = ( globalp1[1] + globalp2[1] ) / 2;
  global[2] = ( globalp1[2] + globalp2[2] ) / 2;

  if (lcoord>0.0 && lcoord<1.0)
  {
    theBndP1 = BNDP2LGM(aBndP0);
    theBndP2 = BNDP2LGM(aBndP1);
    count = 0;
    for (i=0; i<LGM_BNDP_N(theBndP1); i++)
      for (j=0; j<LGM_BNDP_N(theBndP2); j++)
        if (LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j))
        {
          theSurface = LGM_BNDP_SURFACE(theBndP1,i);
          GetLocalKoord(theSurface,global,local);
          ilocal = floor(local[0]);
          ilocal1 = floor(local[1]);
          if(ilocal>ilocal1)
            ilocal = ilocal1;
          slocal[0] = local[0]-ilocal;
          slocal[1] = local[1]-ilocal;
          if( (slocal[0]>=-small) && (slocal[0]<=1.0+small)
              && (slocal[1]>=-small) && (slocal[1]<=1.0+small)
              && (1-slocal[0]-slocal[1]>=-small) && (1-slocal[0]-slocal[1]<=1.0+small) )
          {
            if(count==0)
            {
              size = sizeof(LGM_BNDP)+max*sizeof(struct lgm_bndp_surf);
              theBndP = (LGM_BNDP *)GetFreelistMemory(Heap,size);
              LGM_BNDP_N(theBndP) = 0;
            }
            count++;
            s = LGM_BNDP_SURFACE(theBndP1,i);
            LGM_BNDP_SURFACE(theBndP,LGM_BNDP_N(theBndP)) = LGM_BNDP_SURFACE(theBndP1,i);
            LGM_BNDP_N(theBndP)++;

            LGM_BNDP_LOCAL(theBndP,count-1)[0] = local[0];
            LGM_BNDP_LOCAL(theBndP,count-1)[1] = local[1];

            if( (local[0]<0.0) || (local[1]<0.0) )
              printf("%s \n","ALARM3");
          }
        }
    if (count==0)
      return (NULL);

    return((BNDP *)theBndP);
  }

  return (NULL);
}

/* domain interface function: for description see domain.h */
INT BNDP_Dispose (HEAP *Heap, BNDP *aBndP)
{
  LGM_BNDP *theBndP;
  INT size;

  if (aBndP == NULL)
    return(0);

  theBndP = BNDP2LGM(aBndP);
  size = sizeof(LGM_BNDP) + (LGM_BNDP_N(theBndP)-1)*sizeof(struct lgm_bndp_surf);
  return (PutFreelistMemory(Heap,theBndP,size));
}

/* domain interface function: for description see domain.h */
INT BNDP_SaveBndP (BNDP *aBndP)
{
  INT i;
  LGM_BNDP *theBndP;
  int n;
  double d[2];

  theBndP = BNDP2LGM(aBndP);
  n = LGM_BNDP_N(theBndP);
  if (Bio_Write_mint(1,&n)) return (1);
  for (i=0; i<LGM_BNDP_N(theBndP); i++)
  {
    n = LGM_SURFACE_ID(LGM_BNDP_SURFACE(theBndP,i));
    if (Bio_Write_mint(1,&n)) return (1);
    d[0] = LGM_BNDP_LOCAL(theBndP,i)[0];
    d[1] = LGM_BNDP_LOCAL(theBndP,i)[1];
    if (Bio_Write_mdouble(2,d)) return (1);
  }
  return(0);
}

/* domain interface function: for description see domain.h */
BNDP *BNDP_LoadBndP (BVP *theBVP, HEAP *Heap)
{
  LGM_DOMAIN *theDomain;
  LGM_SURFACE *theSurface;
  int i,n,id;
  double local[2];
  LGM_BNDP *theBndP;

  theDomain = BVP2LGM(theBVP);

  if (Bio_Read_mint(1,&n)) return (NULL);
  theBndP = (LGM_BNDP *)GetFreelistMemory(Heap,sizeof(LGM_BNDP)+(n-1)*sizeof(LGM_BNDP_PSURFACE));
  LGM_BNDP_N(theBndP) = n;
  for (i=0; i<n; i++)
  {
    if (Bio_Read_mint(1,&id)) return (NULL);
    for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
      if (LGM_SURFACE_ID(theSurface)==id) break;
    if (theSurface==NULL) return (NULL);
    if (Bio_Read_mdouble(2,local)) return (NULL);
    LGM_BNDP_SURFACE(theBndP,i) = theSurface;
    LGM_BNDP_LOCAL(theBndP,i)[0] = local[0];
    LGM_BNDP_LOCAL(theBndP,i)[1] = local[1];
  }

  return((BNDP *)theBndP);
}

/* domain interface function: for description see domain.h */
INT BNDS_Global (BNDS *aBndS, DOUBLE *local, DOUBLE *global)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;

  /* global coordinates */
  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);

  Surface_Local2Global (theSurface,global,local);

  return(0);
}

/* domain interface function: for description see domain.h */
INT BNDS_BndCond (BNDS *aBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;
  DOUBLE global[DOM_PARAM_OFFSET];
  DOUBLE global0[3],global1[3],global2[3],new_global[3];
  DOUBLE bnds_local[2],new_local[2];
  DOUBLE loc0[2],loc1[2],loc2[2],loc[2];

  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);
  if (LGM_SURFACE_BNDCOND(theSurface)==NULL)
    return (2);

  loc0[0] =  LGM_BNDS_LOCAL(theBndS,0,0);
  loc0[1] =  LGM_BNDS_LOCAL(theBndS,0,1);
  loc1[0] =  LGM_BNDS_LOCAL(theBndS,1,0);
  loc1[1] =  LGM_BNDS_LOCAL(theBndS,1,1);
  loc2[0] =  LGM_BNDS_LOCAL(theBndS,2,0);
  loc2[1] =  LGM_BNDS_LOCAL(theBndS,2,1);

  Surface_Local2Global (theSurface,global0,loc0);
  Surface_Local2Global (theSurface,global1,loc1);
  Surface_Local2Global (theSurface,global2,loc2);

  global[0] = (1-local[0]-local[1]) * global0[0] + local[0] * global1[0] + local[1] * global2[0];
  global[1] = (1-local[0]-local[1]) * global0[1] + local[0] * global1[1] + local[1] * global2[1];
  global[2] = (1-local[0]-local[1]) * global0[2] + local[0] * global1[2] + local[1] * global2[2];

  GetLocalKoord(theSurface,global,loc);

  if (BNDS_Global(aBndS,loc,new_global))
    return (1);

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
INT BNDS_BndSDesc (BNDS *aBndS, INT *left, INT *right, INT *part)
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

  /*	if(LGM_DEBUG)
                  printf("%s %d %d\n", "left right", *left, *right);*/
  return(0);
}

/* domain interface function: for description see domain.h */
BNDP *BNDS_CreateBndP (HEAP *Heap, BNDS *aBndS, DOUBLE *local)
{
  LGM_BNDS *theBndS;
  LGM_BNDP *theBndP;
  LGM_SURFACE *theSurface;
  DOUBLE loc0[2],loc1[2],loc2[2],loc[2];
  DOUBLE global0[3],global1[3],global2[3],global[3];

  if (local[0]<=0.0 || local[0]>=1.0)
    return (NULL);
  if (local[1]<=0.0 || local[1]>=1.0)
    return (NULL);

  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);

  loc0[0] =  LGM_BNDS_LOCAL(theBndS,0,0);
  loc0[1] =  LGM_BNDS_LOCAL(theBndS,0,1);
  loc1[0] =  LGM_BNDS_LOCAL(theBndS,1,0);
  loc1[1] =  LGM_BNDS_LOCAL(theBndS,1,1);
  loc2[0] =  LGM_BNDS_LOCAL(theBndS,2,0);
  loc2[1] =  LGM_BNDS_LOCAL(theBndS,2,1);

  Surface_Local2Global (theSurface,global0,loc0);
  Surface_Local2Global (theSurface,global1,loc1);
  Surface_Local2Global (theSurface,global2,loc2);

  global[0] = local[0] * global0[0] + local[1] * global1[0] + (1-local[0]-local[1]) * global2[0];
  global[1] = local[0] * global0[1] + local[1] * global1[1] + (1-local[0]-local[1]) * global2[1];
  global[2] = local[0] * global0[2] + local[1] * global1[2] + (1-local[0]-local[1]) * global2[2];

  GetLocalKoord(theSurface,global,loc);

  theBndP = (LGM_BNDP *)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
  LGM_BNDP_N(theBndP) = 1;
  LGM_BNDP_SURFACE(theBndP,0) = theSurface;
  LGM_BNDP_LOCAL(theBndP,0)[0] = loc[0];
  LGM_BNDP_LOCAL(theBndP,0)[1] = loc[1];

  return((BNDP *)theBndP);

}



/* auxiliary function for getting the maximum surface-ID */
INT GetMaximumSurfaceID (LGM_DOMAIN *theDomain)
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
