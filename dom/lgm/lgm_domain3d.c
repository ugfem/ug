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
#include "domain.h"
#include "general.h"
#include "misc.h"
#include "lgm_domain.h"
#include "lgm_load.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define LGM_VECTOR_PRODUCT(A,B,C)       {(C)[0] = (A)[1]*(B)[2] - (A)[2]*(B)[1];\
                                         (C)[1] = (A)[2]*(B)[0] - (A)[0]*(B)[2];\
                                         (C)[2] = (A)[0]*(B)[1] - (A)[1]*(B)[0];}
#define LGM_SCALAR_PRODUCT(A,B,c)       (c) = ((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]);

#define LGM_BUFFERLEN                                           128

#define BVP2LGM(p)                                                      ((LGM_DOMAIN*)(p))
#define BNDP2LGM(p)                                                     ((LGM_BNDP*)(p))
#define BNDS2LGM(p)                                                     ((LGM_BNDS*)(p))

static char buffer[LGM_BUFFERLEN];

static INT SurfaceInfoId;
static INT LineInfoId;
static INT LGM_DEBUG = 0;

/* data for CVS */
/*static char RCS_ID("$Header$",UG_RCS_STRING);
 */
static INT currSubdom, currSurface, currLine;

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
        LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface,k)) = 0;
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

INT PrintDomainInfo (LGM_DOMAIN *aDomain)
{
  INT i,j,k;

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

INT PrintSubdomainInfo (LGM_DOMAIN *aDomain)
{
  INT i,j,k;
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

INT PrintSurfaceInfo (LGM_SURFACE *aSurface)
{
  INT i,j,k;
  DOUBLE local[2],global[3];

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


  return (0);
}

INT Line_Global2Local (LGM_LINE *theLine, DOUBLE *global, INT i)
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

INT PrintLineInfo (LGM_LINE *aLine)
{
  INT i,j,k;
  DOUBLE global[3];
  printf("********* line-info *********\n");
  printf("%s %d\n","LineId: ",LineInfoId);
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
    Line_Global2Local(aLine,global,i);
    printf("%f %f %f\n",global[0],global[1],global[2]);
  }
  return (0);
}

INT PrintMeshInfo (MESH *mesh)
{
  INT i,j,k;
  DOUBLE global[3];
  INT index,n,type[3];
  DOUBLE value[3];

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

/****************************************************************************/
/*D
   BVP_Save - save a BVP

   SYNOPSIS:
   INT BVP_Save (BVP *theBVP, char *name, char argc, char **argv);

   PARAMETERS:
   .  theBVP - BVP structure
   .  name - name of file
   .  argc, argv - command parameters

   DESCRIPTION:
   This function saves a BVP to file named <name>.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BVP_Save (BVP *theBVP, char *name, char *mgname, HEAP *theHeap, INT argc, char **argv)
{
  UserWrite("SORRY: not implemented yet\n");
  return (1);
}

/****************************************************************************/
/*D
   BVP_Check - check consistency of BVP

   SYNOPSIS:
   INT BVP_Check (BVP *aBVP);

   PARAMETERS:
   .  aBVP - BVP structure
   .  CheckResult - 0 if ok, 1 if error detected

   DESCRIPTION:
   This function checks consistency of BVP

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT BVP_Check (BVP *aBVP)
{
  UserWrite("BVP_Check: not implemented\n");

  return (0);
}

static INT Get_NBNDS_Per_Subdomain      (HEAP *Heap, LGM_DOMAIN *theDomain, INT **sizes, DOUBLE h);
static INT DiscretizeLine                       (HEAP *Heap, LGM_LINE *theLine, MESH *theMesh, DOUBLE h);
static INT Get_NBNDP                            (LGM_DOMAIN *theDomain, INT *nBND, DOUBLE h);
static INT DiscretizeDomain             (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h);
static INT DiscretizeSurface                    (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h);
static INT DiscretizeSurfaceold                 (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h);
static INT DiscretizeLines                              (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h);
static INT TransferLines2Mesh                   (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h);
static INT TransferSurfaces2Mesh                (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h);

MESH *BVP_GenerateMesh (HEAP *Heap, BVP *aBVP, INT argc, char **argv)
{
  LGM_DOMAIN *theDomain;
  LGM_SUBDOMAIN *theSubdom;
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  MESH *mesh;
  INT i,j,n,nBNDP_on_Domain;
  INT *BNDS_Per_Subdom, *p;
  float fValue;
  DOUBLE h;
  INT np,old;
  INT id,*newId;
  INT npoints,npsurface;
  old = 0;
  /* read h-option */
  h = 0.0;
  for (i=1; i<argc; i++)
    if (argv[i][0] == 'h')
    {
      if (sscanf(argv[i],"h %f",&fValue) != 1) return (NULL);
      h = fValue;
    }
  if (h<=0.0) return (NULL);

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



  /* prepare for surface-mesh */
  mesh->nBndP = 0;

  /* discretize lines */
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    theLine->ldisc = (LGM_LINEDISC*)GetTmpMem(Heap,sizeof(LGM_LINEDISC)*1);
    if (DiscretizeLine(Heap,theLine,mesh,h))
      return(NULL);
    if (LGM_DEBUG) PrintLineInfo(theLine);
  }

  /* discretize surfaces */
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    LGM_SURFACE_DISC(theSurface) = (LGM_SURFDISC*)GetTmpMem(Heap,sizeof(LGM_SURFDISC)*1);
    if (DiscretizeSurface(Heap,theSurface,mesh,h))
      return(NULL);
    if (LGM_DEBUG) PrintSurfaceInfo(theSurface);
  }

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

  mesh->theBndPs = (BNDP**)GetTmpMem(Heap,sizeof(LGM_BNDP*)*npoints);
  if (mesh->theBndPs == NULL)
    return (NULL);
  /* discretize domain */
  if (DiscretizeDomain(Heap,theDomain,mesh,h))
    return(NULL);

  if (LGM_DEBUG)
    if (PrintMeshInfo(mesh))
      return (NULL);

  /* discretize lines */
  if (TransferLines2Mesh(Heap,theDomain,mesh,h))
    return(NULL);

  if (LGM_DEBUG)
    if (PrintMeshInfo(mesh))
      return (NULL);

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
    if (LGM_DEBUG)
      if (PrintMeshInfo(mesh))
        return (NULL);
  }

  /* print mesh-info */
  if (LGM_DEBUG)
    if (PrintMeshInfo(mesh))
      return (NULL);

  if (LGM_DEBUG) printf("%s\n","und tschuessss");
  return (mesh);
}

static INT InnerPointsPerSurfaceSegment (LGM_SURFACE *theSurface, DOUBLE h, INT i, INT *n)
{
  DOUBLE d;
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
  INT i,*p,n,pn;
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
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  INT i,pn,n;
  DOUBLE d;

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

static INT DiscretizeDomain (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h)
{
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine;
  LGM_TRIANGLE *theTriangle;
  LGM_BNDP *theBndP,*testBndP;
  INT ii, i,j,n,m,id;
  INT *nRef, *newID;
  INT size,a,b,flag;
  DOUBLE local1,local2;
  DOUBLE global[3],local[2];

  /* get number of corners and number of references to them */
  n=0;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    n = MAX(n,LGM_LINE_BEGIN(theLine));
    n = MAX(n,LGM_LINE_END(theLine));
  }

  nRef = (INT *) GetTmpMem(Heap,(n+1)*sizeof(INT));
  if (nRef==NULL)
    return (1);
  newID = (INT *) GetTmpMem(Heap,(n+1)*sizeof(INT));
  if (newID==NULL)
    return (1);
  for (i=0; i<=n; i++)
  {
    nRef[i] = 0;
    newID[i] = -1;
  }

  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    nRef[LGM_LINE_BEGIN(theLine)]++;
    nRef[LGM_LINE_END(theLine)]++;
  }

  /* create BNDP for each corner which is referenced by three or more lines */
  for (i=0; i<=n; i++)
  {
    if (nRef[i]<3) continue;

    size = sizeof(LGM_BNDP)+(nRef[i])*sizeof(struct lgm_bndp_surf);

    theMesh->theBndPs[theMesh->nBndP] = (BNDP*)GetFreelistMemory(Heap,size);
    LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[theMesh->nBndP])) = 0;
    newID[i] = theMesh->nBndP;
    theMesh->nBndP++;
  }
  /* transform BEGIN- and END-id's in LGM_LINEs to list-ids in the 'theMesh->theBndPs'-list */
  /*	for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
          {
                  assert(newID[LGM_LINE_BEGIN(theLine)]!=-1);
                  LGM_LINE_BEGIN(theLine) = newID[LGM_LINE_BEGIN(theLine)];
                  assert(newID[LGM_LINE_END(theLine)]!=-1);
                  LGM_LINE_END(theLine) = newID[LGM_LINE_END(theLine)];
          }*/

  /* set references on BNDPs */
  for (theSurface=FirstSurface(theDomain); theSurface!=NULL; theSurface=NextSurface(theDomain))
  {
    for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    {
      theLine = LGM_SURFACE_LINE(theSurface,i);
      theBndP = BNDP2LGM(theMesh->theBndPs[newID[LGM_LINE_BEGIN(theLine)]]);

      /* get global coordinates */
      global[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[0];
      global[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[1];
      global[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[2];
      /* calculate local coordinates */
      GetLocalKoord(theSurface,global,local);

      flag = -1;
      for(j=0; j<LGM_BNDP_N(theBndP); j++)
        if(LGM_BNDP_SURFACE(theBndP,j)==theSurface)                             /* Surface schon in der Liste vorhanden */
          flag = 1;
      if(flag!=1)
      {
        LGM_BNDP_SURFACE(theBndP,LGM_BNDP_N(theBndP)) = theSurface;
        LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[0] = local[0];
        LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[1] = local[1];
        LGM_BNDP_N(theBndP)++;
      }

      theBndP = BNDP2LGM(theMesh->theBndPs[newID[LGM_LINE_END(theLine)]]);

      /* get global coordinates */
      global[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,LGM_LINE_NPOINT(theLine)-1))[0];
      global[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,LGM_LINE_NPOINT(theLine)-1))[1];
      global[2] = LGM_POINT_POS(LGM_LINE_POINT(theLine,LGM_LINE_NPOINT(theLine)-1))[2];
      /* calculate local coordinates */
      GetLocalKoord(theSurface,global,local);

      flag = -1;
      for(j=0; j<LGM_BNDP_N(theBndP); j++)
        if(LGM_BNDP_SURFACE(theBndP,j)==theSurface)                             /* Surface schon in der Liste vorhanden */
          flag = 1;
      if(flag!=1)
      {
        LGM_BNDP_SURFACE(theBndP,LGM_BNDP_N(theBndP)) = theSurface;
        LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[0] = local[0];
        LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[1] = local[1];
        LGM_BNDP_N(theBndP)++;
      }

    }
  }

  return (0);
}

#define Lenght(vec)             sqrt(vec[0]*vec[0]      \
                                     +vec[1]*vec[1]      \
                                     +vec[2]*vec[2])

#define Cross(vec,vec1,vec2)    vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]; \
  vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]; \
  vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

#define Minus(sol,vec1,vec2)    sol[0] = vec1[0] - vec2[0];     \
  sol[1] = vec1[1] - vec2[1];     \
  sol[2] = vec1[2] - vec2[2];

#define Plus(sol,vec1,vec2)     sol[0] = vec1[0] + vec2[0];     \
  sol[1] = vec1[1] + vec2[1];     \
  sol[2] = vec1[2] + vec2[2];

#define Mult(vec1,vec2)         ( vec1[0] * vec2[0]     \
                                  + vec1[1] * vec2[1]     \
                                  + vec1[2] * vec2[2])

#define Scale(vec,l)            vec[0] = vec[0] * l;            \
  vec[1] = vec[1] * l;            \
  vec[2] = vec[2] * l;

#define Det(a)                  ( - a[0] * a[4] * a[8]  \
                                  + a[0] * a[5] * a[7]  \
                                  + a[3] * a[1] * a[8]  \
                                  - a[3] * a[2] * a[7]  \
                                  - a[6] * a[1] * a[5]  \
                                  + a[6] * a[2] * a[4])

#define InvMatMult(b,c,a)               b[0] =    ( a[5] * a[7] - a[4] * a[8] ) * c[0]          \
                                               + ( a[1] * a[8] - a[2] * a[7] ) * c[1]          \
                                               + ( a[2] * a[4] - a[1] * a[5] ) * c[2];         \
  b[1] =    ( a[3] * a[8] - a[5] * a[6] ) * c[0]          \
         + ( a[2] * a[6] - a[0] * a[8] ) * c[1]          \
         + ( a[0] * a[5] - a[2] * a[3] ) * c[2];         \
  b[2] =    ( a[4] * a[6] - a[3] * a[7] ) * c[0]          \
         + ( a[0] * a[7] - a[1] * a[6] ) * c[1]          \
         + ( a[1] * a[3] - a[0] * a[4] ) * c[2];         \
  b[0] = b[0] / Det(a);                                   \
  b[1] = b[1] / Det(a);                                   \
  b[2] = b[2] / Det(a);

#define Det2d(a)                ( a[0] * a[3] - a[1] * a[2] )

#define InvMatMult2d(b,c,a)             b[0] =    a[3] * c[0]           \
                                               - a[1] * c[1];          \
  b[1] =  - a[2] * c[0]           \
         + a[0] * c[1];          \
  b[0] = b[0] / Det2d(a);         \
  b[1] = b[1] / Det2d(a);


INT GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  INT i,j,test,min,mi;
  DOUBLE *p0,*p1,*p2,e0[3],e1[3],e2[3],eps;
  DOUBLE n0[3],n1[3],n2[3],p[3],np[3],hp[3],l;
  DOUBLE a[9],b[3],c[3];
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
  local[0] = lam[0] + i;
  local[1] = lam[1] + i;

  if( (local[0]<0.0) || (local[1]<0.0) )
  {
    /*		printf("%s %lf %lf\n","WARNING1", local[0], local[1]);*/
    if( (local[0]<0.0) && (local[0]>-0.000001) )
      local[0] = 0.0;
    if( (local[1]<0.0) && (local[1]>-0.000001) )
      local[1] = 0.0;
  }
  return(0);
}

static INT DiscretizeLine (HEAP *Heap, LGM_LINE *theLine, MESH *theMesh, DOUBLE h)
{
  int i,j,k,npoints;
  DOUBLE length_of_line,lenght_of_segment,dist,dist1,dist2,dist3,slocal;

  length_of_line = 0;
  for(i=0; i<LGM_LINE_NPOINT(theLine)-1; i++)
    length_of_line = length_of_line + LGM_POINT_DIST(LGM_LINE_POINT(theLine,i),LGM_LINE_POINT(theLine,i+1));

  /* calculate number of points on the line */
  npoints = floor(length_of_line / h + 0.51) + 1;
  if(npoints<2)
    npoints = 2;
  lenght_of_segment = length_of_line / (npoints - 1);

  /* store npoints in the line-structure */
  LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine)) = npoints;

  /* allocate memery for local coordinates */
  LGM_LINE_LINEDISC(theLine)->local = (DOUBLE*)GetTmpMem(Heap,sizeof(DOUBLE)*(npoints+1));

  for(i=0; i<npoints; i++)
  {
    dist = i * lenght_of_segment;

    dist1 = 0.0;
    for(j=0; j<LGM_LINE_NPOINT(theLine)-1; j++)
    {
      dist1 = dist1 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
      if(dist1>dist)
        break;
    }
    if(j>LGM_LINE_NPOINT(theLine)-1)
      j = LGM_LINE_NPOINT(theLine) - 1;

    dist2 = 0.0;
    for(k=0; k<j; k++)
      dist2 = dist2 + LGM_POINT_DIST(LGM_LINE_POINT(theLine,k),LGM_LINE_POINT(theLine,k+1));
    dist3 = dist - dist2;
    if(j==LGM_LINE_NPOINT(theLine)-1)
      slocal = 0.0;
    else
      slocal = dist3 / LGM_POINT_DIST(LGM_LINE_POINT(theLine,j),LGM_LINE_POINT(theLine,j+1));
    if( (slocal<0.0) && (slocal>-0.0001) )
      slocal = 0.0;
    LGM_LINEDISC_LOCAL(LGM_LINE_LINEDISC(theLine),i) =
      j + slocal;
  }

  UserWriteF(" Line %4d: %4d Linesegments\n",
             LGM_LINE_ID(theLine),
             LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(theLine))-1);

  return (0);
}

static INT TransferLines2Mesh (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h)
{
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine, *aLine;
  LGM_BNDP *theBndPList, *theBndP;
  INT i,j,ni,nl,offset,count,count1,*nRef,size,help;
  DOUBLE local[2],global[3];
  INT loop,ii;

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
              Line_Global2Local(aLine,global,i);
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
static INT DiscretizeLines (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h)
{
  LGM_SURFACE *theSurface;
  LGM_LINE *theLine, *aLine;
  LGM_BNDP *theBndPList, *theBndP;
  INT i,j,ni,nl,offset,count,count1,*nRef,size,help;
  DOUBLE local[2],global[3];
  INT loop,ii;

  offset = theMesh->nBndP;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    theLine->ldisc = (LGM_LINEDISC*)GetTmpMem(Heap,sizeof(LGM_LINEDISC)*1);
    if (DiscretizeLine(Heap,theLine,theMesh,h))
      return(NULL);
    if (LGM_DEBUG) PrintLineInfo(theLine);
  }

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
  theBndPList = (LGM_BNDP *)GetFreelistMemory(Heap,ni*sizeof(LGM_BNDP));
  for (i=0; i<ni; i++)
  {
    theBndP = theBndPList+i;
    /*		theMesh->theBndPs[offset + i] = (BNDP*)theBndP;*/
    theMesh->theBndPs[theMesh->nBndP++] = (BNDP*)theBndP;
    LGM_BNDP_N(theBndP) = 0;
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
              Line_Global2Local(aLine,global,i);
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
  INT i,n,ni,j,k,offset,nils,id,ls_offset;
  LGM_BNDP *theBndPList, *theBndP;
  DOUBLE global[3],globalbndp[3],local[2];
  INT ii,oldline,newline,direction,a,b;
  INT nnew[5],*pi;
  INT *oldId,nsp,aksp;
  INT *newId,nId[3],oldnb,newpoints,oldpoints,size;
  LGM_POINT **ptrlst;
  DOUBLE small;

  /* insert new points into the mesh */
  offset = theMesh->nBndP;
  newpoints = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));
  oldpoints = newpoints;
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
  {
    newpoints = newpoints - ( LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,i))) - 1);
  }
  oldpoints = oldpoints - newpoints;

  theBndPList = (LGM_BNDP *)GetFreelistMemory(Heap,newpoints*sizeof(LGM_BNDP));

  for (i=0; i<newpoints; i++)
  {
    size = sizeof(LGM_BNDP)+(1*sizeof(struct lgm_bndp_surf));
    theMesh->theBndPs[offset + i] = (BNDP*)GetFreelistMemory(Heap,size);
    theMesh->nBndP++;
    LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[offset + i])) = 0;

    local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i+oldpoints,0);
    local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i+oldpoints,1);

    LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[offset + i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[offset + i])))[0] = local[0];
    LGM_BNDP_LOCAL(BNDP2LGM(theMesh->theBndPs[offset + i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[offset + i])))[1] = local[1];
    LGM_BNDP_SURFACE(BNDP2LGM(theMesh->theBndPs[offset + i]),LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[offset + i]))) = theSurface;
    LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[offset + i]))++;
  }


  offset = theMesh->nBndP;

  /* find id's of the trianglecorners */
  newId = (INT*)GetTmpMem(Heap,LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))*sizeof(INT));

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

      /*			if( (global[0]==globalbndp[0])
                               && (global[1]==globalbndp[1])
                               && (global[2]==globalbndp[2]) )
                                      newId[i] = j;*/

      small = 0.1;
      if( (global[0]-globalbndp[0]<small) && (globalbndp[0]-global[0]<small)
          && (global[1]-globalbndp[1]<small) && (globalbndp[1]-global[1]<small)
          && (global[2]-globalbndp[2]<small) && (globalbndp[2]-global[2]<small) )
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

static INT DiscretizeSurface (HEAP *Heap, LGM_SURFACE *theSurface, MESH *theMesh, DOUBLE h)
{
  INT i,n,ni,j,k,offset,nils,id,ls_offset;
  LGM_BNDP *theBndPList, *theBndP;
  DOUBLE global[3],globalbndp[3],local[2];
  INT ii,oldline,newline,direction,a,b;
  INT nnew[5],*pi;
  INT *oldId,nsp,aksp;
  INT *newId,nId[3],oldnb,newpoints,oldpoints,size;
  LGM_POINT **ptrlst;

  nsp = 0;
  for(i=0; i<LGM_SURFACE_NLINE(theSurface); i++)
    nsp = nsp + LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,i)))-1;

  LGM_SURFACE_DISC(theSurface)->local = (DOUBLE **) GetTmpMem(Heap,nsp*nsp*sizeof(DOUBLE*));
  for(i=0; i<nsp*nsp; i++)
    LGM_SURFACE_DISC(theSurface)->local[i] = (DOUBLE *)  GetTmpMem(Heap,3*sizeof(DOUBLE));

  LGM_SURFACE_DISC(theSurface)->triangle = (INT **) GetTmpMem(Heap,nsp*nsp*sizeof(INT*));
  for(i=0; i<nsp*nsp; i++)
    LGM_SURFACE_DISC(theSurface)->triangle[i] = (INT *)  GetTmpMem(Heap,4*sizeof(INT));

  InitSurface();

  LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)) = 0;
  LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)) = 0;

  ptrlst = (LGM_POINT**)GetTmpMem(Heap,LGM_SURFACE_NPOINT(theSurface)*sizeof(LGM_POINT*));

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
    AddGeomElement(nId[0]+1,nId[1]+1,nId[2]+1);
  }

  /* set line-flag to zero */
  /*	for(i=0; i<LGM_SURFACE_NLINE(theSurface);i++)
          LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface, i)) = 0;*/

  /* input of the linediscretisation */
  ii= 0;
  /* find direction for the first line */
  direction = -1;
  for(i=0; i<LGM_SURFACE_NTRIANGLE(theSurface); i++)
  {
    if( (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0))[0]
         ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[0])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0))[1]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[1])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0))[2]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[2])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1))[0]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[0])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1))[1]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[1])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1))[2]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[2]) )
      direction = 1;
    if( (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1))[0]
         ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[0])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1))[1]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[1])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),1))[2]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[2])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2))[0]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[0])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2))[1]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[1])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2))[2]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[2]) )
      direction = 1;
    if( (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2))[0]
         ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[0])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2))[1]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[1])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),2))[2]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),0))[2])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0))[0]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[0])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0))[1]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[1])
        &&      (LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,i),0))[2]
                 ==LGM_POINT_POS(LGM_LINE_POINT(LGM_SURFACE_LINE(theSurface,0),1))[2]) )
      direction = 1;
  }
  /* first line */
  /*	LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface, 0)) = 1;*/
  if(direction==1)
  {
    for(j=0; j<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,0)))-1; j++)
    {
      Line_Global2Local(LGM_SURFACE_LINE(theSurface,0),global,j);
      AddLinePoint(ii+1,global[0],global[1],global[2]);
      ii++;
      GetLocalKoord(theSurface,global,local);
      LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
      LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
      LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
    }
  }
  else
    for(j=LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,0)))-1; j>0; j--)
    {
      Line_Global2Local(LGM_SURFACE_LINE(theSurface,0),global,j);
      AddLinePoint(ii+1,global[0],global[1],global[2]);
      ii++;
      GetLocalKoord(theSurface,global,local);
      LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
      LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
      LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
    }

  oldline = 0;
  /*	direction = 1;*/
  for(i=1; i<LGM_SURFACE_NLINE(theSurface); i++)
  {
    /* search next line */
    for(k=1; k<LGM_SURFACE_NLINE(theSurface); k++)
    {
      if((k!=oldline) /*&& (LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface, k))==0)*/)
      {
        if(direction==1)
        {
          if (LGM_LINE_END(LGM_SURFACE_LINE(theSurface,oldline))
              ==LGM_LINE_BEGIN(LGM_SURFACE_LINE(theSurface,k)))
          {
            newline = k;
            direction = 1;
            for(j=0; j<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,newline)))-1; j++)
            {
              Line_Global2Local(LGM_SURFACE_LINE(theSurface,newline),global,j);
              AddLinePoint(ii+1,global[0],global[1],global[2]);
              ii++;
              GetLocalKoord(theSurface,global,local);
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
              LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
            }
            oldline = newline;
            /*	LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface, oldline)) = 1;*/
            break;
          }
          if (LGM_LINE_END(LGM_SURFACE_LINE(theSurface,oldline))
              ==LGM_LINE_END(LGM_SURFACE_LINE(theSurface,k)))
          {
            newline = k;
            direction = -1;
            for(j=LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,newline)))-1; j>0; j--)
            {
              Line_Global2Local(LGM_SURFACE_LINE(theSurface,newline),global,j);
              AddLinePoint(ii+1,global[0],global[1],global[2]);
              ii++;
              GetLocalKoord(theSurface,global,local);
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
              LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
            }
            oldline = newline;
            /*	LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface, oldline)) = 1;*/
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
            for(j=0; j<LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,newline)))-1; j++)
            {
              Line_Global2Local(LGM_SURFACE_LINE(theSurface,newline),global,j);
              AddLinePoint(ii+1,global[0],global[1],global[2]);
              ii++;
              GetLocalKoord(theSurface,global,local);
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
              LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
            }
            oldline = newline;
            /*		LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface, oldline)) = 1;*/
            break;
          }
          if (LGM_LINE_BEGIN(LGM_SURFACE_LINE(theSurface,oldline))
              ==LGM_LINE_END(LGM_SURFACE_LINE(theSurface,k)))
          {
            newline = k;
            direction = -1;
            for(j=LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,newline)))-1; j>0; j--)
            {
              Line_Global2Local(LGM_SURFACE_LINE(theSurface,newline),global,j);
              AddLinePoint(ii+1,global[0],global[1],global[2]);
              ii++;
              GetLocalKoord(theSurface,global,local);
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
              LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
              LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
            }
            oldline = newline;
            /*	LGM_LINE_FLAG(LGM_SURFACE_LINE(theSurface, oldline)) = 1;*/
            break;
          }
        }
      }
    }
  }
  for(i=1; i<ii; i++)
  {
    AddLineSegment(i,i+1);
  }
  AddLineSegment(ii,1);

  oldnb = LGM_SURFACE_NPOINT(theSurface);

  GenerateSurfaceGrid (theSurface, h, 1,1);

  /* insert new points into the mesh */
  /*	offset = theMesh->nBndP;
          newpoints = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));
          oldpoints = newpoints;
          for(i=0;i<LGM_SURFACE_NLINE(theSurface);i++)
          {
                  newpoints = newpoints - ( LGM_LINEDISC_NPOINT(LGM_LINE_LINEDISC(LGM_SURFACE_LINE(theSurface,i))) - 1);
          }
          oldpoints = oldpoints - newpoints;

          theBndPList = (LGM_BNDP *)GetFreelistMemory(Heap,newpoints*sizeof(LGM_BNDP));

          for (i=0; i<newpoints; i++)
          {
                  theBndP = theBndPList+i;
                  theMesh->theBndPs[offset + i] = (BNDP*)theBndP;
                  theMesh->nBndP++;
                  LGM_BNDP_N(theBndP) = 0;
                  size = sizeof(LGM_BNDP)+(1*sizeof(struct lgm_bndp_surf));
                  theMesh->theBndPs[offset + i] = (BNDP*)GetFreelistMemory(Heap,size);

                  local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i+oldpoints,0);
                  local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i+oldpoints,1);

                  LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[0] = local[0];
                  LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP))[1] = local[1];
                  LGM_BNDP_SURFACE(theBndP,LGM_BNDP_N(theBndP)) = theSurface;
                  LGM_BNDP_N(theBndP)++;
          }


          offset = theMesh->nBndP;*/

  /* find id's of the trianglecorners */
  /*	newId = (INT*)GetTmpMem(Heap,LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))*sizeof(INT));

          for (i=0;i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));i++)
          {
                  local[0] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,0);
                  local[1] = LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),i,1);
                  if(Surface_Local2Global(theSurface,global,local))
                          return(1);

                  for (j=0;j<theMesh->nBndP;j++)
                  {
                          if (BNDP_Global(theMesh->theBndPs[j],globalbndp))
                                  return (1);

                          if( (global[0]==globalbndp[0])
                           && (global[1]==globalbndp[1])
                           && (global[2]==globalbndp[2]) )
                                  newId[i] = j;

                  }
                  for (j=0;j<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));j++)
                          printf("%d %d\n",j,newId[j]);

          }

          for (i=0;i<LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));i++)
                  printf("%d %d\n",i,newId[i]);*/

  UserWriteF(" Surface %4d: %4d Triangles\n",
             LGM_SURFACE_ID(theSurface),
             LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)));

  return (0);
}

INT Surface_Local2Global (LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local)
{
  INT ilocal,ilocal1,id;
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

INT BNDP_Global (BNDP *aBndP, DOUBLE *global)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE slocal[2],*local;
  INT ilocal,ilocal1,id;
  INT a;
  LGM_POINT *test;
  DOUBLE b;
  theBndP = BNDP2LGM(aBndP);
  theSurface = LGM_BNDP_SURFACE(theBndP,0);

  local = LGM_BNDP_LOCAL(theBndP,0);
  Surface_Local2Global (theSurface,global,local);

  return (0);
}
INT BNDP_Globali (BNDP *aBndP, DOUBLE *global, INT i)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE slocal[2],*local;
  INT ilocal,ilocal1,id;
  INT a;
  LGM_POINT *test;
  DOUBLE b;
  theBndP = BNDP2LGM(aBndP);
  theSurface = LGM_BNDP_SURFACE(theBndP,i);

  local = LGM_BNDP_LOCAL(theBndP,i);
  Surface_Local2Global (theSurface,global,local);

  return (0);
}

INT BNDP_BndCond (BNDP *aBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;
  DOUBLE slocal[2];
  INT ilocal,ilocal1;
  DOUBLE global[DIM],*local;

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
    if ((*LGM_SURFACE_BNDCOND (theSurface))(in,value,type))
      return (1);
  }
  else
  if ((*LGM_SURFACE_BNDCOND (theSurface))(global,value,type))
    return (1);

  return (0);
}

INT BNDP_BndPDesc (BNDP *aBndP, INT *move)
{
  LGM_SURFACE *theSurface;
  LGM_BNDP *theBndP;

  theBndP = BNDP2LGM(aBndP);

  if(LGM_BNDP_N(theBndP)==1)
    *move=1;                                            /* Point in the Interior of the surface */
  else
    *move=0;                                            /* Point on the Boundary of the surface */

  *move=0;
  return(0);
}

BNDS *BNDP_CreateBndS (HEAP *Heap, BNDP **aBndP, INT n)
{
  INT i,j,k, l,i0,j0,k0,count, ilocal,ilocal1, ready, found, direction;
  LGM_BNDP *theBndP1, *theBndP2, *theBndP3;
  LGM_SURFACE *theSurface, *s1,*s2,*s3;
  LGM_BNDS *theBndS;
  DOUBLE loc, loc1[2], loc2[2], loc3[2], local[2], slocal[2];
  DOUBLE globalp0[3],globalp1[3],globalp2[3], global[3];
  DOUBLE small, sp;
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

  count = 0;
  for (i=0; i<LGM_BNDP_N(theBndP1); i++)
    for (j=0; j<LGM_BNDP_N(theBndP2); j++)
      for (k=0; k<LGM_BNDP_N(theBndP3); k++)
      {
        if(
          (LGM_BNDP_SURFACE(theBndP1,i)==LGM_BNDP_SURFACE(theBndP2,j))
          &&
          (LGM_BNDP_SURFACE(theBndP2,j)==LGM_BNDP_SURFACE(theBndP3,k)) )
        {
          /* Zusatzabfrage fuer den Fall, dass 3
             Punkte auf zwei Surfaces liegen */
          theSurface =
            LGM_BNDP_SURFACE(theBndP1,i);
          GetLocalKoord(theSurface,global,local);
          ilocal = floor(local[0]);
          ilocal1 = floor(local[1]);
          if(ilocal>ilocal1)
            ilocal = ilocal1;
          slocal[0] = local[0]-ilocal;
          slocal[1] = local[1]-ilocal;
          /* Punkt liegt im Innern des Dreiecks
           */
          if( (slocal[0]>-small) &&
              (slocal[0]<1.0+small)
              && (slocal[1]>-small) &&
              (slocal[1]<1.0+small)
              && (1-slocal[0]-slocal[1]>-small) &&
              (1-slocal[0]-slocal[1]<1.0+small) )
          {
            theSurface =
              LGM_BNDP_SURFACE(theBndP1,i);
            count++;
            i0=i;
            j0=j;
            k0=k;
          }
        }
      }
  /*	if(count>1)
          printf("%d\n", count);*/
  if (count!=1)
    return (NULL);

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
  LGM_BNDS_LOCAL(theBndS,1,0) = LGM_BNDP_LOCAL(theBndP1,j0)[0];
  LGM_BNDS_LOCAL(theBndS,1,1) = LGM_BNDP_LOCAL(theBndP1,j0)[1];
  LGM_BNDS_LOCAL(theBndS,2,0) = LGM_BNDP_LOCAL(theBndP2,k0)[0];
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

  /* global ist der Scherpunkt der Dreiecks,
     dass durch die 3 BNDP's aufgespannt wird.
     Projeziere diesen Punkt auf die Surface und
     finde das zugehoerige Inputdreieck */

  GetLocalKoord(theSurface, global, local);
  ilocal = floor(local[0]);
  ilocal1 = floor(local[1]);
  if(ilocal>ilocal1)
    ilocal = ilocal1;

  /* bestimme den Normalenvektor der des Surfacedreieck */
  A[0] =
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[0]
    -
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[0];
  A[1] =
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[1]
    -
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[1];
  A[2] =
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[2]
    -
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),0))[2];
  B[0] =
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[0]
    -
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[0];
  B[1] =
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[1]
    -
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[1];
  B[2] =
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),2))[2]
    -
    LGM_POINT_POS(LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(theSurface,ilocal),1))[2];

  LGM_VECTOR_PRODUCT(A, B, Surface_NV);

  LGM_SCALAR_PRODUCT(BNDP_NV, Surface_NV, sp);


  if(sp>1)
    LGM_BNDS_N(theBndS) = n;
  else
    LGM_BNDS_N(theBndS) = -n;

  return((BNDS *)theBndS);
}

BNDP *BNDP_CreateBndP (HEAP *Heap, BNDP *aBndP0, BNDP *aBndP1, DOUBLE lcoord)
{
  LGM_BNDP *theBndP1, *theBndP2, *theBndP;
  LGM_SURFACE *theSurface,*s;
  DOUBLE loc[2], loc1[2], loc2[2];
  INT i,j,count,k,size, max, ilocal, ilocal1;
  DOUBLE globalp1[3],globalp2[3],global[3],local[2], slocal[2];
  DOUBLE small;

  small = 0.0001;
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

INT BNDS_Global (BNDS *aBndS, DOUBLE *local, DOUBLE *global)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;
  DOUBLE loc[2], slocal[2];
  INT ilocal,ilocal1;

  /* global coordinates */
  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);

  Surface_Local2Global (theSurface,global,local);

  return(0);
}

INT BNDS_BndCond (BNDS *aBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;
  DOUBLE global[DIM+1];

  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);
  if (LGM_SURFACE_BNDCOND(theSurface)==NULL)
    return (2);
  if (BNDS_Global(aBndS,local,global))
    return (1);
  if (in!=NULL)
  {
    in[0] = global[0];
    in[1] = global[1];
    in[2] = global[2];
    if ((*LGM_SURFACE_BNDCOND (theSurface))(in,value,type))
      return (1);
  }
  else
  if ((*LGM_SURFACE_BNDCOND (theSurface))(global,value,type))
    return (1);

  return (0);
}

INT BNDS_BndSDesc (BNDS *aBndS, INT *left, INT *right)
{
  LGM_BNDS *theBndS;
  LGM_SURFACE *theSurface;

  theBndS = BNDS2LGM(aBndS);
  theSurface = LGM_BNDS_SURFACE(theBndS);

  if(LGM_BNDS_N(theBndS)>0)
  {
    *left = LGM_SURFACE_LEFT(theSurface);
    *right = LGM_SURFACE_RIGHT(theSurface);
  }
  else
  {
    *right = LGM_SURFACE_LEFT(theSurface);
    *left = LGM_SURFACE_RIGHT(theSurface);
  }
  return(0);
}

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
