// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lgm_load.c													*/
/*																			*/
/* Purpose:   load lgm_domain from file										*/
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: klaus@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   04.09.96 begin												*/
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

#include "devices.h"
#include "lgm_domain.h"
#include "lgm_load.h"
#include "lgm_transfer.h"
#ifdef __THREEDIM__
        #include "ansys2lgm.h"
#endif
#include "misc.h"
#include "general.h"
#ifdef ModelP
#include "debug.h"
#include "parallel.h"
#endif

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

typedef int (*ReadDomainProc)(HEAP *theHeap, char *filename, LGM_DOMAIN_INFO *domain_info, INT MarkKey);
typedef int (*ReadSizesProc)(LGM_SIZES *lgm_sizes);
typedef int (*ReadSubDomainProc)(int i, LGM_SUBDOMAIN_INFO *subdom_info);
typedef int (*ReadLinesProc)(int i, LGM_LINE_INFO *line_info);
typedef int (*ReadPointsProc)(LGM_POINT_INFO *lgm_point_info);
typedef int (*ReadMeshProc)(HEAP *theHeap, LGM_MESH_INFO *lgm_mesh_info, INT MarkKey);

#if (LGM_DIM==3)
typedef int (*ReadSurfaceProc)(int i, LGM_SURFACE_INFO *surface_info);
#endif

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

#if (LGM_DIM==2)
LGM_LINE        **LinePtrArray          = NULL;
#endif
#if (LGM_DIM==3)
LGM_SURFACE     **SurfacePtrArray       = NULL;
#endif

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT LGM_DEBUG = 0;

static ReadDomainProc ReadDomain;
static ReadSizesProc ReadSizes;
static ReadSubDomainProc ReadSubDomain;
static ReadLinesProc ReadLines;
static ReadPointsProc ReadPoints;
static ReadMeshProc ReadMesh;

#if (LGM_DIM==3)
static ReadSurfaceProc ReadSurface;
#endif

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*
   LGM_LoadDomain - reads the domain from file

   SYNOPSIS:
   INT LGM_LoadDomain (char *filename, char *name, HEAP *theHeap, INT DomainVarID, INT MarkKey);

   PARAMETERS:
   .  filename - name of file to read
   .  domainname - name of domain

   DESCRIPTION:
   function read a lgm-domain from the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

#if (LGM_DIM==2)

LGM_DOMAIN *LGM_LoadDomain (char *filename, char *name, HEAP *theHeap, INT DomainVarID, INT MarkKey)
{
  LGM_DOMAIN *theDomain;
  LGM_DOMAIN_INFO theDomInfo;
  LGM_SUBDOMAIN_INFO theSubdomInfo;
  LGM_SUBDOMAIN *theSubdom;
  LGM_SIZES lgm_sizes;
  INT i, j, MaxLinePerSubdom,MaxPointPerLine, *intptr, id;
  LGM_LINE **LinePtrList;
  LGM_POINT_INFO *piptr;
  LGM_LINE_INFO theLineInfo;
  char *p;

  /* set transfer functions */
  p = filename + strlen(filename) - 4;
  if (strcmp(p,".lgm")==0 || strcmp(filename,"geometry")==0)
  {
    ReadDomain              = LGM_ReadDomain;
    ReadSizes               = LGM_ReadSizes;
    ReadSubDomain   = LGM_ReadSubDomain;
    ReadLines               = LGM_ReadLines;
    ReadPoints              = LGM_ReadPoints;
    ReadMesh                = NULL;
  }
  else
  {
    UserWrite("ERROR: fimename must end with .lgm or .hgm\n");
    return (NULL);
  }

  /* read general information */
  if ((*ReadDomain)(theHeap,filename,&theDomInfo,MarkKey))
  {
    UserWrite("ERROR in LGM_LoadDomain: ReadDomain failed\n");
    return (NULL);
  }

  /* allocate and initialize the LGM_DOMAIN */
  if (ChangeEnvDir("/LGM_BVP")==NULL) return (NULL);
  theDomain = (LGM_DOMAIN *) MakeEnvItem(name,DomainVarID,sizeof(LGM_DOMAIN)+theDomInfo.nSubDomain*sizeof(void*));
  if (theDomain==NULL) {UserWriteF("cannot create Domain %s\n",name); return (NULL);}
  if (theDomInfo.Dimension != LGM_DIM) {UserWrite("cannot load domain: wrong dimension\n"); return (NULL);}
  LGM_DOMAIN_CONVEX(theDomain)            = theDomInfo.Convex;
  LGM_DOMAIN_RADIUS(theDomain)            = 1.0;
  LGM_DOMAIN_MIDPOINT(theDomain)[0]       = 0.0;
  LGM_DOMAIN_MIDPOINT(theDomain)[1]       = 0.0;
  LGM_DOMAIN_NPOINT(theDomain)            = theDomInfo.nPoint;
  LGM_DOMAIN_NSUBDOM(theDomain)           = theDomInfo.nSubDomain;
  LGM_DOMAIN_DOMDATA(theDomain)           = NULL;       /* to fill later */
  strcpy(LGM_DOMAIN_PROBLEMNAME(theDomain),theDomInfo.ProblemName);
  LGM_DOMAIN_PROBLEM(theDomain)           = NULL;
  LGM_DOMAIN_S2P_PTR(theDomain)           = NULL;

  /* read sizes */
  if ((lgm_sizes.Subdom_nLine=(int*)GetTmpMem(theHeap,sizeof(int)*(theDomInfo.nSubDomain+1),MarkKey)) == NULL) return (NULL);
  if ((lgm_sizes.Polyline_nPoint=(int*)GetTmpMem(theHeap,sizeof(int)*theDomInfo.nPolyline,MarkKey)) == NULL) return (NULL);
  if ((*ReadSizes)(&lgm_sizes))
  {
    UserWrite("ERROR in LGM_LoadDomain: ReadSizes failed\n");
    return (NULL);
  }


  /* prepare for LGM_DOMAIN_SUBDOM and LGM_LINE_INFO */
  MaxLinePerSubdom = 0;
  for (i=1; i<=theDomInfo.nSubDomain; i++) MaxLinePerSubdom = MAX(MaxLinePerSubdom,lgm_sizes.Subdom_nLine[i]);
  if ((theSubdomInfo.LineNumber=(int*)GetTmpMem(theHeap,sizeof(int)*MaxLinePerSubdom,MarkKey))==NULL) return (NULL);
  MaxPointPerLine = 0;
  for (i=0; i<theDomInfo.nPolyline; i++) MaxPointPerLine = MAX(MaxPointPerLine,lgm_sizes.Polyline_nPoint[i]);
  if ((theLineInfo.point=(int*)GetTmpMem(theHeap,sizeof(int)*MaxPointPerLine,MarkKey))==NULL) return (NULL);

  /* allocate lines */
  if ((LinePtrList=(LGM_LINE**)GetFreelistMemory(theHeap,sizeof(LGM_LINE*)*theDomInfo.nPolyline)) == NULL) return (NULL);
        #ifdef ModelP
  if ((LinePtrArray=(LGM_LINE**)GetFreelistMemory(theHeap,sizeof(LGM_LINE*)*theDomInfo.nPolyline)) == NULL) return (NULL);
        #endif
  for (i=0; i<theDomInfo.nPolyline; i++)
  {
    if ((LinePtrList[i]=(LGM_LINE*)GetFreelistMemory(theHeap,sizeof(LGM_LINE)+(lgm_sizes.Polyline_nPoint[i]-1)*sizeof(LGM_POINT))) == NULL) return (NULL);
    if ((*ReadLines)(i,&theLineInfo))
    {
      UserWrite("ERROR in LGM_LoadDomain: ReadLines failed\n");
      return (NULL);
    }
    LGM_LINE_ID(LinePtrList[i])     = i;
    LGM_LINE_NPOINT(LinePtrList[i]) = lgm_sizes.Polyline_nPoint[i];
    LGM_LINE_LEFT(LinePtrList[i]) = theLineInfo.left;
    LGM_LINE_RIGHT(LinePtrList[i]) = theLineInfo.right;
    LGM_LINE_BNDCOND(LinePtrList[i]) = NULL;
    LGM_LINE_BEGIN(LinePtrList[i]) = theLineInfo.point[0];
    LGM_LINE_END(LinePtrList[i]) = theLineInfo.point[lgm_sizes.Polyline_nPoint[i]-1];
    for (j=0; j<lgm_sizes.Polyline_nPoint[i]; j++)
    {
      intptr = (INT*)LGM_LINE_POINT(LinePtrList[i],j);
      *intptr = theLineInfo.point[j];                                                           /* store id of point */
    }
                #ifdef ModelP
    PRINTDEBUG(dom,3,(PFMT "LGM_LoadDomain(): i=%d lineptr=%x\n",me,i,LinePtrList[i]));
    LinePtrArray[i] = LinePtrList[i];
                #endif
  }

  /* allocate and initialize subdomains */
  LGM_DOMAIN_SUBDOM(theDomain,0) = NULL;
  for (i=1; i<=theDomInfo.nSubDomain; i++)
  {
    if ((*ReadSubDomain)(i,&theSubdomInfo))
    {
      UserWrite("ERROR in LGM_LoadDomain: ReadSubDomain failed\n");
      return (NULL);
    }
    if ((theSubdom=(LGM_SUBDOMAIN*)GetFreelistMemory(theHeap,sizeof(LGM_SUBDOMAIN)+(lgm_sizes.Subdom_nLine[i]-1)*sizeof(void*)))==NULL)
    {
      return (NULL);
    }
    strcpy(LGM_SUBDOMAIN_UNIT(theSubdom),theSubdomInfo.Unit);
    LGM_DOMAIN_SUBDOM(theDomain,i)          = theSubdom;
    LGM_SUBDOMAIN_ID(theSubdom)                     = i;
    LGM_SUBDOMAIN_SDDATA(theSubdom)         = NULL;             /* to fill later */
    LGM_SUBDOMAIN_NLINE(theSubdom)          = lgm_sizes.Subdom_nLine[i];
    for (j=0; j<lgm_sizes.Subdom_nLine[i]; j++)
      LGM_SUBDOMAIN_LINE(theSubdom,j) = LinePtrList[theSubdomInfo.LineNumber[j]];
  }

  /* read points */
  if ((piptr=(LGM_POINT_INFO*)GetTmpMem(theHeap,sizeof(LGM_POINT_INFO)*theDomInfo.nPoint,MarkKey)) == NULL) return (NULL);
  if ((*ReadPoints)(piptr))
  {
    UserWrite("ERROR in LGM_LoadDomain: ReadPoints failed\n");
    return (NULL);
  }
  /* set positions of points on the lines */
  for (i=0; i<theDomInfo.nPolyline; i++)
  {
    for (j=0; j<lgm_sizes.Polyline_nPoint[i]; j++)
    {
      intptr = (INT*)LGM_LINE_POINT(LinePtrList[i],j);
      id = intptr[0];
      LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[0] = piptr[id].position[0];
      LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[1] = piptr[id].position[1];
    }
  }

  return (theDomain);
}

INT LGM_LoadMesh (HEAP *theHeap, MESH *theMesh, INT MarkKey)
{
  /* if impossible to read mesh, return 1 */
  if (ReadMesh==NULL) return (1);

  /* do the right thing */
  return (1);
}


#endif

#if (LGM_DIM==3)
LGM_DOMAIN *LGM_LoadDomain (char *filename, char *name, HEAP *theHeap, INT DomainVarID, INT MarkKey)
{
  LGM_DOMAIN *theDomain;
  LGM_DOMAIN_INFO theDomInfo;
  LGM_SUBDOMAIN *theSubdom;
  LGM_SUBDOMAIN_INFO theSubdomInfo;
  LGM_SURFACE_INFO theSurfaceInfo;
  LGM_SIZES lgm_sizes;
  LGM_POINT_INFO *piptr;
  INT i, j, k, l, MaxSurfacePerSubdom, MaxTrianglePerSurface, MaxPointPerLS, MaxLinePerSurface,MaxPointPerSurface,size;
  INT *intlineptr,*intsurfaceptr, id,*pi;
  LGM_LINE **LinePtrList;
  LGM_SURFACE **SurfacePtrList;
  LGM_LINE_INFO theLineInfo;
  char *p;

  /* set transfer functions */
  p = filename + strlen(filename) - 4;
  if (strcmp(p,".lgm")==0 || strcmp(filename,"geometry")==0)
  {
    ReadDomain              = LGM_ReadDomain;
    ReadSizes               = LGM_ReadSizes;
    ReadSubDomain   = LGM_ReadSubDomain;
    ReadSurface             = LGM_ReadSurface;
    ReadLines               = LGM_ReadLines;
    ReadPoints              = LGM_ReadPoints;
    ReadMesh                = NULL;
  }
  else if (strcmp(p,".ans")==0)
  {
    ReadDomain              = LGM_ANSYS_ReadDomain;
    ReadSizes               = LGM_ANSYS_ReadSizes;
    ReadSubDomain           = LGM_ANSYS_ReadSubDomain;
    ReadSurface             = LGM_ANSYS_ReadSurface;
    ReadLines               = LGM_ANSYS_ReadLines;
    ReadPoints              = LGM_ANSYS_ReadPoints;
    ReadMesh                        = LGM_ANSYS_ReadMesh;
  }
  else
  {
    UserWrite("ERROR: fimename must end with .lgm or .ans\n");
    return (NULL);
  }

  /* read general information */
  if ((*ReadDomain)(theHeap,filename,&theDomInfo,MarkKey))
    return (NULL);

  /* allocate and initialize the LGM_DOMAIN */
  if (ChangeEnvDir("/LGM_BVP")==NULL) return (NULL);
  theDomain = (LGM_DOMAIN *) MakeEnvItem(name,DomainVarID,sizeof(LGM_DOMAIN)+(theDomInfo.nSubDomain-1)*sizeof(void*));
  if (theDomInfo.Dimension != LGM_DIM) {UserWrite("cannot load domain: wrong dimension\n"); return (NULL);}
  LGM_DOMAIN_CONVEX(theDomain)            = theDomInfo.Convex;
  LGM_DOMAIN_RADIUS(theDomain)            = 1.0;
  LGM_DOMAIN_MIDPOINT(theDomain)[0]       = 0.0;
  LGM_DOMAIN_MIDPOINT(theDomain)[1]       = 0.0;
  LGM_DOMAIN_MIDPOINT(theDomain)[2]       = 0.0;
  LGM_DOMAIN_NPOINT(theDomain)            = theDomInfo.nPoint;
  LGM_DOMAIN_NSUBDOM(theDomain)           = theDomInfo.nSubDomain;
  LGM_DOMAIN_DOMDATA(theDomain)           = NULL;       /* to fill later */
  strcpy(LGM_DOMAIN_PROBLEMNAME(theDomain),theDomInfo.ProblemName);
  LGM_DOMAIN_PROBLEM(theDomain)           = NULL;

  /* read sizes */
  if ((lgm_sizes.Subdom_nSurf=(int*)GetTmpMem(theHeap,sizeof(int)*(theDomInfo.nSubDomain+1),MarkKey)) == NULL) return (NULL);
  if ((lgm_sizes.Surf_nPolyline=(int*)GetTmpMem(theHeap,sizeof(int)*(theDomInfo.nSurface+1),MarkKey)) == NULL) return (NULL);
  if ((lgm_sizes.Surf_nTriangle=(int*)GetTmpMem(theHeap,sizeof(int)*(theDomInfo.nSurface+1),MarkKey)) == NULL) return (NULL);
  if ((lgm_sizes.Surf_nPoint=(int*)GetTmpMem(theHeap,sizeof(int)*(theDomInfo.nSurface+1),MarkKey)) == NULL) return (NULL);
  if ((lgm_sizes.Polyline_nPoint=(int*)GetTmpMem(theHeap,sizeof(int)*(theDomInfo.nPolyline+1),MarkKey)) == NULL) return (NULL);

  if ((*ReadSizes)(&lgm_sizes))
    return (NULL);

  /* prepare for LGM_DOMAIN_SUBDOM and LGM_LINE_INFO */

  MaxSurfacePerSubdom = 0;
  for (i=1; i<=theDomInfo.nSubDomain; i++) MaxSurfacePerSubdom = MAX(MaxSurfacePerSubdom,lgm_sizes.Subdom_nSurf[i]);
  if ((theSubdomInfo.SurfaceNumber=(int*)GetTmpMem(theHeap,sizeof(int)*MaxSurfacePerSubdom,MarkKey))==NULL)
    return (NULL);

  MaxPointPerLS = 0;
  for (i=0; i<theDomInfo.nPolyline; i++)
    MaxPointPerLS = MAX(MaxPointPerLS,lgm_sizes.Polyline_nPoint[i]);
  if ((theLineInfo.point=(int*)GetTmpMem(theHeap,sizeof(int)*MaxPointPerLS,MarkKey))==NULL)
    return (NULL);

  MaxPointPerSurface = 0;
  for (i=0; i<theDomInfo.nSurface; i++)
    MaxPointPerSurface = MAX(MaxPointPerSurface,lgm_sizes.Surf_nPoint[i]);
  if ((theSurfaceInfo.point=(int*)GetTmpMem(theHeap,sizeof(int)*MaxPointPerSurface,MarkKey))==NULL)
    return (NULL);

  MaxTrianglePerSurface = 0;
  for (i=0; i<theDomInfo.nSurface; i++) MaxTrianglePerSurface = MAX(MaxTrianglePerSurface,lgm_sizes.Surf_nTriangle[i]);
  if ((theSurfaceInfo.Triangle=(LGM_TRIANGLE_INFO*)GetTmpMem(theHeap,sizeof(LGM_TRIANGLE_INFO)*MaxTrianglePerSurface,MarkKey))==NULL)
    return (NULL);

  MaxLinePerSurface = 0;
  for (i=0; i<theDomInfo.nSurface; i++)
    MaxLinePerSurface = MAX(MaxLinePerSurface,lgm_sizes.Surf_nPolyline[i]);
  if ((theSurfaceInfo.line=(int*)GetTmpMem(theHeap,sizeof(int)*MaxLinePerSurface,MarkKey))==NULL)
    return (NULL);

  theSurfaceInfo.length = theDomInfo.nPoint;
  if ((theSurfaceInfo.point_list=(int**)GetTmpMem(theHeap,sizeof(int*)*theDomInfo.nPoint,MarkKey))==NULL)
    return (NULL);
  for(i=0; i<theDomInfo.nPoint; i++)
    if( (theSurfaceInfo.point_list[i] = (int*)GetTmpMem(theHeap,sizeof(int)*MAXTRIANGLES,MarkKey)) == NULL )
      return(NULL);

  /* allocate lines */
  if ((LinePtrList=(LGM_LINE**)GetFreelistMemory(theHeap,sizeof(LGM_LINE*)*theDomInfo.nPolyline)) == NULL)
    return (NULL);
  for (i=0; i<theDomInfo.nPolyline; i++)
  {
    size = sizeof(LGM_LINE) + (lgm_sizes.Polyline_nPoint[i]-1)*sizeof(LGM_POINT);
    if ((LinePtrList[i] = (LGM_LINE*)GetFreelistMemory(theHeap,size)) == NULL)
      return (NULL);

    if ((*ReadLines)(i,&theLineInfo))
      return (NULL);
    LGM_LINE_NPOINT(LinePtrList[i]) = lgm_sizes.Polyline_nPoint[i];
    for (j=0; j<lgm_sizes.Polyline_nPoint[i]; j++)
    {
      intlineptr = (INT*)LGM_LINE_POINT(LinePtrList[i],j);
      *intlineptr = theLineInfo.point[j];
    }

    pi = (INT*)LGM_LINE_POINT(LinePtrList[i],0);
    id = pi[0];
    LGM_LINE_BEGIN(LinePtrList[i]) = id;

    pi = (INT*)LGM_LINE_POINT(LinePtrList[i],lgm_sizes.Polyline_nPoint[i]-1);
    id = pi[0];
    LGM_LINE_END(LinePtrList[i]) = id;
    LGM_LINE_ID(LinePtrList[i])  = i;
  }

  /* allocate surfaces */
  if ((SurfacePtrList=(LGM_SURFACE**)GetFreelistMemory(theHeap,sizeof(LGM_SURFACE*)*theDomInfo.nSurface)) == NULL)
    return (NULL);
        #ifdef ModelP
  if ((SurfacePtrArray=(LGM_SURFACE**)GetFreelistMemory(theHeap,sizeof(LGM_SURFACE*)*theDomInfo.nSurface)) == NULL)
    return (NULL);
        #endif
  for (i=0; i<theDomInfo.nSurface; i++)
  {
    size = sizeof(LGM_SURFACE)+(lgm_sizes.Surf_nPoint[i]-1)*sizeof(LGM_POINT);
    if ((SurfacePtrList[i]=(LGM_SURFACE*)GetFreelistMemory(theHeap,size)) == NULL)
      return (NULL);

    if ((LGM_SURFACE_FPOINT(SurfacePtrList[i])=
           (LGM_POINT*)GetFreelistMemory(theHeap,sizeof(LGM_POINT)*lgm_sizes.Surf_nPoint[i])) == NULL)
      return (NULL);
    if ((LGM_SURFACE_FTRIANGLE(SurfacePtrList[i])
           =(LGM_TRIANGLE*)GetFreelistMemory(theHeap,sizeof(LGM_TRIANGLE)*lgm_sizes.Surf_nTriangle[i])) == NULL)
      return (NULL);

    theSurfaceInfo.nPoint = lgm_sizes.Surf_nPoint[i];
    theSurfaceInfo.nTriangles = lgm_sizes.Surf_nTriangle[i];
    if ((*ReadSurface)(i,&theSurfaceInfo))
      return (NULL);

    LGM_SURFACE_NPOINT(SurfacePtrList[i])           = lgm_sizes.Surf_nPoint[i];
    LGM_SURFACE_ID(SurfacePtrList[i])                       = i;
    LGM_SURFACE_NTRIANGLE(SurfacePtrList[i])        = theSurfaceInfo.nTriangles;
    LGM_SURFACE_LEFT(SurfacePtrList[i])             = theSurfaceInfo.left;
    LGM_SURFACE_RIGHT(SurfacePtrList[i])            = theSurfaceInfo.right;
    LGM_SURFACE_NLINE(SurfacePtrList[i])            = lgm_sizes.Surf_nPolyline[i];
    LGM_SURFACE_DISC(SurfacePtrList[i])                     = NULL;
    LGM_SURFACE_BNDCOND(SurfacePtrList[i])          = NULL;

    /* set references to lines */
    for (j=0; j<LGM_SURFACE_NLINE(SurfacePtrList[i]); j++)
      LGM_SURFACE_LINE(SurfacePtrList[i],j) = LinePtrList[theSurfaceInfo.line[j]];

    /* set point ids */
    for (j=0; j<LGM_SURFACE_NPOINT(SurfacePtrList[i]); j++)
    {
      intsurfaceptr = (INT*)LGM_SURFACE_POINT(SurfacePtrList[i],j);
      intsurfaceptr[0] = theSurfaceInfo.point[j];
    }

    /* set triangle information */
    for (j=0; j<LGM_SURFACE_NTRIANGLE(SurfacePtrList[i]); j++)
    {
      for (k=0; k<3; k++)
      {
        LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(SurfacePtrList[i],j),k) =
          LGM_SURFACE_POINT(SurfacePtrList[i],theSurfaceInfo.Triangle[j].corner[k]);
        LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(SurfacePtrList[i],j),k) =
          theSurfaceInfo.Triangle[j].corner[k];
        LGM_TRIANGLE_NEIGHBOR(LGM_SURFACE_TRIANGLE(SurfacePtrList[i],j),k) =
          theSurfaceInfo.Triangle[j].neighbor[k];
      }
    }
                #ifdef ModelP
    PRINTDEBUG(dom,3,(PFMT "LGM_LoadDomain(): i=%d surfaceptr=%x\n",me,i,SurfacePtrList[i]));
    SurfacePtrArray[i] = SurfacePtrList[i];
                #endif
  }

  /* allocate and initialize subdomains */
  LGM_DOMAIN_SUBDOM(theDomain,0) = NULL;
  for (i=1; i<=theDomInfo.nSubDomain; i++)
  {
    if ((*ReadSubDomain)(i,&theSubdomInfo)) return (NULL);
    if ((theSubdom=(LGM_SUBDOMAIN*)GetFreelistMemory(theHeap,sizeof(LGM_SUBDOMAIN)+
                                                     (lgm_sizes.Subdom_nSurf[i]-1)*sizeof(void*)))==NULL)
      return (NULL);
    strcpy(LGM_SUBDOMAIN_UNIT(theSubdom),theSubdomInfo.Unit);
    LGM_DOMAIN_SUBDOM(theDomain,i)                  = theSubdom;
    LGM_SUBDOMAIN_ID(theSubdom)                     = i;
    LGM_SUBDOMAIN_SDDATA(theSubdom)                 = NULL;             /* to fill later */
    LGM_SUBDOMAIN_NSURFACE(theSubdom)               = lgm_sizes.Subdom_nSurf[i];
    /*		for (j=0; j<lgm_sizes.Subdom_nSurf[i]; j++)
            LGM_SUBDOMAIN_SURFACE(theSubdom,j)	= (LGM_SURFACE*)theSubdomInfo.SurfaceNumber[j];*/	/* store ids of the surfaces */
    for (j=0; j<lgm_sizes.Subdom_nSurf[i]; j++)
      LGM_SUBDOMAIN_SURFACE(theSubdom,j)      = SurfacePtrList[theSubdomInfo.SurfaceNumber[j]];
  }

  /* read points */
  if ((piptr=(LGM_POINT_INFO*)GetFreelistMemory(theHeap,sizeof(LGM_POINT_INFO)*theDomInfo.nPoint)) == NULL)
    return (NULL);

  if ((*ReadPoints)(piptr))
    return (NULL);

  /* set positions of points on the lines */
  for (i=0; i<theDomInfo.nPolyline; i++)
  {
    for (j=0; j<lgm_sizes.Polyline_nPoint[i]; j++)
    {
      intlineptr = (INT*)LGM_LINE_POINT(LinePtrList[i],j);
      id = intlineptr[0];
      LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[0] = piptr[id].position[0];
      LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[1] = piptr[id].position[1];
      LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[2] = piptr[id].position[2];
      if (LGM_DEBUG) printf("%f %f %f\n",LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[0],
                            LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[1],
                            LGM_POINT_POS(LGM_LINE_POINT(LinePtrList[i],j))[2]);
    }
    if (LGM_DEBUG) printf("\n");
  }

  /* the surfaces */
  for (i=0; i<theDomInfo.nSurface; i++)
  {
    for (j=0; j<lgm_sizes.Surf_nPoint[i]; j++)
    {
      intsurfaceptr = (INT*)LGM_SURFACE_POINT(SurfacePtrList[i],j);
      id = intsurfaceptr[0];
      LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],j))[0] = piptr[id].position[0];
      LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],j))[1] = piptr[id].position[1];
      LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],j))[2] = piptr[id].position[2];
      if (LGM_DEBUG) printf("%f %f %f\n",LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],j))[0],
                            LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],j))[1],
                            LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],j))[2]);

    }
    if (LGM_DEBUG) printf("\n");
  }

  for (i=0; i<theDomInfo.nSurface; i++)
  {
    for (j=0; j<LGM_SURFACE_NTRIANGLE(SurfacePtrList[i]); j++)
    {
      for (k=0; k<3; k++)
      {
        for(l=0; l<LGM_SURFACE_NPOINT(SurfacePtrList[i]); l++)
        {
          /* standard input files global id's for the triangles*/
          if(
            (piptr[LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(SurfacePtrList[i],j),k)].position[0]
             ==LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],l))[0])
            &&
            (piptr[LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(SurfacePtrList[i],j),k)].position[1]
             ==LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],l))[1])
            &&
            (piptr[LGM_TRIANGLE_CORNERID(LGM_SURFACE_TRIANGLE(SurfacePtrList[i],j),k)].position[2]
             ==LGM_POINT_POS(LGM_SURFACE_POINT(SurfacePtrList[i],l))[2])
            )
            LGM_TRIANGLE_CORNER(LGM_SURFACE_TRIANGLE(SurfacePtrList[i],j),k) =
              LGM_SURFACE_POINT(SurfacePtrList[i],l);
        }
      }
    }
  }

  return (theDomain);
}

INT LGM_LoadMesh (HEAP *theHeap, MESH *theMesh, INT MarkKey)
{
  LGM_MESH_INFO lgm_mesh_info;
  INT i;

  /* if impossible to read mesh, return 1 */
  if (ReadMesh==NULL) return (1);

  /* do the right thing */
  if ((*ReadMesh)(theHeap,&lgm_mesh_info,MarkKey)) return (1);

  /* copy mesh_info to mesh and create BNDPs */
  theMesh->nBndP                    = lgm_mesh_info.nBndP;
  theMesh->nInnP                     = lgm_mesh_info.nInnP;
  theMesh->Position                 = lgm_mesh_info.InnPosition;
  theMesh->nSubDomains              = lgm_mesh_info.nSubDomains;
  theMesh->nSides                   = lgm_mesh_info.nSides;
  theMesh->Side_corners             = lgm_mesh_info.Side_corners;
  theMesh->Side_corner_ids          = lgm_mesh_info.Side_corner_ids;
  theMesh->nElements                = lgm_mesh_info.nElements;
  theMesh->Element_corners          = lgm_mesh_info.Element_corners;
  theMesh->Element_corner_ids       = lgm_mesh_info.Element_corner_ids;
  theMesh->nbElements               = lgm_mesh_info.nbElements;
  for (i=0; i<lgm_mesh_info.nBndP; i++)
  {
    theMesh->theBndPs[i] = NULL;
  }
  theMesh->ElemSideOnBnd            = NULL;

  return (1);
}

#endif

/****************************************************************************/
/*
   LGM_LoadDomain - reads the domain from file

   SYNOPSIS:
   INT LGM_LoadDomain (char *filemane, char *domainname);

   PARAMETERS:
   .  filename - name of file to read
   .  domainname - name of domain

   DESCRIPTION:
   function read a lgm-domain from the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

INT InitLGMLoad (void)
{
  if (InitLGMTransfer()) return (1);


  return (0);
}
