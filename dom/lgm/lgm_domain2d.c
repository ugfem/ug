// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lgm_domain2d.c												*/
/*																			*/
/* Purpose:   source for lgm_domain                                                                             */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: klaus@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   08.07.96 begin												*/
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
#include "devices.h"
#include "domain.h"
#include "lgm_load.h"
#include "lgm_domain.h"
#include "misc.h"
#include "general.h"
#include "lgm_transfer.h"

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

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

static char buffer[LGM_BUFFERLEN];

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

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*
   Step-functions for lines in LGM_DOMAIN
 */
/****************************************************************************/

static INT currSubdom, currLine;

static INT ResetLineFlags (LGM_DOMAIN *theDomain)
{
  INT i,j;
  LGM_SUBDOMAIN *theSubdom;

  for (i=1; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    theSubdom = LGM_DOMAIN_SUBDOM(theDomain,i);
    for (j=0; j<LGM_SUBDOMAIN_NLINE(theSubdom); j++)
      LGM_LINE_FLAG(LGM_SUBDOMAIN_LINE(theSubdom,j)) = 0;
  }

  return (0);
}

static LGM_LINE *FirstLine (LGM_DOMAIN *theDomain)
{
  LGM_LINE *theLine;

  if (ResetLineFlags(theDomain)) return (NULL);
  theLine = LGM_SUBDOMAIN_LINE(LGM_DOMAIN_SUBDOM(theDomain,1),0);
  LGM_LINE_FLAG(theLine) = 1;
  currLine = 0; currSubdom = 1;
  return (theLine);
}

static LGM_LINE *helpNextLine (LGM_DOMAIN *theDomain)
{
  LGM_SUBDOMAIN *theSubdom;
  LGM_LINE *theLine;

  theSubdom = LGM_DOMAIN_SUBDOM(theDomain,currSubdom);
  if (currLine<LGM_SUBDOMAIN_NLINE(theSubdom)-1)
  {
    currLine++;
    theLine =  LGM_SUBDOMAIN_LINE(theSubdom,currLine);
    return (theLine);
  }
  else if (currSubdom<LGM_DOMAIN_NSUBDOM(theDomain))
  {
    currSubdom++;
    theSubdom = LGM_DOMAIN_SUBDOM(theDomain,currSubdom);
    currLine = 0;
    theLine =  LGM_SUBDOMAIN_LINE(theSubdom,currLine);
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

/****************************************************************************/
/*D
   BVP_Init - create a lgm-domain and return a mesh

   SYNOPSIS:
   INT BVP_Init (char *name, HEAP *Heap, MESH *Mesh);

   PARAMETERS:
   .  filename - name of file
   .  theHeap - heap

   DESCRIPTION:
   function creates a lgm-domain representaion reading the file <filename>.lgm
   and returns a mesh

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.

   SEE ALSO:
   D*/
/****************************************************************************/

INT SetBoundaryCondition (LGM_DOMAIN *theDomain, BndCondProcPtr BndCond)
{
  INT i,k;
  LGM_SUBDOMAIN *theSubdom;
  LGM_LINE *theLine;

  for (i=1; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    theSubdom = LGM_DOMAIN_SUBDOM(theDomain,i);
    for (k=0; k<LGM_SUBDOMAIN_NLINE(theSubdom); k++)
    {
      theLine = LGM_SUBDOMAIN_LINE(theSubdom,k);
      if (LGM_LINE_LEFT(theLine)*LGM_LINE_RIGHT(theLine)!=0)
        LGM_LINE_BNDCOND(theLine) = NULL;
      else
        LGM_LINE_BNDCOND(theLine) = BndCond;
    }
  }

  return (0);
}

/****************************************************************************/
/*D
        SetDomainSize - set the bounding sphere in the domain

   SYNOPSIS:
   INT SetDomainSize (LGM_DOMAIN *theDomain);

   PARAMETERS:
   .  theDomain - domain

   DESCRIPTION:
   function sets the bounding sphere in the domain

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.

   SEE ALSO:
   D*/
/****************************************************************************/

INT SetDomainSize (LGM_DOMAIN *theDomain)
{
  LGM_PROBLEM *theProblem;
  LGM_LINE *theLine;
  DOUBLE min[2], max[2];
  INT i;

  min[0]=min[1]=MAX_C; max[0]=max[1]=-MAX_C;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
    for (i=0; i<LGM_LINE_NPOINT(theLine); i++)
    {
      min[0] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[0],min[0]);
      min[1] = MIN(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[1],min[1]);
      max[0] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[0],max[0]);
      max[1] = MAX(LGM_POINT_POS(LGM_LINE_POINT(theLine,i))[1],max[1]);
    }
  LGM_DOMAIN_MIDPOINT(theDomain)[0] = 0.5*(min[0]+max[0]);
  LGM_DOMAIN_MIDPOINT(theDomain)[1] = 0.5*(min[1]+max[1]);
  LGM_DOMAIN_RADIUS(theDomain) = 0.55*sqrt((max[0]-min[0])*(max[0]-min[0])+(max[1]-min[1])*(max[1]-min[1]));

  theProblem = LGM_DOMAIN_PROBLEM(theDomain);
  if (LGM_PROBLEM_DOMCONFIG(theProblem)!=NULL)
    if ((*LGM_PROBLEM_DOMCONFIG (theProblem))(min,max))
      return (1);

  return (0);
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/* functions called by script commands										*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

static INT PrintMeshInfo (MESH *mesh)
{
  INT i,j,k;
  DOUBLE global[2];

  printf("********* mesh-info *********\n");
  printf("\n");

  printf("nBNPs = %d\n",(int)mesh->nBndP);
  for (i=0; i<mesh->nBndP; i++)
  {
    if (BNDP_Global(mesh->theBndPs[i],global)) return (1);
    printf("    BNP(%d) = %f %f\n",(int)i,(float)global[0],(float)global[1]);
  }
  printf("\n");

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

static INT Get_NBNDS_Per_Subdomain      (HEAP *Heap, LGM_DOMAIN *theDomain, INT **sizes, DOUBLE h);
static INT DiscretizeLine                       (HEAP *Heap, LGM_LINE *theLine, MESH *theMesh, DOUBLE h);
static INT Get_NBNDP                            (LGM_DOMAIN *theDomain, INT *nBND, DOUBLE h);
static INT DiscretizeDomain             (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h);

/* domain interface function: for description see domain.h */
MESH *BVP_GenerateMesh (HEAP *Heap, BVP *aBVP, INT argc, char **argv)
{
  LGM_DOMAIN *theDomain;
  LGM_LINE *theLine;
  MESH *mesh;
  INT i,j,nBNDP_on_Domain;
  INT *BNDS_Per_Subdom, *p;
  float fValue;
  DOUBLE h;

  /* read h-option */
  h = 0.0;
  for (i=1; i<argc; i++)
    if (argv[i][0] == 'H')
    {
      if (sscanf(argv[i],"H %f",&fValue) != 1) return (NULL);
      h = fValue;
    }
  if (h==0.0)
  {
    for (i=1; i<argc; i++)
      if (argv[i][0] == 'h')
      {
        if (sscanf(argv[i],"h %f",&fValue) != 1) return (NULL);
        h = fValue;
      }
  }
  if (h<=0.0) return (NULL);

  /* set LGM_BVP */
  theDomain = BVP2LGM(aBVP);
  if (theDomain==NULL) return (NULL);

  /* allocate mesh */
  mesh = (MESH *) GetTmpMem(Heap,sizeof(MESH));
  if (mesh == NULL) return(NULL);

  /* init mesh: only surface-mesh */
  mesh->nInnP = 0;        mesh->Position = NULL;
  mesh->nElements = NULL; mesh->Element_corners = NULL;   mesh->Element_corner_ids = NULL;
  mesh->nSubDomains = LGM_DOMAIN_NSUBDOM(theDomain);

  /* get heap for surface-mesh-substructures: the subdomain-dependence */
  mesh->nSides = (INT *) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT));
  if (mesh->nSides==NULL) return(NULL);
  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++) mesh->nSides[i] = 0;
  mesh->Side_corners = (INT **) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT*));
  if (mesh->Side_corners == NULL) return (NULL);
  mesh->Side_corner_ids = (INT ***) GetTmpMem(Heap,(LGM_DOMAIN_NSUBDOM(theDomain)+1)*sizeof(INT**));
  if (mesh->Side_corner_ids == NULL) return (NULL);

  /* allocate further substructers */
  if (Get_NBNDS_Per_Subdomain(Heap,theDomain,&BNDS_Per_Subdom,h)) return (NULL);
  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++)
  {
    mesh->Side_corners[i]    = (INT *)  GetTmpMem(Heap,sizeof(INT)*BNDS_Per_Subdom[i]);
    mesh->Side_corner_ids[i] = (INT **) GetTmpMem(Heap,sizeof(INT*)*BNDS_Per_Subdom[i]);
    p = (INT *) GetTmpMem(Heap,2*sizeof(INT)*BNDS_Per_Subdom[i]);
    for (j=0; j<BNDS_Per_Subdom[i]; j++)
      mesh->Side_corner_ids[i][j] = p+2*j;
  }
  if (Get_NBNDP(theDomain,&nBNDP_on_Domain,h)) return (NULL);
  mesh->theBndPs = (BNDP**)GetTmpMem(Heap,sizeof(LGM_BNDP*)*nBNDP_on_Domain);

  /* prepare for surface-mesh */
  mesh->nBndP = 0;

  /* discretize domain */
  if (DiscretizeDomain(Heap,theDomain,mesh,h)) return(NULL);

  /* discretize lines */
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
    if (DiscretizeLine(Heap,theLine,mesh,h)) return(NULL);

  /* print mesh-info */
  /*if (PrintMeshInfo(mesh)) return (NULL);*/

  mesh->VertexLevel = NULL;
  mesh->VertexPrio = NULL;

  return (mesh);
}

/****************************************************************************/
/*																			*/
/* functions concerned with the quality of the surface-mesh					*/
/*																			*/
/****************************************************************************/

static INT InnerPointsPerLineSegment (LGM_LINE *theLine, DOUBLE h, INT i, INT *n)
{
  DOUBLE d;

  if (i<0 || i>LGM_LINE_NPOINT(theLine)-2 || h<=0.0) return (1);
  d = LGM_LINE_POINTDIST(theLine,i,i+1);
  *n = floor(d/h);

  return (0);
}

static INT Get_NBNDS_Per_Subdomain (HEAP *Heap, LGM_DOMAIN *theDomain, INT **sizes, DOUBLE h)
{
  INT i,*p,n,pn;
  LGM_LINE *theLine;

  /* get heap for information */
  p = (INT *) GetTmpMem(Heap,sizeof(INT)*(LGM_DOMAIN_NSUBDOM(theDomain)+1));
  if (p == NULL) return(1);
  *sizes = p;

  /* get information */
  for (i=0; i<=LGM_DOMAIN_NSUBDOM(theDomain); i++) p[i]=0;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    n = 0;
    for (i=0; i<LGM_LINE_NPOINT(theLine)-1; i++)
    {
      if (InnerPointsPerLineSegment(theLine,h,i,&pn)) return (1);
      n += pn;
    }
    p[LGM_LINE_LEFT(theLine)]  += LGM_LINE_NPOINT(theLine)-1+n;
    p[LGM_LINE_RIGHT(theLine)] += LGM_LINE_NPOINT(theLine)-1+n;
  }

  return (0);
}

static INT Get_NBNDP (LGM_DOMAIN *theDomain, INT *nBND, DOUBLE h)
{
  LGM_LINE *theLine;
  INT i,pn,n;

  n = LGM_DOMAIN_NPOINT(theDomain);
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    for (i=0; i<LGM_LINE_NPOINT(theLine)-1; i++)
    {
      if (InnerPointsPerLineSegment(theLine,h,i,&pn)) return (1);
      n += pn;
    }
  }
  *nBND = n;

  return (0);
}

static INT DiscretizeDomain (HEAP *Heap, LGM_DOMAIN *theDomain, MESH *theMesh, DOUBLE h)
{
  LGM_LINE *theLine;
  LGM_BNDP *theBndP;
  INT i,n;
  INT *nRef, *newID;

  /* get number of corners and number of references to them */
  n=0;
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    n = MAX(n,LGM_LINE_BEGIN(theLine));
    n = MAX(n,LGM_LINE_END(theLine));
  }
  nRef = (INT *) GetTmpMem(Heap,(n+1)*sizeof(INT)); if (nRef==NULL) return (1);
  newID = (INT *) GetTmpMem(Heap,(n+1)*sizeof(INT)); if (newID==NULL) return (1);
  for (i=0; i<=n; i++) {nRef[i] = 0; newID[i] = -1;}
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    nRef[LGM_LINE_BEGIN(theLine)]++;
    nRef[LGM_LINE_END(theLine)]++;
  }

  /* check if referenced correctly */
  for (i=0; i<=n; i++)
    if (nRef[i]==1)
    {
      UserWrite("corners not referenced correctly\n");
      return (1);
    }

  /* create BNDP for each corner which is referenced two times or more */
  for (i=0; i<=n; i++)
  {
    if (nRef[i]<2) continue;
    theMesh->theBndPs[theMesh->nBndP] = (BNDP*)GetFreelistMemory(Heap,sizeof(LGM_BNDP)+(nRef[i]-1)*sizeof(LGM_BNDP_PLINE));
    LGM_BNDP_N(BNDP2LGM(theMesh->theBndPs[theMesh->nBndP])) = 0;
    newID[i] = theMesh->nBndP;
    theMesh->nBndP++;
  }

  /* transform BEGIN- and END-id's in LGM_LINEs to list-ids in the 'theMesh->theBndPs'-list */
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    assert(newID[LGM_LINE_BEGIN(theLine)]!=-1);
    LGM_LINE_BEGIN(theLine) = newID[LGM_LINE_BEGIN(theLine)];
    assert(newID[LGM_LINE_END(theLine)]!=-1);
    LGM_LINE_END(theLine) = newID[LGM_LINE_END(theLine)];
  }

  /* set references on BNDPs */
  for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain))
  {
    theBndP = BNDP2LGM(theMesh->theBndPs[LGM_LINE_BEGIN(theLine)]);
    LGM_BNDP_LINE(theBndP,LGM_BNDP_N(theBndP)) = theLine;
    LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP)) = 0.0;
    LGM_BNDP_N(theBndP)++;

    theBndP = BNDP2LGM(theMesh->theBndPs[LGM_LINE_END(theLine)]);
    LGM_BNDP_LINE(theBndP,LGM_BNDP_N(theBndP)) = theLine;
    LGM_BNDP_LOCAL(theBndP,LGM_BNDP_N(theBndP)) = LGM_LINE_NPOINT(theLine)-1.0;
    LGM_BNDP_N(theBndP)++;
  }

  return (0);
}

static INT DiscretizeLine (HEAP *Heap, LGM_LINE *theLine, MESH *theMesh, DOUBLE h)
{
  INT i,ni,j,offset,nils,id,ls_offset;
  LGM_BNDP *theBndPList, *theBndP;

  offset = theMesh->nBndP;
  ni = LGM_LINE_NPOINT(theLine)-2;
  if (LGM_LINE_NPOINT(theLine)<2) return (1);

  /* insert the points off the line */
  theBndPList = (LGM_BNDP *)GetFreelistMemory(Heap,ni*sizeof(LGM_BNDP));
  for (i=0; i<ni; i++)
  {
    theBndP = theBndPList+i;
    theMesh->theBndPs[theMesh->nBndP++] = (BNDP*)theBndP;
    LGM_BNDP_N(theBndP) = 1;
    LGM_BNDP_LINE(theBndP,0) = theLine;
    LGM_BNDP_LOCAL(theBndP,0) = i+1;
  }

  /* run over all line-segments */
  for (i=0; i<=ni; i++)
  {
    /* get nb of inner points of that line-segment */
    if (InnerPointsPerLineSegment(theLine,h,i,&nils)) return (1);

    /* memory for inner points */
    theBndPList = (LGM_BNDP *)GetFreelistMemory(Heap,nils*sizeof(LGM_BNDP));
    ls_offset = theMesh->nBndP;

    /* insert new BNDP in ptr-list of mesh */
    for (j=0; j<nils; j++)
    {
      theBndP = theBndPList+j;
      theMesh->theBndPs[theMesh->nBndP++] = (BNDP*)theBndP;
      LGM_BNDP_N(theBndP) = 1;
      LGM_BNDP_LINE(theBndP,0) = theLine;
      LGM_BNDP_LOCAL(theBndP,0) = i+(j+1.0)/(nils+1.0);
    }

    /* create references to BNDSs */
    for (j=0; j<=nils; j++)
    {
      id = LGM_LINE_LEFT(theLine);
      if (id!=0)
      {
        theMesh->Side_corners[id][theMesh->nSides[id]] = 2;
        if (j==0) theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = (i>0) ? (offset+i-1) : (LGM_LINE_BEGIN(theLine));
        else theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = ls_offset+j-1;
        if (j==nils) theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = (i<ni) ? (offset+i) : (LGM_LINE_END(theLine));
        else theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = ls_offset+j;
        theMesh->nSides[id]++;
      }
      id = LGM_LINE_RIGHT(theLine);
      if (id!=0)
      {
        theMesh->Side_corners[id][theMesh->nSides[id]] = 2;
        if (j==0) theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = (i>0) ? (offset+i-1) : (LGM_LINE_BEGIN(theLine));
        else theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = ls_offset+j-1;
        if (j==nils) theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = (i<ni) ? (offset+i) : (LGM_LINE_END(theLine));
        else theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = ls_offset+j;
        theMesh->nSides[id]++;
      }
    }
  }

  return (0);
}

static INT DiscretizeLine2 (HEAP *Heap, LGM_LINE *theLine, MESH *theMesh, DOUBLE h)
{
  INT i,n,id,offset;
  LGM_BNDP *theBndPList, *theBndP;

  /* number of BNDP on the interior of theLine */
  n = LGM_LINE_NPOINT(theLine)-2;
  if (n<0) return (0);
  offset = theMesh->nBndP;

  /* create, reference and init BNDPs */
  theBndPList = (LGM_BNDP *)GetFreelistMemory(Heap,n*sizeof(LGM_BNDP));
  for (i=0; i<n; i++)
  {
    theBndP = theBndPList+i;
    theMesh->theBndPs[theMesh->nBndP++] = (BNDP*)theBndP;
    LGM_BNDP_N(theBndP) = 1;
    LGM_BNDP_LINE(theBndP,0) = theLine;
    LGM_BNDP_LOCAL(theBndP,0) = i+1;
  }

  /* create, reference and init BNDSs */
  for (i=0; i<=n; i++)
  {
    id = LGM_LINE_LEFT(theLine);
    if (id!=0)
    {
      theMesh->Side_corners[id][theMesh->nSides[id]] = 2;
      if (i==0) theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = LGM_LINE_BEGIN(theLine);
      else theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = offset+i-1;
      if (i==n) theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = LGM_LINE_END(theLine);
      else theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = offset+i;
      theMesh->nSides[id]++;
    }
    if (id == LGM_LINE_RIGHT(theLine)) continue;
    id = LGM_LINE_RIGHT(theLine);
    if (id!=0)
    {
      theMesh->Side_corners[id][theMesh->nSides[id]] = 2;
      if (i==0) theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = LGM_LINE_BEGIN(theLine);
      else theMesh->Side_corner_ids[id][theMesh->nSides[id]][1] = offset+i-1;
      if (i==n) theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = LGM_LINE_END(theLine);
      else theMesh->Side_corner_ids[id][theMesh->nSides[id]][0] = offset+i;
      theMesh->nSides[id]++;
    }
  }

  return (0);
}

/* domain interface function: for description see domain.h */
INT BVP_Check (BVP *aBVP)
{
  INT i,j,k,ret,at_begin,at_left,flags;
  LGM_DOMAIN *theD;
  LGM_SUBDOMAIN *theSD;
  LGM_LINE *theL,*theL2;

  /* cast */
  theD = BVP2LGM(aBVP);
  if (theD==NULL) return (1);

  /* init */
  UserWrite("BVP_Check: ");
  ret = 0;

  /* check subdomains */
  for (i=1; i<=LGM_DOMAIN_NSUBDOM(theD); i++)
  {
    theSD = LGM_DOMAIN_SUBDOM(theD,i);
    if (theSD==NULL)
    {
      if (!ret) UserWrite("\n");
      UserWriteF("Subdomain %d is not referenced from Domain\n",(int)i);
      ret = 1;
      continue;
    }
    for (j=0; j<LGM_SUBDOMAIN_NLINE(theSD); j++)
    {
      theL = LGM_SUBDOMAIN_LINE(theSD,j);
      if (theL==NULL)
      {
        if (!ret) UserWrite("\n");
        UserWriteF("Line %d is not referenced from Subdomain %d\n",(int)LGM_LINE_ID(theL),(int)i);
        ret = 1;
        continue;
      }
      if (LGM_LINE_LEFT(theL)!=i && LGM_LINE_RIGHT(theL)!=i)
      {
        if (!ret) UserWrite("\n");
        UserWriteF("Line %d does not reference Subdomain %d, [%d (left), %d (right)]\n",(int)LGM_LINE_ID(theL),(int)i,(int)LGM_LINE_LEFT(theL),(int)LGM_LINE_RIGHT(theL));
        ret = 1;
        continue;
      }
      if (LGM_LINE_LEFT(theL)==i && LGM_LINE_RIGHT(theL)==i)
      {
        if (!ret) UserWrite("\n");
        UserWriteF("Line %d references Subdomain %d two times [%d (left), %d (right)]\n",(int)LGM_LINE_ID(theL),(int)i,(int)LGM_LINE_LEFT(theL),(int)LGM_LINE_RIGHT(theL));
        ret = 1;
        continue;
      }
      if (LGM_LINE_BEGIN(theL)==LGM_LINE_END(theL))
      {
        if (!ret) UserWrite("\n");
        UserWriteF("Line %d is cyclic\n",(int)LGM_LINE_ID(theL));
        ret = 1;
        continue;
      }
      if (LGM_LINE_LEFT(theL)==i) at_left = 1;
      else at_left = 0;
      for (k=0; k<LGM_SUBDOMAIN_NLINE(theSD); k++)
      {
        if (k==j) continue;
        theL2 = LGM_SUBDOMAIN_LINE(theSD,k);
        if (LGM_LINE_BEGIN(theL2)==LGM_LINE_END(theL))
        {
          at_begin = 1;
          break;
        }
        else if (LGM_LINE_END(theL2)==LGM_LINE_END(theL))
        {
          at_begin = 0;
          break;
        }
      }
      if (k==LGM_SUBDOMAIN_NLINE(theSD))
      {
        if (!ret) UserWrite("\n");
        UserWriteF("Line %d has no successor\n",(int)LGM_LINE_ID(theL));
        ret = 1;
        continue;
      }
      flags = at_left | (at_begin<<1);
      switch (flags)
      {
      case 0 :
        if (LGM_LINE_LEFT(theL2)!=i)
        {
          if (!ret) UserWrite("\n");
          UserWriteF("Subdomain %d: Line %d: right, Succ(inv) %d: right\n",(int)i,(int)LGM_LINE_ID(theL),(int)LGM_LINE_ID(theL2));
          ret = 1;
          continue;
        }
        break;
      case 1 :
        if (LGM_LINE_RIGHT(theL2)!=i)
        {
          if (!ret) UserWrite("\n");
          UserWriteF("Subdomain %d: Line %d: left, Succ(inv) %d: left\n",(int)i,(int)LGM_LINE_ID(theL),(int)LGM_LINE_ID(theL2));
          ret = 1;
          continue;
        }
        break;
      case 2 :
        if (LGM_LINE_RIGHT(theL2)!=i)
        {
          if (!ret) UserWrite("\n");
          UserWriteF("Subdomain %d: Line %d: right, Succ %d: left\n",(int)i,(int)LGM_LINE_ID(theL),(int)LGM_LINE_ID(theL2));
          ret = 1;
          continue;
        }
        break;
      case 3 :
        if (LGM_LINE_LEFT(theL2)!=i)
        {
          if (!ret) UserWrite("\n");
          UserWriteF("Subdomain %d: Line %d: left, Succ %d: right\n",(int)i,(int)LGM_LINE_ID(theL),(int)LGM_LINE_ID(theL2));
          ret = 1;
          continue;
        }
        break;
      }
    }
  }

  if (ret==0) UserWrite("ok\n");

  return (ret);
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/* functions for BNDP														*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/* domain interface function: for description see domain.h */
INT BNDP_Global (BNDP *aBndP, DOUBLE *global)
{
  LGM_LINE *theLine;
  LGM_BNDP *theBndP;
  DOUBLE slocal;
  INT ilocal;

  theBndP = BNDP2LGM(aBndP);
  theLine = LGM_BNDP_LINE(theBndP,0);
  ilocal = floor(LGM_BNDP_LOCAL(theBndP,0));
  slocal = LGM_BNDP_LOCAL(theBndP,0)-ilocal;
  assert(slocal>=0.0);
  assert(ilocal<LGM_LINE_NPOINT(theLine) && ilocal>=0);

  if (ilocal<LGM_LINE_NPOINT(theLine)-1)
  {
    global[0] = (1.0-slocal)*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[0] + slocal*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal+1))[0];
    global[1] = (1.0-slocal)*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[1] + slocal*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal+1))[1];
  }
  else
  {
    assert(slocal==0.0);
    global[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[0];
    global[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[1];
  }

  return(0);
}

/* domain interface function: for description see domain.h */
INT BNDP_BndCond (BNDP *aBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_LINE *theLine;
  LGM_BNDP *theBndP;
  DOUBLE slocal;
  INT ilocal=0;
  DOUBLE global[DOM_PARAM_OFFSET];
  INT id1, id2;

  /* general */
  theBndP = BNDP2LGM(aBndP);
  *n = LGM_BNDP_N(theBndP);
  assert(i>=0 && i<LGM_BNDP_N(theBndP));
  theLine = LGM_BNDP_LINE(theBndP,i);
  if (LGM_LINE_BNDCOND(theLine)==NULL) return (2);

  /* global coordinates */
  ilocal = floor(LGM_BNDP_LOCAL(theBndP,i));
  slocal = LGM_BNDP_LOCAL(theBndP,i)-ilocal;
  assert(slocal>=0.0);
  assert(ilocal<LGM_LINE_NPOINT(theLine) && ilocal>=0);
  if (ilocal<LGM_LINE_NPOINT(theLine)-1)
  {
    global[0] = (1.0-slocal)*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[0] + slocal*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal+1))[0];
    global[1] = (1.0-slocal)*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[1] + slocal*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal+1))[1];
  }
  else
  {
    assert(slocal==0.0);
    global[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[0];
    global[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[1];
  }

  id1 = LGM_LINE_LEFT(theLine);
  id2 = LGM_LINE_RIGHT(theLine);

  /* get values */
  if (in!=NULL)
  {
    in[0] = global[0];
    in[1] = global[1];
    in[DIM] = MAX(id1,id2);
    if ((*LGM_LINE_BNDCOND (theLine))(in,value,type)) return (1);
  }
  else
  {
    global[DIM] = MAX(id1,id2);
    if ((*LGM_LINE_BNDCOND (theLine))(global,value,type)) return (1);
  }

  return (0);
}

/* domain interface function: for description see domain.h */
INT BNDP_BndPDesc (BNDP *aBndP, INT *move, INT *part)
{
  LGM_LINE *theLine;
  LGM_BNDP *theBndP;
  DOUBLE slocal;
  INT ilocal;

  part[0] = 0;
  theBndP = BNDP2LGM(aBndP);
  theLine = LGM_BNDP_LINE(theBndP,0);
  ilocal = floor(LGM_BNDP_LOCAL(theBndP,0));
  slocal = LGM_BNDP_LOCAL(theBndP,0)-ilocal;
  assert(slocal>=0.0);
  assert(ilocal<LGM_LINE_NPOINT(theLine) && ilocal>=0);

  if (ilocal==LGM_LINE_NPOINT(theLine)-1 || ilocal==0) *move=0;
  else *move=1;

  /* HRR_TODO: assign part */
  *part = 0;

  return(0);
}

/* domain interface function: for description see domain.h */
BNDS *BNDP_CreateBndS (HEAP *Heap, BNDP **aBndP, INT n)
{
  INT i,j,i0,j0,count,k;
  LGM_BNDP *theBndP1, *theBndP2;
  LGM_LINE *theLine;
  LGM_BNDS *theBndS;
  DOUBLE loc1, loc2;

  assert(n==2);
  theBndP1 = BNDP2LGM(aBndP[0]); theBndP2 = BNDP2LGM(aBndP[1]);
  count = 0;
  for (i=0; i<LGM_BNDP_N(theBndP1); i++)
    for (j=0; j<LGM_BNDP_N(theBndP2); j++)
      if (LGM_BNDP_LINE(theBndP1,i)==LGM_BNDP_LINE(theBndP2,j))
      {
        theLine = LGM_BNDP_LINE(theBndP1,i);
        count++;
        i0=i; j0=j;
      }
  if (count!=1) return (NULL);
  loc1 = LGM_BNDP_LOCAL(theBndP1,i0);
  loc2 = LGM_BNDP_LOCAL(theBndP2,j0);
  k = (loc1<loc2) ? floor(loc1) : floor(loc2);
  if (loc1-k>1.0 || loc2-k>1.0) return (NULL);
  theBndS = (LGM_BNDS *)GetFreelistMemory(Heap,sizeof(LGM_BNDS));
  LGM_BNDS_LINE(theBndS) = theLine;
  LGM_BNDS_LOCAL(theBndS,0) = LGM_BNDP_LOCAL(theBndP1,i0);
  LGM_BNDS_LOCAL(theBndS,1) = LGM_BNDP_LOCAL(theBndP2,j0);

  return((BNDS *)theBndS);
}

/* domain interface function: for description see domain.h */
BNDP *BNDP_CreateBndP (HEAP *Heap, BNDP *aBndP0, BNDP *aBndP1, DOUBLE lcoord)
{
  LGM_BNDP *theBndP1, *theBndP2, *theBndP;
  LGM_LINE *theLine;
  DOUBLE loc, loc1, loc2;
  INT i,j,count,k;

  if (lcoord>0.0 && lcoord<1.0)
  {
    theBndP1 = BNDP2LGM(aBndP0);
    theBndP2 = BNDP2LGM(aBndP1);
    count = 0;
    for (i=0; i<LGM_BNDP_N(theBndP1); i++)
      for (j=0; j<LGM_BNDP_N(theBndP2); j++)
        if (LGM_BNDP_LINE(theBndP1,i)==LGM_BNDP_LINE(theBndP2,j))
        {
          theLine = LGM_BNDP_LINE(theBndP1,i);
          count++;
          loc1 = LGM_BNDP_LOCAL(theBndP1,i);
          loc2 = LGM_BNDP_LOCAL(theBndP2,j);
          k = (loc1<loc2) ? floor(loc1) : floor(loc2);
          if (loc1-k>1.0 || loc2-k>1.0) return (NULL);
          loc = 0.5*(loc1+loc2);
        }
    if (count!=1) return (NULL);
    theBndP = (LGM_BNDP *)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
    LGM_BNDP_N(theBndP) = 1;
    LGM_BNDP_LINE(theBndP,0) = theLine;
    LGM_BNDP_LOCAL(theBndP,0) = loc;

    return((BNDP *)theBndP);
  }

  return (NULL);
}

/* domain interface function: for description see domain.h */
INT BNDP_Dispose (HEAP *Heap, BNDP *aBndP)
{
  LGM_BNDP *theBndP;
  INT size;

  if (aBndP == NULL) return(0);

  theBndP = BNDP2LGM(aBndP);
  size = sizeof(LGM_BNDP) + (LGM_BNDP_N(theBndP)-1)*sizeof(LGM_BNDP_PLINE);
  return (PutFreelistMemory(Heap,theBndP,size));
}

/* domain interface function: for description see domain.h */
INT BNDP_SaveBndP (BNDP *aBndP)
{
  INT i;
  LGM_BNDP *theBndP;
  int n;
  double d;

  theBndP = BNDP2LGM(aBndP);
  n = LGM_BNDP_N(theBndP);
  if (Bio_Write_mint(1,&n)) return (1);
  for (i=0; i<LGM_BNDP_N(theBndP); i++)
  {
    n = LGM_LINE_ID(LGM_BNDP_LINE(theBndP,i));
    if (Bio_Write_mint(1,&n)) return (1);
    d = LGM_BNDP_LOCAL(theBndP,i);
    if (Bio_Write_mdouble(1,&d)) return (1);
  }

  return(0);
}

/* domain interface function: for description see domain.h */
BNDP *BNDP_LoadBndP (BVP *theBVP, HEAP *Heap)
{
  LGM_DOMAIN *theDomain;
  LGM_LINE *theLine;
  int i,n,id;
  double local;
  LGM_BNDP *theBndP;

  theDomain = BVP2LGM(theBVP);

  if (Bio_Read_mint(1,&n)) return (NULL);
  theBndP = (LGM_BNDP *)GetFreelistMemory(Heap,sizeof(LGM_BNDP)+(n-1)*sizeof(LGM_BNDP_PLINE));
  LGM_BNDP_N(theBndP) = n;
  for (i=0; i<n; i++)
  {
    if (Bio_Read_mint(1,&id)) return (NULL);
    for (theLine=FirstLine(theDomain); theLine!=NULL; theLine=NextLine(theDomain)) if (LGM_LINE_ID(theLine)==id) break;
    if (theLine==NULL) return (NULL);
    if (Bio_Read_mdouble(1,&local)) return (NULL);
    LGM_BNDP_LINE(theBndP,i) = theLine;
    LGM_BNDP_LOCAL(theBndP,i) = local;
  }

  return((BNDP *)theBndP);
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/* functions for BNDS														*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/* domain interface function: for description see domain.h */
INT BNDS_Global (BNDS *aBndS, DOUBLE *local, DOUBLE *global)
{
  LGM_BNDS *theBndS;
  LGM_LINE *theLine;
  DOUBLE loc, slocal;
  INT ilocal;

  theBndS = BNDS2LGM(aBndS);
  theLine = LGM_BNDS_LINE(theBndS);
  loc = (1.0-local[0])*LGM_BNDS_LOCAL(theBndS,0)+local[0]*LGM_BNDS_LOCAL(theBndS,1);

  ilocal = floor(loc);
  slocal = loc-ilocal;
  assert(slocal>=0.0);
  assert(ilocal<LGM_LINE_NPOINT(theLine) && ilocal>=0);

  if (ilocal<LGM_LINE_NPOINT(theLine)-1)
  {
    global[0] = (1.0-slocal)*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[0] + slocal*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal+1))[0];
    global[1] = (1.0-slocal)*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[1] + slocal*LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal+1))[1];
  }
  else
  {
    assert(slocal==0.0);
    global[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[0];
    global[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,ilocal))[1];
  }

  return(0);
}

/* domain interface function: for description see domain.h */
INT BNDS_BndCond (BNDS *aBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value, INT *type)
{
  LGM_BNDS *theBndS;
  LGM_LINE *theLine;
  DOUBLE global[DOM_PARAM_OFFSET+1];
  INT id1, id2;

  theBndS = BNDS2LGM(aBndS);
  theLine = LGM_BNDS_LINE(theBndS);
  if (LGM_LINE_BNDCOND(theLine)==NULL) return (2);
  if (BNDS_Global(aBndS,local,global)) return (1);

  id1 = LGM_LINE_LEFT(theLine);
  id2 = LGM_LINE_RIGHT(theLine);

  if (in!=NULL)
  {
    in[0] = global[0];
    in[1] = global[1];
    in[DIM] = MAX(id1,id2);
    if ((*LGM_LINE_BNDCOND (theLine))(in,value,type)) return (1);
  }
  else
  {
    global[DIM] = MAX(id1,id2);
    if ((*LGM_LINE_BNDCOND (theLine))(global,value,type)) return (1);
  }

  return (0);
}

static INT SideIsCooriented (LGM_BNDS *theBndS)
{
  if (LGM_BNDS_LOCAL(theBndS,1)>LGM_BNDS_LOCAL(theBndS,0))
    return (YES);
  else
    return (NO);
}

/* domain interface function: for description see domain.h */
INT BNDS_BndSDesc (BNDS *aBndS, INT *id, INT *nbid, INT *part)
{
  LGM_BNDS *theBndS;
  LGM_LINE *theLine;
  INT left,right;

  theBndS = BNDS2LGM(aBndS);
  theLine = LGM_BNDS_LINE(theBndS);

  left = LGM_LINE_LEFT(theLine);
  right = LGM_LINE_RIGHT(theLine);

  /* check orientation */
  if (SideIsCooriented(theBndS))
  {
    /* patch and side are co-oriented */
    *id = left;
    *nbid = right;
  }
  else
  {
    /* patch and side are anti-oriented */
    *id = right;
    *nbid = left;
  }

  /* HRR_TODO: assign part */
  *part = 0;

  return(0);
}

/* domain interface function: for description see domain.h */
BNDP *BNDS_CreateBndP (HEAP *Heap, BNDS *aBndS, DOUBLE *local)
{
  LGM_BNDS *theBndS;
  LGM_BNDP *theBndP;
  LGM_LINE *theLine;
  DOUBLE loc;

  if (local[0]<=0.0 || local[0]>=1.0) return (NULL);
  theBndS = BNDS2LGM(aBndS);
  theLine = LGM_BNDS_LINE(theBndS);
  loc = (1.0-local[0])*LGM_BNDS_LOCAL(theBndS,0)+local[0]*LGM_BNDS_LOCAL(theBndS,1);
  theBndP = (LGM_BNDP *)GetFreelistMemory(Heap,sizeof(LGM_BNDP));
  LGM_BNDP_N(theBndP) = 1;
  LGM_BNDP_LINE(theBndP,0) = theLine;
  LGM_BNDP_LOCAL(theBndP,0) = loc;

  return((BNDP *)theBndP);
}

/* domain interface function: for description see domain.h */
INT BVP_Save (BVP *theBVP, char *name, char *mgname, HEAP *theHeap, INT argc, char **argv)
{
  LGM_DOMAIN *lgm_domain;
  LGM_LINE *theLine, **LinePtrList;
  FILE *stream;
  int i,j,npoints,act,nlines;
  LGM_POINT *PointList;
  INT *FlagList;

  /* init */
  lgm_domain = (LGM_DOMAIN*)theBVP;
  if (lgm_domain==NULL) return (1);
  if (lgm_domain->theProblem==NULL) return (1);

  /* open file */
  stream = LGM_WriteOpenFile(name);
  if (stream==NULL) return (1);

  /* save domain-info */
  if (fprintf(stream,"# Domain-Info\n")<0) return (1);
  if (fprintf(stream,"name = %s\n",mgname)<0) return (1);
  if (fprintf(stream,"problemname = %s\n",ENVITEM_NAME(lgm_domain->theProblem))<0) return (1);
  if (fprintf(stream,"convex = %d\n\n",lgm_domain->convex)<0) return (1);

  /* create point lists */
  MarkTmpMem(theHeap);
  npoints = LGM_DOMAIN_NPOINT(lgm_domain);
  PointList = (LGM_POINT *) GetTmpMem(theHeap,npoints*sizeof(LGM_POINT));
  if (PointList==NULL)
  {
    UserWrite("ERROR: cannot allocate memory for PointList\n");
    return (1);
  }
  FlagList = (INT *) GetTmpMem(theHeap,npoints*sizeof(INT));
  if (FlagList==NULL)
  {
    UserWrite("ERROR: cannot allocate memory for FlagList\n");
    ReleaseTmpMem(theHeap);
    return (1);
  }
  for (i=0; i<npoints; i++) FlagList[i]=0;
  nlines=0; for (theLine=FirstLine(lgm_domain); theLine!=NULL; theLine=NextLine(lgm_domain)) nlines++;
  LinePtrList = (LGM_LINE **) GetTmpMem(theHeap,nlines*sizeof(LGM_LINE *));
  if (LinePtrList==NULL)
  {
    UserWrite("ERROR: cannot allocate memory for LineIDList\n");
    ReleaseTmpMem(theHeap);
    return (1);
  }
  for (i=0; i<nlines; i++) LinePtrList[i]=NULL;

  /* first save corners to list */
  for (theLine=FirstLine(lgm_domain); theLine!=NULL; theLine=NextLine(lgm_domain))
  {
    if (FlagList[LGM_LINE_BEGIN(theLine)]==0)
    {
      PointList[LGM_LINE_BEGIN(theLine)].position[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[0];
      PointList[LGM_LINE_BEGIN(theLine)].position[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[1];
      FlagList[LGM_LINE_BEGIN(theLine)] = 1;
    }
    else
    {
      if (PointList[LGM_LINE_BEGIN(theLine)].position[0] != LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[0] ||
          PointList[LGM_LINE_BEGIN(theLine)].position[1] != LGM_POINT_POS(LGM_LINE_POINT(theLine,0))[1]   )
      {
        UserWrite("ERROR: corner of line does not match previous position\n");
        return (1);
      }
    }
    if (FlagList[LGM_LINE_END(theLine)]==0)
    {
      PointList[LGM_LINE_END(theLine)].position[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,LGM_LINE_NPOINT(theLine)-1))[0];
      PointList[LGM_LINE_END(theLine)].position[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,LGM_LINE_NPOINT(theLine)-1))[1];
      FlagList[LGM_LINE_END(theLine)] = 1;
    }
    else
    {
      if (PointList[LGM_LINE_END(theLine)].position[0] != LGM_POINT_POS(LGM_LINE_POINT(theLine,LGM_LINE_NPOINT(theLine)-1))[0] ||
          PointList[LGM_LINE_END(theLine)].position[1] != LGM_POINT_POS(LGM_LINE_POINT(theLine,LGM_LINE_NPOINT(theLine)-1))[1]    )
      {
        UserWrite("ERROR: corner of line does not match previous position\n");
        return (1);
      }
    }
  }

  /* save line info */
  if (fprintf(stream,"# Line-Info\n")<0) return (1);
  act = 0;
  for (theLine=FirstLine(lgm_domain); theLine!=NULL; theLine=NextLine(lgm_domain))
  {
    if (LGM_LINE_ID(theLine)>=nlines || LGM_LINE_ID(theLine)<0)
    {
      UserWrite("ERROR: LineID out of range\n");
      ReleaseTmpMem(theHeap);
      return (1);
    }
    if (LinePtrList[LGM_LINE_ID(theLine)]==NULL)
      LinePtrList[LGM_LINE_ID(theLine)]=theLine;
    else
    {
      UserWrite("ERROR: LineID exists twice\n");
      ReleaseTmpMem(theHeap);
      return (1);
    }
  }
  for (i=0; i<nlines; i++)
  {
    theLine = LinePtrList[i];
    if (theLine==NULL)
    {
      UserWrite("ERROR: LinePtr not set\n");
      ReleaseTmpMem(theHeap);
      return (1);
    }
    if (fprintf(stream,"line %d: left=%d; right=%d; points: %d",(int)LGM_LINE_ID(theLine),(int)LGM_LINE_LEFT(theLine),(int)LGM_LINE_RIGHT(theLine),(int)LGM_LINE_BEGIN(theLine))<0) return (1);
    for (j=1; j<LGM_LINE_NPOINT(theLine)-1; j++)
    {
      while(FlagList[act] && act<npoints) act++;
      if (act>=npoints)
      {
        UserWrite("ERROR in FlagList\n");
        ReleaseTmpMem(theHeap);
        return (1);
      }
      if (fprintf(stream," %d",act)<0) return (1);
      PointList[act].position[0] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[0];
      PointList[act].position[1] = LGM_POINT_POS(LGM_LINE_POINT(theLine,j))[1];
      FlagList[act] = 1;
    }
    if (fprintf(stream," %d;\n",(int)LGM_LINE_END(theLine))<0) return (1);
  }

  /* save points */
  if (fprintf(stream,"\n# Point-Info\n")<0) return (1);
  for (i=0; i<npoints; i++)
  {
    if (FlagList[i]==0)
    {
      UserWrite("ERROR: FlagList-error, not all points are set correctly\n");
      ReleaseTmpMem(theHeap);
      return (1);
    }
    if (fprintf(stream,"%g %g;\n",(float)PointList[i].position[0],(float)PointList[i].position[1])<0)
    {
      UserWrite("ERROR: cannot save points\n");
      ReleaseTmpMem(theHeap);
      return (1);
    }
  }

  /* release tmp mem */
  ReleaseTmpMem(theHeap);

  /* close file */
  if (fclose(stream)==EOF) return (1);

  return (0);
}
