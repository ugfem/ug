// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lgm_transferd.c												*/
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compiler.h"
#include "fileopen.h"
#include "defaults.h"
#include "misc.h"
#include "devices.h"
#include "lgm_domain.h"
#include "lgm_transfer.h"
#include "heaps.h"
#include "general.h"


static FILE *stream;
static INT lgmdomainpathes_set;
static INT LGM_DEBUG = 0;
static HEAP *theHeap;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

static int SkipBTN (void)
{
  int c;

  while (1)
  {
    c = fgetc(stream);
    if (c==EOF) return (1);
    if (c!=' ' && c!='\n' && c!='\t' && c!=';') break;
  }
  if (ungetc(c,stream)==EOF) return (1);

  return (0);
}

static int ReadCommentLine (char *comment)
{
  char buffer[256];

  if (SkipBTN()) return (1);
  if (fscanf(stream,"# %s",buffer)!=1) return (1);
  if (strcmp(comment,buffer)!=0) return (1);

  return (0);
}

/****************************************************************************/
/*
   LGM_ReadDomain - reads general domain information from file

   SYNOPSIS:
   int LGM_ReadDomain (char *filename, LGM_DOMAIN_INFO *domain_info);

   PARAMETERS:
   .  filename - name of file to read
   .  domain_info - information

   DESCRIPTION:
   function read a lgm-domain from the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

static int nSubdomain, nSurface, nLine, nPoint;
static fpos_t filepos;
static fpos_t fileposline;
static fpos_t filepossurface;

int LGM_ReadDomain (HEAP *Heap, char *filename, LGM_DOMAIN_INFO *domain_info)
{
  int i,i0,i1,i2;
  char buffer[256];

  /* store heapptr */
  if (Heap==NULL) return (1);
  theHeap = Heap;

  /* open file */
  if (lgmdomainpathes_set)
    stream = FileOpenUsingSearchPaths(filename,"r","lgmdomainpathes");
  else
    stream = fileopen(filename,"r");
  if (stream==NULL)
  {
    UserWriteF("cannot open file %s\n",filename);
    return(1);
  }

  /********************/
  /* read Domain-Info */
  /********************/

  /* comment */
  if (ReadCommentLine("Domain-Info")) return (1);

  /* domain-name */
  if (SkipBTN()) return (1);
  if (fscanf(stream,"name = %s",buffer)!=1) return (1);
  strcpy(domain_info->Name,buffer);

  /* problem name */
  if (SkipBTN()) return (1);
  if (fscanf(stream,"problemname = %s",buffer)!=1) return (1);
  strcpy(domain_info->ProblemName,buffer);

  domain_info->Dimension = 3;

  /* convex */
  if (SkipBTN()) return (1);
  if (fscanf(stream,"convex = %d",&i)!=1) return (1);
  domain_info->Convex = i;

  /* get number of subdomains, surfaces and lines */
  if (SkipBTN())
    return (1);
  if (ReadCommentLine("Line-Info"))
    return (1);
  if (SkipBTN()) return (1);
  if (fgetpos(stream, &fileposline))
    return (1);
  domain_info->nSubDomain = 0;
  domain_info->nSurface = 0;
  domain_info->nPolyline = 0;
  domain_info->nPoint = 0;

  while (fscanf(stream,"line %d:",&i)==1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"points: %d",&i)!=1)
      return (1);
    while (1)
    {
      if (SkipBTN()) return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
    }
    domain_info->nPolyline++;
    if (SkipBTN())
      return (1);
  }

  if (SkipBTN())
    return (1);
  if (ReadCommentLine("Surface-Info"))
    return (1);
  if (SkipBTN())
    return (1);
  if (fgetpos(stream, &filepossurface))
    return (1);
  while (fscanf(stream,"surface %d:",&i)==1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"left=%d;",&i)!=1)
      return (1);
    domain_info->nSubDomain = MAX(domain_info->nSubDomain,i);
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"right=%d;",&i)!=1)
      return (1);
    domain_info->nSubDomain = MAX(domain_info->nSubDomain,i);

    if (SkipBTN())
      return (1);
    if (fscanf(stream,"points: %d",&i)!=1)
      return (1);
    domain_info->nPoint = MAX(domain_info->nPoint,i);
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
      domain_info->nPoint = MAX(domain_info->nPoint,i);
    }
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"lines: %d",&i)!=1)
      return (1);
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
    }
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"triangles: %d %d %d",&i0, &i1, &i2)!=3)
      return (1);
    while (1)
    {
      if (SkipBTN()) return (1);
      if (fscanf(stream,"%d %d %d",&i0, &i1, &i2)!=3)
        break;
    }
    domain_info->nSurface++;
  }
  domain_info->nPoint++;
  nSubdomain = domain_info->nSubDomain;
  nSurface = domain_info->nSurface;
  nLine = domain_info->nPolyline;
  nPoint = domain_info->nPoint;
  return (0);
}

int LGM_ReadSizes (LGM_SIZES *lgm_sizes)
{
  int i,line_i,surface_i,i0,i1,i2;


  for (i=0; i<=nSubdomain; i++) lgm_sizes->Subdom_nSurf[i] = 0;
  for (i=0; i<=nSurface; i++) lgm_sizes->Surf_nPolyline[i] = 0;
  for (i=0; i<=nSurface; i++) lgm_sizes->Surf_nTriangle[i] = 0;
  for (i=0; i<=nSurface; i++) lgm_sizes->Surf_nPoint[i] = 0;
  for (i=0; i<=nLine; i++) lgm_sizes->Polyline_nPoint[i] = 0;

  line_i=0;
  if (fsetpos(stream, &fileposline))
    return (1);
  while (fscanf(stream,"line %d:",&i)==1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"points: %d",&i)!=1)
      return (1);
    lgm_sizes->Polyline_nPoint[line_i]=1;
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
      lgm_sizes->Polyline_nPoint[line_i]++;
    }
    line_i++;
  }

  if (SkipBTN())
    return (1);
  surface_i=0;
  if (fsetpos(stream, &filepossurface))
    return (1);
  while (fscanf(stream,"surface %d:",&i)==1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"left=%d;",&i)!=1)
      return (1);
    lgm_sizes->Subdom_nSurf[i]++;
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"right=%d;",&i)!=1) return (1);
    lgm_sizes->Subdom_nSurf[i]++;

    if (SkipBTN())
      return (1);
    if (fscanf(stream,"points: %d",&i)!=1) return (1);
    lgm_sizes->Surf_nPoint[surface_i]=1;
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
      lgm_sizes->Surf_nPoint[surface_i]++;
    }

    if (SkipBTN())
      return (1);
    if (fscanf(stream,"lines: %d",&i)!=1) return (1);
    lgm_sizes->Surf_nPolyline[surface_i]=1;
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
      lgm_sizes->Surf_nPolyline[surface_i]++;
    }

    if (SkipBTN())
      return (1);
    if (fscanf(stream,"triangles: %d %d %d;",&i0, &i1, &i2)!=3) return (1);
    lgm_sizes->Surf_nTriangle[surface_i]=1;
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d %d %d",&i0, &i1, &i2)!=3)
        break;
      lgm_sizes->Surf_nTriangle[surface_i]++;
    }
    surface_i++;
  }

  return (0);
}

int LGM_ReadLines (int dummy, LGM_LINE_INFO *line_info)
{
  int i,n;

  if(dummy == 0)
    if (fsetpos(stream, &fileposline))
      return (1);

  /* read line information */
  if (SkipBTN())
    return (1);
  if (fscanf(stream,"line %d:",&i)!=1)
    return (1);

  if (SkipBTN()) return (1);
  if (fscanf(stream,"points: %d",&i)!=1)
    return (1);
  n=0;
  line_info->point[n++] = i;
  while (1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"%d",&i)!=1)
      break;
    line_info->point[n++] = i;
  }

  return (0);
}

int LGM_ReadSubDomain (int subdom_i, LGM_SUBDOMAIN_INFO *subdom_info)
{
  int i,n,surface_i,i1,i2,i3;

  /* read subdomain information */
  if (fsetpos(stream, &filepossurface))
    return (1);
  surface_i=0;
  n=0;
  if (SkipBTN()) return (1);
  while (fscanf(stream,"surface %d:",&i)==1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"left=%d;",&i)!=1)
      return (1);
    if (i==subdom_i)
      subdom_info->SurfaceNumber[n++] = surface_i;
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"right=%d;",&i)!=1) \
      return (1);
    if (i==subdom_i)
      subdom_info->SurfaceNumber[n++] = surface_i;

    if (SkipBTN())
      return (1);
    if (fscanf(stream,"points: %d",&i)!=1) return (1);
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
    }

    if (SkipBTN())
      return (1);
    if (fscanf(stream,"lines: %d",&i)!=1) return (1);
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d",&i)!=1)
        break;
    }

    if (SkipBTN())
      return (1);
    if (fscanf(stream,"triangles: %d %d %d",&i1, &i2, &i3)!=3) return (1);
    while (1)
    {
      if (SkipBTN())
        return (1);
      if (fscanf(stream,"%d %d %d",&i1, &i2, &i3)!=3)
        break;
    }

    surface_i++;
  }

  return (0);
}

int LGM_ReadPoints (LGM_POINT_INFO *lgm_point_info)
{
  int n;
  float f[3];

  if (SkipBTN())
    return (1);
  /* comment */
  if (ReadCommentLine("Point-Info")) return (1);

  n=0;
  while (1)
  {
    if (SkipBTN()) break;
    if (fscanf(stream,"%f %f %f;",f,f+1,f+2)!=3) break;
    lgm_point_info[n].position[0] = f[0];
    lgm_point_info[n].position[1] = f[1];
    lgm_point_info[n].position[2] = f[2];
    n++;
  }

  if (fclose(stream)==EOF) return (1);
  return (0);
}

int LGM_ReadSurface (int dummy, LGM_SURFACE_INFO *surface_info)
{
  int i,k,n,i1,i2,i3;

  if(dummy == 0)
    if (fsetpos(stream, &filepossurface))
      return (1);

  if (fscanf(stream,"surface %d:",&i)!=1)
    return (1);
  if (SkipBTN())
    return (1);
  if (fscanf(stream,"left=%d;",&i)!=1)
    return (1);
  surface_info->left = i;
  if (SkipBTN())
    return (1);
  if (fscanf(stream,"right=%d",&i)!=1)
    return (1);
  surface_info->right = i;
  if (SkipBTN())
    return (1);

  if (SkipBTN())
    return (1);
  if (fscanf(stream,"points: %d",&i)!=1) return (1);
  n=0;
  surface_info->point[n++] = i;
  while (1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"%d",&i)!=1)
      break;
    surface_info->point[n++] = i;
  }

  if (SkipBTN())
    return (1);
  if (fscanf(stream,"lines: %d",&i)!=1) return (1);
  n=0;
  surface_info->line[n++] = i;
  while (1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"%d",&i)!=1)
      break;
    surface_info->line[n++] = i;
  }

  if (fscanf(stream,"triangles: %d %d %d;",&i1,&i2,&i3)!=3)
    return (1);
  k = 0;
  surface_info->Triangle[k].corner[0] = i1;
  surface_info->Triangle[k].corner[1] = i2;
  surface_info->Triangle[k].corner[2] = i3;
  surface_info->Triangle[k].neighbor[0] = 0;
  surface_info->Triangle[k].neighbor[1] = 0;
  surface_info->Triangle[k].neighbor[2] = 0;

  while (1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"%d %d %d;",&i1,&i2,&i3)!=3)
      break;
    k++;
    surface_info->Triangle[k].corner[0] = i1;
    surface_info->Triangle[k].corner[1] = i2;
    surface_info->Triangle[k].corner[2] = i3;
    surface_info->Triangle[k].neighbor[0] = 0;
    surface_info->Triangle[k].neighbor[1] = 0;
    surface_info->Triangle[k].neighbor[2] = 0;
    if (SkipBTN())
      return (1);
  }

  return (0);
}

INT InitLGMTransfer (void)
{
  /* path to dir 'storidge' */
  lgmdomainpathes_set = 0;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"lgmdomainpathes")==0)
    lgmdomainpathes_set = 1;

  return (0);
}
