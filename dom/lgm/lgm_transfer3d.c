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
#include "ugdevices.h"
#include "lgm_domain.h"
#include "lgm_transfer.h"
#include "heaps.h"
#include "fifo.h"
#include "general.h"
#include "ng.h"


static FILE *stream;
static INT lgmdomainpathes_set;
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
    if (c!=' ' && c!='\n' && c!='\r' && c!='\t' && c!=';') break;
  }
  if (ungetc(c,stream)==EOF) return (1);

  return (0);
}


static int SkipEOL (void)
{
  int c;

  while (1)
  {
    c = fgetc(stream);
    if (c==EOF) return (1);
    if (c=='\n') break;
    if (c=='\r') break;
  }

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
static fpos_t UnitInfoFilepos;
static fpos_t fileposline;
static fpos_t filepossurface;

int LGM_ReadDomain (HEAP *Heap, char *filename, LGM_DOMAIN_INFO *domain_info, INT MarkKey)
{
  int i,i0,i1,i2;
  char buffer[256];
  fpos_t filepos_tmp;

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
    UserWriteF("  cannot open file %s\n",filename);
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

  /* if comment: read Unit-Info */
  if (ReadCommentLine("Unit-Info")) return (1);
  if (SkipBTN()) return (1);
  if (fgetpos(stream, &UnitInfoFilepos)) return (1);
  while (fscanf(stream,"unit %d",&i)==1)
    if (SkipEOL()) return (1);

  /* read Subdomain-Info if */
  if (fgetpos(stream, &filepos_tmp)) return (1);
  if (ReadCommentLine("Subdomain-Info"))
  {
    if (fsetpos(stream,&filepos_tmp)) return (1);
  }
  else
  {
    if (SkipBTN()) return (1);
    while (fscanf(stream,"subdomain %s",buffer)==1)
      if (SkipEOL()) return (1);
  }

  /* get number of subdomains, surfaces and lines */
  if (SkipBTN()) return (1);
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
  int i,n, old_i, error, line_id;

  error = 0;
  if(dummy == 0)
    if (fsetpos(stream, &fileposline))
      return (1);

  /* read line information */
  if (SkipBTN())
    return (1);
  if (fscanf(stream,"line %d:",&i)!=1)
    return (1);
  line_id = i;

  if (SkipBTN()) return (1);
  if (fscanf(stream,"points: %d",&i)!=1)
    return (1);
  old_i = i;
  n=0;
  line_info->point[n++] = i;
  while (1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"%d",&i)!=1)
      break;
    if(old_i==i)
      error++;
    line_info->point[n++] = i;
    old_i = i;
  }

  if(error>0)
    UserWriteF("%s %d\n", "Error in Line", line_id);

  return (0);
}

int LGM_ReadSubDomain (int subdom_i, LGM_SUBDOMAIN_INFO *subdom_info)
{
  int i,n,surface_i,i1,i2,i3,found,copy;
  fpos_t filepos_tmp;
  char buffer[256];


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


  /* scan Unit-Info */
  if (fgetpos(stream, &filepos_tmp)) return (1);
  if (fsetpos(stream, &UnitInfoFilepos)) return (1);
  found = 0;
  while(1)
  {
    copy = 0;
    if (fscanf(stream,"%s",buffer)!=1) break;
    if (strcmp(buffer,"unit")!=0) break;
    while (fscanf(stream," %d",&i)==1)
      if (i==subdom_i)
      {
        copy = 1;
        found++;
      }
    if (fscanf(stream,"%s",buffer)!=1) return (1);
    if (copy)
      strcpy(subdom_info->Unit,buffer);
  }
  if (found<1)
  {
    UserWriteF("ERROR: subdomain %d references no unit\n",(int)subdom_i);
    return (1);
  }
  if (found>1)
  {
    UserWriteF("ERROR: subdomain %d references more than 1 unit\n",(int)subdom_i);
    return (1);
  }
  if (fsetpos(stream, &filepos_tmp)) return (1);

  return (0);
}

int LGM_ReadPoints (LGM_POINT_INFO *lgm_point_info)
{
  int n;
  double f[3];

  if (SkipBTN())
    return (1);
  /* comment */
  if (ReadCommentLine("Point-Info")) return (1);

  n=0;
  while (1)
  {
    if (SkipBTN()) break;
    if (fscanf(stream,"%lf %lf %lf;",f,f+1,f+2)!=3) break;
    lgm_point_info[n].position[0] = f[0];
    lgm_point_info[n].position[1] = f[1];
    lgm_point_info[n].position[2] = f[2];
    n++;
  }

  if (fclose(stream)==EOF) return (1);
  return (0);
}

static int Compare_Triangles(LGM_TRIANGLE_INFO* t1, LGM_TRIANGLE_INFO* t2)
{
  int a, b, c, d, k, l;

  for(k=0; k<3; k++)
    for(l=0; l<3; l++)
    {
      a = t1->corner[(k+1)%3];
      b = t1->corner[(k+2)%3];
      c = t2->corner[(l+2)%3];
      d = t2->corner[(l+1)%3];
      if( (a==c) && (b==d) )
        return(-1);
      if( (a==d) && (b==c) )
        return(1);
    }

  return(0);
}

static void Change_Triangle(LGM_TRIANGLE_INFO *tri)
{
  int tmp;

  tmp = tri->corner[0];
  tri->corner[0] = tri->corner[1];
  tri->corner[1] = tmp;

  return;
}

/* OS_CHANGED: faster, no malloc() and even works for non-connected surfaces */
static int OrientateTriangles(LGM_SURFACE_INFO *surface_info, int tri_id, int *tr_used)
{
  HEAP *heap = theHeap;       /* global in this file */
  FIFO shell;
  LGM_TRIANGLE_INFO *triangles, *tri, *nb;
  int err, n, j, nb_id, flip, changed, MarkKey;
  void *buffer;

  /* the 1st triangle is assumed to be already orientated */
  assert(tr_used[tri_id]==1);

  changed = 0;
  triangles = surface_info->Triangle;
  n = surface_info->nTriangles;

  /* init FIFO */
  MarkTmpMem(heap,&MarkKey);
  buffer=(void *)GetTmpMem(heap,sizeof(LGM_TRIANGLE_INFO*)*n,MarkKey);
  assert(buffer!=NULL);
  fifo_init(&shell,buffer,sizeof(LGM_TRIANGLE_INFO*)*n);
  fifo_clear(&shell);

  /* put the first triangle as reference orientation into FIFO */
  err = fifo_in(&shell,(void *)&triangles[tri_id]);
  assert(!err);

  while (!fifo_empty(&shell))
  {
    /* get one triangle from FIFO and orientate its neighbours */
    tri = (LGM_TRIANGLE_INFO*) fifo_out(&shell);
    for(j=0; j<3; j++)
    {
      nb_id = tri->neighbor[j];
      if (nb_id != -1)                   /* neighbour exists */
      {
        if(!tr_used[nb_id])
        {
          nb = &triangles[nb_id];
          flip = Compare_Triangles(tri, nb);
          assert(flip != 0);                               /* Ensure they're neighbours! */
          if(flip==1)
          {
            Change_Triangle(nb);
            changed = 1;
          }
          assert(!fifo_full(&shell));
          err = fifo_in(&shell,(void *)nb);
          assert(!err);
          tr_used[nb_id] = 1;
        }
      }
    }
  }
  ReleaseTmpMem(heap,MarkKey);
  return changed;
}

static int Check_Orientation(LGM_SURFACE_INFO *surface_info, int id)
{
  HEAP *heap = theHeap;       /* global to this file */
  int i, nsfparts, n, changed, MarkKey;
  int *tr_used;

  n = surface_info->nTriangles;

  /* Init a flag-arry for already orientated triangles */
  MarkTmpMem(heap,&MarkKey);
  tr_used=(int *)GetTmpMem(heap,sizeof(int)*n,MarkKey);
  assert(tr_used!=NULL);

  for(i=0; i<n; i++)
    tr_used[i] = 0;

  i = nsfparts = 0;
  do
  {
    /* Take triangle i as reference, try to orientate the others */
    tr_used[i] = 1;
    changed = OrientateTriangles(surface_info,i,tr_used);

    /*
     * Check if all triangles are oriented now - otherwise we have another
     * non-connected part of the surface and we go on with the next non-
     * oriented triangle.
     */
    for (i=0; i<n; i++) if (!tr_used[i]) break;
    nsfparts++;
  } while (i<n);

  if (nsfparts > 1) UserWriteF("\n  * Surface %d has %d non-connected parts. *",
                               id, nsfparts);
  ReleaseTmpMem(heap,MarkKey);
  return changed;
}

static int Search_Neighbours(LGM_SURFACE_INFO *surface_info, int **point_list, int nPoints)
{
  int ni,i,j,k,l,nTriangle,corner_id;
  int a,b,c,d;

  nTriangle = surface_info->nTriangles;

  for(i=0; i<nTriangle; i++)
    for(j=0; j<3; j++)
      surface_info->Triangle[i].neighbor[j] = -1;

  for(i=0; i<nPoints; i++)
  {
    point_list[i][0] = 0;
    for(j=1; j<MAXTRIANGLES; j++)
      point_list[i][j] = -1;
  }
  for(i=0; i<nTriangle; i++)
  {
    for(j=0; j<3; j++)
    {
      corner_id = surface_info->Triangle[i].corner[j];
      point_list[corner_id][++point_list[corner_id][0]] = i;
    }
  }

  for(ni=0; ni<nPoints; ni++)
    for(i=1; i<=point_list[ni][0]; i++)
      for(j=1; j<=point_list[ni][0]; j++)
        if(i!=j)
          for(k=0; k<3; k++)
            for(l=0; l<3; l++)
            {
              a = surface_info->Triangle[point_list[ni][i]].corner[(k+1)%3];
              b = surface_info->Triangle[point_list[ni][i]].corner[(k+2)%3];
              c = surface_info->Triangle[point_list[ni][j]].corner[(l+2)%3];
              d = surface_info->Triangle[point_list[ni][j]].corner[(l+1)%3];
              if( ((a==c) && (b==d)) || ((a==d) && (b==c)) )
              {
                surface_info->Triangle[point_list[ni][i]].neighbor[k] = point_list[ni][j];
              }
            }
  return(0);
}

int LGM_ReadSurface (int dummy, LGM_SURFACE_INFO *surface_info)
{
  int i,k,n,i1,i2,i3,surface_id, error;

  error = 0;
  if(dummy == 0)
    if (fsetpos(stream, &filepossurface))
      return (1);

  if (fscanf(stream,"surface %d:",&i)!=1)
    return (1);
  surface_id = i;
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
  if((i1!=i2)&&(i1!=i3)&&(i2!=i3))
  {
    surface_info->Triangle[k].corner[0] = i1;
    surface_info->Triangle[k].corner[1] = i2;
    surface_info->Triangle[k].corner[2] = i3;
    surface_info->Triangle[k].neighbor[0] = 0;
    surface_info->Triangle[k].neighbor[1] = 0;
    surface_info->Triangle[k].neighbor[2] = 0;
  }
  else
  {
    error++;
    UserWriteF("%s %d %s %d\n", "Error in Surface", surface_id, "; triangle ", k+error);
    k--;
  }
  while (1)
  {
    if (SkipBTN())
      return (1);
    if (fscanf(stream,"%d %d %d;",&i1,&i2,&i3)!=3)
      break;
    if((i1!=i2)&&(i1!=i3)&&(i2!=i3))
    {
      k++;
      surface_info->Triangle[k].corner[0] = i1;
      surface_info->Triangle[k].corner[1] = i2;
      surface_info->Triangle[k].corner[2] = i3;
      surface_info->Triangle[k].neighbor[0] = 0;
      surface_info->Triangle[k].neighbor[1] = 0;
      surface_info->Triangle[k].neighbor[2] = 0;
    }
    else
    {
      error++;
      UserWriteF("%s %d %s %d\n", "Error in Surface", surface_id, "; triangle ", k+error);
    }
    if (SkipBTN())
      return (1);
  }

  surface_info->nTriangles = k+1;
  Search_Neighbours(surface_info, surface_info->point_list, surface_info->length);
  if(Check_Orientation(surface_info, surface_id))
    UserWriteF("Warning: Orientation of input triangles on surface %4d changed.\n",surface_id);

  return (0);
}

/****************************************************************************/
/*D
        LGM_WriteOpenFile - init

   SYNOPSIS:
   FILE *LGM_WriteOpenFile (char* name);

   PARAMETERS:

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

FILE *LGM_WriteOpenFile (char* name)
{
  FILE *stream;

  /* open file */
  if (lgmdomainpathes_set)
    stream = FileOpenUsingSearchPaths(name,"w","lgmdomainpathes");
  else
    stream = fileopen(name,"w");
  if (stream==NULL)
  {
    UserWriteF("cannot open file %s\n",name);
    return(NULL);
  }

  return (stream);
}

INT InitLGMTransfer (void)
{
  /* path to dir 'storidge' */
  lgmdomainpathes_set = 0;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"lgmdomainpathes")==0)
    lgmdomainpathes_set = 1;

  /* init NG */
  if (NG_Init ((int)lgmdomainpathes_set)) return (1);

  return (0);
}
