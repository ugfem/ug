// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lgm_transfer.c												*/
/*																			*/
/* Purpose:   low level read routines from file								*/
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: klaus@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   11.09.96 begin												*/
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
#include "fileopen.h"
#include "ugdevices.h"
#include "defaults.h"
#include "misc.h"
#include "lgm_transfer.h"
#include "general.h"

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

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);

static FILE *stream;
static INT lgmdomainpathes_set;


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

static int SkipEOL (void)
{
  int c;

  while (1)
  {
    c = fgetc(stream);
    if (c==EOF) return (1);
    if (c=='\n') break;
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

static int nSubdomain, nLine;
static fpos_t filepos,filepos2,UnitInfoFilepos;
static HEAP *theHeap;

int LGM_ReadDomain (HEAP *Heap, char *filename, LGM_DOMAIN_INFO *domain_info, INT dummy)
{
  int i;
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

  /* Dimension */
  domain_info->Dimension = 2;

  /* convex */
  if (SkipBTN()) return (1);
  if (fscanf(stream,"convex = %d",&i)!=1) return (1);
  domain_info->Convex = i;

  /********************/
  /* skip unit-Info   */
  /********************/

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

  /* get number of subdomain and line */
  if (SkipBTN()) return (1);
  if (ReadCommentLine("Line-Info")) return (1);
  if (SkipBTN()) return (1);
  if (fgetpos(stream, &filepos)) return (1);
  domain_info->nSubDomain=domain_info->nPolyline=domain_info->nPoint=0;
  while (fscanf(stream,"line %d:",&i)==1)
  {
    if (SkipBTN()) return (1);
    if (fscanf(stream,"left=%d;",&i)!=1) return (1);
    domain_info->nSubDomain = MAX(domain_info->nSubDomain,i);
    if (SkipBTN()) return (1);
    if (fscanf(stream,"right=%d;",&i)!=1) return (1);
    domain_info->nSubDomain = MAX(domain_info->nSubDomain,i);

    if (SkipBTN()) return (1);
    if (fscanf(stream,"points: %d",&i)!=1) return (1);
    domain_info->nPoint = MAX(domain_info->nPoint,i);
    while (1)
    {
      if (SkipBTN()) return (1);
      if (fscanf(stream,"%d",&i)!=1) break;
      domain_info->nPoint = MAX(domain_info->nPoint,i);
    }
    domain_info->nPolyline++;
    if (SkipBTN()) return (1);
  }
  domain_info->nPoint++;
  nSubdomain = domain_info->nSubDomain;
  nLine = domain_info->nPolyline;

  return (0);
}

/****************************************************************************/
/*
   LGM_ReadSizes - reads sizes from file

   SYNOPSIS:
   int LGM_ReadSizes (LGM_SIZES *lgm_sizes);

   PARAMETERS:
   .  lgm_sizes -size information

   DESCRIPTION:
   function read a lgm-sizes from the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int LGM_ReadSizes (LGM_SIZES *lgm_sizes)
{
  int i,line_i;


  for (i=0; i<=nSubdomain; i++) lgm_sizes->Subdom_nLine[i] = 0;
  for (i=0; i<nLine; i++) lgm_sizes->Polyline_nPoint[i] = 0;
  line_i=0;
  if (fsetpos(stream, &filepos)) return (1);
  while (fscanf(stream,"line %d:",&i)==1)
  {
    if (SkipBTN()) return (1);
    if (fscanf(stream,"left=%d;",&i)!=1) return (1);
    lgm_sizes->Subdom_nLine[i]++;
    if (SkipBTN()) return (1);
    if (fscanf(stream,"right=%d;",&i)!=1) return (1);
    lgm_sizes->Subdom_nLine[i]++;

    if (SkipBTN()) return (1);
    if (fscanf(stream,"points: %d",&i)!=1) return (1);
    lgm_sizes->Polyline_nPoint[line_i]=1;
    while (1)
    {
      if (SkipBTN()) return (1);
      if (fscanf(stream,"%d",&i)!=1) break;
      lgm_sizes->Polyline_nPoint[line_i]++;
    }
    line_i++;
  }
  if (fsetpos(stream, &filepos)) return (1);

  return (0);
}

/****************************************************************************/
/*
   LGM_ReadLines - reads line information from file

   SYNOPSIS:
   int LGM_ReadLines (int i, LGM_LINE_INFO *line_info);

   PARAMETERS:
   .  dummy - dummy
   .  line_info - line information

   DESCRIPTION:
   function read a line-information from the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int LGM_ReadLines (int dummy, LGM_LINE_INFO *line_info)
{
  int i,n;

  /* read line information */
  if (SkipBTN()) return (1);
  if (fscanf(stream,"line %d:",&i)!=1) return (1);

  if (SkipBTN()) return (1);
  if (fscanf(stream,"left=%d",&i)!=1) return (1);
  line_info->left = i;

  if (SkipBTN()) return (1);
  if (fscanf(stream,"right=%d",&i)!=1) return (1);
  line_info->right = i;

  if (SkipBTN()) return (1);
  if (fscanf(stream,"points: %d",&i)!=1) return (1);
  n=0;
  line_info->point[n++] = i;
  while (1)
  {
    if (SkipBTN()) return (1);
    if (fscanf(stream,"%d",&i)!=1) break;
    line_info->point[n++] = i;
  }

  return (0);
}

/****************************************************************************/
/*
   LGM_ReadSubDomain - reads line information from file

   SYNOPSIS:
   int LGM_ReadSubDomain (int dummy, LGM_SUBDOMAIN_INFO *subdom_info);

   PARAMETERS:
   .  dummy - dummy
   .  subdom_info - subdomain information

   DESCRIPTION:
   function read a subdomain-information from the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int LGM_ReadSubDomain (int subdom_i, LGM_SUBDOMAIN_INFO *subdom_info)
{
  int i,n,line_i,found,copy;
  char buffer[256];

  /* read subdomain information */
  if (fsetpos(stream, &filepos)) return (1);
  line_i=0;
  n=0;
  if (SkipBTN()) return (1);
  while (fscanf(stream,"line %d:",&i)==1)
  {
    if (SkipBTN()) return (1);
    if (fscanf(stream,"left=%d;",&i)!=1) return (1);
    if (i==subdom_i) subdom_info->LineNumber[n++] = line_i;
    if (SkipBTN()) return (1);
    if (fscanf(stream,"right=%d;",&i)!=1) return (1);
    if (i==subdom_i) subdom_info->LineNumber[n++] = line_i;
    strcpy(subdom_info->Unit,"-");

    if (SkipBTN()) return (1);
    if (fscanf(stream,"points: %d",&i)!=1) return (1);
    while (1)
    {
      if (SkipBTN()) return (1);
      if (fscanf(stream,"%d",&i)!=1) break;
    }
    line_i++;
  }

  /* scan Unit-Info */
  if (fgetpos(stream, &filepos2)) return (1);
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
  if (fsetpos(stream, &filepos2)) return (1);

  return (0);
}

/****************************************************************************/
/*
   LGM_ReadPoints - reads all points from file

   SYNOPSIS:
   int LGM_ReadPoints (LGM_POINT_INFO *lgm_point_info);

   PARAMETERS:
   .  lgm_point_info - points

   DESCRIPTION:
   function reads all points from the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int LGM_ReadPoints (LGM_POINT_INFO *lgm_point_info)
{
  int n;
  float f[2];

  /* comment */
  if (ReadCommentLine("Point-Info")) return (1);

  n=0;
  while (1)
  {
    if (SkipBTN()) break;
    if (fscanf(stream,"%f %f;",f,f+1)!=2) break;
    lgm_point_info[n].position[0] = f[0];
    lgm_point_info[n].position[1] = f[1];
    n++;
  }

  if (fclose(stream)==EOF) return (1);

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

/****************************************************************************/
/*D
   InitLGMTransfer - init

   SYNOPSIS:
   INT InitLGMTransfer (void);

   PARAMETERS:

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

INT InitLGMTransfer (void)
{
  /* path to dir 'storidge' */
  lgmdomainpathes_set = 0;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"lgmdomainpathes")==0)
    lgmdomainpathes_set = 1;

  return (0);
}
