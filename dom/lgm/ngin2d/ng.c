// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gridt.c														*/
/*																			*/
/* Purpose:   grid representation                                           */
/*																			*/
/* Author:	  Klaus Johannsen                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "ng2d.h"
#ifdef __USE_IN_UG__
        #include "general.h"
        #include "fileopen.h"
        #include "defaults.h"
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

#undef NP_Debug

#define NG_NOLEFT_COORD                 -1.0
#define NG_NORIGHT_COORD                12345677890.0
#define NG_HEAPFAULT                    {NG_Print(ERROR_PREFIX "heap-fault\n"); return (1);}

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

extern FILE *ngin;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static int mode;
static int n_bn,n_in,n_el,lineid_max,subdom_max;
static int ng_abort;
static int *Line_npoints;
static LGM_MESH_INFO *Global_Mesh;
static HEAP *Global_Heap;
static Global_MarkKey;

#ifdef __USE_IN_UG__
static lgmdomainpathes_set;
#endif

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

void ngbreak (void);

/****************************************************************************/
/*D
   allocation of grid-objects

   D*/
/****************************************************************************/

/****************************************************************************/
/*

    check functions controlling nodes and elements

 */
/****************************************************************************/

static int CheckElem (NG_ELEMENT *Elem)
{
  int i;

  /* check element */
  if (Elem->subdom<=0) return (1);

  /* check sides */
  switch (Elem->n_c)
  {
  case 3 :
    if (Elem->n_s>3) return (1);
    break;
  case 4 :
    if (Elem->n_s>4) return (1);
    break;
  default :
    return (1);
  }

  return (0);
}

int NP_ElemSideOnBnd (NG_ELEMENT *Elem)
{
  int i,j,esob,c1,c2,ec1,ec2;

  esob=0;
  for (i=0; i<Elem->n_s; i++)
  {
    c1=Elem->side[i].c_id[0];
    c2=Elem->side[i].c_id[1];
    for (j=0; j<Elem->n_c; j++)
    {
      ec1=Elem->c_id[j];
      ec2=Elem->c_id[(j+1)%Elem->n_c];
      if ((c1==ec1 && c2==ec2) || (c1==ec2 && c2==ec1))
        esob|=1<<j;
    }
  }
  return (esob);
}

int OrientateElem (NG_ELEMENT *Elem)
{
  int i;
  double p[4][2];

  for (i=0; i<Elem->n_c; i++)
    if (Elem->c_id[i]<Global_Mesh->nBndP)
    {
      p[i][0]=Global_Mesh->BndPosition[Elem->c_id[i]][0];
      p[i][1]=Global_Mesh->BndPosition[Elem->c_id[i]][1];
    }
    else
    {
      p[i][0]=Global_Mesh->InnPosition[Elem->c_id[i]-Global_Mesh->nBndP][0];
      p[i][1]=Global_Mesh->InnPosition[Elem->c_id[i]-Global_Mesh->nBndP][1];
    }
  for (i=1; i<Elem->n_c; i++)
  {
    p[i][0]-=p[0][0];
    p[i][1]-=p[0][1];
  }
  if (p[1][0]*p[2][1]<p[1][1]*p[2][0])
  {
    i=Elem->c_id[0]; Elem->c_id[0]=Elem->c_id[2]; Elem->c_id[2]=i;
  }

  return(0);
}


/****************************************************************************/
/*

    callback functions for parser

 */
/****************************************************************************/

int PutBndNode (BND_NODE *BndNode)
{
  int i,j,line_id;

  switch (mode)
  {
  case 0 :
    n_bn++;
    for (i=0; i<BndNode->n_lp; i++)
      if (lineid_max<BndNode->lp[i].line_id)
        lineid_max=BndNode->lp[i].line_id;
    break;
  case 1 :
    Global_Mesh->BndP_nLine[n_bn]=BndNode->n_lp;
    if (Global_Mesh->BndP_nLine[n_bn])
    {
      Global_Mesh->BndP_LineID[n_bn]=(int*)NG_MALLOC(Global_Heap,BndNode->n_lp*sizeof(int),Global_MarkKey);
      if (Global_Mesh->BndP_LineID[n_bn]==NULL) return (1);
    }
    else
    {
      Global_Mesh->BndP_LineID[n_bn]=NULL;
    }
    for (i=0; i<BndNode->n_lp; i++)
    {
      Global_Mesh->BndP_LineID[n_bn][i]=BndNode->lp[i].line_id;
      Line_npoints[BndNode->lp[i].line_id]++;
    }
    Global_Mesh->BndP_lcoord[n_bn]=(float*)NG_MALLOC(Global_Heap,BndNode->n_lp*sizeof(float),Global_MarkKey);
    if (Global_Mesh->BndP_lcoord[n_bn]==NULL) return (1);
    for (i=0; i<BndNode->n_lp; i++)
      Global_Mesh->BndP_lcoord[n_bn][i]=BndNode->lp[i].local;
    Global_Mesh->BndPosition[n_bn]=(float*)NG_MALLOC(Global_Heap,2*sizeof(float),Global_MarkKey);
    if (Global_Mesh->BndPosition[n_bn]==NULL) return (1);
    Global_Mesh->BndPosition[n_bn][0]=BndNode->global[0];
    Global_Mesh->BndPosition[n_bn][1]=BndNode->global[1];
    n_bn++;
    break;
  case 2 :
    for (i=0; i<BndNode->n_lp; i++)
    {
      line_id=BndNode->lp[i].line_id;
      Line_npoints[line_id]++;
    }
    break;
  case 3 :
    n_bn++;
    break;
  }

  return (0);
}

int PutInnerNode (INNER_NODE *InnNode)
{
  switch (mode)
  {
  case 0 :
    n_in++;
    break;
  case 1 :
    Global_Mesh->InnPosition[n_in]=(double*)NG_MALLOC(Global_Heap,2*sizeof(double),Global_MarkKey);
    if (Global_Mesh->InnPosition[n_in]==NULL) return (1);
    Global_Mesh->InnPosition[n_in][0]=InnNode->global[0];
    Global_Mesh->InnPosition[n_in][1]=InnNode->global[1];
    n_in++;
    break;
  }

  return (0);
}

int PutElement (NG_ELEMENT *Elem)
{
  int i,j,side;

  switch (mode)
  {
  case 0 :
    if (CheckElem(Elem)) return (1);
    if (Elem->subdom>subdom_max) subdom_max=Elem->subdom;
    break;
  case 1 :
    Global_Mesh->nSides[Elem->subdom]+=Elem->n_s;
    Global_Mesh->nElements[Elem->subdom]++;
    break;
  case 2 :
    if (OrientateElem(Elem)) return (1);
    Global_Mesh->Element_corners[Elem->subdom][Global_Mesh->nElements[Elem->subdom]]=Elem->n_c;
    for (i=0; i<Elem->n_s; i++)
    {
      Global_Mesh->nSides[Elem->subdom]++;
    }
    Global_Mesh->Element_SideOnBnd[Elem->subdom][Global_Mesh->nElements[Elem->subdom]]=NP_ElemSideOnBnd(Elem);
    Global_Mesh->nElements[Elem->subdom]++;
    break;
  case 3 :
    if (OrientateElem(Elem)) return (1);
    for (i=0; i<Elem->n_s; i++)
    {
      side=Global_Mesh->nSides[Elem->subdom];
      for (j=0; j<2; j++)
        Global_Mesh->Side_corner_ids[Elem->subdom][side][j]=Elem->side[i].c_id[j];
      Global_Mesh->nSides[Elem->subdom]++;
    }
    for (i=0; i<Elem->n_c; i++)
      Global_Mesh->Element_corner_ids[Elem->subdom][Global_Mesh->nElements[Elem->subdom]][i]=Elem->c_id[i];
    Global_Mesh->nElements[Elem->subdom]++;
    break;
  }

  return (0);
}

void ngbreak (void)
{
  ng_abort=1;
}

int NG_ReadMesh (char *name, HEAP *Heap, LGM_MESH_INFO *theMesh, int MarkKey)
{
  int i,j,error;
  char ngname[128];
  char *p;
  char *tmp;

  /* init */
  ng_abort=0;
  Global_Mesh=theMesh;
  Global_Heap=Heap;
  Global_MarkKey=MarkKey;

  /* open file */
  strcpy(ngname,name);

  /* check if ngname ends with '.lgm'.  If so, remove '.lgm' */
  tmp = ngname + strlen(ngname) - 4;
  if ( strcmp(tmp,".lgm")==0 )
    ngname[strlen(ngname)-4] = '\0';

  strcat(ngname,".ng");
  NG_FOPEN(ngin,ngname);
  if (ngin==NULL) return (1);

  /* parse cycle 0 */
  mode=0;
  n_bn=n_in=lineid_max=subdom_max=0;
  while(!feof(ngin))
  {
    ngparse();
    if (ng_abort) return (1);
  }

  /* bnd points */
  if (n_bn<=0) {NG_Print(ERROR_PREFIX "nb of bnd points is 0\n"); return (1);}
  theMesh->nBndP=n_bn;
  theMesh->BndP_nLine=(int*)NG_MALLOC(Heap,n_bn*sizeof(int),MarkKey);                                             if (theMesh->BndP_nLine==NULL) NG_HEAPFAULT;
  theMesh->BndP_LineID=(int**)NG_MALLOC(Heap,n_bn*sizeof(int*),MarkKey);                                          if (theMesh->BndP_LineID==NULL) NG_HEAPFAULT;
  theMesh->BndP_lcoord=(float**)NG_MALLOC(Heap,n_bn*sizeof(float*),MarkKey);                                      if (theMesh->BndP_lcoord==NULL) NG_HEAPFAULT;
  theMesh->BndPosition=(float**)NG_MALLOC(Heap,n_bn*sizeof(float*),MarkKey);                                      if (theMesh->BndPosition==NULL) NG_HEAPFAULT;
  theMesh->nbElements=NULL;
  theMesh->Element_SideOnBnd=NULL;
  Line_npoints=(int*)NG_MALLOC(Heap,(lineid_max+1)*sizeof(int),MarkKey);                                          if (Line_npoints==NULL) NG_HEAPFAULT;
  for (i=0; i<=lineid_max; i++) Line_npoints[i]=0;

  /* inner points */
  theMesh->nInnP=n_in;
  if (theMesh->nInnP>0)
  {
    theMesh->InnPosition=(double**)NG_MALLOC(Heap,n_in*sizeof(double*),MarkKey);                    if (theMesh->InnPosition==NULL) NG_HEAPFAULT;
  }
  else
    theMesh->InnPosition=NULL;

  /* elements */
  if (subdom_max<=0) {NG_Print(ERROR_PREFIX "nb of subdomains is 0\n"); return (1);}
  theMesh->nSubDomains=subdom_max;
  theMesh->nSides=(int*)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int),MarkKey);                                       if (theMesh->nSides==NULL) NG_HEAPFAULT;
  for (i=0; i<=subdom_max; i++) theMesh->nSides[i]=0;
  theMesh->Side_corner_ids=(int***)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int**),MarkKey);          if (theMesh->Side_corner_ids==NULL) NG_HEAPFAULT;
  theMesh->nElements=(int*)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int),MarkKey);                            if (theMesh->nElements==NULL) NG_HEAPFAULT;
  for (i=0; i<=subdom_max; i++) theMesh->nElements[i]=0;
  theMesh->Element_corners=(int**)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int*),MarkKey);            if (theMesh->Element_corners==NULL) NG_HEAPFAULT;
  theMesh->Element_SideOnBnd=(int**)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int*),MarkKey);          if (theMesh->Element_SideOnBnd==NULL) NG_HEAPFAULT;
  theMesh->Element_corner_ids=(int***)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int**),MarkKey);       if (theMesh->Element_corner_ids==NULL) NG_HEAPFAULT;

  /* output */
  NG_Print("Parse: [0]");

  /* parse 1 */
  rewind(ngin);
  mode=1;
  n_bn=n_in=0;
  while(!feof(ngin))
  {
    ngparse();
    if (ng_abort) return(1);
  }

  /* bnd points */
  error=0;
  for (i=0; i<=lineid_max; i++)
  {
    if (Line_npoints[i]<=0) {error=1; NG_Print(ERROR_PREFIX "no points on line %d\n",i); continue;}
    Line_npoints[i]=0;
  }
  if (error) return (1);

  /* elements */
  error=0;
  for (i=1; i<=subdom_max; i++)
  {
    if (theMesh->nElements[i]<=0) {error=1; NG_Print(ERROR_PREFIX "no element in subdomain %d\n",i); continue;}
    if (theMesh->nSides[i]<=0) {error=1; NG_Print(ERROR_PREFIX "no side in subdomain %d\n",i); continue;}
    theMesh->Side_corner_ids[i]=(int**)NG_MALLOC(Heap,theMesh->nSides[i]*sizeof(int*),MarkKey);                                             if (theMesh->Side_corner_ids[i]==NULL) NG_HEAPFAULT;
    theMesh->nSides[i]=0;
    theMesh->Element_corners[i]=(int*)NG_MALLOC(Heap,theMesh->nElements[i]*sizeof(int),MarkKey);                                            if (theMesh->Element_corners[i]==NULL) NG_HEAPFAULT;
    theMesh->Element_SideOnBnd[i]=(int*)NG_MALLOC(Heap,theMesh->nElements[i]*sizeof(int),MarkKey);                                          if (theMesh->Element_SideOnBnd[i]==NULL) NG_HEAPFAULT;
    theMesh->Element_corner_ids[i]=(int**)NG_MALLOC(Heap,theMesh->nElements[i]*sizeof(int*),MarkKey);                                       if (theMesh->Element_corner_ids[i]==NULL) NG_HEAPFAULT;
    theMesh->nElements[i]=0;
  }
  if (error) return (1);

  /* output */
  NG_Print(" [1]");

  /* parse 2 */
  rewind(ngin);
  mode=2;
  while(!feof(ngin))
  {
    ngparse();
    if (ng_abort) return(1);
  }

  /* elements */
  for (i=1; i<=subdom_max; i++)
  {
    for (j=0; j<theMesh->nSides[i]; j++)
    {
      theMesh->Side_corner_ids[i][j]=(int*)NG_MALLOC(Heap,2*sizeof(int),MarkKey);
      if (theMesh->Side_corner_ids[i][j]==NULL) NG_HEAPFAULT;
    }
    theMesh->nSides[i]=0;
    for (j=0; j<theMesh->nElements[i]; j++)
    {
      theMesh->Element_corner_ids[i][j]=(int*)NG_MALLOC(Heap,theMesh->Element_corners[i][j]*sizeof(int),MarkKey);
      if (theMesh->Element_corner_ids[i][j]==NULL) NG_HEAPFAULT;
    }
    theMesh->nElements[i]=0;
  }

  NG_Print(" [2]");

  /* parse 3 */
  rewind(ngin);
  mode=3;
  n_bn=0;
  while(!feof(ngin))
  {
    ngparse();
    if (ng_abort) return(1);
  }
  NG_Print(" [3]\n");

  fclose(ngin);
  return (0);
}

#ifdef __USE_IN_UG__

int NG_Init (int domainpathes_set)
{
  lgmdomainpathes_set = domainpathes_set;

  return (0);
}

#else

main ()
{
  LGM_MESH_INFO theMesh;

  NG_ReadMesh ("test",NULL,&theMesh,0);

  return;
}

#endif
