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

#include "ng.h"
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
static float **Line_coords;
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

static int CheckBndNode (BND_NODE *BndNode)
{
  if (BndNode->n_sp<=0) return (1);
  if (BndNode->n_lp>0 && BndNode->n_sp<2) return (1);
  return (0);
}

static int CheckInnNode (INNER_NODE *InnNode)
{
  return (0);
}

static int CheckElem (ELEMENT *Elem)
{
  int i;

  /* check element */
  if (Elem->subdom<=0) return (1);

  /* check faces */
  for (i=0; i<Elem->n_f; i++)
    if (Elem->face[i].n_c!=3)
      return (1);
  switch (Elem->n_c)
  {
  case 4 :
    if (Elem->n_f>4) return (1);
    break;
  case 5 :
    if (Elem->n_f>5) return (1);
    break;
  case 6 :
    if (Elem->n_f>5) return (1);
    break;
  case 8 :
    if (Elem->n_f>6) return (1);
    break;
  default :
    return (1);
  }

  return (0);
}

static int cmp_int (const void *p, const void *q)
{
  int i1,i2;

  i1 = *((int*)p);
  i2 = *((int*)q);
  if (i1<i2) return (-1);
  if (i2<i1) return (1);
  return (0);
}

int NP_ElemSideOnBnd (ELEMENT *Elem)
{
  int i,j,n,esob;
  int cof[4];

  esob=0;
  for (i=0; i<Elem->n_f; i++)
  {
    if (Elem->face[i].n_c!=3) continue;
    for (j=0; j<3; j++)
    {
      for (n=0; n<Elem->n_c; n++)
        if (Elem->face[i].c_id[j]==Elem->c_id[n])
        {
          cof[j]=n;
          break;
        }
      if (n==Elem->n_c) ngbreak();
    }
    qsort((void*)cof,3,sizeof(int),cmp_int);
    switch (Elem->n_c)
    {
    case 4 :
      if (cof[0]==0 && cof[1]==1 && cof[2]==2) esob |= 1;
      if (cof[0]==1 && cof[1]==2 && cof[2]==3) esob |= 2;
      if (cof[0]==0 && cof[1]==2 && cof[2]==3) esob |= 4;
      if (cof[0]==0 && cof[1]==1 && cof[2]==3) esob |= 8;
      break;
    case 5 :
      if (cof[0]==0 && cof[1]==1 && cof[2]==4) esob |= 2;
      if (cof[0]==1 && cof[1]==2 && cof[2]==4) esob |= 4;
      if (cof[0]==2 && cof[1]==3 && cof[2]==4) esob |= 8;
      if (cof[0]==0 && cof[1]==3 && cof[2]==4) esob |= 16;
      break;
    case 6 :
      if (cof[0]==0 && cof[1]==1 && cof[2]==2) esob |= 1;
      if (cof[0]==3 && cof[1]==4 && cof[2]==5) esob |= 16;
      break;
    }
  }

  return (esob);
}

#define NP_SPAT(x,y,z)                  ((x[1])*(y[2])-(x[2])*(y[1]))*z[0] + \
  ((x[2])*(y[0])-(x[0])*(y[2]))*z[1] + \
  ((x[0])*(y[1])-(x[1])*(y[0]))*z[2]

int OrientateElem (ELEMENT *Elem)
{
  int i;
  double sp;
  double p[8][3];

  for (i=0; i<Elem->n_c; i++)
    if (Elem->c_id[i]<Global_Mesh->nBndP)
    {
      p[i][0]=Global_Mesh->BndPosition[Elem->c_id[i]][0];
      p[i][1]=Global_Mesh->BndPosition[Elem->c_id[i]][1];
      p[i][2]=Global_Mesh->BndPosition[Elem->c_id[i]][2];
    }
    else
    {
      p[i][0]=Global_Mesh->InnPosition[Elem->c_id[i]-Global_Mesh->nBndP][0];
      p[i][1]=Global_Mesh->InnPosition[Elem->c_id[i]-Global_Mesh->nBndP][1];
      p[i][2]=Global_Mesh->InnPosition[Elem->c_id[i]-Global_Mesh->nBndP][2];
    }
  for (i=1; i<Elem->n_c; i++)
  {
    p[i][0]-=p[0][0];
    p[i][1]-=p[0][1];
    p[i][2]-=p[0][2];
  }
  switch (Elem->n_c)
  {
  case 4 :
    sp=NP_SPAT(p[1],p[2],p[3]);
    if (sp<0.0)
    {
      i=Elem->c_id[0]; Elem->c_id[0]=Elem->c_id[1]; Elem->c_id[1]=i;
#ifdef NP_Debug
      printf("\nSWITCH TET %d %d %d %d\n",Elem->c_id[0],Elem->c_id[1],Elem->c_id[2],Elem->c_id[3]);
#endif
    }
    break;
  case 5 :
    sp=NP_SPAT(p[1],p[2],p[4]);
    if (sp<0.0)
    {
      i=Elem->c_id[1]; Elem->c_id[1]=Elem->c_id[3]; Elem->c_id[3]=i;
#ifdef NP_Debug
      printf("\nSWITCH PYR %d %d %d %d %d\n",Elem->c_id[0],Elem->c_id[1],Elem->c_id[2],Elem->c_id[3],Elem->c_id[4]);
#endif
    }
    break;
  case 6 :
    sp=NP_SPAT(p[1],p[2],p[3]);
    if (sp<0.0)
    {
      i=Elem->c_id[0]; Elem->c_id[0]=Elem->c_id[1]; Elem->c_id[1]=i;
      i=Elem->c_id[3]; Elem->c_id[3]=Elem->c_id[4]; Elem->c_id[4]=i;
#ifdef NP_Debug
      printf("\nSWITCH PRI %d %d %d %d %d %d\n",Elem->c_id[0],Elem->c_id[1],Elem->c_id[2],Elem->c_id[3],Elem->c_id[4],Elem->c_id[5]);
#endif
    }
    break;
  case 8 :
    sp=NP_SPAT(p[1],p[2],p[4]);
    if (sp<0.0)
    {
      i=Elem->c_id[0]; Elem->c_id[0]=Elem->c_id[2]; Elem->c_id[2]=i;
      i=Elem->c_id[4]; Elem->c_id[4]=Elem->c_id[6]; Elem->c_id[6]=i;
#ifdef NP_Debug
      printf("\nSWITCH HEX %d %d %d %d %d %d %d %d\n",Elem->c_id[0],Elem->c_id[1],Elem->c_id[2],Elem->c_id[3],Elem->c_id[4],Elem->c_id[5],Elem->c_id[6],Elem->c_id[7]);
#endif
    }
    break;
  }

  return (0);
}

/****************************************************************************/
/*

    callback functions for parser

 */
/****************************************************************************/

int PutBndNode (BND_NODE *BndNode)
{
  int i,j,line_id;
  float *fp;

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
    Global_Mesh->BndP_nSurf[n_bn]=BndNode->n_sp;
    if (Global_Mesh->BndP_nLine[n_bn])
    {
      Global_Mesh->BndP_LineID[n_bn]=(int*)NG_MALLOC(Global_Heap,BndNode->n_lp*sizeof(int),Global_MarkKey);
      if (Global_Mesh->BndP_LineID[n_bn]==NULL) return (1);
      Global_Mesh->BndP_lcoord_left[n_bn]=(float*)NG_MALLOC(Global_Heap,BndNode->n_lp*sizeof(float),Global_MarkKey);
      if (Global_Mesh->BndP_lcoord_left[n_bn]==NULL) return (1);
      Global_Mesh->BndP_lcoord_right[n_bn]=(float*)NG_MALLOC(Global_Heap,BndNode->n_lp*sizeof(float),Global_MarkKey);
      if (Global_Mesh->BndP_lcoord_right[n_bn]==NULL) return (1);
    }
    else
    {
      Global_Mesh->BndP_LineID[n_bn]=NULL;
      Global_Mesh->BndP_lcoord_left[n_bn]=NULL;
      Global_Mesh->BndP_lcoord_right[n_bn]=NULL;
    }
    for (i=0; i<BndNode->n_lp; i++)
    {
      Global_Mesh->BndP_LineID[n_bn][i]=BndNode->lp[i].line_id;
      Global_Mesh->BndP_lcoord_left[n_bn][i]=NG_NOLEFT_COORD;
      Global_Mesh->BndP_lcoord_right[n_bn][i]=NG_NORIGHT_COORD;
      Line_npoints[BndNode->lp[i].line_id]++;
    }
    Global_Mesh->BndP_SurfID[n_bn]=(int*)NG_MALLOC(Global_Heap,BndNode->n_sp*sizeof(int),Global_MarkKey);
    if (Global_Mesh->BndP_SurfID[n_bn]==NULL) return (1);
    Global_Mesh->BndP_Cor_TriaID[n_bn]=(int*)NG_MALLOC(Global_Heap,BndNode->n_sp*sizeof(int),Global_MarkKey);
    if (Global_Mesh->BndP_Cor_TriaID[n_bn]==NULL) return (1);
    Global_Mesh->BndP_lcoord[n_bn]=(float**)NG_MALLOC(Global_Heap,BndNode->n_sp*sizeof(float*),Global_MarkKey);
    if (Global_Mesh->BndP_lcoord[n_bn]==NULL) return (1);
    fp=(float*)NG_MALLOC(Global_Heap,2*BndNode->n_sp*sizeof(float),Global_MarkKey);
    if (fp==NULL) return (1);
    for (i=0; i<BndNode->n_sp; i++)
    {
      Global_Mesh->BndP_SurfID[n_bn][i]=BndNode->sp[i].surf_id;
      Global_Mesh->BndP_Cor_TriaID[n_bn][i]=BndNode->sp[i].tri_id;
      Global_Mesh->BndP_lcoord[n_bn][i]=fp;
      Global_Mesh->BndP_lcoord[n_bn][i][0]=1.0-BndNode->sp[i].local[0]-BndNode->sp[i].local[1];
      Global_Mesh->BndP_lcoord[n_bn][i][1]=BndNode->sp[i].local[0];
      fp+=2;
    }
    Global_Mesh->BndPosition[n_bn]=(float*)NG_MALLOC(Global_Heap,3*sizeof(float),Global_MarkKey);
    if (Global_Mesh->BndPosition[n_bn]==NULL) return (1);
    Global_Mesh->BndPosition[n_bn][0]=BndNode->global[0];
    Global_Mesh->BndPosition[n_bn][1]=BndNode->global[1];
    Global_Mesh->BndPosition[n_bn][2]=BndNode->global[2];
    n_bn++;
    break;
  case 2 :
    for (i=0; i<BndNode->n_lp; i++)
    {
      line_id=BndNode->lp[i].line_id;
      Line_coords[line_id][Line_npoints[line_id]]=BndNode->lp[i].local;
      Line_npoints[line_id]++;
    }
    break;
  case 3 :
    for (i=0; i<BndNode->n_lp; i++)
    {
      line_id=BndNode->lp[i].line_id;
      for (j=0; j<Line_npoints[line_id]; j++)
      {
        if (Line_coords[line_id][j]<BndNode->lp[i].local)
          if (Line_coords[line_id][j]>Global_Mesh->BndP_lcoord_left[n_bn][i])
            Global_Mesh->BndP_lcoord_left[n_bn][i]=Line_coords[line_id][j];
        if (Line_coords[line_id][j]>BndNode->lp[i].local)
          if (Line_coords[line_id][j]<Global_Mesh->BndP_lcoord_right[n_bn][i])
            Global_Mesh->BndP_lcoord_right[n_bn][i]=Line_coords[line_id][j];
      }
    }
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
    Global_Mesh->InnPosition[n_in]=(double*)NG_MALLOC(Global_Heap,3*sizeof(double),Global_MarkKey);
    if (Global_Mesh->InnPosition[n_in]==NULL) return (1);
    Global_Mesh->InnPosition[n_in][0]=InnNode->global[0];
    Global_Mesh->InnPosition[n_in][1]=InnNode->global[1];
    Global_Mesh->InnPosition[n_in][2]=InnNode->global[2];
    n_in++;
    break;
  }

  return (0);
}

int PutElement (ELEMENT *Elem)
{
  int i,j,side;

  switch (mode)
  {
  case 0 :
    if (CheckElem(Elem)) return (1);
    if (Elem->subdom>subdom_max) subdom_max=Elem->subdom;
    break;
  case 1 :
    Global_Mesh->nSides[Elem->subdom]+=Elem->n_f;
    Global_Mesh->nElements[Elem->subdom]++;
    break;
  case 2 :
    if (OrientateElem(Elem)) return (1);
    Global_Mesh->Element_corners[Elem->subdom][Global_Mesh->nElements[Elem->subdom]]=Elem->n_c;
    for (i=0; i<Elem->n_f; i++)
    {
      Global_Mesh->Side_corners[Elem->subdom][Global_Mesh->nSides[Elem->subdom]]=Elem->face[i].n_c;
      Global_Mesh->nSides[Elem->subdom]++;
    }
    Global_Mesh->Element_SideOnBnd[Elem->subdom][Global_Mesh->nElements[Elem->subdom]]=NP_ElemSideOnBnd(Elem);
    Global_Mesh->nElements[Elem->subdom]++;
    break;
  case 3 :
    if (OrientateElem(Elem)) return (1);
    for (i=0; i<Elem->n_f; i++)
    {
      side=Global_Mesh->nSides[Elem->subdom];
      for (j=0; j<Global_Mesh->Side_corners[Elem->subdom][side]; j++)
        Global_Mesh->Side_corner_ids[Elem->subdom][side][j]=Elem->face[i].c_id[j];
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

  /* init */
  ng_abort=0;
  Global_Mesh=theMesh;
  Global_Heap=Heap;
  Global_MarkKey=MarkKey;

  /* open file */
  strcpy(ngname,name);
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
  theMesh->BndP_nSurf=(int*)NG_MALLOC(Heap,n_bn*sizeof(int),MarkKey);
  if (theMesh->BndP_nSurf==NULL) NG_HEAPFAULT;
  theMesh->BndP_nLine=(int*)NG_MALLOC(Heap,n_bn*sizeof(int),MarkKey);
  if (theMesh->BndP_nLine==NULL) NG_HEAPFAULT;
  theMesh->BndP_LineID=(int**)NG_MALLOC(Heap,n_bn*sizeof(int*),MarkKey);
  if (theMesh->BndP_LineID==NULL) NG_HEAPFAULT;
  theMesh->BndP_lcoord=(float***)NG_MALLOC(Heap,n_bn*sizeof(float**),MarkKey);
  if (theMesh->BndP_lcoord==NULL) NG_HEAPFAULT;
  theMesh->BndP_SurfID=(int**)NG_MALLOC(Heap,n_bn*sizeof(int*),MarkKey);
  if (theMesh->BndP_SurfID==NULL) NG_HEAPFAULT;
  theMesh->BndP_Cor_TriaID=(int**)NG_MALLOC(Heap,n_bn*sizeof(int*),MarkKey);
  if (theMesh->BndP_Cor_TriaID==NULL) NG_HEAPFAULT;
  theMesh->BndP_lcoord_left=(float**)NG_MALLOC(Heap,n_bn*sizeof(float*),MarkKey);
  if (theMesh->BndP_lcoord_left==NULL) NG_HEAPFAULT;
  theMesh->BndP_lcoord_right=(float**)NG_MALLOC(Heap,n_bn*sizeof(float*),MarkKey);
  if (theMesh->BndP_lcoord_right==NULL) NG_HEAPFAULT;
  theMesh->BndPosition=(float**)NG_MALLOC(Heap,n_bn*sizeof(float*),MarkKey);
  if (theMesh->BndPosition==NULL) NG_HEAPFAULT;
  theMesh->nbElements=NULL;
  theMesh->Element_SideOnBnd=NULL;
  Line_npoints=(int*)NG_MALLOC(Heap,(lineid_max+1)*sizeof(int),MarkKey);
  if (Line_npoints==NULL) NG_HEAPFAULT;
  Line_coords=(float**)NG_MALLOC(Heap,(lineid_max+1)*sizeof(float*),MarkKey);
  if (Line_coords==NULL) NG_HEAPFAULT;
  for (i=0; i<=lineid_max; i++) Line_npoints[i]=0;

  /* inner points */
  theMesh->nInnP=n_in;
  if (theMesh->nInnP>0)
  {
    theMesh->InnPosition=(double**)NG_MALLOC(Heap,n_in*sizeof(double*),MarkKey);
    if (theMesh->InnPosition==NULL) NG_HEAPFAULT;
  }
  else
    theMesh->InnPosition=NULL;

  /* elements */
  if (subdom_max<=0) {NG_Print(ERROR_PREFIX "nb of subdomains is 0\n"); return (1);}
  theMesh->nSubDomains=subdom_max;
  theMesh->nSides=(int*)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int),MarkKey);
  if (theMesh->nSides==NULL) NG_HEAPFAULT;
  for (i=0; i<=subdom_max; i++) theMesh->nSides[i]=0;
  theMesh->Side_corners=(int**)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int*),MarkKey);
  if (theMesh->Side_corners==NULL) NG_HEAPFAULT;
  theMesh->Side_corner_ids=(int***)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int**),MarkKey);
  if (theMesh->Side_corner_ids==NULL) NG_HEAPFAULT;
  theMesh->nElements=(int*)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int),MarkKey);
  if (theMesh->nElements==NULL) NG_HEAPFAULT;
  for (i=0; i<=subdom_max; i++) theMesh->nElements[i]=0;
  theMesh->Element_corners=(int**)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int*),MarkKey);
  if (theMesh->Element_corners==NULL) NG_HEAPFAULT;
  theMesh->Element_SideOnBnd=(int**)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int*),MarkKey);
  if (theMesh->Element_SideOnBnd==NULL) NG_HEAPFAULT;
  theMesh->Element_corner_ids=(int***)NG_MALLOC(Heap,(subdom_max+1)*sizeof(int**),MarkKey);
  if (theMesh->Element_corner_ids==NULL) NG_HEAPFAULT;

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
    Line_coords[i]=(float*)NG_MALLOC(Heap,Line_npoints[i]*sizeof(float),MarkKey);
    if (Line_coords[i]==NULL) NG_HEAPFAULT;
    Line_npoints[i]=0;
  }
  if (error) return (1);

  /* elements */
  error=0;
  for (i=1; i<=subdom_max; i++)
  {
    if (theMesh->nElements[i]<=0) {error=1; NG_Print(ERROR_PREFIX "no element in subdomain %d\n",i); continue;}
    if (theMesh->nSides[i]<=0) {error=1; NG_Print(ERROR_PREFIX "no side in subdomain %d\n",i); continue;}
    theMesh->Side_corners[i]=(int*)NG_MALLOC(Heap,theMesh->nSides[i]*sizeof(int),MarkKey);
    if (theMesh->Side_corners[i]==NULL) NG_HEAPFAULT;
    theMesh->Side_corner_ids[i]=(int**)NG_MALLOC(Heap,theMesh->nSides[i]*sizeof(int*),MarkKey);
    if (theMesh->Side_corner_ids[i]==NULL) NG_HEAPFAULT;
    theMesh->nSides[i]=0;
    theMesh->Element_corners[i]=(int*)NG_MALLOC(Heap,theMesh->nElements[i]*sizeof(int),MarkKey);
    if (theMesh->Element_corners[i]==NULL) NG_HEAPFAULT;
    theMesh->Element_SideOnBnd[i]=(int*)NG_MALLOC(Heap,theMesh->nElements[i]*sizeof(int),MarkKey);
    if (theMesh->Element_SideOnBnd[i]==NULL) NG_HEAPFAULT;
    theMesh->Element_corner_ids[i]=(int**)NG_MALLOC(Heap,theMesh->nElements[i]*sizeof(int*),MarkKey);
    if (theMesh->Element_corner_ids[i]==NULL) NG_HEAPFAULT;
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
      theMesh->Side_corner_ids[i][j]=(int*)NG_MALLOC(Heap,theMesh->Side_corners[i][j]*sizeof(int),MarkKey);
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
