// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      gg3d.c                                                        */
/*                                                                          */
/* Purpose:   interface for the 3d grid generator netgen                            */
/*                                                                          */
/* Author:    Christian Wieners                                             */
/*			  Institut fuer Computeranwendungen III                         */
/*			  Universitaet Stuttgart			                            */
/*			  Pfaffenwaldring 27				                            */
/*			  70569 Stuttgart, Germany			                            */
/*			  email: ug@ica3.uni-stuttgart.de		                        */
/*									                                        */
/* History:   18 March 96 begin, ug version 3.2                             */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "compiler.h"
#include "devices.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"
#include "lgm_domain.h"
#include "domain.h"


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static INT nodeid;
static INT triangleid;
static INT left;
static INT right;
static double h_global;
static int LGM_DEBUG = 0;

static INT ntriangle;

static CoeffProcPtr LOCAL_H;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/*****************************************************************************/

#ifdef _NETGEN
extern int AddGeomPoint (int id, double x, double y, double z);
extern int StartSurfaceNetgen (double h, int smooth, int display, int D);
extern int InitSurfaceNetgen (char * rulefilename);
#endif
extern INT GetLocalKoord(LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local, DOUBLE *n);
extern INT Surface_Local2Global (LGM_SURFACE *theSurface, DOUBLE *global, DOUBLE *local);

static LGM_SURFACE *theSurface;
static HEAP *Heap;
static INT LGM_MarkKey;

int Allocate_Mem_Surfdisc(int npoints, int nelements)
{
  INT i;
  LGM_SURFACE_DISC(theSurface)->local = (DOUBLE **) GetTmpMem(Heap,(npoints+1)*sizeof(DOUBLE*),LGM_MarkKey);
  if(LGM_SURFACE_DISC(theSurface)->local==NULL)
    return(1);
  for(i=0; i<npoints; i++)
  {
    LGM_SURFACE_DISC(theSurface)->local[i] = (DOUBLE *)  GetTmpMem(Heap,3*sizeof(DOUBLE),LGM_MarkKey);
    if(LGM_SURFACE_DISC(theSurface)->local[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }
  LGM_SURFACE_DISC(theSurface)->triangle = (INT **) GetTmpMem(Heap,(nelements+1)*sizeof(INT*), LGM_MarkKey);
  if(LGM_SURFACE_DISC(theSurface)->triangle==NULL)
    return(1);
  for(i=0; i<nelements; i++)
  {
    LGM_SURFACE_DISC(theSurface)->triangle[i] = (INT *)  GetTmpMem(Heap,4*sizeof(INT),LGM_MarkKey);
    if(LGM_SURFACE_DISC(theSurface)->triangle[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }
  LGM_SURFACE_DISC(theSurface)->neighbour = (INT **) GetTmpMem(Heap,(nelements+1)*sizeof(INT*), LGM_MarkKey);
  if(LGM_SURFACE_DISC(theSurface)->neighbour==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<nelements; i++)
  {
    LGM_SURFACE_DISC(theSurface)->neighbour[i] = (INT *)  GetTmpMem(Heap,4*sizeof(INT),LGM_MarkKey);
    if(LGM_SURFACE_DISC(theSurface)->neighbour[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }
  return(0);
}

static INT AddPoint2Netgen (INT id, DOUBLE *global)
{

    #ifdef _NETGEN
  AddGeomPoint (nodeid, (double)global[0],(double)global[1],(double)global[2]);
    #endif

  return(0);
}


int AddInnerNode2ug (double x, double y, double z)
{
  DOUBLE global[3],local[2], n[3];

  global[0] = x;
  global[1] = y;
  global[2] = z;

  n[0] = n[1] = n[2] = 0.0;
  GetLocalKoord(theSurface,global,local, n);

  LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
  LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
  LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;

  global[0] = 0.0;
  global[1] = 0.0;
  global[2] = 0.0;
  Surface_Local2Global(theSurface, global, local);

  if(sqrt( (global[0]-x)*(global[0]-x) + (global[1]-y)*(global[1]-y) + (global[2]-z)*(global[2]-z) ) > 0.001)
  {
    printf("%f %f %f\n",x,y,z);
    printf("%f %f %f\n",global[0], global[1], global[2]);
  }
  if( sqrt( (global[0]-x)*(global[0]-x) + (global[1]-y)*(global[1]-y) + (global[2]-z)*(global[2]-z) ) > 0.00001 )
    printf("%s\n", "Warning in surface.c");
  assert( sqrt( (global[0]-x)*(global[0]-x) + (global[1]-y)*(global[1]-y) + (global[2]-z)*(global[2]-z) ) < 0.01 );

  return(0);
}

int AddSurfaceTriangle2ug (int node0, int node1, int node2)
{
  INT Id[3];

  Id[0] = node0;
  Id[1] = node1;
  Id[2] = node2;

  /*printf("%s %d %d %d\n","outputtriangle from netgen ",Id[0],Id[1],Id[2]);*/

  LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)),0) = Id[0];
  LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)),1) = Id[1];
  LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)),2) = Id[2];
  LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface))++;

  return(0);
}

#define MAX_T 30

static INT Search_Neighbour_Triangle(INT dummy)
{
  INT ni, i, j, k, l, ntriangle, npoint, a, b, c, d;
  INT **point_list, corner_id;

  ntriangle = LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface));
  npoint = LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface));

  for(i=0; i<ntriangle; i++)
    for(j=0; j<3; j++)
      LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), i, j) = -1;

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

  for(i=0; i<npoint; i++)
  {
    point_list[i][0] = 0;
    for(j=1; j<MAX_T; j++)
      point_list[i][j] = -1;
  }

  for(i=0; i<ntriangle; i++)
    for(j=0; j<3; j++)
    {
      corner_id = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),i,j);
      point_list[corner_id][++point_list[corner_id][0]] = i;
    }

  for(ni=0; ni<npoint; ni++)
    for(i=1; i<=point_list[ni][0]; i++)
      for(j=1; j<=point_list[ni][0]; j++)
        if(i!=j)
          for(k=0; k<3; k++)
            for(l=0; l<3; l++)
            {
              a = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),point_list[ni][i],(k+1)%3);
              b = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),point_list[ni][i],(k+2)%3);
              c = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),point_list[ni][j],(l+2)%3);
              d = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),point_list[ni][j],(l+1)%3);
              if( ((a==c)&&(b==d)) )
                LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), point_list[ni][i], k)
                  = point_list[ni][j];
            }


  return(0);
}

static INT Modify_Triangulation(INT dummy)
{
  INT i, j, k, ntriangle, n, m, flag1, flag2, triangle, neighbour, help, id[4], tr_id[3], n_id[3], trn_id[3], nn_id[3];
  INT *modify;

  ntriangle = LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface));
  modify = (INT *)GetTmpMem(Heap,ntriangle*sizeof(INT),LGM_MarkKey);
  if(modify==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<ntriangle; i++)
    modify[i] = 0;

  for(i=0; i<ntriangle; i++)
  {
    n = 0;
    for(j=0; j<3; j++)
      if(LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), i, j)==-1)
        n++;
      else
        flag1 = j;
    if(n==2)
    {
      m = 0;
      triangle = i;
      neighbour = LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, flag1);
      for(j=0; j<3; j++)
      {
        if(LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, j)==i)
          flag2 = j;
        if(LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, j)==-1)
          m++;
      }
      if(m==0)
        modify[neighbour]++;
    }
  }

  for(i=0; i<ntriangle; i++)
  {
    n = 0;
    for(j=0; j<3; j++)
      if(LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), i, j)==-1)
        n++;
      else
        flag1 = j;
    if(n==2)
    {
      m = 0;
      triangle = i;
      neighbour = LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, flag1);
      for(j=0; j<3; j++)
      {
        if(LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, j)==i)
          flag2 = j;
        if(LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, j)==-1)
          m++;
      }
      if((m==0)&&(modify[neighbour]==1))
      {
        /*				printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 0),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 1),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 2),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 0),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 1),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 2),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 0),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 1),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 2),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 0),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 1),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 2));*/

        for(k=0; k<3; k++)
        {
          tr_id[k] = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, (flag1 + 2 + k)%3 );
          n_id[k] = LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, (flag2 + 2 + k)%3 );
          trn_id[k] = LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, (flag1 + 2 + k)%3 );
          nn_id[k] = LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, (flag2 + 2 + k)%3 );
        }
        /*				printf("%d %d %d %d %d %d\n",tr_id[0],tr_id[1],tr_id[2],n_id[0],n_id[1],n_id[2]); */

        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 0) = tr_id[0];
        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 1) = tr_id[1];
        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 2) = n_id[1];
        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 0) = n_id[0];
        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 1) = n_id[1];
        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 2) = tr_id[1];
        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 0) = -1;
        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 1) = -1;
        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 2) = -1;
        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 0) = -1;
        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 1) = -1;
        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 2) = -1;

        /*				printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 0),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 1),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), triangle, 2),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 0),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 1),
                                        LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface), neighbour, 2),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 0),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 1),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), triangle, 2),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 0),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 1),
                                        LGM_SURFDISC_TRIANGLE_NEIGHBOUR(LGM_SURFACE_DISC(theSurface), neighbour, 2));*/
      }

    }
  }

  return(0);
}

INT GenerateSurfaceGrid (HEAP *theHeap, INT MarkKey, LGM_SURFACE *aSurface, DOUBLE h, INT smooth,INT display, INT D)
{
  INT sid,i;
  char rulefilename[128];
  DOUBLE **x;

  theSurface = aSurface;
  Heap = theHeap;
  LGM_MarkKey = MarkKey;

  ntriangle = 0;

  if (GetDefaultValue(DEFAULTSFILENAME,"netgentrianglerules",rulefilename))
    strcpy(rulefilename,rulefilename);

    #ifdef _NETGEN
  if (StartSurfaceNetgen(h,smooth,display, D)) return(1);
    #endif

  Search_Neighbour_Triangle(1);
  Modify_Triangulation(1);

  return(0);
}

int Get_Local_h(double *in, double *out)
{
  (*LOCAL_H)(in, out);
  return(0);
}

INT InitSurface(CoeffProcPtr Coeff)
{
  char rulefilename[128];
  if (GetDefaultValue(DEFAULTSFILENAME,"netgentrianglerules",rulefilename))
    strcpy(rulefilename,"triangle.rls");
  /*	LOCAL_H[0] = Coeff[coeff];*/
  LOCAL_H = Coeff;
    #ifdef _NETGEN
  InitSurfaceNetgen(rulefilename);
    #endif

}
