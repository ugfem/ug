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
#include "gm.h"
#include "evm.h"
#include "misc.h"
#include "gginterface.h"
#include "general.h"
#include "debug.h"
#include "scan.h"

#include "gg3d.h"

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
static MULTIGRID *currMG;
static double h_global;
static INT subdomain;
static INT *transfer, *newId, *oldId, point, element;
static MESH *currmesh;

static INT *nInnP;
static DOUBLE ***Position;
static INT ***nbElement;

static INT nb_boundary_points;
static INT nb_boundary_points_subdom;
static INT nb_inner_points;
static DOUBLE scale[3][3];
static DOUBLE invscale[3][3];
static INT from_sub,  to_sub;
INT **point_list, **elem_list;

static INT GG3_DEBUG = 0;
static INT SAVE;
static CoeffProcPtr Coefficients[8];
static CoeffProcPtr LOCAL_H;

static INT GG3_MarkKey;

#define M_TIMES_M(A,B,C)                           {(C)[0][0] = (A)[0][0]*(B)[0][0]+(A)[0][1]*(B)[1][0]+(A)[0][2]*(B)[2][0];\
                                                    (C)[0][1] = (A)[0][0]*(B)[0][1]+(A)[0][1]*(B)[1][1]+(A)[0][2]*(B)[2][1];\
                                                    (C)[0][2] = (A)[0][0]*(B)[0][2]+(A)[0][1]*(B)[1][2]+(A)[0][2]*(B)[2][2];\
                                                    (C)[1][0] = (A)[1][0]*(B)[0][0]+(A)[1][1]*(B)[1][0]+(A)[1][2]*(B)[2][0];\
                                                    (C)[1][1] = (A)[1][0]*(B)[0][1]+(A)[1][1]*(B)[1][1]+(A)[1][2]*(B)[2][1];\
                                                    (C)[1][2] = (A)[1][0]*(B)[0][2]+(A)[1][1]*(B)[1][2]+(A)[1][2]*(B)[2][2];\
                                                    (C)[2][0] = (A)[2][0]*(B)[0][0]+(A)[2][1]*(B)[1][0]+(A)[2][2]*(B)[2][0];\
                                                    (C)[2][1] = (A)[2][0]*(B)[0][1]+(A)[2][1]*(B)[1][1]+(A)[2][2]*(B)[2][1];\
                                                    (C)[2][2] = (A)[2][0]*(B)[0][2]+(A)[2][1]*(B)[1][2]+(A)[2][2]*(B)[2][2];}

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/*****************************************************************************/

int AllMemInnerPoints(int npoints)
{
  INT i;

  point = 0;
  nInnP[subdomain] = npoints;
  Position[subdomain] = (DOUBLE **) GetTmpMem(MGHEAP(currMG),(npoints+1)*sizeof(DOUBLE*), GG3_MarkKey);
  if(Position==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  for(i=0; i<npoints; i++)
  {
    Position[subdomain][i] = (DOUBLE *) GetTmpMem(MGHEAP(currMG),3*sizeof(DOUBLE), GG3_MarkKey);
    if(Position[subdomain][i]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
  }
  return(0);
}

/****************************************************************************/
/*D
   AddInnerNode - add inner nodes in gg3d

   SYNOPSIS:
   int AddInnerNode (double x, double y, double z);

   PARAMETERS:
   .  x,y,z - coordiates of the node

   DESCRIPTION:
   This function is an interface function for the grid generator.
   It is called in gg3d and adds inner nodes in the node list in ug.

   RETURN VALUE:
   INT
   .n    nodeid if ok
   .n    -1 if error occured
   D*/
/****************************************************************************/

int AddInnerNode (double x, double y, double z)
{
  DOUBLE xc[3], nxc[3];

  PRINTDEBUG(dom,1,(" add inner node %4d %6.3lf %6.3lf %6.3lf\n",
                    nodeid,x,y,z));

  xc[0] = x;
  xc[1] = y;
  xc[2] = z;

  MM_TIMES_V_DIM(invscale, xc, nxc);

  /* write pointposition in the mesh structure */
  Position[subdomain][point][0] = nxc[0];
  Position[subdomain][point][1] = nxc[1];
  Position[subdomain][point][2] = nxc[2];
  point++;

  return(nodeid++);
}

int AllMemElements(int nelements)
{
  INT i;
  char buff[3], name[6];
  FILE *stream;

  if(SAVE)
  {
    name[0] = 'v';
    name[1] = 'o';
    name[2] = 'l';
    sprintf(buff,"%d",subdomain);
    name[3] = buff[0];
    name[4] = buff[1];
    name[5] = buff[2];

    stream = fopen(name,"w+");
    if (stream==NULL)
    {
      printf("%s\n", "cannot open file");
      return(1);
    }

    fprintf(stream, "%s\n", "vol_mesh");

    fprintf(stream,"%d\n",nelements);
    fclose(stream);
  }

  element = 0;
  currmesh->nElements[subdomain] = nelements;
  currmesh->Element_corners[subdomain] = (INT *) GetTmpMem(MGHEAP(currMG),(nelements+1)*sizeof(INT), GG3_MarkKey);
  if(currmesh->Element_corners[subdomain]==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  currmesh->Element_corner_ids[subdomain] = (INT **) GetTmpMem(MGHEAP(currMG),(nelements+1)*sizeof(INT*), GG3_MarkKey);
  if(currmesh->Element_corner_ids[subdomain]==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  return(0);
}

/****************************************************************************/
/*D
   AddElement  - add tetraheron in gg3d

   SYNOPSIS:
   int AddElement (int node0, int node1, int node2, int node3);

   PARAMETERS:
   .  node0,node1,node2,node3,node4,node5 - ids of the corner nodes

   DESCRIPTION:
   This function is an interface function for the grid generator.
   It is called in gg3d and adds elements in the element list in ug.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

int AddElement (int nnodes, int node0, int node1, int node2, int node3, int node4, int node5)
{
  INT Id[6], i;
  char buff[3], name[6];
  FILE *stream;
  ELEMENT *theElement;

  if(nnodes==4)
  {
    Id[0] = node1;
    Id[1] = node0;
    Id[2] = node2;
    Id[3] = node3;
  }
  if(nnodes==5)
  {
    Id[0] = node0;
    Id[1] = node3;
    Id[2] = node2;
    Id[3] = node1;
    Id[4] = node4;
  }
  if(nnodes==6)
  {
    Id[0] = node0;
    Id[1] = node2;
    Id[2] = node1;
    Id[3] = node3;
    Id[4] = node5;
    Id[5] = node4;
  }
  /* write element in the mesh structure */
  currmesh->Element_corner_ids[subdomain][element] = (INT *) GetTmpMem(MGHEAP(currMG),nnodes*sizeof(INT), GG3_MarkKey);
  if(currmesh->Element_corner_ids[subdomain][element]==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  currmesh->Element_corners[subdomain][element] = nnodes;
  for(i=0; i<nnodes; i++)
    if(Id[i]<nb_boundary_points_subdom)
      currmesh->Element_corner_ids[subdomain][element][i] = oldId[Id[i]];
    else
      currmesh->Element_corner_ids[subdomain][element][i] = nb_boundary_points
                                                            + nb_inner_points
                                                            + Id[i]
                                                            - nb_boundary_points_subdom;
  element++;

  return(0);
}

/****************************************************************************/
/*D
   AddBoundaryNode - append a node in the boundary node list of gg3d

   SYNOPSIS:
   static INT AddBoundaryNode (INT nodeid, DOUBLE *global);

   PARAMETERS:
   .  nodeid - id of the node
   .  global - coordiates of the node

   DESCRIPTION:
   This function is an interface function for the grid generator.
   It adds a node in the node list of surface nodes.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

static INT AddBoundaryNode (INT nodeid, DOUBLE *global)
{
  PRINTDEBUG(dom,1,(" add node   %4d    %6.3f %6.3f %6.3f\n",
                    nodeid,global[0],global[1],global[2]));

    #ifdef _NETGEN
  AddSurfaceNode (nodeid,
                  (double)global[0],(double)global[1],(double)global[2]);
    #endif

  return(0);
}

/****************************************************************************/
/*D
   AddBoundaryElement - append an element in the surface element list of gg3d

   SYNOPSIS:
   static INT AddBoundaryElement (INT n, INT *nodelist);

   PARAMETERS:
   .  n - number of corners
   .  nodelist - list or node ids

   DESCRIPTION:
   This function is an interface function for the grid generator.
   It adds a triangle described by its corner ids in the list of surface
   elements.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

static INT AddBoundaryElement (INT n, INT *nodelist, INT prism_flag)
{
  if (n != 3)
    return(1);

  PRINTDEBUG(dom,1,(" add triangle  %4d %4d %4d\n",
                    nodelist[0],nodelist[1],nodelist[2]));

    #ifdef _NETGEN
  AddSurfaceTriangle(nodelist[0],nodelist[1],nodelist[2], prism_flag);
    #endif

  return(0);
}

/* only for debug */
static INT Write_SurfaceMesh(MESH *mesh, MULTIGRID *theMG)
{
  char buff[3], name[10], name1[14];
  INT i, j, k, sid, nnodes, nsides, nodelist[4];
  FILE *stream;
  DOUBLE local[3], global[3], newglobal[3];
  NODE *theNode;
  VERTEX *theVertex;


  /* write mesh-information to file */

  stream = fopen("mesh","w+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }

  fprintf(stream, "%s\n", "mesh");

  fprintf(stream,"%d\n",(int)nb_boundary_points_subdom);
  for (i=0; i<mesh->nBndP; i++)
  {
    if(transfer[i])
    {
      if (BNDP_Global(mesh->theBndPs[i],global))
        return (1);

      MM_TIMES_V_DIM(scale, global, newglobal);

      fprintf(stream,"%lf %lf %lf\n",newglobal[0], newglobal[1], newglobal[2]);
    }
  }

  fprintf(stream,"%d\n",(int)mesh->nSides[subdomain]);

  for (i=subdomain; i<=subdomain; i++)
  {
    for (j=0; j<mesh->nSides[i]; j++)
    {
      for(k=0; k<mesh->Side_corners[i][j]; k++)
        nodelist[k] = newId[mesh->Side_corner_ids[i][j][k]];
      for (k=0; k<mesh->Side_corners[i][j]; k++)
        fprintf(stream,"%d ",nodelist[k]);
      fprintf(stream,"\n");
    }
  }

  fclose(stream);
  return(0);
}

static INT Write_VolumeMesh(MESH *mesh, MULTIGRID *theMG)
{
  char buff[3], name[6], name1[14];
  INT i, j, k, sid, nnodes, nsides, nodelist[6], id[6];
  FILE *stream;
  DOUBLE local[3], global[3];
  NODE *theNode;
  VERTEX *theVertex;

  name[0] = 'v';
  name[1] = 'o';
  name[2] = 'l';
  sprintf(buff,"%d",subdomain);
  name[3] = buff[0];
  name[4] = buff[1];
  name[5] = buff[2];

  /* write mesh-information to file */

  stream = fopen(name,"r+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }

  fprintf(stream, "%s\n", "vol_mesh");

  fprintf(stream,"%d\n",(int)mesh->nElements[subdomain]);

  for (i=subdomain; i<=subdomain; i++)
  {
    for (j=0; j<mesh->nElements[subdomain]; j++)
    {
      for(k=0; k<mesh->Element_corners[i][j]; k++)
        if(mesh->Element_corner_ids[i][j][k]<nb_boundary_points)
          id[k] = newId[mesh->Element_corner_ids[i][j][k]];
        else
          id[k] = mesh->Element_corner_ids[i][j][k]
                  - nb_boundary_points
                  - nb_inner_points
                  + nb_boundary_points_subdom;

      fprintf(stream,"%d\n",mesh->Element_corners[i][j]);
      for(k=0; k<mesh->Element_corners[i][j]; k++)
        fprintf(stream,"%d ",id[k]);
      fprintf(stream,"\n");
    }
  }

  fprintf(stream,"%d\n",(int)nb_boundary_points_subdom+(int)nInnP[subdomain]);
  for (i=0; i<mesh->nBndP; i++)
  {
    if(transfer[i])
    {
      if (BNDP_Global(mesh->theBndPs[i],global))
        return (1);
      fprintf(stream,"%20.16lf %20.16lf %20.16lf\n",global[0], global[1], global[2]);
    }
  }
  for (i=0; i<nInnP[subdomain]; i++)
    fprintf(stream,"%20.16lf %20.16lf %20.16lf\n",Position[subdomain][i][0], Position[subdomain][i][1], Position[subdomain][i][2]);

  fclose(stream);
  return(0);
}

static INT Get_NV(DOUBLE *n, INT i, INT j, INT k)
{
  DOUBLE n0[3], n1[3], n2[3];
  DOUBLE p0[3], p1[3], p2[3];

  p0[0] = Position[subdomain][i][0];
  p0[1] = Position[subdomain][i][1];
  p0[2] = Position[subdomain][i][2];
  p1[0] = Position[subdomain][j][0];
  p1[1] = Position[subdomain][j][1];
  p1[2] = Position[subdomain][j][2];
  p2[0] = Position[subdomain][k][0];
  p2[1] = Position[subdomain][k][1];
  p2[2] = Position[subdomain][k][2];
  V_DIM_SUBTRACT(p2, p0, n0);
  V_DIM_SCALE(1/sqrt(V_DIM_SCAL_PROD(n0, n0)), n0);
  V_DIM_SUBTRACT(p2, p1, n1);
  V_DIM_SCALE(1/sqrt(V_DIM_SCAL_PROD(n1, n1)), n1);

  V_DIM_VECTOR_PRODUCT(n0, n1, n);
  V_DIM_SCALE(1/sqrt(V_DIM_SCAL_PROD(n, n)), n);

  return(0);
}

static INT Get_Ang(DOUBLE *n0, DOUBLE *n1)
{
  DOUBLE s, angle;

  s = V_DIM_SCAL_PROD(n0, n1);
  s=MIN(1,s); s=MAX(-1,s);
  /*	printf("%lf %lf\n", s, acos(s));*/
  /*	if(s>1-0.00001)
                  angle = 3.141592654;
          else
                  if(s<-1+0.00001)
                          angle = 0.0;
                  else*/
  angle = 3.141592654-acos(s);

  return(angle);
}

static INT Angle_of_Element(INT *Id, DOUBLE *max_a, DOUBLE *min_a)
{
  INT i, j, k;
  DOUBLE n0[3], n1[3], n2[3], n3[3], s, angle;

  *max_a = 0.0;
  *min_a = 3.141592654;

  Get_NV(n0, Id[0], Id[1], Id[2]);
  Get_NV(n1, Id[2], Id[1], Id[3]);
  Get_NV(n2, Id[1], Id[0], Id[3]);
  Get_NV(n3, Id[2], Id[3], Id[0]);

  angle = Get_Ang(n0, n1);
  if(*max_a<angle)
    *max_a = angle;
  if(*min_a>angle)
    *min_a = angle;

  angle = Get_Ang(n0, n2);
  if(*max_a<angle)
    *max_a = angle;
  if(*min_a>angle)
    *min_a = angle;

  angle = Get_Ang(n0, n3);
  if(*max_a<angle)
    *max_a = angle;
  if(*min_a>angle)
    *min_a = angle;

  angle = Get_Ang(n1, n2);
  if(*max_a<angle)
    *max_a = angle;
  if(*min_a>angle)
    *min_a = angle;

  angle = Get_Ang(n1, n3);
  if(*max_a<angle)
    *max_a = angle;
  if(*min_a>angle)
    *min_a = angle;

  angle = Get_Ang(n2, n3);
  if(*max_a<angle)
    *max_a = angle;
  if(*min_a>angle)
    *min_a = angle;
  return(0);
}

static INT Read_VolumeMesh(MESH *mesh, MULTIGRID *theMG, INT MarkKey)
{
  char buff[3], name[6], name1[14];
  INT i, j, k, sid, buflen, iv;
  FILE *stream;
  char buffer[256];
  int nelements, npoints, d0, d1, d2, d3, id[6];
  double g0, g1, g2;
  DOUBLE max_angle, min_angle, max_a, min_a;

  buflen = 256;

  name[0] = 'v';
  name[1] = 'o';
  name[2] = 'l';
  sprintf(buff,"%d",subdomain);
  name[3] = buff[0];
  name[4] = buff[1];
  name[5] = buff[2];

  stream = fopen(name,"r+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }

  fscanf(stream,"%s\n", buffer);

  fscanf(stream,"%d\n",&nelements);

  mesh->nElements[subdomain] = nelements;
  mesh->Element_corners[subdomain] = (INT *) GetTmpMem(MGHEAP(currMG),(nelements+1)*sizeof(INT), MarkKey);
  if(mesh->Element_corners[subdomain]==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  mesh->Element_corner_ids[subdomain] = (INT **) GetTmpMem(MGHEAP(currMG),(nelements+1)*sizeof(INT*), MarkKey);
  if(mesh->Element_corner_ids[subdomain]==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  for (i=subdomain; i<=subdomain; i++)
  {
    for (j=0; j<mesh->nElements[subdomain]; j++)
    {
      fscanf(stream,"%d\n", &iv);
      mesh->Element_corners[i][j] = iv;

      for(k=0; k<mesh->Element_corners[i][j]; k++)
      {
        fscanf(stream,"%d ",&iv);
        id[k] = iv;
      }
      fscanf(stream,"\n");

      mesh->Element_corner_ids[subdomain][j] = (INT *) GetTmpMem(MGHEAP(currMG),mesh->Element_corners[i][j]*sizeof(INT), MarkKey);
      if(mesh->Element_corner_ids[subdomain][j]==NULL)
      {
        printf("%s\n", "Not enough memory");
        assert(0);
      }

      for(k=0; k<mesh->Element_corners[i][j]; k++)
      {
        if(id[k]<nb_boundary_points_subdom)
          mesh->Element_corner_ids[subdomain][j][k] = oldId[id[k]];
        else
          mesh->Element_corner_ids[subdomain][j][k] = nb_boundary_points
                                                      + nb_inner_points
                                                      + id[k]
                                                      - nb_boundary_points_subdom;
      }
    }
  }

  fscanf(stream,"%d\n",&npoints);

  nInnP[subdomain] = npoints - nb_boundary_points_subdom;

  Position[subdomain] = (DOUBLE **) GetTmpMem(MGHEAP(currMG),(nInnP[subdomain]+1)*sizeof(DOUBLE*), MarkKey);
  if(Position[subdomain]==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  for(i=0; i<nInnP[subdomain]; i++)
  {
    Position[subdomain][i] = (DOUBLE *) GetTmpMem(MGHEAP(currMG),3*sizeof(DOUBLE), MarkKey);
    if(Position[subdomain][i]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
  }
  for (i=0; i<nb_boundary_points_subdom; i++)
  {
    fscanf(stream,"%lg %lg %lg",&g0, &g1, &g2);
  }

  for (i=0; i<nInnP[subdomain]; i++)
  {
    fscanf(stream,"%lg %lg %lg",&g0, &g1, &g2);
    Position[subdomain][i][0] = g0;
    Position[subdomain][i][1] = g1;
    Position[subdomain][i][2] = g2;
  }
  fclose(stream);

  return(0);
}

static INT MAX_ELE = 75;

static INT Search_Neighbours(MULTIGRID *theMG, MESH *mesh, INT MarkKey)
{
  INT sid, np, i, j, k, l, m, n, a1, a2, a3, b1, b2, b3, npoints, corner_of_elem, flag;
  INT **pointlist, sides_i, sides_j, corners_i, corners_j, corner_i, corner_j, element_type_i, element_type_j;

  nbElement = (INT ***) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT**), MarkKey);
  if(nbElement==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  for(sid=from_sub; sid<=to_sub; sid++)
  {
    nbElement[sid] = (INT **) GetTmpMem(MGHEAP(currMG),(mesh->nElements[sid]+1)*sizeof(INT*), MarkKey);
    if(nbElement[sid]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
    for(i=0; i<mesh->nElements[sid]; i++)
    {
      nbElement[sid][i] = (INT *) GetTmpMem(MGHEAP(currMG),mesh->Element_corners[sid][i]*sizeof(INT), MarkKey);
      if(nbElement[sid][i]==NULL)
      {
        printf("%s\n", "Not enough memory");
        assert(0);
      }
    }
  }

  npoints =  nb_inner_points + nb_boundary_points;
  pointlist = (INT **) GetTmpMem(MGHEAP(currMG),(npoints+1)*sizeof(INT*), MarkKey);
  if(pointlist==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  for(i=0; i<npoints+1; i++)
  {
    pointlist[i] = (INT *) GetTmpMem(MGHEAP(currMG),MAX_ELE*sizeof(INT), MarkKey);
    if(pointlist[i]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }

  }
  if(GG3_DEBUG)
    for(sid=from_sub; sid<=to_sub; sid++)
      for(i=0; i<mesh->nElements[sid]; i++)
      {
        printf("%d ",i);
        for(j=0; j<mesh->Element_corners[sid][i]; j++)
          printf("%d ", mesh->Element_corner_ids[sid][i][j]);
        printf("\n");
      }

  for(sid=from_sub; sid<=to_sub; sid++)
  {
    for(i=0; i<mesh->nElements[sid]; i++)
      for(j=0; j<mesh->Element_corners[sid][i]; j++)
        nbElement[sid][i][j] = -1;
    for(i=0; i<npoints; i++)
      for(j=0; j<MAX_ELE; j++)
        pointlist[i][j] = 0;

    for(i=0; i<mesh->nElements[sid]; i++)
      for(j=0; j<mesh->Element_corners[sid][i]; j++)
      {
        corner_of_elem = mesh->Element_corner_ids[sid][i][j];
        pointlist[corner_of_elem][++pointlist[corner_of_elem][0]] = i;
      }

    if(GG3_DEBUG)
      for(np=0; np<npoints; np++)
      {
        printf("%d %d ", np, pointlist[np][0]);
        for(i=1; i<=pointlist[np][0]; i++)
          printf("%d ", pointlist[np][i]);
        printf("\n");
      }

    for(np=0; np<npoints; np++)
      if(pointlist[np][0]>1)
        for(i=1; i<=pointlist[np][0]; i++)
          for(j=1; j<=pointlist[np][0]; j++)
            if(i!=j)
            {
              element_type_i = mesh->Element_corners[sid][pointlist[np][i]];
              element_type_j = mesh->Element_corners[sid][pointlist[np][j]];
              sides_i = SIDES_OF_REF(element_type_i);
              sides_j = SIDES_OF_REF(element_type_i);

              for(k=0; k<sides_i; k++)
                for(l=0; l<sides_j; l++)
                {
                  corners_i = CORNERS_OF_SIDE_REF(element_type_i,k);
                  corners_j = CORNERS_OF_SIDE_REF(element_type_j,l);

                  if(corners_i==corners_j)
                  {
                    for(m=0; m<corners_i; m++)
                    {
                      flag = 0;
                      for(n=0; n<corners_j; n++)
                      {
                        corner_i = CORNER_OF_SIDE_REF(element_type_i,k,(  n)%corners_i);
                        corner_j = CORNER_OF_SIDE_REF(element_type_j,l,corners_j-1-(m+n)%corners_j);

                        if(
                          (       mesh->Element_corner_ids[sid][pointlist[np][i]][corner_i]
                                  ==  mesh->Element_corner_ids[sid][pointlist[np][j]][corner_j] )
                          )
                          flag++;
                      }
                      if(flag==corners_i)
                      {
                        nbElement[sid][pointlist[np][i]][k] = pointlist[np][j];
                        nbElement[sid][pointlist[np][j]][l] = pointlist[np][i];
                      }
                    }
                    /*	for(m=0;m<corners_i;m++)
                            {
                                    flag = 0;
                                    for(n=0;n<corners_j;n++)
                                    {
                                            corner_i = CORNER_OF_SIDE_REF(element_type_i,k,(  n)%corners_i);
                                            corner_j = CORNER_OF_SIDE_REF(element_type_j,l,(m+n)%corners_j);

                                            if(
                                            (	mesh->Element_corner_ids[sid][pointlist[np][i]][corner_i]
                                            ==  mesh->Element_corner_ids[sid][pointlist[np][j]][corner_j] )
                                            )
                                                    flag++;
                                    }
                                    if(flag==corners_i)
                                    {
                                            nbElement[sid][pointlist[np][i]][k] = pointlist[np][j];
                                            nbElement[sid][pointlist[np][j]][l] = pointlist[np][i];
                                    }
                            }*/
                  }
                }
            }
  }

  if(GG3_DEBUG)
    for(sid=from_sub; sid<=to_sub; sid++)
    {
      for(i=0; i<mesh->nElements[sid]; i++)
      {
        printf("%d%s ", i, ":");
        for(j=0; j<SIDES_OF_REF(mesh->Element_corners[sid][i]); j++)
          printf("%d ", nbElement[sid][i][j]);
        printf("\n");
      }
    }

  return(0);
}

#define Mult(vec1,vec2)         ( vec1[0] * vec2[0]     \
                                  + vec1[1] * vec2[1]     \
                                  + vec1[2] * vec2[2])

#define Lenght(vec)             sqrt(vec[0]*vec[0]      \
                                     +vec[1]*vec[1]  \
                                     +vec[2]*vec[2])

#define Minus(sol,vec1,vec2)    sol[0] = vec1[0] - vec2[0];     \
  sol[1] = vec1[1] - vec2[1];     \
  sol[2] = vec1[2] - vec2[2];

#define Scale(vec,l)            vec[0] = vec[0] * l;            \
  vec[1] = vec[1] * l;            \
  vec[2] = vec[2] * l;

#define Cross(vec,vec1,vec2)    vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]; \
  vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]; \
  vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

int GetNormalVector(double *p0, double *p1, double *p2, double *n)
{
  double n1[3], n2[3], l;

  Minus(n1, p1, p0);
  l = Lenght(n1);
  Scale(n1,1/l);
  Minus(n2, p2, p0);
  l = Lenght(n2);
  Scale(n2,1/l);
  Cross(n, n1, n2);
  l = Lenght(n);
  Scale(n,1/l);

  return(1);
}

static int Same_Plane(  double *p1_0, double *p1_1, double *p1_2,
                        double *p2_0, double *p2_1, double *p2_2)
{
  double n1[3], n2[3], sp;

  GetNormalVector(p1_0, p1_1, p1_2, n1);
  GetNormalVector(p2_0, p2_1, p2_2, n2);

  sp = Mult(n1, n2);

  if(sp<-0.999999)
    return(1);
  else
    return(0);

}

#define MAX_T 30

static INT Check_Volume(MESH *mesh, INT sid)
{
  INT ni, i, j, k, l, ntriangle, npoint, a, b, c, d, e, f, corner_id;
  DOUBLE g1[3], g2[3], g3[3], g4[3], g5[3], g6[3], ng1[3], ng2[3], ng3[3], ng4[3], ng5[3], ng6[3];

  ntriangle = mesh->nSides[sid];
  npoint = mesh->nBndP;

  /* first search neighbour triangle */
  for(i=0; i<ntriangle; i++)
    for(j=0; j<3; j++)
      elem_list[i][j] = -1;

  for(i=0; i<npoint; i++)
  {
    point_list[i][0] = 0;
    for(j=1; j<MAX_T; j++)
      point_list[i][j] = -1;
  }

  for(i=0; i<ntriangle; i++)
    for(j=0; j<3; j++)
    {
      corner_id = mesh->Side_corner_ids[sid][i][j];
      point_list[corner_id][++point_list[corner_id][0]] = i;
    }

  for(ni=0; ni<npoint; ni++)
    for(i=1; i<=point_list[ni][0]; i++)
      for(j=1; j<=point_list[ni][0]; j++)
        if(i!=j)
          for(k=0; k<3; k++)
            for(l=0; l<3; l++)
            {
              a = mesh->Side_corner_ids[sid][point_list[ni][i]][(k+1)%3];
              b = mesh->Side_corner_ids[sid][point_list[ni][i]][(k+2)%3];
              c = mesh->Side_corner_ids[sid][point_list[ni][j]][(l+2)%3];
              d = mesh->Side_corner_ids[sid][point_list[ni][j]][(l+1)%3];
              if( ((a==c)&&(b==d)) )
                elem_list[point_list[ni][i]][k] = point_list[ni][j];
            }

  /* triangles identical */
  for(i=0; i<ntriangle; i++)
  {
    a = mesh->Side_corner_ids[sid][i][0];
    b = mesh->Side_corner_ids[sid][i][1];
    c = mesh->Side_corner_ids[sid][i][2];
    for(j=0; j<3; j++)
    {
      d = mesh->Side_corner_ids[sid][elem_list[i][j]][0];
      e = mesh->Side_corner_ids[sid][elem_list[i][j]][1];
      f = mesh->Side_corner_ids[sid][elem_list[i][j]][2];
      if( ((b==d)&&(a==e)&&(c==f))
          ||      ((b==e)&&(a==f)&&(c==d))
          ||      ((b==f)&&(a==d)&&(c==e)) )
        return(1);
    }
  }
  /* triangles planar */
  for(i=0; i<ntriangle; i++)
  {
    BNDP_Global(mesh->theBndPs[mesh->Side_corner_ids[sid][i][0]],g1);
    BNDP_Global(mesh->theBndPs[mesh->Side_corner_ids[sid][i][1]],g2);
    BNDP_Global(mesh->theBndPs[mesh->Side_corner_ids[sid][i][2]],g3);
    MM_TIMES_V_DIM(scale, g1, ng1);
    MM_TIMES_V_DIM(scale, g2, ng2);
    MM_TIMES_V_DIM(scale, g3, ng3);
    for(j=0; j<3; j++)
    {
      BNDP_Global(mesh->theBndPs[mesh->Side_corner_ids[sid][elem_list[i][j]][0]],g4);
      BNDP_Global(mesh->theBndPs[mesh->Side_corner_ids[sid][elem_list[i][j]][1]],g5);
      BNDP_Global(mesh->theBndPs[mesh->Side_corner_ids[sid][elem_list[i][j]][2]],g6);
      MM_TIMES_V_DIM(scale, g4, ng4);
      MM_TIMES_V_DIM(scale, g5, ng5);
      MM_TIMES_V_DIM(scale, g6, ng6);
      if(Same_Plane(ng1, ng2, ng3, ng4, ng5, ng6))
        return(1);
    }
  }

  return(0);
}

static INT Allocate_Mem(MESH *mesh, INT from, INT to)
{
  int i, j, npoint, nelem;

  npoint = mesh->nBndP;
  nelem = 0;

  for(i=from; i<=to; i++)
    if(nelem < mesh->nSides[i])
      nelem = mesh->nSides[i];

  point_list = (INT **) GetTmpMem(MGHEAP(currMG),(npoint+1)*sizeof(INT*), GG3_MarkKey);
  if(point_list==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<npoint; i++)
  {
    point_list[i] = (INT *) GetTmpMem(MGHEAP(currMG),MAX_T*sizeof(INT), GG3_MarkKey);
    if(point_list[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }
  elem_list = (INT **) GetTmpMem(MGHEAP(currMG),(nelem+1)*sizeof(INT*), GG3_MarkKey);
  if(elem_list==NULL)
  {
    printf("%s\n","Not enough memory");
    assert(0);
  }
  for(i=0; i<nelem; i++)
  {
    elem_list[i] = (INT *) GetTmpMem(MGHEAP(currMG),3*sizeof(INT), GG3_MarkKey);
    if(elem_list[i]==NULL)
    {
      printf("%s\n","Not enough memory");
      assert(0);
    }
  }
  return(1);
}

static INT Write_Domain(MESH *mesh)
{
  INT i, j, nelem, npoint, sid, Id[6];
  FILE *file;
  DOUBLE global[3];

  nelem = 0;
  for (sid=from_sub; sid<=to_sub; sid++)
    nelem = nelem + mesh->nElements[sid];


  file = fopen("domain","w");
  fprintf(file, "%s\n", "volmesh");
  fprintf(file, "%d\n", nelem);

  for (sid=from_sub; sid<=to_sub; sid++)
  {
    for(i=0; i<mesh->nElements[sid]; i++)
    {
      fprintf(file, "%d\n",mesh->Element_corners[sid][i]);
      for(j=0; j<mesh->Element_corners[sid][i]; j++)
        fprintf(file, "%d ", mesh->Element_corner_ids[sid][i][j]);
      fprintf(file, "\n");
    }
  }

  fprintf(file, "%d\n", mesh->nBndP+mesh->nInnP);
  for(i=0; i<mesh->nBndP; i++)
  {
    if (BNDP_Global(mesh->theBndPs[i],global))
      return (1);
    fprintf(file, "%f %f %f\n", global[0], global[1], global[2]);
  }

  for(i=0; i<mesh->nInnP; i++)
  {
    for(j=0; j<3; j++)
      global[j] = mesh->Position[i][j];
    fprintf(file, "%f %f %f\n", global[0], global[1], global[2]);
  }

  return(0);
}

/****************************************************************************/
/*D
   GenerateGrid3d - call the netgen grid generator

   SYNOPSIS:
   INT GenerateGrid3d (MULTIGRID *theMG, MESH *mesh, DOUBLE h, INT smooth, INT coeff, DOUBLE *sc, INT from, INT to);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  mesh - pointer to mesh
   .  h - mesh size
   .  smoothing parameter
   .  coeff - number of coefficientfunction for gridsize
   .  sc - scaling
   .  from - subdomain
   .  to - subdomain
   DESCRIPTION:
   This function is an interface function for the grid generator.
   It calls the netgen grid genarator to compute an inner
   decomposition into tetrahedrons. It loads the surface triangulation
   given by the mesh structure.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

INT GenerateGrid3d (MULTIGRID *theMG, MESH *mesh, DOUBLE h, INT smooth,
                    INT display, INT coeff, INT from, INT to, INT prism, INT save,
                    INT argc, char **argv)
{
  NODE *theNode;
  VERTEX *theVertex;
  INT Scaling,sid,i,j, k, nodelist[6], Id[6], bnds_flag[6], k1, k2;
  char rulefilename[128];
  DOUBLE global[3], newglobal[3], det, vec[3], a1, a2, a3, n1[3], n2[3], n3[3], n[3], m[3], lam1, lam2, lam3, scal;
  DOUBLE T[3][3], T1[3][3], Diag[3][3], dummy[3][3];
  char buff[3], name[6];
  FILE *stream;
  ELEMENT *theElement;
  char buffer[512],scale_name[128];

  GG3_MarkKey = MG_MARK_KEY(theMG);
  SAVE = save;
  from_sub = from;
  to_sub = to;

  if (mesh == NULL)
    return(GM_OK);

  currMG = theMG;
  currmesh = mesh;
  if(h<=0.0)
  {
    /* get Coefficientfunctions */
    if (BVP_SetCoeffFct(MG_BVP(theMG),-1,Coefficients))
      return (0);
    LOCAL_H = Coefficients[coeff];
  }

  nInnP = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT), GG3_MarkKey);
  if(nInnP==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  Position = (DOUBLE ***) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(DOUBLE**), GG3_MarkKey);
  if(Position==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  mesh->nElements = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT), GG3_MarkKey);
  if(mesh->nElements==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  mesh->Element_corners = (INT **) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+2)*sizeof(INT*), GG3_MarkKey);
  if(mesh->Element_corners==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  mesh->Element_corner_ids = (INT ***) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT**), GG3_MarkKey);
  if(mesh->Element_corner_ids==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  transfer = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nBndP+1)*sizeof(INT), GG3_MarkKey);
  if(transfer==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  newId = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nBndP+1)*sizeof(INT), GG3_MarkKey);
  if(newId==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  oldId = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nBndP+1)*sizeof(INT), GG3_MarkKey);
  if(oldId==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }

  nb_boundary_points = theMG->nodeIdCounter;
  nb_inner_points = 0;

  Allocate_Mem(mesh, from, to);
  /* triangulate every subdomain */
  if (ReadArgvChar("Scaling",scale_name,argc,argv)) Scaling=0;
  else Scaling=1;
  for (sid=from; sid<=to; sid++)
  {
    subdomain = sid;

    if (Scaling)
    {
      sprintf(buffer,"%s%s%d",scale_name,"X",(int)sid);
      if (!GetStringValue(buffer,&a1))
      {
        sprintf(buffer,"%s%s%d",scale_name,"Y",(int)sid);
        if (GetStringValue(buffer,&a2)) assert(0);
        sprintf(buffer,"%s%s%d",scale_name,"Z",(int)sid);
        if (GetStringValue(buffer,&a3)) assert(0);
      }
      else
      {
        a1 = 0.0;
        a2 = 0.0;
        a3 = 1.0;
      }
    }
    else
    {
      a1 = 0.0;
      a2 = 0.0;
      a3 = 1.0;
    }
    vec[0] = a1;
    vec[1] = a2;
    vec[2] = a3;

    V_DIM_EUKLIDNORM(vec, lam1);
    V_DIM_SCALE(1/lam1,vec);
    V_DIM_COPY(vec, n1);
    lam2 = 1.0;
    lam3 = 1.0;

    V_DIM_SET(0.0,n2);
    V_DIM_SET(0.0,n3);
    V_DIM_SET(0.0,m);

    n[0] = 1.0; n[1] = 0.0; n[2] = 0.0;
    V_DIM_SCALAR_PRODUCT(n1, n, scal);
    V_DIM_COPY(n1, m);
    V_DIM_SCALE(scal, m);
    V_DIM_SUBTRACT(n, m, n2);
    V_DIM_EUKLIDNORM(n2, scal);
    if(scal<1e-4)
    {
      n[0] = 0.0; n[1] = 1.0; n[2] = 0.0;
      V_DIM_SCALAR_PRODUCT(n1, n, scal);
      V_DIM_COPY(n1, m);
      V_DIM_SCALE(scal, m);
      V_DIM_SUBTRACT(n, m, n2);
      V_DIM_EUKLIDNORM(n2, scal);
      if(scal<1e-4)
      {
        n[0] = 0.0; n[1] = 0.0; n[2] = 1.0;
        V_DIM_SCALAR_PRODUCT(n1, n, scal);
        V_DIM_COPY(n1, m);
        V_DIM_SCALE(scal, m);
        V_DIM_SUBTRACT(n, m, n2);
        V_DIM_EUKLIDNORM(n2, scal);
        if(scal<1e-4)
          assert(0);
      }
    }
    V_DIM_EUKLIDNORM(n2, scal);
    V_DIM_SCALE(1/scal,n2);

    V_DIM_VECTOR_PRODUCT(n1, n2, n3);
    V_DIM_EUKLIDNORM(n3, scal);
    V_DIM_SCALE(1/scal,n3);

    T[0][0] = n1[0]; T[0][1] = n2[0]; T[0][2] = n3[0];
    T[1][0] = n1[1]; T[1][1] = n2[1]; T[1][2] = n3[1];
    T[2][0] = n1[2]; T[2][1] = n2[2]; T[2][2] = n3[2];
    M_DIM_INVERT(T,T1,det);

    Diag[0][0] = lam1; Diag[0][1] = 0.0; Diag[0][2] = 0.0;
    Diag[1][0] = 0.0; Diag[1][1] = lam2; Diag[1][2] = 0.0;
    Diag[2][0] = 0.0; Diag[2][1] = 0.0; Diag[2][2] = lam3;

    M_TIMES_M(T, Diag, dummy);
    M_TIMES_M(dummy, T1, scale);
    M_DIM_INVERT(scale,invscale,det);

    if (GetDefaultValue(DEFAULTSFILENAME,"netgenrules",rulefilename))
      strcpy(rulefilename,"tetra.rls");
                #ifdef _NETGEN
    InitNetgen(rulefilename);
                #else
    PrintErrorMessage('E',"GenerateGrid3d","no netgen");
    return(1);
                #endif

    /* transfer points to netgen */
    for(i=0; i<mesh->nBndP; i++)
    {
      transfer[i] = 0;
      newId[i] = -1;
      oldId[i] = -1;
    }
    k = 0;
    for(i=0; i<mesh->nSides[sid]; i++)
    {
      for(j=0; j<mesh->Side_corners[sid][i]; j++)
      {
        if(transfer[mesh->Side_corner_ids[sid][i][j]]==0)
        {
          transfer[mesh->Side_corner_ids[sid][i][j]] = 1;
        }
      }
    }
    k = 0;
    for(i=0; i<mesh->nBndP; i++)
      if(transfer[i]==1)
      {
        newId[i] = k;
        oldId[k] = i;
        k++;
      }
    nb_boundary_points_subdom = k;
    /* insert boundary-points from mesh-structure to gridgenerator */
    for(i=0; i<mesh->nBndP; i++)
    {
      if(transfer[i])
      {
        if (BNDP_Global(mesh->theBndPs[i],global))
          return (1);
        MM_TIMES_V_DIM(scale, global, newglobal);
        if (AddBoundaryNode (newId[i],newglobal))
          return(1);
      }
    }
    /* transfer surface-triangles to netgen */
    for (i=0; i<mesh->nSides[sid]; i++)
    {
      for(j=0; j<mesh->Side_corners[sid][i]; j++)
        nodelist[j] = newId[mesh->Side_corner_ids[sid][i][j]];
      if (AddBoundaryElement (mesh->Side_corners[sid][i],nodelist, mesh->xy_Side[sid][i]))
        return(1);
    }

    if(SAVE)
    {
      name[0] = 'v';
      name[1] = 'o';
      name[2] = 'l';
      sprintf(buff,"%d",subdomain);
      name[3] = buff[0];
      name[4] = buff[1];
      name[5] = buff[2];

      stream = fopen(name,"r+");
      if (stream==NULL)
      {
        fclose(stream);
        printf("%s %d %s\n", "Subdomain ", subdomain, "not triangulated, do now");

        if(SAVE)
          Write_SurfaceMesh(mesh, theMG);
        if(Check_Volume(mesh, sid))
        {
          UserWriteF("%s %d\n", "Surfaces identical in subdomain", sid);
          UserWriteF("%s\n", "Check Surfaces");
          return(1);
        }
                                #ifdef _NETGEN
        if (StartNetgen(h,smooth,display, prism)) return(1);
                                #endif
        if(SAVE)
          Write_VolumeMesh(mesh, theMG);
      }
      else
      {
        fclose(stream);
        if(GG3_DEBUG) printf("%s %d\n", "Read Subdomain ", subdomain);
        Read_VolumeMesh(mesh, theMG, GG3_MarkKey);
      }
    }
    else
    {
      printf("%s %d %s\n", "Subdomain ", subdomain, "not triangulated, do now");
      if(Check_Volume(mesh, sid))
      {
        UserWriteF("%s %d\n", "Surfaces identical in subdomain", sid);
        UserWriteF("%s\n", "Check Surfaces");
        return(1);
      }
                        #ifdef _NETGEN
      if (StartNetgen(h,smooth,display, prism)) return(1);
                        #endif
    }
    nb_inner_points = nb_inner_points + nInnP[subdomain];
  }

  /* write inner points into the mesh */
  /* future work: change mesh structure */
  mesh->nInnP = 0;
  for (sid=1; sid<=mesh->nSubDomains; sid++)
    if(nInnP[sid]>0)
      mesh->nInnP = mesh->nInnP + nInnP[sid];
  mesh->Position = (DOUBLE **) GetTmpMem(MGHEAP(currMG), (mesh->nInnP+1)*sizeof(DOUBLE*), GG3_MarkKey);
  if(mesh->Position==NULL)
  {
    printf("%s\n", "Not enough memory");
    assert(0);
  }
  for(i=0; i<mesh->nInnP; i++)
  {
    mesh->Position[i] = (DOUBLE *) GetTmpMem(MGHEAP(currMG),3*sizeof(DOUBLE), GG3_MarkKey);
    if(mesh->Position[i]==NULL)
    {
      printf("%s\n", "Not enough memory");
      assert(0);
    }
  }
  k = 0;
  for (sid=from; sid<=to; sid++)
  {
    for(i=0; i<nInnP[sid]; i++)
    {
      for(j=0; j<3; j++)
        mesh->Position[k][j] = Position[sid][i][j];
      k++;
    }
  }

  if(GG3_DEBUG)
    printf("%s\n", "3d-gg fertig");
  /* search neighbor-elements per subdomain */
  /* for InsertElement */
  Search_Neighbours(theMG, mesh, GG3_MarkKey);

  /* insert points into the multigrid */
  for(i=0; i<mesh->nInnP; i++)
  {
    for(j=0; j<3; j++)
      global[j] = mesh->Position[i][j];
    if (InsertInnerNode(GRID_ON_LEVEL(theMG,0),global) == NULL)
      return(-1);
  }

  /* insert elements into the multigrid */
  k1 = k2 = 0;
  for (sid=from; sid<=to; sid++)
  {
    subdomain = sid;
    k1 = k2+1;
    for(i=0; i<mesh->nElements[sid]; i++)
    {
      for(j=0; j<mesh->Element_corners[sid][i]; j++)
      {
        Id[j] = mesh->Element_corner_ids[sid][i][j];
        if(GG3_DEBUG)
          printf("%d ", mesh->Element_corner_ids[sid][i][j]);
        if(nbElement[sid][i][j]==-1)
          bnds_flag[j] = 1;
        else
          bnds_flag[j] = 0;
      }
      if(GG3_DEBUG)
        printf("\n");
      theElement = InsertElementFromIDs(GRID_ON_LEVEL(currMG,0),mesh->Element_corners[sid][i],Id,bnds_flag);
      /*			theElement = InsertElementFromIDs(GRID_ON_LEVEL(currMG,0),mesh->Element_corners[sid][i],Id,NULL);*/
      if (theElement==NULL)
        return(1);
      SETSUBDOMAIN(theElement,subdomain);

    }
    k2 = k2 + mesh->nElements[sid];
    if(GG3_DEBUG)
      printf("%s %d %d %d\n", "subdomain ", sid, k1, k2);
  }
  if(SAVE)
    Write_Domain(mesh);
  return(0);
}

int Get_h (double *in, double *out)
{
  (*LOCAL_H)(in, out);
  return(0);
}
