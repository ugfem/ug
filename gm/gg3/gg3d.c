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

static INT GG3_DEBUG = 0;

static CoeffProcPtr Coefficients[8];
static CoeffProcPtr LOCAL_H;

static INT GG3_MarkKey;

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
  for(i=0; i<npoints; i++)
    Position[subdomain][i] = (DOUBLE *) GetTmpMem(MGHEAP(currMG),3*sizeof(DOUBLE), GG3_MarkKey);
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
  DOUBLE_VECTOR xc;

  PRINTDEBUG(dom,1,(" add inner node %4d %6.3lf %6.3lf %6.3lf\n",
                    nodeid,x,y,z));

  xc[0] = x;
  xc[1] = y;
  xc[2] = z;

  /*	if (InsertInnerNode(GRID_ON_LEVEL(currMG,0),xc) == NULL)
            return(-1);*/

  /* write pointposition in the mesh structure */
  Position[subdomain][point][0] = x;
  Position[subdomain][point][1] = y;
  Position[subdomain][point][2] = z;
  point++;

  return(nodeid++);
}

int AllMemElements(int nelements)
{
  INT i;
  char buff[3], name[6];
  FILE *stream;

  if(GG3_DEBUG)
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
  currmesh->Element_corner_ids[subdomain] = (INT **) GetTmpMem(MGHEAP(currMG),(nelements+1)*sizeof(INT*), GG3_MarkKey);
}

/****************************************************************************/
/*D
   AddTetrahedron  - add tetraheron in gg3d

   SYNOPSIS:
   int AddTetrahedron (int node0, int node1, int node2, int node3);

   PARAMETERS:
   .  node0,node1,node2,node3 - ids of the corner nodes

   DESCRIPTION:
   This function is an interface function for the grid generator.
   It is called in gg3d and adds elements in the element list in ug.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

int AddTetrahedron (int node0, int node1, int node2, int node3)
{
  INT Id[4], i;
  char buff[3], name[6];
  FILE *stream;
  ELEMENT *theElement;

  Id[0] = node1;                /* das gehoert so ! */
  Id[1] = node0;
  Id[2] = node2;
  Id[3] = node3;

  if(GG3_DEBUG)
  {
    name[0] = 'v';
    name[1] = 'o';
    name[2] = 'l';
    sprintf(buff,"%d",subdomain);
    name[3] = buff[0];
    name[4] = buff[1];
    name[5] = buff[2];

    stream = fopen(name,"a+");
    if (stream==NULL)
    {
      printf("%s\n", "cannot open file");
      return(1);
    }

    fprintf(stream,"%d %d %d %d\n",node0, node1, node2, node3);
    fclose(stream);
  }


  PRINTDEBUG(dom,1,(" add element %4d %4d %4d %4d\n",
                    Id[0], Id[1], Id[2], Id[3]));


  /*	if (InsertElementFromIDs(GRID_ON_LEVEL(currMG,0),4,Id) == NULL)
            return(1);*/

  /* write element in the mesh structure */
  currmesh->Element_corner_ids[subdomain][element] = (INT *) GetTmpMem(MGHEAP(currMG),4*sizeof(INT), GG3_MarkKey);

  currmesh->Element_corners[subdomain][element] = 4;
  for(i=0; i<4; i++)
    if(Id[i]<nb_boundary_points_subdom)
      currmesh->Element_corner_ids[subdomain][element][i] = oldId[Id[i]];
    else
      currmesh->Element_corner_ids[subdomain][element][i] = nb_boundary_points
                                                            + nb_inner_points
                                                            + Id[i]
                                                            - nb_boundary_points_subdom;
  /*	printf("%d %d %d %d\n", Id[0], Id[1], Id[2], Id[3]);
          printf("%d %d %d %d\n", currmesh->Element_corner_ids[subdomain][element][0],
                                                          currmesh->Element_corner_ids[subdomain][element][1],
                                                          currmesh->Element_corner_ids[subdomain][element][2],
                                                          currmesh->Element_corner_ids[subdomain][element][3]);*/
  element++;
  /*	theElement = InsertElementFromIDs(GRID_ON_LEVEL(currMG,0),4,Id,NULL);
          if (theElement==NULL)  return(1);*/
  /*SETSUBDOMAIN(theElement,id_from_somewhere);*/

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

static INT AddBoundaryElement (INT n, INT *nodelist)
{
  if (n != 3)
    return(1);

  PRINTDEBUG(dom,1,(" add triangle  %4d %4d %4d\n",
                    nodelist[0],nodelist[1],nodelist[2]));

    #ifdef _NETGEN
  AddSurfaceTriangle(nodelist[0],nodelist[1],nodelist[2]);
    #endif

  return(0);
}

/****************************************************************************/
/*D
   GenerateGrid3d - call the netgen grid generator

   SYNOPSIS:
   INT GenerateGrid3d (MULTIGRID *theMG, MESH *mesh, DOUBLE h, INT smooth, INT coeff);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  mesh - pointer to mesh
   .  h - mesh size
   .  smoothing parameter
   .  coeff - number of coefficientfunction for gridsize
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

INT GenerateGrid3d_OLD (MULTIGRID *theMG, MESH *mesh, DOUBLE h, INT smooth,
                        INT display, INT coeff)
{
  NODE *theNode;
  VERTEX *theVertex;
  INT sid,i;
  char rulefilename[128];
  DOUBLE **x;

  if (mesh == NULL)
    return(GM_OK);

  currMG = theMG;

  if(h<=0.0)
  {
    /* get Coefficientfunctions */
    if (BVP_SetCoeffFct(MG_BVP(theMG),-1,Coefficients))
      return (0);
    LOCAL_H = Coefficients[coeff];
  }

  if (GetDefaultValue(DEFAULTSFILENAME,"netgenrules",rulefilename))
    strcpy(rulefilename,"tetra.rls");
    #ifdef _NETGEN
  InitNetgen(rulefilename);
        #else
  PrintErrorMessage('E',"GenerateGrid3d","no netgen");
  return(1);
    #endif

  IFDEBUG(dom,1)
  x = (DOUBLE **) GetTmpMem(MGHEAP(theMG),mesh->nBndP*sizeof(DOUBLE *),MG_MARK_KEY(theMG));
  for (i=0, theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,0));
       theNode!=NULL; theNode=SUCCN(theNode))
  {
    PRINTDEBUG(dom,0,(" i %d  nid %ld nd %x %x %x \n",
                      i,ID(theNode),theNode,
                      V_BNDP(theVertex),mesh->theBndPs[i]));

    x[i++] = CVECT(MYVERTEX(theNode));
  }
  for (i=0, theNode=LASTNODE(GRID_ON_LEVEL(theMG,0));
       theNode!=NULL; theNode=PREDN(theNode))
  {
    PRINTDEBUG(dom,0,(" i %d  nid %ld nd %x %x %x \n",
                      i,ID(theNode),theNode,
                      V_BNDP(theVertex),mesh->theBndPs[i]));

    x[i++] = CVECT(MYVERTEX(theNode));
  }



  ENDDEBUG

  for (i=0, theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,0));
       theNode!=NULL; theNode=SUCCN(theNode))
  {
    theVertex = MYVERTEX(theNode);
    if (V_BNDP(theVertex) != mesh->theBndPs[i])
      return(1);
    if (AddBoundaryNode (i++,CVECT(theVertex)))
      return(1);
  }
  if (i != mesh->nBndP)
    return(1);

  for (sid=0; sid<=mesh->nSubDomains; sid++)
    for (i=0; i<mesh->nSides[sid]; i++)
    {

      DOUBLE mat[3][3];

      if (AddBoundaryElement (mesh->Side_corners[sid][i],
                              mesh->Side_corner_ids[sid][i]))
        return(1);

      PRINTDEBUG(dom,1,(" sid %ld (%4.2f,%4.2f,%4.2f) (%4.2f,%4.2f,%4.2f) (%4.2f,%4.2f,%4.2f)\n",
                        i,
                        x[mesh->Side_corner_ids[sid][i][0]][0],
                        x[mesh->Side_corner_ids[sid][i][0]][1],
                        x[mesh->Side_corner_ids[sid][i][0]][2],
                        x[mesh->Side_corner_ids[sid][i][1]][0],
                        x[mesh->Side_corner_ids[sid][i][1]][1],
                        x[mesh->Side_corner_ids[sid][i][1]][2],
                        x[mesh->Side_corner_ids[sid][i][2]][0],
                        x[mesh->Side_corner_ids[sid][i][2]][1],
                        x[mesh->Side_corner_ids[sid][i][2]][2]));

    }

    #ifdef _NETGEN
  if (StartNetgen(h,smooth,display)) return(1);
    #endif

  return(0);
}

/* only for debug */
static INT Write_SurfaceMesh(MESH *mesh, MULTIGRID *theMG)
{
  char buff[3], name[10], name1[14];
  INT i, j, k, sid, nnodes, nsides, nodelist[4];
  FILE *stream;
  DOUBLE local[3], global[3];
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
      fprintf(stream,"%lf %lf %lf\n",global[0], global[1], global[2]);
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
  INT i, j, k, sid, nnodes, nsides, nodelist[4], id[4];
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

  stream = fopen(name,"a+");
  if (stream==NULL)
  {
    printf("%s\n", "cannot open file");
    return(1);
  }

  fprintf(stream,"%d\n",(int)nb_boundary_points_subdom+(int)nInnP[subdomain]);
  for (i=0; i<mesh->nBndP; i++)
  {
    if(transfer[i])
    {
      if (BNDP_Global(mesh->theBndPs[i],global))
        return (1);
      fprintf(stream,"%lf %lf %lf\n",global[0], global[1], global[2]);
    }
  }
  for (i=0; i<nInnP[subdomain]; i++)
    fprintf(stream,"%lf %lf %lf\n",Position[subdomain][i][0], Position[subdomain][i][1], Position[subdomain][i][2]);

  /*	fprintf(stream,"%d\n",(int)mesh->nElements[subdomain]);

          for (i=subdomain; i<=subdomain; i++)
          {
                  for (j=0; j<mesh->nElements[subdomain];j++)
                  {
                          for(k=0;k<4;k++)
                                  if(mesh->Element_corner_ids[i][j][k]<nb_boundary_points_subdom)
                                          id[k] = newId[mesh->Element_corner_ids[i][j][k]];
                                  else
                                          id[k] = mesh->Element_corner_ids[i][j][k]
                                                          - nb_boundary_points
                                                          - nb_inner_points
   + nb_boundary_points_subdom;

                                  fprintf(stream,"%d %d %d %d\n",id[0], id[1], id[2], id[3]);
                  }
          }*/

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
  printf("%lf %lf\n", s, acos(s));
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
  INT i, j, k, sid, buflen;
  FILE *stream;
  char buffer[256];
  int nelements, npoints, d0, d1, d2, d3, id[4];
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

  fgets(buffer,buflen,stream);
  fgets(buffer,buflen,stream);
  sscanf(buffer,"%d",&nelements);

  mesh->nElements[subdomain] = nelements;
  mesh->Element_corners[subdomain] = (INT *) GetTmpMem(MGHEAP(currMG),(nelements+1)*sizeof(INT), MarkKey);
  mesh->Element_corner_ids[subdomain] = (INT **) GetTmpMem(MGHEAP(currMG),(nelements+1)*sizeof(INT*), MarkKey);

  for (i=subdomain; i<=subdomain; i++)
  {
    for (j=0; j<mesh->nElements[subdomain]; j++)
    {
      fgets(buffer,buflen,stream);
      sscanf(buffer,"%d %d %d %d",&d0, &d1, &d2, &d3);
      id[0] = d1; id[1] = d0; id[2] = d2; id[3] = d3;

      mesh->Element_corners[subdomain][j] = 4;
      mesh->Element_corner_ids[subdomain][j] = (INT *) GetTmpMem(MGHEAP(currMG),4*sizeof(INT), MarkKey);
      if(mesh->Element_corner_ids[subdomain][j]==NULL)
        return(0);

      for(k=0; k<4; k++)
        if(id[k]<nb_boundary_points_subdom)
          mesh->Element_corner_ids[subdomain][j][k] = oldId[id[k]];
        else
          mesh->Element_corner_ids[subdomain][j][k] = nb_boundary_points
                                                      + nb_inner_points
                                                      + id[k]
                                                      - nb_boundary_points_subdom;

      /*				printf("%d %d %d %d\n", id[0], id[1], id[2], id[3]);*/
      /*				printf("%d %d %d %d %d\n", j, mesh->Element_corner_ids[subdomain][j][0],
                                                                                      mesh->Element_corner_ids[subdomain][j][1],
                                                                                      mesh->Element_corner_ids[subdomain][j][2],
                                                                                      mesh->Element_corner_ids[subdomain][j][3]);*/
    }
  }

  fgets(buffer,buflen,stream);
  sscanf(buffer,"%d",&npoints);

  nInnP[subdomain] = npoints - nb_boundary_points_subdom;

  Position[subdomain] = (DOUBLE **) GetTmpMem(MGHEAP(currMG),(nInnP[subdomain]+1)*sizeof(DOUBLE*), MarkKey);
  if(Position[subdomain]==NULL)
    return(0);
  for(i=0; i<nInnP[subdomain]; i++)
  {
    Position[subdomain][i] = (DOUBLE *) GetTmpMem(MGHEAP(currMG),3*sizeof(DOUBLE), MarkKey);
    if(Position[subdomain][i]==NULL)
      return(0);
  }
  for (i=0; i<nb_boundary_points_subdom; i++)
  {
    fgets(buffer,buflen,stream);
    sscanf(buffer,"%lf %lf %lf",&g0, &g1, &g2);
    /*		Position[subdomain][i][0] = g0;
                    Position[subdomain][i][1] = g1;
                    Position[subdomain][i][2] = g2;*/
  }

  for (i=0; i<nInnP[subdomain]; i++)
  {
    fgets(buffer,buflen,stream);
    sscanf(buffer,"%lf %lf %lf",&g0, &g1, &g2);
    Position[subdomain][i][0] = g0;
    Position[subdomain][i][1] = g1;
    Position[subdomain][i][2] = g2;
  }

  fclose(stream);

  return(0);
}

static INT MAX_ELE = 75;

static INT Search_Tet_Neighbours(MULTIGRID *theMG, MESH *mesh, INT MarkKey)
{
  INT sid, np, i, j, k, l, a1, a2, a3, b1, b2, b3, npoints, corner_of_tet;
  INT **pointlist;

  nbElement = (INT ***) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT**), MarkKey);
  if(nbElement==NULL)
    return(0);
  for(sid=1; sid<=mesh->nSubDomains; sid++)
  {
    nbElement[sid] = (INT **) GetTmpMem(MGHEAP(currMG),(mesh->nElements[sid]+1)*sizeof(INT*), MarkKey);
    if(nbElement[sid]==NULL)
      return(0);
    for(i=0; i<mesh->nElements[sid]; i++)
    {
      nbElement[sid][i] = (INT *) GetTmpMem(MGHEAP(currMG),4*sizeof(INT), MarkKey);
      if(nbElement[sid][i]==NULL)
        return(0);
    }
  }

  npoints =  nb_inner_points + nb_boundary_points;
  pointlist = (INT **) GetTmpMem(MGHEAP(currMG),(npoints+1)*sizeof(INT*), MarkKey);
  for(i=0; i<npoints+1; i++)
    pointlist[i] = (INT *) GetTmpMem(MGHEAP(currMG),MAX_ELE*sizeof(INT), MarkKey);


  for(sid=1; sid<=mesh->nSubDomains; sid++)
  {
    for(i=0; i<mesh->nElements[sid]; i++)
      for(j=0; j<4; j++)
        nbElement[sid][i][j] = -1;
    for(i=0; i<npoints; i++)
      for(j=0; j<MAX_ELE; j++)
        pointlist[i][j] = 0;

    for(i=0; i<mesh->nElements[sid]; i++)
      for(j=0; j<4; j++)
      {
        corner_of_tet = mesh->Element_corner_ids[sid][i][j];
        pointlist[corner_of_tet][++pointlist[corner_of_tet][0]] = i;
      }

    for(np=0; np<npoints; np++)
      if(pointlist[np][0]>1)
        for(i=1; i<=pointlist[np][0]; i++)
          for(j=1; j<=pointlist[np][0]; j++)
            if(i!=j)
              for(k=0; k<4; k++)
                for(l=0; l<4; l++)
                {
                  a1 = mesh->Element_corner_ids[sid][pointlist[np][i]][(k+1)%4];
                  a2 = mesh->Element_corner_ids[sid][pointlist[np][i]][(k+2)%4];
                  a3 = mesh->Element_corner_ids[sid][pointlist[np][i]][(k+3)%4];
                  b1 = mesh->Element_corner_ids[sid][pointlist[np][j]][(l+3)%4];
                  b2 = mesh->Element_corner_ids[sid][pointlist[np][j]][(l+2)%4];
                  b3 = mesh->Element_corner_ids[sid][pointlist[np][j]][(l+1)%4];
                  if( ((a1==b1) && (a2==b2) && (a3==b3))
                      || ((a1==b2) && (a2==b3) && (a3==b1))
                      || ((a1==b3) && (a2==b1) && (a3==b2))
                      || ((a1==b1) && (a2==b3) && (a3==b2))
                      || ((a1==b2) && (a2==b1) && (a3==b3))
                      || ((a1==b3) && (a2==b2) && (a3==b1)) )
                    nbElement[sid][pointlist[np][i]][(k+1)%4] = pointlist[np][j];
                }

  }
  return(0);
}

INT GenerateGrid3d (MULTIGRID *theMG, MESH *mesh, DOUBLE h, INT smooth,
                    INT display, INT coeff)
{
  NODE *theNode;
  VERTEX *theVertex;
  INT sid,i,j, k, nodelist[4], Id[4], bnds_flag[4], k1, k2;
  char rulefilename[128];
  DOUBLE global[3];
  char buff[3], name[6];
  FILE *stream;
  ELEMENT *theElement;

  GG3_MarkKey = MG_MARK_KEY(theMG);

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
  Position = (DOUBLE ***) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(DOUBLE**), GG3_MarkKey);

  mesh->nElements = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT), GG3_MarkKey);
  mesh->Element_corners = (INT **) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT*), GG3_MarkKey);
  mesh->Element_corner_ids = (INT ***) GetTmpMem(MGHEAP(theMG),(mesh->nSubDomains+1)*sizeof(INT**), GG3_MarkKey);

  transfer = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nBndP+1)*sizeof(INT), GG3_MarkKey);
  newId = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nBndP+1)*sizeof(INT), GG3_MarkKey);
  oldId = (INT *) GetTmpMem(MGHEAP(theMG),(mesh->nBndP+1)*sizeof(INT), GG3_MarkKey);

  nb_boundary_points = theMG->nodeIdCounter;
  nb_inner_points = 0;

  /* triangulate every subdomain */
  for (sid=1; sid<=mesh->nSubDomains; sid++)
  {
    subdomain = sid;

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
        if (AddBoundaryNode (newId[i],global))
          return(1);
      }
    }
    /* transfer surface-triangles to netgen */
    for (i=0; i<mesh->nSides[sid]; i++)
    {
      for(j=0; j<mesh->Side_corners[sid][i]; j++)
        nodelist[j] = newId[mesh->Side_corner_ids[sid][i][j]];
      if (AddBoundaryElement (mesh->Side_corners[sid][i],nodelist))
        return(1);
    }

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

      if(GG3_DEBUG)
        Write_SurfaceMesh(mesh, theMG);
                        #ifdef _NETGEN
      if (StartNetgen(h,smooth,display)) return(1);
                        #endif
      if(GG3_DEBUG)
        Write_VolumeMesh(mesh, theMG);
    }
    else
    {
      fclose(stream);
      printf("%s %d\n", "Read Subdomain ", subdomain);
      Read_VolumeMesh(mesh, theMG, GG3_MarkKey);
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
  for(i=0; i<mesh->nInnP; i++)
    mesh->Position[i] = (DOUBLE *) GetTmpMem(MGHEAP(currMG),3*sizeof(DOUBLE), GG3_MarkKey);
  k = 0;
  for (sid=1; sid<=mesh->nSubDomains; sid++)
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
  Search_Tet_Neighbours(theMG, mesh, GG3_MarkKey);

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
  for (sid=1; sid<=mesh->nSubDomains; sid++)
  {
    subdomain = sid;
    k1 = k2+1;
    for(i=0; i<mesh->nElements[sid]; i++)
    {
      for(j=0; j<4; j++)
      {
        Id[j] = mesh->Element_corner_ids[sid][i][j];
        if(nbElement[sid][i][j]==-1)
          bnds_flag[j] = 1;
        else
          bnds_flag[j] = 0;
      }
      /*			printf("%d %d %d %d %d %d %d %d\n", mesh->Element_corner_ids[sid][i][0],
                                                                  mesh->Element_corner_ids[sid][i][1],
                                                                  mesh->Element_corner_ids[sid][i][2],
                                                                  mesh->Element_corner_ids[sid][i][3],
                                                                  nbElement[sid][i][0],
                                                                  nbElement[sid][i][1],
                                                                  nbElement[sid][i][2],
                                                                  nbElement[sid][i][3]);*/
      /*			if (InsertElementFromIDs_New(GRID_ON_LEVEL(theMG,0),mesh->Element_corners[sid][i],Id, bnds_flag, sid) == NULL)
                                      return(1);*/
      theElement = InsertElementFromIDs(GRID_ON_LEVEL(currMG,0),4,Id,bnds_flag);
      if (theElement==NULL)
        return(1);
      SETSUBDOMAIN(theElement,subdomain);

    }
    k2 = k2 + mesh->nElements[sid];
    printf("%s %d %d %d\n", "subdomain ", sid, k1, k2);
  }

  return(0);
}

int Get_h (double *in, double *out)
{
  (*LOCAL_H)(in, out);
  return(0);
}
