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

static CoeffProcPtr Coefficients[8];
static CoeffProcPtr LOCAL_H;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/*****************************************************************************/

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

  if (InsertInnerNode(GRID_ON_LEVEL(currMG,0),xc) == NULL)
    return(-1);

  return(nodeid++);
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
  INT Id[4];
  ELEMENT *theElement;

  PRINTDEBUG(dom,1,(" add element %4d %4d %4d %4d\n",
                    node0,node1,node2,node3));

  Id[0] = node0;
  Id[1] = node1;
  Id[2] = node2;
  Id[3] = node3;

  theElement = InsertElementFromIDs(GRID_ON_LEVEL(currMG,0),4,Id,NULL);
  if (theElement==NULL) return(1);
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

INT GenerateGrid3d (MULTIGRID *theMG, MESH *mesh, DOUBLE h, INT smooth,
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

int Get_h (double *in, double *out)
{
  (*LOCAL_H)(in, out);
  return(0);
}
