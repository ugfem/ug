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
#include "uginterface.h"
#include "general.h"

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

#define DEBUG

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

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

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
  COORD_VECTOR xc;

    #ifdef DEBUG
  UserWriteF(" add inner node %4d %6.3lf %6.3lf %6.3lf\n",nodeid,
             x,y,z);
    #endif

  xc[0] = x;
  xc[1] = y;
  xc[2] = z;

  if (InsertInnerNode(currMG,xc)!=GM_OK)
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

    #ifdef DEBUG
  UserWriteF(" add element %4d %4d %4d %4d\n",node0,node1,node2,node3);
    #endif

  Id[0] = node0;
  Id[1] = node1;
  Id[2] = node2;
  Id[3] = node3;

  if (InsertElementFromIDs(currMG,4,Id)!=GM_OK)
    return(1);

  return(0);
}

/****************************************************************************/
/*D
   AddBoundaryNode - append a node in the boundary node list of gg3d

   SYNOPSIS:
   static INT AddBoundaryNode (INT nodeid, COORD *global);

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

static INT AddBoundaryNode (INT nodeid, COORD *global)
{
    #ifdef DEBUG
  UserWriteF(" add node   %4d    %6.3f %6.3f %6.3f\n",
             nodeid,global[0],global[1],global[2]);
    #endif

    #ifdef NETGENT
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

  if (left > 0)
  {
        #ifdef NETGENT
    AddSurfaceTriangle(nodelist[0],nodelist[1],nodelist[2]);
        #endif
    triangleid++;
        #ifdef DEBUG
    UserWriteF(" add triangle  %4d %4d %4d\n",
               nodelist[0],nodelist[1],nodelist[2]);
        #endif
  }
  if (right > 0)
  {
        #ifdef NETGENT
    AddSurfaceTriangle(nodelist[0],nodelist[2],nodelist[1]);
        #endif
    triangleid++;
        #ifdef DEBUG
    UserWriteF(" add triangle  %4d %4d %4d\n",
               nodelist[0],nodelist[2],nodelist[1]);
        #endif
  }

  return(0);
}

/****************************************************************************/
/*D
   AddBoundaryElements - decompose a strip into triangles

   SYNOPSIS:
   static INT AddBoundaryElements (INT n, INT m,
   INT c0, INT c1, INT c2, INT c3,
   INT s0, INT s1, INT s2, INT s3);

   PARAMETERS:
   .  n,m - stripe with n+1 nodes on the bottom and m+1 nodes on the top
   .  c0,c1,c2,c3 - corner node ids
   .  s0,s1,s2,s3 - side node ids

   DESCRIPTION:
   This function splits a stripe into triangles an calls 'AddBoundaryElements'.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

static INT AddBoundaryElements (INT n, INT m,
                                INT c0, INT c1, INT c2, INT c3,
                                INT s0, INT s1, INT s2, INT s3)
{
  INT nodelist[3];

  if (m < n)
  {
    if (n == 1)
    {
      nodelist[0] = c0;
      nodelist[1] = c2;
      nodelist[2] = c1;
      AddBoundaryElement(3,nodelist);
    }
    else
    {
      nodelist[0] = c0;
      nodelist[1] = c2;
      nodelist[2] = s0;
      AddBoundaryElement(3,nodelist);
      c0 = s0;
      if (s0<s1) s0++;else s0--;
      AddBoundaryElements (n-1,m,c0,c1,c2,c3,s0,s1,s2,s3);
    }
  }
  else
  {
    if (m == 1)
    {
      nodelist[0] = c0;
      nodelist[1] = c2;
      nodelist[2] = c1;
      AddBoundaryElement(3,nodelist);
      nodelist[0] = c1;
      nodelist[1] = c2;
      nodelist[2] = c3;
      AddBoundaryElement(3,nodelist);
    }
    else
    {
      nodelist[0] = c2;
      nodelist[1] = s2;
      nodelist[2] = c0;
      AddBoundaryElement(3,nodelist);
      c2 = s2;
      if (s2<s3) s2++;else s2--;
      AddBoundaryElements (n,m-1,c0,c1,c2,c3,s0,s1,s2,s3);
    }
  }

  return(0);
}

/****************************************************************************/
/*D
   TriangulatePatch - decompose a patch into stripes

   SYNOPSIS:
   static INT TriangulatePatch (DOUBLE h, PATCH *thePatch,
   INT npc, INT *cornerid, COORD local[CORNERS_OF_BND_SEG][DIM-1],
   INT sideid[CORNERS_OF_BND_SEG][2], INT *siden);

   PARAMETERS:
   .  h - maximal length of an edge of the boundary triangles
   .  thePatch - poiter to a patch
   .  npc - number of corner nodes of the patch
   .  cornerid - ids of the corner nodes
   .  local - local coordinates of the corners
   .  sideid - ids of the nodes on the edges of the patch
   .  siden - number of nodes on the edges of the patch

   DESCRIPTION:
   This function splits the patch into stripes and calls
   'AddBoundaryElements' to decompose the stripes into triangles.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

static INT TriangulatePatch (DOUBLE h, PATCH *thePatch, INT npc, INT *cornerid,
                             COORD local[CORNERS_OF_BND_SEG][DIM-1],
                             INT sideid[CORNERS_OF_BND_SEG][2], INT *siden)
{
  INT i,j,k,n;
  COORD lvect[DIM-1],gvect[DIM],gvect1[DIM];
  COORD next_local[CORNERS_OF_BND_SEG][DIM-1];
  DOUBLE lambda,dist,step;
  INT next_cornerid[CORNERS_OF_BND_SEG];
  INT next_sideid[CORNERS_OF_BND_SEG][2];
  INT next_siden[CORNERS_OF_BND_SEG];

  if ((siden[npc-1] > 1) && (siden[1] > 1))
  {
    for (k=2; k<npc; k++)
      for (i=0; i<DIM-1; i++)
        next_local[k][i] = local[k][i];
    next_siden[0] = n;
    next_siden[1] = siden[1] - 1;
    next_siden[2] = siden[2];
    next_siden[npc-1] = siden[npc-1] - 1;
    next_sideid[2][0] = sideid[2][0];
    next_sideid[2][1] = sideid[2][1];
    for (k=2; k<npc; k++)
      next_cornerid[k] = cornerid[k];
    next_cornerid[0] = sideid[npc-1][1];
    next_sideid[npc-1][0] = sideid[npc-1][0];
    if (sideid[npc-1][0] < sideid[npc-1][1])
      next_sideid[npc-1][1] = sideid[npc-1][1] - 1;
    else
      next_sideid[npc-1][1] = sideid[npc-1][1] + 1;
    lambda = (siden[npc-1] - 1.0) / siden[npc-1];
    V2_LINCOMB(lambda,local[0],(1.0-lambda),local[npc-1],
               next_local[0]);
    next_cornerid[1] = sideid[1][0];
    next_sideid[1][1] = sideid[1][1];
    if (sideid[1][0] < sideid[1][1])
      next_sideid[1][0] = sideid[1][0] + 1;
    else
      next_sideid[1][0] = sideid[1][0] - 1;
    lambda = (siden[1] - 1.0) / siden[1];
    V2_LINCOMB(lambda,local[1],(1.0-lambda),local[2],next_local[1]);
    Patch_local2global(thePatch,next_local[0],gvect);
    Patch_local2global(thePatch,next_local[1],gvect1);
    V_DIM_EUKLIDNORM_OF_DIFF(gvect,gvect1,dist);
    next_siden[0] = siden[0];
    if (((INT)(dist/h)) < siden[0])
      next_siden[0] -= 1;
    else if (((INT)(dist/h)) > siden[0])
      next_siden[0] += 1;
    next_sideid[0][0] = nodeid;
    next_sideid[0][1] = nodeid + next_siden[0] - 2;
    step = 1.0 / next_siden[0];
    for (i=1; i<next_siden[0]; i++)
    {
      lambda = i * step;
      V2_LINCOMB(lambda,next_local[1],(1.0-lambda),next_local[0],lvect);
      Patch_local2global(thePatch,lvect,gvect);
      if (AddBoundaryNode(nodeid,gvect))
      {
        Release(MGHEAP(currMG),FROM_TOP);
        return(1);
      }
      if (InsertBoundaryNodeFromPatch(currMG,thePatch,lvect)!=GM_OK)
      {
        Release(MGHEAP(currMG),FROM_TOP);
        return(1);
      }
      nodeid++;
    }

    AddBoundaryElements(siden[0],next_siden[0],
                        cornerid[0],cornerid[1],
                        next_cornerid[0],next_cornerid[1],
                        sideid[0][0],sideid[0][1],
                        next_sideid[0][0],next_sideid[0][1]);

    return(TriangulatePatch(h,thePatch,npc,next_cornerid,
                            next_local,next_sideid,next_siden));
  }
  else if ((siden[npc-1] = 1) && (siden[1] = 1))
  {
    if (npc == 3)
      return(AddBoundaryElements(siden[0],0,
                                 cornerid[0],cornerid[1],cornerid[2],0,
                                 sideid[0][0],sideid[0][1],0,0));
    else
      return(AddBoundaryElements(siden[0],siden[2],
                                 cornerid[0],cornerid[1],
                                 cornerid[3],cornerid[2],
                                 sideid[0][0],sideid[0][1],
                                 sideid[2][1],sideid[2][0]));
  }
  else if (((siden[npc-1] > 1) && (siden[1] = 1)) ||
           ((siden[npc-1] = 1) && (siden[1] > 1))   )
  {
    for (i=0; i<DIM-1; i++)
      next_local[0][i] = local[npc-1][i];
    next_siden[0] = siden[npc-1];
    next_sideid[0][0] = sideid[npc-1][0];
    next_sideid[0][1] = sideid[npc-1][1];
    for (k=1; k<npc; k++)
    {
      for (i=0; i<DIM-1; i++)
        next_local[k][i] = local[k+1][i];
      next_siden[k] = siden[k+1];
      next_sideid[k][0] = sideid[k+1][0];
      next_sideid[k][1] = sideid[k+1][1];
    }

    return(TriangulatePatch(h,thePatch,npc,next_cornerid,
                            next_local,next_sideid,next_siden));
  }


  UserWriteF("TriangulatePatch: not implemented\n");

  return(1);
}

/****************************************************************************/
/*D
   GenerateBoundaryNodes3d - computing of a boundary triangulation

   SYNOPSIS:
   INT GenerateBoundaryNodes3d (MULTIGRID *theMG, DOUBLE h);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  h - maximal length of an edge of the boundary triangles

   DESCRIPTION:
   This function adds node on the edges of all patches with a
   distance of approximately 'h' and calls 'TriangulatePatch'
   to decompose every patch into triangles.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

#define INDEX_IN_LIST(from,to,nc)  (from < to ? from*nc+to : to*nc+from)

INT GenerateBoundaryNodes3d (MULTIGRID *theMG, DOUBLE h)
{
  VERTEX **pv;
  INT i,j,k,n,from,to,nc,npc;
  COORD lvect[DIM-1],gvect[DIM];
  COORD local[CORNERS_OF_BND_SEG][DIM-1],*from_local,*to_local;
  DOUBLE lambda,dist,step;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;
  INT *vlist;
  INT cornerid[CORNERS_OF_BND_SEG];
  INT sideid[CORNERS_OF_BND_SEG][2],siden[CORNERS_OF_BND_SEG];
  char rulefilename[128];

  currMG = theMG;
  theBVP = MG_BVP(theMG);
  if (theBVP==NULL)
  {
    PrintErrorMessage('E',"GenerateBoundaryNodes3d","BVP not found");
    return(1);
  }

  if (BVP_GetBVPDesc(theBVP,&theBVPDesc))
  {
    PrintErrorMessage('E',"CreateMultiGrid","BVP not evaluated");
    return(1);
  }

  pv = theMG->corners;
  nc = BVPD_NCORNERS(theBVPDesc);
  h_global = h;

  if (GetDefaultValue(DEFAULTSFILENAME,"netgenrules",rulefilename))
    strcpy(rulefilename,"tetra.rls");
    #ifdef NETGENT
  InitNetgen(rulefilename);
    #endif
  for (i=0; i<nc; i++)
    AddBoundaryNode(i,CVECT(pv[i]));
  nodeid = nc;
  triangleid = 0;
  vlist = (INT*) GetMem(MGHEAP(theMG),nc*nc*sizeof(INT),FROM_TOP);
  for (i=0; i<nc; i++)
    for (j=0; j<nc; j++)
      vlist[INDEX_IN_LIST(i,j,nc)] = 0;

  for (thePatch=BVP_GetFirstPatch(theBVP); thePatch!=NULL;
       thePatch=BVP_GetNextPatch(theBVP,thePatch))
  {
    Patch_GetPatchDesc(thePatch,&thePatchDesc);
    left = PATCH_LEFT(thePatchDesc);
    right = PATCH_RIGHT(thePatchDesc);
        #ifdef DEBUG
    UserWriteF(" PID %d %d %d\n",PATCH_ID(thePatchDesc),
               PATCH_LEFT(thePatchDesc),PATCH_RIGHT(thePatchDesc));
        #endif
    npc = PATCH_N(thePatchDesc);
    if (npc > 4)
    {
      Release(MGHEAP(theMG),FROM_TOP);
      return(1);
    }
    from = PATCH_CID(thePatchDesc,npc-1);
    from_local =  PATCH_LCVECT(thePatchDesc,npc-1);
    for( k=0; k<npc; k++ )
    {
      to = PATCH_CID(thePatchDesc,k);
      to_local = PATCH_LCVECT(thePatchDesc,k);
      for (i=0; i<DIM-1; i++)
        local[k][i] = from_local[i];
      cornerid[k] = from;
      V_DIM_EUKLIDNORM_OF_DIFF(CVECT(pv[from]),CVECT(pv[to]),
                               dist);
      n = MAX(1,dist/h);
      siden[k] = n;
      sideid[k][0] = vlist[INDEX_IN_LIST(from,to,nc)];
      if (sideid[k][0] > 0)
      {
        if (from < to)
          sideid[k][1] = sideid[k][0] + n - 2;
        else
        {
          sideid[k][1] = sideid[k][0];
          sideid[k][0] = sideid[k][1] + n - 2;
        }
                #ifdef DEBUG
        UserWriteF(" VID %d %d sideid %d %d\n",from,to,
                   sideid[k][0],sideid[k][1]);
                #endif
        from = to;
        from_local =  to_local;
        continue;
      }
      vlist[INDEX_IN_LIST(from,to,nc)] = nodeid;
      if (from < to)
      {
        sideid[k][0] = nodeid;
        sideid[k][1] = nodeid + n - 2;
      }
      else
      {
        sideid[k][1] = nodeid;
        sideid[k][0] = nodeid + n - 2;
      }
                #ifdef DEBUG
      UserWriteF(" VID %d %d sideid %d %d\n",from,to,
                 sideid[k][0],sideid[k][1]);
            #endif
      step = 1.0 / n;
      for (i=1; i<n; i++)
      {
        if (from < to)
          lambda = i * step;
        else
          lambda = (n-i) * step;
        V2_LINCOMB(lambda,to_local,(1.0-lambda),from_local,lvect);
        Patch_local2global(thePatch,lvect,gvect);
        if (AddBoundaryNode(nodeid,gvect))
        {
          Release(MGHEAP(theMG),FROM_TOP);
          return(1);
        }
        if (InsertBoundaryNodeFromPatch(theMG,thePatch,lvect)!=GM_OK)
        {
          Release(MGHEAP(theMG),FROM_TOP);
          return(1);
        }
        nodeid++;
      }
      from = to;
      from_local =  to_local;
    }
    if (TriangulatePatch(h,thePatch,npc,cornerid,
                         local,sideid,siden))
    {
      Release(MGHEAP(theMG),FROM_TOP);
      return(1);
    }
  }

  Release(MGHEAP(theMG),FROM_TOP);

  UserWriteF("bnodes: number of nodes     %d\n", nodeid);
  UserWriteF("        number of triangles %d\n", triangleid);

  return(0);
}

/****************************************************************************/
/*D
   GenerateGrid3d - call the netgen grid generator

   SYNOPSIS:
   INT GenerateGrid3d();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function is an interface function for the grid generator.
   It calls the netgen grid genarator to compute an inner
   decomposition into tetrahedrons. It requires that 'GenerateBoundaryNodes'
   is called before.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

INT GenerateGrid3d(INT smooth)
{
    #ifdef NETGENT
  if (StartNetgen(h_global,smooth)) return(1);
    #endif

  return(0);
}
