// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      smooth.c	                                                        */
/*                                                                          */
/* Purpose:   smoothing a grid                                              */
/*                                                                          */
/* Author:    Bernhard Huurdeman                                            */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   17.07.96 begin, ug version 3.3                                */
/*                                                                          */
/* Remarks:                                                                                                     */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "compiler.h"
#include "gm.h"
#include "evm.h"
#include "misc.h"
#include "shapes.h"
#include "ugm.h"
#include "ugenv.h"
#include "algebra.h"
#include "refine.h"
#include "devices.h"

#define SMALL_LOCAL    2.E-4

#define LOCAL_EQUAL(A,V)     (ABS((A)-(V))< (SMALL_LOCAL))
#define IS_VALUE2(A,V)     (ABS((A)-(V))< (SMALL_C))
#define IS_0_OR_1(V)   (LOCAL_EQUAL(V,0) || LOCAL_EQUAL(V,1))
#define V2_LOCAL_EQUAL(A,B) ((ABS((A)[0]-(B)[0])<SMALL_LOCAL)&&(ABS((A)[1]-(B)[1])<SMALL_LOCAL))

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*
   OneSideMove - calculate displacement of the center point along one local coordinate

   SYNOPSIS
   static DOUBLE OneSideMoveCP(DOUBLE *CenterPVertex, DOUBLE *sideMid,
                         DOUBLE *CenterPointNeEL)

   PARAMETERS:
   . CenterPVertex    - global coordinates of center vertex to move
   . sideMid          - global coordinates of side mid vertex on the side located on the father
                   and neighbour element
   . CenterPointNeEl  - global coordinates of center vertex of the neighbour element

   DESCRIPTION:
   move center point of father element according to the size of the neighbour element

   CP: center point of father element; SM: side mid of side between father and neighbour element;
   CPNE: center point of neighbour element

   father element        neighbour element
   ------------------------------------------------
 |      CP      |SM|            CPNE            |
   ------------------------------------------------
         <--x1---><-----x2------->

   father element: NCP new center point  OCP: old center point
   ------------------------------------
 |        NCP     OCP               |
   ------------------------------------
                   <----x1---------->
   <--y1----><--------y2-------------->

   calculate:
   y1 = 2*x1/(1+(x2/x1)^1/2)  , (2*x1 = y1 + y2)
   moving distance in direction SM-CP in local coordinates:
   LocalMove = 1/2 *(y1/x1) - 1/2

   RETURN VALUE:
   DOUBLE
   .n  LocalMove  - local coordinate of displacement (moving direction is from CenterPVertex to sideMid)
 */
/***************************************************************/

static DOUBLE OneSideMoveCP(DOUBLE *CenterPVertex, DOUBLE *sideMid,
                            DOUBLE *CenterPointNeEL)
{
  DOUBLE x1,x2,yc1;
  DOUBLE LocalMove;

  V_DIM_EUKLIDNORM_OF_DIFF(sideMid,CenterPVertex,x1);
  V_DIM_EUKLIDNORM_OF_DIFF(CenterPointNeEL,sideMid,x2);
  assert(x1!=0 && x2!=0);
  yc1 = 2*x1/(1+sqrt(x2/x1));
  LocalMove = 0.5*(yc1/x1) - 0.5;
  return(LocalMove);
}

/****************************************************************************/
/*
   NewPosCenterNodeCurved - modify displacement of center node according to a curved boundary

   SYNOPSIS
   static INT NewPosCenterNodeCurved(ELEMENT *theElement, DOUBLE *LocCoord)

   PARAMETERS
   .  theElement     - father element of the center node
   .  LocCoord       - new local coordinates of center node of theElement

   DESCRIPTION
   For boundary father elements on a curved boundary it exists one midnode with a 'MOVED' vertex on the
   boundary. In this functions this is taken into account by modifying the displacement of the center node
   calculated before with the function NewPosCenterNode.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

static INT NewPosCenterNodeCurved(ELEMENT *theElement, DOUBLE *LocCoord)
{
  INT i,j,n,nmoved,found;
  ELEMENT *sonElement;
  ELEMENT *SonList[MAX_SONS];
  NODE *theNode[MAX_SIDES_OF_ELEM];
  VERTEX *bndVertex[MAX_SIDES_OF_ELEM],*theVertex;
#ifdef __TWODIM__
  ELEMENT *fatherElement;
  NODE *cornerNode[MAX_CORNERS_OF_SIDE];
  DOUBLE lcorn0[DIM],lcorn1[DIM],*lmid,totalLocal,*theCorners[MAX_CORNERS_OF_ELEM];
  INT ncorn,edge,co0,co1;
#endif

#ifdef __THREEDIM__
  PrintErrorMessage('E',"NewPosCenterNodeCurved","3D not implemented yet");
  return(1);
#endif
  nmoved = 0;
  /*  printf("fatherElement %ld, centerNode %ld, x , y %f %f \n",ID(theElement),ID(centerNode),XC(MYVERTEX(centerNode)),YC(MYVERTEX(centerNode))); */
  /* find boundary sides with moved vertices */
  GetAllSons(theElement, SonList);
  for (i=0; i<NSONS(theElement); i++)
  {
    sonElement = SonList[i];
    for (j=0; j<CORNERS_OF_ELEM(sonElement); j++)
    {
      /* take only midnodes */
      if (NTYPE(CORNER(sonElement,j))!=MID_NODE) continue;
      theVertex = MYVERTEX(CORNER(sonElement,j));
      if (MOVED(theVertex))
      {
        /* already found ? */
        found=FALSE;
        for (n=0; n<nmoved; n++)
          if (theVertex==bndVertex[n]) found=TRUE;
        if (found==FALSE)
        {
          theNode[nmoved] = CORNER(sonElement,j);
          bndVertex[nmoved] = MYVERTEX(theNode[nmoved]);
          nmoved++;
        }
        break;
      }
    }
  }

#ifdef __TWODIM__
  assert(nmoved<=2);

  for (n=0; n<nmoved; n++)
  {
    lmid = LCVECT(bndVertex[n]);
    /* find the two corner nodes of the Element on the boundary */
    edge = ONEDGE(bndVertex[n]);
    fatherElement = VFATHER(bndVertex[n]);
    co0 = CORNER_OF_EDGE(fatherElement,edge,0);
    co1 = CORNER_OF_EDGE(fatherElement,edge,1);
    cornerNode[0] = SONNODE(CORNER(fatherElement,co0));
    cornerNode[1] = SONNODE(CORNER(fatherElement,co1));

    /*      nbn = 0;
          for (theLink=START(theNode[n]);theLink!=0;theLink=NEXT(theLink))
            {
              if (CORNERTYPE(NBNODE(theLink)))
                {
                  cornerNode[nbn] = NBNODE(theLink);
                  nbn++;
                }
            }
          assert(nbn==2); */


    CORNER_COORDINATES(theElement,ncorn,theCorners);
    UG_GlobalToLocal(ncorn,(const DOUBLE **)theCorners,
                     CVECT(MYVERTEX(cornerNode[0])),lcorn0);
    UG_GlobalToLocal(ncorn,(const DOUBLE **)theCorners,
                     CVECT(MYVERTEX(cornerNode[1])),lcorn1);
    found = FALSE;
    for (i=0; i<DIM; i++)
      if (LOCAL_EQUAL(lcorn0[i],lcorn1[i]))
      {
        if (LOCAL_EQUAL(lcorn0[i],0))
        {
          totalLocal = 1-lmid[i];
          LocCoord[i] = totalLocal*LocCoord[i] + lmid[i];
          found = TRUE;
          break;
        }
        else if (LOCAL_EQUAL(lcorn0[i],1))
        {
          totalLocal = lmid[i];
          LocCoord[i] = totalLocal*LocCoord[i];
          found = TRUE;
          break;
        }
      }
  }
  if (found==FALSE)
  {
    printf("NewPosCenterNodeCurved lcorn0: %f %f, lcorn1: %f %f \n",lcorn0[0],lcorn0[1],lcorn1[0],lcorn1[1]);
    printf("center node nacher: xi=%f  eta=%f \n",LocCoord[0],LocCoord[1]);
  }
#endif
  return(0);
}

/****************************************************************************/
/*
   MovedNode - check for 'MOVED' boundary vertices on a father element

   SYNOPSIS
   static INT MovedNode (ELEMENT *theElement)

   PARAMETERS
   .  theElement     - father element to check

   DESCRIPTION
   Search in all son elements of theElement whether there are 'MOVED' vertices existing or not.

   RETURN VALUE
   INT
   .n   0: no 'MOVED' vertices found
   .n   1: at least one 'MOVED' vertex found
 */
/****************************************************************************/

static INT MovedNode (ELEMENT *theElement)
{
  INT i,j;
  ELEMENT *sonElement;
  ELEMENT *SonList[MAX_SONS];
  NODE *theNode;
  VERTEX *theVertex;
  /*  printf("MovedNode? fatherElement %ld \n",ID(theElement));*/
  GetAllSons(theElement, SonList);
  for (i=0; i<NSONS(theElement); i++)
  {
    sonElement = SonList[i];
    for (j=0; j<CORNERS_OF_ELEM(sonElement); j++)
    {
      theNode = CORNER(sonElement,j);
      /* check only midnodes */
      if (NTYPE(theNode)!=MID_NODE) continue;
      theVertex = MYVERTEX(theNode);
      /* printf("Moved Node? %ld, xc yc %f %f\n",ID(theNode),XC(theVertex),YC(theVertex));*/
      /*          if(MOVED(theVertex)) printf("Moved Node! %ld, xc yc %f %f\n",ID(theNode),XC(theVertex),YC(theVertex)); */
      if (MOVED(theVertex)) return(1);
    }
  }
  return(0);
}
/****************************************************************************/
/*
   NewPosCenterNode - calculate new global coordinates for a center node

   SYNOPSIS
   static INT NewPosCenterNode(ELEMENT *fatherElement, ELEMENT *nbElement[MAX_SIDES_OF_ELEM],
                            DOUBLE LimitLocDis, DOUBLE *newPos)

   PARAMETERS
   .  fatherElement     - father element of the center node
   .  nbElement         - neighbouring elements of the father element
   .  LimitLocDis       - maximum displacement in local coordinates
   .  newPos            - new position of the center nodes in global coordinates

   DESCRIPTION
   Calculate new coordinates for the center node of fatherElement according to the size
   of the neighbour elements.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

static INT NewPosCenterNode(ELEMENT *fatherElement, ELEMENT *nbElement[MAX_SIDES_OF_ELEM],
                            DOUBLE LimitLocDis, DOUBLE *newPos)
{
  DOUBLE sideMid[DIM],CenterPoint[DIM];
  DOUBLE *nbCorners[MAX_CORNERS_OF_ELEM],*fatherCorners[MAX_CORNERS_OF_ELEM],CenterPVertex[DIM];
  INT i,j,k,coe,fcorn,numOfSides,nmove[DIM];
  DOUBLE localmove,LocMove[DIM],LocCoord[DIM];
#ifdef __TWODIM__
  DOUBLE lcorn0[DIM],lcorn1[DIM];
  INT idim,found,co0,co1;
#endif

#ifdef __THREEDIM__
  PrintErrorMessage('E',"NewPosCenterNode","3D not implemented yet");
  return(1);
#endif

#ifdef __TWODIM__
  LocCoord[_X_] = 0.5;
  LocCoord[_Y_] = 0.5;
#else
  LocCoord[_X_] = 0.5;
  LocCoord[_Y_] = 0.5;
  LocCoord[_Z_] = 0.5;
#endif

  numOfSides = SIDES_OF_ELEM(fatherElement);
  CORNER_COORDINATES(fatherElement,fcorn,fatherCorners);
  /* calculate old global coordinates of center node according to xi=eta=0.5 */
  LOCAL_TO_GLOBAL(fcorn,fatherCorners,LocCoord,CenterPVertex);
  V_DIM_CLEAR(LocMove);
  V_DIM_CLEAR(nmove);

  /* loop over all neighbour elements */
  for (i=0; i<numOfSides; i++)
  {
    if (!nbElement[i]) continue;

    /* get coordinates of corners of neighbour element */
    CORNER_COORDINATES(nbElement[i],coe,nbCorners);
    /* center point of neigbour element */
    V_DIM_CLEAR(CenterPoint);
    for (j=0; j<coe; j++)
    {
      V_DIM_ADD1(nbCorners[j],CenterPoint);
    }
    V_DIM_SCALE(1.0/coe,CenterPoint);

    /* side mid of father element */
    V_DIM_CLEAR(sideMid);
    for (j=0; j<CORNERS_OF_SIDE(fatherElement,i); j++)
    {
      k = CORNER_OF_SIDE(fatherElement,i,j);
      V_DIM_ADD1(fatherCorners[k],sideMid);
    }
    V_DIM_SCALE(1./CORNERS_OF_SIDE(fatherElement,i),sideMid);

    /* calculate displacement in direction of the neighbour element */
    localmove = OneSideMoveCP(CenterPVertex,sideMid,CenterPoint);
#ifdef __TWODIM__
    co0 = CORNER_OF_SIDE(fatherElement,i,0);
    co1 = CORNER_OF_SIDE(fatherElement,i,1);
    /* determine displacement direction in local coordinates of the father element */
    UG_GlobalToLocal(fcorn,(const DOUBLE **)fatherCorners,
                     CVECT(MYVERTEX(CORNER(fatherElement,co0))),lcorn0);
    UG_GlobalToLocal(fcorn,(const DOUBLE **)fatherCorners,
                     CVECT(MYVERTEX(CORNER(fatherElement,co1))),lcorn1);
    found = FALSE;
    for (idim=0; idim<DIM; idim++)
    {
      if (LOCAL_EQUAL(lcorn0[idim],lcorn1[idim]))
      {
        /* idim is local coordinate to move, now determine direction in local coordinates */
        if (lcorn0[idim]<LocCoord[idim])
        {
          LocMove[idim] -= localmove;
          nmove[idim] ++;
          found = TRUE;
        }
        else
        {
          LocMove[idim] += localmove;
          nmove[idim] ++;
          found = TRUE;
        }
        break;
      }
    }
    if (found==FALSE)
    {
      printf("father corners %f  %f \n",fatherCorners[0][0],fatherCorners[0][1]);
      printf("father corners %f  %f \n",fatherCorners[1][0],fatherCorners[1][1]);
      printf("father corners %f  %f \n",fatherCorners[2][0],fatherCorners[2][1]);
      printf("father corners %f  %f \n",fatherCorners[3][0],fatherCorners[3][1]);
      printf("father corners %f  %f \n",CVECT(MYVERTEX(CORNER(fatherElement,co0)))[0],
             CVECT(MYVERTEX(CORNER(fatherElement,i)))[1]);
      UG_GlobalToLocal(fcorn,(const DOUBLE **)fatherCorners,
                       CVECT(MYVERTEX(CORNER(fatherElement,co0))),lcorn0);
      UG_GlobalToLocal(fcorn,(const DOUBLE **)fatherCorners,
                       CVECT(MYVERTEX(CORNER(fatherElement,co1))),lcorn1);
      printf("NewPosCenterNode lcorn0: %f %f, lcorn1: %f %f \n",lcorn0[0],lcorn0[1],lcorn1[0],lcorn1[1]);
    }
#endif
  }

  /* displacement in local coordinates according to all neighbour elements */
  for (i=0; i<DIM; i++)
    if (nmove[i]) LocMove[i] = LocMove[i]/nmove[i];

  /* new local coordinates of center node */
  V_DIM_ADD1(LocMove,LocCoord);
  /* limit moving distance */
  for (i=0; i<DIM; i++)
  {
    LocCoord[i] = MAX(LocCoord[i],0.5-LimitLocDis);
    LocCoord[i] = MIN(LocCoord[i],0.5+LimitLocDis);
  }

  /* if the father element is part of a curved boundary an additional move in boundary */
  /* direction must be performed                                                       */
  if (OBJT(fatherElement)==BEOBJ)
    if (MovedNode(fatherElement)==1)
      if (NewPosCenterNodeCurved(fatherElement,LocCoord)!=0) return(1);

  /* calculate new global coordinates of center node */
  LOCAL_TO_GLOBAL(fcorn,fatherCorners,LocCoord,newPos);
  return(0);
}

static INT LambdaFromTriangle (ELEMENT *theElement, NODE *cornerNodes[], DOUBLE *lambda)
{
  INT i,j,k,side0,side1,coe,nmove;
  NODE *node0, *node1;
  ELEMENT *nbElement;
  DOUBLE sideMid[DIM],CenterPointElem[DIM],CenterPointNBElem[DIM];
  DOUBLE *Corners[MAX_CORNERS_OF_ELEM],localmove0,localmove1;

  assert(CORNERS_OF_ELEM(theElement)==3);

  /* get coordinates of corners of theElement */
  CORNER_COORDINATES(theElement,coe,Corners);
  /* center point of theElement */
  V_DIM_CLEAR(CenterPointElem);
  for (j=0; j<coe; j++)
  {
    V_DIM_ADD1(Corners[j],CenterPointElem);
  }
  V_DIM_SCALE(1.0/coe,CenterPointElem);

  /* search neighbour element at side0 (S0-S2) and side1 (S1-S2) */
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    node0 = SONNODE(CORNER(theElement,CORNER_OF_SIDE(theElement,i,0)));
    node1 = SONNODE(CORNER(theElement,CORNER_OF_SIDE(theElement,i,1)));
    if ((node0==cornerNodes[0] && node1!=cornerNodes[1]) ||
        (node0!=cornerNodes[1] && node1==cornerNodes[0]))
    {
      side0 = i;
    }
    if ((node0!=cornerNodes[0] && node1==cornerNodes[1]) ||
        (node0==cornerNodes[1] && node1!=cornerNodes[0]))
    {
      side1 = i;
    }
  }
  localmove0 = localmove1 = 0; nmove = 0;
  if ((nbElement = NBELEM(theElement,side0))!=0)
  {
    /* get coordinates of corners of neighbour element */
    CORNER_COORDINATES(nbElement,coe,Corners);
    /* center point of neigbour element */
    V_DIM_CLEAR(CenterPointNBElem);
    for (j=0; j<coe; j++)
    {
      V_DIM_ADD1(Corners[j],CenterPointNBElem);
    }
    V_DIM_SCALE(1.0/coe,CenterPointNBElem);

    /* get coordinates of corners of theElement */
    CORNER_COORDINATES(theElement,coe,Corners);
    /* side mid of side0 */
    V_DIM_CLEAR(sideMid);
    for (j=0; j<CORNERS_OF_SIDE(theElement,side0); j++)
    {
      k = CORNER_OF_SIDE(theElement,side0,j);
      V_DIM_ADD1(Corners[k],sideMid);
    }
    V_DIM_SCALE(1./CORNERS_OF_SIDE(theElement,side0),sideMid);

    localmove0 = OneSideMoveCP(CenterPointElem,sideMid,CenterPointNBElem);
    nmove++;
  }
  if ((nbElement = NBELEM(theElement,side1))!=0)
  {
    /* get coordinates of corners of neighbour element */
    CORNER_COORDINATES(nbElement,coe,Corners);
    /* center point of neigbour element */
    V_DIM_CLEAR(CenterPointNBElem);
    for (j=0; j<coe; j++)
    {
      V_DIM_ADD1(Corners[j],CenterPointNBElem);
    }
    V_DIM_SCALE(1.0/coe,CenterPointNBElem);

    /* get coordinates of corners of theElement */
    CORNER_COORDINATES(theElement,coe,Corners);
    /* side mid of side1 */
    V_DIM_CLEAR(sideMid);
    for (j=0; j<CORNERS_OF_SIDE(theElement,side1); j++)
    {
      k = CORNER_OF_SIDE(theElement,side1,j);
      V_DIM_ADD1(Corners[k],sideMid);
    }
    V_DIM_SCALE(1./CORNERS_OF_SIDE(theElement,side1),sideMid);

    localmove1 = OneSideMoveCP(CenterPointElem,sideMid,CenterPointNBElem);
    nmove++;
  }
  if (nmove==0)
  {
    *lambda = 0.5;
    return(0);
  }
  *lambda = .5 + (localmove1 - localmove0)/nmove;
  return(0);
}

/****************************************************************************/
/*
   LambdaFromQuad - calculate new position of a mid node

   SYNOPSIS
   static INT LambdaFromQuad (ELEMENT *theElement, VERTEX *centerVertex,
                           NODE *sidemidNode,NODE *cornerNodes[], DOUBLE *lambda)

   PARAMETERS
   .  theElement     - father element of the center node
   .  centerVertex   - vertex of the center node
   .  sidemidNode    - mid node on one side of the father element
   .  cornerNodes    - corner nodes on the same side
   .  lambda         - new normalized position of sidemidNode on this side

   DESCRIPTION
   Calculate a new position of the side mid according to the local coordinates of the center node

   ------MS1-------S1
 |                |
 |   theElement   |
 |   (level l-1)  |
 |               OSM
 |                |
 |       CP      NSM
 |                |
   ------MS0-------S0

   CP: center point of theElement;  OSM: old side mid position;
   S1, S0:  corner nodes of side (level=l);
   In addition  MS1 and MS0 are needed when the element is on a curved boundary

   output: lambda (normalized position (NSM-S0)/(S1-S0) )


   RETURN VALUE
   INT
   .n   0: ok
 */
/**********************************************************************************/

static INT LambdaFromQuad (ELEMENT *theElement,VERTEX *centerVertex,
                           NODE *sidemidNode,NODE *cornerNodes[],
                           DOUBLE *lambda)
{
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM],lcorn0[DIM],lcorn1[DIM],lmid0[DIM],lmid1[DIM];
  DOUBLE *LocalCoord;
  INT i,j,k,coe,node_found,coord,curved;
  ELEMENT *sonElement;
  ELEMENT *SonList[MAX_SONS];
  NODE *midNode[2];
  LINK *theLink;

  assert(CORNERS_OF_ELEM(theElement)==4);
  /* local coordinates of center vertex */
  LocalCoord = LCVECT(centerVertex);
  /* local coordinates of corner vertices */
  CORNER_COORDINATES(theElement,coe,CornerPtrs);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,
                   CVECT(MYVERTEX(cornerNodes[0])),lcorn0);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,
                   CVECT(MYVERTEX(cornerNodes[1])),lcorn1);

  /* determine local coordinate of displacement */
  if (LOCAL_EQUAL(lcorn0[0],lcorn1[0]))
    coord = 1;
  else if (LOCAL_EQUAL(lcorn0[1],lcorn1[1]))
    coord = 0;
  else
  {
    printf("LambdaFromQuad lcorn0: %f %f, lcorn1: %f %f \n",lcorn0[0],lcorn0[1],lcorn1[0],lcorn1[1]);
    printf("center node nacher: xi=%f  eta=%f \n",LocalCoord[0],LocalCoord[1]);
    *lambda = 0.5;
    return(0);
  }

  /* determine direction and calculate lambda */
  if (lcorn0[coord]<lcorn1[coord])
    *lambda = LocalCoord[coord];
  else
    *lambda = 1-LocalCoord[coord];

  /* check for boundary element with a moved vertex */
  curved = FALSE;
  if (OBJT(theElement)==BEOBJ)
    if (MovedNode(theElement)) curved = TRUE;

  if (curved==FALSE) return(0);

  /* search for mid nodes linked to the cornerNodes */
  GetAllSons(theElement, SonList);
  for (k=0; k<2; k++)
  {
    node_found=FALSE;
    for (theLink=START(cornerNodes[k]); theLink!=0; theLink=NEXT(theLink))
    {
      for (i=0; i<NSONS(theElement); i++)
      {
        sonElement = SonList[i];
        for (j=0; j<CORNERS_OF_ELEM(sonElement); j++)
        {
          if (NBNODE(theLink)==CORNER(sonElement,j))
          {
            if (NBNODE(theLink)==sidemidNode) continue;
            node_found = TRUE;
            midNode[k] = NBNODE(theLink);
            break;
          }
        }
        if (node_found==TRUE) break;
      }
      if (node_found==TRUE) break;
    }
  }

  /* the two midnodes are not moved vertices, no change for lambda necessary  */
  if (!MOVED(MYVERTEX(midNode[0])) && !MOVED(MYVERTEX(midNode[1]))) return(0);

  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,
                   CVECT(MYVERTEX(midNode[0])),lmid0);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,
                   CVECT(MYVERTEX(midNode[1])),lmid1);
  /* determine local coordinate of displacement */
  if (LOCAL_EQUAL(lcorn0[0],lcorn1[0]))
    coord = 1;
  else if (LOCAL_EQUAL(lcorn0[1],lcorn1[1]))
    coord = 0;
  else
  {
    printf("LambdaFromQuadCurved lcorn0: %f %f, lcorn1: %f %f \n",lcorn0[0],lcorn0[1],lcorn1[0],lcorn1[1]);
    printf("center node nacher: xi=%f  eta=%f \n",LocalCoord[0],LocalCoord[1]);
    *lambda = 0.5;
    return(0);
  }

  *lambda = (LocalCoord[coord]-lmid0[coord])/(lmid1[coord]-lmid0[coord]);
  return(0);
}

static INT DefaultBndElemCenterLocal(NODE *centerNode, DOUBLE *defaultLocal)
{
  ELEMENT *fatherElement;
  DOUBLE global[DIM];
  INT coe,m,j;
  VERTEX *VertexOnEdge[MAX_EDGES_OF_ELEM];
  NODE *theNode;
  EDGE *theEdge;
  DOUBLE fac;
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM];

  if (NTYPE(centerNode)!=CENTER_NODE) return(1);
  fatherElement = VFATHER(MYVERTEX(centerNode));
  if (OBJT(fatherElement)!=BEOBJ) return(1);

  /* check if moved side nodes exist */
  if (MovedNode(fatherElement))
  {
    CORNER_COORDINATES(fatherElement,coe,CornerPtrs);
    m = EDGES_OF_ELEM(fatherElement);
    for (j=0; j<m; j++)
    {
      theEdge=GetEdge(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,j,0)),
                      CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,j,1)));
      assert(theEdge != NULL);
      theNode = MIDNODE(theEdge);
      if (theNode == NULL)
        VertexOnEdge[j] = NULL;
      else
        VertexOnEdge[j] = MYVERTEX(theNode);
    }
    fac = 1.0 / m;
    V_DIM_CLEAR(global);
    for (j=0; j<m; j++)
      if (VertexOnEdge[j] == NULL)
      {
        V_DIM_LINCOMB(1.0,global,0.5*fac,
                      CVECT(MYVERTEX(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,j,0)))),global);
        V_DIM_LINCOMB(1.0,global,0.5*fac,
                      CVECT(MYVERTEX(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,j,1)))),global);
      }
      else
        V_DIM_LINCOMB(1.0,global,fac,CVECT(VertexOnEdge[j]),global);
    UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,global,defaultLocal);
  }
  else
  {
    defaultLocal[0] = defaultLocal[1] = 0.5;
  }
  return(0);
}
#ifdef __TWODIM__
static INT LambdaOrthoBnd2D(const ELEMENT *fatherElement, const INT side, const NODE *MidNode,
                            const NODE *Node0, const NODE *Node1, const NODE *CenterNode,
                            DOUBLE *lambda)
{
  VERTEX *theVertex;
  DOUBLE_VECTOR BndPoint0,BndPoint1,midBndPoint,MPVec,TangVec,NormVec;
  DOUBLE Lambda0, midLambda, Lambda1,bndLambda[DIM-1],veclen ;
  DOUBLE *CenterPoint;
  DOUBLE_VECTOR VecA,VecB,VecC;
  DOUBLE prod0,prod1;
  BNDS *bnds;
  DOUBLE min_fac,area,min_diff;
  INT maxiter,iter,straight,orientation,concave;

  straight=FALSE;
  concave=FALSE;
  /* check if boundary is convex */
  V_DIM_SUBTRACT(CVECT(MYVERTEX(Node1)),CVECT(MYVERTEX(Node0)),VecA);
  V_DIM_SUBTRACT(CVECT(MYVERTEX(MidNode)),CVECT(MYVERTEX(Node0)),VecB);
  V_DIM_SUBTRACT(CVECT(MYVERTEX(CenterNode)),CVECT(MYVERTEX(Node0)),VecC);
  V_DIM_VECTOR_PRODUCT(VecA,VecB,prod0);
  V_DIM_VECTOR_PRODUCT(VecA,VecC,prod1);
  V_DIM_EUKLIDNORM(VecA,veclen);
  /* is the boundary a straight line ? */
  if (fabs(prod0/(veclen*veclen)) < 100*SMALL_C) straight=TRUE;
  if (prod0*prod1 < 0 && straight==FALSE) concave=TRUE;    /* boundary is concave */

  if ((prod1 > 0. && concave==FALSE) || (prod1 < 0. && concave==TRUE))
    orientation = 1;     /* CenterPoint and midNode are left from line Node0->Node1 */
  else
    orientation = 2;     /* CenterPoint and midNode are right from line Node0->Node1 */

  CenterPoint = CVECT(MYVERTEX(CenterNode));
  min_fac = 0.0001;
  maxiter = 50;

  theVertex = MYVERTEX(MidNode);
  if (OBJT(theVertex) != BVOBJ) return(GM_ERROR);

  bnds = ELEM_BNDS(fatherElement,side);
  min_diff = min_fac;

  Lambda0 = 0.;
  Lambda1 = 1.;
  iter=0;

  /*     printf("midnode %ld  ",ID(MidNode));  */
  while (fabs(Lambda1-Lambda0)>fabs(5.*min_diff) && iter<maxiter)
  {
    iter++;
    midLambda = 0.5*(Lambda0+Lambda1);
    do
    {
      bndLambda[0] = midLambda-min_diff;
      bndLambda[0] = MAX(bndLambda[0],0.); bndLambda[0] = MIN(bndLambda[0],1.);
      BNDS_Global(bnds,bndLambda,BndPoint0);
      bndLambda[0]=midLambda+min_diff;
      bndLambda[0] = MAX(bndLambda[0],0.); bndLambda[0] = MIN(bndLambda[0],1.);
      BNDS_Global(bnds,bndLambda,BndPoint1);
      min_diff *=2;
    }
    while (V_DIM_ISEQUAL(BndPoint0,BndPoint1));
    /* reset min_diff */
    min_diff = min_fac;
    bndLambda[0]=midLambda;
    BNDS_Global(bnds,bndLambda,midBndPoint);

    V_DIM_SUBTRACT(BndPoint1,BndPoint0,TangVec);
    V_DIM_EUKLIDNORM(TangVec,veclen);

    V_DIM_SCALE(1./veclen,TangVec);

    /* vector normal to the boundary with direction into domain */
    if (orientation==1)
    {
      NormVec[0] = - TangVec[1]; NormVec[1] = TangVec[0];
    }
    else
    {
      NormVec[0] = TangVec[1]; NormVec[1] = -TangVec[0];
    }
    /* calculate direction MidNode -> CenterNode  */
    V_DIM_SUBTRACT(CenterPoint,midBndPoint,MPVec);
    V_DIM_EUKLIDNORM(MPVec,veclen);
    V_DIM_SCALE(1./veclen,MPVec);

    /* calculate vector product */
    V_DIM_VECTOR_PRODUCT(MPVec,NormVec,area);
    if (area==0.) break;
    if (concave==TRUE)
    {
      if ( (area>0 && orientation==1) || (area<0 && orientation==2))
      {
        Lambda0 = midLambda;
      }
      else
      {
        Lambda1 = midLambda;
      }
    }
    else
    {
      if ( (area>0 && orientation==1) || (area<0 && orientation==2))
      {
        Lambda0 = midLambda;
      }
      else
      {
        Lambda1 = midLambda;
      }
    }
  }

  bndLambda[0] = 0.;
  BNDS_Global(bnds,bndLambda,BndPoint0);

  if (V_DIM_ISEQUAL(CVECT(MYVERTEX(Node0)),BndPoint0))
    *lambda = midLambda;
  else if (V_DIM_ISEQUAL(CVECT(MYVERTEX(Node1)),BndPoint0))
    *lambda = 1.-midLambda;
  else
    return(1);
  if (iter > 40)
    printf ("iter %d, midnode-id %ld, midlambda %f, lambda %f,  area %lf\n",
            iter,ID(MidNode),midLambda,*lambda,area);
  return(0);
}
#endif
/*******************************************************************/
/*                                                                 */
/*                                                                 */
/*   CornerNodes[1]       SideNodes[2]         CornerNodes[0]      */
/*     #----------------------#-------------------#                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |            CenterNode|                   |                */
/*     #----------------------#-------------------#SideNodes[0]    */
/*     |SideNodes[1]          |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     #----------------------#-------------------#                */
/*    BndNodes[1]           BndNode             BndNodes[0]        */
/*                                                                 */
/*******************************************************************/
static INT EqualDistBndElem(ELEMENT *fatherElement, INT side, NODE *CenterNode, DOUBLE *CenterNodeCoord,
                            NODE **SideNode, DOUBLE *SideNodeCoord)
{
  NODE *SideNodes[3], *BndNodes[2], *BndNode, *CornerNodes[2];
  LINK *theLink;
  DOUBLE dist0, dist1, distmax, CenterNodeDist, SideNodeDist, Lambda0, Lambda1, midLambda, bndLambda[DIM_OF_BND];
  DOUBLE *BndPoint, BndPoint0[DIM],BndPoint1[DIM];
  DOUBLE TangVec[DIM],NormVec[DIM],veclen ;
  DOUBLE OldNormVec[DIM],prod,min_diff,min_fac;
  INT i,j,k,nlinks;
  BNDS *bnds;

  /* find side nodes */
  i = nlinks = 0;
  for (theLink=START(CenterNode); theLink!=0; theLink=NEXT(theLink))
  {
    if (OBJT(MYVERTEX(NBNODE(theLink)))==BVOBJ)
      BndNode = NBNODE(theLink);
    else
      SideNodes[i++] = NBNODE(theLink);
    nlinks++;
  }
  assert(nlinks==4);
  assert(i==3);

  /* find boundary nodes of father element */
  j = k = 0;
  for (i=0; i<CORNERS_OF_ELEM(fatherElement); i++)
  {
    if (OBJT(MYVERTEX(CORNER(fatherElement,i)))==BVOBJ)
      BndNodes[j++] = CORNER(fatherElement,i);
    else
      CornerNodes[k++] = CORNER(fatherElement,i);
  }
  assert(j==2);
  assert(k==2);

  /* find lambda of BndNode */
  Lambda0 = 0.;
  Lambda1 = 1.;
  BndPoint = CVECT(MYVERTEX(BndNode));
  bnds = ELEM_BNDS(fatherElement,side);
  do
  {
    midLambda = 0.5*(Lambda0+Lambda1);
    bndLambda[0] = Lambda0;
    BNDS_Global(bnds,bndLambda,BndPoint0);
    bndLambda[0]= midLambda;
    BNDS_Global(bnds,bndLambda,BndPoint1);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint1,BndPoint0,dist1);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint,BndPoint0,dist0);
    if (dist0<dist1)
      Lambda1 = midLambda;
    else
      Lambda0 = midLambda;
  }
  while (!V_DIM_ISEQUAL(BndPoint0,BndPoint1));

  /* calculate new boundary distance of CenterNode */
  CenterNodeDist = distmax = 0;
  for (i=0; i<3; i++)
  {
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[0])),dist0);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[1])),dist1);
    dist0 = MIN(dist0,dist1);
    if (dist0>distmax)
    {
      distmax = dist0;
      SideNode[0] = SideNodes[i];
    }
    CenterNodeDist += dist0;
  }
  CenterNodeDist = (CenterNodeDist-distmax)*0.5;

  /* calculate new boundary distance of SideNode[2] */
  SideNodeDist = 0.;
  for (i=0; i<2; i++)
  {
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CornerNodes[i])),CVECT(MYVERTEX(BndNodes[0])),dist0);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CornerNodes[i])),CVECT(MYVERTEX(BndNodes[1])),dist1);
    dist0 = MIN(dist0,dist1);
    SideNodeDist += dist0;
  }
  SideNodeDist *= 0.5;

  /* calculate tangential vector */
  min_fac = 0.0001;
  min_diff = min_fac;
  do
  {
    bndLambda[0] = midLambda-min_diff;
    bndLambda[0] = MAX(bndLambda[0],0.); bndLambda[0] = MIN(bndLambda[0],1.);
    BNDS_Global(bnds,bndLambda,BndPoint0);
    bndLambda[0]=midLambda+min_diff;
    bndLambda[0] = MAX(bndLambda[0],0.); bndLambda[0] = MIN(bndLambda[0],1.);
    BNDS_Global(bnds,bndLambda,BndPoint1);
    min_diff *=2;
  }
  while (V_DIM_ISEQUAL(BndPoint0,BndPoint1));

  V_DIM_SUBTRACT(BndPoint1,BndPoint0,TangVec);
  V_DIM_EUKLIDNORM(TangVec,veclen);

  V_DIM_SCALE(1./veclen,TangVec);

  /* old normal vector */
  V_DIM_SUBTRACT(CVECT(MYVERTEX(CenterNode)),CVECT(MYVERTEX(BndNode)),OldNormVec);
  /* find correct orientation for normal vector */
  NormVec[0] = - TangVec[1]; NormVec[1] = TangVec[0];
  V_DIM_SCALAR_PRODUCT(OldNormVec,NormVec,prod);
  if (prod<0)
  {
    V_DIM_SCALE(-1.,NormVec);
  }

  V_DIM_LINCOMB(1.,BndPoint,CenterNodeDist,NormVec,CenterNodeCoord);
  V_DIM_LINCOMB(1.,BndPoint,SideNodeDist,NormVec,SideNodeCoord);

  return(0);
}

/*******************************************************************/
/*        fatherElement                                            */
/*        -------------                                            */
/*   CornerNodes[1]       SideNodes[2]         CornerNodes[0]      */
/*     #----------------------#-------------------#                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |            CenterNode|                   |                */
/*     #----------------------#-------------------#SideNodes[0]    */
/*     |SideNodes[1]          |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     #----------------------#-------------------#                */
/*                         SideNodes[3]                            */
/*                                                                 */
/*        boundaryFather                                           */
/*        --------------                                           */
/*                                                                 */
/*     #----------------------#-------------------#                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     #----------------------#-------------------#                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     |                      |                   |                */
/*     #----------------------#-------------------#                */
/*    BndNodes[1]           BndNode             BndNodes[0]        */
/*                                                                 */
/*******************************************************************/
static INT EqualDistInnerElem(ELEMENT *fatherElement, ELEMENT *boundaryFather, INT side, NODE *CenterNode,
                              DOUBLE *CenterNodeCoord, NODE **SideNode, DOUBLE *SideNodeCoord)
{
  NODE *SideNodes[4], *BndNodes[2], *BndNode, *CornerNodes[2], *NearestNode;
  LINK *theLink,*theLink1;
  DOUBLE dist0, dist1, dist2, distmax, distmin, CenterNodeDist, SideNodeDist;
  DOUBLE Lambda0, Lambda1, midLambda, bndLambda[DIM_OF_BND];
  DOUBLE *BndPoint, BndPoint0[DIM],BndPoint1[DIM];
  DOUBLE TangVec[DIM],NormVec[DIM],veclen ;
  DOUBLE OldNormVec[DIM],prod,min_diff,min_fac;
  INT i,j,iter;
  BNDS *bnds;

  /* find side nodes */
  i = 0;
  for (theLink=START(CenterNode); theLink!=0; theLink=NEXT(theLink))
    SideNodes[i++] = NBNODE(theLink);
  assert(i==4);

  /* find boundary nodes of boundary father element  */
  j = 0;
  for (i=0; i<CORNERS_OF_ELEM(boundaryFather); i++)
  {
    if (OBJT(MYVERTEX(CORNER(boundaryFather,i)))==BVOBJ)
      BndNodes[j++] = CORNER(boundaryFather,i);
  }
  assert(j==2);
  /* find son boundary node of boundary father element */
  i = 0;
  for (theLink=START(SONNODE(BndNodes[0])); theLink!=0; theLink=NEXT(theLink))
  {
    for (theLink1=START(SONNODE(BndNodes[1])); theLink1!=0; theLink1=NEXT(theLink1))
    {
      if (NBNODE(theLink)==NBNODE(theLink1))
      {
        BndNode = NBNODE(theLink);
        i++;
      }
    }
  }
  assert(i==1);
  assert(OBJT(MYVERTEX(BndNode))==BVOBJ);

  /* find lambda of BndNode */
  Lambda0 = 0.;
  Lambda1 = 1.;
  iter=0;
  BndPoint = CVECT(MYVERTEX(BndNode));
  bnds = ELEM_BNDS(boundaryFather,side);
  do
  {
    midLambda = 0.5*(Lambda0+Lambda1);
    bndLambda[0] = Lambda0;
    BNDS_Global(bnds,bndLambda,BndPoint0);
    bndLambda[0]= midLambda;
    BNDS_Global(bnds,bndLambda,BndPoint1);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint1,BndPoint0,dist1);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint,BndPoint0,dist0);
    if (dist0<dist1)
      Lambda1 = midLambda;
    else
      Lambda0 = midLambda;
  }
  while (!V_DIM_ISEQUAL(BndPoint0,BndPoint1));

  /* order side nodes */
  /* side node with maximum wall distance */
  distmax = 0;
  for (i=0; i<4; i++)
  {
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[0])),dist0);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[1])),dist1);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNode)),dist2);
    dist0 = MIN(MIN(dist0,dist1),dist2);
    if (dist0>distmax)
    {
      SideNode[0] = SideNodes[i];
      distmax = dist0;
    }
  }
  /* side node with minimum wall distance */
  distmin = MAX_D;
  for (i=0; i<4; i++)
  {
    /*         V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[0])),dist0); */
    /*         V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[1])),dist1); */
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNode)),dist2);
    dist0 = dist2;     /* MIN(MIN(dist0,dist1),dist2); */
    if (dist0<distmin)
    {
      NearestNode = SideNodes[i];
      distmin = dist0;
    }
  }

  /* calculate new boundary distance of CenterNode */
  CenterNodeDist = 0;
  for (i=0; i<4; i++)
  {
    if (SideNodes[i]==NearestNode || SideNodes[i]==SideNode[0]) continue;
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[0])),dist0);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(SideNodes[i])),CVECT(MYVERTEX(BndNodes[1])),dist1);
    CenterNodeDist += MIN(dist0,dist1);
  }
  CenterNodeDist *= 0.5;

  /* calculate new boundary distance of farest SideNode[2] */
  CornerNodes[0] = CORNER(fatherElement,(side+2)%4);
  CornerNodes[1] = CORNER(fatherElement,(side+3)%4);
  SideNodeDist = 0.;
  for (i=0; i<2; i++)
  {
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CornerNodes[i])),CVECT(MYVERTEX(BndNodes[0])),dist0);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CornerNodes[i])),CVECT(MYVERTEX(BndNodes[1])),dist1);
    dist0 = MIN(dist0,dist1);
    SideNodeDist += dist0;
  }
  SideNodeDist *= 0.5;

  /* calculate tangential vector */
  min_fac = 0.0001;
  min_diff = min_fac;
  do
  {
    bndLambda[0] = midLambda-min_diff;
    bndLambda[0] = MAX(bndLambda[0],0.); bndLambda[0] = MIN(bndLambda[0],1.);
    BNDS_Global(bnds,bndLambda,BndPoint0);
    bndLambda[0]=midLambda+min_diff;
    bndLambda[0] = MAX(bndLambda[0],0.); bndLambda[0] = MIN(bndLambda[0],1.);
    BNDS_Global(bnds,bndLambda,BndPoint1);
    min_diff *=2;
  }
  while (V_DIM_ISEQUAL(BndPoint0,BndPoint1));

  V_DIM_SUBTRACT(BndPoint1,BndPoint0,TangVec);
  V_DIM_EUKLIDNORM(TangVec,veclen);

  V_DIM_SCALE(1./veclen,TangVec);

  /* old normal vector */
  V_DIM_SUBTRACT(CVECT(MYVERTEX(CenterNode)),CVECT(MYVERTEX(BndNode)),OldNormVec);
  /* find correct orientation for normal vector */
  NormVec[0] = - TangVec[1]; NormVec[1] = TangVec[0];
  V_DIM_SCALAR_PRODUCT(OldNormVec,NormVec,prod);
  if (prod<0)
  {
    V_DIM_SCALE(-1.,NormVec);
  }

  V_DIM_LINCOMB(1.,BndPoint,CenterNodeDist,NormVec,CenterNodeCoord);
  V_DIM_LINCOMB(1.,BndPoint,SideNodeDist,NormVec,SideNodeCoord);

  return(0);
}


/****************************************************************************/
/*
   SmoothGrid - resize all quadrilaterals and triangles on grid (l) according to the element sizes of
             grid (l-1)

   SYNOPSIS
   static INT SmoothGrid (GRID *theGrid, DOUBLE LimitLocDis, INT *MoveInfo)

   PARAMETERS
   .  theGrid      - resize quadrilaterals and triangles on this grid
   .  LimitLocDis  - maximum displacement of the nodes in local coordinates  (0<LimitLocDis<0.5)
   .  MoveInfo     - give information about moved nodes

   DESCRIPTION
   Resize all quadrilaterals and triangles on grid (l) according to the element sizes of grid (l-1).
   In the first node-loop all center nodes will be modified according to the size of the neighbour elements
   (grid l-1). In the second node-loop all mid nodes will be modified according to the position of the
   (moved) center nodes and the size of the triangles.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

INT SmoothGrid (GRID *theGrid, const DOUBLE LimitLocDis, INT *MoveInfo, const INT ForceLevelSet,
                const INT bnd_num, const INT *bnd)
{
  MULTIGRID *theMG;
  ELEMENT *fatherElement,*nbElement[MAX_SIDES_OF_ELEM];
  ELEMENT *SonList[MAX_SONS];
  VERTEX *theVertex;
  NODE *theNode;
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM];
  DOUBLE newPos[DIM],newLocal[DIM];
  INT i,coe,numOfSides,OnlyRedSons;
#ifdef __TWODIM__
  NODE *theSideNode[1];
  LINK *theLink;
  NODE *CornerNodes[2],*CenterNodes[2],*node0,*node1;
  ELEMENT *oppositeElement,*Level0Father,*boundaryFather;
  DOUBLE lambda,lambda0,lambda1,lambda_old,x1,x2;
  INT ceN,Eside,nlinks,j,type[3],side,sidefound,eqdist;
  DOUBLE dummyValues[10];
  DOUBLE CenterNodeCoord[DIM],SideNodeCoord[DIM];
  static DOUBLE local[] = {0.5};
  INT edge,co0,co1,ortho_used;
#endif

  theMG = MYMG(theGrid);

  /*    move center nodes of quadrilaterals  */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    /* skip node if it is a copy from a lower level */
    if (CORNERTYPE(theNode)) continue;

    /* skip node if it is not a center node */
    if (NTYPE(theNode)!=CENTER_NODE) continue;

    /* use only quadrilaterals */
    theVertex=MYVERTEX(theNode);
    fatherElement = VFATHER(theVertex);

    if (!ForceLevelSet)
    {
      /* move nodes below currentlevel only if irregular or no refinement was performed on
         the father element before */
      if (REFINE(fatherElement)==RED && LEVEL(theNode)<CURRENTLEVEL(theMG))
      {
        OnlyRedSons = TRUE;
        GetAllSons(fatherElement, SonList);
        for (i=0; i<NSONS(fatherElement); i++)
        {
          if (REFINE(SonList[i])!=RED) OnlyRedSons=FALSE;
        }
        if (OnlyRedSons==TRUE)
        {
          /*                  printf("node-id %ld , father-elem %ld wurde uebersprungen \n",ID(theNode),ID(fatherElement)); */
          continue;
        }
      }
    }

#ifdef __TWODIM__
    if (TAG(fatherElement)!=QUADRILATERAL) continue;      /* obsolete */
#else
    PrintErrorMessage('E',"SmoothGrid","3D not implemented yet");
    return(1);
#endif

    /* create list of neighbour elements */
    numOfSides = SIDES_OF_ELEM(fatherElement);
    for (i=0; i<numOfSides; i++)
    {
      if (NBELEM(fatherElement,i))
      {
        nbElement[i] = NBELEM(fatherElement,i);
      }
      else
        nbElement[i] = 0;
    }

    /* calculate new global coordinates of center node*/
    if (NewPosCenterNode(fatherElement,nbElement,LimitLocDis,newPos)!=0) return(1);

    /* calculate new local coordinates of center node */
    CORNER_COORDINATES(fatherElement,coe,CornerPtrs);
    UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs, newPos,newLocal);

    /* move node only when necessary */
    if (V2_LOCAL_EQUAL(LCVECT(MYVERTEX(theNode)),newLocal)!=0) continue;

    /* move center node */
    if (MoveNode(theMG,theNode,newPos)!=0) return(1);
    MoveInfo[0]++;
    /* are there nodes which reached the moving limit ? */
    if ( (LOCAL_EQUAL(XI(theVertex),(0.5+LimitLocDis)) || LOCAL_EQUAL(XI(theVertex),(0.5-LimitLocDis))) ||
         (LOCAL_EQUAL(ETA(theVertex),(0.5+LimitLocDis)) || LOCAL_EQUAL(ETA(theVertex),(0.5-LimitLocDis))))
    {
      printf("center-node %ld, father-elem%ld reached limit\n",ID(theNode),ID(fatherElement));
      MoveInfo[2]++;
    }
  }

  /* search side vertices */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    /* skip node if it is a copy from a lower level */
    if (CORNERTYPE(theNode)) continue;

    /* skip node if it is not a mid node (mid point of an edge) */
    if (NTYPE(theNode)!=MID_NODE) continue;

    /* search for one or two elements on level i-1 */
#ifdef __TWODIM__
    theVertex = MYVERTEX(theNode);

    /* find the two corner nodes on the same edge and the elements on this edge */
    fatherElement = VFATHER(theVertex);
    edge = ONEDGE(theVertex);
    co0 = CORNER_OF_EDGE(fatherElement,edge,0);
    co1 = CORNER_OF_EDGE(fatherElement,edge,1);
    CornerNodes[0] = SONNODE(CORNER(fatherElement,co0));
    CornerNodes[1] = SONNODE(CORNER(fatherElement,co1));

    /* find the center nodes of the neighbouring elements */
    ceN = 0;
    nlinks = 0;
    for (theLink=START(theNode); theLink!=0; theLink=NEXT(theLink))
    {
      /*            if (NTYPE(NBNODE(theLink))==CORNER_NODE)
                  {
                      CornerNodes[coN] = NBNODE(theLink);
                      coN++;
                  }  */
      if (NTYPE(NBNODE(theLink))==CENTER_NODE)
      {
        CenterNodes[ceN] = NBNODE(theLink);
        ceN++;
      }
      nlinks++;
    }

    /* 5 possibilities of element neighbourship relationships */

    ortho_used = FALSE;
    /* quadrilateral is boundary element */
    if (nlinks==3 && ceN==1)
    {
      fatherElement = VFATHER(MYVERTEX(CenterNodes[0]));
      if (OBJT(fatherElement)!=BEOBJ) printf("elem-id %ld, node-id %ld, X %f %f\n",ID(fatherElement),ID(theNode),XC(MYVERTEX(theNode)),YC(MYVERTEX(theNode)));
      if (bnd_num==0 || OBJT(MYVERTEX(theNode))!=BVOBJ)
      {
        if (LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),theNode,CornerNodes,&lambda)!=0)
          continue;
      }
      else
      {
        /* move boundary side nodes of quadrilaterals in order to get orthogonal boundary elements  on
           specified boundaries */
        sidefound = FALSE;
        for (i=0; i<SIDES_OF_ELEM(fatherElement); i++)
        {
          if (!SIDE_ON_BND(fatherElement,i) || INNER_BOUNDARY(fatherElement,i)) continue;
          SideBndCond(fatherElement,i,local,dummyValues,type);
          for (j=0; j<bnd_num; j++)
            if (type[1]==bnd[j])
            {
              sidefound = TRUE;
              side = i;
            }
        }
        if (sidefound==FALSE)
        {
          if (LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),theNode,CornerNodes,&lambda)!=0)
            continue;
        }
        else
        {
          /* apply 'ortho' option for boundary mid-node if possible */
          if (LambdaOrthoBnd2D(fatherElement,side,theNode,(NODE *)NFATHER(CornerNodes[0]),
                               (NODE *)NFATHER(CornerNodes[1]),CenterNodes[0],&lambda)!=0)
          {
            if (LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),theNode,CornerNodes,&lambda)!=0)
              continue;
          }
          else
            ortho_used=TRUE;
        }
      }

    }
    /* triangle is boundary element */
    else if (nlinks==4 && ceN==0)
    {
      fatherElement = VFATHER(MYVERTEX(theNode));
      if (LambdaFromTriangle(fatherElement,CornerNodes,&lambda)!=0) continue;
    }
    /* two quadrilaterals */
    else if (nlinks==4 && ceN==2)
    {
      fatherElement = VFATHER(MYVERTEX(CenterNodes[0]));
      oppositeElement = VFATHER(MYVERTEX(CenterNodes[1]));
      if(LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),theNode,CornerNodes,&lambda0)!=0)
        continue;
      if (LambdaFromQuad(oppositeElement,MYVERTEX(CenterNodes[1]),theNode,CornerNodes,&lambda1)!=0)
        continue;
      lambda = 0.5*(lambda0+lambda1);
    }
    /* two triangles */
    else if (nlinks==6 && ceN==0)
    {
      fatherElement = VFATHER(MYVERTEX(theNode));
      /* opposite element is a triangle */
      for (i=0; i<SIDES_OF_ELEM(fatherElement); i++)
      {
        node0 = SONNODE(CORNER(fatherElement,CORNER_OF_SIDE(fatherElement,i,0)));
        node1 = SONNODE(CORNER(fatherElement,CORNER_OF_SIDE(fatherElement,i,1)));
        if ((node0==CornerNodes[0] && node1==CornerNodes[1]) ||
            (node0==CornerNodes[1] && node1==CornerNodes[0]))
        {
          Eside = i;
          continue;
        }
      }
      if (TAG(fatherElement)!=TRIANGLE)
      {
        lambda0 = 0.5;
      }
      else
      {
        if (LambdaFromTriangle(fatherElement,CornerNodes,&lambda0)!=0) continue;
      }
      oppositeElement = NBELEM(fatherElement,Eside);
      if (TAG(oppositeElement)!=TRIANGLE)
      {
        lambda1 = 0.5;
      }
      else
      {
        if (LambdaFromTriangle(oppositeElement,CornerNodes,&lambda1)!=0) continue;
      }
      lambda = 0.5*(lambda0+lambda1);
    }
    /* one triangle and one quadrilateral */
    else if (nlinks==5 && ceN==1)
    {
      /* quadrilateral */
      fatherElement = VFATHER(MYVERTEX(CenterNodes[0]));
      /* opposite element is a triangle */
      for (i=0; i<SIDES_OF_ELEM(fatherElement); i++)
      {
        node0 = SONNODE(CORNER(fatherElement,CORNER_OF_SIDE(fatherElement,i,0)));
        node1 = SONNODE(CORNER(fatherElement,CORNER_OF_SIDE(fatherElement,i,1)));
        if ((node0==CornerNodes[0] && node1==CornerNodes[1]) ||
            (node0==CornerNodes[1] && node1==CornerNodes[0]))
        {
          Eside = i;
          continue;
        }
      }
      if(LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),theNode,CornerNodes,&lambda0)!=0)
        continue;
      oppositeElement = NBELEM(fatherElement,Eside);
      if (TAG(oppositeElement)!=TRIANGLE)
      {
        lambda1 = 0.5;
      }
      else
      {
        if (LambdaFromTriangle(oppositeElement,CornerNodes,&lambda1)!=0) continue;
      }
      lambda = 0.5*(lambda0+lambda1);
    }
    /* one quadrilateral and one irregular refined quadrilateral */
    else if (ceN==1)
    {
      /* quadrilateral */
      fatherElement = VFATHER(MYVERTEX(CenterNodes[0]));
      if(LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),theNode,CornerNodes,&lambda0)!=0)
        continue;
      lambda = 0.5*(lambda0 + 0.5);
    }
    else
    {
      /* durch adaptives refinement koennen auch noch andere Faelle auftreten */
      printf("Dieser Fall ist bisher noch nicht vorgesehen: ceN=%d, nlinks=%d, MidNode %ld, VaterEl %ld \n"
             ,ceN,nlinks,ID(theNode),ID(fatherElement));
      continue;
    }

    /* apply limits for lambda  */
    lambda = MAX(lambda,0.5-LimitLocDis);
    lambda = MIN(lambda,0.5+LimitLocDis);

    /* move node only when necessary */
    if (OBJT(MYVERTEX(theNode))==BVOBJ)
    {
      if (GetMidNodeParam(theNode,&lambda_old)!=0) return(1);
    }
    else
    {
      /* calculate old lambda */
      V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CornerNodes[1])),CVECT(MYVERTEX(CornerNodes[0])),x1);
      V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(theNode)),CVECT(MYVERTEX(CornerNodes[0])),x2);
      lambda_old = x2/x1;
    }
    if (LOCAL_EQUAL(lambda,lambda_old)) continue;

    /* move boundary node */
    if (MoveMidNode(theMG,theNode,lambda)!=0) return(1);
    MoveInfo[1]++;
    if (ortho_used) MoveInfo[4]++;
    if (LOCAL_EQUAL(lambda,(0.5+LimitLocDis)) || LOCAL_EQUAL(lambda,(0.5-LimitLocDis)) )
    {
      printf("mid-node %ld, father-elem%ld reached limit\n",ID(theNode),ID(fatherElement));
      MoveInfo[3]++;
    }
#else
    PrintErrorMessage('E',"SmoothGrid","3D not implemented yet");
    return(1);
#endif
  }

#ifdef __TWODIM__
  /* move all son nodes of an level-0-boundary-element in order to get 'orthogonal' elements */
  eqdist = TRUE;
  if (eqdist==TRUE && bnd_num > 0)
  {
    /*    move center nodes of quadrilaterals  */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* skip node if it is a copy from a lower level */
      if (CORNERTYPE(theNode)) continue;

      /* skip node if it is not a center node */
      if (NTYPE(theNode)!=CENTER_NODE) continue;

      /* use only boundary father elements  */
      /* the father element on level 0 must be a boundary element */
      theVertex=MYVERTEX(theNode);
      fatherElement = VFATHER(theVertex);

      /* skip node if it is not on a regular refined element */
      if (REFINECLASS(fatherElement)!=RED_CLASS) continue;

      Level0Father=fatherElement;
      for (i=GLEVEL(theGrid)-1; i>0; i--) Level0Father=EFATHER(Level0Father);     /* new */
      if (OBJT(Level0Father)!=BEOBJ) continue;        /* new */
      /*             if (OBJT(fatherElement)!=BEOBJ) continue; */

      if (OBJT(fatherElement)==BEOBJ)
      {
        /* check if father element is on wall boundary */
        sidefound = FALSE;
        for (i=0; i<SIDES_OF_ELEM(fatherElement); i++)
        {
          if (!SIDE_ON_BND(fatherElement,i) || INNER_BOUNDARY(fatherElement,i)) continue;
          SideBndCond(fatherElement,i,local,dummyValues,type);
          for (j=0; j<bnd_num; j++)
            if (type[1]==bnd[j])
            {
              sidefound = TRUE;
              side = i;
            }
        }
        if (sidefound==FALSE) continue;

        if (EqualDistBndElem(fatherElement,side,theNode,CenterNodeCoord,theSideNode,SideNodeCoord)!=0) continue;

        /* move center node */
        if (MoveNode(theMG,theNode,CenterNodeCoord)!=0) return(1);
        /*             MoveInfo[0]++; */
        /* move side node */
        if (MoveNode(theMG,theSideNode[0],SideNodeCoord)!=0) return(1);
        /*             MoveInfo[1]++; */
      }
      else
      {
        /* check if father element is on wall boundary */
        sidefound = FALSE;
        for (i=0; i<SIDES_OF_ELEM(Level0Father); i++)
        {
          if (!SIDE_ON_BND(Level0Father,i) || INNER_BOUNDARY(Level0Father,i)) continue;
          SideBndCond(Level0Father,i,local,dummyValues,type);
          for (j=0; j<bnd_num; j++)
            if (type[1]==bnd[j])
            {
              sidefound = TRUE;
              side = i;
            }
        }
        if (sidefound==FALSE) continue;

        /* find corresponding element at the boundary */
        for (boundaryFather=NBELEM(fatherElement,side); boundaryFather!=0;
             boundaryFather=NBELEM(boundaryFather,side))
        {
          /* otherwise we will never find the correct boundary element */
          if (REFINECLASS(boundaryFather)!=RED_CLASS) break;
          if(OBJT(boundaryFather)==BEOBJ)
          {
            if (EqualDistInnerElem(fatherElement,boundaryFather,side,theNode,
                                   CenterNodeCoord,theSideNode,SideNodeCoord)!=0) continue;
            /* move center node */
            if (MoveNode(theMG,theNode,CenterNodeCoord)!=0) return(1);
            /*             MoveInfo[0]++; */
            /* move side node */
            if (MoveNode(theMG,theSideNode[0],SideNodeCoord)!=0) return(1);
            /*             MoveInfo[1]++; */
          }
        }
      }
    }
  }
#endif
  return(0);
}
/****************************************************************************/
/*
   SmoothGridReset - reset all nodes on surface levels to the default refinement position

   SYNOPSIS
   static INT SmoothGridReset (GRID *theGrid, INT *MoveInfo)

   PARAMETERS
   .  theGrid      - resize quadrilaterals and triangles on this grid
   .  MoveInfo     - give information about reset nodes

   DESCRIPTION
   Reset all center nodes and mid nodes of the quadrilaterals of the grid to the position after
   the grid refinement. This means in local coordinates: center node (0.5,0.5) and mid nodes (0.5,0),
   (1,0.5), (0.5,1), (0,0.5). This works for triangles too.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

INT SmoothGridReset (GRID *theGrid, INT *MoveInfo)
{
  MULTIGRID *theMG;
  ELEMENT *fatherElement;
  VERTEX *theVertex;
  NODE *theNode;
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM],LocalCenter[3]={0.5,0.5,0.5};
  DOUBLE newPos[DIM],defaultLocal[DIM];
  INT coe;
#ifdef __TWODIM__
  LINK *theLink;
  NODE *CornerNodes[2],*CenterNodes[2];
  DOUBLE lambda_old,x1,x2;
  INT ceN,coN;
#endif

  theMG = MYMG(theGrid);
  /* re-move side vertices; this must be performed before re-moving the center nodes because of
     DefaultBndElemCenterLocal */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    /* skip node if it is a copy from a lower level */
    if (CORNERTYPE(theNode)) continue;

    /* skip node if it is not a mid node (mid point of an edge) */
    if (NTYPE(theNode)!=MID_NODE) continue;

#ifdef __TWODIM__
    /* find the two corner nodes on the same edge */
    coN = 0;
    ceN = 0;
    for (theLink=START(theNode); theLink!=0; theLink=NEXT(theLink))
    {
      if (CORNERTYPE(NBNODE(theLink)))
      {
        CornerNodes[coN] = NBNODE(theLink);
        coN++;
      }
      if (NTYPE(NBNODE(theLink))==CENTER_NODE)
      {
        CenterNodes[ceN] = NBNODE(theLink);
        ceN++;
      }
    }
    assert (coN==2);
    if (coN>2) {PrintErrorMessage('E',"SmoothGridReset","cannot reset this topology"); return(1);}

    /* remove node only when necessary */
    if (OBJT(MYVERTEX(theNode))==BVOBJ)
    {
      if (GetMidNodeParam(theNode,&lambda_old)!=0) return(1);
    }
    else
    {
      /* calculate old lambda */
      V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CornerNodes[1])),CVECT(MYVERTEX(CornerNodes[0])),x1);
      V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(theNode)),CVECT(MYVERTEX(CornerNodes[0])),x2);
      lambda_old = x2/x1;
    }
    if (LOCAL_EQUAL(0.5,lambda_old)) continue;

    /* remove mid node */
    if (MoveMidNode(theMG,theNode,0.5)!=0) return(1);
    MoveInfo[1]++;
#else
    PrintErrorMessage('E',"SmoothGridReset","3D not implemented yet");
    return(1);
#endif
  }

  /*    move center nodes of quadrilaterals  */
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    /* skip node if it is a copy from a lower level */
    if (CORNERTYPE(theNode)) continue;

    /* skip node if it is not a center node */
    if (NTYPE(theNode)!=CENTER_NODE) continue;

    /* use only quadrilaterals */
    theVertex=MYVERTEX(theNode);
    fatherElement = VFATHER(theVertex);

#ifdef __TWODIM__
    if (TAG(fatherElement)!=QUADRILATERAL) continue;     /* obsolete */
#else
    PrintErrorMessage('E',"SmoothGridReset","3D not implemented yet");
    return(1);
#endif
    if (OBJT(fatherElement)==BEOBJ)
    {
      if (DefaultBndElemCenterLocal(theNode,defaultLocal)!=0) return(1);
      if (V2_LOCAL_EQUAL(LCVECT(theVertex),defaultLocal)) continue;
      /* calculate new global coordinates */
      CORNER_COORDINATES(fatherElement,coe,CornerPtrs);
      LOCAL_TO_GLOBAL(coe,CornerPtrs,defaultLocal,newPos);
    }
    else
    {
      /* remove center vertices when xi!=0.5 or eta!=0.5 (default in refine) */
      if (V2_LOCAL_EQUAL(LCVECT(theVertex),LocalCenter)) continue;
      /* calculate new global coordinates */
      CORNER_COORDINATES(fatherElement,coe,CornerPtrs);
      LOCAL_TO_GLOBAL(coe,CornerPtrs,LocalCenter,newPos);
    }

    /* re-move center node */
    if (MoveNode(theMG,theNode,newPos)!=0) return(1);
    MoveInfo[0]++;
  }

  return(0);
}

/****************************************************************************/
/*D
   SmoothMultiGrid - Interprete and execute a smooth multigrid command

   SYNOPSIS:
   INT SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  niter - number of iterations to do
   .  bdryFlag - see 'smoothmg' command.

   DESCRIPTION:
   This function smoothes all grid levels of a multigrid hierarchy.
   It processes all grid levels from bottom to top. Within each grid
   level all nodes that are not already in coarser levels are set to
   a new position obtained by the center of gravity of their neighboring
   nodes. The processing from bottom to top level ensures fast convergence
   of the algorithm. Caution! The algorithm may produce undesirable
   results for non-convex domains.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 if error occured.
   D*/
/****************************************************************************/

INT SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag)
{
  INT l,i,n,m,k;
  DOUBLE N;
  GRID *theGrid;
  NODE *node;
  ELEMENT *eptr;
  EDGE *theEdge;
  VERTEX *vptr;
  LINK *lptr;
  DOUBLE *corn[MAX_CORNERS_OF_ELEM],*y,*cvect;
  DOUBLE x[DIM],old_x[DIM];

  if (bdryFlag) {
    PrintErrorMessage('E',"SmoothMultiGrid",
                      "Smoothing boundary nodes not implemented");
    return(GM_ERROR);
  }

  n = niter;
  if (n<=0) n = 1;
  if (n>50) n = 50;

  for (i=0; i<n; i++)
  {
    for (l=0; l<=theMG->topLevel; l++)
    {
      theGrid=GRID_ON_LEVEL(theMG,l);

      /* update global coordinates of new nodes */
      if (l!=0)
        for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node))
          if (!CORNERTYPE(node))
          {
            vptr=MYVERTEX(node);
            if ((OBJT(vptr)!=BVOBJ)||(bdryFlag!=0))
            {
              CORNER_COORDINATES(VFATHER(vptr),m,corn);
              LOCAL_TO_GLOBAL(m,corn,LCVECT(vptr),CVECT(vptr));
            }
          }

      for (node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node))
      {
        /* skip node if it is a copy from a lower level */
        if (l>0)
          if (CORNERTYPE(node))
            continue;
        vptr = MYVERTEX(node);
        /* skip node if it on the boundary */
        if (OBJT(vptr) == BVOBJ)
          continue;
        cvect = CVECT(vptr);
        V_DIM_CLEAR(x);
        N=0.0;
        for (lptr=START(node); lptr!=NULL; lptr=NEXT(lptr))
        {
          y = CVECT(MYVERTEX(NBNODE(lptr)));
          V_DIM_ADD1(y,x);
          N+=1.0;
        }

        V_DIM_SCALESET(1/N,x,cvect);

        /* if there is a father element, change local variables */
        if (l!=0)
        {
          V_DIM_COPY(cvect,old_x);
          eptr = FindFather(vptr);
          if (eptr == NULL)
          {
            PrintErrorMessage('W',"SmoothGrid",
                              "cannot find father element");
            V_DIM_COPY(old_x,cvect);
            return(GM_ERROR);
          }
          else
          {
            CORNER_COORDINATES(eptr,m,corn);
            UG_GlobalToLocal(m,(const DOUBLE **)corn,
                             cvect,LCVECT(vptr));
            for (k=0; k<EDGES_OF_ELEM(eptr); k++) {
              theEdge =
                GetEdge(CORNER(eptr,CORNER_OF_EDGE(eptr,k,0)),
                        CORNER(eptr,CORNER_OF_EDGE(eptr,k,1)));
              if (MIDNODE(theEdge) == node) {
                SETONEDGE(vptr,k);
                break;
              }
            }
            VFATHER(vptr) = eptr;
          }
        }
      }
    }
  }

  return(GM_OK);
}
