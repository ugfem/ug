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
#include "udm.h"

#define SMALL_LOCAL    1.E-4

#define LOCAL_EQUAL(A,V)     (ABS((A)-(V))< (5E-3))
#define IS_VALUE2(A,V)     (ABS((A)-(V))< (SMALL_C))
#define IS_0_OR_1(V)   (LOCAL_EQUAL(V,0) || LOCAL_EQUAL(V,1))
#define V2_LOCAL_EQUAL(A,B) ((ABS((A)[0]-(B)[0])<SMALL_LOCAL)&&(ABS((A)[1]-(B)[1])<SMALL_LOCAL))
#define V3_LOCAL_EQUAL(A,B) ((ABS((A)[0]-(B)[0])<SMALL_LOCAL)&&(ABS((A)[1]-(B)[1])<SMALL_LOCAL)&&(ABS((A)[2]-(B)[2])<SMALL_LOCAL))
#ifdef __TWODIM__
#define V_DIM_LOCAL_EQUAL(A,B)      V2_LOCAL_EQUAL(A,B)
#else
#define V_DIM_LOCAL_EQUAL(A,B)      V3_LOCAL_EQUAL(A,B)
#endif
#define MIN_DIFF 0.05   /* minimum local distance between midnode and cornernode */

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

#define ISTEPS 100
/* calculate lambda as a function of length of boundary edge */
static DOUBLE LambdaOfLengthOfEdge(ELEMENT *theElement, INT edge, DOUBLE LambdaOfLength)
{
  INT n, reverse, step;
  BNDS *bnds;
  DOUBLE_VECTOR BndPoint0, BndPoint1;
  DOUBLE diff, len, len1, len2, lambda[DIM_OF_BND];

  if (OBJT(theElement)!=BEOBJ) return(LambdaOfLength);
  bnds = ELEM_BNDS(theElement,edge);
  if (bnds==NULL) return(LambdaOfLength);    /* nothing to do */

  /* determine orientation of boundary lambda */
  lambda[0] = 0;
  BNDS_Global(bnds,lambda,BndPoint1);
  reverse = TRUE;
  V_DIM_COPY(CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,edge,0)))),BndPoint0);
  if (V_DIM_ISEQUAL(BndPoint0,BndPoint1)) reverse = FALSE;

  /* linear interpolation */
  V_DIM_LINCOMB(1-LambdaOfLength,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,edge,0)))),
                LambdaOfLength,CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_EDGE(theElement,edge,1)))),
                BndPoint0);

  /* boundary interpolation */
  if (reverse==TRUE)
    lambda[0] = 1.-LambdaOfLength;
  else
    lambda[0] = LambdaOfLength;
  BNDS_Global(bnds,lambda,BndPoint1);
  if (reverse==TRUE) printf("reverse: element %d, edge %d \n",ID(theElement),edge);

  V_DIM_EUKLIDNORM_OF_DIFF(BndPoint0,BndPoint1,diff);
  if (diff <= MAX_PAR_DIST) return(LambdaOfLength);    /* boundary is a straight line */

  lambda[0] = 0;
  BNDS_Global(bnds,lambda,BndPoint1);

  /* calculate length */
  len = 0;
  for (n=1; n<=ISTEPS; n++)
  {
    V_DIM_COPY(BndPoint1,BndPoint0);
    lambda[0] = n/(DOUBLE)ISTEPS;
    BNDS_Global(bnds,lambda,BndPoint1);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint1,BndPoint0,diff);
    len = len + diff;
  }
  /* calculate lambda */
  len1 = 0;
  lambda[0] = 0;
  step = 0;
  BNDS_Global(bnds,lambda,BndPoint1);
  for (n=1; n<=ISTEPS; n++)
  {
    V_DIM_COPY(BndPoint1,BndPoint0);
    lambda[0] = n/(DOUBLE)ISTEPS;
    BNDS_Global(bnds,lambda,BndPoint1);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint1,BndPoint0,diff);
    len1 = len1 + diff;
    if (len1/len>=LambdaOfLength) break;
    step = n;
    len2 = len1;
  }

  /* use smaller steps */
  lambda[0] = step/(DOUBLE)ISTEPS;
  BNDS_Global(bnds,lambda,BndPoint1);
  for (n=1; n<=ISTEPS; n++)
  {
    V_DIM_COPY(BndPoint1,BndPoint0);
    lambda[0] = n/(DOUBLE)ISTEPS/(DOUBLE)ISTEPS + step/(DOUBLE)ISTEPS;
    BNDS_Global(bnds,lambda,BndPoint1);
    V_DIM_EUKLIDNORM_OF_DIFF(BndPoint1,BndPoint0,diff);
    len2 = len2 + diff;
    if (len2/len>=LambdaOfLength) break;
  }

  if (reverse==TRUE)
    return(1-lambda[0]);
  else
    return(lambda[0]);
}

/****************************************************************************/
/*
   NewPosCurvedBoundary - modify displacement of center node according to a curved boundary

   SYNOPSIS
   static INT NewPosCurvedBoundary(ELEMENT *theElement, NODE *CenterNode, DOUBLE *MidNodeLambdaNew)

   PARAMETERS
   .  theElement     - father element of the center node
   .  CenterNode     - center node
   . MidNodeLambdaNew  - actual lambda for mid nodes

   DESCRIPTION
   For boundary father elements on a curved boundary it exists one midnode with a 'MOVED' vertex on the
   boundary. In this functions this is taken into account by modifying the displacement of the center node
   calculated before with the function NewPosCenterNode. Additionaly the node opposite to the boundary
   mid node is moved.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

static INT NewPosCurvedBoundary(ELEMENT *theElement, NODE *CenterNode, DOUBLE *MidNodeLambdaNew)
{
  BNDS *bnds;
  DOUBLE_VECTOR bnd_global, bnd_global0, global, CenterDiff, bnd_diff, opp_diff;
  DOUBLE lambda[DIM_OF_BND], *CornerPtrs[MAX_CORNERS_OF_ELEM];
  NODE *theBndNode, *theOppNode;
  EDGE *theEdge;
  INT j, coe, cb0, cb1, co0, co1, moved;
  DOUBLE len_opp, len_bnd, newLambda, *local;

  V_DIM_CLEAR(CenterDiff);
  moved = 0;
  for (j=0; j<EDGES_OF_ELEM(theElement); j++)
  {
    /* find boundary side and get node on boundary edge */
    bnds = ELEM_BNDS(theElement,j);
    if (bnds==NULL) continue;
    cb0 = CORNER_OF_EDGE(theElement,j,0);
    cb1 = CORNER_OF_EDGE(theElement,j,1);
    theEdge=GetEdge(CORNER(theElement,cb0),CORNER(theElement,cb1));
    ASSERT(theEdge != NULL);
    theBndNode = MIDNODE(theEdge);
    if (theBndNode == NULL) continue;
    if (!MOVED(MYVERTEX(theBndNode)) || !OBJT(MYVERTEX(theBndNode))==BVOBJ) continue;

    /* get node on edge opposite to boundary edge */
    co0 = CORNER_OF_EDGE(theElement,OPPOSITE_EDGE(theElement,j),0);
    co1 = CORNER_OF_EDGE(theElement,OPPOSITE_EDGE(theElement,j),1);
    theEdge=GetEdge(CORNER(theElement,co0),CORNER(theElement,co1));
    ASSERT(theEdge != NULL);
    theOppNode = MIDNODE(theEdge);
    if (theOppNode == NULL) continue;
    /* we need the same father element */
    if (VFATHER(MYVERTEX(theBndNode))!=VFATHER(MYVERTEX(theOppNode))) continue;

    /* get lambda of regular refined opposite node and apply it to the boundary node */
    /* only valid for quadrilaterals */
    newLambda = MidNodeLambdaNew[ID(MYVERTEX(theOppNode))];
    lambda[0] = LambdaOfLengthOfEdge(theElement,j,1-newLambda);

    /* calculate global coordinates of (virtual) boundary mid node */
    /* note: (lambda of boundary edge) = 1 - (lambda of opposite edge) */
    BNDS_Global(bnds,lambda,bnd_global);
    /* calculate global coordinates assuming straight boundary */
    V_DIM_LINCOMB(  newLambda,CVECT(MYVERTEX(CORNER(theElement,cb0))),
                    1-newLambda,CVECT(MYVERTEX(CORNER(theElement,cb1))),bnd_global0);

    /* difference between straight and curved boundary at boundary node */
    V_DIM_LINCOMB(1.0,bnd_global,-1.0,bnd_global0,opp_diff);
    /* scale diff according to length of edges */
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CORNER(theElement,cb0))),
                             CVECT(MYVERTEX(CORNER(theElement,cb1))),len_bnd);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CORNER(theElement,co0))),
                             CVECT(MYVERTEX(CORNER(theElement,co1))),len_opp);
    V_DIM_SCALE(0.5*len_opp/len_bnd,opp_diff);    /* half of difference */

    /* current global coordinates of opposite node */
    V_DIM_LINCOMB(1-newLambda,CVECT(MYVERTEX(CORNER(theElement,co0))),
                  newLambda,CVECT(MYVERTEX(CORNER(theElement,co1))),global);
    /* add difference */
    V_DIM_ADD1(opp_diff,global);
    /* update coordinates of opposite node */
    V_DIM_COPY(global,CVECT(MYVERTEX(theOppNode)));
    CORNER_COORDINATES(theElement,coe,CornerPtrs);
    UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,global,LCVECT(MYVERTEX(theOppNode)));

    /* move this node with MoveNode instead with MoveMidNode */
    SETTHEFLAG(theOppNode,1);

    /* local coordinates of boundary node */
    newLambda = MidNodeLambdaNew[ID(MYVERTEX(theBndNode))];
    lambda[0] = LambdaOfLengthOfEdge(theElement,j,newLambda);
    MidNodeLambdaNew[ID(MYVERTEX(theBndNode))] = lambda[0];
    BNDS_Global(bnds,lambda,bnd_global);
    V_DIM_LINCOMB(1-newLambda,CVECT(MYVERTEX(CORNER(theElement,cb0))),
                  newLambda,CVECT(MYVERTEX(CORNER(theElement,cb1))),bnd_global0);
    V_DIM_LINCOMB(1.0,bnd_global,-1.0,bnd_global0,bnd_diff);

    /* calculate displacement of global coordinates of center node w.r.t the
       displacementot the boundary and the opposite node */
    local = LCVECT(MYVERTEX(CenterNode));
    /* only valid for quadrilaterals */
    if (j==0)
    {
      V_DIM_LINCOMB(1.0,CenterDiff,1-local[1],bnd_diff,CenterDiff);
      V_DIM_LINCOMB(1.0,CenterDiff,local[1],opp_diff,CenterDiff);
    }
    else if (j==1)
    {
      V_DIM_LINCOMB(1.0,CenterDiff,local[0],bnd_diff,CenterDiff);
      V_DIM_LINCOMB(1.0,CenterDiff,1-local[0],opp_diff,CenterDiff);
    }
    else if (j==2)
    {
      V_DIM_LINCOMB(1.0,CenterDiff,local[1],bnd_diff,CenterDiff);
      V_DIM_LINCOMB(1.0,CenterDiff,1-local[1],opp_diff,CenterDiff);
    }
    else
    {
      V_DIM_LINCOMB(1.0,CenterDiff,1-local[0],bnd_diff,CenterDiff);
      V_DIM_LINCOMB(1.0,CenterDiff,local[0],opp_diff,CenterDiff);
    }
    moved++;
  }
  /* update center node position */
  if (moved)
  {
    V_DIM_ADD1(CenterDiff,CVECT(MYVERTEX(CenterNode)));
    UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,CVECT(MYVERTEX(CenterNode)),LCVECT(MYVERTEX(CenterNode)));
    SETTHEFLAG(CenterNode,1);
  }
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
   Search on all edges of theElement whether there are 'MOVED' vertices or not.

   RETURN VALUE
   INT
   .n   0: no 'MOVED' vertices found
   .n   1: at least one 'MOVED' vertex found
 */
/****************************************************************************/

static INT MovedNode (ELEMENT *theElement)
{
  INT j;
  NODE *theNode;
  EDGE *theEdge;

  for (j=0; j<EDGES_OF_ELEM(theElement); j++)
  {
    theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),
                    CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
    ASSERT(theEdge != NULL);
    theNode = MIDNODE(theEdge);
    if (theNode == NULL) continue;
    if (MOVED(MYVERTEX(theNode)) && OBJT(MYVERTEX(theNode))==BVOBJ) return(1);
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
   .  newPos            - new position of the center node in global coordinates
   .  newLPos           - new position of the center node in local coordinates

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
                            DOUBLE LimitLocDis, NODE *CenterNode)
{
  DOUBLE *nbCorners[MAX_CORNERS_OF_ELEM],*fatherCorners[MAX_CORNERS_OF_ELEM], localmove;
  INT i,j,k,coe,fcorn,numOfSides,nmove[DIM];
  DOUBLE_VECTOR LocMove, CenterPVertex, sideMid, CenterPoint;
  DOUBLE_VECTOR newLPos;

  V_DIM_SET(0.5,newLPos);

  numOfSides = SIDES_OF_ELEM(fatherElement);
  CORNER_COORDINATES(fatherElement,fcorn,fatherCorners);

  /* calculate old global coordinates of center node according to xi=eta=0.5 */
  LOCAL_TO_GLOBAL(fcorn,fatherCorners,newLPos,CenterPVertex);
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

    /* determine direction in local coordinates */
    /* only valid for quadrilaterals */
    ASSERT(TAG(fatherElement)==QUADRILATERAL);
    if (i==0)
    { LocMove[1] -= localmove; nmove[1]++;}
    else if (i==1)
    { LocMove[0] += localmove; nmove[0]++;}
    else if (i==2)
    { LocMove[1] += localmove; nmove[1]++;}
    else
    { LocMove[0] -= localmove; nmove[0]++;}
  }

  /* displacement in local coordinates according to all neighbour elements */
  for (i=0; i<DIM; i++)
    if (nmove[i]) LocMove[i] = LocMove[i]/nmove[i];

  /* new local coordinates of center node */
  V_DIM_ADD1(LocMove,newLPos);
  /* limit moving distance */
  for (i=0; i<DIM; i++)
  {
    newLPos[i] = MAX(newLPos[i],0.5-LimitLocDis);
    newLPos[i] = MIN(newLPos[i],0.5+LimitLocDis);
  }

  /* update local and global coordinates of center node */
  V_DIM_COPY(newLPos,LCVECT(MYVERTEX(CenterNode)));
  LOCAL_TO_GLOBAL(fcorn,fatherCorners,newLPos,CVECT(MYVERTEX(CenterNode)));

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
                           NODE *cornerNodes[], DOUBLE *lambda)

   PARAMETERS
   .  theElement     - father element of the center node
   .  centerVertex   - vertex of the center node
   .  cornerNodes    - corner nodes on the same side
   .  lambda         - new normalized position of sidemidNode on this side

   DESCRIPTION
   Calculate a new position of the side mid according to the local coordinates of the center node

   ----------------S1
 |                |
 |   theElement   |
 |   (level l-1)  |
 |               OSM
 |                |
 |       CP      NSM
 |                |
   ----------------S0

   CP: center point of theElement;  OSM: old side mid position;
   S1, S0:  corner nodes of side (level=l);

   output: lambda (normalized position (NSM-S0)/(S1-S0) )


   RETURN VALUE
   INT
   .n   0: ok
 */
/**********************************************************************************/

static INT LambdaFromQuad (ELEMENT *theElement,VERTEX *centerVertex,
                           NODE *cornerNodes[], DOUBLE *lambda)
{
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM], *LocalCoord;
  DOUBLE_VECTOR lcorn0, lcorn1;
  INT coe,coord;

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

  return(0);
}
/****************************************************************************/
/*
   DefaultPosCurvedBoundary -  of center node according to a curved boundary

   SYNOPSIS
   static INT DefaultPosCurvedBoundary(ELEMENT *theElement, NODE *CenterNode, INT LambdaMod,
                                    DOUBLE *MidNodeLambdaNew, DOUBLE *MidNodeLambdaOld)

   PARAMETERS
   .  theElement     - father element of the center node
   .  CenterNode     - center node of father element
   .  LambdaMod      - TRUE for 'smoothgrid $b'
   .  MidNodeLambdaNew - new lambda of mid nodes
   .  MidNodeLambdaOld - old lambda of mid nodes

   DESCRIPTION
   Reset the center node and the node on the edge opposite to the boundary to
   their default positions.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

static INT DefaultPosCurvedBoundary(ELEMENT *theElement, NODE *CenterNode, INT LambdaMod,
                                    DOUBLE *MidNodeLambdaNew, DOUBLE *MidNodeLambdaOld)
{
  BNDS *bnds;
  DOUBLE_VECTOR bnd_global, bnd_global0, diff,  global, CenterDiff;
  DOUBLE lambda[DIM_OF_BND], *CornerPtrs[MAX_CORNERS_OF_ELEM], lambda_old;
  NODE *theBndNode, *theOppNode;
  EDGE *theEdge;
  INT j, coe, cb0, cb1, co0, co1, moved;
  DOUBLE len_opp, len_bnd;

  V_DIM_CLEAR(CenterDiff);
  moved = 0;
  for (j=0; j<EDGES_OF_ELEM(theElement); j++)
  {
    /* find boundary side and get node on boundary edge */
    bnds = ELEM_BNDS(theElement,j);
    if (bnds==NULL) continue;
    cb0 = CORNER_OF_EDGE(theElement,j,0);
    cb1 = CORNER_OF_EDGE(theElement,j,1);
    theEdge=GetEdge(CORNER(theElement,cb0),CORNER(theElement,cb1));
    ASSERT(theEdge != NULL);
    theBndNode = MIDNODE(theEdge);
    if (theBndNode == NULL) continue;
    if (!MOVED(MYVERTEX(theBndNode)) || !OBJT(MYVERTEX(theBndNode))==BVOBJ) continue;

    /* get node on edge opposite to boundary edge */
    co0 = CORNER_OF_EDGE(theElement,OPPOSITE_EDGE(theElement,j),0);
    co1 = CORNER_OF_EDGE(theElement,OPPOSITE_EDGE(theElement,j),1);
    theEdge=GetEdge(CORNER(theElement,co0),CORNER(theElement,co1));
    ASSERT(theEdge != NULL);
    theOppNode = MIDNODE(theEdge);
    if (theOppNode == NULL) continue;
    /* we need the same father element */
    if (VFATHER(MYVERTEX(theBndNode))!=VFATHER(MYVERTEX(theOppNode))) continue;

    /* calculate global coordinates of a regular refined boundary mid node */
    if (LambdaMod==TRUE)
    {     /* get lambda of bnd-node and reset center node (only used for smoothgrid $b) */
      if (GetMidNodeParam(theBndNode,&lambda_old)!=0) continue;
      MidNodeLambdaOld[ID(MYVERTEX(theBndNode))] = lambda_old;
      lambda[0] = LambdaOfLengthOfEdge(theElement,j,0.5);
      MidNodeLambdaNew[ID(MYVERTEX(theBndNode))] = lambda[0];
      /* default local coordinates of center node */
      V_DIM_SET(0.5,LCVECT(MYVERTEX(CenterNode)));
      /* calculate new global coordinates */
      CORNER_COORDINATES(theElement,coe,CornerPtrs);
      LOCAL_TO_GLOBAL(coe,CornerPtrs,LCVECT(MYVERTEX(CenterNode)),CVECT(MYVERTEX(CenterNode)));
    }
    else
      lambda[0] = 0.5;
    BNDS_Global(bnds,lambda,bnd_global);
    /* calculate global coordinates assuming straight boundary */
    V_DIM_LINCOMB(0.5,CVECT(MYVERTEX(CORNER(theElement,cb0))),
                  0.5,CVECT(MYVERTEX(CORNER(theElement,cb1))),bnd_global0);

    /* difference between straight and curved boundary at boundary node */
    V_DIM_LINCOMB(1.0,bnd_global,-1.0,bnd_global0,diff);
    /* add diff of boundary node to center node */
    V_DIM_SCALEADD1(0.5,diff,CenterDiff)
    /* scale diff according to length of edges */
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CORNER(theElement,cb0))),
                             CVECT(MYVERTEX(CORNER(theElement,cb1))),len_bnd);
    V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(CORNER(theElement,co0))),
                             CVECT(MYVERTEX(CORNER(theElement,co1))),len_opp);
    V_DIM_SCALE(0.5*len_opp/len_bnd,diff);    /* half of difference */

    /* global coordinates of regular refined opposite node */
    V_DIM_LINCOMB(0.5,CVECT(MYVERTEX(CORNER(theElement,co0))),
                  0.5,CVECT(MYVERTEX(CORNER(theElement,co1))),global);
    /* add difference */
    V_DIM_ADD1(diff,global);
    /* update coordinates of opposite node */
    V_DIM_COPY(global,CVECT(MYVERTEX(theOppNode)));
    CORNER_COORDINATES(theElement,coe,CornerPtrs);
    UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,global,LCVECT(MYVERTEX(theOppNode)));

    /* move this node with MoveNode instead with MoveMidNode */
    SETTHEFLAG(theOppNode,1);

    /* add diff of opposite node to center node */
    V_DIM_SCALEADD1(0.5,diff,CenterDiff)
    moved++;
  }
  /* update center node position */
  if (moved)
  {
    V_DIM_SCALE(1./moved,CenterDiff);
    V_DIM_ADD1(CenterDiff,CVECT(MYVERTEX(CenterNode)));
    UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,CVECT(MYVERTEX(CenterNode)),LCVECT(MYVERTEX(CenterNode)));
    SETTHEFLAG(CenterNode,1);
  }
  return(0);
}

/* calculate position of bnd-midnode in order to get the line bnd-midnode -> center node
   orthogonal to the boundary */
static INT LambdaOrthoBnd2D(const ELEMENT *fatherElement, const INT edge, const NODE *MidNode,
                            DOUBLE *lambda)
{
  VERTEX *theVertex;
  NODE *CenterNode;
  DOUBLE_VECTOR BndPoint0,BndPoint1,midBndPoint,MPVec,TangVec,NormVec;
  DOUBLE Lambda0, midLambda, Lambda1,bndLambda[DIM-1],veclen ;
  DOUBLE *CenterPoint;
  BNDS *bnds;
  EDGE *theEdge;
  LINK *theLink;
  DOUBLE min_fac,area,min_diff;
  INT maxiter,iter;

  /* is midnode on specified edge ? */
  theEdge=GetEdge(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,edge,0)),
                  CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,edge,1)));
  if (theEdge==NULL) return(1);
  if (MidNode!=MIDNODE(theEdge)) return(1);

  /* find center node */
  for (theLink=START(MidNode); theLink!=0; theLink=NEXT(theLink))
    if (NTYPE(NBNODE(theLink))==CENTER_NODE)
    {
      CenterNode = NBNODE(theLink);
      break;
    }

  CenterPoint = CVECT(MYVERTEX(CenterNode));
  min_fac = 0.00001;
  maxiter = 50;

  theVertex = MYVERTEX(MidNode);
  if (OBJT(theVertex) != BVOBJ) return(GM_ERROR);

  bnds = ELEM_BNDS(fatherElement,edge);
  min_diff = min_fac;

  Lambda0 = 0.;
  Lambda1 = 1.;
  iter=0;

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
    NormVec[0] = - TangVec[1]; NormVec[1] = TangVec[0];

    /* calculate direction MidNode -> CenterNode  */
    V_DIM_SUBTRACT(CenterPoint,midBndPoint,MPVec);
    V_DIM_EUKLIDNORM(MPVec,veclen);
    V_DIM_SCALE(1./veclen,MPVec);

    /* calculate vector product */
    V2_VECTOR_PRODUCT(MPVec,NormVec,area);
    if (area==0.) break;
    if (area>0)
      Lambda0 = midLambda;
    else
      Lambda1 = midLambda;
  }
  *lambda = midLambda;

  if (iter > 40)
    printf ("iter %d, midnode-id %ld, midlambda %f, lambda %f,  area %lf\n",
            iter,ID(MidNode),midLambda,*lambda,area);
  return(0);
}
/*******************************************************************/
/*                                                                 */
/*                                                                 */
/*   OppNodes[1]           OppNode                OppNodes[0]      */
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
/*    BndNodes[0]           BndNode             BndNodes[1]        */
/*                                                                 */
/*******************************************************************/
static INT EqualDistBndElem(ELEMENT *fatherElement, INT bnd_side, NODE *CenterNode, DOUBLE *MidNodeLambdaNew)
{
  NODE *BndNodes[2], *BndNode, *OppNodes[2], *OppNode, *SideNodes[2];
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM];
  DOUBLE dist0, dist1, CenterNodeDist, OppNodeDist,  bndLambda[DIM_OF_BND], lambda;
  DOUBLE_VECTOR BndPoint, testLocal;
  DOUBLE NormVec[DIM], VecLen, *CenterLocal, *OppLocal;
  EDGE *theEdge;
  ELEMENT *fatherElement2;
  INT c0,c1, coe;
  BNDS *bnds;


  if (TAG(fatherElement)!=QUADRILATERAL) return(1);
  if (OBJT(fatherElement)!=BEOBJ) return(1);

  /* mid node on boundary side */
  c0 = CORNER_OF_EDGE(fatherElement,bnd_side,0);
  c1 = CORNER_OF_EDGE(fatherElement,bnd_side,1);
  theEdge=GetEdge(CORNER(fatherElement,c0),CORNER(fatherElement,c1));
  if (theEdge==NULL) return(1);
  BndNode = MIDNODE(theEdge);
  if (BndNode==NULL) return(1);
  /* corner nodes on boundary */
  BndNodes[0] = CORNER(fatherElement,c0);
  BndNodes[1] = CORNER(fatherElement,c1);

  /* opposite edge */
  c0 = CORNER_OF_EDGE(fatherElement,OPPOSITE_EDGE(fatherElement,bnd_side),0);
  c1 = CORNER_OF_EDGE(fatherElement,OPPOSITE_EDGE(fatherElement,bnd_side),1);
  theEdge=GetEdge(CORNER(fatherElement,c0),CORNER(fatherElement,c1));
  if (theEdge==NULL) return(1);
  OppNode = MIDNODE(theEdge);
  if (OppNode == NULL) return(1);
  OppNodes[0] = CORNER(fatherElement,c0);
  OppNodes[1] = CORNER(fatherElement,c1);

  /* side nodes */
  coe = EDGES_OF_ELEM(fatherElement);
  c0 = CORNER_OF_EDGE(fatherElement,(bnd_side+1)%coe,0);
  c1 = CORNER_OF_EDGE(fatherElement,(bnd_side+1)%coe,1);
  theEdge=GetEdge(CORNER(fatherElement,c0),CORNER(fatherElement,c1));
  if (theEdge==NULL) return(1);
  SideNodes[0] = MIDNODE(theEdge);
  if (SideNodes[0]==NULL) return(1);
  c0 = CORNER_OF_EDGE(fatherElement,(bnd_side+3)%coe,0);
  c1 = CORNER_OF_EDGE(fatherElement,(bnd_side+3)%coe,1);
  theEdge=GetEdge(CORNER(fatherElement,c0),CORNER(fatherElement,c1));
  if (theEdge==NULL) return(1);
  SideNodes[1] = MIDNODE(theEdge);
  if (SideNodes[1]==NULL) return(1);

  /* get new lambda of BndNode (determined in LambdaOrthoBnd2D) */
  bndLambda[0]= MidNodeLambdaNew[ID(MYVERTEX(BndNode))];
  /* get global coordinates of BndNode */
  bnds = ELEM_BNDS(fatherElement,bnd_side);
  BNDS_Global(bnds,bndLambda,BndPoint);

  /* calculate new boundary distance of OppNode */
  V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(OppNodes[0])),CVECT(MYVERTEX(BndNodes[1])),dist0);
  V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(OppNodes[1])),CVECT(MYVERTEX(BndNodes[0])),dist1);
  lambda = MidNodeLambdaNew[ID(MYVERTEX(OppNode))];
  if (VFATHER(MYVERTEX(OppNode))!=fatherElement) lambda = 1.-lambda;
  OppNodeDist = (1.-lambda)*dist0 + lambda*dist1;

  /* calculate new boundary distance of CenterNode */
  lambda = MidNodeLambdaNew[ID(MYVERTEX(SideNodes[0]))];
  if (VFATHER(MYVERTEX(SideNodes[0]))!=fatherElement) lambda = 1.-lambda;
  dist0 = lambda*dist0;
  lambda = MidNodeLambdaNew[ID(MYVERTEX(SideNodes[1]))];
  if (VFATHER(MYVERTEX(SideNodes[1]))==fatherElement) lambda = 1.-lambda;
  dist1 = lambda*dist1;
  CenterLocal = LCVECT(MYVERTEX(CenterNode));
  /* only valid for quadrilaterals */
  if (bnd_side==0)
    lambda = CenterLocal[0];
  else if (bnd_side==1)
    lambda = CenterLocal[1];
  else if (bnd_side==2)
    lambda = 1.-CenterLocal[0];
  else
    lambda = 1.-CenterLocal[1];
  CenterNodeDist = lambda*dist0 + (1.-lambda)*dist1;

  /* normal vector on wall (already determined in LambdaOrthoBnd2D) */
  V_DIM_SUBTRACT(CVECT(MYVERTEX(CenterNode)),BndPoint,NormVec);
  V_DIM_EUKLIDNORM(NormVec,VecLen);

  /* new OppNode position */
  V_DIM_LINCOMB(1.,BndPoint,OppNodeDist/VecLen,NormVec,CVECT(MYVERTEX(OppNode)));
  /* calculate new local coordinates of OppNode */
  fatherElement2 = VFATHER(MYVERTEX(OppNode));
  CORNER_COORDINATES(fatherElement2,coe,CornerPtrs);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,CVECT(MYVERTEX(OppNode)),
                   LCVECT(MYVERTEX(OppNode)));

  /* apply limits for position of OppNode */
  OppLocal = LCVECT(MYVERTEX(OppNode));
  V_DIM_COPY(OppLocal,testLocal);
  if (bnd_side==0 || bnd_side==2)
    testLocal[0] = MAX(MIN(1-MIN_DIFF,testLocal[0]),MIN_DIFF);
  else
    testLocal[1] = MAX(MIN(1-MIN_DIFF,testLocal[1]),MIN_DIFF);
  if (!V_DIM_ISEQUAL(OppLocal,testLocal))
  {
    V_DIM_COPY(testLocal,OppLocal);
    LOCAL_TO_GLOBAL(coe,CornerPtrs,OppLocal,CVECT(MYVERTEX(OppNode)));
    V_DIM_SUBTRACT(CVECT(MYVERTEX(OppNode)),BndPoint,NormVec);
    V_DIM_EUKLIDNORM(NormVec,VecLen);
  }

  /* new center node position */
  V_DIM_LINCOMB(1.,BndPoint,CenterNodeDist/VecLen,NormVec,CVECT(MYVERTEX(CenterNode)));
  /* calculate new local coordinates of center node */
  CORNER_COORDINATES(fatherElement,coe,CornerPtrs);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,CVECT(MYVERTEX(CenterNode)),
                   LCVECT(MYVERTEX(CenterNode)));

  SETTHEFLAG(OppNode,1);

  return(0);
}

/*******************************************************************/
/*        fatherElement                                            */
/*        -------------                                            */
/*   OppNodes[1]           OppNode                OppNodes[0]      */
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
static INT EqualDistInnerElem(ELEMENT *fatherElement,  ELEMENT *bndElement, INT bnd_side,
                              NODE *CenterNode, DOUBLE *MidNodeLambdaNew)
{
  NODE *BndNodes[2], *BndNode, *OppNodes[2], *OppNode, *SideNodes[2];
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM];
  DOUBLE dist0, dist1, CenterNodeDist, OppNodeDist,  bndLambda[DIM_OF_BND], lambda;
  DOUBLE_VECTOR BndPoint, SidePoint;
  DOUBLE NormVec[DIM], VecLen, *CenterLocal;
  EDGE *theEdge;
  ELEMENT *fatherElement2;
  INT c0,c1, coe;
  BNDS *bnds;


  if (TAG(fatherElement)!=QUADRILATERAL) return(1);
  if (TAG(bndElement)!=QUADRILATERAL) return(1);
  if (OBJT(bndElement)!=BEOBJ) return(1);

  /* mid node on boundary side */
  c0 = CORNER_OF_EDGE(bndElement,bnd_side,0);
  c1 = CORNER_OF_EDGE(bndElement,bnd_side,1);
  theEdge=GetEdge(CORNER(bndElement,c0),CORNER(bndElement,c1));
  if (theEdge==NULL) return(1);
  BndNode = MIDNODE(theEdge);
  if (BndNode==NULL) return(1);
  /* corner nodes on boundary */
  BndNodes[0] = CORNER(bndElement,c0);
  BndNodes[1] = CORNER(bndElement,c1);

  /* opposite edge on father element  */
  c0 = CORNER_OF_EDGE(fatherElement,OPPOSITE_EDGE(fatherElement,bnd_side),0);
  c1 = CORNER_OF_EDGE(fatherElement,OPPOSITE_EDGE(fatherElement,bnd_side),1);
  theEdge=GetEdge(CORNER(fatherElement,c0),CORNER(fatherElement,c1));
  if (theEdge==NULL) return(1);
  OppNode = MIDNODE(theEdge);
  if (OppNode == NULL) return(1);
  OppNodes[0] = CORNER(fatherElement,c0);
  OppNodes[1] = CORNER(fatherElement,c1);

  /* side nodes on father element */
  coe = EDGES_OF_ELEM(fatherElement);
  c0 = CORNER_OF_EDGE(fatherElement,(bnd_side+1)%coe,0);
  c1 = CORNER_OF_EDGE(fatherElement,(bnd_side+1)%coe,1);
  theEdge=GetEdge(CORNER(fatherElement,c0),CORNER(fatherElement,c1));
  if (theEdge==NULL) return(1);
  SideNodes[0] = MIDNODE(theEdge);
  if (SideNodes[0]==NULL) return(1);
  c0 = CORNER_OF_EDGE(fatherElement,(bnd_side+3)%coe,0);
  c1 = CORNER_OF_EDGE(fatherElement,(bnd_side+3)%coe,1);
  theEdge=GetEdge(CORNER(fatherElement,c0),CORNER(fatherElement,c1));
  if (theEdge==NULL) return(1);
  SideNodes[1] = MIDNODE(theEdge);
  if (SideNodes[1]==NULL) return(1);

  /* get new lambda of BndNode (determined in LambdaOrthoBnd2D) */
  bndLambda[0]= MidNodeLambdaNew[ID(MYVERTEX(BndNode))];
  /* get global coordinates of BndNode */
  bnds = ELEM_BNDS(bndElement,bnd_side);
  BNDS_Global(bnds,bndLambda,BndPoint);

  /* calculate new boundary distance of OppNode */
  V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(OppNodes[0])),CVECT(MYVERTEX(BndNodes[1])),dist0);
  V_DIM_EUKLIDNORM_OF_DIFF(CVECT(MYVERTEX(OppNodes[1])),CVECT(MYVERTEX(BndNodes[0])),dist1);
  lambda = MidNodeLambdaNew[ID(MYVERTEX(OppNode))];
  if (VFATHER(MYVERTEX(OppNode))!=fatherElement) lambda = 1.-lambda;
  OppNodeDist = (1.-lambda)*dist0 + lambda*dist1;

  /* calculate new boundary distance of CenterNode */
  lambda = MidNodeLambdaNew[ID(MYVERTEX(SideNodes[0]))];
  if (VFATHER(MYVERTEX(SideNodes[0]))!=fatherElement) lambda = 1.-lambda;
  V_DIM_LINCOMB(1.-lambda,
                CVECT(MYVERTEX(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,(bnd_side+1)%coe,0)))),
                lambda,CVECT(MYVERTEX(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,(bnd_side+1)%coe,1)))),SidePoint);
  V_DIM_EUKLIDNORM_OF_DIFF(SidePoint,CVECT(MYVERTEX(BndNodes[1])),dist0);
  lambda = MidNodeLambdaNew[ID(MYVERTEX(SideNodes[1]))];
  if (VFATHER(MYVERTEX(SideNodes[1]))!=fatherElement) lambda = 1.-lambda;
  V_DIM_LINCOMB(1.-lambda,
                CVECT(MYVERTEX(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,(bnd_side+3)%coe,0)))),
                lambda,CVECT(MYVERTEX(CORNER(fatherElement,CORNER_OF_EDGE(fatherElement,(bnd_side+3)%coe,1)))),SidePoint);
  V_DIM_EUKLIDNORM_OF_DIFF(SidePoint,CVECT(MYVERTEX(BndNodes[0])),dist1);
  CenterLocal = LCVECT(MYVERTEX(CenterNode));
  /* only valid for quadrilaterals */
  if (bnd_side==0)
    lambda = CenterLocal[0];
  else if (bnd_side==1)
    lambda = CenterLocal[1];
  else if (bnd_side==2)
    lambda = 1.-CenterLocal[0];
  else
    lambda = 1.-CenterLocal[1];
  CenterNodeDist = lambda*dist0 + (1.-lambda)*dist1;

  /* normal vector on wall (this is only an approximation) */
  V_DIM_SUBTRACT(CVECT(MYVERTEX(CenterNode)),BndPoint,NormVec);
  V_DIM_EUKLIDNORM(NormVec,VecLen);

  /* new center node position */
  V_DIM_LINCOMB(1.,BndPoint,CenterNodeDist/VecLen,NormVec,CVECT(MYVERTEX(CenterNode)));
  /* calculate new local coordinates of center node */
  CORNER_COORDINATES(fatherElement,coe,CornerPtrs);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,CVECT(MYVERTEX(CenterNode)),
                   LCVECT(MYVERTEX(CenterNode)));

  /* new OppNode position */
  V_DIM_LINCOMB(1.,BndPoint,OppNodeDist/VecLen,NormVec,CVECT(MYVERTEX(OppNode)));
  /* calculate new local coordinates of OppNode */
  fatherElement2 = VFATHER(MYVERTEX(OppNode));
  CORNER_COORDINATES(fatherElement2,coe,CornerPtrs);
  UG_GlobalToLocal(coe,(const DOUBLE **)CornerPtrs,CVECT(MYVERTEX(OppNode)),
                   LCVECT(MYVERTEX(OppNode)));
  SETTHEFLAG(OppNode,1);

  return(0);
}

/* calculate lambda of mid node according to the neighboring center node(s) */
static INT NewLambdaMidNode(NODE *theNode, DOUBLE *lambda)
{
  VERTEX *theVertex;
  ELEMENT *fatherElement, *oppositeElement;
  NODE *CornerNodes[2],*CenterNodes[2],*node0,*node1;
  INT i,edge,co0,co1,ceN,nlinks,Eside;
  LINK *theLink;
  DOUBLE lambda0, lambda1;

  /* find the two corner nodes on the same edge and the elements on this edge */
  theVertex = MYVERTEX(theNode);
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
    if (NTYPE(NBNODE(theLink))==CENTER_NODE)
    {
      CenterNodes[ceN] = NBNODE(theLink);
      ceN++;
    }
    nlinks++;
  }

  /* 5 possibilities of element neighbourship relationships */

  /* quadrilateral is boundary element */
  if (nlinks==3 && ceN==1)
  {
    fatherElement = VFATHER(MYVERTEX(CenterNodes[0]));
    if (LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),CornerNodes,lambda)!=0)
      return(0);
  }
  /* triangle is boundary element */
  else if (nlinks==4 && ceN==0)
  {
    fatherElement = VFATHER(MYVERTEX(theNode));
    if (LambdaFromTriangle(fatherElement,CornerNodes,lambda)!=0) return(0);
  }
  /* two quadrilaterals */
  else if (nlinks==4 && ceN==2)
  {
    fatherElement = VFATHER(MYVERTEX(CenterNodes[0]));
    oppositeElement = VFATHER(MYVERTEX(CenterNodes[1]));
    if(LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),CornerNodes,&lambda0)!=0)
      return(0);
    if (LambdaFromQuad(oppositeElement,MYVERTEX(CenterNodes[1]),CornerNodes,&lambda1)!=0)
      return(0);
    *lambda = 0.5*(lambda0+lambda1);
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
      if (LambdaFromTriangle(fatherElement,CornerNodes,&lambda0)!=0) return(0);
    }
    oppositeElement = NBELEM(fatherElement,Eside);
    if (TAG(oppositeElement)!=TRIANGLE)
    {
      lambda1 = 0.5;
    }
    else
    {
      if (LambdaFromTriangle(oppositeElement,CornerNodes,&lambda1)!=0) return(0);
    }
    *lambda = 0.5*(lambda0+lambda1);
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
    if(LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),CornerNodes,&lambda0)!=0)
      return(0);
    oppositeElement = NBELEM(fatherElement,Eside);
    if (TAG(oppositeElement)!=TRIANGLE)
    {
      lambda1 = 0.5;
    }
    else
    {
      if (LambdaFromTriangle(oppositeElement,CornerNodes,&lambda1)!=0) return(0);
    }
    *lambda = 0.5*(lambda0+lambda1);
  }
  /* one quadrilateral and one irregular refined quadrilateral */
  else if (ceN==1)
  {
    /* quadrilateral */
    fatherElement = VFATHER(MYVERTEX(CenterNodes[0]));
    if(LambdaFromQuad(fatherElement,MYVERTEX(CenterNodes[0]),CornerNodes,&lambda0)!=0)
      return(0);
    *lambda = 0.5*(lambda0 + 0.5);
  }

  return(0);
}
/* update local and global Coordinates and move nodes to their new position */
static INT MoveNodesOnGrid (GRID *theGrid, DOUBLE_VECTOR *VertexCoord, DOUBLE_VECTOR *VertexLCoord,
                            DOUBLE *MidNodeLambdaOld, DOUBLE *MidNodeLambdaNew, DOUBLE LimitLocDis)
{
  NODE *theNode;
  MULTIGRID *theMG;
  VERTEX *theVertex;
  DOUBLE *VertexLocal, *VertexGlobal, lambda_new, lambda_old, *x[MAX_CORNERS_OF_ELEM];
  DOUBLE_VECTOR oldPos,oldLPos, newPos, newLPos;
  INT i,n,MoveInfo[4];

  theMG = MYMG(theGrid);
  for (i=0; i<4; i++) MoveInfo[i] = 0;
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    /* skip node if it is a copy from a lower level */
    if (CORNERTYPE(theNode)) continue;
    theVertex = MYVERTEX(theNode);

    /* reset local and global variables */
    VertexLocal = LCVECT(MYVERTEX(theNode));
    VertexGlobal = CVECT(MYVERTEX(theNode));
    /* old coordinates */
    V_DIM_COPY(VertexCoord[ID(MYVERTEX(theNode))],oldPos);
    V_DIM_COPY(VertexLCoord[ID(MYVERTEX(theNode))],oldLPos);
    /* save new coordinates */
    V_DIM_COPY(VertexLocal,newLPos);
    V_DIM_COPY(VertexGlobal,newPos);
    /* copy old coordinates on vertices */
    V_DIM_COPY(oldLPos,VertexLocal);
    V_DIM_COPY(oldPos,VertexGlobal);

    if (THEFLAG(theNode)==1)
    {
      if (!V2_LOCAL_EQUAL(newLPos,oldLPos))
      {
        if (MoveNode(theMG,theNode,newPos,FALSE)!=0) return(1);
        SETMOVED(theVertex,1);
        if (NTYPE(theNode)==CENTER_NODE) MoveInfo[0]++;
        if (NTYPE(theNode)==MID_NODE) MoveInfo[1]++;
      }
    }
    else if (NTYPE(theNode)==CENTER_NODE)
    {
      if (!V2_LOCAL_EQUAL(newLPos,oldLPos))
      {
        if (MoveNode(theMG,theNode,newPos,FALSE)!=0) return(1);
        SETMOVED(theVertex,1);
        MoveInfo[0]++;
        /* are there nodes which reached the moving limit ? */
        if ( (LOCAL_EQUAL(XI(theVertex),(0.5+LimitLocDis)) ||
              LOCAL_EQUAL(XI(theVertex),(0.5-LimitLocDis))) ||
             (LOCAL_EQUAL(ETA(theVertex),(0.5+LimitLocDis)) ||
              LOCAL_EQUAL(ETA(theVertex),(0.5-LimitLocDis))))
        {
          /*                     printf("center-node %ld, father-elem%ld reached limit\n",ID(theNode),ID(VFATHER(theVertex))); */
          MoveInfo[2]++;
        }

      }
    }
    else if (NTYPE(theNode)==MID_NODE)
    {
      lambda_new = MidNodeLambdaNew[ID(theVertex)];
      lambda_old = MidNodeLambdaOld[ID(theVertex)];

      if (!LOCAL_EQUAL(lambda_new,lambda_old))
      {
        if (MoveMidNode(theMG,theNode,lambda_new,FALSE)!=0) return(1);
        SETMOVED(theVertex,1);
        MoveInfo[1]++;
      }
      if (LOCAL_EQUAL(lambda_new,(0.5+LimitLocDis)) || LOCAL_EQUAL(lambda_new,(0.5-LimitLocDis)) )
      {
        /*                 printf("mid-node %ld, father-elem%ld reached limit\n",ID(theNode),ID(VFATHER(theVertex))); */
        MoveInfo[3]++;
      }

    }
  }
  /* update nodes on higher level */
  for(i=GLEVEL(theGrid)+1; i<=TOPLEVEL(theMG); i++)
    for (theVertex=FIRSTVERTEX(GRID_ON_LEVEL(theMG,i));
         theVertex!=NULL; theVertex=SUCCV(theVertex))
      if ((OBJT(theVertex) != BVOBJ)) {
        CORNER_COORDINATES(VFATHER(theVertex),n,x);
        LOCAL_TO_GLOBAL(n,x,LCVECT(theVertex),CVECT(theVertex));
      }
      else
        MoveBndMidNode(theMG,theVertex);

  UserWriteF(" %d center nodes and %d mid nodes moved on level %d \n",MoveInfo[0],MoveInfo[1],
             GLEVEL(theGrid));
  if (MoveInfo[2]!=0 || MoveInfo[3]!=0)
  {
    UserWriteF("%d center nodes and %d mid nodes reached limit on level %d\n",MoveInfo[2],
               MoveInfo[3],GLEVEL(theGrid));
  }
  return(0);

}

/* TODO: Get boundary-id without boundary condition functions */
/* !!! this routine uses information specified in the boundary condition functions !!! */
static INT ElemOnSpecBnd(ELEMENT *theElement, const INT *bnd, const INT bnd_num, INT *side)
{
  DOUBLE dummyValues[10];
  static DOUBLE local[] = {0.5};
  INT i,j,type[10];


  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    if (!SIDE_ON_BND(theElement,i) || INNER_BOUNDARY(theElement,i)) continue;
    SideBndCond(theElement,i,local,dummyValues,type);
    for (j=0; j<bnd_num; j++)
      if (type[1]==bnd[j])
      {
        *side = i;
        return(TRUE);
      }
  }
  return(FALSE);
}
/****************************************************************************/
/*
   SmoothGrid - resize all quadrilaterals and triangles on specified grids according
             to the element sizes of grid (l-1)

   SYNOPSIS
   static INT SmoothGrid (MULTIGRID *theMG, INT fl, INT tl, const DOUBLE LimitLocDis,
                       const INT bnd_num, const INT *bnd, const INT option)

   PARAMETERS
   .  theMG        - the Multigrid
   .  fl, tl       - resize elements on grids  fl <= level <= tl
   .  LimitLocDis  - maximum displacement of the nodes in local coordinates  (0<LimitLocDis<0.5)
   .  bnd_num      - number of boundaries with special treatment
   .  bnd[]        - apply ortho-option on this boundaries
   .  option       - option type

   DESCRIPTION
   Resize all quadrilaterals and triangles on grid (l) according to the element sizes of grid (l-1).
   In the first node-loop all center nodes will be modified according to the size of the neighbour elements
   (grid l-1). In the second node-loop all mid nodes will be modified according to the position of the
   (moved) center nodes and the size of the triangles. According to the option type, nodes on boundary elements
   will be moved additionaly. In the last loop the calculated new positions of all nodes will be updated using
   the functions MoveNode and MoveMidNode.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

INT SmoothGrid (MULTIGRID *theMG, INT fl, INT tl, const DOUBLE LimitLocDis,
                const INT bnd_num, const INT *bnd, const INT option)
{
  GRID *theGrid;
  ELEMENT *fatherElement, *nbElement[MAX_SIDES_OF_ELEM], *Level0Father, *boundaryFather;
  VERTEX *theVertex;
  NODE *theNode;
  INT i, numOfSides, side, lev;
  DOUBLE_VECTOR *VertexCoord, *VertexLCoord;
  DOUBLE lambda, lambda_old, *MidNodeLambdaOld, *MidNodeLambdaNew;
  INT MarkKey;

  /* allocate temporary memory */
  MarkTmpMem(MGHEAP(theMG),&MarkKey);
  VertexCoord = (DOUBLE_VECTOR *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE_VECTOR),MarkKey);
  VertexLCoord = (DOUBLE_VECTOR *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE_VECTOR),MarkKey);
  MidNodeLambdaOld = (DOUBLE *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE),MarkKey);
  MidNodeLambdaNew = (DOUBLE *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE),MarkKey);

  for (lev=fl; lev<=tl; lev++)
  {
    theGrid=GRID_ON_LEVEL(theMG,lev);
    /* copy current global and local vertex coordinates */
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
    {
      V_DIM_COPY(CVECT(theVertex),VertexCoord[ID(theVertex)]);
      V_DIM_COPY(LCVECT(theVertex),VertexLCoord[ID(theVertex)]);
    }
    /* default values for lambda */
    memset(MidNodeLambdaOld,0.5,VIDCNT(theMG)*sizeof(DOUBLE));
    memset(MidNodeLambdaNew,0.5,VIDCNT(theMG)*sizeof(DOUBLE));

    /* check node flags */
    for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
      if (THEFLAG(theNode))
      {
        PrintErrorMessage('E',"SmoothGrid","node flag already set");
        ReleaseTmpMem(MGHEAP(theMG),MarkKey);
        return(1);
      }

    if (option==3)
      goto option_b;

    /*    move center nodes of quadrilaterals  */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* skip node if it is a copy from a lower level */
      if (CORNERTYPE(theNode)) continue;

      /* skip node if it is not a center node */
      if (NTYPE(theNode)!=CENTER_NODE) continue;

      /* get father element */
      theVertex=MYVERTEX(theNode);
      fatherElement = VFATHER(theVertex);

      /* create list of neighbour elements */
      numOfSides = SIDES_OF_ELEM(fatherElement);
      for (i=0; i<numOfSides; i++)
        nbElement[i] = NBELEM(fatherElement,i);

      /* calculate new global and local coordinates of center node*/
      if (NewPosCenterNode(fatherElement,nbElement,LimitLocDis,theNode)!=0) goto exit;
    }

    /* move mid nodes  */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* skip node if it is a copy from a lower level */
      if (CORNERTYPE(theNode)) continue;

      /* skip node if it is not a mid node (mid point of an edge) */
      if (NTYPE(theNode)!=MID_NODE) continue;

      theVertex = MYVERTEX(theNode);

      /* calculate lambda of mid node according to the neighboring center node(s) */
      if (NewLambdaMidNode(theNode,&lambda)) continue;

      /* apply limits for lambda  */
      lambda = MAX(lambda,0.5-LimitLocDis);
      lambda = MIN(lambda,0.5+LimitLocDis);

      /* get old lambda */
      if (GetMidNodeParam(theNode,&lambda_old)!=0) goto exit;

      MidNodeLambdaOld[ID(MYVERTEX(theNode))] = lambda_old;
      MidNodeLambdaNew[ID(MYVERTEX(theNode))] = lambda;
    }

option_b:
    /* special treatment for center and mid nodes of boundary elements on a curved boundary */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* skip node if it is a copy from a lower level */
      if (CORNERTYPE(theNode)) continue;

      /* skip node if it is not a center node */
      if (NTYPE(theNode)!=CENTER_NODE) continue;

      /* get father element */
      theVertex=MYVERTEX(theNode);
      fatherElement = VFATHER(theVertex);

      /* search for curved boundaries */
      if (OBJT(fatherElement)!=BEOBJ) continue;
      if (!MovedNode(fatherElement)) continue;

      if (option==3)
      {
        /* determine default positions for center node and mid node opposite to boundary */
        if (DefaultPosCurvedBoundary(fatherElement,theNode,TRUE,
                                     MidNodeLambdaNew,MidNodeLambdaOld)) continue;
      }
      else
      {
        /* determine new positions for center node and mid node opposite to boundary */
        if (NewPosCurvedBoundary(fatherElement,theNode,MidNodeLambdaNew)) continue;
      }

    }

    /* move boundary side nodes of quadrilaterals in order to get orthogonal boundary elements  on
       specified boundaries */
    if (option==1 || option==2)       /* search bnd mid nodes  */
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      {
        /* skip node if it is a copy from a lower level */
        if (CORNERTYPE(theNode)) continue;

        /* skip node if it is not a mid node (mid point of an edge) */
        if (NTYPE(theNode)!=MID_NODE) continue;

        /* skip node if its not a boundary node */
        theVertex = MYVERTEX(theNode);
        fatherElement = VFATHER(theVertex);
        if (OBJT(fatherElement)!=BEOBJ) continue;
        if (OBJT(theVertex)!=BVOBJ) continue;

        /* find specified boundary */
        if (!ElemOnSpecBnd(fatherElement,bnd,bnd_num,&side)) continue;

        /* apply 'ortho' option for boundary mid-node if possible */
        if (LambdaOrthoBnd2D(fatherElement,side,theNode,&lambda)) continue;

        /* apply limits for lambda  */
        lambda = MAX(lambda,0.5-LimitLocDis);
        lambda = MIN(lambda,0.5+LimitLocDis);
        MidNodeLambdaNew[ID(MYVERTEX(theNode))] = lambda;
      }

    /* move all son nodes of an level-0-boundary-element in order to get 'orthogonal' elements */
    if (option==2)
    {
      /*    move center nodes of quadrilaterals  */
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      {
        /* skip node if it is a copy from a lower level */
        if (CORNERTYPE(theNode)) continue;

        /* skip node if it is not a center node */
        if (NTYPE(theNode)!=CENTER_NODE) continue;

        theVertex=MYVERTEX(theNode);
        fatherElement = VFATHER(theVertex);

        /* use only boundary father elements  */
        /* if (OBJT(fatherElement)!=BEOBJ) continue; */

        /* skip node if it is not on a regular refined element */
        if (REFINECLASS(fatherElement)!=RED_CLASS) continue;

        /* search father element on level 0 */
        Level0Father=fatherElement;
        for (i=GLEVEL(theGrid)-1; i>0; i--) Level0Father=EFATHER(Level0Father);
        /* the father element on level 0 must be a boundary element */
        if (OBJT(Level0Father)!=BEOBJ) continue;

        if (OBJT(fatherElement)==BEOBJ)
        {
          /* check if father element is on wall boundary */
          if (!ElemOnSpecBnd(fatherElement,bnd,bnd_num,&side)) continue;

          if (EqualDistBndElem(fatherElement,side,theNode,MidNodeLambdaNew)!=0) continue;

        }
        else
        {
          /* check if father element is on wall boundary */
          if (!ElemOnSpecBnd(Level0Father,bnd,bnd_num,&side)) continue;

          /* find corresponding element at the boundary */
          for (boundaryFather=NBELEM(fatherElement,side); boundaryFather!=0;
               boundaryFather=NBELEM(boundaryFather,side))
          {
            /* otherwise we will never find the correct boundary element */
            if (REFINECLASS(boundaryFather)!=RED_CLASS) break;
            if(OBJT(boundaryFather)==BEOBJ)
            {
              if (EqualDistInnerElem(fatherElement,boundaryFather,side,theNode,
                                     MidNodeLambdaNew)!=0) continue;
            }
          }
        }

      }
    }

    /* update local and global Coordinates and move nodes to their new position */
    if (MoveNodesOnGrid(theGrid,VertexCoord,VertexLCoord,MidNodeLambdaOld,MidNodeLambdaNew,
                        LimitLocDis))
    {
      PrintErrorMessage('E',"SmoothGrid","Error in MoveNodesOnGrid");
      goto exit;
    }

    /* reset node flag */
    for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
      SETTHEFLAG(theNode,0);
  }
  ReleaseTmpMem(MGHEAP(theMG),MarkKey);
  return(0);

exit:
  for (lev=fl; lev<=tl; lev++)
  {
    theGrid=GRID_ON_LEVEL(theMG,lev);
    /* reset node flag */
    for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
      SETTHEFLAG(theNode,0);
  }
  ReleaseTmpMem(MGHEAP(theMG),MarkKey);
  return(1);

}
/****************************************************************************/
/*
   SmoothGridReset - reset all nodes on specified grids to the default refinement position

   SYNOPSIS
   static INT SmoothGridReset (MULTIGRID *theMG, INT fl, INT tl)

   PARAMETERS
   .  theMG        - Multigrid
   .  fl, tl       - reset elements of grids on level fl<=level<=tl

   DESCRIPTION
   Reset all center nodes and mid nodes of the grid to the default position.

   RETURN VALUE
   INT
   .n   0: ok
   .n   1: error
 */
/****************************************************************************/

INT SmoothGridReset (MULTIGRID *theMG, INT fl, INT tl)
{
  GRID *theGrid;
  ELEMENT *fatherElement;
  VERTEX *theVertex;
  NODE *theNode;
  DOUBLE *CornerPtrs[MAX_CORNERS_OF_ELEM];
  INT coe,lev;
  DOUBLE_VECTOR *VertexCoord, *VertexLCoord;
  DOUBLE *MidNodeLambdaOld, *MidNodeLambdaNew,lambda_old;
  INT MarkKey;

  /* allocate temporary memory */
  MarkTmpMem(MGHEAP(theMG),&MarkKey);
  VertexCoord = (DOUBLE_VECTOR *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE_VECTOR),MarkKey);
  VertexLCoord = (DOUBLE_VECTOR *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE_VECTOR),MarkKey);
  MidNodeLambdaOld = (DOUBLE *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE),MarkKey);
  MidNodeLambdaNew = (DOUBLE *) GetTmpMem(MGHEAP(theMG),VIDCNT(theMG)*sizeof(DOUBLE),MarkKey);

  for (lev=fl; lev<=tl; lev++)
  {
    theGrid=GRID_ON_LEVEL(theMG,lev);
    /* copy current global and local vertex coordinates */
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
    {
      V_DIM_COPY(CVECT(theVertex),VertexCoord[ID(theVertex)]);
      V_DIM_COPY(LCVECT(theVertex),VertexLCoord[ID(theVertex)]);
    }

    /* check node flags */
    for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
      if (THEFLAG(theNode))
      {
        PrintErrorMessage('E',"SmoothGridReset","node flag already set");
        return(1);
      }

    /*    re-move center nodes of quadrilaterals  */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* skip node if it is a copy from a lower level */
      if (CORNERTYPE(theNode)) continue;

      /* skip node if it is not a center node */
      if (NTYPE(theNode)!=CENTER_NODE) continue;

      /* get father element */
      theVertex=MYVERTEX(theNode);
      fatherElement = VFATHER(theVertex);

      /* default local coordinates of center node */
      V_DIM_SET(0.5,LCVECT(theVertex));
      /* calculate new global coordinates */
      CORNER_COORDINATES(fatherElement,coe,CornerPtrs);
      LOCAL_TO_GLOBAL(coe,CornerPtrs,LCVECT(theVertex),CVECT(theVertex));
    }

    /* re-move mid nodes */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* skip node if it is a copy from a lower level */
      if (CORNERTYPE(theNode)) continue;

      /* skip node if it is not a mid node (mid point of an edge) */
      if (NTYPE(theNode)!=MID_NODE) continue;

      if (GetMidNodeParam(theNode,&lambda_old)!=0) return(1);
      MidNodeLambdaOld[ID(MYVERTEX(theNode))] = lambda_old;
      /* default value for lambda */
      MidNodeLambdaNew[ID(MYVERTEX(theNode))] = 0.5;
    }

    /* special treatment for center and mid nodes of boundary elements on a curved boundary */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* skip node if it is a copy from a lower level */
      if (CORNERTYPE(theNode)) continue;

      /* skip node if it is not a center node */
      if (NTYPE(theNode)!=CENTER_NODE) continue;

      /* get father element */
      theVertex=MYVERTEX(theNode);
      fatherElement = VFATHER(theVertex);

      /* search for curved boundaries */
      if (OBJT(fatherElement)!=BEOBJ) continue;
      if (!MovedNode(fatherElement)) continue;

      /* determine default positions for center node and mid node opposite to boundary */
      if (DefaultPosCurvedBoundary(fatherElement,theNode,FALSE,
                                   MidNodeLambdaNew,MidNodeLambdaOld)) continue;

    }

    /* update local and global Coordinates and move nodes to their new position */
    if (MoveNodesOnGrid(theGrid,VertexCoord,VertexLCoord,MidNodeLambdaOld,MidNodeLambdaNew,
                        0.3))
    {
      PrintErrorMessage('E',"SmoothGridReset","Error in MoveNodesOnGrid");
      return(1);
    }

    /* reset flag */
    for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
      SETTHEFLAG(theNode,0);
  }
  ReleaseTmpMem(MGHEAP(theMG),MarkKey);
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
            PrintErrorMessage('W',"SmoothMultiGrid",
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
