// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugm2d.c														*/
/*																			*/
/* Purpose:   Functions only available in 2D								*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de					                */
/*																			*/
/* History:   13.12.94 begin, ug version 3.0								*/
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

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>

#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"

#include "devices.h"

#include "switch.h"
#include "gm.h"
#include "misc.h"
#include "ugm.h"
#include "ugm2d.h"
#include "shapes2d.h"


#ifdef __THREEDIM__
#error this source file is for 2D ONLY
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

#define SMALL1 0.001


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



/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   FindFather - Find the new father element

   SYNOPSIS:
   static ELEMENT *FindFather(VERTEX *vptr);

   PARAMETERS:
   .  vptr - Pointer to 'VERTEX' whose father element is to be found.

   DESCRIPTION:
   This function finds the new father element of the given vertex.
   It assumes that the  new father is one of the neighbors of the
   old father element. The current father of 'vptr' is not changed.

   RETURN VALUE:
   ELEMENT *
   .n     pointer to an element
   .n     NULL if none or no correct father is found or vertex is level 0
   D*/
/****************************************************************************/

static ELEMENT *FindFather(VERTEX *vptr)
{
  short i;
  COORD x[2],lambda,lambda1,lambda2;
  ELEMENT *eptr,*eptr1;
  ELEMENTSIDE *eside;
  VSEGMENT *vseg;

  if ((eptr=VFATHER(vptr))==NULL)
    return(NULL);

  if (OBJT(vptr)==BVOBJ)
  {
    /* for boundary vertices only a consistency check is made */
    eside=SIDE(eptr,ONEDGE(vptr));
    vseg = VSEG(vptr);

    if (VS_PATCH(vseg)!=ES_PATCH(eside))
      return(NULL);

    /* for higher dimensions (3d) the following must be generalized */
    lambda=LAMBDA(vseg,0);
    lambda1=PARAM(eside,0,0);
    lambda2=PARAM(eside,1,0);
    if ( ((lambda1<=lambda)&&(lambda2>=lambda)) || ((lambda2<=lambda)&&(lambda1>=lambda)) )
      return(eptr);
    else
      return(NULL);
  }

  x[0]=XC(vptr); x[1]=YC(vptr);

  eptr=VFATHER(vptr);
  if (PointInElement(x,eptr))
    return(eptr);

  for (i=0; i<SIDES_OF_ELEM(eptr); i++)
    if (PointInElement(x,(eptr1=NBELEM(eptr,i))))
      return(eptr1);

  return(NULL);
}



/****************************************************************************/
/*D
   Local2Global - Updates the global coordinates of a vertex

   SYNOPSIS:
   static INT Local2Global (MULTIGRID *theMG, VERTEX *vptr);

   PARAMETERS:
   .  theMG - pointer to 'MULTIGRID' structure.
   .  vptr - pointer to a 'VERTEX'

   DESCRIPTION:
   This function updates the global coordinates of a vertex by
   evaluating the local coordinates in the father element of the vertex.
   This function overwrites the current coordinates, it works for interior
   and boundary vertices.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   >0 if an error occured or vertex is on level 0
   D*/
/****************************************************************************/

static INT Local2Global (MULTIGRID *theMG, VERTEX *vptr)
{
  short i,n,side;
  COORD x[DIM],xi,eta;
  COORD lambda1,lambda2,lambdaa,lambdae;
  ELEMENT *eptr;
  VERTEX *vptr1,*vptr2,*vptra,*vptre;
  VSEGMENT *vseg;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;

  if ((eptr=VFATHER(vptr))==NULL)
    return(8040);

  n=TAG(eptr);
  if (OBJT(vptr)==BVOBJ)
  {
    vseg = VSEG(vptr);
    if ((thePatch=VS_PATCH(vseg))==NULL) return(8041);
    if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return(8041);

    side=ONEDGE(vptr);
    vptr1=MYVERTEX(CORNER(eptr,side));
    vptr2=MYVERTEX(CORNER(eptr,(side+1)%n));

    vptra=theMG->corners[PATCH_CID(thePatchDesc,0)]; lambdaa=PATCH_LCVECT(thePatchDesc,0)[0];
    vptre=theMG->corners[PATCH_CID(thePatchDesc,1)]; lambdae=PATCH_LCVECT(thePatchDesc,1)[0];

    if (vptr1==vptra)
      lambda1=lambdaa;
    else
    if (vptr1==vptre)
      lambda1=lambdae;
    else
      lambda1=LAMBDA(VSEG(vptr1),0);

    if (vptr2==vptra)
      lambda2=lambdaa;
    else
    if (vptr2==vptre)
      lambda2=lambdae;
    else
      lambda2=LAMBDA(VSEG(vptr2),0);


    LAMBDA(vseg,0)=(1-ZETA(vseg))*lambda1+ZETA(vseg)*lambda2;
    if (Patch_local2global(thePatch,PVECT(vseg),CVECT(vptr))) return(8041);
  }
  else
  {
    xi=XI(vptr); eta=ETA(vptr);
    x[0]=x[1]=0;
    for (i=0; i<n; i++)
    {
      vptr1=MYVERTEX(CORNER(eptr,i));
      x[0]+=XC(vptr1)*N(n,i,xi,eta);
      x[1]+=YC(vptr1)*N(n,i,xi,eta);
    }

    XC(vptr)=x[0]; YC(vptr)=x[1];
  }

  return(0);
}

/****************************************************************************/
/*D
   Global2Local - Updates the local coordinates of a vertex

   SYNOPSIS:
   static INT Global2Local (MULTIGRID *theMG, VERTEX *vptr);

   PARAMETERS:
   .  theMG - pointer to 'MULTIGRID'
   .  vptr - pointer to a 'VERTEX'

   DESCRIPTION:
   This function updates the local coordinates of a vertex in the
   local coordinate system of the father element. It is assumed that
   the vertex is inside the father element. If not an error is returned.

   RETURN VALUE:
   INT
   .n    0 if ok vertex is level 0
   .n    >0 if error occured (invalid trafo, not inside)
   D*/
/***************************************************************************/

static INT Global2Local (MULTIGRID *theMG, VERTEX *vptr)
{
  short n,side;
  DOUBLE t1x,t2x,t3x,t1y,t2y,t3y,a,b,c,D,xi1,xi2,eta1,eta2;
  DOUBLE x1,x2,y1,y2,x3,y3,x4,y4,x,y;
  DOUBLE a1,a2,a3,a4,b1,b2,b3,b4;
  COORD lambda1,lambda2,lambdaa,lambdae;
  ELEMENT *eptr;
  VERTEX *vptr1,*vptr2,*vptra,*vptre;
  VSEGMENT *vseg;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;

  if ((eptr=VFATHER(vptr))==NULL)
    return(0);

  if (OBJT(vptr)==BVOBJ)
  {
    n=TAG(eptr);
    vseg = VSEG(vptr);
    thePatch=VS_PATCH(vseg);
    if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return (1);
    side=ONEDGE(vptr);
    vptr1=MYVERTEX(CORNER(eptr,side));
    vptr2=MYVERTEX(CORNER(eptr,(side+1)%n));

    vptra=theMG->corners[PATCH_CID(thePatchDesc,0)]; lambdaa=PATCH_LCVECT(thePatchDesc,0)[0];
    vptre=theMG->corners[PATCH_CID(thePatchDesc,1)]; lambdae=PATCH_LCVECT(thePatchDesc,1)[0];

    if (vptr1==vptra)
      lambda1=lambdaa;
    else
    if (vptr1==vptre)
      lambda1=lambdae;
    else
      lambda1=LAMBDA(VSEG(vptr1),0);

    if (vptr2==vptra)
      lambda2=lambdaa;
    else
    if (vptr2==vptre)
      lambda2=lambdae;
    else
      lambda2=LAMBDA(VSEG(vptr2),0);

    /* falls sich Projektion auf Edge als sinnvoller erweist
       x = XC(vptr); y = YC(vptr);
       x1 = XC(vptr1); y1 = YC(vptr1);
       x2 = XC(vptr2); y2 = YC(vptr2);
       a=sqrt(pow(x2-x1,2)+pow(y2-y1,2));
       b=sqrt(pow(x-x1,2)+pow(y-y1,2));
       assert(errno==0);
       if (a*b==0)
            return(8052);
       c=((x-x1)*(x2-x1)+(y-y1)*(y2-y1))/(a*b);
       if ((c<.1)||(c>.9))
            return(8053);
     */

    /* set local boundary coordinate */
    ZETA(vseg)=c=(LAMBDA(vseg,0)-lambda1)/(lambda2-lambda1);

    /* set local coordinates */
    switch(n)
    {
    case TRIANGLE :
      XI(vptr) = ETA(vptr) = 0.0;
      switch (side)
      {
      case 0 : XI(vptr)=c; break;
      case 1 : XI(vptr)=1-c; ETA(vptr)=c; break;
      case 2 : ETA(vptr)=1-c; break;
      }
      break;

    case QUADRILATERAL :
      switch (side)
      {
      case 0 : XI(vptr)=2*c-1; ETA(vptr) = -1;  break;
      case 1 : XI(vptr) = 1;   ETA(vptr)=2*c-1; break;
      case 2 : XI(vptr)=1-2*c; ETA(vptr) = 1;   break;
      case 3 : XI(vptr) = -1;  ETA(vptr)=1-2*c; break;
      }
      break;
    }
  }
  else
  {
    /* compute local coordinates */
    vptr1=MYVERTEX(CORNER(eptr,0)); x1 = XC(vptr1); y1 = YC(vptr1);
    vptr1=MYVERTEX(CORNER(eptr,1)); x2 = XC(vptr1); y2 = YC(vptr1);
    vptr1=MYVERTEX(CORNER(eptr,2)); x3 = XC(vptr1); y3 = YC(vptr1);
    x = XC(vptr);
    y = YC(vptr);

    switch (TAG(eptr))
    {
    case TRIANGLE :
      t1x = x2-x1; t2x = x3-x1; t3x = x1;
      t1y = y2-y1; t2y = y3-y1; t3y = y1;
      D = t1x*t2y-t2x*t1y;
      if (D<0.0)
        return(8051);

      XI(vptr) = (t2y*(x-t3x)-t2x*(y-t3y))/D;
      ETA(vptr) = (-t1y*(x-t3x)+t1x*(y-t3y))/D;
      break;

    case QUADRILATERAL :
      vptr1=MYVERTEX(CORNER(eptr,3)); x4 = XC(vptr1); y4 = YC(vptr1);

      a1=-x+.25*(x1+x2+x3+x4); b1=-y+.25*(y1+y2+y3+y4);
      a2=.25*(-x1+x2+x3-x4); b2=.25*(-y1+y2+y3-y4);
      a3=.25*(-x1-x2+x3+x4); b3=.25*(-y1-y2+y3+y4);
      a4=.25*(x1-x2+x3-x4); b4=.25*(y1-y2+y3-y4);

      c=a2*b1-a1*b2;
      b=a4*b1-a3*b2+a2*b3-a1*b4;
      a=a4*b3-a3*b4;

      if (ABS(a)<SMALL_D)
      {
        if (ABS(b)<SMALL_D)
          return(8051);

        eta1 = eta2 = -c/b;
      }
      else
      {
        D=b*b-4*a*c;
        if (D<0) return(8051);

        eta1=(-b+sqrt(D))/(2*a); eta2=(-b-sqrt(D))/(2*a);
      }

      c=a3*b1-a1*b3;
      b=a4*b1+a3*b2-a2*b3-a1*b4;
      a=a4*b2-a2*b4;

      if (ABS(a)<SMALL_D)
      {
        if (ABS(b)<SMALL_D)
          return(8051);

        xi1 = xi2 = -c/b;
      }
      else
      {
        D=b*b-4*a*c;
        if (D<0) return(8051);

        xi1=(-b+sqrt(D))/(2*a); xi2=(-b-sqrt(D))/(2*a);
      }

      if ((xi1>=-1-SMALL1)&&(xi1<=1+SMALL1)&&(eta1>=-1-SMALL1)&&(eta1<=1+SMALL1))
      {
        XI(vptr) = xi1;
        ETA(vptr) = eta1;
        break;
      }

      if ((xi2>=-1-SMALL1)&&(xi2<=1+SMALL1)&&(eta1>=-1-SMALL1)&&(eta1<=1+SMALL1))
      {
        XI(vptr) = xi2;
        ETA(vptr) = eta1;
        break;
      }

      if ((xi1>=-1-SMALL1)&&(xi1<=1+SMALL1)&&(eta2>=-1-SMALL1)&&(eta2<=1+SMALL1))
      {
        XI(vptr) = xi1;
        ETA(vptr) = eta2;
        break;
      }

      if ((xi2>=-1-SMALL1)&&(xi2<=1+SMALL1)&&(eta2>=-1-SMALL1)&&(eta2<=1+SMALL1))
      {
        XI(vptr) = xi2;
        ETA(vptr) = eta2;
        break;
      }

      return(8051);
      break;
    }
  }

  return(0);
}

/****************************************************************************/
/*D
   MoveInnerNode - Let user enter a new position for an inner node

   SYNOPSIS:
   INT MoveInnerNode (MULTIGRID *theMG, NODE *theNode, COORD *newPos);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  newPos - new position (x,y)

   DESCRIPTION:
   This function moves a given node to a new position. The complete
   multigrid structure is moved hierachically, that all global coordinates
   of nodes are updated in order to reflect the changes on coarser grids.

   `Function only implemented in 2D version !`

   RETURN VALUE:
   INT
   .n   0 when ok
   .n   >0 when error occured.
   D*/
/****************************************************************************/

INT MoveInnerNode (MULTIGRID *theMG, NODE *theNode, COORD *newPos)
{
  GRID *theGrid2;
  int k,k2;
  NODE *theNode2;
  VERTEX *theVertex,*theVertex2;
  ELEMENT *theElement,*oldElement;
  double x,y,oldx,oldy;

  k = LEVEL(theNode);

  /* set k (and theNode) to the level where the node appears the first time */
  while ((theNode2=NFATHER(theNode))!=NULL)
  {
    theNode=theNode2;
    k=k-1;
  }

  theVertex = MYVERTEX(theNode);
  oldx = XC(theVertex); oldy = YC(theVertex);

  if (OBJT(theVertex)!=IVOBJ)
  {
    PrintErrorMessage('E',"MoveInnerNode","no inner node passed");
    return(GM_ERROR);
  }

  x = newPos[0]; y = newPos[1];

  /* set values */
  XC(theVertex) = x; YC(theVertex) = y;
  if (VFATHER(theVertex)!=NULL)
  {
    oldElement=VFATHER(theVertex);
    if ((theElement=FindFather(theVertex))==NULL)
    {
      PrintErrorMessage('E',"MoveInnerNode","No father element! Probably you tried to move the vertex too far!");
      XC(theVertex) = oldx; YC(theVertex) = oldy;
      return(GM_ERROR);
    }
    else
      VFATHER(theVertex)=theElement;

    if (Global2Local(theMG,theVertex)!=0)
    {
      PrintErrorMessage('E',"MoveInnerNode","Error in Global2Local");
      VFATHER(theVertex)=oldElement;
      XC(theVertex) = oldx; YC(theVertex) = oldy;
      return(GM_ERROR);
    }
  }

  /*	now we correct the global coordinates for all levels above, since it is not
          easy to find exactly the vertices whose global coordinates have changed */

  for(k2=k+1; k2<=theMG->topLevel; k2++)
  {
    theGrid2 = theMG->grids[k2];
    for (theVertex2=FIRSTVERTEX(theGrid2); theVertex2!=NULL; theVertex2=SUCCV(theVertex2))
      if (Local2Global(theMG,theVertex2)!=0)
      {
        PrintErrorMessage('E',"MoveInnerNode","Fatal error in correcting global coordinates for higher levels. Grid may be inconsistent from now on.");
        return(GM_ERROR);
      }
  }

  /* OK, done */
  return(GM_OK);
}

/****************************************************************************/
/*D
   MoveBoundaryNode - Let user enter a new position

   SYNOPSIS:
   INT MoveBoundaryNode (MULTIGRID *theMG, NODE *theNode, INT segid,
   COORD *newPos);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theNode - node to move
   .  segid - new boundary segment id
   .  newPos - new parameter lambda

   DESCRIPTION:
   This function moves a boundary node to a new position. The position of
   all nodes on finer grids is updated to reflect these changes.

   `Function only implemented in 2D version!`

   RETURN VALUE:
   INT
   .n    0 when ok
   .n    >0 when error occured.
   D*/
/****************************************************************************/

INT MoveBoundaryNode (MULTIGRID *theMG, NODE *theNode, INT patchid, COORD *newPos)
{
  GRID *theGrid,*theGrid2;
  int i,n,k,k2;
  NODE *theNode2;
  VERTEX *theVertex,*theVertex2;
  ELEMENT *theElement,*oldElement;
  ELEMENTSIDE *theSide;
  double oldx,oldy,l,oldl;
  BVP             *theBVP;
  BVP_DESC theBVPDesc;
  PATCH *thePatch, *oldPatch;
  PATCH_DESC thePatchDesc;

  /* get BVP description */
  theBVP = MG_BVP(theMG);
  if (BVP_GetBVPDesc(theBVP,&theBVPDesc))
  {
    PrintErrorMessage('E',"MoveBoundaryNode","cannot evaluate BVP");
    return(1);
  }

  k = LEVEL(theNode);
  theGrid = GRID_ON_LEVEL(theMG,k);

  /* set k (and theNode) to the level where the node appears the first time */
  while ((theNode2=NFATHER(theNode))!=NULL)
  {
    theNode=theNode2;
    k=k-1;
  }

  theVertex = MYVERTEX(theNode);
  oldx = XC(theVertex); oldy = YC(theVertex);

  if (OBJT(theVertex)!=BVOBJ)
  {
    PrintErrorMessage('E',"MoveBoundaryNode","no boundary node passed");
    return(GM_ERROR);
  }

  if (MOVE(theVertex)==0)
  {
    PrintErrorMessage('W',"MoveBoundaryNode","corners cannot be moved");
    return(GM_ERROR);
  }

  l = oldl = LAMBDA(VSEG(theVertex),0);
  oldPatch = VS_PATCH(VSEG(theVertex));
  if (START(theNode)==NULL)
  {
    if(patchid >= BVPD_NPATCHES(theBVPDesc))
    {
      PrintErrorMessage('E',"MoveBoundaryNode","patchid out of range");
      return(GM_ERROR);
    }
    thePatch = Patch_GetPatchByID(theBVP,patchid);
  }
  else thePatch = VS_PATCH(VSEG(theVertex));

  if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return (1);
  patchid = PATCH_ID(thePatchDesc);

  l = newPos[0];

  if ((l<PATCH_LCVECT(thePatchDesc,0)[0]) || (l>PATCH_LCVECT(thePatchDesc,1)[0]))
  {
    PrintErrorMessage('E',"MoveBoundaryNode","parameter out of range");
    return(GM_ERROR);
  }

  LAMBDA(VSEG(theVertex),0) = l;
  if (Patch_local2global(thePatch,PVECT(VSEG(theVertex)),CVECT(theVertex))) return (1);
  VS_PATCH(VSEG(theVertex)) = thePatch;

  if (VFATHER(theVertex)!=NULL)
  {
    oldElement=VFATHER(theVertex);
    if ((theElement=FindFather(theVertex))==NULL)
    {
      PrintErrorMessage('E',"MoveBoundaryNode","No father element! Probably you have tried to move the vertex too far");
      XC(theVertex) = oldx; YC(theVertex) = oldy; LAMBDA(VSEG(theVertex),0)=oldl; VS_PATCH(VSEG(theVertex)) = oldPatch;
      return(GM_ERROR);
    }
    else
      VFATHER(theVertex)=theElement;

    if (Global2Local(theMG,theVertex)!=0)
    {
      PrintErrorMessage('E',"MoveBoundaryNode","Error in Global2Local!");
      VFATHER(theVertex)=oldElement;
      XC(theVertex) = oldx; YC(theVertex) = oldy; LAMBDA(VSEG(theVertex),0)=oldl;
      return(GM_ERROR);
    }
  }

  /*	now we correct the global coordinates for all new vertices in the levels above, since it is not
          easy to find exactly the vertices whose global coordinates have changed */
  for(k2=k+1; k2<=theMG->topLevel; k2++)
  {
    theGrid2 = theMG->grids[k2];
    for (theVertex2=FIRSTVERTEX(theGrid2); theVertex2!=NULL; theVertex2=SUCCV(theVertex2))
      if (Local2Global(theMG,theVertex2)!=0)
      {
        PrintErrorMessage('E',"MoveBoundaryNode","Fatal error in correcting global coordinates for higher levels. Grid may be inconsistent from now on");
        return(GM_ERROR);
      }
  }

  /* update element sides on this level */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if (OBJT(theElement)==BEOBJ)
    {
      n = SIDES_OF_ELEM(theElement);
      for (i=0; i<n; i++)
      {
        theSide = SIDE(theElement,i);
        if (theSide!=NULL)
        {
          if (MYVERTEX(CORNER(theElement,i))==theVertex)
            PARAM(theSide,0,0) = l;
          if (MYVERTEX(CORNER(theElement,(i+1)%n))==theVertex)
            PARAM(theSide,1,0) = l;
        }
      }
    }

  /* update (for simplicity) all sides on that segment for all levels above */
  for(k2=k+1; k2<=theMG->topLevel; k2++)
  {
    theGrid2 = theMG->grids[k2];
    for (theElement=theGrid2->elements; theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (OBJT(theElement)==BEOBJ)
      {
        n = SIDES_OF_ELEM(theElement);
        for (i=0; i<n; i++)
        {
          theSide = SIDE(theElement,i);
          if (theSide!=NULL)
          {
            theNode2=CORNER(theElement,i);
            if (NFATHER(theNode2)==NULL)
            {
              theVertex2=MYVERTEX(theNode2);
              if (VS_PATCH(VSEG(theVertex2))==thePatch)
                PARAM(theSide,0,0) = LAMBDA(VSEG(theVertex2),0);
            }
            theNode2=CORNER(theElement,(i+1)%n);
            if (NFATHER(theNode2)==NULL)
            {
              theVertex2=MYVERTEX(theNode2);
              if (VS_PATCH(VSEG(theVertex2))==VS_PATCH(VSEG(theVertex)))
                PARAM(theSide,1,0) = LAMBDA(VSEG(theVertex2),0);
            }
          }
        }
      }
    }
  }

  return(GM_OK);
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

   `This function is available in 2D version only!`

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 if error occured.
   D*/
/****************************************************************************/

INT SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag)
{
  int l,i,n;
  double ratio;
  DOUBLE N;
  /* COORD beta,gamma; */
  COORD lambda,lambda0,lambda1,lambda2,lambdaa,lambdae;
  GRID *theGrid;
  NODE *node,*node2;
  ELEMENT *eptr;
  ELEMENTSIDE *eside;
  VERTEX *vptr0,*vptr1,*vptr2,*vptra,*vptre,*vptr;
  LINK *lptr;
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;

  COORD x[2];

  n = niter;
  if (n<=0) n = 1;
  if (n>50) n = 50;

  ratio=.5;

  for (i=0; i<n; i++)
  {
    for (l=0; l<=theMG->topLevel; l++)
    {
      theGrid=theMG->grids[l];

      /* update global coordinates of new nodes */
      if (l!=0)
        for (node=theGrid->firstNode; node!=NULL; node=SUCCN(node))
          if (NFATHER(node)==NULL)
          {
            vptr=MYVERTEX(node);
            if ((OBJT(vptr)!=BVOBJ)||(bdryFlag!=0))
              if (Local2Global(theMG,vptr)!=0)
                return(GM_ERROR);
          }

      for (node=theGrid->firstNode; node!=NULL; node=SUCCN(node))
      {
        /* skip node if it is a copy from a lower level */
        if (NFATHER(node)!=NULL) continue;

        vptr0=MYVERTEX(node);
        if (OBJT(vptr0)==BVOBJ)
        {
          if (bdryFlag==0) continue;                                            /* boundary is not allowed to move */
          if (MOVE(vptr0)==0) continue;                                 /* corner: do not move it */

          /* test if free boundary: since that type is determined by the grid it may only
             be smoothed in certain cases and in a special way */
          thePatch=VS_PATCH(VSEG(vptr0));
          if(Patch_GetPatchDesc(thePatch,&thePatchDesc)) return(1);
          if (PATCH_TYPE(thePatchDesc)==FREE)
            continue;                                           /*free boundary is not allowed to be smoothed */

          /* first find endpoints of vptr0's boundary segment */
          lambda0=LAMBDA(VSEG(vptr0),0);
          vptra=theMG->corners[PATCH_CID(thePatchDesc,0)]; lambdaa=PATCH_LCVECT(thePatchDesc,0)[0];
          vptre=theMG->corners[PATCH_CID(thePatchDesc,1)]; lambdae=PATCH_LCVECT(thePatchDesc,1)[0];

          /* search for the two nearest neighbors on the same segment */
          vptr1=vptr2=NULL;
          for (lptr=START(node); lptr!=NULL; lptr=NEXT(lptr))
            if (EXTRA(MYEDGE(lptr))==0)
            {
              node2=NBNODE(lptr);
              vptr=MYVERTEX(node2);
              if (OBJT(vptr)!=BVOBJ)
                continue;

              /* if on same segment get lambda value */
              if (vptr==vptra)
                lambda=lambdaa;
              else
              if (vptr==vptre)
                lambda=lambdae;
              else
              if (thePatch==VS_PATCH(VSEG(vptr)))
                lambda=LAMBDA(VSEG(vptr),0);
              else
                continue;

              if (lambda>lambda0)
              {
                if (vptr2==NULL)
                {
                  vptr2=vptr;
                  lambda2=lambda;
                }
                else
                if (lambda<lambda2)
                {
                  vptr2=vptr;
                  lambda2=lambda;
                }
              }
              else
              {
                if (vptr1==NULL)
                {
                  vptr1=vptr;
                  lambda1=lambda;
                }
                else
                if (lambda>lambda1)
                {
                  vptr1=vptr;
                  lambda1=lambda;
                }
              }
            }

          if ((vptr1==NULL)||(vptr2==NULL))
            return(GM_ERROR);

          if (PATCH_TYPE(thePatchDesc)==FREE)
          {
            /*	This is only sensible if the free boundary is a line and the endpoints are moved.
                    A more general smoothing would be to use e.g. a quadratic interpolant. */
            XC(vptr0)=.5*(XC(vptr1)+XC(vptr2));
            YC(vptr0)=.5*(YC(vptr1)+YC(vptr2));
            /*	since local and global boundary coordinates are given by the grid itself
                    they are assumed to be correct. */
          }
          else
          {
            LAMBDA(VSEG(vptr0),0)=.5*(lambda1+lambda2);
            /* set global coordinates */
            if (Patch_local2global(thePatch,PVECT(VSEG(vptr0)),CVECT(vptr0))) return (GM_ERROR);

            /* set local boundary coordinates */
            if (Global2Local(theMG,vptr0)!=0)
              return(GM_ERROR);
          }
        }
        else
        {
          x[0]=x[1]=0; N=0;

          for (lptr=START(node); lptr!=NULL; lptr=NEXT(lptr))
          {
            node2=NBNODE(lptr);
            vptr=MYVERTEX(node2);

            if (EXTRA(MYEDGE(lptr))==0)
            {
              x[0]+=XC(vptr);
              x[1]+=YC(vptr);
              N+=1;
            }
            else
            {
              x[0]+=ratio*XC(vptr);
              x[1]+=ratio*YC(vptr);
              N+= ratio;
            }
          }

          XC(vptr0)=x[0]/N; YC(vptr0)=x[1]/N;

          /* if there is a father element, change local variables */
          if (l!=0)
            if ((eptr=FindFather(vptr0))!=NULL)
            {
              VFATHER(vptr0)=eptr;
              if (Global2Local(theMG,vptr0)!=0)
                return(GM_ERROR);
            }
            else
              return(GM_ERROR);
        }
      }

      /* at last, the boundary element sides of this level must be updated! */
      for (eptr=theGrid->elements; eptr!=NULL; eptr=SUCCE(eptr))
        if (OBJT(eptr)==BEOBJ)
        {
          n=SIDES_OF_ELEM(eptr);
          for (i=0; i<n; i++)
            if ((eside=SIDE(eptr,i))!=NULL)
            {
              /* for higher dimensions (3d) the following must be generalized */
              vptr=MYVERTEX(CORNER(eptr,i));
              if (MOVE(vptr)!=0)
                PARAM(eside,0,0)=LAMBDA(VSEG(vptr),0);
              vptr=MYVERTEX(CORNER(eptr,(i+1)%n));
              if (MOVE(vptr)!=0)
                PARAM(eside,1,0)=LAMBDA(VSEG(vptr),0);
            }
        }
    }
  }

  return(GM_OK);
}

/****************************************************************************/
/*D
   CreateAuxEdge - Inserts an 'EDGE' with AUXEDGE flag set

   SYNOPSIS:
   EDGE *CreateAuxEdge (GRID *theGrid, NODE *from, NODE *to);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  from - starting 'NODE'
   .  to - destination 'NODE'

   DESCRIPTION:
   This function inserts a new edge between given nodes.
   This function is useful only in version 2.3 compatibility mode which is
   not documented.

   RETURN VALUE:
   INT
   .n     NULL when alloc failed
   .n     else it returns a pointer to the new edge
   D*/
/****************************************************************************/

EDGE *CreateAuxEdge (GRID *theGrid, NODE *from, NODE *to)
{
  EDGE *theEdge;

  theEdge = CreateEdge(theGrid,from,to);
  if (theEdge==NULL) return(NULL);
  SETEXTRA(theEdge,1);
  SETAUXEDGE(theEdge,1);
  return(theEdge);
}

/****************************************************************************/
/*D
   DisposeAuxEdges - Remove edges with AUXEDGE flag set from theGrid

   SYNOPSIS:
   INT DisposeAuxEdges (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function removes all 'EDGE' objects with AUXEDGE flag set from the
   given grid level. This function is provided for compatibility with
   version 2.3.

   RETURN VALUE:
   INT
   .n    -1 in case of an error
   .n    else it returns the number of edges that have been disposed
   D*/
/****************************************************************************/

INT DisposeAuxEdges (GRID *theGrid)
{
  NODE *theNode;
  LINK *theLink;
  long cnt;

  cnt = 0;

  /* throw away the AUXEDGEs */
  for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
      if (AUXEDGE(MYEDGE(theLink)))
        if (DisposeEdge(theGrid,MYEDGE(theLink)))
          return (-1);
        else
          cnt++;

  if (cnt)
    GSTATUS(theGrid) = 1;

  return (cnt);
}
