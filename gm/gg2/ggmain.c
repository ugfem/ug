// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ggmain.c                                                      */
/*                                                                          */
/* Purpose:   central grid generator functions				                        */
/*                                                                          */
/* Author:    Wolfgang Hoffmann, Henrik Renz-Reichert	                    */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart, Germany										*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   08.03.94 begin, ug version 2.2                                */
/*                15.10.95 implemented in ug31                                  */
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

/* low */
#include "compiler.h"
#include "general.h"
#include "debug.h"
#include "heaps.h"
#include "misc.h"
#include "ugstruct.h"

/* dev */
#include "devices.h"
#ifdef __MPW32__
#include "MacGui.h"
#endif

/* dom */
#include "domain.h"

/* gm */
#include "gm.h"
#include "ggm.h"
#include "ugm.h"
#include "refine.h"
#include "algebra.h"
#include "evm.h"

/* TODO: hierarchy conflict due to UserRead, UserInterrupt */
#include "uginterface.h"


/* gg2 */
#include "ggm.h"
#include "ggaccel.h" /* header file for accelerator */

#include "ggmain.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAXNPOINTS 100          /* static heap size for the points in a square round the new FC          */
#define BARYCENTER      0.1666  /* factor for calculating the barycenter of the new created triangle */
#define ROOT3           1.7321  /* square root (3)                                                                                                       */
#define SMALLDOUBLE     1e-6

/* flags for different checks */
#define CHECKNEAR               (1<<0)
#define CHECKNBCUT              (1<<1)
#define CHECKINSIDE             (1<<2)
#define CHECKINTERSECT  (1<<3)
#define CHECKALL                (CHECKNEAR | CHECKNBCUT | CHECKINSIDE | CHECKINTERSECT)

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

typedef struct {
  INT n;
  INT SubdomainID;
  NODE *theNode[4];
  VERTEX *theVertex[4];
  ELEMENT *theNeighbour[4];
  INT Neighbourside[4];
  ELEMENT *thenewElement;
} ELEMENT_CONTEXT;

struct delaunay_context
{
  DOUBLE midp[2];
  DOUBLE radius;
};

typedef struct delaunay_context DELAUNAY_CONTEXT;

struct element_2d
{
  ELEMENT_CONTEXT elem_context;
  DELAUNAY_CONTEXT del_context;
  struct element_2d *succel;
};

typedef struct element_2d ELEMENT_2D;

struct mesh2d
{
  INT *nElements;
  ELEMENT_2D **start;
  ELEMENT_2D **Element;
};

typedef struct mesh2d MESH2D;

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

static char buffer[256];                /* general purpose text buffer				*/

static INT doupdate; /* if TRUE MakeElement calls UpdateDoc, PlotElement else*/
static INT dostep;              /* if TRUE MakeElement stops after each plot			*/
static INT doanimate;   /* if TRUE MakeElement plots the current element		*/
static INT doedge;      /* if TRUE MakeElement uses accelerator with edgetree*/
static INT doangle;     /* if TRUE MakeElement uses accelerator with angletree		*/
static INT doEdge;      /* if TRUE MakeElement uses smallest edge without accelerator*/
static INT doAngle;     /* if TRUE MakeElement uses smallest angle without accelerator*		*/
static INT doConstDel;  /* if TRUE Constrained Delaunay Triangulation*		*/

static INT ElemID;
static INT equilateral;
static INT plotfront;
static CoeffProcPtr nominal_h;
static DOUBLE searchradius,searchradius2;

static MG_GGDATA *myMGdata;

static GG_PARAM *myPars;

static INT SingleMode;


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

static INT IflObj;
static INT FlObj;
static INT FcObj;

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/*****************************************************************************/

static DOUBLE H_global;

static INT GlobalMeshsize (DOUBLE *in, DOUBLE *out)
{
  /* outvalue is just the constant global meshsize */
  out[0] = H_global;

  return(0);
}

static INT IsPointLeftOfFC (DOUBLE xP,DOUBLE yP, FRONTCOMP *theFC)
{
  VERTEX *theVertex;
  DOUBLE *cvect,xIn,yIn,xOut,yOut,xTest,yTest;
  DOUBLE cosIn,cosOut;
  INT sideIn,sideOut;

  cvect = CVECT(MYVERTEX(FRONTN(theFC)));
  theVertex = MYVERTEX(FRONTN(PREDFC(theFC)));
  xIn = cvect[0] - XC(theVertex);
  yIn = cvect[1] - YC(theVertex);
  theVertex = MYVERTEX(FRONTN(SUCCFC(theFC)));
  xOut = XC(theVertex) - cvect[0];
  yOut = YC(theVertex) - cvect[1];

  xTest = xP - cvect[0];
  yTest = yP - cvect[1];

  /* side.. = FALSE: xTest right of ..,
            = TRUE:  xTest left  of ..   */
  sideIn  = ((yTest*xIn-xTest*yIn)   > SMALLDOUBLE);
  sideOut = ((yTest*xOut-xTest*yOut) > SMALLDOUBLE);

  if (sideIn==sideOut)
  {
    /* in and out lie on the same side of test: both possible for decision */
    if (sideIn)
      return (YES);
    else
      return (NO);
  }
  else
  {
    /* in and out lie on different sides of test: take the closer one for decision */
    cosIn  = (xTest*xIn+yTest*yIn)/sqrt((xIn*xIn+yIn*yIn)*(xTest*xTest+yTest*yTest));
    cosOut = (xTest*xOut+yTest*yOut)/sqrt((xOut*xOut+yOut*yOut)*(xTest*xTest+yTest*yTest));

    if ((-cosIn) - cosOut > SMALLDOUBLE)
    {
      if (sideIn)
        return (YES);
      else
        return (NO);
    }
    else
    {
      if (sideOut)
        return (YES);
      else
        return (NO);
    }
  }
}

/****************************************************************************/
/*                                                                          */
/* Function:  DetermineOrientation                                                      */
/*                                                                          */
/* Purpose:   determines orientation of a FRONTLIST							*/
/*                                                                          */
/* Input:     front list										            */
/*                                                                          */
/* Output:    INT return code see header file                               */
/*                                                                          */
/****************************************************************************/

static INT DetermineOrientation (FRONTLIST *theFL)
{
  FRONTCOMP *theFC;
  DOUBLE np,dx1,dx2,dy1,dy2,*x0,*x1,*x2;
  double anglesum,arg;

  if (NFC(theFL)<3)
  {
    PrintErrorMessage('E',"DetermineOrientation","wrong orientation !!! ");
    return (1);
  }

  theFC = PREDFC(PREDFC(STARTFC(theFL)));
  x0 = CVECT(MYVERTEX(FRONTN(theFC)));
  theFC = SUCCFC(theFC);
  x1 = CVECT(MYVERTEX(FRONTN(theFC)));

  anglesum = 0.0;
  for (theFC=SUCCFC(theFC); theFC!=NULL; theFC=SUCCFC(theFC))
  {
    x2 = CVECT(MYVERTEX(FRONTN(theFC)));

    dx1 = x1[0] - x0[0];
    dy1 = x1[1] - x0[1];
    dx2 = x2[0] - x1[0];
    dy2 = x2[1] - x1[1];

    np = dx1*dy2-dy1*dx2;               /* x1 * normal(x2) */
    arg = (dx1*dx2+dy1*dy2)/sqrt((dx1*dx1+dy1*dy1)*(dx2*dx2+dy2*dy2));
    arg = MIN(1,arg);
    arg = MAX(-1,arg);
    anglesum += SIGNUM(np)*acos(arg);

    if (theFC==LASTFC(theFL))
      break;

    x0 = x1; x1 = x2;
  }

  FLORIENTATION(theFL) = (anglesum>0) ? MATHPOS : MATHNEG;

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  HandleSubdomain                                                           */
/*                                                                          */
/* Purpose:   find associated segments of the subdomain						*/
/*                                                                          */
/* Input:     independent front lists                                                               */
/*                                                                          */
/* Output:                                                                                                  */
/*                                                                          */
/****************************************************************************/

static INT HandleSubdomain (INDEPFRONTLIST *theIFL, NODE **Nodes,
                            INT SubdomainID, INT n, INT **node_id, INT MarkKey)
{
  MULTIGRID *theMG;
  NODE **NodeHandle;
  INT *SideFlag;
  FRONTLIST *theFL;
  INT i,j,id,cnt,orientation,oldside,gcnt;

  theMG = MYMG(MYGRID(theIFL));

  NodeHandle = (NODE **) GetTmpMem(theMG->theHeap,n*sizeof(NODE*),MarkKey);
  if (NodeHandle==NULL) return (1);
  SideFlag = (INT *) GetTmpMem(theMG->theHeap,n*sizeof(INT),MarkKey);
  if (SideFlag==NULL) return (1);
  for (i=0; i<n; i++) SideFlag[i] = 0;

  gcnt = 0;
  cnt = 0;
  orientation = 1;                      /* counterclockwise */
  oldside = -1;
  id = node_id[0][0];
  SideFlag[0] = 1;
  NodeHandle[0] = Nodes[id];
  while (1)
  {
    /* new node */
    PRINTDEBUG(dom,1,("  cnt nid %ld \n",ID(NodeHandle[cnt])));
    for (j=0; j<n; j++)
      if (NodeHandle[cnt] == Nodes[node_id[j][0]] && j!=oldside)
      {
        orientation = 1;
        break;
      }
      else if (NodeHandle[cnt]==Nodes[node_id[j][1]] && j!=oldside)
      {
        orientation = 0;
        break;
      }
    if (j==n) return(1);
    SideFlag[j] = 1;
    cnt++;

    /* check if closed */
    if (node_id[j][orientation] == id)
    {
      theFL = CreateFrontList(theIFL,SubdomainID);
      if (CreateFrontComp (theFL, NULL, cnt, NodeHandle)==NULL)
      {
        PrintErrorMessage('E',"HandleSubdomain","no storage for FC list");
        return (1);
      }
      if (DetermineOrientation (theFL)) return (2);
      gcnt += cnt;
      if (gcnt>=n) break;

      PRINTDEBUG(dom,1,("......\n"));

      for (i=0; i<n; i++) if (SideFlag[i]==0) break;
      if (i==n) return (1);
      cnt = 0;
      id = node_id[i][0];
      NodeHandle[0] = Nodes[id];
      SideFlag[i] = 1;
      j = -1;
    }
    else
      NodeHandle[cnt] = Nodes[node_id[j][orientation]];
    oldside = j;
  }

  PRINTDEBUG(dom,1,("------\n"));
  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  AssembleFrontLists                                                        */
/*                                                                          */
/* Purpose:   assembling of the associated frontlists						*/
/*                                                                          */
/* Input:                                                                                                                   */
/*                                                                          */
/* Output:                                                                                                  */
/*                                                                          */
/****************************************************************************/

static INT AssembleFrontLists (MULTIGRID *theMG, MESH *mesh, INT MarkKey)
{
  GRID *theGrid;
  NODE *theNode,**Nodes;
  INDEPFRONTLIST *theIFL;
  INT numOfSubdomains,SubdomainID,numOfBndP,i;

  theGrid = GRID_ON_LEVEL(theMG,0);
  numOfSubdomains = mesh->nSubDomains;
  numOfBndP = mesh->nBndP;

  Nodes = (NODE **) GetTmpMem(MGHEAP(theMG),numOfBndP*sizeof(NODE *),MarkKey);
  if (Nodes == NULL)
    return(1);

  for (i=0, theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {
    PRINTDEBUG(dom,1,("  nid %ld \n",ID(theNode)));
    if (V_BNDP(MYVERTEX(theNode)) != mesh->theBndPs[i])
      return(1);
    Nodes[i++] = theNode;
  }
  if (i != mesh->nBndP)
    return(1);

  theIFL = CreateIndepFrontList(theGrid);
  if (SingleMode<=0 || SingleMode>numOfSubdomains)
  {
    for (SubdomainID=1; SubdomainID<=numOfSubdomains; SubdomainID++)
      if (HandleSubdomain (theIFL,Nodes,SubdomainID, mesh->nSides[SubdomainID], mesh->Side_corner_ids[SubdomainID],MarkKey))
        return (1);
  }
  else
  {
    SubdomainID = SingleMode;
    if (HandleSubdomain (theIFL,Nodes,SubdomainID, mesh->nSides[SubdomainID], mesh->Side_corner_ids[SubdomainID],MarkKey))
      return (1);
  }

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  ChooseFC					                                    */
/*                                                                          */
/* Purpose:   create the start front components of the advancing front list	*/
/*                                                                          */
/* Input:                                                                                                               */
/*                                                                          */
/* Output:    front component for creating inner node                       */
/*                                                                          */
/****************************************************************************/

static FRONTCOMP *ChooseFCminangle (INDEPFRONTLIST *theIFL, FRONTLIST **myFLHandle)
{
  VERTEX *theVertex;
  FRONTCOMP *theFC,*bestFC;
  FRONTLIST *theFL,*myList;
  DOUBLE xc,yc,px,py,sx,sy,minangle,angle;

  /* return FC with minimal inner angle */

  minangle = MAX_C;
  bestFC = NULL;
  for (theFL=STARTFL(theIFL); theFL!=NULL; theFL=SUCCFL(theFL))
    for (theFC=STARTFC(theFL); theFC!=NULL; theFC=SUCCFC(theFC))
    {
      theVertex = MYVERTEX(FRONTN(theFC));
      xc = XC(theVertex);
      yc = YC(theVertex);
      theVertex = MYVERTEX(FRONTN(PREDFC(theFC)));
      px = xc - XC(theVertex);
      py = yc - YC(theVertex);
      theVertex = MYVERTEX(FRONTN(SUCCFC(theFC)));
      sx = XC(theVertex) - xc;
      sy = YC(theVertex) - yc;

      /* angle > 180 degrees? */
      if ((py*sx-px*sy) > SMALLDOUBLE)                                  /* check (py*sx-px*sy)>0 */
        if (theFC==LASTFC(theFL))
          break;
        else
          continue;

      angle = (px*sx+py*sy)/sqrt((sx*sx+sy*sy)*(px*px+py*py));
      if (angle<minangle)
      {
        minangle = angle;
        bestFC = theFC;
        myList = theFL;
      }

      if (theFC==LASTFC(theFL))
        break;
    }

  *myFLHandle = myList;
  return (bestFC);
}

static FRONTCOMP *ChooseFCminside (INDEPFRONTLIST *theIFL, FRONTLIST **myFLHandle)
{
  VERTEX *theVertex;
  FRONTCOMP *theFC,*shortestFC;
  FRONTLIST *theFL,*myList;
  DOUBLE xc,yc,sxc,syc,mindist2,dist2;

  /* return FC with shortest edge to SUCCFC */

  mindist2 = MAX_C;
  shortestFC = NULL;
  for (theFL=STARTFL(theIFL); theFL!=NULL; theFL=SUCCFL(theFL))
  {
    if ((theFC=LASTFC(theFL))==NULL)
      continue;

    theVertex = MYVERTEX(FRONTN(theFC));
    sxc = XC(theVertex);
    syc = YC(theVertex);

    for (theFC=PREDFC(theFC); theFC!=NULL; theFC=PREDFC(theFC))
    {
      theVertex = MYVERTEX(FRONTN(theFC));
      xc = XC(theVertex);
      yc = YC(theVertex);

      dist2 = (xc-sxc)*(xc-sxc)+(yc-syc)*(yc-syc);
      if (dist2<mindist2)
      {
        mindist2 = dist2;
        shortestFC = theFC;
        myList = theFL;
      }

      if (theFC==STARTFC(theFL))
        break;

      sxc = xc; syc = yc;
    }
  }

  *myFLHandle = myList;
  return (shortestFC);
}

/****************************************************************************/
/*                                                                          */
/* Function:  NearerNeighbour                                                           */
/*                                                                          */
/* Purpose:   checks the minimum dist. of inner nodes and surrounding nodes	*/
/*                                                                          */
/* Input:     front nodes and inner node							        */
/*																			*/
/*                                                                          */
/* Output:    node to be inserted in AFL                                                    */
/*                                                                          */
/****************************************************************************/

static FRONTCOMP *NearerNeighbour (FRONTCOMP *theFC, DOUBLE xt3, DOUBLE yt3)
{
  VERTEX *theVertex;
  DOUBLE xPQ,yPQ,xPR,yPR,xQ,yQ,xR,yR,xP2P,yP2P,xP2Q,yP2Q,xP1P,yP1P,xP1R,yP1R;
  DOUBLE succdistance2,preddistance2,succedgecos,prededgecos,succdenominator,preddenominator;


  /* calculation of the succ point distance to the inner node */
  theVertex = MYVERTEX(FRONTN(SUCCFC(SUCCFC(theFC))));
  xQ = XC(theVertex);
  yQ = YC(theVertex);
  xPQ = xQ - xt3;
  yPQ = yQ - yt3;

  succdistance2 = xPQ*xPQ+yPQ*yPQ;

  /* calculation of the diff. vector between P2 and innerNodeHandle for angle check */
  theVertex = MYVERTEX(FRONTN(SUCCFC(theFC)));
  xP2P = xt3 - XC(theVertex);
  yP2P = yt3 - YC(theVertex);
  xP2Q = xQ - XC(theVertex);
  yP2Q = yQ - YC(theVertex);

  succdenominator = sqrt((xP2Q*xP2Q+yP2Q*yP2Q)*(xP2P*xP2P+yP2P*yP2P));
  /* (fabs(succdenominator)<SMALLDOUBLE) would mean points are "identical" */

  succedgecos= ((xP2Q*xP2P)+(yP2Q*yP2P))/succdenominator;

  /* calculation of the left point distance of innerNodeHandle */
  theVertex = MYVERTEX(FRONTN(PREDFC(theFC)));
  xR = XC(theVertex);
  yR = YC(theVertex);
  xPR = xR - xt3;
  yPR = yR - yt3;

  preddistance2 = xPR*xPR+yPR*yPR;

  /* calculation of the diff. vector between P1 and innerNodeHandle for angle check */
  theVertex = MYVERTEX(FRONTN(theFC));
  xP1P = xt3 - XC(theVertex);
  yP1P = yt3 - YC(theVertex);
  xP1R = xR  - XC(theVertex);
  yP1R = yR  - YC(theVertex);

  preddenominator = sqrt((xP1P*xP1P+yP1P*yP1P)*(xP1R*xP1R+yP1R*yP1R));
  /* (fabs(preddenominator)<SMALLDOUBLE) would mean points are "identical" */

  prededgecos=((xP1P*xP1R)+(yP1P*yP1R))/preddenominator;



  /* find the correct front component proving the distance between the calc. point */
  /* and SUCCFC(SUCCFC(theFC)))... and also the angle between the two edges              */

  if ((succdistance2 < preddistance2) && (succdistance2 < searchradius2))
    return (SUCCFC(SUCCFC(theFC)));

  else if (preddistance2 < searchradius2)
    return (PREDFC(theFC));

  else if (succedgecos > myPars->CheckCos && succedgecos>prededgecos)
    return (SUCCFC(SUCCFC(theFC)));

  else if (prededgecos > myPars->CheckCos)
    return (PREDFC(theFC));

  else
    return (NULL);
}


/****************************************************************************/
/*                                                                          */
/* Function:  CutNeighbour		                                                */
/*                                                                          */
/* Purpose:   checks if the new triangle intersects with the neighbour front*/
/*                                                                          */
/* Input:     front node and inner node								        */
/*																			*/
/*                                                                          */
/* Output:    node to be inserted in AFL                                                    */
/*                                                                          */
/****************************************************************************/

static INT IsBetweenVectors (DOUBLE x1,DOUBLE y1,DOUBLE xb,DOUBLE yb,DOUBLE x2,DOUBLE y2)
{
  if (((xb*y1-yb*x1)<SMALL_C) && ((xb*y2-yb*x2)>SMALL_C))       /* right of 1 and left of 2, SMALL_C means numerical O */
    return (YES);
  else
    return (NO);
}

static FRONTCOMP *CutNeighbour (FRONTCOMP *theFC, DOUBLE xt[3], DOUBLE yt[3])
{
  VERTEX *theVertex;
  DOUBLE xQ,yQ,xR,yR,xP1Q,yP1Q,xP2R,yP2R,xP1P,yP1P,xP2P,yP2P,xP1P2,yP1P2;
  DOUBLE xP1minusP2,yP1minusP2;
  DOUBLE projection1,projection2;
  INT intersectPred,intersectSucc;

  /* determination of the vectors needed */
  theVertex = MYVERTEX(FRONTN(PREDFC(theFC))) ;
  xQ = XC(theVertex);
  yQ = YC(theVertex);
  xP1Q = xQ - xt[0];
  yP1Q = yQ - yt[0];

  theVertex = MYVERTEX(FRONTN(SUCCFC(SUCCFC(theFC))));
  xR = XC(theVertex);
  yR = YC(theVertex);
  xP2R = xR - xt[1];
  yP2R = yR - yt[1];

  xP1P = xt[2] - xt[0];
  yP1P = yt[2] - yt[0];

  xP2P = xt[2] - xt[1];
  yP2P = yt[2] - yt[1];

  xP1P2 = xt[1] - xt[0];
  yP1P2 = yt[1] - yt[0];

  /* vector P1Q between P1P2 and P1P */
  intersectPred = IsBetweenVectors(xP1P2,yP1P2,xP1Q,yP1Q,xP1P,yP1P);

  /* vector P2R between -P1P2 and P2P */
  intersectSucc = IsBetweenVectors(xP2P,yP2P,xP2R,yP2R,-xP1P2,-yP1P2);

  /* choose the right next front component */
  if (intersectPred && !intersectSucc)
    return (PREDFC(theFC));
  else if (!intersectPred && intersectSucc)
    return (SUCCFC(SUCCFC(theFC)));
  else if (intersectPred && intersectSucc)
  {
    /* third case both, P1P and P2P intersect with the AVF */

    xP1minusP2 = xt[0] - xt[1];
    yP1minusP2 = yt[0] - yt[1];

    projection1 = ABS(xP1P*xP1minusP2 + yP1P*yP1minusP2);
    projection2 = ABS(xP2R*(-xP1minusP2) + yP2R*(-yP1minusP2));

    if (projection1 - projection2 > SMALL_C)                                    /* check projection1>projection2 */
      return (SUCCFC(SUCCFC(theFC)));

    else
      return (PREDFC(theFC));

  }
  else return (NULL);
}

/*****************************************************************************/
/*                                                                           */
/* Function:  Point_In_Triangle		                                                 */
/*                                                                           */
/* Purpose:   check if a point of the AF lies inside the new created triangle*/
/*                                                                                           */
/*                                                                           */
/* Input:     coordinates of the AF point and new point						 */
/*																			 */
/*                                                                           */
/* Output:    return (YES) if the point is inside, else (NO)				 */
/*                                                                           */
/*****************************************************************************/

static INT Point_In_Triangle (DOUBLE pt[DIM], DOUBLE x[3], DOUBLE y[3])
{
  INT i,index;
  DOUBLE sx,sy,hlp;
  DOUBLE eps;
  DOUBLE x_eps[3], y_eps[3];

  eps = myPars->epsi;

  for (i=0; i<3; i++)
  {
    x_eps[i] = x[i];
    y_eps[i] = y[i];
  }

  /* Increasing of the search-triangle with eps to find nodes, which lie on its edges */

  /* find x_min */
  if (x_eps[0] <  x_eps[1])
    index = 0;
  else
    index = 1;

  if  (x_eps[index] < x_eps[2] )
    x_eps[index] -= eps;
  else
    x_eps[2] -= eps;

  /* find y_min */
  if (y_eps[0] <  y_eps[1])
    index = 0;
  else
    index = 1;

  if  (y_eps[index] < y_eps[2] )
    y_eps[index] -= eps;
  else
    y_eps[2] -= eps;

  /* find x_max */
  if (x_eps[0] >  x_eps[1])
    index = 0;
  else
    index = 1;

  if  (x_eps[index] > x_eps[2] )
    x_eps[index] += eps;
  else
    x_eps[2] += eps;

  /* find y_max */
  if (y_eps[0] >  y_eps[1])
    index = 0;
  else
    index = 1;

  if  (y_eps[index] > y_eps[2] )
    y_eps[index] += eps;
  else
    y_eps[2] += eps;

  for (i=0; i<3; i++)
  {

    sx = x_eps[(i+1)%3] - x_eps[i];
    sy = y_eps[(i+1)%3] - y_eps[i];

    hlp = ((sy*(pt[0]-x_eps[i])-sx*(pt[1]-y_eps[i]))/(sx*sx+sy*sy));
    if (hlp > SMALLDOUBLE)
      return (NO);
  }
  return (YES);
}

static INT Point_In_Circle (DOUBLE pt[DIM], DOUBLE x, DOUBLE y, DOUBLE searchrad2)
{
  DOUBLE dx,dy;

  dx = pt[0] - x;
  dy = pt[1] - y;

  if (searchrad2-(dx*dx+dy*dy)>SMALLDOUBLE )        /* check of searchrad2>(dx*dx+dy*dy) */
    return (YES);
  else
    return (NO);
}

/****************************************************************************/
/*                                                                          */
/* Function:  AssemblePointsInTriangleOrCircle                              */
/*                                                                          */
/* Purpose:   assemble the adresses of the FC within the square				*/
/*                                                                                          */
/*                                                                          */
/* Input:     Independent Front List										*/
/*																			*/
/*                                                                          */
/* Output:   found points within the triangle                                                           */
/*                                                                          */
/****************************************************************************/

static INT AssemblePointsInTriangleOrCircle (INDEPFRONTLIST *theIFL,
                                             FRONTCOMP *thefoundPoints[MAXNPOINTS], DOUBLE x[3], DOUBLE y[3], DOUBLE radius)
{
  VERTEX *theVertex;
  FRONTLIST *theFL;
  FRONTCOMP *thecompFC;
  DOUBLE rad2;
  INT foundpoints;

  rad2 = radius*radius;

  for (foundpoints=0,theFL=LASTFL(theIFL); theFL != NULL; theFL=PREDFL(theFL))
    for (thecompFC=STARTFC(theFL); thecompFC != NULL; thecompFC=SUCCFC(thecompFC))
    {
      theVertex = MYVERTEX(FRONTN(thecompFC));
      if (    Point_In_Triangle       (CVECT(theVertex),x,y)
              ||      Point_In_Circle (CVECT(theVertex),x[2],y[2],rad2))
      {
        if ( ++foundpoints >= MAXNPOINTS)
        {
          PrintErrorMessage('F',"AssemblePointsInTriangleOrCircle","no. of found points overflows the internal heap size! ");
          return (-1);
        }
        thefoundPoints[foundpoints-1] = thecompFC;
      }
      if (thecompFC==LASTFC(theFL))
        break;
    }

  return (foundpoints);

}

/****************************************************************************/
/*                                                                          */
/* Function:  FCinsideOrNear		                                        */
/*                                                                          */
/* Purpose:   find the nearest front component to the new created FC		*/
/*                                                                                          */
/*                                                                          */
/* Input:     FC within the triangle and coord. of the triangle                         */
/*																			*/
/*                                                                          */
/* Output:   the new FC or the old point                                                                        */
/*                                                                          */
/****************************************************************************/

static INT FCinsideOrNear (INDEPFRONTLIST *theIFL,FRONTCOMP *theFC,FRONTCOMP *theProposedFC, FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS],
                           DOUBLE xt[3],DOUBLE yt[3], FRONTCOMP **newFChandle, INT checkoptions)
{
  VERTEX *theVertex;
  NODE *BaseNode1,*BaseNode2,*PropNode;
  FRONTCOMP *thefoundPoints[MAXNPOINTS],*thenewFC;
  DOUBLE normalx,normaly,relx,rely,xM,yM;
  DOUBLE deltaP,minimaldist;
  INT i,foundpoints;

  /* variables for epsilon begin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  DOUBLE x_eps[3],y_eps[3];
  DOUBLE eps,mz,mn,m,xoffs,yoffs,hilfx,hilfy;
  const DOUBLE smallcoordpos = 10e-6;
  const DOUBLE smallcoordneg = -10e-6;


  /* triangle += epsilon begin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  for (i=0; i<3; i++)
  {
    x_eps[i] = xt[i];
    y_eps[i] = yt[i];
  }

  hilfx = x_eps[2];
  hilfy = y_eps[2];

  eps = myPars->epsi;


  /***********************/
  /* first edge: "basis" */
  /***********************/
  mz = yt[1] - yt[0];
  mn = xt[1] - xt[0];

  if ((mz<smallcoordpos)&&(mz>smallcoordneg))       /* case "-" */
    if (xt[1]>xt[0])             /* subcase "->" */
    {
      x_eps[0] -= eps;
      x_eps[1] += eps;
    }
    else             /* subcase "<-" */
    {
      x_eps[0] += eps;
      x_eps[1] -= eps;
    }
  else if ((mn<smallcoordpos)&&(mn>smallcoordneg))       /* case "|" */
    if (yt[1]>yt[0])             /* subcase "|u" */
    {
      y_eps[0] -= eps;
      y_eps[1] += eps;
    }
    else             /* subcase "|d" */
    {
      y_eps[0] += eps;
      y_eps[1] -= eps;
    }
  else
  {
    m = mz/mn;
    xoffs = eps/sqrt(1+m*m);
    yoffs = m*xoffs;
    if(yoffs<0.0) yoffs = -yoffs;

    if (yt[1]>yt[0])             /* case "/u" or "\u" */
    {
      y_eps[0] -= yoffs;
      y_eps[1] += yoffs;
    }
    else             /* case "/d" or "\d" */
    {
      y_eps[0] += yoffs;
      y_eps[1] -= yoffs;
    }

    if (xt[1]>xt[0])             /* case "/u" or "\d" */
    {
      x_eps[0] -= xoffs;
      x_eps[1] += xoffs;
    }
    else             /* case "/d" or "\u" */
    {
      x_eps[0] += xoffs;
      x_eps[1] -= xoffs;
    }

  }

  /***********************/
  /* second edge: "left" */
  /***********************/
  mz = yt[2] - yt[0];
  mn = xt[2] - xt[0];

  if ((mz<smallcoordpos)&&(mz>smallcoordneg))       /* case "-" */
    if (xt[2]>xt[0])             /* subcase "->" */
    {
      x_eps[2] += eps;
    }
    else             /* subcase "<-" */
    {
      x_eps[2] -= eps;
    }
  else if ((mn<smallcoordpos)&&(mn>smallcoordneg))       /* case "|" */
    if (yt[2]>yt[0])             /* subcase "|u" */
    {
      y_eps[2] += eps;
    }
    else             /* subcase "|d" */
    {
      y_eps[2] -= eps;
    }
  else
  {
    m = mz/mn;
    xoffs = eps/sqrt(1+m*m);
    yoffs = m*xoffs;
    if(yoffs<0.0) yoffs = -yoffs;

    if (yt[2]>yt[0])             /* case "/u" or "\u" */
    {
      y_eps[2] += yoffs;
    }
    else             /* case "/d" or "\d" */
    {
      y_eps[2] -= yoffs;
    }

    if (xt[2]>xt[0])             /* case "/u" or "\d" */
    {
      x_eps[2] += xoffs;
    }
    else             /* case "/d" or "\u" */
    {
      x_eps[2] -= xoffs;
    }

  }



  /***********************/
  /* third edge: "right" */
  /***********************/
  mz = yt[2] - yt[1];
  mn = xt[2] - xt[1];

  if ((mz<smallcoordpos)&&(mz>smallcoordneg))       /* case "-" */
    if (xt[2]>xt[1])             /* subcase "->" */
    {
      hilfx += eps;
    }
    else             /* subcase "<-" */
    {
      hilfx -= eps;
    }
  else if ((mn<smallcoordpos)&&(mn>smallcoordneg))       /* case "|" */
    if (yt[2]>yt[1])             /* subcase "|u" */
    {
      hilfy += eps;
    }
    else             /* subcase "|d" */
    {
      hilfy -= eps;
    }
  else
  {
    m = mz/mn;
    xoffs = eps/sqrt(1+m*m);
    yoffs = m*xoffs;
    if(yoffs<0.0) yoffs = -yoffs;

    if (yt[2]>yt[1])             /* case "/u" or "\u" */
    {
      hilfy += yoffs;
    }
    else             /* case "/d" or "\d" */
    {
      hilfy -= yoffs;
    }

    if (xt[2]>xt[1])             /* case "/u" or "\d" */
    {
      hilfx += xoffs;
    }
    else             /* case "/d" or "\u" */
    {
      hilfx -= xoffs;
    }

  }


  /* the third point of the epsilon-triangle: */
  y_eps[2] = (y_eps[2] + hilfy)/2.0;
  x_eps[2] = (x_eps[2] + hilfx)/2.0;

  /* triangle += epsilon finished !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  if (doAngle || doEdge)
    foundpoints = AccelFCTreeSearch (theIFL,thefoundPoints,theIntersectfoundPoints,xt,yt,searchradius);
  else
    foundpoints = AssemblePointsInTriangleOrCircle (theIFL,thefoundPoints,xt,yt,searchradius);

  *newFChandle = thenewFC = NULL;

  if (foundpoints<0)
    return (1);

  if (foundpoints==0)
    return (0);

  minimaldist = MAX_D;
  normalx =  yt[0] - yt[1];
  normaly = -xt[0] + xt[1];
  xM = 0.5*(xt[0] + xt[1]);
  yM = 0.5*(yt[0] + yt[1]);

  BaseNode1 = FRONTN(theFC);
  BaseNode2 = FRONTN(SUCCFC(theFC));
  PropNode  = (theProposedFC==NULL) ? NULL : FRONTN(theProposedFC);

  for (i=0; i<foundpoints; i++)
  {
    if ((FRONTN(thefoundPoints[i])==BaseNode1) || (FRONTN(thefoundPoints[i])==BaseNode2)
        || (FRONTN(thefoundPoints[i])==PropNode))
      continue;

    /* determination of the distance between the found points and the triangle base */
    theVertex = MYVERTEX(FRONTN(thefoundPoints[i]));
    relx = XC(theVertex) - xt[0];
    rely = YC(theVertex) - yt[0];

    /* are we on the correct side of the base? */
    deltaP = normalx*relx + normaly*rely;

    if (deltaP<0)
      continue;

    if (deltaP -minimaldist < SMALL_C)               /* check deltaP < minimaldist */
      /* take this point only if (xM,yM) lies left of the other front */
      if (IsPointLeftOfFC(xM,yM,thefoundPoints[i]))
      {
        minimaldist = deltaP;
        thenewFC = thefoundPoints[i];
      }
  }

  if (thenewFC!=NULL)
    if (!(checkoptions&CHECKNEAR))                      /* check only inside */
      if (!Point_In_Triangle (CVECT(MYVERTEX(FRONTN(thenewFC))),xt,yt))
        thenewFC = NULL;

  *newFChandle = thenewFC;

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  FrontLineIntersection			                                */
/*                                                                          */
/* Purpose:   find the intersection of the AF with the new triangle			*/
/*                                                                                          */
/*                                                                          */
/* Input:     coordinates of the new point and of the egde points			*/
/*																			*/
/*                                                                          */
/* Output:   the new FC or NULL                                                                                         */
/*                                                                          */
/****************************************************************************/

static FRONTCOMP *FrontLineIntersection (INDEPFRONTLIST *theIFL,FRONTCOMP *theFC,FRONTCOMP *theProposedFC,
                                         FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS], DOUBLE xt[3],DOUBLE yt[3])
{
  NODE *theNode,*BaseNode;
  VERTEX *theVertex;
  FRONTLIST *theFL;
  FRONTCOMP *thecompFC,*thedistFC1;
  DOUBLE xM,yM,xR,yR,xQ,yQ,lambda1,lambda2;
  DOUBLE xMQ,yMQ,xMR,yMR,denominator;
  DOUBLE lambdacomp = MAX_C;
  INT j = 0;

  /* variables for epsilon begin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  INT i;
  DOUBLE x_eps[3],y_eps[3];
  DOUBLE eps,mz,mn,m,xoffs,yoffs,hilfx,hilfy;
  const DOUBLE smallcoordpos = 10e-6;
  const DOUBLE smallcoordneg = -10e-6;



  /* triangle += epsilon begin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  for (i=0; i<3; i++)
  {
    x_eps[i] = xt[i];
    y_eps[i] = yt[i];
  }

  hilfx = x_eps[2];
  hilfy = y_eps[2];

  eps = myPars->epsi;


  /***********************/
  /* first edge: "basis" */
  /***********************/
  mz = yt[1] - yt[0];
  mn = xt[1] - xt[0];

  if ((mz<smallcoordpos)&&(mz>smallcoordneg))       /* case "-" */
    if (xt[1]>xt[0])             /* subcase "->" */
    {
      x_eps[0] -= eps;
      x_eps[1] += eps;
    }
    else             /* subcase "<-" */
    {
      x_eps[0] += eps;
      x_eps[1] -= eps;
    }
  else if ((mn<smallcoordpos)&&(mn>smallcoordneg))       /* case "|" */
    if (yt[1]>yt[0])             /* subcase "|u" */
    {
      y_eps[0] -= eps;
      y_eps[1] += eps;
    }
    else             /* subcase "|d" */
    {
      y_eps[0] += eps;
      y_eps[1] -= eps;
    }
  else
  {
    m = mz/mn;
    xoffs = eps/sqrt(1+m*m);
    yoffs = m*xoffs;
    if(yoffs<0.0) yoffs = -yoffs;

    if (yt[1]>yt[0])             /* case "/u" or "\u" */
    {
      y_eps[0] -= yoffs;
      y_eps[1] += yoffs;
    }
    else             /* case "/d" or "\d" */
    {
      y_eps[0] += yoffs;
      y_eps[1] -= yoffs;
    }

    if (xt[1]>xt[0])             /* case "/u" or "\d" */
    {
      x_eps[0] -= xoffs;
      x_eps[1] += xoffs;
    }
    else             /* case "/d" or "\u" */
    {
      x_eps[0] += xoffs;
      x_eps[1] -= xoffs;
    }

  }

  /***********************/
  /* second edge: "left" */
  /***********************/
  mz = yt[2] - yt[0];
  mn = xt[2] - xt[0];

  if ((mz<smallcoordpos)&&(mz>smallcoordneg))       /* case "-" */
    if (xt[2]>xt[0])             /* subcase "->" */
    {
      x_eps[2] += eps;
    }
    else             /* subcase "<-" */
    {
      x_eps[2] -= eps;
    }
  else if ((mn<smallcoordpos)&&(mn>smallcoordneg))       /* case "|" */
    if (yt[2]>yt[0])             /* subcase "|u" */
    {
      y_eps[2] += eps;
    }
    else             /* subcase "|d" */
    {
      y_eps[2] -= eps;
    }
  else
  {
    m = mz/mn;
    xoffs = eps/sqrt(1+m*m);
    yoffs = m*xoffs;
    if(yoffs<0.0) yoffs = -yoffs;

    if (yt[2]>yt[0])             /* case "/u" or "\u" */
    {
      y_eps[2] += yoffs;
    }
    else             /* case "/d" or "\d" */
    {
      y_eps[2] -= yoffs;
    }

    if (xt[2]>xt[0])             /* case "/u" or "\d" */
    {
      x_eps[2] += xoffs;
    }
    else             /* case "/d" or "\u" */
    {
      x_eps[2] -= xoffs;
    }

  }



  /***********************/
  /* third edge: "right" */
  /***********************/
  mz = yt[2] - yt[1];
  mn = xt[2] - xt[1];

  if ((mz<smallcoordpos)&&(mz>smallcoordneg))       /* case "-" */
    if (xt[2]>xt[1])             /* subcase "->" */
    {
      hilfx += eps;
    }
    else             /* subcase "<-" */
    {
      hilfx -= eps;
    }
  else if ((mn<smallcoordpos)&&(mn>smallcoordneg))       /* case "|" */
    if (yt[2]>yt[1])             /* subcase "|u" */
    {
      hilfy += eps;
    }
    else             /* subcase "|d" */
    {
      hilfy -= eps;
    }
  else
  {
    m = mz/mn;
    xoffs = eps/sqrt(1+m*m);
    yoffs = m*xoffs;
    if(yoffs<0.0) yoffs = -yoffs;

    if (yt[2]>yt[1])             /* case "/u" or "\u" */
    {
      hilfy += yoffs;
    }
    else             /* case "/d" or "\d" */
    {
      hilfy -= yoffs;
    }

    if (xt[2]>xt[1])             /* case "/u" or "\d" */
    {
      hilfx += xoffs;
    }
    else             /* case "/d" or "\u" */
    {
      hilfx -= xoffs;
    }

  }


  /* the third point of the epsilon-triangle: */
  y_eps[2] = (y_eps[2] + hilfy)/2.0;
  x_eps[2] = (x_eps[2] + hilfx)/2.0;

  /* triangle += epsilon finished !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  /* calculation of the vector (xM,yM) to the mid of the edge between P1 and P2 */
  xM = 0.5*(xt[0] + xt[1]);
  yM = 0.5*(yt[0] + yt[1]);

  thedistFC1 = NULL;

  BaseNode = FRONTN(theFC);
  if (doAngle || doEdge)
    while ( theIntersectfoundPoints[j] != NULL )
    {
      if (FRONTN(theIntersectfoundPoints[j])==BaseNode)
      {
        j++;
        continue;                               /* skip the base of our triangle */
      }

      /* calculation of the coordinates of the FC's to proove */
      theVertex = MYVERTEX(FRONTN(theIntersectfoundPoints[j]));
      xR = XC(theVertex);
      yR = YC(theVertex);
      theVertex = MYVERTEX(FRONTN(SUCCFC(theIntersectfoundPoints[j])));
      xQ = XC(theVertex);
      yQ = YC(theVertex);

      /* check if the AF intersects with the triangle */
      denominator = (yR-yQ)*(xt[2]-xM)-(xR-xQ)*(yt[2]-yM);
      if (fabs(denominator)<SMALLDOUBLE)
      {
        j++;
        continue;                               /* lines are parallel */
      }

      if (fabs(xt[2]-xM)<SMALLDOUBLE)
      {
        lambda2 = (xR-xM)/(xR-xQ);
        lambda1 = (yR-yM - lambda2*(yR-yQ)) / (yt[2]-yM);
      }
      else
      {
        lambda2=((yR-yM)*(xt[2]-xM)-(xR-xM)*(yt[2]-yM))/denominator;
        lambda1=(xR-xM-lambda2*(xR-xQ))/(xt[2]-xM);
      }

      /* find the nearest cutting line */
      /* lambda1 parametrizes the cut on the triangle height
         lambda2 parametrizes the cut on the edge from thecompFC to its succ */
      if (lambda1 <= 1.15 && lambda1 >= 0 && lambda2 <= 1.0 && lambda2 >= 0.0)
      {
        if (lambda1 - lambdacomp < SMALLDOUBLE)                           /* check of lambda1 < lambdacomp */
          /* take this edge only if (xM,yM) lies left of the other front */
          if (IsPointLeftOfFC(xM,yM,theIntersectfoundPoints[j]))
          {
            lambdacomp = lambda1;
            thedistFC1 = theIntersectfoundPoints[j];
          }
      }
      ++j;
    }            /* of while */

  else        /* the non accelerated case */

    for (theFL=LASTFL(theIFL); theFL != NULL; theFL=PREDFL(theFL))
    {
      for (thecompFC=STARTFC(theFL); thecompFC != NULL; thecompFC=SUCCFC(thecompFC))
      {
        if (FRONTN(thecompFC)==BaseNode)
          if (thecompFC==LASTFC(theFL))
            break;
          else
            continue;                                           /* skip the base of our triangle */

        /* calculation of the coordinates of the FC's to proove */
        theVertex = MYVERTEX(FRONTN(SUCCFC(thecompFC)));
        xR = XC(theVertex);
        yR = YC(theVertex);
        theVertex = MYVERTEX(FRONTN(thecompFC));
        xQ = XC(theVertex);
        yQ = YC(theVertex);

        /* check if the AF intersects with the triangle */
        denominator = (yR-yQ)*(xt[2]-xM)-(xR-xQ)*(yt[2]-yM);
        if (fabs(denominator)<SMALLDOUBLE)
          if (thecompFC==LASTFC(theFL))
            break;
          else
            continue;                                           /* lines are parallel */

        if (fabs(xt[2]-xM)<SMALLDOUBLE)
        {
          lambda2 = (xR-xM)/(xR-xQ);
          lambda1 = (yR-yM - lambda2*(yR-yQ)) / (yt[2]-yM);
        }
        else
        {
          lambda2=((yR-yM)*(xt[2]-xM)-(xR-xM)*(yt[2]-yM))/denominator;
          lambda1=(xR-xM-lambda2*(xR-xQ))/(xt[2]-xM);
        }

        /* find the nearest cutting line */
        /* lambda1 parametrizes the cut on the triangle height
           lambda2 parametrizes the cut on the edge from thecompFC to its succ */
        if (lambda1 <= 1.15 && lambda1 >= 0 && lambda2 <= 1.0 && lambda2 >= 0.0)
        {
          if (lambda1 - lambdacomp < SMALLDOUBLE)                         /* check of lambda1 < lambdacomp */
            /* take this edge only if (xM,yM) lies left of the other front */
            if (IsPointLeftOfFC(xM,yM,thecompFC))
            {
              lambdacomp = lambda1;
              thedistFC1 = thecompFC;
            }
        }
        if (thecompFC==LASTFC(theFL))
          break;
      }
    }            /* else for */

  if (thedistFC1 != NULL )
  {
    theNode = FRONTN(thedistFC1);

    /* is thedistFC1 a possible alternative */
    if (theNode==FRONTN(PREDFC(theFC)))
      return (NULL);

    if (theNode==FRONTN(theFC))
      return (NULL);

    if (theNode==FRONTN(SUCCFC(theFC)))
      return (NULL);

    if (theProposedFC!=NULL)
    {
      if (theNode==FRONTN(theProposedFC))
        return (NULL);

      /* exclude all edges beginning at theProposedFC */
      if (theNode==FRONTN(PREDFC(theProposedFC)))
        return (NULL);

      /* exclude all edges ending at theProposedFC */
      if (FRONTN(SUCCFC(thedistFC1))==FRONTN(theProposedFC))
        return (NULL);
    }

    /* calculation of the coordinates of the FC */
    /* diagonal comparison of the two possible points */
    theVertex = MYVERTEX(FRONTN(SUCCFC(thedistFC1)));
    xR = XC(theVertex);
    yR = YC(theVertex);
    theVertex = MYVERTEX(FRONTN(thedistFC1));
    xQ = XC(theVertex);
    yQ = YC(theVertex);

    /* calculation of the nearest distance of the two found points to the midpoint of the base */
    xMQ = (xQ - xM); yMQ = (yQ - yM);
    xMR = (xR - xM); yMR = (yR - yM);

    /* find the FC for creating the triangle */

    if ((xMQ*xMQ+yMQ*yMQ)-(xMR*xMR+yMR*yMR) > SMALLDOUBLE)              /* check (xMQ*xMQ+yMQ*yMQ) > (xMR*xMR+yMR*yMR) */
    {
      /* lies R left of P1P2? */
      if ((xMR*(yt[0]-yt[1])-yMR*(xt[0]-xt[1])) > SMALLDOUBLE)
        return (SUCCFC(thedistFC1));                                    /* return R */
      else
        return (thedistFC1);                                                    /* return Q */
    }
    else
    {
      /* lies Q left of P1P2? */
      if ((xMQ*(yt[0]-yt[1])-yMQ*(xt[0]-xt[1])) > SMALLDOUBLE)
        return (thedistFC1);                                                    /* return R */
      else
        return (SUCCFC(thedistFC1));                                    /* return Q */
    }
  }
  else
    return (NULL);

}

/****************************************************************************/
/*                                                                          */
/* Function:  ContainedIn					                                */
/*                                                                          */
/* Purpose:   check wether FL2 is contained in FL1				                        */
/*                                                                          */
/* Input:     FL1,FL2														*/
/*                                                                          */
/* Output:    YES or NO										                                */
/*                                                                          */
/****************************************************************************/

static INT ContainedIn (FRONTLIST *FL1, FRONTLIST *FL2)
{
  FRONTCOMP *testFC,*theFC;
  DOUBLE dx,dy,y1cut,y2cut,x1cut,*testc,*thec,*thesc;
  INT ncut;

  /* check wether one arbitrary point of the second list is contained in the first one */

  testFC = STARTFC(FL2);
  testc = CVECT(MYVERTEX(FRONTN(testFC)));

  ncut = 0;
  for (theFC=STARTFC(FL1); theFC!=NULL; theFC=SUCCFC(theFC))
  {
    /* cut with right half line from testFC to infinity? */
    thec  = CVECT(MYVERTEX(FRONTN(theFC)));
    thesc = CVECT(MYVERTEX(FRONTN(SUCCFC(theFC))));
    y1cut = thec[1]  - testc[1];
    y2cut = thesc[1] - testc[1];

    if ((y1cut*y2cut) < 0)
    {
      /* thec and thesc lie on different sides, do they cut? */
      dx = thesc[0] - thec[0];
      dy = thesc[1] - thec[1];

      if (fabs(dy)>SMALLDOUBLE)
      {
        x1cut = y1cut * dx/dy;

        if ((x1cut+thec[0])-(testc[0]) > SMALLDOUBLE )                          /* check (x1cut+thec[0]) > testc[0] */
          ncut++;
      }
      else
      if (thec[0]-testc[0]>SMALLDOUBLE)
        ncut++;
    }

    if (theFC==LASTFC(FL1))
      break;
  }

  if (FLORIENTATION(FL1)==MATHPOS)
  {
    if (ODD(ncut))
      return (YES);
    else
      return (NO);
  }
  else
  {
    if (EVEN(ncut))
      return (YES);
    else
      return (NO);
  }
}

/****************************************************************************/
/*                                                                          */
/* Function:  RedistributeFLs				                                */
/*                                                                          */
/* Purpose:   Redistribute FLs on the two IFLs					                        */
/*                                                                          */
/* Input:     myList,thenewFL (resulting from merging process)				*/
/*                                                                          */
/* Output:    update of the AF								                                */
/*                                                                          */
/****************************************************************************/

static INT RedistributeFLs (FRONTLIST *myList,FRONTLIST *thenewFL)
{
  INDEPFRONTLIST *myIFL,*thenewIFL;
  FRONTLIST *theFL,*theloopFC;

  myIFL     = MYIFL(myList);
  thenewIFL = MYIFL(thenewFL);

  /* take smaller list for decisions */
  for (theFL=STARTFL(myIFL); theFL!=NULL; theFL=theloopFC)
  {
    theloopFC = SUCCFL(theFL);

    if (theFL==myList)
      continue;

    if (NFC(myList)<NFC(thenewFL))
    {
      if (ContainedIn(myList,theFL))
        continue;
    }
    else
    {
      if (!ContainedIn(thenewFL,theFL))
        continue;
    }

    /* remove theFL from myIFL */
    if (PREDFL(theFL)!=NULL)
      SUCCFL(PREDFL(theFL)) = SUCCFL(theFL);
    else
      STARTFL(myIFL) = SUCCFL(theFL);
    if (SUCCFL(theFL)!=NULL)
      PREDFL(SUCCFL(theFL)) = PREDFL(theFL);
    if (LASTFL(myIFL)==theFL)
      LASTFL(myIFL) = PREDFL(theFL);
    NFL(myIFL)--;

    /* put theFL into thenewIFL */
    MYIFL(theFL) = thenewIFL;
    SUCCFL(theFL) = STARTFL(thenewIFL);
    PREDFL(SUCCFL(theFL)) = theFL;
    PREDFL(theFL) = NULL;
    STARTFL(thenewIFL) = theFL;
    NFL(thenewIFL)++;
  }

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  FrontLineUpDate				                                */
/*                                                                          */
/* Purpose:   update of the AF                                                                                          */
/*                                                                          */
/* Input:     FC and newFC													*/
/*                                                                          */
/* Output:    update of the AF								                                */
/*                                                                          */
/****************************************************************************/

static INT FrontLineUpDate (GRID *theGrid,INDEPFRONTLIST *theIFL,FRONTLIST *myList,FRONTCOMP *theFC,FRONTCOMP *thenewFC, int* FlgForAccel, FRONTCOMP *the_old_succ, FRONTCOMP **disp_FC, FRONTLIST **disp_FL)
{

  MULTIGRID *theMG;
  INDEPFRONTLIST *thenewIFL;
  FRONTLIST *myfoundFL,*thenewFL,*theFL;
  FRONTCOMP  *thenewdoubledFC,*thelistFC,*thesuccFC,*thememsuccFC;
  INT numofFC,totalnumofFC;
  INT n = 3;

  theMG = MYMG(theGrid);

  *FlgForAccel = NORMALCASE;

  /* if FL consists of 3 Fcs, therefore accelerator has to delete them. */
  if (SUCCFC(SUCCFC(theFC))==PREDFC(theFC))
  {
    *FlgForAccel = FINALCASE;
    *disp_FC = NULL;
    *disp_FL = myList;
    return (0);

  }       /* last element of an IFL */


  /* handle simple cases */
  else if (SUCCFC(theFC)==thenewFC)
  {
    *FlgForAccel = NORMALCASE;


    return (0);

  }       /* end of handle simple case */


  else if (SUCCFC(SUCCFC(theFC))==thenewFC)
  {
    /* right */
    *FlgForAccel = RIGHTNEIGHBOURCASE;
    /* Cut Neighbour or Nearer Neighbour  right*/

    /* if FL consists of 3 Fcs, therefore accelerator has to delete them. */
    if (SUCCFC(SUCCFC(theFC))==PREDFC(theFC))
    {
      *FlgForAccel = FINALCASE;

      *disp_FC = SUCCFC(theFC);
      *disp_FL = myList;
      /* no more here: DisposeFrontComp (myList, SUCCFC(theFC)); */
      return (0);

    }             /* end of final case in rightneighbourcase */



    /* the rightneighbourcase - Continue*/

    /* In this case SUCCFC(SUCCFC(theFC)) <--> thenewFC */


    *disp_FC = SUCCFC(theFC);
    *disp_FL = myList;
    /* no more here:DisposeFrontComp (myList, SUCCFC(theFC));*/
    return (0);
  }       /* the end of rightneighbourcase */


  else if (PREDFC(theFC)==thenewFC)
  {
    /* left */
    *FlgForAccel = LEFTNEIGHBOURCASE;
    /* Cut Neighbour or Nearer Neighbour  left*/

    /* if FL consists of 3 Fcs, therefore accelerator has to delete them. */
    if (SUCCFC(SUCCFC(theFC))==PREDFC(theFC))
    {
      *FlgForAccel = FINALCASE;

      *disp_FC = SUCCFC(theFC);
      *disp_FL = myList;
      /*no more here: DisposeFrontComp (myList, SUCCFC(theFC));*/
      return (0);

    }               /* end of final case in leftneighbourcase */

    /* the leftneighbourcase - Continue*/



    *disp_FC = theFC;
    *disp_FL = myList;
    /*no more here: DisposeFrontComp (myList, theFC);*/
    return (0);
  }



  /* find the front List associated to the thenewFC */
  for (theFL=LASTFL(theIFL); theFL!=NULL; theFL=PREDFL(theFL))
  {
    for (thelistFC=STARTFC(theFL); thelistFC!=NULL; thelistFC=SUCCFC(thelistFC))
    {
      if (thenewFC==thelistFC)
      {
        myfoundFL = theFL;
        theFL = STARTFL(theIFL);
        break;
      }

      if (thelistFC==LASTFC(theFL))
        break;
    }
  }

#if (defined __MPW32__ || defined __MWCW__)
  if (plotfront)
  {
    /* plot the old ADVFront */



    UserWrite("front before being updated <RETURN>\n");
    UserRead(buffer);
  }
#endif

  /* duplicate thenewFC in myfoundFL */

  /* if FL consists of 3 Fcs, therefore accelerator has to delete them. */

  *FlgForAccel = ININTERCASE;

  thenewdoubledFC = CreateFrontComp (myfoundFL,thenewFC,1,&FRONTN(thenewFC));

  if (thenewdoubledFC==NULL)
  {
    PrintErrorMessage('E',"FrontLineUpDate","no storage for new FC");
    return (1);
  }

  /* necessary  thenewdoubledFC muss Nb, Sd und Nbsd von thenewFC uebernehmen.  */
  FCNGB(thenewdoubledFC)  = FCNGB(thenewFC);
  FCNGBS(thenewdoubledFC) = FCNGBS(thenewFC);


  /* update of the advancing front */

  /* first case, the FCs are in the same front list */
  if (myfoundFL==myList)
  {
    /* extract the new list beginning with thenewdoubledFC and ending with
       theFC from myList */

    thesuccFC = SUCCFC(theFC);

    thememsuccFC = SUCCFC(theFC);


    SUCCFC(theFC) = thenewdoubledFC;
    PREDFC(thenewdoubledFC) = theFC;
    SUCCFC(thenewFC) = thesuccFC;
    PREDFC(thesuccFC) = thenewFC;

    thenewIFL = CreateIndepFrontList (theGrid);
    if (thenewIFL==NULL)
    {
      PrintErrorMessage('E',"FrontLineUpDate","no storage for new IFL");
      return (1);
    }
    thenewFL = CreateFrontList (thenewIFL,myList->SubdomainID);
    if (thenewFL==NULL)
    {
      PrintErrorMessage('E',"FrontLineUpDate","no storage for new FL");
      return (1);
    }

    STARTFC(thenewFL) = thenewdoubledFC;
    LASTFC(thenewFL)  = theFC;
    /* tell all FCs of the new FL that they belong to the new list newFL */
    for ( thelistFC =  STARTFC(thenewFL); thelistFC != NULL; thelistFC = SUCCFC(thelistFC) )
    {
      MYFL(thelistFC) = thenewFL;

      if ( thelistFC == LASTFC(thenewFL) ) break;
    }

    STARTFC(myList) = thesuccFC;
    LASTFC(myList)  = thenewFC;

    /* update of the number of the front components */
    for (numofFC=0,thelistFC=STARTFC(myList); thelistFC!=NULL; thelistFC=SUCCFC(thelistFC))
    {
      numofFC++;
      if (thelistFC==LASTFC(myList))
        break;
    }

    totalnumofFC = NFC(myList) ;
    NFC(thenewFL) = totalnumofFC - numofFC;
    NFC(myList) = numofFC;

    /* determination of the orientation of the lists */
    if (FLORIENTATION(myList)==MATHPOS)
      FLORIENTATION(thenewFL) = FLORIENTATION(myList);
    else
    {
      /* myList is an 'island': one of the new lists is MATHNEG, the other MATHPOS */
      if (DetermineOrientation (myList))
        return (2);

      if (DetermineOrientation (thenewFL))
        return (2);
    }

    /* redistribute contained FLs on the two IFLs if necessary */
    if (NFL(theIFL)>1)
      if (RedistributeFLs (myList,thenewFL))
        return(2);


  }

  /* second case, the FCs are not in the same front list */
  else
  {
    /* merge the two lists */
    thesuccFC = SUCCFC(theFC);


    SUCCFC(theFC) = thenewdoubledFC;
    PREDFC(thenewdoubledFC) = theFC;

    PREDFC(thesuccFC) = thenewFC;
    SUCCFC(thenewFC) = thesuccFC;

    NFC(myList) += NFC(myfoundFL);

    if (LASTFC(myList)==theFC)
      STARTFC(myList) = SUCCFC(theFC);

    /*???*/	/* tell the former nodes of myfoundFL that they belong now to myList */
    for ( thelistFC =  thenewdoubledFC; thelistFC != NULL; thelistFC = SUCCFC(thelistFC) )
    {
      MYFL(thelistFC) = myList;

      if ( thelistFC == thenewFC )
        break;
    }

    /* indicate that myfoundFL is empty now and dispose it */
    STARTFC(myfoundFL) = NULL;
    LASTFC(myfoundFL) = NULL;

    *disp_FL = myfoundFL;
    /*no more here: DisposeFrontList (myfoundFL);*/
  }

#if (defined __MPW32__ || defined __MWCW__)
  if (plotfront)
  {


    UserWrite("front after updating <RETURN>\n");
    UserRead(buffer);
  }
#endif

  /* endof concerning link informations inf FC-datastructure */

  return (0);

}

/****************************************************************************/
/*                                                                          */
/* Function:  MakeElement                                                               */
/*                                                                          */
/* Purpose:   creating elements												*/
/*                                                                          */
/* Input:     nodes for the element									        */
/*																			*/
/*                                                                          */
/* Output:    return 0 if no error occur							                */
/*                                                                          */
/****************************************************************************/

static INT MakeElement (GRID *theGrid, ELEMENT_CONTEXT* theElementContext)
{
  INT i,n,found,NeighborSide[4];
  NODE *Node[3];
  ELEMENT *theElement,*Neighbor[4];
  INT reply;
  DOUBLE_VECTOR position;
  BVP_DESC theBVPDesc;
  DOUBLE diff;

  n = 3;                /* generate triangle */

  Node[0] = theElementContext->theNode[0];
  Node[1] = theElementContext->theNode[1];
  Node[2] = theElementContext->theNode[2];

  V2_CLEAR(position);
  for (i=0; i<3; i++)
  {
    V2_ADD1(CVECT(MYVERTEX(Node[i])),position);
  }
  V2_SCALE(0.3333333333333333,position);
  if (BVP_SetBVPDesc(MG_BVP(MYMG(theGrid)),&theBVPDesc)) return (1);
  V2_EUKLIDNORM_OF_DIFF(theBVPDesc.midpoint,position,diff);
  if (diff>theBVPDesc.radius)
  {
    UserWrite("\nERROR: trying to create element outside bounding sphere of domain\n");
    return (1);
  }


  /* find neighboring elements */
  found = 0;
  for (i=0; i<n; i++)
  {
    Neighbor[i]     = theElementContext->theNeighbour[i];
    NeighborSide[i] = theElementContext->Neighbourside[i];
  }

  /*
          theElement = InsertElement (theGrid,n,Node,Neighbor,NeighborSide);
   */
  IFDEBUG(gm,2)
  printf("ggmain.c: theMG=%x InsertElement n=%d Node=%x\n",MYMG(theGrid),n,Node);
  printf("ggmain.c: ndelemptrarray=%x\n",MGNDELEMPTRARRAY(MYMG(theGrid)));
  fflush(stdout);
  for (i=0; i<n; i++)
  {
    printf("ggmain.c: i=%d Node=%x id=%d\n",i,Node[i],ID(Node[i]));
    fflush(stdout);
  }
  for (i=0; i<NDELEM_BLKS_MAX; i++)
  {
    printf("i=%d ndblk=%x\n",i,MGNDELEMBLK(MYMG(theGrid),i));
    fflush(stdout);
  }
  ENDDEBUG


    theElement = InsertElement (theGrid,n,Node,NULL,NULL);
  if (theElement == NULL) RETURN(GM_FATAL);

  SETSUBDOMAIN(theElement,theElementContext->SubdomainID);
  theElementContext->thenewElement = theElement;

  /* alternativ O(N*N): InsertElement (MYMG(theGrid),n,Node,NULL,NULL);*/

  /* plot */

        #ifdef ModelP
  if (ID(theElement)==ElemID)
        #else
  if ((ID(theElement)==ElemID) || (UserInterrupt(NULL)))
        #endif
  {
    sprintf(buffer,"zoom current element (ID=%ld) (y/n)? ",ID(theElement));
    UserWrite(buffer);
    UserRead(buffer);


    UserWrite("break here (n/y)? ");
    UserRead(buffer);

    if (buffer[0]=='y')
      return (1);

    UserWrite("step (n/y)? "); UserRead(buffer);
    dostep = (buffer[0]=='y') ? YES : NO;

    UserWrite("update (n/y)? "); UserRead(buffer);
    reply = (buffer[0]=='y') ? YES : NO;
    doupdate = reply;
  }

  if (dostep)
  {
    sprintf(buffer,"current element (ID=%ld)\nstep/zoom & step/cont/break?",ID(theElement));
    UserWrite(buffer); UserRead(buffer);



    if (buffer[0]=='c')
      dostep = doupdate = NO;

    if (buffer[0]=='b')
      return (1);
  }

  return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  CalcNewPoint		                                                */
/*                                                                          */
/* Purpose:   calculate the coordinates of the new point					*/
/*                                                                          */
/* Input:                                                                                                                               */
/*                                                                          */
/* Output:    INT return code see header file                               */
/*                                                                          */
/****************************************************************************/

static INT CalcNewPoint (FRONTLIST *myFL,FRONTCOMP *theFC, DOUBLE xt[3], DOUBLE yt[3])
{
  VERTEX *theVertex;
  DOUBLE pos[2];
  DOUBLE norm,meshsize,hight;
  DOUBLE xP1minusP2,yP1minusP2;

  /* calculation of the coordinates of the inner node for triangulation */
  theVertex = MYVERTEX(FRONTN(theFC));
  xt[0] = XC(theVertex);
  yt[0] = YC(theVertex);
  theVertex = MYVERTEX(FRONTN(SUCCFC(theFC)));
  xt[1] = XC(theVertex);
  yt[1] = YC(theVertex);

  /* calculation of the norm of the perpend. vector of the triangle */
  xP1minusP2 = xt[0] -xt[1];
  yP1minusP2 = yt[0] -yt[1];
  norm = sqrt(xP1minusP2*xP1minusP2+yP1minusP2*yP1minusP2);

  pos[0] = 0.5*((xt[0]+xt[1]) + yP1minusP2*BARYCENTER*ROOT3);
  pos[1] = 0.5*((yt[0]+yt[1]) - xP1minusP2*BARYCENTER*ROOT3);
  (*nominal_h)(pos,&meshsize);
  if (equilateral)
  {
    hight = meshsize*meshsize-0.25*norm*norm;

    /* possibly: hight = MAX(hight,factor*norm*norm); */

    if (hight<0)
    {
      PrintErrorMessage('E',"CalcNewPoint","hight<0");
      return (-3);
    }
    hight = sqrt(hight);
  }
  else
  if(doConstDel)
    hight = meshsize/1.2;
  else
    hight = meshsize;

  xt[2] = 0.5*(xt[0]+xt[1]) + yP1minusP2*hight/norm;
  yt[2] = 0.5*(yt[0]+yt[1]) - xP1minusP2*hight/norm;

  /* definition of the minimal distance of the created node to the surrounding nodes */
  searchradius = myPars->searchconst*meshsize;

  /* comparing the squares avoids sqrt */
  searchradius2 = searchradius*searchradius;

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  CreateOrSelectFC	                                                */
/*                                                                          */
/* Purpose:   find the FC for creating the new triangle						*/
/*                                                                          */
/* Input:                                                                                                                               */
/*                                                                          */
/* Output:    INT return code see header file                               */
/*                                                                          */
/****************************************************************************/

static FRONTCOMP *CreateOrSelectFC (
  GRID *theGrid,
  INDEPFRONTLIST *theIFL,
  FRONTLIST *myList,
  FRONTCOMP *theFC,
  FRONTCOMP *theProposedFC,
  FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS],
  DOUBLE xt[3], DOUBLE yt[3],
  INT checkoptions,
  INT recursiondepth,
  ELEMENT_CONTEXT *theElementContext)
{
  FRONTCOMP *thenewFC;
  NODE *theNode;
  VERTEX *theVertex;
  DOUBLE_VECTOR pos;

  if (recursiondepth > 20)
  {
    PrintErrorMessage('E',"CreateOrSelectFC","recursiondepth > 10 in CreateOrSelectFC");
    sprintf(buffer,"constructing element over %ld %ld\n",ID(FRONTN(theFC)),ID(FRONTN(SUCCFC(theFC))));
    UserWrite(buffer);
    return (NULL);
  }

  if (theProposedFC!=NULL)
  {
    theVertex = MYVERTEX(FRONTN(theProposedFC));
    xt[2] = XC(theVertex);
    yt[2] = YC(theVertex);
  }

  if (checkoptions&CHECKNEAR)
    if ((thenewFC=NearerNeighbour (theFC,xt[2],yt[2]))!=NULL)
    {
      theElementContext->n = 1;
      return (CreateOrSelectFC (theGrid,theIFL,myList,theFC,thenewFC,theIntersectfoundPoints,xt,yt,CHECKINSIDE,++recursiondepth, theElementContext));
    }

  if (checkoptions&CHECKNBCUT)
    if ((thenewFC=CutNeighbour (theFC,xt,yt))!=NULL)
      if (thenewFC!=theProposedFC)
      {
        theElementContext->n = 1;
        return (CreateOrSelectFC (theGrid,theIFL,myList,theFC,thenewFC,theIntersectfoundPoints,xt,yt,CHECKINSIDE,++recursiondepth, theElementContext));
      }


  if (checkoptions&(CHECKNEAR | CHECKINSIDE))
  {
    if (FCinsideOrNear (theIFL,theFC,theProposedFC, theIntersectfoundPoints, xt,yt,&thenewFC,checkoptions))
    {
      theElementContext->n = 4;
      return (NULL);                                            /* error occured */
    }
    else
    if (thenewFC != NULL)                               /* point found   */
    {
      return (CreateOrSelectFC (theGrid,theIFL,myList,theFC,thenewFC,theIntersectfoundPoints,xt,yt,CHECKINSIDE | CHECKINTERSECT,++recursiondepth, theElementContext));
    }
  }

  if (checkoptions&CHECKINTERSECT)
    if ((thenewFC=FrontLineIntersection (theIFL,theFC,theProposedFC, theIntersectfoundPoints, xt,yt))!=NULL)
      if (thenewFC!=theProposedFC)
      {
        theElementContext->n = 4;
        return (CreateOrSelectFC (theGrid,theIFL,myList,theFC,thenewFC,theIntersectfoundPoints, xt,yt,CHECKNBCUT | CHECKINSIDE,++recursiondepth, theElementContext));
      }

  if (theProposedFC!=NULL)
    return (theProposedFC);

  /* create a new FC including node and vertex */
  pos[0] = xt[2];
  pos[1] = yt[2];
  theNode = InsertInnerNode (theGrid,pos);
  if (theNode == NULL)
    return(NULL);

  thenewFC = CreateFrontComp (myList, theFC, 1, &theNode);
  if (thenewFC==NULL)
  {
    PrintErrorMessage('E',"CreateOrSelectFC","no storage for new FC");
    return (NULL);

  }

  return (thenewFC);
}

static INT CheckNewElement(FRONTLIST *theFL,  DOUBLE xt[3], DOUBLE yt[3])
{
  FRONTCOMP *thecompFC;
  VERTEX *theVertex;
  DOUBLE xM, yM, xR, yR, xQ, yQ, lambda1, lambda2;
  DOUBLE denominator;
  DOUBLE lambdacomp = MAX_C;
  INT i;

  xM = 0.5*(xt[0] + xt[1]);
  yM = 0.5*(yt[0] + yt[1]);

  thecompFC=STARTFC(theFL);
  for(i=0; i<theFL->nFrontcomp; i++)
  {
    /* calculation of the coordinates of the FC's to proove */
    theVertex = MYVERTEX(FRONTN(SUCCFC(thecompFC)));
    xR = XC(theVertex);
    yR = YC(theVertex);
    theVertex = MYVERTEX(FRONTN(thecompFC));
    xQ = XC(theVertex);
    yQ = YC(theVertex);

    denominator = (yR-yQ)*(xt[2]-xM)-(xR-xQ)*(yt[2]-yM);
    if (fabs(denominator)<SMALLDOUBLE)
      if (thecompFC==LASTFC(theFL))
        break;
      else
        continue;                               /* lines are parallel */

    if (fabs(xt[2]-xM)<SMALLDOUBLE)
    {
      lambda2 = (xR-xM)/(xR-xQ);
      lambda1 = (yR-yM - lambda2*(yR-yQ)) / (yt[2]-yM);
    }
    else
    {
      lambda2=((yR-yM)*(xt[2]-xM)-(xR-xM)*(yt[2]-yM))/denominator;
      lambda1=(xR-xM-lambda2*(xR-xQ))/(xt[2]-xM);
    }
    if (lambda1 < 0.9999999999 && lambda1 > 0.000000001 && lambda2 < 0.9999999999 && lambda2 > 0.000000001)
      return(0);
    thecompFC=SUCCFC(thecompFC);
  }
  return(1);
}

static INT Cross_Point(DOUBLE p0[2], DOUBLE n0[2], DOUBLE p1[2], DOUBLE n1[2], DOUBLE *cp)
{
  DOUBLE lambda0, lambda1, a, cp0[2], cp1[2];

  a = - n0[0] * n1[1] + n0[1] * n1[0];

  lambda0 = ( n1[1] * (p0[0] - p1[0]) - n1[0] * (p0[1] - p1[1]) ) / a;
  lambda1 = ( n0[1] * (p0[0] - p1[0]) - n0[0] * (p0[1] - p1[1]) ) / a;

  cp0[0] = p0[0] + lambda0 * n0[0];
  cp0[1] = p0[1] + lambda0 * n0[1];

  cp1[0] = p1[0] + lambda1 * n1[0];
  cp1[1] = p1[1] + lambda1 * n1[1];

  /*	if(sqrt( (cp0[0]-cp1[0])*(cp0[0]-cp1[0]) + (cp0[1]-cp1[1])*(cp0[1]-cp1[1]) )>SMALLDOUBLE)
                  printf("%s\n", "ERROR");*/

  cp[0] = cp0[0];
  cp[1] = cp0[1];
  return(0);
}

static DOUBLE Calc_Circumcircusmidpoint(DOUBLE xt[3], DOUBLE yt[3], DOUBLE *midp)
{
  DOUBLE p0[2], p1[2], p2[2], n0[2], n1[2], n2[2], cp01[2], l;

  /* edgemidpoints */
  p0[0] = 0.5 * ( xt[0] + xt[1] );
  p0[1] = 0.5 * ( yt[0] + yt[1] );

  p1[0] = 0.5 * ( xt[1] + xt[2] );
  p1[1] = 0.5 * ( yt[1] + yt[2] );

  p2[0] = 0.5 * ( xt[2] + xt[0] );
  p2[1] = 0.5 * ( yt[2] + yt[0] );

  /* normals */
  n0[0] = - (yt[1] - yt[0]);
  n0[1] =   (xt[1] - xt[0]);
  l = sqrt(n0[0]*n0[0]+n0[1]*n0[1]);
  n0[0] = n0[0] / l;
  n0[1] = n0[1] / l;

  n1[0] = - (yt[2] - yt[1]);
  n1[1] =   (xt[2] - xt[1]);
  l = sqrt(n1[0]*n1[0]+n1[1]*n1[1]);
  n1[0] = n1[0] / l;
  n1[1] = n1[1] / l;

  n2[0] = - (yt[0] - yt[2]);
  n2[1] =   (xt[0] - xt[2]);
  l = sqrt(n2[0]*n2[0]+n2[1]*n2[1]);
  n2[0] = n2[0] / l;
  n2[1] = n2[1] / l;

  Cross_Point(p0, n0, p1, n1, cp01);
  /*	Cross_Point(p1, n1, p2, n2, cp12);
          Cross_Point(p2, n2, p0, n0, cp20);

          if( (sqrt( (cp01[0]-cp12[0])*(cp01[0]-cp12[0]) + (cp01[1]-cp12[1])*(cp01[1]-cp12[1]) )>SMALLDOUBLE) ||
              (sqrt( (cp12[0]-cp20[0])*(cp12[0]-cp20[0]) + (cp12[1]-cp20[1])*(cp12[1]-cp20[1]) )>SMALLDOUBLE) ||
              (sqrt( (cp20[0]-cp01[0])*(cp20[0]-cp01[0]) + (cp20[1]-cp01[1])*(cp20[1]-cp01[1]) )>SMALLDOUBLE))
          {
                  printf("%s\n", "ERROR");
                  midp[0] = MAXDOUBLE;
                  midp[1] = MAXDOUBLE;
          }
          else
          {
                  midp[0] = cp01[0];
                  midp[1] = cp01[1];
          }*/

  midp[0] = cp01[0];
  midp[1] = cp01[1];
  return(0.0);
}

static FRONTCOMP *CreateDelaunayTriangle (
  GRID *theGrid,
  MESH2D mesh_2d,
  INDEPFRONTLIST *theIFL,
  FRONTLIST *theFL,
  FRONTCOMP *theFC,
  DOUBLE xt[3], DOUBLE yt[3])
{
  FRONTCOMP *thenewFC, *thereturnFC;
  NODE *theNode;
  VERTEX *theVertex;
  DOUBLE meshsize, pos[2], circ, min_circ;
  INT i, subdomain, ok;
  DOUBLE midp[2], np_midp[2], np[2];
  DOUBLE np_circ;
  ELEMENT_2D *help;

  min_circ = MAXDOUBLE;
  thereturnFC = NULL;
  subdomain = theFL->SubdomainID;

  pos[0] = 0.5*(xt[0]+xt[1]);
  pos[1] = 0.5*(yt[0]+yt[1]);
  (*nominal_h)(pos,&meshsize);

  /* first: test the new point */
  if(CheckNewElement(theFL, xt, yt))                            /* no crossing with the Front */
  {
    Calc_Circumcircusmidpoint(xt, yt, np_midp);
    np_circ = sqrt( (np_midp[0]-xt[0])*(np_midp[0]-xt[0]) + (np_midp[1]-yt[0])*(np_midp[1]-yt[0]) );
    np[0] = xt[2];
    np[1] = yt[2];
  }
  else
    np_circ = MAXDOUBLE;

  thenewFC=SUCCFC(SUCCFC(theFC));
  for(i=0; i<theFL->nFrontcomp-2; i++)
  {
    theVertex = MYVERTEX(FRONTN(thenewFC));
    xt[2] = XC(theVertex);
    yt[2] = YC(theVertex);

    if( ( ( (xt[1]-xt[0])*(yt[2]-yt[0])-(yt[1]-yt[0])*(xt[2]-xt[0]) )>SMALLDOUBLE ) )
      if(CheckNewElement(theFL, xt, yt))                                        /* no crossing with the Front */
      {
        Calc_Circumcircusmidpoint(xt, yt, midp);
        /*				circ0 = sqrt( (midp[0]-xt[0])*(midp[0]-xt[0]) + (midp[1]-yt[0])*(midp[1]-yt[0]) );
                                        circ1 = sqrt( (midp[0]-xt[1])*(midp[0]-xt[1]) + (midp[1]-yt[1])*(midp[1]-yt[1]) );
                                        circ2 = sqrt( (midp[0]-xt[2])*(midp[0]-xt[2]) + (midp[1]-yt[2])*(midp[1]-yt[2]) );
                                        if( (sqrt((circ0-circ1)*(circ0-circ1))<SMALLDOUBLE) &&
                                            (sqrt((circ1-circ2)*(circ1-circ2))<SMALLDOUBLE) &&
                                                (sqrt((circ2-circ0)*(circ2-circ0))<SMALLDOUBLE) )
                                                circ = circ1;
                                        else
                                                printf("%s\n", "ERROR");*/

        circ = sqrt( (midp[0]-xt[0])*(midp[0]-xt[0]) + (midp[1]-yt[0])*(midp[1]-yt[0]) );
        if(min_circ>circ)
        {
          min_circ = circ;
          thereturnFC = thenewFC;
        }
      }
    thenewFC = SUCCFC(thenewFC);
  }

  if(min_circ>1.5*np_circ)                              /* take the new point */
  {
    /* first check, if new point lies in the circumcircle of a created element */
    ok = 1;
    help = mesh_2d.start[subdomain];
    for(i=0; i<mesh_2d.nElements[subdomain]; i++)
    {
      if(sqrt( (help->del_context.midp[0]-np[0])*(help->del_context.midp[0]-np[0])
               +(help->del_context.midp[1]-np[1])*(help->del_context.midp[1]-np[1]) )
         < help->del_context.radius)
        ok = 0;
      help = help->succel;
    }
    ok = 1;
    if(ok)                                                              /* everything ok */
    {
      pos[0] = np[0];
      pos[1] = np[1];
      xt[2] = np[0];
      yt[2] = np[1];
      theNode = InsertInnerNode (theGrid,pos);
      if (theNode == NULL)
        return(NULL);

      thereturnFC = CreateFrontComp (theFL, theFC, 1, &theNode);
      if (thereturnFC==NULL)
      {
        PrintErrorMessage('E',"CreateOrSelectFC","no storage for new FC");
        return (NULL);
      }
    }
    else
    if(min_circ>meshsize)
    {
      /*	help = mesh_2d.start[subdomain];
              for(i=0;i<mesh_2d.nElements[subdomain];i++)
              {
                      if(sqrt( (help->del_context.midp[0]-np[0])*(help->del_context.midp[0]-np[0])
         +(help->del_context.midp[1]-np[1])*(help->del_context.midp[1]-np[1]) )
                                      < help->del_context.radius)
                      {
                              if(help->succel!=NULL)
                                      help->succel = help->succel->succel;
                              mesh_2d.nElements[subdomain]--;
                      }
                      help = help->succel;
              }*/

    }
  }

  return (thereturnFC);
}


/****************************************************************************/
/*                                                                          */
/* Function:  FillElementContext                                                */
/*                                                                          */
/* Purpose:   create the ElementContext for MakeElement                             */
/*                                                                          */
/* Input:     ELEMENT_CONTEXT* theElementContext, FRONTCOMP* theFC          */
/*                                                                          */
/* Output:    ELEMENT_CONTEXT* theElementContext                            */
/*                                                                          */
/****************************************************************************/

static INT FillElementContext(INT FlgForAccel, ELEMENT_CONTEXT* theElementContext, FRONTCOMP* theFC, FRONTCOMP* thenewFC, FRONTCOMP* the_old_succ)
{
  theElementContext->theNode[0] = FRONTN(theFC);
  theElementContext->theNode[1] = FRONTN(the_old_succ);
  theElementContext->theNode[2] = FRONTN(thenewFC);

  theElementContext->theVertex[0] = (MYVERTEX(FRONTN(theFC)));
  theElementContext->theVertex[1] = (MYVERTEX(FRONTN(the_old_succ)));
  theElementContext->theVertex[2] = (MYVERTEX(FRONTN(thenewFC)));

  theElementContext->SubdomainID = theFC->myFL->SubdomainID;

  switch (FlgForAccel)
  {
  case NORMALCASE :
    theElementContext->theNeighbour[0] = FCNGB(theFC);
    theElementContext->theNeighbour[1] = NULL;
    theElementContext->theNeighbour[2] = NULL;
    theElementContext->Neighbourside[0] = FCNGBS(theFC);
    theElementContext->Neighbourside[1] = -1;
    theElementContext->Neighbourside[2] = -1;
    break;

  case LEFTNEIGHBOURCASE :
    theElementContext->theNeighbour[0] = FCNGB(theFC);
    theElementContext->theNeighbour[1] = NULL;
    theElementContext->theNeighbour[2] = FCNGB(PREDFC(theFC));
    theElementContext->Neighbourside[0] = FCNGBS(theFC);
    theElementContext->Neighbourside[1] = -1;
    theElementContext->Neighbourside[2] = FCNGBS(PREDFC(theFC));
    break;

  case RIGHTNEIGHBOURCASE :
    theElementContext->theNeighbour[0] = FCNGB(theFC);
    theElementContext->theNeighbour[1] = FCNGB(the_old_succ);
    theElementContext->theNeighbour[2] = NULL;
    theElementContext->Neighbourside[0] = FCNGBS(theFC);
    theElementContext->Neighbourside[1] = FCNGBS(the_old_succ);
    theElementContext->Neighbourside[2] = -1;
    break;

  case ININTERCASE :
    theElementContext->theNeighbour[0] = FCNGB(theFC);
    theElementContext->theNeighbour[1] = NULL;
    theElementContext->theNeighbour[2] = NULL;
    theElementContext->Neighbourside[0] = FCNGBS(theFC);
    theElementContext->Neighbourside[1] = -1;
    theElementContext->Neighbourside[2] = -1;
    break;

  case FINALCASE :
    theElementContext->theNeighbour[0] = FCNGB(theFC);
    theElementContext->theNeighbour[1] = FCNGB(the_old_succ);
    theElementContext->theNeighbour[2] = FCNGB(PREDFC(theFC));
    theElementContext->Neighbourside[0] = FCNGBS(theFC);
    theElementContext->Neighbourside[1] = FCNGBS(the_old_succ);
    theElementContext->Neighbourside[2] = FCNGBS(PREDFC(theFC));
    break;

  default :
    return (1);

  }

  if(theElementContext == NULL) return(1);
  return(0);
}

static INT FillDelContext(ELEMENT_CONTEXT* theElementContext, DELAUNAY_CONTEXT* theDelaunayContext)
{
  DOUBLE xt[3], yt[3], midp[2], radius;
  INT i;

  for(i=0; i<3; i++)
  {
    xt[i] = XC(MYVERTEX(theElementContext->theNode[i]));
    yt[i] = YC(MYVERTEX(theElementContext->theNode[i]));
  }

  Calc_Circumcircusmidpoint(xt, yt, midp);

  radius = sqrt( (midp[0]-xt[0])*(midp[0]-xt[0]) + (midp[1]-yt[0])*(midp[1]-yt[0]) );

  theDelaunayContext->midp[0] = midp[0];
  theDelaunayContext->midp[1] = midp[1];
  theDelaunayContext->radius = radius;

  if(theDelaunayContext == NULL)
    return(1);
  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Function:  FrontcomponentUpdate                                              */
/*                                                                          */
/* Purpose:   makes an Update of the frontcomponent informations "Neighbour"*/
/*                        , "Side" and "Neighbourside"									*/
/*                                                                          */
/* Input:    FRONTCOMP* theFC,the_old_succ,thenewFC and the flag describing */
/*                       the concerning case !											*/
/*                                                                          */
/* Output:	  returns 1 if error occurs                                     */
/*																			*/
/****************************************************************************/

static INT FrontcomponentUpdate(INT FlgForAccel, FRONTCOMP* theFC, FRONTCOMP* the_old_succ, FRONTCOMP* thenewFC, ELEMENT_CONTEXT* theElementContext )
{
  /* 1.) theFC : */
  FCNGB(theFC) = theElementContext->thenewElement;
  FCNGBS(theFC) = 2;

  /* 3.) thenewFC : */
  if(FlgForAccel != RIGHTNEIGHBOURCASE)       /*special case*/
  {
    FCNGB(thenewFC) = theElementContext->thenewElement;
    FCNGBS(thenewFC) = 1;
  }

  if(theFC == NULL) return(1);
  return(0);

}

/****************************************************************************/
/*                                                                          */
/* Function:  FL_FC_Disposer                                                            */
/*                                                                          */
/* Purpose:   disposing of frontcomponents and Frontlists                                       */
/*                                                                          */
/* Input:    FRONTCOMP *disp_FC, FRONTLIST *disp_FL							*/
/*                                                                          */
/* Output:	  returns 1 if error occurs                                     */
/*																			*/
/****************************************************************************/

static INT FL_FC_Disposer(FRONTCOMP *disp_FC, FRONTLIST *disp_FL)
{
  if (disp_FC != NULL)
    DisposeFrontComp (disp_FL, disp_FC);
  else if (disp_FL != NULL)
    DisposeFrontList(disp_FL);
  return(0);
}

static INT PrintFront(MULTIGRID *theMG)
{
  INDEPFRONTLIST *theIFL;
  FRONTLIST *theFL;
  FRONTCOMP *theFC;
  VERTEX *theVertex;
  INT i, j, l;

  theIFL=LASTIFL(myMGdata);
  for(i=0; i<myMGdata->nIndepFrontlist; i++)
  {
    theFL=LASTFL(theIFL);
    for(j=0; j<theIFL->nFrontlist; j++)
    {
      theFC=STARTFC(theFL);
      for(l=0; l<theFL->nFrontcomp; l++)
      {
        theVertex = MYVERTEX(FRONTN(theFC));
        printf("%lf %lf\n", XC(theVertex), YC(theVertex));
        theFC = SUCCFC(theFC);
      }
      printf("\n");
      theFL=PREDFL(theFL);
    }
    printf("####\n");
    theIFL=PREDIFL(theIFL);
  }
  printf("#################################################\n");

  return (0);
}

/****************************************************************************/
/*D
   GenerateGrid - Create the advancing frontlists and generates the grid in 2d

   SYNOPSIS:
   INT GenerateGrid (MULTIGRID *theMG, GG_ARG *MyArgs, GG_PARAM *param);

    PARAMETERS:
   .   theMG - pointer to the multigrid
   .   MyArgs - structure for reading grid generator parameters
   .   param -  structure for reading the grid generating control parameters

    DESCRIPTION:
        This function creates automatical the triangular grid in 2d.

    RETURN VALUE:
    INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

static INT debug;

INT GenerateGrid (MULTIGRID *theMG, GG_ARG *MyArgs, GG_PARAM *param, MESH *mesh, CoeffProcPtr coeff, INT Single_Mode, INT display)
{
  GRID *theGrid;
  INDEPFRONTLIST *theIFL,*nextIFL;
  FRONTLIST *myList;
  FRONTCOMP *theFC,*thesuccFC,*thenewFC;
  FRONTCOMP *the_old_succ;
  DOUBLE xt[3],yt[3];
  INT printelem,FlgForAccel,nElement, i, j;
  FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS];
  ELEMENT_CONTEXT theElementContext;
  ELEMENT_2D *elem_context_list;
  MESH2D mesh_2d;
  ELEMENT_2D *elem_2d;
  FRONTCOMP *disp_FC;
  FRONTLIST *disp_FL;
  INT MarkKey=MG_MARK_KEY(theMG);


  SingleMode = Single_Mode;
  nElement = 0;

  SetFlagsfortemporaryGGObjects(IflObj, FlObj, FcObj);

  myMGdata = GetMGdataPointer(theMG);
  myPars = param;

  theGrid = GRID_ON_LEVEL(theMG,0);

  if (coeff == NULL)
  {
    nominal_h = GlobalMeshsize;
    H_global = param->h_global;
  }
  else
    nominal_h = coeff;

  /* check options */
  ElemID = -1;

  PRINTDEBUG(dom,1,("  h %f \n",H_global));

  doanimate       = MyArgs->doanimate;
  doupdate        = MyArgs->doupdate;
  dostep          = MyArgs->dostep;
  plotfront       = MyArgs->plotfront;
  printelem       = MyArgs->printelem;
  equilateral     = MyArgs->equilateral;
  doedge          = MyArgs->doedge;
  doangle         = MyArgs->doangle;
  doEdge          = MyArgs->doEdge;
  doAngle         = MyArgs->doAngle;
  doConstDel  = MyArgs->doConstDel;

  if(!(doConstDel))
    if(((doedge == YES) && ((doangle || doEdge || doAngle) == YES )) ||
       ((doangle == YES) && ((doedge || doEdge || doAngle) == YES )) ||
       ((doEdge == YES) && ((doangle || doedge || doAngle) == YES )) ||
       ((doAngle == YES) && ((doangle || doEdge || doedge) == YES )) ||
       ((doAngle || doangle || doEdge || doedge) == NO ))
    {
      PrintErrorMessage('E',"GenerateGrid","no variable chosen for accelerate or not!");
      return (2);
    }

  /* are there already indep front lists? */
  if (STARTIFL(myMGdata)!=NULL)
  {
    PrintErrorMessage('E',"GenerateGrid","there exist already independent front lists");
    return (4);
  }

  if (AssembleFrontLists (theMG,mesh,MarkKey)!=0)
  {
    return (5);
  }

  /* now we can create all element sides that will be needed by the MakeElement fct. */

  if (doAngle || doEdge)
    if (AccelInit(theGrid, doAngle, doEdge, myPars)!=0) return(1);

  if (doConstDel)
  {
    /* allocate memory for mesh_2d */
    elem_context_list = NULL;
    mesh_2d.nElements = (INT *) GetTmpMem(theMG->theHeap,(theMG->theBVPD.nSubDomains+1)*sizeof(INT),MarkKey);
    if (mesh_2d.nElements == NULL)
      return(NULL);
    mesh_2d.Element = (ELEMENT_2D**) GetTmpMem(theMG->theHeap,(theMG->theBVPD.nSubDomains+1)*sizeof(ELEMENT_2D*),MarkKey);
    if (mesh_2d.Element == NULL)
      return(NULL);
    mesh_2d.start = (ELEMENT_2D**) GetTmpMem(theMG->theHeap,(theMG->theBVPD.nSubDomains+1)*sizeof(ELEMENT_2D*),MarkKey);
    if (mesh_2d.start == NULL)
      return(NULL);
    for(i=1; i<=theMG->theBVPD.nSubDomains; i++)
    {
      mesh_2d.Element[i] = NULL;
      mesh_2d.start[i] = NULL;
      mesh_2d.nElements[i] = 0;
    }
  }



  /*************************************************************************************/
  /*                                                                                     */
  /* creating inner nodes	and vertices for automatically triangulation                             */
  /* loops for the indep. front lists and front lists begin at the end of the lists      */
  /* according to the possibility of creating of new indep. front lists or front lists */
  /*                                                                                     */
  /*************************************************************************************/

  for (theIFL=LASTIFL(myMGdata); theIFL!=NULL; theIFL=nextIFL)
  {
    if (doAngle || doEdge)
      while ((theFC=AccelBaseTreeSearch(&myList)) != NULL)
      {
        theIFL = myList->myIFL;
        the_old_succ = SUCCFC(theFC);
        /* are there only 3 FCs left and lie they not on an inner hole? */
        if (PREDFC(theFC)==SUCCFC(SUCCFC(theFC)) && NFL(theIFL) == 1)
        {
          /* we make this last element and dispose the list */
          /* accelerator final case */
          FlgForAccel = FINALCASE;
          AccelUpdate( theFC, PREDFC(theFC), the_old_succ, FlgForAccel,doAngle,doEdge);

          if (FillElementContext(FlgForAccel, &theElementContext, theFC, PREDFC(theFC), the_old_succ))
            return (1);

          if (MakeElement(theGrid, &theElementContext ))
            return (8);
          if (display>0)
          {
            nElement++;
            if (nElement%display==0)
            {
              if (nElement%(10*display)==0) UserWrite("\n");
              UserWriteF("[%d] ",(int)nElement);
            }
          }
          DisposeFrontList(myList);

          /* in this case "FrontcomponentUpdate(...);" is redundant*/

          continue;
        }

        CalcNewPoint (myList,theFC,xt,yt);

        thesuccFC = SUCCFC(theFC);

        theIntersectfoundPoints[0] = NULL;

        if ((thenewFC=CreateOrSelectFC(theGrid,theIFL,myList,theFC,NULL, theIntersectfoundPoints, xt,yt,CHECKALL,0, &theElementContext))==NULL)
          return (9);

        FlgForAccel = NORMALCASE;

        disp_FC = NULL;
        disp_FL = NULL;

        if (FrontLineUpDate (theGrid,theIFL,myList,theFC,thenewFC,&FlgForAccel, the_old_succ, &disp_FC, &disp_FL))
          return (10);

        if (FillElementContext(FlgForAccel, &theElementContext, theFC, thenewFC, the_old_succ ))
          return (1);

        AccelUpdate( theFC, thenewFC, the_old_succ, FlgForAccel, doAngle, doEdge);

        if (MakeElement(theGrid,&theElementContext))
          return (11);
        if (display>0)
        {
          nElement++;
          if (nElement%display==0)
          {
            if (nElement%(10*display)==0) UserWrite("\n");
            UserWriteF("[%d] ",(int)nElement);
          }
        }

        if (FrontcomponentUpdate(FlgForAccel, theFC, the_old_succ, thenewFC, &theElementContext))
          return (1);

        if(FL_FC_Disposer(disp_FC, disp_FL))
          return (1);



        if (printelem)
        {
          sprintf(buffer,"ELEMID %ld done\n",ID(LASTELEMENT(theGrid)));
          UserWrite(buffer);
        }
      }                   /* while */
    else
    if (doedge)
      while ((theFC=ChooseFCminside(theIFL,&myList)) != NULL)
      {
        theIFL = myList->myIFL;
        the_old_succ = SUCCFC(theFC);
        /* are there only 3 FCs left and lie they not on an inner hole? */
        if (PREDFC(theFC)==SUCCFC(SUCCFC(theFC)) && NFL(theIFL) == 1)
        {
          FlgForAccel = FINALCASE;
          /* we make this last element and dispose the list */
          if (FillElementContext(FlgForAccel, &theElementContext, theFC, PREDFC(theFC), the_old_succ))
            return (1);
          if (MakeElement(theGrid, &theElementContext))
            return (12);
          if (display>0)
          {
            nElement++;
            if (nElement%display==0)
            {
              if (nElement%(10*display)==0) UserWrite("\n");
              UserWriteF("[%d] ",(int)nElement);
            }
          }

          /* in this case "FrontcomponentUpdate(...);" is redundant*/

          DisposeFrontList(myList);
          continue;
        }

        CalcNewPoint (myList,theFC,xt,yt);

        thesuccFC = SUCCFC(theFC);

        theIntersectfoundPoints[0] = NULL;

        if ((thenewFC=CreateOrSelectFC(theGrid,theIFL,myList,theFC,NULL, theIntersectfoundPoints, xt,yt,CHECKALL,0, &theElementContext))==NULL)
          return (13);

        FlgForAccel = NORMALCASE;

        disp_FC = NULL;
        disp_FL = NULL;

        if (FrontLineUpDate (theGrid,theIFL,myList,theFC,thenewFC,&FlgForAccel, the_old_succ, &disp_FC, &disp_FL))
          return (10);

        if (FillElementContext(FlgForAccel, &theElementContext, theFC, thenewFC, the_old_succ))
          return (1);

        if (MakeElement(theGrid, &theElementContext))
          return (15);

        if (display>0)
        {
          nElement++;
          if (nElement%display==0)
          {
            if (nElement%(10*display)==0) UserWrite("\n");
            UserWriteF("[%d] ",(int)nElement);
          }
        }

        if (FrontcomponentUpdate(FlgForAccel, theFC, the_old_succ, thenewFC, &theElementContext))
          return (1);

        debug++;

        if(FL_FC_Disposer(disp_FC, disp_FL))
          return (1);

        if (printelem)
        {
          sprintf(buffer,"ELEMID %ld done\n",ID(LASTELEMENT(theGrid)));
          UserWrite(buffer);
        }
      }                           /* while */
    else
    if (doConstDel)
    {
      while ((theFC=ChooseFCminside(theIFL,&myList)) != NULL)
      {
        theIFL = myList->myIFL;
        the_old_succ = SUCCFC(theFC);
        /* are there only 3 FCs left and lie they not on an inner hole? */
        if (PREDFC(theFC)==SUCCFC(SUCCFC(theFC)) && NFL(theIFL) == 1)
        {
          FlgForAccel = FINALCASE;
          /* we make this last element and dispose the list */

          elem_2d = (ELEMENT_2D *) GetTmpMem(theMG->theHeap,sizeof(ELEMENT_2D),MarkKey);
          if (elem_2d == NULL)
            return(NULL);
          elem_2d->succel = NULL;

          if (FillElementContext(FlgForAccel, &(elem_2d->elem_context), theFC, PREDFC(theFC), the_old_succ))
            return (1);

          if (FillDelContext(&(elem_2d->elem_context), &(elem_2d->del_context)))
            return (1);

          if(mesh_2d.nElements[MYFL(theFC)->SubdomainID]==0)
          {
            mesh_2d.start[MYFL(theFC)->SubdomainID] = elem_2d;
            mesh_2d.Element[MYFL(theFC)->SubdomainID] = elem_2d;
          }
          else
          {
            mesh_2d.Element[MYFL(theFC)->SubdomainID]->succel = elem_2d;
            mesh_2d.Element[MYFL(theFC)->SubdomainID] = mesh_2d.Element[MYFL(theFC)->SubdomainID]->succel;
          }
          mesh_2d.nElements[MYFL(theFC)->SubdomainID]++;

          if (display>0)
          {
            nElement++;
            if (nElement%display==0)
            {
              if (nElement%(10*display)==0) UserWrite("\n");
              UserWriteF("[%d] ",(int)nElement);
            }
          }

          /* in this case "FrontcomponentUpdate(...);" is redundant*/

          DisposeFrontList(myList);
          continue;
        }

        CalcNewPoint (myList,theFC,xt,yt);

        /*	PrintFront(theMG);*/
        thesuccFC = SUCCFC(theFC);

        theIntersectfoundPoints[0] = NULL;

        if ((thenewFC=CreateDelaunayTriangle(theGrid,mesh_2d,theIFL,myList,theFC,xt,yt))==NULL)
          return (13);

        FlgForAccel = NORMALCASE;

        disp_FC = NULL;
        disp_FL = NULL;

        if (FrontLineUpDate (theGrid,theIFL,myList,theFC,thenewFC,&FlgForAccel, the_old_succ, &disp_FC, &disp_FL))
          return (10);

        elem_2d = (ELEMENT_2D *) GetTmpMem(theMG->theHeap,sizeof(ELEMENT_2D),MarkKey);
        if (elem_2d == NULL)
          return(NULL);
        elem_2d->succel = NULL;

        if (FillElementContext(FlgForAccel, &(elem_2d->elem_context), theFC, thenewFC, the_old_succ))
          return (1);

        if (FillDelContext(&(elem_2d->elem_context), &(elem_2d->del_context)))
          return (1);

        if (display>0)
        {
          nElement++;
          if (nElement%display==0)
          {
            if (nElement%(10*display)==0) UserWrite("\n");
            UserWriteF("[%d] ",(int)nElement);
          }
        }

        if (FrontcomponentUpdate(FlgForAccel, theFC, the_old_succ, thenewFC, &(elem_2d->elem_context)))
          return (1);

        if(mesh_2d.nElements[MYFL(theFC)->SubdomainID]==0)
        {
          mesh_2d.start[MYFL(theFC)->SubdomainID] = elem_2d;
          mesh_2d.Element[MYFL(theFC)->SubdomainID] = elem_2d;
        }
        else
        {
          mesh_2d.Element[MYFL(theFC)->SubdomainID]->succel = elem_2d;
          mesh_2d.Element[MYFL(theFC)->SubdomainID] = mesh_2d.Element[MYFL(theFC)->SubdomainID]->succel;
        }
        mesh_2d.nElements[MYFL(theFC)->SubdomainID]++;

        debug++;

        if(FL_FC_Disposer(disp_FC, disp_FL))
          return (1);

        if (printelem)
        {
          sprintf(buffer,"ELEMID %ld done\n",ID(LASTELEMENT(theGrid)));
          UserWrite(buffer);
        }
      }
    }
    else
      while ((theFC=ChooseFCminangle(theIFL,&myList)) != NULL)
      {
        theIFL = myList->myIFL;
        the_old_succ = SUCCFC(theFC);
        /* are there only 3 FCs left and lie they not on an inner hole? */
        if (PREDFC(theFC)==SUCCFC(SUCCFC(theFC)) && NFL(theIFL) == 1)
        {
          FlgForAccel = FINALCASE;

          if (FillElementContext(FlgForAccel, &theElementContext, theFC, PREDFC(theFC), the_old_succ))
            return (1);

          /* we make this last element and dispose the list */
          if (MakeElement(theGrid, &theElementContext ))
            return (16);
          if (display>0)
          {
            nElement++;
            if (nElement%display==0)
            {
              if (nElement%(10*display)==0) UserWrite("\n");
              UserWriteF("[%d] ",(int)nElement);
            }
          }
          DisposeFrontList(myList);
          /* in this case "FrontcomponentUpdate(...);" is redundant*/
          continue;
        }

        CalcNewPoint (myList,theFC,xt,yt);

        thesuccFC = SUCCFC(theFC);

        theIntersectfoundPoints[0] = NULL;

        if ((thenewFC=CreateOrSelectFC(theGrid,theIFL,myList,theFC,NULL, theIntersectfoundPoints, xt,yt,CHECKALL,0, &theElementContext))==NULL)
          return (17);

        FlgForAccel = NORMALCASE;

        disp_FC = NULL;
        disp_FL = NULL;

        if (FrontLineUpDate (theGrid,theIFL,myList,theFC,thenewFC,&FlgForAccel, the_old_succ, &disp_FC, &disp_FL))
          return (10);


        if (FillElementContext(FlgForAccel, &theElementContext, theFC, thenewFC, the_old_succ))
          return (1);

        if (MakeElement(theGrid, &theElementContext))
          return (19);

        if (display>0)
        {
          nElement++;
          if (nElement%display==0)
          {
            if (nElement%(10*display)==0) UserWrite("\n");
            UserWriteF("[%d] ",(int)nElement);
          }
        }
        if (FrontcomponentUpdate(FlgForAccel, theFC, the_old_succ, thenewFC, &theElementContext))
          return (1);

        if(FL_FC_Disposer(disp_FC, disp_FL))
          return (1);

        if (printelem)
        {
          sprintf(buffer,"ELEMID %ld done\n",ID(LASTELEMENT(theGrid)));
          UserWrite(buffer);
        }
      }                                   /* while */


    nextIFL = PREDIFL(theIFL);                          /* remember pred before disposing */
    DisposeIndepFrontList(theIFL);

  }

  if(doConstDel)
  {
    /* transfer elements to theGrid */
    for(i=1; i<=theMG->theBVPD.nSubDomains; i++)
    {
      mesh_2d.Element[i] = mesh_2d.start[i];
      for(j=0; j<mesh_2d.nElements[i]; j++)
      {

        if (MakeElement(theGrid, &(mesh_2d.Element[i]->elem_context)))
          return (15);
        mesh_2d.Element[i] = mesh_2d.Element[i]->succel;
      }
    }
  }

  SetStringValue(":gg:nElem",(double) theGrid->nElem);
  SetStringValue(":gg:nNode",(double) theGrid->nNode);

  if (display>0) UserWrite("\n");

  if (doAngle || doEdge)
    TerminateAccel(theMG, 0);
  return (0);
}

/****************************************************************************/
/*                                                                          */
/* initialization for this source file					                                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* initialization for grid generator library			                                */
/*                                                                          */
/****************************************************************************/

INT InitGG ()
{
  if (MakeStruct(":gg")!=0) return(__LINE__);

  return(0);
}
