// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugrefine3d.c													*/
/*																			*/
/* Purpose:   unstructured 3d-grid refinement (tree version)				*/
/*																			*/
/* Author:	  Juergen Bey													*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 294										*/
/*			  6900 Heidelberg												*/
/*																			*/
/*			  Klaus Johannsen												*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   24.06.92 begin, ug3 version 1.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment ugrefine
#endif

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "heaps.h"
#include "ugenv.h"

#include "devices.h"

#include "compiler.h"
#include "defaults.h"
#include "fileopen.h"
#include "switch.h"
#include "gm.h"
#include "GenerateRules.h"
#include "misc.h"
#include "evm.h"
#include "ugm.h"
#include "ugm3d.h"
#include "algebra.h"
#include "ugrefine.h"
#include "ugrefine3d.h"
#include "simplex.h"
#include "shapes3d.h"
#include "cw.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define RESOLUTION                      20              /* resolution for creating boundary midn*/
#define CENTER_NODE             10              /* position in struct ELEMENTCONTEXT whe*/
/* ptr to center node is stored                 */

#define MINVNCLASS              2                       /* determines copies, dep. on discr. !	*/

#define INNER_EDGE                      1
#define SIDE_EDGE                       2
#define HALF_FATHER_EDGE        3
#define FATHER_EDGE             4

/* macros concerning the MARK/REFINE flags */
#define REF_TYPE_CHANGES(e)             ((REFINE(e)!=MARK(e)) || (REFINECLASS(e)!=MARKCLASS(e)))
#define MARK_BISECT_EDGE(e,i)           ((Rules[MARK(e)].pattern)[i])
#define REFINE_BISECT_EDGE(e,i)         ((Rules[REFINE(e)].pattern)[i])

/* macros defining best refrule, specify exactly one of them !! */
/*#define __SHORTEST_INTERIOR_EDGE__*/
/*#define __MIDDLE_INTERIOR_EDGE__*/
#define __LONGEST_INTERIOR_EDGE__

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef NODE *ELEMENTCONTEXT[MAX_CORNERS_OF_ELEM+NEWCORNERS];

/****************************************************************************/
/*																			*/
/*	the refinement rules													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* edge type  : INNER_EDGE		 = inner edge of the father                             */
/*				SIDE_EDGE		 = inner edge of one side					*/
/*				HALF_FATHER_EDGE = half an edge of the father element		*/
/*				FATHER_EDGE      = edge of the father itself				*/
/*																			*/
/*																			*/
/****************************************************************************/

REFRULE *Rules;
static SHORT *PatternToRefrule;
static FULLREFRULEPTR theFullRefRule;
static int rFlag=GM_REFINE_TRULY_LOCAL; /* type of refine					*/
static INT theBFRRDirID;                                /* env type for BestFullRefRule		*/

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* Function:  ShowRefRule													*/
/*																			*/
/* Purpose:   fill SonList for theElement									*/
/*																			*/
/* Param:	  ELEMENT *theElement, ELEMENT *SonList[MAX_SONS]				*/
/*																			*/
/* return:	  0: ok                                                                                                                 */
/*			  1: error														*/
/*																			*/
/****************************************************************************/

static INT PrintEdgeData (struct edgedata theEdgeData)
{
  char buffer[128];

  sprintf(buffer,"   type=%d from=%d to=%d side=%d\n",(int)theEdgeData.type
          ,(int)theEdgeData.from
          ,(int)theEdgeData.to
          ,(int)theEdgeData.side);
  UserWrite(buffer);
  return(0);
}

static INT PrintSonData(struct sondata theSonData)
{
  char buffer[128];

  sprintf(buffer,"   corners=%d %d %d %d\n",(int)theSonData.corners[0]
          ,(int)theSonData.corners[1]
          ,(int)theSonData.corners[2]
          ,(int)theSonData.corners[3]);
  UserWrite(buffer);
  sprintf(buffer,"   nb=%d %d %d %d\n",(int)theSonData.nb[0]
          ,(int)theSonData.nb[1]
          ,(int)theSonData.nb[2]
          ,(int)theSonData.nb[3]);
  UserWrite(buffer);
  UserWrite("\n");
  return(0);
}

INT ShowRefRule (INT nb)
{
  char buffer[128];
  INT i,j;
  REFRULE *theRule;

  if (nb>241) return (1);

  theRule=&(Rules[nb]);

  /* header */
  UserWrite("\n");
  sprintf(buffer,"refrule %d:\n",nb);
  UserWrite(buffer);

  /* nsons */
  sprintf(buffer,"   nsons=%d\n",(int)theRule->nsons);
  UserWrite(buffer);

  /* pattern */
  UserWrite("   pattern=");
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
  {
    sprintf(buffer,"%d ",(int)theRule->pattern[i]);
    UserWrite(buffer);
  }
  UserWrite("\n");

  /* sonandnode */
  for (i=0; i<NEWCORNERS; i++)
  {
    sprintf(buffer,"   sonandnode[%d][0]=%d\n",i,(int)theRule->sonandnode[i][0]);
    UserWrite(buffer);
    sprintf(buffer,"             [%d][1]=%d\n",i,(int)theRule->sonandnode[i][1]);
    UserWrite(buffer);
  }
  UserWrite("\n");

  /* print edge data */
  UserWrite("edge data\n");
  for (i=0; i<MAXEDGES; i++)
    PrintEdgeData(theRule->edges[i]);
  UserWrite("\n");

  /* print sondata data */
  UserWrite("son data\n");
  for (i=0; i<(int)theRule->nsons; i++)
    PrintSonData(theRule->sons[i]);

  return (0);
}

/****************************************************************************/
/*
   GetSons - Fill SonList for theElement

   SYNOPSIS:
   INT GetSons (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS]);

   PARAMETERS:
   .  theElement -
   .  SonList[MAX_SONS]

   DESCRIPTION:
   This function fills SonList for theElement.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

INT GetSons (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS])
{
  REFRULE *theRule;
  ELEMENT *theSon;
  int SonID,PathPos;

  if (theElement==NULL) return (GM_ERROR);

  for (SonID=0; SonID<SONS_OF_ELEM(theElement); SonID++)
    SonList[SonID] = NULL;

  SonList[0] = SON(theElement,0);

  /* get other sons from path info in rules */
  theRule = &(Rules[REFINE(theElement)]);

  for (SonID=1; SonID<theRule->nsons; SonID++)
  {
    theSon = SonList[0];
    for (PathPos=0; PathPos<PATHDEPTH(theRule->sons[SonID].path); PathPos++)
      theSon = NBELEM(theSon,NEXTSIDE(theRule->sons[SonID].path,PathPos));

    if (theSon==NULL)
      return (GM_ERROR);

    SonList[SonID] = theSon;
  }

  return (GM_OK);
}

/****************************************************************************/
/*																			*/
/* Function:  ComputeCopies                                                                                             */
/*																			*/
/* Purpose:   determine copy elements from node classes                                         */
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  0: ok                                                                                                                 */
/*																			*/
/****************************************************************************/

static int ComputeCopies (GRID *theGrid)
{
  ELEMENT *theElement;
  int i,flag;

  /* set class of all dofs on next level to 0 */
  ClearNextVectorClasses(theGrid);

  /* seed dofs of regularly and irregularly refined elements to 3 */
  flag = 0;
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if (MARK(theElement)!=NOREFRULE && (MARKCLASS(theElement)==RED || MARKCLASS(theElement)==GREEN))
    {
      SeedNextVectorClasses(theElement);
      flag=1;                   /* there is at least one element to be refined */
    }

  /* copy all option or neighborhood */
  if (rFlag==GM_COPY_ALL)
  {
    if (flag)
      for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
        SeedNextVectorClasses(theElement);
  }
  else
  {
    PropagateNextVectorClasses(theGrid);
  }

  /* an element is copied if it has a dof of class 2 and higher */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    if ((MARK(theElement)==NOREFRULE)&&(MaxNextVectorClass(theElement)>=MINVNCLASS))
    {
      SETMARK(theElement,COPY_REFRULE);
      SETMARKCLASS(theElement,YELLOW);
    }

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  RestrictMarks                                                                                             */
/*																			*/
/* Purpose:   restrict refinement marks when going down                                         */
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void RestrictMarks (GRID *theGrid)
{
  ELEMENT *theElement,*SonList[MAX_SONS];
  EDGE *theEdge;
  int myClass,sonClass,myRule;
  int i,j,flag,CondensedPattern;

  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (GetSons(theElement,SonList)!=0) return;
    myRule = REFINE(theElement);
    sonClass = REFINECLASS (theElement);
    myClass = ECLASS(theElement);

    /* if element is not refined anyway */
    if (myRule==NOREFRULE) continue;

    /* copies always start with no refinement */
    if (myClass==YELLOW)
    {
      SETMARK(theElement,NOREFRULE);
      continue;
    }

    /* irregular elements are marked by estimator */
    if (myClass==GREEN) continue;

    /* regular elements with GREEN copies are marked by estimator */
    if (sonClass==YELLOW) continue;

    /* regular elements with GREEN or copy refinement go to no refinement or red refinement */
    if (sonClass==GREEN || sonClass==YELLOW)
    {
      for (i=0; i<NSONS(theElement); i++)
        if (MARK(SonList[i])>0)
        {
          if (MARKCLASS(theElement)==RED)
          {
            /* theElement is marked from outside */
            if (MARK(theElement)!=FULL_REFRULE)
              SETMARK(theElement,REFINE(theElement));
          }
          else
          {
            /* theElement is not marked from outside, so find a regular rule being consistent */
            /* with those neighbors of all sons of theElement which are marked for refine.	  */
            /* this choice will make sure these markes will not be distroyed.				  */
            CondensedPattern = Rules[REFINE(theElement)].pat;
            for (j=0; j<MAX_EDGES_OF_ELEM; j++)
            {
              theEdge=GetEdge(CORNER(theElement,CornerOfEdge[j][0]),CORNER(theElement,CornerOfEdge[j][1]));
              assert (theEdge != NULL);
              if (MIDNODE(theEdge)==NULL)
              {
                theEdge=GetEdge(SONNODE(CORNER(theElement,CornerOfEdge[j][0])),SONNODE(CORNER(theElement,CornerOfEdge[j][1])));
                assert(theEdge != NULL);
                if (ADDPATTERN(theEdge))
                  CondensedPattern |= (1<<j);
              }
            }
            SETMARK(theElement,PatternToRefrule[CondensedPattern]);
            SETMARKCLASS(theElement,RED);
          }
          break;
        }
      continue;
    }

    /* regular elements with regular refinement, are the only ones to coarsen */
    if (REFINECLASS(theElement) == RED)
    {
      SETMARK(theElement,REFINE(theElement));
      SETMARKCLASS(theElement,REFINECLASS(theElement));
    }
    flag = 0;
    for (i=0; i<NSONS(theElement); i++)
      if (!COARSEN(SonList[i]))
      {
        flag = 1;
        break;
      }

    /* at least one son has no coarsen flag set if flag==1 */
    if (flag) continue;
    /* remove refinement */

    SETMARK(theElement,NOREFRULE);
    SETCOARSEN(theElement,0);
  }
}

/****************************************************************************/
/*																			*/
/* Function:  CheckMemoryRequirements										*/
/*																			*/
/* Purpose:   check if there is enough memory for the following refinement	*/
/*																			*/
/* Param:	  MULTIGRID *theMG: pointer to multigrid structure with                 */
/*								computed refinement                                             */
/*																			*/
/* return:	  INT 1: enough memory for the operation						*/
/*			  INT 0: not enough memory for the operation					*/
/*																			*/
/****************************************************************************/

INT CheckMemoryRequirements (MULTIGRID *theMG)
{
  return(1);
}

/****************************************************************************/
/*																			*/
/* Function:  CreateMidNode                                                                                             */
/*																			*/
/* Purpose:   allocate a new node on an edge of an element. Includes vertex */
/*			  best fit boundary coordinates and local coordinates			*/
/*			  insert also links to endpoints of the refined edge			*/
/*																			*/
/* Param:	  ELEMENT *theElement: element to refine						*/
/*			  int edge: edge to refine										*/
/*			  NODE *after: insert new node after that node					*/
/*																			*/
/* return:	  NODE* : pointer to new node									*/
/*			  NULL	: could not allocate									*/
/*																			*/
/****************************************************************************/

static NODE *CreateMidNode (GRID *theGrid,ELEMENT *theElement,int edge,NODE *after)
{
  ELEMENTSIDE   *theSide;
  ELEMENT           *BoundaryElement;
  COORD s,smin;
  COORD x[3],r[3], ropt[3], l[2], dl[2], Inverse[9], Matrix[9];
  COORD_VECTOR HelpVector;
  COORD             *lambda, *lambda1, *lambda2, *cvect;
  int i,ni0,ni1,iopt;
  VERTEX            *theVertex, *v1, *v2;
  VSEGMENT          *vs1,*vs2, *vs;
  NODE              *theNode;
  BNDSEGDESC        *theSeg;

  /* calculate midpoint of edge */
  ni0 = CornerOfEdge[edge][0];
  ni1 = CornerOfEdge[edge][1];
  v1      = MYVERTEX(CORNER(theElement,ni0));
  v2      = MYVERTEX(CORNER(theElement,ni1));

  /* calculate physical position of MidNode */
  V3_LINCOMB(0.5, CVECT(v1), 0.5, CVECT(v2), x);

  /* check for boundary node */
  theVertex = NULL;
  if (OBJT(v1) == BVOBJ && OBJT(v2) == BVOBJ)
  {
    /* We now assume, that the edge is on the boundary if and only if	*/
    /* there is at least one boundary segment containing both vertices. */
    /* This works if we assume that there is no level 0 element                 */
    /* with a boundary side covering more than one boundary segment.	*/
    /* See also "InsertElementCommand" in "plot3d.c".					*/

    for (vs1=VSEG(v1); vs1!=NULL; vs1=NEXTSEG(vs1))
    {
      for (vs2=VSEG(v2); vs2!=NULL; vs2=NEXTSEG(vs2))
      {
        if (BSEGDESC(vs1) == BSEGDESC(vs2))
        {
          if (theVertex == NULL)
          {
            theVertex = CreateBoundaryVertex(theGrid,NULL);
            cvect = CVECT(theVertex);
          }
          if (theVertex == NULL) return(NULL);

          if ((vs = CreateVertexSegment(theGrid,theVertex)) == NULL)
          {
            DisposeVertex(theGrid, theVertex);
            return(NULL);
          }
          theSeg = BSEGDESC(vs) = BSEGDESC(vs1);

          /* find optimal parameters for this segment */
          lambda1 = PVECT(vs1);
          lambda2 = PVECT(vs2);
          lambda  = PVECT(vs);

          /* test midpoint */
          lambda[0] = 0.5*(lambda1[0] + lambda2[0]);
          lambda[1] = 0.5*(lambda1[1] + lambda2[1]);

          (*BNDSEGFUNC (theSeg))(SEGDATA(theSeg),lambda,ropt);

          V3_EUKLIDNORM_OF_DIFF(x,ropt,smin)

          if (smin>MAX_PAR_DIST)                               /* perhaps not the midpoint */
          {
            l[0] = lambda1[0];
            l[1] = lambda1[1];

            assert (RESOLUTION > 0);
            dl[0] = (lambda2[0] - lambda1[0])/((COORD) RESOLUTION);
            dl[1] = (lambda2[1] - lambda1[1])/((COORD) RESOLUTION);

            for (i=0; i<=RESOLUTION; i++)
            {
              (*BNDSEGFUNC (theSeg))(SEGDATA(theSeg),l,r);
              V3_EUKLIDNORM_OF_DIFF(x,r,s)
              if (s<smin)
              {
                smin = s;
                lambda[0] = l[0];
                lambda[1] = l[1];
                V3_COPY(r,ropt);
              }
              l[0]+=dl[0];
              l[1]+=dl[1];
            }
          }

          /* if it is the first vertex segment found fill in geometric data else compare with other vertex segments */
          if (NEXTSEG(vs) == NULL)
          {
            V3_COPY(ropt,cvect);
            SETMOVED(theVertex, smin>MAX_PAR_DIST);
          }
          else
          {
            V3_EUKLIDNORM_OF_DIFF(cvect,ropt,s)
            if (s>MAX_PAR_DIST)
            {
              DisposeVertex(theGrid,theVertex);
              UserWrite("two boundary segments with a common edge are not consistent\n");
              return(NULL);
            }
          }
        }
      }
    }
  }

  if (theVertex == NULL)
  {
    theVertex = CreateInnerVertex(theGrid,NULL);
    if (theVertex == NULL) return(NULL);
    V3_COPY(x,CVECT(theVertex));
  }

  VFATHER(theVertex) = theElement;
  SETONEDGE(theVertex,edge);

  /* create node */
  theNode = CreateNode(theGrid,after);
  if (theNode==NULL)
  {
    DisposeVertex(theGrid,theVertex);
    return(NULL);
  }
  MYVERTEX(theNode) = theVertex;
  NFATHER(theNode) = NULL;
  SETCLASS(theNode,4);

  /* local coordinates for the reference tetrahedron */
  if (!MOVED(theVertex))
  {
    XI(theVertex)  = 0.5*(TRefCoord[ni0][0]+TRefCoord[ni1][0]);
    ETA(theVertex) = 0.5*(TRefCoord[ni0][1]+TRefCoord[ni1][1]);
    NU(theVertex)  = 0.5*(TRefCoord[ni0][2]+TRefCoord[ni1][2]);
  }
  else
  {
    /* set transformation matrix */
    V3_SUBTRACT(CVECT(MYVERTEX(CORNER(theElement,1))), CVECT(MYVERTEX(CORNER(theElement,0))), Matrix)
    V3_SUBTRACT(CVECT(MYVERTEX(CORNER(theElement,2))), CVECT(MYVERTEX(CORNER(theElement,0))), Matrix+3)
    V3_SUBTRACT(CVECT(MYVERTEX(CORNER(theElement,3))), CVECT(MYVERTEX(CORNER(theElement,0))), Matrix+6)

    /* Invert it */
    if (M3_Invert(Inverse, Matrix))
    {
      UserWrite("error: cannot create MidNode on edge of the element\n");
      UserWrite("       the element is degenerated\n");
      return (NULL);
    }
    V3_SUBTRACT(CVECT(theVertex), CVECT(MYVERTEX(CORNER(theElement,0))), HelpVector)
    M3_TIMES_V3(Inverse,HelpVector,LCVECT(theVertex))
  }

  return(theNode);
}

/****************************************************************************/
/*																			*/
/* Function:  GetCurrentContext                                                                                         */
/*																			*/
/* Purpose:   assemble references to objects which interact with the sons	*/
/*			  of the given element, as indicated by REFINE.                                 */
/*			  (i)	 corner nodes											*/
/*			  (ii)	 nodes at midpoints of edges							*/
/*																			*/
/* Param:	  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: context structure to fill                 */
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void GetCurrentContext (ELEMENT *theElement, NODE **theElementContext)
{
  ELEMENT *SonList[MAX_SONS];
  INT i;
  REFRULE *rule;

  if (GetSons(theElement,SonList)!=0) return;

  /* get nodes, can be NULL */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    theElementContext[i] = SONNODE(CORNER(theElement,i));

  /* get midpoints */
  rule = &(Rules[REFINE(theElement)]);
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
  {
    if (rule->pattern[i])
      theElementContext[i+MAX_CORNERS_OF_ELEM] = CORNER(SonList[rule->sonandnode[i][0]],rule->sonandnode[i][1]);
    else
      theElementContext[i+MAX_CORNERS_OF_ELEM] = NULL;
  }

  /* get center node */
  if (rule->sonandnode[6][0] == NO_CENTER_NODE)
    theElementContext[CENTER_NODE] = NULL;
  else
    theElementContext[CENTER_NODE] = CORNER(SonList[rule->sonandnode[6][0]],rule->sonandnode[6][1]);
}


/****************************************************************************/
/*																			*/
/* Function:  UpdateContext                                                                                             */
/*																			*/
/* Purpose:   assemble references to objects which interact with the sons	*/
/*			  of the given element, i.e.									*/
/*			  objects are allocated, kept or deleted as indicated by MARK	*/
/*			  (i)	 corner nodes											*/
/*			  (ii)	 nodes at midpoints of edges							*/
/*																			*/
/* Param:	  GRID *theGrid: grid level of the sons of theElement			*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: context structure to update		*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*			  INT 1: fatal memory error                                                                     */
/*																			*/
/****************************************************************************/

static INT NodeNeededOnlyFor(NODE *theNode,NODE **theElementContext)
{
  LINK *theLink;
  INT i, found;

  for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
  {
    found=0;
    for (i=0; i<CENTER_NODE; i++)
      if (NBNODE(theLink)==theElementContext[i])
        found=1;
    if (!found)
      return (0);
  }

  return (1);
}

static int UpdateContext (GRID *theGrid, ELEMENT *theElement, NODE **theElementContext)
{
  NODE *theNode, *CenterNode;
  ELEMENT *theNeighbor,*theSon;                         /* neighbor and a son of current elem.	*/
  EDGE *theEdge;                                                        /* temporary storage for an edge		*/
  INT i,j,r,Corner0, Corner1,candelete;         /* some integer variables				*/
  NODE **MidNodes;                                                      /* nodes on refined edges				*/
  LINK *theLink;                                                        /* scan through nodes neighbor list     */
  NODE *Node0, *Node1;
  VERTEX *CenterVertex;
  EDGEDATA *edata;
  COORD *x, *xi;
  char buffer[64];


  /* allocate corner nodes if necessary */
  if (MARK(theElement)>0)
    for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    {
      theNode = CORNER(theElement,i);
      if (SONNODE(theNode)==NULL)
      {
        SONNODE(theNode) = CreateNode(theGrid,NULL);
        if (SONNODE(theNode)==NULL) return(1);
        MYVERTEX(SONNODE(theNode)) = MYVERTEX(theNode);
        NFATHER(SONNODE(theNode)) = theNode;
        theElementContext[i] = SONNODE(theNode);
      }
    }

  /* allocate,keep, or delete midpoint nodes */
  /* allocate,keep, or delete corner/corner edges */
  MidNodes = theElementContext+MAX_CORNERS_OF_ELEM;
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
  {
    Corner0 = CornerOfEdge[i][0];
    Corner1 = CornerOfEdge[i][1];
    if (MARK_BISECT_EDGE(theElement,i))
    {
      /* if a corner corner edge exists then delete it */
      if ((theEdge = GetEdge(theElementContext[Corner0],theElementContext[Corner1]))!=NULL)
        DisposeEdge(theGrid,theEdge);

      /* we need a midpoint node */
      if (MidNodes[i]!=NULL) continue;
      Node0 = CORNER(theElement,Corner0);
      Node1 = CORNER(theElement,Corner1);
      if ((theEdge = GetEdge(Node0,Node1))==NULL)
        return (1);
      MidNodes[i] = MIDNODE(theEdge);
      if (MidNodes[i] == NULL)
      {
        MidNodes[i] = CreateMidNode(theGrid,theElement,i,theElementContext[Corner0]);
        if (MidNodes[i]==NULL) return(1);
        MIDNODE(theEdge) = MidNodes[i];
        if ((theEdge=CreateEdge(theGrid,theElementContext[Corner0],MidNodes[i]))==NULL) return(1);
        if ((theEdge=CreateEdge(theGrid,theElementContext[Corner1],MidNodes[i]))==NULL) return(1);
      }
    }
    else
    {
      /* if we need a corner corner edge then allocate it */
      if (MARK(theElement)>0)
      {
        if ((theEdge=CreateEdge(theGrid,theElementContext[Corner0],theElementContext[Corner1]))==NULL) return(1);
      }

      /* we don't need a midpoint node on that edge, lets see if it can be deleted */
      if (MidNodes[i]==NULL) continue;
      candelete = 1;

      /* This midnode can be deleted, if all of it's remaining links */
      /* are to the endpoints of the father edge. In this case all   */
      /* elements sharing that edge have been visited.			   */
      for (theLink=START(MidNodes[i]); theLink!=NULL; theLink=NEXT(theLink))
        if ((NBNODE(theLink)!=theElementContext[Corner0])&&(NBNODE(theLink)!=theElementContext[Corner1]))
        {
          candelete = 0;
          break;
        }

      if (candelete)
      {
        DisposeVertex(theGrid,MYVERTEX(MidNodes[i]));
        DisposeNode(theGrid,MidNodes[i]);
        theEdge = GetEdge(CORNER(theElement,Corner0),CORNER(theElement,Corner1));
        MIDNODE(theEdge) = NULL;
      }
    }
  }

  /* delete corner nodes if possible */
  if (MARK(theElement)==0)
    for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
      if (theElementContext[i]!=NULL)
      {
        if (NodeNeededOnlyFor(theElementContext[i],theElementContext))
        {
          DisposeNode(theGrid,theElementContext[i]);
          SONNODE(CORNER(theElement,i)) = NULL;
        }
      }

  /* allocate/remove center node */
  if ((theElementContext[CENTER_NODE] != NULL) && (Rules[MARK(theElement)].sonandnode[6][0] == NO_CENTER_NODE))
  {
    /* existing center node has to be removed */
    DisposeVertex(theGrid,MYVERTEX(theElementContext[CENTER_NODE]));
    DisposeNode(theGrid,theElementContext[CENTER_NODE]);
    theElementContext[CENTER_NODE] = NULL;
  }
  if ((theElementContext[CENTER_NODE] == NULL) && (Rules[MARK(theElement)].sonandnode[6][0] != NO_CENTER_NODE))
  {
    /* there is no center node, but we need one: so allocate it */
    CenterNode = CreateNode(theGrid,NULL);
    if (CenterNode == NULL) return (1);

    /* allocate center vertex and init local and global position */
    CenterVertex = CreateInnerVertex(theGrid,NULL);
    if (CenterVertex == NULL) return (1);
    x = CVECT(CenterVertex);
    xi = LCVECT(CenterVertex);
    V3_CLEAR(x)
    V3_CLEAR(xi)
    for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    {
      V3_LINCOMB(1.0, x, 1.0, CVECT(MYVERTEX(CORNER(theElement,i))), x);
      V3_LINCOMB(1.0, xi, 1.0, TRefCoord[i], xi);
    }
    V3_SCALE(0.25, x)
    V3_SCALE(0.25, xi)

    /* init ptrs */
    VFATHER(CenterVertex) = theElement;
    TOPNODE(CenterVertex) = CenterNode;
    MYVERTEX(CenterNode) = CenterVertex;

    /* adjust element context */
    theElementContext[CENTER_NODE] = CenterNode;
  }

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  UnrefineElement												*/
/*																			*/
/* Purpose:   remove previous refinement of an element						*/
/*			  (i)	 all interior nodes and edges are deleted				*/
/*			  (ii)	 sons are deleted and references to sons reset to NULL	*/
/*																			*/
/* Param:	  GRID *theGrid: grid level of sons of theElement				*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: current context of element		*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static INT UnrefineElement (GRID *theGrid, ELEMENT *theElement, NODE **theElementContext)
{
  int i,j,s;
  EDGE *theEdge;
  NODE *CenterNode;
  REFRULE *rule;                                        /* current refinement rule of theElement*/
  ELEMENTCONTEXT sonContext;
  ELEMENT *theSon,*SonList[MAX_SONS];
  EDGEDATA *edata;

  /* something to do ? */
  if ((REFINE(theElement)==0)||(theGrid==NULL)) return(0);

  if (GetSons(theElement,SonList)!=0) return(1);

  for (s=0; s<NSONS(theElement); s++)
  {
    theSon = SonList[s];
    SETMARK(theSon,NOREFRULE);
    if (REFINE(theSon)>0)
    {
      GetCurrentContext(theSon,sonContext);
      if (UnrefineElement(theGrid->finer,theSon,sonContext)) return (1);
      UpdateContext(theGrid->finer,theSon,sonContext);
    }
  }

  /* remove connections in neighborhood of sons */
  for (i=0; i<NSONS(theElement); i++)
    DisposeConnectionsInNeighborhood(theGrid,SonList[i]);

  /* remove son elements */
  rule = &(Rules[REFINE(theElement)]);
  for (s=0; s<rule->nsons; s++)
  {
    if (DisposeEdgesFromElement(theGrid,SonList[s])) return (1);
    DisposeElement(theGrid,SonList[s]);
  }


  SETNSONS(theElement,0);
  SET_SON(theElement,0,NULL);

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  MarkForRefinement                                                                                         */
/*																			*/
/* Purpose:   mark an element for refinement								*/
/*																			*/
/* Param:	  ELEMENT *theElement: element to refine						*/
/*			  INT type: type of refinement mark:							*/
/*						REGULAR_MARK										*/
/*																			*/
/* return:	  INT 1: element has been marked								*/
/*				  0: element cannot be marked								*/
/*																			*/
/****************************************************************************/

INT IsAllowedToRefine (const ELEMENT *theElement)
{
  if (REFINECLASS(theElement)==RED) return (FALSE);
  return (TRUE);
}

INT MarkForRefinement (ELEMENT *theElement, INT rule, INT side)
{
  /* regulary refined elements can not be be marked */
  if (REFINECLASS(theElement)==RED) return(GM_ERROR);

  SETCOARSEN(theElement,0);
  switch(rule)
  {
  case (RED) :
    SETMARK(theElement,FULL_REFRULE);
    SETMARKCLASS(theElement,RED);
    break;
  case (COPY) :
    SETMARK(theElement,COPY_REFRULE);
    SETMARKCLASS(theElement,RED);
    break;
  case (NO_REFINEMENT) :
    SETMARK(theElement,NOREFRULE);
    break;
  case (UNREFINE) :
    SETMARK(theElement,NOREFRULE);
    SETCOARSEN(theElement,1);
    break;
  }

  return(GM_OK);
}

/****************************************************************************/
/*
   EstimateHere	- Return true (1) when element can be tagged for refinement

   SYNOPSIS:
   INT EstimateHere (ELEMENT *theElement);

   PARAMETERS:
   .  theElement - element to refine

   DESCRIPTION:
   This function returns true (1) when element can be tagged for refinement.

   RETURN VALUE:
   INT
   .n    false if do not tag element
   .n    true if element can be tagged for refinement.
 */
/****************************************************************************/

INT EstimateHere (ELEMENT *theElement)
{
  return(LEAFELEM(theElement));
}

/****************************************************************************/
/*																*/
/* Function:  GetRefinementMark                                                                                         */
/*																			*/
/* Purpose:   gets rule and variant of refinement							*/
/*																			*/
/* Param:	  ELEMENT *theElement: element to refine						*/
/*			  int *rule: filled with current refinement rule				*/
/*			  int *side: filled with side, if rule is oriented				*/
/*																			*/
/* return:	  int 0: side information valid                                                                 */
/*			  int 1: rule without orientation								*/
/*																			*/
/****************************************************************************/

INT GetRefinementMark (const ELEMENT *theElement, INT *rule, INT *side)
{
  if (MARK(theElement)==FULL_REFRULE) {*rule=RED; return(GM_RULE_WITHOUT_ORIENTATION);}
  switch (MARK(theElement))
  {
  case NOREFRULE :
    *rule=NO_REFINEMENT;
    if (COARSEN(theElement)) *rule = UNREFINE;
    break;
  case COPY_REFRULE : *rule=COPY; break;
  default : *rule=NO_REFINEMENT;  break;
  }
  *side=0;
  return(GM_RULE_WITHOUT_ORIENTATION);
}

/****************************************************************************/
/*																			*/
/* Function:  RefineElement                                                                                             */
/*																			*/
/* Purpose:   refine an element in the given context						*/
/*			  (i)	 corner and midnodes are already allocated				*/
/*			  (ii)	 edges between corner and midnodes are ok				*/
/*			  (iii)  create interior nodes and edges						*/
/*			  (iv)	 create sons and set references to sons                                 */
/*																			*/
/* Param:	  GRID *theGrid: grid level of sons of theElement				*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: current context of element		*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*			  INT 1: fatal memory error                                                                     */
/*																			*/
/****************************************************************************/

static int RefineElement (GRID *theGrid, ELEMENT *theElement, NODE **theElementContext)
{
  INT i,j,l,s,s2,ni,pi,nsi,p,q,side;
  INT points,ni0,ni1,m;
  ELEMENT *theSon,*theNeighbor;
  EDGE *theEdge;
  ELEMENT *SonList[MAX_SONS],*SonList2[MAX_SONS];
  ELEMENTSIDE *oldSide,*newSide;
  VERTEX *myvertex;
  VSEGMENT *vs;
  NODE *mynode;
  COORD *v0,*v1,*v2, vd[3];
  COORD *x;
  INT boundaryelement, found;
  REFRULE *rule, *rule2;
  EDGEDATA *edata;
  SONDATA *sdata, *sdata2;

  /* is something to do ? */
  if (MARK(theElement)==0) return(0);
  rule = &(Rules[MARK(theElement)]);

  /* create interior edges */
  for (s=0; s<MAXEDGES; s++)
  {
    edata = &(rule->edges[s]);
    if (edata->type != INNER_EDGE) continue;
    if ((theEdge=CreateEdge(theGrid,theElementContext[edata->from], theElementContext[edata->to]))==NULL)
      return(1);
  }

  /* create elements */
  for (s=0; s<rule->nsons; s++)
  {
    boundaryelement = 0;
    if (OBJT(theElement) == BEOBJ)
      for (i=0; i<MAX_SIDES_OF_ELEM; i++)
        if ( (side = rule->sons[s].nb[i]) >= FATHER_SIDE_OFFSET )                                  /* exterior side */
          if (SIDE(theElement,side-FATHER_SIDE_OFFSET)!=NULL)                                           /* at the boundary */
          {
            boundaryelement = 1;
            break;
          }
    if (boundaryelement)
      theSon = CreateBoundaryElement(theGrid,NULL,TETRAHEDRON);
    else
      theSon = CreateInnerElement(theGrid,NULL,TETRAHEDRON);
    if (theSon==NULL) return(1);

    /* fill in son data */
    SonList[s] = theSon;
    SETECLASS(theSon,MARKCLASS(theElement));
    SET_EFATHER(theSon,theElement);
    for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
      SET_CORNER(theSon,i,theElementContext[rule->sons[s].corners[i]]);
    for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
      for (j=i+1; j<MAX_CORNERS_OF_ELEM; j++)
      {
        theEdge = CreateEdge(theGrid,CORNER(theSon,i),CORNER(theSon,j));
        assert(theEdge!=NULL);
        if (NO_OF_ELEM(theEdge)<NO_OF_ELEM_MAX-1)
          INC_NO_OF_ELEM(theEdge);
        else
          return (1);
      }
  }
  SETNSONS(theElement,rule->nsons);
  SET_SON(theElement,0,SonList[0]);

  /* create element sides at the boundary */
  if (OBJT(theElement)==BEOBJ)
    for( s=0; s<rule->nsons; s++ )
    {
      if (OBJT(SonList[s]) != BEOBJ) continue;
      for (j=0; j<MAX_SIDES_OF_ELEM; j++)
      {
        SET_SIDE(SonList[s],j,NULL);
        if ((side = rule->sons[s].nb[j]) < FATHER_SIDE_OFFSET) continue;
        side -= FATHER_SIDE_OFFSET;
        if ((oldSide = SIDE(theElement,side)) == NULL) continue;

        newSide = CreateElementSide(theGrid);
        if (newSide==NULL) return(1);
        SET_SIDE(SonList[s],j,newSide);
        SEGDESC(newSide) = SEGDESC(oldSide);

        for (i=0; i<MAX_CORNERS_OF_SIDE; i++)
        {
          ni = CornerOfSide[j][i];                                                                        /* node index of son	 */
          pi = rule->sons[s].corners[ni];                                                         /* point in theElement */

          if (pi<MAX_CORNERS_OF_ELEM)                                                       /* a corner of theElement */
          {
            nsi = CornerOfSideInv[side][pi];                                                             /* position of corner in side numeration */
            assert(nsi != -1);
            PARAM(newSide,i,0) = PARAM(oldSide,nsi,0);
            PARAM(newSide,i,1) = PARAM(oldSide,nsi,1);
          }
          else                                                       /* a midpoint of an edge of theElement */
          {
            myvertex = MYVERTEX(theElementContext[pi]);
            assert (OBJT(myvertex) == BVOBJ);
            assert (VSEG(myvertex) != NULL);

            /* find common boundary segment */
            for( vs=VSEG(myvertex); vs!=NULL; vs = NEXTSEG(vs) )
            {
              if (BSEGDESC(vs) == SEGDESC(oldSide)) break;
            }
            assert(vs!=NULL);
            PARAM(newSide,i,0) =  LAMBDA(vs,0);
            PARAM(newSide,i,1) =  LAMBDA(vs,1);
          }
        }
      }
    }

  /* connect elements */
  for (s=0; s<rule->nsons; s++)
  {
    sdata = &(rule->sons[s]);
    for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    {
      SET_NBELEM(SonList[s],i,NULL);

      /* an interior triangle face */
      if ( (side = sdata->nb[i]) < FATHER_SIDE_OFFSET )
      {
        SET_NBELEM(SonList[s],i,SonList[side]);

        /* dispose doubled side vectors if */
                                #ifdef __SIDEDATA__
        for (l=0; l<MAX_SIDES_OF_ELEM; l++)
          if (rule->sons[side].nb[l]==s)
            break;
        if (DisposeDoubledSideVector(theGrid,SonList[s],i,SonList[side],l))
          return (1);
                                #endif
        continue;
      }

      /* the boundary case */
      if ((OBJT(SonList[s]) == BEOBJ) && (SIDE(SonList[s],i) != NULL)) continue;

      /* check, if neighbor has been refined */
      side -= FATHER_SIDE_OFFSET;
      theNeighbor = NBELEM(theElement,side);
      if(theNeighbor==NULL) continue;

      if (REF_TYPE_CHANGES(theNeighbor) || (REFINE(theNeighbor) == 0 ))
        continue;

      for (l=0; l<MAX_SIDES_OF_ELEM; l++)
        if (NBELEM(theNeighbor,l) == theElement)
          break;
      assert(l<MAX_SIDES_OF_ELEM);

      rule2 = &(Rules[REFINE(theNeighbor)]);
      found = 0;
      if (GetSons(theNeighbor,SonList2)!=0)
        return (1);
      for (s2=0; s2<rule2->nsons; s2++)
      {
        sdata2 = &(rule2->sons[s2]);
        for (j=0; j<MAX_SIDES_OF_ELEM; j++)
        {
          if (sdata2->nb[j] != FATHER_SIDE_OFFSET+l) continue;
          points=0;
          for (p=0; p<MAX_CORNERS_OF_SIDE; p++)
            for (q=0; q<MAX_CORNERS_OF_SIDE; q++)
              if (CORNER(SonList[s],CornerOfSide[i][p]) == CORNER(SonList2[s2],CornerOfSide[j][q]))
              {
                points |= ((1<<p) + (8<<q));
                break;
              }
          if (points == 63)                                                      /* neighbor found */
          {
            /* adjust pointers */
            SET_NBELEM(SonList[s],i,SonList2[s2]);
            SET_NBELEM(SonList2[s2],j,SonList[s]);

            /* dispose doubled side vectors if */
                                                #ifdef __SIDEDATA__
            if (DisposeDoubledSideVector(theGrid,SonList[s],i,SonList2[s2],j))
              return (1);
                                                #endif

            found=1;
            break;
          }
        }
        if (found) break;
      }
      assert (found==1);
    }
  }

  return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  RefineGrid													*/
/*																			*/
/* Purpose:   refine one level of the grid									*/
/*																			*/
/* Param:	  GRID *theGrid: grid level to refine							*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*			  INT 1: fatal memory error                                                                     */
/*																			*/
/****************************************************************************/

static int RefineGrid (GRID *theGrid)
{
  ELEMENT *theElement;
  ELEMENTCONTEXT theContext;
  GRID *fineGrid;
  NODE *theNode;

  fineGrid = theGrid->finer;
  if (fineGrid==NULL) return(1);

  /* refine elements */
  RESETGSTATUS(fineGrid,GRID_CHANGED);
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (REF_TYPE_CHANGES(theElement))
    {
      GetCurrentContext(theElement,theContext);
      if (UnrefineElement(fineGrid,theElement,theContext)) return(1);
      if (UpdateContext(fineGrid,theElement,theContext)!=0) return(1);
      if (RefineElement(fineGrid,theElement,theContext)!=0) return(1);
      SETREFINE(theElement,MARK(theElement));
      SETREFINECLASS(theElement,MARKCLASS(theElement));
      SETGSTATUS(fineGrid,GRID_CHANGED);
    }
  }

  /* reset newref flags */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    SETMARK(theElement,NOREFRULE);
  }

  /* set node class on next level */
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
    if (SONNODE(theNode)!=NULL)
    {
      SETCLASS(SONNODE(theNode),NCLASS(theNode));
      if (NCLASS(theNode)>=2) TOPNODE(MYVERTEX(theNode)) = SONNODE(theNode);
    }

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  YAlignment													*/
/*																			*/
/* Purpose:   compute best full refined refrule for the element                         */
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                                 */
/*																			*/
/* return:	  INT Mark: number of refrule									*/
/*																			*/
/****************************************************************************/

static INT YAlignment (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM], help;
  INT i, flags, imax;
  COORD Dist_0_5, Dist_1_3, Dist_2_4, max;

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* compute differences */
  max=-1.0;
  for (i=0; i<EDGES_OF_ELEM(theElement); i++)
  {
    V3_SUBTRACT(Corners[CornerOfEdge[i][0]], Corners[CornerOfEdge[i][1]], help)
    V3_Normalize(help);
    if (ABS(help[1]+help[2])>max) {imax=i; max=ABS(help[1]+help[2]);}
  }
  V3_EUKLIDNORM_OF_DIFF(MidPoints[0], MidPoints[5], Dist_0_5)
  V3_EUKLIDNORM_OF_DIFF(MidPoints[1], MidPoints[3], Dist_1_3)
  V3_EUKLIDNORM_OF_DIFF(MidPoints[2], MidPoints[4], Dist_2_4)
  switch (imax)
  {
  case 0 : if (Dist_1_3<Dist_2_4) return (FULL_REFRULE_1_3);
    else return (FULL_REFRULE_2_4);
  case 1 : if (Dist_0_5<Dist_2_4) return (FULL_REFRULE_0_5);
    else return (FULL_REFRULE_2_4);
  case 2 : if (Dist_1_3<Dist_0_5) return (FULL_REFRULE_1_3);
    else return (FULL_REFRULE_0_5);
  case 3 : if (Dist_0_5<Dist_2_4) return (FULL_REFRULE_0_5);
    else return (FULL_REFRULE_2_4);
  case 4 : if (Dist_1_3<Dist_0_5) return (FULL_REFRULE_1_3);
    else return (FULL_REFRULE_0_5);
  case 5 : if (Dist_1_3<Dist_2_4) return (FULL_REFRULE_1_3);
    else return (FULL_REFRULE_2_4);
  }
}

/****************************************************************************/
/*																			*/
/* Function:  ShortestInteriorEdge											*/
/*																			*/
/* Purpose:   compute best full refined refrule for the element                         */
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                                 */
/*																			*/
/* return:	  INT Mark: number of refrule									*/
/*																			*/
/****************************************************************************/

static INT ShortestInteriorEdge (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM];
  INT i, flags;
  COORD Dist_0_5, Dist_1_3, Dist_2_4;

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* compute distances */
  V3_EUKLIDNORM_OF_DIFF(MidPoints[0], MidPoints[5], Dist_0_5)
  V3_EUKLIDNORM_OF_DIFF(MidPoints[1], MidPoints[3], Dist_1_3)
  V3_EUKLIDNORM_OF_DIFF(MidPoints[2], MidPoints[4], Dist_2_4)

  /* return best Refrule (shortest interior edge) */
  flags = (Dist_0_5 < Dist_1_3) ;
  flags |= ((Dist_1_3 < Dist_2_4) <<1) ;
  flags |= ((Dist_2_4 < Dist_0_5) <<2) ;
  assert(flags != 7);

  switch(flags)
  {
  case 0 :                                              /* Dist_0_5 = Dist_2_4 = Dist_1_3 */
    return (FULL_REFRULE_0_5);
  case 1 :                                              /* Dist_0_5 < Dist_2_4 < Dist_1_3 */
    return (FULL_REFRULE_0_5);
  case 2 :                                              /* Dist_1_3 < Dist_0_5 < Dist_2_4 */
    return (FULL_REFRULE_1_3);
  case 3 :                                              /* Dist_0_5 < Dist_1_3 < Dist_2_4 */
    return (FULL_REFRULE_0_5);
  case 4 :                                              /* Dist_2_4 < Dist_1_3 < Dist_0_5 */
    return (FULL_REFRULE_2_4);
  case 5 :                                              /* Dist_2_4 < Dist_0_5 < Dist_1_3 */
    return (FULL_REFRULE_2_4);
  case 6 :                                              /* Dist_1_3 < Dist_2_4 < Dist_0_5 */
    return (FULL_REFRULE_1_3);
  }
}

/****************************************************************************/
/*																			*/
/* Function:  MinimalSideAngle													*/
/*																			*/
/* Purpose:   compute best full refined refrule for the element                         */
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                                 */
/*																			*/
/* return:	  INT Mark: number of refrule									*/
/*																			*/
/****************************************************************************/

static INT MinimalSideAngle (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM];
  INT i,j,k,l,imin;
  COORD MaxAngle,Max,Min;

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* try possebilities */
  Min = 190.0;
  for (i=0; i<3; i++)
  {
    j = OppositeEdge[i];
    Corners[2] = MidPoints[i];
    Corners[3] = MidPoints[j];

    Max = 0.0;
    for (k=0; k<2; k++)
    {
      for (l=0; l<2; l++)
        Corners[l] = MidPoints[SideEdgesOfEdge[i][k][l]];
      if (TetMaxSideAngle (Corners,&MaxAngle))
        return (FULL_REFRULE_0_5);
      Max = MAX(Max,MaxAngle);
    }
    for (k=0; k<2; k++)
    {
      for (l=0; l<2; l++)
        Corners[l] = MidPoints[SideEdgesOfEdge[j][k][l]];
      if (TetMaxSideAngle (Corners,&MaxAngle))
        return (FULL_REFRULE_0_5);
      Max = MAX(Max,MaxAngle);
    }
    if (Max<Min)
    {
      Min = Max;
      imin = i;
    }
  }

  switch (imin)
  {
  case 0 :
    return (FULL_REFRULE_0_5);
  case 1 :
    return (FULL_REFRULE_1_3);
  case 2 :
    return (FULL_REFRULE_2_4);
  }
}

/****************************************************************************/
/*																			*/
/* Function:  MinimalSideEntry													*/
/*																			*/
/* Purpose:   compute best full refined refrule for the element                         */
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                                 */
/*																			*/
/* return:	  INT Mark: number of refrule									*/
/*																			*/
/****************************************************************************/

static INT MinimalSideEntry (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM];
  INT i,j,k,l,imin;
  COORD Angle[MAX_EDGES_OF_ELEM],Length[MAX_EDGES_OF_ELEM],Max,Min,help;

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* try possebilities */
  Min = MAX_C;
  for (i=0; i<3; i++)
  {
    j = OppositeEdge[i];
    Corners[2] = MidPoints[i];
    Corners[3] = MidPoints[j];

    Max = 0.0;
    for (k=0; k<2; k++)
    {
      for (l=0; l<2; l++)
        Corners[l] = MidPoints[SideEdgesOfEdge[i][k][l]];
      if (TetAngleAndLength (Corners,Angle,Length))
        return (FULL_REFRULE_0_5);
      for (l=0; l<MAX_EDGES_OF_ELEM; l++)
        if (Angle[l]>PI/2.0)
        {
          help = ABS(Length[l]*(COORD)(cos((double)Angle[l])/sin((double)Angle[l])));
          Max = MAX(Max,help);
        }
    }
    for (k=0; k<2; k++)
    {
      for (l=0; l<2; l++)
        Corners[l] = MidPoints[SideEdgesOfEdge[j][k][l]];
      if (TetAngleAndLength (Corners,Angle,Length))
        return (FULL_REFRULE_0_5);
      for (l=0; l<MAX_EDGES_OF_ELEM; l++)
        if (Angle[l]>PI/2.0)
        {
          help = ABS(Length[l]*(COORD)(cos((double)Angle[l])/sin((double)Angle[l])));
          Max = MAX(Max,help);
        }
    }
    if (Max<Min)
    {
      Min = Max;
      imin = i;
    }
  }

  switch (imin)
  {
  case 0 :
    return (FULL_REFRULE_0_5);
  case 1 :
    return (FULL_REFRULE_1_3);
  case 2 :
    return (FULL_REFRULE_2_4);
  }
}

/****************************************************************************/
/*																			*/
/* Function:  BestLaplaceMMatrix										*/
/*																			*/
/* Purpose:   compute best full refined refrule for the element:		*/
/*			  optimal laplace-disc w.r.t. M-Matrix eigenschaft				*/
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                         */
/*																			*/
/* return:	  INT Mark: number of refrule								*/
/*																			*/
/****************************************************************************/

static INT BestLaplaceMMatrix (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM];
  INT i,j,k,l,imin,TBFR,refrule;
  COORD Angle[MAX_EDGES_OF_ELEM],Length[MAX_EDGES_OF_ELEM],sum,Min;
  char buffer[64];

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* try possebilities */
  Min = MAX_C; imin = -1;
  for (i=0; i<3; i++)
  {
    j = OppositeEdge[i];
    Corners[2] = MidPoints[i];
    Corners[3] = MidPoints[j];

    sum = 0.0;
    for (k=0; k<2; k++)
    {
      for (l=0; l<2; l++)
        Corners[l] = MidPoints[SideEdgesOfEdge[i][k][l]];
      if (TetAngleAndLength (Corners,Angle,Length))
        return (FULL_REFRULE_0_5);
      for (l=0; l<MAX_EDGES_OF_ELEM; l++)
        if (Angle[l]>PI/2.0)
          sum += ABS(Length[l]*(COORD)cos((double)Angle[l])/(COORD)sin((double)Angle[l]));
    }
    for (k=0; k<2; k++)
    {
      for (l=0; l<2; l++)
        Corners[l] = MidPoints[SideEdgesOfEdge[j][k][l]];
      if (TetAngleAndLength (Corners,Angle,Length))
        return (FULL_REFRULE_0_5);
      for (l=0; l<MAX_EDGES_OF_ELEM; l++)
        if (Angle[l]>PI/2.0)
          sum += ABS(Length[l]*(COORD)cos((double)Angle[l])/(COORD)sin((double)Angle[l]));
    }
    if (sum<Min)
    {
      Min = sum;
      imin = i;
    }
  }

  switch (imin)
  {
  case 0 :
    refrule = FULL_REFRULE_0_5;
    break;
  case 1 :
    refrule = FULL_REFRULE_1_3;
    break;
  case 2 :
    refrule = FULL_REFRULE_2_4;
    break;
  }


  TBFR = ShortestInteriorEdge(theElement);
  if (imin == -1)
  {
    refrule = TBFR;
    UserWrite ("#");
  }

  return (refrule);
}

/****************************************************************************/
/*																			*/
/* Function:  MaxPerpendicular											*/
/*																			*/
/* Purpose:   compute best full refined refrule for the element:		*/
/*			  optimal laplace-disc w.r.t. M-Matrix eigenschaft				*/
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                         */
/*																			*/
/* return:	  INT Mark: number of refrule								*/
/*																			*/
/****************************************************************************/

static INT MaxPerpendicular (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM],a,b,c;
  INT i,j,k,l,imin,TBFR,refrule;
  COORD Angle[MAX_EDGES_OF_ELEM],Length[MAX_EDGES_OF_ELEM],sprd,Max;
  char buffer[64];

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* try possebilities */
  Max = -MAX_C; imin = -1;
  for (i=0; i<3; i++)
  {
    j = OppositeEdge[i];

    V3_SUBTRACT(Corners[CornerOfEdge[i][0]],Corners[CornerOfEdge[i][1]],a)
    V3_SUBTRACT(Corners[CornerOfEdge[j][0]],Corners[CornerOfEdge[j][1]],b)
    V3_VECTOR_PRODUCT(a,b,c)
    V3_Normalize(c);
    V3_SUBTRACT(MidPoints[i],MidPoints[j],a)
    V3_Normalize(a);
    V3_SCALAR_PRODUCT(a,c,sprd)
    sprd = ABS(sprd);

    if (sprd>Max)
    {
      Max = sprd;
      imin = i;
    }
  }

  switch (imin)
  {
  case 0 :
    refrule = FULL_REFRULE_0_5;
    break;
  case 1 :
    refrule = FULL_REFRULE_1_3;
    break;
  case 2 :
    refrule = FULL_REFRULE_2_4;
    break;
  }


  TBFR = ShortestInteriorEdge (theElement);
  if (imin == -1)
  {
    refrule = TBFR;
    UserWrite ("#");
  }

  return (refrule);
}

/****************************************************************************/
/*																			*/
/* Function:  MaxRightAngle                                                                                     */
/*																			*/
/* Purpose:   compute best full refined refrule for the element:		*/
/*			  optimal laplace-disc w.r.t. M-Matrix eigenschaft				*/
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                         */
/*																			*/
/* return:	  INT Mark: number of refrule								*/
/*																			*/
/****************************************************************************/

static INT MaxRightAngle (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM],a,b,c;
  INT i,j,k,l,imin,TBFR,refrule;
  COORD Angle[MAX_EDGES_OF_ELEM],Length[MAX_EDGES_OF_ELEM],sprd,Min;
  char buffer[64];

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* try possebilities */
  Min = MAX_C; imin = -1;
  for (i=0; i<3; i++)
  {
    j = OppositeEdge[i];

    V3_SUBTRACT(Corners[CornerOfEdge[i][0]],Corners[CornerOfEdge[i][1]],a)
    V3_Normalize(a);
    V3_SUBTRACT(Corners[CornerOfEdge[j][0]],Corners[CornerOfEdge[j][1]],b)
    V3_Normalize(b);
    V3_SCALAR_PRODUCT(a,b,sprd)
    sprd = ABS(sprd);

    if (sprd<Min)
    {
      Min = sprd;
      imin = i;
    }
  }

  switch (imin)
  {
  case 0 :
    refrule = FULL_REFRULE_0_5;
    break;
  case 1 :
    refrule = FULL_REFRULE_1_3;
    break;
  case 2 :
    refrule = FULL_REFRULE_2_4;
    break;
  }


  TBFR = ShortestInteriorEdge (theElement);
  if (imin == -1)
  {
    refrule = TBFR;
    UserWrite ("#");
  }

  return (refrule);
}

/****************************************************************************/
/*																			*/
/* Function:  MaxArea												*/
/*																			*/
/* Purpose:   compute best full refined refrule for the element:		*/
/*			  optimal laplace-disc w.r.t. M-Matrix eigenschaft				*/
/*																			*/
/* Param:	  ELEMENT *theElement: for that element                                         */
/*																			*/
/* return:	  INT Mark: number of refrule								*/
/*																			*/
/****************************************************************************/

static INT MaxArea (ELEMENT *theElement)
{
  COORD *Corners[MAX_CORNERS_OF_ELEM];
  COORD_VECTOR MidPoints[MAX_EDGES_OF_ELEM],a,b,c;
  INT i,j,k,l,imin,TBFR,refrule;
  COORD Angle[MAX_EDGES_OF_ELEM],Length[MAX_EDGES_OF_ELEM],norm,Max;
  char buffer[64];

  /* get physical position of the corners */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    Corners[i] = CVECT(MYVERTEX(CORNER(theElement,i)));

  /* get physical position of the midpoints of the edges */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    V3_LINCOMB(0.5, Corners[CornerOfEdge[i][0]], 0.5, Corners[CornerOfEdge[i][1]], MidPoints[i]);

  /* try possebilities */
  Max = -MAX_C; imin = -1;
  for (i=0; i<3; i++)
  {
    j = OppositeEdge[i];

    V3_SUBTRACT(Corners[CornerOfEdge[i][0]],Corners[CornerOfEdge[i][1]],a)
    V3_SUBTRACT(Corners[CornerOfEdge[j][0]],Corners[CornerOfEdge[j][1]],b)
    V3_VECTOR_PRODUCT(a,b,c)
    V3_EUKLIDNORM(c,norm)

    if (norm>Max)
    {
      Max = norm;
      imin = i;
    }
  }

  switch (imin)
  {
  case 0 :
    refrule = FULL_REFRULE_0_5;
    break;
  case 1 :
    refrule = FULL_REFRULE_1_3;
    break;
  case 2 :
    refrule = FULL_REFRULE_2_4;
    break;
  }


  TBFR = ShortestInteriorEdge (theElement);
  if (imin == -1)
  {
    refrule = TBFR;
    UserWrite ("#");
  }

  return (refrule);
}

/****************************************************************************/
/*																			*/
/* Function:  CloseGrid                                                                                                         */
/*																			*/
/* Purpose:   compute closure for next level								*/
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  INT >0: elements will be refined								*/
/*			  INT 0: no elements will be refined							*/
/*																			*/
/****************************************************************************/

static int CloseGrid (GRID *theGrid)
{
  ELEMENT *theElement, *NbElement;
  EDGE *MyEdge, *NbEdge;
  INT i, j, cnt;
  INT Mark, MyEdgePattern, MySidePattern, MyEdgeNum;
  INT NbEdgePattern, NbSidePattern, NbEdgeNum, NbSideMask;
  SHORT *myPattern;

  /* reset pattern and used flag on the edges */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    SETUSED(theElement,0);
    for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
      for (j=i+1; j<MAX_CORNERS_OF_ELEM; j++)
      {
        MyEdge=GetEdge(CORNER(theElement,i),CORNER(theElement,j));
        assert (MyEdge != NULL);
        SETPATTERN(MyEdge,0);
        SETADDPATTERN(MyEdge,1);
      }
  }

  /* set pattern on the edges and reset side-and edge pattern field */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (MARKCLASS(theElement)==RED)
    {
      Mark = MARK(theElement);
      myPattern = Rules[Mark].pattern;
      for (i=0; i<MAX_EDGES_OF_ELEM; i++)
        if (myPattern[i])
        {
          MyEdge=GetEdge(CORNER(theElement,CornerOfEdge[i][0]),CORNER(theElement,CornerOfEdge[i][1]));
          if (MyEdge != NULL)
            SETPATTERN(MyEdge,1);
        }
    }
    SETEDGEPATTERN(theElement,0);
    SETSIDEPATTERN(theElement,0);
  }

  /* set pattern (edge and side) on the elements */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    /* make edgepattern consistent with pattern of edges */
    SETUSED(theElement,1);
    MyEdgePattern = 0;
    for (i=MAX_EDGES_OF_ELEM-1; i>=0; i--)
    {
      MyEdge=GetEdge(CORNER(theElement,CornerOfEdge[i][0]),CORNER(theElement,CornerOfEdge[i][1]));
      MyEdgePattern = (MyEdgePattern<<1) | PATTERN(MyEdge);
    }
    SETEDGEPATTERN(theElement,MyEdgePattern);

    /* make SIDEPATTERN consistent with neighbors */
    for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    {
      NbElement = NBELEM(theElement,i);
      if (NbElement == NULL) continue;
      if (!USED(NbElement)) continue;
      /* now edgepattern from theelement and NbElement are in final state */
      for (j=0; j<MAX_SIDES_OF_ELEM; j++)
        if (NBELEM(NbElement,j) == theElement)
          break;

      /* because SETSIDEPATTERN is set to zero, I choose TriSectionEdgeÉÉ[0] */
      MyEdgeNum = TriSectionEdge[MyEdgePattern&CondensedEdgeOfSide[i]][0];
      if (MyEdgeNum == -2) return (0);
      if (MyEdgeNum == -1) continue;

      NbEdgePattern = EDGEPATTERN(NbElement);
      NbEdgeNum = TriSectionEdge[NbEdgePattern&CondensedEdgeOfSide[j]][0];
      if (NbEdgeNum == -2 || NbEdgeNum == -1) return (0);

      if (!( CORNER(theElement,CornerOfEdge[MyEdgeNum][0]) == CORNER(NbElement,CornerOfEdge[NbEdgeNum][0]) &&
             CORNER(theElement,CornerOfEdge[MyEdgeNum][1]) == CORNER(NbElement,CornerOfEdge[NbEdgeNum][1])        ) &&
          !( CORNER(theElement,CornerOfEdge[MyEdgeNum][0]) == CORNER(NbElement,CornerOfEdge[NbEdgeNum][1]) &&
             CORNER(theElement,CornerOfEdge[MyEdgeNum][1]) == CORNER(NbElement,CornerOfEdge[NbEdgeNum][0])        )        )
      {
        NbSidePattern = SIDEPATTERN(NbElement);
        NbSideMask = (1<<j);
        if ( NbSidePattern & NbSideMask )
          NbSidePattern &= ~NbSideMask;
        else
          NbSidePattern |= NbSideMask;
        SETSIDEPATTERN(NbElement,NbSidePattern);
      }
    }
  }

  /* set refinement rules from edge- and sidepattern */
  cnt = 0;
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    MyEdgePattern = EDGEPATTERN(theElement);
    MySidePattern = SIDEPATTERN(theElement);
    Mark = PatternToRefrule[MyEdgePattern | (MySidePattern<<MAX_EDGES_OF_ELEM)];
    assert(Mark != -1);
    if (Mark == FULL_REFRULE)
      Mark = (*theFullRefRule)(theElement);
    if (MARKCLASS(theElement)==RED && Mark==NOREFRULE)
      Mark = COPY_REFRULE;
    if (MARKCLASS(theElement)!=RED && Mark!=NOREFRULE)
      SETMARKCLASS(theElement,GREEN);
    if (Mark)
      cnt++;
    SETMARK(theElement,Mark);
  }

  /* build a green covering around the red elements */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (MARKCLASS(theElement)!=RED) continue;
    for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    {
      NbElement = NBELEM(theElement,i);
      if (NbElement == NULL) continue;
      if ((MARKCLASS(NbElement)==GREEN) || (MARKCLASS(NbElement)==RED)) continue;
      if (MARK(NbElement) == NOREFRULE)
        SETMARK(NbElement,COPY_REFRULE);
      SETMARKCLASS(NbElement,GREEN);
    }
  }

  /* set additional pattern on the edges */
  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (MARKCLASS(theElement)!=RED) continue;
    for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
      for (j=i+1; j<MAX_CORNERS_OF_ELEM; j++)
      {
        MyEdge=GetEdge(CORNER(theElement,i),CORNER(theElement,j));
        assert (MyEdge != NULL);
        SETADDPATTERN(MyEdge,0);
      }
  }

  return(cnt);
}


/****************************************************************************/
/*																			*/
/* Function:  PrintMarks													*/
/*																			*/
/* Purpose:   print marks of all the elements of that grid					*/
/*																			*/
/* Param:	  GRID *theGrid: grid to scan									*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void PrintMarks (GRID *theGrid)
{
  ELEMENT *theElement;
  char buffer[64];

  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    sprintf(buffer,"MARK (ID[ELEM] = %ld) = %ld \n",ID(theElement),MARK(theElement));
    UserWrite(buffer);
  }

}

/****************************************************************************/
/*																			*/
/* Function:  DropMarks                                                                                                         */
/*																			*/
/* Purpose:   drop marks from leafelements to first regular, and reset		*/
/*			  marks on all element above (important for restrict marks)     */
/*																			*/
/* Param:	  MULTIGRID *theMG												*/
/*																			*/
/* return:	  INT 0: ok                                                                                                     */
/*			  INT 1: error													*/
/*																			*/
/****************************************************************************/

static INT DropMarks (MULTIGRID *theMG)
{
  INT k, Mark;
  GRID *theGrid;
  ELEMENT *theElement, *FatherElement;

  for (k=theMG->topLevel; k>0; k--)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
      if (LEAFELEM(theElement) && (MARKCLASS(theElement) == RED) && (ECLASS(theElement) != RED))
      {
        Mark = MARK(theElement);
        FatherElement = theElement;
        while(ECLASS(FatherElement) != RED)
        {
          SETMARK(FatherElement,NOREFRULE);
          FatherElement = EFATHER(FatherElement);
        }
        SETMARK(FatherElement,Mark);
        SETMARKCLASS(FatherElement,RED);
      }
  }
  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  ComputEdgeTypes											*/
/*																			*/
/* Purpose:   refine whole multigrid structure					*/
/*																			*/
/* Param:	  MULTIGRID *theMG: multigrid to refine                                 */
/*																			*/
/* return:	  INT 0: ok                                                                     */
/*			  INT 1: error								*/
/*																			*/
/****************************************************************************/

static INT ResetEdgeNew (GRID *theGrid)
{
  NODE *theNode;
  LINK *theLink;
  EDGE *theEdge;

  if (theGrid==NULL) return(0);
  for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
  {
    for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
      if (ID(NBNODE(theLink))>ID(theNode))
      {
        theEdge = MYEDGE(theLink);
        assert(theEdge!=NULL);
        SETEDGENEW(theEdge,0);
      }
  }

  return (0);

}

static INT ComputEdgeTypes (GRID *theGrid)
{
  ELEMENT *theElement;
  ELEMENTCONTEXT theContext;
  EDGE *theEdge;
  REFRULE *theRule;
  EDGEDATA *eData;
  INT i, j;

  for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (REFINE(theElement)==NOREFRULE) continue;
    theRule = &(Rules[REFINE(theElement)]);
    GetCurrentContext (theElement,theContext);
    for (eData=&(theRule->edges[0]); eData->type!=0; eData++)
    {
      theEdge = GetEdge(theContext[eData->from],theContext[eData->to]);
      assert(theEdge!=NULL);
      SETTAG(theEdge,eData->type);
    }
  }

  return (0);
}


/****************************************************************************/
/*
   RefineMultiGrid - Refine whole multigrid structure

   SYNOPSIS:
   INT RefineMultiGrid (MULTIGRID *theMG, INT flag);

   PARAMETERS:
   .  theMG - multigrid to refine
   .  flag -

   DESCRIPTION:
   This function refines whole multigrid structure.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if out of memory, but data structure as before
   .n   2 if fatal memory error, data structure corrupted.
 */
/****************************************************************************/

INT RefineMultiGrid (MULTIGRID *theMG, INT flag)
{
  int j,k,r;
  int newlevel;
  NODE *theNode;
  GRID *theGrid, *FinerGrid;
  ELEMENT *theElement, *FatherElement;

  theFullRefRule = YAlignment;
  /*theFullRefRule = ShortestInteriorEdge;*/
  rFlag=flag;           /* set global variable */

  /* drop marks to regular elements */
  if (DropMarks(theMG))
    return(GM_ERROR);

  /* prepare algebra (set internal flags correctly */
  PrepareAlgebraModification(theMG);

  /* compute modification of coarser levels from above */
  j = theMG->topLevel;
  for (k=j; k>0; k--)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    CloseGrid(theMG->grids[k]);
    RestrictMarks(theMG->grids[k-1]);
  }

  newlevel = 0;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    if (k<j) FinerGrid = GRID_ON_LEVEL(theMG,k+1);else FinerGrid = NULL;

    /* reset some flags */
    SETMODIFIED(theGrid,0);
    for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode)) SETMODIFIED(theNode,0);

    /* leave only regular marks */
    for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    {
      if ((ECLASS(theElement)==RED) && MARKCLASS(theElement)==RED) continue;
      SETMARK(theElement,NOREFRULE);
    }

    /* determine regular and irregular elements on next level */
    if ((r = CloseGrid(theGrid))<0)
    {
      PrintErrorMessage('E',"RefineMultiGrid","error in CloseGrid");
      return(GM_FATAL);
    }

    ComputeCopies(theGrid);

    /* dispose connections that may be changed on next level, this is determined */
    /* by the neighborhood of elements were MARK != REFINE.                                      */
    /* This will leave some flags where to rebuild connections later			 */
    if (k<j)
    {
      for (theElement=FIRSTELEMENT(FinerGrid); theElement!=NULL; theElement=SUCCE(theElement))
      {
        assert(EFATHER(theElement) != NULL);
        if (REFINE(EFATHER(theElement))!=MARK(EFATHER(theElement)))
          if (DisposeConnectionsInNeighborhood(FinerGrid,theElement)!=GM_OK)
            return(GM_FATAL);
      }
    }

    /* create a new grid level, if at least one element is refined on finest level */
    if ( (r>0) && (k==j) )
    {
      newlevel = 1;
      if (CreateNewLevel(theMG)==NULL)
        return(GM_FATAL);
      FinerGrid = GRID_ON_LEVEL(theMG,j+1);
    }

    /* now really manipulate the next finer level */
    if ( k<j || newlevel )
      if (RefineGrid(theGrid)!=0)
        return(GM_FATAL);

    if ((k<j)||(newlevel))
    {
      /* now rebuild connections in neighborhood of elements which have EBUILDCON set */
      /* This flag has been set either by GridDisposeConnection or by CreateElement	*/
      if (GridCreateConnection(FinerGrid)) return (GM_FATAL);
      /* and compute the vector classes on the new (or changed) level */

      if (rFlag==GM_COPY_ALL)
        for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
          SeedVectorClasses(theElement);
      else
      {
        ClearVectorClasses(FinerGrid);
        for (theElement=FIRSTELEMENT(FinerGrid); theElement!=NULL; theElement=SUCCE(theElement))
          if (ECLASS(theElement)>=IRREGULAR_CLASS) SeedVectorClasses(theElement);
        PropagateVectorClasses(FinerGrid);
      }
    }

  }

  DisposeTopLevel(theMG);
  theMG->currentLevel = theMG->topLevel;

  /* set grid status of grid 0 */
  RESETGSTATUS(theMG->grids[0],GRID_CHANGED);

  return(GM_OK);
}

static int FReadRule (FILE *stream, REFRULE *theRule)
{
  int i;
  int ns,p0,p1,p2,p3,p4,p5,p6,p7,pat;
  int type,from,to,side;
  int c0,c1,c2,c3,n0,n1,n2,n3,pa;
  int sn0,sn1;

  if (fscanf(stream,"%d  %d %d %d %d %d %d  %d",&ns,&p0,&p1,&p2,&p3,&p4,&p5,&pat)!=8) return (1);

  theRule->nsons = ns;
  theRule->pattern[0] = p0;
  theRule->pattern[1] = p1;
  theRule->pattern[2] = p2;
  theRule->pattern[3] = p3;
  theRule->pattern[4] = p4;
  theRule->pattern[5] = p5;
  theRule->pattern[6] = p6;
  theRule->pattern[7] = p7;
  theRule->pat = pat;

  for (i=0; i<MAXEDGES; i++)
  {
    if (fscanf(stream," %d %d %d %d",&type,&from,&to,&side)!=4) return (1);
    theRule->edges[i].type = type;
    theRule->edges[i].from = from;
    theRule->edges[i].to   = to;
    theRule->edges[i].side = side;
  }

  for (i=0; i<MAX_SONS; i++)
  {
    if (fscanf(stream," %d %d %d %d %d %d %d %d %d",&c0,&c1,&c2,&c3,&n0,&n1,&n2,&n3,&pa)!=9) return (1);
    theRule->sons[i].corners[0] = c0;
    theRule->sons[i].corners[1] = c1;
    theRule->sons[i].corners[2] = c2;
    theRule->sons[i].corners[3] = c3;
    theRule->sons[i].nb[0]          = n0;
    theRule->sons[i].nb[1]          = n1;
    theRule->sons[i].nb[2]          = n2;
    theRule->sons[i].nb[3]          = n3;
    theRule->sons[i].path           = pa;
  }

  for (i=0; i<NEWCORNERS; i++)
  {
    if (fscanf(stream," %d %d",&sn0,&sn1)!=2) return (1);
    theRule->sonandnode[i][0] = sn0;
    theRule->sonandnode[i][1] = sn1;
  }

  return (0);
}

/****************************************************************************/
/*
   InitRefine3d - Initialization of ugrefine

   SYNOPSIS:
   INT InitRefine3d (void);

   PARAMETERS:
   .  void

   DESCRIPTON:
   This function initializes ugrefine.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured
 */
/****************************************************************************/

INT InitRefine3d (void)
{
  int nRules, nPatterns, err, i, P2R;
  FILE *stream;
  FULLREFRULE *newFRR;
  char buffer[256];

  /************************************************************************/
  /*																		*/
  /* read refinement rules from file 'RefRules.data'						*/
  /*																		*/
  /************************************************************************/

  /* open file */
  if (GetDefaultValue(DEFAULTSFILENAME,"refrulefile",buffer)==0)
    stream = fileopen(buffer,"r");
  else
    stream = fileopen("RefRules.data","r");
  if (stream==NULL)
  {
    UserWrite("ERROR: could not open file 'RefRules.data'\n");
    fclose(stream);
    return (__LINE__);
  }

  /* read Rules and nPatterns from file */
  if (fscanf(stream,"%d %d\n",&nRules,&nPatterns)!=2)
  {
    UserWrite("ERROR: failed to read Rules and nPatterns\n");
    fclose(stream);
    return (__LINE__);
  }

  /* get storage for Rules */
  Rules = (REFRULE *) malloc(nRules*sizeof(REFRULE));
  if (Rules==NULL)
  {
    UserWrite("ERROR: no storage for Rules\n");
    fclose(stream);
    return (__LINE__);
  }

  /* get storage for PatternToRefrule */
  PatternToRefrule = (SHORT *) malloc(nPatterns*sizeof(SHORT));
  if (PatternToRefrule==NULL)
  {
    UserWrite("ERROR: no storage for PatternToRefrule\n");
    fclose(stream);
    return (__LINE__);
  }

  /* read Rules */
  for (i=0; i<nRules; i++)
    if (FReadRule(stream,Rules+i)) return (__LINE__);

  /* read PatternToRefrule */
  for (i=0; i<nPatterns; i++)
  {
    if (fscanf(stream,"%d",&P2R)!=1) return (__LINE__);
    PatternToRefrule[i] = P2R;
  }

  fclose(stream);

  UserWrite("RefRules installed\n");

  /************************************************************************/
  /*																		*/
  /* install best full refrules											*/
  /*																		*/
  /************************************************************************/

  /* install the /Menu directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitRefine3d","could not changedir to root");
    return(__LINE__);
  }
  theBFRRDirID = GetNewEnvDirID();
  if (MakeEnvItem("best full refrule",theBFRRDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitRefine3d","could not install '/best full refrule' dir");
    return(__LINE__);
  }
  if (ChangeEnvDir("/best full refrule")==NULL)
    return(__LINE__);

  newFRR = (FULLREFRULE *) MakeEnvItem("shortestie",theBFRRDirID,sizeof(FULLREFRULE));
  if (newFRR==NULL)
    return(__LINE__);
  newFRR->theFullRefRule = ShortestInteriorEdge;

  newFRR = (FULLREFRULE *) MakeEnvItem("minangle",theBFRRDirID,sizeof(FULLREFRULE));
  if (newFRR==NULL)
    return(__LINE__);
  newFRR->theFullRefRule = MinimalSideAngle;

  newFRR = (FULLREFRULE *) MakeEnvItem("bestm",theBFRRDirID,sizeof(FULLREFRULE));
  if (newFRR==NULL)
    return(__LINE__);
  newFRR->theFullRefRule = BestLaplaceMMatrix;

  newFRR = (FULLREFRULE *) MakeEnvItem("maxper",theBFRRDirID,sizeof(FULLREFRULE));
  if (newFRR==NULL)
    return(__LINE__);
  newFRR->theFullRefRule = MaxPerpendicular;

  newFRR = (FULLREFRULE *) MakeEnvItem("mra",theBFRRDirID,sizeof(FULLREFRULE));
  if (newFRR==NULL)
    return(__LINE__);
  newFRR->theFullRefRule = MaxRightAngle;

  newFRR = (FULLREFRULE *) MakeEnvItem("maxarea",theBFRRDirID,sizeof(FULLREFRULE));
  if (newFRR==NULL)
    return(__LINE__);
  newFRR->theFullRefRule = MaxArea;

  newFRR = (FULLREFRULE *) MakeEnvItem("minentry",theBFRRDirID,sizeof(FULLREFRULE));
  if (newFRR==NULL)
    return(__LINE__);
  newFRR->theFullRefRule = MinimalSideEntry;

  /* default full refrule */
  theFullRefRule = ShortestInteriorEdge;

  return (GM_OK);
}
