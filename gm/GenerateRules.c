// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	GenerateRules.c                                                                                                 */
/*																			*/
/* Purpose: generate the rules for ug3d                                                                         */
/*																			*/
/* Author:	Henrik Reichert                                                                                                 */
/*			Institut fuer Angewandte Mathematik                                                     */
/*			Universitaet Heidelberg                                                                                 */
/*			Im Neuenheimer Feld 294                                                                                 */
/*			6900 Heidelberg                                                                                                 */
/*																			*/
/* History: 10.9.1993 begin                                                                                             */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "GenerateRules.h"

#include "gm.h"
#include "simplex.h"
#include "ugrefine3d.h"

#include "general.h"

/* only compile for 3D version */
#ifdef __THREEDIM__

/****************************************************************************/
/*																			*/
/*		defines and macros													*/
/*																			*/
/****************************************************************************/

#define PATTERNFILTER   0x3F    /* mask for the first 6 digits of pattern: 00111111     */
#define NOFATHERSIDE    5               /* has to be just greater than MAX_SIDES_OF_ELEM		*/

#define NRULES          242
#define MAXCORNERS      11                      /* max no. of new corners (incl. center)	*/

#define NOTDONE         -1                      /* SHORT has to be signed!					*/
#define DONE             0
#define TOUCHED          1                      /* has been visited by FindPathForNeighbours*/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/*		global variables													*/
/*																			*/
/****************************************************************************/

static EDGEDATA EdgeData[MAXCORNERS][MAXCORNERS] =
{
  {{0,0,0,0},{0,0,1,0},{0,0,2,0},{0,0,3,0},{0,0,4,0},{2,0,5,0},{0,0,6,0},{0,0,7,0},{2,0,8,3},{2,0,9,2},{1,0,10,0}},
  {{0,0,1,0},{0,1,1,0},{0,1,2,0},{0,1,3,0},{0,1,4,0},{0,1,5,0},{2,1,6,0},{2,1,7,3},{0,1,8,0},{2,1,9,1},{1,1,10,0}},
  {{0,0,2,0},{0,1,2,0},{0,2,2,0},{0,2,3,0},{2,2,4,0},{0,2,5,0},{0,2,6,0},{2,2,7,2},{2,2,8,1},{0,2,9,0},{1,2,10,0}},
  {{0,0,3,0},{0,1,3,0},{0,2,3,0},{0,3,3,0},{2,3,4,3},{2,3,5,1},{2,3,6,2},{0,3,7,0},{0,3,8,0},{0,3,9,0},{1,3,10,0}},
  {{0,0,4,0},{0,1,4,0},{2,2,4,0},{2,3,4,3},{0,4,4,0},{2,4,5,0},{2,4,6,0},{2,4,7,3},{2,4,8,3},{1,4,9,0},{1,4,10,0}},
  {{2,0,5,0},{0,1,5,0},{0,2,5,0},{2,3,5,1},{2,4,5,0},{0,5,5,0},{2,5,6,0},{1,5,7,0},{2,5,8,1},{2,5,9,1},{1,5,10,0}},
  {{0,0,6,0},{2,1,6,0},{0,2,6,0},{2,3,6,2},{2,4,6,0},{2,5,6,0},{0,6,6,0},{2,6,7,2},{1,6,8,0},{2,6,9,2},{1,6,10,0}},
  {{0,0,7,0},{2,1,7,3},{2,2,7,2},{0,3,7,0},{2,4,7,3},{1,5,7,0},{2,6,7,2},{0,7,7,0},{2,7,8,3},{2,7,9,2},{1,7,10,0}},
  {{2,0,8,3},{0,1,8,0},{2,2,8,1},{0,3,8,0},{2,4,8,3},{2,5,8,1},{1,6,8,0},{2,7,8,3},{0,8,8,0},{2,8,9,1},{1,8,10,0}},
  {{2,0,9,2},{2,1,9,1},{0,2,9,0},{0,3,9,0},{1,4,9,0},{2,5,9,1},{2,6,9,2},{2,7,9,2},{2,8,9,1},{0,9,9,0},{1,9,10,0}},
  {{1,0,10,0},{1,1,10,0},{1,2,10,0},{1,3,10,0},{1,4,10,0},{1,5,10,0},{1,6,10,0},{1,7,10,0},{1,8,10,0},{1,9,10,0},{0,10,10,0}}
};

static SHORT CornerOfSonSide[MAX_SIDES_OF_ELEM][3]       = {{0,1,2},{1,2,3},{0,2,3},{0,1,3}};
static SHORT IsOnSide[MAX_SIDES_OF_ELEM][MAXCORNERS] = {{1,1,1,0,1,1,1,0,0,0,0},
                                                        {0,1,1,1,0,1,0,0,1,1,0},
                                                        {1,0,1,1,0,0,1,1,0,1,0},
                                                        {1,1,0,1,1,0,0,1,1,0,0}};

static int nRulesTh;                                    /* total number of theoretical refinement rules */
static int nRules;                                              /* total number of generated refinement rules	*/
static int FirstRuleWithEdgePattern;    /* Rule no. of first of several SidePatterns	*/
static REFRULE *RuleList;                               /* list is allocated in main program			*/
static SHORT *PatternToRefrule;                 /* vector is allocated in main program			*/
static int Output;                                              /* if TRUE output of generated rules			*/

/* prototype for new rules without any data */
static REFRULE EmptyRule =
{0,{0,0,0,0,0,0},0,{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},

 {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
          {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
          {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
          {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},

 {{{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
          {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}}};

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/*		functions															*/
/*																			*/
/****************************************************************************/

static SHORT SideWithCorners (SHORT C0,SHORT C1,SHORT C2)
{
  SHORT Side;

  for (Side=0; Side<MAX_SIDES_OF_ELEM; Side++)
    if (IsOnSide[Side][C0] && IsOnSide[Side][C1] &&IsOnSide[Side][C2])
      break;

  if (Side<MAX_SIDES_OF_ELEM)
    return (Side+FATHER_SIDE_OFFSET);
  else
    return (NOFATHERSIDE);
}

/* definition of the rotations of the tetrahedron indices */
static SHORT Ax0Left [MAXCORNERS] = {0,3,1,2,7,8,4,6,9,5,10};
static SHORT Ax0Right[MAXCORNERS] = {0,2,3,1,6,9,7,4,5,8,10};
static SHORT Ax1Left [MAXCORNERS] = {2,1,3,0,5,8,9,6,4,7,10};
static SHORT Ax1Right[MAXCORNERS] = {3,1,0,2,8,4,7,9,5,6,10};
static SHORT Ax2Left [MAXCORNERS] = {3,0,2,1,7,6,9,8,4,5,10};
static SHORT Ax2Right[MAXCORNERS] = {1,3,2,0,8,9,5,4,7,6,10};
static SHORT Ax3Left [MAXCORNERS] = {1,2,0,3,5,6,4,8,9,7,10};
static SHORT Ax3Right[MAXCORNERS] = {2,0,1,3,6,4,5,9,7,8,10};

static int RotateCorners (REFRULE *theRule, SHORT ReplaceBy[MAXCORNERS])
{
  SONDATA *theSon;
  SHORT i,j;

  for (i=0; i<theRule->nsons; i++)
    for (theSon = theRule->sons+i, j=0; j<CORNERS_OF_ELEM; j++)
      theSon->corners[j] = ReplaceBy[theSon->corners[j]];

  return (0);
}

static int PrintRule (REFRULE *theRule, SHORT nEdges)
{
  SHORT i,j;

  printf("nsons %d; pattern %d%d%d%d%d%d; ",(int)theRule->nsons,  (int)theRule->pattern[5],
         (int)theRule->pattern[4],
         (int)theRule->pattern[3],
         (int)theRule->pattern[2],
         (int)theRule->pattern[1],
         (int)theRule->pattern[0]);

  printf("pat %#04x; (nedges %d)\n",(int)theRule->pat,(int)nEdges);

  printf("edges: ");
  for (i=0; i<8; i++)
    printf("{%d,%d,%d,%d}",(int)theRule->edges[i].type,(int)theRule->edges[i].from,
           (int)theRule->edges[i].to,(int)theRule->edges[i].side);
  printf("\n\t");
  for (; i<MAXEDGES; i++)
    printf("{%d,%d,%d,%d}",(int)theRule->edges[i].type,(int)theRule->edges[i].from,
           (int)theRule->edges[i].to,(int)theRule->edges[i].side);

  printf("\nsons:  ");
  for (i=0; i<6; i++)
    printf("{{%d,%d,%d,%d}{%d,%d,%d,%d}}",(int)theRule->sons[i].corners[0],(int)theRule->sons[i].corners[1],
           (int)theRule->sons[i].corners[2],(int)theRule->sons[i].corners[3],
           (int)theRule->sons[i].nb[0],(int)theRule->sons[i].nb[1],
           (int)theRule->sons[i].nb[2],(int)theRule->sons[i].nb[3]);
  printf("\n\t");
  for (; i<SONS; i++)
    printf("{{%d,%d,%d,%d}{%d,%d,%d,%d}}",(int)theRule->sons[i].corners[0],(int)theRule->sons[i].corners[1],
           (int)theRule->sons[i].corners[2],(int)theRule->sons[i].corners[3],
           (int)theRule->sons[i].nb[0],(int)theRule->sons[i].nb[1],
           (int)theRule->sons[i].nb[2],(int)theRule->sons[i].nb[3]);
  for (i=1; i<theRule->nsons; i++)
  {
    printf("\npath to %d: son0",i);
    for (j=0; j<PATHDEPTH(theRule->sons[i].path); j++)
      printf("-->nb%d",(int) NEXTSIDE(theRule->sons[i].path,j));

    printf(" (%d steps)",(int) PATHDEPTH(theRule->sons[i].path));
  }

  printf("\nsonandnode:");
  for (i=0; i<NEWCORNERS; i++)
    printf(" %d->(%d,%d)",(int)(i+CORNERS_OF_ELEM),(int)theRule->sonandnode[i][0],(int)theRule->sonandnode[i][1]);

  printf("\n");

  return (0);
}

static int CopySonsIntoRule (REFRULE *theRule,SONDATA sons[SONS])
{
  SHORT i;

  for (i=0; i<SONS; i++)
    theRule->sons[i] = sons[i];

  return (0);
}

static int CompareEdgesForSides (const void *Edge1, const void *Edge2)
{
  EDGEDATA *theEdge1,*theEdge2;

  theEdge1 = (EDGEDATA *) Edge1;
  theEdge2 = (EDGEDATA *) Edge2;

  /* leave type==0 at end of the list */
  if ((theEdge1->type)==0)
    return (1);
  if ((theEdge2->type)==0)
    return (-1);

  if ((theEdge1->side)==(theEdge2->side))
    return (0);
  else if ((theEdge1->side)>(theEdge2->side))
    return (1);
  else
    return (-1);
}

static int CompareEdges (const void *theEdge1, const void *theEdge2)
{
  if (((EDGEDATA *)theEdge1)->from==((EDGEDATA *)theEdge2)->from)
  {
    if (((EDGEDATA *)theEdge1)->to>((EDGEDATA *)theEdge2)->to)
      return (1);
    else if (((EDGEDATA *)theEdge1)->to<((EDGEDATA *)theEdge2)->to)
      return (-1);
    else
      return (0);
  }
  else if (((EDGEDATA *)theEdge1)->from>((EDGEDATA *)theEdge2)->from)
    return (1);
  else
    return (-1);
}

static int CompareCorners (const void *theCorner1, const void *theCorner2)
{
  SHORT C1,C2;

  C1 = * ((SHORT *) theCorner1);
  C2 = * ((SHORT *) theCorner2);

  assert(C1!=C2);

  if (C1>C2)
    return (1);
  else
    return (-1);
}

static int NewEdgesAreEqual (EDGEDATA *NewEdges0, EDGEDATA *NewEdges1)
{
  SHORT i;

  for (i=0; i<MAXEDGES; i++)
    if ((NewEdges0[i].from!=NewEdges1[i].from)||(NewEdges0[i].to!=NewEdges1[i].to))
      break;

  if (i<MAXEDGES)
    return (0);                 /* EDGEDATA are not equal */
  else
    return (1);                 /* EDGEDATA are equal */
}

static int FindPathForNeighbours (REFRULE *theRule, SHORT myID, SHORT Status[SONS])
{
  SHORT i,nbID;
  SHORT nbPATHDEPTH;
  INT *nbPath;

  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    if (((nbID=theRule->sons[myID].nb[i])<FATHER_SIDE_OFFSET) && (Status[nbID]==NOTDONE))
    {
      nbPath = &(theRule->sons[nbID].path);

      /* copy myPath to nbPath */
      *nbPath = theRule->sons[myID].path;

      /* complete nbPath */
      nbPATHDEPTH = PATHDEPTH(*nbPath);
      SETNEXTSIDE(*nbPath,nbPATHDEPTH,i);
      SETPATHDEPTH(*nbPath,++nbPATHDEPTH);
      Status[nbID] = TOUCHED;
    }

  /* recursive call for TOUCHED sons */
  for (nbID=1; nbID<theRule->nsons; nbID++)
    if (Status[nbID]==TOUCHED)
    {
      Status[nbID] = DONE;
      FindPathForNeighbours (theRule,nbID,Status);
    }

  return (0);
}

/****************************************************************************/
/*																			*/
/*		AddRuleIfNew														*/
/*			called by most of the main functions							*/
/*																			*/
/*			the additional data fields of the refrule are filled and the	*/
/*			rule is stored if not already there                                                     */
/*																			*/
/****************************************************************************/

static int AddRuleIfNew (REFRULE *theRule)
{
  EDGEDATA theEdgeData,*InnerEdge[2];
  SONDATA *SONi,*SONj,*theSon;
  SHORT i,j,nEdges,theEdge,FatherSide,newcorner,corner,son,edge,found,SidePattern;
  SHORT SIDEi,SIDEj,side,SideBit,Status[SONS];

  assert(nRules<NRULES);

  /* search new inner edges (one or two corners > 3) and fill edges field */
  nEdges = 0;
  for (son=0; son<theRule->nsons; son++)
    for (i=0; i<CORNERS_OF_ELEM; i++)
      for (j=i+1; j<CORNERS_OF_ELEM; j++)
        if ((theEdgeData=EdgeData       [theRule->sons[son].corners[i]]
                          [theRule->sons[son].corners[j]]).type)
        {
          /* (From,To) is a new edge: already stored? */
          for (theEdge=0; theEdge<nEdges; theEdge++)
            if ((theRule->edges[theEdge].from==theEdgeData.from)&&(theRule->edges[theEdge].to==theEdgeData.to))
              break;
          if (theEdge>=nEdges)
            theRule->edges[nEdges++] = theEdgeData;
        }

  /* sort edges for corner IDs (necessary for fast comparison of SidePatterns) */
  qsort((void *)theRule->edges,nEdges,sizeof(EDGEDATA),CompareEdges);

  /* if this SidePattern is not already in RuleList: store it */
  for (i=FirstRuleWithEdgePattern; i<nRules; i++)
    if (NewEdgesAreEqual(theRule->edges,RuleList[i].edges))
      break;

  if (i<nRules)
    return (1);                         /* rule is already in RuleList */

  /* sort corners of sons to make searching for neighbour sides easier */
  for (i=0; i<theRule->nsons; i++)
    qsort((void *)theRule->sons[i].corners,4,sizeof(SHORT),CompareCorners);

  /* fill entries of nb field in the SONDATA */
  for (i=0; i<theRule->nsons; i++)
    for (SONi=theRule->sons+i, SIDEi=0; SIDEi<MAX_SIDES_OF_ELEM; SIDEi++)
      if (SONi->nb[SIDEi] == NOTDONE)
      {
        /* outer face? */
        if ((FatherSide=SideWithCorners(SONi->corners[CornerOfSonSide[SIDEi][0]],
                                        SONi->corners[CornerOfSonSide[SIDEi][1]],
                                        SONi->corners[CornerOfSonSide[SIDEi][2]]))!=NOFATHERSIDE)
        {
          SONi->nb[SIDEi] = FatherSide;
          continue;
        }

        for (j=i+1; j<theRule->nsons; j++)
          for (SONj=theRule->sons+j, SIDEj=0; SIDEj<MAX_SIDES_OF_ELEM; SIDEj++)
            if (   (SONi->corners[CornerOfSonSide[SIDEi][0]]==SONj->corners[CornerOfSonSide[SIDEj][0]])
                   && (SONi->corners[CornerOfSonSide[SIDEi][1]]==SONj->corners[CornerOfSonSide[SIDEj][1]])
                   && (SONi->corners[CornerOfSonSide[SIDEi][2]]==SONj->corners[CornerOfSonSide[SIDEj][2]]))
            {
              /* neighbour found */
              SONi->nb[SIDEi] = j;
              SONj->nb[SIDEj] = i;

              j=SONS;
              break;
            }
      }

  /* clear the remaining NOTDONEs */
  for (i=theRule->nsons; i<SONS; i++)
    for (SONi=theRule->sons+i, SIDEi=0; SIDEi<MAX_SIDES_OF_ELEM; SIDEi++)
      SONi->nb[SIDEi] = 0;

  /* fill the sonandnode field */
  for (newcorner=0; newcorner<NEWCORNERS; newcorner++)
    theRule->sonandnode[newcorner][0] = NOTUSED;

  for (newcorner=CORNERS_OF_ELEM; newcorner<(NEWCORNERS+CORNERS_OF_ELEM); newcorner++)
    for (son=0; son<SONS; son++)
      for (theSon=theRule->sons+son, corner=0; corner<CORNERS_OF_ELEM; corner++)
        if (theSon->corners[corner]==newcorner)
        {
          theRule->sonandnode[newcorner-CORNERS_OF_ELEM][0] = son;
          theRule->sonandnode[newcorner-CORNERS_OF_ELEM][1] = corner;

          son = SONS;
          break;
        }

  /* make a consistence check: newcorners in sondata indicate the refined edges */
  for (newcorner=0; newcorner<EDGES; newcorner++)
    if (theRule->sonandnode[newcorner][0]==NOTUSED)
    {
      assert(theRule->pattern[newcorner]==0);
      assert((theRule->pat & (1<<newcorner))==0);
    }
    else
    {
      assert(theRule->pattern[newcorner]==1);
      assert((theRule->pat & (1<<newcorner))>0);
    }

  /* find the SidePattern of theRule */
  SidePattern = 0;
  for (side=0; side<MAX_SIDES_OF_ELEM; side++)
    if (TriSectionEdge[theRule->pat & CondensedEdgeOfSide[side]][0] >= 0)
    {
      /* here we have a side with a trisection edge */
      for (found=0, edge=0; edge<MAXEDGES; edge++)
        if ((theRule->edges[edge].type==2) && (theRule->edges[edge].side==side))
        {
          /* we have an inner edge of side (which has a trisection node) */
          InnerEdge[found++] = theRule->edges+edge;
        }
      assert(found==2);

      if (InnerEdge[0]->from==InnerEdge[1]->from)
        theEdge = InnerEdge[0]->from - CORNERS_OF_ELEM;
      else if (InnerEdge[0]->from==InnerEdge[1]->to)
        theEdge = InnerEdge[0]->from - CORNERS_OF_ELEM;
      else if (InnerEdge[0]->to==InnerEdge[1]->from)
        theEdge = InnerEdge[0]->to - CORNERS_OF_ELEM;
      else if (InnerEdge[0]->to==InnerEdge[1]->to)
        theEdge = InnerEdge[0]->to - CORNERS_OF_ELEM;
      else
        assert(0);

      if (TriSectionEdge[theRule->pat & CondensedEdgeOfSide[side]][0] == theEdge)
        SideBit = 0;
      else if (TriSectionEdge[theRule->pat & CondensedEdgeOfSide[side]][1] == theEdge)
        SideBit = 1;
      else
        assert(0);

      if (SideBit)
        SidePattern |= (1<<side);

      if (Output)
        printf("trisection side is %d, edge %d\n",(int)side,(int)theEdge);
    }

  /* find path from son_1 to son_i */
  Status[0] = DONE;
  for (i=1; i<theRule->nsons; i++)
    Status[i] = NOTDONE;

  FindPathForNeighbours (theRule,0,Status);

  /* put the new rule into the RuleList */
  PatternToRefrule[theRule->pat | (SidePattern<<EDGES)] = nRules;
  RuleList[nRules++] = *theRule;

  /* print data for check purposes */
  if (Output)
  {
    printf("### rule no %3d ###\n",nRules-1);
    PrintRule(RuleList+(nRules-1),nEdges);

    printf("SidePattern is %#04x\n",(int)SidePattern);
    printf("### end of rule ###\n");
  }

  return (0);
}

static int BisectSon (SHORT theSon, SHORT Corner0ID, SHORT Corner1ID, SHORT theEdge, SONDATA *theSD, SHORT nsons)
{
  SHORT i;

  /* copy old son into new son */
  for (i=0; i<CORNERS_OF_ELEM; i++)
    theSD[nsons].corners[i] = theSD[theSon].corners[i];

  theSD[theSon].corners[Corner0ID] = MidNodeOfEdge[theEdge];
  theSD[nsons].corners[Corner1ID]  = MidNodeOfEdge[theEdge];

  return (0);
}

static int BisectEdge (SHORT EdgePattern[EDGES], REFRULE theRule)
{
  REFRULE theNewRule;
  SHORT SonID[2],CornerID0[2],CornerID1[2];
  SHORT theEdge,nEdges,found,Corner0,Corner1,theSon,i,j;

  for (nEdges=0,theEdge=0; theEdge<EDGES; theEdge++)
    if (EdgePattern[theEdge])
    {
      ++nEdges;

      /* refine son(s) containing theEdge by bisection */
      if (Output)
        printf("\trefining edge %d\n",(int)theEdge);

      /* first we need a copy of the original rule */
      theNewRule = theRule;

      found = 0;
      Corner0 = CornerOfEdge[theEdge][0];
      Corner1 = CornerOfEdge[theEdge][1];
      for (theSon=0; theSon<theNewRule.nsons; theSon++)
        for (i=0; i<CORNERS_OF_ELEM; i++)
          if (theNewRule.sons[theSon].corners[i]==Corner0)
            for (j=0; j<CORNERS_OF_ELEM; j++)
              if (theNewRule.sons[theSon].corners[j]==Corner1)
              {
                SonID[found] = theSon;
                CornerID0[found] = i;
                CornerID1[found] = j;
                ++found;
              }

      assert(found==1 || found==2);

      for (--found; found>=0; found--)
        BisectSon(SonID[found],CornerID0[found],CornerID1[found],theEdge,theNewRule.sons,theNewRule.nsons++);

      /* recursively refine the remaining edges: */
      EdgePattern[theEdge] = 0;
      BisectEdge(EdgePattern,theNewRule);
      EdgePattern[theEdge] = 1;
    }

  if (nEdges==0)
    AddRuleIfNew(&theRule);

  return (0);
}


/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		main functions:                                                                                                         */
/*			MakeRule...                                                                                                     */
/*																			*/
/*		the Rules refining n edges are generated							*/
/*																			*/
/****************************************************************************/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/*		MakeRuleByBisection                                                                                             */
/*			generates rules by recursive bisection of the father			*/
/*																			*/
/****************************************************************************/

static int MakeRuleByBisection (SHORT pattern)
{
  REFRULE NewRule;
  SONDATA father = {{0,1,2,3},{NOTDONE,NOTDONE,NOTDONE,NOTDONE}};
  SHORT EdgePattern[EDGES],RefEdgesPerSide[MAX_SIDES_OF_ELEM];
  SHORT i,nTrapezoid;

  NewRule = EmptyRule;
  NewRule.nsons = 1;
  NewRule.sons[0] = father;

  NewRule.pat = pattern;

  /* make extended pattern */
  for (i=0; i<EDGES; i++)
    NewRule.pattern[i] = EdgePattern[i] = (((1<<i)&pattern)>0);

  /* find the number of trapezoids for check purposes */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    RefEdgesPerSide[i] = 0;

  for (i=0; i<EDGES; i++)
    if (EdgePattern[i])
    {
      ++RefEdgesPerSide[(SideWithEdge[i][0])];
      ++RefEdgesPerSide[(SideWithEdge[i][1])];
    }

  for (nTrapezoid=0,i=0; i<MAX_SIDES_OF_ELEM; i++)
    if (RefEdgesPerSide[i]==2)
      ++nTrapezoid;

  printf("   make rule by bisection: %d Rules\n",(int)(1<<nTrapezoid));

  /* refine tetrahedron with this pattern */
  BisectEdge(EdgePattern,NewRule);

  if (Output)
    printf("   %d trapezoids found ==> %d SidePatterns\n",(int)nTrapezoid,(int)(1<<nTrapezoid));

  nRulesTh += 1<<nTrapezoid;

  return (0);
}

static int MakeRule0 (SHORT pattern)
{
  REFRULE CopyRule =
  {1,{0,0,0,0,0,0},0,{{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0}},

   {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
                  {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
                  {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},

   {{{0,1,2,3},{FATHER_SIDE_OFFSET+0,FATHER_SIDE_OFFSET+1,FATHER_SIDE_OFFSET+2,FATHER_SIDE_OFFSET+3},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0}}};
  REFRULE NoRefRule =
  {0,{0,0,0,0,0,0},0,{{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0},{NOTUSED,0}},

   {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
                  {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
                  {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},

   {{{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0},
                  {{0,0,0,0},{0,0,0,0},0}}};

  assert(pattern==0);           /* just to supress 'not used' warning */

  printf("   rule 0: no ref and copy: 2 Rules\n");
  nRulesTh += 2;

  PatternToRefrule[NoRefRule.pat] = nRules;
  RuleList[nRules++] = NoRefRule;
  RuleList[nRules++] = CopyRule;

  return (0);
}


#define IDENTITY                        NULL

#define RULE3_CORNER    1
#define RULE3CORNER_PROTO       2

#define RULE3_SIDE              2

static int MakeRule3 (SHORT pattern)
{
  SONDATA CornerPrototype[RULE3CORNER_PROTO][SONS] =
  {
    {{{0,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{4,6,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,2,3,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,3,4,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,4,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,2,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,4,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,6,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,3,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},

    {{{0,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{4,6,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,2,3,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,3,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,4,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,2,4,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,4,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,3,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,6,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}}
  };
  SONDATA SidePrototype[SONS] =
  {{{0,4,2,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{1,4,2,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{4,2,7,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{7,8,2,3},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}};

  REFRULE theRule;
  SHORT *Rot;
  SHORT theRule3,i,j;

  /* see wether Rule3 is to apply */
  theRule3 = 0;
  switch (pattern)
  {
  /* the edges to be refined share one node */

  case 0x0D :                           /* 3 edges (0,2,3) at node 0: 001101 */
    Rot = IDENTITY; theRule3 = RULE3_CORNER; break;

  case 0x13 :                           /* 3 edges (0,1,4) at node 1: 010011 */
    Rot = Ax3Left; theRule3 = RULE3_CORNER; break;

  case 0x26 :                           /* 3 edges (1,2,5) at node 2: 100110 */
    Rot = Ax3Right; theRule3 = RULE3_CORNER; break;

  case 0x38 :                           /* 3 edges (3,4,5) at node 3: 111000 */
    Rot = Ax1Right; theRule3 = RULE3_CORNER; break;


  /* the edges to be refined lie on one side */

  case 0x07 :                           /* 3 edges (0,1,2) at side 0: 000111 */
    Rot = Ax1Left; theRule3 = RULE3_SIDE; break;

  case 0x32 :                           /* 3 edges (1,4,5) at side 1: 110010 */
    Rot = Ax3Left; theRule3 = RULE3_SIDE; break;

  case 0x2C :                           /* 3 edges (2,3,5) at side 2: 101100 */
    Rot = Ax3Right; theRule3 = RULE3_SIDE; break;

  case 0x19 :                           /* 3 edges (0,3,4) at side 3: 010001 */
    Rot = IDENTITY; theRule3 = RULE3_SIDE; break;
  }

  if (!theRule3)
    return (1);

  /* apply theRule3 */

  /* rotate the Prototype son patterns into the pattern position and copy them into the rule */
  if (theRule3==RULE3_CORNER)
  {
    printf("   rule 3: 3 next to one node: 2 Rules (extra center-node)\n");
    /* nRulesTh += 2;  will be counted by MakeRuleByBisection */

    for (i=0; i<RULE3CORNER_PROTO; i++)
    {
      theRule = EmptyRule;
      theRule.nsons = 9;
      CopySonsIntoRule(&theRule,CornerPrototype[i]);

      if (Rot!=IDENTITY)
        RotateCorners(&theRule,Rot);

      theRule.pat = pattern;

      /* fill extended pattern */
      for (j=0; j<EDGES; j++)
        theRule.pattern[j] = (((1<<j)&pattern)>0);

      if (AddRuleIfNew(&theRule)==1)
        printf("MakeRule3: already there\n");
    }
    return (1);
  }
  else
  {
    printf("   rule 3: 3 on one side: 1 Rule\n");
    nRulesTh += 1;

    theRule = EmptyRule;
    theRule.nsons = 4;

    CopySonsIntoRule(&theRule,SidePrototype);

    if (Rot!=IDENTITY)
      RotateCorners(&theRule,Rot);

    theRule.pat = pattern;

    /* fill extended pattern */
    for (j=0; j<EDGES; j++)
      theRule.pattern[j] = (((1<<j)&pattern)>0);

    if (AddRuleIfNew(&theRule)==1)
      printf("MakeRule3: already there\n");

    return (0);
  }
}

#define RULE4CORNER             1
#define RULE4CORNER_PROTO       4

#define RULE4PLANE                      2
#define RULE4PLANE_PROTO        2

static int MakeRule4 (SHORT pattern)
{
  SONDATA CornerPrototype[RULE4CORNER_PROTO][SONS] =
  {
    {{{0,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{5,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,4,7,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,4,1,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{6,3,7,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{6,3,2,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},
    {{{0,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{5,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,4,7,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,4,1,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,7,3,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,7,6,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},
    {{{0,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{5,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,7,4,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,7,3,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,7,3,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,7,6,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},
    {{{0,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{5,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,7,4,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,7,3,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{6,3,7,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{6,3,2,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}}
  };
  SONDATA PlanePrototype[RULE4PLANE_PROTO][SONS] =
  {
    {{{1,8,5,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,8,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,6,5,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,6,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,5,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,5,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,5,8,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,7,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,7,8,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,2,8,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,2,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},
    {{{1,5,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,5,8,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,7,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,7,8,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,6,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,8,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,6,5,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,6,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,8,5,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{3,8,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,3,5,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,3,7,10},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}}
  };
  SHORT *Rot0,*Rot1;
  REFRULE theRule;
  SHORT i,j,theRule4;


  if ((pattern&0x0D)==0x0D)
    switch (pattern&~0x0D)                              /* corner 0 */
    {
    case (1<<1) : Rot0 = IDENTITY; Rot1 = IDENTITY; theRule4 = RULE4CORNER; break;
    case (1<<4) : Rot0 = IDENTITY; Rot1 = Ax0Left;  theRule4 = RULE4CORNER; break;
    case (1<<5) : Rot0 = IDENTITY; Rot1 = Ax0Right; theRule4 = RULE4CORNER; break;
    }
  else if ((pattern&0x13)==0x13)
    switch (pattern&~0x13)                              /* corner 1 */
    {
    case (1<<2) : Rot0 = Ax3Left; Rot1 = IDENTITY; theRule4 = RULE4CORNER; break;
    case (1<<3) : Rot0 = Ax3Left; Rot1 = Ax1Right; theRule4 = RULE4CORNER; break;
    case (1<<5) : Rot0 = Ax3Left; Rot1 = Ax1Left;  theRule4 = RULE4CORNER; break;
    }
  else if ((pattern&0x26)==0x26)
    switch (pattern&~0x26)                              /* corner 2 */
    {
    case (1<<0) : Rot0 = Ax3Right; Rot1 = IDENTITY; theRule4 = RULE4CORNER; break;
    case (1<<3) : Rot0 = Ax3Right; Rot1 = Ax2Left;  theRule4 = RULE4CORNER; break;
    case (1<<4) : Rot0 = Ax3Right; Rot1 = Ax2Right; theRule4 = RULE4CORNER; break;
    }
  else if ((pattern&0x38)==0x38)
    switch (pattern&~0x38)                              /* corner 3 */
    {
    case (1<<0) : Rot0 = Ax1Right; Rot1 = IDENTITY; theRule4 = RULE4CORNER; break;
    case (1<<1) : Rot0 = Ax1Right; Rot1 = Ax3Left;  theRule4 = RULE4CORNER; break;
    case (1<<2) : Rot0 = Ax1Right; Rot1 = Ax3Right; theRule4 = RULE4CORNER; break;
    }
  else if (pattern==0x1E)                               /* plane 5,6,7,8 (011110) */
  {
    Rot0 = IDENTITY; theRule4 = RULE4PLANE;
  }
  else if (pattern==0x35)                               /* plane 4,6,8,9 (110101) */
  {
    Rot0 = Ax3Left; theRule4 = RULE4PLANE;
  }
  else if (pattern==0x2B)                               /* plane 4,5,7,9 (101011) */
  {
    Rot0 = Ax3Right; theRule4 = RULE4PLANE;
  }
  else
    assert(0);                                                          /* should not happen */

  /* rotate the Prototype son patterns into the pattern position and copy them into the rule */
  if (theRule4==RULE4CORNER)
  {
    printf("   rule 4: 3 next to one node: 4 Rules\n");
    nRulesTh += 4;

    for (i=0; i<RULE4CORNER_PROTO; i++)
    {
      theRule = EmptyRule;
      theRule.nsons = 6;
      CopySonsIntoRule(&theRule,CornerPrototype[i]);

      if (Rot0!=IDENTITY)
        RotateCorners(&theRule,Rot0);
      if (Rot1!=IDENTITY)
        RotateCorners(&theRule,Rot1);

      theRule.pat = pattern;

      /* fill extended pattern */
      for (j=0; j<EDGES; j++)
        theRule.pattern[j] = (((1<<j)&pattern)>0);

      if (AddRuleIfNew(&theRule)==1)
        printf("MakeRule4: already there\n");
    }
    return (0);
  }
  else
  {
    printf("   rule 4: 4 on a plane: 2 Rules\n");
    /* nRulesTh += 2;  will be counted by MakeRuleByBisection */

    for (i=0; i<RULE4PLANE_PROTO; i++)
    {
      theRule = EmptyRule;
      theRule.nsons = 12;
      CopySonsIntoRule(&theRule,PlanePrototype[i]);

      if (Rot0!=IDENTITY)
        RotateCorners(&theRule,Rot0);

      theRule.pat = pattern;

      /* fill extended pattern */
      for (j=0; j<EDGES; j++)
        theRule.pattern[j] = (((1<<j)&pattern)>0);

      if (AddRuleIfNew(&theRule)==1)
        printf("MakeRule4: already there\n");
    }
    return (1);
  }
}

#define RULE5_PROTO     4

static int MakeRule5 (SHORT pattern)
{
  SONDATA CornerPrototype[RULE5_PROTO][SONS] =
  {
    {{{3,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,9,6,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{6,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{5,6,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,6,7,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,6,5,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,5,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},
    {{{3,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,9,6,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{6,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{5,6,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,6,8,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,6,8,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},
    {{{3,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,9,6,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{6,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{5,6,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,6,8,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,6,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,7,6,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}},
    {{{3,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{2,9,6,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{9,7,5,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{9,7,5,6},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,7,6,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,1,7,5},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{1,7,5,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
     {{0,0,0,0},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0}}
  };
  SHORT *Rot0,*Rot1;
  REFRULE theRule;
  SHORT i,j;

  printf("   rule 5: 4 Rules\n");
  nRulesTh += 4;

  switch (~pattern&PATTERNFILTER)
  {
  case (1<<0) : Rot0 = IDENTITY; Rot1 = IDENTITY; break;                /* all but egde 0 (111110) */
  case (1<<1) : Rot0 = Ax3Left;  Rot1 = IDENTITY; break;                /* all but egde 1 (111101) */
  case (1<<2) : Rot0 = Ax3Right; Rot1 = IDENTITY; break;                /* all but egde 2 (111011) */
  case (1<<3) : Rot0 = Ax2Left;  Rot1 = IDENTITY; break;                /* all but egde 3 (110111) */
  case (1<<4) : Rot0 = Ax2Right; Rot1 = IDENTITY; break;                /* all but egde 4 (101111) */
  case (1<<5) : Rot0 = Ax2Left;  Rot1 = Ax1Right; break;                /* all but egde 5 (011111) */

  default : assert (0);                                                                                 /* should not happen	   */
  }

  /* rotate the Prototype son patterns into the pattern position and copy them into the rule */
  for (i=0; i<RULE5_PROTO; i++)
  {
    theRule = EmptyRule;
    theRule.nsons = 7;
    CopySonsIntoRule(&theRule,CornerPrototype[i]);

    if (Rot0!=IDENTITY)
      RotateCorners(&theRule,Rot0);

    if (Rot1!=IDENTITY)
      RotateCorners(&theRule,Rot1);

    theRule.pat = pattern;

    /* fill extended pattern */
    for (j=0; j<EDGES; j++)
      theRule.pattern[j] = (((1<<j)&pattern)>0);

    if (AddRuleIfNew(&theRule)==1)
      printf("MakeRule5: already there\n");
  }

  return (0);
}

static int MakeRule6 (SHORT pattern)
{
  SONDATA FullRefPrototype[SONS] =
  {{{0,4,6,7},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{4,1,5,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{6,5,2,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{7,8,9,3},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{4,6,7,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{4,6,5,8},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{6,5,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{6,7,8,9},{NOTDONE,NOTDONE,NOTDONE,NOTDONE},0},
   {{0,0,0,0},{0,0,0,0},0},
   {{0,0,0,0},{0,0,0,0},0},
   {{0,0,0,0},{0,0,0,0},0},
   {{0,0,0,0},{0,0,0,0},0}};

  REFRULE theRule;
  SHORT j;

  printf("   rule 6: 3 Rules\n");
  nRulesTh += 3;

  /* inner edge 6-8 */
  theRule = EmptyRule;
  theRule.nsons = 8;
  theRule.pat = pattern;
  for (j=0; j<EDGES; j++)
    theRule.pattern[j] = (((1<<j)&pattern)>0);
  CopySonsIntoRule(&theRule,FullRefPrototype);
  AddRuleIfNew(&theRule);

  /* inner edge 4-9 */
  theRule = EmptyRule;
  theRule.nsons = 8;
  theRule.pat = pattern;
  for (j=0; j<EDGES; j++)
    theRule.pattern[j] = (((1<<j)&pattern)>0);
  CopySonsIntoRule(&theRule,FullRefPrototype);
  RotateCorners(&theRule,Ax1Left);
  FirstRuleWithEdgePattern = nRules;
  AddRuleIfNew(&theRule);

  /* inner edge 5-7 */
  theRule = EmptyRule;
  theRule.nsons = 8;
  theRule.pat = pattern;
  for (j=0; j<EDGES; j++)
    theRule.pattern[j] = (((1<<j)&pattern)>0);
  CopySonsIntoRule(&theRule,FullRefPrototype);
  RotateCorners(&theRule,Ax1Right);
  FirstRuleWithEdgePattern = nRules;
  AddRuleIfNew(&theRule);

  /* reset PatternToRefrule[111111] to the first of the FullRefRules */
  PatternToRefrule[pattern] = nRules-3;

  return (0);
}

int main (int argc, char **argv)
{
#ifndef TEST_GR
        #pragma unused (argc,argv)
#endif

  FILE *stream;
  SHORT pattern,nRefEdge,nRE,i;
  int MINnRefEdge,MAXnRefEdge,j,err,SaveList;

  MINnRefEdge = 0;
  MAXnRefEdge = EDGES;

  Output   = 0;
  SaveList = 1;

#ifdef TEST_GR
  if (argc<3)
  {
    printf("usage: GenerateRules <min # of refined edges> <max # of refined edges> [-o] [-s]\n");
    return (1);
  }

  /* scan args for options */
  if (sscanf(argv[1],"%d",&MINnRefEdge)!=1)
    return (1);

  if (sscanf(argv[2],"%d",&MAXnRefEdge)!=1)
    return (1);

  Output = 0;
  for (i=3; i<argc; i++)
    if ((argv[i][0]=='-') && (argv[i][1]=='o'))
    {
      Output = 1;
      break;
    }

  SaveList = 0;
  for (i=3; i<argc; i++)
    if ((argv[i][0]=='-') && (argv[i][1]=='s'))
    {
      SaveList = 1;
      break;
    }

#endif

  if ((MINnRefEdge<0) || (MINnRefEdge>EDGES))
  {
    printf("ERROR: min # of refined edges out of range\n");
    return (1);
  }

  if ((MAXnRefEdge<0) || (MAXnRefEdge>EDGES))
  {
    printf("ERROR: min # of refined edges out of range\n");
    return (1);
  }

  if (MAXnRefEdge<MINnRefEdge)
  {
    printf("ERROR: min # of refined edges > max\n");
    return (1);
  }

  /* get storage for RuleList */
  RuleList = (REFRULE *) malloc(NRULES*sizeof(REFRULE));
  if (RuleList==NULL)
  {
    printf("ERROR: no storage for RuleList\n");
    return (1);
  }

  /* get storage for PatternToRefrule */
  PatternToRefrule = (SHORT *) malloc((1<<(EDGES+MAX_SIDES_OF_ELEM))*sizeof(SHORT));
  if (PatternToRefrule==NULL)
  {
    printf("ERROR: no storage for PatternToRefrule\n");
    return (1);
  }
  for (j=0; j<(1<<(EDGES+MAX_SIDES_OF_ELEM)); j++)
    PatternToRefrule[j] = NOTUSED;

  nRules = nRulesTh = 0;

  for (nRefEdge=MINnRefEdge; nRefEdge<=MAXnRefEdge; nRefEdge++)
  {
    printf("\n### %d refined edges ###\n",(int)nRefEdge);

    for (pattern=0; pattern < (1<<EDGES); pattern++)
    {
      /* count refined edges */
      for (nRE=0, i=0; i<EDGES; i++)
        nRE += ((1<<i)&pattern)>0;

      /* take only pattern with nRE==nRefEdge */
      if (nRE!=nRefEdge) continue;

      /* set begin of possibly identical SidePatterns for AddRuleIfNew () */
      FirstRuleWithEdgePattern = nRules;

      for (i=0; i<EDGES; i++)
        printf("%1d",(((1<<i)&pattern)>0));

      switch (nRefEdge)
      {
      case 0 :
        MakeRule0(pattern);
        break;

      case 1 :
        MakeRuleByBisection(pattern);
        break;

      case 2 :
        MakeRuleByBisection(pattern);
        break;

      case 3 :
        if (MakeRule3(pattern))
          MakeRuleByBisection(pattern);
        break;

      case 4 :
        if (MakeRule4(pattern))
          MakeRuleByBisection(pattern);
        break;

      case 5 :
        MakeRule5(pattern);
        break;

      case 6 :
        MakeRule6(pattern);
        break;
      }
    }
  }

  /* sort edges for sides (necessary for ugrefine) */
  for (i=2; i<nRules; i++)
    qsort((void *)(RuleList[i].edges),MAXEDGES,sizeof(EDGEDATA),CompareEdgesForSides);

  printf("\n#######################\n");
  printf("%d rules generated\n%d rules theoretically\n",nRules,nRulesTh);
  printf("#######################\n");

  if (SaveList)
  {
    /* open file */
    stream = fopen("RefRules.data","w");
    if (stream==NULL)
    {
      printf("ERROR: could not open file 'RefRules.data'\n");
      return (1);
    }

    /* save RuleList and PatternToRefrule to file */
    err = fwrite(&nRules,sizeof(int),1,stream);                                                 /* no. of RefRules		*/
    if (err!=1) goto exit;
    err = fwrite(RuleList,nRules*sizeof(REFRULE),1,stream);                             /* RefRules                     */
    if (err!=1) goto exit;
    j = 1<<(EDGES+MAX_SIDES_OF_ELEM);
    err = fwrite(&j,sizeof(int),1,stream);                                                              /* no. of patterns		*/
    if (err!=1) goto exit;
    err = fwrite(PatternToRefrule,j*sizeof(SHORT),1,stream);                            /* pattern -> RefRules	*/
    if (err!=1) goto exit;
    fclose(stream);

    printf("\nRuleList saved to 'RefRules.data'\n");
  }

  printf("\n######## DONE #########\n");

  return (0);

exit:
  fclose(stream);
  return (1);
}

#endif
