// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugio.c														*/
/*																			*/
/* Purpose:   ug's grid input/output                                        */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de					            */
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment ugio
#endif

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "compiler.h"
#include "fileopen.h"
#include "heaps.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"
#include "bio.h"
#include "misc.h"
#include "ugstruct.h"

#include "devices.h"

#include "parallel.h"
#include "gm.h"
#include "algebra.h"
#include "misc.h"
#include "ugm.h"
#include "ugio.h"
#include "elements.h"
#include "shapes.h"
#include "rm.h"
#include "mgio.h"
#include "fifo.h"

/* include refine because of macros accessed  */
#include "refine.h"
#include "rm.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define BUFFERSIZE                              512     /* size of the general purpose text buff*/

#define HEADER_FMT               "# grid on level 0 for %s\n# saved %s\n# %s\n# %s\n"
#define BN_HEADER_FMT    "\n# boundary nodes\n"
#define BN_FMT               "bn %d"
#define IN_HEADER_FMT    "\n# inner nodes\n"
#define IN_FMT               "in "
#define IE_HEADER_FMT    "\n# elements\n"
#define IE_FMT               "ie "
#define EOL_FMT              ";\n"
#define EOF_FMT              "# end of file\n"

/* size of integer list storing processor-ids of copies of objects,sufficient at least for one cg_element or refinement */
#define ELEMPROCLISTSIZE        2000    /*  n_elem + n_edge + n_node														*/
/*  n_elem = 6+30 (6 h-ghosts, 30 v-ghosts)											*/
/*  n_edge = 12*30 (12 max_edges_of_elem) (30 probably 30 elements per edge)		*/
/*  n_node = 100*8 (8 max_node_of_elem) 100 probably elements per node)				*/
#define PROCLISTSIZE            ELEMPROCLISTSIZE*MAX_SONS

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

static INT gridpaths_set=FALSE;
static MGIO_RR_RULE *rr_rules;
static unsigned short *ProcList, *ActProcListPos;
static int foid,non;
static NODE **nid_n;
static ELEMENT **eid_e;
static int nparfiles;                                     /* nb of parallel files		*/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static INT WriteElementParInfo (ELEMENT *theElement, MGIO_PARINFO *pinfo);

/****************************************************************************/
/*																			*/
/* Function:  MGSetVectorClasses											*/
/*																			*/
/* Purpose:   Returns highest vector class of a dof on next level			*/
/*																			*/
/* Input:	  *theElement													*/
/*																			*/
/* Output:	  INT															*/
/*																			*/
/****************************************************************************/

INT MGSetVectorClasses (MULTIGRID *theMG)
{
  INT i;
  GRID *theGrid;
  ELEMENT *theElement;

  /* set vector classes */
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if (ClearVectorClasses(theGrid)) return (1);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (ECLASS(theElement)!=RED_CLASS && ECLASS(theElement)!=GREEN_CLASS) continue;
      if (SeedVectorClasses(theGrid,theElement)) return (1);
    }
    if (PropagateVectorClasses(theGrid)) return (1);
  }

  /* set NextVectorClasses */
  if (ClearNextVectorClasses(GRID_ON_LEVEL(theMG,TOPLEVEL(theMG)))) return (1);
  for (i=0; i<TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if (ClearNextVectorClasses(theGrid)) return (1);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (NSONS(theElement)==0) continue;
      if (ECLASS(SON(theElement,0))!=RED_CLASS && ECLASS(SON(theElement,0))!=GREEN_CLASS) continue;
      if (SeedNextVectorClasses(theGrid,theElement)) return (1);
    }
    if (PropagateNextVectorClasses(theGrid)) return (1);
  }

  return (0);
}

/****************************************************************************/
/*D
        SaveMultiGrid - Save complete multigrid structure in a text file

        SYNOPSIS:
        INT SaveMultiGrid (MULTIGRID *theMG, char *name, char *comment);

        PARAMETERS:
        .  theMG - pointer to multigrid
        .  name - name of the text file
        .  comment - to be included at beginning of file

        DESCRIPTION:
        This function saves the grid on level 0 to a text file.
        The text file can be used in a script to load the grid.

        RETURN VALUE:
        INT
        .n    0 if ok
        .n    >0 if error occured.
        D*/
/****************************************************************************/

static INT SaveSurfaceGrid  (MULTIGRID *theMG, FILE *stream)
{
  NODE *theNode;
  ELEMENT *theElement;
  VERTEX *theVertex;
  DOUBLE *global;
  char buffer[BUFFERSIZE];
  INT i,id,move,l,tl,part;

  tl = CURRENTLEVEL(theMG);
  for (l=0; l<= tl; l++)
    for (theElement = FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
         theElement != NULL; theElement = SUCCE(theElement))
      if (NSONS(theElement) == 0)
        for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
          ID(MYVERTEX(CORNER(theElement,i))) = 0;

  /* find all boundary nodes witch are no corner nodes */
  fprintf(stream,BN_HEADER_FMT);
  id = 0;
  for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,0));
       theNode!= NULL; theNode=SUCCN(theNode)) {
    theVertex = MYVERTEX(theNode);
    if (OBJT(theVertex) == IVOBJ)
      continue;
    if (BNDP_BndPDesc(V_BNDP(theVertex),&move,&part))
      RETURN(1);
    if (move == 0)
      ID(theVertex) = id++;
  }
  for (l=0; l<= tl; l++)
    for (theElement = FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
         theElement != NULL; theElement = SUCCE(theElement))
      if (NSONS(theElement) == 0)
        for (i=0; i<CORNERS_OF_ELEM(theElement); i++) {
          theVertex = MYVERTEX(CORNER(theElement,i));
          if (OBJT(theVertex) == IVOBJ)
            continue;
          /* skip corner points */
          if (BNDP_BndPDesc(V_BNDP(theVertex),&move,&part))
            RETURN(1);
          if (move == 0)
            continue;
          if (ID(theVertex) > 0)
            continue;
          ID(theVertex) = id++;
          if (BNDP_SaveInsertedBndP(V_BNDP(theVertex),
                                    buffer,BUFFERSIZE))
            RETURN(1);
          fprintf(stream,"%s",buffer);
          fprintf(stream,EOL_FMT);
        }
  /* find all inner nodes */
  fprintf(stream,IN_HEADER_FMT);
  for (l=0; l<= tl; l++)
    for (theElement = FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
         theElement != NULL; theElement = SUCCE(theElement))
      if (NSONS(theElement) == 0)
        for (i=0; i<CORNERS_OF_ELEM(theElement); i++) {
          theVertex = MYVERTEX(CORNER(theElement,i));
          if (OBJT(theVertex) == BVOBJ)
            continue;
          if (ID(theVertex) > 0)
            continue;
          global = CVECT(theVertex);
          fprintf(stream,IN_FMT);
          for (i=0; i<DIM; i++)
            fprintf(stream," %f",global[i]);
          fprintf(stream,EOL_FMT);
          ID(theVertex) = id++;
        }

  /* elements */
  fprintf(stream,IE_HEADER_FMT);
  for (l=0; l<= tl; l++)
    for (theElement = FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
         theElement != NULL; theElement = SUCCE(theElement))
      if (NSONS(theElement) == 0) {
        fprintf(stream,IE_FMT);
        for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
          fprintf(stream," %d",ID(MYVERTEX(CORNER(theElement,i))));
        fprintf(stream,EOL_FMT);
      }

  /* trailer */
  fprintf(stream,EOF_FMT);
  fclose(stream);
  return(GM_OK);
}

static INT SaveMultiGrid_SCR (MULTIGRID *theMG, char *name, char *comment)
{
  FILE *stream;
  GRID *theGrid;
  NODE *theNode;
  ELEMENT *theElement;
  VERTEX *theVertex;
  DOUBLE *global;
  time_t Time;
  char *fmt;
  char buffer[BUFFERSIZE];
  BVP_DESC theBVPDesc;
  INT i,id,move,part;

  if (gridpaths_set)
    /* this way grids are stored to path[0] */
    stream = FileOpenUsingSearchPaths(name,"w","gridpaths");
  else
    stream = fileopen(name,"w");
  if (stream==NULL)
  {
    PrintErrorMessage('E',"SaveMultiGrid","cannot open file");
    RETURN(GM_FILEOPEN_ERROR);
  }

  /* get BVPDesc */
  if (BVP_SetBVPDesc(MG_BVP(theMG),&theBVPDesc))
    RETURN (GM_ERROR);

  /* get time */
  fmt = "%a %b %d %H:%M:%S %Y";
  time(&Time);
  strftime(buffer,BUFFERSIZE,fmt,localtime(&Time));

  /* write header */
  fprintf(stream,HEADER_FMT,BVPD_NAME(&theBVPDesc),buffer,name,comment);

  if (TOPLEVEL(theMG) > 0)
    return(SaveSurfaceGrid(theMG,stream));

  theGrid = GRID_ON_LEVEL(theMG,0);

  /* find all boundary nodes witch are no corner nodes */
  fprintf(stream,BN_HEADER_FMT);
  id = 0;
  for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    theVertex = MYVERTEX(theNode);
    if (OBJT(theVertex) == IVOBJ)
      continue;
    if (BNDP_BndPDesc(V_BNDP(theVertex),&move,&part))
      RETURN(1);
    if (move == 0)
      ID(theNode) = id++;
  }
  for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    theVertex = MYVERTEX(theNode);
    if (OBJT(theVertex) == IVOBJ)
      continue;
    /* skip corner points */
    if (BNDP_BndPDesc(V_BNDP(theVertex),&move,&part))
      RETURN(1);
    if (move == 0)
      continue;
    if (BNDP_SaveInsertedBndP(V_BNDP(theVertex),buffer,BUFFERSIZE))
      RETURN(1);
    fprintf(stream,"%s",buffer);
    fprintf(stream,EOL_FMT);
    ID(theNode) = id++;
  }
  /* find all inner nodes */
  fprintf(stream,IN_HEADER_FMT);
  for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    theVertex = MYVERTEX(theNode);
    if (OBJT(theVertex) == BVOBJ)
      continue;
    global = CVECT(theVertex);
    fprintf(stream,IN_FMT);
    for (i=0; i<DIM; i++)
      fprintf(stream," %f",global[i]);
    fprintf(stream,EOL_FMT);
    ID(theNode) = id++;
  }
  if (id != theGrid->nNode)
    RETURN(1);

  /* elements */
  fprintf(stream,IE_HEADER_FMT);
  for (theElement=FIRSTELEMENT(theGrid); theElement!= NULL;
       theElement=SUCCE(theElement))
  {
    fprintf(stream,IE_FMT);
    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      fprintf(stream," %d",ID(CORNER(theElement,i)));
    fprintf(stream,EOL_FMT);
  }

  /* trailer */
  fprintf(stream,EOF_FMT);
  fclose(stream);
  return(GM_OK);
}

INT OrphanCons(MULTIGRID *theMG)
{
  INT i,j,error,orphan;
  GRID            *theGrid;
  ELEMENT *theElement;
  NODE            *theNode,*FatherNode;
  EDGE            *theEdge;

  error = 0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      orphan = 0;
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      {
        theNode = CORNER(theElement,j);
        switch (NTYPE(theNode))
        {
        case (CORNER_NODE) :
          FatherNode = (NODE *)NFATHER(theNode);
          if (FatherNode == NULL)
          {
            if (EGHOST(theElement)) orphan = 1;
            else if (LEVEL(theElement) != 0)
            {
              error++;
              printf(PFMT "OrphanCons(): ERROR: elem=" EID_FMTX " cornernode=" ID_FMTX " is orphannode\n",
                     me,EID_PRTX(theElement),ID_PRTX(theNode));
            }
          }
          else
            assert (SONNODE(FatherNode) == theNode);
          break;
        case (MID_NODE) :
          theEdge = (EDGE *)NFATHER(theNode);
          if (theEdge == NULL)
          {
            if (EGHOST(theElement)) orphan = 1;
            else if (LEVEL(theElement) != 0)
            {
              error++;
              printf(PFMT "OrphanCons(): ERROR: elem=" EID_FMTX " midnode=" ID_FMTX " is orphannode\n",
                     me,EID_PRTX(theElement),ID_PRTX(theNode));
            }
          }
          else
            assert (MIDNODE(theEdge) == theNode);
          break;
        default :
          break;
        }
      }
      if (orphan)
      {
        ELEMENT *succe          = SUCCE(theElement);
        ELEMENT *theFather      = EFATHER(theElement);

        assert(EGHOST(theElement) || LEVEL(theElement)==0);
        if (theFather != NULL)
        {
          int index = PRIO2INDEX(EPRIO(theElement));
          ELEMENT *Next = NULL;
          ASSERT(index!=-1 && index<2);

          if (SON(theFather,index) == theElement)
          {
            if (succe != NULL)
            {
              if (EFATHER(succe)==theFather)
                if (EPRIO(succe) == EPRIO(theElement))
                {
                  Next = succe;
                }
            }
            SET_SON(theFather,index,Next);
          }

          SETNSONS(theFather,NSONS(theFather)-1);
          SET_EFATHER(theElement,NULL);
          GRID_UNLINK_ELEMENT(theGrid,theElement);
          GRID_LINK_ELEMENT(theGrid,theElement,EPRIO(theElement));
          PRINTDEBUG(gm,1,(PFMT "OrphanCons(): new orphan elem=" EID_FMTX "\n",
                           me,EID_PRTX(theElement)));
        }
      }
    }
  }

  return(error);
}

static INT Write_RefRules (MULTIGRID *theMG, INT *RefRuleOffset)
{
  MGIO_RR_GENERAL rr_general;
  INT i,j,k,t,nRules;
  HEAP *theHeap;
  MGIO_RR_RULE *Refrule;
  REFRULE * ug_refrule;
  struct mgio_sondata *sonData;

  /* init */
  if (theMG==NULL) return (1);
  theHeap = MGHEAP(theMG);

  /* write refrules general */
  nRules = 0; RefRuleOffset[0] = 0;
  for (i=0; i<TAGS; i++)
  {
    nRules += MaxRules[i];
    if (i>0) RefRuleOffset[i] = RefRuleOffset[i-1] +  MaxRules[i-1];
    rr_general.RefRuleOffset[i] = RefRuleOffset[i];
  }
  rr_general.nRules = nRules;
  if (Write_RR_General(&rr_general)) return (1);

  /* write refrules */
  rr_rules = Refrule = (MGIO_RR_RULE *)GetTmpMem(theHeap,nRules*sizeof(MGIO_RR_RULE));
  if (Refrule==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for refrules\n",(int)nRules*sizeof(MGIO_RR_RULE)); return (1);}
  for (t=0; t<TAGS; t++)
  {
    ug_refrule = RefRules[t];
    for (i=0; i<MaxRules[t]; i++)
    {
      Refrule->class = ug_refrule->class;
      Refrule->nsons = ug_refrule->nsons;
      for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
        Refrule->pattern[j] = ug_refrule->pattern[j];
      for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
      {
        Refrule->sonandnode[j][0] = ug_refrule->sonandnode[j][0];
        Refrule->sonandnode[j][1] = ug_refrule->sonandnode[j][1];
      }
      for (j=0; j<Refrule->nsons; j++)
      {
        sonData = &(Refrule->sons[j]);
        sonData->tag = ug_refrule->sons[j].tag;
        for (k=0; k<MGIO_MAX_CORNERS_OF_ELEM; k++)
          sonData->corners[k] = ug_refrule->sons[j].corners[k];
        for (k=0; k<MGIO_MAX_SIDES_OF_ELEM; k++)
          sonData->nb[k] = ug_refrule->sons[j].nb[k];
        sonData->path = ug_refrule->sons[j].path;
      }
      Refrule++;
      ug_refrule++;
    }
  }
  if (Write_RR_Rules(nRules,rr_rules)) return (1);

  return (0);
}

INT GetOrderedSons (ELEMENT *theElement, NODE **NodeContext, ELEMENT **SonList, INT *nmax)
{
  INT i,j,k,l,nfound,found;
  REFRULE *theRule;
  ELEMENT *NonorderedSonList[MAX_SONS];
  NODE *theNode;

  nfound = nmax[0] = 0;
  theRule = RefRules[TAG(theElement)] + REFINE(theElement);
  if (GetAllSons(theElement,NonorderedSonList)) return (1);
  for (i=0; i<NSONS_OF_RULE(theRule); i++)
  {
    found=1;
    for (j=0; j<CORNERS_OF_TAG(SON_TAG_OF_RULE(theRule,i)); j++)
      if (NodeContext[SON_CORNER_OF_RULE(theRule,i,j)] == NULL)
      {
        found=0;
        break;
      }
    if (!found)
    {
      SonList[i] = NULL;
      continue;
    }

    /* identify (hopefully) an element of SonList */
    for (j=0; j<NSONS(theElement); j++)
    {
      found=0;
      for (l=0; l<CORNERS_OF_TAG(SON_TAG_OF_RULE(theRule,i)); l++)
      {
        theNode = NodeContext[SON_CORNER_OF_RULE(theRule,i,l)];
        for (k=0; k<CORNERS_OF_ELEM(NonorderedSonList[j]); k++)
          if (theNode==CORNER(NonorderedSonList[j],k))
          {
            found++;
            break;
          }
      }
      if (found==CORNERS_OF_TAG(SON_TAG_OF_RULE(theRule,i)))
      {
        SonList[i] = NonorderedSonList[j];
        nmax[0] = i+1;
        break;
      }
      else
        SonList[i] = NULL;
    }
  }

  return (0);
}

static INT SetRefinement (ELEMENT *theElement, NODE **NodeContext, ELEMENT *SonList[MAX_SONS], INT nmax, MGIO_REFINEMENT *refinement, INT *RefRuleOffset)
{
  REFRULE *theRule;
  NODE *theNode;
  INT i,j,n,sonRefined,sonex,nex;

  if (nmax==0) return (0);
  refinement->refrule = REFINE(theElement) + RefRuleOffset[TAG(theElement)];
  theRule = RefRules[TAG(theElement)] + REFINE(theElement);

  if (MGIO_PARFILE)
  {
    nex=0;
    for (i=0; i<nmax; i++)
      if (SonList[i] != NULL)
        for (j=0; j<CORNERS_OF_TAG(theRule->sons[i].tag); j++)
          nex |= (1<<theRule->sons[i].corners[j]);
  }
  else
    nex     = ~0;

  /* store new corner ids */
  n=0;
  for (i=0; i<CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement)+SIDES_OF_ELEM(theElement); i++)
    if (NodeContext[i]!=NULL && ((nex>>i)&0x1))
      refinement->newcornerid[n++] = ID(NodeContext[i]);
  if (NodeContext[CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement)]!=NULL)
    refinement->newcornerid[n++] = ID(NodeContext[CENTER_NODE_INDEX(theElement)]);
  refinement->nnewcorners = n;

  /* sons are refined ? */
  sonRefined=sonex=refinement->nbid_ex=0;
  for (i=0,n=0; i<nmax; i++)
  {
    if (SonList[i]==NULL) continue;
    if (REFINE(SonList[i])!=NO_REFINEMENT && NSONS(SonList[i])>0) sonRefined |= (1<<i);
    if (MGIO_PARFILE)
    {
      sonex |= (1<<i);
      if (WriteElementParInfo(SonList[i],&refinement->pinfo[i])) return (1);
      for (j=0; j<SIDES_OF_ELEM(SonList[i]); j++)
        if (NBELEM(SonList[i],j)!=NULL && EFATHER(NBELEM(SonList[i],j))==NULL)
        {
          refinement->nbid_ex |= (1<<i);
          refinement->nbid[i][j] = ID(NBELEM(SonList[i],j));
        }
        else
          refinement->nbid[i][j] = -1;
    }
  }
  refinement->sonref = sonRefined;
  if (MGIO_PARFILE) refinement->sonex = sonex;

  /* not movable at the moment */
  refinement->nmoved = 0;

  /* set refinement class */
  if (nmax>0)
    refinement->refclass = ECLASS(SonList[nmax-1]);

  return (0);
}

static INT CheckNodeContext (ELEMENT *theElement, NODE ** NodeContext)
{
  INT i, mark,nor,noc;

  mark = MARK(theElement);
  for (i=CORNERS_OF_ELEM(theElement); i<CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement); i++)
  {
    if (NODE_OF_RULE(theElement,mark,i-CORNERS_OF_ELEM(theElement))) nor = 1;
    else nor = 0;
    if(NodeContext[i]!=NULL) noc = 1;
    else noc = 0;
    if(nor!=noc)
      if (GetNodeContext(theElement,NodeContext)) return (1);
  }
  for (; i<CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement)+SIDES_OF_ELEM(theElement); i++)
    if (NodeContext[i]!=NULL)
      if (GetNodeContext(theElement,NodeContext)) return (1);

  return (0);
}

static INT SetHierRefinement (ELEMENT *theElement, MGIO_REFINEMENT *refinement, INT *RefRuleOffset)
{
  INT i,nmax;
  ELEMENT *theSon;
  ELEMENT *SonList[MAX_SONS];
  NODE *NodeContext[MAX_NEW_CORNERS_DIM+MAX_CORNERS_OF_ELEM];

  /* sequential part */
  if (REFINE(theElement)==NO_REFINEMENT) return (0);
  if (GetNodeContext(theElement,NodeContext)) return (1);
  if (GetOrderedSons(theElement,NodeContext,SonList,&nmax)) return (1);
  if (nmax==0) return (0);
  if (SetRefinement (theElement,NodeContext,SonList,nmax,refinement,RefRuleOffset)) return (1);
  if (Write_Refinement (refinement,rr_rules)) return (1);

  /* recursive call */
  for (i=0; i<nmax; i++)
  {
    theSon = SonList[i];
    if (theSon==NULL) continue;
    if (REFINE(theSon)!=NO_REFINEMENT)
      if (SetHierRefinement(theSon,refinement,RefRuleOffset)) return (1);
  }

  return (0);
}

static INT nHierElements (ELEMENT *theElement, INT *n)
{
  INT i;
  ELEMENT *SonList[MAX_SONS];

  if (REFINE(theElement)==NO_REFINEMENT) return (0);
  if (GetAllSons(theElement,SonList)) return (1);
  if (SonList[0]==NULL) return (0);
  (*n)++;
  for (i=0; SonList[i]!=NULL; i++)
    if (nHierElements(SonList[i],n)) return (1);

  return (0);
}

static INT WriteCG_Vertices (MULTIGRID *theMG, INT nov)
{
  INT i,j,n;
  MGIO_CG_POINT *cg_point,*cgp;
  HEAP *theHeap;
  VERTEX *theVertex;

  theHeap = MGHEAP(theMG);
  n = nov*(MGIO_CG_POINT_SIZE);
  cg_point = (MGIO_CG_POINT *)GetTmpMem(theHeap,n);
  if (cg_point==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_points\n",(int)n); return (1);}
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theVertex=PFIRSTVERTEX(GRID_ON_LEVEL(theMG,i)); theVertex!=NULL; theVertex=SUCCV(theVertex))
      if (ID(theVertex)<nov)
      {
        assert(USED(theVertex));
        cgp = MGIO_CG_POINT_PS(cg_point,ID(theVertex));
        for (j=0; j<MGIO_DIM; j++)
          cgp->position[j] = CVECT(theVertex)[j];
        if (MGIO_PARFILE)
        {
          cgp->level = LEVEL(theVertex);
          cgp->prio = 0;
        }
      }

  if (Write_CG_Points((int)nov,cg_point)) return (1);

  return (0);
}

#ifdef ModelP
static INT WriteElementParInfo (ELEMENT *theElement, MGIO_PARINFO *pinfo)
{
  INT i,j,k,s,n_max;
  int *pl;
  NODE *theNode;
  EDGE *theEdge;

  n_max = PROCLISTSIZE-(ActProcListPos-ProcList);

  s=0;
  pinfo->prio_elem = DDD_InfoPriority(PARHDRE(theElement));
  pinfo->ncopies_elem = DDD_InfoNCopies(PARHDRE(theElement));
  if (n_max<pinfo->ncopies_elem) RETURN (1);
  if (pinfo->ncopies_elem>0)
  {
    pl = DDD_InfoProcList(PARHDRE(theElement));
    for (i=0,j=2; i<pinfo->ncopies_elem; i++,j+=2)
      ActProcListPos[s++] = pl[j];
  }
  pinfo->e_ident = EGID(theElement);
  for (k=0; k<CORNERS_OF_ELEM(theElement); k++)
  {
    theNode = CORNER(theElement,k);
    pinfo->prio_node[k] = DDD_InfoPriority(PARHDR(theNode));
    pinfo->ncopies_node[k] = DDD_InfoNCopies(PARHDR(theNode));
    if (n_max<pinfo->ncopies_node[k]+s) RETURN (1);
    if (pinfo->ncopies_node[k]>0)
    {
      pl = DDD_InfoProcList(PARHDR(theNode));
      for (i=0,j=2; i<pinfo->ncopies_node[k]; i++,j+=2)
        ActProcListPos[s++] = pl[j];
    }
    pinfo->n_ident[k] = GID(CORNER(theElement,k));
  }
#ifdef __THREEDIM__
  for (k=0; k<EDGES_OF_ELEM(theElement); k++)
  {
    theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,k,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,k,1)));
    pinfo->prio_edge[k] = DDD_InfoPriority(PARHDR(theEdge));
    pinfo->ncopies_edge[k] = DDD_InfoNCopies(PARHDR(theEdge));
    if (n_max<pinfo->ncopies_edge[k]+s) RETURN (1);
    if (pinfo->ncopies_edge[k]>0)
    {
      pl = DDD_InfoProcList(PARHDR(theEdge));
      for (i=0,j=2; i<pinfo->ncopies_edge[k]; i++,j+=2)
        ActProcListPos[s++] = pl[j];
    }
    pinfo->ed_ident[k] = GID(GetEdge(CORNER_OF_EDGE_PTR(theElement,k,0),CORNER_OF_EDGE_PTR(theElement,k,1)));
  }
#endif
  if (s>0) pinfo->proclist = ActProcListPos;
  else pinfo->proclist = NULL;
  ActProcListPos += s;

  return (0);
}
#else
static INT WriteElementParInfo (ELEMENT *theElement, MGIO_PARINFO *pinfo)
{
  return (1);
}
#endif

static INT SaveMultiGrid_SPF (MULTIGRID *theMG, char *name, char * type, char *comment)
{
  GRID *theGrid;
  NODE *theNode;
  ELEMENT *theElement;
  HEAP *theHeap;
  MGIO_MG_GENERAL mg_general;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  MGIO_GE_GENERAL ge_general;
  MGIO_GE_ELEMENT ge_element[TAGS];
  MGIO_CG_GENERAL cg_general;
  MGIO_CG_POINT *cg_point;
  MGIO_CG_ELEMENT *cg_element,*cge;
  MGIO_REFINEMENT *refinement;
  MGIO_BD_GENERAL bd_general;
  MGIO_PARINFO cg_pinfo;
  INT i,j,k,niv,nbv,nie,nbe,n,nhe,hr_max,mode,level,n_inc,n_max,id,foid,non,tl,Err;
  INT RefRuleOffset[TAGS];
  int *vidlist,ftype,error,*cg_ident;
  char *p;
  BNDP **BndPList;
  char filename[NAMESIZE];
  char buf[64];

  /* check */
  if (theMG==NULL) return (1);
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap);

  /* open file */
  nparfiles = procs;
  if (strcmp(type,"dbg")==0) mode = BIO_DEBUG;
  else if (strcmp(type,"asc")==0) mode = BIO_ASCII;
  else if (strcmp(type,"bin")==0) mode = BIO_BIN;
  else return (1);
  strcpy(filename,name);
  sprintf(buf,".ug.mg.");
  strcat(filename,buf);
  strcat(filename,type);
#ifdef ModelP
  error = 0;
  if (me == master)
  {
    if (MGIO_PARFILE)
    {
      ftype = MGIO_filetype(filename);
      if (ftype == FT_FILE)
      {
        error = -1;
      }
      else if (ftype == FT_UNKNOWN)
      {
        if (MGIO_dircreate(filename)) error = -1;
      }
    }
  }
  Broadcast(&error,sizeof(int));
  if (error == -1)
  {
    UserWriteF("SaveMultiGrid_SPF(): error during file/directory creation\n");
    return(1);
  }
  if (MGIO_PARFILE)
  {
    sprintf(buf,"/mg.%04d",(int)me);
    strcat(filename,buf);
  }
#endif
  if (Write_OpenMGFile (filename)) return (1);

  /* write general information */
  theBVP = MG_BVP(theMG);
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc)) return (1);
  mg_general.mode                 = mode;
  mg_general.dim                  = DIM;
  mg_general.magic_cookie = MG_MAGIC_COOKIE(theMG);
  mg_general.heapsize             = MGHEAP(theMG)->size/1024;
  mg_general.nLevel               = TOPLEVEL(theMG) + 1;
  mg_general.nNode = mg_general.nPoint = mg_general.nElement = 0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    mg_general.nNode                += NN(theGrid);
    mg_general.nPoint               += NV(theGrid);
    mg_general.nElement             += NT(theGrid);
  }
  strcpy(mg_general.version,MGIO_VERSION);
  p = GetStringVar (":IDENTIFICATION");
  if (p!=NULL) strcpy(mg_general.ident,p);
  else strcpy(mg_general.ident,"---");
  strcpy(mg_general.DomainName,BVPD_NAME(&theBVPDesc));
  strcpy(mg_general.MultiGridName,MGNAME(theMG));
  strcpy(mg_general.Formatname,ENVITEM_NAME(MGFORMAT(theMG)));
  mg_general.VectorTypes  = 0;

  /* parallel part */
  mg_general.nparfiles = nparfiles;
  mg_general.me = me;

  if (Write_MG_General(&mg_general)) return (1);

  /* write information about general elements */
  ge_general.nGenElement = TAGS;
  if (Write_GE_General(&ge_general)) return (1);
  for (i=0; i<TAGS; i++)
  {
    ge_element[i].tag = i;
    if (element_descriptors[i]!=NULL)
    {
      ge_element[i].nCorner = element_descriptors[i]->corners_of_elem;
      ge_element[i].nEdge = element_descriptors[i]->edges_of_elem;
      ge_element[i].nSide = element_descriptors[i]->sides_of_elem;
      for (j=0; j<MGIO_MAX_EDGES_OF_ELEM; j++)
      {
        ge_element[i].CornerOfEdge[j][0] = element_descriptors[i]->corner_of_edge[j][0];
        ge_element[i].CornerOfEdge[j][1] = element_descriptors[i]->corner_of_edge[j][1];
      }
      for (j=0; j<MGIO_MAX_SIDES_OF_ELEM; j++)
        for (k=0; k<MGIO_MAX_CORNERS_OF_SIDE; k++)
          ge_element[i].CornerOfSide[j][k] = element_descriptors[i]->corner_of_side[j][k];
    }
    else
    {
      ge_element[i].nCorner = 0;
      ge_element[i].nEdge = 0;
      ge_element[i].nSide = 0;
    }
  }
  if (Write_GE_Elements(TAGS,ge_element)) return (1);

  /* write information about refrules used */
  if (Write_RefRules(theMG,RefRuleOffset)) return (1);

  if (MGIO_PARFILE)
  {
    if (i = OrphanCons(theMG))
    {
      printf(PFMT "OrphanCons() returned %d errors\n",me,i);
      fflush(stdout);
      assert(0);
    }

  }

  /* renumber objects */
  if (RenumberMultiGrid (theMG,&nbe,&nie,&nbv,&niv,&foid,&non)) {UserWriteF("ERROR: cannot renumber multigrid\n"); return (1);}

  /* write general information about coarse grid */
  theGrid = GRID_ON_LEVEL(theMG,0);
  cg_general.nPoint = nbv+niv;
  cg_general.nBndPoint = nbv;
  cg_general.nInnerPoint = niv;
  cg_general.nElement = nbe+nie;
  cg_general.nBndElement = nbe;
  cg_general.nInnerElement = nie;
  if (Write_CG_General(&cg_general)) return (1);

  /* write coarse grid points */
  if (WriteCG_Vertices(theMG,nbv+niv)) return (1);

  /* write mapping: node-id --> vertex-id for orphan-nodes, if(MGIO_PARFILE) */
  if (MGIO_PARFILE)
  {
    vidlist = (int*)GetTmpMem(theHeap,(2+non)*sizeof(int));
    vidlist[0] = non;
    vidlist[1] = foid;
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (USED(theNode))
          vidlist[ID(theNode)+2-foid] = ID(MYVERTEX(theNode));
    if (Bio_Write_mint((int)(2+non),vidlist)) return(1);
  }

  /* write orphan elements */
  n = cg_general.nElement; hr_max=0;
  cg_element = (MGIO_CG_ELEMENT*)GetTmpMem(theHeap,n*MGIO_CG_ELEMENT_SIZE);
  if (cg_element==NULL)   {UserWriteF("ERROR: cannot allocate %d bytes for cg_elements\n",(int)n*MGIO_CG_ELEMENT_SIZE); return (1);}
  if (MGIO_PARFILE)
  {
    ActProcListPos = ProcList = (unsigned short*)GetTmpMem(theHeap,PROCLISTSIZE*sizeof(unsigned short));
    if (ProcList==NULL)     {UserWriteF("ERROR: cannot allocate %d bytes for ProcList\n",(int)PROCLISTSIZE*sizeof(int)); return (1);}
  }
  else
    ActProcListPos = ProcList = NULL;
  id=i=0;
  for (level=0; level<=TOPLEVEL(theMG); level++)
    for (theElement = PFIRSTELEMENT(GRID_ON_LEVEL(theMG,level)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      /* entry in cg_element-list */
      cge = MGIO_CG_ELEMENT_PS(cg_element,i);

      /* only orphan elements */
      if (EFATHER(theElement)!=NULL) continue;
      assert(id==ID(theElement));
      assert(id<n);

      if (MGIO_PARFILE)
      {
        /* parallel part */
        cge->level = LEVEL(theElement);
      }

      /* --- */
      cge->ge = TAG(theElement);
      nhe=0;
      if (nHierElements (theElement,&nhe)) return (1);
      hr_max = MAX(hr_max,nhe);
      cge->nhe = nhe;
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
        cge->cornerid[j] = ID(CORNER(theElement,j));
      for (j=0; j<SIDES_OF_ELEM(theElement); j++)
        if (NBELEM(theElement,j)!=NULL)
          cge->nbid[j] = ID(NBELEM(theElement,j));
        else
          cge->nbid[j] = -1;
      cge->subdomain = SUBDOMAIN(theElement);

      /* increment counters */
      i++; id++;                                                                                                        /* id this is the id of theElement, according to RenumberEments (ugm.c) */
    }
  if (Write_CG_Elements((int)i,cg_element)) return (1);

  /* write general bnd information */
  if (Bio_Jump_From ()) return (1);
  bd_general.nBndP = nbv;
  if (Write_BD_General (&bd_general)) return (1);

  /* write bnd information */
  BndPList = (BNDP**)GetTmpMem(theHeap,nbv*sizeof(BNDP*));
  if (BndPList==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for BndPList\n",(int)nbv*sizeof(BNDP*)); return (1);}
  if (procs>1) i=TOPLEVEL(theMG);
  else i=0;
  for (level=0; level<=i; level++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,level)); theNode!=NULL; theNode=SUCCN(theNode))
      if (USED(MYVERTEX(theNode)) && OBJT(MYVERTEX(theNode))==BVOBJ)                            /* see ugm.c: RenumberVertices */
      {
        id = ID(MYVERTEX(theNode));
        if (id<0 || id>=nbv) return (1);
        BndPList[id] = V_BNDP(MYVERTEX(theNode));
        SETUSED(MYVERTEX(theNode),0);
      }
  if (Write_PBndDesc (nbv,BndPList)) return (1);
  if (Bio_Jump_To ()) return (1);

  /* write parinfo of coarse-grid */
  if (MGIO_PARFILE)
  {
    cg_pinfo.proclist = ProcList;
    for (level=0; level<=TOPLEVEL(theMG); level++)
      for (theElement = PFIRSTELEMENT(GRID_ON_LEVEL(theMG,level)); theElement!=NULL; theElement=SUCCE(theElement))
        if (EFATHER(theElement)==NULL)
        {
          if (WriteElementParInfo(theElement,&cg_pinfo)) RETURN (1);
          if (Write_pinfo (TAG(theElement),&cg_pinfo)) RETURN (1);
        }
  }

  /* save refinement */
  refinement = (MGIO_REFINEMENT *)GetTmpMem(theHeap,(hr_max+1)*MGIO_REFINEMENT_SIZE);                   /* size according to procs==1 or procs>1 (see mgio.h) */
  if (refinement==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for refinement\n",(int)hr_max*sizeof(MGIO_REFINEMENT)); return (1);}
  if (procs>1) tl=TOPLEVEL(theMG);
  else tl=0;
  id=0;
  for (level=0; level<=tl; level++)
    for (theElement = PFIRSTELEMENT(GRID_ON_LEVEL(theMG,level)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (EFATHER(theElement)!=NULL) continue;
      assert(id==ID(theElement));
      if (SetHierRefinement (theElement,refinement,RefRuleOffset)) return (1);
      id++;
    }

  /* release tmp mem */
  ReleaseTmpMem(theHeap);

  /* close file */
  if (CloseMGFile ()) return (1);

  /* saved */
  MG_SAVED(theMG) = 1;
  strcpy(MG_FILENAME(theMG),filename);

  return (0);
}

INT SaveMultiGrid (MULTIGRID *theMG, char *name, char *type, char *comment)
{
  char *p;

  /* check name convention */
  p = name + strlen(name) - 4;
  if (strcmp(p,".scr")==0)
  {
    if (SaveMultiGrid_SCR (theMG,name,comment)) return (1);
  }
  else
  {
    if (SaveMultiGrid_SPF (theMG,name,type,comment)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*
   LoadMultiGrid - Load complete multigrid structure from a text file

   SYNOPSIS:
   MULTIGRID *LoadMultiGrid (char *MultigridName, char *FileName,
   char *BVPName, char *format, unsigned long heapSize);

   PARAMETERS:
   .  MultigridName - Name of the new 'MULTIGRID' structure in memory.
   .  FileName - Name of the file to be read.
   .  BVPName - `Name` of the BVP used for the 'MULTIGRID'.
   .  format - `Name` of the 'FORMAT' to be used for the 'MULTIGRID'.
   .  heapSize - Size of the heap in bytes that will be allocated for the 'MULTIGRID'.

   DESCRIPTION:
   This function can read grid files produced with the 'SaveMultiGrid' function.

   RETURN VALUE:
   INT
   .n    NULL if an error occured
   .n    else pointer to new 'MULTIGRID'
 */
/****************************************************************************/

static INT Evaluate_pinfo (GRID *theGrid, ELEMENT *theElement, MGIO_PARINFO *pinfo)
{
  INT i,j,k,s,prio,where,oldwhere,old;
  INT evec,nvec,edvec,svec;
  ELEMENT         *theFather,*After,*Next,*Succe;
  NODE            *theNode;
  EDGE            *theEdge;
  VECTOR          *theVector;


  evec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),ELEMVEC);
  nvec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),NODEVEC);
  edvec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),EDGEVEC);
#ifdef __THREEDIM__
  svec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),SIDEVEC);
#else
  svec = 0;
#endif
  /* this funxtion does not support side vectors                        */
  /* proclist and identificator need to be stored for each  side vector */
  if (svec) assert(0);

  theFather = EFATHER(theElement);
  s = 0;
  if ((prio = pinfo->prio_elem) != PrioMaster)
  {
    old = EPRIO(theElement);
    oldwhere = PRIO2INDEX(old);
    GRID_UNLINK_ELEMENT(theGrid,theElement);
    SETEPRIO(theElement,prio);
    if (theFather != NULL)
    {
      if (theElement == SON(theFather,oldwhere))
      {
        Next = NULL;
        Succe = SUCCE(theElement);
        if (Succe != NULL)
          if (EFATHER(Succe)==theFather && EPRIO(Succe)==old)
            Next = Succe;
        SET_SON(theFather,oldwhere,Next);
      }
      where = PRIO2INDEX(prio);
      After = SON(theFather,where);
      if (After == NULL) SET_SON(theFather,where,theElement);
      GRID_LINKX_ELEMENT(theGrid,theElement,prio,After);
    }
    else
      GRID_LINK_ELEMENT(theGrid,theElement,prio);
    if (evec)
    {
      theVector = EVECTOR(theElement);
      GRID_UNLINK_VECTOR(theGrid,theVector);
      SETPRIO(EVECTOR(theElement),prio);
      GRID_LINK_VECTOR(theGrid,theVector,prio);
    }
  }
  for (i=0; i<pinfo->ncopies_elem; i++)
  {
    DDD_IdentifyNumber(PARHDRE(theElement),pinfo->proclist[s],pinfo->e_ident);
    if (evec)
      DDD_IdentifyNumber(PARHDR(EVECTOR(theElement)),pinfo->proclist[s],pinfo->e_ident);
    s++;
  }

  for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
  {
    theNode = CORNER(theElement,j);
    if (USED(theNode) == 0)
    {
      if ((prio = pinfo->prio_node[j]) != PrioMaster)
      {
        GRID_UNLINK_NODE(theGrid,theNode);
        SETPRIO(theNode,prio);
        GRID_LINK_NODE(theGrid,theNode,prio);
        if (nvec)
        {
          theVector = NVECTOR(theNode);
          GRID_UNLINK_VECTOR(theGrid,theVector);
          SETPRIO(NVECTOR(theNode),prio);
          GRID_LINK_VECTOR(theGrid,theVector,prio);
        }
      }
      /*printf("Evaluate-pinfo():nid=%d prio=%d\n",ID(theNode),prio);fflush(stdout);*/
      for (i=0; i<pinfo->ncopies_node[j]; i++)
      {
        DDD_IdentifyNumber(PARHDR(theNode),pinfo->proclist[s],pinfo->n_ident[j]);
        if (nvec)
          DDD_IdentifyNumber(PARHDR(NVECTOR(theNode)),pinfo->proclist[s],pinfo->n_ident[j]);
        s++;
      }
      SETUSED(theNode,1);
    }
    else
      s += pinfo->ncopies_node[j];
  }
#if (MGIO_DIM==3)
  for (j=0; j<EDGES_OF_ELEM(theElement); j++)
  {
    theEdge = GetEdge(CORNER_OF_EDGE_PTR(theElement,j,0),
                      CORNER_OF_EDGE_PTR(theElement,j,1));
    if (USED(theEdge) == 0)
    {
      if ((prio = pinfo->prio_edge[j]) != PrioMaster)
      {
        SETPRIO(theEdge,prio);
        if (edvec)
        {
          GRID_UNLINK_VECTOR(theGrid,theVector);
          SETPRIO(EDVECTOR(theEdge),prio);
          GRID_LINK_VECTOR(theGrid,theVector,prio);
        }
      }
      for (i=0; i<pinfo->ncopies_edge[j]; i++)
      {
        DDD_IdentifyNumber(PARHDR(theEdge),pinfo->proclist[s],pinfo->ed_ident[j]);
        if (edvec)
          DDD_IdentifyNumber(PARHDR(EDVECTOR(theEdge)),pinfo->proclist[s],pinfo->ed_ident[j]);
        s++;
      }
      SETUSED(theEdge,1);
    }
    else
      s += pinfo->ncopies_edge[j];
  }
#endif

  return(0);
}

static INT IO_GridCons(MULTIGRID *theMG)
{
  INT i,*proclist;
  GRID    *theGrid;
  ELEMENT *theElement;
  VECTOR  *theVector;

  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      proclist = EPROCLIST(theElement);
      while (proclist[0] != -1)
      {
        if (EMASTERPRIO(proclist[1])) PARTITION(theElement) = proclist[0];
        proclist += 2;
      }
    }
    for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
      if (!MASTER(theVector))
        DisposeConnectionFromVector(theGrid,theVector);
    /* spread nodetypes from master to its copies */
    assert(0);
  }

  return(GM_OK);
}

static INT InsertLocalTree (GRID *theGrid, ELEMENT *theElement, MGIO_REFINEMENT *ref, int *RefRuleOffset)
{
  INT i,j,k,r_index,nedge,type,sonRefined,n0,n1,Sons_of_Side,SonSides[MAX_SONS],offset;
  ELEMENT *theSonElem[MAX_SONS];
  ELEMENT *theSonList[MAX_SONS];
  NODE *NodeList[MAX_NEW_CORNERS_DIM+MAX_CORNERS_OF_ELEM];
  NODE *SonNodeList[MAX_CORNERS_OF_ELEM];
  NODE *theNode;
  GRID *upGrid;
  EDGE *theEdge;
  MGIO_RR_RULE *theRule;
  struct mgio_sondata *SonData;
  INT nbside,nex;

  /*PRINTDEBUG(gm,0,(PFMT "InsertLocalTree(): level=%d elem=" EID_FMTX "\n",me,LEVEL(theGrid),EID_PRTX(theElement)));*/

  /* read refinement */
  if (Read_Refinement(ref,rr_rules)) RETURN (1);

  /* init */
  if (ref->refrule==-1) RETURN (1);
  SETREFINE(theElement,ref->refrule-RefRuleOffset[TAG(theElement)]);
  SETREFINECLASS(theElement,ref->refclass);
  SETMARK(theElement,ref->refrule-RefRuleOffset[TAG(theElement)]);
  SETMARKCLASS(theElement,ref->refclass);
  theRule = rr_rules+ref->refrule;
  upGrid = UPGRID(theGrid);

  /* insert nodes */
  if (MGIO_PARFILE)
  {
    nex=0;
    for (i=0; i<theRule->nsons; i++)
      if ((ref->sonex>>i)&1)
        for (j=0; j<CORNERS_OF_TAG(theRule->sons[i].tag); j++)
          nex |= (1<<theRule->sons[i].corners[j]);
  }
  r_index = 0;
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    if (MGIO_PARFILE && !((nex>>i)&1))
    {
      NodeList[i] = NULL;
      continue;
    }
    NodeList[i] = SONNODE(CORNER(theElement,i));
    if (NodeList[i]==NULL)
    {
      if (ref->newcornerid[r_index]-foid>=0 && ref->newcornerid[r_index]-foid<non)
      {
        NodeList[i] = nid_n[ref->newcornerid[r_index]-foid];
        assert(NodeList[i]!=NULL);
        assert(ID(NodeList[i]) == ref->newcornerid[r_index]);
        SONNODE(CORNER(theElement,i)) = NodeList[i];
        SETNFATHER(NodeList[i],(GEOM_OBJECT *)CORNER(theElement,i));
      }
      else
      {
        NodeList[i] = CreateSonNode(upGrid,CORNER(theElement,i));
        if (NodeList[i]==NULL) RETURN (1);
        ID(NodeList[i]) = ref->newcornerid[r_index];
      }
    }
    else if (ID(NodeList[i]) != ref->newcornerid[r_index])
    {
      printf("WARNING: inconsistent CORNER-node-ids(elem=%d): %d != %d\n",ID(theElement),(int)ID(NodeList[i]),(int)ref->newcornerid[r_index]);
      printf("pos: %f %f\n",(float)CVECT(MYVERTEX(NodeList[i]))[0],(float)CVECT(MYVERTEX(NodeList[i]))[1]);
    }
    SETNTYPE(NodeList[i],CORNER_NODE);
    r_index++;
  }
  offset = i;

  nedge = EDGES_OF_ELEM(theElement);
  for (; i<nedge+offset; i++)
  {
    if (MGIO_PARFILE && !((nex>>i)&1))
    {
      NodeList[i] = NULL;
      continue;
    }
    if (theRule->pattern[i-offset]!=1)
    {
      NodeList[i] = NULL;
      continue;
    }
    n0 = CORNER_OF_EDGE(theElement,i-offset,0);
    n1 = CORNER_OF_EDGE(theElement,i-offset,1);
    theEdge = GetEdge(CORNER(theElement,n0),CORNER(theElement,n1));
    if (theEdge==NULL) RETURN (1);
    NodeList[i] = MIDNODE(theEdge);
    if (NodeList[i]==NULL)
    {
      if (ref->newcornerid[r_index]-foid>=0 && ref->newcornerid[r_index]-foid<non)
      {
        NodeList[i] = nid_n[ref->newcornerid[r_index]-foid];
        assert(NodeList[i]!=NULL);
        assert(ID(NodeList[i]) == ref->newcornerid[r_index]);
        MIDNODE(theEdge) = NodeList[i];
        SETNFATHER(NodeList[i],(GEOM_OBJECT *)theEdge);
      }
      else
      {
        NodeList[i] = CreateMidNode(upGrid,theElement,i-offset);
        if (NodeList[i]==NULL) RETURN (1);
        ID(NodeList[i]) = ref->newcornerid[r_index];
        HEAPFAULT(NodeList[i]);
      }
    }
    else if (ID(NodeList[i]) != ref->newcornerid[r_index])
    {
      printf("WARNING: inconsistent EDGE-node-ids(elem=%d): %d != %d\n",ID(theElement),(int)ID(NodeList[i]),(int)ref->newcornerid[r_index]);
      printf("pos: %f %f\n",(float)CVECT(MYVERTEX(NodeList[i]))[0],(float)CVECT(MYVERTEX(NodeList[i]))[1]);
    }
    SETNTYPE(NodeList[i],MID_NODE);
    HEAPFAULT(NodeList[i]);
    r_index++;
  }
  offset = i;

#ifdef __THREEDIM__
  for (; i<SIDES_OF_ELEM(theElement)+offset; i++)
  {
    if (MGIO_PARFILE && !((nex>>i)&1))
    {
      NodeList[i] = NULL;
      continue;
    }
    if (theRule->pattern[i-CORNERS_OF_ELEM(theElement)]!=1)
    {
      NodeList[i] = NULL;
      continue;
    }

    if (ref->newcornerid[r_index]-foid>=0 && ref->newcornerid[r_index]-foid<non)
    {
      NodeList[i] = nid_n[ref->newcornerid[r_index]-foid];
      assert(NodeList[i]!=NULL);
      assert(ID(NodeList[i]) == ref->newcornerid[r_index]);
    }
    else
      NodeList[i] = GetSideNode(theElement,i-offset);
    if (NodeList[i]==NULL)
    {
      NodeList[i] = CreateSideNode(upGrid,theElement,i-offset);
      if (NodeList[i]==NULL) RETURN (1);
      ID(NodeList[i]) = ref->newcornerid[r_index];
    }
    else
      SETONNBSIDE(MYVERTEX(NodeList[i]),i-offset);
    SETNTYPE(NodeList[i],SIDE_NODE);
    r_index++;
  }
  offset = i;
#endif

  if (theRule->pattern[CENTER_NODE_INDEX(theElement)]==1 && (!MGIO_PARFILE || ((nex>>i)&1)))
  {
    if (ref->newcornerid[r_index]-foid>=0 && ref->newcornerid[r_index]-foid<non)
    {
      NodeList[i] = nid_n[ref->newcornerid[r_index]-foid];
      assert(NodeList[i]!=NULL);
      assert(ID(NodeList[i]) == ref->newcornerid[r_index]);
    }
    else
      NodeList[i] = CreateCenterNode(upGrid,theElement);
    if (NodeList[i]==NULL) RETURN (1);
    SETNTYPE(NodeList[i],CENTER_NODE);
    ID(NodeList[i]) = ref->newcornerid[r_index++];
  }
  else
    NodeList[i] = NULL;

  /* insert elements */
  for (i=0; i<theRule->nsons; i++)
  {
    if (MGIO_PARFILE && !((ref->sonex>>i)&1))
    {
      theSonList[i] = NULL;
      continue;
    }
    SonData = theRule->sons+i;
    type = IEOBJ;
    if (OBJT(theElement)==BEOBJ)
      for (j=0; j<SIDES_OF_TAG(SonData->tag); j++)
        if (SonData->nb[j]>=20 && SIDE_ON_BND(theElement,SonData->nb[j]-20))
        {
          type = BEOBJ;
          break;
        }
    for (j=0; j<CORNERS_OF_TAG(SonData->tag); j++)
      SonNodeList[j] = NodeList[SonData->corners[j]];
    theSonList[i] = CreateElement(upGrid,SonData->tag,type,SonNodeList,theElement);
    if (theSonList[i]==NULL) RETURN (1);
    SETECLASS(theSonList[i],ref->refclass);
    SETSUBDOMAIN(theSonList[i],SUBDOMAIN(theElement));
  }

  /* identify objects */
  for (i=0; i<theRule->nsons; i++)
  {
    if (theSonList[i] != NULL)
      Evaluate_pinfo(upGrid,theSonList[i],ref->pinfo+i);
  }

  /* set neighbor relation between sons */
  for (i=0; i<theRule->nsons; i++)
  {
    if (MGIO_PARFILE && !((ref->sonex>>i)&1)) continue;
    SonData = theRule->sons+i;
    for (j=0; j<SIDES_OF_TAG(SonData->tag); j++)
      if (SonData->nb[j]<20 && (!MGIO_PARFILE || ((ref->sonex>>(SonData->nb[j]))&1)))
        SET_NBELEM(theSonList[i],j,theSonList[SonData->nb[j]]);
      else
        SET_NBELEM(theSonList[i],j,NULL);
  }

  /* connect to neighbors */
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    for (k=0,j=0; j<theRule->nsons; j++)
    {
      if (!MGIO_PARFILE || ((ref->sonex>>j)&1))
        theSonElem[k++] = theSonList[j];
    }
    theSonElem[k] = NULL;
    if (Get_Sons_of_ElementSide (theElement,i,&Sons_of_Side,theSonElem,SonSides,0)) RETURN (1);
    if (Connect_Sons_of_ElementSide(UPGRID(theGrid),theElement,i,Sons_of_Side,theSonElem,SonSides,0,1)) RETURN (1);
  }

  /* connect to orphan-neighbors */
  if (MGIO_PARFILE)
    for (i=0; i<theRule->nsons; i++)
      if (((ref->sonex>>i)&1) && ((ref->nbid_ex>>i)&1))
      {
        SonData = theRule->sons+i;
        for (j=0; j<SIDES_OF_TAG(SonData->tag); j++)
          if (ref->nbid[i][j]>=0)
          {
            SET_NBELEM(theSonList[i],j,eid_e[ref->nbid[i][j]]);
            GetNbSideByNodes(eid_e[ref->nbid[i][j]],&nbside,theSonList[i],j);
            SET_NBELEM(eid_e[ref->nbid[i][j]],nbside,theSonList[i]);
          }
      }

  /* jump to the sons ? */
  sonRefined = ref->sonref;

  /* call recoursively */
  for (i=0; i<theRule->nsons; i++)
    if (sonRefined & (1<<i))
    {
      if (InsertLocalTree (upGrid,theSonList[i],ref,RefRuleOffset)) RETURN (1);
    }
    else if (theSonList[i] != NULL)
    {
      SETREFINE(theSonList[i],NO_REFINEMENT);
    }

  return (0);
}

MULTIGRID *LoadMultiGrid (char *MultigridName, char *name, char *type, char *BVPName, char *format, unsigned long heapSize, INT force, INT optimizedIE)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  ELEMENT *theElement,*theNeighbor,*ENext;
  NODE *theNode;
  EDGE *theEdge;
  HEAP *theHeap;
  MGIO_MG_GENERAL mg_general;
  MGIO_GE_GENERAL ge_general;
  MGIO_GE_ELEMENT ge_element[TAGS];
  MGIO_RR_GENERAL rr_general;
  MGIO_CG_GENERAL cg_general;
  MGIO_CG_POINT *cg_point,*cgp;
  MGIO_CG_ELEMENT *cg_element,*cge;
  MGIO_BD_GENERAL bd_general;
  MGIO_PARINFO cg_pinfo;
  MGIO_REFINEMENT *refinement;
  BNDP **BndPList;
  DOUBLE *Positions;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  MESH theMesh;
  char FormatName[NAMESIZE], BndValName[NAMESIZE], MGName[NAMESIZE], filename[NAMESIZE];
  INT i,j,k,*Element_corner_uniq_subdom, *Ecusdp[2],**Enusdp[2],**Ecidusdp[2],
  **Element_corner_ids_uniq_subdom,*Element_corner_ids,max,**Element_nb_uniq_subdom,
  *Element_nb_ids,id,level;
  INT parfile;
  char buf[64];
  int *vidlist,fileprocs;

  /* open file */
  strcpy(filename,name);
  strcat(filename,".ug.mg.");
  strcat(filename,type);
#ifdef ModelP
  if (me == master)
  {
#endif
  nparfiles = 1;
  if (MGIO_filetype(filename) == FT_DIR)
  {
    sprintf(buf,"/mg.%04d",(int)me);
    strcat(filename,buf);
    if (Read_OpenMGFile (filename))                                                                         {nparfiles = -1;}
    else
    if (Read_MG_General(&mg_general))                                                               {CloseMGFile (); nparfiles = -1;}
    nparfiles = mg_general.nparfiles;
  }
  else if(MGIO_filetype(filename) == FT_FILE)
  {
    if (Read_OpenMGFile (filename))                                                                         {nparfiles = -1;}
    else
    if (Read_MG_General(&mg_general))                                                               {CloseMGFile (); nparfiles = -1;}
  }
  else
    nparfiles = -1;
#ifdef ModelP
  Broadcast(&nparfiles,sizeof(int));
}
else
{
  Broadcast(&nparfiles,sizeof(int));
  if (me < nparfiles)
  {
    sprintf(buf,"/mg.%04d",(int)me);
    strcat(filename,buf);
    if (Read_OpenMGFile (filename))                                                                         {nparfiles = -1;}
    else
    if (Read_MG_General(&mg_general))                                                               {CloseMGFile (); nparfiles = -1;}

  }
}
nparfiles = UG_GlobalMinINT(nparfiles);
#endif
  if (nparfiles == -1) return(NULL);

  if (procs>nparfiles)
  {
    Broadcast(&mg_general,sizeof(MGIO_MG_GENERAL));
    if (me < nparfiles)
      mg_general.me = me;
  }
  if (me >= nparfiles)
  {
    /* BVP and format */
    if (BVPName==NULL) strcpy(BndValName,mg_general.DomainName);
    else strcpy(BndValName,BVPName);
    if (MultigridName==NULL) strcpy(MGName,mg_general.MultiGridName);
    else strcpy(MGName,MultigridName);
    if (format==NULL) strcpy(FormatName,mg_general.Formatname);
    else strcpy(FormatName,format);
    if (heapSize==0) heapSize = mg_general.heapsize * 1024;

    /* create a virginenal multigrid on the BVP */
    theMG = CreateMultiGrid(MGName,BndValName,FormatName,heapSize,FALSE);
    if (theMG==NULL)                                                                                                        {UserWrite("ERROR(ugio): cannot create multigrid\n"); CloseMGFile (); return (NULL);}
    MG_MAGIC_COOKIE(theMG) = mg_general.magic_cookie;
    if (DisposeGrid(GRID_ON_LEVEL(theMG,0)))                                                        {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    for (i=0; i<mg_general.nLevel; i++)
      if (CreateNewLevel(theMG,0)==NULL)                                                              {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    return (theMG);
  }
  if (mg_general.dim!=DIM)                                                                                        {UserWrite("ERROR: wrong dimension\n");CloseMGFile (); return (NULL);}
  if (strcmp(mg_general.version,MGIO_VERSION)!=0 && force==0)                                             {UserWrite("ERROR: wrong version\n");CloseMGFile (); return (NULL);}

  /* parallel part */
  if (mg_general.nparfiles>procs)                                                                         {UserWrite("ERROR: too many processors needed\n");CloseMGFile (); return (NULL);}
  assert(mg_general.me == me);

  /* BVP and format */
  if (BVPName==NULL) strcpy(BndValName,mg_general.DomainName);
  else strcpy(BndValName,BVPName);
  if (MultigridName==NULL) strcpy(MGName,mg_general.MultiGridName);
  else strcpy(MGName,MultigridName);
  if (format==NULL) strcpy(FormatName,mg_general.Formatname);
  else strcpy(FormatName,format);
  if (heapSize==0) heapSize = mg_general.heapsize * 1024;

  /* create a virginenal multigrid on the BVP */
  theMG = CreateMultiGrid(MGName,BndValName,FormatName,heapSize,FALSE);
  if (theMG==NULL)                                                                                                        {UserWrite("ERROR(ugio): cannot create multigrid\n"); CloseMGFile (); return (NULL);}
  MG_MAGIC_COOKIE(theMG) = mg_general.magic_cookie;
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap);
  if (DisposeGrid(GRID_ON_LEVEL(theMG,0)))                                                        {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (CreateNewLevel(theMG,0)==NULL)                                                                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  theHeap = MGHEAP(theMG);
  theBVP = MG_BVP(theMG);
  if (theBVP==NULL)                                                                                                       {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc))                                                         {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general element information */
  if (Read_GE_General(&ge_general))                                                                       {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_GE_Elements(TAGS,ge_element))                                                          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read refrule information */
  if (Read_RR_General(&rr_general))                                                                   {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  rr_rules = (MGIO_RR_RULE *)GetTmpMem(theHeap,rr_general.nRules*sizeof(MGIO_RR_RULE));
  if (rr_rules==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for rr_rules\n",(int)rr_general.nRules*sizeof(MGIO_RR_RULE)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_RR_Rules(rr_general.nRules,rr_rules))                                          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general information about coarse grid */
  if (Read_CG_General(&cg_general))                                                                       {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read coarse grid points and elements */
  cg_point = (MGIO_CG_POINT *)GetTmpMem(theHeap,cg_general.nPoint*sizeof(MGIO_CG_POINT));
  if (cg_point==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_points\n",(int)cg_general.nPoint*sizeof(MGIO_CG_POINT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_CG_Points(cg_general.nPoint,cg_point))                                         {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  if (MGIO_PARFILE)
  {
    /* read mapping: node-id --> vertex-id for orphan-nodes */
    if (Bio_Read_mint(1,&non))                                                                              {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    if (Bio_Read_mint(1,&foid))                                                                     {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    vidlist = (int*)GetTmpMem(theHeap,non*sizeof(int));
    if (Bio_Read_mint(non,vidlist))                                                                 {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    for (i=0; i<non; i++) printf("LoadMG(): vidList[%d]=%d\n",i,vidlist[i]);
  }
  else
  {
    foid = 0;
    non = cg_general.nPoint;
  }

  cg_element = (MGIO_CG_ELEMENT *)GetTmpMem(theHeap,cg_general.nElement*sizeof(MGIO_CG_ELEMENT));
  if (cg_element==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_elements\n",(int)cg_general.nElement*sizeof(MGIO_CG_ELEMENT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_CG_Elements(cg_general.nElement,cg_element))                           {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general bnd information */
  if (Bio_Jump (0))                                                                                                       {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_BD_General (&bd_general))                                                                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read bnd points */
  BndPList = (BNDP**)GetTmpMem(theHeap,bd_general.nBndP*sizeof(BNDP*));
  if (BndPList==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for BndPList\n",(int)bd_general.nBndP*sizeof(BNDP*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_PBndDesc (theBVP,theHeap,bd_general.nBndP,BndPList))           {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* create and insert coarse mesh: prepare */
  theMesh.nBndP = cg_general.nBndPoint;
  theMesh.theBndPs = BndPList;
  theMesh.nInnP = cg_general.nInnerPoint;
  if (cg_general.nInnerPoint>0)
  {
    theMesh.Position = (DOUBLE**)GetTmpMem(theHeap,cg_general.nInnerPoint*sizeof(DOUBLE*));
    if (theMesh.Position==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.Position\n",(int)cg_general.nInnerPoint*sizeof(DOUBLE*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    Positions = (DOUBLE*)GetTmpMem(theHeap,MGIO_DIM*cg_general.nInnerPoint*sizeof(DOUBLE));
    if (Positions==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Positions\n",(int)MGIO_DIM*cg_general.nInnerPoint*sizeof(DOUBLE)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }
  if (MGIO_PARFILE)
  {
    theMesh.VertexLevel = (char*)GetTmpMem(theHeap,(cg_general.nBndPoint+cg_general.nInnerPoint)*sizeof(char));
    if (theMesh.VertexLevel==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for VertexLevel\n",(int)(cg_general.nBndPoint+cg_general.nInnerPoint)*sizeof(char)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }
  else
    theMesh.VertexLevel = NULL;
  theMesh.nSubDomains = theBVPDesc.nSubDomains;
  theMesh.nElements = (INT*)GetTmpMem(theHeap,(theMesh.nSubDomains+1)*sizeof(INT));
  if (theMesh.nElements==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.nElements\n",(int)(theMesh.nSubDomains+1)*sizeof(INT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  theMesh.ElementLevel = (char**)GetTmpMem(theHeap,(theMesh.nSubDomains+1)*sizeof(char*));
  if (theMesh.ElementLevel==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.ElementLevel\n",(int)(theMesh.nSubDomains+1)*sizeof(char*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<=theMesh.nSubDomains; i++)
  {
    theMesh.nElements[i] = 0;
    theMesh.ElementLevel[i] = NULL;
  }
  theMesh.nElements[1] = cg_general.nElement;
  theMesh.ElementLevel[1] = (char*)GetTmpMem(theHeap,cg_general.nElement*sizeof(char));
  if (theMesh.ElementLevel[1]==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.ElementLevel[1]\n",(int)cg_general.nElement*sizeof(char)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<cg_general.nInnerPoint; i++)
    theMesh.Position[i] = Positions+MGIO_DIM*i;
  for (i=0; i<cg_general.nInnerPoint; i++)
  {
    cgp = MGIO_CG_POINT_PS(cg_point,cg_general.nBndPoint+i);
    for (j=0; j<MGIO_DIM; j++)
      theMesh.Position[i][j] = cgp->position[j];
  }
  if (MGIO_PARFILE)
    for (i=0; i<cg_general.nBndPoint+cg_general.nInnerPoint; i++)
    {
      cgp = MGIO_CG_POINT_PS(cg_point,i);
      theMesh.VertexLevel[i] = cgp->level;
    }

  /* nb of corners of elements */
  Element_corner_uniq_subdom  = (INT*)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT));
  if (Element_corner_uniq_subdom==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_corner_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  max = 0;
  for (i=0; i<cg_general.nElement; i++)
  {
    cge = MGIO_CG_ELEMENT_PS(cg_element,i);
    Element_corner_uniq_subdom[i] = ge_element[cge->ge].nCorner;
    max = MAX(max,ge_element[cge->ge].nCorner);
    max = MAX(max,ge_element[cge->ge].nSide);
    if (MGIO_PARFILE) theMesh.ElementLevel[1][i] = cge->level;
    else theMesh.ElementLevel[1][i] = 0;
  }
  Ecusdp[1] = Element_corner_uniq_subdom;
  theMesh.Element_corners = Ecusdp;

  /* corners ids of elements */
  Element_corner_ids_uniq_subdom  = (INT**)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT*));
  if (Element_corner_ids_uniq_subdom==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_corner_ids_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  Element_corner_ids  = (INT*)GetTmpMem(theHeap,max*cg_general.nElement*sizeof(INT));
  if (Element_corner_ids==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_corner_ids\n",(int)max*cg_general.nElement*sizeof(INT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<cg_general.nElement; i++)
    Element_corner_ids_uniq_subdom[i] = Element_corner_ids+max*i;
  for (i=0; i<cg_general.nElement; i++)
  {
    cge = MGIO_CG_ELEMENT_PS(cg_element,i);
    if (MGIO_PARFILE)
      for (j=0; j<Element_corner_uniq_subdom[i]; j++)
      {
        Element_corner_ids_uniq_subdom[i][j] = vidlist[cge->cornerid[j]-foid];
        printf("LoadMultiGrid(): cg_elem=%d  cg_nid[%d]=%d foid=%d\n",i,j,cge->cornerid[j],foid);
      }
    else
      for (j=0; j<Element_corner_uniq_subdom[i]; j++)
        Element_corner_ids_uniq_subdom[i][j] = cge->cornerid[j];
  }
  Ecidusdp[1] = Element_corner_ids_uniq_subdom;
  theMesh.Element_corner_ids = Ecidusdp;

  /* nb of elements */
  Element_nb_uniq_subdom  = (INT**)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT*));
  if (Element_nb_uniq_subdom==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_nb_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  Element_nb_ids  = (INT*)GetTmpMem(theHeap,max*cg_general.nElement*sizeof(INT));
  if (Element_nb_ids==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_nb_ids\n",(int)max*cg_general.nElement*sizeof(INT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<cg_general.nElement; i++)
    Element_nb_uniq_subdom[i] = Element_nb_ids+max*i;
  for (i=0; i<cg_general.nElement; i++)
  {
    cge = MGIO_CG_ELEMENT_PS(cg_element,i);
    for (j=0; j<ge_element[cge->ge].nSide; j++)
      Element_nb_uniq_subdom[i][j] = cge->nbid[j];
  }
  Enusdp[1] = Element_nb_uniq_subdom;
  theMesh.nbElements = Enusdp;

  /* create levels */
  for (i=1; i<mg_general.nLevel; i++)
    if (CreateNewLevel(theMG,0)==NULL)                                                              {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* insert coarse mesh */
  if (InsertMesh(theMG,&theMesh))                                                                         {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (FixCoarseGrid (theMG))                                                                                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement = PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      cge = MGIO_CG_ELEMENT_PS(cg_element,ID(theElement));
      SETREFINE(theElement,0);
      SETREFINECLASS(theElement,NO_CLASS);
      SETMARK(theElement,0);
      SETMARKCLASS(theElement,NO_CLASS);
      SETSUBDOMAIN(theElement,cge->subdomain);
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      {
        ID(CORNER(theElement,j)) = cge->cornerid[j];
      }
    }

  i = MG_ELEMUSED | MG_NODEUSED | MG_EDGEUSED | MG_VERTEXUSED |  MG_VECTORUSED;
  ClearMultiGridUsedFlags(theMG,0,TOPLEVEL(theMG),i);

  /* open identification context */
  DDD_IdentifyBegin();

  /* read parinfo of coarse-grid */
  if (MGIO_PARFILE)
  {
    ActProcListPos = ProcList = (unsigned short*)GetTmpMem(theHeap,PROCLISTSIZE*sizeof(unsigned short));
    if (ProcList==NULL)     {UserWriteF("ERROR: cannot allocate %d bytes for ProcList\n",(int)PROCLISTSIZE*sizeof(int)); return (NULL);}

    cg_pinfo.proclist = ProcList;
    for (level=0; level<=TOPLEVEL(theMG); level++)
    {
      theGrid = GRID_ON_LEVEL(theMG,level);
      for (theElement = PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=ENext)
      {
        /* store succe because Evaluate_pinfo() relinks the element list */
        ENext = SUCCE(theElement);
        cge = MGIO_CG_ELEMENT_PS(cg_element,ID(theElement));
        if (Read_pinfo (TAG(theElement),&cg_pinfo)) return (NULL);
        if (Evaluate_pinfo(theGrid,theElement,&cg_pinfo)) return (NULL);
      }
    }
  }

  /* are we ready ? */
  if (mg_general.nLevel==1)
  {
    /* close identification context */
    DDD_IdentifyEnd();

    /* repair inconsistencies */
    if (MGIO_PARFILE)
      if (IO_GridCons(theMG)) return(NULL);

    ReleaseTmpMem(theHeap);
    if (CloseMGFile ())                                                                                             {DisposeMultiGrid(theMG); return (NULL);}

    /* saved */
    MG_SAVED(theMG) = 1;
    strcpy(MG_FILENAME(theMG),filename);

    return (theMG);
  }

  /* list: node-id --> node */
  nid_n = (NODE**)GetTmpMem(theHeap,non*sizeof(NODE*));
  id = foid;
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      assert(ID(theNode)-foid<non);
      nid_n[ID(theNode)-foid] = theNode;
    }

  /* list: elem-id --> elem */
  eid_e = (ELEMENT**)GetTmpMem(theHeap,cg_general.nElement*sizeof(ELEMENT*));
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      assert(ID(theElement)>=0);
      assert(ID(theElement)<cg_general.nElement);
      eid_e[ID(theElement)] = theElement;
    }

  /* read hierarchical elements */
  refinement = (MGIO_REFINEMENT*)malloc(MGIO_REFINEMENT_SIZE);
  if (refinement==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for refinement\n",(int)MGIO_REFINEMENT_SIZE); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  ProcList = (unsigned short*)malloc(PROCLISTSIZE*sizeof(unsigned short));
  if (ProcList==NULL)     {UserWriteF("ERROR: cannot allocate %d bytes for ProcList\n",(int)PROCLISTSIZE*sizeof(int)); return (NULL);}
  for (i=0; i<MAX_SONS; i++) refinement->pinfo[i].proclist = ProcList+i*ELEMPROCLISTSIZE;
  for (j=0; j<cg_general.nElement; j++)
  {
    theElement = eid_e[j];
    theGrid = GRID_ON_LEVEL(theMG,LEVEL(theElement));
    cge = MGIO_CG_ELEMENT_PS(cg_element,j);
    if (cge->nhe==0)
    {
      SETREFINE(theElement,NO_REFINEMENT);
      SETREFINECLASS(theElement,NO_CLASS);
      SETMARK(theElement,NO_REFINEMENT);
      SETMARKCLASS(theElement,NO_CLASS);
      if (LEVEL(theElement)==0) SETECLASS(theElement,RED_CLASS);
#ifndef Debug
      else assert(0);
#endif
      continue;
    }
    if (InsertLocalTree(theGrid,theElement,refinement,rr_general.RefRuleOffset))    {DisposeMultiGrid(theMG); return (NULL);}
  }

  /* close identification context */
  DDD_IdentifyEnd();

  /* repair inconsistencies */
  if (MGIO_PARFILE)
    if (IO_GridCons(theMG)) return(NULL);

  /* postprocess */
  for (i=0; i<TOPLEVEL(theMG); i++)
    for (theElement = FIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      SETMARK(theElement,0);
      SETMARKCLASS(theElement,NO_CLASS);
    }
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);

#ifdef __THREEDIM__
    if (VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),SIDEVEC))
      for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
        for (j=0; j<SIDES_OF_ELEM(theElement); j++)
        {
          theNeighbor = NBELEM(theElement,j);
          if (theNeighbor==NULL) continue;
          for (k=0; k<SIDES_OF_ELEM(theNeighbor); k++)
            if (NBELEM(theNeighbor,k)==theElement)
              break;
          if (k==SIDES_OF_ELEM(theNeighbor))                                                                      {DisposeMultiGrid(theMG); return (NULL);}
          if (DisposeDoubledSideVector (theGrid,theElement,j,theNeighbor,k))      {DisposeMultiGrid(theMG); return (NULL);}
        }
#endif

    if (GridCreateConnection(theGrid))                                                              {DisposeMultiGrid(theMG); return (NULL);}
    ClearVectorClasses(theGrid);
    ClearNextVectorClasses(theGrid);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (ECLASS(theElement)>=GREEN_CLASS)
        SeedVectorClasses(theGrid,theElement);
      if (REFINECLASS(theElement)>=GREEN_CLASS)
        SeedNextVectorClasses(theGrid,theElement);
    }
    PropagateVectorClasses(theGrid);
    PropagateNextVectorClasses(theGrid);
  }

  /* close file */
  ReleaseTmpMem(theHeap);
  if (CloseMGFile ())                                                                                                     {DisposeMultiGrid(theMG); return (NULL);}

  /* saved */
  MG_SAVED(theMG) = 1;
  strcpy(MG_FILENAME(theMG),filename);

  return (theMG);
}

#ifdef __TWODIM__
/****************************************************************************/
/*
   SaveCNomGridAndValues - Save 2d grid & data in cnom format

   SYNOPSIS:
   INT SaveCNomGridAndValues (MULTIGRID *theMG, FILE *stream, char *symname)

   PARAMETERS:
   .  theMG - pointer to multigrid structure
   .  stream - file on which data is written
   .  symname - name of data field

   DESCRIPTION:
   Is called by the CnomCommand.

   RETURN VALUE:
   INT
   .n    NULL if an error occured
   .n    else pointer to new 'MULTIGRID'
 */
/****************************************************************************/

static DOUBLE LocalCoord[2][4][2]=
{ {{ 0, 0},{1, 0},{0,1},{ 0,0}},
  {{-1,-1},{1,-1},{1,1},{-1,1}} };

INT SaveCnomGridAndValues (MULTIGRID *theMG, char *docName, char *plotprocName, char *tag)
{
  ELEMENT *theElement;
  VERTEX *theVertex;
  GRID *theGrid;
  long nv,ne,id;
  int i,j,k,n;
  double min,max,val;
  DOUBLE *CoordOfCornerPtr[8];
  FILE *stream;
  EVALUES *PlotProcInfo;

  if (theMG==NULL)
    return(0);

  if ((PlotProcInfo=GetElementValueEvalProc(plotprocName))==NULL)
  {
    PrintErrorMessage('E',"SaveCnomGridAndValues","can't find ElementValueEvalProc");
    return(1);
  }

  stream = fopen(docName,"w");
  if (stream==NULL)
  {
    PrintErrorMessage('E',"SaveCnomGridAndValues","can't open file");
    return(1);
  }

  if (PlotProcInfo->PreprocessProc!=NULL)
    if ((*PlotProcInfo->PreprocessProc)(NULL,theMG)!=0)
      return(1);

  j=TOPLEVEL(theMG);

  /* count elements and vertices */
  nv = ne = 0;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
    {
      nv++;
      SETUSED(theVertex,0);
    }
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
        ne++;
  }

  /* write header */
  fprintf(stream,">DATA\n");
  fprintf(stream,">TIME(S) 0.0\n");
  fprintf(stream,">NV: %d\n",nv);
  fprintf(stream,">NE: %d\n",ne);

  /* compute min and max */
  min = MAX_D; max = -MAX_D;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
      {
        CORNER_COORDINATES(theElement,n,CoordOfCornerPtr);
        for (i=0; i<n; i++)
        {
          val=(*PlotProcInfo->EvalProc)(theElement, (const DOUBLE **) CoordOfCornerPtr, (DOUBLE *) &(LocalCoord[n-3,i,0]));
          min = MIN(val,min);
          max = MAX(val,max);
        }
      }
  }

  fprintf(stream,">MIN\n");
  fprintf(stream," %s\n",tag);
  fprintf(stream," %15.8LE\n",min);
  fprintf(stream,">MAX\n");
  fprintf(stream," %s\n,tag");
  fprintf(stream," %15.8LE\n",max);
  fprintf(stream,">FIN\n");

  /* write x values now */
  fprintf(stream,">X\n");
  id = 0;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
        for (i=0; i<TAG(theElement); i++)
        {
          theVertex=MYVERTEX(CORNER(theElement,i));
          if (USED(theVertex)) continue;
          fprintf(stream," %15.8lE",(double)XC(theVertex));
          ID(theVertex)=id++;
          if (id%5==0) fprintf(stream,"\n");
          SETUSED(theVertex,1);
        }
  }
  if (id%5!=0) fprintf(stream,"\n");

  /* write y values now */
  fprintf(stream,">Y\n");
  id = 0;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
        for (i=0; i<TAG(theElement); i++)
        {
          theVertex=MYVERTEX(CORNER(theElement,i));
          if (USED(theVertex)==0) continue;
          fprintf(stream," %15.8lE",(double)YC(theVertex));
          id++;
          if (id%5==0) fprintf(stream,"\n");
          SETUSED(theVertex,0);
        }
  }
  if (id%5!=0) fprintf(stream,"\n");

  /* write element list now */
  fprintf(stream,">E\n");
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
      {
        if (TAG(theElement)==3)
          fprintf(stream,"%ld %ld %ld\n",(long)ID(MYVERTEX(CORNER(theElement,0))),(long)ID(MYVERTEX(CORNER(theElement,1))),(long)ID(MYVERTEX(CORNER(theElement,2))));
        else
          fprintf(stream,"%ld %ld %ld %ld\n",ID(MYVERTEX(CORNER(theElement,0))),(long)ID(MYVERTEX(CORNER(theElement,1))),(long)ID(MYVERTEX(CORNER(theElement,2))),(long)ID(MYVERTEX(CORNER(theElement,3))));
      }
  }

  /* write data field now */
  fprintf(stream,">Z\n");
  fprintf(stream," %s\n",tag);
  id = 0;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
      {
        CORNER_COORDINATES(theElement,n,CoordOfCornerPtr);
        for (i=0; i<n; i++)
        {
          theVertex=MYVERTEX(CORNER(theElement,i));
          if (USED(theVertex)) continue;
          val=(*PlotProcInfo->EvalProc)(theElement, (const DOUBLE **) CoordOfCornerPtr, (DOUBLE *) &(LocalCoord[n-3,i,0]));
          fprintf(stream," %15.8lE",val);
          id++;
          if (id%5==0) fprintf(stream,"\n");
          SETUSED(theVertex,1);
        }
      }
  }
  if (id%5!=0) fprintf(stream,"\n");

  fprintf(stream,"<\n");
  fclose(stream);

  return(0);
}
#endif

INT InitUgio ()
{
  /* read gridpaths from defaults file (iff) */
  gridpaths_set = FALSE;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"gridpaths")==0)
    gridpaths_set = TRUE;

  if (MGIO_Init ()) return (1);

  return (0);
}
