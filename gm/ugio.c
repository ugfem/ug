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
#include "ugstruct.h"

#include "devices.h"

#include "switch.h"
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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static INT IO_Get_Sons_of_ElementSide (ELEMENT *theElement, INT side, INT *Sons_of_Side, ELEMENT *SonList[MAX_SONS], INT SonSides[MAX_SONS], INT dummy);

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
  INT i,id,move,l,tl;

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
    if (BNDP_BndPDesc(V_BNDP(theVertex),&move))
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
          if (BNDP_BndPDesc(V_BNDP(theVertex),&move))
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
  INT i,id,move;

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
  fprintf(stream,HEADER_FMT,BVPD_NAME(theBVPDesc),buffer,name,comment);

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
    if (BNDP_BndPDesc(V_BNDP(theVertex),&move))
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
    if (BNDP_BndPDesc(V_BNDP(theVertex),&move))
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

static INT Write_RefRules (MULTIGRID *theMG, INT *RefRuleOffset)
{
  MGIO_RR_GENERAL rr_general;
  INT i,j,k,t,nRules;
  HEAP *theHeap;
  MGIO_RR_RULE *Refrule, *RR;
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
  RR = Refrule = (MGIO_RR_RULE *)GetTmpMem(theHeap,nRules*sizeof(MGIO_RR_RULE));
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
  if (Write_RR_Rules(nRules,RR)) return (1);

  return (0);
}

#ifdef __TWODIM__

static INT SetRefinement (ELEMENT *theElement, MGIO_REFINEMENT *refinement, INT *RefRuleOffset)
{
  REFRULE *theRule;
  NODE *theNode;
  INT i,j,n;
  int sonRefined;
  ELEMENT *SonList[MAX_SONS];

  if (REFINE(theElement)==NO_REFINEMENT) return (1);
  refinement->refrule = REFINE(theElement) + RefRuleOffset[TAG(theElement)];
  theRule = RefRules[TAG(theElement)] + REFINE(theElement);
  refinement->refclass = ECLASS(SON(theElement,0));
  GetSons(theElement,SonList);

  /* store new corner ids */
  for (i=0,n=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    theNode = SONNODE(CORNER(theElement,i));
    if (theNode==NULL) return (1);
    refinement->newcornerid[n++] = ID(theNode);
  }

  for (i=0; i<EDGES_OF_ELEM(theElement)+1; i++)
  {
    if (i<=EDGES_OF_ELEM(theElement)) j=i;
    else j=CENTER_NODE_INDEX(theElement);
    if (theRule->pattern[j]!=1) continue;
    theNode = CORNER(SonList[theRule->sonandnode[j][0]],theRule->sonandnode[j][1]);
    if (theNode==NULL) return (1);
    refinement->newcornerid[n++] = ID(theNode);
  }
  refinement->nnewcorners = n;

  /* sons are refined ? */
  sonRefined=0;
  for (i=0; SonList[i]!=NULL; i++)
    if (REFINE(SonList[i])!=NO_REFINEMENT)
      sonRefined |= (1<<i);
  refinement->sonref = sonRefined;

  /* not movable at the moment */
  refinement->nmoved = 0;

  return (0);
}

static INT SetHierRefinement (ELEMENT *theElement, MGIO_REFINEMENT *refinement, INT *RefRuleOffset)
{
  INT i;
  ELEMENT *SonList[MAX_SONS];
  static MGIO_REFINEMENT *actualPosition;

  if (refinement!=NULL) actualPosition = refinement;

  if (REFINE(theElement)==NO_REFINEMENT) return (1);
  if (SetRefinement (theElement,actualPosition,RefRuleOffset)) return (1);
  actualPosition++;
  GetSons(theElement,SonList);
  for (i=0; SonList[i]!=NULL; i++)
  {
    if (REFINE(SonList[i])==NO_REFINEMENT) continue;
    if (SetHierRefinement(SonList[i],NULL,RefRuleOffset)) return (1);
  }

  return (0);
}

static INT nHierElements (ELEMENT *theElement, INT *n)
{
  INT i;
  ELEMENT *SonList[MAX_SONS];

  if (REFINE(theElement)==NO_REFINEMENT) return (0);
  (*n)++;
  GetSons(theElement,SonList);
  for (i=0; SonList[i]!=NULL; i++)
  {
    if (REFINE(SonList[i])!=NO_REFINEMENT)
      if (nHierElements(SonList[i],n)) return (1);
  }

  return (0);
}

#endif

#ifdef __THREEDIM__

static INT SetRefinement (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS], MGIO_REFINEMENT *refinement, INT *RefRuleOffset)
{
  REFRULE *theRule;
  NODE *theNode;
  INT i,j,n,sonRefined;

  if (REFINE(theElement)==NO_REFINEMENT)
  {
    refinement->refrule = -1;
    refinement->nnewcorners = 0;
    refinement->nmoved = 0;
    refinement->refclass = 0;
    return (0);
  }
  else
    refinement->refrule = REFINE(theElement) + RefRuleOffset[TAG(theElement)];
  theRule = RefRules[TAG(theElement)] + REFINE(theElement);

  /* store new corner ids */
  for (i=0,n=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    theNode = SONNODE(CORNER(theElement,i));
    if (theNode==NULL) return (1);
    refinement->newcornerid[n++] = ID(theNode);
  }
  for (i=0; i<EDGES_OF_ELEM(theElement)+SIDES_OF_ELEM(theElement)+1; i++)
  {
    if (i<=EDGES_OF_ELEM(theElement)+SIDES_OF_ELEM(theElement)) j=i;
    else j=CENTER_NODE_INDEX(theElement);
    if (theRule->pattern[j]!=1) continue;
    theNode = CORNER(SonList[theRule->sonandnode[j][0]],theRule->sonandnode[j][1]);
    if (theNode==NULL) return (1);
    refinement->newcornerid[n++] = ID(theNode);
  }
  refinement->nnewcorners = n;

  /* sons are refined ? */
  sonRefined=0;
  for (i=0,n=0; i<NSONS(theElement); i++)
    if (REFINE(SonList[i])!=NO_REFINEMENT)
      sonRefined |= (1<<i);
  refinement->sonref = sonRefined;

  /* not movable at the moment */
  refinement->nmoved = 0;

  /* set refinement class */
  refinement->refclass = ECLASS(SonList[0]);

  return (0);
}

static INT SetHierRefinement (ELEMENT *theElement, MGIO_REFINEMENT *refinement, INT *RefRuleOffset)
{
  INT i;
  ELEMENT *theSon;
  ELEMENT *SonList[MAX_SONS];
  static MGIO_REFINEMENT *actualPosition;

  if (refinement!=NULL) actualPosition = refinement;

  if (REFINE(theElement)==NO_REFINEMENT) return (1);
  if (GetSons(theElement,SonList)) return (1);
  if (SetRefinement (theElement,SonList,actualPosition,RefRuleOffset)) return (1);
  actualPosition++;
  for (i=0; i<NSONS(theElement); i++)
  {
    theSon = SonList[i];
    if (REFINE(theSon)!=NO_REFINEMENT)
      if (SetHierRefinement(theSon,NULL,RefRuleOffset)) return (1);
  }

  return (0);
}

static INT nHierElements (ELEMENT *theElement, INT *n)
{
  INT i;
  ELEMENT *SonList[MAX_SONS];

  if (REFINE(theElement)==NO_REFINEMENT) return (0);
  (*n)++;
  if (GetSons(theElement,SonList)) return (1);
  for (i=0; i<NSONS(theElement); i++)
    if (REFINE(SonList[i])!=NO_REFINEMENT)
      if (nHierElements(SonList[i],n)) return (1);

  return (0);
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
  MGIO_CG_ELEMENT *cg_element;
  MGIO_REFINEMENT *refinement;
  MGIO_BD_GENERAL bd_general;
  INT i,j,k,niv,nbv,nie,nbe,n,nhe,hr_max,mode;
  INT RefRuleOffset[TAGS];
  char *p;
  BNDP **BndPList;
  char filename[NAMESIZE];

  /* check */
  if (theMG==NULL) return (1);
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap);

  /* open file */
  if (strcmp(type,"dbg")==0) mode = BIO_DEBUG;
  else if (strcmp(type,"asc")==0) mode = BIO_ASCII;
  else if (strcmp(type,"bin")==0) mode = BIO_BIN;
  else return (1);
  strcpy(filename,name);
  strcat(filename,".ug.mg.");
  strcat(filename,type);
  if (Write_OpenMGFile (filename)) return (1);

  /* write general information */
  theBVP = MG_BVP(theMG);
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc)) return (1);
  mg_general.mode                 = mode;
  mg_general.dim                  = DIM;
  mg_general.magic_cookie = MG_MAGIC_COOKIE(theMG);
  mg_general.heapsize             = MGHEAP(theMG)->size/1024;
  mg_general.nLevel               = TOPLEVEL(theMG) + 1;
  mg_general.nNode = mg_general.nPoint = mg_general.nElement;
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
  strcpy(mg_general.DomainName,BVPD_NAME(theBVPDesc));
  strcpy(mg_general.MultiGridName,MGNAME(theMG));
  strcpy(mg_general.Formatname,ENVITEM_NAME(MGFORMAT(theMG)));
  mg_general.VectorTypes  = 0;
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

  /* write general information about coarse grid */
  theGrid = GRID_ON_LEVEL(theMG,0);
  cg_general.nPoint = NV(theGrid);
  niv=nbv=0;
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    if (OBJT(MYVERTEX(theNode))==IVOBJ) niv++;
    else nbv++;
  cg_general.nBndPoint = nbv;
  cg_general.nInnerPoint = niv;
  cg_general.nElement = NT(theGrid);
  nie=nbe=0;
  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (OBJT(theElement)==IEOBJ) nie++;
    else if (OBJT(theElement)==BEOBJ) nbe++;
    else return (1);
  }
  cg_general.nBndElement = nbe;
  cg_general.nInnerElement = nie;
  if (Write_CG_General(&cg_general)) return (1);

  /* write coarse grid points */
  if (RenumberNodeElem (theMG)) return (1);
  theGrid = GRID_ON_LEVEL(theMG,0);
  n = NV(theGrid)*sizeof(MGIO_CG_POINT);
  cg_point = (MGIO_CG_POINT *)GetTmpMem(theHeap,n);
  if (cg_point==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_points\n",(int)n); return (1);}
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    for (j=0; j<MGIO_DIM; j++)
      cg_point[ID(theNode)].position[j] = CVECT(MYVERTEX(theNode))[j];

  if (Write_CG_Points((int)NV(theGrid),cg_point)) return (1);

  /* write coarse grid elements */
  theGrid = GRID_ON_LEVEL(theMG,0);
  n = NT(theGrid); hr_max=0;
  cg_element = (MGIO_CG_ELEMENT*)GetTmpMem(theHeap,n*sizeof(MGIO_CG_ELEMENT));
  if (cg_element==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_elements\n",(int)n*sizeof(MGIO_CG_ELEMENT)); return (1);}
  for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    i = ID(theElement);

    /* coarse grid part */
    cg_element[i].ge = TAG(theElement);
    nhe=0;
    if (nHierElements (theElement,&nhe)) return (1);
    hr_max = MAX(hr_max,nhe);
    cg_element[i].nhe = nhe;
    for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      cg_element[i].cornerid[j] = ID(CORNER(theElement,j));
    for (j=0; j<SIDES_OF_ELEM(theElement); j++)
      if (NBELEM(theElement,j)!=NULL)
        cg_element[i].nbid[j] = ID(NBELEM(theElement,j));
      else
        cg_element[i].nbid[j] = -1;
    cg_element[i].subdomain = SUBDOMAIN(theElement);
  }
  if (Write_CG_Elements((int)n,cg_element)) return (1);

  /* write general bnd information */
  if (Bio_Jump_From ()) return (1);
  bd_general.nBndP = nbv;
  if (Write_BD_General (&bd_general)) return (1);

  /* write bnd information */
  BndPList = (BNDP**)GetTmpMem(theHeap,nbv*sizeof(BNDP*));
  if (BndPList==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for BndPList\n",(int)nbv*sizeof(BNDP*)); return (1);}
  for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,0)); theNode!=NULL; theNode=SUCCN(theNode))
    if (OBJT(MYVERTEX(theNode))==BVOBJ)
    {
      if (ID(theNode) < 0 || ID(theNode) >= nbv) return (1);
      BndPList[ID(theNode)] = V_BNDP(MYVERTEX(theNode));
    }
  if (Write_PBndDesc (nbv,BndPList)) return (1);
  if (Bio_Jump_To ()) return (1);

  /* save refinement */
  refinement = (MGIO_REFINEMENT *)GetTmpMem(theHeap,(hr_max+1)*sizeof(MGIO_REFINEMENT));
  if (refinement==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for refinement\n",(int)hr_max*sizeof(MGIO_REFINEMENT)); return (1);}
  theGrid = GRID_ON_LEVEL(theMG,0);
  for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    i = ID(theElement);
    if (cg_element[i].nhe==0) continue;
    if (SetHierRefinement (theElement,refinement,RefRuleOffset)) return (1);
    if (Write_Refinement (cg_element[i].nhe,refinement)) return (1);
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

#ifdef __THREEDIM__

#define MGIO_PATHDEPTHMASK 0xF0000000
#define MGIO_PATHDEPTHSHIFT 28
#define MGIO_PATHDEPTH(i)                               (((i) & MGIO_PATHDEPTHMASK)>>MGIO_PATHDEPTHSHIFT)

#define MGIO_NEXTSIDEMASK 0x00000003
#define MGIO_NEXTSIDE(i,n)                              (((i) & (MGIO_NEXTSIDEMASK<<(2*(n))))>>(2*(n)))

#define MGIO_NEXTSIDEMASKHEX 0x00000007
#define MGIO_NEXTSIDEHEX(i,n)                                                   (((i) & (MGIO_NEXTSIDEMASKHEX<<(3*(n))))>>(3*(n)))

static INT IO_GetSons (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS])
{
  MGIO_RR_RULE *theRule;
  ELEMENT *theSon;
  struct mgio_sondata *sonData;
  INT SonID,PathPos;

  if (theElement==NULL) return(1);
  theRule = rr_rules+REFINE(theElement);
  for (SonID=0; SonID<MAX_SONS; SonID++) SonList[SonID] = NULL;
  if (NSONS(theElement) == 0) return(0);
  SonList[0] = SON(theElement,0);
  switch (TAG(theElement))
  {
  case (TETRAHEDRON) :
    for (SonID=1; SonID<theRule->nsons; SonID++)
    {
      sonData = &(theRule->sons[SonID]);
      theSon = SonList[0];
      for (PathPos=0; PathPos<MGIO_PATHDEPTH(sonData->path); PathPos++)
        theSon = NBELEM(theSon,MGIO_NEXTSIDE(sonData->path,PathPos));
      if (theSon==NULL) return (1);
      SonList[SonID] = theSon;
    }
    break;

  case (HEXAHEDRON) :
    for (SonID=1; SonID<theRule->nsons; SonID++)
    {
      sonData = &(theRule->sons[SonID]);
      theSon = SonList[0];
      for (PathPos=0; PathPos<MGIO_PATHDEPTH(sonData->path); PathPos++)
        theSon = NBELEM(theSon,MGIO_NEXTSIDEHEX(sonData->path,PathPos));
      if (theSon==NULL) return(1);
      SonList[SonID] = theSon;
    }
    break;
  default :
    return(1);
  }

  return(0);
}

static INT IO_GetSideNode (ELEMENT *theElement, INT side, NODE **theNode, INT *nb_side)
{
  INT nbside,sidenode_index;
  ELEMENT *theNeighbor;
  MGIO_RR_RULE *nbrefrule;
  ELEMENT *SonList[MAX_SONS];

  theNeighbor = NBELEM(theElement,side);
  if (theNeighbor==NULL || SON(theNeighbor,0)==NULL) {*theNode=NULL; return (0);}
  for (nbside=0; nbside<SIDES_OF_ELEM(theNeighbor); nbside++)
    if (NBELEM(theNeighbor,nbside)==theElement)
      break;
  if (nbside==SIDES_OF_ELEM(theNeighbor)) return (1);
  *nb_side = nbside;
  nbrefrule = rr_rules+REFINE(theNeighbor);
  if (IO_GetSons (theNeighbor,SonList)) return (1);
  sidenode_index = EDGES_OF_ELEM(theNeighbor)+nbside;
  *theNode = CORNER(SonList[nbrefrule->sonandnode[sidenode_index][0]],nbrefrule->sonandnode[sidenode_index][1]);
  if (*theNode==NULL) return (1);

  return (0);
}

#endif

static INT IO_Get_Sons_of_ElementSide (ELEMENT *theElement, INT side, INT *Sons_of_Side, ELEMENT *SonList[MAX_SONS], INT SonSides[MAX_SONS], INT dummy)
{
  INT i,j;
  MGIO_RR_RULE *theRule;
  struct mgio_sondata *son;

  /* init */
  *Sons_of_Side = 0;
  theRule = rr_rules + REFINE(theElement);

#ifdef __TWODIM__
  GetSons (theElement,SonList);
#endif
#ifdef __THREEDIM__
  if (IO_GetSons(theElement,SonList) != GM_OK) RETURN(GM_FATAL);
#endif

  for (i=0; i<theRule->nsons; i++)
  {
    son = &(theRule->sons[i]);
    for (j=0; j<SIDES_OF_TAG(son->tag); j++)
      if (son->nb[j]==20+side)
      {
        SonList[*Sons_of_Side] = SonList[i];
        SonSides[*Sons_of_Side] = j;
        (*Sons_of_Side)++;
      }
  }

  /* fill remaining with zero */
  for (i=*Sons_of_Side; i<MAX_SONS; i++) SonList[i] = NULL;

  return(GM_OK);
}

static INT InsertLocalTree (GRID *theGrid, ELEMENT *theElement, MGIO_REFINEMENT *refinement)
{
  INT i,j,l_entry,r_index,nedge,type,sonRefined,n0,n1,Sons_of_Side,Sons_of_Side_List[MAX_SONS];
  ELEMENT *theSonElem[MAX_SONS];
  ELEMENT *SonList[MAX_SONS];
  NODE *NodeList[MAX_NEW_CORNERS_DIM+MAX_CORNERS_OF_ELEM];
  NODE *SonNodeList[MAX_CORNERS_OF_ELEM];
  GRID *upGrid;
  EDGE *theEdge;
  MGIO_RR_RULE *theRule;
  static MGIO_REFINEMENT *ref;
  struct mgio_sondata *SonData;
#       ifdef __THREEDIM__
  INT nbside;
#       endif

  /* init */
  if (refinement!=NULL) ref=refinement;
  if (ref->refrule==-1) return (1);
  SETREFINE(theElement,ref->refrule);
  SETREFINECLASS(theElement,RED_CLASS);
  SETMARK(theElement,ref->refrule);
  SETMARKCLASS(theElement,RED_CLASS);
  theRule = rr_rules+ref->refrule;
  upGrid = UPGRID(theGrid);

  /* insert nodes */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    NodeList[i] = SONNODE(CORNER(theElement,i));
    if (NodeList[i]==NULL)
    {
      NodeList[i] = CreateSonNode(upGrid,CORNER(theElement,i));
      if (NodeList[i]==NULL) return (1);
      ID(NodeList[i]) = ref->newcornerid[i];
    }
  }

  nedge = EDGES_OF_ELEM(theElement);
  l_entry = r_index = CORNERS_OF_ELEM(theElement);
  for (i=0; i<nedge; i++)
  {
    if (theRule->pattern[i]!=1)
    {
      NodeList[l_entry++] = NULL;
      continue;
    }
    n0 = CORNER_OF_EDGE(theElement,i,0);
    n1 = CORNER_OF_EDGE(theElement,i,1);
    theEdge = GetEdge(CORNER(theElement,n0),CORNER(theElement,n1));
    if (theEdge==NULL) return (1);
    NodeList[l_entry] = MIDNODE(theEdge);
    if (NodeList[l_entry]==NULL)
    {
      NodeList[l_entry] = CreateMidNode(upGrid,theElement,i);
      if (NodeList[l_entry]==NULL) return (1);
      ID(NodeList[l_entry]) = ref->newcornerid[r_index];
    }
    l_entry++; r_index++;
  }

#ifdef __THREEDIM__
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    if (theRule->pattern[i+nedge]!=1)
    {
      NodeList[l_entry++] = NULL;
      continue;
    }
    if (IO_GetSideNode(theElement,i,&(NodeList[l_entry]),&nbside)) return (1);
    if (NodeList[l_entry]==NULL)
    {
      NodeList[l_entry] = CreateSideNode(upGrid,theElement,i);
      if (NodeList[l_entry]==NULL) return (1);
      ID(NodeList[l_entry]) = ref->newcornerid[r_index];
    }
    else
      SETONNBSIDE(MYVERTEX(NodeList[l_entry]),nbside);
    l_entry++; r_index++;
  }
#endif

  if (theRule->pattern[CENTER_NODE_INDEX(theElement)]==1)
  {
    NodeList[l_entry] = CreateCenterNode(upGrid,theElement);
    if (NodeList[l_entry]==NULL) return (1);
    ID(NodeList[l_entry]) = ref->newcornerid[r_index];
    l_entry++;
  }

  /* insert elements */
  for (i=0; i<theRule->nsons; i++)
  {
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
    theSonElem[i] = CreateElement(upGrid,SonData->tag,type,SonNodeList,theElement);
    if (theSonElem[i]==NULL) return (1);
    SETECLASS(theSonElem[i],ref->refclass);
    SETSUBDOMAIN(theSonElem[i],SUBDOMAIN(theElement));
  }

  /* set neighbor relation between sons */
  for (i=0; i<theRule->nsons; i++)
  {
    SonData = theRule->sons+i;
    for (j=0; j<SIDES_OF_TAG(SonData->tag); j++)
      if (SonData->nb[j]<20)
        SET_NBELEM(theSonElem[i],j,theSonElem[SonData->nb[j]]);
      else
        SET_NBELEM(theSonElem[i],j,NULL);
  }

  /* connect to neighbors */
  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    if (IO_Get_Sons_of_ElementSide(theElement,i,&Sons_of_Side,SonList,Sons_of_Side_List,1)) return (1);
    if (Connect_Sons_of_ElementSide(UPGRID(theGrid),theElement,i,Sons_of_Side,SonList,Sons_of_Side_List)) return (1);
  }

  /* jump to the sons ? */
  sonRefined = ref->sonref;

  /* one refinement used */
  ref++;

  /* call recoursively */
  for (i=0; i<theRule->nsons; i++)
    if (sonRefined & (1<<i))
    {
      if (InsertLocalTree (upGrid,theSonElem[i],NULL)) return (1);
    }
    else
    {
      SETREFINE(theSonElem[i],NO_REFINEMENT);
    }

  return (0);
}

MULTIGRID *LoadMultiGrid (char *MultigridName, char *name, char *type, char *BVPName, char *format, unsigned long heapSize, DOUBLE_VECTOR global0, DOUBLE_VECTOR global1, DOUBLE_VECTOR global2)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  ELEMENT *theElement;
  HEAP *theHeap;
  MGIO_MG_GENERAL mg_general;
  MGIO_GE_GENERAL ge_general;
  MGIO_GE_ELEMENT ge_element[TAGS];
  MGIO_RR_GENERAL rr_general;
  MGIO_CG_GENERAL cg_general;
  MGIO_CG_POINT *cg_point;
  MGIO_CG_ELEMENT *cg_element;
  MGIO_BD_GENERAL bd_general;
  MGIO_REFINEMENT *refinement;
  BNDP **BndPList;
  DOUBLE *Positions;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  MESH theMesh;
  char FormatName[NAMESIZE], BndValName[NAMESIZE], MGName[NAMESIZE], filename[NAMESIZE];
  INT i,j,*Element_corner_uniq_subdom, *Ecusdp[2],**Enusdp[2],**Ecidusdp[2],**Element_corner_ids_uniq_subdom,*Element_corner_ids,max,**Element_nb_uniq_subdom,*Element_nb_ids;
#       ifdef __THREEDIM__
  INT k;
  ELEMENT *theNeighbor;
  ELEMENT *e, *nbe;
  ELEMENT *e0, *e1, *e2;
  FIFO myfifo;
  void *buffer;
  INT sid;
  INT n;
#       endif

  /* open file */
  strcpy(filename,name);
  strcat(filename,".ug.mg.");
  strcat(filename,type);
  if (Read_OpenMGFile (filename))                                                                         {return (NULL);}
  if (Read_MG_General(&mg_general))                                                                       {CloseMGFile (); return (NULL);}
  if (mg_general.dim!=DIM)                                                                                        {UserWrite("ERROR: wrong dimension\n");CloseMGFile (); return (NULL);}
  if (strcmp(mg_general.version,MGIO_VERSION)!=0)                                         {UserWrite("ERROR: wrong version\n");CloseMGFile (); return (NULL);}

  /* BVP and format */
  if (BVPName==NULL) strcpy(BndValName,mg_general.DomainName);
  else strcpy(BndValName,BVPName);
  if (MultigridName==NULL) strcpy(MGName,mg_general.MultiGridName);
  else strcpy(MGName,MultigridName);
  if (format==NULL) strcpy(FormatName,mg_general.Formatname);
  else strcpy(FormatName,format);
  if (heapSize==0) heapSize = mg_general.heapsize * 1024;

  /* create a virginenal multigrid on the BVP */
  theMG = CreateMultiGrid(MGName,BndValName,FormatName,heapSize);
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
  cg_element = (MGIO_CG_ELEMENT *)GetTmpMem(theHeap,cg_general.nElement*sizeof(MGIO_CG_ELEMENT));
  if (cg_element==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_elements\n",(int)cg_general.nElement*sizeof(MGIO_CG_ELEMENT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_CG_Points(cg_general.nPoint,cg_point))                                         {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
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
  for (i=0; i<cg_general.nInnerPoint; i++)
    theMesh.Position[i] = Positions+MGIO_DIM*i;
  for (i=0; i<cg_general.nInnerPoint; i++)
    for (j=0; j<MGIO_DIM; j++)
      theMesh.Position[i][j] = cg_point[cg_general.nBndPoint+i].position[j];
  theMesh.nSubDomains = theBVPDesc.nSubDomains;
  theMesh.nElements = (INT*)GetTmpMem(theHeap,(theMesh.nSubDomains+1)*sizeof(INT));
  if (theMesh.nElements==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.nElements\n",(int)(theMesh.nSubDomains+1)*sizeof(INT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<=theMesh.nSubDomains; i++) theMesh.nElements[i] = 0;theMesh.nElements[1] = cg_general.nElement;

  /* nb of corners of elements */
  Element_corner_uniq_subdom  = (INT*)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT));
  if (Element_corner_uniq_subdom==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_corner_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  max = 0;
  for (i=0; i<cg_general.nElement; i++)
  {
    Element_corner_uniq_subdom[i] = ge_element[cg_element[i].ge].nCorner;
    max = MAX(max,ge_element[cg_element[i].ge].nCorner);
    max = MAX(max,ge_element[cg_element[i].ge].nSide);
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
    for (j=0; j<Element_corner_uniq_subdom[i]; j++)
      Element_corner_ids_uniq_subdom[i][j] = cg_element[i].cornerid[j];
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
    for (j=0; j<ge_element[cg_element[i].ge].nSide; j++)
      Element_nb_uniq_subdom[i][j] = cg_element[i].nbid[j];
  Enusdp[1] = Element_nb_uniq_subdom;
  theMesh.nbElements = Enusdp;

  /* insert coarse mesh */
  if (InsertMesh(theMG,&theMesh))                                                                         {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  theGrid = GRID_ON_LEVEL(theMG,0);
  for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    SETREFINE(theElement,0);
    SETREFINECLASS(theElement,NO_CLASS);
    SETMARK(theElement,0);
    SETMARKCLASS(theElement,NO_CLASS);
    SETSUBDOMAIN(theElement,cg_element[ID(theElement)].subdomain);
  }

  /* are we ready ? */
  if (mg_general.nLevel==1)
  {
    ReleaseTmpMem(theHeap);
    if (CloseMGFile ())                                                                                             {DisposeMultiGrid(theMG); return (NULL);}

    /* saved */
    MG_SAVED(theMG) = 1;
    strcpy(MG_FILENAME(theMG),filename);

    return (theMG);
  }

  /* create grids */
  for (i=1; i<mg_general.nLevel; i++)
    if (CreateNewLevel(theMG,0)==NULL)                                                              {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read hierarchical elements */
  theGrid = GRID_ON_LEVEL(theMG,0); max=0;
  for (i=0; i<cg_general.nElement; i++) max = MAX(max,cg_element[i].nhe);
  refinement = (MGIO_REFINEMENT*)GetTmpMem(theHeap,max*sizeof(MGIO_REFINEMENT));
  if (refinement==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for refinement\n",(int)max*sizeof(MGIO_REFINEMENT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Set_Get_Sons_of_ElementSideProc(IO_Get_Sons_of_ElementSide))        {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    i = ID(theElement);
    if (cg_element[i].nhe==0) continue;
    if (Read_Refinement (cg_element[i].nhe,refinement))                     {DisposeMultiGrid(theMG); return (NULL);}
    if (InsertLocalTree(theGrid,theElement,refinement))                             {DisposeMultiGrid(theMG); return (NULL);}
  }

  /* postprocess */
  for (i=0; i<TOPLEVEL(theMG); i++)
    for (theElement = FIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (REFINE(theElement)!=NO_REFINEMENT)
        SETREFINE(theElement,REFINE(theElement)-rr_general.RefRuleOffset[TAG(theElement)]);
      SETMARK(theElement,0);
      SETMARKCLASS(theElement,NO_CLASS);
    }
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);

#ifdef __THREEDIM__
    if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
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
