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

#define EXTRACT_RULES           1       /* 1: use er.c, 0: use old version			*/

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
#ifdef ModelP
#include "parallel.h"
#endif
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

#include "er.h"

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
/*#define PROCLISTSIZE      (ELEMPROCLISTSIZE*MAX_SONS * 2) size not enough and change to static variable; Christian Wrobel 980128 */
#define PROCLISTSIZE_VALUE      (ELEMPROCLISTSIZE*MAX_SONS * MAX(5,(int)(2.0+log((double)procs))))
#define PROCLISTSIZE            proc_list_size

/* orphan condition for elements */
#define EORPHAN(e)              (EFATHER(e)==NULL || THEFLAG(e))

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

REP_ERR_FILE;
static INT gridpaths_set=FALSE;
static MGIO_RR_RULE *rr_rules;
INT rr_rule_offsets[TAGS];
static unsigned short *ProcList, *ActProcListPos;
static int foid,non;
static NODE **nid_n;                                            /* mapping: orphan-node-id  -->  orphan node */
static NODE **vid_n;                                            /* mapping: orphan-vertex-id  -->  lowest orphan node */
static INT nov;                                                         /* number of orphan vertices */
static ELEMENT **eid_e;
static int nparfiles;                                           /* nb of parallel files		*/
static int proc_list_size = -1;                         /* hold the computed value for PROCLISTSIZE; initialized with crazy dummy */

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static INT WriteElementParInfo (GRID *theGrid, ELEMENT *theElement, MGIO_PARINFO *pinfo);

/****************************************************************************/
/*D
   RenumberMultiGrid - Recalculate ids in the current order

   SYNOPSIS:
   INT RenumberMultiGrid (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - structure to renumber

   DESCRIPTION:
   This function recalculates ids in the current order.

   RETURN VALUE:
   INT
   .n   0 if ok
   D*/
/****************************************************************************/

/************************************************************************************/
/*                                                                                                                                                              */
/*    Enumeration of nodes from left to right corresponding to 0...n-1                          */
/*    ----------------------------------------------------------------				*/
/*                                                                                                                                                              */
/*                                         orphan nodes												*/
/*                                      |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|							*/
/*   N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N                      */
/*  |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|             */
/*                      ghost nodes						  master nodes						*/
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/*   Mapping should exist for all orphan nodes to vertices. 'foid' is id of first       */
/*   orphan node, 'non' number of orphan nodes. (MGIO_PAR) it looks (consequently): */
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/*                                                                      orphan nodes								*/
/*                                                      |~~~~~~~~~~~~~~~~~~~~|							*/
/*                                   N N N N N N N N N N N N N N N N N N N                      */
/*                                                          |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|            */
/*                                                                                        master nodes						*/
/*                                                                                                                                                              */
/*                                                                                                                                                              */
/************************************************************************************/

static INT RenumberNodes (MULTIGRID *theMG, INT *foid, INT *non)
{
  INT i,nid;
  NODE *theNode;

  nid=0;
  if (procs==1)
  {
    /* ids for all nodes */
    for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,0)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      ID(theNode) = ID(MYVERTEX(theNode));
      nid = MAX(nid,ID(theNode));
    }
    nid++;
    non[0] = nid;
    for (i=1; i<=TOPLEVEL(theMG); i++)
      for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        ID(theNode) = nid++;
    foid[0] = 0;
  }
  else
  {
    /* ids for other nodes */
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (GHOST(theNode) && !USED(theNode))
          ID(theNode) = nid++;
    foid[0] = nid;

    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (GHOST(theNode) && USED(theNode))
          ID(theNode) = nid++;

    /* ids for master nodes */
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (MASTER(theNode) && USED(theNode))
          ID(theNode) = nid++;
    non[0] = nid-foid[0];

    /* ids for master nodes */
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (MASTER(theNode) && !USED(theNode))
          ID(theNode) = nid++;

  }

  return (0);
}

INT RenumberMultiGrid (MULTIGRID *theMG, INT *nboe, INT *nioe, INT *nbov, INT *niov, NODE ***vid_n, INT *foid, INT *non, INT MarkKey)
{
  NODE *theNode;
  ELEMENT *theElement;
  INT i,n_ioe,n_boe,vid,n_iov,n_bov,eid,lfoid,lnon,j;

  /* init used-flags */
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      SETUSED(theNode,0);
      SETUSED(MYVERTEX(theNode),0);
      SETTHEFLAG(MYVERTEX(theNode),0);
    }

  /* renumber elements and set orphan-flags for nodes and vertices */
  eid=n_ioe=n_boe=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
      if (EFATHER(theElement)==NULL || THEFLAG(theElement)==1)
      {
        ID(theElement) = eid++;
        if (OBJT(theElement)==BEOBJ) n_boe++;
        else n_ioe++;
        for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
        {
          SETUSED(CORNER(theElement,j),1);
          SETUSED(MYVERTEX(CORNER(theElement,j)),1);
        }
                                #ifdef ModelP
        assert(i==0 || EGHOST(theElement));
                                #endif
      }

  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
      if (EFATHER(theElement)!=NULL && THEFLAG(theElement)==0)
        ID(theElement) = eid++;
  if (nboe!=NULL) *nboe = n_boe;
  if (nioe!=NULL) *nioe = n_ioe;

  /* renumber vertices */
  vid=n_iov=n_bov=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
      if (THEFLAG(MYVERTEX(theNode))==0 && USED(MYVERTEX(theNode)) && OBJT(MYVERTEX(theNode))==BVOBJ)
      {
        ID(MYVERTEX(theNode)) = vid++;
        n_bov++;
        SETTHEFLAG(MYVERTEX(theNode),1);
      }
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
      if (THEFLAG(MYVERTEX(theNode))==0 && USED(MYVERTEX(theNode)) && OBJT(MYVERTEX(theNode))==IVOBJ)
      {
        ID(MYVERTEX(theNode)) = vid++;
        n_iov++;
        SETTHEFLAG(MYVERTEX(theNode),1);
      }
  if (vid_n!=NULL)
  {
    *vid_n = (NODE**)GetTmpMem(MGHEAP(theMG),(n_iov+n_bov)*sizeof(NODE*),MarkKey);
    for (i=0; i<n_iov+n_bov; i++) (*vid_n)[i] = NULL;
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (USED(theNode))
        {
          assert(ID(MYVERTEX(theNode))<n_iov+n_bov);
          if ((*vid_n)[ID(MYVERTEX(theNode))]!=NULL) continue;
          (*vid_n)[ID(MYVERTEX(theNode))] = theNode;
        }
    IFDEBUG(gm,4)
    for (i=0; i<n_iov+n_bov; i++)
      assert((*vid_n)[i] != NULL);
    ENDDEBUG
  }
  for (i=0; i<=TOPLEVEL(theMG); i++)                                            /* not neccessary for i/o */
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
      if (THEFLAG(MYVERTEX(theNode))==0 && !USED(MYVERTEX(theNode)))
      {
        ID(MYVERTEX(theNode)) = vid++;
        SETTHEFLAG(MYVERTEX(theNode),1);
      }
  if (nbov!=NULL) *nbov = n_bov;
  if (niov!=NULL) *niov = n_iov;

  /* renumber nodes */
  if (RenumberNodes(theMG,&lfoid,&lnon)) return (1);
  if (foid!=NULL) *foid = lfoid;
  if (non!=NULL) *non = lnon;

  return (0);
}

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

  assert(0);    /* this function is corrupted (says Chr. Wieners) and should never be used */
  /* TODO: remove this function */

  /* set vector classes */
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if (ClearVectorClasses(theGrid)) REP_ERR_RETURN(1);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (ECLASS(theElement)!=RED_CLASS && ECLASS(theElement)!=GREEN_CLASS) continue;
      if (SeedVectorClasses(theGrid,theElement)) REP_ERR_RETURN(1);
    }
    if (PropagateVectorClasses(theGrid)) REP_ERR_RETURN(1);
  }

  /* set NextVectorClasses */
  if (ClearNextVectorClasses(GRID_ON_LEVEL(theMG,TOPLEVEL(theMG)))) REP_ERR_RETURN(1);
  for (i=0; i<TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if (ClearNextVectorClasses(theGrid)) REP_ERR_RETURN(1);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (NSONS(theElement)==0) continue;
      if (ECLASS(SON(theElement,0))!=RED_CLASS && ECLASS(SON(theElement,0))!=GREEN_CLASS) continue;
      if (SeedNextVectorClasses(theGrid,theElement)) REP_ERR_RETURN(1);
    }
    if (PropagateNextVectorClasses(theGrid)) REP_ERR_RETURN(1);
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

static INT OrphanCons(MULTIGRID *theMG)
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
      SETTHEFLAG(theElement,0);
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
              PRINTDEBUG(gm,1,(PFMT "OrphanCons(): ERROR: elem=" EID_FMTX " cornernode=" ID_FMTX " is orphannode\n",
                               me,EID_PRTX(theElement),ID_PRTX(theNode)));
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
              PRINTDEBUG(gm,1,(PFMT "OrphanCons(): ERROR: elem=" EID_FMTX " midnode=" ID_FMTX " is orphannode\n",
                               me,EID_PRTX(theElement),ID_PRTX(theNode)));
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

          /* this is an orphan */
          SETTHEFLAG(theElement,1);
          if (0)
          {
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
          }
          PRINTDEBUG(gm,1,(PFMT "OrphanCons(): new orphan elem=" EID_FMTX "\n",
                           me,EID_PRTX(theElement)));
        }
      }
    }
  }

  return(error);
}

static INT Write_RefRules (MULTIGRID *theMG, INT *RefRuleOffset, INT MarkKey)
{
  MGIO_RR_GENERAL rr_general;
  INT i,j,k,t,nRules;
  HEAP *theHeap;
  MGIO_RR_RULE *Refrule;
  REFRULE * ug_refrule;
  struct mgio_sondata *sonData;

  /* init */
  if (theMG==NULL) REP_ERR_RETURN(1);
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
  if (Write_RR_General(&rr_general)) REP_ERR_RETURN(1);

  /* write refrules */
  rr_rules = Refrule = (MGIO_RR_RULE *)GetTmpMem(theHeap,nRules*sizeof(MGIO_RR_RULE),MarkKey);
  if (Refrule==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for refrules\n",(int)nRules*sizeof(MGIO_RR_RULE)); REP_ERR_RETURN(1);}
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
  if (Write_RR_Rules(nRules,rr_rules)) REP_ERR_RETURN(1);

  return (0);
}

static INT SetRefinement (GRID *theGrid, ELEMENT *theElement,
                          NODE **NodeContext, ELEMENT *SonList[MAX_SONS],
                          INT nmax, MGIO_REFINEMENT *refinement,
                          INT *RefRuleOffset)
{
  MGIO_RR_RULE *theRule;
  INT i,j,n,sonRefined,sonex,nex,refined;
  /*ELEMENT *SonSonList[MAX_SONS];*/

  refinement->refclass = REFINECLASS(theElement);
  refinement->refrule = REFINE(theElement) + RefRuleOffset[TAG(theElement)];
  theRule = rr_rules + rr_rule_offsets[TAG(theElement)] + REFINE(theElement);

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
  i = CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement);
  if (NodeContext[i]!=NULL && ((nex>>i)&0x1))
    refinement->newcornerid[n++] = ID(NodeContext[i]);
  refinement->nnewcorners = n;

  /* sons are refined ? */
  sonRefined=sonex=refinement->nbid_ex=0;
  for (i=0,n=0; i<nmax; i++)
  {
    if (SonList[i]==NULL) continue;

    /*
            GetAllSons(SonList[i],SonSonList);

            refined = 0;
            for (j=0; SonSonList[j]!=NULL; j++)
                    if (!EORPHAN(SonSonList[j]))	refined = 1;
            if (REFINE(SonList[i])!=NO_REFINEMENT && refined)		sonRefined |= (1<<i);
     */

    if (REFINE(SonList[i])!=NO_REFINEMENT) sonRefined |= (1<<i);
    if (MGIO_PARFILE)
    {
      sonex |= (1<<i);
      if (WriteElementParInfo(theGrid,SonList[i],&refinement->pinfo[i]))
        REP_ERR_RETURN(1);
      for (j=0; j<SIDES_OF_ELEM(SonList[i]); j++)
        if (NBELEM(SonList[i],j)!=NULL && EORPHAN(NBELEM(SonList[i],j)))
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

  /* write orphan-sonsonson-node info */
  if (MGIO_PARFILE)
  {
    refinement->orphanid_ex = n = 0;
    for (i=0; i<CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement)+SIDES_OF_ELEM(theElement); i++)
    {
      if (NodeContext[i]!=NULL && ((nex>>i)&0x1))
      {
        if (ID(MYVERTEX(NodeContext[i]))<nov
            && LEVEL(NodeContext[i])<LEVEL(vid_n[ID(MYVERTEX(NodeContext[i]))]))
        {
          /* TODO remove the next assert; Christian Wrobel 980127
             it seems that the assert fails for valid load balancings
             consider the following situation:

             level 3:     |ooo*                     * is orphan node
             level 2:         &------|              & is corner node and i points to it; thus i<=CORNERS_OF_ELEM(theElement) is valid
             level 1:         |-------------|
             level 0: |---------------------|       this is also an orphan element but this doesn't influence the current case
           */
          /*assert(i>=CORNERS_OF_ELEM(theElement)); */

          refinement->orphanid[n] = ID(vid_n[ID(MYVERTEX(NodeContext[i]))]);
          refinement->orphanid_ex = 1;
        }
        else
        {
          refinement->orphanid[n] = -1;
        }
        n++;
      }
    }
    i = CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement);
    if (NodeContext[i]!=NULL && ((nex>>i)&0x1))
    {
      if (ID(MYVERTEX(NodeContext[i]))<nov
          && LEVEL(NodeContext[i])<LEVEL(vid_n[ID(MYVERTEX(NodeContext[i]))]))
      {
        refinement->orphanid[n++] = ID(vid_n[ID(MYVERTEX(NodeContext[i]))]);
        refinement->orphanid_ex = 1;
      }
      else
      {
        refinement->orphanid[n++] = -1;
      }
    }
  }

#if (MGIO_DEBUG>0)
  /* write mykey */
  refinement->mykey = KeyForObject((KEY_OBJECT *)theElement);

  /* write myfatherkey */
  if ( THEFLAG(theElement) )
    refinement->myfatherkey = -1;               /* orphan element has no father */
  else
    refinement->myfatherkey = KeyForObject((KEY_OBJECT *)EFATHER(theElement));

  /* write nbkey[] */
  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
    refinement->nbkey[j] = KeyForObject((KEY_OBJECT *)NBELEM(theElement,j));

  n=0;

  refinement->mycorners = CORNERS_OF_ELEM(theElement);

  /* write mycornerkey[], mycornerfatherkey[] and mycornersonkey[] */
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
  {
    refinement->mycornerkey[n] = KeyForObject((KEY_OBJECT *)CORNER(theElement,i));
    refinement->mycornerfatherkey[n] = KeyForObject((KEY_OBJECT *)NFATHER(CORNER(theElement,i)));
    refinement->mycornersonkey[n] = KeyForObject((KEY_OBJECT *)SONNODE(CORNER(theElement,i)));
    n++;
  }

  if (MGIO_PARFILE)
  {
    for (i=0; i<nmax; i++)
    {
      if ((sonex>>i)&1)
      {                         /* son exists */
        refinement->sonskey[i] = KeyForObject((KEY_OBJECT *)SonList[i]);
        for (j=0; j<SIDES_OF_ELEM(SonList[i]); j++)
          refinement->sonsnbkey[i][j] = KeyForObject((KEY_OBJECT *)NBELEM(SonList[i],j));
      }
      else
        refinement->sonskey[i] = -1;
    }
  }
#endif

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
      if (GetNodeContext(theElement,NodeContext)) REP_ERR_RETURN(1);
  }
  for (; i<CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement)+SIDES_OF_ELEM(theElement); i++)
    if (NodeContext[i]!=NULL)
      if (GetNodeContext(theElement,NodeContext)) REP_ERR_RETURN(1);

  return (0);
}

static INT RemoveOrphanSons (ELEMENT **SonList, INT *nmax)
{
  INT i,max;

  if (nmax != NULL)
  {
    for (max=0,i=0; i<*nmax; i++)
    {
      if (SonList[i]!=NULL)
        if (THEFLAG(SonList[i]))
          SonList[i] = NULL;
        else
          max = i+1;
    }
    *nmax = max;
  }
  else
  {
    for (i=0; SonList[i]!=NULL; i++)
      if (SonList[i]!=NULL && THEFLAG(SonList[i])) SonList[i] = NULL;
    for (max=0,i=0; i<MAX_SONS; i++)
      if (SonList[i]!=NULL)
      {
        if (max<i) SonList[max] = SonList[i];
        max++;
      }
  }

  return(0);
}

static INT SetHierRefinement (GRID *theGrid, ELEMENT *theElement, MGIO_REFINEMENT *refinement, INT *RefRuleOffset)
{
  INT i,nmax;
  ELEMENT *theSon;
  ELEMENT *SonList[MAX_SONS];
  NODE *NodeContext[MAX_NEW_CORNERS_DIM+MAX_CORNERS_OF_ELEM];
  MGIO_RR_RULE *theRule;

  /*PRINTDEBUG(gm,0,(PFMT "SetHierRefinement(): level=%d elem=" EID_FMTX "\n",me,LEVEL(theGrid),EID_PRTX(theElement)));*/

  /* sequential part */
  if (REFINE(theElement)==NO_REFINEMENT) return (0);
  if (GetNodeContext(theElement,NodeContext)) REP_ERR_RETURN(1);
  theRule = rr_rules + rr_rule_offsets[TAG(theElement)] + REFINE(theElement);
  if (GetOrderedSons(theElement,theRule,NodeContext,SonList,&nmax)) REP_ERR_RETURN(1);
  if (RemoveOrphanSons(SonList,&nmax)) REP_ERR_RETURN(1);
  if (SetRefinement (theGrid,theElement,NodeContext,SonList,nmax,refinement,RefRuleOffset)) REP_ERR_RETURN(1);
  if (Write_Refinement (refinement,rr_rules)) REP_ERR_RETURN(1);

  /* recursive call */
  for (i=0; i<nmax; i++)
  {
    theSon = SonList[i];
    if (theSon==NULL) continue;
    if (REFINE(theSon)!=NO_REFINEMENT)
      if (SetHierRefinement(theGrid,theSon,refinement,RefRuleOffset))
        REP_ERR_RETURN(1);
  }

  return (0);
}

static INT nRefinements (ELEMENT *theElement, INT *n)
{
  INT i;
  ELEMENT *SonList[MAX_SONS];

  if (REFINE(theElement)==NO_REFINEMENT) return (0);
  if (GetAllSons(theElement,SonList)) REP_ERR_RETURN(1);
  if (RemoveOrphanSons(SonList,NULL)) REP_ERR_RETURN(1);
  (*n)++;
  for (i=0; SonList[i]!=NULL; i++)
    if (nRefinements(SonList[i],n)) REP_ERR_RETURN(1);

  return (0);
}

static INT WriteCG_Vertices (MULTIGRID *theMG, INT MarkKey)
{
  INT i,j,n;
  MGIO_CG_POINT *cg_point,*cgp;
  HEAP *theHeap;
  VERTEX *theVertex;

  theHeap = MGHEAP(theMG);
  n = nov*(MGIO_CG_POINT_SIZE);
  cg_point = (MGIO_CG_POINT *)GetTmpMem(theHeap,n,MarkKey);
  if (cg_point==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_points\n",(int)n); REP_ERR_RETURN(1);}
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

  if (Write_CG_Points((int)nov,cg_point)) REP_ERR_RETURN(1);

  return (0);
}

#ifdef ModelP
static INT WriteElementParInfo (GRID *theGrid,
                                ELEMENT *theElement, MGIO_PARINFO *pinfo)
{
  INT i,j,k,s,n_max;
  int *pl;
  NODE *theNode;
  VERTEX *theVertex;
  EDGE *theEdge;

  memset(pinfo,0,sizeof(MGIO_PARINFO));

  n_max = PROCLISTSIZE-(ActProcListPos-ProcList);

  s=0;
  pinfo->prio_elem = EPRIO(theElement);
  pinfo->ncopies_elem = ENCOPIES(theElement);
  if (n_max<pinfo->ncopies_elem)
  {
    PrintErrorMessage('E',"WriteElementParInfo","increase PROCLISTSIZE in gm/ugio.c\n");
    REP_ERR_RETURN(1);
  }
  if (pinfo->ncopies_elem>0)
  {
    pl = EPROCLIST(theElement);
    for (i=0,j=2; i<pinfo->ncopies_elem; i++,j+=2)
      ActProcListPos[s++] = pl[j];
  }
  pinfo->e_ident = EGID(theElement);
  for (k=0; k<CORNERS_OF_ELEM(theElement); k++)
  {
    theNode = CORNER(theElement,k);
    pinfo->prio_node[k] = PRIO(theNode);
    pinfo->ncopies_node[k] = NCOPIES(theNode);
    if (n_max<pinfo->ncopies_node[k]+s)
    {
      PrintErrorMessage('E',"WriteElementParInfo","increase PROCLISTSIZE in gm/ugio.c\n");
      REP_ERR_RETURN(1);
    }
    if (pinfo->ncopies_node[k]>0)
    {
      pl = PROCLIST(theNode);
      for (i=0,j=2; i<pinfo->ncopies_node[k]; i++,j+=2)
        ActProcListPos[s++] = pl[j];
    }
    pinfo->n_ident[k] = GID(theNode);
  }
  for (k=0; k<CORNERS_OF_ELEM(theElement); k++)
  {
    theVertex = MYVERTEX(CORNER(theElement,k));
    pinfo->prio_vertex[k] = VXPRIO(theVertex);
    pinfo->ncopies_vertex[k] = VXNCOPIES(theVertex);
    if (n_max<pinfo->ncopies_vertex[k]+s)
    {
      PrintErrorMessage('E',"WriteElementParInfo","increase PROCLISTSIZE in gm/ugio.c\n");
      REP_ERR_RETURN(1);
    }
    if (pinfo->ncopies_vertex[k]>0)
    {
      pl = VXPROCLIST(theVertex);
      for (i=0,j=2; i<pinfo->ncopies_vertex[k]; i++,j+=2)
        ActProcListPos[s++] = pl[j];
    }
    pinfo->v_ident[k] = VXGID(theVertex);
  }

#ifdef __TWODIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
  {
    VECTOR *v;
    for (k=0; k<EDGES_OF_ELEM(theElement); k++) {
      theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,k,0)),
                      CORNER(theElement,CORNER_OF_EDGE(theElement,k,1)));
      v = EDVECTOR(theEdge);
      pinfo->prio_edge[k] = PRIO(v);
      pinfo->ncopies_edge[k] = NCOPIES(v);
      if (n_max<pinfo->ncopies_edge[k]+s)
      {
        PrintErrorMessage('E',"WriteElementParInfo","increase PROCLISTSIZE in gm/ugio.c\n");
        REP_ERR_RETURN(1);
      }
      pinfo->ed_ident[k] = GID(v);
      if (pinfo->ncopies_edge[k]>0) {
        pl = PROCLIST(v);
        for (i=0,j=2; i<pinfo->ncopies_edge[k]; i++,j+=2)
          ActProcListPos[s++] = pl[j];
      }
    }
  }
#endif

#ifdef __THREEDIM__
  for (k=0; k<EDGES_OF_ELEM(theElement); k++)
  {
    theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,k,0)),
                      CORNER(theElement,CORNER_OF_EDGE(theElement,k,1)));
    pinfo->prio_edge[k] = PRIO(theEdge);
    pinfo->ncopies_edge[k] = NCOPIES(theEdge);
    if (n_max<pinfo->ncopies_edge[k]+s)
    {
      PrintErrorMessage('E',"WriteElementParInfo","increase PROCLISTSIZE in gm/ugio.c\n");
      REP_ERR_RETURN(1);
    }
    if (pinfo->ncopies_edge[k]>0) {
      pl = PROCLIST(theEdge);
      for (i=0,j=2; i<pinfo->ncopies_edge[k]; i++,j+=2)
        ActProcListPos[s++] = pl[j];
    }
    pinfo->ed_ident[k] = GID(theEdge);
  }
#endif

  if (s>0) pinfo->proclist = ActProcListPos;
  else pinfo->proclist = NULL;
  ActProcListPos += s;

  return (0);
}
#else
static INT WriteElementParInfo (GRID *theGrid, ELEMENT *theElement, MGIO_PARINFO *pinfo)
{
  REP_ERR_RETURN(1);
}
#endif

static INT SaveMultiGrid_SPF (MULTIGRID *theMG, char *name, char *type, char *comment, INT autosave, INT rename)
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
  MGIO_CG_ELEMENT *cg_element,*cge;
  MGIO_REFINEMENT *refinement;
  MGIO_BD_GENERAL bd_general;
  MGIO_PARINFO cg_pinfo;
  INT i,j,k,niv,nbv,nie,nbe,n,nref,hr_max,mode,level,id,foid,non,tl,saved;
  int *vidlist;
  char *p,*f,*s,*l;
  BNDP **BndPList;
  char filename[NAMESIZE];
  char buf[64],itype[10];
  int lastnumber;
  INT MarkKey;
#ifdef ModelP
  int ftype;
  int error;
#endif

  /* check */
  if (theMG==NULL) REP_ERR_RETURN(1);
  if (!MG_COARSE_FIXED(theMG))
  {
    PrintErrorMessage('E',"SaveMultiGrid_SPF","coarse grid not fixed\n");
    REP_ERR_RETURN(1);
  }

  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap,&MarkKey);

  /* something to do ? */
  saved = MG_SAVED(theMG);
#ifdef ModelP
  saved = UG_GlobalMinINT(saved);
  proc_list_size = PROCLISTSIZE_VALUE;
#endif
  if (saved)
  {
    UserWriteF("WARNING: multigrid already saved as %s\n",MG_FILENAME(theMG));
  }
  /* open file */
  nparfiles = procs;
  if (autosave)
  {
    if (name==NULL)
    {
      if (type!=NULL) REP_ERR_RETURN(1);
      strcpy(filename,MG_FILENAME(theMG));
      f = strtok(filename,".");       if (f==NULL) REP_ERR_RETURN(1);
      s = strtok(NULL,".");           if (s==NULL) REP_ERR_RETURN(1);
      l = strtok(NULL,".");           if (l==NULL) REP_ERR_RETURN(1);
      l = strtok(NULL,".");           if (l==NULL) REP_ERR_RETURN(1);
      l = strtok(NULL,".");           if (l==NULL) REP_ERR_RETURN(1);
      strtok(l,"/");
      if (sscanf(s,"%d",&lastnumber)!=1) REP_ERR_RETURN(1);
      if (lastnumber<0) REP_ERR_RETURN(1);
      lastnumber++;
      strcpy(itype,l);
      sprintf(buf,".%06d",lastnumber);
      strcat(filename,buf);
    }
    else
    {
      if (type==NULL) REP_ERR_RETURN(1);
      strcpy(itype,type);
      if (name==NULL) REP_ERR_RETURN(1);
      strcpy(filename,name);
      strcat(filename,".000000");
    }
  }
  else
  {
    if (type==NULL) REP_ERR_RETURN(1);
    strcpy(itype,type);
    if (name==NULL) REP_ERR_RETURN(1);
    strcpy(filename,name);
  }
  if (strcmp(itype,"dbg")==0) mode = BIO_DEBUG;
  else if (strcmp(itype,"asc")==0) mode = BIO_ASCII;
  else if (strcmp(itype,"bin")==0) mode = BIO_BIN;
  else REP_ERR_RETURN(1);
  sprintf(buf,".ug.mg.");
  strcat(filename,buf);
  strcat(filename,itype);
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
    REP_ERR_RETURN(1);
  }
  if (MGIO_PARFILE)
  {
    sprintf(buf,"/mg.%04d",(int)me);
    strcat(filename,buf);
  }
#endif
  if (Write_OpenMGFile (filename,(int)rename)) REP_ERR_RETURN(1);

  /* write general information */
  theBVP = MG_BVP(theMG);
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc)) REP_ERR_RETURN(1);
  mg_general.mode                 = mode;
  mg_general.dim                  = DIM;
  Broadcast(&MG_MAGIC_COOKIE(theMG),sizeof(INT));
  mg_general.magic_cookie = MG_MAGIC_COOKIE(theMG);
  mg_general.heapsize             = MGHEAP(theMG)->size/KBYTE;
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

  if (Write_MG_General(&mg_general)) REP_ERR_RETURN(1);

  /* write information about general elements */
  ge_general.nGenElement = TAGS;
  if (Write_GE_General(&ge_general)) REP_ERR_RETURN(1);
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
  if (Write_GE_Elements(TAGS,ge_element)) REP_ERR_RETURN(1);

  /* write information about refrules used */
  /* TODO (HRR 971209): drop old call of Write_RefRules */
    #if EXTRACT_RULES
  if (NEW_Write_RefRules(theMG,rr_rule_offsets,MarkKey,&rr_rules)) REP_ERR_RETURN(1);
        #else
  if (Write_RefRules(theMG,rr_rule_offsets,MarkKey)) REP_ERR_RETURN(1);
        #endif

  i = OrphanCons(theMG);
  if (i)
  {
    PRINTDEBUG(gm,1,(PFMT "OrphanCons() returned %d errors\n",me,i));
    fflush(stdout);
    assert(0);
  }

  /* renumber objects */
  if (MGIO_PARFILE)
  {
    if (RenumberMultiGrid (theMG,&nbe,&nie,&nbv,&niv,&vid_n,&foid,&non,MarkKey)) {UserWriteF("ERROR: cannot renumber multigrid\n"); REP_ERR_RETURN(1);}
  }
  else
  {
    if (RenumberMultiGrid (theMG,&nbe,&nie,&nbv,&niv,NULL,&foid,&non,MarkKey)) {UserWriteF("ERROR: cannot renumber multigrid\n"); REP_ERR_RETURN(1);}
  }
  nov = nbv+niv;

  /* write general information about coarse grid */
  theGrid = GRID_ON_LEVEL(theMG,0);
  cg_general.nPoint = nbv+niv;
  cg_general.nBndPoint = nbv;
  cg_general.nInnerPoint = niv;
  cg_general.nElement = nbe+nie;
  cg_general.nBndElement = nbe;
  cg_general.nInnerElement = nie;
  if (Write_CG_General(&cg_general)) REP_ERR_RETURN(1);

#ifdef ModelP
  /* this proc has no elements */
  if (cg_general.nElement == 0)
  {
    /* release tmp mem */
    ReleaseTmpMem(theHeap,MarkKey);

    /* close file */
    if (CloseMGFile ()) REP_ERR_RETURN(1);

    /* saved */
    MG_SAVED(theMG) = 1;
    strcpy(MG_FILENAME(theMG),filename);

    return(0);
  }
#endif

  /* write coarse grid points */
  if (WriteCG_Vertices(theMG,MarkKey)) REP_ERR_RETURN(1);

  /* write mapping: node-id --> vertex-id for orphan-nodes, if(MGIO_PARFILE) */
  if (MGIO_PARFILE)
  {
    vidlist = (int*)GetTmpMem(theHeap,(2+non)*sizeof(int),MarkKey);
    vidlist[0] = non;
    vidlist[1] = foid;
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
        if (USED(theNode))
          vidlist[ID(theNode)+2-foid] = ID(MYVERTEX(theNode));
    if (Bio_Write_mint((int)(2+non),vidlist)) REP_ERR_RETURN(1);
  }

  /* write orphan elements */
  n = cg_general.nElement; hr_max=0;
  cg_element = (MGIO_CG_ELEMENT*)GetTmpMem(theHeap,n*MGIO_CG_ELEMENT_SIZE,MarkKey);
  if (cg_element==NULL)   {UserWriteF("ERROR: cannot allocate %d bytes for cg_elements\n",(int)n*MGIO_CG_ELEMENT_SIZE); REP_ERR_RETURN(1);}
  if (MGIO_PARFILE)
  {
    ActProcListPos = ProcList = (unsigned short*)GetTmpMem(theHeap,PROCLISTSIZE*sizeof(unsigned short),MarkKey);
    if (ProcList==NULL)     {UserWriteF("ERROR: cannot allocate %d bytes for ProcList\n",(int)PROCLISTSIZE*sizeof(unsigned short)); REP_ERR_RETURN(1);}
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
      if (!EORPHAN(theElement)) continue;
      assert(id==ID(theElement));
      assert(id<n);

      if (MGIO_PARFILE)
      {
        /* parallel part */
        cge->level = LEVEL(theElement);
      }

      /* --- */
      cge->ge = TAG(theElement);
      nref=0;
      if (nRefinements (theElement,&nref)) REP_ERR_RETURN(1);
      hr_max = MAX(hr_max,nref);
      cge->nref = nref;
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
        cge->cornerid[j] = ID(CORNER(theElement,j));
      for (j=0; j<SIDES_OF_ELEM(theElement); j++)
        if (NBELEM(theElement,j)!=NULL)
          cge->nbid[j] = ID(NBELEM(theElement,j));
        else
          cge->nbid[j] = -1;
      cge->side_on_bnd = 0;
      if (OBJT(theElement)==BEOBJ)
        for (j=0; j<SIDES_OF_ELEM(theElement); j++)
          if (SIDE_ON_BND(theElement,j))
            cge->side_on_bnd |= (1<<j);
      cge->subdomain = SUBDOMAIN(theElement);

#if (MGIO_DEBUG>0)
      /* write debug extension */
      cge->mykey = KeyForObject((KEY_OBJECT *)theElement);
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
        cge->nodekey[j] = KeyForObject((KEY_OBJECT *)CORNER(theElement,j));
      for (j=0; j<SIDES_OF_ELEM(theElement); j++)
        cge->neighborkey[j] = KeyForObject((KEY_OBJECT *)NBELEM(theElement,j));
#endif

      /* increment counters */
      i++; id++;                                                                                                        /* id this is the id of theElement, according to RenumberElements (ugm.c) */
    }
  if (Write_CG_Elements((int)i,cg_element)) REP_ERR_RETURN(1);

  /* write general bnd information */
  if (Bio_Jump_From ()) REP_ERR_RETURN(1);
  bd_general.nBndP = nbv;
  if (Write_BD_General (&bd_general)) REP_ERR_RETURN(1);

  /* write bnd information */
  if (nbv > 0) {
    BndPList = (BNDP**)GetTmpMem(theHeap,nbv*sizeof(BNDP*),MarkKey);
    if (BndPList==NULL) {
      UserWriteF("ERROR: cannot allocate %d bytes for BndPList\n",
                 (int)nbv*sizeof(BNDP*));
      REP_ERR_RETURN(1);
    }
    if (procs>1) i=TOPLEVEL(theMG);
    else i=0;
    for (level=0; level<=i; level++)
      for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,level)); theNode!=NULL; theNode=SUCCN(theNode))
        if (USED(MYVERTEX(theNode)) && OBJT(MYVERTEX(theNode))==BVOBJ) {
          /* see ugm.c: RenumberVertices */
          id = ID(MYVERTEX(theNode));
          if (id<0 || id>=nbv) REP_ERR_RETURN(1);
          BndPList[id] = V_BNDP(MYVERTEX(theNode));
          SETUSED(MYVERTEX(theNode),0);
        }
    if (Write_PBndDesc (nbv,BndPList)) REP_ERR_RETURN(1);
  }
  if (Bio_Jump_To ()) REP_ERR_RETURN(1);

  /* write parinfo of coarse-grid */
  if (MGIO_PARFILE)
  {
    cg_pinfo.proclist = ProcList;
    for (level=0; level<=TOPLEVEL(theMG); level++)
      for (theElement = PFIRSTELEMENT(GRID_ON_LEVEL(theMG,level)); theElement!=NULL; theElement=SUCCE(theElement))
        if (EORPHAN(theElement))
        {
          if (WriteElementParInfo(theGrid, theElement,&cg_pinfo)) REP_ERR_RETURN(1);
          if (Write_pinfo (TAG(theElement),&cg_pinfo)) REP_ERR_RETURN(1);
        }
  }

  /* save refinement */
  refinement = (MGIO_REFINEMENT *)GetTmpMem(theHeap,(hr_max+1)*MGIO_REFINEMENT_SIZE,MarkKey);                   /* size according to procs==1 or procs>1 (see mgio.h) */
  if (refinement==NULL) {UserWriteF("ERROR: cannot allocate %ld bytes for refinement\n",(long)hr_max*MGIO_REFINEMENT_SIZE); REP_ERR_RETURN(1);}
  if (procs>1) tl=TOPLEVEL(theMG);
  else tl=0;

  id=0;
  for (level=0; level<=tl; level++) {
    theGrid = GRID_ON_LEVEL(theMG,level);
    for (theElement = PFIRSTELEMENT(theGrid);
         theElement!=NULL; theElement=SUCCE(theElement)) {
      if (!EORPHAN(theElement)) continue;
      assert(id==ID(theElement));
      if (REFINE(theElement)!=NO_REFINEMENT)
        /* write refinement only if it is neccessary */
        if(SetHierRefinement(theGrid,theElement,refinement,rr_rule_offsets))
          REP_ERR_RETURN(1);
      id++;
    }
  }

  /* release tmp mem */
  ReleaseTmpMem(theHeap,MarkKey);

  /* close file */
  if (CloseMGFile ()) REP_ERR_RETURN(1);

  /* saved */
  MG_SAVED(theMG) = 1;
  strcpy(MG_FILENAME(theMG),filename);

        #if EXTRACT_RULES
  ResetRefineTagsBeyondRuleManager(theMG);
        #endif

  return (0);
}

INT SaveMultiGrid (MULTIGRID *theMG, char *name, char *type, char *comment, INT autosave, INT rename)
{
  /* check name convention */
  if (name==NULL || strcmp(name+strlen(name)-4,".scr")!=0)
  {
    if (SaveMultiGrid_SPF (theMG,name,type,comment,autosave,rename)) REP_ERR_RETURN(1);
  }
  else
  {
    if (SaveMultiGrid_SCR (theMG,name,comment)) REP_ERR_RETURN(1);
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
  INT i,j,s,prio,where,oldwhere,old;
  INT evec,nvec,edvec,svec;
  GRID            *vgrid;
  ELEMENT         *theFather,*After,*Next,*Succe;
  NODE            *theNode;
  VERTEX          *theVertex;
  VECTOR          *theVector;
  EDGE            *theEdge;

  evec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),ELEMVEC);
  nvec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),NODEVEC);
  edvec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),EDGEVEC);
#ifdef __THREEDIM__
  svec = VEC_DEF_IN_OBJ_OF_MG(MYMG(theGrid),SIDEVEC);
#else
  svec = 0;
#endif
  /* this function does not support side vectors                        */
  /* proclist and identificator need to be stored for each  side vector */
  if (svec) assert(0);

  theFather = EFATHER(theElement);
  s = 0;
  if ((prio = pinfo->prio_elem) != PrioMaster)
  {
    old = EPRIO(theElement);
    oldwhere = PRIO2INDEX(old);
    Succe = SUCCE(theElement);
    GRID_UNLINK_ELEMENT(theGrid,theElement);
    SETEPRIO(theElement,prio);
    if (theFather != NULL)
    {
      if (theElement == SON(theFather,oldwhere))
      {
        Next = NULL;
        if (Succe != NULL)
          if (EFATHER(Succe)==theFather && PRIO2INDEX(EPRIO(Succe))==oldwhere)
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
      PRINTDEBUG(gm,1,("Evaluate-pinfo():nid=%d prio=%d\n",ID(theNode),prio);fflush(stdout));
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
  for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
  {
    theVertex = MYVERTEX(CORNER(theElement,j));
    if (USED(theVertex) == 0)
    {
      vgrid = GRID_ON_LEVEL(MYMG(theGrid),LEVEL(theVertex));
      if ((prio = pinfo->prio_vertex[j]) != PrioMaster)
      {
        GRID_UNLINK_VERTEX(vgrid,theVertex);
        SETVXPRIO(theVertex,prio);
        GRID_LINK_VERTEX(vgrid,theVertex,prio);
      }
      PRINTDEBUG(gm,1,("Evaluate-pinfo():vid=%d prio=%d\n",ID(theVertex),prio);fflush(stdout));
      for (i=0; i<pinfo->ncopies_vertex[j]; i++)
      {
        DDD_IdentifyNumber(PARHDRV(theVertex),pinfo->proclist[s],pinfo->v_ident[j]);
        s++;
      }
      SETUSED(theVertex,1);
    }
    else
      s += pinfo->ncopies_vertex[j];
  }

#ifdef __TWODIM__
  if (edvec) {
    for (j=0; j<EDGES_OF_ELEM(theElement); j++) {
      theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),
                      CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
      if (USED(theEdge) == 0) {
        theVector = EDVECTOR(theEdge);
        if ((prio = pinfo->prio_edge[j]) != PrioMaster) {
          GRID_UNLINK_VECTOR(theGrid,theVector);
          SETPRIO(theVector,prio);
          GRID_LINK_VECTOR(theGrid,theVector,prio);
        }
        for (i=0; i<pinfo->ncopies_edge[j]; i++) {
          DDD_IdentifyNumber(PARHDR(theVector),
                             pinfo->proclist[s],pinfo->ed_ident[j]);
          s++;
        }
        SETUSED(theEdge,1);
      }
      else
        s += pinfo->ncopies_edge[j];
    }
  }
#endif

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
          theVector = EVECTOR(theElement);
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

#ifdef ModelP
static int Gather_RefineInfo (DDD_OBJ obj, void *data)
{
  ELEMENT *theElement = (ELEMENT *)obj;

  ((int *)data)[0] = REFINECLASS(theElement);
  ((int *)data)[1] = REFINE(theElement);
  ((int *)data)[2] = MARKCLASS(theElement);
  ((int *)data)[3] = MARK(theElement);

  return(GM_OK);
}

static int Scatter_RefineInfo (DDD_OBJ obj, void *data)
{
  ELEMENT *theElement = (ELEMENT *)obj;

  SETREFINECLASS(theElement,((int *)data)[0]);
  SETREFINE(theElement,((int *)data)[1]);
  SETMARKCLASS(theElement,((int *)data)[2]);
  SETMARK(theElement,((int *)data)[3]);

  return(GM_OK);
}

static INT SpreadRefineInfo(GRID *theGrid)
{
  DDD_IFAOneway(ElementIF,GRID_ATTR(theGrid),IF_FORWARD,4*sizeof(INT),
                Gather_RefineInfo,Scatter_RefineInfo);
  return(GM_OK);
}
#endif

static int Gather_NodeType (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  ((int *)data)[0] = NTYPE(theNode);

  return(GM_OK);
}

static int Scatter_NodeType (DDD_OBJ obj, void *data)
{
  NODE *theNode = (NODE *)obj;

  SETNTYPE(theNode,((int *)data)[0]);

  return(GM_OK);
}

static INT SpreadGridNodeTypes(GRID *theGrid)
{
  DDD_IFAOneway(NodeIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(INT),
                Gather_NodeType,Scatter_NodeType);
  return(GM_OK);
}

static INT IO_GridCons(MULTIGRID *theMG)
{
  INT i,*proclist;
  GRID    *theGrid;
  ELEMENT *theElement;
  VECTOR  *theVector;

  for (i=TOPLEVEL(theMG); i>=0; i--)         /* propagate information top-down */
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
      ASSERT((PARTITION(theElement)==me && EMASTER(theElement))
             || (PARTITION(theElement)!=me && !EMASTER(theElement)));
    }
    for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
      if (!MASTER(theVector))
        DisposeConnectionFromVector(theGrid,theVector);

#ifdef ModelP
    /* spread element refine info */
    if (SpreadRefineInfo(theGrid) != GM_OK) RETURN(GM_FATAL);

    /* spread nodetypes from master to its copies */
    if (SpreadGridNodeTypes(theGrid) != GM_OK) RETURN(GM_FATAL);

    /* repair parallel information */
    ConstructConsistentGrid(theGrid);
#endif
  }

  return(GM_OK);
}


#if (MGIO_DEBUG>0)
static void CheckNeighborElement( ELEMENT *theElement, INT nb_nr, INT expected_nb_key, const char *messagePrefix )
/* Checks wheather the relation between an element and a given neighbor is ok.
   If an error is detected, the programm stopps with 'assert'.
 */
{
  INT k, key;
  ELEMENT *nb_back, *nb_element;

  nb_element = NBELEM(theElement,nb_nr);

  if (nb_element!=NULL)
  {
    key = KeyForObject((KEY_OBJECT *)nb_element);
    if (key!=expected_nb_key)
    {
      if ( EGHOST(theElement) && EGHOST(nb_element) )
      {                         /* if both elements are ghosts, then the neighbor-pointers are
                                   optionally; but if there is one pointer then also the
                                   other must exist.
                                   if at save time none of this pointers have existed,
                                   at load time they may exist and it is ok.
                                 */

        if ( expected_nb_key == -1 )
        {                               /* at save time there were no neighbor pointers
                                           ( i.e. expected_nb_key == -1 ), but now.
                                           ensure that the pointers are symmetric.
                                         */

          /* search for the neighbor pointer from nb_element
             back to theElement */
          nb_back = NULL;
          for (k=0; k<SIDES_OF_ELEM(nb_element); k++)
            if (NBELEM(nb_element,k)==theElement)
            {
              nb_back = theElement;
              break;
            }

          if ( nb_back != theElement )
          {
            printf(PFMT " %s neighbor relation unsymmetric, element: " EID_FMTX "has neighbor[%d]: " EID_FMTX " but not vice versa (this pointer has not been in the coarse grid at save time, but it's ok)\n",
                   me,messagePrefix,EID_PRTX(theElement),nb_nr,EID_PRTX(nb_element));
            assert(0);
          }
          else
          {                                     /* else the neighbor pointer relation is symmetric and hence ok */
            return;
          }
        }
        else                         /* it is an error because there are 2 valid (!=-1) but different keys */
        {
          printf(PFMT " %s neighbor relation with different keys for ghosts, element: " EID_FMTX "has neighbor[%d]: " EID_FMTX "\n",
                 me,messagePrefix,EID_PRTX(theElement),nb_nr,EID_PRTX(nb_element));
          assert(0);
        }
      }
      else
      {
        /* else it is an error because at least one of the 2 considered elements
           is a master and masters must always have pointers to its
           existing neighbors. Thus especially it can't be, that at
           save time expected_nb_key==-1 has held and a key difference
           is in this case always an error.
         */
        printf(PFMT " %s neighbor relation with different keys for master/ghost, element: " EID_FMTX "has neighbor[%d]: " EID_FMTX ":expected key %d does not match\n",
               me,messagePrefix,EID_PRTX(theElement),nb_nr,EID_PRTX(nb_element),expected_nb_key);
        assert(0);
      }
    }
  }
  /* you cannot expect that the neighbors are already existing
      else
      {
              if ( expected_nb_key!=-1 && must_exist )
              {
                      printf(PFMT " %s neighbor relation, element: " EID_FMTX ", neighbor %d doesn't exist but key %d expected\n",me,messagePrefix,EID_PRTX(theElement),nb_nr,ref->mykey);
                      assert(0);
              }
      }
   */
}

static INT CheckLocalElementKeys (ELEMENT *theElement, MGIO_REFINEMENT *ref, INT must_exist)
/* check as many consistencies as possible. Various keys (for neibhbours,
   sons, ...) are stored in the debug mode and are now compared to the
   existing objects in the memory.
   if must_exist == TRUE all sons must exist. */
{
  INT i, j, k, key, nmax;
  NODE *corner_node;
  ELEMENT *SonList[MAX_SONS], *nb_element, *nb_back;

  /* check consistency of element key */
  key = KeyForObject((KEY_OBJECT *)theElement);
  if (key!=ref->mykey)
  {
    printf(PFMT " IO_Loc: element key, element: " EID_FMTX ": exp.key %d does not match\n",me,EID_PRTX(theElement),ref->mykey);
    assert(0);
  }

  /* check consistency of father element key;
     the check is ok also for orphan elements (they have fatherkey -1) */
  key = KeyForObject((KEY_OBJECT *)EFATHER(theElement));
  if (key!=ref->myfatherkey)
  {
    printf(PFMT " IO_Loc: father element key, element: " EID_FMTX ", father element: " EID_FMTX ": exp.key %d does not match\n",me,EID_PRTX(theElement),EID_PRTX(EFATHER(theElement)),ref->mykey);
    assert(0);
  }

  /* check corners */
  if (CORNERS_OF_ELEM(theElement)!=ref->mycorners)
  {
    printf(PFMT " IO_Loc: number of corners, element: " EID_FMTX ": exp.key %d does not matchas %d corners but expected %dh\n",me,EID_PRTX(theElement),CORNERS_OF_ELEM(theElement),ref->mycorners);
    assert(0);
  }
  for (j=0; j<ref->mycorners; j++)
  {
    corner_node = CORNER(theElement,j);
    if (corner_node!=NULL)
    {
      if (ref->mycornerkey[j] != KeyForObject((KEY_OBJECT *)corner_node))
      {
        printf(PFMT " IO_Loc: corner relation, element: " EID_FMTX ", corner %d: " ID_FMTX " :exp.key %d does not match\n",
               me,EID_PRTX(theElement),j,ID_PRTX(corner_node),ref->mycornerkey[j]);
        assert(0);
      }

      /* check father corner */
      if (NFATHER(corner_node)!=NULL)
      {
        if ( ref->mycornerfatherkey[j] != KeyForObject((KEY_OBJECT *)NFATHER(corner_node)) )
        {
          printf(PFMT " IO_Loc: father corner relation, element: " EID_FMTX ", corner[%d]: " ID_FMTX " fathercorner: " ID_FMTX ":exp.key %d does not match\n",
                 me,EID_PRTX(theElement),j,ID_PRTX(corner_node),ID_PRTX(NFATHER(corner_node)),ref->mycornerfatherkey[j]);
          assert(0);
        }
      }
      else
      {
        if ( ref->mycornerfatherkey[j]!=-1 )
        {
          if ( !EORPHAN(theElement) )
          {
            printf(PFMT " IO_Loc: father corner relation, element: " EID_FMTX " corner %d doesn't exist but key %d expected\n",me,EID_PRTX(theElement),j,ref->mycornerfatherkey[j]);
            assert(0);
          }
          /* else for orphan element: at save time there was a father-corner. Such father-corner pointers are
             reconstructed not from the orphan element itself, but from the father-corner which is
             located in a neighbor element-tree; and since it is possible, that this neighbor element-tree
             is not already constructed and the father-corner pointer is still missing. It's OK.
           */
        }
      }

      /* check son corner */
      if (SONNODE(corner_node)!=NULL)
      {
        if ( ref->mycornersonkey[j] != KeyForObject((KEY_OBJECT *)SONNODE(corner_node)) )
        {
          printf(PFMT " IO_Loc: son corner relation, element: " EID_FMTX ", corner[%d]: " ID_FMTX " soncorner: " ID_FMTX ":exp.key %d does not match\n",
                 me,EID_PRTX(theElement),j,ID_PRTX(corner_node),ID_PRTX(SONNODE(corner_node)),ref->mycornersonkey[j]);
          assert(0);
        }
      }
      /* you cannot expect that the sonnodes are already existing
         else
         {
              if ( ref->mycornersonkey[j]!=-1 && must_exist)
              {
                      printf(PFMT " IO_Loc: son corner relation, element: " EID_FMTX ", corner[%d] " ID_FMTX ": son corner doesn't exist but key %d expected\n",me,EID_PRTX(theElement),j,ID_PRTX(corner_node),ref->mycornersonkey[j]);
                      assert(0);
              }
         }
       */
    }
    else
    {
      if ( ref->mycornerkey[j]!=-1 )
      {
        printf(PFMT " IO_Loc: corner relation, element: " EID_FMTX ", corner %d doesn't exist but key %d expected\n",me,EID_PRTX(theElement),j,ref->mycornerkey[j]);
        assert(0);
      }
    }

  }

  /* check nbkey[] */
  for (j=0; j<SIDES_OF_ELEM(theElement); j++)
    CheckNeighborElement( theElement, j, ref->nbkey[j], "IO_Loc:" );

  if (MGIO_PARFILE)
  {
    NODE *NodeContext[MAX_NEW_CORNERS_DIM+MAX_CORNERS_OF_ELEM];
    MGIO_RR_RULE *theRule;

    if (GetNodeContext(theElement,NodeContext)) REP_ERR_RETURN(1);
    theRule = rr_rules + rr_rule_offsets[TAG(theElement)] + REFINE(theElement);
    if (GetOrderedSons(theElement,theRule,NodeContext,SonList,&nmax)) REP_ERR_RETURN(1);

    /* check the sons */
    for (i=0; i<nmax; i++)
    {
      if ( SonList[i]==NULL )
      {
        if (!must_exist)
          continue;

        if (((1<<i) & ref->sonex)!=0)
        {
          printf(PFMT " IO_Loc: sonlist relation, element: " EID_FMTX ", son %d doesn't exist but sonex 0x%x set\n",me,EID_PRTX(theElement),i,ref->sonex);
          assert(0);
        }
        if ( ref->sonskey[i] != -1 )
        {
          printf(PFMT " IO_Loc: sonlist relation, element: " EID_FMTX ", son %d doesn't exist but key %d expected\n",me,EID_PRTX(theElement),i,ref->sonskey[i]);
          assert(0);
        }
        continue;
      }

      if (((1<<i) & ref->sonex)==0)
      {
        /* son exists, but is not set in sonex */
        /* TODO: check whether SonList[i] is really a (orphan) coarse grid element */
        continue;
      }


      /* check consistency of son element key */
      key = KeyForObject((KEY_OBJECT *)SonList[i]);
      if (key!=ref->sonskey[i])
      {
        printf(PFMT " IO_Loc: son element key, element: " EID_FMTX ", son[%d] " EID_FMTX ": exp.key %d does not match\n",me,EID_PRTX(theElement),i,EID_PRTX(SonList[i]),ref->sonskey[i]);
        assert(0);
      }

      for (j=0; j<SIDES_OF_ELEM(SonList[i]); j++)
        CheckNeighborElement( SonList[i], j, ref->sonsnbkey[i][j], "IO_Loc: son" );
    }
  }

  return (0);
}
#else
static INT CheckLocalElementKeys (ELEMENT *el, MGIO_REFINEMENT *ref, INT must_exist)
{
  return (0);
}
#endif


static INT InsertLocalTree (GRID *theGrid, ELEMENT *theElement, MGIO_REFINEMENT *ref, int *RefRuleOffset)
{
  INT i,j,k,r_index,nedge,type,sonRefined,n0,n1,Sons_of_Side,SonSides[MAX_SONS],offset;
  ELEMENT *theSonElem[MAX_SONS];
  ELEMENT *theSonList[MAX_SONS];
  VERTEX *theVertex;
  NODE *NodeList[MAX_NEW_CORNERS_DIM+MAX_CORNERS_OF_ELEM];
  NODE *SonNodeList[MAX_CORNERS_OF_ELEM];
  GRID *upGrid;
  EDGE *theEdge;
  MGIO_RR_RULE *theRule;
  struct mgio_sondata *SonData;
  INT nbside,nex;
#if (MGIO_DEBUG>0)
  MGIO_REFINEMENT ref_copy;
#endif

  /* read refinement */
  if (Read_Refinement(ref,rr_rules)) REP_ERR_RETURN(1);
#if (MGIO_DEBUG>0)
  ref_copy = *ref;      /* copy the read ref. into local buffer for subsequent checks */
#endif

  /*PRINTDEBUG(gm,0,(PFMT "InsertLocalTree(): level=%d elem=" EID_FMTX " REFINECLASS %d\n",me,LEVEL(theGrid),EID_PRTX(theElement),ref->refclass));*/

  /* init */
  if (ref->refrule==-1) REP_ERR_RETURN(1);
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
        if (NodeList[i]==NULL) REP_ERR_RETURN(1);
        ID(NodeList[i]) = ref->newcornerid[r_index];
      }
    }
    else
    {
      ASSERT(ID(NodeList[i]) == ref->newcornerid[r_index]);
    }
    SETNTYPE(NodeList[i],CORNER_NODE);
    r_index++;
  }
  offset = i;

  nedge = EDGES_OF_ELEM(theElement);
  for (; i<nedge+offset; i++)
  {
    theVertex = NULL;
    if (MGIO_PARFILE)
    {
      if (!((nex>>i)&1))
      {
        NodeList[i] = NULL;
        continue;
      }
      else if (ref->orphanid_ex==1)
      {
        if (ref->orphanid[r_index]>=0)
        {
          assert(ref->orphanid[r_index]>=foid && ref->orphanid[r_index]<non+foid);
          theVertex = MYVERTEX(nid_n[ref->orphanid[r_index]-foid]);
        }
      }
    }
    if (theRule->pattern[i-offset]!=1)
    {
      NodeList[i] = NULL;
      continue;
    }
    n0 = CORNER_OF_EDGE(theElement,i-offset,0);
    n1 = CORNER_OF_EDGE(theElement,i-offset,1);
    theEdge = GetEdge(CORNER(theElement,n0),CORNER(theElement,n1));
    if (theEdge==NULL) REP_ERR_RETURN(1);
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
        NodeList[i] = CreateMidNode(upGrid,theElement,theVertex,i-offset);
        if (NodeList[i]==NULL) REP_ERR_RETURN(1);
        ID(NodeList[i]) = ref->newcornerid[r_index];
        HEAPFAULT(NodeList[i]);
      }
    }
    else
    {
      ASSERT(ID(NodeList[i]) == ref->newcornerid[r_index]);
    }
    SETNTYPE(NodeList[i],MID_NODE);
    HEAPFAULT(NodeList[i]);
    r_index++;
  }
  offset = i;

#ifdef __THREEDIM__
  for (; i<SIDES_OF_ELEM(theElement)+offset; i++)
  {
    if (theRule->pattern[i-CORNERS_OF_ELEM(theElement)]!=1)
    {
      NodeList[i] = NULL;
      continue;
    }
    theVertex = NULL;
    if (MGIO_PARFILE)
    {
      if (!((nex>>i)&1))
      {
        NodeList[i] = NULL;
        continue;
      }
      else if (ref->orphanid_ex==1)
      {
        if (ref->orphanid[r_index]>=0)
        {
          assert(ref->orphanid[r_index]>=foid && ref->orphanid[r_index]<non+foid);
          theVertex = MYVERTEX(nid_n[ref->orphanid[r_index]-foid]);
        }
      }
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
      NodeList[i] = CreateSideNode(upGrid,theElement,theVertex,i-offset);
      if (NodeList[i]==NULL) REP_ERR_RETURN(1);
      ID(NodeList[i]) = ref->newcornerid[r_index];
    }
    else
      SETONNBSIDE(MYVERTEX(NodeList[i]),i-offset);
    SETNTYPE(NodeList[i],SIDE_NODE);
    r_index++;
  }
  offset = i;
#endif

  i = CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement);       /* using consistency within rm.c (id of CenterNode in SonData corresponds to CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement)) */
  if (theRule->pattern[CENTER_NODE_INDEX(theElement)]==1 && (!MGIO_PARFILE || ((nex>>(i))&1)))
  {
    theVertex = NULL;
    if (MGIO_PARFILE && ref->orphanid_ex==1 && ref->orphanid[r_index]>=0)
    {
      assert(ref->orphanid[r_index]>=foid && ref->orphanid[r_index]<non+foid);
      theVertex = MYVERTEX(nid_n[ref->orphanid[r_index]-foid]);
    }
    if (ref->newcornerid[r_index]-foid>=0 && ref->newcornerid[r_index]-foid<non)
    {
      NodeList[i] = nid_n[ref->newcornerid[r_index]-foid];
      assert(NodeList[i]!=NULL);
      assert(ID(NodeList[i]) == ref->newcornerid[r_index]);
    }
    else
      NodeList[i] = CreateCenterNode(upGrid,theElement,theVertex);
    if (NodeList[i]==NULL) REP_ERR_RETURN(1);
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
        if (SonData->nb[j]>=MGIO_FATHER_SIDE_OFFSET && SIDE_ON_BND(theElement,SonData->nb[j]-MGIO_FATHER_SIDE_OFFSET))
        {
          type = BEOBJ;
          break;
        }
    for (j=0; j<CORNERS_OF_TAG(SonData->tag); j++)
      SonNodeList[j] = NodeList[SonData->corners[j]];
    theSonList[i] = CreateElement(upGrid,SonData->tag,type,SonNodeList,theElement,1);
    if (theSonList[i]==NULL) REP_ERR_RETURN(1);
    SETECLASS(theSonList[i],ref->refclass);
    SETSUBDOMAIN(theSonList[i],SUBDOMAIN(theElement));
  }

  /* identify objects */
  if (MGIO_PARFILE)
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
      if (SonData->nb[j]<MGIO_FATHER_SIDE_OFFSET && (!MGIO_PARFILE || ((ref->sonex>>(SonData->nb[j]))&1)))
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
    if (Get_Sons_of_ElementSide (theElement,i,&Sons_of_Side,theSonElem,SonSides,0,1)) REP_ERR_RETURN(1);

    if (Connect_Sons_of_ElementSide(UPGRID(theGrid),theElement,i,Sons_of_Side,theSonElem,SonSides,1)) REP_ERR_RETURN(1);
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

#if (MGIO_DEBUG>0)
  if (CheckLocalElementKeys(theElement, &ref_copy, FALSE)==0)
  {
    PRINTDEBUG(gm,4,("InsertLocalTree: preliminary CheckLocalElementKeys ok.\n"));
  }
  else
  {
    PrintErrorMessage('E',"InsertLocalTree","preliminary CheckLocalElementKeys failed\n");
  }
#endif

  /* jump to the sons ? */
  sonRefined = ref->sonref;

  /* call recursively */
  for (i=0; i<theRule->nsons; i++)
  {
    if (sonRefined & (1<<i))
    {
      if (InsertLocalTree (upGrid,theSonList[i],ref,RefRuleOffset)) REP_ERR_RETURN(1);
    }
    else if (theSonList[i] != NULL)
    {
      SETREFINE(theSonList[i],NO_REFINEMENT);
    }

#if (MGIO_DEBUG>0)
    if (CheckLocalElementKeys(theElement, &ref_copy, TRUE)==0)
    {
      PRINTDEBUG(gm,4,("InsertLocalTree: after son %d intermediate CheckLocalElementKeys ok.\n",i));
    }
    else
    {
      PrintErrorMessage('E',"InsertLocalTree","intermediate CheckLocalElementKeys failed\n");
    }
#endif
  }

#if (MGIO_DEBUG>0)
  if (CheckLocalElementKeys(theElement, &ref_copy, TRUE)==0)
  {
    PRINTDEBUG(gm,4,("InsertLocalTree: final CheckLocalElementKeys ok.\n"));
  }
  else
  {
    PrintErrorMessage('E',"InsertLocalTree","final CheckLocalElementKeys failed\n");
  }
#endif

  return (0);
}

#ifdef ModelP
static INT Gather_EClasses (DDD_OBJ obj, void *data)
{
  ELEMENT *p;
  int *d;

  p  = (ELEMENT *)obj;
  d  = (int *)data;
  *d = ECLASS(p);
}

static INT Scatter_EClasses(DDD_OBJ obj, void *data)
{
  ELEMENT *p;
  int *d;

  p = (ELEMENT *)obj;
  d = (int *)data;
  SETECLASS(p,*d);
}

void CommunicateEClasses (MULTIGRID *theMG)
{
  DDD_IFOneway(ElementVHIF,IF_FORWARD,sizeof(int),
               Gather_EClasses, Scatter_EClasses);
  return;
}
#endif

#if (MGIO_DEBUG>0)
static INT CheckCGKeys (INT ne, ELEMENT** eid_e, MGIO_CG_ELEMENT *cg_elem)
{
  INT i,j,k,key;
  MGIO_CG_ELEMENT *cge;
  ELEMENT *theElement, *nb_element, *nb_back;

  for (i=0; i<ne; i++)
  {
    cge = MGIO_CG_ELEMENT_PS(cg_elem,i);
    theElement = eid_e[i];

    /* check consistency of keys */
    key = KeyForObject((KEY_OBJECT *)theElement);
    if (key!=-1 && key!=cge->mykey)
    {
      printf(PFMT " IO_CG: element key, element: " EID_FMTX ": exp.key %d does not match\n",me,EID_PRTX(theElement),cge->mykey);
      assert(0);
    }

    key = KeyForObject((KEY_OBJECT *)EFATHER(theElement));
    if (key!=-1)
    {
      printf(PFMT " IO_CG: father relation, element: " EID_FMTX ", father: " EID_FMTX " ; but no father expected because element is orphan\n",
             me,EID_PRTX(theElement),EID_PRTX(EFATHER(theElement)));

      assert(0);
    }

    for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
    {
      key = KeyForObject((KEY_OBJECT *)CORNER(theElement,j));
      if (key!=-1 && key!=cge->nodekey[j])
      {
        printf(PFMT " IO_CG: corner relation, element: " EID_FMTX ", corner[%d]: " ID_FMTX " :exp.key %d does not match\n",
               me,EID_PRTX(theElement),j,ID_PRTX(CORNER(theElement,j)),cge->nodekey[j]);
        assert(0);
      }
    }

    for (j=0; j<SIDES_OF_ELEM(theElement); j++)
      CheckNeighborElement( theElement, j, cge->neighborkey[j], "IO_CG:" );
  }

  return (0);
}
#else
static INT CheckCGKeys (INT ne, ELEMENT** eid_e, MGIO_CG_ELEMENT *cg_elem)
{
  return (0);
}
#endif

MULTIGRID *LoadMultiGrid (char *MultigridName, char *name, char *type, char *BVPName, char *format, unsigned long heapSize, INT force, INT optimizedIE, INT autosave)
/* Documentation of the intended programm flow resp. communication requirements.
   Functions introducing a global communication (all processors without any exception)
   are:
        UG_GlobalMin* and DDD_XferBegin/End
   the corresponding topics in the following lists are preceeded by "#"
   The sequence of main topics is:
                read coarse grid data
                create all grid levels
                insert coarse mesh

   #CreateAlgebra				(incl. DDD_IF* and
                                                   SetSurfaceClasses incl. UG_GlobalMinINT)
                PrepareAlgebraModification

                DDD_IdentifyBegin
                        read parinfo of coarse-grid
                        read refinements and build local tree (fine grid)
                DDD_IdentifyEnd

                repair inconsistencies:
                        CommunicateEClasses		(incl. DDD_IF*)
   #IO_GridCons				(incl. ConstructConsistentGrid
                                                   incl. DDD_IF* and DDD_XferBegin/End )
                        Propagate*VectorClasses (incl. DDD_IF*)
                        PrepareAlgebraModification
   #SetSurfaceClasses		(incl. UG_GlobalMinINT)

        if you must/will leave this sequence premature ENSURE THE COMMUNICATION PATTERN!!!
        i.e. do all global communications as dummies and
                all element orientated communication (like DDD_IF*) if there exist elements
 */
{
  MULTIGRID *theMG;
  GRID *theGrid;
  ELEMENT *theElement,*ENext;
  NODE *theNode;
  VECTOR *theVector;
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
  INT i,j,*Element_corner_uniq_subdom, *Ecusdp[2],**Enusdp[2],**Ecidusdp[2],
  **Element_corner_ids_uniq_subdom,*Element_corner_ids,max,**Element_nb_uniq_subdom,
  *Element_nb_ids,id,level;
  INT *Element_SideOnBnd_uniq_subdom,*ESoBusdp[2];
  char buf[64],itype[10];
  int *vidlist;
  INT MarkKey;
#ifdef __THREEDIM__
  ELEMENT *theNeighbor;
  INT k;
#endif

  if (autosave)
  {
    if (name==NULL)
    {
      return (NULL);
    }
    else
    {
      if (type==NULL) return (NULL);
      strcpy(itype,type);
      if (name==NULL) return (NULL);
      strcpy(filename,name);
      strcat(filename,".000000");
    }
  }
  else
  {
    if (type==NULL) return (NULL);
    strcpy(itype,type);
    if (name==NULL) return (NULL);
    strcpy(filename,name);
  }
  sprintf(buf,".ug.mg.");
  strcat(filename,buf);
  strcat(filename,itype);

#ifdef ModelP
  proc_list_size = PROCLISTSIZE_VALUE;
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
    if (Read_MG_General(&mg_general))
    {
      CloseMGFile ();
      nparfiles = -1;
#ifdef ModelP
      assert(0);
#endif
    }

    nparfiles = mg_general.nparfiles;
    if (mg_general.nparfiles>procs)                                                                         {UserWrite("ERROR: too many processors needed\n");CloseMGFile (); nparfiles = -1;}
    assert(mg_general.me == me);

  }
  else if(MGIO_filetype(filename) == FT_FILE)
  {
    if (Read_OpenMGFile (filename))                                                                         {nparfiles = -1;}
    else
    if (Read_MG_General(&mg_general))
    {
      CloseMGFile ();
      nparfiles = -1;
#ifdef ModelP
      assert(0);
#endif
    }
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
  if (nparfiles == -1)
  {
    UserWriteF("ERROR in LoadMultiGrid: could not open file %s.\n", filename );
    return(NULL);
  }

  if (procs>nparfiles)
  {
    Broadcast(&mg_general,sizeof(MGIO_MG_GENERAL));
    if (me < nparfiles)
      mg_general.me = me;
  }
  if (mg_general.dim!=DIM) {
    UserWrite("ERROR: wrong dimension\n");
    CloseMGFile ();
    return (NULL);
  }
  if (strcmp(mg_general.version,MGIO_VERSION)!=0 && force==0) {
    UserWrite("ERROR: wrong version\n");
    CloseMGFile ();
    return (NULL);
  }
  /* BVP and format */
  if (BVPName==NULL) strcpy(BndValName,mg_general.DomainName);
  else strcpy(BndValName,BVPName);
  if (MultigridName==NULL) strcpy(MGName,mg_general.MultiGridName);
  else strcpy(MGName,MultigridName);
  if (format==NULL) strcpy(FormatName,mg_general.Formatname);
  else strcpy(FormatName,format);
  if (heapSize==0) heapSize = mg_general.heapsize * KBYTE;

  /* create a virginenal multigrid on the BVP */
  theMG = CreateMultiGrid(MGName,BndValName,FormatName,heapSize,TRUE,FALSE);
  if (theMG==NULL) {
    UserWrite("ERROR(ugio): cannot create multigrid\n");
    CloseMGFile ();
    return (NULL);
  }
  MG_MAGIC_COOKIE(theMG) = mg_general.magic_cookie;

  if (me >= nparfiles)
  {

    if (DisposeGrid(GRID_ON_LEVEL(theMG,0)))        {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    for (i=0; i<mg_general.nLevel; i++)
      if (CreateNewLevel(theMG,0)==NULL)              {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    /* no coarse mesh */

    if (CreateAlgebra (theMG))                                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    if (PrepareAlgebraModification(theMG))          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

        #ifdef ModelP
    DDD_IdentifyBegin();
    /* no elements to insert */
    DDD_IdentifyEnd();

    if (MGIO_PARFILE)
    {
      /* replaced by IO_GridCons because of uniformity; Christian Wrobel 980311. TODO remove:
         for (i=0; i<mg_general.nLevel; i++)
          ConstructConsistentGrid(GRID_ON_LEVEL(theMG,i));
       */
      if (IO_GridCons(theMG))                                 {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    }
                #endif

    if (SetSurfaceClasses (theMG))                          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    if (CloseMGFile ())                                             {DisposeMultiGrid(theMG); return (NULL);}

    /* saved */
    MG_SAVED(theMG) = 1;
    strcpy(MG_FILENAME(theMG),filename);

    return (theMG);

  }             /* else if (me < nparfiles) */
  theHeap = MGHEAP(theMG);
  MarkKey = MG_MARK_KEY(theMG);
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
  rr_rules = (MGIO_RR_RULE *)GetTmpMem(theHeap,rr_general.nRules*sizeof(MGIO_RR_RULE),MarkKey);
  if (rr_rules==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for rr_rules\n",(int)rr_general.nRules*sizeof(MGIO_RR_RULE)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_RR_Rules(rr_general.nRules,rr_rules))                                          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general information about coarse grid */
  if (Read_CG_General(&cg_general))                                                                       {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

#ifdef ModelP
  /* this proc has no elements */
  if (cg_general.nElement == 0)
  {             /* no elements for me; may skip element orientated communication routines */

    ReleaseTmpMem(theHeap,MarkKey);

    /* no coarse mesh */

    /* replaced by IO_GridCons because of uniformity; Christian Wrobel 980311. TODO remove:
       ConstructConsistentGrid(GRID_ON_LEVEL(theMG,0));
     */
    /* create levels */
    for (i=1; i<mg_general.nLevel; i++)
    {
      if (CreateNewLevel(theMG,0)==NULL)      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
      /* replaced by IO_GridCons because of uniformity; Christian Wrobel 980311. TODO remove:
         ConstructConsistentGrid(GRID_ON_LEVEL(theMG,i));
       */
    }

    if (CreateAlgebra (theMG))                                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    if (PrepareAlgebraModification(theMG))          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    DDD_IdentifyBegin();
    /* no elements to insert */
    DDD_IdentifyEnd();

    if (MGIO_PARFILE)
      if (IO_GridCons(theMG))                                 {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}


    if (SetSurfaceClasses (theMG))                          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    if (CloseMGFile ())                                             {DisposeMultiGrid(theMG); return (NULL);}

    /* saved */
    MG_SAVED(theMG) = 1;
    strcpy(MG_FILENAME(theMG),filename);

    return (theMG);
  }
#endif

  /* read coarse grid points and elements */
  cg_point = (MGIO_CG_POINT *)GetTmpMem(theHeap,cg_general.nPoint*sizeof(MGIO_CG_POINT),MarkKey);
  if (cg_point==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_points\n",(int)cg_general.nPoint*sizeof(MGIO_CG_POINT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_CG_Points(cg_general.nPoint,cg_point))                                         {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  if (MGIO_PARFILE)
  {
    /* read mapping: node-id --> vertex-id for orphan-nodes */
    if (Bio_Read_mint(1,&non))                                                                              {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    if (Bio_Read_mint(1,&foid))                                                                     {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    vidlist = (int*)GetTmpMem(theHeap,non*sizeof(int),MarkKey);
    if (Bio_Read_mint(non,vidlist))                                                                 {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
#ifdef Debug
    for (i=0; i<non; i++) PRINTDEBUG(gm,1,("LoadMG(): vidList[%d]=%d\n",i,vidlist[i]));
#endif
  }
  else
  {
    foid = 0;
    non = cg_general.nPoint;
  }

  cg_element = (MGIO_CG_ELEMENT *)GetTmpMem(theHeap,cg_general.nElement*sizeof(MGIO_CG_ELEMENT),MarkKey);
  if (cg_element==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for cg_elements\n",(int)cg_general.nElement*sizeof(MGIO_CG_ELEMENT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_CG_Elements(cg_general.nElement,cg_element))                           {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general bnd information */
  if (Bio_Jump (0))                                                                                                       {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_BD_General (&bd_general))                                                                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read bnd points */
  if (bd_general.nBndP > 0) {
    BndPList = (BNDP**)GetTmpMem(theHeap,bd_general.nBndP*sizeof(BNDP*),MarkKey);
    if (BndPList==NULL) {
      UserWriteF("ERROR: cannot allocate %d bytes for BndPList\n",
                 (int)bd_general.nBndP*sizeof(BNDP*));
      CloseMGFile ();
      DisposeMultiGrid(theMG);
      return (NULL);
    }
    if (Read_PBndDesc (theBVP,theHeap,bd_general.nBndP,BndPList)) {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }

  /* create and insert coarse mesh: prepare */
  theMesh.nBndP = cg_general.nBndPoint;
  theMesh.theBndPs = BndPList;
  theMesh.nInnP = cg_general.nInnerPoint;
  if (cg_general.nInnerPoint>0)
  {
    theMesh.Position = (DOUBLE**)GetTmpMem(theHeap,cg_general.nInnerPoint*sizeof(DOUBLE*),MarkKey);
    if (theMesh.Position==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.Position\n",(int)cg_general.nInnerPoint*sizeof(DOUBLE*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    Positions = (DOUBLE*)GetTmpMem(theHeap,MGIO_DIM*cg_general.nInnerPoint*sizeof(DOUBLE),MarkKey);
    if (Positions==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Positions\n",(int)MGIO_DIM*cg_general.nInnerPoint*sizeof(DOUBLE)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }
  if (MGIO_PARFILE)
  {
    theMesh.VertexLevel = (char*)GetTmpMem(theHeap,(cg_general.nBndPoint+cg_general.nInnerPoint)*sizeof(char),MarkKey);
    if (theMesh.VertexLevel==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for VertexLevel\n",(int)(cg_general.nBndPoint+cg_general.nInnerPoint)*sizeof(char)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }
  else
    theMesh.VertexLevel = NULL;
  theMesh.nSubDomains = theBVPDesc.nSubDomains;
  theMesh.nElements = (INT*)GetTmpMem(theHeap,(theMesh.nSubDomains+1)*sizeof(INT),MarkKey);
  if (theMesh.nElements==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.nElements\n",(int)(theMesh.nSubDomains+1)*sizeof(INT)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  theMesh.ElementLevel = (char**)GetTmpMem(theHeap,(theMesh.nSubDomains+1)*sizeof(char*),MarkKey);
  if (theMesh.ElementLevel==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for theMesh.ElementLevel\n",(int)(theMesh.nSubDomains+1)*sizeof(char*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<=theMesh.nSubDomains; i++)
  {
    theMesh.nElements[i] = 0;
    theMesh.ElementLevel[i] = NULL;
  }
  theMesh.nElements[1] = cg_general.nElement;
  theMesh.ElementLevel[1] = (char*)GetTmpMem(theHeap,cg_general.nElement*sizeof(char),MarkKey);
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
  Element_corner_uniq_subdom  = (INT*)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT),MarkKey);
  if (Element_corner_uniq_subdom==NULL)
  {
    UserWriteF("ERROR: cannot allocate %d bytes for Element_corner_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT));
    CloseMGFile ();
    DisposeMultiGrid(theMG);
    return (NULL);
  }
  Element_SideOnBnd_uniq_subdom  = (INT*)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT),MarkKey);
  if (Element_SideOnBnd_uniq_subdom==NULL)
  {
    UserWriteF("ERROR: cannot allocate %d bytes for Element_SideOnBnd_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT));
    CloseMGFile ();
    DisposeMultiGrid(theMG);
    return (NULL);
  }
  max = 0;
  for (i=0; i<cg_general.nElement; i++)
  {
    cge = MGIO_CG_ELEMENT_PS(cg_element,i);
    Element_corner_uniq_subdom[i] = ge_element[cge->ge].nCorner;
    max = MAX(max,ge_element[cge->ge].nCorner);
    max = MAX(max,ge_element[cge->ge].nSide);
    if (MGIO_PARFILE) theMesh.ElementLevel[1][i] = cge->level;
    else theMesh.ElementLevel[1][i] = 0;
    Element_SideOnBnd_uniq_subdom[i] = cge->side_on_bnd;
  }
  Ecusdp[1] = Element_corner_uniq_subdom;
  theMesh.Element_corners = Ecusdp;
  ESoBusdp[1] = Element_SideOnBnd_uniq_subdom;
  theMesh.ElemSideOnBnd = ESoBusdp;

  /* corners ids of elements */
  Element_corner_ids_uniq_subdom  = (INT**)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT*),MarkKey);
  if (Element_corner_ids_uniq_subdom==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_corner_ids_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  Element_corner_ids  = (INT*)GetTmpMem(theHeap,max*cg_general.nElement*sizeof(INT),MarkKey);
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
        PRINTDEBUG(gm,1,("LoadMultiGrid(): cg_elem=%d  cg_nid[%d]=%d foid=%d\n",i,j,cge->cornerid[j],foid));
      }
    else
      for (j=0; j<Element_corner_uniq_subdom[i]; j++)
        Element_corner_ids_uniq_subdom[i][j] = cge->cornerid[j];
  }
  Ecidusdp[1] = Element_corner_ids_uniq_subdom;
  theMesh.Element_corner_ids = Ecidusdp;

  /* nb of elements */
  Element_nb_uniq_subdom  = (INT**)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT*),MarkKey);
  if (Element_nb_uniq_subdom==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for Element_nb_uniq_subdom\n",(int)cg_general.nElement*sizeof(INT*)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  Element_nb_ids  = (INT*)GetTmpMem(theHeap,max*cg_general.nElement*sizeof(INT),MarkKey);
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
        theNode = CORNER(theElement,j);
        if (OBJT(MYVERTEX(theNode))==BVOBJ) SETNSUBDOM(theNode,0);
        else SETNSUBDOM(theNode,cge->subdomain);
        ID(CORNER(theElement,j)) = cge->cornerid[j];
      }
      for (j=0; j<EDGES_OF_ELEM(theElement); j++)
      {
        SETEDSUBDOM(GetEdge(CORNER_OF_EDGE_PTR(theElement,j,0),CORNER_OF_EDGE_PTR(theElement,j,1)),cge->subdomain);
      }
    }
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement = PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
      if (OBJT(theElement)==BEOBJ)
#ifdef __TWODIM__
        for (j=0; j<EDGES_OF_ELEM(theElement); j++)
          if (EDGE_ON_BND(theElement,j))
            SETEDSUBDOM(GetEdge(CORNER_OF_EDGE_PTR(theElement,j,0),CORNER_OF_EDGE_PTR(theElement,j,1)),0);

#endif
#ifdef __THREEDIM__
        for (j=0; j<SIDES_OF_ELEM(theElement); j++)
          if (SIDE_ON_BND(theElement,j))
            for (k=0; k<EDGES_OF_SIDE(theElement,j); k++)
              SETEDSUBDOM(GetEdge(CORNER_OF_EDGE_PTR(theElement,EDGE_OF_SIDE(theElement,j,k),0),CORNER_OF_EDGE_PTR(theElement,EDGE_OF_SIDE(theElement,j,k),1)),0);

#endif
        /* now CreateAlgebra  is necessary to have the coarse grid nodevectors for DDD identification in Evaluate_pinfo */
        if (CreateAlgebra (theMG))                                                                                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (PrepareAlgebraModification(theMG))                                                          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

  i = MG_ELEMUSED | MG_NODEUSED | MG_EDGEUSED | MG_VERTEXUSED |  MG_VECTORUSED;
  ClearMultiGridUsedFlags(theMG,0,TOPLEVEL(theMG),i);

  /* open identification context */
  DDD_IdentifyBegin();

  /* read parinfo of coarse-grid */
  if (MGIO_PARFILE)
  {
    ActProcListPos = ProcList = (unsigned short*)GetTmpMem(theHeap,PROCLISTSIZE*sizeof(unsigned short),MarkKey);
    if (ProcList==NULL)     {UserWriteF("ERROR: cannot allocate %d bytes for ProcList\n",(int)PROCLISTSIZE*sizeof(unsigned short)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    cg_pinfo.proclist = ProcList;
    for (level=0; level<=TOPLEVEL(theMG); level++)
    {
      theGrid = GRID_ON_LEVEL(theMG,level);
      for (theElement = PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=ENext)
      {
        /* store succe because Evaluate_pinfo() relinks the element list */
        ENext = SUCCE(theElement);
        cge = MGIO_CG_ELEMENT_PS(cg_element,ID(theElement));
        if (Read_pinfo (TAG(theElement),&cg_pinfo))     {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
        if (Evaluate_pinfo(theGrid,theElement,&cg_pinfo))       {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
      }
    }
  }

  /* are we ready ? */
  if (mg_general.nLevel==1)
  {             /* if we are here now, all other processors (wich haven't returned yet) are here TOO
                   because mg_general.nLevel is a global quantity. */

    /* no fine grid elements */

    /* close identification context */
    DDD_IdentifyEnd();

    /* repair inconsistencies */
    if (MGIO_PARFILE)
      if (IO_GridCons(theMG))                                 {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    for (theVector=PFIRSTVECTOR(GRID_ON_LEVEL(theMG,0));
         theVector!=NULL; theVector=SUCCVC(theVector)) {
      SETVCLASS(theVector,3);
      SETVNCLASS(theVector,0);
      SETNEW_DEFECT(theVector,1);
      SETFINE_GRID_DOF(theVector,1);
    }

    if (SetSurfaceClasses (theMG))                          {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}

    ReleaseTmpMem(theHeap,MarkKey);
    if (CloseMGFile ())                                             {DisposeMultiGrid(theMG); return (NULL);}

    /* saved */
    MG_SAVED(theMG) = 1;
    strcpy(MG_FILENAME(theMG),filename);

    return (theMG);
  }

  /* list: node-id --> node */
  nid_n = (NODE**)GetTmpMem(theHeap,non*sizeof(NODE*),MarkKey);
  id = foid;
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theNode=PFIRSTNODE(GRID_ON_LEVEL(theMG,i)); theNode!=NULL; theNode=SUCCN(theNode))
    {
      assert(ID(theNode)-foid<non);
      nid_n[ID(theNode)-foid] = theNode;
    }

  /* list: elem-id --> elem */
  eid_e = (ELEMENT**)GetTmpMem(theHeap,cg_general.nElement*sizeof(ELEMENT*),MarkKey);
  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      assert(ID(theElement)>=0);
      assert(ID(theElement)<cg_general.nElement);
      eid_e[ID(theElement)] = theElement;
    }

  if (CheckCGKeys(cg_general.nElement, eid_e, cg_element)==0)
  {
    PRINTDEBUG(gm,4,("CheckCGKeys ok.\n"));
  }
  else
  {
    PrintErrorMessage('E',"LoadMultiGrid","CheckCGKeys failed\n");
  }

  /* read hierarchical elements */
  refinement = (MGIO_REFINEMENT*)malloc(MGIO_REFINEMENT_SIZE);
  if (refinement==NULL) {UserWriteF("ERROR: cannot allocate %d bytes for refinement\n",(int)MGIO_REFINEMENT_SIZE); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (MGIO_PARFILE)
  {
    ProcList = (unsigned short*)malloc(PROCLISTSIZE*sizeof(unsigned short));
    if (ProcList==NULL)     {UserWriteF("ERROR: cannot allocate %d bytes for ProcList\n",(int)PROCLISTSIZE*sizeof(unsigned short)); CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
    for (i=0; i<MAX_SONS; i++) refinement->pinfo[i].proclist = ProcList+i*ELEMPROCLISTSIZE;
  }
  for (j=0; j<cg_general.nElement; j++)
  {
    theElement = eid_e[j];
    theGrid = GRID_ON_LEVEL(theMG,LEVEL(theElement));
    cge = MGIO_CG_ELEMENT_PS(cg_element,j);
    if (cge->nref==0)
    {
      SETREFINE(theElement,NO_REFINEMENT);
      SETREFINECLASS(theElement,NO_CLASS);
      SETMARK(theElement,NO_REFINEMENT);
      SETMARKCLASS(theElement,NO_CLASS);
      if (LEVEL(theElement)==0) SETECLASS(theElement,RED_CLASS);
#ifdef ModelP
      else assert(EGHOST(theElement));                          /* masters elements must have a father or be on level 0 */
#else
      else assert(0);                           /* all orphans must be on level 0 in sequential mode */
#endif
      continue;
    }
    if (InsertLocalTree(theGrid,theElement,refinement,rr_general.RefRuleOffset)) {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }

  /* close identification context */
  DDD_IdentifyEnd();

  /* repair inconsistencies */
  if (MGIO_PARFILE)
  {
#ifdef ModelP
    CommunicateEClasses(theMG);
#endif
    if (IO_GridCons(theMG))                         {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }

  /* postprocess */
  for (i=0; i<TOPLEVEL(theMG); i++)
    for (theElement = FIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      SETMARK(theElement,0);
      SETMARKCLASS(theElement,NO_CLASS);
      SETEBUILDCON(theElement,1);
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
          if (k==SIDES_OF_ELEM(theNeighbor))      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
          if (DisposeDoubledSideVector (theGrid,theElement,j,theNeighbor,k)) {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
        }
#endif
    if (GridCreateConnection(theGrid))                      {CloseMGFile (); DisposeMultiGrid(theMG); return (NULL);}
  }
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
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
  if (PrepareAlgebraModification(theMG))                                                                                          {DisposeMultiGrid(theMG); return (NULL);}

  /* set DOFs on vectors */
  if (SetSurfaceClasses (theMG))                                                                                                          {DisposeMultiGrid(theMG); return (NULL);}

  /* close file */
  ReleaseTmpMem(theHeap,MarkKey);
  if (CloseMGFile ())                                                                                                                             {DisposeMultiGrid(theMG); return (NULL);}


  /* saved */
  MG_SAVED(theMG) = 1;
  strcpy(MG_FILENAME(theMG),filename);

        #if EXTRACT_RULES
  ResetRefineTagsBeyondRuleManager(theMG);
        #endif

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
    REP_ERR_RETURN(1);
  }

  stream = fopen(docName,"w");
  if (stream==NULL)
  {
    PrintErrorMessage('E',"SaveCnomGridAndValues","can't open file");
    REP_ERR_RETURN(1);
  }

  if (PlotProcInfo->PreprocessProc!=NULL)
    if ((*PlotProcInfo->PreprocessProc)(NULL,theMG)!=0)
      REP_ERR_RETURN(1);

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
          val=(*PlotProcInfo->EvalProc)(theElement, (const DOUBLE **) CoordOfCornerPtr, LOCAL_COORD_OF_ELEM(theElement,i));
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
          val=(*PlotProcInfo->EvalProc)(theElement, (const DOUBLE **) CoordOfCornerPtr,  LOCAL_COORD_OF_ELEM(theElement,i));
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

  if (MGIO_Init ()) return(1);

  return (0);
}
