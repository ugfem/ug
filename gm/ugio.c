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

/* include refine because of macros accessed  */
#include "refine.h"

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

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

static INT RenumberNE (MULTIGRID *theMG)
{
  INT i,nid,eid;
  GRID *theGrid;
  NODE *theNode;
  ELEMENT *theElement;

  nid=eid=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      ID(theElement) = eid++;
    if (i==0)
    {
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        if (OBJT(MYVERTEX(theNode))==BVOBJ)
          ID(theNode) = nid++;
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        if (OBJT(MYVERTEX(theNode))==IVOBJ)
          ID(theNode) = nid++;
    }
    else
      for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        ID(theNode) = nid++;
  }

  return (0);
}

INT SaveMultiGrid_SCR (MULTIGRID *theMG, char *name, char *comment)
{
  FILE *stream;
  GRID *theGrid;
  NODE *theNode;
  ELEMENT *theElement;
  VERTEX *theVertex;
  COORD *global;
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

  if (TOPLEVEL(theMG) > 0)
    PrintErrorMessage('W',"SaveMultiGrid",
                      "only level 0 will be saved");

  /* get BVPDesc */
  if (BVP_SetBVPDesc(MG_BVP(theMG),&theBVPDesc))
    RETURN (GM_ERROR);

  /* get time */
  fmt = "%a %b %d %H:%M:%S %Y";
  time(&Time);
  strftime(buffer,BUFFERSIZE,fmt,localtime(&Time));

  /* write header */
  fprintf(stream,HEADER_FMT,BVPD_NAME(theBVPDesc),buffer,name,comment);

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

INT SaveMultiGrid_SPF (MULTIGRID *theMG, char *name, char *comment)
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
  MGIO_RR_GENERAL rr_general;
  MGIO_CG_GENERAL cg_general;
  MGIO_CG_POINT *cg_point;
  MGIO_CG_ELEMENT *cg_element;
  MGIO_BD_GENERAL bd_general;
  INT i,j,k,niv,nbv,nie,nbe,n,bvi,ivi;
  char *p;
  BNDP **BndPList;

  /* check */
  if (theMG==NULL) return (1);

  /* open file */
  if (Write_OpenFile (name)) return (1);

  /* write general information */
  theBVP = MG_BVP(theMG);
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc)) return (1);
  theGrid = GRID_ON_LEVEL(theMG,0);
  p = name + strlen(name) - 4;
  if (strcmp(p,".asc")==0) mg_general.mode = 0;
  else if (strcmp(p,".bin")==0) mg_general.mode = 1;
  else return (1);
  mg_general.dim                  = DIM;
  mg_general.nLevel               = 1;
  mg_general.nNode                = NN(theGrid);
  mg_general.nPoint               = NV(theGrid);
  mg_general.nElement             = NT(theGrid);
  mg_general.nHierElem    = 0;
  strcpy(mg_general.DomainName,BVPD_NAME(theBVPDesc));
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
  rr_general.nRules = 0;
  if (Write_RR_General(&rr_general)) return (1);

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
  if (RenumberNE (theMG)) return (1);
  theGrid = GRID_ON_LEVEL(theMG,0);
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap);
  n = NV(theGrid)*sizeof(MGIO_CG_POINT);
  cg_point = (MGIO_CG_POINT *)GetTmpMem(theHeap,n);
  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    for (j=0; j<MGIO_DIM; j++)
      cg_point[ID(theNode)].position[j] = CVECT(MYVERTEX(theNode))[j];

  if (Write_CG_Points((int)NV(theGrid),cg_point)) return (1);
  ReleaseTmpMem(theHeap);

  /* write coarse grid elements */
  theGrid = GRID_ON_LEVEL(theMG,0);
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap);
  n = NT(theGrid);
  cg_element = (MGIO_CG_ELEMENT*)GetTmpMem(theHeap,n*sizeof(MGIO_CG_ELEMENT));
  for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    i = ID(theElement);

    cg_element[i].ge = TAG(theElement);
    for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      cg_element[i].cornerid[j] = ID(CORNER(theElement,j));
    for (j=0; j<SIDES_OF_ELEM(theElement); j++)
      if (NBELEM(theElement,j)!=NULL)
        cg_element[i].nbid[j] = ID(NBELEM(theElement,j));
      else
        cg_element[i].nbid[j] = -1;
    cg_element[i].refrule = -1;
    cg_element[i].nmoved = 0;
  }
  if (Write_CG_Elements((int)n,cg_element)) return (1);
  ReleaseTmpMem(theHeap);

  /* write general bnd information */
  bd_general.nBndP = nbv;
  if (Write_BD_General (&bd_general)) return (1);

  /* write bnd information */
  MarkTmpMem(theHeap);
  BndPList = (BNDP**)GetTmpMem(theHeap,nbv*sizeof(BNDP*));
  for (theNode=FIRSTNODE(GRID_ON_LEVEL(theMG,0)); theNode!=NULL; theNode=SUCCN(theNode))
    if (OBJT(MYVERTEX(theNode))==BVOBJ)
    {
      if (ID(theNode) < 0 || ID(theNode) >= nbv) return (1);
      BndPList[ID(theNode)] = V_BNDP(MYVERTEX(theNode));
    }
  if (Write_PBndDesc (nbv,BndPList)) return (1);

  /* release tmp mem */
  ReleaseTmpMem(theHeap);

  /* close file */
  if (CloseFile ()) return (1);

  return (0);
}

INT SaveMultiGrid (MULTIGRID *theMG, char *name, char *comment)
{
  char *p;
  INT mode;

  /* check name convention */
  p = name + strlen(name) - 4;
  if (strcmp(p,".scr")==0)
  {
    if (SaveMultiGrid_SCR (theMG,name,comment)) return (1);
  }
  else if (strcmp(p,".asc")==0 || strcmp(p,".bin")==0)
  {
    if (SaveMultiGrid_SPF (theMG,name,comment)) return (1);
  }
  else
  {
    UserWrite("filename must end with '.scr', '.asc' or '.bin'\n");
    return (1);
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

MULTIGRID *LoadMultiGrid (char *MultigridName, char *FileName, char *BVPName, char *format, unsigned long heapSize)
{
  MULTIGRID *theMG;
  HEAP *theHeap;
  MGIO_MG_GENERAL mg_general;
  MGIO_GE_GENERAL ge_general;
  MGIO_GE_ELEMENT ge_element[TAGS];
  MGIO_RR_GENERAL rr_general;
  MGIO_CG_GENERAL cg_general;
  MGIO_CG_POINT *cg_point;
  MGIO_CG_ELEMENT *cg_element;
  MGIO_BD_GENERAL bd_general;
  BNDP **BndPList;
  COORD **PositionList, *Positions;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  MESH theMesh;
  char FormatName[NAMESIZE], BndValName[NAMESIZE];
  INT i,j,*Element_corner_uniq_subdom, *Ecusdp[2],**Enusdp[2],**Ecidusdp[2],**Element_corner_ids_uniq_subdom,*Element_corner_ids,max,**Element_nb_uniq_subdom,*Element_nb_ids;

  /* open file */
  if (Read_OpenFile (FileName))                                                                           {return (NULL);}
  if (Read_MG_General(&mg_general))                                                                       {CloseFile (); return (NULL);}
  if (mg_general.dim!=DIM)                                                                                        {UserWrite("ERROR: wrong dimension\n");CloseFile (); return (NULL);}

  /* BVP and format */
  if (BVPName==NULL) strcpy(BndValName,mg_general.DomainName);
  else strcpy(BndValName,BVPName);
  if (format==NULL) strcpy(FormatName,mg_general.Formatname);
  else strcpy(FormatName,format);

  /* create a virginenal multigrid on the BVP */
  theMG = CreateMultiGrid(MultigridName,BndValName,FormatName,heapSize);
  if (theMG==NULL)                                                                                                        {CloseFile (); return (NULL);}
  if (DisposeGrid(GRID_ON_LEVEL(theMG,0)))                                                        {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (CreateNewLevel(theMG)==NULL)                                                                        {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  theHeap = MGHEAP(theMG);
  theBVP = MG_BVP(theMG);
  if (theBVP==NULL)                                                                                                       {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (BVP_SetBVPDesc(theBVP,&theBVPDesc))                                                         {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  MarkTmpMem(theHeap);

  /* read general element information */
  if (Read_GE_General(&ge_general))                                                                       {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_GE_Elements(TAGS,ge_element))                                                          {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general refrule information */
  if (Read_RR_General(&rr_general))                                                                   {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general information about coarse grid */
  if (Read_CG_General(&cg_general))                                                                       {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read coarse grid points and elements */
  cg_point = (MGIO_CG_POINT *)GetTmpMem(theHeap,cg_general.nPoint*sizeof(MGIO_CG_POINT));
  cg_element = (MGIO_CG_ELEMENT *)GetTmpMem(theHeap,cg_general.nElement*sizeof(MGIO_CG_ELEMENT));
  if (Read_CG_Points(cg_general.nPoint,cg_point))                                         {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_CG_Elements(cg_general.nElement,cg_element))                           {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read general bnd information */
  if (Read_BD_General (&bd_general))                                                                      {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* read bnd points */
  BndPList = (BNDP**)GetTmpMem(theHeap,bd_general.nBndP*sizeof(BNDP*));
  if (BndPList==NULL)                                                                                                     {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  if (Read_PBndDesc (theBVP,theHeap,bd_general.nBndP,BndPList))           {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* create and insert coarse mesh: prepare */
  theMesh.nBndP = cg_general.nBndPoint;
  theMesh.theBndPs = BndPList;
  theMesh.nInnP = cg_general.nInnerPoint;
  theMesh.Position = (COORD**)GetTmpMem(theHeap,cg_general.nInnerPoint*sizeof(COORD*));
  Positions = (COORD*)GetTmpMem(theHeap,MGIO_DIM*cg_general.nInnerPoint*sizeof(COORD));
  for (i=0; i<cg_general.nInnerPoint; i++)
    theMesh.Position[i] = Positions+MGIO_DIM*i;
  for (i=0; i<cg_general.nInnerPoint; i++)
    for (j=0; j<MGIO_DIM; j++)
      theMesh.Position[i][j] = cg_point[cg_general.nBndPoint+i].position[j];
  theMesh.nSubDomains = theBVPDesc.nSubDomains;
  theMesh.nElements = (INT*)GetTmpMem(theHeap,(theMesh.nSubDomains+1)*sizeof(INT));
  if (theMesh.nElements==NULL)                                                                            {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<=theMesh.nSubDomains; i++) theMesh.nElements[i] = 0;theMesh.nElements[1] = cg_general.nElement;

  /* nb of corners of elements */
  Element_corner_uniq_subdom  = (INT*)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT));
  if (Element_corner_uniq_subdom==NULL)                                                           {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
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
  if (Element_corner_ids_uniq_subdom==NULL)                                                       {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  Element_corner_ids  = (INT*)GetTmpMem(theHeap,max*cg_general.nElement*sizeof(INT));
  if (Element_corner_ids==NULL)                                                                           {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<cg_general.nElement; i++)
    Element_corner_ids_uniq_subdom[i] = Element_corner_ids+max*i;
  for (i=0; i<cg_general.nElement; i++)
    for (j=0; j<Element_corner_uniq_subdom[i]; j++)
      Element_corner_ids_uniq_subdom[i][j] = cg_element[i].cornerid[j];
  Ecidusdp[1] = Element_corner_ids_uniq_subdom;
  theMesh.Element_corner_ids = Ecidusdp;

  /* nb of elements */
  Element_nb_uniq_subdom  = (INT**)GetTmpMem(theHeap,cg_general.nElement*sizeof(INT*));
  if (Element_nb_uniq_subdom==NULL)                                                                       {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  Element_nb_ids  = (INT*)GetTmpMem(theHeap,max*cg_general.nElement*sizeof(INT));
  if (Element_nb_ids==NULL)                                                                                       {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}
  for (i=0; i<cg_general.nElement; i++)
    Element_nb_uniq_subdom[i] = Element_nb_ids+max*i;
  for (i=0; i<cg_general.nElement; i++)
    for (j=0; j<ge_element[cg_element[i].ge].nSide; j++)
      Element_nb_uniq_subdom[i][j] = cg_element[i].nbid[j];
  Enusdp[1] = Element_nb_uniq_subdom;
  theMesh.nbElements = Enusdp;

  /* insert coarse mesh */
  if (InsertMesh(theMG,&theMesh))                                                                         {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

  /* release tmp mem */
  ReleaseTmpMem(theHeap);

  /* close file */
  if (CloseFile ())                                                                                                       {CloseFile (); DisposeMultiGrid(theMG); return (NULL);}

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

static COORD LocalCoord[2][4][2]=
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
  COORD *CoordOfCornerPtr[8];
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
    theGrid = theMG->grids[k];
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
    theGrid = theMG->grids[k];
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
      {
        CORNER_COORDINATES(theElement,n,CoordOfCornerPtr);
        for (i=0; i<n; i++)
        {
          val=(*PlotProcInfo->EvalProc)(theElement, (const COORD **) CoordOfCornerPtr, (COORD *) &(LocalCoord[n-3,i,0]));
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
    theGrid = theMG->grids[k];
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
    theGrid = theMG->grids[k];
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
    theGrid = theMG->grids[k];
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
    theGrid = theMG->grids[k];
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      if ((k==j)||LEAFELEM(theElement))
      {
        CORNER_COORDINATES(theElement,n,CoordOfCornerPtr);
        for (i=0; i<n; i++)
        {
          theVertex=MYVERTEX(CORNER(theElement,i));
          if (USED(theVertex)) continue;
          val=(*PlotProcInfo->EvalProc)(theElement, (const COORD **) CoordOfCornerPtr, (COORD *) &(LocalCoord[n-3,i,0]));
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

  return (0);
}
