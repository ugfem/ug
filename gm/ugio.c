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

INT SaveMultiGrid (MULTIGRID *theMG, char *name, char *comment)
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

MULTIGRID *LoadMultiGrid (char *MultigridName, char *FileName, char *BVPName,
                          char *format, unsigned long heapSize)
{

  PrintErrorMessage('E',"LoadMultiGrid","not supported any more");
  return(NULL);
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
