// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  mgio.c														*/
/*																			*/
/* Purpose:   input/output of loc ref mg	                                                                */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   18.11.96 begin,												*/
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

#include <stdio.h>
#include "mgio.h"

#ifdef __MGIO_USE_IN_UG__

        #include "defaults.h"
        #include "fileopen.h"
        #include "domain.h"

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

#define MGIO_TITLE_LINE                 "#### sparse mg storage format ####\n"
#define MGIO_DATA_LINE                  "######### data following #########\n"
#define MGIO_NEW_LINE                   "\n"
#define MGIO_SEP_LINE                   "##################################\n"

#define MGIO_ASCII                                      0
#define MGIO_BIN                                        1

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef int (*RW_mint_proc)(int n, int *intList);
typedef int (*RW_mfloat_proc)(int n, float *floatList);
typedef int (*RW_string_proc)(char *string);

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

static FILE *stream;                    /* file                                                 */
static int gridpathes_set;              /* pathes used in ug			*/
static int rw_mode;                             /* ASCII or ... (see header)	*/
static char buffer[1024];               /* general purpose buffer		*/
static int intList[100];                /* general purpose integer list */
static float floatList[100];    /* general purpose float list   */

/* low level read/write functions */
static RW_mint_proc Read_mint, Write_mint;
static RW_mfloat_proc Read_mfloat, Write_mfloat;
static RW_string_proc Read_string, Write_string;

#ifdef __MGIO_USE_IN_UG__

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

#endif

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*************                                              *****************/
/*************  low level io functions                      *****************/
/*************                                              *****************/
/****************************************************************************/

static int ASCII_Read_mint (int n, int *intList)
{
  int i;

  for (i=0; i<n; i++)
    if (fscanf(stream,"%d",intList+i)!=1) return (1);
  return (0);
}

static int ASCII_Write_mint (int n, int *intList)
{
  int i;

  for (i=0; i<n; i++)
    if (fprintf(stream,"%d",intList[i])<0) return (1);
  return (0);
}

static int ASCII_Read_mfloat (int n, float *floatList)
{
  int i;
  float fValue;

  for (i=0; i<n; i++)
    if (fscanf(stream,"%f",floatList+i)!=1) return (1);
  return (0);
}

static int ASCII_Write_mfloat (int n, float *floatList)
{
  int i;
  float fValue;

  for (i=0; i<n; i++)
    if (fprintf(stream,"%f",floatList+i)<0) return (1);
  return (0);
}

static int ASCII_Read_string (char *string)
{
  if (fscanf(stream,"%s",string)!=1) return (1);
  return (0);
}

static int ASCII_Write_string (char *string)
{
  if (fprintf(stream,"%s",string)<0) return (1);
  return (0);
}

/****************************************************************************/
/*D
   Read_OpenFile - opens file for reading

   SYNOPSIS:
   int Read_OpenFile (char *filename);

   PARAMETERS:
   .  filename - name of file

   DESCRIPTION:
   opens a file with specified name for reading

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 when error occured.

   SEE ALSO:
   D*/
/****************************************************************************/

int Read_OpenFile (char *filename)
{

#ifdef __MGIO_USE_IN_UG__
  if (gridpathes_set) stream = FileOpenUsingSearchPaths(filename,"r","gridpathes");
  else stream = fileopen(filename,"r");
#else
  stream = fileopen(filename,"r");
#endif

  if (stream==NULL) return (1);

  return (0);
}

/****************************************************************************/
/*D
   Write_OpenFile - opens file for reading

   SYNOPSIS:
   int Write_OpenFile (char *filename);

   PARAMETERS:
   .  filename - name of file

   DESCRIPTION:
   opens a file with specified name for writing

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 when error occured.

   SEE ALSO:
   D*/
/****************************************************************************/

int Write_OpenFile (char *filename)
{

#ifdef __MGIO_USE_IN_UG__
  if (gridpathes_set) stream = FileOpenUsingSearchPaths(filename,"w","gridpathes");
  else stream = fileopen(filename,"w");
#else
  stream = fileopen(filename,"w");
#endif

  if (stream==NULL) return (1);

  return (0);
}

/****************************************************************************/
/*
   Read_MG_General - reads general information about mg

   SYNOPSIS:
   int Read_MG_General (MGIO_MG_GENERAL *mg_general);

   PARAMETERS:
   .  mg_general - general information about mg

   DESCRIPTION:
   function reads general information about the mg

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Read_MG_General (MGIO_MG_GENERAL *mg_general)
{
  if (ASCII_Read_string(buffer)) return (1);if (strcmp(buffer,MGIO_SEP_LINE)!=0) return (1);
  if (ASCII_Read_string(buffer)) return (1);if (strcmp(buffer,MGIO_TITLE_LINE)!=0) return (1);
  if (ASCII_Read_string(buffer)) return (1);if (strcmp(buffer,MGIO_SEP_LINE)!=0) return (1);
  if (ASCII_Read_mint(5,intList)) return (1);
  mg_general->mode                = intList[0];
  mg_general->nLevel              = intList[1];
  mg_general->nNode               = intList[2];
  mg_general->nPoint              = intList[3];
  mg_general->nElement    = intList[4];
  if (ASCII_Read_string(mg_general->DomainName)) return (1);
  if (ASCII_Read_string(mg_general->Formatname)) return (1);
  if (ASCII_Read_mint(1,intList)) return (1);
  mg_general->VectorTypes = intList[0];
  if (ASCII_Read_string(buffer)) return (1);if (strcmp(buffer,MGIO_NEW_LINE)!=0) return (1);
  if (ASCII_Read_string(buffer)) return (1);if (strcmp(buffer,MGIO_SEP_LINE)!=0) return (1);
  if (ASCII_Read_string(buffer)) return (1);if (strcmp(buffer,MGIO_DATA_LINE)!=0) return (1);
  if (ASCII_Read_string(buffer)) return (1);if (strcmp(buffer,MGIO_SEP_LINE)!=0) return (1);

  /* set the read functions */
  switch (mg_general->mode)
  {
  case MGIO_ASCII :
    Read_mint       = ASCII_Read_mint;
    Read_mfloat = ASCII_Read_mfloat;
    Read_string = ASCII_Read_string;
    break;
  case MGIO_BIN :
    Read_mint       = ASCII_Read_mint;
    Read_mfloat = ASCII_Read_mfloat;
    Read_string = ASCII_Read_string;
    break;
  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*
   Write_MG_General - writes general information about mg

   SYNOPSIS:
   int Write_MG_General (MGIO_MG_GENERAL *mg_general);

   PARAMETERS:
   .  mg_general - general information about mg

   DESCRIPTION:
   function writes general information about the mg

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Write_MG_General (MGIO_MG_GENERAL *mg_general)
{
  if (ASCII_Write_string(MGIO_SEP_LINE)) return (1);
  if (ASCII_Write_string(MGIO_TITLE_LINE)) return (1);
  if (ASCII_Write_string(MGIO_SEP_LINE)) return (1);
  intList[0] = mg_general->mode;
  intList[1] = mg_general->nLevel;
  intList[2] = mg_general->nNode;
  intList[3] = mg_general->nPoint;
  intList[4] = mg_general->nElement;
  if (ASCII_Write_mint(5,intList)) return (1);
  if (ASCII_Write_string(mg_general->DomainName)) return (1);
  if (ASCII_Write_string(mg_general->Formatname)) return (1);
  intList[0] = mg_general->VectorTypes;
  if (ASCII_Write_mint(1,intList)) return (1);
  if (ASCII_Write_string(MGIO_NEW_LINE)) return (1);
  if (ASCII_Write_string(MGIO_SEP_LINE)) return (1);
  if (ASCII_Write_string(MGIO_DATA_LINE)) return (1);
  if (ASCII_Write_string(MGIO_SEP_LINE)) return (1);

  /* set the write functions */
  switch (mg_general->mode)
  {
  case MGIO_ASCII :
    Write_mint              = ASCII_Write_mint;
    Write_mfloat    = ASCII_Write_mfloat;
    Write_string    = ASCII_Write_string;
    break;
  case MGIO_BIN :
    Write_mint              = ASCII_Write_mint;
    Write_mfloat    = ASCII_Write_mfloat;
    Write_string    = ASCII_Write_string;
    break;
  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*
   Read_GE_General - reads general information about general elements

   SYNOPSIS:
   int	Read_GE_General (MGIO_GE_GENERAL *ge_general);

   PARAMETERS:
   .  ge_general - general information about general elements

   DESCRIPTION:
   function reads general information about general elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Read_GE_General (MGIO_GE_GENERAL *ge_general)
{
  if ((*Read_mint)(1,intList)) return (1);
  ge_general->nGenElement = intList[0];

  return (0);
}

/****************************************************************************/
/*
   Write_GE_General - writes general information about general elements

   SYNOPSIS:
   int	Write_GE_General (MGIO_GE_GENERAL *ge_general);

   PARAMETERS:
   .  ge_general - general information about general elements

   DESCRIPTION:
   function writes general information about general elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Write_GE_General (MGIO_GE_GENERAL *ge_general)
{
  intList[0] = ge_general->nGenElement;
  if ((*Write_mint)(1,intList)) return (1);

  return (0);
}

/****************************************************************************/
/*
   Read_GE_Elements - reads general elements

   SYNOPSIS:
   int	Read_GE_Elements (int n, MGIO_GE_ELEMENT *ge_element);

   PARAMETERS:
   .  n - nb of general elements to read
   .  ge_element - first general elements

   DESCRIPTION:
   function reads general elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Read_GE_Elements (int n, MGIO_GE_ELEMENT *ge_element)
{
  int i,j,s,m;
  MGIO_GE_ELEMENT *pge;

  pge = ge_element;
  m = 4 + 6*MGIO_MAX_EDGES_OF_ELEM;
  for (i=0; i<n; i++)
  {
    if ((*Read_mint)(m,intList)) return (1);
    s=0;
    pge->tag                = intList[s++];
    pge->nCorner    = intList[s++];
    pge->nEdge              = intList[s++];
    pge->nSide              = intList[s++];
    for (j=0; j<MGIO_MAX_EDGES_OF_ELEM; j++)
    {
      pge->CornerOfEdge[j][0] = intList[s++];
      pge->CornerOfEdge[j][1] = intList[s++];
    }
    for (j=0; j<MGIO_MAX_EDGES_OF_ELEM; j++)
    {
      pge->CornerOfSide[j][0] = intList[s++];
      pge->CornerOfSide[j][1] = intList[s++];
      pge->CornerOfSide[j][2] = intList[s++];
      pge->CornerOfSide[j][3] = intList[s++];
    }
    pge++;
  }

  return (0);
}

/****************************************************************************/
/*
   Write_GE_Elements - writes general elements

   SYNOPSIS:
   int	Write_GE_Elements (int n, MGIO_GE_ELEMENT *ge_element);

   PARAMETERS:
   .  n - nb of general elements to write
   .  ge_element - first general elements

   DESCRIPTION:
   function writes general elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Write_GE_Elements (int n, MGIO_GE_ELEMENT *ge_element)
{
  int i,j,s;
  MGIO_GE_ELEMENT *pge;

  pge = ge_element;
  for (i=0; i<n; i++)
  {
    s=0;
    intList[s++] = pge->tag;
    intList[s++] = pge->nCorner;
    intList[s++] = pge->nEdge;
    intList[s++] = pge->nSide;
    for (j=0; j<MGIO_MAX_EDGES_OF_ELEM; j++)
    {
      intList[s++] = pge->CornerOfEdge[j][0];
      intList[s++] = pge->CornerOfEdge[j][1];
    }
    for (j=0; j<MGIO_MAX_EDGES_OF_ELEM; j++)
    {
      intList[s++] = pge->CornerOfSide[j][0];
      intList[s++] = pge->CornerOfSide[j][1];
      intList[s++] = pge->CornerOfSide[j][2];
      intList[s++] = pge->CornerOfSide[j][3];
    }
    if ((*Write_mint)(s,intList)) return (1);
    pge++;
  }

  return (0);
}

/****************************************************************************/
/*
   Read_RR_General - reads general information about refinement rules

   SYNOPSIS:
   int	Read_RR_General (MGIO_RR_GENERAL *mgio_rr_general);

   PARAMETERS:
   .  mgio_rr_general - general information about refinement rules

   DESCRIPTION:
   function reads general information about refinement rules

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Read_RR_General (MGIO_RR_GENERAL *mgio_rr_general)
{
  if ((*Read_mint)(1,intList)) return (1);
  mgio_rr_general->nRules = intList[0];

  return (0);
}

/****************************************************************************/
/*
   Write_GE_General - writes general information about refinement rules

   SYNOPSIS:
   int	Write_RR_General (MGIO_RR_GENERAL *mgio_rr_general);

   PARAMETERS:
   .  mgio_rr_general - general information about refinement rules

   DESCRIPTION:
   function writes general information about refinement rules

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Write_RR_General (MGIO_RR_GENERAL *mgio_rr_general)
{
  intList[0] = mgio_rr_general->nRules;
  if ((*Write_mint)(1,intList)) return (1);

  return (0);
}

/****************************************************************************/
/*
   Read_RR_Rules - reads refinement rules

   SYNOPSIS:
   int	Read_RR_Rules (int n, MGIO_RR_RULE *rr_rules);

   PARAMETERS:
   .  n - nb of refinement rules to read
   .  rr_rules - first refinement rule

   DESCRIPTION:
   function reads refinement rules

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Read_RR_Rules (int n, MGIO_RR_RULE *rr_rules)
{
  int i,j,k,m,s;
  MGIO_RR_RULE *prr;

  m = 2+MGIO_MAX_NEW_CORNERS+2*MGIO_MAX_NEW_CORNERS+MGIO_MAX_SONS_OF_ELEM*(1+MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_SIDES_OF_ELEM+1);
  prr = rr_rules;
  for (i=0; i<n; i++)
  {
    if ((*Read_mint)(m,intList)) return (1);
    s=0;
    prr->class = intList[s++];
    prr->nsons = intList[s++];
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
      prr->pattern[j] = intList[s++];
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
    {
      prr->sonandnode[j][0] = intList[s++];
      prr->sonandnode[j][1] = intList[s++];
    }
    for (j=0; j<MGIO_MAX_SONS_OF_ELEM; j++)
    {
      prr->sons[j].tag = intList[s++];
      for (k=0; k<MGIO_MAX_CORNERS_OF_ELEM; k++)
        prr->sons[j].corners[k] = intList[s++];
      for (k=0; k<MGIO_MAX_SIDES_OF_ELEM; k++)
        prr->sons[j].nb[k] = intList[s++];
      prr->sons[j].path = intList[s++];
    }
    prr++;
  }

  return (0);
}

/****************************************************************************/
/*
   Write_RR_Rules - writes refinement rules

   SYNOPSIS:
   int	Write_RR_Rules (int n, MGIO_RR_RULE *rr_rules);

   PARAMETERS:
   .  n - nb of refinement rules to write
   .  rr_rules - first refinement rule

   DESCRIPTION:
   function writes refinement rules

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Write_RR_Rules (int n, MGIO_RR_RULE *rr_rules)
{
  int i,j,k,m,s;
  MGIO_RR_RULE *prr;

  m = 2+MGIO_MAX_NEW_CORNERS+2*MGIO_MAX_NEW_CORNERS+MGIO_MAX_SONS_OF_ELEM*(1+MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_SIDES_OF_ELEM+1);
  prr = rr_rules;
  for (i=0; i<n; i++)
  {
    s=0;
    intList[s++] = prr->class;
    intList[s++] = prr->nsons;
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
      intList[s++] = prr->pattern[j];
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
    {
      intList[s++] = prr->sonandnode[j][0];
      intList[s++] = prr->sonandnode[j][1];
    }
    for (j=0; j<MGIO_MAX_SONS_OF_ELEM; j++)
    {
      intList[s++] = prr->sons[j].tag;
      for (k=0; k<MGIO_MAX_CORNERS_OF_ELEM; k++)
        intList[s++] = prr->sons[j].corners[k];
      for (k=0; k<MGIO_MAX_SIDES_OF_ELEM; k++)
        intList[s++] = prr->sons[j].nb[k];
      intList[s++] = prr->sons[j].path;
    }
    if ((*Write_mint)(s,intList)) return (1);
    prr++;
  }

  return (0);
}

/****************************************************************************/
/*
   Read_CG_General - reads general information about coarse grid

   SYNOPSIS:
   int Read_CG_General (MGIO_CG_GENERAL *cg_general);

   PARAMETERS:
   .  cg_general - general information about coarse grid

   DESCRIPTION:
   function reads general information about coarse grid

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Read_CG_General (MGIO_CG_GENERAL *cg_general)
{
  if ((*Read_mint)(6,intList)) return (1);
  cg_general->nPoint                      = intList[0];
  cg_general->nBndPoint           = intList[1];
  cg_general->nInnerPoint         = intList[2];
  cg_general->nElement            = intList[3];
  cg_general->nBndElement         = intList[4];
  cg_general->nInnerElement       = intList[5];

  return (0);
}

/****************************************************************************/
/*
   Write_CG_General - writes general information about coarse grid

   SYNOPSIS:
   int Write_CG_General (MGIO_CG_GENERAL *cg_general);

   PARAMETERS:
   .  cg_general - general information about coarse grid

   DESCRIPTION:
   function writes general information about coarse grid

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Write_CG_General (MGIO_CG_GENERAL *cg_general)
{
  intList[0] = cg_general->nPoint;
  intList[1] = cg_general->nBndPoint;
  intList[2] = cg_general->nInnerPoint;
  intList[3] = cg_general->nElement;
  intList[4] = cg_general->nBndElement;
  intList[5] = cg_general->nInnerElement;
  if ((*Write_mint)(6,intList)) return (1);

  return (0);
}

/****************************************************************************/
/*
   Read_CG_Points - reads coarse grid points

   SYNOPSIS:
   int Read_CG_Points (MGIO_CG_POINT *cg_point);

   PARAMETERS:
   .  n - nb of coarse grid points to read
   .  cg_point - first coarse grid point

   DESCRIPTION:
   function reads coarse grid points

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Read_CG_Points (int n, MGIO_CG_POINT *cg_point)
{
  int i,j,m,s;
  MGIO_CG_POINT *pp;

  m = MGIO_DIM*n;
  pp = cg_point;
  if ((*Read_mfloat)(m,floatList)) return (1);
  s=0;
  for(i=0; i<n; i++)
  {
    for(j=0; j<MGIO_DIM; j++)
      pp[i].position[j] = floatList[s++];
    pp++;
  }

  return (0);
}

/****************************************************************************/
/*
   Write_CG_Points - writes coarse grid points

   SYNOPSIS:
   int Write_CG_Points (MGIO_CG_POINT *cg_point);

   PARAMETERS:
   .  n - nb of coarse grid points to write
   .  cg_point - first coarse grid point

   DESCRIPTION:
   function writes coarse grid points

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Write_CG_Points (int n, MGIO_CG_POINT *cg_point)
{
  int i,j,s;
  MGIO_CG_POINT *pp;

  pp = cg_point;
  s=0;
  for(i=0; i<n; i++)
  {
    for(j=0; j<MGIO_DIM; j++)
      floatList[s++] = pp[i].position[j];
    pp++;
  }
  if ((*Write_mfloat)(s,floatList)) return (1);

  return (0);
}

/****************************************************************************/
/*
   Read_CG_Elements - reads coarse grid elements

   SYNOPSIS:
   int Read_CG_Elements (int n, MGIO_CG_ELEMENT *cg_element);

   PARAMETERS:
   .  n - nb of coarse grid elements to read
   .  cg_element - first coarse grid element

   DESCRIPTION:
   function reads coarse grid elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Read_CG_Elements (int n, MGIO_CG_ELEMENT *cg_element)
{
  int i,j,k,m,s;
  MGIO_CG_ELEMENT *pe;

  m = 1+MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_SIDES_OF_ELEM+2+MGIO_MAX_NEW_CORNERS;
  pe = cg_element;
  for (i=0; i<n; i++)
  {
    if ((*Read_mint)(m,intList)) return (1);
    s=0;
    pe->ge = intList[s++];
    for (j=0; j<MGIO_MAX_CORNERS_OF_ELEM; j++)
      pe->cornerid[j] = intList[s++];
    for (j=0; j<MGIO_MAX_SIDES_OF_ELEM; j++)
      pe->nbid[j] = intList[s++];
    pe->refrule = intList[s++];
    pe->nmoved = intList[s++];
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
      pe->moved[j] = intList[s++];
    m = MGIO_DIM*pe->nmoved;
    if ((*Read_mfloat)(m,floatList)) return (1);
    s=0;
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
    {
      if (!pe->moved[j]) continue;
      for (k=0; k<MGIO_DIM; k++)
        pe->newposition[j].position[k] = floatList[s++];
    }
    pe++;
  }

  return (0);
}

/****************************************************************************/
/*
   Write_CG_Elements - writes coarse grid elements

   SYNOPSIS:
   int Write_CG_Elements (int n, MGIO_CG_ELEMENT *cg_element);

   PARAMETERS:
   .  n - nb of coarse grid elements to write
   .  rr_rules - first coarse grid element

   DESCRIPTION:
   function writes coarse grid elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Write_CG_Elements (int n, MGIO_CG_ELEMENT *cg_element)
{
  int i,j,k,m,s;
  MGIO_CG_ELEMENT *pe;

  m = 1+MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_SIDES_OF_ELEM+2+MGIO_MAX_NEW_CORNERS;
  pe = cg_element;
  for (i=0; i<n; i++)
  {
    s=0;
    intList[s++] = pe->ge;
    for (j=0; j<MGIO_MAX_CORNERS_OF_ELEM; j++)
      intList[s++] = pe->cornerid[j];
    for (j=0; j<MGIO_MAX_SIDES_OF_ELEM; j++)
      intList[s++] = pe->nbid[j];
    intList[s++] = pe->refrule;
    intList[s++] = pe->nmoved;
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
      intList[s++] = pe->moved[j];
    if ((*Write_mint)(s,intList)) return (1);
    s=0;
    for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
    {
      if (!pe->moved[j]) continue;
      for (k=0; k<MGIO_DIM; k++)
        floatList[s++] = pe->newposition[j].position[k];
    }
    if ((*Write_mfloat)(s,floatList)) return (1);
    pe++;
  }

  return (0);
}

/****************************************************************************/
/*
   Read_HE_Refinement - reads hierarchical elements

   SYNOPSIS:
   int Read_HE_Refinement (MGIO_HE_ELEMENT *he_element);

   PARAMETERS:
   .  n - nb of hierarchical elements to read
   .  he_element - hierarchical element

   DESCRIPTION:
   function reads hierarchical elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Read_HE_Refinement (int n, MGIO_HE_ELEMENT *he_element)
{
  return (1);
}

/****************************************************************************/
/*
   Write_HE_Refinement - writes hierarchical elements

   SYNOPSIS:
   int Write_HE_Refinement (int n, MGIO_HE_ELEMENT *he_element);

   PARAMETERS:
   .  n - nb of hierarchical elements to write
   .  he_element - hierarchical element

   DESCRIPTION:
   function writes hierarchical elements

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Write_HE_Refinement (int n, MGIO_HE_ELEMENT *he_element)
{
  return (1);
}

/****************************************************************************/
/*
   Read_BD_General - reads general information about boundary description

   SYNOPSIS:
   int Read_BD_General (MGIO_BD_GENERAL *bd_general);

   PARAMETERS:
   .  bd_general - information about boundary description

   DESCRIPTION:
   function reads general information about boundary description

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Read_BD_General (MGIO_BD_GENERAL *bd_general)
{
  int i,s;

  if ((*Read_mint)(2*MGIO_MAXLEVEL,intList)) return (1);
  s=0;
  for (i=0; i<MGIO_MAXLEVEL; i++)
  {
    bd_general->nBndP[i] = intList[s++];
    bd_general->nBndS[i] = intList[s++];
  }

  return (0);
}

/****************************************************************************/
/*
   Write_BD_General - writes general information about boundary description

   SYNOPSIS:
   int Write_BD_General (MGIO_BD_GENERAL *bd_general);

   PARAMETERS:
   .  bd_general - information about boundary description

   DESCRIPTION:
   function writes general information about boundary description

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int Write_BD_General (MGIO_BD_GENERAL *bd_general)
{
  int i,s;

  s=0;
  for (i=0; i<MGIO_MAXLEVEL; i++)
  {
    intList[s++] = bd_general->nBndP[i];
    intList[s++] = bd_general->nBndS[i];
  }
  if ((*Write_mint)(s,intList)) return (1);

  return (0);
}

/****************************************************************************/
/*
   Read_PBndDesc - reads BNDPs

   SYNOPSIS:
   int Read_PBndDesc (MGIO_HEAP *theHeap, int n, BNDP **BndPList);

   PARAMETERS:
   .  theHeap - heap
   .  n - nb of BndP to read
   .  BndPList - list of ptrs to BndP

   DESCRIPTION:
   function reads BNDPs

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

#ifdef __MGIO_USE_IN_UG__

int Read_PBndDesc (MGIO_HEAP *theHeap, int n, BNDP **BndPList)
{
  int i;

  if (theHeap==NULL) return (1);
  for (i=0; i<n; i++)
  {
    BndPList[i] = BNDP_LoadBndP (theHeap,stream);
    if (BndPList[i]==NULL) return (1);
  }
  return (0);
}

#else

int Read_PBndDesc (MGIO_HEAP *theHeap, int n, BNDP **BndPList)
{
  return (1);
}

#endif

/****************************************************************************/
/*
   Write_PBndDesc - write BNDPs

   SYNOPSIS:
   int Write_PBndDesc (int n, BNDP **BndPList);

   PARAMETERS:
   .  n - nb of BndP to write
   .  BndPList - list of ptrs to BndP

   DESCRIPTION:
   function writes BNDPs

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

#ifdef __MGIO_USE_IN_UG__

int Write_PBndDesc (int n, BNDP **BndPList)
{
  int i;

  for (i=0; i<n; i++)
    if (BNDP_SaveBndP (BndPList[i],stream)) return (1);
  return (0);
}

#else

int Write_PBndDesc (int n, BNDP **BndPList)
{
  return (1);
}

#endif

/****************************************************************************/
/*
   Read_SBndDesc - reads BNDSs

   SYNOPSIS:
   int Read_SBndDesc (MGIO_HEAP *theHeap, int n, BNDS **BndSList);

   PARAMETERS:
   .  theHeap - heap
   .  n - nb of BndP to read
   .  BndSList - list of ptrs to BndS

   DESCRIPTION:
   function reads BNDSs

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

#ifdef __MGIO_USE_IN_UG__

int Read_SBndDesc (MGIO_HEAP *theHeap, int n, BNDP **BndSList)
{
  int i;

  if (theHeap==NULL) return (1);
  for (i=0; i<n; i++)
  {
    BndSList[i] = BNDS_LoadBndS (theHeap,stream);
    if (BndSList[i]==NULL) return (1);
  }
  return (0);
}

#else

int Read_SBndDesc (MGIO_HEAP *theHeap, int n, BNDS **BndSList)
{
  return (1);
}

#endif

/****************************************************************************/
/*
   Write_SBndDesc - write BNDSs

   SYNOPSIS:
   int Write_SBndDesc (int n, BNDS **BndSList);

   PARAMETERS:
   .  n - nb of BndS to write
   .  BndPList - list of ptrs to BndS

   DESCRIPTION:
   function writes BNDSs

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

#ifdef __MGIO_USE_IN_UG__

int Write_SBndDesc (int n, BNDS **BndSList)
{
  int i;

  for (i=0; i<n; i++)
    if (BNDS_SaveBndS (BndSList[i],stream)) return (1);
  return (0);
}

#else

int Write_SBndDesc (int n, BNDS **BndSList)
{
  return (1);
}

#endif

/****************************************************************************/
/*
   CloseFile - close the file

   SYNOPSIS:
   int CloseFile ();

   PARAMETERS:

   DESCRIPTION:
   close the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int CloseFile ()
{
  if (fclose(stream)!=NULL) return (1);
  return (0);
}

/****************************************************************************/
/*
   MGIO_Init - init input/output for mg

   SYNOPSIS:
   int MGIO_Init (void);

   PARAMETERS:

   DESCRIPTION:
   init the i/o of mg

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int MGIO_Init (void)
{

#ifdef __MGIO_USE_IN_UG__

  /* path to grid-dirs */
  gridpathes_set = 0;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"gridpathes")==0)
    gridpathes_set = 1;

#endif

  return (0);
}
