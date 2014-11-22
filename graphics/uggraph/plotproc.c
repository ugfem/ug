// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  plotproc.c													*/
/*																			*/
/* Purpose:   generic plot functions for node vectors                           */
/*																			*/
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   November 23, 1996                                                                         */
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

#include <config.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"

#include "gm.h"
#include "elements.h"
#include "ugenv.h"
#include "evm.h"
#include "ugdevices.h"
#include "shapes.h"
#include "num.h"
#include "cw.h"

#include "plotproc.h"

#include "ppif_namespace.h"

USING_UG_NAMESPACES
USING_PPIF_NAMESPACE

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

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

static INT NodeValueComp;
static INT ElementValueComp;
static INT NodeVectorComp;
static INT ElementVectorComp;
static INT GradientFlag;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/**************************************************************************/
/*
   NodeValue - generic plot funtion

   SYNOPSIS:
   static DOUBLE NodeValue (const ELEMENT *theElement,
   const DOUBLE **CornersCoord, DOUBLE *LocalCoord);

   PARAMETERS:
   .  theElement - pointer to an element
   .  CornersCoord - corner coordinates
   .  LocalCoord - local coordinate

   DESCRIPTION:
   This function plots the values of a node vector, interpolated by the
   standard shape functions.

   RETURN VALUE:
   DOUBLE  value
 */
/*************************************************************************/

static INT PreProcessNodeValue (const char *name, MULTIGRID *theMG)
{
  VECDATA_DESC *theVD;

  theVD = GetVecDataDescByName(theMG,(char *)name);

  if (theVD == NULL) {
    PrintErrorMessage('E',"PreProcessNodeValue","cannot find symbol");
    return (1);
  }

  if (VD_ncmps_in_otype(theVD,NODEVEC)<1)
    return (1);

  NodeValueComp = VD_cmp_of_otype(theVD,NODEVEC,0);

  return (0);
}

static DOUBLE NodeValue (const ELEMENT *theElement, const DOUBLE **CornersCoord, DOUBLE *LocalCoord)
{
  INT i,n;
  DOUBLE phi;

  n = CORNERS_OF_ELEM(theElement);
  phi = 0.0;
  for (i=0; i<n; i++)
    phi += GN(n,i,LocalCoord)*VVALUE(NVECTOR(CORNER(theElement,i)),NodeValueComp);

  return(phi);
}

static INT PreProcessElementValue (const char *name, MULTIGRID *theMG)
{
  VECDATA_DESC *theVD;

  theVD = GetVecDataDescByName(theMG,(char *)name);

  if (theVD == NULL) {
    PrintErrorMessage('E',"PreProcessNodeValue","cannot find symbol");
    return (1);
  }

  if (VD_ncmps_in_otype(theVD,ELEMVEC)<1)
    return (1);

  ElementValueComp = VD_cmp_of_otype(theVD,ELEMVEC,0);

  return (0);
}

static DOUBLE ElementValue (const ELEMENT *theElement, const DOUBLE **CornersCoord, DOUBLE *LocalCoord)
{
  return(VVALUE(EVECTOR(theElement),ElementValueComp));
}

static DOUBLE LevelValue (const ELEMENT *theElement, const DOUBLE **CornersCoord, DOUBLE *LocalCoord)
{
  DOUBLE phi;

  phi = LEVEL(theElement);
  if (phi>0 && ECLASS(theElement)<=GREEN_CLASS) phi -= 1;
  return(phi);
}

/**************************************************************************/
/*
   NodeVector - plots the vector of node values

   SYNOPSIS:
   static void NodeVector (const ELEMENT *theElement,
   const DOUBLE **theCorners, DOUBLE *LocalCoord, DOUBLE *values);

   PARAMETERS:
   .  theElement - pointer to an element
   .  CornersCoord - corner coordinates
   .  LocalCoord - local coordinate
   .  values - plot vector

   DESCRIPTION:
   This function plots the vector of node values.

   RETURN VALUE:
   void
 */
/*************************************************************************/

static INT PreProcessNodeVector (const char *name, MULTIGRID *theMG)
{
  VECDATA_DESC *theVD;
  INT i;

  theVD = GetVecDataDescByName(theMG,(char *)name);
  if (theVD == NULL) {
    PrintErrorMessage('E',"PreProcessNodeVector","cannot find symbol");
    return (1);
  }
  NodeVectorComp = VD_cmp_of_otype(theVD,NODEVEC,0);
  if (VD_ncmps_in_otype(theVD,NODEVEC)<DIM)
    GradientFlag = 1;
  else {
    GradientFlag = 0;
    for (i=1; i<DIM; i++)
      if ((NodeVectorComp+i) != VD_cmp_of_otype(theVD,NODEVEC,i))
        return (1);
  }
  return (0);
}

static void NodeVector (const ELEMENT *theElement, const DOUBLE **theCorners,
                        DOUBLE *LocalCoord, DOUBLE *values)
{
  VECTOR *v;
  INT i,j,n;
  DOUBLE s;
  DOUBLE_VECTOR grad;

  n = CORNERS_OF_ELEM(theElement);
  V_DIM_CLEAR(values);

  if (GradientFlag)
    for (i=0; i<n; i++) {
      v = NVECTOR(CORNER(theElement,i));
      D_GN(n,i,LocalCoord,grad);
      s = VVALUE(v,NodeVectorComp);
      V_DIM_SCALE(s,grad);
      V_DIM_ADD(grad,values,values);
    }
  else
    for (i=0; i<n; i++) {
      v = NVECTOR(CORNER(theElement,i));
      s = GN(n,i,LocalCoord);
      for (j=0; j<DIM; j++)
        values[j] += s * VVALUE(v,NodeVectorComp+j);
    }

  return;
}

static INT PreProcessElementVector (const char *name, MULTIGRID *theMG)
{
  VECDATA_DESC *theVD;

  theVD = GetVecDataDescByName(theMG,(char *)name);

  if (theVD == NULL) {
    PrintErrorMessage('E',"PreProcessNodeValue","cannot find symbol");
    return (1);
  }

  if (VD_ncmps_in_otype(theVD,ELEMVEC)<DIM)
    return (1);

  ElementVectorComp = VD_cmp_of_otype(theVD,ELEMVEC,0);

  return (0);
}

static void ElementVector (const ELEMENT *theElement, const DOUBLE **theCorners,
                           DOUBLE *LocalCoord, DOUBLE *values)
{
  int i;

  for (i = 0; i < DIM; i++)
    values[i] = VVALUE(EVECTOR(theElement), ElementVectorComp+i);
}

/**************************************************************************/
/*
   RefMarks - plot funtion for refinement marks

   SYNOPSIS:
   static DOUBLE RefMarks (const ELEMENT *theElement,
   const DOUBLE **CornersCoord, DOUBLE *LocalCoord);

   PARAMETERS:
   .  theElement - pointer to an element
   .  CornersCoord - corner coordinates
   .  LocalCoord - local coordinate

   DESCRIPTION:
   This function plots the marks returns by 'GetRefienmentMark'.

   RETURN VALUE:
   DOUBLE  1.0  for rule == RED
           0.0  for rule == NO_REFINEMENT
                  -1.0  for rule == COARSE
 */
/*************************************************************************/

static INT PreProcessRefMarks (const char *name, MULTIGRID *theMG)
{
  if (CURRENTLEVEL(theMG) != TOPLEVEL(theMG)) {
    PrintErrorMessage('E',"PreProcessRefMarks",
                      "ref marks can be plotted on toplevel only");
    return (1);
  }
  return (0);
}

static DOUBLE RefMarks (const ELEMENT *theElement,
                        const DOUBLE **CornersCoord, DOUBLE *LocalCoord)
{
  INT rule,side;

  GetRefinementMark((ELEMENT *)theElement,&rule,(void *)&side);

  switch (rule) {
  case RED :           return( 1.0);
  case NO_REFINEMENT : return( 0.0);
  case COARSE :        return(-1.0);
  }
  return(0.0);
}

/**************************************************************************/
/*   ProcID - plot funtion for processor id                               */
/**************************************************************************/

static DOUBLE ProcID(const ELEMENT *e, const DOUBLE **cc, DOUBLE *lc)
{
  return (DOUBLE)me;
}

/**************************************************************************/
/*   SubDomID - plot funtion for subdomain id                             */
/**************************************************************************/

static DOUBLE SubDomID(const ELEMENT *e, const DOUBLE **cc, DOUBLE *lc)
{
  return (DOUBLE)SUBDOMAIN(e);
}

/****************************************************************************/
/*
   InitPlotProc	- Init this file

   SYNOPSIS:
   INT InitPlotProc ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitPlotProc ()
{
  /* install general plot procs */
  if (CreateElementValueEvalProc("nvalue",PreProcessNodeValue,NodeValue) == NULL) return(1);
  if (CreateElementValueEvalProc("evalue",PreProcessElementValue,ElementValue) == NULL) return(1);
  if (CreateElementValueEvalProc("level",NULL,LevelValue) == NULL) return(1);
  if (CreateElementVectorEvalProc("nvector",PreProcessNodeVector,NodeVector,DIM) == NULL) return(1);
  if (CreateElementVectorEvalProc("evector",PreProcessElementVector,ElementVector,DIM) == NULL) return(1);
  if (CreateElementValueEvalProc("refmarks",PreProcessRefMarks,RefMarks) == NULL) return(1);
  if (CreateElementValueEvalProc("procid",NULL,ProcID) == NULL) return(1);
  if (CreateElementValueEvalProc("subdomid",NULL,SubDomID) == NULL) return(1);
  return (0);
}
