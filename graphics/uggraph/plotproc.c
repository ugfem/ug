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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"

#include "gm.h"
#include "ugenv.h"
#include "evm.h"
#include "devices.h"
#include "shapes.h"
#include "num.h"

#include "plotproc.h"

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

static INT nodecomp;

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

static INT PreprocessNodeValue (const char *name, MULTIGRID *theMG)
{
  VECDATA_DESC *theVD;

  theVD = GetVecDataDescByName(theMG,(char *)name);

  if (theVD == NULL) {
    PrintErrorMessage('E',"PreprocessNodeValue","cannot find symbol");
    return (1);
  }

  if (VD_NCMPS_IN_TYPE(theVD,NODEVECTOR)<1)
    return (1);

  nodecomp = VD_CMP_OF_TYPE(theVD,NODEVECTOR,0);

  return (0);
}

static DOUBLE NodeValue (const ELEMENT *theElement,
                         const DOUBLE **CornersCoord, DOUBLE *LocalCoord)
{
  INT i,n;
  DOUBLE phi;

  n = CORNERS_OF_ELEM(theElement);
  phi = 0.0;
  for (i=0; i<n; i++)
    phi += GN(n,i,LocalCoord)*VVALUE(NVECTOR(CORNER(theElement,i)),nodecomp);

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

static INT PreprocessNodeVector (const char *name, MULTIGRID *theMG)
{
  VECDATA_DESC *theVD;
  INT i;

  theVD = GetVecDataDescByName(theMG,(char *)name);

  if (theVD == NULL) {
    PrintErrorMessage('E',"PreprocessNodeVector","cannot find symbol");
    return (1);
  }

  if (VD_NCMPS_IN_TYPE(theVD,NODEVECTOR)<DIM)
    return (1);

  nodecomp = VD_CMP_OF_TYPE(theVD,NODEVECTOR,0);

  for (i=1; i<DIM; i++)
    if ((nodecomp+i) != VD_CMP_OF_TYPE(theVD,NODEVECTOR,i))
      return (1);

  return (0);
}

static void NodeVector (const ELEMENT *theElement, const DOUBLE **theCorners,
                        DOUBLE *LocalCoord, DOUBLE *values)
{
  VECTOR *v;
  INT i,j,n;
  DOUBLE s;

  n = CORNERS_OF_ELEM(theElement);
  V_DIM_CLEAR(values);

  for (i=0; i<n; i++)
  {
    v = NVECTOR(CORNER(theElement,i));
    s = GN(n,i,LocalCoord);
    for (j=0; j<DIM; j++)
      values[j] += s*VVALUE(v,nodecomp+j);
  }

  return;
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

INT InitPlotProc ()
{
  /* install general plot procs */
  if (CreateElementValueEvalProc("nvalue",PreprocessNodeValue,
                                 NodeValue)==NULL)
    return(1);
  if (CreateElementVectorEvalProc("nvector",PreprocessNodeVector,
                                  NodeVector,DIM)==NULL)
    return(1);

  return (0);
}
