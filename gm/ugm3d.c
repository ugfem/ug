// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugm2d.c														*/
/*																			*/
/* Purpose:   Functions only available in 3D								*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   13.12.94 begin, ug version 3.0								*/
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

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"

#include "devices.h"

#include "switch.h"
#include "gm.h"
#include "evm.h"
#include "simplex.h"
#include "ugm.h"
#include "ugm3d.h"


#ifdef __TWODIM__
#error this source file is for 3D ONLY
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

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   DisposeEdgesFromElement - Dispose all edges of element not needed otherwise

   SYNOPSIS:
   INT DisposeEdgesFromElement  (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - grid to remove from
   .  theElement - element whose edges are to be removed

   DESCRIPTION:
   This function disposes all edge of element not needed otherwise.
   It runs through all edges of the element and decreases the 'NO_OF_ELEM' number.
   If this number is zero then the edge is removed.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if edge with zeroin elements found (before decrement).
   D*/
/****************************************************************************/

INT DisposeEdgesFromElement  (GRID *theGrid, ELEMENT *theElement)
{
  INT i,j,ret;
  EDGE *theEdge;

  ret=0;
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    for (j=i+1; j<CORNERS_OF_ELEM(theElement); j++)
    {
      theEdge = GetEdge(CORNER(theElement,i),CORNER(theElement,j));
      if (theEdge!=NULL)
      {
        if (NO_OF_ELEM(theEdge)<1)
          ret=1;
        if (NO_OF_ELEM(theEdge)==1)
          DisposeEdge(theGrid,theEdge);
        else
          DEC_NO_OF_ELEM(theEdge);
      }
    }

  return(ret);
}

/*******************************************************************/
/*
   MoveInnerNode - Assign new position to inner node

   SYNOPSIS:
   INT MoveInnerNode (MULTIGRID *theMG, NODE *theNode, COORD *newPos);

   PARAMETERS:
   .  theMG - pointer to 'MULTIGRID' structure
   .  theNode - pointer to 'NODE' to be moved
   .  newPos - new position

   DESCRIPTION:
   This function

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
 */
/*******************************************************************/

INT MoveInnerNode (MULTIGRID *theMG, NODE *theNode, COORD *newPos)
{
  UserWrite("MoveInnerNode not implemented in 3D version jet\n");
  return (0);
}

/*******************************************************************/
/*
   MoveBoundaryNode - Assign new position to boundary node

   SYNOPSIS:
   INT MoveBoundaryNode (MULTIGRID *theMG, NODE *theNode, INT segid,
   COORD *newPos);

   PARAMETERS:
   .  theMG - pointer to 'MULTIGRID' structure
   .  theNode - pointer to 'NODE' to be moved
   .  segid - new segment id
   .  newPos - new parameter (lambda,mu)

   DESCRIPTION:
   This function not implemented in 3D!

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
 */
/*******************************************************************/

INT MoveBoundaryNode (MULTIGRID *theMG, NODE *theNode, INT segid, COORD *newPos)
{
  UserWrite("MoveBoundaryNode not implemented in 3D version jet\n");
  return (0);
}

/*******************************************************************/
/*
   SmoothMultiGrid -

   SYNOPSIS:
   INT  SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag);

   PARAMETERS:
   .  theMG -
   .  niter -
   .  bdryFlag -

   DESCRIPTION:
   This function

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
 */
/*******************************************************************/

INT  SmoothMultiGrid (MULTIGRID *theMG, INT niter, INT bdryFlag)
{
  UserWrite("SmoothMultiGrid not implemented in 3D version jet\n");
  return (0);
}

/****************************************************************************/
/*
   SetParityOfElement - Set parity of element

   SYNOPSIS:
   INT SetParityOfElement (ELEMENT *theElement);

   PARAMETERS:
   .  theElement -	pointer to 'ELEMENT'

   DESCRIPTION:
   This function sets parity of element.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
 */
/****************************************************************************/

INT SetParityOfElement (ELEMENT *theElement)
{
  INT j;
  ELEMENT *e;
  ELEMENTSIDE *s;
  NODE *n;
  COORD sp,*x[4];
  COORD_VECTOR a,b,c,d;
        #ifdef __SIDEDATA__
  VECTOR *v;
        #endif

  if (TAG(theElement)!=TETRAHEDRON) return (1);
  if (EFATHER(theElement)!=NULL) return (1);
  if (NSONS(theElement)!=0) return (1);

  /* check orientation */
  for (j=0; j<4; j++) x[j] = CVECT(MYVERTEX(CORNER(theElement,j)));
  V3_SUBTRACT(x[CornerOfSide[0][1]],x[CornerOfSide[0][0]],a)
  V3_SUBTRACT(x[CornerOfSide[0][2]],x[CornerOfSide[0][0]],b)
  V3_SUBTRACT(x[OppositeCorner[0]],x[CornerOfSide[0][0]],c)
  V3_VECTOR_PRODUCT(a,b,d)
  V3_SCALAR_PRODUCT(c,d,sp)
  if (sp<0.0) return(0);

  /* invert orientaion */
  n = CORNER(theElement,0);
  SET_CORNER(theElement,0,CORNER(theElement,1));
  SET_CORNER(theElement,1,CORNER(theElement,2));
  SET_CORNER(theElement,2,CORNER(theElement,3));
  SET_CORNER(theElement,3,n);

  e = NBELEM(theElement,0);
  SET_NBELEM(theElement,0,NBELEM(theElement,1));
  SET_NBELEM(theElement,1,NBELEM(theElement,2));
  SET_NBELEM(theElement,2,NBELEM(theElement,3));
  SET_NBELEM(theElement,3,e);

        #ifdef __SIDEDATA__
  v = SVECTOR(theElement,0);
  SET_SVECTOR(theElement,0,SVECTOR(theElement,1));
  SET_SVECTOR(theElement,1,SVECTOR(theElement,2));
  SET_SVECTOR(theElement,2,SVECTOR(theElement,3));
  SET_SVECTOR(theElement,3,v);
        #endif

  if (OBJT(theElement)==IEOBJ) return(0);
  s = SIDE(theElement,0);
  SET_SIDE(theElement,0,SIDE(theElement,1));
  SET_SIDE(theElement,1,SIDE(theElement,2));
  SET_SIDE(theElement,2,SIDE(theElement,3));
  SET_SIDE(theElement,3,s);

  return(0);
}

/****************************************************************************/
/*
   SetParityInGrid - Set parity of elements in grid

   SYNOPSIS:
   INT SetParityInGrid (GRID* theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function sets parity of elements in grid.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
 */
/****************************************************************************/

INT SetParityInGrid (GRID* theGrid)
{
  ELEMENT *theElement;

  if (UPGRID(theGrid)!=NULL) return (1);
  if (DOWNGRID(theGrid)!=NULL) return (1);

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    if (SetParityOfElement(theElement)) return (1);

  return (0);
}

/****************************************************************************/
/*
   CheckParityOfElements - Check parity of elements in multigrid

   SYNOPSIS:
   INT CheckParityOfElements (MULTIGRID* theMG);

   PARAMETERS:
   .  theMG - pointer to multigrid

   DESCRIPTION:
   This function checks parity of elements in multigrid.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
 */
/****************************************************************************/

INT CheckParityOfElements (MULTIGRID* theMG)
{
  INT i,j;
  ELEMENT *theElement;
  COORD sp,*x[4];
  COORD_VECTOR a,b,c,d;
  char buffer[128];

  for (i=0; i<=TOPLEVEL(theMG); i++)
    for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,i)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (TAG(theElement)!=TETRAHEDRON) return (1);

      /* check orientation */
      for (j=0; j<4; j++) x[j] = CVECT(MYVERTEX(CORNER(theElement,j)));
      V3_SUBTRACT(x[CornerOfSide[0][1]],x[CornerOfSide[0][0]],a)
      V3_SUBTRACT(x[CornerOfSide[0][2]],x[CornerOfSide[0][0]],b)
      V3_SUBTRACT(x[OppositeCorner[0]],x[CornerOfSide[0][0]],c)
      V3_VECTOR_PRODUCT(a,b,d)
      V3_SCALAR_PRODUCT(c,d,sp)
      if (sp<0.0) continue;

      sprintf(buffer,"parity wrong: IS = %d\n",(int)ID(theElement));
      UserWrite(buffer);
    }

  return(0);
}
