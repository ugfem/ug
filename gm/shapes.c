// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  shapes.c														*/
/*																			*/
/* Purpose:   shape functions for triangles and quadrilaterals				*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					*/
/*																			*/
/* History:   08.04.92 begin, ug version 2.0								*/
/*			  20.11.94 moved shapes.c from ug to cd folder					*/
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

#include <math.h>
#include <assert.h>

#include "compiler.h"
#include "misc.h"
#include "gm.h"
#include "evm.h"
#include "shapes.h"

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
   GN - General Shape function for nodes

   SYNOPSIS:
   DOUBLE GN (INT n, INT i, COORD local);

   PARAMETERS:
   .  n - number of corners of the element
   .  i - corner number (corner number [0..n-1])
   .  local - local COORDinates

   DESCRIPTION:
   This function finds the linear/bilinear shape functions Ni(local) to approximate the
   solution u in the integration point ip for triangles/quadrilaterals/tetrahedra

   .n   uip(local) = SUM Ni(local)*ui

   where the sum runs over all nodes of the element to which the considered
   ip belongs. The shape function is defined as

   - for all elements who do not have the node i as a corner
   .n   Ni = 0

   - for the elements
   .n   Ni(node i) = 1
   .n   Ni(node k) = 0, if k is not equal i.

   RETURN VALUE:
   DOUBLE
   .n          value
   D*/
/****************************************************************************/

DOUBLE GN (INT n, INT i, COORD *local)
{
#ifdef __TWODIM__
  switch (n)
  {
  case 3 :
    switch (i)
    {
    case 0 : return((DOUBLE)(1-local[0]-local[1]));
    case 1 : return((DOUBLE)local[0]);
    case 2 : return((DOUBLE)local[1]);
    }
  case 4 :
    switch (i)
    {
    case 0 : return((DOUBLE)(0.25*(1-local[0])*(1-local[1])));
    case 1 : return((DOUBLE)(0.25*(1+local[0])*(1-local[1])));
    case 2 : return((DOUBLE)(0.25*(1+local[0])*(1+local[1])));
    case 3 : return((DOUBLE)(0.25*(1-local[0])*(1+local[1])));
    }
  default :
    return (-1.0);
  }
#endif

#ifdef __THREEDIM__
  switch (n)
  {
  case 4 :
    switch (i)
    {
    case 0 : return((DOUBLE)(1.0-local[0]-local[1]-local[2]));
    case 1 : return((DOUBLE)local[0]);
    case 2 : return((DOUBLE)local[1]);
    case 3 : return((DOUBLE)local[2]);
    }
  default :
    return (-1.0);
  }
#endif
}
