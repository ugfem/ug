// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  domain.c														*/
/*																			*/
/* Purpose:   general function bases on domain.h-IF-functions                           */
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
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

/* standard C library */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* low modules */
#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"

/* domain module */
#include "domain.h"

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
   Patch_GetPatchByID - get patch from id

   SYNOPSIS:
   PATCH *Patch_GetPatchByID (BVP *theBVP, INT id);

   PARAMETERS:
   .  theBVP - the BoundaryValueProblem
   .  thePatch - the patch

   DESCRIPTION:
   This function finds the patch from id by running the patch-list

   THE RETURN VALUE:
   PATCH *
   .n      pointer to PATCH
   .n      NULL if not found or error.
   D*/
/****************************************************************************/

PATCH *Patch_GetPatchByID (BVP *theBVP, INT id)
{
  PATCH *thePatch;
  PATCH_DESC thePatchDesc;

  for (thePatch=BVP_GetFirstPatch(theBVP); thePatch!=NULL; thePatch=BVP_GetNextPatch(theBVP,thePatch))
  {
    if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return (NULL);
    if (PATCH_ID(thePatchDesc)==id) break;
  }

  return (thePatch);
}

/****************************************************************************/
/*D
   Patch_GetPatchID - get id from patch from

   SYNOPSIS:
   INT Patch_GetPatchID (PATCH *thePatch);

   PARAMETERS:
   .  thePatch - the patch

   DESCRIPTION:
   This function gets the id of the patch

   THE RETURN VALUE:
   INT
   .n      >= if ok
   .n      -1 if error.
   D*/
/****************************************************************************/

INT Patch_GetPatchID (PATCH *thePatch)
{
  PATCH_DESC thePatchDesc;

  if (Patch_GetPatchDesc(thePatch,&thePatchDesc)) return (-1);
  return(PATCH_ID(thePatchDesc));
}
