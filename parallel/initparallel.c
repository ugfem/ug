// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initparallel.c												*/
/*																			*/
/* Purpose:   call the init routines of the parallel modules				*/
/*																			*/
/* Author:	  Stefan Lang                                                                           */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   980604 start                                                  */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* ANSI includes */
#include "config.h"
#include <stdio.h>

/* low module */
#include "initlow.h"
#include "misc.h"
#include "general.h"
#include "defaults.h"
#include "ugstruct.h"

/* parallelization module */
#ifdef ModelP
#include "parallel.h"
#include "ppif.h"
#include "initparallel.h"
#endif

#include "namespace.h"

/* UG namespaces: */
USING_UG_NAMESPACES

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*D
   InitParallel - Call of the initfunctions for all parallel modules

   SYNOPSIS:
   INT InitParallel (int *argcp, char ***argvp);

   PARAMETERS:
   .  argcp - pointer to argument counter
   .  argvp - pointer to argument vector

   DESCRIPTION:
   This function initializes.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

INT NS_DIM_PREFIX InitParallel (int *argcp, char ***argvp)
{
  INT err;

  /* init ddd module */
  PRINTDEBUG(init,1,("%d:     InitParallel()...\n",me))
  if ((err=InitDDD())!=0)
  {
    printf("ERROR in InitParallel while InitDDD (line %d): called routine line %d\n",(int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  return (0);
}


/****************************************************************************/
/*D
   ExitParallel - Call of the exitfunctions for all parallel modules

   SYNOPSIS:
   INT ExitParallel (void);

   PARAMETERS:

   DESCRIPTION:
   This function exits parallel modules.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

INT NS_DIM_PREFIX ExitParallel (void)
{
  INT err;

  /* exit parallelization module */
  PRINTDEBUG(init,1,("%d:     ExitParallel()...\n",me))
  if ((err=ExitDDD())!=0)
  {
    printf("ERROR in ExitParallel while ExitDDD (line %d): called routine line %d\n",(int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  return (0);
}

#endif
