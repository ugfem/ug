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
#endif

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

INT InitParallel (int *argcp, char ***argvp)
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

        #ifdef __DLB__
  /* init dlb module */
  PRINTDEBUG(init,1,("%d:     InitParallel()...\n",me))
  if ((err=DLB_Init())!=0)
  {
    printf("ERROR in InitParallel while DLB_Init (line %d): called routine line %d\n",(int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }
        #endif

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

INT ExitParallel (void)
{
  INT err;

  /* exit parallelization module */
  /* the following code (ExitDDD) once seemed to crash
     with MPI-PPIF. today it seems to run without problems.
     therefore, we switch it on again, if there are any problems with MPI
     and exiting the program, it may come from here. KB 970527 */

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
