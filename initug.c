// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initug.c														*/
/*																			*/
/* Purpose:   call the init routines of the ug modules						*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   27.02.95 begin, ug version 3.0								*/
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

/* ANSI-C includes */
#include <stdio.h>

/* low module */
#include "initlow.h"
#include "misc.h"
#include "general.h"
#include "defaults.h"

/* parallelization module */
#ifdef ModelP
#include "parallel.h"
#endif

/* devices module */
#include "devices.h"

/* domain module */
#include "domain.h"

/* grid manager module */
#include "initgm.h"
#include "switch.h"

/* numerics module */
#include "initnumerics.h"

/* graph module */
#include "initgraph.h"

/* user interface module */
#include "initui.h"

/* own header */
#include "initug.h"

/* TODO: delete this */
#include "debug.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define UGDEBUGRFILE            "debugfile"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*D
   InitUg - Call of the initfunctions for all the ug modules

   SYNOPSIS:
   INT InitUg (int *argcp, char ***argvp);

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

INT InitUg (int *argcp, char ***argvp)
{
  INT err;
        #if (defined Debug && defined __MWCW__)
  char buffer[256];
  char debugfilename[NAMESIZE];
        #endif

  /* init the low module */
  if ((err=InitLow())!=0)
  {
    printf("ERROR in InitUg while InitLow (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init parallelization module */
        #ifdef ModelP
  PRINTDEBUG(init,1,("%d:     InitParallel()...\n",me))
  if ((err=InitParallel(argcp, argvp))!=0)
  {
    printf("ERROR in InitUg while InitParallel (line %d): called routine line %d\n",(int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }
    #endif

  /* init the devices module */
  if ((err=InitDevices(*argcp,*argvp))!=0)
  {
    printf("ERROR in InitUg while InitDevices (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

        #if (defined Debug && defined __MWCW__)
  if ((GetDefaultValue(DEFAULTSFILENAME,UGDEBUGRFILE,buffer)==0)
      && (sscanf(buffer," %s ",&debugfilename)==1))
  {
    if (SetPrintDebugToFile(debugfilename)!=0)
    {
      printf("ERROR while opening debug file '%s'\n",debugfilename);
      printf ("aborting ug\n");
      return (1);
    }
    UserWriteF("debug info is captured to file '%s'\n",debugfilename);
  }
  else
  {
    SetPrintDebugProc(UserWriteF);
    UserWriteF("debug info is printed to ug's shell window\n");
  }
        #endif

  /* init the domain module */
  if ((err=InitDom())!=0)
  {
    printf("ERROR in InitDom while InitDom (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the gm module */
  if ((err=InitGm())!=0)
  {
    printf("ERROR in InitUg while InitGm (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the numerics module */
  if ((err=InitNumerics())!=0)
  {
    printf("ERROR in InitUg while InitNumerics (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the graph module */
  if ((err=InitGraph())!=0)
  {
    printf("ERROR in InitUg while InitGraph (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the ui module */
  if ((err=InitUi())!=0)
  {
    printf("ERROR in InitUg while InitUi (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  return (0);
}
