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
/*			  email: ug@ica3.uni-stuttgart.de							*/
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
#include "misc.h"
#include "initlow.h"

/* devices module */
#include "devices.h"

/* domain module */
#include "domain.h"

/* grid manager module */
#include "initgm.h"
#include "switch.h"

/* grid generator module */
#include "ggmain.h"

/* numerics module */
#include "initnumerics.h"

/* graph module */
#include "initgraph.h"

/* user interface module */
#include "initui.h"

/* own header */
#include "initug.h"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*D
   InitUg - Call of the initfunctions for all the ug modules

   SYNOPSIS:
   INT InitUg (int argc, char **argv);

   PARAMETERS:
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function initializes.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

INT InitUg (int argc, char **argv)
{
  INT err;

  /* init the low module */
  if ((err=InitLow())!=0)
  {
    printf("ERROR in InitUg while InitLow (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the devices module */
  if ((err=InitDevices(argc,argv))!=0)
  {
    printf("ERROR in InitUg while InitDevices (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }

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

    #ifdef __TWODIM__
  /* init the gg module */
  if ((err=InitGG())!=0)
  {
    printf("ERROR in InitUg while InitGG (line %d): called routine line %d\n",
           (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ug\n");

    return (1);
  }
    #endif

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
