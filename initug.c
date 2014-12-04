// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file initug.c
 * \ingroup ug
 */

/** \addtogroup ug
 *
 * @{
 */


/****************************************************************************/
/*                                                                          */
/* File:      initug.c                                                      */
/*                                                                          */
/* Purpose:   call the init routines of the ug modules                      */
/*                                                                          */
/* Author:    Henrik Rentz-Reichert                                         */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   27.02.95 begin, ug version 3.0                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*    system include files                                                  */
/*    application include files                                             */
/*                                                                          */
/****************************************************************************/

/* ANSI-C includes */
#include <config.h>
#include <stdio.h>
#include <string.h>

/* low module */
#include "initlow.h"
#include "misc.h"
#include "general.h"
#include "defaults.h"
#include "ugstruct.h"

/* parallelization module */
#ifdef ModelP
#include "initparallel.h"
#include "parallel.h"
#include "ppif.h"
USING_PPIF_NAMESPACE
#endif

/* devices module */
#include "dev/ugdevices.h"

/* domain module */
#include "dom/domain.h"

/* grid manager module */
#include "gm/initgm.h"

/* numerics module */
#include "np/initnumerics.h"

/* graphics module */
#include "graphics/graphics.h"

/* user interface module */
#include "ui/initui.h"

/* own header */
#include "initug.h"

/** \todo delete this */
#include "low/debug.h"


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

#define UGDEBUGRFILE            "debugfile"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char
RCS_ID
  ("$Header$",
  UG_RCS_STRING);

/****************************************************************************/
/** \brief Call the init functions for all the ug modules
 *
 * @param argcp - pointer to argument counter
 * @param argvp - pointer to argument vector
 *
 *   This function initializes.
 *
 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error occured. </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitUg (int *argcp, char ***argvp)
{
  INT err;
#ifdef Debug
  char buffer[256];
  char debugfilename[NAMESIZE];
#endif

#ifdef ModelP
  /* init ppif module */
  if ((err = InitPPIF (argcp, argvp)) != PPIF_SUCCESS)
  {
    printf ("ERROR in InitParallel while InitPPIF.\n");
    printf ("aborting ug\n");

    return (1);
  }
#endif

  /* init the low module */
  if ((err = InitLow ()) != 0)
  {
    printf
      ("ERROR in InitUg while InitLow (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init parallelization module */
#ifdef ModelP
  PRINTDEBUG (init, 1, ("%d:     InitParallel()...\n", me))
  if ((err = InitParallel (argcp, argvp)) != 0)
  {
    printf
      ("ERROR in InitUg while InitParallel (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }
#endif

  /* create struct for configuration parameters */
  if (MakeStruct (":conf"))
    return (__LINE__);
  if (SetStringVar ("conf:arch", ARCHNAME))
    return (__LINE__);

  /* set variable for parallel modus */
#ifdef ModelP
  if (SetStringValue ("conf:parallel", 1.0))
    return (__LINE__);
  if (SetStringValue ("conf:procs", (DOUBLE) procs))
    return (__LINE__);
  if (SetStringValue ("conf:me", (DOUBLE) me))
    return (__LINE__);
#else
  if (SetStringValue ("conf:parallel", 0.0))
    return (__LINE__);
  if (SetStringValue ("conf:procs", 1.0))
    return (__LINE__);
  if (SetStringValue ("conf:me", 0.0))
    return (__LINE__);
#endif

  /* init the devices module */
  if ((err = InitDevices (argcp, *argvp)) != 0)
  {
    printf
      ("ERROR in InitUg while InitDevices (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

#ifdef Debug
  {
    int i;
    for (i = 1; i < *argcp; i++)
      if (strncmp ((*argvp)[i], "-dbgfile", 8) == 0)
        break;
    if ((i < *argcp)
        && (GetDefaultValue (DEFAULTSFILENAME, UGDEBUGRFILE, buffer) == 0)
        && (sscanf (buffer, " %s ", debugfilename) == 1))
    {
      if (SetPrintDebugToFile (debugfilename) != 0)
      {
        printf ("ERROR while opening debug file '%s'\n", debugfilename);
        printf ("aborting ug\n");
        return (1);
      }
      UserWriteF ("debug info is captured to file '%s'\n", debugfilename);
    }
    else
    {
      SetPrintDebugProc (printf);
      UserWriteF ("debug info is printed to stdout\n");
    }
  }
#endif

  /* init the domain module */
  if ((err = InitDom ()) != 0)
  {
    printf
      ("ERROR in InitDom while InitDom (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the gm module */
  if ((err = InitGm ()) != 0)
  {
    printf
      ("ERROR in InitUg while InitGm (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the numerics module */
  if ((err = InitNumerics ()) != 0)
  {
    printf
      ("ERROR in InitUg while InitNumerics (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the ui module */
  if ((err = InitUi (argcp[0], argvp[0])) != 0)
  {
    printf
      ("ERROR in InitUg while InitUi (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* init the graphics module */
  if ((err = InitGraphics ()) != 0)
  {
    printf
      ("ERROR in InitUg while InitGraphics (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  return (0);
}


/****************************************************************************/
/** \brief Call of the exitfunctions for all the ug modules
 *
 * This function exits ug. It is called at the end of the CommandLoop.
 * It calls all available exit functions in reverse order of the corresponding
 * calls in InitUg().
 *
 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error occured. </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX
ExitUg (void)
{
  INT err;


  /* exit graphics module */
  PRINTDEBUG (init, 1, ("%d:     ExitGraphics()...\n", me))
  if ((err = ExitGraphics ()) != 0)
  {
    printf
      ("ERROR in ExitUg while ExitGraphics (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* exit ui module */
  PRINTDEBUG (init, 1, ("%d:     ExitUi()...\n", me))
  if ((err = ExitUi ()) != 0)
  {
    printf
      ("ERROR in ExitUg while ExitUi (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* exit gm module */
  PRINTDEBUG (init, 1, ("%d:     ExitGm()...\n", me))
  if ((err = ExitGm ()) != 0)
  {
    printf
      ("ERROR in ExitUg while ExitGm (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  /* exit devices module */
  PRINTDEBUG (init, 1, ("%d:     ExitDevices()...\n", me))
  if ((err = ExitDevices ()) != 0)
  {
    printf
      ("ERROR in ExitUg while ExitDevices (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }


  /* exit parallelization module */
#ifdef ModelP

  /* the following code (ExitParallel) once seemed to crash
     with MPI-PPIF. today it seems to run without problems.
     therefore, we switch it on again, if there are any problems with MPI
     and exiting the program, it may come from here. KB 970527 */

  PRINTDEBUG (init, 1, ("%d:     ExitParallel()...\n", me))
  if ((err = ExitParallel ()) != 0)
  {
    printf
      ("ERROR in ExitUg while ExitParallel (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

#endif

  /* exit low module */
  PRINTDEBUG (init, 1, ("%d:     ExitLow()...\n", me))
  if ((err = ExitLow ()) != 0)
  {
    printf
      ("ERROR in ExitUg while ExitLow (line %d): called routine line %d\n",
      (int) HiWrd (err), (int) LoWrd (err));
    printf ("aborting ug\n");

    return (1);
  }

  return (0);
}

/** @} */
