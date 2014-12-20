// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      tstep.c                                                       */
/*                                                                          */
/* Purpose:   time-step schemes                                             */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   Jan  31, 1998  begin                                          */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <config.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include "ugdevices.h"
#include "general.h"
#include "debug.h"
#include "ugstruct.h"
#include "misc.h"
#include "gm.h"
#include "evm.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "np.h"

#include "nls.h"
#include "ls.h"
#include "nls.h"
#include "assemble.h"
#include "transfer.h"

#include "reinit.h"

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

  REP_ERR_FILE

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* generic routines                                                         */
/*                                                                          */
/****************************************************************************/

INT NS_DIM_PREFIX REINIT_Display (NP_BASE *base)
{
  NP_REINIT *reinit;
  INT i;

  reinit=(NP_REINIT*)base;
  UserWrite("\nreinit configuration:\n");
  for (i=0; i<reinit->n; i++)
    UserWriteF(DISPLAY_NP_FORMAT_SF,reinit->name[i],(float)reinit->parameter[i]);

  return (0);
}

/****************************************************************************/
/*
   InitReinit  - Init this file

   SYNOPSIS:
   INT InitReinit ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitReinit ()
{
  return (0);
}
