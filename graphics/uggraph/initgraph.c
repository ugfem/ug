// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initgraph.c													*/
/*																			*/
/* Purpose:   call the init routines of the graph module					*/
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
#include "compiler.h"
#include "misc.h"
#include "general.h"

/* graph module */
#include "wpm.h"
#include "wop.h"
#include "plotproc.h"

/* own header */
#include "initgraph.h"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*
   InitUGGraph -  Call the inits for the graph module

   SYNOPSIS:
   INT InitUGGraph ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function calls the inits for the graph module.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if some error occured.
 */
/****************************************************************************/

INT InitUGGraph (void)
{
  INT err;

  /* init wpm.c */
  if ((err=InitWPM())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init wop.c */
  if ((err=InitWOP())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init plotproc.c */
  if ((err=InitPlotProc())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if (SetStringValue("Devices:nWindows",0.0))
    return(__LINE__);

  return (0);
}
