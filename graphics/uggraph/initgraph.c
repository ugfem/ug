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
#include "compiler.h"
#include "misc.h"

/* graph module */
#include "wpm.h"
#include "wop.h"

/* own header */
#include "initgraph.h"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*D
   InitGraph -  Call the inits for the graph module

   SYNOPSIS:
   INT InitGraph ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function calls the inits for the graph module.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if some error occured.
   D*/
/****************************************************************************/

INT InitGraph ()
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

  return (0);
}
