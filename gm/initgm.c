// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initgm.c														*/
/*																			*/
/* Purpose:   call the init routines of the grid manager module                         */
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
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
#include "defaults.h"
#include "general.h"

/* gm module */
#include "switch.h"
#include "gm.h"
#include "enrol.h"
#include "algebra.h"
#include "cw.h"
#include "ugm.h"
#include "ugio.h"
#include "elements.h"
#include "refine.h"
#include "rm.h"

#ifdef __TWODIM__
/* grid generator module */
#include "ggmain.h"
#endif

/* own header */
#include "initgm.h"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*
   InitGm - Call the inits for the grid manger module

   SYNOPSIS:
   INT InitGm ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function calls the inits for the grid manger module.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if some error occured.
 */
/****************************************************************************/

INT InitGm ()
{
  INT err;

  /* cw.c */
  if ((err=InitPredefinedControlEntries())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* enrol.c */
  if ((err=InitEnrol())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* algebra.c */
  if ((err=InitAlgebra())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* ugm.c */
  if ((err=InitUGManager())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* ugio.c */
  if ((err=InitUgio())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init evalproc.c */
  if ((err=InitEvalProc())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* rm.c */
  if ((err=InitRuleManager())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
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

  return (0);
}
