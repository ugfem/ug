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
#include "defaults.h"

/* gm module */
#include "switch.h"
#include "gm.h"
#include "enrol.h"
#include "algebra.h"
#include "cw.h"
#include "ugm.h"
#include "ugio.h"
#include "elements.h"
#ifdef __TWODIM__
        #include "ugrefine2d.h"
#endif
#ifdef __THREEDIM__
        #include "ugrefine3d.h"
#endif

/* own header */
#include "initgm.h"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* data for CVS */
static char rcsid[] = "$Header$";

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

  /* elements.c */
  if ((err=InitElementTypes())!=GM_OK)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

        #ifdef __TWODIM__
  /* ugrefine2d.c */
  if ((err=InitRefine2d())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }
        #endif

        #ifdef __THREEDIM__
  /* ugrefine3d.c */
  if ((err=InitRefine3d())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }
        #endif

  return (0);
}
