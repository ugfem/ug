// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initnp.c			                                                                                */
/*																			*/
/* Purpose:   call the init routines of the numerics module		                        */
/*																			*/
/* Author:	  Klaus Johannsen		                                                                                */
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

#include <stdio.h>

#include "np.h"
#include "general.h"
#include "basics.h"
#include "amgtransfer.h"
#include "transfer.h"
#include "assemble.h"
#include "iter.h"
#include "ls.h"
#include "nls.h"
#include "error.h"
#include "fvgeom.h"
#include "udm.h"
#include "formats.h"
#include "dio.h"
#include "newton.h"
#include "ts.h"
#include "bdf.h"
#include "ew.h"
#include "field.h"
#include "db.h"

#include "initnp.h"
#include "numproc.h"
#include "amg_ug.h"

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*D
   InitNumerics	- Call the inits for the grid manger module

   SYNOPSIS:
   INT InitNumerics ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function does initialization of the np subsystem.
   It is called in InitUG.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT InitNumerics ()
{
  INT err;

  /* init procs */
  if ((err=InitNumProcManager())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitTransfer())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitAMGTransfer())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitLinearSolver())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitNewtonSolver())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitAssemble())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitBDFSolver())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitTSolver())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitNonlinearSolver())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitIter())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitBasics())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitError())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitEW())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init finite volumes */
  if ((err=InitFiniteVolumeGeom())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init user data manager */
  if ((err=InitUserDataManager())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }
  if ((err=InitFormats())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init data io */
  if ((err=DIO_Init())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init stochastic field */
  if ((err=InitStochField())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init data base */
  if ((err=InitDb())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init amg solver */
  if ((err=InitAMGSolver())!=0) {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  return (0);
}
