// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initlow.c                                                                                                     */
/*																			*/
/* Purpose:   call the init routines of the low module						*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: henrik@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7007									*/
/*			  fax  : 0049-(0)711-685-7000									*/
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
#include "heaps.h"
#include "ugenv.h"
#include "fileopen.h"
#include "ugstruct.h"

/* own header */
#include "initlow.h"

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*D
   InitLow - Call the inits for the low module

   SYNOPSIS:
   INT InitLow ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function calls the inits for the low module.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

#define DEFAULTENVSIZE  64000   /* size of environment if no default value	*/

INT InitLow ()
{
  INT err;
  char buffer[BUFFSIZE];

  /* keep type for sscanf */
  int heapSize;


  /* init heaps.c */
  if ((err=InitHeaps())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init ugenv.c */
  if (GetDefaultValue(DEFAULTSFILENAME,"envmemory",buffer)==0)
    sscanf(buffer," %d ",&heapSize);
  else
    heapSize = DEFAULTENVSIZE;

  if ((err=InitUgEnv(heapSize))!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init fileopen */
  if ((err=InitFileOpen())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init structs */
  if ((err=InitUgStruct())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  return (0);
}
