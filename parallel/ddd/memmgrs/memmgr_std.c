// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      memmgr_std.c                                                  */
/*                                                                          */
/* Purpose:   basic memory management module                                */
/*            (with standard malloc() calls)                                */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/04/27 kb  begin                                            */
/*            96/01/20 kb  updated to DDD V1.5                              */
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

/* standard C library */
#include <stdlib.h>
#include <stdio.h>

#include "ppif.h"


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


void *memmgr_AllocPMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreePMEM (void *buffer)
{
  free(buffer);
}




void *memmgr_AllocOMEM (size_t size, int ddd_typ, int proc, int attr)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeOMEM (void *buffer, size_t size, int ddd_typ)
{
  free(buffer);
}




void *memmgr_AllocAMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeAMEM (void *buffer)
{
  free(buffer);
}


void *memmgr_AllocTMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeTMEM (void *buffer)
{
  free(buffer);
}


void *memmgr_AllocHMEM (size_t size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeHMEM (void *buffer)
{
  free(buffer);
}


void memmgr_MarkHMEM (void)
{}


void memmgr_ReleaseHMEM (void)
{}



void memmgr_Init (void)
{}
