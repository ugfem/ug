// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      memmgr.c                                                      */
/*                                                                          */
/* Purpose:   memory management module                                      */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/04/27 kb  begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP

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

#include "compiler.h"

#include "parallel.h"
#include "general.h"

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

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


void *memmgr_AllocPMEM (unsigned long size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreePMEM (void *buffer)
{
  free(buffer);
}




void *memmgr_AllocOMEM (size_t size, int ddd_type, int prio, int attr)
{
  void   *buffer;

  buffer = GetFreelistMemory(MGHEAP(dddctrl.currMG), size);

  /*
     printf("%4d: memmgr_AllocOMem: size=%05d ddd_type=%02d prio=%d attr=%d\n",
     me,size,ddd_type,prio,attr);
   */

  return(buffer);
}


void memmgr_FreeOMEM (void *buffer, size_t size, int ddd_type)
{
  /*
     printf("%d: memmgr_FreeOMEM(): buffer=%x, ddd_type=%d\n", me, buffer, ddd_type);
   */

  PutFreelistMemory(MGHEAP(dddctrl.currMG), buffer, size);
}




void *memmgr_AllocAMEM (unsigned long size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeAMEM (void *buffer)
{
  free(buffer);
}


void *memmgr_AllocTMEM (unsigned long size)
{
  void   *buffer;

  buffer = malloc(size);
  return(buffer);
}


void memmgr_FreeTMEM (void *buffer)
{
  free(buffer);
}


void *memmgr_AllocHMEM (unsigned long size)
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
{
  printf("memmgr_Init aufgerufen.\n");
}



#endif /* ModelP */
