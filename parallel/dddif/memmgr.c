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
#include "ppif.h"
#include "ddd.h"
#include "heaps.h"

#include "gm.h"
#include "parallel.h"


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

extern MULTIGRID *DDD_currMG;



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


static HEAP *memmgr_theHeap;
static MEM memmgr_n;

void InitMemMgr(HEAP *theHeap)
{
  memmgr_theHeap = theHeap;
}

void * GetObjMem (DDD_TYPE type, HEAP *theHeap, MEM n, int mode)
{
  INT ds, Size;

  /* TODO: untiges ist wahrscheinlich zu langsam,
                   klappt gar nicht fuer Xfer! */

  memmgr_theHeap = theHeap;        /* wieviele heaps gibs?  */
  /* TODO: check whether this is to implement
          if (type == TypeVector)
          {
                  ds = DDD_currMG->theFormat->VectorSizes[NODEVECTOR];
                  Size = sizeof(VECTOR)-sizeof(DOUBLE)+ds;
                  memmgr_n = Size;
          }
          else
   */
  memmgr_n = n;                   /* kann bis auf vector aus typ berechnet werden */


  return (DDD_ObjGetX(type,memmgr_n));
}

void *memmgr_AllocOMEM (unsigned long size, int id)
{
  void   *buffer;

  /*
     if (id==TypeBElement)
     memmgr_n = sizeof(struct belement);
   */

  printf("%4d: AllocOMem: size=%04d  id=%02d\n",me,size,id);
  buffer = GetMem(memmgr_theHeap, size, FROM_BOTTOM);

  return(buffer);
}


void memmgr_FreeOMEM (void *buffer, int id)
{
        #ifdef Debug
  printf("%d: memmgr_FreeOMEM(): buffer=%x, id=%d\n",me,buffer,id);
        #endif
  /* TODO: delete this */
  return;

  free(buffer);
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
{}

#endif /* ModelP */
