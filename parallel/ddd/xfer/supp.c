// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      supp.c                                                        */
/*                                                                          */
/* Purpose:   support routines for Transfer Module                          */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin (xfer.c)                                   */
/*            95/03/21 kb  added variable sized objects (XferCopyObjX)      */
/*            95/04/05 kb  V1.3: extracted from xfer.c                      */
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


#include "dddi.h"
#include "xfer.h"



/* defined in cmds.c */
extern XICopyObj *theXIAddData;



XFERADDDATA *freeXIAddData;



/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/*
    include templates
 */
#define T XICopyObj
#include "sll.ct"

#define T XIDelCmd
#include "sll.ct"

#define T XIDelObj
#include "sll.ct"

#define T XISetPrio
#include "sll.ct"

#define T XINewCpl
#include "sll.ct"

#define T XIOldCpl
#include "sll.ct"



#define T XIAddCpl
#include "sll.ct"

#define T XIDelCpl
#include "sll.ct"

#define T XIModCpl
#include "sll.ct"



/****************************************************************************/


XFERADDDATA *NewXIAddData (void)
{
  XFERADDDATA *xa;

  if (freeXIAddData==NULL)
  {
    xa = (XFERADDDATA *) AllocTmp(sizeof(XFERADDDATA));
  }
  else
  {
    xa = freeXIAddData;
    freeXIAddData = xa->next;
  }

  xa->next = theXIAddData->add;
  theXIAddData->add = xa;

  return(xa);
}


XFERADDDATA *FreeXIAddData (XFERADDDATA *item)
{
  XFERADDDATA   *next = item->next;

  if (item->sizes!=NULL)
  {
    FreeTmp(item->sizes);
    item->sizes=NULL;
  }

  item->next = freeXIAddData;
  freeXIAddData = item;

  return(next);
}



void FreeAllXIAddData (void)
{
  XICopyObj *xi;

  /* for all XICopyObj-items */
  for(xi=listXICopyObj; xi!=NULL; xi=xi->sll_next)
  {
    XFERADDDATA *xa;

    /* free XferAdd items, if necessary */
    xa = xi->add;
    while (xa!=NULL)
    {
      xa = FreeXIAddData(xa);
    }
  }
}


/****************************************************************************/
