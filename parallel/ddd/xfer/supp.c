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


/****************************************************************************/
/*                                                                          */
/* definition of constants, macros                                          */
/*                                                                          */
/****************************************************************************/

#define ADDDATASEGM_SIZE 256
#define SIZESSEGM_SIZE 2048




/****************************************************************************/
/*                                                                          */
/* data types                                                               */
/*                                                                          */
/****************************************************************************/

/* segment of AddDatas */
typedef struct _AddDataSegm
{
  struct _AddDataSegm *next;
  int nItems;

  XFERADDDATA item[ADDDATASEGM_SIZE];
} AddDataSegm;


/* segment of AddData-Sizes */
typedef struct _SizesSegm
{
  struct _SizesSegm   *next;
  int current;

  int data[SIZESSEGM_SIZE];
} SizesSegm;




/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


static AddDataSegm *segmAddData = NULL;
static SizesSegm   *segmSizes   = NULL;


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
#undef T

#define T XIDelCmd
#include "sll.ct"
#undef T

#define T XIDelObj
#include "sll.ct"
#undef T

#define T XISetPrio
#include "sll.ct"
#undef T

#define T XINewCpl
#include "sll.ct"
#undef T

#define T XIOldCpl
#include "sll.ct"
#undef T



#define T XIAddCpl
#include "sll.ct"
#undef T

#define T XIDelCpl
#include "sll.ct"
#undef T

#define T XIModCpl
#include "sll.ct"
#undef T



/****************************************************************************/


static AddDataSegm *NewAddDataSegm (void)
{
  AddDataSegm *segm;

  segm = (AddDataSegm *) AllocTmp(sizeof(AddDataSegm));

  if (segm==NULL)
  {
    DDD_PrintError('F', 9999, STR_NOMEM " during XferEnd()");
    HARD_EXIT;
  }

  segm->next   = segmAddData;
  segmAddData  = segm;
  segm->nItems = 0;

  return(segm);
}


static void FreeAddDataSegms (void)
{
  AddDataSegm *segm = segmAddData;
  AddDataSegm *next = NULL;

  while (segm!=NULL)
  {
    next = segm->next;
    FreeTmp(segm);

    segm = next;
  }

  segmAddData = NULL;
}


/****************************************************************************/


static SizesSegm *NewSizesSegm (void)
{
  SizesSegm *segm;

  segm = (SizesSegm *) AllocTmp(sizeof(SizesSegm));

  if (segm==NULL)
  {
    DDD_PrintError('F', 9999, STR_NOMEM " during XferEnd()");
    HARD_EXIT;
  }

  segm->next    = segmSizes;
  segmSizes     = segm;
  segm->current = 0;

  return(segm);
}


static void FreeSizesSegms (void)
{
  SizesSegm *segm = segmSizes;
  SizesSegm *next = NULL;

  while (segm!=NULL)
  {
    next = segm->next;
    FreeTmp(segm);

    segm = next;
  }

  segmSizes = NULL;
}


/****************************************************************************/


XFERADDDATA *NewXIAddData (void)
{
  AddDataSegm *segm = segmAddData;
  XFERADDDATA *xa;

  if (segm==NULL || segm->nItems==ADDDATASEGM_SIZE)
  {
    segm = NewAddDataSegm();
  }

  xa = &(segm->item[segm->nItems++]);
  xa->next = theXIAddData->add;
  theXIAddData->add = xa;

  return(xa);
}



void FreeAllXIAddData (void)
{
  FreeAddDataSegms();
  FreeSizesSegms();
}


/****************************************************************************/

int *AddDataAllocSizes (int cnt)
{
  SizesSegm *segm = segmSizes;
  int *pos;

  if (segm==NULL || segm->current+cnt>=SIZESSEGM_SIZE)
  {
    segm = NewSizesSegm();
  }

  pos = segm->data + segm->current;
  segm->current += cnt;

  return(pos);
}


/****************************************************************************/
