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
#include <assert.h>


#include "dddi.h"


/*
        NOTE: all container-classes from ooppcc.h are implemented in this
              source file by setting the following define.
 */
#define ContainerImplementation
/* this is the hardliner version, for debugging only
   #define _CHECKALLOC(ptr)   assert(ptr!=NULL)
 */
#define _CHECKALLOC(ptr)   if (ptr==NULL) return (NULL)


static int TmpMem_kind = TMEM_ANY;

void *xfer_AllocTmp (size_t size)
{
  void *buffer = AllocTmpReq(size, TmpMem_kind);
  return(buffer);
}

void xfer_FreeTmp (void *buffer)
{
  FreeTmpReq(buffer, 0, TmpMem_kind);
}

void xfer_SetTmpMem (int kind)
{
  TmpMem_kind = kind;
}



/* forward declaration */
void *xfer_AllocHeap (size_t);
void xfer_FreeHeap (void *);


#include "xfer.h"



#ifdef XferMemFromHeap
void *xfer_AllocHeap (size_t size)
{
  void *buffer;

  if (xferGlobals.useHeap)
  {
    buffer = AllocHeap(size, xferGlobals.theMarkKey);
  }
  else
  {
    buffer = AllocTmp(size);
  }

  return(buffer);
}

void xfer_FreeHeap (void *buffer)
{
  if (!xferGlobals.useHeap)
  {
    FreeTmp(buffer,0);
  }
  /* else: do nothing for heap-allocated memory */
}

#endif


void *xfer_AllocSend (size_t size)
{
  void *buffer = AllocTmpReq(size, TMEM_ANY);
  return(buffer);
}

void xfer_FreeSend (void *buffer)
{
  FreeTmpReq(buffer, 0, TMEM_ANY);
}



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
/* class member function implementations                                    */
/*                                                                          */
/****************************************************************************/

#define ClassName XICopyObj

/*
        compare-method in order to eliminate double XICopyObj-items.
        merge priorities from similar XICopyObj-items.

        the items are sorted according to key (dest,gid),
        all in ascending order. if dest and gid are equal,
        we merge priorities and get a new priority together with
        the information whether first item wins over second.
        in both cases, we use the new priority for next comparison.

        this implements rule XFER-C1.
 */
int Method(Compare) (ClassPtr item1, ClassPtr item2)
{
  DDD_PRIO newprio;
  int ret;

  if (item1->dest < item2->dest) return(-1);
  if (item1->dest > item2->dest) return(1);

  if (item1->gid < item2->gid) return(-1);
  if (item1->gid > item2->gid) return(1);


  /* items have equal gid and dest, so they are considered as equal. */
  /* however, we must check priority, and patch both items with
     the new priority after merge. */
  ret = PriorityMerge(&theTypeDefs[OBJ_TYPE(item1->hdr)],
                      item1->prio, item2->prio, &newprio);

  item1->prio = newprio;

  if (ret==PRIO_FIRST || ret==PRIO_UNKNOWN)
  {
    /* tell XferInitCopyInfo() that item is rejected */
    item2->prio = PRIO_INVALID;
  }
  else
  {
    /* communicate new priority to XferInitCopyInfo() */
    item2->prio = newprio;
  }

  return(0);
}


void Method(Print) (ParamThis _PRINTPARAMS)
{
  fprintf(fp, "XICopyObj dest=%d gid=%08x prio=%d\n",
          This->dest, This->gid, This->prio);
}

#undef ClassName



/****************************************************************************/

#define ClassName XISetPrio

/*
        compare-method in order to eliminate double XISetPrio-items.
        merge priorities from similar XISetPrio-items.

        the items are sorted according to key (gid),
        all in ascending order. if both gids are equal,
        we merge priorities and get a new priority together with
        the information whether first item wins over second.
        in both cases, we use the new priority for next comparison.

        this implements rule XFER-P1.
 */
int Method(Compare) (ClassPtr item1, ClassPtr item2)
{
  DDD_PRIO newprio;
  int ret;

  /* ascending GID is needed for ExecLocalXI___ */
  if (item1->gid < item2->gid) return(-1);
  if (item1->gid > item2->gid) return(1);


  /* items have equal gid and dest, so they are considered as equal. */
  /* however, we must check priority, and patch both items with
     the new priority after merge. */
  ret = PriorityMerge(&theTypeDefs[OBJ_TYPE(item1->hdr)],
                      item1->prio, item2->prio, &newprio);

  item1->prio = item2->prio = newprio;


  if (ret==PRIO_FIRST || ret==PRIO_UNKNOWN)
  {
    /* tell XferInitCopyInfo() that item is rejected */
    item2->prio = PRIO_INVALID;
  }
  else
  {
    /* communicate new priority to XferInitCopyInfo() */
    item2->prio = newprio;
  }

  return(0);
}


void Method(Print) (ParamThis _PRINTPARAMS)
{
  fprintf(fp, "XISetPrio gid=%08x prio=%d\n", This->gid, This->prio);
}

#undef ClassName


/****************************************************************************/


/*
    include templates
 */
#define T XIDelCmd
#define SLL_WithOrigOrder
#include "sll.ct"
#undef T

#define T XIDelObj
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

  segm = (AddDataSegm *) OO_Allocate(sizeof(AddDataSegm));
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
    OO_Free (segm /*,sizeof(AddDataSegm)*/);

    segm = next;
  }

  segmAddData = NULL;
}


/****************************************************************************/


static SizesSegm *NewSizesSegm (void)
{
  SizesSegm *segm;

  segm = (SizesSegm *) OO_Allocate (sizeof(SizesSegm));
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
    OO_Free (segm /*,sizeof(SizesSegm)*/);

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


/*
        get quantitative resource usage
 */
void GetSizesXIAddData (int *nSegms, int *nItems, size_t *alloc_mem, size_t *used_mem)
{
  size_t allocated=0, used=0;
  int ns=0, ni=0;

  {
    AddDataSegm  *segm;
    for (segm=segmAddData; segm!=NULL; segm=segm->next)
    {
      /* count number of segments and number of items */
      ns++;
      ni+=segm->nItems;

      /* compute memory usage */
      allocated += sizeof(AddDataSegm);
      used += (sizeof(AddDataSegm) -
               (sizeof(XFERADDDATA)*(ADDDATASEGM_SIZE-segm->nItems)));
    }
  }

  /* TODO: add resources for SizesSegm */


  *nSegms    = ns;
  *nItems    = ni;
  *alloc_mem = allocated;
  *used_mem  = used;
}


/****************************************************************************/
