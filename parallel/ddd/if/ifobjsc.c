// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ifobjsc.c                                                     */
/*                                                                          */
/* Purpose:   routines concerning interfaces between processors             */
/*            part 4: routines for creation and management of               */
/*                    object shortcut tables                                */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   96/01/24 kb  copied from ifcreate.c                           */
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
#include "if.h"





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
        convert cpl-IF-table into obj-IF-table
 */
static void IFComputeShortcutTable (DDD_IF ifId)
{
  int nItems = theIF[ifId].nItems;
  COUPLING  **cpls = theIF[ifId].cpl;
  IFObjPtr   *objs = theIF[ifId].obj;
  int i;

  /* mark obj-shortcut-table as valid */
  theIF[ifId].objValid = TRUE;

  if (nItems==0)
    return;

  /* fill in object pointers, this is the 4-fold indirection step */
  for(i=0; i<nItems; i++)
  {
#ifdef CPP_FRONTEND
    // TODO, avoid dirty cast!
    objs[i] = (DDD_Object*)(cpls[i]->obj);
#else
    objs[i] = OBJ_OBJ(cpls[i]->obj);
#endif
    /*
       #ifdef F_FRONTEND
       printf("%4d: Shortcut IF=%d item=%d/%d %08x %d\n",
            me, ifId, i, nItems, OBJ_GID(cpls[i]->obj), objs[i]);
       #endif
     */
  }
}


/****************************************************************************/


/*
        create direct link from ifHead and ifAttr to objects,
        avoid one indirect addressing step across couplings.
        each cpl-entry in an interface has one corresponding obj-entry
 */
void IFCreateObjShortcut (DDD_IF ifId)
{
  COUPLING    **cplarray = theIF[ifId].cpl;
  IFObjPtr     *objarray;
  IF_PROC     *ifHead;

  /* dont create shortcuts for STD_INTERFACE */
  if (ifId==STD_INTERFACE)
    return;

  /* are there any items? */
  if (theIF[ifId].nItems == 0)
    return;

  /* get memory for addresses of objects inside IF */
  objarray = (IFObjPtr *) AllocIF(sizeof(IFObjPtr)*theIF[ifId].nItems);
  if (objarray==NULL) {
    DDD_PrintError('E', 4000, STR_NOMEM " in IFCreateObjShortcut");
    HARD_EXIT;
  }
  theIF[ifId].obj = objarray;

  IFComputeShortcutTable(ifId);


  ForIF(ifId,ifHead)
  {
    IF_ATTR  *ifAttr;

    /* compute pointers to subarrays */
    ifHead->obj    = objarray + (ifHead->cpl    - cplarray);
    ifHead->objAB  = objarray + (ifHead->cplAB  - cplarray);
    ifHead->objBA  = objarray + (ifHead->cplBA  - cplarray);
    ifHead->objABA = objarray + (ifHead->cplABA - cplarray);


    /* compute pointers from ifAttrs to subarrays */
    for(ifAttr=ifHead->ifAttr; ifAttr!=0; ifAttr=ifAttr->next)
    {
      ifAttr->objAB  = objarray + (ifAttr->cplAB  - cplarray);
      ifAttr->objBA  = objarray + (ifAttr->cplBA  - cplarray);
      ifAttr->objABA = objarray + (ifAttr->cplABA - cplarray);
    }
  }
}


/****************************************************************************/


/*
        if object addresses in memory are changed, then the shortcut-tables
        will get invalid. this routine does the invalidation.
 */
void IFInvalidateShortcuts (DDD_TYPE invalid_type)
{
  int i;

  /* test all interfaces */
  for(i=0; i<nIFs; i++)
  {
    if (i==STD_INTERFACE)
      continue;

    if (theIF[i].objValid)
    {
      /* determine whether object belongs to IF */
      if ((1<<invalid_type) & theIF[i].maskO)
      {
        /* yes, invalidate interface */
        theIF[i].objValid = FALSE;
      }
    }
  }
}



/****************************************************************************/


/*
        check if shortcut-table is valid and recompute, if necessary
 */
void IFCheckShortcuts (DDD_IF ifId)
{
  if (ifId==STD_INTERFACE)
    return;

  if (! theIF[ifId].objValid)
  {
    /*
       sprintf(cBuffer, "%04d: IFComputeShortcutTable IF=%d\n", me, ifId);
       DDD_PrintDebug(cBuffer);
     */

    IFComputeShortcutTable(ifId);
  }
}



/****************************************************************************/
