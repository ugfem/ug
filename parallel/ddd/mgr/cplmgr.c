// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      cplmgr.c                                                      */
/*                                                                          */
/* Purpose:   management of couplings                                       */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            94/08/24 kb  added DDD_InfoProcList()                         */
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



#define DebugCoupling 10  /* 10 is off */



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
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



static COUPLING *memlistCpl;



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/

static COUPLING *NewCoupling (void)
{
  COUPLING *cpl;

  if (memlistCpl==NULL)
  {
    cpl = (COUPLING *) AllocCpl(sizeof(COUPLING));
  }
  else
  {
    cpl = memlistCpl;
    memlistCpl = cpl->next;
  }

  nCplItems++;

  return(cpl);
}


static void DisposeCoupling (COUPLING *cpl)
{
  cpl->next = memlistCpl;
  memlistCpl = cpl;

  nCplItems--;
}



/****************************************************************************/
/*                                                                          */
/* Function:  AddCoupling                                                   */
/*                                                                          */
/* Purpose:   get new coupling record and init contents                     */
/*            if coupling is already there, no additional coupling is       */
/*            created. priority is adapted in this case.                    */
/*                                                                          */
/* Input:     hdr: DDD-header of object with new coupling                   */
/*            proc: owner of copy to be registered                          */
/*            prio: priority of copy                                        */
/*                                                                          */
/* Output:    ptr to new cpl record (or old one, if existed before)         */
/*            NULL on error                                                 */
/*                                                                          */
/****************************************************************************/

COUPLING *AddCoupling (DDD_HDR hdr, DDD_PROC proc, DDD_PRIO prio)
{
  COUPLING        *cp, *cp2;
  DDD_HDR oldObj;
  int objIndex;

#       if DebugCoupling<=1
  sprintf(cBuffer, "%4d: AddCoupling %08x proc=%d prio=%d\n",
          me, OBJ_GID(hdr), proc, prio);
  DDD_PrintDebug(cBuffer);
#       endif

  /* find or free position in coupling array */
  objIndex = OBJ_INDEX(hdr);
  if (! HAS_COUPLING(hdr))
  {
    if (nCpls==MAX_CPL) {
      DDD_PrintError('E', 2520, "no more couplings in AddCoupling");
      return(NULL);
    }

    oldObj = theObj[nCpls];

    /* exchange object without coupling and object with coupling */
    theObj[objIndex] = oldObj;
    OBJ_INDEX(oldObj) = objIndex;
    theObj[nCpls]    = hdr;
    OBJ_INDEX(hdr)    = nCpls;

    objIndex = nCpls;
    theCpl[objIndex] = NULL;
    theCplN[objIndex] = 0;
    nCpls++;
  }
  else
  {
    for(cp2=theCpl[objIndex]; cp2!=NULL; cp2=cp2->next)
    {
      if (cp2->proc==proc)
      {
        if (cp2->prio!=prio)
        {
          /* coupling upgrades/downgrades, are they allowed?
                                                  printf("%4d: diff in cpl, %05x old %d-%d new %d-%d\n",
                                                          me,OBJ_GID(hdr),cp2->proc,cp2->prio, proc, prio);
           */
          cp2->prio = prio;
        }
        /*
                                        DDD_PrintError('W', 2600, "coupling already known in AddCoupling");
         */
        return(cp2);
      }
    }
  }

  /* create new coupling record */
  cp = NewCoupling();
  if (cp==NULL) {
    DDD_PrintError('E', 2500, "not enough memory in AddCoupling");
    return(NULL);
  }

  /* init contents */
  cp->obj = hdr;
  cp->proc = proc;
  cp->prio = prio;

  /* insert into theCpl array */
  cp->next = theCpl[objIndex];
  theCpl[objIndex] = cp;
  theCplN[objIndex]++;

  return(cp);
}





/****************************************************************************/
/*                                                                          */
/* Function:  ModCoupling                                                   */
/*                                                                          */
/* Purpose:   find existing coupling record and modify priority             */
/*            this function does coupling upgrade/downgrade without         */
/*            complaining.                                                  */
/*                                                                          */
/* Input:     hdr: DDD-header of object with new coupling                   */
/*            proc: owner of copy to be modified                            */
/*            prio: new priority of copy                                    */
/*                                                                          */
/* Output:    ptr to old cpl record                                         */
/*            NULL on error                                                 */
/*                                                                          */
/****************************************************************************/

COUPLING *ModCoupling (DDD_HDR hdr, DDD_PROC proc, DDD_PRIO prio)
{
  COUPLING        *cp, *cp2;
  DDD_HDR oldObj;
  int objIndex;

#       if DebugCoupling<=1
  sprintf(cBuffer, "%4d: ModCoupling %08x proc=%d prio=%d\n",
          me, OBJ_GID(hdr), proc, prio);
  DDD_PrintDebug(cBuffer);
#       endif

  /* find or free position in coupling array */
  objIndex = OBJ_INDEX(hdr);
  if (! HAS_COUPLING(hdr))
  {
    /* there are no couplings for this object! */
    sprintf(cBuffer, "no couplings for %08x in ModCoupling", OBJ_GID(hdr));
    DDD_PrintError('E', 2530, cBuffer);
    return(NULL);
  }
  else
  {
    /* look if coupling exists and change it */
    for(cp2=theCpl[objIndex]; cp2!=NULL; cp2=cp2->next)
    {
      if (cp2->proc==proc)
      {
        cp2->prio = prio;
        return(cp2);
      }
    }
  }

  /* coupling not found */
  sprintf(cBuffer, "no coupling from %d for %08x in ModCoupling",
          proc, OBJ_GID(hdr));
  DDD_PrintError('E', 2531, cBuffer);
  return(NULL);
}




/****************************************************************************/
/*                                                                          */
/* Function:  DelCoupling                                                   */
/*                                                                          */
/* Purpose:   remove coupling record from object                            */
/*                                                                          */
/* Input:     hdr: DDD-header of object with old coupling                   */
/*            proc: owner of copy to be un-registered                       */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DelCoupling (DDD_HDR hdr, DDD_PROC proc)
{
  COUPLING        *cpl, *cplLast;
  int objIndex;

  objIndex = OBJ_INDEX(hdr);

  if (objIndex<nCpls)
  {
    for(cpl=theCpl[objIndex], cplLast=NULL; cpl!=NULL; cpl=cpl->next)
    {
      if(cpl->proc==proc)
      {
        if (cplLast==NULL) {
          theCpl[objIndex] = cpl->next;
        } else {
          cplLast->next = cpl->next;
        }
#                               if DebugCoupling<=1
        sprintf(cBuffer,"%4d: DelCoupling %07x on proc=%d, now %d cpls\n",
                me, OBJ_GID(hdr), proc, theCplN[objIndex]-1);
        DDD_PrintDebug(cBuffer);
#                               endif

        DisposeCoupling(cpl);

        theCplN[objIndex]--;

        if (theCplN[objIndex]==0) {
          nCpls--;
          OBJ_INDEX(hdr) = nCpls;
          OBJ_INDEX(theObj[nCpls]) = objIndex;
          theObj[objIndex] = theObj[nCpls];
          theObj[nCpls] = hdr;
          theCpl[objIndex] = theCpl[nCpls];
          theCplN[objIndex] = theCplN[nCpls];
        }
        break;
      }
      cplLast = cpl;
    }
  }
}


/****************************************************************************/
/*                                                                          */
/* Function:  DisposeCouplingList                                           */
/*                                                                          */
/* Purpose:   dispose complete coupling list                                */
/*                                                                          */
/* Input:     cpl: first element of coupling list                           */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DisposeCouplingList (COUPLING *cpl)
{
  COUPLING *c, *next;

  c = cpl;
  while (c!=NULL)
  {
    next = c->next;
    DisposeCoupling(c);
    c = next;
  }
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoProcList                                              */
/*                                                                          */
/* Purpose:   return list of couplings of certain object                    */
/*                                                                          */
/* Input:     hdr: DDD-header of object with coupling                       */
/*                                                                          */
/* Output:    pointer to iBuffer, which has been filled with:               */
/*               1) id of calling processor                                 */
/*               2) priority of local object coppy on calling processor     */
/*               3) id of processor which holds a object copy               */
/*               4) priority of copy on that processor                      */
/*               5) 3+4 repeated for each coupling                          */
/*               6) processor number = -1 as end mark                       */
/*                                                                          */
/****************************************************************************/

int *DDD_InfoProcList (DDD_HDR hdr)
{
  COUPLING *cpl;
  int i, objIndex = OBJ_INDEX(hdr);

  /* insert description of own (i.e. local) copy */
  iBuffer[0] = me;
  iBuffer[1] = OBJ_PRIO(hdr);

  i=2;

  /* append descriptions of foreign copies */
  if (objIndex<nCpls) {
    for(cpl=theCpl[objIndex]; cpl!=NULL; cpl=cpl->next, i+=2) {
      iBuffer[i]   = cpl->proc;
      iBuffer[i+1] = cpl->prio;
    }
  }

  /* append end mark */
  iBuffer[i] = -1;

  return(iBuffer);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoProcPrio                                              */
/*                                                                          */
/* Purpose:   return first processor number with a given priority           */
/*                                                                          */
/* Input:     hdr:  DDD-header of object with coupling                      */
/*            prio: priority to search for                                  */
/*                                                                          */
/* Output:    id of processor which holds the object copy with prio         */
/*            (or procs if no such copy exists)                             */
/*                                                                          */
/****************************************************************************/

DDD_PROC DDD_InfoProcPrio (DDD_HDR hdr, DDD_PRIO prio)
{
  COUPLING *cpl;
  int objIndex = OBJ_INDEX(hdr);

  /* append descriptions of foreign copies */
  if (objIndex<nCpls) {
    for(cpl=theCpl[objIndex]; cpl!=NULL; cpl=cpl->next)
    {
      if (cpl->prio == prio)
        return(cpl->proc);
    }
  }

  /* eventually local copy has priority we are looking for */
  if (OBJ_PRIO(hdr)==prio)
    return(me);

  return(procs);
}

int DDD_InfoIsLocal (DDD_HDR hdr)
{
  return(! HAS_COUPLING(hdr));
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoCoupling                                              */
/*                                                                          */
/* Purpose:   displays list of coupling for certain object                  */
/*                                                                          */
/* Input:     hdr: DDD-header of object with coupling                       */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DDD_InfoCoupling (DDD_HDR hdr)
{
  COUPLING *cpl;
  int objIndex = OBJ_INDEX(hdr);

  sprintf(cBuffer, "%4d: InfoCoupling for object %07x (%05d/%05d)\n",
          me, OBJ_GID(hdr), objIndex, nCpls);
  DDD_PrintLine(cBuffer);

  if (objIndex<nCpls) {
    for(cpl=theCpl[objIndex]; cpl!=NULL; cpl=cpl->next) {
      sprintf(cBuffer, "%4d:    cpl %08x proc=%4d prio=%4d\n",
              me, cpl, cpl->proc, cpl->prio);
      DDD_PrintLine(cBuffer);
    }
  }
}



/****************************************************************************/
/*                                                                          */
/* Function:  CplMgrInit and CplMgrExit                                     */
/*                                                                          */
/* Purpose:   init/exit coupling manager                                    */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void ddd_CplMgrInit (void)
{
  memlistCpl = NULL;
}


void ddd_CplMgrExit (void)
{
  /* TODO put freeing of memlist of unused COUPLINGS here */
}
