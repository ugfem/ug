// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      xfer.c                                                        */
/*                                                                          */
/* Purpose:   main module for object transfer                               */
/*            (contains basic functionality used by the rest of the source  */
/*             files in the Xfer-module)                                    */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            95/03/21 kb  added variable sized objects (XferCopyObjX)      */
/*            96/07/03 kb  splitted XferInfo-list into ObjXfer and CplXfer  */
/*            960718 kb  introduced lowcomm-layer (sets of messages)        */
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
#include <string.h>


#include "dddi.h"
#include "xfer.h"


/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/

/* helpful macros for FRONTEND switching, will be #undef'd at EOF */
#ifdef F_FRONTEND
#define _FADR     &
#else
#define _FADR
#endif


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



XFER_GLOBALS xferGlobals;


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)





/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



static int sort_NewOwners (const void *e1, const void *e2)
{
  REGISTER XICopyObj *item1 = *((XICopyObj **)e1);
  REGISTER XICopyObj *item2 = *((XICopyObj **)e2);

  if (item1->gid < item2->gid) return(-1);
  if (item1->gid > item2->gid) return(1);

  if (item1->dest < item2->dest) return(-1);
  if (item1->dest > item2->dest) return(1);

  return(0);
}


/****************************************************************************/

/*
        collect a temporary list of XINewCpl-items.

        for each XICopyObj-command whose destination doesn't have
        an object copy already do: create a set of XINewCpl-items, one
        for every processor which owns a copy of the local object.

        this is an estimate, because without communication the sending
        processor cannot know whether the object copy will be accepted.
        this final information will be transferred in a second pass, as
        soon as the receiver has decided whether he accepts the incoming
        object or not (depending on rules XFER-C2, XFER-C3, XFER-C4,
        XFER-P and XFER-D).
 */

XICopyObj **CplClosureEstimate (XICopyObjPtrArray *arrayItems, int *nRet)
{
  int i, nNewOwners;
  XICopyObj **arrayNewOwners = NULL;
  XICopyObj **items = XICopyObjPtrArray_GetData(arrayItems);
  int n       = XICopyObjPtrArray_GetSize(arrayItems);



  nNewOwners=0;
  for(i=0; i<n; i++)
  {
    REGISTER XICopyObj *xi = items[i];
    REGISTER DDD_PROC dest = xi->dest;              /* destination proc */
    REGISTER COUPLING *cpl, *xicpl = ObjCplList(xi->hdr);
    REGISTER DDD_GID xigid = xi->gid;
    REGISTER DDD_TYPE xitype = OBJ_TYPE(xi->hdr);

    SET_CO_NEWOWNER(xi);

    /* look if there's a coupling for dest */
    for(cpl=xicpl; cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      if (dest==CPL_PROC(cpl))
      {
        /* got one coupling, destination is not a new owner */
        CLEAR_CO_NEWOWNER(xi);

        /* destination proc had a copy before xfer */
        /* check whether priority of that copy will change */
        /*
                                        if (WhichPrioWins(xi->prio, cpl->prio)==1)
                                        {
         */
        /* new prio will win on other proc -> adapt coupling */
        /*
           printf("%4d: XXXCoupling %08x proc=%d prio=%d\n",
           me,xigid,dest,xi->prio);
                                                cpl->prio = xi->prio;
                                        }
         */

        /* this should be the only coupling for that proc.
           leave loop prematurely. */
        break;
      }
    }


    if (CO_NEWOWNER(xi))
    {
      nNewOwners++;

      /* destination proc didn't have a copy before xfer */

      /* inform other owners of local copies (XINewCpl) */
      for(cpl=xicpl; cpl!=NULL; cpl=CPL_NEXT(cpl))
      {
        XINewCpl *xc = NewXINewCpl(SLLNewArgs);
        if (xc==NULL)
          HARD_EXIT;

        xc->to      = CPL_PROC(cpl);                         /* receiver of XINewCpl    */
        NewCpl_SetDest(xc->te,dest);                         /* destination of XICopyObj*/
        NewCpl_SetGid(xc->te,xigid);                         /* the object's gid        */
        NewCpl_SetPrio(xc->te,xi->prio);                         /* new obj's priority  */
        NewCpl_SetType(xc->te,xitype);                           /* the object's type   */
      }

      /* send current couplings (XIOldCpl) to new destination */
      /* note: destination proc can get this information
               multiple times, once for each incoming object
               with same gid (from different senders)  */
      for(cpl=xicpl; cpl!=NULL; cpl=CPL_NEXT(cpl))
      {
        XIOldCpl *xc = NewXIOldCpl(SLLNewArgs);
        if (xc==NULL)
          HARD_EXIT;

        xc->to      = dest;                                    /* receiver of XIOldCpl */
        xc->te.gid  = xigid;                                   /* the object's gid     */
        xc->te.proc = CPL_PROC(cpl);                           /* coupling proc        */
        xc->te.prio = cpl->prio;                               /* coupling priority    */
      }

      /* send one coupling (XIOldCpl) for local copy */
      {
        XIOldCpl *xc = NewXIOldCpl(SLLNewArgs);
        if (xc==NULL)
          HARD_EXIT;

        xc->to      = dest;                                     /* receiver of XIOldCpl */
        xc->te.gid  = xigid;                                    /* the object's gid     */
        xc->te.proc = me;                                       /* coupling proc        */
        xc->te.prio = OBJ_PRIO(xi->hdr);                        /* coupling priority    */
      }
    }
  }

  *nRet = nNewOwners;
  /*
     printf("%4d: nNewOwners=%d\n", me, nNewOwners);
   */

  /* check multiple new-owner-destinations for same gid */
  if (nNewOwners>0)
  {
    int j, k;

    arrayNewOwners =
      (XICopyObj **) OO_Allocate (sizeof(XICopyObj *)* nNewOwners);
    if (arrayNewOwners==NULL)
    {
      DDD_PrintError('E', 6102, STR_NOMEM " in XferEnd()");
      return NULL;
    }

    /* fill pointer array XICopyObj-items marked CO_NEWOWNER */
    for(j=0, k=0; j<n; j++)
    {
      if (CO_NEWOWNER(items[j]))
      {
        arrayNewOwners[k] = items[j];
        k++;
      }
    }

    if (nNewOwners==1)
      return(arrayNewOwners);

    /* sort according to gid (items is sorted according to dest) */
    qsort(arrayNewOwners, nNewOwners, sizeof(XICopyObj *), sort_NewOwners);


    for(j=0; j<nNewOwners-1; j++)
    {
      REGISTER XICopyObj *no1 = arrayNewOwners[j];
      REGISTER DDD_GID gid1 = no1->gid;

      for(k=j+1; k<nNewOwners; k++)
      {
        REGISTER XICopyObj *no2 = arrayNewOwners[k];
        REGISTER DDD_TYPE no2type;

        if (no2->gid != gid1)
          break;

        no2type = OBJ_TYPE(no2->hdr);

        /* inform other new-owners of same obj (also XINewCpl!)    */
        /* tell no1-dest that no2-dest gets a copy with no2->prio  */
        {
          XINewCpl *xc = NewXINewCpl(SLLNewArgs);
          if (xc==NULL)
            HARD_EXIT;

          xc->to      = no1->dest;                                 /* receiver of XINewCpl     */
          NewCpl_SetDest(xc->te,no2->dest);                               /* dest of XICopyObj */
          NewCpl_SetGid(xc->te,gid1);                                     /* the obj's gid     */
          NewCpl_SetPrio(xc->te,no2->prio);                               /* new obj's priority*/
          NewCpl_SetType(xc->te,no2type);                                 /* the obj's type    */
        }
        /* tell no2->dest that no1-dest gets a copy with no1->prio */
        {
          XINewCpl *xc = NewXINewCpl(SLLNewArgs);
          if (xc==NULL)
            HARD_EXIT;

          xc->to      = no2->dest;                                 /* receiver of XINewCpl     */
          NewCpl_SetDest(xc->te,no1->dest);                               /* dest of XICopyObj */
          NewCpl_SetGid(xc->te,gid1);                                     /* the obj's gid     */
          NewCpl_SetPrio(xc->te,no1->prio);                               /* new obj's priority*/
          NewCpl_SetType(xc->te,no2type);                                 /* the obj's type    */
        }
      }
    }
  }

  return(arrayNewOwners);
}


/****************************************************************************/

/*
        auxiliary functions for PrepareObjMsgs()
 */


static void BuildDepDataInfo (XFERMSG *xm, XICopyObj *xi)
{
  XFERADDDATA *xa;
  int ptr, chunks;


  /* count characteristic values for each chunk */
  chunks = ptr = 0;
  for(xa=xi->add; xa!=NULL; xa=xa->next)
  {
    ptr += xa->addNPointers;

    /* add control information size for var-sized AddData-Items */
    if (xa->sizes!=NULL)
      xi->addLen += CEIL(sizeof(int) * xa->addCnt);

    chunks++;
  }


  /* add size of control information */
  if (xi->addLen>0)
    xi->addLen += CEIL(sizeof(int)) + chunks*CEIL(2*sizeof(int));

  /* add to current message size information */
  xm->size      += xi->addLen;
  xm->nPointers += ptr;
}


static XFERMSG *CreateXferMsg (DDD_PROC dest, XFERMSG *lastxm)
{
  XFERMSG *xm;

  xm = (XFERMSG *) OO_Allocate (sizeof(XFERMSG));
  if (xm==NULL)
  {
    DDD_PrintError('E', 6100, STR_NOMEM " in PrepareObjMsgs");
    return NULL;
  }
  xm->nPointers  = 0;
  xm->nObjects   = 0;
  xm->proc = dest;
  xm->size = 0;

  xm->xferObjArray = NULL;
  xm->xferNewCpl   = NULL;
  xm->xferOldCpl   = NULL;
  xm->nObjItems = 0;
  xm->nNewCpl   = 0;
  xm->nOldCpl   = 0;

  xm->next = lastxm;

  return xm;
}



static XFERMSG *AccumXICopyObj (XFERMSG *currxm, int *nMsgs, int *nItems,
                                XICopyObj **items, DDD_PROC dest, int nmax)
{
  XFERMSG *xm;
  int i;

  if (currxm!=NULL && currxm->proc==dest)
  {
    /* there is a XFERMSG with correct processor number -> reuse it */
    xm = currxm;
  }
  else
  {
    /* create new XFERMSG structure */
    xm = CreateXferMsg(dest, currxm);
    (*nMsgs)++;
  }

#       if DebugXfer<=2
  sprintf(cBuffer, "%4d: PrepareObjMsgs, XferMsg proc=%d"
          " nmax=%d\n", me, dest, nmax);
  DDD_PrintDebug(cBuffer);
#       endif


  for (i=0; i<nmax && items[i]->dest==dest; i++)
  {
    REGISTER XICopyObj *xi = items[i];
    DDD_HDR hdr = xi->hdr;
    TYPE_DESC  *desc = &theTypeDefs[OBJ_TYPE(hdr)];

#               if DebugXfer<=0
    sprintf(cBuffer, "%4d: PrepareObjMsgs, proc=%d"
            " i=%d/%d (%08x)\n",
            me, dest, i, nmax, xi->gid);
    DDD_PrintDebug(cBuffer);
#               endif

    /* accumulate xfer-items in message-info */
    xm->nObjects++;

    /* length of object itself, possibly variable  */
    xm->size += CEIL(xi->size);
    xm->nPointers += desc->nPointers;

    if (xi->add != NULL)
      BuildDepDataInfo(xm, xi);
  }

  *nItems = i;
  return xm;
}




static XFERMSG *AccumXINewCpl (XFERMSG *currxm, int *nMsgs, int *nItems,
                               XINewCpl **items, DDD_PROC dest, int nmax)
{
  XFERMSG *xm;
  int i;

  if (currxm!=NULL && currxm->proc==dest)
  {
    /* there is a XFERMSG with correct processor number -> reuse it */
    xm = currxm;
  }
  else
  {
    /* create new XFERMSG structure */
    xm = CreateXferMsg(dest, currxm);
    (*nMsgs)++;
  }

#       if DebugXfer<=2
  sprintf(cBuffer, "%4d: PrepareObjMsgs, XferMsg proc=%d"
          " nmax=%d\n", me, dest, nmax);
  DDD_PrintDebug(cBuffer);
#       endif


  for (i=0; i<nmax && items[i]->to==dest; i++)
#               if DebugXfer<=0
  {
    XINewCpl *xi = items[i];
    sprintf(cBuffer, "%4d: PrepareObjMsgs, proc=%d"
            " i=%d/%d (%08x)\n",
            me, dest, i, nmax, xi->te.gid);
    DDD_PrintDebug(cBuffer);
  }
#               else
    ;
#               endif

  *nItems = i;
  return xm;
}



static XFERMSG *AccumXIOldCpl (XFERMSG *currxm, int *nMsgs, int *nItems,
                               XIOldCpl **items, DDD_PROC dest, int nmax)
{
  XFERMSG *xm;
  int i;

  if (currxm!=NULL && currxm->proc==dest)
  {
    /* there is a XFERMSG with correct processor number -> reuse it */
    xm = currxm;
  }
  else
  {
    /* create new XFERMSG structure */
    xm = CreateXferMsg(dest, currxm);
    (*nMsgs)++;
  }

#       if DebugXfer<=2
  sprintf(cBuffer, "%4d: PrepareObjMsgs, XferMsg proc=%d"
          " nmax=%d\n", me, dest, nmax);
  DDD_PrintDebug(cBuffer);
#       endif


  for (i=0; i<nmax && items[i]->to==dest; i++)
#               if DebugXfer<=0
  {
    XIOldCpl *xi = items[i];
    sprintf(cBuffer, "%4d: PrepareObjMsgs, proc=%d"
            " i=%d/%d (%08x)\n",
            me, dest, i, nmax, xi->te.gid);
    DDD_PrintDebug(cBuffer);
  }
#               else
    ;
#               endif

  *nItems = i;
  return xm;
}



/****************************************************************************/

/*
        prepare messages for phase 1.

        object copies will be sent as well as the estimated
        coupling closure from CplClosureEstimate().
 */

int PrepareObjMsgs (XICopyObjPtrArray *arrayO,
                    XINewCpl **itemsNC, int nNC,
                    XIOldCpl **itemsOC, int nOC,
                    XFERMSG **theMsgs, size_t *memUsage)
{
  XFERMSG    *xm=NULL;
  int iO, iNC, iOC, nMsgs=0;

  XICopyObj  **itemsO = XICopyObjPtrArray_GetData(arrayO);
  int nO       = XICopyObjPtrArray_GetSize(arrayO);


#       if DebugXfer<=3
  printf("%4d: PrepareObjMsgs, nXICopyObj=%d nXINewCpl=%d nXIOldCpl=%d\n",
         me, nO, nNC, nOC);
  fflush(stdout);
#       endif


  /*
          run through both itemsO and itemsNC/itemsOC simultaneously,
          each time a new proc-nr is encountered in one of these
          lists, create a new XFERMSG item.

          (the lists have been sorted according to proc-nr previously.)
   */

  iO=0; iNC=0; iOC=0;
  while (iO<nO || iNC<nNC || iOC<nOC)
  {
    int n;
    DDD_PROC pO = (iO<nO) ? itemsO[iO]->dest : procs;
    DDD_PROC pNC = (iNC<nNC) ? itemsNC[iNC]->to   : procs;
    DDD_PROC pOC = (iOC<nOC) ? itemsOC[iOC]->to   : procs;

    if (pO<=pNC && pO<=pOC && pO<procs)
    {
      xm = AccumXICopyObj(xm, &nMsgs, &n, itemsO+iO, pO, nO-iO);
      xm->xferObjArray = itemsO+iO;
      xm->nObjItems = n;
      iO += n;
    }

    if (pNC<=pO && pNC<=pOC && pNC<procs)
    {
      xm = AccumXINewCpl(xm, &nMsgs, &n, itemsNC+iNC, pNC, nNC-iNC);
      xm->xferNewCpl = itemsNC+iNC;
      xm->nNewCpl = n;
      iNC += n;
    }

    if (pOC<=pO && pOC<=pNC && pOC<procs)
    {
      xm = AccumXIOldCpl(xm, &nMsgs, &n, itemsOC+iOC, pOC, nOC-iOC);
      xm->xferOldCpl = itemsOC+iOC;
      xm->nOldCpl = n;
      iOC += n;
    }

    if (pO==procs) iO = nO;
    if (pNC==procs) iNC = nNC;
    if (pOC==procs) iOC = nOC;
  }
  *theMsgs = xm;


  /* compute brutto message size from netto message size */
  for(xm=*theMsgs; xm!=NULL; xm=xm->next)
  {
    size_t bufSize;
    xm->msg_h = LC_NewSendMsg(xferGlobals.objmsg_t, xm->proc);
    LC_SetTableSize(xm->msg_h, xferGlobals.symtab_id, xm->nPointers);
    LC_SetTableSize(xm->msg_h, xferGlobals.objtab_id, xm->nObjects);
    LC_SetTableSize(xm->msg_h, xferGlobals.newcpl_id, xm->nNewCpl);
    LC_SetTableSize(xm->msg_h, xferGlobals.oldcpl_id, xm->nOldCpl);
    LC_SetChunkSize(xm->msg_h, xferGlobals.objmem_id, xm->size);

    bufSize = LC_MsgFreeze(xm->msg_h);
    *memUsage += bufSize;

    if (DDD_GetOption(OPT_INFO_XFER) & XFER_SHOW_MEMUSAGE)
    {
      sprintf(cBuffer,
              "DDD MESG [%03d]: SHOW_MEM "
              "send msg  dest=%04d size=%010ld\n",
              me, xm->proc, (long)bufSize);
      DDD_PrintLine(cBuffer);
    }
  }


#       if DebugXfer<=3
  printf("%4d: PrepareObjMsgs, nMsgs=%d\n", me, nMsgs);
  fflush(stdout);
#       endif

  return(nMsgs);
}



/****************************************************************************/

/*
        execute SetPrio-commands and create those XIModCpl-items,
        which can be computed without knowledge of information sent by other
        procs during first message phase.
 */
void ExecLocalXISetPrio (
  XISetPrioPtrArray *arrayP,
  XIDelObj  **itemsD, int nD,
  XICopyObj  **itemsNO, int nNO)
{
  int iP, iD, iNO;
  XISetPrio **itemsP = XISetPrioPtrArray_GetData(arrayP);
  int nP       = XISetPrioPtrArray_GetSize(arrayP);

  /*
          execute SetPrio only if no corresponding DelObj exists!
   */
  for(iP=0, iD=0, iNO=0; iP<nP; iP++)
  {
    REGISTER XISetPrio *sp = itemsP[iP];
    REGISTER DDD_HDR hdr = sp->hdr;
    DDD_GID gid      = sp->gid;
    DDD_PRIO newprio  = sp->prio;
    COUPLING   *cpl;

    while ((iD<nD) && (itemsD[iD]->gid<gid))
      iD++;

    /* skip XICopyObj-items until entries for gid found */
    while (iNO<nNO && itemsNO[iNO]->gid<gid)
      iNO++;


    sp->is_valid = (! ((iD<nD) && (itemsD[iD]->gid==gid)));

    if (sp->is_valid)
    {
      /* SetPrio, but _no_ DelObj: execute SetPrio */
      DDD_TYPE typ   = OBJ_TYPE(hdr);
      TYPE_DESC  *desc = &(theTypeDefs[typ]);

      /* call application handler for changing prio of dependent objects */
      if (desc->handlerSETPRIORITY)
      {
        DDD_OBJ obj = HDR2OBJ(hdr,desc);

                                #if defined(C_FRONTEND) || defined(F_FRONTEND)
        desc->handlerSETPRIORITY(_FADR obj, _FADR newprio);
                                #endif
                                #if defined(CPP_FRONTEND)
        CallHandler(desc,SETPRIORITY) (HParam(obj) newprio);
                                #endif
      }

      /* change actual priority to new value */
      OBJ_PRIO(hdr) = newprio;


      /* generate XIModCpl-items */

      /* 1. for all existing couplings */
      for(cpl=ObjCplList(hdr); cpl!=NULL; cpl=CPL_NEXT(cpl))
      {
        XIModCpl *xc = NewXIModCpl(SLLNewArgs);
        if (xc==NULL)
          HARD_EXIT;

        xc->to      = CPL_PROC(cpl);                           /* receiver of XIModCpl  */
        xc->te.gid  = gid;                                     /* the object's gid      */
        xc->te.prio = newprio;                                 /* the object's new prio */
        xc->typ     = typ;                                     /* the object's type     */
      }
      /* 2. for all CopyObj-items with new-owner destinations */
      while (iNO<nNO && itemsNO[iNO]->gid==gid)
      {
        XIModCpl *xc = NewXIModCpl(SLLNewArgs);
        if (xc==NULL)
          HARD_EXIT;

        xc->to      = itemsNO[iNO]->dest;                        /* receiver of XIModCpl */
        xc->te.gid  = gid;                                       /* the object's gid     */
        xc->te.prio = newprio;                                  /* the object's new prio */
        xc->typ     = typ;                                     /* the object's type     */

        iNO++;
      }
    }
    /*
            else: SetPrio _and_ DelObj, SetPrio is invalid,
                  DelObj will be executed lateron
                  (this is rule XFER-M1).
     */
  }

}



/*
        execute local DelObj-commands and create those XIDelCpl-items,
        which can be computed without knowledge of information sent by other
        procs during first message phase.
 */


void ExecLocalXIDelCmd (XIDelCmd  **itemsD, int nD)
{
  int iD;
  XIDelCmd **origD;

  if (nD==0)
    return;

  /* reconstruct original order of DelObj commands */
  origD = (XIDelCmd **) OO_Allocate (sizeof(XIDelCmd *) * nD);
  if (origD==NULL)
  {
    DDD_PrintError('E', 6101, STR_NOMEM " in XferEnd()");
    HARD_EXIT;
  }

  /* copy pointer array and resort it */
  memcpy(origD, itemsD, sizeof(XIDelCmd *) * nD);
  OrigOrderXIDelCmd(origD, nD);


  /* loop in original order (order of Del-cmd issueing) */
  for(iD=0; iD<nD; iD++)
  {
    REGISTER DDD_HDR hdr = origD[iD]->hdr;
    DDD_TYPE typ   = OBJ_TYPE(hdr);
    TYPE_DESC  *desc = &(theTypeDefs[typ]);
    DDD_OBJ obj   = HDR2OBJ(hdr,desc);

                #if defined(C_FRONTEND) || defined(F_FRONTEND)
    /* do deletion */
    if (desc->handlerDELETE)
      desc->handlerDELETE(_FADR obj);
    else
    {
      /* TODO the following three calls should be collected in
         one ObjMgr function */

      /* destruct LDATA and GDATA */
      if (desc->handlerDESTRUCTOR!=NULL)
        desc->handlerDESTRUCTOR(_FADR obj);

      /* HdrDestructor will call ddd_XferRegisterDelete() */
      DDD_HdrDestructor(hdr);
      DDD_ObjDelete(obj, desc->size, typ);
    }
                #endif
                #ifdef CPP_FRONTEND
    // call destructor
    //printf("%4d: calling destructor for %08x\n", me, origD[iD]->obj);
    CallHandler(desc,DESTRUCTOR) (HParamOnly(obj));
                #endif
  }

  OO_Free (origD /*,0*/);
}




void ExecLocalXIDelObj (
  XIDelObj  **itemsD, int nD,
  XICopyObj  **itemsNO, int nNO)
{
  int iD, iNO;


  /* create XIDelCpl for all DelObj-commands (sorted acc. to gid) */
  for(iD=0, iNO=0; iD<nD; iD++)
  {
    DDD_GID gid   = itemsD[iD]->gid;


    /* skip XICopyObj-items until entries for gid found */
    while (iNO<nNO && itemsNO[iNO]->gid<gid)
      iNO++;


    /* generate XIDelCpl-items */
    /* 1. for all existing couplings, has been done during
          ddd_XferRegisterDelete. */

    /* 2. for all CopyObj-items with new-owner destinations */
    while (iNO<nNO && itemsNO[iNO]->gid==gid)
    {
      XIDelCpl *xc = NewXIDelCpl(SLLNewArgs);
      if (xc==NULL)
        HARD_EXIT;

      xc->to      = itemsNO[iNO]->dest;                  /* receiver of XIDelCpl */
      xc->prio    = PRIO_INVALID;                        /* dont remember priority   */
      xc->te.gid  = gid;                                 /* the object's gid     */


      /* we must remember couplings for eventual restoring
         (if this object is received from another proc) */
      xc->next = itemsD[iD]->delcpls;
      itemsD[iD]->delcpls = xc;

      iNO++;
    }
  }
}




/****************************************************************************/


/*
        create those XI???Cpl-items, which require knowledge of information
        sent by other procs during first message phase.
 */
void PropagateCplInfos (
  XISetPrio **itemsP, int nP,
  XIDelObj  **itemsD, int nD,
  TENewCpl  *arrayNC, int nNC)
{
  int iP, iD, iNC;

  /*
          step 1: create XIModCpl-items from SetPrio-cmds
                  (only if no DelObj-items exist)
   */
  for(iP=0, iNC=0; iP<nP; iP++)
  {
    REGISTER XISetPrio *sp = itemsP[iP];

    if (sp->is_valid)
    {
      REGISTER DDD_HDR hdr = sp->hdr;
      DDD_GID gid      = sp->gid;
      DDD_PRIO newprio  = sp->prio;

      /* skip TENewCpl-entries until one for gid found */
      while (iNC<nNC && NewCpl_GetGid(arrayNC[iNC])<gid)
        iNC++;

      /* generate additional XIModCpl-items for all valid NewCpl-items */
      while (iNC<nNC && NewCpl_GetGid(arrayNC[iNC])==gid)
      {
        XIModCpl *xc = NewXIModCpl(SLLNewArgs);
        if (xc==NULL)
          HARD_EXIT;

        /* receiver of XIModCpl */
        xc->to      = NewCpl_GetDest(arrayNC[iNC]);
        xc->te.gid  = gid;                                       /* the object's gid     */
        xc->te.prio = newprio;                                   /* the object's new prio */
        xc->typ     = OBJ_TYPE(hdr);                             /* the object's type     */

        iNC++;
      }
    }
  }



  /*
          step 2: create XIDelCpl-items from DelObj-cmds
   */
  for(iD=0, iNC=0; iD<nD; iD++)
  {
    DDD_GID gid   = itemsD[iD]->gid;

    /* skip TENewCpl-entries until one for gid found */
    while (iNC<nNC && NewCpl_GetGid(arrayNC[iNC])<gid)
      iNC++;

    /* generate additional XIDelCpl-items for all valid NewCpl-items */
    while (iNC<nNC && NewCpl_GetGid(arrayNC[iNC])==gid)
    {
      XIDelCpl *xc = NewXIDelCpl(SLLNewArgs);
      if (xc==NULL)
        HARD_EXIT;

      xc->to      = NewCpl_GetDest(arrayNC[iNC]);                   /* receiver of XIDelCpl */
      xc->prio    = PRIO_INVALID;
      xc->te.gid  = gid;                                 /* the object's gid     */
      /*
         printf("%4d: DelCpl 3      %08x %d %d\n",me,gid,xc->to,xc->prio);
       */

      iNC++;
    }
  }
}





/****************************************************************************/


/*
        this function is called by DDD_HdrDestructor!
 */
void ddd_XferRegisterDelete (DDD_HDR hdr)
{
  COUPLING *cpl;
  XIDelObj *xi;

  /* create new XIDelObj */
  xi      = NewXIDelObj(SLLNewArgs);
  if (xi==NULL)
    HARD_EXIT;

  xi->gid = OBJ_GID(hdr);
  xi->delcpls = NULL;

  /*
          now generate XIDelCpl-items, one for each existing coupling.
          these items serve as notification of this delete operation
          for remote processors with same object.
          these items are also a intermediate storage for the object's
          coupling list, in case the object is received after deletion
          and the coupling list must be restored.
   */
  for(cpl=ObjCplList(hdr); cpl!=NULL; cpl=CPL_NEXT(cpl))
  {
    XIDelCpl *xc = NewXIDelCpl(SLLNewArgs);
    if (xc==NULL)
      HARD_EXIT;

    xc->to      = CPL_PROC(cpl);                 /* receiver of XIDelCpl */
    xc->prio    = cpl->prio;                     /* remember priority    */
    xc->te.gid  = OBJ_GID(hdr);                  /* the object's gid     */

    /* we must remember couplings for eventual restoring
       (if this object is received from another proc) */
    xc->next = xi->delcpls;
    xi->delcpls = xc;
  }
}



/****************************************************************************/

/*
        management functions for XferMode.

        these functions control the mode the xfer-module is
        currently in. this is used for error detection, but
        also for correct detection of coupling inconsistencies
        and recovery.
 */

char *XferModeName (int mode)
{
  switch(mode)
  {
  case XMODE_IDLE : return "idle-mode";
  case XMODE_CMDS : return "commands-mode";
  case XMODE_BUSY : return "busy-mode";
  }
  return "unknown-mode";
}


static void XferSetMode (int mode)
{
  xferGlobals.xferMode = mode;

#       if DebugXfer<=8
  sprintf(cBuffer, "%4d: XferMode=%s.\n",
          me, XferModeName(xferGlobals.xferMode));
  DDD_PrintDebug(cBuffer);
#       endif
}


static int XferSuccMode (int mode)
{
  switch(mode)
  {
  case XMODE_IDLE : return XMODE_CMDS;
  case XMODE_CMDS : return XMODE_BUSY;
  case XMODE_BUSY : return XMODE_IDLE;
  }
  return XMODE_IDLE;
}



int XferMode (void)
{
  return xferGlobals.xferMode;
}


int ddd_XferActive (void)
{
  return xferGlobals.xferMode!=XMODE_IDLE;
}


int XferStepMode (int old)
{
  if (xferGlobals.xferMode!=old)
  {
    sprintf(cBuffer, "wrong xfer-mode (currently in %s, expected %s)",
            XferModeName(xferGlobals.xferMode), XferModeName(old));
    DDD_PrintError('E', 6200, cBuffer);
    return FALSE;
  }

  XferSetMode(XferSuccMode(xferGlobals.xferMode));
  return TRUE;
}


/****************************************************************************/


void ddd_XferInit (void)
{
  /* switch off heap usage, will be switched on during XferBegin/End */
  xferGlobals.useHeap = FALSE;

  /* set kind of TMEM alloc/free requests */
  xfer_SetTmpMem(TMEM_ANY);

  /* init control structures for XferInfo-items in first (?) message */
  xferGlobals.setXICopyObj = New_XICopyObjSet();
  xferGlobals.setXISetPrio = New_XISetPrioSet();
  InitXIDelCmd();
  InitXIDelObj();
  InitXINewCpl();
  InitXIOldCpl();

  /* init control structures for XferInfo-items for second (?) message */
  InitXIDelCpl();
  InitXIModCpl();
  InitXIAddCpl();


  XferSetMode(XMODE_IDLE);

  xferGlobals.objmsg_t = LC_NewMsgType("XferMsg");
  xferGlobals.symtab_id = LC_NewMsgTable("SymTab",
                                         xferGlobals.objmsg_t, sizeof(SYMTAB_ENTRY));
  xferGlobals.objtab_id = LC_NewMsgTable("ObjTab",
                                         xferGlobals.objmsg_t, sizeof(OBJTAB_ENTRY));
  xferGlobals.newcpl_id = LC_NewMsgTable("NewCpl",
                                         xferGlobals.objmsg_t, sizeof(TENewCpl));
  xferGlobals.oldcpl_id = LC_NewMsgTable("OldCpl",
                                         xferGlobals.objmsg_t, sizeof(TEOldCpl));
  xferGlobals.objmem_id = LC_NewMsgChunk("ObjMem",
                                         xferGlobals.objmsg_t);

  /* not used anymore
          xferGlobals.deltab_id =
                  LC_NewMsgTable(xferGlobals.objmsg_t, sizeof(DELTAB_ENTRY));
          xferGlobals.priotab_id =
                  LC_NewMsgTable(xferGlobals.objmsg_t, sizeof(CPLTAB_ENTRY));
   */

  CplMsgInit();
  CmdMsgInit();
}


void ddd_XferExit (void)
{
  /* set kind of TMEM alloc/free requests */
  xfer_SetTmpMem(TMEM_ANY);

  CmdMsgExit();
  CplMsgExit();

  /* TODO data (e.g., lists&trees of XI-items) should be freed!! */
}



/****************************************************************************/

#undef _FADR


/****************************************************************************/
