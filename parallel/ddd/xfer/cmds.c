// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      cmds.c                                                        */
/*                                                                          */
/* Purpose:   DDD-commands for Transfer Module                              */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   931130 kb  begin (xfer.c)                                     */
/*            950321 kb  added variable sized objects (XferCopyObjX)        */
/*            950405 kb  V1.3: extracted from xfer.c                        */
/*            960703 kb  splitted XferInfo-list into ObjXfer and CplXfer    */
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


/*
   #define DebugXfer 2
 */


/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/



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


XICopyObj *theXIAddData;



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


static int sort_XICopyObj (const void *e1, const void *e2)
{
  REGISTER XICopyObj *item1 = *((XICopyObj **)e1);
  REGISTER XICopyObj *item2 = *((XICopyObj **)e2);

  if (item1->dest < item2->dest) return(-1);
  if (item1->dest > item2->dest) return(1);

  if (item1->gid < item2->gid) return(-1);
  if (item1->gid > item2->gid) return(1);

  /* sorting according to priority is not necessary anymore,
     equal items with different priorities will be sorted
     out according to PriorityMerge(). KB 970129
     if (item1->prio < item2->prio) return(-1);
     if (item1->prio > item2->prio) return(1);
   */

  return(0);
}



static int sort_XISetPrio (const void *e1, const void *e2)
{
  REGISTER XISetPrio *item1 = *((XISetPrio **)e1);
  REGISTER XISetPrio *item2 = *((XISetPrio **)e2);

  /* ascending GID is needed for ExecLocalXI___ */
  if (item1->gid < item2->gid) return(-1);
  if (item1->gid > item2->gid) return(1);

  /* sorting according to priority is not necessary anymore,
     equal items with different priorities will be sorted
     out according to PriorityMerge(). KB 970129
     if (item1->prio < item2->prio) return(-1);
     if (item1->prio > item2->prio) return(1);
   */

  return(0);
}



static int sort_XIDelCmd (const void *e1, const void *e2)
{
  REGISTER XIDelCmd *item1 = *((XIDelCmd **)e1);
  REGISTER XIDelCmd *item2 = *((XIDelCmd **)e2);

  /* ascending GID is needed for ExecLocalXIDelCmds */
  if (OBJ_GID(item1->hdr) < OBJ_GID(item2->hdr)) return(-1);
  if (OBJ_GID(item1->hdr) > OBJ_GID(item2->hdr)) return(1);

  return(0);
}


static int sort_XIDelObj (const void *e1, const void *e2)
{
  REGISTER XIDelObj *item1 = *((XIDelObj **)e1);
  REGISTER XIDelObj *item2 = *((XIDelObj **)e2);

  /* ascending GID is needed for ExecLocalXIDelObjs */
  if (item1->gid < item2->gid) return(-1);
  if (item1->gid > item2->gid) return(1);

  return(0);
}


static int sort_XINewCpl (const void *e1, const void *e2)
{
  REGISTER XINewCpl *item1 = *((XINewCpl **)e1);
  REGISTER XINewCpl *item2 = *((XINewCpl **)e2);

  /* receiving processor */
  if (item1->to < item2->to) return(-1);
  if (item1->to > item2->to) return(1);

  return(0);
}



static int sort_XIOldCpl (const void *e1, const void *e2)
{
  REGISTER XIOldCpl *item1 = *((XIOldCpl **)e1);
  REGISTER XIOldCpl *item2 = *((XIOldCpl **)e2);
  REGISTER DDD_GID gid1, gid2;

  /* receiving processor */
  if (item1->to < item2->to) return(-1);
  if (item1->to > item2->to) return(1);

  /* ascending GID is needed for UnpackOldCplTab on receiver side */
  gid1 = item1->te.gid;
  gid2 = item2->te.gid;
  if (gid1 < gid2) return(-1);
  if (gid1 > gid2) return(1);

  return(0);
}



static int sort_XIDelCpl (const void *e1, const void *e2)
{
  REGISTER XIDelCpl *item1 = *((XIDelCpl **)e1);
  REGISTER XIDelCpl *item2 = *((XIDelCpl **)e2);
  REGISTER DDD_GID gid1, gid2;

  /* receiving processor */
  if (item1->to < item2->to) return(-1);
  if (item1->to > item2->to) return(1);

  /* ascending GID is needed by CplMsgUnpack on receiver side */
  gid1 = item1->te.gid;
  gid2 = item2->te.gid;
  if (gid1 < gid2) return(-1);
  if (gid1 > gid2) return(1);

  return(0);
}


static int sort_XIModCpl (const void *e1, const void *e2)
{
  REGISTER XIModCpl *item1 = *((XIModCpl **)e1);
  REGISTER XIModCpl *item2 = *((XIModCpl **)e2);
  REGISTER DDD_GID gid1, gid2;

  /* receiving processor */
  if (item1->to < item2->to) return(-1);
  if (item1->to > item2->to) return(1);

  /* ascending GID is needed by CplMsgUnpack on receiver side */
  gid1 = item1->te.gid;
  gid2 = item2->te.gid;
  if (gid1 < gid2) return(-1);
  if (gid1 > gid2) return(1);

  /* sorting according to priority is not necessary anymore,
     equal items with different priorities will be sorted
     out according to PriorityMerge(). KB 970129
     if (item1->te.prio < item2->te.prio) return(-1);
     if (item1->te.prio > item2->te.prio) return(1);
   */

  return(0);
}


static int sort_XIAddCpl (const void *e1, const void *e2)
{
  REGISTER XIAddCpl *item1 = *((XIAddCpl **)e1);
  REGISTER XIAddCpl *item2 = *((XIAddCpl **)e2);
  REGISTER DDD_GID gid1, gid2;

  /* receiving processor */
  if (item1->to < item2->to) return(-1);
  if (item1->to > item2->to) return(1);

  /* ascending GID is needed by CplMsgUnpack on receiver side */
  gid1 = item1->te.gid;
  gid2 = item2->te.gid;
  if (gid1 < gid2) return(-1);
  if (gid1 > gid2) return(1);

  return(0);
}



/****************************************************************************/



/*
        eliminate double XICopyObj-items.
        merge priorities from similar XICopyObj-items.

        the items have been sorted according to key (dest,gid),
        all in ascending order. if dest or gid are different, then
        at least the first item is relevant. if both are equal,
        we merge priorities and get a new priority together with
        the information whether first item wins over second.
        if first item wins, it is switched into second position and
        the second item (now on first position) is rejected.
        if second item wins, first item is rejected.
        in both cases, we use the new priority for next comparison.

        this implements rule XFER-C1.
 */
static int unify_XICopyObj (XICopyObj **i1p, XICopyObj **i2p)
{
  XICopyObj *i1 = *i1p, *i2 = *i2p;
  DDD_PRIO newprio;
  int ret;

  /* if items are different in gid or dest, take first item */
  if ((i1->hdr != i2->hdr) || (i1->dest != i2->dest))
    return TRUE;

  /* items have equal gid and dest, we must check priority */
  ret = PriorityMerge(&theTypeDefs[OBJ_TYPE(i1->hdr)],
                      i1->prio, i2->prio, &newprio);

  if (ret==PRIO_FIRST || ret==PRIO_UNKNOWN)
  {
    /* i1 is winner, take it, switch it into second position,
       signal rejection of i2 (now on first position).
       use new priority */

    i1->prio = newprio;
    *i1p = i2; *i2p = i1;              /* switch pointers */
  }
  else
  {
    /* i1 lost, i2 is winner. throw away i1, but
       use new priority for next comparison */
    i2->prio = newprio;
  }

  return FALSE;
}



/*
        eliminate double XISetPrio-items.
        merge priorities from similar XISetPrio-items.

        the items have been sorted according to key (gid), in
        ascending order. if gid (i.e., hdr) is different,
        then at least the first item is relevant. if gid is equal,
        we merge priorities and get a new priority together with
        the information whether first item wins over second.
        if first item wins, it is switched into second position and
        the second item (now on first position) is rejected.
        if second item wins, first item is rejected.
        in both cases, we use the new priority for next comparison.

        this implements rule XFER-P1.
 */
static int unify_XISetPrio (XISetPrio **i1p, XISetPrio **i2p)
{
  XISetPrio *i1 = *i1p, *i2 = *i2p;
  DDD_PRIO newprio;
  int ret;

  /* if items are different in gid or dest, take first item */
  if (i1->hdr != i2->hdr)
    return TRUE;

  /* items have equal gid, we must check priority */
  ret = PriorityMerge(&theTypeDefs[OBJ_TYPE(i1->hdr)],
                      i1->prio, i2->prio, &newprio);

  if (ret==PRIO_FIRST || ret==PRIO_UNKNOWN)
  {
    /* i1 is winner, take it, switch it into second position,
       signal rejection of i2 (now on first position).
       use new priority */

    i1->prio = newprio;
    *i1p = i2; *i2p = i1;              /* switch pointers */
  }
  else
  {
    /* i1 lost, i2 is winner. throw away i1, but
       use new priority for next comparison */
    i2->prio = newprio;
  }

  return FALSE;
}



/*
        eliminate double XIDelCmd-items.

        the items have been sorted according to key
    (gid), all in ascending order.
        if gid (i.e., hdr) is equal,
        the item is skipped.
        this implements rule XFER-D1.

        the number of valid items is returned.
 */
static int unify_XIDelCmd (XIDelCmd **i1, XIDelCmd **i2)
{
  return ((*i1)->hdr != (*i2)->hdr);
}



/*
        eliminate double XIModCpl-items.
        merge priorities from similar XIModCpl-items.

        the items have been sorted according to key (to,gid),
        all in ascending order. if to or gid are different, then
        at least the first item is relevant. if both are equal,
        we merge priorities and get a new priority together with
        the information whether first item wins over second.
        if first item wins, it is switched into second position and
        the second item (now on first position) is rejected.
        if second item wins, first item is rejected.
        in both cases, we use the new priority for next comparison.
 */
static int unify_XIModCpl (XIModCpl **i1p, XIModCpl **i2p)
{
  XIModCpl *i1 = *i1p, *i2 = *i2p;
  DDD_PRIO newprio;
  int ret;

  /* if items are different in gid or dest, take first item */
  if ((i1->to != i2->to) || (i1->te.gid != i2->te.gid))
    return TRUE;

  /* items have equal to and gid, we must check priority */
  ret = PriorityMerge(&theTypeDefs[i1->typ],
                      i1->te.prio, i2->te.prio, &newprio);

  if (ret==PRIO_FIRST || ret==PRIO_UNKNOWN)
  {
    /* i1 is winner, take it, switch it into second position,
       signal rejection of i2 (now on first position).
       use new priority */

    i1->te.prio = newprio;
    *i1p = i2; *i2p = i1;              /* switch pointers */
  }
  else
  {
    /* i1 lost, i2 is winner. throw away i1, but
       use new priority for next comparison */
    i2->te.prio = newprio;
  }

  return FALSE;
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferEnd                                                   */
/*                                                                          */
/* Purpose:   end xfer command phase, start real transfer phase             */
/*            (i.e. transaction commit)                                     */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/


#ifdef C_FRONTEND
void DDD_XferEnd (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::XferEnd (void)
#endif
#ifdef F_FRONTEND
void DDD_XferEnd (void)
#endif
{
  XICopyObj   **arrayXICopyObj, **arrayNewOwners;
  int remXICopyObj, nNewOwners;
  XIDelCmd    **arrayXIDelCmd;
  int remXIDelCmd;
  XIDelObj    **arrayXIDelObj;
  XISetPrio   **arrayXISetPrio;
  int remXISetPrio;
  XINewCpl    **arrayXINewCpl;
  XIOldCpl    **arrayXIOldCpl;
  XIDelCpl    **arrayXIDelCpl;
  int remXIDelCpl;
  XIModCpl    **arrayXIModCpl;
  int remXIModCpl;
  XIAddCpl    **arrayXIAddCpl;
  int obsolete, nRecvMsgs, nSendMsgs;
  XFERMSG     *sendMsgs, *sm=0;
  LC_MSGHANDLE *recvMsgs;
  DDD_HDR     *localCplObjs=NULL;


  /* step mode and check whether call to XferEnd is valid */
  if (!XferStepMode(XMODE_CMDS))
  {
    DDD_PrintError('E', 6011, "DDD_XferEnd() aborted");
    HARD_EXIT;
  }


  /*
          COMMUNICATION PHASE 1
   */

  /* create sorted array of XICopyObj-items, and unify it */
  arrayXICopyObj = SortedArrayXICopyObj(sort_XICopyObj);
  remXICopyObj   = UnifyXICopyObj(arrayXICopyObj, unify_XICopyObj);
  obsolete = nXICopyObj-remXICopyObj;

  /* send Cpl-info about new objects to owners of other local copies */
  arrayNewOwners = CplClosureEstimate(
    arrayXICopyObj, remXICopyObj,
    &nNewOwners);

  /* create sorted array of XINewCpl- and XIOldCpl-items.
     TODO. if efficiency is a problem here, use b-tree or similar
           data structure to improve performance. */
  arrayXINewCpl = SortedArrayXINewCpl(sort_XINewCpl);
  arrayXIOldCpl = SortedArrayXIOldCpl(sort_XIOldCpl);


  /* prepare msgs for objects and XINewCpl-items */
  nSendMsgs = PrepareObjMsgs(
    arrayXICopyObj, remXICopyObj,
    arrayXINewCpl, nXINewCpl,
    arrayXIOldCpl, nXIOldCpl,
    &sendMsgs);

  /* init communication topology */
  nRecvMsgs = LC_Connect(xferGlobals.objmsg_t);

  /* build obj msgs on sender side and start send */
  XferPackMsgs(sendMsgs);


  /*
          now messages are in the net, use spare time
   */

  /* create sorted array of XISetPrio-items, and unify it */
  arrayXISetPrio = SortedArrayXISetPrio(sort_XISetPrio);
  remXISetPrio   = UnifyXISetPrio(arrayXISetPrio, unify_XISetPrio);
  obsolete += (nXISetPrio-remXISetPrio);

  /* create sorted array of XIDelCmd-items, and unify it */
  arrayXIDelCmd = SortedArrayXIDelCmd(sort_XIDelCmd);
  remXIDelCmd   = UnifyXIDelCmd(arrayXIDelCmd, unify_XIDelCmd);
  obsolete += (nXIDelCmd-remXIDelCmd);

  /* execute local commands */
  /* NOTE: messages have been build before in order to allow
           deletion of objects. */
  ExecLocalXIDelCmd(arrayXIDelCmd,  remXIDelCmd);

  /* now all XIDelObj-items have been created. these come from:
          1. application->DDD_XferDeleteObj->XIDelCmd->HdrDestructor->
             ->XferRegisterDelete
          2. HANDLER_DELETE->HdrDestructor (for dependent object)->
             ->XferRegisterDelete
   */

  /* create sorted array of XIDelObj-items */
  arrayXIDelObj = SortedArrayXIDelObj(sort_XIDelObj);

  ExecLocalXISetPrio(arrayXISetPrio, remXISetPrio,
                     arrayXIDelObj,  nXIDelObj,
                     arrayNewOwners, nNewOwners);
  ExecLocalXIDelObj(arrayXIDelObj,  nXIDelObj,
                    arrayNewOwners, nNewOwners);


  /* get sorted list of local objects with couplings */
  localCplObjs = LocalCoupledObjectsList();


  if (obsolete>0) {
    if (DDD_GetOption(OPT_INFO_XFER)==OPT_ON)
    {
      sprintf(cBuffer, "DDD MESG [%03d]: %4d from %4d xfer-cmds obsolete.\n",
              me, obsolete, nXICopyObj+nXISetPrio+nXIDelObj);
      DDD_PrintLine(cBuffer);
    }
  }

  /*
          nothing more to do until incoming messages arrive
   */


  /* wait for communication-completion (send AND receive) */
  recvMsgs = LC_Communicate();
  XferUnpack(recvMsgs, nRecvMsgs,
             localCplObjs, nCpls,
             arrayXISetPrio, remXISetPrio,
             arrayXIDelObj, nXIDelObj,
             arrayXICopyObj, remXICopyObj,
             arrayNewOwners, nNewOwners);
  LC_Cleanup();

  /* recreate sorted list of local coupled objects,
     old list might be corrupt due to creation of new objects */
  if (localCplObjs!=NULL) FreeTmp(localCplObjs);
  localCplObjs = LocalCoupledObjectsList();


  /* create sorted array of XIDelCpl-, XIModCpl- and XIAddCpl-items.
     TODO. if efficiency is a problem here, use b-tree or similar
           data structure to improve performance. */
  arrayXIDelCpl = SortedArrayXIDelCpl(sort_XIDelCpl);
  arrayXIModCpl = SortedArrayXIModCpl(sort_XIModCpl);
  arrayXIAddCpl = SortedArrayXIAddCpl(sort_XIAddCpl);


  /* some XIDelCpls have been invalidated by UpdateCoupling(),
     decrease list size to avoid sending them */
  remXIDelCpl = nXIDelCpl;
  while (remXIDelCpl>0 && arrayXIDelCpl[remXIDelCpl-1]->to == procs)
    remXIDelCpl--;

  remXIModCpl   = UnifyXIModCpl(arrayXIModCpl, unify_XIModCpl);

  /*
     printf("%4d: %d XIDelCpls obsolete\n", me, nXIDelCpl-remXIDelCpl);
   */


  /*
          COMMUNICATION PHASE 2
   */

  CommunicateCplMsgs(arrayXIDelCpl, remXIDelCpl,
                     arrayXIModCpl, remXIModCpl,
                     arrayXIAddCpl, nXIAddCpl,
                     localCplObjs, nCpls);


  /*
          free temporary storage
   */

  if (arrayXICopyObj!=NULL) FreeTmp(arrayXICopyObj);
  if (arrayNewOwners!=NULL) FreeTmp(arrayNewOwners);
  FreeAllXIAddData();
  FreeAllXICopyObj();

  if (arrayXISetPrio!=NULL) FreeTmp(arrayXISetPrio);
  FreeAllXISetPrio();

  if (arrayXIDelCmd!=NULL) FreeTmp(arrayXIDelCmd);
  FreeAllXIDelCmd();

  if (arrayXIDelObj!=NULL) FreeTmp(arrayXIDelObj);
  FreeAllXIDelObj();

  if (arrayXINewCpl!=NULL) FreeTmp(arrayXINewCpl);
  FreeAllXINewCpl();

  if (arrayXIOldCpl!=NULL) FreeTmp(arrayXIOldCpl);
  FreeAllXIOldCpl();

  if (arrayXIDelCpl!=NULL) FreeTmp(arrayXIDelCpl);
  FreeAllXIDelCpl();

  if (arrayXIModCpl!=NULL) FreeTmp(arrayXIModCpl);
  FreeAllXIModCpl();

  if (arrayXIAddCpl!=NULL) FreeTmp(arrayXIAddCpl);
  FreeAllXIAddCpl();

  if (localCplObjs!=NULL) FreeTmp(localCplObjs);

  for(; sendMsgs!=NULL; sendMsgs=sm)
  {
    sm = sendMsgs->next;
    FreeTmp(sendMsgs);
  }


#       if DebugXfer<=4
  sprintf(cBuffer,"%4d: XferEnd, before IFAllFromScratch().\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  /* re-create all interfaces and step XMODE */
  IFAllFromScratch();
  XferStepMode(XMODE_BUSY);
}




/* ablage fuer debug-ausgabe

 #	if DebugXfer<=4
                sprintf(cBuffer,"%4d: XferEnd, after XferDeleteObjects(): "
                        "send=%d recv=%d\n", me, nSendMsgs, nRecvMsgs);
                DDD_PrintDebug(cBuffer);
 #	endif
 */


/* aufbewahrungsort fuer timing-direktiven

        STAT_ZEROTIMER;
        STAT_RESET;

        STAT_TIMER(10); STAT_RESET;
        STAT_TIMER(12);
 */


/****************************************************************************/
/*                                                                          */
/* Function:  XferCopyObj                                                   */
/*                                                                          */
/* Purpose:   xfer command: create copy of obj on proc with priority prio   */
/*                                                                          */
/* Input:     hdr:  DDD-header of object to be copied                       */
/*            proc: receiver of object                                      */
/*            prio: local priority of copy                                  */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/



void DDD_XferPrioChange (DDD_HDR hdr, DDD_PRIO prio)
{
  XISetPrio *xi = NewXISetPrio();
  xi->hdr  = hdr;
  xi->gid  = OBJ_GID(hdr);
  xi->prio = prio;

#       if DebugXfer<=2
  sprintf(cBuffer, "%4d: DDD_XferPrioChange %08x, prio=%d\n",
          me, OBJ_GID(hdr), prio);
  DDD_PrintDebug(cBuffer);
#       endif
}



static void XferInitCopyInfo (DDD_HDR hdr,
                              TYPE_DESC *desc,
                              size_t size,
                              DDD_PROC dest,
                              DDD_PRIO prio)
{
  if (!XferActive())
  {
    DDD_PrintError('E', 6012, "Missing DDD_XferBegin(). aborted");
    HARD_EXIT;
  }

  if (dest>=procs)
  {
    sprintf(cBuffer, "cannot transfer %08x to processor %d (procs=%d)",
            OBJ_GID(hdr), dest, procs);
    DDD_PrintError('E', 6003, cBuffer);
    HARD_EXIT;
  }

  if (prio>=MAX_PRIO)
  {
    sprintf(cBuffer, "priority must be less than %d (prio=%d) in xfer-cmd",
            MAX_PRIO, prio);
    DDD_PrintError('E', 6004, cBuffer);
    HARD_EXIT;
  }

  if (dest==me)
  {
    /* XFER-C4: XferCopyObj degrades to SetPrio command */
    XISetPrio *xi = NewXISetPrio();
    xi->hdr  = hdr;
    xi->gid  = OBJ_GID(hdr);
    xi->prio = prio;

    /* nevertheless, call XFERCOPY-handler! */

    /* reset for eventual AddData-calls during handler execution */
    theXIAddData = NULL;

    /* call application handler for xfer of dependent objects */
    if (desc->handlerXFERCOPY)
    {
      DDD_OBJ obj = HDR2OBJ(hdr,desc);
      desc->handlerXFERCOPY(_FADR obj, _FADR dest, _FADR prio);
    }

    /* theXIAddData might be changed during handler execution */
    theXIAddData = NULL;
  }
  else
  {
    /* this is a real transfer to remote proc */
    XICopyObj  *xi = NewXICopyObj();
    xi->hdr  = hdr;
    xi->gid  = OBJ_GID(hdr);
    xi->size = size;
    xi->dest = dest;
    xi->prio = prio;

    xi->add    = NULL;
    xi->addLen = 0;

    /* set XferAddInfo for evtl AddData-calls during handler execution */
    theXIAddData = xi;

    /* call application handler for xfer of dependent objects */
    if (desc->handlerXFERCOPY)
    {
      DDD_OBJ obj = HDR2OBJ(hdr,desc);
      desc->handlerXFERCOPY(_FADR obj, _FADR dest, _FADR prio);
    }

    /* theXIAddData might be changed during handler execution */
    theXIAddData = xi;
  }
}






#ifdef C_FRONTEND
void DDD_XferCopyObj (DDD_HDR hdr, DDD_PROC proc, DDD_PRIO prio)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Object::XferCopyObj (DDD_PROC proc, DDD_PRIO prio)
{
  DDD_HDR hdr = &_hdr;
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
TYPE_DESC *desc =  &(theTypeDefs[OBJ_TYPE(hdr)]);

#       if DebugXfer<=2
sprintf(cBuffer, "%4d: DDD_XferCopyObj %08x, proc=%d prio=%d\n",
        me, OBJ_GID(hdr), proc, prio);
DDD_PrintDebug(cBuffer);
#       endif

XferInitCopyInfo(hdr, desc, desc->size, proc, prio);
}
#endif

#ifdef F_FRONTEND

void DDD_XferCopyObj (DDD_TYPE *type, DDD_OBJ *obj, DDD_PROC *proc,
                      DDD_PRIO *prio)
{
  TYPE_DESC *desc = &(theTypeDefs[*type]);
  DDD_HDR hdr = OBJ2HDR(*obj,desc);

#       if DebugXfer<=2
  sprintf(cBuffer, "%4d: DDD_XferCopyObj %08x(%d), proc=%d prio=%d\n",
          me, OBJ_GID(hdr), *obj, *proc, *prio);
  DDD_PrintDebug(cBuffer);
#       endif

  XferInitCopyInfo(hdr, desc, desc->size, *proc, *prio);
}
#endif



/* XferCopyObj for variable sized objects, 950321 KB */

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)

void DDD_XferCopyObjX (DDD_HDR hdr, DDD_PROC proc, DDD_PRIO prio, size_t size)
{
  TYPE_DESC *desc =  &(theTypeDefs[OBJ_TYPE(hdr)]);

#       if DebugXfer<=2
  sprintf(cBuffer, "%4d: DDD_XferCopyObjX %08x, proc=%d prio=%d size=%d\n",
          me, OBJ_GID(hdr), proc, prio, size);
  DDD_PrintDebug(cBuffer);
#       endif

  if ((desc->size!=size) && (DDD_GetOption(OPT_WARNING_VARSIZE_OBJ)==OPT_ON))
  {
    DDD_PrintError('W', 6001,
                   "object size differs from declared size in DDD_XferCopyObjX");
  }

  if ((desc->size>size) && (DDD_GetOption(OPT_WARNING_SMALLSIZE)==OPT_ON))
  {
    DDD_PrintError('W', 6002,
                   "object size smaller than declared size in DDD_XferCopyObjX");
  }

  XferInitCopyInfo(hdr, desc, size, proc, prio);
}

/* No Fortran Interface. Fortran won't support variable sized object */

#endif


#ifdef C_FRONTEND
void DDD_XferAddData (int cnt, DDD_TYPE typ)
{
  XFERADDDATA *xa;
  TYPE_DESC   *descDepTyp;

#       if DebugXfer<=2
  sprintf(cBuffer, "%4d: DDD_XferAddData cnt=%d typ=%d\n", me, cnt, typ);
  DDD_PrintDebug(cBuffer);
#       endif

  if (theXIAddData==NULL) return;

  xa = NewXIAddData();
  xa->addCnt = cnt;
  xa->addTyp = typ;
  xa->sizes  = NULL;

  if (typ<DDD_USER_DATA || typ>DDD_USER_DATA_MAX)
  {
    /* normal dependent object */
    descDepTyp =  &(theTypeDefs[typ]);

    xa->addLen       = CEIL(descDepTyp->size)   * cnt;
    xa->addNPointers = (descDepTyp->nPointers) * cnt;
  }
  else
  {
    /* stream of bytes, since V1.2 */
    /* many streams, since V1.7.8 */
    xa->addLen       = cnt;
    xa->addNPointers = 0;
  }

  theXIAddData->addLen += xa->addLen;
}



void DDD_XferAddDataX (int cnt, DDD_TYPE typ, size_t *sizes)
{
  XFERADDDATA *xa;
  TYPE_DESC   *descDepTyp;
  int i;

#       if DebugXfer<=2
  sprintf(cBuffer,"%4d: DDD_XferAddData cnt=%d typ=%d\n", me, cnt, typ);
  DDD_PrintDebug(cBuffer);
#       endif

  if (theXIAddData==NULL) return;

  xa = NewXIAddData();
  xa->addCnt = cnt;
  xa->addTyp = typ;

  if (typ<DDD_USER_DATA || typ>DDD_USER_DATA_MAX)
  {
    /* copy sizes array */
    xa->sizes = AddDataAllocSizes(cnt);
    memcpy(xa->sizes, sizes, sizeof(int)*cnt);

    /* normal dependent object */
    descDepTyp =  &(theTypeDefs[typ]);

    xa->addLen = 0;
    for (i=0; i<cnt; i++)
    {
      xa->addLen += CEIL(sizes[i]);
    }
    xa->addNPointers = (descDepTyp->nPointers) * cnt;
  }
  else
  {
    /* stream of bytes, since V1.2 */
    /* many streams, since V1.7.8 */
    xa->addLen       = cnt;
    xa->addNPointers = 0;
  }

  theXIAddData->addLen += xa->addLen;
}
#endif



/****************************************************************************/
/*                                                                          */
/* Function:  XferDeleteObj                                                 */
/*                                                                          */
/* Purpose:   xfer command: delete local object                             */
/*                                                                          */
/* Input:     object to be deleted                                          */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_XferDeleteObj (DDD_HDR hdr)
#endif
#ifdef CPP_FRONTEND
void DDD_Object::XferDeleteObj (void)
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
{
        #ifdef CPP_FRONTEND
  DDD_HDR hdr = &_hdr;
        #endif
  TYPE_DESC *desc =  &(theTypeDefs[OBJ_TYPE(hdr)]);
  XIDelCmd  *dc = NewXIDelCmd();


        #ifdef CPP_FRONTEND
  dc->obj = this;
        #endif

  dc->hdr = hdr;

#       if DebugXfer<=2
  sprintf(cBuffer,"%4d: DDD_XferDeleteObj %08x\n",
          me, OBJ_GID(hdr));
  DDD_PrintDebug(cBuffer);
#       endif


  /* call application handler for deletion of dependent objects */
  if (desc->handlerXFERDELETE!=NULL)
    desc->handlerXFERDELETE(HDR2OBJ(hdr,desc));
}
#endif


#ifdef F_FRONTEND
void DDD_XferDeleteObj (DDD_TYPE *type, DDD_OBJ *obj)
{
  DDD_HDR hdr = OBJ2HDR(*obj,&theTypeDefs[*type]);
  XIDelCmd  *dc = NewXIDelCmd();

  dc->hdr = hdr;


#       if DebugXfer<=2
  sprintf(cBuffer,"%4d: DDD_XferDeleteObj %08x(%d)\n",
          me, OBJ_GID(hdr), *obj);
  DDD_PrintDebug(cBuffer);
#       endif


  /* call application handler for deletion of dependent objects */
  /*
          if (desc->handler[HANDLER_XFERDELETE]!=NULL)
                  desc->handler[HANDLER_XFERDELETE](???);
   */
}
#endif



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferBegin                                                 */
/*                                                                          */
/* Purpose:   start xfer command phase                                      */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_XferBegin (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::XferBegin (void)
#endif
#ifdef F_FRONTEND
void DDD_XferBegin (void)
#endif
{
  theXIAddData = NULL;


  /* step mode and check whether call to XferBegin is valid */
  if (!XferStepMode(XMODE_IDLE))
  {
    DDD_PrintError('E', 6010, "DDD_XferBegin() aborted");
    HARD_EXIT;
  }
}


/****************************************************************************/

#undef _FADR

/****************************************************************************/
