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



/* TODO remove this */
void GetSizesXIAddData (int *, int *, size_t *, size_t *);


/*
        compute and display memory resources used
 */
static void DisplayMemResources (void)
{
  int nSegms=0, nItems=0, nNodes=0;
  size_t memAllocated=0, memUsed=0;

  GetSizesXIAddData(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XIAddData segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);


  XICopyObjSet_GetResources(xferGlobals.setXICopyObj,
                            &nSegms, &nItems, &nNodes, &memAllocated, &memUsed);
  if (nSegms>0) {
    printf("%4d: XferEnd, XICopyObj "
           "segms=%d items=%d nodes=%ld allocated=%ld used=%ld\n",
           me, nSegms, nItems, nNodes, (long)memAllocated, (long)memUsed);
  }


        #ifdef XICOPYOBJ_DETAILED_RESOURCES
  /* this is a different version, split up into BTree and SegmList */
  XICopyObjSegmList_GetResources(xferGlobals.setXICopyObj->list,
                                 &nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XICopyObj segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);

  XICopyObjBTree_GetResources(xferGlobals.setXICopyObj->tree,
                              &nNodes, &nItems, &memAllocated, &memUsed);
  if (nItems>0)
    printf("%4d: XferEnd, XICopyObj nodes=%d items=%d allocated=%ld used=%ld\n",
           me, nNodes, nItems, (long)memAllocated, (long)memUsed);
        #endif


  XISetPrioSet_GetResources(xferGlobals.setXISetPrio,
                            &nSegms, &nItems, &nNodes, &memAllocated, &memUsed);
  if (nSegms>0) {
    printf("%4d: XferEnd, XISetPrio "
           "segms=%d items=%d nodes=%ld allocated=%ld used=%ld\n",
           me, nSegms, nItems, nNodes, (long)memAllocated, (long)memUsed);
  }


  GetSizesXIDelCmd(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XIDelCmd  segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);

  GetSizesXIDelObj(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XIDelObj  segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);

  GetSizesXINewCpl(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XINewCpl  segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);

  GetSizesXIOldCpl(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XIOldCpl  segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);

  GetSizesXIDelCpl(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XIDelCpl  segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);

  GetSizesXIModCpl(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XIModCpl  segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);

  GetSizesXIAddCpl(&nSegms, &nItems, &memAllocated, &memUsed);
  if (nSegms>0)
    printf("%4d: XferEnd, XIAddCpl  segms=%d items=%d allocated=%ld used=%ld\n",
           me, nSegms, nItems, (long)memAllocated, (long)memUsed);


  /*
          sprintf(cBuffer, "%4d: XferEnd, segms=%d items=%d allocated=%ld used=%ld\n",
                  me, nSegms, nItems, (long)memAllocated, (long)memUsed);
          DDD_PrintDebug(cBuffer);
   */
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferEnd                                                   */
/*                                                                          */
/****************************************************************************/

/**
        End of transfer phase.
        This function starts the object transfer process. After a call to
        this function (on all processors) all {\bf Transfer}-commands since
        the last call to \funk{XferBegin} are executed. This involves
        a set of local communications between the processors.
 */

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
  XICopyObjPtrArray *arrayXICopyObj;
  XICopyObj   **arrayNewOwners;
  int nNewOwners;
  XIDelCmd    **arrayXIDelCmd;
  int remXIDelCmd, prunedXIDelCmd;
  XIDelObj    **arrayXIDelObj;
  XISetPrioPtrArray *arrayXISetPrio;
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
  size_t sendMem=0, recvMem=0;
  int DelCmds_were_pruned;

  STAT_SET_MODULE(DDD_MODULE_XFER);
  STAT_ZEROALL;

  /* step mode and check whether call to XferEnd is valid */
  if (!XferStepMode(XMODE_CMDS))
  {
    DDD_PrintError('E', 6011, "DDD_XferEnd() aborted");
    HARD_EXIT;
  }


  /*
          PREPARATION PHASE
   */
  STAT_RESET;
  /* get sorted array of XICopyObj-items */
  arrayXICopyObj = XICopyObjSet_GetArray(xferGlobals.setXICopyObj);
  obsolete = XICopyObjSet_GetNDiscarded(xferGlobals.setXICopyObj);

  /* debugging output, write all XICopyObjs to file
     if (XICopyObjSet_GetNItems(xferGlobals.setXICopyObj)>0)
     {
          FILE *ff = fopen("xfer.dump","w");
          XICopyObjSet_Print(xferGlobals.setXICopyObj, 2, ff);
          fclose(ff);
     }
   */


  /*
          (OPTIONAL) COMMUNICATION PHASE 0
   */
  if (DDD_GetOption(OPT_XFER_PRUNE_DELETE)==OPT_ON)
  {
    /*
            for each XferDelete-Cmd: if there exists at least
                    one XferCopy-cmd with destination=me, then the
                    XferDelete-Cmd is discarded.

            NOTE: the priorities behave like in the specification,
                    i.e., incoming objects with lower priority than the
                    local (deleted) object won't be rejected.
     */
    /* create sorted array of XIDelCmd-items, and unify it */
    /* in case of pruning set to OPT_OFF, this sorting/unifying
       step is done lateron. */
    arrayXIDelCmd = SortedArrayXIDelCmd(sort_XIDelCmd);
    remXIDelCmd   = UnifyXIDelCmd(arrayXIDelCmd, unify_XIDelCmd);
    obsolete += (nXIDelCmd-remXIDelCmd);

    /* do communication and actual pruning */
    prunedXIDelCmd = PruneXIDelCmd(arrayXIDelCmd, remXIDelCmd, arrayXICopyObj);
    obsolete += prunedXIDelCmd;
    remXIDelCmd -= prunedXIDelCmd;

    DelCmds_were_pruned = TRUE;
  }
  else
  {
    DelCmds_were_pruned = FALSE;
  }

  /*
     if (nXIDelCmd>0||remXIDelCmd>0||prunedXIDelCmd>0)
     printf("%4d: XIDelCmd. all=%d rem=%d pruned=%d\n", me,
     nXIDelCmd, remXIDelCmd, prunedXIDelCmd);
   */


  /*
          COMMUNICATION PHASE 1
   */
  STAT_RESET;
  /* send Cpl-info about new objects to owners of other local copies */
  arrayNewOwners = CplClosureEstimate(arrayXICopyObj, &nNewOwners);

  /* create sorted array of XINewCpl- and XIOldCpl-items.
     TODO. if efficiency is a problem here, use b-tree or similar
           data structure to improve performance. */
  arrayXINewCpl = SortedArrayXINewCpl(sort_XINewCpl);
  arrayXIOldCpl = SortedArrayXIOldCpl(sort_XIOldCpl);


  /* prepare msgs for objects and XINewCpl-items */
  nSendMsgs = PrepareObjMsgs(arrayXICopyObj,
                             arrayXINewCpl, nXINewCpl,
                             arrayXIOldCpl, nXIOldCpl,
                             &sendMsgs, &sendMem);


  /*
     DisplayMemResources();
   */

  /* init communication topology */
  nRecvMsgs = LC_Connect(xferGlobals.objmsg_t);
  STAT_TIMER(T_XFER_PREP_MSGS);

  STAT_RESET;
  /* build obj msgs on sender side and start send */
  XferPackMsgs(sendMsgs);
  STAT_TIMER(T_XFER_PACK_SEND);

  /*
          now messages are in the net, use spare time
   */

  /* create sorted array of XISetPrio-items, and unify it */
  STAT_RESET;
  arrayXISetPrio = XISetPrioSet_GetArray(xferGlobals.setXISetPrio);
  obsolete += XISetPrioSet_GetNDiscarded(xferGlobals.setXISetPrio);


  if (!DelCmds_were_pruned)
  {
    /* create sorted array of XIDelCmd-items, and unify it */
    arrayXIDelCmd = SortedArrayXIDelCmd(sort_XIDelCmd);
    remXIDelCmd   = UnifyXIDelCmd(arrayXIDelCmd, unify_XIDelCmd);
    obsolete += (nXIDelCmd-remXIDelCmd);
  }

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

  ExecLocalXISetPrio(arrayXISetPrio,
                     arrayXIDelObj,  nXIDelObj,
                     arrayNewOwners, nNewOwners);
  ExecLocalXIDelObj(arrayXIDelObj,  nXIDelObj,
                    arrayNewOwners, nNewOwners);


  /* get sorted list of local objects with couplings */
  localCplObjs = LocalCoupledObjectsList();

  if (obsolete>0)
  {
    if (DDD_GetOption(OPT_INFO_XFER) & XFER_SHOW_OBSOLETE)
    {
      int all = nXIDelObj+
                XISetPrioSet_GetNItems(xferGlobals.setXISetPrio)+
                XICopyObjSet_GetNItems(xferGlobals.setXICopyObj);

      sprintf(cBuffer, "DDD MESG [%03d]: %4d from %4d xfer-cmds obsolete.\n",
              me, obsolete, all);
      DDD_PrintLine(cBuffer);
    }
  }
  STAT_TIMER(T_XFER_WHILE_COMM);

  /*
          nothing more to do until incoming messages arrive
   */

  /* display information about send-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_XFER) & XFER_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD XFER_SHOW_MSGSALL: ObjMsg.Send\n");
    LC_PrintSendMsgs();
  }


  /* wait for communication-completion (send AND receive) */
  STAT_RESET;
  recvMsgs = LC_Communicate();
  STAT_TIMER(T_XFER_WAIT_RECV);


  /* display information about message buffer sizes */
  if (DDD_GetOption(OPT_INFO_XFER) & XFER_SHOW_MEMUSAGE)
  {
    int k;

    /* sum up sizes of receive mesg buffers */
    for(k=0; k<nRecvMsgs; k++)
    {
      recvMem += LC_GetBufferSize(recvMsgs[k]);
    }

    sprintf(cBuffer,
            "DDD MESG [%03d]: SHOW_MEM "
            "msgs  send=%010ld recv=%010ld all=%010ld\n",
            me, (long)sendMem, (long)recvMem, (long)(sendMem+recvMem));
    DDD_PrintLine(cBuffer);
  }

  /* display information about recv-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_XFER) & XFER_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD XFER_SHOW_MSGSALL: ObjMsg.Recv\n");
    LC_PrintRecvMsgs();
  }

  /* unpack messages */
  STAT_RESET;
  XferUnpack(recvMsgs, nRecvMsgs,
             localCplObjs, NCpl_Get,
             arrayXISetPrio,
             arrayXIDelObj, nXIDelObj,
             arrayXICopyObj,
             arrayNewOwners, nNewOwners);
  LC_Cleanup();
  STAT_TIMER(T_XFER_UNPACK);

  /* recreate sorted list of local coupled objects,
     old list might be corrupt due to creation of new objects */
  STAT_RESET;
  FreeLocalCoupledObjectsList(localCplObjs);
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
  STAT_TIMER(T_XFER_PREP_CPL);

  /*
     printf("%4d: %d XIDelCpls obsolete\n", me, nXIDelCpl-remXIDelCpl);
   */


  /*
          COMMUNICATION PHASE 2
   */

  STAT_RESET;
  CommunicateCplMsgs(arrayXIDelCpl, remXIDelCpl,
                     arrayXIModCpl, remXIModCpl,
                     arrayXIAddCpl, nXIAddCpl,
                     localCplObjs, NCpl_Get);
  STAT_TIMER(T_XFER_CPLMSG);


  /*
          free temporary storage
   */

  XICopyObjPtrArray_Free(arrayXICopyObj);
  XICopyObjSet_Reset(xferGlobals.setXICopyObj);

  if (arrayNewOwners!=NULL) OO_Free (arrayNewOwners,0);
  FreeAllXIAddData();

  XISetPrioPtrArray_Free(arrayXISetPrio);
  XISetPrioSet_Reset(xferGlobals.setXISetPrio);

  if (arrayXIDelCmd!=NULL) OO_Free (arrayXIDelCmd,0);
  FreeAllXIDelCmd();

  if (arrayXIDelObj!=NULL) OO_Free (arrayXIDelObj,0);
  FreeAllXIDelObj();

  if (arrayXINewCpl!=NULL) OO_Free (arrayXINewCpl,0);
  FreeAllXINewCpl();

  if (arrayXIOldCpl!=NULL) OO_Free (arrayXIOldCpl,0);
  FreeAllXIOldCpl();

  if (arrayXIDelCpl!=NULL) OO_Free (arrayXIDelCpl,0);
  FreeAllXIDelCpl();

  if (arrayXIModCpl!=NULL) OO_Free (arrayXIModCpl,0);
  FreeAllXIModCpl();

  if (arrayXIAddCpl!=NULL) OO_Free (arrayXIAddCpl,0);
  FreeAllXIAddCpl();

  FreeLocalCoupledObjectsList(localCplObjs);

  for(; sendMsgs!=NULL; sendMsgs=sm)
  {
    sm = sendMsgs->next;
    OO_Free (sendMsgs,0);
  }

        #ifdef XferMemFromHeap
  xferGlobals.useHeap = FALSE;
  ReleaseHeap(xferGlobals.theMarkKey);
  LC_SetMemMgrDefault();
        #endif

#       if DebugXfer<=4
  sprintf(cBuffer,"%4d: XferEnd, before IFAllFromScratch().\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  /* re-create all interfaces and step XMODE */
  STAT_RESET;
  IFAllFromScratch();
  STAT_TIMER(T_XFER_BUILD_IF);

  XferStepMode(XMODE_BUSY);
}




/* ablage fuer debug-ausgabe

 #	if DebugXfer<=4
                sprintf(cBuffer,"%4d: XferEnd, after XferDeleteObjects(): "
                        "send=%d recv=%d\n", me, nSendMsgs, nRecvMsgs);
                DDD_PrintDebug(cBuffer);
 #	endif
 */



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferPrioChange                                            */
/*                                                                          */
/****************************************************************************/

/**
        Consistent change of a local object's priority during DDD Transfer.
        Local objects which are part of a distributed object must notify
        other copies about local priority changes. This is accomplished
        by issueing \funk{XferPrioChange}-commands during the transfer phase;
        DDD will send appropriate messages to the owner processors of
        the other copies.

        This function is regarded as a {\bf Transfer}-operation due
        to its influence on DDD management information on neighbouring
        processors. Therefore the function has to be issued between
        a starting \funk{XferBegin} and a final \funk{XferEnd} call.

   @param hdr  DDD local object whose priority should be changed.
   @param prio new priority of that local object.
 */



void DDD_XferPrioChange (DDD_HDR hdr, DDD_PRIO prio)
{
  XISetPrio *xi = XISetPrioSet_NewItem(xferGlobals.setXISetPrio);
  xi->hdr  = hdr;
  xi->gid  = OBJ_GID(hdr);
  xi->prio = prio;

  if (! XISetPrioSet_ItemOK(xferGlobals.setXISetPrio))
    return;

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
    XISetPrio *xi = XISetPrioSet_NewItem(xferGlobals.setXISetPrio);
    xi->hdr  = hdr;
    xi->gid  = OBJ_GID(hdr);
    xi->prio = prio;

    if (! XISetPrioSet_ItemOK(xferGlobals.setXISetPrio))
    {
      /* item has been inserted already, don't store it twice. */
      /* even don't call XFERCOPY-handler, this is a real API change! */

      /* xi->prio will be set to PRIO_INVALID if the priority of
         the previously existing XICopyObj-item wins the PriorityMerge
         in the corresponding Compare function (see supp.c). Then,
         we in fact won't need calling the XFERCOPY-handler here,
         because it doesn't give new information. If xi->prio is
         not PRIO_INVALID, then the XICopyObj-item xi wins the merge
         and the XFERCOPY-handler has to be called a second time, now
         with a higher priority. */
      if (xi->prio==PRIO_INVALID)
        return;
    }


    /* although XferCopyObj degrades to SetPrio, call XFERCOPY-handler! */

    /* reset for eventual AddData-calls during handler execution */
    theXIAddData = NULL;

    /* call application handler for xfer of dependent objects */
    if (desc->handlerXFERCOPY)
    {
      DDD_OBJ obj = HDR2OBJ(hdr,desc);

                        #if defined(C_FRONTEND) || defined(F_FRONTEND)
      desc->handlerXFERCOPY(_FADR obj, _FADR dest, _FADR prio);
                        #endif
                        #ifdef CPP_FRONTEND
      CallHandler(desc,XFERCOPY) (HParam(obj) dest, prio);
                        #endif
    }

    /* theXIAddData might be changed during handler execution */
    theXIAddData = NULL;
  }
  else
  {
    /* this is a real transfer to remote proc */
    XICopyObj  *xi = XICopyObjSet_NewItem(xferGlobals.setXICopyObj);
    xi->hdr  = hdr;
    xi->gid  = OBJ_GID(hdr);
    xi->dest = dest;
    xi->prio = prio;

    if (! XICopyObjSet_ItemOK(xferGlobals.setXICopyObj))
    {
      /* item has been inserted already, don't store it twice. */
      /* even don't call XFERCOPY-handler, this is a real API change! */

      /* xi->prio will be set to PRIO_INVALID if the priority of
         the previously existing XICopyObj-item wins the PriorityMerge
         in the corresponding Compare function (see supp.c). Then,
         we in fact won't need calling the XFERCOPY-handler here,
         because it doesn't give new information. If xi->prio is
         not PRIO_INVALID, then the XICopyObj-item xi wins the merge
         and the XFERCOPY-handler has to be called a second time, now
         with a higher priority. */
      if (xi->prio==PRIO_INVALID)
        return;
    }

    xi->size = size;
    xi->add    = NULL;
    xi->addLen = 0;

    /* set XferAddInfo for evtl AddData-calls during handler execution */
    theXIAddData = xi;

    /* call application handler for xfer of dependent objects */
    if (desc->handlerXFERCOPY)
    {
      DDD_OBJ obj = HDR2OBJ(hdr,desc);

                        #if defined(C_FRONTEND) || defined(F_FRONTEND)
      desc->handlerXFERCOPY(_FADR obj, _FADR dest, _FADR prio);
                        #endif
                        #ifdef CPP_FRONTEND
      CallHandler(desc,XFERCOPY) (HParam(obj) dest, prio);
                        #endif
    }

    /* theXIAddData might be changed during handler execution */
    theXIAddData = xi;
  }
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferCopyObj                                               */
/*                                                                          */
/****************************************************************************/

/**
        Transfer-command for copying a local DDD object to another processor.
        After an initial call to \funk{XferBegin}, this function
        creates a copy of one local DDD object on another processor with a certain
        priority. The necessary actions (packing/unpacking of object data, message
        transfer) are executed via the final call to \funk{XferEnd}; therefore
        a whole set of {\bf Transfer}-operations is accumulated.

        Caution: As the original object data is not copied throughout this
        call due to efficiency reasons (transferring a large number of objects
        would result in a huge amount of memory copy operations),
        the object may not be changed or deleted until the actual transfer
        has happened. Otherwise the changes will be sent, too.

        Two different mechanisms allow the transfer of data depending on the
        stated object:

        \begin{itemize}
        \item Right after the function \funk{XferCopyObj} has been called,
        the optional handler #HANDLER_XFERCOPY# is executed by DDD.
        Basically this mechanism allows the transfer of dependent DDD objects
        (or: hierarchies of objects), although arbitrary actions may occur
        inside the handler.
        No relationship between the {\em primary} object and the additional
        objects can be expressed, as many different destination processors
        might be involved.
        %
        \item Via an arbitrary number of additional calls to \funk{XferAddData}
        arrays of {\em data objects} (i.e.~without the usual DDD object header)
        may be sent.
        Using this mechanism the data objects are strongly linked to the
        {\em primary} object; all data objects are transferred to the same
        destination processor.
        \end{itemize}

        After the object copy has been established on the destination processor,
        first the optional handler #HANDLER_LDATACONSTRUCTOR# is called to
        update the LDATA parts of that object (e.g.~locally linked lists).
        Afterwards the optional handler #HANDLER_OBJMKCONS# is called
        to establish consistency for that object (e.g.~backward references on
        the new object copy). For handlers related to function
        \funk{XferAddData} refer to its description.

   @param hdr   DDD local object which has to be copied.
   @param proc  destination processor which will receive the object copy.
   @param prio  DDD priority of new object copy.
 */

#ifdef C_FRONTEND
void DDD_XferCopyObj (DDD_HDR hdr, DDD_PROC proc, DDD_PRIO prio)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Object::XferCopyObj (DDD_PROC proc, DDD_PRIO prio)
{
  DDD_HDR hdr = this;
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




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferCopyObjX                                              */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND)
/**
        Transfer-command for objects of varying sizes.
        This function is an extension of \funk{XferCopyObj}.
        For objects with same DDD type but with variable size in memory,
        one can give the real size as forth parameter to \funk{XferCopyObjX}.
        The DDD Transfer module will use that size value instead of using
        the size computed in \funk{TypeDefine} after the definition of
        the object's DDD type.

   @param hdr   DDD local object which has to be copied.
   @param proc  destination processor which will receive the object copy.
   @param prio  DDD priority of new object copy.
   @param size  real size of local object.
 */
#endif


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



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferAddData                                               */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
/**
        Transfer array of additional data objects with a DDD local object.
        This function transfers an array of additional data objects
        corresponding to one DDD object to the same destination processor.
        Therefore the latest call to \funk{XferCopyObj} defines the
        {\em primary} object and the destination. An arbitrary number of
        \funk{XferAddData}-calls may be issued to transfer various
        data object arrays at once. This serves as a mechanism to send
        object references which cannot be handled by the native DDD
        pointer conversion technique.

        Just before the actual {\bf Transfer}-operation is executed, the (optional)
        handler #HANDLER_XFERGATHER# will be called to fill the data
        object array into some reserved storage space.
        After the {\bf Xfer}-operation the handler #HANDLER_XFERSCATTER# is
        called on the receiving processor to rebuild the data objects.

        As the data objects had to be registered and equipped with a
        usual {\em type\_id} (as with the DDD objects), the
        standard DDD pointer conversion will also take place for the data objects.

   @param cnt  number of data objects in array (i.e.~array length).
   @param typ  DDD type of data objects. All data objects inside one array
        should have the same object type. This object type is defined by
        registering the object structure via \funk{TypeDefine} as usual,
        but without including the DDD object header.
 */

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
  if (xa==NULL)
    HARD_EXIT;

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



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferAddDataX                                              */
/*                                                                          */
/****************************************************************************/

/**
        Transfer array of additional, variable-sized data objects.
        \todo{not documented yet.}
 */

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
  if (xa==NULL)
    HARD_EXIT;

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


#ifdef C_FRONTEND
/**
        Tell application if additional data will be sent.
        If the application issues a \funk{XferCopyObj} command
        with the local processor as the destination processor,
        no object copy will be created. Therefore, also additional
        data objects will not be sent and the corresponding
   #XFER_GATHER#/#XFER_SCATTER# handlers will not be called.
        This function returns a boolean value indicating whether
        the last \funk{XferCopyObj}-command will need the specification
        of additional data objects.

        The application program may use this function in order to
        avoid unnecessary work, \eg, for counting the number of
        additional data objects.

   @return #TRUE# if additional data objects will be gathered, sent
                and scattered; #FALSE# otherwise.
 */
int DDD_XferWithAddData (void)
{
  /* if theXIAddData==NULL, the XferAddData-functions will
     do nothing -> the Gather/Scatter-handlers will not be
     called. */
  return(theXIAddData!=NULL);
}
#endif



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferDeleteObj                                             */
/*                                                                          */
/****************************************************************************/

/**
        Transfer-command for deleting a local DDD object.
        This function is regarded as a {\bf Transfer}-operation due
        to its influence on DDD management information on neighbouring
        processors. Therefore the function has to be issued between
        a starting \funk{XferBegin} and a final \funk{XferEnd} call.

        During the actual {\bf Transfer}-operation all data corresponding
        to that object and the object memory itself will be deleted.

   @param hdr   DDD local object which has to be deleted.
 */

#ifdef C_FRONTEND
void DDD_XferDeleteObj (DDD_HDR hdr)
#endif
#ifdef CPP_FRONTEND
void DDD_Object::XferDeleteObj (void)
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
{
        #ifdef CPP_FRONTEND
  DDD_HDR hdr = this;
        #endif
  TYPE_DESC *desc =  &(theTypeDefs[OBJ_TYPE(hdr)]);
  XIDelCmd  *dc = NewXIDelCmd(SLLNewArgs);

  if (dc==NULL)
    HARD_EXIT;


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
  {
                #if defined(C_FRONTEND) || defined(F_FRONTEND)
    desc->handlerXFERDELETE(HDR2OBJ(hdr,desc));
                #endif
                #ifdef CPP_FRONTEND
    CallHandler(desc,XFERDELETE) (HParamOnly(HDR2OBJ(hdr,desc)));
                #endif
  }
}
#endif


#ifdef F_FRONTEND
void DDD_XferDeleteObj (DDD_TYPE *type, DDD_OBJ *obj)
{
  DDD_HDR hdr = OBJ2HDR(*obj,&theTypeDefs[*type]);
  XIDelCmd  *dc = NewXIDelCmd(SLLNewArgs);

  if (dc==NULL)
    HARD_EXIT;

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
/****************************************************************************/

/**
        Starts transfer phase.
        A call to this function establishes a global transfer operation.
        It should be issued on all processors. After this call an arbitrary
        series of {\bf Xfer}-commands may be issued. The global transfer operation
        is carried out via a \funk{XferEnd} call on each processor.
 */

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


  /* set kind of TMEM alloc/free requests */
  xfer_SetTmpMem(TMEM_XFER);

        #ifdef XferMemFromHeap
  MarkHeap(&xferGlobals.theMarkKey);
  xferGlobals.useHeap = TRUE;
  LC_SetMemMgrRecv(xfer_AllocHeap, NULL);
  LC_SetMemMgrSend(xfer_AllocSend, xfer_FreeSend);
        #endif
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferIsPrunedDelete                                        */
/*                                                                          */
/****************************************************************************/

/**
    Returns information about pruned DDD_XferDeleteObj() command.
    If a \funk{XferDeleteObj} command has been pruned (i.e., option
    OPT_XFER_PRUNE_DELETE is set to OPT_ON and another processor issued
    a \funk{XferCopyObj}, command which sends an object copy to the
    local processor), then this function will return XFER_PRUNED_TRUE,
    otherwise it returns XFER_PRUNED_FALSE. If an error condition
    occurs (e.g., when it is called at the wrong time), the function
    returns XFER_PRUNED_ERROR.


   @param hdr   DDD local object which has to be deleted.
   @return  one of XFER_PRUNED_xxx
 */

#ifdef C_FRONTEND
int DDD_XferIsPrunedDelete (DDD_HDR hdr)
{
  if (XferMode() != XMODE_BUSY)
  {
    return(XFER_PRUNED_ERROR);
  }

  if (OBJ_PRUNED(hdr))
    return(XFER_PRUNED_TRUE);

  return(XFER_PRUNED_FALSE);
}
#endif




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_XferObjIsResent                                           */
/*                                                                          */
/****************************************************************************/

/**
    Returns if an object will receive an additional copy.
    If another processor issued a \funk{XferCopyObj} in order
    to send a further copy of the given local object to the local
    processor, this function will return XFER_RESENT_TRUE.
    Otherwise this function will return XFER_RESENT_FALSE.
    If an error condition occurs (e.g., when it is called at the
    wrong time), the function returns XFER_RESENT_ERROR.

    This function will only work with option OPT_XFER_PRUNE_DELETE
    set to OPT_ON. Otherwise, there will be no kind of communication
    to supply this information.

   @param hdr   DDD local object which has to be deleted.
   @return  one of XFER_RESENT_xxx
 */

#ifdef SUPPORT_RESENT_FLAG

#ifdef C_FRONTEND
int DDD_XferObjIsResent (DDD_HDR hdr)
{
  if (XferMode() != XMODE_BUSY)
  {
    return(XFER_RESENT_ERROR);
  }

  if (DDD_GetOption(OPT_XFER_PRUNE_DELETE)==OPT_OFF)
  {
    return(XFER_RESENT_ERROR);
  }

  if (OBJ_RESENT(hdr))
    return(XFER_RESENT_TRUE);

  return(XFER_RESENT_FALSE);
}
#endif

#endif /* SUPPORT_RESENT_FLAG */


/****************************************************************************/

#undef _FADR

/****************************************************************************/
