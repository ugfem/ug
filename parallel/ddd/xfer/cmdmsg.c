// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      cmdmsg.c                                                      */
/*                                                                          */
/* Purpose:   ddd command transfer:                                         */
/*               send messages with commands to owners of other local       */
/*               object copies. this is used for sending XferCopy commands  */
/*               to owners of local copies, in order to prevent deletion    */
/*               and creation of the same object copy during xfer.          */
/*               (used only if OPT_XFER_PRUNE_DELETE is set to OPT_ON)      */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   970411 kb  created                                            */
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
#include <assert.h>


#include "dddi.h"
#include "xfer.h"





/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

/* CMDMSG: complete description of message on sender side */

typedef struct _CMDMSG
{
  DDD_PROC proc;

  struct _CMDMSG *next;


  DDD_GID   *aUnDelete;
  int nUnDelete;

  /* lowcomm message handle */
  LC_MSGHANDLE msg_h;

} CMDMSG;


/****************************************************************************/
/*                                                                          */
/* variables global to this source file only (static)                       */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



static LC_MSGTYPE cmdmsg_t;
static LC_MSGCOMP undelete_id;



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


void CmdMsgInit (void)
{
  cmdmsg_t = LC_NewMsgType("CmdMsg");
  undelete_id = LC_NewMsgTable("UndelTab", cmdmsg_t, sizeof(DDD_GID));
}


void CmdMsgExit (void)
{}



/****************************************************************************/

static CMDMSG *CreateCmdMsg (DDD_PROC dest, CMDMSG *lastxm)
{
  CMDMSG *xm;

  xm = (CMDMSG *) AllocTmp(sizeof(CMDMSG));
  if (xm==NULL)
  {
    DDD_PrintError('E', 6500, STR_NOMEM " in PrepareCmdMsgs");
    HARD_EXIT;
  }

  xm->aUnDelete = NULL;
  xm->nUnDelete = 0;
  xm->proc = dest;
  xm->next = lastxm;

  return(xm);
}



static int PrepareCmdMsgs (XICopyObj **itemsCO, int nCO, CMDMSG **theMsgs)
{
  CMDMSG    *xm=NULL;
  int j, iCO, markedCO, nMsgs=0;
  DDD_GID   *gids;


  /* set output value in case of early exit */
  *theMsgs = NULL;


  if (nCO==0)
    return(0);

#       if DebugCmdMsg<=3
  printf("%4d: PreparePrune, nCopyObj=%d\n", me, nCO);
  fflush(stdout);
#       endif


  /*
          run through CopyObj table, mark all entries that have
          a coupling with same proc as destination.
   */

  markedCO = 0;
  for(iCO=0; iCO<nCO; iCO++)
  {
    XICopyObj *co = itemsCO[iCO];
    DDD_PROC pCO = co->dest;
    COUPLING *cpl;

    /* run through coupling list of corresponding obj,
       find coupling to destination of XferCopyObj command */
    cpl = ObjCplList(co->hdr);
    while (cpl!=NULL && CPL_PROC(cpl)!=pCO)
      cpl=CPL_NEXT(cpl);

    if (cpl!=NULL)
    {
      /* found coupling -> mark CopyObj */
      SET_CO_SELF(co, 1);
      markedCO++;
    }
    else
    {
      SET_CO_SELF(co, 0);
    }
  }


  if (markedCO==0)
    return(0);


  gids = (DDD_GID *) AllocTmp(sizeof(DDD_GID) * markedCO);
  if (gids==NULL)
  {
    DDD_PrintError('E', 6501, STR_NOMEM " in PrepareCmdMsgs");
    HARD_EXIT;
  }


  /*
          run through CopyObj table again, consider only marked items,
          each time a new proc-nr is encountered in one of these
          tables, create a new CMDMSG item.

          (the lists have been sorted according to proc-nr previously.)
   */
  j=0;
  for(iCO=0; iCO<nCO; iCO++)
  {
    XICopyObj *co = itemsCO[iCO];

    if (CO_SELF(co))
    {
      gids[j] = co->gid;

      if ((xm==NULL) || (xm->proc!=co->dest))
      {
        xm = CreateCmdMsg(co->dest, xm);
        nMsgs++;


        xm->aUnDelete = gids+j;
      }

      xm->nUnDelete++;
      j++;
    }
  }
  *theMsgs = xm;


  /* initiate send messages */
  for(xm=*theMsgs; xm!=NULL; xm=xm->next)
  {
    DDD_GID *array;

    /* create new send message */
    xm->msg_h = LC_NewSendMsg(cmdmsg_t, xm->proc);

    /* init tables inside message */
    LC_SetTableSize(xm->msg_h, undelete_id, xm->nUnDelete);

    /* prepare message for sending away */
    LC_MsgPrepareSend(xm->msg_h);

    array = (DDD_GID *)LC_GetPtr(xm->msg_h, undelete_id);
    memcpy((char *)array,
           (char *)xm->aUnDelete,
           sizeof(DDD_GID)*xm->nUnDelete);
  }

  FreeTmp(gids,0);

  return(nMsgs);
}


static void CmdMsgSend (CMDMSG *theMsgs)
{
  CMDMSG *m;

  for(m=theMsgs; m!=NULL; m=m->next)
  {
    /* schedule message for send */
    LC_MsgSend(m->msg_h);
  }
}



/****************************************************************************/


static int sort_Gids (const void *e1, const void *e2)
{
  if ((*(DDD_GID *)e1) < (*(DDD_GID *)e2)) return(-1);
  if ((*(DDD_GID *)e1) > (*(DDD_GID *)e2)) return(1);

  return(0);
}





static int CmdMsgUnpack (LC_MSGHANDLE *theMsgs, int nRecvMsgs,
                         XIDelCmd  **itemsDC, int nDC)
{
  DDD_GID *unionGidTab;
  int lenGidTab;
  int i, k, jDC, iDC, pos, nPruned;

  lenGidTab = 0;
  for(i=0; i<nRecvMsgs; i++)
  {
    LC_MSGHANDLE xm = theMsgs[i];
    lenGidTab += (int)LC_GetTableLen(xm, undelete_id);
  }

  if (lenGidTab==0)
    return(0);


  unionGidTab = (DDD_GID *) OO_Allocate (sizeof(DDD_GID)*lenGidTab);
  if (unionGidTab==NULL)
  {
    DDD_PrintError('E', 6510, STR_NOMEM " in CmdMsgUnpack");
    HARD_EXIT;
  }


  for(i=0, pos=0; i<nRecvMsgs; i++)
  {
    LC_MSGHANDLE xm = theMsgs[i];
    int len = LC_GetTableLen(xm, undelete_id);

    if (len>0)
    {
      memcpy((char *) (unionGidTab+pos),
             (char *)    LC_GetPtr(xm, undelete_id),
             sizeof(DDD_GID) * len);
      pos += len;
    }
  }
  assert(pos==lenGidTab);

  /* sort GidTab */
  qsort(unionGidTab, lenGidTab, sizeof(DDD_GID), sort_Gids);


        #ifdef SUPPORT_RESENT_FLAG
  {
    DDD_HDR *localCplObjs = LocalCoupledObjectsList();
    int iLCO, nLCO=NCpl_Get;

    /* set RESENT flag for objects which will receive another copy */
    iLCO=0;
    for(k=0; k<lenGidTab; k++)
    {
      while (iLCO<nLCO &&
             (OBJ_GID(localCplObjs[iLCO]) < unionGidTab[k]))
      {
        SET_OBJ_RESENT(localCplObjs[iLCO], 0);

#                               if DebugCmdMsg<=0
        printf("%4d: PruneDelCmds. %08x without resent.\n", me, OBJ_GID(localCplObjs[iLCO]));
        fflush(stdout);
#                               endif

        iLCO++;
      }

      if (iLCO<nLCO && (OBJ_GID(localCplObjs[iLCO]) == unionGidTab[k]))
      {
        SET_OBJ_RESENT(localCplObjs[iLCO], 1);

#                               if DebugCmdMsg<=1
        printf("%4d: PruneDelCmds. %08x will be resent.\n", me, OBJ_GID(localCplObjs[iLCO]));
        fflush(stdout);
#                               endif

        iLCO++;
      }
    }

    /* reset remaining objects' flags */
    while (iLCO<nLCO)
    {
      SET_OBJ_RESENT(localCplObjs[iLCO], 0);

#                       if DebugCmdMsg<=0
      printf("%4d: PruneDelCmds. %08x without resent.\n", me, OBJ_GID(localCplObjs[iLCO]));
      fflush(stdout);
#                       endif

      iLCO++;
    }

    FreeLocalCoupledObjectsList(localCplObjs);
  }
        #endif


  k=0; jDC=0;
  for(iDC=0; iDC<nDC; iDC++)
  {
    DDD_GID gidDC = OBJ_GID(itemsDC[iDC]->hdr);

    while (k<lenGidTab && (unionGidTab[k]<gidDC))
      k++;

    if (k<lenGidTab && (unionGidTab[k]==gidDC))
    {
      /* found a DelCmd-item to prune */
      SET_OBJ_PRUNED(itemsDC[iDC]->hdr, 1);

#                       if DebugCmdMsg<=1
      printf("%4d: PruneDelCmds. pruned %08x\n", me, gidDC);
      fflush(stdout);
#                       endif
    }
    else
    {
      itemsDC[jDC] = itemsDC[iDC];
      jDC++;
    }
  }
  nPruned = nDC-jDC;

#       if DebugCmdMsg<=3
  printf("%4d: PruneDelCmds. nPruned=%d/%d\n", me, nPruned, nDC);
  fflush(stdout);
#       endif


  OO_Free (unionGidTab /*,0*/);

  return(nPruned);
}


/****************************************************************************/


static void CmdMsgDisplay (char *comment, LC_MSGHANDLE xm)
{
  DDD_GID      *theGid;
  char buf[30];
  int i, proc = LC_MsgGetProc(xm);
  int lenGid = (int) LC_GetTableLen(xm, undelete_id);

  sprintf(buf, " %03d-%s-%03d ", me, comment, proc);

  /* get table addresses inside message */
  theGid = (DDD_GID *)    LC_GetPtr(xm, undelete_id);


  sprintf(cBuffer, "%s 04 Gid.size=%05d\n", buf, lenGid);
  DDD_PrintDebug(cBuffer);

  for(i=0; i<lenGid; i++)
  {
    sprintf(cBuffer, "%s 14 gid    %04d - %08x\n",
            buf, i, (unsigned long)theGid[i]);
    DDD_PrintDebug(cBuffer);
  }
}


/****************************************************************************/


/*
        prune superfluous DelCmds. DelCmds are regarded superfluous,
        if there exists another processor which also owns an object copy
        and sends it to destination=me.

        this is implemented by sending the GID of all objects with an
        XferCopyObj-command with destination=p to p, iff p does already
        own a copy of the object. after receiving the messages,
        each processor prunes DelCmds for the GIDs in the message.

        returns:  number of pruned DelCmds.  ( < nDC)
 */

int PruneXIDelCmd (
  XIDelCmd  **itemsDC, int nDC,
  XICopyObjPtrArray *arrayCO)
{
  CMDMSG    *sendMsgs, *sm=0;
  LC_MSGHANDLE *recvMsgs;
  int i, nSendMsgs, nRecvMsgs;
  int nPruned;

  XICopyObj **itemsCO = XICopyObjPtrArray_GetData(arrayCO);
  int nCO       = XICopyObjPtrArray_GetSize(arrayCO);


  /* accumulate messages (one for each partner) */
  nSendMsgs = PrepareCmdMsgs(itemsCO, nCO, &sendMsgs);

#if DebugCmdMsg>2
  if (DDD_GetOption(OPT_DEBUG_XFERMESGS)==OPT_ON)
#endif
  {
    for(sm=sendMsgs; sm!=NULL; sm=sm->next)
    {
      CmdMsgDisplay("PS", sm->msg_h);
    }
  }

  /* init communication topology */
  nRecvMsgs = LC_Connect(cmdmsg_t);

  /* build and send messages */
  CmdMsgSend(sendMsgs);

  /* communicate set of messages (send AND receive) */
  recvMsgs = LC_Communicate();


  nPruned = CmdMsgUnpack(recvMsgs, nRecvMsgs, itemsDC, nDC);


  /*
          for(i=0; i<nRecvMsgs; i++)
          {
     #if DebugCmdMsg>=2
                  if (DDD_GetOption(OPT_DEBUG_XFERMESGS)==OPT_ON)
     #endif
                          CmdMsgDisplay("PR", recvMsgs[i]);
          }
   */


  /* cleanup low-comm layer */
  LC_Cleanup();


  /* free temporary memory */
  for(; sendMsgs!=NULL; sendMsgs=sm)
  {
    sm = sendMsgs->next;
    FreeTmp(sendMsgs,0);
  }

  return(nPruned);
}


/****************************************************************************/
