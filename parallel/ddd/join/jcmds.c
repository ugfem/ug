// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      jcmds.c                                                       */
/*                                                                          */
/* Purpose:   DDD-commands for Join Environment                             */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: birken@ica3.uni-stuttgart.de                           */
/*            phone: 0049-(0)711-685-7007                                   */
/*            fax  : 0049-(0)711-685-7000                                   */
/*                                                                          */
/* History:   980126 kb  begin                                              */
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
#include "join.h"




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


static int sort_NewGid (const void *e1, const void *e2)
{
  REGISTER JIJoin *item1 = *((JIJoin **)e1);
  REGISTER JIJoin *item2 = *((JIJoin **)e2);

  if (item1->new_gid < item2->new_gid) return(-1);
  if (item1->new_gid > item2->new_gid) return(1);

  return(0);
}


/****************************************************************************/


static int sort_Gid (const void *e1, const void *e2)
{
  REGISTER JIPartner *item1 = (JIPartner *)e1;
  REGISTER JIPartner *item2 = (JIPartner *)e2;

  if (OBJ_GID(item1->hdr) < OBJ_GID(item2->hdr)) return(-1);
  if (OBJ_GID(item1->hdr) > OBJ_GID(item2->hdr)) return(1);

  return(0);
}


/****************************************************************************/


/*
        prepare messages for phase 1.

 */

static int PreparePhase1Msgs (JIJoinPtrArray *arrayJoin,
                              JOINMSG1 **theMsgs, size_t *memUsage)
{
  int i, last_i, nMsgs;
  JIJoin   **itemsJ = JIJoinPtrArray_GetData(arrayJoin);
  int nJ       = JIJoinPtrArray_GetSize(arrayJoin);

#       if DebugJoin<=3
  printf("%4d: PreparePhase1Msgs, nJoins=%d\n",
         me, nJ);
  fflush(stdout);
#       endif

  /* init return parameters */
  *theMsgs = NULL;
  *memUsage = 0;


  if (nJ==0)
    /* no messages */
    return(0);


  /* check whether Join objects are really local (without copies) */
  /* and set local GID to invalid (will be set to new value lateron) */
  for(i=0; i<nJ; i++)
  {
    if (ObjHasCpl(itemsJ[i]->hdr))
    {
      sprintf(cBuffer, "cannot join %08x, object already distributed",
              OBJ_GID(itemsJ[i]->hdr));
      DDD_PrintError('E', 7006, cBuffer);
      HARD_EXIT;
    }

    OBJ_GID(itemsJ[i]->hdr) = GID_INVALID;
  }


  /* set local GID to new value */
  for(i=0; i<nJ; i++)
  {
    DDD_GID local_gid = OBJ_GID(itemsJ[i]->hdr);

    /* check for double Joins with different new_gid */
    if (local_gid!=GID_INVALID && local_gid!=itemsJ[i]->new_gid)
    {
      sprintf(cBuffer,
              "several (inconsistent) DDD_JoinObj-commands for local object %08x",
              local_gid);
      DDD_PrintError('E', 7007, cBuffer);
      HARD_EXIT;
    }

    OBJ_GID(itemsJ[i]->hdr) = itemsJ[i]->new_gid;
  }


  nMsgs = 0;
  last_i = i = 0;
  do
  {
    JOINMSG1  *jm;
    size_t bufSize;

    /* skip until dest-processor is different */
    while (i<nJ && (itemsJ[i]->dest == itemsJ[last_i]->dest))
      i++;

    /* create new message */
    jm = (JOINMSG1 *) AllocTmp(sizeof(JOINMSG1));
    if (jm==NULL)
    {
      DDD_PrintError('E', 7900, STR_NOMEM " in PreparePhase1Msgs");
      HARD_EXIT;
    }
    jm->nJoins = i-last_i;
    jm->arrayJoin = &(itemsJ[last_i]);
    jm->dest = itemsJ[last_i]->dest;
    jm->next = *theMsgs;
    *theMsgs = jm;
    nMsgs++;

    /* create new send message */
    jm->msg_h = LC_NewSendMsg(joinGlobals.phase1msg_t, jm->dest);

    /* init table inside message */
    LC_SetTableSize(jm->msg_h, joinGlobals.jointab_id, jm->nJoins);

    /* prepare message for sending away */
    bufSize = LC_MsgPrepareSend(jm->msg_h);
    *memUsage += bufSize;

    if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MEMUSAGE)
    {
      sprintf(cBuffer,
              "DDD MESG [%03d]: SHOW_MEM "
              "send msg phase1   dest=%04d size=%010ld\n",
              me, jm->dest, (long)bufSize);
      DDD_PrintLine(cBuffer);
    }

    last_i = i;

  } while (last_i < nJ);

  return(nMsgs);
}



/****************************************************************************/
/*                                                                          */
/* Function:  PackPhase1Msgs                                                */
/*                                                                          */
/* Purpose:   allocate one message buffer for each outgoing message,        */
/*            fill buffer with message contents and initiate asynchronous   */
/*            send for each message.                                        */
/*                                                                          */
/* Input:     theMsgs: list of message-send-infos                           */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void PackPhase1Msgs (JOINMSG1 *theMsgs)
{
  JOINMSG1 *jm;

  for(jm=theMsgs; jm!=NULL; jm=jm->next)
  {
    TEJoin *theJoinTab;
    int i;

    /* copy data into message */
    theJoinTab = (TEJoin *)LC_GetPtr(jm->msg_h, joinGlobals.jointab_id);
    for(i=0; i<jm->nJoins; i++)
    {
      theJoinTab[i].gid  = jm->arrayJoin[i]->new_gid;
      theJoinTab[i].prio = OBJ_PRIO(jm->arrayJoin[i]->hdr);
    }
    LC_SetTableLen(jm->msg_h, joinGlobals.jointab_id, jm->nJoins);


    /* send away */
    LC_MsgSend(jm->msg_h);
  }
}



/*
        unpack phase1 messages.
 */

static void UnpackPhase1Msgs (LC_MSGHANDLE *theMsgs, int nRecvMsgs,
                              DDD_HDR *localCplObjs, int nLCO,
                              JIPartner **p_joinObjs, int *p_nJoinObjs)
{
  JIPartner *joinObjs;
  int nJoinObjs = 0;
  int m, jo;

  /* init return values */
  *p_joinObjs  = NULL;
  *p_nJoinObjs = 0;


  for(m=0; m<nRecvMsgs; m++)
  {
    LC_MSGHANDLE jm = theMsgs[m];
    TEJoin *theJoin = (TEJoin *) LC_GetPtr(jm, joinGlobals.jointab_id);
    int nJ       = (int) LC_GetTableLen(jm, joinGlobals.jointab_id);
    int i, j;

    nJoinObjs += nJ;

    for(i=0, j=0; i<nJ; i++)
    {
      while ((j<nLCO) && (OBJ_GID(localCplObjs[j]) < theJoin[i].gid))
        j++;

      if ((j<nLCO) && (OBJ_GID(localCplObjs[j])==theJoin[i].gid))
      {
        COUPLING *cpl;

        /* found local object which is join target */
        /* store shortcut to local object */
        theJoin[i].hdr = localCplObjs[j];

        /* generate phase2-JIAddCpl for this object */
        for(cpl=ObjCplList(localCplObjs[j]); cpl!=NULL; cpl=CPL_NEXT(cpl))
        {
          JIAddCpl *ji = JIAddCplSet_NewItem(joinGlobals.setJIAddCpl2);
          ji->dest    = CPL_PROC(cpl);
          ji->te.gid  = theJoin[i].gid;
          ji->te.proc = LC_MsgGetProc(jm);
          ji->te.prio = theJoin[i].prio;

          if (! JIAddCplSet_ItemOK(joinGlobals.setJIAddCpl2))
            continue;

#                                       if DebugJoin<=1
          printf("%4d: Phase1 Join for %08x from %d, "
                 "send AddCpl to %d.\n",
                 me, theJoin[i].gid, ji->te.proc, ji->dest);
#                                       endif

        }

        /* send phase3-JIAddCpl back to Join-proc */
        for(cpl=ObjCplList(localCplObjs[j]); cpl!=NULL; cpl=CPL_NEXT(cpl))
        {
          JIAddCpl *ji = JIAddCplSet_NewItem(joinGlobals.setJIAddCpl3);
          ji->dest    = LC_MsgGetProc(jm);
          ji->te.gid  = OBJ_GID(localCplObjs[j]);
          ji->te.proc = CPL_PROC(cpl);
          ji->te.prio = cpl->prio;

          if (! JIAddCplSet_ItemOK(joinGlobals.setJIAddCpl3))
            continue;
        }
      }
      else
      {
        sprintf(cBuffer, "no object %08x for join from %d",
                theJoin[i].gid, LC_MsgGetProc(jm));
        DDD_PrintError('E', 7300, cBuffer);
        HARD_EXIT;
      }
    }
  }


  /* return immediately if no join-objects have been found */
  if (nJoinObjs==0)
    return;


  /* allocate array of objects, which has been contacted by a join */
  joinObjs = (JIPartner *) AllocTmp(sizeof(JIPartner) * nJoinObjs);
  if (joinObjs==NULL)
  {
    DDD_PrintError('E', 7903, STR_NOMEM " in UnpackPhase1Msgs");
    HARD_EXIT;
  }

  /* set return values */
  *p_joinObjs  = joinObjs;
  *p_nJoinObjs = nJoinObjs;


  /* add one local coupling for each Join */
  for(m=0, jo=0; m<nRecvMsgs; m++)
  {
    LC_MSGHANDLE jm = theMsgs[m];
    TEJoin *theJoin = (TEJoin *) LC_GetPtr(jm, joinGlobals.jointab_id);
    int nJ       = (int) LC_GetTableLen(jm, joinGlobals.jointab_id);
    int i;

    for(i=0; i<nJ; i++)
    {
      AddCoupling(theJoin[i].hdr, LC_MsgGetProc(jm), theJoin[i].prio);

      /* send one phase3-JIAddCpl for symmetric connection */
      {
        JIAddCpl *ji = JIAddCplSet_NewItem(joinGlobals.setJIAddCpl3);
        ji->dest    = LC_MsgGetProc(jm);
        ji->te.gid  = OBJ_GID(theJoin[i].hdr);
        ji->te.proc = me;
        ji->te.prio = OBJ_PRIO(theJoin[i].hdr);

        JIAddCplSet_ItemOK(joinGlobals.setJIAddCpl3);
      }

      joinObjs[jo].hdr  = theJoin[i].hdr;
      joinObjs[jo].proc = LC_MsgGetProc(jm);
      jo++;
    }
  }


  /* sort joinObjs-array according to gid */
  if (nJoinObjs>1)
    qsort(joinObjs, nJoinObjs, sizeof(JIPartner), sort_Gid);
}



/****************************************************************************/

/*
        prepare messages for phase 2.

 */

static int PreparePhase2Msgs (JIAddCplPtrArray *arrayAddCpl,
                              JOINMSG2 **theMsgs, size_t *memUsage)
{
  int i, last_i, nMsgs;
  JIAddCpl **itemsAC = JIAddCplPtrArray_GetData(arrayAddCpl);
  int nAC       = JIAddCplPtrArray_GetSize(arrayAddCpl);

#       if DebugJoin<=3
  printf("%4d: PreparePhase2Msgs, nAddCpls=%d\n",
         me, nAC);
  fflush(stdout);
#       endif

  /* init return parameters */
  *theMsgs = NULL;
  *memUsage = 0;


  if (nAC==0)
    /* no messages */
    return(0);


  nMsgs = 0;
  last_i = i = 0;
  do
  {
    JOINMSG2  *jm;
    size_t bufSize;

    /* skip until dest-processor is different */
    while (i<nAC && (itemsAC[i]->dest == itemsAC[last_i]->dest))
      i++;

    /* create new message */
    jm = (JOINMSG2 *) AllocTmp(sizeof(JOINMSG2));
    if (jm==NULL)
    {
      DDD_PrintError('E', 7901, STR_NOMEM " in PreparePhase2Msgs");
      HARD_EXIT;
    }
    jm->nAddCpls = i-last_i;
    jm->arrayAddCpl = &(itemsAC[last_i]);
    jm->dest = itemsAC[last_i]->dest;
    jm->next = *theMsgs;
    *theMsgs = jm;
    nMsgs++;

    /* create new send message */
    jm->msg_h = LC_NewSendMsg(joinGlobals.phase2msg_t, jm->dest);

    /* init table inside message */
    LC_SetTableSize(jm->msg_h, joinGlobals.addtab_id, jm->nAddCpls);

    /* prepare message for sending away */
    bufSize = LC_MsgPrepareSend(jm->msg_h);
    *memUsage += bufSize;

    if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MEMUSAGE)
    {
      sprintf(cBuffer,
              "DDD MESG [%03d]: SHOW_MEM "
              "send msg phase2   dest=%04d size=%010ld\n",
              me, jm->dest, (long)bufSize);
      DDD_PrintLine(cBuffer);
    }

    last_i = i;

  } while (last_i < nAC);

  return(nMsgs);
}



/****************************************************************************/
/*                                                                          */
/* Function:  PackPhase2Msgs                                                */
/*                                                                          */
/* Purpose:   allocate one message buffer for each outgoing message,        */
/*            fill buffer with message contents and initiate asynchronous   */
/*            send for each message.                                        */
/*                                                                          */
/* Input:     theMsgs: list of message-send-infos                           */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void PackPhase2Msgs (JOINMSG2 *theMsgs)
{
  JOINMSG2 *jm;

  for(jm=theMsgs; jm!=NULL; jm=jm->next)
  {
    TEAddCpl *theAddTab;
    int i;

    /* copy data into message */
    theAddTab = (TEAddCpl *)LC_GetPtr(jm->msg_h, joinGlobals.addtab_id);
    for(i=0; i<jm->nAddCpls; i++)
    {
      /* copy complete TEAddCpl item */
      theAddTab[i] = jm->arrayAddCpl[i]->te;
    }
    LC_SetTableLen(jm->msg_h, joinGlobals.addtab_id, jm->nAddCpls);


    /* send away */
    LC_MsgSend(jm->msg_h);
  }
}



/*
        unpack phase2 messages.
 */

static void UnpackPhase2Msgs (LC_MSGHANDLE *theMsgs2, int nRecvMsgs2,
                              JIPartner *joinObjs, int nJoinObjs,
                              DDD_HDR *localCplObjs, int nLCO)
{
  int m;

  for(m=0; m<nRecvMsgs2; m++)
  {
    LC_MSGHANDLE jm = theMsgs2[m];
    TEAddCpl *theAC = (TEAddCpl *) LC_GetPtr(jm, joinGlobals.addtab_id);
    int nAC      = (int) LC_GetTableLen(jm, joinGlobals.addtab_id);
    int i, j, jo;

    for(i=0, j=0, jo=0; i<nAC; i++)
    {
      while ((j<nLCO) && (OBJ_GID(localCplObjs[j]) < theAC[i].gid))
        j++;

      while ((jo<nJoinObjs) && (OBJ_GID(joinObjs[jo].hdr) < theAC[i].gid))
        jo++;

      if ((j<nLCO) && (OBJ_GID(localCplObjs[j])==theAC[i].gid))
      {
        /* found local object which is AddCpl target */
        AddCoupling(localCplObjs[j], theAC[i].proc, theAC[i].prio);

#                               if DebugJoin<=1
        printf("%4d: Phase2 execute AddCpl(%08x,%d,%d) (from %d).\n",
               me,
               theAC[i].gid, theAC[i].proc, theAC[i].prio,
               LC_MsgGetProc(jm));
#                               endif

        while ((jo<nJoinObjs) && (OBJ_GID(joinObjs[jo].hdr) == theAC[i].gid))
        {
          JIAddCpl *ji = JIAddCplSet_NewItem(joinGlobals.setJIAddCpl3);
          ji->dest    = joinObjs[jo].proc;
          ji->te.gid  = theAC[i].gid;
          ji->te.proc = theAC[i].proc;
          ji->te.prio = theAC[i].prio;
          JIAddCplSet_ItemOK(joinGlobals.setJIAddCpl3);


#                                       if DebugJoin<=1
          printf("%4d: Phase2 forward AddCpl(%08x,%d,%d) to %d.\n",
                 me,
                 theAC[i].gid, theAC[i].proc, theAC[i].prio,
                 ji->dest);
#                                       endif

          jo++;
        }
      }
      else
      {
        /* this should never happen. AddCpl send from invalid proc. */
        assert(0);
      }
    }
  }
}



/****************************************************************************/



/*
        prepare messages for phase 3.

 */

static int PreparePhase3Msgs (JIAddCplPtrArray *arrayAddCpl,
                              JOINMSG3 **theMsgs, size_t *memUsage)
{
  int i, last_i, nMsgs;
  JIAddCpl **itemsAC = JIAddCplPtrArray_GetData(arrayAddCpl);
  int nAC       = JIAddCplPtrArray_GetSize(arrayAddCpl);

#       if DebugJoin<=3
  printf("%4d: PreparePhase3Msgs, nAddCpls=%d\n",
         me, nAC);
  fflush(stdout);
#       endif

  /* init return parameters */
  *theMsgs = NULL;
  *memUsage = 0;


  if (nAC==0)
    /* no messages */
    return(0);


  nMsgs = 0;
  last_i = i = 0;
  do
  {
    JOINMSG3  *jm;
    size_t bufSize;

    /* skip until dest-processor is different */
    while (i<nAC && (itemsAC[i]->dest == itemsAC[last_i]->dest))
      i++;

    /* create new message */
    jm = (JOINMSG3 *) AllocTmp(sizeof(JOINMSG3));
    if (jm==NULL)
    {
      DDD_PrintError('E', 7902, STR_NOMEM " in PreparePhase3Msgs");
      HARD_EXIT;
    }
    jm->nAddCpls = i-last_i;
    jm->arrayAddCpl = &(itemsAC[last_i]);
    jm->dest = itemsAC[last_i]->dest;
    jm->next = *theMsgs;
    *theMsgs = jm;
    nMsgs++;

    /* create new send message */
    jm->msg_h = LC_NewSendMsg(joinGlobals.phase3msg_t, jm->dest);

    /* init table inside message */
    LC_SetTableSize(jm->msg_h, joinGlobals.cpltab_id, jm->nAddCpls);

    /* prepare message for sending away */
    bufSize = LC_MsgPrepareSend(jm->msg_h);
    *memUsage += bufSize;

    if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MEMUSAGE)
    {
      sprintf(cBuffer,
              "DDD MESG [%03d]: SHOW_MEM "
              "send msg phase3   dest=%04d size=%010ld\n",
              me, jm->dest, (long)bufSize);
      DDD_PrintLine(cBuffer);
    }

    last_i = i;

  } while (last_i < nAC);

  return(nMsgs);
}



/****************************************************************************/
/*                                                                          */
/* Function:  PackPhase3Msgs                                                */
/*                                                                          */
/* Purpose:   allocate one message buffer for each outgoing message,        */
/*            fill buffer with message contents and initiate asynchronous   */
/*            send for each message.                                        */
/*                                                                          */
/* Input:     theMsgs: list of message-send-infos                           */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void PackPhase3Msgs (JOINMSG3 *theMsgs)
{
  JOINMSG3 *jm;

  for(jm=theMsgs; jm!=NULL; jm=jm->next)
  {
    TEAddCpl *theAddTab;
    int i;

    /* copy data into message */
    theAddTab = (TEAddCpl *)LC_GetPtr(jm->msg_h, joinGlobals.cpltab_id);
    for(i=0; i<jm->nAddCpls; i++)
    {
      /* copy complete TEAddCpl item */
      theAddTab[i] = jm->arrayAddCpl[i]->te;
    }
    LC_SetTableLen(jm->msg_h, joinGlobals.cpltab_id, jm->nAddCpls);


    /* send away */
    LC_MsgSend(jm->msg_h);
  }
}



/****************************************************************************/


/*
        unpack phase3 messages.
 */

static void UnpackPhase3Msgs (LC_MSGHANDLE *theMsgs, int nRecvMsgs,
                              JIJoinPtrArray *arrayJoin)
{
  JIJoin   **itemsJ = JIJoinPtrArray_GetData(arrayJoin);
  int nJ       = JIJoinPtrArray_GetSize(arrayJoin);
  int m;


  for(m=0; m<nRecvMsgs; m++)
  {
    LC_MSGHANDLE jm = theMsgs[m];
    TEAddCpl *theAC = (TEAddCpl *) LC_GetPtr(jm, joinGlobals.cpltab_id);
    int nAC      = (int) LC_GetTableLen(jm, joinGlobals.cpltab_id);
    int i, j;

    for(i=0, j=0; i<nAC; i++)
    {
      while ((j<nJ) && (OBJ_GID(itemsJ[j]->hdr) < theAC[i].gid))
        j++;

      if ((j<nJ) && (OBJ_GID(itemsJ[j]->hdr) == theAC[i].gid))
      {
        /* found local object which is AddCpl target */
        AddCoupling(itemsJ[j]->hdr, theAC[i].proc, theAC[i].prio);

#                               if DebugJoin<=1
        printf("%4d: Phase3 execute AddCpl(%08x,%d,%d) (from %d).\n",
               me,
               OBJ_GID(itemsJ[j]->hdr), theAC[i].proc, theAC[i].prio,
               LC_MsgGetProc(jm));
#                               endif
      }
      else
      {
        /* this should never happen. AddCpl send for unknown obj. */
        assert(0);
      }
    }
  }
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_JoinEnd                                                   */
/*                                                                          */
/****************************************************************************/

/**
        End of join phase.
        This function starts the actual join process. After a call to
        this function (on all processors) all {\bf Join}-commands since
        the last call to \funk{JoinBegin} are executed. This involves
        a set of local communications between the processors.
 */

#ifdef C_FRONTEND
DDD_RET DDD_JoinEnd (void)
#endif
#ifdef CPP_FRONTEND
DDD_RET DDD_Library::JoinEnd (void)
#endif
#ifdef F_FRONTEND
DDD_RET DDD_JoinEnd (void)
#endif
{
  JIJoinPtrArray   *arrayJIJoin    = NULL;
  JIAddCplPtrArray *arrayJIAddCpl2 = NULL;
  JIAddCplPtrArray *arrayJIAddCpl3 = NULL;
  int obsolete, nRecvMsgs1, nRecvMsgs2, nRecvMsgs3, nSendMsgs;
  JOINMSG1    *sendMsgs1=NULL, *sm1=NULL;
  JOINMSG2    *sendMsgs2=NULL, *sm2=NULL;
  JOINMSG3    *sendMsgs3=NULL, *sm3=NULL;
  LC_MSGHANDLE *recvMsgs1=NULL, *recvMsgs2=NULL, *recvMsgs3=NULL;
  DDD_HDR     *localCplObjs=NULL;
  size_t sendMem=0, recvMem=0;
  JIPartner   *joinObjs = NULL;
  int nJoinObjs;



        #ifdef JoinMemFromHeap
  MarkHeap();
  LC_SetMemMgr(memmgr_AllocTMEM, memmgr_FreeTMEM,
               memmgr_AllocHMEM, NULL);
        #endif

  STAT_SET_MODULE(DDD_MODULE_JOIN);
  STAT_ZEROALL;

  /* step mode and check whether call to JoinEnd is valid */
  if (!JoinStepMode(JMODE_CMDS))
  {
    DDD_PrintError('E', 7011, "DDD_JoinEnd() aborted");
    HARD_EXIT;
  }


  /*
          PREPARATION PHASE
   */
  /* get sorted array of JIJoin-items */
  arrayJIJoin = JIJoinSet_GetArray(joinGlobals.setJIJoin);
  obsolete = JIJoinSet_GetNDiscarded(joinGlobals.setJIJoin);


  /*
          COMMUNICATION PHASE 1
          all processors, where JoinObj-commands have been issued,
          send information about these commands to the target
          processors together with the GID of the objects on the
          target procs and the local priority.
   */
  STAT_RESET;
  /* prepare msgs for JIJoin-items */
  nSendMsgs = PreparePhase1Msgs(arrayJIJoin, &sendMsgs1, &sendMem);
  /* DisplayMemResources(); */

  /* init communication topology */
  nRecvMsgs1 = LC_Connect(joinGlobals.phase1msg_t);
  STAT_TIMER(T_JOIN_PREP_MSGS);

  STAT_RESET;
  /* build phase1 msgs on sender side and start send */
  PackPhase1Msgs(sendMsgs1);
  STAT_TIMER(T_JOIN_PACK_SEND);


  /*
          now messages are in the net, use spare time
   */
  STAT_RESET;
  /* get sorted list of local objects with couplings */
  localCplObjs = LocalCoupledObjectsList();
  if (localCplObjs==NULL && ddd_nCpls>0)
  {
    DDD_PrintError('E', 7020,
                   "Cannot get list of coupled objects in DDD_JoinEnd(). Aborted");
    HARD_EXIT;
  }


  if (obsolete>0)
  {
    if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_OBSOLETE)
    {
      int all = JIJoinSet_GetNItems(joinGlobals.setJIJoin);

      sprintf(cBuffer, "DDD MESG [%03d]: %4d from %4d join-cmds obsolete.\n",
              me, obsolete, all);
      DDD_PrintLine(cBuffer);
    }
  }
  STAT_TIMER(T_JOIN);


  /*
          nothing more to do until incoming messages arrive
   */

  /* display information about send-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD JOIN_SHOW_MSGSALL: Phase1Msg.Send\n");
    LC_PrintSendMsgs();
  }


  /* wait for communication-completion (send AND receive) */
  STAT_RESET;
  recvMsgs1 = LC_Communicate();
  STAT_TIMER(T_JOIN_WAIT_RECV);


  /* display information about message buffer sizes */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MEMUSAGE)
  {
    int k;

    /* sum up sizes of receive mesg buffers */
    for(k=0; k<nRecvMsgs1; k++)
    {
      recvMem += LC_GetBufferSize(recvMsgs1[k]);
    }

    sprintf(cBuffer,
            "DDD MESG [%03d]: SHOW_MEM "
            "msgs  send=%010ld recv=%010ld all=%010ld\n",
            me, (long)sendMem, (long)recvMem, (long)(sendMem+recvMem));
    DDD_PrintLine(cBuffer);
  }

  /* display information about recv-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD JOIN_SHOW_MSGSALL: Phase1Msg.Recv\n");
    LC_PrintRecvMsgs();
  }

  /* unpack messages */
  STAT_RESET;
  UnpackPhase1Msgs(recvMsgs1, nRecvMsgs1, localCplObjs, NCpl_Get,
                   &joinObjs, &nJoinObjs);
  LC_Cleanup();
  STAT_TIMER(T_JOIN_UNPACK);





  /*
          COMMUNICATION PHASE 2
          all processors which received notification of JoinObj-commands
          during phase 1 send AddCpl-requests to all copies of DDD objects,
          for which Joins had been issued remotely.
   */
  /* get sorted array of JIAddCpl-items */
  arrayJIAddCpl2 = JIAddCplSet_GetArray(joinGlobals.setJIAddCpl2);

  STAT_RESET;
  /* prepare msgs for JIAddCpl-items */
  nSendMsgs = PreparePhase2Msgs(arrayJIAddCpl2, &sendMsgs2, &sendMem);
  /* DisplayMemResources(); */

  /* init communication topology */
  nRecvMsgs2 = LC_Connect(joinGlobals.phase2msg_t);
  STAT_TIMER(T_JOIN_PREP_MSGS);

  STAT_RESET;
  /* build phase2 msgs on sender side and start send */
  PackPhase2Msgs(sendMsgs2);
  STAT_TIMER(T_JOIN_PACK_SEND);

  /*
          now messages are in the net, use spare time
   */

  /* reorder Join-commands according to new_gid */
  /* this ordering is needed in UnpackPhase3 */
  if (JIJoinPtrArray_GetSize(arrayJIJoin) > 1)
  {
    qsort(
      JIJoinPtrArray_GetData(arrayJIJoin),
      JIJoinPtrArray_GetSize(arrayJIJoin),
      sizeof(JIJoin *), sort_NewGid);
  }


  /*
          nothing more to do until incoming messages arrive
   */

  /* display information about send-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD JOIN_SHOW_MSGSALL: Phase2Msg.Send\n");
    LC_PrintSendMsgs();
  }


  /* wait for communication-completion (send AND receive) */
  STAT_RESET;
  recvMsgs2 = LC_Communicate();
  STAT_TIMER(T_JOIN_WAIT_RECV);


  /* display information about message buffer sizes */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MEMUSAGE)
  {
    int k;

    /* sum up sizes of receive mesg buffers */
    for(k=0; k<nRecvMsgs2; k++)
    {
      recvMem += LC_GetBufferSize(recvMsgs2[k]);
    }

    sprintf(cBuffer,
            "DDD MESG [%03d]: SHOW_MEM "
            "msgs  send=%010ld recv=%010ld all=%010ld\n",
            me, (long)sendMem, (long)recvMem, (long)(sendMem+recvMem));
    DDD_PrintLine(cBuffer);
  }

  /* display information about recv-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD JOIN_SHOW_MSGSALL: Phase2Msg.Recv\n");
    LC_PrintRecvMsgs();
  }

  /* unpack messages */
  STAT_RESET;
  UnpackPhase2Msgs(recvMsgs2, nRecvMsgs2, joinObjs, nJoinObjs,
                   localCplObjs, NCpl_Get);

  LC_Cleanup();
  STAT_TIMER(T_JOIN_UNPACK);

  for(; sendMsgs2!=NULL; sendMsgs2=sm2)
  {
    sm2 = sendMsgs2->next;
    FreeTmp(sendMsgs2, 0);
  }






  /*
          COMMUNICATION PHASE 3
          all processors which received notification of JoinObj-commands
          during phase 1 send AddCpl-requests to the procs where the
          JoinObj-commands have been issued. One AddCpl-request is sent
          for each cpl in the local objects coupling list. One AddCpl-request
          is sent for each AddCpl-request received during phase 2.
          (i.e., two kinds of AddCpl-requests are send to the processors
          on which the JoinObj-commands have been issued.
   */
  /* get sorted array of JIAddCpl-items */
  arrayJIAddCpl3 = JIAddCplSet_GetArray(joinGlobals.setJIAddCpl3);

  STAT_RESET;
  /* prepare msgs for JIAddCpl-items */
  nSendMsgs = PreparePhase3Msgs(arrayJIAddCpl3, &sendMsgs3, &sendMem);
  /* DisplayMemResources(); */

  /* init communication topology */
  nRecvMsgs3 = LC_Connect(joinGlobals.phase3msg_t);
  STAT_TIMER(T_JOIN_PREP_MSGS);

  STAT_RESET;
  /* build phase3 msgs on sender side and start send */
  PackPhase3Msgs(sendMsgs3);
  STAT_TIMER(T_JOIN_PACK_SEND);

  /*
          now messages are in the net, use spare time
   */
  /* ... */

  /*
          nothing more to do until incoming messages arrive
   */

  /* display information about send-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD JOIN_SHOW_MSGSALL: Phase3Msg.Send\n");
    LC_PrintSendMsgs();
  }


  /* wait for communication-completion (send AND receive) */
  STAT_RESET;
  recvMsgs3 = LC_Communicate();
  STAT_TIMER(T_JOIN_WAIT_RECV);


  /* display information about message buffer sizes */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MEMUSAGE)
  {
    int k;

    /* sum up sizes of receive mesg buffers */
    for(k=0; k<nRecvMsgs3; k++)
    {
      recvMem += LC_GetBufferSize(recvMsgs3[k]);
    }

    sprintf(cBuffer,
            "DDD MESG [%03d]: SHOW_MEM "
            "msgs  send=%010ld recv=%010ld all=%010ld\n",
            me, (long)sendMem, (long)recvMem, (long)(sendMem+recvMem));
    DDD_PrintLine(cBuffer);
  }

  /* display information about recv-messages on lowcomm-level */
  if (DDD_GetOption(OPT_INFO_JOIN) & JOIN_SHOW_MSGSALL)
  {
    DDD_SyncAll();
    if (me==master)
      DDD_PrintLine("DDD JOIN_SHOW_MSGSALL: Phase3Msg.Recv\n");
    LC_PrintRecvMsgs();
  }

  /* unpack messages */
  STAT_RESET;
  UnpackPhase3Msgs(recvMsgs3, nRecvMsgs3, arrayJIJoin);
  LC_Cleanup();
  STAT_TIMER(T_JOIN_UNPACK);

  for(; sendMsgs3!=NULL; sendMsgs3=sm3)
  {
    sm3 = sendMsgs3->next;
    FreeTmp(sendMsgs3, 0);
  }






  /*
          free temporary storage
   */
  JIJoinPtrArray_Free(arrayJIJoin);
  JIJoinSet_Reset(joinGlobals.setJIJoin);

  JIAddCplPtrArray_Free(arrayJIAddCpl2);
  JIAddCplSet_Reset(joinGlobals.setJIAddCpl2);

  JIAddCplPtrArray_Free(arrayJIAddCpl3);
  JIAddCplSet_Reset(joinGlobals.setJIAddCpl3);

  if (localCplObjs!=NULL) FreeTmp(localCplObjs, 0);

  if (joinObjs!=NULL) FreeTmp(joinObjs, 0);

  for(; sendMsgs1!=NULL; sendMsgs1=sm1)
  {
    sm1 = sendMsgs1->next;
    FreeTmp(sendMsgs1, 0);
  }



        #ifdef JoinMemFromHeap
  ReleaseHeap();
  LC_SetMemMgr(memmgr_AllocTMEM, memmgr_FreeTMEM,
               memmgr_AllocTMEM, memmgr_FreeTMEM);
        #endif

#       if DebugJoin<=4
  sprintf(cBuffer,"%4d: JoinEnd, before IFAllFromScratch().\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  /* re-create all interfaces and step JMODE */
  STAT_RESET;
  IFAllFromScratch();
  STAT_TIMER(T_JOIN_BUILD_IF);


  JoinStepMode(JMODE_BUSY);

  return(DDD_RET_OK);
}





/****************************************************************************/
/*                                                                          */
/* Function:  DDD_JoinObj                                                   */
/*                                                                          */
/****************************************************************************/

/**
        Join local object with a distributed object.

        \todoTBC

   @param hdr  DDD local object which should be joined.
 */



void DDD_JoinObj (DDD_HDR hdr, DDD_PROC dest, DDD_GID new_gid)
{
  JIJoin *ji;

  if (!ddd_JoinActive())
  {
    DDD_PrintError('E', 7012, "Missing DDD_JoinBegin(). aborted");
    HARD_EXIT;
  }

  if (dest>=procs)
  {
    sprintf(cBuffer, "cannot join %08x with %08x on processor %d (procs=%d)",
            OBJ_GID(hdr), new_gid, dest, procs);
    DDD_PrintError('E', 7003, cBuffer);
    HARD_EXIT;
  }

  if (dest==me)
  {
    sprintf(cBuffer, "cannot join %08x with myself", OBJ_GID(hdr));
    DDD_PrintError('E', 7004, cBuffer);
    HARD_EXIT;
  }

  if (ObjHasCpl(hdr))
  {
    sprintf(cBuffer, "cannot join %08x, object already distributed",
            OBJ_GID(hdr));
    DDD_PrintError('E', 7005, cBuffer);
    HARD_EXIT;
  }



  ji = JIJoinSet_NewItem(joinGlobals.setJIJoin);
  ji->hdr     = hdr;
  ji->dest    = dest;
  ji->new_gid = new_gid;

  if (! JIJoinSet_ItemOK(joinGlobals.setJIJoin))
    return;

#       if DebugJoin<=2
  sprintf(cBuffer, "%4d: DDD_JoinObj %08x, dest=%d, new_gid=%08x\n",
          me, OBJ_GID(hdr), dest, new_gid);
  DDD_PrintDebug(cBuffer);
#       endif
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_JoinBegin                                                 */
/*                                                                          */
/****************************************************************************/

/**
        Starts join phase.
        A call to this function establishes a global join operation.
        It must be issued on all processors. After this call an arbitrary
        series of {\bf Join}-commands may be issued. The global transfer operation
        is carried out via a \funk{JoinEnd} call on each processor.
 */

#ifdef C_FRONTEND
void DDD_JoinBegin (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::JoinBegin (void)
#endif
#ifdef F_FRONTEND
void DDD_JoinBegin (void)
#endif
{
  /* step mode and check whether call to JoinBegin is valid */
  if (!JoinStepMode(JMODE_IDLE))
  {
    DDD_PrintError('E', 7010, "DDD_JoinBegin() aborted");
    HARD_EXIT;
  }


  /* set kind of TMEM alloc/free requests */
  join_SetTmpMem(TMEM_JOIN);
}


/****************************************************************************/

#undef _FADR

/****************************************************************************/
