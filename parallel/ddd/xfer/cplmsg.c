// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      cplmsg.c                                                      */
/*                                                                          */
/* Purpose:   ddd object transfer:                                          */
/*               last messages in order to create coupling consistency.     */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   960712 kb  created                                            */
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


#include "dddi.h"
#include "xfer.h"


#define DebugCplMsg  10  /* 10 is off */



/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

/* CPLMSG: complete description of message on sender side */

typedef struct _CPLMSG
{
  DDD_PROC proc;

  struct _CPLMSG *next;


  XIDelCpl  **xferDelCpl;
  int nDelCpl;

  XIModCpl  **xferModCpl;
  int nModCpl;

  XIAddCpl  **xferAddCpl;
  int nAddCpl;


  /* lowcomm message handle */
  LC_MSGHANDLE msg_h;

} CPLMSG;


/****************************************************************************/
/*                                                                          */
/* variables global to this source file only (static)                       */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



static LC_MSGTYPE cplmsg_t;
static LC_MSGCOMP delcpl_id, modcpl_id, addcpl_id;



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


void CplMsgInit (void)
{
  cplmsg_t = LC_NewMsgType("CplMsg");
  delcpl_id = LC_NewMsgTable(cplmsg_t, sizeof(TEDelCpl));
  modcpl_id = LC_NewMsgTable(cplmsg_t, sizeof(TEModCpl));
  addcpl_id = LC_NewMsgTable(cplmsg_t, sizeof(TEAddCpl));
}


void CplMsgExit (void)
{}


/****************************************************************************/

static CPLMSG *CreateCplMsg (DDD_PROC dest, CPLMSG *lastxm)
{
  CPLMSG *xm;

  xm = (CPLMSG *) AllocTmp(sizeof(CPLMSG));
  if (xm==NULL)
  {
    DDD_PrintError('E', 6400, STR_NOMEM " in PrepareCplMsgs");
    HARD_EXIT;
  }

  xm->nDelCpl    = 0;
  xm->xferDelCpl = NULL;
  xm->nModCpl    = 0;
  xm->xferModCpl = NULL;
  xm->nAddCpl    = 0;
  xm->xferAddCpl = NULL;
  xm->proc = dest;
  xm->next = lastxm;

  return(xm);
}


static int PrepareCplMsgs (
  XIDelCpl **itemsDC, int nDC,
  XIModCpl **itemsMC, int nMC,
  XIAddCpl **itemsAC, int nAC,
  CPLMSG **theMsgs)
{
  CPLMSG    *xm=NULL;
  int iDC, iMC, iAC, nMsgs=0;

#       if DebugCplMsg<=3
  printf("%4d: PrepareCplMsgs, nXIDelCpl=%d nXIModCpl=%d nXIAddCpl=%d\n",
         me, nDC, nMC, nAC);
  fflush(stdout);
#       endif


  /*
          run through all tables simultaneously,
          each time a new proc-nr is encountered in one of these
          tables, create a new CPLMSG item.

          (the lists have been sorted according to proc-nr previously.)
   */

  iDC=0; iMC=0; iAC=0;
  while (iDC<nDC || iMC<nMC || iAC<nAC)
  {
    DDD_PROC pDC = (iDC<nDC) ? itemsDC[iDC]->to   : procs;
    DDD_PROC pMC = (iMC<nMC) ? itemsMC[iMC]->to   : procs;
    DDD_PROC pAC = (iAC<nAC) ? itemsAC[iAC]->to   : procs;

    /* check DelCpl-items */
    if (pDC<=pMC && pDC<=pAC && pDC<procs)
    {
      int i;

      if (xm==NULL || xm->proc!=pDC)
      {
        xm = CreateCplMsg(pDC, xm);
        nMsgs++;
      }

      xm->xferDelCpl = itemsDC+iDC;
      for(i=iDC; i<nDC && itemsDC[i]->to==pDC; i++)
        ;

      xm->nDelCpl = i-iDC;
      iDC = i;
    }

    /* check ModCpl-items */
    if (pMC<=pDC && pMC<=pAC && pMC<procs)
    {
      int i;

      if (xm==NULL || xm->proc!=pMC)
      {
        xm = CreateCplMsg(pMC, xm);
        nMsgs++;
      }

      xm->xferModCpl = itemsMC+iMC;
      for(i=iMC; i<nMC && itemsMC[i]->to==pMC; i++)
        ;

      xm->nModCpl = i-iMC;
      iMC = i;
    }

    /* check AddCpl-items */
    if (pAC<=pDC && pAC<=pMC && pAC<procs)
    {
      int i;

      if (xm==NULL || xm->proc!=pAC)
      {
        xm = CreateCplMsg(pAC, xm);
        nMsgs++;
      }

      xm->xferAddCpl = itemsAC+iAC;
      for(i=iAC; i<nAC && itemsAC[i]->to==pAC; i++)
        ;

      xm->nAddCpl = i-iAC;
      iAC = i;
    }

    if (pDC==procs) iDC = nDC;
    if (pMC==procs) iMC = nMC;
    if (pAC==procs) iAC = nAC;
  }
  *theMsgs = xm;



  /* initiate send messages */
  for(xm=*theMsgs; xm!=NULL; xm=xm->next)
  {
    /* create new send message */
    xm->msg_h = LC_NewSendMsg(cplmsg_t, xm->proc);

    /* init tables inside message */
    LC_SetTableSize(xm->msg_h, delcpl_id, xm->nDelCpl);
    LC_SetTableSize(xm->msg_h, modcpl_id, xm->nModCpl);
    LC_SetTableSize(xm->msg_h, addcpl_id, xm->nAddCpl);

    /* prepare message for sending away */
    LC_MsgPrepareSend(xm->msg_h);
  }

  return(nMsgs);
}


static void CplMsgSend (CPLMSG *theMsgs)
{
  CPLMSG *m;

  for(m=theMsgs; m!=NULL; m=m->next)
  {
    int i;
    TEDelCpl *arrayDC = (TEDelCpl *)LC_GetPtr(m->msg_h, delcpl_id);
    TEModCpl *arrayMC = (TEModCpl *)LC_GetPtr(m->msg_h, modcpl_id);
    TEAddCpl *arrayAC = (TEAddCpl *)LC_GetPtr(m->msg_h, addcpl_id);

    /* copy data into message */
    for(i=0; i<m->nDelCpl; i++)
    {
      arrayDC[i] = m->xferDelCpl[i]->te;
    }
    for(i=0; i<m->nModCpl; i++)
    {
      arrayMC[i] = m->xferModCpl[i]->te;

#                       ifdef SLL_DebugNew
      {
        XIModCpl *mc = m->xferModCpl[i];

        sprintf(cBuffer, "%4d: send modcpl to %d (%08x, %3d)  %s:%d\n",
                me, m->proc,
                mc->te.gid, mc->te.prio, mc->sll_file, mc->sll_line);
        DDD_PrintDebug(cBuffer);
      }
#                       endif
    }
    for(i=0; i<m->nAddCpl; i++)
    {
      arrayAC[i] = m->xferAddCpl[i]->te;
    }

    /* schedule message for send */
    LC_MsgSend(m->msg_h);
  }
}



/****************************************************************************/


static void CplMsgUnpackSingle (LC_MSGHANDLE xm,
                                DDD_HDR *localCplObjs, int nLCO)
{
  TEDelCpl  *theDelCpl;
  TEModCpl  *theModCpl;
  TEAddCpl  *theAddCpl;
  int i, j, nDelCpl, nModCpl, nAddCpl;
  DDD_PROC proc = LC_MsgGetProc(xm);

  /* get number and address of del-items */
  nDelCpl = (int) LC_GetTableLen(xm, delcpl_id);
  nModCpl = (int) LC_GetTableLen(xm, modcpl_id);
  nAddCpl = (int) LC_GetTableLen(xm, addcpl_id);
  theDelCpl = (TEDelCpl *) LC_GetPtr(xm, delcpl_id);
  theModCpl = (TEModCpl *) LC_GetPtr(xm, modcpl_id);
  theAddCpl = (TEAddCpl *) LC_GetPtr(xm, addcpl_id);


  /* modify couplings according to mod-list */
  for(i=0, j=0; i<nModCpl; i++)
  {
    while ((j<nLCO) && (OBJ_GID(localCplObjs[j]) < theModCpl[i].gid))
      j++;

    if ((j<nLCO) && (OBJ_GID(localCplObjs[j])==theModCpl[i].gid))
    {
      ModCoupling(localCplObjs[j], proc, theModCpl[i].prio);
    }
  }


  /* delete couplings according to del-list */
  for(i=0, j=0; i<nDelCpl; i++)
  {
    while ((j<nLCO) && (OBJ_GID(localCplObjs[j]) < theDelCpl[i].gid))
      j++;

    if ((j<nLCO) && (OBJ_GID(localCplObjs[j])==theDelCpl[i].gid))
    {
      DelCoupling(localCplObjs[j], proc);
    }
  }


  /* add couplings according to add-list */
  for(i=0, j=0; i<nAddCpl; i++)
  {
    while ((j<nLCO) && (OBJ_GID(localCplObjs[j]) < theAddCpl[i].gid))
      j++;

    if ((j<nLCO) && (OBJ_GID(localCplObjs[j])==theAddCpl[i].gid))
    {
      AddCoupling(localCplObjs[j], theAddCpl[i].proc, theAddCpl[i].prio);
    }
  }
}


/****************************************************************************/


static void CplMsgDisplay (char *comment, LC_MSGHANDLE xm)
{
  TEDelCpl     *theDelCpl;
  TEModCpl     *theModCpl;
  TEAddCpl     *theAddCpl;
  char buf[30];
  int i, proc = LC_MsgGetProc(xm);
  int lenDelCpl = (int) LC_GetTableLen(xm, delcpl_id);
  int lenModCpl = (int) LC_GetTableLen(xm, modcpl_id);
  int lenAddCpl = (int) LC_GetTableLen(xm, addcpl_id);

  sprintf(buf, " %03d-%s-%03d ", me, comment, proc);

  /* get table addresses inside message */
  theDelCpl = (TEDelCpl *)    LC_GetPtr(xm, delcpl_id);
  theModCpl = (TEModCpl *)    LC_GetPtr(xm, modcpl_id);
  theAddCpl = (TEAddCpl *)    LC_GetPtr(xm, addcpl_id);


  sprintf(cBuffer, "%s 04 DelCpl.size=%05d\n", buf, lenDelCpl);
  DDD_PrintDebug(cBuffer);
  sprintf(cBuffer, "%s 05 ModCpl.size=%05d\n", buf, lenModCpl);
  DDD_PrintDebug(cBuffer);
  sprintf(cBuffer, "%s 06 AddCpl.size=%05d\n", buf, lenAddCpl);
  DDD_PrintDebug(cBuffer);


  for(i=0; i<lenDelCpl; i++)
  {
    sprintf(cBuffer, "%s 14 delcpl %04d - %08x\n",
            buf, i, theDelCpl[i].gid);
    DDD_PrintDebug(cBuffer);
  }

  for(i=0; i<lenModCpl; i++)
  {
    sprintf(cBuffer, "%s 15 modcpl %04d - %08x %3d\n",
            buf, i, theModCpl[i].gid, theModCpl[i].prio);
    DDD_PrintDebug(cBuffer);
  }

  for(i=0; i<lenAddCpl; i++)
  {
    sprintf(cBuffer, "%s 16 addcpl %04d - %08x %4d %3d\n",
            buf, i, theAddCpl[i].gid, theAddCpl[i].proc, theAddCpl[i].prio);
    DDD_PrintDebug(cBuffer);
  }
}


/****************************************************************************/


/*
        ...

        localCplObjs is a sorted list of all objects with coupling lists.
 */

void CommunicateCplMsgs (
  XIDelCpl **itemsDC, int nDC,
  XIModCpl **itemsMC, int nMC,
  XIAddCpl **itemsAC, int nAC,
  DDD_HDR *localCplObjs, int nLCO)
{
  CPLMSG    *sendMsgs, *sm=0;
  LC_MSGHANDLE *recvMsgs;
  int i, nSendMsgs, nRecvMsgs;


  /* accumulate messages (one for each partner) */
  nSendMsgs = PrepareCplMsgs(itemsDC, nDC,
                             itemsMC, nMC,
                             itemsAC, nAC,
                             &sendMsgs);

  /* init communication topology */
  nRecvMsgs = LC_Connect(cplmsg_t);

  /* build and send messages */
  CplMsgSend(sendMsgs);


#if DebugCplMsg>2
  if (DDD_GetOption(OPT_DEBUG_XFERMESGS)==OPT_ON)
#endif
  {
    for(sm=sendMsgs; sm!=NULL; sm=sm->next)
    {
      CplMsgDisplay("CS", sm->msg_h);
    }
  }


  /* communicate set of messages (send AND receive) */
  recvMsgs = LC_Communicate();


  for(i=0; i<nRecvMsgs; i++)
  {
    CplMsgUnpackSingle(recvMsgs[i], localCplObjs, nLCO);

    /*
                    if (DDD_GetOption(OPT_DEBUG_XFERMESGS)==OPT_ON)
                            CplMsgDisplay("CR", recvMsgs[i]);
     */
  }


  /* cleanup low-comm layer */
  LC_Cleanup();



  /* free temporary memory */
  for(; sendMsgs!=NULL; sendMsgs=sm)
  {
    sm = sendMsgs->next;
    FreeTmp(sendMsgs);
  }
}


/****************************************************************************/
