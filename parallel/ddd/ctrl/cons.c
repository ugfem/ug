// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      cons.c                                                        */
/*                                                                          */
/* Purpose:   consistency checker for ddd structures                        */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/03/21 kb  begin                                            */
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
/*#include "xfer/xfer.h"*/
#include "basic/lowcomm.h"


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
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

typedef struct
{
  DDD_GID gid;
  DDD_TYPE typ;
  DDD_PROC proc;
  DDD_PRIO prio;
} CONS_INFO;



typedef struct _CONSMSG
{
  DDD_PROC proc;

  struct _CONSMSG *next;


  CONS_INFO *consArray;
  int nItems;

  /* lowcomm message handle */
  LC_MSGHANDLE msg_h;

} CONSMSG;



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


static LC_MSGTYPE consmsg_t;
static LC_MSGCOMP constab_id;



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


void ddd_ConsInit (void)
{
  consmsg_t = LC_NewMsgType("ConsCheckMsg");
  constab_id = LC_NewMsgTable(consmsg_t, sizeof(CONS_INFO));
}


void ddd_ConsExit (void)
{}




static int ConsBuildMsgInfos (CONS_INFO *allItems, int nXferItems, CONSMSG **theMsgs)
{
  CONSMSG    *cm, *lastCm;
  int i, lastproc, nMsgs;

  lastproc = -1;
  lastCm = cm = NULL;
  nMsgs = 0;
  for(i=0; i<nXferItems; i++)
  {
    /* eventually create new message item */
    if (allItems[i].proc != lastproc)
    {
      cm = (CONSMSG *) AllocTmp(sizeof(CONSMSG));
      if (cm==NULL)
      {
        DDD_PrintError('E', 9900,
                       "not enough memory in ConsBuildMsgInfos");
        return(0);
      }

      cm->nItems     = 0;
      cm->consArray = &allItems[i];
      cm->proc = allItems[i].proc;
      cm->next = lastCm;
      lastCm = cm;
      lastproc = cm->proc;
      nMsgs++;
    }
    cm->nItems++;
  }
  *theMsgs = cm;


  /* initiate send messages */
  for(cm=*theMsgs; cm!=NULL; cm=cm->next)
  {
    /* create new send message */
    cm->msg_h = LC_NewSendMsg(consmsg_t, cm->proc);

    /* init table inside message */
    LC_SetTableSize(cm->msg_h, constab_id, cm->nItems);

    /* prepare message for sending away */
    LC_MsgPrepareSend(cm->msg_h);
  }

  return(nMsgs);
}



static void ConsSend (CONSMSG *theMsgs)
{
  CONSMSG *cm;

  for(cm=theMsgs; cm!=NULL; cm=cm->next)
  {
    /* copy data into message */
    memcpy(LC_GetPtr(cm->msg_h, constab_id),
           cm->consArray, sizeof(CONS_INFO)*cm->nItems);

    /* send message */
    LC_MsgSend(cm->msg_h);
  }
}



static int ConsCheckSingleMsg (LC_MSGHANDLE xm)
{
  DDD_HDR      *locObjs;
  CONS_INFO    *theCplBuf;
  int i, j, nItems;
  int error_cnt = 0;


  nItems = LC_GetTableLen(xm,constab_id);
  locObjs = LocalObjectsList();
  theCplBuf = (CONS_INFO *) LC_GetPtr(xm, constab_id);


  /*
          sprintf(cBuffer, "%4d: checking message from proc %d (%d items)\n",
                  me, LC_MsgGetProc(xm), nItems);
          DDD_PrintDebug(cBuffer);
   */


  /* test whether there are consistent objects for all couplings */
  for(i=0, j=0; i<nItems; i++)
  {
    while ((j<nObjs) && (OBJ_GID(locObjs[j]) < theCplBuf[i].gid))
      j++;

    if ((j<nObjs) && (OBJ_GID(locObjs[j])==theCplBuf[i].gid))
    {
      if (OBJ_PRIO(locObjs[j])!=theCplBuf[i].prio)
      {
        sprintf(cBuffer, "    DDD-GCC Warning: obj %08x type %d on %d"
                " has prio %d, cpl from %d has prio %d!\n",
                OBJ_GID(locObjs[j]), OBJ_TYPE(locObjs[j]), me, OBJ_PRIO(locObjs[j]),
                LC_MsgGetProc(xm), theCplBuf[i].prio);
        DDD_PrintLine(cBuffer);

        error_cnt++;
      }
    }
    else
    {
      sprintf(cBuffer, "    DDD-GCC Warning: obj %08x type %d on %d for cpl"
              " from %3d missing!\n",
              theCplBuf[i].gid, theCplBuf[i].typ, me, LC_MsgGetProc(xm));
      DDD_PrintLine(cBuffer);

      error_cnt++;
    }
  }

  if (locObjs!=NULL)
    FreeTmp(locObjs);

  return(error_cnt);
}



static int sort_CplBufProc (const void *e1, const void *e2)
{
  CONS_INFO   *ci1, *ci2;

  ci1 = (CONS_INFO *)e1;
  ci2 = (CONS_INFO *)e2;

  if (ci1->proc < ci2->proc) return(-1);
  if (ci1->proc > ci2->proc) return(1);

  if (ci1->gid < ci2->gid) return(-1);
  if (ci1->gid > ci2->gid) return(1);

  return(0);
}


static int ConsCheckGlobalCpl (void)
{
  CONS_INFO *cplBuf;
  COUPLING     *cpl;
  int i, j, lenCplBuf, nRecvMsgs, nSendMsgs;
  CONSMSG      *sendMsgs, *cm=0;
  LC_MSGHANDLE *recvMsgs;
  int error_cnt = 0;


  /* count overall number of couplings */
  for(i=0, lenCplBuf=0; i<nCpls; i++)
    lenCplBuf += theCplN[i];

  /* get storage for messages */
  cplBuf = (CONS_INFO *) AllocTmp(lenCplBuf*sizeof(CONS_INFO));

  /* copy CONS_INFOs into message buffer */
  for(i=0, j=0; i<nCpls; i++)
  {
    for(cpl=theCpl[i]; cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      if ((DDD_PROC)cpl->proc >= procs)
      {
        error_cnt++;
        sprintf(cBuffer, "%4d: DDD-GCC Warning: invalid proc=%d (%08x/%08x)\n",
                me, cpl->proc, OBJ_GID(cpl->obj), OBJ_GID(theObj[i]));
        DDD_PrintLine(cBuffer);
      }
      cplBuf[j].gid  = OBJ_GID(cpl->obj);
      cplBuf[j].typ  = OBJ_TYPE(cpl->obj);
      cplBuf[j].proc = cpl->proc;
      cplBuf[j].prio = cpl->prio;
      j++;
    }
  }

  /* sort couplings */
  qsort(cplBuf, lenCplBuf, sizeof(CONS_INFO), sort_CplBufProc);

  /* accumulate messages (one for each partner); inform receivers */
  nSendMsgs = ConsBuildMsgInfos(cplBuf, lenCplBuf, &sendMsgs);

  /* init communication topology */
  nRecvMsgs = LC_Connect(consmsg_t);

  /* build and send messages */
  ConsSend(sendMsgs);


  /* communicate set of messages (send AND receive) */
  recvMsgs = LC_Communicate();


  /* perform checking of received data */
  for(i=0; i<nRecvMsgs; i++)
  {
    error_cnt += ConsCheckSingleMsg(recvMsgs[i]);
  }


  /* cleanup low-comm layer */
  LC_Cleanup();


  /* free temporary storage */
  if (cplBuf!=NULL)
    FreeTmp(cplBuf);

  for(; sendMsgs!=NULL; sendMsgs=cm)
  {
    cm = sendMsgs->next;
    FreeTmp(sendMsgs);
  }

  return(error_cnt);
}



static int ConsCheckDoubleObj (void)
{
  DDD_HDR      *locObjs;
  int i, error_cnt = 0;

  locObjs = LocalObjectsList();

  for(i=1; i<nObjs; i++)
  {
    if (OBJ_GID(theObj[i-1])==OBJ_GID(theObj[i]))
    {
      error_cnt++;
      sprintf(cBuffer, "    DDD-GCC Warning: obj %08x on %d doubled\n",
              OBJ_GID(locObjs[i]), me);
      DDD_PrintLine(cBuffer);
    }
  }

  if (locObjs!=NULL)
    FreeTmp(locObjs);

  return(error_cnt);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_ConsCheck                                                 */
/*                                                                          */
/* Purpose:   check consistency of ddd structures                           */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    total number of errors (sum of all procs)                     */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND) || defined(F_FRONTEND)
int DDD_ConsCheck (void)
#endif
#ifdef CPP_FRONTEND
int DDD_Library::ConsCheck (void)
#endif
{
  int total_errors=0;

  DDD_Flush();
  Synchronize();
  if (DDD_GetOption(OPT_QUIET_CONSCHECK)==OPT_OFF)
  {
    if (me==master)
      DDD_PrintLine("   DDD-GCC (Global Consistency Check)\n");
  }

  total_errors += ConsCheckDoubleObj();
  total_errors += ConsCheckGlobalCpl();
  total_errors += DDD_CheckInterfaces();


  /* compute sum of errors over all processors */
  total_errors = ddd_GlobalSumInt(total_errors);

  DDD_Flush();
  Synchronize();
  if (DDD_GetOption(OPT_QUIET_CONSCHECK)==OPT_OFF)
  {
    if (me==master)
    {
      sprintf(cBuffer, "   DDD-GCC ready (%d errors)\n", total_errors);
      DDD_PrintLine(cBuffer);
    }
  }

  return(total_errors);
}
