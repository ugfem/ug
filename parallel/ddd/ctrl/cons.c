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
#include <assert.h>

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

/*
   PAIRS:     check existance of object for each coupling
   ALLTOALL:  check if all coupling lists are equal
 */
#define CHECK_CPL_PAIRS
#define CHECK_CPL_ALLTOALL



/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

typedef struct
{
  DDD_GID gid;
  DDD_TYPE typ;
  DDD_PROC dest;
  DDD_PROC proc;
  DDD_PRIO prio;
} CONS_INFO;




typedef struct _CONSMSG
{
  DDD_PROC dest;

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
  constab_id = LC_NewMsgTable("ConsTab", consmsg_t, sizeof(CONS_INFO));
}


void ddd_ConsExit (void)
{}


/****************************************************************************/

static int ConsBuildMsgInfos (CONS_INFO *allItems, int nXferItems, CONSMSG **theMsgs)
{
  CONSMSG    *cm, *lastCm;
  int i, lastdest, nMsgs;

  lastdest = -1;
  lastCm = cm = NULL;
  nMsgs = 0;
  for(i=0; i<nXferItems; i++)
  {
    /* eventually create new message item */
    if (allItems[i].dest != lastdest)
    {
      cm = (CONSMSG *) AllocTmp(sizeof(CONSMSG));
      if (cm==NULL)
      {
        DDD_PrintError('E', 9900, STR_NOMEM " in ConsBuildMsgInfos");
        return(0);
      }

      cm->nItems     = 0;
      cm->consArray = &allItems[i];
      cm->dest = allItems[i].dest;
      cm->next = lastCm;
      lastCm = cm;
      lastdest = cm->dest;
      nMsgs++;
    }
    cm->nItems++;
  }
  *theMsgs = cm;


  /* initiate send messages */
  for(cm=*theMsgs; cm!=NULL; cm=cm->next)
  {
    /* create new send message */
    cm->msg_h = LC_NewSendMsg(consmsg_t, cm->dest);

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



static int ConsCheckSingleMsg (LC_MSGHANDLE xm, DDD_HDR *locObjs)
{
  CONS_INFO    *theCplBuf;
  int i, j, nItems;
  int error_cnt = 0;


  nItems = (int) LC_GetTableLen(xm,constab_id);
  theCplBuf = (CONS_INFO *) LC_GetPtr(xm, constab_id);


  /*
          sprintf(cBuffer, "%4d: checking message from proc %d (%d items)\n",
                  me, LC_MsgGetProc(xm), nItems);
          DDD_PrintDebug(cBuffer);
   */


  /* test whether there are consistent objects for all couplings */
  for(i=0, j=0; i<nItems; i++)
  {
    while ((j<ddd_nObjs) && (OBJ_GID(locObjs[j]) < theCplBuf[i].gid))
      j++;

    if ((j<ddd_nObjs) && (OBJ_GID(locObjs[j])==theCplBuf[i].gid))
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

  return(error_cnt);
}



static int sort_CplBufDest (const void *e1, const void *e2)
{
  CONS_INFO   *ci1, *ci2;

  ci1 = (CONS_INFO *)e1;
  ci2 = (CONS_INFO *)e2;

  if (ci1->dest < ci2->dest) return(-1);
  if (ci1->dest > ci2->dest) return(1);

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
  DDD_HDR      *locObjs = NULL;


  /* count overall number of couplings */
  for(i=0, lenCplBuf=0; i<NCPL_GET; i++)
    lenCplBuf += IdxNCpl(i);

  /* get storage for messages */
  cplBuf = (CONS_INFO *) AllocTmp(lenCplBuf*sizeof(CONS_INFO));

  /* copy CONS_INFOs into message buffer */
  for(i=0, j=0; i<NCPL_GET; i++)
  {
    for(cpl=IdxCplList(i); cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      if ((DDD_PROC)cpl->proc >= procs)
      {
        error_cnt++;
        sprintf(cBuffer, "%4d: DDD-GCC Warning: invalid proc=%d (%08x/%08x)\n",
                me, cpl->proc, OBJ_GID(cpl->obj),
                OBJ_GID(ddd_ObjTable[i])
                );
        DDD_PrintLine(cBuffer);
      }
      cplBuf[j].gid  = OBJ_GID(cpl->obj);
      cplBuf[j].typ  = OBJ_TYPE(cpl->obj);
      cplBuf[j].dest = cpl->proc;
      cplBuf[j].proc = cpl->proc;
      cplBuf[j].prio = cpl->prio;
      j++;
    }
  }
  assert(j==lenCplBuf);

  /* sort couplings */
  qsort(cplBuf, lenCplBuf, sizeof(CONS_INFO), sort_CplBufDest);

  /* accumulate messages (one for each partner); inform receivers */
  nSendMsgs = ConsBuildMsgInfos(cplBuf, lenCplBuf, &sendMsgs);

  /* init communication topology */
  nRecvMsgs = LC_Connect(consmsg_t);

  /* build and send messages */
  ConsSend(sendMsgs);


  /* communicate set of messages (send AND receive) */
  recvMsgs = LC_Communicate();


  /* perform checking of received data */
  if (nRecvMsgs>0) locObjs = LocalObjectsList();
  for(i=0; i<nRecvMsgs; i++)
  {
    error_cnt += ConsCheckSingleMsg(recvMsgs[i], locObjs);
  }
  if (nRecvMsgs>0) FreeTmp(locObjs);



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




/****************************************************************************/


static int Cons2CheckSingleMsg (LC_MSGHANDLE xm, DDD_HDR *locObjs)
{
  CONS_INFO    *theCplBuf;
  int i, inext=0, j, nItems;
  int error_cnt = 0;


  nItems = (int) LC_GetTableLen(xm,constab_id);
  theCplBuf = (CONS_INFO *) LC_GetPtr(xm, constab_id);


  /*
          sprintf(cBuffer, "%4d: checking message from proc %d (%d items)\n",
                  me, LC_MsgGetProc(xm), nItems);
          DDD_PrintDebug(cBuffer);
   */


  /* test whether there are consistent objects for all couplings */
  for(i=0, j=0; i<nItems; i=inext)
  {
    inext = i+1;

    while ((j<ddd_nObjs) && (OBJ_GID(locObjs[j]) < theCplBuf[i].gid))
      j++;

    if ((j<ddd_nObjs) && (OBJ_GID(locObjs[j])==theCplBuf[i].gid))
    {
      if (theCplBuf[i].proc == me)
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
        int i2;
        COUPLING *j2;

        for(j2=ObjCplList(locObjs[j]); j2!=NULL; j2=CPL_NEXT(j2))
        {
          int ifound = -1;

          for(i2=i; i2<nItems && theCplBuf[i2].gid==theCplBuf[i].gid;
              i2++)
          {
            if (theCplBuf[i2].proc==j2->proc)
            {
              ifound = i2;
              break;
            }
          }

          if (ifound==-1)
          {
            sprintf(cBuffer, "    DDD-GCC Warning: obj %08x type %d on %d has cpl"
                    " from%4d, but %d hasn't!\n",
                    theCplBuf[i].gid, theCplBuf[i].typ, me,
                    j2->proc, LC_MsgGetProc(xm));
            DDD_PrintLine(cBuffer);

            error_cnt++;
          }
        }

        for(; inext<nItems && theCplBuf[inext].gid==theCplBuf[i].gid; inext++)
          ;
      }
    }
    /*   this is superfluous and (more important:) wrong!
                    else
                    {
                            sprintf(cBuffer, " X  DDD-GCC Warning: obj %08x type %d on %d for cpl"
                                    " from %3d missing!\n",
                                    theCplBuf[i].gid, theCplBuf[i].typ, me, LC_MsgGetProc(xm));
                            DDD_PrintLine(cBuffer);

                            error_cnt++;
                    }
     */
  }

  /* the next too lines are wrong, bug find by PURIFY. KB 970416.
     the locObjs list is freed in Cons2CheckGlobalCpl!
     if (locObjs!=NULL)
          FreeTmp(locObjs);
   */


  return(error_cnt);
}



static int Cons2CheckGlobalCpl (void)
{
  CONS_INFO *cplBuf;
  COUPLING     *cpl, *cpl2;
  int i, j, lenCplBuf, nRecvMsgs, nSendMsgs;
  CONSMSG      *sendMsgs, *cm=0;
  LC_MSGHANDLE *recvMsgs;
  int error_cnt = 0;
  DDD_HDR      *locObjs = NULL;

  /* count overall number of couplings */
  for(i=0, lenCplBuf=0; i<NCPL_GET; i++)
    lenCplBuf += (IdxNCpl(i) * (IdxNCpl(i)+1));

  /* get storage for messages */
  cplBuf = (CONS_INFO *) AllocTmp(lenCplBuf*sizeof(CONS_INFO));

  /* copy CONS_INFOs into message buffer */
  for(i=0, j=0; i<NCPL_GET; i++)
  {
    for(cpl=IdxCplList(i); cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      cplBuf[j].gid  = OBJ_GID(cpl->obj);
      cplBuf[j].typ  = OBJ_TYPE(cpl->obj);
      cplBuf[j].dest = cpl->proc;
      cplBuf[j].proc = me;
      cplBuf[j].prio = OBJ_PRIO(cpl->obj);
      j++;

      for(cpl2=IdxCplList(i); cpl2!=NULL; cpl2=CPL_NEXT(cpl2))
      {
        cplBuf[j].gid  = OBJ_GID(cpl->obj);
        cplBuf[j].typ  = OBJ_TYPE(cpl->obj);
        cplBuf[j].dest = cpl->proc;
        cplBuf[j].proc = cpl2->proc;
        cplBuf[j].prio = cpl2->prio;
        j++;
      }
    }
  }
  assert(j==lenCplBuf);

  /* sort couplings */
  qsort(cplBuf, lenCplBuf, sizeof(CONS_INFO), sort_CplBufDest);

  /* accumulate messages (one for each partner); inform receivers */
  nSendMsgs = ConsBuildMsgInfos(cplBuf, lenCplBuf, &sendMsgs);

  /* init communication topology */
  nRecvMsgs = LC_Connect(consmsg_t);

  /* build and send messages */
  ConsSend(sendMsgs);

  /* communicate set of messages (send AND receive) */
  recvMsgs = LC_Communicate();


  /* perform checking of received data */
  if (nRecvMsgs>0) locObjs = LocalObjectsList();
  for(i=0; i<nRecvMsgs; i++)
  {
    error_cnt += Cons2CheckSingleMsg(recvMsgs[i], locObjs);
  }
  if (nRecvMsgs>0) FreeTmp(locObjs);


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



/****************************************************************************/

static int ConsCheckDoubleObj (void)
{
  DDD_HDR      *locObjs;
  int i, error_cnt = 0;

  locObjs = LocalObjectsList();

  for(i=1; i<ddd_nObjs; i++)
  {
    if (OBJ_GID(locObjs[i-1])==OBJ_GID(locObjs[i]))
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
/****************************************************************************/

/**
        Check DDD runtime consistency.
        This function performs a combined local/global consistency
        check on the object data structures and interfaces managed by DDD.
        This may be used for debugging purposes; if errors are detected,
        then some understanding of internal DDD structures will be useful.

        The following single aspects will be checked:
        \begin{itemize}
        \item double existence of {\em global ID} numbers in each processor's
                  set of local objects.
        \item consistency of coupling lists and object copies
        \item non-symmetric interfaces between processor pairs
        \item non-symmetric number of items in each interface
        \end{itemize}

   @returns  total number of errors (sum of all procs)
 */

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
#ifdef CHECK_CPL_PAIRS
  total_errors += ConsCheckGlobalCpl();
#endif
#ifdef CHECK_CPL_ALLTOALL
  total_errors += Cons2CheckGlobalCpl();
#endif
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
