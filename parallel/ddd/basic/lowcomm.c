// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      lowcomm.c                                                     */
/*                                                                          */
/* Purpose:   lowlevel communication layer                                  */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   960715 kb  begin                                              */
/*                                                                          */
/*                                                                          */
/*                                                                          */
/* NOTE: THIS MODULE IS NOT IN ITS FINAL STATE. THE INTERFACE DESIGN, IMPL. */
/*       AND DOCUMENTATION NEED SOME REWORKING!!                            */
/*                                                                          */
/*                                                                          */
/* Remarks:                                                                 */
/*            This modules implements a MSG class, which is a hierarchical  */
/*            format for xfer messages. Each MSG is composed of a number of */
/*            Chunks of data, and has the following format:                 */
/*                                                                          */
/*               description                           |  type              */
/*              ---------------------------------------+---------           */
/*               magic number                          |  ULONG             */
/*               #chunks                               |  ULONG             */
/*               offset chunk1 (from beginning of Msg) |  ULONG             */
/*               length chunk1 (in bytes)              |  ULONG             */
/*               nItems chunk1                         |  ULONG             */
/*                 ...                                 |  ...               */
/*               offset chunkN                         |  ULONG             */
/*               length chunkN                         |  ULONG             */
/*               nItems chunkN                         |  ULONG             */
/*               chunk1                                                     */
/*                ...                                                       */
/*               chunkN                                                     */
/*                                                                          */
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
#include "basic/lowcomm.h"
#include "basic/notify.h"


#define DebugLowComm  10  /* 0 is all, 10 is off */


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


/* maximum number of components in a message */
#define MAX_COMPONENTS 8


/* dummy magic number for messages */
#define MAGIC_DUMMY 0x1234


/* number of entries per chunk in message header */
#define HDR_ENTRIES_PER_CHUNK   3

enum CompType {
  CT_NONE,
  CT_TABLE,
  CT_CHUNK
};


/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

typedef struct _COMP_DESC
{
  int type;                    /* type of this component */
  size_t entry_size;           /* size per entry (for tables) */
} COMP_DESC;


typedef struct _MSG_TYPE
{
  char       *name;                      /* textual description of msgtype */
  int nComps;                            /* number of components */
  COMP_DESC comp[MAX_COMPONENTS];        /* component array */

  struct _MSG_TYPE *next;         /* linked list of all message types */
} MSG_TYPE;


typedef struct _CHUNK_DESC
{
  size_t size;           /* size of chunk (in bytes) */
  ULONG entries;         /* number of valid entries (for tables) */
  ULONG offset;          /* offset of chunk in MSG */
} CHUNK_DESC;


typedef struct _MSG_DESC
{
  MSG_TYPE   *msgType;           /* message type of this actual message */

  ULONG magic;                 /* magic number */
  CHUNK_DESC *chunks;          /* array of chunks */


  size_t bufferSize;           /* size of message buffer (in bytes) */
  char    *buffer;             /* address of message buffer */


  struct _MSG_DESC *next;       /* linked list inside Send/Recv-queue */
  DDD_PROC proc;
  msgid msgId;

} MSG_DESC;


typedef struct _TABLE_DESC
{
  size_t entry_size;             /* size of one table entry */
  int nMax;                      /* number of entries in table */
  int nValid;                    /* number of valid entries */
} TABLE_DESC;


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



static MSG_TYPE *LC_MsgTypes;
static MSG_DESC *LC_SendQueue, *LC_RecvQueue;
static int nSends, nRecvs;
static char *theRecvBuffer;
static LC_MSGHANDLE *theRecvArray;
static MSG_DESC *LC_FreeMsgDescs;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



void LowCommInit (void)
{
  LC_SendQueue = NULL;
  LC_RecvQueue = NULL;
  nSends = 0;
  nRecvs = 0;

  LC_MsgTypes = NULL;

  theRecvArray = NULL;

  LC_FreeMsgDescs = NULL;
}


void LowCommExit (void)
{
  /* TODO: free temporary data */
}


/****************************************************************************/

static MSG_DESC *NewMsgDesc (void)
{
  MSG_DESC *md;

  if (LC_FreeMsgDescs!=NULL)
  {
    /* get item from freelist */
    md = LC_FreeMsgDescs;
    LC_FreeMsgDescs = LC_FreeMsgDescs->next;
  }
  else
  {
    /* freelist is empty */
    md = (MSG_DESC *) AllocCom(sizeof(MSG_DESC));
  }

  return(md);
}


static void FreeMsgDesc (MSG_DESC *md)
{
  /* sort into freelist */
  md->next = LC_FreeMsgDescs;
  LC_FreeMsgDescs = md;
}


/****************************************************************************/




LC_MSGHANDLE LC_NewSendMsg (LC_MSGTYPE mt, DDD_PROC dest)
{
  MSG_TYPE *mtyp = (MSG_TYPE *)mt;
  MSG_DESC *msg = NewMsgDesc();;

#       if DebugLowComm<=6
  sprintf(cBuffer, "%4d: LC_NewSendMsg(%s) dest=%d nSends=%d\n",
          me, mtyp->name, dest, nSends+1);
  DDD_PrintDebug(cBuffer);
#       endif


  msg->msgType = mtyp;
  msg->proc = dest;
  msg->bufferSize = 0;

  /* allocate chunks array */
  msg->chunks = (CHUNK_DESC *) AllocTmp(sizeof(CHUNK_DESC)*mtyp->nComps);

  /* enter message into send queue */
  msg->next = LC_SendQueue;
  LC_SendQueue = msg;
  nSends++;

  return((LC_MSGHANDLE) msg);
}



static LC_MSGHANDLE LC_NewRecvMsg (LC_MSGTYPE mt, DDD_PROC source, size_t size)
{
  MSG_TYPE *mtyp = (MSG_TYPE *)mt;
  MSG_DESC *msg = NewMsgDesc();

#       if DebugLowComm<=6
  sprintf(cBuffer, "%4d: LC_NewRecvMsg(%s) source=%d\n",
          me, mtyp->name, source);
  DDD_PrintDebug(cBuffer);
#       endif

  msg->msgType = mtyp;
  msg->proc = source;
  msg->bufferSize = size;

  /* allocate chunks array */
  msg->chunks = (CHUNK_DESC *) AllocTmp(sizeof(CHUNK_DESC)*mtyp->nComps);

  /* enter message into send queue */
  msg->next = LC_RecvQueue;
  LC_RecvQueue = msg;

  return((LC_MSGHANDLE) msg);
}




static void LC_DeleteMsg (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *)msg;

  FreeTmp(md->chunks);
  FreeMsgDesc(md);
}


static void LC_DeleteMsgBuffer (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *)msg;

  FreeMsg(md->buffer);
}



void LC_SetChunkSize (LC_MSGHANDLE msg, LC_MSGCOMP id, size_t size)
{
  MSG_DESC *md = (MSG_DESC *) msg;

  md->chunks[id].size = size;
  md->chunks[id].entries = 1;
}


void LC_SetTableSize (LC_MSGHANDLE msg, LC_MSGCOMP id, ULONG entries)
{
  MSG_DESC *md = (MSG_DESC *) msg;

  md->chunks[id].size = entries * md->msgType->comp[id].entry_size;
  md->chunks[id].entries = entries;
}


/* returns size of message buffer */

size_t LC_MsgPrepareSend (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *) msg;
  ULONG      *hdr;
  int i, j, n = md->msgType->nComps;

  /* compute size of header */
  md->bufferSize  = 2 * sizeof(ULONG);
  md->bufferSize += (n * HDR_ENTRIES_PER_CHUNK * sizeof(ULONG));

  /* compute size and offset for each chunk */
  for(i=0; i<n; i++)
  {
    md->chunks[i].offset = md->bufferSize;
    md->bufferSize += md->chunks[i].size;
  }


  /* allocate buffer for messages */
  md->buffer = (char *) AllocMsg(md->bufferSize);
  if (md->buffer==NULL)
  {
    sprintf(cBuffer, STR_NOMEM " in LC_MsgPrepareSend (size=%ld)",
            (unsigned long)md->bufferSize);
    DDD_PrintError('E', 6600, cBuffer);
    HARD_EXIT;
  }


  /* enter control data into message header */
  hdr = (ULONG *)md->buffer;
  j=0;
  hdr[j++] = MAGIC_DUMMY;         /* magic number */
  hdr[j++] = n;

  /* enter chunk descriptions into message header */
  for(i=0; i<n; i++)
  {
    hdr[j++] = md->chunks[i].offset;
    hdr[j++] = md->chunks[i].size;
    hdr[j++] = md->chunks[i].entries;
  }

  return(md->bufferSize);
}




DDD_PROC LC_MsgGetProc (LC_MSGHANDLE msg)
{
  MSG_DESC *md = (MSG_DESC *)msg;

  return md->proc;
}


void *LC_GetPtr (LC_MSGHANDLE msg, LC_MSGCOMP id)
{
  MSG_DESC *md = (MSG_DESC *)msg;

  return ((void *)(((char *)md->buffer) + md->chunks[id].offset));
}


void LC_SetTableLen (LC_MSGHANDLE msg, LC_MSGCOMP id, ULONG n)
{
  MSG_DESC *md = (MSG_DESC *)msg;
  ULONG *hdr = (ULONG *)md->buffer;

  hdr[HDR_ENTRIES_PER_CHUNK*id+4] = n;
  md->chunks[id].entries = n;
}


ULONG LC_GetTableLen (LC_MSGHANDLE msg, LC_MSGCOMP id)
{
  MSG_DESC *md = (MSG_DESC *)msg;

  return((ULONG)md->chunks[id].entries);
}



void LC_MsgSend (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *)msg;
  int error;

  /* initiate asynchronous send */
  md->msgId = SendASync(VCHAN_TO(md->proc),
                        md->buffer, md->bufferSize, &error);
}




static void LC_MsgRecv (MSG_DESC *md)
{
  int i, j;

  /* get message address */
  ULONG    *hdr = (ULONG *)md->buffer;

  /* get number of chunks */
  int n = hdr[1];

  /* magic number is hdr[0] */
  if (hdr[0]!=MAGIC_DUMMY)
  {
    sprintf(cBuffer,
            "invalid magic number for message from %d in LC_MsgRecv()",
            md->proc);
    DDD_PrintError('E', 6680, cBuffer);
    HARD_EXIT;
  }

  /* number of chunks must be consistent with message type */
  if (n!=md->msgType->nComps)
  {
    sprintf(cBuffer,
            "wrong number of chunks (%d!=%d) in msg from %d in LC_MsgRecv()",
            n, md->msgType->nComps, md->proc);
    DDD_PrintError('E', 6681, cBuffer);
    HARD_EXIT;
  }

  /* get chunk descriptions from message header */
  for(j=2, i=0; i<n; i++)
  {
    md->chunks[i].offset  = hdr[j++];
    md->chunks[i].size    = hdr[j++];
    md->chunks[i].entries = hdr[j++];
  }
}


size_t LC_GetBufferSize (LC_MSGHANDLE msg)
{
  MSG_DESC *md = (MSG_DESC *) msg;

  return(md->bufferSize);
}



/****************************************************************************/
/*                                                                          */
/* Function:  LC_PollSend                                                   */
/*                                                                          */
/* Purpose:   polls all message-sends one time and returns remaining        */
/*            outstanding messages. whenever a message-send has been        */
/*            completed, its message buffer is freed.                       */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    remaining outstanding messages                                */
/*                                                                          */
/****************************************************************************/

static int LC_PollSend (void)
{
  MSG_DESC *md, *prev, *next=0;
  int remaining, error;

  remaining = 0;
  prev = NULL;
  for(md=LC_SendQueue; md!=NULL; md=next)
  {
    next = md->next;

    error = InfoASend(VCHAN_TO(md->proc), md->msgId);
    /* TODO complete error handling */
    if (error==1)
    {
      /* free message and buffer */
      LC_DeleteMsg((LC_MSGHANDLE)md);
      LC_DeleteMsgBuffer((LC_MSGHANDLE)md);

      /* remove from SendQueue */
      if (prev==NULL)
        LC_SendQueue = next;
      else
        prev->next = next;
      nSends--;
    }
    else
    {
      /* we keep this message in SendQueue */
      remaining++;
      prev = md;
    }
  }

        #if     DebugLowComm<=3
  sprintf(cBuffer, "%4d: LC_PollSend, %d msgs remaining\n",me,remaining);
  DDD_PrintDebug(cBuffer);
        #endif

  return(remaining);
}



/****************************************************************************/
/*                                                                          */
/* Function:  LC_PollRecv                                                   */
/*                                                                          */
/* Purpose:   polls all message-recvs one time and returns remaining        */
/*            outstanding messages. this function doesn't free the message  */
/*            buffers.                                                      */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    remaining outstanding messages                                */
/*                                                                          */
/****************************************************************************/

static int LC_PollRecv (void)
{
  MSG_DESC *md;
  int remaining, error;

  remaining = 0;
  for(md=LC_RecvQueue; md!=NULL; md=md->next)
  {
    if (md->msgId!=-1)
    {
      error = InfoARecv(VCHAN_TO(md->proc), md->msgId);
      /* TODO complete error handling */
      if (error==1)
      {
        md->msgId=-1;

        LC_MsgRecv(md);
      }
      else
      {
        remaining++;
      }
    }
  }

        #if     DebugLowComm<=3
  sprintf(cBuffer, "%4d: LC_PollRecv, %d msgs remaining\n",me,remaining);
  DDD_PrintDebug(cBuffer);
        #endif

  return(remaining);
}


/****************************************************************************/


/*
        allocation of receive message buffers.
        NOTE: one big memory block is allocated and used for all
              message buffers.
 */
static void LC_PrepareRecv (void)
{
  MSG_DESC *md;
  size_t sumSize;
  char     *buffer;
  int error;

  /* compute sum of message buffer sizes */
  for(sumSize=0, md=LC_RecvQueue; md!=NULL; md=md->next)
    sumSize += md->bufferSize;


  /* allocate buffer for messages */
  theRecvBuffer = (char *) AllocMsg(sumSize);
  if (theRecvBuffer==NULL)
  {
    DDD_PrintError('E', 6610, STR_NOMEM " in LC_PrepareRecv");
    sprintf(cBuffer, "(size of message buffer: %d)", sumSize);
    DDD_PrintError('E', 6610, cBuffer);
    HARD_EXIT;
  }


  /* initiate receive calls */
  buffer=theRecvBuffer;
  for(md=LC_RecvQueue; md!=NULL; md=md->next)
  {
    md->buffer = buffer;
    buffer += md->bufferSize;

    md->msgId = RecvASync(VCHAN_TO(md->proc),
                          md->buffer, md->bufferSize, &error);
  }
}




/****************************************************************************/


int LC_Connect (LC_MSGTYPE mtyp)
{
  DDD_PROC    *partners = DDD_ProcArray();
  NOTIFY_DESC *msgs = DDD_NotifyBegin(nSends);
  MSG_DESC *md;
  int i, p;

  if (nSends<0 || nSends>procs-1)
  {
    sprintf(cBuffer, "cannot send %d messages (must be less than %d)",
            nSends, procs-1);
    DDD_PrintError('E', 6620, cBuffer);
    HARD_EXIT;
  }

#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Connect(%s) nSends=%d ...\n",
          me, ((MSG_TYPE *)mtyp)->name, nSends);
  DDD_PrintDebug(cBuffer);
#       endif



  /* fill notify array */
  for(i=0, p=0, md=LC_SendQueue; md!=NULL; i++, md=md->next)
  {
    msgs[i].proc = md->proc;
    msgs[i].size = md->bufferSize;

    /* enhance list of communication partners (destinations) */
    partners[p++] = md->proc;

  }


  /* inform message receivers */
  nRecvs = DDD_Notify();

  if (nRecvs<0 || nRecvs>procs-1)
  {
    sprintf(cBuffer, "cannot receive %d messages (must be less than %d)",
            nRecvs, procs-1);
    DDD_PrintError('E', 6620, cBuffer);
    HARD_EXIT;
  }



#       if DebugLowComm<=7
  sprintf(cBuffer, "%4d: LC_Connect() nSends=%d nRecvs=%d\n",
          me, nSends, nRecvs);
  DDD_PrintDebug(cBuffer);
#       endif


  /* create array of receive message handles */
  if (nRecvs>0)
  {
    theRecvArray = (LC_MSGHANDLE *)AllocTmp(sizeof(LC_MSGHANDLE)*nRecvs);
  }


  /* create recv messages from notify array */
  for(i=0; i<nRecvs; i++)
  {
    /* create recv message handle and store it in MSGHANDLE array */
    theRecvArray[i] = LC_NewRecvMsg(mtyp, msgs[i].proc, msgs[i].size);

    /* enhance list of communication partners (sources) */
    partners[p++] = msgs[i].proc;

  }


  DDD_NotifyEnd();


  /* get necessary connections to comm-partners */
  if (p>0)
    DDD_GetChannels(nRecvs+nSends);


#       if DebugLowComm<=5
  DDD_DisplayTopo();
#       endif


  if (nRecvs>0)
    LC_PrepareRecv();


#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Connect() ready\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  return(nRecvs);
}


LC_MSGHANDLE *LC_Communicate (void)
{
  int leftSend, leftRecv;

#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Communicate() ...\n", me);
  DDD_PrintDebug(cBuffer);
#       endif


  /* poll asynchronous send and receives */
  leftSend = nSends;
  leftRecv = nRecvs;
  do {
    if (leftRecv>0)
      leftRecv = LC_PollRecv();

    if (leftSend>0)
      leftSend = LC_PollSend();

  } while (leftRecv>0 || leftSend>0);


#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Communicate() ready\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  return(theRecvArray);
}


void LC_Cleanup (void)
{
  MSG_DESC *md, *next=0;

#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Cleanup() ...\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  /* free recv messages */
  for(md=LC_RecvQueue; md!=NULL; md=next)
  {
    next = md->next;

    /* free message */
    LC_DeleteMsg((LC_MSGHANDLE)md);
  }


  if (nRecvs>0)
  {
    FreeMsg(theRecvBuffer);
    theRecvBuffer=NULL;
  }

  if (theRecvArray!=NULL)
  {
    FreeTmp(theRecvArray);
    theRecvArray=NULL;
  }

  LC_RecvQueue = NULL;
  nRecvs = 0;

#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Cleanup() ready\n", me);
  DDD_PrintDebug(cBuffer);
#       endif
}



/****************************************************************************/

/*
                MSG_TYPE definition functions
 */



LC_MSGTYPE LC_NewMsgType (char *msgname)
{
  MSG_TYPE *mt;

  mt = (MSG_TYPE *) AllocCom(sizeof(MSG_TYPE));
  /* TODO error handling */
  mt->name   = msgname;
  mt->nComps = 0;

  /* insert into linked list of message types */
  mt->next = LC_MsgTypes;
  LC_MsgTypes = mt;

  return((LC_MSGTYPE) mt);
}


LC_MSGCOMP LC_NewMsgChunk (LC_MSGTYPE mt)
{
  MSG_TYPE  *mtyp = (MSG_TYPE *)mt;
  LC_MSGCOMP id = mtyp->nComps++;

  if (id>=MAX_COMPONENTS)
  {
    sprintf(cBuffer, "too many message components (max. %d)",
            MAX_COMPONENTS);
    DDD_PrintError('E', 6630, cBuffer);
    HARD_EXIT;
  }

  mtyp->comp[id].type = CT_CHUNK;

  return(id);
}



LC_MSGCOMP LC_NewMsgTable (LC_MSGTYPE mt, size_t size)
{
  MSG_TYPE  *mtyp = (MSG_TYPE *)mt;
  LC_MSGCOMP id = mtyp->nComps++;

  if (id>=MAX_COMPONENTS)
  {
    sprintf(cBuffer, "too many message components (max. %d)",
            MAX_COMPONENTS);
    DDD_PrintError('E', 6631, cBuffer);
    HARD_EXIT;
  }

  mtyp->comp[id].type = CT_TABLE;
  mtyp->comp[id].entry_size = size;

  return(id);
}



/****************************************************************************/
