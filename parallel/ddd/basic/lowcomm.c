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
/*            971007 kb  reworked                                           */
/*                                                                          */
/* Remarks:                                                                 */
/*            This module provides two basic abstractions:                  */
/*            - sending of messages without explicit receive calls          */
/*            - message types consisting of a set of components, where      */
/*              components are tables (with entries of equal sizes) and     */
/*              raw data chunks.                                            */
/*                                                                          */
/*            The LowComm subsystem uses the Notify-subsystem in order to   */
/*            tell receiving processors that corresponding send-calls had   */
/*            been issued.                                                  */
/*                                                                          */
/*            The structure of each message is:                             */
/*                                                                          */
/*               description                               |  type          */
/*              -------------------------------------------+---------       */
/*               magic number                              |  ULONG         */
/*               #components                               |  ULONG         */
/*               offset component1 (from beginning of Msg) |  ULONG         */
/*               length component1 (in bytes)              |  ULONG         */
/*               nItems component1                         |  ULONG         */
/*                 ...                                     |  ...           */
/*               offset componentN                         |  ULONG         */
/*               length componentN                         |  ULONG         */
/*               nItems componentN                         |  ULONG         */
/*               component1                                                 */
/*                ...                                                       */
/*               componentN                                                 */
/*                                                                          */
/*            The LowComm subsystem is able to handle low-memory situations,*/
/*            where the available memory is not enough for all send- and    */
/*            receive-buffers. See LC_MsgAlloc for details.                 */
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
  char     *name;              /* textual description of component */
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


enum MsgState { MSTATE_NEW, MSTATE_FREEZED, MSTATE_ALLOCATED, MSTATE_COMM, MSTATE_READY };

typedef struct _MSG_DESC
{
  int msgState;                /* message state of this message (one of MsgState) */
  MSG_TYPE   *msgType;         /* message type of this message */

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

static AllocFunc _DefaultAlloc, _SendAlloc, _RecvAlloc;
static FreeFunc _DefaultFree,  _SendFree,  _RecvFree;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/**
        Initiates LowComm subsystem.
        This function has to be called exactly once in order
        to initialize the LowComm subsystem. After a call to
        this function, the functionality of the LowComm can
        be used.

   @param aAllocFunc memory allocation function used as the default
   @param aFreeFunc  memory free function used as the default
 */

void LC_Init (AllocFunc aAllocFunc, FreeFunc aFreeFunc)
{
  LC_SendQueue = NULL;
  LC_RecvQueue = NULL;
  nSends = 0;
  nRecvs = 0;

  LC_MsgTypes = NULL;

  theRecvArray = NULL;

  LC_FreeMsgDescs = NULL;

  _DefaultAlloc = aAllocFunc;
  _DefaultFree  = aFreeFunc;
  LC_SetMemMgrDefault();
}



/**
        Exits LowComm subsystem.
        This function frees memory allocated by the LowComm subsystem
        and shuts down its communication structures.
 */

void LC_Exit (void)
{
  /* TODO: free temporary data */
}




/**
        Customizing memory management for LowComm subsystem.
        Currently this function supports only alloc/free of
        buffers for messages to be sent.
 */
void LC_SetMemMgrSend (AllocFunc aAllocFunc, FreeFunc aFreeFunc)
{
  _SendAlloc = aAllocFunc;
  _SendFree  = aFreeFunc;
}


/**
        Customizing memory management for LowComm subsystem.
        Currently this function supports only alloc/free of
        buffers for messages to be received.
 */
void LC_SetMemMgrRecv (AllocFunc aAllocFunc, FreeFunc aFreeFunc)
{
  _RecvAlloc = aAllocFunc;
  _RecvFree  = aFreeFunc;
}


/**
        Set memory management for LowComm subsystem to its default state.
        Set alloc/free of message buffers to the functions provided to
        \lcfunk{Init}.
 */
void LC_SetMemMgrDefault (void)
{
  _SendAlloc = _DefaultAlloc;
  _SendFree  = _DefaultFree;
  _RecvAlloc = _DefaultAlloc;
  _RecvFree  = _DefaultFree;
}


/****************************************************************************/

/*
        auxiliary functions
 */

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



/*
        this function has internal access only because LowComm initiates
        asynchronous receive calls itself.
 */

static LC_MSGHANDLE LC_NewRecvMsg (LC_MSGTYPE mt, DDD_PROC source, size_t size)
{
  MSG_TYPE *mtyp = (MSG_TYPE *)mt;
  MSG_DESC *msg = NewMsgDesc();

#       if DebugLowComm<=6
  sprintf(cBuffer, "%4d: LC_NewRecvMsg(%s) source=%d\n",
          me, mtyp->name, source);
  DDD_PrintDebug(cBuffer);
#       endif

  msg->msgState = MSTATE_NEW;
  msg->msgType = mtyp;
  msg->proc = source;
  msg->bufferSize = size;

  /* allocate chunks array */
  msg->chunks = (CHUNK_DESC *) AllocTmpReq(sizeof(CHUNK_DESC)*mtyp->nComps, TMEM_LOWCOMM);

  /* enter message into recv queue */
  msg->next = LC_RecvQueue;
  LC_RecvQueue = msg;

  return((LC_MSGHANDLE) msg);
}




static void LC_DeleteMsg (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *)msg;
  size_t size = sizeof(CHUNK_DESC) * md->msgType->nComps;

  FreeTmpReq(md->chunks,size,TMEM_LOWCOMM);
  FreeMsgDesc(md);
}


static void LC_DeleteMsgBuffer (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *)msg;

  if (_SendFree!=NULL)
    (*_SendFree)(md->buffer);
}




static void LC_MsgRecv (MSG_DESC *md)
{
  int i, j;

  /* get message address */
  ULONG    *hdr = (ULONG *)md->buffer;

  /* get number of chunks */
  int n = (int)(hdr[1]);

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
    md->chunks[i].offset  =          hdr[j++];
    md->chunks[i].size    = (size_t)(hdr[j++]);
    md->chunks[i].entries =          hdr[j++];
  }

        #if     DebugLowComm<=2
  sprintf(cBuffer, "%4d: LC_MsgRecv(). from=%d ready.\n",
          me, md->proc);
  DDD_PrintDebug(cBuffer);
        #endif
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
  MSG_DESC *md;
  int remaining, error;

  remaining = 0;
  for(md=LC_SendQueue; md!=NULL; md=md->next)
  {
    if (md->msgState==MSTATE_COMM)
    {
      error = InfoASend(VCHAN_TO(md->proc), md->msgId);
      if (error==-1)
      {
        sprintf(cBuffer,
                "PPIF's InfoASend() failed for send to proc=%d in LowComm",
                md->proc);
        DDD_PrintError('E', 6640, cBuffer);
        HARD_EXIT;
      }

      if (error==1)
      {
        /* free message buffer */
        LC_DeleteMsgBuffer((LC_MSGHANDLE)md);

        md->msgState=MSTATE_READY;
      }
      else
      {
        /* we keep this message in SendQueue */
        remaining++;
      }
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
    if (md->msgState==MSTATE_COMM)
    {
      error = InfoARecv(VCHAN_TO(md->proc), md->msgId);
      if (error==-1)
      {
        sprintf(cBuffer,
                "PPIF's InfoARecv() failed for recv from proc=%d in LowComm",
                md->proc);
        DDD_PrintError('E', 6641, cBuffer);
        HARD_EXIT;
      }

      if (error==1)
      {
        LC_MsgRecv(md);

        md->msgState=MSTATE_READY;
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
/*                                                                          */
/* Function:  LC_FreeSendQueue                                              */
/*                                                                          */
/****************************************************************************/

static void LC_FreeSendQueue (void)
{
  MSG_DESC *md, *next=NULL;

  for(md=LC_SendQueue; md!=NULL; md=next)
  {
    /* the following assertion is too picky. Freeing of
       message queues should be possible in all msgStates. */
    /*assert(md->msgState==MSTATE_READY);*/

    next = md->next;
    LC_DeleteMsg((LC_MSGHANDLE)md);
  }


  LC_SendQueue = NULL;
  nSends = 0;
}


/****************************************************************************/
/*                                                                          */
/* Function:  LC_FreeRecvQueue                                              */
/*                                                                          */
/****************************************************************************/

static void LC_FreeRecvQueue (void)
{
  MSG_DESC *md, *next=NULL;

  for(md=LC_RecvQueue; md!=NULL; md=next)
  {
    /* the following assertion is too picky. Freeing of
       message queues should be possible in all msgStates. */
    /*assert(md->msgState==MSTATE_READY);*/

    next = md->next;
    LC_DeleteMsg((LC_MSGHANDLE)md);
  }


  LC_RecvQueue = NULL;
  nRecvs = 0;
}



/****************************************************************************/


/* LC_MsgFreeze and LC_MsgAlloc are the two parts of LC_MsgPrepareSend(). */

/* returns size of message buffer */

size_t LC_MsgFreeze (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *) msg;
  int i, n = md->msgType->nComps;

  assert(md->msgState==MSTATE_NEW);

  /* compute size of header */
  md->bufferSize  = 2 * sizeof(ULONG);
  md->bufferSize += (n * HDR_ENTRIES_PER_CHUNK * sizeof(ULONG));

  /* compute size and offset for each chunk */
  for(i=0; i<n; i++)
  {
    md->chunks[i].offset = md->bufferSize;
    md->bufferSize += md->chunks[i].size;
  }

  md->msgState=MSTATE_FREEZED;

  return(md->bufferSize);
}



int LC_MsgAlloc (LC_MSGHANDLE msg)
{
  MSG_DESC   *md = (MSG_DESC *) msg;
  ULONG      *hdr;
  int i, j, n = md->msgType->nComps;
  int remaining=1, give_up = FALSE;

  assert(md->msgState==MSTATE_FREEZED);

  /* the following code tries to allocate the message buffer.
     if this fails, the previously started asynchronous sends are
     polled, in order to free their message buffers. if there are
     no remaining async-sends, we give up. */
  do {
    /* allocate buffer for messages */
    md->buffer = (char *) (*_SendAlloc)(md->bufferSize);
    if (md->buffer==NULL)
    {
      if (remaining==0)
        give_up = TRUE;
      else
      {
#                               if DebugLowComm<=7
        sprintf(cBuffer, "%4d: LC_MsgAlloc(%s) detected low memory.\n",
                me, md->msgType->name);
        DDD_PrintDebug(cBuffer);
#                               endif

        /* couldn't get msg-buffer. try to poll previous messages. */
        /* first, poll receives to avoid communication deadlock. */
        LC_PollRecv();

        /* now, try to poll sends and free their message buffers */
        remaining  = LC_PollSend();

#                               if DebugLowComm<=6
        sprintf(cBuffer,
                "%4d: LC_MsgAlloc(%s) preliminary poll, sends_left=%d.\n",
                me, md->msgType->name, remaining);
        DDD_PrintDebug(cBuffer);
#                               endif
      }
    }
  } while (md->buffer==NULL && !give_up);

  if (give_up)
  {
#               if DebugLowComm<=7
    sprintf(cBuffer, "%4d: LC_MsgAlloc(%s) giving up, no memory.\n",
            me, md->msgType->name);
    DDD_PrintDebug(cBuffer);
#               endif
    return(FALSE);
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

  md->msgState=MSTATE_ALLOCATED;

  return(TRUE);
}




/*
        allocation of receive message buffers.
        NOTE: one big memory block is allocated and used for all
              message buffers.
 */
static RETCODE LC_PrepareRecv (void)
{
  MSG_DESC *md;
  size_t sumSize;
  char     *buffer;
  int error;

  /* compute sum of message buffer sizes */
  for(sumSize=0, md=LC_RecvQueue; md!=NULL; md=md->next)
  {
    assert(md->msgState==MSTATE_NEW);
    sumSize += md->bufferSize;
  }


  /* allocate buffer for messages */
  theRecvBuffer = (char *) (*_RecvAlloc)(sumSize);
  if (theRecvBuffer==NULL)
  {
    DDD_PrintError('E', 6610, STR_NOMEM " in LC_PrepareRecv");
    sprintf(cBuffer, "(size of message buffer: %d)", sumSize);
    DDD_PrintError('E', 6610, cBuffer);
    RET_ON_ERROR;
  }


  /* initiate receive calls */
  buffer=theRecvBuffer;
  for(md=LC_RecvQueue; md!=NULL; md=md->next)
  {
    md->buffer = buffer;
    buffer += md->bufferSize;

    md->msgId = RecvASync(VCHAN_TO(md->proc),
                          md->buffer, md->bufferSize, &error);

    md->msgState=MSTATE_COMM;
  }

  RET_ON_OK;
}


/****************************************************************************/

/*
        MSG_TYPE definition functions
 */


/****************************************************************************/
/*                                                                          */
/* Function:  LC_NewMsgType                                                 */
/*                                                                          */
/****************************************************************************/

/**
        Declares new message-type.
        Before messages may be sent and received with the LowComm
        subsystem, at least one {\em message-type} must be defined by
        a global call to this function. Subsequently, calls to
        \lcfunk{NewMsgTable} and \lcfunk{NewMsgChunk} can be used
        in order to define the structure of the new message-type.

        Each message-type in the LowComm subsystem consists of a set
        of {\em message-components}. Possible message-components are:
        {\em tables} (with entries of equal size) and raw {\em data chunks}.
        The set of message-components has the same structure for
        all messages of the same type, but the number of table entries
        and the size of the data chunks differ from message to message.

   @return identifier of new message-type
   @param aName  name of message-type. This string is used for debugging
        and logging output.
 */

LC_MSGTYPE LC_NewMsgType (char *aName)
{
  MSG_TYPE *mt;

  mt = (MSG_TYPE *) AllocCom(sizeof(MSG_TYPE));
  if (mt==NULL)
  {
    DDD_PrintError('E', 6601, STR_NOMEM " in LC_NewMsgType()");
    HARD_EXIT;
  }

  mt->name   = aName;
  mt->nComps = 0;

  /* insert into linked list of message types */
  mt->next = LC_MsgTypes;
  LC_MsgTypes = mt;

  return((LC_MSGTYPE) mt);
}



/****************************************************************************/
/*                                                                          */
/* Function:  LC_NewMsgChunk                                                */
/*                                                                          */
/****************************************************************************/

/**
        Add data chunk to current set of a message-type's message-components.
        This function is called after a previous call to \lcfunk{NewMsgType}
        in order to add a new message-component to the message-type.
        The component added by this function is a chunk of raw data.
        The size of the chunk is not specified here, use \lcfunk{SetChunkSize}
        for specifying the data chunk size for a given (concrete) message.

        See \lcfunk{NewMsgTable} for adding message-tables, which are
        a different kind of message-component.

   @return           identifier of new message-component
   @param  aName     name of new message component
   @param  aMsgType  previously declared message-type
 */

LC_MSGCOMP LC_NewMsgChunk (char *aName, LC_MSGTYPE aMsgType)
{
  MSG_TYPE  *mtyp = (MSG_TYPE *)aMsgType;
  LC_MSGCOMP id = mtyp->nComps++;

  if (id>=MAX_COMPONENTS)
  {
    sprintf(cBuffer, "too many message components (max. %d)",
            MAX_COMPONENTS);
    DDD_PrintError('E', 6630, cBuffer);
    HARD_EXIT;
  }

  mtyp->comp[id].type = CT_CHUNK;
  mtyp->comp[id].name = aName;

  return(id);
}




/****************************************************************************/
/*                                                                          */
/* Function:  LC_NewMsgTable                                                */
/*                                                                          */
/****************************************************************************/

/**
        Add table to current set of a message-type's message-components.
        This function is called after a previous call to \lcfunk{NewMsgType}
        in order to add a new message-component to the message-type.
        The component added by this function is a table of data, where
        each table entry has the same size.
        The overall size of the whole table is not specified here, but only
        the size for one table entry. Use \lcfunk{SetTableSize} for setting
        the number of reserved table entries in a given (concrete) message;
        use \lcfunk{SetTableLen} in order to specify the number of valid
        entries in a given message.

        See \lcfunk{NewMsgChunk} for adding data chunks, which are
        a different kind of message-component.

   @return           identifier of new message-component
   @param  aName     name of new message component
   @param  aMsgType  previously declared message-type
   @param  aSize     size of each table entry (in byte)
 */

LC_MSGCOMP LC_NewMsgTable (char *aName, LC_MSGTYPE aMsgType, size_t aSize)
{
  MSG_TYPE  *mtyp = (MSG_TYPE *)aMsgType;
  LC_MSGCOMP id = mtyp->nComps++;

  if (id>=MAX_COMPONENTS)
  {
    sprintf(cBuffer, "too many message components (max. %d)",
            MAX_COMPONENTS);
    DDD_PrintError('E', 6631, cBuffer);
    HARD_EXIT;
  }

  mtyp->comp[id].type = CT_TABLE;
  mtyp->comp[id].entry_size = aSize;
  mtyp->comp[id].name = aName;

  return(id);
}



/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* Function:  LC_NewSendMsg                                                 */
/*                                                                          */
/****************************************************************************/

/**
        Create new message on sending processor.
        This function creates a new message handle on the sending processor and
        links it into the LowComm send-queue. The message has a given message-type
        and a given destination processor. Before the message is actually sent
        (by calling \lcfunk{MsgSend}), the sizes of the message's components
        must be set (\lcfunk{SetTableSize}, \lcfunk{SetChunkSize}) and the message
        buffer must be prepared (via \lcfunk{MsgPrepareSend}). After that,
        the message's tables and chunks can be filled with data and the message
        sending process can be initiated by \lcfunk{MsgSend}.

   @return          identifier of new message
   @param aMsgType  message-type for new message
   @param aDest     destination processor of new message
 */

LC_MSGHANDLE LC_NewSendMsg (LC_MSGTYPE aMsgType, DDD_PROC aDest)
{
  MSG_TYPE *mtyp = (MSG_TYPE *)aMsgType;
  MSG_DESC *msg = NewMsgDesc();

#       if DebugLowComm<=6
  sprintf(cBuffer, "%4d: LC_NewSendMsg(%s) dest=%d nSends=%d\n",
          me, mtyp->name, aDest, nSends+1);
  DDD_PrintDebug(cBuffer);
#       endif


  msg->msgState = MSTATE_NEW;
  msg->msgType = mtyp;
  msg->proc = aDest;
  msg->bufferSize = 0;

  /* allocate chunks array */
  msg->chunks = (CHUNK_DESC *) AllocTmpReq(sizeof(CHUNK_DESC)*mtyp->nComps, TMEM_LOWCOMM);
  if (msg->chunks==NULL)
  {
    DDD_PrintError('E', 6602, STR_NOMEM " in LC_NewSendMsg()");
    HARD_EXIT;
  }


  /* enter message into send queue */
  msg->next = LC_SendQueue;
  LC_SendQueue = msg;
  nSends++;

  return((LC_MSGHANDLE) msg);
}




void LC_SetChunkSize (LC_MSGHANDLE msg, LC_MSGCOMP id, size_t size)
{
  MSG_DESC *md = (MSG_DESC *) msg;

  assert(md->msgState==MSTATE_NEW);
  assert(id < md->msgType->nComps);

  md->chunks[id].size = size;
  md->chunks[id].entries = 1;
}


void LC_SetTableSize (LC_MSGHANDLE msg, LC_MSGCOMP id, ULONG entries)
{
  MSG_DESC *md = (MSG_DESC *) msg;

  assert(md->msgState==MSTATE_NEW);
  assert(id < md->msgType->nComps);

  md->chunks[id].size = ((int)entries) * md->msgType->comp[id].entry_size;
  md->chunks[id].entries = entries;
}



/****************************************************************************/
/*                                                                          */
/* Function:  LC_MsgPrepareSend                                             */
/*                                                                          */
/****************************************************************************/


/* returns size of message buffer */

size_t LC_MsgPrepareSend (LC_MSGHANDLE msg)
{
  size_t size = LC_MsgFreeze(msg);
  if (! LC_MsgAlloc(msg))
  {
    sprintf(cBuffer, STR_NOMEM " in LC_MsgPrepareSend (size=%ld)",
            (unsigned long)size);
    DDD_PrintError('E', 6600, cBuffer);
    HARD_EXIT;
  }

  return(size);
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

  assert(md->msgState==MSTATE_ALLOCATED);

  /* initiate asynchronous send */
  md->msgId = SendASync(VCHAN_TO(md->proc),
                        md->buffer, md->bufferSize, &error);

  md->msgState=MSTATE_COMM;
}


size_t LC_GetBufferSize (LC_MSGHANDLE msg)
{
  MSG_DESC *md = (MSG_DESC *) msg;

  return(md->bufferSize);
}


/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* Function:  LC_Connect                                                    */
/*                                                                          */
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
  if (nRecvs<0)
  {
    /* some processor raised an exception */
    sprintf(cBuffer,
            "Notify() raised exception #%d in LC_Connect()",
            -nRecvs);
    DDD_PrintError('E', 6624, cBuffer);

    /* automatically call LC_Cleanup() */
    DDD_NotifyEnd();
    LC_Cleanup();

    return(nRecvs);
  }


  if (nRecvs>procs-1)
  {
    sprintf(cBuffer, "cannot receive %d messages (must be less than %d)",
            nRecvs, procs-1);
    DDD_PrintError('E', 6622, cBuffer);
    DDD_NotifyEnd();
    return(EXCEPTION_LOWCOMM_CONNECT);
  }



#       if DebugLowComm<=7
  sprintf(cBuffer, "%4d: LC_Connect() nSends=%d nRecvs=%d\n",
          me, nSends, nRecvs);
  DDD_PrintDebug(cBuffer);
#       endif


  /* create array of receive message handles */
  if (nRecvs>0)
  {
    theRecvArray = (LC_MSGHANDLE *)AllocTmpReq(sizeof(LC_MSGHANDLE)*nRecvs, TMEM_ANY);
    if (theRecvArray==NULL)
    {
      DDD_PrintError('E', 6623, "out of memory in LC_Connect()");
      DDD_NotifyEnd();
      return(EXCEPTION_LOWCOMM_CONNECT);
    }
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
  {
    if (! IS_OK(DDD_GetChannels(nRecvs+nSends)))
    {
      DDD_PrintError('E', 6620, "couldn't get channels in LC_Connect()");
      return(EXCEPTION_LOWCOMM_CONNECT);
    }
  }


#       if DebugLowComm<=5
  DDD_DisplayTopo();
#       endif


  if (nRecvs>0)
  {
    if (! IS_OK(LC_PrepareRecv()))
      return(EXCEPTION_LOWCOMM_CONNECT);
  }


#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Connect() ready\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  return(nRecvs);
}



/****************************************************************************/
/*                                                                          */
/* Function:  LC_Abort                                                      */
/*                                                                          */
/****************************************************************************/

int LC_Abort (int exception)
{
  int retException;

  if (exception>EXCEPTION_LOWCOMM_USER)
  {
    DDD_PrintError('E', 6626,
                   "exception must be <=EXCEPTION_LOWCOMM_USER in LC_Abort()");
    HARD_EXIT;
  }


  DDD_NotifyBegin(exception);

#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Abort() exception=%d ...\n",
          me, exception);
  DDD_PrintDebug(cBuffer);
#       endif


  /* inform message receivers */
  retException = DDD_Notify();

  DDD_NotifyEnd();


#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Abort() ready, exception=%d.\n",
          me, retException);
  DDD_PrintDebug(cBuffer);
#       endif


  /* automatically call LC_Cleanup() */
  LC_Cleanup();

  return(retException);
}



/****************************************************************************/
/*                                                                          */
/* Function:  LC_Communicate                                                */
/*                                                                          */
/****************************************************************************/

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
    if (leftRecv>0) leftRecv = LC_PollRecv();
    if (leftSend>0) leftSend = LC_PollSend();
  } while (leftRecv>0 || leftSend>0);


#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Communicate() ready\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  return(theRecvArray);
}


/****************************************************************************/
/*                                                                          */
/* Function:  LC_Cleanup                                                    */
/*                                                                          */
/****************************************************************************/

void LC_Cleanup (void)
{
#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Cleanup() ...\n", me);
  DDD_PrintDebug(cBuffer);
#       endif

  if (nRecvs>0)
  {
    if (_RecvFree!=NULL)
      (*_RecvFree)(theRecvBuffer);

    theRecvBuffer=NULL;
  }

  if (theRecvArray!=NULL)
  {
    size_t size = sizeof(LC_MSGHANDLE)*nRecvs;
    FreeTmpReq(theRecvArray,size,TMEM_ANY);
    theRecvArray=NULL;
  }

  /* free recv queue */
  LC_FreeRecvQueue();

  /* free send queue */
  LC_FreeSendQueue();


#       if DebugLowComm<=9
  sprintf(cBuffer, "%4d: LC_Cleanup() ready\n", me);
  DDD_PrintDebug(cBuffer);
#       endif
}


/****************************************************************************/

#define LC_COLWIDTH   10
#define LC_DFLTNAME   "<?>"

/* construct name or default name */
#define LC_NAME(n)   (((n)!=NULL) ? (n) : LC_DFLTNAME)


static void LC_PrintMsgList (MSG_DESC *list)
{
  char buf[LC_COLWIDTH*2];
  MSG_DESC *md;
  MSG_TYPE *last_mt=NULL;
  size_t sum, comp_size[MAX_COMPONENTS];
  int i;

  for(md=list; md!=NULL; md=md->next)
  {
    MSG_TYPE *mt = md->msgType;

    if (mt!=last_mt)
    {
      /* msg-type changes, print new header */

      /* first, close part of msg-list with summary */
      if (last_mt!=NULL)
      {
        sprintf(cBuffer, "%4d:        = |", me);
        sum = 0;
        for(i=0; i<last_mt->nComps; i++)
        {
          sprintf(buf, "%9ld", comp_size[i]);
          strcat(cBuffer, buf);

          sum += comp_size[i];                                 /* horizontal sum */
        }
        sprintf(buf, "%9ld\n", sum); strcat(cBuffer, buf);
        DDD_PrintLine(cBuffer);
      }

      /* then, construct header */
      sprintf(cBuffer, "%4d:%9.9s |", me, LC_NAME(mt->name));
      for(i=0; i<mt->nComps; i++)
      {
        if (mt->comp[i].name!=NULL)
          sprintf(buf, "%9.9s", LC_NAME(mt->comp[i].name));
        else
          sprintf(buf, "%9d", i);
        strcat(cBuffer, buf);

        comp_size[i] = 0;
      }
      strcat(cBuffer, "        =\n"); DDD_PrintLine(cBuffer);
      last_mt = mt;
    }

    /* construct info about message components */
    sprintf(cBuffer, "%4d:%9d |", me, md->proc);
    sum = 0;
    for(i=0; i<mt->nComps; i++)
    {
      size_t s = md->chunks[i].size;

      sprintf(buf, "%9ld", s);
      strcat(cBuffer, buf);

      sum          += s;                     /* horizontal sum */
      comp_size[i] += s;                     /* vertical sum */
    }
    sprintf(buf, "%9ld\n", sum); strcat(cBuffer, buf);
    DDD_PrintLine(cBuffer);
  }

  /* close last part of msg-list with summary */
  if (last_mt!=NULL)
  {
    sprintf(cBuffer, "%4d:        = |", me);
    sum = 0;
    for(i=0; i<last_mt->nComps; i++)
    {
      sprintf(buf, "%9ld", comp_size[i]);
      strcat(cBuffer, buf);

      sum += comp_size[i];                     /* horizontal sum */
    }
    sprintf(buf, "%9ld\n", sum); strcat(cBuffer, buf);
    DDD_PrintLine(cBuffer);
  }

}


void LC_PrintSendMsgs (void)
{
  int p;

  for(p=0; p<procs; p++)
  {
    DDD_SyncAll();
    if (p==me)
      LC_PrintMsgList(LC_SendQueue);
  }
  DDD_SyncAll();
}

void LC_PrintRecvMsgs (void)
{
  int p;

  for(p=0; p<procs; p++)
  {
    DDD_SyncAll();
    if (p==me)
      LC_PrintMsgList(LC_RecvQueue);
  }
  DDD_SyncAll();
}


/****************************************************************************/
