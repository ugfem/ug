// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ifuse.c                                                       */
/*                                                                          */
/* Purpose:   routines concerning interfaces between processors             */
/*            part 2: usage of DDD interfaces                               */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            94/03/03 kb  complete rewrite                                 */
/*            94/09/12 kb  IFExchange & IFOneway rewrite, two bugs fixed    */
/*            94/09/21 kb  created from if.c                                */
/*            95/01/13 kb  added range functionality                        */
/*            95/07/26 kb  overlapping of gather/scatter and communication  */
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
#include "if.h"


/****************************************************************************/
/*                                                                          */
/* variables global to this source file (static)                            */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


static int send_mesgs;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/*
        allocate memory for message buffers,
        one for send, one for receive
 */
void IFGetMem (IF_PROC *ifHead, size_t itemSize, int lenIn, int lenOut)
{
  ifHead->lenBufIn  = itemSize * lenIn;
  ifHead->lenBufOut = itemSize * lenOut;

  ifHead->msgBufIn  = (char *) AllocMsg(ifHead->lenBufIn);
  ifHead->msgBufOut = (char *) AllocMsg(ifHead->lenBufOut);
}



/*
        initiate asynchronous receive calls,
        return number of messages to be received
 */
int IFInitComm (DDD_IF ifId, size_t itemSize)
{
  IF_PROC   *ifHead;
  int error;
  int recv_mesgs;


  MarkHeap();

  recv_mesgs = 0;

  /* get memory and initiate receive calls */
  ForIF(ifId,ifHead)
  {
    ifHead->msgIn =
      RecvASync(ifHead->vc,
                ifHead->msgBufIn, ifHead->lenBufIn,
                &error);

    /* TODO error abfrage */
    recv_mesgs++;
  }

  send_mesgs = 0;

  return recv_mesgs;
}



/*
        cleanup memory
 */
void IFExitComm (DDD_IF ifId)
{
  IF_PROC   *ifHead;

  ForIF(ifId,ifHead)
  {
    FreeMsg(ifHead->msgBufIn);
    FreeMsg(ifHead->msgBufOut);
  }

  ReleaseHeap();
}



/*
        initiate single asynchronous send call
 */
void IFInitSend (IF_PROC *ifHead)
{
  int error;

  ifHead->msgOut =
    SendASync(ifHead->vc,
              ifHead->msgBufOut, ifHead->lenBufOut,
              &error);

  /* TODO error abfrage */
  send_mesgs++;
}



/*
        poll asynchronous send calls,
        return if ready
 */
int IFPollSend (DDD_IF ifId)
{
  int try;

  for(try=0; try<MAX_TRIES && send_mesgs>0; try++)
  {
    IF_PROC   *ifHead;

    /* poll send calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgOut!=-1)
      {
        int error = InfoASend(ifHead->vc, ifHead->msgOut);
        /* TODO complete error handling */
        if (error==1) {
          send_mesgs--;
          ifHead->msgOut=-1;
        }
      }
    }
  }

  return(send_mesgs==0);
}



/****************************************************************************/


/*
        do loop over single list of couplings,
        copy object data from/to message buffer

        fast version: uses object pointer shortcut
 */
char *IFCommLoopObj (ComProcPtr LoopProc,
                     DDD_OBJ *obj,
                     char *buffer,
                     size_t itemSize,
                     int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
#ifdef C_FRONTEND
    error = (*LoopProc)(obj[i], buffer);
#else
    error = (*LoopProc)(obj+i, buffer);
#endif
    /* TODO error abfrage */
  }

  return(buffer);
}


/*
        simple variant of above routine. dont communicate,
        but call an application's routine.
 */
void IFExecLoopObj (ExecProcPtr LoopProc, DDD_OBJ *obj, int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++)
  {
#ifdef C_FRONTEND
    error = (*LoopProc)(obj[i]);
#else
    error = (*LoopProc)(obj+i);
#endif
    /* TODO error abfrage */
  }
}




/*
        do loop over single list of couplings,
        copy object data from/to message buffer

        unnecessary indirect addressing!
        (-> CPL -> DDD_HDR.typ -> header offset -> object address)
 */
char *IFCommLoopCpl (ComProcPtr LoopProc,
                     COUPLING **cpl,
                     char *buffer,
                     size_t itemSize,
                     int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
#ifdef C_FRONTEND
    error = (*LoopProc)(OBJ_OBJ(cpl[i]->obj), buffer);
#else
    error = (*LoopProc)((DDD_OBJ *) &(OBJ_INDEX(cpl[i]->obj)), buffer);
#endif
    /* TODO error abfrage */
  }

  return(buffer);
}




/*
        do loop over single list of couplings,
        copy object data from/to message buffer

        extended version: call ComProc with extended parameters

        (necessary) indirect addressing!
        (-> CPL -> DDD_HDR.typ -> header offset -> object address)
 */
char *IFCommLoopCplX (ComProcXPtr LoopProc,
                      COUPLING **cpl,
                      char *buffer,
                      size_t itemSize,
                      int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
#ifdef C_FRONTEND
    error = (*LoopProc)(OBJ_OBJ(cpl[i]->obj),
                        buffer, cpl[i]->proc, cpl[i]->prio);
#else
    error = (*LoopProc)((DDD_OBJ *) &(OBJ_INDEX(cpl[i]->obj)),
                        buffer, (DDD_PROC *) &(cpl[i]->proc),
                        (DDD_PRIO *) &(cpl[i]->prio));
#endif

    /* TODO error abfrage */
  }

  return(buffer);
}


#ifdef C_FRONTEND

/*
        simple variant of above routine. dont communicate,
        but call an application's routine.
 */
void IFExecLoopCplX (ExecProcXPtr LoopProc, COUPLING **cpl, int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++)
  {
    error = (*LoopProc)(OBJ_OBJ(cpl[i]->obj), cpl[i]->proc, cpl[i]->prio);
    /* TODO error abfrage */
  }
}

#endif

/****************************************************************************/
