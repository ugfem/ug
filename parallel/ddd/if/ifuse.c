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
  size_t sizeIn  = itemSize * lenIn;
  size_t sizeOut = itemSize * lenOut;

  BufferCreate(ifHead->bufIn,  sizeIn);
  BufferCreate(ifHead->bufOut, sizeOut);
}



/*
        initiate asynchronous receive calls,
        return number of messages to be received
 */
int IFInitComm (DDD_IF ifId)
{
  IF_PROC   *ifHead;
  int error;
  int recv_mesgs;


  /* MarkHeap(); */

  recv_mesgs = 0;

  /* get memory and initiate receive calls */
  ForIF(ifId,ifHead)
  {
    if (! BufferIsEmpty(ifHead->bufIn))
    {
      ifHead->msgIn =
        RecvASync(ifHead->vc,
                  BufferMem(ifHead->bufIn), BufferLen(ifHead->bufIn),
                  &error);
      if (ifHead->msgIn==0)
      {
        DDD_PrintError('E', 4225, "PPIF's RecvASync() failed in IF-Comm");
        HARD_EXIT;
      }

      recv_mesgs++;
    }
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

  if (DDD_GetOption(OPT_IF_REUSE_BUFFERS) == OPT_OFF)
  {
    ForIF(ifId,ifHead)
    {
      BufferFree(ifHead->bufIn);
      BufferFree(ifHead->bufOut);
    }
  }

  /* ReleaseHeap(); */
}



/*
        initiate single asynchronous send call
 */
void IFInitSend (IF_PROC *ifHead)
{
  int error;

  if (! BufferIsEmpty(ifHead->bufOut))
  {
    ifHead->msgOut =
      SendASync(ifHead->vc,
                BufferMem(ifHead->bufOut), BufferLen(ifHead->bufOut),
                &error);
    if (ifHead->msgOut==0)
    {
      DDD_PrintError('E', 4226, "PPIF's SendASync() failed in IF-Comm");
      HARD_EXIT;
    }

    send_mesgs++;
  }
}



/*
        poll asynchronous send calls,
        return if ready
 */
int IFPollSend (DDD_IF ifId)
{
  unsigned long tries;

  for(tries=0; tries<MAX_TRIES && send_mesgs>0; tries++)
  {
    IF_PROC   *ifHead;

    /* poll send calls */
    ForIF(ifId,ifHead)
    {
      if ((! BufferIsEmpty(ifHead->bufOut)) && ifHead->msgOut!=-1)
      {
        int error = InfoASend(ifHead->vc, ifHead->msgOut);
        if (error==-1)
        {
          sprintf(cBuffer,
                  "PPIF's InfoASend() failed for send to proc=%d in IF-Comm",
                  ifHead->proc);
          DDD_PrintError('E', 4220, cBuffer);
          HARD_EXIT;
        }

        if (error==1)
        {
          send_mesgs--;
          ifHead->msgOut=-1;

                                        #ifdef CtrlTimeoutsDetailed
          printf("%4d: IFCTRL %02d send-completed    to "
                 "%4d after %10ld, size %ld\n",
                 me, ifId, ifHead->proc,
                 (unsigned long)tries,
                 (unsigned long)BufferLen(ifHead->bufOut));
                                        #endif
        }
      }
    }
  }

        #ifdef CtrlTimeouts
  if (send_mesgs==0)
  {
    printf("%4d: IFCTRL %02d send-completed    all after %10ld tries\n",
           me, ifId, (unsigned long)tries);
  }
        #endif

  return(send_mesgs==0);
}



/****************************************************************************/


/*
        do loop over single list of couplings,
        copy object data from/to message buffer

        fast version: uses object pointer shortcut
 */
#if defined(C_FRONTEND) || defined(F_FRONTEND)
char *IFCommLoopObj (ComProcPtr LoopProc,
                     IFObjPtr *obj,
                     char *buffer,
                     size_t itemSize,
                     int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
#if defined(C_FRONTEND)
    error = (*LoopProc)(obj[i], buffer);
#endif
#ifdef F_FRONTEND
    error = (*LoopProc)(obj+i, buffer);
#endif
    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }

  return(buffer);
}
#endif

#ifdef CPP_FRONTEND
char *IFCommLoopObjGather (DDD_GatherScatter& gs,
                           IFObjPtr *obj,
                           char *buffer,
                           size_t itemSize,
                           int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
    error = gs.Gather(obj[i], buffer);
    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }

  return(buffer);
}

char *IFCommLoopObjScatter (DDD_GatherScatter& gs,
                            IFObjPtr *obj,
                            char *buffer,
                            size_t itemSize,
                            int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
    error = gs.Scatter(obj[i], buffer);
    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }

  return(buffer);
}
#endif


#ifndef CPP_FRONTEND  /* for debugging */

/*
        simple variant of above routine. dont communicate,
        but call an application's routine.
 */
void IFExecLoopObj (ExecProcPtr LoopProc, IFObjPtr *obj, int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++)
  {
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    error = (*LoopProc)(obj[i]);
#endif
#ifdef F_FRONTEND
    error = (*LoopProc)(obj+i);
#endif
    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
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
#if defined(C_FRONTEND)
    error = (*LoopProc)(OBJ_OBJ(cpl[i]->obj), buffer);
#endif
#if defined(CPP_FRONTEND)
    // TODO: dirty cast in first argument!
    error = (*LoopProc)((DDD_Object*)(cpl[i]->obj), buffer);
#endif
#ifdef F_FRONTEND
    error = (*LoopProc)((IFObjPtr *) &(OBJ_INDEX(cpl[i]->obj)), buffer);
#endif
    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
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
#if defined(C_FRONTEND)
    error = (*LoopProc)(OBJ_OBJ(cpl[i]->obj),
                        buffer, CPL_PROC(cpl[i]), cpl[i]->prio);
#endif
#if defined(CPP_FRONTEND)
    // TODO: dirty cast in first argument!
    error = (*LoopProc)((DDD_Object*)(cpl[i]->obj),
                        buffer, CPL_PROC(cpl[i]), cpl[i]->prio);
#endif
#ifdef F_FRONTEND
    error = (*LoopProc)((IFObjPtr *) &(OBJ_INDEX(cpl[i]->obj)),
                        buffer, (DDD_PROC *) &(CPL_PROC(cpl[i])),
                        (DDD_PRIO *) &(cpl[i]->prio));
#endif

    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }

  return(buffer);
}


#if defined(C_FRONTEND) || defined(CPP_FRONTEND)

/*
        simple variant of above routine. dont communicate,
        but call an application's routine.
 */
void IFExecLoopCplX (ExecProcXPtr LoopProc, COUPLING **cpl, int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++)
  {
#if defined(C_FRONTEND)
    error = (*LoopProc)(OBJ_OBJ(cpl[i]->obj), CPL_PROC(cpl[i]), cpl[i]->prio);
#endif
#if defined(CPP_FRONTEND)
    // TODO: dirty cast in first argument!
    error = (*LoopProc)((DDD_Object*)(cpl[i]->obj), CPL_PROC(cpl[i]), cpl[i]->prio);
#endif
#ifdef F_FRONTEND
    error = (*LoopProc)((IFObjPtr *) &(OBJ_INDEX(cpl[i]->obj)),
                        (DDD_PROC *) &(CPL_PROC(cpl[i])),
                        (DDD_PRIO *) &(cpl[i]->prio));
#endif
    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }
}

#endif


#endif  /* for debugging */


/****************************************************************************/


/***
        interface loop functions for STD_INTERFACE communication.
        (not DDD_OBJ will be passed as a parameter, but DDD_HDR
        instead).
 ***/


/*
        do loop over single list of couplings,
        copy object data from/to message buffer
 */
char *IFCommHdrLoopCpl (ComProcHdrPtr LoopProc,
                        COUPLING **cpl,
                        char *buffer,
                        size_t itemSize,
                        int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
    error = (*LoopProc)(cpl[i]->obj, buffer);

    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }

  return(buffer);
}


/*
        simple variant of above routine. dont communicate,
        but call an application's routine.
 */
void IFExecHdrLoopCpl (ExecProcHdrPtr LoopProc, COUPLING **cpl, int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++)
  {
    error = (*LoopProc)(cpl[i]->obj);

    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }
}


/*
        do loop over single list of couplings,
        copy object data from/to message buffer

        extended version: call ComProc with extended parameters
 */
char *IFCommHdrLoopCplX (ComProcHdrXPtr LoopProc,
                         COUPLING **cpl,
                         char *buffer,
                         size_t itemSize,
                         int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++, buffer+=itemSize)
  {
    error = (*LoopProc)(cpl[i]->obj,
                        buffer, CPL_PROC(cpl[i]), cpl[i]->prio);

    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }

  return(buffer);
}


/*
        simple variant of above routine. dont communicate,
        but call an application's routine.

        extended version: call ExecProc with extended parameters
 */
void IFExecHdrLoopCplX (ExecProcHdrXPtr LoopProc, COUPLING **cpl, int nItems)
{
  int i, error;

  for(i=0; i<nItems; i++)
  {
    error = (*LoopProc)(cpl[i]->obj, CPL_PROC(cpl[i]), cpl[i]->prio);

    /* TODO: check error-value from IF-LoopProc and issue warning or HARD_EXIT */
  }
}



/****************************************************************************/
