// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ifcmds.c                                                      */
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
/*            96/01/24 kb  added use of object shortcut tables              */
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
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


#ifdef C_FRONTEND
void DDD_IFExecLocal (DDD_IF ifId, ExecProcPtr ExecProc)
{
#else
void DDD_IFExecLocal (DDD_IF *_ifId, ExecProcPtr ExecProc)
{
  DDD_IF ifId     = *_ifId;
#endif
  IF_PROC         *ifHead;

  /* if shortcut tables are invalid -> recompute */
  IFCheckShortcuts(ifId);

  ForIF(ifId,ifHead)
  {
    IFExecLoopObj(ExecProc, ifHead->objBA,  ifHead->nBA);
    IFExecLoopObj(ExecProc, ifHead->objAB,  ifHead->nAB);
    IFExecLoopObj(ExecProc, ifHead->objABA, ifHead->nABA);
  }
}




#ifdef C_FRONTEND
void DDD_IFAExecLocal (DDD_IF ifId, DDD_ATTR attr, ExecProcPtr ExecProc)
{
#else
void DDD_IFAExecLocal (DDD_IF *_ifId, DDD_ATTR *_attr, ExecProcPtr ExecProc)
{
  DDD_IF ifId   = *_ifId;
  DDD_ATTR attr   = *_attr;
#endif
  IF_PROC         *ifHead;

  /* if shortcut tables are invalid -> recompute */
  IFCheckShortcuts(ifId);

  ForIF(ifId,ifHead)
  {
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;

    IFExecLoopObj(ExecProc, ifAttr->objBA,  ifAttr->nBA);
    IFExecLoopObj(ExecProc, ifAttr->objAB,  ifAttr->nAB);
    IFExecLoopObj(ExecProc, ifAttr->objABA, ifAttr->nABA);
  }
}


#ifdef C_FRONTEND
void DDD_IFAExecLocalX (DDD_IF ifId, DDD_ATTR attr, ExecProcXPtr ExecProc)
{
  IF_PROC         *ifHead;

  /* this routine doesn't use shortcuts. */


  ForIF(ifId,ifHead)
  {
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;

    IFExecLoopCplX(ExecProc, ifAttr->cplBA,  ifAttr->nBA);
    IFExecLoopCplX(ExecProc, ifAttr->cplAB,  ifAttr->nAB);
    IFExecLoopCplX(ExecProc, ifAttr->cplABA, ifAttr->nABA);
  }
}

#endif


/****************************************************************************/


#ifdef C_FRONTEND
void DDD_IFExchange (DDD_IF ifId,
                     size_t itemSize,
                     ComProcPtr Gather,
                     ComProcPtr Scatter)
{
#else
void DDD_IFExchange (DDD_IF *_ifId,
                     size_t *_itemSize,
                     ComProcPtr Gather,
                     ComProcPtr Scatter)
{
  size_t itemSize = *_itemSize;
  DDD_IF ifId     = *_ifId;
#endif

  IF_PROC         *ifHead;
  int try;
  int recv_mesgs;

  /* if shortcut tables are invalid -> recompute */
  IFCheckShortcuts(ifId);


  STAT_ZEROTIMER;
  STAT_RESET1;
  STAT_SETCOUNT(20,0);
  STAT_SETCOUNT(21,0);
  STAT_SETCOUNT(22,0);

  /* allocate storage for in and out buffers */
  ForIF(ifId,ifHead)
  {
    IFGetMem(ifHead, itemSize, ifHead->nItems, ifHead->nItems);
    STAT_SETCOUNT(20,STAT_GETCOUNT(20)+ifHead->nItems);
    STAT_SETCOUNT(21,STAT_GETCOUNT(21)+2);
    STAT_SETCOUNT(22,STAT_GETCOUNT(22)+ifHead->lenBufIn+ifHead->lenBufOut);
  }
  STAT_TIMER1(63);

  /* init communication, initiate receives */
  recv_mesgs = IFInitComm(ifId, itemSize);

  /* build messages using gather-handler and send them away */
  ForIF(ifId,ifHead)
  {
    char *buffer;

    /* exchange BA and AB during send */
    buffer = ifHead->msgBufOut;
    buffer = IFCommLoopObj(Gather, ifHead->objBA,  buffer, itemSize, ifHead->nBA);
    buffer = IFCommLoopObj(Gather, ifHead->objAB,  buffer, itemSize, ifHead->nAB);
    buffer = IFCommLoopObj(Gather, ifHead->objABA, buffer, itemSize, ifHead->nABA);

    IFInitSend(ifHead);
  }

  /* poll receives and scatter data */
  for(try=0; try<MAX_TRIES && recv_mesgs>0; try++)
  {
    /* poll receive calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgIn!=-1)
      {
        int error = InfoARecv(ifHead->vc, ifHead->msgIn);
        /* TODO complete error handling */
        if (error==1)
        {
          char *buffer;

          recv_mesgs--;
          ifHead->msgIn=-1;

          /* get data using scatter-handler */
          buffer = ifHead->msgBufIn;
          buffer = IFCommLoopObj(Scatter,
                                 ifHead->objAB,  buffer, itemSize, ifHead->nAB);
          buffer = IFCommLoopObj(Scatter,
                                 ifHead->objBA,  buffer, itemSize, ifHead->nBA);
          buffer = IFCommLoopObj(Scatter,
                                 ifHead->objABA, buffer, itemSize, ifHead->nABA);
        }
      }
    }
  }

  /* finally poll send completion */
  if ((recv_mesgs>0) || (! IFPollSend(ifId)))
  {
    DDD_PrintError('E', 4200, "timeout in DDD_IFExchange");
    /* TODO complete error handling */
  }

  /* free memory */
  IFExitComm(ifId);

  STAT_TIMER1(60);
}



#ifdef C_FRONTEND
void DDD_IFOneway (DDD_IF ifId,
                   int dir,
                   size_t itemSize,
                   ComProcPtr Gather,
                   ComProcPtr Scatter)
{
#else
void DDD_IFOneway (DDD_IF *_ifId,
                   int *_dir,
                   size_t *_itemSize,
                   ComProcPtr Gather,
                   ComProcPtr Scatter)
{
  int dir      = *_dir;
  size_t itemSize = *_itemSize;
  DDD_IF ifId     = *_ifId;
#endif

  IF_PROC         *ifHead;
  int try;
  int recv_mesgs;

  /* if shortcut tables are invalid -> recompute */
  IFCheckShortcuts(ifId);


  STAT_ZEROTIMER;
  STAT_RESET1;
  STAT_SETCOUNT(20,0);
  STAT_SETCOUNT(21,0);
  STAT_SETCOUNT(22,0);

  /* allocate storage for in and out buffers */
  ForIF(ifId,ifHead)
  {
    int nOut, nIn;

    if (dir==IF_FORWARD)
    {
      nOut = ifHead->nAB; nIn = ifHead->nBA;
    }
    else
    {
      nOut = ifHead->nBA; nIn = ifHead->nAB;
    }

    IFGetMem(ifHead, itemSize, nIn+ifHead->nABA, nOut+ifHead->nABA);
    STAT_SETCOUNT(20,STAT_GETCOUNT(20)+ifHead->nItems);
    STAT_SETCOUNT(21,STAT_GETCOUNT(21)+2);
    STAT_SETCOUNT(22,STAT_GETCOUNT(22)+ifHead->lenBufIn+ifHead->lenBufOut);
  }
  STAT_TIMER1(63);

  /* init communication, initiate receives */
  recv_mesgs = IFInitComm(ifId, itemSize);

  /* build messages using gather-handler and send them away */
  ForIF(ifId,ifHead)
  {
    char     *buffer;
    int nOut;
    DDD_OBJ  *objOut;

    if (dir==IF_FORWARD)
    {
      nOut = ifHead->nAB;  objOut = ifHead->objAB;
    }
    else
    {
      nOut = ifHead->nBA;  objOut = ifHead->objBA;
    }

    buffer = IFCommLoopObj(Gather, objOut, ifHead->msgBufOut, itemSize, nOut);
    buffer = IFCommLoopObj(Gather,ifHead->objABA,buffer,itemSize,ifHead->nABA);

    IFInitSend(ifHead);
  }

  /* poll receives and scatter data */
  for(try=0; try<MAX_TRIES && recv_mesgs>0; try++)
  {
    /* poll receive calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgIn!=-1)
      {
        int error = InfoARecv(ifHead->vc, ifHead->msgIn);
        /* TODO complete error handling */
        if (error==1)
        {
          char     *buffer;
          int nIn;
          DDD_OBJ  *objIn;

          recv_mesgs--;
          ifHead->msgIn=-1;

          /* get data using scatter-handler */
          if (dir==IF_FORWARD)
          {
            nIn  = ifHead->nBA;  objIn = ifHead->objBA;
          }
          else
          {
            nIn  = ifHead->nAB;  objIn = ifHead->objAB;
          }

          buffer = IFCommLoopObj(Scatter,
                                 objIn, ifHead->msgBufIn, itemSize, nIn);
          buffer = IFCommLoopObj(Scatter,
                                 ifHead->objABA,buffer,itemSize,ifHead->nABA);
        }
      }
    }
  }

  /* finally poll send completion */
  if ((recv_mesgs>0) || (! IFPollSend(ifId)))
  {
    DDD_PrintError('E', 4201, "timeout in DDD_IFOneway");
    /* TODO complete error handling */
  }

  /* free memory */
  IFExitComm(ifId);

  STAT_TIMER1(60);
}



/****************************************************************************/


void DDD_IFAExchange (DDD_IF ifId,
                      DDD_ATTR attr,
                      size_t itemSize,
                      ComProcPtr Gather,
                      ComProcPtr Scatter)
{
  IF_PROC         *ifHead;
  int try;
  int recv_mesgs;

  /* if shortcut tables are invalid -> recompute */
  IFCheckShortcuts(ifId);


  /* allocate storage for in and out buffers */
  ForIF(ifId,ifHead)
  {
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;

    IFGetMem(ifHead, itemSize, ifAttr->nItems, ifAttr->nItems);
  }

  /* init communication, initiate receives */
  recv_mesgs = IFInitComm(ifId, itemSize);

  /* build messages using gather-handler and send them away */
  ForIF(ifId,ifHead)
  {
    char            *buffer;
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;

    /* exchange BA and AB during send */
    buffer = ifHead->msgBufOut;
    buffer = IFCommLoopObj(Gather,
                           ifAttr->objBA,  buffer, itemSize, ifAttr->nBA);
    buffer = IFCommLoopObj(Gather,
                           ifAttr->objAB,  buffer, itemSize, ifAttr->nAB);
    buffer = IFCommLoopObj(Gather,
                           ifAttr->objABA, buffer, itemSize, ifAttr->nABA);

    IFInitSend(ifHead);
  }

  /* poll receives and scatter data */
  for(try=0; try<MAX_TRIES && recv_mesgs>0; try++)
  {
    /* poll receive calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgIn!=-1)
      {
        int error = InfoARecv(ifHead->vc, ifHead->msgIn);
        /* TODO complete error handling */
        if (error==1)
        {
          char            *buffer;
          IF_ATTR     *ifAttr;

          recv_mesgs--;
          ifHead->msgIn=-1;

          /* get data using scatter-handler */
          /* find ifAttr */
          ifAttr = ifHead->ifAttr;
          while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
          {
            ifAttr = ifAttr->next;
          }
          if (ifAttr==NULL) continue;

          buffer = ifHead->msgBufIn;
          buffer = IFCommLoopObj(Scatter,
                                 ifAttr->objAB,  buffer, itemSize, ifAttr->nAB);
          buffer = IFCommLoopObj(Scatter,
                                 ifAttr->objBA,  buffer, itemSize, ifAttr->nBA);
          buffer = IFCommLoopObj(Scatter,
                                 ifAttr->objABA, buffer, itemSize, ifAttr->nABA);
        }
      }
    }
  }

  /* finally poll send completion */
  if ((recv_mesgs>0) || (! IFPollSend(ifId)))
  {
    DDD_PrintError('E', 4200, "timeout in DDD_IFRExchange");
    /* TODO complete error handling */
  }

  /* free memory */
  IFExitComm(ifId);
}




void DDD_IFAOneway (DDD_IF ifId,
                    DDD_ATTR attr,
                    int dir,
                    size_t itemSize,
                    ComProcPtr Gather,
                    ComProcPtr Scatter)
{
  IF_PROC         *ifHead;
  int try;
  int recv_mesgs;

  /* if shortcut tables are invalid -> recompute */
  IFCheckShortcuts(ifId);


  /* allocate storage for in and out buffers */
  ForIF(ifId,ifHead)
  {
    int nOut, nIn;
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;


    if (dir==IF_FORWARD)
    {
      nOut = ifAttr->nAB; nIn  = ifAttr->nBA;
    }
    else
    {
      nOut = ifAttr->nBA; nIn  = ifAttr->nAB;
    }

    IFGetMem(ifHead, itemSize, nIn+ifAttr->nABA, nOut+ifAttr->nABA);
  }

  /* init communication, initiate receives */
  recv_mesgs = IFInitComm(ifId, itemSize);

  /* build messages using gather-handler and send them away */
  ForIF(ifId,ifHead)
  {
    char            *buffer;
    int nOut;
    DDD_OBJ     *objOut;
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;

    if (dir==IF_FORWARD)
    {
      nOut = ifAttr->nAB;  objOut = ifAttr->objAB;
    }
    else
    {
      nOut = ifAttr->nBA;  objOut = ifAttr->objBA;
    }

    buffer = IFCommLoopObj(Gather, objOut, ifHead->msgBufOut, itemSize, nOut);
    buffer = IFCommLoopObj(Gather,ifAttr->objABA,buffer,itemSize,ifAttr->nABA);

    IFInitSend(ifHead);
  }

  /* poll receives and scatter data */
  for(try=0; try<MAX_TRIES && recv_mesgs>0; try++)
  {
    /* poll receive calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgIn!=-1)
      {
        int error = InfoARecv(ifHead->vc, ifHead->msgIn);
        /* TODO complete error handling */
        if (error==1)
        {
          char            *buffer;
          int nIn;
          DDD_OBJ     *objIn;
          IF_ATTR     *ifAttr;

          recv_mesgs--;
          ifHead->msgIn=-1;

          /* get data using scatter-handler */
          /* find ifAttr */
          ifAttr = ifHead->ifAttr;
          while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
          {
            ifAttr = ifAttr->next;
          }

          if (ifAttr==NULL) continue;

          if (dir==IF_FORWARD)
          {
            nIn  = ifAttr->nBA;  objIn = ifAttr->objBA;
          }
          else
          {
            nIn  = ifAttr->nAB;  objIn = ifAttr->objAB;
          }

          buffer = IFCommLoopObj(Scatter,
                                 objIn, ifHead->msgBufIn, itemSize, nIn);
          buffer = IFCommLoopObj(Scatter,
                                 ifAttr->objABA,buffer,itemSize,ifAttr->nABA);
        }
      }
    }
  }

  /* finally poll send completion */
  if ((recv_mesgs>0) || (! IFPollSend(ifId)))
  {
    DDD_PrintError('E', 4201, "timeout in DDD_IFROneway");
    /* TODO complete error handling */
  }

  /* free memory */
  IFExitComm(ifId);
}



/****************************************************************************/


void DDD_IFExchangeX (DDD_IF ifId,
                      size_t itemSize,
                      ComProcXPtr Gather,
                      ComProcXPtr Scatter)
{
  IF_PROC         *ifHead;
  int try;
  int recv_mesgs;

  /* this routine doesn't use shortcuts. */


  /* allocate storage for in and out buffers */
  ForIF(ifId,ifHead)
  {
    IFGetMem(ifHead, itemSize, ifHead->nItems, ifHead->nItems);
  }

  /* init communication, initiate receives */
  recv_mesgs = IFInitComm(ifId, itemSize);

  /* build messages using gather-handler and send them away */
  ForIF(ifId,ifHead)
  {
    char            *buffer;

    /* exchange BA and AB during send */
    buffer = ifHead->msgBufOut;
    buffer = IFCommLoopCplX(Gather,ifHead->cplBA, buffer,itemSize,ifHead->nBA);
    buffer = IFCommLoopCplX(Gather,ifHead->cplAB, buffer,itemSize,ifHead->nAB);
    buffer = IFCommLoopCplX(Gather,ifHead->cplABA,buffer,itemSize,ifHead->nABA);

    IFInitSend(ifHead);
  }

  /* poll receives and scatter data */
  for(try=0; try<MAX_TRIES && recv_mesgs>0; try++)
  {
    /* poll receive calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgIn!=-1)
      {
        int error = InfoARecv(ifHead->vc, ifHead->msgIn);
        /* TODO complete error handling */
        if (error==1)
        {
          char            *buffer;

          recv_mesgs--;
          ifHead->msgIn=-1;

          /* get data using scatter-handler */
          buffer = ifHead->msgBufIn;
          buffer = IFCommLoopCplX(Scatter,
                                  ifHead->cplAB,  buffer, itemSize, ifHead->nAB);
          buffer = IFCommLoopCplX(Scatter,
                                  ifHead->cplBA,  buffer, itemSize, ifHead->nBA);
          buffer = IFCommLoopCplX(Scatter,
                                  ifHead->cplABA, buffer, itemSize, ifHead->nABA);
        }
      }
    }
  }

  /* finally poll send completion */
  if ((recv_mesgs>0) || (! IFPollSend(ifId)))
  {
    DDD_PrintError('E', 4200, "timeout in DDD_IFExchangeX");
    /* TODO complete error handling */
  }

  /* free memory */
  IFExitComm(ifId);
}



void DDD_IFOnewayX (DDD_IF ifId,
                    int dir,
                    size_t itemSize,
                    ComProcXPtr Gather,
                    ComProcXPtr Scatter)
{
  IF_PROC         *ifHead;
  int recv_mesgs;
  int try;

  /* this routine doesn't use shortcuts. */


  /* allocate storage for in and out buffers */
  ForIF(ifId,ifHead)
  {
    int nOut, nIn;

    if (dir==IF_FORWARD)
    {
      nOut = ifHead->nAB; nIn  = ifHead->nBA;
    }
    else
    {
      nOut = ifHead->nBA; nIn  = ifHead->nAB;
    }

    IFGetMem(ifHead, itemSize, nIn+ifHead->nABA, nOut+ifHead->nABA);
  }

  /* init communication, initiate receives */
  recv_mesgs = IFInitComm(ifId, itemSize);

  /* build messages using gather-handler and send them away */
  ForIF(ifId,ifHead)
  {
    char            *buffer;
    int nOut;
    COUPLING    **cplOut;

    if (dir==IF_FORWARD)
    {
      nOut = ifHead->nAB;  cplOut = ifHead->cplAB;
    }
    else
    {
      nOut = ifHead->nBA;  cplOut = ifHead->cplBA;
    }

    buffer = IFCommLoopCplX(Gather, cplOut,ifHead->msgBufOut,itemSize,nOut);
    buffer = IFCommLoopCplX(Gather,ifHead->cplABA,buffer,itemSize,ifHead->nABA);

    IFInitSend(ifHead);
  }

  /* poll receives and scatter data */
  for(try=0; try<MAX_TRIES && recv_mesgs>0; try++)
  {
    /* poll receive calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgIn!=-1)
      {
        int error = InfoARecv(ifHead->vc, ifHead->msgIn);
        /* TODO complete error handling */
        if (error==1)
        {
          char            *buffer;
          int nIn;
          COUPLING    **cplIn;

          recv_mesgs--;
          ifHead->msgIn=-1;

          /* get data using scatter-handler */
          if (dir==IF_FORWARD)
          {
            nIn  = ifHead->nBA;  cplIn = ifHead->cplBA;
          }
          else
          {
            nIn  = ifHead->nAB;  cplIn = ifHead->cplAB;
          }

          buffer = IFCommLoopCplX(Scatter,
                                  cplIn, ifHead->msgBufIn, itemSize, nIn);
          buffer = IFCommLoopCplX(Scatter,
                                  ifHead->cplABA,buffer,itemSize,ifHead->nABA);
        }
      }
    }
  }

  /* finally poll send completion */
  if ((recv_mesgs>0) || (! IFPollSend(ifId)))
  {
    DDD_PrintError('E', 4201, "timeout in DDD_IFOnewayX");
    /* TODO complete error handling */
  }

  /* free memory */
  IFExitComm(ifId);
}



void DDD_IFAExchangeX (DDD_IF ifId,
                       DDD_ATTR attr,
                       size_t itemSize,
                       ComProcXPtr Gather,
                       ComProcXPtr Scatter)
{
  IF_PROC         *ifHead;
  int try;
  int recv_mesgs;

  /* this routine doesn't use shortcuts. */


  /* allocate storage for in and out buffers */
  ForIF(ifId,ifHead)
  {
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;

    IFGetMem(ifHead, itemSize, ifAttr->nItems, ifAttr->nItems);
  }

  /* init communication, initiate receives */
  recv_mesgs = IFInitComm(ifId, itemSize);

  /* build messages using gather-handler and send them away */
  ForIF(ifId,ifHead)
  {
    char            *buffer;
    IF_ATTR     *ifAttr;

    /* find ifAttr */
    ifAttr = ifHead->ifAttr;
    while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
    {
      ifAttr = ifAttr->next;
    }
    if (ifAttr==NULL) continue;

    /* exchange BA and AB during send */
    buffer = ifHead->msgBufOut;
    buffer = IFCommLoopCplX(Gather,
                            ifAttr->cplBA,  buffer, itemSize, ifAttr->nBA);
    buffer = IFCommLoopCplX(Gather,
                            ifAttr->cplAB,  buffer, itemSize, ifAttr->nAB);
    buffer = IFCommLoopCplX(Gather,
                            ifAttr->cplABA, buffer, itemSize, ifAttr->nABA);

    IFInitSend(ifHead);
  }

  /* poll receives and scatter data */
  for(try=0; try<MAX_TRIES && recv_mesgs>0; try++)
  {
    /* poll receive calls */
    ForIF(ifId,ifHead)
    {
      if (ifHead->msgIn!=-1)
      {
        int error = InfoARecv(ifHead->vc, ifHead->msgIn);
        /* TODO complete error handling */
        if (error==1)
        {
          char            *buffer;
          IF_ATTR     *ifAttr;

          recv_mesgs--;
          ifHead->msgIn=-1;

          /* get data using scatter-handler */
          /* find ifAttr */
          ifAttr = ifHead->ifAttr;
          while ((ifAttr!=NULL) && (ifAttr->attr!=attr))
          {
            ifAttr = ifAttr->next;
          }
          if (ifAttr==NULL) continue;

          buffer = ifHead->msgBufIn;
          buffer = IFCommLoopCplX(Scatter,
                                  ifAttr->cplAB,  buffer, itemSize, ifAttr->nAB);
          buffer = IFCommLoopCplX(Scatter,
                                  ifAttr->cplBA,  buffer, itemSize, ifAttr->nBA);
          buffer = IFCommLoopCplX(Scatter,
                                  ifAttr->cplABA, buffer, itemSize, ifAttr->nABA);
        }
      }
    }
  }

  /* finally poll send completion */
  if ((recv_mesgs>0) || (! IFPollSend(ifId)))
  {
    DDD_PrintError('E', 4200, "timeout in DDD_IFRExchange");
    /* TODO complete error handling */
  }

  /* free memory */
  IFExitComm(ifId);
}

/****************************************************************************/
