// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      notify.c                                                      */
/*                                                                          */
/* Purpose:   notifies destinations for communication with globally         */
/*            unknown topology                                              */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/01/17 kb  begin                                            */
/*            95/04/06 kb  added SpreadNotify                               */
/*            96/07/12 kb  united xxxNotify functions to one DDD_Notify()   */
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
#include <config.h>
#include <stdlib.h>
#include <stdio.h>

#include "dddi.h"
#include "basic/notify.h"

USING_UG_NAMESPACES

/* PPIF namespace: */
USING_PPIF_NAMESPACE

  START_UGDIM_NAMESPACE

#define DebugNotify   10  /* 0 is all, 10 is off */



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


#define MAX_INFOS    ((procs)*(MAX((1+procs),10)))



/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/


enum NotifyTypes {MYSELF,KNOWN,DUMMY,UNKNOWN};


typedef struct {
  short from, to;                       /* source and destination processor */
  unsigned short flag;                  /* one of NotifyTypes */
  size_t size;                          /* message size */
} NOTIFY_INFO;

#define PROC_INVALID_TEMP   -1


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



static NOTIFY_INFO *allInfoBuffer;
static NOTIFY_DESC *theDescs;
static int      *theRouting;
static int maxInfos, lastInfo, nSendDescs;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



void NotifyInit (void)
{
  /* allocate memory */
  theRouting = (int *) AllocFix(procs*sizeof(int));
  if (theRouting==NULL)
  {
    DDD_PrintError('E', 6301, STR_NOMEM " in NotifyInit");
    HARD_EXIT;
  }


  maxInfos = MAX_INFOS;                 /* TODO maximum value, just for testing */


  /* init local array for all Info records */
  allInfoBuffer = (NOTIFY_INFO *) AllocFix(maxInfos*sizeof(NOTIFY_INFO));
  if (allInfoBuffer==NULL)
  {
    DDD_PrintError('E', 6300, STR_NOMEM " in NotifyInit");
    HARD_EXIT;
  }


  /* allocate array of NOTIFY_DESCs */
  if (procs>1)
  {
    theDescs = (NOTIFY_DESC *) AllocTmp(sizeof(NOTIFY_DESC)*(procs-1));
  }
  else
  {
    theDescs = NULL;
  }
}


void NotifyExit (void)
{
  /* free memory */
  FreeFix(theRouting);
  FreeFix(allInfoBuffer);

  if (theDescs!=NULL)
  {
    FreeTmp(theDescs,sizeof(NOTIFY_DESC)*(procs-1));
  }
}


/****************************************************************************/



static int sort_XferInfos (const void *e1, const void *e2)
{
  NOTIFY_INFO *ci1, *ci2;

  ci1 = (NOTIFY_INFO *)e1;
  ci2 = (NOTIFY_INFO *)e2;

  if (ci1->to < ci2->to) return(-1);
  if (ci1->to > ci2->to) return(1);

  if (ci1->from < ci2->from) return(-1);
  if (ci1->from == ci2->from) return(0);
  return(1);
}


static int sort_XferFlags (const void *e1, const void *e2)
{
  NOTIFY_INFO *ci1, *ci2;

  ci1 = (NOTIFY_INFO *)e1;
  ci2 = (NOTIFY_INFO *)e2;

  if (ci1->flag < ci2->flag) return(-1);
  if (ci1->flag == ci2->flag) return(0);
  return(1);
}


static int sort_XferRouting (const void *e1, const void *e2)
{
  NOTIFY_INFO *ci1, *ci2;

  ci1 = (NOTIFY_INFO *)e1;
  ci2 = (NOTIFY_INFO *)e2;

  if (theRouting[ci1->to] < theRouting[ci2->to]) return(-1);
  if (theRouting[ci1->to] == theRouting[ci2->to]) return(0);
  return(1);
}



NOTIFY_INFO *NotifyPrepare (void)
{
  NOTIFY_INFO  *allInfos;

#if     DebugNotify<=4
  printf("%4d:    NotifyPrepare\n", me);
  fflush(stdout);
#endif

  /* init local array for all Info records */
  allInfos = allInfoBuffer;


  /* init local routing array */
  theRouting[me] = -1;


  /* dummy Info if there is no message to be send */
  allInfos[0].from = me;
  allInfos[0].to   = PROC_INVALID_TEMP;
  allInfos[0].size = 0;
  allInfos[0].flag = DUMMY;
  lastInfo = 1;

  return(allInfos);
}



/****************************************************************************/

/*
        If the parameter 'exception' is !=0, this processor invokes a
        global exception, which will cause all processors to abort this
        notify procedure and return the exception code with flipped sign.
        If more processors issue exception codes, the maximum will be
        communicated.
 */

int NotifyTwoWave (NOTIFY_INFO *allInfos, int lastInfo, int exception)
{
  NOTIFY_INFO  *newInfos;
  int l, i, j, n, unknownInfos, myInfos;
  int local_exception = exception;

#if     DebugNotify<=4
  printf("%4d:    NotifyTwoWave, lastInfo=%d\n", me, lastInfo);
  fflush(stdout);
#endif

  /* BOTTOM->TOP WAVE */
  /* get local Info lists from downtree */
  for(l=degree-1; l>=0; l--)
  {
    GetConcentrate(l, &n, sizeof(int));

    if (n<0)
    {
      /* exception from downtree, propagate */
      if (-n > local_exception)
        local_exception = -n;
    }

    if (lastInfo+n >= maxInfos) {
      DDD_PrintError('E', 6321, "msg-info array overflow in NotifyTwoWave");
      local_exception = EXCEPTION_NOTIFY;

      /* receive data, but put it onto dummy position */
      GetConcentrate(l, allInfos, n*sizeof(NOTIFY_INFO));
    }
    else
    {
      if (n>0)
        GetConcentrate(l, &(allInfos[lastInfo]), n*sizeof(NOTIFY_INFO));
    }

    /* construct routing table */
    for(i=0; i<n; i++)
      theRouting[allInfos[lastInfo+i].from] = l;

    if (n>0)
      lastInfo += n;
  }


  if (! local_exception)
  {
    /* determine target direction in tree */
    /* TODO: eventually extra solution for root node!
                     (it knows all flags are MYSELF or KNOWN!)  */
    qsort(allInfos, lastInfo, sizeof(NOTIFY_INFO), sort_XferInfos);
    i = j = 0;
    unknownInfos = lastInfo;
    myInfos = 0;
    while (i<lastInfo && allInfos[j].to==PROC_INVALID_TEMP)
    {
      if (allInfos[j].from==allInfos[i].to)
      {
        allInfos[i].flag = (allInfos[i].to==me) ? MYSELF : KNOWN;
        unknownInfos--;
        if (allInfos[i].to==me)
          myInfos++;
        i++;
      } else {
        if (allInfos[j].from<allInfos[i].to)
          j++;
        else
          i++;
      }
    }
    qsort(allInfos, lastInfo, sizeof(NOTIFY_INFO), sort_XferFlags);


    /* send local Info list uptree, but only unknown Infos */
    newInfos = &allInfos[lastInfo-unknownInfos];
    Concentrate(&unknownInfos, sizeof(int));
    Concentrate(newInfos, unknownInfos*sizeof(NOTIFY_INFO));
    lastInfo -= unknownInfos;

                #if     DebugNotify<=1
    for(i=0; i<unknownInfos; i++)
    {
      printf("%4d:    NotifyTwoWave, "
             "send uptree unknown %d/%d (%d|%d;%d)\n",
             me, i, unknownInfos,
             newInfos[i].to, newInfos[i].from, newInfos[i].size);
    }
                #endif

  }
  else
  {
    /* we have an exception somewhere in the processor tree */
    /* propagate it */
    int neg_exception = -local_exception;
    Concentrate(&neg_exception, sizeof(int));
    /* don't need to send data now */
  }

#if     DebugNotify<=3
  printf("%4d:    NotifyTwoWave, wave 1 ready\n", me);
  fflush(stdout);
#endif



  /* TOP->BOTTOM WAVE */

  /* get Infos local to my subtree from uptree */
  unknownInfos = 0;
  GetSpread(&unknownInfos, sizeof(int));
  if (unknownInfos<0)
  {
    /* exception from downtree, propagate */
    if (-unknownInfos > local_exception)
      local_exception = -unknownInfos;
  }

  if (unknownInfos>0)
  {
    GetSpread(newInfos, unknownInfos*sizeof(NOTIFY_INFO));
    lastInfo += unknownInfos;
  }

  if (! local_exception)
  {
    /* sort Infos according to routing */
    qsort(allInfos, lastInfo, sizeof(NOTIFY_INFO), sort_XferRouting);

                #if     DebugNotify<=1
    for(i=0; i<lastInfo; i++)
    {
      printf("%4d:    NotifyTwoWave, "
             "sorted for routing  %d/%d (%d|%d;%d)\n",
             me, i, lastInfo,
             allInfos[i].to, allInfos[i].from, allInfos[i].size);
    }
                #endif

    /* send relevant Infos downtree */
    i = 0;
    unknownInfos = lastInfo;
    while ((i<unknownInfos)&&(allInfos[i].to==me)) i++;
    lastInfo = i;
    for(l=0; l<degree; l++)
    {
      j = i;
      while ((i<unknownInfos)&&(theRouting[allInfos[i].to]==l)) i++;
      j = i-j;

      Spread(l, &j, sizeof(int));
      if (j>0)
        Spread(l, &allInfos[i-j], j*sizeof(NOTIFY_INFO));
    }


    /* reuse theDescs-array for registering messages to be received */
    for(i=0; i<lastInfo; i++)
    {
      theDescs[i].proc = allInfos[i].from;
      theDescs[i].size = allInfos[i].size;
    }

                #if     DebugNotify<=3
    printf("%4d:    NotifyTwoWave, "
           "wave 2 ready, nRecv=%d\n", me, lastInfo);
    fflush(stdout);
                #endif
  }
  else
  {
    /* we received an exception from uptree, propagate it */
    for(l=0; l<degree; l++)
    {
      int neg_exception = -local_exception;
      Spread(l, &neg_exception, sizeof(int));
      /* dont send any data */
    }

                #if     DebugNotify<=3
    printf("%4d:    NotifyTwoWave, "
           "wave 2 ready, Exception=%d\n", me, local_exception);
    fflush(stdout);
                #endif

    return(-local_exception);
  }

  return(lastInfo);
}



/****************************************************************************/


NOTIFY_DESC *DDD_NotifyBegin (int n)
{
  nSendDescs = n;

  /* allocation of theDescs is done in NotifyInit() */

  if (n>procs-1)
  {
    DDD_PrintError('E', 6340,
                   "more send-messages than other processors in DDD_NotifyBegin");
    return(NULL);
  }

  return(theDescs);
}


void DDD_NotifyEnd (void)
{
  /* free'ing of theDescs is done in NotifyExit() */
}


int DDD_Notify (void)
{
  NOTIFY_INFO  *allInfos;
  int i, nRecvMsgs;

  /* get storage for local info list */
  allInfos = NotifyPrepare();
  if (allInfos==NULL) return(ERROR);

  if (nSendDescs<0)
  {
    /* this processor is trying to send a global notification
       message. this is necessary for communicating fatal error
       conditions to all other processors. */

    sprintf(cBuffer, "proc %d is sending global exception #%d"
            " in DDD_Notify()", me, -nSendDescs);
    DDD_PrintError('W', 6312, cBuffer);

    /* notify partners */
    nRecvMsgs = NotifyTwoWave(allInfos, lastInfo, -nSendDescs);
  }
  else
  {
    /* convert message list to local Info list */
    for(i=0; i<nSendDescs; i++)
    {
                        #if     DebugNotify<=4
      printf("%4d:    Notify send msg #%02d to %3d size=%d\n", me,
             lastInfo, theDescs[i].proc, theDescs[i].size);
                        #endif

      if (theDescs[i].proc==me) {
        sprintf(cBuffer, "proc %d is trying to send message to itself"
                " in DDD_Notify()", me);
        DDD_PrintError('E', 6310, cBuffer);
        return(ERROR);
      }
      if (theDescs[i].proc>=procs) {
        sprintf(cBuffer, "proc %d is trying to send message to proc %d"
                " in DDD_Notify()", me, theDescs[i].proc);
        DDD_PrintError('E', 6311, cBuffer);
        return(ERROR);
      }

      allInfos[lastInfo].from = me;
      allInfos[lastInfo].to   = theDescs[i].proc;
      allInfos[lastInfo].size = theDescs[i].size;
      allInfos[lastInfo].flag = UNKNOWN;
      lastInfo++;
    }

    /* notify partners */
    nRecvMsgs = NotifyTwoWave(allInfos, lastInfo, 0);
  }


#       if      DebugNotify<=4
  for(i=0; i<nRecvMsgs; i++)
  {
    printf("%4d:    Notify recv msg #%02d from %3d size=%d\n", me,
           lastInfo, theDescs[i].proc, theDescs[i].size);
  }
#       endif


  return(nRecvMsgs);
}

/****************************************************************************/

END_UGDIM_NAMESPACE
