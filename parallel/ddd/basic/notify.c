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
#include <stdlib.h>
#include <stdio.h>

#include "dddi.h"
#include "basic/notify.h"


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
  DDD_PROC from, to;               /* source and destination processors */
  size_t size;                     /* message size */
  unsigned short flag;             /* one of NotifyTypes */
} NOTIFY_INFO;



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
    DDD_PrintError('E', 6301, "not enough memory in NotifyInit");
    HARD_EXIT;
  }


  maxInfos = MAX_INFOS;                 /* TODO maximum value, just for testing */


  /* init local array for all Info records */
  allInfoBuffer = (NOTIFY_INFO *) AllocFix(maxInfos*sizeof(NOTIFY_INFO));
  if (allInfoBuffer==NULL)
  {
    DDD_PrintError('E', 6300, "not enough memory in NotifyInit");
    HARD_EXIT;
  }
}


void NotifyExit (void)
{
  /* free memory */
  FreeFix(theRouting);
  FreeFix(allInfoBuffer);
}


/****************************************************************************/



static int sort_XferInfos (const void *e1, const void *e2)
{
  NOTIFY_INFO *ci1, *ci2;

  ci1 = (NOTIFY_INFO *)e1;
  ci2 = (NOTIFY_INFO *)e2;

  /* PROC_INVALID is less than all 'real' processor numbers */
  if (ci1->to==PROC_INVALID && ci2->to!=PROC_INVALID) return(-1);
  if (ci1->to!=PROC_INVALID && ci2->to==PROC_INVALID) return(1);

  if (ci1->to!=PROC_INVALID && ci2->to!=PROC_INVALID)
  {
    if (ci1->to < ci2->to) return(-1);
    if (ci1->to > ci2->to) return(1);
  }

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
  allInfos[0].to   = PROC_INVALID;
  allInfos[0].size = 0;
  allInfos[0].flag = DUMMY;
  lastInfo = 1;

  return(allInfos);
}



/****************************************************************************/

int NotifyTwoWave (NOTIFY_INFO *allInfos, int lastInfo)
{
  NOTIFY_INFO  *newInfos;
  int l, i, j, n, unknownInfos, myInfos;

#if     DebugNotify<=4
  printf("%4d:    NotifyTwoWave, lastInfo=%d\n", me, lastInfo);
  fflush(stdout);
#endif

  /* BOTTOM->TOP WAVE */
  /* get local Info lists from downtree */
  for(l=degree-1; l>=0; l--)
  {
    GetConcentrate(l, &n, sizeof(int));

    /* TODO in error case other processors will hang! */
    if (lastInfo+n >= maxInfos) {
      DDD_PrintError('E', 6321, "msg-info array overflow in NotifyTwoWave");
      return(ERROR);
    }
    GetConcentrate(l, &(allInfos[lastInfo]), n*sizeof(NOTIFY_INFO));

    /* construct routing table */
    for(i=0; i<n; i++)
      theRouting[allInfos[lastInfo+i].from] = l;

    lastInfo += n;
  }


  /* determine target direction in tree */
  /* TODO: eventually extra solution for root node!
       (it knows all flags are MYSELF or KNOWN!)  */
  qsort(allInfos, lastInfo, sizeof(NOTIFY_INFO), sort_XferInfos);
  i = j = 0;
  unknownInfos = lastInfo;
  myInfos = 0;
  while (i<lastInfo && allInfos[j].to==PROC_INVALID)
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
    printf("%4d:    NotifyTwoWave, send uptree unknown %d/%d (%d|%d;%d)\n",
           me, i, unknownInfos,
           newInfos[i].to, newInfos[i].from, newInfos[i].size);
  }
#endif

#if     DebugNotify<=3
  printf("%4d:    NotifyTwoWave, wave 1 ready\n", me);
  fflush(stdout);
#endif



  /* TOP->BOTTOM WAVE */

  /* get Infos local to my subtree from uptree */
  unknownInfos = 0;
  GetSpread(&unknownInfos, sizeof(int));

  if (unknownInfos>0)
    GetSpread(newInfos, unknownInfos*sizeof(NOTIFY_INFO));
  lastInfo += unknownInfos;

  /* sort Infos according to routing */
  qsort(allInfos, lastInfo, sizeof(NOTIFY_INFO), sort_XferRouting);

#if     DebugNotify<=1
  for(i=0; i<lastInfo; i++)
  {
    printf("%4d:    NotifyTwoWave, sorted for routing  %d/%d (%d|%d;%d)\n",
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
  printf("%4d:    NotifyTwoWave, wave 2 ready, nRecv=%d\n", me, lastInfo);
  fflush(stdout);
#endif

  return(lastInfo);
}



/****************************************************************************/


NOTIFY_DESC *DDD_NotifyBegin (int n)
{
  nSendDescs = n;
  if (procs>1)
  {
    theDescs = (NOTIFY_DESC *) AllocTmp(sizeof(NOTIFY_DESC)*(procs-1));
  }
  else
  {
    theDescs = NULL;
  }

  if (n>procs-1)
  {
    DDD_PrintError('E', 6340,
                   "more send-messages than other processors in DDD_NotifyAlloc");
    return(NULL);
  }

  return(theDescs);
}


void DDD_NotifyEnd (void)
{
  if (theDescs!=NULL)
  {
    FreeTmp(theDescs);
  }
}


int DDD_Notify (void)
{
  NOTIFY_INFO  *allInfos;
  int i, nRecvMsgs;

  /* get storage for local info list */
  allInfos = NotifyPrepare();
  if (allInfos==NULL) return(ERROR);

  /* convert message list to local Info list */
  for(i=0; i<nSendDescs; i++)
  {
#               if      DebugNotify<=4
    printf("%4d:    Notify send msg #%02d to %3d size=%d\n", me,
           lastInfo, theDescs[i].proc, theDescs[i].size);
#               endif

    allInfos[lastInfo].from = me;
    allInfos[lastInfo].to   = theDescs[i].proc;
    allInfos[lastInfo].size = theDescs[i].size;
    allInfos[lastInfo].flag = UNKNOWN;
    lastInfo++;
  }


  /* notify partners */
  nRecvMsgs = NotifyTwoWave(allInfos, lastInfo);

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
