// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      stat.c                                                        */
/*                                                                          */
/* Purpose:   DDD statistical evaluation                                    */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   95/01/16 kb  begin                                            */
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
#include <stdarg.h>

#include "dddi.h"


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
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/* Storage for statistical data (compile with option -DStatistics)
   Caution: statistics should be used via macros in dddi.h! */

#ifdef Statistics
STATISTICS statistics;
#endif


char *DDD_StatXfer[] =
{
  /* DDD_XferEnd */
  /* 00 */ "XferThirdParty",
  /* 01 */ "Create sorted XferItem-array",
  /* 02 */ "XferElimDoubleXferInfos",
  /* 03 */ "Get array of new objects",
  /* 04 */ "XferBuildMsgInfos",
  /* 05 */ "XferNotify",
  /* 06 */ "Get comm.-connections",
  /* 07 */ "XferInitRecv",
  /* 08 */ "XferPackMsgs",
  /* 09 */ "XferDeleteObjects",
  /* 10 */ "XferReceiveMsgs/XferPollSend Loop",
  /* 11 */ "Free TMEM",
  /* 12 */ "IFAllFromScratch",

  /* XferReceiveMsgs */
  /* 13 */ "XUnp: PollReceiveMsgs",
  /* 14 */ "XUnp: XferUnpack",

  /* 15 */ "XferPollSend",
  /* 16 */ "XPac: Allocate send-buffers",
  /* 17 */ "XPac: Pack all messages",
  /* 18 */ "Xfer",
  /* 19 */ "Xfer",
  /* 20 */ "XUSM",
  /* 21 */ "XUSM",
  /* 22 */ "XUSM",
  /* 23 */ "XUSM",
  /* 24 */ "XUSM",
  /* 25 */ "XUSM",
  /* 26 */ "XUSM",
  /* 27 */ "XUnp: qsort GlobalTabs",
  /* 28 */ "XUnp: XferUnifyDel",
  /* 29 */ "XUnp: XferUnifyCpl",
  /* 30 */ "XUnp: XferDelLocalCpl",
  /* 31 */ "XUnp: Unpack all messages",

  /* XferPackSingleMsg */
  /* 32 */ "XPSM: Copy object into message",
  /* 33 */ "XPSM: BuildSymTab",
  /* 34 */ "XPSM: qsort SymTab",
  /* 35 */ "XPSM: qsort ObjTab",
  /* 36 */ "XPSM: qsort CplTab",
  /* 37 */ "XPSM: Replace ptrs by SymIndex",

  /* 38 */ "Xfer",
  /* 39 */ "Xfer",
  /* 40 */ "Xfer",
  /* 41 */ "Xfer",
  /* 42 */ "Xfer",
  /* 43 */ "Xfer",
  /* 44 */ "Xfer",
  /* 45 */ "Xfer",
  /* 46 */ "Xfer",
  /* 47 */ "Xfer",
  /* 48 */ "Xfer",
  /* 49 */ "Xfer",
  /* 50 */ "Xfer",
  /* 51 */ "Xfer",
  /* 52 */ "Xfer",
  /* 53 */ "Xfer",
  /* 54 */ "Xfer",
  /* 55 */ "Xfer",
  /* 56 */ "Xfer",
  /* 57 */ "Xfer",
  /* 58 */ "Xfer",
  /* 59 */ "Xfer",
  /* 60 */ "IF: overall communication",
  /* 61 */ "IF: first phase",
  /* 62 */ "IF: second phase",
  /* 63 */ "IF: third phase"
};


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


double DDD_StatClock (int index)
{
  return(STAT_GETTIMER(index));
}

int DDD_StatCount (int index)
{
  return(STAT_GETCOUNT(index));
}
