// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      join.c                                                        */
/*                                                                          */
/* Purpose:   main module for object join                                   */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: birken@ica3.uni-stuttgart.de                           */
/*            phone: 0049-(0)711-685-7007                                   */
/*            fax  : 0049-(0)711-685-7000                                   */
/*                                                                          */
/* History:   980122 kb  begin                                              */
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



/*
        NOTE: all container-classes from ooppcc.h are implemented in this
              source file by setting the following define.
 */
#define ContainerImplementation
#define _CHECKALLOC(ptr)   assert(ptr!=NULL)


static int TmpMem_kind = TMEM_ANY;

static void *join_AllocTmp (size_t size)
{
  return AllocTmpReq(size, TmpMem_kind);
}

static void join_FreeTmp (void *buffer)
{
  FreeTmpReq(buffer, 0, TmpMem_kind);
}

void join_SetTmpMem (int kind)
{
  TmpMem_kind = kind;
}



#include "join.h"



/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/

/* helpful macros for FRONTEND switching, will be #undef'd at EOF */
#ifdef F_FRONTEND
#define _FADR     &
#else
#define _FADR
#endif


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



JOIN_GLOBALS joinGlobals;


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* class member function implementations                                    */
/*                                                                          */
/****************************************************************************/

#define ClassName JIJoin

/*
        compare-method in order to eliminate double JIJoin-items.

        the items are sorted according to key (dest,remote_gid),
        all in ascending order.
 */
int Method(Compare) (ClassPtr item1, ClassPtr item2)
{
  int ret;

  if (item1->dest < item2->dest) return(-1);
  if (item1->dest > item2->dest) return(1);

  if (item1->new_gid < item2->new_gid) return(-1);
  if (item1->new_gid > item2->new_gid) return(1);


  /* items have equal gid and dest, so they are considered as equal. */
  return(0);
}


void Method(Print) (ParamThis _PRINTPARAMS)
{
  fprintf(fp, "JIJoin local_gid=%08x dest=%d new_gid=%08x\n",
          OBJ_GID(This->hdr), This->dest, This->new_gid);
}

#undef ClassName




#define ClassName JIAddCpl

/*
        compare-method in order to eliminate double JIAddCpl-items.

        the items are sorted according to key (dest,gid),
        all in ascending order.
 */
int Method(Compare) (ClassPtr item1, ClassPtr item2)
{
  int ret;

  if (item1->dest < item2->dest) return(-1);
  if (item1->dest > item2->dest) return(1);

  if (item1->te.gid < item2->te.gid) return(-1);
  if (item1->te.gid > item2->te.gid) return(1);

  if (item1->te.proc < item2->te.proc) return(-1);
  if (item1->te.proc > item2->te.proc) return(1);


  /* items have equal gid and dest, so they are considered as equal. */
  return(0);
}


void Method(Print) (ParamThis _PRINTPARAMS)
{
  fprintf(fp, "JIAddCpl gid=%08x dest=%d proc=%d prio=%d\n",
          This->te.gid, This->dest, This->te.proc, This->te.prio);
}

#undef ClassName




/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/*
        management functions for JoinMode.

        these functions control the mode the join-module is
        currently in. this is used for error detection, but
        also for correct detection of coupling inconsistencies
        and recovery.
 */

char *JoinModeName (int mode)
{
  switch(mode)
  {
  case JMODE_IDLE : return "idle-mode";
  case JMODE_CMDS : return "commands-mode";
  case JMODE_BUSY : return "busy-mode";
  }
  return "unknown-mode";
}


static void JoinSetMode (int mode)
{
  joinGlobals.joinMode = mode;

#       if DebugJoin<=8
  sprintf(cBuffer, "%4d: JoinMode=%s.\n",
          me, JoinModeName(joinGlobals.joinMode));
  DDD_PrintDebug(cBuffer);
#       endif
}


static int JoinSuccMode (int mode)
{
  switch(mode)
  {
  case JMODE_IDLE : return JMODE_CMDS;
  case JMODE_CMDS : return JMODE_BUSY;
  case JMODE_BUSY : return JMODE_IDLE;
  }
  return JMODE_IDLE;
}



int JoinMode (void)
{
  return joinGlobals.joinMode;
}


int ddd_JoinActive (void)
{
  return joinGlobals.joinMode!=JMODE_IDLE;
}


int JoinStepMode (int old)
{
  if (joinGlobals.joinMode!=old)
  {
    sprintf(cBuffer, "wrong join-mode (currently in %s, expected %s)",
            JoinModeName(joinGlobals.joinMode), JoinModeName(old));
    DDD_PrintError('E', 7200, cBuffer);
    return FALSE;
  }

  JoinSetMode(JoinSuccMode(joinGlobals.joinMode));
  return TRUE;
}


/****************************************************************************/


void ddd_JoinInit (void)
{
  /* set kind of TMEM alloc/free requests */
  join_SetTmpMem(TMEM_ANY);

  /* init control structures for JoinInfo-items in messages */
  joinGlobals.setJIJoin    = New_JIJoinSet();
  joinGlobals.setJIAddCpl2 = New_JIAddCplSet();
  joinGlobals.setJIAddCpl3 = New_JIAddCplSet();

  JoinSetMode(JMODE_IDLE);

  joinGlobals.phase1msg_t = LC_NewMsgType("Join1Msg");
  joinGlobals.jointab_id = LC_NewMsgTable("GidTab",
                                          joinGlobals.phase1msg_t, sizeof(TEJoin));

  joinGlobals.phase2msg_t = LC_NewMsgType("Join2Msg");
  joinGlobals.addtab_id = LC_NewMsgTable("AddCplTab",
                                         joinGlobals.phase2msg_t, sizeof(TEAddCpl));

  joinGlobals.phase3msg_t = LC_NewMsgType("Join3Msg");
  joinGlobals.cpltab_id = LC_NewMsgTable("AddCplTab",
                                         joinGlobals.phase3msg_t, sizeof(TEAddCpl));
}


void ddd_JoinExit (void)
{
  /* set kind of TMEM alloc/free requests */
  join_SetTmpMem(TMEM_ANY);


  /* TODO data (e.g., lists&trees of JI-items) should be freed!! */
}



/****************************************************************************/

#undef _FADR


/****************************************************************************/
