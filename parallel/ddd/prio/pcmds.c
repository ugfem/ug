// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pcmds.c                                                       */
/*                                                                          */
/* Purpose:   DDD-commands for Prio Environment                             */
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
/* History:   980720 kb  begin                                              */
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


#define DebugPrio     10   /* 0 is all, 10 is off */


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


/* overall mode of prio-environment */
enum PrioMode {
  PMODE_IDLE = 0,            /* waiting for next DDD_PrioBegin() */
  PMODE_CMDS,                /* after DDD_PrioBegin(), before DDD_PrioEnd() */
  PMODE_BUSY                 /* during DDD_PrioEnd() */
};





/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/* PRIO_GLOBALS: global data for prio module                                */
/****************************************************************************/

typedef struct _PRIO_GLOBALS
{
  /* mode of prio module */
  int prioMode;

} PRIO_GLOBALS;






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



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/* one instance of PRIO_GLOBALS */
static PRIO_GLOBALS prioGlobals;



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/*
        management functions for PrioMode.

        these functions control the mode the prio-module is
        currently in. this is used for error detection, but
        also for correct detection of coupling inconsistencies
        and recovery.
 */

char *PrioModeName (int mode)
{
  switch(mode)
  {
  case PMODE_IDLE : return "idle-mode";
  case PMODE_CMDS : return "commands-mode";
  case PMODE_BUSY : return "busy-mode";
  }
  return "unknown-mode";
}


static void PrioSetMode (int mode)
{
  prioGlobals.prioMode = mode;

#       if DebugPrio<=8
  sprintf(cBuffer, "%4d: PrioMode=%s.\n",
          me, PrioModeName(prioGlobals.prioMode));
  DDD_PrintDebug(cBuffer);
#       endif
}


static int PrioSuccMode (int mode)
{
  switch(mode)
  {
  case PMODE_IDLE : return PMODE_CMDS;
  case PMODE_CMDS : return PMODE_BUSY;
  case PMODE_BUSY : return PMODE_IDLE;
  }
  return PMODE_IDLE;
}



static int PrioMode (void)
{
  return prioGlobals.prioMode;
}


int ddd_PrioActive (void)
{
  return prioGlobals.prioMode!=PMODE_IDLE;
}


static int PrioStepMode (int old)
{
  if (prioGlobals.prioMode!=old)
  {
    sprintf(cBuffer, "wrong prio-mode (currently in %s, expected %s)",
            PrioModeName(prioGlobals.prioMode), PrioModeName(old));
    DDD_PrintError('E', 8200, cBuffer);
    return FALSE;
  }

  PrioSetMode(PrioSuccMode(prioGlobals.prioMode));
  return TRUE;
}


/****************************************************************************/


void ddd_PrioInit (void)
{
  PrioSetMode(PMODE_IDLE);
}


void ddd_PrioExit (void)
{}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_PrioChange                                                */
/*                                                                          */
/****************************************************************************/

/**
        Consistent change of a local object's priority during DDD Prio Environment.
        Local objects which are part of a distributed object must notify
        other copies about local priority changes.
        DDD will send appropriate messages to the owner processors of
        the other copies.

        This function is regarded as a {\bf Prio}-operation due
        to its influence on DDD management information on neighbouring
        processors. Therefore the function has to be issued between
        a starting \funk{PrioBegin} and a final \funk{PrioEnd} call.

   @param hdr  DDD local object whose priority should be changed.
   @param prio new priority for this local object.
 */

void DDD_PrioChange (DDD_HDR hdr, DDD_PRIO prio)
{
  DDD_PRIO old_prio = OBJ_PRIO(hdr);


  if (!ddd_PrioActive())
  {
    DDD_PrintError('E', 8030, "Missing DDD_PrioBegin(). aborted");
    HARD_EXIT;
  }


  /* change priority of object directly, for local objects this
     is all what we need. */
  {
    DDD_PRIO newprio;
    PriorityMerge(&theTypeDefs[OBJ_TYPE(hdr)],
                  OBJ_PRIO(hdr), prio, &newprio);
    OBJ_PRIO(hdr) = newprio;
  }

  /* handle distributed objects
     if (ObjHasCpl(hdr))
     {
          nothing to do here:
          for distributed objects, we will communicate the prio
          via standard interface later.
     }
   */


#       if DebugPrio<=2
  sprintf(cBuffer, "%4d: DDD_PrioChange %08x, old_prio=%d. new_prio=%d\n",
          me, OBJ_GID(hdr), old_prio, OBJ_PRIO(hdr));
  DDD_PrintDebug(cBuffer);
#       endif
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_PrioEnd                                                   */
/*                                                                          */
/****************************************************************************/

static int GatherPrio (DDD_HDR obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  *((DDD_PRIO *)data) = OBJ_PRIO(obj);
  return(0);
}

static int ScatterPrio (DDD_HDR obj, void *data, DDD_PROC proc, DDD_PRIO prio)
{
  DDD_PRIO real_prio = *((DDD_PRIO *)data);

  /* if prio on other processor has been changed, we adapt it here. */
  if (real_prio!=prio)
    ModCoupling(obj, proc, real_prio);

  return(0);
}


/**
        End of PrioEnvironment.
        This function starts the actual process of changing priorities.
        After a call to this function (on all processors) all
        \funk{PrioChange}-commands since the last call to \funk{PrioBegin}
        are executed. This involves a set of interface communications
        between the processors.
 */

#ifdef C_FRONTEND
DDD_RET DDD_PrioEnd (void)
#endif
#ifdef CPP_FRONTEND
DDD_RET DDD_Library::PrioEnd (void)
#endif
#ifdef F_FRONTEND
DDD_RET DDD_PrioEnd (void)
#endif
{
  /* step mode and check whether call to PrioEnd is valid */
  if (!PrioStepMode(PMODE_CMDS))
  {
    DDD_PrintError('E', 8011, "DDD_PrioEnd() aborted");
    HARD_EXIT;
  }


  ddd_StdIFExchangeX(sizeof(DDD_PRIO), GatherPrio, ScatterPrio);

  /*
          free temporary storage
   */
  STAT_RESET;
  IFAllFromScratch();
  STAT_TIMER(T_PRIO_BUILD_IF);


  PrioStepMode(PMODE_BUSY);

  return(DDD_RET_OK);
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_PrioBegin                                                 */
/*                                                                          */
/****************************************************************************/

/**
        Starts a PrioEnvironment.
        A call to this function establishes a global operation of changing
        priorities.  It must be issued on all processors. After this call an
        arbitrary series of \funk{PrioChange}-commands may be issued. The
        global transfer operation is carried out via a \funk{PrioEnd} call on
        each processor.
 */

#ifdef C_FRONTEND
void DDD_PrioBegin (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::PrioBegin (void)
#endif
#ifdef F_FRONTEND
void DDD_PrioBegin (void)
#endif
{
  /* step mode and check whether call to JoinBegin is valid */
  if (!PrioStepMode(PMODE_IDLE))
  {
    DDD_PrintError('E', 8010, "DDD_PrioBegin() aborted");
    HARD_EXIT;
  }
}



/****************************************************************************/

#undef _FADR

/****************************************************************************/
