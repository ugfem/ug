// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ddd.c                                                         */
/*                                                                          */
/* Purpose:   distributed dynamic data module                               */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/22 kb  begin                                            */
/*            94/09/13 kb  added DDD_Status                                 */
/*            95/11/03 kb  complete rewrite of StructRegister code          */
/*            95/11/16 kb  moved DDD_TYPE-definition code to typemgr.c      */
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
#include <string.h>

#include "dddi.h"
#include "basic/notify.h"
#include "basic/lowcomm.h"



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define  BUFFER_SIZE_FACTOR   3
#define  MIN_BUFFER_SIZE     256


/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

DDD_HDR theObj[MAX_OBJ];
int nObjs;

COUPLING   *theCpl[MAX_CPL];
int theCplN[MAX_CPL];
int nCpls;
int nCplItems;

int theIdCount;             /* local unique ID count */

int        *iBuffer;        /* general bufferspace, integer */
char       *cBuffer;        /* general bufferspace, integer */

int theOptions[OPT_END];


#define ddd_SetOption(o,v)      theOptions[o]=(v)


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_Init                                                      */
/*                                                                          */
/* Purpose:   initialize the DDD library                                    */
/*                                                                          */
/* Input:     argcp: pointer to argc (the applications parameter count)     */
/*            argvp: pointer to argv (the applications parameter list)      */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void DDD_Init (int *argcp, char ***argvp)
#endif
#ifdef CPP_FRONTEND
DDD_Library::DDD_Library (int *argcp, char ***argvp)
#endif
{
#ifdef CPP_FRONTEND
  // check existence of another instance of DDD_Library
  if (_instance!=0)
  {
    DDD_PrintError('E', 1021,
                   "construction of two instances of DDD_Library is not allowed");
    HARD_EXIT;
  }
#endif

  int buffsize;

  /* init lineout-interface to stdout */
  DDD_UserLineOutFunction = NULL;

  /* if first arg is NULL, we assume that PPIF has been initialized elsewhere */
  if (argcp!=NULL)
  {
    /* init PPIF */
    if (InitPPIF(argcp, argvp) != PPIF_SUCCESS)
    {
      DDD_PrintError('E', 1005, "PPIF initialization failed");
      HARD_EXIT;
    }
  }


  /*
     printf("%4d: process_id=%d\n", me, getpid());
   */


  /* check max. number of procs (limited by GID construction) */
  if (procs>MAX_PROCS) {
    DDD_PrintError('E', 1010,
                   "too many processors, cannot construct global IDs in DDD_Init");
    HARD_EXIT;
  }

  /* compute size for general buffer */
  buffsize = (procs+1)*(sizeof(int)*BUFFER_SIZE_FACTOR);
  if (buffsize<MIN_BUFFER_SIZE)
  {
    buffsize = MIN_BUFFER_SIZE;
  }

  /* get bufferspace */
  iBuffer = (int *)AllocFix(buffsize);
  if (iBuffer==NULL) {
    DDD_PrintError('E', 1000, "not enough memory in DDD_Init");
    HARD_EXIT;
  }
  /* overlay with other buffers */
  cBuffer = (char *)iBuffer;

  /* init all DDD components */
  NotifyInit();
  LowCommInit();
  ddd_StatInit();
  ddd_TypeMgrInit();
  ddd_CplMgrInit();
  ddd_TopoInit();
  ddd_IdentInit();
  ddd_IFInit();
  ddd_XferInit();
  ddd_ConsInit();

  /* reset all global counters */
  nObjs  = 0;
  nCpls  = 0;
  nCplItems  = 0;
  theIdCount = 1;        /* start with 1, for debugging reasons */

  /* set options on default values */
  ddd_SetOption(OPT_WARNING_VARSIZE_OBJ,   OPT_ON);
  ddd_SetOption(OPT_WARNING_SMALLSIZE,     OPT_ON);
  ddd_SetOption(OPT_WARNING_PRIOCHANGE,    OPT_ON);
  ddd_SetOption(OPT_WARNING_DESTRUCT_HDR,  OPT_ON);
  ddd_SetOption(OPT_DEBUG_XFERMESGS,       OPT_OFF);
  ddd_SetOption(OPT_QUIET_CONSCHECK,       OPT_OFF);
  ddd_SetOption(OPT_IDENTIFY_MODE,         IDMODE_LISTS);
  ddd_SetOption(OPT_WARNING_REF_COLLISION, OPT_ON);
  ddd_SetOption(OPT_INFO_XFER,             XFER_SHOW_NONE);
  ddd_SetOption(OPT_WARNING_OLDSTYLE,      OPT_ON);
  ddd_SetOption(OPT_INFO_IF_WITH_ATTR,     OPT_OFF);
  ddd_SetOption(OPT_XFER_PRUNE_DELETE,     OPT_OFF);
  ddd_SetOption(OPT_IF_REUSE_BUFFERS,      OPT_OFF);
  ddd_SetOption(OPT_IF_CREATE_EXPLICIT,    OPT_OFF);

#ifdef CPP_FRONTEND
  // remember pointer to singleton
  _instance = this;
#endif
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_Exit                                                      */
/*                                                                          */
/* Purpose:   exit from DDD library                                         */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void DDD_Exit (void)
#endif
#ifdef CPP_FRONTEND
DDD_Library::~DDD_Library (void)
#endif
{
  /* free bufferspace */
  FreeFix(iBuffer);

  /* close up all DDD components */
  ddd_ConsExit();
  ddd_XferExit();
  ddd_IFExit();
  ddd_IdentExit();
  ddd_TopoExit();
  ddd_CplMgrExit();
  ddd_TypeMgrExit();
  ddd_StatExit();
  LowCommExit();
  NotifyExit();

  /* exit PPIF */
  ExitPPIF();

#ifdef CPP_FRONTEND
  _instance = 0;
#endif
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_Status                                                    */
/*                                                                          */
/* Purpose:   print out background information about DDD library status     */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void DDD_Status (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::Status (void)
#endif
{
  sprintf(cBuffer, "| DDD_Status for proc=%03d, DDD-Version %s\n", me,
          DDD_VERSION);
  DDD_PrintLine(cBuffer);
  sprintf(cBuffer, "|\n|     MAX_ELEMDESC = %4d\n", MAX_ELEMDESC);
  sprintf(cBuffer, "|     MAX_TYPEDESC = %4d\n", MAX_TYPEDESC);
  sprintf(cBuffer, "|     MAX_PROCS    = %4d\n", MAX_PROCS);
  sprintf(cBuffer, "|     MAX_PRIO     = %4d\n", MAX_PRIO);
  DDD_PrintLine(cBuffer);
  sprintf(cBuffer, "|\n|     MAX_OBJ = %8d  MAX_CPL = %8d\n", MAX_OBJ, MAX_CPL);
  DDD_PrintLine(cBuffer);

  sprintf(cBuffer, "|     nObjs   = %8d  nCpls   = %8d  nCplItems = %8d\n",
          nObjs, nCpls, nCplItems);
  DDD_PrintLine(cBuffer);
  DDD_PrintLine("|\n|     Timeouts:\n");
  sprintf(cBuffer, "|        IFComm:  %12ld\n", (unsigned long)MAX_TRIES);
  DDD_PrintLine(cBuffer);

  sprintf(cBuffer, "|\n|     Compile-Time Options: ");

#       ifdef Statistics
  strcat(cBuffer, "Statistics ");
#       endif

  strcat(cBuffer, "\n");
  DDD_PrintLine(cBuffer);
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_LineOutRegister                                           */
/*                                                                          */
/* Purpose:   redirect DDD text output                                      */
/*                                                                          */
/* Input:     func: function which takes output string as an argument       */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_LineOutRegister (void (*func)(char *s))
#endif
#ifdef CPP_FRONTEND
void DDD_Library::LineOutRegister (void (*func)(char *))
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
{
  DDD_UserLineOutFunction = func;
}
#endif




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_SetOption                                                 */
/*                                                                          */
/* Purpose:   set DDD runtime options                                       */
/*                                                                          */
/* Input:     option:  OptionType of option to be set                       */
/*            val:     value of option                                      */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_SetOption (DDD_OPTION option, int val)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Library::SetOption (DDD_OPTION option, int val)
{
#endif
#ifdef F_FRONTEND
void DDD_SetOption (DDD_OPTION *_option, int *_val)
{
  DDD_OPTION option = *_option;
  int val = *_val;
#endif
if (option>=OPT_END)
{
  DDD_PrintError('E', 1999, "invalid DDD_OPTION in DDD_SetOption()");
  return;
}

ddd_SetOption(option, val);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_GetOption                                                 */
/*                                                                          */
/* Purpose:   get DDD runtime options                                       */
/*                                                                          */
/* Input:     option:  OptionType of option to get                          */
/*                                                                          */
/* Output:    value of option                                               */
/*                                                                          */
/****************************************************************************/

int DDD_GetOption (DDD_OPTION option)
{
  if (option>=OPT_END)
  {
    DDD_PrintError('E', 1999, "invalid DDD_OPTION in DDD_GetOption()");
    return 0;
  }

  return theOptions[option];
}

/****************************************************************************/

/*
        transparent access to global variables from PPIF
 */

#ifdef C_FRONTEND
DDD_PROC DDD_InfoMe (void)
{
  return me;
}
#endif
#ifdef CPP_FRONTEND
DDD_PROC DDD_Library::InfoMe (void)
{
  return me;
}
#endif


#ifdef C_FRONTEND
DDD_PROC DDD_InfoMaster (void)
{
  return master;
}
#endif
#ifdef CPP_FRONTEND
DDD_PROC DDD_Library::InfoMaster (void)
{
  return master;
}
#endif


#ifdef C_FRONTEND
DDD_PROC DDD_InfoProcs (void)
{
  return procs;
}
#endif
#ifdef CPP_FRONTEND
DDD_PROC DDD_Library::InfoProcs (void)
{
  return procs;
}
#endif


/****************************************************************************/

#ifdef CPP_FRONTEND
/*
        implementation for DDD_Library class
 */


// pointer to single instance of DDD_Library
DDD_Library* DDD_Library::_instance = 0;


DDD_Library* DDD_Library::Instance (void)
{
  if (_instance==0)
  {
    DDD_PrintError('E', 1020, "no instance of DDD_Library exists");
    HARD_EXIT;
  }

  return _instance;
}


#endif


/****************************************************************************/
