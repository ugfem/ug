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

void DDD_Init (int *argcp, char ***argvp)
{
  int buffsize;

  /* init lineout-interface to stdout */
  DDD_UserLineOutFunction = NULL;

  /* init PPIF */
  InitPPIF(argcp, argvp);

  /*
     printf("%4d: process_id=%d\n", me, getpid());
   */

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
    return;
  }
  /* overlay with other buffers */
  cBuffer = (char *)iBuffer;

  /* init all DDD components */
  NotifyInit();
  LowCommInit();
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
  DDD_SetOption(OPT_WARNING_VARSIZE_OBJ, OPT_ON);
  DDD_SetOption(OPT_WARNING_SMALLSIZE, OPT_ON);
  DDD_SetOption(OPT_WARNING_PRIOCHANGE, OPT_ON);
  DDD_SetOption(OPT_WARNING_DESTRUCT_HDR, OPT_ON);
  DDD_SetOption(OPT_DEBUG_XFERMESGS, OPT_OFF);
  DDD_SetOption(OPT_QUIET_CONSCHECK, OPT_OFF);
  DDD_SetOption(OPT_IDENTIFY_MODE, IDMODE_LISTS);
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

void DDD_Exit (void)
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
  LowCommExit();
  NotifyExit();

  /* exit PPIF */
  ExitPPIF();
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

void DDD_Status (void)
{
  sprintf(cBuffer, "| DDD_Status for proc=%03d, DDD-Version %s\n", me,
          DDD_VERSION);
  DDD_PrintLine(cBuffer);
  sprintf(cBuffer, "|\n|     MAX_ELEMDESC = %3d\n", MAX_ELEMDESC);
  sprintf(cBuffer, "|     MAX_TYPEDESC = %3d\n", MAX_TYPEDESC);
  DDD_PrintLine(cBuffer);
  sprintf(cBuffer, "|\n|     MAX_OBJ = %8d  MAX_CPL = %8d\n", MAX_OBJ, MAX_CPL);
  DDD_PrintLine(cBuffer);

  sprintf(cBuffer, "|     nObjs   = %8d  nCpls   = %8d  nCplItems = %8d\n",
          nObjs, nCpls, nCplItems);
  DDD_PrintLine(cBuffer);
  DDD_PrintLine("|\n|     Timeouts:\n");
  sprintf(cBuffer, "|        IFComm:  %9d\n", MAX_TRIES);
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

void DDD_LineOutRegister (void (*func)(char *s))
{
  DDD_UserLineOutFunction = func;
}




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

void DDD_SetOption (DDD_OPTION option, int val)
{
  if (option<0 || option>=OPT_END)
  {
    DDD_PrintError('E', 1999, "invalid DDD_OPTION in DDD_SetOption()");
    return;
  }

  theOptions[option] = val;
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
  if (option<0 || option>=OPT_END)
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

DDD_PROC DDD_InfoMe (void)
{
  return me;
}


DDD_PROC DDD_InfoMaster (void)
{
  return master;
}


DDD_PROC DDD_InfoProcs (void)
{
  return procs;
}


/****************************************************************************/
