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

/* standard C library  */
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

DDD_HDR   *ddd_ObjTable;
int ddd_ObjTabSize;
int ddd_nObjs;

COUPLING **ddd_CplTable;
short     *ddd_NCplTable;
int ddd_CplTabSize;
int ddd_nCpls;
int nCplItems;


int        *iBuffer;        /* general bufferspace, integer */
char       *cBuffer;        /* general bufferspace, integer */

int theOptions[OPT_END];


#define ddd_SetOption(o,v)      theOptions[o]=(v)


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


static void *LowComm_DefaultAlloc (size_t s)
{
  return memmgr_AllocTMEM(s, TMEM_LOWCOMM);
}

static void LowComm_DefaultFree (void *buffer)
{
  memmgr_FreeTMEM(buffer, TMEM_LOWCOMM);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_Init                                                      */
/*                                                                          */
/****************************************************************************/

/**
        Initialisation of the DDD library.
        This function has to be called before any other function
        of the DDD library is called. It initializes the underlying
        PPIF-library, sets all DDD options to their default values
        and initiates all DDD subsystems.

        As some of the memory handler calls will be initiated during
        the execution of this function, the memory manager has to be
        initialized before calling \funk{Init}.

        Note: Not the actual arguments of the {\em main()}-function
        are passed here (\ie, {\em argc} and {\em argv}), but {\em pointers}
        to these arguments (\ie, {\em argcp} and {\em argvp}). This is
        due to manipulation of the argument lists by certain underlying
        message-passing-implementations (\eg, Argonne MPI).
        The parameters are not needed for initialisation of the DDD
        library itself, but will be passed to the PPIF-libraries
        initialisation function \ppiffunk{InitPPIF}. If the first
        parameter {\em argcp} is #NULL#, DDD assumes that the PPIF-library
        has been initialized explicitly by the application program.

   @param  argcp      pointer to argc (the application's parameter count)
   @param  argvp      pointer to argv (the application's parameter list)
 */

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
  LC_Init(LowComm_DefaultAlloc, LowComm_DefaultFree);
  ddd_StatInit();
  ddd_TypeMgrInit();
  ddd_ObjMgrInit();
  ddd_CplMgrInit();
  ddd_TopoInit();
  ddd_IdentInit();
  ddd_IFInit();
  ddd_XferInit();
  ddd_JoinInit();
  ddd_ConsInit();

  /* reset all global counters */
  ddd_nObjs  = 0;
  NCpl_Init;
  nCplItems  = 0;

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
  ddd_SetOption(OPT_INFO_JOIN,             JOIN_SHOW_NONE);
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
/****************************************************************************/

/**
        Clean-up of the DDD library.
        This function frees memory previously allocated by DDD and finally
        finishes up the PPIF library. After the call to \funk{Exit}
        further usage of the DDD library is no longer possible during this program
        run.

        The clean-up of the memory manager should happen afterwards and is left
        to the DDD application programmer.
 */

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
  ddd_JoinExit();
  ddd_XferExit();
  ddd_IFExit();
  ddd_IdentExit();
  ddd_TopoExit();
  ddd_CplMgrExit();
  ddd_ObjMgrExit();
  ddd_TypeMgrExit();
  ddd_StatExit();
  LC_Exit();
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
/****************************************************************************/

/**
        Show global status information.
        This function displays information concerning both
        the compile-time parameters of the DDD-library and some important
        runtime-variables. Overview of compile time parameters that will
        be displayed:

        \begin{tabular}{l|l}
        Parameter       & Description\\ \hline
        DDD-Version     & current library version number\\ ##
   #MAX_TYPEDESC#  & maximum number of #DDD_TYPE# IDs\\ ##
   #MAX_OBJ#       & maximum number of DDD-objects on one processor\\ ##
   #MAX_CPL#       & maximum number of couplings on one processor\\ ##
        \end{tabular}
 */

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
#ifdef WithFullObjectTable
  sprintf(cBuffer, "|\n|     MAX_OBJ = %8d  MAX_CPL = %8d\n",
          ddd_ObjTabSize, ddd_CplTabSize);
#else
  sprintf(cBuffer, "|\n|     MAX_CPL = %8d\n", ddd_CplTabSize);
#endif
  DDD_PrintLine(cBuffer);

  sprintf(cBuffer, "|     nObjs   = %8d  nCpls   = %8d  nCplItems = %8d\n",
          ddd_nObjs, NCpl_Get, nCplItems);
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
/****************************************************************************/

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
/**
        Redirect text output.
        This function sets the DDD-textport to a given handler function.
        The handler should be declared as follows:

   #void func(char *line_of_text)#

        Instead of printing text for error, debugging and info messages
        directly to {\em standard output}, DDD will redirect all output
        one line at a time and send it to the handler {\em func}.
        This can be used to send each processor's output into a separate file.

   @param  func  handler function which should be used for text redirection
 */
#endif

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
/****************************************************************************/

/**
        Set a DDD-option to a given value.
        The current behaviour of the DDD library can be configured
        at runtime by setting a variety of options to given values.
        For each option, there is a default setting and a set of
        possible values. See \Ref{DDD Options} for a description
        of all possible options with their default settings and
        meaning.

   @param option   DDD option specifier
   @param value    option value, possible values depend on option specifier
 */

#ifdef C_FRONTEND
void DDD_SetOption (DDD_OPTION option, int value)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Library::SetOption (DDD_OPTION option, int value)
{
#endif
#ifdef F_FRONTEND
void DDD_SetOption (DDD_OPTION *_option, int *_value)
{
  DDD_OPTION option = *_option;
  int value = *_value;
#endif
if (option>=OPT_END)
{
  DDD_PrintError('E', 1999, "invalid DDD_OPTION in DDD_SetOption()");
  return;
}

ddd_SetOption(option, value);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_GetOption (not exported)                                  */
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

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
/**
        Get local processor number.
        This function returns the local processor number, which is an
        integer number between $0$ and \funk{InfoProcs}.

   @return  local processor number
 */
#endif

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



#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
/**
        Get master processor number.
        This function returns the processor number of the
        master processor. Processor numbers are always
        integer numbers between 0 and \funk{InfoProcs}.
        Usually, processor 0 is the master processor.

   @return  master processor number
 */
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


#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
/**
        Get total number of processors.

   @return  total number of processors
 */
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
