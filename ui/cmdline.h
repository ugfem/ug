// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  cmdline.h                                                                                                     */
/*																			*/
/* Purpose:   defines the command abstract data type						*/
/*																			*/
/* Author:	  Peter Bastian,												*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/*																			*/
/* History:   18.02.92 begin, ug version 2.0								*/
/*			  05 Sep 1992, split cmd.c into cmdint.c and commands.c                 */
/*			  17.12.94 ug 3.0												*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __CMDLINE__
#define __CMDLINE__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* return values for commands */
#define DONE                    0               /* command exectuted properly				*/
#define OKCODE                  0               /* command exectuted properly				*/
#define QUITCODE                12345   /* indicates end of program to calling if	*/
#define PARAMERRORCODE  3               /* not enough parameters					*/
#define CMDERRORCODE    4               /* other error condition					*/
#define INTERRUPTCODE   5               /* cmd terminated by user interrupt			*/
#define FATAL                   9999    /* fatal error, quit program				*/

#define MAXOPTIONS              100             /* maximum number of options for a command	*/
#define OPTIONBUFFERLEN 16384   /* length of option buffer					*/

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef INT (*CommandProcPtr)(INT,char **);

typedef struct {                                /* executable command variable				*/
  ENVVAR v;                                             /* this is an environment variable			*/
  CommandProcPtr cmdProc;               /* function to be executed					*/
} COMMAND ;

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT      InitCmdline            (void);

/* command creation and execution */
COMMAND *GetFirstCommand        (void);
COMMAND *GetNextCommand         (const COMMAND *cmd);
COMMAND *SearchUgCmd            (const char *cmdName);
COMMAND *CreateCommand          (const char *name, CommandProcPtr cmdProc);
COMMAND *GetCommand             (const char *name);
COMMAND *ReplaceCommand         (const char *name, CommandProcPtr cmdProc);
INT      ExecCommand            (char *cmdLine);

#endif
