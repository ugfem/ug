// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  cmdline.c                                                                                                     */
/*																			*/
/* Purpose:   command structure and execution								*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   17.12.94 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "compiler.h"
#include "debug.h"
#include "misc.h"
#include "ugenv.h"
#include "cmdline.h"
#include "devices.h"
#include "general.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define OPTIONSEP               "$"     /* option character                                             */

#define WHITESPACE              " \t\n" /* for skipping of blanks					*/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT theMenuDirID;                        /* env type for Menu dir				*/
static INT theCommandVarID;             /* env type for Command vars			*/

/* for ExecCommand */
static char optionBuffer[OPTIONBUFFERLEN]; /* buffer to assemble options	*/
static char *options[MAXOPTIONS];               /* array of pointers to strings         */
static INT optionCount=0;                               /* number of options incl. cmd name */

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   CreateCommand - Register the commands exported by this module in the environ

   SYNOPSIS:
   COMMAND *CreateCommand (const char *name, CommandProcPtr cmdProc);

   PARAMETERS:
   .  name - Name of the command
   .  cmdProc - Pointer to a function of type 'CommandProcPtr'

   DESCRIPTION:
   This function registers a new command that can be executed from the UG
   shell. This process is described in detail on the page 'commands'.

   SEE ALSO:
   'commands'.

   RETURN VALUE:
   COMMAND *
   .n   pointer to new 'COMMAND' structure if o.k.
   .n   NULL pointer in case of an error.
   D*/
/****************************************************************************/

COMMAND *CreateCommand (const char *name, CommandProcPtr cmdProc)
{
  COMMAND *newCommand;

  /* change to Menu directory */
  if (ChangeEnvDir("/Menu")==NULL)
    return (NULL);

  /* allocate structure */
  newCommand = (COMMAND *) MakeEnvItem(name,theCommandVarID,sizeof(COMMAND));
  if (newCommand==NULL) return(NULL);

  /* fill data */
  newCommand->cmdProc = cmdProc;

  return(newCommand);
}

/****************************************************************************/
/*D
   GetCommand - Return pointer to command structure with name name

   SYNOPSIS:
   COMMAND *GetCommand (const char *name);

   PARAMETERS:
   .  name - Name of the command.

   DESCRIPTION:
   This function returns pointer to command structure with name 'name'. If
   it does not exist a 'NULL' pointer is returned.

   RETURN VALUE:
   COMMAND *
   .n      POINTER to 'COMMAND' structure
   .n      NULL if not found
   D*/
/****************************************************************************/

COMMAND *GetCommand (const char *name)
{
  COMMAND *cmd;

  if (ChangeEnvDir("/Menu")==NULL) return(NULL);
  cmd = (COMMAND *) SearchEnv(name,".",theCommandVarID,theMenuDirID);
  return(cmd);
}

/****************************************************************************/
/*D
   GetFirstCommand - Return pointer to first command structure of /Menu dir

   SYNOPSIS:
   COMMAND *GetFirstCommand ();

   PARAMETERS:
   .  void - no arguments

   DESCRIPTION:
   This function returns pointer to first command structure of the '/Menu'
   environment directory.

   RETURN VALUE:
   COMMAND *
   .n      POINTER to 'COMMAND' structure
   .n     NULL if not found
   D*/
/****************************************************************************/

COMMAND *GetFirstCommand ()
{
  ENVITEM *cmd;

  if ((cmd=(ENVITEM*)ChangeEnvDir("/Menu")) == NULL) return (NULL);

  for (cmd=ENVITEM_DOWN(cmd); cmd!=NULL; cmd=NEXT_ENVITEM(cmd))
    if (ENVITEM_TYPE(cmd) == theCommandVarID)
      return ((COMMAND*)cmd);
  return (NULL);
}

/****************************************************************************/
/*D
   GetNextCommand - Return pointer to command structure of /Menu dir following cmd

   SYNOPSIS:
   COMMAND *GetNextCommand (const COMMAND *cmd);

   PARAMETERS:
   .  cmd - pointer to a 'COMMAND'

   DESCRIPTION:
   This function returns a pointer to the next command after 'cmd' in the
   list of commands.

   RETURN VALUE:
   COMMAND *
   .n      POINTER to
   .n      NULL if not found
   D*/
/****************************************************************************/

COMMAND *GetNextCommand (const COMMAND *cmd)
{
  ENVITEM *nextCmd;

  for (nextCmd=NEXT_ENVITEM(cmd); nextCmd!=NULL; nextCmd=NEXT_ENVITEM(nextCmd))
    if (ENVITEM_TYPE(nextCmd) == theCommandVarID)
      return ((COMMAND*)nextCmd);
  return (NULL);
}
/**********************************************************************************/
/*
   Str1inStr2 -

   SYNOPSIS:
   static int Str1inStr2 (const char *name1, const char *name2);

   PARAMETERS:
   .  name1 -
   .  name2 -

   DESCRIPTION:
   This function

   RETURN VALUE
   int
   .n     0 if ok
   .n     1 if error occured.
 */
/***************************************************************************/
static int Str1inStr2 (const char *name1, const char *name2)
{
  do
    if (*name1=='\0')
      return (1);
  while ((*name2!='\0') && (tolower(*(name1++)) == tolower(*(name2++))));

  return (0);
}
/****************************************************************************/
/*D
   SearchUgCmd - Find a UG command by (part of) name

   SYNOPSIS:
   COMMAND *SearchUgCmd (const char *cmdName);

   PARAMETERS:
   .  cmdName - name of command

   DESCRIPTION:
   This function searches for a command with name 'name'. The difference
   to 'GetCommand' is that the name need not be complete. Only the first few
   characters must be supplied that make the name unique.

   RETURN VALUE:
   COMMAND *
   .n      pointer to 'COMMAND' if found and unique
   .n      NULL if not found or ambigous
   D*/
/********************************************************************************/
COMMAND *SearchUgCmd (const char *cmdName)
{
  ENVDIR *currentDir;
  COMMAND *Cmd;
  ENVITEM *theItem;
  char text[128];

  if (ChangeEnvDir("/Menu")==NULL)
  {
    UserWrite("ERROR: could not ChangeDir to /Menu\n");
    return(NULL);
  }

  currentDir = GetCurrentDir();

  /* loop in current directory */
  theItem = currentDir->down;
  Cmd = NULL;
  while (theItem!=NULL)
  {
    if (theItem->v.type == theCommandVarID)
      if (Str1inStr2(cmdName,theItem->v.name))
        if (Cmd!=NULL)
        {
          sprintf(text," ambiguos: %s\n",cmdName);
          UserWrite(text);
          return (NULL);
        }
        else
          Cmd = (COMMAND *) theItem;
    theItem = theItem->v.next;
  }

  return (Cmd);
}

/****************************************************************************/
/*D
   ReplaceCommand - Change an existing command or create a new one

   SYNOPSIS:
   COMMAND *ReplaceCommand (const char *name, CommandProcPtr cmdProc)

   PARAMETERS:
   .  name - name of a command
   .  cmdProc - new command function

   DESCRIPTION:
   This function changes the execution function of an existing command
   or creates a new one in case it does not exist already.

   RETURN VALUE:
   COMMAND *
   .n     pointer to new or replaced command
   .n     NULL in case of an error
   D*/
/****************************************************************************/

COMMAND *ReplaceCommand (const char *name, CommandProcPtr cmdProc)
{
  COMMAND *theCommand;

  /* change to menu directory */
  if (ChangeEnvDir("/Menu")==NULL)
    return (NULL);

  /* allocate new command only if command does not exist */
  if ((theCommand=GetCommand(name))==NULL)
    if ((theCommand = (COMMAND *) MakeEnvItem(name,theCommandVarID,sizeof(COMMAND)))==NULL)
      return(NULL);

  /* fill data */
  theCommand->cmdProc = cmdProc;

  return(theCommand);
}


/****************************************************************************/
/*D
   ExecCommand - Processes the command line  and execute command

   SYNOPSIS:
   INT ExecCommand (char *cmdLine)

   PARAMETERS:
   .  cmdLine - character string containing a complete command line

   DESCRIPTION:
   This function processes the command line, i.e. it constructs the 'argc'
   'argv' arrays and calls the execution function of the command.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT ExecCommand (char *cmdLine)
{
  char *s,*token,commandstr[NAMESIZE];
  int i,j,error;
  COMMAND *commandItem;
        #ifdef ModelP
  char cmd[OPTIONBUFFERLEN],*cmdptr;
        #endif

  optionCount = 0;
  s = optionBuffer;
        #ifdef ModelP
  strncpy(cmd,cmdLine,OPTIONBUFFERLEN);
  cmdptr = cmd;
        #endif
  token = strtok(cmdLine,OPTIONSEP);
  while (token!=NULL)
  {
    if (optionCount>=MAXOPTIONS)
    {
      PrintErrorMessage('E',"ExecCommand","too many options");
      return(8410);                     /* too many options */
    }

    strcpy(s,token);
    options[optionCount++] = s;
    s = s+(strlen(token)+1);
    token = strtok(NULL,OPTIONSEP);
  }

  if (optionCount==0)
    return (1);

  /* strip trailing blanks from arguments */
  for (i=0; i<optionCount; i++)
  {
    s = options[i];
    if (strlen(s)>0)
      for (j=strlen(s)-1; strchr(WHITESPACE,(int) s[j])!=NULL; j--) s[j] = '\0';
  }

  if (sscanf(options[0],expandfmt(CONCAT3("%",NAMELENSTR,"[a-zA-Z_0-9]")),commandstr)!=1)
    return (2);

  commandItem = GetCommand(commandstr);
  if (commandItem==NULL) return(1);
  else
  {
        #ifdef ModelP
    /* TODO: introduces a mask sign to hide option sign $, the $ has to be masked */
    /*		 context sensitve in the command interpreter						  */
    /* if string which is assigned in SetCommand has options need to send  whole string */
    /* in options[0], this excludes the only Option $r									*/
    if (strcmp(commandstr,"set") == 0)
    {
      if (optionCount>1 && strcmp(options[1],"r")!=0)
      {
        optionCount = 1;
        error=(*commandItem->cmdProc)(optionCount,&cmdptr);
        return(error);
      }
    }
        #endif
  }

  if (strcmp("reperr",options[0])!=0) REP_ERR_RESET;

  error=(*commandItem->cmdProc)(optionCount,options);
  if (error==PARAMERRORCODE)
    UserWrite("ERROR: invalid parameters\n");
  if (error!=OKCODE && error!=QUITCODE)
    UserWrite("ERROR in command execution\n");

  return(error);
}

INT InitCmdline ()
{
  /* install the /Menu directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitCmdline","could not changedir to root");
    return(__LINE__);
  }
  theMenuDirID = GetNewEnvDirID();
  if (MakeEnvItem("Menu",theMenuDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitCmdline","could not install '/Menu' dir");
    return(__LINE__);
  }
  theCommandVarID = GetNewEnvVarID();

  return (0);
}
