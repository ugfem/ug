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

#include <config.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctype.h>

#include "ugtypes.h"
#include "debug.h"
#include "misc.h"
#include "ugenv.h"
#include "cmdline.h"
#include "ugdevices.h"
#include "general.h"

USING_UG_NAMESPACES

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
#define COMMENT_CHAR    '#'             /* ignore rest of line after #				*/

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
/** \brief Register the commands exported by this module in the environ

   \param name - Name of the command
   \param cmdProc - Pointer to a function of type 'CommandProcPtr'

   This function registers a new command that can be executed from the UG
   shell. This process is described in detail on the page 'commands'.

   \sa
   'commands'.

   \return <ul>
   <li> pointer to new 'COMMAND' structure if o.k </li>
   <li> NULL pointer in case of an error </li>
   </ul>
 */
/****************************************************************************/

COMMAND * NS_DIM_PREFIX CreateCommand (const char *name, CommandProcPtr cmdProc)
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
/** \brief Return pointer to command structure with name name

   \param name - Name of the command.

   This function returns pointer to command structure with name 'name'. If
   it does not exist a 'NULL' pointer is returned.

   \return <ul>
   <li> POINTER to 'COMMAND' structure </li>
   <li> NULL if not found </li>
   </ul>
 */
/****************************************************************************/

COMMAND * NS_DIM_PREFIX GetCommand (const char *name)
{
  COMMAND *cmd;

  if (ChangeEnvDir("/Menu")==NULL) return(NULL);
  cmd = (COMMAND *) SearchEnv(name,".",theCommandVarID,theMenuDirID);
  return(cmd);
}

/****************************************************************************/
/** \brief Return pointer to first command structure of /Menu dir

   This function returns pointer to first command structure of the '/Menu'
   environment directory.

   \return <ul>
   <li> POINTER to 'COMMAND' structure </li>
   <li> NULL if not found </li>
   </ul>
 */
/****************************************************************************/

COMMAND * NS_DIM_PREFIX GetFirstCommand ()
{
  ENVITEM *cmd;

  if ((cmd=(ENVITEM*)ChangeEnvDir("/Menu")) == NULL) return (NULL);

  for (cmd=ENVITEM_DOWN(cmd); cmd!=NULL; cmd=NEXT_ENVITEM(cmd))
    if (ENVITEM_TYPE(cmd) == theCommandVarID)
      return ((COMMAND*)cmd);
  return (NULL);
}

/****************************************************************************/
/** \brief Return pointer to command structure of /Menu dir following cmd

   \param cmd - pointer to a 'COMMAND'

   This function returns a pointer to the next command after 'cmd' in the
   list of commands.

   \return <ul>
   <li> POINTER to next command </li>
   <li> NULL if not found </li>
   </ul>
 */
/****************************************************************************/

COMMAND * NS_DIM_PREFIX GetNextCommand (const COMMAND *cmd)
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
   \param name1 -
   \param name2 -

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
/** \brief Find a UG command by (part of) name

   \param cmdName - name of command

   This function searches for a command with name 'name'. The difference
   to 'GetCommand' is that the name need not be complete. Only the first few
   characters must be supplied that make the name unique.

   \return <ul>
   <li> pointer to 'COMMAND' if found and unique </li>
   <li> NULL if not found or ambiguous </li>
   </ul>
 */
/********************************************************************************/

COMMAND * NS_DIM_PREFIX SearchUgCmd (const char *cmdName)
{
  ENVDIR *currentDir;
  COMMAND *Cmd;
  ENVITEM *theItem;

  if (ChangeEnvDir("/Menu")==NULL)
  {
    UserWrite("ERROR: could not ChangeDir to /Menu\n");
    return(NULL);
  }

  currentDir = GetCurrentDir();

  /* loop in current directory */
  Cmd = NULL;
  for (theItem=currentDir->down; theItem!=NULL; theItem=theItem->v.next)
    if (theItem->v.type == theCommandVarID)
      if (strcmp(cmdName,ENVITEM_NAME(theItem))==0)
        return ((COMMAND *) theItem);
      else if (Str1inStr2(cmdName,ENVITEM_NAME(theItem)))
      {
        if (Cmd!=NULL)
        {
          UserWriteF(" '%s' ambiguos:\n",cmdName);
          UserWriteF("      %s\n",ENVITEM_NAME(Cmd));
          UserWriteF("      %s\n",ENVITEM_NAME(theItem));
          while ((theItem = theItem->v.next)!=NULL)
            if (Str1inStr2(cmdName,ENVITEM_NAME(theItem)))
              UserWriteF("      %s\n",ENVITEM_NAME(theItem));
          return (NULL);
        }
        else
          Cmd = (COMMAND *) theItem;
      }

  return (Cmd);
}

/****************************************************************************/
/** \brief Change an existing command or create a new one

   \param name - name of a command
   \param cmdProc - new command function

   This function changes the execution function of an existing command
   or creates a new one in case it does not exist already.

   \return </ul>
   <li> pointer to new or replaced command </li>
   <li> NULL in case of an error </li>
   </ul>
 */
/****************************************************************************/

COMMAND * NS_DIM_PREFIX ReplaceCommand (const char *name, CommandProcPtr cmdProc)
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
/** \brief Processes the command line  and execute command

   \param cmdLine - character string containing a complete command line

   This function processes the command line, i.e. it constructs the 'argc'
   'argv' arrays and calls the execution function of the command.

   \return <ul>
   <li> 0 if ok </li>
   <li> 1 if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX ExecCommand (char *cmdLine)
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

  /* strip comments from arguments */
  for (i=0; i<optionCount; i++)
    if ((s=strchr(options[i],COMMENT_CHAR))!=NULL)
      *s = '\0';

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
#ifdef Debug
  IFDEBUG(ui,0)
  if ((error==OKCODE) && REP_ERR_ENCOUNTERED && strcmp("reperr",options[0])!=0)
  {
    PrintErrorMessageF('E',"ExecCommand","Huh??? %s returns OKCODE but rep err encountered",ENVITEM_NAME(commandItem));
    PrintRepErrStack(printf);
    fflush(stdout);
    return (FATAL);
  }
  ENDDEBUG
#endif

  return(error);
}

INT NS_DIM_PREFIX InitCmdline ()
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
