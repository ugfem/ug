// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  helpmsg.c                                                                                                     */
/*																			*/
/* Purpose:   implements a flexible online help system						*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							*/
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

#include "helpmsg.h"
#include "cmdline.h"
#include "num.h"
#include "devices.h"
#include "misc.h"
#include "fileopen.h"
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

/* for help functions */
#define BUFFERSIZE              256
#define LONGBUFFSIZE    1024

#define HELPSEP                 '-'     /* sperator of help pages					*/
#define KEWORDCHAR              '>'     /* seperator for keyword list				*/
#define FILENAMESEP     " \t"   /* token seperators for help filenames		*/
#define MAXHELPFILES    8               /* max number of helpfiles					*/
#define NOTFOUNDMESSAGE "sorry, no message file available\n"

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

/* general purpose text buffer */
static char longbuff[LONGBUFFSIZE];     /* this way Mac output is way faster*/
static char buffer[BUFFERSIZE];

/* for help functions */
static FILE *HelpFile[MAXHELPFILES];    /* the help messages files			*/
static int NHelpFiles;                                  /* no of open helpfiles                         */

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   PrintHelp - Print entry from help file

   SYNOPSIS:
   INT PrintHelp (const char *HelpFor,int mode, const char *addText);

   PARAMETERS:
   .  HelpFor - command name or keyword
   .  mode - operation mode (see below)
   .  addText - additional text to be printed

   DESCRIPTION:
   This command processes the helpfiles declared in the 'defaults'
   file of the application and prints the corresponding message
   for the command if 'mode' is 'HELPITEM' or all messages that
   have the corresponding keyword, if 'mode' is 'KEYWORD'.

   RETURN VALUE:
   INT
   .n 'HELP_OK'
   .n 'HELP_STRING_EMPTY'
   .n 'HELP_NOT_FOUND'
   .n 'HELP_STRING_TOO_LONG'

   D*/
/****************************************************************************/

INT PrintHelp (const char *HelpFor,int mode, const char *addText)
{
  char *s,HelpItem[64],helpfor[BUFFERSIZE];
  INT i,len,found;

  if (strlen(HelpFor)==0)
    return (HELP_STRING_EMPTY);
  if (strlen(HelpFor)>BUFFERSIZE-1)
    return (HELP_STRING_TOO_LONG);

  /* case INSENSITIVE: convert HelpFor to lower case */
  strcpy(helpfor,HelpFor);
  s = helpfor;
  while ((*s=tolower(*(s)))!='\0') s++;

  if (mode==KEYWORD)
  {
    /* print first lines of all help pages with helpfor in keyword list
       (not necessary entire word) */
    found = 0;
    for (i=0; i<NHelpFiles; i++)
    {
      if (HelpFile[i]==NULL) continue;

      rewind(HelpFile[i]);

      while (fgets(buffer,BUFFERSIZE-1,HelpFile[i])!=NULL)
        if (buffer[0]==HELPSEP)
        {
          /* old version: look only in keywords
             s = strchr(buffer,KEWORDCHAR);
             if ((s!=NULL)&&(strstr(s,helpfor)!=NULL)) */
          if (strstr(buffer,helpfor)!=NULL)
            /* matching: get and print next line */
            if ((fgets(buffer,BUFFERSIZE-1,HelpFile[i])!=NULL)&&(buffer[0]!=HELPSEP))
            {
              found++;
              sprintf(longbuff,"help %s",buffer);
              UserWrite(longbuff);
            }
        }
    }

    if (found)
      return (HELP_OK);
    else
      return (HELP_NOT_FOUND);
  }
  else
  {
    /* print help page with help item helpfor */
    longbuff[0] = '\0';
    len = 0;
    for (i=0; i<NHelpFiles; i++)
    {
      if (HelpFile[i]==NULL) continue;

      rewind(HelpFile[i]);

      while (fgets(buffer,BUFFERSIZE-1,HelpFile[i])!=NULL)
        if (buffer[0]==HELPSEP)
          if ((sscanf(buffer+1,"%s",HelpItem)==1)&&(strcmp(HelpItem,helpfor)==0))
          {
            while ((fgets(buffer,BUFFERSIZE-1,HelpFile[i])!=NULL)&&(buffer[0]!=HELPSEP))
              if ((len+=strlen(buffer))<LONGBUFFSIZE-1)
                strcat(longbuff,buffer);
              else
              {
                UserWrite(longbuff);
                longbuff[0] = '\0';
                len = 0;
              }

            UserWrite(longbuff);
            if (addText!=NULL) {UserWrite(addText); UserWrite("\n");}

            return (HELP_OK);
          }
    }
  }

  /* help not found */
  if (addText!=NULL)
  {
    UserWrite(addText);
    UserWrite("\n");
  }

  return (HELP_NOT_FOUND);
}

/****************************************************************************/
/*D
   CheckHelp - Check if all commands have a help message

   SYNOPSIS:
   INT CheckHelp (void);

   PARAMETERS:
   .  void - no parameters

   DESCRIPTION:
   This command checks whether all UG commands have a corresponding
   entry in any of the help files. All commands without entry
   are lsted.

   RETURN VALUE:
   INT
   .n 0 all commands have entry in help file
   .n 1 one or more helps missing
   D*/
/****************************************************************************/

INT CheckHelp ()
{
  COMMAND *theCmd;
  NP_TYPE *theNP;
  char HelpItem[128],cmdname[NAMESIZE],npname[NAMESIZE],*s;
  int i,found,rv;

  /* loop commands */
  UserWrite("checking commands...\n");
  rv = 0;
  for (theCmd=GetFirstCommand(); theCmd!=NULL; theCmd=GetNextCommand(theCmd))
  {
    found = FALSE;
    strcpy(cmdname,ENVITEM_NAME(theCmd));

    /* case INSENSITIVE: convert the command name to lower case */
    s = cmdname;
    while ((*s=tolower(*(s)))!='\0') s++;

    for (i=0; i<NHelpFiles; i++)
    {
      if (HelpFile[i]==NULL) continue;

      rewind(HelpFile[i]);

      while (fgets(buffer,BUFFERSIZE-1,HelpFile[i])!=NULL)
        if (buffer[0]==HELPSEP)
          if ((sscanf(buffer+1,"%s",HelpItem)==1)&&(strcmp(HelpItem,cmdname)==0))
          {
            found = TRUE;
            break;
          }

      if (found)
        break;
    }
    if (!found)
    {
      rv = 1;
      sprintf(buffer,"no help found for '%s'\n",ENVITEM_NAME(theCmd));
      UserWrite(buffer);
    }
  }
  if (rv)
    UserWrite("for all other commands on-line help is available\n\n");
  else
    UserWrite("for all commands on-line help is available\n\n");

  /* loop num proc types */
  UserWrite("checking num proc types...\n");
  rv = 0;
  for (theNP=GetFirstNumProcType(); theNP!=NULL; theNP=GetNextNumProcType(theNP))
  {
    found = FALSE;
    strcpy(npname,ENVITEM_NAME(theNP));

    /* case INSENSITIVE: convert the command name to lower case */
    s = npname;
    while ((*s=tolower(*(s)))!='\0') s++;

    for (i=0; i<NHelpFiles; i++)
    {
      if (HelpFile[i]==NULL) continue;

      rewind(HelpFile[i]);

      while (fgets(buffer,BUFFERSIZE-1,HelpFile[i])!=NULL)
        if (buffer[0]==HELPSEP)
          if ((sscanf(buffer+1,"%s",HelpItem)==1)&&(strcmp(HelpItem,npname)==0))
          {
            found = TRUE;
            break;
          }

      if (found)
        break;
    }
    if (!found)
    {
      rv = 1;
      sprintf(buffer,"no help found for '%s'\n",ENVITEM_NAME(theNP));
      UserWrite(buffer);
    }
  }
  if (rv)
    UserWrite("for all other num proc types on-line help is available\n");
  else
    UserWrite("for all num proc types on-line help is available\n");

  return (rv);
}

/****************************************************************************/
/*																			*/
/* Function:  InitHelpMsg													*/
/*																			*/
/* Purpose:   init this module (open help messages file)					*/
/*																			*/
/* Input:	  char *helpfile: name of the help messages file				*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*			  INT>0: could not open that file								*/
/*																			*/
/****************************************************************************/

INT InitHelpMsg (char *helpfiles)
{
  char *token,*FileName[MAXHELPFILES];
  int i;

  /* help files */
  NHelpFiles = 0;
  token = strtok(helpfiles,FILENAMESEP);
  while (token!=NULL)
  {
    if (NHelpFiles>=MAXHELPFILES)
      return (1);                       /* too many items */

    FileName[NHelpFiles++] = token;
    token = strtok(NULL,FILENAMESEP);
  }

  for (i=0; i<NHelpFiles; i++)
    HelpFile[i] = fileopen(FileName[i],"r");

  return(0);
}
