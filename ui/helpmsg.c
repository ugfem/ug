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
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   17.12.94 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#define KEEP_HELPFILES_OPEN

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
#include "defaults.h"
#include "cmdline.h"
#include "ugdevices.h"
#include "misc.h"
#include "fileopen.h"
#include "general.h"
#include "debug.h"

#ifdef ModelP
#include "parallel.h"
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

/* for help functions */
#define BUFFERSIZE                      256
#define LONGBUFFSIZE            1024

#define HELPFILE_LIST           "lib/ugdata/helpfile.list"

#define HELPSEP                         '-'     /* sperator of help pages				*/
#define KEWORDCHAR                      '>'     /* seperator for keyword list			*/
#define FILENAMESEP             " \t\n" /* token seperators for help filenames	*/
#define MAXHELPFILES            50              /* max number of helpfiles				*/

#define REPLACE_DOT                     " "
#define REPLACE_DOT_N           "  "
#define REPLACE_TILDE           ' '
#define INDENT_VERBATIM         ":   "
#define TAB_SIZE                        4

#define DOC_TEXT_BEGIN(s)       ((s[0]=='/') && (s[1]=='*') && (s[2]=='D'))
#define DOC_TEXT_END(s)         ((s[0]=='D') && (s[1]=='*') && (s[2]=='/'))
#define BEGIN_VERBATIM(s)       ((s[0]=='.') && (s[1]=='v') && (s[2]=='b'))
#define END_VERBATIM(s)         ((s[0]=='.') && (s[1]=='v') && (s[2]=='e'))

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
static char buffer2[BUFFERSIZE];
static char LowerBuffer[BUFFERSIZE];

/* for help functions */
static char *HelpFile[MAXHELPFILES];    /* the help messages files			*/
static FILE *HelpFilePtr[MAXHELPFILES]; /* the help messages files			*/
static int NHelpFiles;                                  /* no of open helpfiles                         */

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static char *ToLower (char *s)
{
  char *p;

  p = LowerBuffer;
  while (*s!='\0')
    *p++ = tolower(*s++);

  *p = '\0';

  return (LowerBuffer);
}

static void WriteFormatted (const char *text)
{
  char buffer[LONGBUFFSIZE];
  INT tp,bp;
  static INT verbatim = FALSE;

  tp = bp = 0;
  buffer[0] = '\0';

  /* consider begin of line */
  if (verbatim)
  {
    if (END_VERBATIM(text))
    {
      verbatim = FALSE;

      /* end verbatim: skip this line */
      return;
    }
    strcat(buffer,INDENT_VERBATIM);
    bp = strlen(INDENT_VERBATIM);
  }
  else if (text[0]=='.')
  {
    if (text[1]=='n')
    {
      strcat(buffer,REPLACE_DOT_N);
      tp = 2;
      bp = strlen(REPLACE_DOT_N);
    }
    else if (BEGIN_VERBATIM(text))
    {
      verbatim = TRUE;

      /* begin verbatim: skip this line */
      return;
    }
    else
    {
      strcat(buffer,REPLACE_DOT);
      tp = 1;
      bp = strlen(REPLACE_DOT);
    }
  }

  /* copy text, replace tilde and replace tabs by runs of spaces */
  while (text[tp]!='\0')
  {
    switch (text[tp])
    {
    case '\t' :
      do
        buffer[bp++] = ' ';
      while (bp%TAB_SIZE);
      break;
    case '~' :
      buffer[bp++] = REPLACE_TILDE;
      break;
    default :
      buffer[bp++] = text[tp];
    }
    tp++;
  }

  /* terminate string */
  buffer[bp] = '\0';

  UserWrite(buffer);
}

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

        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

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
      FILE *hfile =
#                               ifdef KEEP_HELPFILES_OPEN
        HelpFilePtr[i];
#                               else
        fileopen(HelpFile[i],"r");
#                               endif

      if (hfile==NULL) continue;
      rewind(hfile);

      while (fgets(buffer,BUFFERSIZE-1,hfile)!=NULL)
        if (DOC_TEXT_BEGIN(buffer))
        {
          /* move on to first line */
          do
            if (fgets(buffer,BUFFERSIZE-1,hfile)==NULL)
            {
#                                                       ifndef KEEP_HELPFILES_OPEN
              fclose(hfile);
#                                                       endif
              REP_ERR_RETURN(1);
            }
          while (sscanf(buffer,"%s",HelpItem)!=1);

          PRINTDEBUG(ui,1,("PrintHelp: first line is '%s'\n",buffer));

          if ((sscanf(ToLower(buffer),"%s",HelpItem)==1)&&(strstr(HelpItem,helpfor)!=NULL))
          {
            /* matching: print line */
            found++;
            WriteFormatted(buffer);
          }
          else
          {
            /* search KEYWORDS line */
            while ((fgets(buffer2,BUFFERSIZE-1,hfile)!=NULL) && !DOC_TEXT_END(buffer2))
              if (strstr(buffer2,"KEYWORDS")!=NULL)
              {
                PRINTDEBUG(ui,1,("PrintHelp: keywords line is '%s'\n",buffer2));

                /* move on to next line and check it */
                if (fgets(buffer2,BUFFERSIZE-1,hfile)==NULL)
                {
#                                                                       ifndef KEEP_HELPFILES_OPEN
                  fclose(hfile);
#                                                                       endif
                  REP_ERR_RETURN(1);
                }
                if (strstr(ToLower(buffer2),helpfor)!=NULL)
                {
                  /* got it: print first line */
                  found++;
                  WriteFormatted(buffer);
                }
                break;
              }
          }
          /* move on to end of doc text */
          do
            if (DOC_TEXT_END(buffer))
              break;
          while (fgets(buffer,BUFFERSIZE-1,hfile)!=NULL);
        }
#                       ifndef KEEP_HELPFILES_OPEN
      fclose(hfile);
#                       endif
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
      FILE *hfile =
#                               ifdef KEEP_HELPFILES_OPEN
        HelpFilePtr[i];
#                               else
        fileopen(HelpFile[i],"r");
#                               endif

      if (hfile==NULL) continue;
      rewind(hfile);

      while (fgets(buffer,BUFFERSIZE-1,hfile)!=NULL)
        if (DOC_TEXT_BEGIN(buffer))
        {
          /* move on to first line */
          do
            if (fgets(buffer,BUFFERSIZE-1,hfile)==NULL)
            {
#                                                       ifndef KEEP_HELPFILES_OPEN
              fclose(hfile);
#                                                       endif
              REP_ERR_RETURN(1);
            }
          while (sscanf(buffer,"%s",HelpItem)!=1);

          /* scan first word */
          if ((sscanf(ToLower(buffer),"%s",HelpItem)==1)&&(strcmp(HelpItem,helpfor)==0))
          {
            /* print help text including first line */
            do
              WriteFormatted(buffer);
            while ((fgets(buffer,BUFFERSIZE-1,hfile)!=NULL) && !DOC_TEXT_END(buffer));

            if (addText!=NULL)
              UserWriteF("%s\n",addText);

#                                               ifndef KEEP_HELPFILES_OPEN
            fclose(hfile);
#                                               endif
            return (HELP_OK);
          }
        }
#                       ifndef KEEP_HELPFILES_OPEN
      fclose(hfile);
#                       endif
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
   are listed.

   RETURN VALUE:
   INT
   .n 0 all commands have entry in help file
   .n 1 one or more helps missing
   D*/
/****************************************************************************/

INT CheckHelp ()
{
  COMMAND *theCmd;
  char HelpItem[128],cmdname[NAMESIZE],*s;
  int i,found,rv;

        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

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
      FILE *hfile =
#                               ifdef KEEP_HELPFILES_OPEN
        HelpFilePtr[i];
#                               else
        fileopen(HelpFile[i],"r");
#                               endif

      if (hfile==NULL) continue;
      rewind(hfile);

      while (fgets(buffer,BUFFERSIZE-1,hfile)!=NULL)
        if (DOC_TEXT_BEGIN(buffer))
        {
          /* move on to next line */
          if (fgets(buffer,BUFFERSIZE-1,hfile)==NULL)
          {
#                                               ifndef KEEP_HELPFILES_OPEN
            fclose(hfile);
#                                               endif
            REP_ERR_RETURN(1);
          }

          /* scan first word */
          if ((sscanf(ToLower(buffer),"%s",HelpItem)==1)&&(strcmp(HelpItem,cmdname)==0))
          {
            found = TRUE;
            break;
          }
        }
#                       ifndef KEEP_HELPFILES_OPEN
      fclose(hfile);
#                       endif
      if (found)
        break;
    }
    if (!found)
    {
      if (rv==0)
        UserWrite("no help found for:\n");

      rv = 1;
      UserWriteF("    '%s'\n",ENVITEM_NAME(theCmd));
    }
  }
  if (rv)
    UserWrite("for all other commands on-line help is available\n\n");
  else
    UserWrite("for all commands on-line help is available\n\n");


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

INT InitHelpMsg (void)
{
  FILE *fp;
  char *token,buffer[BUFFSIZE+64],path2ug[64],fname[64],*s;

        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

  NHelpFiles = 0;

  /* get application help files */
  if (GetDefaultValue(DEFAULTSFILENAME,"helpfiles",buffer))
  {
    PrintErrorMessageF('W',"InitHelpMsg","could not read 'helpfiles' in defaults file '%s'",DEFAULTSFILENAME);
  }
  else
  {
    /* open help files */
    token = strtok(buffer,FILENAMESEP);
    while (token!=NULL)
    {
      FILE *hfile;

      if (NHelpFiles>=MAXHELPFILES)
      {
        PrintErrorMessage('W',"InitHelpMsg","max number of helpfiles exceeded");
        return (__LINE__);
      }

      if ((hfile=fileopen(token,"r"))==NULL)
      {
        PrintErrorMessageF('W',"InitHelpMsg","could not open file '%s'",token);
      }
      else
      {
#                               ifdef KEEP_HELPFILES_OPEN
        HelpFilePtr[NHelpFiles] = hfile;
#                               else
        fclose(hfile);
#                               endif
        HelpFile[NHelpFiles] = StrDup(token);
        if (HelpFile[NHelpFiles]!=NULL)
          NHelpFiles++;
      }

      token = strtok(NULL,FILENAMESEP);
    }
  }
  /* get path from application to UG/ug */
  if (GetDefaultValue(DEFAULTSFILENAME,"path2ug",buffer)!=0)
  {
    PrintErrorMessageF('W',"InitHelpMsg","could not read 'path2ug' in defaults file '%s'",DEFAULTSFILENAME);
    return (__LINE__);
  }
  if (ExpandCShellVars(buffer)==NULL)
  {
    PrintErrorMessageF('W',"InitHelpMsg","could not expand shell variables in 'path2ug' of defaults file '%s'",DEFAULTSFILENAME);
    return (__LINE__);
  }
  if (sscanf(buffer,"%s",path2ug)!=1)
  {
    PrintErrorMessageF('W',"InitHelpMsg","no path specified in 'path2ug' in defaults file '%s'",DEFAULTSFILENAME);
    return (__LINE__);
  }

  /* now read ug's helpfiles and add path2ug */
  strcpy(buffer,path2ug);
  strcat(buffer,HELPFILE_LIST);
  fp = fileopen(buffer,"r");
  if (fp==NULL)
  {
    PrintErrorMessageF('W',"InitHelpMsg","could not open ug's helpfile list '%s'",buffer);
    return (__LINE__);
  }

  s = buffer + strlen(path2ug);
  while (fgets(s,BUFFERSIZE-1,fp)!=NULL)
  {
    FILE *hfile;

    if (NHelpFiles>=MAXHELPFILES)
    {
      PrintErrorMessage('W',"InitHelpMsg","max number of helpfiles exceeded");
      return (__LINE__);
    }

    if (sscanf(buffer,"%s",fname)!=1)
    {
      PrintErrorMessageF('W',"InitHelpMsg","could not scan filename from '%s'",buffer);
      continue;
    }

    if ((hfile=fileopen(fname,"r"))==NULL)
    {
      PrintErrorMessageF('W',"InitHelpMsg","could not open file '%s'",fname);
      continue;
    }
#               ifdef KEEP_HELPFILES_OPEN
    HelpFilePtr[NHelpFiles] = hfile;
#               else
    fclose(hfile);
#               endif
    HelpFile[NHelpFiles] = StrDup(fname);
    if (HelpFile[NHelpFiles]!=NULL)
      NHelpFiles++;
  }
  fclose (fp);

  return(0);
}
