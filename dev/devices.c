// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  devices.c                                                                                                     */
/*																			*/
/* Purpose:   Initialization and hardware independent part of devices		*/
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

/* low module */
#include "compiler.h"
#include "fileopen.h"
#include "misc.h"
#include "ugenv.h"
#include "defaults.h"
#include "debug.h"
#include "ugstruct.h"
#include "general.h"

/* dev module */
#include "devices.h"
#include "initdev.h"

/* dddif module */
#ifdef ModelP
#include "ppif.h"
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

/* the mute level is set by the MuteCommand and used for output control		*/
/* convention: 0 is default, <0 produces less output, >0 produces more outpu*/
static INT mutelevel=0;

static FILE *logFile=NULL;                                              /* log file pointer             */
static OUTPUTDEVICE *defaultOuputDevice=NULL;   /* console graphical outp.	*/

static INT theOutputDevDirID;                   /* env type for Output Device dir	*/
static INT theOutputDevVarID;                   /* env type for Output Device vars	*/

static char ToolName[nboftools][NAMESIZE]; /* array of tool names			*/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   CreateOutputDevice - allocate OUTPUTDEVICE structure in environment

   SYNOPSIS:
   OUTPUTDEVICE *CreateOutputDevice (char *name);

   PARAMETERS:
   .  name - name of the device

   DESCRIPTION:
   This function allocates OUTPUTDEVICE structure in environment

   RETURN VALUE:
   OUTPUTDEVICE *
   .n    pointer to requested structure
   .n    NULL if operation failed.
   D*/
/****************************************************************************/

OUTPUTDEVICE *CreateOutputDevice (char *name)
{
  OUTPUTDEVICE *dev;

  if (ChangeEnvDir("/Output Devices")==NULL)
    return (NULL);

  if ((dev=(OUTPUTDEVICE *)MakeEnvItem(name,theOutputDevVarID,sizeof(OUTPUTDEVICE))) == NULL )
  {
    printf("error: cannot create output device %s\n",name);
    return(NULL);
  }

  return(dev);
}

/****************************************************************************/
/*D
   GetOutputDevice - search an OUTPUTDEVICE structure with name in the environment

   SYNOPSIS:
   OUTPUTDEVICE *GetOutputDevice (const char *name);

   PARAMETERS:
   .  name - name of the device

   DESCRIPTION:
   This function returns an OUTPUTDEVICE structure with name in the environment.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n    pointer to requested structure
   .n    NULL if operation failed.
   D*/
/****************************************************************************/

OUTPUTDEVICE *GetOutputDevice (const char *name)
{
  return((OUTPUTDEVICE *) SearchEnv(name,"/Output Devices",theOutputDevVarID,theOutputDevDirID));
}

/****************************************************************************/
/*D
   GetDefaultOutputDevice - return OUTPUTDEVICE structure for screen output

   SYNOPSIS:
   OUTPUTDEVICE *GetDefaultOutputDevice (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function returns the OUTPUTDEVICE structure for screen output.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n     pointer to
   .n     NULL if no screen output device available.
   D*/
/****************************************************************************/

OUTPUTDEVICE *GetDefaultOutputDevice (void)
{
  return(defaultOuputDevice);
}

/****************************************************************************/
/*D
   OpenLogFile - open a log file

   SYNOPSIS:
   INT OpenLogFile (const char *name);

   PARAMETERS:
   .  name -

   DESCRIPTION:
   This function opens a log file where all output to 'UserWrite', 'UserWriteF'
   and 'PrintErrorMessage' is protocoled.

   RETURN VALUE:
   INT
   .n    0 if operation ok
   .n    1 if a file is already open
   .n    2 if file open failed.
   D*/
/****************************************************************************/

INT OpenLogFile (const char *name)
{
  char logpath[BUFFSIZE];

  if (logFile!=NULL) return(1);

  /* get path to logfile directory */
  if (GetDefaultValue(DEFAULTSFILENAME,"logfilesdir",logpath)==0)
    logFile = FileOpenUsingSearchPath(name,"w",logpath);
  else
    logFile = fileopen(name,"w");

  if (logFile==NULL) return(2);

  return(0);
}

/****************************************************************************/
/*D
   CloseLogFile	- Close log file

   SYNOPSIS:
   INT CloseLogFile (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function closes the log file.

   RETURN VALUE:
   INT
   .n    0 if operation ok
   .n    1 if an error occured.
   D*/
/****************************************************************************/

INT CloseLogFile (void)
{
  if (logFile==NULL) return(1);

  fclose(logFile);
  logFile = NULL;
  return(0);
}

/****************************************************************************/
/*D
   SetLogFile - set log file ptr

   SYNOPSIS:
   INT SetLogFile (FILE *file);

   PARAMETERS:
   .  file -

   DESCRIPTION:
   This function et log file ptr.

   RETURN VALUE:
   INT
   .n    0 if operation ok
   .n    1 if an error occured.
   D*/
/****************************************************************************/

INT SetLogFile (FILE *file)
{
  logFile = file;
  return(0);
}

/****************************************************************************/
/*D
   WriteLogFile	- Write to (open) log file

   SYNOPSIS:
   INT WriteLogFile (const char *text);

   PARAMETERS:
   .  text -

   DESCRIPTION:
   This function writes to (open) log file.

   RETURN VALUE:
   INT
   .n    0 if operation ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT WriteLogFile (const char *text)
{
  if (logFile==NULL) return(1);

  fputs(text,logFile);
        #ifdef Debug
  fflush(logFile);
        #endif

  return(0);
}

/****************************************************************************/
/*D
   UserWrite - Write a string to shell window

   SYNOPSIS:
   void UserWrite (const char *s);

   PARAMETERS:
   .  s -

   DESCRIPTION:
   This function writes a string to shell window with log file mechanism.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UserWrite (const char *s)
{
        #ifdef ModelP
  if (me==master)
  {
        #endif

  if (mutelevel>-1000)
    WriteString(s);
  if (logFile!=NULL) {
    fputs(s,logFile);
                #ifdef Debug
    fflush(logFile);
                #endif
  }

        #ifdef ModelP
}
else
{
  PRINTDEBUG(ui,1,("%d: %s\n", me,s))
  IFDEBUG(ui,0)
  if (logFile!=NULL) {
    fputs(s,logFile);
    fflush(logFile);
  }
  ENDDEBUG
}
        #endif
}

/****************************************************************************/
/*D
   UserWriteF - write a formated string to shell window

   SYNOPSIS:
   int UserWriteF (const char *format, ...);

   PARAMETERS:
   .  format -
   .  ... - list of arguments for format string

   DESCRIPTION:
   This function writes a formated string to shell
   window with log file  mechanism.

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

int UserWriteF (const char *format, ...)
{
  char buffer[256];
  va_list args;

  /* initialize args */
  va_start(args,format);

  vsprintf(buffer,format,args);

        #ifdef ModelP
  if (me==master)
  {
        #endif

  if (mutelevel>-1000)
    WriteString(buffer);

        #ifdef ModelP
}
else
{
  PRINTDEBUG(ui,1,("%d: %s\n", me,buffer))
  IFDEBUG(ui,0)
  if (logFile!=NULL) {
    fputs(buffer,logFile);
    fflush(logFile);
  }
  ENDDEBUG
}
        #endif

  if (logFile!=NULL) fputs(buffer,logFile);

  /* garbage collection */
  va_end(args);

  return (0);
}

/****************************************************************************/
/*D
   PrintErrorMessage - Formatted error output (also to log file)

   SYNOPSIS:
   void PrintErrorMessage (char type, const char *procName, const char *text);

   PARAMETERS:
   .  type - 'W','E','F'
   .  procName - name  of procedure where error occured
   .  text -  additional explanation

   DESCRIPTION:
   This function formats error output (also to log file).

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void PrintErrorMessage (char type, const char *procName, const char *text)
{
  char buffer[256];
  char classText[32];

  switch (type)
  {
  case 'W' :
    strcpy(classText,"WARNING");
    break;
  case 'E' :
    strcpy(classText,"ERROR");
    break;

  case 'F' :
    strcpy(classText,"FATAL");
    break;

  default :
    strcpy(classText,"USERERROR");
    break;
  }
  sprintf(buffer,"%s in %.20s: %.200s\n",classText,procName,text);
  UserWrite(buffer);
}

/****************************************************************************/
/*D
   SetMuteLevel - set mute level for verbosing level

   SYNOPSIS:
   void SetMuteLevel (INT mute);

   PARAMETERS:
   .  mute -

   DESCRIPTION:
   This function sets the mute level for verbosing level.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void SetMuteLevel (INT mute)
{
  mutelevel = mute;
}

/****************************************************************************/
/*D
   GetMuteLevel - return mute level for verbosing level

   SYNOPSIS:
   void SetMuteLevel (INT mute);

   PARAMETERS:
   .  mute -

   DESCRIPTION:
   This function return the mute level for verbosing level.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

INT GetMuteLevel (void)
{
  return (mutelevel);
}

/****************************************************************************/
/*D
   SetToolName - set a tools name

   SYNOPSIS:
   INT SetToolName (INT tool, const char *name)

   PARAMETERS:
   .  tool - tool id
   .  name - name string

   DESCRIPTION:
   This function sets a tools name.

   RETURN VALUE:
   INT
   .n   0 ok
   .n           1 wrong tool id
   .n   2 name too long
   D*/
/****************************************************************************/

INT SetToolName (INT tool, const char *name)
{
  if (tool>nboftools)
    return (1);
  if (strlen(name)>NAMELEN)
    return (2);

  strcpy(ToolName[tool],name);

  return (0);
}

/****************************************************************************/
/*D
   GetToolName - get a tools name

   SYNOPSIS:
   const char *GetToolName (INT tool)

   PARAMETERS:
   .  tool - tool id

   DESCRIPTION:
   This function returns a tools name.

   RETURN VALUE:
   INT
   .n   NULL error
   .n      else string pointer
   D*/
/****************************************************************************/

const char *GetToolName (INT tool)
{
  if (tool>nboftools)
    return (NULL);

  return (ToolName[tool]);
}

/****************************************************************************/
/*D
   InitDevices - Initialize all devices at startup

   SYNOPSIS:
   INT InitDevices (int argc, char **argv);

   PARAMETERS:
   .  argc - argument counter
   .  argv - command line parameters

   DESCRIPTION:
   This function initializes all devices at startup.
   It must be extended when an output device is added.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if some error occured.
   D*/
/****************************************************************************/


INT InitDevices (int argc, char **argv)
{
  ENVDIR *DevDir;
  ENVITEM *dev;
  INT error,i;
  char sv[32];

  /* init tool names */
  strcpy(ToolName[arrowTool],arrowToolName);
  strcpy(ToolName[crossTool],crossToolName);
  strcpy(ToolName[choiceTool],choiceToolName);
  strcpy(ToolName[circleTool],circleToolName);
  strcpy(ToolName[handTool],handToolName);
  strcpy(ToolName[heartTool],heartToolName);
  strcpy(ToolName[gnoedelTool],gnoedelToolName);

  /* install the /Output Devices directory */
  if (ChangeEnvDir("/")==NULL)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }
  theOutputDevDirID = GetNewEnvDirID();
  if ((DevDir=(ENVDIR*)MakeEnvItem("Output Devices",theOutputDevDirID,sizeof(ENVDIR)))==NULL)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }
  theOutputDevVarID = GetNewEnvVarID();

  /* init screen device */
        #ifdef ModelP
  if (me == master)
  {
        #endif
  defaultOuputDevice = InitScreen(argc,argv,&error);
  if (error) return(1);
        #ifdef ModelP
}
        #endif

  /* init metafile device */
  if (InitMeta()!=0)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }

  /* init postscript device */
  if (InitPostScript()!=0)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }

  /* create struct and fill stringvars */
  if (MakeStruct(":Devices")!=0)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }
  for (i=0, dev=ENVDIR_DOWN(DevDir); dev!=NULL; i++, dev=NEXT_ENVITEM(dev))
  {
    sprintf(sv,":Devices:device%d",(int)i);
    if (SetStringVar(sv,ENVITEM_NAME(dev))!=0)
    {
      SetHiWrd(error,__LINE__);
      return (error);
    }
  }
  if (SetStringValue(":Devices:nDevices",i)!=0)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }

  return(0);
}
