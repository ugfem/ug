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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

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
   UgSetPalette - set color/black-white/gray palette

   SYNOPSIS:
   void UgSetPalette (OUTPUTDEVICE *dev, INT palette);

   PARAMETERS:
   .  dev - outputdevice
   .  palette - color/bw/gray

   DESCRIPTION:
   This function sets a color/black-white/gray palette.

   RETURN VALUE:
   INT
   D*/
/****************************************************************************/

INT UgSetPalette (OUTPUTDEVICE *dev, INT palette)
{
  short red[256],green[256],blue[256];
  short i;

        #ifdef ModelP
  if (me != master)
    return (1);
        #endif

  if (dev==NULL)
    return (1);

  switch (palette)
  {
  case COLOR_PALETTE :
  {
    short res = 63;
    short delta = 4;
    short max = 252;
    short r,g,b,j;

    /* fixed colors */
    i = 0;
    red[i] = 255; green[i] = 255; blue[i++] = 255;                      /* 0 = white */
    red[i] = 255; green[i] = 0      ; blue[i++] = 255;                          /* 1 = magenta */

    /* color spectrum */
    r = g = 0; b = max;
    red[i] = r; green[i] = g; blue[i++] = b;                                    /* 2 = blue */

    /* blue to cyan */
    for (j=0; j<res; j++)
    {
      g += delta;
      red[i] = r; green[i] = g; blue[i++] = b;
    }                                                                           /* 65 = cyan */
    /* cyan to green */
    for (j=0; j<res; j++)
    {
      b -= delta;
      red[i] = r; green[i] = g; blue[i++] = b;
    }                                                                           /* 128 = green */
    /* green to yellow */
    for (j=0; j<res; j++)
    {
      r += delta;
      red[i] = r; green[i] = g; blue[i++] = b;
    }                                                                           /* 191 = yellow */
    /* yellow to rot */
    for (j=0; j<res; j++)
    {
      g -= delta;
      red[i] = r; green[i] = g; blue[i++] = b;
    }                                                                           /* 254 = red */
    red[i] = 0; green[i] = 0  ; blue[i++] = 0;                                  /* 255 = black */
    break;
  }
  case BLACK_WHITE_PALETTE :
    red[0] = green[0] = blue[0] = 0;                                    /* white */
    for (i=1; i<256; i++)
      red[i] = green[i] = blue[i] = 1;                                  /* black */
    break;
  case GRAY_PALETTE :
    for (i=0; i<256; i++)
      red[i] = green[i] = blue[i] = i;                                  /* gray  */
    break;
  default :
    return(1);
  }

  (*dev->SetNewPalette)(0,256,red,green,blue);

  return (0);
}

/****************************************************************************/
/*D
   OpenLogFile - open a log file

   SYNOPSIS:
   INT OpenLogFile (const char *name, int rename);

   PARAMETERS:
   .  name -
   .  rename -

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

INT OpenLogFile (const char *name, int rename)
{
  char logpath[BUFFSIZE];

  if (logFile!=NULL) return(1);

  /* get path to logfile directory */
  if (GetDefaultValue(DEFAULTSFILENAME,"logfilesdir",logpath)==0)
    logFile = FileOpenUsingSearchPath_r(name,"w",logpath,rename);
  else
    logFile = fileopen_r(name,"w",rename);

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

  if ( fputs(text,logFile) < 0 )
  {
    UserWrite( "ERROR in writing logfile\n" );
                #ifdef Debug
    printf( "ERROR in writing logfile\n" );
    fflush(logFile);
                #endif
    return 1;
  }
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
    if ( fputs(s,logFile) < 0 )
    {
      UserWrite( "ERROR in writing logfile\n" );
                        #ifdef Debug
      printf( "ERROR in writing logfile\n" );
      fflush(logFile);
                        #endif
    }
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
    if ( fputs(s,logFile) < 0 )
    {
      UserWrite( "ERROR in writing logfile\n" );
      printf( "ERROR in writing logfile\n" );
      fflush(logFile);
    }
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
    if ( fputs(buffer,logFile) < 0 )
    {
      UserWrite( "ERROR in writing logfile\n" );
      printf( "ERROR in writing logfile\n" );
    }
    fflush(logFile);
  }
  ENDDEBUG
}
        #endif

  if (logFile!=NULL) {
    if ( fputs(buffer,logFile) < 0 )
    {
      UserWrite( "ERROR in writing logfile\n" );
                        #ifdef Debug
      printf( "ERROR in writing logfile\n" );
      fflush(logFile);
                        #endif
      va_end(args);
      return 1;
    }
                #ifdef Debug
    fflush(logFile);
                #endif
  }

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

   SEE ALSO:
   'PrintErrorMessageF'
   D*/
/****************************************************************************/

void PrintErrorMessage (char type, const char *procName, const char *text)
{
  char classText[32];
  INT oldmutelevel;

  oldmutelevel = mutelevel;
  switch (type)
  {
  case 'W' :
    strcpy(classText,"WARNING");
    break;

  case 'E' :
    strcpy(classText,"ERROR");
    mutelevel = 0;
    break;

  case 'F' :
    strcpy(classText,"FATAL");
    mutelevel = 0;
    break;

  default :
    strcpy(classText,"USERERROR");
    break;
  }
  UserWriteF("%s in %.20s: %.200s\n",classText,procName,text);
  mutelevel = oldmutelevel;
}

/****************************************************************************/
/*D
   PrintErrorMessageF - Formatted error output with fotmatted message (also to log file)

   SYNOPSIS:
   void PrintErrorMessageF (char type, const char *procName, const char *format, ...);

   PARAMETERS:
   .  type - 'W','E','F'
   .  procName - name  of procedure where error occured
   .  format -  additional formatted explanation (like printf)

   DESCRIPTION:
   This function formats error output (also to log file).
   After expanding message 'PrintErrorMessage' is called.

   RETURN VALUE:
   void

   SEE ALSO:
   'PrintErrorMessage'
   D*/
/****************************************************************************/

void PrintErrorMessageF (char type, const char *procName, const char *format, ...)
{
  char buffer[256];
  va_list args;

  /* initialize args */
  va_start(args,format);

  vsprintf(buffer,format,args);

  PrintErrorMessage(type,procName,buffer);

  /* garbage collection */
  va_end(args);
}

/****************************************************************************/
/*D
   SetMuteLevel - set mute level for verbosing level

   SYNOPSIS:
   void SetMuteLevel (INT mute);

   PARAMETERS:
   .  mute - indicator of amount of output

   DESCRIPTION:
   This function sets the mute level for verbosing level.

   CONVENTION:
   'mute <= -1' cancels the echoing of `ug`-commands, 'mute >= 0'
   restores the echoing (which is the default state). 'mute <= -1000' suppresses
   also the output of `ug` commands.

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
   INT GetMuteLevel (void);

   PARAMETERS:
   .  void

   CONVENTION:
   'mute <= -1' cancels the echoing of `ug`-commands, 'mute >= 0'
   restores the echoing (which is the default state). 'mute <= -1000' suppresses
   also the output of `ug` commands.

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
   InitDevices - Initialize all devices at startup

   SYNOPSIS:
   INT InitDevices (int *argcp, char **argv);

   PARAMETERS:
   .  argcp - pointer to argument counter
   .  argv  - command line parameters

   DESCRIPTION:
   This function initializes all devices at startup.
   It must be extended when an output device is added.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if some error occured.
   D*/
/****************************************************************************/


INT InitDevices (int *argcp, char **argv)
{
  ENVDIR *DevDir;
  ENVITEM *dev;
  INT error=0,i,screen;
  char sv[32];
#       ifdef ModelP
  INT with_defaultOuputDevice;
#       endif
  char buffer[256];
  int ival;

  /* get default mutelevel from defaults file */
  if (GetDefaultValue(DEFAULTSFILENAME,"mutelevel",buffer)==0) {
    sscanf(buffer," %d ",&ival);
    SetMuteLevel ((INT) ival);
  }

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
  if (me == master) {
        #endif

  defaultOuputDevice = InitScreen(argcp,argv,&error);
  if (error) return(1);

        #ifdef ModelP
  /* send number of command line arguments after InitScreen() */
  Broadcast(argcp, sizeof(int));

  with_defaultOuputDevice = (defaultOuputDevice!=NULL);
}
else {
  int i, new_argc;

  /* get number of command line arguments after InitScreen() */
  Broadcast(&new_argc, sizeof(int));

  /* if number has been reduced, remove first arg from command line */
  while (new_argc<*argcp)
  {
    for(i=1; i < (*argcp)-1; i++)
    {
      argv[i] = argv[i+1];
    }
    if (*argcp > 1) (*argcp)--;
  }
}

Broadcast(&with_defaultOuputDevice, sizeof(int));
if (with_defaultOuputDevice)
{
  if (me!=master)
  {
    defaultOuputDevice = malloc(sizeof(OUTPUTDEVICE));
    /* TODO:  set function pointers to NULL */
  }

  Broadcast((void *)defaultOuputDevice,sizeof(OUTPUTDEVICE));
}
else
{
  if (me!=master)
    defaultOuputDevice = NULL;
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
  screen=0;
  for (i=0, dev=ENVDIR_DOWN(DevDir); dev!=NULL; i++, dev=NEXT_ENVITEM(dev))
  {
    sprintf(sv,":Devices:device%d",(int)i);
    if (SetStringVar(sv,ENVITEM_NAME(dev))!=0)
    {
      SetHiWrd(error,__LINE__);
      return (error);
    }
    if (strcmp(ENVITEM_NAME(dev),"screen")==0) screen=1;
  }
  if (SetStringValue(":Devices:nDevices",i)!=0)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }
  if (SetStringValue(":Devices:Screen",screen)!=0)
  {
    SetHiWrd(error,__LINE__);
    return (error);
  }

  return(0);
}



INT ExitDevices (void)
{
  /* clean up screen device */
  ExitScreen();

  return(0);         /* no error */
}
