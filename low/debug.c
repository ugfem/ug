// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  debug.c                                                                                                               */
/*																			*/
/* Purpose:   ug internal debugger functions								*/
/*																			*/
/* Author:	  Stefan Lang                                                                                   */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: stefan@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   10.07.95 begin                                                                            */
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

#include "general.h"

#ifdef Debug
#define compile_debug

#include <stdio.h>
#include <stdarg.h>

#include "devices.h"
#include "fileopen.h"
#include "debug.h"

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

int Debuginit           =       0;
int Debugdddif              =       0;      /* temporary setting for debugging ModelP */
int Debugdev                =       0;
int Debugdom                =       0;
int Debuggm                 =       0;
int Debuggraph              =       0;
int Debuggui        =   0;
int Debuglow                =       0;
int Debugmachines   =       0;
int Debugnp         =       0;
int Debugui                 =       0;
int Debugappl               =       0;
int Debugpclib              =       0;
int Debugtime               =       0;

int rep_err_count;
int rep_err_line[REP_ERR_MAX];
const char  *rep_err_file[REP_ERR_MAX];

/* for timings */
int __debug_time_count;
double __debug_time[DEBUG_TIME_MAX];
int __debug_time_line[DEBUG_TIME_MAX];
const char  *__debug_time_file[DEBUG_TIME_MAX];

/* from dddif/ppif.h */
extern int me, master;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static PrintfProcPtr printdebug=printf;
static FILE                                     *debugfile=NULL;
static char                             *debugfilename;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* Function:  PrintDebug													*/
/*																			*/
/* Purpose:   Print debugging information to ugshell or logfile				*/
/*																			*/
/* Input:     arguments are passed to PrintDebug in a printf like manner.   */
/*			  char *format: format, which contains the debugging info	        */
/*              ...     :	list of arguments for format string				*/
/*																			*/
/* Output:    void															*/
/*																			*/
/****************************************************************************/

int PrintDebug (const char *format, ...)
{
  char buffer[4096];
  va_list args;

  va_start(args,format);

  vsprintf(buffer,format,args);
        #ifdef ModelP
  if (me==master) {
        #endif

  /* use specific debug function for displaying */
  (*printdebug)(buffer);
  WriteLogFile((const char *)buffer);

        #ifdef ModelP
}
else
{
  printf(buffer);
  fflush(stdout);
  WriteLogFile(buffer);
}

        #endif

  va_end(args);
  return (0);
}

void SetPrintDebugProc (PrintfProcPtr print)
{
  printdebug = print;
}

int PrintDebugToFile (const char *format, ...)
{
  va_list args;

  /* initialize args */
  va_start(args,format);

  vfprintf(debugfile,format,args);
  fflush(debugfile);

  /* garbage collection */
  va_end(args);

  return (0);
}

int SetPrintDebugToFile (const char *fname)
{
  if (debugfile!=NULL)
    return (1);
  if ((debugfile=fileopen(fname,"w"))==NULL)
    return (1);

  debugfilename = StrDup(fname);
  SetPrintDebugProc(PrintDebugToFile);

  return (0);
}

int PostprocessDebugFile (const char *newname)
{
#       ifndef ModelP
  char c;

  if (debugfile==NULL)
    return (1);
  if (debugfilename==NULL)
    return (1);
  if (fclose(debugfile))
    return (1);
  if ((debugfile=fileopen(debugfilename,"r"))==NULL)
    return (1);
  if ((c=getc(debugfile))==EOF)
  {
    /* remove empty file */
    if (fclose(debugfile))
      return (1);
    if (remove(debugfilename))
      return (1);
  }
  else if (newname!=NULL)
  {
    /* rename nonemty file */
    if (fclose(debugfile))
      return (1);
    remove(newname);
    if (rename(debugfilename,newname))
      return (1);
  }
#       endif
  return (0);
}

INT PrintRepErrStack (PrintfProcPtr print)
{
  if (rep_err_count==0)
    print("no errors are reported\n");
  else
  {
    INT i;

    print("reported errors are:\n\n");

    for (i=0; i<rep_err_count; i++)
      print("%2d: File: %20s, Line: %5d\n",i,rep_err_file[i],rep_err_line[i]);
  }

  return (0);
}

/* TODO: delete this */
/*

   static int InitDebug()
   {
        SetPrintDebugProc(printf);
        return(0);
   }

   main()
   {
        char string[10]= "Hallo";
        int n=1234;
        InitDebug();
        PrintDebug("This is a string:%s\n",string);
        PrintDebug("This is a integer:%d\n",n);
   }
 */

#endif /* Debug */
