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

#ifdef Debug

#include <stdio.h>
#include <stdarg.h>


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

typedef int (*PrintDebugProcPtr)(const char *, ...);

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

int Debuginit           =       0;
int Debugdddif              =       6;      /* temporary setting for debugging ModelP */
int Debugdev                =       0;
int Debuggm                 =       0;
int Debuggraph              =       0;
int Debuglow                =       0;
int Debugmachines   =       0;
int Debugnumerics   =       0;
int Debugui                 =       2;      /* temporary setting for debugging ModelP */


/* from dddif/ppif.h */
extern int me, master;


/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

PrintDebugProcPtr printdebug=printf;

/* data for CVS */
static char rcsid[] = "$Header$";

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

void PrintDebug (const char *format, ...)
{
  char buffer[256];
  va_list args;

  va_start(args,format);

  vsprintf(buffer,format,args);
        #ifdef ModelP
  if (me==master) {
        #endif

  /* use specific debug function for displaying */
  (*printdebug)(buffer);

        #ifdef ModelP
}
else
{
  printf(buffer);
  fflush(stdout);
}

        #endif

  va_end(args);
}

void SetPrintDebugProc (PrintDebugProcPtr print)
{
  printdebug = print;
}


int InitDebug()
{
  SetPrintDebugProc(printf);
  return(0);
}


/* TODO: delete this */
/*
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
