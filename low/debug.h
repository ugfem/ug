// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  debug.h														*/
/*																			*/
/* Purpose:   header file for ug internal debugger							*/
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
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __DEBUG__
#define __DEBUG__

#ifndef __GENERAL__
#include "general.h"
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

#ifdef Debug
#include <assert.h>

#define CAT(x,y)                x ## y
#define IFDEBUG(m,n)    if (CAT(Debug, m) >=(n)) {
#define PRINTDEBUG(m,n,s) IFDEBUG(m,n) PrintDebug s; ENDDEBUG
#define ENDDEBUG  }
#define RETURN(rcode)   {INT rc; rc = rcode; assert(!rc); return (rc);}
#define HEAPFAULT(ptr)  assert(((int *)ptr)[1]!=-1);
#define ASSERT(exp)             assert(exp)
#else
#define IFDEBUG(m,n)    if (1==0) {
#define ENDDEBUG        }
#define PRINTDEBUG(m,n,s) /* no debugging */
#define RETURN(rcode)   return (rcode)
#define HEAPFAULT(ptr)
#define ASSERT(exp)
#endif

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef int (*PrintDebugProcPtr)(const char *, ...);

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

#if (defined Debug && !defined compile_debug)

extern int Debuginit;
extern int Debugdddif;
extern int Debugdev;
extern int Debuggm;
extern int Debuggraph;
extern int Debuglow;
extern int Debugdom;
extern int Debugmachines;
extern int Debugnumerics;
extern int Debugui;

#endif

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

void SetPrintDebugProc          (PrintDebugProcPtr print);
void PrintDebug                         (const char *format, ...);
int  PrintDebugToFile           (const char *format, ...);
int  SetPrintDebugToFile        (const char *fname);

#endif
