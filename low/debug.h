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
#include "misc.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define REP_ERR_MAX             10

/* if HEAPCHECK is true p is a pointer to a zombie object */
#define HEAPCHECK(ptr) (((int *)ptr)[1] == -1)

#ifdef Debug
#include <assert.h>

#define IFDEBUG(m,n)    if (Debug ## m >=(n)) {
#define PRINTDEBUG(m,n,s) IFDEBUG(m,n) PrintDebug s; ENDDEBUG
#define ENDDEBUG  }
#define RETURN(rcode)   {INT rc; rc = rcode; assert(!rc); return (rc);}
#define HEAPFAULT(ptr)  assert(((int *)ptr)[1]!=-1);
#define ASSERT(exp)             assert(exp)

#define REP_ERR_INC             {rep_err_line[rep_err_count] = __LINE__;  \
                                 rep_err_file[rep_err_count] = this_file; \
                                 rep_err_count = (rep_err_count+1)%REP_ERR_MAX;}
#define REP_ERR_RETURN(err)             { REP_ERR_INC  return (err);}
#define REP_ERR_RESET                   rep_err_count = 0;
#define REP_ERR_FILE                    static char *this_file=__FILE__
#else
#define IFDEBUG(m,n)    if (1==0) {
#define ENDDEBUG        }
#define PRINTDEBUG(m,n,s) /* no debugging */
#define RETURN(rcode)   return (rcode)
#define HEAPFAULT(ptr)
#define ASSERT(exp)

#define REP_ERR_RETURN(err)             return (err);
#define REP_ERR_INC
#define REP_ERR_RESET
#define REP_ERR_FILE
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
extern int Debugnp;
extern int Debugui;

/* for reporting of erros (using the REP_ERR_RETURN-macro) */
extern int rep_err_count;
extern int rep_err_line[REP_ERR_MAX];
extern const char  *rep_err_file[REP_ERR_MAX];

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
