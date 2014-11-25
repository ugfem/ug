// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  cmdint.h														*/
/*																			*/
/* Purpose:   command interpreter header file								*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*																			*/
/* History:   18.02.92 begin, ug version 2.0								*/
/*			  05 Sep 1992, split cmd.c into cmdint.c and commands.c                 */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#include <cstdio>

#ifndef __CMDINT__
#define __CMDINT__


#include "ugtypes.h"
#include "ugenv.h"
#include "cmdline.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define PROMPT          "> "            /* prompt string printed on each input line */

/* return values for InterpretCommand */
#define NOTFOUNDCODE    2               /* command not found						*/
#define NOSPACECODE     3               /* not enough space to execute properly         */
#define CMDINTBUFSIZE   32000


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

extern INT cmdintbufsize;

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

void    CommandLoop                     (int argc, char **argv);
#ifdef ModelP
void    ParCommandLoop                  (char *inpLine);
#endif
void    SetDoneFlag                     (void);
void    ResetDoneFlag                   (void);
int             GetDoneFlag                     (void);
FILE   *FOpenScript                             (const char *script, const char *mode);

INT     InterpretCommand                (const char *cmds);

INT     InitCommandInterpreter  (INT argc, char **argv);
INT     ExitCommandInterpreter       ();

END_UGDIM_NAMESPACE

#endif
