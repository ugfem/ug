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

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __CMDINT__
#define __CMDINT__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

#ifndef __CMDLINE__
#include "cmdline.h"
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

#define PROMPT          "> "            /* prompt string printed on each input line */

/* return values for InterpretCommand */
#define NOTFOUNDCODE    2               /* command not found						*/
#define NOSPACECODE     3               /* not enough space to execute properly         */
#define MAXCMDSIZE  2048



/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

void    CommandLoop                     (int argc, char **argv);
void    SetDoneFlag                     (void);

INT     InterpretCommand                (char *cmds);

INT     InitCommandInterpreter  (void);

#endif
