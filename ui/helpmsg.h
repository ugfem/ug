// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  helpmsg.h                                                                                                     */
/*																			*/
/* Purpose:   implements the access to online help message files			*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/*																			*/
/* History:   18.02.92 begin, ug version 2.0								*/
/*			  02.02.95 ug 3.0												*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __HELPMSG__
#define __HELPMSG__


#ifndef __COMPILER__
#include "compiler.h"
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

/* modes for PrintHelp() */
#define HELPITEM                                0       /* searching for help name				*/
#define KEYWORD                                 1       /* searching for keyword				*/

/* return codes of PrintHelp */
#define HELP_OK                                 0
#define HELP_STRING_EMPTY               1
#define HELP_NOT_FOUND                  2
#define HELP_STRING_TOO_LONG    3

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* online help mechenism */
INT     InitHelpMsg                     (char *helpfiles);
INT     PrintHelp                               (const char *HelpFor,int mode, const char *addText);
INT     CheckHelp                               (void);

#endif
