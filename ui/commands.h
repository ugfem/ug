// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  commands.h													*/
/*																			*/
/* Purpose:   command for ug command line interface                                             */
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

#ifndef __COMMANDS__
#define __COMMANDS__

#include <stdio.h>


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DEVICESH__
#include "devices.h"
#endif

#ifndef __CMDLINE__
#include "cmdline.h"
#endif

#ifndef __GM__
#include "gm.h"
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

#define NO_OPTION_CHECK(argc,argv)      if (argc>1)                                                     \
  {                                                                               \
    UserWrite("don't specify arguments with "); \
    UserWrite(argv[0]);                             \
    UserWrite("\n");                                        \
    return (CMDERRORCODE);                          \
  }

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

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

INT              InitCommands                   ();

FILE            *GetProtocolFile                ();

MULTIGRID       *GetCurrentMultigrid    ();
INT              SetCurrentMultigrid    (MULTIGRID *theMG);

#endif
