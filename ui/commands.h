// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      commands.h                                                    */
/*                                                                          */
/* Purpose:   command for ug command line interface                         */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*                                                                          */
/* History:   18.02.92 begin, ug version 2.0                                */
/*            05 Sep 1992, split cmd.c into cmdint.c and commands.c         */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __COMMANDS__
#define __COMMANDS__

#include <cstdio>


#include "ugtypes.h"
#include "gm.h"
#include "namespace.h"

START_UGDIM_NAMESPACE


/** This method is not static because it is needed in DUNE */
INT NewCommand(INT argc, char **argv);

/** This method is not static because it is needed in DUNE */
INT ConfigureCommand (INT argc, char **argv);

/** This method is not static because it is needed in DUNE */
INT LBCommand (INT argc, char **argv);


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*    compile time constants defining static data size (i.e. arrays)        */
/*    other constants                                                       */
/*    macros                                                                */
/*                                                                          */
/****************************************************************************/

#define NO_OPTION_CHECK(argc,argv)      if (argc>1)                                                     \
  {                                                                               \
    UserWrite("don't specify arguments with "); \
    UserWrite(argv[0]);                             \
    UserWrite("\n");                                        \
    return (CMDERRORCODE);                          \
  }

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT              InitCommands                   (void);

FILE            *GetProtocolFile                (void);
INT                      QualityElement                 (MULTIGRID *theMG, ELEMENT *theElement);

MULTIGRID       *GetCurrentMultigrid    (void);
INT              SetCurrentMultigrid    (MULTIGRID *theMG);

END_UGDIM_NAMESPACE

#endif
