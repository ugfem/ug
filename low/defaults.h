// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      defaults.h                                                    */
/*                                                                          */
/* Purpose:   implements access to defaults file                            */
/*                                                                          */
/* Author:      Peter Bastian,                                              */
/*              Institut fuer Computeranwendungen III                       */
/*              Universitaet Stuttgart                                      */
/*              Pfaffenwaldring 27                                          */
/*              70569 Stuttgart                                             */
/*              email: ug@ica3.uni-stuttgart.de                          */
/*                                                                          */
/* History:   18.02.92 begin, ug version 2.0                                */
/*              05 Sep 1992, split cmd.c into cmdint.c and commands.c       */
/*              17.12.94 ug 3.0                                             */
/*                                                                          */
/* Revision:  06.09.95                                                      */
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

#ifndef __DEFAULTS__
#define __DEFAULTS__


#ifndef __COMPILER__
#include "compiler.h"
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

#define DEFAULTSFILENAME            "defaults"
#define UGRC_NAME                           ".ugrc"

#define BUFFSIZE    256             /* max length of name string    */
#define BUFFLEN     255             /* BUFFSIZE-1                   */
#define BUFFLENSTR "255"            /* BUFFSIZE-1 as string         */

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* access to default file */
INT     GetDefaultValue                         (const char *filename, const char *name, char *value);
INT     GetLocalizedDefaultValue        (const char *filename, const char *name, char *value);

#endif
