// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  scan.c	                                                                                                */
/*																			*/
/* Purpose:   tools for reading script arguments                                */
/*																			*/
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   November 23, 1996                                                                         */
/*			  low part of former np/udm/scan.c, 15.5.97						*/
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "debug.h"
#include "ugenv.h"
#include "scan.h"
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

#define OPTIONLEN                       32
#define OPTIONLENSTR            "31"
#define VALUELEN                        64
#define VALUELENSTR                     "63"

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

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   ReadArgvDOUBLE - Read command strings

   SYNOPSIS:
   INT ReadArgvDOUBLE (const char *name, DOUBLE *a, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  a - DOUBLE value
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns a DOUBLE value.

   RETURN VALUE:
   INT
   .n    0 if the argument was found and a DOUBLE value could be read
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvDOUBLE (const char *name, DOUBLE *a, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  double value;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %lf",option,&value)!=2)
        continue;
      if (strcmp(option,name) == 0)
      {
        a[0] = value;
        return(0);
      }
    }

  return (1);
}

/****************************************************************************/
/*D
   ReadArgvINT - Read command strings

   SYNOPSIS:
   INT ReadArgvINT (const char *name, INT *j, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  j - integer value
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns an integer value in 'j'.

   RETURN VALUE:
   INT
   .n    0 if the argument was found and a value could be read
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvINT (const char *name, INT *j, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  int value;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %d",option,&value)!=2)
        continue;
      if (strcmp(option,name) == 0)
      {
        j[0] = value;
        return(0);
      }
    }

  return (1);
}

/****************************************************************************/
/*D
   ReadArgvChar - Read command strings

   SYNOPSIS:
   INT ReadArgvChar (const char *name, char *buffer, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  buffer - string
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns an string value in 'buffer'.

   RETURN VALUE:
   INT
   .n    0 if the argument was found and a value could be read
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvChar (const char *name, char *buffer, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  char value[VALUELEN];

  strcpy(buffer,"");
  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0]) {
      if (sscanf(argv[i],
                 expandfmt(CONCAT5("%",OPTIONLENSTR,"[a-zA-Z0-9_] %",
                                   VALUELENSTR,"[ -~]")),option,value)!=2)
        continue;
      if (strcmp(option,name) == 0) {
        strcpy(buffer,value);
        return(0);
      }
    }

  return (1);
}

/****************************************************************************/
/*D
   ReadArgvMEM - Read command strings for (memory) size specification

   SYNOPSIS:
   INT ReadArgvMEM (const char *name, MEM *mem_size, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  mem_size - integer value
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns an MEM value in 'mem_size'.
   It converts a (memory)size specification from String to
   type MEM (an integer type).
   The size specification contains an integer number followed by an
   optional unit specifier:
      G for gigabyte
      M for megabyte
      K for kilobyte
   (also the lower case char's are recognized).

   EXAMPLE:
      "10M" is converted to 10485760 (10 mega byte).

   RETURN VALUE:
   INT
   .n    0 ok (the argument was found and a value could be read)
   .n    1 integer could not be read or invalid unit specifier

   SEE ALSO:
      MEM
   D*/
/****************************************************************************/

INT ReadArgvMEM (const char *name, MEM *mem_size, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN],size_input[20];

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %s",option,size_input)!=2)
        continue;
      if (strcmp(option,name) == 0)
      {
        switch(ReadMemSizeFromString( size_input, mem_size ))
        {
        case 0 : return(0);
        case 1 : PrintErrorMessage( 'E', "ReadArgvMEM", "invalid integer number read" );
          return(1);
        case 2 : PrintErrorMessage( 'E', "ReadArgvMEM", "invalid unit specifier read" );
          return(1);
        default : assert(0);
        }
      }
    }

  return (1);
}

/****************************************************************************/
/*D
   ReadArgvOption - Read command strings

   SYNOPSIS:
   INT ReadArgvOption (const char *name, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns an integer value.

   RETURN VALUE:
   INT
   .n    0 if the option is not set
   .n    n if an integer n is given with the option
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvOption (const char *name, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  int value;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %d",option,&value) == 2)
        if (strcmp(option,name) == 0)
          return(value);
      if (strcmp(argv[i],name) == 0)
        return(1);
    }

  return (0);
}
