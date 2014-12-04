// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      scan.c                                                        */
/*                                                                          */
/* Purpose:   tools for reading script arguments                            */
/*                                                                          */
/*                                                                          */
/* Author:    Christian Wieners                                             */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/* email:     ug@ica3.uni-stuttgart.de                                      */
/*                                                                          */
/* History:   November 23, 1996                                             */
/*            low part of former np/udm/scan.c, 15.5.97                     */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/* system include files                                                     */
/* application include files                                                */
/*                                                                          */
/****************************************************************************/

#include <config.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "debug.h"
#include "ugenv.h"
#include "misc.h"
#include "scan.h"

USING_UG_NAMESPACE

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
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/* in the corresponding include file!)                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

  REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/** \brief Read command line doubles

   \param name - name of the argument
   \param a - DOUBLE value
   \param argc - argument counter
   \param argv - argument vector

   This function reads command strings and returns a DOUBLE value.

   \return <ul>
   <li>  0 if the argument was found and a DOUBLE value could be read </li>
   <li>  1 else. </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX ReadArgvDOUBLE (const char *name, DOUBLE *a, INT argc, char **argv)
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
/** \brief Read command line integers

   \param name - name of the argument
   \param j - integer value
   \param argc - argument counter
   \param argv - argument vector

   This function reads command strings and returns an integer value in 'j'.

   \return <ul>
   <li>  0 if the argument was found and a value could be read </li>
   <li>  1 else. </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX ReadArgvINT (const char *name, INT *j, INT argc, char **argv)
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
/** \brief Read DOUBLE and INT from an option

   \param name - name of the argument
   \param a - DOUBLE value
   \param m - INT value
   \param argc - argument counter
   \param argv - argument vector

   This function reads command strings and sets a DOUBLE value and an
   optional INT-Value.

   \return <ul>
   <li>  0 - !!!! no argument was found !!!! </li>
   <li>  1,2 - number of arguments </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX ReadArgvDOUBLE_INT (const char *name, DOUBLE *a, INT *j, INT argc, char **argv)
{
  INT i,k;
  char option[OPTIONLEN];
  double value;
  int ivalue;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if ((k=sscanf(argv[i],"%s %lf %d",option,&value,&ivalue))<2)
        continue;
      if (strcmp(option,name) == 0)
      {
        *a = value;
        if (k==3) *j = ivalue;else *j=0;
        return(k-1);
      }
    }

  return (0);
}

/****************************************************************************/
/** \brief Read command strings

   \param name - name of the argument
   \param buffer - string
   \param argc - argument counter
   \param argv - argument vector

   This function reads command strings and returns an string value in 'buffer'.

   \return <ul>
   <li>  0 if the argument was found and a value could be read </li>
   <li>  1 else. </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX ReadArgvChar (const char *name, char *buffer, INT argc, char **argv)
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
/** \brief Read command strings for (memory) size specification

   \param name - name of the argument
   \param mem_size - integer value
   \param argc - argument counter
   \param argv - argument vector

   This function reads command strings and returns an MEM value in 'mem_size'.
   It converts a (memory)size specification from String to
   type MEM (an integer type).
   The size specification contains an integer number followed by an
   optional unit specifier:
      G for gigabyte
      M for megabyte
      K for kilobyte
   (also the lower case chars are recognized).

   EXAMPLE:
      "10M" is converted to 10485760 (10 mega byte).

   \return <ul>
   <li>  0 ok (the argument was found and a value could be read) </li>
   <li>  1 integer could not be read or invalid unit specifier </li>
   </ul>

   \sa
      MEM
 */
/****************************************************************************/

INT NS_PREFIX ReadArgvMEM (const char *name, MEM *mem_size, INT argc, char **argv)
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
        case 1 : return(1);
        case 2 : return(1);
        default : ASSERT(0);
        }
      }
    }

  return (1);
}

/****************************************************************************/
/** \brief Read command line options

   \param name - name of the argument
   \param argc - argument counter
   \param argv - argument vector

   This function reads command strings and returns an integer value.

   \return <ul>
   <li>  0 if the option is not set </li>
   <li>  n if an integer n is given with the option </li>
   <li>  1 else. </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX ReadArgvOption (const char *name, INT argc, char **argv)
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
