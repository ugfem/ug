// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      defaults.c                                                    */
/*                                                                          */
/* Purpose:   access to defaults file                                       */
/*                                                                          */
/* Author:      Peter Bastian                                               */
/*              Institut fuer Computeranwendungen III                       */
/*              Universitaet Stuttgart                                      */
/*              Pfaffenwaldring 27                                          */
/*              70569 Stuttgart                                             */
/*              email: ug@ica3.uni-stuttgart.de                             */
/*                                                                          */
/* History:   17.12.94 begin, ug version 3.0                                */
/*                                                                          */
/* Revision:  06.09.95                                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*              system include files                                        */
/*              application include files                                   */
/*                                                                          */
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "fileopen.h"
#include "general.h"
#include "misc.h"
#include "defaults.h"
#ifdef ModelP
#include "ppif.h"
#endif


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

#ifdef ModelP
static char *defaults_buffer=NULL;
#endif


/****************************************************************************/
/*D
   GetDefaultValue - Provide access to defaults file

   SYNOPSIS:
   INT GetDefaultValue (const char *filename, const char *name, char *value);

   PARAMETERS:
   .  filename - pointer to char (const)
   .  name -     pointer to char (const)
   .  value -    pointer to char

   DESCRIPTION:
   This function provides access to defaults file. When 'ug' is started
   this function reads the defaults on file in order to set some
   parameters of 'ug' in advance.

   RETURN VALUE:
   INT
   .n    0 if OK
   .n    1 if error in opening or reading defaults file
   D*/
/****************************************************************************/

#define NAMESIZE    32
#define NAMELEN     31
#define NAMELENSTR    "31"

INT GetDefaultValue (const char *filename, const char *name, char *value)
{
  FILE *defaultsFile;
  char Name[NAMESIZE], buffer[BUFFSIZE];

        #ifdef ModelP
  int file_ok;
  size_t fsize,actsize;
  char *curr;

  if (defaults_buffer==NULL)
  {
    /* first call to GetDefaultValue, we must read the defaults file */

    /* get filesize and broadcast it */
    if (me==master) fsize = filesize(filename);
    Broadcast(&fsize, sizeof(fsize));
    if (fsize==0) return(1);

    /* get buffer for file */
    defaults_buffer = (char *) malloc(fsize+1);
    assert(defaults_buffer!=NULL);

    /* open defaults file */
    if (me==master) {
      defaultsFile = fileopen(filename,"r");
      file_ok = (defaultsFile!=NULL);
    }
    Broadcast(&file_ok, sizeof(file_ok));
    if (!file_ok) {
      free(defaults_buffer);
      defaults_buffer = NULL;
      return(1);
    }


    /* read file into buffer */
    if (me==master) {
      actsize = fread(defaults_buffer, 1, fsize, defaultsFile);
      fclose(defaultsFile);

      /* set end mark */
      defaults_buffer[actsize] = 0;
    }

    /* broadcast buffer */
    Broadcast(defaults_buffer, fsize);
  }

  curr = defaults_buffer;
  while (curr!=NULL && *curr!=0)
  {
    if (sscanf(curr,
               expandfmt(CONCAT5(" %",NAMELENSTR,"[0-9a-zA-Z_] %",
                                 BUFFLENSTR,"[ -~]")), Name,value) ==2)
    {
      if (strcmp(Name,name)==0) return(0);
    }

    /* get next line */
    curr = strchr(curr, '\n');
    if (curr!=NULL) curr++;
  }

  return(1);


        #else

  /* sequential version */
  defaultsFile = fileopen(filename,"r");
  if (defaultsFile==NULL) return(1);

  rewind(defaultsFile);
  while (fgets(buffer,255,defaultsFile)!=NULL)
  {
    if (sscanf(buffer,
               expandfmt(CONCAT5(" %",NAMELENSTR,"[0-9a-zA-Z_] %",
                                 BUFFLENSTR,"[ -~]")),Name,value) ==2)
    {
      if (strcmp(Name,name)==0) {
        fclose(defaultsFile);
        return(0);
      }
    }
  }
  fclose(defaultsFile);
  return(1);

        #endif
}
