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
#include "debug.h"
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
   GetLocalizedDefaultValue - Provide access to defaults file

   SYNOPSIS:
   INT GetLocalizedDefaultValue (const char *filename, const char *name, char *value);

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

   SEE ALSO:
   GetDefaultValue
   D*/
/****************************************************************************/

#define NAMESIZE    32
#define NAMELEN     31
#define NAMELENSTR    "31"

INT GetLocalizedDefaultValue (const char *filename, const char *name, char *value)
{
  FILE *defaultsFile;
  char Name[NAMESIZE];

        #ifdef ModelP
  int file_ok;
  size_t fsize,actsize;
  char *curr;
  static char *buffered_filename;

  if (defaults_buffer==NULL)
  {
    /* first call to GetLocalizedDefaultValue, we must read the defaults file */
    if (filename==NULL)
      return(1);

    buffered_filename = StrDup(filename);

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
  else
  {
    if (filename!=NULL && strcmp(buffered_filename,filename)!=0)
      return(1);
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
  char buffer[BUFFSIZE];

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

   SEE ALSO:
   GetLocalizedDefaultValue
   D*/
/****************************************************************************/


#define MAX_PATH_LEN            1024
enum {PATH_TOO_LONG = 1, COULD_NOT_STAT};

static INT GetPathedUGRCValue (const char *path, const char *name, char *value, int *err)
{
  char full_path[MAX_PATH_LEN];

  if (strlen(path)+strlen(UGRC_NAME)+2 >= MAX_PATH_LEN)
    return PATH_TOO_LONG;

  strcpy(full_path,path);
  AppendTrailingSlash(full_path);
  strcat(full_path,UGRC_NAME);

  PRINTDEBUG(low,2,("GetPathedUGRCValue: trying: '%s'\n",full_path));
  if (filetype(full_path)==FT_FILE)
  {
    *err = GetLocalizedDefaultValue(full_path,name,value);
    return 0;
  }
  else
  {
    PRINTDEBUG(low,2,("GetPathedUGRCValue: could not stat: '%s'\n",full_path));
    return COULD_NOT_STAT;
  }
}

INT GetDefaultValue (const char *filename, const char *name, char *value)
{
        #ifdef ModelP
  static int already_called = FALSE;

  if (already_called)
    return GetLocalizedDefaultValue(NULL,name,value);
  else
    already_called = TRUE;
        #endif

  PRINTDEBUG(low,2,("GetDefaultValue\n"));

  if (strchr(filename,'/')!=NULL)
  {
    PRINTDEBUG(low,2,("GetDefaultValue: GetLocalizedDefaultValue called directly with: '%s'\n",filename));
    return GetLocalizedDefaultValue(filename,name,value);
  }
  else if (strcmp(filename,DEFAULTSFILENAME)==0)
  {
    if (filetype(filename)==FT_FILE)
    {
      PRINTDEBUG(low,2,("GetDefaultValue: GetLocalizedDefaultValue called with defaults: '%s'\n",filename));
      return GetLocalizedDefaultValue(filename,name,value);
    }
  }
  else if (strcmp(filename,UGRC_NAME)!=0)
    ASSERT(FALSE);                      /* try GetLocalizedDefaultValue */

  /* localize defaults file */
  {
    int err;

    const char *path = getenv("HOME");
    if (path!=NULL)
      if (GetPathedUGRCValue(path,name,value,&err)==0)
      {
        PRINTDEBUG(low,2,("GetDefaultValue: GetPathedUGRCValue called for HOME='%s', err=%d\n",path,err));
        return err;
      }
    path = getenv("UGROOT");
    if (path!=NULL)
    {
      char data_path[MAX_PATH_LEN];
      strcpy(data_path,path);
      AppendTrailingSlash(data_path);
      strcat(data_path,"lib/ugdata");
      if (GetPathedUGRCValue(data_path,name,value,&err)==0)
      {
        PRINTDEBUG(low,2,("GetDefaultValue: GetPathedUGRCValue called for DATA='%s', err=%d\n",data_path,err));
        return err;
      }
    }
  }
  return 1;
}
