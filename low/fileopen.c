// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  fileopen.c													*/
/*																			*/
/* Purpose:   definition of a fopen fct. that accepts UNIX-style pathnames	*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: henrik@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7007									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   02.02.95 new for ug version 3.0								*/
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

/* standard C library */
#include <stdio.h>
#include <string.h>

/* low module */
#include "compiler.h"
#include "defaults.h"
#include "ugenv.h"

#include "fileopen.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAXPATHLENGTH           256
#define MAXPATHS                        16

#define SEPERATOR                       " \t"

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef char PATH[MAXPATHLENGTH];

typedef struct
{

  /* env item */
  ENVVAR v;

  INT nPaths;                           /* number of paths stored						*/
  PATH path[1];                         /* begin of path list							*/

} PATHS;

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static INT thePathsDirID;
static INT thePathsVarID;

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* Function:  MakePathsItem													*/
/*																			*/
/* Purpose:   create a Paths environment item								*/
/*																			*/
/* Input:	  name of the Paths, number of Paths	                                                */
/*																			*/
/* Output:	  PATHS *: pointer to the Paths struct							*/
/*			  NULL: if an error occured                                                                     */
/*																			*/
/****************************************************************************/

static PATHS *MakePathsItem (const char *name, INT nPaths)
{
  if (ChangeEnvDir("/Paths") == NULL) return (NULL);
  if (strlen(name)>=NAMESIZE || strlen(name)<=1) return (NULL);

  return ((PATHS *) MakeEnvItem(name,thePathsVarID,sizeof(PATHS)+(nPaths-1)*sizeof(PATH)));
}

/****************************************************************************/
/*																			*/
/* Function:  GetPaths														*/
/*																			*/
/* Purpose:   find the Paths environment item with name                                         */
/*																			*/
/* Input:	  name of the Paths to find                                                                     */
/*																			*/
/* Output:	  PATHS *: pointer to the multigrid struct						*/
/*			  NULL: if an error occured                                                                     */
/*																			*/
/****************************************************************************/

static PATHS *GetPaths (const char *name)
{
  return ((PATHS *) SearchEnv(name,"/Paths",thePathsVarID,thePathsDirID));
}

/****************************************************************************/
/*D
        fileopen - like ANSI-C 'fopen' but convert UNIX-style paths to machine format

        SYNOPSIS:
        FILE *fileopen (const char *fname, const char *mode)

    PARAMETERS:
   .   fname - filename with path convention in UNIX-style
   .   mode - see ANSI-C 'fopen'

        DESCRIPTION:
        Convert UNIX-style paths to machine format. If '__MWCW__' or '__MPW32__'
        are define (Macintosh compilers) the UNIX-sytle path is converted to
        Macintosh-style. Then standard 'fopen' is called.

        For other platforms 'fileopen' is just mapped to 'fopen'.

        RETURN VALUE:
        FILE *
   .n   file ptr (NULL if error)

        SEE ALSO:
        fopen
   D*/
/****************************************************************************/

#if (defined __MWCW__) || (defined __MPW32__)

/* Macintosh computers */

FILE *fileopen (const char *fname, const char *mode)
{
  char fullpath[MAXPATHLENGTH];
  int pos;

  if (*fname=='/')
    /* root not defined on Macintosh computers */
    return (NULL);
  if (*fname=='~')
    /* home directory not defined on Macintosh computers */
    return (NULL);

  /* something to convert? */
  if (strchr(fname,'/')==NULL)
    return (fopen(fname,mode));

  pos = 0;
  while ((*fname!='\0') && (pos<MAXPATHLENGTH-2))
    switch (fname[0])
    {
    case '.' :
      if ((fname[1]=='.')&&(fname[2]=='/'))
      {
        /* "../" */

        /* if the path starts with "../" we interpret "./../", i.e. "::" */
        if (pos==0)
          fullpath[pos++] = ':';

        fullpath[pos++] = ':';
        fname += 3;
      }
      else if (fname[1]=='/')
      {
        /* "./" --> ":" */
        fullpath[pos++] = ':';
        fname += 2;
      }
      else
        fullpath[pos++] = *(fname++);
      break;

    case '/' :
      /* "/" --> ":" */
      fullpath[pos++] = ':';
      while (*fname=='/') fname++;
      break;

    default :
      fullpath[pos++] = *(fname++);
    }

  if (pos>=MAXPATHLENGTH)
    /* filename too long */
    return (NULL);

  /* 0-terminate string */
  fullpath[pos] = '\0';

  return (fopen(fullpath,mode));
}

#else

/* UNIX machines */

FILE *fileopen (const char *fname, const char *mode)
{
  return (fopen(fname,mode));
}

#endif

/****************************************************************************/
/*D
        ReadSearchingPaths - read searching paths from a defaults file

        SYNOPSIS:
        INT ReadSearchingPaths (const char *filename, const char *paths)

    PARAMETERS:
   .   filename - search paths in a defaults file with this name (most likely to be
                        `the` --> 'defaults' file, use 'DEFAULTSFILENAME')
   .   paths - name of the paths item looked for in a defaults file

        DESCRIPTION:
        From a defaultsfile using --> 'GetDefaultValue' the specified paths item
        is read containing one or more paths sperated by blanks (tab or space).
        The paths are stored in an environment item with that same name in the
        environment directory '/Paths'. The function --> 'FileOpenUsingSearchPaths'
        is looking up the paths to be tryed there.

        RETURN VALUE:
        INT
   .n   1: failed to 'GetDefaultValue'
   .n   2: more than 'MAXPATHS' specified in the defaults file
   .n   3: failed to 'MakePathsItem'
   .n   0: ok

        SEE ALSO:
        FileOpenUsingSearchPaths
   D*/
/****************************************************************************/

INT ReadSearchingPaths (const char *filename, const char *paths)
{
  PATHS *thePaths;
  INT i,nPaths;
  char *Path[MAXPATHS];
  char *token,buffer[BUFFLEN];

  if (GetDefaultValue(filename,paths,buffer)!=0)
    return (1);

  /* get Paths */
  nPaths = 0;
  token = strtok(buffer,SEPERATOR);
  while (token!=NULL)
  {
    if (nPaths>=MAXPATHS)
      return (2);                       /* too many paths */

    Path[nPaths++] = token;
    token = strtok(NULL,SEPERATOR);
  }

  /* create env item */
  if ((thePaths=MakePathsItem(paths,nPaths))==NULL)
    return (3);

  /* fill data */
  thePaths->nPaths = nPaths;
  for (i=0; i<nPaths; i++)
    strcpy(thePaths->path[i],Path[i]);

  return (0);
}

/****************************************************************************/
/*D
        FileOpenUsingSearchPaths - open file searching in the directories specified
                        in the environment item '/Paths/<paths>'

        SYNOPSIS:
        FILE *FileOpenUsingSearchPaths (const char *fname, const char *mode, const char *paths)

    PARAMETERS:
   .   fname - file name to be opened
   .   mode - see ANSI-C 'fopen'
   .   paths - try paths specified in the environment item '/Paths/<paths> which was
                        set by --> 'ReadSearchingPaths'

        DESCRIPTION:
        The functions trys to open the file with 'filename' using one by one the
        paths specified in the environment item '/Paths/<paths> which was
        set by --> 'ReadSearchingPaths'. It is used in several places in ug (all paths
        are read from the standard --> 'defaults' file)":"

   .n   'srciptpaths' is used by the interpreter for script execution
   .n   'gridpaths' is used by ugio to read grids from (they are stored in the
   .n   first path

        RETURN VALUE:
        FILE *
   .n   pointer to file opened, 'NULL' if error

        SEE ALSO:
        ReadSearchingPaths, fileopen
   D*/
/****************************************************************************/

FILE *FileOpenUsingSearchPaths (const char *fname, const char *mode, const char *paths)
{
  PATHS *thePaths;
  FILE *theFile;
  char fullname[MAXPATHLENGTH];
  INT i,fnamelen;

  fnamelen = strlen(fname);

  if ((thePaths=GetPaths(paths))==NULL)
    return (NULL);

  for (i=0; i<thePaths->nPaths; i++)
  {
    if (strlen(thePaths->path[i])+fnamelen>MAXPATHLENGTH)
      return (NULL);

    strcpy(fullname,thePaths->path[i]);
    strcat(fullname,"/");
    strcat(fullname,fname);

    if ((theFile=fileopen(fullname,mode))!=NULL)
      return (theFile);
  }

  return (NULL);
}

/****************************************************************************/
/*D
        FileOpenUsingSearchPath - try to open a file in the specified path

        SYNOPSIS:
        FILE *FileOpenUsingSearchPath (const char *fname, const char *mode, const char *path)

    PARAMETERS:
   .   fname - open file with this name
   .   mode - see ANSI-C 'fopen'
   .   path - path to which fname is to be appended

        DESCRIPTION:
        Try to open a file in the specified path.

        RETURN VALUE:
        FILE *
   .n   pointer to file opened, 'NULL' if error

        SEE ALSO:
        FileOpenUsingSearchPaths, fileopen
   D*/
/****************************************************************************/

FILE *FileOpenUsingSearchPath (const char *fname, const char *mode, const char *path)
{
  FILE *theFile;
  char fullname[MAXPATHLENGTH];

  if (strlen(path)+strlen(fname)>MAXPATHLENGTH)
    return (NULL);

  strcpy(fullname,path);
  strcat(fullname,"/");
  strcat(fullname,fname);

  if ((theFile=fileopen(fullname,mode))!=NULL)
    return (theFile);

  return (NULL);
}

/****************************************************************************/
/*D
        InitFileOpen - init 'fileopen.c'

        SYNOPSIS:
        INT InitFileOpen ()

    PARAMETERS:
    --

        DESCRIPTION:
        An environment directory '/Paths' is created where the paths read by
        'ReadSearchingPaths' are stored.

        RETURN VALUE:
        INT
   .n   __LINE__: could not create '/Paths'
   .n           0: ok
   D*/
/****************************************************************************/

INT InitFileOpen ()
{
  /* install the /Paths directory */
  if (ChangeEnvDir("/")==NULL)
    return(__LINE__);

  thePathsDirID = GetNewEnvDirID();
  if (MakeEnvItem("Paths",thePathsDirID,sizeof(ENVDIR))==NULL)
    return(__LINE__);

  thePathsVarID = GetNewEnvVarID();

  return (0);
}
