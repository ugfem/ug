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
/*			  email: ug@ica3.uni-stuttgart.de							    */
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <assert.h>

/* first compiler header for __MACINTOSH__ definition iff */
#include "compiler.h"

/* includes for filesize(), filetype() */
#ifdef __MACINTOSH__
#include <unistd.h>
#include <stat.h>
#include <Files.h>
/* NB: On Macs the structs of <types.h> are defined locally in <stat.h> */
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif


/* low module */
/* avoid assertion in REP_ERR_RETURN for ModelP and GUI */
#if defined(ModelP) && defined(__GUI__)
#define PARALLEL
#undef ModelP
#endif

#include "debug.h"

#ifdef PARALLEL
#define ModelP
#undef PARALLEL
#endif

#include "defaults.h"
#include "general.h"
#include "ugenv.h"
#include "ugdevices.h"

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

#define MAX_PATH_LEN            1024
#define BASE_PATH_SIZE          512

#ifndef __MACINTOSH__
        #define ConvertUNIX_2_MachinePath(fname)                fname
        #define ConvertMachine_2_UNIXPath(fname)                fname
#endif

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

  INT nPaths;                                   /* number of paths stored						*/
  PATH path[1];                         /* begin of path list							*/

} PATHS;

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static INT thePathsDirID;
static INT thePathsVarID;


static char OldBasePath[BASE_PATH_SIZE];
static char BasePath[BASE_PATH_SIZE] = "./";

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

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
        ConvertUNIX_2_MachinePath - convert UNIX-style file paths to machine format

        SYNOPSIS:
        const char *ConvertUNIX_2_MachinePath (const char *fname)


    PARAMETERS:
   .   fname - filename with path convention in UNIX-style

        DESCRIPTION:
        Convert UNIX-style paths to machine format. If '__MACINTOSH__'
        is defined (Apple Macintosh computer) the UNIX-sytle path is converted to
        Macintosh-style.

        For other platforms 'ConvertUNIX_2_MachinePath' returns 'fname' (macro).

        RETURN VALUE:
        char *
   .n   converted file name (NULL if error)
   D*/
/****************************************************************************/

#ifdef __MACINTOSH__

/* Macintosh computers */
static const char *ConvertUNIX_2_MachinePath (const char *fname)
{
  static char fullpath[MAXPATHLENGTH];
  int pos;

  if (*fname=='/')
    /* root not defined on Macintosh computers */
    return (NULL);
  if (*fname=='~')
    /* home directory not defined on Macintosh computers */
    return (NULL);

  /* something to convert? */
  if (strchr(fname,'/')==NULL)
    return (fname);

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
        if (pos)
        {
          /* eat ./ */
          fname += 2;
          continue;
        }

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

  return (fullpath);
}

static const char* ConvertMachine_2_UNIXPath (const char* fname)
{
  static char fullpath[MAXPATHLENGTH];
  int pos;

  /* something to convert? */
  if (strchr(fname,':')==NULL)
    return (fname);

  pos = 0;
  while ((*fname!='\0') && (pos<MAXPATHLENGTH-2))
    switch (fname[0])
    {
    case ':' :
      if (fname[1]==':')
      {
        /* "::" */

        fullpath[pos++] = '.';
        fullpath[pos++] = '.';
        fullpath[pos++] = '/';
        fname++;
      }
      else
      {
        /* ":" */
        if (fullpath[pos]!='/')
          fullpath[pos++] = '/';
        fname++;
      }
      break;

    default :
      fullpath[pos++] = *(fname++);
    }

  if (pos>=MAXPATHLENGTH)
    /* filename too long */
    return (NULL);

  /* 0-terminate string */
  fullpath[pos] = '\0';

  return (fullpath);
}

#endif

const char *BasedConvertedFilename (const char *fname)
/* NOTE: once a filename has passed through BasedConvertedFilename() it is forbidden to
                 call BasedConvertedFilename() a second time for this filename due to
                 static result string.
 */
{
  PRINTDEBUG(low,2,("BasedConvertedFilename: fname= '%s'\n",fname));
  if (fname[0]!='/' && fname[0]!='~')                   /* use BasePath only if no absolute path specified */
  {
    static char based_filename[MAXPATHLENGTH];

    assert(fname!=based_filename);              /* avoid that the result of a previous call to BasedConvertedFilename() is fed back;
                                                                                   in this case the new result would interfere with the given input name and corrupt the result
                                                                                   (for example only ./ would be returned)*/
    strcpy(based_filename,BasePath);
    strcat(based_filename,fname);
    SimplifyPath(based_filename);
    PRINTDEBUG(low,1,("BasedConvertedFilename: based_filename= '%s'\n",based_filename));
    return ConvertUNIX_2_MachinePath(based_filename);
  }
  else
  {
    PRINTDEBUG(low,1,("BasedConvertedFilename: filename not based= '%s'\n",fname));
    return ConvertUNIX_2_MachinePath(fname);
  }
}

static int rename_if_necessary( const char *fname, int do_rename)
{
  FILE *f;

  if (do_rename && (f=fopen(fname,"r"))!=NULL)
  {
    time_t Time;
    struct stat fstat;
    char new_fname[128];

    fclose(f);

    strcpy(new_fname,fname);
    strcat(new_fname,".");

    if (stat(fname, &fstat)<0)
      return(1);
    Time = fstat.st_mtime;
    strftime(new_fname+strlen(fname)+1,64,"%y%m%d%H%M%S",localtime(&Time));

    if (rename(fname,new_fname)!=0)
      return(1);
  }
  return (0);           /* ok */
}


/****************************************************************************/
/*D
        mkdir_r - create a directory incl. renaming option

        SYNOPSIS:
        int mkdir_r (const char *fname, mode_t mode, int do_rename)

    PARAMETERS:
   .   fname - directory name with path convention in UNIX-style; may not be passed through BasedConvertedFilename()
   .   mode - creation mode in UNIX-style (see mkdir(2))
   .   do_rename - if TRUE an already existing subdirectory will be renamed

        DESCRIPTION:
    This function creates an directory and renames an already existing one
        instead of overwriting it if specified.
        fname may not be passed through BasedConvertedFilename().

        RETURN VALUE:
        int
   .n   0 sucessfull completion
   .n      != 0 error occured

        SEE ALSO:
        mkdir(2), fopen_r
   D*/
/****************************************************************************/

int mkdir_r (const char *fname, mode_t mode, int do_rename)
{
  const char *converted_name = BasedConvertedFilename(fname);

  if (do_rename)
  {
    if (rename_if_necessary( converted_name, do_rename)!=0)
      return (1);

    return mkdir(converted_name,mode);
  }
  else
  {
    switch (filetype(fname))                    /* filetype needs an NOT BasedConvertedFilename'ed filename */
    {
    case FT_UNKNOWN :                           /* file doesn't exist, thus create it */
      return mkdir(converted_name,mode);

    case FT_DIR :
      return 0;                                         /* OK, directory exists already */

    case FT_FILE :
      UserWriteF("mkdir_r(): file %s exists already as ordinary file; can't create directory with same name.\n",converted_name);
      return 1;

    case FT_LINK :
      UserWriteF("mkdir_r(): file %s exists already as a link; can't create directory with same name.\n",converted_name);
      return 1;

    default :
      UserWriteF("mkdir_r(): unknown file type %d for file %s\n",filetype(fname),converted_name);
      return 1;
    }
  }
}

/****************************************************************************/
/*D
        fopen_r - create a file incl. renaming option

        SYNOPSIS:
        int fopen_r (const char *fname, const char *mode, int do_rename)

    PARAMETERS:
   .   fname - file name with path convention in UNIX-style
   .   mode - file opening mode in UNIX-style (see fopen(2))
   .   do_rename - if TRUE an already existing file will be renamed

        DESCRIPTION:
    This function opens a file and renames an already existing one
        instead of overwriting it if specified.

        RETURN VALUE:
        int
   .n   0 sucessfull completion
   .n      != 0 error occured

        SEE ALSO:
        fopen(2), mkdir_r
   D*/
/****************************************************************************/

FILE *fopen_r (const char *fname, const char *mode, int do_rename)
{
  if (rename_if_necessary( fname, do_rename)!=0)
    return (NULL);

  return fopen(fname,mode);
}

/****************************************************************************/
/*D
        filesize - get size of a file with given name

        SYNOPSIS:
        size_t filesize (const char *fname)

    PARAMETERS:
   .   fname - filename with path convention in UNIX-style

        DESCRIPTION:
    This function returns the size of the given file or 0 if an error occurs.

        RETURN VALUE:
        size_t
   .n      file size (0 if error)

        SEE ALSO:
        fopen, fileopen, filetype
   D*/
/****************************************************************************/

size_t filesize (const char *fname)
{
  struct stat fstat;

  /* get (Unix) file descriptor */
  PRINTDEBUG(low,1,("filesize\n"));
  if (stat(BasedConvertedFilename(fname), &fstat)<0)
    return(0);

  return((size_t)fstat.st_size);
}

/****************************************************************************/
/*D
        filetype - get type of a file with given name

        SYNOPSIS:
        int filetype (const char *fname)

    PARAMETERS:
   .   fname - filename with path convention in UNIX-style; it may not be passed through BasedConvertedFilename()

        DESCRIPTION:
    This functon returns the type of the given file
    or FT_UNKNOWN if an error occurs.
        fname may not be passed through BasedConvertedFilename().

        RETURN VALUE:
        int
   .n      file type (one of FT_UNKNOWN, FT_FILE, FT_DIR, FT_LINK)

        SEE ALSO:
        fopen, fileopen, filesize
   D*/
/****************************************************************************/

int filetype (const char *fname)
{
  struct stat fstat;
  int r;

  /* get Unix file descriptor */
  PRINTDEBUG(low,1,("filetype\n"));
  if ((r=stat(BasedConvertedFilename(fname), &fstat))<0)
    return(FT_UNKNOWN);

        #ifdef __CC__
  switch (fstat.st_mode & _S_IFMT)
  {
  case _S_IFREG :   return FT_FILE;
  case _S_IFDIR :   return FT_DIR;
#ifdef S_IFLNK
  case _S_IFLNK :   return FT_LINK;
#endif
  }
#else
  switch (fstat.st_mode & S_IFMT)
  {
  case S_IFREG :   return FT_FILE;
  case S_IFDIR :   return FT_DIR;
#ifdef S_IFLNK
  case S_IFLNK :   return FT_LINK;
#endif
  }
#endif
  return(FT_UNKNOWN);
}

/****************************************************************************/
/*D
        DirWalk - loop the names of files in a directory and call
                                a ProcessFileProc for each

        SYNOPSIS:
        INT DirWalk (const char *dir, ProcessFileProc fcn)

    PARAMETERS:
   .   dir - UNIX-style path
   .   fcn - will be called with UNIX-style full path for each file in dir

        DESCRIPTION:
    This functon loops the names of files in a directory and calls
        a user specified ProcessFileProc for each.

        RETURN VALUE:
        INT
   .n      0: success
   .n      PATH_INVALID: dir is no path to a directory
   .n      NAME_TOO_LONG: full file path exceeds buffer length
   D*/
/****************************************************************************/

#ifdef ModelP
#define PARALLEL
#undef ModelP
#endif

INT DirWalk (const char *dir, ProcessFileProc fcn)
{

  /* encapsulate implementation dependent stuff for DirWalk */
#if defined __HP__ || __SGI__ || __T3E__ || __PARAGON__ || __DEC__ || __SUN__
        #include <dirent.h>
  typedef struct dirent DIRENT;
        #define D_NAME(d)               ((d)->d_name)

  const char *cb_dir = BasedConvertedFilename(dir);
  DIR *dfd = opendir(cb_dir);
  int ft = filetype(dir);

  PRINTDEBUG(low,1,("DirWalk: dir   = '%s'\n",dir));
  PRINTDEBUG(low,1,("DirWalk: cb_dir= '%s'\n",cb_dir));

  if (ft!=FT_DIR)
    REP_ERR_RETURN (PATH_NO_DIR)
    if (!dfd)
      REP_ERR_RETURN (PATH_INVALID)
      else
      {
        DIRENT *dp;
        while ((dp=readdir(dfd))!=NULL)
        {
          char name[MAX_PATH_LEN];

          if (strcmp(D_NAME(dp),".")==0
              ||
              strcmp(D_NAME(dp),"..")==0)
            /* skip current and up dir */
            continue;

          if (strlen(dir)+strlen(D_NAME(dp))+2 > sizeof(name))
            REP_ERR_RETURN (NAME_TOO_LONG)
            else
            {
              strcpy(name,dir);
              strcat(name,D_NAME(dp));
              (*fcn)(name);
            }

        }
      }
  closedir(dfd);
  return 0;

#elif macintosh

#define pstrlen(s)              ((s)[0])                /* length of pascal string */

  const char *cb_dir = BasedConvertedFilename(dir);
  CInfoPBRec cipbr;                                                     /* local parameter block */
  HFileInfo       *fpb = (HFileInfo *)&cipbr;           /* two pointers */
  DirInfo *dpb = (DirInfo *) &cipbr;
  HVolumeParam hpb;
  OSErr err;
  short idx;
  char dirname[1024];
  char parentdir[1024];
  short volID;
  long dirID;
  int len;

  /* omit leading single ":" (which means "./") */
  if (cb_dir[0]==':' && cb_dir[1]!=':')
    cb_dir++;

  /* copy dir to dirname and delete trailing ':' (but not "::") */
  /* dir is const, thats why we don't modify it directly */
  strcpy(dirname,cb_dir);
  len = strlen(dirname);
  if (len>1 && dirname[len-1]==':' && dirname[len-2]!=':')
    dirname[len-1] = '\0';

  /* first, get volume reference number via PBHGetVInfo... */
  hpb.ioNamePtr = c2pstr((char *)dirname);
  hpb.ioVRefNum = 0;                                            /* use pathname, not ref num */
  hpb.ioVolIndex = -1;                                  /* don't use volume index either */
  err = PBHGetVInfo( (HParmBlkPtr)&hpb, FALSE );
  if (err) {
    if (err==nsvErr)             /* No such volume error */
      REP_ERR_RETURN (VOLUME_NOT_FOUND);
    if (err==paramErr)             /* No such volume error */
      REP_ERR_RETURN (NO_DEFAULT_VOLUME);
  }
  volID = hpb.ioVRefNum;

  /* next, get directory ID via PBGetCatInfo... */
  fpb->ioVRefNum = 0;
  /* partial pathname, and buffer to receive name */
  strcpy(dirname,cb_dir);
  if (len>1 && dirname[len-1]==':' && dirname[len-2]!=':')
    dirname[len-1] = '\0';
  fpb->ioNamePtr = c2pstr((char *)dirname);
  fpb->ioDirID = 0;                                             /* search from working directory */
  fpb->ioFDirIndex = 0;                                 /* gimme info about the named directory */
  err = PBGetCatInfo( &cipbr, FALSE );      /* get the catalog info */
  if (err)        {
    /* No such volume */
    if (err == nsvErr) REP_ERR_RETURN (VOLUME_NOT_FOUND);
    /* I/O error */
    if (err == ioErr) REP_ERR_RETURN (IO_ERROR);
    /* Bad Filename */
    if (err == bdNamErr) REP_ERR_RETURN (BAD_FILENAME);
    /* File not found */
    if (err == fnfErr) REP_ERR_RETURN (FILE_NOT_FOUND);
    /* No default volume */
    if (err == paramErr) REP_ERR_RETURN (NO_DEFAULT_VOLUME);
    /* Directory not found or incomplete pathname */
    if (err == dirNFErr) REP_ERR_RETURN (DIR_NOT_FOUND);
    /* User does not have the correct access */
    if (err == afpAccessDenied) REP_ERR_RETURN (ACCESS_ERROR);
    /* Directory not found or incomplete pathname */
    if (err == afpObjectTypeErr) REP_ERR_RETURN (DIR_NOT_FOUND);
  }
  dirID = dpb->ioDrDirID;               /* get directory ID for subsequent searching */

  /* return if path is not a directory */
  if (!(dpb->ioFlAttrib & ioDirMask)) REP_ERR_RETURN (PATH_NO_DIR)

    /* save parent directory */
    strcpy(parentdir,p2cstr(fpb->ioNamePtr));

  /* now loop through files using PBGetCatInfo... */
  fpb->ioVRefNum = volID;
  for( idx=1; TRUE; idx++) {
    char name[MAX_PATH_LEN];

    fpb->ioDirID = dirID;                       /* set ioDirID on each loop */
    fpb->ioFDirIndex = idx;                     /* index of entry to return */

    err = PBGetCatInfo( &cipbr, FALSE );
    if (err) break;                                     /* exit when no more entries */

    if (strlen(parentdir)+pstrlen(fpb->ioNamePtr)+2 > sizeof(name))
      REP_ERR_RETURN (NAME_TOO_LONG)
      else
      {
        if (parentdir[0]!='\0')
        {
          /* parentdir is not empty */
          strcpy(name,parentdir);
          if (name[strlen(name)-1]!=':')
            strcat(name,":");
          strcat(name,p2cstr(fpb->ioNamePtr));
        }
        else
          strcpy(name,p2cstr(fpb->ioNamePtr));

        /* hrr 980901: subtract BasePath first? */
        (*fcn)(ConvertMachine_2_UNIXPath(name));
      }

  }

  return 0;

#endif

  printf("fileopen.c: DirWalk() not implemented for architecture: %s\n",ARCHNAME);

  REP_ERR_RETURN (NOT_IMPLEMENTED);
}

#ifdef PARALLEL
#define ModelP
#undef PARALLEL
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
  {
    strcpy(thePaths->path[i],Path[i]);
    AppendTrailingSlash(thePaths->path[i]);
  }

  return (0);
}

/****************************************************************************/
/*D
        DirCreateUsingSearchPaths - create a directory searching in the directories specified
                        in the environment item '/Paths/<paths>'

        SYNOPSIS:
        int DirCreateUsingSearchPaths (const char *fname, const char *paths);

        PARAMETERS:
   .   fname - subdirectory name to be created
   .   paths - try paths specified in the environment item '/Paths/<paths> which was
                        set by --> 'ReadSearchingPaths'

        DESCRIPTION:
        The functions trys to create a directory with 'filename' using one by one the
        paths specified in the environment item '/Paths/<paths> which was
        set by --> 'ReadSearchingPaths'. It is used in several places in ug (all paths
        are read from the standard --> 'defaults' file)":"

   .n   'srciptpaths' is used by the interpreter for script execution
   .n   'gridpaths' is used by ugio to read grids from (they are stored in the
   .n   first path

        RETURN VALUE:
        int
   .n   0 sucessfull completion
   .n      != 0 error occured

        SEE ALSO:
        DirCreateUsingSearchPaths_r, mkdir(2)
   D*/
/****************************************************************************/

int DirCreateUsingSearchPaths (const char *fname, const char *paths)
{
  return DirCreateUsingSearchPaths_r ( fname, paths, FALSE);            /* no renaming */
}

/****************************************************************************/
/*D
        DirCreateUsingSearchPaths_r - create a subdirectory searching in the directories specified
                        in the environment item '/Paths/<paths>' incl. renaming option

        SYNOPSIS:
        int DirCreateUsingSearchPaths (const char *fname, const char *paths, int rename);

        PARAMETERS:
   .   fname - subdirectory name to be created
   .   paths - try paths specified in the environment item '/Paths/<paths> which was
                        set by --> 'ReadSearchingPaths'
   .   rename - if TRUE an already existing subdirectory will be renamed

        DESCRIPTION:
        The functions trys to create a subdirectory with 'filename' using one by one the
        paths specified in the environment item '/Paths/<paths> which was
        set by --> 'ReadSearchingPaths'. It is used in several places in ug (all paths
        are read from the standard --> 'defaults' file)":"

   .n   'srciptpaths' is used by the interpreter for script execution
   .n   'gridpaths' is used by ugio to read grids from (they are stored in the
   .n   first path

        RETURN VALUE:
        int
   .n   0 sucessfull completion
   .n      != 0 error occured

        SEE ALSO:
        DirCreateUsingSearchPaths, mkdir(2)
   D*/
/****************************************************************************/

int DirCreateUsingSearchPaths_r (const char *fname, const char *paths, int rename)
{
  PATHS *thePaths;
  FILE *parentDir;

  char fullname[MAXPATHLENGTH];
  INT i,fnamelen,error;
  mode_t mode;

  fnamelen = strlen(fname);
        #ifndef __MACINTOSH__
  mode = S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP;
        #else
  mode = 0;       /* ignored on Macintosh */
        #endif

  PRINTDEBUG(low,1,("DirCreateUsingSearchPaths\n"));
  if (paths == NULL)
  {
    if ((error=mkdir_r(fname,mode,rename))!=0)
      return (1);
    return (0);
  }

  if ((thePaths=GetPaths(paths))==NULL)
    return (1);

  for (i=0; i<thePaths->nPaths; i++)
  {
    /* test whether parent directory exists */
    if ( (parentDir=fopen(thePaths->path[i],"r")) == NULL )
      continue;                         /* this parent directory doesn't exist; try the next one */
    if( (error=fclose(parentDir)) != 0 )
      return (1);

    if (strlen(thePaths->path[i])+fnamelen>MAXPATHLENGTH)
      return (1);

    strcpy(fullname,thePaths->path[i]);
    strcat(fullname,fname);

    if ((error=mkdir_r(fullname,mode,rename))!=0)
      return (1);
    return (0);                 /* subdirectory created sucessfully */
  }
  return (1);
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
        FileOpenUsingSearchPaths_r, FileOpenUsingSearchPath, FileOpenUsingSearchPath_r, GetPaths, ReadSearchingPaths, fileopen
   D*/
/****************************************************************************/

FILE *FileOpenUsingSearchPaths (const char *fname, const char *mode, const char *paths)
{
  return FileOpenUsingSearchPaths_r( fname, mode, paths, FALSE );       /* no renaming */
}

/****************************************************************************/
/*D
        FileOpenUsingSearchPaths_r - open file searching in the directories specified
                        in the environment item '/Paths/<paths>' incl. renaming option

        SYNOPSIS:
        FILE *FileOpenUsingSearchPaths_r (const char *fname, const char *mode, const char *paths, int rename)

    PARAMETERS:
   .   fname - file name to be opened
   .   mode - see ANSI-C 'fopen'
   .   paths - try paths specified in the environment item '/Paths/<paths> which was
                        set by --> 'ReadSearchingPaths'
   .   rename - if TRUE an already existing file will be renamed (wise only for writing)

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
        FileOpenUsingSearchPaths, FileOpenUsingSearchPath, FileOpenUsingSearchPath_r, GetPaths, ReadSearchingPaths, fileopen
   D*/
/****************************************************************************/

FILE *FileOpenUsingSearchPaths_r (const char *fname, const char *mode, const char *paths, int rename)
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
    strcat(fullname,fname);

    if ((theFile=fileopen_r(fullname,mode,rename))!=NULL)
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
        FileOpenUsingSearchPath_r, FileOpenUsingSearchPaths, FileOpenUsingSearchPaths_r, fileopen
   D*/
/****************************************************************************/

FILE *FileOpenUsingSearchPath (const char *fname, const char *mode, const char *path)
{
  return FileOpenUsingSearchPath_r( fname, mode, path, FALSE );         /* no renaming */
}

/****************************************************************************/
/*D
        FileOpenUsingSearchPath_r - try to open a file in the specified path and renaming option

        SYNOPSIS:
        FILE *FileOpenUsingSearchPath_r (const char *fname, const char *mode, const char *path, int rename)

    PARAMETERS:
   .   fname - open file with this name
   .   mode - see ANSI-C 'fopen'
   .   path - path to which fname is to be appended
   .   rename - if TRUE an already existing file will be renamed (wise only for writing)

        DESCRIPTION:
        Try to open a file in the specified path.

        RETURN VALUE:
        FILE *
   .n   pointer to file opened, 'NULL' if error

        SEE ALSO:
        FileOpenUsingSearchPath, FileOpenUsingSearchPaths, FileOpenUsingSearchPaths_r, fileopen
   D*/
/****************************************************************************/

FILE *FileOpenUsingSearchPath_r (const char *fname, const char *mode, const char *path, int rename)
{
  FILE *theFile;
  char fullname[MAXPATHLENGTH];

  if (strlen(path)+strlen(fname)>MAXPATHLENGTH)
    return (NULL);

  strcpy(fullname,path);
  strcat(fullname,fname);

  if ((theFile=fileopen_r(fullname,mode,rename))!=NULL)
    return (theFile);

  return (NULL);
}

/****************************************************************************/
/*D
        FileTypeUsingSearchPaths - give type of file searching in the
            directories specified in the environment item '/Paths/<paths>'

        SYNOPSIS:
        int FileTypeUsingSearchPaths (const char *fname, const char *paths)

    PARAMETERS:
   .   fname - file name to be opened
   .   paths - try paths specified in the environment item '/Paths/<paths> which was
                        set by --> 'ReadSearchingPaths'

        DESCRIPTION:
        The functions trys to determine the file type of the file named
        'filename' using one by one the paths specified in the environment
        item '/Paths/<paths> which was set by --> 'ReadSearchingPaths'.
        It is used in several places in ug (all paths are read from the
        standard --> 'defaults' file)":"

   .n   'srciptpaths' is used by the interpreter for script execution
   .n   'gridpaths' is used by ugio to read grids from (they are stored in the
   .n   first path)

        RETURN VALUE:
        int
   .n      file type (one of FT_UNKNOWN, FT_FILE, FT_DIR, FT_LINK)

        SEE ALSO:
        ReadSearchingPaths, filetype
   D*/
/****************************************************************************/

int FileTypeUsingSearchPaths (const char *fname, const char *paths)
{
  PATHS *thePaths;
  int ftype;
  char fullname[MAXPATHLENGTH];
  INT i,fnamelen;

  fnamelen = strlen(fname);

  if ((thePaths=GetPaths(paths))==NULL)
    return (FT_UNKNOWN);

  for (i=0; i<thePaths->nPaths; i++)
  {
    if (strlen(thePaths->path[i])+fnamelen>MAXPATHLENGTH)
      return (FT_UNKNOWN);

    strcpy(fullname,thePaths->path[i]);
    strcat(fullname,fname);

    if ((ftype=filetype(fullname))!=FT_UNKNOWN)
      return (ftype);
  }

  return (FT_UNKNOWN);
}

int AppendTrailingSlash (char *path)
{
  if (path[0]!='\0' && path[strlen(path)-1]!='/')
  {
    strcat(path,"/");
    return YES;
  }
  return NO;
}

char *SimplifyPath (char *path)
{
  const char *pf;
  char       *pt;

  PRINTDEBUG(low,2,("SimplifyPath: original path= '%s'\n",path));

  /* cancel ./ (not first one) */
  pf = pt = strchr(path,'/');
  if (pf!=NULL)
  {
    while (*pf)
    {
      if (pf[0]=='.' && pf[1]=='/')
        if (*(pf-1)=='/')
        {
          /* eat ./ */
          pf += 2;
          continue;
        }
      if (pt!=pf)
        *pt = *pf;
      pf++;
      pt++;
    }
    *pt = '\0';
  }

  PRINTDEBUG(low,2,("SimplifyPath: path= '%s'\n",path));

  /* cancel ../ where possible */
  pf = pt = path;
  for (; *pf; pf++,pt++)
  {
    if (pf[0]=='.' && pf[1]=='.' && pf[2]=='/')
      if (pf==path || *(pf-1)=='/')
      {
        char *pd = pt-1;

        while (pd>path)
          if (*(--pd)=='/')
            break;
        if (*pd=='/' && !(pd[0]=='/' && pd[1]=='.' && pd[2]=='.' && pd[3]=='/'))
        {
          /* eat ../ and reset pt */
          pf += 2;
          pt = pd;
          continue;
        }
      }
    *pt = *pf;
  }
  *pt = '\0';

  return path;
}

const char *GetBasePath (void)
{
  return BasePath;
}

const char *SetBasePath (const char *path)
{
  strcpy(OldBasePath,path);
  strcpy(BasePath,path);
  AppendTrailingSlash(BasePath);

  PRINTDEBUG(low,1,("SetBasePath: BasePath= '%s'\n",BasePath));

  return OldBasePath;
}

const char *AddBasePath (const char *path)
{
  strcpy(OldBasePath,path);
  strcat(BasePath,path);
  AppendTrailingSlash(BasePath);

  SimplifyPath(BasePath);

  PRINTDEBUG(low,1,("AddBasePath: BasePath= '%s'\n",BasePath));

  return OldBasePath;
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
