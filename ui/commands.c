// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  commands.c													*/
/*																			*/
/* Purpose:   definition of all dimension independent commands of ug		*/
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*			  02.02.95 begin, ug version 3.0								*/
/*																			*/
/* Remarks:   for dimension dependent commands see commands2d/3d.c			*/
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment cmd
#endif

/****************************************************************************/
/*																			*/
/*		defines to exclude functions										*/
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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/* low module */
#include "compiler.h"
#include "defaults.h"
#include "misc.h"
#include "ugstruct.h"
#include "fileopen.h"
#include "ugenv.h"

/* devices module */
#include "devices.h"

/* grid generatormodule */
#include "ggm.h"
#include "ggmain.h"

/* grid manager module */
#include "gm.h"
#include "evm.h"
#include "formats.h"

/* numerics module */
#include "num.h"

/* graph module */
#include "wpm.h"
#include "wop.h"

/* user interface module */
#include "uginterface.h"
#include "ugstruct.h"
#include "cmdint.h"
#include "cmdline.h"
#include "helpmsg.h"

/* own header */
#include "commands.h"

/* includes for 3D application */
#ifdef __THREEDIM__
#include "ugm3d.h"
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

#define BUFFERSIZE                              512     /* size of the general purpose text buff*/

#define WHITESPACE                              " \t"

#define LONGSTRSIZE                     256 /* size of some strings                             */
#define LONGSTRLEN                              255 /* len of some strings                                      */
#define LONGSTRLENSTR                   "255" /* LONGSTRLEN as string                           */

/* for ProtoOnCommand */
#define NORENAME_PROTO                  0
#define APPEND_PROTO                    1
#define RENAME_PROTO                    2
#define TRYRENAME_PROTO                 3
#define MAXPATHLENGTH                   255
#define MAXRENAMECHAR                   'z'

/* for the .list commands */
#define DO_ID                                   1
#define DO_SELECTION                    2
#define DO_ALL                                  3

/* for MarkCommand */
#define MARK_ALL                                1
#define MARK_COARSEN                    2
#define MARK_ID                                 3
#define MARK_SELECTION                  4
#define NO_SIDE_SPECIFIED               -1
#define NO_RULE_SPECIFIED               -1
#define NO_OF_RULES                     10

/* for save command */
#define NO_COMMENT                               "no comment"

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct MarkRule
{
  char *RuleName;                                       /* what you type in the mark cmdline	*/
  INT RuleId;                                           /* corresponding rule ID for refine     */
};

typedef struct MarkRule MARKRULE;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static MULTIGRID *currMG=NULL;                  /* the current multigrid			*/

static NUM_PROC *currNumProc=NULL;              /* current numerical procedure		*/
static char buffer[BUFFERSIZE];         /* general purpose text buffer		*/

static FILE     *protocolFile=NULL;     /* for protocol commands			*/

static MARKRULE myMR[NO_OF_RULES]=      /* name and ID of available rules	*/
{{"red",        RED},
 {"no",         NO_REFINEMENT},
 {"blue",       BLUE},
 {"copy",       COPY},
 {"bi_1",       BISECTION_1},
 {"bi_2q",      BISECTION_2_Q},
 {"bi_2t1", BISECTION_2_T1},
 {"bi_2t2", BISECTION_2_T2},
 {"bi_3",       BISECTION_3},
 {"coarse", UNREFINE}};

static clock_t Time0;                                   /* time offset for readclock		*/

static char userPath[1024];                     /* environment path for ls,cd		*/

static INT untitledCounter=0;                   /* counter for untitled multigrids	*/

/* some variables to transfer info between QualityCommand and QualityElement*/
static DOUBLE min,max,themin,themax,minangle,maxangle;
static INT lessopt,greateropt;
static char mintext[32],maxtext[32],minmaxtext[32];

/* counters for windows and pictures */
static INT wincounter=1;
static INT piccounter=1;

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*D
   GetCurrentMultigrid - return a pointer to the current multigrid

   SYNOPSIS:
   MULTIGRID *GetCurrentMultigrid ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function returns a pionter to the current multigrid.

   RETURN VALUE:
   MULTIGRID *
   .n     pointer to multigrid
   .n     NULL if there is no current multigrid.
   D*/
/****************************************************************************/

MULTIGRID *GetCurrentMultigrid ()
{
  return (currMG);
}

/****************************************************************************/
/*D
   SetCurrentMultigrid - set the current multigrid if it is valid

   SYNOPSIS:
   INT SetCurrentMultigrid (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - pointer to multigrid

   DESCRIPTION:
   This function sets the current multigrid if it is valid, i. e.
   the function checks whether 'theMG' acually points to a multigrid.
   It can be NULL only if no multigrid is open.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if theMG is not in the multigrid list
   D*/
/****************************************************************************/

INT SetCurrentMultigrid (MULTIGRID *theMG)
{
  MULTIGRID *mg;

  mg = GetFirstMultigrid();
  if (mg==theMG)
  {
    /* possibly NULL */
    currMG = theMG;
    return (0);
  }

  for (; mg!=NULL; mg=GetNextMultigrid(mg))
    if (mg==theMG)
    {
      /* never NULL */
      currMG = theMG;
      return (0);
    }

  return (1);
}

/****************************************************************************/
/*D
   GetCurrentNumProc - return a pointer to the current numproc

   SYNOPSIS:
   static NUM_PROC *GetCurrentNumProc ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function returns a pointer to the current numproc.

   RETURN VALUE:
   NUM_PROC *
   .n      pointer to numproc
   .n      NULL if there is no current numproc.
   D*/
/****************************************************************************/

static NUM_PROC *GetCurrentNumProc ()
{
  return (currNumProc);
}

/****************************************************************************/
/*D
   SetCurrentNumProc -	Set the current NumProc if it is valid

   SYNOPSIS:
   static INT SetCurrentNumProc (NUM_PROC *theNumProc);

   PARAMETERS:
   .  theNumProc - pointer to multigrid

   DESCRIPTION:
   This function sets the current NumProc if it is valid, i. e.
   the function checks whether 'theNumProc' acually points to a numproc.
   It can be NULL only if no numproc is defined.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if theNumProc is not in the numproc list
   D*/
/****************************************************************************/

static INT SetCurrentNumProc (NUM_PROC *theNumProc)
{
  NUM_PROC *np;

  np = GetFirstNumProc();
  if (np==theNumProc)
  {
    /* possibly NULL */
    currNumProc = theNumProc;
    return (0);
  }

  for (; np!=NULL; np=GetNextNumProc(np))
    if (np==theNumProc)
    {
      /* never NULL */
      currNumProc = theNumProc;
      return (0);
    }

  return (1);
}

/****************************************************************************/
/*D
   quit - quit command

   DESCRIPTION:
   This command quits the program and closes the shell.

   'quit'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   QuitCommand - Quit programm

   SYNOPSIS:
   static INT QuitCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function quits programm (checks open files, documents, ... ).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT QuitCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  return(QUITCODE);
}

/****************************************************************************/
/*D
   mute - set mutelevel

   DESCRIPTION:
   This command sets a mutelevel.
   The default value is 0 and all skript lines will be printed on the shell.
   This will be suppressed by mutelevel -1.
   Smaller muteleveles should reduce the output further.

   'mute <value>'
   .   <value> - integer which gives the mutelevel

   REMARK:
   Formally, this is not an ug command, 'mute' is checked in
   'InterpretString'.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   help - search for for a command and prints the help

   DESCRIPTION:
   This command searches for a command and prints the help.
   It gets online help for commands etc.

   help [[<helpitem>] $k]

   .   no~option             - this is  equivalent to 'help help'
   .   <helpitem>            - print help for <helpitem> (string)
   .   $k                    - search for keyword <helpitem>

   EXAMPLE:
   'help plot $k'

   prints a list of all commands which are relevant for plotting
   (openwindow, setview, zoom ...)

   D*/
/****************************************************************************/

/****************************************************************************/
/*
   HelpCommand - Search for for a command and prints the help

   SYNOPSIS:
   static INT HelpCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function searches for a command and prints the help.
   It gets online help for commands etc.

   help [[<helpitem>] $k]

   .   no option             - this is  equivalent to 'help help'
   .   <helpitem>            - print help for <helpitem> (string)
   .   $k                    - search for keyword <helpitem>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT HelpCommand (INT argc, char **argv)
{
  INT i,res,mode,rv;
  COMMAND *Cmd;
  char buf[NAMESIZE];

  mode = HELPITEM;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'k' :
      mode = KEYWORD;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("help",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* get HelpFor string */
  res = sscanf(argv[0],expandfmt(CONCAT3("help %",NAMELENSTR,"[0-9a-zA-Z_]")),buf);
  if (res==1)
  {
    rv = PrintHelp(buf,mode,NULL);
    if (rv!=HELP_OK)
    {
      /* another try: ug command buf existing? */
      Cmd = SearchUgCmd(buf);
      if (Cmd!=NULL)
        rv = PrintHelp(ENVITEM_NAME(Cmd),mode,NULL);
    }
  }
  else
    rv = PrintHelp("help",mode,NULL);

  switch (rv)
  {
  case HELP_OK :
    return (OKCODE);

  case HELP_NOT_FOUND :
    sprintf(buffer," no help entry found for '%s'\n",buf);
    UserWrite(buffer);
    return (OKCODE);

  default :
    PrintErrorMessage('E',"help","(unknown)");
  }

  return (CMDERRORCODE);
}

/****************************************************************************/
/*D
   checkhelp - Check wether all commands in /menu have a help item

   DESCRIPTION:
   This function checks wether for all commands in /menu a help item exists.
   It also checks wether for all num proc types a help item exists.

   It prints all commands and num proc types for which help does NOT exist.

   It calls the funtion 'CheckHelp'.

   'checkhelp'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   CheckHelpCommand - Check wether all commands in /menu have a help item

   SYNOPSIS:
   static INT CheckHelpCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function checks wether for all commands in /menu a help item exists.
   It also checks wether for all num proc types a help item exists.

   It prints all commands and num proc types for which help does NOT exist.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT CheckHelpCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  CheckHelp();

  return (OKCODE);
}

/****************************************************************************/
/*D
   cmfn - creates a metafile name

   DESCRIPTION:
   This command creates a metafile names.
   It creates a string containing the name of a metafile name for
   animation by xugv.

   'cmfn <name> <var>'

   .  <name> - first part of the metafile names
   .  <var> - the contents of var will be appended to the name

   EXAMPLE:
   .vb
   frame="
    cmfn film step;
    openwindow 0 0 820 420 $d meta $n @film;
    openpicture $s 10 10 800 400 $n framepic;
    setplotobject EScalar $e S2 $m COLOR $d 0 $f 0.0 $t 1.0;
    setview;
    zoom 0.4;
    plot;
    closewindow;
   ";
   step = 0;
   steps = 100;
   @frame;
   repeat {
    print "STEP ", step;
    @mysolve;
    step=step+1;
    @frame;
    if (step==steps) break;
   }
   .ve

   This runs 'mysolve' 100 times and
   creates metafiles 'film.0000', 'film.0001', 'film.0002', ... 'film.0100'.
   D*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  CreateMetafileNameCommand										*/
/*																			*/
/* Purpose:   create a string containing the name of a metafile name for	*/
/*			  animation by xugv												*/
/*																			*/
/* Input:	  INT argc: number of arguments (incl. its own name                     */
/*			  char **argv: array of strings giving the arguments			*/
/*																			*/
/* Output:	  INT 0: everything ok											*/
/*			  INT >0: error                                                                                                 */
/*																			*/
/****************************************************************************/

static INT CreateMetafileNameCommand (INT argc, char **argv)
{
  INT res,frame;
  char name[LONGSTRSIZE];
  char fullname[LONGSTRSIZE];

  res = sscanf(argv[0],expandfmt(CONCAT3(" cmfn %",LONGSTRLENSTR,"[0-9:.a-zA-Z_] %255[0-9:.a-zA-Z_]")),name,buffer);
  if (res!=2) return(CMDERRORCODE);

  if (GetStringValueInt(buffer,&frame)) return(CMDERRORCODE);

  sprintf(fullname,"%s.%04d",name,frame);

  if (SetStringVar(name,fullname)) return(CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   readclock - printing execution time

   DESCRIPTION:
   This command is for measuring the time used.
   It prints the execution time since the last 'resetclock' to
   string variable ':CLOCK'.

   'readclock'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ReadClockCommand - For measuring the time used

   SYNOPSIS:
   static INT ReadClockCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function is for measuring the time used.
   It reads all time clock.

   .  readclock       - prints the execution time since the last 'resetclock' to
                     string variable ':CLOCK'

   `Warning:` be careful measuring long times because depending on the
   implementation the wrap around time may be short.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ReadClockCommand (INT argc, char **argv)
{
  DOUBLE Time;

  NO_OPTION_CHECK(argc,argv);

  Time = (clock()-Time0)/((DOUBLE) CLOCKS_PER_SEC);

  if (SetStringValue(":CLOCK",Time)!=0)
  {
    PrintErrorMessage('E',"readclock","could not get string variable :CLOCK");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   resetclock - starting the time mesuring

   DESCRIPTION:
   This command starts the time mesuring.
   It sets the global variable 'Time0' to zero.

   'resetclock'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ResetClockCommand - For measuring the time used

   SYNOPSIS:
   static INT ResetClockCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function is for measuring the time used.
   It resets all time clock.

   .  resetclock      - reset the execution time clock to zero

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ResetClockCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  Time0 = clock();

  return (OKCODE);
}

/****************************************************************************/
/*D
   InitClock() - starting the time mesuring

   SYNOPSIS:
   static INT InitClock();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function starts the time mesuring.
   It sets the global variable 'Time0' to zero.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT InitClock()
{

  Time0 = clock();

  return(0);
}

/****************************************************************************/
/*D
   date - prints the date

   DESCRIPTION:
   This command prints the date on the shell resp.
   writes it in the string variable ':date'.

   'date [$s] [$S]'

   .  no~option - print the date on the shell
   .  $s   -  put in the string variable ':date'.
   .  $S   -  use short format of the form yy.mm.dd

   D*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  DateCommand													*/
/*																			*/
/* Purpose:   prints date and time					                                                */
/*																			*/
/* Input:	  INT argc: number of arguments (incl. its own name                     */
/*			  char **argv: array of strings giving the arguments			*/
/*																			*/
/* Output:	  INT return code see header file								*/
/*																			*/
/****************************************************************************/

static INT DateCommand (INT argc, char **argv)
{
  time_t Time;
  char *fmt;
  INT i,svopt;

  /* check options */
  svopt = FALSE;
  fmt = "%a %b %d %H:%M:%S %Y";
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      svopt = TRUE;
      break;

    case 'S' :
      fmt = "%y.%m.%d";;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("date",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  time(&Time);
  strftime(buffer,BUFFERSIZE,fmt,localtime(&Time));

  if (svopt)
    SetStringVar(":date",buffer);
  else
    UserWriteF("%s\n",buffer);

  return (OKCODE);
}

/****************************************************************************/
/*D
   ls - lists the content of an environment directory.

   DESCRIPTION:
   This command lists the content of an environment directory.

   'ls [<path>]'

   .  no~option - lists the content of the current directory.
   .  <path>  - contains the relative or absolute path in UNIX-style
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ListEnvCommand - Navigate through the environment tree

   SYNOPSIS:
   static INT ListEnvCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function navigates through the environment tree.
   It lists environment directory.
   It uses the function 'ChangeEnvDir'.

   ls [<path>]

   .  <path>                  - <path> contains the relative or absolute path in UNIX-style
   .n                           if <path> is omitted: ls <the current directory>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ListEnvCommand (INT argc, char **argv)
{
  ENVDIR *currentDir;
  ENVITEM *theItem;
  char *s;
  INT i;

  NO_OPTION_CHECK(argc,argv);

  /* load current directory */
  currentDir = ChangeEnvDir(userPath);
  if (currentDir==NULL)
  {
    /* directory is invalid -> change to root directory */
    strcpy(userPath,DIRSEP);
    currentDir = ChangeEnvDir(userPath);
    if (currentDir == NULL)
      return (CMDERRORCODE);
  }

  /* strip ' '*ls' '* */
  s = strchr(argv[0],'l');
  strcpy(buffer,s);
  i = 2;
  while ((buffer[i]!='\0') && (strchr(WHITESPACE,buffer[i])!=NULL)) i++;
  s = buffer+i;

  /* pathname is now in s: change dir to path */
  if (strlen(s)>0)
    currentDir = ChangeEnvDir(s);

  if (currentDir==NULL)
  {
    PrintErrorMessage('E',"ls","invalid path as argument");
    return (CMDERRORCODE);
  }

  theItem = currentDir->down;
  while (theItem!=NULL)
  {
    UserWrite(theItem->v.name);
    if (theItem->v.type%2==0)
      UserWrite("\n");
    else
      UserWrite("*\n");
    theItem = theItem->v.next;
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   cd - change the current environment directory.

   DESCRIPTION:
   This command changes the current environment directory.
   It uses the function 'ChangeEnvDir'.

   'cd [<path>]'

   .  <path> - <path> contains the relative or absolute path in UNIX-style
   .  no~option - cd to root (cd /)
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ChangeEnvCommand - Navigate through the environment tree

   SYNOPSIS:
   static INT ChangeEnvCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function navigates through the environment tree.
   It changes environment directory to specified path.

   cd [<path>]

   .  <path>                  - <path> contains the relative or absolute path in UNIX-style
   .                          - if <path> is omitted: cd to root (cd /)
 */
/****************************************************************************/

static INT ChangeEnvCommand (INT argc, char **argv)
{
  ENVDIR *currentDir;
  char *s;
  INT i;

  NO_OPTION_CHECK(argc,argv);

  /* load current directory */
  currentDir = ChangeEnvDir(userPath);
  if (currentDir==NULL)
  {
    /* directory is invalid -> change to root directory */
    strcpy(userPath,DIRSEP);
    currentDir = ChangeEnvDir(userPath);
    if (currentDir == NULL)
      return (CMDERRORCODE);
  }

  /* strip ' '*cd' '* */
  s = strchr(argv[0],'c');
  strcpy(buffer,s);
  i = 2;
  while ((buffer[i]!='\0') && (strchr(WHITESPACE,buffer[i])!=NULL)) i++;
  s = buffer+i;

  /* pathname is now in buffer */
  currentDir = ChangeEnvDir(s);
  if (currentDir==NULL)
  {
    PrintErrorMessage('E',"cd","invalid path as argument");
    return (CMDERRORCODE);
  }
  GetPathName(userPath);
  UserWrite(userPath);
  UserWrite("\n");

  return (OKCODE);
}

/****************************************************************************/
/*D
   pwd - print the current environment

   DESCRIPTION:
   This command print the current environment on the shell.
   It uses the function 'CangeEnvDir'.

   'pwd'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   PrintEnvDirCommand - Navigate through the environment tree

   SYNOPSIS:
   static INT PrintEnvDirCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function navigates through the environment tree.
   It prints environment working directory.

   pwd

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PrintEnvDirCommand (INT argc, char **argv)
{
  ENVDIR *currentDir;

  NO_OPTION_CHECK(argc,argv);

  /* load current directory */
  currentDir = ChangeEnvDir(userPath);
  if (currentDir==NULL)
  {
    /* directory is invalid */
    strcpy(userPath,DIRSEP);
    currentDir = ChangeEnvDir(userPath);
    if (currentDir == NULL)
      return (CMDERRORCODE);
  }

  GetPathName(userPath);
  UserWrite(userPath);
  UserWrite("\n");

  return(OKCODE);
}

/****************************************************************************/
/*D
   envinfo - print total size and used memory

   DESCRIPTION:
   This command prints total size and used memory on the shell.

   'envinfo'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   EnvInfoCommand - Print total size and used memory

   SYNOPSIS:
   static INT EnvInfoCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function prints used and size of environment heap.
   It calls the function 'EnvHeapInfo'.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EnvInfoCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  EnvHeapInfo(buffer);

  UserWrite(buffer);

  return (OKCODE);
}

/****************************************************************************/
/*
   SetCommand - Set (or print) a string variable struct

   SYNOPSIS:
   static INT SetCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function sets (or prints) a string variable struct.
   It sets or prints the contents of a struct or struct directory
   using the functions 'PrintStructContents', 'SetStringVar',
   'PrintCurrentStructContents'.

   set {<struct> <value>} | {[<structdir> | <struct>] [$r]}

   .  <struct> <value>              - assign <value> (just a string of arbitraray length) to <struct>
   .  [<structdir> | <struct>] [$r] - display contents of <struct> or <structdir> .n                               (default: current struct dir)

   . $r                             - specifies the directory, its contents is listed recursively
   .n                               (see also: help structpath)


   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   set - set (or print) a string variable struct

   DESCRIPTION:
   This command sets (or prints) a string variable struct.
   It sets or prints the contents of a struct or struct directory.

   'set {<struct> <value>} | {[<structdir> | <struct>] [$r]}'

   .  <struct>~<value>   - assign <value> (just a string of arbitraray length) to <struct>
   .  [<structdir>|<struct>]~[$r] - display contents of <struct> or <structdir>

   .n                               (default: current struct dir)

   . $r        - specifies the directory, its contents is listed recursively


   SEE ALSO:
   'structpath'
   D*/
/****************************************************************************/

static INT SetCommand (INT argc, char **argv)
{
  char name[LONGSTRSIZE], *namePtr;
  INT flag,ropt,i,res,rv;

  res = sscanf(argv[0],expandfmt(CONCAT3(" set %",LONGSTRLENSTR,"[0-9:.a-zA-Z_] %255[ -~]")),name,buffer);

  /* check options */
  ropt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'r' :
      if (res>1)
      {
        PrintHelp("set",HELPITEM," (the r option applies not with setting a value)");
        return (PARAMERRORCODE);
      }
      ropt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("set",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* print or set struct */
  switch (res)
  {
  case 1 :
    namePtr=&(name[0]);
    do
    {
      rv=PrintStructContents(namePtr,buffer,BUFFERSIZE,ropt);
      if ((rv!=0)&&(rv!=4))
      {
        PrintErrorMessage('E',"set","structure not found or bad structure");
        return (CMDERRORCODE);
      }
      UserWrite(buffer);
      namePtr=NULL;
    }
    while (rv==4);
    break;

  case 2 :
    rv =  SetStringVar(name,buffer);
    if (rv!=0)
    {
      PrintErrorMessage('E',"set","could not allocate variable");
      return (CMDERRORCODE);
    }

    break;

  default :
    flag=1;
    do
    {
      rv=PrintCurrentStructContents(flag,buffer,BUFFERSIZE,ropt);
      if ((rv!=0)&&(rv!=4))
      {
        PrintErrorMessage('E',"set","structure not found or bad structure");
        return (CMDERRORCODE);
      }
      UserWrite(buffer);
      flag=0;
    }
    while (rv==4);
    break;
  }

  if (rv==0)
    return (OKCODE);
  else
    return (CMDERRORCODE);
}

/****************************************************************************/
/*
   DeleteVariableCommand - Delete an existing variable

   SYNOPSIS:
   static INT DeleteVariableCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function deletes an existing variable.

   dv <variable name>

   .  <variable name>         - <variable name> consists of a complete path related to the
   .n                         current struct dir or the structure root directory in the environment
   .n                         (see also: help structpath)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT DeleteVariableCommand (INT argc, char **argv)
{
  INT res;
  char name[LONGSTRSIZE];

  NO_OPTION_CHECK(argc,argv);

  res = sscanf(argv[0],expandfmt(CONCAT3(" dv %",LONGSTRLENSTR,"[0-9:.a-zA-Z_]")),name);

  if (res!=1)
  {
    PrintHelp("dv",HELPITEM," (could not read name of variable)");
    return(PARAMERRORCODE);
  }

  if (argc!=1)
  {
    PrintHelp("dv",HELPITEM,NULL);
    return(PARAMERRORCODE);
  }

  if (DeleteVariable(name)!=0)
  {
    PrintErrorMessage('E',"dv","could not delete variable");
    return (CMDERRORCODE);
  }
  else
    return(DONE);
}

/****************************************************************************/
/*D
   ms  - create a structure

   DESCRIPTION:
   This commands creates a new string variable struct.
   It calls the function 'MakeStruct'.

   'ms <structdir>'

   .  <structdir>             - the <structdir> consists of a complete path related to the
                             current struct dir or the string variable root in the environment

   SEE ALSO:
   'structpath'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   MakeStructCommand - Create a structure

   SYNOPSIS:
   static INT MakeStructCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function creates a structure.
   It creates a new string variable struct.

   ms <structdir>

   .  <structdir>             - the <structdir> consists of a complete path related to the
                             current struct dir or the string variable root in the environment
                             (see also: help structpath)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MakeStructCommand (INT argc, char **argv)
{
  INT res;
  char name[LONGSTRSIZE];

  NO_OPTION_CHECK(argc,argv);

  res = sscanf(argv[0],expandfmt(CONCAT3(" ms %",LONGSTRLENSTR,"[0-9:.a-zA-Z_]")),name);

  if (res!=1)
  {
    PrintHelp("ms",HELPITEM," (could not read name of struct)");
    return(PARAMERRORCODE);
  }

  if (MakeStruct(name)!=0)
    return (CMDERRORCODE);
  else
    return (OKCODE);
}

/****************************************************************************/
/*D
   cs  - change to a struct

   DESCRIPTION:
   This commands changes to a struct.
   It calls the function 'ChangeStructDir'.

   'cs <structdir>'

   .  <structdir>             - the <structdir> consists of a complete path related to the
                             current struct dir or the string variable root in the environment
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ChangeStructCommand - Create a structure

   SYNOPSIS:
   static INT ChangeStructCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function changes a structure.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ChangeStructCommand (INT argc, char **argv)
{
  char *s;
  INT i;

  NO_OPTION_CHECK(argc,argv);

  /* strip ' '*cs' '* */
  s = strchr(argv[0],'c');
  strcpy(buffer,s);
  i = 2;
  while ((buffer[i]!='\0') && (strchr(WHITESPACE,buffer[i])!=NULL)) i++;
  s = buffer+i;

  /* s is now the pathname */
  if (ChangeStructDir(s)==NULL)
  {
    PrintErrorMessage('E',"cs","invalid path as argument");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   pws  - prints the current struct path

   DESCRIPTION:
   This commands calls the function 'GetStructPathName' and
   puts the result on the screen.

   'pws'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   PrintWorkStructCommand - prints the path of a struct

   SYNOPSIS:
   static INT PrintWorkStructCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function prints the path of a struct.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PrintWorkStructCommand (INT argc, char **argv)
{
  char structPath[1024];

  NO_OPTION_CHECK(argc,argv);

  GetStructPathName(structPath, 1024);
  UserWrite(structPath);
  UserWrite("\n");

  return(OKCODE);
}

/****************************************************************************/
/*D
   ds  - deletes a struct

   DESCRIPTION:
   This commands calls the function 'DeleteStruct' to remove a struct.

   'ds <structdir>'

   .  <structdir>             - the <structdir> consists of a complete path related to the
   .n                         current struct dir or the string variable root in the environment
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   DeleteStructCommand - Delete an existing structure

   SYNOPSIS:
   static INT DeleteStructCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function deletes an existing structure.
   It deletes an existing struct directory.

   ds <structdir>

   .  <structdir>             - the <structdir> consists of a complete path related to the
   .n                         current struct dir or the string variable root in the environment
   .n                         (see also: help structpath)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT DeleteStructCommand (INT argc, char **argv)
{
  INT res;
  char name[LONGSTRSIZE];

  NO_OPTION_CHECK(argc,argv);

  res = sscanf(argv[0],expandfmt(CONCAT3(" ds %",LONGSTRLENSTR,"[0-9:.a-zA-Z_]")),name);

  if (res!=1)
  {
    PrintHelp("ds",HELPITEM," (could not read name of struct)");
    return(PARAMERRORCODE);
  }

  if (argc!=1)
  {
    PrintHelp("ds",HELPITEM,NULL);
    return(PARAMERRORCODE);
  }

  if (DeleteStruct(name)!=0)
  {
    PrintErrorMessage('E',"ds","could not delete structure");
    return (CMDERRORCODE);
  }
  else
    return(DONE);
}

/****************************************************************************/
/*D
   protocol - print strings to the protocol file

   DESCRIPTION:
   This command prints strings to protocol file.
   It writes formatted output to the open protocol file.

   'protocol {$i[ ]<verbatim text> | $n[ ]<verbatim text> | $t[ ]<verbatim text> | $f}*'

   .vb
    $\i   append <verbatim text> to protocol file
    $\n   write a line feed and append <verbatim text> to protocol file
    $\t   write a tab and append <verbatim text> to protocol file
          NOTE: the first space (if there) following the option character is skipped

    $\f   flush the file buffer
   .ve

   EXAMPLE:
   .vb
   x = exp(1);
   protoOn exp.proto;
   protocol $\i the value of exp(1) is $\t @x;
   protocol $\n you can use $s in protocol;
   protoOff
   .ve

   Then, the file 'exp.proto' will consists of the string
   .vb
   "the value of exp(1) is\t2.7182818\nyou can use $s in protocol"
   .ve

   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ProtocolCommand - Print strings to protocol file

   SYNOPSIS:
   static INT ProtocolCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function prints strings to protocol file.
   It writes formatted ouput to the open protocol file.

   protocol {$i[ ]<verbatim text> | $n[ ]<verbatim text> | $t[ ]<verbatim text> | $f}*

   .   $\i                     - append <verbatim text> to protocol file
   .   $\n                     - write a line feed and append <verbatim text> to protocol file
   .   $\t                     - write a tab and append <verbatim text> to protocol file
   .n                         NOTE: the first space (if there) following the option character is skipped

   .   $\f                     - flush the file buffer

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ProtocolCommand (INT argc, char **argv)
{
  INT i,from;

  if (protocolFile==NULL)
  {
    PrintErrorMessage('E',"protocol","no protocol file open!");
    return (CMDERRORCODE);
  }

  for (i=1; i<argc; i++)
  {
    if (argv[i][0]!='\\')
    {
      PrintErrorMessage('E',"protocol","protocol options have to begin with %");
      return (PARAMERRORCODE);
    }
    from = (argv[i][2]==' ') ? 3 : 2;
    switch (argv[i][1])
    {
    case 'i' :
      fprintf(protocolFile,"%s",(argv[i])+from);
      break;

    case 't' :
      fprintf(protocolFile,"\t%s",(argv[i])+from);
      break;

    case 'n' :
      fprintf(protocolFile,"\n%s",(argv[i])+from);
      break;

    case 'f' :
      fflush(protocolFile);
      continue;

    default :
      sprintf(buffer," (unknown option '%s')",argv[i]);
      PrintHelp("protocol",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }
    /* write options not followed by a \ */
    while ((i+1<argc) && (argv[i+1][0]!='\\'))
      fprintf(protocolFile," $%s",(argv[++i]));
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   OpenProto - Open protocol file where specially formatted output is saved

   SYNOPSIS:
   static INT OpenProto (char *name, INT mode);

   PARAMETERS:
   .  name - name of the protocol file
   .  mode - APPEND_PROTO, RENAME_PROTO, TRYRENAME_PROTO, or NORENAME_PROTO

   DESCRIPTION:
   This function opens protocol file where specially formatted output is saved.
   It opens a protocol file for specially formatted output to file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT OpenProto (char *name, INT mode)
{
  char realname[MAXPATHLENGTH],fullname[MAXPATHLENGTH],*pos;
  INT pathlen;
  char c;

  pathlen = 0;
  if (GetDefaultValue(DEFAULTSFILENAME,"protocoldir",fullname)==0)
  {
    pathlen = strlen(fullname);
    strcat(fullname,name);
  }
  else
    strcpy(fullname,name);

  if (protocolFile!=NULL)
  {
    fclose(protocolFile);
    protocolFile = NULL;
    PrintErrorMessage('W',"OpenProto","open protocol file closed!!\n");
  }

  if (mode==APPEND_PROTO)
  {
    protocolFile = fileopen(fullname,"a");
    if (protocolFile==NULL)
      return (1);
    else
      return (0);
  }

  strcpy(realname,fullname);

  if ((mode==RENAME_PROTO)||(mode==TRYRENAME_PROTO))
  {
    /* while file with realname exists */
    c = 'a';
    while ((protocolFile=fileopen(realname,"r"))!=NULL)
    {
      fclose(protocolFile);
      protocolFile = NULL;

      if (c<=MAXRENAMECHAR)
      {
        /* try new name */
        strcpy(realname,fullname);
        if (strchr(name,'.')!=NULL)
        {
          if ((pos=strrchr(realname,'.'))!=NULL)
          {
            /* place a char before .ext in name.ext */
            *pos++ = c++;
            *pos = '\0';
            pos = strrchr(fullname,'.');
            strcat(realname,pos);
          }
        }
        else
        {
          /* place a char after name */
          pos = realname + strlen(realname);
          *pos++ = c++;
          *pos = '\0';
        }
      }
      else if (mode==RENAME_PROTO)
      {
        sprintf(buffer,"could't find a new name for '%s'",fullname);
        PrintErrorMessage('E',"OpenProto",buffer);
        return (1);
      }
      else
        break;
    }
  }

  protocolFile = fileopen(realname,"w");
  if (protocolFile==NULL)
    return (1);

  SetStringVar(":protofilename",realname+pathlen);

  if (strcmp(realname+pathlen,name)!=0)
  {
    sprintf(buffer,"opened protcol file '%s' (instead of '%s')",realname+pathlen,name);
    PrintErrorMessage('W',"OpenProto",buffer);
  }

  return(0);
}

/****************************************************************************/
/*D
   protoOn - open protocol file where specially formatted output is saved

   DESCRIPTION:
   This command opens protocol file where specially formatted output is saved.

   'protoOn <filename> [$r[!] | $a]'

   .   <filename>  - name of the protocol file
   .    $r!        - if a file named <filename> exist already, rename it to <filename>.saved
   .n                          break if the renaming fails

   .   $r                     - like above but proceed even if renaming fails
   .   $a                     - append to existing file named <filename>
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ProtoOnCommand - Open protocol file where specially formatted output is saved

   SYNOPSIS:
   static INT ProtoOnCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function opens protocol file where specially formatted output is saved.
   It writes formatted ouput to the open protocol file.

   protoOn

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ProtoOnCommand (INT argc, char **argv)
{
  static char protoFileName[NAMESIZE];
  INT res,i,RenameMode;

  /* get document name */
  protoFileName[0] = '\0';
  res = sscanf(argv[0],expandfmt(CONCAT3(" protoOn %",NAMELENSTR,"[ -~]")),protoFileName);
  if (res!=1)
  {
    PrintHelp("protoOn",HELPITEM," (filename not found)");
    return (PARAMERRORCODE);
  }

  /* check options */
  RenameMode = NORENAME_PROTO;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (RenameMode!=NORENAME_PROTO)
      {
        PrintErrorMessage('E',"protoOn","specify either $r or $a");
        return (PARAMERRORCODE);
      }
      RenameMode = APPEND_PROTO;
      break;

    case 'r' :
      if (RenameMode!=NORENAME_PROTO)
      {
        PrintErrorMessage('E',"protoOn","specify either $r or $a");
        return (PARAMERRORCODE);
      }
      if (argv[i][1]=='!')
        RenameMode = RENAME_PROTO;
      else
        RenameMode = TRYRENAME_PROTO;
      break;

    default :
      sprintf(buffer," (unknown option '%s')",argv[i]);
      PrintHelp("protoOn",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (OpenProto(protoFileName,RenameMode)>0)
  {
    PrintErrorMessage('E',"protoOn","could not open protocol file");
    return(CMDERRORCODE);
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   protoOff - close protocol file

   DESCRIPTION:
   This command closes the protocol file.

   'protoOff'

   SEE ALSO:
   'protoOn'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ProtoOffCommand - Close protocol file

   SYNOPSIS:
   static INT ProtoOffCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function closes protocol file.
   It closes the protocol file (see the 'protoOn' command).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ProtoOffCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  if (protocolFile==NULL)
  {
    PrintErrorMessage('E',"protoOff","no protocol file open");
    return(PARAMERRORCODE);
  }

  fclose(protocolFile);

  protocolFile = NULL;

  return(OKCODE);
}

FILE *GetProtocolFile (void)
{
  return (protocolFile);
}

/****************************************************************************/
/*D
   logon - open protocol file where all output is saved

   DESCRIPTION:
   This command opens protocol file where all is saved.

   'logon <logfilename>'

   .   <filename>  - name of logfile
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   LogOnCommand - open protocol file where all output is saved

   SYNOPSIS:
   static INT LogOnCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function opens protocol file where all output is saved.
   It directs total shell output also to a logfile.

   logon <logfilename>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT LogOnCommand (INT argc, char **argv)
{
  char logfile[NAMESIZE];
  INT i,rv,popt;

  /* check options */
  popt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'p' :
      if (protocolFile==NULL)
      {
        PrintErrorMessage('E',"logon","no protocol file open");
        return (PARAMERRORCODE);
      }
      popt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("logon",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (popt)
  {
    SetLogFile(protocolFile);
    WriteLogFile("\nbeginlog\n");
    return (OKCODE);
  }

  /* get logfile name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" logon %",NAMELENSTR,"[ -~]")),logfile)!=1)
  {
    PrintErrorMessage('E',"logon","could not read name of logfile");
    return(PARAMERRORCODE);
  }

  rv = OpenLogFile(logfile);
  switch (rv)
  {
  case 0 :
    return (OKCODE);

  case 1 :
    PrintErrorMessage('E',"logon","logfile already open");
    break;

  case 2 :
    PrintErrorMessage('E',"logon","could not open logfile");
    break;

  default :
    PrintErrorMessage('E',"logon","(unknown)");
  }

  return(CMDERRORCODE);
}

/****************************************************************************/
/*D
   logoff - close logfile

   DESCRIPTION:
   This command closes the logfile.

   'logoff'

   SEE ALSO:
   'logon'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   LogOffCommand - Close protocol file

   SYNOPSIS:
   static INT LogOffCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function closes protocol file.
   It switches off the logging mechanism.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT LogOffCommand (INT argc, char **argv)
{
  INT i,popt;

  /* check options */
  popt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'p' :
      if (protocolFile==NULL)
      {
        PrintErrorMessage('E',"logoff","no protocol file open");
        return (PARAMERRORCODE);
      }
      popt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("logon",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (popt)
  {
    WriteLogFile("\nendlog\n");
    SetLogFile(NULL);
    return (OKCODE);
  }

  if (CloseLogFile()!=0)
    PrintErrorMessage('W',"logoff","no logfile open");

  return(OKCODE);
}

/****************************************************************************/
/*D
   new - allocate a new multigrid

   DESCRIPTION:
   This command allocates a new multigrid, using the function 'CreateMultiGrid'.
   It allocates heap and a new multigrid structure.
   The specification of the problem and the domain must be supplied by
   the user with the functions 'CreateProblem' and 'CreateDomain'.
   It also creates the corner vertices and nodes of the domain.

   'new [<mgname>] $d <domain> $p <problem> $f <format> $h <heapsize>'

   .  <mgname>               - the name of the multigrid (default is 'untitled-<nb>')
   .  $d~<domain>            - one of the enroled domains
   .  $p~<problem>           - one of the problems enroled for <domain>
   .  $f~<format>            - one of the enroled formats matching with <problem>
   .  $h~<heapsize>          - the heapsize to be allocated

   EXAMPLE:
   'new $d unit square $p TestProblem $f nc $h 30000000;'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   NewCommand - Allocate a new multigrid

   SYNOPSIS:
   static INT NewCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function allocates a new multigrid.
   It allocates heap and new multigrid structure.

   'new [<mgname>] $d <domain> $p <problem> $f <format> $h <heapsize>'
   .  <mgname>               - the name of the multigrid (default is 'untitled-<nb>')
   .  $d <domain>            - one of the enroled domains
   .  $p <problem>           - one of the problems enroled for <domain>
   .  $f <format>            - one of the enroled formats matching with <problem>
   .  $h <heapsize>          - the heapsize to be allocated

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT NewCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char Multigrid[NAMESIZE],BVPName[NAMESIZE],Format[NAMESIZE];
  unsigned long heapSize;
  INT i,bopt,fopt,hopt;

  /* get multigrid name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" new %",NAMELENSTR,"[ -~]")),Multigrid)!=1) || (strlen(Multigrid)==0))
    sprintf(Multigrid,"untitled-%d",(int)untitledCounter++);

  /* get problem, domain and format */
  heapSize = 0;
  bopt = fopt = hopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'b' :
      if (sscanf(argv[i],expandfmt(CONCAT3("b %",NAMELENSTR,"[ -~]")),BVPName)!=1)
      {
        PrintHelp("new",HELPITEM," (cannot read BndValProblem specification)");
        return(PARAMERRORCODE);
      }
      bopt = TRUE;
      break;

    case 'f' :
      if (sscanf(argv[i],expandfmt(CONCAT3("f %",NAMELENSTR,"[ -~]")),Format)!=1)
      {
        PrintHelp("new",HELPITEM," (cannot read format specification)");
        return(PARAMERRORCODE);
      }
      fopt = TRUE;
      break;

    case 'h' :
      if (sscanf(argv[i],"h %lu",&heapSize)!=1)
      {
        PrintHelp("new",HELPITEM," (cannot read heapsize specification)");
        return(PARAMERRORCODE);
      }
      hopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("new",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (!(bopt && fopt && hopt))
  {
    PrintHelp("new",HELPITEM," (the d, p, f and h arguments are mandatory)");
    return(PARAMERRORCODE);
  }

  /* allocate the multigrid structure */
  theMG = CreateMultiGrid(Multigrid,BVPName,Format,heapSize);
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"new","could not create multigrid");
    return(CMDERRORCODE);
  }

  currMG = theMG;

  return(OKCODE);
}

/****************************************************************************/
/*D
   open - load a new multigrid from a data file

   DESCRIPTION:
   This command loads a new multigrid, using the function 'LoadMultiGrid'.
   Usually, this file should be generated by the 'save' command.
   It allocates the heap and a new multigrid structure.
   The specification of the problem and the domain must be supplied by
   the user with the functions 'CreateProblem' and 'CreateDomain'.

   'open [<mgname>] $d <domain> $p <problem> $f <format> $h <heapsize>'

   .  <mgname>               - the name of the multigrid
   .  $d~<domain>            - one of the enroled domains
   .  $p~<problem>           - one of the problems enroled for <domain>
   .  $f~<format>            - one of the enroled formats matching with <problem>
   .  $h~<heapsize>          - the heapsize to be allocated

   SEE ALSO:
   'new', 'save'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   OpenCommand - load a new multigrid from a data file

   SYNOPSIS:
   static INT OpenCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function loads a new multigrid from a data file.
   It allocates heap and new multigrid structure.

   open <mgname> [$d <domain> [$p <problem>]] [$f <format>] $h <heapsize>

   .   <mgname>               - the name of the multigrid
   .   $d <domain>            - one of the enroled domains
   .   $p <problem>           - one of the problems enroled for <domain>
   .   $f <format>            - one of the enroled formats matching with <problem>
   .   $h <heapsize>          - the heapsize to be allocated

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT OpenCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char theMultigrid[NAMESIZE],Format[NAMESIZE];
  char *theBVP,*theFormat;
  unsigned long heapSize;
  INT i;

  /* get multigrid name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" open %",NAMELENSTR,"[ -~]")),theMultigrid)!=1) || (strlen(theMultigrid)==0))
  {
    PrintErrorMessage('E',"open","specify the name of the multigrid to open");
    return (PARAMERRORCODE);
  }

  /* get problem, domain and format */
  theBVP = theFormat = NULL;
  heapSize = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'b' :
      if (sscanf(argv[i],expandfmt(CONCAT3("b %",NAMELENSTR,"[ -~]")),theBVP)!=1)
      {
        PrintHelp("new",HELPITEM," (cannot read BndValProblem specification)");
        return(PARAMERRORCODE);
      }
      break;

    case 'f' :
      if (sscanf(argv[i],expandfmt(CONCAT3("f %",NAMELENSTR,"[ -~]")),Format)!=1)
      {
        PrintHelp("open",HELPITEM," (cannot read format specification)");
        return(PARAMERRORCODE);
      }
      theFormat = Format;
      break;

    case 'h' :
      if (sscanf(argv[i],"h %lu",&heapSize)!=1)
      {
        PrintHelp("open",HELPITEM," (cannot read heapsize specification)");
        return(PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("open",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (heapSize==0)
  {
    PrintErrorMessage('E',"open","heapsize not specified");
    return(CMDERRORCODE);
  }

  /* allocate the multigrid structure */
  theMG = LoadMultiGrid(theMultigrid,theMultigrid,theBVP,theFormat,heapSize);
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"open","could not open multigrid");
    return(CMDERRORCODE);
  }

  currMG = theMG;

  return(OKCODE);
}

/****************************************************************************/
/*D
   close - close current multigrid

   DESCRIPTION:
   This command closes the current (or all) open multigrid(s),
   frees their heaps and closes all the pictures belonging to them,
   calling 'DisposeMultiGrid' and 'DisposePicture'.

   'close [$a]'

   .   $a                     - close all multigrids
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   CloseCommand	- Close current multigrid

   SYNOPSIS:
   static INT CloseCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function closes current multigrid.
   It closes the current (or all) open multigrid(s) and frees their heaps.

   close [$a]
   .   $a                     - close all multigrids

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT CloseCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  UGWINDOW *theWin;
  PICTURE *thePic,*NextPic,*currPic;
  INT i,closeonlyfirst;

  closeonlyfirst = TRUE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      closeonlyfirst = FALSE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("close",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  i = 0;
  do
  {
    /* get multigrid */
    theMG = currMG;
    if (theMG==NULL)
    {
      if (i==0)
      {
        PrintErrorMessage('W',"close","no open multigrid");
        return (OKCODE);
      }
      closeonlyfirst = FALSE;
      break;
    }

    currPic = GetCurrentPicture();

    /* close pictures first */
    for (theWin=GetFirstUgWindow(); theWin!=NULL; theWin=GetNextUgWindow(theWin))
      for (thePic=GetFirstPicture(theWin); thePic!=NULL; thePic=NextPic)
      {
        NextPic = GetNextPicture(thePic);

        if (PIC_MG(thePic)==theMG)
        {
          if (thePic==currPic)
            SetCurrentPicture(NULL);
          if (DisposePicture(thePic)!=0)
          {
            PrintErrorMessage('E',"closewindow","could not close a picture of that window");
            return (CMDERRORCODE);
          }
        }
      }

    if (DisposeMultiGrid(theMG)!=0)
    {
      PrintErrorMessage('E',"close","closing the mg failed");
      return (CMDERRORCODE);
    }
    i++;

    currMG = GetFirstMultigrid();
  }
  while (!closeonlyfirst);

  return(OKCODE);
}

/****************************************************************************/
/*D
   save - save a multigrid structure in a file

   DESCRIPTION:
   This command writes the current multigrid structure in a file.

   'save [<name>] [$c <comment>]'

   .  <name>                  - name to save with (default is the mgname)
   .  $c <comment>            - optionally specify a comment string

   SEE ALSO:
   'open'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SaveCommand - Save multigrid structure in a file

   SYNOPSIS:
   static INT SaveCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function saves multigrid structure in a file.
   It saves the current multigrid.

   save [<name>] [$c <comment>]

   .  <name>                  - name to save with (default is the mgname)
   .  $c <comment>            - optionally specify a comment string

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SaveCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char Name[NAMESIZE],Comment[LONGSTRSIZE];
  INT i;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"save","no open multigrid");
    return (CMDERRORCODE);
  }

  /* scan name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" save %",NAMELENSTR,"[ -~]")),Name)!=1)
    strcpy(Name,ENVITEM_NAME(theMG));

  /* check options */
  strcpy(Comment,NO_COMMENT);
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'c' :
      if (sscanf(argv[i],expandfmt(CONCAT3(" c %",LONGSTRLENSTR,"[ -~]")),Comment)!=1)
      {
        PrintErrorMessage('E',"save","couldn't read the comment string");
        return (PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("save",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  SaveMultiGrid(theMG,Name,Comment);

  return(OKCODE);
}

/****************************************************************************/
/*D
   level - select another level

   DESCRIPTION:
   This command selects another level or lists an info.
   It changes the working (current) level of the current multigrid.

   level <level> | + | -

   .  <level> - go to level <level>
   .  +       - go to the next finer level
   .  -       - go to the next coarser level
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   LevelCommand	- Select another level or list info

   SYNOPSIS:
   static INT LevelCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function selects another level or list info.
   It changes the working (current) level of the current multigrid.

   level <level> | + | -

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT LevelCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

  /* following variables: keep type for sscanf */
  int l;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"level","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check parameters */
  if (sscanf(argv[0]," level %d",&l)==1)
  {
    if ((l<0) || (l>TOPLEVEL(theMG)))
    {
      PrintErrorMessage('E',"level","level out of range");
      return (PARAMERRORCODE);
    }
    else
      CURRENTLEVEL(theMG) = l;
  }
  else if (strchr(argv[0],'+')!=NULL)
  {
    if (CURRENTLEVEL(theMG)==TOPLEVEL(theMG))
    {
      PrintErrorMessage('W',"level","already on TOPLEVEL");
      return (OKCODE);
    }
    else
      CURRENTLEVEL(theMG)++;
  }
  else if (strchr(argv[0],'-')!=NULL)
  {
    if (CURRENTLEVEL(theMG)==0)
    {
      PrintErrorMessage('W',"level","already on level 0");
      return (OKCODE);
    }
    else
      CURRENTLEVEL(theMG)--;
  }
  else
  {
    PrintErrorMessage('E',"level","specify <level>, + or - with the level command");
    return (CMDERRORCODE);
  }

  UserWriteF("  current level is %d (top level %d)\n",CURRENTLEVEL(theMG),TOPLEVEL(theMG));

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return(OKCODE);
}

/****************************************************************************/
/*D
   renumber - reassign the object IDs in the multigrid


   DESCRIPTION:
   This command reassigns the object IDs in the multigrid
   subsequently to fill the gaps, calling the function 'RenumberMultiGrid'.

   'renumber'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   RenumberMGCommand - Reassign the object IDs in the multigrid

   SYNOPSIS:
   static INT RenumberMGCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function reassigns the object IDs in the multigrid
   subsequently to fill the gaps.
   It assigns consecutive id's to all objects.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT RenumberMGCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"renumber","no open multigrid");
    return (CMDERRORCODE);
  }

  if (RenumberMultiGrid(theMG)!=GM_OK)
  {
    PrintErrorMessage('E',"renumber","renumbering of the mg failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   wplist - list information on all windows and pictures

   DESCRIPTION:
   This command lists information on all windows and pictures, calling
   the functions 'ListWindowPictureHeader', 'ListPicture' and
   'ListUGWindow'.

   'wplist'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   WindowPictureListCommand - List info of all windows and pictures

   SYNOPSIS:
   static INT WindowPictureListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists info of all windows and pictures.
   It lists all windows with their pictures.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT WindowPictureListCommand (INT argc, char **argv)
{
  UGWINDOW *currUgW, *theUgW;
  PICTURE *currPic, *thePic;

  NO_OPTION_CHECK(argc,argv);

  currUgW = GetCurrentUgWindow();
  currPic = GetCurrentPicture();

  ListWindowPictureHeader();
  for (theUgW=GetFirstUgWindow(); theUgW!=NULL; theUgW=GetNextUgWindow(theUgW))
  {
    ListUgWindow(theUgW,(theUgW==currUgW));
    for (thePic=GetFirstPicture(theUgW); thePic!=NULL; thePic=GetNextPicture(thePic))
      ListPicture(thePic,(thePic==currPic));
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   mglist - list information on all multigrids

   DESCRIPTION:
   This command lists information on all multigrids, calling
   the functions 'ListMultiGridHeader' and 'ListMultiGrid'.

   'mglist [$l]'

   .  $l - long format for additional information
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   MGListCommand - List info all multigrids

   SYNOPSIS:
   static INT MGListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists info all multigrids.
   It lists all open multigrids.

   mglist [$l]
   .  $l                     - list long format

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MGListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG,*theCurrMG;
  INT i,longformat;

  theCurrMG = GetCurrentMultigrid();
  if (theCurrMG==NULL)
  {
    PrintErrorMessage('W',"mglist","no multigrid open\n");
    return (OKCODE);
  }

  longformat = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'l' :
      longformat = TRUE;
      break;

    default :
      sprintf(buffer," (unknown option '%s')",argv[i]);
      PrintHelp("mglist",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  ListMultiGridHeader(longformat);

  for (theMG=GetFirstMultigrid(); theMG!=NULL; theMG=GetNextMultigrid(theMG))
    ListMultiGrid(theMG,(theMG==theCurrMG),longformat);

  return (OKCODE);
}

/****************************************************************************/
/*D
   glist - list information on the current multigrid

   DESCRIPTION:
   This command lists information on the current multigrid, calling
   the function 'ListGrids'.

   'glist'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   GListCommand	- List info for the current multigrid

   SYNOPSIS:
   static INT GListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists info for the current multigrid.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT GListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    UserWrite("no multigrid open\n");
    return (OKCODE);
  }

  ListGrids(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   nlist - list information on specified nodes

   DESCRIPTION:
   This command lists information on specified nodes, calling
   the functions 'ListNodeRange' and 'ListNodeSelection'.

   'nlist $s | {$i <fromID> [<toID>]} [$d] [$b] [$n] [$v] [$a]'

   .  $s                     - list info for the selected nodes
   .  $i                     - list info for nodes with an ID in the range <fromID> through <toID>                           # if <fromID> is omitted only the node with <fromID> is listed
   .  $d                     - up to version 2.3 ONLY: list also user data space
   .  $b                     - print additional info for boundary nodes
   .  $n                     - list also neighbours of each node
   .  $v                     - print extended info (verbose mode)
   .  $a                     - list all nodes
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   NListCommand - List info for specified nodes

   SYNOPSIS:
   static INT NListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists info for specified nodes.
   It prints info for the specified nodes.

   nlist $s | {$i <fromID> [<toID>]} [$d] [$b] [$n] [$v] [$a]

   .  $s                     - list info for the selected nodes
   .  $i                     - list info for nodes with an ID in the range <fromID> through <toID>                           # if <fromID> is omitted only the node with <fromID> is listed
   .  $d                     - up to version 2.3 ONLY: list also user data space
   .  $b                     - print additional info for boundary nodes
   .  $n                     - list also neighbours of each node
   .  $v                     - print extended info (verbose mode)
   .  $a                     - list all nodes

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT NListCommand (INT argc, char **argv)
{

  MULTIGRID *theMG;
  INT i,fromN,toN,res,mode,dataopt,boundaryopt,neighbouropt,verboseopt;

  /* following variables: keep type for sscanf */
  long f,t;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"nlist","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  dataopt = boundaryopt = neighbouropt = verboseopt = mode = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'i' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"nlist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ID;
      res = sscanf(argv[i]," i %ld %ld",&f,&t);
      fromN = f;
      toN   = t;
      if (res<1)
      {
        PrintErrorMessage('E',"nlist","specify at least one id with the i option");
        return (PARAMERRORCODE);
      }
      else if (res==1)
        toN = fromN;
      else if (fromN>toN)
      {
        PrintErrorMessage('E',"nlist","from ID > to ID");
        return (PARAMERRORCODE);
      }
      break;

    case 's' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"nlist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_SELECTION;
      break;

    case 'a' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"nlist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ALL;
      break;

    case 'd' :
      dataopt = TRUE;
      break;

    case 'b' :
      boundaryopt = TRUE;
      break;

    case 'n' :
      neighbouropt = TRUE;
      break;

    case 'v' :
      verboseopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("nlist",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  switch (mode)
  {
  case DO_ID :
    ListNodeRange(theMG,fromN,toN,dataopt,boundaryopt,neighbouropt,verboseopt);
    break;

  case DO_ALL :
    ListNodeRange(theMG,0,MAX_I,dataopt,boundaryopt,neighbouropt,verboseopt);
    break;

  case DO_SELECTION :
    ListNodeSelection(theMG,dataopt,boundaryopt,neighbouropt,verboseopt);
    break;

  default :
    PrintErrorMessage('E',"nlist","specify either the a, s or i option");
    return (PARAMERRORCODE);
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   elist - list information on specified elements

   DESCRIPTION:
   This command lists information on spezified elements, calling
   the functions 'ListElementRange' and 'ListElementSelection'.

   'elist $s | {$i <fromID> [<toID>]} [$d] [$b] [$n] [$v] [$a]'

   .  $s                     - list info for the selected elements
   .  $i                     - list info for elements with an ID in the range <fromID> through <toID>
   .n                        - if <fromID> is omitted only the element with <fromID> is listed

   .  $d                     - up to version 2.3 ONLY: list also user data space
   .  $b                     - print additional info for boundary elements
   .  $n                     - list also neighbours of each element
   .  $v                     - print extended info (verbose mode)
   .  $l                     - list only elements of current level
   .  $a                     - list all elements
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   EListCommand	- List info for specified elements

   SYNOPSIS:
   static INT EListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists info for specified elements.
   It prints info for the specified elements.

   elist $s | {$i <fromID> [<toID>]} [$d] [$b] [$n] [$v] [$a]

   .  $s                     - list info for the selected elements
   .  $i                     - list info for elements with an ID in the range <fromID> through <toID>
   .n                        if <fromID> is omitted only the element with <fromID> is listed

   .  $d                     - up to version 2.3 ONLY: list also user data space
   .  $b                     - print additional info for boundary elements
   .  $n                     - list also neighbours of each element
   .  $v                     - print extended info (verbose mode)
   .  $l                     - list only elements of current level
   .  $a                     - list all elements

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT EListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,fromE,toE,res,mode,dataopt,boundaryopt,neighbouropt,verboseopt,levelopt;

  /* following variables: keep type for sscanf */
  long f,t;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"elist","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  dataopt = boundaryopt = neighbouropt = verboseopt = levelopt = mode = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'i' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"elist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ID;
      res = sscanf(argv[i]," i %ld %ld",&f,&t);
      fromE = f;
      toE   = t;
      if (res<1)
      {
        PrintErrorMessage('E',"elist","specify at least one id with the i option");
        return (PARAMERRORCODE);
      }
      else if (res==1)
        toE = fromE;
      else if (fromE>toE)
      {
        PrintErrorMessage('E',"elist","from ID > to ID");
        return (PARAMERRORCODE);
      }
      break;

    case 's' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"elist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_SELECTION;
      break;

    case 'a' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"elist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ALL;
      break;

    case 'l' :
      levelopt = TRUE;
      break;

    case 'd' :
      dataopt = TRUE;
      break;

    case 'b' :
      boundaryopt = TRUE;
      break;

    case 'n' :
      neighbouropt = TRUE;
      break;

    case 'v' :
      verboseopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("elist",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  switch (mode)
  {
  case DO_ID :
    ListElementRange(theMG,fromE,toE,dataopt,boundaryopt,neighbouropt,verboseopt,levelopt);
    break;

  case DO_ALL :
    ListElementRange(theMG,0,MAX_I,dataopt,boundaryopt,neighbouropt,verboseopt,levelopt);
    break;

  case DO_SELECTION :
    ListElementSelection(theMG,dataopt,boundaryopt,neighbouropt,verboseopt);
    break;

  default :
    PrintErrorMessage('E',"elist","specify either the a, s or i option");
    return (PARAMERRORCODE);
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   slist - list information on all selected nodes and elements

   DESCRIPTION:
   This command lists information on selected nodes and elements, calling
   the functions 'ListNodeSelection' and 'ListElementSelection'.

   'slist [$d] [$b] [$n] [$v]'

   .   $d                     - up to version 2.3 ONLY: list also user data space
   .   $b                     - print additional info for boundary nodes/elements
   .   $n                     - list also neighbours of each node/element
   .   $v                     - print extended info (verbose mode)

   SEE ALSO:
   'select'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SelectionListCommand - List all nodes/elements from selection buffer

   SYNOPSIS:
   static INT SelectionListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists all nodes/elements from selection buffer.
   It lists the contents of the selction buffer (elist/nlist format resp.).

   slist [$d] [$b] [$n] [$v]
   .   $d                     - up to version 2.3 ONLY: list also user data space
   .   $b                     - print additional info for boundary nodes/elements
   .   $n                     - list also neighbours of each node/element
   .   $v                     - print extended info (verbose mode)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SelectionListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,dataopt,boundaryopt,neighbouropt,verboseopt;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"slist","no open multigrid");
    return (CMDERRORCODE);
  }

  if (SELECTIONSIZE(theMG)==0)
  {
    PrintErrorMessage('W',"slist","nothing selected");
    return (OKCODE);
  }

  /* check options */
  dataopt = boundaryopt = neighbouropt = verboseopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'd' :
      dataopt = TRUE;
      break;

    case 'b' :
      boundaryopt = TRUE;
      break;

    case 'n' :
      neighbouropt = TRUE;
      break;

    case 'v' :
      verboseopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("slist",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  switch (SELECTIONMODE(theMG))
  {
  case elementSelection :
    ListElementSelection(theMG,dataopt,boundaryopt,neighbouropt,verboseopt);
    break;

  case nodeSelection :
    ListNodeSelection(theMG,dataopt,boundaryopt,neighbouropt,verboseopt);
    break;

  default :
    PrintErrorMessage('W',"slist","selectionmode ???");
    return (PARAMERRORCODE);
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   vmlist - list information on specified vectors and matrices

   DESCRIPTION:
   This command lists information on specified vectors and matrices, calling
   the functions 'ListVectorRange' and 'ListVectorSelection'.

   'vmlist $s | {$i <fromID> [<toID>]} [$m] [$d] [$a] [$l <f> <t>]'

   .  $s                     - list info for the selected vectors
   .  $i                     - list info for vectors with an ID in the range <fromID> through <toID>
   .n                        if <fromID> is omitted only the vector with <fromID> is listed

   .  $m                     - list also the associated matrix entries
   .  $d                     - list also the user data
   .  $a                     - list all vectors
   .  $l <f> <t>             - process levels f <= l <= t
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   VMListCommand - List info for specified vectors (and matrices)

   SYNOPSIS:
   static INT VMListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists info for specified vectors (and matrices).
   It prints info for the specified vectors (possibly including associated matrix entries).

   vmlist $s | {$i <fromID> [<toID>]} [$m] [$d] [$a] [$l <f> <t>]

   .  $s                     - list info for the selected vectors
   .  $i                     - list info for vectors with an ID in the range <fromID> through <toID>
   .n                        if <fromID> is omitted only the vector with <fromID> is listed

   .  $m                     - list also the associated matrix entries
   .  $d                     - list also the user data
   .  $a                     - list all vectors
   .  $l <f> <t>             - process levels f <= l <= t

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT VMListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,fl,tl,fromV,toV,res,mode,dataopt,matrixopt;

  /* following variables: keep type for sscanf */
  long f,t;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"vmlist","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  dataopt = matrixopt = mode = FALSE;
  fl = tl = CURRENTLEVEL(theMG);
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'i' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"vmlist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ID;
      res = sscanf(argv[i]," i %ld %ld",&f,&t);
      fromV = f;
      toV   = t;
      if (res<1)
      {
        PrintErrorMessage('E',"vmlist","specify at least one id with the i option");
        return (PARAMERRORCODE);
      }
      else if (res==1)
        toV = fromV;
      else if (fromV>toV)
      {
        PrintErrorMessage('E',"vmlist","from ID > to ID");
        return (PARAMERRORCODE);
      }
      break;

    case 's' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"vmlist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_SELECTION;
      break;

    case 'a' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"vmlist","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ALL;
      break;

    case 'l' :
      res = sscanf(argv[i]," l %ld %ld",&f,&t);
      fl = f;
      tl = t;
      if (res!=2)
      {
        PrintErrorMessage('E',"vmlist","specify from and to level with the l option");
        return (PARAMERRORCODE);
      }
      else if (fl>tl)
      {
        PrintErrorMessage('E',"vmlist","from level > to level");
        return (PARAMERRORCODE);
      }
      break;

    case 'd' :
      dataopt = TRUE;
      break;

    case 'm' :
      matrixopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("vmlist",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  switch (mode)
  {
  case DO_ID :
    ListVectorRange(theMG,fl,tl,fromV,toV,matrixopt,dataopt);
    break;

  case DO_ALL :
    ListVectorRange(theMG,fl,tl,0,MAX_I,matrixopt,dataopt);
    break;

  case DO_SELECTION :
    if (SELECTIONMODE(theMG)==elementSelection)
      ListVectorOfElementSelection(theMG,matrixopt,dataopt);
    else
      ListVectorSelection(theMG,matrixopt,dataopt);
    break;

  default :
    PrintErrorMessage('E',"vmlist","specify either the a, s or i option");
    return (PARAMERRORCODE);
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   in - insert an inner node and vertex

   DESCRIPTION:
   This command inserts an inner node and the corresponding vertex
   into a multigrid with only level 0, calling the function
   'InsertInnerNode'.

   'in <x> <y> [<z>]'

   .  <x>~<y>~[<z>]            - specify as much coordinates as the space has dimensions
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   InsertInnerNodeCommand - Insert an inner node+vertex

   SYNOPSIS:
   static INT InsertInnerNodeCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function inserts an inner node+vertex.
   It inserts an inner node into a multigrid with only level 0.

   in <x> <y> <z>

   .  <x> <y> <z>            - specify as much coordinates as the space has dimensions

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InsertInnerNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  COORD xc[DIM];
  INT i;

  /* following variables: keep type for sscanf */
  float x[DIM_MAX];

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"in","no open multigrid");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"in %f %f %f",x,x+1,x+2)!=DIM)
  {
    sprintf(buffer,"specify %d coordinates for an inner node",(int)DIM);
    PrintErrorMessage('E',"in",buffer);
    return (PARAMERRORCODE);
  }
  for (i=0; i<DIM; i++)
    xc[i] = x[i];

  /* NB: toplevel=0 is checked by InsertInnerNode() */
  if (InsertInnerNode(theMG,xc)!=GM_OK)
  {
    PrintErrorMessage('E',"in","inserting an inner node failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   bn - insert a boundary node and vertex

   DESCRIPTION:
   This command inserts an boundary node and the corresponding vertex
   into a multigrid with only level 0, calling the function
   'InsertBoubdaryNode'.

   'bn <Id> <s> [<t>]'

   .  <Id>                   - insert a boundary node on the boundary segment with <Id>
   .  <s>~[<t>]                - specify as much boundary segment coordinates as the boundary has dimensions
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   InsertBoundaryNodeCommand - Insert a boundary node+vertex

   SYNOPSIS:
   static INT InsertBoundaryNodeCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function inserts a boundary node+vertex.
   It inserts a boundary node into a multigrid with only level 0.

   bn <Id> <s> <t>

   .  <Id>                   - insert a boundary node on the boundary segment with <Id>
   .  <s> <t>                - specify as much boundary segment coordinates as the boundary has dimensions

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InsertBoundaryNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  COORD xc[DIM_OF_BND];
  INT i;

  /* following variables: keep type for sscanf */
  float x[DIM_OF_BND_MAX];
  int segid;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"bn","no open multigrid");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"bn %d %f %f",&segid,x,x+1)!=1+DIM_OF_BND)
  {
    sprintf(buffer,"specify %d coordinates for a boundary node",(int)DIM);
    PrintErrorMessage('E',"bn",buffer);
    return (PARAMERRORCODE);
  }
  for (i=0; i<DIM_OF_BND; i++)
    xc[i] = x[i];

  /* NB: toplevel=0 is checked by InsertBoundaryNode() */
  if (InsertBoundaryNode(theMG,segid,xc)!=GM_OK)
  {
    PrintErrorMessage('E',"bn","inserting a boundary node failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   deln - delete a node and vertex

   DESCRIPTION:
   This command deletes a node and the corresponding vertex
   of the current multigrid, calling the function 'DeleteNode'.

   'deln <Id> | $s'

   .  <Id>                   - ID of the node to be deleted
   .  $s                     - delete ALL nodes from the selection
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   DeleteNodeCommand - Delete a node+vertex

   SYNOPSIS:
   static INT DeleteNodeCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function deletes a node+vertex.
   It deletes a node from the multigrid.

   deln <Id> | $s

   .  <Id>                   - ID of the node to be deleted
   .  $s                     - delete ALL nodes from the selection

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT DeleteNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT sopt,i;

  /* following variables: keep type for sscanf */
  int id;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"deln","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  sopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      sopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("deln",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (sopt)
  {
    if (SELECTIONMODE(theMG)==nodeSelection)
      for (i=0; i<SELECTIONSIZE(theMG); i++)
        if (DeleteNode(theMG,(NODE *)SELECTIONOBJECT(theMG,i))!=GM_OK)
        {
          PrintErrorMessage('E',"deln","deleting the node failed");
          return (CMDERRORCODE);
        }
    ClearSelection(theMG);
    InvalidatePicturesOfMG(theMG);
    InvalidateUgWindowsOfMG(theMG);

    return (OKCODE);
  }

  if (sscanf(argv[0],"deln %d",&id)!=1)
  {
    PrintErrorMessage('E',"deln","specify the ID of the node to be deleted");
    return (PARAMERRORCODE);
  }

  /* NB: toplevel=0 is checked by DeleteNode() */
  if (DeleteNodeWithID(theMG,id)!=GM_OK)
  {
    PrintErrorMessage('E',"deln","deleting the node failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   move - move a node and vertex

   DESCRIPTION:
   This command moves a node and the corresponding vertex
   of the current multigrid to a new position,
   calling the functions 'MoveInnerNode' and 'MoveBoundaryNode'.

   'move {<Id> | $s} {$i <x> <y> [<z>] | $b <SegId> <s> [<t>]}'

   .  <Id>                   - Id of the node to be moved
   .  $i~<x>~<y>~[<z>]         - specify as much coordinates as the space has dimensions
   .  $b~<Id>~<s>~[<t>]        - specify as much boundary segment coordinates as the boundary has dimensions
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   MoveNodeCommand - Move a node of the current multigrid

   SYNOPSIS:
   static INT MoveNodeCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function moves a node of the current multigrid.
   It moves a node of the multigrid to a new postion.

   move {<Id> | $s} {$i <x> <y> <z> | $b <SegId> <s> <t>}

   .  <Id>                   - Id of the node to be moved
   .  $i <x> <y> <z>         - specify as much coordinates as the space has dimensions
   .  $b <Id> <s> <t>        - specify as much boundary segment coordinates as the boundary has dimensions

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MoveNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  NODE *theNode;
  COORD xc[DIM];
  INT type,i,j,level;

  /* following variables: keep type for sscanf */
  float x[DIM_MAX];
  int id,segid;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"move","no open multigrid");
    return (CMDERRORCODE);
  }

  theNode = NULL;

  /* check parameters */
  if (sscanf(argv[0],"move %d",&id)==1)
  {
    /* search node */
    for (level=0; level<=TOPLEVEL(theMG); level++)
      if ((theNode=FindNodeFromId(GRID_ON_LEVEL(theMG,level),id))!=NULL)
        break;
    if (theNode==NULL)
    {
      sprintf(buffer,"node with ID %ld not found",(long)id);
      PrintErrorMessage('E',"move",buffer);
      return (CMDERRORCODE);
    }
  }

  /* check options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      if (SELECTIONMODE(theMG)==nodeSelection)
      {
        PrintErrorMessage('E',"move","there is no node in the selection");
        return (PARAMERRORCODE);
      }
      if (SELECTIONSIZE(theMG)!=1)
      {
        PrintErrorMessage('E',"move","there is more than one node in the selection");
        return (PARAMERRORCODE);
      }
      theNode = (NODE *)SELECTIONOBJECT(theMG,0);
      break;

    case 'i' :
      if (OBJT(MYVERTEX(theNode))!=IVOBJ)
      {
        sprintf(buffer,"node with ID %ld is no inner node",(long)id);
        PrintErrorMessage('E',"move",buffer);
        return (CMDERRORCODE);
      }
      type = IVOBJ;
      if (sscanf(argv[i],"i %f %f %f",x,x+1,x+2)!=DIM)
      {
        sprintf(buffer,"specify %d new coordinates for an inner node",(int)DIM);
        PrintErrorMessage('E',"move",buffer);
        return (PARAMERRORCODE);
      }
      for (j=0; j<DIM; j++)
        xc[j] = x[j];
      break;

    case 'b' :
      if (OBJT(MYVERTEX(theNode))!=BVOBJ)
      {
        sprintf(buffer,"node with ID %ld is no boundary node",(long)id);
        PrintErrorMessage('E',"move",buffer);
        return (CMDERRORCODE);
      }
      type = BVOBJ;
      if (sscanf(argv[i],"b %d %f %f",&segid,x,x+1)!=1+DIM_OF_BND)
      {
        sprintf(buffer,"specify the segment if and %d new coordinates for a boundary node",(int)DIM_OF_BND);
        PrintErrorMessage('E',"move",buffer);
        return (PARAMERRORCODE);
      }
      for (j=0; j<DIM_OF_BND; j++)
        xc[j] = x[j];
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("move",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (theNode==NULL)
  {
    PrintErrorMessage('E',"move","you have to either specify\n"
                      "the ID of the node to move or the s option");
    return (PARAMERRORCODE);
  }

  if (type==IVOBJ)
  {
    if (MoveInnerNode(theMG,theNode,xc)!=GM_OK)
    {
      PrintErrorMessage('E',"move","failed moving the node");
      return (CMDERRORCODE);
    }
  }
  else if (type==BVOBJ)
  {
    if (MoveBoundaryNode(theMG,theNode,segid,xc)!=GM_OK)
    {
      PrintErrorMessage('E',"move","failed moving the node");
      return (CMDERRORCODE);
    }
  }
  else
  {
    PrintHelp("move",HELPITEM," (either i or b option is mandatory)");
    return (PARAMERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   ie - insert an element

   DESCRIPTION:
   This command inserts an element into a multigrid with only level 0,
   calling the function 'InsertElement'.


   'ie {<Id>}+ | $s'

   .  {<Id>}+                - specify at least three (2d) or four (3d) corner nodes, the corresponding (unique) element will be created
   .  $s - taking selected nodes
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   InsertElementCommand - Insert an element

   SYNOPSIS:
   static INT InsertElementCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function inserts an element.
   It inserts an element for specified corner nodes.

   ie {<Id>}+ | $s

   .  {<Id>}+                - specify at least three corner nodes, the coresponding (unique) element will
   .n                        - be created

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InsertElementCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  NODE *theNode,*theNodes[MAX_CORNERS_OF_ELEM];
  char *token,*vstr;
  INT i,nNodes,Id[MAX_CORNERS_OF_ELEM];

  /* following variables: keep type for sscanf */
  int id;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"ie","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  nNodes = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      if (SELECTIONMODE(theMG)==nodeSelection)
        for (i=0; i<SELECTIONSIZE(theMG); i++)
        {
          theNode = (NODE *)SELECTIONOBJECT(theMG,i);
          if (nNodes>=MAX_CORNERS_OF_ELEM)
          {
            PrintErrorMessage('E',"ie","too many nodes are in the selection");
            return (CMDERRORCODE);
          }
          theNodes[nNodes++] = theNode;
        }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("ie",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* got the nodes via s option? */
  if (nNodes>0)
  {
    if (InsertElement(theMG,nNodes,theNodes)!=GM_OK)
    {
      PrintErrorMessage('E',"ie","inserting the element failed");
      return (CMDERRORCODE);
    }
    else
    {
      InvalidatePicturesOfMG(theMG);
      InvalidateUgWindowsOfMG(theMG);
      return (OKCODE);
    }
  }

  /* set vstr after 'ie' */
  if ((vstr=strchr(argv[0],'e'))!=NULL)
    ++vstr;
  else
    return (CMDERRORCODE);

  /* we need the string split into tokens */
  nNodes = 0;
  token = strtok(vstr,WHITESPACE);
  while (token!=NULL)
  {
    if (nNodes>=MAX_CORNERS_OF_ELEM)
    {
      sprintf(buffer,"specify at most %d id's",(int)MAX_CORNERS_OF_ELEM);
      PrintErrorMessage('E',"ie",buffer);
      return (PARAMERRORCODE);                                  /* too many items */
    }
    if (sscanf(token," %d",&id)!=1)
    {
      sprintf(buffer,"could not read the id of corner no %d",(int)i);
      PrintErrorMessage('E',"ie",buffer);
      return (PARAMERRORCODE);                                  /* too many items */
    }
    Id[nNodes++] = id;
    token = strtok(NULL,WHITESPACE);
  }

  /* NB: toplevel=0 is checked by InsertElementFromIDs() */
  if (InsertElementFromIDs(theMG,nNodes,Id)!=GM_OK)
  {
    PrintErrorMessage('E',"ie","inserting the element failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   dele - delete an element

   DESCRIPTION:
   This command deletes the specified  element of a multigrid
   with only level 0, calling the function 'DeleteElement'.

   'dele <Id> | $s'

   .  <Id>                   - ID of the element to be deleted
   .  $s                     - delete all elements from the selection
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   DeleteElementCommand - Delete an element

   SYNOPSIS:
   static INT DeleteElementCommand (INT argc, char **argv)

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function deletes an element.
   It deletes an element from the multigrid including edges not deeded anymore.

   dele <Id> | $s

   .  <Id>                   - ID of the element to be deleted
   .  $s                     - delete ALL elements from the selection

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT DeleteElementCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,sopt;

  /* following variables: keep type for sscanf */
  int id;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"dele","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  sopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      sopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("dele",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (sopt)
  {
    if (SELECTIONMODE(theMG)==elementSelection)
      for (i=0; i<SELECTIONSIZE(theMG); i++)
        if (DeleteElement(theMG,(ELEMENT *)SELECTIONOBJECT(theMG,i))!=GM_OK)
        {
          PrintErrorMessage('E',"dele","deleting the element failed");
          return (CMDERRORCODE);
        }
    ClearSelection(theMG);
    InvalidatePicturesOfMG(theMG);
    InvalidateUgWindowsOfMG(theMG);

    return (OKCODE);
  }

  if (sscanf(argv[0],"dele %d",&id)!=1)
  {
    PrintErrorMessage('E',"dele","specify the ID of the element to be deleted");
    return (PARAMERRORCODE);
  }

  /* NB: toplevel=0 is checked by DeleteElement() */
  if (DeleteElementWithID(theMG,id)!=GM_OK)
  {
    PrintErrorMessage('E',"dele","deleting the element failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   refine - refine the current multigrid

   DESCRIPTION:
   This command refines the multigrid according to the refinement marks
   set in the elements, calling the function 'RefineMultiGrid'.

   'refine [$g]'

   .  no~option - only local refinement
   .  $g - copy nonrefined regions to new level
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   RefineCommand - Refine the current multigrid

   SYNOPSIS:
   static INT RefineCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function refines the current multigrid.
   It refines the multigrid according to the refinement marks set in the elements.

   refine [$g]

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT RefineCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,mode,rv;
  char buffer[128];
  EVECTOR *theElemEvalDirection;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"refine","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  theElemEvalDirection = NULL;
  mode = GM_REFINE_TRULY_LOCAL;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'g' :
      mode = GM_COPY_ALL;
      break;
#ifdef __THREEDIM__
    case 'a' :
      if (sscanf(argv[i],"a %s",buffer)==1)
        theElemEvalDirection = GetElementVectorEvalProc(buffer);
      if (theElemEvalDirection==NULL)
        UserWrite("direction eval fct not found: taking shortest interior edge\n");
      break;
#endif
    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("refine",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  rv = RefineMultiGrid(theMG,mode,theElemEvalDirection);

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  switch (rv)
  {
  case GM_OK :
    sprintf(buffer," %s refined\n",ENVITEM_NAME(theMG));
    UserWrite(buffer);
    SetStringVar(":errno","0");
    return (OKCODE);

  case GM_ERROR :
    PrintErrorMessage('E',"refine","could not refine, data structure still ok");
    SetStringVar(":errno","1");
    return (CMDERRORCODE);

  case GM_FATAL :
    PrintErrorMessage('F',"refine","could not refine, data structure inconsistent\n");
    SetStringVar(":errno","1");
    return (CMDERRORCODE);

  default :
    PrintErrorMessage('E',"refine","unknown error in refine");
    SetStringVar(":errno","1");
    return (CMDERRORCODE);
  }
}

/****************************************************************************/
/*D
   mark - mark elements with refinement type

   DESCRIPTION:
   This command marks elements with refinement type,
   calling the function 'MarkForRefinement'.

   'mark [$h | {[<rule> [<side>]] [$a | $i <Id> | $s]} | $c]'

   .  <rule>                 - specify a refinement rule ("red" is default)
   .  <side>                 - has to be specified if the corresponding rule can be applied in several orientations
   .  $a                     - refine all (leave) elements
   .  $i <Id>                - refine the element with <Id>
   .  $s                     - refine all elements from the current selection
   .  $h                     - show available rules
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   MarkCommand - Mark element with refinement type

   SYNOPSIS:
   static INT MarkCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function marks element with refinement type.
   It marks one (or several) elements for refinement.

   mark [$h | {[<rule> [<side>]] [$a | $i <Id> | $s]} | $c]

   .  <rule>                 - specify a refinement rule ("red" is default)
   .  <side>                 - has to be specified if the corresponding rule can be applied in several orientations
   .  $a                     - refine all (leave) elements
   .  $i <Id>                - refine the element with <Id>
   .  $s                     - refine all elements from the current selection
   .  $h                     - show available rules

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MarkCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  ELEMENT *theElement;
  char rulename[32];
  INT i,j,l,mode,rv,Rule;
  long nmarked;

  /* following variables: keep type for sscanf */
  int id,Side;


  /* first check help option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='h')
    {
      UserWrite("the following rules are available:\n");
      for (j=0; j<NO_OF_RULES; j++)
      {
        UserWrite(myMR[j].RuleName);
        UserWrite("\n");
      }
      return (OKCODE);
    }

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"mark","no open multigrid");
    return (CMDERRORCODE);
  }

  /* scan parameters */
  rv = sscanf(argv[0],"mark %31[redbluecopycoarsnoi_123qt] %d",rulename,&Side);
  if (rv<1)
  {
    /* set the default rule */
    strcpy(rulename,myMR[0].RuleName);
    Rule = myMR[0].RuleId;
    Side = NO_SIDE_SPECIFIED;
  }
  else
  {
    Rule = NO_RULE_SPECIFIED;
    for (i=0; i<NO_OF_RULES; i++)
      if (strcmp(rulename,myMR[i].RuleName)==0)
      {
        Rule = myMR[i].RuleId;
        break;
      }

    if (Rule==NO_RULE_SPECIFIED)
    {
      sprintf(buffer,"unknown rule '%s'",rulename);
      PrintErrorMessage('E',"mark",buffer);
      return (PARAMERRORCODE);
    }

    if (rv!=2)
      Side = NO_SIDE_SPECIFIED;
  }

  /* check options */
  mode = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"mark","specify only one option of a, i, s");
        return (PARAMERRORCODE);
      }
      mode = MARK_ALL;
      break;

    case 'i' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"mark","specify only one option of a, i, s");
        return (PARAMERRORCODE);
      }
      mode = MARK_ID;
      if (sscanf(argv[i],"i %d ", &id)!=1)
      {
        PrintErrorMessage('E',"mark","cannot scan id");
        return (PARAMERRORCODE);
      }
      break;

    case 's' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"mark","specify only one option of a, i, s");
        return (PARAMERRORCODE);
      }
      mode = MARK_SELECTION;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("mark",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (mode==FALSE)
  {
    PrintErrorMessage('E',"mark","specify exactly one option of a, i, s");
    return (PARAMERRORCODE);
  }

  /* print rule and side */
  if (Side==NO_SIDE_SPECIFIED)
    sprintf(buffer,"   using rule %s (no side given)\n",rulename);
  else
    sprintf(buffer,"   using rule %s, side %d\n",rulename,(int)Side);
  UserWrite(buffer);

  nmarked = 0;
  switch (mode)
  {
  case MARK_ALL :
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l)); theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
        {
          if ((rv = MarkForRefinement(theElement,Rule,Side))!=0)
          {
            l = TOPLEVEL(theMG);
            break;
          }
          else
            nmarked++;
        }
    break;

  case MARK_ID :
    for (l=0; l<=TOPLEVEL(theMG); l++)
      if ((theElement=FindElementFromId(GRID_ON_LEVEL(theMG,l),id))!=NULL)
        break;
      else
        nmarked++;

    if (theElement == NULL)
    {
      PrintErrorMessage('E',"mark","no element with this id");
      return (PARAMERRORCODE);
    }

    if (EstimateHere(theElement))
      rv = MarkForRefinement(theElement,Rule,Side);
    break;

  case MARK_SELECTION :
    if (SELECTIONMODE(theMG)==elementSelection)
      for (i=0; i<SELECTIONSIZE(theMG); i++)
      {
        theElement = (ELEMENT *)SELECTIONOBJECT(theMG,i);
        if (EstimateHere(theElement))
        {
          if ((rv = MarkForRefinement(theElement,Rule,Side))!=0)
            break;
          else
            nmarked++;
        }
      }
    break;
  }

  sprintf(buffer," %ld elements marked for refinement\n",nmarked);
  UserWrite(buffer);

  if (rv)
  {
    sprintf(buffer,"rule could not be applied for element with ID %ld, nothing marked",ID(theElement));
    PrintErrorMessage('W',"mark",buffer);
    return (CMDERRORCODE);
  }
  else
    return(OKCODE);
}

/****************************************************************************/
/*D
   smooth - invoke hierarchical multigrid smoother

   DESCRIPTION:
   This command invokes hierarchical multigrid smoother,
   calling the function 'SmoothMultiGrid'.

   'smooth <nIt> [$b]'

   .    <nIt>   - number of iterations
   .    $b      - also smooth boundary nodes
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SmoothMGCommand - Invoke hierarchical multigrid smoother

   SYNOPSIS:
   static INT SmoothMGCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function invokes hierarchical multigrid smoother.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SmoothMGCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,smbdry;

  /* following variables: keep type for sscanf */
  int nit;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"smooth","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check pararmeters */
  if (sscanf(argv[0],"smooth %d",&nit)!=1)
  {
    PrintHelp("smooth",HELPITEM," (specify number of iterations)");
    return (PARAMERRORCODE);
  }

  /* check options */
  smbdry = GM_KEEP_BOUNDARY_NODES;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'b' :
      smbdry = GM_MOVE_BOUNDARY_NODES;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("move",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (SmoothMultiGrid(theMG,nit,smbdry)!=GM_OK)
  {
    PrintErrorMessage('E',"smooth","failed smoothing the multigrid");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   ordernodes - order the nodes lexicographically according to the specified directions


   DESCRIPTION:
   This command orders the nodes according to the user provided dependencies.
   It orders the nodes of the current multigrid, calling the function
   'OrderNodesInGrid' on all levels.

   If specified the links are ordered in the corresponding order.

   'ordernodes ur|ul|dr|dl|ru|rd|lu|ld' [$l <level>] [$L]

   . [$l~<level>] - only on level <level>
   .n      u=up, d=down, r=right, l=left

    EXAMPLE:
        'ordernodes rd $l2'

        Order nodes of grid level 2 lexicographically in horizontal lines from
    left to right and the lines vertical from top down.
   D*/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Function:  OrderNodesCommand                                             */
/*                                                                          */
/* Purpose:   reorder nodes in lexicographical order                        */
/*                                                                          */
/* Input:     INT argc: number of arguments (incl. its own name)            */
/*            char **argv: array of strings giving the arguments            */
/*                                                                          */
/* Output:    INT return code see header file                               */
/*                                                                          */
/****************************************************************************/

static INT OrderNodesCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  INT i,res,level,fromLevel,toLevel;
  INT sign[DIM],order[DIM],xused,yused,zused,error,AlsoOrderLinks;
  char ord[3];

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"ordernodes","no open multigrid");
    return (CMDERRORCODE);
  }
  fromLevel = 0;
  toLevel   = TOPLEVEL(theMG);

  /* read ordering directions */
        #ifdef __TWODIM__
  res = sscanf(argv[0],expandfmt("ordernodes %2[rlud]"),ord);
        #else
  res = sscanf(argv[0],expandfmt("ordernodes %3[rlbfud]"),ord);
        #endif
  if (res!=1)
  {
    PrintHelp("ordernodes",HELPITEM," (could not read order type)");
    return(PARAMERRORCODE);
  }
  if (strlen(ord)!=DIM)
  {
    PrintHelp("ordernodes",HELPITEM," (specify DIM chars out of 'rlud' or 'rlbfud' resp.)");
    return(PARAMERRORCODE);
  }
  error = xused = yused = zused = FALSE;
  for (i=0; i<DIM; i++)
    switch (ord[i])
    {
    case 'r' :
      if (xused) error = TRUE;
      xused = TRUE;
      order[i] = _X_; sign[i] =  1; break;
    case 'l' :
      if (xused) error = TRUE;
      xused = TRUE;
      order[i] = _X_; sign[i] = -1; break;

                        #ifdef __TWODIM__
    case 'u' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] =  1; break;
    case 'd' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] = -1; break;
                        #else
    case 'b' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] =  1; break;
    case 'f' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] = -1; break;

    case 'u' :
      if (zused) error = TRUE;
      zused = TRUE;
      order[i] = _Z_; sign[i] =  1; break;
    case 'd' :
      if (zused) error = TRUE;
      zused = TRUE;
      order[i] = _Z_; sign[i] = -1; break;
                        #endif
    }
  if (error)
  {
    PrintHelp("ordernodes",HELPITEM," (bad combination of 'rludr' or 'rlbfud' resp.)");
    return(PARAMERRORCODE);
  }

  /* check options */
  AlsoOrderLinks = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'l' :
      res = sscanf(argv[i],"l %d",&level);
      if (res!=1)
      {
        PrintErrorMessage('E',"ordernodes","could not read level");
        return(PARAMERRORCODE);
      }
      if ((level>=fromLevel)&&(level<=toLevel))
        fromLevel = toLevel = level;
      else
      {
        PrintErrorMessage('E',"ordernodes","level out of range");
        return(PARAMERRORCODE);
      }
      break;

    case 'L' :
      AlsoOrderLinks = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("ordernodes",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* first we renumber the multigrid (to have node-IDs coinciding with lists) */
  if (RenumberMultiGrid(theMG)!=GM_OK)
  {
    PrintErrorMessage('E',"ordernodes","renumbering of the mg failed");
    return (CMDERRORCODE);
  }

  /* now we reorder the nodes on the specified levels */
  for (level=fromLevel; level<=toLevel; level++)
  {
    theGrid = GRID_ON_LEVEL(theMG,level);

    sprintf(buffer," [%d:",level);
    UserWrite(buffer);

    if (OrderNodesInGrid(theGrid,order,sign,AlsoOrderLinks)!=GM_OK)
    {
      PrintErrorMessage('E',"ordernodes","OrderNodesInGrid failed");
      return (CMDERRORCODE);
    }
    UserWrite("o]");
  }

  UserWrite("\n");

  return (OKCODE);
}


/****************************************************************************/
/*D
   lexorderv - order the vectors lexicographically


   DESCRIPTION:
   This command orders the vectors lexicographically according to the user
   specified directions.
   It orders the vectors of the current multigrid, calling the function
   'LexOrderVectorsInGrid'.

   'lexorderv '

   .  $m FFCCLL | FCFCLL     - possible modes are FFCCLL or FCFCLL
   .  $d <dep-proc>          - the ordering algorithm uses this dependency procedure...
   .  $o <dep-proc options>  - ...and passes these options to it
   .  $a                     - order all levels of the current multigrid
   D*/
/****************************************************************************/

static INT LexOrderVectorsCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  INT i,res,level,fromLevel,toLevel;
  INT sign[DIM],order[DIM],which,xused,yused,zused,error,AlsoOrderMatrices,SpecialTreatSkipVecs;
  char ord[3];

  theMG = GetCurrentMultigrid();
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"lexorderv","no open multigrid");
    return (CMDERRORCODE);
  }
  fromLevel = 0;
  toLevel   = TOPLEVEL(theMG);

  /* read ordering directions */
        #ifdef __TWODIM__
  res = sscanf(argv[0],expandfmt("lexorderv %2[rlud]"),ord);
        #else
  res = sscanf(argv[0],expandfmt("lexorderv %3[rlbfud]"),ord);
        #endif
  if (res!=1)
  {
    PrintHelp("lexorderv",HELPITEM," (could not read order type)");
    return(PARAMERRORCODE);
  }
  if (strlen(ord)!=DIM)
  {
    PrintHelp("lexorderv",HELPITEM," (specify DIM chars out of 'rlud' or 'rlbfud' resp.)");
    return(PARAMERRORCODE);
  }
  error = xused = yused = zused = FALSE;
  for (i=0; i<DIM; i++)
    switch (ord[i])
    {
    case 'r' :
      if (xused) error = TRUE;
      xused = TRUE;
      order[i] = _X_; sign[i] =  1; break;
    case 'l' :
      if (xused) error = TRUE;
      xused = TRUE;
      order[i] = _X_; sign[i] = -1; break;

                        #ifdef __TWODIM__
    case 'u' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] =  1; break;
    case 'd' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] = -1; break;
                        #else
    case 'b' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] =  1; break;
    case 'f' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] = -1; break;

    case 'u' :
      if (zused) error = TRUE;
      zused = TRUE;
      order[i] = _Z_; sign[i] =  1; break;
    case 'd' :
      if (zused) error = TRUE;
      zused = TRUE;
      order[i] = _Z_; sign[i] = -1; break;
                        #endif
    }
  if (error)
  {
    PrintHelp("lexorderv",HELPITEM," (bad combination of 'rludr' or 'rlbfud' resp.)");
    return(PARAMERRORCODE);
  }

  /* check options */
  AlsoOrderMatrices = SpecialTreatSkipVecs = FALSE;
  /* which = GM_TAKE_SKIP | GM_TAKE_NONSKIP;*/
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'l' :
      res = sscanf(argv[i],"l %d",&level);
      if (res!=1)
      {
        PrintErrorMessage('E',"lexorderv","could not read level");
        return(PARAMERRORCODE);
      }
      if ((level>=fromLevel)&&(level<=toLevel))
        fromLevel = toLevel = level;
      else
      {
        PrintErrorMessage('E',"lexorderv","level out of range");
        return(PARAMERRORCODE);
      }
      break;

    case 'm' :
      AlsoOrderMatrices = TRUE;
      break;

    /*			case 'w':
                                    which = 0;
                                    if (strchr(argv[i],'s')!=NULL)
                                            which |= GM_TAKE_SKIP;
                                    if (strchr(argv[i],'n')!=NULL)
                                            which |= GM_TAKE_NONSKIP;
                                    break; */

    case 's' :
      if              (strchr(argv[i],'<')!=NULL)
        SpecialTreatSkipVecs = GM_PUT_AT_BEGIN;
      else if (strchr(argv[i],'>')!=NULL)
        SpecialTreatSkipVecs = GM_PUT_AT_END;
      else if (strchr(argv[i],'0')!=NULL)
        SpecialTreatSkipVecs = FALSE;
      else
      {
        PrintErrorMessage('E',"lexorderv","use < or > with s-option");
        return(PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("lexorderv",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* now we reorder the vectors on the specified levels */
  for (level=fromLevel; level<=toLevel; level++)
  {
    theGrid = GRID_ON_LEVEL(theMG,level);

    sprintf(buffer," [%d:",level);
    UserWrite(buffer);

    if (LexOrderVectorsInGrid(theGrid,order,sign,SpecialTreatSkipVecs,AlsoOrderMatrices)!=GM_OK)
    {
      PrintErrorMessage('E',"lexorderv","LexOrderVectorsInGrid failed");
      return (CMDERRORCODE);
    }
    UserWrite("ov]");
  }

  UserWrite("\n");

  return (OKCODE);
}

/****************************************************************************/
/*D
   shellorderv - order the vectors shell by shell


   DESCRIPTION:
   This command orders the vectors of the current level of the current
   multigrid in shells starting from a seed.

   'shellorderv f | l | s'

   .  f - take first vector as seed
   .  l - take last vector as seed
   .  s - take selected vector as seed

   D*/
/****************************************************************************/

static INT ShellOrderVectorsCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  VECTOR *seed;
  char option;

  NO_OPTION_CHECK(argc,argv);

  theMG = GetCurrentMultigrid();
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"shellorderv","no open multigrid");
    return (CMDERRORCODE);
  }
  theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));

  if (sscanf(argv[0],"shellorderv %c",&option)!=1)
  {
    PrintErrorMessage('E',"shellorderv","specify f, l or s");
    return (CMDERRORCODE);
  }

  switch (option)
  {
  case 'f' :
    seed = FIRSTVECTOR(theGrid);
    break;
  case 'l' :
    seed = LASTVECTOR(theGrid);
    break;
  case 's' :
    if (SELECTIONMODE(theMG)!=vectorSelection)
    {
      PrintErrorMessage('E',"shellorderv","no vector selection");
      return (CMDERRORCODE);
    }
    if (SELECTIONSIZE(theMG)!=1)
    {
      PrintErrorMessage('E',"shellorderv","select ONE vector");
      return (CMDERRORCODE);
    }
    seed = (VECTOR *)SELECTIONOBJECT(theMG,0);
    break;
  default :
    PrintErrorMessage('E',"shellorderv","specify f, l or s");
    return (CMDERRORCODE);
  }

  if (ShellOrderVectors(theGrid,seed)!=0)
  {
    PrintErrorMessage('E',"shellorderv","ShellOrderVectors failed");
    return (CMDERRORCODE);
  }
  else
    return (OKCODE);
}

/****************************************************************************/
/*D
   orderv - order the vectors according to the user provided dependencies


   DESCRIPTION:
   This command orders the vectors according to the user provided dependencies.
   It orders the vectors of the current multigrid, calling the function
   'OrderVectors'.

   'orderv $m FFCCLL | FCFCLL $d <dep-proc> $o <dep-proc options> $c <find-cut-proc> [$a]'

   .  $m FFCCLL | FCFCLL     - possible modes are FFCCLL or FCFCLL
   .  $d <dep-proc>          - the ordering algorithm uses this dependency procedure...
   .  $o <dep-proc options>  - ...and passes these options to it
   .  $a                     - order all levels of the current multigrid
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   OrderVectorsCommand - Order the vectors according to the user provided dependencies

   SYNOPSIS:
   static INT OrderVectorsCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function orders the vectors according to the user provided dependencies.
   It orders the vectors of the current multigrid.

   orderv $m FFCCLL | FCFCLL $d <dep-proc> $o dep-proc options> [$a]

   .  $m FFCCLL | FCFCLL     - possible modes are FFCCLL or FCFCLL
   .  $d <dep-proc>          - the ordering algorithm uses this dependency procedure...
   .  $o <dep-proc options>  - ...and passes these options to it
   .  $a                     - order all levels of the current multigrid

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT OrderVectorsCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char modestr[7];              /* <-- with change of size also change format string in sscanf! */
  char *dep,*dep_opt,*cut;
  INT mode,i,levels,PutSkipFirst,SkipPat;
  int iValue;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"orderv","no open multigrid");
    return (CMDERRORCODE);
  }

  levels           = GM_CURRENT_LEVEL;
  mode             = FALSE;
  dep              = dep_opt = cut = NULL;
  PutSkipFirst = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'm' :
      if (sscanf(argv[i],"m %6[FCL]",modestr)!=1)
      {
        PrintHelp("orderv",HELPITEM," (could not read the mode)");
        return (PARAMERRORCODE);
      }
      if (strcmp(modestr,"FCFCLL")==0)
        mode = GM_FCFCLL;
      else if (strcmp(modestr,"FFCCLL")==0)
        mode = GM_FFCCLL;
      else if (strcmp(modestr,"FFLCLC")==0)
        mode = GM_FFLCLC;
      else
      {
        PrintHelp("orderv",HELPITEM," (you have to specify FFCCLL or FCFCLL as mode)");
        return (PARAMERRORCODE);
      }
      break;

    case 'd' :
      /* skip leading blanks */
      dep = argv[i]+1;
      while ((*dep!='\0') && (strchr(WHITESPACE,*dep)!=NULL))
        dep++;
      break;

    case 'o' :
      /* skip leading blanks */
      dep_opt = argv[i]+1;
      while ((*dep_opt!='\0') && (strchr(WHITESPACE,*dep_opt)!=NULL))
        dep_opt++;
      break;

    case 'c' :
      /* skip leading blanks */
      cut = argv[i]+1;
      while ((*cut!='\0') && (strchr(WHITESPACE,*cut)!=NULL))
        cut++;
      break;

    case 's' :
      PutSkipFirst = TRUE;
      if (sscanf(argv[i],"s %x",&iValue)!=1)
      {
        PrintErrorMessage('E',"orderv","could not read skip pattern");
        return(PARAMERRORCODE);
      }
      SkipPat = iValue;
      break;

    case 'a' :
      levels = GM_ALL_LEVELS;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("orderv",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (mode==FALSE)
  {
    PrintErrorMessage('E',"orderv","the m option is mandatory");
    return (PARAMERRORCODE);
  }

  if (dep==NULL)
  {
    PrintErrorMessage('E',"orderv","the d option is mandatory");
    return (PARAMERRORCODE);
  }

  if (dep_opt==NULL)
  {
    PrintErrorMessage('E',"orderv","the o option is mandatory");
    return (PARAMERRORCODE);
  }

  if (OrderVectors(theMG,levels,mode,PutSkipFirst,SkipPat,dep,dep_opt,cut)!=GM_OK)
  {
    PrintErrorMessage('E',"orderv","order vectors failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   setindex - set vector index in ascending order

   DESCRIPTION:
   'setindex' sets the vector index in ascending order.
   D*/
/****************************************************************************/

static INT SetIndexCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"orderv","no open multigrid");
    return (CMDERRORCODE);
  }
  theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));

  if (l_setindex(theGrid)!=0)
  {
    PrintErrorMessage('E',"setindex","l_setindex failed");
    return (CMDERRORCODE);
  }
  else
    return (OKCODE);
}

/****************************************************************************/
/*D
   find - find (and select) a node (element) from a given position

   DESCRIPTION:
   This function finds (and selects) a node (element) from a given position,
   where some tolerance can be specified.
   It finds a node (element) on the current level of the current multigrid,
   using the functions
   'FindNodeFromPosition', 'FindElementFromPosition'
   'AddNodeToSelection', 'RemoveNodeFromSelection',
   'AddElementToSelection' and 'RemoveElementFromSelection'

   'find <x> <y> <z> {$n <tol> | $e} [$s]'

   .  <x>~<y>~<z>            - specify as much coordinates as the space has dimensions
   .  $n~<tol>               - find a node maching the position with tolerance <tol>
   .  $e                     - find an element maching the position
   .  $s                     - add the selected node (element) to the selection buffer
   .n                           (if not specified the node is just listed)
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   FindCommand - Find (& select) a node (element) from a given position (+tol)

   SYNOPSIS:
   static INT FindCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function finds (& selects) a node (element) from a given position (+tol).
   It finds a node (element) on the current level of the current multigrid.

   find <x> <y> <z> {$n <tol> | $v <tol> | $e} [$s]

   .  <x> <y> <z>            - specify as much coordinates as the space has dimensions
   .  $n <tol>               - find a node maching the position with tolerance <tol>
   .  $v <tol>               - find a vector maching the position with tolerance <tol>
   .  $e                     - find an element maching the position
   .  $s                     - add the selected node (element) to the selection buffer
                          - (if not specified the node is just listed)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT FindCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  NODE *theNode;
  VECTOR *theVector;
  ELEMENT *theElement;
  COORD xc[DIM],tolc[DIM];
  INT i,j,select,isNode,isElement,isVector;

  /* following variables: keep type for sscanf */
  float x[DIM_MAX],tol;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"find","no open multigrid");
    return (CMDERRORCODE);
  }
  theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));

  /* check pararmeters */
  if (sscanf(argv[0],"find %f %f %f",x,x+1,x+2)!=DIM)
  {
    PrintHelp("find",HELPITEM," (could not get coordinates)");
    return (PARAMERRORCODE);
  }
  for (i=0; i<DIM; i++)
    xc[i] = x[i];

  /* check options */
  select = isNode = isElement = isVector = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'n' :
      if (sscanf(argv[i],"n %f",&tol)!=1)
      {
        PrintHelp("find",HELPITEM," (could not read tolerance)");
        return (PARAMERRORCODE);
      }
      for (j=0; j<DIM; j++)
        tolc[j] = tol;
      theNode = FindNodeFromPosition(theGrid,xc,tolc);
      if (theNode==NULL)
      {
        PrintErrorMessage('W',"find","no node is matching");
        return (OKCODE);
      }
      isNode = TRUE;
      break;

    case 'v' :
      if (sscanf(argv[i],"v %f",&tol)!=1)
      {
        PrintHelp("find",HELPITEM," (could not read tolerance)");
        return (PARAMERRORCODE);
      }
      for (j=0; j<DIM; j++)
        tolc[j] = tol;
      theVector = FindVectorFromPosition(theGrid,xc,tolc);
      if (theVector==NULL)
      {
        PrintErrorMessage('W',"find","no vector is matching");
        return (OKCODE);
      }
      isVector = TRUE;
      break;

    case 'e' :
      theElement = FindElementFromPosition(theGrid,xc);
      if (theElement==NULL)
      {
        PrintErrorMessage('W',"find","no element is matching");
        return (OKCODE);
      }
      isElement = TRUE;
      break;

    case 's' :
      select = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("find",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (select)
  {
    if (isNode)
      if (AddNodeToSelection(theMG,theNode)!=GM_OK)
      {
        PrintErrorMessage('E',"find","selecting the node failed");
        return (CMDERRORCODE);
      }

    if (isVector)
      if (AddVectorToSelection(theMG,theVector)!=GM_OK)
      {
        PrintErrorMessage('E',"find","selecting the vector failed");
        return (CMDERRORCODE);
      }

    if (isElement)
      if (AddElementToSelection(theMG,theElement)!=GM_OK)
      {
        PrintErrorMessage('E',"find","selecting the element failed");
        return (CMDERRORCODE);
      }
  }
  else
  {
    if (isNode)
      ListNode(theMG,theNode,FALSE,FALSE,FALSE,FALSE);

    if (isVector)
      ListVector(theMG,theVector,FALSE,FALSE);

    if (isElement)
      ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   select - select a node or element from a given position

   DESCRIPTION:
   This function finds (and selects) a node (element) from a given position,
   where some tolerance can be specified.
   It adds/removes nodes/elements from the selection buffer of the current
   multigrid, using the functions
   'FindNodeFromId', 'FindElementFromId'
   'AddNodeToSelection', 'RemoveNodeFromSelection',
   'AddElementToSelection' and 'RemoveElementFromSelection'

   'select $c | $n {+|-} <Id> | $e {+|-} <Id>'

   .  $c                     - clear the selection buffer
   .  $n {+|-} <Id>          - add (+) or remove (-) the node with <Id> to (from) the selection buffer
   .  $e {+|-} <Id>          - add (+) or remove (-) the element with <Id> to (from) the selection buffer
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SelectCommand - Find (and selct) a node (element) from a given position (+tol

   SYNOPSIS:
   static INT SelectCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function finds (& selects) a node (element) from a given position (+tol).
   It adds/removes nodes/elements from the selction buffer of the current multigrid.

   select $c | $n {+|-} <Id> | $e {+|-} <Id>

   .  $c                     - clear the selection buffer
   .  $n~{+|-}~<Id>          - add (+) or remove (-) the node with <Id> to (from) the selection buffer
   .  $e~{+|-}~<Id>          - add (+) or remove (-) the element with <Id> to (from) the selection buffer

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SelectCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  NODE *theNode;
  ELEMENT *theElement;
  INT i,level;
  char c;

  /* following variables: keep type for sscanf */
  int id;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"select","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'c' :
      ClearSelection(theMG);
      break;

    case 'n' :
      if (sscanf(argv[i],"n %c %d",&c,&id)!=2)
      {
        PrintErrorMessage('E',"select","could not get +/- or ID");
        return (PARAMERRORCODE);
      }
      if (c=='+')
      {
        /* search node */
        theNode = NULL;
        for (level=0; level<=TOPLEVEL(theMG); level++)
          if ((theNode=FindNodeFromId(GRID_ON_LEVEL(theMG,level),id))!=NULL)
            break;
        if (theNode==NULL)
        {
          sprintf(buffer,"node with ID %ld not found",(long)id);
          PrintErrorMessage('E',"select",buffer);
          return (CMDERRORCODE);
        }
        if (AddNodeToSelection(theMG,theNode)!=GM_OK)
        {
          PrintErrorMessage('E',"select","selecting the node failed");
          return (CMDERRORCODE);
        }
      }
      else if (c=='-')
      {
        if (SELECTIONMODE(theMG)==nodeSelection)
          for (i=0; i<SELECTIONSIZE(theMG); i++)
          {
            theNode = (NODE *)SELECTIONOBJECT(theMG,i);
            if (ID(theNode)==id)
              break;
          }
        if (theNode==NULL)
        {
          sprintf(buffer,"node with ID %ld is not in selection",(long)id);
          PrintErrorMessage('E',"select",buffer);
          return (CMDERRORCODE);
        }
        if (RemoveNodeFromSelection(theMG,theNode)!=GM_OK)
        {
          PrintErrorMessage('E',"select","removing the node failed");
          return (CMDERRORCODE);
        }
      }
      else
      {
        PrintErrorMessage('E',"select","specify + or - with n option");
        return (PARAMERRORCODE);
      }
      break;

    case 'e' :
      if (sscanf(argv[i],"e %c %d",&c,&id)!=2)
      {
        PrintErrorMessage('E',"select","could not get +/- or ID");
        return (PARAMERRORCODE);
      }
      if (c=='+')
      {
        /* search element */
        theElement = NULL;
        for (level=0; level<=TOPLEVEL(theMG); level++)
          if ((theElement=FindElementFromId(GRID_ON_LEVEL(theMG,level),id))!=NULL)
            break;
        if (theElement==NULL)
        {
          sprintf(buffer,"element with ID %ld not found",(long)id);
          PrintErrorMessage('E',"select",buffer);
          return (CMDERRORCODE);
        }
        if (AddElementToSelection(theMG,theElement)!=GM_OK)
        {
          PrintErrorMessage('E',"select","selecting the element failed");
          return (CMDERRORCODE);
        }
      }
      else if (c=='-')
      {
        if (SELECTIONMODE(theMG)==elementSelection)
          for (i=0; i<SELECTIONSIZE(theMG); i++)
          {
            theElement = (ELEMENT *)SELECTIONOBJECT(theMG,i);
            if (ID(theElement)==id)
              break;
          }
        if (theElement==NULL)
        {
          sprintf(buffer,"element with ID %ld is not in selection",(long)id);
          PrintErrorMessage('E',"select",buffer);
          return (CMDERRORCODE);
        }
        if (RemoveElementFromSelection(theMG,theElement)!=GM_OK)
        {
          PrintErrorMessage('E',"select","removing the element failed");
          return (CMDERRORCODE);
        }
      }
      else
      {
        PrintErrorMessage('E',"select","specify + or - with n option");
        return (PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("select",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  return (OKCODE);
}

/****************************************************************************/
/*D
   extracon - display or delete extra connections

   DESCRIPTION:
   This command displays or deleteâ extra connections.

   'extracon [$d]'

   .  $c - also check the connections
   D*/
/****************************************************************************/

static INT ExtraConnectionCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  VECTOR *vec;
  MATRIX *mat;
  INT delete,i,nextra;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"extracon","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  delete = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'd' :
      delete = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("extracon",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }
  theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));

  /* count extra connections on current level */
  nextra = 0;
  for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
    for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
      if (CEXTRA(MMYCON(mat)))
        nextra++;
  nextra /= 2;                  /* have been counted twice */

  UserWriteF("%d extra connections on level %d (total %d)\n",(int)nextra,(int)CURRENTLEVEL(theMG),(int)NC(theGrid));

  SetStringValue(":extraconratio",nextra/((DOUBLE)NC(theGrid)));

  if (delete)
  {
    if (DisposeExtraConnections(theGrid)!=GM_OK)
    {
      PrintErrorMessage('E',"extracon","deleting extra connections failed");
      return (CMDERRORCODE);
    }
    UserWrite("...deleted\n");
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   check - check consistency of the data structure

   DESCRIPTION:
   This command checks consistency of the data structure, using
   the functions 'CheckGrid' and 'CheckConnections'.

   'check [$c]'

   .  $c - also check the connections
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   CheckCommand - Check consistency of the data structure

   SYNOPSIS:
   static INT CheckCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function checks consistency of the data structure.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT CheckCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  INT checkconn,level,err,i;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"check","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  checkconn = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'c' :
      checkconn = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("check",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  err = 0;
  for (level=0; level<=TOPLEVEL(theMG); level++)
  {
    theGrid = GRID_ON_LEVEL(theMG,level);
    sprintf(buffer,"[%d:",(int)level);
    UserWrite(buffer);
    if (CheckGrid(theGrid)!=GM_OK)
    {
      err++;
      UserWrite(" grid bad");
    }
    else if ( (checkconn) && (CheckConnections(theGrid)!=GM_OK))
    {
      err++;
      UserWrite(" connections bad");
    }
    else
      UserWrite("ok");

    UserWrite("] ");
  }
  UserWrite("\n");

  if (err)
    return (CMDERRORCODE);
  else
    return (OKCODE);
}

/****************************************************************************/
/*D
   quality - calculate minimal and maximal angle of specified elements

   DESCRIPTION:
   This command calculates minimal and maximal angle of the specified elements
   and lists elements with angle < or > given angles.
   It calls the functions 'QualityElement'.

   'quality $a | $s | {$i <fromID> [<toID>]} [$< <angle>] [$> <angle>]'

   .    $a                     - check angles of all elements in the multigrid
   .    $s                     - check angles of the selected elements
   .    $i                     - check angles of elements with an ID in the range <fromID> through <toID>
   .n                           if <fromID> is omitted only the element with <fromID> is listed

   .    $<~<angle>             - print info for all elements the minangle of which is < <angle>
   .    $>~<angle>             - print info for all elements the maxangle of which is > <angle>

     (angles in degree 0-360)
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   QualityElement - calculate minimal and maximal angle of an element

   SYNOPSIS:
   static void QualityElement (MULTIGRID *theMG, ELEMENT *theElement);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  theElement - pointer to Element

   DESCRIPTION:
   This function calculates the minimal and the maximal angle of elements
   and lists elements with angle < or > given angles.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

static void QualityElement (MULTIGRID *theMG, ELEMENT *theElement)
{
  min = 360.0;
  max = 0.0;

  MinMaxAngle(theElement,&min,&max);
  minangle = MIN(min,minangle);
  maxangle = MAX(max,maxangle);
  if ((lessopt && (min<themin)) && (greateropt && (max>themax)))
  {
    UserWrite(minmaxtext);
    ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
  }
  else if (lessopt && (min<themin))
  {
    UserWrite(mintext);
    ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
  }
  else if (lessopt && (min<themin))
  {
    UserWrite(mintext);
    ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
  }
  return;
}

/****************************************************************************/
/*
   QualityCommand - Calculate min and max angle of elements

   SYNOPSIS:
   static INT QualityCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function calculates min and max angle of elements
   and lists elements with angle < or > given angles.
   It checks the min and max angles of the elements in the grid.

   quality $a | $s | {$i <fromID> [<toID>]} [$< <angle>] [$> <angle>]

   .  $a                     - check angles of all elements in the multigrid
   .  $s                     - check angles of the selected elements
   .  $i                     - check angles of elements with an ID in the range <fromID> through <toID>
   .n                        if <fromID> is omitted only the element with <fromID> is listed

   .  $< <angle>             - print info for all elements the minangle of which is < <angle>
   .  $> <angle>             - print info for all elements the maxangle of which is > <angle>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT QualityCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  ELEMENT *theElement;
  INT i,fromE,toE,res,mode;

  /* following variables: keep type for sscanf */
  float angle;
  long f,t;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"quality","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  lessopt = greateropt = mode = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"quality","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ALL;
      break;

    case 'i' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"quality","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_ID;
      res = sscanf(argv[i]," i %ld %ld",&f,&t);
      fromE = f;
      toE   = t;
      if (res<1)
      {
        PrintErrorMessage('E',"quality","specify at least one id with the i option");
        return (PARAMERRORCODE);
      }
      else if (res==1)
        toE = fromE;
      else if (fromE>toE)
      {
        PrintErrorMessage('E',"quality","from ID > to ID");
        return (PARAMERRORCODE);
      }
      break;

    case '<' :
      lessopt = TRUE;
      if (sscanf(argv[i],"< %f",&angle)!=1)
      {
        PrintErrorMessage('E',"quality","could not get angle of < option");
        return (CMDERRORCODE);
      }
      themin = angle;
      break;

    case '>' :
      greateropt = TRUE;
      if (sscanf(argv[i],"> %f",&angle)!=1)
      {
        PrintErrorMessage('E',"quality","could not get angle of > option");
        return (CMDERRORCODE);
      }
      themax = angle;
      break;

    case 's' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"quality","specify either the a, s or i option");
        return (PARAMERRORCODE);
      }
      mode = DO_SELECTION;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("quality",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  sprintf(mintext," < %g: ",(float)themin);
  sprintf(maxtext," > %g: ",(float)themax);
  sprintf(minmaxtext," < %g and > %g: ",(float)themin,(float)themax);

  minangle =      MAX_D;
  maxangle = -MAX_D;

  switch (mode)
  {
  case DO_ID :
    for (theGrid=GRID_ON_LEVEL(theMG,0); theGrid!=NULL; theGrid=UPGRID(theGrid))
      for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
        if ((ID(theElement)>=fromE) && (ID(theElement)<=toE))
          QualityElement(theMG,theElement);
    break;

  case DO_ALL :
    for (theGrid=GRID_ON_LEVEL(theMG,0); theGrid!=NULL; theGrid=UPGRID(theGrid))
      for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
        QualityElement(theMG,theElement);
    break;

  case DO_SELECTION :
    if (SELECTIONMODE(theMG)==elementSelection)
      for (i=0; i<SELECTIONSIZE(theMG); i++)
        QualityElement(theMG,(ELEMENT *)SELECTIONOBJECT(theMG,i));
    break;

  default :
    PrintErrorMessage('E',"quality","specify one option of a, s or i");
    return (PARAMERRORCODE);
  }

  sprintf(buffer," min angle = %12.4f\n max angle = %12.4f\n",(float)minangle,(float)maxangle);
  UserWrite(buffer);

  return(OKCODE);
}

/****************************************************************************/
/*D
   bnodes - generate boundary nodes

   DESCRIPTION:
   This command generates boundary nodes.
   It reads the environ variables ':gg:RelRasterSize', ':gg:h_global',
   ':gg:searchconst', ':gg:angle', ':gg:epsi'.

   'bnodes'

   SEE ALSO:
   'makegrid'
   D*/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Function:  BnodesCommand                                                             */
/*                                                                          */
/* Purpose:   allocate a new document record with bndry node generation		*/
/*                                                                          */
/* Input:     INT argc: number of arguments (incl. its own name             */
/*            char **argv: array of strings giving the arguments            */
/*                                                                          */
/* Output:    INT return code see header file                               */
/*                                                                          */
/****************************************************************************/

#ifdef __TWODIM__

static INT BnodesCommand  (INT argc, char **argv)
{
  MULTIGRID *theMG;
  DOUBLE RelRasterSize,h_global;
  INT i,msizecoeffno;

  /* get current multigrid */
  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"bnodes","no open multigrid");
    return (CMDERRORCODE);
  }
  /* check options (no option up to now) */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("bnodes",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* read RelRasterSize */
  if(GetStringDOUBLEInRange(":gg:RelRasterSize",
                            10e-10, 10e+10, &RelRasterSize) != 0)
    return(PARAMERRORCODE);

  h_global = 1.0;
  if(GetStringINTInRange   (":gg:meshsizecoeffno",
                            -1     , 100    , &msizecoeffno) != 0)
    msizecoeffno = -1;

  if (msizecoeffno == -1)
    if(GetStringDOUBLEInRange(":gg:h_global",
                              10e-10, 10e+10, &h_global) != 0)
      return(PARAMERRORCODE);

  if (GenerateBnodes(theMG,RelRasterSize,h_global,msizecoeffno) != 0)
  {
    PrintErrorMessage('E',"bnodes","execution failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

#endif

/****************************************************************************/
/*D
   makegrid - generate grid

   DESCRIPTION:
   This command generates the grid. First, the command bnodes must be called.
   It reads the environ variables ':gg:RelRasterSize', ':gg:h_global',
   ':gg:searchconst', ':gg:angle', ':gg:epsi'.

   'makegrid ${W|w|K|k} [$E]'

   .  ${W|w|K|k} - W resp. K are using the accellerator,
   W resp. w use the angle criterion,
   K resp. k use the edge criterion
   .  $E - grid generator tries to create eqilateral triangles

   .n default: isosceles triangles

   EXAMPLE:
   .vb
   :gg:RelRasterSize	= 0.005;
   :gg:h_global		= 0.005;
   :gg:searchconst		= 0.002;
   :gg:angle			= 15;
   :gg:epsi			= 0.000625;

   bnodes;
   makegrid $E $W;
   .ve
   D*/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Function:  MakeGridCommand                                                   */
/*                                                                          */
/* Purpose:   create the advancing frontlists and generates the grid		*/
/*                                                                          */
/* Input:     INT argc: number of arguments (incl. its own name             */
/*            char **argv: array of strings giving the arguments            */
/*                                                                          */
/* Output:    INT return code see header file                               */
/*                                                                          */
/****************************************************************************/

#ifdef __TWODIM__

static INT MakeGridCommand  (INT argc, char **argv)
{
  long ElemID;
  MULTIGRID *theMG;
  DOUBLE angle;
  INT i;
  GG_ARG args;
  GG_PARAM params;

  /* get current multigrid */
  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"makegrid","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  args.doanimate = args.doupdate = args.dostep = args.equilateral = args.plotfront
                                                                      = args.printelem = args.doedge = args.doangle = args.doEdge = args.doAngle =  NO;
  ElemID = -1;


  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      args.doanimate = YES; break;

    case 'u' :
      args.doupdate = YES; break;

    case 's' :
      args.dostep = YES; break;

    case 'f' :
      args.plotfront = YES; break;

    case 'p' :
      args.printelem = YES; break;

    case 'E' :
      args.equilateral = YES; break;

    case 'e' :
      if (sscanf(argv[i],"e %ld",&ElemID)!=1)
      {
        PrintHelp("makegrid",HELPITEM," (could not read <element id>)");
        return (PARAMERRORCODE);
      }
      break;

    case 'k' :
      args.doedge = YES;
      break;

    case 'w' :
      args.doangle = YES;
      break;

    case 'K' :
      args.doEdge = YES;
      break;

    case 'W' :
      args.doAngle = YES;
      break;

    default :
      sprintf(buffer," (unknown option '%s')",argv[i]);
      PrintHelp("makegrid",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if(GetStringDOUBLEInRange(":gg:angle", 10e-10, 10e+10, &angle) != 0)
    return(PARAMERRORCODE);
  if(GetStringDOUBLEInRange(":gg:epsi", 10e-10, 10e+10, &(params.epsi))!= 0)
    return(PARAMERRORCODE);
  if(GetStringDOUBLEInRange(":gg:searchconst", 10e-10, 10e+10,
                            &(params.searchconst)) != 0)
    return(PARAMERRORCODE);
  if(GetStringINTInRange   (":gg:meshsizecoeffno",
                            -1      , 100   , &(params.msizecoeffno)) != 0)
    params.msizecoeffno = -1;

  if (params.msizecoeffno == -1)
    if(GetStringDOUBLEInRange(":gg:h_global",
                              10e-10, 10e+10, &(params.h_global)) != 0)
      return(PARAMERRORCODE);

  params.CheckCos = cos(angle*PI/180.0);

  if (GenerateGrid(theMG, &args, &params) != 0)
  {
    PrintErrorMessage('E',"makegrid","execution failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

#endif

/****************************************************************************/
/*D
   screensize - print the size of the monitor screen in pixels

   DESCRIPTION:
   This command prints the size of the monitor screen in pixels.
   It prints the size in pixels of the screen (if there) on the shell
   and in the variables ':screensize:width', ':screensize:height'

   'screensize'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ScreenSizeCommand - Print the size of the monitor screen in pixels

   SYNOPSIS:
   static INT ScreenSizeCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function prints the size of the monitor screen in pixels.
   It prints the size in pixels of the screen (if there) ==> max window size.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ScreenSizeCommand (INT argc, char **argv)
{
  INT size[2];

  NO_OPTION_CHECK(argc,argv);

  if (GetScreenSize(size)==FALSE)
  {
    PrintErrorMessage('W',"screensize","there is no monitor");
    return (OKCODE);
  }

  sprintf(buffer," screen width: %d, screen height: %d\n",(int)size[0], (int)size[1]);
  UserWrite(buffer);

  if (    (SetStringValue(":screensize:width",size[0])!=0)
          ||      (SetStringValue(":screensize:height",size[1])!=0))
  {
    PrintErrorMessage('E',"screensize","could not set :screensize:width or :screensize:height");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   openwindow - open a new window

   DESCRIPTION:
   This command opens an ug-window on an outputdevice
   (this will be the current window then).
   It calls the function 'CreateUGWindow'.


   'openwindow <h> <v> <dh> <dv> [$d <output device>] [$n <window name>]'

   .  <h>~<v>                - the lower left corner of the plotting region in the --> standardRefSys
   .  <dh>~<dv>              - the width and height resp. of the plotting region of the window
   .  $d~<output~device>     - specify the name of an output device (default: screen)
   .  $n~<window~name>       - optionally you can specify the window name
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   OpenWindowCommand - Open a new window

   SYNOPSIS:
   static INT OpenWindowCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function opens a new window.
   It opens an ug-window on an outputdevice (this will be the current window then).

   openwindow <h> <v> <dh> <dv> [$d <output device>] [$n <window name>]

   .  <h> <v>                - the lower left corner of the plotting region in the --> standardRefSys
   .  <dh> <dv>              - the width and height resp. of the plotting region of the window
   .  $d <output device>     - specify the name of an output device (default: screen)
   .  $n <window name>       - optionally you can specify the window name

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT OpenWindowCommand (INT argc, char **argv)
{
  OUTPUTDEVICE *theOutDev;
  UGWINDOW *theWin;
  char devname[NAMESIZE],winname[NAMESIZE];
  INT i;

  /* following variables: keep type for sscanf */
  int x,y,w,h;


  /* scan parameters */
  if (sscanf(argv[0],"openwindow %d %d %d %d",&x,&y,&w,&h)!=4)
  {
    PrintHelp("openwindow",HELPITEM," could not get all mandatory parameters");
    return (PARAMERRORCODE);
  }

  /* check options */
  theOutDev  = GetDefaultOutputDevice();
  winname[0] = '\0';
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'd' :
      if (sscanf(argv[i],expandfmt(CONCAT3("d %",NAMELENSTR,"[a-zA-Z0-9_-]")),devname)!=1)
      {
        PrintErrorMessage('E',"openwindow","specify device name with d option");
        return (PARAMERRORCODE);
      }
      if ((theOutDev=GetOutputDevice(devname))==NULL)
      {
        sprintf(buffer,"there is no device named '%s'",devname);
        PrintErrorMessage('E',"openwindow",buffer);
        return (PARAMERRORCODE);
      }
      break;

    case 'n' :
      if (sscanf(argv[i],expandfmt(CONCAT3("n %",NAMELENSTR,"[a-zA-Z0-9_.-]")),winname)!=1)
      {
        PrintErrorMessage('E',"openwindow","specify window name with n option");
        return (PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("openwindow",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (strlen(winname)==0)
    sprintf(winname,"window%d",(int) wincounter++);

  if (theOutDev==NULL) {
    PrintErrorMessage('E',"openwindow","no output device");
    return (PARAMERRORCODE);
  }

  if ((theWin=CreateUgWindow(theOutDev,winname,x,y,w,h))==NULL)
  {
    PrintErrorMessage('E',"openwindow","failed to open a window");
    return (CMDERRORCODE);
  }

  SetCurrentUgWindow(theWin);

  return (OKCODE);
}

/****************************************************************************/
/*D
   closewindow - close the current ug-window

   DESCRIPTION:
   This command closes one (or all) ug-window(s)
   (including the pictures residing there, of course),
   calling the functions 'DisposeMultiGrid' and 'DisposePicture'.

   'closewindow [$n <window name> | $a]'

   .  $n~<window~name>       - close the window with the specified name
   .n                        (default: the current window)

   .  $a                     - close all open windows
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   CloseWindowCommand - Close the current ug-window

   SYNOPSIS:
   static INT CloseWindowCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function closes the current ug-window.
   It closes one (or all) ug-window(s) (including the pictures residing there, of course).

   closewindow [$n <window name> | $a]

   .  $n <window name>       - close the window with the specified name
   .n                        (default: the current window)

   .  $a                     - close all open windows

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT CloseWindowCommand (INT argc, char **argv)
{
  UGWINDOW *theWin, *currWin;
  PICTURE *thePic,*currPic;
  char winname[NAMESIZE];
  INT i,aopt;

  currWin = GetCurrentUgWindow();
  theWin = currWin;

  /* check options */
  aopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'n' :
      if (sscanf(argv[i],expandfmt(CONCAT3("n %",NAMELENSTR,"[a-zA-Z0-9_.]")),winname)!=1)
      {
        PrintErrorMessage('E',"closewindow","specify a window name with n option");
        return (PARAMERRORCODE);
      }
      if ((theWin=GetUgWindow(winname))==NULL)
      {
        sprintf(buffer,"there is no window named '%s'",winname);
        PrintErrorMessage('W',"closewindow",buffer);
        return (OKCODE);
      }
      break;

    case 'a' :
      aopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("closewindow",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (theWin==NULL)
  {
    PrintErrorMessage('W',"closewindow","there's no current window");
    return (OKCODE);
  }
  if (theWin == currWin) currWin = NULL;
  currPic = GetCurrentPicture();

  if (aopt)
  {
    while((theWin=GetFirstUgWindow())!=NULL)
    {
      /* first dispose all pictures of that window */
      while((thePic=GetFirstPicture(theWin))!=NULL)
      {
        if (thePic==currPic)
          SetCurrentPicture(NULL);
        if (DisposePicture(thePic))
        {
          PrintErrorMessage('E',"closewindow","could not close a picture of that window");
          return (CMDERRORCODE);
        }
      }

      /* now dispose the window itself */
      if (DisposeUgWindow(theWin))
      {
        PrintErrorMessage('E',"closewindow","could not close the window");
        return (CMDERRORCODE);
      }
    }

    currWin = GetFirstUgWindow();

    SetCurrentUgWindow(currWin);
    if (currWin!=NULL)
      SetCurrentPicture(GetFirstPicture(currWin));
    else
      SetCurrentPicture(NULL);

    return (OKCODE);
  }

  /* first dispose all pictures of that window */
  while ((thePic=GetFirstPicture(theWin))!=NULL)
  {
    if (thePic==currPic)
      SetCurrentPicture(NULL);
    if (DisposePicture(thePic))
    {
      PrintErrorMessage('E',"closewindow","could not close a picture of that window");
      return (CMDERRORCODE);
    }
  }

  /* now dispose the window itself */
  if (DisposeUgWindow(theWin)!=0)
  {
    PrintErrorMessage('E',"closewindow","could not close the window");
    return (CMDERRORCODE);
  }

  currWin = GetFirstUgWindow();
  SetCurrentUgWindow(currWin);
  if (currWin!=NULL)
    SetCurrentPicture(GetFirstPicture(currWin));
  else
    SetCurrentPicture(NULL);

  return (OKCODE);
}

/****************************************************************************/
/*D
   setcurrwindow - set the current window

   DESCRIPTION:
   This command makes a window the current window.
   It calls the function 'GetUGWindow'.

   'setcurrwindow <window name>'

   .  <window~name> - name of a window

   SEE ALSO:
   'openwindow'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SetCurrentWindowCommand - Make a window the current window

   SYNOPSIS:
   static INT SetCurrentWindowCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function makes a window the current window.
   It makes a window the current window.

   setcurrwindow <window name>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SetCurrentWindowCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  char winname[NAMESIZE];

  NO_OPTION_CHECK(argc,argv);

  /* get window name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" setcurrwindow %",NAMELENSTR,"[a-zA-Z0-9_]")),winname)!=1)
  {
    PrintHelp("setcurrwindow",HELPITEM," (specify a window name)");
    return(PARAMERRORCODE);
  }

  theWin = GetUgWindow(winname);

  if (theWin==NULL)
  {
    PrintErrorMessage('E',"setcurrwindow","no window with this name open");
    return (CMDERRORCODE);
  }

  SetCurrentUgWindow(theWin);

  return (OKCODE);
}

/****************************************************************************/
/*D
   drawtext - draw text in a ug window

   DESCRIPTION:
   This command draws text into a ug window.

   'drawtext <xpos> <ypos> <text> [$w <window name>] [$c] [$i]'

   .  <xpos>                 - x-coordinate in pixels (origin is the lower left corner)
   .  <ypos>                 - y-coordinate in pixels (origin is the lower left corner)
   .  <text>                 - text to draw
   .  $w~<window~name>       - draw text into this window (default: current window)
   .  $i                     - draw text inverse
   .  $c                     - center text at <xpos> <ypos>
   .  $s~<size>              - text size
   D*/
/****************************************************************************/

static INT DrawTextCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  char winname[NAMESIZE],text[NAMESIZE];
  COORD_POINT pos;
  INT i,invopt,centeropt,size;

  theWin = GetCurrentUgWindow();
  if (theWin==NULL)
  {
    PrintErrorMessage('E',"drawtext","there's no window draw text");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],expandfmt(CONCAT3("drawtext %f %f %",NAMELENSTR,"[ -~]")),&pos.x,&pos.y,text)!=3)
  {
    PrintErrorMessage('E',"drawtext","specify position with two integers and then the text");
    return (CMDERRORCODE);
  }

  /* check options */
  invopt = centeropt = FALSE;
  size   = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'w' :
      if (sscanf(argv[i],expandfmt(CONCAT3("w %",NAMELENSTR,"[a-zA-Z0-9_]")),winname)!=1)
      {
        PrintErrorMessage('E',"drawtext","specify a window name with w option");
        return (PARAMERRORCODE);
      }
      if ((theWin=GetUgWindow(winname))==NULL)
      {
        sprintf(buffer,"there is no window named '%s'",winname);
        PrintErrorMessage('E',"drawtext",buffer);
        return (PARAMERRORCODE);
      }
      break;

    case 's' :
      if (sscanf(argv[i],"s %d",&size)!=1)
      {
        PrintErrorMessage('E',"drawtext","specify a size with s option");
        return (PARAMERRORCODE);
      }
      break;

    case 'i' :
      invopt = TRUE;
      break;

    case 'c' :
      centeropt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("drawtext",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  DrawWindowText(theWin,pos,text,size,centeropt,invopt);

  return (OKCODE);
}

/****************************************************************************/
/*D
   openpicture - open a new picture

   DESCRIPTION:
   This command opens a picture on a window
   (these will be the current window and picture then).
   It calls the function 'CreatePicture'.

   'openpicture [$w <window name>] [$s <h> <v> <dh> <dv>] [$n <picture name>]'

   .  $w~<window~name>       - open a picture on this window (default: current window)
   .  $s~<h>~<v>~<dh>~<dv>   - specify the location and size in the --> standardRefSys with the origin located in the lower left corner of the parent window
   .n                        (default: picture size = parent window size)

   .  $n~<picture~name>      - optionally you can specify the picture name
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   OpenPictureCommand - open a new picture

   SYNOPSIS:
   static INT OpenPictureCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function opens a new window.
   It opens a picture on a window (these will be the current window and picture then).

   openpicture [$w <window name>] [$s <h> <v> <dh> <dv>] [$n <picture name>]

   .  $w~<window~name>       - open a picture on this window (default: current window)
   .  $s~<h>~<v>~<dh>~<dv>   - specify the location and size in the --> standardRefSys with the origin located in the lower left corner of the parent window
   .n                        (default: picture size = parent window size)

   .  $n~<picture name>      - optionally you can specify the picture name

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT OpenPictureCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  PICTURE *thePic;
  char picname[NAMESIZE],winname[NAMESIZE];
  INT local_LL[2],local_UR[2];
  INT i,sopt;

  /* following variables: keep type for sscanf */
  int h,v,dh,dv;


  theWin = GetCurrentUgWindow();
  if (theWin==NULL)
  {
    PrintErrorMessage('E',"openpicture","there's no window to open a picture on");
    return (CMDERRORCODE);
  }

  /* check options */
  picname[0] = '\0';
  sopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'n' :
      if (sscanf(argv[i],expandfmt(CONCAT3("n %",NAMELENSTR,"[a-zA-Z0-9_]")),picname)!=1)
      {
        PrintErrorMessage('E',"openpicture","specify a picture name with n option");
        return (PARAMERRORCODE);
      }
      break;

    case 'w' :
      if (sscanf(argv[i],expandfmt(CONCAT3("w %",NAMELENSTR,"[a-zA-Z0-9_]")),winname)!=1)
      {
        PrintErrorMessage('E',"openpicture","specify a window name with w option");
        return (PARAMERRORCODE);
      }
      if ((theWin=GetUgWindow(winname))==NULL)
      {
        sprintf(buffer,"there is no window named '%s'",winname);
        PrintErrorMessage('E',"openpicture",buffer);
        return (PARAMERRORCODE);
      }
      break;

    case 's' :
      sopt = TRUE;
      if (sscanf(argv[i],"s %d %d %d %d",&h,&v,&dh,&dv)!=4)
      {
        PrintErrorMessage('E',"openpicture","specify h, v, dh, dv with s option");
        return (PARAMERRORCODE);
      }
      local_LL[0] = h;
      local_LL[1] = v;
      local_UR[0] = h+dh;
      local_UR[1] = v+dv;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("openpicture",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (strlen(picname)==0)
    sprintf(picname,"picture%d",(int) piccounter++);

  if (!sopt)
  {
    /* default size is the window's size */
    local_LL[0] = 0;
    local_LL[1] = 0;
    local_UR[0] = ABS(UGW_LUR(theWin)[0]-UGW_LLL(theWin)[0]);
    local_UR[1] = ABS(UGW_LUR(theWin)[1]-UGW_LLL(theWin)[1]);
  }

  if ((thePic=CreatePicture(picname,theWin,local_LL,local_UR))==NULL)
  {
    PrintErrorMessage('E',"openpicture","failed to open a picture");
    return (CMDERRORCODE);
  }

  SetCurrentPicture(thePic);

  return (OKCODE);
}

/****************************************************************************/
/*D
   closepicture - close a picture

   DESCRIPTION:
   This command closes one (or all) picture(s) on a window.
   It calls the function 'DeletePicture'.

   'closepicture [$a | {$w <window name> {<picture name> | $a}}]'

   .  $w~<window~name>        - close a picture of this window
   (default: the current picture)
   .  $a                      - close all pictures of the current window
   .  {<picture~name>~|~$a}   - close the picture with the specified name or all pictures of that window
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ClosePictureCommand - Close the current view

   SYNOPSIS:
   static INT ClosePictureCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function closes the current view.
   It closes one (or all) picture(s) on a window.

   closepicture [$a | {$w <window name> {<picture name> | $a}}]

   .  default                 - the current picture
   .  $a                      - close all pictures of the current window
   .  $w <window name>        - close a picture of this window
   .  {<picture name> | $a}   - close the picture with the specified name or all pictures of that window

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ClosePictureCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  PICTURE *thePic,*NextPic;
  char picname[NAMESIZE],winname[NAMESIZE];
  INT i,aopt,wopt;

  theWin = GetCurrentUgWindow();
  if (theWin==NULL)
  {
    PrintErrorMessage('W',"closepicture","there's no open window");
    return (OKCODE);
  }

  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('W',"closepicture","there's no picture to dispose");
    return (OKCODE);
  }

  /* check options */
  aopt = wopt = FALSE;
  picname[0] = '\0';
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'w' :
      wopt = TRUE;
      if (sscanf(argv[i],expandfmt(CONCAT5("w %",NAMELENSTR,"[a-zA-Z0-9_] %",NAMELENSTR,"[a-zA-Z0-9_]")),winname,picname)<1)
      {
        PrintErrorMessage('E',"closepicture","specify a window name with w option");
        return (PARAMERRORCODE);
      }
      if ((theWin=GetUgWindow(winname))==NULL)
      {
        sprintf(buffer,"there is no window named '%s'",winname);
        PrintErrorMessage('E',"closepicture",buffer);
        return (PARAMERRORCODE);
      }
      break;

    case 'a' :
      aopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("closepicture",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (wopt)
  {
    if (!aopt && (strlen(picname)==0))
    {
      PrintErrorMessage('E',"closepicture","specify picture name or $a with window name");
      return (PARAMERRORCODE);
    }
    else if (strlen(picname)==0)
      if ((thePic=GetUgPicture(theWin,picname))==NULL)
      {
        sprintf(buffer,"there is no picture named '%s'",picname);
        PrintErrorMessage('E',"closepicture",buffer);
        return (PARAMERRORCODE);
      }
  }

  if (aopt)
  {
    for (thePic=GetFirstPicture(theWin); thePic!=NULL; thePic=NextPic)
    {
      NextPic = GetNextPicture(thePic);
      if (DisposePicture(thePic)!=0)
      {
        PrintErrorMessage('E',"closepicture","could not close the picture");
        return (CMDERRORCODE);
      }
    }
    SetCurrentPicture(NULL);

    return (OKCODE);
  }

  SetCurrentPicture(NULL);
  if (ErasePicture(thePic)) return (CMDERRORCODE);

  if (DisposePicture(thePic)!=0)
  {
    PrintErrorMessage('E',"closepicture","could not close the picture");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   setcurrpicture - make a picture the current picture

   DESCRIPTION:
   This command makes a picture the current picture.

   'setcurrpicture <picture name> [$w <window name>]'

   .  <picture~name>       - name of the picture
   .  $w~<window~name>       - picture resides in this window (default: current window)
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SetCurrentPictureCommand - make a picture the current picture

   SYNOPSIS:
   static INT SetCurrentPictureCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function makes a picture the current picture.
   It makes a picture the current picture.

   'setcurrpicture <picture name> [$w <window name>]'

   .  $w<window~name>       - picture resides in this window (default: current window)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SetCurrentPictureCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  PICTURE *thePic;
  char picname[NAMESIZE], winname[NAMESIZE];
  INT i;

  theWin = GetCurrentUgWindow();
  if (theWin==NULL)
  {
    PrintErrorMessage('E',"setcurrpicture","there's no open window (and therefore no picture)");
    return (CMDERRORCODE);
  }

  /* get picture name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" setcurrpicture %",NAMELENSTR,"[a-zA-Z0-9_]")),picname)!=1)
  {
    PrintHelp("setcurrpicture",HELPITEM," (specify a picture name)");
    return(PARAMERRORCODE);
  }

  /* check options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'w' :
      if (sscanf(argv[i],expandfmt(CONCAT3("w %",NAMELENSTR,"[a-zA-Z0-9_]")),winname)!=1)
      {
        PrintErrorMessage('E',"setcurrpicture","specify a window name with w option");
        return (PARAMERRORCODE);
      }
      if ((theWin=GetUgWindow(winname))==NULL)
      {
        sprintf(buffer,"there is no window named '%s'",winname);
        PrintErrorMessage('E',"setcurrpicture",buffer);
        return (PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("setcurrpicture",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  thePic = GetUgPicture(theWin,picname);

  if (thePic==NULL)
  {
    PrintErrorMessage('E',"setcurrpicture","no picture with this name open");
    return (CMDERRORCODE);
  }

  SetCurrentPicture(thePic);

  return (OKCODE);
}

/****************************************************************************/
/*D
   clearpicture - clear current picture

   DESCRIPTION:
   This command clears current picture.
   It calls the function 'ErasePicture'.

   'clearpicture'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ClearPictureCommand - clear current picture

   SYNOPSIS:
   static INT ClearPictureCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function clears current picture.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ClearPictureCommand (INT argc, char **argv)
{
  PICTURE *thePicture;

  NO_OPTION_CHECK(argc,argv);

  thePicture = GetCurrentPicture();
  if (thePicture==NULL)
  {
    UserWrite("WARNING: there is no current picture\n");
    return (OKCODE);
  }
  ErasePicture(thePicture);
  DrawPictureFrame(thePicture,WOP_ACTIVE);

  /* picture has changed */
  if (InvalidatePicture(thePicture))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   picframe - toggle framing of pictures

   DESCRIPTION:
   This command toggles the framing of pictures.

   SYMTAX:
   'picframe 0|1'
   D*/
/****************************************************************************/

static INT PicFrameCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  if (strchr(argv[0],'0')!=NULL)
    SetDoFramePicture (NO);
  else
    SetDoFramePicture (YES);

  return (OKCODE);
}
/****************************************************************************/
/*D
   setview - specifies the view on the object

   DESCRIPTION:
   This command specifies or changes the observer view of the object
   of the current picture.
   It calls the function 'SetView'.

    in 2D: 'setview [$i] [$t <x> <y>] [$x  <x> <y>]'

    in 3D: 'setview [$i] [$o <x> <y> <z> $t <x> <y> <z>] [$x <x> <y> [<z>]] [$p < | =]'

   .n                         all coordinates have to be given in physical coordinates

   .   $i                     - return to default settings first
   .   $o~<x>~<y>~<z>         - 3D objects ONLY: specify the observer stand
   .   $p~<~|~=               - 3D objects ONLY: choose central (<) or parallel (=) perspective
   .   $t~<x>~<y>~[<z>]       - specify the target point in the viewplane
   .n                         (NB: the viewplane is then defined to be normal to the line observer-target)

   .    $x~<x>~<y>~[<z>]      - define an x-axis in the viewplane (which will be to the right in the picture)
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SetViewCommand - Specify or change the observer's view

   SYNOPSIS:
   static INT SetViewCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function specifies or changes the observer's view of
   the object of the current picture.
   It specifies or change the observer's view of the object of the current picture.

    in 2D: 'setview [$i] [$t <x> <y>] [$x  <x> <y>]'

    in 3D: 'setview [$i] [$o <x> <y> <z> $t <x> <y> <z>] [$x <x> <y> [<z>]] [$p < | =]'

   .n                         all coordinates have to be given in physical coordinates

   .   $i                     - return to default settings first
   .   $o <x> <y> <z>         - 3D objects ONLY: specify the observer stand
   .   $p < | =               - 3D objects ONLY: choose central (<) or parallel (=) perspective
   .   $t <x> <y> [<z>]       - specify the target point in the viewplane
   .n                         (NB: the viewplane is then defined to be normal to the line observer-target)

   .    $x <x> <y> [<z>]      - define an x-axis in the viewplane (which will be to the right in the picture)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SetViewCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  VIEWEDOBJ *theViewedObj;
  COORD *viewPoint,*targetPoint,*xAxis;
  COORD vP[3],tP[3],xA[3];
  INT *perspective;
  INT per;
  INT i,veclen;

  /* following variables: keep type for sscanf */
  float x[3];

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"setview","there's no current picture");
    return (CMDERRORCODE);
  }

  theViewedObj = PIC_VO(thePic);
  veclen = (VO_DIM(theViewedObj)==TYPE_2D) ? 2 : 3;

  /* check options */
  viewPoint = targetPoint = xAxis = NULL;
  perspective = NULL;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'o' :
      if (VO_DIM(theViewedObj)!=TYPE_3D)
      {
        PrintErrorMessage('E',"setview","the o option applies ONLY with 3D objects");
        return (PARAMERRORCODE);
      }
      if (sscanf(argv[i],"o %f %f %f",x,x+1,x+2)!=veclen)
      {
        sprintf(buffer,"o option: %d coordinates required for a %dD object",(int)veclen,(int)veclen);
        PrintErrorMessage('E',"setview",buffer);
        return (PARAMERRORCODE);
      }
      for (i=0; i<veclen; i++)
        vP[i] = x[i];
      viewPoint = vP;
      break;

    case 't' :
      if (sscanf(argv[i],"t %f %f %f",x,x+1,x+2)!=veclen)
      {
        sprintf(buffer,"t option: %d coordinates required for a %dD object",(int)veclen,(int)veclen);
        PrintErrorMessage('E',"setview",buffer);
        return (PARAMERRORCODE);
      }
      for (i=0; i<veclen; i++)
        tP[i] = x[i];
      targetPoint = tP;
      break;

    case 'x' :
      if (sscanf(argv[i],"x %f %f %f",x,x+1,x+2)!=veclen)
      {
        sprintf(buffer,"x option: %d coordinates required for a %dD object",(int)veclen,(int)veclen);
        PrintErrorMessage('E',"setview",buffer);
        return (PARAMERRORCODE);
      }
      for (i=0; i<veclen; i++)
        xA[i] = x[i];
      xAxis = xA;
      break;

    case 'p' :
      if (VO_DIM(theViewedObj)!=TYPE_3D)
      {
        PrintErrorMessage('E',"setview","the p option applies ONLY with 3D objects");
        return (PARAMERRORCODE);
      }
      if ((strchr(argv[i],'<')!=NULL) && (strchr(argv[i],'=')!=NULL))
      {
        PrintErrorMessage('E',"setview","specify EITHER < OR = for the perspective");
        return (PARAMERRORCODE);
      }
      else if (strchr(argv[i],'<')!=NULL)
        per = YES;
      else if (strchr(argv[i],'=')!=NULL)
        per = NO;
      else
      {
        PrintErrorMessage('E',"setview","specify AT LEAST < OR = for the perspective");
        return (PARAMERRORCODE);
      }
      perspective = &per;
      break;

    case 'i' :
      VO_STATUS(theViewedObj) = NOT_INIT;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("setview",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (SetView(thePic,viewPoint,targetPoint,xAxis,perspective)!=0)
  {
    PrintErrorMessage('E',"setview","error during SetView");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   vdisplay - display view of current picture

   DESCRIPTION:
   This command displays view of current picture.
   It calls the function 'DisplayViewOfDisplayedObject'.

   'vdisplay'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   DisplayViewCommand - Display view of current picture

   SYNOPSIS:
   static INT DisplayViewCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function displays view of current picture.
   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT DisplayViewCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"vdisplay","there's no current picture");
    return (CMDERRORCODE);
  }

  if (DisplayViewOfViewedObject(thePic))
  {
    PrintErrorMessage('E',"vdisplay","error during DisplayView");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   walk - let the observer walk relative to the viewRefSys

   DESCRIPTION:
   This command lets the observer walk relative to the viewRefSys
   in the current picture.
   It calls the function 'walk'.

   'walk <x> <y> [<z>]'

   .   <x>~<y>~[<z>] - coordinates
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   WalkCommand - let the observer walk relative to the viewRefSys

   SYNOPSIS:
   static INT WalkCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lets the observer walk relative to the viewRefSys
   in the current picture.

   walk <x> <y> [<z>]

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT WalkCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  VIEWEDOBJ *theViewedObj;
  COORD dx[3];
  INT i,veclen;

  /* following variables: keep type for sscanf */
  float x[3];


  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"walk","there's no current picture");
    return (CMDERRORCODE);
  }

  theViewedObj = PIC_VO(thePic);
  veclen = (VO_DIM(theViewedObj)==TYPE_2D) ? 2 : 3;

  if (sscanf(argv[0],"walk %f %f %f",x,x+1,x+2)!=veclen)
  {
    sprintf(buffer,"%d coordinates required for a %dD object",(int)veclen,(int)veclen);
    PrintErrorMessage('E',"walk",buffer);
    return (PARAMERRORCODE);
  }
  for (i=0; i<veclen; i++)
    dx[i] = x[i];

  if (Walk(thePic,dx)!=0)
  {
    PrintErrorMessage('E',"walk","error during Walk");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   walkaround - let the observer walk on a sphere around the target point

   DESCRIPTION:
   This command lets the observer walk on a sphere around the target point.
   It calls the function 'RunAroundTargetPoint'.

   'walkaround <viewplane angle> <rotation angle>'

   .  <viewplane~angle>      - this angle runs in the view plane math pos from the x-axis and defines
   .n                        - together with the target-observer direction a plane

   .  <rotation~angle>       - the observer will be rotated around the target point in the above plane

   (angles in degree 0 - 360)
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   WalkAroundCommand - let the observer walk on a sphere around the target point

   SYNOPSIS:
   static INT WalkAroundCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lets the observer walk on a sphere around
   the target point in the current picture. (3D ONLY)
   It lets the observer walk on a sphere around the target point in the current picture.

   walkaround <viewplane angle> <rotation angle>

   .  <viewplane~angle>      - this angle runs in the view plane math pos from the x-axis and defines
   .n                        - together with the target-observer direction a plane

   .  <rotation~angle>       - the observer will be rotated around the target point in the above plane

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT WalkAroundCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  VIEWEDOBJ *theViewedObj;

  /* following variables: keep type for sscanf */
  float dirAngle,Angle;


  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"walkaround","there's no current picture");
    return (CMDERRORCODE);
  }

  theViewedObj = PIC_VO(thePic);

  if (VO_DIM(theViewedObj)!=TYPE_3D)
  {
    PrintErrorMessage('E',"walkaround","walkaround only possible for 3D objects");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"walkaround %f %f",&dirAngle,&Angle)!=2)
  {
    PrintErrorMessage('E',"walkaround","2 angles required");
    return (PARAMERRORCODE);
  }

  /* transform degree into arclength */
  dirAngle *= PI/180.0;
  Angle    *= PI/180.0;

  if (RunAroundTargetPoint(thePic,dirAngle,Angle)!=0)
  {
    PrintErrorMessage('E',"walkaround","error during WalkAroundTargetPoint");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   zoom - zoom the current picture

   DESCRIPTION:
   This command zooms the current picture.
   the zoom factor is always relative to the current setting.
   It calls the function 'Zoom'.

   'zoom <factor>'

   .       <factor>              - values < 1 magnify picture

   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ZoomCommand - zoom the current picture

   SYNOPSIS:
   static INT ZoomCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function zooms the current picture.

   zoom <factor>

   .  zoom <factor>              - values < 1 magnify picture

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ZoomCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  /* following variables: keep type for sscanf */
  float factor;


  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"zoom","there's no current picture");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"zoom %f",&factor)!=1)
  {
    PrintErrorMessage('E',"zoom","zoom factor required");
    return (PARAMERRORCODE);
  }

  if (Zoom(thePic,factor)!=0)
  {
    PrintErrorMessage('E',"zoom","error during Zoom");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   drag - drag the projection plane of the current picture

   DESCRIPTION:
   This command drags the projection plane of the current picture relative
   to its x-axis.
   It calls the function 'DragProjectionPlane'.

   drag <dx> <dy>

   .  <dx>~<dy>  - displacement vector

   D*/
/****************************************************************************/

/****************************************************************************/
/*
   DragCommand - drag the projection plane of the current picture

   SYNOPSIS:
   static INT DragCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function drags the projection plane of the current picture
   relative to its x-axis.

   drag <dx> <dy>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT DragCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  /* following variables: keep type for sscanf */
  float dx,dy;


  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"drag","there's no current picture");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"drag %f %f",&dx,&dy)!=2)
  {
    PrintErrorMessage('E',"drag","dx, dy required");
    return (PARAMERRORCODE);
  }

  if (DragProjectionPlane(thePic,dx,dy)!=0)
  {
    PrintErrorMessage('E',"drag","error during DragProjectionPlane");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   rotate - rotate the projection plane of the current picture

   DESCRIPTION:
   This command rotates the projection plane of the current picture around
   the target point.
   It calls the function 'RotateProjectionPlane'.

   'rotate <angle>'

   .  <angle>                - <angle> runs in the view plane math pos from the x-axis
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   RotateCommand - Rotate the projection plane of the current picture

   SYNOPSIS:
   static INT RotateCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function rotates the projection plane of the current picture
   around the target point.

   rotate <angle>

   .  <angle>                - <angle> runs in the view plane math pos from the x-axis

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT RotateCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  /* following variables: keep type for sscanf */
  float angle;


  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"rotate","there's no current picture");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"rotate %f",&angle)!=1)
  {
    PrintErrorMessage('E',"rotate","angle required");
    return (PARAMERRORCODE);
  }

  /* transform degree into arclength */
  angle *= PI/180.0;

  if (RotateProjectionPlane(thePic,angle)!=0)
  {
    PrintErrorMessage('E',"rotate","error during RotateProjectionPlane");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   setplotobject - set plotting specification

   DESCRIPTION:
   This command specifies the object which will be plotted in the current
   picture and associates the current multigrid with it (if not done yet) .
   It calls the function 'SpecifyPlotObjOfViewedObject'.

   'setplotobject [<object type name>] [$a] ... '

   .    <object~type~name>   - possible are: 'Grid', 'EScalar', 'EVector'
   .    $a - (for 3d plot objects) moves the observer to a place outside
                   the bounding sphere of the plot object


   The remaining options depend on which object you specified.

   2D:

   'Grid [$c {0|1}] [$e {0|1}] [$n {0|1}] [$w {0|1}]'

   .    $c~{0|1}      - plot colored: no/yes
   .    $e~{0|1}      - plot ElemIDs: no/yes
   .    $n~{0|1}      - plot NodeIDs: no/yes
   .    $w~{c/i/r/a}  - which elements to plot:
   .n                                        c  copies and up
   .n                                        i  irregular and up
   .n                                        r  regular and up
   .n                                        a  all

    'EScalar $e <ElemEvalProc> [$m {COLOR | CONTOURS_EQ}] [$f <fromValue>] [$t <toValue>] [$n <nContours>]'

   .    $e~<ElemEvalProc>               - name of element scalar evaluation procedure
   .    $d~<depth>                      - depth of plot
   .    $m~{COLOR|CONTOURS_EQ}        - mode: COLOR-plot or CONTOUR-plot
   .    $f~<fromValue>~$t~<toValue>     - range [fromValue,toValue]
   .    $n~<nContours>                  - number of contours

   'EVector $e <ElemEvalProc> [$c {0|1}] [$t <toValue>] [$r <rastersize>] [$l <cutlength>]'

   .    $e~<ElemEvalProc>               - name of element vector evaluation procedure
   .    $c~{0|1}                        - cut vectors if to long
   .    $t~<toValue>                    - range: [0,toValue]
   .    $r~<rastersize>                 - physical meshsize of rasterpoints where
                                     vectors are plotted
   .    $l~<cutlength>                  - cutlength in units of rastersize

   3D:

   'Grid [$c {0|1}] [$w {0|1}] [$P <x> <y> <z>] [$N <x> <y> <z>]'

   .    $c~{0|1}                        - plot colored: no/yes
   .    $w~{c/i/r/a}                    - which elements to plot:
   .n                                        c  copies and up
   .n                                        i  irregular and up
   .n                                        r  regular and up
   .n                                        a  all

   .    $P~<x>~<y>~<z>                  - a point on the cut plane
   .    $N~<x>~<y>~<z>                  - the normal of the cut plane


   'EScalar $e <ElemEvalProc> [$m {COLOR | CONTOURS_EQ}] [$f <fromValue>] [$t <toValue>] [$n <nContours>] [$P <x> <y> <z>] [$N <x> <y> <z>]'

   .    $e~<ElemEvalProc>               - name of element scalar evaluation procedure
   .    $d~<depth>                      - depth of plot
   .    $m~{COLOR|CONTOURS_EQ}          - mode: COLOR-plot or CONTOUR-plot
   .    $f~<fromValue>~$t~<toValue>     - range [fromValue,toValue]
   .    $n~<nContours>                  - number of contours
   .    $P~<x>~<y>~<z>                  - a point on the cut plane
   .    $N~<x>~<y>~<z>                  - the normal of the cut plane


    'EVector $e <ElemEvalProc> [$c {0|1}] [$t <toValue>] [$r <rastersize>] [$P <x> <y> <z>] [$N <x> <y> <z>]'

   .    $e~<ElemEvalProc>               - name of element vector evaluation procedure
   .    $c~{0|1}                        - cut vectors if to long
   .    $t~<toValue>                    - range: [0,toValue]
   .    $r~<rastersize>                 - physical meshsize of rasterpoints where
   .    $P~<x>~<y>~<z>                  - a point on the cut plane
   .    $N~<x>~<y>~<z>                  - the normal of the cut plane
                                    - vectors are plotted

   EXAMPLE:
   .vb
   openpicture $n p0 $s 5 5  300 300;
   setplotobject Grid $b 1 $n 0 $e 0 $c 1 $w a;
   setview;
   zoom 0.75;

   openpicture $n p1 $s 5 310 300 300;
   setplotobject EScalar $e u_sol $d 1;
   setview;
   zoom 0.75;

   openpicture $n p2 $s 5 615 300 300;
   setplotobject EScalar $e v_sol $d 1;
   .ve
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SetPlotObjectCommand - Specify the object

   SYNOPSIS:
   static INT SetPlotObjectCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function specifies the object which will be plotted in the current
   picture and associates the current multigrid with it (if not done yet) .
   It specifies the object which will be plotted in the current picture and associates
   the current multigrid with it (if not done yet).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SetPlotObjectCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  MULTIGRID *theMG;
  char potname[NAMESIZE],*thePlotObjTypeName;

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"setplotobject","there's no current picture");
    return (CMDERRORCODE);
  }

  /* get plot object name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" setplotobject %",NAMELENSTR,"[a-zA-Z0-9_]")),potname)==1)
  {
    thePlotObjTypeName = potname;

    /* current multigrid */
    theMG = currMG;
    if (theMG==NULL)
    {
      PrintErrorMessage('E',"setplotobject","no current multigrid\n");
      return (CMDERRORCODE);
    }
  }
  else
  {
    thePlotObjTypeName = NULL;
    theMG = NULL;
  }

  if (theMG!=NULL)
  {
    sprintf(buffer," picture '%s' and multigrid '%s' coupled\n",ENVITEM_NAME(thePic),ENVITEM_NAME(theMG));
    UserWrite(buffer);
  }

  if (SpecifyPlotObjOfViewedObject(thePic,theMG,thePlotObjTypeName,argc,argv)!=0)
  {
    PrintErrorMessage('E',"setplotobject","error during SpecifyPlotObjOfViewedObject");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);


  return (OKCODE);
}

/****************************************************************************/
/*D
   polist - print the specifications of the object

   DESCRIPTION:
   This command prints the specifications of the object defined in the
   current picture.
   It calls the function 'DisplayObject'.

   'polist'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   PlotObjectListCommand -  print the specifications of the object

   SYNOPSIS:
   static INT PlotObjectListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function prints the specifications of the object defined in the current picture.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PlotObjectListCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('W',"listplotobject","there's no current picture");
    return (OKCODE);
  }

  if (DisplayObject(PIC_PO(thePic))!=0)
  {
    PrintErrorMessage('E',"listplotobject","error during DisplayPlotObjOfViewedObject");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   plot - plot an object

   DESCRIPTION:
   This command plots the object of the current picture according to
   its specifications.
   It calls the function 'WorkOnPicture'.

   'plot'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   PlotCommand - Plot the object of the current picture

   SYNOPSIS:
   static INT PlotCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function plots the object of the current picture according to its specifications.
   It prints the specifications of the object defined in the current picture.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PlotCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  WORK myWork,*theWork;

  NO_OPTION_CHECK(argc,argv);

  theWork = &myWork;

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"plot","there's no current picture");
    return (CMDERRORCODE);
  }

  /* fill work struct */
  W_ID(theWork) = DRAW_WORK;

  if (WorkOnPicture(thePic,theWork)!=0)
  {
    PrintErrorMessage('E',"plot","error during WorkOnPicture");
    return (CMDERRORCODE);
  }

  /* picture is current */
  DrawPictureFrame(thePic,WOP_ACTIVE);

  return (OKCODE);
}

/****************************************************************************/
/*D
   findrange - find the range

   DESCRIPTION:
   This command computes the range of  values to be plotted for the object
   in the current picture. The result will be printed and
   stored into :findrange:min and :findrange:max.
   It calls the function 'WorkOnPicture'.

   'findrange [$z <zoom factor>] [$s] [$p]'

   .  $z~<zoom factor>       - zoom the range by this factor
   .  $s                     - symmetrize the range (min=max)
   .  $p                     - store range in plot object

   SEE ALSO:
   'plot'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   FindRangeCommand -  Plot the object of the current picture

   SYNOPSIS:
   static INT FindRangeCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function plots the object of the current picture according to its specifications.
   It finds the range of values to be plotted for the object in the current picture.
   The result will be printed and stored into :findrange:min and :findrange:max.

   findrange [$z <zoom factor>] [$s] [$p]

   .  $z <zoom factor>       - zoom the range by this factor
   .  $s                     - symmetrize the range (min=max)
   .  $p                     - store range in plot object

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT FindRangeCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  WORK myWork,*theWork;
  INT i,sym,put;
  DOUBLE min,max;

  /* following variables: keep type for sscanf */
  float zoom;

  theWork = &myWork;

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"findrange","there's no current picture");
    return (CMDERRORCODE);
  }

  /* check options */
  put = sym = NO;
  zoom = 1.0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      sym = YES;
      break;

    case 'p' :
      put = YES;
      break;

    case 'z' :
      if (sscanf(argv[i],"z %f",&zoom)!=1)
      {
        PrintErrorMessage('E',"findrange","specify a zoom factor with z option");
        return (PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("findrange",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* fill work struct */
  W_ID(theWork) = FINDRANGE_WORK;

  W_FINDRANGE_WORK(theWork)->symmetric = sym;
  W_FINDRANGE_WORK(theWork)->zoom          = zoom;
  W_FINDRANGE_WORK(theWork)->put           = put;

  if (WorkOnPicture(thePic,theWork)!=0)
  {
    PrintErrorMessage('E',"findrange","error during WorkOnPicture");
    return (CMDERRORCODE);
  }

  /* read out min and max */
  min = W_FINDRANGE_WORK(theWork)->min;
  max = W_FINDRANGE_WORK(theWork)->max;

  sprintf(buffer," FR_min = %10.3g\n FR_max = %10.3g\n",(float)min,(float)max);
  UserWrite(buffer);

  if (    (SetStringValue(":findrange:min",min)!=0)
          ||      (SetStringValue(":findrange:max",max)!=0))
  {
    PrintErrorMessage('E',"findrange","could not set :findrange:min or :findrange:max");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   setcurrmg - set a current multigrid

   DESCRIPTION:
   This command sets a current multigrid.

   'setcurrmg <mgname>'

   .    <mgname>  - name of the open multigrid which will be the current one
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SetCurrentMultigridCommand - Make a multigrid the current multigrid

   SYNOPSIS:
   static INT SetCurrentMultigridCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function makes a multigrid the current multigrid.
   It specifies the object which will be plotted in the current picture and associates
   the current multigrid with it (if not done yet).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SetCurrentMultigridCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char mgname[NAMESIZE];

  NO_OPTION_CHECK(argc,argv);

  /* get multigrid name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" setcurrmg %",NAMELENSTR,"[ -~]")),mgname)!=1)
  {
    PrintHelp("setcurrmg",HELPITEM," (specify current multigrid name)");
    return(PARAMERRORCODE);
  }

  theMG = GetMultigrid(mgname);

  if (theMG==NULL)
  {
    PrintErrorMessage('E',"setcurrmg","no multigrid with this name open");
    return (CMDERRORCODE);
  }

  if (SetCurrentMultigrid(theMG)!=0)
    return (CMDERRORCODE);

  return(OKCODE);
}

/****************************************************************************/
/*D
   updateDoc - updating the windows and pictures of the current multigrid

   DESCRIPTION:
   This command runs 'InvalidatePicturesOfMG' and
   'InvalidateUgWindowsOfMG'.
   If the refresh state is on, the pictures will be replotted.

   'updateDoc'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   UpdateDocumentCommand - Redraw pictures of current multigrid

   SYNOPSIS:
   static INT UpdateDocumentCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function redraws pictures of current multigrid.
   It redraws all pictures of the current multigrid.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT UpdateDocumentCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  InvalidatePicturesOfMG(currMG);
  InvalidateUgWindowsOfMG(currMG);

  return (OKCODE);
}

/****************************************************************************/
/*D
   clear - assign a value to a symbolic vector or matrix

   DESCRIPTION:
   This function sets the value of a symbol.
   It clears or assign a constant value to a user defined symbol.

   'clear <symbol name> [$a] [$u] [$v <value>]'

   .  $a                     - from level 0 through current level (default: current level only)
   .  $u                     - consider the skip fields
   .  $v~<value>             - assign this value (instead of 0.0)

   SEE ALSO:
   'cv', 'cm'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ClearCommand - Set the value of a symbol

   SYNOPSIS:
   static INT ClearCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function sets the value of a symbol.
   It clears or assign a constant value to a user defined symbol.

   clear <symbol name> [$a] [$s] [$v <value>]

   .  $a                     - from level 0 through current level (default: current level only)
   .  $s                     - consider the skip fields
   .  $v <value>             - assign this value (instead of 0.0)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ClearCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  SYMBOL *sym;
  INT i,fl,tl,skip;
  float value;

  theMG = GetCurrentMultigrid();
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"clear","no current multigrid");
    return(CMDERRORCODE);
  }

  /* check options */
  fl = tl = CURRENTLEVEL(theMG);
  skip = FALSE;
  value = 0.0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      fl = 0;
      break;

    case 's' :
      skip = TRUE;
      break;

    case 'v' :
      if (sscanf(argv[i],"v %f",&value)!=1)
      {
        PrintErrorMessage('E',"clear","could not read value");
        return(CMDERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("clear",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if ((sym = ReadVecSymbolOfFormat(ENVITEM_NAME(MGFORMAT(theMG)),
                                   "clear",argc,argv))==NULL)
  {
    PrintErrorMessage('E',"clear","could not read symbol");
    return (PARAMERRORCODE);
  }

  if (skip)
  {
    if (a_dsetnonskip(theMG,fl,tl,SYM_VEC_DESC(sym),EVERY_CLASS,value)
        !=NUM_OK)
      return (CMDERRORCODE);
  }
  else
  {
    if (a_dset(theMG,fl,tl,SYM_VEC_DESC(sym),EVERY_CLASS,value)!=NUM_OK)
      return (CMDERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   reinit - reinitialize a problem

   DESCRIPTION:
   This command reinitializes the problem with the user defined reinit of the
   problem

   'reinit [$d <domain> $p <problem>] [reinit options]'

   .  no~option			  - reinit the problem of the current multigrid
   .  $d~<domain>            - domain to initialize (default is the domain of the current mg)
   .  $p~<problem>           - problem to initialize (default is the problem of the current mg)
   .  reinit~options		  - these options will be passed to the problem defined reinit procedure
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   reinit -

   SYNOPSIS:
   static INT ReInitCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function reinitializes the problem.
   It reinitializes a problem.

   reinit [$d <domain> $p <problem>]

   .  $d <domain>            - domain to initialize (default is the domain of the current mg)
   .  $p <problem>           - problem to initialize (default is the problem of the current mg)

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ReInitCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  INT i,bopt;
  char BVPName[NAMESIZE];

  /* check options */
  bopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'b' :
      if (sscanf(argv[i],expandfmt(CONCAT3("b %",NAMELENSTR,"[0-9a-zA-Z/_ ]")),BVPName)!=1)
      {
        PrintErrorMessage('E',"reinit","could not read BndValProblem string");
        return(PARAMERRORCODE);
      }
      bopt = TRUE;
      break;

      /* no default because param list is passed to reinit function */
    }

  /* set the document */
  if (bopt)
  {
    theBVP = GetBVP(BVPName);
    if(theBVP==NULL)
    {
      sprintf(buffer,"could not interpret '%s' as a BVP name",BVPName);
      PrintErrorMessage('E',"reinit",buffer);
      return (CMDERRORCODE);
    }
  }
  else
  {
    theMG = currMG;
    if (theMG==NULL)
    {
      PrintErrorMessage('E',"reinit","no open multigrid (specify problem and domain instead)");
      return (CMDERRORCODE);
    }
    theBVP = MG_BVP(theMG);
  }
  if (BVP_GetBVPDesc(theBVP,&theBVPDesc)) return (CMDERRORCODE);

  /* reconfigure the problem */
  if (BVPD_CONFIG(theBVPDesc)==NULL)
  {
    PrintErrorMessage('E',"reinit","the problem has no reinit\n");
    return(CMDERRORCODE);
  }
  (*BVPD_CONFIG (theBVPDesc))(argc,argv);

  return(OKCODE);
}

/****************************************************************************/
/*D
   npexecute - execute a NumProc

   DESCRIPTION:
   This command executes a NumProc.
   It calls the function 'ExecuteNumProc'.

   'npexecute [<num proc name>] <argument list to be passed>'

   . <num~proc~name> - name of an existing NumProc

   SEE ALSO
   'numerics', 'NUM_PROC'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ExecuteNumProcCommand - Solve a numerical problem

   SYNOPSIS:
   static INT ExecuteNumProcCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function solves a numerical problem.
   It executes a num proc.

   npexecute [<num proc name>] <argument list to be passed>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ExecuteNumProcCommand (INT argc, char **argv)
{
  char theNumProcName[NAMESIZE];
  NUM_PROC *theNumProc;
  MULTIGRID *theMG;
  INT err;

  /* get NumProc */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" npexecute %",NAMELENSTR,"[ -~]")),theNumProcName)!=1) || (strlen(theNumProcName)==0))
  {
    theNumProc = GetCurrentNumProc();
    if (theNumProc == NULL)
    {
      PrintErrorMessage('E',"npexecute","there is no current numerical procedure");
      return (CMDERRORCODE);
    }
  }
  else
  {
    theNumProc = GetNumProcFromName(theNumProcName);
    if (theNumProc == NULL)
    {
      PrintErrorMessage('E',"npexecute","cannot find specified numerical procedure");
      return (CMDERRORCODE);
    }
  }

  theMG = GetCurrentMultigrid();
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"npexecute","there is no current multigrid\n");
    return (CMDERRORCODE);
  }

  if ((err=ExecuteNumProc(theNumProc,theMG,argc,argv))!=0)
  {
    sprintf(buffer,"execution of '%s' failed (error code %d)",theNumProcName,err);
    PrintErrorMessage('E',"npexecute",buffer);
    return (CMDERRORCODE);
  }

  return(OKCODE);
}


/****************************************************************************/
/*D
   nplist - list of all NumProcs created

   DESCRIPTION:
   This command list all NumProcs and displays their status.
   It calls the function 'ListNumProc'.

   'nplist'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   NumProcListCommand - List NumProc

   SYNOPSIS:
   static INT NumProcListCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists solvers.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT NumProcListCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);
  if (ListNumProc(GetCurrentNumProc())) return (CMDERRORCODE);
  return(OKCODE);
}

/****************************************************************************/
/*D
   npdisplay - display a NumProc

   DESCRIPTION:
   This command displays a NumProc.
   It calls the function 'DisplayNumProc'.

   'npdisplay [<num proc name>]'

   . <num~proc~name> - name of an existing NumProc (default is the current NumProc)

   D*/
/****************************************************************************/

/****************************************************************************/
/*
   NumProcDisplayCommand - display a NumProc

   SYNOPSIS:
   static INT NumProcDisplayCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This command displays current num proc.

   npdisplay

   . <num~proc~name> - name of an existing NumProc

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT NumProcDisplayCommand (INT argc, char **argv)
{
  char theNumProcName[NAMESIZE];
  NUM_PROC *theNumProc;

  if (argc!=1)
  {
    PrintErrorMessage('W',"npdisplay","don't specify arguments with npdisplay");
    return (OKCODE);
  }

  /* get NumProc */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" npdisplay %",NAMELENSTR,"[ -~]")),theNumProcName)!=1) || (strlen(theNumProcName)==0))
  {
    theNumProc = GetCurrentNumProc();
    if (theNumProc == NULL)
    {
      PrintErrorMessage('W',"npdisplay","there is no current numproc");
      return (OKCODE);
    }
  }
  else
  {
    theNumProc = GetNumProcFromName(theNumProcName);
    if (theNumProc == NULL)
    {
      PrintErrorMessage('W',"npdisplay","cannot find specified numerical procedure");
      return (OKCODE);
    }
  }

  /* display NumProc */
  if (DisplayNumProc(theNumProc))
  {
    PrintErrorMessage('E',"npdisplay","error during DisplayNumProc procedure\n");
    return (CMDERRORCODE);
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   npcreate - creating a NumProc

   DESCRIPTION:
   This command creates a NumProc.
   It calls the function 'CreateNumProc'.

   'npcreate <num proc name> $t <num proc type> $f <format>'

   . <num~proc~name> - name of the new NumProc
   . $t~<num~proc~type> - name of an existing NumProcType
   . $f~<format> - name of an existing format

   SEE ALSO
   'numerics', 'NUM_PROC',  'NP_TYPE_BASIC', 'NP_TYPE_SMOOTHER',
   'NP_TYPE_SOLVER', 'NP_TYPE_ASSEMBLE'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   NumProcCreateCommand - Create NumProc

   SYNOPSIS:
   static INT NumProcCreateCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function creates NumProc.
   It creates a num proc from a num proc type.

   npcreate <num proc name> $t <num proc type> $f <format>

   . <num~proc~name> - name of an existing NumProc

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT NumProcCreateCommand (INT argc, char **argv)
{
  NP_TYPE *theNumProcType;
  NUM_PROC *theNumProc;
  FORMAT *theFormat;
  char theNumProcTypeName[NAMESIZE], theFormatName[NAMESIZE], theNumProcName[NAMESIZE], buffer[128];
  INT i, topt, fopt;

  /* get NumProc name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" npcreate %",NAMELENSTR,"[ -~]")),theNumProcName)!=1) || (strlen(theNumProcName)==0))
  {
    PrintErrorMessage('E',"npcreate","specify the name of the theNumProcName to create");
    return (PARAMERRORCODE);
  }

  /* get format and NumProcType names */
  topt = fopt = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 't' :
      if (sscanf(argv[i],expandfmt(CONCAT3("t %",NAMELENSTR,"[ -~]")),theNumProcTypeName)!=1)
      {
        PrintHelp("npcreate",HELPITEM," (cannot read theNumProcType name)");
        return(PARAMERRORCODE);
      }
      topt = 1;
      break;
    case 'f' :
      if (sscanf(argv[i],expandfmt(CONCAT3("f %",NAMELENSTR,"[ -~]")),theFormatName)!=1)
      {
        PrintHelp("socreate",HELPITEM," (cannot read format name)");
        return(PARAMERRORCODE);
      }
      fopt = 1;
      break;
    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("npcreate",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }
  if (topt+fopt!=2)
  {
    PrintHelp("npcreate",HELPITEM," (specify n, t, f option)");
    return(PARAMERRORCODE);
  }

  /* find NumProcType and format */
  theNumProcType = GetNumProcTypeFromName(theNumProcTypeName);
  if (theNumProcType==NULL)
  {
    PrintHelp("npcreate",HELPITEM," (cannot find specified  NumProcType)");
    return(PARAMERRORCODE);
  }
  theFormat = GetFormat(theFormatName);
  if (theFormat==NULL)
  {
    PrintHelp("npcreate",HELPITEM," (cannot find specified  Format)");
    return(PARAMERRORCODE);
  }

  /* create NumProc */
  if ((theNumProc=CreateNumProc(theNumProcName,theNumProcType,theFormat))==NULL)
  {
    PrintHelp("npcreate",HELPITEM," (create solver failed)");
    return(PARAMERRORCODE);
  }

  if (SetCurrentNumProc(theNumProc))
    return(PARAMERRORCODE);

  return(OKCODE);
}

/****************************************************************************/
/*D
   npinit - inizialize a NumProc

   DESCRIPTION:
   This command inizializes a NumProc.
   It calls the function 'SetNumProc'.

   'npinit <argument list to be passed>'

   SEE ALSO:
   'numerics', 'NUM_PROC'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   NumProcInitCommand - Init NumProc

   SYNOPSIS:
   static INT NumProcInitCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function inits NumProc.
   It initialises current numproc.

   npinit <argument list to be passed>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT NumProcInitCommand (INT argc, char **argv)
{
  NUM_PROC *theNumProc;

  theNumProc = GetCurrentNumProc();
  if (theNumProc == NULL)
  {
    PrintErrorMessage('W',"npinit","there is no current NumProc");
    return (OKCODE);
  }

  if (SetNumProc(theNumProc,argc,argv))
  {
    PrintErrorMessage('W',"npinit","error during SetNumProc procedure");
    return (CMDERRORCODE);
  }


  return(OKCODE);
}

/****************************************************************************/
/*D
   scnp - make a NumProc the current NumProc

   DESCRIPTION:
   This command makes a NumProc the current NumProc.
   It sets current num proc.
   It calls the function 'GetNumProcFromName'.

   'scnp <num proc name>'

   . <num~proc~name> - name of an existing NumProc
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SetCurrentNumProcCommand - Make a NumProc the current NumProc

   SYNOPSIS:
   static INT SetCurrentNumProcCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function makes a NumProc the current NumProc.
   It sets current num proc.

   scnp <num proc name>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SetCurrentNumProcCommand (INT argc, char **argv)
{
  NUM_PROC *theNumProc;
  char npname[NAMESIZE];

  NO_OPTION_CHECK(argc,argv);

  /* get solver name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" scnp %",NAMELENSTR,"[ -~]")),npname)!=1)
  {
    PrintHelp("scnp",HELPITEM," (specify current NumProc name)");
    return(PARAMERRORCODE);
  }

  theNumProc = GetNumProcFromName(npname);

  if (theNumProc==NULL)
  {
    PrintErrorMessage('E',"scnp","no NumProc with this name exists");
    return (CMDERRORCODE);
  }

  if (SetCurrentNumProc(theNumProc))
    return (CMDERRORCODE);

  return(OKCODE);
}

/****************************************************************************/
/*D
   symlist - list contents of vector and matrix symbols

   DESCRIPTION:
   This command lists the contents of vector and matrix symbols.

   'scnp <num proc name>'

   . <num~proc~name> - name of an existing NumProc
   D*/
/****************************************************************************/

static INT SymListCommand (INT argc, char **argv)
{
  FORMAT *fmt;
  SYMBOL *sym;
  char format[NAMESIZE],name[NAMESIZE];
  INT found,res;

  if (argc!=2)
  {
    PrintErrorMessage('E',"symlist","specify one option with symlist");
    return (PARAMERRORCODE);
  }

  switch (argv[1][0])
  {
  case 's' :
    if ((res=sscanf(argv[1],"s %s %s",name,format))<1)
    {
      PrintErrorMessage('E',"symlist","specify symbol name with s option");
      return (PARAMERRORCODE);
    }
    if (res==1)
    {
      /* loop all formats */
      found = FALSE;
      for (fmt=GetFirstFormat(); fmt!=NULL; fmt=GetNextFormat(fmt))
        if ((sym=GetSymbol(ENVITEM_NAME(fmt),name))!=NULL)
        {
          found = TRUE;
          DisplaySymbol(sym);
        }
      if (!found) UserWriteF("no symbol '%s' in any format\n",name);
      return (OKCODE);
    }
    if (format[0]=='*')
      if (currMG!=NULL)
        strcpy(format,ENVITEM_NAME(MGFORMAT(currMG)));
      else
      {
        PrintErrorMessage('E',"symlist","no current multigrid");
        return (PARAMERRORCODE);
      }
    /* search in format */
    if ((sym=GetSymbol(format,name))!=NULL)
      DisplaySymbol(sym);
    else
      UserWriteF("no symbol '%s' in format '%s'%s\n",name,format,(format[0]=='*') ? " of current multigrid" : "");
    return (OKCODE);

  case 'V' :
    res = sscanf(argv[1],"V %s",format);
    if (res!=1)
    {
      /* loop all formats */
      for (fmt=GetFirstFormat(); fmt!=NULL; fmt=GetNextFormat(fmt))
      {
        UserWriteF("scanning format '%s'\n",ENVITEM_NAME(fmt));
        found = FALSE;
        for (sym=GetFirstSymbol(ENVITEM_NAME(fmt)); sym!=NULL; sym=GetNextSymbol(sym))
          if (SYM_IS_VEC(sym))
          {
            found = TRUE;
            DisplaySymbol(sym);
          }
        if (!found) UserWriteF("   >>> no vector symbol in format '%s'\n",ENVITEM_NAME(fmt));
      }
      return (OKCODE);
    }
    if (format[0]=='*')
      if (currMG!=NULL)
        strcpy(format,ENVITEM_NAME(MGFORMAT(currMG)));
      else
      {
        PrintErrorMessage('E',"symlist","no current multigrid");
        return (PARAMERRORCODE);
      }
    /* search in format */
    found = FALSE;
    for (sym=GetFirstSymbol(format); sym!=NULL; sym=GetNextSymbol(sym))
      if (SYM_IS_VEC(sym))
      {
        found = TRUE;
        DisplaySymbol(sym);
      }
    if (!found) UserWriteF("no vector symbol in format '%s'\n",format,(format[0]=='*') ? " of current multigrid" : "");
    return (OKCODE);

  case 'M' :
    res = sscanf(argv[1],"M %s",format);
    if (res!=1)
    {
      /* loop all formats */
      for (fmt=GetFirstFormat(); fmt!=NULL; fmt=GetNextFormat(fmt))
      {
        UserWriteF("scanning format '%s'\n",ENVITEM_NAME(fmt));
        found = FALSE;
        for (sym=GetFirstSymbol(ENVITEM_NAME(fmt)); sym!=NULL; sym=GetNextSymbol(sym))
          if (SYM_IS_MAT(sym))
          {
            found = TRUE;
            DisplaySymbol(sym);
          }
        if (!found) UserWriteF("   >>> no matrix symbol in format '%s'\n",ENVITEM_NAME(fmt));
      }
      return (OKCODE);
    }
    if (format[0]=='*')
      if (currMG!=NULL)
        strcpy(format,ENVITEM_NAME(MGFORMAT(currMG)));
      else
      {
        PrintErrorMessage('E',"symlist","no current multigrid");
        return (PARAMERRORCODE);
      }
    /* search in format */
    found = FALSE;
    for (sym=GetFirstSymbol(format); sym!=NULL; sym=GetNextSymbol(sym))
      if (SYM_IS_MAT(sym))
      {
        found = TRUE;
        DisplaySymbol(sym);
      }
    if (!found) UserWriteF("no matrix symbol in format '%s'\n",format,(format[0]=='*') ? " of current multigrid" : "");
    return (OKCODE);

  default :
    sprintf(buffer,"(invalid option '%s')",argv[1]);
    PrintHelp("symlist",HELPITEM,buffer);
    return (PARAMERRORCODE);
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
        newformat - init a format and allocate symbols

        SYNTAX:
   .vb
        newformat <format_name> [$V <vec_size>: {<vec_sym_name>}+ [$comp <comp_names> {$sub <sub_sym_name> <comps>}*]]
                            [$M <mat_size>: {<mat_sym_name>}+ [$comp <comp_names> {$sub <sub_sym_name> <matsizes+comps>}*]
                            [$d <mtype> <depth>]

      $V            # enrol storage for vector data: the specified size per type is multiplied
 # by the number of <vec_sym_name> specified. SYMBOLs are allocated with the
 # specified names, also scalar SYMBOLs are allocated

         <vec_size> := {<tp><size> }+, where <tp> is a char of
                                               n for NODE
                                               k for EDGE
                                               e for ELEM
                                               s for SIDE (3D only)
                       and <size> is a integer


      $comp         # optionally specify names for the components
         <comp_names> := one char per component of the symbol (without spaces)


      $sub          # if comp names have been specified optionally specify sub symbols for
 # each of the previously allocated symbols by sppecifying their comps

         <comps> := list of comps as specified in <comp_names>


          $M            # enrol storage for matrix data: the specified size per type is multiplied
 # by the number of <mat_sym_name> specified. SYMBOLs are allocated with the
 # specified names, also scalar SYMBOLs are allocated

         <mat_size> := {<rtp><rsize>x<ctp><csize> }+, where r and c refer to rows and columns of the matrix


      $comp         # optionally specify names for the components
         <comp_names> := TWO chars per component of the symbol (without spaces)


         <matsizes+comps> := {<rsize>x<csize> <comp_names>}+


      $d <mtype> <depth>  # set connection depth (default 0)

         <mtype> := <rtp>x<ctp>
   .ve

        DESCRIPTION:
        The 'newformat' command enrols a format for multigrid user data and creates symbols for the format.
        It also creates scalar and possibly other sub symbols.

        EXAMPLE:
   .vb
    newformat nsr                                                 # format name
              $V k2 e1: sol rhs tmp cor                           # 2 components per symbol and edge vector
 # 1 component per symbol and element vector
                        $comp uvp                                 # components are named u,v,p
                              $sub vel u v                        # create a sub symbol for each symbol with components u,v
              $M k2xk2 k2xe1 e1xk2: MAT                           # 2x2 components in edge-edge matrices
 # 2x1 components in edge-element matrices
 # 1x2 components in element-edge matrices
                        $comp uuuvvuvvupvppupv                    # components are named uu,uv,vu,vv...
                              $sub mom 2x2 uu uv vu vv 2x1 up vp  # create a sub symbol for each symbol with
 # components uu,uv,vu,vv in edge-edge matrices and
 # components up,vp in edge-element matrices
                              $sub velp 2x1 up vp
                              $sub off 2x1 up vp 1x2 pu pv
              $M k2xk2 k2xe1 e1xk2 e1xe1: LU
                        $comp uuuvvuvvupvppupvpp
              $d exe 1;                                           # set connection depth for element-element matrices 1
   .ve
        This format description will result in
   .       - 2*4 DOUBLEs in edge vectors
   .       - 1*4 DOUBLEs in element vectors
   .       - 4*1+4*1 DOUBLEs in edge-edge matrices
   .       - 2*1+2*1 DOUBLEs in edge-element matrices
   .       - 2*1+2*1 DOUBLEs in element-edge matrices
   .       - 1*1 DOUBLEs in element-element matrices
   and element-element matrices are connected with depth 1.

   The 'SYMBOL's created are
   .vb
   sol
   usol
   vsol
   psol
   velsol
   rhs
   urhs
   vrhs
   prhs
   velrhs
   tmp
   utmp
   vtmp
   ptmp
   veltmp
   cor
   ucor
   vcor
   pcor
   velcor
   MAT
   MATuu
   MATuv
   MATvu
   MATvv
   MATup
   MATvp
   MATpu
   MATpv
   MATmom
   MAToff
   LU
   LUuu
   LUuv
   LUvu
   LUvv
   LUup
   LUvp
   LUpu
   LUpv
   LUpp
   LUvelp
   .ve
    (To see this list enter 'ls /Formats/nsr' into the shell window after executing the above command.)

    SEE ALSO:
    'symlist'
   D*/
/****************************************************************************/

static INT CreateFormatCommand (INT argc, char **argv)
{
  INT err;

  err = CreateFormatCmd(argc,argv);

  switch (err)
  {
  case 0 : return (OKCODE);
  case 1 : PrintHelp("newformat",HELPITEM,NULL);
    return (PARAMERRORCODE);
  default : return (CMDERRORCODE);
  }
}

/****************************************************************************/
/*D
        setpf -  command to change current settings of the data
                                        listing functions of a format

        SYNTAX:
        setpf <format_name> [$V{0 | {+|-} {<vecsym_name>}+}]+ [$M{0 | {+|-} {<matsym_name>}+}]+

        DESCRIPTION:
        For a format previously enroled by use of the 'newformat' command the 'setpf'
        command specifies the symbols that are displayed when 'vmlist' is called with
        '$d' (list vector data) or '$d $m' (list vector and matrix data).

        EXAMPLE:
   .vb
        newformat ns $V n3: sol rhs tmp cor $M n3xn3: MAT LU;

        open grid $f ns $h 1000000;

 # print sol, rhs vector data and MAT, LU matrix data of vector with index 10
        setpf ns $V 0 $M 0 $V + sol rhs $M + MAT LU;
        vmlist $i 10 $m $d;

 # print sol vector data and MAT matrix data of vector with index 12
        setpf ns $V - rhs $M - LU;
        vmlist $i 12 $m $d;
   .ve

        SEE ALSO:
        showpf
   D*/
/****************************************************************************/

static INT SetPrintingFormatCommand (INT argc, char **argv)
{
  INT err;
  char format[NAMESIZE];

  if (sscanf(argv[0],"setpf %s",format)!=1)
  {
    if (currMG==NULL)
    {
      PrintErrorMessage('E',"setpf","no format specified and no current mg");
      return (PARAMERRORCODE);
    }
    strcpy(format,ENVITEM_NAME(MGFORMAT(currMG)));
  }

  err = SetPrintingFormatCmd(format,argc,argv);

  switch (err)
  {
  case 0 : return (OKCODE);
  case 1 : PrintHelp("setpf",HELPITEM,NULL);
    return (PARAMERRORCODE);
  default : return (CMDERRORCODE);
  }
}

/****************************************************************************/
/*D
        showpf - command to display current settings of data
                                        listing functions

        SYNTAX:
        showpf

        DESCRIPTION:
        ...

        SEE ALSO:
        setpf
   D*/
/****************************************************************************/

static INT ShowPrintingFormatCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  DisplayPrintingFormat();

  return (OKCODE);
}

/****************************************************************************/
/*D
   setkey - associate a command key with a ug command

   DESCRIPTION:
   This command associates a command key with a ug command.
   It calls the function 'SetCmdKey'.

   'setkey $<command char> $"<command sequence>"*'

   .  $k~<command~char>      - specifiy a single character which will be the command key
   .  $"<command~sequence>"  - give an arbitrary sequence of statements which is to be executed when the command key is pressed

   EXAMPLE:
   'setkey $r $"mark $a" $"refine";'

   Typing <alt> and then 'r' refines all elements on a UNIX system,
   <cmd>+<r> on a Macintosh.

   NOTICE:
   When the command sequence contains a @ character, the following token is interpreted
   as string variable and expanded as usual instantaneously when the command key is created.
   If the sequence contains a ? instead of the @, the string variable will be expanded with its current
   contents at execution time of the command sequence.

   The special character ? is limited however to the 'setkey' command.

   SEE ALSO:
   'setkey', 'delkey'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   SetCommandKeyCommand - associate a command key with a ug command

   SYNOPSIS:
   static INT SetCommandKeyCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function associates a command key with a ug command.
   It creates a command key and assign executable statements to it.

   setkey $<command char> $"<command sequence>"

   .  $k <command char>      - specifiy a single character which will be the command key
   .  $"<command sequence>"  - give an arbitrary sequence of statements which is to be
   .n                           executed when the command key is pressed

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT SetCommandKeyCommand (INT argc, char **argv)
{
  INT i,j,odd;
  char *ChatPtr;
  char CmdBuffer[INPUTBUFFERLEN];

  /* check input */
  if (argc < 3 )
    return(CMDERRORCODE);

  if (strlen(argv[1])!=1)
  {
    PrintErrorMessage('E',"setkey","only one character for cmd key");
    return (PARAMERRORCODE);
  }

  /* store input */
  ChatPtr = CmdBuffer;
  for (i=2; i<argc; i++)
  {
    *ChatPtr = '$';
    ChatPtr++;
    strcpy(ChatPtr,argv[i]);
    ChatPtr += strlen(argv[i]);
  }

  /* check input */
  if ((argv[2][0] != '\"') || (argv[argc-1][strlen(argv[argc-1])-1] != '\"'))
    return (CMDERRORCODE);
  j=0;
  for (i=0; i<strlen(CmdBuffer); i++)
    if (CmdBuffer[i] == '\"') j++;
  if (j%2 == 1)
    return (CMDERRORCODE);

  /* reconstruct command */
  odd = 0;
  for (i=0; i<strlen(CmdBuffer); i++)
  {
    if (CmdBuffer[i] != '\"') continue;
    odd = 1-odd;
    if (odd)
    {
      if (CmdBuffer[i-1] != '$')
        return (CMDERRORCODE);
      if (i==1)
        CmdBuffer[0] = ' ';
      else
        CmdBuffer[i-1] = ';';
    }
    CmdBuffer[i]   = ' ';
  }

  /* create cmd key */
  if (SetCmdKey(argv[1][0],CmdBuffer)!=0)
  {
    PrintErrorMessage('E',"setkey","cannot create cmd key");
    return(CMDERRORCODE);
  }

  return(OKCODE);
}

/****************************************************************************/
/*D
   delkey - delete an existing command key

   DESCRIPTION:
   This command deletes an existing command key.
   It calls the function 'DelCmdKey'.

   'delkey $all | $<command char>'

   .  $all                   - delete all command keys allocated before
   .  $<command~char>        - delete only the command key associated with <command char>

   SEE ALSO:
   'setkey', 'delkey'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   DeleteCommandKeyCommand - delete an existing command key

   SYNOPSIS:
   static INT DeleteCommandKeyCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function deletes an existing command key.
   It removes a command key allocated before.

   delkey $all | $<command char>

   .  $all                   - delete all command keys allocated before
   .  $<command char>        - delete only the command key associated with <command char>

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT DeleteCommandKeyCommand (INT argc, char **argv)
{
  /* check input */
  if (argc!=2)
  {
    PrintHelp("delkey",HELPITEM," (give exactly one argument)");
    return(CMDERRORCODE);
  }

  /* remove key(s) */
  if (strcmp(argv[1],"all")==0)
  {
    if (DelAllCmdKeys()!=0)
    {
      PrintErrorMessage('E',"delkey","failed deleting all cmd keys");
      return(CMDERRORCODE);
    }
  }
  else
  {
    if (DelCmdKey(argv[1][0])!=0)
    {
      PrintErrorMessage('E',"delkey","failed deleting cmd key");
      return(CMDERRORCODE);
    }
  }

  return (OKCODE);
}

/****************************************************************************/
/*D
   key - list all existing command keys

   DESCRIPTION:
   This command lists all existing command keys.
   It calls the function 'ListCmdKeys'.

   'key'

   SEE ALSO:
   'setkey', 'delkey'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   ListCommandKeysCommand - List all existing command keys

   SYNOPSIS:
   static INT ListCommandKeysCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function lists all existing command keys.
   It lists all command keys defined together with their contents.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ListCommandKeysCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  ListCmdKeys();

  return (OKCODE);
}

/****************************************************************************/
/*D
   refreshon - sets the refresh state on

   DESCRIPTION:
   This command sets the refresh state on: the pictures on the screen
   device will be updated instantaneously.

   'refreshon'

    SEE ALSO:
        'InvalidatePicturesOfMG', 'InvalidateUgWindowsOfMG'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   RefreshOnCommand - Set refresh status TRUE

   SYNOPSIS:
   static INT RefreshOnCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function sets refresh status TRUE.
   It pictures on the screen device will be updated instantaneously.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT RefreshOnCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  SetRefreshState(ON);
  return(OKCODE);
}

/****************************************************************************/
/*D
   refreshoff - sets the refresh state off

   DESCRIPTION:
   This command sets the refresh state off.

   'refreshoff'

   SEE ALSO:
   'refreshon'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   RefreshOffCommand - Set refresh status FALSE

   SYNOPSIS:
   static INT RefreshOffCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function sets refresh status FALSE.
   It pictures on the screen device will be updated ONLY explicitely.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/


static INT RefreshOffCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  SetRefreshState(OFF);
  return(OKCODE);
}

/****************************************************************************/
/*
   InitFindRange -

   SYNOPSIS:
   static INT InitFindRange ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function defines ':findrange'.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InitFindRange ()
{
  /* install a struct for findrange */
  if (MakeStruct(":findrange")!=0)
    return (1);
  return (0);
}

/****************************************************************************/
/*
   InitScreenSize -

   SYNOPSIS:
   static INT InitScreenSize ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function prints the size in pixels of the screen (if there) ==> max window size.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InitScreenSize ()
{
  /* install a struct for screensize */
  if (MakeStruct(":screensize")!=0)
    return (1);
  return (0);
}

#ifdef __THREEDIM__

/****************************************************************************/
/*D
   checkparity - check parity of elements (3d only)

   DESCRIPTION:
   This command sets the refresh status FALSE.
   It calls the function 'CheckParityOfElements'.

   'checkparity'
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   CheckParityCommand - Set refresh status FALSE

   SYNOPSIS:
   static INT CheckParityCommand (INT argc, char **argv);

   PARAMETERS:
   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   DESCRIPTION:
   This function sets refresh status FALSE.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT CheckParityCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

  /* current multigrid */
  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"checkparity","no current multigrid\n");
    return (CMDERRORCODE);
  }

  if (CheckParityOfElements(theMG)) return (CMDERRORCODE);

  return(OKCODE);
}

#endif

/****************************************************************************/
/*D
   InitCommands - Initialization of the commands

   SYNOPSIS:
   INT InitCommands ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function does initialization of all ug-commands, using
   'CreateCommand'.
   It initializes the 'clock', 'findrange' and 'screensize' comand.

   SEE ALSO:
   commands

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if error occured.
   D*/
/****************************************************************************/

INT InitCommands ()
{
  /* general commands */
  if (CreateCommand("quit",                       QuitCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("help",                       HelpCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("checkhelp",          CheckHelpCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("readclock",          ReadClockCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("resetclock",         ResetClockCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("date",                       DateCommand                                     )==NULL) return (__LINE__);

  /* commands for environement management */
  if (CreateCommand("cd",                         ChangeEnvCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("ls",                         ListEnvCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("pwd",                        PrintEnvDirCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("envinfo",            EnvInfoCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("set",                        SetCommand                                              )==NULL) return (__LINE__);
  if (CreateCommand("dv",                         DeleteVariableCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("ms",                         MakeStructCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("cs",                         ChangeStructCommand                             )==NULL) return (__LINE__);
  if (CreateCommand("pws",                        PrintWorkStructCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("ds",                         DeleteStructCommand                     )==NULL) return (__LINE__);


  /* commands for protocol and logfile output */
  if (CreateCommand("protoOn",            ProtoOnCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("protoOff",           ProtoOffCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("protocol",           ProtocolCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("logon",                      LogOnCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("logoff",             LogOffCommand                                   )==NULL) return (__LINE__);


  /* commands for grid management */
  if (CreateCommand("setcurrmg",          SetCurrentMultigridCommand              )==NULL) return (__LINE__);
  if (CreateCommand("new",                        NewCommand                                              )==NULL) return (__LINE__);
  if (CreateCommand("open",                       OpenCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("close",                      CloseCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("save",                       SaveCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("level",                      LevelCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("renumber",           RenumberMGCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("smooth",             SmoothMGCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("ordernodes",         OrderNodesCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("lexorderv",          LexOrderVectorsCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("orderv",             OrderVectorsCommand                     )==NULL) return (__LINE__);
  if (CreateCommand("shellorderv",        ShellOrderVectorsCommand                )==NULL) return (__LINE__);
  if (CreateCommand("setindex",           SetIndexCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("extracon",           ExtraConnectionCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("check",                      CheckCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("in",                         InsertInnerNodeCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("bn",                         InsertBoundaryNodeCommand               )==NULL) return (__LINE__);
  if (CreateCommand("deln",                       DeleteNodeCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("move",                       MoveNodeCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("ie",                         InsertElementCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("dele",                       DeleteElementCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("refine",             RefineCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("mark",                       MarkCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("find",                       FindCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("select",             SelectCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("wplist",             WindowPictureListCommand                )==NULL) return (__LINE__);
  if (CreateCommand("mglist",             MGListCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("glist",                      GListCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("nlist",                      NListCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("elist",                      EListCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("slist",                      SelectionListCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("vmlist",             VMListCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("quality",            QualityCommand                                  )==NULL) return (__LINE__);
    #ifdef __TWODIM__
  if (CreateCommand("bnodes",                 BnodesCommand                                       )==NULL) return (__LINE__);
  if (CreateCommand("makegrid",           MakeGridCommand                                 )==NULL) return (__LINE__);
    #endif

  /* commands for window and picture management */
  if (CreateCommand("screensize",         ScreenSizeCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("openwindow",         OpenWindowCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("closewindow",        CloseWindowCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("setcurrwindow",      SetCurrentWindowCommand                 )==NULL) return (__LINE__);
  if (CreateCommand("drawtext",           DrawTextCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("openpicture",        OpenPictureCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("closepicture",       ClosePictureCommand                     )==NULL) return (__LINE__);
  if (CreateCommand("clearpicture",       ClearPictureCommand                     )==NULL) return (__LINE__);
  if (CreateCommand("picframe",           PicFrameCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("setcurrpicture", SetCurrentPictureCommand            )==NULL) return (__LINE__);
  if (CreateCommand("setview",            SetViewCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("vdisplay",           DisplayViewCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("walk",                       WalkCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("walkaround",         WalkAroundCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("zoom",                       ZoomCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("drag",                       DragCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("rotate",             RotateCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("setplotobject",      SetPlotObjectCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("polist",             PlotObjectListCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("plot",                       PlotCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("findrange",          FindRangeCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("updateDoc",          UpdateDocumentCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("cmfn",                       CreateMetafileNameCommand               )==NULL) return (__LINE__);

  /* commands for problem management */
  if (CreateCommand("clear",                      ClearCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("reinit",             ReInitCommand                                   )==NULL) return (__LINE__);

  /* commands for NumProc management */
  if (CreateCommand("npexecute",          ExecuteNumProcCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("nplist",                     NumProcListCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("npdisplay",          NumProcDisplayCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("npcreate",           NumProcCreateCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("npinit",                     NumProcInitCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("scnp",                       SetCurrentNumProcCommand                )==NULL) return (__LINE__);

  /* symbols */
  if (CreateCommand("symlist",            SymListCommand                                  )==NULL) return (__LINE__);

  /* formats */
  if (CreateCommand("newformat",          CreateFormatCommand                     )==NULL) return (__LINE__);
  if (CreateCommand("setpf",                      SetPrintingFormatCommand        )==NULL) return (__LINE__);
  if (CreateCommand("showpf",             ShowPrintingFormatCommand       )==NULL) return (__LINE__);

  /* miscellaneous commands */
  if (CreateCommand("setkey",             SetCommandKeyCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("delkey",             DeleteCommandKeyCommand                 )==NULL) return (__LINE__);
  if (CreateCommand("keylist",            ListCommandKeysCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("refreshon",          RefreshOnCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("refreshoff",         RefreshOffCommand                               )==NULL) return (__LINE__);
        #ifdef __THREEDIM__
  if (CreateCommand("checkparity",        CheckParityCommand                              )==NULL) return (__LINE__);
        #endif

  if (InitClock()                 !=0) return (__LINE__);
  if (InitFindRange()     !=0) return (__LINE__);
  if (InitScreenSize()    !=0) return (__LINE__);

  return(0);
}
