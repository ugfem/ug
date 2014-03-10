// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file commands.c
 * \ingroup ui
 */

/** \addtogroup ui
 *
 * @{
 */

/****************************************************************************/
/*                                                                          */
/* File:      commands.c                                                    */
/*                                                                          */
/* Purpose:   definition of all dimension independent commands of ug        */
/*                                                                          */
/* Author:    Henrik Rentz-Reichert                                         */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   29.01.92 begin, ug version 2.0                                */
/*            02.02.95 begin, ug version 3.0                                */
/*                                                                          */
/* Remarks:   for dimension dependent commands see commands2d/3d.c          */
/*                                                                          */
/****************************************************************************/

#ifdef __MPW32__
#pragma segment cmd
#endif

/****************************************************************************/
/*                                                                          */
/* defines to exclude functions                                             */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/* system include files                                                     */
/* application include files                                                */
/*                                                                          */
/****************************************************************************/

/* standard C library */
#include "config.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/* low module */
#include "ugtypes.h"
#include "architecture.h"
#include "ugtime.h"
#include "initug.h"
#include "defaults.h"
#include "misc.h"
#include "ugstruct.h"
#include "fileopen.h"
#include "ugenv.h"
#include "debug.h"
#include "heaps.h"              /* for MEM declaration */
#include "general.h"

/* devices module */
#include "ugdevices.h"

/* grid manager module */
#include "gm.h"
#include "elements.h"
#include "cw.h"
#include "pargm.h"
#include "rm.h"
#include "evm.h"
#include "ugm.h"
#include "algebra.h"
#include "shapes.h"
#include "mgio.h"

/* grid generator module */
#ifdef __TWODIM__
#include "gm/gg2/ggm.h"
#include "gm/gg2/ggmain.h"
#endif
#if defined __THREEDIM__ && defined _NETGEN
#include "gm/gg3/gg3d.h"
#endif

/* numerics module */
#include "np.h"
#include "ugblas.h"
#include "formats.h"
#include "disctools.h"
#include "data_io.h"
#include "npcheck.h"
#include "udm.h"
#include "fvgeom.h"

/* graph module */
#include "wpm.h"
#include "wop.h"
#include "graph.h"
#include "graphics/grape/connectuggrape.h"
#ifdef _COVISE
#include "graphics/covise/coviseif.h"
#endif

/* user interface module */
#include "uginterface.h"
#include "ugstruct.h"
#include "cmdint.h"
#include "cmdline.h"
#include "helpmsg.h"
#include "mmio.h"
#include "dio.h"

#ifdef ModelP
#include "parallel.h"
#endif

#include "ppif_namespace.h"

/* own header */
#include "commands.h"


USING_UG_NAMESPACES
  USING_PPIF_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*    compile time constants defining static data size (i.e. arrays)        */
/*    other constants                                                       */
/*    macros                                                                */
/*                                                                          */
/****************************************************************************/

/** \brief Size of the general purpose text buffer*/
#define BUFFERSIZE                              512

#define WHITESPACE                              " \t"

/** \brief Size of some strings */
#define LONGSTRSIZE                     256
/** \brief Length of some strings                                       */
#define LONGSTRLEN                              255
/** \brief LONGSTRLEN as string                                 */
#define LONGSTRLENSTR                   "255"

/** @name For ProtoOnCommand */
/*@{*/
#define NORENAME_PROTO                  0
#define APPEND_PROTO                    1
#define RENAME_PROTO                    2
#define TRYRENAME_PROTO                 3
#define MAXPATHLENGTH                   255
#define MAXRENAMECHAR                   'z'
/*@}*/

/** @name For the .list commands */
/*@{*/
#define DO_ID                                   1
#define DO_SELECTION                    2
#define DO_ALL                                  3
/*@}*/

/** @name For MarkCommand */
/*@{*/
#define MARK_ALL                                1
#define AI_MARK_ALL             256
#define MARK_COARSEN                    2
#define MARK_ID                                 3
#define MARK_SELECTION                  4
#define NO_SIDE_SPECIFIED               -1
#define NO_RULE_SPECIFIED               -1
#define NO_OF_RULES                     64
/*@}*/

/* for anisotropic refinement */

/* for save command */
#define NO_COMMENT                               "no comment"

/** @name For array commands */
/*@{*/
#define AR_NVAR_MAX                     10
#define AR_NVAR(p)                      ((p)->nVar)
#define AR_VARDIM(p,i)          ((p)->VarDim[i])
#define AR_DATA(p,i)            ((p)->data[i])
/*@}*/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/* in the corresponding include file!)                                      */
/*                                                                          */
/****************************************************************************/

struct MarkRule
{
  const char *RuleName;         /*!< what you type in the mark cmdline	*/
  INT RuleId;           /*!< corresponding rule ID for refine   */
};

typedef struct MarkRule MARKRULE;

typedef struct
{
  /** \brief Fields for environment directory */
  ENVVAR v;

  /* data */
  INT nVar;
  INT VarDim[AR_NVAR_MAX];
  DOUBLE data[1];

} ARRAY;

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static MULTIGRID *currMG=NULL;                  /*!< The current multigrid			*/

static NP_BASE *currNumProc=NULL;               /*!< Current numerical procedure		*/

static char buffer[BUFFERSIZE];         /*!< General purpose text buffer		*/

static FILE     *protocolFile=NULL;     /*!< For protocol commands			*/

/** \brief Name and ID of available rules	*/
static MARKRULE myMR[NO_OF_RULES]=     {{"red",        RED},
                                        {"no",         NO_REFINEMENT},
                                        {"blue",       BLUE},
                                        {"copy",       COPY},
                                         #ifdef __TWODIM__
                                        {"bi_1",       BISECTION_1},
                                        {"bi_2q",      BISECTION_2_Q},
                                        {"bi_2t1", BISECTION_2_T1},
                                        {"bi_2t2", BISECTION_2_T2},
                                        {"bi_3",       BISECTION_3},
                                         #endif
                                         #ifdef __THREEDIM__
                                        /* rules for tetrahedra */
                                         #ifndef TET_RULESET
                                        {"tet2hex",TETRA_RED_HEX},
                                        {"pri2hex",PRISM_RED_HEX},
                                         #endif
                                        /* rules for prisms */
                                        {"pri_quadsect",PRISM_QUADSECT},
                                        {"pri_bisect_hex0",PRISM_BISECT_HEX0},
                                        {"pri_bisect_hex1",PRISM_BISECT_HEX1},
                                        {"pri_bisect_hex2",PRISM_BISECT_HEX2},
                                        {"pri_rot_l",PRISM_ROTATE_LEFT},
                                        {"pri_rot_r",PRISM_ROTATE_RGHT},
                                        {"pri_quadsect_eins",PRISM_QUADSECT_HEXPRI0},
                                        /* rules for tetrahedra */
                                        {"hex_bisect_eins",HEX_BISECT_0_1},
                                        {"hex_bisect_zwei",HEX_BISECT_0_2},
                                        {"hex_bisect_drei",HEX_BISECT_0_3},
                                        {"hex_trisect_eins",HEX_TRISECT_0},
                                        {"hex_trisect_fuenf",HEX_TRISECT_5},
                                        {"hex_quadsect_null",HEX_QUADSECT_0},
                                        {"hex_quadsect_eins",HEX_QUADSECT_1},
                                        {"hex_quadsect_zwei",HEX_QUADSECT_2},
                                        {"hex_bisect_vier",HEX_BISECT_HEXPRI0},
                                        {"hex_bisect_fuenf",HEX_BISECT_HEXPRI1},
                                         #endif
                                        {"coarse", COARSE}};

static DOUBLE Time0;                            /*!< Time offset for readclock		*/

static char userPath[1024];             /*!< Environment path for ls,cd		*/

static INT untitledCounter=0;   /*!< Counter for untitled multigrids	*/

/** @name Some variables to transfer info between QualityCommand and QualityElement*/
/*@{*/
static DOUBLE min,max,themin,themax,minangle,maxangle;
static INT lessopt,greateropt,selectopt;
static char mintext[32],maxtext[32],minmaxtext[32];
/*@}*/

/** @name Counters for windows and pictures */
/*@{*/
static INT wincounter=1;
static INT piccounter=1;
/*@}*/

/** @name Stuff for the array commands */
/*@{*/
static INT theArrayDirID;
static INT theArrayVarID;
static INT arraypathes_set=FALSE;
/*@}*/

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

#if defined(CAD) && defined(__THREEDIM__)
MULTIGRID *ConvertCADGrid  (char *theFormat, char *CADOutputFileName,unsigned long heapSize);
#endif


/****************************************************************************/
/** \brief Return a pointer to the current multigrid
 *
 * This function returns a pionter to the current multigrid.
 *
 * @return <ul>
 *    <li> pointer to multigrid </li>
 *    <li> NULL if there is no current multigrid. </li>
 * </ul>
 */
/****************************************************************************/

MULTIGRID * NS_DIM_PREFIX GetCurrentMultigrid (void)
{
  return (currMG);
}

/****************************************************************************/
/** \brief Set the current multigrid if it is valid
 *
 * @param theMG pointer to multigrid
 *
 * This function sets the current multigrid if it is valid, i. e.
 * the function checks whether 'theMG' acually points to a multigrid.
 * It can be NULL only if no multigrid is open.
 *
 * @result <ul>
 *    <li> 0 if ok </li>
 *    <li> 1 if theMG is not in the multigrid list </li>
 * </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX SetCurrentMultigrid (MULTIGRID *theMG)
{
  MULTIGRID *mg;

  if (ResetPrintingFormat())
    REP_ERR_RETURN(CMDERRORCODE);

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


/** \brief Implementation of \ref quit. */
static INT QuitCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);
  SetDoneFlag();
  return(QUITCODE);
}



/** \brief Implementation of \ref exitug. */
static INT ExitUgCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  ExitUg();

  exit(0);

  return 0;
}


/** \brief Implementation of \ref help. */
static INT HelpCommand (INT argc, char **argv)
{
  INT i,res,mode,rv;
  COMMAND *Cmd;
  char buf[NAMESIZE];

        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

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
      UserWrite("no help found\n"
                "maybe a command matches...\n");
      Cmd = SearchUgCmd(buf);
      if (Cmd!=NULL)
        rv = PrintHelp(ENVITEM_NAME(Cmd),mode,NULL);
    }
  }
  else
    rv = PrintHelp("help",HELPITEM,NULL);

  switch (rv)
  {
  case HELP_OK :
    return (OKCODE);

  case HELP_NOT_FOUND :
    UserWriteF(" no help entry found for '%s'\n",buf);
    return (OKCODE);

  default :
    PrintErrorMessage('E',"help","(unknown)");
  }

  return (CMDERRORCODE);
}


/** \brief Implementation of \ref checkhelp. */
static INT CheckHelpCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

  CheckHelp();

  return (OKCODE);
}

/** \brief Implementation of \ref cmfn. */
static INT CreateMetafileNameCommand (INT argc, char **argv)
{
  INT res,i,nopt;
  int frame;
  char name[LONGSTRSIZE],varname[LONGSTRSIZE];
  char fullname[LONGSTRSIZE];
  char *ext;

  nopt = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'n' :
      if (sscanf(argv[i],expandfmt(CONCAT3("n %",NAMELENSTR,"[ -~]")),varname)!=1)
      {
        PrintErrorMessage('E',"cmfn","can't read varname");
        return(PARAMERRORCODE);
      }
      nopt = 1;
      break;
    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      return (PARAMERRORCODE);
    }

  res = sscanf(argv[0],expandfmt(CONCAT3(" cmfn %",LONGSTRLENSTR,"[0-9:.a-zA-Z_] %255[0-9:.a-zA-Z_]")),name,buffer);
  if (res!=2) return(CMDERRORCODE);

  if (GetStringValueInt(buffer,&frame)) return(CMDERRORCODE);
  ext = GetStringVar("EXT");

  if (ext==NULL)
    sprintf(fullname,"%s.%04d",name,frame);
  else
    sprintf(fullname,"%s.%04d.%s",name,frame,ext);

  if (nopt)
  {
    if (SetStringVar(varname,fullname)) return(CMDERRORCODE);
  }
  else
  {
    if (SetStringVar(name,fullname)) return(CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref readclock. */
static INT ReadClockCommand (INT argc, char **argv)
{
  DOUBLE Time;

  NO_OPTION_CHECK(argc,argv);

  Time = ARCH_DIFF_TIMER(CURRENT_TIME_LONG,Time0);

  if (SetStringValue(":CLOCK",Time)!=0) {
    PrintErrorMessage('E',"readclock",
                      "could not get string variable :CLOCK");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}

/** \brief Implementation of \ref resetclock. */
static INT ResetClockCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  Time0 = CURRENT_TIME_LONG;

  return (OKCODE);
}

/****************************************************************************/
/** \brief Starting the time mesuring
 *
 * This function starts the time mesuring.
 * It sets the global variable 'Time0' to zero.
 *
 * @return <ul>
 *    <li> 0 if ok </li>
 *    <li> 1 if error occured. </li>
 * </ul>
 */
/****************************************************************************/

static INT InitClock(void)
{
  Time0 = CURRENT_TIME_LONG;

  return(0);
}


/** \brief Implementation of \ref date. */
static INT DateCommand (INT argc, char **argv)
{
  time_t Time;
  const char *fmt;
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


/** \brief Implementation of \ref ls. */
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


/** \brief Implementation of \ref cd. */
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
  if (strlen(buffer)==0)
  {
    /* empty path: change to root directory */
    strcpy(userPath,DIRSEP);
    currentDir = ChangeEnvDir(userPath);
    if (currentDir == NULL)
      return (CMDERRORCODE);
    else
      return (OKCODE);
  }
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



/** \brief Implementation of \ref pwd. */
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


/** \brief Implementation of \ref envinfo. */
static INT EnvInfoCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  EnvHeapInfo(buffer);

  UserWrite(buffer);

  return (OKCODE);
}


/** \brief Implementation of \ref set. */
static INT SetCommand (INT argc, char **argv)
{
  char name[LONGSTRSIZE], *namePtr;
  INT flag,ropt,i,res,rv;
        #ifdef ModelP
  /* allocate buffer for definition of long string variables */
  char *buffer;

  /* alloc command buffer */
  if ((buffer=(char *)malloc(cmdintbufsize))==NULL)
  {
    PrintErrorMessage('F',"SetCommand","could not allocate buffer");
    return(__LINE__);
  }

  /* if there are problems with expandfmt print format string !
          printf("SetCommand(): formatstring=%s\n",
                  expandfmt(CONCAT3(" set %",LONGSTRLENSTR,"[0-9:.a-zA-Z_] %16384[]\t\n -~]")));
   */
  res = sscanf(argv[0],expandfmt(CONCAT3(" set %",LONGSTRLENSTR,"[0-9:.a-zA-Z_] %16384[]\t\n -~]")),name,buffer);
        #ifdef __CC__
  res = sscanf(argv[0],expandfmt(CONCAT3(" set %",LONGSTRLENSTR,"[0-9:.a-zA-Z_] %16384[^]\t\n -~]")),name,buffer);
        #endif
        #else
  res = sscanf(argv[0],expandfmt(CONCAT3(" set %",LONGSTRLENSTR,"[0-9:.a-zA-Z_] %255[ -~]")),name,buffer);
        #endif /* ModelP */

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

        #ifdef ModelP
  free(buffer);
        #endif

  if (rv==0)
    return (OKCODE);
  else
    return (CMDERRORCODE);
}

/** \brief Implementation of \ref dv. */
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


/** \brief Implementation of \ref ms. */
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


/** \brief Implementation of \ref cs. */
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


/** \brief Implementation of \ref pws. */
static INT PrintWorkStructCommand (INT argc, char **argv)
{
  char structPath[1024];

  NO_OPTION_CHECK(argc,argv);

  GetStructPathName(structPath, 1024);
  UserWrite(structPath);
  UserWrite("\n");

  return(OKCODE);
}


/** \brief Implementation of \ref ds. */
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


#define PROTOCOL_SEP '%'
/** \brief Implementation of \ref protocol. */
static INT ProtocolCommand (INT argc, char **argv)
{
  INT i,from;

        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

  if (protocolFile==NULL)
  {
    PrintErrorMessage('E',"protocol","no protocol file open!");
    return (CMDERRORCODE);
  }

  for (i=1; i<argc; i++)
  {
    if (argv[i][0]!=PROTOCOL_SEP)
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
    while ((i+1<argc) && (argv[i+1][0]!=PROTOCOL_SEP))
      fprintf(protocolFile," $%s",(argv[++i]));
  }

  return (OKCODE);
}

/****************************************************************************/
/** \brief Open protocol file where specially formatted output is saved
 *
 * @param name - name of the protocol file
 * @param mode - APPEND_PROTO, RENAME_PROTO, TRYRENAME_PROTO, or NORENAME_PROTO
 *
 * This function opens protocol file where specially formatted output is saved.
 * It opens a protocol file for specially formatted output to file.
 *
 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error occured. </li>
 * </ul>
 */
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
        PrintErrorMessageF('E',"OpenProto","could't find a new name for '%s'",fullname);
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
    PrintErrorMessageF('W',"OpenProto","opened protcol file '%s' (instead of '%s')",realname+pathlen,name);
  }

  return(0);
}


/** \brief Implementation of \ref protoOn. */
static INT ProtoOnCommand (INT argc, char **argv)
{
  static char protoFileName[NAMESIZE];
  INT res,i,RenameMode;

        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

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


/** \brief Implementation of \ref protoOff. */
static INT ProtoOffCommand (INT argc, char **argv)
{
        #ifdef ModelP
  if (me != master) return(OKCODE);
        #endif

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

/****************************************************************************/
/** \brief Return pointer to current protocol file
 *
 * This function returns a pointer to the current protocol file (NULL if not open).
 *
 * @return <ul>
 *    <li> file ptr if ok </li>
 *    <li> NULL if no protocol file open </li>
 * </ul>
 */
/****************************************************************************/

FILE* NS_DIM_PREFIX GetProtocolFile ()
{
  return (protocolFile);
}


/** \brief Implementation of \ref logon. */
static INT LogOnCommand (INT argc, char **argv)
{
  char logfile[NAMESIZE];
  INT i,rv,popt,res,pext,meext,rename;
  int ropt;

  /* check options */
  popt = pext = meext = rename = FALSE;
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

    case 'e' :
                                #ifdef ModelP
      pext = TRUE;
                                #endif
      break;

    case 'a' :
                                #ifdef ModelP
      meext = TRUE;
                                #endif
      break;

    case 'f' :
      CloseLogFile();
      break;

    case 'r' :
      res = sscanf(argv[i]," r %d",&ropt);
      if (res==0 || (res==1 && ropt==1)) rename = TRUE;
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
        #ifdef ModelP
  if (pext == TRUE)
  {
    sprintf(logfile,"%s.p%04d",logfile,procs);
  }
  if (meext == TRUE)
  {
    sprintf(logfile,"%s.%04d",logfile,me);
  }
  else if (me != master)
    return (OKCODE);
        #endif

  rv = OpenLogFile(logfile,rename);
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


/** \brief Implementation of \ref logoff. */
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

#ifdef __TWODIM__


/** \brief Implementation of \ref cnom. */
static INT CnomCommand (INT argc, char **argv)
{
  char docName[32],plotprocName[NAMESIZE],tagName[NAMESIZE];
  int i,flag;

  if (currMG==NULL)
  {
    PrintErrorMessage('E',"cnom","no multigrid active");
    return(CMDERRORCODE);
  }

  /* get document name */
  docName[0] = (char) 0;
  sscanf(argv[0]," cnom %31[ -~]",docName);
  if (strlen(docName)==0)
  {
    PrintErrorMessage('E',"cnom","no document name");
    return(PARAMERRORCODE);
  }
  if (argc!=2)
  {
    PrintErrorMessage('E',"cnom","specify only one argument with cnom");
    PrintHelp("cnom",HELPITEM,buffer);
    return(PARAMERRORCODE);
  }

  flag=0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'p' :
      if (sscanf(argv[i],expandfmt(CONCAT3("p %",NAMELENSTR,"[ -~]")),plotprocName)!=1)
      {
        PrintErrorMessage('E',"cnom","can't read plotprocName");
        return(PARAMERRORCODE);
      }
      flag=flag|1;
      break;
    case 't' :
      if (sscanf(argv[i],expandfmt(CONCAT3("t %",NAMELENSTR,"[ -~]")),tagName)!=1)
      {
        PrintErrorMessage('E',"cnom","can't read tagName");
        return(PARAMERRORCODE);
      }
      flag=flag|2;
      break;
    default :
      flag=flag|4;
      break;
    }

  if (flag!=3)
  {
    PrintHelp("cnom",HELPITEM,buffer);
    return (PARAMERRORCODE);
  }


  return(SaveCnomGridAndValues(currMG,docName,plotprocName,tagName));
}
#endif


/** \brief Implementation of \ref configure. */
INT NS_DIM_PREFIX ConfigureCommand (INT argc, char **argv)
{
  BVP *theBVP;
  BVP_DESC theBVPDesc;
  char BVPName[NAMESIZE];

  /* get BVP name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" configure %",NAMELENSTR,"[ -~]")),BVPName)!=1) || (strlen(BVPName)==0))
  {
    PrintHelp("configure",HELPITEM," (cannot read BndValProblem specification)");
    return(PARAMERRORCODE);
  }

  theBVP = BVP_GetByName(BVPName);
  if (theBVP == NULL)
  {
    PrintHelp("configure",HELPITEM," (cannot read BndValProblem specification)");
    return(PARAMERRORCODE);
  }

  if (BVP_SetBVPDesc(theBVP,&theBVPDesc))
    return (CMDERRORCODE);

  if (BVPD_CONFIG(&theBVPDesc)!=NULL)
    if ((*BVPD_CONFIG(&theBVPDesc))(argc,argv))
    {
      PrintErrorMessage('E',"configure"," (could not configure BVP)");
      return(CMDERRORCODE);
    }

  return(OKCODE);
}


/** \brief Implementation of \ref close. */
static INT CloseCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  UGWINDOW *theWin;
  PICTURE *thePic,*NextPic,*currPic;
  INT i,closeonlyfirst;

  if (ResetPrintingFormat())
    REP_ERR_RETURN(CMDERRORCODE);

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


/** \brief Implementation of \ref new. */
INT NS_DIM_PREFIX NewCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char Multigrid[NAMESIZE],BVPName[NAMESIZE],Format[NAMESIZE];
  MEM heapSize;
  INT i,bopt,fopt,hopt,IEopt,emptyGrid;

  /* get multigrid name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" new %",NAMELENSTR,"[ -~]")),Multigrid)!=1) || (strlen(Multigrid)==0))
    sprintf(Multigrid,"untitled-%d",(int)untitledCounter++);

  theMG = GetMultigrid(Multigrid);
  if ((theMG != NULL) && (theMG == currMG)) CloseCommand(0,NULL);

  /* get problem, domain and format */
  heapSize = 0;
  bopt = fopt = hopt = FALSE;
  IEopt = TRUE;
  emptyGrid = FALSE;
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

    case 'n' :
      IEopt = FALSE;
      break;

    case 'e' :
      emptyGrid = TRUE;
      break;

    case 'h' :
      if (ReadMemSizeFromString(argv[i]+1,&heapSize)!=0)           /* skip leading 'h' in argv */
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
  theMG = CreateMultiGrid(Multigrid,BVPName,Format,heapSize,IEopt,!emptyGrid);
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"new","could not create multigrid");
    return(CMDERRORCODE);
  }

  currMG = theMG;

  return(OKCODE);
}


/** \brief Implementation of \ref open. */
static INT OpenCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char Multigrid[NAMESIZE],File[NAMESIZE],BVPName[NAMESIZE],Format[NAMESIZE],type[NAMESIZE];
  char *theBVP,*theFormat,*theMGName;
  MEM heapSize;
  INT i,force,IEopt,autosave,try_load,fqn,mgpathes_set_old;

  /* get multigrid name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" open %",NAMELENSTR,"[ -~]")),File)!=1) || (strlen(File)==0))
  {
    PrintErrorMessage('E',"open","specify the name of the file to open");
    return (PARAMERRORCODE);
  }

  /* get problem and format */
  strcpy(type,"asc");
  theBVP = theFormat = theMGName = NULL;
  heapSize = force = autosave = fqn = 0;
  try_load = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      autosave = 1;
      break;

    case 'b' :
      if (sscanf(argv[i],expandfmt(CONCAT3("b %",NAMELENSTR,"[ -~]")),BVPName)!=1)
      {
        PrintHelp("open",HELPITEM," (cannot read BndValProblem specification)");
        return(PARAMERRORCODE);
      }
      theBVP = BVPName;
      break;

    case 'f' :
      if (sscanf(argv[i],expandfmt(CONCAT3("f %",NAMELENSTR,"[ -~]")),Format)!=1)
      {
        PrintHelp("open",HELPITEM," (cannot read format specification)");
        return(PARAMERRORCODE);
      }
      theFormat = Format;
      break;

    case 'F' :
      force = 1;
      break;

    case 'n' :
      IEopt = FALSE;
      break;

    case 'm' :
      if (sscanf(argv[i],expandfmt(CONCAT3("m %",NAMELENSTR,"[ -~]")),Multigrid)!=1)
      {
        PrintHelp("open",HELPITEM," (cannot read multigrid specification)");
        return(PARAMERRORCODE);
      }
      theMGName = Multigrid;
      break;

    case 't' :
      if (strncmp(argv[i],"try",3)==0)
      {
        try_load = TRUE;
        break;
      }
      if (sscanf(argv[i],expandfmt(CONCAT3("t %",NAMELENSTR,"[ -~]")),type)!=1)
      {
        PrintHelp("open",HELPITEM," (cannot read type specification)");
        return(PARAMERRORCODE);
      }
      break;

    case 'h' :
      if (ReadMemSizeFromString(argv[i]+1,&heapSize)!=0)                           /* skip leading 'h' in argv */
      {
        PrintHelp("open",HELPITEM," (cannot read heapsize specification)");
        return(PARAMERRORCODE);
      }
      break;

    case 'z' :
      fqn = 1;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("open",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (fqn) {
    mgpathes_set_old = mgpathes_set;
    mgpathes_set = 0;
  }

  /* allocate the multigrid structure */
  theMG = LoadMultiGrid(theMGName,File,type,theBVP,theFormat,
                        heapSize,force,IEopt,autosave);

  if (fqn) mgpathes_set = mgpathes_set_old;

  if (theMG==NULL)
  {
    PrintErrorMessage('E',"open","could not open multigrid");
    if (try_load)
      return(CMDERRORCODE);
    else
      RETURN(CMDERRORCODE);
  }
  currMG = theMG;

  return(OKCODE);
}


/** \brief Implementation of \ref save. */
static INT SaveCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char Name[NAMESIZE],type[NAMESIZE],Comment[LONGSTRSIZE];
  INT i,autosave,rename,res;
  int ropt;

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
  autosave=rename=0;
  strcpy(Comment,NO_COMMENT);
  strcpy(type,"asc");
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

    case 't' :
      if (sscanf(argv[i],expandfmt(CONCAT3("t %",NAMELENSTR,"[ -~]")),type)!=1)
      {
        PrintHelp("open",HELPITEM," (cannot read type specification)");
        return(PARAMERRORCODE);
      }
      break;

    case 'a' :
      autosave=1;
      break;

    case 'r' :
      res = sscanf(argv[i]," r %d",&ropt);
      if (res==0 || (res==1 && ropt==1)) rename = 1;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("save",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (SaveMultiGrid(theMG,Name,type,Comment,autosave,rename)) return (CMDERRORCODE);

  return(OKCODE);
}


/** \brief Implementation of \ref savedomain. */
static INT SaveDomainCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char Name[NAMESIZE];

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"savedomain","no open multigrid");
    return (CMDERRORCODE);
  }

  /* scan name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" savedomain %",NAMELENSTR,"[ -~]")),Name)!=1)
    strcpy(Name,BVPD_NAME(MG_BVPD(theMG)));

  if (BVP_Save(MG_BVP(theMG),Name,ENVITEM_NAME(theMG),MGHEAP(theMG),argc,argv)) return (CMDERRORCODE);

  return(OKCODE);
}



#ifdef ModelP
static INT comm_comp, comm_box_comp;

static int Gather_Scalar (DDD_OBJ obj, void *data)
{
  NODE *pv = (NODE *)obj;
  DOUBLE *d = (DOUBLE *)data;

  d[0] = VVALUE(NVECTOR(pv),comm_comp);
  d[1] = VVALUE(NVECTOR(pv),comm_box_comp);

  return (NUM_OK);
}

static int Scatter_Scalar (DDD_OBJ obj, void *data)
{
  NODE *pv = (NODE *)obj;
  DOUBLE *d = (DOUBLE *)data;

  VVALUE(NVECTOR(pv),comm_comp) += d[0];
  VVALUE(NVECTOR(pv),comm_box_comp) += d[1];

  return (NUM_OK);
}

INT l_sum_scalar (GRID *g, INT comp, INT box_comp)
{
  comm_comp = comp; comm_box_comp = box_comp;
  DDD_IFAExchange(BorderNodeSymmIF, GRID_ATTR(g), 2*sizeof(DOUBLE),
                  Gather_Scalar, Scatter_Scalar);
  return (NUM_OK);
}

static int Gather_Vector (DDD_OBJ obj, void *data)
{
  NODE *pv = (NODE *)obj;
  DOUBLE *d = (DOUBLE *)data;
  INT i;

  for (i=0; i<DIM; i++)
    d[i] = VVALUE(NVECTOR(pv),comm_comp+i);
  d[DIM] = VVALUE(NVECTOR(pv),comm_box_comp);

  return (NUM_OK);
}

static int Scatter_Vector (DDD_OBJ obj, void *data)
{
  NODE *pv = (NODE *)obj;
  DOUBLE *d = (DOUBLE *)data;
  INT i;

  for (i=0; i<DIM; i++)
    VVALUE(NVECTOR(pv),comm_comp+i) += d[i];
  VVALUE(NVECTOR(pv),comm_box_comp) += d[DIM];

  return (NUM_OK);
}

INT l_sum_vector (GRID *g, INT comp, INT box_comp)
{
  comm_comp = comp; comm_box_comp = box_comp;
  DDD_IFAExchange(BorderNodeSymmIF, GRID_ATTR(g), (DIM+1)*sizeof(DOUBLE),
                  Gather_Vector, Scatter_Vector);
  return (NUM_OK);
}
#endif


INT AverageScalar (MULTIGRID *mg, EVALUES *eval, char *eval_name, VECDATA_DESC *vecdesc)
{
  INT box_compo,compo,n,k,i,error,j;
  GRID *g;
  NODE *theNode;
  VECTOR *theVector;
  SHORT NCmpInType[NVECTYPES];
  VECDATA_DESC *box_vd=NULL;
  PreprocessingProcPtr preProc;
  ElementEvalProcPtr evalProc;
  ELEMENT *el;
  DOUBLE value;
  DOUBLE *CornersCoord[MAX_CORNERS_OF_ELEM];
  DOUBLE LocalCoord[DIM];
  DOUBLE local[DIM];
  FVElementGeometry Geo;
  SubControlVolume *scv;

  /* get component */
  compo = VD_ncmp_cmpptr_of_otype(vecdesc,NODEVEC,&n)[0];
  assert(n>0);

  /* reset average to zero */
  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    g = GRID_ON_LEVEL(mg,k);
    for (theNode=FIRSTNODE(g); theNode!= NULL; theNode=SUCCN(theNode))
    {
      theVector = NVECTOR(theNode);
      VVALUE(theVector,compo) = 0.0;
    }
  }

  /* allocate vecdata desc for box volume */
  for (k=0; k<NVECTYPES; k++) NCmpInType[k]=0;
  NCmpInType[NODEVEC]=1;       /* one component in node */
  if (AllocVDfromNCmp(mg,0,TOPLEVEL(mg),NCmpInType,NULL,&box_vd)) return(1);

  /* get component of box volume */
  box_compo = VD_ncmp_cmpptr_of_otype(box_vd,NODEVEC,&n)[0];

  /* reset box volume to zero */
  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    g = GRID_ON_LEVEL(mg,k);
    for (theNode=FIRSTNODE(g); theNode!= NULL; theNode=SUCCN(theNode))
    {
      theVector = NVECTOR(theNode);
      VVALUE(theVector,box_compo) = 0.0;
    }
  }

  /* prepare plot function */
  preProc = eval->PreprocessProc;
  if (preProc!=NULL) preProc(eval_name,mg);
  evalProc = eval->EvalProc;

  /* run through all levels, elements, corners, ... */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      /* get control volume geometry */
      error = EvaluateFVGeometry(el,&Geo);

      /* process all corners of element */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        /* sub control volume */
        scv = Geo.scv+i;

        /* get value in corner i */
        for (j=0; j<CORNERS_OF_ELEM(el); j++)
          CornersCoord[j] = CVECT(MYVERTEX(CORNER(el,j)));                       /* x,y,z of corners */
        LocalCornerCoordinates(DIM,TAG(el),i,local);
        for (j=0; j<DIM; j++) LocalCoord[j] = local[j];
        value = evalProc(el,(const DOUBLE **)CornersCoord,LocalCoord);

        /* accumulate value weighted with box volume */
        theVector = NVECTOR(CORNER(el,i));
        VVALUE(theVector,compo) += scv->volume*value;

        /* ... and accumulate box volume */
        VVALUE(theVector,box_compo) += scv->volume;
      }
    }

  /* if it is parallel, then sum values at interfaces */
        #ifdef ModelP
  for (k=0; k<=TOPLEVEL(mg); k++)
    l_sum_scalar(GRID_ON_LEVEL(mg,k),compo,box_compo);
        #endif

  /* finally, divide by box volume */
  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    g = GRID_ON_LEVEL(mg,k);
    for (theNode=FIRSTNODE(g); theNode!= NULL; theNode=SUCCN(theNode))
    {
      theVector = NVECTOR(theNode);
      VVALUE(theVector,compo) =
        VVALUE(theVector,compo)/VVALUE(theVector,box_compo);
    }
  }

  /* and free box volume */
  FreeVD(mg,0,TOPLEVEL(mg),box_vd);

  return(0);
}


INT AverageVector (MULTIGRID *mg, EVECTOR *eval, char *eval_name, VECDATA_DESC *vecdesc)
{
  INT box_compo,compo,n,k,i,error,j;
  GRID *g;
  NODE *theNode;
  VECTOR *theVector;
  SHORT NCmpInType[NVECTYPES];
  VECDATA_DESC *box_vd=NULL;
  PreprocessingProcPtr preProc;
  ElementVectorProcPtr evalProc;
  ELEMENT *el;
  DOUBLE value[DIM];
  DOUBLE *CornersCoord[MAX_CORNERS_OF_ELEM];
  DOUBLE LocalCoord[DIM];
  DOUBLE local[DIM];
  FVElementGeometry Geo;
  SubControlVolume *scv;

  /* get component */
  compo = VD_ncmp_cmpptr_of_otype(vecdesc,NODEVEC,&n)[0];
  assert(n==DIM);
  for (j=1; j<DIM; j++)
    if (VD_ncmp_cmpptr_of_otype(vecdesc,NODEVEC,&n)[j]!=compo+j)
    {
      UserWrite("can only handle consecutive components!\n");
      return(1);
    }

  /* reset average to zero */
  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    g = GRID_ON_LEVEL(mg,k);
    for (theNode=FIRSTNODE(g); theNode!= NULL; theNode=SUCCN(theNode))
    {
      theVector = NVECTOR(theNode);
      for (j=0; j<DIM; j++) VVALUE(theVector,compo+j) = 0.0;
    }
  }

  /* allocate vecdata desc for box volume */
  for (k=0; k<NVECTYPES; k++) NCmpInType[k]=0;
  NCmpInType[NODEVEC]=1;       /* one component in node */
  if (AllocVDfromNCmp(mg,0,TOPLEVEL(mg),NCmpInType,NULL,&box_vd)) return(1);

  /* get component of box volume */
  box_compo = VD_ncmp_cmpptr_of_otype(box_vd,NODEVEC,&n)[0];

  /* reset box volume to zero */
  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    g = GRID_ON_LEVEL(mg,k);
    for (theNode=FIRSTNODE(g); theNode!= NULL; theNode=SUCCN(theNode))
    {
      theVector = NVECTOR(theNode);
      VVALUE(theVector,box_compo) = 0.0;
    }
  }

  /* prepare plot function */
  preProc = eval->PreprocessProc;
  if (preProc!=NULL) preProc(eval_name,mg);
  evalProc = eval->EvalProc;

  /* run through all levels, elements, corners, ... */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      /* get control volume geometry */
      error = EvaluateFVGeometry(el,&Geo);

      /* process all corners of element */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        /* sub control volume */
        scv = Geo.scv+i;

        /* get value in corner i */
        for (j=0; j<CORNERS_OF_ELEM(el); j++)
          CornersCoord[j] = CVECT(MYVERTEX(CORNER(el,j)));                       /* x,y,z of corners */
        LocalCornerCoordinates(DIM,TAG(el),i,local);
        for (j=0; j<DIM; j++) LocalCoord[j] = local[j];
        evalProc(el,(const DOUBLE **)CornersCoord,LocalCoord,value);

        /* accumulate value weighted with box volume */
        theVector = NVECTOR(CORNER(el,i));
        for (j=0; j<DIM; j++) VVALUE(theVector,compo+j) += scv->volume*value[j];

        /* ... and accumulate box volume */
        VVALUE(theVector,box_compo) += scv->volume;
      }
    }

  /* if it is parallel, then sum values at interfaces */
        #ifdef ModelP
  for (k=0; k<=TOPLEVEL(mg); k++)
    l_sum_vector(GRID_ON_LEVEL(mg,k),compo,box_compo);
        #endif

  /* finally, divide by box volume */
  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    g = GRID_ON_LEVEL(mg,k);
    for (theNode=FIRSTNODE(g); theNode!= NULL; theNode=SUCCN(theNode))
    {
      theVector = NVECTOR(theNode);
      for (j=0; j<DIM; j++)
        VVALUE(theVector,compo+j) =
          VVALUE(theVector,compo+j)/VVALUE(theVector,box_compo);
    }
  }

  /* and free box volume */
  FreeVD(mg,0,TOPLEVEL(mg),box_vd);

  return(0);
}

#define MAXVARIABLES    10
/** \brief Implementation of \ref average. */
static INT AverageCommand (INT argc, char **argv)
{
  MULTIGRID *mg;
  INT ns;
  INT nv;
  EVALUES *es[MAXVARIABLES];
  char es_name[MAXVARIABLES][NAMESIZE];
  EVECTOR *ev[MAXVARIABLES];
  char ev_name[MAXVARIABLES][NAMESIZE];
  char s[NAMESIZE];
  VECDATA_DESC *vd,*vd2;
  SHORT NCmpInType[NVECTYPES];
  INT i,k;

  /* get current multigrid */
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"average","no multigrid open\n");
    return (OKCODE);
  }

  /* scan input for $ns \<scalar eval proc> [$s \<symbol>] or $nv \<scalar eval proc> [$s \<symbol>]*/
  ns=nv=0;       /* no plot functions yet */
  for(i=1; i<argc; i++)
  {
    if (strncmp(argv[i],"ns",2)==0) {
      if (ns>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"average:","too many scalar variables specified\n");
        break;
      }
      sscanf(argv[i],"ns %s", s);
      es[ns] = GetElementValueEvalProc(s);
      if (es[ns]==NULL)
      {
        PrintErrorMessageF('E',"average:","could not find scalar eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(es_name[ns],s);
        i++;
      }
      else
        strcpy(es_name[ns],es[ns]->v.name);
      ns++;
      continue;
    }

    if (strncmp(argv[i],"nv",2)==0) {
      if (nv>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"average:","too many vector variables specified\n");
        break;
      }
      sscanf(argv[i],"nv %s", s);
      ev[nv] = GetElementVectorEvalProc(s);
      if (ev[nv]==NULL)
      {
        PrintErrorMessageF('E',"average:","could not find vector eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(ev_name[nv],s);
        i++;
      }
      else
        strcpy(ev_name[nv],ev[nv]->v.name);
      nv++;
      continue;
    }

  }

  /* scalar functions with one component in node .. */
  for (k=0; k<NVECTYPES; k++) NCmpInType[k]=0;
  NCmpInType[NODEVEC]=1;
  for (i=0; i<ns; i++)
  {
    /* allocate free VECDATA_DESC */
    vd=NULL;
    if (AllocVDfromNCmp(mg,0,TOPLEVEL(mg),NCmpInType,NULL,&vd)) return(1);

    /* if name exist then it must be that of vd */
    vd2=GetVecDataDescByName(mg,es[i]->v.name);
    if (vd2!=vd && vd2!=NULL)
    {
      UserWrite(es[i]->v.name); UserWrite(": name exists already, skipping\n");
      FreeVD(mg,0,TOPLEVEL(mg),vd);
      return(1);
      continue;
    }

    /* rename vd to name of plot proc */
    strcpy(vd->v.name,es[i]->v.name);
    UserWrite(es[i]->v.name); UserWrite(": created\n");

    /* ... and average it */
    if (AverageScalar(mg,es[i],es_name[i],vd)) return(1);
  }

  /* vector functions with DIM components in node .. */
  for (k=0; k<NVECTYPES; k++) NCmpInType[k]=0;
  NCmpInType[NODEVEC]=DIM;
  for (i=0; i<nv; i++)
  {
    /* allocate free VECDATA_DESC */
    vd=NULL;
    if (AllocVDfromNCmp(mg,0,TOPLEVEL(mg),NCmpInType,NULL,&vd)) return(1);

    /* if name exist then it must be that of vd */
    vd2=GetVecDataDescByName(mg,ev[i]->v.name);
    if (vd2!=vd && vd2!=NULL)
    {
      UserWrite(ev[i]->v.name); UserWrite(": name exists already, skipping\n");
      FreeVD(mg,0,TOPLEVEL(mg),vd);
      return(1);
      continue;
    }

    /* rename vd to name of plot proc */
    strcpy(vd->v.name,ev[i]->v.name);
    UserWrite(ev[i]->v.name); UserWrite(": created\n");

    /* ... and average it */
    if (AverageVector(mg,ev[i],ev_name[i],vd)) return(1);
  }

  return(0);
}

/** \brief Implementation of \ref freeaverage. */
static INT FreeAverageCommand (INT argc, char **argv)
{
  MULTIGRID *mg;
  INT ns;
  INT nv;
  EVALUES *es[MAXVARIABLES];
  char es_name[MAXVARIABLES][NAMESIZE];
  EVECTOR *ev[MAXVARIABLES];
  char ev_name[MAXVARIABLES][NAMESIZE];
  char s[NAMESIZE];
  VECDATA_DESC *vd;
  INT i;

  /* get current multigrid */
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"average","no multigrid open\n");
    return (OKCODE);
  }

  /* scan input for $ns <scalar eval proc> [$s <symbol>] or $nv <scalar eval proc> [$s <symbol>]*/
  ns=nv=0;       /* no plot functions yet */
  for(i=1; i<argc; i++)
  {
    if (strncmp(argv[i],"ns",2)==0) {
      if (ns>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"freeaverage:","too many scalar variables specified\n");
        break;
      }
      sscanf(argv[i],"ns %s", s);
      es[ns] = GetElementValueEvalProc(s);
      if (es[ns]==NULL)
      {
        PrintErrorMessageF('E',"freeaverage:","could not find scalar eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(es_name[ns],s);
        i++;
      }
      else
        strcpy(es_name[ns],es[ns]->v.name);

      /* get VECDATA_DESC */
      vd = GetVecDataDescByName(mg,es[ns]->v.name);
      if (vd==NULL)
      {
        UserWrite(es[ns]->v.name); UserWrite(": VECDATA_DESC not found\n");
        continue;
      }

      /* free it */
      FreeVD(mg,0,TOPLEVEL(mg),vd);
      UserWrite(vd->v.name); UserWrite(": freed\n");

      ns++;
      continue;
    }

    if (strncmp(argv[i],"nv",2)==0) {
      if (nv>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"freeaverage:","too many vector variables specified\n");
        break;
      }
      sscanf(argv[i],"nv %s", s);
      ev[nv] = GetElementVectorEvalProc(s);
      if (ev[nv]==NULL)
      {
        PrintErrorMessageF('E',"freeaverage:","could not find vector eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(ev_name[nv],s);
        i++;
      }
      else
        strcpy(ev_name[nv],ev[nv]->v.name);

      /* get VECDATA_DESC */
      vd = GetVecDataDescByName(mg,ev[nv]->v.name);
      if (vd==NULL)
      {
        UserWrite(ev[nv]->v.name); UserWrite(": VECDATA_DESC not found\n");
        continue;
      }

      /* free it */
      FreeVD(mg,0,TOPLEVEL(mg),vd);
      UserWrite(vd->v.name); UserWrite(": freed\n");

      nv++;
      continue;
    }

  }

  return(0);
}



static INT ReadSaveDataInput (MULTIGRID *theMG, INT argc, char **argv, const char *VDSym, char EvalChar, VECDATA_DESC **theVD, EVALUES **theEVal, EVECTOR **theEVec)
{
  INT i;

  /* init */
  *theVD = NULL; *theEVal = NULL; *theEVec = NULL;

  /* read VecDesc */
  for (i=1; i<argc; i++)
    if (argv[i][0]==VDSym[0])
    {
      if (sscanf(argv[i]+1," %s",buffer)!=1) break;
      if (strlen(buffer)>=NAMESIZE) break;
      *theVD = GetVecDataDescByName(theMG,buffer);
      if (*theVD!=NULL) return (1);
    }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    if (argv[i][0]==EvalChar)
    {
      if (sscanf(argv[i]+1," %s",buffer)!=1) break;
      if (strlen(buffer)>=NAMESIZE) break;
      *theEVal = GetElementValueEvalProc(buffer);
      if (*theEVal!=NULL) return (2);
      *theEVec = GetElementVectorEvalProc(buffer);
      if (*theEVec!=NULL) return (3);
    }

  return (0);
}

#define NM_MAX          100
/** \brief Implementation of \ref savedata. */
static INT SaveDataCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char FileName[NAMESIZE],type[NAMESIZE],mname[NAMESIZE];
  VECDATA_DESC *theVDList[NM_MAX];
  char NameList[5][NAMESIZE];
  char **Names;
  char *NamePtr[5];
  EVALUES *theEValues[5];
  EVECTOR *theEVector[5];
  INT i,j,n,ret,number,rename,nm,save_without_mg;
  int iValue;
  double Value[3];
  DOUBLE t[3];

  theMG = currMG;
  if (theMG==NULL) {PrintErrorMessage('E',"savedata","no open multigrid"); return (CMDERRORCODE);}

  /* scan filename */
  if (sscanf(argv[0],expandfmt(CONCAT3(" savedata %",NAMELENSTR,"[ -~]")),FileName)!=1) { PrintErrorMessage('E',"save","cannot read filename"); return (CMDERRORCODE);}

  strcpy(type,"asc");
  number = -1;
  t[0]=t[1]=t[2]=-1.0; nm=0;
  rename=0;
  save_without_mg=0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'm' :
      ret=sscanf(argv[i]+1," %s %d",mname,&iValue);
      if (ret!=2)
      {
        PrintHelp("savedata",HELPITEM," (multiple vector specification)");
        return(PARAMERRORCODE);
      }
      nm=iValue;
      if (nm<1 || nm>NM_MAX)
      {
        PrintHelp("savedata",HELPITEM," (multiple vector number out of range [0,xxx])");
        return(PARAMERRORCODE);
      }
      break;
    case 't' :
      if (sscanf(argv[i],expandfmt(CONCAT3("t %",NAMELENSTR,"[ -~]")),type)!=1)
      {
        PrintHelp("savedata",HELPITEM," (cannot read type specification)");
        return(PARAMERRORCODE);
      }
      break;

    case 'n' :
      if (sscanf(argv[i],"n %d",&iValue)!=1)
      {
        PrintHelp("savedata",HELPITEM," (cannot read number specification)");
        return(PARAMERRORCODE);
      }
      number = iValue;
      if (number<0 || number > 999999)
      {
        PrintHelp("savedata",HELPITEM," (number out of range [0,9999999])");
        return(PARAMERRORCODE);
      }
      break;

    case 'T' :
      ret = sscanf(argv[i],"T %lf %lf %lf",Value,Value+1,Value+2);
      if (ret<1 || ret>3)
      {
        PrintHelp("savedata",HELPITEM," (cannot read TIME specification)");
        return(PARAMERRORCODE);
      }
      for (j=0; j<ret; j++)
        t[j] = Value[j];
      if (t[0]<0.0)
      {
        PrintHelp("savedata",HELPITEM," (TIME out of range ]-inf, 0.0[)");
        return(PARAMERRORCODE);
      }
      break;
    case 'r' :
      ret = sscanf(argv[i]," r %d",&iValue);
      if (ret==0 || (ret==1 && iValue==1)) rename = 1;
      break;
    case 'p' :
      save_without_mg=1;
      break;
    }
  if (((t[0]<0.0) && number>=0) || ((t[0]>=0.0) && number<0))
  {
    PrintHelp("savedata",HELPITEM," (specify both or none the options 'n' and 'T')");
    return(PARAMERRORCODE);
  }

  /* get input */
  n=0;
  if (nm<=0)
  {
    ret = ReadSaveDataInput (theMG,argc,argv,"a",'A',theVDList+0,theEValues+0,theEVector+0);        if (ret) n++;
    ret = ReadSaveDataInput (theMG,argc,argv,"b",'B',theVDList+1,theEValues+1,theEVector+1);        if (ret) n++;
    ret = ReadSaveDataInput (theMG,argc,argv,"c",'C',theVDList+2,theEValues+2,theEVector+2);        if (ret) n++;
    ret = ReadSaveDataInput (theMG,argc,argv,"d",'D',theVDList+3,theEValues+3,theEVector+3);        if (ret) n++;
    ret = ReadSaveDataInput (theMG,argc,argv,"e",'E',theVDList+4,theEValues+4,theEVector+4);        if (ret) n++;
  }
  else
  {
    for (i=0; i<nm; i++)
    {
      sprintf(buffer,"%s%d",mname,i);
      theVDList[i]=GetVecDataDescByName(theMG,buffer);
      if (theVDList[i]==NULL) return (PARAMERRORCODE);
    }
    n=nm;
  }

  /* read names */
  Names=NULL;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'N' :
      if (sscanf(argv[i],"N %s %s %s %s %s",NameList[0],NameList[1],NameList[2],NameList[3],NameList[4])!=n) return (PARAMERRORCODE);
      Names=NamePtr;
      for (j=0; j<5; j++) NamePtr[j]=NameList[j];
      break;
    }

  if (n<=0) return (PARAMERRORCODE);
  if (SaveData(theMG,FileName,rename,save_without_mg,type,number,t[0],t[1],t[2],n,theVDList,theEValues,theEVector,Names)) return (PARAMERRORCODE);

  return(OKCODE);
}


/** \brief Implementation of \ref loaddata. */
static INT LoadDataCommand (INT argc, char **argv)
{
  char FileName[NAMESIZE],type[NAMESIZE],mname[NAMESIZE];
  VECDATA_DESC *theVDList[NM_MAX];
  INT i,m,n,number,force,ret,nm,fqn,datapathes_set_old;
  int iValue;
  INT read_in_currmg=0;
  MEM heapSize;

  /* scan filename */
  if (sscanf(argv[0],expandfmt(CONCAT3(" loaddata %",NAMELENSTR,"[ -~]")),FileName)!=1) { PrintErrorMessage('E',"save","cannot read filename"); return (CMDERRORCODE);}

  strcpy(type,"asc");
  heapSize=0;
  force=0;
  fqn=0;
  number = -1;
  nm=0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'm' :
      ret=sscanf(argv[i]+1," %s %d",mname,&iValue);
      if (ret!=2)
      {
        PrintHelp("savedata",HELPITEM," (multiple vector specification)");
        return(PARAMERRORCODE);
      }
      nm=iValue;
      if (nm<1 || nm>NM_MAX)
      {
        PrintHelp("savedata",HELPITEM," (multiple vector number out of range [0,xxx])");
        return(PARAMERRORCODE);
      }
      break;

    case 't' :
      if (sscanf(argv[i],expandfmt(CONCAT3("t %",NAMELENSTR,"[ -~]")),type)!=1)
      {
        PrintHelp("loaddata",HELPITEM," (cannot read type specification)");
        return(PARAMERRORCODE);
      }
      break;

    case 'f' :
      force=1;
      break;

    case 'n' :
      if (sscanf(argv[i],"n %d",&iValue)!=1)
      {
        PrintHelp("loaddata",HELPITEM," (cannot read number specification)");
        return(PARAMERRORCODE);
      }
      number = iValue;
      if (number<0 || number > 999999)
      {
        PrintHelp("loaddata",HELPITEM," (number out of range [0,999999])");
        return(PARAMERRORCODE);
      }
      break;

    case 'h' :
      if (ReadMemSizeFromString(argv[i]+1,&heapSize)!=0)                           /* skip leading 'h' in argv */
      {
        PrintHelp("new",HELPITEM," (cannot read heapsize specification)");
        return(PARAMERRORCODE);
      }
      break;

    case 'r' :
      read_in_currmg=1;
      break;

    case 'z' :
      fqn = 1;
      break;
    }

  /* consistency */
  if (read_in_currmg) force=0;

  if (fqn) {
    datapathes_set_old = datapathes_set;
    datapathes_set = 0;
  }

  /* open or reopen multigrid */
  if (force)
  {
    currMG = OpenMGFromDataFile(currMG,number,type,FileName,heapSize);
    if (currMG==NULL)
    {
      PrintErrorMessage('E',"loaddata","cannot open multigrid");
      return (CMDERRORCODE);
    }
  }
  if (currMG==NULL) {PrintErrorMessage('E',"loaddata","no open multigrid"); return (CMDERRORCODE);}

  /* get vecdatadesc */
  n=0;
  if (nm==0)
  {
    theVDList[n++]=ReadArgvVecDesc(currMG,"a",argc,argv);
    theVDList[n++]=ReadArgvVecDesc(currMG,"b",argc,argv);
    theVDList[n++]=ReadArgvVecDesc(currMG,"c",argc,argv);
    theVDList[n++]=ReadArgvVecDesc(currMG,"d",argc,argv);
    theVDList[n++]=ReadArgvVecDesc(currMG,"e",argc,argv);
  }
  else
  {
    for (i=0; i<nm; i++)
    {
      sprintf(buffer,"%s%d",mname,i);
      theVDList[i]=GetVecDataDescByName(currMG,buffer);
      if (theVDList[i]==NULL)
      {
        theVDList[i]=CreateVecDescOfTemplate (currMG,buffer,NULL);
        if (theVDList[i]==NULL) return (CMDERRORCODE);
      }
    }
    n=nm;
  }
  m=0;
  for (i=0; i<n; i++)
    if (theVDList[i]!=NULL)
      m = i+1;

  if (m<=0) return (PARAMERRORCODE);
  if (read_in_currmg) {
    /* renumber multigrid corresponding to saved order */
    if (RenumberMultiGrid(currMG,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0)!=GM_OK) {
      PrintErrorMessage('E',"loaddata","renumbering of the mg failed");
      return (CMDERRORCODE);
    }
  }
  if (LoadData(currMG,FileName,type,number,m,theVDList)) return (CMDERRORCODE);

  if (fqn) datapathes_set = datapathes_set_old;

  return(OKCODE);
}


/** \brief Implementation of \ref changemc. */
static INT ChangeMagicCookieCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  int iValue;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"changemc","no open multigrid");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0]," changemc %d",&iValue)!=1)
  {
    PrintErrorMessage('E',"changemc","cannot read magic-cookie");
    return (CMDERRORCODE);
  }
  MG_MAGIC_COOKIE(theMG) = iValue;
  return(OKCODE);
}


/** \brief Implementation of \ref level. */
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
    if ((l<BOTTOMLEVEL(theMG)) || (l>TOPLEVEL(theMG)))
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
    if (CURRENTLEVEL(theMG)==BOTTOMLEVEL(theMG))
    {
      PrintErrorMessage('W',"level","already on BOTTOMLEVEL");
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

  UserWriteF("  current level is %d (bottom level %d, top level %d)\n",
             CURRENTLEVEL(theMG),BOTTOMLEVEL(theMG),TOPLEVEL(theMG));

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return(OKCODE);
}


/** \brief Implementation of \ref renumber. */
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

  if (RenumberMultiGrid(theMG,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0)!=GM_OK)
  {
    PrintErrorMessage('E',"renumber","renumbering of the mg failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref wplist. */
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


/** \brief Implementation of \ref mglist. */
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

  longformat = TRUE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 's' :
      longformat = FALSE;
      break;

    case 'l' :
      /* old syntax */
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


/** \brief Implementation of \ref glist. */
static INT GListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: GListCommand(): me not in Context,"\
                     " no listing of grid\n",me))
    return(OKCODE);
  }
        #endif

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


/** \brief Implementation of \ref nlist. */
static INT NListCommand (INT argc, char **argv)
{

  MULTIGRID *theMG;
  INT i,fromN,toN,res,mode,idopt,dataopt,boundaryopt,neighbouropt,verboseopt;
  char buff[32];

  /* following variables: keep type for sscanf */
  long f,t;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: NListCommand(): me not in Context," \
                     " no listing of nodes\n",me))
    return(OKCODE);
  }
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"nlist","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  idopt = LV_ID;
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

#ifdef ModelP
    case 'g' :
      mode = DO_ID;
      idopt = LV_GID;
      res = sscanf(argv[i]," g %s",buff);
      fromN = toN = (DDD_GID) strtol(buff, NULL, 0);
      break;
#endif
    case 'k' :
      mode = DO_ID;
      idopt = LV_KEY;
      res = sscanf(argv[i]," k %s",buff);
      fromN = toN = (INT) strtol(buff, NULL, 0);
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
    ListNodeRange(theMG,fromN,toN,idopt,dataopt,boundaryopt,neighbouropt,verboseopt);
    break;

  case DO_ALL :
    ListNodeRange(theMG,0,MAX_I,idopt,dataopt,boundaryopt,neighbouropt,verboseopt);
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


/** \brief Implementation of \ref elist. */
static INT EListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,fromE,toE,res,mode,idopt,dataopt,boundaryopt,neighbouropt,verboseopt,levelopt;
  char buff[32];

  /* following variables: keep type for sscanf */
  long f,t;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: EListCommand(): me not in Context," \
                     " no listing of elements\n",me))
    return(OKCODE);
  }
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"elist","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  idopt = LV_ID;
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

#ifdef ModelP
    case 'g' :
      mode = DO_ID;
      idopt = LV_GID;
      res = sscanf(argv[i]," g %s",buff);
      fromE = toE = (DDD_GID) strtol(buff, NULL, 0);
      break;
#endif
    case 'k' :
      mode = DO_ID;
      idopt = LV_KEY;
      res = sscanf(argv[i]," k %s",buff);
      fromE = toE = (INT) strtol(buff, NULL, 0);
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
    ListElementRange(theMG,fromE,toE,idopt,dataopt,boundaryopt,neighbouropt,verboseopt,levelopt);
    break;

  case DO_ALL :
    ListElementRange(theMG,0,MAX_I,idopt,dataopt,boundaryopt,neighbouropt,verboseopt,levelopt);
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


/** \brief Implementation of \ref slist. */
static INT SelectionListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,dataopt,boundaryopt,neighbouropt,verboseopt;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: SelectionListCommand(): me not in Context,"\
                     " no listing of selection\n",me))
    return(OKCODE);
  }
        #endif

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

  case vectorSelection :
    UserWrite("sorry, this service is not available for vector selections\n");
    break;

  default :
    PrintErrorMessage('W',"slist","selectionmode ???");
    return (PARAMERRORCODE);
  }

  return(OKCODE);
}


/** \brief Implementation of \ref rlist. */
static INT RuleListCommand (INT argc, char **argv)
{
  INT i,allopt,rn,rv,tag;
  char etype[32];

  rn = -1;
  allopt = FALSE;

  /* check options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      allopt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("rlist",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  /* scan parameters */
  if (allopt == FALSE)
    rv = sscanf(argv[0],"rlist %31[triquatethexa] %d",etype,&rn);
  else
    rv = sscanf(argv[0],"rlist %31[triaquadtetrahexa]",etype);

  tag = -1;
        #ifdef __TWODIM__
  if (strcmp("tri",etype)==0) tag = TRIANGLE;
  if (strcmp("qua",etype)==0) tag = QUADRILATERAL;
        #endif
        #ifdef __THREEDIM__
  if (strcmp("tet",etype)==0) tag = TETRAHEDRON;
  if (strcmp("hex",etype)==0) tag = HEXAHEDRON;
        #endif

  if (tag==-1)
  {
    PrintErrorMessage('E',"rlist","wrong element type");
    return (CMDERRORCODE);
  }

  if (rn==-1 && allopt==FALSE || rn>=0 && allopt==TRUE)
  {
    PrintErrorMessage('E',"rlist","specify rulenumber OR $a option!");
    return (CMDERRORCODE);
  }

  if (allopt == TRUE)
    for (i=0; i<MaxRules[tag]; i++)
      ShowRefRule(tag,i);
  else
    ShowRefRule(tag,rn);

  return(OKCODE);
}


/** \todo Please doc me! */
static INT PrintValueCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *theVD;
  double val;
  char name[NAMESIZE];
  char value[VALUELEN];
  int found = FALSE;
  int idx;
  int n;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"printvalue","no open multigrid");
    return (CMDERRORCODE);
  }
  if (sscanf(argv[0],"printvalue %s %d",name,&n)!=2)
  {
    PrintErrorMessage('E',"printvalue","could not scan vec desc and selection number");
    return (PARAMERRORCODE);
  }
  theVD = GetVecDataDescByName(theMG,name);
  if (theVD==NULL)
  {
    PrintErrorMessageF('E',"printvalue","vec desc '%s' not found",name);
    return (PARAMERRORCODE);
  }

  if ((SELECTIONMODE(theMG)==vectorSelection) && (SELECTIONSIZE(theMG)>n))
  {
    VECTOR *vec = (VECTOR *)SELECTIONOBJECT(theMG,n);
    if (VD_ISDEF_IN_TYPE(theVD,VTYPE(vec)))
    {
      val = VVALUE(vec,VD_CMP_OF_TYPE(theVD,VTYPE(vec),0));
      idx = VINDEX(vec);
      found = TRUE;
    }
  }
  if (found)
    sprintf(buffer,"%.10e",val);
  else
    sprintf(buffer,"---");

  UserWriteF("value 0 of %s in vec %d = %s\n",name,idx,buffer);
  if (ReadArgvChar("s",value,argc,argv)==0)
    if (SetStringVar(value,buffer))
    {
      PrintErrorMessageF('E',"printvalue","coul not write onto string var '%s'",value);
      return (PARAMERRORCODE);
    }

  return OKCODE;
}


/** \brief Implementation of \ref vmlist. */
static INT VMListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  FORMAT *theFormat;
  GRID *theGrid;
  INT i,j,fl,tl,fromV,toV,res,mode,idopt,dataopt,matrixopt,vclass,vnclass,datatypes,modifiers;
  VECDATA_DESC *theVD;
  MATDATA_DESC *theMD;
  char value[VALUELEN];
  char buff[32];
  /* following variables: keep type for sscanf */
  long f,t;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: VMListCommand(): me not in Context," \
                     " no listing of VM\n",me))
    return(OKCODE);
  }
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"vmlist","no open multigrid");
    return (CMDERRORCODE);
  }
  theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));
  theFormat = MGFORMAT(theMG);

  if (ReadArgvINT("vclass",&vclass,argc,argv))
    vclass = 3;
  if (ReadArgvINT("vnclass",&vnclass,argc,argv))
    vnclass = 3;

  if (ReadArgvChar("vmlist",value,argc,argv) == 0) {
    theVD = GetVecDataDescByName(theMG,value);
    if (theVD != NULL) {
      if (ReadArgvOption("S",argc,argv))
        PrintSVector(theMG,theVD);
            #ifdef __INTERPOLATION_MATRIX__
      else if (ReadArgvOption("I",argc,argv))
        PrintIMatrix(theGrid,theVD,vclass,vnclass);
            #endif
      else if (ReadArgvOption("s",argc,argv))
      {
        /* get selection list */
        if ((SELECTIONMODE(theMG)==vectorSelection) && (SELECTIONSIZE(theMG)>=1))
        {
          VECTOR **vlist = (VECTOR**)malloc((SELECTIONSIZE(theMG)+1)*sizeof(VECTOR*));
          if (vlist!=NULL)
          {
            int i;
            for (i=0; i<SELECTIONSIZE(theMG); i++)
              vlist[i] = (VECTOR *)SELECTIONOBJECT(theMG,i);
            vlist[SELECTIONSIZE(theMG)] = NULL;
            PrintVectorListX((const VECTOR **)vlist,theVD,vclass,vnclass,UserWriteF);
            free(vlist);
          }
        }
      }
      else
        PrintVector(theGrid,theVD,vclass,vnclass);
      return(OKCODE);
    }
    theMD = GetMatDataDescByName(theMG,value);
    if (theMD != NULL) {
      if (ReadArgvOption("T",argc,argv))
        PrintTMatrix(theGrid,theMD,vclass,vnclass);
      else if (ReadArgvOption("D",argc,argv))
        PrintDiagMatrix(theGrid,theMD,vclass,vnclass);
      else
        PrintMatrix(theGrid,theMD,vclass,vnclass);
      return(OKCODE);
    }
  }
  modifiers = LV_MOD_DEFAULT;
  if (ReadArgvINT("skip",&j,argc,argv)==0)
    if (j)
      SET_FLAG(modifiers,LV_SKIP);
    else
      CLEAR_FLAG(modifiers,LV_SKIP);
  if (ReadArgvINT("pos",&j,argc,argv)==0)
    if (j)
      SET_FLAG(modifiers,LV_POS);
    else
      CLEAR_FLAG(modifiers,LV_POS);
  if (ReadArgvINT("obj",&j,argc,argv)==0)
    if (j)
      SET_FLAG(modifiers,LV_VO_INFO);
    else
      CLEAR_FLAG(modifiers,LV_VO_INFO);

  /* check options */
  datatypes = 0;
  idopt = LV_ID;
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

#ifdef ModelP
    case 'g' :
      mode = DO_ID;
      idopt = LV_GID;
      res = sscanf(argv[i]," g %s",buff);
      fromV = toV = (DDD_GID) strtol(buff, NULL, 0);
      break;
#endif
    case 'k' :
      mode = DO_ID;
      idopt = LV_KEY;
      res = sscanf(argv[i]," k %s",buff);
      fromV = toV = (INT) strtol(buff, NULL, 0);
      break;

    case 's' :
      if (strncmp(argv[i],"skip",4)==0)
        /* handled by ReadArgvINT */
        break;
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

    case 't' :
      for (j=0; j<NVECTYPES; j++)
        if (FMT_S_VEC_TP(theFormat,j)>0)
          if (strchr(argv[i]+1,FMT_VTYPE_NAME(theFormat,j))!=NULL)
            datatypes |= BITWISE_TYPE(j);
      break;

    case 'm' :
      matrixopt = TRUE;
      break;

    case 'z' :
      matrixopt = -TRUE;
      break;

    case 'p' :
    case 'o' :
      /* handled by ReadArgvINT */
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("vmlist",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }
  if (datatypes==0)
    /* default */
    for (j=0; j<NVECTYPES; j++)
      datatypes |= BITWISE_TYPE(j);

  switch (mode)
  {
  case DO_ID :
    ListVectorRange(theMG,fl,tl,fromV,toV,idopt,matrixopt,dataopt,datatypes,modifiers);
    break;

  case DO_ALL :
    ListVectorRange(theMG,fl,tl,0,MAX_I,idopt,matrixopt,dataopt,datatypes,modifiers);
    break;

  case DO_SELECTION :
    if (SELECTIONMODE(theMG)==elementSelection)
      ListVectorOfElementSelection(theMG,matrixopt,dataopt,modifiers);
    else
      ListVectorSelection(theMG,matrixopt,dataopt,modifiers);
    break;

  default :
    PrintErrorMessage('E',"vmlist","specify either the a, s or i option");
    return (PARAMERRORCODE);
  }

  return(OKCODE);
}


/** \todo Please doc me! */
static int ReadMatrixDimensions (char *name, int *n, int *na)
{
  int i;

  FILE *stream = fileopen(name,"r");
  if (stream == NULL) return(1);
  fscanf(stream," %d\n",n);
  for (i=0; i<=*n; i++)
    fscanf(stream," %d ",na);
  fclose(stream);

  return(0);
}

/** \todo Please doc me! */
static int ReadMatrix (char *name, int n, int *ia, int *ja, double *a)
{
  int i;

  FILE *stream = fileopen(name,"r");
  if (stream == NULL) return(1);
  fscanf(stream," %d\n",&i);
  if (i != n) return(1);
  for (i=0; i<=n; i++)
    fscanf(stream," %d ",ia+i);
  fscanf(stream,"\n");
  for (i=0; i<ia[n]; i++)
    fscanf(stream," %d ",ja+i);
  fscanf(stream,"\n");
  for (i=0; i<ia[n]; i++)
    fscanf(stream," %lf ",a+i);
  fscanf(stream,"\n");
  fclose(stream);

  return(0);
}

/** \todo Please doc me! */
static int WriteMatrix (char *name, int n, int *ia, int *ja, double *a)
{
  int i;

  FILE *stream = fileopen(name,"w");
  if (stream == NULL) return(1);
  fprintf(stream," %d\n",n);
  for (i=0; i<=n; i++)
    fprintf(stream," %d ",ia[i]);
  fprintf(stream,"\n");
  for (i=0; i<ia[n]; i++)
    fprintf(stream," %d ",ja[i]);
  fprintf(stream,"\n");
  for (i=0; i<ia[n]; i++)
    fprintf(stream," %f ",a[i]);
  fprintf(stream,"\n");
  fclose(stream);

  return(0);
}

/** \todo Please doc me! */
static int WriteMatrixfmt (char *name, int n, int *ia, int *ja, double *a,
                           int inc)
{
  int i;

  FILE *stream = fileopen(name,"w");
  if (stream == NULL) return(1);
  fprintf(stream,"%d %d",n,ia[n]+inc);
  for (i=0; i<=n; i++) {
    if ((i%10) == 0) fprintf(stream,"\n");
    fprintf(stream,"%6d",ia[i]+inc);
  }
  for (i=0; i<ia[n]; i++) {
    if ((i%3) == 0) fprintf(stream,"\n");
    fprintf(stream,"%6d %18.9f",ja[i]+inc,a[i]);
  }
  fprintf(stream,"\n");
  fclose(stream);

  return(0);
}

/** \brief Implementation of \ref convert. */
static INT ConvertCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  HEAP *theHeap;
  MATDATA_DESC *A;
  int n,nn,*ia,*ja,inc;
  double *a,*r;
  char name[32];
  INT i,j,MarkKey,symmetric,ncomp;

  theMG = GetCurrentMultigrid();
  if (theMG==NULL) {
    PrintErrorMessage('E',"convert","no current multigrid");
    return(CMDERRORCODE);
  }
  theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));
  if ((A = ReadArgvMatDesc(theMG,"convert",argc,argv))==NULL) {
    PrintErrorMessage('E',"convert","could not read symbol");
    return (PARAMERRORCODE);
  }
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap,&MarkKey);
  symmetric = ReadArgvOption("symmetric",argc,argv);
  inc = ReadArgvOption("inc",argc,argv);
  if (ReadArgvINT("ncomp",&ncomp,argc,argv))
    ncomp = 1;
  if (ReadArgvChar("r",name,argc,argv) == 0) {
    if (ReadMatrixDimensions(name,&n,&nn)) {
      PrintErrorMessage('E',"convert",
                        "could not read matrix dimensions");
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
      return(CMDERRORCODE);
    }
    ia = (int *)GetTmpMem(theHeap,sizeof(int) * (n+1),MarkKey);
    a = (double *)GetTmpMem(theHeap,sizeof(double) * nn,MarkKey);
    ja = (int *)GetTmpMem(theHeap,sizeof(int) * nn,MarkKey);
    if ((ia == NULL) || (a == NULL) || (ja == NULL)) {
      PrintErrorMessage('E',"convert",
                        "could not allocate memory");
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
      return(CMDERRORCODE);
    }
    if (ReadMatrix(name,n,ia,ja,a)) {
      PrintErrorMessage('E',"convert","could write matrix");
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
      return(CMDERRORCODE);
    }
  }
  else if (ConvertMatrix(theGrid,MGHEAP(theMG),MarkKey,A,symmetric,
                         &n,&ia,&ja,&a)) {
    PrintErrorMessage('E',"convert","could not read matrix");
    ReleaseTmpMem(MGHEAP(theMG),MarkKey);
    return(CMDERRORCODE);
  }
  if (ReadArgvChar("f",name,argc,argv) == 0)
    if (ReadArgvOption("fmt",argc,argv)) {
      if (WriteMatrixfmt(name,n,ia,ja,a,inc)) {
        PrintErrorMessage('E',"convert","could write matrix");
        ReleaseTmpMem(MGHEAP(theMG),MarkKey);
        return(CMDERRORCODE);
      }
    }
    else if (WriteMatrix(name,n,ia,ja,a)) {
      PrintErrorMessage('E',"convert","could write matrix");
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
      return(CMDERRORCODE);
    }
  if (ReadArgvOption("p",argc,argv)) {
    r = (DOUBLE *)GetTmpMem(MGHEAP(theMG),sizeof(DOUBLE) * n,MarkKey);
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++)
        r[j] = 0.0;
      for (j=ia[i]; j<ia[i+1]; j++)
        r[ja[j]] = a[j];
      for (j=0; j<n; j++)
        UserWriteF("%8.4f",r[j]);
      UserWrite("\n");
    }
  }
  ReleaseTmpMem(MGHEAP(theMG),MarkKey);
  return (OKCODE);
}



/** \brief Implementation of \ref in. */
static INT InsertInnerNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  DOUBLE xc[DIM];

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"in","no open multigrid");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"in %lf %lf %lf",xc,xc+1,xc+2)!=DIM)
  {
    PrintErrorMessageF('E',"in","specify %d coordinates for an inner node",(int)DIM);
    return (PARAMERRORCODE);
  }

  /* NB: toplevel=0 is checked by InsertInnerNode() */
  if (InsertInnerNode(GRID_ON_LEVEL(theMG,0),xc)==NULL)
  {
    PrintErrorMessage('E',"in","inserting an inner node failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}

/** \todo Please doc me! */
static INT NGInsertInnerNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  DOUBLE xc[DIM];
  static int n;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"in","no open multigrid");
    return (CMDERRORCODE);
  }

  UserWriteF("# IPoint %d\n",n);
  n++;
  UserWriteF("# %s\n",argv[0]);
  if (sscanf(argv[0],"ngin %lf %lf %lf",xc,xc+1,xc+2)!=DIM)
  {
    PrintErrorMessageF('E',"in","specify %d coordinates for an inner node",(int)DIM);
    return (PARAMERRORCODE);
  }

  UserWriteF("I %lf %lf %lf;\n",xc[0],xc[1],xc[2]);

  return (OKCODE);
}


/** \brief Implementation of \ref bn. */
static INT InsertBoundaryNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  BNDP *bndp;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"bn","no open multigrid");
    return (CMDERRORCODE);
  }

  bndp = BVP_InsertBndP (MGHEAP(theMG),MG_BVP(theMG),argc,argv);
  if (bndp == NULL)
  {
    PrintErrorMessage('E',"bn","inserting a boundary point failed");
    return (CMDERRORCODE);
  }

  if (InsertBoundaryNode(GRID_ON_LEVEL(theMG,0),bndp)==NULL)
  {
    PrintErrorMessage('E',"bn","inserting a boundary node failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}


/** \todo Please doc me! */
static INT NGInsertBoundaryNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  BNDP *bndp;
  static int i;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"ngbn","no open multigrid");
    return (CMDERRORCODE);
  }

  UserWriteF("# BPoint %d \n",i);
  /* this works only for LGM domain and does no real insertion !!! */
  bndp = BVP_InsertBndP (MGHEAP(theMG),MG_BVP(theMG),argc,argv);
  if (bndp == NULL)
  {
    i++;
    return (OKCODE);
  }
  return (CMDERRORCODE);
}


/** \brief Implementation of \ref gn. */
static INT InsertGlobalNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  BNDP *bndp;
  int i;
  int ropt = FALSE;
  int err = OKCODE;
  DOUBLE resolution;
  INT n = 2;
  INT my_argc = 0;
  char **my_argv;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"gn","no open multigrid");
    return (CMDERRORCODE);
  }

  /* assemble command line for bn */
  if (ReadArgvDOUBLE("r",&resolution,argc,argv)==0)
  {
    ropt = TRUE;
    n++;
  }
  my_argv = (char**)malloc(n*sizeof(char*));
  if (my_argv==NULL)
    return CMDERRORCODE;
  my_argv[my_argc] = StrDup(argv[0]);
  if (my_argv[my_argc]==NULL)
  {
    err = CMDERRORCODE;
    goto Exit_gn;
  }
  my_argv[my_argc][0] = 'b';                    /* gn --> bn */
  my_argc++;

  my_argv[my_argc] = StrDup("g");
  if (my_argv[my_argc]==NULL)
  {
    err = CMDERRORCODE;
    goto Exit_gn;
  }
  my_argc++;

  if (ropt)
  {
    char r_opt[64];
    sprintf(r_opt,"$r %g",resolution);
    my_argv[my_argc] = StrDup(r_opt);
    if (my_argv[my_argc]==NULL)
    {
      err = CMDERRORCODE;
      goto Exit_gn;
    }
    my_argc++;
  }
  ASSERT(n==my_argc);

  /* try inserting a boundary node */
  bndp = BVP_InsertBndP (MGHEAP(theMG),MG_BVP(theMG),my_argc,my_argv);
  if (bndp == NULL)
  {
    double x[DIM_MAX];
    DOUBLE xc[DIM];

    /* failed: try inserting an inner node */
    if (sscanf(argv[0],"gn %lf %lf %lf",x,x+1,x+2)!=DIM)
    {
      PrintErrorMessageF('E',"gn","specify %d global coordinates",(int)DIM);
      err = PARAMERRORCODE;
      goto Exit_gn;
    }
    for (i=0; i<DIM; i++)
      xc[i] = x[i];

    /* NB: toplevel=0 is checked by InsertInnerNode() */
    if (InsertInnerNode(GRID_ON_LEVEL(theMG,0),xc)==NULL)
    {
      PrintErrorMessage('E',"gn","inserting an inner node failed");
      err = CMDERRORCODE;
      goto Exit_gn;
    }
    UserWrite("  ### gn: inserted a in\n");
  }
  else if (InsertBoundaryNode(GRID_ON_LEVEL(theMG,0),bndp)==NULL)
  {
    PrintErrorMessage('E',"gn","inserting a boundary node failed");
    err =  CMDERRORCODE;
    goto Exit_gn;
  }
  else
    UserWrite("  ### gn: inserted a bn\n");

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

Exit_gn:
  for (i=0; i<my_argc; i++)
    if (my_argv[i]!=NULL)
      free(my_argv[i]);
  free(my_argv);

  return err;
}


/** \brief Implementation of \ref deln. */
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
        if (DeleteNode(GRID_ON_LEVEL(theMG,0),(NODE *)SELECTIONOBJECT(theMG,i))!=GM_OK)
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
  if (DeleteNodeWithID(GRID_ON_LEVEL(theMG,0),id)!=GM_OK)
  {
    PrintErrorMessage('E',"deln","deleting the node failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}


/** \brief Implementation of \ref move. */
static INT MoveNodeCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VERTEX *myVertex;
  NODE *theNode;
  DOUBLE xc[DIM];
  INT type,i,j,level,relative;

  /* following variables: keep type for sscanf */
  double x[DIM_MAX];
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
      PrintErrorMessageF('E',"move","node with ID %ld not found",(long)id);
      return (CMDERRORCODE);
    }
  }

  /* check options */
  relative = FALSE;
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
        PrintErrorMessageF('E',"move","node with ID %ld is no inner node",(long)id);
        return (CMDERRORCODE);
      }
      type = IVOBJ;
      if (sscanf(argv[i],"i %lf %lf %lf",x,x+1,x+2)!=DIM)
      {
        PrintErrorMessageF('E',"move","specify %d new coordinates for an inner node",(int)DIM);
        return (PARAMERRORCODE);
      }
      for (j=0; j<DIM; j++)
        xc[j] = x[j];
      break;

    case 'b' :
      if (OBJT(MYVERTEX(theNode))!=BVOBJ)
      {
        PrintErrorMessageF('E',"move","node with ID %ld is no boundary node",(long)id);
        return (CMDERRORCODE);
      }
      type = BVOBJ;
      if (sscanf(argv[i],"b %d %lf %lf",&segid,x,x+1)!=1+DIM_OF_BND)
      {
        PrintErrorMessageF('E',"move","specify the segment if and %d new coordinates for a boundary node",(int)DIM_OF_BND);
        return (PARAMERRORCODE);
      }
      for (j=0; j<DIM_OF_BND; j++)
        xc[j] = x[j];
      break;

    case 'r' :
      relative = TRUE;
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

  myVertex = MYVERTEX(theNode);
  if (type==IVOBJ)
  {
    if (relative)
    {
      /* move relative to old position */
      for (j=0; j<DIM; j++)
        xc[j] += CVECT(myVertex)[j];
    }
    if (MoveNode(theMG,theNode,xc,TRUE)!=GM_OK)
    {
      PrintErrorMessage('E',"move","failed moving the node");
      return (CMDERRORCODE);
    }
  }
  else if (type==BVOBJ)
  {
    PrintErrorMessage('E',"move","moving boundary nodes not implemented yet");
    return (CMDERRORCODE);
  }
  else
  {
    PrintHelp("move",HELPITEM," (either i or b option is mandatory)");
    return (PARAMERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);

  return (OKCODE);
}


/** \brief Implementation of \ref ie. */
static INT InsertElementCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  NODE *theNode,*theNodes[MAX_CORNERS_OF_ELEM];
  char *token,*vstr;
  INT i,nNodes,Id[MAX_CORNERS_OF_ELEM];

  /* following variables: keep type for sscanf */
  int id;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

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
        for (nNodes=0; nNodes<SELECTIONSIZE(theMG); nNodes++)
        {
          theNode = (NODE *)SELECTIONOBJECT(theMG,nNodes);
          if (nNodes>=MAX_CORNERS_OF_ELEM)
          {
            PrintErrorMessage('E',"ie","too many nodes are in the selection");
            return (CMDERRORCODE);
          }
          theNodes[nNodes] = theNode;
        }
      else
      {
        PrintErrorMessage('E',"ie","objects other than nodes are in the selection");
        return (PARAMERRORCODE);
      }
      if (nNodes==0)
      {
        PrintErrorMessage('E',"ie","no nodes are in the selection");
        return (PARAMERRORCODE);
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
    if (InsertElement(GRID_ON_LEVEL(theMG,0),nNodes,theNodes,NULL,NULL,NULL)==NULL)
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
      PrintErrorMessageF('E',"ie","specify at most %d id's",(int)MAX_CORNERS_OF_ELEM);
      return (PARAMERRORCODE);                                  /* too many items */
    }
    if (sscanf(token," %d",&id)!=1)
    {
      PrintErrorMessageF('E',"ie","could not read the id of corner no %d",(int)i);
      return (PARAMERRORCODE);                                  /* too many items */
    }
    Id[nNodes++] = id;
    token = strtok(NULL,WHITESPACE);
  }

  /* NB: toplevel=0 is checked by InsertElementFromIDs() */
  if (InsertElementFromIDs(GRID_ON_LEVEL(theMG,0),nNodes,Id,NULL)==NULL)
  {
    PrintErrorMessage('E',"ie","inserting the element failed");
    return (CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}

/** \todo Please doc me! */
static INT NGInsertElementCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char *token,*vstr;
  INT i,bf,nNodes,Id[MAX_CORNERS_OF_ELEM];
  static int n;

  /* following variables: keep type for sscanf */
  int id;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"ngie","no open multigrid");
    return (CMDERRORCODE);
  }

  /* set vstr after 'ngie' */
  if ((vstr=strchr(argv[0],'e'))!=NULL)
    ++vstr;
  else
    return (CMDERRORCODE);

  UserWriteF("# %s\n",argv[0]);
  UserWriteF("# element %d\n",n);
  n++;
  UserWriteF("E ");

  /* we need the string split into tokens */
  nNodes = bf = 0;
  token = strtok(vstr,WHITESPACE);
  i = 0;
  while (token!=NULL)
  {
    /* separator for boundary faces */
    if (strcmp(token,"F") == 0)
    {
      UserWriteF("\n");
      bf = 1;
      goto NEXTTOKEN;
    }

    /* read next boundary face */
    if (bf > 0)
    {

      if (sscanf(token," %d",&id)!=1)
      {
        PrintErrorMessageF('E',"ngie","could not read the id of boundary face no %d",(int)bf);
        return (PARAMERRORCODE);                                        /* too many items */
      }
      UserWriteF("	F");
      switch (nNodes)
      {
      case 4 :
      case 5 :
      case 6 :
        UserWriteF("ngie: elementtype = %d not implemented!\n",nNodes);
        break;
      case 8 :
      {
        INT n;
        /* assert(id<SIDES_OF_ELEM_TAG(nNodes));*/

        for (n=0; n<CORNERS_OF_SIDE_TAG(7,id); n++)
        {
          UserWriteF(" %d",Id[CORNER_OF_SIDE_TAG(7,id,n)]);
        }
        UserWriteF("\n");
      }
      break;
      default :
        assert(0);
      }

      bf++;
      goto NEXTTOKEN;
    }

    /* read next node */
    if (nNodes>=MAX_CORNERS_OF_ELEM)
    {
      PrintErrorMessageF('E',"ngie","specify at most %d id's",(int)MAX_CORNERS_OF_ELEM);
      return (PARAMERRORCODE);                                  /* too many items */
    }
    if (sscanf(token," %d",&id)!=1)
    {
      PrintErrorMessageF('E',"ngie","could not read the id of corner no %d",(int)nNodes);
      return (PARAMERRORCODE);                                  /* too many items */
    }

    /* first id is subdomain */
    if (i > 0)
    {
      Id[nNodes] = id;
      nNodes++;
    }
    UserWriteF(" %d",id);

NEXTTOKEN:
    token = strtok(NULL,WHITESPACE);
    i++;
  }

  UserWriteF(";\n");

  return (OKCODE);
}


/** \brief Implementation of \ref dele. */
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


/** \brief Implementation of \ref adapt. */
static INT AdaptCommand (INT argc, char **argv)
{
  MULTIGRID       *theMG;
  INT i,mode,mark,rv;
  INT seq;
  INT mgtest;
  EVECTOR         *theElemEvalDirection;

        #ifdef ModelP
  if (!CONTEXT(me))
  {
    PRINTDEBUG(ui,0,("%2d: AdaptCommand(): me not in Context,"
                     " grid not refined\n",me))
    return (OKCODE);
  }
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"adapt","no open multigrid");
    return (CMDERRORCODE);
  }

  /* init defaults */
  seq = GM_REFINE_PARALLEL;
  mgtest = GM_REFINE_NOHEAPTEST;

  /* check options */
  theElemEvalDirection = NULL;
  mode = GM_REFINE_TRULY_LOCAL;
  mark = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      mark = MARK_ALL;
      break;

            #ifdef __THREEDIM__
    case 'd' :
      if (sscanf(argv[i],"a %s",buffer)==1)
        theElemEvalDirection = GetElementVectorEvalProc(buffer);
      if (theElemEvalDirection==NULL)
        UserWrite("direction eval fct not found: taking shortest interior edge\n");
      break;
                        #endif

    case 'g' :
      mode = mode | GM_COPY_ALL;
      break;

    case 'h' :
      mode = mode | GM_REFINE_NOT_CLOSED;
      break;

    case 's' :
      seq = GM_REFINE_SEQUENTIAL;
      break;

    case 't' :
      mgtest = GM_REFINE_HEAPTEST;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("refine",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

        #ifdef ModelP
  /* currently only this is supported in parallel */
  if (0)
    if (procs > 1)
    {
      mark = MARK_ALL;
      mode = mode | GM_REFINE_NOT_CLOSED;
    }
        #endif

  if (mark == MARK_ALL)
  {
    INT l,nmarked;
    ELEMENT *theElement;

    nmarked = 0;

    for (l=TOPLEVEL(theMG); l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l)); theElement!=NULL; theElement=SUCCE(theElement))
      {
        if (EstimateHere(theElement))
          if ((rv = MarkForRefinement(theElement,RED,0))!=0)
          {
            l = TOPLEVEL(theMG);
            break;
          }
          else
            nmarked++;
      }
    UserWriteF("%d: %d elements marked for regular refinement\n",me,nmarked);
  }

  /* get velocity */
        #ifdef __THREEDIM__
  SetAlignmentPtr (theMG, theElemEvalDirection);
        #endif

  rv = AdaptMultiGrid(theMG,mode,seq,mgtest);

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  switch (rv)
  {
  case GM_OK :
    UserWriteF(" %s refined\n",ENVITEM_NAME(theMG));
    SetStringVar(":errno","0");
    return (OKCODE);

  case GM_COARSE_NOT_FIXED :
    PrintErrorMessage('E',"refine","do 'fixcoarsegrid' first and then refine!");
    SetStringVar(":errno","1");
    return (CMDERRORCODE);

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

/** \brief Implementation of \ref fixcoarsegrid. */
static INT FixCoarseGridCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

  theMG = currMG;

  PRINTDEBUG(ui,2,("%d: FixCoarseGrid currMG %x fixed %d\n",
                   me,theMG,MG_COARSE_FIXED(theMG)));

  if (theMG==NULL) {
    PrintErrorMessage('E',"fixcoarsegrid","no open multigrid");
    return (CMDERRORCODE);
  }
  if (FixCoarseGrid(theMG)) return (CMDERRORCODE);

  PRINTDEBUG(ui,2,("%d: FixCoarseGrid currMG %x fixed %d\n",
                   me,theMG,MG_COARSE_FIXED(theMG)));

  return(OKCODE);
}


/** \brief Implementation of \ref collapse. */
static INT CollapseCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

  theMG = currMG;

  if (theMG==NULL) {
    PrintErrorMessage('E',"collapse","no open multigrid");
    return (CMDERRORCODE);
  }
  if (Collapse(theMG)) return (CMDERRORCODE);

  return(OKCODE);
}


/** \brief Implementation of \ref mark. */
static INT MarkCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  ELEMENT *theElement;
  char rulename[32];
  INT i,j,l,mode,rv,sid;
  enum RefinementRule Rule;
  DOUBLE_VECTOR global;
  DOUBLE x,y;
  long nmarked;
#       ifdef __THREEDIM__
  DOUBLE X,Y;
  DOUBLE z,Z;
#       endif

  /* following variables: keep type for sscanf */
  int id,idfrom,idto;
  INT Side;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%d: MarkCommand() me not in Context," \
                     " nothing marked\n",me))
    return (OKCODE);
  }
        #endif

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"mark","no open multigrid");
    return (CMDERRORCODE);
  }

  /* first check help option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='h')
    {
      UserWrite("the following rules are available:\n");
      for (j=0; j<NO_OF_RULES && myMR[j].RuleName!=NULL; j++)
      {
        UserWrite(myMR[j].RuleName);
        UserWrite("\n");
      }
      return (OKCODE);
    }

  /* scan parameters */
  /*rv = sscanf(argv[0],"mark %31[redbluecopycoarsnoi_123qtpritet2hex] %d",rulename,&Side);*/
  rv = sscanf(argv[0],"mark %31[a-z_0-9] %d",rulename,&Side);
  if (rv<1)
  {
    /* set the default rule */
    strcpy(rulename,myMR[0].RuleName);
    Rule = (enum RefinementRule)myMR[0].RuleId;
    Side = NO_SIDE_SPECIFIED;
  }
  else
  {
    Rule = (enum RefinementRule)NO_RULE_SPECIFIED;
    for (i=0; i<NO_OF_RULES; i++)
      if (strcmp(rulename,myMR[i].RuleName)==0)
      {
        Rule = (enum RefinementRule)myMR[i].RuleId;
        break;
      }

    if (Rule==NO_RULE_SPECIFIED)
    {
      PrintErrorMessageF('E',"mark","unknown rule '%s'",rulename);
      return (PARAMERRORCODE);
    }

    if (rv!=2)
      Side = NO_SIDE_SPECIFIED;
  }

  if (ReadArgvOption("c",argc, argv))
  {
    for (i=0; i<=TOPLEVEL(theMG); i++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,i));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          MarkForRefinement(theElement,NO_REFINEMENT,0);

    UserWrite("all refinement marks removed\n");

    return(OKCODE);
  }

  if (ReadArgvDOUBLE("x",&x,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if (XC(MYVERTEX(CORNER(theElement,j))) < x)
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in x < %f marked for refinement\n",
               (float) x);

    return(OKCODE);
  }

  if (ReadArgvDOUBLE("X",&x,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if (XC(MYVERTEX(CORNER(theElement,j))) > x)
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in x > %f marked for refinement\n",
               (float) x);

    return(OKCODE);
  }

  if (ReadArgvDOUBLE("y",&y,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if (YC(MYVERTEX(CORNER(theElement,j))) < y)
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in y < %f marked for refinement\n",
               (float) y);

    return(OKCODE);
  }

  if (ReadArgvDOUBLE("Y",&y,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if (YC(MYVERTEX(CORNER(theElement,j))) > y)
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in y > %f marked for refinement\n",
               (float) y);

    return(OKCODE);
  }

  if (ReadArgvDOUBLE("stripes",&x,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
        {
          INT flag=1;
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
          {
            DOUBLE xc=YC(MYVERTEX(CORNER(theElement, j)));
            xc=fmod(xc,(4*x));
            if ((xc<0.9*x)||(xc>2.1*x)) flag=0;
          }
          if(flag)
            MarkForRefinement(theElement,Rule,0);
        }
    UserWriteF("stripes %f\n", (float) x);

    return(OKCODE);
  }

  if (ReadArgvINT("S",&sid,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if (SUBDOMAIN(theElement) == sid)
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in subdomain %d marked for refinement\n",sid);

    return(OKCODE);
  }

#ifdef __THREEDIM__
  if ((ReadArgvDOUBLE("x0",&x,argc, argv)==0) &&
      (ReadArgvDOUBLE("x1",&X,argc, argv)==0) &&
      (ReadArgvDOUBLE("y0",&y,argc, argv)==0) &&
      (ReadArgvDOUBLE("y1",&Y,argc, argv)==0) &&
      (ReadArgvDOUBLE("z0",&z,argc, argv)==0) &&
      (ReadArgvDOUBLE("z1",&Z,argc, argv)==0))
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if ((XC(MYVERTEX(CORNER(theElement,j))) < X) &&
                (XC(MYVERTEX(CORNER(theElement,j))) > x) &&
                (YC(MYVERTEX(CORNER(theElement,j))) < Y) &&
                (YC(MYVERTEX(CORNER(theElement,j))) > y) &&
                (ZC(MYVERTEX(CORNER(theElement,j))) < Z) &&
                (ZC(MYVERTEX(CORNER(theElement,j))) > z))
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in box marked for refinement\n");

    return(OKCODE);
  }

  if (ReadArgvDOUBLE("z",&z,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if (ZC(MYVERTEX(CORNER(theElement,j))) < z)
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in z < %f marked for refinement\n",
               (float) z);

    return(OKCODE);
  }

  if (ReadArgvDOUBLE("Z",&z,argc, argv)==0)
  {
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement))
        if (EstimateHere(theElement))
          for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            if (ZC(MYVERTEX(CORNER(theElement,j))) > z)
              MarkForRefinement(theElement,Rule,0);

    UserWriteF("all elements in z > %f marked for refinement\n",
               (float) z);

    return(OKCODE);
  }

#endif

  if (ReadArgvPosition("pos",argc,argv,global)==0)
  {
    DOUBLE r;

    if (ReadArgvDOUBLE("r",&r,argc, argv)==0) {
      for (l=0; l<=TOPLEVEL(theMG); l++)
        for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
             theElement!=NULL; theElement=SUCCE(theElement))
          if (EstimateHere(theElement))
            for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
            {
              DOUBLE dist;

              V_DIM_EUKLIDNORM_OF_DIFF(global,
                                       CVECT(MYVERTEX(CORNER(theElement,j))),dist);
              if (dist <= r) {
                MarkForRefinement(theElement,Rule,0);
                break;
              }
            }
      UserWriteF("all elements in |x - p|  < %f marked for refinement\n",
                 (float) r);

      return(OKCODE);
    }
    theElement = FindElementOnSurface(theMG,global);
    if (theElement != NULL) {
      MarkForRefinement(theElement,Rule,0);
        #ifdef ModelP
      j = (INT) UG_GlobalSumDOUBLE(1.0);
      i = DDD_InfoGlobalId(PARHDRE(theElement));
    }
    else {
      j = (INT) UG_GlobalSumDOUBLE(0.0);
      i = -1;
    }
    if (j == 0)
      return(PARAMERRORCODE);

    for (l=0; l<j; l++) {
      rv = UG_GlobalMaxINT(i);
      UserWriteF("element GID %08x marked for refinement\n",rv);
      if (rv == i) i = -1;
    }
        #else
      UserWriteF("element %d marked for refinement\n",ID(theElement));
    }
    else
      return(PARAMERRORCODE);
            #endif

    return (OKCODE);
  }

  /* check options */
  mode = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"mark","specify only one option of a, b, i, s");
        return (PARAMERRORCODE);
      }
      mode = MARK_ALL;
      break;

    case 'i' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"mark","specify only one option of a, b, i, s");
        return (PARAMERRORCODE);
      }
      mode = MARK_ID;
      j = sscanf(argv[i],"i %d %d", &idfrom, &idto);
      if (j<1 || j>2)
      {
        PrintErrorMessage('E',"mark","cannot scan id(s)");
        return (PARAMERRORCODE);
      }
      if (j == 1) idto = idfrom;
      break;

    case 's' :
      if (mode!=FALSE)
      {
        PrintErrorMessage('E',"mark","specify only one option of a, b, i, s");
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
    PrintErrorMessage('E',"mark","specify exactly one option of a, b, i, s");
    return (PARAMERRORCODE);
  }

  /* print rule and side */
  if (Side==NO_SIDE_SPECIFIED)
    UserWriteF("   using rule %s (no side given)\n",rulename);
  else
    UserWriteF("   using rule %s, side %d\n",rulename,(int)Side);

  nmarked = 0;
  switch (mode)
  {
  case MARK_ALL :
    for (l=0; l<=TOPLEVEL(theMG); l++)
      for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,l));
           theElement!=NULL; theElement=SUCCE(theElement)) {
        if (EstimateHere(theElement))
          if ((rv = MarkForRefinement(theElement,
                                      Rule,Side))!=0)
          {
            l = TOPLEVEL(theMG);
            break;
          }
          else
            nmarked++;
      }
    break;

  case MARK_ID :
    theElement = NULL;
    for (id=idfrom; id<=idto; id++) {
      for (l=0; l<=TOPLEVEL(theMG); l++)
        if ((theElement=FindElementFromId(GRID_ON_LEVEL(theMG,l),id))!=NULL)
          break;

      if (theElement==NULL)
      {
        PrintErrorMessageF('W',"mark","element with ID %ld could not be found, nothing marked",id);
      }

      if (EstimateHere(theElement))
        if ((rv = MarkForRefinement(theElement,
                                    Rule,Side))!=0)
          break;
        else
          nmarked++;
    }
    break;

  case MARK_SELECTION :
    if (SELECTIONMODE(theMG)==elementSelection)
      for (i=0; i<SELECTIONSIZE(theMG); i++)
      {
        theElement = (ELEMENT *)SELECTIONOBJECT(theMG,i);
        if (EstimateHere(theElement))
        {
          if ((rv = MarkForRefinement(theElement,
                                      Rule,Side))!=0)
            break;
          else
            nmarked++;
        }
      }
    break;
  }

        #ifdef ModelP
  nmarked = UG_GlobalSumINT(nmarked);
        #endif
  UserWriteF(" %ld elements marked for refinement\n",nmarked);

  if (rv && theElement)
  {
    PrintErrorMessageF('W',"mark","rule could not be applied for element with ID %ld, nothing marked",ID(theElement));
    return (CMDERRORCODE);
  }
  else
    return(OKCODE);
}


/** \brief Implementation of \ref smooth. */
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
      smbdry = 1;
      break;

    case 'n' :
      smbdry = 2;
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


/** \brief Implementation of \ref smoothgrid. */
static INT SmoothGridCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  DOUBLE LimitLocDis;
  INT i,GridReset,lowLevel;
  INT bnd_num,bnd[22],fl,option,dummy;
  float fValue;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"smoothgrid","no current multigrid");
    return(CMDERRORCODE);
  }
  if (CURRENTLEVEL(theMG)==0)
  {
    PrintErrorMessage('E',"smoothgrid","cannot smooth grid on level 0");
    return(CMDERRORCODE);
  }
#ifdef __THREEDIM__
  PrintErrorMessage('E',"smoothgrid","3D not implemented yet");
  return(CMDERRORCODE);
#endif

  /* options */
  lowLevel = CURRENTLEVEL(theMG);
  option = 0;
  GridReset = FALSE;
  LimitLocDis = 0.3;
  bnd_num = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'r' :
      if (strcmp(argv[i],"reset")==0)
        GridReset=TRUE;
      else
      {
        PrintErrorMessageF('E',"smoothgrid","(invalid option '%s')",argv[i]);
        return (PARAMERRORCODE);
      }
      break;
    case 'l' :
      if (sscanf(argv[i],"limit %f",&fValue)!=1)
      {
        PrintErrorMessageF('E',"smoothgrid","(invalid option '%s')",argv[i]);
        return (PARAMERRORCODE);
      }
      LimitLocDis = (DOUBLE) fValue;
      if (LimitLocDis>=0.5 || LimitLocDis <=0)
      {
        PrintErrorMessage('E',"smoothgrid","specify a local limit between 0 and 0.5 (default 0.3)");
        return (PARAMERRORCODE);
      }
      break;

    case 'f' :
      if (sscanf(argv[i],"f %d",&fl)!=1)
      {
        PrintErrorMessageF('E',"smoothgrid","(invalid option '%s')",argv[i]);
        return (PARAMERRORCODE);
      }
      lowLevel = fl;
      break;
    case 'o' :
      if (strstr(argv[i],"ortho0")!=NULL)
      {
        if ((bnd_num=sscanf(argv[i],"ortho0 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                            &bnd[0],&bnd[1],&bnd[2],&bnd[3],&bnd[4],&bnd[5],&bnd[6],&bnd[7],&bnd[8],&bnd[9],
                            &bnd[10],&bnd[11],&bnd[12],&bnd[13],&bnd[14],&bnd[15],&bnd[16],&bnd[17],&bnd[18],
                            &bnd[19],&bnd[20],&bnd[21],&dummy))<1)
        {
          PrintErrorMessage('E',"smoothgrid","specify at least one boundary-id with 'ortho0' option");
          return (PARAMERRORCODE);
        }
        if (option!=0)
        {
          PrintErrorMessage('E',"smoothgrid","specify either $b, $ortho0 or $ortho1 option");
          return (PARAMERRORCODE);
        }
        option = 1;
      }
      else if (strstr(argv[i],"ortho1")!=NULL)
      {
        if ((bnd_num=sscanf(argv[i],"ortho1 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                            &bnd[0],&bnd[1],&bnd[2],&bnd[3],&bnd[4],&bnd[5],&bnd[6],&bnd[7],&bnd[8],&bnd[9],
                            &bnd[10],&bnd[11],&bnd[12],&bnd[13],&bnd[14],&bnd[15],&bnd[16],&bnd[17],&bnd[18],
                            &bnd[19],&bnd[20],&bnd[21],&dummy))<1)
        {
          PrintErrorMessage('E',"smoothgrid","specify at least one boundary-id with 'ortho1' option");
          return (PARAMERRORCODE);
        }
        if (option!=0)
        {
          PrintErrorMessage('E',"smoothgrid","specify either $b, $ortho0 or $ortho1 option");
          return (PARAMERRORCODE);
        }
        option = 2;
      }
      else
      {
        PrintErrorMessageF('E',"smoothgrid","(invalid option '%s')",argv[i]);
        return (PARAMERRORCODE);
      }
      if (bnd_num>21)
      {
        PrintErrorMessage('E',"smoothgrid","cannot process more than 9 boundaries with 'ortho' option");
        return (PARAMERRORCODE);
      }
      break;

    case 'b' :
      if (option!=0)
      {
        PrintErrorMessage('E',"smoothgrid","specify either $b, $ortho0 or $ortho1 option");
        return (PARAMERRORCODE);
      }
      option = 3;
      break;

    case 's' :
      break;
    default :
      PrintErrorMessageF('E',"smoothgrid","(invalid option '%s')",argv[i]);
      return (PARAMERRORCODE);
    }

  if (ReadArgvOption("spline",argc,argv))
  {
    if (option==0)
      option = 5;
    else if (option==1)
      option = 6;
    else if (option==2)
      option = 7;
  }
  if (ReadArgvOption("spline0",argc,argv))
  {
    option = 4;
  }
  UserWriteF("option = %d\n",option);
  if (GridReset==TRUE)
  {
    if (SmoothGridReset(theMG,lowLevel,CURRENTLEVEL(theMG))!=0) return(CMDERRORCODE);
  }
  else
  {
    lowLevel = MIN(FULLREFINELEVEL(theMG),lowLevel);
    lowLevel = MAX(lowLevel,1);
    if (SmoothGrid(theMG,lowLevel,CURRENTLEVEL(theMG),LimitLocDis,bnd_num,bnd,option)!=0) return(CMDERRORCODE);
  }

  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

  return (OKCODE);
}


/** \brief Implementation of \ref ordernodes. */
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
  if (RenumberMultiGrid(theMG,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0)!=GM_OK)
  {
    PrintErrorMessage('E',"ordernodes","renumbering of the mg failed");
    return (CMDERRORCODE);
  }

  /* now we reorder the nodes on the specified levels */
  for (level=fromLevel; level<=toLevel; level++)
  {
    theGrid = GRID_ON_LEVEL(theMG,level);

    UserWriteF(" [%d:",level);

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



/** \brief Implementation of \ref lexorderv. */
static INT LexOrderVectorsCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  INT i,res,level,fromLevel,toLevel;
  INT sign[DIM],order[DIM],which,xused,yused,zused,rused,pused,error,AlsoOrderMatrices,SpecialTreatSkipVecs,mode;
  char ord[3];

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"lexorderv","no open multigrid");
    return (CMDERRORCODE);
  }
  fromLevel = 0;
  toLevel   = TOPLEVEL(theMG);

  /* read ordering directions */
        #ifdef __TWODIM__
  res = sscanf(argv[0],expandfmt("lexorderv %2[rludIOPN]"),ord);
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
    PrintHelp("lexorderv",HELPITEM," (specify DIM chars out of 'rlud', 'IOPN' or 'rlbfud' resp.)");
    return(PARAMERRORCODE);
  }
  error = xused = yused = zused = rused = pused = FALSE;
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

    case 'u' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] =  1; break;
    case 'd' :
      if (yused) error = TRUE;
      yused = TRUE;
      order[i] = _Y_; sign[i] = -1; break;

#ifdef __THREEDIM__
    case 'b' :
      if (zused) error = TRUE;
      zused = TRUE;
      order[i] = _Z_; sign[i] =  1; break;
    case 'f' :
      if (zused) error = TRUE;
      zused = TRUE;
      order[i] = _Z_; sign[i] = -1; break;
#endif

                        #ifdef __TWODIM__

    /* polar coordiante directions */
    case 'I' :                          /* capital i */
      if (rused) error = TRUE;
      rused = TRUE;
      order[i] = 0; sign[i] =  1; break;

    case 'O' :
      if (rused) error = TRUE;
      rused = TRUE;
      order[i] = 0; sign[i] = -1; break;

    case 'P' :
      if (pused) error = TRUE;
      pused = TRUE;
      order[i] = 1; sign[i] =  1; break;

    case 'N' :
      if (pused) error = TRUE;
      pused = TRUE;
      order[i] = 1; sign[i] = -1; break;
                        #endif
    }
  if (error)
  {
    PrintHelp("lexorderv",HELPITEM," (bad combination of 'rludr' or 'rlbfud' resp.)");
    return(PARAMERRORCODE);
  }
  mode = OV_CARTES;
  if (rused || pused)
    if (!(rused && pused))
    {
      PrintHelp("lexorderv",HELPITEM," (bad combination of cartesian/polar direction)");
      return(PARAMERRORCODE);
    }
    else
      mode = OV_POLAR;

  /* check options */
  AlsoOrderMatrices = SpecialTreatSkipVecs = FALSE;
  which = GM_TAKE_SKIP | GM_TAKE_NONSKIP;
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

    case 'w' :
      which = 0;
      if (strchr(argv[i],'s')!=NULL)
        which |= GM_TAKE_SKIP;
      if (strchr(argv[i],'n')!=NULL)
        which |= GM_TAKE_NONSKIP;
      break;

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

    UserWriteF(" [%d:",level);

    if (LexOrderVectorsInGrid(theGrid,mode,order,sign,which,SpecialTreatSkipVecs,AlsoOrderMatrices)!=GM_OK)
    {
      PrintErrorMessage('E',"lexorderv","LexOrderVectorsInGrid failed");
      return (CMDERRORCODE);
    }
    UserWrite("ov]");
  }

  UserWrite("\n");

  return (OKCODE);
}


/** \brief Implementation of \ref shellorderv. */
static INT ShellOrderVectorsCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  VECTOR *seed;
  char option;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
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
  {
    l_setindex(theGrid);
    return (OKCODE);
  }
}


/** \brief Implementation of \ref orderv. */
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
      else if (strcmp(modestr,"FFLLCC")==0)
        mode = GM_FFLLCC;
      else if (strcmp(modestr,"FFLCLC")==0)
        mode = GM_FFLCLC;
      else if (strcmp(modestr,"CCFFLL")==0)
        mode = GM_CCFFLL;
      else
      {
        PrintHelp("orderv",HELPITEM," (you have to specify FFLLCC, FFLCLC, CCFFLL or FCFCLL as mode)");
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
    UserWrite("WARNING: no depency specified\n");
    if (dep_opt!=NULL)
    {
      UserWrite("WARNING: ignore specified options for dependency\n");
      dep_opt=NULL;
    }
  }

  if (dep!=NULL && dep_opt==NULL)
  {
    PrintErrorMessage('E',"orderv","the o option is mandatory if dopt specified");
    return (PARAMERRORCODE);
  }

  if (OrderVectors(theMG,levels,mode,PutSkipFirst,SkipPat,dep,dep_opt,cut)!=GM_OK)
  {
    PrintErrorMessage('E',"orderv","order vectors failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref revvecorder. */
static INT RevertVecOrderCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,from,to,l;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"revvecorder","no open multigrid");
    return (CMDERRORCODE);
  }

  from = to = CURRENTLEVEL(theMG);
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      from = 0;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("revvecorder",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  for (l=from; l<=to; l++)
  {
    RevertVecOrder(GRID_ON_LEVEL(theMG,l));
    UserWriteF(" [%d:rev]",l);
  }
  UserWrite("\n");

  return (OKCODE);
}


/** \brief Implementation of \ref lineorderv. */
static INT LineOrderVectorsCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char *dep,*dep_opt,*cut;
  INT i,levels,verboselevel;
  int iValue;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"lineorderv","no open multigrid");
    return (CMDERRORCODE);
  }

  levels           = GM_CURRENT_LEVEL;
  dep              = dep_opt = cut = NULL;
  verboselevel = 0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
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

    case 'a' :
      levels = GM_ALL_LEVELS;
      break;

    case 'v' :
      if (sscanf(argv[i],"v %d",&iValue)!=1)
      {
        PrintErrorMessage('E',"lineorderv","specify integer with v option");
        return (CMDERRORCODE);
      }
      verboselevel = iValue;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("lineorderv",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (dep==NULL)
  {
    PrintErrorMessage('E',"lineorderv","the d option is mandatory");
    return (PARAMERRORCODE);
  }

  if (dep_opt==NULL)
  {
    PrintErrorMessage('E',"lineorderv","the o option is mandatory");
    return (PARAMERRORCODE);
  }

  if (LineOrderVectors(theMG,levels,dep,dep_opt,cut,verboselevel)!=GM_OK)
  {
    PrintErrorMessage('E',"lineorderv","order vectors failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref setindex. */
static INT SetIndexCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"setindex","no open multigrid");
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


/** \brief Implementation of \ref find. */
static INT FindCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  NODE *theNode;
  VECTOR *theVector;
  ELEMENT *theElement;
  DOUBLE xc[DIM],tolc[DIM];
  INT i,j,select,isNode,isElement,isVector;

  /* following variables: keep type for sscanf */
  double x[DIM_MAX],tol;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"find","no open multigrid");
    return (CMDERRORCODE);
  }
  theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));

  /* check pararmeters */
  if (sscanf(argv[0],"find %lf %lf %lf",x,x+1,x+2)!=DIM)
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
      if (sscanf(argv[i],"n %lf",&tol)!=1)
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
        return (CMDERRORCODE);
      }
      isNode = TRUE;
      break;

    case 'v' :
      if (sscanf(argv[i],"v %lf",&tol)!=1)
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
        return (CMDERRORCODE);
      }
      isVector = TRUE;
      break;

    case 'e' :
      theElement = FindElementFromPosition(theGrid,xc);
      if (theElement==NULL)
      {
        PrintErrorMessage('W',"find","no element is matching");
        return (CMDERRORCODE);
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
      ListVector(theMG,theVector,FALSE,FALSE,LV_MOD_DEFAULT);

    if (isElement)
      ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref select. */
static INT SelectCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  NODE *theNode;
  ELEMENT *theElement;
  VECTOR *theVector;
  INT i,level;
  char c;

  /* following variables: keep type for sscanf */
  int id;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: SelectCommand(): me not in Context," \
                     " no selection of elements\n",me))
    return(OKCODE);
  }
        #endif

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
          PrintErrorMessageF('E',"select","node with ID %ld not found",(long)id);
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
          PrintErrorMessageF('E',"select","node with ID %ld is not in selection",(long)id);
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
          PrintErrorMessageF('E',"select","element with ID %ld not found",(long)id);
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
          PrintErrorMessageF('E',"select","element with ID %ld is not in selection",(long)id);
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

    case 'v' :
      if (sscanf(argv[i],"v %c %d",&c,&id)!=2)
      {
        PrintErrorMessage('E',"select","could not get +/- or ID");
        return (PARAMERRORCODE);
      }
      if (c=='+')
      {
        /* search vector */
        theVector = NULL;
        for (level=0; level<=TOPLEVEL(theMG); level++)
          if ((theVector=FindVectorFromIndex(GRID_ON_LEVEL(theMG,level),id))!=NULL)
            break;
        if (theVector==NULL)
        {
          PrintErrorMessageF('E',"select","vector with ID %ld not found",(long)id);
          return (CMDERRORCODE);
        }
        if (AddVectorToSelection(theMG,theVector)!=GM_OK)
        {
          PrintErrorMessage('E',"select","selecting the vector failed");
          return (CMDERRORCODE);
        }
      }
      else if (c=='-')
      {
        if (SELECTIONMODE(theMG)==vectorSelection)
          for (i=0; i<SELECTIONSIZE(theMG); i++)
          {
            theVector = (VECTOR *)SELECTIONOBJECT(theMG,i);
            if (ID(theVector)==id)
              break;
          }
        if (theVector==NULL)
        {
          PrintErrorMessageF('E',"select","vector with ID %ld is not in selection",(long)id);
          return (CMDERRORCODE);
        }
        if (RemoveVectorFromSelection(theMG,theVector)!=GM_OK)
        {
          PrintErrorMessage('E',"select","removing the vector failed");
          return (CMDERRORCODE);
        }
      }
      else
      {
        PrintErrorMessage('E',"select","specify + or - with n option");
        return (PARAMERRORCODE);
      }
      break;

    case 'i' :
      if (SELECTIONSIZE(theMG)==0)
        UserWrite("nothing selected\n");
      else switch (SELECTIONMODE(theMG))
        {
        case elementSelection :
          UserWriteF("%d elements selected (use for example 'elist $s')\n",SELECTIONSIZE(theMG));
          break;
        case nodeSelection :
          UserWriteF("%d nodes selected (use for example 'nlist $s')\n",SELECTIONSIZE(theMG));
          break;
        case vectorSelection :
          UserWriteF("%d vectors selected (use for example 'vmlist $s')\n",SELECTIONSIZE(theMG));
          break;
        default :
          UserWrite("unknown selection type\n");
        }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("select",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  return (OKCODE);
}


/** \brief Implementation of \ref extracon. */
static INT ExtraConnectionCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  VECTOR *vec;
  MATRIX *mat;
  INT Delete,i,nextra,nc;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"extracon","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  Delete = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'd' :
      Delete = TRUE;
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

  nc = NC(theGrid);
        #ifdef ModelP
  nextra = UG_GlobalSumINT(nextra);
  nc = UG_GlobalSumINT(nc);
        #endif

  UserWriteF("%d extra connections on level %d (total %d)\n",
             (int)nextra,(int)CURRENTLEVEL(theMG),(int)NC(theGrid));

  SetStringValue(":extraconratio",nextra/((DOUBLE)nc));

  if (Delete)
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



/** \brief Implementation of \ref check. */
static INT CheckCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  INT checkgeom,checkalgebra,checklists,checkbvp,checknp;
        #ifdef ModelP
  INT checkif;
        #endif
  INT level,err,i;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"check","no open multigrid");
    return (CMDERRORCODE);
  }

  /* set default options */
  checkgeom = TRUE;
  checkalgebra = checklists = checkbvp = checknp = FALSE;
        #ifdef ModelP
  checkif = FALSE;
        #endif

  /* read options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      checkgeom = checkalgebra = checklists = checknp = TRUE;
                                #ifdef ModelP
      checkif = TRUE;
                                #endif
      break;

    case 'g' :
      checkgeom = TRUE;
      break;

    case 'c' :
      checkalgebra = TRUE;
      break;

    case 'l' :
      checklists = TRUE;
      break;

                        #ifdef ModelP
    case 'i' :
      checkif = TRUE;
      break;
                        #endif

    case 'b' :
      checkbvp = TRUE;
      break;

    case 'n' :
      checknp = TRUE;
      break;

    case 'w' :
      ListAllCWsOfAllObjectTypes(UserWriteF);
      break;

    default :
      if (!checknp) {
        sprintf(buffer,"(invalid option '%s')",argv[i]);
        PrintHelp("check",HELPITEM,buffer);
        return (PARAMERRORCODE);
      }
    }
  err = 0;

  /* check BVP if */
  if (checkbvp==TRUE)
    if (BVP_Check (MG_BVP(theMG)))
      err++;

  for (level=0; level<=TOPLEVEL(theMG); level++)
  {
    theGrid = GRID_ON_LEVEL(theMG,level);
    UserWriteF("[%d:",(int)level);

                #ifndef ModelP
    if (CheckGrid(theGrid,checkgeom,checkalgebra,checklists)!=GM_OK)
                #else
    if (CheckGrid(theGrid,checkgeom,checkalgebra,checklists,checkif)!=GM_OK)
                #endif
      err++;

    UserWrite("]\n");
  }
  UserWrite("\n");

  if (checknp)
    if (CheckNP(theMG,argc,argv))
      err++;

  if (err)
    return (CMDERRORCODE);
  else
    return (OKCODE);
}



/****************************************************************************/
/** \brief Calculate minimal and maximal angle of an element
 *
 * @param theMG - pointer to multigrid
 * @param theElement - pointer to Element
 *
   This function calculates the minimal and the maximal angle of elements
   and lists elements with angle \< or \> given angles.


 */
/****************************************************************************/

INT NS_DIM_PREFIX QualityElement (MULTIGRID *theMG, ELEMENT *theElement)
{
  INT error;

  min = 360.0;
  max = 0.0;

  if ((error=MinMaxAngle(theElement,&min,&max))!=GM_OK)
    return(error);
  minangle = MIN(min,minangle);
  maxangle = MAX(max,maxangle);
  if ((lessopt && (min<themin)) && (greateropt && (max>themax)))
  {
    UserWrite(minmaxtext);
    ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
    if (selectopt) AddElementToSelection(theMG,theElement);
  }
  else if (lessopt && (min<themin))
  {
    UserWrite(mintext);
    ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
    if (selectopt) AddElementToSelection(theMG,theElement);
  }
  else if (greateropt && (max>themax))
  {
    UserWrite(maxtext);
    ListElement(theMG,theElement,FALSE,FALSE,FALSE,FALSE);
    if (selectopt) AddElementToSelection(theMG,theElement);
  }

  return(0);
}

/** \brief Implementation of \ref quality. */
static INT QualityCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *theGrid;
  ELEMENT *theElement;
  INT error,i,fromE,toE,res,mode;

  /* following variables: keep type for sscanf */
  double angle;
  long f,t;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"quality","no open multigrid");
    return (CMDERRORCODE);
  }

  /* check options */
  lessopt = greateropt = selectopt = mode = FALSE;
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
      if (sscanf(argv[i],"< %lf",&angle)!=1)
      {
        PrintErrorMessage('E',"quality","could not get angle of < option");
        return (CMDERRORCODE);
      }
      themin = angle;
      break;

    case '>' :
      greateropt = TRUE;
      if (sscanf(argv[i],"> %lf",&angle)!=1)
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

    case 'S' :
      selectopt = TRUE;
      ClearSelection(theMG);
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

  error=0;
  switch (mode)
  {
  case DO_ID :
    for (theGrid=GRID_ON_LEVEL(theMG,0); theGrid!=NULL; theGrid=UPGRID(theGrid))
      for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
        if ((ID(theElement)>=fromE) && (ID(theElement)<=toE))
          if ((error=QualityElement(theMG,theElement))!=0)
            break;
    break;

  case DO_ALL :
    for (theGrid=GRID_ON_LEVEL(theMG,0); theGrid!=NULL; theGrid=UPGRID(theGrid))
      for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
        if ((error=QualityElement(theMG,theElement))!=0)
          break;
    break;

  case DO_SELECTION :
    if (SELECTIONMODE(theMG)==elementSelection)
      for (i=0; i<SELECTIONSIZE(theMG); i++)
        if ((error=QualityElement(theMG,(ELEMENT *)SELECTIONOBJECT(theMG,i)))!=0)
          break;
    break;

  default :
    PrintErrorMessage('E',"quality","specify one option of a, s or i");
    return (PARAMERRORCODE);
  }

  if (error)
  {
    PrintErrorMessage('E',"quality","error in QualityElement/MinMaxAngle");
    return (CMDERRORCODE);
  }

  UserWriteF(" min angle = %20.12f\n max angle = %20.12f\n",(float)minangle,(float)maxangle);

  return(OKCODE);
}

/** \brief Implementation of \ref fiflel. */
#ifdef __THREEDIM__
static INT FindFlippedElementsCommand(INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT verbose;

  theMG = GetCurrentMultigrid();
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"fiflel","no current multigrid");
    return (CMDERRORCODE);
  }

  /* verbose mode */
  verbose = ReadArgvOption("v",argc,argv);

  if(FindFlippedElements(theMG,verbose))
    return (CMDERRORCODE);

  return(OKCODE);
}
#endif

/** \brief Implementation of \ref makegrid. */
static INT MakeGridCommand  (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT i,Single_Mode,display;
  MESH *mesh;
  INT MarkKey;
#ifdef __TWODIM__
  CoeffProcPtr coeff;
  GG_ARG args;
  GG_PARAM params;
  INT smooth;
  long ElemID,m;
  int iValue;
  double tmp;
#endif
#if defined __THREEDIM__ && defined _NETGEN
  INT smooth, from, to, prism, save;
  DOUBLE h,vol_h;
  INT coeff;
#ifdef ModelP
  INT mprocs;
#endif
#endif

  /* get current multigrid */
  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"makegrid","no open multigrid");
    return (CMDERRORCODE);
  }
        #if defined ModelP && !defined __THREEDIM__
  if (me!=master)
  {
    if (FixCoarseGrid(theMG)) return (CMDERRORCODE);
    return (OKCODE);
  }
        #endif
  if ((CURRENTLEVEL(theMG)!=0)||(TOPLEVEL(theMG)!=0))
  {
    PrintErrorMessage('E',"MakeGridCommand","only for a multigrid with exactly one level a grid can be generated");
    RETURN(GM_ERROR);
  }
  MarkKey = MG_MARK_KEY(theMG);
  if (MG_COARSE_FIXED(theMG)) {
    MG_COARSE_FIXED(theMG) = FALSE;
    MarkTmpMem(MGHEAP(theMG),&MarkKey);
    MG_MARK_KEY(theMG) = MarkKey;
    if ((MGNDELEMPTRARRAY(theMG) =
           (ELEMENT***)GetTmpMem(MGHEAP(theMG),NDELEM_BLKS_MAX*sizeof(ELEMENT**),MarkKey))==NULL)
    {
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
      PrintErrorMessage('E',"makegrid","ERROR: could not allocate memory from the MGHeap");
      return (CMDERRORCODE);
    }
    for (i=0; i<NDELEM_BLKS_MAX; i++) MGNDELEMBLK(theMG,i) = NULL;
  }

  Single_Mode = 0;
  display = 0;

  /* check options */
#ifdef __TWODIM__
  args.doanimate = args.doupdate = args.dostep = args.equilateral = args.plotfront = args.printelem = args.doangle = args.doEdge = args.doAngle = args.doConstDel =  NO;
  args.doedge = YES;
  ElemID = -1;
#endif

  if (DisposeGrid(GRID_ON_LEVEL(theMG,0)))
  {
    UserWriteF("makegrid: cannot dispose coarse grid\n");
    DisposeMultiGrid(theMG);
    return (CMDERRORCODE);
  }

  if (CreateNewLevel(theMG,0)==NULL)
  {
    UserWriteF("makegrid: cannot create new level\n");
    DisposeMultiGrid(theMG);
    return (CMDERRORCODE);
  }
  mesh = BVP_GenerateMesh (MGHEAP(theMG),MG_BVP(theMG),argc,argv,MarkKey);
  if (mesh == NULL)
  {
    UserWriteF("makegrid: cannot generate boundary mesh\n");
    ReleaseTmpMem(MGHEAP(theMG),MarkKey);
    return (CMDERRORCODE);
  }
  else
    InsertMesh(theMG,mesh);

  if (mesh->nElements == NULL)
  {
#ifdef __TWODIM__
    coeff = NULL;
    params.h_global = .0;
    params.CheckCos = cos(PI/18.0);
    params.searchconst = 0.2;
    smooth = 5;
    for (i=1; i<argc; i++)
      switch (argv[i][0])
      {
      case 'a' :
        args.doanimate = YES;
        break;
      case 'u' :
        args.doupdate = YES;
        break;
      case 's' :
        args.dostep = YES;
        break;
      case 'f' :
        args.plotfront = YES;
        break;
      case 'p' :
        args.printelem = YES;
        break;
      case 'E' :
        args.equilateral = YES;
        break;
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
      case 'C' :
        args.doConstDel = YES;
        args.doedge = NO;
        break;
      case 'm' :
        if (sscanf(argv[i],"m %ld",&m)!=1)
        {
          PrintHelp("makegrid",HELPITEM," (could not read <element id>)");
          return (PARAMERRORCODE);
        }
        coeff = MG_GetCoeffFct(theMG,m);
        break;
      case 'h' :
        if (sscanf(argv[i],"h %lf",&tmp)!=1)
        {
          PrintHelp("makegrid",HELPITEM," (could not read <element id>)");
          return (PARAMERRORCODE);
        }
        if (tmp > 0) params.h_global = tmp;
        break;
      case 'A' :
        if (sscanf(argv[i],"A %lf",&tmp)!=1)
        {
          PrintHelp("makegrid",HELPITEM," (could not read <element id>)");
          return (PARAMERRORCODE);
        }
        if ((tmp > 0) && (tmp < 90)) params.CheckCos = cos(tmp*PI/180.0);
        break;
      case 'S' :
        if (sscanf(argv[i],"S %lf",&tmp)!=1)
        {
          PrintHelp("makegrid",HELPITEM," (could not read <element id>)");
          return (PARAMERRORCODE);
        }
        if ((tmp > 0) && (tmp < 1.0)) params.searchconst = tmp;
        break;
      case 'd' :
        if (sscanf(argv[i],"d %d",&iValue)!=1) break;
        Single_Mode = iValue;
        break;
      case 'D' :
        if (sscanf(argv[i],"D %d",&iValue)!=1) break;
        display = iValue;
        break;
      case 'g' :
        if (sscanf(argv[i],"g %d",&iValue)!=1) break;
        smooth = iValue;
        break;
      default :
        break;
      }
    params.epsi = params.h_global * 0.125;

    if (GenerateGrid(theMG, &args, &params, mesh, coeff, Single_Mode, display))
    {
      PrintErrorMessage('E',"makegrid","execution failed");
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
      return (CMDERRORCODE);
    }
    if (SmoothMultiGrid(theMG,smooth,GM_KEEP_BOUNDARY_NODES)!=GM_OK)
    {
      PrintErrorMessage('E',"makegrid","failed smoothing the multigrid");
      return (CMDERRORCODE);
    }
    if (CheckOrientationInGrid (GRID_ON_LEVEL(theMG,0)))
    {
      PrintErrorMessage('E',"makegrid","orientation wrong");
      return (CMDERRORCODE);
    }
#endif

#if defined __THREEDIM__ && defined _NETGEN
    if (ReadArgvINT("s",&smooth,argc,argv)) smooth = 0;
    if (ReadArgvDOUBLE("h",&h,argc,argv)) h = 1.0;
    if(h<0)
      if (ReadArgvINT("c",&coeff,argc,argv)) coeff = 0;

    if (ReadArgvINT("s",&save,argc,argv))
      save = 0;
    if (ReadArgvINT("p",&prism,argc,argv))
      prism = 0;
    if (ReadArgvINT("f",&from,argc,argv))
      from = 1;
    if(from<1)
      from = 1;
    if (ReadArgvINT("t",&to,argc,argv))
      to = theMG->theBVPD.nSubDomains;
    if(to>theMG->theBVPD.nSubDomains)
      to = theMG->theBVPD.nSubDomains;

    if (ReadArgvDOUBLE("v",&vol_h,argc,argv))
      vol_h = h;

                #ifdef ModelP && defined __THREEDIM__ && defined _NETGEN
    if (ReadArgvINT("x",&mprocs,argc,argv))
    {
      /* sequential */
      mprocs = 1;
    }
    else
    {
      INT nsub,ndiv,i,j;

      mprocs = MAX(1,MIN(procs,mprocs));
      nsub = theMG->theBVPD.nSubDomains/mprocs;
      ndiv = theMG->theBVPD.nSubDomains%mprocs;

      if (me<ndiv)
      {
        nsub++;
        from = nsub*me+1;
        to   = MIN(theMG->theBVPD.nSubDomains,nsub*me+nsub);
      }
      else
      {
        from = nsub*me+1+ndiv;
        to   = MIN(theMG->theBVPD.nSubDomains,nsub*me+nsub+ndiv);
      }

      UserWriteF("%5d: Parallel Triangulation of %d subdomains:\n",me,theMG->theBVPD.nSubDomains);
      UserWriteF("%5d: mprocs=%d nsub=%d from=%d to=%d\n",me,mprocs,nsub,from,to);
    }
                #endif

    if (GenerateGrid3d(theMG,mesh,vol_h,smooth,ReadArgvOption("d",argc,argv),coeff, from, to, prism, save, argc, argv))
    {
      PrintErrorMessage('E',"makegrid","execution failed");
      ReleaseTmpMem(MGHEAP(theMG),MarkKey);
      return (CMDERRORCODE);
    }
#endif
  }

  if (FixCoarseGrid(theMG)) return (CMDERRORCODE);
  InvalidatePicturesOfMG(theMG);
  InvalidateUgWindowsOfMG(theMG);

        #if defined ModelP && defined __THREEDIM__ && defined _NETGEN
  if (mprocs > 1)
    if (IdentifySDGrid(theMG,mprocs)) return (CMDERRORCODE);
        #endif

  return (OKCODE);
}


/** \brief Implementation of \ref status. */
static INT StatusCommand  (INT argc, char **argv)
{
  MULTIGRID       *theMG;
  INT i,grid,green,load,verbose;

  /* get current multigrid */
  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"status command","no open multigrid");
    return (CMDERRORCODE);
  }
  grid = green = load = 0;
  verbose = 1;

  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      grid = 1;
      green = 1;
                                #ifdef ModelP
      load = 1;
                                #endif
      break;
    case 'g' :
      green = 1;
      break;
                        #ifdef ModelP
    case 'l' :
      load = 1;
      sscanf(argv[i],"l %d",&load);
      break;
                        #endif
    case 'm' :
      grid = 1;
      break;
    default :
      break;
    }

  if (MultiGridStatus(theMG,grid,green,load,verbose) != 0)
  {
    PrintErrorMessage('E',"GridStatus()","execution failed");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}


#if defined(CAD) && defined(__THREEDIM__)
/** \brief Implementation of \ref cadconvert. */
static INT CADGridConvertCommand(INT argc, char **argv)
{
  MULTIGRID *theMG;
  char CADOutputFileName[NAMESIZE];

  char Format[NAMESIZE];
  char *theFormat;

  unsigned long heapSize;
  INT i,hopt;

  /* get CADfile name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" cadconvert %",NAMELENSTR,"[ -~]")),CADOutputFileName)!=1) || (strlen(CADOutputFileName)==0))
    sprintf(CADOutputFileName,"untitled-%d",(int)untitledCounter++);

  /* get problem, domain and format */
  theFormat = NULL;
  heapSize = 0;
  hopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
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

  if (!(hopt))
  {
    PrintHelp("new",HELPITEM," (the d, p, f and h arguments are mandatory)");
    return(PARAMERRORCODE);
  }

  /* check options */
  theMG = ConvertCADGrid(theFormat, CADOutputFileName, heapSize);
  if (theMG == NULL)
  {
    PrintErrorMessage('E',"cadconvert","execution failed");
    return (CMDERRORCODE);
  }

  if (SetCurrentMultigrid(theMG)!=0)
    return (CMDERRORCODE);



  return (OKCODE);
}
#endif


/** \brief Implementation of \ref grape. */
static INT CallGrapeCommand (INT argc, char **argv)
{
  MULTIGRID *theCurrMG;

  /* see if multigrid exists */
  theCurrMG = currMG;
  if (theCurrMG==NULL)
  {
    UserWrite("cannot call grape without multigrid\n");
    return (CMDERRORCODE);
  }

  /* call grape */
  if (CallGrape(theCurrMG)) return (CMDERRORCODE);

  return (OKCODE);
}


#ifdef _COVISE

/** \brief Implementation of \ref covise. */
static INT CoviseCommand (INT argc, char **argv)
{
  MULTIGRID *theCurrMG;
  char hostname[NAMESIZE];

  /* see if multigrid exists */
  theCurrMG = currMG;
  if (theCurrMG==NULL)
  {
    UserWrite("cannot start covise without multigrid\n");
    return (CMDERRORCODE);
  }


  if (sscanf(argv[0],"covise %s",hostname)!=1)
  {
    UserWrite("covise: specify hostname to connect to!\n");
    return(PARAMERRORCODE);
  }

  /* call covise */
  if (ConnectCovise(theCurrMG, hostname)) return (CMDERRORCODE);

  return (OKCODE);
}
#endif



/** \brief Implementation of \ref screensize. */
static INT ScreenSizeCommand (INT argc, char **argv)
{
  INT size[2];

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  NO_OPTION_CHECK(argc,argv);

  if (GetScreenSize(size)==FALSE)
  {
    PrintErrorMessage('W',"screensize","there is no monitor");
    return (OKCODE);
  }

  UserWriteF(" screen width: %d, screen height: %d\n",(int)size[0], (int)size[1]);

  if (    (SetStringValue(":screensize:width",size[0])!=0)
          ||      (SetStringValue(":screensize:height",size[1])!=0))
  {
    PrintErrorMessage('E',"screensize","could not set :screensize:width or :screensize:height");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref openwindow. */
static INT OpenWindowCommand (INT argc, char **argv)
{
  OUTPUTDEVICE *theOutDev;
  UGWINDOW *theWin;
  char devname[NAMESIZE],winname[NAMESIZE];
  INT i,rename,res;
  int ropt;

  /* following variables: keep type for sscanf */
  int x,y,w,h;
        #ifdef ModelP
  if (!CONTEXT(me))
    PRINTDEBUG(ui,0,("%2d: OpenWindowCommand(): me not in Context," \
                     "  no window structure allocated\n",me))
    /* Note: it is forbidden that a part of the processors leave this function premature
       since in the subsequent CreateUgWindow communication occures, whis requires
       all processors to participate. Thus ensure that all or none processor leaves.
       Christian Wrobel 990329 */
    for(i=0; i<procs; i++)
      if(!CONTEXT(i))                   /* each procesor has the same context; thus no communication is required */
        return(CMDERRORCODE);
        #endif

  /* scan parameters */
  if (sscanf(argv[0],"openwindow %d %d %d %d",&x,&y,&w,&h)!=4)
  {
    PrintHelp("openwindow",HELPITEM," could not get all mandatory parameters");
    return (PARAMERRORCODE);
  }

  /* check options */
  theOutDev  = GetDefaultOutputDevice();
  winname[0] = '\0';
  rename = 0;
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
        PrintErrorMessageF('E',"openwindow","there is no device named '%s'",devname);
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

    case 'r' :
      res = sscanf(argv[i]," r %d",&ropt);
      if (res==0 || (res==1 && ropt==1)) rename = 1;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("openwindow",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (strlen(winname)==0)
    sprintf(winname,"window%d",(int) wincounter++);


        #ifdef ModelP
  if (me == master)
        #endif
  if (theOutDev==NULL) {
    PrintErrorMessage('E',"openwindow","no output device");
    return (PARAMERRORCODE);
  }

  if ((theWin=CreateUgWindow(theOutDev,winname,rename,x,y,w,h))==NULL)
  {
    PrintErrorMessage('E',"openwindow","failed to open a window");
    return (CMDERRORCODE);
  }

  SetCurrentUgWindow(theWin);

  return (OKCODE);
}


/** \brief Implementation of \ref closewindow. */
static INT CloseWindowCommand (INT argc, char **argv)
{
  UGWINDOW *theWin, *currWin;
  PICTURE *thePic,*currPic;
  char winname[NAMESIZE];
  INT i,aopt;

        #ifdef ModelP
  /*if (me!=master) return (OKCODE); since OpenWindow allocates structures for all procs, CloseWindow must deallocate for all those (and not only for master). Christian Wrobel 990326 */
  if (!CONTEXT(me))
    PRINTDEBUG(ui,0,("%2d: CloseWindowCommand(): me not in Context," \
                     "  no window structure closed\n",me))
    /* Note: it is forbidden that a part of the processors leave this function premature
       since all allocated structures must be deallocated.
       Thus ensure that all or none processor leaves.
       Christian Wrobel 990329 */
    for(i=0; i<procs; i++)
      if(!CONTEXT(i))                   /* each procesor has the same context; thus no communication is required */
        return(CMDERRORCODE);
        #endif

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
        PrintErrorMessageF('W',"closewindow","there is no window named '%s'",winname);
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


/** \brief Implementation of \ref setcurrwindow. */
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


/** \brief Implementation of \ref drawtext. */
static INT DrawTextCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  char winname[NAMESIZE],text[NAMESIZE];
  COORD_POINT pos;
  INT i,mode,centeropt,size;
  double x,y;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  theWin = GetCurrentUgWindow();
  if (theWin==NULL)
  {
    PrintErrorMessage('E',"drawtext","there's no window to draw text");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],expandfmt(CONCAT3("drawtext %lf %lf %",NAMELENSTR,"[ -~]")),&x,&y,text)!=3)
  {
    PrintErrorMessage('E',"drawtext","specify position with two integers and then the text");
    return (CMDERRORCODE);
  }
  pos.x = x; pos.y = y;

  /* check options */
  centeropt = FALSE;
  mode = TEXT_REGULAR;
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
        PrintErrorMessageF('E',"drawtext","there is no window named '%s'",winname);
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

    case 'm' :
      if (strstr(argv[i],"reg")!=NULL)
        mode = TEXT_REGULAR;
      else if (strstr(argv[i],"inv")!=NULL)
        mode = TEXT_INVERSE;
      else if (strstr(argv[i],"ind")!=NULL)
        mode = TEXT_INDEXED;
      break;

    case 'c' :
      centeropt = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("drawtext",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  DrawWindowText(theWin,pos,text,size,centeropt,mode);

  return (OKCODE);
}


/** \brief Implementation of \ref openpicture. */
static INT OpenPictureCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  PICTURE *thePic;
  char picname[NAMESIZE],winname[NAMESIZE];
  INT local_LL[2],local_UR[2];
  INT i,sopt;

  /* following variables: keep type for sscanf */
  int h,v,dh,dv;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: OpenPictureCommand(): me not in Context,"\
                     " no picture stucture allocated\n",me))
    return(OKCODE);
  }
        #endif

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
        PrintErrorMessageF('E',"openpicture","there is no window named '%s'",winname);
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


/** \brief Implementation of \ref openppic. */
static INT OpenPlacedPicturesCommand (INT argc, char **argv)
{
  INT i,qopt,ropt,nPic,sopt,wopt,res,rename;
  PLACEMENT_TASK task;
  int iValue,v,h,dv,dh;
  OUTPUTDEVICE *theOutDev;
  char devname[NAMESIZE],qarray[NAMESIZE],rarray[NAMESIZE],buffer[NAMESIZE];
  UGWINDOW *theWin;

  /* get number of pictures */
  if (sscanf(argv[0],"openppic %d",&iValue)!=1)
  {
    PrintErrorMessage('E',"openppic","specify number of pictures with n option");
    return (PARAMERRORCODE);
  }
  nPic = iValue;

  /* check options */
  theOutDev  = GetDefaultOutputDevice();
  rename=sopt=wopt=qopt=ropt=0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'w' :
      wopt=1;
      if (sscanf(argv[i],expandfmt(CONCAT3("w %",NAMELENSTR,"[a-zA-Z0-9_.]")),task.win_name)!=1)
      {
        PrintErrorMessage('E',"openppic","specify a window name with w option");
        return (PARAMERRORCODE);
      }
      break;

    case 'q' :
      qopt=1;
      if (sscanf(argv[i],expandfmt(CONCAT3("q %",NAMELENSTR,"[a-zA-Z0-9_:]")),qarray)!=1)
      {
        PrintErrorMessage('E',"openppic","specify an array name with q option");
        return (PARAMERRORCODE);
      }
      break;

    case 'r' :
      ropt=1;
      if (sscanf(argv[i],expandfmt(CONCAT3("r %",NAMELENSTR,"[a-zA-Z0-9_:]")),rarray)!=1)
      {
        PrintErrorMessage('E',"openppic","specify an array name with r option");
        return (PARAMERRORCODE);
      }
      break;

    case 'R' :
      res = sscanf(argv[i]," R %d",&iValue);
      if (res==0 || (res==1 && iValue==1)) rename = 1;
      break;

    case 's' :
      sopt = 1;
      if (sscanf(argv[i],"s %d %d %d %d",&h,&v,&dh,&dv)!=4)
      {
        PrintErrorMessage('E',"openpicture","specify h, v, dh, dv with s option");
        return (PARAMERRORCODE);
      }
      task.winLL[0] = h;
      task.winLL[1] = v;
      task.winUR[0] = h+dh;
      task.winUR[1] = v+dv;
      break;

    case 'd' :
      if (sscanf(argv[i],expandfmt(CONCAT3("d %",NAMELENSTR,"[a-zA-Z0-9_-]")),devname)!=1)
      {
        PrintErrorMessage('E',"openppic","specify device name with d option");
        return (PARAMERRORCODE);
      }
      if ((theOutDev=GetOutputDevice(devname))==NULL)
      {
        PrintErrorMessageF('E',"openppic","there is no device named '%s'",devname);
        return (PARAMERRORCODE);
      }
      break;

    default :
      PrintErrorMessage('E',"openppic","unknown option");
      return (PARAMERRORCODE);
    }

  /* check if size initialized */
  if (!sopt)
  {
    PrintErrorMessage('E',"openppic","size not specified");
    return (PARAMERRORCODE);
  }

  /* check if window name initialized */
  if (!wopt)
  {
    PrintErrorMessage('E',"openppic","window name not specified");
    return (PARAMERRORCODE);
  }

  /* check if aspect-ratio-array name initialized */
  if (!qopt)
  {
    PrintErrorMessage('E',"openppic","q-array name not specified");
    return (PARAMERRORCODE);
  }

  /* check if ratio-array name initialized */
  if (!ropt)
  {
    PrintErrorMessage('E',"openppic","r-array name not specified");
    return (PARAMERRORCODE);
  }

  task.n = nPic;

  /* set names of pictues */
  for (i=0; i<nPic; i++)
  {
    sprintf(task.pic_name[i],"pic_%d",(int)i);
    sprintf(buffer,"%s%d",qarray,i);
    if (GetStringValueDouble (buffer,&(task.aspect_ratio[i])))
    {
      PrintErrorMessage('E',"openppic","q-array entry not found");
      return (PARAMERRORCODE);
    }
    sprintf(buffer,"%s%d",rarray,i);
    if (GetStringValueDouble (buffer,&(task.rel_size[i])))
    {
      PrintErrorMessage('E',"openppic","r-array entry not found");
      return (PARAMERRORCODE);
    }
  }

  /* check device */
  if (theOutDev==NULL)
  {
    PrintErrorMessage('E',"openppic","cannot find outputdevice");
    return (PARAMERRORCODE);
  }

  /* place pictures */
  theWin = OpenPlacedPictures(theOutDev,&task,rename);
  if (theWin==NULL) return (PARAMERRORCODE);
  SetCurrentUgWindow(theWin);

  return (OKCODE);
}


/** \brief Implementation of \ref closepicture. */
static INT ClosePictureCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  PICTURE *thePic,*NextPic;
  char picname[NAMESIZE],winname[NAMESIZE];
  INT i,aopt,wopt;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

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
        PrintErrorMessageF('E',"closepicture","there is no window named '%s'",winname);
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
        PrintErrorMessageF('E',"closepicture","there is no picture named '%s'",picname);
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

  if (DisposePicture(thePic)!=0)
  {
    PrintErrorMessage('E',"closepicture","could not close the picture");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref setcurrpicture. */
static INT SetCurrentPictureCommand (INT argc, char **argv)
{
  UGWINDOW *theWin;
  PICTURE *thePic;
  char picname[NAMESIZE], winname[NAMESIZE];
  INT i;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: OpenPictureCommand(): me not in Context,"\
                     " no picture stucture allocated\n",me))
    return(OKCODE);
  }
        #endif

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
        PrintErrorMessageF('E',"setcurrpicture","there is no window named '%s'",winname);
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


/** \brief Implementation of \ref picwin. */
static INT PictureWindowCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('W',"picwin","there's no picture to move");
    return (OKCODE);
  }
  if (ErasePicture(thePic)) return (CMDERRORCODE);
  if (MovePictureToNewWindow(thePic))
  {
    PrintErrorMessage('E',"picwin","failed to create a new window for the picture");
    return (CMDERRORCODE);
  }

  SetCurrentUgWindow(PIC_UGW(thePic));
  SetCurrentPicture(thePic);

  return (OKCODE);
}


/** \brief Implementation of \ref clearpicture. */
static INT ClearPictureCommand (INT argc, char **argv)
{
  PICTURE *thePicture;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

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


/** \brief Implementation of \ref picframe. */
static INT PicFrameCommand (INT argc, char **argv)
{
        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  NO_OPTION_CHECK(argc,argv);

  if (strchr(argv[0],'0')!=NULL)
    SetDoFramePicture (NO);
  else
    SetDoFramePicture (YES);

  return (OKCODE);
}


/** \brief Implementation of \ref setview. */
static INT SetViewCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  VIEWEDOBJ *theViewedObj;
  DOUBLE *viewPoint,*targetPoint,*xAxis;
  DOUBLE vP[3],tP[3],xA[3];
  DOUBLE PlanePoint[3],PlaneNormal[3];
  DOUBLE *CutPoint,*CutNormal,*scaleptr, scale[3];
  INT *perspective;
  INT per;
  INT i,j,veclen,res,RemoveCut;

  /* following variables: keep type for sscanf */
  double x[3],help[3];

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: SetViewCommand(): me not in Context,"\
                     " view not specified",me))
    return(OKCODE);
  }
        #endif

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
  CutPoint = CutNormal = NULL;
  RemoveCut = NO;
  scaleptr = NULL;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'o' :
      if (VO_DIM(theViewedObj)!=TYPE_3D)
      {
        PrintErrorMessage('E',"setview","the o option applies ONLY with 3D objects");
        return (PARAMERRORCODE);
      }
      if (sscanf(argv[i],"o %lf %lf %lf",x,x+1,x+2)!=veclen)
      {
        PrintErrorMessageF('E',"setview","o option: %d coordinates required for a %dD object",(int)veclen,(int)veclen);
        return (PARAMERRORCODE);
      }
      for (j=0; j<veclen; j++)
        vP[j] = x[j];
      viewPoint = vP;
      break;

    case 't' :
      if (sscanf(argv[i],"t %lf %lf %lf",x,x+1,x+2)!=veclen)
      {
        PrintErrorMessageF('E',"setview","t option: %d coordinates required for a %dD object",(int)veclen,(int)veclen);
        return (PARAMERRORCODE);
      }
      for (j=0; j<veclen; j++)
        tP[j] = x[j];
      targetPoint = tP;
      break;

    case 's' :
      if (sscanf(argv[i],"s %lf %lf %lf",x,x+1,x+2)!=veclen)
      {
        PrintErrorMessageF('E',"setview","s option: %d scalings required for a %dD object",(int)veclen,(int)veclen);
        return (PARAMERRORCODE);
      }
      for (j=0; j<veclen; j++)
        scale[j] = x[j];
      scaleptr = scale;
      break;

    case 'x' :
      if (sscanf(argv[i],"x %lf %lf %lf",x,x+1,x+2)!=veclen)
      {
        PrintErrorMessageF('E',"setview","x option: %d coordinates required for a %dD object",(int)veclen,(int)veclen);
        return (PARAMERRORCODE);
      }
      for (j=0; j<veclen; j++)
        xA[j] = x[j];
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


    /* capitals for cut definition */

    case 'C' :
      if (!PO_USESCUT(PIC_PO(thePic)))
      {
        PrintErrorMessage('E',"setview","plot object does not use a cut");
        return(PARAMERRORCODE);
      }

      /* set default cut */
      V3_COPY(PO_MIDPOINT(PIC_PO(thePic)),PlanePoint);
      V3_CLEAR(PlaneNormal);
      CutPoint = PlanePoint;
      CutNormal = PlaneNormal;
      break;

    case 'R' :
      if (!PO_USESCUT(PIC_PO(thePic)))
      {
        PrintErrorMessage('E',"setview","plot object does not use a cut");
        return(PARAMERRORCODE);
      }

      RemoveCut = YES;
      break;

    case 'P' :
      if (!PO_USESCUT(PIC_PO(thePic)))
      {
        PrintErrorMessage('E',"setview","plot object does not use a cut");
        return(PARAMERRORCODE);
      }
      res = sscanf(argv[i],"P %lg %lg %lg",help, help+1, help+2);
      if (res!=3)
      {
        PrintErrorMessage('E',"setview","specify three values for cut plane point");
        return(PARAMERRORCODE);
      }
      V3_COPY(help,PlanePoint);
      CutPoint = PlanePoint;
      break;

    case 'N' :
      if (!PO_USESCUT(PIC_PO(thePic)))
      {
        PrintErrorMessage('E',"setview","plot object does not use a cut");
        return(PARAMERRORCODE);
      }
      res = sscanf(argv[i],"N %lg %lg %lg",help, help+1, help+2);
      if (res!=3)
      {
        PrintErrorMessage('E',"setview","specify three values for cut normal point");
        return(PARAMERRORCODE);
      }
      V3_COPY(help,PlaneNormal);
      CutNormal = PlaneNormal;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("setview",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (SetView(thePic,viewPoint,targetPoint,xAxis,perspective,RemoveCut,CutPoint,CutNormal,scaleptr)!=0)
  {
    PrintErrorMessage('E',"setview","error during SetView");
    return (CMDERRORCODE);
  }

  /* picture has changed */
  if (InvalidatePicture(thePic))
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref vdisplay. */
static INT DisplayViewCommand (INT argc, char **argv)
{
  PICTURE *thePic;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"vdisplay","there's no current picture");
    return (CMDERRORCODE);
  }

  /* check options */
  switch (argc)
  {
  case 1 :
    if (DisplayViewOfViewedObject(thePic))
    {
      PrintErrorMessage('E',"vdisplay","error during DisplayView");
      return (CMDERRORCODE);
    }
    break;
  case 2 :
    if (argv[1][0]!='s')
    {
      sprintf(buffer,"(invalid option '%s')",argv[1]);
      PrintHelp("vdisplay",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }
    PrintViewSettings(thePic);
    break;
  default :
    PrintErrorMessage('E',"vdisplay","too many options");
    return (CMDERRORCODE);
  }

  return(OKCODE);
}



/** \brief Implementation of \ref cpview. */
static INT CopyViewCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  INT i,all,cut;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"cpview","there's no current picture");
    return (CMDERRORCODE);
  }
  all = cut = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      all = TRUE;
      break;

    case 'c' :
      cut = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("cpview",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (CopyView(thePic,all,cut))
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref walk. */
static INT WalkCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  VIEWEDOBJ *theViewedObj;
  DOUBLE dx[3];
  INT i,veclen;

  /* following variables: keep type for sscanf */
  double x[3];

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

  if (sscanf(argv[0],"walk %lf %lf %lf",x,x+1,x+2)!=veclen)
  {
    PrintErrorMessageF('E',"walk","%d coordinates required for a %dD object",(int)veclen,(int)veclen);
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
/** \brief Implementation of \ref walkaround

   WalkAroundCommand - let the observer walk on a sphere around the target point

 * @param argc - number of arguments (incl. its own name)
 * @param argv - array of strings giving the arguments

   This function lets the observer walk on a sphere around
   the target point in the current picture. (3D ONLY)
   It lets the observer walk on a sphere around the target point in the current picture.

   walkaround \<viewplane angle> \<rotation angle>

   .  \<viewplane~angle>      - this angle runs in the view plane math pos from the x-axis and defines
   .n                        - together with the target-observer direction a plane

   .  \<rotation~angle>       - the observer will be rotated around the target point in the above plane

   RETURN VALUE:
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT WalkAroundCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  VIEWEDOBJ *theViewedObj;

  /* following variables: keep type for sscanf */
  double dirAngle,Angle;

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

  if (sscanf(argv[0],"walkaround %lf %lf",&dirAngle,&Angle)!=2)
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

/** \brief Implementation of \ref zoom. */
static INT ZoomCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  /* following variables: keep type for sscanf */
  double factor;

  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"zoom","there's no current picture");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"zoom %lf",&factor)!=1)
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


/** \brief Implementation of \ref drag. */
static INT DragCommand (INT argc, char **argv)
{
  PICTURE *thePic;

  /* following variables: keep type for sscanf */
  double dx,dy;

  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"drag","there's no current picture");
    return (CMDERRORCODE);
  }

  if (sscanf(argv[0],"drag %lf %lf",&dx,&dy)!=2)
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


/** \brief Implementation of \ref rotate. */
static INT RotateCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  VIEWEDOBJ *theVO;
  DOUBLE ex,ey,norm;

  /* following variables: keep type for sscanf */
  double angle;

  NO_OPTION_CHECK(argc,argv);

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"rotate","there's no current picture");
    return (CMDERRORCODE);
  }
  theVO=PIC_VO(thePic);

  if (sscanf(argv[0],"rotate %lf",&angle)!=1)
  {
    V_DIM_EUKLIDNORM(VO_PXD(theVO),norm);
    if (norm==0.0) return(CMDERRORCODE);
    ex=VO_PXD(theVO)[DIM-1]/norm;
    V_DIM_EUKLIDNORM(VO_PYD(theVO),norm);
    if (norm==0.0) return(CMDERRORCODE);
    ey=VO_PYD(theVO)[DIM-1]/norm;
    if (ex==0.0 && ey==0.0) return(CMDERRORCODE);
    angle=-atan2(ex,ey);
    if (ey*cos(angle)<ex*sin(angle))
      angle+=PI;
  }
  else
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


/** \brief Implementation of \ref textfac. */
static INT TextFacCommand (INT argc, char **argv)
{
  double Value;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  NO_OPTION_CHECK(argc,argv);

  if (sscanf(argv[0],"textfac %lf",&Value)!=1)
  {
    PrintErrorMessage('E',"textfac","specify a factor");
    return (PARAMERRORCODE);
  }

  SetTextFactor(Value);

  InvalidatePicturesOfMG(currMG);

  return (OKCODE);
}


/** \brief Implementation of \ref linefac. */
static INT LineFacCommand (INT argc, char **argv)
{
  double Value;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  NO_OPTION_CHECK(argc,argv);

  if (sscanf(argv[0],"linefac %lf",&Value)!=1)
  {
    PrintErrorMessage('E',"linefac","specify a factor");
    return (PARAMERRORCODE);
  }

  SetLineFactor(Value);

  InvalidatePicturesOfMG(currMG);

  return (OKCODE);
}



/** \brief Implementation of \ref setplotobject. */
static INT SetPlotObjectCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  MULTIGRID *theMG;
  char potname[NAMESIZE],*thePlotObjTypeName;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: SetPlotObjectCommand(): me not in Context," \
                     " plot object not specified\n",me))
    return(OKCODE);
  }
        #endif

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
    UserWriteF(" picture '%s' and multigrid '%s' coupled\n",ENVITEM_NAME(thePic),ENVITEM_NAME(theMG));

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


/** \brief Implementation of \ref polist. */
static INT PlotObjectListCommand (INT argc, char **argv)
{
  PICTURE *thePic;

        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

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


/** \brief Implementation of \ref plot. */
static INT PlotCommand (INT argc, char **argv)
{
  UGWINDOW *theUgW;
  PICTURE *thePic,*currPic;
  INT i,OrderStrategy,all,bullet,noframe;
  DOUBLE zOffsetFactor;

  /* scan for options */

  OrderStrategy = 0;
  all = bullet = noframe = NO;
  zOffsetFactor = 1.0;

  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'o' :
      sscanf(argv[i],"o %d", &OrderStrategy);
      break;
    case 'a' :
      all = YES;
      break;
    case 'b' :
      bullet = YES;
      sscanf(argv[i],"b %lf", &zOffsetFactor);
      break;
    case 'n' :
      noframe = YES;
      break;

    default :
      break;
    }
  if (SetOrderStrategy(OrderStrategy)!=0)
  {
    PrintErrorMessage('E',"plot","invalid order mode");
    return (CMDERRORCODE);
  }

  if (all)
  {
    currPic = GetCurrentPicture();

    /* plot all pictures in all windows */
    for (theUgW=GetFirstUgWindow(); theUgW!=NULL; theUgW=GetNextUgWindow(theUgW))
      for (thePic=GetFirstPicture(theUgW); thePic!=NULL; thePic=GetNextPicture(thePic))
      {
        if (bullet) {
          if (BulletDrawUgPicture(thePic, zOffsetFactor)!=0)
          {
            PrintErrorMessage('E',"plot","error during WorkOnPicture");
            return (CMDERRORCODE);
          }
        }
        else {
          if (DrawUgPicture(thePic)!=0)
          {
            PrintErrorMessage('E',"plot","error during WorkOnPicture");
            return (CMDERRORCODE);
          }
        }

                                #ifdef ModelP
        if (me == master)
                                #endif

        if (thePic==currPic)
          DrawPictureFrame(thePic,WOP_ACTIVE);
        else
          DrawPictureFrame(thePic,WOP_NOT_ACTIVE);
      }

    return (OKCODE);
  }

  /* current picture */
  thePic = GetCurrentPicture();
  if (thePic==NULL)
  {
    PrintErrorMessage('E',"plot","there's no current picture");
    return (CMDERRORCODE);
  }

  if (bullet) {
    if (BulletDrawUgPicture(thePic, zOffsetFactor)!=0)
    {
      PrintErrorMessage('E',"plot","error during WorkOnPicture");
      return (CMDERRORCODE);
    }
  }
  else {
    if (DrawUgPicture(thePic)!=0)
    {
      PrintErrorMessage('E',"plot","error during WorkOnPicture");
      return (CMDERRORCODE);
    }
  }

        #ifdef ModelP
  if (me == master)
        #endif

  /* picture is current */
  if (noframe==NO)
    DrawPictureFrame(thePic,WOP_ACTIVE);

  return (OKCODE);
}


/** \brief Implementation of \ref findrange. */
static INT FindRangeCommand (INT argc, char **argv)
{
  PICTURE *thePic;
  WORK myWork,*theWork;
  INT i,sym,put;
  DOUBLE min,max;

  /* following variables: keep type for sscanf */
  double zoom;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: FindRangeCommand(): me not in Context,"\
                     " range found\n",me))
    return(OKCODE);
  }
        #endif

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
      if (sscanf(argv[i],"z %lf",&zoom)!=1)
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
        #ifdef ModelP
  Broadcast(&min,sizeof(double));
  W_FINDRANGE_WORK(theWork)->min = min;
        #endif

  max = W_FINDRANGE_WORK(theWork)->max;
        #ifdef ModelP
  Broadcast(&max,sizeof(double));
  W_FINDRANGE_WORK(theWork)->max = max;
        #endif

  UserWriteF(" FR_min = %20.16e\n FR_max = %20.16e\n",min,max);

  if (put)
    if (InvalidatePicture(thePic))
      return (CMDERRORCODE);

  if (    (SetStringValue(":findrange:min",min)!=0)
          ||      (SetStringValue(":findrange:max",max)!=0))
  {
    PrintErrorMessage('E',"findrange","could not set :findrange:min or :findrange:max");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref setcurrmg. */
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


/** \brief Implementation of \ref updateDoc. */
static INT UpdateDocumentCommand (INT argc, char **argv)
{
        #ifdef ModelP
  if (me!=master) return (OKCODE);
        #endif

  NO_OPTION_CHECK(argc,argv);

  InvalidatePicturesOfMG(currMG);
  InvalidateUgWindowsOfMG(currMG);

  return (OKCODE);
}


/** \brief Implementation of \ref rotmode. */
static INT RotModeCommand (INT argc, char **argv)
{
  INT mode;

  NO_OPTION_CHECK(argc,argv);

  if (strchr(argv[0],'E')!=NULL)
    mode = ROTMODE_EULER;
  else if (strchr(argv[0],'S')!=NULL)
    mode = ROTMODE_SPHERE;
  else
  {
    PrintHelp("rotmode",HELPITEM," (specify Euler or Sphere)");
    return(PARAMERRORCODE);
  }

  SetRotMode(mode);

  return (OKCODE);
}


/** \brief Implementation of \ref setpalette. */
static INT SetPaletteCommand (INT argc, char **argv)
{
  OUTPUTDEVICE *theOutDev;
  INT i,palette;
  char plt,devname[NAMESIZE];

  if (sscanf(argv[0],"setpalette %c",&plt)!=1)
  {
    PrintHelp("setpalette",HELPITEM," (specify c|bw|g)");
    return(PARAMERRORCODE);
  }
  switch (plt)
  {
  case 'c' : palette = COLOR_PALETTE; break;
  case 'b' : palette = BLACK_WHITE_PALETTE; break;
  case 'g' : palette = GRAY_PALETTE; break;
  default :
    PrintHelp("setpalette",HELPITEM," (specify c|bw|g)");
    return(PARAMERRORCODE);
  }

  /* check options */
  theOutDev  = GetDefaultOutputDevice();
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'd' :
      if (sscanf(argv[i],expandfmt(CONCAT3("d %",NAMELENSTR,"[a-zA-Z0-9_-]")),devname)!=1)
      {
        PrintErrorMessage('E',"setpalette","specify device name with d option");
        return (PARAMERRORCODE);
      }
      if ((theOutDev=GetOutputDevice(devname))==NULL)
      {
        PrintErrorMessageF('E',"setpalette","there is no device named '%s'",devname);
        return (PARAMERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("setpalette",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  if (UgSetPalette(theOutDev,palette))
    REP_ERR_RETURN(CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref clear. */
static INT ClearCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *theVD;
  VECTOR *v;
  INT i,l,fl,tl,n,skip,xflag;
  int j;
  double value;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"clear","no current multigrid");
    return(CMDERRORCODE);
  }
  theVD = ReadArgvVecDesc(theMG,"clear",argc,argv);
  if (theVD == NULL) {
    PrintErrorMessage('E',"clear","could not read data descriptor");
    return (PARAMERRORCODE);
  }
  if (ReadArgvOption("d",argc,argv)) {
    for (i=theMG->bottomLevel; i<=TOPLEVEL(theMG); i++)
      ClearVecskipFlags(GRID_ON_LEVEL(theMG,i),theVD);
    return (OKCODE);
  }
  if (ReadArgvOption("r",argc,argv)) {
    i = CURRENTLEVEL(theMG);
    l_dsetrandom(GRID_ON_LEVEL(theMG,i),theVD,EVERY_CLASS,1.0);
    if (ReadArgvOption("d",argc,argv))
      ClearDirichletValues(GRID_ON_LEVEL(theMG,i),theVD);
    return (OKCODE);
  }
  /* check options */
  fl = tl = CURRENTLEVEL(theMG);
  skip = FALSE;
  xflag = -1;
  value = 0.0;
  j = -1;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      fl = 0;
      break;

    case 's' :
      skip = TRUE;
      break;

    case 'x' :
      xflag = 0;
      break;

    case 'y' :
      xflag = 1;
      break;

    case 'z' :
      xflag = 2;
      break;

    case 'i' :
      if (sscanf(argv[i],"i %d",&j)!=1)
      {
        PrintErrorMessage('E',"clear","could not read value");
        return(CMDERRORCODE);
      }
      break;

    case 'v' :
      if (sscanf(argv[i],"v %lf",&value)!=1)
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
  if (j >= 0) {
    for (v = FIRSTVECTOR(GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG)));
         v != NULL; v = SUCCVC(v)) {
      n = VD_NCMPS_IN_TYPE(theVD,VTYPE(v));
      if (j < n) {
        VVALUE(v,VD_CMP_OF_TYPE(theVD,VTYPE(v),j)) = value;
        return (OKCODE);
      }
      j -= n;
    }
    return (CMDERRORCODE);
  }
  if (xflag != -1) {
    for (l=fl; l<=tl; l++)
      for (v=FIRSTVECTOR(GRID_ON_LEVEL(theMG,l)); v!=NULL; v=SUCCVC(v))
      {
        DOUBLE_VECTOR pos;

        if (VD_NCMPS_IN_TYPE(theVD,VTYPE(v)) == 0) continue;
        if (VectorPosition(v,pos)) continue;
        VVALUE(v,VD_CMP_OF_TYPE(theVD,VTYPE(v),0)) = pos[xflag];
      }
    return (OKCODE);
  }
  if (skip) {
    if (a_dsetnonskip(theMG,fl,tl,theVD,EVERY_CLASS,value)
        !=NUM_OK)
      return (CMDERRORCODE);
  }
  else {
    if (dset(theMG,fl,tl,ALL_VECTORS,theVD,value)!=NUM_OK)
      return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref makevdsub. */
static INT MakeVDsubCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *theVD,*subVD;
  VEC_TEMPLATE *vt;
  INT sub;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"makevdsub","no current multigrid");
    return(CMDERRORCODE);
  }

  theVD = ReadArgvVecDescX(theMG,"makevdsub",argc,argv,NO);

  if (theVD == NULL) {
    PrintErrorMessage('E',"makevdsub","could not read data descriptor");
    return (PARAMERRORCODE);
  }
  vt = ReadArgvVecTemplateSub(MGFORMAT(theMG),"sub",argc,argv,&sub);
  if (vt==NULL)
    REP_ERR_RETURN(PARAMERRORCODE);

  if (VDsubDescFromVT(theVD,vt,sub,&subVD))
    REP_ERR_RETURN(CMDERRORCODE);

  UserWriteF("sub descriptor '%s' for '%s' created\n",ENVITEM_NAME(subVD),ENVITEM_NAME(theVD));

  return (OKCODE);
}


/** \brief Implementation of \ref mflops. */
static INT MFLOPSCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *g;
  INT i,l;
  INT ncomp;
  VECDATA_DESC *x,*y;
  MATDATA_DESC *A;
  VECTOR *v;
  MATRIX *mat;
  INT n,m,loop;
  DOUBLE nop;
  VEC_SCALAR scal;
  DOUBLE time_ddot, time_matmul;

  /* get MG */
  theMG = GetCurrentMultigrid();
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"value","no current multigrid");
    return(CMDERRORCODE);
  }
  l = CURRENTLEVEL(theMG);
  g = GRID_ON_LEVEL(theMG,l);

  /* read arguments, allocate two vectors and a matrix */
  A = ReadArgvMatDesc(theMG,"A",argc,argv);
  x = ReadArgvVecDesc(theMG,"x",argc,argv);
  y = ReadArgvVecDesc(theMG,"y",argc,argv);
  if (x == NULL)
  {
    PrintErrorMessage('E',"x","could not read symbol");
    return (PARAMERRORCODE);
  }
  if (AllocVDFromVD(theMG,l,l,x,&y))
    return (CMDERRORCODE);
  if (AllocMDFromVD(theMG,l,l,x,x,&A))
    return (CMDERRORCODE);
  if (ReadArgvINT("loop",&loop,argc,argv))
  {
    loop = 100;
  }

  /* gather statistics */
  n = m = 0;
  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    n++;
    for (mat=VSTART(v); mat!=NULL; mat = MNEXT(mat))
      m++;
  }
  ncomp = VD_ncmps_in_otype(x,NODEVEC);
  if ((ncomp <= 0) || (VD_NCOMP(x) != ncomp)) {
    PrintErrorMessage('E',"mflops","only for NODEVEC");
    return (PARAMERRORCODE);
  }

  /* initialize */
  dset(theMG,l,l,ALL_VECTORS,x,1.0);
  dset(theMG,l,l,ALL_VECTORS,y,1.0);
  dmatset(theMG,l,l,ALL_VECTORS,A,1.0);

  /* loop */
  time_ddot = CURRENT_TIME;
  for (i=1; i<=loop; i++)
    ddot(theMG,l,l,ALL_VECTORS,x,x,scal);
  time_ddot = CURRENT_TIME - time_ddot;

  time_matmul = CURRENT_TIME;
  for (i=1; i<=loop; i++)
    dmatmul(theMG,l,l,ALL_VECTORS,y,A,x);
  time_matmul = CURRENT_TIME - time_matmul;

  if (FreeMD(theMG,l,l,A)) REP_ERR_RETURN(CMDERRORCODE);
  if (FreeVD(theMG,l,l,y)) REP_ERR_RETURN(CMDERRORCODE);

  nop = 2*n*ncomp*loop;
  UserWriteF("DDOT t=%12.4E op=%12.4E MFLOPs=%12.6f\n",
             (double)time_ddot,(double)nop,
             (double)0.000001*nop/time_ddot);
  nop = m*ncomp*ncomp*2*loop;
  UserWriteF("MMUL t=%12.4E op=%12.4E MFLOPs=%12.6f\n",
             (double)time_matmul,(double)nop,
             (double)0.000001*nop/time_matmul);

  return (OKCODE);
}


/** \brief Implementation of \ref rand. */
static INT RandCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  GRID *g;
  VECDATA_DESC *theVD;
  INT i,fl,tl,skip;
  double from_value,to_value;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"rand","no current multigrid");
    return(CMDERRORCODE);
  }

  /* check options */
  fl = tl = CURRENTLEVEL(theMG);
  skip = 0;
  from_value = 0.0;
  to_value = 1.0;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      fl = 0;
      break;

    case 's' :
      skip = 1;
      break;

    case 'f' :
      if (sscanf(argv[i],"f %lf",&from_value)!=1)
      {
        PrintErrorMessage('E',"rand","could not read from value");
        return(CMDERRORCODE);
      }
      break;

    case 't' :
      if (sscanf(argv[i],"t %lf",&to_value)!=1)
      {
        PrintErrorMessage('E',"rand","could not read to value");
        return(CMDERRORCODE);
      }
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("rand",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }

  theVD = ReadArgvVecDesc(theMG,"rand",argc,argv);

  if (theVD == NULL)
  {
    PrintErrorMessage('E',"rand","could not read data descriptor");
    return (PARAMERRORCODE);
  }

  for (i=fl; i<=tl; i++)
  {
    g = GRID_ON_LEVEL(theMG,i);
    if (l_dsetrandom2(g,theVD,EVERY_CLASS,(DOUBLE)from_value,(DOUBLE)to_value,skip)) return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref copy. */
static INT CopyCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *from,*to;
  INT fl,tl;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"copy","no current multigrid");
    return(CMDERRORCODE);
  }
  fl = tl = CURRENTLEVEL(theMG);
  if (argc<3 || argc>4)
  {
    PrintErrorMessage('E',"copy","specify exactly the f and t option");
    return(PARAMERRORCODE);
  }

  from = ReadArgvVecDescX(theMG,"f",argc,argv,NO);
  to = ReadArgvVecDesc(theMG,"t",argc,argv);

  if (from == NULL) {
    PrintErrorMessage('E',"copy","could not read 'f' symbol");
    return (PARAMERRORCODE);
  }
  if (to == NULL) {
    PrintErrorMessage('E',"copy","could not read 't' symbol");
    return (PARAMERRORCODE);
  }

  if (ReadArgvOption("a",argc,argv)) fl = 0;

  if (dcopy(theMG,fl,tl,ALL_VECTORS,to,from) != NUM_OK)
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref add. */
static INT AddCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *x,*y;
  INT fl,tl;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"copy","no current multigrid");
    return(CMDERRORCODE);
  }
  fl = tl = CURRENTLEVEL(theMG);
  if (argc<3 || argc>4)
  {
    PrintErrorMessage('E',"copy","specify exactly the f and t option");
    return(PARAMERRORCODE);
  }

  x = ReadArgvVecDesc(theMG,"x",argc,argv);
  y = ReadArgvVecDesc(theMG,"y",argc,argv);

  if (x == NULL) {
    PrintErrorMessage('E',"copy","could not read 'f' symbol");
    return (PARAMERRORCODE);
  }
  if (y == NULL) {
    PrintErrorMessage('E',"copy","could not read 't' symbol");
    return (PARAMERRORCODE);
  }
  if (ReadArgvOption("a",argc,argv)) fl = 0;

  if (dadd(theMG,fl,tl,ALL_VECTORS,x,y) != NUM_OK)
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref sub. */
static INT SubCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *x,*y;
  INT fl,tl;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"copy","no current multigrid");
    return(CMDERRORCODE);
  }
  fl = tl = CURRENTLEVEL(theMG);
  if (argc<3 || argc>4)
  {
    PrintErrorMessage('E',"copy","specify exactly the f and t option");
    return(PARAMERRORCODE);
  }

  x = ReadArgvVecDesc(theMG,"x",argc,argv);
  y = ReadArgvVecDesc(theMG,"y",argc,argv);

  if (x == NULL) {
    PrintErrorMessage('E',"copy","could not read 'f' symbol");
    return (PARAMERRORCODE);
  }
  if (y == NULL) {
    PrintErrorMessage('E',"copy","could not read 't' symbol");
    return (PARAMERRORCODE);
  }
  if (ReadArgvOption("a",argc,argv)) fl = 0;

  if (dsub(theMG,fl,tl,ALL_VECTORS,x,y) != NUM_OK)
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref homotopy. */
static INT HomotopyCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *x,*y;
  DOUBLE mu;
  DOUBLE v[MAX_VEC_COMP];
  INT i;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"homotopy","no current multigrid");
    return(CMDERRORCODE);
  }

  x = ReadArgvVecDesc(theMG,"x",argc,argv);
  if (x == NULL) {
    PrintErrorMessage('E',"homotopy","could not read 'x' symbol");
    return (PARAMERRORCODE);
  }
  y = ReadArgvVecDesc(theMG,"y",argc,argv);
  if (y == NULL) {
    PrintErrorMessage('E',"homotopy","could not read 'y' symbol");
    return (PARAMERRORCODE);
  }

  if (ReadArgvDOUBLE("v",&mu,argc,argv))
    return (PARAMERRORCODE);

  if (ReadArgvOption("a",argc,argv))
  {
    for (i=0; i<VD_NCOMP(x); i++)
      v[i] = 1.0 - mu;
    if (a_dscale(theMG,0,CURRENTLEVEL(theMG),x,EVERY_CLASS,v)!=NUM_OK)
      return (CMDERRORCODE);
    for (i=0; i<VD_NCOMP(x); i++)
      v[i] = mu;
    if (a_daxpy(theMG,0,CURRENTLEVEL(theMG),x,EVERY_CLASS,v,y)!=NUM_OK)
      return (CMDERRORCODE);
  }
  else
  {
    for (i=0; i<VD_NCOMP(x); i++)
      v[i] = 1.0 - mu;
    if (l_dscale(GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG)),
                 x,EVERY_CLASS,v)!=NUM_OK)
      return (CMDERRORCODE);
    for (i=0; i<VD_NCOMP(x); i++)
      v[i] = mu;
    if (l_daxpy(GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG)),
                x,EVERY_CLASS,v,y)!=NUM_OK)
      return (CMDERRORCODE);
  }

  return (OKCODE);
}


/** \brief Implementation of \ref interpolate. */
static INT InterpolateCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *theVD;
  INT lev,currlev;

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"interpolate","no current multigrid");
    return(CMDERRORCODE);
  }

  theVD = ReadArgvVecDescX(theMG,"interpolate",argc,argv,NO);

  if (theVD == NULL) {
    PrintErrorMessage('E',"interpolate","could not read symbol");
    return (PARAMERRORCODE);
  }

  currlev = CURRENTLEVEL(theMG);
  for (lev=1; lev<=currlev; lev++)
    if (StandardInterpolateNewVectors(GRID_ON_LEVEL(theMG,lev),theVD)!=NUM_OK)
      return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref reinit. */
static INT ReInitCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  BVP *theBVP;
  BVP_DESC theBVPDesc,*theBVPD;
  INT i,bopt;
  char BVPName[NAMESIZE];

  /* check options */
  bopt = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'b' :
      if (argv[i][1]!=' ') continue;
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
    theBVP = BVP_GetByName(BVPName);
    if(theBVP==NULL)
    {
      PrintErrorMessageF('E',"reinit","could not interpret '%s' as a BVP name",BVPName);
      return (CMDERRORCODE);
    }
    if (BVP_SetBVPDesc(theBVP,&theBVPDesc)) return (CMDERRORCODE);
    theBVPD = &theBVPDesc;
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
    theBVPD = MG_BVPD(theMG);
  }

  if (BVPD_CONFIG(theBVPD)!=NULL)
    if ((*BVPD_CONFIG (theBVPD))(argc,argv)!=0)
      return (CMDERRORCODE);

  return(OKCODE);
}

/* see formats.c for the man page */

/** \brief Implementation of \ref newformat. */
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


/** \brief Implementation of \ref delformat. */
static INT DeleteFormatCommand (INT argc, char **argv)
{
  INT err;
  char fmtname[NAMESIZE];

  NO_OPTION_CHECK(argc,argv);

  if (sscanf(argv[0],"delformat %s",fmtname)!=1)
  {
    PrintErrorMessage('E',"delformat","specify format to delete");
    return (PARAMERRORCODE);
  }

  err = RemoveFormatWithSubs(fmtname);

  switch (err)
  {
  case GM_OK : return (OKCODE);
  default : return (CMDERRORCODE);
  }
}


/** \brief Implementation of \ref setpf. */
static INT SetPrintingFormatCommand (INT argc, char **argv)
{
  INT err;
  MULTIGRID *theMG;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"setpf","there is no current multigrid\n");
    return (CMDERRORCODE);
  }
  err = SetPrintingFormatCmd(theMG,argc,argv);

  switch (err)
  {
  case 0 : return (OKCODE);
  case 1 : PrintHelp("setpf",HELPITEM,NULL);
    return (PARAMERRORCODE);
  default : return (CMDERRORCODE);
  }
}


/** \brief Implementation of \ref showpf. */
static INT ShowPrintingFormatCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  DisplayPrintingFormat();

  return (OKCODE);
}

/****************************************************************************/
/** \brief Return a pointer to the current numproc

   This function returns a pointer to the current numproc.

   RETURN VALUE:
   .n      pointer to numproc
   .n      NULL if there is no current numproc.
 */
/****************************************************************************/

static NP_BASE *GetCurrentNumProc (void)
{
  return (currNumProc);
}

/****************************************************************************/
/** \brief Set the current NumProc if it is valid
 *
 * @param theNumProc - pointer to multigrid

   This function sets the current NumProc if it is valid, i. e.
   the function checks whether 'theNumProc' acually points to a numproc.
   It can be NULL only if no numproc is defined.

 * @return <ul>
   .n    0 if ok
   .n    1 if theNumProc is not in the numproc list
 */
/****************************************************************************/

static INT SetCurrentNumProc (NP_BASE *theNumProc)
{
  currNumProc = theNumProc;
  return (0);
}


/** \brief Implementation of \ref npexecute. */
static INT ExecuteNumProcCommand (INT argc, char **argv)
{
  char theNumProcName[NAMESIZE];
  NP_BASE *theNumProc;
  MULTIGRID *theMG;
  INT err;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"npexecute","there is no current multigrid\n");
    return (CMDERRORCODE);
  }

  /* get NumProc */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" npexecute %",NAMELENSTR,"[ -~]")),theNumProcName)!=1) || (strlen(theNumProcName)==0))
  {
    theNumProc = GetCurrentNumProc();
    if (theNumProc == NULL)
    {
      PrintErrorMessage('E',"npexecute",
                        "there is no current numerical procedure");
      return (CMDERRORCODE);
    }
  }
  else
  {
    theNumProc = GetNumProcByName (theMG,theNumProcName,"");
    if (theNumProc == NULL)
    {
      PrintErrorMessage('E',"npexecute",
                        "cannot find specified numerical procedure");
      return (CMDERRORCODE);
    }
  }
  if (theNumProc->status != NP_EXECUTABLE) {
    PrintErrorMessage('E',"npexecute",
                      "the num proc is not executable");
    return (CMDERRORCODE);
  }
  if ((err=((*theNumProc->Execute)(theNumProc,argc,argv)))!=0)
  {
    PrintErrorMessageF('E',"npexecute","execution of '%s' failed (error code %d)",theNumProcName,err);
    return (CMDERRORCODE);
  }

  return(OKCODE);
}


/** \brief Implementation of \ref npdisplay. */
static INT NumProcDisplayCommand (INT argc, char **argv)
{
  NP_BASE *theNumProc;
  MULTIGRID *theMG;
  INT i,All,Class,err;
  char theNumProcName[NAMESIZE],ClassName[NAMESIZE];

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"npdisplay","there is no current multigrid\n");
    return (CMDERRORCODE);
  }

  /* check options */
  All = Class = FALSE;
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'a' :
      All = TRUE;
      break;

    case 'c' :
      if (sscanf(argv[i],expandfmt(CONCAT3("c %",NAMELENSTR,"[ -~]")),ClassName)!=1)
      {
        PrintErrorMessage('W',"npdisplay","no class specified\n");
        UserWrite("enroled classes are:\n");
        if (MGListNPClasses(theMG))
          return (CMDERRORCODE);
        return (OKCODE);
      }
      Class = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintHelp("npdisplay",HELPITEM,buffer);
      return (PARAMERRORCODE);
    }
  if (All && Class)
  {
    PrintErrorMessage('E',"npdisplay","a and c option are mutually exclusive");
    return (CMDERRORCODE);
  }
  if (Class)
  {
    if (MGListNPsOfClass(theMG,ClassName))
      return (CMDERRORCODE);
    return (OKCODE);
  }
  if (All)
  {
    if (MGListAllNPs(theMG))
      return (CMDERRORCODE);
    return (OKCODE);
  }

  /* get NumProc */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" npdisplay %",NAMELENSTR,"[ -~]")),theNumProcName)!=1) || (strlen(theNumProcName)==0))
  {
    theNumProc = GetCurrentNumProc();
    if (theNumProc == NULL)
    {
      PrintErrorMessage('E',"npdisplay","there is no current numerical procedure");
      return (CMDERRORCODE);
    }
  }
  else
  {
    theNumProc = GetNumProcByName (theMG,theNumProcName,"");
    if (theNumProc == NULL)
    {
      PrintErrorMessage('E',"npdisplay","cannot find specified numerical procedure");
      return (CMDERRORCODE);
    }
  }
  if ((err=ListNumProc(theNumProc))!=0)
  {
    PrintErrorMessageF('E',"npdisplay","execution of '%s' failed (error code %d)",theNumProcName,err);
    return (CMDERRORCODE);
  }

  return(OKCODE);
}


/** \brief Implementation of \ref npcreate. */
static INT NumProcCreateCommand (INT argc, char **argv)
{
  char theNumProcName[NAMESIZE];
  char ConstructName[NAMESIZE];
  NP_BASE *theNumProc;
  MULTIGRID *theMG;
  INT err,ignore;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"npexecute","there is no current multigrid\n");
    return (CMDERRORCODE);
  }

  /* get NumProc name */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" npcreate %",NAMELENSTR,"[ -~]")),theNumProcName)!=1) || (strlen(theNumProcName)==0))
  {
    PrintErrorMessage('E',"npcreate","specify the name of the theNumProcName to create");
    return (PARAMERRORCODE);
  }
  if (ReadArgvChar("c",ConstructName,argc,argv))
  {
    PrintErrorMessage('E',"npcreate","specify the name of the constructor");
    return (PARAMERRORCODE);
  }
  ignore=ReadArgvOption("i",argc,argv);
  if (ignore)
  {
    theNumProc = GetNumProcByName(theMG,theNumProcName,"");
  }
  if (!ignore || theNumProc==NULL)
    if ((err=CreateObject(theMG,theNumProcName,ConstructName))!=0)
    {
      UserWriteF("creating of '%s' failed (error code %d)\n",theNumProcName,err);
      return (CMDERRORCODE);
    }
  theNumProc = GetNumProcByName(theMG,theNumProcName,"");
  if (SetCurrentNumProc(theNumProc)) return(PARAMERRORCODE);

  return(OKCODE);
}


/** \brief Implementation of \ref npinit. */
static INT NumProcInitCommand (INT argc, char **argv)
{
  char theNumProcName[NAMESIZE];
  NP_BASE *theNumProc;
  MULTIGRID *theMG;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"npinit","there is no current multigrid\n");
    return (CMDERRORCODE);
  }

  /* get NumProc */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" npinit %",NAMELENSTR,"[ -~]")),theNumProcName)!=1) || (strlen(theNumProcName)==0))
  {
    theNumProc = GetCurrentNumProc();
    if (theNumProc == NULL)
    {
      PrintErrorMessage('E',"npinit","there is no current numerical procedure");
      return (CMDERRORCODE);
    }
  }
  else
  {
    theNumProc = GetNumProcByName (theMG,theNumProcName,"");
    if (theNumProc == NULL)
    {
      sprintf(buffer,"cannot find specified numerical procedure '%s'",theNumProcName);
      PrintErrorMessage('E',"npinit",buffer);
      return (CMDERRORCODE);
    }
  }
  theNumProc->status = (*theNumProc->Init)(theNumProc,argc,argv);
  switch (theNumProc->status) {
  case NP_NOT_INIT :
    UserWriteF("num proc %s has status NOT_INIT\n",theNumProcName);
    return (CMDERRORCODE);
  case NP_NOT_ACTIVE :
    UserWriteF("num proc %s has status NOT_ACTIVE\n",theNumProcName);
    return (CMDERRORCODE);
  case NP_ACTIVE :
    UserWriteF("num proc %s has status ACTIVE\n",theNumProcName);
    break;
  case NP_EXECUTABLE :
    UserWriteF("num proc %s has status EXECUTABLE\n",theNumProcName);
    break;
  default :
    PrintErrorMessage('E',"npinit","unknown status");
    return (CMDERRORCODE);
  }

  return(OKCODE);
}


/** \brief Implementation of \ref scnp. */
static INT SetCurrentNumProcCommand (INT argc, char **argv)
{
  NP_BASE *theNumProc;
  MULTIGRID *theMG;
  char npname[NAMESIZE];

  NO_OPTION_CHECK(argc,argv);

  /* get solver name */
  if (sscanf(argv[0],expandfmt(CONCAT3(" scnp %",NAMELENSTR,"[ -~]")),npname)!=1)
  {
    PrintHelp("scnp",HELPITEM," (specify current NumProc name)");
    return(PARAMERRORCODE);
  }
  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"npexecute","there is no current multigrid\n");
    return (CMDERRORCODE);
  }
  theNumProc = GetNumProcByName (theMG,npname,"");
  if (theNumProc == NULL)
  {
    PrintErrorMessage('E',"npexecute","cannot find specified numerical procedure");
    return (CMDERRORCODE);
  }

  if (SetCurrentNumProc(theNumProc))
    return (CMDERRORCODE);

  return(OKCODE);
}


/** \brief Implementation of \ref createvector. */
static INT CreateVecDescCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char name[NAMESIZE];

  if (ReadArgvChar("m",name,argc,argv))
    theMG = currMG;
  else
    theMG = GetMultigrid(name);
  if (theMG==NULL) {
    PrintErrorMessage('E',"createvector","no current multigrid");
    return(CMDERRORCODE);
  }
  if (CreateVecDescCmd(theMG,argc,argv))
    return(CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref creatematrix. */
static INT CreateMatDescCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char name[NAMESIZE];

  if (ReadArgvChar("m",name,argc,argv))
    theMG = currMG;
  else
    theMG = GetMultigrid(name);
  if (theMG==NULL) {
    PrintErrorMessage('E',"creatematrix","no current multigrid");
    return(CMDERRORCODE);
  }
  if (CreateMatDescCmd(theMG,argc,argv))
    return(CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref freematrix. */
static INT FreeMatDescCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  char name[NAMESIZE];

  if (ReadArgvChar("m",name,argc,argv))
    theMG = currMG;
  else
    theMG = GetMultigrid(name);
  if (theMG==NULL) {
    PrintErrorMessage('E',"freematrix","no current multigrid");
    return(CMDERRORCODE);
  }
  if (FreeMatDescCmd(theMG,argc,argv))
    return(CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref symlist. */
static INT SymListCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  VECDATA_DESC *vd;
  MATDATA_DESC *md;
  char name[NAMESIZE];
  INT i,res,modifiers;

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"npinit","there is no current multigrid\n");
    return (CMDERRORCODE);
  }

  /* first get modifiers */
  modifiers = 0;
  if (ReadArgvOption("scal",argc,argv))
    SET_FLAG(modifiers,SCAL_PROP);
  if (ReadArgvOption("alloc",argc,argv))
    SET_FLAG(modifiers,ALLOC_STAT);

  /* scan V and M options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'V' :
      res = sscanf(argv[1],"V %s",name);
      if (res!=1)
      {
        /* print all vectors */
        for (vd = GetFirstVector(theMG); vd != NULL; vd = GetNextVector(vd))
        {
          DisplayVecDataDesc(vd,modifiers,buffer);
          UserWrite(buffer);
        }
        return (OKCODE);
      }
      vd = GetVecDataDescByName(theMG,name);
      if (vd!=NULL)
      {
        DisplayVecDataDesc(vd,modifiers,buffer);
        UserWrite(buffer);
      }
      break;

    case 'M' :
      res = sscanf(argv[1],"M %s",name);
      if (res!=1)
      {
        /* print all matrices */
        for (md = GetFirstMatrix(theMG); md != NULL; md = GetNextMatrix(md))
        {
          DisplayMatDataDesc(md,buffer);
          UserWrite(buffer);
        }
        return (OKCODE);
      }
      md = GetMatDataDescByName(theMG,name);
      if (md!=NULL)
      {
        DisplayMatDataDesc(md,buffer);
        UserWrite(buffer);
        return (OKCODE);
      }
      break;
    }

  return (OKCODE);
}


/** \brief Implementation of \ref setkey. */
static INT SetCommandKeyCommand (INT argc, char **argv)
{
  INT opt,begin,i,j,odd,ShowBar;
  char *ChatPtr;
  char CmdBuffer[INPUTBUFFERLEN],comment[KEY_COMMENT_SIZE];

  /* check input */
  if (argc < 3 )
    return(CMDERRORCODE);

  if (strlen(argv[1])!=1)
  {
    PrintErrorMessage('E',"setkey","only one character for cmd key");
    return (PARAMERRORCODE);
  }

  begin = 2;

  /* comment given? */
  comment[0] = '\0';
  if (argv[begin][0]=='c')
  {
    if (sscanf(argv[begin],expandfmt(CONCAT3("c %",KEY_COMMENT_LEN_STR,"[ -~]")),comment)!=1)
    {
      PrintErrorMessage('E',"setkey","could not read comment");
      return (PARAMERRORCODE);
    }
    begin++;
  }

  /* show bar before in keylist? */
  ShowBar = FALSE;
  if (argv[begin][0]=='-')
  {
    ShowBar = TRUE;
    begin++;
  }

  /* store input */
  ChatPtr = CmdBuffer;
  for (opt=begin; opt<argc; opt++)
  {
    *ChatPtr = '$';
    ChatPtr++;
    strcpy(ChatPtr,argv[opt]);
    ChatPtr += strlen(argv[opt]);
  }

  /* check input */
  if ((argv[begin][0] != '\"') || (argv[argc-1][strlen(argv[argc-1])-1] != '\"'))
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
  if (SetCmdKey(argv[1][0],comment,ShowBar,CmdBuffer)!=0)
  {
    PrintErrorMessage('E',"setkey","cannot create cmd key");
    return(CMDERRORCODE);
  }

  return(OKCODE);
}


/** \brief Implementation of \ref delkey. */
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


/** \brief Implementation of \ref keylist. */
static INT ListCommandKeysCommand (INT argc, char **argv)
{
  if (argc>2)
  {
    PrintErrorMessage('E',"setkey","max of one option exceeded");
    return (PARAMERRORCODE);
  }
  if (argc==2)
    if (argv[1][0]=='l')
    {
      ListCmdKeys(YES);
      return (OKCODE);
    }

  ListCmdKeys(NO);

  return (OKCODE);
}


/** \brief Implementation of \ref refreshon. */
static INT RefreshOnCommand (INT argc, char **argv)
{
  DOUBLE factor = 1.0;

  if (argc >= 2 && argv[1][0] == 'b') {
    sscanf(argv[1],"b %lf", &factor);
    SetRefreshState(ON, YES, factor);
  }
  else
    SetRefreshState(ON, NO, factor);

  return(OKCODE);
}


/** \brief Implementation of \ref refreshoff. */
static INT RefreshOffCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  SetRefreshState(OFF, NO, 1.0);
  return(OKCODE);
}

/****************************************************************************/
/** \brief Implementation of \ref ???

   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   This function tests interesting machine parameters.

   RETURN VALUE:
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT MachineTestCommand (INT argc, char **argv)
{
  /* determine parameters of memory */
        #ifdef ModelP
  if (me == master)
        #endif
  MemoryParameters();

        #ifdef ModelP
  /* determine performance of network */
  /*	NetworkPerformance(); */
        #endif

  return (OKCODE);
}

/****************************************************************************/
/** \brief Implementation of \ref system

   SystemCommand - calls system routine

   .  argc - number of arguments (incl. its own name)
   .  argv - array of strings giving the arguments

   This function calls a system routine.

   'system \<system~command>'

   .   \<system~command>	- this command string is passed to the system

   KEYWORDS:
   system, cshell, UNIX
 */
/****************************************************************************/

static INT SystemCommand (INT argc, char **argv)
{
  char *p;

  if (strlen(argv[0])<8) return (PARAMERRORCODE);
  p = argv[0]+7;

  printf("system \n%s\n",p);

  if (system(p)==-1)
    UserWrite("system-error\n");

  return (OKCODE);
}


/** \brief Implementation of \ref resetCEstat. */
static INT ResetCEstatCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  ResetCEstatistics();

  return (OKCODE);
}


/** \brief Implementation of \ref printCEstat. */
static INT PrintCEstatCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  PrintCEstatistics();

  return (OKCODE);
}


/** \brief Implementation of \ref heapstat. */
static INT HeapStatCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: HeapStatCommand(): me not in Context,"\
                     " no heap stat\n",me))
    return(OKCODE);
  }
        #endif

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    UserWrite("no multigrid open\n");
    return (OKCODE);
  }

  HeapStat((const HEAP *)MGHEAP(theMG));

  return (OKCODE);
}



/** \brief Implementation of \ref getheapused. */
static INT GetHeapUsedCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  INT used;

        #ifdef ModelP
  if (!CONTEXT(me)) {
    PRINTDEBUG(ui,0,("%2d: GetHeapUsedCommand(): me not in Context,"\
                     " no heap info\n",me))
    return(OKCODE);
  }
        #endif

  NO_OPTION_CHECK(argc,argv);

  theMG = currMG;
  if (theMG==NULL)
  {
    UserWrite("no multigrid open\n");
    return (OKCODE);
  }

  used = (INT)HeapUsed(MGHEAP(theMG));

        #ifdef ModelP
  used = UG_GlobalMaxINT(used);
        #endif

  if (SetStringValue(":HEAPUSED",used)!=0) {
    PrintErrorMessage('E',"getheapused","could not get string variable :HEAPUSED");
    return (CMDERRORCODE);
  }

  return (OKCODE);
}
/****************************************************************************/
/** \brief Create struct where findrange stores results (min and max)

   This function creates the struct ':findrange'.

   RETURN VALUE:
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InitFindRange (void)
{
  /* install a struct for findrange */
  if (MakeStruct(":findrange")!=0)
    return (1);
  return (0);
}

/****************************************************************************/
/** \brief Create struct :screensize

   This function creates the struct _screensize for the 'screensize' command.

   RETURN VALUE:
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT InitScreenSize (void)
{
  /* install a struct for screensize */
  if (MakeStruct(":screensize")!=0)
    return (1);
  return (0);
}


/** \brief Implementation of \ref lb. */
INT NS_DIM_PREFIX LBCommand (INT argc, char **argv)
{
                #ifndef ModelP
  /* dummy command in seriell version */
  return(OKCODE);
                #endif

                #ifdef ModelP
  INT res,cmd_error,maxlevel,i;
  int minlevel,cluster_depth,threshold,Const,n,c,
      strategy,eigen,loc,dims,weights,coarse,mode,iopt;
  char levelarg[32];
  MULTIGRID *theMG;

  theMG = currMG;

  if (theMG == NULL)
  {
    UserWrite("LBCommand: no open multigrid\n");
    return(OKCODE);
  }

  if (procs==1) return(OKCODE);

  /* defaults */
  minlevel                = 1;
  cluster_depth   = 2;
  maxlevel                = MAXLEVEL;
  threshold               = 1;
  Const                   = 1;
  n                               = 500;
  c                               = 10;
  strategy                = 0;
  eigen                   = 0;
  loc                             = 0;
  dims                    = 1;
  weights                 = 0;
  coarse                  = 0;
  mode                    = 0;
  iopt                    = 10000;

  /* scan lb strategy*/
  res = sscanf(argv[0]," lb %d", &strategy);
  if (res > 1)
  {
    UserWriteF("lb [<strategy>] [$c <minlevel>] [$d <depth>] [$f <maxlevel>] [$e <minelem>] \n");
    UserWriteF("default lb 0 $c 1 $d 2 $e 1\n");
    return(OKCODE);
  }

  /* parse options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'c' :
      sscanf(argv[i],"c %d",&minlevel);
      break;

    case 'd' :
      UserWriteF("lb: depth parameter skipped\n");
      break;

    case 'f' :
      sscanf(argv[i],"f %d",&maxlevel);
      break;

    case 'e' :
      UserWriteF("lb: minelem parameter skipped\n");
      break;

    default :
      UserWriteF("lb [<strategy>] [$c <minlevel>] [$d <depth>] [$f <maxlevel>] [$e <minelem>]\n");
      UserWriteF("default lb 0 $c 1 $d 2\n");
      break;
    }

  /* check for parameter consistency! */
  cmd_error = 0;

  if ((minlevel<0)||(minlevel>TOPLEVEL(theMG)))
  {
    UserWriteF("Choose <minlevel>: 0-%d (toplevel)\n",TOPLEVEL(theMG));
    cmd_error = 1;
  }

  if (maxlevel < minlevel)
  {
    UserWriteF("Choose <maxlevel>: %d-%d (MAXLEVEL)\n",minlevel,MAXLEVEL);
    cmd_error = 1;
  }
  if (maxlevel < TOPLEVEL(theMG))
  {
    UserWriteF("%s: maxlevel reached: no loadbalancing done!\n"
               "    maxlevel=%d toplevel=%d\n",argv[0],maxlevel,TOPLEVEL(theMG));
    return(OKCODE);
  }

  if (Const < 1)
  {
    UserWriteF("Choose <minelem> > 0\n");
    cmd_error = 1;
  }

  if (strategy==1 || strategy==2)
  {
    if (strategy == 1)
    {
      coarse = 50;
      loc = 1;
    }
    eigen = 1;
  }

  if (strategy>2 || strategy<6)
  {
    weights = 1;
  }

  if (cmd_error) return(CMDERRORCODE);

  sprintf(levelarg,"%d",minlevel);
  lbs(levelarg, theMG);

  return(OKCODE);
                #endif
}


#ifdef ModelP
/** \brief Implementation of \ref lbs. */
static INT LBSCommand (INT argc, char **argv)
{
  MULTIGRID *theCurrMG;

  theCurrMG = currMG;
  if (theCurrMG==NULL)
  {
    PrintErrorMessage('W',"mglist","no multigrid open\n");
    return (OKCODE);
  }

  if (argc==2)
    lbs(argv[1], theCurrMG);
  else
    lbs("0", theCurrMG);

  return(OKCODE);
}


/** \brief Implementation of \ref context. */
static INT ContextCommand (INT argc, char **argv)
{
  INT proc = INT_MAX;
  INT flag_all, flag_empty, flag_invert;


  flag_all   = ReadArgvOption("a", argc, argv);
  flag_empty = ReadArgvOption("e", argc, argv);
  flag_invert= ReadArgvOption("i", argc, argv);

  ReadArgvINT("context", &proc, argc, argv);
  if (proc<0 || proc>=procs)
  {
    if (proc!=INT_MAX && me==0)
      UserWriteF("context: invalid processor id (procs=%d)\n", procs);
  }
  else
  {
    /* switch context for proc on/off */
    CONTEXT(proc) = 1-CONTEXT(proc);
  }

  if (proc==INT_MAX)
  {
    if (flag_all && !flag_empty)
    {
      int p; for(p=0; p<procs; p++) CONTEXT(p) = 1;
    }
    if (flag_empty && !flag_all)
    {
      int p; for(p=0; p<procs; p++) CONTEXT(p) = 0;
    }
    if (flag_empty && flag_all)
    {
      if (me==0)
        UserWriteF("context: invalid option combination\n");
    }

    if (flag_invert)
    {
      int p; for(p=0; p<procs; p++) CONTEXT(p) = 1-CONTEXT(p);
    }
  }

  ddd_DisplayContext();
  return(OKCODE);
}


/** \brief Implementation of \ref pstat. */
static INT PStatCommand (INT argc, char **argv)
{
  if (argc!=2)
    return (CMDERRORCODE);

  ddd_pstat(argv[1]);
  return(OKCODE);
}

#ifdef USE_FAMG
static INT pamgCheckCommand (INT argc, char **argv)
{
  MULTIGRID *theMG;
  int level;

  theMG = currMG;
  level = 0;
  if (argv[1][0]=='l')
  {
    if (sscanf(argv[1],"l %d",&level)!=1)
      return (CMDERRORCODE);
  }

  if (pamgCheckDo( theMG, level ) > 0)
    return(CMDERRORCODE);

  return (OKCODE);
}
#endif

#endif /* ModelP */


#ifdef Debug
/** \brief Implementation of \ref debug. */
static INT DebugCommand (INT argc, char **argv)
{
  char *module;
  int l;

  if (argc<2 || argc>3)
  {
    UserWriteF("usage: debug $<module> [$<level>]\n");
    return (CMDERRORCODE);
  }

  if (argc==3)
  {
    if              (strcmp("init",argv[1])==0) Debuginit               = atoi(argv[2]);
    else if (strcmp("dddif",argv[1])==0) Debugdddif              = atoi(argv[2]);
    else if (strcmp("dev",argv[1])==0) Debugdev                = atoi(argv[2]);
    else if (strcmp("dom",argv[1])==0) Debugdom                = atoi(argv[2]);
    else if (strcmp("gm",argv[1])==0) Debuggm                 = atoi(argv[2]);
    else if (strcmp("graph",argv[1])==0) Debuggraph              = atoi(argv[2]);
    else if (strcmp("low",argv[1])==0) Debuglow                = atoi(argv[2]);
    else if (strcmp("machines",argv[1])==0) Debugmachines   = atoi(argv[2]);
    else if (strcmp("np",argv[1])==0) Debugnp             = atoi(argv[2]);
    else if (strcmp("ui",argv[1])==0) Debugui                 = atoi(argv[2]);
    else if (strcmp("time",argv[1])==0) Debugtime               = atoi(argv[2]);
    else if (strcmp("pclib",argv[1])==0) Debugpclib              = atoi(argv[2]);
    else if (strcmp("appl",argv[1])==0) Debugappl               = atoi(argv[2]);
    else
    {
      UserWriteF("no debug variable for module %s found!\n",argv[1]);
      return (CMDERRORCODE);
    }
    UserWriteF("set debuglevel for module %s to %d\n",argv[1],atoi(argv[2]));
  }
  else
  {
    if              (strcmp("init",argv[1])==0)             {module="init";         l=Debuginit;}
    else if (strcmp("dddif",argv[1])==0)    {module="dddif";        l=Debugdddif;}
    else if (strcmp("dev",argv[1])==0)              {module="dev";          l=Debugdev;}
    else if (strcmp("dom",argv[1])==0)              {module="dom";          l=Debugdom;}
    else if (strcmp("gm",argv[1])==0)               {module="gm";           l=Debuggm;}
    else if (strcmp("graph",argv[1])==0)    {module="graph";        l=Debuggraph;}
    else if (strcmp("low",argv[1])==0)              {module="low";          l=Debuglow;}
    else if (strcmp("machines",argv[1])==0) {module="machines";     l=Debugmachines;}
    else if (strcmp("np",argv[1])==0)           {module="np";           l=Debugnp;}
    else if (strcmp("ui",argv[1])==0)               {module="ui";           l=Debugui;}
    else if (strcmp("time",argv[1])==0)             {module="time";         l=Debugtime;}
    else if (strcmp("pclib",argv[1])==0)    {module="pclib";        l=Debugpclib;}
    else if (strcmp("appl",argv[1])==0)             {module="appl";         l=Debugappl;}
    else
    {
      UserWriteF("no debug variable for module %s found!\n",argv[1]);
      return (CMDERRORCODE);
    }

    UserWriteF("debuglevel for module %s is %d\n",module,l);
  }
  return(OKCODE);
}
#endif


#ifdef Debug
/** \brief Implementation of \ref trace. */
static INT TraceCommand (INT argc, char **argv)
{
  int i;

  for (i=1; i<argc; i++)
    STR_SWITCH(argv[i])
    STR_CASE("blas")
    int n;
  if (sscanf(argv[i],"blas %d",&n)==1)
    TraceUGBlas(n);
  else
    TraceUGBlas(TRBL_PARAMS);
  STR_BREAK

  STR_DEFAULT
  sprintf(buffer,"(invalid option '%s')",argv[i]);
  PrintHelp("trace",HELPITEM,buffer);
  return (PARAMERRORCODE);
  STR_BREAK

    STR_SWITCH_END

  return (OKCODE);
}
#endif


#ifdef Debug
/** \brief Implementation of \ref reperr. */
static INT RepErrCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  PrintRepErrStack(UserWriteF);

  return (OKCODE);
}
#endif



#ifdef Debug
/** \brief Implementation of \ref timing. */
static INT TimingCommand (INT argc, char **argv)
{
  INT i;

  if (ReadArgvOption("r",argc,argv)) {
    DEBUG_TIME_RESET;
    return (OKCODE);
  }
  if (__debug_time_count==0)
    UserWrite("no timing\n");
  else
  {
    UserWrite("timing:\n\n");
    for (i=0; i<__debug_time_count; i++) {
      UserWriteF("%2d: File:%15s, Line:%5d elapsed time%10.4f",
                 i,__debug_time_file[i],__debug_time_line[i],
                 __debug_time[i]-__debug_time[0]);
      if (i > 0) UserWriteF(" diff%8.4f",
                            __debug_time[i]-__debug_time[i-1]);
      UserWriteF("\n");
    }
  }
  return (OKCODE);
}
#endif


/** \brief Implementation of \ref showconfig. */
static INT ShowConfigCommand (INT argc, char **argv)
{
  NO_OPTION_CHECK(argc,argv);

  UserWrite("Configuration of this programm:\n");

#ifdef _2
  UserWrite("    Dimension:    2\n");
#elif defined _3
  UserWrite("    Dimension:    3\n");
#else
  UserWrite("    Dimension:    unknown\n");
#endif

#ifdef ModelP
  UserWrite("    Model:        parallel\n");
#else
  UserWrite("    Model:        sequential\n");
#endif

  UserWriteF("   Architecture: %s\n",ARCHNAME);

#ifdef Debug
  UserWrite("    Debugging:    ON\n");
#elif defined NDEBUG
  UserWrite("    Debugging:    OFF\n");
#else
  UserWrite("    Debugging:    unknown\n");
#endif

#ifdef RIF_SOCKETS
  UserWrite("    remote:       ON\n");
#else
  UserWrite("    remote:       OFF\n");
#endif

  return (OKCODE);
}



static INT ClearArray (ARRAY *theAR)
{
  INT i, size;

  size = 1;
  for (i=0; i<AR_NVAR(theAR); i++)
    size *= AR_VARDIM(theAR,i);
  for (i=0; i<size; i++)
    AR_DATA(theAR,i) = 0.0;

  return (0);
}

/****************************************************************************/
/** \brief Allocate a new array structure

   .  name - name under which the array is allocated in '/Array'
   .  nVar - number of dimensions of the data field
   .  VarDim - extension of the data field in each dimension

   Allocate a new array structure in the directory '/Array' and
   allocate the data field. The maximum number of dimensions is
   'AR_NVAR_MAX'.

   RETURN VALUE:
   .n     pointer to allocated array
   .n     NULL on error.
 */
/****************************************************************************/

static ARRAY *CreateArray (char *name, INT nVar, INT *VarDim)
{
  INT i, size;
  ARRAY *theAR;

  if (nVar<1 || nVar>AR_NVAR_MAX) return (NULL);

  /* change to directory */
  if (ChangeEnvDir("/Array")==NULL)
    return(NULL);

  /* get size */
  size = sizeof(DOUBLE);
  for (i=0; i<nVar; i++)
    size *= VarDim[i];
  size += sizeof(ARRAY) - sizeof(DOUBLE);

  /* allocate structure */
  theAR = (ARRAY*)MakeEnvItem (name,theArrayVarID,size);
  if (theAR==NULL) return (NULL);

  /* init structure */
  ENVITEM_LOCKED(theAR) = 0;
  AR_NVAR(theAR) = nVar;
  for (i=0; i<nVar; i++)
    AR_VARDIM(theAR,i) = VarDim[i];

  if (ClearArray(theAR)) return (NULL);

  return (theAR);
}

/****************************************************************************/
/** \brief Set one single entry of the data field of the array

   .  theAR - array structure to work on
   .  Point - specify the coordinate of the entry in each dimension
   .  value - value to be stored

   Set one single entry of the data field of the array to the given value.

   RETURN VALUE:
   .n    0 always
 */
/****************************************************************************/

static INT WriteArray (ARRAY *theAR, INT *Point, DOUBLE value)
{
  INT i, pos;

  pos = Point[AR_NVAR(theAR)-1];
  for (i=AR_NVAR(theAR)-2; i>=0; i--)
    pos = Point[i] + AR_VARDIM(theAR,i)*pos;
  AR_DATA(theAR,pos) = value;

  return (0);
}

/****************************************************************************/
/** \brief Read one single entry of the data field of the array

   .  theAR - array structure to work on
   .  Point - specify the coordinate of the entry in each dimension
   .  value - read value

   Read one single entry of the data field of the array.

   RETURN VALUE:
   .n    0 always
 */
/****************************************************************************/

static INT ReadArray (ARRAY *theAR, INT *Point, DOUBLE *value)
{
  INT i, pos;

  pos = Point[AR_NVAR(theAR)-1];
  for (i=AR_NVAR(theAR)-2; i>=0; i--)
    pos = Point[i] + AR_VARDIM(theAR,i)*pos;
  value[0] = AR_DATA(theAR,pos);

  return (0);
}


/** \brief Implementation of \ref crar. */
static INT CreateArrayCommand (INT argc, char **argv)
{
  INT i, nVar, VarDim[AR_NVAR_MAX];
  int iValue;
  char name[128];

  nVar = argc-2;
  if (nVar<1 || nVar>AR_NVAR_MAX)
    return (CMDERRORCODE);
  if (argv[1][0]=='n')
  {
    if (sscanf(argv[1],"n %s",name)!=1)
      return (CMDERRORCODE);
  }
  for (i=0; i<nVar; i++)
  {
    if (sscanf(argv[i+2],"%d",&iValue)!=1)
      return (CMDERRORCODE);
    if (iValue<1)
      return (CMDERRORCODE);
    VarDim[i] = iValue;
  }

  /* create Array */
  if (CreateArray(name,nVar,VarDim)==NULL)
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref dear. */
static INT DeleteArrayCommand (INT argc, char **argv)
{
  char name[128];
  ARRAY *theAR;

  if (argv[1][0]=='n')
  {
    if (sscanf(argv[1],"n %s",name)!=1)
      return (CMDERRORCODE);
  }

  /* find array */
  if (ChangeEnvDir("/Array")==NULL)
  {
    PrintErrorMessage('F',"DeleteArrayCommand","could not changedir to /Array");
    return(CMDERRORCODE);
  }
  theAR = (ARRAY *)SearchEnv(name,".",theArrayVarID,SEARCHALL);
  if (theAR==NULL)
    return (CMDERRORCODE);

  /* delete Array */
  if (RemoveEnvItem((ENVITEM *)theAR))
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref saar. */
static INT SaveArrayCommand (INT argc, char **argv)
{
  INT i, size;
  char name[128];
  ARRAY *theAR;
  FILE *stream;

  if (argv[1][0]=='n')
  {
    if (sscanf(argv[1],"n %s",name)!=1)
      return (CMDERRORCODE);
  }

  /* find array */
  if (ChangeEnvDir("/Array")==NULL)
  {
    PrintErrorMessage('F',"SaveArrayCommand","could not changedir to /Array");
    return(CMDERRORCODE);
  }
  theAR = (ARRAY *)SearchEnv(name,".",theArrayVarID,SEARCHALL);
  if (theAR==NULL)
    return (CMDERRORCODE);

  /* save Array */
  strcat(name,".array");
  if (arraypathes_set)
    stream = FileOpenUsingSearchPaths(name,"w","arraypathes");
  else
    stream = fileopen(name,"w");
  if (stream==NULL)
  {
    PrintErrorMessage('E',"SaveArrayCommand","cannot open file");
    return(CMDERRORCODE);
  }

  /* store */
  if (fwrite((void*)(&(theAR->nVar)),sizeof(INT),1,stream)!=1) return(CMDERRORCODE);
  if (fwrite((void*)theAR->VarDim,sizeof(INT),AR_NVAR(theAR),stream)!=AR_NVAR(theAR)) return(CMDERRORCODE);
  size = 1;
  for (i=0; i<AR_NVAR(theAR); i++)
    size *= AR_VARDIM(theAR,i);
  if (fwrite((void*)(theAR->data),sizeof(DOUBLE),size,stream)!=size) return(CMDERRORCODE);
  if (fclose(stream)) return(CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref loar. */
static INT LoadArrayCommand (INT argc, char **argv)
{
  INT i, size, nVar, VarDim[AR_NVAR_MAX];
  char name[128], filename[128];
  ARRAY *theAR;
  FILE *stream;

  if (argv[1][0]=='n')
  {
    if (sscanf(argv[1],"n %s",name)!=1)
      return (CMDERRORCODE);
  }

  strcpy(filename,name);
  strcat(filename,".array");
  if (arraypathes_set)
    stream = FileOpenUsingSearchPaths(filename,"r","arraypathes");
  else
    stream = fileopen(filename,"r");
  if (stream==NULL)
  {
    PrintErrorMessage('E',"LoadArrayCommand","cannot open file");
    return(CMDERRORCODE);
  }
  if (fread((void*)(&nVar),sizeof(INT),1,stream)!=1)
    return(CMDERRORCODE);
  if (nVar>AR_NVAR_MAX)
    return (CMDERRORCODE);
  if (fread((void*)VarDim,sizeof(INT),nVar,stream)!=nVar)
    return(CMDERRORCODE);
  theAR = CreateArray (name,nVar,VarDim);
  if (theAR==NULL)
    return(CMDERRORCODE);
  size = 1;
  for (i=0; i<AR_NVAR(theAR); i++)
    size *= AR_VARDIM(theAR,i);
  if (fread((void*)(theAR->data),sizeof(DOUBLE),size,stream)!=size)
    return(CMDERRORCODE);
  if (fclose(stream))
    return(CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref wrar. */
static INT WriteArrayCommand (INT argc, char **argv)
{
  INT i, Point[AR_NVAR_MAX];
  int iValue;
  double Value;
  char name[128];
  ARRAY *theAR;

  if (argv[1][0]=='n')
  {
    if (sscanf(argv[1],"n %s",name)!=1)
      return (CMDERRORCODE);
  }
  if (ChangeEnvDir("/Array")==NULL)
  {
    PrintErrorMessage('F',"WriteArrayCommand","could not changedir to /Array");
    return(CMDERRORCODE);
  }
  theAR = (ARRAY *)SearchEnv(name,".",theArrayVarID,SEARCHALL);
  if (theAR==NULL)
    return (CMDERRORCODE);

  if (AR_NVAR(theAR) != argc-3)
    return (CMDERRORCODE);
  for (i=0; i<AR_NVAR(theAR); i++)
  {
    if (sscanf(argv[i+2],"%d",&iValue)!=1)
      return (CMDERRORCODE);
    if (iValue<0 || iValue>=AR_VARDIM(theAR,i))
    {
      PrintErrorMessage( 'E', "WriteArrayCommand", "Index Range Error" );
      return (CMDERRORCODE);
    }
    Point[i] = iValue;
  }

  /* write */
  if (sscanf(argv[argc-1],"v %lf",&Value)!=1)
    return (CMDERRORCODE);
  if (WriteArray(theAR,Point,(DOUBLE)Value))
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref rear. */
static INT ReadArrayCommand (INT argc, char **argv)
{
  INT i, Point[AR_NVAR_MAX];
  int iValue;
  DOUBLE value;
  char name[128];
  ARRAY *theAR;

  if (argv[1][0]=='n')
  {
    if (sscanf(argv[1],"n %s",name)!=1)
      return (CMDERRORCODE);
  }
  if (ChangeEnvDir("/Array")==NULL)
  {
    PrintErrorMessage('F',"ReadArrayCommand","could not changedir to /Array");
    return(CMDERRORCODE);
  }
  theAR = (ARRAY *)SearchEnv(name,".",theArrayVarID,SEARCHALL);
  if (theAR==NULL)
    return (CMDERRORCODE);

  if (AR_NVAR(theAR) != argc-2)
    return (CMDERRORCODE);
  for (i=0; i<AR_NVAR(theAR); i++)
  {
    if (sscanf(argv[i+2],"%d",&iValue)!=1)
      return (CMDERRORCODE);
    if (iValue<0 || iValue>=AR_VARDIM(theAR,i))
    {
      PrintErrorMessage( 'E', "ReadArrayCommand", "Index Range Error" );
      return (CMDERRORCODE);
    }
    Point[i] = iValue;
  }

  /* read */
  if (ReadArray(theAR,Point,&value))
    return (CMDERRORCODE);
  if (SetStringValue(":ARRAY_VALUE",(double)value))
    return (CMDERRORCODE);

  return (OKCODE);
}


/** \brief Implementation of \ref clar. */
static INT ClearArrayCommand (INT argc, char **argv)
{
  char name[128];
  ARRAY *theAR;

  if (argv[1][0]=='n')
  {
    if (sscanf(argv[1],"n %s",name)!=1)
      return (CMDERRORCODE);
  }
  if (ChangeEnvDir("/Array")==NULL)
  {
    PrintErrorMessage('F',"ClearArrayCommand","could not changedir to /Array");
    return(CMDERRORCODE);
  }
  theAR = (ARRAY *)SearchEnv(name,".",theArrayVarID,SEARCHALL);
  if (theAR==NULL)
    return (CMDERRORCODE);

  if (ClearArray(theAR))
    return (CMDERRORCODE);

  return (OKCODE);
}

/****************************************************************************/
/** \brief Initialization of the array commands

   This function does initialization of the ug-commands concerning arrays.

   \sa   array, crar, dear, wrar, rear, saar, loar, clar

   @return
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

static INT InitArray (void)
{
  /* install the /Array directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitArray","could not changedir to root");
    return(__LINE__);
  }
  theArrayDirID = GetNewEnvDirID();
  if (MakeEnvItem("Array",theArrayDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitArray","could not install '/Array' dir");
    return(__LINE__);
  }
  theArrayVarID = GetNewEnvVarID();

  /* path to dir for 'array' files */
  arraypathes_set = FALSE;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"arraypathes")==0)
    arraypathes_set = TRUE;

  return (0);
}



/** \brief Implementation of \ref dumpalg. */
static INT DumpAlgCommand(INT argc, char **argv)
{
  INT level, comp;
  VECTOR *v;
  MULTIGRID *theMG;
  GRID *theGrid;
  VECDATA_DESC *v_desc;
  char buffer[1024];

  theMG = currMG;
  if (theMG==NULL)
  {
    PrintErrorMessage('E',"dumpalg","no open multigrid");
    return (CMDERRORCODE);
  }

  v_desc = ReadArgvVecDesc(theMG,"v",argc,argv);
  if (v_desc == NULL)
  {
    PrintErrorMessage('E',"dumpalg","wrong vector specification");
    return (CMDERRORCODE);
  }
  UserWriteF(DISPLAY_NP_FORMAT_SS,"vector displayed",ENVITEM_NAME(v_desc));
  DisplayVecDataDesc(v_desc,~0,buffer);

  for (level=0; level<=TOPLEVEL(theMG); level++)
  {
    theGrid = GRID_ON_LEVEL(theMG,level);

#ifdef ModelP
    for( v=PFIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v) )
#else
    for( v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v) )
#endif
    {
      printf( "Vec key=%d level=%d type=%d pe=%d fine=%d new_def=%d ",
              KeyForObject((KEY_OBJECT*)v),level,VTYPE(v),me,
              FINE_GRID_DOF(v),NEW_DEFECT(v) );
      for( comp=0; comp<VD_NCMPS_IN_TYPE(v_desc,VTYPE(v)); comp++ )
        printf(" %g ",comp,VVALUE(v,VD_CMP_OF_TYPE(v_desc,VTYPE(v),comp)));
      printf("\n");
    }
  }

  /* aus nstools.c
          {
                  MATRIX *m;
              DisplayMatDataDesc(m,buffer);
              fprintf(fptr,"%s",buffer);

              fprintf(fptr,"Matrix:\n");
              for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
              {
                  for (mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
                  {
                      fprintf(fptr,"%d\t%d\n",(int)VINDEX(vec),(int)VINDEX(MDEST(mat)));
                      fprintf(fptr,"f\tt");
                      for (i=0; i<MD_ROWS_IN_RT_CT(m,VTYPE(vec),VTYPE(MDEST(mat))); i++)
                      {
                          for (j=0; j<MD_COLS_IN_RT_CT(m,VTYPE(vec),VTYPE(MDEST(mat))); j+
     +)
                              fprintf(fptr,"\t%.9e",(double)MVALUE(mat,MD_IJ_CMP_OF_RT_CT
          (m,VTYPE(vec),VTYPE(MDEST(mat)),i,j)));
                          fprintf(fptr,"\n\t");
                      }
                      fprintf(fptr,"\n");
                  }
              }
          }
   */

  return(OKCODE);
}

/****************************************************************************/
/* Quick Hack for periodic boundaries                                       */
/****************************************************************************/

#ifdef __PERIODIC_BOUNDARY__
static INT ListPeriodicPosCommand (INT argc, char **argv)
{
  MULTIGRID *theMG = currMG;
  DOUBLE_VECTOR pos;

  if (theMG == NULL)
  {
    UserWrite("ListPeriodicPos: no open multigrid\n");
    return(OKCODE);
  }

#ifdef __THREEDIM__
  if (sscanf(argv[0],"lppos %lf %lf %lf",pos,pos+1,pos+2) != 3)
#else
  if (sscanf(argv[0],"lppos %lf %lf",pos,pos+1) != 2)
#endif
  {
    if (me == master)
      UserWriteF("ListPeriodicPos wrong number of coords\n");
  }

  if (MG_ListPeriodicPos(theMG,0,TOPLEVEL(theMG),pos))
    REP_ERR_RETURN (CMDERRORCODE);

  return (OKCODE);
}

static INT MakePeriodicCommand (INT argc, char **argv)
{
  MULTIGRID *theMG = currMG;

  if (theMG == NULL)
  {
    UserWrite("MakePeriodic: no open multigrid\n");
    return(OKCODE);
  }

  if (MG_GeometricToPeriodic(theMG,0,TOPLEVEL(theMG)))
    REP_ERR_RETURN (CMDERRORCODE);

  return (OKCODE);
}
#endif


/****************************************************************************/
/** \brief Initialization of the commands

   This function does initialization of all ug-commands, using
   'CreateCommand'.
   It initializes 'clock', 'findrange', 'screensize' and 'array'
   commands.

   SEE ALSO:
   commands

   RETURN VALUE:
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitCommands ()
{
  /* quick hack */
#ifdef __PERIODIC_BOUNDARY__
  if (CreateCommand("makeperiodic",       MakePeriodicCommand                             )==NULL) return (__LINE__);
  if (CreateCommand("lppos",                      ListPeriodicPosCommand                  )==NULL) return (__LINE__);
#endif

  /* general commands */
  if (CreateCommand("quit",                       QuitCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("exitug",                     ExitUgCommand                                   )==NULL) return (__LINE__);
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
        #ifdef __TWODIM__
  if (CreateCommand("cnom",                       CnomCommand                                             )==NULL) return (__LINE__);
        #endif

  /* commands for grid management */
  if (CreateCommand("configure",          ConfigureCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("setcurrmg",          SetCurrentMultigridCommand              )==NULL) return (__LINE__);
  if (CreateCommand("new",                        NewCommand                                              )==NULL) return (__LINE__);
  if (CreateCommand("open",                       OpenCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("close",                      CloseCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("save",                       SaveCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("savedomain",         SaveDomainCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("savedata",           SaveDataCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("loaddata",           LoadDataCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("changemc",           ChangeMagicCookieCommand                )==NULL) return (__LINE__);
  if (CreateCommand("level",                      LevelCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("average",        AverageCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("freeaverage",    FreeAverageCommand              )==NULL) return (__LINE__);
  if (CreateCommand("renumber",           RenumberMGCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("smooth",             SmoothMGCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("smoothgrid",     SmoothGridCommand               )==NULL) return(__LINE__);
  if (CreateCommand("ordernodes",         OrderNodesCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("lexorderv",          LexOrderVectorsCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("orderv",             OrderVectorsCommand                     )==NULL) return (__LINE__);
  if (CreateCommand("lineorderv",         LineOrderVectorsCommand                 )==NULL) return (__LINE__);
  if (CreateCommand("revvecorder",        RevertVecOrderCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("shellorderv",        ShellOrderVectorsCommand                )==NULL) return (__LINE__);
  if (CreateCommand("setindex",           SetIndexCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("extracon",           ExtraConnectionCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("check",                      CheckCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("in",                         InsertInnerNodeCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("ngin",                       NGInsertInnerNodeCommand                )==NULL) return (__LINE__);
  if (CreateCommand("bn",                         InsertBoundaryNodeCommand               )==NULL) return (__LINE__);
  if (CreateCommand("ngbn",                       NGInsertBoundaryNodeCommand             )==NULL) return (__LINE__);
  if (CreateCommand("gn",                         InsertGlobalNodeCommand                 )==NULL) return (__LINE__);
  if (CreateCommand("deln",                       DeleteNodeCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("move",                       MoveNodeCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("ie",                         InsertElementCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("ngie",                       NGInsertElementCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("dele",                       DeleteElementCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("refine",             AdaptCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("adapt",                      AdaptCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("fixcoarsegrid",      FixCoarseGridCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("collapse",           CollapseCommand                         )==NULL) return (__LINE__);
  if (CreateCommand("mark",                       MarkCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("find",                       FindCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("select",             SelectCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("wplist",             WindowPictureListCommand                )==NULL) return (__LINE__);
  if (CreateCommand("mglist",             MGListCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("glist",                      GListCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("nlist",                      NListCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("elist",                      EListCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("slist",                      SelectionListCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("rlist",                      RuleListCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("printvalue",         PrintValueCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("vmlist",             VMListCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("convert",        ConvertCommand                  )==NULL) return(__LINE__);
  if (CreateCommand("quality",            QualityCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("makegrid",           MakeGridCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("status",                     StatusCommand                                   )==NULL) return (__LINE__);
#ifdef __THREEDIM__
  if (CreateCommand("fiflel",                     FindFlippedElementsCommand              )==NULL) return (__LINE__);
#endif

#if defined(CAD) && defined(__THREEDIM__)
  if (CreateCommand("cadconvert",     CADGridConvertCommand           )==NULL) return (__LINE__);
#endif

  /* commands for grape */
  if (CreateCommand("grape",                      CallGrapeCommand                                )==NULL) return (__LINE__);
#ifdef _COVISE
  /* commands for covise */
  if (CreateCommand("covise",                     CoviseCommand                                   )==NULL) return (__LINE__);
#endif

  /* commands for window and picture management */
  if (CreateCommand("screensize",         ScreenSizeCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("openwindow",         OpenWindowCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("openppic",       OpenPlacedPicturesCommand           )==NULL) return (__LINE__);
  if (CreateCommand("closewindow",        CloseWindowCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("setcurrwindow",      SetCurrentWindowCommand                 )==NULL) return (__LINE__);
  if (CreateCommand("drawtext",           DrawTextCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("openpicture",        OpenPictureCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("closepicture",       ClosePictureCommand                     )==NULL) return (__LINE__);
  if (CreateCommand("clearpicture",       ClearPictureCommand                     )==NULL) return (__LINE__);
  if (CreateCommand("picframe",           PicFrameCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("setcurrpicture", SetCurrentPictureCommand            )==NULL) return (__LINE__);
  if (CreateCommand("picwin",                     PictureWindowCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("setview",            SetViewCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("cpview",                     CopyViewCommand                                 )==NULL) return (__LINE__);
  if (CreateCommand("vdisplay",           DisplayViewCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("walk",                       WalkCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("walkaround",         WalkAroundCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("zoom",                       ZoomCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("drag",                       DragCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("rotate",             RotateCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("textfac",            TextFacCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("linefac",            LineFacCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("setplotobject",      SetPlotObjectCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("polist",             PlotObjectListCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("plot",                       PlotCommand                                     )==NULL) return (__LINE__);
  if (CreateCommand("findrange",          FindRangeCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("updateDoc",          UpdateDocumentCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("rotmode",            RotModeCommand                                  )==NULL) return (__LINE__);
  if (CreateCommand("cmfn",                       CreateMetafileNameCommand               )==NULL) return (__LINE__);
  if (CreateCommand("setpalette",         SetPaletteCommand                               )==NULL) return (__LINE__);

  /* commands for problem management */
  if (CreateCommand("reinit",             ReInitCommand                                   )==NULL) return (__LINE__);

  /* commands for NumProc management */
  if (CreateCommand("npexecute",          ExecuteNumProcCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("npdisplay",          NumProcDisplayCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("npcreate",           NumProcCreateCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("npinit",                     NumProcInitCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("scnp",                       SetCurrentNumProcCommand                )==NULL) return (__LINE__);

  /* vectors and matrices */
  if (CreateCommand("clear",                      ClearCommand                                    )==NULL) return (__LINE__);
  if (CreateCommand("mflops",         MFLOPSCommand                   )==NULL) return (__LINE__);
  if (CreateCommand("makevdsub",      MakeVDsubCommand                )==NULL) return (__LINE__);

  if (CreateCommand("rand",                       RandCommand                                             )==NULL) return (__LINE__);
  if (CreateCommand("copy",                       CopyCommand                                             )==NULL) return (__LINE__);
  if (CreateCommand("add",                        AddCommand                                              )==NULL) return (__LINE__);
  if (CreateCommand("sub",                        SubCommand                                              )==NULL) return (__LINE__);
  if (CreateCommand("homotopy",       HomotopyCommand                 )==NULL) return(__LINE__);
  if (CreateCommand("interpolate",        InterpolateCommand                              )==NULL) return (__LINE__);

  /* formats */
  if (CreateCommand("newformat",          CreateFormatCommand                             )==NULL) return (__LINE__);
  if (CreateCommand("delformat",          DeleteFormatCommand                             )==NULL) return (__LINE__);
  if (CreateCommand("showpf",             ShowPrintingFormatCommand               )==NULL) return (__LINE__);
  if (CreateCommand("setpf",                      SetPrintingFormatCommand                )==NULL) return (__LINE__);
  if (CreateCommand("createvector",   CreateVecDescCommand            )==NULL) return (__LINE__);
  if (CreateCommand("creatematrix",   CreateMatDescCommand            )==NULL) return (__LINE__);
  if (CreateCommand("freematrix",     FreeMatDescCommand              )==NULL) return (__LINE__);
  if (CreateCommand("symlist",            SymListCommand                      )==NULL) return (__LINE__);

  /* miscellaneous commands */
  if (CreateCommand("setkey",             SetCommandKeyCommand                    )==NULL) return (__LINE__);
  if (CreateCommand("delkey",             DeleteCommandKeyCommand                 )==NULL) return (__LINE__);
  if (CreateCommand("keylist",            ListCommandKeysCommand                  )==NULL) return (__LINE__);
  if (CreateCommand("refreshon",          RefreshOnCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("refreshoff",         RefreshOffCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("machinetest",        MachineTestCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("system",                     SystemCommand                                   )==NULL) return (__LINE__);
  if (CreateCommand("resetCEstat",        ResetCEstatCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("printCEstat",        PrintCEstatCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("heapstat",           HeapStatCommand                             )==NULL) return (__LINE__);
  if (CreateCommand("getheapused",        GetHeapUsedCommand                          )==NULL) return (__LINE__);

  /* commands for debugging */
        #ifdef Debug
  if (CreateCommand("debug",                      DebugCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("trace",                      TraceCommand                                )==NULL) return (__LINE__);
  if (CreateCommand("reperr",             RepErrCommand                               )==NULL) return (__LINE__);
  if (CreateCommand("timing",             TimingCommand                               )==NULL) return (__LINE__);
        #endif
  if (CreateCommand("showconfig",         ShowConfigCommand                           )==NULL) return (__LINE__);

#ifdef ModelP
  /* commands for parallel version */
  if (CreateCommand("lb",                         LBCommand                                               )==NULL) return (__LINE__);
  if (CreateCommand("ptest",                      LBSCommand                                      )==NULL) return (__LINE__);
  if (CreateCommand("lbs",                        LBSCommand                                      )==NULL) return (__LINE__);
  if (CreateCommand("context",            ContextCommand                              )==NULL) return (__LINE__);
  if (CreateCommand("pstat",                      PStatCommand                                )==NULL) return (__LINE__);

#ifdef USE_FAMG
  if (CreateCommand("pamgcheck",      pamgCheckCommand                )==NULL) return(__LINE__);
#endif

#endif /* ModelP */

  /* array commands */
  if (CreateCommand("crar",               CreateArrayCommand                                      )==NULL) return (__LINE__);
  if (CreateCommand("dear",               DeleteArrayCommand                                      )==NULL) return (__LINE__);
  if (CreateCommand("saar",               SaveArrayCommand                                        )==NULL) return (__LINE__);
  if (CreateCommand("loar",               LoadArrayCommand                                        )==NULL) return (__LINE__);
  if (CreateCommand("wrar",               WriteArrayCommand                                       )==NULL) return (__LINE__);
  if (CreateCommand("rear",               ReadArrayCommand                                        )==NULL) return (__LINE__);
  if (CreateCommand("clar",               ClearArrayCommand                                       )==NULL) return (__LINE__);

  if (CreateCommand("dumpalg",            DumpAlgCommand                                  )==NULL) return (__LINE__);

  if (InitClock()                 !=0) return (__LINE__);
  if (InitFindRange()     !=0) return (__LINE__);
  if (InitScreenSize()    !=0) return (__LINE__);
  if (InitArray()                 !=0) return (__LINE__);

  return(0);
}

/** @} */
