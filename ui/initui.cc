// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initui.c														*/
/*																			*/
/* Purpose:   init file to init user interface								*/
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de				                        */
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
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

#include <config.h>

/* ANSI-C includes */
#include <cstdio>
#include <cstring>

/* low module */
#include "misc.h"
#include "defaults.h"
#include "general.h"

/* user interface module */
#include "uginterface.h"
#include "helpmsg.h"
#include "cmdint.h"
#include "ugstruct.h"
#include "commands.h"
#include "mmio.h"
#include "tecplot.h"
#include "avs.h"
#include "dataexplorer.h"
#include "fieldio.h"

/* own header */
#include "initui.h"

USING_UG_NAMESPACES

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* Function:  InitUi														*/
/*																			*/
/* Purpose:   Init ui module												*/
/*																			*/
/* Input:	  int argc, char **argv: of main-call							*/
/*																			*/
/* Output: INT 0: ok														*/
/*			   1: init failed												*/
/*																			*/
/****************************************************************************/

INT NS_DIM_PREFIX InitUi (INT argc, char **argv)
{
  INT err;

  /* init ug interface */
  if ((err=InitUgInterface())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init cmdline */
  if ((err=InitCmdline())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init help mechanism */
  if ((err=InitHelpMsg())!=0)
    PrintErrorMessage('W',"InitUi","help mechanism not working properly");

  /* init command interpreter */
  if ((err=InitCommandInterpreter(argc,argv))!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* init commands of ug's command line interface */
  if ((err=InitCommands())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* tecplot output */
  if ((err=InitTecplot())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* avs output */
  if ((err=InitAVS())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* DataExplorer output */
  if ((err=InitDataExplorer())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* matrix market input/output */
  if ((err=InitMMIO())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* fieldio module */
  if ((err=InitFieldIO())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* return to application */
  return(0);
}


INT NS_DIM_PREFIX ExitUi()
{
  INT err;

  /* shut down command interpreter */
  if ((err=ExitCommandInterpreter())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  return 0;
}
