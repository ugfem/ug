// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  graphics.c                                                                                                    */
/*																			*/
/* Purpose:   common functions for graphical interfaces                         */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken                                                                     */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   971216 begin													*/
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
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <stdarg.h>

/* low module */
#include "ugtypes.h"
/*
   #include "fileopen.h"
   #include "misc.h"
   #include "ugenv.h"
   #include "defaults.h"
   #include "debug.h"
   #include "ugstruct.h"
 */
#include "general.h"

/* gm module */
#include "gm.h"


#include "initgraph.h"
#ifdef _COVISE
#include "coviseif.h"
#endif
#ifdef _PV3
#include "pv3if.h"
#endif

/* own header */
#include "graphics.h"

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
/*D
   InitGraphics - Initialize all graphical interfaces at startup

   SYNOPSIS:
   INT InitGraphics (void);


   DESCRIPTION:
   This function initializes all graphical interfaces at startup.
   It must be extended when a new graphical interfaces is added.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if some error occured.
   D*/
/****************************************************************************/


INT NS_DIM_PREFIX InitGraphics (void)
{
  INT error;

  /* init UG-graphics */
  error = InitUGGraph();
  if (error!=0)
    return(error);

#ifdef _COVISE
  /* init Covise interface */
  error = InitCoviseIF();
  if (error!=0)
    return(error);
#endif

#ifdef _PV3
  error = InitPV3();
  if (error!=0)
    return(error);
#endif

  return(0);         /* no error */
}



INT NS_DIM_PREFIX ExitGraphics (void)
{
        #ifdef _COVISE
  INT error;

  /* close Covise interface */
  error = ExitCoviseIF();
  if (error!=0)
    return(error);
        #endif

  return(0);         /* no error */
}
