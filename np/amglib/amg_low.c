// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_low.c														*/
/*																			*/
/* Purpose:   handlers, memory management for amg							*/
/*																			*/
/* Author:	  Peter Bastian					                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   28 Jan 1996 Begin												*/
/*            02 Apr 1996 new memory allocation strategy					*/
/*            30 Sep 1997 redesign											*/
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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "amg_header.h"
#include "amg_low.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#undef DEBUG

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

/* RCS_ID
   $Header$
 */

/* user installable print handler */
static AMG_PrintFuncPtr AMG_UserPrint=NULL;
static AMG_MallocFuncPtr AMG_UserMalloc=NULL;

/* log file */
FILE *outFile=NULL;

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   AMG_InstallPrintHandler - install user print function for amg package

   SYNOPSIS:
   typedef void (*AMG_PrintFuncPtr) (char *);
   int AMG_InstallPrintHandler (AMG_PrintFuncPtr print);

   PARAMETERS:
   .  print - Function of type AMG_PrintFuncPtr.

   DESCRIPTION:
   This function allows you to define your own print function, i.e. a function
   that writes a string to the console. This is necessary since printf
   is usually not available in windowed environments. If no handler is ever
   installed then AMG will use the standard fputs function to write out strings.

   RETURN VALUE:
   .n AMG_OK

   D*/
/****************************************************************************/

int AMG_InstallPrintHandler (AMG_PrintFuncPtr print)
{
  AMG_UserPrint = print;
  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_Print - function to write a string

   SYNOPSIS:
   int AMG_Print (char *s);
   int AMG_RedirectToFile (char *name)
   int AMG_RedirectToScreen (void)

   PARAMETERS:
   .  s - string to print.

   DESCRIPTION:
   In order to enable proper working of amg within other codes
   all diagnostic output should be done through these functions.
   The usual way is to use sprintf to provide
   formatted output into a buffer which is then written with AMG_Print
   to the screen. The other two functions allow to redirect the output to
   a text file and back to the screen.

   RETURN VALUE:
   .n AMG_OK

   D*/
/****************************************************************************/

int AMG_Print (char *s)
{
  if (outFile!=NULL)
  {
    fputs(s,outFile);
    return(AMG_OK);
  }
  if (AMG_UserPrint==NULL)
    fputs(s,stdout);
  else
    (*AMG_UserPrint)(s);
  return(AMG_OK);
}

int AMG_RedirectToFile (char *name)
{
  if (outFile!=NULL) return(AMG_OK);
  outFile = fopen(name,"w");
  if (outFile==NULL) return(AMG_FATAL);
  return(AMG_OK);
}

int AMG_RedirectToScreen (void)
{
  fclose(outFile);
  outFile = NULL;
  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_InstallMallocHandler - install user malloc function

   SYNOPSIS:
   typedef void (*AMG_MallocFuncPtr) (size_t n);
   int AMG_InstallMallocHandler (AMG_MallocFuncPtr mall);

   PARAMETERS:
   .  mall - pointer to user definable malloc function.

   DESCRIPTION:
   This allows to redefine the memory allocation function used by amg
   in order to be able to use amg within other programs having
   their own memory management. All memory allocation is usually static
   within amg, i.e. memory is allocated but never freed again.

   RETURN VALUE:
   .n AMG_OK
   .n AMG_FATAL
   D*/
/****************************************************************************/

int AMG_InstallMallocHandler (AMG_MallocFuncPtr mall)
{
  AMG_UserMalloc = mall;
  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_Malloc - malloc function for amg

   SYNOPSIS:
   void *AMG_Malloc (size_t nbytes);

   PARAMETERS:
   .  size_t - nbytes.

   DESCRIPTION:
   All dynamic memory allocation in amg is done thru this function.
   It allows redefinition.

   RETURN VALUE:
   .n pointer to memory or NULL

   D*/
/****************************************************************************/

void *AMG_Malloc (size_t nbytes)
{
  if (AMG_UserMalloc==NULL)
    return(malloc(nbytes));
  else
    return(AMG_UserMalloc(nbytes));
}
