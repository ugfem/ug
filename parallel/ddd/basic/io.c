// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      io.c                                                              */
/*                                                                          */
/* Purpose:   routines for I/O used by DDD                                      */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*                                                                          */
/* History:   92/01/29 PrintErrorMessage by Peter Bastian                   */
/*            95/03/21 kb  added PrintString()                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "compiler.h"
#include "ppif.h"
#include "include/dddio.h"
#include "dddi.h"


void (*DDD_UserLineOutFunction)(char *s);



/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_PrintLine                                                 */
/*                                                                          */
/* Purpose:   print interface, all output lines will be routed to here      */
/*                                                                          */
/* Input:     char *s: string to be printed                                 */
/*                                                                          */
/* Output:    none                                                              */
/*                                                                          */
/****************************************************************************/

void DDD_PrintLine (char *s)
{
  /* newline character will be included in s */

  if (DDD_UserLineOutFunction!=NULL)
  {
    DDD_UserLineOutFunction(s);
  }
  else
  {
    PrintHostMessage(s);
  }
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_Flush                                                     */
/*                                                                          */
/* Purpose:   flush output device                                           */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    none                                                              */
/*                                                                          */
/****************************************************************************/

void DDD_Flush (void)
{
  fflush(stdout);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_SyncAll                                                   */
/*                                                                          */
/* Purpose:   flush output devices and synchronize                          */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    none                                                              */
/*                                                                          */
/****************************************************************************/

void DDD_SyncAll (void)
{
  DDD_Flush();
  Synchronize();
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_PrintDebug                                                */
/*                                                                          */
/* Purpose:   print interface for debug output                              */
/*                                                                          */
/* Input:     char *s: string to be printed                                 */
/*                                                                          */
/* Output:    none                                                              */
/*                                                                          */
/****************************************************************************/

void DDD_PrintDebug (char *s)
{
  /* newline character will be included in s */

  DDD_PrintLine(s);

  DDD_Flush();
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_PrintError                                                */
/*                                                                          */
/* Purpose:   print formatted error message on user screen                  */
/*                                                                          */
/* Input:     char class: 'W' Warning, 'E' Error, 'F' Fatal                 */
/*            int errno: error number                                       */
/*            char *text: error message text                                */
/*                                                                          */
/* Output:    none                                                              */
/*                                                                          */
/****************************************************************************/

void DDD_PrintError (char error_class, int error_no, char *text)
{
  char buffer[256];
  char classText[32];

  switch (error_class)
  {
  case 'W' :
    strcpy(classText,"WARNING");
    break;

  case 'E' :
    strcpy(classText,"ERROR");
    break;

  case 'F' :
    strcpy(classText,"FATAL");
    break;

  default :
    strcpy(classText,"USER");
    break;
  }
  sprintf(buffer,"DDD [%03d] %s %05d: %s\n",me,classText,error_no,text);
  DDD_PrintLine(buffer);
}
