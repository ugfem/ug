// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      notify.h                                                      */
/*                                                                          */
/* Purpose:   header file for notification module                           */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/*                                                                          */
/* History:   95/04/04 kb  created                                          */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD_NOTIFY_H__
#define __DDD_NOTIFY_H__


#include "include/ddd.h"


/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

typedef struct _NOTIFY_DESC
{
  DDD_PROC proc;
  size_t size;

} NOTIFY_DESC;


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


void         NotifyInit (void);
void         NotifyExit (void);

NOTIFY_DESC *DDD_NotifyBegin (int);
void         DDD_NotifyEnd (void);
int          DDD_Notify (void);


#endif
