// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/*            DDD V1.1                                                      */
/*                                                                          */
/* File:      dddaddon.h                                                    */
/*                                                                          */
/* Purpose:   header file for ddd additional functionality                  */
/*              i.e.:                                                       */
/*                    1) Statistical evaluation and performance measures    */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/*                                                                          */
/* History:   95/01/18 kb  begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __DDDADDON__
#define __DDDADDON__

#include "dddresources.h"


#ifdef __cplusplus
extern "C" {
#endif


/****************************************************************************/
/*                                                                          */
/*  definition of exported global variables                                 */
/*                                                                          */
/****************************************************************************/


/* description for StatClock entries after Xfer measurement */
extern char *DDD_StatXfer[];



/****************************************************************************/
/*                                                                          */
/*  declaration of DDD-addon functional interface                           */
/*                                                                          */
/****************************************************************************/


/*
        Statistic Evaluation and Performance Measurements
 */
double    DDD_StatClock (int index);
int       DDD_StatCount (int index);


#ifdef __cplusplus
}
#endif

#endif
