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
/*  definition of constants                                                 */
/*                                                                          */
/****************************************************************************/


enum {
  DDD_MODULE_MGR,
  DDD_MODULE_XFER,
  DDD_MODULE_IDENT,
  DDD_MODULE_IF,
  DDD_MODULES
};


/* for DDD_MODULE_IDENT */
/* times */
enum {
  T_PREPARE,
  T_PREPARE_SORT,
  T_QSORT_TUPEL,
  T_RESOLVE_DEP,
  T_QSORT_LOI,
  T_BUILD_GRAPH,
  T_CONSTRUCT_ARRAY,
  T_COMM_AND_IDENT,
  T_BUILD_IF
};
/* numbers */
enum {
  N_PARTNERS
};


/* for DDD_MODULE_IF */
/* times */
enum {
  T_CREATE_COLLECT,
  T_CREATE_SORT,
  T_CREATE_BUILD,
  T_CREATE_SHORTCUT,
  T_CREATE_COMM
};
/* numbers */
/*enum {
   };
 */



/* for DDD_MODULE_XFER */
/* times */
enum {
  T_XFER_PREP_CMDS,
  T_XFER_PREP_MSGS,
  T_XFER_PACK_SEND,
  T_XFER_WHILE_COMM,
  T_XFER_WAIT_RECV,
  T_XFER_UNPACK,
  T_XFER_PREP_CPL,
  T_XFER_CPLMSG,
  T_XFER_BUILD_IF
};
/* numbers */
/*enum {
   };
 */


/****************************************************************************/
/*                                                                          */
/*  definition of exported global variables                                 */
/*                                                                          */
/****************************************************************************/




/****************************************************************************/
/*                                                                          */
/*  declaration of DDD-addon functional interface                           */
/*                                                                          */
/****************************************************************************/


/*
        Statistic Evaluation and Performance Measurements
 */
double    DDD_StatClock (int module, int index);
long      DDD_StatCount (int module, int index);
char *    DDD_StatClockDesc (int module, int index);
char *    DDD_StatCountDesc (int module, int index);


#ifdef __cplusplus
}
#endif

#endif
