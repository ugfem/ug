// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      dddstr.h                                                      */
/*                                                                          */
/* Purpose:   ddd string constants                                          */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/*                                                                          */
/* History:   970903 kb  begin                                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */


#ifndef __DDDSTR_H__
#define __DDDSTR_H__


/****************************************************************************/
/*                                                                          */
/* string constants for error codes                                         */
/*                                                                          */
/****************************************************************************/


/* error strings for Identify */

#define ERR_ID_NOTIFY_FAILED \
  "Notify failed in Ident-ConsCheck"
#define ERR_ID_WRONG_MODE \
  "wrong Ident-mode (currently in %s, expected %s)"
#define ERR_ID_SAME_TUPEL \
  "same identification tupel for objects %08x and %08x"
#define ERR_ID_OBJ_CYCLE \
  "IdentifyObject-cycle, objects %08x and %08x"
#define ERR_ID_NOMEM_RESOLV \
  STR_NOMEM " in ResolveDependencies"
#define ERR_ID_UNKNOWN_OPT \
  "unknown OPT_IDENTIFY_MODE"
#define ERR_ID_NOMEM_SORT \
  STR_NOMEM " in IdentifySort"
#define ERR_ID_DIFF_IDENT \
  "Identify: no Ident-calls from proc %d, expected %d"
#define ERR_ID_DIFF_N_IDENT \
  "Identify: %d Ident-calls from proc %d, expected %d"
#define ERR_ID_DIFF_N_OBJECTS \
  "Identify: %d identified objects from proc %d, expected %d"
#define ERR_ID_ERRORS \
  "found errors in IdentifyEnd()"
#define ERR_ID_OK \
  "Ident-ConsCheck level 0: ok"
#define ERR_ID_ABORT_END \
  "DDD_IdentifyEnd() aborted"
#define ERR_ID_NOMEM_IDENT_END \
  STR_NOMEM " in DDD_IdentifyEnd"
#define ERR_ID_INCONS_TUPELS \
  "inconsistent tupels, gid %08x on %d, gid %08x on %d, in IdentifyEnd()"
#define ERR_ID_CANT_RECV \
  "couldn't receive message from %d in IdentifyEnd()"
#define ERR_ID_NO_BEGIN \
  "Missing DDD_IdentifyBegin(), aborted"
#define ERR_ID_NOT_WITH_ME \
  "cannot identify %08x with myself"
#define ERR_ID_NOT_WITH_PROC \
  "cannot identify %08x with processor %d"
#define ERR_ID_NOMEM_IDENTRY \
  STR_NOMEM "in IdentifyIdEntry"
#define ERR_ID_NOMEM_IDNUMBER \
  STR_NOMEM " in DDD_IdentifyNumber"
#define ERR_ID_NOMEM_IDSTRING \
  STR_NOMEM "in DDD_IdentifyString"
#define ERR_ID_NOMEM_IDOBJ \
  STR_NOMEM " in DDD_IdentifyObject"
#define ERR_ID_ABORT_BEGIN \
  "DDD_IdentifyBegin() aborted."





/****************************************************************************/

#endif
