// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pfile.h                                                       */
/*                                                                          */
/* Purpose:   parallel file output	                                                        */
/*                                                                          */
/* Author:	  Peter Bastian                                                                                     */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
/*																			*/
/* History:   29.06.95 begin, ug version 3.0								*/
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __PFILE__
#define __PFILE__

#include <limits.h>
#include "compiler.h"

#define PFILE_BUFFER_SIZE       16384
#define PFILE_MAX_TREE          32
#define PFILE_MAX_INT       INT_MAX

typedef struct {
  INT nchars;                                                           /* number of characters in buffer	*/
  INT first_key, last_key;                              /* for sorted output				*/
} PFILE_STATE ;

typedef struct {
  FILE *stream;                                                 /* the output file on master		*/
  PFILE_STATE local;                                    /* the local state					*/
  INT valid_state[PFILE_MAX_TREE+1];    /* valid flags for states			*/
  PFILE_STATE state[PFILE_MAX_TREE+1];  /* all downtree states	and own		*/
  INT robin;                                                            /* round robin for minimum search	*/
        #ifdef ModelP
  char buffer[PFILE_BUFFER_SIZE];               /* local buffer						*/
  char buffer2[PFILE_BUFFER_SIZE];              /* local buffer						*/
        #endif
} PFILE ;

typedef struct {
  INT nints;
  INT nfloats;
  INT nbytes;
  INT first_key, last_key;
} PFILE_STATE_BIN;

typedef struct {
  FILE *stream;
  PFILE_STATE_BIN local;
  INT valid_state[PFILE_MAX_TREE+1];
  PFILE_STATE_BIN state[PFILE_MAX_TREE+1];
  INT robin;
#ifdef ModelP
  INT buf_INT[PFILE_BUFFER_SIZE];
  INT buf_INT2[PFILE_BUFFER_SIZE];
  FLOAT buf_FLOAT[PFILE_BUFFER_SIZE];
  FLOAT buf_FLOAT2[PFILE_BUFFER_SIZE];
  unsigned char buf_BYTE[PFILE_BUFFER_SIZE];
  unsigned char buf_BYTE2[PFILE_BUFFER_SIZE];
#endif
} PFILE_BIN ;

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

PFILE *pfile_open               (char *name);
INT    pfile_master_puts        (PFILE *pf, char *s);
INT    pfile_puts               (PFILE *pf, char *s);
INT    pfile_tagged_puts        (PFILE *pf, char *s, INT key);
INT    pfile_sync               (PFILE *pf);
INT    pfile_close              (PFILE *pf);
PFILE_BIN *pfile_open_bin           (char *name);
INT        pfile_tagged_write_INT   (PFILE_BIN *pf, INT *values, int n, INT key);
INT        pfile_tagged_write_FLOAT (PFILE_BIN *pf, FLOAT *values, int n, INT key);
/* TODO: pfile_tagged_write_BYTE works only in seq Mode !!! */
INT        pfile_tagged_write_BYTE  (PFILE_BIN *pf, unsigned char *values, int n, INT key);
INT        pfile_sync_bin           (PFILE_BIN *pf);
INT        pfile_close_bin          (PFILE_BIN *pf);

#endif
