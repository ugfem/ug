// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      fifo.h                                                        */
/*                                                                            */
/* Purpose:   header file for general purpose fifo                            */
/*                                                                            */
/* Author:      Peter Bastian                                                 */
/*              Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen    */
/*              Universitaet Heidelberg                                        */
/*              Im Neuenheimer Feld 368                                        */
/*              6900 Heidelberg                                                */
/*                                                                            */
/* History:   30.01.92 begin, ug version 2.0                                */
/*                                                                            */
/* Revision:  07.09.95                                                      */
/*                                                                            */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                            */
/* auto include mechanism and other include files                            */
/*                                                                            */
/****************************************************************************/

#ifndef __FIFO__
#define __FIFO__


#ifndef __COMPILER__
#include "compiler.h"
#endif

/****************************************************************************/
/*                                                                            */
/* data structures exported by the corresponding source file                */
/*                                                                            */
/****************************************************************************/

typedef struct {
  INT start,end,size,used;
  void **elements;
} FIFO ;


/****************************************************************************/
/*                                                                            */
/* function declarations                                                    */
/*                                                                            */
/****************************************************************************/

INT     fifo_init    (FIFO *myfifo, void *buffer, INT size);
void    fifo_clear    (FIFO *myfifo);
INT     fifo_empty    (const FIFO *myfifo);
INT     fifo_full    (const FIFO *myfifo);
INT     fifo_in     (FIFO *myfifo, void *newelement);
void    *fifo_out    (FIFO *myfifo);

#endif
