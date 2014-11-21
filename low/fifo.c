// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief Functions to handle a fifo (first in, first out) structure
 */

/****************************************************************************/
/*	                                                                        */
/* File:      fifo.c                                                        */
/*                                                                          */
/* Purpose:   general purpose first in first out queue                      */
/*                                                                          */
/* Author:      Peter Bastian                                               */
/*              Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen */
/*              Universitaet Heidelberg                                     */
/*              Im Neuenheimer Feld 368                                     */
/*              6900 Heidelberg                                             */
/*                                                                          */
/* History:   30.01.92 begin, ug version 2.0                                */
/*                                                                          */
/* Revision:  07.09.95                                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*              system include files                                        */
/*              application include files                                   */
/*                                                                          */
/****************************************************************************/

#include <config.h>
#include <stdio.h>

#include "ugtypes.h"
#include "general.h"
#include "fifo.h"

USING_UG_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/** \brief Initialize fifo data structure

   \param myfifo - pointer to a fifo record
   \param buffer - pointer to a memory area for the fifo
   \param size -   size of the buffer in bytes

   This function initializes 'fifo' (first in, first out) data structure.

   RETURN VALUE:
   .n     size of fifo record
   .n     0 if size lower or equal to zero
 */
/****************************************************************************/

INT NS_PREFIX fifo_init (FIFO *myfifo, void *buffer, INT size)
{
  myfifo->size = size / sizeof(void *);
  if (myfifo->size<=0) return(0);
  myfifo->elements = (void **) buffer;
  myfifo->start = myfifo->end = myfifo->used = 0;
  return(myfifo->size);
}

/****************************************************************************/
/** \brief Reset a previously initialized FIFO structure

   \param myfifo - pointer to FIFO structure

   This function resets a previously initialized FIFO structure.
 */
/****************************************************************************/

void NS_PREFIX fifo_clear (FIFO *myfifo)
{
  myfifo->start = myfifo->end = myfifo->used = 0;
}

/****************************************************************************/
/** \brief Test if fifo is empty

   \param myfifo - pointer to fifo structure

   This function tests if fifo is empty.

   RETURN VALUE:
   .n     0 if fifo is not empty
   .n     1 if fifo is empty
 */
/****************************************************************************/

INT NS_PREFIX fifo_empty (const FIFO *myfifo)
{
  if (myfifo->used==0)
    return(1);
  else
    return(0);
}

/****************************************************************************/
/** \brief Test if fifo is full

   \param myfifo - pointer to FIFO structure

   This function tests if fifo is full.

   RETURN VALUE:
   .n     0 if fifo is not empty
   .n     1 if fifo is empty
 */
/****************************************************************************/

INT NS_PREFIX fifo_full (const FIFO *myfifo)
{
  if (myfifo->used==myfifo->size)
    return(1);
  else
    return(0);
}

/****************************************************************************/
/** \brief Insert an element in the fifo

   \param myfifo - pointer to FIFO structure
   \param newelement - pointer to the new element (any type!)

   This function inserts an element in the fifo.

   RETURN VALUE:
   .n     1 if error occurred.
   .n     0 if OK
 */
/****************************************************************************/

INT NS_PREFIX fifo_in (FIFO *myfifo, void *newelement)
{
  if (myfifo->used<myfifo->size)
  {
    (myfifo->elements)[myfifo->end] = newelement;
    myfifo->end = (myfifo->end+1)%myfifo->size;
    myfifo->used++;
    return(0);
  }
  else
    return(1);
}

/****************************************************************************/
/** \brief Get an element from the fifo

   \param myfifo - pointer to fifo structure

   This function gets an element from the fifo.

   RETURN VALUE:
   .n        pointer to the element to take out from the fifo
   .n        NULL if fifo is empty
 */
/****************************************************************************/

void *NS_PREFIX fifo_out (FIFO *myfifo)
{
  INT i;

  if (myfifo->used==0)
    return(NULL);
  else
  {
    i = myfifo->start;
    myfifo->start = (myfifo->start+1)%myfifo->size;
    myfifo->used--;
    return (myfifo->elements[i]);
  }
}
