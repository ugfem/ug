// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup gm The Grid Manager
 * \ingroup ug
 */
/*! \file initgm.h
 * \ingroup gm
 */

/****************************************************************************/
/*                                                                                                                                                      */
/* File:          initgm.h                                                                                                      */
/*                                                                                                                                                      */
/* Purpose:   call the init routines of the grid manager module                         */
/*                        (header file)                                                                                                 */
/*                                                                                                                                                      */
/* Author:        Henrik Rentz-Reichert                                                                                 */
/*                        Institut fuer Computeranwendungen III                                                 */
/*                        Universitaet Stuttgart                                                                                */
/*                        Pfaffenwaldring 27                                                                                    */
/*                        70569 Stuttgart                                                                                               */
/*                        email: ug@ica3.uni-stuttgart.de                                                           */
/*                                                                                                                                                      */
/* History:   27.02.95 begin, ug version 3.0                                                            */
/*                                                                                                                                                      */
/* Remarks:                                                                                                                             */
/*                                                                                                                                                      */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

#ifndef __INITGM__
#define __INITGM__

#include "compiler.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* initialisation of the gm module */
INT InitGm ();

/* Clean up of the gm module */
INT ExitGm ();

END_NAMESPACE

#endif
