// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup low Low-Level Stuff
 * \ingroup ug
 *
 *
 * @{
 *
 */
/**
 * \file initlow.h
 */
/****************************************************************************/
/*                                                                          */
/* File:      initlow.h                                                     */
/*                                                                          */
/* Purpose:   call the init routines of the low module (header file)        */
/*                                                                          */
/* Author:    Henrik Rentz-Reichert                                         */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   27.02.95 begin, ug version 3.0                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

#ifndef __INITLOW__
#define __INITLOW__

#include "compiler.h"


#include "namespace.h"

START_UG_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/** \brief Initialisation of the low module */
INT InitLow (void);

INT ExitLow();

END_NAMESPACE

/** @} */
#endif
