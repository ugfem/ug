// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup np Numerical Procedures
 * \ingroup ug
 */
/*! \file initnp.h
 * \ingroup np
 */

/****************************************************************************/
/*                                                                                                                                                      */
/* File:          initnumerics.h                                                                                                */
/*                                                                                                                                                      */
/* Purpose:   call the init routines of the numerics module                                     */
/*                        (header file)                                                                                                 */
/*                                                                                                                                                      */
/* Author:        Klaus Johannsen                                                                                               */
/*                        Institut fuer Computeranwendungen III                                                 */
/*                        Universitaet Stuttgart                                                                                */
/*                        Pfaffenwaldring 27                                                                                    */
/*                        70569 Stuttgart                                                                                               */
/*                        email: ug@ica3.uni-stuttgart.de                                                           */
/*                                                                                                                                                      */
/* History:   24.04.95 begin, ug version 3.0                                                            */
/*                                                                                                                                                      */
/* Remarks:                                                                                                                             */
/*                                                                                                                                                      */
/****************************************************************************/


/* RCS_ID
   $Header$
 */
#ifndef __INITNP__
#define __INITNP__

#include "compiler.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* initialisation of the numerics module */
INT InitNumerics();


END_UGDIM_NAMESPACE

#endif
