// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file ugio.h
 * \ingroup gm
 */


/****************************************************************************/
/*                                                                                                                                                      */
/* File:          ugio.h                                                                                                                */
/*                                                                                                                                                      */
/* Purpose:   ug input/output header file                                       */
/*                                                                                                                                                      */
/* Author:        Peter Bastian                                                                                                 */
/*                        Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*                        Universitaet Heidelberg                                                                               */
/*                        Im Neuenheimer Feld 368                                                                               */
/*                        6900 Heidelberg                                                                                               */
/*                        internet: ug@ica3.uni-stuttgart.de                                            */
/*                                                                                                                                                      */
/* History:   15.04.92 begin, ug version 2.0                                                            */
/*                                                                                                                                                      */
/* Remarks:                                                                                                                             */
/*                                                                                                                                                      */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                                                                                                      */
/* auto include mechanism and other include files                                                       */
/*                                                                                                                                                      */
/****************************************************************************/

#ifndef __UGIO__
#define __UGIO__

#include "compiler.h"
#include "gm.h"

/**************************************************/
/* A namespace for the c++ version                */
/**************************************************/
#ifdef __cplusplus
#ifdef __TWODIM__
namespace UG2d {
#else
namespace UG3d {
#endif
#endif


/****************************************************************************/
/*                                                                                                                                                      */
/* function declarations                                                                                                        */
/*                                                                                                                                                      */
/****************************************************************************/

INT InitUgio ();

#ifdef __cplusplus
}  /* namespace UG{2|3}d */
#endif

#endif
