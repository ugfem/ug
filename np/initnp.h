// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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

#include "compiler.h"

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

/* initialisation of the numerics module */
INT InitNumerics (void);


#ifdef __cplusplus
}  /* namespace UG{2|3}d */
#endif
