// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup ug The UG Kernel
 */
/*! \file initug.h
 * \ingroup ug
 */


/****************************************************************************/
/*                                                                          */
/* File:      initug.h                                                      */
/*                                                                          */
/* Purpose:   call the init routines of the ug modules (header file)        */
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

#ifndef __INITUG__
#define __INITUG__

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
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/** \brief Initialisation of the ug library */
INT InitUg (int *argcp, char ***argvp);

/** \brief Finalisation of the ug library */
INT ExitUg (void);



#ifdef __cplusplus
}  /* namespace UG{2|3}d */
#endif

#endif
