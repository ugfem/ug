// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup gg2 Another Grid Generator
 * \ingroup gm
 */
/*! \file ggmain.h
 * \ingroup gg2
 */

/****************************************************************************/
/*                                                                          */
/* File:      ggmain.h                                                      */
/*                                                                          */
/* Purpose:   header file for ggmain                                        */
/*                                                                          */
/* Author:    Wolfgang Hoffmann, Henrik Renz-Reichert                       */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart, Germany                                      */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   08.03.94 begin, ug version 2.2                                */
/*            15.10.95 implemented in ug31                                  */
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

#ifndef __GGMAIN__
#define __GGMAIN__

#include "compiler.h"
#include "gm.h"
#include "ggm.h"

#include "namespace.h"

START_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT GenerateBnodes  (MULTIGRID *theMG, DOUBLE RelRasterSize,
                     DOUBLE h_global, INT meshsizecoeffno);
INT GenerateGrid (MULTIGRID *theMG, GG_ARG *MyArgs, GG_PARAM *param, MESH *mesh, CoeffProcPtr coeff, INT Single_Mode, INT display);
INT InitGG (void);

END_NAMESPACE

#endif
