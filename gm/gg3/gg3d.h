// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup gg3 A 3d Grid Generator
 * \ingroup gm
 */
/*! \file gg3d.h
 * \ingroup gg3
 */

/****************************************************************************/
/*                                                                          */
/* File:      gg3d.h                                                        */
/*                                                                          */
/* Purpose:   header file for the 3d grid generator                         */
/*                                                                          */
/* Author:    Christian Wieners                                             */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart, Germany                                      */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   18 March 96 begin, ug version 3.2                             */
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

#ifndef __GG3D__
#define __GG3D__

#include "compiler.h"
#include "gm.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT GenerateGrid3d (MULTIGRID *theMG, MESH *mesh, DOUBLE h, INT smooth,
                    INT display, INT coeff, INT from, INT to, INT prism, INT save,
                    INT argc, char **argv);
int Get_Local_h(double *in, double *out);

END_NAMESPACE

#endif
