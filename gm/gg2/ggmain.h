// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ggmain.h                                                      */
/*                                                                          */
/* Purpose:   header file for ggmain			                                                */
/*                                                                          */
/* Author:    Wolfgang Hoffmann, Henrik Renz-Reichert	                    */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart, Germany										*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   08.03.94 begin, ug version 2.2                                */
/*                15.10.95 implemented in ug31                                  */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __GGMAIN__
#define __GGMAIN__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __GM__
#include "gm.h"
#endif

#ifndef __GGM__
#include "ggm.h"
#endif

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT GenerateBnodes  (MULTIGRID *theMG, COORD RelRasterSize,
                     DOUBLE h_global, INT meshsizecoeffno);
INT GenerateGrid (MULTIGRID *theMG, GG_ARG *MyArgs, GG_PARAM *param, MESH *mesh, CoeffProcPtr coeff);
INT InitGG (void);

#endif
