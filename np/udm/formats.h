// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      formats.h                                                     */
/*                                                                          */
/* Purpose:   header file for format definition				                */
/*                                                                          */
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: henrik@ica3.uni-stuttgart.de							*/
/*			  fon: 0049-(0)711-685-7007										*/
/*			  fax: 0049-(0)711-685-7000										*/
/*																			*/
/* History:   27.03.95 begin, ug version 3.0								*/
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __FORMATS__
#define __FORMATS__

#include "gm.h"
#include "np.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT DisplayPrintingFormat                               (void);
INT SetPrintingFormatCmd                                (const MULTIGRID *mg, INT argc, char **argv);

VECDATA_DESC *CreateVecDescOfTemplate   (MULTIGRID *theMG,
                                         const char *name, const char *template);
MATDATA_DESC *CreateMatDescOfTemplate   (MULTIGRID *theMG,
                                         const char *name, const char *template);

INT VDsubDescFromVT                                             (const VECDATA_DESC *vd, const char *tplt, const char *sub, VECDATA_DESC **subvd);
INT MDsubDescFromVT                                             (const MATDATA_DESC *md, const char *tplt, const char *sub, MATDATA_DESC **submd);
INT VDinterfaceDesc                                             (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi);

INT CreateFormatCmd                                             (INT argc, char **argv);
INT CreateVecDescCmd                            (MULTIGRID *theMG, INT argc, char **argv);
INT CreateMatDescCmd                                (MULTIGRID *theMG, INT argc, char **argv);

INT InitFormats                                                 (void);

#endif
