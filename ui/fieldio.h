// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      fieldio.h                                                     */
/*                                                                          */
/* Purpose:   Field I/O commands                                            */
/*                                                                          */
/* Author:    Michael Lampe                                                 */
/*            IWR - Technische Simulation                                   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            Email: Michael.Lampe@iwr.uni-heidelberg.de                    */
/*                                                                          */
/* History:   20011002 begin, ug3.8                                         */
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

#ifndef __FIELDIO__
#define __FIELDIO__

#include "compiler.h"

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT InitFieldIO(void);
INT LoadFieldCommand(INT argc, char **argv);

#endif
