// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amgtransfer.h                                                                                                 */
/*																			*/
/* Purpose:   initialization for algebraic multigrid	                    */
/*                                                                          */
/* Author:	  Nicolas Neuss                                                                                     */
/*			  Institut fuer Angewandte Mathematik                           */
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 294										*/
/*			  69120 Heidelberg												*/
/*			  email: neuss@iwr.uni-heidelberg.de			                        */
/*																			*/
/* History:   1994-1995 in old ug2.0							            */
/*            May 1997  in new ug3.7                                        */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __AMGTRANSFER__
#define __AMGTRANSFER__

#include "np.h"
#include "transfer.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT AMGTransferInit       (NP_BASE *np, INT argc , char **argv);
INT AMGTransferDisplay    (NP_BASE *theNP);
INT AMGTransferPreProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *result);
INT InitAMGTransfer       (void);

#endif
