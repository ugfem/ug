// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugrefine2d.h													*/
/*																			*/
/* Purpose:   unstructured grid refinement header file						*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					*/
/*																			*/
/* History:   09.03.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UGREFINE2D__
#define __UGREFINE2D__

#include "compiler.h"

#ifndef __GM__
#include "gm.h"
#endif

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
/* macros																	*/
/*																			*/
/****************************************************************************/

#define LEAFELEM(e)      (NSONS(e)==0)
#define LEAFNODE(n)      (SONNODE(n)==NULL)
#define REGULAR(e)               (ECLASS(e)==REGULAR_CLASS)
#define IRREGULAR(e)     (ECLASS(e)==IRREGULAR_CLASS)

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT                      InitRefine2d           (void);

#endif
