// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  data_io.h														*/
/*																			*/
/* Purpose:   header file for data_io.c			                                                        */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   16.12.96 begin												*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/* switch */
#define __DIO_USE_IN_UG__

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __DATA_IO__
#define __DATA_IO__

#include <stdio.h>


/****************************************************************************/
/*																			*/
/* configuration of interface                                                                                           */
/*																			*/
/****************************************************************************/

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

INT LoadData (MULTIGRID *theMG, char *FileName, char *type, INT number, INT n, VECDATA_DESC **theVDList);
INT SaveData (MULTIGRID *theMG, char *FileName, INT rename, char *type, INT number, DOUBLE time, DOUBLE dt, DOUBLE ndt, INT n, VECDATA_DESC **theVDList, EVALUES **theEVal, EVECTOR **theEVec, char **NameList);

#endif
