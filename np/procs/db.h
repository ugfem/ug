// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  db.h                                                                                                          */
/*																			*/
/* Purpose:   header data base interface                                                                */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Sep 11, 1997 begin                                                                */
/*																			*/
/* Remarks:   not finished!                                                                     */
/*																			*/
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __DB__
#define __DB__

#include "np.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define DATA_BASE_CLASS_NAME "data_base"

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct np_data_base {
  NP_BASE base;                              /* inherits base class             */

  /* functions */
  INT (*PreProcess)
    (struct np_data_base *,                  /* pointer to (derived) object     */
    INT *);                                      /* result                          */
  INT (*GetList)
    (struct np_data_base *,                  /* pointer to (derived) object     */
    DOUBLE **,                                   /* ptr to list                     */
    INT *,                                       /* size of list                    */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_data_base *,                  /* pointer to (derived) object     */
    INT *);                                      /* result                          */
};
typedef struct np_data_base NP_DATA_BASE;

typedef INT (*PreProcessDataBase)(NP_DATA_BASE *, INT *);
typedef INT (*GetList)(NP_DATA_BASE *, DOUBLE **List, INT *, INT *);
typedef INT (*PostProcessDataBase)(NP_DATA_BASE *, INT *);

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

INT     InitDb (void);

INT DB_Init (NP_BASE *theNP, INT argc, char **argv);
INT DB_Display (NP_BASE *theNP);


#endif
