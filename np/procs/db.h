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

#define ORDERED_LIST_CLASS_NAME "ordered_list"

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct np_ordered_list {
  NP_BASE base;                              /* inherits base class             */

  /* functions */
  INT (*PreProcess)
    (struct np_ordered_list *,               /* pointer to (derived) object     */
    INT *);                                      /* result                          */
  INT (*GetListEntry_Index)
    (struct np_ordered_list *,               /* pointer to (derived) object     */
    INT,                                 /* index of entry                  */
    DOUBLE *,                                    /* list entry                      */
    INT *);                                      /* result                          */
  INT (*GetListEntry_NextHigherEntry)
    (struct np_ordered_list *,               /* pointer to (derived) object     */
    DOUBLE,                              /* value                           */
    DOUBLE *,                                    /* list entry                      */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_ordered_list *,               /* pointer to (derived) object     */
    INT *);                                      /* result                          */
};
typedef struct np_ordered_list NP_ORDERED_LIST;

typedef INT (*PreProcessDataBase)(NP_ORDERED_LIST *, INT *);
typedef INT (*GetEntry_Index)(NP_ORDERED_LIST *, INT , DOUBLE *, INT *);
typedef INT (*GetEntry_NextHigherEntry)(NP_ORDERED_LIST *, DOUBLE, DOUBLE *, INT *);
typedef INT (*PostProcessDataBase)(NP_ORDERED_LIST *, INT *);

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

INT InitDb (void);

#endif
