// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  iter.h                                                                                                        */
/*																			*/
/* Purpose:   definition of the iter num proc type			                        */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   November 29, 1996                                                                         */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __ITER__
#define __ITER__

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

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct np_iter {

  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  VECDATA_DESC *c;                       /* correction                      */
  VECDATA_DESC *b;                       /* defect                          */
  MATDATA_DESC *A;                       /* matrix                          */

  /* functions */
  INT (*PreProcess)
    (struct np_iter *,                       /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* correction vector               */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *,                                       /* baselevel used by iter          */
    INT *);                                      /* result                          */
  INT (*Iter)
    (struct np_iter *,                       /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* correction vector               */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_iter *,                       /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* correction vector               */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_iter NP_ITER;

typedef INT (*PreProcessIterProcPtr) \
  (NP_ITER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *, INT *);
typedef INT (*IterProcPtr) \
  (NP_ITER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*PostProcessIterProcPtr) \
  (NP_ITER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for iter num procs */
INT NPIterInit (NP_ITER *theNP, INT argc , char **argv);

/* generic display function for iter num procs */
INT NPIterDisplay (NP_ITER *theNP);

/* generic execute function for iter num procs */
INT NPIterExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard iter num proc type */
INT InitIter (void);

#endif
