// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  order.h                                                                                                       */
/*																			*/
/* Purpose:   order num proc header                                                                     */
/*																			*/
/* Author:	  Klaus Johannsen                                                                           */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Nov 18, 1997 begin                                                                */
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

#ifndef __ORDER__
#define __ORDER__

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

#define ORDER_CLASS_NAME "order"

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct np_order {
  NP_BASE base;                              /* inherits base class             */

  /* data (optional, necessary for calling the generic execute routine)   */
  MATDATA_DESC *A;                       /* matrix symbol                   */

  /* functions */
  INT (*Order)
    (struct np_order *,                      /* pointer to (derived) object     */
    INT level,                           /* level to order                  */
    MATDATA_DESC *A,                     /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_order NP_ORDER;

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

INT     InitOrder (void);

#endif
