// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  project.h                                                                                                     */
/*																			*/
/* Purpose:   projection header                                                                         */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Nov 5, 1997 begin                                                                 */
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

#ifndef __PROJECT__
#define __PROJECT__

#include "np.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  project the solution of an eigenvalue-calculation in a space          */
/*		  where the eigenvalues aren't zero.                                                            */
/*		                                                                                                                                */
/*																			*/
/****************************************************************************/

#define PROJECT_CLASS_NAME "project"

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct np_project {
  NP_BASE base;                              /* inherits base class             */

  /* data (optional, necessary for calling the generic execute routine)   */
  VECDATA_DESC *x;                           /* vector                          */

  /* functions */
  INT (*PreProcess)
    (struct np_project *,                    /* pointer to (derived) object     */
    INT *);                                      /* result                          */
  INT (*Project)
    (struct np_project *,                    /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* from level                      */
    VECDATA_DESC *,                              /* vector                          */
    INT *);                                      /* result                          */
  INT (*ProjectionVector)
    (struct np_project *,                    /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* from level                      */
    INT,                                         /* index                           */
    VECDATA_DESC *,                              /* vector                          */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_project *,                    /* pointer to (derived) object     */
    INT *);                                      /* result                          */

  INT dim;                                   /* dimension reduction             */
                                             /* (linear case)                   */

};
typedef struct np_project NP_PROJECT;

typedef INT (*PreProcessProject)(NP_PROJECT *, INT *);
typedef INT (*Project)(NP_PROJECT *, INT, INT, VECDATA_DESC *,INT *);
typedef INT (*ProjectionVector)(NP_PROJECT *, INT, INT, INT, VECDATA_DESC *,INT *);
typedef INT (*PostProcessProject)(NP_PROJECT *, INT *);

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

INT Project_Init (NP_PROJECT *theNP, INT argc, char **argv);
INT Project_Display (NP_PROJECT *theNP);


INT     InitProject (void);

#endif
