// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  error.h                                                                                                       */
/*																			*/
/* Purpose:   header for indicator.c										*/
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de					                */
/*																			*/
/* History:   Sep 4, 1996, ug version 3.4                                                               */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __ERROR__
#define __ERROR__

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
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct np_error {

  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  VECDATA_DESC *x;                       /* solution                        */
  VECDATA_DESC *o;                       /* old solution                    */

  /* functions */
  INT (*PreProcess)
    (struct np_error *,                      /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
  INT (*Error)
    (struct np_error *,                      /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*TimeError)
    (struct np_error *,                      /* pointer to (derived) object     */
    INT,                                         /* level                           */
    DOUBLE,                                      /* time                            */
    DOUBLE *,                                    /* time step                       */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* old solution vector             */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_error *,                      /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_error NP_ERROR;

typedef INT (*PreProcessErrorProcPtr)                                       \
  (NP_ERROR *, INT, INT *);
typedef INT (*ErrorProcPtr)                                                 \
  (NP_ERROR *, INT, VECDATA_DESC *, VECDATA_DESC *, INT *);
typedef INT (*TimeErrorProcPtr)                                             \
  (NP_ERROR *, INT, DOUBLE, DOUBLE *, VECDATA_DESC *, VECDATA_DESC *, INT *);
typedef INT (*PostProcessErrorProcPtr)                                      \
  (NP_ERROR *, INT, INT *);

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

INT SurfaceIndicator (MULTIGRID *theMG, VECDATA_DESC *theVD,
                      DOUBLE refine, DOUBLE coarse, INT project,
                      INT from, INT to, INT clear);

/* generic init function for Error num procs */
INT NPErrorInit (NP_ERROR *theNP, INT argc , char **argv);

/* generic display function for Error num procs */
INT NPErrorDisplay (NP_ERROR *theNP);

/* generic execute function for Error num procs */
INT NPErrorExecute (NP_BASE *theNP, INT argc , char **argv);

INT     InitError (void);

#endif
