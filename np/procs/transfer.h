// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  transfer.h                                                                                                    */
/*																			*/
/* Purpose:   definition of the transfer num proc type						*/
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

#ifndef __TRANSFER__
#define __TRANSFER__

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

#define TRANSFER_CLASS_NAME "transfer"

/* access macros for NP_TRANSFER */
#define NPTR_x(p)                               (((NP_TRANSFER*)(p))->x)
#define NPTR_c(p)                               (((NP_TRANSFER*)(p))->c)
#define NPTR_b(p)                               (((NP_TRANSFER*)(p))->b)
#define NPTR_A(p)                               (((NP_TRANSFER*)(p))->A)
#define NPTR_DAMP(p)                    (((NP_TRANSFER*)(p))->damp)

#define NPTR_PRE(p)                             (((NP_TRANSFER*)(p))->PreProcess)
#define NPTR_PREPROJ(p)                 (((NP_TRANSFER*)(p))->PreProcessProject)
#define NPTR_INTCOR(p)                  (((NP_TRANSFER*)(p))->InterpolateCorrection)
#define NPTR_RESTRICT(p)                (((NP_TRANSFER*)(p))->RestrictDefect)
#define NPTR_INTNEW(p)                  (((NP_TRANSFER*)(p))->InterpolateNewVectors)
#define NPTR_PROJSOL(p)                 (((NP_TRANSFER*)(p))->ProjectSolution)
#define NPTR_ADPTCOR(p)                 (((NP_TRANSFER*)(p))->AdaptCorrection)
#define NPTR_POST(p)                    (((NP_TRANSFER*)(p))->PostProcess)

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct np_transfer {

  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  VECDATA_DESC *x;                       /* solution                        */
  VECDATA_DESC *c;                       /* correction                      */
  VECDATA_DESC *b;                       /* defect                          */
  MATDATA_DESC *A;                       /* matrix                          */
  VEC_SCALAR damp;                           /* damping factor                  */

  /* functions */
  INT (*PreProcess)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PreProcessProject)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *,                                       /* baselevel                       */
    INT *);                                      /* result                          */
  INT (*InterpolateCorrection)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* destination vector              */
    VECDATA_DESC *,                              /* source vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    VEC_SCALAR,                                  /* damping factor                  */
    INT *);                                      /* result                          */
  INT (*RestrictDefect)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* destination vector              */
    VECDATA_DESC *,                              /* source vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    VEC_SCALAR,                                  /* damping factor                  */
    INT *);                                      /* result                          */
  INT (*InterpolateNewVectors)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*ProjectSolution)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*AdaptCorrection)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* correction vector               */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PostProcessProject)
    (struct np_transfer *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    INT *);                                      /* result                          */
};
typedef struct np_transfer NP_TRANSFER;

typedef INT (*PreProcessTransferProcPtr)                                    \
  (NP_TRANSFER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *,       \
  INT *, INT *);
typedef INT (*InterpolateCorrectionProcPtr)                                 \
  (NP_TRANSFER *, INT, VECDATA_DESC *, VECDATA_DESC *, DOUBLE *, INT *);
typedef INT (*RestrictDefectProcPtr)                                        \
  (NP_TRANSFER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, DOUBLE *, INT *);
typedef INT (*InterpolateSolutionProcPtr)                                   \
  (NP_TRANSFER *, INT, VECDATA_DESC *, INT *);
typedef INT (*ProjectSolutionProcPtr)                                       \
  (NP_TRANSFER *, INT, VECDATA_DESC *, INT *);
typedef INT (*PostProcessTransferProcPtr)                                   \
  (NP_TRANSFER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for transfer num procs */
INT NPTransferInit (NP_TRANSFER *theNP, INT argc , char **argv);

/* generic display function for transfer num procs */
INT NPTransferDisplay (NP_TRANSFER *theNP);

/* generic execute function for transfer num procs */
INT NPTransferExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard transfer num proc type */
INT InitTransfer (void);

#endif
