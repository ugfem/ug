// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  assemble.h                                                                                                    */
/*																			*/
/* Purpose:   definition of the assemble num proc type                                  */
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

#ifndef __ASSEMBLE__
#define __ASSEMBLE__

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

struct np_assemble {

  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  VECDATA_DESC *x;                       /* solution                        */
  VECDATA_DESC *c;                       /* correction                      */
  VECDATA_DESC *b;                       /* defect                          */
  MATDATA_DESC *A;                       /* matrix                          */

  /* functions */
  INT (*PreProcess)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*Assemble)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*AssembleSolution)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*AssembleDefect)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*AssembleMatrix)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                                          /* current solution	(initial)	*/
    VECDATA_DESC *,                                          /* defect for current solution     */
    VECDATA_DESC *,                                          /* correction to be computed               */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_assemble NP_ASSEMBLE;

typedef INT (*PreProcessAssembleProcPtr)                                    \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*AssembleProcPtr)                                              \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*AssembleSolutionProcPtr)                                      \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, INT *);
typedef INT (*AssembleDefectProcPtr)                                        \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*AssembleMatrixProcPtr)                                        \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *,       \
  MATDATA_DESC *, INT *);
typedef INT (*PostProcessAssembleProcPtr)                                   \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for Assemble num procs */
INT NPAssembleInit (NP_ASSEMBLE *theNP, INT argc , char **argv);

/* generic display function for Assemble num procs */
INT NPAssembleDisplay (NP_ASSEMBLE *theNP);

/* generic execute function for Assemble num procs */
INT NPAssembleExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard Assemble num proc type
   INT InitAssemble (void); */

#endif
