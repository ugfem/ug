// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  eiter.h                                                                                                       */
/*																			*/
/* Purpose:   definition of the extended iter num proc type			                */
/*																			*/
/* Author:    Klaus Johannsen                                               */
/*            IWR/Technische Simulation                                     */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: klaus.johannsen@iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   22.07.02 begin                                                */
/*																			*/
/* Remarks:                                                                                                                             */
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

#ifndef __EITER__
#define __EITER__

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

#define EITER_CLASS_NAME "ext_iter"

/* access macros */
#define NPEIT_A(p)                      (((NP_EITER*)(p))->A)
#define NPEIT_b(p)                      (((NP_EITER*)(p))->b)
#define NPEIT_c(p)                      (((NP_EITER*)(p))->c)
#define NPEIT_PREP(p)           (((NP_EITER*)(p))->PreProcess)
#define NPEIT_ITER(p)           (((NP_EITER*)(p))->Iter)
#define NPEIT_POST(p)           (((NP_EITER*)(p))->PostProcess)

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct np_eiter {

  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  EVECDATA_DESC *c;                      /* correction                      */
  EVECDATA_DESC *b;                      /* defect                          */
  EMATDATA_DESC *A;                      /* matrix                          */

  /* functions */
  INT (*PreProcess)
    (struct np_eiter *,                      /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector					*/
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    INT *,                                       /* baselevel used by iter          */
    INT *);                                      /* result                          */
  INT (*Iter)
    (struct np_eiter *,                      /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* correction vector               */
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_eiter *,                      /* pointer to (derived) object     */
    INT,                                         /* level                           */
    EVECDATA_DESC *,                             /* solution vector					*/
    EVECDATA_DESC *,                             /* defect vector                   */
    EMATDATA_DESC *,                             /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_eiter NP_EITER;

typedef INT (*PreProcessEIterProcPtr)(NP_EITER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, INT *, INT *);
typedef INT (*EIterProcPtr)(NP_EITER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, INT *);
typedef INT (*PostProcessEIterProcPtr)(NP_EITER *, INT, EVECDATA_DESC *, EVECDATA_DESC *, EMATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for iter num procs */
INT NPEIterInit (NP_EITER *theNP, INT argc , char **argv);

/* generic display function for iter num procs */
INT NPEIterDisplay (NP_EITER *theNP);

/* generic execute function for iter num procs */
INT NPEIterExecute (NP_BASE *theNP, INT argc , char **argv);

/* create standard iter num proc type */
INT InitEIter (void);

#endif
