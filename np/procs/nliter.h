// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      nliter.c                                                          */
/*                                                                          */
/* Purpose:   nonlinear iteration num procs                                 */
/*                                                                          */
/* Author:    Gabriele Beddies                                              */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   30.07.97 begin, ug version 3.8                                */
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
#include "assemble.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define NL_ITER_CLASS_NAME "nliter"

/* access macros */
#define NPINL_A(p)                      (((NP_NL_ITER*)(p))->A)
#define NPINL_x(p)                      (((NP_NL_ITER*)(p))->x)
#define NPINL_b(p)                      (((NP_NL_ITER*)(p))->b)

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct np_nl_iter {

  NP_BASE base;                              /* inherits base class         */

  /* data (optinal, necessary for calling the generic execute routine)*/
  VECDATA_DESC *b;                       /* defect                      */
  VECDATA_DESC *x;                       /* solution                    */
  MATDATA_DESC *A;                       /* matrix                      */
  NP_NL_ASSEMBLE *Assemble;              /* the assemble numproc		*/

  /* functions */
  INT (*PreProcess)
    (struct np_nl_iter *,            /* pointer to (derived) object     */
    INT,                                     /* level                           */
    VECDATA_DESC *,                          /* defect vector                   */
    VECDATA_DESC *,                          /* solution vector		            */
    MATDATA_DESC *,                          /* matrix                          */
    INT *,                                   /* baselevel used by iter          */
    INT *);                                  /* result                          */
  INT (*NLIter)
    (struct np_nl_iter *,                /* pointer to (derived) object     */
    INT,                                     /* level                           */
    VECDATA_DESC *,                          /* solution vector		            */
    VECDATA_DESC *,                          /* defect vector                   */
    MATDATA_DESC *,                          /* matrix                          */
    NP_NL_ASSEMBLE *,                            /* the assemble numproc			*/
    INT *);                              /* result                          */
  INT (*PostProcess)
    (struct np_nl_iter *,                /* pointer to (derived) object     */
    INT,                                     /* level                           */
    VECDATA_DESC *,                          /* solution vector					*/
    VECDATA_DESC *,                          /* defect vector                   */
    MATDATA_DESC *,                          /* matrix                          */
    INT *);                                  /* result                          */
};
typedef struct np_nl_iter NP_NL_ITER;

typedef INT (*PreProcessNLIterProcPtr) \
  (NP_NL_ITER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *, INT *);
typedef INT (*NLIterProcPtr) \
  (NP_NL_ITER *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, NP_NL_ASSEMBLE *, INT *);
typedef INT (*PostProcessNLIterProcPtr) \
  (NP_NL_ITER *, INT, VECDATA_DESC *, MATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for iter num procs */
INT NPNLIterInit (NP_NL_ITER *theNP, INT argc , char **argv);

/* generic display function for iter num procs */
INT NPNLIterDisplay (NP_NL_ITER *theNP);

/* generic execute function for iter num procs */
INT NPNLIterExecute (NP_BASE *theNP, INT argc , char **argv);

/* create nonlinear iter num proc type */
INT InitNLIter (void);

#endif
