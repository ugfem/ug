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

#define ASSEMBLE_CLASS_NAME     "assemble"
#define NL_ASSEMBLE_CLASS_NAME  "nlass"

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* linear assemble interface												*/
/*																			*/
/****************************************************************************/

struct np_assemble {

  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  VECDATA_DESC *x;                       /* solution                        */
  VECDATA_DESC *b;                       /* defect                          */
  MATDATA_DESC *A;                       /* matrix                          */

  /* functions */
  INT (*PreProcess)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* rhs vector                          */
    MATDATA_DESC *,                              /* matrix                                                      */
    INT *);                                      /* result                          */
  INT (*Assemble)
    (struct np_assemble *,                   /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                                          /* current solution	(initial)	*/
    VECDATA_DESC *,                                          /* right hand side                         */
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
typedef INT (*PostProcessAssembleProcPtr)                                   \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* local assemble interface													*/
/*																			*/
/****************************************************************************/

struct np_local_assemble {

  NP_ASSEMBLE assemble;                      /* inherits assemble class         */

  /* data */
  INT galerkin;                              /* Galerkin assembling             */

  /* functions */
  INT (*PreProcess)
    (struct np_local_assemble *,             /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    DOUBLE **,                           /* local solution                  */
    DOUBLE **,                           /* local defect                    */
    DOUBLE **,                           /* local matrix                    */
    INT **,                              /* local vecskip                   */
    INT *);                                      /* result                          */
  INT (*AssembleLocal)
    (ELEMENT *,                              /* pointer to an element           */
    INT *);                                      /* result                          */
  INT (*AssembleLocalDefect)
    (ELEMENT *,                              /* pointer to an element           */
    INT *);                                      /* result                          */
  INT (*AssembleLocalMatrix)
    (ELEMENT *,                              /* pointer to an element           */
    INT *);                                      /* result                          */
  INT (*PostMatrix)
    (struct np_local_assemble *,             /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_local_assemble *,             /* pointer to (derived) object     */
    INT,                                         /* level                           */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_local_assemble NP_LOCAL_ASSEMBLE;

typedef INT (*PreProcessLocalAssembleProcPtr)                               \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *,       \
  DOUBLE **, DOUBLE **, DOUBLE **, INT **, INT *);
typedef INT (*AssembleLocalProcPtr)                                         \
  (ELEMENT *, INT *);
typedef INT (*AssembleLocalDefectProcPtr)                                   \
  (ELEMENT *, INT *);
typedef INT (*AssembleLocalMatrixProcPtr)                                   \
  (ELEMENT *, INT *);
typedef INT (*PostMatrixLocalAssembleProcPtr)                               \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*PostProcessLocalAssembleProcPtr)                              \
  (NP_ASSEMBLE *, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);


/****************************************************************************/
/*																			*/
/* nonlinear assemble interface												*/
/*																			*/
/****************************************************************************/

struct np_nl_assemble {

  NP_BASE base;                              /* inherits base class             */

  /* data (optinal, necessary for calling the generic execute routine)    */
  VECDATA_DESC *x;                       /* solution                        */
  VECDATA_DESC *c;                       /* correction                      */
  VECDATA_DESC *b;                       /* defect                          */
  MATDATA_DESC *A;                       /* matrix                          */

  /* functions */
  INT (*PreProcess)
    (struct np_nl_assemble *,                /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*NLAssembleSolution)
    (struct np_nl_assemble *,                /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    VECDATA_DESC *,                              /* solution vector                 */
    INT *);                                      /* result                          */
  INT (*NLAssembleDefect)
    (struct np_nl_assemble *,                /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*NLAssembleMatrix)
    (struct np_nl_assemble *,                /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    VECDATA_DESC *,                                          /* current solution	(initial)	*/
    VECDATA_DESC *,                                          /* defect for current solution     */
    VECDATA_DESC *,                                          /* correction to be computed               */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*PostProcess)
    (struct np_nl_assemble *,                /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    VECDATA_DESC *,                              /* solution vector                 */
    VECDATA_DESC *,                              /* defect vector                   */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
};
typedef struct np_nl_assemble NP_NL_ASSEMBLE;

typedef INT (*PreProcessNLAssembleProcPtr)                                   \
  (NP_NL_ASSEMBLE *, INT, INT, VECDATA_DESC *, INT *);
typedef INT (*NLAssembleSolutionProcPtr)                                     \
  (NP_NL_ASSEMBLE *, INT, INT, VECDATA_DESC *, INT *);
typedef INT (*NLAssembleDefectProcPtr)                                       \
  (NP_NL_ASSEMBLE *, INT, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*NLAssembleMatrixProcPtr)                                       \
  (NP_NL_ASSEMBLE *, INT, INT, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *,\
  MATDATA_DESC *, INT *);
typedef INT (*PostProcessNLAssembleProcPtr)                                  \
  (NP_NL_ASSEMBLE *, INT, INT, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init function for Assemble num procs */
INT NPAssembleInit (NP_BASE *theNP, INT argc , char **argv);

/* generic display function for Assemble num procs */
INT NPAssembleDisplay (NP_BASE *theNP);

/* generic execute function for Assemble num procs */
INT NPAssembleExecute (NP_BASE *theNP, INT argc , char **argv);

/* generic init function for LocalAssemble num procs */
INT NPLocalAssembleInit (NP_LOCAL_ASSEMBLE *theNP, INT argc , char **argv);

/* generic display function for LocalAssemble num procs */
INT NPLocalAssembleDisplay (NP_LOCAL_ASSEMBLE *theNP);

/* modification of the matrix for Dirichlet values */
INT NPLocalAssemblePostMatrix (NP_LOCAL_ASSEMBLE *theNP, INT level,
                               VECDATA_DESC *x,
                               VECDATA_DESC *b, MATDATA_DESC *A, INT *result);

/* generic construction of NP_ASSEMBLE from NP_LOCAL_ASSEMBLE */
INT NPLocalAssembleConstruct (NP_ASSEMBLE *theNP);

/* generic init function for NLAssemble num procs */
INT NPNLAssembleInit (NP_BASE *theNP, INT argc , char **argv);

/* generic display function for NLAssemble num procs */
INT NPNLAssembleDisplay (NP_BASE *theNP);

/* generic execute function for NLAssemble num procs */
INT NPNLAssembleExecute (NP_BASE *theNP, INT argc , char **argv);

#endif
