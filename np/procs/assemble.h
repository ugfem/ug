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
#define T_ASSEMBLE_CLASS_NAME   "tass"

/* access macros for NP_ASSEMBLE */
#define NPASS_x(p)                              (((NP_ASSEMBLE*)(p))->x)
#define NPASS_c(p)                              (((NP_ASSEMBLE*)(p))->c)
#define NPASS_b(p)                              (((NP_ASSEMBLE*)(p))->b)
#define NPASS_A(p)                              (((NP_ASSEMBLE*)(p))->A)

#define NPASS_PRE(p)                    (((NP_ASSEMBLE*)(p))->PreProcess)
#define NPASS_ASSSOL(p)                 (((NP_ASSEMBLE*)(p))->AssembleSolution)
#define NPASS_ASSDEF(p)                 (((NP_ASSEMBLE*)(p))->AssembleDefect)
#define NPASS_ASSMAT(p)                 (((NP_ASSEMBLE*)(p))->AssembleMatrix)
#define NPASS_ASS(p)                    (((NP_ASSEMBLE*)(p))->Assemble)
#define NPASS_POST(p)                   (((NP_ASSEMBLE*)(p))->PostProcess)

/* access macros for NP_LOCAL_ASSEMBLE */
#define NPLOC_GALERKIN(p)               (((NP_LOCAL_ASSEMBLE*)(p))->galerkin)
#define NPLOC_PRE(p)                    (((NP_LOCAL_ASSEMBLE*)(p))->PreProcess)
#define NPLOC_ASS(p)                    (((NP_LOCAL_ASSEMBLE*)(p))->AssembleLocal)
#define NPLOC_ASSDEF(p)                 (((NP_LOCAL_ASSEMBLE*)(p))->AssembleLocalDefect)
#define NPLOC_ASSMAT(p)                 (((NP_LOCAL_ASSEMBLE*)(p))->AssembleLocalMatrix)
#define NPLOC_POSTMAT(p)                (((NP_LOCAL_ASSEMBLE*)(p))->PostMatrix)
#define NPLOC_POST(p)                   (((NP_LOCAL_ASSEMBLE*)(p))->PostProcess)

/* access macros for NP_NL_ASSEMBLE */
#define NPANL_x(p)                              (((NP_NL_ASSEMBLE*)(p))->x)
#define NPANL_c(p)                              (((NP_NL_ASSEMBLE*)(p))->c)
#define NPANL_b(p)                              (((NP_NL_ASSEMBLE*)(p))->b)
#define NPANL_A(p)                              (((NP_NL_ASSEMBLE*)(p))->A)

#define NPANL_PRE(p)                    (((NP_NL_ASSEMBLE*)(p))->PreProcess)
#define NPANL_ASSDEF(p)                 (((NP_NL_ASSEMBLE*)(p))->NLAssembleDefect)
#define NPANL_ASSSOL(p)                 (((NP_NL_ASSEMBLE*)(p))->NLAssembleSolution)
#define NPANL_ASSMAT(p)                 (((NP_NL_ASSEMBLE*)(p))->NLAssembleMatrix)
#define NPANL_POST(p)                   (((NP_NL_ASSEMBLE*)(p))->PostProcess)

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
/* time-dependent assemble interface										*/
/*																			*/
/****************************************************************************/

struct np_t_assemble {

  NP_BASE base;                              /* inherits base class             */

  /* functions */
  INT (*TAssemblePreProcess)                     /* call at begin of timestep	        */
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    DOUBLE,                                                              /* time t_k+1						*/
    DOUBLE,                                                              /* time t_k						*/
    DOUBLE,                                                              /* time t_k-1						*/
    VECDATA_DESC *,                              /* (unknown) solution at t_k+1         */
    VECDATA_DESC *,                              /* solution vector at t_k          */
    VECDATA_DESC *,                              /* solution vector at t_k-1        */
    INT *);                                      /* result                          */
  INT (*TAssembleInitial)                        /* set initial values				*/
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    DOUBLE,                                                              /* time value t					*/
    VECDATA_DESC *,                              /* solution vector at time t       */
    INT *);                                      /* result                          */
  INT (*TAssembleSolution)               /* set dirichlet conditions in sol.*/
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    DOUBLE,                                                              /* time value t					*/
    VECDATA_DESC *,                              /* solution vector at time t       */
    INT *);                                      /* result                          */
  INT (*TAssembleDefect)                     /* accumulate to defect vector		*/
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    DOUBLE,                                                              /* time value t					*/
    DOUBLE,                                                              /* scaling for m-term: s_m			*/
    DOUBLE,                                                              /* scaling for a-term: s_a			*/
    VECDATA_DESC *,                              /* solution vector y               */
    VECDATA_DESC *,                              /* accumulate s_m*m(t,y)+s_a*a(t,y)*/
    MATDATA_DESC *,                              /* matrix may be handy for Picard  */
    INT *);                                      /* result                          */
  INT (*TAssembleMatrix)                         /* compute linearization (Jacobian)*/
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    DOUBLE,                                                              /* time value t					*/
    DOUBLE,                                                              /* scaling for a-term: s_a	(s_m=1!)*/
    VECDATA_DESC *,                                          /* current sol (linearization pt)  */
    VECDATA_DESC *,                                          /* defect for current solution     */
    VECDATA_DESC *,                                          /* correction to be computed               */
    MATDATA_DESC *,                              /* matrix                          */
    INT *);                                      /* result                          */
  INT (*TAssemblePostProcess)                /* call after solution t_k+1 known */
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    DOUBLE,                                                              /* time t_k+1						*/
    DOUBLE,                                                              /* time t_k						*/
    DOUBLE,                                                              /* time t_k-1						*/
    VECDATA_DESC *,                              /* solution t_k+1 (just computed!)	*/
    VECDATA_DESC *,                              /* solution vector at t_k          */
    VECDATA_DESC *,                              /* solution vector at t_k-1        */
    INT *);                                      /* result                          */
};
typedef struct np_t_assemble NP_T_ASSEMBLE;

typedef INT (*TAssemblePreProcessProcPtr)                                     \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *, INT *);
typedef INT (*TAssembleInitialProcPtr)                                       \
  (NP_T_ASSEMBLE *, INT, INT, VECDATA_DESC *, INT *);
typedef INT (*TAssembleSolutionProcPtr)                                      \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, VECDATA_DESC *, INT *);
typedef INT (*TAssembleDefectProcPtr)                                        \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*TAssembleMatrixProcPtr)                                        \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*TAssemblePostProcessProcPtr)                                    \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *, INT *);

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

/* default functions for time--dependent assembly */
INT NPTAssembleInit     (NP_BASE *theNP, INT argc , char **argv);
INT NPTAssembleDisplay  (NP_BASE *theNP);
INT NPTAssembleExecute  (NP_BASE *theNP, INT argc , char **argv);

#endif
