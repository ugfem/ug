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


/* RCS_ID
   $Header$
 */

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
#define NL_PARTASS_CLASS_NAME   "nlpass"
#define T_PARTASS_CLASS_NAME    "tpass"

enum PP_ACTIONS {
  PARTASS_UNKNOWN,
  PARTASS_DEFECT  = 1<<0,
  PARTASS_MATRIX  = 1<<1
};

#define PARTASS_DEF_AND_MAT     (PARTASS_DEFECT | PARTASS_MATRIX)

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

/* access macros for NP_T_ASSEMBLE */
#define NPAT_x(p)                               (((NP_T_ASSEMBLE*)(p))->x)
#define NPAT_c(p)                               (((NP_T_ASSEMBLE*)(p))->c)
#define NPAT_b(p)                               (((NP_T_ASSEMBLE*)(p))->b)
#define NPAT_A(p)                               (((NP_T_ASSEMBLE*)(p))->A)

#define NPAT_PRE(p)                             (((NP_T_ASSEMBLE*)(p))->TAssemblePreProcess)
#define NPAT_INITIAL(p)                 (((NP_T_ASSEMBLE*)(p))->TAssembleInitial)
#define NPAT_ASSDEF(p)                  (((NP_T_ASSEMBLE*)(p))->TAssembleDefect)
#define NPAT_ASSSOL(p)                  (((NP_T_ASSEMBLE*)(p))->TAssembleSolution)
#define NPAT_ASSMAT(p)                  (((NP_T_ASSEMBLE*)(p))->TAssembleMatrix)
#define NPAT_POST(p)                    (((NP_T_ASSEMBLE*)(p))->TAssemblePostProcess)
#define NPAT_FINAL(p)                   (((NP_T_ASSEMBLE*)(p))->TAssembleFinal)

/* PARTASS_PARAMS access macros */
#define PP_ACTION(p)                    ((p)->action)
#define PP_SCALE_A(p)                   ((p)->s_a)
#define PP_SCALE_M(p)                   ((p)->s_m)
#define PP_TIME(p)                              ((p)->time)
#define PP_DELTA_T(p)                   ((p)->dt)
#define PP_OLD_DELTA_T(p)               ((p)->dt_old)
#define PP_TIMEDEP(p)                   ((p)->dt!=0.0)
#define PP_ASS_PART(p)                  ((p)->ass_part)
#define PP_SKIP(p)                              ((p)->partskip)
#define PP_CO_SKIP(p)                   ((p)->co_partskip)
#define PP_MD_A(p)                              ((p)->MD_A)
#define PP_MD_A_glob(p)                 ((p)->MD_A_glob)
#define PP_VD_s(p)                              ((p)->VD_s)
#define PP_VD_s_glob(p)                 ((p)->VD_s_glob)
#define PP_VD_s_i(p)                    ((p)->VD_s_i)
#define PP_VD_s_co(p)                   ((p)->VD_s_co)
#define PP_VD_s_ico(p)                  ((p)->VD_s_ico)
#define PP_VD_o(p)                              ((p)->VD_o)
#define PP_VD_o_glob(p)                 ((p)->VD_o_glob)
#define PP_VD_c(p)                              ((p)->VD_c)
#define PP_VD_c_glob(p)                 ((p)->VD_c_glob)
#define PP_VD_r(p)                              ((p)->VD_r)
#define PP_VD_r_glob(p)                 ((p)->VD_r_glob)
#define PP_VD_gridvel(p)                ((p)->VD_gridvel)

/* NP_NL_PARTASS access macros */
#define NPPNL_t(p)                              (((NP_NL_PARTASS*)(p))->t)
#define NPPNL_s(p)                              (((NP_NL_PARTASS*)(p))->s)
#define NPPNL_x(p)                              (((NP_NL_PARTASS*)(p))->x)
#define NPPNL_c(p)                              (((NP_NL_PARTASS*)(p))->c)
#define NPPNL_b(p)                              (((NP_NL_PARTASS*)(p))->b)
#define NPPNL_g(p)                              (((NP_NL_PARTASS*)(p))->g)
#define NPPNL_A(p)                              (((NP_NL_PARTASS*)(p))->A)

#define NPPNL_PRE(p)                    (((NP_NL_PARTASS*)(p))->NLPpreprocess)
#define NPPNL_ASSSOL(p)                 (((NP_NL_PARTASS*)(p))->NLPassembleSolution)
#define NPPNL_ASS(p)                    (((NP_NL_PARTASS*)(p))->NLPassemble)
#define NPPNL_POST(p)                   (((NP_NL_PARTASS*)(p))->NLPpostprocess)

/* NP_T_PARTASS access macros */
#define NPPT_t(p)                               (((NP_T_PARTASS*)(p))->t)
#define NPPT_s(p)                               (((NP_T_PARTASS*)(p))->s)

#define NPPT_INITIAL(p)                 (((NP_T_PARTASS*)(p))->TPassembleInitial)
#define NPPT_PRE(p)                     (((NP_T_PARTASS*)(p))->TPpreprocess)
#define NPPT_ASSSOL(p)                  (((NP_T_PARTASS*)(p))->TPassembleSolution)
#define NPPT_ASS(p)                     (((NP_T_PARTASS*)(p))->TPassemble)
#define NPPT_POST(p)                    (((NP_T_PARTASS*)(p))->TPpostprocess)
#define NPPT_FINAL(p)                   (((NP_T_PARTASS*)(p))->TPassembleFinal)

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
  INT (*NLNAssembleMatrix)
    (struct np_nl_assemble *,                /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    NODE *,                              /* pointer to node                 */
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
typedef INT (*NLNAssembleMatrixProcPtr)                                       \
  (NP_NL_ASSEMBLE *, INT, INT, NODE *, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *,\
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
  INT (*TNAssembleMatrix)                        /* compute linearization (Jacobian)*/
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    NODE *,                              /* pointer to node                 */
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
  INT (*TAssembleFinal)                                  /* call after finishing integration*/
    (struct np_t_assemble *,                 /* pointer to (derived) object     */
    INT,                                         /* from level                      */
    INT,                                         /* to level                        */
    INT *);                                      /* result                          */
};
typedef struct np_t_assemble NP_T_ASSEMBLE;

typedef INT (*TAssemblePreProcessProcPtr)                                     \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *, INT *);
typedef INT (*TAssembleInitialProcPtr)                                       \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, VECDATA_DESC *, INT *);
typedef INT (*TAssembleSolutionProcPtr)                                      \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, VECDATA_DESC *, INT *);
typedef INT (*TAssembleDefectProcPtr)                                        \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*TAssembleMatrixProcPtr)                                        \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*TNAssembleMatrixProcPtr)                                        \
  (NP_T_ASSEMBLE *, INT, INT, NODE *, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *, MATDATA_DESC *, INT *);
typedef INT (*TAssemblePostProcessProcPtr)                                    \
  (NP_T_ASSEMBLE *, INT, INT, DOUBLE, DOUBLE, DOUBLE, VECDATA_DESC *, VECDATA_DESC *, VECDATA_DESC *, INT *);

/****************************************************************************/
/*																			*/
/* nonlinear assemble interface												*/
/*																			*/
/****************************************************************************/

typedef struct {

  INT action;                                                           /* def/mat (enum PP_ACTIONS)	*/
  DOUBLE s_a;                                                           /* scale stiffness mat			*/
  DOUBLE s_m;                                                           /* scale mass matrix			*/
  DOUBLE time;                                                          /* time                                                 */
  DOUBLE dt;                                                            /* time step					*/
  DOUBLE dt_old;                                                        /* last time step				*/
  INT ass_part;                                                         /* assemble part only			*/
  INT partskip[NVECTYPES];                              /* own skip flag positions		*/
  INT co_partskip[NVECTYPES];                           /* co skip flag positions		*/
  MATDATA_DESC   *MD_A;                                         /* stiffness matrix                     */
  MATDATA_DESC   *MD_A_glob;
  VECDATA_DESC   *VD_s;                                         /* solution                                     */
  VECDATA_DESC   *VD_s_glob;
  VECDATA_DESC   *VD_s_i;
  VECDATA_DESC   *VD_s_co;
  VECDATA_DESC   *VD_s_ico;
  VECDATA_DESC   *VD_o;                                         /* last time step sol			*/
  VECDATA_DESC   *VD_o_glob;
  VECDATA_DESC   *VD_c;                                         /* correction					*/
  VECDATA_DESC   *VD_c_glob;
  VECDATA_DESC   *VD_r;                                         /* right hand side				*/
  VECDATA_DESC   *VD_r_glob;
  VECDATA_DESC   *VD_gridvel;                                   /* grid velocity				*/

} PARTASS_PARAMS;

struct np_nl_partass {

  NP_BASE base;                                                 /* inherits base class				*/

  /* data (optional, necessary for calling the generic execute routine)	*/
  VEC_TEMPLATE *t;                                              /* template matching x				*/
  INT s;                                                        /* sub vec for own part                         */
  VECDATA_DESC *x;                                              /* solution                                             */
  VECDATA_DESC *c;                                              /* correction						*/
  VECDATA_DESC *b;                                              /* defect							*/
  VECDATA_DESC *g;                                              /* grid velocity					*/
  MATDATA_DESC *A;                                              /* matrix							*/

  /* functions */
  INT (*NLPpreprocess)(
    struct np_nl_partass *,                             /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*NLPassembleSolution)(
    struct np_nl_partass *,                             /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*NLPassemble)(
    struct np_nl_partass *,                             /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*NLPpostprocess)(
    struct np_nl_partass *,                             /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
};
typedef struct np_nl_partass NP_NL_PARTASS;

/****************************************************************************/
/*																			*/
/* time-dependent part assemble interface									*/
/*																			*/
/****************************************************************************/

struct np_t_partass {

  NP_BASE base;                                                 /* inherits base class				*/

  /* data (optional, necessary for calling the generic execute routine)	*/
  VEC_TEMPLATE *t;                                              /* template matching x				*/
  INT s;                                                        /* sub vec for own part                         */

  /* functions */
  INT (*TPassembleInitial)(
    struct np_t_partass *,                              /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*TPpreprocess)(
    struct np_t_partass *,                              /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*TPassembleSolution)(
    struct np_t_partass *,                              /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*TPassemble)(
    struct np_t_partass *,                              /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*TPpostprocess)(
    struct np_t_partass *,                              /* pointer to (derived) object		*/
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
  INT (*TPassembleFinal)                                /* call after finishing integration*/
    (struct np_t_partass *,                             /* pointer to (derived) object     */
    INT fl,                                                             /* from level						*/
    INT tl,                                                             /* to level                                             */
    PARTASS_PARAMS *pp,                                 /* part assemble parameters             */
    INT *                                                               /* result							*/
    );
};
typedef struct np_t_partass NP_T_PARTASS;

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* generic init/display/execute functions for Assemble num procs */
INT NPAssembleInit                              (NP_BASE *theNP, INT argc , char **argv);
INT NPAssembleDisplay                   (NP_BASE *theNP);
INT NPAssembleExecute                   (NP_BASE *theNP, INT argc , char **argv);

/* generic init/display function for LocalAssemble num procs */
INT NPLocalAssembleInit                 (NP_LOCAL_ASSEMBLE *theNP, INT argc , char **argv);
INT NPLocalAssembleDisplay              (NP_LOCAL_ASSEMBLE *theNP);

/* modification of the matrix for Dirichlet values */
INT NPLocalAssemblePostMatrix   (NP_LOCAL_ASSEMBLE *theNP, INT level,
                                 VECDATA_DESC *x,
                                 VECDATA_DESC *b, MATDATA_DESC *A, INT *result);

/* generic construction of NP_ASSEMBLE from NP_LOCAL_ASSEMBLE */
INT NPLocalAssembleConstruct    (NP_ASSEMBLE *theNP);

/* generic init/display/execute functions for NLAssemble num procs */
INT NPNLAssembleInit                    (NP_BASE *theNP, INT argc , char **argv);
INT NPNLAssembleDisplay                 (NP_BASE *theNP);
INT NPNLAssembleExecute                 (NP_BASE *theNP, INT argc , char **argv);

/* generic init/display/execute functions for time--dependent assembly */
INT NPTAssembleInit                     (NP_BASE *theNP, INT argc , char **argv);
INT NPTAssembleDisplay                  (NP_BASE *theNP);
INT NPTAssembleExecute                  (NP_BASE *theNP, INT argc , char **argv);

void DefaultPartassParams               (PARTASS_PARAMS *pp);
INT SetPartassParams                    (PARTASS_PARAMS *pp, const VEC_TEMPLATE *vt, INT sub,
                                         DOUBLE s_a, DOUBLE s_m, DOUBLE t, DOUBLE dt, DOUBLE dt_old,
                                         VECDATA_DESC *s, VECDATA_DESC *r, VECDATA_DESC *o,
                                         VECDATA_DESC *c, VECDATA_DESC *g, MATDATA_DESC *A);

/* generic init/display/execute functions for NP_NL_PARTASS num procs */
INT NPNLPartAssInit                     (NP_BASE *theNP, INT argc, char **argv);
INT NPNLPartAssDisplay                  (NP_BASE *theNP);
INT NPNLPartAssExecute                  (NP_BASE *theNP, INT argc, char **argv);

/* generic init/display/execute functions for NP_T_PARTASS num procs */
INT NPTPartAssInit                              (NP_BASE *theNP, INT argc, char **argv);
INT NPTPartAssDisplay                   (NP_BASE *theNP);
INT NPTPartAssExecute                   (NP_BASE *theNP, INT argc, char **argv);

const char *pp_action2str               (const PARTASS_PARAMS *pp);

/* init tis file */
INT InitAssemble (void);

#endif
