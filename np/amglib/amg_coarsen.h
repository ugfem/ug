// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amg_coarsen.h													*/
/*																			*/
/* Purpose:   algebraic multigrid coarse grid setup							*/
/*																			*/
/* Author:	  Peter Bastian                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   29 FEB 1996 Begin												*/
/*			  01 OKT 1997 redesign											*/
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

#ifndef __AMG_COARSEN__
#define __AMG_COARSEN__

#include "amg_sp.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define AMG_MAX_COMP                    5               /* max number of components in syste*/

#define AMG_MAX_LEVELS              32          /* max no of levels in hierarchy	*/
#define AMG_MAX_FRONT                   1024    /* max number of front nodes		*/
#define AMG_MAX_CLUSTER                 256             /* max size of a cluster			*/
#define AMG_MAX_STACK                   256             /* size of seed stack				*/
#define AMG_MAX_ROW                             512             /* max number of nonzeros in row	*/

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef struct {                                        /* parameters of coarsening strategy	*/
  int verbose;                                          /* be verbose							*/
  double alpha;                                         /* "relative" strongness                                */
  double beta;                                          /* "absolute" strongness				*/
  int mincluster;                                       /* minimal cluster size					*/
  int maxcluster;                                       /* maximum cluster size					*/
  int maxdistance;                                      /* maximum distance in one cluster		*/
  int maxconnectivity;                          /* limit for connectivity                               */
  int coarsentarget;                                    /* quit coarsening if nodes reached		*/
  int depthtarget;                                      /* create at most so many levels		*/
  double coarsenrate;                                   /* quit if coarsening is too slow		*/
  int major;                                                    /* use major component strategy if >=0  */
} AMG_CoarsenContext;

typedef struct {
  /* graph data structure */
  int n;                                                        /* number of nodes (this is n from A)	*/
  int e;                                                        /* number of edges (nonzeros from A)	*/
  int *ra,*ja;                                          /* connectivity from A					*/
  int *ca;                                                      /* cluster array (fine to coarse map)	*/
  char *na;                                                     /* node array (flags per node)			*/
  char *la;                                                     /* link array							*/
  float *da;                                                    /* damping array with automatic damping */
  int clusters;                                         /* total number of clusters build		*/
  int conclusters;                                      /* number of nonisolated clusters		*/
  int system_as_scalar;                         /* copied from matrix					*/
} AMG_GRAPH;

#define AMG_GRAPH_N(p)                          (p)->n
#define AMG_GRAPH_E(p)                          (p)->e
#define AMG_GRAPH_SAS(p)                        (p)->system_as_scalar
#define AMG_GRAPH_RA(p)                         (p)->ra
#define AMG_GRAPH_JA(p)                         (p)->ja
#define AMG_GRAPH_CA(p)                         (p)->ca
#define AMG_GRAPH_NA(p)                         (p)->na
#define AMG_GRAPH_LA(p)                         (p)->la
#define AMG_GRAPH_DA(p)                         (p)->da


/****************************************************************************/
/*																			*/
/* functions																*/
/*																			*/
/****************************************************************************/

int AMG_BuildHierarchy (AMG_CoarsenContext *cc, AMG_MATRIX *A,
                        AMG_MATRIX *H[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS]);

#endif
