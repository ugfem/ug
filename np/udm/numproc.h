// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  numproc.h                                                                                                     */
/*																			*/
/* Purpose:   definition of the basic num proc type                                 */
/*																			*/
/* Author:	  Peter Bastian                                                                                 */
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

#ifndef __NUMPROC__
#define __NUMPROC__

#include "ugenv.h"
#include "gm.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* status for num procs  */
#define NP_NOT_INIT                                     0
#define NP_NOT_ACTIVE                           1
#define NP_ACTIVE                                       2
#define NP_EXECUTABLE                           3
#define NP_PDT_SIZE                32

/* macros for NP_BASE access */
#define NP_MG(p)                                (((NP_BASE*)(p))->mg)
#define NP_STATUS(p)                    (((NP_BASE*)(p))->status)
#define NP_INIT(p)                              (((NP_BASE*)(p))->Init)
#define NP_DISPLAY(p)                   (((NP_BASE*)(p))->Display)
#define NP_EXECUTE(p)                   (((NP_BASE*)(p))->Execute)

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

struct np_base {

  /* data */
  ENVVAR v;                                                     /* is an environment variable		*/
  MULTIGRID *mg;                                                /* associated multigrid				*/
  INT status;                                                   /* has a status, NO type and size...*/

  /* functions */
  INT (*Init)(struct np_base *, INT, char **);    /* initializing routine   */
  INT (*Display)(struct np_base *);                      /* Display routine */
  INT (*Execute)(struct np_base *, INT, char **);        /* Execute routine */
};
typedef struct np_base NP_BASE;

typedef INT (*InitNumProcProcPtr)(NP_BASE *, INT, char **);
typedef INT (*DisplayNumProcProcPtr)(NP_BASE *);
typedef INT (*ExecuteNumProcProcPtr)(NP_BASE *, INT, char **);

typedef INT (*ConstructorProcPtr)(NP_BASE *);

typedef struct
{
  ENVVAR v;                                                     /* class name							*/
  INT size;                             /* size of the object                   */
  ConstructorProcPtr Construct;         /* construct one object of this class   */
} NP_CONSTRUCTOR;

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

INT CreateClass (char *classname, INT size, ConstructorProcPtr Construct);
NP_CONSTRUCTOR *GetConstructor (const char *classname);
INT CreateObject (MULTIGRID *theMG, char *objectname, char *classname);

NP_BASE *GetNumProcByName (MULTIGRID *theMG, char *objectname, char *classname);

INT InitNumProcManager (void);

#endif
