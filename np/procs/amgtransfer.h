// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  amgtransfer.h                                                                                                 */
/*																			*/
/* Purpose:   initialization for algebraic multigrid	                    */
/*                                                                          */
/* Author:	  Nicolas Neuss                                                                                     */
/*			  Institut fuer Angewandte Mathematik                           */
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 294										*/
/*			  69120 Heidelberg												*/
/*			  email: neuss@iwr.uni-heidelberg.de			                        */
/*																			*/
/* History:   1994-1995 in old ug2.0							            */
/*            May 1997  in new ug3.7                                        */
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

#ifndef __AMGTRANSFER__
#define __AMGTRANSFER__

#include "np.h"
#include "transfer.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define SELECTION_AMG 1
#define CLUSTER_AMG 2

#define DISPLAY_NP_AMG_STRING "Level  NVectors   NMatrices  NIMats (to finer level)\n"
#define DISPLAY_NP_AMG_FORMAT "%3d   %8d   %8d   %8d\n"

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/* a data type for returning the status of the coarsening procedure         */
typedef struct {
  INT nVects;
  INT nMats;
  INT nIMats;
} CLRESULT;

typedef struct {
  INT error_code;                           /* error code                       */
  INT nLevels;                              /* number of created CG levels      */
  CLRESULT clres[MAXLEVEL];                 /* data for each level              */
} CRESULT;

typedef INT (*MarkConnectionsProcPtr)(GRID *, MATDATA_DESC *, DOUBLE, INT);
typedef INT (*CoarsenProcPtr)(GRID *);
typedef INT (*SetupIRMatProcPtr)(GRID *, MATDATA_DESC *, MATDATA_DESC *);
typedef INT (*SetupCGMatProcPtr)(GRID *, MATDATA_DESC *, MATDATA_DESC *, INT);

typedef struct
{
  NP_TRANSFER transfer;
  INT display;                               /* display modus                   */

  INT AMGtype;                               /* type of AMG                     */
  MarkConnectionsProcPtr MarkStrong;         /* mark strong connections         */
  DOUBLE thetaS;                             /* parameter                       */
  INT compS;                                 /* (vector) component to be used   */

  CoarsenProcPtr Coarsen;                    /* the coarsening routine          */

  SetupIRMatProcPtr SetupIR;                 /* setup interpolation/restriction */

  SetupCGMatProcPtr SetupCG;                 /* setup coarse grid matrix        */
  INT CMtype;                                /* Bits 0:symm & 1:R=Inj & 2:P=Inj */

  MarkConnectionsProcPtr MarkKeep;           /* mark connections to keep        */
  DOUBLE thetaK;                             /* parameter                       */
  INT compK;                                 /* (vector) component to be used   */
  INT sparsenFlag;                           /* if set, lump to diagonal        */

  INT reorderFlag;                           /* ordering of fine grid points    */

  INT fgcstep;                                   /* do fine grid correction step in */
  /* Reusken/Wagner scheme           */
  VECDATA_DESC *p;                           /* for fgcstep                     */

  INT vectLimit;                             /* stop if vects<vectLimit         */
  INT matLimit;                              /* stop if matrices<matLimit       */
  DOUBLE bandLimit;                          /* stop if matrices/vects>bandLimit*/
  DOUBLE vRedLimit;                          /* stop if vectReduction<vRedLimit */
  DOUBLE mRedLimit;                          /* stop if matReduction<mRedLimit  */
  INT levelLimit;                            /* stop if -level>levelLimit       */
  INT aggLimit;                              /* agglomerate to one processor    */
  /* if level <= aggLimit            */
  INT agglevel;                          /* agglomerated bottom level       */

  INT explicitFlag;                          /* clear only by npexecute         */
  INT hold;                                  /* no clear in postprocess         */

  INT symmIR;                                /* internal: 1 if R=I^t, 0 else    */

} NP_AMG_TRANSFER;

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

INT AMGTransferInit       (NP_BASE *np, INT argc , char **argv);
INT AMGTransferDisplay    (NP_BASE *theNP);
INT AMGTransferPreProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *result);
INT AMGTransferExecute    (NP_BASE *theNP, INT argc , char **argv);

INT InitAMGTransfer       (void);

#endif
