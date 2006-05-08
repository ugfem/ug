// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      amgtransfer.h                                                 */
/*                                                                          */
/* Purpose:   initialization for algebraic multigrid                        */
/*                                                                          */
/* Author:    Nicolas Neuss                                                 */
/*            Institut fuer Angewandte Mathematik                           */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 294                                       */
/*            69120 Heidelberg                                              */
/*            email: neuss@iwr.uni-heidelberg.de                            */
/*                                                                          */
/* History:   1994-1995 in old ug2.0                                        */
/*            May 1997  in new ug3.7                                        */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __AMGTRANSFER__
#define __AMGTRANSFER__

#include "np.h"
#include "transfer.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*    compile time constants defining static data size (i.e. arrays)        */
/*    other constants                                                       */
/*    macros                                                                */
/*                                                                          */
/****************************************************************************/

#define SELECTION_AMG 1
#define CLUSTER_AMG 2
#define FAMG 3

#define DISPLAY_NP_AMG_STRING "Level  NVectors   NMatrices  NIMats (to finer level)\n"
#define DISPLAY_NP_AMG_FORMAT "%3d   %8d   %8d   %8d\n"

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/** \brief Data type for returning the status of the coarsening procedure */
typedef struct {
  INT nVects;
  INT nMats;
  INT nIMats;
} CLRESULT;

typedef struct {

  /** \brief Error code */
  INT error_code;

  /** \brief Number of created CG levels */
  INT nLevels;

  /** \brief Data for each level */
  CLRESULT clres[MAXLEVEL];
} CRESULT;

typedef INT (*MarkConnectionsProcPtr)(GRID *, MATDATA_DESC *, DOUBLE, INT);
typedef INT (*CoarsenProcPtr)(GRID *);
typedef INT (*SetupIRMatProcPtr)(GRID *, MATDATA_DESC *, MATDATA_DESC *);
typedef INT (*SetupCGMatProcPtr)(GRID *, MATDATA_DESC *, MATDATA_DESC *, INT);

typedef struct
{
  NP_TRANSFER transfer;

  /** \brief Display mode */
  INT display;

  /** \brief Type of AMG                     */
  INT AMGtype;

  /** \brief Mark strong connections         */
  MarkConnectionsProcPtr MarkStrong;

  /** \brief Parameter                       */
  DOUBLE thetaS;

  /** \brief (vector) component to be used   */
  INT compS;

  /** \brief The coarsening routine          */
  CoarsenProcPtr Coarsen;

  /** \brief Setup interpolation/restriction */
  SetupIRMatProcPtr SetupIR;

  /** \brief Setup coarse grid matrix        */
  SetupCGMatProcPtr SetupCG;

  /** \brief Bits 0:symm & 1:R=Inj & 2:P=Inj */
  INT CMtype;

  /** \brief Mark connections to keep        */
  MarkConnectionsProcPtr MarkKeep;

  /** \brief Parameter                       */
  DOUBLE thetaK;

  /** \brief (vector) component to be used   */
  INT compK;

  /** \brief If set, lump to diagonal        */
  INT sparsenFlag;

  /** \brief Ordering of fine grid points    */
  INT reorderFlag;

  /** \brief Transform defect in RS scheme   */
  INT transformdef;

  /** \brief Do fine grid correction step in Reusken/Wagner scheme */
  INT fgcstep;

  /** \brief For fgcstep                     */
  VECDATA_DESC *p;

  /** \brief Stop if vects<vectLimit         */
  INT vectLimit;

  /** \brief Stop if matrices<matLimit       */
  INT matLimit;

  /** \brief Stop if matrices/vects>bandLimit*/
  DOUBLE bandLimit;

  /** \brief Stop if vectReduction<vRedLimit */
  DOUBLE vRedLimit;

  /** \brief Stop if matReduction<mRedLimit  */
  DOUBLE mRedLimit;

  /** \brief Stop if -level>levelLimit       */
  INT levelLimit;

  /** \brief Agglomerate to one processor if level <= aggLimit   */
  INT aggLimit;

  /** \brief Agglomerated bottom level       */
  INT agglevel;

  /** \brief Clear only by npexecute         */
  INT explicitFlag;

  /** \brief No clear in postprocess         */
  INT hold;

  /** \brief Internal: 1 if R=I^t, 0 else    */
  INT symmIR;

} NP_AMG_TRANSFER;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT AMGTransferInit       (NP_BASE *np, INT argc , char **argv);
INT AMGTransferDisplay    (NP_BASE *theNP);
INT AMGTransferPreProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
                           VECDATA_DESC *x, VECDATA_DESC *b,
                           MATDATA_DESC *A, INT *result);
INT AMGTransferExecute    (NP_BASE *theNP, INT argc , char **argv);

INT InitAMGTransfer       (void);
INT AMGTransferConstruct  (NP_BASE *theNP);

END_UGDIM_NAMESPACE

#endif
