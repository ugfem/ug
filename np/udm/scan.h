// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  scan.h                                                                                                        */
/*																			*/
/* Purpose:   header file for scanning routines for npinit calls            */
/*																			*/
/* Author:	  Christian Wieners                                                                     */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   December 8, 1996                                                                  */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __SCAN__
#define __SCAN__

#include "udm.h"
#include "numproc.h"

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
/* general numerics defines													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* macros concerned with data descriptors and symbols						*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* macros concerned with solving											*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* structures concerned with symbolic user data management					*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* reading VECDESC and MATDESC                                              */

VECDATA_DESC *ReadArgvVecDesc      (MULTIGRID *theMG, char *name,
                                    INT argc, char **argv);
MATDATA_DESC *ReadArgvMatDesc      (MULTIGRID *theMG, char *name,
                                    INT argc, char **argv);

NP_BASE      *ReadArgvNumProc      (MULTIGRID *theMG, char *name, char *class,
                                    INT argc, char **argv);

/* for reading dampinf factors etc. */
INT ReadVecTypeINTs             (char *str, INT n, INT nINT[MAXVECTORS], INT theINTs[][MAXVECTORS]);
INT ReadVecTypeDOUBLEs  (char *str, INT n, INT nDOUBLE[MAXVECTORS], DOUBLE theDOUBLEs[][MAXVECTORS]);
INT ReadVecTypeOrder    (char *str, INT n, INT MaxPerType, INT *nOrder, INT theOrder[]);

/* tools for VEC_SCALAR                                                     */
INT sc_read          (VEC_SCALAR x, const VECDATA_DESC *theVD, const char *name,
                      INT argc, char **argv);
INT sc_disp      (VEC_SCALAR x, const VECDATA_DESC *theVD, const char *name);
INT sc_cmp           (VEC_SCALAR x, const VEC_SCALAR y, const VECDATA_DESC *theVD);
INT sc_mul           (VEC_SCALAR x, const VEC_SCALAR y, const VEC_SCALAR z,
                      const VECDATA_DESC *theVD);
INT sc_mul_check (VEC_SCALAR x, const VEC_SCALAR y, const VEC_SCALAR z,
                  const VECDATA_DESC *theVD);

/* scanning argument lists                                                  */
INT ReadArgvDOUBLE (char *name, DOUBLE *a, INT argc, char **argv);
INT ReadArgvINT (char *name, INT *j, INT argc, char **argv);
INT ReadArgvChar (char *name, char *buffer, INT argc, char **argv);
INT ReadArgvDisplay (INT argc, char **argv);
INT ReadArgvOption (char *name, INT argc, char **argv);
INT ReadArgvPosition (char *name, INT argc, char **argv, DOUBLE *pos);

#endif
