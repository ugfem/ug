// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  npscan.h                                                                                                      */
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
/*			  np part of former np/udm/scan.h, 15.5.97						*/
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

#ifndef __NPSCAN__
#define __NPSCAN__

#include "scan.h"
#include "udm.h"
#include "numproc.h"
#include "formats.h"

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

INT ReadArgvPosition                                    (const char *name, INT argc, char **argv, DOUBLE *pos);

VECDATA_DESC *ReadArgvVecDescX                  (MULTIGRID *theMG, const char *name,
                                                 INT argc, char **argv, INT CreateIfNonExistent);
VEC_TEMPLATE *ReadArgvVecTemplate               (const FORMAT *fmt, const char *name,
                                                 INT argc, char **argv);
VEC_TEMPLATE *ReadArgvVecTemplateSub    (const FORMAT *fmt, const char *name,
                                         INT argc, char **argv, INT *sub);
MAT_TEMPLATE *ReadArgvMatTemplateSub    (const FORMAT *fmt, const char *name,
                                         INT argc, char **argv, INT *sub);
MATDATA_DESC *ReadArgvMatDescX                  (MULTIGRID *theMG, const char *name,
                                                 INT argc, char **argv, INT CreateIfNonExistent);

NP_BASE      *ReadArgvNumProc                   (MULTIGRID *theMG, const char *name, const char *cls,
                                                 INT argc, char **argv);

#define ReadArgvVecDesc(mg,n,ac,av)             ReadArgvVecDescX(mg,n,ac,av,YES)
#define ReadArgvMatDesc(mg,n,ac,av)             ReadArgvMatDescX(mg,n,ac,av,YES)

/* for reading damping factors etc. */
INT ReadVecTypeINTs             (const FORMAT *fmt, char *str, INT n, INT nINT[MAXVECTORS], INT theINTs[][MAXVECTORS]);
INT ReadVecTypeDOUBLEs  (const FORMAT *fmt, char *str, INT n, INT nDOUBLE[MAXVECTORS], DOUBLE theDOUBLEs[][MAXVECTORS]);
INT ReadVecTypeOrder    (const FORMAT *fmt, char *str, INT n, INT MaxPerType, INT *nOrder, INT theOrder[]);
INT ReadVecTypeNUMPROCs (const MULTIGRID *theMG, char *str, char *class_name, INT n, INT nNUMPROC[MAXVECTORS], NP_BASE *theNUMPROCs[][MAXVECTORS]);

/* tools for VEC_SCALAR */
INT sc_read          (VEC_SCALAR x, const FORMAT *fmt, const VECDATA_DESC *theVD, const char *name, INT argc, char **argv);
INT sc_disp      (VEC_SCALAR x, const VECDATA_DESC *theVD, const char *name);
INT sc_cmp           (VEC_SCALAR x, const VEC_SCALAR y, const VECDATA_DESC *theVD);
INT sc_eq            (VEC_SCALAR x, const VEC_SCALAR y, DOUBLE ac, const VECDATA_DESC *theVD);
INT sc_mul           (VEC_SCALAR x, const VEC_SCALAR y, const VEC_SCALAR z,
                      const VECDATA_DESC *theVD);
INT sc_mul_check (VEC_SCALAR x, const VEC_SCALAR y, const VEC_SCALAR z,
                  const VECDATA_DESC *theVD);

/* scanning argument lists                                                  */
INT ReadArgvDisplay (INT argc, char **argv);

#endif
