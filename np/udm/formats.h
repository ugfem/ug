// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      formats.h                                                     */
/*                                                                          */
/* Purpose:   header file for format definition				                */
/*                                                                          */
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: henrik@ica3.uni-stuttgart.de							*/
/*			  fon: 0049-(0)711-685-7007										*/
/*			  fax: 0049-(0)711-685-7000										*/
/*																			*/
/* History:   27.03.95 begin, ug version 3.0								*/
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

#ifndef __FORMATS__
#define __FORMATS__

#include "gm.h"
#include "udm.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* limits for XDATA_DESC handling */
#define MAX_SUB                         10

#define V_COMP_NAMES            (MAX_VEC_COMP*NVECTYPES)
#define M_COMP_NAMES            (2*V_COMP_NAMES*V_COMP_NAMES)

/* macros for SUBVEC */
#define SUBV_NAME(s)            ((s)->Name)
#define SUBV_NCOMPS(s)          ((s)->Comp)
#define SUBV_NCOMP(s,tp)        ((s)->Comp[tp])
#define SUBV_COMP(s,tp,i)       ((s)->Comps[tp][i])

/* macros for SUBMAT */
#define SUBM_NAME(s)            ((s)->Name)
#define SUBM_RCOMPS(s)          ((s)->RComp)
#define SUBM_CCOMPS(s)          ((s)->CComp)
#define SUBM_RCOMP(s,tp)        ((s)->RComp[tp])
#define SUBM_CCOMP(s,tp)        ((s)->CComp[tp])
#define SUBM_COMP(s,tp,i)       ((s)->Comps[tp][i])

/* macros for VEC_TEMPLATE */
#define VT_COMPS(vt)            ((vt)->Comp)
#define VT_COMP(vt,tp)          ((vt)->Comp[tp])
#define VT_COMPNAMES(vt)        ((vt)->CompNames)
#define VT_COMPNAME(vt,i)       ((vt)->CompNames[i])
#define VT_NID(vt)                      ((vt)->nId)
#define VT_IDENT_PTR(vt)        ((vt)->Ident)
#define VT_IDENT(vt,i)          ((vt)->Ident[i])
#define VT_SUB(vt,i)            ((vt)->SubVec[i])
#define VT_NSUB(vt)                     ((vt)->nsub)

/* macros for MAT_TEMPLATE */
#define MT_RCOMPS(mt)           ((mt)->RComp)
#define MT_RCOMP(mt,tp)         ((mt)->RComp[tp])
#define MT_CCOMPS(mt)           ((mt)->CComp)
#define MT_CCOMP(mt,tp)         ((mt)->CComp[tp])
#define MT_COMPNAMES(mt)        ((mt)->CompNames)
#define MT_COMPNAME(mt,i)       ((mt)->CompNames[i])
#define MT_SUB(mt,i)            ((mt)->SubMat[i])
#define MT_NSUB(mt)                     ((mt)->nsub)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/* sub vector of vector template (components form a subset of the template)	*/
typedef struct {

  char Name[NAMESIZE];                                          /* prefix for sub vector name	*/
  SHORT Comp[NVECTYPES];                                        /* number of comps per type		*/
  SHORT Comps[NVECTYPES][MAX_VEC_COMP];         /* subsequent comps rel to tplt	*/

} SUBVEC;

/* vector template specifying number of comps per type and comp names */
typedef struct {

  ENVITEM v;                                                                    /* environment item				*/

  SHORT Comp[NVECTYPES];                                        /* number of comps per type		*/
  char CompNames[V_COMP_NAMES];                         /* comp names (one char each)	*/

  SHORT nId;                                                            /* number of comps after ident	*/
  SHORT Ident[V_COMP_NAMES];                            /* identification table			*/

  SHORT nsub;                                                           /* number of sub vectors		*/
  SUBVEC  *SubVec[MAX_SUB];                                     /* pointers to sub vectors		*/

} VEC_TEMPLATE;

/* sub matrix of matrix template (components form a subset of the template)	*/
typedef struct {

  char Name[NAMESIZE];                                          /* prefix for sub matrix name	*/
  SHORT RComp[NMATTYPES];                                       /* number of row comps per type	*/
  SHORT CComp[NMATTYPES];                                       /* number of col comps per type	*/
  SHORT Comps[NMATTYPES][MAX_MAT_COMP];         /* subsequent comps rel to tplt	*/

} SUBMAT;

/* matrix template specifying number of row/col comps per type and comp names */
typedef struct {

  ENVITEM v;                                                                    /* environment item				*/

  SHORT RComp[NMATTYPES];                                       /* number of comps per type		*/
  SHORT CComp[NMATTYPES];                                       /* number of col comps per type	*/
  char CompNames[M_COMP_NAMES];                         /* comp names (two chars each)	*/

  SHORT nsub;                                                           /* number of sub matrices		*/
  SUBMAT  *SubMat[MAX_SUB];                                     /* pointers to sub matrices		*/

} MAT_TEMPLATE;

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

INT DisplayPrintingFormat                               (void);
INT SetPrintingFormatCmd                                (const MULTIGRID *mg, INT argc, char **argv);
INT ResetPrintingFormat                                 (void);

VEC_TEMPLATE *GetVectorTemplate                 (const FORMAT *theFmt, const char *template);
MAT_TEMPLATE *GetMatrixTemplate                 (const FORMAT *theFmt, const char *template);

VECDATA_DESC *CreateVecDescOfTemplate   (MULTIGRID *theMG,
                                         const char *name, const char *template);
MATDATA_DESC *CreateMatDescOfTemplate   (MULTIGRID *theMG,
                                         const char *name, const char *template);

INT VDmatchesVT                                                 (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt);
INT MDmatchesMT                                                 (const MATDATA_DESC *md, const MAT_TEMPLATE *mt);
INT MDmatchesVT                                                 (const MATDATA_DESC *md, const VEC_TEMPLATE *vt);
INT MDmatchesVTxVT                                              (const MATDATA_DESC *md, const VEC_TEMPLATE *rvt, const VEC_TEMPLATE *cvt);

INT VDsubDescFromVT                                             (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt, INT sub, VECDATA_DESC **subvd);
INT VDsubDescFromVS                                             (const VECDATA_DESC *vd, const SUBVEC *subv, VECDATA_DESC **subvd);
INT MDsubDescFromVT                                             (const MATDATA_DESC *md, const VEC_TEMPLATE *vt, INT sub, MATDATA_DESC **submd);
INT MDsubDescFromVTxVT                                  (const MATDATA_DESC *md, const VEC_TEMPLATE *rvt, INT rsub,
                                                         const VEC_TEMPLATE *cvt, INT csub,
                                                         MATDATA_DESC **submd);
INT MDsubDescFromMT                                             (const MATDATA_DESC *md, const MAT_TEMPLATE *mt, INT sub, MATDATA_DESC **submd);

INT CreateFormatCmd                                             (INT argc, char **argv);
INT RemoveFormatWithSubs                                (const char *name);

INT CreateVecDescCmd                            (MULTIGRID *theMG, INT argc, char **argv);
INT CreateMatDescCmd                                (MULTIGRID *theMG, INT argc, char **argv);

INT InitFormats                                                 (void);

#endif
