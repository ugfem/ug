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
#define MAX_SUB                         5

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
#define VF_COMPS(vt)            ((vt)->Comp)
#define VF_COMP(vt,tp)          ((vt)->Comp[tp])
#define VF_COMPNAMES(vt)        ((vt)->CompNames)
#define VF_COMPNAME(vt,i)       ((vt)->CompNames[i])
#define VF_SUB(vt,i)            ((vt)->SubVec[i])
#define VF_NSUB(vt)                     ((vt)->nsub)

/* macros for MAT_TEMPLATE */
#define MF_RCOMPS(mt)           ((mt)->RComp)
#define MF_RCOMP(mt,tp)         ((mt)->RComp[tp])
#define MF_CCOMPS(mt)           ((mt)->CComp)
#define MF_CCOMP(mt,tp)         ((mt)->CComp[tp])
#define MF_COMPNAMES(mt)        ((mt)->CompNames)
#define MF_COMPNAME(mt,i)       ((mt)->CompNames[i])
#define MF_SUB(mt,i)            ((mt)->SubMat[i])
#define MF_NSUB(mt)                     ((mt)->nsub)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef struct {

  char Name[NAMESIZE];
  SHORT Comp[NVECTYPES];
  SHORT Comps[NVECTYPES][MAX_VEC_COMP];

} SUBVEC;

typedef struct {

  char Name[NAMESIZE];
  SHORT RComp[NMATTYPES];
  SHORT CComp[NMATTYPES];
  SHORT Comps[NMATTYPES][MAX_MAT_COMP];

} SUBMAT;

typedef struct {

  ENVITEM v;

  SHORT Comp[NVECTYPES];
  char CompNames[MAX_VEC_COMP];

  SHORT nsub;
  SUBVEC  *SubVec[MAX_SUB];

} VEC_TEMPLATE;

typedef struct {

  ENVITEM v;

  SHORT RComp[NMATTYPES];
  SHORT CComp[NMATTYPES];
  char CompNames[2*MAX_MAT_COMP];

  SHORT nsub;
  SUBMAT  *SubMat[MAX_SUB];

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

VEC_TEMPLATE *GetVectorTemplate                 (const FORMAT *theFmt, const char *template);
MAT_TEMPLATE *GetMatrixTemplate                 (const FORMAT *theFmt, const char *template);

VECDATA_DESC *CreateVecDescOfTemplate   (MULTIGRID *theMG,
                                         const char *name, const char *template);
MATDATA_DESC *CreateMatDescOfTemplate   (MULTIGRID *theMG,
                                         const char *name, const char *template);

INT VDsubDescFromVT                                             (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt, INT sub, VECDATA_DESC **subvd);
INT MDsubDescFromVT                                             (const MATDATA_DESC *md, const VEC_TEMPLATE *vt, INT sub, MATDATA_DESC **submd);

INT CreateFormatCmd                                             (INT argc, char **argv);
INT CreateVecDescCmd                            (MULTIGRID *theMG, INT argc, char **argv);
INT CreateMatDescCmd                                (MULTIGRID *theMG, INT argc, char **argv);

INT InitFormats                                                 (void);

#endif
