// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  udm.h															*/
/*																			*/
/* Purpose:   user data manager (header file)								*/
/*																			*/
/* Author:	  Peter Bastian													*/
/*																			*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/*																			*/
/* History:   02.12.96 begin, ug version 3.4								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UDM__
#define __UDM__

#include "compiler.h"
#include "gm.h"

/****************************************************************************/
/*																			*/
/* macros concerned with data descriptors                                                       */
/*																			*/
/****************************************************************************/

#define NVECTYPES                               MAXVECTORS
#define NMATTYPES                               (MAXVECTORS*MAXVECTORS)
#define MAX_SINGLE_VEC_COMP              9      /* max nb of comp in one TYPE                   */
#define MAX_SINGLE_MAT_COMP             81      /* max nb of comp in one TYPE		    */
#define MAX_VEC_COMP                     9      /* max nb of comp in one VECDATA_DESC	*/
#define MAX_MAT_COMP                    81      /* max nb of comp in one VECDATA_DESC	*/

#define NVECOFFSETS                             (NVECTYPES+1)
/* for offset component in VECDATA_DESC	*/
#define NMATOFFSETS                             (NMATTYPES+1)

#define DEFAULT_NAMES "uvwzpqrst"   /* of size MAX_VEC_COMP                 */

/* VECDATA_DESC */
#define VD_ISDEF_IN_TYPE(vd,tp)             (VD_NCMPS_IN_TYPE(vd,tp)>0)
#define VD_NCMPPTR(vd)                              ((vd)->NCmpInType)
#define VD_NCMPS_IN_TYPE(vd,tp)             (VD_NCMPPTR(vd)[tp])
#define VD_CMP_OF_TYPE(vd,tp,i)             ((vd)->CmpsInType[tp][i])
#define VD_CMPPTR_OF_TYPE(vd,tp)            ((vd)->CmpsInType[tp])

#define VD_IS_SCALAR(vd)                    ((vd)->IsScalar)
#define VD_SCALCMP(vd)                                          ((vd)->ScalComp)
#define VD_SCALTYPEMASK(vd)                                     ((vd)->ScalTypeMask)
#define VD_OFFSETPTR(vd)                    ((vd)->offset)
#define VD_OFFSET(vd,tp)                    (VD_OFFSETPTR(vd)[tp])
#define VD_NCOMP(vd)                        (VD_OFFSETPTR(vd)[NVECTYPES])

/* MATDATA_DESC */
#define MTP(rt,ct)                          ((rt)*NVECTYPES+(ct))
#define MCMP(row,col,ncol)                  ((row)*(ncol)+col)
#define MTYPE_RT(mtp)                                           ((mtp)/NVECTYPES)
#define MTYPE_CT(mtp)                                           ((mtp)%NVECTYPES)
#define MCMP_I(mc,ncol)                                         ((mc)/(ncol))
#define MCMP_J(mc,ncol)                                         ((mc)%(ncol))

#define MD_ISDEF_IN_MTYPE(md,mtp)           (MD_ROWS_IN_MTYPE(md,mtp)>0)
#define MD_ISDEF_IN_RT_CT(md,rt,ct)         MD_ISDEF_IN_MTYPE(md,MTP(rt,ct))
#define MD_ROWPTR(md)                               ((md)->RowsInType)
#define MD_ROWS_IN_MTYPE(md,mtp)            (MD_ROWPTR(md)[mtp])
#define MD_ROWS_IN_RT_CT(md,rt,ct)          MD_ROWS_IN_MTYPE(md,MTP(rt,ct))
#define MD_COLPTR(md)                               ((md)->ColsInType)
#define MD_COLS_IN_MTYPE(md,mtp)            (MD_COLPTR(md)[mtp])
#define MD_COLS_IN_RT_CT(md,rt,ct)          MD_COLS_IN_MTYPE(md,MTP(rt,ct))
#define MD_NCMPS_IN_MTYPE(md,mtp)           MD_ROWS_IN_MTYPE(md,mtp)*MD_COLS_IN_MTYPE(md,mtp)
#define MD_NCMPS_IN_RT_CT(md,rt,ct)         MD_NCMPS_IN_MTYPE(md,MTP(rt,ct))
#define MD_MCMPPTR_OF_MTYPE(md,mtp)         ((md)->CmpsInType[mtp])
#define MD_MCMP_OF_MTYPE(md,mtp,i)          ((md)->CmpsInType[mtp][i])
#define MD_MCMPPTR_OF_RT_CT(md,rt,ct)       ((md)->CmpsInType[MTP(rt,ct)])
#define MD_MCMP_OF_RT_CT(md,rt,ct,i)        MD_MCMP_OF_MTYPE(md,MTP(rt,ct),i)
#define MD_IJ_CMP_OF_MTYPE(md,mtp,i,j)      MD_MCMP_OF_MTYPE(md,mtp,MCMP(i,j,MD_COLS_IN_MTYPE(md,mtp)))
#define MD_IJ_CMP_OF_RT_CT(md,rt,ct,i,j)    MD_IJ_CMP_OF_MTYPE(md,MTP(rt,ct),i,j)

#define MD_IS_SCALAR(md)                    ((md)->IsScalar)
#define MD_SCALCMP(md)                                          ((md)->ScalComp)
#define MD_SCAL_RTYPEMASK(md)                           ((md)->ScalRowTypeMask)
#define MD_SCAL_CTYPEMASK(md)                           ((md)->ScalColTypeMask)
#define MD_OFFSETPTR(md)                        ((md)->offset)
#define MD_MTYPE_OFFSET(md,mtp)             (MD_OFFSETPTR(md)[mtp])
#define MD_RT_CT_OFFSET(md,mtp)             (MD_MTYPE_OFFSET(md,MTP(rt,ct)))
#define MD_SUCC_COMP(md)                    ((md)->SuccComp)

/* VEC_SCALAR */
#define VS_CMP_AT_OFFSET(vs,os,i)                       (vs[(os)+(i)])
#define VS_CMP_OF_TYPE(vs,vd,tp,i)                      VS_CMP_AT_OFFSET(vs,VD_OFFSET(vd,tp),i)

#define VM_COMP_NAMEPTR(p)                 ((p)->compNames)
#define VM_COMP_NAME(p,i)                  (VM_COMP_NAMEPTR(p)[i])
#define VM_COMPPTR(p)                      ((p)->Components)
#define VM_LOCKED(p)                       ((p)->locked)

/****************************************************************************/
/*																			*/
/* data structures                                                                                                                      */
/*																			*/
/****************************************************************************/

typedef struct {

  /* fields for environment list variable */
  ENVVAR v;

  INT locked;                          /* locked for dynamic allocation         */
  char compNames[MAX_VEC_COMP];    /* names for symbol components           */
  SHORT NCmpInType[NVECTYPES];     /* number of components of a vector      */
                                   /* per type                              */
  SHORT *CmpsInType[NVECTYPES];    /* pointer to SHORT vector containing    */
  /*    the components                     */
  /* redundant (but frequently used) information                          */
  SHORT IsScalar;                  /* TRUE if desc is scalar:               */
                                   /*  same settings in all types           */
  SHORT ScalComp;                  /* location of scalar component          */
  INT ScalTypeMask;                /* mask for used vectypes                */
  SHORT offset[NVECOFFSETS];       /* offsets for VEC_SCALARs               */
  SHORT Components[1];                 /* memory for component mapping	        */

} VECDATA_DESC;

typedef struct {

  ENVVAR v;

  INT locked;                          /* locked for dynamic allocation         */
  char compNames[2*MAX_MAT_COMP];   /* names for symbol components          */
  SHORT RowsInType[NMATTYPES];          /* number of rows of a matrix per type  */
  SHORT ColsInType[NMATTYPES];          /* number of columns of a matrix        */
                                        /* per type                             */
  SHORT *CmpsInType[NMATTYPES];         /* pointer to SHORT vector containing   */
                                        /* the components                       */
  /* redundant (but frequently used) information                          */
  SHORT IsScalar;                       /* TRUE if desc is scalar:              */
  /* same settings in all types           */
  SHORT SuccComp;                   /* successive components                */
  SHORT ScalComp;                       /* location of scalar component         */
  INT ScalRowTypeMask;                  /* mask for used vectypes in rows       */
  INT ScalColTypeMask;                  /* mask for used vectypes in cols       */
  SHORT offset[NMATOFFSETS];            /* offsets for what ever you need it    */
  SHORT Components[1];                  /* memory for component mapping	        */

} MATDATA_DESC;


typedef DOUBLE VEC_SCALAR[MAX_VEC_COMP];

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

VECDATA_DESC *GetFirstVector (MULTIGRID *theMG);
VECDATA_DESC *GetNextVector (VECDATA_DESC *vd);
MATDATA_DESC *GetFirstMatrix (MULTIGRID *theMG);
MATDATA_DESC *GetNextMatrix (MATDATA_DESC *md);

VECDATA_DESC *CreateVecDesc (MULTIGRID *theMG, const char *name, const char *compNames,
                             const SHORT *NCmpInType);
MATDATA_DESC *CreateMatDesc (MULTIGRID *theMG, const char *name, const char *compNames,
                             const SHORT *RowsInType, const SHORT *ColsInType);
VECDATA_DESC *CreateSubVecDesc (MULTIGRID *theMG, const VECDATA_DESC *theVD, const char *name,
                                const SHORT *NCmpInType, const SHORT *Comps, const char *CompNames);
MATDATA_DESC *CreateSubMatDesc (MULTIGRID *theMG, const MATDATA_DESC *theMD,
                                const char *name, const SHORT *RowsInType,
                                const SHORT *ColsInType, const SHORT *Comps, const char *CompNames);

INT FillRedundantComponentsOfVD (VECDATA_DESC *vd);
INT FillRedundantComponentsOfMD (MATDATA_DESC *md);

INT DisplayVecDataDesc (const VECDATA_DESC *vd, char *buffer);
INT DisplayMatDataDesc (const MATDATA_DESC *md, char *buffer);

VECDATA_DESC *GetVecDataDescByName (const MULTIGRID *theMG, char *name);
MATDATA_DESC *GetMatDataDescByName (const MULTIGRID *theMG, char *name);

/* allocation of vector descriptors */
INT AllocVDfromNCmp (MULTIGRID *theMG, INT fl, INT tl,
                     const SHORT *NCmpInType, const char *compNames, VECDATA_DESC **new_desc);
INT AllocVDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                   const VECDATA_DESC *template_desc, VECDATA_DESC **new_desc);

/* allocation of matrix descriptors */
INT AllocMDFromMRowMCol (MULTIGRID *theMG, INT fl, INT tl,
                         const SHORT *RowsInType,const SHORT *ColsInType,const char *compNames,
                         MATDATA_DESC **new_desc);
INT AllocMDFromVRowVCol (MULTIGRID *theMG, INT fl, INT tl,
                         const SHORT *RowsInType, const SHORT *ColsInType, MATDATA_DESC **new_desc);
INT AllocMDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                   const VECDATA_DESC *x, const VECDATA_DESC *y, MATDATA_DESC **new_desc);
INT AllocMDFromMD (MULTIGRID *theMG, INT fl, INT tl,
                   const MATDATA_DESC *template_desc, MATDATA_DESC **new_desc);

/* freeing of vector and matrix descriptors */
INT FreeVD        (MULTIGRID *theMG, INT fl, INT tl, VECDATA_DESC *x);
INT FreeMD        (MULTIGRID *theMG, INT fl, INT tl, MATDATA_DESC *A);

INT ConstructVecOffsets (const SHORT *NCmpInType, SHORT *offset);
INT ConstructMatOffsets (const SHORT *RowsInType, const SHORT *ColsInType, SHORT *offset);


/* init user data manager */
INT InitUserDataManager (void);

#endif
