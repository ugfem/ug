// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      udm.h                                                         */
/*                                                                          */
/* Purpose:   user data manager (header file)                               */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*                                                                          */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/*                                                                          */
/* History:   02.12.96 begin, ug version 3.4                                */
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

#ifndef __UDM__
#define __UDM__

#include "compiler.h"
#include "gm.h"
#include "sm.h"

#include "namespace.h"

START_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* macros concerned with data descriptors                                   */
/*                                                                          */
/****************************************************************************/

#define NVECTYPES                       MAXVECTORS
#define NMATTYPES                       MAXCONNECTIONS
#define NMATTYPES_NORMAL        MAXMATRICES

#define MTP(rt,ct)          ((rt)*NVECTYPES+(ct))
#define DMTP(rt)            (NMATTYPES_NORMAL+rt)
#define MTYPE_RT(mtp)       (((mtp)<NMATTYPES_NORMAL) ? (mtp)/NVECTYPES : (mtp)%NVECTYPES)
#define MTYPE_CT(mtp)           ((mtp)%NVECTYPES)

/** \brief Max nb of vec comps in one TYPE      */
#define MAX_SINGLE_VEC_COMP             40
/** \brief Max nb of mat comps in one TYPE              */
#define MAX_SINGLE_MAT_COMP       1600
/** \brief Max nb of comps in one VECDATA_DESC  */
#define MAX_VEC_COMP                40
/** \brief Max nb of comps in one MATDATA_DESC  */
#define MAX_MAT_COMP              7000
/** \brief Max #(comp) in one MATDATA_DESC       */
#define MAX_MAT_COMP_TOTAL        7000

#define NVECOFFSETS                             (NVECTYPES+1)

/** \brief For offset component in VECDATA_DESC */
#define NMATOFFSETS                             (NMATTYPES+1)

#define DEFAULT_NAMES "uvwzpabcdefghijklmnoPQRSTUVWXYZ123456789" /* of size MAX_VEC_COMP */

/** \brief No identification of components */
#define NO_IDENT                        -1

/** \brief Full template rather than sub */
#define FULL_TPLT                       -1

#define GENERATED_NAMES_SEPERATOR               "_"

/** \brief Defines for getting object type specific information from XXXDATA_DESCs      */
enum VECTOR_DATA_IN_OBJ
{
  STRICT,
  NON_STRICT
};

/** \brief Modifier flags for DisplayXXXDataDesc */
enum DISP_DATA_DESC_MODIF
{
  DEFAULT                 = (1<<0),
  ALLOC_STAT              = (1<<1),
  SCAL_PROP               = (1<<2)
};

/** @name Defines for VECDATA_DESC */
/*@{*/
#define VD_MG(vd)                                                       ((vd)->mg)
#define VD_ISDEF_IN_TYPE(vd,tp)             (VD_NCMPS_IN_TYPE(vd,tp)>0)
#define VD_NCMPPTR(vd)                              ((vd)->NCmpInType)
#define VD_NCMPS_IN_TYPE(vd,tp)             (VD_NCMPPTR(vd)[tp])
#define VD_CMP_OF_TYPE(vd,tp,i)             ((vd)->CmpsInType[tp][i])
#define VD_CMPPTR_OF_TYPE(vd,tp)            ((vd)->CmpsInType[tp])

#define VD_NID(vd)                                                      ((vd)->nId)
#define VD_IDENT_PTR(vd)                                        ((vd)->Ident)
#define VD_IDENT(vd,i)                                          ((vd)->Ident[i])

#define VD_DATA_TYPES(vd)                                       ((vd)->datatypes)
#define VD_OBJ_USED(vd)                                         ((vd)->objused)

#define VD_IS_SCALAR(vd)                    ((vd)->IsScalar)
#define VD_SCALCMP(vd)                                          ((vd)->ScalComp)
#define VD_SCALTYPEMASK(vd)                                     ((vd)->ScalTypeMask)
#define VD_OFFSETPTR(vd)                    ((vd)->offset)
#define VD_OFFSET(vd,tp)                    (VD_OFFSETPTR(vd)[tp])
#define VD_NCOMP(vd)                        (VD_OFFSETPTR(vd)[NVECTYPES])
#define VD_MINTYPE(vd)                      ((vd)->mintype)
#define VD_MAXTYPE(vd)                      ((vd)->maxtype)
#define VD_SUCC_COMP(vd)                    ((vd)->SuccComp)
/*@}*/


/** @name Macros for MATDATA_DESC */
/*@{*/
#define MCMP(row,col,ncol)                  ((row)*(ncol)+col)
#define MCMP_I(mc,ncol)                                         ((mc)/(ncol))
#define MCMP_J(mc,ncol)                                         ((mc)%(ncol))

#define MD_MG(md)                                                       ((md)->mg)
#define MD_ISDEF_IN_MTYPE(md,mtp)           (MD_ROWS_IN_MTYPE(md,mtp)>0)
#define MD_ISDEF_IN_RT_CT(md,rt,ct)         MD_ISDEF_IN_MTYPE(md,MTP(rt,ct))
#define MD_ROWPTR(md)                               ((md)->RowsInType)
#define MD_ROWS_IN_MTYPE(md,mtp)            (MD_ROWPTR(md)[mtp])
#define MD_ROWS_IN_RT_CT(md,rt,ct)          MD_ROWS_IN_MTYPE(md,MTP(rt,ct))
#define MD_COLPTR(md)                               ((md)->ColsInType)
#define MD_COLS_IN_MTYPE(md,mtp)            (MD_COLPTR(md)[mtp])
#define MD_COLS_IN_RT_CT(md,rt,ct)          MD_COLS_IN_MTYPE(md,MTP(rt,ct))
#define MD_NCMPS_IN_MTYPE(md,mtp)           (MD_ROWS_IN_MTYPE(md,mtp)*MD_COLS_IN_MTYPE(md,mtp))
#define MD_NCMPS_IN_RT_CT(md,rt,ct)         MD_NCMPS_IN_MTYPE(md,MTP(rt,ct))
#define MD_MCMPPTR(md)                      ((md)->CmpsInType)
#define MD_MCMPPTR_OF_MTYPE(md,mtp)         ((md)->CmpsInType[mtp])
#define MD_MCMP_OF_MTYPE(md,mtp,i)          ((md)->CmpsInType[mtp][i])
#define MD_MCMPPTR_OF_RT_CT(md,rt,ct)       ((md)->CmpsInType[MTP(rt,ct)])
#define MD_MCMP_OF_RT_CT(md,rt,ct,i)        MD_MCMP_OF_MTYPE(md,MTP(rt,ct),i)
#define MD_IJ_CMP_OF_MTYPE(md,mtp,i,j)      MD_MCMP_OF_MTYPE(md,mtp,MCMP(i,j,MD_COLS_IN_MTYPE(md,mtp)))
#define MD_IJ_CMP_OF_RT_CT(md,rt,ct,i,j)    MD_IJ_CMP_OF_MTYPE(md,MTP(rt,ct),i,j)

#define MD_ROW_DATA_TYPES(md)                           ((md)->rowdatatypes)
#define MD_COL_DATA_TYPES(md)                           ((md)->coldatatypes)
#define MD_ROW_OBJ_USED(md)                                     ((md)->rowobjused)
#define MD_COL_OBJ_USED(md)                                     ((md)->colobjused)
#define MD_SM(md,i)                                 ((md)->sm[i])
#define MD_SMP(md)                                  ((md)->sm)

#define MD_IS_SPARSE(md)                    ((md)->IsSparse)
#define MD_IS_SCALAR(md)                    ((md)->IsScalar)
#define MD_SCALCMP(md)                                          ((md)->ScalComp)
#define MD_SCAL_RTYPEMASK(md)                           ((md)->ScalRowTypeMask)
#define MD_SCAL_CTYPEMASK(md)                           ((md)->ScalColTypeMask)
#define MD_OFFSETPTR(md)                        ((md)->offset)
#define MD_MTYPE_OFFSET(md,mtp)             (MD_OFFSETPTR(md)[mtp])
#define MD_RT_CT_OFFSET(md,mtp)             (MD_MTYPE_OFFSET(md,MTP(rt,ct)))
#define MD_SUCC_COMP(md)                    ((md)->SuccComp)
/*@}*/

/** @name Macros for VEC_SCALAR */
/*@{*/
#define VS_CMP_AT_OFFSET(vs,os,i)                       (vs[(os)+(i)])
#define VS_CMP_OF_TYPE(vs,vd,tp,i)                      VS_CMP_AT_OFFSET(vs,VD_OFFSET(vd,tp),i)

#define VM_COMP_NAMEPTR(p)                 ((p)->compNames)
#define VM_COMP_NAME(p,i)                  (VM_COMP_NAMEPTR(p)[i])
#define VM_COMPPTR(p)                      ((p)->Components)
/*@}*/

/* we remove this for security reasons: please use function calls
   #define VM_LOCKED(p)                       ((p)->locked)
 */
#define EVM_LOCKED(p)                      ((p)->locked)

/** @name Swapping part interface data */
/*@{*/
#define SPID_NVD_MAX            4
#define SPID_NMD_MAX            2
#define SPID_FORTH                      69
#define SPID_BACK                       96

#define SPID_NVD(p)                     ((p)->nvd)
#define SPID_VD(p,i)            ((p)->vd[i])
#define SPID_VDI(p,i)           ((p)->vdi[i])
#define SPID_NMD(p)                     ((p)->nmd)
#define SPID_MD(p,i)            ((p)->md[i])
#define SPID_MDI(p,i)           ((p)->mdi[i])
/*@}*/

/** @name Flags for consistent status and collect status */
/*@{*/
#define READ_VEC_CONS_FLAG(p,vt,i)      READ_FLAG((p)->data_status.VecConsistentStatus[vt][(i)/32],1<<((i)%32))
#define SET_VEC_CONS__FLAG(p,vt,i)      SET_FLAG((p)->data_status.VecConsistentStatus[vt][(i)/32],1<<((i)%32))
#define CLEAR_VEC_CONS__FLAG(p,vt,i)    CLEAR_FLAG((p)->data_status.VecConsistentStatus[vt][(i)/32],1<<((i)%32))

#define READ_VEC_COLLECT_FLAG(p,vt,i)   READ_FLAG((p)->data_status.VecCollectStatus[vt][(i)/32],1<<((i)%32))
#define SET_VEC_COLLECT__FLAG(p,vt,i)   SET_FLAG((p)->data_status.VecCollectStatus[vt][(i)/32],1<<((i)%32))
#define CLEAR_VEC_COLLECT__FLAG(p,vt,i) CLEAR_FLAG((p)->data_status.VecCollectStatus[vt][(i)/32],1<<((i)%32))
/*@}*/

/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

/** \brief A set of degrees of freedom associated to a geometrical object
 */
typedef struct {

  /** \brief Fields for environment list variable */
  ENVVAR v;

  /** \brief Locked for dynamic allocation */
  SHORT locked;

  /** \brief Associated multigrid */
  MULTIGRID *mg;

  /** \brief Names for symbol components */
  char compNames[MAX_VEC_COMP];

  /** \brief Number of components of a vector per type */
  SHORT NCmpInType[NVECTYPES];

  /** \brief Pointer to SHORT vector containing the components */
  SHORT *CmpsInType[NVECTYPES];

  /* redundant (but frequently used) information */
  /** \brief TRUE if desc is scalar: same settings in all types */
  SHORT IsScalar;

  /** \brief Successive components */
  SHORT SuccComp;

  /** \brief Location of scalar component */
  SHORT ScalComp;

  /** \brief Mask for used vectypes */
  SHORT ScalTypeMask;

  /** \brief Offsets for VEC_SCALARs */
  SHORT offset[NVECOFFSETS];

  /** \brief Compact form of vtypes (bitwise) */
  SHORT datatypes;

  /** \brief Compact form of otypes (bitwise) */
  SHORT objused;

  /** \brief Minimal used type */
  SHORT mintype;

  /** \brief Maximal used type */
  SHORT maxtype;

  /** \brief Number of comps after ident */
  SHORT nId;

  /** \brief Identification table */
  SHORT *Ident;

  /** \brief Memory for component mapping */
  SHORT Components[1];

} VECDATA_DESC;

typedef struct {

  /** \brief Inheritance from environment variable class */
  ENVVAR v;

  /** \brief Locked for dynamic allocation        */
  SHORT locked;

  /** \brief Associated multigrid */
  MULTIGRID *mg;

  /** \brief Names for symbol components          */
  char compNames[2*MAX_MAT_COMP];

  /** \brief Number of rows of a matrix per type  */
  SHORT RowsInType[NMATTYPES];

  /** \brief Number of columns of a matrix per type */
  SHORT ColsInType[NMATTYPES];

  /** \brief Pointer to SHORT vector containing the components */
  SHORT *CmpsInType[NMATTYPES];

  /** \brief Pointers to sm form, if not full */
  SPARSE_MATRIX *sm[NMATTYPES];

  /* redundant (but frequently used) information                          */
  /** \brief TRUE if sparse form should be used   */
  SHORT IsSparse;

  /** \brief TRUE if desc is scalar: same settings in all types             */
  SHORT IsScalar;

  /** \brief Successive components                */
  SHORT SuccComp;

  /** \brief Location of scalar component         */
  SHORT ScalComp;

  /** \brief Mask for used vectypes in rows       */
  SHORT ScalRowTypeMask;

  /** \brief Mask for used vectypes in cols       */
  SHORT ScalColTypeMask;

  /** \brief Offsets for what ever you need it    */
  SHORT offset[NMATOFFSETS];

  /** \brief Compact form of row vtypes (bitwise)     */
  SHORT rowdatatypes;

  /** \brief Compact form of col vtypes (bitwise)     */
  SHORT coldatatypes;

  /** \brief Compact form of row otypes (bitwise) */
  SHORT rowobjused;

  /** \brief Compact form of col otypes (bitwise) */
  SHORT colobjused;

  /** \brief Memory for component mapping         */
  SHORT Components[1];

} MATDATA_DESC;


typedef DOUBLE VEC_SCALAR[MAX_VEC_COMP];

typedef struct {

  /** \brief Number of vector data descriptors               */
  INT nvd;

  /** \brief Vector data descriptors                                 */
  VECDATA_DESC *vd[SPID_NVD_MAX];

  /** \brief Vector data interface descriptors               */
  VECDATA_DESC *vdi[SPID_NVD_MAX];

  /** \brief Number of matrix data descriptors               */
  INT nmd;

  /** \brief Matrix data descriptors                                 */
  MATDATA_DESC *md[SPID_NMD_MAX];

  /** \brief Matrix data interface descriptors               */
  MATDATA_DESC *mdi[SPID_NMD_MAX];

} SPID_DESC;


#define EXTENSION_MAX                   10

#define EVDD_E(v,l,i)                           ((v)->e[l][i])
#define EVDD_E_PTR(v,l)                         (&EVDD_E(v,l,0))
#define EVD_NCOMP(v)                (VD_NCOMP((v)->vd)+(v)->n)

typedef struct {

  /** \brief Fields for environment list variable */
  ENVVAR v;

  /** \brief Locked for dynamic allocation        */
  SHORT locked;

  /** \brief Size of extension                    */
  INT n;

  /** \brief Vector descriptor                    */
  VECDATA_DESC *vd;

  /** \brief Extension                         */
  DOUBLE e[MAXLEVEL][EXTENSION_MAX];

} EVECDATA_DESC;

#define EMDD_EE(m,l,i)                          ((m)->ee[l][i])
#define EMDD_EE_PTR(m,l)                        (&EMDD_EE(m,l,0))

typedef struct {

  /** \brief Fields for environment list variable */
  ENVVAR v;

  /** \brief Locked for dynamic allocation        */
  SHORT locked;

  /** \brief Size of extension                    */
  INT n;

  /** \brief 11-block                     */
  MATDATA_DESC *mm;

  /** \brief 12-block                     */
  VECDATA_DESC *me[EXTENSION_MAX];

  /** \brief 21-block                     */
  VECDATA_DESC *em[EXTENSION_MAX];

  /** \brief 22-block           */
  DOUBLE ee[MAXLEVEL][EXTENSION_MAX*EXTENSION_MAX];

} EMATDATA_DESC;

typedef DOUBLE EVEC_SCALAR[MAX_VEC_COMP+EXTENSION_MAX];



/****************************************************************************/
/*                                                                          */
/* definition of exported data structures                                   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported functions                                         */
/*                                                                          */
/****************************************************************************/

VECDATA_DESC *GetFirstVector (MULTIGRID *theMG);
VECDATA_DESC *GetNextVector  (VECDATA_DESC *vd);
MATDATA_DESC *GetFirstMatrix (MULTIGRID *theMG);
MATDATA_DESC *GetNextMatrix  (MATDATA_DESC *md);

VECDATA_DESC *CreateVecDesc (MULTIGRID *theMG, const char *name, const char *compNames,
                             const SHORT *NCmpInType, SHORT nId, SHORT *Ident);
MATDATA_DESC *CreateMatDesc (MULTIGRID *theMG, const char *name, const char *compNames,
                             const SHORT *RowsInType, const SHORT *ColsInType,
                             SHORT **CmpsInType);
VECDATA_DESC *CreateSubVecDesc (MULTIGRID *theMG, const char *name,
                                const SHORT *NCmpInType, const SHORT *Comps, const char *CompNames);
MATDATA_DESC *CreateSubMatDesc (MULTIGRID *theMG, const char *name, const char *CompNames,
                                const SHORT *RowsInType, const SHORT *ColsInType, SHORT **CmpsInType);
VECDATA_DESC *CombineVecDesc (MULTIGRID *theMG, const char *name, const VECDATA_DESC **theVDs,
                              const INT nrOfVDs);

INT VDequal (const VECDATA_DESC *vd0, const VECDATA_DESC *vd1);

INT FillRedundantComponentsOfVD (VECDATA_DESC *vd);
INT FillRedundantComponentsOfMD (MATDATA_DESC *md);

INT DisplayVecDataDesc (const VECDATA_DESC *vd, INT modifiers, char *buffer);
INT DisplayMatDataDesc (const MATDATA_DESC *md, char *buffer);

VECDATA_DESC *GetVecDataDescByName (const MULTIGRID *theMG, char *name);
MATDATA_DESC *GetMatDataDescByName (const MULTIGRID *theMG, char *name);

INT CompMatDesc (const MATDATA_DESC *md, const SHORT *RowsInType,
                 const SHORT *ColsInType, SHORT *const*CmpsInType);

/** @name Allocation of vector descriptors */
/*@{*/
INT AllocVDfromNCmp (MULTIGRID *theMG, INT fl, INT tl,
                     const SHORT *NCmpInType, const char *compNames, VECDATA_DESC **new_desc);
INT AllocVDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                   const VECDATA_DESC *template_desc, VECDATA_DESC **new_desc);
/*@}*/

/** @name Allocation of extended vector descriptors */
/*@{*/
INT AllocEVDForVD (MULTIGRID *theMG, const VECDATA_DESC *vd, INT n, EVECDATA_DESC **new_desc);
INT AllocEVDFromEVD (MULTIGRID *theMG, INT fl, INT tl, const EVECDATA_DESC *vd, EVECDATA_DESC **new_desc);
/*@}*/

/** @name Allocation of matrix descriptors */
/*@{*/
INT AllocMDFromMRowMCol (MULTIGRID *theMG, INT fl, INT tl,
                         const SHORT *RowsInType,const SHORT *ColsInType,const char *compNames,
                         MATDATA_DESC **new_desc);
INT AllocMDFromVRowVCol (MULTIGRID *theMG, INT fl, INT tl,
                         const SHORT *RowsInType, const SHORT *ColsInType, MATDATA_DESC **new_desc);
INT AllocMDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                   const VECDATA_DESC *x, const VECDATA_DESC *y, MATDATA_DESC **new_desc);
INT AllocMDFromMD (MULTIGRID *theMG, INT fl, INT tl,
                   MATDATA_DESC *template_desc, MATDATA_DESC **new_desc);
/*@}*/

/** @name Allocation of extended matrix descriptors */
/*@{*/
INT AllocEMDFromEVD (MULTIGRID *theMG, INT fl, INT tl, const EVECDATA_DESC *x, const EVECDATA_DESC *y, EMATDATA_DESC **new_desc);
INT AllocEMDForMD (MULTIGRID *theMG, const MATDATA_DESC *md, INT n, EMATDATA_DESC **new_desc);
/*@}*/

/** @name Locking of vector and matrix descriptors */
/*@{*/
INT LockVD (MULTIGRID *theMG, VECDATA_DESC *vd);
INT LockMD (MATDATA_DESC *md);
INT UnlockMD (MATDATA_DESC *md);
/*@}*/

/** @name ??? */
/*@{*/
INT TransmitLockStatusVD (const VECDATA_DESC *vd, VECDATA_DESC *svd);
INT TransmitLockStatusMD (const MATDATA_DESC *md, MATDATA_DESC *smd);
/*@}*/

/** @name Freeing of vector and matrix descriptors */
/*@{*/
INT FreeVD        (MULTIGRID *theMG, INT fl, INT tl, VECDATA_DESC *x);
INT FreeMD        (MULTIGRID *theMG, INT fl, INT tl, MATDATA_DESC *A);
/*@}*/

/** @name Freeing of extended vector and matrix descriptors */
/*@{*/
INT FreeEVD       (MULTIGRID *theMG, INT fl, INT tl, EVECDATA_DESC *vd);
INT FreeEMD       (MULTIGRID *theMG, INT fl, INT tl, EMATDATA_DESC *md);
/*@}*/

/** \brief Interpolate allocation on new level */
INT InterpolateVDAllocation (MULTIGRID *theMG, VECDATA_DESC *vd);

/** @name Disposing of vector and matrix descriptors */
/*@{*/
INT DisposeVD     (VECDATA_DESC *vd);
INT DisposeMD     (MATDATA_DESC *md);
/*@}*/

/** @name Constructing part interface descriptors */
/*@{*/
INT VDinterfaceDesc                             (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi);
INT VDinterfaceCoDesc                   (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi);
INT VDCoDesc                                    (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi);
INT MDinterfaceDesc                             (const MATDATA_DESC *md, const MATDATA_DESC *mds, MATDATA_DESC **mdi);
INT MDinterfaceCoCoupleDesc             (const MATDATA_DESC *md, const MATDATA_DESC *mds, MATDATA_DESC **mdi);
/*@}*/

INT ConstructVecOffsets         (const SHORT *NCmpInType, SHORT *offset);
INT ConstructMatOffsets         (const SHORT *RowsInType, const SHORT *ColsInType, SHORT *offset);
INT ConstructMatOffsetsAlt      (const SHORT *CmpsInType, SHORT *offset);

/* swapping data on part interfaces */
/*@{*/
INT SwapPartInterfaceData       (INT fl, INT tl, SPID_DESC *spid, INT direction);
INT SwapPartSkipflags           (INT fl, INT tl, const VECDATA_DESC *vdg, const VECDATA_DESC *vdi, INT direction);
/*@}*/

/****************************************************************************/
/*      getting object type specific information from XXXDATA_DESCs         */
/****************************************************************************/

/** @name vtypes and object types */
/*@{*/
INT             GetUniqueOTypeOfVType                                           (const FORMAT *fmt, INT vtype);
INT             GetUniquePartOfVType                                            (const MULTIGRID *mg, INT vtype);
INT             IsVDdefinedInAllObjects                                         (const MULTIGRID *mg, const VECDATA_DESC *vd, INT obj_flags);
INT             FillCompsForOType                                                       (const FORMAT *fmt, INT otype, INT n, SHORT cmps[]);
/*@}*/

/** @name VECDATA_DESCs and object type */
/*@{*/
INT             VD_ncmps_in_otype_mod                                           (const VECDATA_DESC *vd, INT otype, INT mode);
INT             VD_cmp_of_otype_mod                                                     (const VECDATA_DESC *vd, INT otype, INT i, INT mode);
SHORT  *VD_ncmp_cmpptr_of_otype_mod                             (const VECDATA_DESC *vd, INT otype, INT *ncomp, INT mode);
INT             VDusesVOTypeOnly                                                        (const VECDATA_DESC *vd, INT votype);

#define VD_ncmps_in_otype(vd,ot)                                        VD_ncmps_in_otype_mod(vd,ot,STRICT)
#define VD_cmp_of_otype(vd,ot,i)                                        VD_cmp_of_otype_mod(vd,ot,i,STRICT)
#define VD_cmpptr_of_otype(vd,ot)                                       VD_ncmp_cmpptr_of_otype_mod(vd,ot,NULL,STRICT)
#define VD_cmpptr_of_otype_mod(vd,ot,mo)                        VD_ncmp_cmpptr_of_otype_mod(vd,ot,NULL,mo)
#define VD_ncmp_cmpptr_of_otype(vd,ot,nc)                       VD_ncmp_cmpptr_of_otype_mod(vd,ot,nc,STRICT)
/*@}*/

/** @name MATDATA_DESCs and object type */
/*@{*/
INT MD_rows_in_ro_co_mod                   (const MATDATA_DESC *md, INT rowobj, INT colobj, INT mode);
INT MD_cols_in_ro_co_mod                   (const MATDATA_DESC *md, INT rowobj, INT colobj, INT mode);
INT MD_rows_cols_in_ro_co_mod             (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nr, INT *nc, INT mode);
INT MD_mcmp_of_ro_co_mod                   (const MATDATA_DESC *md, INT rowobj, INT colobj, INT i, INT mode);
SHORT *MD_nr_nc_mcmpptr_of_ro_co_mod      (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nrow, INT *ncol, INT mode);
INT    MDusesVOTypeOnly                     (const MATDATA_DESC *md, INT votype);

#define MD_rows_in_ro_co(md,ro,co)                                      MD_rows_in_ro_co_mod(md,ro,co,STRICT)
#define MD_cols_in_ro_co(md,ro,co)                                      MD_cols_in_ro_co_mod(md,ro,co,STRICT)
#define MD_rows_cols_in_ro_co(md,ro,co,nr,nc)           MD_rows_cols_in_ro_co_mod(md,ro,co,nr,nc,STRICT)
#define MD_mcmp_of_ro_co(md,ro,co,i)                            MD_mcmp_of_ro_co_mod(md,ro,co,i,STRICT)
#define MD_mcmpptr_of_ro_co(md,ro,co)                           MD_nr_nc_mcmpptr_of_ro_co_mod(md,ro,co,NULL,NULL,STRICT)
#define MD_mcmpptr_of_ro_co_mod(md,ro,co,mo)            MD_nr_nc_mcmpptr_of_ro_co_mod(md,ro,co,NULL,NULL,mo)
#define MD_nr_nc_mcmpptr_of_ro_co(md,ro,co,nr,nc)       MD_nr_nc_mcmpptr_of_ro_co_mod(md,ro,co,nr,nc,STRICT)
/*@}*/

/** \brief Init user data manager */
INT InitUserDataManager (void);

END_NAMESPACE

#endif
