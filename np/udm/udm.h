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


/* RCS_ID
   $Header$
 */

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
#define VD_MG(vd)                                                       ((vd)->mg)
#define VD_ISDEF_IN_TYPE(vd,tp)             (VD_NCMPS_IN_TYPE(vd,tp)>0)
#define VD_NCMPPTR(vd)                              ((vd)->NCmpInType)
#define VD_NCMPS_IN_TYPE(vd,tp)             (VD_NCMPPTR(vd)[tp])
#define VD_CMP_OF_TYPE(vd,tp,i)             ((vd)->CmpsInType[tp][i])
#define VD_CMPPTR_OF_TYPE(vd,tp)            ((vd)->CmpsInType[tp])

#define VD_DATA_TYPES(vd)                                       ((vd)->datatypes)
#define VD_OBJ_USED(vd)                                         ((vd)->objused)

#define VD_IS_SCALAR(vd)                    ((vd)->IsScalar)
#define VD_SCALCMP(vd)                                          ((vd)->ScalComp)
#define VD_SCALTYPEMASK(vd)                                     ((vd)->ScalTypeMask)
#define VD_OFFSETPTR(vd)                    ((vd)->offset)
#define VD_OFFSET(vd,tp)                    (VD_OFFSETPTR(vd)[tp])
#define VD_NCOMP(vd)                        (VD_OFFSETPTR(vd)[NVECTYPES])
#define VD_MAXTYPE(vd)                      ((vd)->maxtype)

/* MATDATA_DESC */
#define MTP(rt,ct)                          ((rt)*NVECTYPES+(ct))
#define MCMP(row,col,ncol)                  ((row)*(ncol)+col)
#define MTYPE_RT(mtp)                                           ((mtp)/NVECTYPES)
#define MTYPE_CT(mtp)                                           ((mtp)%NVECTYPES)
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

/* swapping part interface data */
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

/* flags for constistent status and collect status */

#define READ_VEC_CONS_FLAG(p,vt,i)      READ_FLAG((p)->data_status.VecConsistentStatus[vt][(i)/32],1<<((i)%32))
#define SET_VEC_CONS__FLAG(p,vt,i)      SET_FLAG((p)->data_status.VecConsistentStatus[vt][(i)/32],1<<((i)%32))
#define CLEAR_VEC_CONS__FLAG(p,vt,i)    CLEAR_FLAG((p)->data_status.VecConsistentStatus[vt][(i)/32],1<<((i)%32))

#define READ_VEC_COLLECT_FLAG(p,vt,i)   READ_FLAG((p)->data_status.VecCollectStatus[vt][(i)/32],1<<((i)%32))
#define SET_VEC_COLLECT__FLAG(p,vt,i)   SET_FLAG((p)->data_status.VecCollectStatus[vt][(i)/32],1<<((i)%32))
#define CLEAR_VEC_COLLECT__FLAG(p,vt,i) CLEAR_FLAG((p)->data_status.VecCollectStatus[vt][(i)/32],1<<((i)%32))

/****************************************************************************/
/*																			*/
/* data structures                                                                                                                      */
/*																			*/
/****************************************************************************/

typedef struct {

  /* fields for environment list variable */
  ENVVAR v;

  INT locked;                          /* locked for dynamic allocation         */
  MULTIGRID *mg;                                   /* associated multigrid					*/
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

  INT datatypes;                                   /* compact form of vtypes (bitwise)		*/
  INT objused;                                     /* compact form of otypes (bitwise)		*/
  INT maxtype;                                     /* maximal used type                         */

  SHORT Components[1];                 /* memory for component mapping	        */

} VECDATA_DESC;

typedef struct {

  ENVVAR v;

  INT locked;                           /* locked for dynamic allocation        */
  MULTIGRID *mg;                                    /* associated multigrid					*/
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

  INT rowdatatypes;                                /* compact form of row vtypes (bitwise)	*/
  INT coldatatypes;                                /* compact form of col vtypes (bitwise)	*/
  INT rowobjused;                                  /* compact form of row otypes (bitwise)	*/
  INT colobjused;                                  /* compact form of col otypes (bitwise)	*/

  SHORT Components[1];                  /* memory for component mapping	        */

} MATDATA_DESC;


typedef DOUBLE VEC_SCALAR[MAX_VEC_COMP];

/* special const pointers */
typedef const VECDATA_DESC *CONST_VECDATA_DESC_PTR;
typedef const MATDATA_DESC *CONST_MATDATA_DESC_PTR;

typedef struct {

  INT nvd;                                                      /* number of vec data descriptors		*/
  VECDATA_DESC *vd[SPID_NVD_MAX];       /* vec data descriptors					*/
  VECDATA_DESC *vdi[SPID_NVD_MAX];      /* vec data interface descriptors		*/

  INT nmd;                                                      /* number of mat data descriptors		*/
  MATDATA_DESC *md[SPID_NMD_MAX];       /* mat data descriptors					*/
  MATDATA_DESC *mdi[SPID_NMD_MAX];      /* mat data interface descriptors		*/

} SPID_DESC;

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
VECDATA_DESC *CreateSubVecDesc (MULTIGRID *theMG, const char *name,
                                const SHORT *NCmpInType, const SHORT *Comps, const char *CompNames);
MATDATA_DESC *CreateSubMatDesc (MULTIGRID *theMG,
                                const char *name, const SHORT *RowsInType,
                                const SHORT *ColsInType, const SHORT *Comps, const char *CompNames);
VECDATA_DESC *CombineVecDesc (MULTIGRID *theMG, const char *name, const VECDATA_DESC **theVDs,
                              const INT nrOfVDs);

INT VDequal (const VECDATA_DESC *vd0, const VECDATA_DESC *vd1);

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

/* disposing of vector and matrix descriptors */
INT DisposeVD     (VECDATA_DESC *vd);
INT DisposeMD     (MATDATA_DESC *md);

/* constructing part interface descriptors */
INT VDinterfaceDesc                                             (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi);
INT VDinterfaceCoDesc                                   (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi);
INT MDinterfaceDesc                                             (const MATDATA_DESC *md, const MATDATA_DESC *mds, MATDATA_DESC **mdi);

INT ConstructVecOffsets         (const SHORT *NCmpInType, SHORT *offset);
INT ConstructMatOffsets         (const SHORT *RowsInType, const SHORT *ColsInType, SHORT *offset);
INT ConstructMatOffsetsAlt      (const SHORT *CmpsInType, SHORT *offset);

/* swapping data on part interfaces */
INT SwapPartInterfaceData       (INT fl, INT tl, SPID_DESC *spid, INT direction);
INT SwapPartSkipflags           (INT fl, INT tl, const VECDATA_DESC *vdg, const VECDATA_DESC *vdi, INT direction);

/****************************************************************************/
/*	getting object type specific information from XXXDATA_DESCs
 */

/* vtypes and object types */
INT             GetUniqueOTypeOfVType           (const FORMAT *fmt, INT vtype);
INT             GetUniquePartOfVType            (const MULTIGRID *mg, INT vtype);
INT             FillCompsForOType                       (const FORMAT *fmt, INT otype, INT n, SHORT cmps[]);

/* VECDATA_DESCs and object type */
INT             VD_ncmps_in_otype                       (const VECDATA_DESC *vd, INT otype);
INT             VD_cmp_of_otype                         (const VECDATA_DESC *vd, INT otype, INT i);
#define VD_cmpptr_of_otype(vd,ot)       VD_ncmp_cmpptr_of_otype(vd,ot,NULL)
SHORT   *VD_ncmp_cmpptr_of_otype        (const VECDATA_DESC *vd, INT otype, INT *ncmp);
INT             VDusesVOTypeOnly                        (const VECDATA_DESC *vd, INT votype);

/* MATDATA_DESCs and object type */
INT             MD_rows_in_ro_co                        (const MATDATA_DESC *md, INT rowobj, INT colobj);
INT             MD_cols_in_ro_co                        (const MATDATA_DESC *md, INT rowobj, INT colobj);
INT             MD_rows_cols_in_ro_co           (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nr, INT *nc);
INT             MD_mcmp_of_ro_co                        (const MATDATA_DESC *md, INT rowobj, INT colobj, INT i);
#define MD_mcmpptr_of_ro_co(md,ro,co)   MD_nr_nc_mcmpptr_of_ro_co(md,ro,co,NULL,NULL)
SHORT   *MD_nr_nc_mcmpptr_of_ro_co      (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nrow, INT *ncol);
INT             MDusesVOTypeOnly                        (const MATDATA_DESC *md, INT votype);

/* init user data manager */
INT InitUserDataManager (void);

#endif
