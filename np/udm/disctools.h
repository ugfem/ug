// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      disctools.h                                                   */
/*                                                                          */
/* Purpose:   tools for assembling (header file)		                        */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   Nov 27 95                                                                                 */
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

#ifndef __DISCTOOLS__
#define __DISCTOOLS__

#include "compiler.h"
#include "gm.h"
#include "np.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#ifdef __TWODIM__
#define MAX_NODAL_VECTORS       8
#define MAX_NODAL_VALUES        20
#define MAX_BND_VECTORS         4
#endif

#ifdef __THREEDIM__
#define MAX_NODAL_VECTORS       20
#define MAX_NODAL_VALUES        68
#define MAX_BND_VECTORS         8
#endif

#define MAXVD                           10
#define MAXMD                           5

#define MVMD_NVD(p)                     ((p)->nvd)
#define MVMD_NMD(p)                     ((p)->nmd)
#define MVMD_VD(p,j)            ((p)->vdlist[j])
#define MVMD_MD(p,j)            ((p)->mdlist[j])

#define MVMD_TYPES(p)           ((p)->types)
#define MVMD_TYPE(p,i)          (MVMD_TYPES(p)[i])
#define MVMD_DATATYPES(p)       ((p)->datatypes)
#define MVMD_OBJTYPES(p)        ((p)->objtypes)
#define MVMD_VDSUBSEQ(p,i)      ((p)->vdsubseq[i])
#define MVMD_MDSUBSEQ(p,i)      ((p)->mdsubseq[i])
#define MVMD_M_OF_1_ONLY(p)     ((p)->MatOfFirstVecOnly)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef struct {

  /* to be filled before call of PrepareElementMultipleVMPtrs */
  INT nvd;
  const VECDATA_DESC *vdlist[MAXVD];
  INT nmd;
  const MATDATA_DESC *mdlist[MAXMD];

  /* filled by PrepareElementMultipleVMPtrs */
  SHORT types[NVECTYPES];
  INT datatypes;
  INT objtypes;
  INT vdsubseq[MAXVD];
  INT mdsubseq[MAXMD];
  INT MatOfFirstVecOnly;                /* default FALSE, may be changed after call	*/
  /* of PrepareElementMultipleVMPtrs			*/

} MVM_DESC;

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

CoeffProcPtr MG_GetCoeffFct     (MULTIGRID *theMG, INT n);
UserProcPtr MG_GetUserFct       (MULTIGRID *theMG, INT n);
INT GetElementVertices          (ELEMENT *theElement, DOUBLE **x);
INT GetAllVectorsOfElementOfType(ELEMENT *theElement, VECTOR **vec,
                                 const VECDATA_DESC *theVD);
INT GetAllVectorsOfElementsideOfType (ELEMENT *theElement, INT side,
                                      VECTOR **vec,
                                      const VECDATA_DESC *theVD);
INT GetElementsideIndices       (ELEMENT *theElement, INT side,
                                 const VECDATA_DESC *theTVD, INT *index);

/* getting (pointers to) data corresponding to a (fixed, s.b.) set of XXXDATA_DESCs */
INT GetElementVPtrs             (ELEMENT *theElement, const VECDATA_DESC *theTVD,
                                 DOUBLE **vptr);
INT GetElementVValues           (ELEMENT *theElement,
                                 const VECDATA_DESC *theVD, DOUBLE *value);
INT AddElementVValues           (ELEMENT *theElement,
                                 const VECDATA_DESC *theVD, DOUBLE *value);
INT GetVlistVValues             (INT cnt, VECTOR **theVec,
                                 const VECDATA_DESC *theVD, DOUBLE *value);
INT AddVlistVValues             (INT cnt, VECTOR **theVec,
                                 const VECDATA_DESC *theVD, DOUBLE *value);
INT SetVlistVValues             (INT cnt, VECTOR **theVec,
                                 const VECDATA_DESC *theVD, DOUBLE *value);
INT GetElementVPtrsVecskip      (ELEMENT *theElement, const VECDATA_DESC *theTVD,
                                 DOUBLE **vptr, INT *vecskip);
INT GetElementNewVPtrs          (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                 DOUBLE **vptr, INT *newField);
INT GetElementMPtrs                             (ELEMENT *theElement,
                                                 const MATDATA_DESC *theTMD, DOUBLE **mptr);
INT GetVlistMValues             (INT cnt, VECTOR **theVec,
                                 const MATDATA_DESC *theMD, DOUBLE *value);
INT AddVlistMValues             (GRID *theGrid, INT cnt, VECTOR **theVec,
                                 const MATDATA_DESC *theMD, DOUBLE *value);
INT GetElementVMPtrs            (ELEMENT *theElement,
                                 const VECDATA_DESC *theTVD, const MATDATA_DESC *theTMD,
                                 DOUBLE **vptr, DOUBLE **mptr);
INT GetElementVVMPtrs           (ELEMENT *theElement, const VECDATA_DESC *theTVD1,
                                 const VECDATA_DESC *theTVD2, const MATDATA_DESC *theTMD,
                                 DOUBLE **vptr1, DOUBLE **vptr2, DOUBLE **mptr,
                                 INT *vecskip);

/* getting pointers to data corresponding to a variable set of XXXDATA_DESCs */
/* CAUTION: these fcts. follow other conventions than GetElementVVMPtrs... */
INT PrepareElementMultipleVMPtrs (MVM_DESC *mvmd);
INT GetElementMultipleVMPtrs    (ELEMENT *elem, const MVM_DESC *mvmd,
                                 DOUBLE **vptrlist[MAXVD],
                                 DOUBLE **mptrlist[MAXMD],
                                 INT *vecskip, INT *vtype, INT nvec[MAXVD]);
INT PrepareBndVecMultipleVMPtrs (GRID *theGrid, MVM_DESC *mvmd);
INT GetBndVecMultipleVMPtrs             (const MVM_DESC *mvmd,
                                         INT *cnt,
                                         VECTOR *VecList[],
                                         DOUBLE **vptrlist[MAXVD],
                                         DOUBLE **mptrlist[MAXMD],
                                         INT *vecskip, INT *vtype, INT nvec[MAXVD], INT *end);
INT ResetBndVecMultipleVMPtrs   (void);
INT FinishBndVecMultipleVMPtrs  (void);

/* skip flags */
INT ComputePartVecskip                  (const VECDATA_DESC *vd, const VECDATA_DESC *vds, INT typeskip[NVECTYPES], INT co_typeskip[NVECTYPES]);
INT ClearPartVecskipFlags               (GRID *theGrid, const INT typeskip[NVECTYPES]);
INT ClearVecskipFlags           (GRID *theGrid, const VECDATA_DESC *theVD);
INT GetElementDirichletFlags    (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                 INT *vecskip);
INT SetElementDirichletFlags    (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                 INT *vecskip);

/* modifications in Dirichlet dofs */
INT ModifyDirichletMatrix               (GRID *theGrid, const MATDATA_DESC *Mat);
INT ModifyDirichletDefect               (GRID *theGrid, const VECDATA_DESC *Def);
INT ClearDirichletValues                (GRID *theGrid, VECDATA_DESC *x);
INT AssembleDirichletBoundary   (GRID *theGrid, const MATDATA_DESC *Mat,
                                 const VECDATA_DESC *Sol, const VECDATA_DESC *Rhs);
INT AssembleTotalDirichletBoundary (GRID *theGrid, const MATDATA_DESC *Mat,
                                    const VECDATA_DESC *Sol, const VECDATA_DESC *Rhs);

/* display data */
INT PrintVector (GRID *g, VECDATA_DESC *X, INT vclass, INT vnclass);
INT PrintSVector (MULTIGRID *mg, VECDATA_DESC *X);
INT PrintMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass);
INT PrintTMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass);
INT PrintDiagMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass);
INT PrintIMatrix (GRID *g, VECDATA_DESC *V, INT vclass, INT vnclass);

#endif
