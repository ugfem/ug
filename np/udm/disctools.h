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

#define MAX_NODAL_VECTORS 20
#define MAX_NODAL_VALUES  60

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

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
INT GetElementsideIndices       (ELEMENT *theElement, INT side,
                                 const VECDATA_DESC *theTVD, INT *index);
INT GetElementVPtrs             (ELEMENT *theElement, const VECDATA_DESC *theTVD,
                                 DOUBLE **vptr);
INT GetElementVPtrsVecskip      (ELEMENT *theElement, const VECDATA_DESC *theTVD,
                                 DOUBLE **vptr, INT *vecskip);
INT GetElementNewVPtrs          (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                 DOUBLE **vptr, INT *new);
INT GetElementMPtrs                             (ELEMENT *theElement, const MATDATA_DESC *theTMD, DOUBLE **mptr);
INT GetElementVMPtrs            (ELEMENT *theElement,
                                 const VECDATA_DESC *theTVD, const MATDATA_DESC *theTMD,
                                 DOUBLE **vptr, DOUBLE **mptr);
INT GetElementVVMPtrs           (ELEMENT *theElement, const VECDATA_DESC *theTVD1,
                                 const VECDATA_DESC *theTVD2, const MATDATA_DESC *theTMD,
                                 DOUBLE **vptr1, DOUBLE **vptr2, DOUBLE **mptr,
                                 INT *vecskip);
INT ClearVecskipFlags           (GRID *theGrid);
INT GetElementDirichletFlags    (ELEMENT *theElement, const VECDATA_DESC *theTVD,
                                 INT *vecskip);
INT SetElementDirichletFlags    (ELEMENT *theElement, const VECDATA_DESC *theTVD,
                                 INT *vecskip);
INT AssembleDirichletBoundary   (GRID *theGrid, const MATDATA_DESC *Mat,
                                 const VECDATA_DESC *Sol, const VECDATA_DESC *Rhs);
INT AssembleTotalDirichletBoundary (GRID *theGrid, const MATDATA_DESC *Mat,
                                    const VECDATA_DESC *Sol, const VECDATA_DESC *Rhs);

INT PrintVector (GRID *g, VECDATA_DESC *X, INT vclass, INT vnclass);
INT PrintMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass);
INT PrintIMatrix (GRID *g, VECDATA_DESC *V, INT vclass, INT vnclass);


#endif
