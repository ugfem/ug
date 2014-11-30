// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      disctools.c                                                   */
/*                                                                          */
/* Purpose:   tools for assembling                                          */
/*                                                                          */
/* Author:    Christian Wieners                                             */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/* email:     ug@ica3.uni-stuttgart.de                                      */
/*                                                                          */
/* History:   Nov 27 95 begin                                               */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <config.h>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "gm.h"       /* for data structure               */
#include "evm.h"      /* for data structure               */
#include "general.h"
#include "ugdevices.h"

#include "ugblas.h"
#include "disctools.h"
#ifdef ModelP
#include "parallel.h"   /* for PRIO */
#endif

USING_UG_NAMESPACES

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

  REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/** \brief Get function pointer

   \param theMG - pointer to a multigrid
   \param n - number of coefficient function

   This function returns a pointer to the nth coefficient function
   of the multigrid.

   RETURN VALUE:
   .n    pointer to function
   .n    NULL if n is too large
 */
/****************************************************************************/

CoeffProcPtr NS_DIM_PREFIX MG_GetCoeffFct (const MULTIGRID *theMG, INT n)
{
  BVP *myBVP;
  const BVP_DESC *BVP_desc;
  CoeffProcPtr cpp;

  myBVP = MG_BVP(theMG);
  BVP_desc = MG_BVPD(theMG);
  cpp = NULL;
  if ((n >= 0) && (n < BVPD_NCOEFFF(BVP_desc)))
    BVP_SetCoeffFct       (myBVP,n,&cpp);
  return(cpp);
}

/****************************************************************************/
/** \brief Get function pointer

   \param theMG - pointer to a multigrid
   \param n - number of a user function

   This function returns a pointer to the nth user function
   of the multigrid.

   RETURN VALUE:
   .n    pointer to function
   .n    NULL if n is too large
 */
/****************************************************************************/

UserProcPtr NS_DIM_PREFIX MG_GetUserFct (MULTIGRID *theMG, INT n)
{
  BVP *myBVP;
  BVP_DESC *BVP_desc;
  UserProcPtr upp;

  myBVP = MG_BVP(theMG);
  BVP_desc = MG_BVPD(theMG);
  upp = NULL;
  if ((n >= 0) && (n < BVPD_NUSERF(BVP_desc)))
    BVP_SetUserFct (myBVP,n,&upp);
  return(upp);
}

/****************************************************************************/
/** \brief Get vector list for an element
 *
 * @param theElement pointer to an element
   \param vec - vector list
   \param theVD - vector descriptor

   This function gets a list of vectors corresponding to an element.
   It uses GetVectorsOfDataTypesInObjects (which should be preferred).

   RETURN VALUE:
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetAllVectorsOfElementOfType (ELEMENT *theElement, VECTOR **vec,
                                                const VECDATA_DESC *theVD)
{
  INT cnt;

  if (GetVectorsOfDataTypesInObjects(theElement,VD_DATA_TYPES(theVD),VD_OBJ_USED(theVD),&cnt,vec))
    return (-1);

  return (cnt);
}

/****************************************************************************/
/** \brief Get vector list for an element side
 *
 * @param theElement pointer to an element
   \param side the number of the element side
   \param vec - vector list
   \param theVD - vector descriptor

   This function gets a list of vectors corresponding to an element.
   It uses GetVectorsOfDataTypesInObjects (which should be preferred).

   RETURN VALUE:
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetAllVectorsOfElementsideOfType (ELEMENT *theElement, INT side,
                                                    VECTOR **vec,
                                                    const VECDATA_DESC *theVD)
{
  VECTOR *v,*vv[MAX_CORNERS_OF_ELEM];
  INT i,n;
  INT obj = VD_OBJ_USED(theVD);
  INT cnt = 0;

  if (obj & BITWISE_TYPE(NODEVEC)) {
    if (GetVectorsOfNodes(theElement,&n,vv) != GM_OK)
      return(GM_ERROR);
    for (i=0; i<CORNERS_OF_SIDE(theElement,side); i++) {
      v = vv[CORNER_OF_SIDE(theElement,side,i)];
      if (VD_NCMPS_IN_TYPE(theVD,VTYPE(v)))
        vec[cnt++] = v;
    }
  }
  if (obj & BITWISE_TYPE(EDGEVEC)) {
    if (GetVectorsOfEdges(theElement,&n,vv) != GM_OK)
      return(GM_ERROR);
    for (i=0; i<EDGES_OF_SIDE(theElement,side); i++) {
      v = vv[EDGE_OF_SIDE(theElement,side,i)];
      if (VD_NCMPS_IN_TYPE(theVD,VTYPE(v)))
        vec[cnt++] = v;
    }
  }
  if (obj & BITWISE_TYPE(ELEMVEC)) {
    if (GetVectorsOfElement(theElement,&n,vec+cnt) != GM_OK)
      return(GM_ERROR);
    if (VD_NCMPS_IN_TYPE(theVD,VTYPE(vec[cnt])))
      cnt++;
  }
    #ifdef __THREEDIM__
  if (obj & BITWISE_TYPE(SIDEVEC)) {
    if (GetVectorsOfSides(theElement,&n,vec+cnt) != GM_OK)
      return(GM_ERROR);
    if (VD_NCMPS_IN_TYPE(theVD,VTYPE(vec[cnt])))
      cnt++;
  }
    #endif

  return (cnt);
}



/****************************************************************************/
/** \brief Compute vector indices of an element side

   \param theElement - pointer to an element
   \param side - element side
   \param theVD - type vector descriptor
   \param index - index list

   This function gets the indices of the vector pointers of an element
   corresponding to a side.

   RETURN VALUE:
   .n    number of components of the element side
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementsideIndices (ELEMENT *theElement, INT side,
                                         const VECDATA_DESC *theVD, INT *index)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  INT vncomp;
  INT i,j,k,l,m,cnt,votype;
  INT itype[NVECTYPES];

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  for (i=0; i<NVECTYPES; i++)
    itype[i] = 0;

  m = 0;
  k = 0;
  for (i=0; i<cnt; i++)
  {
    votype = VOTYPE(theVec[i]);
    vncomp = VD_NCMPS_IN_TYPE(theVD,VTYPE(theVec[i]));
    if (votype == NODEVEC)
      if (itype[votype] == 0)
        for (l=0; l<CORNERS_OF_SIDE(theElement,side); l++)
          for (j=0; j<vncomp; j++)
            index[k++] = m + CORNER_OF_SIDE(theElement,side,l)*vncomp + j;
    if (votype == EDGEVEC)
              #ifdef __TWODIM__
      if (itype[votype] == side)
        for (j=0; j<vncomp; j++)
          index[k++] = m + j;
          #else
      if (itype[EDGEVEC] == 0)
        for (l=0; l<EDGES_OF_SIDE(theElement,side); l++)
          for (j=0; j<vncomp; j++)
            index[k++] = m + EDGE_OF_SIDE(theElement,side,l)*vncomp + j;
          #endif
    if (votype == SIDEVEC)
              #ifdef __TWODIM__
      if (itype[votype] == side)
        for (j=0; j<vncomp; j++)
          index[k++] = m + j;
          #else
      if (itype[votype] == side)
        for (j=0; j<vncomp; j++)
          index[k++] = m + j;
          #endif
    itype[votype] += 1;
    m += vncomp;
  }

  return (k);
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for vectors

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param vptr - pointer to double values

   This function gets all local vector pointers corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementVPtrs (ELEMENT *theElement, const VECDATA_DESC *theVD, DOUBLE **vptr)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  INT i,j,m,cnt,vtype;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype = VTYPE(theVec[i]);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++)
      vptr[m++] = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,j));
  }

  return (m);
}

/****************************************************************************/
/** \brief Get list of DOUBLE values for vectors

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param value - pointer to double values

   This function gets all local vector values corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementVValues (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                     DOUBLE *value)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  DOUBLE *vptr;
  INT i,j,m,cnt,vtype;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype = VTYPE(theVec[i]);
    vptr = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,0));
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++) {
      value[m++] = *vptr;
      vptr++;
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief Add list of DOUBLE values for vectors

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param value - pointer to double values

   This function adds all local vector values corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX AddElementVValues (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                     DOUBLE *value)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  DOUBLE *vptr;
  INT i,j,m,cnt,vtype;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype = VTYPE(theVec[i]);
    vptr = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,0));
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++) {
      *vptr += value[m++];
      vptr++;
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief Get list of vecskip flags

   \param cnt - number of vectors
   \param theVec - vector list
   \param theVD - type vector descriptor
   \param vecskip - set 1 for DIRICHLET boundary, 0 else

   This function gets vecskip flags for a vector set.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetVlistVecskip (INT cnt, VECTOR **theVec, const VECDATA_DESC *theVD,
                                   INT *vecskip)
{
  INT i,j,m,vtype;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype = VTYPE(theVec[i]);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++) {
      vecskip[m++] = ((VECSKIP(theVec[i]) & (1<<j))!=0);
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief Set list of vecskip flags

   \param cnt - number of vectors
   \param theVec - vector list
   \param theVD - type vector descriptor
   \param vecskip - set 1 for DIRICHLET boundary, 0 else

   This function sets vecskip flags for a vector set.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX SetVlistVecskip (INT cnt, VECTOR **theVec, const VECDATA_DESC *theVD,
                                   INT *vecskip)
{
  INT i,j,m,vtype;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype = VTYPE(theVec[i]);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++) {
      if (vecskip[m++] == 1)
        VECSKIP(theVec[i]) |= (1<<j);
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief Get list of DOUBLE values for vectors

   \param cnt - number of vectors
   \param theVec - vector list
   \param theVD - type vector descriptor
   \param value - pointer to double values

   This function gets all local vector values corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetVlistVValues (INT cnt, VECTOR **theVec,
                                   const VECDATA_DESC *theVD, DOUBLE *value)
{
  DOUBLE *vptr;
  INT i,j,m,vtype;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype = VTYPE(theVec[i]);
    vptr = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,0));
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++) {
      value[m++] = *vptr;
      vptr++;
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief Add list of DOUBLE values for vectors

   \param cnt - number of vectors
   \param theVec - vector list
   \param theVD - type vector descriptor
   \param value - pointer to double values

   This function adds all local vector values corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX AddVlistVValues (INT cnt, VECTOR **theVec,
                                   const VECDATA_DESC *theVD, DOUBLE *value)
{
  DOUBLE *vptr;
  INT i,j,m,vtype;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype = VTYPE(theVec[i]);
    vptr = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,0));
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++) {
      *vptr += value[m++];
      vptr++;
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief Set list of DOUBLE values for vectors

   \param cnt - number of vectors
   \param theVec - vector list
   \param theVD - type vector descriptor
   \param value - pointer to double values

   This function sets all local vector values corresponding to a vector list.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX SetVlistVValues (INT cnt, VECTOR **theVec,
                                   const VECDATA_DESC *theVD, DOUBLE *value)
{
  DOUBLE *vptr;
  INT i,j,m,vtype;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype = VTYPE(theVec[i]);
    vptr = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,0));
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++) {
      *vptr = value[m++];
      vptr++;
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for vectors

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param vptr - pointer to double values
   \param vecskip - set 1 for DIRICHLET boundary, 0 else

   This function gets all local vector pointers corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementVPtrsVecskip (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                          DOUBLE **vptr, INT *vecskip)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  INT i,j,m,cnt,vtype;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype = VTYPE(theVec[i]);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++)
    {
      vptr[m] = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,j));
      vecskip[m++] = ((VECSKIP(theVec[i]) & (1<<j))!=0);
    }
  }

  return (m);
}

/****************************************************************************/
/** \brief ???

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param vptr - pointer to double values
   \param newField - set 1 for new vectors, 0 else

   This function gets all local vector pointers corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    0 if there is no new vector
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementNewVPtrs (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                      DOUBLE **vptr, INT *newField)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  INT i,j,m,cnt,vtype,found;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  found = 0;
  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype = VTYPE(theVec[i]);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++)
    {
      vptr[m] = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype,j));
      if ((newField[m++] = VNEW(theVec[i])) > 0)
        found++;
    }
  }

  if (found)
    return (m);
  else
    return(0);
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for matrices

   \param theElement - pointer to an element
   \param md - matrix data descriptor
   \param mptr - pointer to double values corresponding to the local stiffness matrix

   This function gets all local matrix pointers corresponding to an element.

   RETURN VALUE:
   .n    number of dofs
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementMPtrs (ELEMENT *theElement, const MATDATA_DESC *md,
                                   DOUBLE **mptr)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  MATRIX *theMatrix;
  INT vncomp[MAX_NODAL_VECTORS];
  INT vtype[MAX_NODAL_VECTORS];
  SHORT types[NVECTYPES];
  INT i,j,k,l,m,m1,m2,cnt;

  for (i=0; i<NVECTYPES; i++)
    types[i] = MD_ISDEF_IN_RT_CT(md,i,i);

  if (GetVectorsOfDataTypesInObjects(theElement,MD_ROW_DATA_TYPES(md),MD_ROW_OBJ_USED(md),&cnt,theVec)!=GM_OK)
    return (-1);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = MD_ROWS_IN_RT_CT(md,vtype[i],vtype[i]);
    m += vncomp[i];
  }

  m1 = 0;
  for (i=0; i<cnt; i++)
  {
    theMatrix = VSTART(theVec[i]);
    for (k=0; k<vncomp[i]; k++)
      for (l=0; l<vncomp[i]; l++)
        mptr[(m1+k)*m+m1+l] =
          MVALUEPTR(theMatrix,
                    MD_MCMP_OF_RT_CT(md,vtype[i],vtype[i],k*vncomp[i]+l));
    m2 = 0;
    for (j=0; j<i; j++)
    {
      if ((theMatrix = GetMatrix(theVec[i],theVec[j]))==NULL)
        return (-1);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m1+k)*m+m2+l] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(md,vtype[i],vtype[j],k*vncomp[j]+l));
      theMatrix = MADJ(theMatrix);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m2+l)*m+m1+k] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(md,vtype[i],vtype[j],l*vncomp[i]+k));
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

/****************************************************************************/
/** \brief Get list of DOUBLE values for matrices

   \param cnt - number of vectors
   \param theVec - vector list
   \param theMD - type matrix descriptor
   \param value - pointer to double values

   This function get all local matrix values corresponding to a vector list.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetVlistMValues (INT cnt, VECTOR **theVec,
                                   const MATDATA_DESC *theMD, DOUBLE *value)
{
  MATRIX *theMatrix;
  INT vncomp[MAX_NODAL_VECTORS];
  INT vtype[MAX_NODAL_VECTORS];
  const SHORT *Comp[MAX_NODAL_VECTORS][MAX_NODAL_VECTORS];
  INT i,j,k,l,m,m1,m2;
  DOUBLE *mptr;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = MD_ROWS_IN_RT_CT(theMD,vtype[i],vtype[i]);
    m += vncomp[i];
  }
  for (i=0; i<cnt; i++)
    for (j=0; j<cnt; j++)
      Comp[i][j] = MD_MCMPPTR_OF_MTYPE(theMD,MTP(vtype[i],vtype[j]));

  m1 = 0;
  for (i=0; i<cnt; i++) {
    theMatrix = VSTART(theVec[i]);
    mptr = MVALUEPTR(theMatrix,0);
    for (k=0; k<vncomp[i]; k++)
      for (l=0; l<vncomp[i]; l++)
        value[(m1+k)*m+m1+l] =
          mptr[Comp[i][i][k*vncomp[i]+l]];
    m2 = 0;
    for (j=0; j<i; j++) {
      theMatrix = GetMatrix(theVec[i],theVec[j]);
      if (theMatrix == NULL) {
        for (k=0; k<vncomp[i]; k++)
          for (l=0; l<vncomp[j]; l++)
            value[(m1+k)*m+m2+l] = value[(m2+l)*m+m1+k] = 0.0;
      }
      else {
        mptr = MVALUEPTR(theMatrix,0);
        for (k=0; k<vncomp[i]; k++)
          for (l=0; l<vncomp[j]; l++)
            value[(m1+k)*m+m2+l] =
              mptr[Comp[i][j][k*vncomp[j]+l]];
        mptr = MVALUEPTR(MADJ(theMatrix),0);
        for (k=0; k<vncomp[i]; k++)
          for (l=0; l<vncomp[j]; l++)
            value[(m2+l)*m+m1+k] =
              mptr[Comp[i][j][l*vncomp[i]+k]];
      }
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

/****************************************************************************/
/** \brief Add list of DOUBLE values for matrices

   \param theGrid - pointer to grid
   \param cnt - number of vectors
   \param theVec - vector list
   \param theMD - type matrix descriptor
   \param value - pointer to double values

   This function adds all local matrix values corresponding to a vector list.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX AddVlistMValues (GRID *theGrid, INT cnt, VECTOR **theVec,
                                   const MATDATA_DESC *theMD, DOUBLE *value)
{
  MATRIX *theMatrix;
  INT vncomp[MAX_NODAL_VECTORS];
  INT vtype[MAX_NODAL_VECTORS];
  const SHORT *Comp[MAX_NODAL_VECTORS][MAX_NODAL_VECTORS],*comp;
  INT i,j,k,l,m,m1,m2;
  DOUBLE *mptr;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = MD_ROWS_IN_RT_CT(theMD,vtype[i],vtype[i]);
    m += vncomp[i];
  }
  for (i=0; i<cnt; i++)
    for (j=0; j<cnt; j++)
      Comp[i][j] = MD_MCMPPTR_OF_MTYPE(theMD,MTP(vtype[i],vtype[j]));

  if (MD_SUCC_COMP(theMD)) {
    m1 = 0;
    for (i=0; i<cnt; i++) {
      theMatrix = VSTART(theVec[i]);
      mptr = MVALUEPTR(theMatrix,Comp[i][i][0]);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[i]; l++)
          mptr[k*vncomp[i]+l]
            += value[(m1+k)*m+m1+l];
      m2 = 0;
      for (j=0; j<i; j++) {
        GET_MATRIX(theVec[i],theVec[j],theMatrix);
        if (theMatrix == NULL) {
          theMatrix =
            CreateExtraConnection(theGrid,theVec[i],theVec[j]);
          if (theMatrix == NULL)
            return (-1);
        }
        mptr = MVALUEPTR(theMatrix,Comp[i][j][0]);
        for (k=0; k<vncomp[i]; k++)
          for (l=0; l<vncomp[j]; l++)
            mptr[k*vncomp[j]+l] += value[(m1+k)*m+m2+l];
        mptr = MVALUEPTR(MADJ(theMatrix),Comp[j][i][0]);
        for (k=0; k<vncomp[i]; k++)
          for (l=0; l<vncomp[j]; l++)
            mptr[l*vncomp[i]+k] += value[(m2+l)*m+m1+k];
        m2 += vncomp[j];
      }
      m1 += vncomp[i];
    }
    return (m);
  }
  m1 = 0;
  for (i=0; i<cnt; i++) {
    theMatrix = VSTART(theVec[i]);
    mptr = MVALUEPTR(theMatrix,0);
    comp = Comp[i][i];
    for (k=0; k<vncomp[i]; k++)
      for (l=0; l<vncomp[i]; l++)
        mptr[comp[k*vncomp[i]+l]]
          += value[(m1+k)*m+m1+l];
    m2 = 0;
    for (j=0; j<i; j++) {
      GET_MATRIX(theVec[i],theVec[j],theMatrix);
      if (theMatrix == NULL) {
        theMatrix = CreateExtraConnection(theGrid,theVec[i],theVec[j]);
        if (theMatrix == NULL)
          return (-1);
      }
      mptr = MVALUEPTR(theMatrix,0);
      comp = Comp[i][j];
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[comp[k*vncomp[j]+l]] += value[(m1+k)*m+m2+l];
      mptr = MVALUEPTR(MADJ(theMatrix),0);
      comp = Comp[j][i];
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[comp[l*vncomp[i]+k]] += value[(m2+l)*m+m1+k];
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

INT NS_DIM_PREFIX SetVlistMValues (GRID *theGrid, INT cnt, VECTOR **theVec,
                                   const MATDATA_DESC *theMD, DOUBLE *value)
{
  MATRIX *theMatrix;
  INT vncomp[MAX_NODAL_VECTORS];
  INT vtype[MAX_NODAL_VECTORS];
  const SHORT *Comp[MAX_NODAL_VECTORS][MAX_NODAL_VECTORS],*comp;
  INT i,j,k,l,m,m1,m2;
  DOUBLE *mptr;

  m = 0;
  for (i=0; i<cnt; i++) {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = MD_ROWS_IN_RT_CT(theMD,vtype[i],vtype[i]);
    m += vncomp[i];
  }
  for (i=0; i<cnt; i++)
    for (j=0; j<cnt; j++)
      Comp[i][j] = MD_MCMPPTR_OF_MTYPE(theMD,MTP(vtype[i],vtype[j]));

  if (MD_SUCC_COMP(theMD)) {
    m1 = 0;
    for (i=0; i<cnt; i++) {
      theMatrix = VSTART(theVec[i]);
      mptr = MVALUEPTR(theMatrix,Comp[i][i][0]);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[i]; l++)
          mptr[k*vncomp[i]+l]
            = value[(m1+k)*m+m1+l];
      m2 = 0;
      for (j=0; j<i; j++) {
        GET_MATRIX(theVec[i],theVec[j],theMatrix);
        if (theMatrix == NULL) {
          theMatrix =
            CreateExtraConnection(theGrid,theVec[i],theVec[j]);
          if (theMatrix == NULL)
            return (-1);
        }
        mptr = MVALUEPTR(theMatrix,Comp[i][j][0]);
        for (k=0; k<vncomp[i]; k++)
          for (l=0; l<vncomp[j]; l++)
            mptr[k*vncomp[j]+l] = value[(m1+k)*m+m2+l];
        mptr = MVALUEPTR(MADJ(theMatrix),Comp[j][i][0]);
        for (k=0; k<vncomp[i]; k++)
          for (l=0; l<vncomp[j]; l++)
            mptr[l*vncomp[i]+k] = value[(m2+l)*m+m1+k];
        m2 += vncomp[j];
      }
      m1 += vncomp[i];
    }
    return (m);
  }
  m1 = 0;
  for (i=0; i<cnt; i++) {
    theMatrix = VSTART(theVec[i]);
    mptr = MVALUEPTR(theMatrix,0);
    comp = Comp[i][i];
    for (k=0; k<vncomp[i]; k++)
      for (l=0; l<vncomp[i]; l++)
        mptr[comp[k*vncomp[i]+l]]
          = value[(m1+k)*m+m1+l];
    m2 = 0;
    for (j=0; j<i; j++) {
      GET_MATRIX(theVec[i],theVec[j],theMatrix);
      if (theMatrix == NULL) {
        theMatrix = CreateExtraConnection(theGrid,theVec[i],theVec[j]);
        if (theMatrix == NULL)
          return (-1);
      }
      mptr = MVALUEPTR(theMatrix,0);
      comp = Comp[i][j];
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[comp[k*vncomp[j]+l]] = value[(m1+k)*m+m2+l];
      mptr = MVALUEPTR(MADJ(theMatrix),0);
      comp = Comp[j][i];
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[comp[l*vncomp[i]+k]] = value[(m2+l)*m+m1+k];
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for vectors and matrices

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param md - type matrix descriptor
   \param vptr - pointer to double values corresponding to the local right hand side
   \param mptr - pointer to double values corresponding to the local stiffness matrix

   This function gets all local vector pointers corresponding to an element
   and the connecting matrix pointers.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementVMPtrs (ELEMENT *theElement,
                                    const VECDATA_DESC *theVD, const MATDATA_DESC *md,
                                    DOUBLE **vptr, DOUBLE **mptr)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  MATRIX *theMatrix;
  INT vncomp[MAX_NODAL_VECTORS];
  INT vtype[MAX_NODAL_VECTORS];
  INT i,j,k,l,m,m1,m2,cnt;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = VD_NCMPS_IN_TYPE (theVD,vtype[i]);
    for (j=0; j<vncomp[i]; j++)
      vptr[m++] = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD,vtype[i],j));
  }

  m1 = 0;
  for (i=0; i<cnt; i++)
  {
    theMatrix = VSTART(theVec[i]);
    for (k=0; k<vncomp[i]; k++)
      for (l=0; l<vncomp[i]; l++)
        mptr[(m1+k)*m+m1+l] =
          MVALUEPTR(theMatrix,
                    MD_MCMP_OF_RT_CT(md,vtype[i],vtype[i],k*vncomp[i]+l));
    m2 = 0;
    for (j=0; j<i; j++)
    {
      if ((theMatrix = GetMatrix(theVec[i],theVec[j]))==NULL)
        return (-1);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m1+k)*m+m2+l] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(md,vtype[i],vtype[j],k*vncomp[j]+l));
      theMatrix = MADJ(theMatrix);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m2+l)*m+m1+k] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(md,vtype[i],vtype[j],l*vncomp[i]+k));
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for vectors and matrices

   \param theElement - pointer to an element
   \param theVD1 - type vector descriptor
   \param theVD2 - type vector descriptor
   \param md - type matrix descriptor
   \param vptr1 - pointer to double values corresponding to the local right hand side
   \param vptr2 - pointer to double values corresponding to the local right hand side
   \param mptr - pointer to double values corresponding to the local stiffness matrix
   \param vecskip - set 1 for DIRICHLET boundary, 0 else

   This function gets all local vector pointers corresponding to an element
   and the connecting matrix pointers.

   RETURN VALUE:
   .n    total number of components if ok
   .n    -1: error in GetAllVectorsOfElementOfType
   .n    -2: vecdata descriptors of different size
   .n    -3: could not get matrix
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementVVMPtrs (ELEMENT *theElement, const VECDATA_DESC *theVD1,
                                     const VECDATA_DESC *theVD2, const MATDATA_DESC *md,
                                     DOUBLE **vptr1, DOUBLE **vptr2, DOUBLE **mptr,
                                     INT *vecskip)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  MATRIX *theMatrix;
  INT vncomp[MAX_NODAL_VECTORS];
  INT vtype[MAX_NODAL_VECTORS];
  INT i,j,k,l,m,m1,m2,cnt;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD1);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    REP_ERR_RETURN(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = VD_NCMPS_IN_TYPE (theVD1,vtype[i]);
    if (vncomp[i] != VD_NCMPS_IN_TYPE (theVD2,vtype[i]))
      REP_ERR_RETURN (-2);
    for (j=0; j<vncomp[i]; j++)
    {
      vptr1[m] = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD1,vtype[i],j));
      vptr2[m] = VVALUEPTR(theVec[i],VD_CMP_OF_TYPE(theVD2,vtype[i],j));
      vecskip[m++] = ((VECSKIP(theVec[i]) & (1<<j))!=0);
    }
  }

  m1 = 0;
  for (i=0; i<cnt; i++)
  {
    theMatrix = VSTART(theVec[i]);
    for (k=0; k<vncomp[i]; k++)
      for (l=0; l<vncomp[i]; l++)
        mptr[(m1+k)*m+m1+l] =
          MVALUEPTR(theMatrix,
                    MD_MCMP_OF_RT_CT(md,vtype[i],vtype[i],k*vncomp[i]+l));
    m2 = 0;
    for (j=0; j<i; j++)
    {
      if ((theMatrix = GetMatrix(theVec[i],theVec[j]))==NULL)
        REP_ERR_RETURN (-3);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m1+k)*m+m2+l] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(md,vtype[i],vtype[j],k*vncomp[j]+l));
      theMatrix = MADJ(theMatrix);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m2+l)*m+m1+k] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(md,vtype[i],vtype[j],l*vncomp[i]+k));
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

/****************************************************************************/
/** \brief Prepare execution of GetMultipleVMPtrs

   \param mvmd - multiple VM data structure, partially to be filled before call

   This function prepares the execution of GetMultipleVMPtrs in a
   loop. It has to be called once before running the loop. The types needed are
   evaluated and it is checked, whether the components described by the
   XXXDATA_DESCs passed are in subsequent order. In this case ONLY a pointer to
   the first component per type is returned. The pointer can be incremented by
   the user.

   RETURN VALUE:
   .n    0

   SEE ALSO:
   GetMultipleVMPtrs
 */
/****************************************************************************/

static INT PrepareMultipleVMPtrs (MVM_DESC *mvmd)
{
  FORMAT *fmt;
  INT tp,dt,ot,j;

  /* get format */
  if (MVMD_NVD(mvmd)>0)
    fmt = MGFORMAT(VD_MG(MVMD_VD(mvmd,0)));
  else if (MVMD_NMD(mvmd)>0)
    fmt = MGFORMAT(MD_MG(MVMD_MD(mvmd,0)));
  else
    /* no XXXDATA_DESCs defined at all */
    REP_ERR_RETURN (1);

  dt = ot = 0;
  for (j=0; j<MVMD_NVD(mvmd); j++)
  {
    MVMD_VDSUBSEQ(mvmd,j) = VD_SUCC_COMP(MVMD_VD(mvmd,j));
    dt |= VD_DATA_TYPES(MVMD_VD(mvmd,j));
    ot |= VD_OBJ_USED(MVMD_VD(mvmd,j));
  }
  for (j=0; j<MVMD_NMD(mvmd); j++)
  {
    MVMD_MDSUBSEQ(mvmd,j) = MD_SUCC_COMP(MVMD_MD(mvmd,j));
    dt |= MD_ROW_DATA_TYPES(MVMD_MD(mvmd,j));
    dt |= MD_COL_DATA_TYPES(MVMD_MD(mvmd,j));
    ot |= MD_ROW_OBJ_USED(MVMD_MD(mvmd,j));
    ot |= MD_COL_OBJ_USED(MVMD_MD(mvmd,j));
  }
  MVMD_DATATYPES(mvmd) = dt;
  MVMD_OBJTYPES(mvmd)  = ot;

  /* set vtypes used */
  for (tp=0; tp<NVECTYPES; tp++)
    if (READ_FLAG(dt,BITWISE_TYPE(tp)))
      MVMD_TYPE(mvmd,tp) = TRUE;
    else
      MVMD_TYPE(mvmd,tp) = FALSE;

  MVMD_M_OF_1_ONLY(mvmd) = FALSE;

  return (0);
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for vectors and matrices

   \param[in] mvmd - data filled by PrepareElementMultipleVMPtrs
   \param[in] cnt - number of vectors to extract pointers from
   \param[in] VecList - list of vectors to extract pointers from

   \param[out] vptrlist - pointer to lists of double values corresponding the  VECDATA_DESC-list
   \param[out] mptrlist - pointer to lists of double values corresponding the  MATDATA_DESC-list
   \param[out] vecskip - set 1 for DIRICHLET boundary, 0 else (ordering corresponds to the first
                                VECDATA_DESC)
   \param[out] vtype - types of vectors collected for first vd
   \param[out] nvec - number of vectors involved from per vd of mvmd

   This functions returns pointers to the data fields described in a VECDATA_DESC-list
   and a MATDATA_DESC-list passed in the mvmd argument. Before call of this function
   in a loop PrepareMultipleVMPtrs has to be called with mvmd.

   In the case that components are ordered subsequently only a pointer to the first
   DOUBLE value is returned. The pointers can be incremented by the user.

   RETURN VALUE:
   .n    0 if ok
   .n    >0: error

   SEE ALSO:
   PrepareMultipleVMPtrs
 */
/****************************************************************************/

static INT GetMultipleVMPtrs (const MVM_DESC *mvmd, INT cnt, VECTOR *VecList[],
                              DOUBLE **vptrlist[MAXVD],
                              DOUBLE **mptrlist[MAXMD],
                              INT *vecskip, INT *vtype, INT nvec[MAXVD])
{
  VECTOR *rv,*cv;
  MATRIX *mat;
  INT i,j,k,l,nskip,rt,ct;
  INT vc[MAXVD],mc[MAXMD];

  for (i=0; i<MAXVD; i++)
    nvec[i] = 0;

  nskip = 0;
  for (l=0; l<MVMD_NVD(mvmd); l++) vc[l] = 0;
  for (l=0; l<MVMD_NMD(mvmd); l++) mc[l] = 0;
  for (i=0; i<cnt; i++)
  {
    rv = VecList[i];
    rt = VTYPE(rv);

    /* read skip flags */
    if (MVMD_NVD(mvmd)>1)
      for (k=0; k<VD_NCMPS_IN_TYPE(MVMD_VD(mvmd,0),rt); k++)
        vecskip[nskip++] = ((VECSKIP(rv) & (1<<k))!=0);

    /* get VECDATA_DESC ptrs */
    for (l=0; l<MVMD_NVD(mvmd); l++)
      if (VD_ISDEF_IN_TYPE(MVMD_VD(mvmd,l),rt))
      {
        if (l==0)
          vtype[nvec[l]] = rt;
        nvec[l]++;
        if (MVMD_VDSUBSEQ(mvmd,l))
        {
          /* by convention only return ptr to first component */
          /* TODO: possibly increase speed by introducing offset pointers
                           VVALUEPTR(rv,VD_CMP_OF_TYPE(MVMD_VD(mvmd,l),rt,0));
                           which are computed by PrepareElementMultipleVMPtrs for
                           all vector types */
          vptrlist[l][vc[l]++] = VVALUEPTR(rv,VD_CMP_OF_TYPE(MVMD_VD(mvmd,l),rt,0));
        }
        else
        {
          /* fill DOUBLE pointers, subsequently all needed in rv */
          for (k=0; k<VD_NCMPS_IN_TYPE(MVMD_VD(mvmd,l),rt); k++)
            vptrlist[l][vc[l]++] = VVALUEPTR(rv,VD_CMP_OF_TYPE(MVMD_VD(mvmd,l),rt,k));
        }
      }

    if (MVMD_NMD(mvmd)<=0)
      continue;

    if (MVMD_M_OF_1_ONLY(mvmd) && (i>0))
      continue;

    /* now connections */

    /* diag first */
    mat = VSTART(rv);
    for (l=0; l<MVMD_NMD(mvmd); l++)
      if (MD_ISDEF_IN_RT_CT(MVMD_MD(mvmd,l),rt,rt))
        if (MVMD_MDSUBSEQ(mvmd,l))
        {
          /* by convention only return ptr to first component */
          mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),rt,rt,0));
        }
        else
        {
          /* fill DOUBLE pointers, subsequently all needed in mat */
          for (k=0; k<MD_NCMPS_IN_RT_CT(MVMD_MD(mvmd,l),rt,rt); k++)
            mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),rt,rt,k));
        }

    /* off diag follows */
    if (MVMD_M_OF_1_ONLY(mvmd))
      for (j=1; j<cnt; j++)
      {
        cv = VecList[j];
        ct = VTYPE(cv);

        mat = GetMatrix(rv,cv);
        if (mat==NULL) REP_ERR_RETURN(1);

        for (l=0; l<MVMD_NMD(mvmd); l++)
          if (MD_ISDEF_IN_RT_CT(MVMD_MD(mvmd,l),rt,ct))
            if (MVMD_MDSUBSEQ(mvmd,l))
            {
              /* by convention only return ptr to first component */
              mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),rt,ct,0));
            }
            else
            {
              /* fill DOUBLE pointers, subsequently all needed in mat */
              for (k=0; k<MD_NCMPS_IN_RT_CT(MVMD_MD(mvmd,l),rt,ct); k++)
                mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),rt,ct,k));
            }

        /* no adjoint matrix */
      }
    else
      for (j=0; j<i; j++)
      {
        cv = VecList[j];
        ct = VTYPE(cv);

        mat = GetMatrix(rv,cv);
        if (mat==NULL) REP_ERR_RETURN(1);

        for (l=0; l<MVMD_NMD(mvmd); l++)
          if (MD_ISDEF_IN_RT_CT(MVMD_MD(mvmd,l),rt,ct))
            if (MVMD_MDSUBSEQ(mvmd,l))
            {
              /* by convention only return ptr to first component */
              mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),rt,ct,0));
            }
            else
            {
              /* fill DOUBLE pointers, subsequently all needed in mat */
              for (k=0; k<MD_NCMPS_IN_RT_CT(MVMD_MD(mvmd,l),rt,ct); k++)
                mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),rt,ct,k));
            }

        /* adjoint matrix */
        mat = MADJ(mat);
        for (l=0; l<MVMD_NMD(mvmd); l++)
          if (MD_ISDEF_IN_RT_CT(MVMD_MD(mvmd,l),ct,rt))
            if (MVMD_MDSUBSEQ(mvmd,l))
            {
              /* by convention only return ptr to first component */
              mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),ct,rt,0));
            }
            else
            {
              /* fill DOUBLE pointers, subsequently all needed in mat */
              for (k=0; k<MD_NCMPS_IN_RT_CT(MVMD_MD(mvmd,l),ct,rt); k++)
                mptrlist[l][mc[l]++] = MVALUEPTR(mat,MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,l),ct,rt,k));
            }
      }
  }

  return (0);
}

/****************************************************************************/
/** \brief Prepare execution of GetElementMultipleVMPtrs

   \param mvmd - multiple VM data structure, partially to be filled before call

   This function prepares the execution of GetElementMultipleVMPtrs in an element
   loop. It has to be called once before running the loop.
   For further description see 'PrepareMultipleVMPtrs'.

   RETURN VALUE:
   .n    0

   SEE ALSO:
   PrepareMultipleVMPtrs, GetElementMultipleVMPtrs
 */
/****************************************************************************/

INT NS_DIM_PREFIX PrepareElementMultipleVMPtrs (MVM_DESC *mvmd)
{
  return (PrepareMultipleVMPtrs(mvmd));
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for vectors and matrices

   \param elem - pointer to an element
   \param mvmd - data filled by PrepareElementMultipleVMPtrs
   \param vptrlist - pointer to lists of double values corresponding the  VECDATA_DESC-list
   \param mptrlist - pointer to lists of double values corresponding the  MATDATA_DESC-list
   \param vecskip - set 1 for DIRICHLET boundary, 0 else (ordering corresponds to the first
                                VECDATA_DESC)
   \param vtype - types of vectors collected for first vd
   \param nvec - number of vectors involved from this element

   This functions returns pointers to the data fields described in a VECDATA_DESC-list
   and a MATDATA_DESC-list passed in the mvmd argument. Before call of this function
   in an element loop PrepareElementMultipleVMPtrs has to be called with mvmd.

   In the case that components are ordered subsequently only a pointer to the first
   DOUBLE value is returned. The pointers can be incremented by the user.

   RETURN VALUE:
   .n    0: ok
   .n    >0: error

   SEE ALSO:
   GetMultipleVMPtrs,PrepareElementMultipleVMPtrs
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementMultipleVMPtrs (ELEMENT *elem, const MVM_DESC *mvmd,
                                            DOUBLE **vptrlist[MAXVD],
                                            DOUBLE **mptrlist[MAXMD],
                                            INT *vecskip, INT *vtype, INT nvec[MAXVD])
{
  VECTOR *VecList[MAX_NODAL_VECTORS];
  INT cnt;


  if (GetVectorsOfDataTypesInObjects(elem,MVMD_DATATYPES(mvmd),MVMD_OBJTYPES(mvmd),&cnt,VecList)!=GM_OK)
    REP_ERR_RETURN (1);

  return (GetMultipleVMPtrs(mvmd,cnt,VecList,vptrlist,mptrlist,vecskip,vtype,nvec));
}

/****************************************************************************/
/** \brief Prepare execution of GetBndVecMultipleVMPtrs

   \param theGrid - grid level
   \param mvmd - multiple VM data structure, partially to be filled before call

   This function prepares the execution of GetBndVecMultipleVMPtrs in an element
   loop. It has to be called once before running the loop.
   For further description see 'PrepareMultipleVMPtrs'.

   RETURN VALUE:
   INT
   .n    0

   SEE ALSO:
   PrepareMultipleVMPtrs, GetBndVecMultipleVMPtrs
 */
/****************************************************************************/

INT NS_DIM_PREFIX PrepareBndVecMultipleVMPtrs (GRID *theGrid, MVM_DESC *mvmd)
{
  INT MaxListLen;

  if (PrepareMultipleVMPtrs(mvmd))
    REP_ERR_RETURN (1);

  if (MVMD_OBJTYPES(mvmd)!=BITWISE_TYPE(NODEVEC))
    REP_ERR_RETURN (1);

  if (PrepareGetBoundaryNeighbourVectors(theGrid,&MaxListLen))
    REP_ERR_RETURN (1);

  if (MaxListLen>MAX_BND_VECTORS)
    REP_ERR_RETURN (1);

  /* index has been used by PrepareGetBoundaryNeighbourVectors */
  l_setindex(theGrid);

  return (0);
}

/****************************************************************************/
/** \brief Get list of DOUBLE pointers for vectors and matrices

   \param mvmd - data filled by PrepareBndVecMultipleVMPtrs
   \param cnt - items in the vector list
   \param VecList - list of local boundary vector neighborhood
   \param vptrlist - pointer to lists of double values corresponding the  VECDATA_DESC-list
   \param mptrlist - pointer to lists of double values corresponding the  MATDATA_DESC-list
   \param vecskip - set 1 for DIRICHLET boundary, 0 else (ordering corresponds to the first
                                VECDATA_DESC)
   \param vtype - types of vectors collected for first vd
   \param nvec - number of vectors involved from this element

   This functions returns pointers to the data fields described in a VECDATA_DESC-list
   and a MATDATA_DESC-list passed in the mvmd argument. Before call of this function
   in an element loop PrepareBndVecMultipleVMPtrs has to be called with mvmd.

   In the case that components are ordered subsequently only a pointer to the first
   DOUBLE value is returned. The pointers can be incremented by the user.

   RETURN VALUE:
   .n    0: ok
   .n    >0: error

   SEE ALSO:
   GetMultipleVMPtrs,PrepareBndVecMultipleVMPtrs
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetBndVecMultipleVMPtrs (const MVM_DESC *mvmd,
                                           INT *cnt,
                                           VECTOR *VecList[],
                                           DOUBLE **vptrlist[MAXVD],
                                           DOUBLE **mptrlist[MAXMD],
                                           INT *vecskip, INT *vtype, INT nvec[MAXVD])
{
  if (GetBoundaryNeighbourVectors(MVMD_DATATYPES(mvmd),MVMD_OBJTYPES(mvmd),cnt,VecList)!=GM_OK)
    REP_ERR_RETURN (1);

  return (GetMultipleVMPtrs(mvmd,*cnt,VecList,vptrlist,mptrlist,vecskip,vtype,nvec));
}

INT NS_DIM_PREFIX ResetBndVecMultipleVMPtrs (void)
{
  if (ResetGetBoundaryNeighbourVectors())
    REP_ERR_RETURN (1);
  return (0);
}

INT NS_DIM_PREFIX FinishBndVecMultipleVMPtrs (void)
{
  if (FinishBoundaryNeighbourVectors()!=GM_OK)
    REP_ERR_RETURN (1);
  return (0);
}

/****************************************************************************/
/** \brief Set all vecskip flags to 0

   \param theGrid - pointer to a grid
   \param theVD - type vector descriptor

   This function sets the vecskip flags for all vectors of a grid to 0.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX ClearVecskipFlags (GRID *theGrid, const VECDATA_DESC *theVD)
{
  VECTOR *theVector;
  INT j,n;

  for (theVector=PFIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
  {
    n = VD_NCMPS_IN_TYPE(theVD,VTYPE(theVector));
    for (j=0; j<n; j++)
      VECSKIP(theVector) =  CLEAR_FLAG(VECSKIP(theVector),(1<<j));
  }

  return(0);
}

/****************************************************************************/
/** \brief Set bits where skip component of vectors have to be cleared

   \param vd - descriptor of the global data
   \param vds - sub descriptor of vd describing the part to be cleared
   \param n - number of components at interface
   \param typeskip - bits are set where comps of vds relative to vd
   \param co_n - number of co-components at interface
   \param co_typeskip - bits are set where not comps of vds relative to vd

   This function sets bits where skip the component of vectors have to be cleared.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX ComputePartVecskip (const VECDATA_DESC *vd, const VECDATA_DESC *vds,
                                      INT typeskip[NVECTYPES],
                                      INT co_typeskip[NVECTYPES])
{
  INT tp,n,ns,i,j,cmp;

  for (tp=0; tp<NVECTYPES; tp++)
  {
    typeskip[tp] = co_typeskip[tp] = 0;
    if (VD_ISDEF_IN_TYPE(vds,tp))
    {
      if (!VD_ISDEF_IN_TYPE(vd,tp))
        REP_ERR_RETURN (1);

      n  = VD_NCMPS_IN_TYPE(vd,tp);
      ns = VD_NCMPS_IN_TYPE(vds,tp);
      if (ns<n)
      {
        /* set all bits of comps in vds relative to vd */
        for (i=0; i<n; i++)
        {
          cmp = VD_CMP_OF_TYPE(vd,tp,i);
          for (j=0; j<ns; j++)
            if (VD_CMP_OF_TYPE(vds,tp,j)==cmp)
              break;
          if (j<ns)
            /* cmp contained in vds */
            typeskip[tp] |= 1<<i;
          else
            co_typeskip[tp] |= 1<<i;
        }
      }
      else if (ns==n)
      {
        for (i=0; i<n; i++)
          typeskip[tp] |= 1<<i;
        co_typeskip[tp] = 0;
      }
      else
        /* vd does not contain vds */
        REP_ERR_RETURN (1);
    }
  }
  return (0);
}

/****************************************************************************/
/** \brief Set vecskip bits to 0 in part of the domain

   \param theGrid - pointer to a grid
   \param typeskip - clear flags where bits are set

   This function resets the vecskip flags in all vectors of a grid at positions
   where flags are set in typeskip.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX ClearPartVecskipFlags (GRID *theGrid, const INT typeskip[NVECTYPES])
{
  VECTOR *theVector;
  INT pattern[NVECTYPES];
  INT tp;

  for (tp=0; tp<NVECTYPES; tp++)
    pattern[tp] = ~typeskip[tp];

  /* compute skip positions of each type */
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
    VECSKIP(theVector) &= pattern[VTYPE(theVector)];

  return(0);
}

/****************************************************************************/
/** \brief Get vecskip entries

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param vecskip - returns 1 for DIRICHLET boundary, 0 else

   This function gets the vecskip flags for all local vectors
   corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetElementDirichletFlags (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                            INT *vecskip)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  INT i,j,m,cnt,vtype;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype = VTYPE(theVec[i]);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++)
      vecskip[m++] = ((VECSKIP(theVec[i]) & (1<<j))!=0);
  }

  return (m);
}

/****************************************************************************/
/** \brief Set vecskip entries

   \param theElement - pointer to an element
   \param theVD - type vector descriptor
   \param vecskip - set 1 for DIRICHLET boundary, 0 else

   This function sets the vecskip flags for all local vectors
   corresponding to an element.

   RETURN VALUE:
   .n    number of components
   .n    -1 if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX SetElementDirichletFlags (ELEMENT *theElement, const VECDATA_DESC *theVD,
                                            INT *vecskip)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  INT i,j,m,cnt,vtype;

  cnt = GetAllVectorsOfElementOfType(theElement,theVec,theVD);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype = VTYPE(theVec[i]);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++)
      if (vecskip[m++] == 1)
        VECSKIP(theVec[i]) |= (1<<j);
  }

  return (m);
}

/****************************************************************************/
/** \brief Sets values with Dirichlet flag to 0

   \param theGrid - pointer to a grid
   \param x - type vector descriptor

   This function sets all components of x where the vecskip flag is set
   to zero.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX ClearDirichletValues (GRID *theGrid, VECDATA_DESC *x)
{
  VECTOR *theVector;
  INT j,type,ncomp;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL;
       theVector=SUCCVC(theVector)) {
    type = VTYPE(theVector);
    ncomp = VD_NCMPS_IN_TYPE (x,type);
    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(theVector) & (1<<j))
        VVALUE(theVector,VD_CMP_OF_TYPE(x,type,j)) = 0.0;
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Modifies matrix entries and right hand side
   for Dirichlet boundary

   \param theGrid - pointer to a grid
   \param Mat - type matrix descriptor for the stiffness matix
   \param Sol - type vector descriptor for the solution vector
   \param Rhs - type vector descriptor for the right hand side

   This function sets all components of Rhs where the vecskip flag is set
   to the Sol value and sets the corresponding matrix row to the unit vector.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX AssembleDirichletBoundary (GRID *theGrid, const MATDATA_DESC *Mat,
                                             const VECDATA_DESC *Sol, const VECDATA_DESC *Rhs)
{
  VECTOR *theVector;
  MATRIX *theMatrix;
  INT i,j,comp1,comp2,ncomp,dcomp,type,dtype,cnt;

  cnt = 0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL;
       theVector=SUCCVC(theVector))
  {
    type = VTYPE(theVector);
    ncomp = VD_NCMPS_IN_TYPE (Sol,type);
    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(theVector) & (1<<j))
      {
        cnt++;
        comp1 = VD_CMP_OF_TYPE(Sol,type,j);
        comp2 = VD_CMP_OF_TYPE(Rhs,type,j);
        VVALUE(theVector,comp2) =     VVALUE(theVector,comp1);
        theMatrix = VSTART(theVector);
        for (i=j*ncomp; i<(j+1)*ncomp; i++)
          MVALUE(theMatrix,MD_MCMP_OF_RT_CT(Mat,type,type,i)) = 0.0;
        MVALUE(theMatrix,MD_MCMP_OF_RT_CT(Mat,type,type,j+j*ncomp))=1.0;
        for (theMatrix=MNEXT(theMatrix); theMatrix!=NULL;
             theMatrix=MNEXT(theMatrix))
        {
          dtype = MDESTTYPE(theMatrix);
          dcomp = VD_NCMPS_IN_TYPE (Sol,dtype);
          if (dcomp == 0) continue;
          for (i=j*dcomp; i<(j+1)*dcomp; i++)
            MVALUE(theMatrix,MD_MCMP_OF_RT_CT(Mat,type,dtype,i)) = 0.0;
        }
      }
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Modifies matrix entries for Dirichlet values
   for Dirichlet boundary

   \param theGrid - pointer to a grid
   \param Mat - type matrix descriptor for the stiffness matix

   This function sets all matrix rows to the unit vector where the vecskip flag is set.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX ModifyDirichletMatrix (GRID *theGrid, const MATDATA_DESC *Mat)
{
  VECTOR *theVector;
  MATRIX *theMatrix;
  INT i,j,ncomp,dcomp,type,dtype;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
  {
    type = VTYPE(theVector);
    ncomp = MD_ROWS_IN_RT_CT (Mat,type,type);
    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(theVector) & (1<<j))
      {
        theMatrix = VSTART(theVector);
        for (i=j*ncomp; i<(j+1)*ncomp; i++)
          MVALUE(theMatrix,MD_MCMP_OF_RT_CT(Mat,type,type,i)) = 0.0;
        MVALUE(theMatrix,MD_MCMP_OF_RT_CT(Mat,type,type,j+j*ncomp)) = 1.0;

        for (theMatrix=MNEXT(theMatrix); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
        {
          dtype = MDESTTYPE(theMatrix);
          dcomp = MD_COLS_IN_RT_CT (Mat,type,dtype);
          if (dcomp == 0) continue;
          for (i=j*dcomp; i<(j+1)*dcomp; i++)
            MVALUE(theMatrix,MD_MCMP_OF_RT_CT(Mat,type,dtype,i)) = 0.0;
        }
      }
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Set defect to zero for Dirichlet boundary

   \param theGrid - pointer to a grid
   \param Cor - type vector descriptor for the correction vector

   This function sets all components of Cor to zero where the vecskip flag is set.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX ModifyDirichletDefect (GRID *theGrid, const VECDATA_DESC *Cor)
{
  VECTOR *theVec;
  INT j,comp1,ncomp,type;

  for (theVec=FIRSTVECTOR(theGrid); theVec!= NULL; theVec=SUCCVC(theVec))
  {
    type = VTYPE(theVec);
    ncomp = VD_NCMPS_IN_TYPE (Cor,type);
    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(theVec) & (1<<j))
      {
        comp1 = VD_CMP_OF_TYPE(Cor,type,j);
        VVALUE(theVec,comp1) = 0.;
      }
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Modifies matrix entries and right hand side
   for Dirichlet boundary

   \param theGrid - pointer to a grid
   \param Mat - type matrix descriptor for the stiffness matix
   \param Sol - type vector descriptor for the solution vector
   \param Rhs - type vector descriptor for the right hand side

   This function sets all components of Rhs where the vecskip flag is set
   to the Sol value and sets the corresponding matrix row and column to the unit vector.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX AssembleTotalDirichletBoundary (GRID *theGrid, const MATDATA_DESC *Mat,
                                                  const VECDATA_DESC *Sol,
                                                  const VECDATA_DESC *Rhs)
{
  VECTOR *theVector;
  INT i,j;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL;
       theVector=SUCCVC(theVector))
  {
    INT type = VTYPE(theVector);
    INT ncomp = VD_NCMPS_IN_TYPE (Sol,type);

    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(theVector) & (1<<j))
      {
        MATRIX *theMatrix = VSTART(theVector);
        INT comp1 = VD_CMP_OF_TYPE(Sol,type,j);
        DOUBLE s = VVALUE(theVector,comp1);

        VVALUE(theVector,VD_CMP_OF_TYPE(Rhs,type,j)) = 0.0;
        for (i=0; i<ncomp; i++)
          if ((i != j) && (!(VECSKIP(theVector) & (1<<i))))
            VVALUE(theVector,VD_CMP_OF_TYPE(Rhs,type,i)) -=
              s * MVALUE(theMatrix,
                         MD_MCMP_OF_RT_CT(Mat,type,type,i*ncomp+j));
        for (i=0; i<ncomp; i++) {
          MVALUE(theMatrix,
                 MD_MCMP_OF_RT_CT(Mat,type,type,i*ncomp+j)) = 0.0;
          MVALUE(theMatrix,
                 MD_MCMP_OF_RT_CT(Mat,type,type,j*ncomp+i)) = 0.0;
        }
        MVALUE(theMatrix,
               MD_MCMP_OF_RT_CT(Mat,type,type,j+j*ncomp)) = 1.0;
        for (theMatrix=MNEXT(theMatrix); theMatrix!=NULL;
             theMatrix=MNEXT(theMatrix))
        {
          VECTOR *dest = MDEST(theMatrix);
          INT dtype = MDESTTYPE(theMatrix);
          INT dcomp = VD_NCMPS_IN_TYPE (Sol,dtype);

          if (dcomp == 0) continue;
          for (i=0; i<dcomp; i++) {
            if (!(VECSKIP(dest) & (1<<i)))
              VVALUE(dest,VD_CMP_OF_TYPE(Rhs,dtype,i)) -=
                s * MVALUE(MADJ(theMatrix),
                           MD_MCMP_OF_RT_CT(Mat,dtype,type,i*ncomp+j));
            MVALUE(theMatrix,
                   MD_MCMP_OF_RT_CT(Mat,type,dtype,j*dcomp+i))
              = 0.0;
            MVALUE(MADJ(theMatrix),
                   MD_MCMP_OF_RT_CT(Mat,dtype,type,i*ncomp+j))
              = 0.0;
          }
        }
      }
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Converts a matrix in a sparse format

   \param theGrid - pointer to a grid
   \param theHeap - pointer to heap
   \param MarkKey - mark in temp mem
   \param A - pointer to matrix descriptor
   \param symmetric - if symmetric, store only upper part
   \param pn - number of lines
   \param pia - diagonal indices
   \param pja - row indices
   \param pa - matrix entries

   This function converts a matrix in a sparse format.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX ConvertMatrix (GRID *theGrid, HEAP *theHeap, INT MarkKey,
                                 MATDATA_DESC *A, INT symmetric,
                                 int *pn, int **pia, int **pja, double **pa)
{
  VECTOR *v;
  MATRIX *m;
  INT rtype,ctype,rcomp,ccomp,i,j,k,n,nn;
  int *ia,*ja;
  double *a;
  SHORT *Mcomp;

  n = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(A,rtype,rtype);
    VINDEX(v) = n;
    n += rcomp;
  }

  n = nn = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(A,rtype,rtype);
    for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
      if (symmetric)
        if (VINDEX(MDEST(m)) > n) continue;
      ctype = MDESTTYPE(m);
      ccomp = MD_COLS_IN_RT_CT(A,rtype,ctype);
      if (ccomp == 0) continue;
      nn += rcomp * ccomp;
    }
    n += rcomp;
  }
  ia = (int *)GetTmpMem(theHeap,sizeof(int) * (n+1),MarkKey);
  a = (double *)GetTmpMem(theHeap,sizeof(double) * nn,MarkKey);
  ja = (int *)GetTmpMem(theHeap,sizeof(int) * nn,MarkKey);
  if (ia == NULL || a == NULL || ja == NULL)
    return(NUM_ERROR);
  n = nn = 0;
  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(A,rtype,rtype);
    for (i=0; i<rcomp; i++) {
      ia[n++] = nn;
      for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
        k = VINDEX(MDEST(m));
        ctype = MDESTTYPE(m);
        ccomp = MD_COLS_IN_RT_CT(A,rtype,ctype);
        if (ccomp == 0) continue;
        Mcomp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        for (j=0; j<ccomp; j++) {
          if (symmetric)
            if (k >= n) continue;
          a[nn] = MVALUE(m,Mcomp[i*ccomp+j]);
          ja[nn] = k++;
          nn++;
        }
      }
    }
  }
  ia[n] = nn;
  pn[0] = n;
  pia[0] = ia;
  pja[0] = ja;
  pa[0] = a;

  return(NUM_OK);
}


static void PrintSingleVectorX (const VECTOR *v, const VECDATA_DESC *X, INT vclass, INT vnclass, PrintfProcPtr Printf, INT *info)
{
  char buffer[256];
  DOUBLE_VECTOR pos;
  INT comp,ncomp,i,j;

  if (VCLASS(v) > vclass) return;
  if (VNCLASS(v) > vnclass) return;
  ncomp = VD_NCMPS_IN_TYPE(X,VTYPE(v));
  if (ncomp == 0) return;
  /* Check if there is an object associated with the vector. */
  i = 0;
  if (VOBJECT(v) != NULL) {
    VectorPosition(v,pos);
    i += sprintf(buffer,"x=%5.2f y=%5.2f ",pos[0],pos[1]);
    if (DIM == 3)
      i += sprintf(buffer+i,"z=%5.2f ",pos[2]);
  }
  else {
    *info = TRUE;
    i += sprintf(buffer,"                ");
    if (DIM == 3)
      i += sprintf(buffer+i,"        ");
  }
  for (j=0; j<ncomp; j++)
  {
    comp = VD_CMP_OF_TYPE(X,VTYPE(v),j);
    i += sprintf(buffer+i,"u[%d]=%15.8f ",j,VVALUE(v,comp));
  }
  i += sprintf(buffer+i,"   cl %d %d sk ",VCLASS(v),VNCLASS(v));
  for (j=0; j<ncomp; j++)
    i += sprintf(buffer+i,"%d ",((VECSKIP(v) & (1<<j))!=0));
  i += sprintf(buffer+i,"n %d t %d o %d\n",VNEW(v),VTYPE(v),VOTYPE(v));
  Printf(buffer);

        #ifdef Debug
  if (Printf!=PrintDebug)
    PRINTDEBUG(np,1,("%d: %s",me,buffer));
        #endif

  return;
}

INT NS_DIM_PREFIX PrintVectorListX (const VECTOR *vlist[], const VECDATA_DESC *X, INT vclass, INT vnclass, PrintfProcPtr Printf)
{
  INT info=FALSE;

  for (; *vlist!= NULL; ++vlist)
    PrintSingleVectorX(*vlist,X,vclass,vnclass,Printf,&info);

  if (info) Printf("NOTE: Geometrical information not available for some vectors.\n");

  return(NUM_OK);
}

INT NS_DIM_PREFIX PrintVectorX (const GRID *g, const VECDATA_DESC *X, INT vclass, INT vnclass, PrintfProcPtr Printf)
{
  const VECTOR *v;
  INT info=FALSE;

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
    PrintSingleVectorX(v,X,vclass,vnclass,Printf,&info);

  if (info) Printf("NOTE: Geometrical information not available for some vectors.\n");

  return(NUM_OK);
}

/****************************************************************************/
/** \brief Print a vector list

   \param g - pointer to a grid
   \param X - pointer to vector descriptor
   \param vclass - class number
   \param vnclass - class number

   This function prints the values of X for all vectors with class smaller
   or equal to vclass and next class smaller or equal to vnclass.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX PrintVector (GRID *g, const VECDATA_DESC *X, INT vclass, INT vnclass)
{
  return (PrintVectorX(g,X,vclass,vnclass,UserWriteF));
}

INT NS_DIM_PREFIX PrintSVector (MULTIGRID *mg, VECDATA_DESC *X)
{
  VECTOR *v;
  DOUBLE_VECTOR pos;
  INT vtype,comp,ncomp,j,lev;

  for (vtype=0; vtype<NVECTYPES; vtype++) {
    ncomp = VD_NCMPS_IN_TYPE(X,vtype);
    if (ncomp == 0) continue;
    comp = VD_CMP_OF_TYPE(X,vtype,0);
    S_BELOW_VLOOP__TYPE(lev,0,CURRENTLEVEL(mg),v,mg,vtype) {
      VectorPosition(v,pos);
      UserWriteF("x=%5.2f y=%5.2f ",pos[0],pos[1]);
      if (DIM == 3)
        UserWriteF("z=%5.2f ",pos[2]);
      for (j=0; j<ncomp; j++)
        UserWriteF("u[%d]=%15.8f ",j,VVALUE(v,comp+j));
      UserWriteF("   cl %d %d sk ",VCLASS(v),VNCLASS(v));
      for (j=0; j<ncomp; j++)
        UserWriteF("%d ",((VECSKIP(v) & (1<<j))!=0));
      UserWriteF("\n");
    }
    S_FINE_VLOOP__TYPE(CURRENTLEVEL(mg),v,mg,vtype) {
      VectorPosition(v,pos);
      UserWriteF("x=%5.2f y=%5.2f ",pos[0],pos[1]);
      if (DIM == 3)
        UserWriteF("z=%5.2f ",pos[2]);
      for (j=0; j<ncomp; j++)
        UserWriteF("u[%d]=%15.8f ",j,VVALUE(v,comp+j));
      UserWriteF("   cl %d %d sk ",VCLASS(v),VNCLASS(v));
      for (j=0; j<ncomp; j++)
        UserWriteF("%d ",((VECSKIP(v) & (1<<j))!=0));
      UserWriteF("\n");
    }
  }

  return(NUM_OK);
}

/****************************************************************************/
/** \brief Print a matrix list

   \param g - pointer to a grid
   \param M - pointer to matrix descriptor
   \param vclass - class number
   \param vnclass - class number

   This function prints the values of M of the matrix list
   for all vectors with class smaller
   or equal to vclass and next class smaller or equal to vnclass.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX PrintMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass)
{
  VECTOR *v;
  MATRIX *m;
  INT Mcomp,rcomp,ccomp,i,j,rtype,ctype;

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    if (VCLASS(v) > vclass) continue;
    if (VNCLASS(v) > vnclass) continue;
    rtype = VTYPE(v);
    rcomp = MD_ROWS_IN_RT_CT(Mat,rtype,rtype);
    for (i=0; i<rcomp; i++)
    {
      for (m=VSTART(v); m!=NULL; m = MNEXT(m))
      {
        ctype = MDESTTYPE(m);
        ccomp = MD_COLS_IN_RT_CT(Mat,rtype,ctype);
        if (ccomp == 0) continue;
        if (rcomp != MD_ROWS_IN_RT_CT(Mat,rtype,ctype))
          UserWrite("wrong type\n");
        Mcomp = MD_MCMP_OF_RT_CT(Mat,rtype,ctype,i*ccomp);
        for (j=0; j<ccomp; j++)
          UserWriteF("%16.8e ",MVALUE(m,Mcomp+j));
      }
      UserWrite("\n");
    }
  }

  return(NUM_OK);
}

INT NS_DIM_PREFIX PrintTMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass)
{
  VECTOR *v;
  MATRIX *m;
  INT Mcomp,rcomp,ccomp,i,j,rtype,ctype;

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    if (VCLASS(v) > vclass) continue;
    if (VNCLASS(v) > vnclass) continue;
    rtype = VTYPE(v);
    ccomp = MD_COLS_IN_RT_CT(Mat,rtype,rtype);
    for (i=0; i<ccomp; i++)
    {
      for (m=VSTART(v); m!=NULL; m = MNEXT(m))
      {
        ctype = MDESTTYPE(MADJ(m));
        rcomp = MD_ROWS_IN_RT_CT(Mat,rtype,ctype);
        Mcomp = MD_MCMP_OF_RT_CT(Mat,rtype,ctype,0);
        for (j=0; j<rcomp; j++)
          UserWriteF("%4.2f ",MVALUE(MADJ(m),Mcomp+j*ccomp+i));
      }
      UserWrite("\n");
    }
  }

  return(NUM_OK);
}

INT NS_DIM_PREFIX PrintDiagMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass)
{
  char buffer[256];
  VECTOR *v;
  DOUBLE_VECTOR pos;
  INT info=FALSE;
  MATRIX *m;
  INT Mcomp,ccomp,i,j,rtype;

  for (v=PFIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    if (VCLASS(v) > vclass) continue;
    if (VNCLASS(v) > vnclass) continue;
    rtype = VTYPE(v);
    ccomp = MD_COLS_IN_RT_CT(Mat,rtype,rtype);
    if (ccomp == 0) continue;
    m = VSTART(v);
    Mcomp = MD_MCMP_OF_RT_CT(Mat,rtype,rtype,0);
    /* Check if there is an object associated with the vector. */
    i = 0;
    if (VOBJECT(v) != NULL) {
      VectorPosition(v,pos);
      i += sprintf(buffer,"x=%5.2f y=%5.2f ",pos[0],pos[1]);
      if (DIM == 3)
        i += sprintf(buffer+i,"z=%5.2f ",pos[2]);
#ifdef ModelP
      i += sprintf(buffer+i,"l %d p %d ",GLEVEL(g),PRIO(v));
#endif
    }
    else {
      info = TRUE;
      i += sprintf(buffer,"                ");
      if (DIM == 3)
        i += sprintf(buffer+i,"        ");
#ifdef ModelP
      i += sprintf(buffer+i,"l %d p %d ",GLEVEL(g),-1);
#endif
    }
    for (j=0; j<ccomp; j++)
      i += sprintf(buffer+i,"d[%d]=%15.8f ",j,
                   MVALUE(m,Mcomp+j*ccomp+j));
    i += sprintf(buffer+i,"\n");
    UserWrite(buffer);

    PRINTDEBUG(np,1,("%d: %s",me,buffer));
  }


  if (info)
    UserWrite("NOTE: "
              "Geometrical information not available for some vectors.\n");

  return(NUM_OK);
}

#ifdef __INTERPOLATION_MATRIX__
/****************************************************************************/
/** \brief Print the interpolation matrix list

   \param g - pointer to a grid
   \param V - pointer to vector descriptor
   \param vclass - class number
   \param vnclass - class number

   This function prints the values of the interpolation matrix list
   for all vectors with class smaller
   or equal to vclass and next class smaller or equal to vnclass.

   RETURN VALUE:
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX PrintIMatrix (GRID *g, VECDATA_DESC *V, INT vclass, INT vnclass)
{
  VECTOR *v;
  MATRIX *m;
  INT Mcomp,rcomp,ccomp,i,j;

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    if (VCLASS(v) > vclass) continue;
    if (VNCLASS(v) > vnclass) continue;
    rcomp = VD_NCMPS_IN_TYPE(V,VTYPE(v));
    for (i=0; i<rcomp; i++)
    {
      for (m=VISTART(v); m!=NULL; m = MNEXT(m))
      {
        ccomp = VD_NCMPS_IN_TYPE(V,MDESTTYPE(m));
        Mcomp = i * ccomp;
        for (j=0; j<ccomp; j++)
          UserWriteF("%+5.3f ",MVALUE(m,Mcomp+j));
      }
      UserWrite("\n");
    }
  }

  return(NUM_OK);
}
#endif  /* __INTERPOLATION_MATRIX__ */
