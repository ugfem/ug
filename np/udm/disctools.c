// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      disctools.c                                                   */
/*                                                                          */
/* Purpose:   tools for assembling                                          */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Nov 27 95 begin                                                                           */
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "gm.h"       /* for data structure               */
#include "evm.h"      /* for data structure               */
#include "general.h"
#include "devices.h"

#include "ugblas.h"
#include "disctools.h"

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
/*D
   MG_GetCoeffFct - get function pointer

   SYNOPSIS:
   CoeffProcPtr MG_GetCoeffFct (MULTIGRID *theMG, INT n);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  n - number of coefficient function

   DESCRIPTION:
   This function returns a pointer to the nth coefficient function
   of the multigrid.

   RETURN VALUE:
   CoeffProcPtr
   .n    pointer to function
   .n    NULL if n is too large
   D*/
/****************************************************************************/

CoeffProcPtr MG_GetCoeffFct (MULTIGRID *theMG, INT n)
{
  BVP *myBVP;
  BVP_DESC *BVP_desc;
  CoeffProcPtr cpp;

  myBVP = MG_BVP(theMG);
  BVP_desc = MG_BVPD(theMG);
  cpp = NULL;
  if ((n >= 0) && (n < BVPD_NCOEFFF(BVP_desc)))
    BVP_SetCoeffFct       (myBVP,n,&cpp);
  return(cpp);
}

/****************************************************************************/
/*D
   MG_GetUserFct - get function pointer

   SYNOPSIS:
   UserProcPtr MG_GetUserFct (MULTIGRID *theMG, INT n);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  n - number of a user function

   DESCRIPTION:
   This function returns a pointer to the nth user function
   of the multigrid.

   RETURN VALUE:
   CoeffProcPtr
   .n    pointer to function
   .n    NULL if n is too large
   D*/
/****************************************************************************/

UserProcPtr MG_GetUserFct (MULTIGRID *theMG, INT n)
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
/*D
   GetAllVectorsOfElementOfType - get vector list

   SYNOPSIS:
   INT GetAllVectorsOfElementOfType (ELEMENT *theElement, VECTOR **vec,
   VECDATA_DESC *theVD);

   PARAMETERS:
   .  theElement - pointer to an element
   .  vec - vector list
   .  theVD - vector descriptor

   DESCRIPTION:
   This function gets a list of vectors corresponding to an element.
   It uses GetVectorsOfDataTypesInObjects (which should be preferred).

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetAllVectorsOfElementOfType (ELEMENT *theElement, VECTOR **vec,
                                  const VECDATA_DESC *theVD)
{
  INT cnt;

  if (GetVectorsOfDataTypesInObjects(theElement,VD_DATA_TYPES(theVD),VD_OBJ_USED(theVD),&cnt,vec))
    return (-1);

  return (cnt);
}

INT GetAllVectorsOfElementsideOfType (ELEMENT *theElement, INT side,
                                      VECTOR **vec,
                                      const VECDATA_DESC *theVD)
{
  VECTOR *v[MAX_NODAL_VECTORS];
  INT i,cnt,cnt1;

  if (GetVectorsOfDataTypesInObjects(theElement,VD_DATA_TYPES(theVD),
                                     VD_OBJ_USED(theVD),&cnt,v))
    return (-1);

  cnt = cnt1 = 0;
  if (VD_NCMPS_IN_TYPE(theVD,NODEVEC)) {
    for (i=0; i<CORNERS_OF_SIDE(theElement,side); i++)
      vec[cnt++] = v[CORNER_OF_SIDE(theElement,side,i)];
    cnt1 += CORNERS_OF_ELEM(theElement);
  }
  if (VD_NCMPS_IN_TYPE(theVD,EDGEVEC)) {
    for (i=0; i<EDGES_OF_SIDE(theElement,side); i++)
      vec[cnt++] = v[cnt1+EDGE_OF_SIDE(theElement,side,i)];
    cnt1 += EDGES_OF_ELEM(theElement);
  }
    #ifdef __THREEDIM__
  if (VD_NCMPS_IN_TYPE(theVD,SIDEVEC))
    vec[cnt++] = v[cnt1];
    #endif

  return (cnt);
}

/****************************************************************************/
/*D
   GetElementsideIndices - compute vector indices of an element side

   SYNOPSIS:
   INT GetElementsideIndices (ELEMENT *theElement, INT side,
   VECDATA_DESC *theVD, INT *index);

   PARAMETERS:
   .  theElement - pointer to an element
   .  side - element side
   .  theVD - type vector descriptor
   .  index - index list

   DESCRIPTION:
   This function gets the indices of the vector pointers of an element
   corresponding to a side.

   RETURN VALUE:
   INT
   .n    number of components of the element side
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetElementsideIndices (ELEMENT *theElement, INT side,
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
/*D
   GetElementVPtrs - get list of DOUBLE pointers for vectors

   SYNOPSIS:
   INT GetElementVPtrs (ELEMENT *theElement, VECDATA_DESC *theVD,
   DOUBLE **vptr);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  vptr - pointer to double values

   DESCRIPTION:
   This function gets all local vector pointers corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetElementVPtrs (ELEMENT *theElement, const VECDATA_DESC *theVD, DOUBLE **vptr)
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
/*D
   GetElementVValue - get list of DOUBLE values for vectors

   SYNOPSIS:
   INT GetElementVValues (ELEMENT *theElement, VECDATA_DESC *theVD,
   DOUBLE *value);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  value - pointer to double values

   DESCRIPTION:
   This function gets all local vector values corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetElementVValues (ELEMENT *theElement, const VECDATA_DESC *theVD,
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
/*D
   AddElementVValue - add list of DOUBLE values for vectors

   SYNOPSIS:
   INT AddElementVValues (ELEMENT *theElement, VECDATA_DESC *theVD,
   DOUBLE *value);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  value - pointer to double values

   DESCRIPTION:
   This function adds all local vector values corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT AddElementVValues (ELEMENT *theElement, const VECDATA_DESC *theVD,
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
/*D
   GetVlistVValue - get list of DOUBLE values for vectors

   SYNOPSIS:
   INT GetVlistVValues (INT cnt, VECTOR **theVec,
   const VECDATA_DESC *theVD, DOUBLE *value);

   PARAMETERS:
   .  cnt - number of vectors
   .  theVec - vector list
   .  theVD - type vector descriptor
   .  value - pointer to double values

   DESCRIPTION:
   This function gets all local vector values corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetVlistVValues (INT cnt, VECTOR **theVec,
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
/*D
   AddVlistVValues - get list of DOUBLE values for vectors

   SYNOPSIS:
   INT AddVlistVValues (INT cnt, VECTOR **theVec,
   const VECDATA_DESC *theVD, DOUBLE *value);

   PARAMETERS:
   .  cnt - number of vectors
   .  theVec - vector list
   .  theVD - type vector descriptor
   .  value - pointer to double values

   DESCRIPTION:
   This function adds all local vector values corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT AddVlistVValues (INT cnt, VECTOR **theVec,
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
/*D
   SetVlistVValues - set list of DOUBLE values for vectors

   SYNOPSIS:
   INT SetVlistVValues (INT cnt, VECTOR **theVec,
   const VECDATA_DESC *theVD, DOUBLE *value);

   PARAMETERS:
   .  cnt - number of vectors
   .  theVec - vector list
   .  theVD - type vector descriptor
   .  value - pointer to double values

   DESCRIPTION:
   This function sets all local vector values corresponding to a vector list.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT SetVlistVValues (INT cnt, VECTOR **theVec,
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
/*D
   GetElementVPtrsVecskip - get list of DOUBLE pointers for vectors

   SYNOPSIS:
   INT GetElementVPtrsVecskip (ELEMENT *theElement, VECDATA_DESC *theVD,
   DOUBLE **vptr, INT *vecskip);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  vptr - pointer to double values
   .  vecskip - set 1 for DIRICHLET boundary, 0 else

   DESCRIPTION:
   This function gets all local vector pointers corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetElementVPtrsVecskip (ELEMENT *theElement, const VECDATA_DESC *theVD,
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
/*D
   GetElementNewVPtrs - get list of DOUBLE pointers for vectors

   SYNOPSIS:
   INT GetElementNewVPtrs (ELEMENT *theElement, const VECDATA_DESC *theVD,
   DOUBLE **vptr, INT *newField);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  vptr - pointer to double values
   .  newField - set 1 for new vectors, 0 else

   DESCRIPTION:
   This function gets all local vector pointers corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    0 if there is no new vector
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetElementNewVPtrs (ELEMENT *theElement, const VECDATA_DESC *theVD,
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
/*D
   GetElementMPtrs - get list of DOUBLE pointers for matrices

   SYNOPSIS:
   INT GetElementMPtrs (ELEMENT *theElement, const MATDATA_DESC *md,
   DOUBLE **mptr);

   PARAMETERS:
   .  theElement - pointer to an element
   .  md - matrix data descriptor
   .  mptr - pointer to double values corresponding to the local stiffness matrix

   DESCRIPTION:
   This function gets all local matrix pointers corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of dofs
   D*/
/****************************************************************************/

INT GetElementMPtrs (ELEMENT *theElement, const MATDATA_DESC *md,
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
/*D
   GetVlistMValues - get list of DOUBLE values for matrices

   SYNOPSIS:
   INT GetVlistMValues (INT cnt, VECTOR **theVec,
   const MATDATA_DESC *theMD, DOUBLE *value);

   PARAMETERS:
   .  cnt - number of vectors
   .  theVec - vector list
   .  theMD - type matrix descriptor
   .  value - pointer to double values

   DESCRIPTION:
   This function get all local matrix values corresponding to a vector list.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetVlistMValues (INT cnt, VECTOR **theVec,
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
/*D
   AddVlistMValues - add list of DOUBLE values for matrices

   SYNOPSIS:
   INT AddVlistMValues (GRID *theGrid, INT cnt, VECTOR **theVec,
   const MATDATA_DESC *theMD, DOUBLE *value);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  cnt - number of vectors
   .  theVec - vector list
   .  theMD - type matrix descriptor
   .  value - pointer to double values

   DESCRIPTION:
   This function adds all local matrix values corresponding to a vector list.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT AddVlistMValues (GRID *theGrid, INT cnt, VECTOR **theVec,
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
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[comp[l*vncomp[i]+k]] += value[(m2+l)*m+m1+k];
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

/****************************************************************************/
/*D
   GetElementVMPtrs - get list of DOUBLE pointers for vectors and matrices

   SYNOPSIS:
   INT GetElementVMPtrs (ELEMENT *theElement,
   VECDATA_DESC *theVD, MATDATA_DESC *md,
   DOUBLE **vptr, DOUBLE **mptr);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  md - type matrix descriptor
   .  vptr - pointer to double values corresponding to the local right hand side
   .  mptr - pointer to double values corresponding to the local stiffness matrix

   DESCRIPTION:
   This function gets all local vector pointers corresponding to an element
   and the connecting matrix pointers.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetElementVMPtrs (ELEMENT *theElement,
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
/*D
   GetElementVVMPtrs - get list of DOUBLE pointers for vectors and matrices

   SYNOPSIS:
   INT GetElementVVMPtrs (ELEMENT *theElement, VECDATA_DESC *theVD1,
   VECDATA_DESC *theVD2, MATDATA_DESC *md,
   DOUBLE **vptr1, DOUBLE **vptr2, DOUBLE **mptr, INT *vecskip);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD1 - type vector descriptor
   .  theVD2 - type vector descriptor
   .  md - type matrix descriptor
   .  vptr1 - pointer to double values corresponding to the local right hand side
   .  vptr2 - pointer to double values corresponding to the local right hand side
   .  mptr - pointer to double values corresponding to the local stiffness matrix
   .  vecskip - set 1 for DIRICHLET boundary, 0 else

   DESCRIPTION:
   This function gets all local vector pointers corresponding to an element
   and the connecting matrix pointers.

   RETURN VALUE:
   INT
   .n    total number of components if ok
   .n    -1: error in GetAllVectorsOfElementOfType
   .n    -2: vecdata descriptors of different size
   .n    -3: could not get matrix
   D*/
/****************************************************************************/

INT GetElementVVMPtrs (ELEMENT *theElement, const VECDATA_DESC *theVD1,
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
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = VD_NCMPS_IN_TYPE (theVD1,vtype[i]);
    if (vncomp[i] != VD_NCMPS_IN_TYPE (theVD2,vtype[i]))
      return (-2);
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
        return (-3);
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
/*D
   PrepareMultipleVMPtrs - prepare execution of GetMultipleVMPtrs

   SYNOPSIS:
   INT PrepareMultipleVMPtrs (MVM_DESC *mvmd)

   PARAMETERS:
   .  mvmd - multiple VM data structure, partially to be filled before call

   DESCRIPTION:
   This function prepares the execution of GetMultipleVMPtrs in a
   loop. It has to be called once before running the loop. The types needed are
   evaluated and it is checked, whether the components described by the
   XXXDATA_DESCs passed are in subsequent order. In this case ONLY a pointer to
   the first component per type is returned. The pointer can be incremented by
   the user.

   RETURN VALUE:
   INT
   .n    0

   SEE ALSO:
   GetMultipleVMPtrs
   D*/
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
/*D
   GetMultipleVMPtrs - get list of DOUBLE pointers for vectors and matrices

   SYNOPSIS:
   INT GetMultipleVMPtrs (const MVM_DESC *mvmd, INT cnt, VECTOR *VecList[],
                                                                                         DOUBLE **vptrlist[MAXVD],
                                                                                         DOUBLE **mptrlist[MAXMD],
                                                                                         INT *vecskip, INT *vtype, INT nvec[MAXVD])


   PARAMETERS:
   input variables:~
   .  mvmd - data filled by PrepareElementMultipleVMPtrs
   .  cnt - number of vectors to extract pointers from
   .  VecList - list of vectors to extract pointers from

   output variables:~
   .  vptrlist - pointer to lists of double values corresponding the  VECDATA_DESC-list
   .  mptrlist - pointer to lists of double values corresponding the  MATDATA_DESC-list
   .  vecskip - set 1 for DIRICHLET boundary, 0 else (ordering corresponds to the first
                                VECDATA_DESC)
   .  vtype - types of vectors collected for first vd
   .  nvec - number of vectors involved from per vd of mvmd

   DESCRIPTION:
   This functions returns pointers to the data fields described in a VECDATA_DESC-list
   and a MATDATA_DESC-list passed in the mvmd argument. Before call of this function
   in a loop PrepareMultipleVMPtrs has to be called with mvmd.

   In the case that components are ordered subsequently only a pointer to the first
   DOUBLE value is returned. The pointers can be incremented by the user.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0: error

   SEE ALSO:
   PrepareMultipleVMPtrs
   D*/
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
/*D
   PrepareElementMultipleVMPtrs - prepare execution of GetElementMultipleVMPtrs

   SYNOPSIS:
   INT PrepareElementMultipleVMPtrs (MVM_DESC *mvmd)

   PARAMETERS:
   .  mvmd - multiple VM data structure, partially to be filled before call

   DESCRIPTION:
   This function prepares the execution of GetElementMultipleVMPtrs in an element
   loop. It has to be called once before running the loop.
   For further description see 'PrepareMultipleVMPtrs'.

   RETURN VALUE:
   INT
   .n    0

   SEE ALSO:
   PrepareMultipleVMPtrs, GetElementMultipleVMPtrs
   D*/
/****************************************************************************/

INT PrepareElementMultipleVMPtrs (MVM_DESC *mvmd)
{
  return (PrepareMultipleVMPtrs(mvmd));
}

/****************************************************************************/
/*D
   GetElementMultipleVMPtrs - get list of DOUBLE pointers for vectors and matrices

   SYNOPSIS:
   INT GetElementMultipleVMPtrs (ELEMENT *elem, const MVM_DESC *mvmd,
                                                                                         DOUBLE **vptrlist[MAXVD],
                                                                                         DOUBLE **mptrlist[MAXMD],
                                                                                         INT *vecskip, INT *vtype, INT nvec[MAXVD])


   PARAMETERS:
   .  elem - pointer to an element
   .  mvmd - data filled by PrepareElementMultipleVMPtrs
   .  vptrlist - pointer to lists of double values corresponding the  VECDATA_DESC-list
   .  mptrlist - pointer to lists of double values corresponding the  MATDATA_DESC-list
   .  vecskip - set 1 for DIRICHLET boundary, 0 else (ordering corresponds to the first
                                VECDATA_DESC)
   .  vtype - types of vectors collected for first vd
   .  nvec - number of vectors involved from this element

   DESCRIPTION:
   This functions returns pointers to the data fields described in a VECDATA_DESC-list
   and a MATDATA_DESC-list passed in the mvmd argument. Before call of this function
   in an element loop PrepareElementMultipleVMPtrs has to be called with mvmd.

   In the case that components are ordered subsequently only a pointer to the first
   DOUBLE value is returned. The pointers can be incremented by the user.

   RETURN VALUE:
   INT
   .n    0: ok
   .n    >0: error

   SEE ALSO:
   GetMultipleVMPtrs,PrepareElementMultipleVMPtrs
   D*/
/****************************************************************************/

INT GetElementMultipleVMPtrs (ELEMENT *elem, const MVM_DESC *mvmd,
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
/*D
   PrepareBndVecMultipleVMPtrs - prepare execution of GetBndVecMultipleVMPtrs

   SYNOPSIS:
   INT PrepareBndVecMultipleVMPtrs (GRID *theGrid, MVM_DESC *mvmd)

   PARAMETERS:
   .  theGrid - grid level
   .  mvmd - multiple VM data structure, partially to be filled before call

   DESCRIPTION:
   This function prepares the execution of GetBndVecMultipleVMPtrs in an element
   loop. It has to be called once before running the loop.
   For further description see 'PrepareMultipleVMPtrs'.

   RETURN VALUE:
   INT
   .n    0

   SEE ALSO:
   PrepareMultipleVMPtrs, GetBndVecMultipleVMPtrs
   D*/
/****************************************************************************/

INT PrepareBndVecMultipleVMPtrs (GRID *theGrid, MVM_DESC *mvmd)
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
/*D
   GetBndVecMultipleVMPtrs - get list of DOUBLE pointers for vectors and matrices

   SYNOPSIS:
   INT GetBndVecMultipleVMPtrs (const MVM_DESC *mvmd,
                                                                                         INT *cnt,
                                                                                         VECTOR *VecList[],
                                                                                         DOUBLE **vptrlist[MAXVD],
                                                                                         DOUBLE **mptrlist[MAXMD],
                                                                                         INT *vecskip, INT *vtype, INT nvec[MAXVD])


   PARAMETERS:
   .  mvmd - data filled by PrepareBndVecMultipleVMPtrs
   .  cnt - items in the vector list
   .  VecList - list of local boundary vector neighborhood
   .  vptrlist - pointer to lists of double values corresponding the  VECDATA_DESC-list
   .  mptrlist - pointer to lists of double values corresponding the  MATDATA_DESC-list
   .  vecskip - set 1 for DIRICHLET boundary, 0 else (ordering corresponds to the first
                                VECDATA_DESC)
   .  vtype - types of vectors collected for first vd
   .  nvec - number of vectors involved from this element

   DESCRIPTION:
   This functions returns pointers to the data fields described in a VECDATA_DESC-list
   and a MATDATA_DESC-list passed in the mvmd argument. Before call of this function
   in an element loop PrepareBndVecMultipleVMPtrs has to be called with mvmd.

   In the case that components are ordered subsequently only a pointer to the first
   DOUBLE value is returned. The pointers can be incremented by the user.

   RETURN VALUE:
   INT
   .n    0: ok
   .n    >0: error

   SEE ALSO:
   GetMultipleVMPtrs,PrepareBndVecMultipleVMPtrs
   D*/
/****************************************************************************/

INT GetBndVecMultipleVMPtrs (const MVM_DESC *mvmd,
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

INT ResetBndVecMultipleVMPtrs (void)
{
  if (ResetGetBoundaryNeighbourVectors())
    REP_ERR_RETURN (1);
  return (0);
}

INT FinishBndVecMultipleVMPtrs (void)
{
  if (FinishBoundaryNeighbourVectors()!=GM_OK)
    REP_ERR_RETURN (1);
  return (0);
}

/****************************************************************************/
/*D
   ClearVecskipFlags - set all vecskip flags to 0

   SYNOPSIS:
   INT ClearVecskipFlags (GRID *theGrid, VECDATA_DESC *theVD);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  theVD - type vector descriptor

   DESCRIPTION:
   This function sets the vecskip flags for all vectors of a grid to 0.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ClearVecskipFlags (GRID *theGrid, const VECDATA_DESC *theVD)
{
  VECTOR *theVector;
  INT j,n;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
  {
    n = VD_NCMPS_IN_TYPE(theVD,VTYPE(theVector));
    for (j=0; j<n; j++)
      VECSKIP(theVector) =  CLEAR_FLAG(VECSKIP(theVector),(1<<j));
  }

  return(0);
}

/****************************************************************************/
/*D
   ComputePartVecskip - set bits where skip component of vectors have to be cleared

   SYNOPSIS:
   INT ComputePartVecskip (const VECDATA_DESC *vd, const VECDATA_DESC *vds,
                                                                INT n[NVECTYPES], INT typeskip[NVECTYPES],
                                                                INT co_n[NVECTYPES], INT co_typeskip[NVECTYPES])

   PARAMETERS:
   .  vd - descriptor of the global data
   .  vds - sub descriptor of vd describing the part to be cleared
   .  n - number of components at interface
   .  typeskip - bits are set where comps of vds relative to vd
   .  co_n - number of co-components at interface
   .  co_typeskip - bits are set where not comps of vds relative to vd

   DESCRIPTION:
   This function sets bits where skip the component of vectors have to be cleared.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ComputePartVecskip (const VECDATA_DESC *vd, const VECDATA_DESC *vds,
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
/*D
   ClearPartVecskipFlags - set vecskip bits to 0 in part of the domain

   SYNOPSIS:
   INT ClearPartVecskipFlags (GRID *theGrid, const INT typeskip[NVECTYPES])

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  typeskip - clear flags where bits are set

   DESCRIPTION:
   This function resets the vecskip flags in all vectors of a grid at positions
   where flags are set in typeskip.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ClearPartVecskipFlags (GRID *theGrid, const INT typeskip[NVECTYPES])
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
/*D
   GetElementDirichletFlags - get vecskip entries

   SYNOPSIS:
   INT GetElementDirichletFlags (ELEMENT *theElement, VECDATA_DESC *theVD,
   INT *vecskip);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  vecskip - returns 1 for DIRICHLET boundary, 0 else

   DESCRIPTION:
   This function gets the vecskip flags for all local vectors
   corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetElementDirichletFlags (ELEMENT *theElement, const VECDATA_DESC *theVD,
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
/*D
   SetElementDirichletFlags - set vecskip entries

   SYNOPSIS:
   INT SetElementDirichletFlags (ELEMENT *theElement, VECDATA_DESC *theVD,
   INT *vecskip);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  vecskip - set 1 for DIRICHLET boundary, 0 else

   DESCRIPTION:
   This function sets the vecskip flags for all local vectors
   corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT SetElementDirichletFlags (ELEMENT *theElement, const VECDATA_DESC *theVD,
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
/*D
   ClearDirichletValues - sets values with Dirichlet flag to 0

   SYNOPSIS:
   INT ClearDirichletValues (GRID *theGrid, VECDATA_DESC *x);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  x - type vector descriptor

   DESCRIPTION:
   This function sets all components of x where the vecskip flag is set
   to zero.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ClearDirichletValues (GRID *theGrid, VECDATA_DESC *x)
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
/*D
   AssembleDirichletBoundary - modifies matrix entries and right hand side
   for Dirichlet boundary

   SYNOPSIS:
   INT AssembleDirichletBoundary (GRID *theGrid, MATDATA_DESC *Mat,
   VECDATA_DESC *Sol, VECDATA_DESC *Rhs);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  Mat - type matrix descriptor for the stiffness matix
   .  Sol - type vector descriptor for the solution vector
   .  Rhs - type vector descriptor for the right hand side

   DESCRIPTION:
   This function sets all components of Rhs where the vecskip flag is set
   to the Sol value and sets the corresponding matrix row to the unit vector.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT AssembleDirichletBoundary (GRID *theGrid, const MATDATA_DESC *Mat,
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
/*D
   ModifyDirichletMatrix - modifies matrix entries for Dirichlet values
   for Dirichlet boundary

   SYNOPSIS:
   INT ModifyDirichletMatrix (GRID *theGrid, const MATDATA_DESC *Mat);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  Mat - type matrix descriptor for the stiffness matix

   DESCRIPTION:
   This function sets all matrix rows to the unit vector where the vecskip flag is set.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ModifyDirichletMatrix (GRID *theGrid, const MATDATA_DESC *Mat)
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
/*D
   ModifyDirichletDefect - set defect to zero for Dirichlet boundary

   SYNOPSIS:
   INT ModifyDirichletDefect (GRID *theGrid, VECDATA_DESC *Cor);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  Cor - type vector descriptor for the correction vector

   DESCRIPTION:
   This function sets all components of Cor to zero where the vecskip flag is set.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ModifyDirichletDefect (GRID *theGrid, const VECDATA_DESC *Cor)
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
/*D
   AssembleTotalDirichletBoundary - modifies matrix entries and right hand side
   for Dirichlet boundary

   SYNOPSIS:
   INT AssembleTotalDirichletBoundary (GRID *theGrid, MATDATA_DESC *Mat,
   VECDATA_DESC *Sol, VECDATA_DESC *Rhs);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  Mat - type matrix descriptor for the stiffness matix
   .  Sol - type vector descriptor for the solution vector
   .  Rhs - type vector descriptor for the right hand side

   DESCRIPTION:
   This function sets all components of Rhs where the vecskip flag is set
   to the Sol value and sets the corresponding matrix row and column to the unit vector.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT AssembleTotalDirichletBoundary (GRID *theGrid, const MATDATA_DESC *Mat,
                                    const VECDATA_DESC *Sol,
                                    const VECDATA_DESC *Rhs)
{
  VECTOR *theVector;
  MATRIX *theMatrix;
  INT i,j,comp1,comp2,ncomp,dcomp,type,dtype;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL;
       theVector=SUCCVC(theVector)) {
    type = VTYPE(theVector);
    ncomp = VD_NCMPS_IN_TYPE (Sol,type);
    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(theVector) & (1<<j)) {
        comp1 = VD_CMP_OF_TYPE(Sol,type,j);
        comp2 = VD_CMP_OF_TYPE(Rhs,type,j);
        VVALUE(theVector,comp2) =       VVALUE(theVector,comp1);
        theMatrix = VSTART(theVector);
        for (i=0; i<ncomp; i++) {
          MVALUE(theMatrix,
                 MD_MCMP_OF_RT_CT(Mat,type,type,i*ncomp+j)) = 0.0;
          MVALUE(theMatrix,
                 MD_MCMP_OF_RT_CT(Mat,type,type,j*ncomp+i)) = 0.0;
        }
        MVALUE(theMatrix,
               MD_MCMP_OF_RT_CT(Mat,type,type,j+j*ncomp)) = 1.0;
        for (theMatrix=MNEXT(theMatrix); theMatrix!=NULL;
             theMatrix=MNEXT(theMatrix)) {
          dtype = MDESTTYPE(theMatrix);
          dcomp = VD_NCMPS_IN_TYPE (Sol,dtype);
          if (dcomp == 0) continue;
          for (i=0; i<dcomp; i++) {
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
/*D
   ConvertMatrix - converts a matrix in a sparce format

   SYNOPSIS:
   INT ConvertMatrix (GRID *theGrid, HEAP *theHeap, INT MarkKey,
   MATDATA_DESC *A, INT symmetric, int *pn, int **pia, int **pja, double **pa);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  theHeap - pointer to heap
   .  MarkKey - mark in temp mem
   .  A - pointer to matrix descriptor
   .  symmetric - if symmetric, store only upper part
   .  pn - number of lines
   .  pia - diagonal indices
   .  pja - row indices
   .  pa - matrix entries

   DESCRIPTION:
   This function converts a matrix in a sparce format.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ConvertMatrix (GRID *theGrid, HEAP *theHeap, INT MarkKey,
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

/****************************************************************************/
/*D
   PrintVector - print a vector list

   SYNOPSIS:
   INT PrintVector (GRID *g, VECDATA_DESC *X, INT vclass, INT vnclass);

   PARAMETERS:
   .  g - pointer to a grid
   .  X - pointer to vector descriptor
   .  vclass - class number
   .  vnclass - class number

   DESCRIPTION:
   This function prints the values of X for all vectors with class smaller
   or equal to vclass and next class smaller or equal to vnclass.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT PrintVectorX (const GRID *g, const VECDATA_DESC *X, INT vclass, INT vnclass, PrintfProcPtr Printf)
{
  char buffer[256];
  const VECTOR *v;
  DOUBLE_VECTOR pos;
  INT comp,ncomp,i,j;
  INT info=FALSE;

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    if (VCLASS(v) > vclass) continue;
    if (VNCLASS(v) > vnclass) continue;
    ncomp = VD_NCMPS_IN_TYPE(X,VTYPE(v));
    if (ncomp == 0) continue;
    comp = VD_CMP_OF_TYPE(X,VTYPE(v),0);
    /* Check if there is an object associated with the vector. */
    i = 0;
    if (VOBJECT(v) != NULL) {
      VectorPosition(v,pos);
      i += sprintf(buffer,"x=%5.2f y=%5.2f ",pos[0],pos[1]);
      if (DIM == 3)
        i += sprintf(buffer+i,"z=%5.2f ",pos[2]);
    }
    else {
      info = TRUE;
      i += sprintf(buffer,"                ");
      if (DIM == 3)
        i += sprintf(buffer+i,"        ");
    }
    for (j=0; j<ncomp; j++)
      i += sprintf(buffer+i,"u[%d]=%15.8lf ",j,VVALUE(v,comp+j));
    i += sprintf(buffer+i,"   cl %d %d sk ",VCLASS(v),VNCLASS(v));
    for (j=0; j<ncomp; j++)
      i += sprintf(buffer+i,"%d ",((VECSKIP(v) & (1<<j))!=0));
    i += sprintf(buffer+i,"n %d t %d o %d\n",VNEW(v),VTYPE(v),VOTYPE(v));
    Printf(buffer);

    if (Printf!=PrintDebug)
      PRINTDEBUG(np,1,("%d: %s",me,buffer));

  }

  if (info) Printf("NOTE: Geometrical information not available for some vectors.\n");

  return(NUM_OK);
}

INT PrintVector (GRID *g, VECDATA_DESC *X, INT vclass, INT vnclass)
{
  return (PrintVectorX(g,X,vclass,vnclass,UserWriteF));
}

INT PrintSVector (MULTIGRID *mg, VECDATA_DESC *X)
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
        UserWriteF("u[%d]=%15.8lf ",j,VVALUE(v,comp+j));
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
        UserWriteF("u[%d]=%15.8lf ",j,VVALUE(v,comp+j));
      UserWriteF("   cl %d %d sk ",VCLASS(v),VNCLASS(v));
      for (j=0; j<ncomp; j++)
        UserWriteF("%d ",((VECSKIP(v) & (1<<j))!=0));
      UserWriteF("\n");
    }
  }

  return(NUM_OK);
}

/****************************************************************************/
/*D
   PrintMatrix - print a matrix list

   SYNOPSIS:
   INT PrintMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass);

   PARAMETERS:
   .  g - pointer to a grid
   .  M - pointer to matrix descriptor
   .  vclass - class number
   .  vnclass - class number

   DESCRIPTION:
   This function prints the values of M of the matrix list
   for all vectors with class smaller
   or equal to vclass and next class smaller or equal to vnclass.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT PrintMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass)
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
          UserWriteF("%4.2lf ",MVALUE(m,Mcomp+j));
      }
      UserWrite("\n");
    }
  }

  return(NUM_OK);
}

INT PrintTMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass)
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
          UserWriteF("%4.2lf ",MVALUE(MADJ(m),Mcomp+j*ccomp+i));
      }
      UserWrite("\n");
    }
  }

  return(NUM_OK);
}

INT PrintDiagMatrix (GRID *g, MATDATA_DESC *Mat, INT vclass, INT vnclass)
{
  char buffer[256];
  VECTOR *v;
  DOUBLE_VECTOR pos;
  INT info=FALSE;
  MATRIX *m;
  INT Mcomp,ccomp,i,j,rtype;

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
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
    }
    else {
      info = TRUE;
      i += sprintf(buffer,"                ");
      if (DIM == 3)
        i += sprintf(buffer+i,"        ");
    }
    for (j=0; j<ccomp; j++)
      i += sprintf(buffer+i,"d[%d]=%15.8lf ",j,
                   MVALUE(m,Mcomp+j*ccomp+j));
    i += sprintf(buffer+i,"\n");
    UserWrite(buffer);
  }

  if (info)
    UserWrite("NOTE: "
              "Geometrical information not available for some vectors.\n");

  return(NUM_OK);
}

#ifdef __INTERPOLATION_MATRIX__
/****************************************************************************/
/*D
   PrintIMatrix - print the interpolation matrix list

   SYNOPSIS:
   INT PrintIMatrix (GRID *g, VECDATA_DESC *V, INT vclass, INT vnclass);

   PARAMETERS:
   .  g - pointer to a grid
   .  V - pointer to vector descriptor
   .  vclass - class number
   .  vnclass - class number

   DESCRIPTION:
   This function prints the values of the interpolation matrix list
   for all vectors with class smaller
   or equal to vclass and next class smaller or equal to vnclass.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT PrintIMatrix (GRID *g, VECDATA_DESC *V, INT vclass, INT vnclass)
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
          UserWriteF("%4.2lf ",MVALUE(m,i+j*rcomp));
      }
      UserWrite("\n");
    }
  }

  return(NUM_OK);
}
#endif  /* __INTERPOLATION_MATRIX__ */
