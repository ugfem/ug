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
   INT GetElementVMPtrs (ELEMENT *theElement, const MATDATA_DESC *md,
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
    theMatrix = START(theVec[i]);
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
    theMatrix = START(theVec[i]);
    mptr = MVALUEPTR(theMatrix,0);
    for (k=0; k<vncomp[i]; k++)
      for (l=0; l<vncomp[i]; l++)
        value[(m1+k)*m+m1+l] =
          mptr[Comp[i][i][k*vncomp[i]+l]];
    m2 = 0;
    for (j=0; j<i; j++) {
      GET_MATRIX(theVec[i],theVec[j],theMatrix);
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
      theMatrix = START(theVec[i]);
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
    theMatrix = START(theVec[i]);
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
    theMatrix = START(theVec[i]);
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
    theMatrix = START(theVec[i]);
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
   PrepareElementMultipleVMPtrs - prepare execution of GetElementMultipleVMPtrs

   SYNOPSIS:
   INT PrepareElementMultipleVMPtrs (MVM_DESC *mvmd)

   PARAMETERS:
   .  MVM_DESC - multiple VM data structure, partially to be filled before call

   DESCRIPTION:
   This function prepares the execution of GetElementMultipleVMPtrs in an element
   loop. It has to be called once before running the loop. The types needed are
   evaluated and it is checked, whether the components described by the
   XXXDATA_DESCs passed are in subsequent order. In this case ONLY a pointer to
   the first component per type is returned. The pointer can be incremented by
   the user.

   RETURN VALUE:
   INT
   .n    0

   SEE ALSO:
   GetElementMultipleVMPtrs
   D*/
/****************************************************************************/

INT PrepareElementMultipleVMPtrs (MVM_DESC *mvmd)
{
  FORMAT *fmt;
  INT tp,ctp,j,k,n,def;

  /* get format */
  if (MVMD_NVD(mvmd)>0)
    fmt = MGFORMAT(VD_MG(MVMD_VD(mvmd,0)));
  else if (MVMD_NMD(mvmd)>0)
    fmt = MGFORMAT(MD_MG(MVMD_MD(mvmd,0)));
  else
    /* no XXXDATA_DESCs defined at all */
    REP_ERR_RETURN (1);

  for (j=0; j<MVMD_NVD(mvmd); j++) MVMD_VDSUBSEQ(mvmd,j) = TRUE;
  for (j=0; j<MVMD_NMD(mvmd); j++) MVMD_MDSUBSEQ(mvmd,j) = TRUE;

  for (tp=0; tp<NVECTYPES; tp++)
  {
    def = FALSE;
    for (j=0; j<MVMD_NVD(mvmd); j++)
      if (VD_ISDEF_IN_TYPE(MVMD_VD(mvmd,j),tp))
      {
        def = TRUE;
        if (MVMD_VDSUBSEQ(mvmd,j))
        {
          /* check whether components are arranged subsequently in VECTORs of tp */
          n = VD_NCMPS_IN_TYPE(MVMD_VD(mvmd,j),tp)-1;
          for (k=0; k<n; k++)
            if (VD_CMP_OF_TYPE(MVMD_VD(mvmd,j),tp,k+1)!=VD_CMP_OF_TYPE(MVMD_VD(mvmd,j),tp,k)+1)
            {
              MVMD_VDSUBSEQ(mvmd,j) = FALSE;
              break;
            }
        }
      }
    for (j=0; j<MVMD_NMD(mvmd); j++)
      if (MD_ISDEF_IN_RT_CT(MVMD_MD(mvmd,j),tp,tp))
      {
        def = TRUE;
        if (MVMD_VDSUBSEQ(mvmd,j))
          for (ctp=0; ctp<NVECTYPES; ctp++)
          {
            /* check whether components are arranged subsequently in MATRIXs of tp,ctp */
            n = MD_ROWS_IN_RT_CT(MVMD_MD(mvmd,j),tp,ctp)*MD_COLS_IN_RT_CT(MVMD_MD(mvmd,j),tp,ctp)
                -1;
            for (k=0; k<n; k++)
              if (MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,j),tp,ctp,k+1)!=MD_MCMP_OF_RT_CT(MVMD_MD(mvmd,j),tp,ctp,k)+1)
              {
                MVMD_VDSUBSEQ(mvmd,j) = FALSE;
                break;
              }
          }
      }

    MVMD_TYPE(mvmd,tp) = def;
  }

  /* fill data and object types */
  MVMD_DATATYPES(mvmd) = MVMD_OBJTYPES(mvmd) = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (MVMD_TYPE(mvmd,tp))
    {
      MVMD_DATATYPES(mvmd) |= BITWISE_TYPE(tp);
      MVMD_OBJTYPES(mvmd)  |= FMT_T2O(fmt,tp);
    }

  return (0);
}

/****************************************************************************/
/*D
   GetElementMultipleVMPtrs - get list of DOUBLE pointers for vectors and matrices

   SYNOPSIS:
   INT GetElementMultipleVMPtrs (ELEMENT *elem, const MVM_DESC *mvmd,
                                                                                         DOUBLE **vptrlist[MAXVD],
                                                                                         DOUBLE **mptrlist[MAXMD],
                                                                                         INT *vecskip, INT *nvec)


   PARAMETERS:
   .  elem - pointer to an element
   .  mvmd - data filled by PrepareElementMultipleVMPtrs
   .  vptrlist - pointer to lists of double values corresponding the  VECDATA_DESC-list
   .  mptrlist - pointer to lists of double values corresponding the  MATDATA_DESC-list
   .  vptr1 - pointer to double values corresponding to the local right hand side
   .  vptr2 - pointer to double values corresponding to the local right hand side
   .  vecskip - set 1 for DIRICHLET boundary, 0 else (ordering corresponds to the first
                                VECDATA_DESC)
   .  nvec - number of vectors involved from this element

   DESCRIPTION:
   This functions returns pointers to the data fields described in a VECDATA_DESC-list
   and a MATDATA_DESC-list passed in the mvmd argument. Before call of this function
   in an element loop PrepareElementMultipleVMPtrs has to be called with mvmd.

   In the case that components are ordered subsequently only a pointer to the first
   DOUBLE value is returned. The pointers can be incremented by the user.

   The non-subsequent case is NOT IMPLEMENTED yet.

   RETURN VALUE:
   INT
   .n    total number of components if ok
   .n    -1: error in GetAllVectorsOfElementOfType
   .n    -2: vecdata descriptors of different size
   .n    -3: could not get matrix

   SEE ALSO:
   PrepareElementMultipleVMPtrs
   D*/
/****************************************************************************/

INT GetElementMultipleVMPtrs (ELEMENT *elem, const MVM_DESC *mvmd,
                              DOUBLE **vptrlist[MAXVD],
                              DOUBLE **mptrlist[MAXMD],
                              INT *vecskip, INT *nvec)
{
  VECTOR *theVec[MAX_NODAL_VECTORS],*rv,*cv;
  MATRIX *mat;
  INT i,j,k,l,nskip,cnt,rt,ct;
  INT vc[MAXVD],mc[MAXMD];


  if (GetVectorsOfDataTypesInObjects(elem,MVMD_DATATYPES(mvmd),MVMD_OBJTYPES(mvmd),&cnt,theVec)!=GM_OK)
    return (-1);

  *nvec = cnt;

  nskip = 0;
  for (l=0; l<MVMD_NVD(mvmd); l++) vc[l] = 0;
  for (l=0; l<MVMD_NMD(mvmd); l++) mc[l] = 0;
  for (i=0; i<cnt; i++)
  {
    rv = theVec[i];
    rt = VTYPE(rv);

    /* read skip flags */
    if (MVMD_NVD(mvmd)>1)
      for (k=0; k<VD_NCMPS_IN_TYPE(MVMD_VD(mvmd,0),rt); k++)
        vecskip[nskip++] = ((VECSKIP(rv) & (1<<k))!=0);

    /* get VECDATA_DESC ptrs */
    for (l=0; l<MVMD_NVD(mvmd); l++)
      if (VD_ISDEF_IN_TYPE(MVMD_VD(mvmd,l),rt))
        if (MVMD_VDSUBSEQ(mvmd,l))
        {
          /* by convention only return ptr to first component */
          /* TODO: possibly increase speed by introducing offset pointers
                           VVALUEPTR(rv,VD_CMP_OF_TYPE(MVMD_VD(mvmd,l),rt,0));
                           which are computed by PrepareElementMultipleVMPtrs for
                           all element types */
          vptrlist[l][vc[l]++] = VVALUEPTR(rv,VD_CMP_OF_TYPE(MVMD_VD(mvmd,l),rt,0));
        }
        else
        {
          /* fill DOUBLE pointers, subsequently all needed in rv */
          for (k=0; k<VD_NCMPS_IN_TYPE(MVMD_VD(mvmd,l),rt); k++)
            vptrlist[l][vc[l]++] = VVALUEPTR(rv,VD_CMP_OF_TYPE(MVMD_VD(mvmd,l),rt,k));
        }

    if (MVMD_NMD(mvmd)<=0)
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
    for (j=0; j<i; j++)
    {
      cv = theVec[j];
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

INT ClearVecskipFlags (GRID *theGrid, VECDATA_DESC *theVD)
{
  VECTOR *theVector;
  INT j,vtype;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL;
       theVector=SUCCVC(theVector)) {
    vtype = VTYPE(theVector);
    for (j=0; j<VD_NCMPS_IN_TYPE (theVD,vtype); j++)
      VECSKIP(theVector) =  (VECSKIP(theVector) & (~(1<<j)) | (0<<j));
  }
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
  INT i,j,comp1,comp2,ncomp,dcomp,type,dtype;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL;
       theVector=SUCCVC(theVector))
  {
    type = VTYPE(theVector);
    ncomp = VD_NCMPS_IN_TYPE (Sol,type);
    if (ncomp == 0) continue;
    for (j=0; j<ncomp; j++)
      if (VECSKIP(theVector) & (1<<j))
      {
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

INT PrintVector (GRID *g, VECDATA_DESC *X, INT vclass, INT vnclass)
{
  char buffer[256];
  VECTOR *v;
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
    UserWrite(buffer);
  }

  if (info) UserWrite("NOTE: Geometrical information not available for some vectors.\n");

  return(NUM_OK);
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
          UserWriteF("%4.2lf ",MVALUE(m,Mcomp+j));
      }
      UserWrite("\n");
    }
  }

  return(NUM_OK);
}
#endif  /* __INTERPOLATION_MATRIX__ */
