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

#include "switch.h"
#include "gm.h"       /* for data structure               */
#include "evm.h"      /* for data structure               */
#include "general.h"
#include "devices.h"

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
  BVP_DESC BVP_desc;
  CoeffProcPtr cpp;

  myBVP = MG_BVP(theMG);
  BVP_SetBVPDesc (myBVP,&BVP_desc);
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
  BVP_DESC BVP_desc;
  UserProcPtr upp;

  myBVP = MG_BVP(theMG);
  BVP_SetBVPDesc (myBVP,&BVP_desc);
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

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetAllVectorsOfElementOfType (ELEMENT *theElement, VECTOR **vec,
                                  const VECDATA_DESC *theVD)
{
  INT i;
  INT cnt;

  cnt = 0;
  if (VD_NCMPS_IN_TYPE(theVD,NODEVECTOR)>0)
  {
    if (GetVectorsOfNodes(theElement,&i,vec) != GM_OK)
      return(-1);
    cnt += i;
  }
  if (VD_NCMPS_IN_TYPE(theVD,EDGEVECTOR)>0)
  {
    if (GetVectorsOfEdges(theElement,&i,vec+cnt) != GM_OK)
      return(-1);
    cnt += i;
  }
  if (VD_NCMPS_IN_TYPE(theVD,ELEMVECTOR)>0)
  {
    if (GetVectorsOfElement(theElement,&i,vec+cnt) != GM_OK)
      return(-1);
    cnt += i;
  }
    #ifdef __THREEDIM__
  if (VD_NCMPS_IN_TYPE(theVD,SIDEVECTOR)>0)
  {
    if (GetVectorsOfSides(theElement,&i,vec+cnt) != GM_OK)
      return(-1);
    cnt += i;
  }
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
  INT i,j,k,l,m,cnt,vtype;
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
    vtype = VTYPE(theVec[i]);
    vncomp = VD_NCMPS_IN_TYPE(theVD,vtype);
    if (vtype == NODEVECTOR)
      if (itype[vtype] == 0)
        for (l=0; l<CORNERS_OF_SIDE(theElement,side); l++)
          for (j=0; j<vncomp; j++)
            index[k++] = m + CORNER_OF_SIDE(theElement,side,l)*vncomp + j;
    if (vtype == EDGEVECTOR)
              #ifdef __TWODIM__
      if (itype[vtype] == side)
        for (j=0; j<vncomp; j++)
          index[k++] = m + j;
          #else
      if (itype[EDGEVECTOR] == 0)
        for (l=0; l<EDGES_OF_SIDE(theElement,side); l++)
          for (j=0; j<vncomp; j++)
            index[k++] = m + EDGE_OF_SIDE(theElement,side,l)*vncomp + j;
          #endif
    if (vtype == SIDEVECTOR)
              #ifdef __TWODIM__
      if (itype[vtype] == side)
        for (j=0; j<vncomp; j++)
          index[k++] = m + j;
          #else
      if (itype[vtype] == side)
        for (j=0; j<vncomp; j++)
          index[k++] = m + j;
          #endif
    itype[vtype] += 1;
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
   DOUBLE **vptr, INT *new);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  vptr - pointer to double values
   .  new - set 1 for new vectors, 0 else

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
                        DOUBLE **vptr, INT *new)
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
      if ((new[m++] = VNEW(theVec[i])) > 0)
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
   INT GetElementVMPtrs (ELEMENT *theElement, const MATDATA_DESC *theTMD,
                                          DOUBLE **mptr);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theTMD - matrix data descriptor
   .  mptr - pointer to double values corresponding to the local stiffness matrix

   DESCRIPTION:
   This function gets all local matrix pointers corresponding to an element.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT GetElementMPtrs (ELEMENT *theElement, const MATDATA_DESC *theTMD,
                     DOUBLE **mptr)
{
  VECTOR *theVec[MAX_NODAL_VECTORS];
  MATRIX *theMatrix;
  INT vncomp[MAX_NODAL_VECTORS];
  INT vtype[MAX_NODAL_VECTORS];
  INT types[NVECTYPES];
  INT i,j,k,l,m,m1,m2,cnt;

  for (i=0; i<NVECTYPES; i++)
    types[i] = MD_ISDEF_IN_RT_CT(theTMD,i,i);

  if (GetVectorsOfTypes(theElement,types,&cnt,theVec)!=GM_OK)
    return (-1);

  if (cnt > MAX_NODAL_VECTORS || cnt < 1)
    return(-1);

  m = 0;
  for (i=0; i<cnt; i++)
  {
    vtype[i] = VTYPE(theVec[i]);
    vncomp[i] = MD_ROWS_IN_RT_CT(theTMD,vtype[i],vtype[i]);
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
                    MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[i],k*vncomp[i]+l));
    m2 = 0;
    for (j=0; j<i; j++)
    {
      if ((theMatrix = GetMatrix(theVec[i],theVec[j]))==NULL)
        return (-1);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m1+k)*m+m2+l] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[j],k*vncomp[j]+l));
      theMatrix = MADJ(theMatrix);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m2+l)*m+m1+k] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[j],l*vncomp[i]+k));
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
   VECDATA_DESC *theVD, MATDATA_DESC *theTMD,
   DOUBLE **vptr, DOUBLE **mptr);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD - type vector descriptor
   .  theTMD - type matrix descriptor
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
                      const VECDATA_DESC *theVD, const MATDATA_DESC *theTMD,
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
                    MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[i],k*vncomp[i]+l));
    m2 = 0;
    for (j=0; j<i; j++)
    {
      if ((theMatrix = GetMatrix(theVec[i],theVec[j]))==NULL)
        return (-1);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m1+k)*m+m2+l] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[j],k*vncomp[j]+l));
      theMatrix = MADJ(theMatrix);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m2+l)*m+m1+k] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[j],l*vncomp[i]+k));
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
   VECDATA_DESC *theVD2, MATDATA_DESC *theTMD,
   DOUBLE **vptr1, DOUBLE **vptr2, DOUBLE **mptr, INT *vecskip);

   PARAMETERS:
   .  theElement - pointer to an element
   .  theVD1 - type vector descriptor
   .  theVD2 - type vector descriptor
   .  theTMD - type matrix descriptor
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
                       const VECDATA_DESC *theVD2, const MATDATA_DESC *theTMD,
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
                    MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[i],k*vncomp[i]+l));
    m2 = 0;
    for (j=0; j<i; j++)
    {
      if ((theMatrix = GetMatrix(theVec[i],theVec[j]))==NULL)
        return (-3);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m1+k)*m+m2+l] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[j],k*vncomp[j]+l));
      theMatrix = MADJ(theMatrix);
      for (k=0; k<vncomp[i]; k++)
        for (l=0; l<vncomp[j]; l++)
          mptr[(m2+l)*m+m1+k] =
            MVALUEPTR(theMatrix,
                      MD_MCMP_OF_RT_CT(theTMD,vtype[i],vtype[j],l*vncomp[i]+k));
      m2 += vncomp[j];
    }
    m1 += vncomp[i];
  }

  return (m);
}

/****************************************************************************/
/*D
   ClearVecskipFlags - set all vecskip flags to 0

   SYNOPSIS:
   INT ClearVecskipFlags (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to a grid

   DESCRIPTION:
   This function sets the vecskip flags for all vectors of a grid to 0.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured
   D*/
/****************************************************************************/

INT ClearVecskipFlags (GRID *theGrid)
{
  VECTOR *theVector;

  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL;
       theVector=SUCCVC(theVector))
    VECSKIP(theVector) = 0;

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
            MVALUE(theMatrix,MD_MCMP_OF_RT_CT(Mat,type,dtype,i)) =  MVALUE(MADJ(theMatrix),MD_MCMP_OF_RT_CT(Mat,type,dtype,i)) = 0.0;
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
  VECTOR *v;
  DOUBLE pos[DIM];
  INT comp,ncomp,j;

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    if (VCLASS(v) > vclass) continue;
    if (VNCLASS(v) > vnclass) continue;
    ncomp = VD_NCMPS_IN_TYPE(X,VTYPE(v));
    if (ncomp == 0) continue;
    comp = VD_CMP_OF_TYPE(X,VTYPE(v),0);
    VectorPosition(v,pos);
    UserWriteF("x=%5.2f y=%5.2f ",pos[0],pos[1]);
    if (DIM == 3)
      UserWriteF("z=%5.2f ",pos[2]);
    for (j=0; j<ncomp; j++)
      UserWriteF("u[%d]=%15.8lf ",j,VVALUE(v,comp+j));
    UserWriteF("   cl %d %d sk ",VCLASS(v),VNCLASS(v));
    for (j=0; j<ncomp; j++)
      UserWriteF("%d ",((VECSKIP(v) & (1<<j))!=0));
    UserWrite("\n");
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
        if (rcomp != MD_ROWS_IN_RT_CT(Mat,rtype,ctype))
          UserWrite("wrong type\n");
        ccomp = MD_COLS_IN_RT_CT(Mat,rtype,ctype);
        Mcomp = MD_MCMP_OF_RT_CT(Mat,rtype,ctype,i*ccomp);
        for (j=0; j<ccomp; j++)
          UserWriteF("%4.2lf ",MVALUE(m,Mcomp+j));
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
